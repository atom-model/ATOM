/*
 * Ocean General Circulation Modell (OGCM) applied to laminar flow
 * program for the computation of geo-atmospherical circulating flows in a spherical shell
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz (roger.grundmann@web.de)
*/
#include "cHydrosphereModel.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "BC_Hyd.h"
#include "BC_Bath_Hyd.h"
#include "BC_Thermohalin.h"
#include "Accuracy_Hyd.h"
#include "PostProcess_Hyd.h"
#include "Pressure_Hyd.h"
#include "MinMax_Hyd.h"
#include "Results_Hyd.h"
#include "Utils.h"

#include "Config.h"
#include "tinyxml2.h"
#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

// Earth's radius is r_earth = 6731 km compares to 6.731 [/]
// for 6 km expansion of the area of circulation compares to 0.02 [/] with 40 steps of size 0.0005 
// Definition of meridional and longitudinal step sizes 
// for example: dthe = the_Grad / pi180 = 1.125 / 57.3 = 0.01963
// maximum velocity on the sea surface  w_max = 0.29 [/] compares to 0.21 m/s = 0.78 km/h as annual mean 
// mean velocity for sea level is 0.5 to 1 km/h compares to 0.14 to 0.28 m/s
// maximum temperature of earth's surface at equator t_max = 1.1355 compares to 37° C compares to 310 K
// maximum temperature of earth's surface at equator t_max = 1.0974 compares to 27° C compares to 300 K
// minimum temperature at the poles t_pol = .7803 compares to -60° C compares to 213.15 K
// minimum temperature in the deep ocean t_deep_ocean = 1.0146 compares to 4° C compares to 277.15 K
// temperature t_0 = 1.000 compares to 0° C compares to 273,15 K
// temperature t_0 = 0.003661 compares to 1° C compares to 1 K
// 1 PSU (Practical Salt Unit) = 1 g/kg, means g of salt per kg sweet water
// mass of water compares to 1.0, rate of salt compares to 0.0346
// c_0 compares to the total mass for mean salinity of 34.6 psu or dimensionsless 1.
// for c = 0.9249 compares to a salinity of 32.0 psu
// for c = 0.9682 compares to a salinity of 33.5 psu
// for c = 1.0000 compares to a salinity of 34.6 psu
// for c = 1.0983 compares to a salinity of 38.0 psu
const float cHydrosphereModel::dr = 0.025;  // 0.025 x 40 = 1.0 compares to 200m : 40 = 5m for 1 radial step
const float cHydrosphereModel::dt = 0.000001; // time step satisfies the CFL condition
const float cHydrosphereModel::pi180 = 180./ M_PI;  // pi180 = 57.3
const float cHydrosphereModel::the_degree = 1.;   // compares to 1° step size laterally
const float cHydrosphereModel::phi_degree = 1.;  // compares to 1° step size longitudinally
const float cHydrosphereModel::dthe = the_degree / pi180; // dthe = the_degree / pi180 = 1.0 / 57.3 = 0.01745, 180 * .01745 = 3.141
const float cHydrosphereModel::dphi = phi_degree / pi180; // dphi = phi_degree / pi180 = 1.0 / 57.3 = 0.01745, 360 * .01745 = 6.282

cHydrosphereModel::cHydrosphereModel():
    has_printed_welcome_msg(false){
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.
    PythonStream::OverrideCout();
    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);
    // set default configuration
    SetDefaultConfig();
}

cHydrosphereModel::~cHydrosphereModel(){}

#include "cHydrosphereDefaults.cpp.inc"

void cHydrosphereModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    try{
        if (err){
            doc.PrintError();
            throw std::invalid_argument(std::string("unable to load config file:  ") + filename);
        }
        XMLElement *atom = doc.FirstChildElement("atom"), *elem_common = NULL, *elem_hydrosphere = NULL;
        if (!atom){
            throw std::invalid_argument(std::string("Failed to find the 'atom' element in config file: ") + filename);
        }else{
            elem_common = atom->FirstChildElement("common");
            if(!elem_common){
                throw std::invalid_argument(std::string(
                    "Failed to find the 'common' element in 'atom' element in config file: ") + filename);
            }
            elem_hydrosphere = atom->FirstChildElement("hydrosphere");
            if (!elem_hydrosphere){
                throw std::invalid_argument(std::string(
                    "Failed to find the 'hydrosphere' element in 'atom' element in config file: ") + filename);
            }
        }
        #include "HydrosphereLoadConfig.cpp.inc"
    }catch(const std::exception &exc){
        std::cerr << exc.what() << std::endl;
        abort();
    }
}

void cHydrosphereModel::RunTimeSlice(int Ma){
    m_current_time = Ma;
    reset_arrays();
    mkdir(output_path.c_str(), 0777);
    int Ma_max = 300;   // parabolic temperature distribution 300 Ma back
    int Ma_max_half = 150;  // half of time scale
    cout.precision(6);
    cout.setf(ios::fixed);
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
    int switch_2D = 0;
    double emin = epsres * 100.;
    // radial expansion of the computational field for the computation of initial values
    int i_max = im - 1;   // corresponds to sea level
    int i_beg = 0;  // compares to an ocean depth of 1000 m with step size 25m, location of the thermocline
    string Name_v_w_Transfer_File;
    stringstream ssName_v_w_Transfer_File;
    string Name_SurfaceTemperature_File  = temperature_file;
    string Name_SurfaceSalinity_File = salinity_file;
    if(Ma != 0 && use_earthbyte_reconstruction){
        Name_SurfaceTemperature_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";
        Name_SurfaceSalinity_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Salinity.xyz";
        struct stat info;
        if(stat(Name_SurfaceSalinity_File.c_str(), &info) != 0){
            std::string cmd_str = "python " + reconstruction_script_path 
                + " " + std::to_string(Ma - time_step) + " " + std::to_string(Ma) 
                + " " + output_path + " " + BathymetrySuffix +" hyd";
            int ret = system(cmd_str.c_str());
            std::cout << " reconstruction script returned: " << ret << std::endl;
        }
    }
    string bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    if(!has_printed_welcome_msg)
        print_welcome_msg();
    PostProcess_Hydrosphere  read_Transfer(im, jm, km, output_path);
    read_Transfer.Atmosphere_TransferFile_read(bathymetry_name, v, w, t, 
        p_dyn, Evaporation_Dalton, Precipitation);
    cout << "***** time slice for the Oceanic Global Circulation Modell (OGCM) is:    Ma = " << Ma << " million years" 
        << endl << endl;
    cout << "***** bathymetry/topography given y the x-y-z data set:    " << bathymetry_name.c_str() << endl << endl;
    BC_Bathymetry_Hydrosphere  depth(im, jm, km);
    depth.BC_SeaGround(bathymetry_path + "/" + bathymetry_name, L_hyd, h, Bathymetry);
    BC_Hydrosphere  boundary(im, jm, km);
    Pressure_Hyd  startPressure(im, jm, km, dr, dthe, dphi);
    Results_Hyd  calculate_MSL(im, jm, km);
    BC_Thermohalin  oceanflow(im, jm, km, i_beg, i_max, Ma, Ma_max, Ma_max_half, 
        dr, g, r_0_water, u_0, p_0, t_0, c_0, cp_w, 
        L_hyd, t_average, t_paleo_max, t_equator, t_pole, output_path);
//    goto Printout;
    if(use_NASA_velocity){
        read_IC(velocity_v_file, v.x[im-1], jm, km);
        read_IC(velocity_w_file, w.x[im-1], jm, km);    
        read_IC(Name_SurfaceTemperature_File, t.x[im-1], jm, km);
        read_IC(Name_SurfaceSalinity_File, c.x[im-1], jm, km);
    }
    iter_cnt_3d = -1;
    if(debug) save_data();
    iter_cnt_3d++;
//    goto Printout;
    IC_EkmanSpiral();
//    goto Printout;
//    IC_u_WestEastCoast();
//    IC_Equatorial_Currents();
//    if(Ma <= 41)  IC_CircumPolar_Current(); // Drake passage closed 41 Ma ago
//    oceanflow.BC_Temperature_Salinity(h, t, c, p_dyn );
    oceanflow.BC_Temperature_Salinity(h, t, c, p_dyn, Evaporation_Dalton, Precipitation, salinity_evaporation, c_fix );
//    BC_Pressure_Density();
//    goto Printout;
    store_intermediate_data_2D();
    store_intermediate_data_3D();
    calculate_MSL.land_oceanFraction (h);
    /*******************************************   start of pressure and velocity iterations ******************************/
    iter_cnt = 1;
    /** ::::::::::::::::::::::::::::::::::::::   begin of 2D loop for initial surface conditions: if (switch_2D == 0)   ::::::**/
    if(switch_2D != 1){
        for(int pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++){
            for(int velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++){
/*
                cout << endl << endl;
                cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
                cout << " 2D OGCM iterational process" << endl;
                cout << " max total iteration number nm = " << nm << endl << endl;

                cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, \
                    maximum and current number of velocity iterations, maximum and current number of pressure iterations " 
                    << endl << endl << " Ma = " << Ma << "     n = " << iter_cnt << "    velocity_iter_max_2D = " << 
                    velocity_iter_max_2D << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " 
                    << pressure_iter_max_2D << "    pressure_iter_2D = " << pressure_iter_2D << endl;
*/
                boundary.RB_theta(t, u, v, w, p_dyn, c);
                boundary.RB_phi(t, u, v, w, p_dyn, c);
                oceanflow.Value_Limitation_Hyd(h, u, v, w, p_dyn, t, c);
                run_data(); 
                solveRungeKutta_2D_Hydrosphere();
                store_intermediate_data_2D();
                iter_cnt++;
            } // end of velocity loop_2D: if(velocity_iter_2D > velocity_iter_max_2D)
            startPressure.computePressure_2D(u_0, r_0_water, rad, the, p_dyn, p_dynn, h, aux_v, aux_w);
            if(iter_cnt > nm){
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        } // end of pressure loop_2D: if (pressure_iter_2D > pressure_iter_max_2D) 
        iter_cnt = 1;
        switch_2D = 1;
        emin = epsres * 100.;
    }// end of 2D loop for initial surface conditions: if (switch_2D == 0)   ::::::::::
    cout << endl << endl;
    iter_cnt_3d = 0;
    if(debug) save_data();
    for(int pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++){
        for(int velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++){
            cout << endl << endl;
            cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            cout << " 3D OGCM iterational process" << endl;
            cout << " max total iteration number nm = " << nm << endl << endl;
            cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum \
                and current number of velocity iterations, maximum and current number of pressure iterations " << endl 
                << endl << " Ma = " << Ma << "     n = " << iter_cnt << "    velocity_iter_max = " << velocity_iter_max << 
                "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << 
                "    pressure_iter = " << pressure_iter << endl;
            boundary.RB_radius(dr, rad, t, u, v, w, p_dyn, c);
            boundary.RB_theta(t, u, v, w, p_dyn, c);
            boundary.RB_phi(t, u, v, w, p_dyn, c);
            depth.BC_SolidGround(h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn);
            BC_Pressure_Density();
            oceanflow.Value_Limitation_Hyd(h, u, v, w, p_dyn, t, c);
            solveRungeKutta_3D_Hydrosphere(); 
            print_min_max();
            run_data(); 
            store_intermediate_data_3D();
            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        } // end of velocity loop_3D: if(velocity_iter > velocity_iter_max)
        startPressure.computePressure_3D(u_0, r_0_water, rad, the, p_dyn, p_dynn, h, aux_u, aux_v, aux_w);

        if(debug && pressure_iter % checkpoint == 0){
            write_file(bathymetry_name, output_path);
        }
        if(iter_cnt > nm){
            cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
            break;
        }
    }// end of pressure loop_3D: if(pressure_iter > pressure_iter_max)   :::::::::::
//    Printout:
    cout << endl << endl;
    write_file(bathymetry_name, output_path, true);
    iter_cnt_3d++;
    save_data();
    cout << endl << "***** end of the Hydrosphere General Circulation Modell (OGCM) *****" << endl << endl;
    if(emin <= epsres)  cout << "***** steady solution reached! *****" << endl;
}


void cHydrosphereModel::reset_arrays(){
    rad.initArray_1D(im, 1.); // radial coordinate direction
    the.initArray_1D(jm, 2.); // lateral coordinate direction
    phi.initArray_1D(km, 3.); // longitudinal coordinate direction
    aux_grad_v.initArray_1D(km, 4.); // auxilliar array
    aux_grad_w.initArray_1D(km, 5.); // auxilliar array

    Bathymetry.initArray_2D(jm, km, 0.); // Bathymetry in m
    Upwelling.initArray_2D(jm, km, 0.); // upwelling
    Downwelling.initArray_2D(jm, km, 0.); // downwelling
    EkmanPumping.initArray_2D(jm, km, 0.); // 2D bottom water summed up in a vertical column
    SaltFinger.initArray_2D(jm, km, 0.);      // salt bulge of higher density
    SaltDiffusion.initArray_2D(jm, km, 0.);   // salt bulge of lower density
    Salt_total.initArray_2D(jm, km, 0.);     // rate of salt summed up in a vertical column
    BuoyancyForce_2D.initArray_2D(jm, km, 0.); // radiation balance at the surface
    salinity_evaporation.initArray_2D(jm, km, 0.); // additional salinity by evaporation
    Evaporation_Dalton.initArray_2D(jm, km, 0.); // evaporation by Dalton in [mm/d]
    Precipitation.initArray_2D(jm, km, 0.); // areas of higher precipitation
    c_fix.initArray_2D(jm, km, 0.); // local surface salinity fixed for iterations

    h.initArray(im, jm, km, 0.); // bathymetry, depth from sea level
    t.initArray(im, jm, km, 1.); // temperature
    u.initArray(im, jm, km, 0.); // u-component velocity component in r-direction
    v.initArray(im, jm, km, 0.); // v-component velocity component in theta-direction
    w.initArray(im, jm, km, 0.); // w-component velocity component in phi-direction
    c.initArray(im, jm, km, 1.); // salinity

    tn.initArray(im, jm, km, 1.); // temperature new
    un.initArray(im, jm, km, 0.); // u-velocity component in r-direction new
    vn.initArray(im, jm, km, 0.); // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, 0.); // w-velocity component in phi-direction new
    cn.initArray(im, jm, km, 1.); // salinity new

    p_dyn.initArray(im, jm, km, 1.); // dynamic pressure
    p_dynn.initArray(im, jm, km, 1.); // dynamic pressure new
    p_stat.initArray(im, jm, km, 1.); // static pressure

    rhs_t.initArray(im, jm, km, 0.); // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.); // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.); // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.); // auxilliar field RHS w-velocity component
    rhs_c.initArray(im, jm, km, 0.); // auxilliar field RHS water vapour

    aux_u.initArray(im, jm, km, 0.); // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.); // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.); // auxilliar field w-velocity component

    Salt_Finger.initArray(im, jm, km, 0.); // salt bulge of higher density
    Salt_Diffusion.initArray(im, jm, km, 0.); // salt bulge of lowerer density and temperature
    Salt_Balance.initArray(im, jm, km, 0.); // +/- salt balance

    r_water.initArray(im, jm, km, 0.); // water density as function of pressure
    r_salt_water.initArray(im, jm, km, 0.); // salt water density as function of pressure and temperature
    BuoyancyForce_3D.initArray(im, jm, km, 0.); // 3D buoyancy force
}

void cHydrosphereModel::Run(){
    mkdir(output_path.c_str(), 0777);
    cout << "Output is being written to " << output_path << "\n";
    // write out the config for reproducibility
    // disabled for now
    // std::stringstream output_config_path;
    // output_config_path << output_path << "/config_hyd.xml";
    // WriteConfig(output_config_path.str().c_str());
    if(verbose){
        cout << endl << endl << endl;
        cout << "***** Hydrosphere General Circulation Model (OGCM) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 1 additional transport equations to describe the salinity" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** salinity is part of the Boussinesq approximation" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz (roger.grundmann@web.de)" << endl << endl;
        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
    }
    for(int i = time_start; i <= time_end; i+=time_step){
        RunTimeSlice(i);
    }
    cout << endl << "***** end of the Hydrosphere General Circulation Modell (OGCM) *****" << endl << endl;
    cout << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-hydrospheric circulation *****";
    cout << "\n\n\n\n";
}

void cHydrosphereModel::write_file(std::string &bathymetry_name, string& filepath, bool is_final_result){
    PostProcess_Hydrosphere  write_File(im, jm, km, filepath);
    int j_longal = 75;
    write_File.paraview_vtk_longal(bathymetry_name, j_longal, iter_cnt-1, u_0, c_0, r_0_water, h, p_dyn, p_stat, r_water, 
        r_salt_water, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance);
    int k_zonal = 185;
    write_File.paraview_vtk_zonal(bathymetry_name, k_zonal, iter_cnt-1, u_0, c_0, r_0_water, h, p_dyn, p_stat, r_water, r_salt_water, 
            t, u, v, w, c, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance);
    int i_radial = 40;
    write_File.paraview_vtk_radial(bathymetry_name, i_radial, iter_cnt-1, u_0, t_0, c_0, r_0_water, h, p_dyn, p_stat, r_water, 
        r_salt_water, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance, 
        Upwelling, Downwelling, SaltFinger, SaltDiffusion, BuoyancyForce_2D, EkmanPumping, Evaporation_Dalton, 
        Precipitation, Bathymetry);
    if(paraview_panorama_vts_flag){
        write_File.paraview_panorama_vts(bathymetry_name, iter_cnt-1, u_0, c_0, r_0_water, h, t, p_dyn, p_stat, r_water, r_salt_water, 
            u, v, w, c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance);
    }
    PostProcess_Hydrosphere  ppa(im, jm, km, output_path);
    Hydrosphere_PlotData(bathymetry_name, (is_final_result ? -1 : iter_cnt-1));
}


void  cHydrosphereModel::save_data(){
    struct stat info;
    string path = output_path + "/bin_data/";
    if(stat(path.c_str(), &info) != 0){
        mkdir(path.c_str(), 0777);
    }
    std::ostringstream ss;
    if(iter_cnt_3d == pressure_iter_max * velocity_iter_max+1)
        ss << "_time_" << m_current_time << "_iter_n";
    else
        ss << "_time_" << m_current_time << "_iter_" << iter_cnt_3d;
    std::string postfix_str = ss.str();
    Array  v_t(im, jm, km, 0), w_t(im, jm, km, 0), c_t(im, jm, km, 0), t_t(im, jm, km, 0);
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                c_t.x[i][j][k] = c.x[i][j][k] * 35.;
                v_t.x[i][j][k] = v.x[i][j][k] * u_0;
                w_t.x[i][j][k] = w.x[i][j][k] * u_0;
                t_t.x[i][j][k] = t.x[i][j][k] * 273.15 - 273.15;
            }
        }
    }

    for(int i=im-1; i>19; i--){
        v_t.save(path + std::string("hyd_v")+postfix_str, i);
        w_t.save(path + std::string("hyd_w")+postfix_str, i);
    }

    c_t.save(path + std::string("hyd_s")+postfix_str, im-1);
    t_t.save(path + std::string("hyd_t")+postfix_str, im-1);
}

void cHydrosphereModel::store_intermediate_data_3D(float coeff){
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                un.x[i][j][k] = coeff * u.x[i][j][k];
                vn.x[i][j][k] = coeff * v.x[i][j][k];
                wn.x[i][j][k] = coeff * w.x[i][j][k];
                p_dynn.x[i][j][k] = coeff * p_dyn.x[i][j][k];
                tn.x[i][j][k] = coeff * t.x[i][j][k];
                cn.x[i][j][k] = coeff * c.x[i][j][k];
            }
        }
    }
}

void cHydrosphereModel::store_intermediate_data_2D(float coeff){   
    for(int j = 0; j < jm; j++){   
        for(int k = 0; k < km; k++){   
            vn.x[im-1][j][k] = coeff * v.x[im-1][j][k];
            wn.x[im-1][j][k] = coeff * w.x[im-1][j][k]; 
            p_dynn.x[im-1][j][k] = coeff * p_dyn.x[im-1][j][k];
        }
    }
}

