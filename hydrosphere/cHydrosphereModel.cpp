/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * program for the computation of geo-atmospherical circulating flows in a spherical shell
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)
*/

#include <fenv.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>    
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include "Accuracy_Hyd.h"
#include "Utils.h"
#include "Config.h"
#include "AtomMath.h"
#include "cHydrosphereModel.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

cHydrosphereModel* cHydrosphereModel::m_model = NULL;

// Earth's radius is r_earth = 6731 km compares to 6.731 [/]
// for 6 km expansion of the area of circulation compares to 0.02 [/] with 40 steps of size 0.0005 
// Definition of meridional and longitudinal step sizes 
// for example: dthe = the_Grad/pi180 = 1.125/57.3 = 0.01963
// maximum velocity on the sea surface  w_max = 0.29 [/] compares to 0.21 m/s = 0.78 km/h as annual mean 
// mean velocity for sea level is 0.5 to 1 km/h compares to 0.14 to 0.28 m/s
// maximum temperature of earth's surface at equator t_max = 1.1355 compares to 37° C compares to 310 K
// maximum temperature of earth's surface at equator t_max = 1.0974 compares to 27° C compares to 300 K
// minimum temperature at the poles t_pol = .7803 compares to -60° C compares to 213.15 K
// minimum temperature in the deep ocean t_deep_ocean = 1.0146 compares to 4° C compares to 277.15 K
// temperature t_0 = 1.000 compares to 0° C compares to 273,15 K
// temperature t_0 = 0.003661 compares to 1° C compares to 1 K
// 1 PSU(Practical Salt Unit) = 1 g/kg, means g of salt per kg sweet water
// mass of water compares to 1.0, rate of salt compares to 0.0346
// c_0 compares to the total mass for mean salinity of 34.6 psu or dimensionsless 1.
// for c = 0.9249 compares to a salinity of 32.0 psu
// for c = 0.9682 compares to a salinity of 33.5 psu
// for c = 1.0000 compares to a salinity of 34.6 psu
// for c = 1.0983 compares to a salinity of 38.0 psu
const float cHydrosphereModel::dr = 0.025;  // 0.025 x 40 = 1.0 compares to 200m : 40 = 5m for 1 radial step
//const float cHydrosphereModel::dr = - 0.025;  // 0.025 x 40 = 1.0 compares to 200m : 40 = 5m for 1 radial step
//const float cHydrosphereModel::dt = 0.00001; // time step satisfies the CFL condition
//const float cHydrosphereModel::dt = 0.0001; // time step satisfies the CFL condition
//const float cHydrosphereModel::dt = 0.00005; // time step satisfies the CFL condition
const float cHydrosphereModel::dt = 0.00001; // time step satisfies the CFL condition
const float cHydrosphereModel::pi180 = 180./M_PI;  // pi180 = 57.3
const float cHydrosphereModel::the_degree = 1.;   // compares to 1° step size laterally
const float cHydrosphereModel::phi_degree = 1.;  // compares to 1° step size longitudinally
const float cHydrosphereModel::dthe = the_degree/pi180; // dthe = the_degree/pi180 = 1.0/57.3 = 0.01745, 180 * .01745 = 3.141
const float cHydrosphereModel::dphi = phi_degree/pi180; // dphi = phi_degree/pi180 = 1.0/57.3 = 0.01745, 360 * .01745 = 6.282
    
const double cHydrosphereModel::the0 = 0.;             // North Pole
const double cHydrosphereModel::phi0 = 0.;             // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cHydrosphereModel::r0 = - 1.; // non-dimensional
//const double cHydrosphereModel::r0 = 6731000.; // in m

cHydrosphereModel::cHydrosphereModel():
    i_bathymetry(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    has_printed_welcome_msg(false){
    if(PythonStream::is_enable()){
        backup = std::cout.rdbuf();
        std::cout.rdbuf(&ps);
    }
    signal(SIGINT, exit);
    SetDefaultConfig();
    m_model = this;
    rad.initArray_1D(im, 0);
    the.initArray_1D(jm, 0);
    phi.initArray_1D(km, 0);
}

cHydrosphereModel::~cHydrosphereModel(){
    if(PythonStream::is_enable()){
        std::cout.rdbuf(backup);
    }
    m_model = NULL;
    logger().close();
}

#include "cHydrosphereDefaults.cpp.inc"

void cHydrosphereModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    try{
        if(err){
            doc.PrintError();
            throw std::invalid_argument(std::string("unable to load config file:  ") + filename);
        }
        XMLElement *atom = doc.FirstChildElement("atom"), *elem_common = NULL, *elem_hydrosphere = NULL;
        if(!atom){
            throw std::invalid_argument(std::string("Failed to find the 'atom' element in config file: ") + filename);
        }else{
            elem_common = atom->FirstChildElement("common");
            if(!elem_common){
                throw std::invalid_argument(std::string(
                    "Failed to find the 'common' element in 'atom' element in config file: ") + filename);
            }
            elem_hydrosphere = atom->FirstChildElement("hydrosphere");
            if(!elem_hydrosphere){
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
/*
*
*/
void cHydrosphereModel::RunTimeSlice(int Ma){
    if(!is_temperature_curve_loaded()) 
        load_temperature_curve();
    reset_arrays();
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
    m_current_time = m_time_list.insert(float(Ma)).first;
    mkdir(output_path.c_str(), 0777);
    Ma = int(round(*get_current_time()));
    cout.precision(6);
    cout.setf(ios::fixed);
    string Name_v_w_Transfer_File;
    stringstream ssName_v_w_Transfer_File;
    string Name_SurfaceTemperature_File = temperature_file;
    string Name_SurfaceNASATemperature_File  = temperature_file;
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
    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    if(!has_printed_welcome_msg)
        print_welcome_msg();
    Atmosphere_v_w_Transfer(bathymetry_name);
    cout << "***** time slice for the Oceanic Global Circulation Modell(OGCM) is:    Ma = " << Ma << " million years" 
        << endl << endl;
    cout << "***** bathymetry/topography given by the x-y-z data set:    " << bathymetry_name.c_str() << endl << endl;
    init_bathymetry(bathymetry_path + "/" + bathymetry_name);
//    goto Printout;
    if(use_NASA_velocity){
        read_IC(velocity_v_file, v.x[im-1], jm, km);
        read_IC(velocity_w_file, w.x[im-1], jm, km);    
        read_IC(Name_SurfaceTemperature_File, t.x[im-1], jm, km);
        read_IC(Name_SurfaceNASATemperature_File, temperature_NASA.y, jm, km);
        read_IC(Name_SurfaceSalinity_File, c.x[im-1], jm, km);
    }
    iter_cnt = 1;
    iter_cnt_3d = -1;
    if(debug) save_data();
    iter_cnt_3d++;
    land_oceanFraction();
    EkmanSpiral();
    Value_Limitation_Hyd();
//    goto Printout;
//    IC_u_WestEastCoast();
//    IC_Equatorial_Currents();
//    if(Ma <= 41)  IC_CircumPolar_Current(); // Drake passage closed 41 Ma ago
    init_temperature();
    IC_t_WestEastCoast();
//    goto Printout;
    fft_gaussian_filter_3d(t,1);
    init_salinity();
//    goto Printout;
    PresStat_SaltWaterDens();
    SalinityEvaporation();
    init_dynamic_pressure();
    fft_gaussian_filter_3d(c,1);
//    store_intermediate_data_2D(1.);
    store_intermediate_data_3D(1.);
//    run_2D_loop();
//    fft_gaussian_filter_3d(u,1);
//    fft_gaussian_filter_3d(v,1);
//    fft_gaussian_filter_3d(w,1);
//    goto Printout;
    run_data_hyd(); 
    print_min_max_hyd();
    run_3D_loop();
//    goto Printout;
    cout << endl << endl;
/*
    Printout:
    run_data_hyd(); 
    print_min_max_hyd();
    write_file(bathymetry_name, output_path, true);
*/
    cout << endl << endl;
    iter_cnt_3d++;
    save_data();
    cout << endl << "***** end of the Hydrosphere General Circulation Modell(OGCM) *****" << endl << endl;
}
/*
*
*/
void cHydrosphereModel::reset_arrays(){
    rad.initArray_1D(im, 1.); // radial coordinate direction
    the.initArray_1D(jm, 2.); // lateral coordinate direction
    phi.initArray_1D(km, 3.); // longitudinal coordinate direction
    aux_grad_v.initArray_1D(im, 4.); // auxilliar array
    aux_grad_w.initArray_1D(im, 5.); // auxilliar array

    Bathymetry.initArray_2D(jm, km, 0.); // Bathymetry in m
    Upwelling.initArray_2D(jm, km, 0.); // upwelling
    Downwelling.initArray_2D(jm, km, 0.); // downwelling
    EkmanPumping.initArray_2D(jm, km, 0.); // 2D Ekman pumping vertical velocity summed up in a vertical column
    SaltFinger.initArray_2D(jm, km, 0.);      // salt bulge of higher density
    SaltDiffusion.initArray_2D(jm, km, 0.);   // salt bulge of lower density
    Salt_total.initArray_2D(jm, km, 0.);     // rate of salt summed up in a vertical column
    BuoyancyForce_2D.initArray_2D(jm, km, 0.); // radiation balance at the surface
    salinity_evaporation.initArray_2D(jm, km, 0.); // additional salinity by evaporation
    Evaporation_Dalton.initArray_2D(jm, km, 0.); // evaporation by Dalton in [mm/d]
    Evaporation_Penman.initArray_2D(jm, km, 0.); // evaporation by Penman in [mm/d]
    Precipitation.initArray_2D(jm, km, 0.); // areas of higher precipitation
    temperature_NASA.initArray_2D(jm, km, 0.); // surface temperature from NASA
    c_fix.initArray_2D(jm, km, 0.); // local surface salinity fixed for iterations
    v_wind.initArray_2D(jm, km, 0.); // v-component of surface wind
    w_wind.initArray_2D(jm, km, 0.); // w-component of surface wind

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

    p_dyn.initArray(im, jm, km, 0.); // dynamic pressure
    p_dynn.initArray(im, jm, km, 0.); // dynamic pressure new
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

    r_water.initArray(im, jm, km, r_0_water); // water density as function of pressure
    r_salt_water.initArray(im, jm, km, r_0_water); // salt water density as function of pressure and temperature
    BuoyancyForce.initArray(im, jm, km, 0.); // 3D buoyancy force
    CoriolisForce.initArray(im, jm, km, 0.); // Coriolis force terms
    PressureGradientForce.initArray(im, jm, km, 0.); // Force caused by normal pressure gradient
}
/*
*
*/
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
        cout << "***** Hydrosphere General Circulation Model(OGCM) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 1 additional transport equations to describe the salinity" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** salinity is part of the Boussinesq approximation" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)" << endl << endl;
        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
    }
    for(int i = time_start; i <= time_end; i+=time_step){
        RunTimeSlice(i);
    }
    cout << endl << "***** end of the Hydrosphere General Circulation Modell(OGCM) *****" << endl << endl;
    cout << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-hydrospheric circulation *****";
    cout << "\n\n\n\n";
}
/*
*
*/
void cHydrosphereModel::write_file(std::string &bathymetry_name, 
    std::string &output_path, bool is_final_result){
    int i_radial = 40;
    paraview_vtk_radial(bathymetry_name, i_radial, iter_cnt-1);
    int j_longal = 75;
    paraview_vtk_longal(bathymetry_name, j_longal, iter_cnt-1);
    int k_zonal = 185;
    paraview_vtk_zonal(bathymetry_name, k_zonal, iter_cnt-1);
    if(paraview_panorama_vts_flag){
        paraview_panorama_vts(bathymetry_name, iter_cnt-1);
    }
    Hydrosphere_PlotData(bathymetry_name,(is_final_result ? -1 : iter_cnt-1));
}
/*
*
*/
void  cHydrosphereModel::save_data(){
    struct stat info;
    string path = output_path + "/bin_data/";
    if(stat(path.c_str(), &info) != 0){
        mkdir(path.c_str(), 0777);
    }
    std::ostringstream ss;
    if(iter_cnt_3d == pressure_iter_max * velocity_iter_max+1)
        ss << "_time_" << (int)(*get_current_time()) << "_iter_n";

    else
        ss << "_time_" << (int)(*get_current_time()) << "_iter_" << iter_cnt_3d;
    std::string postfix_str = ss.str();
    Array  v_t(im, jm, km, 0), w_t(im, jm, km, 0), c_t(im, jm, km, 0), t_t(im, jm, km, 0);
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                c_t.x[i][j][k] = c.x[i][j][k] * c_0;
                v_t.x[i][j][k] = v.x[i][j][k] * u_0;
                w_t.x[i][j][k] = w.x[i][j][k] * u_0;
                t_t.x[i][j][k] = t.x[i][j][k] * t_0 - t_0;
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
/*
*
*/
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
/*
*
*/
void cHydrosphereModel::store_intermediate_data_2D(float coeff){   
    for(int j = 0; j < jm; j++){   
        for(int k = 0; k < km; k++){   
            vn.x[im-1][j][k] = coeff * v.x[im-1][j][k];
            wn.x[im-1][j][k] = coeff * w.x[im-1][j][k]; 
            p_dynn.x[im-1][j][k] = coeff * p_dyn.x[im-1][j][k];
        }
    }
}
/*
*
*/
void  cHydrosphereModel::run_2D_loop(){
    int switch_2D = 0;    
    iter_cnt = 1;
    if(switch_2D != 1){
        for(int pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++){
            for(int velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++){

                cout << endl << endl;
                cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
                cout << " 2D OGCM iterational process" << endl;
                cout << " max total iteration number nm = " << nm << endl << endl;
                cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, \
                    maximum and current number of velocity iterations, maximum and current number of pressure iterations " 
                    << endl << endl << " Ma = " << (int)*get_current_time() << "     n = " << iter_cnt << "    velocity_iter_max_2D = " << 
                    velocity_iter_max_2D << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " 
                    << pressure_iter_max_2D << "    pressure_iter_2D = " << pressure_iter_2D << endl;
                BC_theta();
                BC_phi();
                BC_SolidGround();
                solveRungeKutta_2D_Hydrosphere();
                store_intermediate_data_2D();
//                Value_Limitation_Hyd();
                iter_cnt++;
            } // end of velocity loop
            computePressure_2D();
            print_min_max_hyd();
            run_data_hyd(); 
            store_intermediate_data_2D();
            if(iter_cnt > nm){
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        } // end of pressure loop
    }// end of 2D loop
    return;
}
/*
*
*/
void  cHydrosphereModel::run_3D_loop(){
cout << endl << "      run_3D_loop hyd" << endl;
    iter_cnt = 1;
    iter_cnt_3d = 0;
    for(int pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++){
        for(int velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++){
            cout << endl << endl;
            cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            cout << " 3D OGCM iterational process" << endl;
            cout << " max total iteration number nm = " << nm << endl << endl;
            cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum \
                and current number of velocity iterations, maximum and current number of pressure iterations " << endl 
                << endl << " Ma = " << (int)*get_current_time() << "     n = " << iter_cnt << "    velocity_iter_max = " << velocity_iter_max << 
                "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << 
                "    pressure_iter = " << pressure_iter << endl;
            BC_radius();
            BC_theta();
            BC_phi();
            BC_SolidGround();
//            EkmanSpiral();
            if(velocity_iter % 2 == 0){
//                SalinityEvaporation();
                PresStat_SaltWaterDens();
            }
            solveRungeKutta_3D_Hydrosphere(); 
            run_data_hyd();
            Value_Limitation_Hyd();
            print_min_max_hyd();
            store_intermediate_data_3D(1.);
            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        } // end of velocity loop
        computePressure_3D();
/*
    p_dyn.printArray(im, jm, km);
    aux_u.printArray(im, jm, km);
    aux_v.printArray(im, jm, km);
    aux_w.printArray(im, jm, km);
    rhs_u.printArray(im, jm, km);
    rhs_v.printArray(im, jm, km);
    rhs_w.printArray(im, jm, km);
*/
        if(pressure_iter % checkpoint == 0){
            write_file(bathymetry_name, output_path, true);
        }
        if(iter_cnt > nm){
            cout << "       nm = " << nm 
                << "     .....     maximum number of iterations   nm   reached!" 
                << endl;
            break;
        }
    }// end of pressure loop
cout << endl << "      run_3D_loop hyd ended" << endl;
}
/*
*
*/
void cHydrosphereModel::EkmanSpiral(){
    //Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
    //a further turning downwards until the end of the shear layer such that finally 90° of turning are reached
//    double wind_water_vel_ratio = 0.03;
    // initial conditions for v and w velocity components at the sea surface
    // ocean surface velocity is about 3% of the wind velocity at the surface
    // u_0 for the hydrosphere is (0.03 * u_0) for the atmosphere
    // surface wind vector driving the Ekman spiral in the Ekman layer
    // northern and southern hemisphere
    cout.precision(8);
    int i_Ekman = 0;
    int i_max = im-1;
    double a = 0.; // a in 1/m
    double Ekman_angle = 45./pi180;
    double sinthe = 0.;
    double rm = 0.;
    double rmsinthe = 0;
    double a_z = 0.;
    double exp_a_z = 0.;
    double sin_a_z = 0;
    double cos_a_z = 0;
    double alfa = 0;
    double angle = 0;
    double CD = 2.6e-3;  // drag coefficient in ./.
    double U_10 = 0.;  // wind velocity 10 m above sea surface directed to the north in m/s
    double V_0 = 0.;  // water velocity at surface shifted by Ekman angle in m/s
    double T_yz = 0.;  // wind stress in v-direction ( y, j ) in kg/(m*s2)
    double f = 0.;  // Coriolis parameter in 1/s
    double Az = 0.;  // constant vertical eddy viscosity in m2/s
    double D_E = 0.;  // Ekman layer depth in m
    double D_E_op = 0.;  // Ekman layer depth opposite to the surface wind in m
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){ // in m/s
            v_wind.y[j][k] = v.x[im-1][j][k];
            w_wind.y[j][k] = w.x[im-1][j][k];
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 1; j < jm-1; j++){
            v_wind.y[90][k] = 0.;
            if(w_wind.y[j][k] == 0.) w_wind.y[j][k] = 1.e-8;
            alfa = atan(fabs(v_wind.y[j][k]/w_wind.y[j][k]));
//            sinthe = sin(fabs((90. - (double)j) * M_PI/180.));
            if(j <= 90)
                sinthe = sin(the.z[j]);
            if(j > 90){
                int j_rev = 90 - fabs(j - 90);
                sinthe = sin(the.z[j_rev]);
            }
            if(sinthe == 0.)  sinthe = 1.e-5;
            U_10 = sqrt(v_wind.y[j][k] * v_wind.y[j][k] 
                + w_wind.y[j][k] * w_wind.y[j][k]); // in m/s, dimensional surface wind velocity U_10 in m/s
                // original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 139, eq. 9.16
            D_E = 7.6/sqrt(sinthe) * U_10; // in m
            i_Ekman = int(D_E/L_hyd * (double)i_max);
            T_yz = r_air * CD * U_10 * U_10;  // in kg/(m*s*s))
            V_0 = 0.0127 * U_10/sqrt(sinthe);  // in m/s
            f = 2. * omega * fabs(sinthe);  // in 1/s
            Az = pow((T_yz/(r_0_water * V_0)),2)/f;  // in m*m/s
            a = sqrt(f/(2. * Az));  // in 1/m
            D_E_op = sqrt(2. * Az * M_PI * M_PI/f);  // in m
/*
    if((j == 75) &&(k == 180)) cout << endl << "Ekman-Layer north" << endl
        << "   j = " << j << "   k = " << k << endl
        << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
        << "   v_wind = " << v_wind.y[j][k] 
        << "   w_wind = " << w_wind.y[j][k] << endl
        << "   v = " << v.x[im-1][j][k] 
        << "   w = " << w.x[im-1][j][k] << endl
        << "   Ekman_angle = " << Ekman_angle
        << "   Ekman_angle_deg = " << Ekman_angle * pi180 << endl
        << "   alfa = " << alfa
        << "   alfa_deg = " << alfa * pi180 << endl
        << "   U_10 = " << U_10 
        << "   V_0 = " << V_0 << endl
        << "   D_E = " << D_E << "   D_E_op = " << D_E_op << endl
        << "   T_yz = " << T_yz << endl
        << "   f = " << f
        << "   Az = " << Az
        << "   a = " << a
        << "   i_Ekman = " << i_Ekman 
        << endl << endl;

    if((j == 105) &&(k == 180)) cout << endl << "Ekman-Layer south" << endl
        << "   j = " << j << "   k = " << k << endl
        << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
        << "   v_wind = " << v_wind.y[j][k] 
        << "   w_wind = " << w_wind.y[j][k] << endl
        << "   v = " << v.x[im-1][j][k] 
        << "   w = " << w.x[im-1][j][k] << endl
        << "   Ekman_angle = " << Ekman_angle
        << "   Ekman_angle_deg = " << Ekman_angle * pi180 << endl
        << "   alfa = " << alfa
        << "   alfa_deg = " << alfa * pi180 << endl
        << "   U_10 = " << U_10 
        << "   V_0 = " << V_0 << endl
        << "   D_E = " << D_E << "   D_E_op = " << D_E_op << endl
        << "   T_yz = " << T_yz << endl
        << "   f = " << f
        << "   Az = " << Az
        << "   a = " << a
        << "   i_Ekman = " << i_Ekman 
        << endl << endl;
*/
            if(j <= (jm-1)/2){
                if((w_wind.y[j][k] <= 0.)&&(v_wind.y[j][k] >= 0.)){
                    if((alfa >= 0.)&&(alfa <= 45./pi180)){  // section I (j = 83 - 90) (7°N - 0°N)
                        angle = 180./pi180 - (alfa - Ekman_angle);
                    }
                    if((alfa > 45./pi180)&&(alfa <= 90./pi180)){  // section II (j = 76 - 83) (14°N - 7°N)
                        angle = 180./pi180 - (alfa - Ekman_angle);
                    }
                }
                if((w_wind.y[j][k] >= 0.)&&(v_wind.y[j][k] >= 0.)){
                    if((alfa <= 90./pi180)&&(alfa >= 45./pi180)){  // section III (j = 67 - 76) (23°N - 14°N)
                        angle = alfa + Ekman_angle;
                    }
                    if((alfa < 45./pi180)&&(alfa >= 0./pi180)){  // section IV (j = 49 - 67) (41°N - 23°N)
                        angle = alfa + Ekman_angle;
                    }
                }
                if((w_wind.y[j][k] >= 0.)&&(v_wind.y[j][k] <= 0.)){  // section V (j = 33 - 49) (57°N - 41°N)
                    if((alfa >= 0.)&&(alfa <= 45./pi180)){
                        angle = alfa + Ekman_angle;
                    }
                }
/*
// overcomes the singular behaviour around the equator
                if((j>=88)&&(j<=90)){
                    v_wind.y[j][k] = v_wind.y[87][k];
                    w_wind.y[j][k] = w_wind.y[87][k];
                }
*/
//              original laws in Robert H. Stewart, Introduction to Physical Oceanography, p. 138, eq. 9.11a/b
                for(int i = 0; i < im; i++){
//              original laws in Robert H. Stewart, Introduction to Physical Oceanography, p. 137, eq. 9.9a/b
                    a_z = a * rad.z[i] * L_hyd;
                    exp_a_z = exp(a_z);
                    sin_a_z = sin(angle + a_z);
//                    double angle_backturn = atan(angle + a_z) - angle;
//                    double angle_new = angle - angle_backturn;
                    sin_a_z = sin(angle + a_z);
                    cos_a_z = cos(angle + a_z);
//                    sin_a_z = sin(angle_new + a_z);
//                    cos_a_z = cos(angle_new + a_z);
                    v.x[i][j][k] = V_0 * exp_a_z * sin_a_z; 
                    w.x[i][j][k] = V_0 * exp_a_z * cos_a_z;
/*
    if((j == 75) &&(k == 180)) cout << "north" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   rad = " << rad.z[i] << endl
        << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
        << "   a = " << a << "   a_z = " << a_z << "   exp_a_z = " << exp_a_z 
        << "   sin_a_z = " << sin_a_z << "   cos_a_z = " << cos_a_z << endl
        << "   alfa = " << alfa
        << "   angle = " << angle << endl
        << "   alfa_deg = " << alfa * pi180
        << "   angle_deg = " << angle * pi180 << endl
        << "   Ekman = " << i_Ekman 
        << "   D_E = " << D_E << "   D_E_op = " << D_E_op << endl
        << "   T_yz = " << T_yz
        << "   U_10 = " << U_10
        << "   V_0 = " << V_0 << endl
        << "   v_wind = " << v_wind.y[j][k] 
        << "   w_wind = " << w_wind.y[j][k] << endl
        << "   v = " << v.x[i][j][k] 
        << "   w = " << w.x[i][j][k] << endl
        << "   f = " << f
        << "   Az = " << Az
        << "   a = " << a
        << endl;
*/
                }
            }
            if(j > (jm-1)/2){
                if((w_wind.y[j][k] <= 0.)&&(v_wind.y[j][k] <= 0.)){
                    if((alfa >= 0.)&&(alfa <= 45./pi180)){  // section I (j = 173 - 180) (7°S - 0°S)
                        angle = 180./pi180 - (alfa - Ekman_angle);
                    }
                    if((alfa > 45./pi180)&&(alfa <= 90./pi180)){  // section II (j = 166 - 173) (14°S - 7°S)
                        angle = 180./pi180 - (alfa - Ekman_angle);
                    }
                }
                if((w_wind.y[j][k] >= 0.)&&(v_wind.y[j][k] <= 0.)){
                    if((alfa > 45./pi180)&&(alfa <= 90./pi180)){  // section III (j = 157 - 166) (23°S - 14°S)
                        angle = alfa + Ekman_angle;
                    }
                    if((alfa >= 0.)&&(alfa <= 45./pi180)){  // section IV (j = 139 - 157) (41°S - 23°S)
                        angle = alfa + Ekman_angle;
                    }
                }
                if((w_wind.y[j][k] >= 0.)&&(v_wind.y[j][k] >= 0.)){
                    if((alfa < 45./pi180)&&(alfa >= 0./pi180)){  // section V (j = 123 - 139) (57°S - 41°S)
                        angle = alfa + Ekman_angle;
                    }
                }
/*
// overcomes the singular behaviour around the equator
                if((j>90)&&(j<=92)){
                    v_wind.y[j][k] = - v_wind.y[87][k];
                    w_wind.y[j][k] = w_wind.y[87][k];
                    }
*/
//              original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 138, eq. 9.11a/b
                for(int i = 0; i < im; i++){
//              original laws in Robert H. Stewart, Introduction to Physical Oceanography, p. 137, eq. 9.9a/b
                    a_z = a * rad.z[i] * L_hyd;
                    exp_a_z = exp(a_z);
                    sin_a_z = - sin(angle + a_z);
                    cos_a_z = cos(angle + a_z);
                    v.x[i][j][k] = V_0 * exp_a_z * sin_a_z; 
                    w.x[i][j][k] = V_0 * exp_a_z * cos_a_z;
/*
    if((j == 105) &&(k == 180)) cout << "south" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   rad = " << rad.z[i] << endl
        << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
        << "   a = " << a << "   a_z = " << a_z << "   exp_a_z = " << exp_a_z 
        << "   sin_a_z = " << sin_a_z << "   cos_a_z = " << cos_a_z << endl
        << "   alfa = " << alfa
        << "   angle = " << angle << endl
        << "   alfa_deg = " << alfa * pi180
        << "   angle_deg = " << angle * pi180 << endl
        << "   Ekman = " << i_Ekman 
        << "   D_E = " << D_E << "   D_E_op = " << D_E_op << endl
        << "   T_yz = " << T_yz
        << "   U_10 = " << U_10
        << "   V_0 = " << V_0 << endl
        << "   v_wind = " << v_wind.y[j][k] 
        << "   w_wind = " << w_wind.y[j][k] << endl
        << "   v = " << v.x[i][j][k] 
        << "   w = " << w.x[i][j][k] << endl
        << "   f = " << f
        << "   Az = " << Az
        << "   a = " << a
        << endl;
*/
                }
            }
        }
    }
// Ekman layer computation, variable EkmanPumping is equal to the vertical velocity u at the surface
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            u.x[0][j][k] = 0.;
            v.x[0][j][k] = 0.;
            w.x[0][j][k] = 0.;
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            aux_grad_v.z[0] = v.x[0][j][k];
            aux_grad_w.z[0] = w.x[0][j][k];
            for(int i = 1; i < im; i++){
                if(is_land(h, i, j, k)){
                    aux_grad_v.z[i] = 0.;
                    aux_grad_w.z[i] = 0.;
                }else{
                    aux_grad_v.z[i] = v.x[i][j][k];
                    aux_grad_w.z[i] = w.x[i][j][k];
                }
/*
                if((i_max-i_Ekman) % 2 == 0){
                    aux_v.x[im-1][j][k] = simpson(0, i, dr_rm, aux_grad_v);
                    aux_w.x[im-1][j][k] = simpson(0, i, dr_rm, aux_grad_w);
                }else  cout << "   i_max-i_Ekman  must be an even number to use the Simpson integration method" << endl;
*/
                double dr_rm = dr * L_hyd/(double)(im-1);
                aux_v.x[i][j][k] = trapezoidal(0, i, dr_rm, aux_grad_v);
                aux_w.x[i][j][k] = trapezoidal(0, i, dr_rm, aux_grad_w);
//                aux_v.x[i][j][k] = rectangular(0, i, dr_rm, aux_grad_v);
//                aux_w.x[i][j][k] = rectangular(0, i, dr_rm, aux_grad_w);
                if(is_land(h, i, j, k)){
                    aux_v.x[i][j][k] = 0.;
                    aux_w.x[i][j][k] = 0.;
                }
            }
        }
    }
     for(int j = 1; j < jm-1; j++){
        if(j <= 90)
            sinthe = sin(the.z[j]);
        if(j > 90){
            int j_rev = 90 - fabs(j - 90);
            sinthe = sin(the.z[j_rev]);
        }
        if(sinthe == 0.)  sinthe = 1.e-5;
        for(int k = 1; k < km-1; k++){
            for(int i = 1; i < im-1; i++){
                rm = rad.z[i] * L_hyd;
                rmsinthe = rm * sinthe;
                u.x[i][j][k] =
                    - ((aux_v.x[i][j+1][k] - aux_v.x[i][j-1][k])
                    /(2. * rm * dthe) 
                    + (aux_w.x[i][j][k+1] - aux_w.x[i][j][k-1])
                    /(2. * rmsinthe * dphi));
                u.x[i][j][k] = - u.x[i][j][k];  //rm is negative
/*
                u.x[i][j][k] = u.x[i-1][j][k] - dr * 
                    ((v.x[i-1][j+1][k] - v.x[i-1][j-1][k])
                    /(2. * rm * dthe) 
                    + (w.x[i-1][j][k+1] - w.x[i-1][j][k-1])
                    /(2. * rmsinthe * dphi));
*/
                if(is_land(h, i, j, k))  u.x[i][j][k] = 0.;
                if((j>=88)&&(j<=90)){
                    u.x[i][j][k] = u.x[i][87][k];
                }
                if((j>90)&&(j<=92)){
                    u.x[i][j][k] = u.x[i][87][k];
                }
                double residuum = (u.x[i][j][k] - u.x[i-1][j][k])/dr 
                    + ((v.x[i-1][j+1][k] - v.x[i-1][j-1][k])
                    /(2. * rm * dthe) 
                    + (w.x[i-1][j][k+1] - w.x[i-1][j][k-1])
                    /(2. * rmsinthe * dphi));
/*
    if((j == 75) &&(k == 180)) cout << "north pumping" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   aux_v = " << aux_v.x[i][j][k] 
        << "   aux_w = " << aux_w.x[i][j][k] << endl
        << "   u = " << u.x[i][j][k] 
        << "   v = " << v.x[i][j][k] 
        << "   w = " << w.x[i][j][k] << endl
        << "   rm = " << rm
        << "   dr = " << dr
        << "   residuum = " << residuum << endl;
    if((j == 105) &&(k == 180)) cout << "south pumping" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   aux_v = " << aux_v.x[i][j][k] 
        << "   aux_w = " << aux_w.x[i][j][k] << endl
        << "   u = " << u.x[i][j][k] 
        << "   v = " << v.x[i][j][k] 
        << "   w = " << w.x[i][j][k] << endl
        << "   rm = " << rm
        << "   dr = " << dr
        << "   residuum = " << residuum << endl;
*/
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            u.x[im-1][j][k] = c43 * u.x[im-2][j][k] - c13 * u.x[im-3][j][k];
//            u.x[im-1][j][k] = 0.;
            for(int i = 0; i < im; i++){
                    u.x[i][j][k] = u.x[i][j][k]/u_0;
                    v.x[i][j][k] = v.x[i][j][k]/u_0;
                    w.x[i][j][k] = w.x[i][j][k]/u_0;
                v.x[i][90][k] = 0.;
                if(is_land(h, i, j, k)){
                    u.x[i][j][k] = 0.;
                    v.x[i][j][k] = 0.;
                    w.x[i][j][k] = 0.;
                }
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::init_bathymetry(const string &bathymetry_file){
//  default adjustment, h must be 0 everywhere
    h.initArray(im, jm, km, 0.);
    ifstream ifile(bathymetry_file);
    if (!ifile.is_open()){
        cerr << "ERROR: could not open bathymetry file " << bathymetry_file << "\n";
        abort();
    }
    double lon, lat, depth;
    int j, k;
    for(j = 0; j < jm && !ifile.eof(); j++){
        for(k = 0; k < km && !ifile.eof(); k++){
            depth = 999; // in case the height is NaN
            ifile >> lon >> lat >> depth;
            if( depth > 0.){
                depth = 0;
            }
            int i_boden = (im-1) + int(floor(depth/L_hyd * (im-1)));
            i_bathymetry[j][k] = i_boden;
            Bathymetry.y[j][k] = - depth;
            for(int i = 0; i <= i_boden; i++){
                h.x[i][j][k] = 1.;
            }
            if(ifile.fail()){
                ifile.clear();
                std::string tmp;
                std::getline(ifile, tmp);
                logger() << "bad data in topography at: " << lon << " " << lat << " " << tmp << std::endl;
            }
        }
    }
    if(j != jm || k != km ){
        std::cerr << "wrong topography file size! aborting..."<<std::endl;
        abort();
    }
//  rewrite bathymetric data from -180° - 0° - +180° to 0°- 360°
    for(int j = 0; j < jm; j++){
        move_data(Bathymetry.y[j], km);
        for(int i = 0; i < im; i++){
            move_data(h.x[i][j], km);
            move_data(i_bathymetry[j], km);
        }
    }
//  reduction and smoothing of peaks and needles in the bathymetry
    for(int i = 0; i < im; i++){
            for(int j = 1; j < jm-1; j++){
        for(int k = 1; k < km-1; k++){
                if((is_land(h, i, j, k))
                     &&((is_water(h, i, j-1, k))
                     &&(is_water(h, i, j+1, k)))){
                    h.x[i][j][k] = 0.;
                }
                if((is_land(h, i, j, k))
                     &&((is_water(h, i, j, k-1))
                     &&(is_water(h, i, j, k+1)))){
                    h.x[i][j][k] = 0.;
                }
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::init_temperature(){
cout << endl << "      init_temperature" << endl;
    // Lenton_etal_COPSE_time_temp, constant paleo mean temperature, added to the surface initial temperature
    // difference between mean temperature (Ma) and mean temperature (previous Ma) == t_paleo_add
    std::map<float, float> pole_temp_map;  // Stein/Rüdiger/Parish linear pole temperature (Ma) distribution
    int Ma = *get_current_time();
    double t_equator = pow((rad_equator/sigma),0.25)/t_0;  // surface temperature at the equator t_equator = 1.0976 compares to 28.0°C = to 299.81 K
    double t_pole = pow((rad_pole/sigma),0.25)/t_0;  // surface temperature at the poles t_pole = 0.9436 compares to -15.4°C = to 250.25 K
    double t_eff = t_pole - t_equator;
    double get_pole_temperature(int Ma, const std::map<float, float> &pole_temp_map);
    load_map_from_file(pole_temperature_file, pole_temp_map); 
    float pole_temperature = (1.0 + get_pole_temperature(*get_current_time(), 
        pole_temp_map)/t_0) - t_pole;
//    if(pole_temperature <= t_pole)  pole_temperature = t_pole;
    double t_paleo_add = 0.0; 
    if(!is_first_time_slice()){
        if((NASATemperature != 0)&&(*get_current_time() > 0))  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - get_mean_temperature_from_curve(*get_previous_time());
        if(NASATemperature == 0)  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - t_average;
        t_paleo_add /= t_0; // non-dimensional 
    }

    cout.precision(3);
    const char* time_slice_comment = "      time slice of Paleo-OGCM:";
    const char* time_slice_number = " Ma = ";
    const char* time_slice_unit = " million years";
    cout << endl << setiosflags(ios::left) << setw(55) << setfill('.') 
        << time_slice_comment << resetiosflags(ios::left) << setw (6) 
        << fixed << setfill(' ') << time_slice_number << setw (3) << Ma 
        << setw(12) << time_slice_unit << endl << endl;
    const char* temperature_comment = "      temperature increase at paleo times: ";
    const char* temperature_gain = " t increase";
    const char* temperature_pole_comment = "      pole temperature increase: ";
    const char* temperature_gain_pole = " t pole increase";
    const char* temperature_modern = "      mean temperature at modern times: ";
    const char* temperature_paleo = "      mean temperature at paleo times: ";
    const char* temperature_average = " t modern";
    const char* temperature_average_pal = " t paleo";
    const char* temperature_unit =  "°C ";
    cout << endl << setiosflags(ios::left) 
        << setw(55) << setfill('.') 
        << temperature_comment << resetiosflags(ios::left) << setw(13) 
        << temperature_gain << " = " << setw(7) << setfill(' ') 
        << t_paleo_add * t_0 << setw(5) << temperature_unit << endl 

        << setiosflags(ios::left) 
        << setw(52) << setfill('.') 
        << temperature_pole_comment << resetiosflags(ios::left) << setw(13) 
        << temperature_gain_pole << " = " << setw(7) << setfill(' ') 
        << pole_temperature * t_0 << setw(5) << temperature_unit << endl 

        << setw(55) << setfill('.') 
        << setiosflags(ios::left) << temperature_modern 
        << resetiosflags(ios::left) << setw(13) << temperature_average 
        << " = " << setw(7) << setfill(' ') << t_average << setw(5) 
        << temperature_unit << endl 

        << setw(55) << setfill('.') 
        << setiosflags(ios::left) << temperature_paleo << resetiosflags(ios::left) 
        << setw(13) << temperature_average_pal << " = "  << setw(7) 
        << setfill(' ') << t_average + t_paleo_add * t_0 << setw(5) 
        << temperature_unit << endl;
    // temperatur distribution at a prescribed sun position
    // sun_position_lat = 60,    position of sun j = 120 means 30°S, j = 60 means 30°N
    // sun_position_lon = 180, position of sun k = 180 means 0° or 180° E (Greenwich, zero meridian)
    // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature (linear equation + parabola)
/*
    if((*get_current_time() > 0) && (sun == 1)){
        double j_par = sun_position_lat; // position of maximum temperature, sun position
        j_par = j_par + declination; // angle of sun axis, declination = 23,4°
        double j_pol = jm-1;
        double j_par_f = (double)j_par;
        double j_pol_f = (double)j_pol;
        double aa = (t_equator - t_pole)/(((j_par_f * j_par_f) 
            - (j_pol_f * j_pol_f)) - 2. * j_par_f * (j_par_f - j_pol_f));
        double bb = - 2. * aa * j_par_f;
        double cc = t_equator + aa * j_par_f * j_par_f;
        double j_d = sqrt ((cc - t_pole)/aa);
        double dd = 2. * aa * j_d + bb;
        double e = t_pole;
        // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature (linear equation + parabola)
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                double d_j = (double)j;
                if(d_j <= j_d){
                    t.x[im-1][j][k] = dd * d_j + e + t_paleo_add;
                }
                if(d_j > j_d){
                    t.x[im-1][j][k] = aa * d_j * d_j + bb * d_j 
                        + cc + t_paleo_add;
                }
            }
        }
        // longitudinally variable temperature distribution from west to east in parabolic form
        // communicates the impression of local sun radiation on the southern hemisphere
        double k_par = sun_position_lon;  // position of the sun at constant longitude
        double k_pol = km - 1;
        double t_360 = (t_0 + 5.)/t_0;
        for(int j = 0; j < jm; j++){
            double jm_temp_asym = t.x[im-1][j][20];//transfer of zonal constant temperature into aa 1D-temperature field
            for(int k = 0; k < km; k++){
                double k_par_f = (double)k_par;
                double k_pol_f = (double)k_pol;
                double d_k = (double) k;
                aa = (jm_temp_asym - t_360)/(((k_par_f * k_par_f) 
                    - (k_pol_f * k_pol_f)) - 2. * k_par_f * 
                    (k_par_f - k_pol_f));
                bb = - 2. * aa * k_par_f;
                cc = jm_temp_asym + aa * k_par_f * k_par_f;
                t.x[im-1][j][k] = aa * d_k * d_k + bb * d_k + cc;
            }
        }
    }// temperatur distribution at aa prescribed sun position
*/
    // pole temperature adjustment, combination of linear time dependent functions 
    // Stein/Rüdiger/Parish locally constant pole temperature
    // difference between pole temperature (Ma) and pole temperature (previous Ma)
    double d_j_half = (double)(jm-1)/2.0;
    float t_pole_diff_ocean = 0.0;
    float t_pole_diff_land = 0.0;
    //the t_pole_diff_ocean should be the difference between this time slice and the previous one
    if(!is_first_time_slice()){
        t_pole_diff_ocean = get_pole_temperature(*get_current_time(), pole_temp_map)  // no differences between t_pole_diff_ocean and t_pole_diff_land
            - get_pole_temperature(*get_previous_time(), pole_temp_map);
        t_pole_diff_land = get_pole_temperature(*get_current_time(), pole_temp_map)
            - get_pole_temperature(*get_previous_time(), pole_temp_map);
    }
    // in °C, constant local pole temperature as function of Ma for hothouse climates 
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double d_j = (double)j;
            if(NASATemperature == 0){  // parabolic ocean surface temperature assumed
                t.x[im-1][j][k] = t_eff 
                    * parabola((double)d_j/(double)d_j_half) + t_pole;
                t.x[im-1][j][k] += t_paleo_add + m_model->t_land
                    + pole_temperature * fabs(parabola((double)d_j
                    /(double)d_j_half) + 1.0);
                if(t.x[im-1][j][k] <= t_pole_salt)  t.x[im-1][j][k] = t_pole_salt;
                if(is_land(h, im-1, j, k)){
                    t.x[im-1][j][k] += m_model->t_land;
                }
/*
    cout.precision(9);
    if(k == 180) 
        cout << j << "     " << pole_temperature * t_0
            << "      " << pole_temperature * fabs(parabola((double)d_j
                /(double)d_j_half) + 1.0) * t_0
            << "     " << t_pole * t_0 - t_0
            << "     " << t.x[im-1][j][k] * t_0 - t_0
            << endl;
*/
            }else{  // if(NASATemperature == 1) ocean surface temperature based on NASA temperature distribution
                // transported for later time slices Ma by use_earthbyte_reconstruction
                if(is_land (h, im-1, j, k)){  // on land a parabolic distribution assumed, no NASA based data transportable
                    if(*get_current_time() > 0){
//                        t.x[im-1][j][k] = (t_eff * parabola(d_j/d_j_half) 
//                            + pole_temperature) + t_paleo_add + m_model->t_land;
                        t.x[im-1][j][k] += t_paleo_add + m_model->t_land
                            + t_pole_diff_land * fabs(parabola((double)d_j
                            /(double)d_j_half) + 1.0)/t_0;
                        if(t.x[im-1][j][k] <= t_pole_salt)  t.x[im-1][j][k] = t_pole_salt;
                        // land surface temperature increased by mean t_paleo_add
                        // and by a zonally equatorwards decreasing temperature difference is added
                        // Stein/Rüdiger/Parish pole temperature decreasing equator wards
                    }
                    if(*get_current_time() == 0)
                        t.x[im-1][j][k] = temperature_NASA.y[j][k];  // initial temperature by NASA for Ma=0
                }else{ // if the location is ocean
                    if(*get_current_time() > 0){        
                        // ocean surface temperature increased by mean t_paleo_add
                        // and by a zonally equatorwards decreasing temperature difference is added
                        // Stein/Rüdiger/Parish pole temperature decreasing equator wards
                        t.x[im-1][j][k] += t_paleo_add 
                            + t_pole_diff_ocean * fabs(parabola((double)d_j 
                           /(double)d_j_half) + 1.0)/t_0;
                        if(t.x[im-1][j][k] <= t_pole_salt)  t.x[im-1][j][k] = t_pole_salt;
                    }
                    if(*get_current_time() == 0)
                        t.x[im-1][j][k] = temperature_NASA.y[j][k];  // initial temperature by NASA for Ma=0
                }
            }// else(NASATemperature == 1)
        }// for j
    }// for k
    int i_beg = 0; // 200m depth  
    double tm_tbeg = 0.0;
    double d_i_max =(double)i_max;
    double d_i_beg =(double)i_beg;
    double d_i = 0.0;
// distribution of t and c with increasing depth
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            tm_tbeg = (t.x[im-1][j][k] - 1.) 
                /(d_i_max * d_i_max - d_i_beg * d_i_beg);
            for(int i = i_beg; i < im-1; i++){
                    d_i = (double)i;
                    t.x[i][j][k] = 1. + tm_tbeg 
                        * (d_i * d_i - d_i_beg * d_i_beg);// parabolic approach
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_sea_mount = i_bathymetry[j][k];
            for(int i = 0; i < im; i++){
                if((is_land(h, i, j, k))&&(i < i_sea_mount)){
                    t.x[i][j][k] = 1.;
                }
                if((is_land(h, 0, j, k))&&(i == 0)){
                    t.x[0][j][k] = t.x[i_sea_mount][j][k];
                }
            }
        }
    }
cout << "      init_temperature ended" << endl;
}
/*
*
*/
void cHydrosphereModel::init_salinity(){
cout << endl << "      init_salinity" << endl;
    // Lenton_etal_COPSE_time_temp, constant paleo mean temperature, added to the surface initial temperature
    // difference between mean temperature (Ma) and mean temperature (previous Ma) == t_paleo_add
    int Ma = *get_current_time();
    double get_pole_temperature(int Ma, const std::map<float, float> &pole_temp_map);
    double t_paleo_add = 0; 
    if(!is_first_time_slice()){
        if((NASATemperature != 0)&&(*get_current_time() > 0))  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - get_mean_temperature_from_curve(*get_previous_time());
        if(NASATemperature == 0)  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - t_average;
        t_paleo_add /= t_0; // non-dimensional 
    }

    double c_average = (t_average + 346.)/10.;// in psu, relation taken from "Ocean Circulation, The Open University"
    double c_paleo = (t_average + t_paleo_add * t_0 + 346.)/10.;// in psu
    c_paleo = c_paleo - c_average;
    if(Ma == 0)  c_paleo = 0.;
    cout.precision(3);
    string temperature_comment = "      temperature increase at paleo times: ";
    string temperature_gain = " t increase";
    string temperature_unit =  "°C ";
    string salinity_comment = "      salinity increase at paleo times: ";
    string salinity_gain = " salinity increase";
    string salinity_modern = "      mean salinity at modern times: ";
    string salinity_paleo = "      mean salinity at paleo times: ";
    string salinity_average = " salinity modern";
    string salinity_average_pal = " salinity paleo";
    string salinity_unit =  "psu ";
    cout << endl << setiosflags(ios::left) << setw(55) << setfill('.') 
        << temperature_comment << resetiosflags(ios::left) << setw(13) 
        << temperature_gain << " = " << setw(7) << setfill(' ') 
        << t_paleo_add * t_0 << setw(5) << temperature_unit << endl
        << setiosflags (ios::left) << setw (50) << setfill ('.')
        << salinity_comment << resetiosflags(ios::left) << setw(12) 
        << salinity_gain << " = " << setw(7) << setfill(' ') << c_paleo 
        << setw(5) << salinity_unit << endl << setw(52) << setfill('.') 
        << setiosflags(ios::left) << salinity_modern << resetiosflags(ios::left) 
        << setw(13) << salinity_average  << " = " << setw(7) << setfill(' ') 
        << c_average << setw(5) << salinity_unit << endl << setw(53) 
        << setfill('.') << setiosflags(ios::left) << salinity_paleo 
        << resetiosflags(ios::left) << setw(13) << salinity_average_pal 
        << " = " << setw(7) << setfill(' ') << c_average + c_paleo 
        << setw(5) << salinity_unit << endl;

    double C_p = 999.83;
    double beta_p = 0.808;
    double r_0_saltwater = 1027.0;  // in kg/m³

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            if(is_water(h, im-1, j, k)){
                double t_Celsius = t.x[im-1][j][k] * t_0 - t_0;
//                if(t_Celsius <= t_pole * t_0 - t_0)
//                    t_Celsius = t_pole * t_0 - t_0;
//                c.x[im-1][j][k] = ((t_Celsius + 346.)/10.)/c_0;
//                if(c.x[im-1][j][k] <= 1.) c.x[im-1][j][k] = 1.;
                double alfa_t_p = 0.0708 *(1.0 + 0.068 * t_Celsius);
                double gamma_t_p = 0.003 *(1.0 - 0.012 * t_Celsius);
                c.x[im-1][j][k] = ((r_0_saltwater - C_p  //in kg/m³ approximation by Gill (saltwater.pdf)
                    + (alfa_t_p + 35.0 * gamma_t_p)       // salt water density assumed as constant
                    * t_Celsius)/(beta_p + gamma_t_p * t_Celsius))/c_0;
                if(c.x[im-1][j][k] <= 1.) c.x[im-1][j][k] = 1.0;


            }else{
                c.x[im-1][j][k] = 1.0;
            }
        }
    }
    int i_beg = 0; // 200m depth  
//    double tm_tbeg = 0.0;
    double cm_cbeg = 0.0;
    double d_i_max =(double)i_max;
    double d_i_beg =(double)i_beg;
    double d_i = 0.0;
// distribution of t and c with increasing depth
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
//            tm_tbeg = (t.x[im-1][j][k] - 1.0) 
//                /(d_i_max * d_i_max - d_i_beg * d_i_beg);
//            if(t.x[im-1][j][k] <= t_pole)
//                t.x[im-1][j][k] = t_pole;
            cm_cbeg = (c.x[im-1][j][k] - 1.0) 
                /(d_i_max * d_i_max - d_i_beg * d_i_beg);
            for(int i = i_beg; i < im-1; i++){
                    d_i = (double)i;
//                    t.x[i][j][k] = 1. + tm_tbeg 
//                        * (d_i * d_i - d_i_beg * d_i_beg);// parabolic approach
//                    if(t.x[i][j][k] <= t_pole)
//                        t.x[i][j][k] = t_pole;
                    c.x[i][j][k] = 1. + cm_cbeg 
                        * (d_i * d_i - d_i_beg * d_i_beg);// parabolic approach
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_sea_mount = i_bathymetry[j][k];
            for(int i = 0; i < im; i++){
                if((is_land(h, i, j, k))&&(i < i_sea_mount)){
                    c.x[i][j][k] = 1.0;
                }
                if((is_land(h, 0, j, k))&&(i == 0)){
                    c.x[0][j][k] = c.x[i_sea_mount][j][k];
                }
            }
        }
    }
cout << "      init_salinity ended" << endl;
}
/*
*
*/
void cHydrosphereModel::init_dynamic_pressure(){
cout << endl << "      init_dynamic_pressure" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                p_dyn.x[i][j][k] = 0.5 * r_salt_water.x[i][j][k] 
                    * sqrt(((u.x[i][j][k] * u.x[i][j][k]) 
                    + (v.x[i][j][k] * v.x[i][j][k]) 
                    + (w.x[i][j][k] * w.x[i][j][k]))/3.0) * u_0;
                if(is_land(h, i, j, k))
                    p_dyn.x[i][j][k] = 0.0;
            }
        }
    }
cout << "      init_dynamic_pressure ended" << endl;
}
/*
*
*/
void cHydrosphereModel::BC_SolidGround(){
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                if(is_land(h, i, j, k)){
                    p_dyn.x[i][j][k] = 0.0;
                    t.x[i][j][k] = 1.0;
                    c.x[i][j][k] = 1.0;
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;
                    r_water.x[i][j][k] = r_0_water;
                    r_salt_water.x[i][j][k] = r_0_water;
                }
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::IC_t_WestEastCoast(){
// initial conditions for temperature at 
// the sea surface close to west coasts
// transition between coast flows and open sea flows included
// upwelling along west coasts causes a temperature reduction
    double t_neg = 0.;
    int k_mid_temp = 15; //  maximun extension of temperature management from the west coasts
//    int k_mid_temp = 10; //  maximun extension of temperature management from the west coasts
    int k_temp = 0; //  maximun extension of temperature management from the west coasts
    double temp_inter = 0.0;
    double temp_red = 0.95;
    for(int j = 0; j < jm; j++){
        for(int k = 1; k < km-1; k++){
            if((is_air(h, im-1, j, k-1))&&(is_land(h, im-1, j, k))){
                while(k >= 1){
                    if(is_air(h, im-1, j, k))  break;
                    t_neg = temp_red * t.x[im-1][j][k];
                    if(k <= k_mid_temp) k_temp = k_mid_temp - k;
                    else k_temp = k_mid_temp;
                    for(int l = 0; l <= k_temp; l++){
                        temp_inter = t.x[im-1][j][k-l];
                        t.x[im-1][j][k-l] = (t.x[im-1][j][k-k_mid_temp] - t_neg) 
                           /(double)k_mid_temp * (double)l + t_neg; // decreasing temperature leaving west coasts
                            // upwelling along west coasts lowers temperature
                        if(is_land(h, im-1, j, k-l))  t.x[im-1][j][k-l] = temp_inter;
                        if(t.x[im-1][j][k-l] <= t_pole_salt)  t.x[im-1][j][k-l] = t_pole_salt;
                    }
                    break;
                } // end while
            } // end if
        } // end k
    } // end j
    int j_var_north = 86;
    int j_var_south = 94;
//    int k_mid_eq = 120; // maximun extension from west coasts for equatorial lowered initial temperatures
    int k_mid_eq = 100; // maximun extension from west coasts for equatorial lowered initial temperatures
    for(int j = j_var_north; j <= j_var_south; j++){  // zonal extension of equatorial lowered initial temperatures
        for(int k = 1; k < km-1; k++){
            if((is_air(h, im-1, j, k-1))&&(is_land(h, im-1, j, k))){
                k_temp = k_mid_eq;
                t_neg = temp_red * t.x[im-1][j][k];  // reduction of the original parabolic temperature approach due to upwelling with lower temperatures
                for(int l = 0; l <= k_temp; l++){
                    temp_inter = t.x[im-1][j][k-l];
                    t.x[im-1][j][k-l] = (t.x[im-1][j][k-k_mid_eq] - t_neg) 
                       /(double)k_mid_eq * (double)l + t_neg; // decreasing temperature leaving west coasts
                        // upwelling regions along equatorial currents beginning on west coasts lower the local temperature
                    if(is_land(h, im-1, j, k-l))  t.x[im-1][j][k-l] = temp_inter;
                    if(t.x[im-1][j][k-l] <= t_pole_salt)  t.x[im-1][j][k-l] = t_pole_salt;
                }
            } // end if
        } // end k
    } // end j
}
/*
*
*/
void cHydrosphereModel::IC_u_WestEastCoast(){
// initial conditions for the u velocity component at the sea surface close to east or west coasts
    int i_beg = 20;                                                     // == 100m depth
    int i_half = i_beg + 10;                                            // location of u-max
    int j_half =(jm-1)/2;
    int i_max = im-1;
    double d_i_half =(double)i_half;
// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included
// northern hemisphere: east coast
//    k_grad = 6;                                                       // extension of velocity change
    int k_grad = 20;                                                    // extension of velocity change
    int k_water = 0;                                                    // on water closest to coast
    int k_sequel = 1;                                                   // on solid ground
    int m = 0;
    double d_i = 0.;
    for(int j = 0; j <= j_half; j++){                                   // outer loop: latitude
        for(int k = 0; k < km; k++){                                    // inner loop: longitude
            if(is_land(h, i_half, j, k)) k_sequel = 0;                  // if solid ground: k_sequel = 0
            if((is_water(h, i_half, j, k))&&(k_sequel == 0)) 
                k_water = 0;                                            // if water and and k_sequel = 0 then is water closest to coast
            else k_water = 1;                                           // somewhere on water
            if((is_water(h, i_half, j, k))&&(k_water == 0)){            // if water is closest to coast, change of velocity components begins
                for(int l = 0; l < k_grad; l++){                        // extension of change, sign change in v-velocity and distribution of u-velocity with depth
                    if(k+l > km-1) break;
                    for(int i = i_beg; i <= i_half; i++){               // loop in radial direction, extension for u-velocity component, downwelling here
                        m = i_max -(i - i_beg);
                        d_i =(double)i;
                        u.x[i][j][k+l] = - d_i/d_i_half/(double)(l+1);  // increase with depth, decrease with distance from coast
                        u.x[m][j][k+l] = - d_i/d_i_half/(double)(l+1);  // decrease with depth, decrease with distance from coast
                    }
                }
                k_sequel = 1;                                           // looking for another east coast
            }
        }                                                               // end of longitudinal loop
        k_water = 0;                                                    // starting at another latitude
    }                                                                   // end of latitudinal loop
// southern hemisphere: east coast
    k_water = 0;
    k_sequel = 1;
    for(int j = j_half+1; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_land(h, i_half, j, k)) k_sequel = 0;
            if((is_water(h, i_half, j, k))&&(k_sequel == 0)) 
                k_water = 0;
            else k_water = 1;
            if((is_water(h, i_half, j, k))&&(k_water == 0)){
                for(int l = 0; l < k_grad; l++){
                    if(k+l > km-1) break; 
                    for(int i = i_beg; i <= i_half; i++){
                        m = i_max -(i - i_beg);
                        d_i =(double)i;
                        u.x[i][j][k+l] = - d_i/d_i_half/(double)(l+1);  // increase with depth, decrease with distance from coast
                        u.x[m][j][k+l] = - d_i/d_i_half/(double)(l+1);  // decrease with depth, decrease with distance from coast
                    }
                }
                k_sequel = 1;
            }
        }
        k_water = 0;
    }
// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included
// northern hemisphere: west coast
//    k_grad = 6;                                                       // extension of velocity change
    k_water = 0;                                                        // somewhere on water
    int flip = 0;                                                           // somewhere on water
    for(int j = 0; j <= j_half; j++){                                        // outer loop: latitude
        for(int k = 0; k < km; k++){                                    // inner loop: longitude
            if(is_water(h, i_half, j, k)){                              // if somewhere on water
                k_water = 0;                                            // somewhere on water: k_water = 0
                flip = 0;                                               // somewhere on water: flip = 0
            }
            else k_water = 1;                                           // first time on land
            if((flip == 0)&&(k_water == 1)){                           // on water closest to land
                for(int l = k; l >(k - k_grad + 1); l--){               // backward extention of velocity change: nothing changes
                    if(l < 0) break;
                    for(int i = i_beg; i <= i_half; i++){               // loop in radial direction, extension for u -velocity component, downwelling here
                        m = i_max -(i - i_beg);
                        d_i =(double)i;
                        u.x[i][j][l] = + d_i/d_i_half/((double)(k - l + 1));  // increase with depth, decrease with distance from coast
                        u.x[m][j][l] = + d_i/d_i_half/((double)(k - l + 1));  // decrease with depth, decrease with distance from coast
                    }
                }
                flip = 1;
            }
        }
        flip = 0;
    }
// southern hemisphere: west coast
    k_water = 0;
    flip = 0;
    for(int j = j_half+1; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_water(h, i_half, j, k)){
                k_water = 0;
                flip = 0;
            }
            else k_water = 1;
            if((flip == 0)&&(k_water == 1)){
                for(int l = k; l >(k - k_grad + 1); l--){
                    if(l < 0) break;    
                    for(int i = i_beg; i <= i_half; i++){
                        m = i_max -(i - i_beg);
                        d_i =(double)i;
                        u.x[i][j][l] = + d_i/d_i_half/((double)(k - l + 1));
                        u.x[m][j][l] = + d_i/d_i_half/((double)(k - l + 1));
                    }
                }
                flip = 1;
            }
        }
        flip = 0;
    }
    for(int i = 0; i < i_beg; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                u.x[i][j][k] = 0.;
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::IC_Equatorial_Currents(){
// currents along the equator
// equatorial undercurrent - Cromwell flow, EUC
// equatorial intermediate current, EIC
// nothern and southern equatorial subsurface counter-currents, NSCC und SSCC
// nothern and southern equatorial counter-currents, NECC und SECC
    double IC_water = 1.;  // no dimension,(average velocity compares to   u_0 * IC_water = 0,25)
/*
// one grid step compares to a depth of 25 m for L_hyd = 1000m   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    i_EIC_u = i_beg; // 1000m depth  // westward equatorial intermediate current, max 0.1 m/s
    i_EIC_o = 28; // 0m depth
    i_SCC_u = 8; // 800m depth  // eastward subsurface counter current, max 0.05 m/s
    i_SCC_o = 28; // 300m depth
    i_ECC_u = 28; // 400m depth  // eastward equatorial counter current, max 0.2 m/s
    i_ECC_o = im; // 0m depth
    i_EUC_u = 32; // 200m depth  // eastward equatorial under current(Cromwell current), max 0.8 m/s
    i_EUC_o = 36; // 100m depth
*/
// one grid step compares to a depth of 5 m for L_hyd = 200m   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    int i_beg = 0; // 200m depth  
    int i_EIC_u = i_beg; // 200m depth  // westward equatorial intermediate current, max 0.1 m/s
    int i_EIC_o = 0; // 0m depth
    int i_SCC_u = 0; // 800m depth  // eastward subsurface counter current, max 0.05 m/s
    int i_SCC_o = 0; // 300m depth
    int i_ECC_u = 0; // 200m depth  // eastward equatorial counter current, max 0.2 m/s
    int i_ECC_o = im; // 0m depth
    int i_EUC_u = 0; // 200m depth  // eastward equatorial under current(Cromwell current), max 0.8 m/s
    int i_EUC_o = 20; // 100m depth
// equatorial currents and counter-currents
//  §§§§§§§§§§§§§§§§§§§§§§§§§   valid for all paleo-ocean constellations along the equator   §§§§§§§§§§§§§§§§§§§§§§§§§§
    int j_half =(jm -1)/2;
    int k_beg = 0;
    int k_end = 0;
// extention of land and ocean areas
    for(int k = 1; k < km-1; k++){
        if(k < km){
            if((is_water(h, im-1, j_half, k)) 
                &&(is_land(h, im-1, j_half, k+1))){
                k_end = k;
            }
            if((is_land(h, im-1, j_half, k)) 
                &&(is_water(h, im-1, j_half, k+1))){
                k_beg = k;
            }
            if(((is_water(h, im-1, j_half, k)) 
                &&(is_water(h, im-1, j_half, k+1))) 
                &&(k == km-2)){
                k_end = k;
            }
// equatorial northern counter-current(NECC, i=im-1 until i=im compares to 0 until -200m depth)
// equatorial northern counter-current(from j=83 until j=88 compares to 3°N until 8°N)
            for(int i = i_ECC_u; i < i_ECC_o; i++){
                for(int j =83; j < 88; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = + .01 * IC_water 
                                *(double)(i - i_ECC_u) 
                               /(double)(i_ECC_o - i_ECC_u);
                        }
                    }
                }
            }
// equatorial southern counter-current(SECC, i=im-1 until i=im compares to 0 until -200m depth)
// equatorial southern counter-current(from j=93 until j=96 compares to 3°S until 6°S)
            for(int i = i_ECC_u; i < i_ECC_o; i++){
                for(int j = 93; j < 98; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = + .01 * IC_water 
                                *(double)(i - i_ECC_u) 
                               /(double)(i_ECC_o - i_ECC_u);
                        }
                    }
                }
            }
// equatorial undercurrent - Cromwell current(EUC, i=im-2 until i=im-1 compares to -100 until -200m depth)
// equatorial undercurrent - Cromwell current(from j=87 until j=93 compares to 3°N until 3°S)
            for(int i = i_EUC_u; i < i_EUC_o; i++){
                for(int j = 87; j < 94; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = + .4 * IC_water;
                        }
                    }
                }
            }
// equatorial intermediate current(EIC, i=im-4 until i=im-2 compares to -300 until -1000m depth)
// equatorial intermediate current(from j=88 until j=92 compares to 2°N until 2°S)
            for(int i = i_EIC_u; i < i_EIC_o; i++){
                for(int j = 88; j < 93; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = - .05 * IC_water 
                                *(double)(i - i_EIC_u) 
                               /(double)(i_EIC_o - i_EIC_u);
                        }
                    }
                }
            }
// equatorial northern and southern subsurface counter-current
// equatorial northern subsurface counter-current(NSCC, i=im-3 until i=im-2 compares to -300 until -800m depth)
// equatorial northern subsurface counter-current(from j=86 until j=87 compares to 3°N until 4°N)
            for(int i = i_SCC_u; i < i_SCC_o; i++){
                for(int j = 86; j < 88; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = .025 * IC_water 
                                *(double)(i - i_SCC_u) 
                               /(double)(i_SCC_o - i_SCC_u);
                        }
                    }
                }
            }
// equatorial southern subsurface counter-current(SSCC, i=im-3 until i=im-2 compares to -300 until -800m depth)
// equatorial southern subsurface counter-current(from j=93 until j=94 compares to 3°S until 4°S)
            for(int i = i_SCC_u; i < i_SCC_o; i++){
                for(int j = 93; j < 95; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = .025 * IC_water 
                                *(double)(i - i_SCC_u) 
                               /(double)(i_SCC_o - i_SCC_u);
                        }
                    }
                }
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::IC_CircumPolar_Current(){ 
// south polar sea
// antarctic circumpolar current(-5000m deep)(from j=147 until j=152 compares to 57°S until 62°S,
// from k=0 until k=km compares to 0° until 360°)
    int i_beg = 0;
    for(int i = i_beg; i < im; i++){
        for(int j = 147; j < 153; j++){
            for(int k = 0; k < km; k++){
                if(is_water(h, i, j, k)){
//                  c.x[i][j][k] = 1.;
                    w.x[i][j][k] = 2. * u_0;  // 0.5 m/s
                }
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::BC_Surface_Temperature_NASA
   (const string &Name_SurfaceTemperature_File){
    cout.precision(3);
    cout.setf(ios::fixed);
    ifstream Name_SurfaceTemperature_File_Read;
    Name_SurfaceTemperature_File_Read.open(Name_SurfaceTemperature_File);
    if(!Name_SurfaceTemperature_File_Read.is_open()){
        cerr << "ERROR: could not open SurfaceTemperature_File file at " << Name_SurfaceTemperature_File << "\n";
        abort();
    }
    int j = 0;
    int k = 0;
    double dummy_1, dummy_2, dummy_3;
    while((k < km) &&(!Name_SurfaceTemperature_File_Read.eof())){
        while(j < jm){
            Name_SurfaceTemperature_File_Read >> dummy_1;
            Name_SurfaceTemperature_File_Read >> dummy_2;
            Name_SurfaceTemperature_File_Read >> dummy_3;
            t.x[im-1][j][k] =(dummy_3 + 273.15)/273.15;
            j++;
        }
        j = 0;
        k++;
    }
    Name_SurfaceTemperature_File_Read.close();
    for(int j = 0; j < jm; j++){
        for(int k = 1; k < km-1; k++){
            if(k == 180) t.x[im-1][j][k] =(t.x[im-1][j][k + 1] 
                + t.x[im-1][j][k - 1]) * .5;
        }
    }
}
/*
*
*/
void cHydrosphereModel::BC_Surface_Salinity_NASA(const string &Name_SurfaceSalinity_File){
    // initial conditions for the salinity at the sea surface
    streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;
    cout.precision(3);
    cout.setf(ios::fixed);
    ifstream Name_SurfaceSalinity_File_Read;
    Name_SurfaceSalinity_File_Read.open(Name_SurfaceSalinity_File);
    if(!Name_SurfaceSalinity_File_Read.is_open()){
        cerr << "ERROR: could not open SurfaceSalinity_File file at " << Name_SurfaceSalinity_File << "\n";
        abort();
    }
    int j = 0;
    int k = 0;
    double dummy_1, dummy_2, dummy_3;
    while((k < km) &&(!Name_SurfaceSalinity_File_Read.eof())){
        while(j < jm){
            Name_SurfaceSalinity_File_Read >> dummy_1;
            Name_SurfaceSalinity_File_Read >> dummy_2;
            Name_SurfaceSalinity_File_Read >> dummy_3;
            if(dummy_3 < 0.) dummy_3 = 1.;
            else  c.x[im-1][j][k] = dummy_3/c_0;
            j++;
        }
        j = 0;
        k++;
    }
    Name_SurfaceSalinity_File_Read.close();
}
/*
*
*/
void cHydrosphereModel::load_temperature_curve(){
    load_map_from_file(temperature_curve_file, m_temperature_curve);
}
/*
*
*/
float cHydrosphereModel::get_mean_temperature_from_curve(float time) const{
    if(time<m_temperature_curve.begin()->first 
        || time>(--m_temperature_curve.end())->first){
        std::cout << "Input time out of range: " <<time<< std::endl;    
        return NAN;
    }
    if(m_temperature_curve.size()<2){
        std::cout << "No enough data in m_temperature_curve  map" << std::endl;
        return NAN;
    }
    map<float, float>::const_iterator upper=m_temperature_curve.begin(), 
        bottom=++m_temperature_curve.begin(); 
    for(map<float, float>::const_iterator it = m_temperature_curve.begin();
            it != m_temperature_curve.end(); ++it){
        if(time < it->first){
            bottom = it;
            break;
        }else{
            upper = it;
        }
    }
    //std::cout << upper->first << " " << bottom->first << std::endl;
    return upper->second + (time - upper->first) 
       /(bottom->first - upper->first) 
        * (bottom->second - upper->second);
}
