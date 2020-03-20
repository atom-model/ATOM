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
#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "Accuracy_Hyd.h"
#include "Utils.h"
#include "Config.h"
#include "tinyxml2.h"
#include "PythonStream.h"
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
//const float cHydrosphereModel::dt = 0.00001; // time step satisfies the CFL condition
const float cHydrosphereModel::dt = 0.0001; // time step satisfies the CFL condition
const float cHydrosphereModel::pi180 = 180./M_PI;  // pi180 = 57.3
const float cHydrosphereModel::the_degree = 1.;   // compares to 1° step size laterally
const float cHydrosphereModel::phi_degree = 1.;  // compares to 1° step size longitudinally
const float cHydrosphereModel::dthe = the_degree/pi180; // dthe = the_degree/pi180 = 1.0/57.3 = 0.01745, 180 * .01745 = 3.141
const float cHydrosphereModel::dphi = phi_degree/pi180; // dphi = phi_degree/pi180 = 1.0/57.3 = 0.01745, 360 * .01745 = 6.282

cHydrosphereModel::cHydrosphereModel():
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
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
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
    reset_arrays();
    m_current_time = m_time_list.insert(float(Ma)).first;
    mkdir(output_path.c_str(), 0777);
    Ma = int(round(*get_current_time()));
    cout.precision(6);
    cout.setf(ios::fixed);
    string Name_v_w_Transfer_File;
    stringstream ssName_v_w_Transfer_File;
    string Name_SurfaceTemperature_File = temperature_file;
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
    Atmosphere_TransferFile_read(bathymetry_name);
    cout << "***** time slice for the Oceanic Global Circulation Modell(OGCM) is:    Ma = " << Ma << " million years" 
        << endl << endl;
    cout << "***** bathymetry/topography given by the x-y-z data set:    " << bathymetry_name.c_str() << endl << endl;
    BC_SeaGround(bathymetry_path + "/" + bathymetry_name);
    if(use_NASA_velocity){
        read_IC(velocity_v_file, v.x[im-1], jm, km);
        read_IC(velocity_w_file, w.x[im-1], jm, km);    
        read_IC(Name_SurfaceTemperature_File, t.x[im-1], jm, km);
        read_IC(Name_SurfaceSalinity_File, c.x[im-1], jm, km);
    }
    iter_cnt_3d = -1;
    if(debug) save_data();
    iter_cnt_3d++;
    land_oceanFraction();
    EkmanSpiral();
    Value_Limitation_Hyd();
//    IC_u_WestEastCoast();
//    IC_Equatorial_Currents();
//    if(Ma <= 41)  IC_CircumPolar_Current(); // Drake passage closed 41 Ma ago
    BC_Temperature_Salinity();
//    goto Printout;
    PresStat_SaltWaterDens();
    store_intermediate_data_2D();
    store_intermediate_data_3D();
    run_2D_loop();
    cout << endl << endl;
    run_3D_loop();
    cout << endl << endl;
    if(debug) save_data();
/*
Printout:
//run_data_hyd(); 
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

    r_water.initArray(im, jm, km, r_0_water); // water density as function of pressure
    r_salt_water.initArray(im, jm, km, r_0_water); // salt water density as function of pressure and temperature
    BuoyancyForce_3D.initArray(im, jm, km, 0.); // 3D buoyancy force
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
    int j_longal = 75;
    paraview_vtk_longal(bathymetry_name, j_longal, iter_cnt-1);
    int k_zonal = 185;
    paraview_vtk_zonal(bathymetry_name, k_zonal, iter_cnt-1);
    int i_radial = 40;
    paraview_vtk_radial(bathymetry_name, i_radial, iter_cnt-1);
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
            if(velocity_iter % 2 == 0){
                PresStat_SaltWaterDens();
            }
            solveRungeKutta_3D_Hydrosphere(); 
            run_data_hyd();
            Value_Limitation_Hyd();
            print_min_max_hyd();
            store_intermediate_data_3D();
            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        } // end of velocity loop
        computePressure_3D();
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
    return;
}
/*
*
*/
void cHydrosphereModel::EkmanSpiral(){
    //Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
    //a further turning downwards until the end of the shear layer such that finally 90° of turning are reached
    float Ekman_angle = 45./pi180;
    // initial conditions for v and w velocity components at the sea surface
    // ocean surface velocity is about 3% of the wind velocity at the surface
    // u_0 for the hydrosphere is (0.03 * u_0) for the atmosphere
    // surface wind vector driving the Ekman spiral in the Ekman layer
    // northern and southern hemisphere
    int i_Ekman = 0;
    int i_max = im-1;
    double sinthe = 0.;
    double rm = 0.;
    double vel_magnit = 0.;
    double vel_mag = 0.;
    double i_Ekman_layer = 0.;
    double coeff = 0.;
    double gam_z = 0.;
    double exp_gam_z = 0.;
    double sin_gam_z = 0;
    double cos_gam_z = 0;
    double alfa = 0;
    float angle = 0;
    for(int j = 1; j < jm-1; j++){
        sinthe = sin(fabs((90. -(double)j) * M_PI/180.));
        if(sinthe == 0.)  sinthe = 1.e-5;
        for(int k = 1; k < km-1; k++){
            vel_magnit = sqrt(v.x[im-1][j][k] * v.x[im-1][j][k] 
                + w.x[im-1][j][k] * w.x[im-1][j][k]); // in m/s
            // dimensional surface wind velocity U_10 in m/s
            // original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 139, eq. 9.16
            i_Ekman_layer = 3. * 7.6/sqrt(sinthe) * vel_magnit; // dimensional surface wind velocity U_10 in m/s
            i_Ekman_layer = 1.5 * i_Ekman_layer;            // assumed depth at i_Ekman = 3 x pi/gam, velocity opposite direction plus pi to surface velocity
            coeff = i_Ekman_layer/L_hyd;
            if(coeff >= 1.) coeff = .623;
            i_Ekman = (im-1) * (1.-coeff);
/*
    if((j == 90) &&(k == 180)) cout << endl << "   coeff = " << coeff 
        << "   j = " << j << "   k = " << k 
        << "   sqrt(sinthe) = " << sqrt(sinthe) 
        << "   i_Ekman_layer = " << i_Ekman_layer 
        << "   u_0 = " << u_0 
        << "   i_Ekman = " << i_Ekman 
        << "   v = " << v.x[im-1][j][k] 
        << "   w = " << w.x[im-1][j][k] 
        << "   vel_magnit = " << vel_magnit 
        << endl << endl;
*/
            v.x[im-1][j][k] = v.x[im-1][j][k]/u_0;
            w.x[im-1][j][k] = w.x[im-1][j][k]/u_0;
// overcomes the singular behaviour around the equator
            if((j>=88)&&(j<=90)){
                v.x[im-1][j][k] = v.x[im-1][87][k];
                w.x[im-1][j][k] = w.x[im-1][87][k];
            }
            if((j>90)&&(j<=92)){
                v.x[im-1][j][k] = - v.x[im-1][87][k];
                w.x[im-1][j][k] = w.x[im-1][87][k];
            }
            v.x[im-1][90][k] = 0.;
            if(w.x[im-1][j][k] == 0.) w.x[im-1][j][k] = 1.e-8;
            alfa = atan(fabs(v.x[im-1][j][k]/w.x[im-1][j][k]));
            if(j <=(jm-1)/2){
                if((w.x[im-1][j][k] >= 0.)&&(v.x[im-1][j][k] >= 0.))
                    angle = alfa + Ekman_angle;
                if((w.x[im-1][j][k] <= 0.)&&(v.x[im-1][j][k] >= 0.))
                    angle = M_PI - alfa + Ekman_angle;
                if((w.x[im-1][j][k] >= 0.)&&(v.x[im-1][j][k] <= 0.))
                    angle = 2. * M_PI + alfa + Ekman_angle;
                if((w.x[im-1][j][k] <= 0.)&&(v.x[im-1][j][k] <= 0.))
                    angle = 2. * M_PI - alfa + Ekman_angle;
//              original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 138, eq. 9.11a/b
                vel_mag = - sqrt(v.x[im-1][j][k] * v.x[im-1][j][k] 
                          + w.x[im-1][j][k] * w.x[im-1][j][k]);

                for(int i = im-2; i >= i_Ekman; i--){
//              original laws in Robert H. Stewart, Introduction to Physical Oceanography, p. 137, eq. 9.9a/b
                    gam_z = - 3. * M_PI * ((double)i-(double)i_max) 
                        /((double)i_Ekman-(double)i_max);
                    exp_gam_z = exp(gam_z);
                    sin_gam_z = - sin(angle + gam_z);
                    cos_gam_z = cos(angle + gam_z);
                    v.x[i][j][k] = vel_mag * exp_gam_z * sin_gam_z; 
                    w.x[i][j][k] = - vel_mag * exp_gam_z * cos_gam_z;
/*
    if((j == 89) &&(k == 180)) cout << "north" << "   i = " << i 
        << "   j = " << j << "   k = " << k 
        << "   gam_z = " << gam_z << "   exp_gam_z = " << exp_gam_z 
        << "   sin_gam_z = " << sin_gam_z << "   cos_gam_z = " << cos_gam_z 
        << "   i_Ekman = " << i_Ekman 
        << "   i_Ekman_layer = " << i_Ekman_layer 
        << "   v = " << v.x[i][j][k] 
        << "   w = " << w.x[i][j][k] 
        << "   vel_mag = " << vel_mag 
        << "   vel_magnit = " << vel_magnit << endl;
*/
                }
            }else{
                if((w.x[im-1][j][k] >= 0.)&&(v.x[im-1][j][k] <= 0.))
                    angle = M_PI -(2. * M_PI - alfa - Ekman_angle);
                if((w.x[im-1][j][k] <= 0.)&&(v.x[im-1][j][k] <= 0.))
                    angle = M_PI -(M_PI + alfa - Ekman_angle);
                if((w.x[im-1][j][k] >= 0.)&&(v.x[im-1][j][k] >= 0.))
                    angle = M_PI -(alfa - Ekman_angle);
                if((w.x[im-1][j][k] <= 0.)&&(v.x[im-1][j][k] >= 0.))
                    angle = M_PI -(M_PI - alfa - Ekman_angle);
//              original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 138, eq. 9.11a/b
                vel_mag = - sqrt(v.x[im-1][j][k] * v.x[im-1][j][k] 
                          + w.x[im-1][j][k] * w.x[im-1][j][k]);

                for(int i = im-2; i >= i_Ekman; i--){
//              original laws in Robert H. Stewart, Introduction to Physical Oceanography, p. 137, eq. 9.9a/b
                    gam_z = - 3. * M_PI * ((double)i-(double)i_max) 
                        /((double)i_Ekman-(double)i_max);
                    exp_gam_z = exp(gam_z);
                    sin_gam_z = - sin(angle + gam_z);
                    cos_gam_z = cos(angle + gam_z);
                    v.x[i][j][k] = vel_mag * exp_gam_z * sin_gam_z; 
                    w.x[i][j][k] = vel_mag * exp_gam_z * cos_gam_z;
/*
    if((j == 91) &&(k == 180)) cout << "south" << "   i = " << i 
        << "   j = " << j << "   k = " << k 
        << "   gam_z = " << gam_z << "   exp_gam_z = " << exp_gam_z 
        << "   sin_gam_z = " << sin_gam_z << "   cos_gam_z = " << cos_gam_z 
        << "   i_Ekman = " << i_Ekman 
        << "   i_Ekman_layer = " << i_Ekman_layer 
        << "   v = " << v.x[i][j][k] 
        << "   w = " << w.x[i][j][k] 
        << "   vel_mag = " << vel_mag 
        << "   vel_magnit = " << vel_magnit << endl;
*/
                }
            }
        }
    }
// vertical velocity u calculated by the continuity equation
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            u.x[0][j][k] = 0.;
            for(int i = 0; i < im-1; i++){
                rm = rad.z[i];
                u.x[i+1][j][k] = u.x[i][j][k]-dr*((v.x[i][j+1][k] 
                    -v.x[i][j-1][k])/(2.*rm*dthe)+(w.x[i][j][k+1] 
                    -w.x[i][j][k-1])/(2.*rm*sinthe*dphi));
/*
                if((j>=88)&&(j<=90)){
                    u.x[i+1][j][k] = u.x[i+1][87][k];
                }
                if((j>90)&&(j<=92)){
                    u.x[i+1][j][k] = u.x[i+1][87][k];
                }
*/
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
//            u.x[im-1][j][k] = u.x[im-4][j][k] - 3. * u.x[im-3][j][k] + 3. * u.x[im-2][j][k];        // extrapolation
            u.x[im-1][j][k] = 0.;
            v.x[im-1][j][k] = v.x[im-4][j][k] - 3. * v.x[im-3][j][k] + 3. * v.x[im-2][j][k];        // extrapolation
            w.x[im-1][j][k] = w.x[im-4][j][k] - 3. * w.x[im-3][j][k] + 3. * w.x[im-2][j][k];        // extrapolation
            for(int i = 0; i <= im-1; i++){
                v.x[i][90][k] = 0.;
                if(is_land(h, i, j, k)){
                    u.x[i][j][k] = 0.;
                    v.x[i][j][k] = 0.;
                    w.x[i][j][k] = 0.;
                }
            }
        }
    }
//    BC_radius();
    BC_theta();
    BC_phi();
}
/*
*
*/
void cHydrosphereModel::BC_SeaGround(const string &bathymetry_file){
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
            int i_boden = (im-1) + int(floor(depth / L_hyd * (im - 1)));
            Bathymetry.y[j][k] = -depth;
            for(int i = 0; i <= i_boden; i++){
                h.x[i][j][k] = 1.;
            }
            if(ifile.fail()){
                ifile.clear();
                std::string tmp;
                std::getline(ifile, tmp);
                logger() << "bad data in topography at: " << lon << " " << lat << " " << tmp << std::endl;
            }
            //logger() << lon << " " << lat << " " << h.x[0][j][k] << std::endl;            
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
        }
    }
//  reduction and smoothing of peaks and needles in the bathymetry
    double h_center = 0.;
    for(int k = 2; k < km-2; k++){
        for(int j = 2; j < jm-2; j++){
            for(int i = 2; i <= im-1; i++){
                if((is_land(h, i, j, k)) && ((is_water(h, i, j+1, k)) 
                    && (is_water(h, i, j-1, k))) && ((is_water(h, i, j, k+1)) 
                    && (is_water(h, i, j, k-1)))){
                    h_center = h.x[i][j][k];
                    h.x[i][j][k] = 0.;
                 }
/*
                if((h_center == 1.) && ((is_water(h, i, j-2, k)) 
                    && (is_water(h, i, j, k+2)))){
                    h.x[i][j-1][k+1] = 0.;
                }
                if((h_center == 1.) && ((is_water(h, i, j-2, k)) 
                    && (is_water(h, i, j, k-2)))){ 
                    h.x[i][j-1][k-1] = 0.;
                }
                if((h_center == 1.) && ((is_water(h, i, j+2, k)) 
                    && (is_water(h, i, j, k-2)))){
                    h.x[i][j+1][k-1] = 0.;
                }
                if((h_center == 1.) && ((is_water(h, i, j+2, k)) 
                    && (is_water(h, i, j, k+2)))){ 
                    h.x[i][j+1][k+1] = 0.;
                }
*/
/*
                if((h_center == 1.) && ((is_land(h, i, j, k)) 
                    && (is_water(h, i, j, k+2)))){
                    h.x[i][j][k] = h.x[i][j][k+1] = 0.;
                }
                if((h_center == 1.) && ((is_land(h, i, j, k)) 
                    && (is_water(h, i, j+2, k)))){ 
                    h.x[i][j][k] = h.x[i][j+1][k] = 0.;
                }
*/
                if((i >= im-2) && (h_center == 1.) && (is_water(h, im-3, j, k))){
                    h.x[im-1][j][k] = 0.;
                }
/*
                if((h_center == 1.) && ((is_water(h, i, j+2, k)) 
                    && (is_water(h, i, j, k+2)))){ 
                    h.x[i][j+1][k+1] = 0.;
                }
*/
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::BC_SolidGround(){
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                if(is_land(h, i, j, k)){
                    p_dyn.x[i][j][k] = 0.;
                    t.x[i][j][k] = 1.;
                    c.x[i][j][k] = c_0;
                    u.x[i][j][k] = 0.;
                    v.x[i][j][k] = 0.;
                    w.x[i][j][k] = 0.;
                    r_water.x[i][j][k] = r_0_water;
                    r_salt_water.x[i][j][k] = r_0_water;
                }
            }
        }
    }
}

