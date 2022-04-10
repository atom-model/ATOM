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

#include "Utils.h"
#include "Config.h"
#include "AtomMath.h"
#include "cHydrosphereModel.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

#define dxdthe_a(X) \
    (h_d_j * (X->x[i][j+1][k] - X->x[i][j-1][k])/(2.0 * dthe))
#define dxdthe_b(X) \
    (h_d_j * (- 3.0 * X->x[i][j][k] + 4.0 * X->x[i][j+1][k] - X->x[i][j+2][k])/(2.0 * dthe))
#define dxdthe_c(X) \
    (- h_d_j * (- 3.0 * X->x[i][j][k] + 4.0 * X->x[i][j-1][k] - X->x[i][j-2][k])/(2.0 * dthe))
#define dxdphi_a(X) \
    (h_d_k * (X->x[i][j][k+1] - X->x[i][j][k-1])/(2.0 * dphi))
#define dxdphi_b(X) \
    (h_d_k * (- 3.0 * X->x[i][j][k] + 4.0 * X->x[i][j][k+1] - X->x[i][j][k+2])/(2.0 * dphi))
#define dxdphi_c(X) \
    (- h_d_k * (- 3.0 * X->x[i][j][k] + 4.0 * X->x[i][j][k-1] - X->x[i][j][k-2])/(2.0 * dphi))

extern std::vector<std::vector<double> > m_node_weights;

cHydrosphereModel* cHydrosphereModel::m_model = NULL;

const double cHydrosphereModel::pi180 = 180.0/M_PI;  // pi180 = 57.3

const double cHydrosphereModel::the_degree = 1.0;   // compares to 1° step size laterally
const double cHydrosphereModel::phi_degree = 1.0;  // compares to 1° step size longitudinally

const double cHydrosphereModel::dr = 0.025;  // 0.025 x 40 = 1.0 compares to 200m : 40 = 5m for 1 radial step
const double cHydrosphereModel::dt = 0.001; // time step satisfies the CFL condition
//const double cHydrosphereModel::dt = 0.0001; // time step satisfies the CFL condition

const double cHydrosphereModel::dthe = the_degree/pi180; // dthe = the_degree/pi180 = 1.0/57.3 = 0.01745, 180 * .01745 = 3.141
const double cHydrosphereModel::dphi = phi_degree/pi180; // dphi = phi_degree/pi180 = 1.0/57.3 = 0.01745, 360 * .01745 = 6.282
    
const double cHydrosphereModel::the0 = 0.0;             // North Pole
const double cHydrosphereModel::phi0 = 0.0;             // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cHydrosphereModel::r0 = 1.0; // non-dimensional

//const double cHydrosphereModel::residuum_ref_hyd = 1.0e-4;     // criterium to finish global iterations

cHydrosphereModel::cHydrosphereModel():
    i_bathymetry(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    has_printed_welcome_msg(false){
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.
    if(PythonStream::is_enable()){
        backup = std::cout.rdbuf();
        std::cout.rdbuf(&ps);
    }
    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);
    // set default configuration
    SetDefaultConfig();
    m_model = this;
    //  Coordinate system in form of a spherical shell
    //  rad for r-direction bottom to the surface of the earth, the for lateral and phi for longitudinal direction
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
            throw std::invalid_argument(std::string("unable to load config file:  ") 
                + filename);
        }
        XMLElement *atom = doc.FirstChildElement("atom"), *elem_common = NULL, 
            *elem_hydrosphere = NULL;
        if(!atom){
            throw std::invalid_argument(std::string
                ("Failed to find the 'atom' element in config file: ") 
                + filename);
        }else{
            elem_common = atom->FirstChildElement("common");
            if(!elem_common){
                throw std::invalid_argument(std::string
                    ("Failed to find the 'common' element in 'atom' element in config file: ") 
                    + filename);
            }
            elem_hydrosphere = atom->FirstChildElement("hydrosphere");
            if(!elem_hydrosphere){
                throw std::invalid_argument(std::string(
                    "Failed to find the 'hydrosphere' element in 'atom' element in config file: ") 
                    + filename);
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
    int start = RunStart(hy);

    m_current_time = m_time_list.insert(int(Ma)).first;

    if(debug){
        feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO); //not platform independent, bad, very bad, I know
    }
    if(debug) save_data();

    use_presets(use_earthbyte_reconstruction, 
        use_NASA_temperature, use_NASA_velocity);

    iter_cnt = 1;
    iter_cnt_3d = -1;
    iter_cnt_3d++;
    velocity_iter = 0;
    pressure_iter = 0;

    reset_arrays();

    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);

    read_Hydrosphere_Surface_Data(Ma);

    land_oceanFraction(im-1, jm, km, h);

    init_EkmanSpiral();

//    goto Printout;

    init_dynamic_pressure();

    fft_gaussian_filter_3d(u,1);
    fft_gaussian_filter_3d(v,1);
    fft_gaussian_filter_3d(w,1);

//    IC_u_WestEastCoast();
//    IC_Equatorial_Currents();
//    if(Ma <= 41)  IC_CircumPolar_Current(); // Drake passage closed 41 Ma ago
    init_temperature(Ma);

//    goto Printout;
    
    if(!use_NASA_temperature) 
        IC_t_WestEastCoast();

    init_salinity();

    fft_gaussian_filter_3d(t,1);
//    fft_gaussian_filter_3d(t, 3, direction_k);
    fft_gaussian_filter_3d(c,1);

    store_intermediate_data_2D(1.0);
    store_intermediate_data_3D(1.0);

    PresStat_SaltWaterDens();
    SalinityEvaporation();

    store_intermediate_data_2D(1.0);
    store_intermediate_data_3D(1.0);
    cout << endl << endl;
    run_3D_loop(); // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
    cout << endl << endl;
/*
    Printout:
    run_data_hyd(); 
    print_min_max_hyd();
    write_file(bathymetry_name, output_path, true);
*/
    iter_cnt_3d++;
    save_data();
    if(debug){
        fedisableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO); //not platform independent(bad, very bad, I know)
    }
    RunEnd(hy, Ma, start);
    print_final_remarks();
    return;
}
/*
*
*/
void cHydrosphereModel::reset_arrays(){
    cout << endl << "      OGCM: reset_arrays" << endl;
    rad.initArray_1D(im, 1.0); // radial coordinate direction
    the.initArray_1D(jm, 2.0); // lateral coordinate direction
    phi.initArray_1D(km, 3.0); // longitudinal coordinate direction
    aux_grad_v.initArray_1D(im, 4.0); // auxilliar array
    aux_grad_w.initArray_1D(im, 5.0); // auxilliar array

    Bathymetry.initArray_2D(jm, km, 0.0); // Bathymetry in m
    Upwelling.initArray_2D(jm, km, 0.0); // upwelling
    Downwelling.initArray_2D(jm, km, 0.0); // downwelling
    EkmanPumping.initArray_2D(jm, km, 0.0); // 2D Ekman pumping vertical velocity summed up in a vertical column
    SaltFinger.initArray_2D(jm, km, 0.0);      // salt bulge of higher density
    SaltDiffusion.initArray_2D(jm, km, 0.0);   // salt bulge of lower density
    Salt_total.initArray_2D(jm, km, 0.0);     // rate of salt summed up in a vertical column
    salinity_evaporation.initArray_2D(jm, km, 0.0); // additional salinity by evaporation
    Evaporation_Dalton.initArray_2D(jm, km, 0.0); // evaporation by Dalton in [mm/d]
    Evaporation_Penman.initArray_2D(jm, km, 0.0); // evaporation by Penman in [mm/d]
    Precipitation.initArray_2D(jm, km, 0.0); // areas of higher precipitation
    precipitation_NASA.initArray_2D(jm, km, 0.0); // surface precipitation from NASA
    temperature_NASA.initArray_2D(jm, km, 0.0); // surface temperature from NASA
    temp_reconst.initArray_2D(jm, km, 0.0); // surface temperature from reconstruction tool
    c_fix.initArray_2D(jm, km, 0.0); // local surface salinity fixed for iterations
    v_wind.initArray_2D(jm, km, 0.0); // v-component of surface wind
    w_wind.initArray_2D(jm, km, 0.0); // w-component of surface wind
    temp_landscape.initArray_2D(jm, km, 0.0); // landscape temperature

    h.initArray(im, jm, km, 0.0); // bathymetry, depth from sea level
    t.initArray(im, jm, km, 1.0); // temperature
    u.initArray(im, jm, km, 0.0); // u-component velocity component in r-direction
    v.initArray(im, jm, km, 0.0); // v-component velocity component in theta-direction
    w.initArray(im, jm, km, 0.0); // w-component velocity component in phi-direction
    c.initArray(im, jm, km, 1.0); // salinity

    tn.initArray(im, jm, km, 1.0); // temperature new
    un.initArray(im, jm, km, 0.0); // u-velocity component in r-direction new
    vn.initArray(im, jm, km, 0.0); // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, 0.0); // w-velocity component in phi-direction new
    cn.initArray(im, jm, km, 1.0); // salinity new

    p_stat.initArray(im, jm, km, 0.0); // dynamic pressure
    p_hydro.initArray(im, jm, km, 1.0); // static pressure

    rhs_t.initArray(im, jm, km, 0.0); // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.0); // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.0); // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.0); // auxilliar field RHS w-velocity component
    rhs_c.initArray(im, jm, km, 0.0); // auxilliar field RHS water vapour

    aux_u.initArray(im, jm, km, 0.0); // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.0); // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.0); // auxilliar field w-velocity component

    Salt_Finger.initArray(im, jm, km, 0.0); // salt bulge of higher density
    Salt_Diffusion.initArray(im, jm, km, 0.0); // salt bulge of lowerer density and temperature
    Salt_Balance.initArray(im, jm, km, 0.0); // +/- salt balance

    r_water.initArray(im, jm, km, r_0_water); // water density as function of pressure
    r_salt_water.initArray(im, jm, km, r_0_water); // salt water density as function of pressure and temperature
    BuoyancyForce.initArray(im, jm, km, 0.0); // 3D buoyancy force
    CoriolisForce.initArray(im, jm, km, 0.0); // Coriolis force terms
    PressureGradientForce.initArray(im, jm, km, 0.0); // Force caused by normal pressure gradient
    cout << "      OGCM: reset_arrays ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::Run(){
    std::time_t Run_start;
    struct tm * timeinfo_begin;
    std::time(&Run_start);
    timeinfo_begin = std::localtime(&Run_start);
    std::cout << std::endl << std::endl;
    std::cout << " ... OGCM: time and date at run time begin:   " 
        << std::asctime(timeinfo_begin);
    std::cout << std::endl << std::endl;

    mkdir(output_path.c_str(), 0777);
    if(!has_printed_welcome_msg)  print_welcome_msg();
    for(int i = time_start; i <= time_end; i+=time_step){
        RunTimeSlice(i);
    }
    print_final_remarks();

    std::time_t Run_end;
    struct tm * timeinfo_end;
    std::time(&Run_end);
    timeinfo_end = std::localtime(&Run_end);
    std::cout << " ... OGCM: time and date at run time end:   " 
        << std::asctime(timeinfo_end) << std::endl;
    int Run_total = Run_end - Run_start;
    int Run_total_minutes = Run_total/60;
    int Run_total_hours = Run_total_minutes/60;
    std::cout <<  " ... OGCM: computer time needed for Ma = " << time_start 
        << " to Ma = " << time_end << " in " << time_step << " Ma steps:     "
        << Run_total << " seconds" << std::endl
        << " ... compares to:" << std::endl << std::endl
        << setw(30) << setfill(' ') << Run_total_hours 
        << " hours" << std::endl
        << setw(30) << setfill(' ') << Run_total_minutes 
        << " minutes" << std::endl
        << setw(30) << setfill(' ') << Run_total%60 
        << " seconds" << std::endl << std::endl
        << std::endl;
    return;
}
/*
*
*/
void cHydrosphereModel::write_file(std::string &bathymetry_name, 
    std::string &output_path, bool is_final_result){
    cout << endl << "      OGCM: write_file" << endl;
    int i_radial = 40;
    paraview_vtk_radial(bathymetry_name, i_radial, iter_cnt-1);
    int j_longal = 75;
    paraview_vtk_longal(bathymetry_name, j_longal, iter_cnt-1);
    int k_zonal = 185;
    paraview_vtk_zonal(bathymetry_name, k_zonal, iter_cnt-1);
    if(paraview_panorama_vts_flag){
        paraview_panorama_vts(bathymetry_name, iter_cnt-1);
    }
    HydrospherePlotData(bathymetry_name,(is_final_result ? -1 : iter_cnt-1));
    cout << endl << "      OGCM: write_file ended" << endl;
    return;
}
/*
*
*/
void  cHydrosphereModel::save_data(){
    cout << endl << "      OGCM: save_data" << endl;
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
    cout << "      OGCM: save_data ended" << endl;
    return;
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
                tn.x[i][j][k] = coeff * t.x[i][j][k];
                cn.x[i][j][k] = coeff * c.x[i][j][k];
            }
        }
    }
    return;
}
/*
*
*/
void cHydrosphereModel::store_intermediate_data_2D(float coeff){   
    for(int j = 0; j < jm; j++){   
        for(int k = 0; k < km; k++){   
            vn.x[im-1][j][k] = coeff * v.x[im-1][j][k];
            wn.x[im-1][j][k] = coeff * w.x[im-1][j][k]; 
        }
    }
    return;
}
/*
*
*/
void cHydrosphereModel::run_2D_loop(){
cout << endl << "      OGCM: run_2D_loop hyd ................................." << endl;
    int switch_2D = 0;
    iter_cnt = 1;
    if(switch_2D != 1){
        for(pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++){
            for(velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++){
                print_loop_2D_headings();
                BC_theta();
                BC_phi();
                BC_SolidGround();
                if(velocity_iter_2D % 2 == 0){
                     computePressure_2D();
                }
                solveRungeKutta_2D_Hydrosphere();
                ValueLimitationHyd();
                store_intermediate_data_2D(1.0);
                print_min_max_hyd();
                run_data_hyd(); 
                iter_cnt++;
            }  //  ::::::   end of velocity loop_2D: if (velocity_iter_2D > velocity_iter_max_2D)   ::::::::::::::::::::::

            if(pressure_iter_2D % checkpoint == 0){
            cout << endl << "      OGCM: write_file in run_2D_loop hyd ......................." << endl;
                write_file(bathymetry_name, output_path, false);
            }

            if(iter_cnt > nm){
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        } // :::::::::::::::::::   end of pressure loop_2D: if (pressure_iter_2D > pressure_iter_max_2D)   ::::::::::
    } // ::::::::   end of 2D loop for initial surface conditions: if (switch_2D == 0)   :::::::::::::::::::::::::::::
    cout << endl << "      OGCM: run_2D_loop hyd ended ................................" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::run_3D_loop(){
cout << endl << "      OGCM: run_3D_loop hyd ..........................." << endl;
    iter_cnt = 1;
    iter_cnt_3d = 0;
    double residuum_loop = 0.0;
    if(iter_cnt_3d == 0){
        residuum_old = 1.0;
        residuum_loop = 1.0;
    }
    if(debug){
        save_data();
    }
    for(pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++){
        for(velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++){
            print_loop_3D_headings();
            BC_radius();
            BC_theta();
            BC_phi();
            BC_SolidGround();
            if(velocity_iter % 2 == 0){
                computePressure_3D();
                residuum_loop = find_residuum_hyd();
                if(fabs(residuum_loop/residuum_old - 1.0) <= eps_residuum){
                    cout << endl << "      OGCM: write_file in run_3D_loop hyd, termination due to convergence .... relative error less than eps_residuum ..................." << endl;
                    cout << endl << "      residuum_loop = " << residuum_loop
                                 << "   residuum_old = " << residuum_old
                                 << "   eps_residuum = " << eps_residuum
                                 << "   absolute error = " 
                                 << fabs(residuum_old - residuum_loop)
                                 << "   relative error = " 
                                 << fabs(residuum_loop/residuum_old - 1.0) << endl;
                    write_file(bathymetry_name, output_path, false);
                    goto endloop;
                }else  cout << "      OGCM: write_file in find_residuum_hyd, relative error is too high ......................." << endl << endl;
                if(iter_cnt > nm){
                    cout << "       nm = " << nm 
                        << "     .....     maximum number of iterations   nm   reached!" 
                        << endl;
                    goto endloop;
                }
            }
/*
            if(velocity_iter % 4 == 0){
                PresStat_SaltWaterDens();
                SalinityEvaporation();
            }
*/
            solveRungeKutta_3D_Hydrosphere(); 
            ValueLimitationHyd();
            residuum_old = residuum_loop;
            store_intermediate_data_3D(1.0);
            print_min_max_hyd();
            residuum_old = residuum_loop;
            run_data_hyd();
            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        } // end of velocity loop
        if(pressure_iter % checkpoint == 0){
            cout << endl << "      OGCM: write_file in run_3D_loop hyd ......................." << endl;
            write_file(bathymetry_name, output_path, false);
        }
    } // end of pressure loop
    endloop:
    cout << endl << "      OGCM: run_3D_loop hyd ended ..........................." << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::init_EkmanSpiral(){
cout << endl << "      OGCM: init_EkmanSpiral" << endl;
    //Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
    //a further turning downwards until the end of the shear layer such that finally 90° of turning are reached
//    double wind_water_vel_ratio = 0.03;
    // initial conditions for v and w velocity components at the sea surface
    // ocean surface velocity is about 3% of the wind velocity at the surface
    // u_0 for the hydrosphere is (0.03 * u_0) for the atmosphere
    // surface wind vector driving the Ekman spiral in the Ekman layer
    // northern and southern hemisphere
    std::vector<double>radius(im, 0);
    cout.precision(8);
    double a = 0.0; // in 1/m
    double Ekman_angle = 45.0/pi180;
    double sinthe = 0.0;
    double rm = 0.0;
    double rmsinthe = 0.0;
    double a_z = 0.0;
    double exp_a_z = 0.0;
    double sin_a_z = 0.0;
    double cos_a_z = 0.0;
    double alfa = 0.0;
    double angle = 0.0;
    double CD = 2.6e-3;  // drag coefficient in ./.
    double U_10 = 0.0;  // wind velocity 10 m above sea surface directed to the north in m/s
    double V_0 = 0.0;  // water velocity at surface shifted by Ekman angle in m/s
    double T_yz = 0.0;  // wind stress in v-direction ( y, j ) in kg/(m*s2)
    double f = 0.0;  // Coriolis parameter in 1/s
    double Az = 0.0;  // constant vertical eddy viscosity in m2/s
//    double D_E = 0.0;  // Ekman layer depth in m
//    double D_E_op = 0.0;  // Ekman layer depth opposite to the surface wind in m
    radius[0] = 1.0;
    for(int i = 1; i < im; i++){
        radius[i] = radius[i-1] - dr;
    }
    radius[im-1] = 0.0;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){           // non-dimensional
            v_wind.y[j][k] = v.x[im-1][j][k];
            w_wind.y[j][k] = w.x[im-1][j][k];

            if(is_land(h, im-1, j, k)){
                v_wind.y[j][k] = 0.0; 
                w_wind.y[j][k] = 0.0;
            }
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
//    for(int k = 0; k < km; k++){
//        for(int j = 0; j < jm; j++){
            v_wind.y[90][k] = 0.0;
//            if(w_wind.y[j][k] == 0.0) w_wind.y[j][k] = 1.0e-8;
            alfa = atan(fabs(v_wind.y[j][k]/w_wind.y[j][k]));
//            if(w_wind.y[j][k] == 0.0) alfa = 0.0;
            sinthe = sin(the.z[j]);
            if(sinthe == 0.0) sinthe = 1.0e-5;
            U_10 = sqrt(v_wind.y[j][k] * v_wind.y[j][k] 
                + w_wind.y[j][k] * w_wind.y[j][k]); // in m/s, dimensional surface wind velocity U_10 in m/s
            // original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 139, eq. 9.16
            T_yz = r_air * CD * U_10 * U_10;  // in kg/(m*s*s))
            V_0 = 0.0127 * U_10/sqrt(sinthe);  // in m/s
            f = 2.0 * omega * fabs(sinthe);  // in 1/s
            Az = pow((T_yz/(r_0_water * V_0)), 2)/f;  // in m*m/s
            a = sqrt(f/(2.0 * Az));  // in 1/m
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
//        << "   D_E = " << D_E << "   D_E_op = " << D_E_op << endl
        << "   T_yz = " << T_yz << endl
        << "   f = " << f
        << "   Az = " << Az
        << "   a = " << a
//        << "   i_Ekman = " << i_Ekman 
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
//        << "   D_E = " << D_E << "   D_E_op = " << D_E_op << endl
        << "   T_yz = " << T_yz << endl
        << "   f = " << f
        << "   Az = " << Az
        << "   a = " << a
//        << "   i_Ekman = " << i_Ekman 
        << endl << endl;
*/
            if(j <= (jm-1)/2){
                if((w_wind.y[j][k] <= 0.0)&&(v_wind.y[j][k] >= 0.0)){
                    if((alfa >= 0.0)&&(alfa <= 45.0/pi180)){  // section I (j = 83 - 90) (7°N - 0°N)
                        angle = 180.0/pi180 - (alfa - Ekman_angle);
                    }
                    if((alfa > 45.0/pi180)&&(alfa <= 90.0/pi180)){  // section II (j = 76 - 83) (14°N - 7°N)
                        angle = 180.0/pi180 - (alfa - Ekman_angle);
                    }
                }
                if((w_wind.y[j][k] >= 0.0)&&(v_wind.y[j][k] >= 0.0)){
                    if((alfa <= 90.0/pi180)&&(alfa >= 45.0/pi180)){  // section III (j = 67 - 76) (23°N - 14°N)
                        angle = alfa + Ekman_angle;
                    }
                    if((alfa < 45.0/pi180)&&(alfa >= 0.0/pi180)){  // section IV (j = 49 - 67) (41°N - 23°N)
                        angle = alfa + Ekman_angle;
                    }
                }
                if((w_wind.y[j][k] >= 0.0)&&(v_wind.y[j][k] <= 0.0)){  // section V (j = 33 - 49) (57°N - 41°N)
                    if((alfa >= 0.0)&&(alfa <= 45.0/pi180)){
                        angle = alfa + Ekman_angle;
                    }
                }
//              original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 137, eq. 9.9a/b, p. 138, eq. 9.11a/b
                for(int i = 0; i < im; i++){
                    a_z = - a * radius[i] * L_hyd * dr;
                    exp_a_z = exp(a_z);
                    sin_a_z = sin(angle + a_z);
                    cos_a_z = cos(angle + a_z);
                    v.x[i][j][k] = V_0 * exp_a_z * sin_a_z; 
                    w.x[i][j][k] = V_0 * exp_a_z * cos_a_z;
                    v.x[i][90][k] = 0.0;
                    if(is_land(h, i, j, k)){
                        v.x[i][j][k] = 0.0; 
                        w.x[i][j][k] = 0.0;
                    }
/*
//    if((j == 75) &&(k == 180)) cout << "north" << endl
    if((j == 1) &&(k == 180)) cout << "north" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   rad = " << radius[i] << endl
        << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
        << "   a = " << a << "   a_z = " << a_z << "   exp_a_z = " << exp_a_z 
        << "   sin_a_z = " << sin_a_z << "   cos_a_z = " << cos_a_z << endl
        << "   alfa = " << alfa
        << "   angle = " << angle << endl
        << "   alfa_deg = " << alfa * pi180
        << "   angle_deg = " << angle * pi180 << endl
//        << "   Ekman = " << i_Ekman 
//        << "   D_E = " << D_E << "   D_E_op = " << D_E_op << endl
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
                } // i end north
            }
            if(j > (jm-1)/2){
                if((w_wind.y[j][k] <= 0.0)&&(v_wind.y[j][k] <= 0.0)){
                    if((alfa >= 0.0)&&(alfa <= 45.0/pi180)){  // section I (j = 173 - 180) (7°S - 0°S)
                        angle = 180.0/pi180 - (alfa - Ekman_angle);
                    }
                    if((alfa > 45.0/pi180)&&(alfa <= 90.0/pi180)){  // section II (j = 166 - 173) (14°S - 7°S)
                        angle = 180.0/pi180 - (alfa - Ekman_angle);
                    }
                }
                if((w_wind.y[j][k] >= 0.0)&&(v_wind.y[j][k] <= 0.0)){
                    if((alfa > 45.0/pi180)&&(alfa <= 90.0/pi180)){  // section III (j = 157 - 166) (23°S - 14°S)
                        angle = alfa + Ekman_angle;
                    }
                    if((alfa >= 0.0)&&(alfa <= 45.0/pi180)){  // section IV (j = 139 - 157) (41°S - 23°S)
                        angle = alfa + Ekman_angle;
                    }
                }
                if((w_wind.y[j][k] >= 0.0)&&(v_wind.y[j][k] >= 0.0)){
                    if((alfa < 45.0/pi180)&&(alfa >= 0.0/pi180)){  // section V (j = 123 - 139) (57°S - 41°S)
                        angle = alfa + Ekman_angle;
                    }
                }
//              original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 137, eq. 9.9a/b, p. 138, eq. 9.11a/b
                for(int i = 0; i < im; i++){
                    a_z = - a * radius[i] * L_hyd * dr;
                    exp_a_z = exp(a_z);
                    sin_a_z = - sin(angle + a_z);
                    cos_a_z = cos(angle + a_z);
                    v.x[i][j][k] = V_0 * exp_a_z * sin_a_z; 
                    w.x[i][j][k] = V_0 * exp_a_z * cos_a_z;
                    v.x[i][90][k] = 0.0;
                    if(is_land(h, i, j, k)){
                        v.x[i][j][k] = 0.0; 
                        w.x[i][j][k] = 0.0;
                    }

/*
    if((j == 105) &&(k == 180)) cout << "south" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   rad = " << radius[i] << endl
        << "   sinthe = " << sinthe << "   sqrt(sinthe) = " << sqrt(sinthe) << endl
        << "   a = " << a << "   a_z = " << a_z << "   exp_a_z = " << exp_a_z 
        << "   sin_a_z = " << sin_a_z << "   cos_a_z = " << cos_a_z << endl
        << "   alfa = " << alfa
        << "   angle = " << angle << endl
        << "   alfa_deg = " << alfa * pi180
        << "   angle_deg = " << angle * pi180 << endl
//        << "   Ekman = " << i_Ekman 
//        << "   D_E = " << D_E << "   D_E_op = " << D_E_op << endl
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
                } // i end south
            }
        }
    }

    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
//            v.x[i][0][k] = 0.0;
//            w.x[i][0][k] = 0.0;
            v.x[i][0][k] = c43 * v.x[i][1][k] - c13 * v.x[i][2][k];
            w.x[i][0][k] = c43 * w.x[i][1][k] - c13 * w.x[i][2][k];
//            v.x[i][jm-1][k] = 0.0;
//            w.x[i][jm-1][k] = 0.0;
            v.x[i][jm-1][k] = c43 * v.x[i][jm-2][k] - c13 * v.x[i][jm-3][k];
            w.x[i][jm-1][k] = c43 * w.x[i][jm-2][k] - c13 * w.x[i][jm-3][k];
        }
    }

    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            v.x[i][j][0] = c43 * v.x[i][j][1] - c13 * v.x[i][j][2];
            w.x[i][j][0] = c43 * w.x[i][j][1] - c13 * w.x[i][j][2];
            v.x[i][j][km-1] = c43 * v.x[i][j][km-2] - c13 * v.x[i][j][km-3];
            w.x[i][j][km-1] = c43 * w.x[i][j][km-2] - c13 * w.x[i][j][km-3];
            v.x[i][j][0] = v.x[i][j][km-1] 
                = (v.x[i][j][0] + v.x[i][j][km-1])/2.0;
            w.x[i][j][0] = w.x[i][j][km-1] 
                = (w.x[i][j][0] + w.x[i][j][km-1])/2.0;
        }
    }

    for(int j = 1; j < jm-1; j++){
        double dr_rm = dr * L_hyd/(double)(im-1);
        sinthe = sin(the.z[j]);
        if(sinthe == 0.0) sinthe = 1.0e-5;
        for(int k = 1; k <= km-1; k++){
                u.x[0][j][k] = 0.0;
            for(int i = 1; i <= im-1; i++){
                rm = - L_hyd * (1.0 - rad.z[i]);
                rmsinthe = rm * sinthe;
                double h_d_j = 1.0;
                double h_d_k = 1.0;
                std::vector<Array*> arrays_1{&v, &w};
                enum array_index_1{i_v_1, i_w_1, last_array_index_1};
                std::vector<double> dxdthe_vals(last_array_index_1), 
                                    dxdphi_vals(last_array_index_1);
                bool the_flag = false, phi_flag = false;
                if((j >= 1)&&(j < jm-2)){
                    if((is_land(h, i, j, k))
                        &&((is_water(h, i, j+1, k)) 
                        &&(is_water(h, i, j+2, k)))){
                        for(int n=0; n<last_array_index_1; n++)
                            dxdthe_vals[n] = dxdthe_b(arrays_1[n]);
                        the_flag = true;
                    }
                    if((is_land(h, i, j, k))
                        &&(is_water(h, i, j-1, k)) 
                        &&(is_water(h, i, j-2, k))){
                        for(int n=0; n<last_array_index_1; n++)
                            dxdthe_vals[n] = dxdthe_c(arrays_1[n]);
                        the_flag = true;
                    }
                }
                if((k >= 1)&&(k < km-2)){
                    if((is_land(h, i, j, k))
                        &&(is_water(h, i, j, k+1)) 
                        &&(is_water(h, i, j, k+2))){
                        for(int n=0; n<last_array_index_1; n++)
                            dxdphi_vals[n] = dxdphi_b(arrays_1[n]);
                        phi_flag = true;
                    }
                    if((is_land(h, i, j, k))
                        &&(is_water(h, i, j, k-1)) 
                        &&(is_water(h, i, j, k-2))){
                        for(int n=0; n<last_array_index_1; n++)
                            dxdphi_vals[n] = dxdphi_c(arrays_1[n]);
                        phi_flag = true;
                    }
                }
                for(int n=0; n<last_array_index_1; n++){
                    if(!the_flag) dxdthe_vals[n] = dxdthe_a(arrays_1[n]);
                    if(!phi_flag) dxdphi_vals[n] = dxdphi_a(arrays_1[n]);
                }
                double dvdthe = dxdthe_vals[i_v_1];
                double dwdphi = dxdphi_vals[i_w_1];

                u.x[i][j][k] = u.x[i-1][j][k] 
                    - dr_rm * (dvdthe/rm + dwdphi/rmsinthe);
/*
                double residuum = (u.x[i][j][k] - u.x[i-1][j][k])/dr 
                    + ((v.x[i-1][j+1][k] - v.x[i-1][j-1][k])
                    /(2. * rm * dthe) 
                    + (w.x[i-1][j][k+1] - w.x[i-1][j][k-1])
                    /(2. * rmsinthe * dphi));

                if((j == 75) &&(k == 180)) cout << "north pumping" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                    << "   u = " << u.x[i][j][k] 
                    << "   v = " << v.x[i][j][k] 
                    << "   w = " << w.x[i][j][k] << endl
                    << "   rm = " << rm
                    << "   dr = " << dr
                    << "   residuum = " << residuum << endl;
                if((j == 105) &&(k == 180)) cout << "south pumping" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                    << "   u = " << u.x[i][j][k] 
                    << "   v = " << v.x[i][j][k] 
                    << "   w = " << w.x[i][j][k] << endl
                    << "   rm = " << rm
                    << "   dr = " << dr
                    << "   residuum = " << residuum << endl;
*/
            }  // end i
        }  // end k
    }  // end j

/*
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
//            u.x[0][j][k] = c43 * u.x[1][j][k] - c13 * u.x[2][j][k];
            u.x[im-1][j][k] = u.x[im-4][j][k] 
                - 3. * u.x[im-3][j][k] + 3. * u.x[im-2][j][k];  // extrapolation
            u.x[0][j][k] = 0.0;
//            u.x[im-1][j][k] = c43 * u.x[im-2][j][k] - c13 * u.x[im-3][j][k];
//            u.x[im-1][j][k] = 0.0;
        }
    }
*/

    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            u.x[i][0][k] = c43 * u.x[i][1][k] - c13 * u.x[i][2][k];
//            u.x[i][0][k] = 0.0;
            u.x[i][jm-1][k] = c43 * u.x[i][jm-2][k] - c13 * u.x[i][jm-3][k];
//            u.x[i][jm-1][k] = 0.0;
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            u.x[i][j][0] = c43 * u.x[i][j][1] - c13 * u.x[i][j][2];
            u.x[i][j][0] = u.x[i][j][km-1] 
                = (u.x[i][j][0] + u.x[i][j][km-1])/2.0;
        }
    }

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                u.x[i][j][k] = u.x[i][j][k]/u_0;
//                u.x[i][j][k] = 0.0;
                v.x[i][j][k] = v.x[i][j][k]/u_0;
                w.x[i][j][k] = w.x[i][j][k]/u_0;

                if((j>=89)&&(j<=90)){
                    u.x[i][j][k] = u.x[i][88][k];
                }
                if((j>90)&&(j<=91)){
                    u.x[i][j][k] = u.x[i][88][k];
                }
                if(is_land(h, i, j, k)){
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;
                }
            }
        }
    }
    cout << "      OGCM: init_EkmanSpiral ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::init_bathymetry(const string &bathymetry_file){
    cout << endl << "      OGCM: init_bathymetry" << endl;
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
                    &&(is_water(h, i, j+1, k)))
                    &&((is_water(h, i, j, k-1))
                    &&(is_water(h, i, j, k+1)))){
                    h.x[i][j][k] = 0.0;
                }
            }
        }
    }
    cout << "      OGCM: init_bathymetry ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::init_temperature(int Ma){
cout << endl << "      OGCM: init_temperature" << endl;
    double t_paleo_add = 0.0; 
    double t_pole_add = 0.0; 
    double t_global_mean_exp = 0.0;

    // correction of surface temperature around 180°E due to bad data around 180E
    if(is_first_time_slice()){
        int k_half = (km -1)/2;
        for(int j = 0; j < jm; j++){
            t.x[im-1][j][k_half] = (t.x[im-1][j][k_half+1] 
                + t.x[im-1][j][k_half-1])/2.0;
            temperature_NASA.y[j][k_half] = (temperature_NASA.y[j][k_half+1] +
                temperature_NASA.y[j][k_half-1])/2.0;
        }
    }
    // if use_earthbyte_reconstruction temperature in °C converted to non-dimensional
    if(Ma != 0 && use_earthbyte_reconstruction){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                t.x[im-1][j][k] = (t.x[im-1][j][k] + t_0)/t_0; // non-dimensional, (reconstructed temperature in °C)
                temp_reconst.y[j][k] = t.x[im-1][j][k]; // non-dimensional
            }
        }
    }
    if((is_first_time_slice())&&(use_earthbyte_reconstruction)&&(!use_NASA_temperature)){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                t.x[im-1][j][k] = (temperature_NASA.y[j][k] + t_0)/t_0; // non-dimensional
            }
        }
    }
    t_average = 
        get_temperatures_from_curve(0, m_equat_temperature_curve);
    t_pole_modern = 
        get_temperatures_from_curve(0, m_pole_temperature_curve);
    t_global_mean_exp = 
        get_temperatures_from_curve(*get_current_time(), 
        m_global_temperature_curve);
    if(is_first_time_slice()){
        t_global_mean_exp = 
            get_temperatures_from_curve(0, 
                m_global_temperature_curve);
        t_global_mean = GetMean_2D(jm, km, temperature_NASA);
    }
    if(!is_first_time_slice()){                                         // temperature difference between adjacent time steps
        t_paleo_add =                                                   // at poles
            get_temperatures_from_curve(*get_current_time(), 
            m_equat_temperature_curve)
            - get_temperatures_from_curve(*get_previous_time(), 
            m_equat_temperature_curve);
        t_pole_add =                                                    // at equator
            get_temperatures_from_curve(*get_current_time(), 
            m_pole_temperature_curve)
            - get_temperatures_from_curve(*get_previous_time(), 
            m_pole_temperature_curve);
        t_paleo_add /= t_0; // non-dimensional 
        t_pole_add /= t_0; // non-dimensional 
    }
    t_paleo_total += t_paleo_add;
    t_pole_total += t_pole_add;

    cout.precision(3);
    const char* time_slice_comment = "      time slice of Paleo-OGCM: ";
    const char* time_slice_number = " Ma = ";
    const char* time_slice_unit = " million years";

    cout << setiosflags(ios::left) << setw(58) << setfill('.') 
        << time_slice_comment << resetiosflags(ios::left) << setw (6) 
        << fixed << setfill(' ') << time_slice_number << setw (3) << Ma 
        << setw(12) << time_slice_unit << endl;

    const char* temperature_equat_comment = "      equat temperature increase: ";  // t global increase
    const char* temperature_equat_gain = " t equat increase";

    const char* temperature_pole_comment = "      pole temperature increase: ";  // t pole increase
    const char* temperature_pole_gain = " t pole increase";

    const char* temperature_equat = "      equat temperature at paleo times: ";  // t equatorial
    const char* temperature_average_pal = " t paleo";

    const char* temperature_pole = "      pole temperature at paleo times: ";  // t polar
    const char* temperature_pole_pal = " t pole";

    const char* temperature_mean = "      mean temperature at paleo times: ";  // t global mean ATOM
    const char* temperature_mean_pal = " t global mean";

    const char* temperature_mean_exp = "      exp mean temperature at paleo times: ";  // t global mean experimental
    const char* temperature_mean_exp_pal = " t global mean exp";

    const char* temperature_equat_modern = "      equat temperature at modern times: ";  // t modern equatorial
    const char* temperature_average_equat = " t modern equat";

    const char* temperature_pole_modern = "      pole temperature at modern times: ";  // t modern pole
    const char* temperature_average_pole = " t modern pole";

    const char* temperature_unit =  "°C ";

    cout << setiosflags(ios::left) 
        << setw(51) << setfill('.')  // t global increase
        << temperature_equat_comment << resetiosflags(ios::left) << setw(13) 
        << temperature_equat_gain << " = " << setw(7) << setfill(' ') 
        << t_paleo_add * t_0 << setw(5) << temperature_unit << endl 

        << setiosflags(ios::left)  // t pole increase
        << setw(52) << setfill('.') 
        << temperature_pole_comment << resetiosflags(ios::left) << setw(13) 
        << temperature_pole_gain << " = " << setw(7) << setfill(' ') 
        << t_pole_add * t_0 << setw(5) << temperature_unit << endl 

        << setw(55) << setfill('.')  // t equatorial
        << setiosflags(ios::left) << temperature_equat << resetiosflags(ios::left) 
        << setw(13) << temperature_average_pal << " = "  << setw(7) 
        << setfill(' ') << ((t_average + t_0) + t_paleo_total * t_0) - t_0 << setw(5) 
        << temperature_unit << endl

        << setw(55) << setfill('.')  // t polar
        << setiosflags(ios::left) << temperature_pole << resetiosflags(ios::left) 
        << setw(13) << temperature_pole_pal << " = "  << setw(7) 
        << setfill(' ') << ((t_pole_modern + t_0) + t_pole_total * t_0) - t_0 << setw(5)
        << temperature_unit << endl

        << setw(54) << setfill('.')  // t global mean ATOM
        << setiosflags(ios::left) << temperature_mean << resetiosflags(ios::left) 
        << setw(13) << temperature_mean_pal << " = "  << setw(7) 
        << setfill(' ') << t_global_mean << setw(5)
        << temperature_unit << endl

        << setw(50) << setfill('.')  // t global mean experimental
        << setiosflags(ios::left) << temperature_mean_exp << resetiosflags(ios::left) 
        << setw(13) << temperature_mean_exp_pal << " = "  << setw(7) 
        << setfill(' ') << t_global_mean_exp << setw(5)
        << temperature_unit << endl

        << setw(53) << setfill('.')  // t modern equatorial
        << setiosflags(ios::left) << temperature_equat_modern 
        << resetiosflags(ios::left) << setw(13) << temperature_average_equat 
        << " = " << setw(7) << setfill(' ') << t_average << setw(5) 
        << temperature_unit << endl

        << setw(54) << setfill('.')  // t modern pole
        << setiosflags(ios::left) << temperature_pole_modern 
        << resetiosflags(ios::left) << setw(13) << temperature_average_pole 
        << " = " << setw(7) << setfill(' ') << t_pole_modern << setw(5) 
        << temperature_unit << endl;

    double t_corr = 0.0;
    double d_j_half = (double)(jm-1)/2.0;
    double t_equator = (get_temperatures_from_curve(*get_current_time(), 
                m_equat_temperature_curve) + t_0)/t_0;
    double t_pole = (get_temperatures_from_curve(*get_current_time(), 
                m_pole_temperature_curve) + t_0)/t_0;
    double t_eff = t_pole - t_equator;
    double t_eff_rec = t_pole_add - t_paleo_add;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double d_j = (double)j;
            t_corr = (t_eff_rec                                         // average temperature distributions in time at poles and equator corrects parabolic assumption
                * parabola((double)d_j/(double)d_j_half)                // parabolic correction function from pole to pole, values t_corr to be added  
                + t_pole_add);                                          // on any pole to pole reconstructed initial temperature assumtion
            if((!use_NASA_temperature)&&(!use_earthbyte_reconstruction)){  // parabolic surface land and ocean temperature assumed
                t.x[im-1][j][k] = t_eff 
                    * parabola((double)d_j/(double)d_j_half) + t_pole;  // parabolic temperature assumption
                t.x[im-1][j][k] = t.x[im-1][j][k] + t_corr;             // based on the parabolic temperature assumption
            }
            if((use_NASA_temperature)&&(use_earthbyte_reconstruction)){ // NASA surface temperature assumed and corrected by known local values
                if(*get_current_time() == 0){
                    t.x[im-1][j][k] = (temperature_NASA.y[j][k] + t_0)/t_0;  // initial temperature by NASA for Ma=0, non-dimensional
                }else{
                    t.x[im-1][j][k] = t.x[im-1][j][k] + t_corr;         // based on the reconstructed temperature
                }
                if(*get_current_time() >= Ma_switch){                   // parabolic temperature distribution starting at Ma_switch
                    t.x[im-1][j][k] = t_eff 
                        * parabola((double)d_j/(double)d_j_half) + t_pole;
                    t.x[im-1][j][k] = t.x[im-1][j][k] + t_corr;         // based on the parabolic temperature assumption
                }
            }
            temp_landscape.y[j][k] = t.x[im-1][j][k] * t_0 - t_0; // in °C
        }// for j
    }// for k
    t_global_mean = GetMean_2D(jm, km, temp_landscape);
    int i_beg = 0; // 200m depth  
    double tm_tbeg = 0.0;
    int i_max = im - 1;   // corresponds to sea level
    double d_i_max =(double)i_max;
    double d_i_beg =(double)i_beg;
    double d_i = 0.0;
// distribution of t with increasing depth
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            tm_tbeg = (t.x[im-1][j][k] - 1.0) 
                /(d_i_max * d_i_max - d_i_beg * d_i_beg);
            for(int i = i_beg; i < im-1; i++){
                    d_i = (double)i;
                    t.x[i][j][k] = 1.0 + tm_tbeg 
                        * (d_i * d_i - d_i_beg * d_i_beg);// parabolic approach
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(is_land(h, i, j, k)){
                    t.x[i][j][k] = 1.0;
                }
            }
        }
    }
    cout << "      OGCM: init_temperature ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::init_salinity(){
cout << endl << "      OGCM: init_salinity" << endl;
    // Lenton_etal_COPSE_time_temp, constant paleo mean temperature, added to the surface initial temperature
    // difference between mean temperature (Ma) and mean temperature (previous Ma) == t_paleo_add
    double t_paleo_add = 0; 
    if(Ma > 0){
        if((use_NASA_temperature)&&(*get_current_time() > 0))  
            t_paleo_add = 
                get_temperatures_from_curve(*get_current_time(), 
                m_global_temperature_curve)
                - get_temperatures_from_curve(*get_previous_time(), 
                m_global_temperature_curve);
        if(!use_NASA_temperature) 
            t_paleo_add = get_temperatures_from_curve(*get_current_time(), 
                m_global_temperature_curve) - t_average;
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
    int i_max = im - 1;   // corresponds to sea level
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
    cout << "      OGCM: init_salinity ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::init_dynamic_pressure(){
cout << endl << "      OGCM: init_dynamic_pressure" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
//                p_stat.x[i][j][k] = 0.5 * r_salt_water.x[i][j][k] 
//                p_stat.x[i][j][k] = 1e-2 * 0.5 * r_0_water // hPa
                p_stat.x[i][j][k] = 1e-2 * 0.5 * r_0_saltwater // hPa
                    * sqrt(((u.x[i][j][k] * u.x[i][j][k]) 
                    + (v.x[i][j][k] * v.x[i][j][k]) 
                    + (w.x[i][j][k] * w.x[i][j][k]))/3.0) * u_0;
                if(is_land(h, i, j, k))
                    p_stat.x[i][j][k] = 0.0;
            }
        }
    }
    cout << "      OGCM: init_dynamic_pressure ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::ValueLimitationHyd(){
cout << endl << "      Value_Limitation_Hyd" << endl;
// the limiting values depend on local singular behaviour
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){

                if(u.x[i][j][k] >= 0.06/u_0)  u.x[i][j][k] = 0.06/u_0; // non-dimensional
                if(u.x[i][j][k] <= - 0.06/u_0)  u.x[i][j][k] = - 0.06/u_0; // non-dimensional

                if(v.x[i][j][k] >= 0.06/u_0)  v.x[i][j][k] = 0.06/u_0;
                if(v.x[i][j][k] <= - 0.06/u_0)  v.x[i][j][k] = - 0.06/u_0;
                if(w.x[i][j][k] >= 0.06/u_0)  w.x[i][j][k] = 0.06/u_0;
                if(w.x[i][j][k] <= - 0.06/u_0)  w.x[i][j][k] = - 0.06/u_0;

                if(p_stat.x[i][j][k] >= 0.3)  p_stat.x[i][j][k] = 0.3; // hPa
                if(p_stat.x[i][j][k] < 0.0) p_stat.x[i][j][k] = 0.0;
                if(t.x[i][j][k] >= 1.147)  t.x[i][j][k] = 1.147; // 40.15 °C
                if(t.x[i][j][k] <= 0.9927)  t.x[i][j][k] = 0.9927;// -1.0 °C
                if(c.x[i][j][k] * c_0 >= 50.)  c.x[i][j][k] = 1.4451;  // 50.0 psu
                if(c.x[i][j][k] * c_0 <= 32.)  c.x[i][j][k] = 0.9249;      // 32.0 psu
            }
        }
    }
cout << "      Value_Limitation_Hyd ended" << endl;
}
/*
*
*/
void cHydrosphereModel::BC_SolidGround(){
    cout << endl << "      OGCM: BC_SolidGround" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                if(is_land(h, i, j, k)){
                    p_stat.x[i][j][k] = 0.0;
                    t.x[i][j][k] = 1.0;
                    c.x[i][j][k] = 1.0;
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;
                    un.x[i][j][k] = 0.0;
                    vn.x[i][j][k] = 0.0;
                    wn.x[i][j][k] = 0.0;
                    r_water.x[i][j][k] = r_0_water;
                    r_salt_water.x[i][j][k] = r_0_saltwater;
                }
            }
        }
    }
    cout << "      OGCM: BC_SolidGround" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::IC_t_WestEastCoast(){
    cout << endl << "      OGCM: IC_t_WestEastCoast" << endl;
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
            if((is_water(h, im-1, j, k-1))&&(is_land(h, im-1, j, k))){
                while(k >= 1){
                    if(is_water(h, im-1, j, k))  break;
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
            if((is_water(h, im-1, j, k-1))&&(is_land(h, im-1, j, k))){
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
    cout << "      OGCM: IC_t_WestEastCoast" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::IC_u_WestEastCoast(){
    cout << endl << "      OGCM: IC_u_WestEastCoast" << endl;
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
    cout << "      OGCM: IC_u_WestEastCoast" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::IC_Equatorial_Currents(){
    cout << endl << "      OGCM: IC_Equatorial_Currents" << endl;
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
    cout << "      OGCM: IC_Equatorial_Currents ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::IC_CircumPolar_Current(){ 
    cout << endl << "      OGCM: IC_CircumPolar_Current" << endl;
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
    cout << "      OGCM: IC_CircumPolar_Current ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::load_global_temperature_curve(){
    load_map_from_file(temperature_global_file, m_global_temperature_curve);
    return;
}
/*
*
*/
void cHydrosphereModel::load_equat_temperature_curve(){
    load_map_from_file(temperature_equat_file, m_equat_temperature_curve);
    return;
}
/*
*
*/
void cHydrosphereModel::load_pole_temperature_curve(){
    load_map_from_file(temperature_pole_file, m_pole_temperature_curve);
    return;
}
/*
*
*/
float cHydrosphereModel::get_temperatures_from_curve(float time, 
    std::map<float, float>& m) const{
    if(time < m.begin()->first 
        || time > (--m.end())->first){
        std::cout << "Input time out of range: " << time << std::endl;    
        return NAN;
    }
    if(m.size() < 2){
        std::cout << "No enough data in map m" << std::endl;
        return NAN;
    }
    map<float, float>::const_iterator upper = m.begin(), 
        bottom = ++m.begin(); 
    for(map<float, float>::const_iterator it = m.begin();
            it != m.end(); ++it){
        if(time < it->first){
            bottom = it;
            break;
        }else{
            upper = it;
        }
    }
/*
    std::cout << "   get_temperatures_from_curve" << std::endl;
    std::cout << "   Ma ->   " << upper->first << " " 
        << bottom->first << std::endl;
    std::cout << "   temp-range ->   "<< upper->second 
        << " " << bottom->second << std::endl;
    std::cout << "   temp-interpolation ->   " 
        << upper->second + (time - upper->first) 
       /(bottom->first - upper->first) 
        * (bottom->second - upper->second) << std::endl << std::endl;
*/
    return upper->second + (time - upper->first) 
       /(bottom->first - upper->first) 
        * (bottom->second - upper->second);
}
/*
*
*/
double cHydrosphereModel::find_residuum_hyd(){
    cout << endl << "      OGCM: find_residuum_hyd" << endl;
    bool residum_found = false;
    double residuum_hyd = 0.0;
    double maxValue = 0.0;
    double dudr = 0.0;
    double dvdthe = 0.0;
    double dwdphi = 0.0;
    int i_error = 0;
    int j_error = 0;
    int k_error = 0;  
//    for(int j = 1; j < jm-1; j++){
    for(int j = 75; j <= 75; j++){
        double sinthe = sin(the.z[j]);
        if(sinthe == 0.0) sinthe = 1.0e-5;
//        for(int k = 1; k < km-1; k++){
        for(int k = 180; k <= 180; k++){
            for(int i = 1; i < im-3; i++){
                double rm = rad.z[i];
                double rmsinthe = rm * sinthe;
                if((is_water(h, i, j, k)&&(is_water(h, i-1, j, k)))
                    &&(is_water(h, i, j, k)&&(is_water(h, i, j-1, k))
                    &&(is_water(h, i, j+1, k)))
                    &&(is_water(h, i, j, k)&&(is_water(h, i, j, k-1))
                    &&(is_water(h, i, j, k+1)))){
                    dudr = (u.x[i+1][j][k] - u.x[i-1][j][k])/dr;
                    dvdthe = (v.x[i][j+1][k] - v.x[i][j-1][k])
                        /(2.0 * rm * dthe);
                    dwdphi = (w.x[i][j][k+1] - w.x[i][j][k-1])
                        /(2.0 * rmsinthe * dphi);
                    residuum_hyd = fabs(dudr + dvdthe + dwdphi);
                    if(residuum_hyd > maxValue){
                        maxValue = residuum_hyd;
                        i_error = i;
                        j_error = j;
                        k_error = k;
                    }
/*
                    cout.precision(8);
                    if((j == 75) &&(k == 180)) cout << "residuum_hyd" << endl
                        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                        << "   u = " << u.x[i][j][k] 
                        << "   v = " << v.x[i][j][k] 
                        << "   w = " << w.x[i][j][k] << endl
                        << "   dudr = " << dudr
                        << "   dvdthe = " << dvdthe
                        << "   dwdphi = " << dwdphi << endl
                        << "   residuum_old = " << residuum_old
                        << "   residuum_hyd = " << residuum_hyd
                        << "   maxValue = " << maxValue
                        << "   eps_residuum = " << eps_residuum << endl
                        << "   i_error = " << i_error << "   j_error = " << j_error << "   k_error = " << k_error  << endl
                        << "   relative error = " 
                        << fabs(residuum_hyd/residuum_old - 1.0) << endl
                        << "   absolute error = " 
                        << fabs(residuum_old - maxValue) << endl << endl;
*/
                }  // end if loop
            }  //  end i
        }  // end k
    }  // end j
                    cout.precision(8);
                    if((residuum_old - maxValue) > 0.0){
                        cout << endl << "      OGCM: write_file in find_residuum_hyd, absolute error declining ......................." << endl;
                        cout << endl << "      residuum_hyd = " << maxValue
                                     << "      residuum_old = " << residuum_old
                                     << "      eps_residuum = " << eps_residuum << endl
                                     << "      i_error = " << i_error << "   j_error = " << j_error << "   k_error = " << k_error  << endl
                                     << "      relative error = " 
                                     << fabs(residuum_hyd/residuum_old - 1.0) << endl
                                     << "      absolute error = " 
                                     << fabs(residuum_old - maxValue) << endl;
                        residum_found = true;
                    }
    if(residum_found == false){  
        cout << "      OGCM: write_file in find_residuum_hyd, absolute error is too high ......................." << endl;
        cout << endl << "      residuum_hyd = " << maxValue
             << "      residuum_old = " << residuum_old
             << "      eps_residuum = " << eps_residuum << endl
             << "      i_error = " << i_error 
             << "      j_error = " << j_error 
             << "      k_error = " << k_error  << endl
             << "      relative error = " 
             << fabs(maxValue/residuum_old - 1.0) << endl
             << "      absolute error = " 
             << fabs(residuum_old - maxValue) << endl;
    }
    return maxValue;
}
/*
*
*/

