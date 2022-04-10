/*
 * Atmosphere General Circulation Modell(AGCM) applied to laminar flow
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
#include "cAtmosphereModel.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

extern std::vector<std::vector<double> > m_node_weights;

cAtmosphereModel* cAtmosphereModel::m_model = NULL;

const double cAtmosphereModel::pi180 = 180./ M_PI;      // pi180 = 57.3

const double cAtmosphereModel::the_degree = 1.;         // compares to 1° step size laterally
const double cAtmosphereModel::phi_degree = 1.;         // compares to 1° step size longitudinally

const double cAtmosphereModel::dr = 0.025;    // 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
//const double cAtmosphereModel::dt = 0.00001;  // time step coincides with the CFL condition
const double cAtmosphereModel::dt = 0.0001;  // time step coincides with the CFL condition
//const double cAtmosphereModel::dt = 0.001;  // time step crashes

const double cAtmosphereModel::dthe = the_degree/pi180; 
const double cAtmosphereModel::dphi = phi_degree/pi180;
    
const double cAtmosphereModel::the0 = 0.;             // North Pole
const double cAtmosphereModel::phi0 = 0.;             // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cAtmosphereModel::r0 = 1.0; // non-dimensional
//const double cAtmosphereModel::r0 = 0.0; // non-dimensional

const double cAtmosphereModel::residuum_ref_atm = 1.0e-2;     // criterium to finish global iterations

cAtmosphereModel::cAtmosphereModel():
    i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
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
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.initArray_1D(im, 0); // radial coordinate direction
    the.initArray_1D(jm, 0); // lateral coordinate direction
    phi.initArray_1D(km, 0); // longitudinal coordinate direction
}

cAtmosphereModel::~cAtmosphereModel(){
    if(PythonStream::is_enable()){
        std::cout.rdbuf(backup);
    }
    m_model = NULL;
    logger().close();
}
 
#include "cAtmosphereDefaults.cpp.inc"

void cAtmosphereModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile (filename);
    try{
        if(err){
            doc.PrintError();
            throw std::invalid_argument(std::string("unable to load config file:  ") 
                + filename);
        }
        XMLElement *atom = doc.FirstChildElement("atom"), *elem_common = NULL, 
            *elem_atmosphere = NULL;
        if (!atom){
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
            elem_atmosphere = atom->FirstChildElement("atmosphere");
            if (!elem_atmosphere) {
                throw std::invalid_argument(std::string
                    ("Failed to find the 'atmosphere' element in 'atom' element in config file: ") 
                    + filename);
            }
        }
        #include "AtmosphereLoadConfig.cpp.inc"
    }catch(const std::exception &exc){
        std::cerr << exc.what() << std::endl;
        abort();
    }
}
/*
*
*/
void cAtmosphereModel::RunTimeSlice(int Ma){
    int start = AtomUtils::RunStart(at);

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

    init_layer_heights();
    init_tropopause_layers();

    read_Atmosphere_Surface_Data(Ma);  // reading barthymetry data and NASA measurements

    if(!use_NASA_velocity){
        init_velocities();             // construction of zonal initial velocities from measurements
    }
    if(use_NASA_velocity)
        IC_vwt_WestEastCoast();        // horizontal velocity initial condition and temperature adjustments along coasts

    init_dynamic_pressure();           // initialization of static pressure data by dynamic pressure data to prevent zero values

    fft_gaussian_filter_3d(u,1);
    fft_gaussian_filter_3d(v,1);
    fft_gaussian_filter_3d(w,1);

    BC_SolidGround();                  // description of values inside of mountain areas

    init_temperature(Ma);              // initialization of temperature, hydrostatic pressure and density of dry air, reconstruction of potential surface values

//    goto Printout;

//    t.printArray("AGCM", im, jm, km);
//    temp_reconst.printArray_2D(jm, km);

    store_intermediate_data_2D(1.0);
    store_intermediate_data_3D(1.0);

    land_oceanFraction(0, jm, km, h);  // ratio of land to sea surface of the various time slices

    if(!use_NASA_temperature) 
        IC_t_WestEastCoast();          // horizontal temperature adjustments along west coasts and equator line due to upwelling

//    fft_gaussian_filter_3d(t, 3, direction_k);
    fft_gaussian_filter_3d(t,4);
    fft_gaussian_filter_3d(p_hydro,4);
    fft_gaussian_filter_3d(r_dry,4);

    init_water_vapour();               // initialisation of water vapour, cloud water and cloud ice based on the temperature profile

    PressureDensity();                 // recalculation of hydrostatic pressure and density of dry and humid air

    SaturationAdjustment();            // based on the initial distribution, recomputation of the cloud water and cloud ice formation in case of saturated water vapour detected

    fft_gaussian_filter_3d(c,1);
    fft_gaussian_filter_3d(cloud,1);
    fft_gaussian_filter_3d(ice,1);

    RadiationMultiLayer();             // incoming short wave and outgoing long wave radiation in dependence on atmosphere's emissivity changes initial distribution of temperature 

    switch(CategoryIceScheme){
        case 0: ZeroCategoryIceScheme();     // development of rain and snow fall, water vapour and cloud water
                break;
        case 1: OneCategoryIceScheme();      // development of rain and snow fall, water vapour and cloud water
                break;
        case 2: TwoCategoryIceScheme();      // development of rain and snow fall, water vapour, cloud water and cloud ice
                break;
        case 3: ThreeCategoryIceScheme();    // development of rain and snow fall, water vapour, cloud water, cloud ice and graupel
                break;
    }

//    MoistConvectionMidL();             // precipitation due to local mid-level moist convection by buoyancy effects
//    MoistConvectionShall();            // precipitation due to local shallow moist convection by buoyancy effects

    StandAtm_DewPoint_HumidRel();      // International Standard Atmosphere temperature profile, dew point temperature, relative humidity profile
    WaterVapourEvaporation();          // correction of surface water vapour by evaporation 
    MassStreamfunction();              // mass stream function
//    LatentSensibleHeat();              // latent and sensible heat
    init_co2(Ma);                      // greenhouse gas co2 as function of temperature

    ValueLimitationAtm();              // value limitation prevents local formation of NANs

    store_intermediate_data_2D(1.0);
    store_intermediate_data_3D(1.0);
    cout << endl << endl;
    run_3D_loop();                     // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
    cout << endl << endl;
/*
    Printout:
    run_data_atm(); 
    print_min_max_atm();
    write_file(bathymetry_name, output_path, true); // printing files for ParaView, AtmosphereDataTransfer and AtmospherePlotData
*/
    iter_cnt_3d++;
    save_data();    
    if(debug){
        fedisableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO); //not platform independent(bad, very bad, I know)
    }
    RunEnd(at, Ma, start);
    print_final_remarks();
    return;
}
/*
*
*/
void cAtmosphereModel::run_2D_loop(){
cout << endl << "      AGCM: run_2D_loop atm ................................." << endl;
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
                solveRungeKutta_2D_Atmosphere();
                ValueLimitationAtm();
                store_intermediate_data_2D(1.0);
                print_min_max_atm();
                run_data_atm();
                iter_cnt++;
            }  //  ::::::   end of velocity loop_2D: if (velocity_iter_2D > velocity_iter_max_2D)   ::::::::::::::::::::::
            if(pressure_iter_2D % checkpoint == 0){
            cout << endl << "      AGCM: write_file in run_2D_loop atm ......................." << endl;
                write_file(bathymetry_name, output_path, false);
            }
            if(iter_cnt > nm){
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        } // :::::::::::::::::::   end of pressure loop_2D: if (pressure_iter_2D > pressure_iter_max_2D)   ::::::::::
    } // ::::::::   end of 2D loop for initial surface conditions: if (switch_2D == 0)   :::::::::::::::::::::::::::::
    cout << endl << "      AGCM: run_2D_loop atm ended ................................" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::run_3D_loop(){ 
cout << endl << "      AGCM: run_3D_loop atm ..........................." << endl;
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
                residuum_loop = find_residuum_atm();
                if(fabs(residuum_loop/residuum_old - 1.0) <= eps_residuum){
                    cout << endl << "      AGCM: write_file in run_3D_loop atm, termination due to convergence .... relative error less than eps_residuum ..................." << endl;
                    cout << endl << "      residuum_loop = " << residuum_loop
                                 << "      residuum_old = " << residuum_old
                                 << "      eps_residuum = " << eps_residuum
                                 << "      absolute error = " 
                                 << fabs(residuum_old - residuum_loop)
                                 << "      relative error = " 
                                 << fabs(residuum_loop/residuum_old - 1.0) << endl;
                    write_file(bathymetry_name, output_path, false);
                    goto endloop;
                }else  cout << "      AGCM: write_file in find_residuum_atm, relative error is too high ......................." << endl << endl;
                if(iter_cnt > nm){
                    cout << "       nm = " << nm 
                        << "     .....     maximum number of iterations   nm   reached!" 
                        << endl;
                    goto endloop;
                }
            }
/*
            if(velocity_iter % 4 == 0){
//                SaturationAdjustment();
                RadiationMultiLayer(); 
                PressureDensity();
                WaterVapourEvaporation();
                TwoCategoryIceScheme(); 
                MoistConvectionShall();
                MoistConvectionMidL();
                StandAtm_DewPoint_HumidRel();
                LatentSensibleHeat(); 
            }
*/
            solveRungeKutta_3D_Atmosphere();
            residuum_old = residuum_loop;
            ValueLimitationAtm();
            store_intermediate_data_3D(1.0);
            print_min_max_atm();
            run_data_atm();
            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        } // end of velocity loop
        if(pressure_iter % checkpoint == 0){
            cout << endl << "      AGCM: write_file in run_3D_loop atm ......................." << endl;
            write_file(bathymetry_name, output_path, false);
        }
    } // end of pressure loop
    endloop:
    cout << endl << "      AGCM: run_3D_loop atm ended ..........................." << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::Run(){
    std::time_t Run_start;
    struct tm * timeinfo_begin;
    std::time(&Run_start);
    timeinfo_begin = std::localtime(&Run_start);
    std::cout << std::endl << std::endl;
    std::cout << " ... AGCM: time and date at run time begin:   " 
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
    std::cout << " ... AGCM: time and date at run time end:   " 
        << std::asctime(timeinfo_end) << std::endl;
    int Run_total = Run_end - Run_start;
    int Run_total_minutes = Run_total/60;
    int Run_total_hours = Run_total_minutes/60;
    std::cout <<  " ... AGCM: computer time needed for Ma = " << time_start 
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
void cAtmosphereModel::reset_arrays(){
    cout << endl << "      AGCM: reset_arrays" << endl;
    rad.initArray_1D(im, 1.0); // radial coordinate direction
    the.initArray_1D(jm, 2.0); // lateral coordinate direction
    phi.initArray_1D(km, 3.0); // longitudinal coordinate direction
    aux_grad_v.initArray_1D(im, 4.0); // auxilliar array

    Topography.initArray_2D(jm, km, 0.0); // topography
    Vegetation.initArray_2D(jm, km, 0.0); // vegetation via precipitation
    Precipitation.initArray_2D(jm, km, 0.0); // areas of higher precipitation
    precipitable_water.initArray_2D(jm, km, 0.0); // areas of precipitable water in the air
    precipitation_NASA.initArray_2D(jm, km, 0.0); // surface precipitation from NASA
    temperature_NASA.initArray_2D(jm, km, 0.0); // surface temperature from NASA
    velocity_v_NASA.initArray_2D(jm, km, 0.0); // surface v-velocity from NASA
    velocity_w_NASA.initArray_2D(jm, km, 0.0); // surface w-velocity from NASA
    temp_reconst.initArray_2D(jm, km, 0.0); // surface temperature from reconstruction tool
    temp_landscape.initArray_2D(jm, km, 0.0); // landscape temperature
    p_stat_landscape.initArray_2D(jm, km, 0.0); // landscape static pressure
    r_dry_landscape.initArray_2D(jm, km, 0.0); // landscape dry air density
    r_humid_landscape.initArray_2D(jm, km, 0.0); // landscape humid air density
    albedo.initArray_2D(jm, km, 0.0); // albedo = reflectivity
    epsilon_2D.initArray_2D(jm, km, 0.0); // epsilon = absorptivity
    Q_radiation.initArray_2D(jm, km, 0.0); // heat from the radiation balance in [W/m2]
    Q_Evaporation.initArray_2D(jm, km, 0.0); // evaporation heat of water by Kuttler
    Q_latent.initArray_2D(jm, km, 0.0); // latent heat from bottom values by the energy transport equation
    Q_sensible.initArray_2D(jm, km, 0.0); // sensible heat from bottom values by the energy transport equation
    Q_bottom.initArray_2D(jm, km, 0.0); // difference by Q_Radiation - Q_latent - Q_sensible
    vapour_evaporation.initArray_2D(jm, km, 0.0); // additional water vapour by evaporation
    Evaporation_Dalton.initArray_2D(jm, km, 0.0); // evaporation by Dalton in [mm/d]
    Evaporation_Penman.initArray_2D(jm, km, 0.0); // evaporation by Dalton in [mm/d]
    co2_total.initArray_2D(jm, km, 0.0); // areas of higher co2 concentration
    dew_point_temperature.initArray_2D(jm, km, 0.0); // dew point temperature
    condensation_level.initArray_2D(jm, km, 0.0); // areas of higher co2 concentration // local condensation level
    c_fix.initArray_2D(jm, km, 0.0); // local surface water vapour fixed for iterations
    i_Base.initArray_2D(jm, km, 0.0); // locations of the cloud base
    i_LFS.initArray_2D(jm, km, 0.0); // locations of the cloud top

    h.initArray(im, jm, km, 0.0); // bathymetry, depth from sea level
    t.initArray(im, jm, km, 1.0); // temperature
    u.initArray(im, jm, km, 0.0); // u-component velocity component in r-direction
    v.initArray(im, jm, km, 0.0); // v-component velocity component in theta-direction
    w.initArray(im, jm, km, 0.0); // w-component velocity component in phi-direction
    c.initArray(im, jm, km, 0.0); // water vapour
    cloud.initArray(im, jm, km, 0.0); // cloud water
    gr.initArray(im, jm, km, 0.0); // cloud graupel
    ice.initArray(im, jm, km, 0.0); // cloud ice
    cloudiness.initArray(im, jm, km, 0.0); // cloudiness, N in literature
    co2.initArray(im, jm, km, 1.0); // CO2

    tn.initArray(im, jm, km, 1.0); // temperature new
    un.initArray(im, jm, km, 0.0); // u-velocity component in r-direction new
    vn.initArray(im, jm, km, 0.0); // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, 0.0); // w-velocity component in phi-direction new
    cn.initArray(im, jm, km, 1.0); // water vapour new
    cloudn.initArray(im, jm, km, 0.0); // cloud water new
    icen.initArray(im, jm, km, 0.0); // cloud ice new
    grn.initArray(im, jm, km, 0.0); // cloud graupel new
    co2n.initArray(im, jm, km, 1.0); // CO2 new

    p_stat.initArray(im, jm, km, 0.0); // dynamic pressure
    p_hydro.initArray(im, jm, km, 1.0); // static pressure

    u_stream.initArray(im, jm, km, 0.0); // u-velocity by mass stream function
    stream.initArray(im, jm, km, 0.0); // mass streamfunction

    rhs_t.initArray(im, jm, km, 0.0); // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.0); // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.0); // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.0); // auxilliar field RHS w-velocity component
    rhs_c.initArray(im, jm, km, 0.0); // auxilliar field RHS water vapour
    rhs_cloud.initArray(im, jm, km, 0.0); // auxilliar field RHS cloud water
    rhs_ice.initArray(im, jm, km, 0.0); // auxilliar field RHS cloud ice
    rhs_g.initArray(im, jm, km, 0.0); // auxilliar field RHS cloud graupel
    rhs_co2.initArray(im, jm, km, 0.0); // auxilliar field RHS CO2

    aux_u.initArray(im, jm, km, 0.0); // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.0); // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.0); // auxilliar field w-velocity component
    aux_t.initArray(im, jm, km, 0.0); // auxilliar field t

    Q_Latent.initArray(im, jm, km, 0.0); // latent heat
    Q_Sensible.initArray(im, jm, km, 0.0); // sensible heat
    BuoyancyForce.initArray(im, jm, km, 0.0); // buoyancy force, Boussinesque approximation
    CoriolisForce.initArray(im, jm, km, 0.0); // Coriolis force terms
    PressureGradientForce.initArray(im, jm, km, 0.0); // Force caused by normal pressure gradient
    epsilon.initArray(im, jm, km, 0.0); // emissivity/ absorptivity
    radiation.initArray(im, jm, km, 0.0); // radiation
    P_rain.initArray(im, jm, km, 0.0); // rain precipitation mass rate
    P_snow.initArray(im, jm, km, 0.0); // snow precipitation mass rate
    P_graupel.initArray(im, jm, km, 0.0); // graupel precipitation mass rate
    P_conv_midl.initArray(im, jm, km, 0.0); // rain formation by mid-level cloud convection
    P_conv_shall.initArray(im, jm, km, 0.0); // rain formation by shallow cloud convection
    TempStand.initArray(im, jm, km, 0.0); // International Standard Atmosphere (ISA)
    TempDewPoint.initArray(im, jm, km, 0.0); // Dew Point Temperature
    HumidityRel.initArray(im, jm, km, 100.0); // relative humidity
    S_v.initArray(im, jm, km, 0.0); // water vapour mass rate
    S_c.initArray(im, jm, km, 0.0); // cloud water mass rate
    S_i.initArray(im, jm, km, 0.0); // cloud ice mass rate
    S_r.initArray(im, jm, km, 0.0); // rain mass rate
    S_s.initArray(im, jm, km, 0.0); // snow mass rate
    S_g.initArray(im, jm, km, 0.0); // graupel mass rate
    S_c_c.initArray(im, jm, km, 0.0); // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique
    M_u.initArray(im, jm, km, 0.0); // moist convection within the updraft
    M_d.initArray(im, jm, km, 0.0); // moist convection within the downdraft
    MC_t.initArray(im, jm, km, 0.0); // moist convection  acting on dry static energy
    MC_q.initArray(im, jm, km, 0.0); // moist convection acting on water vapour development
    MC_v.initArray(im, jm, km, 0.0); // moist convection acting on v-velocity component
    MC_w.initArray(im, jm, km, 0.0); // moist convection acting on w-velocity component
    r_dry.initArray(im, jm, km, 0.0); // density of dry air
    r_humid.initArray(im, jm, km, 0.0); // density of humid air
    g_p.initArray(im, jm, km, 0.0); // conversion cloud droplets to raindrops
    c_u.initArray(im, jm, km, 0.0); // condensation in the updraft
    e_d.initArray(im, jm, km, 0.0); // evaporation of precipitation in the downdraft
    e_l.initArray(im, jm, km, 0.0); // evaporation of cloud water in the environment
    e_p.initArray(im, jm, km, 0.0); // evaporation of precipitation below cloud base
    s.initArray(im, jm, km, 1.0); // dry static energy
    s_u.initArray(im, jm, km, 1.0); // dry static energy in the updraft
    s_d.initArray(im, jm, km, 1.0); // dry static energy in the downdraft
    u_u.initArray(im, jm, km, 0.0); // u-velocity component in the updraft
    u_d.initArray(im, jm, km, 0.0); // u-velocity component in the downdraft
    v_u.initArray(im, jm, km, 0.0); // v-velocity component in the updraft
    v_d.initArray(im, jm, km, 0.0); // v-velocity component in the downdraft
    w_u.initArray(im, jm, km, 0.0); // w-velocity component in the updraft
    w_d.initArray(im, jm, km, 0.0); // w-velocity component in the downdraft
    q_v_u.initArray(im, jm, km, 0.0); // water vapour in the updraft
    q_v_d.initArray(im, jm, km, 0.0); // water vapour in the downdraft
    q_c_u.initArray(im, jm, km, 0.0); // cloud water in the updraft
    E_u.initArray(im, jm, km, 0.0); // moist entrainment in the updraft
    D_u.initArray(im, jm, km, 0.0); // moist detrainment in the updraft
    E_d.initArray(im, jm, km, 0.0); // moist entrainment in the downdraft
    D_d.initArray(im, jm, km, 0.0); // moist detrainment in the downdraft

    for(auto &i : i_topography)
        std::fill(i.begin(), i.end(), 0);
    cout << "      AGCM: reset_arrays ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::write_file(std::string &bathymetry_name, 
    std::string &output_path, bool is_final_result){
    cout << endl << "      AGCM: write_file" << endl;
    int i_radial = 0;
//    int i_radial = 10;
/*
    double t_eff_tropo = t_tropopause_pole - t_tropopause_equator;
    if(i_radial == 0){
        for(int j = 0; j < jm; j++){
            double temp_tropopause = t_eff_tropo * parabola((double)j
                /((double)(jm-1)/2.0)) + t_tropopause_pole;   //temperature at tropopause     
            for(int k = 0; k < km; k++){
                int i_mount = i_topography[j][k];
//                int i_trop = get_tropopause_layer(j);
                int i_trop = im-1;
                double t_mount_top = (temp_tropopause - t.x[0][j][k]) *
                    (get_layer_height(i_mount) 
                   /get_layer_height(i_trop)) + t.x[0][j][k]; //temperature at mountain top
                double c_mount_top = (c_tropopause - c.x[0][j][k]) *
                    (get_layer_height(i_mount) 
                   /get_layer_height(i_trop)) + c.x[0][j][k]; //water vapour at mountain top
                t.x[i_radial][j][k] = t_mount_top;
                c.x[i_radial][j][k] = c_mount_top;
                co2.x[i_radial][j][k] = co2.x[i_mount][j][k];
            }
        }
    }
*/
    AtmosphereDataTransfer(bathymetry_name);
    AtmospherePlotData(bathymetry_name, (is_final_result ? -1 : iter_cnt-1));
    paraview_vtk_radial(bathymetry_name, Ma, i_radial, iter_cnt-1); 
    int j_longal = 62;          // Mount Everest/Himalaya
//    int j_longal = 90;          // Pacific center
    paraview_vtk_longal(bathymetry_name, j_longal, iter_cnt-1); 
    int k_zonal = 87;           // Mount Everest/Himalaya
//    int k_zonal = 180;          // Pacific center
    paraview_vtk_zonal(bathymetry_name, k_zonal, iter_cnt-1); 
    if(paraview_panorama_vts_flag){ //This function creates a large file. Use a flag to control if it is wanted.
        paraview_panorama_vts (bathymetry_name, iter_cnt-1); 
    }
    cout << endl << "      AGCM: write_file ended " << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::load_global_temperature_curve(){
    load_map_from_file(temperature_global_file, m_global_temperature_curve);
/*
    cout << "   m_global_temperature_curve" << endl;
    for(auto &printout : m_global_temperature_curve){
        cout << printout.first << " ..... " << printout.second << endl;
    }
*/
    return;
}
/*
*
*/
void cAtmosphereModel::load_equat_temperature_curve(){
    load_map_from_file(temperature_equat_file, m_equat_temperature_curve);
/*
    cout << "   m_equat_temperature_curve" << endl;
    for(auto &printout : m_equat_temperature_curve){
        cout << printout.first << " ..... " << printout.second << endl;
    }
*/
    return;
}
/*
*
*/
void cAtmosphereModel::load_pole_temperature_curve(){
    load_map_from_file(temperature_pole_file, m_pole_temperature_curve);
/*
    cout << "   m_pole_temperature_curve" << endl;
    for(auto &printout : m_pole_temperature_curve){
        cout << printout.first << " ..... " << printout.second << endl;
    }
*/
    return;
}
/*
*
*/
float cAtmosphereModel::get_temperatures_from_curve(float time, 
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
void cAtmosphereModel::init_water_vapour(){
    cout << endl << "      AGCM: init_water_vapour" << endl;
    // initial and boundary conditions of water vapour on water and land surfaces
    // a value of 0.04 stands for 40 g/kg, g water vapour per kg dry air
    // water vapour contents computed by Clausius-Clapeyron-formula
    // water vapour decreases approaching tropopause
    double E = 0.0;
    double q_rain = 0.0;
    double t_u = 0.0;
    double c_land_min = c_tropopause;
//    double c_land_eff = c_land_min - c_land;  // coefficient for the vertical decrease of water vapour over land
    double c_ocean_min = c_tropopause;
//    double c_ocean_eff = c_ocean_min - c_ocean;  // coefficient for the vertical decrease of water vapour over ocean
    c_land_red = std::vector<double>(im, c_land_min);
    c_ocean_red = std::vector<double>(im, c_ocean_min);
//    double d_i_max = (double)(im-1);
    for(int j = 0; j < jm; j++){
//        int i_trop = get_tropopause_layer(j);
//        int i_trop = im-1;
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            t_u = t.x[i_mount][j][k] * t_0;
            if(t_u >= t_0)
                E = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
            else
                E = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
            if(is_air(h, 0, j, k))
                c.x[i_mount][j][k] = c_ocean * ep * E
                    /(p_hydro.x[i_mount][j][k] - E); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
            if(is_land(h, 0, j, k))
                c.x[i_mount][j][k] = c_land * ep * E
                    /(p_hydro.x[i_mount][j][k] - E); // relativ water vapour contents on land surface reduced by factor in kg/kg
            for(int i = i_mount; i < im; i++){
/*
                double d_i = (double)i;
                c_land_red[i] = c_land - sqrt(c_land_eff * c_land_eff
                    * d_i/d_i_max);  // radially parabolic decrease with height
                c_ocean_red[i] = c_ocean - sqrt(c_ocean_eff * c_ocean_eff
                    * d_i/d_i_max);  // radially parabolic decrease with height
                c_land_red[i] = c_land * (1.0 - d_i/d_i_max);  // linear decrease with height
                c_ocean_red[i] = c_ocean * (1.0 - d_i/d_i_max);  // linear decrease with height
*/
                c_land_red[i] = c_land;  // constant, best approach
                c_ocean_red[i] = c_ocean;  // constant, best approach
                t_u = t.x[i][j][k] * t_0;
                if(t_u >= t_0)
                    E = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                else
                    E = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
                if(is_air(h, 0, j, k))
                    c.x[i][j][k] = c_ocean_red[i] * ep * E
                        /(p_hydro.x[i][j][k] - E); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                if(is_land(h, 0, j, k))
                    c.x[i][j][k] = c_land_red[i] * ep * E
                        /(p_hydro.x[i][j][k] - E); // relativ water vapour contents on land surface reduced by factor in kg/kg
/*
                cout.precision(5);
                cout.setf(ios::fixed);
                if((j == 60)&&(k == 87))  cout << endl
                    << "  WaterVapour" << endl 
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  i_mount = " << i_mount << endl
                    << "  p_hydro = " << p_hydro.x[i][j][k]
                    << "  c = " << c.x[i][j][k] * 1e3 << endl
                    << "  c_land_red = " << c_land_red[i] 
                    << "  c_ocean_red = " << c_ocean_red[i]  << endl;
*/
            } // end i
        }// end k
    }// end j
/*
*       initial distribution of cloud water and ice water
*/
    // Humility_critical approach from the textbook: Numerical Weather Prediction, Jean Coiffier, p. 211
    double x_cloud = 0.0;
    double h_T = 0.0;
    double cloud_loc_equator = 22.0;
    double cloud_loc_pole = 18.0;
//     double alfa_s = 30.0;  // controles the amount of cloud water and cloud ice
     double alfa_s = 15.0;  // controles the amount of cloud water and cloud ice
//     double alfa_s = 1.0;  // controles the amount of cloud water and cloud ice
    cloud_loc = std::vector<double>(jm, cloud_loc_pole); // radial location of cloud water maximum
    double cloud_loc_eff = cloud_loc_pole - cloud_loc_equator;  // coefficient for the zonal parabolic cloudwater extention
    double d_j_half = (double)(jm-1)/2.0;
    for(int j = 0; j < jm; j++){
        double d_j = (double)j;
        cloud_loc[j] = cloud_loc_eff 
            * parabola((double)d_j/(double)d_j_half) + cloud_loc_pole;
    }
    for(int j = 0; j < jm; j++){
        double d_i_max_cloud_loc = cloud_loc[j];
        for(int k = 0; k < km; k++){
            int i_trop = im-1;
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i <= i_trop; i++){
                double d_i = (double)i;
                x_cloud = d_i/d_i_max_cloud_loc;
                t_u = t.x[i][j][k] * t_0;
                if(t_u >= t_0)  E = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                else            E = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
                q_rain = ep * E/(p_hydro.x[i][j][k] - E); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                double del_q_ls = c.x[i][j][k] - q_rain
                    * Humility_critical(x_cloud, 1.0, 0.8);
                if(del_q_ls < 0.0)  del_q_ls = 0.0;
                cloud.x[i][j][k] = c.x[i][j][k] 
                    * (1.0 - exp(- alfa_s * del_q_ls/c.x[i][j][k]));
                double det_T_0 = t_0 - 37.0;  //  = 236.15 K = -37 °C
                if(t_u < t_0)
                    h_T = 1.0 - exp(- 0.5 * pow((t_u - t_0)/det_T_0, 2.0));
                else  h_T = 0.0;
                cloud.x[i][j][k] = cloud.x[i][j][k] * (1.0 - h_T);
                ice.x[i][j][k] = cloud.x[i][j][k] * h_T;
                gr.x[i][j][k] = cloud.x[i][j][k] * h_T;  // to be adjusted
                if(c.x[i][j][k]/q_rain < 1.0)
                    cloudiness.x[i][j][k] = pow(c.x[i][j][k]/q_rain, 0.25)
                        * (1.0 - exp(- (100.0 * cloud.x[i][j][k]/(1.0 - h_T))
                        /pow((q_rain - c.x[i][j][k]), 0.49)));
                else  cloudiness.x[i][j][k] = 1.0;  // overcast layer
/*
                cout.precision(8);
                cout.setf(ios::fixed);
                if((j == 90)&&(k == 180)) cout << endl
                    << "  cloudwater and cloudice    °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k << endl
                    << "   cloud_loc = " << cloud_loc[j]
                    << "   cloud_loc_equator = " << cloud_loc_equator 
                    << "   cloud_loc_pole = " << cloud_loc_pole << endl
                    << "   x_cloud = " << x_cloud << endl
                    << "   alfa_s = " << alfa_s << endl
                    << "   Humility_critical = " << Humility_critical(x_cloud, 1.0, 0.8) << endl
                    << "   del_q_ls = " << del_q_ls * 1e3 << endl
                    << "   h_T = " << h_T << endl
                    << "   q_rain = " << q_rain * 1e3
                    << "   c = " << c.x[i][j][k] * 1e3
                    << "   cloud = " << cloud.x[i][j][k] * 1e3
                    << "   ice = " << ice.x[i][j][k] * 1e3
                    << "   gr = " << gr.x[i][j][k] * 1e3
                    << "   cloudiness = " << cloudiness.x[i][j][k] << endl << endl;
*/
            } // end i cloud
        }// end k
    }// end j
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
                if(is_land(h, 0, j, k)){
                    c.x[0][j][k] = c.x[i_mount][j][k];
                    cloud.x[0][j][k] = cloud.x[i_mount][j][k];
                    ice.x[0][j][k] = ice.x[i_mount][j][k];
                    gr.x[0][j][k] = gr.x[i_mount][j][k];
                }
            for(int i = i_mount; i >= 0; i--){
                if(t.x[i][j][k] * t_0 <= t_00){
                    c.x[i][j][k] = 0.0;
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k] = 0.0;
                    gr.x[i][j][k] = 0.0;
                }
            }
        }
    }
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                c.x[i][j][k] = c.x[i][j][k] // cloud water and ice added to water vapour
                    + cloud.x[i][j][k] + ice.x[i][j][k] + gr.x[i][j][k];
            }
        }
    }
    cout << "      AGCM: init_water_vapour ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::init_temperature(int Ma){
    cout << endl << "      AGCM: init_temperature" << endl;
    double t_paleo_add = 0.0; 
    double t_pole_add = 0.0; 
    double t_global_mean_exp = 0.0;
    // correction of surface temperature around 180°E due to bad data around 180E
    if(is_first_time_slice()){
        int k_half = (km -1)/2;
        for(int j = 0; j < jm; j++){
            t.x[0][j][k_half] = (t.x[0][j][k_half+1] 
                + t.x[0][j][k_half-1])/2.0;
            temperature_NASA.y[j][k_half] = (temperature_NASA.y[j][k_half+1] +
                temperature_NASA.y[j][k_half-1])/2.0;
        }
    }
    // if use_earthbyte_reconstruction temperature in °C converted to non-dimensional
    if((Ma != 0)&&(use_earthbyte_reconstruction)){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                t.x[0][j][k] = (t.x[0][j][k] + t_0)/t_0; // non-dimensional, (reconstructed temperature in °C)
                temp_reconst.y[j][k] = t.x[0][j][k]; // non-dimensional
            }
        }
    }
    if((is_first_time_slice())&&(use_earthbyte_reconstruction)&&(!use_NASA_temperature)){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                t.x[0][j][k] = (temperature_NASA.y[j][k] + t_0)/t_0; // non-dimensional
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
    const char* time_slice_comment = "      time slice of Paleo-AGCM: ";
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
        << setw(13) << temperature_average_pal << " = " << setw(7) 
        << setfill(' ') << ((t_average + t_0) + t_paleo_total * t_0) - t_0 << setw(5) 
        << temperature_unit << endl

        << setw(55) << setfill('.')  // t polar
        << setiosflags(ios::left) << temperature_pole << resetiosflags(ios::left) 
        << setw(13) << temperature_pole_pal << " = " << setw(7) 
        << setfill(' ') << ((t_pole_modern + t_0) + t_pole_total * t_0) - t_0 << setw(5)
        << temperature_unit << endl

        << setw(54) << setfill('.')  // t global mean ATOM
        << setiosflags(ios::left) << temperature_mean << resetiosflags(ios::left) 
        << setw(13) << temperature_mean_pal << " = " << setw(7) 
        << setfill(' ') << t_global_mean << setw(5)
        << temperature_unit << endl

        << setw(50) << setfill('.')  // t global mean experimental
        << setiosflags(ios::left) << temperature_mean_exp << resetiosflags(ios::left) 
        << setw(13) << temperature_mean_exp_pal << " = " << setw(7) 
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
/*
    // temperatur distribution at a prescribed sun position
    // sun_position_lat = 60,    position of sun j = 120 means 30°S, j = 60 means 30°N
    // sun_position_lon = 180, position of sun k = 180 means 0° or 180° E (Greenwich, zero meridian)
    // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature (linear equation + parabola)
    if((*get_current_time() > 0)&&(sun == 1)){
        double j_par = sun_position_lat; // position of maximum temperature, sun position
        j_par = j_par + declination; // angle of sun axis, declination = 23,4°
        double j_pol = jm-1;
        double j_par_f = (double)j_par;
        double j_pol_f = (double)j_pol;
        double aa = (t_equator - t_pole)/(((j_par_f * j_par_f) 
            - (j_pol_f * j_pol_f)) - 2.0 * j_par_f * (j_par_f - j_pol_f));
        double bb = - 2.0 * aa * j_par_f;
        double cc = t_equator + aa * j_par_f * j_par_f;
        double j_d = sqrt ((cc - t_pole)/aa);
        double dd = 2.0 * aa * j_d + bb;
        double e = t_pole;
        // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature (linear equation + parabola)
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                double d_j = (double)j;
                if(d_j <= j_d){
                    t.x[0][j][k] = dd * d_j + e + t_paleo_add;
                }
                if(d_j > j_d){
                    t.x[0][j][k] = aa * d_j * d_j + bb * d_j 
                        + cc + t_paleo_add;
                }
            }
        }
        // longitudinally variable temperature distribution from west to east in parabolic form
        // communicates the impression of local sun radiation on the southern hemisphere
        double k_par = sun_position_lon;  // position of the sun at constant longitude
        double k_pol = km - 1;
        double t_360 = (t_0 + 5.0)/t_0;
        for(int j = 0; j < jm; j++){
            double jm_temp_asym = t.x[0][j][20];//transfer of zonal constant temperature into aa 1D-temperature field
            for(int k = 0; k < km; k++){
                double k_par_f = (double)k_par;
                double k_pol_f = (double)k_pol;
                double d_k = (double) k;
                aa = (jm_temp_asym - t_360)/(((k_par_f * k_par_f) 
                    - (k_pol_f * k_pol_f)) - 2. * k_par_f * 
                    (k_par_f - k_pol_f));
                bb = - 2.0 * aa * k_par_f;
                cc = jm_temp_asym + aa * k_par_f * k_par_f;
                t.x[0][j][k] = aa * d_k * d_k + bb * d_k + cc;
            }
        }
    }// temperatur distribution at aa prescribed sun position
*/
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
                t.x[0][j][k] = t_eff 
                    * parabola((double)d_j/(double)d_j_half) + t_pole;  // parabolic temperature assumption
                t.x[0][j][k] = t.x[0][j][k] + t_corr;                   // based on the parabolic temperature assumption
                if(is_land(h, 0, j, k)){
                    t.x[0][j][k] += m_model->t_land;
                }
            }
            if((use_NASA_temperature)&&(use_earthbyte_reconstruction)){ // NASA surface temperature assumed and corrected by known local values
                if(*get_current_time() == 0){
                    t.x[0][j][k] = (temperature_NASA.y[j][k] + t_0)/t_0;// initial temperature by NASA for Ma=0, non-dimensional
                }else{
                    t.x[0][j][k] = t.x[0][j][k] + t_corr;                // based on the reconstructed temperature
                    if(is_land(h, 0, j, k)){
                        t.x[0][j][k] += m_model->t_land;
                    }
                }
                if(*get_current_time() >= Ma_switch){                   // parabolic temperature distribution starting at Ma_switch
                    t.x[0][j][k] = t_eff 
                        * parabola((double)d_j/(double)d_j_half) + t_pole;
                    t.x[0][j][k] = t.x[0][j][k] + t_corr;               // based on the parabolic temperature assumption
                    if(is_land(h, 0, j, k)){
                        t.x[0][j][k] += m_model->t_land;
                    }
                }
            }
        }// for j
    }// for k
//    t.printArray("AGCM", im, jm, km);
//    temperature_NASA.printArray_2D(jm, km);
//    int i_potential = 0;
    double beta = 42.0; // in K, COSMO
    double t_u = 0.0;
//    double height_mount = 0.0;
//    double height_potential = 0.0;
    double t_potential = 0.0;
    double p_hydro_potential = 0.0;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            int i_mount = i_topography[j][k];
//            height_mount = get_layer_height(i_mount);
            if(*get_current_time() == 0){
                if(is_land(h, 0, j, k)){
                    t_u = t.x[0][j][k] * t_0;  // 3D-temperature profile on land and ocean surfaces projected to zero level in K
                    temp_landscape.y[j][k] = temperature_NASA.y[j][k] + t_0; // in K
                    // International Standard Atmosphere (ISA) with constant lapse rate, dry-adiabatic lifting
                    for(int i = 0; i < im; i++){
                        double height = get_layer_height(i);
                        t.x[i][j][k] = t_u 
                            * sqrt(1.0 - (2.0 * beta * g * height)
                            /(R_Air * t_u * t_u)); // COSMO in K
                        p_hydro.x[i][j][k] = p_0 * exp(- t_u/beta 
                            * (1.0 - sqrt(1.0 - (2.0 * beta * g * height)
                            /(R_Air * t_u * t_u))));  // reference static pressure given in hPa, by COSMO
/*
                        if(i == i_mount){
                            i_potential = i;
                            height_potential = height;
                            t_potential = t.x[i][j][k];
                            p_hydro_potential = p_hydro.x[i][j][k];
                        }
*/
                    } // end i
                    t.x[0][j][k] = t_potential/t_0 // potential temperature, landscape temperature projected to see level
                        * pow(p_0/p_hydro_potential, R_Air/cp_l); // non-dimensional
                    p_hydro.x[0][j][k] = p_hydro_potential // potential hydro static pressure, landscape pressure projected to see level
                        * pow((t.x[0][j][k] * t_0)/t_potential, cp_l/R_Air);
                    temp_reconst.y[j][k] = t.x[0][j][k] * t_0; //  in K
                    temp_landscape.y[j][k] = t.x[i_mount][j][k] - t_0; // in °C
                    p_stat_landscape.y[j][k] = p_hydro.x[i_mount][j][k]; // in hPa
/*
                    cout.precision(6);
                    if((j == 62)&&(k == 87))  cout << endl 
                        << "  printout in init_temperature (potential)  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
                        << "  Ma = " << (int)*get_current_time() << endl
                        << "  j = " << j << "  k = " << k << endl 
                        << "  i_mount = " << i_mount
                        << "  i_potential = " << i_potential
                        << "  height_mount = " << height_mount
                        << "  height_potential = " << height_potential << endl << endl

                        << "  r_air = " << r_air
                        << "  p_0 = " << p_0
                        << "  t_0 = " << t_0
                        << "  t_land = " << t_land * t_0 << endl << endl

                        << "  r_dry_potential = " << r_dry_potential
                        << "  r_dry_sl = " << r_dry.x[0][j][k] << endl << endl

                        << "  p_hydro_potential = " << p_hydro_potential
                        << "  p_hydro_sl = " << p_hydro.x[0][j][k] << endl << endl

                        << "  t_potential = " << t_potential - t_0
                        << "  t_sl = " << t.x[0][j][k] * t_0 - t_0 << endl << endl;
*/
                }else{  // is_water
                    t_u = t.x[0][j][k] * t_0;  // 3D-temperature profile on land and ocean surfaces projected to zero level in K
                    p_hydro.x[0][j][k] = 1e-2 * (r_air * R_Air * t_u);  // reference static pressure given in hPa, by gas equation
                    temp_reconst.y[j][k] = t_u; //  in K
                    temp_landscape.y[j][k] = t_u - t_0; // in °C
                }
            // Initial conditions proposed by COSMO, dry-adiabatic lifting
                for(int i = 0; i < im; i++){
                    double height = get_layer_height(i);
                    t.x[i][j][k] = temp_reconst.y[j][k] 
                        * sqrt(1.0 - (2.0 * beta * g * height)
                        /(R_Air * temp_reconst.y[j][k] 
                        * temp_reconst.y[j][k])); // COSMO in K
//                t.x[i][j][k] = temp_reconst.y[j][k] - 6.5 * height/1000.0; // International Standard Atmosphere (ISA)
                    p_hydro.x[i][j][k] = p_hydro.x[0][j][k] 
                        * exp(- temp_reconst.y[j][k]/beta 
                        * (1.0 - sqrt(1.0 - (2.0 * beta *g * height)
                        /(R_Air * temp_reconst.y[j][k] 
                        * temp_reconst.y[j][k]))));  // reference static pressure given in hPa, by COSMO
                } // end i
            }else{  // (*get_current_time() > 0)
                t_u = t.x[0][j][k] * t_0;  // in K
                p_hydro.x[0][j][k] = 1e-2 * (r_air * R_Air * t_u);  // reference static pressure given in hPa, by gas equation
//                r_dry.x[0][j][k] = 1e2 * p_hydro.x[0][j][k]
//                    /(R_Air * t_u); // in kg/m3
                temp_reconst.y[j][k] = t_u; //  in K
            // Initial conditions proposed by COSMO, dry-adiabatic lifting
                for(int i = 0; i < im; i++){
                    double height = get_layer_height(i);
                    t.x[i][j][k] = temp_reconst.y[j][k] 
                        * sqrt(1.0 - (2.0 * beta * g * height)
                        /(R_Air * temp_reconst.y[j][k] 
                        * temp_reconst.y[j][k])); // COSMO in K
//                t.x[i][j][k] = temp_reconst.y[j][k] - 6.5 * height/1000.0; // International Standard Atmosphere (ISA)
                    p_hydro.x[i][j][k] = p_hydro.x[0][j][k] 
                        * exp(- temp_reconst.y[j][k]/beta 
                        * (1.0 - sqrt(1.0 - (2.0 * beta *g * height)
                        /(R_Air * temp_reconst.y[j][k] 
                        * temp_reconst.y[j][k]))));  // reference static pressure given in hPa, by COSMO
                } // end i
                temp_landscape.y[j][k] = t.x[i_mount][j][k] - t_0; // in °C
                p_stat_landscape.y[j][k] = p_hydro.x[i_mount][j][k]; // in hPa
/*
                cout.precision(6);
                if((j == 62)&&(k == 87))  cout << endl 
                    << "  printout in init_temperature   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  j = " << j << "  k = " << k << endl 

                    << "  r_air = " << r_air
                    << "  p_0 = " << p_0
                    << "  t_0 = " << t_0
                    << "  t_land = " << t_land * t_0 << endl << endl

                    << "  i_mount = " << i_mount
                    << "  height_mount = " << get_layer_height(i_mount) << endl << endl

                    << "  r_dry_landscape = " << r_dry_landscape.y[j][k]
                    << "  r_dry = " << r_dry.x[0][j][k] << endl << endl

                    << "  p_hydro_landscape = " << p_stat_landscape.y[j][k]
                    << "  p_hydro = " << p_hydro.x[0][j][k] << endl << endl

                    << "  t_landscape = " << temp_landscape.y[j][k]
                    << "  t_reconstruct = " << temp_reconst.y[j][k] - t_0
                    << "  t = " << t.x[0][j][k] - t_0 << endl << endl;
*/
            }  // (*get_current_time() > 0)
        } // end k
    } // end j

    t_global_mean = GetMean_2D(jm, km, temp_landscape);
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            temp_reconst.y[j][k] = t.x[0][j][k] - t_0; // in °C
            for(int i = 0; i < im; i++){
                t.x[i][j][k] = t.x[i][j][k]/t_0;  // non-dimensional
            }
        }
    }
    cout << "      AGCM: init_temperature ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::init_dynamic_pressure(){
    cout << endl << "      AGCM: init_dynamic_pressure" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                p_stat.x[i][j][k] = 1e-2 * 0.5 * r_air // hPa
                    * sqrt(((u.x[i][j][k] * u.x[i][j][k]) 
                    + (v.x[i][j][k] * v.x[i][j][k]) 
                    + (w.x[i][j][k] * w.x[i][j][k]))/3.0) * u_0;
                if(is_land(h, i, j, k))
                    p_stat.x[i][j][k] = 0.0;
                if(p_stat.x[i][j][k] < 0.0) p_stat.x[i][j][k] = 0.0;
            }
        }
    }
    cout << "      AGCM: init_dynamic_pressure ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::init_topography(const string &topo_filename){
    cout << endl << "      AGCM: init_topography" << endl;
    ifstream ifile(topo_filename);
    if(! ifile.is_open()){
        std::cerr << "ERROR: could not open Name_Bathymetry_File file: " 
            <<  topo_filename << std::endl;
        abort();
    }
    double lon, lat, height;
    int j, k;
    for(j = 0; j < jm && !ifile.eof(); j++){
        for(k = 0; k < km && !ifile.eof(); k++){
            height = -999.0; // in case the height is NaN
            ifile >> lon >> lat >> height;
            if(!(height > 0.0)){
                h.x[0][j][k] = Topography.y[j][k] = 0.0;
            }else{
                Topography.y[j][k] = height;
                for(int i = 0; i < im; i++){
                    if(height > get_layer_height(i)){
                        h.x[i][j][k] = 1.0;
                    }else{
                        i_topography[j][k] = i-1;
                        break;
                    }
                }
            }
            if(ifile.fail()){
                ifile.clear();
                std::string tmp;
                std::getline(ifile, tmp);
                logger() << "bad data in topography at: " << lon << " " 
                    << lat << " " << tmp << std::endl;
            }   
        }   
    }
//  reduction and smoothing of peaks and needles in the topography
    for(int i = 0; i < im; i++){
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){

                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j-1, k))
                     &&(is_air(h, i, j+1, k)))){
                    h.x[i][j][k] = 0.0;
                    Topography.y[j][k] = h.x[i][j][k];
                }
                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j, k-1))
                     &&(is_air(h, i, j, k+1)))){
                    h.x[i][j][k] = 0.0;
                    Topography.y[j][k] = h.x[i][j][k];
                }
/*
                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j-1, k))
                     &&(is_air(h, i, j+1, k)))
                     &&((is_air(h, i, j, k-1))
                     &&(is_air(h, i, j, k+1)))){
                    h.x[i][j][k] = 0.0;
                    Topography.y[j][k] = h.x[i][j][k];
                }
*/
            }
        }
    }
    for(int i = 1; i < im; i++){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                if((is_air(h, i, j, k))&&(is_land(h, i-1, j, k))){
                    i_topography[j][k] = i;
                }
            }
        }
    }
    if(j != jm || k != km){
       std::cerr << "wrong topography file size! aborting..."<<std::endl;
        abort();
    }
    // rewriting bathymetrical data from -180° _ 0° _ +180° coordinate system to 0°- 360°
    for(int j = 0; j < jm; j++){
        move_data(Topography.y[j], km);
        move_data(i_topography[j], km);
        for(int i = 0; i < im; i++){
            move_data(h.x[i][j], km);
        }
    }
    cout << "      AGCM: init_topography ended" << endl;
    return;
}

/*
*
*/
void cAtmosphereModel::init_co2(int Ma){
    cout << endl << "      AGCM: init_co2" << endl;
    // initial and boundary conditions of CO2 content on water and land surfaces
    double t_paleo_add = 0.0; 
    if(Ma > 0){
        if((use_NASA_temperature)&&(*get_current_time() > 0))  
            t_paleo_add = 
                get_temperatures_from_curve(*get_current_time(), 
                m_global_temperature_curve)
                - get_temperatures_from_curve(*get_previous_time(), 
                m_global_temperature_curve);
        if(!use_NASA_temperature) 
//            t_paleo_add = get_temperatures_from_curve(*get_current_time(), 
//                m_global_temperature_curve) - t_average;
            t_paleo_add = 
                get_temperatures_from_curve(*get_current_time(), 
                m_global_temperature_curve)
                - get_temperatures_from_curve(*get_previous_time(), 
                m_global_temperature_curve);
    }
    co2_paleo = 3.2886 * pow ((t_paleo_add + t_average), 2) - 32.8859 *
        (t_paleo_add + t_average) + 102.2148;  // in ppm
    co2_average = 3.2886 * pow (t_average, 2) - 32.8859 * t_average + 102.2148;  // in ppm
    co2_paleo = co2_paleo - co2_average;
    cout.precision(3);
    const char* temperature_comment = "      temperature increase at paleo times: ";
    const char* temperature_gain = " t increase";
    const char* temperature_unit =  "°C ";
    const char* co_comment = "      co2 increase at paleo times: ";
    const char* co_gain = " co2 increase";
    const char* co_modern = "      mean co2 at modern times: ";
    const char* co_paleo_str = "      mean co2 at paleo times: ";
    const char* co_average_str = " co2 modern";
    const char* co_average_pal = " co2 paleo";
    const char* co_unit =  "ppm ";
    cout << setiosflags(ios::left) << setw(55) << setfill('.') 
        << temperature_comment << resetiosflags(ios::left) << setw(13) 
        << temperature_gain << " = " << setw(7) << setfill(' ') 
        << t_paleo_add << setw(5) << temperature_unit << endl
        << setiosflags (ios::left) << setw (55) << setfill ('.')
        << co_comment << resetiosflags (ios::left)         << setw (12) << co_gain << " = "
        << setw (7) << setfill (' ') << co2_paleo << setw (5) << co_unit << 
        endl << setw (55) << setfill ('.')  << setiosflags (ios::left) << co_modern
        << resetiosflags (ios::left) << setw (13) << co_average_str  << " = "
        << setw (7)  << setfill (' ') << co2_average << setw (5) << co_unit 
        << endl << setw (55) << setfill ('.')  << setiosflags (ios::left)
        << co_paleo_str << resetiosflags (ios::left) << setw (13) << co_average_pal
        << " = "  << setw (7)  << setfill (' ') << co2_average + co2_paleo
        << setw (5) << co_unit << endl;
    co2_paleo = co2_paleo/co2_0;
    co2_land = co2_land/co2_0;
    co2_ocean = co2_ocean/co2_0;
    co2_tropopause = co2_tropopause/co2_0;
    double emittancy_total = 0.423; // in W/m²
    double coeff_em = 5.6697e-8; // in W/(m² K)
    double delta_T = 0.02; // in K
    // CO2-content as initial solution
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
//            int i_mount = i_topography[j][k];
            int i_mount = 0;
            if(is_air(h, i_mount, j, k)){
                     // reciprocal formula for the temperature increase by co2, 
                     // delta_T = emittancy_total*ln(co2)/(4*coeff_em * t³)
                     // original formula by Nasif Nahle Sabag delta T(co2)
                co2.x[i_mount][j][k] = exp (4. * delta_T * coeff_em 
                     * pow((t.x[i_mount][j][k] * t_0), 3)/emittancy_total)
                     + co2_paleo + co2_ocean;
            }
            if(is_land(h, i_mount, j, k)){
                     // reciprocal formula for the temperature increase by co2, 
                     // delta_T = emittancy_total*ln(co2)/(4*coeff_em * t³)
                     // original formula by Nasif Nahle Sabag delta T(co2)
                co2.x[i_mount][j][k] = exp (4. * delta_T * coeff_em 
                     * pow((t.x[i_mount][j][k] * t_0), 3)/emittancy_total)
                     + co2_paleo + co2_land - co2_vegetation
                     * Vegetation.y[j][k]/co2_0;
            }
        }
    }
    // co2 distribution decreasing approaching tropopause, above no co2
    for(int j = 0; j < jm; j++){
//        int i_trop = get_tropopause_layer(j);
        int i_trop = im-1;
        for(int k = 0; k < km; k++){
//         int i_mount = i_topography[j][k];
            int i_mount = 0;
            for(int i = 0; i < im; i++){
                if(i <= i_trop){
                    // radial distribution approximated by a parabola
//                    double x = get_layer_height(i) 
//                       /get_layer_height(i_trop); 
                    double x = (double)i/(double)i_trop; 
                    co2.x[i][j][k] = parabola_interp(co2_tropopause, 
                        co2.x[i_mount][j][k], x); 
                }else  co2.x[i][j][k] = co2_tropopause;
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = 0; i < im; i++){
                if((is_land(h, i, j, k))&&(i < i_mount)){
                    co2.x[i][j][k] = co2.x[i_mount][j][k];
                }
                if((is_land(h, 0, j, k))&&(i == 0)){
                    co2.x[0][j][k] = co2.x[i_mount][j][k];
                }
            }
        }
    }
    cout << "      AGCM: init_co2 ended" << endl;
    return;
}
/*
*
*/
void  cAtmosphereModel::save_data(){
    cout << endl << "      AGCM: save_data" << endl;
    struct stat info;
    string path = output_path + "/bin_data/";
    if(stat(path.c_str(), &info) != 0){
        mkdir(path.c_str(), 0777);
    }
    bool last_iter = false;
    std::ostringstream ss;
    if(iter_cnt_3d == pressure_iter_max * velocity_iter_max+1){
        last_iter = true;
    }
    if(last_iter){
        ss << "_time_" << (int)(*get_current_time()) << "_iter_n";
    }else{
        ss << "_time_" << (int)(*get_current_time()) << "_iter_" << iter_cnt_3d;
    }
    std::string postfix_str = ss.str();
    Array t_t(im, jm, km, 0),  v_t(im, jm, km, 0), w_t(im, jm, km, 0), 
        m_t(im, jm, km, 0), u_t(im, jm, km, 0);
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                //m_t.x[i][j][k] = sqrt (pow (v.x[i][j][k] * u_0, 2) + pow (w.x[i][j][k] * u_0, 2));
                t_t.x[i][j][k] = t.x[i][j][k] * t_0 - t_0;
                v_t.x[i][j][k] = v.x[i][j][k] * u_0;
                w_t.x[i][j][k] = w.x[i][j][k] * u_0;
                if(last_iter){
                    if(fabs(v_t.x[i][j][k]) > 0 && !(fabs(v_t.x[0][j][k]) > 0))
                        v_t.x[0][j][k] = v_t.x[i][j][k];
                    if(fabs(w_t.x[i][j][k]) > 0 && !(fabs(w_t.x[0][j][k]) > 0))
                        w_t.x[0][j][k] = w_t.x[i][j][k];
                }
            }
        }
    }
    t_t.save(path + string("atm_t") + postfix_str, 0);
    t_t.save(path + string("atm_t") + postfix_str, 1);
    v_t.save(path + string("atm_v") + postfix_str, 0);
    w_t.save(path + string("atm_w") + postfix_str, 0);
    h.save(path + string("atm_h") + postfix_str, 0);
    Precipitation.save(path + string("atm_p") + postfix_str);
    cout << "      AGCM: save_data ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::ValueLimitationAtm(){
    cout << endl << "      ValueLimitationAtm" << endl;
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            if(Precipitation.y[j][k] >= 25.0)  Precipitation.y[j][k] = 25.0;
            if(Precipitation.y[j][k] <= 0.0)  Precipitation.y[j][k] = 0.0;
            for(int i = 0; i < im; i++){
//                if(u.x[i][j][k] >= 0.106)  u.x[i][j][k] = 0.106;
//                if(u.x[i][j][k] <= - 0.106)  u.x[i][j][k] = - 0.106;
                if(u.x[i][j][k] >= 0.5)  u.x[i][j][k] = 0.5;
                if(u.x[i][j][k] <= - 0.5)  u.x[i][j][k] = - 0.5;
                if(v.x[i][j][k] >= 0.5)  v.x[i][j][k] = 0.5;
                if(v.x[i][j][k] <= - 0.5)  v.x[i][j][k] = - 0.5;
                if(w.x[i][j][k] >= 6.25)  w.x[i][j][k] = 6.25;
                if(w.x[i][j][k] <= - 2.0)  w.x[i][j][k] = - 2.0;
                if(t.x[i][j][k] >= 1.165)  t.x[i][j][k] = 1.165;  // == 45.0 °C
                if(t.x[i][j][k] <= - 0.78)  t.x[i][j][k] = - 0.78;  // == 59.82 °C
                if(c.x[i][j][k] >= 0.04)  c.x[i][j][k] = 0.04;
                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] >= 0.02)  cloud.x[i][j][k] = 0.02;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] >= 0.01)  ice.x[i][j][k] = 0.01;
                if(ice.x[i][j][k] < 0.0)  ice.x[i][j][k] = 0.0;
                double t_u = t.x[i][j][k] * t_0;
                if(t_u <= t_00){
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k] = 0.0;
                    gr.x[i][j][k] = 0.0;
                }
//                if(P_rain.x[i][j][k] >= 10.0/8.64e4)  P_rain.x[i][j][k] = 10.0/8.64e4;
//                if(P_snow.x[i][j][k] >= 1.0/8.64e4)  P_snow.x[i][j][k] = 1.0/8.64e4;
//                if(P_graupel.x[i][j][k] >= 1.0/8.64e4)  P_graupel.x[i][j][k] = 1.0/8.64e4;
//                if(P_conv_midl.x[i][j][k] >= 10.0/8.64e4)  P_conv_midl.x[i][j][k] = 10.0/8.64e4;
//                if(P_conv_shall.x[i][j][k] >= 10.0/8.64e4)  P_conv_shall.x[i][j][k] = 10.0/8.64e4;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                if(P_graupel.x[i][j][k] < 0.0)  P_graupel.x[i][j][k] = 0.0;
                if(P_conv_midl.x[i][j][k] < 0.0)  P_conv_midl.x[i][j][k] = 0.0;
                if(P_conv_shall.x[i][j][k] < 0.0)  P_conv_shall.x[i][j][k] = 0.0;
                if(co2.x[i][j][k] >= 5.36)  co2.x[i][j][k] = 5.36;
                if(co2.x[i][j][k] <= 1.0)  co2.x[i][j][k] = 1.0;
                if(is_land(h, i, j, k))  p_stat.x[i][j][k] = 0.0;
                if(p_stat.x[i][j][k] < 0.0) p_stat.x[i][j][k] = 0.0;
            }
        }
    }
    cout << "      ValueLimitationAtm ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::BC_SolidGround(){
    cout << endl << "      AGCM: BC_SolidGround" << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(is_land(h, i, j, k)){
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;
                    un.x[i][j][k] = 0.0;
                    vn.x[i][j][k] = 0.0;
                    wn.x[i][j][k] = 0.0;
//                    t.x[i][j][k] = 1.0;  // = 273.15 K
//                    c.x[i][j][k] = 0.0; 
//                    cloud.x[i][j][k] = 0.0;
//                    ice.x[i][j][k] = 0.0;
//                    P_rain.x[i][j][k] = 0.0;
//                    P_snow.x[i][j][k] = 0.0;
//                    co2.x[i][j][k] = 1.0;  // = 280 ppm
                    p_stat.x[i][j][k] = 0.0;
                } // is_land
            } // i
        } // k
    } // j
    cout << "      AGCM: BC_SolidGround ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::store_intermediate_data_2D(float coeff){
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            vn.x[0][j][k] = coeff * v.x[0][j][k];
            wn.x[0][j][k] = coeff * w.x[0][j][k];
        }
    }
    return;
}
/*
*
*/
void cAtmosphereModel::store_intermediate_data_3D(float coeff){ 
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                un.x[i][j][k] = coeff * u.x[i][j][k];
                vn.x[i][j][k] = coeff * v.x[i][j][k];
                wn.x[i][j][k] = coeff * w.x[i][j][k];
                tn.x[i][j][k] = coeff * t.x[i][j][k];
                cn.x[i][j][k] = coeff * c.x[i][j][k];
                cloudn.x[i][j][k] = coeff * cloud.x[i][j][k];
                icen.x[i][j][k] = coeff * ice.x[i][j][k];
                grn.x[i][j][k] = coeff * gr.x[i][j][k];
                co2n.x[i][j][k] = coeff * co2.x[i][j][k];
            }
        }
    }
    return;
}
/*
*
*/
void cAtmosphereModel::adjust_temperature_IC(double** t, int jm, int km){
    for(int k=0; k < km; k++){
        for(int j=0; j < jm; j++){
//            t[j][k] = temperature_NASA.y[j][k] = (t[j][k] + t_0)/t_0;
            t[j][k] = (t[j][k] + t_0)/t_0;
        }
    }
    // correction of surface temperature around 180°E
    int k_half = (km-1)/2;
    for(int j = 0; j < jm; j++){
        t[j][k_half] = (t[j][k_half+1] + t[j][k_half-1])/2.;
        temperature_NASA.y[j][k_half] = (temperature_NASA.y[j][k_half+1] +
            temperature_NASA.y[j][k_half-1])/2.;
    }
    return;
}
/*
*
*/
void cAtmosphereModel::check_data(Array& a, Array&an, const std::string& name){
    float t_diff_min, t_diff_max, t_diff_mean, t_min, t_max;
    for(int i = 0; i < im; i++){
        t_diff_min = t_diff_max = t_diff_mean = 0; 
        t_min = t_max = a.x[i][0][0];
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                float t_diff = fabs(an.x[i][j][k]-a.x[i][j][k]);
                float tt = an.x[i][j][k];
                if(tt < t_min) t_min = tt;
                if(tt > t_max) t_max = tt;
                if(t_diff < t_diff_min) t_diff_min = t_diff;
                if(t_diff > t_diff_max) t_diff_max = t_diff;
                t_diff_mean += t_diff;
            }
        }
/*
        logger() << "layer: " << i << std::endl;
        logger() << name << " min: " << t_min << std::endl;
        logger() << name << " max: " << t_max << std::endl;
        logger() << name << " diff min: " << t_diff_min << std::endl;
        logger() << name << " diff max: " << t_diff_max << std::endl;
        logger() << name << " diff mean: " << t_diff_mean << "   " 
            << (double)t_diff_mean/(jm*km) << std::endl;
*/
    }
    return;
}
/*
*
*/
void cAtmosphereModel::check_data(){
    check_data(t,tn,"t");
    check_data(u,un,"u");
    check_data(v,vn,"v");
    check_data(w,wn,"w");
    check_data(c,cn,"c");
    check_data(cloud,cloudn,"cloud");
    check_data(ice,icen,"ice");
    check_data(co2,co2n,"c02");
    return;
}
/*
*
*/
void cAtmosphereModel::init_tropopause_layers(){
    tropopause_layers = std::vector<double>(jm, tropopause_pole);
    cout << endl << "      AGCM: init_tropopause_layers" << endl;
// Versiera di Agnesi approach, two inflection points
    int i_max = im-1;
    int j_max = jm-1;
    int j_half = j_max/2;
    double coeff_pole = 6e3; // coefficient meets tropopause_pole value of 8 km
    for(int j=j_half; j>=0; j--){
        double x = coeff_pole * (double)(j_half-j)/(double)j_half;
        tropopause_layers[j] = Agnesi(tropopause_equator, x); 
        tropopause_layers[j] = round(tropopause_layers[j] 
           /L_atm * (double)i_max); 
/*
    cout << "   j = " << j << "   x = " << x 
        << "   Agnesi = " << Agnesi(tropopause_equator, x) 
        << "   tropopause_equator = " << tropopause_equator 
        << "   tropopause_layers = " << round(tropopause_layers[j]) << endl;
*/
    }
    for(int j=j_max; j>j_half; j--){
        tropopause_layers[j] = tropopause_layers[j_max-j];
    }
/*
    if(debug){
        for(int j=0; j<jm; j++){
            std::cout << tropopause_layers[j] << " ";
        }
        std::cout << std::endl;
    }
*/
    cout << "      AGCM: init_tropopause_layers ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::IC_vwt_WestEastCoast(){
    cout << endl << "      AGCM: IC_vwt_WestEastCoast" << endl;
// initial conditions for v and w velocity components and temperature at 
// the sea surface close to east or west coasts
// transition between coast flows and open sea flows included
// upwelling along west coasts causes a temperature reduction
/** ........................... east coast extension of v- and w-velocity decrease due to upwelling ........................... **/
//    int i_max = 20; // maximun height of the coast near velocity and temperature changes
//    int i_max = 10; // maximun height of the coast near velocity and temperature changes
    int i_max = 5; // maximun height of the coast near velocity and temperature changes
    double d_i_max = get_layer_height(i_max);
    double d_i = 0.;
    double v_neg = 0.;
    double w_neg = 0.;
    double t_neg = 0.;
    int k_mid = 30; //  maximun extension of velocity management from the east coasts
    int k_beg = 1;
    for(int j = 0; j < jm; j++){
        for(int k = 1; k < km-1; k++){
            if((is_air(h, 0, j, k))&&(is_land(h, 0, j, k-1))){
                while(k >= k_beg){
                    if(is_land(h, 0, j, k+1))  break;
                    if(is_air(h, 0, j, k-2))  break;
                    if(is_air(h, 0, j, k-3))  break;
                    v_neg = - v.x[0][j][k]; // negative values for v to form global vortices
                    for(int l = 0; l <= k_mid; l++){
                        if(k+l > km-1)  break;
                        v.x[0][j][k+l] = (v.x[0][j][k+k_mid] - v_neg) 
                           /(double)k_mid * (double)l + v_neg;
                        for(int i = 1; i <= i_max; i++){ // development with height, decreasing influence
                            d_i = get_layer_height(i);
                            v.x[i][j][k+l] = (v.x[i_max][j][k+l] - v.x[0][j][k+l]) 
                               /d_i_max * d_i + v.x[0][j][k+l];
                        }
                        if(is_land(h, 0, j, k+l))  v.x[0][j][k+l] = 0.;
                    }
                    w_neg = 0.; // on land w = 0
                    for(int l = 0; l <= k_mid; l++){
                        if(k+l > km-1)  break;
                        w.x[0][j][k+l] = (w.x[0][j][k+k_mid] - w_neg) 
                           /(double)k_mid * (double)l + w_neg; // decreasing w approaching east coasts
                        for(int i = 1; i <= i_max; i++){ // development with height, decreasing influence
                            d_i = get_layer_height(i);
                            w.x[i][j][k+l] = (w.x[i_max][j][k+l] - w.x[0][j][k+l]) 
                               /d_i_max * d_i + w.x[0][j][k+l];
                        }
                        if(is_land(h, 0, j, k+l))  w.x[0][j][k+l] = 0.;
                    }
                    k_beg = k + k_mid;
                    break;
                } // end while
            } // end if
        } // end k
        k_beg = 1;
    } // end j
/** ........................... west coast extension of w-velocity decrease due to upwelling ........................... **/
    k_mid = 50; //  maximun extension of velocity management from the west coasts
    int k_mid_temp = 15; //  maximun extension of temperature management from the west coasts
//    int k_mid_temp = 10; //  maximun extension of temperature management from the west coasts
    int k_temp = 0; //  maximun extension of temperature management from the west coasts
    double vel_inter = 0.;
    double temp_inter = 0.;
//    double temp_red = 0.95;
    double temp_red = 0.98;
    for(int j = 0; j < jm; j++){
        for(int k = 1; k < km-1; k++){
            if((is_air(h, 0, j, k-1))&&(is_land(h, 0, j, k))){
                while(k >= 1){
                    if(is_air(h, 0, j, k))  break;
                    w_neg = w.x[0][j][k];
                    if(k<=k_mid) k_temp = k_mid - k;
                    else k_temp = k_mid;
                    for(int l = 0; l <= k_temp; l++){
                        vel_inter = w.x[0][j][k-l];
                        w.x[0][j][k-l] = (w.x[0][j][k-k_mid] - w_neg) 
                           /(double)k_mid * (double)l + w_neg; // decreasing w approaching west coasts
                        for(int i = 1; i <= i_max; i++){ // development with height, decreasing influence
                            d_i = get_layer_height(i);
                            w.x[i][j][k-l] = (w.x[i_max][j][k-l] - w.x[0][j][k-l]) 
                               /d_i_max * d_i + w.x[0][j][k-l];
                        }
                        if(is_land(h, 0, j, k-l))  w.x[0][j][k-l] = vel_inter;
                    }
/** ........................... west coast extension of temperature decrease ........................... **/
                    t_neg = temp_red * t.x[0][j][k];
                    if(k<=k_mid_temp) k_temp = k_mid_temp - k;
                    else k_temp = k_mid_temp;
                    for(int l = 0; l <= k_temp; l++){
                        temp_inter = t.x[0][j][k-l];
                        t.x[0][j][k-l] = (t.x[0][j][k-k_mid_temp] - t_neg) 
                           /(double)k_mid_temp * (double)l + t_neg; // decreasing temperature leaving west coasts
                            // upwelling along west coasts lowers temperature
                        for(int i = 1; i <= i_max; i++){ // development with height, decreasing influence
                            d_i = get_layer_height(i);
                            t.x[i][j][k-l] = (t.x[i_max][j][k-l] - t.x[0][j][k-l]) 
                               /d_i_max * d_i + t.x[0][j][k-l];
                        }
                        if(is_land(h, 0, j, k-l))  t.x[0][j][k-l] = temp_inter;
                    }
                    break;
                } // end while
            } // end if
        } // end k
    } // end j

/** ........................... west coast equatorial extension of temperature decrease ........................... **/

    int j_var_north = 86;
    int j_var_south = 94;
//    int k_mid_eq = 100; // maximun extension from west coasts for equatorial lowered initial temperatures
    int k_mid_eq = 70; // maximun extension from west coasts for equatorial lowered initial temperatures
    for(int j = j_var_north; j <= j_var_south; j++){  // zonal extension of equatorial lowered initial temperatures
//        for(int k = 1; k < km-1; k++){
        for(int k = 1; k < km-2; k++){
            if((is_air(h, 0, j, k-1))&&(is_land(h, 0, j, k))){
                k_temp = k_mid_eq;
                t_neg = temp_red * t.x[0][j][k];  // reduction of the original parabolic temperature approach due to upwelling with lower temperatures
                for(int l = 0; l <= k_temp; l++){
                    temp_inter = t.x[0][j][k-l];
                    t.x[0][j][k-l] = (t.x[0][j][k-k_mid_eq] - t_neg) 
                       /(double)k_mid_eq * (double)l + t_neg; // decreasing temperature leaving west coasts
                        // upwelling regions along equatorial currents beginning on west coasts lower the local temperature
                    for(int i = 1; i <= i_max; i++){
                        d_i = get_layer_height(i);
                        t.x[i][j][k-l] = (t.x[i_max][j][k-l] - t.x[0][j][k-l]) 
                           /d_i_max * d_i + t.x[0][j][k-l]; // development with height, decreasing influence
                    }
                    if(is_land(h, 0, j, k-l))  t.x[0][j][k-l] = temp_inter;
                }
            } // end if
        } // end k
    } // end j
    cout << "      AGCM: IC_vwt_WestEastCoast ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::IC_t_WestEastCoast(){
    cout << endl << "      AGCM: IC_t_WestEastCoast" << endl;
// initial conditions for temperature at 
// the sea surface close to west coasts
// transition between coast flows and open sea flows included
// upwelling along west coasts causes a temperature reduction
    double t_neg = 0.;
    int k_mid_temp = 15; //  maximun extension of temperature management from the west coasts
//    int k_mid_temp = 10; //  maximun extension of temperature management from the west coasts
    int k_temp = 0; //  maximun extension of temperature management from the west coasts
    double temp_inter = 0.;
    double temp_red = 0.95;
    for(int j = 0; j < jm; j++){
        for(int k = 1; k < km-1; k++){
            if((is_air(h, 0, j, k-1))&&(is_land(h, 0, j, k))){
                while(k >= 1){
                    if(is_air(h, 0, j, k))  break;
                    t_neg = temp_red * t.x[0][j][k];
                    if(k <= k_mid_temp) k_temp = k_mid_temp - k;
                    else k_temp = k_mid_temp;
                    for(int l = 0; l <= k_temp; l++){
                        temp_inter = t.x[0][j][k-l];
                        t.x[0][j][k-l] = (t.x[0][j][k-k_mid_temp] - t_neg) 
                           /(double)k_mid_temp * (double)l + t_neg; // decreasing temperature leaving west coasts
                            // upwelling along west coasts lowers temperature
                        if(is_land(h, 0, j, k-l))  t.x[0][j][k-l] = temp_inter;
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
            if((is_air(h, 0, j, k-1))&&(is_land(h, 0, j, k))){
                k_temp = k_mid_eq;
                t_neg = temp_red * t.x[0][j][k];  // reduction of the original parabolic temperature approach due to upwelling with lower temperatures
                for(int l = 0; l <= k_temp; l++){
                    temp_inter = t.x[0][j][k-l];
                    t.x[0][j][k-l] = (t.x[0][j][k-k_mid_eq] - t_neg) 
                       /(double)k_mid_eq * (double)l + t_neg; // decreasing temperature leaving west coasts
                        // upwelling regions along equatorial currents beginning on west coasts lower the local temperature
                    if(is_land(h, 0, j, k-l))  t.x[0][j][k-l] = temp_inter;
                }
            } // end if
        } // end k
    } // end j
    cout << "      AGCM: IC_t_WestEastCoast" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::MassStreamfunction(){
    cout << endl << "      AGCM: MassStreamfunction" << endl;
    double r0 = 6731000.0; // in m Earth radius
    double sinthe = 0.0;  // sin instead of cos function due to coordinate system
    double coeff_mass = 0.0;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            sinthe = sin(the.z[j]);
            coeff_mass = 2.0 * M_PI * r0 * sinthe;
            stream.x[im-1][j][k] = 0.0;
            for(int i = im-2; i >= 0; i--){
                double dp = get_layer_height(i+1) 
                    - get_layer_height(i);
                stream.x[i][j][k] = stream.x[i+1][j][k] 
                    + coeff_mass * r_humid.x[i][j][k] * v.x[i][j][k] 
                    * u_0 * dp;
                if(is_land(h, i, j, k))  stream.x[i][j][k] = 0.0;
/*
                cout.precision(8);
                if((j == 75) &&(k == 180)) cout << "north mass streamfunction" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                    << "   stream = " << stream.x[i][j][k] * 1.0e-10 << endl
                    << "   u = " << u.x[i][j][k] 
                    << "   v = " << v.x[i][j][k] 
                    << "   w = " << w.x[i][j][k] << endl
                    << "   costhe = " << costhe
                    << "   dp = " << dp
                    << "   coeff_mass = " << coeff_mass << endl;
                if((j == 105) &&(k == 180)) cout << "south mass streamfunction" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                    << "   stream = " << stream.x[i][j][k] * 1.0e-10 << endl
                    << "   u = " << u.x[i][j][k] 
                    << "   v = " << v.x[i][j][k] 
                    << "   w = " << w.x[i][j][k] << endl
                    << "   costhe = " << costhe
                    << "   dp = " << dp
                    << "   coeff_mass = " << coeff_mass << endl;
*/
            }  // end i
        }  // end j
    }  // end k
    for(int i = 1; i < im; i++){
        for(int k = 0; k < km; k++){
            for(int j = 1; j < jm-1; j++){
                sinthe = sin(the.z[j]);
                double coeff_mass_u = - g/(2.0 * M_PI * r0 * r0 * sinthe);
                u_stream.x[i][j][k] = coeff_mass_u 
                    * (stream.x[i][j+1][k] - stream.x[i][j-1][k])
                    /(2.0 * dthe);
                if(is_land(h, i, j, k))  u_stream.x[i][j][k] = 0.0;
            }
        }
    }

    cout << "      AGCM: MassStreamfunction ended" << endl;
    return;
}
/*
*
*/
double cAtmosphereModel::find_residuum_atm(){
    cout << endl << "      AGCM: find_residuum_atm" << endl;
    bool residum_found = false;
    double residuum_atm = 0.0;
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
                double exp_rm = 1.0/(rm + 1.0);
                double rmsinthe = rm * sinthe;
                if((is_air(h, i, j, k)&&(is_air(h, i-1, j, k)))
                    &&(is_air(h, i, j, k)&&(is_air(h, i, j-1, k))
                    &&(is_air(h, i, j+1, k)))
                    &&(is_air(h, i, j, k)&&(is_air(h, i, j, k-1))
                    &&(is_air(h, i, j, k+1)))){
                    dudr = (u.x[i+1][j][k] - u.x[i-1][j][k])/dr * exp_rm;
                    dvdthe = (v.x[i][j+1][k] - v.x[i][j-1][k])
                        /(2.0 * rm * dthe);
                    dwdphi = (w.x[i][j][k+1] - w.x[i][j][k-1])
                        /(2.0 * rmsinthe * dphi);
                    residuum_atm = fabs(dudr + dvdthe + dwdphi);
                    if(residuum_atm > maxValue){
                        maxValue = residuum_atm;
                        i_error = i;
                        j_error = j;
                        k_error = k;
                    }
/*
                    cout.precision(8);
                    if((j == 75) &&(k == 180)) cout << "residuum_atm" << endl
                        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                        << "   u = " << u.x[i][j][k] 
                        << "   v = " << v.x[i][j][k] 
                        << "   w = " << w.x[i][j][k] << endl
                        << "   dudr = " << dudr
                        << "   dvdthe = " << dvdthe
                        << "   dwdphi = " << dwdphi << endl
                        << "   residuum_old = " << residuum_old
                        << "   residuum_atm = " << residuum_atm
                        << "   maxValue = " << maxValue
                        << "   eps_residuum = " << eps_residuum << endl
                        << "   i_error = " << i_error << "   j_error = " << j_error << "   k_error = " << k_error  << endl
                        << "   relative error = " 
                        << fabs(residuum_atm/residuum_old - 1.0) << endl
                        << "   absolute error = " 
                        << fabs(residuum_old - maxValue) << endl << endl;
*/
                }  // end if loop
            }  //  end i
        }  // end k
    }  // end j
                    cout.precision(8);
                    if((residuum_old - maxValue) > 0.0){
                        cout << endl << "      AGCM: write_file in find_residuum_atm, absolute error declining ......................." << endl;
                        cout << endl << "      residuum_atm = " << maxValue
                                     << "   residuum_old = " << residuum_old
                                     << "   eps_residuum = " << eps_residuum << endl
                                     << "   i_error = " << i_error << "   j_error = " << j_error << "   k_error = " << k_error  << endl
                                     << "   relative error = " 
                                     << fabs(residuum_atm/residuum_old - 1.0) << endl
                                     << "   absolute error = " 
                                     << fabs(residuum_old - maxValue) << endl;
                        residum_found = true;
                    }
    if(residum_found == false){  
        cout << "      AGCM: write_file in find_residuum_atm, absolute error is too high ......................." << endl;
        cout << endl << "      residuum_atm = " << maxValue
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


