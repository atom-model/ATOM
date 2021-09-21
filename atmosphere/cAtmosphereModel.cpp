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

cAtmosphereModel* cAtmosphereModel::m_model = NULL;

const double cAtmosphereModel::pi180 = 180./ M_PI;      // pi180 = 57.3

const double cAtmosphereModel::the_degree = 1.;         // compares to 1° step size laterally
const double cAtmosphereModel::phi_degree = 1.;         // compares to 1° step size longitudinally

// dthe = the_degree/pi180 = 1.0/57.3 = 0.01745, 180 * .01745 = 3.141
const double cAtmosphereModel::dthe = the_degree/pi180; 
// dphi = phi_degree/pi180 = 1.0/57.3 = 0.01745, 360 * .01745 = 6.282
const double cAtmosphereModel::dphi = phi_degree/pi180;
    
const double cAtmosphereModel::dr = 0.025;    // 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
//const double cAtmosphereModel::dt = 0.00001;  // time step coincides with the CFL condition
const double cAtmosphereModel::dt = 0.00005;  // time step coincides with the CFL condition
    
const double cAtmosphereModel::the0 = 0.;             // North Pole
const double cAtmosphereModel::phi0 = 0.;             // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cAtmosphereModel::r0 = 1.; // non-dimensional
//const double cAtmosphereModel::r0 = 6731000.; // in m

cAtmosphereModel::cAtmosphereModel():
    i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    is_node_weights_initialised(false), 
    has_welcome_msg_printed(false){
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.
    //PythonStream::OverrideCout();
    if(PythonStream::is_enable()){
        backup = std::cout.rdbuf();
        std::cout.rdbuf(&ps);
    }
    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);
    // set default configuration
    SetDefaultConfig();
    coeff_mmWS = r_air/r_water_vapour; // coeff_mmWS = 1.2041/0.0094 [kg/m³/kg/m³] = 128,0827 [/]
    m_model = this;
    //  Coordinate system in form of a spherical shell
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.initArray_1D(im, 0); // radial coordinate direction
    the.initArray_1D(jm, 0); // lateral coordinate direction
    phi.initArray_1D(km, 0); // longitudinal coordinate direction
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
//    rad.printArray_1D(im);
//    the.printArray_1D(jm);
//    phi.printArray_1D(km);
    init_layer_heights();
    init_tropopause_layers();
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
        if (err){
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
    if(debug){
        feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO); //not platform independent, bad, very bad, I know
    }
    if(!is_temperature_curve_loaded()) 
        load_temperature_curve();
    reset_arrays();    
    m_current_time = m_time_list.insert(float(Ma)).first;
    struct stat info;
    if(stat(output_path.c_str(), &info) != 0){
        mkdir(output_path.c_str(), 0777);
    }
    //Prepare the temperature and precipitation data file
    string Name_SurfaceTemperature_File  = temperature_file;
    string Name_SurfaceNASATemperature_File  = temperature_file;
    string Name_SurfacePrecipitation_File = precipitation_file;
    if(Ma != 0 && use_earthbyte_reconstruction){
        Name_SurfaceTemperature_File = output_path + "/" 
            + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";
//        Name_SurfacePrecipitation_File = output_path + "/" 
//            + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";    
        velocity_v_file = output_path + "/" + std::to_string(Ma) 
            + "Ma_Reconstructed_wind_v.xyz";
        velocity_w_file = output_path + "/" + std::to_string(Ma) 
            + "Ma_Reconstructed_wind_w.xyz";
        if(stat(Name_SurfaceTemperature_File.c_str(), &info) != 0 || 
//           stat(Name_SurfacePrecipitation_File.c_str(), &info) != 0 ||
           stat(velocity_v_file.c_str(), &info) != 0 ||
           stat(velocity_w_file.c_str(), &info) != 0){
               std::string cmd_str = "python " + reconstruction_script_path 
                   + " " + std::to_string(Ma - time_step) + " " 
                   + std::to_string(Ma) + " " + output_path + " " 
                   + BathymetrySuffix + " atm";
               int ret = system(cmd_str.c_str());
               std::cout << " reconstruction script returned: " 
                   << ret << std::endl;
        } 
    }
    if(!has_welcome_msg_printed)  print_welcome_msg();
    //  initialization of the bathymetry/topography
    //  topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    init_topography(bathymetry_path + "/" + bathymetry_name);
    if(use_NASA_velocity){
        read_IC(velocity_v_file, v.x[0], jm, km);
        read_IC(velocity_w_file, w.x[0], jm, km);    
    }
    if(use_NASA_temperature) 
        read_IC(Name_SurfaceTemperature_File, t.x[0], jm, km);
    read_IC(Name_SurfacePrecipitation_File, precipitation_NASA.y, jm, km);
    read_IC(Name_SurfaceNASATemperature_File, temperature_NASA.y, jm, km);

//t.printArray(im, jm, km);

/*
*
*/
    iter_cnt = 1;
    iter_cnt_3d = -1;
    if(debug) save_data();
    iter_cnt_3d++;
    init_velocities();

    fft_gaussian_filter_3d(u,1);
    fft_gaussian_filter_3d(v,1);
    fft_gaussian_filter_3d(w,1);

    init_temperature(); // initialisation of temperature by NASA measurements or parabolic assumption, vertical extension based on moist adiabatic process

    fft_gaussian_filter_3d(t,1);

    if(use_NASA_velocity)
        IC_vwt_WestEastCoast();
    if(!use_NASA_temperature) 
        IC_t_WestEastCoast();

    init_water_vapour(); // initialisation of water vapour based on the temperature profile
//    goto Printout;
    PressureDensity(); // static pressure ==> barometric pressure without velocity gradients

    fft_gaussian_filter_3d(p_stat,1);
    fft_gaussian_filter_3d(r_dry,1);
    fft_gaussian_filter_3d(r_humid,1);

    init_dynamic_pressure(); // pressure gradient force by velocity gradients, dynamic pressure
    SaturationAdjustment(); // cloud water and cloud ice formation in case of saturated water vapour 

    fft_gaussian_filter_3d(t,1);

    fft_gaussian_filter_3d(c,1);
    fft_gaussian_filter_3d(cloud,1);
    fft_gaussian_filter_3d(ice,1);

    RadiationMultiLayer(); // incoming short wave and outgoing longwave radiation in dependence of atmosphere's emissivity ==> temperature 
    TwoCategoryIceScheme(); // development of rain and snow fall
    MoistConvection(); // precipitation due to local moist convection by buoyancy effects

    StandAtm_DewPoint_HumidRel(); // International Standard Atmospher temperature profile, dew point temperature, relative humidity profile
    WaterVapourEvaporation(); // correction of surface water vapour by evaporation 
    MassStreamfunction(); // mass stream function
    VegetationDistribution(); // areas of vegetation on land equator to boreal regions
    init_co2(); // greenhouse gas co2 as function of temperature

    ValueLimitationAtm(); // value limitation prevents local formation of NANs
    print_min_max_atm(); // min/max values of 3D variables printed
    run_data_atm(); // additional results printed

    BC_radius(); // vertical boundary conditions along computational grid
    BC_theta(); // zonal, latitudinal boundary conditions along computational grid
    BC_phi(); // longitudinal boundary conditions along computational grid

    goto Printout;
//    store_intermediate_data_2D(1.); // 2D data stored 
    store_intermediate_data_3D(1.); // 3D data stored 
//    run_2D_loop(); // iterational 2D loop to solve variables in 4-step Runge-Kutta time scheme, not applied
    cout << endl << endl;
    run_3D_loop(); // iterational 3D loop to solve variables in 4-step Runge-Kutta time scheme
    cout << endl << endl;

    Printout:
//    print_min_max_atm();
//    run_data_atm(); 
    write_file(bathymetry_name, output_path, true); // printing files for ParaView and additional data transfer

    iter_cnt_3d++;
    save_data();    
    if(debug){
        fedisableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO); //not platform independent(bad, very bad, I know)
    }
    return;
}
/*
*
*/
void cAtmosphereModel::run_2D_loop(){
    int switch_2D = 0;    
    iter_cnt = 1;
    int Ma = int(round(*get_current_time())); 
    if(switch_2D != 1){
        for(int pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++){
            for(int velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++){

                cout << endl << endl;
                cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
                cout << " 2D AGCM iterational process" << endl;
                cout << " max total iteration number nm = " << nm << endl << endl;
                cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, maximum "
                    << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl 
                    << endl << " Ma = " << Ma << "     n = " << iter_cnt << "    velocity_iter_max_2D = " << velocity_iter_max_2D
                    << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " << pressure_iter_max_2D << 
                    "    pressure_iter_2D = " << pressure_iter_2D << endl;

                BC_theta();
                BC_phi();
//                BC_SolidGround();
                solveRungeKutta_2D_Atmosphere();
                ValueLimitationAtm();
                print_min_max_atm();
                store_intermediate_data_2D();
                iter_cnt++;
            }  //  ::::::   end of velocity loop_2D: if (velocity_iter_2D > velocity_iter_max_2D)   ::::::::::::::::::::::
            computePressure_2D();
            run_data_atm();
            if(iter_cnt > nm){
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        } // :::::::::::::::::::   end of pressure loop_2D: if (pressure_iter_2D > pressure_iter_max_2D)   ::::::::::
    } // ::::::::   end of 2D loop for initial surface conditions: if (switch_2D == 0)   :::::::::::::::::::::::::::::
    return;
}
/*
*
*/
void cAtmosphereModel::run_3D_loop(){ 
cout << endl << "      run_3D_loop atm" << endl;
    iter_cnt = 1;
    iter_cnt_3d = 0;
    if(debug){
        save_data();
    }
    for(int pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++){
        for(int velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++){
            if(debug){ Array tmp = (t-1)*t_0; tmp.inspect();}
            cout << endl << endl;
            cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            cout << " 3D AGCM iterational process" << endl;
            cout << " max total iteration number nm = " << nm << endl << endl;
            cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum "
                << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl << 
                endl << " Ma = " << (int)*get_current_time() << "     n = " << iter_cnt << "    velocity_iter_max = " << velocity_iter_max << 
                "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << 
                "    pressure_iter = " << pressure_iter << endl;
            BC_radius();
            BC_theta();
            BC_phi();
//            BC_SolidGround();
/*
            if(velocity_iter % 2 == 0){
                SaturationAdjustment();
                RadiationMultiLayer(); 
//                PressureDensity();
//                WaterVapourEvaporation();
//                TwoCategoryIceScheme(); 
//                MoistConvection();
//                StandAtm_DewPoint_HumidRel();
            }
*/
            solveRungeKutta_3D_Atmosphere();
//            fft_gaussian_filter_3d(t,1);
//            fft_gaussian_filter_3d(u,1);
//            fft_gaussian_filter_3d(v,1);
//            fft_gaussian_filter_3d(w,1);
            ValueLimitationAtm();
            store_intermediate_data_3D(1.);
//            LatentHeat(); 
            print_min_max_atm();
//            VegetationDistribution();
            run_data_atm();
            if(debug)check_data(); 
            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        }//end of velocity loop
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
    }//end of pressure loop
    cout << endl << "      run_3D_loop atm ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::Run(){
    auto start_time = std::chrono::system_clock::now();
    std::time_t start_time_t = std::chrono::system_clock::to_time_t(start_time);
    logger() << "Start Time:" << std::ctime(&start_time_t) << std::endl;
    mkdir(output_path.c_str(), 0777);
    if(!has_welcome_msg_printed)  print_welcome_msg();
    for(int i = time_start; i <= time_end; i+=time_step){
        RunTimeSlice(i);
    }
    print_final_remarks();
    auto end_time = std::chrono::system_clock::now();
    std::time_t end_time_t = std::chrono::system_clock::to_time_t(end_time);
    logger() << "End Time:" << std::ctime(&end_time_t) << std::endl;
    return;
}
/*
*
*/
void cAtmosphereModel::reset_arrays(){
    rad.initArray_1D(im, 1.); // radial coordinate direction
    the.initArray_1D(jm, 2.); // lateral coordinate direction
    phi.initArray_1D(km, 3.); // longitudinal coordinate direction
    aux_grad_v.initArray_1D(im, 4.); // auxilliar array

    Topography.initArray_2D(jm, km, 0.); // topography
    Vegetation.initArray_2D(jm, km, 0.); // vegetation via precipitation
    Precipitation.initArray_2D(jm, km, 0.); // areas of higher precipitation
    precipitable_water.initArray_2D(jm, km, 0.); // areas of precipitable water in the air
    precipitation_NASA.initArray_2D(jm, km, 0.); // surface precipitation from NASA
    temperature_NASA.initArray_2D(jm, km, 0.); // surface temperature from NASA
    temp_NASA.initArray_2D(jm, km, 0.); // surface temperature from NASA for print function
    temp_reconst.initArray_2D(jm, km, 0.); // surface temperature from reconstruction tool
    albedo.initArray_2D(jm, km, 0.); // albedo = reflectivity
    epsilon_2D.initArray_2D(jm, km, 0.); // epsilon = absorptivity
    Q_radiation.initArray_2D(jm, km, 0.); // heat from the radiation balance in [W/m2]
    Q_Evaporation.initArray_2D(jm, km, 0.); // evaporation heat of water by Kuttler
//    Q_latent.initArray_2D(jm, km, 0.); // latent heat from bottom values by the energy transport equation
//    Q_sensible.initArray_2D(jm, km, 0.); // sensible heat from bottom values by the energy transport equation
    Q_bottom.initArray_2D(jm, km, 0.); // difference by Q_Radiation - Q_latent - Q_sensible
    vapour_evaporation.initArray_2D(jm, km, 0.); // additional water vapour by evaporation
    Evaporation_Dalton.initArray_2D(jm, km, 0.); // evaporation by Dalton in [mm/d]
    Evaporation_Penman.initArray_2D(jm, km, 0.); // evaporation by Dalton in [mm/d]
    co2_total.initArray_2D(jm, km, 0.); // areas of higher co2 concentration
    dew_point_temperature.initArray_2D(jm, km, 0.); // dew point temperature
    condensation_level.initArray_2D(jm, km, 0.); // areas of higher co2 concentration // local condensation level
    c_fix.initArray_2D(jm, km, 0.); // local surface water vapour fixed for iterations
    i_Base.initArray_2D(jm, km, 0.); // locations of the cloud base
    i_LFS.initArray_2D(jm, km, 0.); // locations of the cloud top

    h.initArray(im, jm, km, 0.); // bathymetry, depth from sea level
    t.initArray(im, jm, km, 1.); // temperature
    u.initArray(im, jm, km, 0.); // u-component velocity component in r-direction
    v.initArray(im, jm, km, 0.); // v-component velocity component in theta-direction
    w.initArray(im, jm, km, 0.); // w-component velocity component in phi-direction
    c.initArray(im, jm, km, 0.); // water vapour
    cloud.initArray(im, jm, km, 0.); // cloud water
    ice.initArray(im, jm, km, 0.); // cloud ice
    cloudiness.initArray(im, jm, km, 0.); // cloudiness, N in literature
    co2.initArray(im, jm, km, 1.); // CO2

    tn.initArray(im, jm, km, 1.); // temperature new
    un.initArray(im, jm, km, 0.); // u-velocity component in r-direction new
    vn.initArray(im, jm, km, 0.); // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, 0.); // w-velocity component in phi-direction new
    cn.initArray(im, jm, km, 1.); // water vapour new
    cloudn.initArray(im, jm, km, 0.); // cloud water new
    icen.initArray(im, jm, km, 0.); // cloud ice new
    co2n.initArray(im, jm, km, 1.); // CO2 new

    p_dyn.initArray(im, jm, km, 0.); // dynamic pressure
    p_dynn.initArray(im, jm, km, 0.); // dynamic pressure
    p_stat.initArray(im, jm, km, 1.); // static pressure

    u_stream.initArray(im, jm, km, 0.); // u-velocity by mass stream function
    stream.initArray(im, jm, km, 0.); // mass streamfunction

    rhs_t.initArray(im, jm, km, 0.); // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.); // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.); // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.); // auxilliar field RHS w-velocity component
    rhs_c.initArray(im, jm, km, 0.); // auxilliar field RHS water vapour
    rhs_cloud.initArray(im, jm, km, 0.); // auxilliar field RHS cloud water
    rhs_ice.initArray(im, jm, km, 0.); // auxilliar field RHS cloud ice
    rhs_co2.initArray(im, jm, km, 0.); // auxilliar field RHS CO2

    aux_u.initArray(im, jm, km, 0.); // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.); // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.); // auxilliar field w-velocity component
    aux_t.initArray(im, jm, km, 0.); // auxilliar field t

    Q_Latent.initArray(im, jm, km, 0.); // latent heat
    Q_Sensible.initArray(im, jm, km, 0.); // sensible heat
    BuoyancyForce.initArray(im, jm, km, 0.); // buoyancy force, Boussinesque approximation
    CoriolisForce.initArray(im, jm, km, 0.); // Coriolis force terms
    PressureGradientForce.initArray(im, jm, km, 0.); // Force caused by normal pressure gradient
    epsilon.initArray(im, jm, km, 0.); // emissivity/ absorptivity
    radiation.initArray(im, jm, km, 0.); // radiation
    P_rain.initArray(im, jm, km, 0.); // rain precipitation mass rate
    P_snow.initArray(im, jm, km, 0.); // snow precipitation mass rate
    P_conv.initArray(im, jm, km, 0.); // rain formation by cloud convection
    TempStand.initArray(im, jm, km, 0.); // International Standard Atmosphere (ISA)
    TempDewPoint.initArray(im, jm, km, 0.); // Dew Point Temperature
    HumidityRel.initArray(im, jm, km, 100.); // relative humidity
    S_v.initArray(im, jm, km, 0.); // water vapour mass rate due to category two ice scheme
    S_c.initArray(im, jm, km, 0.); // cloud water mass rate due to category two ice scheme
    S_i.initArray(im, jm, km, 0.); // cloud ice mass rate due to category two ice scheme
    S_r.initArray(im, jm, km, 0.); // rain mass rate due to category two ice scheme
    S_s.initArray(im, jm, km, 0.); // snow mass rate due to category two ice scheme
    S_c_c.initArray(im, jm, km, 0.); // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique
    M_u.initArray(im, jm, km, 0.); // moist convection within the updraft
    M_d.initArray(im, jm, km, 0.); // moist convection within the downdraft
    MC_t.initArray(im, jm, km, 0.); // moist convection  acting on dry static energy
    MC_q.initArray(im, jm, km, 0.); // moist convection acting on water vapour development
    MC_v.initArray(im, jm, km, 0.); // moist convection acting on v-velocity component
    MC_w.initArray(im, jm, km, 0.); // moist convection acting on w-velocity component
    r_dry.initArray(im, jm, km, 0.); // density of dry air
    r_humid.initArray(im, jm, km, 0.); // density of humid air
    g_p.initArray(im, jm, km, 0.); // conversion cloud droplets to raindrops
    c_u.initArray(im, jm, km, 0.); // condensation in the updraft
    e_d.initArray(im, jm, km, 0.); // evaporation of precipitation in the downdraft
    e_l.initArray(im, jm, km, 0.); // evaporation of cloud water in the environment
    e_p.initArray(im, jm, km, 0.); // evaporation of precipitation below cloud base
    s.initArray(im, jm, km, 1.); // dry static energy
    s_u.initArray(im, jm, km, 1.); // dry static energy in the updraft
    s_d.initArray(im, jm, km, 1.); // dry static energy in the downdraft
    u_u.initArray(im, jm, km, 0.); // u-velocity component in the updraft
    u_d.initArray(im, jm, km, 0.); // u-velocity component in the downdraft
    v_u.initArray(im, jm, km, 0.); // v-velocity component in the updraft
    v_d.initArray(im, jm, km, 0.); // v-velocity component in the downdraft
    w_u.initArray(im, jm, km, 0.); // w-velocity component in the updraft
    w_d.initArray(im, jm, km, 0.); // w-velocity component in the downdraft
    q_v_u.initArray(im, jm, km, 0.); // water vapour in the updraft
    q_v_d.initArray(im, jm, km, 0.); // water vapour in the downdraft
    q_c_u.initArray(im, jm, km, 0.); // cloud water in the updraft
    E_u.initArray(im, jm, km, 0.); // moist entrainment in the updraft
    D_u.initArray(im, jm, km, 0.); // moist detrainment in the updraft
    E_d.initArray(im, jm, km, 0.); // moist entrainment in the downdraft
    D_d.initArray(im, jm, km, 0.); // moist detrainment in the downdraft

    for(auto &i : i_topography)
        std::fill(i.begin(), i.end(), 0);
    return;
}
/*
*
*/
void cAtmosphereModel::write_file(std::string &bathymetry_name, 
    std::string &output_path, bool is_final_result){
    cout << endl << "      write_file" << endl;
    int Ma = int(round(*get_current_time()));
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
    Atmosphere_v_w_Transfer(bathymetry_name);
    Atmosphere_PlotData(bathymetry_name, (is_final_result ? -1 : iter_cnt-1));
    paraview_vtk_radial(bathymetry_name, Ma, i_radial, iter_cnt-1); 
    int j_longal = 62;          // Mount Everest/Himalaya
//    int j_longal = 22;          // Mount Everest/Himalaya
//    int j_longal = 90;          // Pacific center
    paraview_vtk_longal(bathymetry_name, j_longal, iter_cnt-1); 
    int k_zonal = 87;           // Mount Everest/Himalaya
//    int k_zonal = 328;           // Mount Everest/Himalaya
//    int k_zonal = 180;          // Pacific center
    paraview_vtk_zonal(bathymetry_name, k_zonal, iter_cnt-1); 
    if(paraview_panorama_vts_flag){ //This function creates a large file. Use a flag to control if it is wanted.
        paraview_panorama_vts (bathymetry_name, iter_cnt-1); 
    }
    cout << "      write_file ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::load_temperature_curve(){
    load_map_from_file(temperature_curve_file, m_temperature_curve);
    return;
}
/*
*
*/
float cAtmosphereModel::get_mean_temperature_from_curve(float time) const{
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
/*
*
*/
double get_pole_temperature(int Ma, int Ma_1, int Ma_2, double t_1, double t_2){
    return (t_2 - t_1)/(double) (Ma_2 - Ma_1) * (double) (Ma - Ma_1) + t_1;
}
/*
*
*/
double get_pole_temperature(int Ma, const std::map<float, float> &pole_temp_map){
    assert(pole_temp_map.size()>1);
    std::pair<int, double> up = *pole_temp_map.begin(), bottom = *++pole_temp_map.begin();
    // when Ma out of boundary
    if(Ma <= pole_temp_map.begin()->first){
        return pole_temp_map.begin()->second; 
    }else if(Ma > (--pole_temp_map.end())->first){
        return (--pole_temp_map.end())->second;
    }
    for(const auto& pair : pole_temp_map){
        if(pair.first>=Ma){
            bottom = pair;
            break;
        }else{
            up = pair;
        }
    }
    return get_pole_temperature(Ma, up.first, bottom.first, up.second, bottom.second);
}
/*
*
*/
float cAtmosphereModel::calculate_mean_temperature(const Array& temp){
    if(!is_node_weights_initialised){
        calculate_node_weights();
        is_node_weights_initialised = true;
    }
    double ret=0., weight=0.;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (t.x[0][j][k]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret += temp.x[0][j][k] * m_node_weights[j][k];
            weight += m_node_weights[j][k];
        }
    }
    return (ret/weight-1)*t_0;
}
/*
*
*/
void cAtmosphereModel::calculate_node_weights(){
    //use cosine of latitude as weights for now
    //longitudes: 0-360(km) latitudes: 90-(-90)(jm)
    double weight = 0.;
    m_node_weights.clear();
    for(int i=0; i<jm; i++){
        if(i<=90){
            weight = cos((90-i) * M_PI/180.0);
        }else{
            weight = cos((i-90) * M_PI/180.0);
        }
        m_node_weights.push_back(std::vector<double>());
        m_node_weights[i].resize(km, weight);
    }
    return;
}
/*
*
*/
void cAtmosphereModel::restrain_temperature(){
    for(int j=0;j<jm;j++){
        for(int k=0; k<km; k++){
            if(t.x[0][j][k] - 1 > 0){
                t.x[0][j][k] -= exp((t.x[0][j][k] - 1) * t_0/4 - 10) * 6/t_0;
            }
        }
    }
    double tmp_1 = get_mean_temperature_from_curve(*get_current_time());
    double tmp_2 = calculate_mean_temperature();
    double diff = tmp_2 - tmp_1;
    for(int j=0;j<jm;j++){
        for(int k=0; k<km; k++){
            t.x[0][j][k] -= diff/t_0;
            if(t.x[0][j][k] > (1 + 38/t_0)){
                t.x[0][j][k] = 1 + 38/t_0;//don't allow temperature to exceed 38 degrees.
            }
        }
    }
    return;
}
/*
*
*/
void cAtmosphereModel::init_water_vapour(){
    cout << endl << "      init_water_vapour" << endl;
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
                    /(p_stat.x[i_mount][j][k] - E); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
            if(is_land(h, 0, j, k))
                c.x[i_mount][j][k] = c_land * ep * E
                    /(p_stat.x[i_mount][j][k] - E); // relativ water vapour contents on land surface reduced by factor in kg/kg
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
                        /(p_stat.x[i][j][k] - E); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                if(is_land(h, 0, j, k))
                    c.x[i][j][k] = c_land_red[i] * ep * E
                        /(p_stat.x[i][j][k] - E); // relativ water vapour contents on land surface reduced by factor in kg/kg
/*
                cout.precision(5);
                cout.setf(ios::fixed);
                if((j == 60)&&(k == 87))  cout << endl
                    << "  WaterVapour" << endl 
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  i_mount = " << i_mount << endl
                    << "  p_stat = " << p_stat.x[i][j][k]
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
//    double alfa_s = 40.0;  // controles the amount of cloud water and cloud ice
     double alfa_s = 10.0;  // controles the amount of cloud water and cloud ice
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
//            int i_trop = get_tropopause_layer(j);
            int i_trop = im-1;
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i <= i_trop; i++){
                double d_i = (double)i;
                x_cloud = d_i/d_i_max_cloud_loc;
                t_u = t.x[i][j][k] * t_0;
                if(t_u >= t_0)
                    E = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                else
                    E = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
                if(is_air(h, i, j, k))
                    q_rain = ep * E/(p_stat.x[i][j][k] - E); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                if(is_land(h, i, j, k))
                    q_rain = ep * E/(p_stat.x[i][j][k] - E); // relativ water vapour contents on land surface reduced by factor in kg/kg
                double del_q_ls = c.x[i][j][k] - q_rain
                    * Humility_critical(x_cloud, 1.0, 0.8);
                if(del_q_ls < 0.0)  del_q_ls = 0.0;
                cloud.x[i][j][k] = c.x[i][j][k] 
                    * (1.0 - exp(- alfa_s * del_q_ls/c.x[i][j][k]));
                double det_T_0 = t_0 - 37.0;  //  = 236.15 K = -37 °C
                if(t_u < t_0)
                    h_T = 1.0 - exp(- 0.5 * pow((t_u - t_0)/det_T_0, 2.0));
                else  h_T = 0.0;
                ice.x[i][j][k] = cloud.x[i][j][k] * h_T;
                cloud.x[i][j][k] = cloud.x[i][j][k] * (1.0 - h_T);
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
                    << "   cloudiness = " << cloudiness.x[i][j][k] << endl << endl;
*/
            } // end i cloud
        }// end k
    }// end j
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if(is_land(h, i, j, k)){
                    c.x[i][j][k] = c.x[i_mount][j][k];
                    cloud.x[i][j][k] = cloud.x[i_mount][j][k];
                    ice.x[i][j][k] = ice.x[i_mount][j][k];
                }
                if(t.x[i][j][k] * t_0 <= t_00){
                    c.x[i][j][k] = 0.0;
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k] = 0.0;
                }
            }
        }
    }
/*
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                c.x[i][j][k] = c.x[i][j][k] // cloud water and ice added to water vapour
                    + cloud.x[i][j][k] + ice.x[i][j][k];
            }
        }
    }
*/
    cout << "      init_water_vapour ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::init_temperature(){
    cout << endl << "      init_temperature" << endl;
    // Lenton_etal_COPSE_time_temp, constant paleo mean temperature, added to the surface initial temperature
    // difference between mean temperature (Ma) and mean temperature (previous Ma) == t_paleo_add
    std::map<float, float> pole_temp_map;  // Stein/Rüdiger/Parish linear pole temperature (Ma) distribution
//    std::vector<double> lapse_rate(im, 0);
    int Ma = *get_current_time();
    double t_equator = pow((rad_equator/sigma),0.25)/t_0;  // surface temperature at the equator t_equator = 1.0976 compares to 28.0°C = to 299.81 K
    double t_pole = pow((rad_pole/sigma),0.25)/t_0;  // surface temperature at the poles t_pole = 0.9436 compares to -15.4°C = to 250.25 K
    double t_eff = t_pole - t_equator;
    double get_pole_temperature(int Ma, const std::map<float, float> &pole_temp_map);
    load_map_from_file(pole_temperature_file, pole_temp_map); 
    float pole_temperature = 
        (1.0 + get_pole_temperature(*get_current_time(), pole_temp_map)/t_0) 
        - t_pole;
    // Lenton_etal_COPSE_time_temp, constant paleo mean temperature, added to the surface initial temperature
    double t_paleo_add = 0.0; 
    if(!is_first_time_slice()){
        if((use_NASA_temperature)&&(*get_current_time() > 0))  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - get_mean_temperature_from_curve(*get_previous_time());
        if(!use_NASA_temperature)  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - t_average;
        t_paleo_add /= t_0; // non-dimensional 
    }
    // correction of surface temperature around 180°E due to bad data around 180E
    if(is_first_time_slice()){
        int k_half = (km -1)/2;
        for(int j = 0; j < jm; j++){
            t.x[0][j][k_half] = (t.x[0][j][k_half+1] 
                + t.x[0][j][k_half-1])/2.;
            temperature_NASA.y[j][k_half] = (temperature_NASA.y[j][k_half+1] +
                temperature_NASA.y[j][k_half-1])/2.;
        }
    }
    // if use_earthbyte_reconstruction temperature in °C converted to non-dimensional
    if(Ma != 0 && use_earthbyte_reconstruction){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                t.x[0][j][k] = (t.x[0][j][k] + t_0)/t_0;
                temp_reconst.y[j][k] = t.x[0][j][k];
            }
        }
    }
    cout.precision(3);
    const char* time_slice_comment = "      time slice of Paleo-AGCM: ";
    const char* time_slice_number = " Ma = ";
    const char* time_slice_unit = " million years";
    cout << setiosflags(ios::left) << setw(58) << setfill('.') 
        << time_slice_comment << resetiosflags(ios::left) << setw (6) 
        << fixed << setfill(' ') << time_slice_number << setw (3) << Ma 
        << setw(12) << time_slice_unit << endl;
    const char* temperature_comment = "      temperature increase at paleo times: ";
    const char* temperature_gain = " t increase";
    const char* temperature_pole_comment = "      pole temperature increase: ";
    const char* temperature_gain_pole = " t pole increase";
    const char* temperature_modern = "      mean temperature at modern times: ";
    const char* temperature_paleo = "      mean temperature at paleo times: ";
    const char* temperature_average = " t modern";
    const char* temperature_average_pal = " t paleo";
    const char* temperature_unit =  "°C ";
    cout << setiosflags(ios::left) 
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
    // pole temperature adjustment, combination of linear time dependent functions 
    // Stein/Rüdiger/Parish locally constant pole temperature
    // difference between pole temperature (Ma) and pole temperature (previous Ma)
    double d_j_half = (double)(jm-1)/2.0;
    float t_pole_diff_ocean = 0.0;
    float t_pole_diff_land = 0.0;
    //the t_pole_diff_ocean should be the difference between this time slice and the previous one
    if(!is_first_time_slice()){  // no differences between t_pole_diff_ocean and t_pole_diff_land
        t_pole_diff_ocean = get_pole_temperature(*get_current_time(), pole_temp_map)
            - get_pole_temperature(*get_previous_time(), pole_temp_map);
        t_pole_diff_land = get_pole_temperature(*get_current_time(), pole_temp_map)
            - get_pole_temperature(*get_previous_time(), pole_temp_map);
    }
    // in °C, constant local pole temperature as function of Ma for hothouse climates 
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double d_j = (double)j;
            if(!use_NASA_temperature){  // parabolic ocean surface temperature assumed
                t.x[0][j][k] = t_eff 
                    * parabola((double)d_j/(double)d_j_half) + t_pole;
                t.x[0][j][k] += t_paleo_add + m_model->t_land
                    + pole_temperature * fabs(parabola((double)d_j
                    /(double)d_j_half) + 1.0);
                if(is_land(h, 0, j, k)){
                    t.x[0][j][k] += m_model->t_land;
                }
            }else{  // if(use_NASA_temperature) ocean surface temperature based on NASA temperature distribution
                // transported for later time slices Ma by use_earthbyte_reconstruction
                if(is_land (h, 0, j, k)){  // on land a parabolic distribution assumed
                    if(*get_current_time() > 0){
                        t.x[0][j][k] += t_paleo_add + m_model->t_land
                            + t_pole_diff_land * fabs(parabola((double)d_j
                            /(double)d_j_half) + 1.0)/t_0;
                        // land surface temperature increased by mean t_paleo_add
                        // and by a zonally equator wards decreasing temperature difference
                        // Stein/Rüdiger/Parish pole temperature decreasing equator wards
                    }
                    if(*get_current_time() == 0)
                        t.x[0][j][k] = (temperature_NASA.y[j][k] + t_0)/t_0;  // initial temperature by NASA for Ma=0, non-dimensional
                }else{ // if the location is ocean
                    if(*get_current_time() > 0){
                        // ocean surface temperature increased by mean t_paleo_add
                        // and by a zonally equator wards decreasing temperature difference
                        // Stein/Rüdiger/Parish pole temperature decreasing equator wards
                        t.x[0][j][k] += t_paleo_add 
                            + t_pole_diff_ocean * fabs(parabola((double)d_j 
                           /(double)d_j_half) + 1.0)/t_0;
                    }
                    if(*get_current_time() == 0)
                        t.x[0][j][k] = (temperature_NASA.y[j][k] + t_0)/t_0;  // initial temperature by NASA for Ma=0, non-dimensional
                }
            }// else(use_NASA_temperature)
        }// for j
    }// for k
    // zonal temperature along tropopause
    //use "linear temperature decay" to generate temperature data for layers between mountain top and tropopause
    //use "mountain top temperature" for the layers below mountain top
    //use "tropopause tempeature" for the layers above tropopause
    // temperature approaching the tropopause, above constant temperature following the International Standard Atmosphere (ISA)
/*
// increases temperature considerably
    double a_s = 0.0;  // paper by Ali Varmaghani, saturated-adiabatic lifting
    double a_s0 = 3.5; // paper by Ali Varmaghani, saturated-adiabatic lifting
    double a_sm = 8.0; // paper by Ali Varmaghani, saturated-adiabatic lifting
    double height_m = get_layer_height(im-1);
*/
//    for(int i = 0; i < im; i++){
//        lapse_rate[i] = 6.5; //  == 6.5 K/km, dry-adiabatic lapes rate
//    }
    double beta = 42.0; // in K, COSMO
    double t_0_20 = t_0 + 20.0;  // reference temperature in K
    double p_SL = p_0;  // reference pressure in hPa, COSMO
    double t_u = 0.0;
    double t_0_parabola = 0.0;
    for(int k = 0; k < km; k++){
        double t_u_pole = t.x[0][0][k] * t_0;
        double t_u_equator = t.x[0][90][k] * t_0;
        for(int j = 0; j < jm; j++){
            int i_mount = i_topography[j][k];
            double height = get_layer_height(i_mount);
            t_u = t.x[0][j][k] * t_0;  // NASA-temperature at height but projected to zero level in K
            double t_eff = t_u_pole - t_u_equator;
            t_0_parabola = (t_eff * parabola((double)j
                /((double)(jm-1)/2.0)) + t_u_pole);
            if(is_land(h, 0, j, k)){
                r_dry.x[0][j][k] = 1e2 * p_SL/(R_Air * t_0_parabola); // in kg/m3, COSMO
                p_stat.x[0][j][k] = 1e-2 * (r_dry.x[0][j][k] * R_Air 
                    * t_0_parabola);  // reference static pressure given in hPa, by gas equation
                t.x[0][j][k] = t_0_parabola;
                p_stat.x[i_mount][j][k] = p_stat.x[0][j][k] // potential static pressure expected at height by barometric height formula,
                    * exp(- t_u/beta * (1.0 - sqrt(1.0 
                    - (2.0 * beta * g * height)
                    /(R_Air * t_u * t_u)))); // COSMO
                t.x[i_mount][j][k] = t_u; // COSMO
            }else{
                t.x[0][j][k] = t_u;
                r_dry.x[0][j][k] = 1e2 * p_SL/(R_Air * t.x[0][j][k]); // in kg/m3, COSMO
                p_stat.x[0][j][k] = 1e-2 * (r_dry.x[0][j][k] * R_Air 
                    * t.x[0][j][k]);  // reference static pressure given in hPa, by gas equation
            }

            cout.precision(6);
            if((j == 90)&&(k == 180))  cout << endl 
                << "  printout in init_temperature   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
                << "  j = " << j << "  k = " << k << endl 
                << "  i_mount = " << i_mount
                << "  height = " << height << endl
                << "  eta = " << cp_l/cv_l << endl
                << "  r_air = " << r_air
                << "  r_dry_zero = " << r_dry.x[0][j][k] << endl
                << "  p_0 = " << p_0 << endl
                << "  p_stat_0_ref = " << p_stat.x[0][j][k]
                << "  p_stat_t_u_height = " << p_stat.x[i_mount][j][k]
                << "  p_stat_zero = " << p_stat.x[0][j][k] << endl
                << "  t_land = " << t_land * t_0
//                << "  lapse_rate_heigt = " << lapse_rate[i_mount] << endl
                << "  t_0 = " << t_0 - t_0 << endl
                << "  t_0_parabola = " << t_0_parabola - t_0
                << "  t_u_pole = " << t_u_pole - t_0
                << "  t_u_equator = " << t_u_equator - t_0 << endl
                << "  t_0_ref = " << t_0_20 - t_0
                << "  t_u_height = " << t.x[i_mount][j][k] - t_0
                << "  t_zero = " << t.x[0][j][k] - t_0 << endl << endl;

            // International Standard Atmosphere (ISA) with constant lapse rate, dry-adiabatic lifting
            for(int i = 0; i < im; i++){
                double height = get_layer_height(i);
                p_stat.x[i][j][k] = p_stat.x[0][j][k] // potential static pressure expected at height by barometric height formula,
                    * exp(- t.x[0][j][k]/beta * (1.0 - sqrt(1.0 
                    - (2.0 * beta * g * height)
                    /(R_Air * t.x[0][j][k] * t.x[0][j][k])))); // COSMO
                t.x[i][j][k] = t.x[0][j][k] * sqrt(1.0 
                    - (2.0 * beta * g * height)
                    /(R_Air * t.x[0][j][k] * t.x[0][j][k])); // COSMO
                r_dry.x[i][j][k] = 1e2 * p_stat.x[i][j][k]/(R_Air * t.x[i][j][k]); // in kg/m3, COSMO

/*
                cout.precision(6);
                if((j == 62)&&(k == 87))  cout << endl 
                    << "  lapse_rate in init_temperature" << endl
                    << "  i = " << i
                    << "  j = " << j << "  k = " << k << endl
                    << "  height = " << height << endl
//                    << "  lapse_rate = " << lapse_rate[i]
                    << "  p_stat_height = " << p_stat.x[i][j][k]
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl;
*/

/*
            // paper by Ali Varmaghani "An Analytical Formula for Potential Water Vapor in an Atmosphere of Constant Lapse Rate", resarchgate
            // increases temperature considerably
                a_s = (a_sm - a_s0) * get_layer_height(i)/height_m + a_s0;
                t.x[i][j][k] = t.x[i_mount][j][k] 
                    - a_s * get_layer_height(i)/1000.0/t_0;

*/
/*
// virtual temperatur increase due to water vapour, Kraus, Die Atmosphäre der Erde, p. 76
// only applicable when water vapour is known
            double E = 0.0;
            double t_u = t.x[i][j][k] * t_0;
            if(t_u >= t_0)
                E = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
            else
                E = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
            if(is_air(h, 0, j, k))
                c.x[i][j][k] = c_ocean * ep * E
                    /(p_stat.x[i][j][k] - E); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
            if(is_land(h, 0, j, k))
                c.x[i][j][k] = c_land * ep * E
                    /(p_stat.x[i][j][k] - E); // relativ water vapour contents on land surface reduced by factor in kg/kg
            t.x[i][j][k] += t.x[i][j][k] 
                * (1.0 + 0.61 * c.x[i][j][k] * 1e-3)/t_0; // virtual temperatur increase due to water vapour
*/
            } // end i
        } // end k
    } // end j
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                t.x[i][j][k] = t.x[i][j][k]/t_0;
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, 0, j, k))&&(i < i_mount)){
                    t.x[i][j][k] = t.x[i_mount][j][k];
                    r_dry.x[i][j][k] = r_dry.x[i_mount][j][k];
                    p_stat.x[i][j][k] = p_stat.x[i_mount][j][k];
                }
            }
        }
    }
    cout << "      init_temperature ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::init_dynamic_pressure(){
    cout << endl << "      init_dynamic_pressure" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                p_dyn.x[i][j][k] = 0.5 * r_humid.x[i][j][k] 
                    * sqrt(((u.x[i][j][k] * u.x[i][j][k]) 
                    + (v.x[i][j][k] * v.x[i][j][k]) 
                    + (w.x[i][j][k] * w.x[i][j][k]))/3.) * u_0;
                if(is_land(h, i, j, k))
                    p_dyn.x[i][j][k] = 0.;
            }
        }
    }
    cout << "      init_dynamic_pressure ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::init_topography(const string &topo_filename){
    cout << endl << "      init_topography" << endl;
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
            height = -999; // in case the height is NaN
            ifile >> lon >> lat >> height;
            if(! (height > 0)){
                h.x[0][j][k] = Topography.y[j][k] = 0;
            }else{
                Topography.y[j][k] = height;
                for(int i = 0; i < im; i++){
                    if(height > get_layer_height(i)){
                        h.x[i][j][k] = 1;
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
            //logger() << lon << " " << lat << " " << h.x[0][j][k] << std::endl;            
        }   
    }
    for(int i = 1; i < im; i++){
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j-1, k))
                     &&(is_air(h, i, j+1, k)))){
                    h.x[i][j][k] = 0.;
                    Topography.y[j][k] = h.x[i][j][k];
//                    i_topography[j][k] = i-1;
//                    break;
                }
                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j, k-1))
                     &&(is_air(h, i, j, k+1)))){
                    h.x[i][j][k] = 0.;
                    Topography.y[j][k] = h.x[i][j][k];
//                    i_topography[j][k] = i-1;
//                    break;
                }
                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j-1, k))
                     &&(is_air(h, i, j+1, k)))
                     &&((is_air(h, i, j, k-1))
                     &&(is_air(h, i, j, k+1)))){
                    h.x[i][j][k] = 0.;
                    Topography.y[j][k] = h.x[i][j][k];
//                    i_topography[j][k] = i-1;
//                    break;
                }
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
    cout << "      init_topography ended" << endl;
    return;
}

/*
*
*/
void cAtmosphereModel::init_co2(){
    cout << endl << "      init_co2" << endl;
    // initial and boundary conditions of CO2 content on water and land surfaces
    // Lenton_etal_COPSE_time_temp, constant paleo mean temperature, added to the surface initial temperature
    // difference between mean temperature (Ma) and mean temperature (previous Ma) == t_paleo_add
    int Ma = *get_current_time();
    double get_pole_temperature(int Ma, const std::map<float, float> &pole_temp_map);
    double t_paleo_add = 0.0; 
    if(Ma > 0){
        if((use_NASA_temperature)&&(*get_current_time() > 0))  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - get_mean_temperature_from_curve(*get_previous_time());
        if(!use_NASA_temperature) 
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - t_average;
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
                     // co2_paleo taken over from Ruddiman, p 86, effect of co2 on global temperature
                co2.x[i_mount][j][k] = exp (4. * delta_T * coeff_em 
                     * pow((t.x[i_mount][j][k] * t_0), 3)/emittancy_total)
                     + co2_paleo + co2_ocean;
            }
            if(is_land(h, i_mount, j, k)){
                     // reciprocal formula for the temperature increase by co2, 
                     // delta_T = emittancy_total*ln(co2)/(4*coeff_em * t³)
                     // original formula by Nasif Nahle Sabag delta T(co2)
                     // co2_paleo taken over from Ruddiman, p 86, effect of co2 on global temperature
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
    cout << "      init_co2 ended" << endl;
    return;
}
/*
*
*/
void  cAtmosphereModel::save_data(){
    cout << endl << "      save_data" << endl;
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
    cout << "      save_data ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::BC_SolidGround(){
    cout << endl << "      BC_SolidGround" << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = im-2; i >= 0; i--){
                if(is_land(h, i, j, k)){
                    u.x[i][j][k] = 0.;
                    v.x[i][j][k] = 0.;
                    w.x[i][j][k] = 0.;
//                    t.x[i][j][k] = 1.;  // = 273.15 K
//                    c.x[i][j][k] = 0.; 
//                    cloud.x[i][j][k] = 0.;
//                    ice.x[i][j][k] = 0.;
//                    P_rain.x[i][j][k] = 0.;
//                    P_snow.x[i][j][k] = 0.;
//                    co2.x[i][j][k] = 1.;  // = 280 ppm
//                    p_dyn.x[i][j][k] = 0.;
                } // is_land
            } // i
        } // k
    } // j
    cout << "      BC_SolidGround ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::VegetationDistribution(){
    // vegetation as function of precipitation, timberline and temperature
    cout << endl << "      VegetationDistribution" << endl;
    double max_Precipitation = 0.0;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(Precipitation.y[j][k] > max_Precipitation){
                max_Precipitation = Precipitation.y[j][k];
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            if(max_Precipitation > 0 && is_land(h, 0, j, k) 
                && !(get_layer_height(i_mount) > 4400.0) // above the vegetation line of 4.4km no vegetation
                && !((t.x[0][j][k] * t_0 - t_0) < - 40.0)){ //  vegetation >= -40°C
                Vegetation.y[j][k] = Precipitation.y[j][k] 
                    /max_Precipitation; // actual vegetation areas
            }else{
                Vegetation.y[j][k] = 0.;
            }
        }
    }
    cout << "      VegetationDistribution ended" << endl;
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
            p_dynn.x[0][j][k] = coeff * p_dyn.x[0][j][k];
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
                p_dynn.x[i][j][k] = coeff * p_dyn.x[i][j][k];
                tn.x[i][j][k] = coeff * t.x[i][j][k];
                cn.x[i][j][k] = coeff * c.x[i][j][k];
                cloudn.x[i][j][k] = coeff * cloud.x[i][j][k];
                icen.x[i][j][k] = coeff * ice.x[i][j][k];
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
    cout << endl << "      init_tropopause_layers" << endl;
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
    cout << "      init_tropopause_layers ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::IC_vwt_WestEastCoast(){
    cout << endl << "      IC_vwt_WestEastCoast" << endl;
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
    cout << "      IC_vwt_WestEastCoast ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::IC_t_WestEastCoast(){
    cout << endl << "      IC_t_WestEastCoast" << endl;
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
    cout << "      IC_t_WestEastCoast ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::MassStreamfunction(){
    cout << endl << "      MassStreamfunction" << endl;
    double r0 = 6731000.; // in m Earth radius
    double costhe = 0.;
    double coeff_mass = 0.;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            costhe = cos(the.z[j]);
            coeff_mass = 2. * M_PI * r0 * costhe/g;
            aux_grad_v.z[0] = p_stat.x[0][j][k];
            aux_v.x[0][j][k] = v.x[0][j][k] * aux_grad_v.z[0];
            for(int i = 1; i < im; i++){
                double height = get_layer_height(i-1);
                double step_meter = get_layer_height(i) - height;
                double dp = r_humid.x[i][j][k] * g * step_meter;
                if(is_land(h, i, j, k)){
                    aux_grad_v.z[i] = 0.;
                }else{
                    aux_grad_v.z[i] = p_stat.x[i][j][k];
                }
                aux_v.x[i][j][k] = v.x[i][j][k] * trapezoidal(0, i, dp, aux_grad_v);
                stream.x[i][j][k] = coeff_mass * aux_v.x[i][j][k];  // in (kg/s)
//                stream.x[i][j][k] = 1e-10 * stream.x[i][j][k];  // in (kg/s)/1e10
/*
                cout.precision(8);
                if((j == 75) &&(k == 180)) cout << "north mass streamfunction" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                    << "   aux_v = " << aux_v.x[i][j][k] 
                    << "   stream = " << stream.x[i][j][k] << endl
                    << "   u = " << u.x[i][j][k] 
                    << "   v = " << v.x[i][j][k] 
                    << "   w = " << w.x[i][j][k] << endl
                    << "   costhe = " << costhe
                    << "   step_meter = " << step_meter
                    << "   dp = " << dp
                    << "   coeff_mass = " << coeff_mass << endl;
                if((j == 105) &&(k == 180)) cout << "south mass streamfunction" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k  << endl
                    << "   aux_v = " << aux_v.x[i][j][k] 
                    << "   stream = " << stream.x[i][j][k] << endl
                    << "   u = " << u.x[i][j][k] 
                    << "   v = " << v.x[i][j][k] 
                    << "   w = " << w.x[i][j][k] << endl
                    << "   costhe = " << costhe
                    << "   step_meter = " << step_meter
                    << "   dp = " << dp
                    << "   coeff_mass = " << coeff_mass << endl;
*/
            }  // end i
        }  // end j
    }  // end k
    for(int i = 1; i < im; i++){
//        double rm = rad.z[i];
        double rm = get_layer_height(i);
        for(int k = 0; k < km; k++){
            for(int j = 1; j < jm-1; j++){
                costhe = cos(the.z[j]);
                double coeff_mass_u = - g/(2. * M_PI * r0 * r0 * costhe);
                u_stream.x[i][j][k] = coeff_mass_u 
                    * (stream.x[i][j+1][k] - stream.x[i][j-1][k])
                    /(2. * rm * dthe);
            }
        }
    }
    cout << "      MassStreamfunction ended" << endl;
    return;
}


