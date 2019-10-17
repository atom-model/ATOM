#include "cAtmosphereModel.h"

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

#include "Accuracy_Atm.h"
#include "RHS_Atm.h"
#include "RungeKutta_Atm.h"
#include "Utils.h"
#include "Config.h"
#include "AtomMath.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

cAtmosphereModel* cAtmosphereModel::m_model = NULL;

const double cAtmosphereModel::pi180 = 180./ M_PI;      // pi180 = 57.3

const double cAtmosphereModel::the_degree = 1.;         // compares to 1° step size laterally
const double cAtmosphereModel::phi_degree = 1.;         // compares to 1° step size longitudinally

// dthe = the_degree / pi180 = 1.0 / 57.3 = 0.01745, 180 * .01745 = 3.141
const double cAtmosphereModel::dthe = the_degree / pi180; 
// dphi = phi_degree / pi180 = 1.0 / 57.3 = 0.01745, 360 * .01745 = 6.282
const double cAtmosphereModel::dphi = phi_degree / pi180;
    
const double cAtmosphereModel::dr = 0.025;    // 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
const double cAtmosphereModel::dt = 0.00001;  // time step coincides with the CFL condition
    
const double cAtmosphereModel::the0 = 0.;             // North Pole
const double cAtmosphereModel::phi0 = 0.;             // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cAtmosphereModel::r0 = 1.; 

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
    coeff_mmWS = r_air / r_water_vapour; // coeff_mmWS = 1.2041 / 0.0094 [kg/m³ / kg/m³] = 128,0827 [/]
    m_model = this;
    //  Coordinate system in form of a spherical shell
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.initArray_1D(im, 0); // radial coordinate direction
    the.initArray_1D(jm, 0); // lateral coordinate direction
    phi.initArray_1D(km, 0); // longitudinal coordinate direction
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
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
//        Name_SurfacePrecipitation_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";    
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
    if(!has_welcome_msg_printed)
        print_welcome_msg();
    //  initialization of the bathymetry/topography
    //  topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    init_topography(bathymetry_path + "/" + bathymetry_name);
    if(use_NASA_velocity){
        read_IC(velocity_v_file, v.x[0], jm, km);
        read_IC(velocity_w_file, w.x[0], jm, km);    
    }
//    read_IC(Name_SurfaceTemperature_File, t.x[0], jm, km);
    read_IC(Name_SurfacePrecipitation_File, precipitation_NASA.y, jm, km);
    read_IC(Name_SurfaceNASATemperature_File, temperature_NASA.y, jm, km);
    iter_cnt_3d = -1;
    if(debug) save_data();
    iter_cnt_3d++;
    init_velocities();
//    near_wall_values();
//    adjust_temperature_IC(t.x[0], jm, km);
    init_temperature();
    init_water_vapour();
//    goto Printout;
    store_intermediate_data_3D();
    BC_Pressure();
    init_co2();
//    goto Printout;
    Ice_Water_Saturation_Adjustment();
    BC_Radiation_multi_layer(); 
//    goto Printout;
    store_intermediate_data_2D();
    store_intermediate_data_3D();
    run_2D_loop();
    cout << endl << endl;
    run_3D_loop();
    cout << endl << endl;
//    restrain_temperature();
//    Printout:
//    write_file(bathymetry_name, output_path, true);
    iter_cnt_3d++;
    save_data();    
    if(debug){
        fedisableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO); //not platform independent(bad, very bad, I know)
    }
}



void cAtmosphereModel::Run(){
    auto start_time = std::chrono::system_clock::now();
    std::time_t start_time_t = std::chrono::system_clock::to_time_t(start_time);
    logger() << "Start Time:" << std::ctime(&start_time_t) << std::endl;
    mkdir(output_path.c_str(), 0777);
    print_welcome_msg();
    for(int i = time_start; i <= time_end; i+=time_step){
        RunTimeSlice(i);
    }
    print_final_remarks();
    auto end_time = std::chrono::system_clock::now();
    std::time_t end_time_t = std::chrono::system_clock::to_time_t(end_time);
    logger() << "End Time:" << std::ctime(&end_time_t) << std::endl;
}


void cAtmosphereModel::reset_arrays(){
    Topography.initArray_2D(jm, km, 0.); // topography
    Vegetation.initArray_2D(jm, km, 0.); // vegetation via precipitation
    Precipitation.initArray_2D(jm, km, 0.); // areas of higher precipitation
    precipitable_water.initArray_2D(jm, km, 0.); // areas of precipitable water in the air
    precipitation_NASA.initArray_2D(jm, km, 0.); // surface precipitation from NASA
    radiation_surface.initArray_2D(jm, km, 0.); // direct sun radiation, short wave
    temperature_NASA.initArray_2D(jm, km, 0.); // surface temperature from NASA
    temp_NASA.initArray_2D(jm, km, 0.); // surface temperature from NASA for print function
    albedo.initArray_2D(jm, km, 0.); // albedo = reflectivity
    epsilon.initArray_2D(jm, km, 0.); // epsilon = absorptivity
    Q_radiation.initArray_2D(jm, km, 0.); // heat from the radiation balance in [W/m2]
    Q_Evaporation.initArray_2D(jm, km, 0.); // evaporation heat of water by Kuttler
    Q_latent.initArray_2D(jm, km, 0.); // latent heat from bottom values by the energy transport equation
    Q_sensible.initArray_2D(jm, km, 0.); // sensible heat from bottom values by the energy transport equation
    Q_bottom.initArray_2D(jm, km, 0.); // difference by Q_Radiation - Q_latent - Q_sensible
    vapour_evaporation.initArray_2D(jm, km, 0.); // additional water vapour by evaporation
    Evaporation_Dalton.initArray_2D(jm, km, 0.); // evaporation by Dalton in [mm/d]
    co2_total.initArray_2D(jm, km, 0.); // areas of higher co2 concentration
    dew_point_temperature.initArray_2D(jm, km, 0.); // dew point temperature
    condensation_level.initArray_2D(jm, km, 0.); // areas of higher co2 concentration // local condensation level
    c_fix.initArray_2D(jm, km, 0.); // local surface water vapour fixed for iterations
    h.initArray(im, jm, km, 0.); // bathymetry, depth from sea level
    t.initArray(im, jm, km, 1.); // temperature
    u.initArray(im, jm, km, 0.); // u-component velocity component in r-direction
    v.initArray(im, jm, km, 0.); // v-component velocity component in theta-direction
    w.initArray(im, jm, km, 0.); // w-component velocity component in phi-direction
    c.initArray(im, jm, km, 1.); // water vapour
    cloud.initArray(im, jm, km, 0.); // cloud water
    ice.initArray(im, jm, km, 0.); // cloud ice
    co2.initArray(im, jm, km, 1.); // CO2
    tn.initArray(im, jm, km, 1.); // temperature new
    un.initArray(im, jm, km, 0.); // u-velocity component in r-direction new
    vn.initArray(im, jm, km, 0.); // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, 0.); // w-velocity component in phi-direction new
    cn.initArray(im, jm, km, 1.); // water vapour new
    cloudn.initArray(im, jm, km, 0.); // cloud water new
    icen.initArray(im, jm, km, 0.); // cloud ice new
    co2n.initArray(im, jm, km, 1.); // CO2 new
    p_dyn.initArray(im, jm, km, 1.); // dynamic pressure
    p_dynn.initArray(im, jm, km, 1.); // dynamic pressure
    p_stat.initArray(im, jm, km, 1.); // static pressure
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
    Q_Latent.initArray(im, jm, km, 0.); // latent heat
    Q_Sensible.initArray(im, jm, km, 0.); // sensible heat
    BuoyancyForce.initArray(im, jm, km, 0.); // buoyancy force, Boussinesque approximation
    epsilon_3D.initArray(im, jm, km, 0.); // emissivity/ absorptivity
    radiation_3D.initArray(im, jm, km, 0.); // radiation
    P_rain.initArray(im, jm, km, 0.); // rain precipitation mass rate
    P_snow.initArray(im, jm, km, 0.); // snow precipitation mass rate
    P_conv.initArray(im, jm, km, 0.); // rain formation by cloud convection
    S_v.initArray(im, jm, km, 0.); // water vapour mass rate due to category two ice scheme
    S_c.initArray(im, jm, km, 0.); // cloud water mass rate due to category two ice scheme
    S_i.initArray(im, jm, km, 0.); // cloud ice mass rate due to category two ice scheme
    S_r.initArray(im, jm, km, 0.); // rain mass rate due to category two ice scheme
    S_s.initArray(im, jm, km, 0.); // snow mass rate due to category two ice scheme
    S_c_c.initArray(im, jm, km, 0.); // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique
    M_u.initArray(im, jm, km, 0.); // moist convection within the updraft
    M_d.initArray(im, jm, km, 0.); // moist convection within the downdraft
    MC_s.initArray(im, jm, km, 0.); // moist convection  acting on dry static energy
    MC_q.initArray(im, jm, km, 0.); // moist convection acting on water vapour development
    MC_v.initArray(im, jm, km, 0.); // moist convection acting on v-velocity component
    MC_w.initArray(im, jm, km, 0.); // moist convection acting on w-velocity component
    for(auto &i : i_topography)
        std::fill(i.begin(), i.end(), 0);
}


void cAtmosphereModel::write_file(std::string &bathymetry_name, 
    std::string &output_path, bool is_final_result){
    int Ma = int(round(*get_current_time()));
    int i_radial = 0;
    double t_eff_tropo = t_tropopause_pole - t_tropopause;
    if(i_radial == 0){
        for(int j = 0; j < jm; j++){
            double temp_tropopause = t_eff_tropo * parabola((double)j
                /((double)(jm-1)/2.0)) + t_tropopause_pole;   //temperature at tropopause     
            for(int k = 0; k < km; k++){
                int i_mount = i_topography[j][k];
                int i_trop = get_tropopause_layer(j);
                double t_mount_top = (temp_tropopause - t.x[0][j][k]) *
                    (get_layer_height(i_mount) 
                    / get_layer_height(i_trop)) + t.x[0][j][k]; //temperature at mountain top
                double c_mount_top = (c_tropopause - c.x[0][j][k]) *
                    (get_layer_height(i_mount) 
                    / get_layer_height(i_trop)) + c.x[0][j][k]; //temperature at mountain top
                t.x[i_radial][j][k] = t_mount_top;
                c.x[i_radial][j][k] = c_mount_top;
                co2.x[i_radial][j][k] = co2.x[i_mount][j][k];
            }
        }
    }
    paraview_vtk_radial(bathymetry_name, Ma, i_radial, iter_cnt-1); 
    int j_longal = 62;          // Mount Everest/Himalaya
    paraview_vtk_longal(bathymetry_name, j_longal, iter_cnt-1); 
    int k_zonal = 87;           // Mount Everest/Himalaya
    paraview_vtk_zonal(bathymetry_name, k_zonal, iter_cnt-1); 
    if(paraview_panorama_vts_flag){ //This function creates a large file. Use a flag to control if it is wanted.
        paraview_panorama_vts (bathymetry_name, iter_cnt-1); 
    }
    Value_Limitation_Atm();
    Atmosphere_v_w_Transfer(bathymetry_name);
    Atmosphere_PlotData(bathymetry_name, (is_final_result ? -1 : iter_cnt-1));
}

void cAtmosphereModel::run_2D_loop(){
    int switch_2D = 0;    
    iter_cnt = 1;
//    int Ma = int(round(*get_current_time())); 
    if(switch_2D != 1){
        for(int pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++){
            for(int velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++){
/*
                cout << endl << endl;
                cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
                cout << " 2D AGCM iterational process" << endl;
                cout << " max total iteration number nm = " << nm << endl << endl;
                cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, maximum "
                    << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl 
                    << endl << " Ma = " << Ma << "     n = " << iter_cnt << "    velocity_iter_max_2D = " << velocity_iter_max_2D
                    << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " << pressure_iter_max_2D << 
                    "    pressure_iter_2D = " << pressure_iter_2D << endl;
*/
                BC_theta();
                BC_phi();
                Value_Limitation_Atm();
                BC_SolidGround();
                solveRungeKutta_2D_Atmosphere();
                store_intermediate_data_2D();
                iter_cnt++;
            }  //  ::::::   end of velocity loop_2D: if (velocity_iter_2D > velocity_iter_max_2D)   ::::::::::::::::::::::
            computePressure_2D();
            if(iter_cnt > nm){
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        } // :::::::::::::::::::   end of pressure loop_2D: if (pressure_iter_2D > pressure_iter_max_2D)   ::::::::::
    } // ::::::::   end of 2D loop for initial surface conditions: if (switch_2D == 0)   :::::::::::::::::::::::::::::
}


void cAtmosphereModel::run_3D_loop(){ 

    cout << "    begin run_3D_loop    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

    iter_cnt = 1;
    iter_cnt_3d = 0;
    if(debug){
        save_data();
        //chessboard_grid(t.x[0], 30, 30, jm, km);
    }
    for(int pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++){
        for(int velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++){

    cout << "    in run_3D_loop    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

            if(debug){ Array tmp = (t-1)*t_0; tmp.inspect();}
            cout << endl << endl;
            cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            cout << " 3D AGCM iterational process" << endl;
            cout << " max total iteration number nm = " << nm << endl << endl;
            cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum "
                << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl << 
                endl << " Ma = " << *get_current_time() << "     n = " << iter_cnt << "    velocity_iter_max = " << velocity_iter_max << 
                "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << 
                "    pressure_iter = " << pressure_iter << endl;
            BC_radius();
            BC_theta();
            BC_phi();
            BC_SolidGround();
            if(velocity_iter % 2 == 0){
                Ice_Water_Saturation_Adjustment();
//            Moist_Convection(); // work in progress
                Two_Category_Ice_Scheme(); 
                Value_Limitation_Atm();
//                if(pressure_iter == 2)  WaterVapourEvaporation();
                WaterVapourEvaporation();
                init_co2();
            }
            solveRungeKutta_3D_Atmosphere();
            Value_Limitation_Atm();
            store_intermediate_data_3D();
            if(debug)check_data(); 
            BC_Radiation_multi_layer(); 
            Latent_Heat(); 
            print_min_max_values();
            vegetation_distribution();
            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        }//end of velocity loop
        computePressure_3D();
//        if(debug && pressure_iter % checkpoint == 0){
        if(pressure_iter % checkpoint == 0){
            write_file(bathymetry_name, output_path);
        }
        if(iter_cnt > nm){
            cout << "       nm = " << nm 
                << "     .....     maximum number of iterations   nm   reached!" 
                << endl;
            break;
        }
    }//end of pressure loop
}


/*
*
*/
void cAtmosphereModel::load_temperature_curve(){
    load_map_from_file(temperature_curve_file, m_temperature_curve);
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
        / (bottom->first - upper->first) 
        * (bottom->second - upper->second);
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
            weight = cos((90-i) * M_PI / 180.0);
        }else{
            weight = cos((i-90) * M_PI / 180.0);
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
                t.x[0][j][k] -= exp((t.x[0][j][k] - 1) * t_0 / 4 - 10) * 6 / t_0;
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
}

/*
*
*/
void cAtmosphereModel::init_water_vapour(){
    // initial and boundary conditions of water vapour on water and land surfaces
    // value 0.04 stands for the maximum value of 40 g/kg, g water vapour per kg dry air
    // water vapour contents computed by Clausius-Clapeyron-formula
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
//            int i_mount = get_surface_layer(j, k);
            int i_mount = 0.;
            c.x[i_mount][j][k] = hp * ep * exp (17.0809 
                * (t.x[i_mount][j][k] * t_0 - t_0) / (234.175 
                + (t.x[i_mount][j][k] * t_0 - t_0))) / ((r_air * R_Air 
                * t.x[i_mount][j][k] * t_0) * .01);
                // saturation of relative water vapour in kg/kg
            if(is_air(h, 0, j, k)){
                c.x[i_mount][j][k] = c_ocean * c.x[i_mount][j][k];
                // relativ water vapour contents on ocean surface reduced by factor
            }
            else if(is_land(h, 0, j, k)){
                c.x[i_mount][j][k] = c_land * c.x[i_mount][j][k];
            }
        }
    }
    // water vapour distribution decreasing approaching tropopause
    for(int j = 0; j < jm; j++){
//        int i_trop = get_tropopause_layer(j);
        int i_trop = im-1;
        for(int k = 0; k < km; k++){
//            int i_mount = get_surface_layer(j, k);
            int i_mount = 0.;
            for(int i = 0; i < im; i++){
                if(i <= i_trop){
                    if(i>i_mount){
                        double x = get_layer_height(i) 
                            / get_layer_height(i_trop); 
                        c.x[i][j][k] = parabola_interp(c_tropopause, 
                            c.x[i_mount][j][k], x); 
                    }
                }
            } // end i
        }// end k
    }// end j
    // cloud water and cloud ice distribution formed by the curve Versiera (Witch) of Agnesi, algebraic curve of 3. order
    for(int j = 0; j < jm; j++){
//        int i_trop = get_tropopause_layer(j);
        int i_trop = im-1;
        for(int k = 0; k < km; k++){
//            int i_mount = get_surface_layer(j, k);
            int i_mount = 0.;
            double i_max = (double)(im-1);
            double cloud_x_max = .01; // scaling of the maximum x value in the linear function x(i) 
            double ice_x_max = .01; // scaling of the maximum x value in the linear function x(i)
            int cloud_Agnesi = 20; // radial location of cloud water maximum
            int ice_Agnesi = 33; // radial location of cloud ice maximum
            double cloud_max = c.x[i_mount][j][k] / 3.; // maximum cloud water based on the surface water vapour
            double ice_max = c.x[i_mount][j][k] / 14.; // maximum ice water based on the surface water vapour
            for(int i = 0; i < im; i++){
                if(i <= i_trop){
                    if(i>i_mount){
                        double x_cloud = ( (double)i - cloud_Agnesi ) 
                            / ( i_max - cloud_Agnesi ) * cloud_x_max;
                        double x_ice = ( (double)i - ice_Agnesi ) 
                            / ( ice_Agnesi - ice_max ) * ice_x_max;
                        cloud.x[i][j][k] = pow(cloud_max, 3.) 
                            / ( pow(cloud_max, 2.) + pow(x_cloud, 2.)); 
                        ice.x[i][j][k] = pow(ice_max, 3.) 
                            / ( pow(ice_max, 2.) + pow(x_ice, 2.)); 
                    }
                }
            } // end i
        }// end k
    }// end j
    float c43 = 4./3.;
    float c13 = 1./3.;
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            cloud.x[im-1][j][k] = c43 * cloud.x[im-2][j][k] -
                c13 * cloud.x[im-3][j][k];
            ice.x[im-1][j][k] = c43 * ice.x[im-2][j][k] -
                c13 * ice.x[im-3][j][k];
            if(is_land (h, 0, j, k)){  
                cloud.x[0][j][k] = 0.;
                ice.x[0][j][k] = 0.;
            }else{
                cloud.x[0][j][k] = c43 * cloud.x[1][j][k] -
                    c13 * cloud.x[2][j][k];
                ice.x[0][j][k] = c43 * ice.x[1][j][k] -
                    c13 * ice.x[2][j][k];
            }
        }
    }
}


/*
*
*/
void cAtmosphereModel::near_wall_values(){
    // all variables need initial values to represent a flow field 
    // as close as possible to the final results due to the frictional 
    // effects of walls, the first near-wall grid point receives reduced 
    // values of the velocity variables u, v and w to form a boundary layer
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if((j >= 2) && (j < jm - 3)){
                    if((is_land(h, i, j, k)) && ((is_air(h, i, j+1, k)) 
                        && (is_air(h, i, j+2, k)))){
                            u.x[i][j+2][k] = u.x[i][j+2][k] / 2.;
                            u.x[i][j+1][k] = u.x[i][j+2][k] / 4.;
                            v.x[i][j+2][k] = v.x[i][j+2][k] / 2.;
                            v.x[i][j+1][k] = v.x[i][j+2][k] / 4.;
                            w.x[i][j+2][k] = w.x[i][j+2][k] / 2.;
                            w.x[i][j+1][k] = w.x[i][j+2][k] / 4.;
                    }
                    if((is_land(h, i, j, k)) && (is_air(h, i, j-1, k)) 
                        && (is_air(h, i, j-2, k))){
                            u.x[i][j-1][k] = u.x[i][j-2][k] / 4.;
                            u.x[i][j-2][k] = u.x[i][j-2][k] / 2.;
                            v.x[i][j-2][k] = v.x[i][j-2][k] / 2.;
                            v.x[i][j-1][k] = v.x[i][j-2][k] / 4.;
                            w.x[i][j-2][k] = w.x[i][j-2][k] / 2.;
                            w.x[i][j-1][k] = w.x[i][j-2][k] / 4.;
                    }
                    if(((is_land(h, i, j, k)) 
                            && ((is_air(h, i, j+1, k)) && (is_land(h, i, j+2, k)))) 
                            || ((j == jm - 2) && ((is_air(h, i, j, k)) 
                            && (is_land(h, i, j+1, k))))){
                                u.x[i][j+2][k] = u.x[i][j+2][k] / 2.;
                                u.x[i][j+1][k] = u.x[i][j+2][k] / 4.;
                                v.x[i][j+2][k] = v.x[i][j+2][k] / 2.;
                                v.x[i][j+1][k] = v.x[i][j+2][k] / 4.;
                                w.x[i][j+2][k] = w.x[i][j+2][k] / 2.;
                                w.x[i][j+1][k] = w.x[i][j+2][k] / 4.;
                    }
                    if(((is_land(h, i, j, k)) 
                            && ((is_air(h, i, j-1, k)) && (is_land(h, i, j-2, k)))) 
                            || ((j == 1) && ((is_land(h, i, j, k)) 
                            && (is_air(h, i, j-1, k))))){
                                u.x[i][j-2][k] = u.x[i][j-2][ k ] / 2.;
                                u.x[i][j-1][k] = u.x[i][j-2][ k ] / 4.;
                                v.x[i][j-2][k] = v.x[i][j-2][ k ] / 2.;
                                v.x[i][j-1][k] = v.x[i][j-2][ k ] / 4.;
                                w.x[i][j-2][k] = w.x[i][j-2][ k ] / 2.;
                                w.x[i][j-1][k] = w.x[i][j-2][ k ] / 4.;
                    }
                }
                if((k >= 2) && (k < km - 3)){
                    if((is_land(h, i, j, k)) && (is_air(h, i, j, k+1)) 
                        && (is_air(h, i, j, k+2))){
                            u.x[i][j][k+2] = u.x[i][j][k+2] / 2.;
                            u.x[i][j][k+1] = u.x[i][j][k+2] / 4.;
                            v.x[i][j][k+2] = v.x[i][j][k+2] / 2.;
                            v.x[i][j][k+1] = v.x[i][j][k+2] / 4.;
                            w.x[i][j][k+2] = w.x[i][j][k+2] / 2.;
                            w.x[i][j][k+1] = w.x[i][j][k+2] / 4.;
                    }
                    if((is_land(h, i, j, k)) && (is_air(h, i, j, k-1)) 
                        && (is_air(h, i, j, k-2))){
                            u.x[i][j][k-2] = u.x[i][j][k-2] / 2.;
                            u.x[i][j][k-1] = u.x[i][j][k-2] / 4.;
                            v.x[i][j][k-2] = v.x[i][j][k-2] / 2.;
                            v.x[i][j][k-1] = v.x[i][j][k-2] / 4.;
                            w.x[i][j][k-2] = w.x[i][j][k-2] / 2.;
                            w.x[i][j][k-1] = w.x[i][j][k-2] / 4.;
                    }
                    if(((is_land(h, i, j, k)) && ((is_air(h, i, j, k+1)) 
                        && (is_land(h, i, j, k+2)))) || ((k == km - 2)
                        && ((is_air(h, i, j, k)) && (is_land(h, i, j, k+1))))){
                            u.x[i][j][k+2] = u.x[i][j][k+2] / 2.;
                            u.x[i][j][k+1] = u.x[i][j][k+2] / 4.;
                            v.x[i][j][k+2] = v.x[i][j][k+2] / 2.;
                            v.x[i][j][k+1] = v.x[i][j][k+2] / 4.;
                            w.x[i][j][k+2] = w.x[i][j][k+2] / 2.;
                            w.x[i][j][k+1] = w.x[i][j][k+2] / 4.;
                    }
                    if(((is_land(h, i, j, k)) 
                            && ((is_air(h, i, j, k-1)) && (is_land(h, i, j, k-2)))) 
                            || ((k == 1) && ((is_land(h, i, j, k)) 
                            && (is_air(h, i, j, k-1))))){
                                u.x[i][j][k-2] = u.x[i][j][k-2] / 2.;
                                u.x[i][j][k-1] = u.x[i][j][k-2] / 4.;
                                v.x[i][j][k-2] = v.x[i][j][k-2] / 2.;
                                v.x[i][j][k-1] = v.x[i][j][k-2] / 4.;
                                w.x[i][j][k-2] = w.x[i][j][k-2] / 2.;
                                w.x[i][j][k-1] = w.x[i][j][k-2] / 4.;
                    }
                }
            } // end i
        }// end k
    }// end j
}

/*
*
*/
void cAtmosphereModel::init_topography(const string &topo_filename){
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
}

/*
*
*/
void cAtmosphereModel::init_co2(){
    // initial and boundary conditions of CO2 content on water and land surfaces
    // parabolic CO2 content distribution from pole to pole accepted
    // CO2-distribution by Ruddiman approximated by a parabola
    co2_paleo = 3.2886 * pow ((t_paleo + t_average), 2) - 32.8859 *
        (t_paleo + t_average) + 102.2148;  // in ppm
    co2_average = 3.2886 * pow (t_average, 2) - 32.8859 * t_average + 102.2148;  // in ppm
    co2_paleo = co2_paleo - co2_average;
    cout.precision(3);
    const char* co_comment = "      co2 increase at paleo times: ";
    const char* co_gain = " co2 increase";
    const char* co_modern = "      mean co2 at modern times: ";
    const char* co_paleo_str = "      mean co2 at paleo times: ";
    const char* co_average_str = " co2 modern";
    const char* co_average_pal = " co2 paleo";
    const char* co_unit =  "ppm ";
    cout << endl << setiosflags (ios::left) << setw (55) << setfill ('.') <<
        co_comment << resetiosflags (ios::left)         << setw (12) << co_gain << " = "
        << setw (7) << setfill (' ') << co2_paleo << setw (5) << co_unit << 
        endl << setw (55) << setfill ('.')  << setiosflags (ios::left) << co_modern
        << resetiosflags (ios::left) << setw (13) << co_average_str  << " = "
        << setw (7)  << setfill (' ') << co2_average << setw (5) << co_unit 
        << endl << setw (55) << setfill ('.')  << setiosflags (ios::left)
        << co_paleo_str << resetiosflags (ios::left) << setw (13) << co_average_pal
        << " = "  << setw (7)  << setfill (' ') << co2_average + co2_paleo
        << setw (5) << co_unit << endl;
    cout << endl;
    co2_paleo = co2_paleo / co2_0;
    co2_land = co2_land / co2_0;
    co2_ocean = co2_ocean / co2_0;
    co2_tropopause = co2_tropopause / co2_0;
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
                     * pow((t.x[i_mount][j][k] * t_0), 3) / emittancy_total)
                     + co2_paleo + co2_ocean;
            }
            if(is_land(h, i_mount, j, k)){
                     // reciprocal formula for the temperature increase by co2, 
                     // delta_T = emittancy_total*ln(co2)/(4*coeff_em * t³)
                     // original formula by Nasif Nahle Sabag delta T(co2)
                     // co2_paleo taken over from Ruddiman, p 86, effect of co2 on global temperature
                co2.x[i_mount][j][k] = exp (4. * delta_T * coeff_em 
                     * pow((t.x[i_mount][j][k] * t_0), 3) / emittancy_total)
                     + co2_paleo + co2_land - co2_vegetation
                     * Vegetation.y[j][k] / co2_0;
            }
        }
    }
    if(iter_cnt_3d > -1){
        // co2 distribution decreasing approaching tropopause, above no co2
        for(int j = 0; j < jm; j++){
            int i_trop = get_tropopause_layer(j);
            for(int k = 0; k < km; k++){
//             int i_mount = i_topography[j][k];
                int i_mount = 0;
                for(int i = 0; i < im; i++){
                    if(i <= i_trop){
                        // radial distribution approximated by a parabola
                        double x = get_layer_height(i) 
                            / get_layer_height(i_trop); 
                        co2.x[i][j][k] = parabola_interp(co2_tropopause, 
                            co2.x[i_mount][j][k], x); 
                    }else  co2.x[i][j][k] = co2_tropopause;
                }
            }
        }
    }
    float c43 = 4./3.;
    float c13 = 1./3.;
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            co2.x[im-1][j][k] = c43 * co2.x[im-2][j][k] -
                c13 * co2.x[im-3][j][k];
            if(is_land (h, 0, j, k)){  
                co2.x[0][j][k] = c43 * co2.x[1][j][k] -
                    c13 * co2.x[2][j][k];
            }
        }
    }
}
/*
*
*/
void  cAtmosphereModel::save_data(){
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
}

/*
*
*/
void cAtmosphereModel::BC_SolidGround(){
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = im-2; i >= 0; i--){
                if(is_land(h, i, j, k)){
                    u.x[i][j][k] = 0.;
                    v.x[i][j][k] = 0.;
                    w.x[i][j][k] = 0.;
                    //t.x[i][j][k] = 1.;  // = 273.15 K
                    //c.x[i][j][k] = c_tropopause;  // = 1 g/kg water vapour
                    //c.x[i][j][k] = 0.; 
                    cloud.x[i][j][k] = 0.;
                    ice.x[i][j][k] = 0.;
//                    co2.x[i][j][k] = 1.;  // = 280 ppm
                    p_dyn.x[i][j][k] = 0.;
                }// is_land
            } // i
        } // k
    } // j
}

/*
*
*/
void cAtmosphereModel::vegetation_distribution(){
    // vegetation as function of precipitation, timberline and temperature
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            if(max_Precipitation > 0 && is_land(h, 0, j, k) 
                && !( get_layer_height(i_mount) > 2400.) // above the timberline no vegeation
                && !( ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) < - 10. ) ){ //  vegetation >= -10°C
                Vegetation.y[j][k] = Precipitation.y[j][k] 
                    / max_Precipitation; // actual vegetation areas
            }else{
                Vegetation.y[j][k] = 0.;
            }
        }
    }
}

void cAtmosphereModel::store_intermediate_data_2D(float coeff){
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            vn.x[0][j][k] = coeff * v.x[0][j][k];
            wn.x[0][j][k] = coeff * w.x[0][j][k];
            p_dynn.x[0][j][k] = coeff * p_dyn.x[0][j][k];
        }
    }
}

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
}

void cAtmosphereModel::adjust_temperature_IC(double** t, int jm, int km){
    for(int k=0; k < km; k++){
        for(int j=0; j < jm; j++){
            t[j][k] = temperature_NASA.y[j][k] = (t[j][k] + t_0) / t_0;
        }
    }
    // correction of surface temperature around 180°E
    int k_half = (km-1)/2;
    for(int j = 0; j < jm; j++){
        t[j][k_half] = (t[j][k_half+1] + t[j][k_half-1]) / 2.;
        temperature_NASA.y[j][k_half] = (temperature_NASA.y[j][k_half+1] +
            temperature_NASA.y[j][k_half-1]) / 2.;
    }
}

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
        logger() << "layer: " << i << std::endl;
        logger() << name << " min: " << t_min << std::endl;
        logger() << name << " max: " << t_max << std::endl;
        logger() << name << " diff min: " << t_diff_min << std::endl;
        logger() << name << " diff max: " << t_diff_max << std::endl;
        logger() << name << " diff mean: " << t_diff_mean << "   " 
            << (double)t_diff_mean / (jm*km) << std::endl;
    }
}


void cAtmosphereModel::check_data(){
    check_data(t,tn,"t");
    check_data(u,un,"u");
    check_data(v,vn,"v");
    check_data(w,wn,"w");
    check_data(c,cn,"c");
    check_data(cloud,cloudn,"cloud");
    check_data(ice,icen,"ice");
    check_data(co2,co2n,"c02");
}


void cAtmosphereModel::init_tropopause_layers(){
    tropopause_layers = std::vector<int>(jm, tropopause_pole);
// Versiera di Agnesi approach, two inflection points
    int j_max = jm-1;
    int j_half = j_max/2;
    int coeff_trop = 8;
    // scaling of the maximum x value in the linear function x(j) to get tropopause_pole
    for(int j=j_half; j<=jm-1; j++){
        int x = coeff_trop * (double)j / (double)j_half;
        tropopause_layers[j] = Agnesi(x, tropopause_equator); 
    }
    for(int j=0; j<j_half; j++){
        tropopause_layers[j] = tropopause_layers[j_max-j];
    }
    if(debug){
        for(int j=0; j<jm; j++){
            std::cout << tropopause_layers[j] << " ";
        }
        std::cout << std::endl;
    }
}

