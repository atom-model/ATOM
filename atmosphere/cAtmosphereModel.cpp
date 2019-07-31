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
//const double cAtmosphereModel::dt = 0.0001;
    
const double cAtmosphereModel::the0 = 0.;             // North Pole
const double cAtmosphereModel::phi0 = 0.;             // zero meridian in Greenwich

//earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km 
const double cAtmosphereModel::r0 = 1.; 

cAtmosphereModel::cAtmosphereModel() :
    i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    is_node_weights_initialised(false), 
    has_welcome_msg_printed(false)
{
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.
    //PythonStream::OverrideCout();
    if(PythonStream::is_enable())
    {
        backup = std::cout.rdbuf();
        std::cout.rdbuf(&ps);
    }
    
    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

    // set default configuration
    SetDefaultConfig();

    coeff_mmWS = r_air / r_water_vapour; // coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]

    m_model = this;

    //  Coordinate system in form of a spherical shell
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.initArray_1D(im, 0); // radial coordinate direction
    the.initArray_1D(jm, 0); // lateral coordinate direction
    phi.initArray_1D(km, 0); // longitudinal coordinate direction
    rad.Coordinates ( im, r0, dr );
    the.Coordinates ( jm, the0, dthe );
    phi.Coordinates ( km, phi0, dphi );

    init_layer_heights();
}

cAtmosphereModel::~cAtmosphereModel() {
    if(PythonStream::is_enable()){
        std::cout.rdbuf(backup);
    }

    m_model = NULL;
    logger().close();
}
 
#include "cAtmosphereDefaults.cpp.inc"

void cAtmosphereModel::LoadConfig ( const char *filename ) 
{
    XMLDocument doc;
    XMLError err = doc.LoadFile ( filename );
    try{
        if (err) {
            doc.PrintError();
            throw std::invalid_argument(std::string("unable to load config file:  ") + filename);
        }

        XMLElement *atom = doc.FirstChildElement("atom"), *elem_common = NULL, *elem_atmosphere = NULL;
        if (!atom) {
            throw std::invalid_argument(std::string("Failed to find the 'atom' element in config file: ") + filename);
        }else{
            elem_common = atom->FirstChildElement( "common" );
            if(!elem_common){
                throw std::invalid_argument(std::string("Failed to find the 'common' element in 'atom' element in config file: ") + filename);
            }
            elem_atmosphere = atom->FirstChildElement( "atmosphere" );
            if (!elem_atmosphere) {
                throw std::invalid_argument(std::string("Failed to find the 'atmosphere' element in 'atom' element in config file: ") + filename);
            }
        }

        #include "AtmosphereLoadConfig.cpp.inc"

    }catch(const std::exception &exc){
        std::cerr << exc.what() << std::endl;
        abort();
    }
}



void cAtmosphereModel::RunTimeSlice ( int Ma )
{
    if(debug){
        feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO); //not platform independent, bad, very bad, I know
    }

    if(!is_temperature_curve_loaded()) 
        load_temperature_curve();

    reset_arrays();    

    m_current_time = m_time_list.insert(float(Ma)).first;

    struct stat info;
    if( stat( output_path.c_str(), &info ) != 0 ){
        mkdir(output_path.c_str(), 0777);
    }


    //Prepare the temperature and precipitation data file
    string Name_SurfaceTemperature_File  = temperature_file;
    string Name_SurfacePrecipitation_File = precipitation_file;

    if(Ma != 0 && use_earthbyte_reconstruction){
        Name_SurfaceTemperature_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";
        Name_SurfacePrecipitation_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";    
        velocity_v_file = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_wind_v.xyz";
        velocity_w_file = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_wind_w.xyz";
    
        if( stat( Name_SurfaceTemperature_File.c_str(), &info ) != 0 || 
            stat( Name_SurfacePrecipitation_File.c_str(), &info ) != 0 ||
            stat( velocity_v_file.c_str(), &info ) != 0 ||
            stat( velocity_w_file.c_str(), &info ) != 0 )
        {
            std::string cmd_str = "python " + reconstruction_script_path + " " + std::to_string(Ma - time_step) + " " + 
                std::to_string(Ma) + " " + output_path + " " + BathymetrySuffix + " atm";
            int ret = system(cmd_str.c_str());
            std::cout << " reconstruction script returned: " << ret << std::endl;
        } 
    }

    if(!has_welcome_msg_printed)
        print_welcome_msg();

    //  initialization of the bathymetry/topography
    //  topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    init_topography(bathymetry_path + "/" + bathymetry_name);

    read_IC(velocity_v_file, v.x[0], jm, km);
    read_IC(velocity_w_file, w.x[0], jm, km);    
    read_IC(Name_SurfaceTemperature_File, t.x[0], jm, km);
    read_IC(Name_SurfacePrecipitation_File, Precipitation.y, jm, km);

//    goto Printout;

    iter_cnt_3d = -1;
    if(debug) save_data();
    iter_cnt_3d++;

    //  class element for the initial conditions for u-v-w-velocity components
    //circulation.IC_CellStructure ( h, u, v, w );
    init_velocities();

//    goto Printout;

    //IC_v_w_WestEastCoast();//adjust east coast velocities.

    adjust_temperature_IC(t.x[0], jm, km);
    init_temperature();

//    goto Printout;

    //  class element for the surface pressure computed by surface temperature with gas equation
    BC_Pressure();

    //parabolic water vapour distribution from pol to pol, maximum water vapour volume at equator
    //circulation.BC_WaterVapour ( h, p_stat, t, c, v, w );
    init_water_vapour();

//    goto Printout;

    //  class element for the parabolic CO2 distribution from pol to pol, maximum CO2 volume at equator
    //circulation.BC_CO2 ( Vegetation, h, t, p_dyn, co2 );
    init_co2();

    // class element for the surface temperature computation by radiation flux density
    if ( RadiationModel == 1 ){
        BC_Radiation_multi_layer(); 
    }

//    goto Printout;

    // class element for the storing of velocity components, pressure and temperature for iteration start
    store_intermediate_data_2D();
    store_intermediate_data_3D();

    run_2D_loop();
    
    cout << endl << endl;

    run_3D_loop();

    cout << endl << endl;

    restrain_temperature();


//    Printout:

    //write the ouput files
    write_file(bathymetry_name, output_path, true);

    iter_cnt_3d++;
    save_data();    

    //  final remarks
    cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;

    if(debug){
        fedisableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO); //not platform independent(bad, very bad, I know)
    }
}



void cAtmosphereModel::Run() 
{
    auto start_time = std::chrono::system_clock::now();
    std::time_t start_time_t = std::chrono::system_clock::to_time_t(start_time);
    logger() << "Start Time:" << std::ctime(&start_time_t) << std::endl;

    mkdir(output_path.c_str(), 0777);

    print_welcome_msg();

    for(int i = time_start; i <= time_end; i+=time_step)
    {
        RunTimeSlice(i);
    }

    print_final_remarks();

    auto end_time = std::chrono::system_clock::now();
    std::time_t end_time_t = std::chrono::system_clock::to_time_t(end_time);
    logger() << "End Time:" << std::ctime(&end_time_t) << std::endl;
}


void cAtmosphereModel::reset_arrays()
{
    // 2D arrays
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
    Evaporation_Penman.initArray_2D(jm, km, 0.); // evaporation by Penman in [mm/d]

    co2_total.initArray_2D(jm, km, 0.); // areas of higher co2 concentration

    // 3D arrays
    h.initArray(im, jm, km, 0.); // bathymetry, depth from sea level

    t.initArray(im, jm, km, ta); // temperature
    u.initArray(im, jm, km, ua); // u-component velocity component in r-direction
    v.initArray(im, jm, km, va); // v-component velocity component in theta-direction
    w.initArray(im, jm, km, wa); // w-component velocity component in phi-direction
    c.initArray(im, jm, km, ca); // water vapour
    cloud.initArray(im, jm, km, 0.); // cloud water
    ice.initArray(im, jm, km, 0.); // cloud ice
    co2.initArray(im, jm, km, coa); // CO2

    tn.initArray(im, jm, km, ta); // temperature new
    un.initArray(im, jm, km, ua); // u-velocity component in r-direction new
    vn.initArray(im, jm, km, va); // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, wa); // w-velocity component in phi-direction new
    cn.initArray(im, jm, km, ca); // water vapour new
    cloudn.initArray(im, jm, km, 0.); // cloud water new
    icen.initArray(im, jm, km, 0.); // cloud ice new
    co2n.initArray(im, jm, km, coa); // CO2 new

    p_dyn.initArray(im, jm, km, pa); // dynamic pressure
    p_dynn.initArray(im, jm, km, pa); // dynamic pressure
    p_stat.initArray(im, jm, km, pa); // static pressure

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
    S_v.initArray(im, jm, km, 0.); // water vapour mass rate due to category two ice scheme
    S_c.initArray(im, jm, km, 0.); // cloud water mass rate due to category two ice scheme
    S_i.initArray(im, jm, km, 0.); // cloud ice mass rate due to category two ice scheme
    S_r.initArray(im, jm, km, 0.); // rain mass rate due to category two ice scheme
    S_s.initArray(im, jm, km, 0.); // snow mass rate due to category two ice scheme
    S_c_c.initArray(im, jm, km, 0.); // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique

    for (auto &i : i_topography)
        std::fill(i.begin(), i.end(), 0);
}


void cAtmosphereModel::write_file(std::string &bathymetry_name, std::string &output_path, bool is_final_result){
    int Ma = int(round(*get_current_time()));
    //  printout in ParaView files and sequel files
    //  writing of data in ParaView files
    //  radial data along constant hight above ground
    int i_radial = 0;
    //  int i_radial = 10;
    paraview_vtk_radial ( bathymetry_name, Ma, i_radial, iter_cnt-1 ); 

    //  londitudinal data along constant latitudes
    int j_longal = 62;          // Mount Everest/Himalaya
    paraview_vtk_longal ( bathymetry_name, j_longal, iter_cnt-1 ); 

    int k_zonal = 87;           // Mount Everest/Himalaya
    paraview_vtk_zonal ( bathymetry_name, k_zonal, iter_cnt-1 ); 

    //  3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
    if(paraview_panorama_vts_flag) //This function creates a large file. Use a flag to control if it is wanted.
    {
        paraview_panorama_vts ( bathymetry_name, iter_cnt-1 ); 
    }

    Value_Limitation_Atm();
    
    //  writing of v-w-data in the v_w_transfer file
    Atmosphere_v_w_Transfer ( bathymetry_name );
    Atmosphere_PlotData ( bathymetry_name, (is_final_result ? -1 : iter_cnt-1) );
}

void cAtmosphereModel::run_2D_loop(){
    int switch_2D = 0;    
    iter_cnt = 1;
    int Ma = int(round(*get_current_time())); 

    // ::::::::::: :::::::::::::::::::::::   begin of 2D loop for initial surface conditions: if ( switch_2D == 0 )   ::::
    if ( switch_2D != 1 )
    {
        // **************   iteration of initial conditions on the surface for the correction of flows close to coasts   **
        // **************   start of pressure and velocity iterations for the 2D iterational process   ********************
        // ::::::::::::::   begin of pressure loop_2D : if ( pressure_iter_2D > pressure_iter_max_2D )   ::::::::::::::::::
        for ( int pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++)
        {
            // ::::::::   begin of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::::::
            for ( int velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++)
            {

                cout << endl << endl;
                cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
                cout << " 2D AGCM iterational process" << endl;
                cout << " max total iteration number nm = " << nm << endl << endl;

                cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, maximum "
                    << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl 
                    << endl << " Ma = " << Ma << "     n = " << iter_cnt << "    velocity_iter_max_2D = " << velocity_iter_max_2D
                    << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " << pressure_iter_max_2D << 
                    "    pressure_iter_2D = " << pressure_iter_2D << endl;

                //  class BC_Atmosphaere for the geometry of a shell of a sphere
                BC_theta();
                BC_phi();
                
                Value_Limitation_Atm( );

                BC_SolidGround(); 
                
                //  class RungeKutta for the solution of the differential equations describing the flow properties
                solveRungeKutta_2D_Atmosphere();
                
                store_intermediate_data_2D();

                iter_cnt++;
            }
            //  ::::::   end of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::::::::::::::


            //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
            computePressure_2D();

            // limit of the computation in the sense of time steps
            if ( iter_cnt > nm )
            {
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        }
        // :::::::::::::::::::   end of pressure loop_2D: if ( pressure_iter_2D > pressure_iter_max_2D )   ::::::::::
    }
    // ::::::::   end of 2D loop for initial surface conditions: if ( switch_2D == 0 )   :::::::::::::::::::::::::::::
}


void cAtmosphereModel::run_3D_loop(){ 
    iter_cnt = 1;
    iter_cnt_3d = 0;

    if(debug)
    {   
        save_data();
        //chessboard_grid(t.x[0], 30, 30, jm, km);
    }

    store_intermediate_data_3D();
    
    for ( int pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++ )
    {
        for ( int velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++ )
        {
            if(debug){ Array tmp = (t-1)*t_0; tmp.inspect(); }
            
            cout << endl << endl;
            cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            cout << " 3D AGCM iterational process" << endl;
            cout << " max total iteration number nm = " << nm << endl << endl;

            cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum "
                << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl << 
                endl << " Ma = " << *get_current_time() << "     n = " << iter_cnt << "    velocity_iter_max = " << velocity_iter_max << 
                "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << 
                "    pressure_iter = " << pressure_iter << endl;

            //  class BC_Atmosphaere for the geometry of a shell of a sphere
            BC_radius();
            BC_theta();
            BC_phi();

            //Ice_Water_Saturation_Adjustment, distribution of cloud ice and cloud water dependent on water vapour amount and temperature
            if ( velocity_iter % 2 == 0 ){
                Ice_Water_Saturation_Adjustment();
            }

            Value_Limitation_Atm();

            BC_SolidGround();

            BC_Evaporation(); 
            
            // class RungeKutta for the solution of the differential equations describing the flow properties
            solveRungeKutta_3D_Atmosphere();

            if(debug)check_data(); 
            
            Value_Limitation_Atm();

            // class element for the surface temperature computation by radiation flux density
            if ( RadiationModel == 1 ){
                BC_Radiation_multi_layer(); 
            }
            
            //  class element for the initial conditions the latent heat
            Latent_Heat(); 

            print_min_max_values();

            //  computation of vegetation areas
            vegetationDistribution();

            //  composition of results
            run_MSL_data(); 

            //  Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
            //  resulting the precipitation distribution formed of rain and snow
            if ( velocity_iter % 2 == 0){
                Two_Category_Ice_Scheme(); 
//                Moist_Convection(); 
            }

            store_intermediate_data_3D();
            
            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        }//end of velocity loop
        
        //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
        computePressure_3D();
        
        if( debug && pressure_iter % checkpoint == 0 ){
            write_file(bathymetry_name, output_path);
        }

        //  limit of the computation in the sense of time steps
        if ( iter_cnt > nm )
        {
            cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
            break;
        }
    }//end of pressure loop
}


/*
*
*/
void cAtmosphereModel::load_temperature_curve()
{
    load_map_from_file(temperature_curve_file, m_temperature_curve);
}

/*
*
*/
float cAtmosphereModel::get_mean_temperature_from_curve(float time) const
{
    if(time<m_temperature_curve.begin()->first || time>(--m_temperature_curve.end())->first){
        std::cout << "Input time out of range: " <<time<< std::endl;    
        return NAN;
    }
    if(m_temperature_curve.size()<2){
        std::cout << "No enough data in m_temperature_curve  map" << std::endl;
        return NAN;
    }
    map<float, float >::const_iterator upper=m_temperature_curve.begin(), bottom=++m_temperature_curve.begin(); 
    for(map<float, float >::const_iterator it = m_temperature_curve.begin();
            it != m_temperature_curve.end(); ++it)
    {
        if(time < it->first){
            bottom = it;
            break;
        }else{
            upper = it;
        }
    }
    //std::cout << upper->first << " " << bottom->first << std::endl;
    return upper->second + (time - upper->first) / (bottom->first - upper->first) * (bottom->second - upper->second);
}

/*
*
*/
float cAtmosphereModel::calculate_mean_temperature(const Array& temp) 
{
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
void cAtmosphereModel::calculate_node_weights()
{
    //use cosine of latitude as weights for now
    //longitudes: 0-360(km) latitudes: 90-(-90)(jm)
    double weight = 0.;
    m_node_weights.clear();
    for(int i=0; i<jm; i++){
        if(i<=90){
            weight = cos((90-i) * M_PI / 180.0 );
        }else{
            weight = cos((i-90) * M_PI / 180.0 );
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
    // parabolic water vapour distribution from pole to pole accepted

    // maximum water vapour content on water surface at equator c_equator = 1.04 compares to 0.04 volume parts
    // minimum water vapour at tropopause c_tropopause = 0.0 compares to 0.0 volume parts
    // value 0.04 stands for the maximum value of 40 g/kg, g water vapour per kg dry air

    // water vapour contents computed by Clausius-Clapeyron-formula
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            int i_mount = get_surface_layer(j, k);
            c.x[ i_mount ][ j ][ k ] = hp * ep * exp ( 17.0809 * ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 +
                    ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) ) ) / ( ( r_air * R_Air * t.x[ i_mount ][ j ][ k ] * t_0 ) * .01 );
                // saturation of relative water vapour in kg/kg

            if ( is_air ( h, 0, j, k ) ){
                c.x[ i_mount ][ j ][ k ] = c_ocean * c.x[ i_mount ][ j ][ k ];
                // relativ water vapour contents on ocean surface reduced by factor
            }
            else if ( is_land ( h, 0, j, k ) ){
                c.x[ i_mount ][ j ][ k ] = c_land * c.x[ i_mount ][ j ][ k ];
            }
        }
    }

    // water vapour distribution decreasing approaching tropopause
    for ( int j = 0; j < jm; j++ ){
        int i_trop = get_tropopause_layer(j);
        for ( int k = 0; k < km; k++ ){
            int i_mount = get_surface_layer(j, k);

            for ( int i = 0; i < im; i++ ){
                if ( i < i_trop ){
                    if(i>i_mount){
                        double x = (get_layer_height(i) - get_layer_height(i_mount)) / 
                            (get_layer_height(i_trop) - get_layer_height(i_mount));
                        c.x[ i ][ j ][ k ] = parabola_interp(c_tropopause, c.x[ i_mount ][ j ][ k ], x); 
                    }else{
                        c.x[ i ][ j ][ k ] = c.x[ i_mount ][ j ][ k ];
                    }
                }else{
                    c.x[ i ][ j ][ k ] = c_tropopause;
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
    if ( ! ifile.is_open()) {
        std::cerr << "ERROR: could not open Name_Bathymetry_File file: " <<  topo_filename << std::endl;
        abort();
    }

    double lon, lat, height;
    int j, k;
    for (j = 0; j < jm && !ifile.eof(); j++) {
        for (k = 0; k < km && !ifile.eof(); k++) {
            height = -999; // in case the height is NaN
            ifile >> lon >> lat >> height;
            if ( ! (height > 0) ){
                h.x[ 0 ][ j ][ k ] = Topography.y[ j ][ k ] = 0;
            }else{
                Topography.y[ j ][ k ] = height;
                for ( int i = 0; i < im; i++ ){
                    if(height > get_layer_height(i)){
                        h.x[ i ][ j ][ k ] = 1;
                    }else{
                        i_topography[ j ][ k ] = i-1;
                        break;
                    }   
                }   
            }   
            if(ifile.fail()){
                ifile.clear();
                std::string tmp;
                std::getline(ifile, tmp);
                logger() << "bad data in topography at: " << lon << " " << lat << " " << tmp << std::endl;
            }   
            //logger() << lon << " " << lat << " " << h.x[ 0 ][ j ][ k ] << std::endl;            
        }   
    }   
    if(j != jm || k != km ){
       std::cerr << "wrong topography file size! aborting..."<<std::endl;
        abort();
    }

    // rewriting bathymetrical data from -180° _ 0° _ +180° coordinate system to 0°- 360°
    for ( int j = 0; j < jm; j++ ){
        move_data(Topography.y[ j ], km);
        move_data(i_topography[ j ], km);
        for ( int i = 0; i < im; i++ ){
            move_data(h.x[ i ][ j ], km);
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
    co2_paleo = 3.2886 * pow ( ( t_paleo + t_average ), 2 ) - 32.8859 *
        ( t_paleo + t_average ) + 102.2148;  // in ppm
    co2_average = 3.2886 * pow ( t_average, 2 ) - 32.8859 * t_average + 102.2148;  // in ppm
    co2_paleo = co2_paleo - co2_average;

    cout.precision ( 3 );

    const char* co_comment = "      co2 increase at paleo times: ";
    const char* co_gain = " co2 increase";
    const char* co_modern = "      mean co2 at modern times: ";
    const char* co_paleo_str = "      mean co2 at paleo times: ";
    const char* co_average_str = " co2 modern";
    const char* co_average_pal = " co2 paleo";
    const char* co_unit =  "ppm ";

    cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) <<
        co_comment << resetiosflags ( ios::left )         << setw ( 12 ) << co_gain << " = "
        << setw ( 7 ) << setfill ( ' ' ) << co2_paleo << setw ( 5 ) << co_unit << 
        endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) << co_modern
        << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_str  << " = "
        << setw ( 7 )  << setfill ( ' ' ) << co2_average << setw ( 5 ) << co_unit 
        << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left )
        << co_paleo_str << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_pal
        << " = "  << setw ( 7 )  << setfill ( ' ' ) << co2_average + co2_paleo
        << setw ( 5 ) << co_unit << endl;
    cout << endl;

    co2_equator = co2_equator / co2_0;
    co2_pole = co2_pole / co2_0;
    co2_paleo = co2_paleo / co2_0;
    co2_land = co2_land / co2_0;
    co2_ocean = co2_ocean / co2_0;
    co2_tropopause = co2_tropopause / co2_0;
    double emittancy_total = 0.423; // in W/m²
    double coeff_em = 5.6697e-8; // in W/(m² K)
    double delta_T = 0.02; // in K
    // CO2-content as initial solution
    for( int k = 0; k < km; k++ ){
        for( int j = 0; j < jm; j++ ){
//            i_mount = i_topography[ j ][ k ];
            int i_mount = 0;
            if( is_air ( h, i_mount, j, k ) ){
                co2.x[ i_mount ][ j ][ k ] = exp ( 4. * delta_T * coeff_em 
                     * pow( ( t.x[ i_mount ][ j ][ k ] * t_0 ), 3 ) / emittancy_total )
                     + co2_paleo + co2_ocean;
                     // reciprocal formula for the temperature increase by co2, 
                     // Temp_co2_add by Nasif Nahle Sabag in PostProcess_Atm.cpp
                     // taken over from Ruddiman, p 86, effect of co2 on global temperature
//                co2.x[ i_mount ][ j ][ k ] = ( 3.2886 * pow ( ( t.x[ i_mount ][ j ][ k ] 
//                    * t_0 - t_0 ), 2 ) - 32.8859 * ( t.x[ i_mount ][ j ][ k ] 
//                    * t_0 - t_0 ) + 102.2148 + co2_ocean ) / co2_0;  // non-dimensional
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_average = " << t_average << "   t = " << t.x[ i_mount ][ 0 ][ k ] * t_0 - t_0 << "   co2_eff = " << ( co2_pole - co2_equator ) * co2_0 << "   co2_pole = " << co2_pole * co2_0 << "   co2_equator = " << co2_equator * co2_0 << "   co2_paleo = " << co2_paleo * co2_0 << "   co2_average = " << co2_average << "   co2 = " << co2.x[ i_mount ][ 0 ][ k ] * co2_0 << endl;
*/
            }
            if( is_land ( h, i_mount, j, k ) ){
                co2.x[ i_mount ][ j ][ k ] = exp ( 4. * delta_T * coeff_em 
                     * pow( ( t.x[ i_mount ][ j ][ k ] * t_0 ), 3 ) / emittancy_total )
                     + co2_paleo + co2_land + co2_vegetation
                    * Vegetation.y[ j ][ k ] / co2_0;
                     // reciprocal formula for the temperature increase by co2, 
                     // Temp_co2_add by Nasif Nahle Sabag in PostProcess_Atm.cpp
                     // taken over from Ruddiman, p 86, effect of co2 on global temperature
//                co2.x[ i_mount ][ j ][ k ] = ( 3.2886 * pow ( ( t.x[ i_mount ][ j ][ k ] 
//                    * t_0 - t_0 ), 2 ) - 32.8859 * ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) 
//                    + 102.2148 + co2_land - co2_vegetation * Vegetation.y[ j ][ k ] ) / co2_0;  // non-dimensional
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    if ( ( j == 90 ) && ( k == 30 ) ) cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_average = " << t_average << "   t = " << t.x[ i_mount ][ 0 ][ k ] * t_0 - t_0 << "   co2_eff = " << ( co2_pole - co2_equator ) * co2_0 << "   co2_pole = " << co2_pole * co2_0 << "   co2_equator = " << co2_equator * co2_0 << "   co2_paleo = " << co2_paleo * co2_0 << "   co2_average = " << co2_average << "   co2 = " << co2.x[ i_mount ][ 0 ][ k ] * co2_0 << "   co2_vegetation = " << co2_vegetation << "   Vegetation = " << Vegetation.y[ j ][ k ] << endl;
*/
            }
        }
    }
    // co2 distribution decreasing approaching tropopause, above no co2
    for( int j = 0; j < jm; j++ ){
        int i_trop = get_tropopause_layer(j);
        for( int k = 0; k < km; k++ ){
//             double i_mount = i_trop;
            int i_mount = 0;
            for( int i = 0; i < im; i++ ){
                if( i <= i_trop ){
                    co2.x[ i ][ j ][ k ] = co2.x[ i_mount ][ j ][ k ] 
                        - ( co2_tropopause - co2.x[ i_mount ][ j ][ k ] ) 
                        * ( get_layer_height(i) / get_layer_height(i_trop) 
                        * ( get_layer_height(i) / get_layer_height(i_trop) - 2. ) );
                        // radial distribution approximated by a parabola
                }
                else  co2.x[ i ][ j ][ k ] = co2_tropopause;
            }
        }
    }
}


/*
*
*/
void cAtmosphereModel::init_velocities(){
    // boundary condition for the velocity components in the circulation cells

    // latest version by Grotjahn ( Global Atmospheric Circulations, 1993 )
    // default for the velocity components u, v, and w as initial conditions

    // velocities given in m/s, 1 m/s compares to 3.6 km/h, non-dimensionalized by u_0 at the end of this class element
    // do not change the velocity initial conditions !!

    // initial velocity components in the northern and southern
    // Pole, Ferrel and Hadley cells

    // equator ( at j=90 compares to 0° latitude )
    // u-component up to tropopause and back on half distance (  i = 20 )
    double ua_00 = 1.;  // in m/s compares to 3.6 km/h, non-dimensionalized by u_0 at the end of this function
    double va_equator_SL =  0.000;
    double va_equator_Tropopause = 0.000;
    double wa_equator_SL = - 1.;
    double wa_equator_Tropopause = - 7.5;

    //equator
    init_u(u,90,90,ua_00);//lat: 0
    
    init_v_or_w(v,90,90,va_equator_Tropopause,va_equator_SL);
    init_v_or_w(w,90,90,wa_equator_Tropopause,wa_equator_SL);
    
    //polar cell
    double ua_90 = 0.5;
    double va_Polar_SL = 0.;
    double va_Polar_Tropopause = 0.;
    double va_Polar_SL_75 = .5;
    double va_Polar_Tropopause_75 = - 1.;
    double wa_Polar_SL = - 0.01;
    double wa_Polar_Tropopause = 0.;

    //northern polar cell
    init_u(u, 0, 0, ua_90); //lat: 90

    init_v_or_w(v,0,30,va_Polar_Tropopause,va_Polar_SL); //lat: 90-60
    init_v_or_w(w,0,30,wa_Polar_Tropopause,wa_Polar_SL); //lat: 90-60
    init_v_or_w(v,15,15,va_Polar_Tropopause_75,va_Polar_SL_75); //lat: 75

    //southern polar cell
    init_u(u, 180, 180, ua_90); //lat: -90

    init_v_or_w(v,150,180,va_Polar_Tropopause,va_Polar_SL); //lat: -90-(-60)
    init_v_or_w(w,150,180,wa_Polar_Tropopause,wa_Polar_SL); //lat: -90-(-60)
    init_v_or_w(v,165,165,va_Polar_Tropopause_75,va_Polar_SL_75); //lat: -75

    //Ferrel cell
    double ua_60 = 0.5;
    double va_Ferrel_SL = 0.5;
    double va_Ferrel_Tropopause = 1.;
    double va_Ferrel_SL_45 = - 0.1;
    double va_Ferrel_Tropopause_45 = 1.;
    double wa_Ferrel_SL = -0.2;    // subpolar jet
    double wa_Ferrel_Tropopause = 10.;           

    //northern Ferrel cell
    init_u(u, 30, 30, ua_60); //lat: 60 

    init_v_or_w(v,30,30,va_Ferrel_Tropopause,va_Ferrel_SL); //lat: 60
    init_v_or_w(w,30,30,wa_Ferrel_Tropopause,wa_Ferrel_SL); //lat: 60
    init_v_or_w(v,45,45,va_Ferrel_Tropopause_45,va_Ferrel_SL_45); //lat: 45   

    //southern Ferrel cell
    init_u(u, 150, 150, ua_60); //lat: -60 

    init_v_or_w(v,150,150,va_Ferrel_Tropopause,va_Ferrel_SL); //lat: -60
    init_v_or_w(w,150,150,wa_Ferrel_Tropopause,wa_Ferrel_SL); //lat: -60
    init_v_or_w(v,135,135,va_Ferrel_Tropopause_45,va_Ferrel_SL_45); //lat: -45   
 
    // Hadley cell
    double ua_30 = 1.;
    double va_Hadley_SL = .25;
    double va_Hadley_Tropopause = - 1.;
    double va_Hadley_SL_15 = 1.;
    double va_Hadley_Tropopause_15 = - 1.;
    double wa_Hadley_SL = 1.;            // at surface
    double wa_Hadley_Tropopause = 30.;  // subtropic jet in m/s compares to 108 km/h

    //northern Hadley cell
    init_u(u, 60, 60, ua_30); //lat: 30 

    init_v_or_w(v,60,60,va_Hadley_Tropopause,va_Hadley_SL); //lat: 30
    init_v_or_w(w,60,60,wa_Hadley_Tropopause,wa_Hadley_SL); //lat: 30
    init_v_or_w(v,75,75,va_Hadley_Tropopause_15,va_Hadley_SL_15); //lat: 15   

    //southern Hadley cell
    init_u(u, 120, 120, ua_30); //lat: -30 

    init_v_or_w(v,120,120,va_Hadley_Tropopause,va_Hadley_SL); //lat: -30
    init_v_or_w(w,120,120,wa_Hadley_Tropopause,wa_Hadley_SL); //lat: -30
    init_v_or_w(v,105,105,va_Hadley_Tropopause_15,va_Hadley_SL_15); //lat: -15 

    // forming diagonals 
    //northen hemisphere
    form_diagonals(u, 0, 30);
    form_diagonals(w, 0, 30);
    form_diagonals(v, 0, 15);
    form_diagonals(v, 15, 30);

    form_diagonals(u, 30, 60);
    form_diagonals(w, 30, 60);
    form_diagonals(v, 30, 45);
    form_diagonals(v, 45, 60);

    form_diagonals(u, 60, 90);
    form_diagonals(w, 60, 90);
    form_diagonals(v, 60, 75);
    form_diagonals(v, 75, 90);

    //southen hemisphere
    form_diagonals(u, 90, 120);
    form_diagonals(w, 90, 120);
    form_diagonals(v, 90, 105);
    form_diagonals(v, 105, 120);

    form_diagonals(u, 120, 150);
    form_diagonals(w, 120, 150);
    form_diagonals(v, 120, 135);
    form_diagonals(v, 135, 150);

    form_diagonals(u, 150, 180);
    form_diagonals(w, 150, 180);
    form_diagonals(v, 150, 165);
    form_diagonals(v, 165, 180);

    //change the direction for southen hemisphere
    /*for ( int i = 0; i < im; i++ ){
        for ( int j = 91; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                v.x[ i ][ j ][ k ] = - v.x[ i ][ j ][ k ];
            }
        }
    }*/

    //smoothing transitions from cell to cell
    /*smooth_transition(u,v,w,60); 
    smooth_transition(u,v,w,30);    
    smooth_transition(u,v,w,75);
    smooth_transition(u,v,w,45);

    smooth_transition(u,v,w,90); 

    smooth_transition(u,v,w,120);
    smooth_transition(u,v,w,150);
    smooth_transition(u,v,w,105);
    smooth_transition(u,v,w,135);
    */
    // non dimensionalization by u_0
    for ( int i = 0; i < im; i++ ){
        for ( int k = 0; k < km; k++ ){
            for ( int j = 0; j < jm; j++ ){
                if ( is_land ( h, i, j, k ) )     
                    u.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] = 0.;
                else{
                    u.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ] / u_0;
                    v.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] / u_0;
                    w.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] / u_0;
                }
            }
        }
    }
}
/*
*
*/
void cAtmosphereModel::smooth_transition(Array &u, Array &v, Array &w, int lat){
    int start = lat-3, end = lat+3;
    for ( int i = 0; i < im; i++ ){
        for ( int k = 0; k < km; k++ ){
            for ( int j = start; j <= end; j++ ){
                u.x[ i ][ j ][ k ] = ( u.x[ i ][ end ][ k ] - u.x[ i ][ start ][ k ] ) *
                    ( ( double ) ( j - start )  / ( end - start ) ) + u.x[ i ][ start ][ k ];
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ end ][ k ] - v.x[ i ][ start ][ k ] ) *
                    ( ( double ) ( j - start )  / ( end - start ) ) + v.x[ i ][ start ][ k ];
                w.x[ i ][ j ][ k ] = ( w.x[ i ][ end ][ k ] - w.x[ i ][ start ][ k ] ) *
                    ( ( double ) ( j - start )  / ( end - start ) ) + w.x[ i ][ start ][ k ];
            }
        }
    }
}

/*
*
*/
void cAtmosphereModel::form_diagonals(Array &a, int start, int end){
    for ( int k = 0; k < km; k++ ){
        for ( int j = start; j < end; j++ ){
            for ( int i = 1; i < im; i++ ){
                a.x[ i ][ j ][ k ] = ( a.x[ i ][ end ][ k ] - a.x[ i ][ start ][ k ] ) *
                    ( j - start ) / (double)(end - start) + a.x[ i ][ start ][ k ];
            }
        }
    }
}

/*
*
*/
void  cAtmosphereModel::init_u(Array &u, int lat_1, int lat_2, double coefficient){
    for(int j = lat_1; j < lat_2+1; j++ ){
        int tropopause_layer = get_tropopause_layer(j);             
        double tropopause_height = get_layer_height(tropopause_layer);
        for( int k = 0; k < km; k++ ){
            for( int i = 0; i < tropopause_layer; i++ ){
                u.x[ i ][ j ][ k ] = -coefficient * 
                    parabola_interp(-1, 0, get_layer_height(i)*2/tropopause_height);
            }
        }
    }
}

/*
*
*/
void  cAtmosphereModel::init_v_or_w(Array &v_or_w, int lat_1, int lat_2, double coeff_trop, double coeff_sl){
    for(int j = lat_1; j < lat_2+1; j++ ){
        int tropopause_layer = get_tropopause_layer(j);
        double tropopause_height = get_layer_height(tropopause_layer);
        for( int k = 0; k < km; k++ ){
            if(is_ocean_surface(h, 0, j, k))
            {
                coeff_sl = v_or_w.x[0][j][k];
            }
            for( int i = 0; i < tropopause_layer; i++ ){
                v_or_w.x[ i ][ j ][ k ] = ( coeff_trop - coeff_sl ) *
                    get_layer_height(i)/tropopause_height + coeff_sl;
            }
        }
    }
    init_v_or_w_above_tropopause(v_or_w, lat_1, lat_2, coeff_trop);
}

/*
*
*/
void  cAtmosphereModel::init_v_or_w_above_tropopause(Array &v_or_w, int lat_1, int lat_2, double coeff){
    for(int j = lat_1; j < lat_2+1; j++ ){
        int tropopause_layer = get_tropopause_layer(j);
        if(tropopause_layer >= im-1) return;
        double tropopause_height = get_layer_height(tropopause_layer);
        for( int k = 0; k < km; k++ ){
            for( int i = tropopause_layer; i < im; i++ ){
                v_or_w.x[ i ][ j ][ k ] = coeff * (get_layer_height(im-1) - get_layer_height(i)) / 
                    (get_layer_height(im-1) - tropopause_height);
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
    if( stat( path.c_str(), &info ) != 0 ){
        mkdir(path.c_str(), 0777);
    }
    bool last_iter = false;
    std::ostringstream ss;
    if(iter_cnt_3d == pressure_iter_max * velocity_iter_max + 1){
        last_iter = true;
    }
    if(last_iter)
    {
        ss << "_time_" << (int)(*get_current_time()) << "_iter_n";
    }
    else
    {
        ss << "_time_" << (int)(*get_current_time()) << "_iter_" << iter_cnt_3d;
    }
    std::string postfix_str = ss.str();

    Array t_t(im, jm, km, 0),  v_t(im, jm, km, 0), w_t(im, jm, km, 0), m_t(im, jm, km, 0), u_t(im, jm, km, 0);
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                //m_t.x[ i ][ j ][ k ] = sqrt ( pow ( v.x[ i ][ j ][ k ] * u_0, 2 ) + pow ( w.x[ i ][ j ][ k ] * u_0, 2 ) );
                t_t.x[ i ][ j ][ k ] = t.x[ i ][ j ][ k ] * t_0 - t_0;
                v_t.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] * u_0;
                w_t.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] * u_0;

                if(last_iter){
                    if(fabs(v_t.x[ i ][ j ][ k ]) > 0 && !(fabs(v_t.x[ 0 ][ j ][ k ]) > 0))
                        v_t.x[ 0 ][ j ][ k ] = v_t.x[ i ][ j ][ k ];
                    if(fabs(w_t.x[ i ][ j ][ k ]) > 0 && !(fabs(w_t.x[ 0 ][ j ][ k ]) > 0))
                        w_t.x[ 0 ][ j ][ k ] = w_t.x[ i ][ j ][ k ];
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
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            for ( int i = im-2; i >= 0; i-- ){
                if ( is_land ( h, i, j, k ) ){
                    u.x[ i ][ j ][ k ] = 0.;
                    v.x[ i ][ j ][ k ] = 0.;
                    w.x[ i ][ j ][ k ] = 0.;
                    //t.x[ i ][ j ][ k ] = 1.;  // = 273.15 K
                    //c.x[ i ][ j ][ k ] = c_tropopause;  // = 1 g/kg water vapour
                    //c.x[ i ][ j ][ k ] = 0.; 
                    cloud.x[ i ][ j ][ k ] = 0.;
                    ice.x[ i ][ j ][ k ] = 0.;
//                    co2.x[ i ][ j ][ k ] = 1.;  // = 280 ppm
                    p_dyn.x[ i ][ j ][ k ] = 0.;
                }// is_land
            } // i
        } // k
    } // j
}

/*
*
*/
void cAtmosphereModel::vegetationDistribution(){
    // description or vegetation areas following the local dimensionsles values of precipitation, maximum value is 1
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( max_Precipitation > 0 && is_land( h, 0, j, k ) && !(t.x[ 0 ][ j ][ k ] < 1.) ){
                Vegetation.y[ j ][ k ] = Precipitation.y[ j ][ k ] / max_Precipitation; // actual vegetation areas
            }else{
                Vegetation.y[ j ][ k ] = 0.;
            }
        }
    }
}

void cAtmosphereModel::store_intermediate_data_2D(float coeff)
{
    for ( int j = 0; j < jm; j++ )
    {
        for ( int k = 0; k < km; k++ )
        {
            vn.x[ 0 ][ j ][ k ] = coeff * v.x[ 0 ][ j ][ k ];
            wn.x[ 0 ][ j ][ k ] = coeff * w.x[ 0 ][ j ][ k ];
            p_dynn.x[ 0 ][ j ][ k ] = coeff * p_dyn.x[ 0 ][ j ][ k ];
        }
    }
}

void cAtmosphereModel::store_intermediate_data_3D(float coeff)
{ 
    for ( int i = 0; i < im; i++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int k = 0; k < km; k++ )
            {
                un.x[ i ][ j ][ k ] = coeff * u.x[ i ][ j ][ k ];
                vn.x[ i ][ j ][ k ] = coeff * v.x[ i ][ j ][ k ];
                wn.x[ i ][ j ][ k ] = coeff * w.x[ i ][ j ][ k ];
                p_dynn.x[ i ][ j ][ k ] = coeff * p_dyn.x[ i ][ j ][ k ];
                tn.x[ i ][ j ][ k ] = coeff * t.x[ i ][ j ][ k ];
                cn.x[ i ][ j ][ k ] = coeff * c.x[ i ][ j ][ k ];
                cloudn.x[ i ][ j ][ k ] = coeff * cloud.x[ i ][ j ][ k ];
                icen.x[ i ][ j ][ k ] = coeff * ice.x[ i ][ j ][ k ];
                co2n.x[ i ][ j ][ k ] = coeff * co2.x[ i ][ j ][ k ];

            }
        }
    }
}

void cAtmosphereModel::adjust_temperature_IC(double** t, int jm, int km)
{
    for(int k=0; k < km; k++ ){
        for(int j=0; j < jm; j++ ){
            t[ j ][ k ] = temperature_NASA.y[ j ][ k ] = ( t[ j ][ k ] + t_0 ) / t_0;
        }
    }

    // correction of surface temperature around 180°E
    int k_half = ( km -1 ) / 2;
    for ( int j = 0; j < jm; j++ ){
        t[ j ][ k_half ] = ( t[ j ][ k_half + 1 ] + t[ j ][ k_half - 1 ] ) / 2.;
        temperature_NASA.y[ j ][ k_half ] = ( temperature_NASA.y[ j ][ k_half + 1 ] +
            temperature_NASA.y[ j ][ k_half - 1 ] ) / 2.;
    }

}

void cAtmosphereModel::check_data(Array& a, Array&an, const std::string& name){
    float t_diff_min, t_diff_max, t_diff_mean, t_min, t_max;
    for ( int i = 0; i < im; i++ )
    {
        t_diff_min = t_diff_max = t_diff_mean = 0; 
        t_min = t_max = a.x[i][0][0];
        for ( int j = 0; j < jm; j++ )
        {
            for ( int k = 0; k < km; k++ )
            {
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
        logger() << name << " diff mean: " << t_diff_mean << "   " << (double)t_diff_mean / (jm*km) << std::endl;
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
