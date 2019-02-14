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

#include "BC_Atm.h"
#include "BC_Bath_Atm.h"
#include "BC_Thermo.h"
#include "Accuracy_Atm.h"
#include "RHS_Atm.h"
#include "RungeKutta_Atm.h"
#include "PostProcess_Atm.h"
#include "Pressure_Atm.h"
#include "Results_Atm.h"
#include "MinMax_Atm.h"
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

cAtmosphereModel::cAtmosphereModel() :
    i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    im_tropopause(NULL),
    is_node_weights_initialised(false), 
    old_arrays_3d {&u,  &v,  &w,  &t,  &p_dyn,  &c,  &cloud,  &ice,  &co2 },
    new_arrays_3d {&un, &vn, &wn, &tn, &p_dynn, &cn, &cloudn, &icen, &co2n},
    old_arrays_2d {&v,  &w,  &p_dyn }, 
    new_arrays_2d {&vn, &wn, &p_dynn},
    residuum_2d(1, jm, km),
    residuum_3d(im, jm, km)
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

    im_tropopause = new int [ jm ];// location of the tropopaus

    emin = epsres * 100.;
    
    m_model = this;

    load_temperature_curve();

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

    delete [] im_tropopause;
    im_tropopause = NULL;
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
//    logger() << "RunTimeSlice: " << Ma << " Ma"<< std::endl <<std::endl;

    reset_arrays();    

    m_current_time = m_time_list.insert(float(Ma)).first;

    struct stat info;
    if( stat( output_path.c_str(), &info ) != 0 ){
        mkdir(output_path.c_str(), 0777);
    }
    // maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km )
    // maximum number of overall iterations ( n )
    // maximum number of inner velocity loop iterations ( velocity_iter_max )
    // maximum number of outer pressure loop iterations ( pressure_iter_max )

    cout.precision ( 6 );
    cout.setf ( ios::fixed );

    //  initial values for the number of computed steps and the time
    double t_cretaceous = 0.;

    //Prepare the temperature and precipitation data file
    string Name_SurfaceTemperature_File  = temperature_file;
    string Name_SurfacePrecipitation_File = precipitation_file;

    if(Ma != 0 && use_earthbyte_reconstruction){
        Name_SurfaceTemperature_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";
        Name_SurfacePrecipitation_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";    
    
        if( stat( Name_SurfaceTemperature_File.c_str(), &info ) != 0 || 
            stat( Name_SurfacePrecipitation_File.c_str(), &info ) != 0 )
        {
            std::string cmd_str = "python " + reconstruction_script_path + " " + std::to_string(Ma - time_step) + " " + 
                std::to_string(Ma) + " " + output_path + " " + BathymetrySuffix + " atm";
            int ret = system(cmd_str.c_str());
            std::cout << " reconstruction script returned: " << ret << std::endl;
        } 
    }

    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    bathymetry_filepath = bathymetry_path + "/" + bathymetry_name;

    cout << "\n   Output is being written to " << output_path << "\n";
    cout << "   Ma = " << Ma << "\n";
    cout << "   bathymetry_path = " << bathymetry_path << "\n";
    cout << "   bathymetry_filepath = " << bathymetry_filepath << "\n\n";

    if (verbose) {
        cout << endl << endl << endl;
        cout << "***** Atmosphere General Circulation Model ( AGCM ) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 4 additional transport equations to describe the water vapour, cloud water, cloud ice and co2 concentration" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** water vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code ( German Weather Forecast )" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** co2 concentration appears in the absorptivity of the radiation models" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )" << endl << endl;

        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
    }

    //  initialization of the bathymetry/topography

    //  class BC_Bathymetry_Atmosphere for the geometrical boundary condition of the computational area
    BC_Bathymetry_Atmosphere LandArea(this, NASATemperature, im, jm, km, co2_vegetation, co2_land, co2_ocean);

    //  topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
    //LandArea.BC_MountainSurface ( bathymetry_filepath, Topography, h );
    init_topography(bathymetry_filepath);
    //  class element for the computation of the ratio ocean to land areas, also supply and removal of CO2 on land, ocean and by vegetation
    LandArea.land_oceanFraction ( h );

    //  class calls for the solution of the flow properties

    //  class BC_Atmosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
    BC_Atmosphere  boundary ( im, jm, km, t_tropopause );

    //  class RHS_Atmosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
    //  class RungeKutta_Atmosphere for the explicit solution of the Navier-Stokes equations

    //  class Results_MSL_Atm to compute and show results on the mean sea level, MSL
    Results_MSL_Atm  calculate_MSL ( im, jm, km, sun, g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_equator, lv, ls, 
                                     cp_l, L_atm, dt, dr, dthe, dphi, r_air, R_Air, r_water_vapour, R_WaterVapour, 
                                     co2_vegetation, co2_ocean, co2_land, gam, t_pole, t_cretaceous, t_average );

    //  class Pressure for the subsequent computation of the pressure by a separate Euler equation
    Pressure_Atm  startPressure ( im, jm, km, dr, dthe, dphi );

    //  class BC_Thermo for the initial and boundary conditions of the flow properties
    BC_Thermo  circulation (this, im, jm, km, h ); 

//    t.printArray ( im, jm, km );

    //  class element calls for the preparation of initial conditions for the flow properties

    //  class element for the tropopause location as a parabolic distribution from pole to pole 
    circulation.TropopauseLocation ();

    //  class element for the initial conditions for u-v-w-velocity components
    //circulation.IC_CellStructure ( h, u, v, w );
    init_velocities();

    //  class element for the surface temperature from NASA for comparison
    //  if ( Ma == 0 ) circulation.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, temperature_NASA, t );
    circulation.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, temperature_NASA, t );

    //  class element for the surface precipitation from NASA for comparison
    circulation.BC_Surface_Precipitation_NASA ( Name_SurfacePrecipitation_File, precipitation_NASA );

    //  class element for the parabolic temperature distribution from pol to pol, maximum temperature at equator
    init_temperature();

    //  class element for the correction of the temperature initial distribution around coasts
    if ( ( NASATemperature == 1 ) && ( Ma > 0 ) && !use_earthbyte_reconstruction) 
    {
        circulation.IC_Temperature_WestEastCoast ( h, t );
    }

    //  class element for the surface pressure computed by surface temperature with gas equation
    circulation.BC_Pressure ( p_stat, p_dyn, t, h );

    //parabolic water vapour distribution from pol to pol, maximum water vapour volume at equator
    //circulation.BC_WaterVapour ( h, p_stat, t, c, v, w );
    init_water_vapour();

    //  class element for the parabolic CO2 distribution from pol to pol, maximum CO2 volume at equator
    circulation.BC_CO2 ( Vegetation, h, t, p_dyn, co2 );

    // class element for the surface temperature computation by radiation flux density
    if ( RadiationModel == 1 ){
        BC_Radiation_multi_layer(); 
    }

    // class element for the storing of velocity components, pressure and temperature for iteration start
    move_data_to_new_arrays(im, jm, km, 1., old_arrays_3d, new_arrays_3d);
    move_data_to_new_arrays(jm, km, 1., old_arrays_2d, new_arrays_2d);



    // ***********************************   start of pressure and velocity iterations ***********************************

    run_2D_loop(boundary, LandArea, startPressure, circulation);
    
    cout << endl << endl;

    run_3D_loop( boundary, LandArea, startPressure, calculate_MSL, circulation);

    cout << endl << endl;

    restrain_temperature();

    //write the ouput files
    write_file(bathymetry_name, output_path, true);

    iter_cnt_3d++;
    save_data();    

    //  final remarks
    cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;
    if ( emin <= epsres ){
        cout << "***** steady solution reached! *****" << endl;
    }

    if(debug){
        fedisableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO); //not platform independent(bad, very bad, I know)
    }
}



void cAtmosphereModel::Run() 
{
    auto start_time = std::chrono::system_clock::now();
    std::time_t start_time_t = std::chrono::system_clock::to_time_t(start_time);
    logger() << "Start Time:" << std::ctime(&start_time_t) << std::endl;
    //char buffer[80];
    //strftime(buffer, 80, "%c", localtime(&start_time_t));
    //timestamp = string(buffer);
    //std::cout << timestamp << std::endl;

    mkdir(output_path.c_str(), 0777);

    cout << std::endl << "Output is being written to " << output_path << std::endl << std::endl;

    if (verbose) {
        cout << endl << endl << endl;
        cout << "***** Atmosphere General Circulation Model ( AGCM ) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 4 additional transport equations to describe the water vapour, cloud water, cloud ice and co2 concentration" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** water vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code ( German Weather Forecast )" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** co2 concentration appears in the absorptivity of the radiation models" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )" << endl << endl;

        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
    }

    for(int i = time_start; i <= time_end; i+=time_step)
    {
        RunTimeSlice(i);
    }

    //  final remarks
    cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
    cout << "\n\n\n\n";

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



void cAtmosphereModel::print_min_max_values()
{
    MinMax_Atm min_max_3d( im, jm, km );

    //  searching of maximum and minimum values of temperature
    min_max_3d.searchMinMax_3D( " max 3D temperature ", " min 3D temperature ", "°C", t, h, 273.15, 
                                [](double i)->double{return i - 273.15;},
                                true );

    //  searching of maximum and minimum values of u-component
    min_max_3d.searchMinMax_3D ( " max 3D u-component ", " min 3D u-component ", "m/s", u, h, u_0);

    //  searching of maximum and minimum values of v-component
    min_max_3d.searchMinMax_3D ( " max 3D v-component ", " min 3D v-component ", "m/s", v, h, u_0 );

    //  searching of maximum and minimum values of w-component
    min_max_3d.searchMinMax_3D ( " max 3D w-component ", " min 3D w-component ", "m/s", w, h, u_0 );

    //  searching of maximum and minimum values of dynamic pressure
    min_max_3d.searchMinMax_3D ( " max 3D pressure dynamic ", " min 3D pressure dynamic ", "hPa", p_dyn, h, 0.768 ); // 0.768 = 0.01 * r_air *u_0*u_0 in hPa

    //  searching of maximum and minimum values of static pressure
    min_max_3d.searchMinMax_3D ( " max 3D pressure static ", " min 3D pressure static ", "hPa", p_stat, h );

    cout << endl << " energies in the three dimensional space: " << endl << endl;

    //  searching of maximum and minimum values of radiation_3D
    min_max_3d.searchMinMax_3D ( " max 3D radiation ",  " min 3D radiation ",  "W/m2", radiation_3D, h );

    //  searching of maximum and minimum values of sensible heat
    min_max_3d.searchMinMax_3D ( " max 3D sensible heat ", " min 3D sensible heat ", "W/m2", Q_Sensible, h );

    //  searching of maximum and minimum values of latency
    min_max_3d.searchMinMax_3D ( " max 3D latent heat ", " min 3D latent heat ", "W/m2", Q_Latent, h );

    cout << endl << " greenhouse gases: " << endl << endl;

    //  searching of maximum and minimum values of water vapour
    min_max_3d.searchMinMax_3D ( " max 3D water vapour ",  " min 3D water vapour ", "g/kg", c, h, 1000. );

    //  searching of maximum and minimum values of cloud water
    min_max_3d.searchMinMax_3D ( " max 3D cloud water ", " min 3D cloud water ", "g/kg", cloud, h, 1000. );

    //  searching of maximum and minimum values of cloud ice
    min_max_3d.searchMinMax_3D ( " max 3D cloud ice ", " min 3D cloud ice ", "g/kg", ice, h, 1000. );

    //  searching of maximum and minimum values of rain precipitation
    min_max_3d.searchMinMax_3D ( " max 3D rain ", " min 3D rain ", "mm/d", P_rain, h, 8.46e4 );

    //  searching of maximum and minimum values of snow precipitation
    min_max_3d.searchMinMax_3D (  " max 3D snow ", " min 3D snow ", "mm/d", P_snow, h, 8.46e4 );

    //  searching of maximum and minimum values of co2
    min_max_3d.searchMinMax_3D ( " max 3D co2 ", " min 3D co2 ", "ppm", co2, h, 280. );

    //  searching of maximum and minimum values of epsilon
    min_max_3d.searchMinMax_3D ( " max 3D epsilon ",  " min 3D epsilon ", "%", epsilon_3D, h );

    //  searching of maximum and minimum values of buoyancy force
    min_max_3d.searchMinMax_3D (  " max 3D buoyancy force ", " min 3D buoyancy force ", "kN/m2", BuoyancyForce, h );



    // 2D-fields

    cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" << endl << 
        " results based on two dimensional considerations of the problem" << endl;

    cout << endl << " co2 distribution row-wise: " << endl << endl;

    MinMax_Atm  min_max_2d ( jm, km );

    //  searching of maximum and minimum values of co2 total
    min_max_2d.searchMinMax_2D ( " max co2_total ", " min co2_total ", " ppm ", co2_total, h, 280. );

    cout << endl << " precipitation: " << endl << endl;

    //  searching of maximum and minimum values of precipitation
    min_max_2d.searchMinMax_2D (  " max precipitation ", " min precipitation ", "mm/d", Precipitation, h, 1. );
    max_Precipitation = min_max_2d.out_maxValue (  );

    //  searching of maximum and minimum values of precipitable water
    min_max_2d.searchMinMax_2D ( " max precipitable water ", " min precipitable water ", "mm", 
                                 precipitable_water, h, 1. );

    cout << endl << " energies at see level without convection influence: " << endl << endl;

    //  searching of maximum and minimum values of radiation
    min_max_2d.searchMinMax_2D ( " max 2D Q radiation ", " min 2D Q radiation ",  "W/m2", Q_radiation, h );

    //  searching of maximum and minimum values of latent energy
    min_max_2d.searchMinMax_2D ( " max 2D Q latent ", " min 2D Q latent ", "W/m2", Q_latent, h );

    //  searching of maximum and minimum values of sensible energy
    min_max_2d.searchMinMax_2D ( " max 2D Q sensible ", " min 2D Q sensible ", "W/m2", Q_sensible, h );

    //  searching of maximum and minimum values of bottom heat
    min_max_2d.searchMinMax_2D ( " max 2D Q bottom ", " min 2D Q bottom heat ", "W/m2", Q_bottom, h );

    cout << endl << " secondary data: " << endl << endl;

    //  searching of maximum and minimum values of Evaporation
    min_max_2d.searchMinMax_2D ( " max heat Evaporation ", " min heat Evaporation ", " kJ/kg", Q_Evaporation, h );

    //  searching of maximum and minimum values of Evaporation by Dalton
    min_max_2d.searchMinMax_2D ( " max Evaporation Dalton ", " min Evaporation Dalton ", "mm/d", Evaporation_Dalton, h );

    //  searching of maximum and minimum values of Evaporation by Penman
    min_max_2d.searchMinMax_2D ( " max Evaporation Penman ", " min Evaporation Penman ", "mm/d", Evaporation_Penman, h );

    cout << endl << " properties of the atmosphere at the surface: " << endl << endl;

    //  searching of maximum and minimum values of albedo
    min_max_2d.searchMinMax_2D (  " max 2D albedo ", " min 2D albedo ", "%", albedo, h );

    //  searching of maximum and minimum values of epsilon
    min_max_2d.searchMinMax_2D ( " max 2D epsilon ", " min 2D epsilon ", "%", epsilon, h );

    //  searching of maximum and minimum values of topography
    min_max_2d.searchMinMax_2D ( " max 2D topography ", " min 2D topography ", "m", Topography, h );
}


void cAtmosphereModel::write_file(std::string &bathymetry_name, std::string &output_path, bool is_final_result){
    int Ma = int(round(*get_current_time()));
    //  Printout:

    //  printout in ParaView files and sequel files

    //  class PostProcess_Atmosphaere for the printing of results
    PostProcess_Atmosphere write_File ( im, jm, km, output_path );

    //  writing of data in ParaView files
    //  radial data along constant hight above ground
    int i_radial = 0;
    //  int i_radial = 10;
/*    write_File.paraview_vtk_radial ( bathymetry_name, Ma, i_radial, iter_cnt-1, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, 
                                     BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, radiation_3D, 
                                     Q_Latent, Q_Sensible, epsilon_3D, P_rain, P_snow, precipitable_water, Q_bottom, 
                                     Q_radiation, Q_latent, Q_sensible, Evaporation_Penman, Evaporation_Dalton, 
                                     Q_Evaporation, temperature_NASA, precipitation_NASA, Vegetation, albedo, epsilon, 
                                     Precipitation, Topography, temp_NASA );

    //  londitudinal data along constant latitudes
    int j_longal = 62;          // Mount Everest/Himalaya
    write_File.paraview_vtk_longal ( bathymetry_name, j_longal, iter_cnt-1, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, 
                                     BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Q_Latent, 
                                     Q_Sensible, epsilon_3D, P_rain, P_snow );

    int k_zonal = 87;           // Mount Everest/Himalaya
    write_File.paraview_vtk_zonal ( bathymetry_name, k_zonal, iter_cnt-1, hp, ep, R_Air, g, L_atm, u_0, t_0, p_0, 
                                    r_air, c_0, co2_0, 
                                    h, p_dyn, p_stat, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, 
                                    Q_Latent, Q_Sensible, radiation_3D, epsilon_3D, P_rain, P_snow, S_v, S_c, S_i, S_r, 
                                    S_s, S_c_c );

    //  3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
    if(paraview_panorama_vts) //This function creates a large file. Use a flag to control if it is wanted.
    {
        write_File.paraview_panorama_vts ( bathymetry_name, iter_cnt-1, u_0, t_0, p_0, r_air, c_0, co2_0, h, t, p_dyn, p_stat, 
                                           BuoyancyForce, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Q_Latent, 
                                           Q_Sensible, epsilon_3D, P_rain, P_snow );
    }
*/
    //  writing of v-w-data in the v_w_transfer file
    PostProcess_Atmosphere ppa ( im, jm, km, output_path );
    ppa.Atmosphere_v_w_Transfer ( bathymetry_name, u_0, v, w, t, p_dyn, Evaporation_Dalton, Precipitation );
    ppa.Atmosphere_PlotData ( bathymetry_name, (is_final_result ? -1 : iter_cnt-1), u_0, t_0, h, v, w, t, c, 
                              Precipitation, precipitable_water );

    if(debug){
        ppa.save(output_path+"/residuum_"+std::to_string(iter_cnt-1)+".dat", 
                std::vector<std::string>{"residuum"},
                std::vector<Vector3D<>* >{&residuum_3d},
                1);
    }
}

void cAtmosphereModel::run_2D_loop( BC_Atmosphere &boundary,
                                    BC_Bathymetry_Atmosphere &LandArea, 
                                    Pressure_Atm &startPressure, BC_Thermo &circulation){
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
                boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
                boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
                
                //  old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
                Accuracy_Atm        min_Residuum_2D ( im, jm, km, dthe, dphi );
                double residuum_old = std::get<0>(min_Residuum_2D.residuumQuery_2D ( rad, the, v, w , residuum_2d));

                circulation.Value_Limitation_Atm ( h, u, v, w, p_dyn, t, c, cloud, ice, co2 );

                LandArea.BC_SolidGround ( RadiationModel, Ma, g, hp, ep, r_air, R_Air, t_0, c_0, t_land, t_cretaceous, 
                                          t_equator, t_pole, 
                                          t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, 
                                          pa, gam, sigma, h, u, v, w, t, p_dyn, c, cloud, ice, co2, 
                                          radiation_3D, Vegetation );
                
                //  class RungeKutta for the solution of the differential equations describing the flow properties
                solveRungeKutta_2D_Atmosphere();
                
                //  new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
                double residuum = std::get<0>(min_Residuum_2D.residuumQuery_2D ( rad, the, v, w, residuum_2d ));

                emin = fabs ( ( residuum - residuum_old ) / residuum_old );
                
                //  state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
                min_Residuum_2D.steadyQuery_2D ( v, vn, w, wn, p_dyn, p_dynn );

                move_data_to_new_arrays(jm, km, 1., old_arrays_2d, new_arrays_2d);

                iter_cnt++;
            }
            //  ::::::   end of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::::::::::::::


            //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
            startPressure.computePressure_2D ( u_0, r_air, rad, the, p_dyn, p_dynn, h, aux_v, aux_w );

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


void cAtmosphereModel::run_3D_loop( BC_Atmosphere &boundary,
                                    BC_Bathymetry_Atmosphere &LandArea,
                                    Pressure_Atm &startPressure, Results_MSL_Atm &calculate_MSL,                  
                                    BC_Thermo &circulation){
    
    iter_cnt = 1;
    iter_cnt_3d = 0;
    emin = epsres * 100.;

    int Ma = int(round(*get_current_time()));

    move_data_to_new_arrays(im, jm, km, 1., old_arrays_3d, new_arrays_3d);

    /** ::::::::::::::   begin of 3D pressure loop : if ( pressure_iter > pressure_iter_max )   :::::::::::::::: **/
    for ( int pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++ )
    {
        /** ::::::::::::   begin of 3D velocity loop : if ( velocity_iter > velocity_iter_max )   ::::::::::::::::::: **/
        for ( int velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++ )
        {
            Array tmp = (t-1)*t_0;
            tmp.inspect();
            //  query to realize zero divergence of the continuity equation ( div c = 0 )
            cout << endl << endl;
            cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            cout << " 3D AGCM iterational process" << endl;
            cout << " max total iteration number nm = " << nm << endl << endl;

            cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum "
                << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl << 
                endl << " Ma = " << Ma << "     n = " << iter_cnt << "    velocity_iter_max = " << velocity_iter_max << 
                "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << 
                "    pressure_iter = " << pressure_iter << endl;

            //  old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
            Accuracy_Atm        min_Residuum ( im, jm, km, dr, dthe, dphi );
            double residuum_old = std::get<0>(min_Residuum.residuumQuery_3D ( rad, the, u, v, w, residuum_3d ));
            
            //logger() <<  residuum_3d(1, 30, 150) << " residuum_mchin" <<Ma<<std::endl;
            
            //  class BC_Atmosphaere for the geometry of a shell of a sphere
            boundary.BC_radius ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
            boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
            boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

            //Ice_Water_Saturation_Adjustment, distribution of cloud ice and cloud water dependent on water vapour amount and temperature
            if ( velocity_iter % 2 == 0 ){
                Ice_Water_Saturation_Adjustment();
            }

            circulation.Value_Limitation_Atm ( h, u, v, w, p_dyn, t, c, cloud, ice, co2 );

            LandArea.BC_SolidGround ( RadiationModel, Ma, g, hp, ep, r_air, R_Air, t_0, c_0, t_land, t_cretaceous, t_equator, 
                                      t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, 
                                      co2_tropopause, pa, gam, sigma, h, u, v, w, t, p_dyn, c, cloud, 
                                      ice, co2, radiation_3D, Vegetation );
            
            // class RungeKutta for the solution of the differential equations describing the flow properties
            solveRungeKutta_3D_Atmosphere(); 
            
            circulation.Value_Limitation_Atm ( h, u, v, w, p_dyn, t, c, cloud, ice, co2 );

            // class element for the surface temperature computation by radiation flux density
            if ( RadiationModel == 1 ){
                BC_Radiation_multi_layer(); 
            }
            //  new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
            double residuum = std::get<0>(min_Residuum.residuumQuery_3D ( rad, the, u, v, w, residuum_3d ));

            emin = fabs ( ( residuum - residuum_old ) / residuum_old );

            //  statements on the convergence und iterational process
            min_Residuum.steadyQuery_3D ( u, un, v, vn, w, wn, t, tn, c, cn, cloud, cloudn, ice, icen, co2, co2n, p_dyn, 
                                          p_dynn, L_atm);

            // 3D_fields

            //  class element for the initial conditions the latent heat
            circulation.Latent_Heat ( rad, the, phi, h, t, tn, u, v, w, p_dyn, p_stat, c, ice, Q_Latent, Q_Sensible, 
                                      radiation_3D, Q_radiation, Q_latent, Q_sensible, Q_bottom );

            print_min_max_values();

            //  computation of vegetation areas
            LandArea.vegetationDistribution ( max_Precipitation, Precipitation, Vegetation, t, h );


            //  composition of results
            run_MSL_data(); 

            //  Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
            //  resulting the precipitation distribution formed of rain and snow
            if ( velocity_iter % 2 == 0){
                circulation.Two_Category_Ice_Scheme ( h, c, t, p_stat, 
                                                      cloud, ice, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );
            }

            move_data_to_new_arrays(im, jm, km, 1., old_arrays_3d, new_arrays_3d);
            
            iter_cnt++;
            iter_cnt_3d++;
            save_data();
        }
        /**  ::::::::::::   end of velocity loop_3D: if ( velocity_iter > velocity_iter_max )   :::::::::::::::::::::::::::: **/
        
        //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
        startPressure.computePressure_3D ( u_0, r_air, rad, the, p_dyn, p_dynn, h, aux_u, aux_v, aux_w );
/*
        //  Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
        //  resulting the precipitation formed of rain and snow
        if(iter_cnt >= 2){
            circulation.Two_Category_Ice_Scheme ( h, c, t, p_stat, cloud, ice, 
                                              P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );
        }
*/
        //logger() << fabs ( p_dyn.x[ 20 ][ 30 ][ 150 ] - p_dynn.x[ 20 ][ 30 ][ 150 ] ) << " pressure_mchin" <<Ma<<std::endl;
        //logger() << std::get<0>(max_diff( im, jm, km, p_dyn, p_dynn)) << " pressure_max_diff" <<Ma<<std::endl;

        if( pressure_iter % checkpoint == 0 ){
            write_file(bathymetry_name, output_path);
        }

        //  limit of the computation in the sense of time steps
        if ( iter_cnt > nm )
        {
            cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
            break;
        }
    }
    /**  :::::   end of pressure loop_3D: if ( pressure_iter > pressure_iter_max )   ::::::::::::::::::::::::::::: **/
}


/*
*
*/
void cAtmosphereModel::load_temperature_curve()
{
    std::string line;
    std::ifstream f(temperature_curve_file);
    
    if (!f.is_open()){
        std::cout << "error while opening file: "<< temperature_curve_file << std::endl;
    }

    float time=0.,temperature=0;
    while(getline(f, line)) {
        std::stringstream(line) >> time >> temperature;
        m_temperature_curve.insert(std::pair<float,float>(time, temperature));
        //std::cout << time <<"  " <<temperature<< std::endl;
    }
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
    double tmp_1 = get_mean_temperature_from_curve(*get_current_time());
    double tmp_2 = calculate_mean_temperature();
    double diff = tmp_2 - tmp_1;
    for(int j=0;j<jm;j++){
        for(int k=0; k<km; k++){
            t.x[0][j][k] -= diff/t_0;
            if(t.x[0][j][k] > (1+40/t_0)){
                t.x[0][j][k] = 1 + 40/t_0;//don't allow temperature to exceed 40 degrees.
                //logger() << "temperature is restraint " << (t.x[0][j][k] - 1)*t_0 << std::endl;
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

    double t_u = 0.;  // not air but ocean surface temperature should be used
    double sat_difference = 0.;  // not air but ocean surface temperature
    // should be involved in water vapour saturation difference, it is not the saturation deficit
    double Dalton_Evaporation = 0;  // ocean surface evaporation

    // water vapour contents computed by Clausius-Clapeyron-formula
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            int i_mount = get_surface_layer(j, k);

            if ( is_air ( h, 0, j, k ) ){
                c.x[ i_mount ][ j ][ k ] = hp * ep *exp ( 17.0809 * ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + 
                    ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) ) ) / ( ( r_air * R_Air * t.x[ i_mount ][ j ][ k ] * t_0 ) * .01 );
                // saturation of relative water vapour in kg/kg
                c.x[ i_mount ][ j ][ k ] = c_ocean * c.x[ i_mount ][ j ][ k ];
                // relativ water vapour contents on ocean surface reduced by factor

                p_stat.x[ i_mount ][ j ][ k ] = ( r_air * R_Air * t.x[ i_mount ][ j ][ k ] * t_0 ) * .01;  // given in hPa
                t_u = t.x[ i_mount ][ j ][ k ] * t_0; // in K
                double r_dry = 100. * p_stat.x[ i_mount ][ j ][ k ] / ( R_Air * t.x[ i_mount ][ j ][ k ] * t_0 );
                double r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i_mount ][ j ][ k ] );
                double e = c.x[ i_mount ][ j ][ k ] * p_stat.x[ i_mount ][ j ][ k ] / ep;  // water vapour pressure in hPa
                double E = hp * exp_func ( t_u, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                double dt_dim = L_atm / u_0 * dt;// dimensional time step of the system in s == 0.02 s
                double dr_dim = dr * L_atm;  // = 0.025 * 16000 = 400

                sat_difference = ( E - e );  // saturation difference in hPa/K
                Dalton_Evaporation = 8.46e-4 * C_Dalton ( u_0, v.x[ i_mount ][ j ][ k ], w.x[ i_mount ][ j ][ k ] ) *
                    sat_difference * dt_dim / ( r_humid * dr_dim ) * 24.;  // mm/h in mm/d
                // Evaporation_Dalton given in mm/d,  recalculation in kg/kg water vapour/ dry air by
                // c_ev = 8.46e-4 * Dalton_Evaporation * dt_dim / ( r_humid * dr_dim ) [ mm/d * m³ * s / ( kg * m ) * ( kg/kg_air ) ] == [ kg/kg_air ]
                // [ mm/s ] == [ kg / ( m² * s ) ]

                c.x[ i_mount ][ j ][ k ] = Dalton_Evaporation + c.x[ i_mount ][ j ][ k ]; // kg/kg_air
            }
            else if ( is_land ( h, 0, j, k ) ){
                c.x[ i_mount ][ j ][ k ] = hp * ep * exp ( 17.0809 * ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + 
                    ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) ) ) / ( ( r_air * R_Air * t.x[ i_mount ][ j ][ k ] * t_0 ) * .01 );
                c.x[ i_mount ][ j ][ k ] = c_land * c.x[ i_mount ][ j ][ k ];
                // relativ water vapour contents on land reduced by factor
                // Dalton_Evaporation = 8.46e-4 * C_Dalton ( u_0, v.x[ i_mount + 1 ][ j ][ k ], w.x[ i_mount + 1 ][ j ][ k ] ) *
                // sat_difference * dt_dim / ( r_humid * dr_dim ) * 24.;  // mm/h in mm/d
                c.x[ i_mount ][ j ][ k ] = Dalton_Evaporation + c.x[ i_mount ][ j ][ k ]; // kg/kg_air
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
void cAtmosphereModel::init_topography(string &topo_filename){
    // default adjustment, h must be 0 everywhere
    h.initArray(im, jm, km, 0.);

    // reading data from file Name_Bathymetry_File_Read
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
    double co2_cretaceous = 3.2886 * pow ( ( t_cretaceous + t_average ), 2 ) - 32.8859 *
        ( t_cretaceous + t_average ) + 102.2148;  // in ppm
    double co2_average = 3.2886 * pow ( t_average, 2 ) - 32.8859 * t_average + 102.2148;  // in ppm
    co2_cretaceous = co2_cretaceous - co2_average;

    //use parabolic distribution from pole to pole to initialize co2  
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            int i_mount = get_surface_layer(j, k);

            if ( is_water ( h, 0, j, k ) ){
                co2.x[ 0 ][ j ][ k ] = parabola_interp(co2_equator, co2_pole, 2*j/(jm-1)) + co2_cretaceous + co2_ocean; 
            }else{
                co2.x[ i_mount ][ j ][ k ] = parabola_interp(co2_equator, co2_pole, 2*j/(jm-1)) + 
                    + co2_cretaceous + co2_land - co2_vegetation * Vegetation.y[ j ][ k ];
            }
            co2.x[ 0 ][ j ][ k ] /= co2_0;// non-dimensional
        }
    }

    // co2 distribution decreasing approaching tropopause, above no co2
    for ( int j = 0; j < jm; j++ ){
        int i_trop = get_tropopause_layer(j);
        for ( int k = 0; k < km; k++ ){
            int i_mount = get_surface_layer(j, k);
            
            for ( int i = 1; i < im; i++ ){
                if ( i < i_trop ){
                    if(i>i_mount){
                        double x = (get_layer_height(i) - get_layer_height(i_mount)) / 
                            (get_layer_height(i_trop) - get_layer_height(i_mount));
                        co2.x[ i ][ j ][ k ] = parabola_interp(co2_tropopause, co2.x[ i_mount ][ j ][ k ], x); 
                    }else{
                        c.x[ i ][ j ][ k ] = c.x[ i_mount ][ j ][ k ];
                    }
                }else{
                    c.x[ i ][ j ][ k ] = c_tropopause;
                }
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
    double ua_90 = - 0.5;
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
    double ua_30 = - 1.;
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
    for ( int i = 0; i < im; i++ ){
        for ( int j = 91; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                v.x[ i ][ j ][ k ] = - v.x[ i ][ j ][ k ];
            }
        }
    }

    //smoothing transitions from cell to cell
    smooth_transition(u,v,w,60); 
    smooth_transition(u,v,w,30);    
    smooth_transition(u,v,w,75);
    smooth_transition(u,v,w,45);

    smooth_transition(u,v,w,90); 

    smooth_transition(u,v,w,120);
    smooth_transition(u,v,w,150);
    smooth_transition(u,v,w,105);
    smooth_transition(u,v,w,135);

    // non dimensionalization by u_0
    for ( int i = 0; i < im; i++ ){
        for ( int k = 0; k < km; k++ ){
            for ( int j = 0; j < jm; j++ ){
                u.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ] / u_0;
                v.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] / u_0;
                w.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] / u_0;
                if ( is_land ( h, i, j, k ) )     u.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] = 0.;
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
            for ( int i = 0; i < im; i++ ){
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

void cAtmosphereModel::save_array(const string& fn, const Array& a){
    for(int i=0; i<im; i++){
        std::ofstream os(fn + "_" + to_string(i) + ".bin", std::ios::binary | std::ios::out);
        for(int j=jm-1; j>=0; j--){
            os.write(reinterpret_cast<const char*>(a.x[i][j]+(km/2)), std::streamsize((km/2)*sizeof(double)));
            os.write(reinterpret_cast<const char*>(a.x[i][j]), std::streamsize((km/2+1)*sizeof(double)));
        }
        os.close();
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
    std::ostringstream ss;
    ss << "_" << (int)(*get_current_time()) << "_" << iter_cnt_3d /*<< "_" << timestamp << ".bin"*/;
    std::string postfix_str = ss.str();

    Array t_t(im, jm, km, 0),  v_t(im, jm, km, 0), w_t(im, jm, km, 0), m_t(im, jm, km, 0), u_t(im, jm, km, 0);
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                m_t.x[ i ][ j ][ k ] = sqrt ( pow ( v.x[ i ][ j ][ k ] * u_0, 2 ) + pow ( w.x[ i ][ j ][ k ] * u_0, 2 ) );
                t_t.x[ i ][ j ][ k ] = t.x[ i ][ j ][ k ] * t_0 - t_0;
                v_t.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] * u_0;
                w_t.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] * u_0;
                u_t.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ] * u_0;
            }
        }
    }

    save_array(path + string("t") + postfix_str, t_t);
    save_array(path + string("v") + postfix_str, v_t);
    save_array(path + string("w") + postfix_str, w_t);
    save_array(path + string("u") + postfix_str, u_t);
    save_array(path + string("m") + postfix_str, m_t);
    save_array(path + string("h") + postfix_str, h);
    save_array(path + string("c") + postfix_str, c);
    save_array(path + string("p") + postfix_str, p_dyn);
}
