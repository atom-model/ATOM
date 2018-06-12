#include "cAtmosphereModel.h"

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

#include "BC_Atm.h"
#include "BC_Bath_Atm.h"
#include "BC_Thermo.h"
#include "Accuracy_Atm.h"
#include "RHS_Atm.h"
#include "RungeKutta_Atm.h"
#include "PostProcess_Atm.h"
#include "Pressure_Atm.h"
#include "Restore_Atm.h"
#include "Results_Atm.h"
#include "MinMax_Atm.h"

#include "Config.h"

using namespace std;
using namespace tinyxml2;

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

cAtmosphereModel::cAtmosphereModel() {
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
}

cAtmosphereModel::~cAtmosphereModel() {
    if(PythonStream::is_enable())
    {    
        std::cout.rdbuf(backup);
    }

    delete [] im_tropopause;
}
 
#include "cAtmosphereDefaults.cpp.inc"

void cAtmosphereModel::LoadConfig ( const char *filename ) 
{
    XMLDocument doc;
    XMLError err = doc.LoadFile ( filename );

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

}



void cAtmosphereModel::RunTimeSlice ( int Ma )
{
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

    //  Coordinate system in form of a spherical shell
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.Coordinates ( im, r0, dr );
    the.Coordinates ( jm, the0, dthe );
    phi.Coordinates ( km, phi0, dphi );


    //  initial values for the number of computed steps and the time
    double t_cretaceous = 0.;
    double co2_cretaceous = 0.;


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

    string bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    string bathymetry_filepath = bathymetry_path + "/" + bathymetry_name;

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
    BC_Bathymetry_Atmosphere LandArea ( NASATemperature, im, jm, km, co2_vegetation, co2_land, co2_ocean );

    //  topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
    LandArea.BC_MountainSurface ( bathymetry_filepath, L_atm, Topography, value_top, h, aux_w );

    //  class element for the computation of the ratio ocean to land areas, also supply and removal of CO2 on land, ocean and by vegetation
    LandArea.land_oceanFraction ( h );

    //  class calls for the solution of the flow properties

    //  class BC_Atmosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
    BC_Atmosphere  boundary ( im, jm, km, t_tropopause );

    //  class RHS_Atmosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
    RHS_Atmosphere  prepare ( im, jm, km, dt, dr, dthe, dphi, re, sc_WaterVapour, sc_CO2, g, pr, 
                              WaterVapour, Buoyancy, CO2, gam, sigma, lamda );
    RHS_Atmosphere  prepare_2D ( jm, km, dthe, dphi, re );

    //  class RungeKutta_Atmosphere for the explicit solution of the Navier-Stokes equations
    RungeKutta_Atmosphere  result ( im, jm, km, dt, dr, dphi, dthe );

    //  class Results_MSL_Atm to compute and show results on the mean sea level, MSL
    Results_MSL_Atm  calculate_MSL ( im, jm, km, sun, g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_equator, lv, ls, 
                                     cp_l, L_atm, dt, dr, dthe, dphi, r_air, R_Air, r_water, r_water_vapour, R_WaterVapour, 
                                     co2_vegetation, co2_ocean, co2_land, gam, t_pole, t_cretaceous, t_average );

    //  class Pressure for the subsequent computation of the pressure by a separate Euler equation
    Pressure_Atm  startPressure ( im, jm, km, dr, dthe, dphi );


    int Ma_prev;
    double t_cretaceous_prev;

    if( is_first_time_slice() ){
        Ma_prev = int(round(*get_current_time()));
        t_cretaceous_prev = 0.;
    }else{
        Ma_prev = int(round(*get_previous_time()));
    }



    //  class BC_Thermo for the initial and boundary conditions of the flow properties
    BC_Thermo  circulation ( output_path, im, jm, km, tropopause_equator, tropopause_pole, RadiationModel, NASATemperature, 
                             sun, declination, sun_position_lat, sun_position_lon, Ma, Ma_prev, Ma_max, Ma_max_half, dt, dr, 
                             dthe, dphi, g, ep, hp, u_0, p_0, t_0, c_0, sigma, lv, ls, cp_l, L_atm, r_air, R_Air, 
                             r_water_vapour, R_WaterVapour, co2_0, co2_cretaceous, co2_vegetation, co2_ocean, co2_land, 
                             co2_factor, c_tropopause, co2_tropopause, c_ocean, c_land, t_average, co2_average, co2_equator, 
                             co2_pole, t_cretaceous, t_cretaceous_prev, t_cretaceous_max, t_land, t_tropopause, t_equator, 
                             t_pole, gam, epsilon_equator, epsilon_pole, epsilon_tropopause, albedo_equator, albedo_pole, 
                             rad_equator, rad_pole );

    //  class Restore to restore the iterational values from new to old
    Restore_Atm oldnew( im, jm, km );

    //  class element calls for the preparation of initial conditions for the flow properties

    //  class element for the tropopause location as a parabolic distribution from pole to pole 
    circulation.TropopauseLocation ( im_tropopause );

    //  class element for the surface temperature from NASA for comparison
    //  if ( Ma == 0 ) circulation.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, temperature_NASA, t );
    circulation.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, temperature_NASA, t );

    //  class element for the surface temperature based on NASA temperature for progressing timeslices
    //  if ( ( Ma > 0 ) && ( NASATemperature == 1 ) ) circulation.BC_NASAbasedSurfTempRead ( Ma_prev, t_cretaceous_prev, t, c, cloud, ice );

    //  class element for the surface precipitation from NASA for comparison
    circulation.BC_Surface_Precipitation_NASA ( Name_SurfacePrecipitation_File, precipitation_NASA );

    //  class element for the parabolic temperature distribution from pol to pol, maximum temperature at equator
    circulation.BC_Temperature ( im_tropopause, t_cretaceous, t_cretaceous_prev, temperature_NASA, h, t, p_dyn, p_stat );
    //  t_cretaceous = circulation.out_t_cretaceous (  );

    //  class element for the correction of the temperature initial distribution around coasts
    if ( ( NASATemperature == 1 ) && ( Ma > 0 ) && !use_earthbyte_reconstruction) 
    {
        circulation.IC_Temperature_WestEastCoast ( h, t );
    }

    //  class element for the surface pressure computed by surface temperature with gas equation
    circulation.BC_Pressure ( p_stat, p_dyn, t, h );

    //  parabolic water vapour distribution from pol to pol, maximum water vapour volume at equator
    circulation.BC_WaterVapour ( im_tropopause, t_cretaceous, h, t, c );

    //  class element for the parabolic CO2 distribution from pol to pol, maximum CO2 volume at equator
    circulation.BC_CO2 ( im_tropopause, t_cretaceous, Vegetation, h, t, p_dyn, co2 );
    //  co2_cretaceous = circulation.out_co2 (  );

    // class element for the surface temperature computation by radiation flux density
    if ( RadiationModel == 1 ){
        circulation.BC_Radiation_multi_layer ( im_tropopause, t_cretaceous, CO2, albedo, epsilon, precipitable_water, 
                                               radiation_surface, Q_radiation, Q_latent, Q_sensible, Q_bottom, co2_total, 
                                               p_stat, t, c, h, epsilon_3D, radiation_3D, cloud, ice, co2 );
    }
    //  class element for the initial conditions for u-v-w-velocity components
    circulation.IC_CellStructure ( im_tropopause, h, u, v, w );

    // class element for the storing of velocity components, pressure and temperature for iteration start
    //  oldnew.restoreOldNew_3D(.9, u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, p_dynn, cn, cloudn, icen, co2n);
    //  oldnew.restoreOldNew_2D(.9, v, w, p_dyn, p_dynn, vn, wn);



    // ***********************************   start of pressure and velocity iterations ***********************************

    double emin = epsres * 100.;

    run_2D_loop(boundary, result, LandArea, prepare_2D, startPressure);
    
    cout << endl << endl;

    run_3D_loop( boundary, result, LandArea, prepare, startPressure, calculate_MSL, circulation);

    cout << endl << endl;

    n--;

    //write the ouput files
    write_file(bathymetry_name, output_path);

    if ( NASATemperature == 1 && !use_earthbyte_reconstruction ) 
    {
        t_cretaceous_prev = t_cretaceous;
        circulation.BC_NASAbasedSurfTempWrite ( Ma, t_cretaceous_prev, t, c, cloud, ice );
    }

    //  final remarks
    cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;
    if ( emin <= epsres ){
        cout << "***** steady solution reached! *****" << endl;
    }
}



void cAtmosphereModel::Run() 
{
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
}


void cAtmosphereModel::reset_arrays()
{
    // reset of arrays to the initial value
    // 1D arrays
    rad.initArray_1D(im, 1.); // radial coordinate direction
    the.initArray_1D(jm, 2.); // lateral coordinate direction
    phi.initArray_1D(km, 3.); // longitudinal coordinate direction

    // 2D arrays
    Topography.initArray_2D(jm, km, 0.); // topography
    value_top.initArray_2D(jm, km, 0.); // auxiliar topography

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
    min_max_3d.searchMinMax_3D ( " max 3D pressure dynamic ", " min 3D pressure dynamic ", "hPa", p_dyn, h, 0.01 );

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
    min_max_3d.searchMinMax_3D ( " max 3D rain ", " min 3D rain ", "g/kg", P_rain, h, 1000. );

    //  searching of maximum and minimum values of snow precipitation
    min_max_3d.searchMinMax_3D (  " max 3D snow ", " min 3D snow ", "g/kg", P_snow, h, 1000. );

    //  searching of maximum and minimum values of co2
    min_max_3d.searchMinMax_3D ( " max 3D co2 ", " min 3D co2 ", "ppm", co2, h, 280. );

    //  searching of maximum and minimum values of epsilon
    min_max_3d.searchMinMax_3D ( " max 3D epsilon ",  " min 3D epsilon ", "%", epsilon_3D, h );

    //  searching of maximum and minimum values of buoyancy force
    min_max_3d.searchMinMax_3D (  " max 3D buoyancy force ", " min 3D buoyancy force ", "N", BuoyancyForce, h );



    // 2D-fields

    cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" << endl << 
        " results based on two dimensional considerations of the problem" << endl;

    cout << endl << " co2 distribution row-wise: " << endl << endl;

    MinMax_Atm  min_max_2d ( jm, km );

    //  searching of maximum and minimum values of co2 total
    min_max_2d.searchMinMax_2D ( " max co2_total ", " min co2_total ", " ppm ", co2_total, h, 280. );

    cout << endl << " precipitation: " << endl << endl;

    //  searching of maximum and minimum values of precipitation
    min_max_2d.searchMinMax_2D (  " max precipitation ", " min precipitation ", "g/kg", Precipitation, h, 1.0/coeff_mmWS );
    max_Precipitation = min_max_2d.out_maxValue (  );

    //  searching of maximum and minimum values of precipitable water
    min_max_2d.searchMinMax_2D ( " max precipitable water ", " min precipitable water ", "g/kg", 
                                 precipitable_water, h, 1.0/coeff_mmWS  );

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


void cAtmosphereModel::write_file(std::string &bathymetry_name, std::string &output_path){
    int Ma = int(round(*get_current_time()));
    //  Printout:

    //  printout in ParaView files and sequel files

    //  class PostProcess_Atmosphaere for the printing of results
    PostProcess_Atmosphere write_File ( im, jm, km, output_path );

    //  writing of data in ParaView files
    //  radial data along constant hight above ground
    int i_radial = 0;
    //  int i_radial = 10;
    write_File.paraview_vtk_radial ( bathymetry_name, Ma, i_radial, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, 
                                     BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, radiation_3D, 
                                     Q_Latent, Q_Sensible, epsilon_3D, P_rain, P_snow, precipitable_water, Q_bottom, 
                                     Q_radiation, Q_latent, Q_sensible, Evaporation_Penman, Evaporation_Dalton, 
                                     Q_Evaporation, temperature_NASA, precipitation_NASA, Vegetation, albedo, epsilon, 
                                     Precipitation, Topography, temp_NASA );

    //  londitudinal data along constant latitudes
    int j_longal = 62;          // Mount Everest/Himalaya
    write_File.paraview_vtk_longal ( bathymetry_name, j_longal, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, 
                                     BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Q_Latent, 
                                     Q_Sensible, epsilon_3D, P_rain, P_snow );

    int k_zonal = 87;           // Mount Everest/Himalaya
    write_File.paraview_vtk_zonal ( bathymetry_name, k_zonal, n, hp, ep, R_Air, g, L_atm, u_0, t_0, p_0, r_air, c_0, co2_0, 
                                    h, p_dyn, p_stat, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, 
                                    Q_Latent, Q_Sensible, radiation_3D, epsilon_3D, P_rain, P_snow, S_v, S_c, S_i, S_r, 
                                    S_s, S_c_c );

    //  3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
    if(paraview_panorama_vts) //This function creates a large file. Use a flag to control if it is wanted.
    {
        write_File.paraview_panorama_vts ( bathymetry_name, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, t, p_dyn, p_stat, 
                                           BuoyancyForce, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Q_Latent, 
                                           Q_Sensible, epsilon_3D, P_rain, P_snow );
    }

    //  writing of v-w-data in the v_w_transfer file
    PostProcess_Atmosphere ppa ( im, jm, km, output_path );
    ppa.Atmosphere_v_w_Transfer ( bathymetry_name, u_0, v, w, t, p_dyn, Evaporation_Dalton, Precipitation );
    ppa.Atmosphere_PlotData ( bathymetry_name, u_0, t_0, h, v, w, t, c, Precipitation, precipitable_water );
}

void cAtmosphereModel::run_2D_loop( BC_Atmosphere &boundary, RungeKutta_Atmosphere &result,
                                    BC_Bathymetry_Atmosphere &LandArea, RHS_Atmosphere &prepare_2D, 
                                    Pressure_Atm &startPressure){
    int switch_2D = 0;    

    n = 1;
    int velocity_iter_2D = 0;
    int pressure_iter_2D = 0;

    double emin = epsres * 100.;
    int pressure_plus_2D = 1;

    double residuum;
    double residuum_old = 0.;

    int j_res = 0, k_res = 0;

    int Ma = int(round(*get_current_time())); 

    // ::::::::::: :::::::::::::::::::::::   begin of 2D loop for initial surface conditions: if ( switch_2D == 0 )   ::::
    if ( switch_2D != 1 )
    {
        // **************   iteration of initial conditions on the surface for the correction of flows close to coasts   **
        // **************   start of pressure and velocity iterations for the 2D iterational process   ********************
        // ::::::::::::::   begin of pressure loop_2D : if ( pressure_iter_2D > pressure_iter_max_2D )   ::::::::::::::::::
        for ( pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++)
        {
            // ::::::::   begin of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::::::
            for (velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++)
            {

                cout << endl << endl;
                cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
                cout << " 2D AGCM iterational process" << endl;
                cout << " max total iteration number nm = " << nm << endl << endl;

                cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, maximum "
                    << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl 
                    << endl << " Ma = " << Ma << "     n = " << n << "    velocity_iter_max_2D = " << velocity_iter_max_2D << 
                    "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " << pressure_iter_max_2D << 
                    "    pressure_iter_2D = " << pressure_iter_2D << endl;

                //  class BC_Atmosphaere for the geometry of a shell of a sphere
                boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
                boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

                //  old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
                Accuracy_Atm        min_Residuum_old_2D ( im, jm, km, dthe, dphi );
                min_Residuum_old_2D.residuumQuery_2D ( rad, the, v, w );
                emin = min_Residuum_old_2D.out_min (  );

                residuum_old = emin;

                //  class RungeKutta for the solution of the differential equations describing the flow properties
                result.solveRungeKutta_2D_Atmosphere ( prepare_2D, n, r_air, u_0, p_0, L_atm, rad, the, rhs_v, rhs_w, h, v, w, 
                                                       p_dyn, vn, wn, p_dynn, aux_v, aux_w );

                //  class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of 
                //  the continents and the ocean ground
                LandArea.BC_SolidGround ( RadiationModel, Ma, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, 
                                          t_equator, t_pole, 
                                          t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, 
                                          co2_cretaceous, pa, gam, sigma, h, u, v, w, t, p_dyn, c, cloud, ice, co2, 
                                          radiation_3D, Vegetation );

                //  new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
                Accuracy_Atm        min_Residuum_2D ( im, jm, km, dthe, dphi );
                min_Residuum_2D.residuumQuery_2D ( rad, the, v, w );
                emin = min_Residuum_2D.out_min (  );
                j_res = min_Residuum_2D.out_j_res (  );
                k_res = min_Residuum_2D.out_k_res (  );

                residuum = emin;
                emin = fabs ( ( residuum - residuum_old ) / residuum_old );

                //  state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
                Accuracy_Atm        min_Stationary_2D ( n, nm, Ma, im, jm, km, emin, j_res, k_res, velocity_iter_2D, 
                                                        pressure_iter_2D, velocity_iter_max_2D, pressure_iter_max_2D );
                min_Stationary_2D.steadyQuery_2D ( v, vn, w, wn, p_dyn, p_dynn );

                n++;
            }
            //  ::::::   end of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::::::::::::::


            //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
            if ( pressure_iter_2D == pressure_plus_2D )
            {
                startPressure.computePressure_2D ( r_air, rad, the, p_dyn, p_dynn, h, rhs_v, rhs_w, aux_v, aux_w );
                pressure_plus_2D = pressure_plus_2D + 1;
            }



            // limit of the computation in the sense of time steps
            if ( n > nm )
            {
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        }
        // :::::::::::::::::::   end of pressure loop_2D: if ( pressure_iter_2D > pressure_iter_max_2D )   ::::::::::
        n = 1;
        emin = epsres * 100.;
    }
    // ::::::::   end of 2D loop for initial surface conditions: if ( switch_2D == 0 )   :::::::::::::::::::::::::::::
}


void cAtmosphereModel::run_3D_loop( BC_Atmosphere &boundary, RungeKutta_Atmosphere &result,
                                    BC_Bathymetry_Atmosphere &LandArea, RHS_Atmosphere &prepare,
                                    Pressure_Atm &startPressure, Results_MSL_Atm &calculate_MSL,                  
                                    BC_Thermo &circulation){
    
    n = 1;

    double emin = epsres * 100.;
    int pressure_plus_3D = 1;

    double residuum;
    double residuum_old = 0.;

    int j_res = 0, k_res = 0;

    int Ma = int(round(*get_current_time()));
    int velocity_n = 1;

    int velocity_iter = 0;
    int pressure_iter = 0;

    // ::::::::::::::   begin of 3D pressure loop : if ( pressure_iter > pressure_iter_max )   ::::::::::::::::
    for ( pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++ )
    {
        // ::::::::::::   begin of 3D velocity loop : if ( velocity_iter > velocity_iter_max )   :::::::::::::::::::
        for ( velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++ )
        {
            //  query to realize zero divergence of the continuity equation ( div c = 0 )
            cout << endl << endl;
            cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            cout << " 3D AGCM iterational process" << endl;
            cout << " max total iteration number nm = " << nm << endl << endl;

            cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum "
                << "and current number of velocity iterations, maximum and current number of pressure iterations " << endl << 
                endl << " Ma = " << Ma << "     n = " << n << "    velocity_iter_max = " << velocity_iter_max << 
                "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << 
                "    pressure_iter = " << pressure_iter << endl;

            //  old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
            Accuracy_Atm        min_Residuum_old ( im, jm, km, dr, dthe, dphi );
            min_Residuum_old.residuumQuery_3D ( rad, the, u, v, w );
            emin = min_Residuum_old.out_min (  );

            residuum_old = emin;

            //  class BC_Atmosphaere for the geometry of a shell of a sphere
            boundary.BC_radius ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
            boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
            boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

            //Ice_Water_Saturation_Adjustment, distribution of cloud ice and cloud water dependent on water vapour amount and temperature
            if ( velocity_iter == velocity_n ){
                circulation.Ice_Water_Saturation_Adjustment ( im_tropopause, n, velocity_iter_max, RadiationModel, h, c, 
                                                              cn, cloud, cloudn, ice, icen, t, p_stat, S_c_c );
            }
            // class RungeKutta for the solution of the differential equations describing the flow properties
            result.solveRungeKutta_3D_Atmosphere ( prepare, n, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, 
                                                   r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, 
                                                   the, phi, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_cloud, rhs_ice, rhs_co2, 
                                                   h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, p_dynn, 
                                                   cn, cloudn, icen, co2n, aux_u, aux_v, aux_w, Q_Latent, BuoyancyForce, 
                                                   Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c, Topography, 
                                                   Evaporation_Dalton, Precipitation );

            //  class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of 
            //  the continents and the ocean ground
            LandArea.BC_SolidGround ( RadiationModel, Ma, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, 
                                      t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, 
                                      co2_tropopause, co2_cretaceous, pa, gam, sigma, h, u, v, w, t, p_dyn, c, cloud, 
                                      ice, co2, radiation_3D, Vegetation );

            // class element for the surface temperature computation by radiation flux density
            if ( RadiationModel == 1 ){
                 circulation.BC_Radiation_multi_layer ( im_tropopause, t_cretaceous, CO2, albedo, epsilon, precipitable_water, 
                                                        radiation_surface, Q_radiation, Q_latent, Q_sensible, Q_bottom, co2_total, 
                                                        p_stat, t, c, h, epsilon_3D, radiation_3D, cloud, ice, co2 );
            }
            //  new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
            Accuracy_Atm      min_Residuum ( im, jm, km, dr, dthe, dphi );
            min_Residuum.residuumQuery_3D ( rad, the, u, v, w );
            emin = min_Residuum.out_min (  );
            int i_res = min_Residuum.out_i_res (  );
            j_res = min_Residuum.out_j_res (  );
            k_res = min_Residuum.out_k_res (  );

            residuum = emin;
            emin = fabs ( ( residuum - residuum_old ) / residuum_old );

            //  statements on the convergence und iterational process
            Accuracy_Atm        min_Stationary ( n, nm, Ma, im, jm, km, emin, i_res, j_res, k_res, velocity_iter, pressure_iter, 
                                                 velocity_iter_max, pressure_iter_max, L_atm );
            min_Stationary.steadyQuery_3D ( u, un, v, vn, w, wn, t, tn, c, cn, cloud, cloudn, ice, icen, co2, co2n, p_dyn, p_dynn );

            // 3D_fields

            //  class element for the initial conditions the latent heat
            circulation.Latent_Heat ( rad, the, phi, h, t, tn, u, v, w, p_dyn, p_stat, c, ice, Q_Latent, Q_Sensible, 
                                      radiation_3D, Q_radiation, Q_latent, Q_sensible, Q_bottom );

            print_min_max_values();

            //  computation of vegetation areas
            LandArea.vegetationDistribution ( max_Precipitation, Precipitation, Vegetation, t, h );


            //  composition of results
            calculate_MSL.run_MSL_data ( n, velocity_iter_max, RadiationModel, t_cretaceous, rad, the, phi, h, c, cn, 
                                         co2, co2n, t, tn, p_dyn, p_stat, BuoyancyForce, u, v, w, Q_Latent, Q_Sensible, 
                                         radiation_3D, cloud, cloudn, ice, icen, P_rain, P_snow, aux_u, aux_v, aux_w, 
                                         temperature_NASA, precipitation_NASA, precipitable_water, Q_radiation, Q_Evaporation, 
                                         Q_latent, Q_sensible, Q_bottom, Evaporation_Penman, Evaporation_Dalton, Vegetation, 
                                         albedo, co2_total, Precipitation, S_v, S_c, S_i, S_r, S_s, S_c_c );



            //  Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation 
            //  distribution formed of rain and snow
            if ( velocity_iter == velocity_n )
            {
                circulation.Two_Category_Ice_Scheme ( n, velocity_iter_max, RadiationModel, t_cretaceous, h, c, t, p_stat, 
                                                      cloud, ice, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );
                velocity_n = velocity_iter + 2;
            }
            n++;
        }
        //  ::::::::::::   end of velocity loop_3D: if ( velocity_iter > velocity_iter_max )   ::::::::::::::::::::::::::::

        
        //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
        if ( pressure_iter == pressure_plus_3D )
        {
            startPressure.computePressure_3D ( r_air, rad, the, p_dyn, p_dynn, h, rhs_u, rhs_v, rhs_w, aux_u, aux_v, aux_w );
            pressure_plus_3D = pressure_plus_3D + 1;
        }

        //  Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation distribution 
        //  formed of rain and snow
        circulation.Two_Category_Ice_Scheme ( n, velocity_iter_max, RadiationModel, t_cretaceous, h, c, t, p_stat, cloud, ice, 
                                              P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );


        //  limit of the computation in the sense of time steps
        if ( n > nm )
        {
            cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
            break;
        }
    }
    //  :::::   end of pressure loop_3D: if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::

    if ( velocity_iter == velocity_iter_max ){
        cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of \
                maximum time steps ***** \n\n" << endl;
    }
}
