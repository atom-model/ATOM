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

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
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
#include "tinyxml2.h"
#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;

cAtmosphereModel::cAtmosphereModel() {
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.

    PythonStream::OverrideCout();

    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

    // set default configuration
    SetDefaultConfig();

}

cAtmosphereModel::~cAtmosphereModel() { }
 
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
    m_current_time = m_time_list.insert(float(Ma)).first;

    struct stat info;
    if( stat( output_path.c_str(), &info ) != 0 ){
        mkdir(output_path.c_str(), 0777);
    }
// maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km )
// maximum number of overall iterations ( n )
// maximum number of inner velocity loop iterations ( velocity_iter_max )
// maximum number of outer pressure loop iterations ( pressure_iter_max )

	const int im = 41, jm = 181, km = 361, nm = 200;

	int j_res = 0, k_res = 0;

	constexpr double pi180 = 180./ M_PI;								// pi180 = 57.3
	constexpr double the_degree = 1.;										// compares to 1° step size laterally
	constexpr double phi_degree = 1.;										// compares to 1° step size longitudinally

	double dthe = the_degree / pi180;										// dthe = the_degree / pi180 = 1.0 / 57.3 = 0.01745, 180 * .01745 = 3.141
	double dphi = phi_degree / pi180;										// dphi = phi_degree / pi180 = 1.0 / 57.3 = 0.01745, 360 * .01745 = 6.282
	double dr = 0.025;																	// 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
	double dt = 0.00001;																// time step coincides with the CFL condition

	const double the0 = 0.;														// North Pole
	const double phi0 = 0.;														// zero meridian in Greenwich
	const double r0 = 1.;															// earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km

	const double coeff_mmWS = r_air / r_water_vapour;			// coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]

	int *im_tropopause = new int [ jm ];									// location of the tropopause

//  class Array for 1-D, 2-D and 3-D field declarations

// 1D arrays
	Array_1D rad(im, 1.); // radial coordinate direction
	Array_1D the(jm, 2.); // lateral coordinate direction
	Array_1D phi(km, 3.); // longitudinal coordinate direction

// 2D arrays
	Array_2D Topography(jm, km, 0.); // topography
	Array_2D value_top(jm, km, 0.); // auxiliar topography

	Array_2D Vegetation(jm, km, 0.); // vegetation via precipitation

	Array_2D Precipitation(jm, km, 0.); // areas of higher precipitation
	Array_2D precipitable_water(jm, km, 0.); // areas of precipitable water in the air
	Array_2D precipitation_NASA(jm, km, 0.); // surface precipitation from NASA

	Array_2D radiation_surface(jm, km, 0.); // direct sun radiation, short wave

	Array_2D temperature_NASA(jm, km, 0.); // surface temperature from NASA
	Array_2D temp_NASA(jm, km, 0.); // surface temperature from NASA for print function

	Array_2D albedo(jm, km, 0.); // albedo = reflectivity
	Array_2D epsilon(jm, km, 0.); // epsilon = absorptivity

	Array_2D Q_radiation(jm, km, 0.); // heat from the radiation balance in [W/m2]
	Array_2D Q_Evaporation(jm, km, 0.); // evaporation heat of water by Kuttler
	Array_2D Q_latent(jm, km, 0.); // latent heat from bottom values by the energy transport equation
	Array_2D Q_sensible(jm, km, 0.); // sensible heat from bottom values by the energy transport equation
	Array_2D Q_bottom(jm, km, 0.); // difference by Q_radiation - Q_latent - Q_sensible

	Array_2D Evaporation_Dalton(jm, km, 0.); // evaporation by Dalton in [mm/d]
	Array_2D Evaporation_Penman(jm, km, 0.); // evaporation by Penman in [mm/d]

	Array_2D co2_total(jm, km, 0.); // areas of higher co2 concentration

// 3D arrays
	Array h(im, jm, km, 0.); // bathymetry, depth from sea level

	Array t(im, jm, km, ta); // temperature
	Array u(im, jm, km, ua); // u-component velocity component in r-direction
	Array v(im, jm, km, va); // v-component velocity component in theta-direction
	Array w(im, jm, km, wa); // w-component velocity component in phi-direction
	Array c(im, jm, km, ca); // water vapour
	Array cloud(im, jm, km, 0.); // cloud water
	Array ice(im, jm, km, 0.); // cloud ice
	Array co2(im, jm, km, coa); // CO2

	Array tn(im, jm, km, ta); // temperature new
	Array un(im, jm, km, ua); // u-velocity component in r-direction new
	Array vn(im, jm, km, va); // v-velocity component in theta-direction new
	Array wn(im, jm, km, wa); // w-velocity component in phi-direction new
	Array cn(im, jm, km, ca); // water vapour new
	Array cloudn(im, jm, km, 0.); // cloud water new
	Array icen(im, jm, km, 0.); // cloud ice new
	Array co2n(im, jm, km, coa); // CO2 new

	Array p_dyn(im, jm, km, pa); // dynamic pressure
	Array p_dynn(im, jm, km, pa); // dynamic pressure
	Array p_stat(im, jm, km, pa); // static pressure

	Array rhs_t(im, jm, km, 0.); // auxilliar field RHS temperature
	Array rhs_u(im, jm, km, 0.); // auxilliar field RHS u-velocity component
	Array rhs_v(im, jm, km, 0.); // auxilliar field RHS v-velocity component
	Array rhs_w(im, jm, km, 0.); // auxilliar field RHS w-velocity component
	Array rhs_c(im, jm, km, 0.); // auxilliar field RHS water vapour
	Array rhs_cloud(im, jm, km, 0.); // auxilliar field RHS cloud water
	Array rhs_ice(im, jm, km, 0.); // auxilliar field RHS cloud ice
	Array rhs_co2(im, jm, km, 0.); // auxilliar field RHS CO2

	Array aux_u(im, jm, km, 0.); // auxilliar field u-velocity component
	Array aux_v(im, jm, km, 0.); // auxilliar field v-velocity component
	Array aux_w(im, jm, km, 0.); // auxilliar field w-velocity component

	Array Q_Latent(im, jm, km, 0.); // latent heat
	Array Q_Sensible(im, jm, km, 0.); // sensible heat
	Array BuoyancyForce(im, jm, km, 0.); // buoyancy force, Boussinesque approximation
	Array epsilon_3D(im, jm, km, 0.); // emissivity/ absorptivity
	Array radiation_3D(im, jm, km, 0.); // radiation

	Array P_rain(im, jm, km, 0.); // rain precipitation mass rate
	Array P_snow(im, jm, km, 0.); // snow precipitation mass rate
	Array S_v(im, jm, km, 0.); // water vapour mass rate due to category two ice scheme
	Array S_c(im, jm, km, 0.); // cloud water mass rate due to category two ice scheme
	Array S_i(im, jm, km, 0.); // cloud ice mass rate due to category two ice scheme
	Array S_r(im, jm, km, 0.); // rain mass rate due to category two ice scheme
	Array S_s(im, jm, km, 0.); // snow mass rate due to category two ice scheme
	Array S_c_c(im, jm, km, 0.); // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique

//	cout << endl << "§§§§§§§§§§§§§§§§§§ after BC_Surface_Temperature_NASA §§§§§§§§§§§§§§§§§§§§§§§§" << endl;
//	cout << endl << " ***** printout of 3D-field temperature ***** " << endl << endl;
//	t.printArray( im, jm, km );

//	cout << endl << " ***** printout of 2D-field vegetation ***** " << endl << endl;
//	Vegetation.printArray_2D( jm, km );

//	cout << endl << " ***** printout of 1D-field radius ***** " << endl << endl;
//	rad.printArray_1D( im );

	cout.precision ( 6 );
	cout.setf ( ios::fixed );

//	Coordinate system in form of a spherical shell
//	rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
	rad.Coordinates ( im, r0, dr );
	the.Coordinates ( jm, the0, dthe );
	phi.Coordinates ( km, phi0, dphi );


//	initial values for the number of computed steps and the time
	int n = 1;
	int velocity_iter = 0;
	int velocity_iter_2D = 0;
	int pressure_iter = 0;
	int pressure_iter_2D = 0;
	int switch_2D = 0;
	double residuum;
	double residuum_old = 0.;
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


//	initialization of the bathymetry/topography

//	class BC_Bathymetry_Atmosphere for the geometrical boundary condition of the computational area
	BC_Bathymetry_Atmosphere		LandArea ( NASATemperature, im, jm, km, co2_vegetation, co2_land, co2_ocean );

//	topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
	LandArea.BC_MountainSurface ( bathymetry_filepath, L_atm, Topography, value_top, h, aux_w );

//	class element for the computation of the ratio ocean to land areas, also supply and removal of CO2 on land, ocean and by vegetation
	LandArea.land_oceanFraction ( h );


//	class calls for the solution of the flow properties

//	class BC_Atmosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
	BC_Atmosphere							boundary ( im, jm, km, t_tropopause );

//	class RHS_Atmosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
	RHS_Atmosphere						prepare ( im, jm, km, dt, dr, dthe, dphi, re, sc_WaterVapour, sc_CO2, g, pr, WaterVapour, Buoyancy, CO2, gam, sigma, lamda );
	RHS_Atmosphere						prepare_2D ( jm, km, dthe, dphi, re );

//	class RungeKutta_Atmosphere for the explicit solution of the Navier-Stokes equations
	RungeKutta_Atmosphere			result ( im, jm, km, dt, dr, dphi, dthe );

//	class Results_MSL_Atm to compute and show results on the mean sea level, MSL
	Results_MSL_Atm						calculate_MSL ( im, jm, km, sun, g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_equator, lv, ls, cp_l, L_atm, dt, dr, dthe, dphi, r_air, R_Air, r_water, r_water_vapour, R_WaterVapour, co2_vegetation, co2_ocean, co2_land, gam, t_pole, t_cretaceous, t_average );

//	class Pressure for the subsequent computation of the pressure by a separate Euler equation
	Pressure_Atm										startPressure ( im, jm, km, dr, dthe, dphi );


	int Ma_prev;
	double t_cretaceous_prev;

	if( is_first_time_slice() ){
		Ma_prev = int(round(*get_current_time()));
		t_cretaceous_prev = 0.;
	}else{
		Ma_prev = int(round(*get_previous_time()));
	}



//	class BC_Thermo for the initial and boundary conditions of the flow properties
	BC_Thermo									circulation ( output_path, im, jm, km, tropopause_equator, tropopause_pole, RadiationModel, NASATemperature, sun, declination, sun_position_lat, sun_position_lon, Ma, Ma_prev, Ma_max, Ma_max_half, dt, dr, dthe, dphi, g, ep, hp, u_0, p_0, t_0, c_0, sigma, lv, ls, cp_l, L_atm, r_air, R_Air, r_water_vapour, R_WaterVapour, co2_0, co2_cretaceous, co2_vegetation, co2_ocean, co2_land, co2_factor, c_tropopause, co2_tropopause, c_ocean, c_land, t_average, co2_average, co2_equator, co2_pole, t_cretaceous, t_cretaceous_prev, t_cretaceous_max, t_land, t_tropopause, t_equator, t_pole, gam, epsilon_equator, epsilon_pole, epsilon_tropopause, albedo_equator, albedo_pole, rad_equator, rad_pole );

//	class Restore to restore the iterational values from new to old
	Restore_Atm								oldnew( im, jm, km );


//	class element calls for the preparation of initial conditions for the flow properties

//	class element for the tropopause location as a parabolic distribution from pole to pole 
	circulation.TropopauseLocation ( im_tropopause );

//  class element for the surface temperature from NASA for comparison
//	if ( Ma == 0 ) circulation.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, temperature_NASA, t );
	circulation.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, temperature_NASA, t );

//  class element for the surface temperature based on NASA temperature for progressing timeslices
//	if ( ( Ma > 0 ) && ( NASATemperature == 1 ) ) circulation.BC_NASAbasedSurfTempRead ( Ma_prev, t_cretaceous_prev, t, c, cloud, ice );

//  class element for the surface precipitation from NASA for comparison
	circulation.BC_Surface_Precipitation_NASA ( Name_SurfacePrecipitation_File, precipitation_NASA );

//  class element for the parabolic temperature distribution from pol to pol, maximum temperature at equator
	circulation.BC_Temperature ( im_tropopause, t_cretaceous, t_cretaceous_prev, temperature_NASA, h, t, p_dyn, p_stat );
//	t_cretaceous = circulation.out_t_cretaceous (  );

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
//	co2_cretaceous = circulation.out_co2 (  );

// class element for the surface temperature computation by radiation flux density
	if ( RadiationModel == 1 ) circulation.BC_Radiation_multi_layer ( im_tropopause, t_cretaceous, n, CO2, albedo, epsilon, precipitable_water, radiation_surface, Q_radiation, Q_latent, Q_sensible, Q_bottom, co2_total, p_stat, t, c, h, epsilon_3D, radiation_3D, cloud, ice, co2 );

// 	class element for the initial conditions for u-v-w-velocity components
	circulation.IC_CellStructure ( im_tropopause, h, u, v, w );

// class element for the storing of velocity components, pressure and temperature for iteration start
//	oldnew.restoreOldNew_3D(.9, u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, p_dynn, cn, cloudn, icen, co2n);
//	oldnew.restoreOldNew_2D(.9, v, w, p_dyn, p_dynn, vn, wn);



// ***********************************   start of pressure and velocity iterations ************************************************************************

	double emin = epsres * 100.;
	int pressure_plus_2D = 1;
	int pressure_plus_3D = 1;

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of 2D loop for initial surface conditions: if ( switch_2D == 0 )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	if ( switch_2D != 1 )
	{
// **********************************   iteration of initial conditions on the surface for the correction of flows close to coasts   *********************************
// **********************************   start of pressure and velocity iterations for the 2D iterational process   *********************************
// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of pressure loop_2D : if ( pressure_iter_2D > pressure_iter_max_2D )   :::::::::::::::::::::::::::::::::::::::::::
		for ( pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++)
		{
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
			for (velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++)
			{

				cout << endl << endl;
				cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
				cout << " 2D AGCM iterational process" << endl;
				cout << " max total iteration number nm = " << nm << endl << endl;

				cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations " << endl << endl << " Ma = " << Ma << "     n = " << n << "    velocity_iter_max_2D = " << velocity_iter_max_2D << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " << pressure_iter_max_2D << "    pressure_iter_2D = " << pressure_iter_2D << endl;

//	class BC_Atmosphaere for the geometry of a shell of a sphere
				boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
				boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

//	old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
				Accuracy_Atm		min_Residuum_old_2D ( im, jm, km, dthe, dphi );
				min_Residuum_old_2D.residuumQuery_2D ( rad, the, v, w );
				emin = min_Residuum_old_2D.out_min (  );

				residuum_old = emin;

//	class RungeKutta for the solution of the differential equations describing the flow properties
				result.solveRungeKutta_2D_Atmosphere ( prepare_2D, n, r_air, u_0, p_0, L_atm, rad, the, rhs_v, rhs_w, h, v, w, p_dyn, vn, wn, p_dynn, aux_v, aux_w );

//	class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
				LandArea.BC_SolidGround ( RadiationModel, Ma, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, co2_cretaceous, pa, gam, sigma, h, u, v, w, t, p_dyn, c, cloud, ice, co2, radiation_3D, Vegetation );

//	new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
				Accuracy_Atm		min_Residuum_2D ( im, jm, km, dthe, dphi );
				min_Residuum_2D.residuumQuery_2D ( rad, the, v, w );
				emin = min_Residuum_2D.out_min (  );
				j_res = min_Residuum_2D.out_j_res (  );
				k_res = min_Residuum_2D.out_k_res (  );

				residuum = emin;
				emin = fabs ( ( residuum - residuum_old ) / residuum_old );

//	state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
				Accuracy_Atm		min_Stationary_2D ( n, nm, Ma, im, jm, km, emin, j_res, k_res, velocity_iter_2D, pressure_iter_2D, velocity_iter_max_2D, pressure_iter_max_2D );
				min_Stationary_2D.steadyQuery_2D ( v, vn, w, wn, p_dyn, p_dynn );

				n++;
			}

//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



//  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
		if ( pressure_iter_2D == pressure_plus_2D )
		{
			startPressure.computePressure_2D ( r_air, rad, the, p_dyn, p_dynn, h, rhs_v, rhs_w, aux_v, aux_w );
			pressure_plus_2D = pressure_plus_2D + 1;
		}



//		limit of the computation in the sense of time steps
			if ( n > nm )
			{
				cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
				break;
			}
		}
		n = 1;
		switch_2D = 1;						// 2D calculations finished
		emin = epsres * 100.;
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of pressure loop_2D: if ( pressure_iter_2D > pressure_iter_max_2D )   :::::::::::::::::::::::::::::::::::::::::::
	}

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of 2D loop for initial surface conditions: if ( switch_2D == 0 )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	cout << endl << endl;


	int velocity_n = 1;

//	goto Printout;

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of 3D pressure loop : if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::
	for ( pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++ )
	{
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of 3D velocity loop : if ( velocity_iter > velocity_iter_max )   :::::::::::::::::::::::::::::::::::::::::::
		for ( velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++ )
		{
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of zero-divergence loop: while ( min >= epsres )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//	query to realize zero divergence of the continuity equation ( div c = 0 )
			cout << endl << endl;
			cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
			cout << " 3D AGCM iterational process" << endl;
			cout << " max total iteration number nm = " << nm << endl << endl;

			cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations " << endl << endl << " Ma = " << Ma << "     n = " << n << "    velocity_iter_max = " << velocity_iter_max << "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << "    pressure_iter = " << pressure_iter << endl;

//	old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy_Atm		min_Residuum_old ( im, jm, km, dr, dthe, dphi );
			min_Residuum_old.residuumQuery_3D ( rad, the, u, v, w );
			emin = min_Residuum_old.out_min (  );

			residuum_old = emin;

//	class BC_Atmosphaere for the geometry of a shell of a sphere
			boundary.BC_radius ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
			boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
			boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

//		Ice_Water_Saturation_Adjustment, distribution of cloud ice and cloud water dependent on water vapour amount and temperature
			if ( velocity_iter == velocity_n )		circulation.Ice_Water_Saturation_Adjustment ( im_tropopause, n, velocity_iter_max, RadiationModel, h, c, cn, cloud, cloudn, ice, icen, t, p_stat, S_c_c );

// 		class RungeKutta for the solution of the differential equations describing the flow properties
			result.solveRungeKutta_3D_Atmosphere ( prepare, n, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_cloud, rhs_ice, rhs_co2, h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, p_dynn, cn, cloudn, icen, co2n, aux_u, aux_v, aux_w, Q_Latent, BuoyancyForce, Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c, Topography, Evaporation_Dalton, Precipitation );

//	class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
			LandArea.BC_SolidGround ( RadiationModel, Ma, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, co2_cretaceous, pa, gam, sigma, h, u, v, w, t, p_dyn, c, cloud, ice, co2, radiation_3D, Vegetation );

// class element for the surface temperature computation by radiation flux density
			if ( RadiationModel == 1 )			circulation.BC_Radiation_multi_layer ( im_tropopause, t_cretaceous, n, CO2, albedo, epsilon, precipitable_water, radiation_surface, Q_radiation, Q_latent, Q_sensible, Q_bottom, co2_total, p_stat, t, c, h, epsilon_3D, radiation_3D, cloud, ice, co2 );

//	new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy_Atm	  min_Residuum ( im, jm, km, dr, dthe, dphi );
			min_Residuum.residuumQuery_3D ( rad, the, u, v, w );
			emin = min_Residuum.out_min (  );
			int i_res = min_Residuum.out_i_res (  );
			j_res = min_Residuum.out_j_res (  );
			k_res = min_Residuum.out_k_res (  );

			residuum = emin;
			emin = fabs ( ( residuum - residuum_old ) / residuum_old );

//	statements on the convergence und iterational process
			Accuracy_Atm		min_Stationary ( n, nm, Ma, im, jm, km, emin, i_res, j_res, k_res, velocity_iter, pressure_iter, velocity_iter_max, pressure_iter_max, L_atm );
			min_Stationary.steadyQuery_3D ( u, un, v, vn, w, wn, t, tn, c, cn, cloud, cloudn, ice, icen, co2, co2n, p_dyn, p_dynn );

// 3D_fields

//  class element for the initial conditions the latent heat
			circulation.Latent_Heat ( rad, the, phi, h, t, tn, u, v, w, p_dyn, p_stat, c, ice, Q_Latent, Q_Sensible, radiation_3D, Q_radiation, Q_latent, Q_sensible, Q_bottom );

//	searching of maximum and minimum values of temperature
			string str_max_temperature = " max 3D temperature ", str_min_temperature = " min 3D temperature ", str_unit_temperature = "C";
			MinMax_Atm		minmaxTemperature ( im, jm, km, u_0 );
			minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

//	searching of maximum and minimum values of u-component
			string str_max_u = " max 3D u-component ", str_min_u = " min 3D u-component ", str_unit_u = "m/s";
			MinMax_Atm		minmax_u ( im, jm, km, u_0 );
			minmax_u.searchMinMax_3D ( str_max_u, str_min_u, str_unit_u, u, h );

//	searching of maximum and minimum values of v-component
			string str_max_v = " max 3D v-component ", str_min_v = " min 3D v-component ", str_unit_v = "m/s";
			MinMax_Atm		minmax_v ( im, jm, km, u_0 );
			minmax_v.searchMinMax_3D ( str_max_v, str_min_v, str_unit_v, v, h );

//	searching of maximum and minimum values of w-component
			string str_max_w = " max 3D w-component ", str_min_w = " min 3D w-component ", str_unit_w = "m/s";
			MinMax_Atm		minmax_w ( im, jm, km, u_0 );
			minmax_w.searchMinMax_3D ( str_max_w, str_min_w, str_unit_w, w, h );

//	searching of maximum and minimum values of dynamic pressure
			string str_max_pressure = " max 3D pressure dynamic ", str_min_pressure = " min 3D pressure dynamic ", str_unit_pressure = "hPa";
			MinMax_Atm	minmaxPressure ( im, jm, km, u_0 );
			minmaxPressure.searchMinMax_3D ( str_max_pressure, str_min_pressure, str_unit_pressure, p_dyn, h );

//	searching of maximum and minimum values of static pressure
			string str_max_pressure_stat = " max 3D pressure static ", str_min_pressure_stat = " min 3D pressure static ", str_unit_pressure_stat = "hPa";
			MinMax_Atm		minmaxPressure_stat ( im, jm, km, u_0 );
			minmaxPressure_stat.searchMinMax_3D ( str_max_pressure_stat, str_min_pressure_stat, str_unit_pressure_stat, p_stat, h );

			cout << endl << " energies in the three dimensional space: " << endl << endl;

//	searching of maximum and minimum values of radiation_3D
			string str_max_radiation_3D = " max 3D radiation ", str_min_radiation_3D = " min 3D radiation ", str_unit_radiation_3D = "W/m2";
			MinMax_Atm		minmaxRadiation ( im, jm, km, u_0 );
			minmaxRadiation.searchMinMax_3D ( str_max_radiation_3D, str_min_radiation_3D, str_unit_radiation_3D, radiation_3D, h );

//	searching of maximum and minimum values of sensible heat
			string str_max_Q_Sensible = " max 3D sensible heat ", str_min_Q_Sensible = " min 3D sensible heat ", str_unit_Q_Sensible = "W/m2";
			MinMax_Atm	minmaxQ_Sensible ( im, jm, km, u_0 );
			minmaxQ_Sensible.searchMinMax_3D ( str_max_Q_Sensible, str_min_Q_Sensible, str_unit_Q_Sensible, Q_Sensible, h );

//	searching of maximum and minimum values of latency
			string str_max_latency = " max 3D latent heat ", str_min_latency = " min 3D latent heat ", str_unit_latency = "W/m2";
			MinMax_Atm		minmaxQ_Latent ( im, jm, km, u_0 );
			minmaxQ_Latent.searchMinMax_3D ( str_max_latency, str_min_latency, str_unit_latency, Q_Latent, h );

			cout << endl << " greenhouse gases: " << endl << endl;

//	searching of maximum and minimum values of water vapour
			string str_max_water_vapour = " max 3D water vapour ", str_min_water_vapour = " min 3D water vapour ", str_unit_water_vapour = "g/kg";
			MinMax_Atm		minmaxWaterVapour ( im, jm, km, u_0 );
			minmaxWaterVapour.searchMinMax_3D ( str_max_water_vapour, str_min_water_vapour, str_unit_water_vapour, c, h );

//	searching of maximum and minimum values of cloud water
			string str_max_cloud_water = " max 3D cloud water ", str_min_cloud_water = " min 3D cloud water ", str_unit_cloud_water = "g/kg";
			MinMax_Atm		minmaxCloudWater ( im, jm, km, u_0 );
			minmaxCloudWater.searchMinMax_3D ( str_max_cloud_water, str_min_cloud_water, str_unit_cloud_water, cloud, h );

//	searching of maximum and minimum values of cloud ice
			string str_max_cloud_ice = " max 3D cloud ice ", str_min_cloud_ice = " min 3D cloud ice ", str_unit_cloud_ice = "g/kg";
			MinMax_Atm		minmaxCloudIce ( im, jm, km, u_0 );
			minmaxCloudIce.searchMinMax_3D ( str_max_cloud_ice, str_min_cloud_ice, str_unit_cloud_ice, ice, h );

//	searching of maximum and minimum values of rain precipitation
			string str_max_P_rain = " max 3D rain ", str_min_P_rain = " min 3D rain ", str_unit_P_rain = "g/kg";
			MinMax_Atm		minmaxPRain ( im, jm, km, u_0 );
			minmaxPRain.searchMinMax_3D ( str_max_P_rain, str_min_P_rain, str_unit_P_rain, P_rain, h );

//	searching of maximum and minimum values of snow precipitation
			string str_max_P_snow = " max 3D snow ", str_min_P_snow = " min 3D snow ", str_unit_P_snow = "g/kg";
			MinMax_Atm		minmaxPSnow ( im, jm, km, u_0 );
			minmaxPSnow.searchMinMax_3D ( str_max_P_snow, str_min_P_snow, str_unit_P_snow, P_snow, h );

//	searching of maximum and minimum values of co2
			string str_max_co2 = " max 3D co2 ", str_min_co2 = " min 3D co2 ", str_unit_co2 = "ppm";
			MinMax_Atm		minmaxCO2 ( im, jm, km, u_0 );
			minmaxCO2.searchMinMax_3D ( str_max_co2, str_min_co2, str_unit_co2, co2, h );

//	searching of maximum and minimum values of epsilon
			string str_max_epsilon = " max 3D epsilon ", str_min_epsilon = " min 3D epsilon ", str_unit_epsilon = "%";
			MinMax_Atm		minmaxEpsilon_3D ( im, jm, km, u_0 );
			minmaxEpsilon_3D.searchMinMax_3D ( str_max_epsilon, str_min_epsilon, str_unit_epsilon, epsilon_3D, h );

//	searching of maximum and minimum values of buoyancy force
			string str_max_buoyancy_force = " max 3D buoyancy force ", str_min_buoyancy_force = " min 3D buoyancy force ", str_unit_buoyancy_force = "N";
			MinMax_Atm		minmaxBuoyancyForce ( im, jm, km, u_0 );
			minmaxBuoyancyForce.searchMinMax_3D ( str_max_buoyancy_force, str_min_buoyancy_force, str_unit_buoyancy_force, BuoyancyForce, h );



// 2D-fields

//	searching of maximum and minimum values of co2 total
			cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" << endl << " results based on two dimensional considerations of the problem" << endl;

			cout << endl << " co2 distribution row-wise: " << endl << endl;

			string str_max_co2_total = " max co2_total ", str_min_co2_total = " min co2_total ", str_unit_co2_total = " / ";
			MinMax_Atm	minmaxCO2_total ( jm, km, coeff_mmWS );
			minmaxCO2_total.searchMinMax_2D ( str_max_co2_total, str_min_co2_total, str_unit_co2_total, co2_total, h );
//			max_CO2_total = minmaxCO2_total.out_maxValue (  );

			cout << endl << " precipitation: " << endl << endl;

//	searching of maximum and minimum values of precipitation
			string str_max_precipitation = " max precipitation ", str_min_precipitation = " min precipitation ", str_unit_precipitation = "mm";
			MinMax_Atm		minmaxPrecipitation ( jm, km, coeff_mmWS );
			minmaxPrecipitation.searchMinMax_2D ( str_max_precipitation, str_min_precipitation, str_unit_precipitation, Precipitation, h );
			double max_Precipitation = minmaxPrecipitation.out_maxValue (  );

//	searching of maximum and minimum values of precipitable water
			string str_max_precipitable_water = " max precipitable water ", str_min_precipitable_water = " min precipitable water ", str_unit_precipitable_water = "mm";
			MinMax_Atm		minmaxPrecipitable_water ( jm, km, coeff_mmWS );
			minmaxPrecipitable_water.searchMinMax_2D ( str_max_precipitable_water, str_min_precipitable_water, str_unit_precipitable_water, precipitable_water, h );

			cout << endl << " energies at see level without convection influence: " << endl << endl;

//	searching of maximum and minimum values of radiation
			string str_max_Q_radiation = " max 2D Q radiation ", str_min_Q_radiation = " min 2D Q radiation ", str_unit_Q_radiation = "W/m2";
			MinMax_Atm		minmaxQ_radiation ( jm, km, coeff_mmWS );
			minmaxQ_radiation.searchMinMax_2D ( str_max_Q_radiation, str_min_Q_radiation, str_unit_Q_radiation, Q_radiation, h );

//	searching of maximum and minimum values of latent energy
			string str_max_Q_latent = " max 2D Q latent ", str_min_Q_latent = " min 2D Q latent ", str_unit_Q_latent = "W/m2";
			MinMax_Atm		minmaxQ_latent ( jm, km, coeff_mmWS );
			minmaxQ_latent.searchMinMax_2D ( str_max_Q_latent, str_min_Q_latent, str_unit_Q_latent, Q_latent, h );

//	searching of maximum and minimum values of sensible energy
			string str_max_Q_sensible = " max 2D Q sensible ", str_min_Q_sensible = " min 2D Q sensible ", str_unit_Q_sensible = "W/m2";
			MinMax_Atm	minmaxQ_sensible ( jm, km, coeff_mmWS );
			minmaxQ_sensible.searchMinMax_2D ( str_max_Q_sensible, str_min_Q_sensible, str_unit_Q_sensible, Q_sensible, h );

//	searching of maximum and minimum values of bottom heat

			string str_max_Q_bottom = " max 2D Q bottom ", str_min_Q_bottom = " min 2D Q bottom heat ", str_unit_Q_bottom = "W/m2";
			MinMax_Atm		minmaxQ_bottom ( jm, km, coeff_mmWS );
			minmaxQ_bottom.searchMinMax_2D ( str_max_Q_bottom, str_min_Q_bottom, str_unit_Q_bottom, Q_bottom, h );


			cout << endl << " secondary data: " << endl << endl;


//	searching of maximum and minimum values of Evaporation
			string str_max_heat_t_Evaporation = " max heat Evaporation ", str_min_heat_t_Evaporation = " min heat Evaporation ", str_unit_heat_t_Evaporation = " kJ/kg";
			MinMax_Atm		minmaxQ_t_Evaporation ( jm, km, coeff_mmWS );
			minmaxQ_t_Evaporation.searchMinMax_2D ( str_max_heat_t_Evaporation, str_min_heat_t_Evaporation, str_unit_heat_t_Evaporation, Q_Evaporation, h );


//	searching of maximum and minimum values of Evaporation by Dalton
			string str_max_t_Evaporation_Dalton = " max Evaporation Dalton ", str_min_t_Evaporation_Dalton = " min Evaporation Dalton ", str_unit_t_Evaporation_Dalton = "mm/d";
			MinMax_Atm		minmaxt_Evaporation_Dalton ( jm, km, coeff_mmWS );
			minmaxt_Evaporation_Dalton.searchMinMax_2D ( str_max_t_Evaporation_Dalton, str_min_t_Evaporation_Dalton, str_unit_t_Evaporation_Dalton, Evaporation_Dalton, h );

//	searching of maximum and minimum values of Evaporation by Penman
			string str_max_t_Evaporation_Penman = " max Evaporation Penman ", str_min_t_Evaporation_Penman = " min Evaporation Penman ", str_unit_t_Evaporation_Penman = "mm/d";
			MinMax_Atm		minmaxt_Evaporation_Penman ( jm, km, coeff_mmWS );
			minmaxt_Evaporation_Penman.searchMinMax_2D ( str_max_t_Evaporation_Penman, str_min_t_Evaporation_Penman, str_unit_t_Evaporation_Penman, Evaporation_Penman, h );

			cout << endl << " properties of the atmosphere at the surface: " << endl << endl;

//	searching of maximum and minimum values of albedo
			string str_max_albedo = " max 2D albedo ", str_min_albedo = " min 2D albedo ", str_unit_albedo = "%";
			MinMax_Atm		minmaxAlbedo ( jm, km, coeff_mmWS );
			minmaxAlbedo.searchMinMax_2D ( str_max_albedo, str_min_albedo, str_unit_albedo, albedo, h );

//	searching of maximum and minimum values of epsilon
			string str_max_epsilon_2D = " max 2D epsilon ", str_min_epsilon_2D = " min 2D epsilon ", str_unit_epsilon_2D = "%";
			MinMax_Atm	minmaxEpsilon ( jm, km, coeff_mmWS );
			minmaxEpsilon.searchMinMax_2D ( str_max_epsilon_2D, str_min_epsilon_2D, str_unit_epsilon_2D, epsilon, h );

//	searching of maximum and minimum values of topography
			string str_max_topography = " max 2D topography ", str_min_topography = " min 2D topography ", str_unit_topography = "m";
			MinMax_Atm	minmaxTopography ( jm, km, coeff_mmWS );
			minmaxTopography.searchMinMax_2D ( str_max_topography, str_min_topography, str_unit_topography, Topography, h );

//	computation of vegetation areas
			LandArea.vegetationDistribution ( max_Precipitation, Precipitation, Vegetation, t, h );


//	composition of results
			calculate_MSL.run_MSL_data ( n, velocity_iter_max, RadiationModel, t_cretaceous, rad, the, phi, h, c, cn, co2, co2n, t, tn, p_dyn, p_stat, BuoyancyForce, u, v, w, Q_Latent, Q_Sensible, radiation_3D, cloud, cloudn, ice, icen, P_rain, P_snow, aux_u, aux_v, aux_w, temperature_NASA, precipitation_NASA, precipitable_water, Q_radiation, Q_Evaporation, Q_latent, Q_sensible, Q_bottom, Evaporation_Penman, Evaporation_Dalton, Vegetation, albedo, co2_total, Precipitation, S_v, S_c, S_i, S_r, S_s, S_c_c );

//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//	Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation distribution formed of rain and snow
			if ( velocity_iter == velocity_n )
			{
				circulation.Two_Category_Ice_Scheme ( n, velocity_iter_max, RadiationModel, t_cretaceous, h, c, t, p_stat, cloud, ice, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );
				velocity_n = velocity_iter + 2;
			}
			n++;
		}
//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of velocity loop_3D: if ( velocity_iter > velocity_iter_max )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
		if ( pressure_iter == pressure_plus_3D )
		{
			startPressure.computePressure_3D ( r_air, rad, the, p_dyn, p_dynn, h, rhs_u, rhs_v, rhs_w, aux_u, aux_v, aux_w );
			pressure_plus_3D = pressure_plus_3D + 1;
		}

//	Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation distribution formed of rain and snow
	circulation.Two_Category_Ice_Scheme ( n, velocity_iter_max, RadiationModel, t_cretaceous, h, c, t, p_stat, cloud, ice, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );

//	goto Printout;


//	limit of the computation in the sense of time steps
		if ( n > nm )
		{
			cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
			break;
		}
	}
//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of pressure loop_3D: if ( pressure_iter > pressure_iter_max )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	cout << endl << endl;

	n--;

//	Printout:

//	printout in ParaView files and sequel files

//	class PostProcess_Atmosphaere for the printing of results
	PostProcess_Atmosphere write_File ( im, jm, km, output_path );

//	writing of data in ParaView files
//	radial data along constant hight above ground
	int i_radial = 0;
//	int i_radial = 10;
	write_File.paraview_vtk_radial ( bathymetry_name, Ma, i_radial, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, radiation_3D, Q_Latent, Q_Sensible, epsilon_3D, P_rain, P_snow, precipitable_water, Q_bottom, Q_radiation, Q_latent, Q_sensible, Evaporation_Penman, Evaporation_Dalton, Q_Evaporation, temperature_NASA, precipitation_NASA, Vegetation, albedo, epsilon, Precipitation, Topography, temp_NASA );

//	londitudinal data along constant latitudes
	int j_longal = 62;			// Mount Everest/Himalaya
	write_File.paraview_vtk_longal ( bathymetry_name, j_longal, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Q_Latent, Q_Sensible, epsilon_3D, P_rain, P_snow );

	int k_zonal = 87;			// Mount Everest/Himalaya
	write_File.paraview_vtk_zonal ( bathymetry_name, k_zonal, n, hp, ep, R_Air, g, L_atm, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Q_Latent, Q_Sensible, radiation_3D, epsilon_3D, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );

//	3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
	if(paraview_panorama_vts) //This function creates a large file. Use a flag to control if it is wanted.
    {
        write_File.paraview_panorama_vts ( bathymetry_name, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, t, p_dyn, p_stat, BuoyancyForce, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Q_Latent, Q_Sensible, epsilon_3D, P_rain, P_snow );
    }

//	writing of v-w-data in the v_w_transfer file
	PostProcess_Atmosphere ppa ( im, jm, km, output_path );
	ppa.Atmosphere_v_w_Transfer ( bathymetry_name, u_0, v, w, t, p_dyn, Evaporation_Dalton, Precipitation );
	ppa.Atmosphere_PlotData ( bathymetry_name, u_0, t_0, h, v, w, t, c, Precipitation, precipitable_water );

	t_cretaceous_prev = t_cretaceous;

	if ( NASATemperature == 1 && !use_earthbyte_reconstruction ) 
    {
        circulation.BC_NASAbasedSurfTempWrite ( Ma, t_cretaceous_prev, t, c, cloud, ice );
    }

// reset of results to the initial value
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Vegetation.y[ j ][ k ] = 0.;

			precipitable_water.y[ j ][ k ] = 0.;
			Precipitation.y[ j ][ k ] = 0.;
			precipitation_NASA.y[ j ][ k ] = 0.;

			temperature_NASA.y[ j ][ k ] = 0.;

			radiation_surface.y[ j ][ k ] = 0.;

			albedo.y[ j ][ k ] = 0.;
			epsilon.y[ j ][ k ] = 0.;

			Q_radiation.y[ j ][ k ] = 0.;
			Q_Evaporation.y[ j ][ k ] = 0.;
			Q_latent.y[ j ][ k ] = 0.;
			Q_sensible.y[ j ][ k ] = 0.;
			Q_bottom.y[ j ][ k ] = 0.;

			Evaporation_Dalton.y[ j ][ k ] = 0.;
			Evaporation_Penman.y[ j ][ k ] = 0.;

			co2_total.y[ j ][ k ] = 0.;


			for ( int i = 0; i < im; i++ )
			{
				h.x[ i ][ j ][ k ] = 0.;
				t.x[ i ][ j ][ k ] = 0.;
				u.x[ i ][ j ][ k ] = 0.;
				v.x[ i ][ j ][ k ] = 0.;
				w.x[ i ][ j ][ k ] = 0.;
				c.x[ i ][ j ][ k ] = 0.;
				cloud.x[ i ][ j ][ k ] = 0.;
				ice.x[ i ][ j ][ k ] = 0.;
				co2.x[ i ][ j ][ k ] = 0.;

				tn.x[ i ][ j ][ k ] = 0.;
				un.x[ i ][ j ][ k ] = 0.;
				vn.x[ i ][ j ][ k ] = 0.;
				wn.x[ i ][ j ][ k ] = 0.;
				cn.x[ i ][ j ][ k ] = 0.;
				cloudn.x[ i ][ j ][ k ] = 0.;
				icen.x[ i ][ j ][ k ] = 0.;
				co2n.x[ i ][ j ][ k ] = 0.;

				p_dyn.x[ i ][ j ][ k ] = 0.;
				p_dynn.x[ i ][ j ][ k ] = 0.;
				p_stat.x[ i ][ j ][ k ] = 0.;

				P_rain.x[ i ][ j ][ k ] = 0.;
				P_snow.x[ i ][ j ][ k ] = 0.;
				S_c_c.x[ i ][ j ][ k ] = 0.;
				S_v.x[ i ][ j ][ k ] = 0.;
				S_c.x[ i ][ j ][ k ] = 0.;
				S_i.x[ i ][ j ][ k ] = 0.;
				S_r.x[ i ][ j ][ k ] = 0.;
				S_s.x[ i ][ j ][ k ] = 0.;

				rhs_t.x[ i ][ j ][ k ] = 0.;
				rhs_u.x[ i ][ j ][ k ] = 0.;
				rhs_v.x[ i ][ j ][ k ] = 0.;
				rhs_w.x[ i ][ j ][ k ] = 0.;
				rhs_c.x[ i ][ j ][ k ] = 0.;
				rhs_cloud.x[ i ][ j ][ k ] = 0.;
				rhs_ice.x[ i ][ j ][ k ] = 0.;
				rhs_co2.x[ i ][ j ][ k ] = 0.;

				aux_u.x[ i ][ j ][ k ] = 0.;
				aux_v.x[ i ][ j ][ k ] = 0.;
				aux_w.x[ i ][ j ][ k ] = 0.;

				Q_Latent.x[ i ][ j ][ k ] = 0.;
				Q_Sensible.x[ i ][ j ][ k ] = 0.;
				BuoyancyForce.x[ i ][ j ][ k ] = 0.;
				epsilon_3D.x[ i ][ j ][ k ] = 0.;
				radiation_3D.x[ i ][ j ][ k ] = 0.;
			}
		}
	}


//  final remarks
	cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;

	if ( velocity_iter == velocity_iter_max )   cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of maximum time steps ***** \n\n" << endl;

	if ( emin <= epsres )		cout << "***** steady solution reached! *****" << endl;
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


