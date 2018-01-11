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
#include "tinyxml2.h"
#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;

const double cAtmosphereModel::pi180 = 180./ M_PI;      // pi180 = 57.3

const double cAtmosphereModel::the_degree = 1.;         // compares to 1° step size laterally
const double cAtmosphereModel::phi_degree = 1.;         // compares to 1° step size longitudinally

const double cAtmosphereModel::dthe = the_degree / pi180; // dthe = the_degree / pi180 = 1.0 / 57.3 = 0.01745, 180 * .01745 = 3.141
const double cAtmosphereModel::dphi = phi_degree / pi180; // dphi = phi_degree / pi180 = 1.0 / 57.3 = 0.01745, 360 * .01745 = 6.282
    
const double cAtmosphereModel::dr = 0.025;    // 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
const double cAtmosphereModel::dt = 0.00001;  // time step coincides with the CFL condition
    
const double cAtmosphereModel::the0 = 0.;             // North Pole
const double cAtmosphereModel::phi0 = 0.;             // zero meridian in Greenwich
const double cAtmosphereModel::r0 = 1.;       /* earth's radius is r_earth = 6731 km, here it is assumed to be infinity, 
                                                circumference of the earth 40074 km */

cAtmosphereModel::cAtmosphereModel():
    j_res(0), 
    k_res(0),
    im_tropopause(new int [ jm ]),                        
    rad(im, 1.),
    the(jm, 2.),
    phi(km, 3.),
    Topography(jm, km, 0.),
    value_top(jm, km, 0.),
    Vegetation(jm, km, 0.),
    LatentHeat(jm, km, 0.),
    Condensation(jm, km, 0.),
    Evaporation(jm, km, 0.),
    Precipitation(jm, km, 0.),
    precipitable_water(jm, km, 0.),
    precipitation_NASA(jm, km, 0.),
    Ik(jm, km, 0.),
    Radiation_Balance(jm, km, 0.),
    temperature_NASA(jm, km, 0.),
    albedo(jm, km, 0.),
    epsilon(jm, km, 0.),
    Q_Radiation(jm, km, 0.),
    Q_Evaporation(jm, km, 0.),
    Q_latent(jm, km, 0.),
    Q_sensible(jm, km, 0.),
    Q_bottom(jm, km, 0.),
    Evaporation_Haude(jm, km, 0.),
    Evaporation_Penman(jm, km, 0.),
    co2_total(jm, km, 0.),
    h(im, jm, km, 0.),
    t(im, jm, km, 0.),    
    u(im, jm, km, 0.),
    v(im, jm, km, 0.),
    w(im, jm, km, 0.),
    c(im, jm, km, 0.),
    cloud(im, jm, km, 0.),
    ice(im, jm, km, 0.),
    co2(im, jm, km, 0.),
    tn(im, jm, km, 0.),
    un(im, jm, km, 0.),
    vn(im, jm, km, 0.),
    wn(im, jm, km, 0.),
    cn(im, jm, km, 0.),
    cloudn(im, jm, km, 0.),
    icen(im, jm, km, 0.),
    co2n(im, jm, km, 0.),
    p_dyn(im, jm, km, 0.),
    p_dynn(im, jm, km, 0.),
    p_stat(im, jm, km, 0.),
    rhs_t(im, jm, km, 0.),
    rhs_u(im, jm, km, 0.),
    rhs_v(im, jm, km, 0.),
    rhs_w(im, jm, km, 0.),
    rhs_p(im, jm, km, 0.),
    rhs_c(im, jm, km, 0.),
    rhs_cloud(im, jm, km, 0.),
    rhs_ice(im, jm, km, 0.),
    rhs_co2(im, jm, km, 0.),
    aux_u(im, jm, km, 0.),
    aux_v(im, jm, km, 0.),
    aux_w(im, jm, km, 0.),
    Latency(im, jm, km, 0.),
    Q_Sensible(im, jm, km, 0.),
    IceLayer(im, jm, km, 0.),
    t_cond_3D(im, jm, km, 0.),
    t_evap_3D(im, jm, km, 0.),
    BuoyancyForce(im, jm, km, 0.),
    epsilon_3D(im, jm, km, 0.),
    radiation_3D(im, jm, km, 0.),
    P_rain(im, jm, km, 0.),
    P_snow(im, jm, km, 0.),
    S_v(im, jm, km, 0.),
    S_c(im, jm, km, 0.),
    S_i(im, jm, km, 0.),
    S_r(im, jm, km, 0.),
    S_s(im, jm, km, 0.),
    S_c_c(im, jm, km, 0.)
{
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.

    PythonStream::OverrideCout();

    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

    // set default configuration
    SetDefaultConfig();

    coeff_mmWS = r_air / r_water_vapour; // coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]
    
    t=Array(im, jm, km, ta); // temperature
    u=Array(im, jm, km, ua); // u-component velocity component in r-direction
    v=Array(im, jm, km, va); // v-component velocity component in theta-direction
    w=Array(im, jm, km, wa); // w-component velocity component in phi-direction
    c=Array(im, jm, km, ca); // water vapour
    co2=Array(im, jm, km, coa); // CO2

    tn=Array(im, jm, km, ta); // temperature new
    un=Array(im, jm, km, ua); // u-velocity component in r-direction new
    vn=Array(im, jm, km, va); // v-velocity component in theta-direction new
    wn=Array(im, jm, km, wa); // w-velocity component in phi-direction new
    cn=Array(im, jm, km, ca); // water vapour new
    co2n=Array(im, jm, km, coa); // CO2 new

    p_dyn=Array(im, jm, km, pa); // dynamic pressure
    p_dynn=Array(im, jm, km, pa); // dynamic pressure
    p_stat=Array(im, jm, km, pa); // static pressure


//#include "PythonStream.h"
//	PythonStream buff;
//	std::ostream RogerStream ( &buff );
//	RogerStream << 1000;
//	RogerStream.flush();
}


cAtmosphereModel::~cAtmosphereModel() {
    delete [] im_tropopause;
}
 
#include "cAtmosphereDefaults.cpp.inc"


void cAtmosphereModel::LoadConfig ( const char *filename ) {
	XMLDocument doc;
	XMLError err = doc.LoadFile ( filename );

	if (err) {
		doc.PrintError();
		throw std::invalid_argument("   couldn't load config file inside cAtmosphereModel");
	}


	XMLElement *atom = doc.FirstChildElement("atom");
	if (!atom) {
		return;
	}


	XMLElement* elem_common = doc.FirstChildElement( "atom" )->FirstChildElement( "common" );
	if (!elem_common) {
		return;
	}


	XMLElement* elem_atmosphere = doc.FirstChildElement( "atom" )->FirstChildElement( "atmosphere" );
	if (!elem_atmosphere) {
		return;
	}

#include "AtmosphereLoadConfig.cpp.inc"

}

void cAtmosphereModel::RunTimeSlice ( int Ma )
{
// maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km )
// maximum number of overall iterations ( n )
// maximum number of inner velocity loop iterations ( velocity_iter_max )
// maximum number of outer pressure loop iterations ( pressure_iter_max )

	mkdir(output_path.c_str(), 0777);

	Array IceLayer(im, jm, km, 0.); // ice shield
	Array t_cond_3D(im, jm, km, 0.); // condensation temperature
	Array t_evap_3D(im, jm, km, 0.); // evaporation temperature
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
	int n_pres = 2;
	int velocity_iter = 0;
	int velocity_iter_2D = 0;
	int pressure_iter = 0;
	int pressure_iter_2D = 0;
	int switch_2D = 0;
	double residuum;
	double residuum_old = 0.;

//	radial expansion of the computational field for the computation of initial values
	int i_max = 32;		 // corresponds to about 12.8 km above sea level, maximum hight of the tropopause at equator
	int i_beg = 16;		 // corresponds to about 6.4 km above sea level, maximum hight of the tropopause at poles

//	naming a file to read the surface temperature by NASA of the modern world
	string Name_NASAbasedSurfaceTemperature_File = output_path + "/NASAbasedSurfaceTemperature.xyz";
	if(Ma == 0){
        	Name_NASAbasedSurfaceTemperature_File =  "../data/NASA_based_SurfaceTemperature.xyz";
     	}

//	naming a file to read the surface temperature by NASA of the modern world
	string Name_SurfaceTemperature_File = temperature_file;	
 	if(Ma != 0 && use_earthbyte_reconstruction){
        	Name_SurfaceTemperature_File = output_path + "/" + std::to_string(Ma) + "Ma_SurfaceTemperature.xyz";
     	}

// naming a file to read the surface precipitation by NASA
	string Name_SurfacePrecipitation_File = precipitation_file;	
	if(Ma != 0 && use_earthbyte_reconstruction){
		Name_SurfacePrecipitation_File = output_path + "/" + std::to_string(Ma) + "Ma_SurfacePrecipitation.xyz";
 	}

	string bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
	string bathymetry_filepath = bathymetry_path + "/" + bathymetry_name;

	cout << "\n   Output is being written to " << output_path << "\n";
	cout << "   Ma = " << Ma << "\n";
	cout << "   bathymetry_path = " << bathymetry_path << "\n";
	cout << "   bathymetry_filepath = " << bathymetry_filepath << "\n\n";
	cout << "   Name_SurfaceTemperature_File = " << Name_SurfaceTemperature_File << std::endl;
	cout << "   Name_SurfacePrecipitation_File = " << Name_SurfacePrecipitation_File << std::endl;
	
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
	RHS_Atmosphere						prepare ( im, jm, km, dt, dr, dthe, dphi, re, ec, sc_WaterVapour, sc_CO2, g, pr, omega, coriolis, centrifugal, WaterVapour, buoyancy, CO2, gam, sigma, lamda );
	RHS_Atmosphere						prepare_2D ( jm, km, dthe, dphi, re, omega, coriolis, centrifugal );

//	class RungeKutta_Atmosphere for the explicit solution of the Navier-Stokes equations
	RungeKutta_Atmosphere			result ( im, jm, km, dt, dr, dphi, dthe );

//	class Results_MSL_Atm to compute and show results on the mean sea level, MSL
	Results_MSL_Atm						calculate_MSL ( im, jm, km, sun, g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_equator, lv, ls, cp_l, L_atm, dt, dr, dthe, dphi, r_air, R_Air, r_water, r_water_vapour, R_WaterVapour, co2_vegetation, co2_ocean, co2_land, gam, t_pole, t_cretaceous );

//	class Pressure for the subsequent computation of the pressure by a separate Euler equation
	Pressure_Atm										startPressure ( im, jm, km, dr, dthe, dphi );

//	class BC_Thermo for the initial and boundary conditions of the flow properties
	BC_Thermo									circulation ( output_path, im, jm, km, i_beg, i_max, RadiationModel, NASATemperature, sun, declination, sun_position_lat, sun_position_lon, Ma, Ma_max, Ma_max_half, dt, dr, dthe, dphi, g, ep, hp, u_0, p_0, t_0, c_0, sigma, lv, ls, cp_l, L_atm, r_air, R_Air, r_water_vapour, R_WaterVapour, co2_0, co2_cretaceous, co2_vegetation, co2_ocean, co2_land, c_tropopause, co2_tropopause, c_ocean, c_land, t_average, co2_average, co2_pole, t_cretaceous, t_cret_cor, t_cretaceous_max, t_land, t_tropopause, t_equator, t_pole, gam, epsilon_equator, epsilon_pole, epsilon_tropopause, albedo_equator, albedo_pole, ik_equator, ik_pole );

//	class Restore to restore the iterational values from new to old
	Restore_Atm								oldnew( im, jm, km );


//	class element calls for the preparation of initial conditions for the flow properties

//	class element for the tropopause location as a parabolic distribution from pole to pole 
	circulation.TropopauseLocation ( im_tropopause );

//  class element for the surface temperature from NASA for comparison
	if ( Ma == 0 || use_earthbyte_reconstruction){ 
            circulation.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, temperature_NASA, t );
	}
	

//  class element for the surface temperature based on NASA temperature for progressing timeslices
//	if ( ( Ma > 0 ) && ( NASATemperature == 1 ) ) circulation.BC_NASAbasedSurfTempRead ( Name_NASAbasedSurfaceTemperature_File, t_cretaceous, t_cret_cor, t, c, cloud, ice );

//  class element for the surface precipitation from NASA for comparison
	if ( Ma == 0 || use_earthbyte_reconstruction){ 
		circulation.BC_Surface_Precipitation_NASA ( Name_SurfacePrecipitation_File, precipitation_NASA );
        }
//  class element for the parabolic temperature distribution from pol to pol, maximum temperature at equator
	circulation.BC_Temperature ( im_tropopause, temperature_NASA, h, t, p_dyn, p_stat );
	t_cretaceous = circulation.out_t_cretaceous (  );



//  class element for the correction of the temperature initial distribution around coasts
    if ( ( NASATemperature == 1 ) && ( Ma > 0 ) && !use_earthbyte_reconstruction ){ 
        circulation.IC_Temperature_WestEastCoast ( h, t );
    }

//  class element for the surface pressure computed by surface temperature with gas equation
	circulation.BC_Pressure ( p_stat, p_dyn, t, h );

//  parabolic water vapour distribution from pol to pol, maximum water vapour volume at equator
	circulation.BC_WaterVapour ( h, t, c );

//  class element for the parabolic CO2 distribution from pol to pol, maximum CO2 volume at equator
	circulation.BC_CO2 ( Vegetation, h, t, p_dyn, co2 );
	co2_cretaceous = circulation.out_co2 (  );

// class element for the surface temperature computation by radiation flux density
    if ( RadiationModel == 1 ){ 
        circulation.BC_Radiation_multi_layer(albedo, Ik, p_stat,t, c, h.to_Int3DArray(), 
                                             epsilon_3D, radiation_3D, cloud, ice );
    }

// 	class element for the initial conditions for u-v-w-velocity components
	circulation.IC_CellStructure ( im_tropopause, h, u, v, w );

// class element for the storing of velocity components, pressure and temperature for iteration start
	oldnew.restoreOldNew_3D(.9, u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, p_dynn, cn, cloudn, icen, co2n);
	oldnew.restoreOldNew_2D(.9, v, w, p_dyn, p_dynn, vn, wn);



// ***********************************   start of pressure and velocity iterations ************************************************************************

	double emin = epsres * 100.;

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

				cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations " << endl << endl << " Ma = " << Ma << "     n = " << n << "     n_pres = " << n_pres << "    velocity_iter_max_2D = " << velocity_iter_max_2D << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " << pressure_iter_max_2D << "    pressure_iter_2D = " << pressure_iter_2D << endl;

//	class BC_Atmosphaere for the geometry of a shell of a sphere
				boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
				boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

//	old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
				Accuracy_Atm		min_Residuum_old_2D ( im, jm, km, dthe, dphi );
				min_Residuum_old_2D.residuumQuery_2D ( rad, the, v, w );
				emin = min_Residuum_old_2D.out_min (  );

				residuum_old = emin;

//	class RungeKutta for the solution of the differential equations describing the flow properties
				result.solveRungeKutta_2D_Atmosphere ( prepare_2D, n, n_pres, r_air, u_0, p_0, L_atm, rad, the, rhs_v, rhs_w, rhs_p, h, v, w, p_dyn, vn, wn, p_dynn, aux_v, aux_w );

//	class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
				LandArea.BC_SolidGround ( RadiationModel, Ma, i_max, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, co2_cretaceous, pa, gam, sigma, h, u, v, w, t, p_dyn, c, cloud, ice, co2, radiation_3D, Vegetation );

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

				oldnew.restoreOldNew_2D ( 1., v, w, p_dyn, p_dynn, vn, wn );

				n++;
			}

//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
		if ( pressure_iter_2D == 1 ) 			startPressure.computePressure_2D ( pa, rad, the, p_dyn, p_dynn, h, rhs_v, rhs_w, aux_v, aux_w );


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

	n_pres = 2;
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

//  restoring the velocity component and the temperature for the new time step
			oldnew.restoreOldNew_3D ( 1., u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, p_dynn, cn, cloudn, icen, co2n );

			cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations " << endl << endl << " Ma = " << Ma << "     n = " << n << "     n_pres = " << n_pres << "    velocity_iter_max = " << velocity_iter_max << "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << "    pressure_iter = " << pressure_iter << endl;

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
			result.solveRungeKutta_3D_Atmosphere ( prepare, n, n_pres, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, rhs_t, rhs_u, rhs_v, rhs_w, rhs_p, rhs_c, rhs_cloud, rhs_ice, rhs_co2, h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, p_dynn, cn, cloudn, icen, co2n, aux_u, aux_v, aux_w, Latency, t_cond_3D, t_evap_3D, IceLayer, BuoyancyForce, Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c, Topography );

//	class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
			LandArea.BC_SolidGround ( RadiationModel, Ma, i_max, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, co2_cretaceous, pa, gam, sigma, h, u, v, w, t, p_dyn, c, cloud, ice, co2, radiation_3D, Vegetation );

// class element for the surface temperature computation by radiation flux density
            if ( RadiationModel == 1 ){
                circulation.BC_Radiation_multi_layer(albedo, Ik, p_stat, t, c, 
                                                     h.to_Int3DArray(), epsilon_3D, radiation_3D, 
                                                     cloud, ice );
            }


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
			circulation.Latent_Heat ( rad, the, phi, h, t, tn, u, v, w, p_dyn, p_stat, c, Latency, Q_Sensible, t_cond_3D, t_evap_3D, radiation_3D );

//	searching of maximum and minimum values of temperature
			string str_max_temperature = " max 3D temperature ", str_min_temperature = " min 3D temperature ", str_unit_temperature = "C";
			MinMax_Atm		minmaxTemperature ( im, jm, km );
			minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

//	searching of maximum and minimum values of u-component
			string str_max_u = " max 3D u-component ", str_min_u = " min 3D u-component ", str_unit_u = "m/s";
			MinMax_Atm		minmax_u ( im, jm, km );
			minmax_u.searchMinMax_3D ( str_max_u, str_min_u, str_unit_u, u, h );

//	searching of maximum and minimum values of v-component
			string str_max_v = " max 3D v-component ", str_min_v = " min 3D v-component ", str_unit_v = "m/s";
			MinMax_Atm		minmax_v ( im, jm, km );
			minmax_v.searchMinMax_3D ( str_max_v, str_min_v, str_unit_v, v, h );

//	searching of maximum and minimum values of w-component
			string str_max_w = " max 3D w-component ", str_min_w = " min 3D w-component ", str_unit_w = "m/s";
			MinMax_Atm		minmax_w ( im, jm, km );
			minmax_w.searchMinMax_3D ( str_max_w, str_min_w, str_unit_w, w, h );

//	searching of maximum and minimum values of dynamic pressure
			string str_max_pressure = " max 3D pressure dynamic ", str_min_pressure = " min 3D pressure dynamic ", str_unit_pressure = "hPa";
			MinMax_Atm	minmaxPressure ( im, jm, km );
			minmaxPressure.searchMinMax_3D ( str_max_pressure, str_min_pressure, str_unit_pressure, p_dyn, h );

//	searching of maximum and minimum values of static pressure
			string str_max_pressure_stat = " max 3D pressure static ", str_min_pressure_stat = " min 3D pressure static ", str_unit_pressure_stat = "hPa";
			MinMax_Atm		minmaxPressure_stat ( im, jm, km );
			minmaxPressure_stat.searchMinMax_3D ( str_max_pressure_stat, str_min_pressure_stat, str_unit_pressure_stat, p_stat, h );

			cout << endl << " energies in the three dimensional space: " << endl << endl;

//	searching of maximum and minimum values of radiation_3D
			string str_max_radiation_3D = " max 3D radiation ", str_min_radiation_3D = " min 3D radiation ", str_unit_radiation_3D = "W/m2";
			MinMax_Atm		minmaxLatency ( im, jm, km );
			minmaxLatency.searchMinMax_3D ( str_max_radiation_3D, str_min_radiation_3D, str_unit_radiation_3D, radiation_3D, h );

//	searching of maximum and minimum values of sensible heat
			string str_max_Q_Sensible = " max 3D sensible heat ", str_min_Q_Sensible = " min 3D sensible heat ", str_unit_Q_Sensible = "W/m2";
			MinMax_Atm	minmaxQ_Sensible ( im, jm, km );
			minmaxQ_Sensible.searchMinMax_3D ( str_max_Q_Sensible, str_min_Q_Sensible, str_unit_Q_Sensible, Q_Sensible, h );

//	searching of maximum and minimum values of latency
			string str_max_latency = " max 3D latent heat ", str_min_latency = " min 3D latent heat ", str_unit_latency = "W/m2";
			MinMax_Atm		minmaxRadiation ( im, jm, km );
			minmaxRadiation.searchMinMax_3D ( str_max_latency, str_min_latency, str_unit_latency, Latency, h );

//	searching of maximum and minimum values of t_cond_3D
			string str_max_t_cond_3D = " max 3D condensation temp ", str_min_t_cond_3D = " min 3D condensation temp ", str_unit_t_cond_3D = "C";
			MinMax_Atm		minmaxt_cond_3D ( im, jm, km );
			minmaxt_cond_3D.searchMinMax_3D ( str_max_t_cond_3D, str_min_t_cond_3D, str_unit_t_cond_3D, t_cond_3D, h );

//	searching of maximum and minimum values of t_evap_3D
			string str_max_t_evap_3D = " max 3D evaporation temp ", str_min_t_evap_3D = " min 3D evaporation temp ", str_unit_t_evap_3D = "C";
			MinMax_Atm		minmaxt_evap_3D ( im, jm, km );
			minmaxt_evap_3D.searchMinMax_3D ( str_max_t_evap_3D, str_min_t_evap_3D, str_unit_t_evap_3D, t_evap_3D, h );

			cout << endl << " greenhouse gases: " << endl << endl;

//	searching of maximum and minimum values of water vapour
			string str_max_water_vapour = " max 3D water vapour ", str_min_water_vapour = " min 3D water vapour ", str_unit_water_vapour = "g/kg";
			MinMax_Atm		minmaxWaterVapour ( im, jm, km );
			minmaxWaterVapour.searchMinMax_3D ( str_max_water_vapour, str_min_water_vapour, str_unit_water_vapour, c, h );

//	searching of maximum and minimum values of cloud water
			string str_max_cloud_water = " max 3D cloud water ", str_min_cloud_water = " min 3D cloud water ", str_unit_cloud_water = "g/kg";
			MinMax_Atm		minmaxCloudWater ( im, jm, km );
			minmaxCloudWater.searchMinMax_3D ( str_max_cloud_water, str_min_cloud_water, str_unit_cloud_water, cloud, h );

//	searching of maximum and minimum values of cloud ice
			string str_max_cloud_ice = " max 3D cloud ice ", str_min_cloud_ice = " min 3D cloud ice ", str_unit_cloud_ice = "g/kg";
			MinMax_Atm		minmaxCloudIce ( im, jm, km );
			minmaxCloudIce.searchMinMax_3D ( str_max_cloud_ice, str_min_cloud_ice, str_unit_cloud_ice, ice, h );

//	searching of maximum and minimum values of rain precipitation
			string str_max_P_rain = " max 3D rain ", str_min_P_rain = " min 3D rain ", str_unit_P_rain = "g/kg";
			MinMax_Atm		minmaxPRain ( im, jm, km );
			minmaxPRain.searchMinMax_3D ( str_max_P_rain, str_min_P_rain, str_unit_P_rain, P_rain, h );

//	searching of maximum and minimum values of snow precipitation
			string str_max_P_snow = " max 3D snow ", str_min_P_snow = " min 3D snow ", str_unit_P_snow = "g/kg";
			MinMax_Atm		minmaxPSnow ( im, jm, km );
			minmaxPSnow.searchMinMax_3D ( str_max_P_snow, str_min_P_snow, str_unit_P_snow, P_snow, h );

//	searching of maximum and minimum values of co2
			string str_max_co2 = " max 3D co2 ", str_min_co2 = " min 3D co2 ", str_unit_co2 = "ppm";
			MinMax_Atm		minmaxCO2 ( im, jm, km );
			minmaxCO2.searchMinMax_3D ( str_max_co2, str_min_co2, str_unit_co2, co2, h );

//	searching of maximum and minimum values of epsilon
			string str_max_epsilon = " max 3D epsilon ", str_min_epsilon = " min 3D epsilon ", str_unit_epsilon = "%";
			MinMax_Atm		minmaxEpsilon_3D ( im, jm, km );
			minmaxEpsilon_3D.searchMinMax_3D ( str_max_epsilon, str_min_epsilon, str_unit_epsilon, epsilon_3D, h );

//	searching of maximum and minimum values of buoyancy force
			string str_max_buoyancy_force = " max 3D buoyancy force ", str_min_buoyancy_force = " min 3D buoyancy force ", str_unit_buoyancy_force = "N/m2";
			MinMax_Atm		minmaxBuoyancyForce ( im, jm, km );
			minmaxBuoyancyForce.searchMinMax_3D ( str_max_buoyancy_force, str_min_buoyancy_force, str_unit_buoyancy_force, BuoyancyForce, h );



// 2D-fields

//	searching of maximum and minimum values of co2 total
			cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" << endl << " results based on two dimensional considerations of the problem" << endl;

			cout << endl << " co2 distribution columnwise: " << endl << endl;

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

/*
//	searching of maximum and minimum values of NASA precipitation
			if ( Ma == 0 )
			{
				string str_max_precipitation = " max precipitation_NASA ", str_min_precipitation = " min precipitation_NASA ", str_unit_precipitation = "mm";
				MinMax	minmaxPrecipitation_NASA ( jm, km, coeff_mmWS );
				minmaxPrecipitation_NASA.searchMinMax_2D ( str_max_precipitation, str_min_precipitation, str_unit_precipitation, precipitation_NASA, h );
			}
*/

//	searching of maximum and minimum values of precipitable water
			string str_max_precipitable_water = " max precipitable water ", str_min_precipitable_water = " min precipitable water ", str_unit_precipitable_water = "mm";
			MinMax_Atm		minmaxPrecipitable_water ( jm, km, coeff_mmWS );
			minmaxPrecipitable_water.searchMinMax_2D ( str_max_precipitable_water, str_min_precipitable_water, str_unit_precipitable_water, precipitable_water, h );

			cout << endl << " energies at see level without convection influence: " << endl << endl;

//	searching of maximum and minimum values of radiation balance
//	string str_max_Radiation_Balance = " max radiation balance ", str_min_Radiation_Balance = " min radiation balance ", str_unit_Radiation_Balance = "W/m2";
//	MinMax	minmaxRadiation_Balance ( jm, km, coeff_mmWS );
//	minmaxRadiation_Balance.searchMinMax_2D ( str_max_Radiation_Balance, str_min_Radiation_Balance, str_unit_Radiation_Balance, Radiation_Balance, h );
//	min_Radiation_Balance = minmaxRadiation_Balance.out_minValue (  );

//	searching of maximum and minimum values of radiation
			string str_max_Q_Radiation = " max 2D Q radiation ", str_min_Q_Radiation = " min 2D Q radiation ", str_unit_Q_Radiation = "W/m2";
			MinMax_Atm		minmaxQ_Radiation ( jm, km, coeff_mmWS );
			minmaxQ_Radiation.searchMinMax_2D ( str_max_Q_Radiation, str_min_Q_Radiation, str_unit_Q_Radiation, Q_Radiation, h );

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


//	searching of maximum and minimum values of latent heat
			string str_max_LatentHeat = " max 2D latent heat ", str_min_LatentHeat = " min 2D latent heat ", str_unit_LatentHeat = "W/m2";
			MinMax_Atm		minmaxLatentHeat_2D ( jm, km, coeff_mmWS );
			minmaxLatentHeat_2D.searchMinMax_2D ( str_max_LatentHeat, str_min_LatentHeat, str_unit_LatentHeat, LatentHeat, h );

//	searching of maximum and minimum values of t_cond_2D
			string str_max_t_cond = " max 2D Condensation ", str_min_t_cond = " min 2D Condensation ", str_unit_t_cond = "W/m2";
			MinMax_Atm		minmaxt_cond_2D ( jm, km, coeff_mmWS );
			minmaxt_cond_2D.searchMinMax_2D ( str_max_t_cond, str_min_t_cond, str_unit_t_cond, Condensation, h );


//	searching of maximum and minimum values of t_evap_2D
			string str_max_t_evap = " max 2D Evaporation ", str_min_t_evap = " min 2D Evaporation ", str_unit_t_evap = "W/m2";
			MinMax_Atm		minmaxt_evap_2D ( jm, km, coeff_mmWS );
			minmaxt_evap_2D.searchMinMax_2D ( str_max_t_evap, str_min_t_evap, str_unit_t_evap, Evaporation, h );

			cout << endl << " secondary data: " << endl << endl;

/*
//	searching of maximum and minimum values of Evaporation
			string str_max_heat_t_Evaporation = " max heat Evaporation ", str_min_heat_t_Evaporation = " min heat Evaporation ", str_unit_heat_t_Evaporation = " W/m2";
			MinMax		minmaxQ_t_Evaporation ( jm, km, coeff_mmWS );
			minmaxQ_t_Evaporation.searchMinMax_2D ( str_max_heat_t_Evaporation, str_min_heat_t_Evaporation, str_unit_heat_t_Evaporation, Q_Evaporation, h );
*/
/*
//	searching of maximum and minimum values of Evaporation by Haude
			string str_max_t_Evaporation_Haude = " max Evaporation Haude ", str_min_t_Evaporation_Haude = " min Evaporation Haude ", str_unit_t_Evaporation_Haude = "mm/d";
			MinMax		minmaxt_Evaporation_Haude ( jm, km, coeff_mmWS );
			minmaxt_Evaporation_Haude.searchMinMax_2D ( str_max_t_Evaporation_Haude, str_min_t_Evaporation_Haude, str_unit_t_Evaporation_Haude, Evaporation_Haude, h );
*/
//	searching of maximum and minimum values of Evaporation by Penman
			string str_max_t_Evaporation_Penman = " max Evaporation Penman ", str_min_t_Evaporation_Penman = " min Evaporation Penman ", str_unit_t_Evaporation_Penman = "mm/d";
			MinMax_Atm		minmaxt_Evaporation_Penman ( jm, km, coeff_mmWS );
			minmaxt_Evaporation_Penman.searchMinMax_2D ( str_max_t_Evaporation_Penman, str_min_t_Evaporation_Penman, str_unit_t_Evaporation_Penman, Evaporation_Penman, h );

			cout << endl << " properties of the atmosphere at the surface: " << endl << endl;

//	searching of maximum and minimum values of albedo
			string str_max_albedo = " max 2D albedo ", str_min_albedo = " min 2D albedo ", str_unit_albedo = "%";
			MinMax_Atm		minmaxAlbedo ( jm, km, coeff_mmWS );
			minmaxAlbedo.searchMinMax_2D ( str_max_albedo, str_min_albedo, str_unit_albedo, albedo, h );
/*
//	searching of maximum and minimum values of epsilon
			string str_max_epsilon = " max 2D epsilon ", str_min_epsilon = " min 2D epsilon ", str_unit_epsilon = "%";
			MinMax	minmaxEpsilon ( jm, km, coeff_mmWS );
			minmaxEpsilon.searchMinMax_2D ( str_max_epsilon, str_min_epsilon, str_unit_epsilon, epsilon, h );
*/

//	computation of vegetation areas
			LandArea.vegetationDistribution ( max_Precipitation, Precipitation, Vegetation, t, h );


//	composition of results
			calculate_MSL.run_MSL_data ( n, velocity_iter_max, RadiationModel, rad, the, phi, h, c, cn, co2, co2n, t, tn, p_dyn, p_stat, BuoyancyForce, u, v, w, Latency, Q_Sensible, radiation_3D, t_cond_3D, t_evap_3D, cloud, cloudn, ice, icen, P_rain, P_snow, aux_u, aux_v, aux_w, precipitation_NASA, Evaporation, Condensation, LatentHeat, precipitable_water, Q_Radiation, Q_Evaporation, Q_latent, Q_sensible, Q_bottom, Evaporation_Penman, Evaporation_Haude, Vegetation, Radiation_Balance, albedo, co2_total, Precipitation, S_v, S_c, S_i, S_r, S_s, S_c_c );


//  restoring the velocity component and the temperature for the new time step
			oldnew.restoreOldNew_3D ( 1., u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, p_dynn, cn, cloudn, icen, co2n );
//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//	Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation distribution formed of rain and snow
			if ( velocity_iter == velocity_n )
			{
				circulation.Two_Category_Ice_Scheme ( n, velocity_iter_max, RadiationModel, h, c, t, p_stat, cloud, ice, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );
				velocity_n = velocity_iter + 2;
			}
			n++;
		}
//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of velocity loop_3D: if ( velocity_iter > velocity_iter_max )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
		if ( pressure_iter == 1 ) 			startPressure.computePressure_3D ( pa, rad, the, p_dyn, p_dynn, h, rhs_u, rhs_v, rhs_w, aux_u, aux_v, aux_w );



//	Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation distribution formed of rain and snow
	circulation.Two_Category_Ice_Scheme ( n, velocity_iter_max, RadiationModel, h, c, t, p_stat, cloud, ice, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );


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


//	printout in ParaView files and sequel files

//	class PostProcess_Atmosphaere for the printing of results
	PostProcess_Atmosphere write_File ( im, jm, km, output_path );

//	writing of data in ParaView files
//	radial data along constant hight above ground
	int i_radial = 0;
	write_File.paraview_vtk_radial ( bathymetry_name, i_radial, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, t_cond_3D, t_evap_3D , BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain, P_snow, Evaporation, Condensation, precipitable_water, Q_bottom, Radiation_Balance, Q_Radiation, Q_latent, Q_sensible, Evaporation_Penman, Evaporation_Haude, Q_Evaporation, temperature_NASA, precipitation_NASA, Vegetation, albedo, epsilon, Precipitation );

//	londitudinal data along constant latitudes
	int j_longal = 62;			// Mount Everest/Himalaya
	write_File.paraview_vtk_longal ( bathymetry_name, j_longal, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, t_cond_3D, t_evap_3D, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain, P_snow );

	int k_zonal = 87;			// Mount Everest/Himalaya
	write_File.paraview_vtk_zonal ( bathymetry_name, k_zonal, n, hp, ep, R_Air, g, L_atm, u_0, t_0, p_0, r_air, c_0, co2_0, h, p_dyn, p_stat, t_cond_3D, t_evap_3D, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, radiation_3D, epsilon_3D, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c );

//	3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
	write_File.paraview_panorama_vts ( bathymetry_name, n, u_0, t_0, p_0, r_air, c_0, co2_0, h, t, p_dyn, p_stat, BuoyancyForce, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain, P_snow );

//	writing of v-w-data in the v_w_transfer file
	PostProcess_Atmosphere ppa ( im, jm, km, output_path );
	ppa.Atmosphere_v_w_Transfer ( bathymetry_name, v, w, p_dyn );
	ppa.Atmosphere_PlotData ( bathymetry_name, u_0, t_0, h, v, w, t, c, Precipitation, precipitable_water );

	t_cret_cor = t_cretaceous;

	if ( NASATemperature == 1 && !use_earthbyte_reconstruction){ 
		circulation.BC_NASAbasedSurfTempWrite ( Name_NASAbasedSurfaceTemperature_File, t_cretaceous, t_cret_cor, t, c, cloud, ice );
	}

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of pressure loop: if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::

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

			LatentHeat.y[ j ][ k ] = 0.;
			Condensation.y[ j ][ k ] = 0.;
			Evaporation.y[ j ][ k ] = 0.;

			Ik.y[ j ][ k ] = 0.;
			Radiation_Balance.y[ j ][ k ] = 0.;

			albedo.y[ j ][ k ] = 0.;
			epsilon.y[ j ][ k ] = 0.;

			Q_Radiation.y[ j ][ k ] = 0.;
			Q_Evaporation.y[ j ][ k ] = 0.;
			Q_latent.y[ j ][ k ] = 0.;
			Q_sensible.y[ j ][ k ] = 0.;
			Q_bottom.y[ j ][ k ] = 0.;

			Evaporation_Haude.y[ j ][ k ] = 0.;
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

				Latency.x[ i ][ j ][ k ] = 0.;
				Q_Sensible.x[ i ][ j ][ k ] = 0.;
				IceLayer.x[ i ][ j ][ k ] = 0.;
				t_cond_3D.x[ i ][ j ][ k ] = 0.;
				t_evap_3D.x[ i ][ j ][ k ] = 0.;
				BuoyancyForce.x[ i ][ j ][ k ] = 0.;
				epsilon_3D.x[ i ][ j ][ k ] = 0.;
				radiation_3D.x[ i ][ j ][ k ] = 0.;
			}
		}
	}


//   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//  final remarks
	cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;

	if ( velocity_iter == velocity_iter_max )   cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of maximum time steps ***** \n\n" << endl;

	if ( emin <= epsres )		cout << "***** steady solution reached! *****" << endl;
}




void cAtmosphereModel::Run() {
// create the output dir
	mkdir(output_path.c_str(), 0777);

	cout << "Output is being written to " << output_path << "\n";

// write out the config for reproducibility
// disabled for now
// std::stringstream output_config_path;
// output_config_path << output_path << "/config_atm.xml";
// WriteConfig(output_config_path.str().c_str());


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

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// time slices to be run after actualizing
	int i_time_slice_max = 15;
	int *time_slice = new int [ i_time_slice_max ];	 // time slices in Ma

	time_slice [ 0 ] = 0;								   // Golonka Bathymetry and Topography
	time_slice [ 1 ] = 10;
	time_slice [ 2 ] = 20;
	time_slice [ 3 ] = 30;
	time_slice [ 4 ] = 40;
	time_slice [ 5 ] = 50;
	time_slice [ 6 ] = 60;
	time_slice [ 7 ] = 70;
	time_slice [ 8 ] = 80;
	time_slice [ 9 ] = 90;
	time_slice [ 10 ] = 100;
	time_slice [ 11 ] = 110;
	time_slice [ 12 ] = 120;
	time_slice [ 13 ] = 130;
	time_slice [ 14 ] = 140;

//  choice of the time slice to be computed
	for (int i_time_slice = 0; i_time_slice < i_time_slice_max; i_time_slice++) {
		int Ma = time_slice[i_time_slice];
		cout << "Ma = " << Ma << "\n";

		RunTimeSlice(Ma);
	}
//   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//  final remarks
	cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;

	cout << endl;
	cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
	cout << "\n\n\n\n";

	delete [ ] time_slice;

}

void cAtmosphereModel::PrintMaxMinValues() {
/*//  searching of maximum and minimum values of temperature
            string str_max_temperature = " max 3D temperature ", str_min_temperature = " min 3D temperature ", str_unit_temperature = "C";
            MinMax_Atm      minmaxTemperature ( im, jm, km );
            minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

//  searching of maximum and minimum values of u-component
            string str_max_u = " max 3D u-component ", str_min_u = " min 3D u-component ", str_unit_u = "m/s";
            MinMax_Atm      minmax_u ( im, jm, km );
            minmax_u.searchMinMax_3D ( str_max_u, str_min_u, str_unit_u, u, h );

//  searching of maximum and minimum values of v-component
            string str_max_v = " max 3D v-component ", str_min_v = " min 3D v-component ", str_unit_v = "m/s";
            MinMax_Atm      minmax_v ( im, jm, km );
            minmax_v.searchMinMax_3D ( str_max_v, str_min_v, str_unit_v, v, h );

//  searching of maximum and minimum values of w-component
            string str_max_w = " max 3D w-component ", str_min_w = " min 3D w-component ", str_unit_w = "m/s";
            MinMax_Atm      minmax_w ( im, jm, km );
            minmax_w.searchMinMax_3D ( str_max_w, str_min_w, str_unit_w, w, h );

//  searching of maximum and minimum values of dynamic pressure
            string str_max_pressure = " max 3D pressure dynamic ", str_min_pressure = " min 3D pressure dynamic ", str_unit_pressure = "hPa";
            MinMax_Atm  minmaxPressure ( im, jm, km );
            minmaxPressure.searchMinMax_3D ( str_max_pressure, str_min_pressure, str_unit_pressure, p_dyn, h );

//  searching of maximum and minimum values of static pressure
            string str_max_pressure_stat = " max 3D pressure static ", str_min_pressure_stat = " min 3D pressure static ", str_unit_pressure_stat = "hPa";
            MinMax_Atm      minmaxPressure_stat ( im, jm, km );
            minmaxPressure_stat.searchMinMax_3D ( str_max_pressure_stat, str_min_pressure_stat, str_unit_pressure_stat, p_stat, h );

            cout << endl << " energies in the three dimensional space: " << endl << endl;

//  searching of maximum and minimum values of radiation_3D
            string str_max_radiation_3D = " max 3D radiation ", str_min_radiation_3D = " min 3D radiation ", str_unit_radiation_3D = "W/m2";
            MinMax_Atm      minmaxLatency ( im, jm, km );
            minmaxLatency.searchMinMax_3D ( str_max_radiation_3D, str_min_radiation_3D, str_unit_radiation_3D, radiation_3D, h );

//  searching of maximum and minimum values of sensible heat
            string str_max_Q_Sensible = " max 3D sensible heat ", str_min_Q_Sensible = " min 3D sensible heat ", str_unit_Q_Sensible = "W/m2";
            MinMax_Atm  minmaxQ_Sensible ( im, jm, km );
            minmaxQ_Sensible.searchMinMax_3D ( str_max_Q_Sensible, str_min_Q_Sensible, str_unit_Q_Sensible, Q_Sensible, h );
*/}
