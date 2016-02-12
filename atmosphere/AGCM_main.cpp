/*
 * Atmosphere General Circulation Model ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include <iomanip>
//#include <netcdf.h>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "BC_Atm.h"
#include "BC_Bath_Atm.h"
#include "BC_Thermo.h"
#include "Accuracy_Atm.h"
#include "RHS_Atm.h"
#include "RungeKutta_Atm.h"
#include "Print_Atm.h"
#include "PostProcess_Atm.h"
#include "Pressure_Atm.h"
#include "Restore_Atm.h"
#include "Results_Atm.h"
#include "MinMax_Atm.h"
//#include "File_NetCDF.h"


using namespace std;


// "SequelFile"				Default for a sequel file to be written
// "SequelFile 1"			a sequel file will be written
// "SequelFile 0"			a sequel file will not be written

// "omega"					Default for the rotation of the earth by omega = 7.29e-5 in radians for Coriolis und centrifugal forces
// "omega 1"				with earth rotation, Coriolis and centrifugal forces are different of zero
// "omega 0"				without earth rotation, Coriolis and centrifugal forces are different of zero

// "coriolis"					Computation with and without a rate of Coriolis force
// "coriolis 1"				with inclusion of the Coriolis force
// "coriolis 0"				without inclusion of the Coriolis force

// "centrifugal"				Computation with and without a rate of centrifugal force
// "centrifugal 1"			with inclusion of the centrifugal force
// "centrifugal 0"			without inclusion of the centrifugal force

// "WaterVapour"			Computation with and without a concentration equation for water vapor
// "WaterVapour 1"		with inclusion of water vapour
// "WaterVapour 0"		without inclusion of water vapour

// "CO2"						Computation with and without a concentration equation for co2
// "CO2 1"					with inclusion of co2
// "CO2 0"					without inclusion of co2

// "buoyancy"				Computation with and without a buoyancy term following the Boussinesq-Approximation
// "buoyancy 1"			with inclusion of a buoyancy term
// "buoyancy 0"			without inclusion of a buoyancy term

// "sun"							Computation with and without an averaged sun radiation
// "sun 1"						with inclusion of radiation
// "sun 0"						without inclusion of radiation

// "IceShield"				Computation with and without ice shields
// "IceShield 1"				with inclusion of ice shields
// "IceShield 0"				without inclusion of ice shields

// Earth's radius is r_earth = 6731 km compares to 6.731 [ / ]
// for 20 km expansion of the area of circulation compares to 0.02 [ / ] with 40 steps of size 0.0005 

// Definition of meridional and longitudinal step sizes 
// for example: dthe = the_Grad / pi180 = 1.125 / 57.3 = 0.01963

// maximum velocity in the subtropical jet  w_max = 1.2 [ / ] compares to 30 m/s = 108 km/h as annual mean 

// minimum temperature at the tropopause t_min = .789 compares to -60° C compares to 213 K
// temperature t_0 = 1.0 compares to 0° C compares to 273,15 K
// temperature t_0 = 0.9561 compares to -12° C compares to 261,15 K
// temperature difference from equator to pole   18°C compares to  t_delta = 0.0665  compares to  18 K

// rate of CO2 compares to 600Gt in the atmosphere ( 0.02% ), 100Gt consumed by vegetation ( 0.00333% ), ratio 1/6 = 0.166667 to be checked!

// maximum rate of water vapour for sea level at the equator c_equator = 40 g/kg  = 0.040 kg/kg ( mass of water vapour / mass of air 




int main ( int argc, char *argv[ ] )
{
// maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km ), 
// maximum number of overall iterations ( n ),
// maximum number of inner velocity loop iterations ( velocity_iter_max ),
// maximum number of outer pressure loop iterations ( pressure_iter_max )

	int im = 41, jm = 181, km = 361, nm = 200, velocity_iter_max = 10, pressure_iter_max = 5;
	int velocity_iter_max_2D = 10, pressure_iter_max_2D = 5;

	int n, i_radial, j_longal, k_zonal, i_max, i_beg;
	int velocity_iter, pressure_iter, pressure_iter_aux, velocity_iter_2D, pressure_iter_2D;
	int i_res, j_res, k_res;
	int i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_co2, j_co2, k_co2, i_p, j_p, k_p;
	int Ma, i_time_slice, i_time_slice_max;
	int switch_2D;

	double time;
	double residuum, residuum_old, min, min_u, min_v, min_w, min_t, min_c, min_co2, min_p;
	double max_Precipitation;

	const int SequelFile = 0;									// sequel file will be written
	const double coriolis = 1.;									// computation with Coriolis force
	const double centrifugal = 1.;							// computation with centrifugal force
	const double WaterVapour = 1.;						// computation with water vapour
	const double buoyancy = 1.;								// computation with buoyancy
	const double CO2 = 1.;									// computation with CO2

	const int declination = 0;									// position of sun axis, today 23,4° 
																				// 21.12.: -23,4°, am 21.3. und 23.9.: 0°
																				// 21.6.: +23,4°, in between sin form
	const int sun_position_lat = 0;							// position of sun j = 120 means 30°S, j = 60 means 30°N
	const int sun_position_lon = 180;						// position of sun k = 180 means 0° or 180° O ( Greenwich, zero meridian )

	const int  sun = 0;												// while no variable sun position wanted
	const int  IceShield = 0;										// while no ice shields wanted

	const int Ma_max = 300;									// parabolic temperature distribution 300 Ma (From Ruddiman)
	const int Ma_max_half = 150;							// half of time scale

	const double L_atm = 20000.;							// extension of the atmosphere shell in m
	const double dt = 0.0001;									// time step coincides with the CFL condition
	const double dr = 0.0005;									// compares to 500m hight
	const double the_Grad = 1.;								// compares to 1° step size laterally
	const double phi_Grad = 1.;								// compares to 1° step size longitudinally
	const double pi180 = 180./M_PI;						// pi180 = 57.3
	const double dthe = the_Grad / pi180, dphi = phi_Grad / pi180;//dthe = the_Grad / pi180 = 1.125 / 57.3 = 0.01963
	const double epsres = 0.00001;						// accuracy of relative and absolute errors

	const double re = 1000.;									// Reynolds number: ratio viscous to inertia forces, Re = u * L / nue
	const double ec = .00044;								// Eckert number: ratio kinetic energy to enthalpy, Ec  = u² / cp T
	const double sc_WaterVapour = .6;					// Schmidt number of water vapour
	const double sc_CO2 = .96;								// Schmidt number of CO2
	const double pr = .7179;									// Prandtl number of air for laminar flows
	const double gr = 9.8066;									// gravitational acceleration of the earth
	const double omega = 7.29e-5;						// rotation number of the earth
	const double ep = .623;									// ratio of the gas constants of dry air to water vapour
	const double hp = 6.1078;								// water vapour pressure at T = 0°C: E = 6.1 hPa
	const double u_0 = 15.;									// maximum value of velocity in m/s
	const double p_0 = 1013.25;							// pressure at sea level 1000 hPa
	const double t_0 = 273.15;								// temperature in K compare to 0°C
	const double ik = 1366.;									// solar constant in W/m2
	const double sigma = 5.670280e-8;					// Stefan-Boltzmann constant W/( m2*K4 )
	const double albedo = .3;									// capability of reflection of short wave radiation
	const double R_Air = 287.1;								// specific gas constant of air in J/( kg*K ))
	const double R_WaterVapour = 461.6;			// specific gas constant of water vapour in J/( kg*K ))
	const double R_co2 = 188.91;							// specific gas constant of CO2 in J/( kg*K ))
	const double lv = 2.5e6;									// specific latent evaporation heat ( condensation heat ) in J/kg
	const double ls = 2.83e6;									// specific latent vaporisation heat ( sublimation heat ) in J/kg
	const double cp_l = 1870.;								// specific heat capacity of water vapour at constant pressure and 20°C in J/( kg K )
	const double r_0_air = 1.2041;							// density of air in kg/m3 at 20°C
	const double r_0_water_vapour = 0.0094;		// density of water vapour in kg/m3 at 25°C and dewpoint temperature at 10°C
	const double r_0_co2 = 0.0019767;					// density of CO2 in kg/m3 at 25°C
	const double c_0 = .035;									// maximum value of water vapour in kg / kg
	const double co2_0 = 280.;								// maximum value of CO2 in ppm

	const double r0 = 6.731; 									// Earth's radius is r_earth = 6731 km compares to 6.731 [ / ]
	const double the0 = 0.;										// North Pole
	const double phi0 = 0.;										// zero meridian in Greenwich

	const double ua = 0.;										// initial velocity component in r-direction
	const double va = 0.;										// initial velocity component in theta-direction
	const double wa = 0.;										// initial velocity component in phi-direction
	const double pa = 0.;										// initial value for the pressure field
	const double ca = 0.;										// value 1.0 stands for the maximum value of 35 g/kg water vapour
	const double tau = 1.;										// initial value for the temperature field,1.0 compares to 0° C compares to 273.15 K
	const double tao = .78;										// initial value for the temperature field in the tropopause, tao = 0.78 compares to -60.093°C compares to 213,057 K
	const double t_Boussinesq = 1.07;					// compares to 20°C, threshhold temperature for the Boussinesq-approximation concerning bouyancy effect
	const double t_cretaceous_max = 10.;				// maximum add of mean temperature in °C during cretaceous times
	const double co2a = 280.;								// 280 ppm CO2 at preindustrial times
	const double coeff_mmWS = r_0_air / r_0_water_vapour / 1000.;	// air density [kg/m3] / water vapour density [kg/m3]   coeff_mmWS = 0.12766

	double t_Average = 15.;									// mean temperature of the modern earth
	double t_equator = 1.119;									// temperature t_0 = 1.119 compares to 32.5° C compares to 305.65 K
//	double t_equator = 1.139;									// temperature t_0 = 1.139 compares to 38.0° C compares to 311.15 K
//	double t_equator = 1.165;									// temperature t_0 = 1.139 compares to 45.0° C compares to 318.15 K
	double t_pole = .8;												// temperature at the poles t_pole = 0.8 compares to -54.63°C compares to 218.52 K
	double t_tropopause = .78;								// temperature in the tropopause
	double t_land_plus = .007322;							// temperature increase on land ( 1°C compares to t_land_plus = 0.003661 )

	double c_tropopause = 0.;									// minimum water vapour at tropopause c_tropopause = 0.01429 compares to 0.05 g/kg
//	double c_land_minus = .9;									// water vapour reduction on land ( 90% )
//	double c_ocean_minus = 1.;								// water vapour reduction on sea surface ( 100% )
	double c_land_minus = 1.3;								// water vapour reduction on land ( 90% )
	double c_ocean_minus = 1.4;							// water vapour reduction on sea surface ( 100% )

	double co2_Average = 280.;								// rate of CO2 at preindustrial times
	double co2_equator = 280.;								// maximum rate of CO2 at sea level at equator, 1. compares to 280 ppm
	double co2_tropopause = 280.;							// minimum rate of CO2 at tropopause  0 ppm
	double co2_pole = 260.;									// percentage of CO2 of the sea surface
	double co2_vegetation = 0.16667;					// value compares to 100/600Gt per year on the global surface for the purpose of testing
	double co2_ocean = 0.1;									// value compares to 0.6/600Gt per year on the sea surface for the purpose of testing
	double co2_land = 0.033;									// value compares to 0.2/600Gt per year on land for the purpose of testing

	double albedo_extra = .3;									// global albedo extraterrestric, capability of reflection
	double epsilon_atmos = .77;								// capability of emissions in the atmosphere


// time slices to be run after actualizing 
	i_time_slice_max = 15;
	int *time_slice = new int [ i_time_slice_max ]; 	// time slices in Ma

	time_slice [ 0 ] = 0;												// Golonka Bathymetry and Topography
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


// 	class Array for 1-D, 2-D and 3-D field declarations
// 1D arrays
	Array_1D	rad ( im, 0., r0, dr );							// radial coordinate direction
	Array_1D	the ( jm, 0., the0, dthe );					// lateral coordinate direction
	Array_1D	phi ( km, 0., phi0, dphi );					// longitudinal coordinate direction

// 2D arrays
	Array_2D	Vegetation ( jm, km, 0. );					// vegetation via precipitation
	Array_2D	Evaporation ( jm, km, 0. );				// areas of higher evaporation
	Array_2D	IceAir ( jm, km, 0. );							// areas of higher IceAir
	Array_2D	Condensation ( jm, km, 0. );			// areas of highercondensation
	Array_2D	Precipitation ( jm, km, 0. );				// areas of precipitation
	Array_2D	precipitable_water ( jm, km, 0. );		// areas of total water in the air
	Array_2D	t_j ( jm, km, 0. );								// auxilliar field for the surface temperature
	Array_2D	c_j ( jm, km, 0. );								// auxilliar field for water vapour  at the surface
	Array_2D	p_j ( jm, km, 0. );								// auxilliar field for pressure  at the surface
	Array_2D	co2_j ( jm, km, 0. );							// auxilliar field for rate of CO2 at the surface
	Array_2D	Ice_Balance ( jm, km, 0. );				// rate of the ice shield
	Array_2D	Ice_Balance_add ( jm, km, 0. );		// addition of all ice layers
	Array_2D	Ik ( jm, km, 0. );								// radiation in the atmosphere
	Array_2D	Radiation_Balance_atm ( jm, km, 0. );// radiation balance of the atmosphere
	Array_2D	Radiation_Balance_bot ( jm, km, 0. );// radiation balance of the ground
	Array_2D	temp_eff_atm ( jm, km, 0. );			// effektive temperature in the atmosphere
	Array_2D	temp_eff_bot ( jm, km, 0. );				// effektive temperature on the ground
	Array_2D	Q_Evaporation ( jm, km, 0. );			// evaporation heat of water by Kuttler
	Array_2D	Q_latent ( jm, km, 0. );					// latent heat from bottom values by the energy transport equation
	Array_2D	Q_sensible ( jm, km, 0. );				// sensible heat from bottom values by the energy transport equation
	Array_2D	Q_diff ( jm, km, 0. );						// difference by Q_Balance_Radiation - Q_latent - Q_sensible
	Array_2D	Evaporation_Haude ( jm, km, 0. );	// evaporation by Haude in [mm/d]
	Array_2D	Evaporation_Penman ( jm, km, 0. );// evaporation by Penman in [mm/d]
	Array_2D	Q_Balance_Radiation ( jm, km, 0. );// heat from the radiation balance in [W/m2]
	Array_2D	precipitation_j ( jm, km, 0. );			// surface precipitation from NASA
	Array_2D	Water_super ( jm, km, 0. );				// supercooled rain
	Array_2D	Water ( jm, km, 0. );						// rain water
	Array_2D	aux_2D_v ( jm, km, 0. );					// auxilliar field v
	Array_2D	aux_2D_w ( jm, km, 0. );					// auxilliar field w

// 3D arrays
	Array	t ( im, jm, km, tau );									// temperature
	Array	u ( im, jm, km, ua );									// u-component velocity component in r-direction
	Array	v ( im, jm, km, va );									// v-component velocity component in theta-direction
	Array	w ( im, jm, km, wa );								// w-component velocity component in phi-direction
	Array	p ( im, jm, km, pa );									// pressure
	Array	c ( im, jm, km, ca );									// salt
	Array	co2 ( im, jm, km, co2a );							// CO2
	Array	tn ( im, jm, km, tau );								// temperature new
	Array	un ( im, jm, km, ua );								// u-component velocity component in r-direction new
	Array	vn ( im, jm, km, va );								// v-component velocity component in theta-direction new
	Array	wn ( im, jm, km, wa );								// w-component velocity component in phi-direction new
	Array	pn ( im, jm, km, pa );								// pressure new
	Array	cn ( im, jm, km, ca );								// salt new
	Array	co2n ( im, jm, km, co2a );						// CO2 new
	Array	h ( im, jm, km, 0. );									// bathymetry, depth from sea level
	Array	rhs_t ( im, jm, km, 0. );								// auxilliar field RHS temperature
	Array	rhs_u ( im, jm, km, 0. );								// auxilliar field RHS
	Array	rhs_v ( im, jm, km, 0. );								// auxilliar field RHS
	Array	rhs_w ( im, jm, km, 0. );								// auxilliar field RHS
	Array	rhs_p ( im, jm, km, 0. );								// auxilliar field RHS pressure
	Array	rhs_c ( im, jm, km, 0. );								// auxilliar field RHS salt
	Array	rhs_co2 ( im, jm, km, 0. );							// auxilliar field RHS CO2
	Array	aux_u ( im, jm, km, 0. );							// auxilliar field u
	Array	aux_v ( im, jm, km, 0. );							// auxilliar field v
	Array	aux_w ( im, jm, km, 0. );							// auxilliar field w
	Array	Latency ( im, jm, km, 0. );						// latent heat
	Array	Rain ( im, jm, km, 0. );								// rain
	Array	Ice ( im, jm, km, 0. );								// ice deposit
	Array	Rain_super ( im, jm, km, 0. );					// Rain_super occurrence
	Array	IceLayer ( im, jm, km, 0. );						// ice shield
	Array	Condensation_3D ( im, jm, km, 0. );		// 3D condensation
	Array	Evaporation_3D ( im, jm, km, 0. );			// 3D Evaporation


/*
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] = 0.;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] = 0.;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] = 0.;
			}
		}
	}
*/



//	cout << " ***** printout of 3D-fields ***** " << endl;
//	v.printArray();
//	cout << " ***** printout of 2D-fields ***** " << endl;
//	Vegetation.printArray_2D();
//	cout << " ***** printout of 1D-fields ***** " << endl;
//	rad.printArray_1D();

	cout.precision ( 6 );
	cout.setf ( ios::fixed );


// 	initial values for the number of computed steps and the time
	n = 0;
	time = dt;
	velocity_iter = 1;
	pressure_iter = 1;
	velocity_iter_2D = 1;
	pressure_iter_2D = 1;
	switch_2D = 0;
	residuum = residuum_old = 0.;

// radial expansion of the computational field for the computation of initial values
	i_max = 32;			// corresponds to about 16 km above sea level, maximum hight of the tropopause at equator
	i_beg = 16;			// corresponds to about 8 km above sea level, maximum hight of the tropopause at poles

// time slice to start with is the modern world
	Ma = 0;
	i_time_slice = 0;


// choice of the time slice by Ma and by author
	string Name_Bathymetry_File;
	stringstream My;

// naming a file to read the surface temperature of the modern world
	string Name_SurfaceTemperature_File; 
	stringstream ssNameSurfaceTemperature;
	ssNameSurfaceTemperature << "SurfaceTemperature.xyz";
	Name_SurfaceTemperature_File = ssNameSurfaceTemperature.str();

// naming a file to read the surface precipitation by NASA
	string Name_SurfacePrecipitation_File; 
	stringstream ssNameSurfacePrecipitation;
	ssNameSurfacePrecipitation << "SurfacePrecipitation_NASA.xyz";
	Name_SurfacePrecipitation_File = ssNameSurfacePrecipitation.str();

// naming a possibly available sequel file
	string Name_Sequel_File;
	stringstream ssNameSequel;

// naming the output netCDF-file
	string Name_netCDF_File;
	stringstream ssNameNetCDF;


	cout << "\n\n\n\n";
	cout << "***** Atmosphere General Circulation Model ( AGCM ) applied to laminar flow" << endl;
	cout << "***** Program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
	cout << "***** Finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
	cout << "***** with 2 additional transport equations to describe the water vapour and co2 concentration" << endl;
	cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations" << endl << endl;

	cout << "***** original program name:  " << __FILE__ << endl;
	cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;


//	Coordinate system in form of a spherical shell
//	rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
	rad.Coordinates ( );
	the.Coordinates ( );
	phi.Coordinates ( );


// 	reading of the surface temperature file if available
	FILE *SurfaceTemperature;
	SurfaceTemperature = fopen ( Name_SurfaceTemperature_File.c_str(), "r" );

// 	if not available, prepare initial conditions, otherwise skip
	if ( SurfaceTemperature != NULL )
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: exists!" << endl;
	}
	else
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: could not be read!" << endl << endl;
		cout << endl << endl;
	}


// 	reading of the surface precipitation file if available
	FILE *SurfacePrecipitation;
	SurfacePrecipitation = fopen ( Name_SurfacePrecipitation_File.c_str(), "r" );

// 	if not available, prepare initial conditions, otherwise skip
	if ( SurfacePrecipitation != NULL )
	{
		cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: exists!" << endl;
	}
	else
	{
		cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: could not be read!" << endl << endl;
		cout << endl << endl;
		n++;
	}



// 	class calls for the following solution of the flow properties

//  class RHS_Atmosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
	RHS_Atmosphere		prepare ( im, jm, km, dt, dr, dthe, dphi, re, ec, sc_WaterVapour, sc_CO2, gr, pr, omega, coriolis, centrifugal, WaterVapour, buoyancy, CO2 );

//  class RungeKutta_Atmosphere for the explicit solution of the Navier-Stokes equations
	RungeKutta_Atmosphere		result ( im, jm, km, dt );

//  class Pressure for the subsequent computation of the pressure by a separat Euler equation
	Pressure		startPressure ( im, jm, km, dr, dthe, dphi );

// 	class BC_Thermo for the initial and boundary conditions of the flow properties
	BC_Thermo		circulation ( im, jm, km, t, c, aux_v, aux_w );

// 	class BC_Atmosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
	BC_Atmosphere		boundary ( im, jm, km );

//	class Restore to restore the iterational values from new to old
	Restore		oldnew( im, jm, km );

//	class Results_MSL_Atm to compute and show results on the mean sea level, MSL
	Results_MSL_Atm		calculate_MSL ( im, jm, km, sun, ep, hp, p_0, t_0, c_0, sigma, albedo, lv, cp_l, L_atm, dr, r_0_air, R_Air, r_0_water_vapour, R_WaterVapour, co2_vegetation, co2_ocean, co2_land );

//	class File_NetCDF to write results in the format of a netCDF-file
//	File_NetCDF		printoutNetCDF ( im, jm, km );





// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//	choice of the time slice to be computed
	time_slice_sequel:

// choice of the time slice by Ma and by author
	if ( Ma == 0 )
	{
		Name_Bathymetry_File = "0Ma_etopo.xyz";
		Name_Sequel_File = "[0Ma_etopo.xyz]_Sequel_Atm.seq";
		Name_netCDF_File = "[0Ma_etopo.xyz]_atmosphere.nc";
	}
	else 
	{
		n = 1;
		My << time_slice [ i_time_slice ] << "Ma_Golonka.xyz";
		Name_Bathymetry_File = My.str();
		My.str("");
		My.ignore(My.rdbuf()->in_avail());
	}


	ssNameSequel << "[" << Name_Bathymetry_File << "]_Sequel_Atm.seq";
	Name_Sequel_File = ssNameSequel.str();
	ssNameSequel.str("");
	ssNameSequel.ignore(ssNameSequel.rdbuf()->in_avail());

	ssNameNetCDF << Name_Bathymetry_File << "_atmosphere.nc";
	Name_netCDF_File = ssNameNetCDF.str();
	ssNameNetCDF.str("");
	ssNameNetCDF.ignore(ssNameNetCDF.rdbuf()->in_avail());


// 	reading of the sequel file if available
	FILE *Atmosphere_read;
	Atmosphere_read = fopen ( Name_Sequel_File.c_str(), "r" );

// 	if not available, skip
	if ( Atmosphere_read == NULL )
	{
		cout << "***** file ::::: " << Name_Sequel_File << " ::::: not yet exists!" << endl;
	}
	else
	{
		PostProcess_Atmosphere		read_File ( im, jm, km );
		read_File.Atmosphere_SequelFile_read ( Name_Bathymetry_File, n, time, rad, the, phi, h, t, u, v, w, c, co2, tn, un, vn, wn, cn, co2n, aux_u, aux_v, aux_w, t_j, c_j );
		cout << "***** file ::::: " << Name_Sequel_File << " ::::: could be read!" << endl << endl;
		cout << "***** Atmosphere_SequelFile_read in AGCM_main.cpp:   n = " << n << "  time = " << time << endl << endl;
		n++;
	}



// 	reading of the bathymetry file
	FILE *Atmosphere_Bathymetry;
	Atmosphere_Bathymetry = fopen ( Name_Bathymetry_File.c_str(), "r" );

	if ( Atmosphere_Bathymetry != NULL )
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: exists!" << endl << endl;
	}
	else
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: could not be read!" << endl << endl;
		n++;
	}


	cout << "***** time slice for the Atmospheric Global Circulation Modell ( AGCM ) is:    Ma = " << Ma << " million years" << endl << endl;
	cout << "***** bathymetry/topography given by the x-y-z data set:    " << Name_Bathymetry_File.c_str() << endl << endl;





// 	class BC_Bathymetry_Atmosphere for the geometrical boundary condition of the computational area
	BC_Bathymetry_Atmosphere		LandArea ( im, jm, km );

// 	class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
	LandArea.BC_MountainSurface ( Name_Bathymetry_File, L_atm, h, aux_w );

//  configuration of the initial and boundary conditions for the temperature, CO2 und water vapour on land and ocean surface

//  surface temperature from World Ocean Atlas 2009 given as boundary condition
	if ( Ma == 0 ) circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t_j, t );

//  surface precipitation from NASA for comparison
	if ( Ma == 0 ) circulation.BC_Surface_Precipitation ( Name_SurfacePrecipitation_File, precipitation_j );

//  surface pressure computed by surface temperature with gas equation
//	circulation.BC_Pressure ( r_0_air, R_Air, p_0, t_0, p_j, t_j, p, t, h );

//  parabolic temperature distribution from pol to pol, maximum temperature at equator
	circulation.BC_Temperature ( i_max, i_beg, Ma, Ma_max, Ma_max_half, sun_position_lat, sun_position_lon, declination, sun, ep, hp, t_0, p_0, t_land_plus, t_cretaceous_max, t_Average, co2_Average, t_equator, t_pole, t_tropopause, t_j, c_j, h, t, p );

//  computation of radiation at variable sun position
	if ( sun == 1 ) LandArea.BC_Radiation ( t_0, ik, sigma, albedo_extra, epsilon_atmos, Precipitation, Ik, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, t_j, t );

//	computation of ice shield following the theorie by Milankowitsch
	if ( IceShield == 1 ) LandArea.BC_IceShield ( Ma, t_0, h, t, c, IceLayer, Ice_Balance, Ice_Balance_add );

//  parabolic water vapour distribution from pol to pol, maximum water vapour volume at equator
	circulation.BC_WaterVapour ( i_max, i_beg, ep, hp, t_0, c_0, p_0, c_land_minus, c_ocean_minus, c_tropopause, c_j, h, t, p, c, t_j );

//  parabolic CO2 distribution from pol to pol, maximum CO2 volume at equator
	circulation.BC_CO2 ( i_max, i_beg, co2_0, co2_Average, co2_equator, co2_pole, co2_tropopause, co2_vegetation, co2_ocean, co2_land, co2_j, Vegetation, h, t, p, co2 );

// 	initial conditions for u-v-w-velocity components
	if ( Ma == 0 ) circulation.IC_CellStructure ( i_max, i_beg, u, v, w );

//	initial conditions for v and w velocity components at the sea surface close to east or west coasts, to close gyres
//	if ( Ma == 0 ) circulation.IC_v_w_WestEastCoast ( h, u, v, w );
//	circulation.IC_v_w_WestEastCoast ( h, u, v, w );

//	initial conditions for v and w velocity components at the sea surface close to east or west coasts, to close gyres
//	if ( Ma == 0 ) oceanflow.IC_v_w_Smoothing ( h, v, w, t, c );
//	circulation.IC_v_w_Smoothing ( velocity_iter, h, u, v, w, t, c );


//	storing of velocity components, pressure and temperature for iteration start
	oldnew.restoreOldNew ( .9, u, v, w, t, p, c, co2, un, vn, wn, tn, pn, cn, co2n );
	oldnew.restoreOldNew_2D ( .9, v, w, p, vn, wn, pn );

// computation of the ratio ocean to land areas, also supply and removal of CO2 on land, ocean and by vegetation
	calculate_MSL.land_oceanFraction ( h );

//		goto Print_commands;														// only initial conditions wihtout iterations


// ******************************************   start of pressure and velocity iterations ************************************************************************


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of pressure loop : if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::

Pressure_loop:

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of velocity loop: while ( min >= epsres )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	min = epsres * 20.;
	min_u = min_v = min_w = min_t = min_c = min_p = epsres * 3.;
	velocity_iter = 0;
	velocity_iter_2D = 0;


//	query to realize zero divergence of the continuity equation ( div c = 0 )
	while ( min >= epsres )
	{
//		limit of the computation in the sense of time steps
		n++;
		if ( n > nm )
		{
			cout << endl; 
			cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
			cout << endl;
			break;
		}

//		limit of maximum number of iterations ( velocity_iter_max )
		velocity_iter++;
		if ( velocity_iter > velocity_iter_max )
		{
			n--;
			velocity_iter--;
			break;
		}

		 if ( switch_2D == 1 ) goto process_3D;


// **********************************   start of pressure and velocity iterations for the 2D iterational process   *********************************


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of pressure loop_2D : if ( pressure_iter_2D > pressure_iter_max_2D )   :::::::::::::::::::::::::::::::::::::::::::

Pressure_loop_2D:

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of velocity loop_2D: while ( min >= epsres )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		while ( velocity_iter_2D <= velocity_iter_max_2D )
		{

//		limit of maximum number of iterations ( velocity_iter_max_2D )
			velocity_iter_2D++;
			if ( velocity_iter_2D > velocity_iter_max_2D )
			{
				goto Pressure_iteration_2D;
			}

//		class BC_Atmosphaere for the geometry of a shell of a sphere
			boundary.BC_theta ( t_tropopause, c_tropopause, t, u, v, w, p, c, co2, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Rain, Ice );
			boundary.BC_phi ( t, u, v, w, p, c, co2, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Rain, Ice );

// 		old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy		min_Residuum_old_2D ( n, im, jm, km, dr, dthe, dphi );
			min_Residuum_old_2D.residuumQuery_2D ( j_res, k_res, min, rad, the, v, w );
			residuum_old = min;

//		class RungeKutta for the solution of the differential equations describing the flow properties
			result.solveRungeKutta_2D_Atmosphere ( prepare, rad, the, phi, rhs_v, rhs_w, h, v, w, p, vn, wn, aux_v, aux_w );

//		state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
			Accuracy		min_Stationary_2D ( n, im, jm, km, dr, dthe, dphi );
			min_Stationary_2D.steadyQuery_2D ( j_v, k_v, j_w, k_w, j_p, k_p, min_v, min_w, min_p, v, vn, w, wn, p, pn );

// 		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy		min_Residuum_2D ( n, im, jm, km, dr, dthe, dphi );
			min_Residuum_2D.residuumQuery_2D ( j_res, k_res, min, rad, the, v, w );
			residuum = min;
			min = fabs ( ( residuum - residuum_old ) / residuum_old );

//		statements on the convergence und iterational process
			Accuracy		printout_2D ( im, Ma, n, velocity_iter_2D, pressure_iter_2D, min, L_atm );
			printout_2D.iterationPrintout_2D ( nm, velocity_iter_max_2D, pressure_iter_max_2D, j_res, k_res, j_v, k_v, j_w, k_w, j_p, k_p, min, min_v, min_w, min_p );

			oldnew.restoreOldNew_2D ( 1., v, w, p, vn, wn, pn );

		}



//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop_2D: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Pressure_iteration_2D:

//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
	startPressure.computePressure_2D ( pa, rad, the, p, h, rhs_v, rhs_w, aux_v, aux_w );

//	statements on the convergence und iterational process
	pressure_iter_2D++;
	velocity_iter_2D = 0;

	if ( pressure_iter_2D >= pressure_iter_max_2D + 1 ) 
	{
		switch_2D = 1;
		goto process_3D;
	}
	else goto Pressure_loop_2D;


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of pressure loop_2D: if ( pressure_iter_2D > pressure_iter_max_2D )   :::::::::::::::::::::::::::::::::::::::::::




	process_3D:

		if ( min >= epsres )
		{
			time = time + dt;
		}

// 		old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
		Accuracy		min_Residuum_old ( n, im, jm, km, dr, dthe, dphi );
		min_Residuum_old.residuumQuery ( i_res, j_res, k_res, min, rad, the, u, v, w );
		residuum_old = min;

// 		class BC_Atmosphaere for the geometry of a shell of a sphere
		boundary.BC_radius ( tao, tau, pa, ca, co2a, dr, t_tropopause, c_tropopause, rad, co2_vegetation, co2_ocean, co2_land, Vegetation, h, t, u, v, w, p, c, co2, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Rain, Ice );
		boundary.BC_theta ( t_tropopause, c_tropopause, t, u, v, w, p, c, co2, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Rain, Ice );
		boundary.BC_phi    ( t, u, v, w, p, c, co2, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Rain, Ice );
//		boundary.BC_NST_control_3D ( dr, dthe, dphi, re, mue_air, mue_water, h, u, v, w, t, p, c, co2, aux_u, aux_v, aux_w, rad, the );

// initial condition for u-v-w-velocity components on west/east coasts
//		if ( n <= 3 )    circulation.IC_v_w_WestEastCoast ( h, u, v, w );
//		circulation.IC_v_w_WestEastCoast ( h, u, v, w );


// 		class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
		LandArea.BC_SolidGround ( Ma, hp, ep, c_land_minus, co2_vegetation, t_0, tau, p_0, pa, h, t, u, v, w, p, c, co2, tn, un, vn, wn, pn, cn, co2n, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, rhs_co2, t_j, c_j, co2_j, Vegetation );

// 		class RungeKutta for the solution of the differential equations describing the flow properties
		result.solveRungeKutta_Atmosphere ( prepare, lv, ls, ep, hp, u_0, t_0, t_Boussinesq, c_0, co2_0, p_0, r_0_air, r_0_water_vapour, r_0_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_co2, h, t, u, v, w, p, c, co2, tn, un, vn, wn, cn, co2n, aux_u, aux_v, aux_w, Latency, Rain, Ice, Rain_super, IceLayer );

// 		class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
		LandArea.BC_SolidGround ( Ma, hp, ep, c_land_minus, co2_vegetation, t_0, tau, p_0, pa, h, t, u, v, w, p, c, co2, tn, un, vn, wn, pn, cn, co2n, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, rhs_co2, t_j, c_j, co2_j, Vegetation );

//		state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
		Accuracy		min_Stationary ( n, im, jm, km, dr, dthe, dphi );
		min_Stationary.steadyQuery ( i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_co2, j_co2, k_co2, i_p, j_p, k_p, min_u, min_v, min_w, min_t, min_c, min_co2, min_p, u, un, v, vn, w, wn, t, tn, c, cn, co2, co2n, p, pn );


// 		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
		Accuracy		min_Residuum ( n, im, jm, km, dr, dthe, dphi );
		min_Residuum.residuumQuery ( i_res, j_res, k_res, min, rad, the, u, v, w );
		residuum = min;
		min = fabs ( ( residuum - residuum_old ) / residuum_old );


//		statements on the convergence und iterational process
		Accuracy		printout ( im, Ma, n, velocity_iter, pressure_iter, min, L_atm );
		printout.iterationPrintout ( nm, velocity_iter_max, pressure_iter_max, i_res, j_res, k_res, i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_co2, j_co2, k_co2, i_p, j_p, k_p, min_u, min_v, min_w, min_t, min_c, min_co2, min_p );


//		searching of maximum and minimum values of temperature
		string str_max_temperature = " max temperature ", str_min_temperature = " min temperature ", str_unit_temperature = "°C";
		MinMax		minmaxTemperature ( im, jm, km );
		minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

//		searching of maximum and minimum values of water vapour
		string str_max_water_vapour = " max water vapour ", str_min_water_vapour = " min water vapour ", str_unit_water_vapour = "g/kg";
		MinMax		minmaxWaterVapour ( im, jm, km );
		minmaxWaterVapour.searchMinMax_3D ( str_max_water_vapour, str_min_water_vapour, str_unit_water_vapour, c, h );

//		searching of maximum and minimum values of co2
		string str_max_co2 = " max co2 ", str_min_co2 = " min co2 ", str_unit_co2 = "ppm";
		MinMax		minmaxCO2 ( im, jm, km );
		minmaxCO2.searchMinMax_3D ( str_max_co2, str_min_co2, str_unit_co2, co2, h );


//		total precipitation as the sum on rain and snow along a normally extended virtual column
//		in r-direction water vapour volume above the saturation vapour pressure is added
		calculate_MSL.run_MSL_data ( h, c, t, p, u, Rain, Rain_super, Ice, Latency, Condensation_3D, Evaporation_3D, Precipitation, Water, Water_super, IceAir, Evaporation, Condensation, precipitable_water, Q_Balance_Radiation, Q_Evaporation, Q_latent, Q_sensible, Q_diff, Evaporation_Penman, Evaporation_Haude, t_j, c_j );


//		searching of maximum and minimum values of precipitation
		string str_max_precipitation = " max precipitation ", str_min_precipitation = " min precipitation ", str_unit_precipitation = "mm";
		MinMax		minmaxPrecipitation ( jm, km, coeff_mmWS );
		minmaxPrecipitation.searchMinMax ( str_max_precipitation, str_min_precipitation, str_unit_precipitation, Precipitation, h );
		max_Precipitation = minmaxPrecipitation.out_maxValue (  );


//		description of vegetation areas following the local dimensionles precipitation values, maximum value now is 1
		calculate_MSL.vegetationDistribution ( max_Precipitation, Precipitation, Vegetation, t, h );


//		searching of maximum and minimum values of precipitable water
		string str_max_precipitable_water = " max precipitable water ", str_min_precipitable_water = " min precipitable water ", str_unit_precipitable_water = "mm";
		MinMax		minmaxPrecipitable_water ( jm, km, coeff_mmWS );
		minmaxPrecipitable_water.searchMinMax ( str_max_precipitable_water, str_min_precipitable_water, str_unit_precipitable_water, precipitable_water, h );

//		searching of maximum and minimum values of radiation balance
		string str_max_Q_Balance_Radiation = " max Q Balance Radiation ", str_min_Q_Balance_Radiation = " min Q Balance Radiation ", str_unit_Q_Balance_Radiation = " W/m2";
		MinMax		minmaxQ_Balance_Radiation ( jm, km, coeff_mmWS );
		minmaxQ_Balance_Radiation.searchMinMax ( str_max_Q_Balance_Radiation, str_min_Q_Balance_Radiation, str_unit_Q_Balance_Radiation, Q_Balance_Radiation, h );

//		searching of maximum and minimum values of latent energy
		string str_max_Q_latent = " max Q latent ", str_min_Q_latent = " min Q latent ", str_unit_Q_latent = " W/m2";
		MinMax		minmaxQ_latent ( jm, km, coeff_mmWS );
		minmaxQ_latent.searchMinMax ( str_max_Q_latent, str_min_Q_latent, str_unit_Q_latent, Q_latent, h );

//		searching of maximum and minimum values of sensible energy
		string str_max_Q_sensible = " max Q sensible ", str_min_Q_sensible = " min Q sensible ", str_unit_Q_sensible = " W/m2";
		MinMax		minmaxQ_sensible ( jm, km, coeff_mmWS );
		minmaxQ_sensible.searchMinMax ( str_max_Q_sensible, str_min_Q_sensible, str_unit_Q_sensible, Q_sensible, h );

//		searching of maximum and minimum values of difference in energy
		string str_max_bottom_heat = " max bottom heat ", str_min_bottom_heat = " min bottom heat ", str_unit_bottom_heat = " W/m2";
		MinMax		minmaxQ_diff ( jm, km, coeff_mmWS );
		minmaxQ_diff.searchMinMax ( str_max_bottom_heat, str_min_bottom_heat, str_unit_bottom_heat, Q_diff, h );

//		searching of maximum and minimum values of evaporation
		string str_max_heat_evaporation = " max heat evaporation ", str_min_heat_evaporation = " min heat evaporation ", str_unit_heat_evaporation = " W/m2";
		MinMax		minmaxQ_Evaporation ( jm, km, coeff_mmWS );
		minmaxQ_Evaporation.searchMinMax ( str_max_heat_evaporation, str_min_heat_evaporation, str_unit_heat_evaporation, Q_Evaporation, h );

//		searching of maximum and minimum values of evaporation by Haude
		string str_max_evaporation_Haude = " max evaporation Haude ", str_min_evaporation_Haude = " min evaporation Haude ", str_unit_evaporation_Haude = "mm/d";
		MinMax		minmaxEvaporation_Haude ( jm, km, coeff_mmWS );
		minmaxEvaporation_Haude.searchMinMax ( str_max_evaporation_Haude, str_min_evaporation_Haude, str_unit_evaporation_Haude, Evaporation_Haude, h );

//		searching of maximum and minimum values of evaporation by Penman
		string str_max_evaporation_Penman = " max evaporation Penman ", str_min_evaporation_Penman = " min evaporation Penman ", str_unit_evaporation_Penman = "mm/d";
		MinMax		minmaxEvaporation_Penman ( jm, km, coeff_mmWS );
		minmaxEvaporation_Penman.searchMinMax ( str_max_evaporation_Penman, str_min_evaporation_Penman, str_unit_evaporation_Penman, Evaporation_Penman, h );


// printout of results at certain positions
		calculate_MSL.show_MSL_data ( h, c, t, p, u, Rain, Ice, Latency, Precipitation, IceAir, Evaporation, Condensation, precipitable_water, Q_Balance_Radiation, Q_Evaporation, Q_latent, Q_sensible, Q_diff, Evaporation_Penman, Evaporation_Haude );


//	restoring the velocity component and the temperature for the new time step
		oldnew.restoreOldNew ( 1., u, v, w, t, p, c, co2, un, vn, wn, tn, pn, cn, co2n );
	}


//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
	startPressure.computePressure ( pa, rad, the, p, h, rhs_u, rhs_v, rhs_w, aux_u, aux_v, aux_w );

//	statements on the convergence und iterational process
	pressure_iter++;
	velocity_iter = 0;

	if ( pressure_iter >= pressure_iter_max + 1 ) 
	{
		goto Print_commands;
	}
	else 
	{
		velocity_iter_2D = velocity_iter_max_2D;
		pressure_iter_2D = pressure_iter_max_2D;
		switch_2D = 1;
		goto Pressure_loop;
	}

                                            //  end of loop: 			if ( pressure_iter > pressure_iter_max )


Print_commands:

//	class Print_Atmosphaere for the printout of results
	Print_Atmosphere		printout ( im, jm, km, nm, n, time );

//	printout in ParaView files, netCDF files and sequel files

	pressure_iter_aux = pressure_iter - 1;

//	results written in netCDF format
//	printoutNetCDF.aus_NetCDF( Name_netCDF_File, v, w, h, Precipitation, precipitable_water );

//	class PostProcess_Atmosphaere for the printing of results
	PostProcess_Atmosphere		write_File ( im, jm, km );

//	writing of data in ParaView files
//	londitudinal data along constant latitudes
	j_longal = 75;
	write_File.paraview_vtk_longal ( Name_Bathymetry_File, j_longal, pressure_iter_aux, h, p, t, u, v, w, c, co2, aux_u, aux_v, aux_w, Latency, Rain, Ice, Rain_super, IceLayer );

	k_zonal = 185;
	write_File.paraview_vtk_zonal ( Name_Bathymetry_File, k_zonal, pressure_iter_aux, h, p, t, u, v, w, c, co2, aux_u, aux_v, aux_w, Latency, Rain, Ice, Rain_super, Condensation_3D, Evaporation_3D );

//	radial data along constant hight above ground
	i_radial = 0;
	write_File.paraview_vtk_radial ( Name_Bathymetry_File, i_radial, pressure_iter_aux, h, p, t, u, v, w, c, co2, aux_u, aux_v, aux_w, Latency, Rain, Ice, Rain_super, IceLayer, Precipitation, Evaporation, IceAir, Condensation, precipitable_water, Q_diff, Q_Balance_Radiation, Q_latent, Q_sensible, Evaporation_Penman, Evaporation_Haude, Q_Evaporation, precipitation_j, Water_super, Water );

//	3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
	write_File.paraview_panorama_vts ( Name_Bathymetry_File, pressure_iter_aux, h, t, p, u, v, w, c, co2, aux_u, aux_v, aux_w, Latency, Rain, Ice, Rain_super, IceLayer );

//	3-dimensional data in spherical coordinate system for a streamline pattern in a shell of a sphere
//	write_File.paraview_vts ( Name_Bathymetry_File, n, rad, the, phi, h, t, p, u, v, w, c, co2, aux_u, aux_v, aux_w, Latency, Rain, Ice, Rain_super, IceLayer );


//	writing of sequential data for the sequel file
	if ( SequelFile == 1 )
	{
		write_File.Atmosphere_SequelFile_write ( Name_Bathymetry_File, n, time, rad, the, phi, h, t, u, v, w, c, co2, tn, un, vn, wn, cn, co2n, aux_u, aux_v, aux_w, t_j, c_j  );
	}

//	writing of v-w-data in the v_w_transfer file
		PostProcess_Atmosphere		write_v_w_Transfer_File ( im, jm, km );
		write_v_w_Transfer_File.Atmosphere_v_w_Transfer ( Name_Bathymetry_File, v, w, p );


//	writing of plot data in the PlotData file
		PostProcess_Atmosphere		write_PlotData_File ( im, jm, km );
		write_PlotData_File.Atmosphere_PlotData ( Name_Bathymetry_File, u_0, t_0, v, w, t, c, Precipitation, precipitable_water );




//	statements on the convergence und iterational process
	velocity_iter = 0;
	n--;

	if ( pressure_iter > pressure_iter_max ) goto time_slice_change;
	else goto Pressure_loop;

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of pressure loop: if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::



	time_slice_change:

//	choice of the next time slice after nm iterations reached
	time = dt;
	velocity_iter = 1;
	pressure_iter = 1;
	velocity_iter_2D = 1;
	pressure_iter_2D = 1;
	switch_2D = 0;

	i_time_slice++;
	Ma = time_slice [ i_time_slice ];
	if ( Ma > 0 ) goto finish;
	if ( i_time_slice >= i_time_slice_max ) goto finish;
	else goto time_slice_sequel;

//   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	finish:

// delete temporary arrays
	delete [ ] time_slice;

// 	final remarks
	cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;

	if ( velocity_iter == velocity_iter_max )	cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of maximum time steps ***** \n\n" << endl;

	if ( min <= epsres )		cout << "***** steady solution reached! *****" << endl;

	cout << endl;
	cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
	cout << "\n\n\n\n";

	return 0;
}
