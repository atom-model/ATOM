/*
 * Atmosphere General Circulation Model ( AGCM ) applied to laminar flow
 * program for the computation of geo-atmospherical circulating flows in a spherical shell
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 4 additional transport equations to describe the water vapour, cloud water, cloud ice and co2 concentration
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * multi-layer and two-layer radiation model for the computation of the surface temperature
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * water vapour distribution given by Clausius-Claperon equation for the partial pressure
 * water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model
 * two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code ( German Weather Forecast )
 * rain and snow precipitation solved by column equilibrium applying the diagnostic equations
 * co2 concentration appears in the absorptivity of the radiation models
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <new>
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
#include "File_NetCDF_Atm.h"

// #include "memcheck.h"
// #include "valgrind.h"


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
// "CO2 1"						with inclusion of co2
// "CO2 0"						without inclusion of co2

// "buoyancy"				Computation with and without a buoyancy term following the Boussinesq-Approximation
// "buoyancy 1"				with inclusion of a buoyancy term
// "buoyancy 0"				without inclusion of a buoyancy term

// "sun"							Computation with and without an averaged sun radiation
// "sun 1"						with inclusion of radiation
// "sun 0"						without inclusion of radiation

// "IceShield"					Computation with and without ice shields
// "IceShield 1"				with inclusion of ice shields
// "IceShield 0"				without inclusion of ice shields

// "RadiationModel"		Computation of temperature with radiation model
// "RadiationModel 1"		with a constant parabolic temperature distribution
// "RadiationModel 2"		with a two-layer radiation model
// "RadiationModel 3"		with a multi-layer radiation model

// Earth's radius is r_earth = 6731 km, equatorial circumference of the earth 40074 km
// 40074 km / 360 deg = 111.32 km
// 40074 km / 180 deg = 55.66 km

// for 20 km expansion of the area of circulation compares to 0.04 [ / ] with 40 steps of size 0.001 
// 20000 m / 40 steps = 500 m

// Definition of meridional and longitudinal step sizes 
// for example: dthe = the_degree / pi180 = 1. / 57.3 = 0.01745

// velocities given in m/s, 1 m/s compares to 3.6 km/h
// maximum velocity in the subtropical jet  w_max = 30 m/s compares to 108 km/h as annual mean 

// minimum temperature at the tropopause t_min = .789 compares to -60° C compares to 213 K
// temperature t_0 = 1.0 compares to 0° C compares to 273,15 K
// temperature t_0 = 0.9561 compares to -12° C compares to 261,15 K
// temperature difference from equator to pole   18°C compares to  t_delta = 0.0665  compares to  18 K

// water vapour c = 35 g/kg, means g water vapour per kg dry air
// density of humid air is r = r_air * ( 1 + c ) / ( 1 + R_WaterVapour / R_Air * c )

// rate of CO2 compares to 600Gt in the atmosphere ( 0.02% ), 100Gt consumed by vegetation ( 0.00333% ), ratio 1/6 = 0.166667 to be checked!





int main ( int argc, char *argv[ ] )
{
	bool verbose = false;

// maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km )
// maximum number of overall iterations ( n )
// maximum number of inner velocity loop iterations ( velocity_iter_max )
// maximum number of outer pressure loop iterations ( pressure_iter_max )

	int im = 41, jm = 181, km = 361, nm = 200, velocity_iter_max = 2, pressure_iter_max = 2;
	int velocity_iter_max_2D = 2, pressure_iter_max_2D = 2;

	int n = 0, i_radial =0, j_longal = 0, k_zonal = 0, i_max = 0, i_beg = 0;
	int velocity_iter = 0, pressure_iter = 0, pressure_iter_aux = 0, velocity_iter_2D = 0, pressure_iter_2D = 0;
//	int velocity_iter, pressure_iter, pressure_iter_aux;
	int i_res = 0, j_res = 0, k_res = 0;
	int Ma = 0, i_time_slice = 0, i_time_slice_max = 0;
	int switch_2D = 0;

	double time = 0;
	double residuum = 0, residuum_old = 0, min = 0;
//	double max_Precipitation, max_precipitable_water, max_CO2_total;
	double max_Precipitation = 0;

	int SequelFile = 0;										// sequel file will not be written

	double coriolis = 1.;									// computation with Coriolis force
	double centrifugal = 1.;							// computation with centrifugal force
	double WaterVapour = 1.;						// computation with water vapour
	double buoyancy = 1.;								// computation with buoyancy
	double CO2 = 1.;										// computation with CO2

	double epsres = 0.00001;							// accuracy of relative and absolute errors

	int  sun = 0;												// while no variable sun position wanted
	int RadiationModel = 3;								// surface temperature computation by a radiation model
	int  IceShield = 0;										// while no ice shields wanted

	int declination = 0;									// position of sun axis, today 23,4°, 21.12.: -23,4°, am 21.3. und 23.9.: 0°, 21.6.: +23,4°, in between sin form
	int sun_position_lat = 60;							// position of sun j = 120 means 30°S, j = 60 means 30°N
	int sun_position_lon = 180;						// position of sun k = 180 means 0° or 180° E ( Greenwich, zero meridian )

	int Ma_max = 300;									// parabolic temperature distribution 300 Ma (From Ruddiman)
	int Ma_max_half = 150;								// half of time scale

	double pi180 = 180./M_PI;						// pi180 = 57.3

	double L_atm = 20000.;							// extension of the atmosphere shell in m, 20000 m / 40 steps = 500 m
	double dt = 0.0001;									// time step coincides with the CFL condition

	double dr = 0.0005;									// compares to 500 m hight, 0.0005 * 40 = .02 * 1000 km = 20 km
	double the_degree = 1.;							// compares to 1° step size laterally
	double phi_degree = 1.;							// compares to 1° step size longitudinally
	double dthe = the_degree / pi180;			//dthe = the_degree / pi180 = 1.0 / 57.3 = 0.01745, 180 * .01745 = 3.141
	double dphi = phi_degree / pi180;				//dphi = phi_degree / pi180 = 1.0 / 57.3 = 0.01745, 360 * .01745 = 6.282

	double the0 = 0.;										// North Pole
	double phi0 = 0.;										// zero meridian in Greenwich
	double r0 = 6.731; 									// earth's radius is r_earth = 6731 km compares to 6.731 [ / ] * 1000 km, circumference of the earth 40074 km

	double ik = 1366.;									// solar constant in W/m2
	double sigma = 5.670280e-8;					// Stefan-Boltzmann constant W/( m²*K4 )
	double albedo_extra = .15;						// capability of reflection of short wave radiation, global albedo_extra extraterrestric
	double epsilon_extra = .71;						// capability of emissions in the atmosphere

	double re = 1000.;									// Reynolds number: ratio viscous to inertia forces, Re = u * L / nue
	double ec = .00044;									// Eckert number: ratio kinetic energy to enthalpy, Ec = u² / cp T
	double sc_WaterVapour = .6;					// Schmidt number of water vapour, Sc = nue / D
	double sc_CO2 = .96;								// Schmidt number of CO2
	double pr = .7179;									// Prandtl number of air for laminar flows
	double g = 9.8066;									// gravitational acceleration of the earth
	double omega = 7.29e-5;							// rotation number of the earth
	double ep = .623;									// ratio of the gas constants of dry air to water vapour [ / ]
	double hp = 6.1078;									// water vapour pressure at T = 0°C: E = 6.1 hPa
	double R_Air = 287.1;								// specific gas constant of air in J/( kg*K ))
	double R_WaterVapour = 461.6;				// specific gas constant of water vapour in J/( kg*K ))
	double R_co2 = 188.91;							// specific gas constant of CO2 in J/( kg*4.5K ))
	double lv = 2.52e6;									// specific latent Evaporation heat ( Condensation heat ) in J/kg
	double ls = 2.83e6;								 	// specific latent vaporisation heat ( sublimation heat ) in J/kg
	double cp_l = 1005.;								// specific heat capacity of dry air at constant pressure and 20°C in J/( kg K )
	double lambda = .0262;							// heat transfer coefficient of air in W/m² K )
	double r_air = 1.2041;							// density of dry air in kg/m³ at 20°C
	double r_water = 1000.;						// density of water in kg/m³ at 20°C
	double r_water_vapour = 0.0094;			// density of saturated water vapour in kg/m³ at 10°C
	double r_co2 = 0.0019767;					// density of CO2 in kg/m³ at 25°C
	double gam = .65;									// constant slope of temperature 	gam = 0.65 K/100 m

	double u_0 = 15.;										// maximum value of velocity in 15 m/s compares to 54 km/h
	double p_0 = 1013.25;								// pressure at sea level in hPa
	double t_0 = 273.15;								// temperature in K compare to 0°C
	double c_0 = .035;									// maximum value of water vapour in kg / kg
	double co2_0 = 280.;								// maximum value of CO2 in ppm at preindustrial times

	double ua = 0.;										// initial velocity component in r-direction
	double va = 0.;										// initial velocity component in theta-direction
	double wa = 0.;										// initial velocity component in phi-direction
	double pa = 0.;										// initial value for the pressure field
	double ca = 0.;										// value 1.0 stands for the maximum value of 35 g/kg water vapour
	double ta = 1.;										// initial value for the temperature field, 1.0 compares to 0° C compares to 273.15 K
	double coa = 1.;									// initial value of co2 = 1.0 compares to 280 ppm in preindustrial times

	double t_cretaceous_max = 10.;				// maximum add of mean temperature in °C during cretaceous times
	double t_cretaceous = 0.;							// value at modern times
	double coeff_mmWS = r_air / r_water_vapour;	// coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]

	double radiation_ocean = 40.;					// increase of radiation at equator in W/m²
	double radiation_pole = - 40.;					// negative amount of radiation at poles in W/m²
	double radiation_equator = 100.;				// positive amount of radiation at equator in W/m²

	double t_average = 15.;							// mean temperature of the modern earth
	double t_equator = 1.1103;						// temperature t_0 = 1.1103 compares to 30.13° C compares to 303.28 K
	double t_pole = .8;									// temperature at the poles t_pole = 0.8 compares to -54.63°C compares to 218.52 K
	double t_tropopause = .78;						// temperature in the tropopause, t = 0.78 compares to -60.093°C compares to 213,057 K
//	double t_land = .003661;							// temperature increase on land by 1°C ( 1°C compares to t_land = 0.003661 )
	double t_land = .018305;							// temperature increase on land by 5°C ( 5°C compares to t_land = 0.018305 )

	double c_tropopause = 0.;						// minimum water vapour at tropopause c_tropopause = 0.01429 compares to 0.05 g/kg
	double c_land = .4;									// water vapour reduction on land ( 90% of the saturation value )
	double c_ocean = .5;								// water vapour reduction on sea surface ( 100% of the saturation value )

	double co2_average = 372.;						// rate of CO2 at preindustrial times
	double co2_equator = 330.;						// maximum rate of CO2 at sea level at equator, 1. compares to 330 ppm
	double co2_tropopause = 0.;					// minimum rate CO2 at tropopause 0 ppm
	double co2_pole = 305.;							// maximum rate of CO2 of the sea surface at poles
	double co2_cretaceous = 0.;						// value at modern times
	double co2_vegetation = 3.;						// value compares to 100/600Gt per year on the global surface by vegetation
	double co2_ocean = 0.;							// value compares to 0.6/600Gt per year on the sea surface
	double co2_land = 3.;								// value compares to 0.2/600Gt per year on land

	int *im_tropopause = new int [ jm ];			// location of the tropopause

// time slices to be run after actualizing 
	i_time_slice_max = 15;
	int *time_slice = new int [ i_time_slice_max ]; 	// time slices in Ma


	time_slice [ 0 ] = 0;									// Golonka Bathymetry and Topography
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
	Array_1D	rad;																			// radial coordinate direction
	rad.initArray_1D ( im, 1. );
	Array_1D	the;																			// lateral coordinate direction
	the.initArray_1D ( jm, 2. );

	Array_1D	phi;																			// longitudinal coordinate direction
	phi.initArray_1D ( km, 3. );


// 2D arrays
	Array_2D	Vegetation;																// vegetation via precipitation
	Vegetation.initArray_2D ( jm, km, 0. );

	Array_2D	LatentHeat;																// areas of higher latent heat
	LatentHeat.initArray_2D ( jm, km, 0. );
	Array_2D	Condensation;															// areas of higher condensation
	Condensation.initArray_2D ( jm, km, 0. );
	Array_2D	Evaporation;															// areas of higher evaporation
	Evaporation.initArray_2D ( jm, km, 0. );

	Array_2D	Precipitation;															// areas of higher precipitation
	Precipitation.initArray_2D ( jm, km, 0. );
	Array_2D	precipitable_water;													// areas of precipitable water in the air
	precipitable_water.initArray_2D ( jm, km, 0. );
	Array_2D	precipitation_NASA;												// surface precipitation from NASA
	precipitation_NASA.initArray_2D ( jm, km, 0. );

	Array_2D	Ice_Balance;															// rate of the ice shield
	Ice_Balance.initArray_2D ( jm, km, 0. );
	Array_2D	Ice_Balance_add;													// addition of all ice layers
	Ice_Balance_add.initArray_2D ( jm, km, 0. );

	Array_2D	Ik;																			// direct sun radiation, short wave
	Ik.initArray_2D ( jm, km, 0. );
	Array_2D	Radiation_Balance;													// radiation balance at the surface
	Radiation_Balance.initArray_2D ( jm, km, 0. );
	Array_2D	Radiation_Balance_par;											// radiation balance at the surface parabolic distribution
	Radiation_Balance_par.initArray_2D ( jm, km, 0. );
	Array_2D	Radiation_Balance_atm;										// radiation balance of the atmosphere
	Radiation_Balance_atm.initArray_2D ( jm, km, 0. );
	Array_2D	Radiation_Balance_bot;											// radiation balance of the ground
	Radiation_Balance_bot.initArray_2D ( jm, km, 0. );

	Array_2D	t_j;																			// 2D temperature at the surface
	t_j.initArray_2D ( jm, km, 0. );
	Array_2D	temp_eff_atm;														// effektive temperature in the atmosphere
	temp_eff_atm.initArray_2D ( jm, km, 0. );
	Array_2D	temp_eff_bot;															// effektive temperature on the ground
	temp_eff_bot.initArray_2D ( jm, km, 0. );
	Array_2D	temp_rad;																// effektive temperature based on radiation
	temp_rad.initArray_2D ( jm, km, 0. );

	Array_2D	albedo;																	// albedo = reflectivity
	albedo.initArray_2D ( jm, km, 0. );
	Array_2D	epsilon;																	// epsilon = absorptivity
	epsilon.initArray_2D ( jm, km, 0. );

	Array_2D	Q_Radiation;															// heat from the radiation balance in [W/m2]
	Q_Radiation.initArray_2D ( jm, km, 0. );
	Array_2D	Q_Evaporation;														// evaporation heat of water by Kuttler
	Q_Evaporation.initArray_2D ( jm, km, 0. );
	Array_2D	Q_latent;																	// latent heat from bottom values by the energy transport equation
	Q_latent.initArray_2D ( jm, km, 0. );
	Array_2D	Q_sensible;																// sensible heat from bottom values by the energy transport equation
	Q_sensible.initArray_2D ( jm, km, 0. );
	Array_2D	Q_bottom;																// difference by Q_Radiation - Q_latent - Q_sensible
	Q_bottom.initArray_2D ( jm, km, 0. );

	Array_2D	Evaporation_Haude;												// evaporation by Haude in [mm/d]
	Evaporation_Haude.initArray_2D ( jm, km, 0. );
	Array_2D	Evaporation_Penman;											// evaporation by Penman in [mm/d]
	Evaporation_Penman.initArray_2D ( jm, km, 0. );

	Array_2D	MaxCloud;																// maximum cloud water mass in a vertical column
	MaxCloud.initArray_2D ( jm, km, 0. );
	Array_2D	MaxIce;																	// maximum cloud ice mass in a vertical column
	MaxIce.initArray_2D ( jm, km, 0. );

	Array_2D	co2_total;																// areas of higher co2 concentration
	co2_total.initArray_2D ( jm, km, 0. );

	Array_2D	aux_2D_v;																// auxilliar field v
	aux_2D_v.initArray_2D ( jm, km, 0. );
	Array_2D	aux_2D_w;																// auxilliar field w
	aux_2D_w.initArray_2D ( jm, km, 0. );


// 3D arrays
	Array	h;																					// bathymetry, depth from sea level
	h.initArray ( im, jm, km, 0. );

	Array	t;																					// temperature
	t.initArray ( im, jm, km, ta );
	Array	u;																					// u-component velocity component in r-direction
	u.initArray ( im, jm, km, ua );
	Array	v;																					// v-component velocity component in theta-direction
	v.initArray ( im, jm, km, va );
	Array	w;																				// w-component velocity component in phi-direction
	w.initArray ( im, jm, km, wa );
	Array	c;																					// water vapour
	c.initArray ( im, jm, km, ca );
	Array	cloud;																			// cloud water
	cloud.initArray ( im, jm, km, 0. );
	Array	ice;																				// cloud ice
	ice.initArray ( im, jm, km, 0. );
	Array	co2;																				// CO2
	co2.initArray ( im, jm, km, coa );

	Array	tn;																					// temperature new
	tn.initArray ( im, jm, km, ta );
	Array	un;																					// u-velocity component in r-direction new
	un.initArray ( im, jm, km, ua );
	Array	vn;																					// v-velocity component in theta-direction new
	vn.initArray ( im, jm, km, va );
	Array	wn;																				// w-velocity component in phi-direction new
	wn.initArray ( im, jm, km, wa );
	Array	cn;																					// water vapour new
	cn.initArray ( im, jm, km, ca );
	Array	cloudn;																			// cloud water new
	cloudn.initArray ( im, jm, km, 0. );
	Array	icen;																				// cloud ice new
	icen.initArray ( im, jm, km, 0. );
	Array	co2n;																				// CO2 new
	co2n.initArray ( im, jm, km, coa );

	Array	p_dyn;																			// dynamic pressure
	p_dyn.initArray ( im, jm, km, pa );
	Array	p_stat;																			// static pressure
	p_stat.initArray ( im, jm, km, pa );

	Array	rhs_t;																			// auxilliar field RHS temperature
	rhs_t.initArray ( im, jm, km, 0. );
	Array	rhs_u;																			// auxilliar field RHS u-velocity component
	rhs_u.initArray ( im, jm, km, 0. );
	Array	rhs_v;																			// auxilliar field RHS v-velocity component
	rhs_v.initArray ( im, jm, km, 0. );
	Array	rhs_w;																			// auxilliar field RHS w-velocity component
	rhs_w.initArray ( im, jm, km, 0. );
	Array	rhs_c;																			// auxilliar field RHS water vapour
	rhs_c.initArray ( im, jm, km, 0. );
	Array	rhs_cloud;																	// auxilliar field RHS cloud water
	rhs_cloud.initArray ( im, jm, km, 0. );
	Array	rhs_ice;																		// auxilliar field RHS cloud ice
	rhs_ice.initArray ( im, jm, km, 0. );
	Array	rhs_co2;																		// auxilliar field RHS CO2
	rhs_co2.initArray ( im, jm, km, 0. );

	Array	aux_u;																			// auxilliar field u-velocity component
	aux_u.initArray ( im, jm, km, 0. );
	Array	aux_v;																			// auxilliar field v-velocity component
	aux_v.initArray ( im, jm, km, 0. );
	Array	aux_w;																		// auxilliar field w-velocity component
	aux_w.initArray ( im, jm, km, 0. );
	Array	aux_p;																			// auxilliar field p
	aux_p.initArray ( im, jm, km, pa );

	Array	Latency;																		// latent heat
	Latency.initArray ( im, jm, km, 0. );
	Array	Q_Sensible;																	// sensible heat
	Q_Sensible.initArray ( im, jm, km, 0. );
	Array	IceLayer;																		// ice shield
	IceLayer.initArray ( im, jm, km, 0. );
	Array	t_cond_3D;																	// condensation temperature
	t_cond_3D.initArray ( im, jm, km, 0. );
	Array	t_evap_3D;																	// evaporation temperature
	t_evap_3D.initArray ( im, jm, km, 0. );
	Array	BuoyancyForce;															// buoyancy force, Boussinesque approximation
	BuoyancyForce.initArray ( im, jm, km, 0. );
	Array	epsilon_3D;																	// emissivity/ absorptivity
	epsilon_3D.initArray ( im, jm, km, 0. );
	Array	radiation_3D;																// radiation
	radiation_3D.initArray ( im, jm, km, 0. );

	Array	P_rain;																			// rain precipitation mass rate
	P_rain.initArray ( im, jm, km, 0. );
	Array	P_snow;																		// snow precipitation mass rate
	P_snow.initArray ( im, jm, km, 0. );
	Array	S_v;																				// water vapour mass rate due to category two ice scheme
	S_v.initArray ( im, jm, km, 0. );
	Array	S_c;																				// cloud water mass rate due to category two ice scheme
	S_c.initArray ( im, jm, km, 0. );
	Array	S_i;																				// cloud ice mass rate due to category two ice scheme
	S_i.initArray ( im, jm, km, 0. );
	Array	S_r;																				// rain mass rate due to category two ice scheme
	S_r.initArray ( im, jm, km, 0. );
	Array	S_s;																				// snow mass rate due to category two ice scheme
	S_s.initArray ( im, jm, km, 0. );



//	cout << endl << " ***** printout of 3D-field temperature ***** " << endl << endl;
//	t.printArray();

//	cout << endl << " ***** printout of 2D-field vegetation ***** " << endl << endl;
//	Vegetation.printArray_2D();

//	cout << endl << " ***** printout of 1D-field radius ***** " << endl << endl;
//	rad.printArray_1D();

	cout.precision ( 6 );
	cout.setf ( ios::fixed );


//	Coordinate system in form of a spherical shell
//	rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
	rad.Coordinates ( im, r0, dr );
	the.Coordinates ( jm, the0, dthe );
	phi.Coordinates ( km, phi0, dphi );

//	VALGRIND_CHECK_VALUE_IS_DEFINED ( rad );

//	coeff_mmWS = r_air / r_water_vapour;	// coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]

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

// time slice to start with the modern world
	Ma = 0;
	i_time_slice = 0;
//	Ma = 50;
//	i_time_slice = 5;


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



// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//	choice of the time slice to be computed
	time_slice_sequel:

// choice of the time slice by Ma and by author ( etopo for modern times, Galonka for time slices )
	if ( Ma == 0 )
	{
		Name_Bathymetry_File = "0Ma_etopo.xyz";
		Name_Sequel_File = "[0Ma_etopo.xyz]_Sequel_Atm.seq";
		Name_netCDF_File = "[0Ma_etopo.xyz]_atmosphere.nc";
	} 
	else 
	{
		n = 0;
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
		read_File.Atmosphere_SequelFile_read ( Name_Bathymetry_File, n, time, rad, the, phi, h, t, u, v, w, c, co2, tn, un, vn, wn, cn, co2n );
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



//  initialization of the bathymetry/topography

// 	class BC_Bathymetry_Atmosphere for the geometrical boundary condition of the computational area
	BC_Bathymetry_Atmosphere		LandArea ( im, jm, km, co2_vegetation, co2_land, co2_ocean );

// 	topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
	LandArea.BC_MountainSurface ( Name_Bathymetry_File, L_atm, h, aux_w );

//	computation of ice shield following the theorie by Milankowitsch
	if ( IceShield == 1 ) LandArea.BC_IceShield ( Ma, t_0, h, t, c, IceLayer, Ice_Balance, Ice_Balance_add );




// 	class calls for the solution of the flow properties

// 	class BC_Atmosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
	BC_Atmosphere		boundary ( im, jm, km, t_tropopause );

//  class RHS_Atmosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
	RHS_Atmosphere		prepare ( im, jm, km, dt, dr, dthe, dphi, re, ec, sc_WaterVapour, sc_CO2, g, pr, omega, coriolis, centrifugal, WaterVapour, buoyancy, CO2, gam, sigma, lambda );
	RHS_Atmosphere		prepare_2D ( jm, km, dthe, dphi, re, omega, coriolis, centrifugal );

//  class RungeKutta_Atmosphere for the explicit solution of the Navier-Stokes equations
	RungeKutta_Atmosphere		result ( im, jm, km, dt, dr, dphi, dthe );

//  class Pressure for the subsequent computation of the pressure by a separat Euler equation
	Pressure		startPressure ( im, jm, km, dr, dthe, dphi );

//	class Restore to restore the iterational values from new to old
	Restore		oldnew( im, jm, km );

//	class Results_MSL_Atm to compute and show results on the mean sea level, MSL
	Results_MSL_Atm		calculate_MSL ( im, jm, km, sun, g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_extra, lv, ls, cp_l, L_atm, dt, dr, dthe, dphi, r_air, R_Air, r_water, r_water_vapour, R_WaterVapour, co2_vegetation, co2_ocean, co2_land, gam, t_pole, t_cretaceous );

//	class File_NetCDF to write results in the format of a netCDF-file
	File_NetCDF		printoutNetCDF ( im, jm, km );



//  configuration of the initial and boundary conditions for the temperature, CO2 und water vapour on land and ocean surfaces

// 	class BC_Thermo for the initial and boundary conditions of the flow properties
	BC_Thermo		circulation ( im, jm, km, i_beg, i_max, RadiationModel, sun, declination, sun_position_lat, sun_position_lon, Ma, Ma_max, Ma_max_half, dr, dthe, dphi, g, ep, hp, u_0, p_0, t_0, c_0, sigma, albedo_extra, epsilon_extra, lv, cp_l, L_atm, r_air, R_Air, r_water_vapour, R_WaterVapour, co2_0, co2_cretaceous, co2_vegetation, co2_ocean, co2_land, ik, c_tropopause, co2_tropopause, c_ocean, c_land, t_average, co2_average, co2_pole, t_cretaceous, t_cretaceous_max, radiation_ocean, radiation_pole, radiation_equator, t_land, t_tropopause, t_equator, t_pole, gam, false, verbose );


//  class element for the tropopause location as a parabolic distribution from pole to pole 
	circulation.TropopauseLocation ( im_tropopause );

//  class element for the surface temperature from World Ocean Atlas 2009 given as boundary condition
//	if ( Ma == 0 ) circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );
//	circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );

//  class element for the surface precipitation from NASA for comparison
//	if ( Ma == 0 ) circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );
//	circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );

//  class element for the parabolic temperature distribution from pol to pol, maximum temperature at equator
	circulation.BC_Temperature ( h, t, p_dyn, p_stat );
	t_cretaceous = circulation.out_temperature (  );

//  class element for the surface pressure computed by surface temperature with gas equation
	circulation.BC_Pressure ( p_stat, t, h );

//  parabolic water vapour distribution from pol to pol, maximum water vapour volume at equator
	circulation.BC_WaterVapour ( h, t, p_stat, c, cloud, ice, P_rain, P_snow, precipitation_NASA );

//	circulation.IC_Cloud_Ice ( h, c, t, p_stat, cloud, ice );


//  class element for the parabolic CO2 distribution from pol to pol, maximum CO2 volume at equator
	circulation.BC_CO2 ( Vegetation, h, t, p_dyn, co2 );
	co2_cretaceous = circulation.out_co2 (  );

// 	class element for the initial conditions for u-v-w-velocity components
	circulation.IC_CellStructure ( im_tropopause, u, v, w );

// class element for the surface temperature computation by radiation flux density
	if ( RadiationModel == 3 ) circulation.BC_Radiation_multi_layer ( n, t_j, albedo, epsilon, precipitable_water, Ik, Q_Radiation, Radiation_Balance, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, temp_rad, Q_latent, Q_sensible, Q_bottom, co2_total, p_stat, t, c, h, epsilon_3D, radiation_3D, cloud, ice );

// class element for the surface temperature computation by radiation flux density
//	if ( RadiationModel >= 2 ) circulation.BC_Radiation_2D_layer ( im_tropopause, albedo, epsilon, precipitable_water, Ik, Q_Radiation, Radiation_Balance, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, temp_rad, Q_latent, Q_sensible, Q_bottom, co2_total, t, c, h, epsilon_3D, radiation_3D );

//  class element for the parabolic radiation balance distribution from pole to pole, maximum radiation balance amount at equator
	if ( RadiationModel == 1 ) circulation.BC_Radiation_parabolic ( Radiation_Balance_par, h );

//	class element for the initial conditions the latent heat
//	circulation.Latent_Heat ( rad, the, phi, h, t, tn, u, v, w, p_dyn, p_stat, c, Latency, Q_Sensible, t_cond_3D, t_evap_3D, radiation_3D );

//	class element for the storing of velocity components, pressure and temperature for iteration start
	oldnew.restoreOldNew_3D ( .9, u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, aux_p, cn, cloudn, icen, co2n );
	oldnew.restoreOldNew_2D ( .9, v, w, p_dyn, aux_p, vn, wn );

// class element for the computation of the ratio ocean to land areas, also supply and removal of CO2 on land, ocean and by vegetation
	LandArea.land_oceanFraction ( h );



/*
	cout << endl << " ***** printout of 3D-field water vapour ***** " << endl << endl;
	c.printArray();

	cout << endl << " ***** printout of 3D-field water cloud ***** " << endl << endl;
	cloud.printArray();

	cout << endl << " ***** printout of 3D-field water ice ***** " << endl << endl;
	ice.printArray();

	cout << endl << " ***** printout of 3D-field temperature ***** " << endl << endl;
	t.printArray();

	cout << endl << " ***** printout of 3D-field epsilon ***** " << endl << endl;
	epsilon_3D.printArray();

	cout << endl << " ***** printout of 3D-field radiation_3D ***** " << endl << endl;
	radiation_3D.printArray();

	cout << endl << " ***** printout of 3D-field p_stat ***** " << endl << endl;
	p_stat.printArray();

	goto Print_commands;
*/


// ******************************************   start of pressure and velocity iterations ************************************************************************


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of pressure loop : if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::

Pressure_loop:

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of velocity loop: while ( min >= epsres )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//	min = min_u = min_v = min_w = min_t = min_c = min_p = epsres * 3.;
	min = epsres * 3.;
	velocity_iter = 0;
	velocity_iter_2D = 0;


//	query to realize zero divergence of the continuity equation ( div c = 0 )
	while ( min >= epsres )
	{
//		limit of the computation in the sense of time stepscoeff_P
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

		cout << endl << "  present state of the computation " << endl << "  current time slice, number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations " << endl << endl << "  Ma = " << Ma << "     n = " << n << "    velocity_iter_max = " << velocity_iter_max << "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << "    pressure_iter = " << pressure_iter << endl;

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
			boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
			boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

// 		old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy		min_Residuum_old_2D ( im, jm, km, dthe, dphi );
			min_Residuum_old_2D.residuumQuery_2D ( rad, the, v, w );
			min = min_Residuum_old_2D.out_min (  );

			residuum_old = min;

//		class RungeKutta for the solution of the differential equations describing the flow properties
			result.solveRungeKutta_2D_Atmosphere ( prepare_2D, rad, the, rhs_v, rhs_w, h, v, w, p_dyn, vn, wn, aux_v, aux_w );

// 		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy		min_Residuum_2D ( im, jm, km, dthe, dphi );
			min_Residuum_2D.residuumQuery_2D ( rad, the, v, w );
			min = min_Residuum_2D.out_min (  );
			j_res = min_Residuum_2D.out_j_res (  );
			k_res = min_Residuum_2D.out_k_res (  );

			residuum = min;
			min = fabs ( ( residuum - residuum_old ) / residuum_old );

//		state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
			Accuracy		min_Stationary_2D ( n, nm, Ma, im, jm, km, min, j_res, k_res, velocity_iter_2D, pressure_iter_2D, velocity_iter_max_2D, pressure_iter_max_2D );
			min_Stationary_2D.steadyQuery_2D ( v, vn, w, wn, p_dyn, aux_p );

			oldnew.restoreOldNew_2D ( 1., v, w, p_dyn, aux_p, vn, wn );


		}



//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop_2D: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Pressure_iteration_2D:

//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
	startPressure.computePressure_2D ( pa, rad, the, p_dyn, h, rhs_v, rhs_w, aux_v, aux_w, aux_p );		// 2D pressure computation causes a pressure jump in radial direction along coast lines, 3D treatment needed later, 2D velocities are though corrected

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
		Accuracy		min_Residuum_old ( im, jm, km, dr, dthe, dphi );
		min_Residuum_old.residuumQuery_3D ( rad, the, u, v, w );
		min = min_Residuum_old.out_min (  );

		residuum_old = min;

// 		class BC_Atmosphaere for the geometry of a shell of a sphere
		boundary.BC_radius ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
		boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
		boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

// 		class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
		LandArea.BC_SolidGround ( RadiationModel, i_max, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, co2_cretaceous, pa, gam, h, u, v, w, t, p_dyn, c, cloud, ice, co2, Vegetation );

// 		class RungeKutta for the solution of the differential equations describing the flow properties
		result.solveRungeKutta_3D_Atmosphere ( prepare, n, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_cloud, rhs_ice, rhs_co2, h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, cn, cloudn, icen, co2n, aux_u, aux_v, aux_w, Latency, t_cond_3D, t_evap_3D, IceLayer, BuoyancyForce, Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s );

// 		class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
		LandArea.BC_SolidGround ( RadiationModel, i_max, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, co2_cretaceous, pa, gam, h, u, v, w, t, p_dyn, c, cloud,ice, co2, Vegetation );

// class element for the surface temperature computation by radiation flux density
		if ( RadiationModel == 3 ) 					circulation.BC_Radiation_multi_layer ( n, t_j, albedo, epsilon, precipitable_water, Ik, Q_Radiation, Radiation_Balance, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, temp_rad, Q_latent, Q_sensible, Q_bottom, co2_total, p_stat, t, c, h, epsilon_3D, radiation_3D, cloud, ice );

// class element for the surface temperature computation by radiation flux density
//		if ( RadiationModel >= 2 ) 					circulation.BC_Radiation_2D_layer ( im_tropopause, albedo, epsilon, precipitable_water, Ik, Q_Radiation, Radiation_Balance, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, temp_rad, Q_latent, Q_sensible, Q_bottom, co2_total, t, c, h, epsilon_3D, radiation_3D );

// 		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
		Accuracy		min_Residuum ( im, jm, km, dr, dthe, dphi );
		min_Residuum.residuumQuery_3D ( rad, the, u, v, w );
		min = min_Residuum.out_min (  );
		i_res = min_Residuum.out_i_res (  );
		j_res = min_Residuum.out_j_res (  );
		k_res = min_Residuum.out_k_res (  );

		residuum = min;
		min = fabs ( ( residuum - residuum_old ) / residuum_old );

//		statements on the convergence und iterational process
		Accuracy		min_Stationary ( n, nm, Ma, im, jm, km, min, i_res, j_res, k_res, velocity_iter, pressure_iter, velocity_iter_max, pressure_iter_max, L_atm );
		min_Stationary.steadyQuery_3D ( u, un, v, vn, w, wn, t, tn, c, cn, cloud, cloudn, ice, icen, co2, co2n, p_dyn, aux_p );
/*
	cout << endl << " ***** printout of 3D-field water vapour ***** " << endl << endl;
	c.printArray();

	cout << endl << " ***** printout of 3D-field water cloud ***** " << endl << endl;
	cloud.printArray();

	cout << endl << " ***** printout of 3D-field water ice ***** " << endl << endl;
	ice.printArray();

	cout << endl << " ***** printout of 3D-field temperature ***** " << endl << endl;
	t.printArray();

	cout << endl << " ***** printout of 3D-field epsilon ***** " << endl << endl;
	epsilon_3D.printArray();

	cout << endl << " ***** printout of 3D-field radiation_3D ***** " << endl << endl;
	radiation_3D.printArray();

	cout << endl << " ***** printout of 3D-field p_stat ***** " << endl << endl;
	p_stat.printArray();

	cout << endl << " ***** printout of 3D-field p_dyn ***** " << endl << endl;
	p_dyn.printArray();

	cout << endl << " ***** printout of 2D-field precipitable_water ***** " << endl << endl;
	precipitable_water.printArray_2D();

	cout << endl << " ***** printout of 2D-field Precipitation ***** " << endl << endl;
	Precipitation.printArray_2D();

	cout << endl << " ***** printout of 3D-field P_rain ***** " << endl << endl;
	P_rain.printArray();

	cout << endl << " ***** printout of 3D-field P_snow ***** " << endl << endl;
	P_snow.printArray();

	cout << endl << " ***** printout of 3D-field co2 ***** " << endl << endl;
	co2.printArray();
*/


//	goto Print_commands;

//	if ( n >= 3 ) goto Print_commands;

// 3D_fields

//		searching of maximum and minimum values of temperature
		string str_max_temperature = " max 3D temperature ", str_min_temperature = " min 3D temperature ", str_unit_temperature = "C";
		MinMax		minmaxTemperature ( im, jm, km );
		minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

//		searching of maximum and minimum values of dynamic pressure
		string str_max_pressure = " max 3D pressure dynamic ", str_min_pressure = " min 3D pressure dynamic ", str_unit_pressure = "hPa";
		MinMax		minmaxPressure ( im, jm, km );
		minmaxPressure.searchMinMax_3D ( str_max_pressure, str_min_pressure, str_unit_pressure, p_dyn, h );

//		searching of maximum and minimum values of static pressure
		string str_max_pressure_stat = " max 3D pressure static ", str_min_pressure_stat = " min 3D pressure static ", str_unit_pressure_stat = "hPa";
		MinMax		minmaxPressure_stat ( im, jm, km );
		minmaxPressure_stat.searchMinMax_3D ( str_max_pressure_stat, str_min_pressure_stat, str_unit_pressure_stat, p_stat, h );

		cout << endl << " energies in the three dimensional space: " << endl << endl;

//		searching of maximum and minimum values of radiation_3D
		string str_max_radiation_3D = " max 3D radiation ", str_min_radiation_3D = " min 3D radiation ", str_unit_radiation_3D = "W/m2";
		MinMax		minmaxLatency ( im, jm, km );
		minmaxLatency.searchMinMax_3D ( str_max_radiation_3D, str_min_radiation_3D, str_unit_radiation_3D, radiation_3D, h );

//		searching of maximum and minimum values of sensible heat
		string str_max_Q_Sensible = " max 3D sensible heat ", str_min_Q_Sensible = " min 3D sensible heat ", str_unit_Q_Sensible = "W/m2";
		MinMax		minmaxQ_Sensible ( im, jm, km );
		minmaxQ_Sensible.searchMinMax_3D ( str_max_Q_Sensible, str_min_Q_Sensible, str_unit_Q_Sensible, Q_Sensible, h );

//		searching of maximum and minimum values of latency
		string str_max_latency = " max 3D latent heat ", str_min_latency = " min 3D latent heat ", str_unit_latency = "W/m2";
		MinMax		minmaxRadiation ( im, jm, km );
		minmaxRadiation.searchMinMax_3D ( str_max_latency, str_min_latency, str_unit_latency, Latency, h );

//		searching of maximum and minimum values of t_cond_3D
		string str_max_t_cond_3D = " max 3D condensation temp ", str_min_t_cond_3D = " min 3D condensation temp ", str_unit_t_cond_3D = "C";
		MinMax		minmaxt_cond_3D ( im, jm, km );
		minmaxt_cond_3D.searchMinMax_3D ( str_max_t_cond_3D, str_min_t_cond_3D, str_unit_t_cond_3D, t_cond_3D, h );

//		searching of maximum and minimum values of t_evap_3D
		string str_max_t_evap_3D = " max 3D evaporation temp ", str_min_t_evap_3D = " min 3D evaporation temp ", str_unit_t_evap_3D = "C";
		MinMax		minmaxt_evap_3D ( im, jm, km );
		minmaxt_evap_3D.searchMinMax_3D ( str_max_t_evap_3D, str_min_t_evap_3D, str_unit_t_evap_3D, t_evap_3D, h );

		cout << endl << " greenhouse gases: " << endl << endl;

//		searching of maximum and minimum values of water vapour
		string str_max_water_vapour = " max 3D water vapour ", str_min_water_vapour = " min 3D water vapour ", str_unit_water_vapour = "g/kg";
		MinMax		minmaxWaterVapour ( im, jm, km );
		minmaxWaterVapour.searchMinMax_3D ( str_max_water_vapour, str_min_water_vapour, str_unit_water_vapour, c, h );

//		searching of maximum and minimum values of cloud water
		string str_max_cloud_water = " max 3D cloud water ", str_min_cloud_water = " min 3D cloud water ", str_unit_cloud_water = "g/kg";
		MinMax		minmaxCloudWater ( im, jm, km );
		minmaxCloudWater.searchMinMax_3D ( str_max_cloud_water, str_min_cloud_water, str_unit_cloud_water, cloud, h );

//		searching of maximum and minimum values of cloud ice
		string str_max_cloud_ice = " max 3D cloud ice ", str_min_cloud_ice = " min 3D cloud ice ", str_unit_cloud_ice = "g/kg";
		MinMax		minmaxCloudIce ( im, jm, km );
		minmaxCloudIce.searchMinMax_3D ( str_max_cloud_ice, str_min_cloud_ice, str_unit_cloud_ice, ice, h );

//		searching of maximum and minimum values of rain precipitation
		string str_max_P_rain = " max 3D rain ", str_min_P_rain = " min 3D rain ", str_unit_P_rain = "g/kg";
		MinMax		minmaxPRain ( im, jm, km );
		minmaxPRain.searchMinMax_3D ( str_max_P_rain, str_min_P_rain, str_unit_P_rain, P_rain, h );

//		searching of maximum and minimum values of snow precipitation
		string str_max_P_snow = " max 3D snow ", str_min_P_snow = " min 3D snow ", str_unit_P_snow = "g/kg";
		MinMax		minmaxPSnow ( im, jm, km );
		minmaxPSnow.searchMinMax_3D ( str_max_P_snow, str_min_P_snow, str_unit_P_snow, P_snow, h );

//		searching of maximum and minimum values of co2
		string str_max_co2 = " max 3D co2 ", str_min_co2 = " min 3D co2 ", str_unit_co2 = "ppm";
		MinMax		minmaxCO2 ( im, jm, km );
		minmaxCO2.searchMinMax_3D ( str_max_co2, str_min_co2, str_unit_co2, co2, h );

//		searching of maximum and minimum values of epsilon
		string str_max_epsilon = " max 3D epsilon ", str_min_epsilon = " min 3D epsilon ", str_unit_epsilon = "%";
		MinMax		minmaxEpsilon_3D ( im, jm, km );
		minmaxEpsilon_3D.searchMinMax_3D ( str_max_epsilon, str_min_epsilon, str_unit_epsilon, epsilon_3D, h );

//		searching of maximum and minimum values of buoyancy force
		string str_max_buoyancy_force = " max 3D buoyancy force ", str_min_buoyancy_force = " min 3D buoyancy force ", str_unit_buoyancy_force = "N/m2";
		MinMax		minmaxBuoyancyForce ( im, jm, km );
		minmaxBuoyancyForce.searchMinMax_3D ( str_max_buoyancy_force, str_min_buoyancy_force, str_unit_buoyancy_force, BuoyancyForce, h );



// 2D-fields

//		searching of maximum and minimum values of co2 total
		cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" << endl << " results based on two dimensional considerations of the problem" << endl;

		cout << endl << " co2 distribution columnwise: " << endl << endl;

		string str_max_co2_total = " max co2_total ", str_min_co2_total = " min co2_total ", str_unit_co2_total = " / ";
		MinMax		minmaxCO2_total ( jm, km, coeff_mmWS );
		minmaxCO2_total.searchMinMax_2D ( str_max_co2_total, str_min_co2_total, str_unit_co2_total, co2_total, h );
//		max_CO2_total = minmaxCO2_total.out_maxValue (  );

		cout << endl << " precipitation: " << endl << endl;

//		searching of maximum and minimum values of precipitation
		string str_max_precipitation = " max precipitation ", str_min_precipitation = " min precipitation ", str_unit_precipitation = "mm";
		MinMax		minmaxPrecipitation ( jm, km, coeff_mmWS );
		minmaxPrecipitation.searchMinMax_2D ( str_max_precipitation, str_min_precipitation, str_unit_precipitation, Precipitation, h );
		max_Precipitation = minmaxPrecipitation.out_maxValue (  );

/*
//		searching of maximum and minimum values of NASA precipitation
		if ( Ma == 0 )
		{
			string str_max_precipitation = " max precipitation_NASA ", str_min_precipitation = " min precipitation_NASA ", str_unit_precipitation = "mm";
			MinMax		minmaxPrecipitation_NASA ( jm, km, coeff_mmWS );
			minmaxPrecipitation_NASA.searchMinMax_2D ( str_max_precipitation, str_min_precipitation, str_unit_precipitation, precipitation_NASA, h );
		}
*/
//		searching of maximum and minimum values of precipitable water
		string str_max_precipitable_water = " max precipitable water ", str_min_precipitable_water = " min precipitable water ", str_unit_precipitable_water = "mm";
		MinMax		minmaxPrecipitable_water ( jm, km, coeff_mmWS );
		minmaxPrecipitable_water.searchMinMax_2D ( str_max_precipitable_water, str_min_precipitable_water, str_unit_precipitable_water, precipitable_water, h );
//		max_precipitable_water = minmaxPrecipitable_water.out_minValue (  );

		cout << endl << " energies at see level without convection influence: " << endl << endl;

//		searching of maximum and minimum values of radiation balance
//		string str_max_Radiation_Balance = " max radiation balance ", str_min_Radiation_Balance = " min radiation balance ", str_unit_Radiation_Balance = "W/m2";
//		MinMax		minmaxRadiation_Balance ( jm, km, coeff_mmWS );
//		minmaxRadiation_Balance.searchMinMax_2D ( str_max_Radiation_Balance, str_min_Radiation_Balance, str_unit_Radiation_Balance, Radiation_Balance, h );
//		min_Radiation_Balance = minmaxRadiation_Balance.out_minValue (  );

//		searching of maximum and minimum values of radiation
		string str_max_Q_Radiation = " max 2D Q radiation ", str_min_Q_Radiation = " min 2D Q radiation ", str_unit_Q_Radiation = "W/m2";
		MinMax		minmaxQ_Radiation ( jm, km, coeff_mmWS );
		minmaxQ_Radiation.searchMinMax_2D ( str_max_Q_Radiation, str_min_Q_Radiation, str_unit_Q_Radiation, Q_Radiation, h );
//		min_Q_Radiation = minmaxQ_Radiation.out_minValue (  );

//		searching of maximum and minimum values of latent energy
		string str_max_Q_latent = " max 2D Q latent ", str_min_Q_latent = " min 2D Q latent ", str_unit_Q_latent = "W/m2";
		MinMax		minmaxQ_latent ( jm, km, coeff_mmWS );
		minmaxQ_latent.searchMinMax_2D ( str_max_Q_latent, str_min_Q_latent, str_unit_Q_latent, Q_latent, h );

//		searching of maximum and minimum values of sensible energy
		string str_max_Q_sensible = " max 2D Q sensible ", str_min_Q_sensible = " min 2D Q sensible ", str_unit_Q_sensible = "W/m2";
		MinMax		minmaxQ_sensible ( jm, km, coeff_mmWS );
		minmaxQ_sensible.searchMinMax_2D ( str_max_Q_sensible, str_min_Q_sensible, str_unit_Q_sensible, Q_sensible, h );

//		searching of maximum and minimum values of bottom heat

		string str_max_Q_bottom = " max 2D Q bottom ", str_min_Q_bottom = " min 2D Q bottom heat ", str_unit_Q_bottom = "W/m2";
		MinMax		minmaxQ_bottom ( jm, km, coeff_mmWS );
		minmaxQ_bottom.searchMinMax_2D ( str_max_Q_bottom, str_min_Q_bottom, str_unit_Q_bottom, Q_bottom, h );


//		searching of maximum and minimum values of latent heat
		string str_max_LatentHeat = " max 2D latent heat ", str_min_LatentHeat = " min 2D latent heat ", str_unit_LatentHeat = "W/m2";
		MinMax		minmaxLatentHeat_2D ( jm, km, coeff_mmWS );
		minmaxLatentHeat_2D.searchMinMax_2D ( str_max_LatentHeat, str_min_LatentHeat, str_unit_LatentHeat, LatentHeat, h );

//		searching of maximum and minimum values of t_cond_2D
		string str_max_t_cond = " max 2D Condensation ", str_min_t_cond = " min 2D Condensation ", str_unit_t_cond = "W/m2";
		MinMax		minmaxt_cond_2D ( jm, km, coeff_mmWS );
		minmaxt_cond_2D.searchMinMax_2D ( str_max_t_cond, str_min_t_cond, str_unit_t_cond, Condensation, h );


//		searching of maximum and minimum values of t_evap_2D
		string str_max_t_evap = " max 2D Evaporation ", str_min_t_evap = " min 2D Evaporation ", str_unit_t_evap = "W/m2";
		MinMax		minmaxt_evap_2D ( jm, km, coeff_mmWS );
		minmaxt_evap_2D.searchMinMax_2D ( str_max_t_evap, str_min_t_evap, str_unit_t_evap, Evaporation, h );

		cout << endl << " secondary data: " << endl << endl;


/*
//		searching of maximum and minimum values of Evaporation
		string str_max_heat_t_Evaporation = " max heat Evaporation ", str_min_heat_t_Evaporation = " min heat Evaporation ", str_unit_heat_t_Evaporation = " W/m2";
		MinMax		minmaxQ_t_Evaporation ( jm, km, coeff_mmWS );
		minmaxQ_t_Evaporation.searchMinMax_2D ( str_max_heat_t_Evaporation, str_min_heat_t_Evaporation, str_unit_heat_t_Evaporation, Q_Evaporation, h );
*/
/*
//		searching of maximum and minimum values of Evaporation by Haude
		string str_max_t_Evaporation_Haude = " max Evaporation Haude ", str_min_t_Evaporation_Haude = " min Evaporation Haude ", str_unit_t_Evaporation_Haude = "mm/d";
		MinMax		minmaxt_Evaporation_Haude ( jm, km, coeff_mmWS );
		minmaxt_Evaporation_Haude.searchMinMax_2D ( str_max_t_Evaporation_Haude, str_min_t_Evaporation_Haude, str_unit_t_Evaporation_Haude, Evaporation_Haude, h );
*/
//		searching of maximum and minimum values of Evaporation by Penman
		string str_max_t_Evaporation_Penman = " max Evaporation Penman ", str_min_t_Evaporation_Penman = " min Evaporation Penman ", str_unit_t_Evaporation_Penman = "mm/d";
		MinMax		minmaxt_Evaporation_Penman ( jm, km, coeff_mmWS );
		minmaxt_Evaporation_Penman.searchMinMax_2D ( str_max_t_Evaporation_Penman, str_min_t_Evaporation_Penman, str_unit_t_Evaporation_Penman, Evaporation_Penman, h );

		cout << endl << " properties of the atmosphere at the surface: " << endl << endl;

//		searching of maximum and minimum values of albedo
		string str_max_albedo = " max 2D albedo ", str_min_albedo = " min 2D albedo ", str_unit_albedo = "%";
		MinMax		minmaxAlbedo ( jm, km, coeff_mmWS );
		minmaxAlbedo.searchMinMax_2D ( str_max_albedo, str_min_albedo, str_unit_albedo, albedo, h );
/*
//		searching of maximum and minimum values of epsilon
		string str_max_epsilon = " max 2D epsilon ", str_min_epsilon = " min 2D epsilon ", str_unit_epsilon = "%";
		MinMax		minmaxEpsilon ( jm, km, coeff_mmWS );
		minmaxEpsilon.searchMinMax_2D ( str_max_epsilon, str_min_epsilon, str_unit_epsilon, epsilon, h );
*/

//		computation of vegetation areas
		LandArea.vegetationDistribution ( max_Precipitation, Precipitation, Vegetation, t, h );


//	class element for the initial conditions the latent heat
		circulation.Latent_Heat ( rad, the, phi, h, t, tn, u, v, w, p_dyn, p_stat, c, Latency, Q_Sensible, t_cond_3D, t_evap_3D, radiation_3D );


//		composition of results
		calculate_MSL.run_MSL_data ( n, velocity_iter_max, RadiationModel, rad, the, phi, h, c, cn, co2, co2n, t, tn, p_dyn, p_stat, BuoyancyForce, u, v, w, Latency, Q_Sensible, radiation_3D, t_cond_3D, t_evap_3D, cloud, cloudn, ice, icen, P_rain, P_snow, aux_u, aux_v, aux_w, precipitation_NASA, Evaporation, Condensation, LatentHeat, precipitable_water, Q_Radiation, Q_Evaporation, Q_latent, Q_sensible, Q_bottom, Evaporation_Penman, Evaporation_Haude, Vegetation, Radiation_Balance, Radiation_Balance_par, Radiation_Balance_bot, albedo, co2_total, Precipitation, S_v, S_c, S_i, S_r, S_s );



// printout of results at certain positions
//		calculate_MSL.show_MSL_data ( h, c, t, p_dyn, u, Latency, Q_Sensible, t_cond_3D, t_evap_3D, precipitation_NASA, Evaporation, Condensation, precipitable_water, Q_Radiation, Q_Evaporation, Q_latent, Q_sensible, Q_bottom, Evaporation_Penman, Evaporation_Haude );

//	restoring the velocity component and the temperature for the new time step
		oldnew.restoreOldNew_3D ( 1., u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, aux_p, cn, cloudn, icen, co2n );

	}


//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
	startPressure.computePressure_3D ( pa, rad, the, p_dyn, h, rhs_u, rhs_v, rhs_w, aux_u, aux_v, aux_w, aux_p );

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
	printoutNetCDF.out_NetCDF( Name_netCDF_File, v, w, h, Precipitation, precipitable_water );


//	class PostProcess_Atmosphaere for the printing of results
	PostProcess_Atmosphere		write_File ( im, jm, km );


//	writing of data in ParaView files

//	radial data along constant hight above ground
	i_radial = 0;
	write_File.paraview_vtk_radial ( Name_Bathymetry_File, i_radial, pressure_iter_aux, u_0, t_0, p_0, r_air, c_0, co2_0, radiation_equator, h, p_dyn, p_stat, t_cond_3D, t_evap_3D , BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain, P_snow, Evaporation, Condensation, precipitable_water, Q_bottom, Radiation_Balance, Q_Radiation, Q_latent, Q_sensible, Evaporation_Penman, Evaporation_Haude, Q_Evaporation, precipitation_NASA, Vegetation, albedo, epsilon, Precipitation );

//	londitudinal data along constant latitudes
	j_longal = 75;
	write_File.paraview_vtk_longal ( Name_Bathymetry_File, j_longal, pressure_iter_aux, u_0, t_0, p_0, r_air, c_0, co2_0, radiation_equator, h, p_dyn, p_stat, t_cond_3D, t_evap_3D, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain, P_snow );

	k_zonal = 145;
	write_File.paraview_vtk_zonal ( Name_Bathymetry_File, k_zonal, pressure_iter_aux, u_0, t_0, p_0, r_air, c_0, co2_0, radiation_equator, h, p_dyn, p_stat, t_cond_3D, t_evap_3D, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, radiation_3D, epsilon_3D, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s );

//	3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
//	write_File.paraview_panorama_vts ( Name_Bathymetry_File, pressure_iter_aux, u_0, t_0, p_0, r_air, c_0, co2_0, h, t, p_dyn, p_stat, BuoyancyForce, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain );

//	3-dimensional data in spherical coordinate system for a streamline pattern in a shell of a sphere
//	write_File.paraview_vts ( Name_Bathymetry_File, n, rad, the, phi, h, t, p_dyn, u, v, w, c, co2, aux_u, aux_v, aux_w, Latency, Rain, Ice, Rain_super, IceLayer );


//	writing of sequential data for the sequel file
	if ( SequelFile == 1 )
	{
		write_File.Atmosphere_SequelFile_write ( Name_Bathymetry_File, n, time, rad, the, phi, h, t, u, v, w, c, co2, tn, un, vn, wn, cn, co2n );
	}

//	writing of v-w-data in the v_w_transfer file
		PostProcess_Atmosphere		write_v_w_Transfer_File ( im, jm, km );
		write_v_w_Transfer_File.Atmosphere_v_w_Transfer ( Name_Bathymetry_File, v, w, p_dyn );


//	writing of plot data in the PlotData file
		PostProcess_Atmosphere		write_PlotData_File ( im, jm, km );
		write_PlotData_File.Atmosphere_PlotData ( Name_Bathymetry_File, u_0, t_0, v, w, t, c, Precipitation, precipitable_water );



//	statements on the convergence und iterational process
	velocity_iter = 0;
	n--;

	if ( pressure_iter <= pressure_iter_max )
	{
		goto Pressure_loop;
	}

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of pressure loop: if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::



	// time_slice_change:


// reset of results to the initial value
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				precipitable_water.y[ j ][ k ] = 0.;
				Precipitation.y[ j ][ k ] = 0.;
				precipitation_NASA.y[ j ][ k ] = 0.;
				Evaporation_Penman.y[ j ][ k ] = 0.;

				for ( int i = 0; i < im; i++ )
				{
					h.x[ i ][ j ][ k ] = 0.;

					t.x[ i ][ j ][ k ] = 0.;
					c.x[ i ][ j ][ k ] = 0.;
					cloud.x[ i ][ j ][ k ] = 0.;
					ice.x[ i ][ j ][ k ] = 0.;
					co2.x[ i ][ j ][ k ] = 0.;
					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
					tn.x[ i ][ j ][ k ] = 0.;
					cn.x[ i ][ j ][ k ] = 0.;

					cloudn.x[ i ][ j ][ k ] = 0.;
					icen.x[ i ][ j ][ k ] = 0.;
					co2n.x[ i ][ j ][ k ] = 0.;
					un.x[ i ][ j ][ k ] = 0.;
					vn.x[ i ][ j ][ k ] = 0.;
					wn.x[ i ][ j ][ k ] = 0.;

					p_dyn.x[ i ][ j ][ k ] = 0.;
					p_stat.x[ i ][ j ][ k ] = 0.;

					P_rain.x[ i ][ j ][ k ] = 0.;
					P_snow.x[ i ][ j ][ k ] = 0.;

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
					aux_p.x[ i ][ j ][ k ] = 0.;

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





//	choice of the next time slice after nm iterations reached
	time = dt;
	velocity_iter = 1;
	pressure_iter = 1;
	velocity_iter_2D = 1;
	pressure_iter_2D = 1;
	switch_2D = 0;

	i_time_slice++;
	Ma = time_slice [ i_time_slice ];

	if ( i_time_slice >= i_time_slice_max ) goto finish;
	else goto time_slice_sequel;

//   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	finish:

// delete temporary arrays
	delete [ ] time_slice;
	delete [ ] im_tropopause;

// 	final remarks
	cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;

	if ( velocity_iter == velocity_iter_max )	cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of maximum time steps ***** \n\n" << endl;

	if ( min <= epsres )		cout << "***** steady solution reached! *****" << endl;

	cout << endl;
	cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
	cout << "\n\n\n\n";

	return 0;
}
