/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-hydrological circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with an additional transport equation to describe the salt concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <netcdf.h>
#include <vector>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "BC_Hyd.h"
#include "BC_Bath_Hyd.h"
#include "IC_Thermohalin.h"
#include "Accuracy_Hyd.h"
#include "RHS_Hyd.h"
#include "RungeKutta_Hyd.h"
#include "Print_Hyd.h"
#include "PostProcess_Hyd.h"
#include "Pressure_Hyd.h"
#include "Restore_Hyd.h"
#include "MinMax_Hyd.h"
#include "Results_Hyd.h"
#include "File_NetCDF.h"


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

// "salt"							Computation with and without a concentration equation for water vapor
// "salt 1"						with inclusion of water vapour
// "salt 0"						without inclusion of water vapour


// Earth's radius is r_earth = 6731 km compares to 6.731 [ / ]
// for 6 km expansion of the area of circulation compares to 0.02 [ / ] with 40 steps of size 0.0005 

// Definition of meridional and longitudinal step sizes 
// for example: dthe = the_Grad / pi180 = 1.125 / 57.3 = 0.01963

// maximum velocity on the sea surface  w_max = 0.29 [ / ] compares to 0.21 m/s = 0.78 km/h as annual mean 

// maximum temperature of earth's surface at equator t_max = 1.1355 compares to 37° C compares to 310 K
// maximum temperature of earth's surface at equator t_max = 1.0974 compares to 27° C compares to 300 K
// minimum temperature at the poles t_pol = .7803 compares to -60° C compares to 213.15 K
// minimum temperature in the deep ocean t_deep_ocean = 1.0146 compares to 4° C compares to 277.15 K
// temperature t_0 = 1.000 compares to 0° C compares to 273,15 K
// temperature t_0 = 1.073 compares to 20° C compares to 293,15 K
// temperature t_0 = 0.003661 compares to 1° C compares to 1 K

// mean velocity for sea level is 0.5 to 1 km/h compares to 0.14 to 0.28 m/s

// mass of water compares to 1.0, rate of salt compares to 0.035,
// c_0 compares to the total mass for mean salinity of 35.0 psu or dimensionsless 1.
// for c = 0.8857 compares to a salinity of 31.0 psu
// for c = 0.9686 compares to a salinity of 33.9 psu
// for c = 1.0000 compares to a salinity of 35.0 psu
// for c = 1.0571 compares to a salinity of 37.0 psu




int main ( int argc, char *argv[ ] )
{
// maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km ), 
// maximum number of overall iterations ( n ),
// maximum number of inner velocity loop iterations ( velocity_iter_max ),
// maximum number of outer pressure loop iterations ( pressure_iter_max )

	int im = 41, jm = 181, km = 361, nm = 200, velocity_iter_max = 2, pressure_iter_max = 2;
	int velocity_iter_max_2D = 10, pressure_iter_max_2D =5;

	int n, i_radial, j_longal, k_zonal;
	int velocity_iter, pressure_iter, pressure_iter_aux, velocity_iter_2D, pressure_iter_2D;
	int i_res, j_res, k_res;
	int i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_p, j_p, k_p;
	int Ma, i_time_slice, i_time_slice_max;
	int switch_2D;

	double time, residuum_old, residuum;
	double min, min_u, min_v, min_w, min_t, min_c, min_p;

	const int SequelFile = 0;									// sequel file will be written
	const double coriolis = 1.;									// computation with Coriolis force
	const double centrifugal = 1.;							// computation with centrifugal force
	const double salt = 1.;										// computation with salt concentration

	const int Ma_max = 300;									// parabolic temperature distribution 300 Ma back
	const int Ma_max_half = 150;							// half of time scale

	const double L_hyd = 6000.;							// extension of the hydrosphere shell in m
	const double dt = 0.0001;									// time step coincides with the CFL condition
	const double dr = 0.0005;									// compares to 150m depth
	const double the_Grad = 1.;								// compares to 1° step size laterally
	const double phi_Grad = 1.;								// compares to 1° step size longitudinally
	const double pi180 = 180./M_PI;						// pi180 = 57.3
	const double dthe = the_Grad / pi180, dphi = phi_Grad / pi180;//dthe = the_Grad / pi180 = 1.125 / 57.3 = 0.01963

	const double re = 1000.;									// Reynolds number
	const double ec = .0001;									// Eckert number
	const double sc = 10.;										// Schmidt number for salt water
	const double pr = 6.957;									// Prandtl number for water
	const double gr = 9.8066;									// gravitational acceleration of the earth
	const double cp_w = 4182.;								// specific heat capacity of water at constant pressure and 20°C in J/( kg K )
	const double omega = 7.29e-5;						// rotation number of the earth
	const double p_0 = 1013.25;							// pressure at sea level
	const double t_0 = 273.15;								// temperature in K compares to 0°C
	const double c_0 = 35.0;									// rate of salt in psu as maximum value
	const double u_0 = 0.45;									// maximum value of velocity in m/s
	const double r_0_water = 1026.0;					// density of salt water in kg/m3
//	const double mue_air = 17.1;							// dynamic viscosity of air in muePa * s at 20°C
//	const double mue_water = 1000.;					// dynamic viscosity of water in muePa * s at 20°C
	const double epsres = 0.0005;							// accuracy for relative and absolute errors

	const double ua = 0.;										// initial velocity component in r-direction
	const double va = 0.;										// initial velocity component in theta-direction
	const double wa = 0.;										// initial velocity component in phi-direction
	const double pa = 0.;										// initial value for the pressure field
	const double ta = 1.0146;									// compares to 4°C, threshhold temperature for the Boussinesq-approximation concerning bouyancy effect
	const double ca = 0.9686;								// c = 0.9686 compares to a salinity of 33.9 psu, mean value, ca corresponds to ta
	const double t_Cretaceous_max = 10.;			// maximum add of mean temperature during cretaceous
	const double c_Boussinesq = 1.;						// compares to 20°C, threshhold temperature for the Boussinesq-approximation concerning bouyancy effect
	const double r0 = 6.731; 									// Earth's radius is r_earth = 6731 km compares to 6.731 [ / ]
	const double the0 = 0.;										// North Pole
	const double phi0 = 0.;										// zero meridian in Greenwich
	const double t_Average = 15.;							// mean temperature of the modern earth
	const double t_equator = 1.1355;						// temperature t_0 = 1.1355 compares to 37° C compares to 310 K
	const double t_pole = .7803;							// temperature at the poles t_pole = 0.8 compares to -60°C compares to 213.15 K


// 	classe Array for 1-D, 2-D und 3-D field declarations

// time slices to be run after actualizing 
	i_time_slice_max = 15;
	int *time_slice = new int [ i_time_slice_max ]; 	// time slices in Ma

	time_slice [ 0 ] = 0;							// Golonka-Bathymetry
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



	Array_1D	rad ( im, 0., r0, dr );							// radial coordinate direction
	Array_1D	the ( jm, 0., the0, dthe );					// lateral coordinate direction
	Array_1D	phi ( km, 0., phi0, dphi );					// longitudinal coordinate direction

	Array_2D	Upwelling ( jm, km, 0. );					// upwelling
	Array_2D	Downwelling ( jm, km, 0. );				// downwelling
	Array_2D	SaltFinger ( jm, km, 0. );					// salt bulge of higher density
	Array_2D	SaltDiffusion ( jm, km, 0. );				// salt diffusion of lower density
	Array_2D	Salt_total ( jm, km, 0. );					// rate of salt summed up in a vertical column
	Array_2D	BottomWater ( jm, km, 0. );				// BottomWater
	Array_2D	t_j ( jm, km, ta );								// auxilliar field for temperature at the sea surface
	Array_2D	c_j ( jm, km, 0. );								// auxilliar field for salinity at the sea surface
	Array_2D	p_j ( jm, km, 0. );								// auxilliar field for pressure  at the surface
	Array_2D	aux_2D_v ( jm, km, 0. );					// auxilliar field v
	Array_2D	aux_2D_w ( jm, km, 0. );					// auxilliar field w

	Array	t ( im, jm, km, ta );									// temperature
	Array	u ( im, jm, km, ua );									// u-component velocity component in r-direction
	Array	v ( im, jm, km, va );									// v-component velocity component in theta-direction
	Array	w ( im, jm, km, wa );								// w-component velocity component in phi-direction
	Array	p ( im, jm, km, pa );									// pressure
	Array	c ( im, jm, km, ca );									// salt
	Array	tn ( im, jm, km, ta );								// temperature new
	Array	un ( im, jm, km, ua );								// u-component velocity component in r-direction new
	Array	vn ( im, jm, km, va );								// v-component velocity component in theta-direction new
	Array	wn ( im, jm, km, wa );								// w-component velocity component in phi-direction new
	Array	pn ( im, jm, km, pa );								// pressure new
	Array	cn ( im, jm, km, ca );								// salt new
	Array	h ( im, jm, km, 0. );									// bathymetry, depth from sea level
	Array	rhs_t ( im, jm, km, 0. );								// auxilliar field RHS temperature
	Array	rhs_u ( im, jm, km, 0. );								// auxilliar field RHS
	Array	rhs_v ( im, jm, km, 0. );								// auxilliar field RHS
	Array	rhs_w ( im, jm, km, 0. );								// auxilliar field RHS
	Array	rhs_p ( im, jm, km, 0. );								// auxilliar field RHS pressure
	Array	rhs_c ( im, jm, km, 0. );								// auxilliar field RHS salt
	Array	aux_u ( im, jm, km, 0. );							// auxilliar field u
	Array	aux_v ( im, jm, km, 0. );							// auxilliar field v
	Array	aux_w ( im, jm, km, 0. );							// auxilliar field w
	Array	Salt_Finger ( im, jm, km, 0. );					// salt bulge of higher density
	Array	Salt_Diffusion ( im, jm, km, 0. );				// salt diffusion of lower densityl
	Array	Salt_Balance ( im, jm, km, 0. );				// +/- salt balance


//	cout << " ***** printout of 3D-fields ***** " << endl;
//	p.printArray();
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

// time slice to start with is the modern world
	i_time_slice = 0;
	Ma = 0;


// naming the bathymetry file by the choice of the time slice by Ma and by author
	string Name_Bathymetry_File;
	stringstream My;

// naming a transfer file for the v- and w-vlocity component of the atmosphere at sea level
	string Name_v_w_Transfer_File;
	stringstream ssName_v_w_Transfer_File;

// naming a possibly available sequel file
	string Name_Sequel_File; 
	stringstream ssNameSequel;

// naming the output netCDF-file
	string Name_netCDF_File;
	stringstream ssNameNetCDF;

// naming a file to read the surface temperature of the modern world
	string Name_SurfaceTemperature_File; 
	stringstream ssNameSurfaceTemperature;
	ssNameSurfaceTemperature << "SurfaceTemperature.xyz";
	Name_SurfaceTemperature_File = ssNameSurfaceTemperature.str();

// naming a file to read the surface salinity of the modern world
	string Name_SurfaceSalinity_File; 
	stringstream ssNameSurfaceSalinity;
	ssNameSurfaceSalinity << "SurfaceSalinity.xyz";
	Name_SurfaceSalinity_File = ssNameSurfaceSalinity.str();


	cout << "\n\n\n\n";
	cout << " Ocean General Circulation Modell ( OGCM ) applied to laminar flow" << endl;
	cout << " Program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
	cout << " Finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
	cout << " with 1 additional transport equations to describe the salt concentration" << endl;
	cout << " 4. order Runge-Kutta scheme to solve 2. order differential equations" << endl << endl;

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
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: exists" << endl;
	}
	else
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: could not be read" << endl << endl;
		cout << endl << endl;
	}


// 	reading of the surface salinity file if available
	FILE *SurfaceSalinity;
	SurfaceSalinity = fopen ( Name_SurfaceSalinity_File.c_str(), "r" );

// 	if not available, prepare initial conditions, otherwise skip

	if ( SurfaceSalinity != NULL )
	{
		cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: exists" << endl << endl;
	}
	else
	{
		cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: could not be read" << endl << endl;
		cout << endl << endl;
	}




// class RHS_Hydrosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
	RHS_Hydrosphere		prepare ( im, jm, km, dt, dr, dthe, dphi, re, ec, sc, gr, pr, omega, coriolis, centrifugal, salt );

// class RungeKutta_Hydrosphere for the explicit solution of the Navier-Stokes equations
	RungeKutta_Hydrosphere		result ( n, im, jm, km, dt );

// class Pressure for the subsequent computation of the pressure by a separat Euler equation
	Pressure		startPressure ( im, jm, km, dr, dthe, dphi );

// class IC_Thermohalin for the initial and boundary conditions of the flow properties
	IC_Thermohalin		oceanflow ( im, jm, km );

// class BC_Hydrosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
	BC_Hydrosphere		boundary ( im, jm, km );

// class Results_MSL_Hyd to compute and show results on the mean sea level, MSL
	Restore		oldnew( im, jm, km );

// class Results_MSL_Hyd to compute and show results on the mean sea level, MSL
	Results_MSL_Hyd		calculate_MSL ( im, jm, km );

// class File_NetCDF to write results in the format of a netCDF-file
//	File_NetCDF		printoutNetCDF ( im, jm, km );

/*
//  constant hydrostatic pressure gradient
	double dpdr_hyd;
	dpdr_hyd = - r_0_water * gr * L_hyd / ( 100. * p_0 );
*/



// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



//	choice of the time slice to be computed
	time_slice_sequel:

// choice of the time slice by Ma and by author
// choice of the transfer file for the surface velocity components
	if ( Ma == 0 )
	{
		Name_Bathymetry_File = "0Ma_etopo.xyz";
		Name_v_w_Transfer_File = "[0Ma_etopo.xyz]_Transfer_Atm.vw";
		Name_Sequel_File = "[0Ma_etopo.xyz]_Sequel_Hyd.seq";
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

	ssName_v_w_Transfer_File << "[" << Name_Bathymetry_File << "]_Transfer_Atm.vw";
	Name_v_w_Transfer_File = ssName_v_w_Transfer_File.str();
	ssName_v_w_Transfer_File.str("");
	ssName_v_w_Transfer_File.ignore(ssName_v_w_Transfer_File.rdbuf()->in_avail());

	ssNameSequel << "[" << Name_Bathymetry_File << "]_Sequel_Hyd.seq";
	Name_Sequel_File = ssNameSequel.str();
	ssNameSequel.str("");
	ssNameSequel.ignore(ssNameSequel.rdbuf()->in_avail());

	ssNameNetCDF << Name_Bathymetry_File << "_atmosphere.nc";
	Name_netCDF_File = ssNameNetCDF.str();
	ssNameNetCDF.str("");
	ssNameNetCDF.ignore(ssNameNetCDF.rdbuf()->in_avail());


// 	reading a transfer file for the v-w-velocity components at sea level if available 
	FILE *Transfer_read;
	Transfer_read = fopen ( Name_v_w_Transfer_File.c_str(), "r" );

// 	reading of the sequel file if available
	FILE *Hydrosphere_read;
	Hydrosphere_read = fopen ( Name_Sequel_File.c_str(), "r" );


// 	if not available, skip
	if ( Transfer_read == NULL )
	{
		cout << "***** file ::::: " << Name_v_w_Transfer_File << " ::::: not yet exists" << endl << endl;
	}
	else
	{
//	class PostProcess for data transport, read and write
		PostProcess_Hydrosphere		read_Transfer ( im, jm, km );
		read_Transfer.Atmosphere_TransferFile_read ( Name_Bathymetry_File, v, w, p );

		cout << "***** file ::::: " << Name_v_w_Transfer_File << " ::::: could be read" << endl << endl;
	}


// 	if not available, skip
	if ( Hydrosphere_read == NULL )
	{
		cout << "***** file ::::: " << Name_Sequel_File << " ::::: not yet exists" << endl << endl;
	}
	else
	{
//	class PostProcess for data transport, read and write
		PostProcess_Hydrosphere		read_File ( im, jm, km );
		read_File.Hydrosphere_SequelFile_read ( Name_Bathymetry_File, n, pressure_iter, time, rad, the, phi, h, t, u, v, w, c, tn, un, vn, wn, cn, aux_u, aux_v, aux_w, t_j, c_j );

		cout << "***** file ::::: " << Name_Sequel_File << " ::::: could be read" << endl << endl;
		cout << "\n\n***** Hydrosphere_SequelFile_read in OGCM_main.cpp:   n = " << n << "   iter_BC = " << pressure_iter << "   time = " << time << endl << endl << endl;
		cout << endl << endl;
	}


	cout << endl << endl;
	cout << "***** time slice for the Oceanic Global Circulation Modell ( OGCM ) is:    Ma = " << Ma << " million years" << endl << endl;
	cout << "***** bathymetry/topography given by the x-y-z data set:    " << Name_Bathymetry_File.c_str() << endl << endl;


// 	reading of the bathymetry file if available
	FILE *Hydrosphere_Bathymetry;
	Hydrosphere_Bathymetry = fopen ( Name_Bathymetry_File.c_str(), "r" );

// 	if not available, prepare initial conditions, otherwise skip

	if ( Hydrosphere_Bathymetry != NULL )
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: exists" << endl << endl;
	}
	else
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: could not be read" << endl << endl;
		cout << endl << endl;
	}



//	initial conditions for u-v-w-velocity components and salinity

// class BC_Bathymetry_Hydrosphere for the geometrical boundary condition of the computational area
	BC_Bathymetry_Hydrosphere		depth ( im, jm, km );

// 	class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
	depth.BC_SeaGround ( Name_Bathymetry_File, h, aux_w );

//	surface temperature from World Ocean Atlas 2009 given as boundary condition
	if ( Ma == 0 ) oceanflow.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t_j, t );

//	surface salinity from World Ocean Atlas 2009 given as boundary condition
	if ( Ma == 0 ) oceanflow.BC_Surface_Salinity ( Name_SurfaceSalinity_File, c_j, c );

//	surface pressure computed by surface temperature with gas equation
//	oceanflow.BC_Surface_Pressure ( pa, gr, r_0_water, p_0, t_0, p_j, t_j, h, p, t );

//	salinity distribution as initial condition in 3 dimensions
	oceanflow.BC_Temperature_Salinity ( Ma, Ma_max, Ma_max_half, t_0, p_0, c_0, t_Cretaceous_max, t_Average, t_equator, t_pole, ua, va, wa, ta, ca, pa, t_j, c_j, p_j, h, t, c, p );

// Restore from new to old values to fill the gaps from time slice to time slice
//	if ( Ma > 0 ) depth.IC_SeaGroundGaps ( h, u, v, w, t, p, c, un, vn, wn, tn, pn, cn );

//	import of surface v- and w-velocity components from AGCM, surface velocity reduced to 3% of the wind velocity
	oceanflow.IC_v_w_Atmosphere ( h, u, v, w );

//	initial conditions for v and w velocity components at the sea surface, reduction of velocity with increasing depth for the purpose of the Ekman spiral
//	oceanflow.IC_v_w_Ekman ( h, v, w );

//	initial conditions for v and w velocity components at the sea surface close to east or west coasts, to close gyres
	oceanflow.IC_v_w_WestEastCoast ( h, u, v, w, c );

//	initial conditions for u-, v- and w-velocity components in deep flows, assumptions for thermohaline transveyor belt
//	if ( Ma == 0 ) oceanflow.IC_DeepWater ( h, u, v, w, c );
//	oceanflow.IC_DeepWater ( h, u, v, w, c );

//	currents along the equator
//	if ( Ma == 0 ) oceanflow.IC_EquatorialCurrents ( h, u, v, w );
//	oceanflow.IC_EquatorialCurrents ( h, u, v, w );

//	antarctic circumpolar current
//	if ( Ma == 0 ) oceanflow.IC_South_Polar_Sea ( h, u, v, w, c );
//	oceanflow.IC_South_Polar_Sea ( h, u, v, w, c );

//	arctic currents
//	if ( Ma == 0 ) oceanflow.IC_Nord_Polar_Meer ( h, u, v, w );
//	oceanflow.IC_Nord_Polar_Meer ( h, u, v, w );

//	Atlantic ocean currents
//	if ( Ma == 0 ) oceanflow.IC_Atlantischer_Ozean ( h, u, v, w, c );

//	Pacific ocean currents
//	if ( Ma == 0 ) oceanflow.IC_Pazifischer_Ozean ( h, u, v, w );

//	Indic ocean currents
//	if ( Ma == 0 ) oceanflow.IC_Indischer_Ozean ( h, u, v, w );

//	initial conditions for v and w velocity components at the sea surface close to east or west coasts, to close gyres
//	if ( Ma == 0 ) oceanflow.IC_v_w_Smoothing ( h, v, w, t, c );
//	oceanflow.IC_v_w_Smoothing ( iter, h, u, v, w, t, c );



//	storing of velocity components, pressure and temperature for iteration start
//	if ( n > 0 )   	oldnew.restoreOldNew ( .99, u, v, w, t, p, c, un, vn, wn, tn, pn, cn );
//	else    	  oldnew.restoreOldNew_2D ( .99, v, w, p, vn, wn, pn );
	oldnew.restoreOldNew ( .99, u, v, w, t, p, c, un, vn, wn, tn, pn, cn );
	oldnew.restoreOldNew_2D ( .99, v, w, p, vn, wn, pn );

	cout << " |||||||||||||||||||||||||||||||||      after restoreOldNew, before begin of pressure loop     ||||||||||||||||||||||||||||||||||||| " << endl;
	cout << " ***** printout of 3D-fields ***** " << endl;
	w.printArray();


// computation of the ratio ocean to land areas
	calculate_MSL.land_oceanFraction ( h );



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

	cout << " |||||||||||||||||||||||||||||||||      begin of velocity loop_2D     ||||||||||||||||||||||||||||||||||||| " << endl;
	cout << " ***** printout of 3D-fields ***** " << endl;
	w.printArray();


//		class BC_Atmosphaere for the geometry of a shell of a sphere
			boundary.RB_theta ( ca, ta, pa, t, u, v, w, p, c, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, aux_u, aux_v, aux_w, h, Salt_Finger, Salt_Diffusion, Salt_Balance );
			boundary.RB_phi ( t, u, v, w, p, c, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, aux_u, aux_v, aux_w, h, Salt_Finger, Salt_Diffusion, Salt_Balance );
//			if ( velocity_iter_2D == 1 )  boundary.BC_NST_control_2D ( dr, dthe, dphi, re, mue_air, mue_water, h, v, w, p, aux_2D_v, aux_2D_w, rad, the );
//			boundary.BC_NST_control_2D ( dr, dthe, dphi, re, mue_air, mue_water, h, v, w, p, aux_2D_v, aux_2D_w, rad, the );

//		pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
			startPressure.computePressure_2D ( pa, rad, the, p, h, rhs_v, rhs_w, aux_v, aux_w );

// 		old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy		min_Residuum_old_2D ( n, im, jm, km, dr, dthe, dphi );
			min_Residuum_old_2D.residuumQuery_2D ( j_res, k_res, min, rad, the, v, w );
			residuum_old = min;

//		class RungeKutta for the solution of the differential equations describing the flow properties
			result.solveRungeKutta_2D_Hydrosphere ( prepare, rad, the, phi, rhs_v, rhs_w, h, v, w, p, vn, wn, aux_v, aux_w );

//		state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
			Accuracy		min_Stationary_2D ( n, im, jm, km, dr, dthe, dphi );
			min_Stationary_2D.steadyQuery_2D ( i_v, j_v, k_v, i_w, j_w, k_w, i_p, j_p, k_p, min_v, min_w, min_p, v, vn, w, wn, p, pn );

// 		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy		min_Residuum_2D ( n, im, jm, km, dr, dthe, dphi );
			min_Residuum_2D.residuumQuery_2D ( j_res, k_res, min, rad, the, v, w );
			residuum = min;
			min = fabs ( ( residuum - residuum_old ) / residuum_old );


//		statements on the convergence und iterational process
			Accuracy		printout_2D ( im, Ma, n, velocity_iter_2D, pressure_iter_2D, min, L_hyd );
			printout_2D.iterationPrintout_2D ( nm, velocity_iter_max_2D, pressure_iter_max_2D, j_res, k_res, i_v, j_v, k_v, i_w, j_w, k_w, i_p, j_p, k_p, min_v, min_w, min_p );

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

	cout << " |||||||||||||||||||||||||||||||||      end of velocity loop_2D     ||||||||||||||||||||||||||||||||||||| " << endl;
	cout << " ***** printout of 3D-fields ***** " << endl;
	w.printArray();






	process_3D:

		if ( min >= epsres )
		{
			time = time + dt;
		}

// 		old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
		Accuracy		min_Residuum_old ( n, im, jm, km, dr, dthe, dphi );
		min_Residuum_old.residuumQuery ( i_res, j_res, k_res, min, rad, the, u, v, w );
		residuum_old = min;

// 		class RB_Hydrosphaere for the geometry of a shell of a sphere
		boundary.RB_radius ( ca, ta, pa, dr, rad, t, u, v, w, p, c, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, aux_u, aux_v, aux_w, h, Salt_Finger, Salt_Diffusion, Salt_Balance );
		boundary.RB_theta ( ca, ta, pa, t, u, v, w, p, c, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, aux_u, aux_v, aux_w, h, Salt_Finger, Salt_Diffusion, Salt_Balance );
		boundary.RB_phi ( t, u, v, w, p, c, rhs_u, rhs_v, rhs_w, rhs_t, rhs_c, aux_u, aux_v, aux_w, h, Salt_Finger, Salt_Diffusion, Salt_Balance );
//		if ( velocity_iter == 1 )  boundary.BC_NST_control_3D ( dr, dthe, dphi, re, mue_air, mue_water, h, u, v, w, t, p, c, co2, aux_u, aux_v, aux_w, rad, the );
//		boundary.BC_NST_control_3D ( dr, dthe, dphi, re, mue_air, mue_water, h, u, v, w, t, p, c, aux_u, aux_v, aux_w, rad, the );

// 		class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
		depth.BC_SolidGround ( ca, ta, pa, h, t, u, v, w, p, c, tn, un, vn, wn, pn, cn, t_j, c_j, p_j );

// 		class RungeKutta for the solution of the differential equations describing the flow properties
		result.solveRungeKutta_Hydrosphere ( prepare, L_hyd, gr, cp_w, c_Boussinesq, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, t, u, v, w, p, c, tn, un, vn, wn, cn, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Salt_Balance );

//	initial conditions for v and w velocity components at the sea surface close to east or west coasts, to close gyres
//	if ( Ma == 0 ) oceanflow.IC_v_w_Smoothing ( h, v, w, t, c );
//	oceanflow.IC_v_w_Smoothing ( iter, h, u, v, w, t, c );


// 		class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
		depth.BC_SolidGround ( ca, ta, pa, h, t, u, v, w, p, c, tn, un, vn, wn, pn, cn, t_j, c_j, p_j );


//		state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
		Accuracy		min_Stationary ( n, im, jm, km, dr, dthe, dphi );
		min_Stationary.steadyQuery ( i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_p, j_p, k_p, min_u, min_v, min_w, min_t, min_c, min_p, u, un, v, vn, w, wn, t, tn, c, cn, p, pn );

		Accuracy		min_Residuum ( n, im, jm, km, dr, dthe, dphi );
		min_Residuum.residuumQuery ( i_res, j_res, k_res, min, rad, the, u, v, w );
		residuum = min;
		min = fabs ( ( residuum - residuum_old ) / residuum_old );

//		statements on the convergence und iterational process
		Accuracy		printout ( im, Ma, n, velocity_iter, pressure_iter, min, L_hyd );
		printout.iterationPrintout ( nm, velocity_iter_max, pressure_iter_max, i_res, j_res, k_res, i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_p, j_p, k_p, min_u, min_v, min_w, min_t, min_c, min_p );





//	searching of maximum and minimum values of salt concentration
		string str_max_salt_concentration = " max salt concentration ", str_min_salt_concentration = " min salt concentration ", str_unit_salt_concentration = "psu";
		MinMax		minmaxSalt ( im, jm, km, c_0, L_hyd );
		minmaxSalt.searchMinMax_3D ( str_max_salt_concentration, str_min_salt_concentration, str_unit_salt_concentration, c );


//	total salinity as the sum along a normally extended virtual column
//	in r-direction salinity above the average value is added
		calculate_MSL.run_MSL_data ( c_0, h, u, v, w, c, Salt_Finger, Salt_Diffusion, Upwelling, Downwelling, SaltFinger, SaltDiffusion, Salt_total, BottomWater );


//	searching of maximum and minimum values of total salt volume in a column
		string str_max_salt_total = " max salt total ", str_min_salt_total = " min salt total ", str_unit_salt_total = "psu";
		MinMax		minmaxSalt_total ( jm, km, c_0 );
		minmaxSalt_total.searchMinMax ( str_max_salt_total, str_min_salt_total, str_unit_salt_total, Salt_total );

//	searching of maximum and minimum values of salt finger volume in a column
		string str_max_salt_finger = " max salt finger ", str_min_salt_finger = " min salt finger ", str_unit_salt_finger = "psu";
		MinMax		minmaxSalt_finger ( jm, km, c_0 );
		minmaxSalt_finger.searchMinMax ( str_max_salt_finger, str_min_salt_finger, str_unit_salt_finger, SaltFinger );

//	searching of maximum and minimum values of salt diffusion volume in a column
		string str_max_salt_diffusion = " max salt diffusion ", str_min_salt_diffusion = " min salt diffusion ", str_unit_salt_diffusion = "psu";
		MinMax		minmaxSalt_diffusion ( jm, km, c_0 );
		minmaxSalt_diffusion.searchMinMax ( str_max_salt_diffusion, str_min_salt_diffusion, str_unit_salt_diffusion, SaltDiffusion );

//	searching of maximum and minimum values of upwelling volume in a column
		string str_max_upwelling = " max upwelling ", str_min_upwelling = " min upwelling ", str_unit_upwelling = "m/s";
		MinMax		minmaxUpwelling ( jm, km, c_0 );
		minmaxUpwelling.searchMinMax ( str_max_upwelling, str_min_upwelling, str_unit_upwelling, Upwelling );

//	searching of maximum and minimum values of upwelling volume in a column
		string str_max_downwelling = " max downwelling ", str_min_downwelling = " min downwelling ", str_unit_downwelling = "m/s";
		MinMax		minmaxDownwelling ( jm, km, c_0 );
		minmaxDownwelling.searchMinMax ( str_max_downwelling, str_min_downwelling, str_unit_downwelling, Downwelling );

//	searching of maximum and minimum values of bottom water volume in a column
		string str_max_bottom_water = " max bottom water ", str_min_bottom_water = " min bottom water ", str_unit_bottom_water = "m/s";
		MinMax		minmaxBottom_water ( jm, km, c_0 );
		minmaxBottom_water.searchMinMax ( str_max_bottom_water, str_min_bottom_water, str_unit_bottom_water, BottomWater );


//	total values of properties in a normally extended virtual column are added
		calculate_MSL.show_MSL_data ( c_0, h, c, t, p, u, Upwelling, Downwelling, BottomWater, SaltFinger, SaltDiffusion, Salt_total );


//	restoring the velocity component and the temperature for the new time step
		oldnew.restoreOldNew ( 1., u, v, w, t, p, c, un, vn, wn, tn, pn, cn );
	}


//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	cout << " |||||||||||||||||||||||||||||||||      end of loop: while, before printout    ||||||||||||||||||||||||||||||||||||| " << endl;
	cout << " ***** printout of 3D-fields ***** " << endl;
	w.printArray();



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

//	class Print_Hydrosphaere for the printout of results
	Print_Hydrosphere		printout ( im, jm, km, nm, n, time );

//	printout in ParaView files, netCDF files and sequel files

	pressure_iter_aux = pressure_iter - 1;

//	results written in netCDF format
//	printoutNetCDF.out_NetCDF( Name_netCDF_File, v, w, h, Upwelling, Downwelling, BottomWater );

//	class PostProcess_Hydrosphaere for the printing of results
	PostProcess_Hydrosphere		write_File ( im, jm, km );


	j_longal = 75;
	write_File.paraview_vtk_longal ( Name_Bathymetry_File, j_longal, pressure_iter_aux, h, p, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, Salt_Balance );

//	zonal data along constant longitudes
	k_zonal = 185;
	write_File.paraview_vtk_zonal ( Name_Bathymetry_File, k_zonal, pressure_iter_aux, h, p, t, u, v, w, c, Salt_Finger, Salt_Diffusion, Salt_Balance );

//	radial data along constant hight above ground
	i_radial = 40;
	write_File.paraview_vtk_radial ( Name_Bathymetry_File, i_radial, pressure_iter_aux, h, p, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, Salt_Balance, Upwelling, Downwelling, SaltFinger, SaltDiffusion, BottomWater );

//	3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
	write_File.paraview_panorama_vts ( Name_Bathymetry_File, pressure_iter_aux, h, t, p, u, v, w, c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Salt_Balance );

//	3-dimensional data in spherical coordinate system for a streamline pattern in a shell of a sphere
//	write_File.paraview_vts ( Name_Bathymetry_File, n, rad, the, phi, h, t, p, u, v, w, c, rhs_u, rhs_v, rhs_w, rhs_c, rhs_p, rhs_t, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Salt_Balance );


//	writing of sequential data for the sequel file
	if ( SequelFile == 1 )
	{
		write_File.Hydrosphere_SequelFile_write ( Name_Bathymetry_File, n, pressure_iter, time, rad, the, phi, h, t, u, v, w, c, tn, un, vn, wn, cn, aux_u, aux_v, aux_w, t_j, c_j  );
	}


//	writing of plot data in the PlotData file
	PostProcess_Hydrosphere		write_PlotData_File ( im, jm, km );
	write_PlotData_File.Hydrosphere_PlotData ( Name_Bathymetry_File, v, w, t, c, BottomWater, Upwelling, Downwelling );




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
	if ( i_time_slice >= i_time_slice_max ) goto finish;
	else goto time_slice_sequel;

//   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	finish:

// delete temporary arrays
	delete [ ] time_slice;

// 	final remarks
	cout << endl << "***** end of the Oceanic General Circulation Modell ( OGCM ) *****" << endl << endl;

	if ( velocity_iter == velocity_iter_max )	cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of maximum time steps ***** \n\n" << endl;

	if ( min <= epsres )		cout << "***** steady solution reached! *****" << endl;

	cout << endl;
	cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
	cout << "\n\n\n\n";

	return 0;
}
