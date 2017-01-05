/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * program for the computation of geo-atmospherical circulating flows in a spherical shell
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <vector>
#include <netcdf.h>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "BC_Hyd.h"
#include "BC_Bath_Hyd.h"
#include "BC_Thermohalin.h"
#include "Accuracy_Hyd.h"
#include "RHS_Hyd.h"
#include "RungeKutta_Hyd.h"
#include "Print_Hyd.h"
#include "PostProcess_Hyd.h"
#include "Pressure_Hyd.h"
#include "Restore_Hyd.h"
#include "MinMax_Hyd.h"
#include "Results_Hyd.h"
#include "File_NetCDF_Hyd.h"

#include "tinyxml2.h"

#include "cHydrosphereModel.h"
#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;

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


// Earth's radius is r_earth = 6731 km compares to 6.731 [ / ]
// for 6 km expansion of the area of circulation compares to 0.02 [ / ] with 40 steps of size 0.0005 

// Definition of meridional and longitudinal step sizes 
// for example: dthe = the_Grad / pi180 = 1.125 / 57.3 = 0.01963

// maximum velocity on the sea surface  w_max = 0.29 [ / ] compares to 0.21 m/s = 0.78 km/h as annual mean 
// mean velocity for sea level is 0.5 to 1 km/h compares to 0.14 to 0.28 m/s

// maximum temperature of earth's surface at equator t_max = 1.1355 compares to 37° C compares to 310 K
// maximum temperature of earth's surface at equator t_max = 1.0974 compares to 27° C compares to 300 K
// minimum temperature at the poles t_pol = .7803 compares to -60° C compares to 213.15 K
// minimum temperature in the deep ocean t_deep_ocean = 1.0146 compares to 4° C compares to 277.15 K
// temperature t_0 = 1.000 compares to 0° C compares to 273,15 K
// temperature t_0 = 1.073 compares to 20° C compares to 293,15 K
// temperature t_0 = 0.003661 compares to 1° C compares to 1 K

// 1 PSU ( Practical Salt Unit ) = 1 g/kg, means g of salt per kg sweet water
// mass of water compares to 1.0, rate of salt compares to 0.0346
// c_0 compares to the total mass for mean salinity of 34.6 psu or dimensionsless 1.
// for c = 0.9249 compares to a salinity of 32.0 psu
// for c = 0.9682 compares to a salinity of 33.5 psu
// for c = 1.0000 compares to a salinity of 34.6 psu
// for c = 1.0983 compares to a salinity of 38.0 psu




cHydrosphereModel::cHydrosphereModel() {
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.
    PythonStream::OverrideCout();

    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

   	// TODO: set default configuration
}

cHydrosphereModel::~cHydrosphereModel() { }

void cHydrosphereModel::LoadConfig(const char *filename) {
	// FIXME: not implemented
	cout << "FIXME: not implemented\n";
}

void cHydrosphereModel::Run() {
// maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km ), 
// maximum number of overall iterations ( n ),
// maximum number of inner velocity loop iterations ( velocity_iter_max ),
// maximum number of outer pressure loop iterations ( pressure_iter_max )

	int im = 41, jm = 181, km = 361, nm = 200, velocity_iter_max = 5, pressure_iter_max = 5;
	int velocity_iter_max_2D = 5, pressure_iter_max_2D = 5;
	int n = 0, i_radial = 0, j_longal = 0, k_zonal = 0;
	int velocity_iter = 0, pressure_iter = 0, pressure_iter_aux = 0, velocity_iter_2D = 0, pressure_iter_2D = 0;
	int i_res = 0, j_res = 0, k_res = 0;
	int Ma = 0, i_time_slice = 0, i_time_slice_max = 0, i_max = 0, i_beg = 0;
	int switch_2D = 0;

	double time = 0., residuum_old = 0., residuum = 0., min = 0.;

	const int SequelFile = 0;							// sequel file will be written
	double coriolis = 1.;									// computation with Coriolis force
	double centrifugal = 1.;							// computation with centrifugal force
	double buoyancy = 1.;								// computation with buoyancy

	int Ma_max = 300;									// parabolic temperature distribution 300 Ma back
	int Ma_max_half = 150;								// half of time scale

//	double L_hyd = 10000.;							// extension of the hydrosphere shell in m, assumption of maximum depth of sea 10000 m compares to 40 steps times 250 m
	double L_hyd = 6000.;							// extension of the hydrosphere shell in m, assumption of maximum depth of sea 6000 m compares to 40 steps times 150 m

	double dt = 0.0001;									// time step coincides with the CFL condition
	double dr = 0.0005;									// compares to 150m depth
	double the_Grad = 1.;								// compares to 1° step size laterally
	double phi_Grad = 1.;								// compares to 1° step size longitudinally
	double pi180 = 180./M_PI;						// pi180 = 57.3
	double dthe = the_Grad / pi180, dphi = phi_Grad / pi180;//dthe = the_Grad / pi180 = 1.125 / 57.3 = 0.01963

	double re = 1000.;									// Reynolds numberc_0
	double ec = .0001;									// Eckert number
	double sc = 10.;										// Schmidt number for salt water
	double pr = 6.957;									// Prandtl number for water
	double g = 9.8066;									// gravitational acceleration of the earth
	double cp_w = 4182.;								// specific heat capacity of water at constant pressure and 20°C in J/( kg K )
	double omega = 7.29e-5;							// rotation number of the earth
	double p_0 = 1013.25;								// pressure at sea level in hPa
	double t_0 = 273.15;								// temperature in K compares to 0°C
	double c_0 = 34.6;									// rate of salt in psu at temperature t_0
	double u_0 = 0.45;									// maximum value of velocity in m/s
	double r_0_water = 1026.0;						// reference density of salt water in kg/m3
	double epsres = 0.0005;							// accuracy for relative and absolute errors0,988571429

	double ua = 0.;										// initial velocity component in r-direction
	double va = 0.;										// initial velocity component in theta-direction
	double wa = 0.;										// initial velocity component in phi-direction
	double pa = 0.;										// initial value for the pressure field
	double ta = 1.01464;								// compares to 4°C
	double ca = 1.01156;								// c = 1.01156 compares to a salinity of 35.0 psu, mean value, ca corresponds to ta = 1.01464  ( = 4°C )
	double ca_max = 1.0983;							// c = 1.0983 compares to a salinity of 38.00 psu  used for deep flow initialization
	double t_cretaceous_max = 10.;				// maximum add of mean temperature during cretaceous
//	double r0 = 6.731; 									// Earth's radius is r_earth = 6731 km compares to 6.731 [ / ]
	double r0 = 1.; 										// Earth's radius is r_earth = 6731 km compares to 6.731 [ / ]
	double the0 = 0.;										// North Pole
	double phi0 = 0.;										// zero meridian in Greenwich
	double t_average = 15.;							// mean temperature of the modern earth
	double t_equator = 1.1355;						// temperature t_0 = 1.1355 compares to 37° C compares to 310 K
	double t_pole = 1.0146;							// compares to 4°C, threshhold temperature for the Boussinesq-approximation concerning bouyancy effect


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


// 	class Array for 1-D, 2-D and 3-D field declarations

// 1D arrays
	Array_1D		rad;																		// radial coordinate direction
	rad.initArray_1D ( im, 1. );
	Array_1D		the;																		// lateral coordinate direction
	the.initArray_1D ( jm, 2. );
	Array_1D		phi;																		// longitudinal coordinate direction
	phi.initArray_1D ( km, 3. );


// 2D arrays
	Array_2D		Upwelling;															// upwelling
	Upwelling.initArray_2D ( jm, km, 0. );
	Array_2D		Downwelling;														// downwelling
	Downwelling.initArray_2D ( jm, km, 0. );
	Array_2D		BottomWater;														// 2D bottom water summed up in a vertical column
	BottomWater.initArray_2D ( jm, km, 0. );

	Array_2D		SaltFinger;															// salt bulge of higher density
	SaltFinger.initArray_2D ( jm, km, 0. );
	Array_2D		SaltDiffusion;														// salt bulge of lower density
	SaltDiffusion.initArray_2D ( jm, km, 0. );
	Array_2D		Salt_total;															// rate of salt summed up in a vertical column
	Salt_total.initArray_2D ( jm, km, 0. );

	Array_2D		BuoyancyForce_2D;											// radiation balance at the surface
	BuoyancyForce_2D.initArray_2D ( jm, km, 0. );


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

	Array	tn;																				// temperature new
	tn.initArray ( im, jm, km, ta );
	Array	un;																				// u-velocity component in r-direction new
	un.initArray ( im, jm, km, ua );
	Array	vn;																				// v-velocity component in theta-direction new
	vn.initArray ( im, jm, km, va );
	Array	wn;																				// w-velocity component in phi-direction new
	wn.initArray ( im, jm, km, wa );
	Array	cn;																				// water vapour new
	cn.initArray ( im, jm, km, ca );

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

	Array	aux_u;																			// auxilliar field u-velocity component
	aux_u.initArray ( im, jm, km, 0. );
	Array	aux_v;																			// auxilliar field v-velocity component
	aux_v.initArray ( im, jm, km, 0. );
	Array	aux_w;																		// auxilliar field w-velocity component
	aux_w.initArray ( im, jm, km, 0. );
	Array	aux_p;																			// auxilliar field p
	aux_p.initArray ( im, jm, km, pa );

	Array	Salt_Finger;																	// salt bulge of higher density
	Salt_Finger.initArray ( im, jm, km, 0. );
	Array	Salt_Diffusion;															// salt bulge of lowerer density
	Salt_Diffusion.initArray ( im, jm, km, 0. );
	Array	Salt_Balance;																// +/- salt balance
	Salt_Balance.initArray ( im, jm, km, 0. );

	Array	BuoyancyForce_3D;													// 3D buoyancy force
	BuoyancyForce_3D.initArray ( im, jm, km, 0. );


//	cout << " ***** printout of 3D-fields ***** " << endl;
//	t.printArray( im, jm, km );

//	cout << " ***** printout of 2D-fields ***** " << endl;
//	BottomWatert.printArray_2D( jm, km );

//	cout << " ***** printout of 1D-fields ***** " << endl;
//	rad.printArray_1D( im );

	cout.precision ( 6 );
	cout.setf ( ios::fixed );

//	Coordinate system in form of a spherical shell
//	rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
	rad.Coordinates ( im, r0, dr );
	the.Coordinates ( jm, the0, dthe );
	phi.Coordinates ( km, phi0, dphi );

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
	i_max = im - 1;																		// corresponds to sea level
//	i_beg = 27;																				// compares to an ocean depth of 1950 m, location of the thermocline
	i_beg = 33;																				// compares to an ocean depth of 1050 m, location of the thermocline

// time slice to start with is the modern world
	Ma = 0;
	i_time_slice = 0;


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


	cout << endl << endl << endl;
	cout << "***** Hydrosphere General Circulation Model ( OGCM ) applied to laminar flow" << endl;
	cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
	cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
	cout << "***** with 1 additional transport equations to describe the salinity" << endl;
	cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
	cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
	cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
	cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
	cout << "***** salinity is part of the Boussinesq approximation" << endl;
	cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )" << endl << endl;

	cout << "***** original program name:  " << __FILE__ << endl;
	cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;





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
		n++;
	}



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
		read_Transfer.Atmosphere_TransferFile_read ( Name_Bathymetry_File, v, w, p_dyn );

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
		read_File.Hydrosphere_SequelFile_read ( Name_Bathymetry_File, n, pressure_iter, time, rad, the, phi, h, t, u, v, w, c, tn, un, vn, wn, cn, aux_u, aux_v, aux_w );

		cout << "***** file ::::: " << Name_Sequel_File << " ::::: could be read" << endl << endl;
		cout << "\n\n***** Hydrosphere_SequelFile_read in OGCM_main.cpp:   n = " << n << "   iter_BC = " << pressure_iter << "   time = " << time << endl << endl << endl;
		cout << endl << endl;
		n++;
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
		n++;
	}



//	initial conditions for u-v-w-velocity components and salinity

// class BC_Bathymetry_Hydrosphere for the geometrical boundary condition of the computational area
	BC_Bathymetry_Hydrosphere		depth ( im, jm, km );

// 	class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
	depth.BC_SeaGround ( Name_Bathymetry_File, L_hyd, h, aux_w );



// class BC_Hydrosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
	BC_Hydrosphere		boundary ( im, jm, km );

// class RHS_Hydrosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
	RHS_Hydrosphere		prepare ( im, jm, km, r0, dt, dr, dthe, dphi, re, ec, sc, g, pr, omega, coriolis, centrifugal, buoyancy );

// class RungeKutta_Hydrosphere for the explicit solution of the Navier-Stokes equations
	RungeKutta_Hydrosphere		result ( n, im, jm, km, dt );

// class Pressure for the subsequent computation of the pressure by a separat Euler equation
	Pressure_Hyd		startPressure ( im, jm, km, dr, dthe, dphi );

// class Restore_Hyd to restore between iterational results
	Restore_Hyd		oldnew( im, jm, km );

// class Results_MSL_Hyd to compute and show results on the mean sea level, MSL
	Results_Hyd		calculate_MSL ( im, jm, km );

// class File_NetCDF to write results in the format of a netCDF-file
//	File_NetCDF		printoutNetCDF ( im, jm, km );



// class BC_Thermohalin for the initial and boundary conditions of the flow properties
	BC_Thermohalin		oceanflow ( im, jm, km, i_beg, i_max, Ma, Ma_max, Ma_max_half, dr, g, r_0_water, ua, va, wa, ta, ca, ca_max, pa, u_0, p_0, t_0, c_0, cp_w, L_hyd, t_average, t_cretaceous_max, t_equator, t_pole );

//	surface temperature from World Ocean Atlas 2009 given as boundary condition
//	if ( Ma == 0 ) oceanflow.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );

//	surface salinity from World Ocean Atlas 2009 given as boundary condition
//	if ( Ma == 0 ) oceanflow.BC_Surface_Salinity ( Name_SurfaceSalinity_File, c );

//	surface pressure computed by surface temperature with gas equation
//	oceanflow.BC_Surface_Pressure ( pa, g, r_0_water, p_0, t_0, p_j, t_j, h, p_dyn, t );

//	import of surface v- and w-velocity components from AGCM, surface velocity reduced to 3% of the wind velocity
	oceanflow.IC_v_w_Atmosphere ( h, u, v, w );

//	salinity distribution as initial condition in 3 dimensions
	oceanflow.BC_Temperature_Salinity ( h, t, c, p_dyn );


//  surface pressure computed by surface temperature with gas equation
	oceanflow.BC_Pressure ( p_stat, t, h );

//	initial conditions for v and w velocity components at the sea surface, reduction of velocity with increasing depth for the purpose of the Ekman spiral
//	oceanflow.IC_v_w_Ekman ( h, v, w );

//	initial conditions for v and w velocity components at the sea surface close to east or west coasts, to close gyres
	oceanflow.IC_v_w_WestEastCoast ( h, u, v, w, c );

//	initial conditions for u-, v- and w-velocity components in deep flows, assumptions for thermohaline transveyor belt
//	if ( Ma == 0 ) oceanflow.IC_DeepWater ( h, u, v, w, c );
//	oceanflow.IC_DeepWater ( h, u, v, w, c );

//	currents along the equator
//	if ( Ma == 0 ) oceanflow.IC_EquatorialCurrents ( h, u, v, w );
	oceanflow.IC_EquatorialCurrents ( h, u, v, w );

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


//	goto Print_commands;														// only initial conditions wihtout iterations


//	storing of velocity components, pressure and temperature for iteration start
	oldnew.restoreOldNew_3D ( .9, u, v, w, t, p_dyn, c, un, vn, wn, tn, aux_p, cn );
	oldnew.restoreOldNew_2D ( .9, v, w, p_dyn, vn, wn, aux_p );

// computation of the ratio ocean to land areas
	calculate_MSL.land_oceanFraction ( h );


// ******************************************   start of pressure and velocity iterations ************************************************************************


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of pressure loop : if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::

Pressure_loop:

// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of velocity loop: while ( min >= epsres )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	min = epsres * 20.;
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
			boundary.RB_theta ( ca, ta, pa, t, u, v, w, p_dyn, c );
			boundary.RB_phi ( t, u, v, w, p_dyn, c );

// 		old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy_Hyd		min_Residuum_old_2D ( im, jm, km, dthe, dphi );
			min_Residuum_old_2D.residuumQuery_2D ( rad, the, v, w );
			min = min_Residuum_old_2D.out_min (  );

			residuum_old = min;

//		class RungeKutta for the solution of the differential equations describing the flow properties
			result.solveRungeKutta_2D_Hydrosphere ( prepare, rad, the, phi, rhs_v, rhs_w, h, v, w, p_dyn, vn, wn, aux_v, aux_w );

// 		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
			Accuracy_Hyd		min_Residuum_2D ( im, jm, km, dthe, dphi );
			min_Residuum_2D.residuumQuery_2D ( rad, the, v, w );
			min = min_Residuum_2D.out_min (  );
			j_res = min_Residuum_2D.out_j_res (  );
			k_res = min_Residuum_2D.out_k_res (  );

			residuum = min;
			min = fabs ( ( residuum - residuum_old ) / residuum_old );


//		state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
			Accuracy_Hyd		min_Stationary_2D ( n, nm, Ma, im, jm, km, min, j_res, k_res, velocity_iter_2D, pressure_iter_2D, velocity_iter_max_2D, pressure_iter_max_2D );
			min_Stationary_2D.steadyQuery_2D ( h, v, vn, w, wn, p_dyn, aux_p );

			oldnew.restoreOldNew_2D ( 1., v, w, p_dyn, vn, wn, aux_p );

		}

 
//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop_2D: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Pressure_iteration_2D:

//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
	startPressure.computePressure_2D ( pa, rad, the, p_dyn, aux_p, h, rhs_v, rhs_w, aux_v, aux_w );		// 2D pressure computation causes a pressure jump in radial direction along coast lines, 3D treatment needed later, 2D velocities are though corrected

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
		Accuracy_Hyd		min_Residuum_old ( im, jm, km, dr, dthe, dphi );
		min_Residuum_old.residuumQuery_3D ( rad, the, u, v, w );
		min = min_Residuum_old.out_min (  );

		residuum_old = min;

// 		class RB_Hydrosphaere for the geometry of a shell of a sphere
		boundary.RB_radius ( ca, ta, pa, dr, rad, t, u, v, w, p_dyn, c );
		boundary.RB_theta ( ca, ta, pa, t, u, v, w, p_dyn, c );
		boundary.RB_phi ( t, u, v, w, p_dyn, c );

// 		class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
		depth.BC_SolidGround ( ca, ta, pa, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, aux_p, cn );

// 		class RungeKutta for the solution of the differential equations describing the flow properties
		result.solveRungeKutta_3D_Hydrosphere ( prepare, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, t, u, v, w, p_dyn, c, tn, un, vn, wn, cn, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

// 		class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
		depth.BC_SolidGround ( ca, ta, pa, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, aux_p, cn );

// 		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
		Accuracy_Hyd		min_Residuum ( im, jm, km, dr, dthe, dphi );
		min_Residuum.residuumQuery_3D ( rad, the, u, v, w );
		min = min_Residuum.out_min (  );
		i_res = min_Residuum.out_i_res (  );
		j_res = min_Residuum.out_j_res (  );
		k_res = min_Residuum.out_k_res (  );

		residuum = min;
		min = fabs ( ( residuum - residuum_old ) / residuum_old );

//		statements on the convergence und iterational process
		Accuracy_Hyd		min_Stationary ( n, nm, Ma, im, jm, km, min, i_res, j_res, k_res, velocity_iter, pressure_iter, velocity_iter_max, pressure_iter_max, L_hyd );
		min_Stationary.steadyQuery_3D ( u, un, v, vn, w, wn, t, tn, c, cn, p_dyn, aux_p );



// 3D_fields

//		searching of maximum and minimum values of temperature
		string str_max_temperature = " max temperature ", str_min_temperature = " min temperature ", str_unit_temperature = "C";
		MinMax_Hyd		minmaxTemperature ( im, jm, km, c_0, L_hyd );
		minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

//		searching of maximum and minimum values of pressure
		string str_max_pressure = " max pressure dynamic ", str_min_pressure = " min pressure dynamic ", str_unit_pressure = "hPa";
		MinMax_Hyd		minmaxPressure ( im, jm, km, c_0, L_hyd );
		minmaxPressure.searchMinMax_3D ( str_max_pressure, str_min_pressure, str_unit_pressure, p_dyn, h );

//		searching of maximum and minimum values of static pressure
		string str_max_pressure_stat = " max pressure static ", str_min_pressure_stat = " min pressure static ", str_unit_pressure_stat = "bar";
		MinMax_Hyd		minmaxPressure_stat ( im, jm, km, c_0, L_hyd );
		minmaxPressure_stat.searchMinMax_3D ( str_max_pressure_stat, str_min_pressure_stat, str_unit_pressure_stat, p_stat, h );

		cout << endl << " salinity based results in the three dimensional space: " << endl << endl;

//	searching of maximum and minimum values of salt concentration
		string str_max_salt_concentration = " max salt concentration ", str_min_salt_concentration = " min salt concentration ", str_unit_salt_concentration = "psu";
		MinMax_Hyd		minmaxSalt ( im, jm, km, c_0, L_hyd );
		minmaxSalt.searchMinMax_3D ( str_max_salt_concentration, str_min_salt_concentration, str_unit_salt_concentration, c, h );

//	searching of maximum and minimum values of salt balance
		string str_max_salt_balance = " max salt balance ", str_min_salt_balance = " min salt balance ", str_unit_salt_balance = "psu";
		MinMax_Hyd		minmaxSaltBalance ( im, jm, km, c_0, L_hyd );
		minmaxSaltBalance.searchMinMax_3D ( str_max_salt_balance, str_min_salt_balance, str_unit_salt_balance, Salt_Balance, h );

//	searching of maximum and minimum values of salt finger
		string str_max_salt_finger = " max salt finger ", str_min_salt_finger = " min salt finger ", str_unit_salt_finger = "psu";
		MinMax_Hyd		minmaxSaltFinger ( im, jm, km, c_0, L_hyd );
		minmaxSaltFinger.searchMinMax_3D ( str_max_salt_finger, str_min_salt_finger, str_unit_salt_finger, Salt_Finger, h );

//	searching of maximum and minimum values of salt diffusion
		string str_max_salt_diffusion = " max salt diffusion ", str_min_salt_diffusion = " min salt diffusion ", str_unit_salt_diffusion = "psu";
		MinMax_Hyd		minmaxSaltDiffusion ( im, jm, km, c_0, L_hyd );
		minmaxSaltDiffusion.searchMinMax_3D ( str_max_salt_diffusion, str_min_salt_diffusion, str_unit_salt_diffusion, Salt_Diffusion, h );

//	searching of maximum and minimum values of buoyancy force
		string str_max_BuoyancyForce_3D = " max buoyancy force ", str_min_BuoyancyForce_3D = " min buoyancy force ", str_unit_BuoyancyForce_3D = "N";
		MinMax_Hyd		minmaxBuoyancyForce_3D ( im, jm, km, c_0, L_hyd );
		minmaxBuoyancyForce_3D.searchMinMax_3D ( str_max_BuoyancyForce_3D, str_min_BuoyancyForce_3D, str_unit_BuoyancyForce_3D, BuoyancyForce_3D, h );




// 2D_fields

//	searching of maximum and minimum values of total salt volume in a column
		string str_max_salt_total = " max salt total ", str_min_salt_total = " min salt total ", str_unit_salt_total = "psu";
		MinMax_Hyd		minmaxSalt_total ( jm, km, c_0 );
		minmaxSalt_total.searchMinMax_2D ( str_max_salt_total, str_min_salt_total, str_unit_salt_total, Salt_total, h );

//	searching of maximum and minimum values of salt finger volume in a column
		string str_max_Salt_Finger = " max Salt_Finger ", str_min_Salt_Finger = " min Salt_Finger ", str_unit_Salt_Finger = "psu";
		MinMax_Hyd		minmaxSalt_finger ( jm, km, c_0 );
		minmaxSalt_finger.searchMinMax_2D ( str_max_Salt_Finger, str_min_Salt_Finger, str_unit_Salt_Finger, SaltFinger, h );

//	searching of maximum and minimum values of salt diffusion volume in a column
		string str_max_Salt_Diffusion = " max Salt_Diffusion ", str_min_Salt_Diffusion = " min Salt_Diffusion ", str_unit_Salt_Diffusion = "psu";
		MinMax_Hyd		minmaxSalt_diffusion ( jm, km, c_0 );
		minmaxSalt_diffusion.searchMinMax_2D ( str_max_Salt_Diffusion, str_min_Salt_Diffusion, str_unit_Salt_Diffusion, SaltDiffusion, h );

//	searching of maximum and minimum values of salt diffusion volume in a column
		string str_max_BuoyancyForce_2D = " max BuoyancyForce_2D ", str_min_BuoyancyForce_2D = " min BuoyancyForce_2D ", str_unit_BuoyancyForce_2D = "N";
		MinMax_Hyd		minmaxBuoyancyForce_2D ( jm, km, c_0 );
		minmaxBuoyancyForce_2D.searchMinMax_2D ( str_max_BuoyancyForce_2D, str_min_BuoyancyForce_2D, str_unit_BuoyancyForce_2D, BuoyancyForce_2D, h );

		cout << endl << " deep currents averaged for a two dimensional plane: " << endl << endl;

//	searching of maximum and minimum values of upwelling volume in a column
		string str_max_upwelling = " max upwelling ", str_min_upwelling = " min upwelling ", str_unit_upwelling = "m/s";
		MinMax_Hyd		minmaxUpwelling ( jm, km, c_0 );
		minmaxUpwelling.searchMinMax_2D ( str_max_upwelling, str_min_upwelling, str_unit_upwelling, Upwelling, h );

//	searching of maximum and minimum values of downwelling volume in a column
		string str_max_downwelling = " max downwelling ", str_min_downwelling = " min downwelling ", str_unit_downwelling = "m/s";
		MinMax_Hyd		minmaxDownwelling ( jm, km, c_0 );
		minmaxDownwelling.searchMinMax_2D ( str_max_downwelling, str_min_downwelling, str_unit_downwelling, Downwelling, h );

//	searching of maximum and minimum values of bottom water volume in a column
		string str_max_bottom_water = " max bottom water ", str_min_bottom_water = " min bottom water ", str_unit_bottom_water = "m/s";
		MinMax_Hyd		minmaxBottom_water ( jm, km, c_0 );
		minmaxBottom_water.searchMinMax_2D ( str_max_bottom_water, str_min_bottom_water, str_unit_bottom_water, BottomWater, h );


//		composition of results
		calculate_MSL.run_data ( u_0, c_0, h, u, v, w, c, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Upwelling, Downwelling, SaltFinger, SaltDiffusion, BuoyancyForce_2D, Salt_total, BottomWater );


//	restoring the velocity component and the temperature for the new time step
		oldnew.restoreOldNew_3D ( 1., u, v, w, t, p_dyn, c, un, vn, wn, tn, aux_p, cn );
	}


//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
	startPressure.computePressure_3D ( pa, rad, the, p_dyn, aux_p, h, rhs_u, rhs_v, rhs_w, aux_u, aux_v, aux_w );

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
	write_File.paraview_vtk_longal ( Name_Bathymetry_File, j_longal, pressure_iter_aux, h, p_dyn, p_stat, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

//	zonal data along constant longitudes
	k_zonal = 185;
	write_File.paraview_vtk_zonal ( Name_Bathymetry_File, k_zonal, pressure_iter_aux, h, p_dyn, p_stat, t, u, v, w, c, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

//	radial data along constant hight above ground
	i_radial = 40;
	write_File.paraview_vtk_radial ( Name_Bathymetry_File, i_radial, pressure_iter_aux, h, p_dyn, p_stat, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance, Upwelling, Downwelling, SaltFinger, SaltDiffusion, BuoyancyForce_2D, BottomWater );

//	3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
//	write_File.paraview_panorama_vts ( Name_Bathymetry_File, pressure_iter_aux, h, t, p_dyn, p_stat, u, v, w, c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

//	3-dimensional data in spherical coordinate system for a streamline pattern in a shell of a sphere
//	write_File.paraview_vts ( Name_Bathymetry_File, n, rad, the, phi, h, t, p_dyn, p_stat, u, v, w, c, rhs_u, rhs_v, rhs_w, rhs_c, rhs_p, rhs_t, aux_u, aux_v, aux_w, Salt_Finger, BuoyancyForce_3D, Salt_Balance );


//	writing of sequential data for the sequel file
	if ( SequelFile == 1 )
	{
		write_File.Hydrosphere_SequelFile_write ( Name_Bathymetry_File, n, pressure_iter, time, rad, the, phi, h, t, u, v, w, c, tn, un, vn, wn, cn, aux_u, aux_v, aux_w );
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


// reset of results to the initial value
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Downwelling.y[ j ][ k ] = 0.;
				Upwelling.y[ j ][ k ] = 0.;
				BottomWater.y[ j ][ k ] = 0.;

				SaltFinger.y[ j ][ k ] = 0.;
				SaltDiffusion.y[ j ][ k ] = 0.;
				Salt_total.y[ j ][ k ] = 0.;

				BuoyancyForce_2D.y[ j ][ k ] = 0.;

				for ( int i = 0; i < im; i++ )
				{
					h.x[ i ][ j ][ k ] = 0.;

					t.x[ i ][ j ][ k ] = 0.;
					c.x[ i ][ j ][ k ] = 0.;
					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
					tn.x[ i ][ j ][ k ] = 0.;
					cn.x[ i ][ j ][ k ] = 0.;
					un.x[ i ][ j ][ k ] = 0.;
					vn.x[ i ][ j ][ k ] = 0.;
					wn.x[ i ][ j ][ k ] = 0.;

					p_dyn.x[ i ][ j ][ k ] = 0.;
					p_stat.x[ i ][ j ][ k ] = 0.;

					rhs_t.x[ i ][ j ][ k ] = 0.;
					rhs_u.x[ i ][ j ][ k ] = 0.;
					rhs_v.x[ i ][ j ][ k ] = 0.;
					rhs_w.x[ i ][ j ][ k ] = 0.;
					rhs_c.x[ i ][ j ][ k ] = 0.;
					aux_u.x[ i ][ j ][ k ] = 0.;
					aux_v.x[ i ][ j ][ k ] = 0.;
					aux_w.x[ i ][ j ][ k ] = 0.;
					aux_p.x[ i ][ j ][ k ] = 0.;

					Salt_Finger.x[ i ][ j ][ k ] = 0.;
					Salt_Diffusion.x[ i ][ j ][ k ] = 0.;
					Salt_Balance.x[ i ][ j ][ k ] = 0.;

					BuoyancyForce_3D.x[ i ][ j ][ k ] = 0.;
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

// delete temporary files
	delete [ ] time_slice;

// 	final remarks
	cout << endl << "***** end of the Oceanic General Circulation Modell ( OGCM ) *****" << endl << endl;

	if ( velocity_iter == velocity_iter_max )	cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of maximum time steps ***** \n\n" << endl;

	if ( min <= epsres )		cout << "***** steady solution reached! *****" << endl;

	cout << endl;
	cout << "***** end of object oriented C++ program for the computation of 3D-oceanic circulation *****";
	cout << "\n\n\n\n";
}
