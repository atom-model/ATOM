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
#include "cHydrosphereModel.h"

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
#include "BC_Hyd.h"
#include "BC_Bath_Hyd.h"
#include "BC_Thermohalin.h"
#include "Accuracy_Hyd.h"
#include "RHS_Hyd.h"
#include "RungeKutta_Hyd.h"
#include "PostProcess_Hyd.h"
#include "Pressure_Hyd.h"
#include "Restore_Hyd.h"
#include "MinMax_Hyd.h"
#include "Results_Hyd.h"

#include "Config.h"
#include "tinyxml2.h"
#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;

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

// set default configuration
	SetDefaultConfig();
}




cHydrosphereModel::~cHydrosphereModel() { }

#include "cHydrosphereDefaults.cpp.inc"



void cHydrosphereModel::LoadConfig(const char *filename) {
	XMLDocument doc;
	XMLError err = doc.LoadFile(filename);


	if (err) {
		doc.PrintError();
		throw std::invalid_argument("   couldn't load config file inside cHydrosphereModel");
	}


	XMLElement *atom = doc.FirstChildElement("atom");
	if (!atom) {
		return;
	}


	XMLElement* elem_common = doc.FirstChildElement( "atom" )->FirstChildElement( "common" );
	if (!elem_common) {
		return;
	}


	XMLElement* elem_hydrosphere = doc.FirstChildElement( "atom" )->FirstChildElement( "hydrosphere" );
	if (!elem_hydrosphere) {
		return;
	}


#include "HydrosphereLoadConfig.cpp.inc"
}




void cHydrosphereModel::RunTimeSlice(int Ma)
{
// maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km ), 
// maximum number of overall iterations ( n ),
// maximum number of inner velocity loop iterations ( velocity_iter_max ),
// maximum number of outer pressure loop iterations ( pressure_iter_max )


    mkdir(output_path.c_str(), 0777);

	const int im = 41, jm = 181, km = 361, nm = 200;

	int j_res = 0.0, k_res = 0.0;

	int Ma_max = 300;																	// parabolic temperature distribution 300 Ma back
	int Ma_max_half = 150;															// half of time scale

	constexpr double pi180 = 180./ M_PI;								// pi180 = 57.3
	constexpr double the_degree = 1.;										// compares to 1° step size laterally
	constexpr double phi_degree = 1.;										// compares to 1° step size longitudinally

	double dthe = the_degree / pi180;										// dthe = the_degree / pi180 = 1.0 / 57.3 = 0.01745, 180 * .01745 = 3.141
	double dphi = phi_degree / pi180;										// dphi = phi_degree / pi180 = 1.0 / 57.3 = 0.01745, 360 * .01745 = 6.282
	double dr = 0.025;																	// 0.025 x 40 = 1.0 compares to 16 km : 40 = 150 m for 1 radial step
	double dt = 0.0001;																// time step satisfies the CFL condition

	const double the0 = 0.;														// North Pole
	const double phi0 = 0.;														// zero meridian in Greenwich
	const double r0 = 1.;															// earth's radius is r_earth = 6731 km, here it is assumed to be infinity, circumference of the earth 40074 km

// 1D arrays
	Array_1D rad(im, 1.); // radial coordinate direction
	Array_1D the(jm, 2.); // lateral coordinate direction
	Array_1D phi(km, 3.); // longitudinal coordinate direction

// 2D arrays
	Array_2D Upwelling(jm, km, 0.); // upwelling
	Array_2D Downwelling(jm, km, 0.); // downwelling
	Array_2D BottomWater(jm, km, 0.); // 2D bottom water summed up in a vertical column

	Array_2D SaltFinger(jm, km, 0.);											// salt bulge of higher density
	Array_2D SaltDiffusion(jm, km, 0.);										// salt bulge of lower density
	Array_2D Salt_total(jm, km, 0.);											// rate of salt summed up in a vertical column

	Array_2D BuoyancyForce_2D(jm, km, 0.); // radiation balance at the surface

// 3D arrays
	Array h(im, jm, km, 0.); // bathymetry, depth from sea level

	Array t(im, jm, km, ta); // temperature
	Array u(im, jm, km, ua); // u-component velocity component in r-direction
	Array v(im, jm, km, va); // v-component velocity component in theta-direction
	Array w(im, jm, km, wa); // w-component velocity component in phi-direction
	Array c(im, jm, km, ca); // water vapour

	Array tn(im, jm, km, ta); // temperature new
	Array un(im, jm, km, ua); // u-velocity component in r-direction new
	Array vn(im, jm, km, va); // v-velocity component in theta-direction new
	Array wn(im, jm, km, wa); // w-velocity component in phi-direction new
	Array cn(im, jm, km, ca); // water vapour new

	Array p_dyn(im, jm, km, pa); // dynamic pressure
	Array p_dynn(im, jm, km, pa); // dynamic pressure new
	Array p_stat(im, jm, km, pa); // static pressure

	Array rhs_t(im, jm, km, 0.); // auxilliar field RHS temperature
	Array rhs_u(im, jm, km, 0.); // auxilliar field RHS u-velocity component
	Array rhs_v(im, jm, km, 0.); // auxilliar field RHS v-velocity component
	Array rhs_w(im, jm, km, 0.); // auxilliar field RHS w-velocity component
	Array rhs_c(im, jm, km, 0.); // auxilliar field RHS water vapour
	Array rhs_p(im, jm, km, 0.); // auxilliar field RHS water vapour

	Array aux_u(im, jm, km, 0.); // auxilliar field u-velocity component
	Array aux_v(im, jm, km, 0.); // auxilliar field v-velocity component
	Array aux_w(im, jm, km, 0.); // auxilliar field w-velocity component

	Array Salt_Finger(im, jm, km, 0.); // salt bulge of higher density
	Array Salt_Diffusion(im, jm, km, 0.); // salt bulge of lowerer density
	Array Salt_Balance(im, jm, km, 0.); // +/- salt balance

	Array BuoyancyForce_3D(im, jm, km, 0.); // 3D buoyancy force

	cout.precision ( 6 );
	cout.setf ( ios::fixed );

//	Coordinate system in form of a spherical shell
//	rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
	rad.Coordinates ( im, r0, dr );
	the.Coordinates ( jm, the0, dthe );
	phi.Coordinates ( km, phi0, dphi );


//	cout << endl << "§§§§§§§§§§§§§§§§§§ after BC_Surface_Temperature_NASA §§§§§§§§§§§§§§§§§§§§§§§§" << endl;
//	cout << endl << " ***** printout of 3D-field temperature ***** " << endl << endl;
//	t.printArray( im, jm, km );

//	cout << endl << " ***** printout of 2D-field vegetation ***** " << endl << endl;
//	Vegetation.printArray_2D( jm, km );

//	cout << endl << " ***** printout of 1D-field radius ***** " << endl << endl;
//	rad.printArray_1D( im );

//	initial values for the number of computed steps and the time
	int n = 1;
	int velocity_iter = 0;
	int velocity_iter_2D = 0;
	int pressure_iter = 0;
	int pressure_iter_2D = 0;
	int switch_2D = 0;
	double residuum;
	double residuum_old = 0.;

// radial expansion of the computational field for the computation of initial values
	int i_max = im - 1;																	// corresponds to sea level
	int i_beg = 33;																		// compares to an ocean depth of 1050 m, location of the thermocline

// naming a transfer file for the v- and w-vlocity component of the atmosphere at sea level
	string Name_v_w_Transfer_File;
	stringstream ssName_v_w_Transfer_File;

// naming a file to read the surface temperature of the modern world
	string Name_SurfaceTemperature_File; 
	stringstream ssNameSurfaceTemperature;
	ssNameSurfaceTemperature << "../data/" << "SurfaceTemperature_NASA.xyz";
	Name_SurfaceTemperature_File = ssNameSurfaceTemperature.str();

// naming a file to read the surface salinity of the modern world
	string Name_SurfaceSalinity_File; 
	stringstream ssNameSurfaceSalinity;
	ssNameSurfaceSalinity << "../data/" << "SurfaceSalinity_NASA.xyz";
	Name_SurfaceSalinity_File = ssNameSurfaceSalinity.str();


    string bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
   string bathymetry_filepath = bathymetry_path + "/" + bathymetry_name;
    string input_path = output_path;

    cout << "\n   Input is being read from " << input_path << "\n";
    cout << "\n   Output is being written to " << output_path << "\n";
    cout << "   Ma = " << Ma << "\n";
    cout << "   bathymetry_path = " << bathymetry_path << "\n";
    cout << "   bathymetry_filepath = " << bathymetry_filepath << "\n\n";


	if (verbose) {
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
	}



//	class PostProcess for data transport, read and write
	PostProcess_Hydrosphere		read_Transfer ( im, jm, km, input_path, output_path );
	read_Transfer.Atmosphere_TransferFile_read ( bathymetry_name, v, w, t, p_dyn );


	cout << "***** time slice for the Oceanic Global Circulation Modell ( OGCM ) is:    Ma = " << Ma << " million years" << endl << endl;
	cout << "***** bathymetry/topography given by the x-y-z data set:    " << bathymetry_name.c_str() << endl << endl;


// class BC_Bathymetry_Hydrosphere for the geometrical boundary condition of the computational area
	BC_Bathymetry_Hydrosphere		depth ( im, jm, km );

// 	class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
	depth.BC_SeaGround(bathymetry_path, bathymetry_name, L_hyd, h, aux_w);

// class BC_Hydrosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
	BC_Hydrosphere		boundary ( im, jm, km );

// class RHS_Hydrosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
	RHS_Hydrosphere		prepare ( im, jm, km, r0, dt, dr, dthe, dphi, re, ec, sc, g, pr, Buoyancy );
	RHS_Hydrosphere		prepare_2D ( jm, km, dthe, dphi, re );

// class RungeKutta_Hydrosphere for the explicit solution of the Navier-Stokes equations
	RungeKutta_Hydrosphere		result ( n, im, jm, km, dt );

// class Pressure for the subsequent computation of the pressure by a separat Euler equation
	Pressure_Hyd		startPressure ( im, jm, km, dr, dthe, dphi );

// class Restore to restore between iterational results
	Restore_Hyd		oldnew( im, jm, km );

// class Results_MSL_Hyd to compute and show results on the mean sea level, MSL
	Results_Hyd		calculate_MSL ( im, jm, km );

// class File_NetCDF to write results in the format of a netCDF-file
// File_NetCDF_Hyd		printoutNetCDF ( im, jm, km );

// class BC_Thermohalin for the initial and boundary conditions of the flow properties
	BC_Thermohalin		oceanflow ( im, jm, km, i_beg, i_max, Ma, Ma_max, Ma_max_half, dr, g, r_0_water, ua, va, wa, ta, ca, pa, u_0, p_0, t_0, c_0, cp_w, L_hyd, t_average, t_cretaceous_max, t_equator, t_pole, input_path );

//	surface temperature from World Ocean Atlas 2009 given as boundary condition
	if ( Ma == 0 ) oceanflow.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, t );

//	surface salinity from World Ocean Atlas 2009 given as boundary condition
//	if ( Ma == 0 ) oceanflow.BC_Surface_Salinity_NASA ( Name_SurfaceSalinity_File, c );

//	import of surface v- and w-velocity components from AGCM, surface velocity reduced to 3% of the wind velocity
	oceanflow.IC_v_w_Atmosphere ( h, u, v, w );

//	salinity distribution as initial condition in 3 dimensions
	oceanflow.BC_Temperature_Salinity ( h, t, c, p_dyn );

//  surface pressure computed by surface temperature with gas equation
	oceanflow.BC_Pressure ( p_stat, t, h );

//	initial conditions for v and w velocity components at the sea surface close to east or west coasts, to close gyres
	oceanflow.IC_v_w_WestEastCoast ( h, u, v, w, c );

//	storing of velocity components, pressure and temperature for iteration start
	oldnew.restoreOldNew_3D(.9, u, v, w, t, p_dyn, c, un, vn, wn, tn, p_dynn, cn);
	oldnew.restoreOldNew_2D(.9, v, w, p_dyn, p_dynn, vn, wn);

// computation of the ratio ocean to land areas
	calculate_MSL.land_oceanFraction ( h );


// ******************************************   start of pressure and velocity iterations ************************************************************************

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
				cout << " 2D OGCM iterational process" << endl;
				cout << " max total iteration number nm = " << nm << endl << endl;

				cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations " << endl << endl << " Ma = " << Ma << "     n = " << n << "    velocity_iter_max_2D = " << velocity_iter_max_2D << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " << pressure_iter_max_2D << "    pressure_iter_2D = " << pressure_iter_2D << endl;

//		class BC_Atmosphaere for the geometry of a shell of a sphere
				boundary.RB_theta ( ca, ta, pa, t, u, v, w, p_dyn, c );
				boundary.RB_phi ( t, u, v, w, p_dyn, c );

//		old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( emin )
				Accuracy_Hyd		min_Residuum_old_2D ( im, jm, km, dthe, dphi );
				min_Residuum_old_2D.residuumQuery_2D ( rad, the, v, w );
				emin = min_Residuum_old_2D.out_min (  );

				residuum_old = emin;

//		class RungeKutta for the solution of the differential equations describing the flow properties
				result.solveRungeKutta_2D_Hydrosphere ( prepare_2D, n, rad, the, phi, rhs_v, rhs_w, rhs_p, h, v, w, p_dyn, vn, wn, p_dynn, aux_v, aux_w );

//		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( emin )
				Accuracy_Hyd		min_Residuum_2D ( im, jm, km, dthe, dphi );
				min_Residuum_2D.residuumQuery_2D ( rad, the, v, w );
				emin = min_Residuum_2D.out_min (  );
				j_res = min_Residuum_2D.out_j_res (  );
				k_res = min_Residuum_2D.out_k_res (  );

				residuum = emin;
				emin = fabs ( ( residuum - residuum_old ) / residuum_old );

//		state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
				Accuracy_Hyd		min_Stationary_2D ( n, nm, Ma, im, jm, km, emin, j_res, k_res, velocity_iter_2D, pressure_iter_2D, velocity_iter_max_2D, pressure_iter_max_2D );
				min_Stationary_2D.steadyQuery_2D ( h, v, vn, w, wn, p_dyn, p_dynn );

				oldnew.restoreOldNew_2D(1., v, w, p_dyn, p_dynn, vn, wn);

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
			oldnew.restoreOldNew_3D(1., u, v, w, t, p_dyn, c, un, vn, wn, tn, p_dynn, cn);

			cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations " << endl << endl << " Ma = " << Ma << "     n = " << n << "    velocity_iter_max = " << velocity_iter_max << "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << "    pressure_iter = " << pressure_iter << endl;

//		old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( emin )
			Accuracy_Hyd		min_Residuum_old ( im, jm, km, dr, dthe, dphi );
			min_Residuum_old.residuumQuery_3D ( rad, the, u, v, w );
			emin = min_Residuum_old.out_min (  );

			residuum_old = emin;

//		class RB_Hydrosphaere for the geometry of a shell of a sphere
			boundary.RB_radius ( ca, ta, pa, dr, rad, t, u, v, w, p_dyn, c );
			boundary.RB_theta ( ca, ta, pa, t, u, v, w, p_dyn, c );
			boundary.RB_phi ( t, u, v, w, p_dyn, c );

//		class RungeKutta for the solution of the differential equations describing the flow properties
			result.solveRungeKutta_3D_Hydrosphere ( prepare, n, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, rhs_t, rhs_u, rhs_v, rhs_w, rhs_p, rhs_c, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance, p_stat );

//		class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
			depth.BC_SolidGround ( ca, ta, pa, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn );

//		new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( emin )
			Accuracy_Hyd		min_Residuum ( im, jm, km, dr, dthe, dphi );
			min_Residuum.residuumQuery_3D ( rad, the, u, v, w );
			emin = min_Residuum.out_min (  );
			int i_res = min_Residuum.out_i_res (  );
			j_res = min_Residuum.out_j_res (  );
			k_res = min_Residuum.out_k_res (  );

			residuum = emin;
			emin = fabs ( ( residuum - residuum_old ) / residuum_old );

//		statements on the convergence und iterational process
			Accuracy_Hyd		min_Stationary ( n, nm, Ma, im, jm, km, emin, i_res, j_res, k_res, velocity_iter, pressure_iter, velocity_iter_max, pressure_iter_max, L_hyd );
			min_Stationary.steadyQuery_3D ( u, un, v, vn, w, wn, t, tn, c, cn, p_dyn, p_dynn );

// 3D_fields

//		searching of maximum and minimum values of temperature
			string str_max_temperature = " max temperature ", str_min_temperature = " min temperature ", str_unit_temperature = "C";
			MinMax_Hyd		minmaxTemperature ( im, jm, km, c_0, L_hyd );
			minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

//	searching of maximum and minimum values of u-component
			string str_max_u = " max 3D u-component ", str_min_u = " min 3D u-component ", str_unit_u = "m/s";
			MinMax_Hyd		minmax_u ( im, jm, km, c_0, L_hyd );
			minmax_u.searchMinMax_3D ( str_max_u, str_min_u, str_unit_u, u, h );

//	searching of maximum and minimum values of v-component
			string str_max_v = " max 3D v-component ", str_min_v = " min 3D v-component ", str_unit_v = "m/s";
			MinMax_Hyd		minmax_v ( im, jm, km, c_0, L_hyd );
			minmax_v.searchMinMax_3D ( str_max_v, str_min_v, str_unit_v, v, h );

//	searching of maximum and minimum values of w-component
			string str_max_w = " max 3D w-component ", str_min_w = " min 3D w-component ", str_unit_w = "m/s";
			MinMax_Hyd		minmax_w ( im, jm, km, c_0, L_hyd );
			minmax_w.searchMinMax_3D ( str_max_w, str_min_w, str_unit_w, w, h );

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
			string str_max_BuoyancyForce_2D = " max BuoyancyForce_2D ", str_min_BuoyancyForce_2D = " min BuoyancyForce_2D ", str_unit_BuoyancyForce_2D = "N/m2";
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

//  restoring the velocity component and the temperature for the new time step
			oldnew.restoreOldNew_3D(1., u, v, w, t, p_dyn, c, un, vn, wn, tn, p_dynn, cn);
//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

			n++;
		}
//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of velocity loop_3D: if ( velocity_iter > velocity_iter_max )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//	pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
		if ( pressure_iter == 1 ) 			startPressure.computePressure_3D ( pa, rad, the, p_dyn, p_dynn, h, rhs_u, rhs_v, rhs_w, aux_u, aux_v, aux_w );


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


//	printout in ParaView and plot files

//	class PostProcess_Hydrosphaere for the printing of results
	PostProcess_Hydrosphere		write_File ( im, jm, km, input_path, output_path );

	int j_longal = 75;
	write_File.paraview_vtk_longal ( bathymetry_name, j_longal, n, u_0, r_0_water, h, p_dyn, p_stat, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

//	zonal data along constant longitudes
	int k_zonal = 185;
	write_File.paraview_vtk_zonal ( bathymetry_name, k_zonal, n, u_0, r_0_water, h, p_dyn, p_stat, t, u, v, w, c, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

//	radial data along constant hight above ground
	int i_radial = 40;
	write_File.paraview_vtk_radial ( bathymetry_name, i_radial, n, u_0, r_0_water, h, p_dyn, p_stat, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance, Upwelling, Downwelling, SaltFinger, SaltDiffusion, BuoyancyForce_2D, BottomWater );

//	3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
	write_File.paraview_panorama_vts ( bathymetry_name, n, u_0, r_0_water, h, t, p_dyn, p_stat, u, v, w, c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

//	writing of plot data in the PlotData file
	PostProcess_Hydrosphere		ppa ( im, jm, km, input_path, output_path );
	ppa.Hydrosphere_PlotData ( bathymetry_name, u_0, h, v, w, t, c, BottomWater, Upwelling, Downwelling );


// reset of results to the initial value
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Upwelling.y[ j ][ k ] = 0.;
			Downwelling.y[ j ][ k ] = 0.;
			BottomWater.y[ j ][ k ] = 0.;

			BuoyancyForce_2D.y[ j ][ k ] = 0.;

			SaltFinger.y[ j ][ k ] = 0.;
			SaltDiffusion.y[ j ][ k ] = 0.;
			Salt_total.y[ j ][ k ] = 0.;

			for ( int i = 0; i < im; i++ )
			{
				h.x[ i ][ j ][ k ] = 0.;

				c.x[ i ][ j ][ k ] = 0.;
				u.x[ i ][ j ][ k ] = 0.;
				v.x[ i ][ j ][ k ] = 0.;
				w.x[ i ][ j ][ k ] = 0.;

				un.x[ i ][ j ][ k ] = 0.;
				vn.x[ i ][ j ][ k ] = 0.;
				wn.x[ i ][ j ][ k ] = 0.;
				cn.x[ i ][ j ][ k ] = 0.;

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

				Salt_Finger.x[ i ][ j ][ k ] = 0.;
				Salt_Diffusion.x[ i ][ j ][ k ] = 0.;
				Salt_Balance.x[ i ][ j ][ k ] = 0.;

				BuoyancyForce_3D.x[ i ][ j ][ k ] = 0.;
			}
		}
	}

	//   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	//  final remarks
    cout << endl << "***** end of the Hydrosphere General Circulation Modell ( OGCM ) *****" << endl << endl;

    if ( velocity_iter == velocity_iter_max )   cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of maximum time steps ***** \n\n" << endl;

    if ( emin <= epsres )        cout << "***** steady solution reached! *****" << endl;
}





void cHydrosphereModel::Run() {
    // create the output dir
//	output_path = output_path + "-" + std::to_string ( Ma );
    mkdir(output_path.c_str(), 0777);

    cout << "Output is being written to " << output_path << "\n";

    // write out the config for reproducibility
    // disabled for now
    // std::stringstream output_config_path;
    // output_config_path << output_path << "/config_hyd.xml";
    // WriteConfig(output_config_path.str().c_str());


	if (verbose) {
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
	}




	// time slices to be run after actualizing 
	int i_time_slice_max = 15;
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

	for (int i_time_slice = 0; i_time_slice < i_time_slice_max; i_time_slice++) {
		int Ma = time_slice[i_time_slice];
		RunTimeSlice(Ma);
	}
//   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//  final remarks
    cout << endl << "***** end of the Hydrosphere General Circulation Modell ( OGCM ) *****" << endl << endl;

    cout << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-hydrospheric circulation *****";
    cout << "\n\n\n\n";

	delete [ ] time_slice;
}
