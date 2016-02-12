/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the salt concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to read and prepare the bathymetric and topografic data
*/


#include <iostream>

#include "Array.h"
#include "Array_2D.h"

#ifndef _BC_BATHYMETRY_HYDROSPHERE_
#define _BC_BATHYMETRY_HYDROSPHERE_

using namespace std;



class BC_Bathymetry_Hydrosphere
{
	private:
		int i, j, k, im, jm, km, i_boden, l;
		int k_water, k_sequel, flip;

		double dummy_1, dummy_2, dummy_3, h_max;

	public:

//		BC_Bathymetry_Hydrosphere ( const string &, int, int, int, Array &, Array & );
		BC_Bathymetry_Hydrosphere ( int, int, int );
		~BC_Bathymetry_Hydrosphere();

		void BC_SeaGround ( const string &, Array &, Array & );

		void BC_SolidGround ( double, double, double, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D & );

		void IC_SeaGroundGaps ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );
};
#endif
