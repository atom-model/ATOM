/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 3D arrays
*/


#include <iostream>

#ifndef _ARRAY_
#define _ARRAY_

using namespace std;


class Array
{
	private:
		int im, jm, km;
		// double aa;

	public:
		double ***x;

		Array ( );

		~Array ( );

		void printArray ( int, int, int );
		void initArray ( int, int, int, double );
};
#endif
