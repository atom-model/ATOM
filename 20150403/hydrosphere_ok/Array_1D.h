/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 1D arrays
*/


#include <iostream>

#ifndef _ARRAY_1D_
#define _ARRAY_1D_

using namespace std;


class Array_1D
{
	private:
		int im, jm, km;
		double cc, z0, dz;

	public:
		double *z;

		Array_1D ( int, double, double, double );
		~Array_1D();

		void printArray_1D();

		void Coordinates ();
};
#endif
