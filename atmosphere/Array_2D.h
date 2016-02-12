/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 2D arrays
*/


#include <iostream>

#ifndef _ARRAY_2D_
#define _ARRAY_2D_

using namespace std;


class Array_2D
{
	private:
		int jm, km;
		double bb;

	public:
		double **y;

		Array_2D ( int, int, double );
		~Array_2D();

		void printArray_2D();
};
#endif
