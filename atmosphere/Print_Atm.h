/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to print results on screen
*/


#include <iostream>

#include "Array.h"
#include "Array_1D.h"

#ifndef _PRINT_ATMOSPHERE_
#define _PRINT_ATMOSPHERE_

using namespace std;



class Print_Atmosphere
{
	private:
		int im, jm, km, nm, n;
		double time;
		
	public:
		Print_Atmosphere ( int, int, int, int, int, double );
		~Print_Atmosphere();


		void printData ( Array &, Array &, Array &, Array &, Array &, Array & );

		void printIntro ( Array_1D &, Array_1D &, Array_1D & );
};
#endif
