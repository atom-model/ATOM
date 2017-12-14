/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to restore the old by new values inside the iterational processes
*/

#include <iostream>
#include "Array.h"

#ifndef _RESTORE_
#define _RESTORE_

using namespace std;



class Restore_Atm
{
	private:
		int im, jm, km;
		// double coeff;

	public:
		Restore_Atm ( int, int, int );
		~Restore_Atm ();

		void restoreOldNew_3D ( double, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void restoreOldNew_2D ( double, Array &, Array &, Array &, Array &, Array &, Array & );
};
#endif

