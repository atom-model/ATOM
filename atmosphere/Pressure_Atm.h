/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to compute the pressure independent of the other variables
*/


#include <iostream>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "BC_Thermo.h"

#ifndef _PRESSURE_
#define _PRESSURE_

using namespace std;

class Pressure_Atm
{
	private:
		int im, jm, km;

		double dr, dthe, dphi, c43, c13;

	public:
		Pressure_Atm ( int, int, int, double, double, double );
		~Pressure_Atm ();

		void computePressure_3D ( BC_Thermo &circulation, double, double, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array & );

		void computePressure_2D ( BC_Thermo &circulation, double, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array & );
};
#endif
