/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
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

#ifndef _PRESSURE_
#define _PRESSURE_

using namespace std;



class Pressure_Hyd
{
	private:
		int im, jm, km;
		double dr, dthe, dphi, dr2, dthe2, dphi2;
		double rm, rm2, sinthe, sinthe2, costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2, rm2dthe2;
		double drhs_udr, drhs_vdthe, drhs_wdphi;
		double denom, num1, num2, num3;

	public:
		Pressure_Hyd ( int, int, int, double, double, double );
		~Pressure_Hyd ();


		void computePressure_3D ( double, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void computePressure_2D ( double, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

};
#endif
