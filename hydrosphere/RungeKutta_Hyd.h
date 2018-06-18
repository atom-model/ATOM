/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include <iostream>
#include "Array.h"
#include "Array_1D.h"
#include "RHS_Hyd.h"

#ifndef _RUNGEKUTTA_HYDROSPHERE_
#define _RUNGEKUTTA_HYDROSPHERE_

using namespace std;


class RungeKutta_Hydrosphere
{
	private:
		int n, im, jm, km;

		double dt, kt1, ku1, kv1, kw1, kc1, kp1, kt2, ku2, kv2, kw2, kc2, kp2, kt3, ku3, kv3, kw3, kc3, kp3, kt4, ku4, kv4, kw4, kc4, kp4;

	public:
		RungeKutta_Hydrosphere ( int, int, int, int, double );
		~RungeKutta_Hydrosphere ();

		void solveRungeKutta_3D_Hydrosphere ( RHS_Hydrosphere &, int &, double, double, double, double, double, double, double, double, double, double, Array_1D &, Array_1D &, Array_1D &, Array_2D &, Array_2D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D & );

		void solveRungeKutta_2D_Hydrosphere ( RHS_Hydrosphere &, int &, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

};
#endif
