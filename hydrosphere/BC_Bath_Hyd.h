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
		int im, jm, km;

	public:

		BC_Bathymetry_Hydrosphere ( int im, int jm, int km );
		~BC_Bathymetry_Hydrosphere();

        void BC_SeaGround(const string &bathymetry_file, double L_hyd, Array &h, Array_2D &Bathymetry);

        void BC_SolidGround ( double ca, double ta, double pa, Array &h, Array &t, Array &u, Array &v, 
            Array &w, Array &p_dyn, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &p_dynn, Array &cn);
};
#endif
