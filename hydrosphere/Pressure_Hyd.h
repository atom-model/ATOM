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

#ifndef _PRESSURE_
#define _PRESSURE_

using namespace std;

class Pressure_Hyd
{
    private:
        int im, jm, km;

        double dr, dthe, dphi, c43, c13;

    public:
        Pressure_Hyd ( int, int, int, double, double, double );
        ~Pressure_Hyd ();

        void computePressure_3D ( double u_0, double r_0_water, Array_1D &rad, Array_1D &the,
                 Array &p_dyn, Array &p_dynn, Array &h, Array &aux_u, Array &aux_v, Array &aux_w );

        void computePressure_2D ( double u_0, double r_0_water, Array_1D &rad, Array_1D &the,
                 Array &p_dyn, Array &p_dynn, Array &h, Array &aux_v, Array &aux_w );
};
#endif
