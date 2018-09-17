/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to prepare the coordinate system
*/

#include <iostream>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"

#ifndef _BC_HYDROSPHERE_
#define _BC_HYDROSPHERE_

using namespace std;

class BC_Hydrosphere{
private:
    int im, jm, km;
    double c13, c43;

public:
    BC_Hydrosphere ( int, int, int );
    ~BC_Hydrosphere();


    void RB_radius ( double, double, double, double, Array_1D &, Array &, Array &,
        Array &, Array &, Array &, Array & );

    void RB_theta ( double, double, double, Array &, Array &, Array &, Array &,
        Array &, Array & );

    void RB_phi ( Array &, Array &, Array &, Array &, Array &, Array & );
};
#endif
