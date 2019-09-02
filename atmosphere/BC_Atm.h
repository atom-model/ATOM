/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the coordinate system
*/

#include <iostream>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"

#ifndef _BC_ATMOSPHAERE_
#define _BC_ATMOSPHERE_

using namespace std;

class BC_Atmosphere{
    private:
        int im, jm, km;
        double c13, c43;
        double t_tropopause;
    public:
        BC_Atmosphere(int, int, int, double);
        ~BC_Atmosphere();
        void BC_radius(Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &);
        void BC_theta(Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &);
        void BC_phi(Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &);
};
#endif
