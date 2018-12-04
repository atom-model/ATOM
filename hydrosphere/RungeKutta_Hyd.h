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

class RungeKutta_Hydrosphere{
     private:
        int im, jm, km;

        double dt, kt1, ku1, kv1, kw1, kc1, kt2, ku2, kv2, kw2, kc2,
                     kt3, ku3, kv3, kw3, kc3, kt4, ku4, kv4, kw4, kc4;

    public:
        RungeKutta_Hydrosphere ( int, int, int, double );
         ~RungeKutta_Hydrosphere ();

        void solveRungeKutta_3D_Hydrosphere ( RHS_Hydrosphere &prepare,
                   int &iter_cnt, double L_hyd, double g, double cp_w, double u_0, double t_0, double c_0,
                   double r_0_water, double ta, double pa, double ca, Array_1D &rad, Array_1D &the, Array_1D &phi,
                   Array_2D &Evaporation_Dalton, Array_2D &Precipitation, Array &h, Array &rhs_t,
                   Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c,
                   Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c,
                   Array &tn, Array &un, Array &vn, Array &wn, Array &p_dynn, Array &cn,
                   Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion,
                   Array &Buoyancy_Force, Array &Salt_Balance, Array &p_stat,
                   Array &r_water, Array &r_salt_water, Array_2D &Bathymetry );

        void solveRungeKutta_2D_Hydrosphere ( RHS_Hydrosphere &, int &,
                    double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &,
                    Array &, Array &, Array &, Array &, Array &, Array &, Array & );
};
#endif
