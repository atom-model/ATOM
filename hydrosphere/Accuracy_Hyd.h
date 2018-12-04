/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to surveil the accuracy of the iterations
*/

#include <iostream>

#include "Array.h"
#include "Array_1D.h"

#ifndef _ACCURACY_
#define _ACCURACY_

using namespace std;

class Accuracy_Hyd{
    private:
        int n, nm, im, jm, km, velocity_iter_2D, pressure_iter_2D, velocity_iter, pressure_iter,
            velocity_iter_max, pressure_iter_max, velocity_iter_max_2D, pressure_iter_max_2D, Ma;
        int i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_p, j_p, k_p;
        int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;
        int i_res, j_res, k_res;

        double dr, dthe, dphi;
        double sinthe, costhe, rmsinthe;
        double dudr, dvdthe, dwdphi;
        double residuum, max_u, max_v, max_w, max_t, max_c, max_p;
        double Value, L_hyd;
        double min, min_u, min_v, min_w, min_t, min_c, min_p;

        string name_Value;
        string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading;

    public:
        Accuracy_Hyd ( int, int, int, double, double );
        Accuracy_Hyd ( int, int, int, double, double, double );
        Accuracy_Hyd ( int, int, int, int, int, int, double, int, int, int, int, int, int );
        Accuracy_Hyd ( int, int, int, int, int, int, double, int, int, int, int, int, int, int, double );

        ~Accuracy_Hyd ();

        double residuumQuery_2D ( Array_1D &, Array_1D &, Array &, Array & );
        double residuumQuery_3D ( double L_hyd, Array_1D &, Array_1D &, Array &, Array &, Array & );

        double steadyQuery_2D ( Array &, Array &, Array &, Array &, Array &, Array &, Array & );
        double steadyQuery_3D ( double L_hyd, Array_1D &rad, Array &u, Array &un, Array &v, Array &vn,
                                Array &w, Array &wn, Array &t, Array &tn, Array &c, Array &cn,
                                Array &p_dyn, Array &p_dynn );

        double out_min () const;
        int out_i_res () const;
        int out_j_res () const;
        int out_k_res () const;
};
#endif
