/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to produce resulting data on mean sea level
*/

#include <iostream>
#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"

#ifndef _Results_Hyd_
#define _Results_Hyd_

using namespace std;

class Results_Hyd{
    private:
        int im, jm, km, h_land, h_ocean, h_point_max;
        int j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;

        double ozean_land, dr, dthe, dphi;
        double Value_1, Value_2, Value_3, Value_4, Value_5, Value_6;
        double c43, c13;
        double **aux_v, **aux_w, *aux_grad_v, *aux_grad_w;

        string name_Value_1, name_Value_2, name_Value_3, name_Value_4,
            name_Value_5, name_Value_6, name_unit_ms, name_unit_psu;
        string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading;


    public:
        Results_Hyd ( int, int, int );

        ~Results_Hyd (  );

        void run_data ( int i_beg, double dr, double dthe, double dphi, double L_hyd, double u_0,
                        double c_0, Array_1D &rad, Array_1D &the, Array &h, Array &u, Array &v,
                        Array &w, Array &c, Array &Salt_Balance, Array &Salt_Finger,
                        Array &Salt_Diffusion, Array &Buoyancy_Force_3D, Array_2D &Upwelling,
                        Array_2D &Downwelling, Array_2D &SaltFinger, Array_2D &SaltDiffusion,
                        Array_2D &BuoyancyForce_2D, Array_2D &Salt_total,
                        Array_2D &EkmanPumping, Array_1D &aux_grad_v, Array_1D &aux_grad_w, 
                        Array &aux_v, Array &aux_w );

        void land_oceanFraction ( Array & );
};
#endif
