/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to prepare the boundary and initial conditions for diverse variables
*/

#include <iostream>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"

#ifndef _BC_THERMOHALIN_
#define _BC_THERMOHALIN_

using namespace std;

class BC_Thermohalin{
    private:
        int j, k, im, jm, km, m, i_beg, j_max, i_max, j_half, i_bottom, i_deep, i_middle,
            i_EIC_o, i_EIC_u, i_SCC_o, i_SCC_u, i_ECC_o, i_ECC_u, i_EUC_o, i_EUC_u;
        int i_half, j_beg, j_end, j_run, j_step, k_beg, k_end, k_run, k_step, k_exp, j_z,
            j_n, k_z, k_n, k_w;
        int k_a, k_b, flip, k_grad;
        int k_water, k_sequel;
        int Ma, Ma_max, Ma_max_half;
        int j1, j2, j3, jn, jd, k1, k2, k3, kn, kd;

        double dummy_1, dummy_2, dummy_3, IC_water, water_wind, c_0,
            Ekman_angle, vel_magnitude, alfa, beta, angle, Ekman_angle_add, Ekman, pi180;
        double t_equator, t_pole, d_i, d_i_half, d_i_middle, d_i_beg, d_i_max,
            d_j, d_j_half, d_j_max, t_coeff, c_average, c_cretaceous, p_0, t_0;
        double v_grad, t_Celsius, u_max;
        double t_cretaceous, t_cretaceous_max, t_cretaceous_coeff;
        double rg;
        double dr, g, r_0_water, u_0, cp_w, L_hyd, t_average;

        string time_slice_comment, time_slice_number, time_slice_unit;
        string temperature_comment, temperature_gain, temperature_modern,
            temperature_average, temperature_unit, temperature_cretaceous,
            temperature_average_cret;
        string salinity_comment, salinity_gain, salinity_modern, salinity_average,
        salinity_unit, salinity_cretaceous, salinity_average_cret;

        string input_path;

    public:

        BC_Thermohalin (int, int, int, int , int , int, int, int, double, double, double,
            double, double, double, double, double, double,  double, double, double,
            double, const string &);

        ~BC_Thermohalin();

        void IC_EkmanSpiral ( Array_1D &, Array_1D &, Array &, Array &, Array & );

//        void IC_u_WestEastCoast ( Array_1D &rad, Array &h, Array &u, Array &v, Array &w, Array &un, Array &vn, Array &wn );

        void BC_Temperature_Salinity ( Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

        void BC_Surface_Temperature_NASA ( const string &, Array & );

        void BC_Surface_Salinity_NASA ( const string &, Array & );

        void Value_Limitation_Hyd ( Array &, Array &, Array &, Array &, Array &, Array &, Array & );

        void Pressure_Limitation_Hyd ( Array &, Array & );
};
#endif
