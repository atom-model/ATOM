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

#ifndef _IC_THERMOHALIN_
#define _IC_THERMOHALIN_

using namespace std;



class IC_Thermohalin
{
	private:
		int i, j, k, im, jm, km, k_half, m, i_beg, j_max, i_max, j_half, i_bottom, i_deep, i_middle, i_EIC_o, i_EIC_u, i_SCC_o, i_SCC_u, i_ECC_o, i_ECC_u;
		int i_half, j_beg, j_end, j_run, j_step, k_beg, k_end, k_run, k_step, k_exp, j_z, j_n, k_z, k_n, k_w;
		int k_a, k_b, flip, k_grad;
		int k_water, k_sequel;
		int n_smooth;
		int j1, j2, j3, jn, jd, k1, k2, k3, kn, kd;

		double dummy_1, dummy_2, dummy_3, IC_water, water_wind, Ekman_angle, vel_magnitude, alfa, beta, angle, Ekman_angle_add, Ekman, pi180;
		double d_i, d_i_half, d_i_max, d_j, d_j_half, d_j_max, t_coeff;
		double v_grad, t_Celsius, ca_max;
		double t_Cretaceous, t_Cretaceous_coeff;
		double c_Average, c_Cretaceous, c_diff;
		double rg, jmkm, v_sum, w_sum, t_sum, c_sum;

		string time_slice_comment, time_slice_number, time_slice_unit;
		string temperature_comment, temperature_gain, temperature_modern, temperature_average, temperature_unit;
		string salinity_comment, salinity_gain, salinity_modern, salinity_average, salinity_unit;

	public:

		IC_Thermohalin ( int, int, int );
		~IC_Thermohalin();


		void IC_South_Polar_Sea ( Array &, Array &, Array &, Array &, Array & );

		void IC_DeepWater ( double, Array &, Array &, Array &, Array &, Array & );

		void IC_v_w_Atmosphere ( Array &, Array &, Array &, Array & );

		void IC_v_w_Smoothing ( int, Array &, Array &, Array &, Array &, Array &, Array & );

		void IC_v_w_WestEastCoast ( Array &, Array &, Array &, Array &, Array & );

		void IC_v_w_Ekman ( Array &, Array &, Array & );

		void BC_Temperature_Salinity ( int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, Array_2D &, Array_2D &, Array_2D &, Array &, Array &, Array &, Array & );

		void IC_Atlantischer_Ozean ( Array &h, Array &u, Array &v, Array &w, Array &c );

		void IC_Indischer_Ozean ( Array &h, Array &u, Array &v, Array &w );

		void IC_Pazifischer_Ozean ( Array &h, Array &u, Array &v, Array &w );

		void IC_Nord_Polar_Meer ( Array &h, Array &u, Array &v, Array &w );

		void IC_EquatorialCurrents ( Array &h, Array &u, Array &v, Array &w );

		void BC_Surface_Temperature ( const string &, Array_2D &, Array & );

		void BC_Surface_Salinity ( const string &, Array_2D &, Array & );

		void BC_Surface_Pressure ( double, double const, double const, double const, double const, Array_2D &, Array_2D &, Array &, Array &, Array & );
};
#endif
