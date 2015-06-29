/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
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

#ifndef _BC_THERMO_
#define _BC_THERMO_

using namespace std;

class BC_Thermo
{
	private:
		int i, j, k, im, jm, km, k_half, m, j_half, i_half, i_max, j_max, i_beg, i_tropopause, im_1, i_land, l;
		int j_aeq, j_pol_n, j_pol_s, j_pol_v_n, j_pol_v_s, j_fer_n, j_fer_s, j_fer_v_n, j_fer_v_s, j_had_n, j_had_s, j_had_v_n, j_had_v_s;
		int j_pac_had_n_end, j_pac_had_s_end, k_pac_w, k_pac_w_end, k_pac_e;
		int j_n, j_s;
		int *im_tropopause;
		int Ma_max, Ma_max_half;
		int j_par, j_pol, k_par, k_pol;
		int k_a, k_b, flip, k_grad, i_middle;
		int k_water, k_sequel;
		int n_smooth;

		double dummy_1, dummy_2, dummy_3, t_equator, t_tropopause, t_coeff, t_pole, c_equator, c_tropopause, c_coeff, c_pol, p_equator, p_tropopause, p_coeff, p_pol;
		double co2_equator, co2_tropopause, co2_coeff, co2_pol, co2_vegetation, co2_ocean, co2_land, co2_cretaceous;
		double ca, ua_00, ua_30, ua_60, ua_90;
		double va_Hadley_Tropopause, va_Hadley_Tropopause_15, va_Ferrel_Tropopause, va_Ferrel_Tropopause_45, va_Polar_Tropopause, va_Polar_Tropopause_75, va_Hadley_SL, va_Hadley_SL_15, va_Ferrel_SL, va_Ferrel_SL_45, va_Polar_SL, va_Polar_SL_75;
		double wa_Ferrel_Tropopause, wa_Polar_Tropopause, wa_Ferrel_SL, wa_Polar_SL;
		double wa_Hadley_SL, wa_Hadley_Tropopause;
		double wa_equator_Tropopause, wa_equator_SL, va_equator_Tropopause, va_equator_SL;
		double d_i, d_i_max, d_i_half, d_j, d_j_half, d_j_max, d_k, i_coeff_plus, i_coeff_minus, i_trop, pi180, d_j_w;
		double d_j_5n, d_j_15n, d_j_45n, d_j_75n, d_j_30n, d_j_5s, d_j_15s, d_j_45s, d_j_75s, d_j_30s, d_j_60n, d_j_60s, d_j_90n, d_j_90s, d_diff;
		double c13, c43, ep, hp, t_0, p_0, factor_velocity;
		double max_u, max_v, max_w, residuum_u, residuum_v, residuum_w;
		double t_Cretaceous, t_Cretaceous_coeff, t_360;
		double j_par_f, j_pol_f, a, b, cc, d, e, j_d, t_d, k_par_f, k_pol_f;
		double *jm_temp_asym;
		double water_wind, t_pol;
		double rR, rg, jmkm, u_sum, v_sum, w_sum, t_sum, c_sum;

		char Temperature_West[50], Salinity_West[50], Temperature_East[50], Salinity_East[50];
 
		string time_slice_comment, time_slice_number, time_slice_unit;
		string temperature_comment, temperature_gain, temperature_modern, temperature_average, temperature_unit;
		string co2_comment, co2_gain, co2_modern, co2_average, co2_unit;

	public:
		BC_Thermo ( int, int, int, Array &, Array &, Array &, Array & );
		~BC_Thermo();


		void IC_CellStructure ( Array &, Array &, Array & );

		void IC_v_w_WestEastCoast ( Array &, Array &, Array &, Array & );

		void BC_Temperature ( int, int, int, int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, Array_2D &, Array_2D &, Array &, Array &, Array &, Array & );

		void BC_WaterVapour ( int, double, double, double, double, double, double, double, Array_2D &, Array &, Array &, Array &, Array &, Array &, Array_2D & );

		void BC_CO2 ( int, double, double, double, double, double, double, double, double, Array_2D &, Array_2D &, Array &, Array &, Array &, Array & );

		void BC_Tropopause ( int, double, double, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_Surface_Temperature ( const string &, Array_2D &, Array & );

		void BC_Surface_Precipitation ( const string &, Array_2D & );

		void BC_Pressure ( double const, double const, double const, double const, Array_2D &, Array_2D &, Array &, Array &, Array & );
};
#endif
