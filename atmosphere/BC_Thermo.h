/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
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
		int i, j, k, im, jm, km, k_half, m, j_half, i_half, i_max, j_max, k_max, i_beg, im_1, i_land, l;
		int j_aeq, j_pol_n, j_pol_s, j_pol_v_n, j_pol_v_s, j_fer_n, j_fer_s, j_fer_v_n, j_fer_v_s, j_had_n, j_had_s, j_had_v_n, j_had_v_s;
		int j_pac_had_n_end, j_pac_had_s_end, k_pac_w, k_pac_w_end, k_pac_e;
		int j_n, j_s;
		int *im_tropopause;
		int sun, Ma, Ma_max, Ma_max_half;
		int j_par, j_pol, k_par, k_pol;
		int k_a, k_b, flip, k_grad, k_grad_init, i_middle;
		int k_water, k_sequel;
		int n_smooth;
		int j_r, k_r, j_sun;
		int RadiationModel, sun_position_lat, sun_position_lon, declination;

		double d_k_half, d_k_max, t_eff_earth, r, rad_bal_minus; 
		double dummy_1, dummy_2, dummy_3;
		double ca, ua_00, ua_30, ua_60, ua_90, u_coeff, u_end;
		double va_Hadley_Tropopause, va_Hadley_Tropopause_15, va_Ferrel_Tropopause, va_Ferrel_Tropopause_45, va_Polar_Tropopause, va_Polar_Tropopause_75, va_Hadley_SL, va_Hadley_SL_15, va_Ferrel_SL, va_Ferrel_SL_45, va_Polar_SL, va_Polar_SL_75;
		double wa_Ferrel_Tropopause, wa_Polar_Tropopause, wa_Ferrel_SL, wa_Polar_SL;
		double wa_Hadley_SL, wa_Hadley_Tropopause;
		double wa_equator_Tropopause, wa_equator_SL, va_equator_Tropopause, va_equator_SL, trop_coeff;
		double d_i, d_i_max, d_i_half, d_j, d_j_half, d_j_max, d_k, pi180, d_j_w;
		double d_j_5n, d_j_15n, d_j_45n, d_j_75n, d_j_30n, d_j_5s, d_j_15s, d_j_45s, d_j_75s, d_j_30s, d_j_60n, d_j_60s, d_j_90n, d_j_90s, d_diff;
		double factor_velocity;
		double max_u, max_v, max_w, residuum_u, residuum_v, residuum_w;
		double t_cretaceous, t_cretaceous_coeff, t_cretaceous_max, t_360;
		double j_par_f, j_pol_f, cc, e, j_d, t_d, k_par_f, k_pol_f;
		double *jm_temp_asym, *alfa, *beta;
		double water_wind, t_pol;
		double rR, rg, jmkm, u_sum, v_sum, w_sum, t_sum, c_sum;
		double radiation_ocean_coeff, radiation_land_coeff, radiation_ocean;
		double D, H, G, R_short, Q_short, AG, A, R_long, Q_long, Q_total, Q_lat, Q_sen, Q_bot, Q_rad, t_rad;
		double g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_extra, cp_l, r_0_air, ozean_land, p_baro, L_atm, c13, c43;
		double R_Air, r_h, r_0_water_vapour, R_WaterVapour, precipitablewater_average, precipitation_average, precipitation_NASA_average;
		double ik, eps, epsilon_extra, c_ocean, c_land, t_average, co2_average, co2_pole, gam, t_Ik;
		double radiation_pole, radiation_equator, t_land, t_land_corr;
		double albedo_coeff, albedo_equator, albedo_pole, epsilon_tropopause;
		double ik_coeff, ik_equator, ik_pole;
		double aa, bb, dd, f;
		double *AA, *BB, **CC, CCC, DDD;
		double epsilon_coeff, epsilon_pole, epsilon_average, rad;
		double e_h, a_h, p_h, q_h, t_tau_h, t_Celsius, t_Celsius_ni, t_Celsius_pi, t_Celsius_nj, t_Celsius_pj, t_Celsius_nk, t_Celsius_pk, dp_hdr, dp_hdthe, dp_hdphi;
		double sinthe, sinthe2, kro, lv, ls, coeff_lv, coeff_ls, p_0_vapour, r_0, t_Boussinesq, t_Boussinesq_diff, dqdr, dqdthe, dqdphi, dpdr_c, dpdthe_c, dpdphi_c, dpdr_co2, dpdthe_co2, dpdphi_co2;
		double dudr, dudthe, dudphi, dvdr, dvdthe, dvdphi, dwdr, dwdthe, dwdphi, drdr, drdthe, drdphi, E_dEdr_Rain, E_dEdr_Rain_super, E_dEdr_Ice, E_dEdthe_Rain, E_dEdthe_Rain_super, E_dEdthe_Ice, E_dEdphi_Rain, E_dEdphi_Rain_super, E_dEdphi_Ice, dRaindr, dRaindthe, dRaindphi, dRain_superdr, dRain_superdthe, dRain_superdphi, dIcedr, dIcedthe, dIcedphi;
		double dt, dr, dthe, dphi, dt2, dr2, dthe2, dphi2, rm2;
		double rm, costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2, rmtanthe;
		double E, E_Rain_SL, E_Rain, E_Rain_super, E_Ice, q_Rain, q_Rain_super, q_Ice, E_Rain_super_SL, E_Ice_SL, q_Rain_SL, q_Rain_super_SL, q_Ice_SL, e_1, eps_1, eps_ad, K, e_0;
		double c12, c32, c42, t_Celsius_0, t_Celsius_1, t_Celsius_2;

		char Temperature_West[50], Salinity_West[50], Temperature_East[50], Salinity_East[50];
 
		string time_slice_comment, time_slice_number, time_slice_unit;
		string temperature_comment, temperature_gain, temperature_modern, temperature_average, temperature_unit, temperature_cretaceous, temperature_average_cret;
		string co2_comment, co2_gain, co2_modern, co2_av, co2_unit, co2_cretaceous_str, co2_average_cret, co2_average_str;


	public:
		double t_equator, t_tropopause, t_coeff, t_pole, c_equator, c_tropopause, c_coeff, c_pol, p_equator, p_tropopause, p_coeff, p_pol;
		double co2_equator, co2_tropopause, co2_coeff, co2_pol, co2_vegetation, co2_ocean, co2_land, co2_cretaceous;

		BC_Thermo ( int, int, int, int, int, int, int, int, int, int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double );
		~BC_Thermo();


		void IC_CellStructure ( int *, Array &, Array &, Array & );

		void IC_WestEastCoast ( double, Array &, Array &, Array &, Array &, Array & );

		void BC_Temperature ( int *, Array &, Array &, Array &, Array & );

		void TropopauseLocation ( int * );

		void BC_Radiation_two_layer ( int *, double, double, double, double, double, double, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_Radiation_multi_layer ( int *, double, double, double, double, double, double, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_Radiation_parabolic ( Array_2D &, Array & );

		void BC_WaterVapour ( int *, Array &, Array &, Array &, Array &, Array_2D & );

		void BC_CO2 ( int *, Array_2D &, Array &, Array &, Array &, Array & );

		void BC_Surface_Temperature ( const string &, Array & );

		void BC_Surface_Precipitation ( const string &, Array_2D & );

		void BC_Pressure ( int *, Array &, Array &, Array & );

		void Latent_Heat ( double, double, double, double, double, double, double, double, double, double, double, double, double, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		double out_temperature (  ) const;

		double out_co2 (  ) const;

};
#endif
