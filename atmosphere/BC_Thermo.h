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
		int i, j, k, im, jm, km, k_half, j_half, i_half, i_max, j_max, k_max, i_beg, im_1, i_land, iter_rad;
		int j_aeq, j_pol_n, j_pol_s, j_pol_v_n, j_pol_v_s, j_fer_n, j_fer_s, j_fer_v_n, j_fer_v_s, j_had_n, j_had_s, j_had_v_n, j_had_v_s;
		int j_had_n_end, j_had_s_end, k_w, k_w_end, k_e;
		int j_n, j_s;
		int sun, Ma, Ma_max, Ma_max_half;
		int j_par, j_pol, k_par, k_pol;
		int k_a, k_b, flip, k_grad, k_grad_init, i_middle;
		int k_water, k_sequel;
		int n_smooth;
		int j_r, k_r, j_sun;
		int RadiationModel, sun_position_lat, sun_position_lon, declination;
		int iter_prec; 
		int *im_tropopause;

		double d_k_half, d_k_max, t_eff_earth, r, rad_bal_minus; 
		double dummy_1, dummy_2, dummy_3;
		double ca, ua_00, ua_30, ua_60, ua_90, u_co2_eff, u_end;
		double va_Hadley_Tropopause, va_Hadley_Tropopause_15, va_Ferrel_Tropopause, va_Ferrel_Tropopause_45, va_Polar_Tropopause, va_Polar_Tropopause_75, va_Hadley_SL, va_Hadley_SL_15, va_Ferrel_SL, va_Ferrel_SL_45, va_Polar_SL, va_Polar_SL_75;
		double wa_Ferrel_Tropopause, wa_Polar_Tropopause, wa_Ferrel_SL, wa_Polar_SL;
		double wa_Hadley_SL, wa_Hadley_Tropopause;
		double wa_equator_Tropopause, wa_equator_SL, va_equator_Tropopause, va_equator_SL, trop_co2_eff;
		double d_i, d_i_max, d_i_half, d_j, d_j_half, d_j_max, d_k, pi180, d_j_w;
		double d_j_5n, d_j_15n, d_j_45n, d_j_75n, d_j_30n, d_j_5s, d_j_15s, d_j_45s, d_j_75s, d_j_30s, d_j_60n, d_j_60s, d_j_90n, d_j_90s, d_diff;
		double max_u, max_v, max_w;
		double t_cretaceous, t_cretaceous_co2_eff, t_cretaceous_max, t_360;
		double j_par_f, j_pol_f, cc, e, j_d, t_d, k_par_f, k_pol_f;
		double water_wind, t_pol;
		double rR, rg, jmkm, u_sum, v_sum, w_sum, t_sum, c_sum;
		double radiation_ocean_co2_eff, radiation_land_co2_eff, radiation_ocean;
		double D, H, G, R_short, Q_short, AG, A, R_long, Q_long, Q_total, Q_lat, Q_sen, Q_bot, Q_rad, t_rad;
		double g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_extra, cp_l, r_air, ozean_land, L_atm, c13, c43;
		double R_Air, r_h, r_water_vapour, R_WaterVapour, precipitablewater_average, precipitation_average, precipitation_NASA_average;
		double ik, eps, epsilon_extra, c_ocean, c_land, t_average, co2_average, co2_pole, gam, t_Ik, Ik_loss, Ik_tot;
		double radiation_pole, radiation_equator, t_land;
		double albedo_co2_eff, albedo_equator, albedo_pole, epsilon_tropopause;
		double ik_co2_eff, ik_equator, ik_pole;
		double aa, bb, dd, f;
		double epsilon_co2_eff, epsilon_pole, epsilon_average;
		double e_h, a_h, p_h, q_h, t_tau_h, t_Celsius, t_Celsius_ni, t_Celsius_pi, t_Celsius_nj, t_Celsius_pj, t_Celsius_nk, t_Celsius_pk, dp_hdr, dp_hdthe, dp_hdphi;
		double sinthe, sinthe2, lv, ls, coeff_lv, coeff_ls, r_0, dqdr, dqdthe, dqdphi, dpdr_c, dpdthe_c, dpdphi_c, dpdr_co2_, dpdthe_co2_, dpdphi_co2_;
		double dudr, dudthe, dudphi, dvdr, dvdthe, dvdphi, dwdr, dwdthe, dwdphi, drdr, drdthe, drdphi, E_dEdr_Rain, E_dEdr_Rain_super, E_dEdr_Ice, E_dEdthe_Rain, E_dEdthe_Rain_super, E_dEdthe_Ice, E_dEdphi_Rain, E_dEdphi_Rain_super, E_dEdphi_Ice, dRaindr, dRaindthe, dRaindphi, dRain_superdr, dRain_superdthe, dRain_superdphi, dIcedr, dIcedthe, dIcedphi;
		double dt, dr, dthe, dphi, dt2, dr2, dthe2, dphi2, rm2;
		double rm, costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2;
		double E, E_Rain_SL, E_Rain, E_Rain_super, E_Ice, q_Rain, q_Rain_super, q_Ice;
		double c12, c32, c42, t_Celsius_0, t_Celsius_1, t_Celsius_2;
		double TK, rad_lon_terrestic, rad_lon_back;
		double epsilon_tropo, epsilon_co2_eff_max;
		double r_dry, r_humid, p_SL, t_SL;
		double t_1, t_2, t_u, t_Celsius_SL, t_dew, t_dew_SL, T, T_nue, T_tilda_h, T_t_in, T_t_in0, T_t_in1, q_T, q_Rain_n, e_SL, a_SL, h_level, i_level, h_h, sat_deficit, RF_e;
		double q_v_b, q_c_b, q_i_b, q_v_hyp, q_v_hyp_n, CND, DEP, d_q_v, d_q_c, d_q_i, d_t, q_Ice_n;
		double t_equator, t_tropopause, t_co2_eff, t_pole, c_equator, c_tropopause;
		double co2_equator, co2_tropopause, co_co2_eff, co_pol, co2_vegetation, co2_ocean, co2_land, co2_cretaceous;
		double *AA, **CC, CCC, DDD, *cloud_max;
		double *jm_temp_asym, *alfa, *beta;
 
		string time_slice_co2_mment, time_slice_number, time_slice_unit;
		string temperature_co2_mment, temperature_gain, temperature_modern, temperature_average, temperature_unit, temperature_cretaceous, temperature_average_cret;
		string co_co2_mment, co_gain, co_modern, co_av, co_unit, co_cretaceous_str, co_average_cret, co_average_str;


	public:
		BC_Thermo ( int, int, int, int, int, int, int, int, int, int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double );
		~BC_Thermo();


		void IC_CellStructure ( int *, Array &, Array &, Array & );

		void BC_Temperature ( Array &, Array &, Array &, Array & );

		void TropopauseLocation ( int * );

		void BC_Radiation_2D_layer ( Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array &, Array &, Array &, Array &, Array & );

		void BC_Radiation_multi_layer ( int, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_Radiation_parabolic ( Array_2D &, Array & );

		void BC_WaterVapour ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D & );

		void IC_Cloud_Ice ( Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_CO2 ( Array_2D &, Array &, Array &, Array &, Array & );

		void BC_Surface_Temperature ( const string &, Array & );

		void BC_Surface_Precipitation ( const string &, Array_2D & );

		void BC_Pressure ( Array &, Array &, Array & );

		void Latent_Heat ( Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		double out_temperature (  ) const;

		double out_co2 (  ) const;

};
#endif
