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
		int i, j, k, im, jm, km, k_half, j_half, i_half, i_max, j_max, k_max, i_beg, im_1, i_land, iter_rad, ll;
		int j_aeq, j_pol_n, j_pol_s, j_pol_v_n, j_pol_v_s, j_fer_n, j_fer_s, j_fer_v_n, j_fer_v_s, j_had_n, j_had_s, j_had_v_n, j_had_v_s;
		int j_had_n_end, j_had_s_end, k_w, k_w_end, k_e;
		int j_n, j_s;
		int sun, Ma, Ma_max, Ma_max_half;
		int j_par, j_pol, k_par, k_pol;
		int k_a, k_b, flip, j_grad, k_grad, k_grad_init, i_middle;
		int j_water, k_water, j_sequel, k_sequel;
		int n_smooth;
		int j_r, k_r, j_sun;
		int RadiationModel, sun_position_lat, sun_position_lon, declination, NASATemperature;
		int iter_prec; 
		int *im_tropopause;

		double d_k_half, d_k_max; 
		double dummy_1, dummy_2, dummy_3;
		double ca, ua_00, ua_30, ua_60, ua_90;
		double va_Hadley_Tropopause, va_Hadley_Tropopause_15, va_Ferrel_Tropopause, va_Ferrel_Tropopause_45, va_Polar_Tropopause, va_Polar_Tropopause_75, va_Hadley_SL, va_Hadley_SL_15, va_Ferrel_SL, va_Ferrel_SL_45, va_Polar_SL, va_Polar_SL_75;
		double wa_Ferrel_Tropopause, wa_Polar_Tropopause, wa_Ferrel_SL, wa_Polar_SL;
		double wa_Hadley_SL, wa_Hadley_Tropopause;
		double wa_equator_Tropopause, wa_equator_SL, va_equator_Tropopause, va_equator_SL, trop_co2_eff;
		double d_i, d_i_max, d_i_half, d_j, d_j_half, d_j_max, d_k, pi180, d_j_w;
		double d_j_5n, d_j_15n, d_j_45n, d_j_75n, d_j_30n, d_j_5s, d_j_15s, d_j_45s, d_j_75s, d_j_30s, d_j_60n, d_j_60s, d_j_90n, d_j_90s, d_diff;
		double t_cretaceous, t_cretaceous_co2_eff, t_cretaceous_max, t_cret_cor, t_360;
		double j_par_f, j_pol_f, e, j_d, t_dd, k_par_f, k_pol_f;
		double g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, cp_l, r_air, L_atm, c13, c43;
		double R_Air, r_h, r_water_vapour, R_WaterVapour, precipitablewater_average, precipitation_average, precipitation_NASA_average;
		double eps, c_ocean, t_land, c_land, c_coeff, t_average, co2_average, co2_pole, gam, t_Ik, Ik_loss, Ik_tot;
		double albedo_co2_eff, albedo_equator, albedo_pole;
		double ik_co2_eff, ik_equator, ik_pole;
		double aa, bb, cc, dd, f;
		double epsilon_co2_eff_2D, epsilon_co2_eff, epsilon_pole, epsilon_equator, epsilon_tropopause, epsilon_co2_eff_max;

		double e_h, a_h, p_h, q_h, t_tau_h, t_Celsius, t_Celsius_ni, t_Celsius_pi, t_Celsius_nj, t_Celsius_pj, t_Celsius_nk, t_Celsius_pk, dp_hdr, dp_hdthe, dp_hdphi;
		double sinthe, sinthe2, lv, ls, coeff_lv, coeff_ls, coeff_L_atm_u_0, r_0;
		double dt, dt_dim, dt_rain_dim, dt_snow_dim, dr, dthe, dphi, dt2, dr2, dthe2, dphi2, rm2;
		double rm, costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2;
		double E, E_Rain_SL, E_Rain, E_Rain_super, E_Ice, q_Rain, q_Rain_super, q_Ice;
		double c12, c32, c42, t_Celsius_it, t_Celsius_0, t_Celsius_1, t_Celsius_2;
		double TK, rad_lon_terrestic, rad_lon_back;
		double r_dry, r_humid, p_SL, t_SL, exp_pressure, hight;
		double t_u, t_Celsius_SL, t_dew, t_dew_SL, T, T_nue, T_it, q_T, q_Rain_n;
		double q_v_b, q_c_b, q_i_b, q_v_hyp, CND, DEP, d_q_v, d_q_c, d_q_i, d_t, q_Ice_n;
		double t_equator, t_tropopause, t_co2_eff, t_pole, c_equator, c_tropopause, coeff_mmWS;
		double co2_equator, co2_tropopause, co_co2_eff, co_pol, co2_vegetation, co2_ocean, co2_land, co2_cretaceous;
		double *AA, **CC, CCC, DDD;
		double *jm_temp_asym, *alfa, *beta;

		double S_c_c, S_au, S_nuc, S_ac, S_rim, S_shed, S_ev, S_dep, S_i_dep, S_melt, S_if_frz, S_cf_frz, S_r_frz, S_c_frz, S_c_au, S_i_au, S_d_au, S_agg, S_i_cri, S_r_cri, S_s_dep, S_i_melt, S_s_melt;
		double S_v, S_c, S_i, S_r, S_s, P_rain_n, P_snow_n, P_rain_n_o, P_snow_n_o;
		double a_if, c_ac, c_rim, bet_ev, alf_melt, bet_melt, bet_dep, bet_s_dep, alf_if, alf_cf, a_s_melt, b_s_melt, E_cf, N_cf, N_cf_0, N_cf_0_surf, N_cf_0_6km;
		double tau_r, tau_s, t_1, t_00, t_m1, t_m2, t_r_frz, c_r_frz, alf_ev, alf_dep, a_d, b_u, alf_1, alf_2, p_ps, bet_p, p_t_in, E_Rain_t_in, q_Rain_t_in;
		double a_mc, a_mv, a_m, coeff_P, a_i_m, a_s_m, N_r_0, N_s_0;
		double N_i_0, N_i, t_nuc, t_d, t_hn, m_i, m_i_0, m_i_max, m_s_0, c_i_dep, c_c_au, c_i_au, c_agg, c_i_cri, c_r_cri, c_s_dep, c_s_melt;
		double u_av, v_av, w_av, coeff_L, coeff_Q;

 
 
		string time_slice_comment, time_slice_number, time_slice_unit, output_path;
		string temperature_comment, temperature_gain, temperature_modern, temperature_average, temperature_unit, temperature_cretaceous, temperature_average_cret;
		string co_comment, co_gain, co_modern, co_av, co_unit, co_cretaceous_str, co_average_cret, co_average_str;


	public:
		BC_Thermo ( string &, int, int, int, int, int, int, int, int, int, int, int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double );
		~BC_Thermo();


		void IC_CellStructure ( int *, Array &, Array &, Array &, Array & );

		void BC_Temperature ( int *, Array_2D &, Array &, Array &, Array &, Array & );

		void TropopauseLocation ( int * );

        void BC_Radiation_2D_layer(Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, 
                                   Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, 
                                   Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, 
                                   Array &, Array &, Array &, Array &, Array & );

        void BC_Radiation_multi_layer(int*, int, Array_2D &, Array_2D &, Array_2D &, Array_2D &, 
                                      Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, 
                                      Array_2D &, Array &, Array &, Array &, const Int3DArray &h, 
                                      Array &, Array &, Array &, Array & );

		void BC_WaterVapour ( Array &, Array &, Array & );

		void BC_CloudWaterIce ( Array &, Array &, Array &, Array & );

		void Ice_Water_Saturation_Adjustment ( int *, int, int, int, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void Two_Category_Ice_Scheme ( int, int, int, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_CO2 ( Array_2D &, Array &, Array &, Array &, Array & );

		void BC_NASAbasedSurfTempRead ( const string &, double &, double &, Array &, Array &, Array &, Array & );

		void BC_NASAbasedSurfTempWrite ( const string &, double &, double &, Array &, Array &, Array &, Array & );

		void BC_Surface_Temperature_NASA ( const string &, Array_2D &, Array & );

		void BC_Surface_Precipitation_NASA ( const string &, Array_2D & );

		void BC_Pressure ( Array &, Array &, Array &, Array & );

		void Latent_Heat ( Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void IC_Temperature_WestEastCoast ( Array &, Array & );

		double exp_func ( double &, const double &, const double & );

		double cloud_ice ( const double &, int, const int & );

		double out_t_cretaceous (  ) const;

		double out_t_cret_cor (  ) const;

		double out_co2 (  ) const;

};
#endif
