/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
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

#ifndef Results_MSL_Atm_Atm
#define Results_MSL_Atm_Atm

using namespace std;

class Results_MSL_Atm
{
	private:
		int i, j, k, l, im, jm, km, sun, h_point_max, h_land, h_ocean;
		int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg, max_ch, min_ch;
		int iter_prec;
		int *i_T, i_T_max, i_T_a, i_T_b, i_b, i_LFS, i_parc;

		double E, E_Rain_SL, E_Rain, E_Rain_super, E_Ice, q_Rain, q_Rain_n, q_Rain_super, q_Ice, q_Ice_n, E_Rain_super_SL, E_Ice_SL, q_Rain_SL, q_Rain_super_SL, q_Ice_SL;
		double e, a, e_SL, a_SL, p_SL, q_SL, t_dew_SL, t_Celsius_SL, q_h_p, q_T;
		double e_h, a_h, p_h, q_h, t_dew, t_Celsius, t_Celsius_1, c_grad, t_grad, t_denom, Delta, E_a, gamma, g, gam;
		double i_level, h_level, h_h, sat_deficit, RF_e, RF_c, Evap_Haude, Evaporation_Penman_average, Evaporation_Haude_average;
		double ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_extra, lv, ls, cp_l, r_0_air, dt, dr, ozean_land, p_baro, L_atm, c13, c43;
		double R_Air, r_h, r_0_water, r_0_water_vapour, R_WaterVapour, precipitablewater_average, precipitation_average, precipitation_NASA_average;
		double coeff_mmWS, coeff_lv, coeff_ls, coeff_Diffusion_latent, coeff_Diffusion_sensibel, f_Haude, co2_vegetation, co2_ocean, co2_land, prec, r_humid, r_dry, r_dry_0;
		double Value_1, Value_2, Value_3, Value_4, Value_5, Value_6, Value_7, Value_8, Value_10, Value_12, Value_13, Value_14, Value_15, Value_16, Value_17, Value_18, Value_19;
		double rm, rmsinthe, sinthe, dthe, dphi, dtdthe, dtdphi;
		double alf_ev, alf_dep, eps_T;
		double  S_c_c, S_au, S_nuc, S_ac, S_rim, S_shed, S_ev, S_dep, S_i_dep, S_melt, S_if_frz, S_cf_frz, S_frz, S_c_frz, S_c_au, S_i_au, S_d_au, S_agg, S_i_cri, S_r_cri, S_s_dep, S_i_melt, S_s_melt;
		double S_v, S_c, S_i, S_r, S_s;
		double a_if, c_ac, c_rim, bet_ev, alf_melt, bet_melt, bet_dep, bet_s_dep, alf_if, alf_cf, a_s_melt, b_s_melt, E_cf, N_cf, N_cf_0, N_cf_0_surf, N_cf_0_6km;
		double tau_r, tau_s, t_1, t_2, t_m1, t_m2;
		double a_mc, a_mv, a_m, coeff_P, a_i_m, a_s_m, N_r_0, N_s_0;
		double r_q_r, r_q_c, r_q_i, r_q_s, r_q_v;
		double D_r, D_s, D_i, m_r, m_s, v_r_0, v_s_0, v_rp_T, v_sp_T;
		double N_i_0, N_i, t_nuc, t_d, t_hn, m_i, m_i_0, m_i_max, m_s_0, c_i_dep, c_c_au, c_i_au, c_agg, c_i_cri, c_r_cri, c_s_dep, c_s_melt;
		double MC_s, MC_q_v, MC_v, MC_w, s, cu, C_p, u_b, u_d1, w_d1, w_d, v_d1, v_d, q_v_d1, q_v_d, s_d1, s_d, u_u1, w_u1, w_u, v_u1, v_u, q_v_u1, q_v_u, q_v_b, q_i_b, q_c_u1, q_c_u, q_c_d, q_c_b, s_u1, s_u, M_d, E_d, D_d, M_u, E_u, D_u, r_dry_parc, r_humid_parc, a_w, b_w, a_i, b_i;
		double a_u, a_d, b_u, alf_1, alf_2, p_ps, bet_p, eps_u, delta_i_c, K_p, del_u, cloud_u, dt_u, t_u, buoyancy_check, u_u1_check, u_d1_check, sat_check, u_max_check, t_b, s_b, c_b, cloud_b, v_b, w_b, gam_d, dcdr, dcdthe, dcdphi, eps_d, del_d, T, T_nue, T_tilda_h, t_Celsius_ni, t_Celsius_pi, E_dEdr_Rain, E_dEdr_Ice, T_p, T_n, p_h_p, p_h_n, p_v_sw_p, p_v_sw_n, q_v_sw_p, q_v_sw_n, dq_v_swdT, p_v_0, p_v_hyp, q_v_hyp, q_v_hyp_n, p_v_sw, q_v_sw, dq_v_sidT, p_v_si, q_v_si, Rd_Rv, Rd_Rv_1, T_t_0, T_b_w, T_t_0_T_b_w, dp_v_swdT, T_b_i, T_t_0_T_b_i, dp_v_sidT, dp_dT, p_denom, p_nom;
		double *e_d, *e_l, *e_p, *g_p, *c_u, *r_humid_u, *r_humid_u_parc, *r_dry_u, *r_dry_u_parc;
		double T_t_00, T_t_01, DEP, CND,d_q_v, d_q_c, d_q_i, d_t, p_t_0, t_Celsius_0, E_Rain_t_0, q_Rain_t_0;

		string name_Value_1, name_Value_2, name_Value_3, name_Value_4, name_Value_5, name_Value_6, name_Value_7, name_Value_8, name_Value_9, name_Value_10, name_Value_11, name_Value_12, name_Value_13, name_Value_14, name_Value_15, name_Value_16, name_Value_17, name_Value_18, name_Value_19, name_Value_20, name_Value_21, name_unit_wm2, name_unit_mm, name_unit_mmd, name_unit_mma;
		string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading, heading_Dresden, heading_Sydney, heading_Equator;

	public:
		Results_MSL_Atm ( int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double ); 


		~Results_MSL_Atm (  );

		void run_MSL_data ( int, int, int, double, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void show_MSL_data ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

};
#endif
