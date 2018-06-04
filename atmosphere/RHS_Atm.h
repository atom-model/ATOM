/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include <iostream>
#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"

#ifndef _RHS_ATMOSPHERE_
#define _RHS_ATMOSPHERE_

using namespace std;


class RHS_Atmosphere
{
	private:
		int i, j, k, im, jm, km;
		int i_level, h_check_i, h_check_j, h_check_k;
		int i_T, i_b;
		int **im_tropopause;

		double zeit, dt, dr, dthe, dphi;
		double re, pr, sc_WaterVapour, sc_CO2, g, WaterVapour, Buoyancy, CO2, lambda;
		double E_Rain_SL, E_Rain, E_Rain_super, E_Ice;
		double e_SL, a_SL, p_SL, q_SL, t_tau_SL, t_Celsius_SL, exp_pressure, t_u;
		double e_h, a_h, p_h, q_h, t_tau_h, t_Celsius, dp_hdr, dp_hdthe, dp_hdphi;
		double h_level, h_h, sat_Deficit, RF_e, RF_c, Evaporation_Haude;
		double a2, dt2, dr2, dthe2, dphi2, rm2;
		double kr1, kr2, kthe1, kthe2, rm, c43, c13;
		double tanthe, kpr, kpthe, kpphi, kphi;
		double RS_WaterVapour_Balance_Rain, RS_WaterVapour_Balance_Ice, RS_co2_Balance_co2;
		double RS_buoyancy_Energy, RS_buoyancy_Momentum, RS_buoyancy_Water_Vapour; 
		double RS_LatentHeat_Energy_Rain, RS_LatentHeat_Energy_Rain_super, RS_LatentHeat_Energy_Ice, Rain_aux, Rain_super_aux, Ice_aux;
		double q_i_plus, q_i_minus, q_j_plus, q_j_minus, q_k_plus, q_k_minus, q_Rain, q_Rain_super, q_Ice;
		double sinthe, sinthe2, kro, lv, ls, ep, hp, coeff_trans, coeff_energy, coeff_buoy, p_0, p_in_vapour, r_0, t_0, cp_l, R_Air, R_WaterVapour, dqdr, dqdthe, dqdphi, dpdr_c, dpdthe_c, dpdphi_c, dpdr_co2_, dpdthe_co2_, dpdphi_co2_;
		double costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2, rmtanthe;
		double resra, resre, dummy_1, dummy_2, dummy_3, k_Force;
		double dudr, dudthe, dudphi, dvdr, dvdthe, dvdphi, dwdr, dwdthe, dwdphi, drdr, drdthe, drdphi, E_dEdr_Rain, E_dEdr_Rain_super, E_dEdr_Ice, E_dEdthe_Rain, E_dEdthe_Rain_super, E_dEdthe_Ice, E_dEdphi_Rain, E_dEdphi_Rain_super, E_dEdphi_Ice, dRaindr, dRaindthe, dRaindphi, dRain_superdr, dRain_superdthe, dRain_superdphi, dIcedr, dIcedthe, dIcedphi;
		double d2udr2, d2udthe2, d2udphi2, d2vdr2, d2vdthe2, d2vdphi2, d2wdr2, d2wdthe2, d2wdphi2, d2rdr2, d2rdthe2, d2rdphi2;
		double d2udrdthe, d2udrdphi, d2vdrdthe, d2vdrdphi, d2wdrdthe, d2wdrdphi, d2rdrdthe, d2rdrdphi, d2udthedphi, d2vdthedphi, d2wdthedphi, d2rdthedphi;
		double d2udrdt, d2vdrdt, d2wdrdt, d2rdrdt, d2udthedt, d2vdthedt, d2wdthedt, d2rdthedt, d2udphidt, d2vdphidt, d2wdphidt, d2rdphidt;
		double dudt, dvdt, dwdt, drdt, dtdr, dpdr, dpdr_atm, dtdthe, dpdthe, dtdphi, dpdphi, d2tdr2, d2tdthe2, d2tdphi2;
		double dcdr, dcdthe, dcdphi, d2cdr2, d2cdthe2, d2cdphi2, c_0;
		double dclouddr, dclouddthe, dclouddphi, d2clouddr2, d2clouddthe2, d2clouddphi2;
		double dicedr, dicedthe, dicedphi, d2icedr2, d2icedthe2, d2icedphi2;
		double dcodr, dcodthe, dcodphi, d2codr2, d2codthe2, d2codphi2, co2_0;
		double h_0_0, h_0_i, h_c_i, h_d_i, h_0_j, h_c_j, h_d_j, h_0_k, h_c_k, h_d_k, cc; 
		double d_i, d_i_max, t_tropopause, j_half, d_j_half, trop_co2_eff, d_j, i_max, i_beg, gam, sigma;
		double alf_ev, alf_dep, eps_T;
		double  S_c_c, S_au, S_nuc, S_ac, S_rim, S_shed, S_ev, S_dep, S_i_dep, S_melt, S_if_frz, S_cf_frz, S_r_frz, S_c_frz, S_c_au, S_i_au, S_d_au, S_agg, S_i_cri, S_r_cri, S_s_dep, S_i_melt, S_s_melt;
		double S_v, S_c, S_i, S_r, S_s;
		double a_if, c_ac, c_rim, bet_ev, alf_melt, bet_melt, bet_dep, bet_s_dep, alf_if, alf_cf, a_s_melt, b_s_melt, E_cf, N_cf, N_cf_0, N_cf_0_surf, N_cf_0_6km;
		double tau_r, tau_s, t_1, t_2, t_m1, t_m2, t_r_frz, c_r_frz;
		double a_mc, a_mv, a_m, a_i_m, a_s_m, N_r_0, N_s_0;
		double *cloud_max;
		double r_q_r, r_q_c, r_q_i, r_q_s, r_q_v, r_dry, r_humid;
		double D_r, D_s, D_i, m_r, m_s, v_r_0, v_s_0, v_rp_T, v_sp_T, R_r;
		double c78, hight, dist;
		double N_i_0, N_i, t_nuc, t_d, t_hn, m_i, m_i_0, m_i_max, m_s_0, c_i_dep, c_c_au, c_i_au, c_agg, c_i_cri, c_r_cri, c_s_dep, c_s_melt;
		double MC_s, MC_q_v, MC_v, MC_w, s, cu, e_l, e_p, g_p, C_p, u_d1, w_d1, w_d, v_d1, v_d, q_v_d1, q_v_d, e_d, s_d1, s_d, u_u1, w_u1, w_u, v_u1, v_u, q_v_u1, q_v_u, q_c_u1, q_c_u, s_u1, s_u, M_d, E_d, D_d, M_u, E_u, D_u;
		double b_u, alf_1, alf_2, p_ps, bet_p, eps_u, delta_i_c, K_p, del_u, cloud_u, c_u, p_t_in, t_Celsius_0, E_Rain_t_in, q_Rain_t_in;

	public:
		RHS_Atmosphere ( int, int, double, double, double );
		RHS_Atmosphere ( int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double );
		~RHS_Atmosphere ();

		void RK_RHS_3D_Atmosphere ( int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D & );

		void RK_RHS_2D_Atmosphere ( int, int, double, double, double, double, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );
};
#endif
