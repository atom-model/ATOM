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
		int im, jm, km;
		int h_check_i, h_check_j, h_check_k;

		double dt, dr, dthe, dphi;
		double re, pr, ec, sc_WaterVapour, sc_CO2, g, omega, coriolis, centrifugal, WaterVapour, buoyancy, CO2, lambda;
		double E_Rain;
		double E_Ice;
		double p_h, q_h;
		double t_Celsius;
		double dr2, dthe2, dphi2, rm2;
		double rm, c43, c13;
		double RS_Coriolis_Energy, RS_centrifugal_Energy;
		double RS_buoyancy_Momentum, RS_buoyancy_Water_Vapour; 
		double RS_Coriolis_Momentum_rad, RS_Coriolis_Momentum_the, RS_Coriolis_Momentum_phi, RS_centrifugal_Momentum_rad, RS_centrifugal_Momentum_the;
		double q_Rain;
		double q_Ice;
		double sinthe, sinthe2;
		double coeff_lv, coeff_ls, coeff_lf;
		double costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2;
		double k_Force;
		double dudr, dudthe, dudphi, dvdr, dvdthe, dvdphi, dwdr, dwdthe, dwdphi;
		double d2udr2, d2udthe2, d2udphi2, d2vdr2, d2vdthe2, d2vdphi2, d2wdr2, d2wdthe2, d2wdphi2;
		double dtdr, dpdr;
		double dtdthe, dpdthe, dtdphi, dpdphi, d2tdr2, d2tdthe2, d2tdphi2;
		double dcdr, dcdthe, dcdphi, d2cdr2, d2cdthe2, d2cdphi2;
		double dclouddr, dclouddthe, dclouddphi, d2clouddr2, d2clouddthe2, d2clouddphi2;
		double dicedr, dicedthe, dicedphi, d2icedr2, d2icedthe2, d2icedphi2;
		double dcodr, dcodthe, dcodphi, d2codr2, d2codthe2, d2codphi2;
		double h_0_i, h_c_i, h_d_i, h_0_j, h_c_j, h_d_j, h_0_k, h_c_k, h_d_k, cc; 
		double gam, sigma;
		double alf_ev;
		double S_c_c;
		double S_nuc, S_ac, S_rim, S_shed, S_ev;
		double S_i_dep;
		double S_if_frz, S_cf_frz;
		double S_c_frz, S_c_au, S_i_au, S_d_au, S_agg, S_i_cri, S_r_cri, S_s_dep, S_i_melt, S_s_melt;
		double a_if, c_ac, c_rim, bet_ev, alf_melt, bet_melt, bet_dep, bet_s_dep, alf_if, alf_cf, a_s_melt, b_s_melt, E_cf, N_cf, N_cf_0, N_cf_0_surf, N_cf_0_6km;
		double tau_r, tau_s, t_1, t_2, t_m1, t_m2;
		double a_mc, a_mv;
		double a_i_m, a_s_m, N_r_0, N_s_0;
		double r_dry, r_humid;
		double N_i_0, N_i, t_nuc, t_d, t_hn, m_i, m_i_0, m_i_max, m_s_0, c_i_dep, c_c_au, c_i_au, c_agg, c_i_cri, c_r_cri, c_s_dep, c_s_melt;
		double b_u, alf_1, alf_2, p_ps, bet_p;
		double p_t_in, t_Celsius_0, E_Rain_t_in, q_Rain_t_in;

	public:
		RHS_Atmosphere ( int, int, double, double, double, double, double, double );
		RHS_Atmosphere ( int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double );
		~RHS_Atmosphere ();

		void RK_RHS_3D_Atmosphere ( int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void RK_RHS_2D_Atmosphere ( int, int, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );
};
#endif
