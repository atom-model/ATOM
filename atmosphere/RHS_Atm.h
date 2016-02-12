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

#ifndef _RHS_ATMOSPHERE_
#define _RHS_ATMOSPHERE_

using namespace std;



class RHS_Atmosphere
{
	private:
		int i, j, k, im, jm, km;
		int i_level, h_check_i, h_check_j, h_check_k;

		double zeit, dt, dr, dthe, dphi;
		double re, pr, ec, sc_WaterVapour, sc_CO2, gr, omega, coriolis, centrifugal, WaterVapour, buoyancy, CO2;
		double E_Rain_SL, E_Rain, E_Rain_super, E_Ice;
		double e_SL, a_SL, p_SL, q_SL, t_tau_SL, t_Celsius_SL;
		double e_h, a_h, p_h, q_h, t_tau_h, t_Celsius;
		double h_level, h_h, sat_Deficit, RF_e, RF_c, Evaporation_Haude;
		double a2, dt2, dr2, dthe2, dphi2, rm2;
		double kr1, kr2, kthe1, kthe2, rm, c43, c13;
		double tanthe, kpr, kpthe, kpphi, kphi;
		double RS_Coriolis_Energy, RS_centrifugal_Energy;
		double RS_WaterVapour_Balance_Rain, RS_WaterVapour_Balance_Ice, RS_co2_Balance_co2;
		double RS_buoyancy_Energy, RS_buoyancy_Momentum, RS_buoyancy_Water_Vapour; 
		double RS_LatentHeat_Energy_Rain, RS_LatentHeat_Energy_Rain_super, RS_LatentHeat_Energy_Ice, Rain_aux, Rain_super_aux, Ice_aux;
		double RS_Coriolis_Momentum_rad, RS_Coriolis_Momentum_the, RS_Coriolis_Momentum_phi, RS_centrifugal_Momentum_rad, RS_centrifugal_Momentum_the;
		double q_i_plus, q_i_minus, q_j_plus, q_j_minus, q_k_plus, q_k_minus, q_Rain, q_Rain_super, q_Ice;
		double sinthe, sinthe2, kro, lv, ls, ep, hp, coeff_lv, coeff_ls, p_0, p_0_vapour, r_0, t_0, t_Boussinesq, cp_l, R_Air, R_WaterVapour, dqdr, dqdthe, dqdphi, dpdr_c, dpdthe_c, dpdphi_c, dpdr_co2, dpdthe_co2, dpdphi_co2;
		double costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2, rmtanthe;
		double resra, resre, dummy_1, dummy_2, dummy_3, k_Force;
		double dudr, dudthe, dudphi, dvdr, dvdthe, dvdphi, dwdr, dwdthe, dwdphi, drdr, drdthe, drdphi;
		double d2udr2, d2udthe2, d2udphi2, d2vdr2, d2vdthe2, d2vdphi2, d2wdr2, d2wdthe2, d2wdphi2, d2rdr2, d2rdthe2, d2rdphi2;
		double d2udrdthe, d2udrdphi, d2vdrdthe, d2vdrdphi, d2wdrdthe, d2wdrdphi, d2rdrdthe, d2rdrdphi, d2udthedphi, d2vdthedphi, d2wdthedphi, d2rdthedphi;
		double d2udrdt, d2vdrdt, d2wdrdt, d2rdrdt, d2udthedt, d2vdthedt, d2wdthedt, d2rdthedt, d2udphidt, d2vdphidt, d2wdphidt, d2rdphidt;
		double dudt, dvdt, dwdt, drdt, dtdr, dpdr, dpdr_atm, dtdthe, dpdthe, dtdphi, dpdphi, d2tdr2, d2tdthe2, d2tdphi2;
		double dcdr, dcdthe, dcdphi, d2cdr2, d2cdthe2, d2cdphi2, c_0;
		double dco2dr, dco2dthe, dco2dphi, d2co2dr2, d2co2dthe2, d2co2dphi2, co2_0;
		double h_0_i, h_c_i, h_d_i, h_0_j, h_c_j, h_d_j, h_0_k, h_c_k, h_d_k, cc; 

	public:
		RHS_Atmosphere ( int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double );
		~RHS_Atmosphere ();


		void RK_RHS_Atmosphere ( int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void RK_RHS_2D_Atmosphere ( int, int, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );
};
#endif
