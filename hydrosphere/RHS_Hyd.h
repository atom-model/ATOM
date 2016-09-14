/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include <iostream>
#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"

#ifndef _RHS_HYDROSPHERE_
#define _RHS_HYDROSPHERE_

using namespace std;


class RHS_Hydrosphere
{
	private:
		int i, j, k, im, jm, km;
		int h_check_i, h_check_j, h_check_k;

		double dt, dr, dthe, dphi, r0;
		double re, pr, ec, sc, g, omega, coriolis, centrifugal, c_salt;
		double a2, dt2, dr2, dthe2, dphi2, rm2;
		double kr1, kr2, kthe1, kthe2, rm, rm_1, rm_2, r, y, y2, d_u, d_l, sig_k, beta_k;
		double sinphi, tanthe, kpr, kpthe, kpphi, kphi;
		double time, c_0, c_TS, t_0, t_Celsius, c_Boussinesq;
		double sinthe, sinthe2, kro, c43, c13, k_Force;
		double costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2, rmtanthe;
		double resra, resre;
		double dudr, dudthe, dudphi, dvdr, dvdthe, dvdphi, dwdr, dwdthe, dwdphi, drdr, drdthe, drdphi;
		double d2udr2, d2udthe2, d2udphi2, d2vdr2, d2vdthe2, d2vdphi2, d2wdr2, d2wdthe2, d2wdphi2, d2rdr2, d2rdthe2, d2rdphi2;
		double d2udrdthe, d2udrdphi, d2vdrdthe, d2vdrdphi, d2wdrdthe, d2wdrdphi, d2rdrdthe, d2rdrdphi, d2udthedphi, d2vdthedphi, d2wdthedphi, d2rdthedphi;
		double d2udrdt, d2vdrdt, d2wdrdt, d2rdrdt, d2udthedt, d2vdthedt, d2wdthedt, d2rdthedt, d2udphidt, d2vdphidt, d2wdphidt, d2rdphidt;
		double dudt, dvdt, dwdt, drdt, dtdr, dpdr, dtdthe, dpdthe, dtdphi, dpdphi, d2tdr2, d2tdthe2, d2tdphi2;
		double dcdr, dcdthe, dcdphi, d2cdr2, d2cdthe2, d2cdphi2;
		double RS_Salt_Balance, RS_Salt_Energy, RS_Coriolis_Energy, RS_Centrifugal_Energy;
		double RS_Salt_Momentum, RS_Coriolis_Momentum_rad, RS_Coriolis_Momentum_the, RS_Coriolis_Momentum_phi, RS_Centrifugal_Momentum_rad, RS_Centrifugal_Momentum_the;
		double buoyancy, RS_buoyancy_Energy, RS_buoyancy_Momentum, RS_buoyancy_Water_Vapour; 
		double h_0_i, h_c_i, h_d_i, h_0_j, h_c_j, h_d_j, h_0_k, h_c_k, h_d_k, cc; 


	public:
		RHS_Hydrosphere ( int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double ); 
		~RHS_Hydrosphere ();

		void Pressure_RHS_Hydrosphere ( double, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );


		void RK_RHS_3D_Hydrosphere ( int, int, int, double, double, double, double, double, double, double, double, double, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void RK_RHS_2D_Hydrosphere ( int, int, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

};
#endif
