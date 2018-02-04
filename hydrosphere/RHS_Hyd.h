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
		int im, jm, km;
		int h_check_i, h_check_j, h_check_k;

		double dt, dr, dthe, dphi, r0;
		double re, pr, ec, sc, g;
		double dr2, dthe2, dphi2, rm2;
		double rm;
		double c_0, t_0, c_Boussinesq, r_fluid, r_salt_water, coeff_density, coeff_buoy;
		double sinthe, sinthe2, c43, c13, k_Force;
		double costhe, cotthe, rmsinthe, rm2sinthe, rm2sinthe2;
		double dudr, dudthe, dudphi, dvdr, dvdthe, dvdphi, dwdr, dwdthe, dwdphi;
		double d2udr2, d2udthe2, d2udphi2, d2vdr2, d2vdthe2, d2vdphi2, d2wdr2, d2wdthe2, d2wdphi2;
		double dtdr, dpdr, dtdthe, dpdthe, dtdphi, dpdphi, d2tdr2, d2tdthe2, d2tdphi2;
		double dcdr, dcdthe, dcdphi, d2cdr2, d2cdthe2, d2cdphi2;
		double RS_Salt_Balance, RS_Salt_Energy, RS_Coriolis_Energy, RS_Centrifugal_Energy;
		double RS_Coriolis_Momentum_rad, RS_Coriolis_Momentum_the, RS_Coriolis_Momentum_phi, RS_Centrifugal_Momentum_rad, RS_Centrifugal_Momentum_the;
		double buoyancy, RS_buoyancy_Energy, RS_buoyancy_Momentum;
		double h_0_i, h_c_i, h_d_i, h_0_j, h_c_j, h_d_j, h_0_k, h_c_k, h_d_k, cc;


	public:
		RHS_Hydrosphere ( int, int, double, double, double );

		RHS_Hydrosphere ( int, int, int, double, double, double, double, double, double, double, double, double, double, double );
		~RHS_Hydrosphere ();

		void RK_RHS_3D_Hydrosphere ( int, int, int, double, double, double, double, double, double, double, double, double, double, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void RK_RHS_2D_Hydrosphere ( int, int, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

};
#endif
