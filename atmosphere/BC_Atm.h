/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the coordinate system
*/

#include <iostream>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"

#ifndef _BC_ATMOSPHAERE_
#define _BC_ATMOSPHERE_

using namespace std;



class BC_Atmosphere
{
	private:
		int im, jm, km;

		double c13, c43;
		double dr, dr2, dthe, dthe2, dphi, dphi2, rm, rm2, sinthe, sinthe2, costhe, rmsinthe, rm2sinthe, rm2sinthe2, cotthe;
		double LHS_u, RHS_u, LHS_v, RHS_v, LHS_w, RHS_w, mue;
		double dudr, dudthe, dudphi, dvdr, dvdthe, dvdphi, dwdr, dwdthe, dwdphi, drdr, drdthe, drdphi;
		double d2udr2, d2udthe2, d2udphi2, d2vdr2, d2vdthe2, d2vdphi2, d2wdr2, d2wdthe2, d2wdphi2, d2rdr2, d2rdthe2, d2rdphi2;
		double d2udrdthe, d2udrdphi, d2vdrdthe, d2vdrdphi, d2wdrdthe, d2wdrdphi, d2rdrdthe, d2rdrdphi, d2udthedphi, d2vdthedphi, d2wdthedphi, d2rdthedphi;
		double d2udrdt, d2vdrdt, d2wdrdt, d2rdrdt, d2udthedt, d2vdthedt, d2wdthedt, d2rdthedt, d2udphidt, d2vdphidt, d2wdphidt, d2rdphidt;
		double dudt, dvdt, dwdt, drdt, dtdr, dpdr, dpdr_atm, dtdthe, dpdthe, dtdphi, dpdphi, d2tdr2, d2tdthe2, d2tdphi2;
		double dcdr, dcdthe, dcdphi, d2cdr2, d2cdthe2, d2cdphi2;
		double dco2dr, dco2dthe, dco2dphi, d2co2dr2, d2co2dthe2, d2co2dphi2;

	public:
		BC_Atmosphere ( int, int, int );
		~BC_Atmosphere();


		void BC_radius ( double, double, double, double, double, double, Array_1D &, double, double, double, Array_2D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_theta ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_phi ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void BC_NST_control_2D ( double, double, double, double, double, double, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_1D &, Array_1D & );

		void BC_NST_control_3D ( double, double, double, double, double, double, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_1D &, Array_1D & );
};
#endif
