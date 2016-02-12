/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to read and prepare the bathymetric and topografic data
*/


#include <iostream>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"

#ifndef _BC_BATHYMETRIE_ATMOSPHERE_
#define _BC_BATHYMETRIE_ATMOSPHERE_

using namespace std;



class BC_Bathymetry_Atmosphere
{
	private:

		int i, j, k, im, jm, km, l, ll, i_SL, i_land, i_max;
		int i_anf, i_end, j_anf, j_end, j_half, j_max, k_half, k_max, k_anf, k_end;
		int Ma, Ma_max, Ma_max_halb, i_Ice, i_Ice_lauf;
		int j_par, j_pol, k_par, k_pol;
		int j_r, k_r, j_sun;
		int *im_tropopause;

		double dummy_1, dummy_2, dummy_3, L_atm, c13, c43;
		double d_i, d_i_max, d_j, d_j_half, d_j_max, d_k_half, d_k_max, t_0, t_coeff, t_pole, t_tropopause, t_equator, t_12, p_0, c_tropopause, c_equator, c_coeff, c_pole, ep, hp, p_h, trop_coeff; 
		double t_delta, t_cretaceous_max, t_Cretaceous_coeff, t_Cretaceous, t_360, pi180;
		double Akkumulation_1, Akkumulation_2, Ablation, Ice_Balance_add_diff, min, max;
		double t_equi, t_equi_Celsius, t_plus, t_plus_Celsius, t_minus, t_minus_Celsius, t_minuss, t_minuss_Celsius, t_pluss, t_pluss_Celsius;
		double Hoehe_equi, Hoehe_equi_km, Hoehe_delta;
		double Ice_Hoehe, t_eff_earth, r;
		double co2_equator, co2_tropopause, co2_coeff, co2_pol, co2_vegetation, co2_ocean, co2_land, co2_Cretaceous;
		double j_par_f, j_pol_f, a, b, cc, d, e, j_d, t_d, k_par_f, k_pol_f, d_k;
		double *jm_temp_asym;

		string Daten_NW;


	public:

		BC_Bathymetry_Atmosphere ( int, int, int );
		~BC_Bathymetry_Atmosphere();

		void BC_MountainSurface ( const string &, double, Array &, Array & );

		void BC_IceShield ( int, double, Array &, Array &, Array &, Array &, Array_2D &, Array_2D & );

		void BC_SolidGround ( int, double, double, double, double, double, double, double, double, double, double, double, double, double, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D & );

};
#endif
