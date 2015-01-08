/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to produce resulting data on mean sea level
*/


//#include <iostream>

#include "Array.h"
#include "Array_2D.h"

#ifndef Results_MSL_Atm_Atm
#define Results_MSL_Atm_Atm

using namespace std;

class Results_MSL_Atm
{
	private:
		int i, j, k, im, jm, km, l, sun, h_point_max, h_land, h_ocean;
		int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;

		double t_Celsius_SL, e_SL, a_SL, c_grad, t_grad, t_denom, Delta, E_Rain_SL, E_Rain, sat_deficit, E_a, gamma;
		double ep, hp, p_0, t_0, c_0, sigma, albedo, lv, cp_l, r_0_air, dr, ozean_land, p_baro, L_atm;
		double p_h, e_h, a_h, R_Air, r_h, r_0_water_vapour, R_WaterVapour;
		double coeff_mmWS, coeff_Diffusion_latent, coeff_Diffusion_sensibel, f_Haude, co2_vegetation, co2_ocean, co2_land;
		double Value_1, Value_2, Value_3, Value_4, Value_5, Value_6;

		Array				hc ( int, int, int, double ), Rai ( int, int, int, double ), Rai_su ( int, int, int, double ), Ic ( int, int, int, double ), Lat ( int, int, int, double ), cond_3D ( int, int, int, double ), evap_3D ( int, int, int, double ), cc ( int, int, int, double ), tc ( int, int, int, double ), pc ( int, int, int, double ), uc ( int, int, int, double );

		string name_Value_1, name_Value_2, name_Value_3, name_Value_4, name_Value_5, name_Value_6, name_unit_wm2, name_unit_mm;
		string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading;

	public:
		Results_MSL_Atm ( int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double ); 

		Array_2D		prec ( int, int, double ), ice ( int, int, double ), rai ( int, int, double ), rai_su ( int, int, double ), evap ( int, int, double ), cond ( int, int, double ), veg ( int, int, double ), prec_wat ( int, int, double ), Q_Bal ( int, int, double ), Q_evap ( int, int, double ), Q_lat ( int, int, double ), Q_sen ( int, int, double ), Q_dif ( int, int, double ), evap_Pen ( int, int, double ), evap_Hau ( int, int, double );

		~Results_MSL_Atm (  );

		void run_MSL_data ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void show_MSL_data ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void land_oceanFraction ( Array & );

		void vegetationDistribution ( double, Array_2D &, Array_2D &, Array &, Array & );
};
#endif
