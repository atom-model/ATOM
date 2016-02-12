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
#include "Array_2D.h"

#ifndef Results_MSL_Atm_Atm
#define Results_MSL_Atm_Atm

using namespace std;

class Results_MSL_Atm
{
	private:
		int i, j, k, l, im, jm, km, sun, h_point_max, h_land, h_ocean;
		int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;

		double E, E_Rain_SL, E_Rain, E_Rain_super, E_Ice, q_Rain, q_Rain_super, q_Ice, E_Rain_super_SL, E_Ice_SL, q_Rain_SL, q_Rain_super_SL, q_Ice_SL;
		double e, a, e_SL, a_SL, p_SL, q_SL, t_dew_SL, t_Celsius_SL;
		double e_h, a_h, p_h, q_h, t_dew, t_Celsius, c_grad, t_grad, t_denom, Delta, E_a, gamma, g, gam;
		double i_level, h_level, h_h, sat_deficit, RF_e, RF_c, Evap_Haude, Evaporation_Penman_average, Evaporation_Haude_average;

		double ep, hp, u_0, p_0, t_0, c_0, sigma, albedo_extra, lv, cp_l, r_0_air, dr, ozean_land, p_baro, L_atm, c13, c43;
		double R_Air, r_h, r_0_water_vapour, R_WaterVapour, precipitablewater_average, precipitation_average, precipitation_NASA_average;
		double coeff_mmWS, coeff_lv, coeff_Diffusion_latent, coeff_Diffusion_sensibel, f_Haude, co2_vegetation, co2_ocean, co2_land;
		double Value_1, Value_2, Value_3, Value_4, Value_5, Value_6, Value_7, Value_8, Value_10, Value_12, Value_13, Value_14, Value_15, Value_16, Value_17, Value_18, Value_19;

		string name_Value_1, name_Value_2, name_Value_3, name_Value_4, name_Value_5, name_Value_6, name_Value_7, name_Value_8, name_Value_9, name_Value_10, name_Value_11, name_Value_12, name_Value_13, name_Value_14, name_Value_15, name_Value_16, name_Value_17, name_Value_18, name_Value_19, name_Value_20, name_Value_21, name_unit_wm2, name_unit_mm, name_unit_mmd, name_unit_mma;
		string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading, heading_Dresden, heading_Sydney, heading_Equator;

	public:
		Results_MSL_Atm ( int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double ); 


		~Results_MSL_Atm (  );

		void run_MSL_data ( int, double, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void show_MSL_data ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void land_oceanFraction ( Array & );

		void vegetationDistribution ( double, Array_2D &, Array_2D &, Array &, Array & );
};
#endif
