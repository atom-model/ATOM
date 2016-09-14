/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
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

#ifndef _Results_Hyd_
#define _Results_Hyd_

using namespace std;


class Results_Hyd
{
	private:
		int im, jm, km, sun, h_land, h_ocean, h_ground, h_point_max;
		int j_max_Salt_total, k_max_Salt_total, j_min_Salt_total, k_min_Salt_total;
		int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;

		double co2_vegetation, co2_ocean, co2_land, ozean_land;
		double h_max_salt, h_min_salt, max_Salt_total, min_Salt_total;
		double Value_1, Value_2, Value_3, Value_4, Value_5, Value_6;
		double c43, c13;

		string name_Value_1, name_Value_2, name_Value_3, name_Value_4, name_Value_5, name_Value_6, name_unit_ms, name_unit_psu;
		string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading;

	public:
		Results_Hyd ( int, int, int ); 

		~Results_Hyd (  );

		void run_data ( double, double, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void land_oceanFraction ( Array & );

//		void show_data ( double, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

};
#endif

