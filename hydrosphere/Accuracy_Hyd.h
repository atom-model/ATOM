/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to surveil the accuracy of the iterations
*/


#include <iostream>

#include "Array.h"
#include "Array_1D.h"

#ifndef _ACCURACY_
#define _ACCURACY_

using namespace std;



class Accuracy
{
	private:
		int n, im, jm, km, l_anf, l_end, velocity_iter, pressure_iter, Ma;
		int i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_p, j_p, k_p;
		int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg, choice;

		double dr, dthe, dphi;
		double sinthe, costhe, rmsinthe;
		double dudr, dvdthe, dwdphi;
		double residuum, max_u, max_v, max_w, max_t, max_c, max_p;
		double Value, L_hyd;

		string name_Value;
		string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading;

	public:
		int i_res, j_res, k_res;
		double min, min_u, min_v, min_w, min_t, min_c, min_p;

		Accuracy ( int, int, int, int, double, double, double ); 
		Accuracy ( int, int, int, int, int, double, double ); 

		~Accuracy ();

		double residuumQuery_3D ( int &, int &, int &, double &, Array_1D &, Array_1D &, Array &, Array &, Array & );

		double steadyQuery_3D ( int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, double &, double &, double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void iterationPrintout_3D ( int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, double &, double &, double &, double &, double &, double & );

		double residuumQuery_2D ( int &, int &, double &, Array_1D &, Array_1D &, Array &, Array & );

		double steadyQuery_2D ( int &, int &, int &, int &, int &, int &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void iterationPrintout_2D ( int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &, double &, double &, double &, double & );

};
#endif

