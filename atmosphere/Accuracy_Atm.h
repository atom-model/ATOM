/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to surveil the accuracy of the iterations
*/



#include "Array.h"
#include "Array_1D.h"

#ifndef _ACCURACY_
#define _ACCURACY_

using namespace std;



class Accuracy_Atm
{
	private:
		int n, nm, im, jm, km;
		int velocity_iter_2D, pressure_iter_2D, velocity_iter, pressure_iter, velocity_iter_max, pressure_iter_max, velocity_iter_max_2D, pressure_iter_max_2D, Ma;
		int i_res, j_res, k_res;
		int i_u, j_u, k_u, i_v, j_v, k_v, i_w, j_w, k_w, i_t, j_t, k_t, i_c, j_c, k_c, i_cloud, j_cloud, k_cloud, i_ice, j_ice, k_ice, i_co2, j_co2, k_co2, i_p, j_p, k_p;
		int i_loc, j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;

		double dr, dthe, dphi;
		double sinthe, costhe, rmsinthe;
		double dudr, dvdthe, dwdphi;
		double residuum, max_u, max_v, max_w, max_t, max_c, max_cloud, max_ice, max_p, max_co2;
		double min, min_u, min_v, min_w, min_t, min_c, min_cloud, min_ice, min_co2, min_p;
		double Value, L_atm;

		string name_Value;
		string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading;

	public:
		Accuracy_Atm ( int, int, int, double, double );
		Accuracy_Atm ( int, int, int, double, double, double );
        Accuracy_Atm ( int n, int nm, int Ma, int im, int jm, int km, double min, int j_res, int k_res );
        Accuracy_Atm ( int n, int nm, int Ma, int im, int jm, int km, double min, int i_res, int j_res, int k_res, double L_atm );

		~Accuracy_Atm ();


		double residuumQuery_2D ( Array_1D &, Array_1D &, Array &, Array & );
		double residuumQuery_3D ( Array_1D &, Array_1D &, Array &, Array &, Array & );

		double steadyQuery_2D ( Array &, Array &, Array &, Array &, Array &, Array & );
		double steadyQuery_3D ( Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

        void print(const string& name, double value, int j, int k) const;        
        void print(const string& name, double value, int i, int j, int k) const;

		double out_min (  ) const;
		int out_i_res (  ) const;
		int out_j_res (  ) const;
		int out_k_res (  ) const;

};
#endif
