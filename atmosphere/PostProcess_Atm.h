/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/


#include <iostream>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"

#ifndef _POSTPROCESS_ATMOSPHERE_
#define _POSTPROCESS_ATMOSPHERE_

using namespace std;


class PostProcess_Atmosphere
{
	private:
		int im, jm, km, jr, kr, a, i_radial;
		int i_max, j_max, k_max, l;

		double rotu, rotv, rotw;

		double max_u, max_v, max_w, max_t, max_c, max_co2, max_cloud, max_ice, max_P_rain, max_P_snow, max_P_conv, max_M_u, max_M_d, max_p_dyn, max_p_stat, max_Rain, max_Rain_super, max_Ice, max_Latency, max_Q_Sensible, max_Precipitation, max_albedo, max_epsilon, max_t_Evaporation, max_t_Condensation, max_t_evap_3D, max_t_cond_3D, max_precipitable_water, max_IceAir, max_Q_bottom, max_Q_latent, max_Q_sensible, max_t_Evaporation_Penman, max_t_Evaporation_Haude, max_Radiation_Balance, max_Q_Radiation, max_Q_t_Evaporation, max_precipitation_NASA, max_Water, max_Water_super, max_Vegetation,  max_IceLayer, max_buoyancy_force, max_radiation_3D;


	public:
		PostProcess_Atmosphere ( int, int, int );
		~PostProcess_Atmosphere();


		void paraview_vts ( const string &, int &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_panorama_vts ( const string &, int &, double &, double &, double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_vtk_zonal ( const string &, int &, int &, double &, double &, double &, double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_vtk_longal ( const string &, int &, int &, double &, double &, double &, double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_vtk_radial ( const string &, int &, int &, double &, double &, double &,double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void Atmosphere_SequelFile_write ( const string &, int &, double &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void Atmosphere_SequelFile_read ( const string &, int &, double &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void Atmosphere_v_w_Transfer ( const string &, Array &, Array &, Array & );

		void Atmosphere_PlotData ( const string &, double, double, Array &, Array &, Array &, Array &, Array_2D &, Array_2D & );
};
#endif
