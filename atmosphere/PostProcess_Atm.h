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

		double max_u, max_v, max_w, max_t, max_c, max_co2, max_p, max_Rain, max_Rain_super, max_Ice, max_Latency, max_Precipitation, max_Evaporation, max_Condensation, max_Evaporation_3D, max_Condensation_3D, max_precipitable_water, max_IceAir, max_Q_bottom, max_Q_latent, max_Q_sensible, max_Evaporation_Penman, max_Evaporation_Haude, max_Radiation_Balance, max_Q_Evaporation, max_precipitation_j, max_Water, max_Water_super, max_Vegetation,  max_IceLayer;


	public:
		PostProcess_Atmosphere ( int, int, int );
		~PostProcess_Atmosphere();


		void paraview_vts ( const string &, int &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_panorama_vts ( const string &, int &, const double &, const double &, const double &, const double &, const double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_vtk_zonal ( const string &, int &, int &, const double &, const double &, const double &,const double &,const double &, const double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_vtk_longal ( const string &, int &, int &, const double &, const double &, const double &, const double &, const double &, const double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_vtk_radial ( const string &, int &, int &, const double &, const double &,const double &, const double &, const double &, const double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void Atmosphere_SequelFile_write ( const string &, int &, double &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D & );

		void Atmosphere_SequelFile_read ( const string &, int &, double &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D & );

		void Atmosphere_v_w_Transfer ( const string &, Array &, Array &, Array & );

		void Atmosphere_PlotData ( const string &, double, double, Array &, Array &, Array &, Array &, Array_2D &, Array_2D & );
};
#endif
