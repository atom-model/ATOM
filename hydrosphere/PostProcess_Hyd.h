/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
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

#ifndef _POSTPROCESS_HYDROSPHERE_
#define _POSTPROCESS_HYDROSPHERE_

using namespace std;


class PostProcess_Hydrosphere
{
	private:
		int im, jm, km;
		int i_max, j_max, k_max;

		double max_p_dyn;

		string input_path, output_path;


	public:
		PostProcess_Hydrosphere (int, int, int, const string &input_path, const string &output_path);
		~PostProcess_Hydrosphere();


		void Hydrosphere_SequelFile_write ( const string &, int &, int &, double &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void Hydrosphere_SequelFile_read ( const string &, int &, int &, double &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void Atmosphere_TransferFile_read ( const string &, Array &, Array &, Array & );

		void paraview_vts ( const string &, int &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_panorama_vts ( const string &, int &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_vtk_radial ( const string &, int &, int &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

		void paraview_vtk_longal ( const string &, int &, int &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void paraview_vtk_zonal ( const string &, int &, int &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

		void Hydrosphere_PlotData ( const string &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D & );

};
#endif
