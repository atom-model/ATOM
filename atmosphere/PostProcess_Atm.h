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

class PostProcess_Atmosphere {
private:
    int im, jm, km;
    // int i_max, j_max, k_max;
    // double max_u, max_v, max_w, max_t, max_c, max_co2, max_cloud, max_ice, max_P_rain, max_P_snow, max_P_conv, max_M_u, max_M_d, max_p_dyn, max_p_stat, max_Rain, max_Rain_super, max_Ice, max_Latency, max_Q_Sensible, max_Precipitation, max_t_Evaporation, max_t_Condensation, max_t_evap_3D, max_t_cond_3D, max_precipitable_water, max_IceAir, max_Q_bottom, max_Q_latent, max_Q_sensible, max_t_Evaporation_Penman, max_t_Evaporation_Haude, max_Q_Radiation, max_Q_t_Evaporation, max_precipitation_NASA, max_Water, max_Water_super, max_Vegetation,  max_IceLayer, max_buoyancy_force;

	double t_u, T, p_SL, p_h, hp, ep, R_Air, r_dry, g, L_atm, E_Rain, E_Ice, q_Rain, q_Ice;
	double t_Cel, co_1, co_2, q_m, r_m_i;

    string output_path;

    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f);
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f);
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f);
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f);
    void dump_array(const string &name, Array &a, double multiplier, ofstream &f);

public:
    PostProcess_Atmosphere(int, int, int, string &output_path);
    ~PostProcess_Atmosphere();

    void paraview_vts ( string &, int &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );
    void paraview_panorama_vts ( string &, int &, double &, double &, double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );
    void paraview_vtk_zonal ( string &, int &, int &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &,Array &, Array &, Array & );
    void paraview_vtk_longal ( string &, int &, int &, double &, double &, double &, double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &,Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );
    void paraview_vtk_radial ( string &, int &, int &, double &, double &, double &,double &, double &, double &, double &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &,Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );
    void Atmosphere_v_w_Transfer ( string &, Array &, Array &, Array & );
    void Atmosphere_PlotData ( string &, double, double, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D & );

	double exp_func ( double &, const double &, const double & );
};
#endif
