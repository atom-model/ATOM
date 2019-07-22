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
    int im, jm, km, i_mount;

    double t_u, T, p_SL, p_h, hp, ep, R_Air, r_dry, g, L_atm, E_Rain, E_Ice, q_Rain, q_Ice, zeta;

    string output_path;

    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f);
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f);
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f);
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f);
    void dump_array(const string &name, Array &a, double multiplier, ofstream &f);

public:
    PostProcess_Atmosphere( int, int, int, string &output_path );
    ~PostProcess_Atmosphere();

    void paraview_panorama_vts (string &Name_Bathymetry_File, 
        int n, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, 
        double &co2_0, Array &h, Array &t, Array &p_dyn, Array &p_stat, Array &BuoyancyForce, 
        Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, 
        Array &aux_u, Array &aux_v, Array &aux_w, Array &Q_Latent, Array &Q_Sensible, 
        Array &epsilon_3D, Array &P_rain, Array &P_snow, Array &P_conv );

    void paraview_vtk_zonal ( string &Name_Bathymetry_File, 
        int k_zonal, int n, double &gam, double &hp, double &ep, double &R_Air, 
        double &g, double &L_atm, double &u_0, double &t_0, double &p_0, double &r_air, 
        double &c_0, double &co2_0, Array_1D &rad, Array &h, Array &p_dyn, 
        Array &p_stat, Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, 
        Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, 
        Array &aux_w, Array &Q_Latent, Array &Q_Sensible, Array &radiation_3D, 
        Array &epsilon_3D, Array &P_rain, Array &P_snow, Array &P_conv, Array &S_v, 
        Array &S_c, Array &S_i, Array &S_r, Array &S_s, Array &S_c_c, Array &M_u, 
        Array &M_d, Array &MC_s, Array &MC_q, Array &MC_v, Array &MC_w );

    void paraview_vtk_longal (string &Name_Bathymetry_File, 
        int j_longal, int n, double &L_atm, double &u_0, double &t_0, double &p_0, double &r_air, 
        double &c_0, double &co2_0, Array_1D &rad, Array &h, Array &p_dyn, Array &p_stat, 
        Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, Array &c, 
        Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, 
        Array &aux_w, Array &Q_Latent, Array &Q_Sensible, Array &epsilon_3D, 
        Array &P_rain, Array &P_snow, Array &P_conv, Array &M_u, Array &M_d, 
        Array &MC_s, Array &MC_q, Array &MC_v, Array &MC_w );

    void paraview_vtk_radial ( string &Name_Bathymetry_File, 
        int Ma, int i_radial, int n, double &u_0, double &t_0, double &p_0, double &r_air, 
        double &c_0, double &co2_0, Array &h, Array &p_dyn, Array &p_stat, 
        Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, Array &c, 
        Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, 
        Array &aux_w, Array &radiation_3D, Array &Q_Latent, Array &Q_Sensible, 
        Array &epsilon_3D, Array &P_rain, Array &P_snow, Array &P_conv, Array &M_u, Array &M_d, 
        Array &MC_s, Array &MC_q, Array &MC_v, Array &MC_w, Array_2D &precipitable_water, 
        Array_2D &Q_bottom, Array_2D &Q_radiation, Array_2D &Q_latent, Array_2D &Q_sensible, 
        Array_2D &Evaporation_Penman, Array_2D &Evaporation_Dalton, Array_2D &Q_Evaporation, 
        Array_2D &temperature_NASA, Array_2D &precipitation_NASA, Array_2D &Vegetation, 
        Array_2D &albedo, Array_2D &epsilon, Array_2D &Precipitation, 
        Array_2D &Topography, Array_2D &temp_NASA, Array_2D &vapour_evaporation, 
        Array_2D &radiation_surface );

    void Atmosphere_v_w_Transfer ( string &, double, Array &, Array &, Array &, Array &,
                                    Array_2D &, Array_2D & );

    void Atmosphere_PlotData ( string &Name_Bathymetry_File, int iter_cnt, double u_0, double t_0,
        Array &h, Array &v, Array &w, Array &t, Array &c, Array_2D &Precipitation, Array_2D &precipitable_water );

    double exp_func ( double &, const double &, const double & );

    void save( const string &filename, const std::vector<string> &field_names,
               const std::vector<Vector3D<>* > &data, unsigned layer=0 );
};
#endif
