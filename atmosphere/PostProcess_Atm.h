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

    double t_u, T, p_SL, p_h, hp, ep, R_Air, r_dry, g, L_atm, E_Rain, E_Ice, q_Rain, q_Ice;

    string output_path;

    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f);
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f);
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f);
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f);
    void dump_array(const string &name, Array &a, double multiplier, ofstream &f);

public:
    PostProcess_Atmosphere(int, int, int, string &output_path);
    ~PostProcess_Atmosphere();

    void paraview_vts ( string &, int iter_cnt, Array_1D &, Array_1D &, Array_1D &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array & );

    void paraview_panorama_vts ( string &, int iter_cnt, double &, double &, double &,
                                    double &, double &, double &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

    void paraview_vtk_zonal ( string &, int, int iter_cnt, double &, double &, double &,
                                    double &, double &, double &, double &, double &, double &,
                                    double &, double &, Array &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &,Array & );

    void paraview_vtk_longal ( string &, int, int iter_cnt, double &, double &, double &,
                                    double &, double &, double &, Array &, Array &, Array &, Array &,
                                    Array &,Array &, Array &, Array &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array & );

    void paraview_vtk_radial ( string &, int, int, int iter_cnt, double &, double &, double &,
                                    double &, double &, double &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &,
                                    Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &,
                                    Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &,
                                    Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &,
                                    Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D & );

    void Atmosphere_v_w_Transfer ( string &, double, Array &, Array &, Array &, Array &,
                                    Array_2D &, Array_2D & );

    void Atmosphere_PlotData ( string &Name_Bathymetry_File, int iter_cnt, double u_0, double t_0,
        Array &h, Array &v, Array &w, Array &t, Array &c, Array_2D &Precipitation, Array_2D &precipitable_water );

    double exp_func ( double &, const double &, const double & );

    void save( const string &filename, const std::vector<string> &field_names,
               const std::vector<Vector3D<>* > &data, unsigned layer=0 );
};
#endif
