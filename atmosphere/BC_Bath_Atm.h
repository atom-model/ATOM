/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to read and prepare the bathymetric and topografic data
*/


#include <iostream>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"

#ifndef _BC_BATHYMETRIE_ATMOSPHERE_
#define _BC_BATHYMETRIE_ATMOSPHERE_

using namespace std;

class BC_Bathymetry_Atmosphere
{
    private:

        int im, jm, km;
        int j_half, j_max;
        int i_h, i_Ice_lauf;
        int NASATemperature;

        double d_i, d_i_max, d_j, d_j_half, d_j_max;
        double t_co2_eff;
        double Akkumulation_1, Akkumulation_2, Ablation, Ice_Balance_add_diff, min, max;
        double t_equi, t_equi_Celsius, t_plus, t_plus_Celsius, t_minus, t_minus_Celsius,
            t_minuss, t_minuss_Celsius, t_pluss, t_pluss_Celsius;
        double Hoehe_equi, Hoehe_equi_km, Hoehe_delta;
        double Ice_Hoehe;
        double co_co2_eff;
        double co2_vegetation, co2_ocean, co2_land;
        double h_point_max, h_land, h_ocean, ozean_land;
        
        BC_Bathymetry_Atmosphere(){};

    public:
        BC_Bathymetry_Atmosphere ( int, int, int, int, double, double, double );
        ~BC_Bathymetry_Atmosphere();

        void BC_MountainSurface( string &topo_filename, double L_atm, Array_1D &rad, Array_2D &Topography, Array &h );

        void BC_IceShield ( int, double, Array &, Array &, Array &, Array &, Array_2D &, Array_2D & );

        void BC_SolidGround ( int RadiationModel, int Ma, double g, double hp, double ep,
                                            double r_air, double R_Air, double t_0, double c_0, double t_land,
                                            double t_cretaceous, double t_equator, double t_pole,
                                            double t_tropopause, double c_land, double c_tropopause,
                                            double co2_0, double co2_equator, double co2_pole,
                                            double co2_tropopause, double pa, double gam, double sigma,
                                            Array &h, Array &u, Array &v, Array &w, Array &t, Array &p_dyn,
                                            Array &c, Array &cloud, Array &ice, Array &co2, Array &radiation_3D,
                                            Array_2D &Vegetation );

        void vegetationDistribution ( double, Array_2D &, Array_2D &, Array &, Array & );

        void land_oceanFraction ( Array & );
};
#endif
