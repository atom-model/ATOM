/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include <iostream>
#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "BC_Thermo.h"

#ifndef _RHS_ATMOSPHERE_
#define _RHS_ATMOSPHERE_

class cAtmosphereModel;

using namespace std;


class RHS_Atmosphere
{
    private:
        cAtmosphereModel* m_model;
        int im, jm, km;

        double dt, dr, dthe, dphi, zeta;
        double re, sc_WaterVapour, sc_CO2, g, pr, gam, WaterVapour, Buoyancy, CO2, sigma;

    public:
        RHS_Atmosphere ( int, int, double, double, double );
        RHS_Atmosphere ( cAtmosphereModel* model, int, int, int, double, double, double, double, double, double,
                                       double, double, double, double, double, double, double, double, double );
        ~RHS_Atmosphere ();

        void RK_RHS_3D_Atmosphere ( int n, int i, int j, int k, double lv, double ls, double ep,
                                            double hp, double u_0, double t_0, double c_0, double co2_0,
                                            double p_0, double r_air, double r_water_vapour, double r_co2,
                                            double L_atm, double cp_l, double R_Air, double R_WaterVapour,
                                            double R_co2, Array_1D &rad, Array_1D &the, Array_1D &phi,
                                            Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn,
                                            Array &p_stat, Array &c, Array &cloud, Array &ice, Array &co2,
                                            Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c,
                                            Array &rhs_cloud, Array &rhs_ice, Array &rhs_co2, Array &aux_u,
                                            Array &aux_v, Array &aux_w, Array &Q_Latent, Array &BuoyancyForce,
                                            Array &Q_Sensible, Array &P_rain, Array &P_snow, Array &S_v,
                                            Array &S_c, Array &S_i, Array &S_r, Array &S_s, Array &S_c_c,
                                            Array_2D &Topography, Array_2D &Evaporation_Dalton,
                                            Array_2D &Precipitation );


        void RK_RHS_2D_Atmosphere ( int j, int k, double r_air, double u_0, double p_0, double L_atm,
                                            Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &v, Array &w,
                                            Array &p_dyn, Array &rhs_v, Array &rhs_w, Array &aux_v, Array &aux_w );
};
#endif
