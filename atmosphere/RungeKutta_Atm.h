/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include <iostream>
#include "Array.h"
#include "Array_1D.h"
#include "RHS_Atm.h"
#include "BC_Thermo.h"

#ifndef _RUNGEKUTTA_ATMOSPHERE_
#define _RUNGEKUTTA_ATMOSPHERE_

using namespace std;


class RungeKutta_Atmosphere
{
    private:
        int im, jm, km;

        double dt, dr, dphi, dthe, kt1, ku1, kv1, kw1, kp1, kc1, kcloud1, kice1, kco1,
                     kt2, ku2, kv2, kw2, kp2, kc2, kcloud2, kice2, kco2, kt3, ku3, kv3, kw3,
                     kp3, kc3, kcloud3, kice3, kco3, kt4, ku4, kv4, kw4, kp4, kc4, kcloud4, kice4, kco4;

    public:
        RungeKutta_Atmosphere ( int, int, int, double, double, double, double );
        ~RungeKutta_Atmosphere ();


        void solveRungeKutta_3D_Atmosphere ( RHS_Atmosphere &prepare,
            int &n, double lv, double ls, double ep, double hp, double u_0,
            double t_0, double c_0, double co2_0, double p_0, double r_air,
            double r_water_vapour, double r_co2, double L_atm, double cp_l,
            double R_Air, double R_WaterVapour, double R_co2,
            Array_1D &rad, Array_1D &the, Array_1D &phi, Array &rhs_t, Array &rhs_u,
            Array &rhs_v, Array &rhs_w, Array &rhs_c, Array &rhs_cloud, Array &rhs_ice,
            Array &rhs_co2, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn,
            Array &p_stat, Array &c, Array &cloud, Array &ice, Array &co2, Array &tn, Array &un,
            Array &vn, Array &wn, Array &p_dynn, Array &cn, Array &cloudn, Array &icen,
            Array &co2n, Array &aux_u, Array &aux_v, Array &aux_w, Array &Q_Latent,
            Array &BuoyancyForce, Array &Q_Sensible, Array &P_rain, Array &P_snow,
            Array &MC_s,Array &MC_q,Array &MC_v,Array &MC_w,
            Array &S_v, Array &S_c, Array &S_i, Array &S_r, Array &S_s, Array &S_c_c, Array &radiation_3D,
            Array_2D &Topography, Array_2D &Evaporation_Dalton,
            Array_2D &Precipitation, Array_2D &albedo );

        void solveRungeKutta_2D_Atmosphere ( RHS_Atmosphere &prepare_2D,
            int &n, double r_air, double u_0, double p_0, double L_atm, 
            Array_1D &rad, Array_1D &the, Array_1D &phi, Array &rhs_v, 
            Array &rhs_w, Array &h, Array &v, Array &w, Array &p_dyn, 
            Array &vn, Array &wn, Array &p_dynn, Array &aux_v, Array &aux_w );
};
#endif
