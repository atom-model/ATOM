/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include <iostream>
#include <cmath>
#include "RungeKutta_Hyd.h"
#include "RHS_Hyd.h"

using namespace std;


RungeKutta_Hydrosphere::RungeKutta_Hydrosphere ( int im, int jm, int km, double dt )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> dt = dt;
}

RungeKutta_Hydrosphere::~RungeKutta_Hydrosphere () {}



void RungeKutta_Hydrosphere::solveRungeKutta_3D_Hydrosphere ( RHS_Hydrosphere &prepare, int &n, double L_hyd, double g, double cp_w, double u_0, double t_0, double c_0, double r_0_water, double ta, double pa, double ca, Array_1D &rad, Array_1D &the, Array_1D &phi, Array_2D &Evaporation_Dalton, Array_2D &Precipitation, Array &h, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &p_dynn, Array &cn, Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion, Array &Buoyancy_Force, Array &Salt_Balance, Array &p_stat, Array &r_water, Array &r_salt_water, Array_2D &Bathymetry )
{
//  3D volume iterations
// Runge-Kutta 4. order for u, v and w component, temperature and salt concentration

    for ( int i = 1; i < im-1; i++ )
    {
        for ( int j = 1; j < jm-1; j++ )
        {
            for ( int k = 1; k < km-1; k++ )
            {
                tn.x[ i ][ j ][ k ] = t.x[ i ][ j ][ k ];
                un.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ];
                vn.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ];
                wn.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ];
                cn.x[ i ][ j ][ k ] =  c.x[ i ][ j ][ k ];


// Runge-Kutta 4. order for k1 step ( dt )
                prepare.RK_RHS_3D_Hydrosphere ( i, j, k, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Buoyancy_Force, Salt_Balance, p_stat, r_water, r_salt_water, Evaporation_Dalton, Precipitation, Bathymetry );

                kt1 = rhs_t.x[ i ][ j ][ k ];
                ku1 = rhs_u.x[ i ][ j ][ k ];
                kv1 = rhs_v.x[ i ][ j ][ k ];
                kw1 = rhs_w.x[ i ][ j ][ k ];
                kc1 = rhs_c.x[ i ][ j ][ k ];

                t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt1 * .5 * dt;
                u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku1 * .5 * dt;
                v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv1 * .5 * dt;
                w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw1 * .5 * dt;
                c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc1 * .5 * dt;

// Runge-Kutta 4. order for k2 step ( dt )
                prepare.RK_RHS_3D_Hydrosphere ( i, j, k, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Buoyancy_Force, Salt_Balance, p_stat, r_water, r_salt_water, Evaporation_Dalton, Precipitation, Bathymetry );

                kt2 = rhs_t.x[ i ][ j ][ k ];
                ku2 = rhs_u.x[ i ][ j ][ k ];
                kv2 = rhs_v.x[ i ][ j ][ k ];
                kw2 = rhs_w.x[ i ][ j ][ k ];
                kc2 = rhs_c.x[ i ][ j ][ k ];

                t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt2 * .5 * dt;
                u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku2 * .5 * dt;
                v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv2 * .5 * dt;
                w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw2 * .5 * dt;
                c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc2 * .5 * dt;

    // Runge-Kutta 4. order for k3 step ( dt )
                prepare.RK_RHS_3D_Hydrosphere ( i, j, k, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Buoyancy_Force, Salt_Balance, p_stat, r_water, r_salt_water, Evaporation_Dalton, Precipitation, Bathymetry );

                kt3 = rhs_t.x[ i ][ j ][ k ];
                ku3 = rhs_u.x[ i ][ j ][ k ];
                kv3 = rhs_v.x[ i ][ j ][ k ];
                kw3 =rhs_w.x[ i ][ j ][ k ];
                kc3 = rhs_c.x[ i ][ j ][ k ];

                t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt3 * dt;
                u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku3 * dt;
                v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv3 * dt;
                w.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] + kw3 * dt;
                c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc3 * dt;

    // Runge-Kutta 4. order for k4 step ( dt )
                prepare.RK_RHS_3D_Hydrosphere ( i, j, k, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Buoyancy_Force, Salt_Balance, p_stat, r_water, r_salt_water, Evaporation_Dalton, Precipitation, Bathymetry );

                kt4 = rhs_t.x[ i ][ j ][ k ];
                ku4 = rhs_u.x[ i ][ j ][ k ];
                kv4 = rhs_v.x[ i ][ j ][ k ];
                kw4 = rhs_w.x[ i ][ j ][ k ];
                kc4 = rhs_c.x[ i ][ j ][ k ];

                t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + dt * ( kt1 + 2. * kt2 + 2. * kt3 + kt4 ) / 6.;
                u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + dt * ( ku1 + 2. * ku2 + 2. * ku3 + ku4 ) / 6.;
                v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + dt * ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
                w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + dt * ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
                c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + dt * ( kc1 + 2. * kc2 + 2. * kc3 + kc4 ) / 6.;

                tn.x[ i ][ j ][ k ] = t.x[ i ][ j ][ k ];
                un.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ];
                vn.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ];
                wn.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ];
                cn.x[ i ][ j ][ k ] =  c.x[ i ][ j ][ k ];
            }
        }
    }
}





void RungeKutta_Hydrosphere::solveRungeKutta_2D_Hydrosphere ( RHS_Hydrosphere &prepare_2D, int &n, double r_0_water, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &rhs_v, Array &rhs_w, Array &h, Array &v, Array &w, Array &p_dyn, Array &vn, Array &wn, Array &p_dynn, Array &aux_v, Array &aux_w )
{
//  2D surface iterations
// Runge-Kutta 4. order for u, v and w component, temperature and salt concentration

    cout.precision ( 9 );
    cout.setf ( ios::fixed );

    for ( int j = 1; j < jm-1; j++ )
    {
        for ( int k = 1; k < km-1; k++ )
        {

            vn.x[ im-1 ][ j ][ k ] = v.x[ im-1 ][ j ][ k ];
            wn.x[ im-1 ][ j ][ k ] = w.x[ im-1 ][ j ][ k ];

// Runge-Kutta 4. order for k1 step ( dt )
            prepare_2D.RK_RHS_2D_Hydrosphere ( j, k, r_0_water, rad, the, phi, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, aux_v, aux_w );

            kv1 = rhs_v.x[ im-1 ][ j ][ k ];
            kw1 = rhs_w.x[ im-1 ][ j ][ k ];

            v.x[ im-1 ][ j ][ k ] = vn.x[ im-1 ][ j ][ k ] + kv1 * .5 * dt;
            w.x[ im-1 ][ j ][ k ] = wn.x[ im-1 ][ j ][ k ] + kw1 * .5 * dt;

    // Runge-Kutta 4. order for k2 step ( dt )
            prepare_2D.RK_RHS_2D_Hydrosphere ( j, k, r_0_water, rad, the, phi, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, aux_v, aux_w );

            kv2 = rhs_v.x[ im-1 ][ j ][ k ];
            kw2 = rhs_w.x[ im-1 ][ j ][ k ];

            v.x[ im-1 ][ j ][ k ] = vn.x[ im-1 ][ j ][ k ] + kv2 * .5 * dt;
            w.x[ im-1 ][ j ][ k ] = wn.x[ im-1 ][ j ][ k ] + kw2 * .5 * dt;

        // Runge-Kutta 4. order for k3 step ( dt )
            prepare_2D.RK_RHS_2D_Hydrosphere ( j, k, r_0_water, rad, the, phi, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, aux_v, aux_w );

            kv3 = rhs_v.x[ im-1 ][ j ][ k ];
            kw3 = rhs_w.x[ im-1 ][ j ][ k ];

            v.x[ im-1 ][ j ][ k ] = vn.x[ im-1 ][ j ][ k ] + kv3 * dt;
            w.x[ im-1 ][ j ][ k ] = wn.x[ im-1 ][ j ][ k ] + kw3 * dt;

        // Runge-Kutta 4. order for k4 step ( dt )
            prepare_2D.RK_RHS_2D_Hydrosphere ( j, k, r_0_water, rad, the, phi, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, aux_v, aux_w );

            kv4 = rhs_v.x[ im-1 ][ j ][ k ];
            kw4 = rhs_w.x[ im-1 ][ j ][ k ];

            v.x[ im-1 ][ j ][ k ] = vn.x[ im-1 ][ j ][ k ] + dt * ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
            w.x[ im-1 ][ j ][ k ] = wn.x[ im-1 ][ j ][ k ] + dt * ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;

            vn.x[ im-1 ][ j ][ k ] = v.x[ im-1 ][ j ][ k ];
            wn.x[ im-1 ][ j ][ k ] = w.x[ im-1 ][ j ][ k ];
        }
    }
}
