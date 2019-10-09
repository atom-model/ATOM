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
#include "cHydrosphereModel.h"

using namespace std;

void cHydrosphereModel::solveRungeKutta_3D_Hydrosphere()
{
    float kt1, ku1, kv1, kw1, kc1;
    float kt2, ku2, kv2, kw2, kc2;
    float kt3, ku3, kv3, kw3, kc3;
    float kt4, ku4, kv4, kw4, kc4;
    //  3D volume iterations
    // Runge-Kutta 4. order for u, v and w component, temperature and salt concentration
    for ( int i = 1; i < im-1; i++ )
    {
        for ( int j = 1; j < jm-1; j++ )
        {
            for ( int k = 1; k < km-1; k++ )
            {
                // Runge-Kutta 4. order for k1 step ( dt )
                RK_RHS_3D_Hydrosphere(i,j,k);

                kt1 = rhs_t.x[i][j][k];
                ku1 = rhs_u.x[i][j][k];
                kv1 = rhs_v.x[i][j][k];
                kw1 = rhs_w.x[i][j][k];
                kc1 = rhs_c.x[i][j][k];

                t.x[i][j][k] = tn.x[i][j][k] + kt1 * .5 * dt;
                u.x[i][j][k] = un.x[i][j][k] + ku1 * .5 * dt;
                v.x[i][j][k] = vn.x[i][j][k] + kv1 * .5 * dt;
                w.x[i][j][k] = wn.x[i][j][k] + kw1 * .5 * dt;
                c.x[i][j][k] = cn.x[i][j][k] + kc1 * .5 * dt;

                // Runge-Kutta 4. order for k2 step ( dt )
                RK_RHS_3D_Hydrosphere( i, j, k);

                kt2 = rhs_t.x[i][j][k];
                ku2 = rhs_u.x[i][j][k];
                kv2 = rhs_v.x[i][j][k];
                kw2 = rhs_w.x[i][j][k];
                kc2 = rhs_c.x[i][j][k];

                t.x[i][j][k] = tn.x[i][j][k] + kt2 * .5 * dt;
                u.x[i][j][k] = un.x[i][j][k] + ku2 * .5 * dt;
                v.x[i][j][k] = vn.x[i][j][k] + kv2 * .5 * dt;
                w.x[i][j][k] = wn.x[i][j][k] + kw2 * .5 * dt;
                c.x[i][j][k] = cn.x[i][j][k] + kc2 * .5 * dt;

                // Runge-Kutta 4. order for k3 step ( dt )
                RK_RHS_3D_Hydrosphere(i,j,k);

                kt3 = rhs_t.x[i][j][k];
                ku3 = rhs_u.x[i][j][k];
                kv3 = rhs_v.x[i][j][k];
                kw3 = rhs_w.x[i][j][k];
                kc3 = rhs_c.x[i][j][k];

                t.x[i][j][k] = tn.x[i][j][k] + kt3 * dt;
                u.x[i][j][k] = un.x[i][j][k] + ku3 * dt;
                v.x[i][j][k] = vn.x[i][j][k] + kv3 * dt;
                w.x[i][j][k] = wn.x[i][j][k] + kw3 * dt;
                c.x[i][j][k] = cn.x[i][j][k] + kc3 * dt;

                // Runge-Kutta 4. order for k4 step ( dt )
                RK_RHS_3D_Hydrosphere( i, j, k);

                kt4 = rhs_t.x[i][j][k];
                ku4 = rhs_u.x[i][j][k];
                kv4 = rhs_v.x[i][j][k];
                kw4 = rhs_w.x[i][j][k];
                kc4 = rhs_c.x[i][j][k];

                t.x[i][j][k] = tn.x[i][j][k] + dt * ( kt1 + 2. * kt2 + 2. * kt3 + kt4 ) / 6.;
                u.x[i][j][k] = un.x[i][j][k] + dt * ( ku1 + 2. * ku2 + 2. * ku3 + ku4 ) / 6.;
                v.x[i][j][k] = vn.x[i][j][k] + dt * ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
                w.x[i][j][k] = wn.x[i][j][k] + dt * ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
                c.x[i][j][k] = cn.x[i][j][k] + dt * ( kc1 + 2. * kc2 + 2. * kc3 + kc4 ) / 6.;
            }
        }
    }
}


void cHydrosphereModel::solveRungeKutta_2D_Hydrosphere()
{
    float kv1, kw1;
    float kv2, kw2;
    float kv3, kw3;
    float kv4, kw4;
    //  2D surface iterations
    // Runge-Kutta 4. order for u, v and w component, temperature and salt concentration
    for ( int j = 1; j < jm-1; j++ )
    {
        for ( int k = 1; k < km-1; k++ )
        {
            // Runge-Kutta 4. order for k1 step ( dt )
            RK_RHS_2D_Hydrosphere( j, k);

            kv1 = rhs_v.x[im-1][j][k];
            kw1 = rhs_w.x[im-1][j][k];

            v.x[im-1][j][k] = vn.x[im-1][j][k] + kv1 * .5 * dt;
            w.x[im-1][j][k] = wn.x[im-1][j][k] + kw1 * .5 * dt;

            // Runge-Kutta 4. order for k2 step ( dt )
            RK_RHS_2D_Hydrosphere( j, k);

            kv2 = rhs_v.x[im-1][j][k];
            kw2 = rhs_w.x[im-1][j][k];

            v.x[im-1][j][k] = vn.x[im-1][j][k] + kv2 * .5 * dt;
            w.x[im-1][j][k] = wn.x[im-1][j][k] + kw2 * .5 * dt;

            // Runge-Kutta 4. order for k3 step ( dt )
            RK_RHS_2D_Hydrosphere( j, k);

            kv3 = rhs_v.x[im-1][j][k];
            kw3 = rhs_w.x[im-1][j][k];

            v.x[im-1][j][k] = vn.x[im-1][j][k] + kv3 * dt;
            w.x[im-1][j][k] = wn.x[im-1][j][k] + kw3 * dt;

            // Runge-Kutta 4. order for k4 step ( dt )
            RK_RHS_2D_Hydrosphere( j, k);

            kv4 = rhs_v.x[im-1][j][k];
            kw4 = rhs_w.x[im-1][j][k];

            v.x[im-1][j][k] = vn.x[im-1][j][k] + dt * ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
            w.x[im-1][j][k] = wn.x[im-1][j][k] + dt * ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
        }
    }
}
