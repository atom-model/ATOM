/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include <iostream>
#include <cmath>
#include "cAtmosphereModel.h"
#include "RungeKutta_Atm.h"
#include "RHS_Atm.h"

using namespace std;

void cAtmosphereModel::solveRungeKutta_3D_Atmosphere(){
    double kt1, ku1, kv1, kw1, kc1, kcloud1, kice1, kco1,
           kt2, ku2, kv2, kw2, kc2, kcloud2, kice2, kco2, 
           kt3, ku3, kv3, kw3, kc3, kcloud3, kice3, kco3, 
           kt4, ku4, kv4, kw4, kc4, kcloud4, kice4, kco4;
    // Runge-Kutta 4. order for u, v and w component, temperature, water vapour and co2 content
    for ( int i = 1; i < im-1; i++ ){
        for ( int j = 1; j < jm-1; j++ ){
            for ( int k = 1; k < km-1; k++ ){
                // Runge-Kutta 4. order for k1 step ( dt )
                RK_RHS_3D_Atmosphere(i, j, k);

                kt1 = rhs_t.x[ i ][ j ][ k ];
                ku1 = rhs_u.x[ i ][ j ][ k ];
                kv1 = rhs_v.x[ i ][ j ][ k ];
                kw1 = rhs_w.x[ i ][ j ][ k ];
                kc1 = rhs_c.x[ i ][ j ][ k ];
                kcloud1 = rhs_cloud.x[ i ][ j ][ k ];
                kice1 = rhs_ice.x[ i ][ j ][ k ];
                kco1 = rhs_co2.x[ i ][ j ][ k ];

                t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt1 * .5 * dt;
                u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku1 * .5 * dt;
                v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv1 * .5 * dt;
                w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw1 * .5 * dt;
                c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc1 * .5 * dt;
                cloud.x[ i ][ j ][ k ] = cloudn.x[ i ][ j ][ k ] + kcloud1 * .5 * dt;
                ice.x[ i ][ j ][ k ] = icen.x[ i ][ j ][ k ] + kice1 * .5 * dt;
                co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco1 * .5 * dt;

                // Runge-Kutta 4. order for k2 step ( dt )
                RK_RHS_3D_Atmosphere(i, j, k);

                kt2 = rhs_t.x[ i ][ j ][ k ];
                ku2 = rhs_u.x[ i ][ j ][ k ];
                kv2 = rhs_v.x[ i ][ j ][ k ];
                kw2 = rhs_w.x[ i ][ j ][ k ];
                kc2 = rhs_c.x[ i ][ j ][ k ];
                kcloud2 = rhs_cloud.x[ i ][ j ][ k ];
                kice2 = rhs_ice.x[ i ][ j ][ k ];
                kco2 = rhs_co2.x[ i ][ j ][ k ];

                t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt2 * .5 * dt;
                u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku2 * .5 * dt;
                v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv2 * .5 * dt;
                w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw2 * .5 * dt;
                c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc2 * .5 * dt;
                cloud.x[ i ][ j ][ k ] = cloudn.x[ i ][ j ][ k ] + kcloud2 * .5 * dt;
                ice.x[ i ][ j ][ k ] = icen.x[ i ][ j ][ k ] + kice2 * .5 * dt;
                co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco2 * .5 * dt;

                // Runge-Kutta 4. order for k3 step ( dt )
                RK_RHS_3D_Atmosphere(i, j, k);

                kt3 = rhs_t.x[ i ][ j ][ k ];
                ku3 = rhs_u.x[ i ][ j ][ k ];
                kv3 = rhs_v.x[ i ][ j ][ k ];
                kw3 = rhs_w.x[ i ][ j ][ k ];
                kc3 = rhs_c.x[ i ][ j ][ k ];
                kcloud3 = rhs_cloud.x[ i ][ j ][ k ];
                kice3 = rhs_ice.x[ i ][ j ][ k ];
                kco3 = rhs_co2.x[ i ][ j ][ k ];

                t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt3 * dt;
                u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku3 * dt;
                v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv3 * dt;
                w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw3 * dt;
                c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc3 * dt;
                cloud.x[ i ][ j ][ k ] = cloudn.x[ i ][ j ][ k ] + kcloud3 * dt;
                ice.x[ i ][ j ][ k ] = icen.x[ i ][ j ][ k ] + kice3 * dt;
                co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco3 * dt;

                // Runge-Kutta 4. order for k4 step ( dt )
                RK_RHS_3D_Atmosphere(i, j, k);

                kt4 = rhs_t.x[ i ][ j ][ k ];
                ku4 = rhs_u.x[ i ][ j ][ k ];
                kv4 = rhs_v.x[ i ][ j ][ k ];
                kw4 = rhs_w.x[ i ][ j ][ k ];
                kc4 = rhs_c.x[ i ][ j ][ k ];
                kcloud4 = rhs_cloud.x[ i ][ j ][ k ];
                kice4 = rhs_ice.x[ i ][ j ][ k ];
                kco4 = rhs_co2.x[ i ][ j ][ k ];

                t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + dt * ( kt1 + 2. * kt2 + 2. * kt3 + kt4 ) / 6.;
                u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + dt * ( ku1 + 2. * ku2 + 2. * ku3 + ku4 ) / 6.;
                v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + dt * ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
                w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + dt * ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
                c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + dt * ( kc1 + 2. * kc2 + 2. * kc3 + kc4 ) / 6.;
                cloud.x[ i ][ j ][ k ] = cloudn.x[ i ][ j ][ k ] + dt * ( kcloud1 + 2. * kcloud2 + 2. * kcloud3 + kcloud4 ) / 6.;
                ice.x[ i ][ j ][ k ] = icen.x[ i ][ j ][ k ] + dt * ( kice1 + 2. * kice2 + 2. * kice3 + kice4 ) / 6.;
                co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + dt *( kco1 + 2. * kco2 + 2. * kco3 + kco4 ) / 6.;
            }
        }
    }
}

void cAtmosphereModel::solveRungeKutta_2D_Atmosphere(){
// Runge-Kutta 4. order for u, v and w component, temperature, water vapour and co2 content
//  2D surface iterations
// Runge-Kutta 4. order for u, v and w component, temperature and salt concentration
    double kv1, kw1, kv2, kw2, kv3, kw3, kw4, kv4;
    for ( int j = 1; j < jm-1; j++ ){
        for ( int k = 1; k < km-1; k++ ){

            // Runge-Kutta 4. order for k1 step ( dt )
            RK_RHS_2D_Atmosphere( j, k);

            kv1 = rhs_v.x[ 0 ][ j ][ k ];
            kw1 =rhs_w.x[ 0 ][ j ][ k ];

            v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv1 * .5 * dt;
            w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw1 * .5 * dt;

            // Runge-Kutta 4. order for k2 step ( dt )
            RK_RHS_2D_Atmosphere( j, k);

            kv2 = rhs_v.x[ 0 ][ j ][ k ];
            kw2 = rhs_w.x[ 0 ][ j ][ k ];

            v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv2 * .5 * dt;
            w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw2 * .5 * dt;

            // Runge-Kutta 4. order for k3 step ( dt )
            RK_RHS_2D_Atmosphere( j, k);

            kv3 = rhs_v.x[ 0 ][ j ][ k ];
            kw3 = rhs_w.x[ 0 ][ j ][ k ];

            v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv3 * dt;
            w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw3 * dt;

            // Runge-Kutta 4. order for k4 step ( dt )
            RK_RHS_2D_Atmosphere( j, k);

            kv4 = rhs_v.x[ 0 ][ j ][ k ];
            kw4 =rhs_w.x[ 0 ][ j ][ k ];

            v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + dt * ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
            w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + dt * ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
        }
    }
}
