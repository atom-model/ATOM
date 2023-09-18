/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include <iostream>
#include <cmath>
#include "cAtmosphereModel.h"

using namespace std;

void cAtmosphereModel::solveRungeKutta_3D_Atmosphere(){
cout << endl << " .................... solveRungeKutta_3D_Atmosphere begin" << endl;
    double kt1, ku1, kv1, kw1, kc1, kcloud1, kice1, kg1, kco1,
           kt2, ku2, kv2, kw2, kc2, kcloud2, kice2, kg2, kco2,
           kt3, ku3, kv3, kw3, kc3, kcloud3, kice3, kg3, kco3,
           kt4, ku4, kv4, kw4, kc4, kcloud4, kice4, kg4, kco4;
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            for(int k = 1; k < km-1; k++){
                RK_RHS_3D_Atmosphere(i, j, k);
                kt1 = dt * rhs_t.x[i][j][k];
                ku1 = dt * rhs_u.x[i][j][k];
                kv1 = dt * rhs_v.x[i][j][k];
                kw1 = dt * rhs_w.x[i][j][k];
                kc1 = dt * rhs_c.x[i][j][k];
                kcloud1 = dt * rhs_cloud.x[i][j][k];
                kice1 = dt * rhs_ice.x[i][j][k];
                kg1 = dt * rhs_g.x[i][j][k];
                kco1 = dt * rhs_co2.x[i][j][k];
                t.x[i][j][k] = tn.x[i][j][k] + kt1 * 0.5;
                u.x[i][j][k] = un.x[i][j][k] + ku1 * 0.5;
                v.x[i][j][k] = vn.x[i][j][k] + kv1 * 0.5;
                w.x[i][j][k] = wn.x[i][j][k] + kw1 * 0.5;
                c.x[i][j][k] = cn.x[i][j][k] + kc1 * 0.5;
                cloud.x[i][j][k] = cloudn.x[i][j][k] + kcloud1 * 0.5;
                ice.x[i][j][k] = icen.x[i][j][k] + kice1 * 0.5;
                gr.x[i][j][k] = grn.x[i][j][k] + kg1 * 0.5;
                co2.x[i][j][k] = co2n.x[i][j][k] + kco1 * 0.5;

                RK_RHS_3D_Atmosphere(i, j, k);
                kt2 = dt * rhs_t.x[i][j][k];
                ku2 = dt * rhs_u.x[i][j][k];
                kv2 = dt * rhs_v.x[i][j][k];
                kw2 = dt * rhs_w.x[i][j][k];
                kc2 = dt * rhs_c.x[i][j][k];
                kcloud2 = dt * rhs_cloud.x[i][j][k];
                kice2 = dt * rhs_ice.x[i][j][k];
                kg2 = dt * rhs_g.x[i][j][k];
                kco2 = dt * rhs_co2.x[i][j][k];
                t.x[i][j][k] = tn.x[i][j][k] + kt2 * 0.5;
                u.x[i][j][k] = un.x[i][j][k] + ku2 * 0.5;
                v.x[i][j][k] = vn.x[i][j][k] + kv2 * 0.5;
                w.x[i][j][k] = wn.x[i][j][k] + kw2 * 0.5;
                c.x[i][j][k] = cn.x[i][j][k] + kc2 * 0.5;
                cloud.x[i][j][k] = cloudn.x[i][j][k] + kcloud2 * 0.5;
                ice.x[i][j][k] = icen.x[i][j][k] + kice2 * 0.5;
                gr.x[i][j][k] = grn.x[i][j][k] + kg2 * 0.5;
                co2.x[i][j][k] = co2n.x[i][j][k] + kco2 * 0.5;

                RK_RHS_3D_Atmosphere(i, j, k);
                kt3 = dt * rhs_t.x[i][j][k];
                ku3 = dt * rhs_u.x[i][j][k];
                kv3 = dt * rhs_v.x[i][j][k];
                kw3 = dt * rhs_w.x[i][j][k];
                kc3 = dt * rhs_c.x[i][j][k];
                kcloud3 = dt * rhs_cloud.x[i][j][k];
                kice3 = dt * rhs_ice.x[i][j][k];
                kg3 = dt * rhs_g.x[i][j][k];
                kco3 = dt * rhs_co2.x[i][j][k];
                t.x[i][j][k] = tn.x[i][j][k] + kt3;
                u.x[i][j][k] = un.x[i][j][k] + ku3;
                v.x[i][j][k] = vn.x[i][j][k] + kv3;
                w.x[i][j][k] = wn.x[i][j][k] + kw3;
                c.x[i][j][k] = cn.x[i][j][k] + kc3;
                cloud.x[i][j][k] = cloudn.x[i][j][k] + kcloud3;
                ice.x[i][j][k] = icen.x[i][j][k] + kice3;
                gr.x[i][j][k] = grn.x[i][j][k] + kg3;
                co2.x[i][j][k] = co2n.x[i][j][k] + kco3;

                RK_RHS_3D_Atmosphere(i, j, k);
                kt4 = dt * rhs_t.x[i][j][k];
                ku4 = dt * rhs_u.x[i][j][k];
                kv4 = dt * rhs_v.x[i][j][k];
                kw4 = dt * rhs_w.x[i][j][k];
                kc4 = dt * rhs_c.x[i][j][k];
                kcloud4 = dt * rhs_cloud.x[i][j][k];
                kice4 = dt * rhs_ice.x[i][j][k];
                kg4 = dt * rhs_g.x[i][j][k];
                kco4 = dt * rhs_co2.x[i][j][k];
                t.x[i][j][k] = tn.x[i][j][k] 
                    + (kt1 + 2.0 * kt2 + 2.0 * kt3 + kt4)/6.0;
                u.x[i][j][k] = un.x[i][j][k] 
                    + (ku1 + 2.0 * ku2 + 2.0 * ku3 + ku4)/6.0;
                v.x[i][j][k] = vn.x[i][j][k] 
                    + (kv1 + 2.0 * kv2 + 2.0 * kv3 + kv4)/6.0;
                w.x[i][j][k] = wn.x[i][j][k] 
                    + (kw1 + 2.0 * kw2 + 2.0 * kw3 + kw4)/6.0;
                c.x[i][j][k] = cn.x[i][j][k] 
                    + (kc1 + 2.0 * kc2 + 2.0 * kc3 + kc4)/6.0;
                cloud.x[i][j][k] = cloudn.x[i][j][k] 
                    + (kcloud1 + 2.0 * kcloud2 + 2.0 * kcloud3 + kcloud4)/6.0;
                ice.x[i][j][k] = icen.x[i][j][k] 
                    + (kice1 + 2.0 * kice2 + 2.0 * kice3 + kice4)/6.0;
                gr.x[i][j][k] = grn.x[i][j][k] 
                    + (kg1 + 2.0 * kg2 + 2.0 * kg3 + kg4)/6.0;
                co2.x[i][j][k] = co2n.x[i][j][k] 
                    + (kco1 + 2.0 * kco2 + 2.0 * kco3 + kco4)/6.0;
            }
        }
    }
cout << " .................... solveRungeKutta_3D_Atmosphere end" << endl;
    return;
}
