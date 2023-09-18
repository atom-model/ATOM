/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <Utils.h>
#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "cHydrosphereModel.h"
#include "AtomMath.h"

using namespace std;
using namespace AtomUtils;

void cHydrosphereModel::solveRungeKutta_3D_Hydrosphere(){
cout << endl << " .................... solveRungeKutta_3D_Hydosphere begin" << endl;
    double kt1, ku1, kv1, kw1, kc1;
    double kt2, ku2, kv2, kw2, kc2;
    double kt3, ku3, kv3, kw3, kc3;
    double kt4, ku4, kv4, kw4, kc4;
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            for(int k = 1; k < km-1; k++){
                RK_RHS_3D_Hydrosphere(i,j,k);
                kt1 = dt * rhs_t.x[i][j][k];
                ku1 = dt * rhs_u.x[i][j][k];
                kv1 = dt * rhs_v.x[i][j][k];
                kw1 = dt * rhs_w.x[i][j][k];
                kc1 = dt * rhs_c.x[i][j][k];
                t.x[i][j][k] = tn.x[i][j][k] + kt1 * 0.5;
                u.x[i][j][k] = un.x[i][j][k] + ku1 * 0.5;
                v.x[i][j][k] = vn.x[i][j][k] + kv1 * 0.5;
                w.x[i][j][k] = wn.x[i][j][k] + kw1 * 0.5;
                c.x[i][j][k] = cn.x[i][j][k] + kc1 * 0.5;
                RK_RHS_3D_Hydrosphere(i, j, k);
                kt2 = dt * rhs_t.x[i][j][k];
                ku2 = dt * rhs_u.x[i][j][k];
                kv2 = dt * rhs_v.x[i][j][k];
                kw2 = dt * rhs_w.x[i][j][k];
                kc2 = dt * rhs_c.x[i][j][k];
                t.x[i][j][k] = tn.x[i][j][k] + kt2 * 0.5;
                u.x[i][j][k] = un.x[i][j][k] + ku2 * 0.5;
                v.x[i][j][k] = vn.x[i][j][k] + kv2 * 0.5;
                w.x[i][j][k] = wn.x[i][j][k] + kw2 * 0.5;
                c.x[i][j][k] = cn.x[i][j][k] + kc2 * 0.5;
                RK_RHS_3D_Hydrosphere(i,j,k);
                kt3 = dt * rhs_t.x[i][j][k];
                ku3 = dt * rhs_u.x[i][j][k];
                kv3 = dt * rhs_v.x[i][j][k];
                kw3 = dt * rhs_w.x[i][j][k];
                kc3 = dt * rhs_c.x[i][j][k];
                t.x[i][j][k] = tn.x[i][j][k] + kt3;
                u.x[i][j][k] = un.x[i][j][k] + ku3;
                v.x[i][j][k] = vn.x[i][j][k] + kv3;
                w.x[i][j][k] = wn.x[i][j][k] + kw3;
                c.x[i][j][k] = cn.x[i][j][k] + kc3;
                RK_RHS_3D_Hydrosphere(i, j, k);
                kt4 = dt * rhs_t.x[i][j][k];
                ku4 = dt * rhs_u.x[i][j][k];
                kv4 = dt * rhs_v.x[i][j][k];
                kw4 = dt * rhs_w.x[i][j][k];
                kc4 = dt * rhs_c.x[i][j][k];
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
            }
        }
    }
cout << endl << " .................... solveRungeKutta_3D_Hydrosphere end" << endl;
    return;
}
/*
*
*/
