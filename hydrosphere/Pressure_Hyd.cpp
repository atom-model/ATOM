/*
 * Atmosphere General Circulation Modell(AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <cmath>
#include "Array.h"
#include "Array_2D.h"
#include "Utils.h"
#include "cHydrosphereModel.h"

using namespace std;
using namespace AtomUtils;

void cHydrosphereModel::computePressure_3D(){
    cout << endl << "      OGCM: computePressure_3D" << endl;
    for(int j = 1; j < jm-1; j++){
        for(int k = 1; k < km-1; k++){
            aux_u.x[0][j][k] = c43 * aux_u.x[1][j][k] - c13 * aux_u.x[2][j][k];
            aux_u.x[im-1][j][k] = c43 * aux_u.x[im-2][j][k] - c13 * aux_u.x[im-3][j][k];
            aux_v.x[0][j][k] = c43 * aux_v.x[1][j][k] - c13 * aux_v.x[2][j][k];
            aux_v.x[im-1][j][k] = c43 * aux_v.x[im-2][j][k] - c13 * aux_v.x[im-3][j][k];
            aux_w.x[0][j][k] = c43 * aux_w.x[1][j][k] - c13 * aux_w.x[2][j][k];
            aux_w.x[im-1][j][k] = c43 * aux_w.x[im-2][j][k] - c13 * aux_w.x[im-3][j][k];
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            aux_u.x[i][0][k] = c43 * aux_u.x[i][1][k] - c13 * aux_u.x[i][2][k];
            aux_u.x[i][jm-1][k] = c43 * aux_u.x[i][jm-2][k] - c13 * aux_u.x[i][jm-3][k];
            aux_v.x[i][0][k] = c43 * aux_v.x[i][1][k] - c13 * aux_v.x[i][2][k];
            aux_v.x[i][jm-1][k] = c43 * aux_v.x[i][jm-2][k] - c13 * aux_v.x[i][jm-3][k];
            aux_w.x[i][0][k] = c43 * aux_w.x[i][1][k] - c13 * aux_w.x[i][2][k];
            aux_w.x[i][jm-1][k] = c43 * aux_w.x[i][jm-2][k] - c13 * aux_w.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            aux_u.x[i][j][0] = c43 * aux_u.x[i][j][1] - c13 * aux_u.x[i][j][2];
            aux_u.x[i][j][km-1] = c43 * aux_u.x[i][j][km-2] - c13 * aux_u.x[i][j][km-3];
            aux_u.x[i][j][0] = aux_u.x[i][j][km-1] = (aux_u.x[i][j][0] + aux_u.x[i][j][km-1])/2.0;
            aux_v.x[i][j][0] = c43 * aux_v.x[i][j][1] - c13 * aux_v.x[i][j][2];
            aux_v.x[i][j][km-1] = c43 * aux_v.x[i][j][km-2] - c13 * aux_v.x[i][j][km-3];
            aux_v.x[i][j][0] = aux_v.x[i][j][km-1] = (aux_v.x[i][j][0] + aux_v.x[i][j][km-1])/2.0;
            aux_w.x[i][j][0] = c43 * aux_w.x[i][j][1] - c13 * aux_w.x[i][j][2];
            aux_w.x[i][j][km-1] = c43 * aux_w.x[i][j][km-2] - c13 * aux_w.x[i][j][km-3];
            aux_w.x[i][j][0] = aux_w.x[i][j][km-1] = (aux_w.x[i][j][0] + aux_w.x[i][j][km-1])/2.0;
        }
    }

    double rm = 0.0;
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double sinthe = 0.0;
    double rmsinthe = 0.0;
    double denom = 0.0;
    double num1 = 0.0;
    double num2 = 0.0;
    double num3 = 0.0;
    double daux_udr = 0.0;
    double daux_vdthe = 0.0;
    double daux_wdphi = 0.0;

    for(int i = 1; i < im-1; i++){
        rm = rad.z[i];
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            if(sinthe == 0.0) sinthe = 1.0e-5;
            rmsinthe = rm * sinthe;
            denom = 2.0/dr2 + 2.0/(rm * dthe2) + 2.0/(rmsinthe * dphi2);
            num1 = 1.0/dr2;
            num2 = 1.0/(rm * dthe2);
            num3 = 1.0/(rmsinthe * dphi2);

            for(int k = 1; k < km-1; k++){
                daux_udr = 
                    (aux_u.x[i+1][j][k] - aux_u.x[i-1][j][k])
                    /(2.0 * dr);
                daux_vdthe = 
                    (aux_v.x[i][j+1][k] - aux_v.x[i][j-1][k])
                    /(2.0 * dthe * rm);
                daux_wdphi = 
                    (aux_w.x[i][j][k+1] - aux_w.x[i][j][k-1])
                    /(2.0 * dphi * rmsinthe);

                p_dyn.x[i][j][k] = 
                     ((p_dyn.x[i+1][j][k] + p_dyn.x[i-1][j][k]) * num1 
                    + (p_dyn.x[i][j+1][k] + p_dyn.x[i][j-1][k]) * num2 
                    + (p_dyn.x[i][j][k+1] + p_dyn.x[i][j][k-1]) * num3 
                    - daux_udr - daux_vdthe - daux_wdphi)/denom;
            }  // end k
        }  // end j
    }  // end i

    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){

            p_dyn.x[i][0][k] = c43 * p_dyn.x[i][1][k] 
                - c13 * p_dyn.x[i][2][k];
            p_dyn.x[i][jm-1][k] = c43 * p_dyn.x[i][jm-2][k] 
                - c13 * p_dyn.x[i][jm-3][k];
        }
    }

    for(int k = 0; k < km; k++){
        for(int j = 1; j < jm-1; j++){
            p_dyn.x[0][j][k] = p_dyn.x[3][j][k] 
                - 3.0 * p_dyn.x[2][j][k] 
                + 3.0 * p_dyn.x[1][j][k];  // extrapolation
            p_dyn.x[im-1][j][k] = p_dyn.x[im-4][j][k] 
                - 3.0 * p_dyn.x[im-3][j][k] 
                + 3.0 * p_dyn.x[im-2][j][k];  // extrapolation
        }
    }
    for(int i = 1; i < im-1; i++){
        for(int j = 1; j < jm-1; j++){
            p_dyn.x[i][j][0] = c43 * p_dyn.x[i][j][1] 
                - c13 * p_dyn.x[i][j][2];
            p_dyn.x[i][j][km-1] = c43 * p_dyn.x[i][j][km-2] 
                - c13 * p_dyn.x[i][j][km-3];
            p_dyn.x[i][j][0] = p_dyn.x[i][j][km-1] =
               (p_dyn.x[i][j][0] + p_dyn.x[i][j][km-1])/2.0;
        }
    }
    cout << "      OGCM: computePressure_3D ended" << endl;
    return;
}
/*
*
*/
