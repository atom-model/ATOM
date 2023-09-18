/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/
#include <iostream>
#include <cmath>
#include "cAtmosphereModel.h"
#include "Utils.h"
using namespace std;
using namespace AtomUtils;

void cAtmosphereModel::BC_radius(){
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
/*
            t.x[0][j][k] = c43 * t.x[1][j][k] - c13 * t.x[2][j][k];  // von Neumann
            u.x[0][j][k] = c43 * u.x[1][j][k] - c13 * u.x[2][j][k];
            v.x[0][j][k] = c43 * v.x[1][j][k] - c13 * v.x[2][j][k];
            w.x[0][j][k] = c43 * w.x[1][j][k] - c13 * w.x[2][j][k];
*/

            t.x[0][j][k] = t.x[3][j][k] 
                - 3.0 * t.x[2][j][k] + 3.0 * t.x[1][j][k];  // extrapolation
            u.x[0][j][k] = u.x[3][j][k] 
                - 3.0 * u.x[2][j][k] + 3.0 * u.x[1][j][k];  // extrapolation
            v.x[0][j][k] = v.x[3][j][k] 
                - 3.0 * v.x[2][j][k] + 3.0 * v.x[1][j][k];  // extrapolation
            w.x[0][j][k] = w.x[3][j][k] 
                - 3.0 * w.x[2][j][k] + 3.0 * w.x[1][j][k];  // extrapolation
/*
            t.x[im-1][j][k] = t.x[im-4][j][k] 
                - 3.0 * t.x[im-3][j][k] + 3.0 * t.x[im-2][j][k];  // extrapolation
            u.x[im-1][j][k] = u.x[im-4][j][k] 
                - 3.0 * u.x[im-3][j][k] + 3.0 * u.x[im-2][j][k];  // extrapolation
            v.x[im-1][j][k] = v.x[im-4][j][k] 
                - 3.0 * v.x[im-3][j][k] + 3.0 * v.x[im-2][j][k];  // extrapolation
            w.x[im-1][j][k] = w.x[im-4][j][k] 
                 - 3.0 * w.x[im-3][j][k] + 3.0 * w.x[im-2][j][k];  // extrapolation
*/
            u.x[im-1][j][k] = 0.0;
            v.x[im-1][j][k] = 0.0;
            w.x[im-1][j][k] = 0.0;
/*
            c.x[im-1][j][k] = c.x[im-4][j][k] 
                - 3.0 * c.x[im-3][j][k] + 3.0 * c.x[im-2][j][k];  // extrapolation
            cloud.x[im-1][j][k] = cloud.x[im-4][j][k] 
                - 3.0 * cloud.x[im-3][j][k] + 3.0 * cloud.x[im-2][j][k];  // extrapolation
            ice.x[im-1][j][k] = ice.x[im-4][j][k] 
                - 3.0 * ice.x[im-3][j][k] + 3.0 * ice.x[im-2][j][k];  // extrapolation
            co2.x[im-1][j][k] = co2.x[im-4][j][k] 
                - 3.0 * co2.x[im-3][j][k] + 3.0 * co2.x[im-2][j][k];  // extrapolation
*/
            c.x[0][j][k] = c.x[3][j][k] 
                - 3.0 * c.x[2][j][k] + 3.0 * c.x[1][j][k];  // extrapolation
            cloud.x[0][j][k] = cloud.x[3][j][k] 
                - 3.0 * cloud.x[2][j][k] + 3.0 * cloud.x[1][j][k];  // extrapolation
            ice.x[0][j][k] = ice.x[3][j][k] 
                - 3.0 * ice.x[2][j][k] + 3.0 * ice.x[1][j][k];  // extrapolation
            gr.x[0][j][k] = gr.x[3][j][k] 
                - 3.0 * gr.x[2][j][k] + 3.0 * gr.x[1][j][k];  // extrapolation
            co2.x[0][j][k] = co2.x[3][j][k] 
                - 3.0 * co2.x[2][j][k] + 3.0 * co2.x[1][j][k];  // extrapolation
/*
            c.x[0][j][k] = c43 * c.x[1][j][k] - c13 * c.x[2][j][k];  // von Neumann
            cloud.x[0][j][k] = c43 * cloud.x[1][j][k] - c13 * cloud.x[2][j][k];
            ice.x[0][j][k] = c43 * ice.x[1][j][k] - c13 * ice.x[2][j][k];
            co2.x[0][j][k] = c43 * co2.x[1][j][k] - c13 * co2.x[2][j][k];
*/
            c.x[im-1][j][k] = c43 * c.x[im-2][j][k] - c13 * c.x[im-3][j][k];  // von Neumann
            cloud.x[im-1][j][k] = c43 * cloud.x[im-2][j][k] - c13 * cloud.x[im-3][j][k];
            ice.x[im-1][j][k] = c43 * ice.x[im-2][j][k] - c13 * ice.x[im-3][j][k];
            gr.x[im-1][j][k] = c43 * gr.x[im-2][j][k] - c13 * gr.x[im-3][j][k];
            co2.x[im-1][j][k] = c43 * co2.x[im-2][j][k] - c13 * co2.x[im-3][j][k];

        }
    }
    return;
}
/*
*
*/
void cAtmosphereModel::BC_theta(){
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            t.x[i][0][k] = c43 * t.x[i][1][k] - c13 * t.x[i][2][k];
            t.x[i][jm-1][k] = c43 * t.x[i][jm-2][k] - c13 * t.x[i][jm-3][k];
/*
            u.x[i][0][k] = c43 * u.x[i][1][k] - c13 * u.x[i][2][k];
            u.x[i][jm-1][k] = c43 * u.x[i][jm-2][k] - c13 * u.x[i][jm-3][k];

            v.x[i][0][k] = c43 * v.x[i][1][k] - c13 * v.x[i][2][k];
            v.x[i][jm-1][k] = c43 * v.x[i][jm-2][k] - c13 * v.x[i][jm-3][k];

            w.x[i][0][k] = c43 * w.x[i][1][k] - c13 * w.x[i][2][k];
            w.x[i][jm-1][k] = c43 * w.x[i][jm-2][k] - c13 * w.x[i][jm-3][k];
*/
            r_dry.x[i][0][k] = c43 * r_dry.x[i][1][k] - c13 * r_dry.x[i][2][k];
            r_dry.x[i][jm-1][k] = c43 * r_dry.x[i][jm-2][k] - c13 * r_dry.x[i][jm-3][k];

            r_humid.x[i][0][k] = c43 * r_humid.x[i][1][k] - c13 * r_humid.x[i][2][k];
            r_humid.x[i][jm-1][k] = c43 * r_humid.x[i][jm-2][k] - c13 * r_humid.x[i][jm-3][k];

            c.x[i][0][k] = c43 * c.x[i][1][k] - c13 * c.x[i][2][k];
            c.x[i][jm-1][k] = c43 * c.x[i][jm-2][k] - c13 * c.x[i][jm-3][k];

            cloud.x[i][0][k] = c43 * cloud.x[i][1][k] - c13 * cloud.x[i][2][k];
            cloud.x[i][jm-1][k] = c43 * cloud.x[i][jm-2][k] - c13 * cloud.x[i][jm-3][k];

            ice.x[i][0][k] = c43 * ice.x[i][1][k] - c13 * ice.x[i][2][k];
            ice.x[i][jm-1][k] = c43 * ice.x[i][jm-2][k] - c13 * ice.x[i][jm-3][k];

            gr.x[i][0][k] = c43 * gr.x[i][1][k] - c13 * gr.x[i][2][k];
            gr.x[i][jm-1][k] = c43 * gr.x[i][jm-2][k] - c13 * gr.x[i][jm-3][k];

            co2.x[i][0][k] = c43 * co2.x[i][1][k] - c13 * co2.x[i][2][k];
            co2.x[i][jm-1][k] = c43 * co2.x[i][jm-2][k] - c13 * co2.x[i][jm-3][k];
        }
    }
    return;
}
/*
*
*/
void cAtmosphereModel::BC_phi(){
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            t.x[i][j][0] = c43 * t.x[i][j][1] - c13 * t.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            t.x[i][j][km-1] = c43 * t.x[i][j][km-2] - c13 * t.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            t.x[i][j][0] = t.x[i][j][1] - 3.0 * t.x[i][j][2] + 3.0 * t.x[i][j][3];  // extrapolation
//            t.x[i][j][km-1] = t.x[im-4][j][km-4] - 3.0 * t.x[im-3][j][km-3] + 3.0 * t.x[im-2][j][km-2];  // extrapolation
            t.x[i][j][0] = t.x[i][j][km-1] = (t.x[i][j][0] + t.x[i][j][km-1])/2.0;

            u.x[i][j][0] = c43 * u.x[i][j][1] - c13 * u.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            u.x[i][j][km-1] = c43 * u.x[i][j][km-2] - c13 * u.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            u.x[i][j][0] = u.x[i][j][1] - 3.0 * u.x[i][j][2] + 3.0 * u.x[i][j][3];  // extrapolation
//            u.x[i][j][km-1] = u.x[im-4][j][km-4] - 3.0 * u.x[im-3][j][km-3] + 3.0 * u.x[im-2][j][km-2];  // extrapolation
            u.x[i][j][0] = u.x[i][j][km-1] = (u.x[i][j][0] + u.x[i][j][km-1])/2.0;

            v.x[i][j][0] = c43 * v.x[i][j][1] - c13 * v.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            v.x[i][j][km-1] = c43 * v.x[i][j][km-2] - c13 * v.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            v.x[i][j][0] = v.x[i][j][1] - 3.0 * v.x[i][j][2] + 3.0 * v.x[i][j][3];  // extrapolation
//            v.x[i][j][km-1] = v.x[im-4][j][km-4] - 3.0 * v.x[im-3][j][km-3] + 3.0 * v.x[im-2][j][km-2];  // extrapolation
            v.x[i][j][0] = v.x[i][j][km-1] = (v.x[i][j][0] + v.x[i][j][km-1])/2.0;

            w.x[i][j][0] = c43 * w.x[i][j][1] - c13 * w.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            w.x[i][j][km-1] = c43 * w.x[i][j][km-2] - c13 * w.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            w.x[i][j][0] = w.x[i][j][1] - 3.0 * w.x[i][j][2] + 3.0 * w.x[i][j][3];  // extrapolation
//            w.x[i][j][km-1] = w.x[im-4][j][km-4] - 3.0 * w.x[im-3][j][km-3] + 3.0 * w.x[im-2][j][km-2];  // extrapolation
            w.x[i][j][0] = w.x[i][j][km-1] = (w.x[i][j][0] + w.x[i][j][km-1])/2.0;

            c.x[i][j][0] = c43 * c.x[i][j][1] - c13 * c.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            c.x[i][j][km-1] = c43 * c.x[i][j][km-2] - c13 * c.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            c.x[i][j][0] = c.x[i][j][1] - 3.0 * c.x[i][j][2] + 3.0 * c.x[i][j][3];  // extrapolation
//            c.x[i][j][km-1] = c.x[im-4][j][km-4] - 3.0 * c.x[im-3][j][km-3] + 3.0 * c.x[im-2][j][km-2];  // extrapolation
            c.x[i][j][0] = c.x[i][j][km-1] = (c.x[i][j][0] + c.x[i][j][km-1])/2.0;

            cloud.x[i][j][0] = c43 * cloud.x[i][j][1] - c13 * cloud.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            cloud.x[i][j][km-1] = c43 * cloud.x[i][j][km-2] - c13 * cloud.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            cloud.x[i][j][0] = cloud.x[i][j][1] - 3.0 * cloud.x[i][j][2] + 3.0 * cloud.x[i][j][3];  // extrapolation
//            cloud.x[i][j][km-1] = cloud.x[im-4][j][km-4] - 3.0 * cloud.x[im-3][j][km-3] + 3.0 * cloud.x[im-2][j][km-2];  // extrapolation
            cloud.x[i][j][0] = cloud.x[i][j][km-1] = (cloud.x[i][j][0] + cloud.x[i][j][km-1])/2.0;

            ice.x[i][j][0] = c43 * ice.x[i][j][1] - c13 * ice.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            ice.x[i][j][km-1] = c43 * ice.x[i][j][km-2] - c13 * ice.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            ice.x[i][j][0] = ice.x[i][j][1] - 3.0 * ice.x[i][j][2] + 3.0 * ice.x[i][j][3];  // extrapolation
//            ice.x[i][j][km-1] = ice.x[im-4][j][km-4] - 3.0 * ice.x[im-3][j][km-3] + 3.0 * ice.x[im-2][j][km-2];  // extrapolation
            ice.x[i][j][0] = ice.x[i][j][km-1] = (ice.x[i][j][0] + ice.x[i][j][km-1])/2.0;

            gr.x[i][j][0] = c43 * gr.x[i][j][1] - c13 * gr.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            gr.x[i][j][km-1] = c43 * gr.x[i][j][km-2] - c13 * gr.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            gr.x[i][j][0] = gr.x[i][j][1] - 3.0 * gr.x[i][j][2] + 3.0 * gr.x[i][j][3];  // extrapolation
//            gr.x[i][j][km-1] = gr.x[im-4][j][km-4] - 3.0 * gr.x[im-3][j][km-3] + 3.0 * gr.x[im-2][j][km-2];  // extrapolation
            gr.x[i][j][0] = gr.x[i][j][km-1] = (gr.x[i][j][0] + gr.x[i][j][km-1])/2.0;

            co2.x[i][j][0] = c43 * co2.x[i][j][1] - c13 * co2.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            co2.x[i][j][km-1] = c43 * co2.x[i][j][km-2] - c13 * co2.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            co2.x[i][j][0] = co2.x[i][j][1] - 3.0 * co2.x[i][j][2] + 3.0 * co2.x[i][j][3];  // extrapolation
//            co2.x[i][j][km-1] = co2.x[im-4][j][km-4] - 3.0 * co2.x[im-3][j][km-3] + 3.0 * co2.x[im-2][j][km-2];  // extrapolation
            co2.x[i][j][0] = co2.x[i][j][km-1] = (co2.x[i][j][0] + co2.x[i][j][km-1])/2.0;

            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] - c13 * P_rain.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] - c13 * P_rain.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            P_rain.x[i][j][0] = P_rain.x[i][j][1] - 3.0 * P_rain.x[i][j][2] + 3.0 * P_rain.x[i][j][3];  // extrapolation
//            P_rain.x[i][j][km-1] = P_rain.x[im-4][j][km-4] - 3.0 * P_rain.x[im-3][j][km-3] + 3.0 * P_rain.x[im-2][j][km-2];  // extrapolation
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] = (P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/2.0;

            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] - c13 * P_snow.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] - c13 * P_snow.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            P_snow.x[i][j][0] = P_snow.x[i][j][1] - 3.0 * P_snow.x[i][j][2] + 3.0 * P_snow.x[i][j][3];  // extrapolation
//            P_snow.x[i][j][km-1] = P_snow.x[im-4][j][km-4] - 3.0 * P_snow.x[im-3][j][km-3] + 3.0 * P_snow.x[im-2][j][km-2];  // extrapolation
            P_snow.x[i][j][0] = P_snow.x[i][j][km-1] = (P_snow.x[i][j][0] + P_snow.x[i][j][km-1])/2.0;

            P_graupel.x[i][j][0] = c43 * P_graupel.x[i][j][1] - c13 * P_graupel.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_graupel.x[i][j][km-1] = c43 * P_graupel.x[i][j][km-2] - c13 * P_graupel.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
//            P_graupel.x[i][j][0] = P_graupel.x[i][j][1] - 3.0 * P_graupel.x[i][j][2] + 3.0 * P_graupel.x[i][j][3];  // extrapolation
//            P_graupel.x[i][j][km-1] = P_graupel.x[im-4][j][km-4] - 3.0 * P_graupel.x[im-3][j][km-3] + 3.0 * P_graupel.x[im-2][j][km-2];  // extrapolation
            P_graupel.x[i][j][0] = P_graupel.x[i][j][km-1] = (P_graupel.x[i][j][0] + P_graupel.x[i][j][km-1])/2.0;
        }
    }
    return;
}
/*
*
*/
void cAtmosphereModel::BC_SolidGround(){
//    cout << endl << "      AGCM: BC_SolidGround" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                if(is_land(h, i, j, k)){
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;
                    un.x[i][j][k] = 0.0;
                    vn.x[i][j][k] = 0.0;
                    wn.x[i][j][k] = 0.0;
                }
            }
        }
    }
//    cout << "      AGCM: BC_SolidGround ended" << endl;
    return;
}
