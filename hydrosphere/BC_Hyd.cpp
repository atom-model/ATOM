/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the coordinate system
*/


#include <iostream>
#include <cmath>
#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "cHydrosphereModel.h"
#include "Utils.h"
using namespace std;
using namespace AtomUtils;

void cHydrosphereModel::BC_radius(){
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){

            t.x[0][j][k] = t.x[3][j][k] 
                - 3.0 * t.x[2][j][k] + 3.0 * t.x[1][j][k];  // extrapolation

            u.x[0][j][k] = u.x[3][j][k] 
                - 3.0 * u.x[2][j][k] + 3.0 * u.x[1][j][k];  // extrapolation
/*
            v.x[0][j][k] = v.x[3][j][k] 
                - 3.0 * v.x[2][j][k] + 3.0 * v.x[1][j][k];  // extrapolation
            w.x[0][j][k] = w.x[3][j][k] 
                - 3.0 * w.x[2][j][k] + 3.0 * w.x[1][j][k];  // extrapolation
*/
//            u.x[0][j][k] = 0.0;
            v.x[0][j][k] = 0.0;
            w.x[0][j][k] = 0.0;

            c.x[0][j][k] = c.x[3][j][k] 
                - 3.0 * c.x[2][j][k] + 3.0 * c.x[1][j][k];  // extrapolation

//            t.x[im-1][j][k] = t.x[im-4][j][k] 
//                - 3.0 * t.x[im-3][j][k] + 3.0 * t.x[im-2][j][k];  // extrapolation

            u.x[im-1][j][k] = u.x[im-4][j][k] 
                - 3.0 * u.x[im-3][j][k] + 3.0 * u.x[im-2][j][k];  // extrapolation

            v.x[im-1][j][k] = v.x[im-4][j][k] 
                - 3.0 * v.x[im-3][j][k] + 3.0 * v.x[im-2][j][k];  // extrapolation
            w.x[im-1][j][k] = w.x[im-4][j][k] 
                - 3.0 * w.x[im-3][j][k] + 3.0 * w.x[im-2][j][k];  // extrapolation

            c.x[im-1][j][k] = c.x[im-4][j][k] 
                - 3.0 * c.x[im-3][j][k] + 3.0 * c.x[im-2][j][k];  // extrapolation

/*
            t.x[im-1][j][k] = c43 * t.x[im-2][j][k] - c13 * t.x[im-3][j][k];  // von Neumann
            u.x[im-1][j][k] = c43 * u.x[im-2][j][k] - c13 * u.x[im-3][j][k];  // von Neumann
            v.x[im-1][j][k] = c43 * v.x[im-2][j][k] - c13 * v.x[im-3][j][k];  // von Neumann
            w.x[im-1][j][k] = c43 * w.x[im-2][j][k] - c13 * w.x[im-3][j][k];  // von Neumann
            c.x[im-1][j][k] = c43 * c.x[im-2][j][k] - c13 * c.x[im-3][j][k];  // von Neumann
*/
        }
    }
    return;
}
/*
*
*/
void cHydrosphereModel::BC_theta(){
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
/*
            t.x[i][0][k] = c43 * t.x[i][1][k] - c13 * t.x[i][2][k];
            u.x[i][0][k] = c43 * u.x[i][1][k] - c13 * u.x[i][2][k];
            v.x[i][0][k] = c43 * v.x[i][1][k] - c13 * v.x[i][2][k];
            w.x[i][0][k] = c43 * w.x[i][1][k] - c13 * w.x[i][2][k];
            c.x[i][0][k] = c43 * c.x[i][1][k] - c13 * c.x[i][2][k];

            t.x[i][jm-1][k] = c43 * t.x[i][jm-2][k] - c13 * t.x[i][jm-3][k];
            u.x[i][jm-1][k] = c43 * u.x[i][jm-2][k] - c13 * u.x[i][jm-3][k];
            v.x[i][jm-1][k] = c43 * v.x[i][jm-2][k] - c13 * v.x[i][jm-3][k];
            w.x[i][jm-1][k] = c43 * w.x[i][jm-2][k] - c13 * w.x[i][jm-3][k];
            c.x[i][jm-1][k] = c43 * c.x[i][jm-2][k] - c13 * c.x[i][jm-3][k];
*/
/*
            t.x[i][0][k] = t.x[i][3][k] 
                - 3.0 * t.x[i][2][k] + 3.0 * t.x[i][1][k];  // extrapolation
            u.x[i][0][k] = u.x[i][3][k] 
                - 3.0 * u.x[i][2][k] + 3.0 * u.x[i][1][k];  // extrapolation
            v.x[i][0][k] = v.x[i][3][k] 
                - 3.0 * v.x[i][2][k] + 3.0 * v.x[i][1][k];  // extrapolation
            w.x[i][0][k] = w.x[i][3][k] 
                - 3.0 * w.x[i][2][k] + 3.0 * w.x[i][1][k];  // extrapolation
            c.x[i][0][k] = c.x[i][3][k] 
                - 3.0 * c.x[i][2][k] + 3.0 * c.x[i][1][k];  // extrapolation

            t.x[i][jm-1][k] = t.x[i][jm-4][k] 
                - 3.0 * t.x[i][jm-3][k] + 3.0 * t.x[i][jm-2][k];  // extrapolation
            u.x[i][jm-1][k] = u.x[i][jm-4][k] 
                - 3.0 * u.x[i][jm-3][k] + 3.0 * u.x[i][jm-2][k];  // extrapolation
            v.x[i][jm-1][k] = v.x[i][jm-4][k] 
                - 3.0 * v.x[i][jm-3][k] + 3.0 * v.x[i][jm-2][k];  // extrapolation
            w.x[i][jm-1][k] = w.x[i][jm-4][k] 
                - 3.0 * w.x[i][jm-3][k] + 3.0 * w.x[i][jm-2][k];  // extrapolation
            c.x[i][jm-1][k] = c.x[i][jm-4][k] 
                - 3.0 * c.x[i][jm-3][k] + 3.0 * c.x[i][jm-2][k];  // extrapolation
*/

            t.x[i][0][k] = t.x[i][0][182];  // along k all values are the same north pole value
            t.x[i][jm-1][k] = t.x[i][jm-1][182];  // along k all values are the same south pole value

            u.x[i][0][k] = u.x[i][0][182];  // along k all values are the same north pole value
            u.x[i][jm-1][k] = u.x[i][jm-1][182];  // along k all values are the same south pole value

            v.x[i][0][k] = v.x[i][0][182];  // along k all values are the same north pole value
            v.x[i][jm-1][k] = v.x[i][jm-1][182];  // along k all values are the same south pole value

            w.x[i][0][k] = w.x[i][0][182];  // along k all values are the same north pole value
            w.x[i][jm-1][k] = w.x[i][jm-1][182];  // along k all values are the same south pole value

            c.x[i][0][k] = c.x[i][0][182];  // along k all values are the same north pole value
            c.x[i][jm-1][k] = c.x[i][jm-1][182];  // along k all values are the same south pole value

/*
            u.x[i][0][k] = 0.0;  // along k all values are the same north pole value
            u.x[i][jm-1][k] = 0.0;  // along k all values are the same south pole value

            v.x[i][0][k] = 0.0;  // along k all values are the same north pole value
            v.x[i][jm-1][k] = 0.0;  // along k all values are the same south pole value

            w.x[i][0][k] = 0.0;  // along k all values are the same north pole value
            w.x[i][jm-1][k] = 0.0;  // along k all values are the same south pole value
*/
        }
    }
    return;
}
/*
*
*/
void cHydrosphereModel::BC_phi(){
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            t.x[i][j][0] = c43 * t.x[i][j][1] - c13 * t.x[i][j][2];
            t.x[i][j][km-1] = c43 * t.x[i][j][km-2] - c13 * t.x[i][j][km-3];
            t.x[i][j][0] = t.x[i][j][km-1] = (t.x[i][j][0] + t.x[i][j][km-1])/2.0;

            u.x[i][j][0] = c43 * u.x[i][j][1] - c13 * u.x[i][j][2];
            u.x[i][j][km-1] = c43 * u.x[i][j][km-2] - c13 * u.x[i][j][km-3];
            u.x[i][j][0] = u.x[i][j][km-1] = (u.x[i][j][0] + u.x[i][j][km-1])/2.0;

            v.x[i][j][0] = c43 * v.x[i][j][1] - c13 * v.x[i][j][2];
            v.x[i][j][km-1] = c43 * v.x[i][j][km-2] - c13 * v.x[i][j][km-3];
            v.x[i][j][0] = v.x[i][j][km-1] = (v.x[i][j][0] + v.x[i][j][km-1])/2.0;

            w.x[i][j][0] = c43 * w.x[i][j][1] - c13 * w.x[i][j][2];
            w.x[i][j][km-1] = c43 * w.x[i][j][km-2] - c13 * w.x[i][j][km-3];
            w.x[i][j][0] = w.x[i][j][km-1] = (w.x[i][j][0] + w.x[i][j][km-1])/2.0;

            c.x[i][j][0] = c43 * c.x[i][j][1] - c13 * c.x[i][j][2];
            c.x[i][j][km-1] = c43 * c.x[i][j][km-2] - c13 * c.x[i][j][km-3];
            c.x[i][j][0] = c.x[i][j][km-1] = (c.x[i][j][0] + c.x[i][j][km-1])/2.0;
        }
    }
    return;
}
/*
*
*/
void cHydrosphereModel::BC_SolidGround(){
    cout << endl << "      OGCM: BC_SolidGround" << endl;
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
                    r_water.x[i][j][k] = r_0_water;
                    r_salt_water.x[i][j][k] = r_0_saltwater;
                }
            }
        }
    }
    cout << "      OGCM: BC_SolidGround" << endl;
    return;
}

