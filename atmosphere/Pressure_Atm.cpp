/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
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
#include "cAtmosphereModel.h"

using namespace std;
using namespace AtomUtils;

void cAtmosphereModel::computePressure_3D(){
    cout << endl << "      AGCM: computePressure_3D" << endl;
    for(int j = 1; j < jm-1; j++){
        for(int k = 1; k < km-1; k++){
            aux_u.x[0][j][k] = c43 * aux_u.x[1][j][k] - c13 * aux_u.x[2][j][k];
            aux_u.x[im-1][j][k] = c43 * aux_u.x[im-2][j][k] - c13 * aux_u.x[im-3][j][k];
            aux_v.x[0][j][k] = c43 * aux_v.x[1][j][k] - c13 * aux_v.x[2][j][k];
            aux_v.x[im-1][j][k] = c43 * aux_v.x[im-2][j][k] - c13 * aux_v.x[im-3][j][k];
            aux_w.x[0][j][k] = c43 * aux_w.x[1][j][k] - c13 * aux_w.x[2][j][k];
            aux_w.x[im-1][j][k] = c43 * aux_w.x[im-2][j][k] - c13 * aux_w.x[im-3][j][k];

            rhs_u.x[0][j][k] = c43 * rhs_u.x[1][j][k] - c13 * rhs_u.x[2][j][k];
            rhs_u.x[im-1][j][k] = c43 * rhs_u.x[im-2][j][k] - c13 * rhs_u.x[im-3][j][k];
            rhs_v.x[0][j][k] = c43 * rhs_v.x[1][j][k] - c13 * rhs_v.x[2][j][k];
            rhs_v.x[im-1][j][k] = c43 * rhs_v.x[im-2][j][k] - c13 * rhs_v.x[im-3][j][k];
            rhs_w.x[0][j][k] = c43 * rhs_w.x[1][j][k] - c13 * rhs_w.x[2][j][k];
            rhs_w.x[im-1][j][k] = c43 * rhs_w.x[im-2][j][k] - c13 * rhs_w.x[im-3][j][k];
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

            rhs_u.x[i][0][k] = c43 * rhs_u.x[i][1][k] - c13 * rhs_u.x[i][2][k];
            rhs_u.x[i][jm-1][k] = c43 * rhs_u.x[i][jm-2][k] - c13 * rhs_u.x[i][jm-3][k];
            rhs_v.x[i][0][k] = c43 * rhs_v.x[i][1][k] - c13 * rhs_v.x[i][2][k];
            rhs_v.x[i][jm-1][k] = c43 * rhs_v.x[i][jm-2][k] - c13 * rhs_v.x[i][jm-3][k];
            rhs_w.x[i][0][k] = c43 * rhs_w.x[i][1][k] - c13 * rhs_w.x[i][2][k];
            rhs_w.x[i][jm-1][k] = c43 * rhs_w.x[i][jm-2][k] - c13 * rhs_w.x[i][jm-3][k];
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

            rhs_u.x[i][j][0] = c43 * rhs_u.x[i][j][1] - c13 * rhs_u.x[i][j][2];
            rhs_u.x[i][j][km-1] = c43 * rhs_u.x[i][j][km-2] - c13 * rhs_u.x[i][j][km-3];
            rhs_u.x[i][j][0] = rhs_u.x[i][j][km-1] = (rhs_u.x[i][j][0] + rhs_u.x[i][j][km-1])/2.0;
            rhs_v.x[i][j][0] = c43 * rhs_v.x[i][j][1] - c13 * rhs_v.x[i][j][2];
            rhs_v.x[i][j][km-1] = c43 * rhs_v.x[i][j][km-2] - c13 * rhs_v.x[i][j][km-3];
            rhs_v.x[i][j][0] = rhs_v.x[i][j][km-1] = (rhs_v.x[i][j][0] + rhs_v.x[i][j][km-1])/2.0;
            rhs_w.x[i][j][0] = c43 * rhs_w.x[i][j][1] - c13 * rhs_w.x[i][j][2];
            rhs_w.x[i][j][km-1] = c43 * rhs_w.x[i][j][km-2] - c13 * rhs_w.x[i][j][km-3];
            rhs_w.x[i][j][0] = rhs_w.x[i][j][km-1] = (rhs_w.x[i][j][0] + rhs_w.x[i][j][km-1])/2.0;
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                if(is_land(h, i, j, k)){
                    rhs_u.x[i][j][k] = 0.0;
                    rhs_v.x[i][j][k] = 0.0;
                    rhs_w.x[i][j][k] = 0.0;
                    aux_u.x[i][j][k] = 0.0;
                    aux_v.x[i][j][k] = 0.0;
                    aux_w.x[i][j][k] = 0.0;
                }
            }
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
    double drhs_udr = 0.0;
    double drhs_vdthe = 0.0;
    double drhs_wdphi = 0.0;
    for(int i = 1; i < im-1; i++){
        rm = rad.z[i];
        double exp_rm = 1.0/(rm + 1.0);
        double exp_2_rm = 1.0/pow((rm + 1.0),2);
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            if(sinthe == 0.0) sinthe = 1.0e-5;
            rmsinthe = rm * sinthe;
            denom = 2.0/dr2 * exp_2_rm + 2.0/(rm * dthe2) + 2.0/(rmsinthe * dphi2);
            num1 = 1.0/dr2 * exp_2_rm;
            num2 = 1.0/(rm * dthe2);
            num3 = 1.0/(rmsinthe * dphi2);
            for(int k = 1; k < km-1; k++){
                daux_udr = (aux_u.x[i+1][j][k] - aux_u.x[i-1][j][k])
                    /(2.0 * dr) * exp_rm;
                drhs_udr = (rhs_u.x[i+1][j][k] - rhs_u.x[i-1][j][k])
                    /(2.0 * dr) * exp_rm;

                if(i <= im-3){
                    if((is_land(h, i, j, k))
                        &&(is_air(h, i+1, j, k))){        
                        daux_udr = (- 3.0 * aux_u.x[i][j][k] 
                            + 4.0 * aux_u.x[i+1][j][k] 
                            - aux_u.x[i+2][j][k])/(2.0 * dr) * exp_rm;
                        drhs_udr = (- 3.0 * rhs_u.x[i][j][k] 
                            + 4.0 * rhs_u.x[i+1][j][k] 
                            - rhs_u.x[i+2][j][k])/(2.0 * dr) * exp_rm;
                    }
                }else{
                    daux_udr = (aux_u.x[i+1][j][k] - aux_u.x[i][j][k])
                        /dr * exp_rm;
                    drhs_udr = (rhs_u.x[i+1][j][k] - rhs_u.x[i][j][k])
                        /dr * exp_rm;
                }

                daux_vdthe = (aux_v.x[i][j+1][k] 
                    - aux_v.x[i][j-1][k])/(2.0 * dthe * rm);
                drhs_vdthe = (rhs_v.x[i][j+1][k] 
                    - rhs_v.x[i][j-1][k])/(2.0 * dthe * rm);

                if((is_land(h, i, j, k))
                    &&(is_air(h, i, j+1, k))){
                    daux_vdthe = (aux_v.x[i][j+1][k] 
                        - aux_v.x[i][j][k])/(dthe * rm);
                    drhs_vdthe = (rhs_v.x[i][j+1][k] 
                        - rhs_v.x[i][j][k])/(dthe * rm);
                    }
                if((is_land(h, i, j, k))
                    &&(is_air(h, i, j-1, k))){
                    daux_vdthe = (aux_v.x[i][j-1][k] 
                        - aux_v.x[i][j][k])/(dthe * rm);
                    drhs_vdthe = (rhs_v.x[i][j-1][k] 
                        - rhs_v.x[i][j][k])/(dthe * rm);
                    }

                if((j >= 2)&&(j <= jm-3)){
                    if((is_land(h, i, j, k))
                        &&(is_air(h, i, j+1, k)) 
                        &&(is_air(h, i, j+2, k))){
                        daux_vdthe = (- 3.0 * aux_v.x[i][j][k] 
                            + 4.0 * aux_v.x[i][j+1][k] 
                            - aux_v.x[i][j+2][k])/(2.0 * dthe * rm);
                        drhs_vdthe = (- 3.0 * rhs_v.x[i][j][k] 
                            + 4.0 * rhs_v.x[i][j+1][k] 
                            - rhs_v.x[i][j+2][k])/(2.0 * dthe * rm);
                    }
                    if((is_land(h, i, j, k))
                        &&(is_air(h, i, j-1, k)) 
                        &&(is_air(h, i, j-2, k))){
                        daux_vdthe = - (- 3.0 * aux_v.x[i][j][k] 
                            + 4.0 * aux_v.x[i][j-1][k] 
                            - aux_v.x[i][j-2][k])/(2.0 * dthe * rm);
                        drhs_vdthe = - (- 3.0 * rhs_v.x[i][j][k] 
                            + 4.0 * rhs_v.x[i][j-1][k] 
                            - rhs_v.x[i][j-2][k])/(2.0 * dthe * rm);
                    }
                }

                daux_wdphi = (aux_w.x[i][j][k+1] 
                    - aux_w.x[i][j][k-1])/(2.0 * dphi * rmsinthe);
                drhs_wdphi = (rhs_w.x[i][j][k+1] 
                    - rhs_w.x[i][j][k-1])/(2.0 * dphi * rmsinthe);

                if((is_land(h, i, j, k))
                    &&(is_air(h, i, j, k+1))){
                    daux_wdphi = (aux_w.x[i][j][k+1] 
                        - aux_w.x[i][j][k])/(dphi * rmsinthe);
                    drhs_wdphi = (rhs_w.x[i][j][k+1] 
                        - rhs_w.x[i][j][k])/(dphi * rmsinthe);
                    }
                if((is_land(h, i, j, k))
                    &&(is_air(h, i, j, k-1))){
                    daux_wdphi = (aux_w.x[i][j][k-1] 
                        - aux_w.x[i][j][k])/(dphi * rmsinthe);
                    drhs_wdphi = (rhs_w.x[i][j][k-1] 
                        - rhs_w.x[i][j][k])/(dphi * rmsinthe);
                    }

                if((k >= 2)&&(k <= km-3)){
                    if((is_land(h, i, j, k))
                        &&(is_air(h, i, j, k+1)) 
                        &&(is_air(h, i, j, k+2))){
                        daux_wdphi = (- 3.0 * aux_w.x[i][j][k] 
                            + 4.0 * aux_w.x[i][j][k+1] 
                            - aux_w.x[i][j][k+2])/(2.0 * rmsinthe * dphi);
                        drhs_wdphi = (- 3.0 * rhs_w.x[i][j][k] 
                            + 4.0 * rhs_w.x[i][j][k+1] 
                            - rhs_w.x[i][j][k+2])/(2.0 * rmsinthe * dphi);
                    }
                    if((is_land(h, i, j, k))
                        &&(is_air(h, i, j, k-1)) 
                        &&(is_air(h, i, j, k-2))){
                        daux_wdphi = - (- 3.0 * aux_w.x[i][j][k] 
                            + 4.0 * aux_w.x[i][j][k-1] 
                            - aux_w.x[i][j][k-2])/(2.0 * rmsinthe * dphi);
                        drhs_wdphi = - (- 3.0 * rhs_w.x[i][j][k] 
                            + 4.0 * rhs_w.x[i][j][k-1] 
                            - rhs_w.x[i][j][k-2])/(2.0 * rmsinthe * dphi);
                    }
                }

                if(((j <= 10)||(j <= jm-11))||((k <= 10)||(k <= km-11))){ // prevents oscillations along pole regions
                    p_stat.x[i][j][k] = 
                         ((p_stat.x[i+1][j][k] + p_stat.x[i-1][j][k]) * num1 
                        + (p_stat.x[i][j+1][k] + p_stat.x[i][j-1][k]) * num2 
                        + (p_stat.x[i][j][k+1] + p_stat.x[i][j][k-1]) * num3 
                        )/denom;
                }else{
                    p_stat.x[i][j][k] = 
                         ((p_stat.x[i+1][j][k] + p_stat.x[i-1][j][k]) * num1 
                        + (p_stat.x[i][j+1][k] + p_stat.x[i][j-1][k]) * num2 
                        + (p_stat.x[i][j][k+1] + p_stat.x[i][j][k-1]) * num3 
                        - (daux_udr - drhs_udr) 
                        - (daux_vdthe - drhs_vdthe)
                        - (daux_wdphi - drhs_wdphi)
                        )/denom;
                }
                if(is_land(h, i, j, k))  p_stat.x[i][j][k] = 0.0;
                if(p_stat.x[i][j][k] < 0.0) p_stat.x[i][j][k] = 0.0;
/*
    cout.precision(12);
    if((j == 75) &&(k == 180)) cout << "northern hemisphere" << endl
//    if((j == 1) &&(k == 180)) cout << "northern hemisphere" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   dt = " << dt 
        << "   dr = " << dr 
        << "   dthe = " << dthe 
        << "   dphi = " << dphi << endl
        << "   rm > rad = " << rm 
        << "   sinthe = " << sinthe
        << "   rmsinthe = " << rmsinthe << endl
        << "   num1 = " << num1
        << "   num2 = " << num2
        << "   num3 = " << num3
        << "   denom = " << denom << endl

        << "   aux_u = " << aux_u.x[i][j][k]
        << "   aux_v = " << aux_v.x[0][j][k]
        << "   aux_w = " << aux_w.x[0][j][k] << endl

        << "   rhs_u = " << rhs_u.x[i][j][k]
        << "   rhs_v = " << rhs_v.x[0][j][k]
        << "   rhs_w = " << rhs_w.x[0][j][k] << endl

        << "   daux_udr = " << daux_udr
        << "   drhs_udr = " << drhs_udr << endl
        << "   daux_vdthe = " << daux_vdthe
        << "   drhs_vdthe = " << drhs_vdthe << endl
        << "   daux_wdphi = " << daux_wdphi
        << "   drhs_wdphi = " << drhs_wdphi << endl

        << "   p_dyni-1 = " << p_stat.x[i-1][j][k]
        << "   p_dyni   = " << p_stat.x[i][j][k]
        << "   p_dyni+1 = " << p_stat.x[i+1][j][k] << endl

        << "   p_dynj-1 = " << p_stat.x[0][j-1][k]
        << "   p_dynj   = " << p_stat.x[0][j][k]
        << "   p_dynj+1 = " << p_stat.x[0][j+1][k] << endl

        << "   p_dynk-1 = " << p_stat.x[0][j][k-1]
        << "   p_dynk   = " << p_stat.x[0][j][k]
        << "   p_dynk+1 = " << p_stat.x[0][j][k+1] << endl

        << "   ui-1 = " << u.x[i-1][j][k]
        << "   ui   = " << u.x[i][j][k]
        << "   ui+1 = " << u.x[i+1][j][k] << endl

        << "   uj-1 = " << u.x[i][j-1][k]
        << "   uj   = " << u.x[i][j][k]
        << "   uj+1 = " << u.x[i][j+1][k] << endl

        << "   uk-1 = " << u.x[i][j][k-1]
        << "   uk   = " << u.x[i][j][k]
        << "   uk+1 = " << u.x[i][j][k+1] << endl << endl

        << "   vi-1 = " << v.x[i-1][j][k]
        << "   vi   = " << v.x[i][j][k]
        << "   vi+1 = " << v.x[i+1][j][k] << endl

        << "   vj-1 = " << v.x[i][j-1][k]
        << "   vj   = " << v.x[i][j][k]
        << "   vj+1 = " << v.x[i][j+1][k] << endl

        << "   vk-1 = " << v.x[i][j][k-1]
        << "   vk   = " << v.x[i][j][k]
        << "   vk+1 = " << v.x[i][j][k+1] << endl << endl

        << "   wi-1 = " << w.x[i-1][j][k]
        << "   wi   = " << w.x[i][j][k]
        << "   wi+1 = " << w.x[i+1][j][k] << endl

        << "   wj-1 = " << w.x[i][j-1][k]
        << "   wj   = " << w.x[i][j][k]
        << "   wj+1 = " << w.x[i][j+1][k] << endl

        << "   wk-1 = " << w.x[i][j][k-1]
        << "   wk   = " << w.x[i][j][k]
        << "   wk+1 = " << w.x[i][j][k+1] << endl << endl;
*/
            }  // end k
        }  // end j
    }  // end i



    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            p_stat.x[0][j][k] = p_stat.x[3][j][k] 
                - 3. * p_stat.x[2][j][k] 
                + 3. * p_stat.x[1][j][k];  // extrapolation
            p_stat.x[im-1][j][k] = p_stat.x[im-4][j][k] 
                - 3.0 * p_stat.x[im-3][j][k] 
                + 3.0 * p_stat.x[im-2][j][k];
            if(is_land(h, 0, j, k))  p_stat.x[0][j][k] = 0.0;
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            p_stat.x[i][0][k] = c43 * p_stat.x[i][1][k] -
                c13 * p_stat.x[i][2][k];
            p_stat.x[i][jm-1][k] = c43 * p_stat.x[i][jm-2][k] -
                c13 * p_stat.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            p_stat.x[i][j][0] = c43 * p_stat.x[i][j][1] -
                c13 * p_stat.x[i][j][2];
            p_stat.x[i][j][km-1] = c43 * p_stat.x[i][j][km-2] -
                c13 * p_stat.x[i][j][km-3];
            p_stat.x[i][j][0] = p_stat.x[i][j][km-1] =
               (p_stat.x[i][j][0] + p_stat.x[i][j][km-1])/2.;
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im-1; i++){
                if(p_stat.x[i][j][k] >= 0.15) p_stat.x[i][j][k] = 0.15;
                if(is_land(h, i, j, k)){
                    p_stat.x[i][j][k] = 0.0;
                }
            }
        }
    }
    cout << "      AGCM: computePressure_3D ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::computePressure_2D(){
    cout << endl << "      AGCM: computePressure_2D" << endl;
    for(int k = 0; k < km; k++){
        aux_v.x[0][0][k] = c43 * aux_v.x[0][1][k] - c13 * aux_v.x[0][2][k];
        aux_v.x[0][jm-1][k] = c43 * aux_v.x[0][jm-2][k] - c13 * aux_v.x[0][jm-3][k];
        aux_w.x[0][0][k] = c43 * aux_w.x[0][1][k] - c13 * aux_w.x[0][2][k];
        aux_w.x[0][jm-1][k] = c43 * aux_w.x[0][jm-2][k] - c13 * aux_w.x[0][jm-3][k];

        rhs_v.x[0][0][k] = c43 * rhs_v.x[0][1][k] - c13 * rhs_v.x[0][2][k];
        rhs_v.x[0][jm-1][k] = c43 * rhs_v.x[0][jm-2][k] - c13 * rhs_v.x[0][jm-3][k];
        rhs_w.x[0][0][k] = c43 * rhs_w.x[0][1][k] - c13 * rhs_w.x[0][2][k];
        rhs_w.x[0][jm-1][k] = c43 * rhs_w.x[0][jm-2][k] - c13 * rhs_w.x[0][jm-3][k];
    }
    for(int j = 0; j < jm; j++){
        aux_v.x[0][j][0] = c43 * aux_v.x[0][j][1] - c13 * aux_v.x[0][j][2];
        aux_v.x[0][j][km-1] = c43 * aux_v.x[0][j][km-2] - c13 * aux_v.x[0][j][km-3];
        aux_v.x[0][j][0] = aux_v.x[0][j][km-1] = (aux_v.x[0][j][0] + aux_v.x[0][j][km-1])/2.0;
        aux_w.x[0][j][0] = c43 * aux_w.x[0][j][1] - c13 * aux_w.x[0][j][2];
        aux_w.x[0][j][km-1] = c43 * aux_w.x[0][j][km-2] - c13 * aux_w.x[0][j][km-3];
        aux_w.x[0][j][0] = aux_w.x[0][j][km-1] = (aux_w.x[0][j][0] + aux_w.x[0][j][km-1])/2.0;

        rhs_v.x[0][j][0] = c43 * rhs_v.x[0][j][1] - c13 * rhs_v.x[0][j][2];
        rhs_v.x[0][j][km-1] = c43 * rhs_v.x[0][j][km-2] - c13 * rhs_v.x[0][j][km-3];
        rhs_v.x[0][j][0] = rhs_v.x[0][j][km-1] = (rhs_v.x[0][j][0] + rhs_v.x[0][j][km-1])/2.0;
        rhs_w.x[0][j][0] = c43 * rhs_w.x[0][j][1] - c13 * rhs_w.x[0][j][2];
        rhs_w.x[0][j][km-1] = c43 * rhs_w.x[0][j][km-2] - c13 * rhs_w.x[0][j][km-3];
        rhs_w.x[0][j][0] = rhs_w.x[0][j][km-1] = (rhs_w.x[0][j][0] + rhs_w.x[0][j][km-1])/2.0;
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_land(h, 0, j, k)){
                rhs_v.x[0][j][k] = 0.0;
                rhs_w.x[0][j][k] = 0.0;
                aux_v.x[0][j][k] = 0.0;
                aux_w.x[0][j][k] = 0.0;
            }
        }
    }

    double rm = 0;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double sinthe = 0;
    double rmsinthe = 0;
    double denom = 0;
    double num2 = 0;
    double num3 = 0;
    double daux_vdthe = 0;
    double daux_wdphi = 0;
    double drhs_vdthe = 0;
    double drhs_wdphi = 0;

    rm = rad.z[0];
    for(int j = 1; j < jm-1; j++){
        sinthe = sin(the.z[j]);
        if(sinthe == 0.0) sinthe = 1.0e-5;
        rmsinthe = rm * sinthe;
        denom = 2.0/(rm * dthe2) + 2.0/(rmsinthe * dphi2);
        num2 = 1.0/(rm * dthe2);
        num3 = 1.0/(rmsinthe * dphi2);
        for(int k = 1; k < km-1; k++){
            daux_vdthe = (aux_v.x[0][j+1][k] 
                - aux_v.x[0][j-1][k])/(2.0 * dthe * rm);
            drhs_vdthe = (rhs_v.x[0][j+1][k] 
                - rhs_v.x[0][j-1][k])/(2.0 * dthe * rm);

            if((is_land(h, 0, j, k))
                &&(is_air(h, 0, j+1, k))){
                daux_vdthe = (aux_v.x[0][j+1][k] 
                    - aux_v.x[0][j][k])/(dthe * rm);
                drhs_vdthe = (rhs_v.x[0][j+1][k] 
                    - rhs_v.x[0][j][k])/(dthe * rm);
                }
            if((is_land(h, 0, j, k))
                &&(is_air(h, 0, j-1, k))){
                daux_vdthe = (aux_v.x[0][j-1][k] 
                    - aux_v.x[0][j][k])/(dthe * rm);
                drhs_vdthe = (rhs_v.x[0][j-1][k] 
                    - rhs_v.x[0][j][k])/(dthe * rm);
                }

            if((j >= 1)&&(j <= jm-2)){
                if((is_land(h, 0, j, k))
                    &&(is_air(h, 0, j+1, k)) 
                    &&(is_air(h, 0, j+2, k))){
                    daux_vdthe = (- 3.0 * aux_v.x[0][j][k] 
                        + 4.0 * aux_v.x[0][j+1][k] 
                        - aux_v.x[0][j + 2][k])/(2.0 * dthe * rm);
                    drhs_vdthe = (- 3.0 * rhs_v.x[0][j][k] 
                        + 4.0 * rhs_v.x[0][j+1][k] 
                        - rhs_v.x[0][j + 2][k])/(2.0 * dthe * rm);
                }
                if((is_land(h, 0, j, k))
                    &&(is_air(h, 0, j-1, k)) 
                    &&(is_air(h, 0, j-2, k))){
                    daux_vdthe = - (- 3.0 * aux_v.x[0][j][k] 
                        + 4.0 * aux_v.x[0][j-1][k] 
                        - aux_v.x[0][j-2][k])/(2.0 * dthe * rm);
                    drhs_vdthe = - (- 3.0 * rhs_v.x[0][j][k] 
                        + 4.0 * rhs_v.x[0][j-1][k] 
                        - rhs_v.x[0][j-2][k])/(2.0 * dthe * rm);
                }
            }

            daux_wdphi =(aux_w.x[0][j][k+1] 
                - aux_w.x[0][j][k-1])/(2.0 * dphi * rmsinthe);
            drhs_wdphi =(rhs_w.x[0][j][k+1] 
                - rhs_w.x[0][j][k-1])/(2.0 * dphi * rmsinthe);

            if((is_land(h, 0, j, k))
                &&(is_air(h, 0, j, k+1))){
                daux_wdphi = (aux_w.x[0][j][k+1] 
                    - aux_w.x[0][j][k])/(dphi * rmsinthe);
                drhs_wdphi = (rhs_w.x[0][j][k+1] 
                    - rhs_w.x[0][j][k])/(dphi * rmsinthe);
                }
            if((is_land(h, 0, j, k))
                &&(is_air(h, 0, j, k-1))){
                daux_wdphi = (aux_w.x[0][j][k-1] 
                    - aux_w.x[0][j][k])/(dphi * rmsinthe);
                drhs_wdphi = - (rhs_w.x[0][j][k-1] 
                    - rhs_w.x[0][j][k])/(dphi * rmsinthe);
                }

            if((k >= 1)&&(k <= km - 2)){
                if((is_land(h, 0, j, k))
                    &&(is_air(h, 0, j, k+1)) 
                    &&(is_air(h, 0, j, k+2))){
                    daux_wdphi = (- 3.0 * aux_w.x[0][j][k] 
                    + 4.0 * aux_w.x[0][j][k+1] 
                    - aux_w.x[0][j][k + 2])/(2.0 * rmsinthe * dphi);
                    drhs_wdphi = (- 3.0 * rhs_w.x[0][j][k] 
                    + 4.0 * rhs_w.x[0][j][k+1] 
                    - rhs_w.x[0][j][k + 2])/(2.0 * rmsinthe * dphi);
                }
                if((is_land(h, 0, j, k))
                    &&(is_air(h, 0, j, k-1)) 
                    &&(is_air(h, 0, j, k-2))){
                    daux_wdphi = - (- 3.0 * aux_w.x[0][j][k] 
                        + 4.0 * aux_w.x[0][j][k-1] 
                        - aux_w.x[0][j][k-2])/(2.0 * rmsinthe * dphi);
                    drhs_wdphi = - (- 3.0 * rhs_w.x[0][j][k] 
                        + 4.0 * rhs_w.x[0][j][k-1] 
                        - rhs_w.x[0][j][k-2])/(2.0 * rmsinthe * dphi);
                    }
            }
            if(((j <= 10)||(j <= jm-11))||((k <= 10)||(k <= km-11))){ // prevents oscillations along pole regions
                p_stat.x[0][j][k] =
                    ((p_stat.x[0][j+1][k] + p_stat.x[0][j-1][k]) * num2 
                    + (p_stat.x[0][j][k+1] + p_stat.x[0][j][k-1]) * num3 
                    )/denom;
            }else{
                p_stat.x[0][j][k] =
                    ((p_stat.x[0][j+1][k] + p_stat.x[0][j-1][k]) * num2 
                    + (p_stat.x[0][j][k+1] + p_stat.x[0][j][k-1]) * num3 
                    - (daux_vdthe - drhs_vdthe) 
                    - (daux_wdphi - drhs_wdphi)
                    )/denom;
            }
/*
            p_stat.x[0][j][k] =
                ((p_stat.x[0][j+1][k] + p_stat.x[0][j-1][k]) * num2 
                + (p_stat.x[0][j][k+1] + p_stat.x[0][j][k-1]) * num3 
                - (daux_vdthe - drhs_vdthe) 
                - (daux_wdphi - drhs_wdphi))/denom;
*/
            if(is_land(h, 0, j, k))  p_stat.x[0][j][k] = 0.0;
                if(p_stat.x[0][j][k] < 0.0) p_stat.x[0][j][k] = 0.0;


/*
    cout.precision(12);
//    if((j == 75) &&(k == 180)) cout << "northern hemisphere" << endl
    if((j == 1) &&(k == 180)) cout << "northern hemisphere" << endl
        << "   0 = " << 0 << "   j = " << j << "   k = " << k  << endl
        << "   dt = " << dt 
        << "   dthe = " << dthe 
        << "   dphi = " << dphi << endl
        << "   rm > rad = " << rm 
        << "   sinthe = " << sinthe
        << "   rmsinthe = " << rmsinthe << endl
        << "   num2 = " << num2
        << "   num3 = " << num3
        << "   denom = " << denom << endl
        << "   daux_vdthe = " << daux_vdthe
        << "   drhs_vdthe = " << drhs_vdthe << endl
        << "   daux_wdphi = " << daux_wdphi
        << "   drhs_wdphi = " << drhs_wdphi << endl

        << "   p_dynj-1 = " << p_stat.x[0][j-1][k]

        << "   p_dynj   = " << p_stat.x[0][j][k]

        << "   p_dynj+1 = " << p_stat.x[0][j+1][k] << endl

        << "   p_dynk-1 = " << p_stat.x[0][j][k-1]

        << "   p_dynk   = " << p_stat.x[0][j][k]

        << "   p_dynk+1 = " << p_stat.x[0][j][k+1] << endl

        << "   aux_v = " << aux_v.x[0][j][k]
        << "   aux_w = " << aux_w.x[0][j][k] << endl

        << "   rhs_v = " << rhs_v.x[0][j][k]
        << "   rhs_w = " << rhs_w.x[0][j][k] << endl << endl;
*/
            }
        }
    for(int k = 0; k < km; k++){
        p_stat.x[0][0][k] = c43 * p_stat.x[0][1][k] -
            c13 * p_stat.x[0][2][k];
        p_stat.x[0][jm-1][k] = c43 * p_stat.x[0][jm-2][k] -
            c13 * p_stat.x[0][jm-3][k];
    }
    for(int j = 0; j < jm; j++){
        p_stat.x[0][j][0] = c43 * p_stat.x[0][j][1] -
            c13 * p_stat.x[0][j][2];
        p_stat.x[0][j][km-1] = c43 * p_stat.x[0][j][km-2] -
            c13 * p_stat.x[0][j][km-3];
        p_stat.x[0][j][0] = p_stat.x[0][j][km-1] =
           (p_stat.x[0][j][0] + p_stat.x[0][j][km-1])/2.;
    }
    return;
}
