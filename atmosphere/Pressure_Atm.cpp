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
    double coeff_p = p_0/(r_air*u_0*u_0);
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
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
            aux_u.x[i][j][0] = aux_u.x[i][j][km-1] = (aux_u.x[i][j][0] + aux_u.x[i][j][km-1])/2.;
            aux_v.x[i][j][0] = c43 * aux_v.x[i][j][1] - c13 * aux_v.x[i][j][2];
            aux_v.x[i][j][km-1] = c43 * aux_v.x[i][j][km-2] - c13 * aux_v.x[i][j][km-3];
            aux_v.x[i][j][0] = aux_v.x[i][j][km-1] = (aux_v.x[i][j][0] + aux_v.x[i][j][km-1])/2.;
            aux_w.x[i][j][0] = c43 * aux_w.x[i][j][1] - c13 * aux_w.x[i][j][2];
            aux_w.x[i][j][km-1] = c43 * aux_w.x[i][j][km-2] - c13 * aux_w.x[i][j][km-3];
            aux_w.x[i][j][0] = aux_w.x[i][j][km-1] = (aux_w.x[i][j][0] + aux_w.x[i][j][km-1])/2.;

            rhs_u.x[i][j][0] = c43 * rhs_u.x[i][j][1] - c13 * rhs_u.x[i][j][2];
            rhs_u.x[i][j][km-1] = c43 * rhs_u.x[i][j][km-2] - c13 * rhs_u.x[i][j][km-3];
            rhs_u.x[i][j][0] = rhs_u.x[i][j][km-1] = (rhs_u.x[i][j][0] + rhs_u.x[i][j][km-1])/2.;
            rhs_v.x[i][j][0] = c43 * rhs_v.x[i][j][1] - c13 * rhs_v.x[i][j][2];
            rhs_v.x[i][j][km-1] = c43 * rhs_v.x[i][j][km-2] - c13 * rhs_v.x[i][j][km-3];
            rhs_v.x[i][j][0] = rhs_v.x[i][j][km-1] = (rhs_v.x[i][j][0] + rhs_v.x[i][j][km-1])/2.;
            rhs_w.x[i][j][0] = c43 * rhs_w.x[i][j][1] - c13 * rhs_w.x[i][j][2];
            rhs_w.x[i][j][km-1] = c43 * rhs_w.x[i][j][km-2] - c13 * rhs_w.x[i][j][km-3];
            rhs_w.x[i][j][0] = rhs_w.x[i][j][km-1] = (rhs_w.x[i][j][0] + rhs_w.x[i][j][km-1])/2.;
        }
    }
    double rm = 0;
    double rm2 = 0;
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double sinthe = 0;
    double rmsinthe = 0;
    double rm2sinthe2 = 0;
    double denom = 0;
    double num1 = 0;
    double num2 = 0;
    double num3 = 0;
    double daux_udr = 0;
    double daux_vdthe = 0;
    double daux_wdphi = 0;
    double drhs_udr = 0;
    double drhs_vdthe = 0;
    double drhs_wdphi = 0;
    for(int i = 1; i < im-1; i++){
        rm = rad.z[i];
        rm2 = rm * rm;
        double zeta = 3.715;
        double exp_rm = 1./exp(zeta * rm);
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            rmsinthe = rm * sinthe;
            rm2sinthe2 = rmsinthe * rmsinthe;
            denom = 2./dr2 + 2./(rm2 * dthe2) + 2./(rm2sinthe2 * dphi2);
            num1 = 1./dr2;
            num2 = 1./(rm2 * dthe2);
            num3 = 1./(rm2sinthe2 * dphi2);
            for(int k = 1; k < km-1; k++){
                daux_udr = (aux_u.x[i+1][j][k] - aux_u.x[i-1][j][k])/(2. * dr) * exp_rm;
                if(i <= im - 3){
                    if((is_land(h, i, j, k)) && (is_water(h, i+1, j, k))){        
                        daux_udr = (- 3. * aux_u.x[i][j][k] 
                            + 4. * aux_u.x[i+1][j][k] 
                            - aux_u.x[i + 2][j][k])/(2. * dr) * exp_rm;
                        drhs_udr = (- 3. * rhs_u.x[i][j][k] 
                            + 4. * rhs_u.x[i+1][j][k] 
                            - rhs_u.x[i + 2][j][k])/(2. * dr) * exp_rm;
                        }
                }else{
                    daux_udr = (aux_u.x[i+1][j][k] - aux_u.x[i][j][k])/dr * exp_rm;
                    drhs_udr = (rhs_u.x[i+1][j][k] - rhs_u.x[i][j][k])/dr * exp_rm;
                }
                daux_vdthe = (aux_v.x[im-1][j+1][k] - aux_v.x[im-1][j-1][k])
                    /(2. * dthe * rm);
                drhs_vdthe = (rhs_v.x[im-1][j+1][k] - rhs_v.x[im-1][j-1][k])
                    /(2. * dthe * rm);
                if((is_land(h, i, j, k)) && (is_water(h, i, j+1, k))){
                    daux_vdthe = (aux_v.x[im-1][j+1][k] 
                        - aux_v.x[im-1][j][k])/(dthe * rm);
                    drhs_vdthe = (rhs_v.x[im-1][j+1][k] 
                        - rhs_v.x[im-1][j][k])/(dthe * rm);
                    }
                if((is_land(h, i, j, k)) && (is_water(h, i, j-1, k))){
                    daux_vdthe = (aux_v.x[im-1][j-1][k] 
                        - aux_v.x[im-1][j][k])/(dthe * rm);
                    drhs_vdthe = (rhs_v.x[im-1][j-1][k] 
                        - rhs_v.x[im-1][j][k])/(dthe * rm);
                    }
                if((j >= 2) &&(j <= jm - 3)){
                    if((is_land(h, i, j, k)) && (is_water(h, i, j+1, k)) 
                        &&(is_water(h, i, j+2, k))){
                        daux_vdthe = (- 3. * aux_v.x[im-1][j][k] 
                            + 4. * aux_v.x[im-1][j+1][k] 
                            - aux_v.x[im-1][j + 2][k])/(2. * dthe * rm);
                        drhs_vdthe = (- 3. * rhs_v.x[im-1][j][k] 
                            + 4. * rhs_v.x[im-1][j+1][k] 
                            - rhs_v.x[im-1][j + 2][k])/(2. * dthe * rm);
                        }
                    if((is_land(h, i, j, k)) && (is_water(h, i, j-1, k)) 
                        &&(is_water(h, i, j-2, k))){
                        daux_vdthe = (- 3. * aux_v.x[im-1][j][k] 
                            + 4. * aux_v.x[im-1][j-1][k] 
                            - aux_v.x[im-1][j-2][k])/(2. * dthe * rm);
                        drhs_vdthe = (- 3. * rhs_v.x[im-1][j][k] 
                            + 4. * rhs_v.x[im-1][j-1][k] 
                            - rhs_v.x[im-1][j-2][k])/(2. * dthe * rm);
                        }
                }
                daux_wdphi =(aux_w.x[im-1][j][k+1] 
                    - aux_w.x[im-1][j][k-1])/(2. * dphi * rmsinthe);
                drhs_wdphi =(rhs_w.x[im-1][j][k+1] 
                    - rhs_w.x[im-1][j][k-1])/(2. * dphi * rmsinthe);
                if((is_land(h, i, j, k)) && (is_water(h, i, j, k+1))){
                    daux_wdphi = (aux_w.x[im-1][j][k+1] 
                        - aux_w.x[im-1][j][k])/(dphi * rmsinthe);
                    drhs_wdphi = (rhs_w.x[im-1][j][k+1] 
                        - rhs_w.x[im-1][j][k])/(dphi * rmsinthe);
                    }
                if((is_land(h, i, j, k)) && (is_water(h, i, j, k-1))){
                    daux_wdphi = (aux_w.x[im-1][j][k-1] 
                        - aux_w.x[im-1][j][k])/(dphi * rmsinthe);
                    drhs_wdphi = (rhs_w.x[im-1][j][k-1] 
                        - rhs_w.x[im-1][j][k])/(dphi * rmsinthe);
                    }
                if((k >= 2) &&(k <= km - 3)){
                    if((is_land(h, i, j, k)) && (is_water(h, i, j, k+1)) 
                        &&(is_water(h, i, j, k+2))){
                        daux_wdphi = (- 3. * aux_w.x[im-1][j][k] 
                        + 4. * aux_w.x[im-1][j][k+1] 
                        - aux_w.x[im-1][j][k + 2])/(2. * rmsinthe * dphi);
                        drhs_wdphi = (- 3. * rhs_w.x[im-1][j][k] 
                        + 4. * rhs_w.x[im-1][j][k+1] 
                        - rhs_w.x[im-1][j][k + 2])/(2. * rmsinthe * dphi);
                    }
                    if((is_land(h, i, j, k)) && (is_water(h, i, j, k-1)) 
                        &&(is_water(h, i, j, k-2))){
                        daux_wdphi = (- 3. * aux_w.x[im-1][j][k] 
                            + 4. * aux_w.x[im-1][j][k - 1] 
                            - aux_w.x[im-1][j][k-2])/(2. * rmsinthe * dphi);
                        drhs_wdphi = (- 3. * rhs_w.x[im-1][j][k] 
                            + 4. * rhs_w.x[im-1][j][k - 1] 
                            - rhs_w.x[im-1][j][k-2])/(2. * rmsinthe * dphi);
                        }
                }
                p_dyn.x[i][j][k] = ((p_dyn.x[i+1][j][k] 
                    + p_dyn.x[i-1][j][k]) * num1 + (p_dyn.x[i][j+1][k] 
                    + p_dyn.x[i][j-1][k]) * num2 + (p_dyn.x[i][j][k+1] 
                    + p_dyn.x[i][j][k-1]) * num3 
                    - ((daux_udr - drhs_udr)/coeff_p 
                    + (daux_vdthe - drhs_vdthe) * rm/coeff_p 
                    + (daux_wdphi - drhs_wdphi) * rmsinthe/coeff_p))/denom;
            }
        }
    }
}
/*
*
*/
void cAtmosphereModel::computePressure_2D(){
    double coeff_p = p_0/(r_air*u_0*u_0);
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
        aux_v.x[0][j][0] = aux_v.x[0][j][km-1] = (aux_v.x[0][j][0] + aux_v.x[0][j][km-1])/2.;
        aux_w.x[0][j][0] = c43 * aux_w.x[0][j][1] - c13 * aux_w.x[0][j][2];
        aux_w.x[0][j][km-1] = c43 * aux_w.x[0][j][km-2] - c13 * aux_w.x[0][j][km-3];
        aux_w.x[0][j][0] = aux_w.x[0][j][km-1] = (aux_w.x[0][j][0] + aux_w.x[0][j][km-1])/2.;

        rhs_v.x[0][j][0] = c43 * rhs_v.x[0][j][1] - c13 * rhs_v.x[0][j][2];
        rhs_v.x[0][j][km-1] = c43 * rhs_v.x[0][j][km-2] - c13 * rhs_v.x[0][j][km-3];
        rhs_v.x[0][j][0] = rhs_v.x[0][j][km-1] = (rhs_v.x[0][j][0] + rhs_v.x[0][j][km-1])/2.;
        rhs_w.x[0][j][0] = c43 * rhs_w.x[0][j][1] - c13 * rhs_w.x[0][j][2];
        rhs_w.x[0][j][km-1] = c43 * rhs_w.x[0][j][km-2] - c13 * rhs_w.x[0][j][km-3];
        rhs_w.x[0][j][0] = rhs_w.x[0][j][km-1] = (rhs_w.x[0][j][0] + rhs_w.x[0][j][km-1])/2.;
    }
    double rm= 0;
    double rm2= 0;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double sinthe = 0;
    double rmsinthe = 0;
    double rm2sinthe2 = 0;
    double denom = 0;
    double num2 = 0;
    double num3 = 0;
    double daux_vdthe = 0;
    double daux_wdphi = 0;
    double drhs_vdthe = 0;
    double drhs_wdphi = 0;
    rm = rad.z[0];
    rm2 = rm * rm;
    for(int j = 1; j < jm-1; j++){
        sinthe = sin(the.z[j]);
        rmsinthe = rm * sinthe;
        rm2sinthe2 = rmsinthe * rmsinthe;
        denom = 2./(rm2 * dthe2) + 2./(rm2sinthe2 * dphi2);
        num2 = 1./(rm2 * dthe2);
        num3 = 1./(rm2sinthe2 * dphi2);
        for(int k = 1; k < km-1; k++){
            daux_vdthe = (aux_v.x[0][j+1][k] - aux_v.x[0][j-1][k])
                /(2. * dthe * rm);
            drhs_vdthe = (rhs_v.x[0][j+1][k] - rhs_v.x[0][j-1][k])
                /(2. * dthe * rm);
            if(is_land(h, 0, j, k) && is_water(h, 0, j+1, k)){
                if((j < jm-2) && is_water(h, 0, j+2, k)){
                    daux_vdthe = (- 3. * aux_v.x[0][j][k] 
                        + 4. * aux_v.x[0][j+1][k] - aux_v.x[0][j+2][k])
                        /(2. * dthe * rm);
                    drhs_vdthe = (- 3. * rhs_v.x[0][j][k] 
                        + 4. * rhs_v.x[0][j+1][k] - rhs_v.x[0][j+2][k])
                        /(2. * dthe * rm);
                }else{
                    daux_vdthe = (aux_v.x[0][j+1][k] 
                        - aux_v.x[0][j][k])/(dthe * rm);
                    drhs_vdthe = (rhs_v.x[0][j+1][k] 
                        - rhs_v.x[0][j][k])/(dthe * rm);
                }
            }
            if(is_land(h, 0, j, k) && is_water(h, 0, j-1, k)){
                if((j > 1)  && is_water(h, 0, j-2, k)){
                    daux_vdthe = (- 3. * aux_v.x[0][j][k] 
                        + 4. * aux_v.x[0][j-1][k] - aux_v.x[0][j-2][k])
                        /(2. * dthe * rm);
                    drhs_vdthe = (- 3. * rhs_v.x[0][j][k] 
                        + 4. * rhs_v.x[0][j-1][k] - rhs_v.x[0][j-2][k])
                        /(2. * dthe * rm);
                }else{
                    daux_vdthe = (aux_v.x[0][j-1][k] 
                        - aux_v.x[0][j][k])/(dthe * rm);
                    drhs_vdthe = (rhs_v.x[0][j-1][k] 
                        - rhs_v.x[0][j][k])/(dthe * rm);
                }
            }
            daux_wdphi = (aux_w.x[0][j][k+1] - aux_w.x[0][j][k-1])
                /(2. * dphi * rmsinthe);
            drhs_wdphi = (rhs_w.x[0][j][k+1] - rhs_w.x[0][j][k-1])
                /(2. * dphi * rmsinthe);
            if(is_land(h, 0, j, k) && is_water(h, 0, j, k+1)){
                if(k < km-2 && is_water(h, 0, j, k+2)){
                    daux_wdphi = (- 3. * aux_w.x[0][j][k] 
                        + 4. * aux_w.x[0][j][k+1] - aux_w.x[0][j][k+2])
                        /(2. * rmsinthe * dphi);
                    drhs_wdphi = (- 3. * rhs_w.x[0][j][k] 
                        + 4. * rhs_w.x[0][j][k+1] - rhs_w.x[0][j][k+2])
                        /(2. * rmsinthe * dphi);
                }else{
                    daux_wdphi = (aux_w.x[0][j][k+1] 
                        - aux_w.x[0][j][k])/(dphi * rmsinthe);
                    drhs_wdphi = (rhs_w.x[0][j][k+1] 
                        - rhs_w.x[0][j][k])/(dphi * rmsinthe);
                }
            }
            if(is_land(h, 0, j, k) && is_water(h, 0, j, k-1)){
                if(k >= 2 && is_water(h, 0, j, k-2)){
                    daux_wdphi = (- 3. * aux_w.x[0][j][k] 
                        + 4. * aux_w.x[0][j][k-1] - aux_w.x[0][j][k-2])
                        /(2. * rmsinthe * dphi);
                    drhs_wdphi = (- 3. * rhs_w.x[0][j][k] 
                        + 4. * rhs_w.x[0][j][k-1] - rhs_w.x[0][j][k-2])
                        /(2. * rmsinthe * dphi);
                }else{
                    daux_wdphi = (aux_w.x[0][j][k-1]  
                        - aux_w.x[0][j][k])/(dphi * rmsinthe);
                    drhs_wdphi = (rhs_w.x[0][j][k-1]  
                        - rhs_w.x[0][j][k])/(dphi * rmsinthe);
                }
            }
            p_dyn.x[0][j][k] = ((p_dyn.x[0][j+1][k] 
                + p_dyn.x[0][j-1][k]) * num2 + (p_dyn.x[0][j][k+1] 
                + p_dyn.x[0][j][k-1]) * num3 
                - ((daux_vdthe  - drhs_vdthe) * rm/coeff_p 
                + (daux_wdphi - drhs_wdphi) * rmsinthe/coeff_p))/denom;
        }
    }
}
