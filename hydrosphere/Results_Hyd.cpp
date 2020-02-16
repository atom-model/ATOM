/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "cHydrosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

void cHydrosphereModel::run_data_hyd(){
    float c43 = 4./3.;
    float c13 = 1./3.;
// total upwelling as sum on normal velocity component values in a virtual vertical column
    int i_half =(im-1)/2;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Upwelling.y[j][k] = 0.;  // Upwelling
            Downwelling.y[j][k] = 0.;  // Downwelling
            EkmanPumping.y[j][k] = 0.;  // Ekman pumping
            SaltFinger.y[j][k] = 0.;  // SaltFinger
            SaltDiffusion.y[j][k] = 0.;  // SaltFinger
            BuoyancyForce_2D.y[j][k] = 0.; // Saltdiffusion
            Salt_total.y[j][k] = 0.;  // total Salt
        }
    }
// Ekman layer computation, variable EkmanPumping is equal to the vertical velocity u at the surface
    double rmsinthe = 0.;
    int i_max = im-1;
    int i_beg = 0;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
        double sinthe = sin(fabs((90. -(double)j) * M_PI/180.));
        if(sinthe == 0.)  sinthe = 1.e-5;
            for(int i = 1; i < im; i++){
                double vel_magnit = sqrt(v.x[im-1][j][k] * v.x[im-1][j][k] 
                    + w.x[im-1][j][k] * w.x[im-1][j][k]); // in m/s
                // dimensional surface wind velocity U_10 in m/s
                // original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 139, eq. 9.16
                // assumed depth at i_Ekman = 3 x pi/gam, velocity opposite direction plus pi to surface velocity
                double i_Ekman_layer = 3. * 7.6/sqrt(sinthe) * vel_magnit; // dimensional surface wind velocity U_10 in m/s
                i_Ekman_layer = 1.5 * i_Ekman_layer;  // factor gives an Ekman depth of 200m
                double coeff = i_Ekman_layer/L_hyd;
                if(coeff >= 1.) coeff = .623;
                int i_Ekman = (im-1) * (1.-coeff);
                if(is_land(h, i-1, j, k) && is_land(h, i, j, k)) i_beg = i;
                else  i_beg = i_Ekman;
            } 
            for(int i = 0; i < im; i++){
                if(is_land(h, i, j, k)){
                    aux_grad_v.z[i] = 0.;
                    aux_grad_w.z[i] = 0.;
                }else{
                    i_beg = i;
                    aux_grad_v.z[i] = v.x[i][j][k];
                    aux_grad_w.z[i] = w.x[i][j][k];
                }
            }
/*
            if((i_max-i_Ekman) % 2 == 0){
                aux_v.x[im-1][j][k] = simpson(i_beg, i_max, dr, aux_grad_v);
                aux_w.x[im-1][j][k] = simpson(i_beg, i_max, dr, aux_grad_w);
            }
            else  cout << "   i_max-i_Ekman  must be an even number to use the Simpson integration method" << endl;
*/
            aux_v.x[im-1][j][k] = trapezoidal(i_beg, i_max, dr, aux_grad_v);
            aux_w.x[im-1][j][k] = trapezoidal(i_beg, i_max, dr, aux_grad_w);
//                aux_v.x[im-1][j][k] = rectangular(i_beg, i_max, dr, aux_grad_v);
//                aux_w.x[im-1][j][k] = rectangular(i_beg, i_max, dr, aux_grad_w);
            if(is_land(h, im-1, j, k)){
                aux_v.x[im-1][j][k] = 0.;
                aux_w.x[im-1][j][k] = 0.;
            }
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            rmsinthe = rad.z[im-1] * sin(the.z[j]);
            EkmanPumping.y[j][k] = -((aux_v.x[im-1][j+1][k] 
                - aux_v.x[im-1][j-1][k])/(2. * rad.z[im-1] * dthe) 
                + (aux_w.x[im-1][j][k+1] - aux_w.x[im-1][j][k-1])
                /(2. * rmsinthe * dphi)) * u_0;
            if(is_land(h, im-1, j, k)){
                EkmanPumping.y[j][k] = 0.;
            }
            if(EkmanPumping.y[j][k] > 0.)  Upwelling.y[j][k] = EkmanPumping.y[j][k];
            else  Upwelling.y[j][k] = 0.;
            if(EkmanPumping.y[j][k] < 0.)  Downwelling.y[j][k] = EkmanPumping.y[j][k];
            else  Downwelling.y[j][k] = 0.;
        }
    }
    for(int k = 0; k < km; k++){
        EkmanPumping.y[0][k] = c43 * EkmanPumping.y[1][k] -
            c13 * EkmanPumping.y[2][k];
        EkmanPumping.y[jm-1][k] = c43 * EkmanPumping.y[jm-2][k] -
            c13 * EkmanPumping.y[jm-3][k];
        Upwelling.y[0][k] = c43 * Upwelling.y[1][k] -
            c13 * Upwelling.y[2][k];
        Upwelling.y[jm-1][k] = c43 * Upwelling.y[jm-2][k] -
            c13 * Upwelling.y[jm-3][k];
        Downwelling.y[0][k] = c43 * Downwelling.y[1][k] -
            c13 * Downwelling.y[2][k];
        Downwelling.y[jm-1][k] = c43 * Downwelling.y[jm-2][k] -
            c13 * Downwelling.y[jm-3][k];
    }
    for(int j = 0; j < jm; j++){
        EkmanPumping.y[j][0] = c43 * EkmanPumping.y[j][1] - c13 * EkmanPumping.y[j][2];
        EkmanPumping.y[j][km-1] = c43 * EkmanPumping.y[j][km-2] -
            c13 * EkmanPumping.y[j][km-3];
        EkmanPumping.y[j][0] = EkmanPumping.y[j][km-1] =(EkmanPumping.y[j][0] +
            EkmanPumping.y[j][km-1])/2.;
        Upwelling.y[j][0] = c43 * Upwelling.y[j][1] - c13 * Upwelling.y[j][2];
        Upwelling.y[j][km-1] = c43 * Upwelling.y[j][km-2] -
            c13 * Upwelling.y[j][km-3];
        Upwelling.y[j][0] = Upwelling.y[j][km-1] =(Upwelling.y[j][0] +
            Upwelling.y[j][km-1])/2.;
        Downwelling.y[j][0] = c43 * Downwelling.y[j][1] - c13 * Downwelling.y[j][2];
        Downwelling.y[j][km-1] = c43 * Downwelling.y[j][km-2] -
            c13 * Downwelling.y[j][km-3];
        Downwelling.y[j][0] = Downwelling.y[j][km-1] =(Downwelling.y[j][0] +
            Downwelling.y[j][km-1])/2.;
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = i_half; i < im; i++){
                if(is_water( h, i, j, k)){
                    SaltFinger.y[j][k] += Salt_Finger.x[i][j][k];
                    SaltDiffusion.y[j][k] += Salt_Diffusion.x[i][j][k];
                    BuoyancyForce_2D.y[j][k] += BuoyancyForce_3D.x[i][j][k];
                    Salt_total.y[j][k] += c.x[i][j][k];
                }
            }
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            BuoyancyForce_3D.x[0][j][k] = c43 * BuoyancyForce_3D.x[1][j][k] -
                c13 * BuoyancyForce_3D.x[2][j][k];
            BuoyancyForce_3D.x[im-1][j][k] = c43 * BuoyancyForce_3D.x[im-2][j][k] -
                c13 * BuoyancyForce_3D.x[im-3][j][k];
            Salt_Finger.x[0][j][k] = c43 * Salt_Finger.x[1][j][k] -
                c13 * Salt_Finger.x[2][j][k];
            Salt_Finger.x[im-1][j][k] = c43 * Salt_Finger.x[im-2][j][k] -
                c13 * Salt_Finger.x[im-3][j][k];
            Salt_Diffusion.x[0][j][k] = c43 * Salt_Diffusion.x[1][j][k] -
                c13 * Salt_Diffusion.x[2][j][k];
            Salt_Diffusion.x[im-1][j][k] = c43 * Salt_Diffusion.x[im-2][j][k] -
                c13 * Salt_Diffusion.x[im-3][j][k];
            Salt_Balance.x[0][j][k] = c43 * Salt_Balance.x[1][j][k] -
                c13 * Salt_Balance.x[2][j][k];
            Salt_Balance.x[im-1][j][k] = c43 * Salt_Balance.x[im-2][j][k] -
                c13 * Salt_Balance.x[im-3][j][k];
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            BuoyancyForce_3D.x[i][0][k] = c43 * BuoyancyForce_3D.x[i][1][k] -
                c13 * BuoyancyForce_3D.x[i][2][k];
            BuoyancyForce_3D.x[i][jm-1][k] = c43 * BuoyancyForce_3D.x[i][jm-2][k] -
                c13 * BuoyancyForce_3D.x[i][jm-3][k];
            Salt_Finger.x[i][0][k] = c43 * Salt_Finger.x[i][1][k] -
                c13 * Salt_Finger.x[i][2][k];
            Salt_Finger.x[i][jm-1][k] = c43 * Salt_Finger.x[i][jm-2][k] -
                c13 * Salt_Finger.x[i][jm-3][k];
            Salt_Diffusion.x[i][0][k] = c43 * Salt_Diffusion.x[i][1][k] -
                 c13 * Salt_Diffusion.x[i][2][k];
            Salt_Diffusion.x[i][jm-1][k] = c43 * Salt_Diffusion.x[i][jm-2][k] -
                c13 * Salt_Diffusion.x[i][jm-3][k];
            Salt_Balance.x[i][0][k] = c43 * Salt_Balance.x[i][1][k] -
                c13 * Salt_Balance.x[i][2][k];
            Salt_Balance.x[i][jm-1][k] = c43 * Salt_Balance.x[i][jm-2][k] -
                c13 * Salt_Balance.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            BuoyancyForce_3D.x[i][j][0] = c43 * BuoyancyForce_3D.x[i][j][1] -
                c13 * BuoyancyForce_3D.x[i][j][2];
            BuoyancyForce_3D.x[i][j][km-1] = c43 * BuoyancyForce_3D.x[i][j][km-2] -
                c13 * BuoyancyForce_3D.x[i][j][km-3];
            BuoyancyForce_3D.x[i][j][0] = BuoyancyForce_3D.x[i][j][km-1] =
               (BuoyancyForce_3D.x[i][j][0] + BuoyancyForce_3D.x[i][j][km-1])/2.;
            Salt_Finger.x[i][j][0] = c43 * Salt_Finger.x[i][j][1] -
                c13 * Salt_Finger.x[i][j][2];
            Salt_Finger.x[i][j][km-1] = c43 * Salt_Finger.x[i][j][km-2] -
                c13 * Salt_Finger.x[i][j][km-3];
            Salt_Finger.x[i][j][0] = Salt_Finger.x[i][j][km-1] =(Salt_Finger.x[i][j][0] +
                 Salt_Finger.x[i][j][km-1])/2.;
            Salt_Diffusion.x[i][j][0] = c43 * Salt_Diffusion.x[i][j][1] -
                c13 * Salt_Diffusion.x[i][j][2];
            Salt_Diffusion.x[i][j][km-1] = c43 * Salt_Diffusion.x[i][j][km-2] -
                c13 * Salt_Diffusion.x[i][j][km-3];
            Salt_Diffusion.x[i][j][0] = Salt_Diffusion.x[i][j][km-1] =
               (Salt_Diffusion.x[i][j][0] + Salt_Diffusion.x[i][j][km-1])/2.;
            Salt_Balance.x[i][j][0] = c43 * Salt_Balance.x[i][j][1] -
                c13 * Salt_Balance.x[i][j][2];
            Salt_Balance.x[i][j][km-1] = c43 * Salt_Balance.x[i][j][km-2] -
                c13 * Salt_Balance.x[i][j][km-3];
            Salt_Balance.x[i][j][0] = Salt_Balance.x[i][j][km-1] =
               (Salt_Balance.x[i][j][0] + Salt_Balance.x[i][j][km-1])/2.;
        }
    }
}
/*
*
*/
void cHydrosphereModel::land_oceanFraction(){
    double h_point_max = (jm-1) * (km-1);
    double h_land = 0;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_land( h, im-1, j, k))  h_land++;
        }
    }
    double h_ocean = h_point_max - h_land;
    double ocean_land =(double)h_ocean/(double)h_land;
    cout.precision(3);
    cout << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.')
        << "      total number of points at constant hight " << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << h_point_max << endl << setiosflags(ios::left)
        << setw(50) << setfill('.') << "      number of points on the ocean surface " << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ') << h_ocean << endl
        << setiosflags(ios::left) << setw(50) << setfill('.') << "      number of points on the land surface "
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ') << h_land
        << endl << setiosflags(ios::left) << setw(50) << setfill('.') << "      ocean/land ratio "
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << ocean_land << endl << endl;
    cout << endl;
}


