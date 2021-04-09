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
// total upwelling as sum on normal velocity component values in a virtual vertical column
    int i_half = (im-1)/2;
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
    double rm = 0.;
    double sinthe = 0.;
    double rmsinthe = 0.;
    int i_max = im-1;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            aux_grad_v.z[0] = r_salt_water.x[0][j][k] * v.x[0][j][k] * u_0;
            aux_grad_w.z[0] = r_salt_water.x[0][j][k] * w.x[0][j][k] * u_0;
            for(int i = 0; i <= i_max; i++){
                if(is_land(h, i, j, k)){
                    aux_grad_v.z[i] = 0.;
                    aux_grad_w.z[i] = 0.;
                }else{
                    aux_grad_v.z[i] = r_salt_water.x[i][j][k] * v.x[i][j][k] * u_0;
                    aux_grad_w.z[i] = r_salt_water.x[i][j][k] * w.x[i][j][k] * u_0;
                }
/*
    if((j == 75) &&(k == 180)) cout << endl << "north result preparation I" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k << endl
        << "   u_0 = " << u_0
        << "   r_salt_water = " << r_salt_water.x[i][j][k] << endl
        << "   aux_grad_v = " << aux_grad_v.z[i]
        << "   aux_grad_w = " << aux_grad_w.z[i] << endl;
*/
            }
/*
            if((i_max-0) % 2 == 0){
                aux_v.x[i_max][j][k] = simpson(0, i_max, dr, aux_grad_v);
                aux_w.x[i_max][j][k] = simpson(0, i_max, dr, aux_grad_w);
            }else  cout << "   i_max-0  must be an even number to use the Simpson integration method" << endl;
*/
            double dr_rm = dr * L_hyd/(double)(im-1);
            aux_v.x[i_max][j][k] = trapezoidal(0, i_max, dr_rm, aux_grad_v);
            aux_w.x[i_max][j][k] = trapezoidal(0, i_max, dr_rm, aux_grad_w);

//            aux_v.x[i_max][j][k] = trapezoidal(0, i_max, dr, aux_grad_v);
//            aux_w.x[i_max][j][k] = trapezoidal(0, i_max, dr, aux_grad_w);
//                aux_v.x[i_max][j][k] = rectangular(0, i_max, dr, aux_grad_v);
//                aux_w.x[i_max][j][k] = rectangular(0, i_max, dr, aux_grad_w);
            if(is_land(h, i_max, j, k)){
                aux_v.x[i_max][j][k] = 0.;
                aux_w.x[i_max][j][k] = 0.;
            }
/*
    if((j == 75) &&(k == 180)) cout << endl << "north result preparation II" << endl
        << "   i_max = " << i_max << "   j = " << j << "   k = " << k << endl
        << "   dr_rm = " << dr_rm << endl
        << "   dr = " << dr << endl
        << "   aux_v = " << aux_v.x[i_max][j][k]
        << "   aux_w = " << aux_w.x[i_max][j][k] << endl;
*/
        }
    }
    double coeff_pumping = 1./(365 * 24 * 60 * 60); // = 3.171e-8 produces EkmanPumping in m/a
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            if(j <= 90)
                sinthe = sin(the.z[j]);
            if(j > 90){
                int j_rev = 90 - fabs(j - 90);
                sinthe = sin(the.z[j_rev]);
            }
            if(sinthe == 0.)  sinthe = 1.e-5;
            rm = rad.z[i_max] * L_hyd;
            rmsinthe = rm * sin(the.z[j]);
            EkmanPumping.y[j][k] = - coeff_pumping * 
                ((aux_v.x[i_max][j+1][k] - aux_v.x[i_max][j-1][k])
                /(2. * rm * dthe) 
                + (aux_w.x[i_max][j][k+1] - aux_w.x[i_max][j][k-1])
                /(2. * rmsinthe * dphi));
            if(is_land(h, i_max, j, k)){
                EkmanPumping.y[j][k] = 0.;
            }
            if(EkmanPumping.y[j][k] > 0.)  Upwelling.y[j][k] = EkmanPumping.y[j][k];
            else  Upwelling.y[j][k] = 0.;
            if(EkmanPumping.y[j][k] < 0.)  Downwelling.y[j][k] = EkmanPumping.y[j][k];
            else  Downwelling.y[j][k] = 0.;
/*
    if((j == 75) &&(k == 180)) cout << endl << "north result pumping" << endl
        << "   j = " << j << "   k = " << k  << endl
        << "   aux_grad_v = " << aux_grad_v.z[i_max]
        << "   aux_grad_w = " << aux_grad_w.z[i_max] << endl
        << "   aux_v = " << aux_v.x[i_max][j][k] 
        << "   aux_w = " << aux_w.x[i_max][j][k] << endl
        << "   r_0_water = " << r_0_water << endl
        << "   EkmanPumping = " << EkmanPumping.y[j][k] << endl
        << "   u = " << u.x[i_max][j][k] 
        << "   v = " << v.x[i_max][j][k] 
        << "   w = " << w.x[i_max][j][k] << endl;
*/
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
                    BuoyancyForce_2D.y[j][k] += BuoyancyForce.x[i][j][k];
                    Salt_total.y[j][k] += c.x[i][j][k];
                }
            }
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
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
/*
            BuoyancyForce.x[0][j][k] = c43 * BuoyancyForce.x[1][j][k] -
                c13 * BuoyancyForce.x[2][j][k];
            BuoyancyForce.x[im-1][j][k] = c43 * BuoyancyForce.x[im-2][j][k] -
                c13 * BuoyancyForce.x[im-3][j][k];
            CoriolisForce.x[0][j][k] = c43 * CoriolisForce.x[1][j][k] 
                - c13 * CoriolisForce.x[2][j][k];
            CoriolisForce.x[im-1][j][k] = c43 * CoriolisForce.x[im-2][j][k] 
                - c13 * CoriolisForce.x[im-3][j][k];
            PressureGradientForce.x[0][j][k] = c43 * PressureGradientForce.x[1][j][k] 
                - c13 * PressureGradientForce.x[2][j][k];
            PressureGradientForce.x[im-1][j][k] = c43 * PressureGradientForce.x[im-2][j][k] 
                - c13 * PressureGradientForce.x[im-3][j][k];
*/
            BuoyancyForce.x[0][j][k] = BuoyancyForce.x[3][j][k] 
                - 3. * BuoyancyForce.x[2][j][k] + 3. * BuoyancyForce.x[1][j][k];  // extrapolation
            BuoyancyForce.x[im-1][j][k] = BuoyancyForce.x[im-4][j][k] 
                - 3. * BuoyancyForce.x[im-3][j][k] + 3. * BuoyancyForce.x[im-2][j][k];  // extrapolation

            CoriolisForce.x[0][j][k] = CoriolisForce.x[3][j][k] 
                - 3. * CoriolisForce.x[2][j][k] + 3. * CoriolisForce.x[1][j][k];  // extrapolation
            CoriolisForce.x[im-1][j][k] = CoriolisForce.x[im-4][j][k] 
                - 3. * CoriolisForce.x[im-3][j][k] + 3. * CoriolisForce.x[im-2][j][k];  // extrapolation

            PressureGradientForce.x[0][j][k] = PressureGradientForce.x[3][j][k] 
                - 3. * PressureGradientForce.x[2][j][k] + 3. * PressureGradientForce.x[1][j][k];  // extrapolation
            PressureGradientForce.x[im-1][j][k] = PressureGradientForce.x[im-4][j][k] 
                - 3. * PressureGradientForce.x[im-3][j][k] + 3. * PressureGradientForce.x[im-2][j][k];  // extrapolation
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
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
            BuoyancyForce.x[i][0][k] = c43 * BuoyancyForce.x[i][1][k] -
                c13 * BuoyancyForce.x[i][2][k];
            BuoyancyForce.x[i][jm-1][k] = c43 * BuoyancyForce.x[i][jm-2][k] -
                c13 * BuoyancyForce.x[i][jm-3][k];
            CoriolisForce.x[i][0][k] = c43 * CoriolisForce.x[i][1][k] 
                - c13 * CoriolisForce.x[i][2][k];
            CoriolisForce.x[i][jm-1][k] = c43 * CoriolisForce.x[i][jm-2][k]       
                - c13 * CoriolisForce.x[i][jm-3][k];
            PressureGradientForce.x[i][0][k] = c43 * PressureGradientForce.x[i][1][k] 
                - c13 * PressureGradientForce.x[i][2][k];
            PressureGradientForce.x[i][jm-1][k] = c43 * PressureGradientForce.x[i][jm-2][k]       
                - c13 * PressureGradientForce.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
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
            BuoyancyForce.x[i][j][0] = c43 * BuoyancyForce.x[i][j][1] -
                c13 * BuoyancyForce.x[i][j][2];
            BuoyancyForce.x[i][j][km-1] = c43 * BuoyancyForce.x[i][j][km-2] -
                c13 * BuoyancyForce.x[i][j][km-3];
            BuoyancyForce.x[i][j][0] = BuoyancyForce.x[i][j][km-1] =
               (BuoyancyForce.x[i][j][0] + BuoyancyForce.x[i][j][km-1])/2.;
            CoriolisForce.x[i][j][0] = c43 * CoriolisForce.x[i][j][1] 
                - c13 * CoriolisForce.x[i][j][2];
            CoriolisForce.x[i][j][km-1] = c43 * CoriolisForce.x[i][j][km-2] 
                - c13 * CoriolisForce.x[i][j][km-3];
            CoriolisForce.x[i][j][0] = CoriolisForce.x[i][j][km-1] =
               (CoriolisForce.x[i][j][0] + CoriolisForce.x[i][j][km-1])/ 2.;
            PressureGradientForce.x[i][j][0] = c43 * PressureGradientForce.x[i][j][1] 
                - c13 * PressureGradientForce.x[i][j][2];
            PressureGradientForce.x[i][j][km-1] = c43 * PressureGradientForce.x[i][j][km-2] 
                - c13 * PressureGradientForce.x[i][j][km-3];
            PressureGradientForce.x[i][j][0] = PressureGradientForce.x[i][j][km-1] =
               (PressureGradientForce.x[i][j][0] + PressureGradientForce.x[i][j][km-1])/ 2.;
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


