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

#include "cHydrosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

#define dxdthe_a(X) \
    (h_d_j * (X->x[i][j+1][k] - X->x[i][j-1][k])/(2. * dthe))
#define dxdthe_b(X) \
    (h_d_j * (- 3. * X->x[i][j][k] + 4. * X->x[i][j+1][k] - X->x[i][j+2][k])/(2. * dthe))
#define dxdthe_c(X) \
    (- h_d_j * (- 3. * X->x[i][j][k] + 4. * X->x[i][j-1][k] - X->x[i][j-2][k])/(2. * dthe))
#define dxdphi_a(X) \
    (h_d_k * (X->x[i][j][k+1] - X->x[i][j][k-1])/(2. * dphi))
#define dxdphi_b(X) \
    (h_d_k * (- 3. * X->x[i][j][k] + 4. * X->x[i][j][k+1] - X->x[i][j][k+2])/(2. * dphi))
#define dxdphi_c(X) \
    (- h_d_k * (- 3. * X->x[i][j][k] + 4. * X->x[i][j][k-1] - X->x[i][j][k-2])/(2. * dphi))

void cHydrosphereModel::run_data_hyd(){
// total upwelling as sum on normal velocity component values in a virtual vertical column
    int i_max = im-1;
    int i_half = (im-1)/2;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Upwelling.y[j][k] = 0.;  // Upwelling
            Downwelling.y[j][k] = 0.;  // Downwelling
            EkmanPumping.y[j][k] = 0.;  // Ekman pumping
            SaltFinger.y[j][k] = 0.;  // SaltFinger
            SaltDiffusion.y[j][k] = 0.;  // SaltFinger
            Salt_total.y[j][k] = 0.;  // total Salt
        }
    }
// Ekman layer computation, variable EkmanPumping is equal to the vertical velocity u at the surface
    double rm = 0.;
    double rmsinthe = 0.;
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
            double dr_rm = dr * L_hyd/(double)i_max;
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
    for(int k = 0; k < km; k++){
        aux_v.x[i_max][0][k] = c43 * aux_v.x[i_max][1][k] - c13 * aux_v.x[i_max][2][k];
        aux_v.x[i_max][jm-1][k] = c43 * aux_v.x[i_max][jm-2][k] - c13 * aux_v.x[i_max][jm-3][k];
        aux_w.x[i_max][0][k] = c43 * aux_w.x[i_max][1][k] - c13 * aux_w.x[i_max][2][k];
        aux_w.x[i_max][jm-1][k] = c43 * aux_w.x[i_max][jm-2][k] - c13 * aux_w.x[i_max][jm-3][k];
    }
    for(int j = 0; j < jm; j++){
        aux_v.x[i_max][j][0] = c43 * aux_v.x[i_max][j][1] - c13 * aux_v.x[i_max][j][2];
        aux_v.x[i_max][j][km-1] = c43 * aux_v.x[i_max][j][km-2] - c13 * aux_v.x[i_max][j][km-3];
        aux_v.x[i_max][j][0] = aux_v.x[i_max][j][km-1] = (aux_v.x[i_max][j][0] + aux_v.x[i_max][j][km-1])/2.;
        aux_w.x[i_max][j][0] = c43 * aux_w.x[i_max][j][1] - c13 * aux_w.x[i_max][j][2];
        aux_w.x[i_max][j][km-1] = c43 * aux_w.x[i_max][j][km-2] - c13 * aux_w.x[i_max][j][km-3];
        aux_w.x[i_max][j][0] = aux_w.x[i_max][j][km-1] = (aux_w.x[i_max][j][0] + aux_w.x[i_max][j][km-1])/2.;
    }
    int i = i_max; // for the defined macros
    double coeff_pumping = 864.0; // produces EkmanPumping from m/s into cm/d
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            double sinthe = sin(the.z[j]);
            rm = L_hyd;
            rmsinthe = rm * sinthe;
            double h_d_j = 0.0;
            double h_d_k = 0.0;
            if(is_water(h, im-1, j, k)){
                h_d_j = 1.0;
                h_d_k = 1.0;
           }
           if(is_land(h, im-1, j, k)){
               h_d_j = 0.0;
               h_d_k = 0.0;
            }
            std::vector<Array*> arrays_1{&aux_v, &aux_w};
            enum array_index_1{i_v_1, i_w_1, last_array_index_1};
            std::vector<double> dxdthe_vals(last_array_index_1), 
                                dxdphi_vals(last_array_index_1);
            bool the_flag = false, phi_flag = false;
            if((j >= 2)&&(j < jm-3)){
                if((is_land(h, 0, j, k))
                    &&((is_air(h, 0, j+1, k)) 
                    &&(is_air(h, 0, j+2, k)))){
                    for(int n=0; n<last_array_index_1; n++)
                        dxdthe_vals[n] = dxdthe_b(arrays_1[n]);
                    the_flag = true;
                }
                if((is_land(h, 0, j, k))
                    &&(is_air(h, 0, j-1, k)) 
                    &&(is_air(h, 0, j-2, k))){
                    for(int n=0; n<last_array_index_1; n++)
                        dxdthe_vals[n] = dxdthe_c(arrays_1[n]);
                    the_flag = true;
                }
            }
            if((k >= 2)&&(k < km-3)){
                if((is_land(h, 0, j, k))
                    &&(is_air(h, 0, j, k+1)) 
                    &&(is_air(h, 0, j, k+2))){
                    for(int n=0; n<last_array_index_1; n++)
                        dxdphi_vals[n] = dxdphi_b(arrays_1[n]);
                    phi_flag = true;
                }
                if((is_land(h, 0, j, k))
                    &&(is_air(h, 0, j, k-1)) 
                    &&(is_air(h, 0, j, k-2))){
                    for(int n=0; n<last_array_index_1; n++)
                        dxdphi_vals[n] = dxdphi_c(arrays_1[n]);
                    phi_flag = true;
                }
            }
            for(int n=0; n<last_array_index_1; n++){
                if(!the_flag) dxdthe_vals[n] = dxdthe_a(arrays_1[n]);
                if(!phi_flag) dxdphi_vals[n] = dxdphi_a(arrays_1[n]);
            }
            double dvdthe = dxdthe_vals[i_v_1];
            double dwdphi = dxdphi_vals[i_w_1];
            EkmanPumping.y[j][k] = - coeff_pumping 
                /r_salt_water.x[i_max][j][k] 
                * (dvdthe/rm + dwdphi/rmsinthe);
/*
            EkmanPumping.y[j][k] = - coeff_pumping // no refinement around walls
                /r_salt_water.x[i_max][j][k] 
                * ((aux_v.x[i_max][j+1][k] - aux_v.x[i_max][j-1][k])
                /(2. * rm * dthe) 
                + (aux_w.x[i_max][j][k+1] - aux_w.x[i_max][j][k-1])
                /(2. * rmsinthe * dphi));
*/
/*
            if((j>=88)&&(j<=90)){
                EkmanPumping.y[j][k] = EkmanPumping.y[87][k];
            }
            if((j>90)&&(j<=92)){
                EkmanPumping.y[j][k] = EkmanPumping.y[87][k];
            }
*/
/*
            if(EkmanPumping.y[j][k] >= 50.0)  
                EkmanPumping.y[j][k] = 50.0; // in cm/d
            if(EkmanPumping.y[j][k] < 50.0)  
                EkmanPumping.y[j][k] = -50.0; // in cm/d
*/
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
                    Salt_total.y[j][k] += c.x[i][j][k];
                }
            }
        }
    }
    double coeff_buoy = r_0_water * (u_0 * u_0)/L_hyd; // coefficient for bouancy term = 0.2871
    double coeff_Coriolis = r_0_water * u_0; // coefficient for Coriolis term = 239.28
    double coriolis = 1.0;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Salt_Finger.x[0][j][k] = c43 * Salt_Finger.x[1][j][k] -
                c13 * Salt_Finger.x[2][j][k];
//            Salt_Finger.x[i_max][j][k] = c43 * Salt_Finger.x[im-2][j][k] -
//                c13 * Salt_Finger.x[im-3][j][k];
            Salt_Diffusion.x[0][j][k] = c43 * Salt_Diffusion.x[1][j][k] -
                c13 * Salt_Diffusion.x[2][j][k];
//            Salt_Diffusion.x[i_max][j][k] = c43 * Salt_Diffusion.x[im-2][j][k] -
//                c13 * Salt_Diffusion.x[im-3][j][k];
            Salt_Balance.x[0][j][k] = c43 * Salt_Balance.x[1][j][k] -
                c13 * Salt_Balance.x[2][j][k];
//            Salt_Balance.x[i_max][j][k] = c43 * Salt_Balance.x[im-2][j][k] -
//                c13 * Salt_Balance.x[im-3][j][k];
            Salt_Finger.x[im-1][j][k] = 
                Salt_Finger.x[im-4][j][k] 
                - 3.0 * Salt_Finger.x[im-3][j][k] 
                + 3.0 * Salt_Finger.x[im-2][j][k];  // extrapolation
            Salt_Diffusion.x[im-1][j][k] = 
                Salt_Diffusion.x[im-4][j][k] 
                - 3.0 * Salt_Diffusion.x[im-3][j][k] 
                + 3.0 * Salt_Diffusion.x[im-2][j][k];  // extrapolation
            Salt_Balance.x[im-1][j][k] = 
                Salt_Balance.x[im-4][j][k] 
                - 3.0 * Salt_Balance.x[im-3][j][k] 
                + 3.0 * Salt_Balance.x[im-2][j][k];  // extrapolation
            BuoyancyForce.x[0][j][k] = 
                BuoyancyForce.x[3][j][k] 
                - 3.0 * BuoyancyForce.x[2][j][k] 
                + 3.0 * BuoyancyForce.x[1][j][k];  // extrapolation
            CoriolisForce.x[0][j][k] = 
                CoriolisForce.x[3][j][k] 
                - 3.0 * CoriolisForce.x[2][j][k] 
                + 3.0 * CoriolisForce.x[1][j][k];  // extrapolation
            PressureGradientForce.x[0][j][k] = 
                PressureGradientForce.x[3][j][k] 
                - 3.0 * PressureGradientForce.x[2][j][k] 
                + 3.0 * PressureGradientForce.x[1][j][k];  // extrapolation

            double dpdr = p_dyn.x[i_max][j][k] = p_dyn.x[im-4][j][k] 
                - 3. * p_dyn.x[im-3][j][k] + 3. * p_dyn.x[im-2][j][k];  // extrapolation
            double sinthe_coriolis = sin(the.z[j]);
            double costhe = cos(the.z[j]);
            double coriolis_rad = coriolis * 2.0 * omega
                * costhe * w.x[i_max][j][k];
            double coriolis_the = - coriolis * 2.0 * omega
                * sinthe_coriolis * w.x[i_max][j][k];
            double coriolis_phi = coriolis * 2.0 * omega
                * (sinthe_coriolis * v.x[i_max][j][k] 
                - costhe * u.x[i_max][j][k]);
            CoriolisForce.x[i_max][j][k] = coeff_Coriolis 
                * sqrt((pow (coriolis_rad,2) 
                + pow (coriolis_the,2) 
                + pow (coriolis_phi,2))/3.0);
            BuoyancyForce.x[i_max][j][k] = buoyancy 
                * r_0_water * (t.x[i_max][j][k] - 1.0) * g;
            PressureGradientForce.x[i_max][j][k] = - coeff_buoy * dpdr;
            if(is_land(h, i_max, j, k)){
                BuoyancyForce.x[i_max][j][k] = 0.0;
                PressureGradientForce.x[i_max][j][k] = 0.0;
                CoriolisForce.x[i_max][j][k] = 0.0;
            }
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
void cHydrosphereModel::print_min_max_hyd(){
    cout << endl << " flow properties: " << endl << endl;
    searchMinMax_3D(" max temperature ", " min temperature ", 
        " deg", t, 273.15, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max u-component ", " min u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max v-component ", " min v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max w-component ", " min w-component ", "m/s", w, u_0);
    searchMinMax_3D(" max pressure dynamic ", " min pressure dynamic ", "hPa", p_dyn, 1.0);
    searchMinMax_3D(" max pressure static ", " min pressure static ", "bar", p_stat, 1.0);
    searchMinMax_3D(" max water density ", " min water density ", "kg/m3", r_water, 1.0);
    searchMinMax_3D(" max salt water density ", " min salt water density ", "kg/m3", r_salt_water, 1.0);
    cout << endl << " salinity based results in the three dimensional space: " << endl << endl;
    searchMinMax_3D(" max salt concentration ", " min salt concentration ", "psu", c, c_0);
    searchMinMax_3D(" max salt balance ", " min salt balance ", "kg/m3", Salt_Balance, 1.0);
    searchMinMax_3D(" max salt finger ", " min salt finger ", "kg/m3", Salt_Finger, 1.0);
    searchMinMax_3D(" max salt diffusion ", " min salt diffusion ", "kg/m3", Salt_Diffusion, 1.0);
    cout << endl << " forces per unit volume: " << endl << endl;
    searchMinMax_3D(" max pressure force ", " min pressure force ", "N", PressureGradientForce, 1.0);
    searchMinMax_3D(" max buoyancy force ", " min buoyancy force ", "N", BuoyancyForce, 1.0);
    searchMinMax_3D(" max Coriolis force ", " min Coriolis force ", "N", CoriolisForce, 1.0);
    cout << endl << " salt concentration averaged for the two dimensional surface plane: " << endl << endl;
    searchMinMax_2D(" max salt total ", " min salt total ", "psu", Salt_total, c_0);
    searchMinMax_2D(" max Salt_Finger ", " min Salt_Finger ", "kg/m3", SaltFinger, 1.0);
    searchMinMax_2D(" max Salt_Diffusion ", " min Salt_Diffusion ", "kg/m3", SaltDiffusion, 1.0);
    cout << endl << " deep currents averaged for a two dimensional plane: " << endl << endl;
    searchMinMax_2D(" max EkmanPumping ", " min EkmanPumping ", "cm/d", EkmanPumping, 1.0);
    searchMinMax_2D(" max upwelling ", " min upwelling ", "cm/d", Upwelling, 1.0);
    searchMinMax_2D(" max downwelling ", " min downwelling ", "cm/d", Downwelling, 1.0);
    searchMinMax_2D(" max bathymetry ", " min bathymetry ", "m", Bathymetry, 1.0);
}
/*
 * 
*/

