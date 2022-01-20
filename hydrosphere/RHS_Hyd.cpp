/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include <iostream>
#include <cmath>
#include "cHydrosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

#define dxdr_a(X) \
    (h_d_i * (X->x[i+1][j][k] - X->x[i-1][j][k])/(2.0 * dr))
#define dxdr_b(X) \
    (h_d_i * (- 3.0 * X->x[i][j][k] + 4.0 * X->x[i+1][j][k] - X->x[i+2][j][k])/(2.0 * dr))
#define dxdr_d(X) \
    (h_d_i * (X->x[i+1][j][k] - X->x[i][j][k])/dthe)

#define d2xdr2_a(X) \
    (h_d_i * (X->x[i+1][j][k] - 2.0 * X->x[i][j][k] + X->x[i-1][j][k])/dr2)
#define d2xdr2_b(X) \
    (h_d_i * (X->x[i][j][k] - 2.0 * X->x[i+1][j][k] + X->x[i+2][j][k])/dr2)

#define dxdthe_a(X) \
    ((X->x[i][j+1][k] - X->x[i][j-1][k])/(2.0 * dthe))
#define dxdthe_b(X) \
    ((- 3.0 * X->x[i][j][k] + 4.0 * X->x[i][j+1][k] - X->x[i][j+2][k])/(2.0 * dthe))
#define dxdthe_c(X) \
    (- (- 3.0 * X->x[i][j][k] + 4.0 * X->x[i][j-1][k] - X->x[i][j-2][k])/(2.0 * dthe))
#define dxdthe_d(X) \
    ((X->x[i][j+1][k] - X->x[i][j][k])/dthe)
#define dxdthe_e(X) \
    ((X->x[i][j-1][k] - X->x[i][j][k])/dthe)

#define d2xdthe2_a(X) \
    ((X->x[i][j+1][k] - 2.0 * X->x[i][j][k] + X->x[i][j-1][k])/dthe2)
#define d2xdthe2_b(X) \
    ((X->x[i][j][k] - 2.0 * X->x[i][j+1][k] + X->x[i][j+2][k])/dthe2)
#define d2xdthe2_c(X) \
    (- (X->x[i][j][k] - 2.0 * X->x[i][j-1][k] + X->x[i][j-2][k])/dthe2)

#define dxdphi_a(X) \
    ((X->x[i][j][k+1] - X->x[i][j][k-1])/(2.0 * dphi))
#define dxdphi_b(X) \
    ((- 3.0 * X->x[i][j][k] + 4.0 * X->x[i][j][k+1] - X->x[i][j][k+2])/(2.0 * dphi))
#define dxdphi_c(X) \
    (- (- 3.0 * X->x[i][j][k] + 4.0 * X->x[i][j][k-1] - X->x[i][j][k-2])/(2.0 * dphi))
#define dxdphi_d(X) \
    ((X->x[i][j][k+1] - X->x[i][j][k])/dphi)
#define dxdphi_e(X) \
    ((X->x[i][j][k-1] - X->x[i][j][k])/dphi)

#define d2xdphi2_a(X) \
    ((X->x[i][j][k+1] - 2.0 * X->x[i][j][k] + X->x[i][j][k-1])/dphi2)
#define d2xdphi2_b(X) \
    ((X->x[i][j][k] - 2.0 * X->x[i][j][k+1] + X->x[i][j][k+2])/dphi2)
#define d2xdphi2_c(X) \
    (- (X->x[i][j][k] - 2.0 * X->x[i][j][k-1] + X->x[i][j][k-2])/dphi2)
/*
*
*/
void cHydrosphereModel::RK_RHS_3D_Hydrosphere(int i, int j, int k){
//  (Reinout vander Meulen, The immersed Boundary Method for the Incompressible Navier-Stokes Equations)
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double rm = rad.z[i];
    double rm2 = rm * rm;
    double sinthe = sin(the.z[j]);
    if(sinthe == 0.0) sinthe = 1.0e-5;
    double sinthe2 = sinthe * sinthe;
    double costhe = cos(the.z[j]);
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;
    double h_0_i = 0.0, h_d_i = 0.0;
    double topo_step = L_hyd/(double)(im-1);
    double height = L_hyd - (double)i * topo_step;
    double topo_diff = fabs(height - Bathymetry.y[j][k]);
    if((topo_diff <= topo_step)
        &&((is_water(h, i, j, k)) 
        &&(is_land(h, i-1, j, k)))){
//        h_0_i = 0.5 *(acos(topo_diff * 3.14/L_hyd) + 1.);   // cosine distribution function
        h_d_i = 1.0 - topo_diff/topo_step;   // hat distribution function
        h_0_i = 0.0;
    }
        if(is_water(h, i, j, k)){
        h_0_i = 0.0;
        h_d_i = 1.0;
    }
    if(is_land(h, i, j, k)){
        h_0_i = 1.0;
        h_d_i = 0.0;
    }
    std::vector<Array*> arrays_1{&u, &v, &w, &t, &p_dyn, &c};
    std::vector<Array*> arrays_2{&u, &v, &w, &t, &c};
    enum array_index_1{i_u_1, i_v_1, i_w_1, i_t_1, i_p_1, i_c_1, 
        last_array_index_1};
    enum array_index_2{i_u_2, i_v_2, i_w_2, i_t_2, i_c_2,
        last_array_index_2};
    std::vector<double> dxdr_vals(last_array_index_1), 
                        dxdthe_vals(last_array_index_1), 
                        dxdphi_vals(last_array_index_1),
                        d2xdr2_vals(last_array_index_2),
                        d2xdthe2_vals(last_array_index_2),
                        d2xdphi2_vals(last_array_index_2);
    bool r_flag = false, the_flag = false, phi_flag = false;
    if(i < im-2){
        if((is_land(h, i, j, k))
            &&(is_water(h, i+1, j, k))){
            for(int n=0; n<last_array_index_1; n++)
                dxdr_vals[n] = dxdr_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdr2_vals[n] = d2xdr2_b(arrays_2[n]);
            r_flag = true;
        }
    }
    if(i == im-2){
        if((is_land(h, i, j, k)) 
            &&(is_water(h, i+1, j, k))){
            for(int n=0; n<last_array_index_1; n++)
                dxdr_vals[n] = dxdr_d(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdr2_vals[n] = 0.0;
            r_flag = true;
        }
    }
    if((j >= 1)&&(j < jm-2)){
        if((is_land(h, i, j, k))
            &&((is_water(h, i, j+1, k)) 
            &&(is_water(h, i, j+2, k)))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_b(arrays_2[n]);
            the_flag = true;
        }
        if((is_land(h, i, j, k))
            &&(is_water(h, i, j-1, k)) 
            &&(is_water(h, i, j-2, k))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_c(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_c(arrays_2[n]);
            the_flag = true;
        }
        if(((is_land(h, i, j, k)) 
            &&((is_water(h, i, j+1, k))
            &&(is_land(h, i, j+2, k)))) 
            ){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_d(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = 0.0;
            the_flag = true;
        }
        if(((is_land(h, i, j, k)) 
            &&((is_water(h, i, j-1, k))
            &&(is_land(h, i, j-2, k)))) 
            ){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_e(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = 0.0;
            the_flag = true;
        }
    }
    if((k >= 1)&&(k < km-2)){
        if((is_land(h, i, j, k))
            &&(is_water(h, i, j, k+1)) 
            &&(is_water(h, i, j, k+2))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_b(arrays_2[n]);
            phi_flag = true;
        }
        if((is_land(h, i, j, k))
            &&(is_water(h, i, j, k-1)) 
            &&(is_water(h, i, j, k-2))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_c(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_c(arrays_2[n]);
            phi_flag = true;
        }

        if(((is_land(h, i, j, k))
            &&((is_water(h, i, j, k+1)) 
            &&(is_land(h, i, j, k+2))))
            ){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_d(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = 0.0;
            phi_flag = true;
        }
        if(((is_land(h, i, j, k)) 
            &&((is_water(h, i, j, k-1))
            &&(is_land(h, i, j, k-2)))) 
            ){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_e(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = 0.0;
            phi_flag = true;
        }
    }

    for(int n=0; n<last_array_index_1; n++){
        if(!r_flag) dxdr_vals[n] = dxdr_a(arrays_1[n]);
        if(!the_flag) dxdthe_vals[n] = dxdthe_a(arrays_1[n]);
        if(!phi_flag) dxdphi_vals[n] = dxdphi_a(arrays_1[n]);
    }
    for(int n=0; n<last_array_index_2; n++){
        if(!r_flag) d2xdr2_vals[n] = d2xdr2_a(arrays_2[n]);
        if(!the_flag) d2xdthe2_vals[n] = d2xdthe2_a(arrays_2[n]);
        if(!phi_flag) d2xdphi2_vals[n] = d2xdphi2_a(arrays_2[n]);
    }
    double dudr = dxdr_vals[i_u_1], dvdr = dxdr_vals[i_v_1], 
           dwdr = dxdr_vals[i_w_1], dtdr = dxdr_vals[i_t_1],
           dpdr = dxdr_vals[i_p_1], dcdr = dxdr_vals[i_c_1];
    double dudthe = dxdthe_vals[i_u_1], dvdthe = dxdthe_vals[i_v_1], 
           dwdthe = dxdthe_vals[i_w_1], dtdthe = dxdthe_vals[i_t_1],
           dpdthe = dxdthe_vals[i_p_1], dcdthe = dxdthe_vals[i_c_1];
    double dudphi = dxdphi_vals[i_u_1], dvdphi = dxdphi_vals[i_v_1], 
           dwdphi = dxdphi_vals[i_w_1], dtdphi = dxdphi_vals[i_t_1],
           dpdphi = dxdphi_vals[i_p_1], dcdphi = dxdphi_vals[i_c_1];

    double d2udr2 = d2xdr2_vals[i_u_2], d2vdr2 = d2xdr2_vals[i_v_2], 
           d2wdr2 = d2xdr2_vals[i_w_2], d2tdr2 = d2xdr2_vals[i_t_2],
           d2cdr2 = d2xdr2_vals[i_c_2];
    double d2udthe2 = d2xdthe2_vals[i_u_2], d2vdthe2 = d2xdthe2_vals[i_v_2], 
           d2wdthe2 = d2xdthe2_vals[i_w_2], d2tdthe2 = d2xdthe2_vals[i_t_2],
           d2cdthe2 = d2xdthe2_vals[i_c_2];
    double d2udphi2 = d2xdphi2_vals[i_u_2], d2vdphi2 = d2xdphi2_vals[i_v_2], 
           d2wdphi2 = d2xdphi2_vals[i_w_2], d2tdphi2 = d2xdphi2_vals[i_t_2],
           d2cdphi2 = d2xdphi2_vals[i_c_2];

    // buoyancy effects by salt water density changes
    double drodc = 0.7;    // gradient given in kg/m³/m
    double salt_water_ref = r_water.x[i][j][k] + drodc * c.x[i][j][k] * c_0;
                          // common linear approach for salt water based on fresh water
    double coeff_g_p = 1.0; // coefficient allows the buoyancy term
    double coeff_energy_p = u_0 * u_0/(cp_w * t_0); // coefficient for the source term = 2.33e-4
    double coeff_buoy = r_0_water * (u_0 * u_0)/L_hyd; // coefficient for bouancy term = 0.2871
    double coeff_Coriolis = r_0_water * u_0; // coefficient for Coriolis term = 239.28
    // Boussineq-approximation for the buoyancy force caused by salinity, higher salinity causes negative buoyancy
    Salt_Balance.x[i][j][k] = salt_water_ref - r_salt_water.x[i][j][k]; // difference of salinity compared to average
    if(Salt_Balance.x[i][j][k] < 0.0){
        Salt_Diffusion.x[i][j][k] = Salt_Balance.x[i][j][k]; // for negativ salinity balance, higher than salt_water_ref
        Salt_Finger.x[i][j][k] = 0.0;
    }else{
        Salt_Finger.x[i][j][k] = Salt_Balance.x[i][j][k]; // for positiv salinity balance, lower than salt_water_ref
        Salt_Diffusion.x[i][j][k] = 0.0;
    }
    double coriolis = 1.0;
    double coriolis_rad = h_d_i * coriolis * 2.0 * omega
        * costhe * w.x[i][j][k];
    double coriolis_the = - coriolis * 2.0 * omega
        * sinthe * w.x[i][j][k];
    double coriolis_phi = coriolis * 2.0 * omega
        * (sinthe * v.x[i][j][k] - costhe * u.x[i][j][k]);
    CoriolisForce.x[i][j][k] = coeff_Coriolis * sqrt((pow (coriolis_rad,2) 
        + pow (coriolis_the, 2) + pow (coriolis_phi,2))/3.0);
    BuoyancyForce.x[i][j][k] = buoyancy * r_0_water 
        * (t.x[i][j][k] - 1.0) * g;
    PressureGradientForce.x[i][j][k] = - coeff_buoy * dpdr;
    if(is_land(h, i, j, k)){
        BuoyancyForce.x[i][j][k] = 0.0;
        PressureGradientForce.x[i][j][k] = 0.0;
        CoriolisForce.x[i][j][k] = 0.0;
    }

    if(((j <= 10)||(j <= jm-11))||((k <= 10)||(k <= km-11))){ // prevents oscillations along pole regions (Euler-equations)
        rhs_t.x[i][j][k] = - (u.x[i][j][k] * dtdr + v.x[i][j][k] * dtdthe/rm 
            + w.x[i][j][k] * dtdphi/rmsinthe) 
            + (dtdr * 2.0/rm 
            + dtdthe * costhe/rm2sinthe)/(re * pr)
            + coeff_energy_p * (u.x[i][j][k] * dpdr 
            + v.x[i][j][k]/rm * dpdthe 
            + w.x[i][j][k]/rmsinthe * dpdphi)
            - h_0_i * t.x[i][j][k]/dr2;
        rhs_u.x[i][j][k] = - (u.x[i][j][k] * dudr + v.x[i][j][k] * dudthe/rm 
            + w.x[i][j][k] * dudphi/rmsinthe) 
            + coeff_g_p * (t.x[i][j][k] - 1.0) * g
            - dpdr
            + ( h_d_i * 2.0 * u.x[i][j][k]/rm2 
            + 4.0 * dudr/rm + dudthe * costhe/rm2sinthe 
            )/re 
            + coriolis_rad
            - h_0_i * u.x[i][j][k]/dr2;
        rhs_v.x[i][j][k] = - (u.x[i][j][k] * dvdr + v.x[i][j][k] * dvdthe/rm 
            + w.x[i][j][k] * dvdphi/rmsinthe) 
            - dpdthe/rm
            + (dvdr * 2.0/rm 
            + dvdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * v.x[i][j][k] 
            + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re 
            + coriolis_the
            - h_0_i * v.x[i][j][k]/dr2;
        rhs_w.x[i][j][k] = - (u.x[i][j][k] * dwdr + v.x[i][j][k] * dwdthe/rm 
            + w.x[i][j][k] * dwdphi/rmsinthe) 
            - dpdphi/rmsinthe
            + (dwdr * 2.0/rm 
            + dwdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * w.x[i][j][k] 
            + 2.0 * dudphi/rm2sinthe 
//            + dvdphi * 2.0 * costhe/rm2sinthe2
            )/re 
            + coriolis_phi
            - h_0_i * w.x[i][j][k]/dr2;
        rhs_c.x[i][j][k] = - (u.x[i][j][k] * dcdr + v.x[i][j][k] * dcdthe/rm 
            + w.x[i][j][k] * dcdphi/rmsinthe)
            - h_0_i * c.x[i][j][k]/dr2;
    }

    else{
        rhs_t.x[i][j][k] = - (u.x[i][j][k] * dtdr + v.x[i][j][k] * dtdthe/rm  // Navier-Stokes equations
            + w.x[i][j][k] * dtdphi/rmsinthe) 
            + (d2tdr2 + dtdr * 2.0/rm + d2tdthe2/rm2 
            + dtdthe * costhe/rm2sinthe + d2tdphi2/rm2sinthe2)/(re * pr)
            + coeff_energy_p * (u.x[i][j][k] * dpdr 
            + v.x[i][j][k]/rm * dpdthe 
            + w.x[i][j][k]/rmsinthe * dpdphi)
            - h_0_i * t.x[i][j][k]/dr2;
        rhs_u.x[i][j][k] = - (u.x[i][j][k] * dudr + v.x[i][j][k] * dudthe/rm 
            + w.x[i][j][k] * dudphi/rmsinthe) 
            + coeff_g_p * (t.x[i][j][k] - 1.0) * g
            - dpdr
            + (d2udr2 + h_d_i * 2.0 * u.x[i][j][k]/rm2 + d2udthe2/rm2 
            + 4.0 * dudr/rm + dudthe * costhe/rm2sinthe 
            + d2udphi2/rm2sinthe2)/re 
            + coriolis_rad
            - h_0_i * u.x[i][j][k]/dr2;
        rhs_v.x[i][j][k] = - (u.x[i][j][k] * dvdr + v.x[i][j][k] * dvdthe/rm 
            + w.x[i][j][k] * dvdphi/rmsinthe) 
            - dpdthe/rm
            + (d2vdr2 + dvdr * 2.0/rm + d2vdthe2/rm2 
            + dvdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * v.x[i][j][k] 
            + d2vdphi2/rm2sinthe2 + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re 
            + coriolis_the
            - h_0_i * v.x[i][j][k]/dr2;
        rhs_w.x[i][j][k] = - (u.x[i][j][k] * dwdr + v.x[i][j][k] * dwdthe/rm 
            + w.x[i][j][k] * dwdphi/rmsinthe) 
            - dpdphi/rmsinthe
            + (d2wdr2 + dwdr * 2.0/rm + d2wdthe2/rm2 
            + dwdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * w.x[i][j][k] 
            + d2wdphi2/rm2sinthe2 + 2.0 * dudphi/rm2sinthe 
            + dvdphi * 2.0 * costhe/rm2sinthe2)/re 
            + coriolis_phi
            - h_0_i * w.x[i][j][k]/dr2;
        rhs_c.x[i][j][k] = - (u.x[i][j][k] * dcdr + v.x[i][j][k] * dcdthe/rm 
            + w.x[i][j][k] * dcdphi/rmsinthe) 
            + (d2cdr2 + dcdr * 2.0/rm 
            + d2cdthe2/rm2 
            + dcdthe * costhe/rm2sinthe 
            + d2cdphi2/rm2sinthe2)/(sc * re)
            - h_0_i * c.x[i][j][k]/dr2;
    }

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe/rm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi/rmsinthe;

    if(is_land(h, i, j, k)){
        u.x[i][j][k] = 0.0;
        v.x[i][j][k] = 0.0;
        w.x[i][j][k] = 0.0;
        p_dyn.x[i][j][k] = 0.0;
        rhs_u.x[i][j][k] = 0.0;
        rhs_v.x[i][j][k] = 0.0;
        rhs_w.x[i][j][k] = 0.0;
        rhs_t.x[i][j][k] = 0.0;
        rhs_c.x[i][j][k] = 0.0;
        aux_u.x[i][j][k] = 0.0;
        aux_v.x[i][j][k] = 0.0;
        aux_w.x[i][j][k] = 0.0;
    }
/*
    double residuum = (u.x[i][j][k] - u.x[i-1][j][k])/dr 
        + ((v.x[i-1][j+1][k] - v.x[i-1][j-1][k])
        /(2. * rm * dthe) 
        + (w.x[i-1][j][k+1] - w.x[i-1][j][k-1])
        /(2. * rmsinthe * dphi));

    cout.precision(8);
    cout << "Ocean code" << endl
        << "northern hemisphere" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   dt = " << dt 
        << "   dr = " << dr 
        << "   dthe = " << dthe 
        << "   dphi = " << dphi
        << "   dr2 = " << dr2 
        << "   dthe2 = " << dthe2 
        << "   dphi2 = " << dphi2 << endl
        << "   rm > rad = " << rm 
        << "   sinthe = " << sinthe
        << "   sinthe2 = " << sinthe2
        << "   costhe = " << costhe
        << "   rmsinthe = " << rmsinthe
        << "   rm2sinthe = " << rm2sinthe
        << "   rm2sinthe2 = " << rm2sinthe2 << endl
        << "   topo_step = " << topo_step
        << "   height = " << height
        << "   Bathymetry = " << Bathymetry.y[j][k]
        << "   topo_diff = " << topo_diff << endl
        << "   h_0_i = " << h_0_i
        << "   h_d_i = " << h_d_i << endl
        << "   residuum = " << residuum << endl
        << "   dpdr = " << dpdr
        << "   dtdr = " << dtdr
        << "   dudr = " << dudr
        << "   dvdr = " << dvdr
        << "   dwdr = " << dwdr
        << "   dcdr = " << dcdr << endl
        << "   dpdthe = " << dpdthe
        << "   dtdthe = " << dtdthe
        << "   dudthe = " << dudthe
        << "   dvdthe = " << dvdthe
        << "   dwdthe = " << dwdthe
        << "   dcdthe = " << dcdthe << endl
        << "   dpdphi = " << dpdphi
        << "   dtdphi = " << dtdphi
        << "   dudphi = " << dudphi
        << "   dvdphi = " << dvdphi
        << "   dwdphi = " << dwdphi
        << "   dcdphi = " << dcdphi << endl

        << "   ti-1 = " << t.x[i-1][j][k]
        << "   ui-1 = " << u.x[i-1][j][k]
        << "   vi-1 = " << v.x[i-1][j][k]
        << "   wi-1 = " << w.x[i-1][j][k]
        << "   ci-1 = " << t.x[i-1][j][k] << endl

        << "   ti   = " << t.x[i][j][k]
        << "   ui   = " << u.x[i][j][k]
        << "   vi   = " << v.x[i][j][k]
        << "   wi  = " << w.x[i][j][k]
        << "   ci   = " << t.x[i][j][k] << endl

        << "   ti+1 = " << t.x[i+1][j][k]
        << "   ui+1 = " << u.x[i+1][j][k]
        << "   vi+1 = " << v.x[i+1][j][k]
        << "   wi+1 = " << w.x[i+1][j][k]
        << "   ci+1 = " << t.x[i+1][j][k] << endl

        << "   tj-1 = " << t.x[i][j-1][k]
        << "   uj-1 = " << u.x[i][j-1][k]
        << "   vj-1 = " << v.x[i][j-1][k]
        << "   wj-1 = " << w.x[i][j-1][k]
        << "   cj-1 = " << t.x[i][j-1][k] << endl

        << "   tj   = " << t.x[i][j][k]
        << "   uj   = " << u.x[i][j][k]
        << "   vj   = " << v.x[i][j][k]
        << "   wj   = " << w.x[i][j][k]
        << "   cj   = " << t.x[i][j][k] << endl

        << "   tj+1 = " << t.x[i][j+1][k]
        << "   uj+1 = " << u.x[i][j+1][k]
        << "   vj+1 = " << v.x[i][j+1][k]
        << "   wj+1 = " << w.x[i][j+1][k]
        << "   cj+1 = " << t.x[i][j+1][k] << endl

        << "   tk-1 = " << t.x[i][j][k-1]
        << "   uk-1 = " << u.x[i][j][k-1]
        << "   vk-1 = " << v.x[i][j][k-1]
        << "   wk-1 = " << w.x[i][j][k-1]
        << "   ck-1 = " << t.x[i][j][k-1] << endl

        << "   tk   = " << t.x[i][j][k]
        << "   uk   = " << u.x[i][j][k]
        << "   vk   = " << v.x[i][j][k]
        << "   wk   = " << w.x[i][j][k]
        << "   ck   = " << t.x[i][j][k] << endl

        << "   tk+1 = " << t.x[i][j][k+1]
        << "   uk+1 = " << u.x[i][j][k+1]
        << "   vk+1 = " << v.x[i][j][k+1]
        << "   wk+1 = " << w.x[i][j][k+1]
        << "   ck+1 = " << t.x[i][j][k+1] << endl

        << "   tn = " << tn.x[i][j][k]
        << "   un = " << un.x[i][j][k]
        << "   vn = " << vn.x[i][j][k]
        << "   wn = " << wn.x[i][j][k]
        << "   cn = " << cn.x[i][j][k] << endl

        << "   rhs_t = " << rhs_t.x[i][j][k]
        << "   rhs_u = " << rhs_u.x[i][j][k]
        << "   rhs_v = " << rhs_v.x[i][j][k]
        << "   rhs_w = " << rhs_w.x[i][j][k]
        << "   rhs_c = " << rhs_t.x[i][j][k] << endl;
*/
}
/*
*
*/
void cHydrosphereModel::RK_RHS_2D_Hydrosphere(int i, int j, int k){
//  (Reinout vander Meulen, The immersed Boundary Method for the Incompressible Navier-Stokes Equations)
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double rm = rad.z[i];
    double rm2 = rm * rm;
    double sinthe = sin(the.z[j]);
    if(sinthe == 0.0) sinthe = 1.0e-5;
    double sinthe2 = sinthe * sinthe;
    double costhe = cos(the.z[j]);
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;
    double h_0_j = 0.0, h_d_j = 0.0;
    double h_0_k = 0.0, h_d_k = 0.0;
    if(is_water(h, im-1, j, k)){
        h_0_j = 0.0;
        h_0_k = 0.0;
        h_d_j = 1.0;
        h_d_k = 1.0;
    }
    if(is_land(h, im-1, j, k)){
        h_0_j = 1.0;
        h_0_k = 1.0;
        h_d_j = 0.0;
        h_d_k = 0.0;
    }
    std::vector<Array*> arrays_1{&u, &v, &w, &p_dyn};
    std::vector<Array*> arrays_2{&v, &w};
    enum array_index_1{i_u_1, i_v_1, i_w_1, i_p_1, last_array_index_1};
    enum array_index_2{i_v_2, i_w_2, last_array_index_2};
    std::vector<double> dxdthe_vals(last_array_index_1), 
                        dxdphi_vals(last_array_index_1),
                        d2xdthe2_vals(last_array_index_2),
                        d2xdphi2_vals(last_array_index_2);
    bool the_flag = false, phi_flag = false;

    if((j >= 1)&&(j < jm-2)){
        if((is_land(h, im-1, j, k))
            &&((is_water(h, im-1, j+1, k)) 
            &&(is_water(h, im-1, j+2, k)))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_b(arrays_2[n]);
            the_flag = true;
        }
        if((is_land(h, im-1, j, k))
            &&(is_water(h, im-1, j-1, k)) 
            &&(is_water(h, im-1, j-2, k))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_c(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_c(arrays_2[n]);
            the_flag = true;
        }
    }
    if((k >= 1)&&(k < km-2)){
        if((is_land(h, im-1, j, k))
            &&(is_water(h, im-1, j, k+1)) 
            &&(is_water(h, im-1, j, k+2))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_b(arrays_2[n]);
            phi_flag = true;
        }
        if((is_land(h, im-1, j, k))
            &&(is_water(h, im-1, j, k-1)) 
            &&(is_water(h, im-1, j, k-2))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_c(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_c(arrays_2[n]);
            phi_flag = true;
        }
    }

    for(int n=0; n<last_array_index_1; n++){
        if(!the_flag) dxdthe_vals[n] = dxdthe_a(arrays_1[n]);
        if(!phi_flag) dxdphi_vals[n] = dxdphi_a(arrays_1[n]);
    }
    for(int n=0; n<last_array_index_2; n++){
        if(!the_flag) d2xdthe2_vals[n] = d2xdthe2_a(arrays_2[n]);
        if(!phi_flag) d2xdphi2_vals[n] = d2xdphi2_a(arrays_2[n]);
    }
    double dudthe = dxdthe_vals[i_u_1], dvdthe = dxdthe_vals[i_v_1], 
           dwdthe = dxdthe_vals[i_w_1], dpdthe = dxdthe_vals[i_p_1];
    double dudphi = dxdphi_vals[i_u_1], dvdphi = dxdphi_vals[i_v_1], 
           dwdphi = dxdphi_vals[i_w_1], dpdphi = dxdphi_vals[i_p_1];
    double d2vdthe2 = d2xdthe2_vals[i_v_2], d2wdthe2 = d2xdthe2_vals[i_w_2];
    double d2vdphi2 = d2xdphi2_vals[i_v_2], d2wdphi2 = d2xdphi2_vals[i_w_2];

    if(((j <= 10)||(j <= jm-11))||((k <= 10)||(k <= km-11))){ // prevents oscillations along pole regions
        rhs_v.x[im-1][j][k] = - (v.x[im-1][j][k] * dvdthe/rm 
            + w.x[im-1][j][k] * dvdphi/rmsinthe) 
            - dpdthe/rm 
            + (dvdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_j * v.x[im-1][j][k] 
            + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re
            - h_0_j * v.x[im-1][j][k]/dthe2;                  // RHS v-velocity Euler equations without friction terms
        rhs_w.x[im-1][j][k] = - (v.x[im-1][j][k] * dwdthe/rm 
            + w.x[im-1][j][k] * dwdphi/rmsinthe) 
            - dpdphi/rmsinthe 
            + (dwdthe/rm2sinthe  * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_k * w.x[im-1][j][k] 
            + 2.0 * dudphi/rm2sinthe 
//            + dvdphi * 2.0 * costhe/rm2sinthe2
            )/re
            - h_0_k * w.x[im-1][j][k]/dphi2;                  // RHS w-velocity Euler equations without friction terms
    }else{
        rhs_v.x[im-1][j][k] = - (v.x[im-1][j][k] * dvdthe/rm 
            + w.x[im-1][j][k] * dvdphi/rmsinthe) 
            - dpdthe/rm 
            + (d2vdthe2/rm2 + dvdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_j * v.x[im-1][j][k] 
            + d2vdphi2/rm2sinthe2 + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re
            - h_0_j * v.x[im-1][j][k]/dthe2;                  // RHS v-velocity Navier-Stokes equations with friction terms

        rhs_w.x[im-1][j][k] = - (v.x[im-1][j][k] * dwdthe/rm 
            + w.x[im-1][j][k] * dwdphi/rmsinthe) 
            - dpdphi/rmsinthe 
            + (d2wdthe2/rm2 + dwdthe/rm2sinthe  * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_k * w.x[im-1][j][k] 
            + d2wdphi2/rm2sinthe2 + 2.0 * dudphi/rm2sinthe 
            + dvdphi * 2.0 * costhe/rm2sinthe2
            )/re
            - h_0_k * w.x[im-1][j][k]/dphi2;                  // RHS w-velocity Navier-Stokes equations with friction terms
    }

    aux_v.x[im-1][j][k] = rhs_v.x[im-1][j][k] + dpdthe/rm;
    aux_w.x[im-1][j][k] = rhs_w.x[im-1][j][k] + dpdphi/rmsinthe;

    if(is_land(h, im-1, j, k)){
        rhs_v.x[im-1][j][k] = 0.0; 
        rhs_w.x[im-1][j][k] = 0.0;
        aux_v.x[im-1][j][k] = 0.0; 
        aux_w.x[im-1][j][k] = 0.0;
    }
/*
    cout.precision(8);
//    if((j == 90) &&(k == 180)) cout << "Ocean code" << endl
        << "northern hemisphere" << endl
        << "   im-1 = " << im-1 << "   j = " << j << "   k = " << k  << endl
        << "   dt = " << dt 
        << "   dr = " << dr 
        << "   dthe = " << dthe 
        << "   dphi = " << dphi
        << "   dthe2 = " << dthe2 
        << "   dphi2 = " << dphi2 << endl
        << "   rm > rad = " << rm 
        << "   sinthe = " << sinthe
        << "   sinthe2 = " << sinthe2
        << "   costhe = " << costhe
        << "   rmsinthe = " << rmsinthe
        << "   rm2sinthe = " << rm2sinthe
        << "   rm2sinthe2 = " << rm2sinthe2 << endl
        << "   Bathymetry = " << Bathymetry.y[j][k] << endl
        << "   dudthe = " << dudthe
        << "   dvdthe = " << dvdthe
        << "   dwdthe = " << dwdthe << endl
        << "   d2vdthe2 = " << d2vdthe2
        << "   d2wdthe2 = " << d2wdthe2 << endl
        << "   dpdthe = " << dpdthe << endl
        << "   dudphi = " << dudphi
        << "   dvdphi = " << dvdphi
        << "   dwdphi = " << dwdphi << endl
        << "   d2vdphi2 = " << d2vdphi2
        << "   d2wdphi2 = " << dwdphi << endl
        << "   dpdphi = " << dpdphi << endl
        << "   1. order the = " << - (v.x[im-1][j][k] * dvdthe/rm 
            + w.x[im-1][j][k] * dvdphi/rmsinthe) << endl
        << "   1.order phi = " << - (v.x[im-1][j][k] * dwdthe/rm 
            + w.x[im-1][j][k] * dwdphi/rmsinthe) << endl
        << "   2. order the = " << + (d2vdthe2/rm2 + dvdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_j * v.x[im-1][j][k] 
            + d2vdphi2/rm2sinthe2 + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re << endl
        << "   2.order phi = " << + (d2wdthe2/rm2 + dwdthe/rm2sinthe  * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_k * w.x[im-1][j][k] 
            + d2wdphi2/rm2sinthe2 + 2. * dudphi/rm2sinthe 
            + dvdphi * 2.0 * costhe/rm2sinthe2)/re << endl

        << "   uj-1 = " << u.x[im-1][j-1][k]
        << "   vj-1 = " << v.x[im-1][j-1][k]
        << "   wj-1 = " << w.x[im-1][j-1][k]
        << "   pj-1 = " << p_dyn.x[im-1][j-1][k] << endl

        << "   uj   = " << u.x[im-1][j][k]
        << "   vj   = " << v.x[im-1][j][k]
        << "   wj   = " << w.x[im-1][j][k]
        << "   pj   = " << p_dyn.x[im-1][j][k] << endl

        << "   uj+1 = " << u.x[im-1][j+1][k]
        << "   vj+1 = " << v.x[im-1][j+1][k]
        << "   wj+1 = " << w.x[im-1][j+1][k]
        << "   pj+1 = " << p_dyn.x[im-1][j+1][k] << endl

        << "   uk-1 = " << u.x[im-1][j][k-1]
        << "   vk-1 = " << v.x[im-1][j][k-1]
        << "   wk-1 = " << w.x[im-1][j][k-1]
        << "   pk-1 = " << p_dyn.x[im-1][j][k-1] << endl

        << "   uk   = " << u.x[im-1][j][k]
        << "   vk   = " << v.x[im-1][j][k]
        << "   wk   = " << w.x[im-1][j][k]
        << "   pk   = " << p_dyn.x[im-1][j][k] << endl

        << "   uk+1 = " << u.x[im-1][j][k+1]
        << "   vk+1 = " << v.x[im-1][j][k+1]
        << "   wk+1 = " << w.x[im-1][j][k+1]
        << "   pk+1 = " << p_dyn.x[im-1][j][k+1] << endl

        << "   un = " << un.x[im-1][j][k]
        << "   vn = " << vn.x[im-1][j][k]
        << "   wn = " << wn.x[im-1][j][k] << endl

        << "   rhs_v = " << rhs_v.x[im-1][j][k]
        << "   rhs_w = " << rhs_w.x[im-1][j][k] << endl

        << "   aux_v = " << aux_v.x[im-1][j][k]
        << "   aux_w = " << aux_w.x[im-1][j][k] << endl;
*/
    return;
}
