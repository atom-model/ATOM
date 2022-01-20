/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include <iostream>
#include <cmath>
#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

#define dxdr_a(X) \
    (h_d_i * (X->x[i+1][j][k] - X->x[i-1][j][k])/(2.0 * dr) * exp_rm)
#define dxdr_b(X) \
    (h_d_i * (- 3.0 * X->x[i][j][k] + 4.0 * X->x[i+1][j][k] - X->x[i+2][j][k])/(2.0 * dr) * exp_rm)
#define dxdr_d(X) \
    (h_d_i * (X->x[i+1][j][k] - X->x[i][j][k])/dr * exp_rm)

#define d2xdr2_a(X) \
    (h_d_i * (X->x[i+1][j][k] - 2.0 * X->x[i][j][k] + X->x[i-1][j][k])/dr2 * exp_2_rm)
#define d2xdr2_b(X) \
    (h_d_i * (X->x[i][j][k] - 2.0 * X->x[i+1][j][k] + X->x[i+2][j][k])/dr2 * exp_2_rm)

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
void cAtmosphereModel::RK_RHS_3D_Atmosphere(int i, int j, int k){
//  (Reinout vander Meulen, The immersed Boundary Method for the Incompressible Navier-Stokes Equations)
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double rm = rad.z[i];
    double rm2 = rm * rm;
    double exp_rm = 1.0/(rm + 1.0);
    double exp_2_rm = 1.0/pow((rm + 1.0),2);
    double sinthe = sin(the.z[j]);
    if(sinthe == 0.0) sinthe = 1.0e-5;
    double sinthe2 = sinthe * sinthe;
    double costhe = cos(the.z[j]);
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;
    double h_0_i = 0.0, h_d_i = 0.0;
    double topo_step = get_layer_height(i) - get_layer_height(i-1);
    double height = get_layer_height(i);
    double topo_diff = fabs(height - Topography.y[j][k]);
    if((topo_diff <= topo_step)
        &&((is_air(h, i, j, k)) 
        &&(is_land(h, i-1, j, k)))){
//        h_0_i = (0.5 * (acos(topo_diff * 3.14/L_atm) + 1.0));   // cosine distribution function
        h_d_i = 1.0 - topo_diff/topo_step;   // hat distribution function
        h_0_i = 0.0;
    }
    if(is_air(h, i, j, k)){
        h_0_i = 0.0;
        h_d_i = 1.0;
    }
    if(is_land(h, i, j, k)){
        h_0_i = 1.0;
        h_d_i = 0.0;
    }
    std::vector<Array*> arrays_1{&u, &v, &w, &t, &p_dyn, &c, &cloud, &ice, &co2};
    std::vector<Array*> arrays_2{&u, &v, &w, &t, &c, &cloud, &ice, &co2};
    enum array_index_1{i_u_1, i_v_1, i_w_1, i_t_1, i_p_1, i_c_1, 
        i_cloud_1, i_ice_1, i_co_1, last_array_index_1};
    enum array_index_2{i_u_2, i_v_2, i_w_2, i_t_2, i_c_2, i_cloud_2, 
        i_ice_2, i_co_2, last_array_index_2};
    std::vector<double> dxdr_vals(last_array_index_1), 
                        dxdthe_vals(last_array_index_1), 
                        dxdphi_vals(last_array_index_1),
                        d2xdr2_vals(last_array_index_2),
                        d2xdthe2_vals(last_array_index_2),
                        d2xdphi2_vals(last_array_index_2);
    bool r_flag = false, the_flag = false, phi_flag = false;
    if(i < im-2){
        if((is_land(h, i, j, k))
            &&(is_air(h, i+1, j, k))){
            for(int n=0; n<last_array_index_1; n++)
                dxdr_vals[n] = dxdr_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdr2_vals[n] = d2xdr2_b(arrays_2[n]);
            r_flag = true;
        }
    }
    if(i == im-2){
        if((is_land(h, i, j, k))
            &&(is_air(h, i+1, j, k))){
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
           dpdr = dxdr_vals[i_p_1], dcdr = dxdr_vals[i_c_1], 
           dclouddr = dxdr_vals[i_cloud_1], dicedr = dxdr_vals[i_ice_1],
           dcodr = dxdr_vals[i_co_1];
    double dudthe = dxdthe_vals[i_u_1], dvdthe = dxdthe_vals[i_v_1], 
           dwdthe = dxdthe_vals[i_w_1], dtdthe = dxdthe_vals[i_t_1],
           dpdthe = dxdthe_vals[i_p_1], dcdthe = dxdthe_vals[i_c_1], 
           dclouddthe = dxdthe_vals[i_cloud_1], dicedthe = dxdthe_vals[i_ice_1],
           dcodthe = dxdthe_vals[i_co_1];
    double dudphi = dxdphi_vals[i_u_1], dvdphi = dxdphi_vals[i_v_1], 
           dwdphi = dxdphi_vals[i_w_1], dtdphi = dxdphi_vals[i_t_1],
           dpdphi = dxdphi_vals[i_p_1], dcdphi = dxdphi_vals[i_c_1], 
           dclouddphi = dxdphi_vals[i_cloud_1], dicedphi = dxdphi_vals[i_ice_1],
           dcodphi = dxdphi_vals[i_co_1];
    double d2udr2 = d2xdr2_vals[i_u_2], d2vdr2 = d2xdr2_vals[i_v_2], 
           d2wdr2 = d2xdr2_vals[i_w_2], d2tdr2 = d2xdr2_vals[i_t_2],
           d2cdr2 = d2xdr2_vals[i_c_2], d2clouddr2 = d2xdr2_vals[i_cloud_2], 
           d2icedr2 = d2xdr2_vals[i_ice_2],
           d2codr2 = d2xdr2_vals[i_co_2];
    double d2udthe2 = d2xdthe2_vals[i_u_2], d2vdthe2 = d2xdthe2_vals[i_v_2], 
           d2wdthe2 = d2xdthe2_vals[i_w_2], d2tdthe2 = d2xdthe2_vals[i_t_2],
           d2cdthe2 = d2xdthe2_vals[i_c_2], d2clouddthe2 = d2xdthe2_vals[i_cloud_2], 
           d2icedthe2 = d2xdthe2_vals[i_ice_2],
           d2codthe2 = d2xdthe2_vals[i_co_2];
    double d2udphi2 = d2xdphi2_vals[i_u_2], d2vdphi2 = d2xdphi2_vals[i_v_2], 
           d2wdphi2 = d2xdphi2_vals[i_w_2], d2tdphi2 = d2xdphi2_vals[i_t_2],
           d2cdphi2 = d2xdphi2_vals[i_c_2], d2clouddphi2 = d2xdphi2_vals[i_cloud_2], 
           d2icedphi2 = d2xdphi2_vals[i_ice_2],
           d2codphi2 = d2xdphi2_vals[i_co_2];

    double coeff_g_p = 1.0; // coefficient allows the buoyancy term
    double coeff_energy_p = u_0 * u_0/(cp_l * t_0); // coefficient for the source terms = 2.33e-4 (Eckert-number)
    double coeff_energy = 1.0/(c_0 * cp_l * t_0); // coefficient for the source terms = 1.041e-4
    double coeff_buoy = r_air * (u_0 * u_0)/L_atm; // coefficient for bouancy term = 4.96e-3
    double coeff_trans = 1.0; // coefficient for the water vapour term = 1.
    double coeff_MC_vel = L_atm/(u_0 * u_0); // coefficient for MC_v and MC_w term = 250.0
    double coeff_MC_q = L_atm/(u_0 * c_0); // coefficient for MC_q term = 57142.9
    double coeff_MC_t = L_atm/(u_0 * t_0); // coefficient for MC_t term = 7.322
    double coriolis = 1.0;
    double coeff_Coriolis = r_air * u_0; // coefficient for Coriolis term = 239.28
    double coriolis_rad = h_d_i * coriolis * 2.0 * omega
        * costhe * w.x[i][j][k];
    double coriolis_the = - coriolis * 2.0 * omega
        * sinthe * w.x[i][j][k];
    double coriolis_phi = coriolis * 2.0 * omega
        * (sinthe * v.x[i][j][k] - costhe * u.x[i][j][k]);
    CoriolisForce.x[i][j][k] = coeff_Coriolis * sqrt((pow (coriolis_rad,2) 
        + pow (coriolis_the, 2) + pow (coriolis_phi,2))/3.0);
    BuoyancyForce.x[i][j][k] = buoyancy * r_air 
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
            + (dtdr * 2./rm 
            + dtdthe * costhe/rm2sinthe)/(re * pr)
            + coeff_MC_t * MC_t.x[i][j][k]
            + coeff_energy_p * (u.x[i][j][k] * dpdr 
            + v.x[i][j][k]/rm * dpdthe 
            + w.x[i][j][k]/rmsinthe * dpdphi)
            + coeff_energy * (S_c.x[i][j][k] + S_r.x[i][j][k]) 
                * lv * r_humid.x[i][j][k]
            + coeff_energy * (S_i.x[i][j][k] + S_s.x[i][j][k]) 
                * ls * r_humid.x[i][j][k]
            - h_0_i * t.x[i][j][k]/dr2 * exp_2_rm;
        rhs_u.x[i][j][k] = - (u.x[i][j][k] * dudr + v.x[i][j][k] * dudthe/rm 
            + w.x[i][j][k] * dudphi/rmsinthe) 
            + coeff_g_p * (t.x[i][j][k] - 1.0) * g
            - dpdr
            + (h_d_i * 2.0 * u.x[i][j][k]/rm2
            + 4.0 * dudr/rm + dudthe * costhe/rm2sinthe 
            )/re
            + coriolis_rad
            - h_0_i * u.x[i][j][k]/dr2 * exp_2_rm;
        rhs_v.x[i][j][k] = - (u.x[i][j][k] * dvdr + v.x[i][j][k] * dvdthe/rm
            + w.x[i][j][k] * dvdphi/rmsinthe)
            - dpdthe/rm
            + (dvdr * 2.0/rm
            + dvdthe/rm2sinthe * costhe
            - (1.0 + costhe/sinthe2)/rm2 * v.x[i][j][k] 
            + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re
            + coriolis_the
            + coeff_MC_vel * MC_v.x[i][j][k]
            - h_0_i * v.x[i][j][k]/dr2 * exp_2_rm;
        rhs_w.x[i][j][k] = - (u.x[i][j][k] * dwdr + v.x[i][j][k] * dwdthe/rm
            + w.x[i][j][k] * dwdphi/rmsinthe)
            - dpdphi/rmsinthe
            + (dwdr * 2.0/rm
            + dwdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * w.x[i][j][k]
            + 2.0 * dudphi/rm2sinthe 
            + dvdphi * 2.0 * costhe/rm2sinthe2
            )/re
            + coriolis_phi
            + coeff_MC_vel * MC_w.x[i][j][k]
            - h_0_i * w.x[i][j][k]/dr2 * exp_2_rm;
        rhs_c.x[i][j][k] = - (u.x[i][j][k] * dcdr + v.x[i][j][k] * dcdthe/rm
            + w.x[i][j][k] * dcdphi/rmsinthe) 
            + (dcdr * 2./rm 
            + dcdthe * costhe/rm2sinthe 
            )/(sc_WaterVapour * re)
            + coeff_trans * S_v.x[i][j][k] * r_humid.x[i][j][k]
            + coeff_MC_q * MC_q.x[i][j][k]
            - h_0_i * c.x[i][j][k]/dr2 * exp_2_rm;
        rhs_cloud.x[i][j][k] = - (u.x[i][j][k] * dclouddr 
            + v.x[i][j][k] * dclouddthe/rm
            + w.x[i][j][k] * dclouddphi/rmsinthe) 
            + (dclouddr * 2.0/rm
            + dclouddthe * costhe/rm2sinthe 
            )/(sc_WaterVapour * re)
            + coeff_trans * S_c.x[i][j][k] * r_humid.x[i][j][k]
            - h_0_i * cloud.x[i][j][k]/dr2 * exp_2_rm;
        rhs_ice.x[i][j][k] = - (u.x[i][j][k] * dicedr 
            + v.x[i][j][k] * dicedthe/rm
            + w.x[i][j][k] * dicedphi/rmsinthe) 
            + (dicedr * 2./rm
            + dicedthe * costhe/rm2sinthe 
            )/(sc_WaterVapour * re)
            + coeff_trans * S_i.x[i][j][k] * r_humid.x[i][j][k]
            - h_0_i * ice.x[i][j][k]/dr2 * exp_2_rm;
        rhs_co2.x[i][j][k] = - (u.x[i][j][k] * dcodr 
            + v.x[i][j][k] * dcodthe/rm
            + w.x[i][j][k] * dcodphi/rmsinthe) 
            + (dcodr * 2.0/rm
            + dcodthe * costhe/rm2sinthe 
            )/(sc_CO2 * re)
            - h_0_i * co2.x[i][j][k]/dr2 * exp_2_rm;
    }else{
        rhs_t.x[i][j][k] = - (u.x[i][j][k] * dtdr + v.x[i][j][k] * dtdthe/rm // Navier-Stokes equations
            + w.x[i][j][k] * dtdphi/rmsinthe) 
            + (d2tdr2 + dtdr * 2./rm + d2tdthe2/rm2
            + dtdthe * costhe/rm2sinthe + d2tdphi2/rm2sinthe2)/(re * pr)
            + coeff_MC_t * MC_t.x[i][j][k]
            + coeff_energy_p * (u.x[i][j][k] * dpdr 
            + v.x[i][j][k]/rm * dpdthe 
            + w.x[i][j][k]/rmsinthe * dpdphi)
            + coeff_energy * (S_c.x[i][j][k] + S_r.x[i][j][k]) 
                * lv * r_humid.x[i][j][k]
            + coeff_energy * (S_i.x[i][j][k] + S_s.x[i][j][k]) 
                * ls * r_humid.x[i][j][k]
            - h_0_i * t.x[i][j][k]/dr2 * exp_2_rm;
        rhs_u.x[i][j][k] = - (u.x[i][j][k] * dudr + v.x[i][j][k] * dudthe/rm 
            + w.x[i][j][k] * dudphi/rmsinthe) 
            + coeff_g_p * (t.x[i][j][k] - 1.0) * g
            - dpdr
            + (d2udr2 + h_d_i * 2.0 * u.x[i][j][k]/rm2 + d2udthe2/rm2
            + 4.0 * dudr/rm + dudthe * costhe/rm2sinthe 
            + d2udphi2/rm2sinthe2)/re
            + coriolis_rad
            - h_0_i * u.x[i][j][k]/dr2 * exp_2_rm;
        rhs_v.x[i][j][k] = - (u.x[i][j][k] * dvdr + v.x[i][j][k] * dvdthe/rm
            + w.x[i][j][k] * dvdphi/rmsinthe)
            - dpdthe/rm
            + (d2vdr2 + dvdr * 2.0/rm + d2vdthe2/rm2 
            + dvdthe/rm2sinthe * costhe
            - (1.0 + costhe/sinthe2)/rm2 * v.x[i][j][k] 
            + d2vdphi2/rm2sinthe2 + 2.0 * dudthe/rm2 
            - dwdphi * 2. * costhe/rm2sinthe2)/re
            + coriolis_the
            + coeff_MC_vel * MC_v.x[i][j][k]
            - h_0_i * v.x[i][j][k]/dr2 * exp_2_rm;
        rhs_w.x[i][j][k] = - (u.x[i][j][k] * dwdr + v.x[i][j][k] * dwdthe/rm
            + w.x[i][j][k] * dwdphi/rmsinthe)
            - dpdphi/rmsinthe
            + (d2wdr2 + dwdr * 2.0/rm + d2wdthe2/rm2
            + dwdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * w.x[i][j][k]
            + d2wdphi2/rm2sinthe2 + 2.0 * dudphi/rm2sinthe 
            + dvdphi * 2. * costhe/rm2sinthe2)/re
            + coriolis_phi
            + coeff_MC_vel * MC_w.x[i][j][k]
            - h_0_i * w.x[i][j][k]/dr2 * exp_2_rm;
        rhs_c.x[i][j][k] = - (u.x[i][j][k] * dcdr + v.x[i][j][k] * dcdthe/rm
            + w.x[i][j][k] * dcdphi/rmsinthe) 
            + (d2cdr2 + dcdr * 2./rm 
            + d2cdthe2/rm2
            + dcdthe * costhe/rm2sinthe 
            + d2cdphi2/rm2sinthe2)/(sc_WaterVapour * re)
            + coeff_trans * S_v.x[i][j][k] * r_humid.x[i][j][k]
            + coeff_MC_q * MC_q.x[i][j][k]
            - h_0_i * c.x[i][j][k]/dr2 * exp_2_rm;
        rhs_cloud.x[i][j][k] = - (u.x[i][j][k] * dclouddr 
            + v.x[i][j][k] * dclouddthe/rm
            + w.x[i][j][k] * dclouddphi/rmsinthe) 
            + (d2clouddr2 + dclouddr * 2.0/rm + d2clouddthe2/rm2
            + dclouddthe * costhe/rm2sinthe 
            + d2clouddphi2/rm2sinthe2)/(sc_WaterVapour * re)
            + coeff_trans * S_c.x[i][j][k] * r_humid.x[i][j][k]
            - h_0_i * cloud.x[i][j][k]/dr2 * exp_2_rm;
        rhs_ice.x[i][j][k] = - (u.x[i][j][k] * dicedr 
            + v.x[i][j][k] * dicedthe/rm
            + w.x[i][j][k] * dicedphi/rmsinthe) 
            + (d2icedr2 + dicedr * 2./rm + d2icedthe2/rm2
            + dicedthe * costhe/rm2sinthe 
            + d2icedphi2/rm2sinthe2)/(sc_WaterVapour * re)
            + coeff_trans * S_i.x[i][j][k] * r_humid.x[i][j][k]
            - h_0_i * ice.x[i][j][k]/dr2 * exp_2_rm;
        rhs_co2.x[i][j][k] = - (u.x[i][j][k] * dcodr 
            + v.x[i][j][k] * dcodthe/rm
            + w.x[i][j][k] * dcodphi/rmsinthe) 
            + (d2codr2 + dcodr * 2.0/rm + d2codthe2/rm2
            + dcodthe * costhe/rm2sinthe 
            + d2codphi2/rm2sinthe2)/(sc_CO2 * re)
            - h_0_i * co2.x[i][j][k]/dr2 * exp_2_rm;
    }
    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe/rm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi/rmsinthe;
    if(is_land(h, i, j, k)){
        rhs_u.x[i][j][k] = 0.0;
        rhs_v.x[i][j][k] = 0.0;
        rhs_w.x[i][j][k] = 0.0;
        rhs_t.x[i][j][k] = 0.0;
        rhs_c.x[i][j][k] = 0.0;
        rhs_cloud.x[i][j][k] = 0.0;
        rhs_ice.x[i][j][k] = 0.0;
        rhs_co2.x[i][j][k] = 0.0;
        aux_u.x[i][j][k] = 0.0;
        aux_v.x[i][j][k] = 0.0;
        aux_w.x[i][j][k] = 0.0;
    }
/*
    cout.precision(8);
//    if((j == 75) &&(k == 180)) cout << "Atmosphere code" << endl
    if((j == 1) &&(k == 180)) cout << "Atmosphere code" << endl
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
        << "   Topography = " << Topography.y[j][k]
        << "   topo_diff = " << topo_diff << endl
        << "   coriolis = " << coriolis
        << "   omega = " << omega
        << "   coriolis_rad = " << coriolis_rad
        << "   coriolis_the = " << coriolis_the
        << "   coriolis_phi = " << coriolis_phi
        << "   CoriolisForce = " << CoriolisForce.x[i][j][k] << endl
        << "   h_0_i = " << h_0_i
        << "   h_d_i = " << h_d_i << endl
//        << "   residuum = " << residuum << endl
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
        << "   rhs_c = " << rhs_t.x[i][j][k] << endl

        << "   aux_t = " << rhs_t.x[i][j][k]
        << "   aux_u = " << rhs_u.x[i][j][k]
        << "   aux_v = " << rhs_v.x[i][j][k]
        << "   aux_w = " << rhs_w.x[i][j][k]
        << "   aux_c = " << rhs_t.x[i][j][k] << endl;
*/
    return;
}
/*
*
*/
void cAtmosphereModel::RK_RHS_2D_Atmosphere(int i, int j, int k){
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
    if(is_air(h, 0, j, k)){
        h_0_j = 0.0;
        h_0_k = 0.0;
        h_d_j = 1.0;
        h_d_k = 1.0;
    }
    if(is_land(h, 0, j, k)){
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
        if((is_land(h, 0, j, k))
            &&((is_air(h, 0, j+1, k)) 
            &&(is_air(h, 0, j+2, k)))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_b(arrays_2[n]);
            the_flag = true;
        }
        if((is_land(h, 0, j, k))
            &&(is_air(h, 0, j-1, k)) 
            &&(is_air(h, 0, j-2, k))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_c(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_c(arrays_2[n]);
            the_flag = true;
        }
        if(((is_land(h, 0, j, k)) 
            &&((is_air(h, 0, j+1, k))
            &&(is_land(h, 0, j+2, k))))||((j == jm-2)
            &&((is_air(h, 0, j, k))
            &&(is_land(h, 0, j+1, k))))){
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_b(arrays_2[n]);
            the_flag = true;
        }
        if(((is_land(h, 0, j, k)) 
            &&((is_air(h, 0, j-1, k))
            &&(is_land(h, 0, j-2, k)))) 
            ||((j == 1)&&((is_land(h, 0, j, k)) 
            &&(is_air(h, 0, j-1, k))))){
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_c(arrays_2[n]);
            the_flag = true;
        }
    }
    if((k >= 1)&&(k < km-2)){
        if((is_land(h, 0, j, k))
            &&(is_air(h, 0, j, k+1)) 
            &&(is_air(h, 0, j, k+2))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_b(arrays_2[n]);
            phi_flag = true;
        }
        if((is_land(h, 0, j, k))
            &&(is_air(h, 0, j, k-1)) 
            &&(is_air(h, 0, j, k-2))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_c(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_c(arrays_2[n]);
            phi_flag = true;
        }
        if(((is_land(h, 0, j, k))
            &&((is_air(h, 0, j, k+1)) 
            &&(is_land(h, 0, j, k+2))))||((k == km-2)
            &&((is_air(h, 0, j, k))
            &&(is_land(h, 0, j, k+1))))){
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_b(arrays_2[n]);
            phi_flag = true;
        }
        if(((is_land(h, 0, j, k)) 
            &&((is_air(h, 0, j, k-1))
            &&(is_land(h, 0, j, k-2)))) 
            ||((k == 1)&&((is_land(h, 0, j, k)) 
            &&(is_air(h, 0, j, k-1))))){
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
        rhs_v.x[0][j][k] = - (v.x[0][j][k] * dvdthe/rm 
            + w.x[0][j][k] * dvdphi/rmsinthe) 
            - dpdthe/rm 
            + (dvdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_j * v.x[0][j][k] 
            + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re
            - h_0_j * v.x[0][j][k]/dthe2;                  // RHS v-velocity Euler equations without friction terms
        rhs_w.x[0][j][k] = - (v.x[0][j][k] * dwdthe/rm 
            + w.x[0][j][k] * dwdphi/rmsinthe) 
            - dpdphi/rmsinthe 
            + (dwdthe/rm2sinthe  * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_k * w.x[0][j][k] 
            + 2.0 * dudphi/rm2sinthe 
            + dvdphi * 2.0 * costhe/rm2sinthe2)/re
            - h_0_k * w.x[0][j][k]/dphi2;                  // RHS w-velocity Euler equations without friction terms
    }else{
        rhs_v.x[0][j][k] = - (v.x[0][j][k] * dvdthe/rm 
            + w.x[0][j][k] * dvdphi/rmsinthe) 
            - dpdthe/rm 
            + (d2vdthe2/rm2 + dvdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_j * v.x[0][j][k] 
            + d2vdphi2/rm2sinthe2 + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re
            - h_0_j * v.x[0][j][k]/dthe2;                  // RHS v-velocity Navier-Stokes equations with friction terms

        rhs_w.x[0][j][k] = - (v.x[0][j][k] * dwdthe/rm 
            + w.x[0][j][k] * dwdphi/rmsinthe) 
            - dpdphi/rmsinthe 
            + (d2wdthe2/rm2 + dwdthe/rm2sinthe  * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_k * w.x[0][j][k] 
            + d2wdphi2/rm2sinthe2 + 2.0 * dudphi/rm2sinthe 
            + dvdphi * 2.0 * costhe/rm2sinthe2)/re
            - h_0_k * w.x[0][j][k]/dphi2;                  // RHS w-velocity Navier-Stokes equations with friction terms
    }
 
    aux_v.x[0][j][k] = rhs_v.x[0][j][k] + dpdthe/rm;
    aux_w.x[0][j][k] = rhs_w.x[0][j][k] + dpdphi/rmsinthe;

    if(is_land(h, 0, j, k)){
        aux_u.x[0][j][k] = 0.0;
        aux_v.x[0][j][k] = 0.0;
        aux_w.x[0][j][k] = 0.0;
    }

    cout.precision(8);
    if((j == 90) &&(k == 180)) cout << "Atmosphere code" << endl
        << "northern hemisphere" << endl
        << "   0 = " << 0 << "   j = " << j << "   k = " << k  << endl
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
        << "   Topography = " << Topography.y[j][k] << endl
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
        << "   1. order the = " << - (v.x[0][j][k] * dvdthe/rm 
            + w.x[0][j][k] * dvdphi/rmsinthe) << endl
        << "   1.order phi = " << - (v.x[0][j][k] * dwdthe/rm 
            + w.x[0][j][k] * dwdphi/rmsinthe) << endl
        << "   2. order the = " << + (d2vdthe2/rm2 + dvdthe/rm2sinthe * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_j * v.x[0][j][k] 
            + d2vdphi2/rm2sinthe2 + 2.0 * dudthe/rm2 
            - dwdphi * 2.0 * costhe/rm2sinthe2)/re << endl
        << "   2.order phi = " << + (d2wdthe2/rm2 + dwdthe/rm2sinthe  * costhe 
            - (1.0 + costhe/sinthe2)/rm2 * h_d_k * w.x[0][j][k] 
            + d2wdphi2/rm2sinthe2 + 2. * dudphi/rm2sinthe 
            + dvdphi * 2.0 * costhe/rm2sinthe2)/re << endl

        << "   uj-1 = " << u.x[0][j-1][k]
        << "   vj-1 = " << v.x[0][j-1][k]
        << "   wj-1 = " << w.x[0][j-1][k]
        << "   pj-1 = " << p_dyn.x[0][j-1][k] << endl

        << "   uj   = " << u.x[0][j][k]
        << "   vj   = " << v.x[0][j][k]
        << "   wj   = " << w.x[0][j][k]
        << "   pj   = " << p_dyn.x[0][j][k] << endl

        << "   uj+1 = " << u.x[0][j+1][k]
        << "   vj+1 = " << v.x[0][j+1][k]
        << "   wj+1 = " << w.x[0][j+1][k]
        << "   pj+1 = " << p_dyn.x[0][j+1][k] << endl

        << "   uk-1 = " << u.x[0][j][k-1]
        << "   vk-1 = " << v.x[0][j][k-1]
        << "   wk-1 = " << w.x[0][j][k-1]
        << "   pk-1 = " << p_dyn.x[0][j][k-1] << endl

        << "   uk   = " << u.x[0][j][k]
        << "   vk   = " << v.x[0][j][k]
        << "   wk   = " << w.x[0][j][k]
        << "   pk   = " << p_dyn.x[0][j][k] << endl

        << "   uk+1 = " << u.x[0][j][k+1]
        << "   vk+1 = " << v.x[0][j][k+1]
        << "   wk+1 = " << w.x[0][j][k+1]
        << "   pk+1 = " << p_dyn.x[0][j][k+1] << endl

        << "   un = " << un.x[0][j][k]
        << "   vn = " << vn.x[0][j][k]
        << "   wn = " << wn.x[0][j][k] << endl

        << "   rhs_v = " << rhs_v.x[0][j][k]
        << "   rhs_w = " << rhs_w.x[0][j][k] << endl

        << "   aux_v = " << aux_v.x[0][j][k]
        << "   aux_w = " << aux_w.x[0][j][k] << endl << endl;

    return;
}


