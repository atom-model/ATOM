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
    (h_d_i * (X->x[i+1][j][k] - X->x[i-1][j][k])/(2. * dr) * exp_rm)
#define dxdr_b(X) \
    (h_d_i * (- 3. * X->x[i][j][k] + 4. * X->x[i+1][j][k] - X->x[i+2][j][k])/(2. * dr) * exp_rm)
#define d2xdr2_a(X) \
    (h_d_i * (X->x[i+1][j][k] - 2. * X->x[i][j][k] + X->x[i-1][j][k])/dr2 * exp_2_rm)
#define d2xdr2_b(X) \
    (h_d_i * (X->x[i][j][k] - 2. * X->x[i+1][j][k] + X->x[i+2][j][k])/dr2 * exp_2_rm)
#define dxdthe_a(X) \
    (h_d_j * (X->x[i][j+1][k] - X->x[i][j-1][k])/(2. * dthe))
#define dxdthe_b(X) \
    (h_d_j * (- 3. * X->x[i][j][k] + 4. * X->x[i][j+1][k] - X->x[i][j+2][k])/(2. * dthe))
#define dxdthe_c(X) \
    (h_d_j * (- 3. * X->x[i][j][k] + 4. * X->x[i][j-1][k] - X->x[i][j-2][k])/(2. * dthe))
#define dxdthe_d(X) \
    (h_d_j * (X->x[i][j+1][k] - X->x[i][j][k])/dthe)
#define dxdthe_e(X) \
    (h_d_j * (X->x[i][j-1][k] - X->x[i][j][k])/dthe)
#define d2xdthe2_a(X) \
    (h_d_j * (X->x[i][j+1][k] - 2. * X->x[i][j][k] + X->x[i][j-1][k])/dthe2)
#define d2xdthe2_b(X) \
    (h_d_j * (X->x[i][j][k] - 2. * X->x[i][j+1][k] + X->x[i][j+2][k])/dthe2)
#define d2xdthe2_c(X) \
    (h_d_j * (X->x[i][j][k] - 2. * X->x[i][j-1][k] + X->x[i][j-2][k])/dthe2)
#define dxdphi_a(X) \
    (h_d_k * (X->x[i][j][k+1] - X->x[i][j][k-1])/(2. * dphi))
#define dxdphi_b(X) \
    (h_d_k * (- 3. * X->x[i][j][k] + 4. * X->x[i][j][k+1] - X->x[i][j][k+2])/(2. * dphi))
#define dxdphi_c(X) \
    (h_d_k * (- 3. * X->x[i][j][k] + 4. * X->x[i][j][k-1] - X->x[i][j][k-2])/(2. * dphi))
#define dxdphi_d(X) \
    (h_d_k * (X->x[i][j][k+1] - X->x[i][j][k])/dphi)
#define dxdphi_e(X) \
    (h_d_k * (X->x[i][j][k-1] - X->x[i][j][k])/dphi)
#define d2xdphi2_a(X) \
    (h_d_k * (X->x[i][j][k+1] - 2. * X->x[i][j][k] + X->x[i][j][k-1])/dphi2)
#define d2xdphi2_b(X) \
    (h_d_k * (X->x[i][j][k] - 2. * X->x[i][j][k+1] + X->x[i][j][k+2])/dphi2)
#define d2xdphi2_c(X) \
    (h_d_k * (X->x[i][j][k] - 2. * X->x[i][j][k-1] + X->x[i][j][k-2])/dphi2)
/*
*
*/
void cAtmosphereModel::RK_RHS_3D_Atmosphere(int i, int j, int k){
    rad.Coordinates(im, r0, dr);
    double cc = - 1.0;  // factor leads to better results 
//  (Reinout vander Meulen, The immersed Boundary Method for the Incompressible Navier-Stokes Equations)
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double rm = rad.z[i];
    double rm2 = rm * rm;
    double exp_rm = 1./(rm + 1.0);
    double exp_2_rm = 1./pow((rm + 1.0),2);
    double sinthe = sin(the.z[j]);
    double sinthe2 = sinthe * sinthe;
    double costhe = cos(the.z[j]);
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;
    double dist = 0.0;
    double h_0_i = 0.0, h_d_i = 0.0;
    double h_0_j = 0.0, h_d_j = 0.0;
    double h_0_k = 0.0, h_d_k = 0.0;
    double topo_step = get_layer_height(i) - get_layer_height(i-1);
    double height = get_layer_height(i);
    double dist_coeff = 0.0;
//    double dist_coeff = 0.5;
    double topo_diff = fabs(height - Topography.y[j][k]);
    if((topo_diff <= topo_step)&&((is_air(h, i, j, k)) 
        &&(is_land(h, i-1, j, k)))){
        double h_0_i = (0.5 * (acos(topo_diff * 3.14/L_atm) + 1.0));   // cosine distribution function, better results for benchmark case
//        h_0_i = topo_diff/topo_step;  // hat distribution function
        h_d_i = 1.0 - h_0_i; 
    }else{
        h_0_i = 0.;
        h_d_i = 1.0; 
    }
    if((is_air(h, i, j, k))&&(is_land(h, i, j+1, k))){ 
        dist = dist_coeff * dthe;
        h_0_j = dist/dthe;
        h_d_j = 1. - h_0_j; 
    }
    if((is_air(h, i, j, k))&&(is_land(h, i, j-1, k))){
        dist = dist_coeff * dthe;
        h_0_j = dist/dthe;
        h_d_j = 1. - h_0_j; 
    }
    if((is_air(h, i, j, k))&&(is_land(h, i, j, k+1))){
        dist = dist_coeff * dphi;
        h_0_k = dist/dphi;
        h_d_k = 1. - h_0_k; 
    }
    if((is_air(h, i, j, k))&&(is_land(h, i, j, k-1))){
        dist = dist_coeff * dphi;
        h_0_k = dist/dphi;
        h_d_k = 1. - h_0_k; 
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
    bool r_flag = false , the_flag = false, phi_flag = false;
    if(i < im - 2){
        if((is_land(h, i, j, k))&&(is_air(h, i+1, j, k))){
            for(int n=0; n<last_array_index_1; n++)
                dxdr_vals[n] = dxdr_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdr2_vals[n] = d2xdr2_b(arrays_2[n]);
            r_flag = true;
        }
    }
    if((j >= 2)&&(j < jm - 3)){
        if((is_land(h, i, j, k))&&((is_air(h, i, j+1, k)) 
           &&(is_air(h, i, j+2, k)))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_b(arrays_2[n]);
            the_flag = true;
        }
        if((is_land(h, i, j, k))&&(is_air(h, i, j-1, k)) 
           &&(is_air(h, i, j-2, k))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_c(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = d2xdthe2_c(arrays_2[n]);
            the_flag = true;
        }
        if(((is_land(h, i, j, k)) 
               &&((is_air(h, i, j+1, k))&&(is_land(h, i, j+2, k)))) 
               ||((j == jm - 2)
               &&((is_air(h, i, j, k))&&(is_land(h, i, j+1, k))))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_d(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = 0.;
            the_flag = true;
        }
        if(((is_land(h, i, j, k)) 
               &&((is_air(h, i, j-1, k))&&(is_land(h, i, j-2, k)))) 
               ||((j == 1)&&((is_land(h, i, j, k)) 
               &&(is_air(h, i, j-1, k))))){
            for(int n=0; n<last_array_index_1; n++)
                dxdthe_vals[n] = dxdthe_e(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdthe2_vals[n] = 0.;
            the_flag = true;
        }
    }
    if((k >= 2)&&(k < km - 3)){
        if((is_land(h, i, j, k))&&(is_air(h, i, j, k+1)) 
           &&(is_air(h, i, j, k+2))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_b(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_b(arrays_2[n]);
            phi_flag = true;
        }
        if((is_land(h, i, j, k))&&(is_air(h, i, j, k-1)) 
           &&(is_air(h, i, j, k-2))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_c(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = d2xdphi2_c(arrays_2[n]);
            phi_flag = true;
        }
        if(((is_land(h, i, j, k))&&((is_air(h, i, j, k+1)) 
           &&(is_land(h, i, j, k+2))))||((k == km - 2)
           &&((is_air(h, i, j, k))&&(is_land(h, i, j, k+1))))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_d(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = 0.;
            phi_flag = true;
        }
        if(((is_land(h, i, j, k)) 
               &&((is_air(h, i, j, k-1))&&(is_land(h, i, j, k-2)))) 
               ||((k == 1)&&((is_land(h, i, j, k)) 
               &&(is_air(h, i, j, k-1))))){
            for(int n=0; n<last_array_index_1; n++)
                dxdphi_vals[n] = dxdphi_e(arrays_1[n]);
            for(int n=0; n<last_array_index_2; n++)
                d2xdphi2_vals[n] = 0.;
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
    double dudr = dxdr_vals[i_u_1], dvdr = dxdr_vals[i_v_1], dwdr = dxdr_vals[i_w_1], dtdr = dxdr_vals[i_t_1],
           dpdr = dxdr_vals[i_p_1], dcdr = dxdr_vals[i_c_1], dclouddr = dxdr_vals[i_cloud_1], dicedr = dxdr_vals[i_ice_1],
           dcodr = dxdr_vals[i_co_1];
    double dudthe = dxdthe_vals[i_u_1], dvdthe = dxdthe_vals[i_v_1], dwdthe = dxdthe_vals[i_w_1], dtdthe = dxdthe_vals[i_t_1],
           dpdthe = dxdthe_vals[i_p_1], dcdthe = dxdthe_vals[i_c_1], dclouddthe = dxdthe_vals[i_cloud_1], dicedthe = dxdthe_vals[i_ice_1],
           dcodthe = dxdthe_vals[i_co_1];
    double dudphi = dxdphi_vals[i_u_1], dvdphi = dxdphi_vals[i_v_1], dwdphi = dxdphi_vals[i_w_1], dtdphi = dxdphi_vals[i_t_1],
           dpdphi = dxdphi_vals[i_p_1], dcdphi = dxdphi_vals[i_c_1], dclouddphi = dxdphi_vals[i_cloud_1], dicedphi = dxdphi_vals[i_ice_1],
           dcodphi = dxdphi_vals[i_co_1];
    double d2udr2 = d2xdr2_vals[i_u_2], d2vdr2 = d2xdr2_vals[i_v_2], d2wdr2 = d2xdr2_vals[i_w_2], d2tdr2 = d2xdr2_vals[i_t_2],
           d2cdr2 = d2xdr2_vals[i_c_2], d2clouddr2 = d2xdr2_vals[i_cloud_2], d2icedr2 = d2xdr2_vals[i_ice_2],
           d2codr2 = d2xdr2_vals[i_co_2];
    double d2udthe2 = d2xdthe2_vals[i_u_2], d2vdthe2 = d2xdthe2_vals[i_v_2], d2wdthe2 = d2xdthe2_vals[i_w_2], d2tdthe2 = d2xdthe2_vals[i_t_2],
           d2cdthe2 = d2xdthe2_vals[i_c_2], d2clouddthe2 = d2xdthe2_vals[i_cloud_2], d2icedthe2 = d2xdthe2_vals[i_ice_2],
           d2codthe2 = d2xdthe2_vals[i_co_2];
    double d2udphi2 = d2xdphi2_vals[i_u_2], d2vdphi2 = d2xdphi2_vals[i_v_2], d2wdphi2 = d2xdphi2_vals[i_w_2], d2tdphi2 = d2xdphi2_vals[i_t_2],
           d2cdphi2 = d2xdphi2_vals[i_c_2], d2clouddphi2 = d2xdphi2_vals[i_cloud_2], d2icedphi2 = d2xdphi2_vals[i_ice_2],
           d2codphi2 = d2xdphi2_vals[i_co_2];

    double coeff_p = 1./r_air; // coefficient allows the dynamic pressure term
    double coeff_g_p = 1.; // coefficient allows the buoyancy term
    double coeff_energy_p = u_0 * u_0/(cp_l * t_0); // coefficient for the source terms = 2.33e-4 (Eckert-number)
    double coeff_energy = 1./(c_0 * cp_l * t_0); // coefficient for the source terms = 1.041e-4
    double coeff_buoy = r_air * (u_0 * u_0)/L_atm; // coefficient for bouancy term = 4.96e-3
    double coeff_trans = 1.; // coefficient for the water vapour term = 1.
    double coeff_MC_vel = L_atm/(u_0 * u_0); // coefficient for MC_v and MC_w term = 250.0
    double coeff_MC_q = L_atm/(u_0 * c_0); // coefficient for MC_q term = 57142.8
    double coeff_MC_t = L_atm/(u_0 * t_0); // coefficient for MC_t term = 7.322
    double coriolis = 1.;
    double coriolis_rad = coriolis * 2. * omega
        * costhe * w.x[i][j][k];
    double coriolis_the = - coriolis * 2. * omega
        * sinthe * w.x[i][j][k];
    double coriolis_phi = coriolis * 2. * omega
        * (sinthe * v.x[i][j][k] - costhe * u.x[i][j][k]);
    CoriolisForce.x[i][j][k] = sqrt((pow (coriolis_rad, 2) 
        + pow (coriolis_the, 2) + pow (coriolis_phi, 2))/3.);
    BuoyancyForce.x[i][j][k] = buoyancy * coeff_buoy * t.x[i][j][k]/t_0 * g;
    PressureGradientForce.x[i][j][k] = - coeff_buoy * coeff_p * dpdr;

    rhs_t.x[i][j][k] = - (u.x[i][j][k] * dtdr + v.x[i][j][k] * dtdthe/rm
        + w.x[i][j][k] * dtdphi/rmsinthe) + (d2tdr2 + dtdr * 2./rm + d2tdthe2/rm2
        + dtdthe * costhe/rm2sinthe + d2tdphi2/rm2sinthe2)/(re * pr)
        + coeff_MC_t * MC_t.x[i][j][k]
        + coeff_energy_p * (u.x[i][j][k] * dpdr 
        + v.x[i][j][k]/rm * dpdthe 
        + w.x[i][j][k]/rmsinthe * dpdphi)
        + coeff_energy * (S_c.x[i][j][k] + S_r.x[i][j][k]) 
            * lv * r_humid.x[i][j][k]
        + coeff_energy * (S_i.x[i][j][k] + S_s.x[i][j][k]) 
            * ls * r_humid.x[i][j][k]
        + cc * h_0_i * t.x[i][j][k]/dr2 * exp_2_rm;
    rhs_u.x[i][j][k] = - (u.x[i][j][k] * dudr + v.x[i][j][k] * dudthe/rm 
        + w.x[i][j][k] * dudphi/rmsinthe) 
        + coeff_g_p * t.x[i][j][k] * 1.
        - coeff_p * dpdr
        + (d2udr2 + h_d_i * 2. * u.x[i][j][k]/rm2 + d2udthe2/rm2
        + 4. * dudr/rm + dudthe * costhe/rm2sinthe + d2udphi2/rm2sinthe2)/re
        + coriolis_rad
        + cc * h_0_i * u.x[i][j][k]/dr2 * exp_2_rm;
    rhs_v.x[i][j][k] = - (u.x[i][j][k] * dvdr + v.x[i][j][k] * dvdthe/rm
        + w.x[i][j][k] * dvdphi/rmsinthe)
        - coeff_p * dpdthe/rm
        + (d2vdr2 + dvdr * 2./rm + d2vdthe2/rm2 
        + dvdthe/rm2sinthe * costhe
        - (1. + costhe * costhe/(rm * sinthe2)) * h_d_j * v.x[i][j][k] 
        + d2vdphi2/rm2sinthe2 + 2. * dudthe/rm2 
        - dwdphi * 2. * costhe/rm2sinthe2)/re
        + coriolis_the
        + coeff_MC_vel * MC_v.x[i][j][k]
        + cc * h_0_i * v.x[i][j][k]/dr2 * exp_2_rm;
    rhs_w.x[i][j][k] = - (u.x[i][j][k] * dwdr + v.x[i][j][k] * dwdthe/rm
        + w.x[i][j][k] * dwdphi/rmsinthe)
        - coeff_p * dpdphi/rmsinthe
        + (d2wdr2 + dwdr * 2./rm + d2wdthe2/rm2
        + dwdthe/rm2sinthe  * costhe - (1. + costhe 
        * costhe/(rm * sinthe2)) * h_d_k * w.x[i][j][k]
        + d2wdphi2/rm2sinthe2 + 2. * dudphi/rm2sinthe 
        + dvdphi * 2. * costhe/rm2sinthe2)/re
        + coriolis_phi
        + coeff_MC_vel * MC_w.x[i][j][k]
        + cc * h_0_i * w.x[i][j][k]/dr2 * exp_2_rm;
    rhs_c.x[i][j][k] = - (u.x[i][j][k] * dcdr + v.x[i][j][k] * dcdthe/rm
        + w.x[i][j][k] * dcdphi/rmsinthe) 
        + (d2cdr2 + dcdr * 2./rm + d2cdthe2/rm2
        + dcdthe * costhe/rm2sinthe 
        + d2cdphi2/rm2sinthe2)/(sc_WaterVapour * re)
        + coeff_trans * S_v.x[i][j][k] * r_humid.x[i][j][k]
        + coeff_MC_q * MC_q.x[i][j][k]
        + cc * h_0_i * c.x[i][j][k]/dr2 * exp_2_rm;
    rhs_cloud.x[i][j][k] = - (u.x[i][j][k] * dclouddr 
        + v.x[i][j][k] * dclouddthe/rm
        + w.x[i][j][k] * dclouddphi/rmsinthe) 
        + (d2clouddr2 + dclouddr * 2./rm + d2clouddthe2/rm2
        + dclouddthe * costhe/rm2sinthe 
        + d2clouddphi2/rm2sinthe2)/(sc_WaterVapour * re)
        + coeff_trans * S_c.x[i][j][k] * r_humid.x[i][j][k]
        + cc * h_0_i * cloud.x[i][j][k]/dr2 * exp_rm;
    rhs_ice.x[i][j][k] = - (u.x[i][j][k] * dicedr 
        + v.x[i][j][k] * dicedthe/rm
        + w.x[i][j][k] * dicedphi/rmsinthe) 
        + (d2icedr2 + dicedr * 2./rm + d2icedthe2/rm2
        + dicedthe * costhe/rm2sinthe 
        + d2icedphi2/rm2sinthe2)/(sc_WaterVapour * re)
        + coeff_trans * S_i.x[i][j][k] * r_humid.x[i][j][k]
        + cc * h_0_i * ice.x[i][j][k]/dr2 * exp_2_rm;
    rhs_co2.x[i][j][k] = - (u.x[i][j][k] * dcodr 
        + v.x[i][j][k] * dcodthe/rm
        + w.x[i][j][k] * dcodphi/rmsinthe) 
        + (d2codr2 + dcodr * 2./rm + d2codthe2/rm2
        + dcodthe * costhe/rm2sinthe 
        + d2codphi2/rm2sinthe2)/(sc_CO2 * re)
        + cc * h_0_i * co2.x[i][j][k]/dr2 * exp_2_rm;
    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + coeff_p * dpdr;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + coeff_p * dpdthe/rm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + coeff_p * dpdphi/rmsinthe;
    if(is_land(h, i, j, k)){
        rhs_u.x[i][j][k] = rhs_v.x[i][j][k] = rhs_w.x[i][j][k] = 
            rhs_t.x[i][j][k] = rhs_c.x[i][j][k] = rhs_cloud.x[i][j][k] = 
            rhs_ice.x[i][j][k] = rhs_co2.x[i][j][k] = 0.;
        aux_u.x[i][j][k] = aux_v.x[i][j][k] = aux_w.x[i][j][k] = 0.;
    }
/*
    cout.precision(8);
    if((j == 75) &&(k == 180)) cout << "Atmosphere code" << endl
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
        << "   rhs_c = " << rhs_t.x[i][j][k] << endl;
*/
    return;
}
/*
*
*/
void cAtmosphereModel::RK_RHS_2D_Atmosphere(int j, int k){
    rad.Coordinates(im, r0, dr);
    double cc = - 1.;  // factor leads to better results (adapted method)
//  (Reinout vander Meulen, The immersed Boundary Method for the Incompressible Navier-Stokes Equations)
//    double coeff_p = 100. * p_0/(r_air*u_0*u_0);
    double coeff_p = L_atm/(r_air * u_0 * u_0); // coefficient for the pressure terms = 207.624
//    double coeff_p = 1.;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
//    double rm = log(get_layer_height(0) + 1.);
//    double rm = log(get_layer_height(0)/L_atm + 1.);
    double rm = log(rad.z[0] + 1.);
    double rm2 = rm * rm;
    double sinthe = sin(the.z[j]);
    double sinthe2 = sinthe * sinthe;
    double costhe = cos(the.z[j]);
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;
    double dist = 0;
    double h_0_j = 0, h_d_j = 0;
    double h_0_k = 0, h_d_k = 0;
    if(is_air(h, 0, j, k)){
        h_0_j = h_0_k = 0.;
        h_d_j = h_d_k = 1.; 
    }
    if(is_land(h, 0, j, k)){
        h_0_j = h_0_k = 1.;
        h_d_j = h_d_k = 0.; 
    }
    double dist_coeff = .5;
    if((is_air(h, 0, j, k))&&(is_land(h, 0, j+1, k))){ 
        dist = dist_coeff * dthe;
        h_0_j = dist/dthe;
        h_d_j = 1. - h_0_j; 
    }
    if((is_air(h, 0, j, k))&&(is_land(h, 0, j-1, k))){
        dist = dist_coeff * dthe;
        h_0_j = dist/dthe;
        h_d_j = 1. - h_0_j; 
    }
    if((is_air(h, 0, j, k))&&(is_land(h, 0, j, k+1))){
        dist = dist_coeff * dphi;
        h_0_k = dist/dphi;
        h_d_k = 1. - h_0_k; 
    }
    if((is_air(h, 0, j, k))&&(is_land(h, 0, j, k-1))){
        dist = dist_coeff * dphi;
        h_0_k = dist/dphi;
        h_d_k = 1. - h_0_k; 
    }
    double dvdthe = h_d_j * (v.x[0][j+1][k] - v.x[0][j-1][k]) 
        /(2. * dthe);
    double dwdthe = h_d_j * (w.x[0][j+1][k] - w.x[0][j-1][k]) 
        /(2. * dthe);
    double dpdthe = h_d_j * (p_dyn.x[0][j+1][k] - p_dyn.x[0][j-1][k]) 
        /(2. * dthe);
    double dvdphi = h_d_k * (v.x[0][j][k+1] - v.x[0][j][k-1]) 
        /(2. * dphi);
    double dwdphi = h_d_k * (w.x[0][j][k+1] - w.x[0][j][k-1]) 
        /(2. * dphi);
    double dpdphi = h_d_k * (p_dyn.x[0][j][k+1] - p_dyn.x[0][j][k-1]) 
        /(2. * dphi);
    double d2vdthe2 = h_d_j *  (v.x[0][j+1][k] - 2. * v.x[0][j][k] 
        + v.x[0][j-1][k])/dthe2;
    double d2wdthe2 = h_d_j * (w.x[0][j+1][k] - 2. * w.x[0][j][k] 
        + w.x[0][j-1][k])/dthe2;
    double d2vdphi2 = h_d_k * (v.x[0][j][k+1] - 2. * v.x[0][j][k] 
        + v.x[0][j][k-1])/dphi2;
    double d2wdphi2 = h_d_k * (w.x[0][j][k+1] - 2. * w.x[0][j][k] 
        + w.x[0][j][k-1])/dphi2;
    if((j >= 2)&&(j < jm-3)){
        if((is_land(h, 0, j, k))&&((is_air(h, 0, j+1, k)&&(is_air(h, 0, j+2, k))))){
            dvdthe = h_d_j * (- 3. * v.x[0][j][k] + 4. * v.x[0][j+1][k] - v.x[0][j+2][k]) 
                       /(2. * dthe);
            dwdthe = h_d_j * (- 3. * w.x[0][j][k] + 4. * w.x[0][j+1][k] - w.x[0][j+2][k]) 
                       /(2. * dthe);
            dpdthe = h_d_j * (- 3. * p_dyn.x[0][j][k] + 4. * p_dyn.x[0][j+1][k] - 
                        p_dyn.x[0][j+2][k])/(2. * dthe);
            d2vdthe2 = h_d_j * (2. * v.x[0][j][k] - 2. * v.x[0][j+1][k]+v.x[0][j+2][k])/dthe2;
            d2wdthe2 = h_d_j * (2. * w.x[0][j][k] - 2. * w.x[0][j+1][k]+w.x[0][j+2][k])/dthe2;
        }
        if((is_land(h, 0, j, k))&&(is_air(h, 0, j-1, k))&&(is_air(h, 0, j-2, k))){
            dvdthe = h_d_j * (- 3. * v.x[0][j][k] + 4. * v.x[0][j-1][k] 
                - v.x[0][j-2][k])/(2. * dthe);
            dwdthe = h_d_j * (- 3. * w.x[0][j][k] + 4. * w.x[0][j-1][k] 
                 - w.x[0][j-2][k])/(2. * dthe);
            dpdthe = h_d_j * (- 3. * p_dyn.x[0][j][k] 
                + 4. * p_dyn.x[0][j-1][k] - p_dyn.x[0][j-2][k]) 
                   /(2. * dthe);
            d2vdthe2 = h_d_j * (2. * v.x[0][j][k] - 2. * v.x[0][j-1][k] 
                + v.x[0][j-2][k])/dthe2;
            d2wdthe2 = h_d_j * (2. * w.x[0][j][k] - 2. * w.x[0][j-1][k] 
                + w.x[0][j-2][k])/dthe2;
        }
        if(((is_land(h, 0, j, k)) 
            &&((is_air(h, 0, j+1, k))&&(is_land(h, 0, j+2, k)))) 
            ||((j == jm-2)
            &&((is_air(h, 0, j, k))&&(is_land(h, 0, j+1, k))))){
            dvdthe = h_d_j * (v.x[0][j+1][k] - v.x[0][j][k])/dthe;
            dwdthe = h_d_j * (w.x[0][j+1][k] - w.x[0][j][k])/dthe;
            d2vdthe2 = d2wdthe2 = 0.;
        }
        if(((is_land(h, 0, j, k)) 
            &&((is_air(h, 0, j-1, k))&&(is_land(h, 0, j-2, k)))) 
            ||((j == 1)&&((is_land(h, 0, j, k))&&(is_air(h, 0, j-1, k))))){
            dvdthe = h_d_j * (v.x[0][j-1][k] - v.x[0][j][k])/dthe;
            dwdthe = h_d_j * (w.x[0][j-1][k] - w.x[0][j][k])/dthe;
            d2vdthe2 = d2wdthe2 = 0.;
        }
    }
    if((k >= 2)&&(k < km-3)){
            if((is_land(h, 0, j, k))&&(is_air(h, 0, j, k+1)) 
            &&(is_air(h, 0, j, k+2))){
            dvdphi = h_d_k * (- 3. * v.x[0][j][k] + 4. * v.x[0][j][k+1] 
                - v.x[0][j][k+2])/(2. * dphi);
            dwdphi = h_d_k * (- 3. * w.x[0][j][k] + 4. * w.x[0][j][k+1] 
                - w.x[0][j][k+2])/(2. * dphi);
            dpdphi = h_d_k * (- 3. * p_dyn.x[0][j][k] 
                + 4. * p_dyn.x[0][j][k+1] - p_dyn.x[0][j][k+2]) 
               /(2. * dphi);
            d2vdthe2 = h_d_k * (2. * v.x[0][j][k] - 2. * v.x[0][j][k+1] 
                + v.x[0][j][k+2])/dphi2;
            d2wdthe2 = h_d_k * (2. * w.x[0][j][k] - 2. * w.x[0][j][k+1] 
                + w.x[0][j][k+2])/dphi2;
        }
        if((is_land(h, 0, j, k))&&(is_air(h, 0, j, k-1))&&(is_air(h, 0, j, k-2))){
            dvdphi = h_d_k * (- 3. * v.x[0][j][k] + 4. * v.x[0][j][k-1] 
                - v.x[0][j][k-2])/(2. * dphi);
            dwdphi = h_d_k * (- 3. * w.x[0][j][k] + 4. * w.x[0][j][k-1] 
                - w.x[0][j][k-2])/(2. * dphi);
            dpdphi = h_d_k * (- 3. * p_dyn.x[0][j][k] 
                + 4. * p_dyn.x[0][j][k-1] - p_dyn.x[0][j][k-2]) 
               /(2. * dphi);
            d2vdthe2 = h_d_k * (2. * v.x[0][j][k] - 2. * v.x[0][j][k-1] 
                + v.x[0][j][k-2])/dphi2;
            d2wdthe2 = h_d_k * (2. * w.x[0][j][k] - 2. * w.x[0][j][k-1] 
                + w.x[0][j][k-2])/dphi2;
        }
        if(((is_land(h, 0, j, k))
            &&((is_air(h, 0, j, k+1))&&(is_land(h, 0, j, k+2))))||((k == km-2)
            &&((is_air(h, 0, j, k))&&(is_land(h, 0, j, k+1))))){
            dvdphi = h_d_k * (v.x[0][j][k+1] - v.x[0][j][k])/dphi;
            dwdphi = h_d_k * (w.x[0][j][k+1] - w.x[0][j][k])/dphi;
            d2vdphi2 = d2wdphi2 = 0.;    
        }
        if(((is_land(h, 0, j, k)) 
            &&((is_air(h, 0, j, k-1))&&(is_land(h, 0, j, k-2))))||((k == 1)
            &&((is_land(h, 0, j, k))&&(is_air(h, 0, j, k-1))))){
            dvdphi = h_d_k * (v.x[0][j][k-1] - v.x[0][j][k])/dphi;
            dwdphi = h_d_k * (w.x[0][j][k-1] - w.x[0][j][k])/dphi;
            d2vdphi2 = d2wdphi2 = 0.;
        }
    }
    rhs_v.x[0][j][k] = - (v.x[0][j][k] * dvdthe/rm 
        + w.x[0][j][k] * dvdphi/rmsinthe) 
        - coeff_p * dpdthe/rm 
        - (d2vdthe2/rm2 + dvdthe/rm2sinthe * costhe - (1. + costhe * costhe/sinthe2) 
        * h_d_j * v.x[0][j][k] + d2vdphi2/rm2sinthe2 - dwdphi 
        * 2. * costhe/rm2sinthe2)/re 
        + cc * h_0_j * v.x[0][j][k]/dthe2;
    rhs_w.x[0][j][k] = - (v.x[0][j][k] * dwdthe/rm 
        +  w.x[0][j][k] * dwdphi/rmsinthe) 
        - coeff_p * dpdphi/rmsinthe 
        + (d2wdthe2/rm2 + dwdthe/rm2sinthe  * costhe - (1. + costhe * costhe/sinthe2) 
        * h_d_k * w.x[0][j][k] + d2wdphi2/rm2sinthe2 
        + dvdphi * 2. * costhe/rm2sinthe2)/re
        + cc * h_0_k * w.x[0][j][k]/dphi2;
    aux_v.x[0][j][k] = rhs_v.x[0][j][k] + coeff_p * dpdthe/rm;
    aux_w.x[0][j][k] = rhs_w.x[0][j][k] + coeff_p * dpdphi/rmsinthe;

    if(is_land(h, 0, j, k)){
        aux_u.x[0][j][k] = aux_v.x[0][j][k] = aux_w.x[0][j][k] = 0.;
    }

    return;
}


