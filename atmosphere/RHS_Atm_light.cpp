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

#define dxdr_a(X, dx) \
    ((X->x[i+1][j][k] - X->x[i-1][j][k])/(2.0 * dx))
#define dxdthe_a(X, dx) \
    ((X->x[i][j+1][k] - X->x[i][j-1][k])/(2.0 * dx))
#define dxdphi_a(X, dx) \
    ((X->x[i][j][k+1] - X->x[i][j][k-1])/(2.0 * dx))

void cAtmosphereModel::RK_RHS_3D_Atmosphere(int i, int j, int k){
    double rm = rad.z[i];
    double exp_rm = 1.0/(rm + 1.0);
    double sinthe = sin(the.z[j]);
    if(sinthe == 0.0) sinthe = 1.0e-5;
    double costhe = cos(the.z[j]);
    if(j > 90) costhe = - costhe;
    double rmsinthe = rm * sinthe;

    std::vector<Array*> arrays_1{&u, &v, &w, &t, &p_dyn, &c, &cloud, 
        &ice, &gr, &co2};

    enum array_index_1{i_u_1, i_v_1, i_w_1, i_t_1, i_p_1, i_c_1, 
        i_cloud_1, i_ice_1, i_g_1, i_co_1, 
        last_array_index_1};

    std::vector<double> dxdr_vals(last_array_index_1), 
                        dxdthe_vals(last_array_index_1), 
                        dxdphi_vals(last_array_index_1);

// field gradients
    for(int n=0; n<last_array_index_1; n++){
        dxdr_vals[n] = dxdr_a(arrays_1[n], dr); // 2. order accurate
        dxdthe_vals[n] = dxdthe_a(arrays_1[n], dthe); // 2. order accurate
        dxdphi_vals[n] = dxdphi_a(arrays_1[n], dphi); // 2. order accurate
    }

    double dudr = dxdr_vals[i_u_1], dvdr = dxdr_vals[i_v_1], 
           dwdr = dxdr_vals[i_w_1], dtdr = dxdr_vals[i_t_1],
           dpdr = dxdr_vals[i_p_1], dcdr = dxdr_vals[i_c_1], 
           dclouddr = dxdr_vals[i_cloud_1], dicedr = dxdr_vals[i_ice_1],
           dgdr = dxdr_vals[i_g_1], dcodr = dxdr_vals[i_co_1];

    double dudthe = dxdthe_vals[i_u_1], dvdthe = dxdthe_vals[i_v_1], 
           dwdthe = dxdthe_vals[i_w_1], dtdthe = dxdthe_vals[i_t_1],
           dpdthe = dxdthe_vals[i_p_1], dcdthe = dxdthe_vals[i_c_1], 
           dclouddthe = dxdthe_vals[i_cloud_1], dicedthe = dxdthe_vals[i_ice_1],
           dgdthe = dxdthe_vals[i_g_1], dcodthe = dxdthe_vals[i_co_1];

    double dudphi = dxdphi_vals[i_u_1], dvdphi = dxdphi_vals[i_v_1], 
           dwdphi = dxdphi_vals[i_w_1], dtdphi = dxdphi_vals[i_t_1],
           dpdphi = dxdphi_vals[i_p_1], dcdphi = dxdphi_vals[i_c_1], 
           dclouddphi = dxdphi_vals[i_cloud_1], dicedphi = dxdphi_vals[i_ice_1],
           dgdphi = dxdphi_vals[i_g_1], dcodphi = dxdphi_vals[i_co_1];

    double buoyancy = 1.0;

    double coeff_buoy = r_air * g; // coefficient allows the buoyancy term
    double coeff_energy_p = u_0 * u_0/(cp_l * t_0); // coefficient for the source terms = 2.33e-4 (Eckert-number)
    double coeff_energy = 1.0/(c_0 * cp_l * t_0); // coefficient for the source terms = 1.041e-4
    double coeff_u_p = 0.5 * r_air * u_0 * u_0/L_atm; // coefficient for pressure force term = 0.002408
    double coeff_trans = 1.0; // coefficient for the water vapour term = 1.
    double coeff_MC_vel = L_atm/(u_0 * u_0); // coefficient for MC_v and MC_w term = 250.0
    double coeff_MC_q = L_atm/(u_0 * c_0); // coefficient for MC_q term = 57142.9
    double coeff_MC_t = L_atm/(u_0 * t_0); // coefficient for MC_t term = 7.322

// influence of the Coriolis force
    double coriolis = 1.0;
    double coeff_Coriolis = r_air * u_0 * omega; // coefficient for Coriolis term = 7.024e-4
    double coriolis_rad = 2.0 * costhe * w.x[i][j][k];
    double coriolis_the = - 2.0 * sinthe * w.x[i][j][k];
    double coriolis_phi = 2.0 * (sinthe * v.x[i][j][k] 
        - costhe * u.x[i][j][k]);

// acting forces
    CoriolisForce.x[i][j][k] = coriolis * coeff_Coriolis 
        * sqrt((pow(coriolis_rad, 2) 
        + pow(coriolis_the, 2) 
        + pow(coriolis_phi, 2))/3.0);
    BuoyancyForce.x[i][j][k] = buoyancy * coeff_buoy 
        * (t.x[i][j][k] - 1.0);
    PressureGradientForce.x[i][j][k] = - coeff_u_p 
        * sqrt((pow((dpdr * exp_rm), 2) 
        + pow(dpdthe/rm, 2) 
        + pow(dpdphi/rmsinthe, 2))/3.0);

//    double t_u = t.x[i][j][k] * t_0; // in K
//    double E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
//    double E_Ice = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
//    double q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
//    double q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
    double coeff_L = r_air * c_0 * lv * u_0/L_atm; // coefficient for Q_Latent = 53.0964
/*
    double Q_Latent_Ice = 0.0;
//    double coeff_S = lamda * t_0/L_atm; // coefficient for Q_Sensible = 0.0004473
    if(c.x[i][j][k] >= 0.85 * q_Rain){
        Q_Latent.x[i][j][k] = (u.x[i][j][k] * dclouddr 
        + v.x[i][j][k]/rm * dclouddthe 
        + w.x[i][j][k]/rmsinthe * dclouddphi);
    }else Q_Latent.x[i][j][k] = 0.0;
    if(c.x[i][j][k] >= 0.85 * q_Ice){
        Q_Latent_Ice = (u.x[i][j][k] * dicedr 
        + v.x[i][j][k]/rm * dicedthe 
        + w.x[i][j][k]/rmsinthe * dicedphi);
    }else Q_Latent_Ice = 0.0;

    Q_Latent.x[i][j][k] = coeff_L 
        * (Q_Latent.x[i][j][k] + Q_Latent_Ice);
    Q_Latent.x[i][j][k] = 0.0; 
*/
/*
    Q_Sensible.x[i][j][k] = coeff_S 
        * (d2tdr2 + dtdr * 2./rm + d2tdthe2/rm2
        + dtdthe * costhe/rm2sinthe + d2tdphi2/rm2sinthe2);  // sensible heat in [W/m2] from energy transport equation
*/
    if(is_land(h, i, j, k)){
//        BuoyancyForce.x[i][j][k] = 0.0;
//        PressureGradientForce.x[i][j][k] = 0.0;
        CoriolisForce.x[i][j][k] = 0.0;
//        Q_Latent.x[i][j][k] = 0.0;
//        Q_Sensible.x[i][j][k] = 0.0;
    }

// transport terms in the Euler-equations
    double pressure_t = 0.0;
    double transport_t = 0.0;
    double transport_u = 0.0;
    double transport_v = 0.0;
    double transport_w = 0.0;
    double transport_c = 0.0;
    double transport_cloud = 0.0;
    double transport_ice = 0.0;
    double transport_g = 0.0;
    double transport_co2 = 0.0;

    pressure_t = coeff_energy_p * (u.x[i][j][k] * dpdr * exp_rm
        + v.x[i][j][k] * dpdthe/rm 
        + w.x[i][j][k] * dpdphi/rmsinthe);

    transport_t = u.x[i][j][k] * dtdr * exp_rm + v.x[i][j][k] * dtdthe/rm
        + w.x[i][j][k] * dtdphi/rmsinthe; 
    transport_u = u.x[i][j][k] * dudr * exp_rm + v.x[i][j][k] * dudthe/rm 
        + w.x[i][j][k] * dudphi/rmsinthe;
    transport_v = u.x[i][j][k] * dvdr * exp_rm + v.x[i][j][k] * dvdthe/rm
        + w.x[i][j][k] * dvdphi/rmsinthe;
    transport_w = u.x[i][j][k] * dwdr * exp_rm + v.x[i][j][k] * dwdthe/rm
        + w.x[i][j][k] * dwdphi/rmsinthe;
    transport_c = u.x[i][j][k] * dcdr * exp_rm + v.x[i][j][k] * dcdthe/rm
        + w.x[i][j][k] * dcdphi/rmsinthe;
    transport_cloud = u.x[i][j][k] * dclouddr * exp_rm + v.x[i][j][k] * dclouddthe/rm
        + w.x[i][j][k] * dclouddphi/rmsinthe;
    transport_ice = u.x[i][j][k] * dicedr * exp_rm + v.x[i][j][k] * dicedthe/rm
        + w.x[i][j][k] * dicedphi/rmsinthe;
    transport_g = u.x[i][j][k] * dgdr * exp_rm + v.x[i][j][k] * dgdthe/rm
        + w.x[i][j][k] * dgdphi/rmsinthe;
    transport_co2 = u.x[i][j][k] * dcodr * exp_rm + v.x[i][j][k] * dcodthe/rm
        + w.x[i][j][k] * dcodphi/rmsinthe;

// right hand sides for Runge-Kutta solution scheme
    rhs_t.x[i][j][k] = 
        + pressure_t
        - transport_t 
        + coeff_MC_t * MC_t.x[i][j][k]
        + coeff_energy * (S_c.x[i][j][k] + S_r.x[i][j][k]) 
            * lv * r_humid.x[i][j][k]
        + coeff_energy * (S_i.x[i][j][k] + S_s.x[i][j][k] + S_g.x[i][j][k]) 
            * ls * r_humid.x[i][j][k]
        + Q_Latent.x[i][j][k]/coeff_L;

    rhs_u.x[i][j][k] = 
        - dpdr * exp_rm
        - transport_u 
        + buoyancy * (t.x[i][j][k] - 1.0)
        + coriolis * coriolis_rad;

    rhs_v.x[i][j][k] = 
        - dpdthe/rm
        - transport_v
        + coriolis * coriolis_the
        + coeff_MC_vel * MC_v.x[i][j][k];

    rhs_w.x[i][j][k] = 
        - dpdphi/rmsinthe
        - transport_w
        + coriolis * coriolis_phi
        + coeff_MC_vel * MC_w.x[i][j][k];

    rhs_c.x[i][j][k] = 
        - transport_c
        + coeff_trans * S_v.x[i][j][k] * r_humid.x[i][j][k]
        + coeff_MC_q * MC_q.x[i][j][k];

    rhs_cloud.x[i][j][k] = 
        - transport_cloud
        + coeff_trans * S_c.x[i][j][k] * r_humid.x[i][j][k];

    rhs_ice.x[i][j][k] = 
        - transport_ice
        + coeff_trans * S_i.x[i][j][k] * r_humid.x[i][j][k];

    rhs_g.x[i][j][k] = 
        - transport_g
        + coeff_trans * S_g.x[i][j][k] * r_humid.x[i][j][k];

    rhs_co2.x[i][j][k] = 
        - transport_co2;

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr * exp_rm;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe/rm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi/rmsinthe;
/*
    cout.precision(8);
    if((i == 1)&&(j == 70)&&(k == 180)) cout << " RHS_atm turbulent code" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   dt = " << dt 
        << "   dr = " << dr 
        << "   dthe = " << dthe 
        << "   dphi = " << dphi << endl

        << "   rm = " << rm 
        << "   rm2 = " << rm2 
        << "   rmsinthe = " << rmsinthe 
        << "   costhe = " << costhe << endl 

        << "   topo_step = " << topo_step
        << "   height = " << height
        << "   Topography = " << Topography.y[j][k]
        << "   topo_diff = " << topo_diff << endl << endl

        << "   h_0_i = " << h_0_i
        << "   h_d_i = " << h_d_i << endl << endl

        << "   t = " << t.x[i][j][k]
        << "   u = " << u.x[i][j][k]
        << "   v = " << v.x[i][j][k]
        << "   w = " << w.x[i][j][k]
        << "   p_dyn = " << p_dyn.x[i][j][k] << endl

        << "   tn = " << tn.x[i][j][k]
        << "   un = " << un.x[i][j][k]
        << "   vn = " << vn.x[i][j][k]
        << "   wn = " << wn.x[i][j][k]

        << "   dpdr = " << dpdr
        << "   dtdr = " << dtdr
        << "   dudr = " << dudr
        << "   dvdr = " << dvdr
        << "   dwdr = " << dwdr

        << "   dpdthe = " << dpdthe
        << "   dtdthe = " << dtdthe
        << "   dudthe = " << dudthe
        << "   dvdthe = " << dvdthe
        << "   dwdthe = " << dwdthe

        << "   dpdphi = " << dpdphi
        << "   dtdphi = " << dtdphi
        << "   dudphi = " << dudphi
        << "   dvdphi = " << dvdphi
        << "   dwdphi = " << dwdphi
        
        << "   rhs_u = " << rhs_u.x[i][j][k]
        << "   rhs_v = " << rhs_v.x[i][j][k]
        << "   rhs_w = " << rhs_w.x[i][j][k]
        << "   rhs_t = " << rhs_t.x[i][j][k] << endl

        << "   aux_u = " << aux_u.x[i][j][k]
        << "   aux_v = " << aux_v.x[i][j][k]
        << "   aux_w = " << aux_w.x[i][j][k] << endl << endl

        << "   transport_u = " << transport_u
        << "   transport_v = " << transport_v
        << "   transport_w = " << transport_w;
*/
    return;
}
/*
*
*/


