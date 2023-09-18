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

#define dxdr_a(X, dx) \
    ((X->x[i+1][j][k] - X->x[i-1][j][k])/(2.0 * dx))
#define dxdthe_a(X, dx) \
    ((X->x[i][j+1][k] - X->x[i][j-1][k])/(2.0 * dx))
#define dxdphi_a(X, dx) \
    ((X->x[i][j][k+1] - X->x[i][j][k-1])/(2.0 * dx))

void cHydrosphereModel::RK_RHS_3D_Hydrosphere(int i, int j, int k){
    double rm = rad.z[i];
    double sinthe = sin(the.z[j]);
    if(sinthe == 0.0) sinthe = 1.0e-5;
    double costhe = cos(the.z[j]);
    if(j > 90) costhe = - costhe;
    double rmsinthe = rm * sinthe;

    std::vector<Array*> arrays_1{&u, &v, &w, &t, &p_dyn, &c};

    enum array_index_1{i_u_1, i_v_1, i_w_1, i_t_1, i_p_1, i_c_1, 
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
           dpdr = dxdr_vals[i_p_1], dcdr = dxdr_vals[i_c_1];

    double dudthe = dxdthe_vals[i_u_1], dvdthe = dxdthe_vals[i_v_1], 
           dwdthe = dxdthe_vals[i_w_1], dtdthe = dxdthe_vals[i_t_1],
           dpdthe = dxdthe_vals[i_p_1], dcdthe = dxdthe_vals[i_c_1];

    double dudphi = dxdphi_vals[i_u_1], dvdphi = dxdphi_vals[i_v_1], 
           dwdphi = dxdphi_vals[i_w_1], dtdphi = dxdphi_vals[i_t_1],
           dpdphi = dxdphi_vals[i_p_1], dcdphi = dxdphi_vals[i_c_1];

    double buoyancy = 0.0;

    double coeff_buoy = r_0_water * g; // coefficient allows the buoyancy term
    double coeff_energy_p = u_0 * u_0/(cp_w * t_0); // coefficient for the source term = 2.33e-4
    double coeff_u_p = 0.5 * r_0_water * u_0 * u_0/L_hyd; // coefficient for pressure force term = 0.1489

// Boussineq-approximation for the buoyancy force caused by salinity, higher salinity causes negative buoyancy
// buoyancy effects by salt water density changes
    double drodc = 0.7;    // gradient given in (kg/m³)/m
    double salt_water_ref = r_water.x[i][j][k] 
        + drodc * L_hyd/(double)(im-1) * (double)((im-1) - i);  // linear approximation for the first 200m depth
    Salt_Balance.x[i][j][k] = salt_water_ref - r_salt_water.x[i][j][k]; // difference of salinity compared to average
    if(Salt_Balance.x[i][j][k] < 0.0){
        Salt_Diffusion.x[i][j][k] = Salt_Balance.x[i][j][k]; // for negativ salinity balance, higher than salt_water_ref
        Salt_Finger.x[i][j][k] = 0.0;
    }else{
        Salt_Finger.x[i][j][k] = Salt_Balance.x[i][j][k]; // for positiv salinity balance, lower than salt_water_ref
        Salt_Diffusion.x[i][j][k] = 0.0;
    }

// influence of the Coriolis force
    double coriolis = 1.0;
    double coeff_Coriolis = r_0_water * u_0 * omega; // coefficient for Coriolis term = 0.01797
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
        * sqrt((pow(dpdr, 2) 
        + pow(dpdthe/rm, 2) 
        + pow(dpdphi/rmsinthe, 2))/3.0);
    if(is_land(h, i, j, k)){
        BuoyancyForce.x[i][j][k] = 0.0;
        PressureGradientForce.x[i][j][k] = 0.0;
        CoriolisForce.x[i][j][k] = 0.0;
    }

// transport terms in the Euler-equations
    double pressure_t = 0.0;
    double transport_t = 0.0;
    double transport_u = 0.0;
    double transport_v = 0.0;
    double transport_w = 0.0;
    double transport_c = 0.0;

    pressure_t = coeff_energy_p * (u.x[i][j][k] * dpdr
        + v.x[i][j][k] * dpdthe/rm 
        + w.x[i][j][k] * dpdphi/rmsinthe);

    transport_t = u.x[i][j][k] * dtdr + v.x[i][j][k] * dtdthe/rm
        + w.x[i][j][k] * dtdphi/rmsinthe; 
    transport_u = u.x[i][j][k] * dudr + v.x[i][j][k] * dudthe/rm 
        + w.x[i][j][k] * dudphi/rmsinthe;
    transport_v = u.x[i][j][k] * dvdr + v.x[i][j][k] * dvdthe/rm
        + w.x[i][j][k] * dvdphi/rmsinthe;
    transport_w = u.x[i][j][k] * dwdr + v.x[i][j][k] * dwdthe/rm
        + w.x[i][j][k] * dwdphi/rmsinthe;
    transport_c = u.x[i][j][k] * dcdr + v.x[i][j][k] * dcdthe/rm
        + w.x[i][j][k] * dcdphi/rmsinthe;

// right hand sides for Runge-Kutta solution scheme
    rhs_t.x[i][j][k] = 
        + pressure_t
        - transport_t; 

    rhs_u.x[i][j][k] = 
        - dpdr
        - transport_u 
        + buoyancy * (t.x[i][j][k] - 1.0)
        + coriolis * coriolis_rad;

    rhs_v.x[i][j][k] = 
        - dpdthe/rm
        - transport_v
        + coriolis * coriolis_the;

    rhs_w.x[i][j][k] = 
        - dpdphi/rmsinthe
        - transport_w
        + coriolis * coriolis_phi;

    rhs_c.x[i][j][k] = 
        - transport_c;

    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe/rm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi/rmsinthe;
/*
    cout.precision(8);
    cout << "Ocean code" << endl
        << "northern hemisphere" << endl
        << "   i = " << i << "   j = " << j << "   k = " << k  << endl
        << "   dt = " << dt 
        << "   dr = " << dr 
        << "   dthe = " << dthe 
        << "   dphi = " << dphi
        << "   rm > rad = " << rm 
        << "   sinthe = " << sinthe
        << "   costhe = " << costhe
        << "   rmsinthe = " << rmsinthe
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

