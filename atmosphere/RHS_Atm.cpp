/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
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

//
#define dxdr_a(X) \
    (h_d_i * ( X->x[ i+1 ][ j ][ k ] - X->x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * exp_rm)
#define dxdr_b(X) \
    (h_d_i * ( - 3. * X->x[ i ][ j ][ k ] + 4. * X->x[ i + 1 ][ j ][ k ] - X->x[ i + 2 ][ j ][ k ] ) / ( 2. * dr ) * exp_rm)

//
#define d2xdr2_a(X) \
    (h_d_i * ( X->x[ i+1 ][ j ][ k ] - 2. * X->x[ i ][ j ][ k ] + X->x[ i-1 ][ j ][ k ] ) / dr2 * exp_2_rm \
    - h_d_i * ( X->x[ i+1 ][ j ][ k ] - X->x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * exp_2_rm_2)
#define d2xdr2_b(X) \
    (h_d_i * ( X->x[ i ][ j ][ k ] - 2. * X->x[ i + 1 ][ j ][ k ] + X->x[ i + 2 ][ j ][ k ] ) / dr2 * exp_2_rm \
    - h_d_i * ( - 3. * X->x[ i ][ j ][ k ] + 4. * X->x[ i + 1 ][ j ][ k ] - X->x[ i + 2 ][ j ][ k ] ) / ( 2. * dr ) * exp_2_rm_2)

//
#define dxdthe_a(X) \
    (h_d_j * ( X->x[ i ][ j+1 ][ k ] - X->x[ i ][ j-1 ][ k ] ) / ( 2. * dthe ))
#define dxdthe_b(X) \
    (h_d_j * ( - 3. * X->x[ i ][ j ][ k ] + 4. * X->x[ i ][ j + 1 ][ k ] - X->x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe ))
#define dxdthe_c(X) \
    (h_d_j * ( - 3. * X->x[ i ][ j ][ k ] + 4. * X->x[ i ][ j - 1 ][ k ] - X->x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe ))
#define dxdthe_d(X) \
    (h_d_j * ( X->x[ i ][ j + 1 ][ k ] - X->x[ i ][ j ][ k ] ) / dthe)
#define dxdthe_e(X) \
    (h_d_j * ( X->x[ i ][ j - 1 ][ k ] - X->x[ i ][ j ][ k ] ) / dthe)

//
#define d2xdthe2_a(X) \
    (h_d_j * ( X->x[ i ][ j+1 ][ k ] - 2. * X->x[ i ][ j ][ k ] + X->x[ i ][ j-1 ][ k ] ) / dthe2)
#define d2xdthe2_b(X) \
    (h_d_j * ( X->x[ i ][ j ][ k ] - 2. * X->x[ i ][ j + 1 ][ k ] + X->x[ i ][ j + 2 ][ k ] ) / dthe2)
#define d2xdthe2_c(X) \
    (h_d_j * ( X->x[ i ][ j ][ k ] - 2. * X->x[ i ][ j - 1 ][ k ] + X->x[ i ][ j - 2 ][ k ] ) / dthe2)

//
#define dxdphi_a(X) \
    (h_d_k * ( X->x[ i ][ j ][ k+1 ] - X->x[ i ][ j ][ k-1 ] ) / ( 2. * dphi ))
#define dxdphi_b(X) \
    (h_d_k * ( - 3. * X->x[ i ][ j ][ k ] + 4. * X->x[ i ][ j ][ k + 1 ] - X->x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi ))
#define dxdphi_c(X) \
    (h_d_k * ( - 3. * X->x[ i ][ j ][ k ] + 4. * X->x[ i ][ j ][ k - 1 ] - X->x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi ))
#define dxdphi_d(X) \
    (h_d_k * ( X->x[ i ][ j ][ k + 1 ] - X->x[ i ][ j ][ k ] ) / dphi)
#define dxdphi_e(X) \
    (h_d_k * ( X->x[ i ][ j ][ k - 1 ] - X->x[ i ][ j ][ k ] ) / dphi)

//
#define d2xdphi2_a(X) \
    (h_d_k * ( X->x[ i ][ j ][ k+1 ] - 2. * X->x[ i ][ j ][ k ] + X->x[ i ][ j ][ k-1 ] ) / dphi2)
#define d2xdphi2_b(X) \
    (h_d_k * ( X->x[ i ][ j ][ k ] - 2. * X->x[ i ][ j ][ k + 1 ] + X->x[ i ][ j ][ k + 2 ] ) / dphi2)
#define d2xdphi2_c(X) \
    (h_d_k * ( X->x[ i ][ j ][ k ] - 2. * X->x[ i ][ j ][ k - 1 ] + X->x[ i ][ j ][ k - 2 ] ) / dphi2)

void cAtmosphereModel::RK_RHS_3D_Atmosphere(int i, int j, int k) 
{
    double k_Force = 1.;// factor for acceleration of convergence processes inside the immersed boundary conditions
    double cc = 1.;
    double dist_coeff = 1.;

    // 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )
    // collection of coefficients9
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;

    // collection of coefficients
    double rm = rad.z[ i ];
    double rm2 = rm * rm;

    // collection of coefficients for coordinate stretching
    double zeta = 3.715;
    double exp_rm = 1. / exp( zeta * rm );
    double exp_2_rm = 1. / exp( 2. * zeta * rm );
    double exp_2_rm_2 = 2. / exp( 2. * zeta * rm );

    // collection of coefficients
    double sinthe = sin( the.z[ j ] );
    double sinthe2 = sinthe * sinthe;
    double costhe = cos( the.z[ j ] );
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;

    //  3D volume iterations in case 1. and 2. order derivatives at walls are needed >>>>>>>>>>>>>>>>>>>>>>>> 
    // only in positive r-direction above ground 
    double topo_step = L_atm / ( double ) ( im-1 );
    double hight = ( double ) i * topo_step;
    double topo_diff = hight - Topography.y[ j ][ k ];
    double h_0_i = topo_diff / topo_step;  // hat distribution function
//    double h_0_i = cc * ( .5 * ( acos ( topo_diff * 3.14 / L_atm ) + 1. ) );   // cosine distribution function, better results for benchmark case

    double h_d_i = 0, h_0_0 = 0;
    if ( ( topo_diff < topo_step ) && ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i-1, j, k ) ) ) ){
        h_0_0 = 1. - h_0_i;
        h_d_i = cc * ( 1. - h_0_0 ); 
    }
    if ( ( topo_diff == topo_step ) || ( is_air ( h, i, j, k ) ) ){
        h_0_i = 1.;
        h_d_i = cc * ( 1. - h_0_i ); 
    }

    // 3D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // only in positive the-direction along northerly boundaries 

    double dist = 0, h_0_j = 0, h_d_j = 0;
    if ( ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j-1, k ) ) ) || 
          ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j+1, k ) ) ) ){
        dist = dist_coeff * dthe;
        h_0_j = dist / dthe;
        h_d_j = cc * ( 1. - h_0_j ); 
    }else{
        h_0_j = 0.;
        h_d_j = cc * ( 1. - h_0_j ); 
    }

    double h_0_k = 0, h_d_k = 0;     
    if ( ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j, k-1 ) ) ) || 
          ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j, k+1 ) ) ) ){
        dist = dist_coeff * dphi;
        h_0_k = dist / dphi;
        h_d_k = cc * ( 1. - h_0_k ); 
    }else{
        h_0_k = 0.; 
        h_d_j = cc * ( 1. - h_0_k ); 
    }

    std::vector<Array*> arrays{&u, &v, &w, &t, &p_dyn, &c, &cloud, &ice, &co2};
    enum array_index {i_u, i_v, i_w, i_t, i_p, i_c, i_cloud, i_ice, i_co, last_array_index};
    std::vector<double> dxdr_vals(last_array_index), 
                        dxdthe_vals(last_array_index), 
                        dxdphi_vals(last_array_index),
                        d2xdr2_vals(last_array_index),
                        d2xdthe2_vals(last_array_index),
                        d2xdphi2_vals(last_array_index);

    // finite differences 1. and 2. order in the free flow field with no contact to solid surfaces
    for(int n=0; n<last_array_index; n++){
        // 1. order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
        dxdr_vals[n] = dxdr_a(arrays[n]);
        dxdthe_vals[n] = dxdthe_a(arrays[n]);
        dxdphi_vals[n] = dxdphi_a(arrays[n]);
        // 2. order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
        d2xdr2_vals[n] = d2xdr2_a(arrays[n]);
        d2xdthe2_vals[n] = d2xdthe2_a(arrays[n]);
        d2xdphi2_vals[n] = d2xdphi2_a(arrays[n]); 
    }

    if ( i < im - 2 ){
        if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i+1, j, k ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdr_vals[n] = dxdr_b(arrays[n]);
                d2xdr2_vals[n] = d2xdr2_b(arrays[n]);
            }
        }
    }

    if ( ( j >= 2 ) && ( j < jm - 3 ) ){
        // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral southward direction
        if ( ( is_land ( h, i, j, k ) ) && ( ( is_air ( h, i, j+1, k ) ) && ( is_air ( h, i, j+2, k ) ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdthe_vals[n] = dxdthe_b(arrays[n]);
                d2xdthe2_vals[n] = d2xdthe2_b(arrays[n]);
            }            
        }
        // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral northward direction
        if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j-1, k ) ) && ( is_air ( h, i, j-2, k ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdthe_vals[n] = dxdthe_c(arrays[n]);
                d2xdthe2_vals[n] = d2xdthe2_c(arrays[n]);
            }
        }

        // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral southward direction
        if ( ( ( is_land ( h, i, j, k ) ) 
                && ( ( is_air ( h, i, j+1, k ) ) && ( is_land ( h, i, j+2, k ) ) ) ) 
                || ( ( j == jm - 2 )
                && ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j+1, k ) ) ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdthe_vals[n] = dxdthe_d(arrays[n]);
                d2xdthe2_vals[n] = 0.;
            }
        }

        // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral northward direction
        if ( ( ( is_land ( h, i, j, k ) ) 
                && ( ( is_air ( h, i, j-1, k ) ) && ( is_land ( h, i, j-2, k ) ) ) ) 
                || ( ( j == 1 )
                && ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j-1, k ) ) ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdthe_vals[n] = dxdthe_e(arrays[n]);
                d2xdthe2_vals[n] = 0.;
            }            
        }
    }


    if ( ( k >= 2 ) && ( k < km - 3 ) ){
        // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in longitudial eastward direction
        if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k+1 ) ) && ( is_air ( h, i, j, k+2 ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdphi_vals[n] = dxdphi_b(arrays[n]);
                d2xdphi2_vals[n] = d2xdphi2_b(arrays[n]);
            }
        }
        
        // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in longitudial westward direction
        if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k-1 ) ) && ( is_air ( h, i, j, k-2 ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdphi_vals[n] = dxdphi_c(arrays[n]);
                d2xdphi2_vals[n] = d2xdphi2_c(arrays[n]);
            }
        }

        // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral eastward direction
        if ( ( ( is_land ( h, i, j, k ) ) 
                && ( ( is_air ( h, i, j, k+1 ) ) && ( is_land ( h, i, j, k+2 ) ) ) ) 
                || ( ( k == km - 2 )
                && ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j, k+1 ) ) ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdphi_vals[n] = dxdphi_e(arrays[n]);
                d2xdphi2_vals[n] = 0.;
            }
        }

        // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral westward direction
        if ( ( ( is_land ( h, i, j, k ) ) 
                && ( ( is_air ( h, i, j, k-1 ) ) && ( is_land ( h, i, j, k-2 ) ) ) ) 
                || ( ( k == 1 )
                && ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k-1 ) ) ) ) ){
            for(int n=0; n<last_array_index; n++){
                dxdphi_vals[n] = dxdphi_e(arrays[n]);
                d2xdphi2_vals[n] = 0.;
            }
        }
    }

    double dudr = dxdr_vals[i_u], dvdr = dxdr_vals[i_v], dwdr = dxdr_vals[i_w], dtdr = dxdr_vals[i_t],
           dpdr = dxdr_vals[i_p], dcdr = dxdr_vals[i_c], dclouddr = dxdr_vals[i_cloud], dicedr = dxdr_vals[i_ice],
           dcodr = dxdr_vals[i_co];
    
    double dudthe = dxdthe_vals[i_u], dvdthe = dxdthe_vals[i_v], dwdthe = dxdthe_vals[i_w], dtdthe = dxdthe_vals[i_t],
           dpdthe = dxdthe_vals[i_p], dcdthe = dxdthe_vals[i_c], dclouddthe = dxdthe_vals[i_cloud], dicedthe = dxdthe_vals[i_ice],
           dcodthe = dxdthe_vals[i_co];
    
    double dudphi = dxdphi_vals[i_u], dvdphi = dxdphi_vals[i_v], dwdphi = dxdphi_vals[i_w], dtdphi = dxdphi_vals[i_t],
           dpdphi = dxdphi_vals[i_p], dcdphi = dxdphi_vals[i_c], dclouddphi = dxdphi_vals[i_cloud], dicedphi = dxdphi_vals[i_ice],
           dcodphi = dxdphi_vals[i_co];

    double d2udr2 = d2xdr2_vals[i_u], d2vdr2 = d2xdr2_vals[i_v], d2wdr2 = d2xdr2_vals[i_w], d2tdr2 = d2xdr2_vals[i_t],
           d2cdr2 = d2xdr2_vals[i_c], d2clouddr2 = d2xdr2_vals[i_cloud], d2icedr2 = d2xdr2_vals[i_ice],
           d2codr2 = d2xdr2_vals[i_co];

    double d2udthe2 = d2xdthe2_vals[i_u], d2vdthe2 = d2xdthe2_vals[i_v], d2wdthe2 = d2xdthe2_vals[i_w], d2tdthe2 = d2xdthe2_vals[i_t],
           d2cdthe2 = d2xdthe2_vals[i_c], d2clouddthe2 = d2xdthe2_vals[i_cloud], d2icedthe2 = d2xdthe2_vals[i_ice],
           d2codthe2 = d2xdthe2_vals[i_co];

    double d2udphi2 = d2xdphi2_vals[i_u], d2vdphi2 = d2xdphi2_vals[i_v], d2wdphi2 = d2xdphi2_vals[i_w], d2tdphi2 = d2xdphi2_vals[i_t],
           d2cdphi2 = d2xdphi2_vals[i_c], d2clouddphi2 = d2xdphi2_vals[i_cloud], d2icedphi2 = d2xdphi2_vals[i_ice],
           d2codphi2 = d2xdphi2_vals[i_co];


    double exp_pressure = g / ( 1.e-2 * gam * R_Air );
    double t_u = t.x[ i ][ j ][ k ] * t_0;  // in K
    double p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );     // given in hPa
    hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
    double p_h = pow ( ( ( t.x[ 0 ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) / 
                         ( t.x[ 0 ][ j ][ k ] * t_0 ) ), exp_pressure ) * p_SL;
    double r_dry = 100. * p_h / ( R_Air * t_u );

    // collection of coefficients
    double coeff_energy = L_atm / ( cp_l * t_0 * u_0 );  // coefficient for the source terms = .00729
    double coeff_buoy = r_air * u_0 * u_0 / L_atm; // coefficient for bouancy term = 0.00482
    double coeff_trans = L_atm / u_0;   // coefficient for the concentration terms = 2000.
    double vapour_evaporation = 0.;
    double vapour_surface = 0.;
    double evap_precip = 0.;
    double coeff_vapour = 1.1574e-5 * L_atm / u_0;// 1.1574e-3 == mm/d to m/s, == .01234

    // Boussineq-approximation for the buoyancy force caused by humid air lighter than dry air
    double r_humid = r_dry * ( 1. + c.x[ i ][ j ][ k ] ) / ( 1. + R_WaterVapour / R_Air * c.x[ i ][ j ][ k ] );

    double RS_buoyancy_Momentum = Buoyancy * ( r_humid - r_dry ) / r_dry * g; // any humid air is less dense than dry air

    BuoyancyForce.x[ i ][ j ][ k ] = - RS_buoyancy_Momentum * coeff_buoy * 1000.;// dimension as pressure in kN/m2

    if ( is_land ( h, i, j, k ) ){
        BuoyancyForce.x[ i ][ j ][ k ] = 0.;
    }

    if ( i == 1 ){
        evap_precip = Evaporation_Dalton.y[ j ][ k ] - Precipitation.y[ j ][ k ];
        if ( evap_precip >= 6. ){
            evap_precip = 6.;     // vapour gradient causes values too high at shelf corners
        }
        if ( evap_precip <= - 6. ){
            evap_precip = - 6.;         // vapour gradient causes values too high at shelf corners
        }
        vapour_surface = r_humid * ( - 3. * c.x[ 0 ][ j ][ k ] + 4. * c.x[ 1 ][ j ][ k ] - c.x[ 3 ][ j ][ k ] ) / 
                         ( 2. * dr ) * ( 1. - 2. * c.x[ 0 ][ j ][ k ] ) * evap_precip;     // 2. ord.
//      vapour_surface = r_humid * ( c.x[ 0 ][ j ][ k ] - c.x[ 1 ][ j ][ k ] ) / dr * ( 1. - 2. * c.x[ 0 ][ j ][ k ] ) * evap_precip;

        vapour_evaporation = + coeff_vapour * vapour_surface;
        if ( is_land ( h, i, j, k ) ){
            vapour_evaporation = 0.;
        }
    }else{
        vapour_evaporation = 0.;
    }

    // Right Hand Side of the time derivative ot temperature, pressure, water vapour concentration and velocity components
    rhs_t.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dtdr + v.x[ i ][ j ][ k ] * dtdthe / rm
            + w.x[ i ][ j ][ k ] * dtdphi / rmsinthe ) + ( d2tdr2 + dtdr * 2. / rm + d2tdthe2 / rm2
            + dtdthe * costhe / rm2sinthe + d2tdphi2 / rm2sinthe2 ) / ( re * pr )
            + coeff_energy * ( S_c.x[ i ][ j ][ k ] + S_r.x[ i ][ j ][ k ] )
            + coeff_energy * ( S_i.x[ i ][ j ][ k ] + S_s.x[ i ][ j ][ k ] );
//            - h_0_i * t.x[ i ][ j ][ k ] * k_Force / dr2;

    rhs_u.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dudr + v.x[ i ][ j ][ k ] * dudthe / rm 
            + w.x[ i ][ j ][ k ] * dudphi / rmsinthe )
            - h_d_i * dpdr / r_air + ( d2udr2 + h_d_i * 2. * u.x[ i ][ j ][ k ] / rm2 + d2udthe2 / rm2
            + 4. * dudr / rm + dudthe * costhe / rm2sinthe + d2udphi2 / rm2sinthe2 ) / re
            - RS_buoyancy_Momentum
            - h_0_i * u.x[ i ][ j ][ k ] * k_Force / dr2;

    rhs_v.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dvdr + v.x[ i ][ j ][ k ] * dvdthe / rm
            + w.x[ i ][ j ][ k ] * dvdphi / rmsinthe ) +
            - h_d_j * dpdthe / rm / r_air + ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
            - ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * h_d_j * v.x[ i ][ j ][ k ] + d2vdphi2 / rm2sinthe2
            + 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
            - h_0_j * v.x[ i ][ j ][ k ] * k_Force / dthe2;

    rhs_w.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dwdr + v.x[ i ][ j ][ k ] * dwdthe / rm
           + w.x[ i ][ j ][ k ] * dwdphi / rmsinthe ) +
            - h_d_k * dpdphi / rmsinthe / r_air + ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2
            + dwdthe / rm2sinthe  * costhe - ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * h_d_k * w.x[ i ][ j ][ k ]
            + d2wdphi2 / rm2sinthe2 + 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
            - h_0_k * w.x[ i ][ j ][ k ] * k_Force / dphi2;

    rhs_c.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dcdr + v.x[ i ][ j ][ k ] * dcdthe / rm
            + w.x[ i ][ j ][ k ] * dcdphi / rmsinthe ) + ( d2cdr2 + dcdr * 2. / rm + d2cdthe2 / rm2
            + dcdthe * costhe / rm2sinthe + d2cdphi2 / rm2sinthe2 ) / ( sc_WaterVapour * re )
            + S_v.x[ i ][ j ][ k ] * coeff_trans
            + vapour_evaporation;
//            - h_0_i * c.x[ i ][ j ][ k ] * k_Force / dr2;

    rhs_cloud.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dclouddr + v.x[ i ][ j ][ k ] * dclouddthe / rm
            + w.x[ i ][ j ][ k ] * dclouddphi / rmsinthe ) + ( d2clouddr2 + dclouddr * 2. / rm + d2clouddthe2 / rm2
            + dclouddthe * costhe / rm2sinthe + d2clouddphi2 / rm2sinthe2 ) / ( sc_WaterVapour * re )
            + S_c.x[ i ][ j ][ k ] * coeff_trans
            - h_0_i * cloud.x[ i ][ j ][ k ] * k_Force / dr2;

    rhs_ice.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dicedr + v.x[ i ][ j ][ k ] * dicedthe / rm
            + w.x[ i ][ j ][ k ] * dicedphi / rmsinthe ) + ( d2icedr2 + dicedr * 2. / rm + d2icedthe2 / rm2
            + dicedthe * costhe / rm2sinthe + d2icedphi2 / rm2sinthe2 ) / ( sc_WaterVapour * re )
            + S_i.x[ i ][ j ][ k ] * coeff_trans
            - h_0_i * ice.x[ i ][ j ][ k ] * k_Force / dr2;

    rhs_co2.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dcodr + v.x[ i ][ j ][ k ] * dcodthe / rm
            + w.x[ i ][ j ][ k ] * dcodphi / rmsinthe ) + ( d2codr2 + dcodr * 2. / rm + d2codthe2 / rm2
            + dcodthe * costhe / rm2sinthe + d2codphi2 / rm2sinthe2 ) / ( sc_CO2 * re )
            - h_0_i * co2.x[ i ][ j ][ k ] * k_Force / dr2;

    // for the Poisson equation to solve for the pressure, pressure gradient substracted from the above RHS
    aux_u.x[ i ][ j ][ k ] = rhs_u.x[ i ][ j ][ k ] + h_d_i * dpdr / r_air;
    aux_v.x[ i ][ j ][ k ] = rhs_v.x[ i ][ j ][ k ] + h_d_j * dpdthe / rm / r_air;
    aux_w.x[ i ][ j ][ k ] = rhs_w.x[ i ][ j ][ k ] + h_d_k * dpdphi / rmsinthe / r_air;

    if ( is_land ( h, i, j, k ) ){
        aux_u.x[ i ][ j ][ k ] = aux_v.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k ] = 0.;
    }
}




void cAtmosphereModel::RK_RHS_2D_Atmosphere(int j, int k){
    //  2D surface iterations
    double k_Force = 1.;// factor for acceleration of convergence processes inside the immersed boundary conditions
    double cc = 1.;
    double dist_coeff = 1.;

    // collection of coefficients
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double rm = rad.z[ 0 ];
    double rm2 = rm * rm;

    // collection of coefficients
    double sinthe = sin( the.z[ j ] );
    double sinthe2 = sinthe * sinthe;
    double costhe = cos( the.z[ j ] );
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;
    double dist = 0, h_0_j = 0, h_d_j = 0;
    // 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>
    // only in positive the-direction along northerly and southerly boundaries 
    if ( ( ( is_air ( h, 0, j, k ) ) && ( is_land ( h, 0, j+1, k ) ) ) || 
          ( ( is_air ( h, 0, j, k ) ) && ( is_land ( h, 0, j-1, k ) ) ) ){
        dist = dist_coeff * dthe;
        h_0_j = dist / dthe;
        h_d_j = cc * ( 1. - h_0_j ); 
    }else{
        h_0_j = 0.; 
        h_d_j = cc * ( 1. - h_0_j ); 
    }

    double h_0_k = 0, h_d_k = 0;
    // 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>
    // only in positive phi-direction on westerly and easterly boundaries 
    if ( ( ( is_air ( h, 0, j, k ) ) && ( is_land ( h, 0, j, k+1 ) ) ) || 
          ( ( is_air ( h, 0, j, k ) ) && ( is_land ( h, 0, j, k-1 ) ) ) ){
        dist = dist_coeff * dphi;
        h_0_k = dist / dphi;
        h_d_k = cc * ( 1. - h_0_k ); 
    }else{
        h_0_k = 0.; 
        h_d_k = cc * ( 1. - h_0_k ); 
    }

    double dvdthe = h_d_j * ( v.x[ 0 ][ j+1 ][ k ] - v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );
    double dwdthe = h_d_j * ( w.x[ 0 ][ j+1 ][ k ] - w.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );
    double dpdthe = h_d_j * ( p_dyn.x[ 0 ][ j+1 ][ k ] - p_dyn.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );

    double dvdphi = h_d_k * ( v.x[ 0 ][ j ][ k+1 ] - v.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
    double dwdphi = h_d_k * ( w.x[ 0 ][ j ][ k+1 ] - w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
    double dpdphi = h_d_k * ( p_dyn.x[ 0 ][ j ][ k+1 ] - p_dyn.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );

    double d2vdthe2 = h_d_j *  ( v.x[ 0 ][ j+1 ][ k ] - 2. * v.x[ 0 ][ j ][ k ] + v.x[ 0 ][ j-1 ][ k ] ) / dthe2;
    double d2wdthe2 = h_d_j * ( w.x[ 0 ][ j+1 ][ k ] - 2. * w.x[ 0 ][ j ][ k ] + w.x[ 0 ][ j-1 ][ k ] ) / dthe2;

    double d2vdphi2 = h_d_k * ( v.x[ 0 ][ j ][ k+1 ] - 2. * v.x[ 0 ][ j ][ k ] + v.x[ 0 ][ j ][ k-1 ] ) / dphi2;
    double d2wdphi2 = h_d_k * ( w.x[ 0 ][ j ][ k+1 ] - 2. * w.x[ 0 ][ j ][ k ] + w.x[ 0 ][ j ][ k-1 ] ) / dphi2;

if ( ( j >= 2 ) && ( j < jm - 3 ) ){
        if ( ( is_land ( h, 0, j, k ) ) 
                && ( ( is_air ( h, 0, j+1, k ) ) 
                && ( is_air ( h, 0, j+2, k ) ) ) ){
            dvdthe = h_d_j * ( - 3. * v.x[ 0 ][ j ][ k ] + 4. * v.x[ 0 ][ j + 1 ][ k ] - v.x[ 0 ][ j + 2 ][ k ] ) 
                        / ( 2. * dthe );
            dwdthe = h_d_j * ( - 3. * w.x[ 0 ][ j ][ k ] + 4. * w.x[ 0 ][ j + 1 ][ k ] - w.x[ 0 ][ j + 2 ][ k ] ) 
                        / ( 2. * dthe );
            dpdthe = h_d_j * ( - 3. * p_dyn.x[ 0 ][ j ][ k ] + 4. * p_dyn.x[ 0 ][ j + 1 ][ k ] - 
                        p_dyn.x[ 0 ][ j + 2 ][ k ] ) / ( 2. * dthe );

            d2vdthe2 = h_d_j * ( 2. * v.x[ 0 ][ j ][ k ] - 2. * v.x[ 0 ][ j + 1 ][ k ] + v.x[ 0 ][ j + 2 ][ k ] ) / dthe2;
            d2wdthe2 = h_d_j * ( 2. * w.x[ 0 ][ j ][ k ] - 2. * w.x[ 0 ][ j + 1 ][ k ] + w.x[ 0 ][ j + 2 ][ k ] ) / dthe2;
        }
        if ( ( is_land ( h, 0, j, k ) ) 
                && ( is_air ( h, 0, j+1, k ) ) ){
            dvdthe = h_d_j * ( v.x[ 0 ][ j + 1 ][ k ] - v.x[ 0 ][ j ][ k ] ) / dthe;
            dwdthe = h_d_j * ( w.x[ 0 ][ j + 1 ][ k ] - w.x[ 0 ][ j ][ k ] ) / dthe;
            dpdthe = h_d_j * ( p_dyn.x[ 0 ][ j + 1 ][ k ] - p_dyn.x[ 0 ][ j ][ k ] ) / dthe;

            d2vdthe2 = d2wdthe2 = 0.;
        }
        if ( ( is_land ( h, 0, j, k ) ) 
            && ( is_air ( h, 0, j-1, k ) ) 
            && ( is_air ( h, 0, j-2, k ) ) ){
            dvdthe = h_d_j * ( - 3. * v.x[ 0 ][ j ][ k ] + 4. * v.x[ 0 ][ j - 1 ][ k ] - v.x[ 0 ][ j - 2 ][ k ] ) 
                        / ( 2. * dthe );
            dwdthe = h_d_j * ( - 3. * w.x[ 0 ][ j ][ k ] + 4. * w.x[ 0 ][ j - 1 ][ k ] - w.x[ 0 ][ j - 2 ][ k ] ) 
                        / ( 2. * dthe );
            dpdthe = h_d_j * ( - 3. * p_dyn.x[ 0 ][ j ][ k ] + 4. * p_dyn.x[ 0 ][ j - 1 ][ k ] - 
                        p_dyn.x[ 0 ][ j - 2 ][ k ] ) / ( 2. * dthe );

            d2vdthe2 = h_d_j * ( 2. * v.x[ 0 ][ j ][ k ] - 2. * v.x[ 0 ][ j - 1 ][ k ] + v.x[ 0 ][ j - 2 ][ k ] ) / dthe2;
            d2wdthe2 = h_d_j * ( 2. * w.x[ 0 ][ j ][ k ] - 2. * w.x[ 0 ][ j - 1 ][ k ] + w.x[ 0 ][ j - 2 ][ k ] ) / dthe2;
        }
        if ( ( is_land ( h, 0, j, k ) ) 
            && ( is_air ( h, 0, j-1, k ) ) ){
            dvdthe = h_d_j * ( v.x[ 0 ][ j ][ k ] - v.x[ 0 ][ j - 1 ][ k ] ) / dthe;
            dwdthe = h_d_j * ( w.x[ 0 ][ j ][ k ] - w.x[ 0 ][ j - 1 ][ k ] ) / dthe;
            dpdthe = h_d_j * ( p_dyn.x[ 0 ][ j ][ k ] - p_dyn.x[ 0 ][ j - 1 ][ k ] ) / dthe;

            d2vdthe2 = d2wdthe2 = 0.;
        }
         d2vdthe2 = d2wdthe2 = 0.;
    }

    if ( ( k >= 2 ) && ( k < km - 3 ) ){
        if ( ( is_land ( h, 0, j, k ) ) 
            && ( is_air ( h, 0, j, k+1 ) ) 
            && ( is_air ( h, 0, j, k+2 ) ) ){
            dvdphi = h_d_k * ( - 3. * v.x[ 0 ][ j ][ k ] + 4. * v.x[ 0 ][ j ][ k + 1 ] - v.x[ 0 ][ j ][ k + 2 ] ) 
                        / ( 2. * dphi );
            dwdphi = h_d_k * ( - 3. * w.x[ 0 ][ j ][ k ] + 4. * w.x[ 0 ][ j ][ k + 1 ] - w.x[ 0 ][ j ][ k + 2 ] ) 
                        / ( 2. * dphi );
            dpdphi = h_d_k * ( - 3. * p_dyn.x[ 0 ][ j ][ k ] + 4. * p_dyn.x[ 0 ][ j ][ k + 1 ] 
                        - p_dyn.x[ 0 ][ j ][ k + 2 ] ) / ( 2. * dphi );

            d2vdthe2 = h_d_k * ( 2. * v.x[ 0 ][ j ][ k ] - 2. * v.x[ 0 ][ j ][ k + 1 ] + v.x[ 0 ][ j ][ k + 2 ] ) / dphi2;
            d2wdthe2 = h_d_k * ( 2. * w.x[ 0 ][ j ][ k ] - 2. * w.x[ 0 ][ j ][ k + 1 ] + w.x[ 0 ][ j ][ k + 2 ] ) / dphi2;
        }
        if ( ( is_land ( h, 0, j, k ) ) 
            && ( is_air ( h, 0, j, k+1 ) ) ){
            dvdphi = h_d_k * ( v.x[ 0 ][ j ][ k + 1 ] - v.x[ 0 ][ j ][ k ] ) / dphi;
            dwdphi = h_d_k * ( w.x[ 0 ][ j ][ k + 1 ] - w.x[ 0 ][ j ][ k ] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[ 0 ][ j ][ k + 1 ] - p_dyn.x[ 0 ][ j ][ k ] ) / dphi;

            d2vdphi2 = d2wdphi2 = 0.;
        }
        if ( ( is_land ( h, 0, j, k ) ) 
            && ( is_air ( h, 0, j, k-1 ) ) 
            && ( is_air ( h, 0, j, k-2 ) ) ){
            dvdphi = h_d_k * ( - 3. * v.x[ 0 ][ j ][ k ] + 4. * v.x[ 0 ][ j ][ k - 1 ] - v.x[ 0 ][ j ][ k - 2 ] ) 
                        / ( 2. * dphi );
            dwdphi = h_d_k * ( - 3. * w.x[ 0 ][ j ][ k ] + 4. * w.x[ 0 ][ j ][ k - 1 ] - w.x[ 0 ][ j ][ k - 2 ] ) 
                        / ( 2. * dphi );
            dpdphi = h_d_k * ( - 3. * p_dyn.x[ 0 ][ j ][ k ] + 4. * p_dyn.x[ 0 ][ j ][ k - 1 ] 
                        - p_dyn.x[ 0 ][ j ][ k - 2 ] ) / ( 2. * dphi );

            d2vdthe2 = h_d_k * ( 2. * v.x[ 0 ][ j ][ k ] - 2. * v.x[ 0 ][ j ][ k - 1 ] + v.x[ 0 ][ j ][ k - 2 ] ) / dphi2;
            d2wdthe2 = h_d_k * ( 2. * w.x[ 0 ][ j ][ k ] - 2. * w.x[ 0 ][ j ][ k - 1 ] + w.x[ 0 ][ j ][ k - 2 ] ) / dphi2;
        }
        if ( ( is_land ( h, 0, j, k ) ) 
            && ( is_air ( h, 0, j, k-1 ) ) ){
            dvdphi = h_d_k * ( v.x[ 0 ][ j ][ k ] - v.x[ 0 ][ j ][ k - 1 ] ) / dphi;
            dwdphi = h_d_k * ( w.x[ 0 ][ j ][ k ] - w.x[ 0 ][ j ][ k - 1 ] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[ 0 ][ j ][ k ] - p_dyn.x[ 0 ][ j ][ k - 1 ] ) / dphi;

            d2vdphi2 = d2wdphi2 = 0.;
        }
        d2vdphi2 = d2wdphi2 = 0.;
    }else{
        if ( ( is_land ( h, 0, j, k ) ) && ( is_air ( h, 0, j, k+1 ) ) ){
            dvdphi = h_d_k * ( v.x[ 0 ][ j ][ k + 1 ] - v.x[ 0 ][ j ][ k ] ) / dphi;
            dwdphi = h_d_k * ( w.x[ 0 ][ j ][ k + 1 ] - w.x[ 0 ][ j ][ k ] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[ 0 ][ j ][ k + 1 ] - p_dyn.x[ 0 ][ j ][ k ] ) / dphi;
        }
        if ( ( is_air ( h, 0, j, k ) ) && ( is_land ( h, 0, j, k-1 ) ) ){
            dvdphi = h_d_k * ( v.x[ 0 ][ j ][ k ] - v.x[ 0 ][ j ][ k - 1 ] ) / dphi;
            dwdphi = h_d_k * ( w.x[ 0 ][ j ][ k ] - w.x[ 0 ][ j ][ k - 1 ] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[ 0 ][ j ][ k ] - p_dyn.x[ 0 ][ j ][ k - 1 ] ) / dphi;
        }
        d2vdthe2 = d2wdthe2 = 0.;
        d2vdphi2 = d2wdphi2 = 0.;
    }

    rhs_v.x[ 0 ][ j ][ k ] = - ( v.x[ 0 ][ j ][ k ] * dvdthe / rm + w.x[ 0 ][ j ][ k ] * dvdphi / rmsinthe ) +
                - h_d_j * dpdthe / rm / r_air - ( d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
                - ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[ 0 ][ j ][ k ]
                + d2vdphi2 / rm2sinthe2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
                - h_0_j * v.x[ 0 ][ j ][ k ] * k_Force / dthe2;

    rhs_w.x[ 0 ][ j ][ k ] = - ( v.x[ 0 ][ j ][ k ] * dwdthe / rm +  w.x[ 0 ][ j ][ k ] * dwdphi / rmsinthe ) +
                - h_d_k * dpdphi / rmsinthe / r_air + ( d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
                - ( 1. + costhe * costhe / sinthe2 ) * h_d_k * w.x[ 0 ][ j ][ k ]
                + d2wdphi2 / rm2sinthe2 + dvdphi * 2. * costhe / rm2sinthe2 ) / re
                - h_0_k * w.x[ 0 ][ j ][ k ] * k_Force / dphi2;

    aux_v.x[ 0 ][ j ][ k ] = rhs_v.x[ 0 ][ j ][ k ] + h_d_j * dpdthe / rm / r_air;
    aux_w.x[ 0 ][ j ][ k ] = rhs_w.x[ 0 ][ j ][ k ] + h_d_k * dpdphi / rmsinthe / r_air;
}


