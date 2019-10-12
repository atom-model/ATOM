/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
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


void cHydrosphereModel::RK_RHS_3D_Hydrosphere(int i, int j, int k)
{
    // collection of coefficients for phase transformation

    double k_Force = 1.;// factor for accelleration of convergence processes inside the immersed boundary conditions
    double cc = 1.;
    double dist_coeff = 1.;
    double h_d_i = 0;

    // 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )

    // collection of coefficients
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;

    // collection of coefficients
    double rm = rad.z[i];
    double rm2 = rm * rm;

    // collection of coefficients
    double sinthe = sin( the.z[j] );
    double sinthe2 = sinthe * sinthe;
    double costhe = cos( the.z[j] );
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;

    //  3D volume iterations in case 1. and 2. order derivatives at walls are needed >>>>>>>>>>>>>>>>>>>>>>>> 
    // only in positive r-direction above ground 

    double topo_step = L_hyd / ( double ) ( im-1 );
    double height = ( double ) i * topo_step;
    double topo_diff = fabs ( height - Bathymetry.y[j][k] );
    double h_0_i = topo_diff / topo_step;  // hat distribution function
    // cosine distribution function, better results for benchmark case
    //    double h_0_i = cc * ( .5 * ( acos ( topo_diff * 3.14 / L_hyd ) + 1. ) ); 

    if ( ( topo_diff <= topo_step ) && ( ( is_water ( h, i, j, k ) ) 
        && ( is_land ( h, i-1, j, k ) ) ) ){
        h_d_i = cc * ( 1. - h_0_i ); 
    }

    if ( ( is_water ( h, i, j, k ) ) && ( is_water ( h, i-1, j, k ) ) ){
        h_0_i = 0.;
        h_d_i = 1.;
    }

    if ( ( is_land ( h, i, j, k ) ) && ( is_land ( h, i-1, j, k ) ) ){
        h_0_i = 1.;
        h_d_i = 0.;
    }

    // 3D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // only in positive the-direction along northerly and southerly boundaries 
    double dist = 0, h_0_j = 0, h_d_j = 0;
    if ( ( ( is_water ( h, i, j, k ) ) && ( is_land ( h, i, j-1, k ) ) ) || 
          ( ( is_water ( h, i, j, k ) ) && ( is_land ( h, i, j+1, k ) ) ) ){
        dist = dist_coeff * dthe;
        h_0_j = dist / dthe;
        h_d_j = cc * ( 1. - h_0_j ); 
    }else{
        h_0_j = 0.;
        h_d_j = cc * ( 1. - h_0_j ); 
    }

    double h_0_k = 0, h_d_k = 0;     
    if ( ( ( is_water ( h, i, j, k ) ) && ( is_land ( h, i, j, k-1 ) ) ) || 
          ( ( is_water ( h, i, j, k ) ) && ( is_land ( h, i, j, k+1 ) ) ) ){
        dist = dist_coeff * dphi;
        h_0_k = dist / dphi;
        h_d_k = cc * ( 1. - h_0_k ); 
    }else{
        h_0_k = 0.; 
        h_d_k = cc * ( 1. - h_0_k ); 
    }


    // buoyancy effects by salt water density changes
    double drodc = .7;    // gradient given in kg/m³/m
    double salt_water_ref = r_water.x[i][j][k] + drodc * c.x[i][j][k];
                          // common linear approach for salt water based on fresh water

    double coeff_buoy = L_hyd / ( u_0 * u_0 );
                          // coefficient for the buoyancy term == 16.0

    // Boussineq-approximation for the buoyancy force caused by salinity, higher salinity causes negative buoyancy
    double RS_buoyancy_Momentum = Buoyancy * coeff_buoy * g * ( r_salt_water.x[i][j][k] - salt_water_ref )
                                  / salt_water_ref;  // buoyancy based on water density, salt water is heavier than fresh water 
    //    double RS_buoyancy_Momentum = 0.;  // test case

    BuoyancyForce_3D.x[i][j][k] = RS_buoyancy_Momentum / coeff_buoy;
                                        // dimension as pressure in kN/m2
    Salt_Balance.x[i][j][k] = salt_water_ref - r_salt_water.x[i][j][k];
                                        // difference of salinity compared to average

    if ( Salt_Balance.x[i][j][k] < 0. ){
        Salt_Diffusion.x[i][j][k] = Salt_Balance.x[i][j][k];
                                        // for negativ salinity balance, higher than reference
        Salt_Finger.x[i][j][k] = 0.;
    }else{
        Salt_Finger.x[i][j][k] = Salt_Balance.x[i][j][k];
                                        // for positiv salinity balance, lower than reference
        Salt_Diffusion.x[i][j][k] = 0.;
    }
    if ( is_land( h, i, j, k) ){
        Salt_Balance.x[i][j][k] = 0.;
        Salt_Finger.x[i][j][k] = 0.;
        Salt_Diffusion.x[i][j][k] = 0.;
        BuoyancyForce_3D.x[i][j][k] = 0.;
        r_water.x[i][j][k] = 1007.;
        r_salt_water.x[i][j][k] = 1032.;
    }

    // 2. order derivative for temperature, pressure, salt concentrations and velocity components

    // computation of initial and boundary conditions for the v and w velocity component
    // computation at the surface

    //  3D volume iterations

    // 1st order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
    double dudr = h_d_i * ( u.x[i+1][j][k] - u.x[i-1][j][k] ) / ( 2. * dr );
    double dvdr = h_d_i * ( v.x[i+1][j][k] - v.x[i-1][j][k] ) / ( 2. * dr );
    double dwdr = h_d_i * ( w.x[i+1][j][k] - w.x[i-1][j][k] ) / ( 2. * dr );
    double dtdr = h_d_i * ( t.x[i+1][j][k] - t.x[i-1][j][k] ) / ( 2. * dr );
    double dpdr = h_d_i * ( p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k] ) / ( 2. * dr );
    double dcdr = h_d_i * ( c.x[i+1][j][k] - c.x[i-1][j][k] ) / ( 2. * dr );

    double dudthe = h_d_j * ( u.x[i][j+1][k] - u.x[i][j-1][k] ) / ( 2. * dthe );
    double dvdthe = h_d_j * ( v.x[i][j+1][k] - v.x[i][j-1][k] ) / ( 2. * dthe );
    double dwdthe = h_d_j * ( w.x[i][j+1][k] - w.x[i][j-1][k] ) / ( 2. * dthe );
    double dtdthe = h_d_j * ( t.x[i][j+1][k] - t.x[i][j-1][k] ) / ( 2. * dthe );
    double dpdthe = h_d_j * ( p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k] ) / ( 2. * dthe );
    double dcdthe = h_d_j * ( c.x[i][j+1][k] - c.x[i][j-1][k] ) / ( 2. * dthe );

    double dudphi = h_d_k * ( u.x[i][j][k+1] - u.x[i][j][k-1] ) / ( 2. * dphi );
    double dvdphi = h_d_k * ( v.x[i][j][k+1] - v.x[i][j][k-1] ) / ( 2. * dphi );
    double dwdphi = h_d_k * ( w.x[i][j][k+1] - w.x[i][j][k-1] ) / ( 2. * dphi );
    double dtdphi = h_d_k * ( t.x[i][j][k+1] - t.x[i][j][k-1] ) / ( 2. * dphi );
    double dpdphi = h_d_k * ( p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1] ) / ( 2. * dphi );
    double dcdphi = h_d_k * ( c.x[i][j][k+1] - c.x[i][j][k-1] ) / ( 2. * dphi );

    // 2nd order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
    double d2udr2 = h_d_i * ( u.x[i+1][j][k] - 2. * u.x[i][j][k] + u.x[i-1][j][k] ) / dr2;
    double d2vdr2 = h_d_i * ( v.x[i+1][j][k] - 2. * v.x[i][j][k] + v.x[i-1][j][k] ) / dr2; 
    double d2wdr2 = h_d_i * ( w.x[i+1][j][k] - 2. * w.x[i][j][k] + w.x[i-1][j][k] ) / dr2; 
    double d2tdr2 = h_d_i * ( t.x[i+1][j][k] - 2. * t.x[i][j][k] + t.x[i-1][j][k] ) / dr2; 
    double d2cdr2 = h_d_i * ( c.x[i+1][j][k] - 2. * c.x[i][j][k] + c.x[i-1][j][k] ) / dr2; 

    double d2udthe2 = h_d_j * ( u.x[i][j+1][k] - 2. * u.x[i][j][k] + u.x[i][j-1][k] ) / dthe2;
    double d2vdthe2 = h_d_j * ( v.x[i][j+1][k] - 2. * v.x[i][j][k] + v.x[i][j-1][k] ) / dthe2;
    double d2wdthe2 = h_d_j * ( w.x[i][j+1][k] - 2. * w.x[i][j][k] + w.x[i][j-1][k] ) / dthe2;
    double d2tdthe2 = h_d_j * ( t.x[i][j+1][k] - 2. * t.x[i][j][k] + t.x[i][j-1][k] ) / dthe2;
    double d2cdthe2 = h_d_j * ( c.x[i][j+1][k] - 2. * c.x[i][j][k] + c.x[i][j-1][k] ) / dthe2;

    double d2udphi2 = h_d_k * ( u.x[i][j][k+1] - 2. * u.x[i][j][k] + u.x[i][j][k-1] ) / dphi2;
    double d2vdphi2 = h_d_k * ( v.x[i][j][k+1] - 2. * v.x[i][j][k] + v.x[i][j][k-1] ) / dphi2;
    double d2wdphi2 = h_d_k * ( w.x[i][j][k+1] - 2. * w.x[i][j][k] + w.x[i][j][k-1] ) / dphi2;
    double d2tdphi2 = h_d_k * ( t.x[i][j][k+1] - 2. * t.x[i][j][k] + t.x[i][j][k-1] ) / dphi2;
    double d2cdphi2 = h_d_k * ( c.x[i][j][k+1] - 2. * c.x[i][j][k] + c.x[i][j][k-1] ) / dphi2;


    // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in radial direction
    if ( i < im - 2 ){
        if ( ( is_land ( h, i, j, k ) ) && ( ( is_water ( h, i+1, j, k ) ) && ( is_water ( h, i+2, j, k ) ) ) ){
            dudr = h_d_i * ( - 3. * u.x[i][j][k] + 4. * u.x[i + 1][j][k] - u.x[i + 2][j][k] ) / ( 2. * dr );
            dvdr = h_d_i * ( - 3. * v.x[i][j][k] + 4. * v.x[i + 1][j][k] - v.x[i + 2][j][k] ) / ( 2. * dr );
            dwdr = h_d_i * ( - 3. * w.x[i][j][k] + 4. * w.x[i + 1][j][k] - w.x[i + 2][j][k] ) / ( 2. * dr );
            dtdr = h_d_i * ( - 3. * t.x[i][j][k] + 4. * t.x[i + 1][j][k] - t.x[i + 2][j][k] ) / ( 2. * dr );
            dpdr = h_d_i * ( - 3. * p_dyn.x[i][j][k] + 4. * p_dyn.x[i + 1][j][k] - p_dyn.x[i + 2][j][k] ) / ( 2. * dr );
            dcdr = h_d_i * ( - 3. * u.x[i][j][k] + 4. * u.x[i + 1][j][k] - u.x[i + 2][j][k] ) / ( 2. * dr );

            d2udr2 = h_d_i * ( u.x[i][j][k] - 2. * u.x[i + 1][j][k] + u.x[i + 2][j][k] ) / dr2; 
            d2vdr2 = h_d_i * ( v.x[i][j][k] - 2. * v.x[i + 1][j][k] + v.x[i + 2][j][k] ) / dr2; 
            d2wdr2 = h_d_i * ( w.x[i][j][k] - 2. * w.x[i + 1][j][k] + w.x[i + 2][j][k] ) / dr2; 
            d2tdr2 = h_d_i * ( t.x[i][j][k] - 2. * t.x[i + 1][j][k] + t.x[i + 2][j][k] ) / dr2; 
            d2cdr2 = h_d_i * ( u.x[i][j][k] - 2. * u.x[i + 1][j][k] + u.x[i + 2][j][k] ) / dr2; 
        }
    }

    // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral southward direction
    if ( ( j >= 2 ) && ( j <= jm - 3 ) ){
        if ( ( is_land ( h, i, j, k ) ) && ( ( is_water ( h, i, j+1, k ) ) && ( is_water ( h, i, j+2, k ) ) ) ){
            dudthe = h_d_j * ( - 3. * u.x[i][j][k] + 4. * u.x[i][j + 1][k] - u.x[i][j + 2][k] ) / ( 2. * dthe );
            dvdthe = h_d_j * ( - 3. * v.x[i][j][k] + 4. * v.x[i][j + 1][k] - v.x[i][j + 2][k] ) / ( 2. * dthe );
            dwdthe = h_d_j * ( - 3. * w.x[i][j][k] + 4. * w.x[i][j + 1][k] - w.x[i][j + 2][k] ) / ( 2. * dthe );
            dtdthe = h_d_j * ( - 3. * t.x[i][j][k] + 4. * t.x[i][j + 1][k] - t.x[i][j + 2][k] ) / ( 2. * dthe );
            dpdthe = h_d_j * ( - 3. * p_dyn.x[i][j][k] + 4. * p_dyn.x[i][j + 1][k] - p_dyn.x[i][j + 2][k] ) / ( 2. * dthe );
            dcdthe = h_d_j * ( - 3. * c.x[i][j][k] + 4. * c.x[i][j + 1][k] - c.x[i][j + 2][k] ) / ( 2. * dthe );

            d2udthe2 = h_d_j * ( u.x[i][j][k] - 2. * u.x[i][j + 1][k] + u.x[i][j + 2][k] ) / dthe2;
            d2vdthe2 = h_d_j * ( v.x[i][j][k] - 2. * v.x[i][j + 1][k] + v.x[i][j + 2][k] ) / dthe2;
            d2wdthe2 = h_d_j * ( w.x[i][j][k] - 2. * w.x[i][j + 1][k] + w.x[i][j + 2][k] ) / dthe2;
            d2tdthe2 = h_d_j * ( t.x[i][j][k] - 2. * t.x[i][j + 1][k] + t.x[i][j + 2][k] ) / dthe2;
            d2cdthe2 = h_d_j * ( c.x[i][j][k] - 2. * c.x[i][j + 1][k] + c.x[i][j + 2][k] ) / dthe2;
        }

        // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral northward direction
        if ( ( is_land ( h, i, j, k ) ) && ( is_water ( h, i, j-1, k ) ) && ( is_water ( h, i, j-2, k ) ) ){
            dudthe = h_d_j * ( - 3. * u.x[i][j][k] + 4. * u.x[i][j - 1][k] - u.x[i][j - 2][k] ) / ( 2. * dthe );
            dvdthe = h_d_j * ( - 3. * v.x[i][j][k] + 4. * v.x[i][j - 1][k] - v.x[i][j - 2][k] ) / ( 2. * dthe );
            dwdthe = h_d_j * ( - 3. * w.x[i][j][k] + 4. * w.x[i][j - 1][k] - w.x[i][j - 2][k] ) / ( 2. * dthe );
            dtdthe = h_d_j * ( - 3. * t.x[i][j][k] + 4. * t.x[i][j - 1][k] - t.x[i][j - 2][k] ) / ( 2. * dthe );
            dpdthe = h_d_j * ( - 3. * p_dyn.x[i][j][k] + 4. * p_dyn.x[i][j - 1][k] - p_dyn.x[i][j - 2][k] ) / ( 2. * dthe );
            dcdthe = h_d_j * ( - 3. * c.x[i][j][k] + 4. * c.x[i][j - 1][k] - c.x[i][j - 2][k] ) / ( 2. * dthe );

            d2udthe2 = h_d_j * ( u.x[i][j][k] - 2. * u.x[i][j - 1][k] + u.x[i][j - 2][k] ) / dthe2;
            d2vdthe2 = h_d_j * ( v.x[i][j][k] - 2. * v.x[i][j - 1][k] + v.x[i][j - 2][k] ) / dthe2;
            d2wdthe2 = h_d_j * ( w.x[i][j][k] - 2. * w.x[i][j - 1][k] + w.x[i][j - 2][k] ) / dthe2;
            d2tdthe2 = h_d_j * ( t.x[i][j][k] - 2. * t.x[i][j - 1][k] + t.x[i][j - 2][k] ) / dthe2;
            d2cdthe2 = h_d_j * ( c.x[i][j][k] - 2. * c.x[i][j - 1][k] + c.x[i][j - 2][k] ) / dthe2;
        }

        // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral southward direction
        if ( ( ( is_land ( h, i, j, k ) ) 
            && ( ( is_water ( h, i, j+1, k ) ) && ( is_land ( h, i, j+2, k ) ) ) ) 
            || ( ( j == jm - 2 )
            && ( ( is_water ( h, i, j, k ) ) && ( is_land ( h, i, j+1, k ) ) ) ) ){
            dudthe = h_d_j * ( u.x[i][j + 1][k] - u.x[i][j][k] ) / dthe;
            dvdthe = h_d_j * ( v.x[i][j + 1][k] - v.x[i][j][k] ) / dthe;
            dwdthe = h_d_j * ( w.x[i][j + 1][k] - w.x[i][j][k] ) / dthe;
            dtdthe = h_d_j * ( t.x[i][j + 1][k] - t.x[i][j][k] ) / dthe;
            dpdthe = h_d_j * ( p_dyn.x[i][j + 1][k] - p_dyn.x[i][j][k] ) / dthe;
            dcdthe = h_d_j * ( c.x[i][j + 1][k] - c.x[i][j][k] ) / dthe;

            d2udthe2 = d2vdthe2 = d2wdthe2 = d2tdthe2 = d2cdthe2 = 0.;
        }

        // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral northward direction
        if ( ( ( is_land ( h, i, j, k ) ) 
            && ( ( is_water ( h, i, j-1, k ) ) && ( is_land ( h, i, j-2, k ) ) ) ) 
            || ( ( j == 1 )
            && ( ( is_land ( h, i, j, k ) ) && ( is_water ( h, i, j-1, k ) ) ) ) ){
            dudthe = h_d_j * ( u.x[i][j - 1][k] - u.x[i][j][k] ) / dthe;
            dvdthe = h_d_j * ( v.x[i][j - 1][k] - v.x[i][j][k] ) / dthe;
            dwdthe = h_d_j * ( w.x[i][j - 1][k] - w.x[i][j][k] ) / dthe;
            dtdthe = h_d_j * ( t.x[i][j - 1][k] - t.x[i][j][k] ) / dthe;
            dpdthe = h_d_j * ( p_dyn.x[i][j - 1][k] - p_dyn.x[i][j][k] ) / dthe;
            dcdthe = h_d_j * ( c.x[i][j - 1][k] - c.x[i][j][k] ) / dthe;

            d2udthe2 = d2vdthe2 = d2wdthe2 = d2tdthe2 = d2cdthe2 = 0.;
        }
    }


    // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in longitudial eastward direction
    if ( ( k >= 2 ) && ( k <= km - 3 ) ){
        if ( ( is_land ( h, i, j, k ) ) && ( is_water ( h, i, j, k+1 ) ) && ( is_water ( h, i, j, k+2 ) ) ){
            dudphi = h_d_k * ( - 3. * u.x[i][j][k] + 4. * u.x[i][j][k + 1] - u.x[i][j][k + 2] ) / ( 2. * dphi );
            dvdphi = h_d_k * ( - 3. * v.x[i][j][k] + 4. * v.x[i][j][k + 1] - v.x[i][j][k + 2] ) / ( 2. * dphi );
            dwdphi = h_d_k * ( - 3. * w.x[i][j][k] + 4. * w.x[i][j][k + 1] - w.x[i][j][k + 2] ) / ( 2. * dphi );
            dtdphi = h_d_k * ( - 3. * t.x[i][j][k] + 4. * t.x[i][j][k + 1] - t.x[i][j][k + 2] ) / ( 2. * dphi );
            dpdphi = h_d_k * ( - 3. * p_dyn.x[i][j][k] + 4. * p_dyn.x[i][j][k + 1] - p_dyn.x[i][j][k + 2] ) / ( 2. * dphi );
            dcdphi = h_d_k * ( - 3. * c.x[i][j][k] + 4. * c.x[i][j][k + 1] - c.x[i][j][k + 2] ) / ( 2. * dphi );

            d2udphi2 = h_d_k * ( u.x[i][j][k] - 2. * u.x[i][j][k + 1] + u.x[i][j][k + 2] ) / dphi2;
            d2vdphi2 = h_d_k * ( v.x[i][j][k] - 2. * v.x[i][j][k + 1] + v.x[i][j][k + 2] ) / dphi2;
            d2wdphi2 = h_d_k * ( w.x[i][j][k] - 2. * w.x[i][j][k + 1] + w.x[i][j][k + 2] ) / dphi2;
            d2tdphi2 = h_d_k * ( t.x[i][j][k] - 2. * t.x[i][j][k + 1] + t.x[i][j][k + 2] ) / dphi2;
            d2cdphi2 = h_d_k * ( c.x[i][j][k] - 2. * c.x[i][j][k + 1] + c.x[i][j][k + 2] ) / dphi2;
        }

        // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in longitudial westward direction
        if ( ( is_land ( h, i, j, k ) ) && ( is_water ( h, i, j, k-1 ) ) && ( is_water ( h, i, j, k-2 ) ) ){
            dudphi = h_d_k * ( - 3. * u.x[i][j][k] + 4. * u.x[i][j][k - 1] - u.x[i][j][k - 2] ) / ( 2. * dphi );
            dvdphi = h_d_k * ( - 3. * v.x[i][j][k] + 4. * v.x[i][j][k - 1] - v.x[i][j][k - 2] ) / ( 2. * dphi );
            dwdphi = h_d_k * ( - 3. * w.x[i][j][k] + 4. * w.x[i][j][k - 1] - w.x[i][j][k - 2] ) / ( 2. * dphi );
            dtdphi = h_d_k * ( - 3. * t.x[i][j][k] + 4. * t.x[i][j][k - 1] - t.x[i][j][k - 2] ) / ( 2. * dphi );
            dpdphi = h_d_k * ( - 3. * p_dyn.x[i][j][k] + 4. * p_dyn.x[i][j][k - 1] - p_dyn.x[i][j][k - 2] ) / ( 2. * dphi );
            dcdphi = h_d_k * ( - 3. * c.x[i][j][k] + 4. * c.x[i][j][k - 1] - c.x[i][j][k - 2] ) / ( 2. * dphi );

            d2udphi2 = h_d_k * ( u.x[i][j][k] - 2. * u.x[i][j][k - 1] + u.x[i][j][k - 2] ) / dphi2;
            d2vdphi2 = h_d_k * ( v.x[i][j][k] - 2. * v.x[i][j][k - 1] + v.x[i][j][k - 2] ) / dphi2;
            d2wdphi2 = h_d_k * ( w.x[i][j][k] - 2. * w.x[i][j][k - 1] + w.x[i][j][k - 2] ) / dphi2;
            d2tdphi2 = h_d_k * ( t.x[i][j][k] - 2. * t.x[i][j][k - 1] + t.x[i][j][k - 2] ) / dphi2;
            d2cdphi2 = h_d_k * ( c.x[i][j][k] - 2. * c.x[i][j][k - 1] + c.x[i][j][k - 2] ) / dphi2;
        }

        // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral eastward direction
        if ( ( ( is_land ( h, i, j, k ) ) 
            && ( ( is_water ( h, i, j, k+1 ) ) && ( is_land ( h, i, j, k+2 ) ) ) ) 
            || ( ( k == km - 2 )
            && ( ( is_water ( h, i, j, k ) ) && ( is_land ( h, i, j, k+1 ) ) ) ) ){
            dudphi = h_d_k * ( u.x[i][j][k + 1] - u.x[i][j][k] ) / dphi;
            dvdphi = h_d_k * ( v.x[i][j][k + 1] - v.x[i][j][k] ) / dphi;
            dwdphi = h_d_k * ( w.x[i][j][k + 1] - w.x[i][j][k] ) / dphi;
            dtdphi = h_d_k * ( t.x[i][j][k + 1] - t.x[i][j][k] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[i][j][k + 1] - p_dyn.x[i][j][k] ) / dphi;
            dcdphi = h_d_k * ( c.x[i][j][k + 1] - c.x[i][j][k] ) / dphi;

            d2udphi2 = d2vdphi2 = d2wdphi2 = d2tdphi2 = d2cdphi2 = 0.;
        }

        // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral westward direction
        if ( ( ( is_land ( h, i, j, k ) ) 
            && ( ( is_water ( h, i, j, k-1 ) ) && ( is_land ( h, i, j, k-2 ) ) ) ) 
            || ( ( k == 1 )
            && ( ( is_land ( h, i, j, k ) ) && ( is_water ( h, i, j, k-1 ) ) ) ) ){
            dudphi = h_d_k * ( u.x[i][j][k - 1] - u.x[i][j][k] ) / dphi;
            dvdphi = h_d_k * ( v.x[i][j][k - 1] - v.x[i][j][k] ) / dphi;
            dwdphi = h_d_k * ( w.x[i][j][k - 1] - w.x[i][j][k] ) / dphi;
            dtdphi = h_d_k * ( t.x[i][j][k - 1] - t.x[i][j][k] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[i][j][k - 1] - p_dyn.x[i][j][k] ) / dphi;
            dcdphi = h_d_k * ( c.x[i][j][k - 1] - c.x[i][j][k] ) / dphi;

            d2udphi2 = d2vdphi2 = d2wdphi2 = d2tdphi2 = d2cdphi2 = 0.;
        }
    }


    // Right Hand Side of the time derivative ot temperature, pressure, salt concentration and velocity components
    //  3D volume iterations

    rhs_t.x[i][j][k] = - ( u.x[i][j][k] * dtdr + v.x[i][j][k] * dtdthe / rm 
            + w.x[i][j][k] * dtdphi / rmsinthe ) + ( d2tdr2 + dtdr * 2. / rm + d2tdthe2 / rm2 
            + dtdthe * costhe / rm2sinthe + d2tdphi2 / rm2sinthe2 ) / ( re * pr );
            //- h_0_i * t.x[i][j][k] * k_Force / dr2;

    rhs_u.x[i][j][k] = - ( u.x[i][j][k] * dudr + v.x[i][j][k] * dudthe / rm 
            + w.x[i][j][k] * dudphi / rmsinthe ) 
            + dpdr / salt_water_ref + ( d2udr2 + h_d_i * 2. * u.x[i][j][k] / rm2 + d2udthe2 / rm2 
            + 4. * dudr / rm + dudthe * costhe / rm2sinthe + d2udphi2 / rm2sinthe2 ) / re 
            + RS_buoyancy_Momentum 
            - h_0_i * u.x[i][j][k] * k_Force / dthe2;

    rhs_v.x[i][j][k] = - ( u.x[i][j][k] * dvdr + v.x[i][j][k] * dvdthe / rm 
            + w.x[i][j][k] * dvdphi / rmsinthe ) 
            - dpdthe / rm / salt_water_ref + ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe 
            - ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[i][j][k] + d2vdphi2 / rm2sinthe2 
            + 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re 
            - h_0_j * v.x[i][j][k] * k_Force / dthe2;

    rhs_w.x[i][j][k] = - ( u.x[i][j][k] * dwdr + v.x[i][j][k] * dwdthe / rm 
            + w.x[i][j][k] * dwdphi / rmsinthe ) 
            - dpdphi / rmsinthe / salt_water_ref + ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 
            + dwdthe / rm2sinthe  * costhe - ( 1. + costhe * costhe / sinthe2 ) * h_d_k * w.x[i][j][k] 
            + d2wdphi2 / rm2sinthe2 + 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re 
            - h_0_k * w.x[i][j][k] * k_Force / dphi2;

    rhs_c.x[i][j][k] = - ( u.x[i][j][k] * dcdr + v.x[i][j][k] * dcdthe / rm 
            + w.x[i][j][k] * dcdphi / rmsinthe ) + ( d2cdr2 + dcdr * 2. / rm + d2cdthe2 / rm2 
            + dcdthe * costhe / rm2sinthe + d2cdphi2 / rm2sinthe2 ) / ( sc * re ); 
            //- h_0_i * c.x[i][j][k] * k_Force / dr2;

    // for the Poisson equation to solve for the pressure, pressure gradient sbstracted from the RHS
    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + h_d_i * dpdr / salt_water_ref;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + h_d_j * dpdthe / rm / salt_water_ref;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + h_d_k * dpdphi / rmsinthe / salt_water_ref;

    if ( is_land( h, i, j, k) ){
        aux_u.x[i][j][k] = aux_v.x[i][j][k] = aux_w.x[i][j][k] = 0.;
    }
}

void cHydrosphereModel::RK_RHS_2D_Hydrosphere( int j, int k)
{

    //  2D surface iterations
    double k_Force = 1.;// factor for accelleration of convergence processes inside the immersed boundary conditions
    double cc = 1.;
    double dist_coeff = 1.;

    // collection of coefficients
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;

    double rm = rad.z[im-1];
    double rm2 = rm * rm;

    // collection of coefficients
    double sinthe = sin( the.z[j] );
    double sinthe2 = sinthe * sinthe;
    double costhe = cos( the.z[j] );
    double rmsinthe = rm * sinthe;
    double rm2sinthe = rm2 * sinthe;
    double rm2sinthe2 = rm2 * sinthe2;

    double dist = 0, h_0_j = 0, h_d_j = 0;

    // 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>
    // only in positive the-direction along northerly and southerly boundaries 
    if ( ( ( is_land ( h, im-1, j, k ) ) && ( is_land ( h, im-1, j+1, k ) ) ) || 
          ( ( is_land ( h, im-1, j, k ) ) && ( is_land ( h, im-1, j-1, k ) ) ) ){
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
    if ( ( ( is_land ( h, im-1, j, k ) ) && ( is_land ( h, im-1, j, k+1 ) ) ) || 
          ( ( is_land ( h, im-1, j, k ) ) && ( is_land ( h, im-1, j, k-1 ) ) ) ){
        dist = dist_coeff * dphi;
        h_0_k = dist / dphi;
        h_d_k = cc * ( 1. - h_0_k ); 
    }else{
        h_0_k = 0.; 
        h_d_k = cc * ( 1. - h_0_k ); 
    }

    double dvdthe = h_d_j * ( v.x[im-1][j+1][k] - v.x[im-1][j-1][k] ) / ( 2. * dthe );
    double dwdthe = h_d_j * ( w.x[im-1][j+1][k] - w.x[im-1][j-1][k] ) / ( 2. * dthe );
    double dpdthe = h_d_j * ( p_dyn.x[im-1][j+1][k] - p_dyn.x[im-1][j-1][k] ) / ( 2. * dthe );

    double dvdphi = h_d_k * ( v.x[im-1][j][k+1] - v.x[im-1][j][k-1] ) / ( 2. * dphi );
    double dwdphi = h_d_k * ( w.x[im-1][j][k+1] - w.x[im-1][j][k-1] ) / ( 2. * dphi );
    double dpdphi = h_d_k * ( p_dyn.x[im-1][j][k+1] - p_dyn.x[im-1][j][k-1] ) / ( 2. * dphi );

    double d2vdthe2 = h_d_j * ( v.x[im-1][j+1][k] - 2. * v.x[im-1][j][k] + v.x[im-1][j-1][k] ) / dthe2;
    double d2wdthe2 = h_d_j * ( w.x[im-1][j+1][k] - 2. * w.x[im-1][j][k] + w.x[im-1][j-1][k] ) / dthe2;

    double d2vdphi2 = h_d_k * ( v.x[im-1][j][k+1] - 2. * v.x[im-1][j][k] + v.x[im-1][j][k-1] ) / dphi2;
    double d2wdphi2 = h_d_k * ( w.x[im-1][j][k+1] - 2. * w.x[im-1][j][k] + w.x[im-1][j][k-1] ) / dphi2;

    if ( ( j >= 2 ) && ( j < jm - 3 ) ){
        if ( ( is_land( h, im-1, j, k) ) 
                && ( ( is_land( h, im-1, j+1, k ) ) 
                && ( is_land( h, im-1, j+2, k ) ) ) ){
            dvdthe = h_d_j * ( - 3. * v.x[im-1][j][k] + 4. * v.x[im-1][j + 1][k] - v.x[im-1][j + 2][k] ) 
                        / ( 2. * dthe );
            dwdthe = h_d_j * ( - 3. * w.x[im-1][j][k] + 4. * w.x[im-1][j + 1][k] - w.x[im-1][j + 2][k] ) 
                        / ( 2. * dthe );
            dpdthe = h_d_j * ( - 3. * p_dyn.x[im-1][j][k] + 4. * p_dyn.x[im-1][j + 1][k] - 
                        p_dyn.x[im-1][j + 2][k] ) / ( 2. * dthe );

            d2vdthe2 = h_d_j * ( 2. * v.x[im-1][j][k] - 2. * v.x[im-1][j + 1][k] + v.x[im-1][j + 2][k] ) / dthe2;
            d2wdthe2 = h_d_j * ( 2. * w.x[im-1][j][k] - 2. * w.x[im-1][j + 1][k] + w.x[im-1][j + 2][k] ) / dthe2;
        }
        if ( ( is_land( h, im-1, j, k) ) 
                && ( is_land( h, im-1, j+1, k ) ) ){
            dvdthe = h_d_j * ( v.x[im-1][j + 1][k] - v.x[im-1][j][k] ) / dthe;
            dwdthe = h_d_j * ( w.x[im-1][j + 1][k] - w.x[im-1][j][k] ) / dthe;
            dpdthe = h_d_j * ( p_dyn.x[im-1][j + 1][k] - p_dyn.x[im-1][j][k] ) / dthe;

            d2vdthe2 = d2wdthe2 = 0.;
        }
        if ( ( is_land( h, im-1, j, k) ) 
            && ( is_land( h, im-1, j-1, k ) ) 
            && ( is_land( h, im-1, j-2, k ) ) ){
            dvdthe = h_d_j * ( - 3. * v.x[im-1][j][k] + 4. * v.x[im-1][j - 1][k] - v.x[im-1][j - 2][k] ) 
                        / ( 2. * dthe );
            dwdthe = h_d_j * ( - 3. * w.x[im-1][j][k] + 4. * w.x[im-1][j - 1][k] - w.x[im-1][j - 2][k] ) 
                        / ( 2. * dthe );
            dpdthe = h_d_j * ( - 3. * p_dyn.x[im-1][j][k] + 4. * p_dyn.x[im-1][j - 1][k] - 
                        p_dyn.x[im-1][j - 2][k] ) / ( 2. * dthe );

            d2vdthe2 = h_d_j * ( 2. * v.x[im-1][j][k] - 2. * v.x[im-1][j - 1][k] + v.x[im-1][j - 2][k] ) / dthe2;
            d2wdthe2 = h_d_j * ( 2. * w.x[im-1][j][k] - 2. * w.x[im-1][j - 1][k] + w.x[im-1][j - 2][k] ) / dthe2;
        }
        if ( ( is_land( h, im-1, j, k) ) 
            && ( is_land( h, im-1, j-1, k ) ) ){
            dvdthe = h_d_j * ( v.x[im-1][j][k] - v.x[im-1][j - 1][k] ) / dthe;
            dwdthe = h_d_j * ( w.x[im-1][j][k] - w.x[im-1][j - 1][k] ) / dthe;
            dpdthe = h_d_j * ( p_dyn.x[im-1][j][k] - p_dyn.x[im-1][j - 1][k] ) / dthe;

            d2vdthe2 = d2wdthe2 = 0.;
        }
        d2vdthe2 = d2wdthe2 = 0.;
    }


    if ( ( k >= 2 ) && ( k < km - 3 ) ){
        if ( ( is_land( h, im-1, j, k) ) 
            && ( is_land( h, im-1, j, k+1) ) 
            && ( is_land( h, im-1, j, k+2 ) ) ){
            dvdphi = h_d_k * ( - 3. * v.x[im-1][j][k] + 4. * v.x[im-1][j][k + 1] - v.x[im-1][j][k + 2] ) 
                        / ( 2. * dphi );
            dwdphi = h_d_k * ( - 3. * w.x[im-1][j][k] + 4. * w.x[im-1][j][k + 1] - w.x[im-1][j][k + 2] ) 
                        / ( 2. * dphi );
            dpdphi = h_d_k * ( - 3. * p_dyn.x[im-1][j][k] + 4. * p_dyn.x[im-1][j][k + 1] 
                        - p_dyn.x[im-1][j][k + 2] ) / ( 2. * dphi );

            d2vdthe2 = h_d_k * ( 2. * v.x[im-1][j][k] - 2. * v.x[im-1][j][k + 1] + v.x[im-1][j][k + 2] ) / dphi2;
            d2wdthe2 = h_d_k * ( 2. * w.x[im-1][j][k] - 2. * w.x[im-1][j][k + 1] + w.x[im-1][j][k + 2] ) / dphi2;
        }
        if ( ( is_land( h, im-1, j, k) ) 
            && ( is_land( h, im-1, j, k+1) ) ){
            dvdphi = h_d_k * ( v.x[im-1][j][k + 1] - v.x[im-1][j][k] ) / dphi;
            dwdphi = h_d_k * ( w.x[im-1][j][k + 1] - w.x[im-1][j][k] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[im-1][j][k + 1] - p_dyn.x[im-1][j][k] ) / dphi;

            d2vdphi2 = d2wdphi2 = 0.;
        }
        if ( ( is_land( h, im-1, j, k) ) 
            && ( is_land( h, im-1, j, k-1) ) 
            && ( is_land( h, im-1, j, k-2) ) ){
            dvdphi = h_d_k * ( - 3. * v.x[im-1][j][k] + 4. * v.x[im-1][j][k - 1] - v.x[im-1][j][k - 2] ) 
                        / ( 2. * dphi );
            dwdphi = h_d_k * ( - 3. * w.x[im-1][j][k] + 4. * w.x[im-1][j][k - 1] - w.x[im-1][j][k - 2] ) 
                        / ( 2. * dphi );
            dpdphi = h_d_k * ( - 3. * p_dyn.x[im-1][j][k] + 4. * p_dyn.x[im-1][j][k - 1] 
                        - p_dyn.x[im-1][j][k - 2] ) / ( 2. * dphi );

            d2vdthe2 = h_d_k * ( 2. * v.x[im-1][j][k] - 2. * v.x[im-1][j][k - 1] + v.x[im-1][j][k - 2] ) / dphi2;
            d2wdthe2 = h_d_k * ( 2. * w.x[im-1][j][k] - 2. * w.x[im-1][j][k - 1] + w.x[im-1][j][k - 2] ) / dphi2;
        }
        if ( ( is_land( h, im-1, j, k) ) 
            && ( is_land( h, im-1, j, k-1) ) ){
            dvdphi = h_d_k * ( v.x[im-1][j][k] - v.x[im-1][j][k - 1] ) / dphi;
            dwdphi = h_d_k * ( w.x[im-1][j][k] - w.x[im-1][j][k - 1] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[im-1][j][k] - p_dyn.x[im-1][j][k - 1] ) / dphi;

            d2vdphi2 = d2wdphi2 = 0.;
        }
        d2vdphi2 = d2wdphi2 = 0.;
    }else{
        if ( ( is_land( h, im-1, j, k) ) && ( is_land( h, im-1, j, k+1) ) ){
            dvdphi = h_d_k * ( v.x[im-1][j][k + 1] - v.x[im-1][j][k] ) / dphi;
            dwdphi = h_d_k * ( w.x[im-1][j][k + 1] - w.x[im-1][j][k] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[im-1][j][k + 1] - p_dyn.x[im-1][j][k] ) / dphi;
        }
        if ( ( h.x[im-1][j][k] == 0. ) && ( is_land( h, im-1, j, k-1) ) ){
            dvdphi = h_d_k * ( v.x[im-1][j][k] - v.x[im-1][j][k - 1] ) / dphi;
            dwdphi = h_d_k * ( w.x[im-1][j][k] - w.x[im-1][j][k - 1] ) / dphi;
            dpdphi = h_d_k * ( p_dyn.x[im-1][j][k] - p_dyn.x[im-1][j][k - 1] ) / dphi;
        }
        d2vdthe2 = d2wdthe2 = 0.;
        d2vdphi2 = d2wdphi2 = 0.;
    }

    rhs_v.x[im-1][j][k] = - ( v.x[im-1][j][k] * dvdthe / rm + w.x[im-1][j][k] * dvdphi / rmsinthe ) +
                - h_d_j * dpdthe / rm / r_0_water - ( d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
                - ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[im-1][j][k] + d2vdphi2 / rm2sinthe2 
                - dwdphi * 2. * costhe / rm2sinthe2 ) / re
                - h_0_j * v.x[im-1][j][k] * k_Force / dthe2;

    rhs_w.x[im-1][j][k] = - ( v.x[im-1][j][k] * dwdthe / rm +  w.x[im-1][j][k] * dwdphi / rmsinthe ) +
                - h_d_k * dpdphi / rmsinthe / r_0_water + ( d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
                - ( 1. + costhe * costhe / sinthe2 ) * h_d_k * w.x[im-1][j][k] + d2wdphi2 / rm2sinthe2 
                + dvdphi * 2. * costhe / rm2sinthe2 ) / re
                - h_0_k * w.x[im-1][j][k] * k_Force / dphi2;

    aux_v.x[im-1][j][k] = rhs_v.x[im-1][j][k] + h_d_j * dpdthe / rm / r_0_water;
    aux_w.x[im-1][j][k] = rhs_w.x[im-1][j][k] + h_d_k * dpdphi / rmsinthe / r_0_water;
}





