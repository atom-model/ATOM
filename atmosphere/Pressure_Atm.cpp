/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to compute the pressure independent of the other variables
*/

#include <iostream>
#include <cmath>

#include "Array.h"
#include "Array_2D.h"
#include "MinMax_Atm.h"
#include "Utils.h"
#include "cAtmosphereModel.h"

using namespace std;
using namespace AtomUtils;

void cAtmosphereModel::computePressure_3D (){
    const  int c43 = 4./3., c13 = 1./3.;
    // boundary conditions for the r-direction, loop index i

    logger() << "enter $$$$$$$$$$$$$$$$$ Pressure_Atm::computePressure_3D: p_dyn: " << p_dyn.max() * u_0 * u_0 * r_air *.01 << std::endl;

    // boundary conditions for the r-direction, loop index i
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            aux_u.x[ 0 ][ j ][ k ] = 0.;
            aux_u.x[ im-1 ][ j ][ k ] = 0.;
            aux_v.x[ 0 ][ j ][ k ] = c43 * aux_v.x[ 1 ][ j ][ k ] - c13 * aux_v.x[ 2 ][ j ][ k ];
            aux_v.x[ im-1 ][ j ][ k ] = c43 * aux_v.x[ im-2 ][ j ][ k ] - c13 * aux_v.x[ im-3 ][ j ][ k ];
            aux_w.x[ 0 ][ j ][ k ] = c43 * aux_w.x[ 1 ][ j ][ k ] - c13 * aux_w.x[ 2 ][ j ][ k ];
            aux_w.x[ im-1 ][ j ][ k ] = c43 * aux_w.x[ im-2 ][ j ][ k ] - c13 * aux_w.x[ im-3 ][ j ][ k ];
        }
    }

    double zeta = 3.715;
    // collection of coefficients for coordinate stretching
    double exp_rm = 0.;
    double exp_2_rm = 0.;
    double exp_2_rm_2 = 0.;

    // boundary conditions for the the-direction, loop index j
    for ( int k = 0; k < km; k++ ){
        for ( int i = 0; i < im; i++ ){
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
            aux_u.x[ i ][ 0 ][ k ] = c43 * aux_u.x[ i ][ 1 ][ k ] - c13 * aux_u.x[ i ][ 2 ][ k ];
            aux_u.x[ i ][ jm-1 ][ k ] = c43 * aux_u.x[ i ][ jm-2 ][ k ] - c13 * aux_u.x[ i ][ jm-3 ][ k ];
            aux_v.x[ i ][ 0 ][ k ] = c43 * aux_v.x[ i ][ 1 ][ k ] - c13 * aux_v.x[ i ][ 2 ][ k ];
            aux_v.x[ i ][ jm-1 ][ k ] = c43 * aux_v.x[ i ][ jm-2 ][ k ] - c13 * aux_v.x[ i ][ jm-3 ][ k ];
            aux_w.x[ i ][ 0 ][ k ] = c43 * aux_w.x[ i ][ 1 ][ k ] - c13 * aux_w.x[ i ][ 2 ][ k ];
            aux_w.x[ i ][ jm-1 ][ k ] = c43 * aux_w.x[ i ][ jm-2 ][ k ] - c13 * aux_w.x[ i ][ jm-3 ][ k ];
        }
    }

// boundary conditions for the phi-direction, loop index k
    for ( int i = 0; i < im; i++ ){
        for ( int j = 0; j < jm; j++ ){
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
            aux_u.x[ i ][ j ][ 0 ] = c43 * aux_u.x[ i ][ j ][ 1 ] - c13 * aux_u.x[ i ][ j ][ 2 ];
            aux_u.x[ i ][ j ][ km-1 ] = c43 * aux_u.x[ i ][ j ][ km-2 ] - c13 * aux_u.x[ i ][ j ][ km-3 ];
            aux_u.x[ i ][ j ][ 0 ] = aux_u.x[ i ][ j ][ km-1 ] = ( aux_u.x[ i ][ j ][ 0 ] + aux_u.x[ i ][ j ][ km-1 ] ) / 2.;
            aux_v.x[ i ][ j ][ 0 ] = c43 * aux_v.x[ i ][ j ][ 1 ] - c13 * aux_v.x[ i ][ j ][ 2 ];
            aux_v.x[ i ][ j ][ km-1 ] = c43 * aux_v.x[ i ][ j ][ km-2 ] - c13 * aux_v.x[ i ][ j ][ km-3 ];
            aux_v.x[ i ][ j ][ 0 ] = aux_v.x[ i ][ j ][ km-1 ] = ( aux_v.x[ i ][ j ][ 0 ] + aux_v.x[ i ][ j ][ km-1 ] ) / 2.;
            aux_w.x[ i ][ j ][ 0 ] = c43 * aux_w.x[ i ][ j ][ 1 ] - c13 * aux_w.x[ i ][ j ][ 2 ];
            aux_w.x[ i ][ j ][ km-1 ] = c43 * aux_w.x[ i ][ j ][ km-2 ] - c13 * aux_w.x[ i ][ j ][ km-3 ];
            aux_w.x[ i ][ j ][ 0 ] = aux_w.x[ i ][ j ][ km-1 ] = ( aux_w.x[ i ][ j ][ 0 ] + aux_w.x[ i ][ j ][ km-1 ] ) / 2.;
        }
    }

    double rm= 0;
    double rm2= 0;
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double sinthe = 0;
    double rmsinthe = 0;
    double rm2sinthe2 = 0;
    double denom = 0;
    double num1 = 0;
    double num2 = 0;
    double num3 = 0;
    double drhs_udr = 0;
    double drhs_vdthe = 0;
    double drhs_wdphi = 0;

// Pressure using Euler equation ( 2. derivative of pressure added to the Poisson-right-hand-side )
    for ( int i = 1; i < im-1; i++ ){
        rm = rad.z[ i ];
        rm2 = rm * rm;

// collection of coefficients for coordinate stretching
        exp_rm = 1. / exp( zeta * rm );
        exp_2_rm = 1. / exp( 2. * zeta * rm );
        exp_2_rm_2 = 2. / exp( 2. * zeta * rm );

        for ( int j = 1; j < jm-1; j++ ){
            sinthe = sin( the.z[ j ] );
            rmsinthe = rm * sinthe;
            rm2sinthe2 = rmsinthe * rmsinthe;
            denom = 2. / dr2 * exp_2_rm + 2. / ( rm2 * dthe2 ) + 2. / ( rm2sinthe2 * dphi2 );
            num1 = 1. / dr2 * exp_2_rm;            
            num2 = 1. / ( rm2 * dthe2 );
            num3 = 1. / ( rm2sinthe2 * dphi2 );

// determining RHS values around mountain surface
            for ( int k = 1; k < km-1; k++ ){
// determining RHS-derivatives around mountain surfaces
                drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * exp_rm;
                if ( i <= im - 3 ){
                    if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i+1, j, k ) ) )
                            drhs_udr = ( - 3. * aux_u.x[ i ][ j ][ k ] + 4. * aux_u.x[ i + 1 ][ j ][ k ] - aux_u.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr ) * exp_rm;
                }else  drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i ][ j ][ k ] ) / dr * exp_rm;

// gradients of RHS terms at mountain sides 2.order accurate in the-direction
                drhs_vdthe = ( aux_v.x[ i ][ j+1 ][ k ] - aux_v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe * rm );
                if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j+1, k ) ) )
                                    drhs_vdthe = ( aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j ][ k ] ) / ( dthe * rm );
                if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j-1, k ) ) )
                                    drhs_vdthe = ( aux_v.x[ i ][ j - 1 ][ k ] - aux_v.x[ i ][ j ][ k ] ) / ( dthe * rm);
                if ( ( j >= 2 ) && ( j <= jm - 3 ) ){
                    if ( ( is_land ( h, i, j, k ) ) && ( ( is_air ( h, i, j+1, k ) ) && ( is_air ( h, i, j+2, k ) ) ) )
                        drhs_vdthe = ( - 3. * aux_v.x[ i ][ j ][ k ] + 4. * aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe * rm );
                    if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j-1, k ) ) && ( is_air ( h, i, j-2, k ) ) )
                        drhs_vdthe = ( - 3. * aux_v.x[ i ][ j ][ k ] + 4. * aux_v.x[ i ][ j - 1 ][ k ] - aux_v.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe * rm );
                }

// gradients of RHS terms at mountain sides 2.order accurate in phi-direction
                drhs_wdphi = ( aux_w.x[ i ][ j ][ k+1 ] - aux_w.x[ i ][ j ][ k-1 ] ) / ( 2. * rmsinthe * dphi );
                if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k+1 ) ) )
                                    drhs_wdphi = ( aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k ] ) / ( rmsinthe * dphi );
                if ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j, k-1 ) ) )
                                    drhs_wdphi = ( aux_w.x[ i ][ j ][ k - 1 ] - aux_w.x[ i ][ j ][ k ] ) / ( rmsinthe * dphi );
                if ( ( k >= 2 ) && ( k <= km - 3 ) ){
                    if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k+1 ) ) && ( is_air ( h, i, j, k+2 ) ) )
                        drhs_wdphi = ( - 3. * aux_w.x[ i ][ j ][ k ] + 4. * aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k + 2 ] ) / ( 2. * rmsinthe * dphi );
                    if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k-1 ) ) && ( is_air ( h, i, j, k-2 ) ) )
                        drhs_wdphi = ( - 3. * aux_w.x[ i ][ j ][ k ] + 4. * aux_w.x[ i ][ j ][ k - 1 ] - aux_w.x[ i ][ j ][ k - 2 ] ) / ( 2. * rmsinthe * dphi );
                }

// explicit pressure computation based on the Poisson equation
                p_dyn.x[ i ][ j ][ k ] = ( ( p_dynn.x[ i+1 ][ j ][ k ] + p_dynn.x[ i-1 ][ j ][ k ] ) * num1 * exp_2_rm_2
                                                    + ( p_dynn.x[ i ][ j+1 ][ k ] + p_dynn.x[ i ][ j-1 ][ k ] ) * num2
                                                    + ( p_dynn.x[ i ][ j ][ k+1 ] + p_dynn.x[ i ][ j ][ k-1 ] ) * num3 
                                                    - r_air * ( drhs_udr + drhs_vdthe + drhs_wdphi ) ) / denom;
                if ( is_land ( h, i, j, k ) )  p_dyn.x[ i ][ j ][ k ] = .0;
            }
        }
    }

// boundary conditions for the r-direction, loop index i
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            p_dyn.x[ 0 ][ j ][ k ] = p_dyn.x[ 3 ][ j ][ k ] - 3. * p_dyn.x[ 2 ][ j ][ k ] + 3. * p_dyn.x[ 1 ][ j ][ k ];  // extrapolation
            p_dyn.x[ im-1 ][ j ][ k ] = c43 * p_dyn.x[ im-2 ][ j ][ k ] - c13 * p_dyn.x[ im-3 ][ j ][ k ];        // Neumann
        }
    }

// boundary conditions for the phi-direction, loop index k
    for ( int i = 0; i < im; i++ ){
        for ( int j = 0; j < jm; j++ ){
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
            p_dyn.x[ i ][ j ][ 0 ] = c43 * p_dyn.x[ i ][ j ][ 1 ] - c13 * p_dyn.x[ i ][ j ][ 2 ];
            p_dyn.x[ i ][ j ][ km-1 ] = c43 * p_dyn.x[ i ][ j ][ km-2 ] - c13 * p_dyn.x[ i ][ j ][ km-3 ];
            p_dyn.x[ i ][ j ][ 0 ] = p_dyn.x[ i ][ j ][ km-1 ] = ( p_dyn.x[ i ][ j ][ 0 ] + p_dyn.x[ i ][ j ][ km-1 ] ) / 2.;
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int i = 0; i < im; i++ ){
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
            p_dyn.x[ i ][ 0 ][ k ] = c43 * p_dyn.x[ i ][ 1 ][ k ] - c13 * p_dyn.x[ i ][ 2 ][ k ];
            p_dyn.x[ i ][ jm-1 ][ k ] = c43 * p_dyn.x[ i ][ jm-2 ][ k ] - c13 * p_dyn.x[ i ][ jm-3 ][ k ];
        }
    }


    logger() << "exit $$$$$$$$$$$$$$$ Pressure_Atm::computePressure_3D: p_dyn: " << p_dyn.max() * u_0 * u_0 * r_air *.01 << std::endl << std::endl;

//    circulation.Pressure_Limitation_Atm ( p_dyn, p_dynn );
}

/*
*
*/
void cAtmosphereModel::computePressure_2D(){
    const  int c43 = 4./3., c13 = 1./3.;
    logger() << "enter &&&&&&&&&&&& computePressure_2D: p_dyn: " << p_dyn.max() * u_0 * u_0 * r_air *.01 << std::endl;

    // Pressure using Euler equation ( 2. derivative of pressure added to the Poisson-right-hand-side )
    // boundary conditions for the the-direction, loop index j
    for ( int k = 0; k < km; k++ ){
        // zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
        aux_v.x[ 0 ][ 0 ][ k ] = c43 * aux_v.x[ 0 ][ 1 ][ k ] - c13 * aux_v.x[ 0 ][ 2 ][ k ];
        aux_v.x[ 0 ][ jm-1 ][ k ] = c43 * aux_v.x[ 0 ][ jm-2 ][ k ] - c13 * aux_v.x[ 0 ][ jm-3 ][ k ];
    }

    // boundary conditions for the phi-direction, loop index k
    for ( int j = 0; j < jm; j++ ){
        // zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
        aux_v.x[ 0 ][ j ][ 0 ] = c43 * aux_v.x[ 0 ][ j ][ 1 ] - c13 * aux_v.x[ 0 ][ j ][ 2 ];
        aux_v.x[ 0 ][ j ][ km-1 ] = c43 * aux_v.x[ 0 ][ j ][ km-2 ] - c13 * aux_v.x[ 0 ][ j ][ km-3 ];
        aux_v.x[ 0 ][ j ][ 0 ] = aux_v.x[ 0 ][ j ][ km-1 ] = ( aux_v.x[ 0 ][ j ][ 0 ] + aux_v.x[ 0 ][ j ][ km-1 ] ) / 2.;
    }


    double rm= 0;
    double rm2= 0;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double sinthe = 0;
    double rmsinthe = 0;
    double rm2sinthe2 = 0;
    double denom = 0;
    double num2 = 0;
    double num3 = 0;
    double drhs_vdthe = 0;
    double drhs_wdphi = 0;

    rm = rad.z[ 0 ];
    rm2 = rm * rm;
    for ( int j = 1; j < jm-1; j++ ){
        sinthe = sin( the.z[ j ] );
        rmsinthe = rm * sinthe;
        rm2sinthe2 = rmsinthe * rmsinthe;
        denom = 2. / ( rm2 * dthe2 ) + 2. / ( rm2sinthe2 * dphi2 );
        num2 = 1. / ( rm2 * dthe2 );
        num3 = 1. / ( rm2sinthe2 * dphi2 );
        for ( int k = 1; k < km-1; k++ ){
            // gradients of RHS terms at mountain sides 2.order accurate in the-direction
            drhs_vdthe = ( aux_v.x[ 0 ][ j+1 ][ k ] - aux_v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe * rm );
            if ( is_land( h, 0, j, k ) && is_air( h, 0, j+1, k) ){
                if ( j < jm-2 && is_air( h, 0, j+2, k) ){
                    drhs_vdthe = ( - 3. * aux_v.x[ 0 ][ j ][ k ] + 4. * aux_v.x[ 0 ][ j+1 ][ k ] -
                           aux_v.x[ 0 ][ j+2 ][ k ] ) / ( 2. * dthe * rm );
                }else{
                    drhs_vdthe = ( aux_v.x[ 0 ][ j+1 ][ k ] - aux_v.x[ 0 ][ j ][ k ] ) / ( dthe * rm );
                }
            }
            if ( is_land( h, 0, j, k ) && is_air( h, 0,  j-1,  k) ){
                if ( j > 1 && is_air( h, 0, j-2,  k) ){
                    drhs_vdthe = ( - 3. * aux_v.x[ 0 ][ j ][ k ] + 4. * aux_v.x[ 0 ][ j-1 ][ k ] -
                            aux_v.x[ 0 ][ j-2 ][ k ] ) / ( 2. * dthe * rm );
                }else{
                    drhs_vdthe = ( aux_v.x[ 0 ][ j-1 ][ k ] - aux_v.x[ 0 ][ j ][ k ] ) / ( dthe * rm );
                }
            }

            // gradients of RHS terms at mountain sides 2.order accurate in phi-direction
            drhs_wdphi = ( aux_w.x[ 0 ][ j ][ k+1 ] - aux_w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe );
            if ( is_land( h, 0,  j, k ) && is_air( h, 0, j, k+1) ){
                if ( k < km-2 && is_air(h, 0, j, k+2 )){
                    drhs_wdphi = ( - 3. * aux_w.x[ 0 ][ j ][ k ] + 4. * aux_w.x[ 0 ][ j ][ k+1 ] -
                            aux_w.x[ 0 ][ j ][ k+2 ] ) / ( 2. * rmsinthe * dphi );
                }else{
                    drhs_wdphi = ( aux_w.x[ 0 ][ j ][ k+1 ] - aux_w.x[ 0 ][ j ][ k ] ) / ( dphi * rmsinthe );
                }
            }
            if ( is_land( h, 0, j,  k ) && is_air( h, 0, j, k-1 ) ){
                if ( k >= 2 && is_air( h, 0, j, k-2 ) ){
                    drhs_wdphi = ( - 3. * aux_w.x[ 0 ][ j ][ k ] + 4. * aux_w.x[ 0 ][ j ][ k-1 ] -
                            aux_w.x[ 0 ][ j ][ k-2 ] ) / ( 2. * rmsinthe * dphi );
                }else{
                    drhs_wdphi = ( aux_w.x[ 0 ][ j ][ k-1 ] - aux_w.x[ 0 ][ j ][ k ] ) / ( dphi * rmsinthe );
                }
            }
            p_dyn.x[ 0 ][ j ][ k ] = ( ( p_dynn.x[ 0 ][ j+1 ][ k ] + p_dynn.x[ 0 ][ j-1 ][ k ] ) * num2
                                                + ( p_dynn.x[ 0 ][ j ][ k+1 ] + p_dynn.x[ 0 ][ j ][ k-1 ] ) * num3 
                                                - r_air * ( drhs_vdthe + drhs_wdphi ) ) / denom;
        }
    }

    // boundary conditions for the the-direction, loop index j
    for ( int k = 0; k < km; k++ ){
        // zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
        p_dyn.x[ 0 ][ 0 ][ k ] = 0.;
        p_dyn.x[ 0 ][ jm-1 ][ k ] = 0.;
    }

    // boundary conditions for the phi-direction, loop index k
    for ( int j = 1; j < jm - 1; j++ ){
        // zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
        p_dyn.x[ 0 ][ j ][ 0 ] = c43 * p_dyn.x[ 0 ][ j ][ 1 ] - c13 * p_dyn.x[ 0 ][ j ][ 2 ];
        p_dyn.x[ 0 ][ j ][ km-1 ] = c43 * p_dyn.x[ 0 ][ j ][ km-2 ] - c13 * p_dyn.x[ 0 ][ j ][ km-3 ];
        p_dyn.x[ 0 ][ j ][ 0 ] = p_dyn.x[ 0 ][ j ][ km-1 ] = ( p_dyn.x[ 0 ][ j ][ 0 ] + p_dyn.x[ 0 ][ j ][ km-1 ] ) / 2.;
    }

    logger() << "exit &&&&&&&&&&&&&&& computePressure_2D: p_dyn: " << p_dyn.max() * u_0 * u_0 * r_air *.01 << std::endl << std::endl;

    //    circulation.Pressure_Limitation_Atm ( p_dyn, p_dynn );
}


