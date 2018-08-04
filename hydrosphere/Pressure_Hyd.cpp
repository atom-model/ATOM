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

#include "Pressure_Hyd.h"
#include "Array.h"
#include "Array_2D.h"
#include "MinMax_Atm.h"
#include "Accuracy_Hyd.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


Pressure_Hyd::Pressure_Hyd ( int im, int jm, int km, double dr, double dthe, double dphi )
{
    this-> im = im;
    this-> jm = jm;
    this-> km = km;
    this-> dr = dr;
    this-> dthe = dthe;
    this-> dphi = dphi;


    c43 = 4./3.;
    c13 = 1./3.;
}


Pressure_Hyd::~Pressure_Hyd (){}



void Pressure_Hyd::computePressure_3D ( BC_Thermohalin &oceanflow, double u_0, double r_0_water,
                                 Array_1D &rad, Array_1D &the, Array &p_dyn, Array &p_dynn,
                                 Array &h, Array &aux_u, Array &aux_v, Array &aux_w )
{
// boundary conditions for the r-direction, loop index i

    logger() << "enter computePressure_3D: p_dyn: " << p_dyn.max() * u_0 * u_0 * r_0_water *.01 << std::endl;
    logger() << "enter computePressure_3D: p_dynn: " << p_dynn.max() * u_0 * u_0 * r_0_water *.01 << std::endl << std::endl;

// boundary conditions for the r-direction, loop index i
    for ( int j = 0; j < jm; j++ )
    {
        for ( int k = 0; k < km; k++ )
        {
            aux_u.x[ 0 ][ j ][ k ] = 0.;
            aux_u.x[ im-1 ][ j ][ k ] = 0.;

            aux_v.x[ 0 ][ j ][ k ] = c43 * aux_v.x[ 1 ][ j ][ k ] - c13 * aux_v.x[ 2 ][ j ][ k ];
            aux_v.x[ im-1 ][ j ][ k ] = c43 * aux_v.x[ im-2 ][ j ][ k ] - c13 * aux_v.x[ im-3 ][ j ][ k ];

            aux_w.x[ 0 ][ j ][ k ] = c43 * aux_w.x[ 1 ][ j ][ k ] - c13 * aux_w.x[ 2 ][ j ][ k ];
            aux_w.x[ im-1 ][ j ][ k ] = c43 * aux_w.x[ im-2 ][ j ][ k ] - c13 * aux_w.x[ im-3 ][ j ][ k ];
        }
    }

// boundary conditions for the the-direction, loop index j
    for ( int k = 0; k < km; k++ )
    {
        for ( int i = 0; i < im; i++ )
        {
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
    for ( int i = 0; i < im; i++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
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
    for ( int i = 1; i < im-1; i++ )
    {
        rm = rad.z[ i ];
        rm2 = rm * rm;

        for ( int j = 1; j < jm-1; j++ )
        {
            sinthe = sin( the.z[ j ] );
            rmsinthe = rm * sinthe;
            rm2sinthe2 = rmsinthe * rmsinthe;
            denom = 2. / dr2 + 2. / ( rm2 * dthe2 ) + 2. / ( rm2sinthe2 * dphi2 );
            num1 = 1. / dr2;
            num2 = 1. / ( rm2 * dthe2 );
            num3 = 1. / ( rm2sinthe2 * dphi2 );

// determining RHS values around mountain surface
            for ( int k = 1; k < km-1; k++ )
            {
// determining RHS-derivatives around mountain surfaces
                drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );

                if ( i <= im - 3 )
                {
                    if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i + 1 ][ j ][ k ] == 0. ) )
                            drhs_udr = ( - 3. * aux_u.x[ i ][ j ][ k ] + 4. * aux_u.x[ i + 1 ][ j ][ k ] - aux_u.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );
                }
                else     drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i ][ j ][ k ] ) / dr;



// gradients of RHS terms at mountain sides 2.order accurate in the-direction
            drhs_vdthe = ( aux_v.x[ im-1 ][ j+1 ][ k ] - aux_v.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe * rm );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j + 1 ][ k ] == 0. ) )
                        drhs_vdthe = ( aux_v.x[ im-1 ][ j + 1 ][ k ] - aux_v.x[ im-1 ][ j ][ k ] ) / ( dthe * rm );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j + 1 ][ k ] == 0. ) && ( h.x[ im-1 ][ j + 2 ][ k ] == 0. ) )
                        drhs_vdthe = ( - 3. * aux_v.x[ im-1 ][ j ][ k ] + 4. * aux_v.x[ im-1 ][ j + 1 ][ k ] - aux_v.x[ im-1 ][ j + 2 ][ k ] ) / ( 2. * dthe * rm );


            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j - 1 ][ k ] == 0. ) && ( h.x[ im-1 ][ j - 2 ][ k ] == 0. ) )
                        drhs_vdthe = ( aux_v.x[ im-1 ][ j - 1 ][ k ] - aux_v.x[ im-1 ][ j ][ k ] ) / ( dthe * rm );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j - 1 ][ k ] == 0. ) )
                        drhs_vdthe = ( - 3. * aux_v.x[ im-1 ][ j ][ k ] + 4. * aux_v.x[ im-1 ][ j - 1 ][ k ] - aux_v.x[ im-1 ][ j - 2 ][ k ] ) / ( 2. * dthe * rm );


            if ( j == 1 )  drhs_vdthe = ( - 3. * aux_v.x[ im-1 ][ j ][ k ] + 4. * aux_v.x[ im-1 ][ j + 1 ][ k ] - aux_v.x[ im-1 ][ j + 2 ][ k ] ) / ( 2. * dthe * rm );

            if ( j == jm - 2 )  drhs_vdthe = ( - 3. * aux_v.x[ im-1 ][ j ][ k ] + 4. * aux_v.x[ im-1 ][ j - 1 ][ k ] - aux_v.x[ im-1 ][ j - 2 ][ k ] ) / ( 2. * dthe * rm );


// gradients of RHS terms at mountain sides 2.order accurate in phi-direction
            drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k+1 ] - aux_w.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k + 1 ] == 0. ) )
                        drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k + 1 ] - aux_w.x[ im-1 ][ j ][ k ] ) / ( dphi * rmsinthe );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k + 1 ] == 0. ) && ( h.x[ im-1 ][ j ][ k + 2 ] == 0. ) )
                        drhs_wdphi = ( - 3. * aux_w.x[ im-1 ][ j ][ k ] + 4. * aux_w.x[ im-1 ][ j ][ k + 1 ] - aux_w.x[ im-1 ][ j ][ k + 2 ] ) / ( 2. * rmsinthe * dphi );


            if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h.x[ im-1 ][ j ][ k - 1 ] == 1. ) )
                        drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k - 1 ] - aux_w.x[ im-1 ][ j ][ k ] ) / ( dphi * rmsinthe );

            if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h.x[ im-1 ][ j ][ k - 1 ] == 1. ) && ( h.x[ im-1 ][ j ][ k - 2 ] == 1. ) )
                        drhs_wdphi = ( - 3. * aux_w.x[ im-1 ][ j ][ k ] + 4. * aux_w.x[ im-1 ][ j ][ k - 1 ] - aux_w.x[ im-1 ][ j ][ k - 2 ] ) / ( 2. * rmsinthe * dphi );


            if ( k == 2 )  drhs_wdphi = ( - 3. * aux_w.x[ im-1 ][ j ][ k ] + 4. * aux_w.x[ im-1 ][ j ][ k + 1 ] - aux_w.x[ im-1 ][ j ][ k + 2 ] ) / ( 2. * rmsinthe * dphi );

            if ( k == km - 2 )  drhs_wdphi = ( - 3. * aux_w.x[ im-1 ][ j ][ k ] + 4. * aux_w.x[ im-1 ][ j ][ k - 1 ] - aux_w.x[ im-1 ][ j ][ k - 2 ] ) / ( 2. * rmsinthe * dphi );


// explicit pressure computation based on the Poisson equation
                p_dyn.x[ i ][ j ][ k ] = ( ( p_dynn.x[ i+1 ][ j ][ k ] + p_dynn.x[ i-1 ][ j ][ k ] ) * num1
                                                    + ( p_dynn.x[ i ][ j+1 ][ k ] + p_dynn.x[ i ][ j-1 ][ k ] ) * num2
                                                    + ( p_dynn.x[ i ][ j ][ k+1 ] + p_dynn.x[ i ][ j ][ k-1 ] ) * num3 
                                                    - r_0_water * ( drhs_udr + drhs_vdthe + drhs_wdphi ) ) / denom;

                if ( h.x[ i ][ j ][ k ] == 1. )                        p_dyn.x[ i ][ j ][ k ] = .0;
            }
        }
    }


// boundary conditions for the r-direction, loop index i
    for ( int j = 0; j < jm; j++ )
    {
        for ( int k = 0; k < km; k++ )
        {
            p_dyn.x[ 0 ][ j ][ k ] = p_dyn.x[ 3 ][ j ][ k ] - 3. * p_dyn.x[ 2 ][ j ][ k ] + 3. * p_dyn.x[ 1 ][ j ][ k ];  // extrapolation
            p_dyn.x[ im-1 ][ j ][ k ] = c43 * p_dyn.x[ im-2 ][ j ][ k ] - c13 * p_dyn.x[ im-3 ][ j ][ k ];        // Neumann
        }
    }

// boundary conditions for the phi-direction, loop index k
    for ( int i = 0; i < im; i++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
            p_dyn.x[ i ][ j ][ 0 ] = c43 * p_dyn.x[ i ][ j ][ 1 ] - c13 * p_dyn.x[ i ][ j ][ 2 ];
            p_dyn.x[ i ][ j ][ km-1 ] = c43 * p_dyn.x[ i ][ j ][ km-2 ] - c13 * p_dyn.x[ i ][ j ][ km-3 ];
            p_dyn.x[ i ][ j ][ 0 ] = p_dyn.x[ i ][ j ][ km-1 ] = ( p_dyn.x[ i ][ j ][ 0 ] + p_dyn.x[ i ][ j ][ km-1 ] ) / 2.;
        }
    }

    for ( int k = 0; k < km; k++ )
    {
        for ( int i = 0; i < im; i++ )
        {
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
            p_dyn.x[ i ][ 0 ][ k ] = c43 * p_dyn.x[ i ][ 1 ][ k ] - c13 * p_dyn.x[ i ][ 2 ][ k ];
            p_dyn.x[ i ][ jm-1 ][ k ] = c43 * p_dyn.x[ i ][ jm-2 ][ k ] - c13 * p_dyn.x[ i ][ jm-3 ][ k ];
        }
    }


    logger() << "exit computePressure_3D: p_dyn: " << p_dyn.max() * u_0 * u_0 * r_0_water *.01 << std::endl;
    logger() << "exit computePressure_3D: p_dynn: " << p_dynn.max() * u_0 * u_0 * r_0_water *.01 << std::endl << std::endl;

//    oceanflow.Pressure_Limitation_Hyd ( p_dyn, p_dynn );
}









void Pressure_Hyd::computePressure_2D ( BC_Thermohalin &oceanflow, double u_0, double r_0_water,
                                 Array_1D &rad, Array_1D &the, Array &p_dyn, Array &p_dynn,
                                 Array &h, Array &aux_v, Array &aux_w )
{
// Pressure using Euler equation ( 2. derivative of pressure added to the Poisson-right-hand-side )
    logger() << "enter computePressure_2D: p_dyn: " << p_dyn.max() * u_0 * u_0 * r_0_water *.01 << std::endl;
    logger() << "enter computePressure_2D: p_dynn: " << p_dynn.max() * u_0 * u_0 * r_0_water *.01 << std::endl << std::endl;

// boundary conditions for the the-direction, loop index j
    for ( int k = 0; k < km; k++ )
    {
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
        aux_v.x[ im-1 ][ 0 ][ k ] = c43 * aux_v.x[ im-1 ][ 1 ][ k ] - c13 * aux_v.x[ im-1 ][ 2 ][ k ];
        aux_v.x[ im-1 ][ jm-1 ][ k ] = c43 * aux_v.x[ im-1 ][ jm-2 ][ k ] - c13 * aux_v.x[ im-1 ][ jm-3 ][ k ];
    }

// boundary conditions for the phi-direction, loop index k
    for ( int j = 0; j < jm; j++ )
    {
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
        aux_v.x[ im-1 ][ j ][ 0 ] = c43 * aux_v.x[ im-1 ][ j ][ 1 ] - c13 * aux_v.x[ im-1 ][ j ][ 2 ];
        aux_v.x[ im-1 ][ j ][ km-1 ] = c43 * aux_v.x[ im-1 ][ j ][ km-2 ] - c13 * aux_v.x[ im-1 ][ j ][ km-3 ];
        aux_v.x[ im-1 ][ j ][ 0 ] = aux_v.x[ im-1 ][ j ][ km-1 ] = ( aux_v.x[ im-1 ][ j ][ 0 ] + aux_v.x[ im-1 ][ j ][ km-1 ] ) / 2.;
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


    rm = rad.z[ im-1 ];
    rm2 = rm * rm;

    for ( int j = 1; j < jm-1; j++ )
    {
        sinthe = sin( the.z[ j ] );
        rmsinthe = rm * sinthe;
        rm2sinthe2 = rmsinthe * rmsinthe;
        denom = 2. / ( rm2 * dthe2 ) + 2. / ( rm2sinthe2 * dphi2 );
        num2 = 1. / ( rm2 * dthe2 );
        num3 = 1. / ( rm2sinthe2 * dphi2 );

        for ( int k = 1; k < km-1; k++ )
        {
// gradients of RHS terms at mountain sides 2.order accurate in the-direction
            drhs_vdthe = ( aux_v.x[ im-1 ][ j+1 ][ k ] - aux_v.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe * rm );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j+1 ][ k ] == 0. ) )
                        drhs_vdthe = ( aux_v.x[ im-1 ][ j+1 ][ k ] - aux_v.x[ im-1 ][ j ][ k ] ) / ( dthe * rm );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j+1 ][ k ] == 0. ) && ( h.x[ im-1 ][ j+2 ][ k ] == 0. ) )
                        drhs_vdthe = ( - 3. * aux_v.x[ im-1 ][ j ][ k ] + 4. * aux_v.x[ im-1 ][ j+1 ][ k ] -
                        aux_v.x[ im-1 ][ j+2 ][ k ] ) / ( 2. * dthe * rm );


            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j-1 ][ k ] == 0. ) && ( h.x[ im-1 ][ j-2 ][ k ] == 0. ) )
                        drhs_vdthe = ( aux_v.x[ im-1 ][ j - 1 ][ k ] - aux_v.x[ im-1 ][ j ][ k ] ) / ( dthe * rm );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j-1 ][ k ] == 0. ) )
                        drhs_vdthe = ( - 3. * aux_v.x[ im-1 ][ j ][ k ] + 4. * aux_v.x[ im-1 ][ j-1 ][ k ] -
                        aux_v.x[ im-1 ][ j-2 ][ k ] ) / ( 2. * dthe * rm );


            if ( j == 1 )  drhs_vdthe = ( - 3. * aux_v.x[ im-1 ][ j ][ k ] + 4. * aux_v.x[ im-1 ][ j+1 ][ k ] -
                                                        aux_v.x[ im-1 ][ j+2 ][ k ] ) / ( 2. * dthe * rm );

            if ( j == jm-2 )  drhs_vdthe = ( - 3. * aux_v.x[ im-1 ][ j ][ k ] + 4. * aux_v.x[ im-1 ][ j-1 ][ k ] -
                                                        aux_v.x[ im-1 ][ j-2 ][ k ] ) / ( 2. * dthe * rm );


// gradients of RHS terms at mountain sides 2.order accurate in phi-direction
            drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k+1 ] - aux_w.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k+1 ] == 0. ) )
                        drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k+1 ] - aux_w.x[ im-1 ][ j ][ k ] ) / ( dphi * rmsinthe );

            if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k+1 ] == 0. ) && ( h.x[ im-1 ][ j ][ k+2 ] == 0. ) )
                        drhs_wdphi = ( - 3. * aux_w.x[ im-1 ][ j ][ k ] + 4. * aux_w.x[ im-1 ][ j ][ k+1 ] -
                                                 aux_w.x[ im-1 ][ j ][ k+2 ] ) / ( 2. * rmsinthe * dphi );


            if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h.x[ im-1 ][ j ][ k-1 ] == 1. ) )
                        drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k-1 ] - aux_w.x[ im-1 ][ j ][ k ] ) / ( dphi * rmsinthe );

            if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h.x[ im-1 ][ j ][ k-1 ] == 1. ) && ( h.x[ im-1 ][ j ][ k-2 ] == 1. ) )
                        drhs_wdphi = ( - 3. * aux_w.x[ im-1 ][ j ][ k ] + 4. * aux_w.x[ im-1 ][ j ][ k-1 ] -
                                                 aux_w.x[ im-1 ][ j ][ k-2 ] ) / ( 2. * rmsinthe * dphi );


            if ( k == 2 )  drhs_wdphi = ( - 3. * aux_w.x[ im-1 ][ j ][ k ] + 4. * aux_w.x[ im-1 ][ j ][ k+1 ] -
                                                             aux_w.x[ im-1 ][ j ][ k+2 ] ) / ( 2. * rmsinthe * dphi );

            if ( k == km-2 )  drhs_wdphi = ( - 3. * aux_w.x[ im-1 ][ j ][ k ] + 4. * aux_w.x[ im-1 ][ j ][ k-1 ] -
                                                                 aux_w.x[ im-1 ][ j ][ k-2 ] ) / ( 2. * rmsinthe * dphi );


            p_dyn.x[ im-1 ][ j ][ k ] = ( ( p_dynn.x[ im-1 ][ j+1 ][ k ] + p_dynn.x[ im-1 ][ j-1 ][ k ] ) * num2
                                                + ( p_dynn.x[ im-1 ][ j ][ k+1 ] + p_dynn.x[ im-1 ][ j ][ k-1 ] ) * num3 
                                                - r_0_water * ( drhs_vdthe + drhs_wdphi ) ) / denom;

            if ( h.x[ im-1 ][ j ][ k ] == 1. )                        p_dyn.x[ im-1 ][ j ][ k ] = .0;
        }
    }


// boundary conditions for the the-direction, loop index j
    for ( int k = 0; k < km; k++ )
    {
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
        p_dyn.x[ im-1 ][ im-1 ][ k ] = 0.;
        p_dyn.x[ 0 ][ jm-1 ][ k ] = 0.;
    }

// boundary conditions for the phi-direction, loop index k
    for ( int j = 1; j < jm - 1; j++ )
    {
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
        p_dyn.x[ im-1 ][ j ][ 0 ] = c43 * p_dyn.x[ im-1 ][ j ][ 1 ] - c13 * p_dyn.x[ im-1 ][ j ][ 2 ];
        p_dyn.x[ im-1 ][ j ][ km-1 ] = c43 * p_dyn.x[ im-1 ][ j ][ km-2 ] - c13 * p_dyn.x[ im-1 ][ j ][ km-3 ];
        p_dyn.x[ im-1 ][ j ][ 0 ] = p_dyn.x[ im-1 ][ j ][ km-1 ] = ( p_dyn.x[ im-1 ][ j ][ 0 ] + p_dyn.x[ im-1 ][ j ][ km-1 ] ) / 2.;
    }


//    oceanflow.Pressure_Limitation_Hyd ( p_dyn, p_dynn );
    logger() << "exit computePressure_2D: p_dyn: " << p_dyn.max() * u_0 * u_0 * r_0_water *.01 << std::endl;
    logger() << "exit computePressure_2D: p_dynn: " << p_dynn.max() * u_0 * u_0 * r_0_water *.01 << std::endl << std::endl;
}
