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

using namespace std;


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




void Pressure_Hyd::computePressure_3D ( double r_0_water, double pa, Array_1D &rad, Array_1D &the, Array &p_dyn, Array &p_dynn, Array &h, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &aux_u, Array &aux_v, Array &aux_w )
{
// boundary conditions for the r-direction, loop index i
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			rhs_u.x[ 0 ][ j ][ k ] = 0.;				// Dirichlet
			rhs_u.x[ im-1 ][ j ][ k ] = 0.;		// Dirichlet

			rhs_v.x[ 0 ][ j ][ k ] = c43 * rhs_v.x[ 1 ][ j ][ k ] - c13 * rhs_v.x[ 2 ][ j ][ k ];									// von Neumann
			rhs_v.x[ im-1 ][ j ][ k ] = c43 * rhs_v.x[ im-2 ][ j ][ k ] - c13 * rhs_v.x[ im-3 ][ j ][ k ];		// Neumann

			rhs_w.x[ 0 ][ j ][ k ] = c43 * rhs_w.x[ 1 ][ j ][ k ] - c13 * rhs_w.x[ 2 ][ j ][ k ];									// von Neumann
			rhs_w.x[ im-1 ][ j ][ k ] = c43 * rhs_w.x[ im-2 ][ j ][ k ] - c13 * rhs_w.x[ im-3 ][ j ][ k ];		// Neumann

			p_dyn.x[ 0 ][ j ][ k ] = p_dyn.x[ 3 ][ j ][ k ] - 3. * p_dyn.x[ 2 ][ j ][ k ] + 3. * p_dyn.x[ 1 ][ j ][ k ];		// extrapolation
			p_dyn.x[ im-1 ][ j ][ k ] = p_dyn.x[ im-4 ][ j ][ k ] - 3. * p_dyn.x[ im-3 ][ j ][ k ] + 3. * p_dyn.x[ im-2 ][ j ][ k ];		// extrapolation
		}
	}

// boundary conditions for the the-direction, loop index j
	for ( int k = 0; k < km; k++ )
	{
		for ( int i = 0; i < im; i++ )
		{
			rhs_u.x[ i ][ 0 ][ k ] = c43 * rhs_u.x[ i ][ 1 ][ k ] - c13 * rhs_u.x[ i ][ 2 ][ k ];
			rhs_u.x[ i ][ jm-1 ][ k ] = c43 * rhs_u.x[ i ][ jm-2 ][ k ] - c13 * rhs_u.x[ i ][ jm-3 ][ k ];

			rhs_v.x[ i ][ 0 ][ k ] = c43 * rhs_v.x[ i ][ 1 ][ k ] - c13 * rhs_v.x[ i ][ 2 ][ k ];
			rhs_v.x[ i ][ jm-1 ][ k ] = c43 * rhs_v.x[ i ][ jm-2 ][ k ] - c13 * rhs_v.x[ i ][ jm-3 ][ k ];

			rhs_w.x[ i ][ 0 ][ k ] = c43 * rhs_w.x[ i ][ 1 ][ k ] - c13 * rhs_w.x[ i ][ 2 ][ k ];
			rhs_w.x[ i ][ jm-1 ][ k ] = c43 * rhs_w.x[ i ][ jm-2 ][ k ] - c13 * rhs_w.x[ i ][ jm-3 ][ k ];

			p_dyn.x[ i ][ 0 ][ k ] = c43 * p_dyn.x[ i ][ 1 ][ k ] - c13 * p_dyn.x[ i ][ 2 ][ k ];
			p_dyn.x[ i ][ jm-1 ][ k ] = c43 * p_dyn.x[ i ][ jm-2 ][ k ] - c13 * p_dyn.x[ i ][ jm-3 ][ k ];
		}
	}

// boundary conditions for the phi-direction, loop index k
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			rhs_u.x[ i ][ j ][ 0 ] = c43 * rhs_u.x[ i ][ j ][ 1 ] - c13 * rhs_u.x[ i ][ j ][ 2 ];
			rhs_u.x[ i ][ j ][ km-1 ] = c43 * rhs_u.x[ i ][ j ][ km-2 ] - c13 * rhs_u.x[ i ][ j ][ km-3 ];
			rhs_u.x[ i ][ j ][ 0 ] = rhs_u.x[ i ][ j ][ km-1 ] = ( rhs_u.x[ i ][ j ][ 0 ] + rhs_u.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_v.x[ i ][ j ][ 0 ] = c43 * rhs_v.x[ i ][ j ][ 1 ] - c13 * rhs_v.x[ i ][ j ][ 2 ];
			rhs_v.x[ i ][ j ][ km-1 ] = c43 * rhs_v.x[ i ][ j ][ km-2 ] - c13 * rhs_v.x[ i ][ j ][ km-3 ];
			rhs_v.x[ i ][ j ][ 0 ] = rhs_v.x[ i ][ j ][ km-1 ] = ( rhs_v.x[ i ][ j ][ 0 ] + rhs_v.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_w.x[ i ][ j ][ 0 ] = c43 * rhs_w.x[ i ][ j ][ 1 ] - c13 * rhs_w.x[ i ][ j ][ 2 ];
			rhs_w.x[ i ][ j ][ km-1 ] = c43 * rhs_w.x[ i ][ j ][ km-2 ] - c13 * rhs_w.x[ i ][ j ][ km-3 ];
			rhs_w.x[ i ][ j ][ 0 ] = rhs_w.x[ i ][ j ][ km-1 ] = ( rhs_w.x[ i ][ j ][ 0 ] + rhs_w.x[ i ][ j ][ km-1 ] ) / 2.;

			p_dyn.x[ i ][ j ][ 0 ] = c43 * p_dyn.x[ i ][ j ][ 1 ] - c13 * p_dyn.x[ i ][ j ][ 2 ];
			p_dyn.x[ i ][ j ][ km-1 ] = c43 * p_dyn.x[ i ][ j ][ km-2 ] - c13 * p_dyn.x[ i ][ j ][ km-3 ];
			p_dyn.x[ i ][ j ][ 0 ] = p_dyn.x[ i ][ j ][ km-1 ] = ( p_dyn.x[ i ][ j ][ 0 ] + p_dyn.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}

// boundary conditions for the r-direction, loop index i
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			aux_u.x[ 0 ][ j ][ k ] = 0.;									// von Neumann
			aux_u.x[ im-1 ][ j ][ k ] = 0.;		// Neumann

			aux_v.x[ 0 ][ j ][ k ] = c43 * aux_v.x[ 1 ][ j ][ k ] - c13 * aux_v.x[ 2 ][ j ][ k ];									// von Neumann
			aux_v.x[ im-1 ][ j ][ k ] = c43 * aux_v.x[ im-2 ][ j ][ k ] - c13 * aux_v.x[ im-3 ][ j ][ k ];		// Neumann

			aux_w.x[ 0 ][ j ][ k ] = c43 * aux_w.x[ 1 ][ j ][ k ] - c13 * aux_w.x[ 2 ][ j ][ k ];									// von Neumann
			aux_w.x[ im-1 ][ j ][ k ] = c43 * aux_w.x[ im-2 ][ j ][ k ] - c13 * aux_w.x[ im-3 ][ j ][ k ];		// Neumann
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


	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

	iter_prec = 0;
	while ( iter_prec <= 20 )																// iter_prec may be varied
	{
		iter_prec = iter_prec + 1;

// Pressure using Euler equation ( 2. derivative of pressure added to the Poisson-right-hand-side )
		for ( int i = 1; i < im-1; i++ )
		{
			rm = rad.z[ i ];
			rm2 = rm * rm;

			for ( int j = 1; j < jm-1; j++ )
			{
				sinthe = sin( the.z[ j ] );
				sinthe2 = sinthe * sinthe;
				costhe = cos( the.z[ j ] );
				cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
				rmsinthe = rm * sinthe;
				denom = 2. / dr2 + 2. / ( rm * dthe2 ) + 2. / ( rmsinthe * dphi2 );
				num1 = 1. / dr2;
				num2 = 1. / ( rm * dthe2 );
				num3 = 1. / ( rmsinthe * dphi );

// gradients of RHS terms at mountain summits 2.order accurate in r-direction
				for ( int k = 1; k < km-1; k++ )
				{
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j + 1 ][ k ] == 0. ) )			aux_v.x[ i ][ j ][ k ] = aux_v.x[ i ][ j + 1 ][ k ];
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) )				aux_v.x[ i ][ j ][ k ] = aux_v.x[ i ][ j - 1 ][ k ];

					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) )			aux_w.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k + 1 ];
					if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h.x[ i ][ j ][ k - 1 ] == 1. ) )				aux_w.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k - 1 ];

					if ( ( j >= 2 ) && ( j < jm - 3 ) )
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( ( h.x[ i ][ j + 1 ][ k ] == 0. ) && ( h.x[ i ][ j + 2 ][ k ] == 0. ) ) )
						{
							aux_v.x[ i ][ j ][ k ] = c43 * aux_v.x[ i ][ j + 1 ][ k ] - c13 * aux_v.x[ i ][ j + 2 ][ k ];
						}
						else		aux_v.x[ i ][ j ][ k ] = aux_v.x[ i ][ j + 1 ][ k ];

						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) && ( h.x[ i ][ j - 2 ][ k ] == 0. ) )
						{
							aux_v.x[ i ][ j ][ k ] = c43 * aux_v.x[ i ][ j - 1 ][ k ] - c13 * aux_v.x[ i ][ j - 2 ][ k ];
						}
						else		aux_v.x[ i ][ j ][ k ] = aux_v.x[ i ][ j - 1 ][ k ];
					}

					if ( ( k >= 2 ) && ( k < km - 3 ) )
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) && ( h.x[ i ][ j ][ k + 2 ] == 0. ) )
						{
							aux_w.x[ i ][ j ][ k ] = c43 * aux_w.x[ i ][ j ][ k + 1 ] - c13 * aux_w.x[ i ][ j ][ k + 2 ];
						}
						else		aux_w.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k + 1 ];

						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k - 1 ] == 0. ) && ( h.x[ i ][ j ][ k - 2 ] == 0. ) )
						{
							aux_w.x[ i ][ j ][ k ] = c43 * aux_w.x[ i ][ j ][ k - 1 ] - c13 * aux_w.x[ i ][ j ][ k - 2 ];
						}
						else		aux_w.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k - 1 ];
					}

// determining RHS-derivatives around mountain surfaces
					drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );

					if ( i < im - 2 )
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i + 1 ][ j ][ k ] == 0. ) )			drhs_udr = ( - 3. * aux_u.x[ i ][ j ][ k ] + 4. * aux_u.x[ i + 1 ][ j ][ k ] - aux_u.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );			// 2. order accurate
					}
					else 		drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i ][ j ][ k ] ) / dr;

// gradients of RHS terms at mountain sides 2.order accurate in the-direction
					drhs_vdthe = ( aux_v.x[ i ][ j+1 ][ k ] - aux_v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );

					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j + 1 ][ k ] == 0. ) )			drhs_vdthe = ( aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j ][ k ] ) / dthe ;					// 1. order accurate
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) )				drhs_vdthe = ( aux_v.x[ i ][ j ][ k ] - aux_v.x[ i ][ j - 1 ][ k ] ) / dthe;					// 1. order accurate

					if ( ( j >= 2 ) && ( j < jm - 3 ) )
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( ( h.x[ i ][ j + 1 ][ k ] == 0. ) && ( h.x[ i ][ j + 2 ][ k ] == 0. ) ) )
						{
							drhs_vdthe = ( - 3. * aux_v.x[ i ][ j ][ k ] + 4. * aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
						}
						else 			drhs_vdthe = ( aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j ][ k ] ) / dthe;

						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) && ( h.x[ i ][ j - 2 ][ k ] == 0. ) )
						{
							drhs_vdthe = ( - 3. * aux_v.x[ i ][ j ][ k ] + 4. * aux_v.x[ i ][ j - 1 ][ k ] - aux_v.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
						}
						else			drhs_vdthe = ( aux_v.x[ i ][ j ][ k ] - aux_v.x[ i ][ j - 1 ][ k ] ) / dthe;
					}

// gradients of RHS terms at mountain sides 2.order accurate in phi-direction
					drhs_wdphi = ( aux_w.x[ i ][ j ][ k+1 ] - aux_w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );

					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) )			drhs_wdphi = ( aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k ] ) / dphi;					// 1. order accurate
					if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h.x[ i ][ j ][ k - 1 ] == 1. ) )				drhs_wdphi = ( aux_w.x[ i ][ j ][ k ] - aux_w.x[ i ][ j ][ k - 1 ] ) / dphi;					// 1. order accurate

					if ( ( k >= 2 ) && ( k < km - 3 ) )
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) && ( h.x[ i ][ j ][ k + 2 ] == 0. ) )
						{
							drhs_wdphi = ( - 3. * aux_w.x[ i ][ j ][ k ] + 4. * aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
						}
						else 			drhs_wdphi = ( aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k ] ) / dphi;

						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k - 1 ] == 0. ) && ( h.x[ i ][ j ][ k - 2 ] == 0. ) )
						{
							drhs_wdphi = ( - 3. * aux_w.x[ i ][ j ][ k ] + 4. * aux_w.x[ i ][ j ][ k - 1 ] - aux_w.x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
						}
						else			drhs_wdphi = ( aux_w.x[ i ][ j ][ k ] - aux_w.x[ i ][ j ][ k - 1 ] ) / dphi;
					}


// explicit pressure computation based on the Poisson equation
						p_dynn.x[ i ][ j ][ k ] = ( ( p_dyn.x[ i+1 ][ j ][ k ] + p_dyn.x[ i-1 ][ j ][ k ] ) * num1
															+ ( p_dyn.x[ i ][ j+1 ][ k ] + p_dyn.x[ i ][ j-1 ][ k ] ) * num2
															+ ( p_dyn.x[ i ][ j ][ k+1 ] + p_dyn.x[ i ][ j ][ k-1 ] ) * num3 
															+ r_0_water * ( drhs_udr + drhs_vdthe + drhs_wdphi ) ) / denom;

						if ( h.x[ i ][ j ][ k ] == 1. )				p_dynn.x[ i ][ j ][ k ] = .0;
				}
			}
		}


// boundary conditions for the r-direction, loop index i
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				p_dynn.x[ 0 ][ j ][ k ] = p_dynn.x[ 3 ][ j ][ k ] - 3. * p_dynn.x[ 2 ][ j ][ k ] + 3. * p_dynn.x[ 1 ][ j ][ k ];		// extrapolation
				p_dynn.x[ im-1 ][ j ][ k ] = c43 * p_dynn.x[ im-2 ][ j ][ k ] - c13 * p_dynn.x[ im-3 ][ j ][ k ];		// Neumann
			}
		}

// boundary conditions for the phi-direction, loop index k
		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
				p_dynn.x[ i ][ j ][ 0 ] = c43 * p_dynn.x[ i ][ j ][ 1 ] - c13 * p_dynn.x[ i ][ j ][ 2 ];
				p_dynn.x[ i ][ j ][ km-1 ] = c43 * p_dynn.x[ i ][ j ][ km-2 ] - c13 * p_dynn.x[ i ][ j ][ km-3 ];
				p_dynn.x[ i ][ j ][ 0 ] = p_dynn.x[ i ][ j ][ km-1 ] = ( p_dynn.x[ i ][ j ][ 0 ] + p_dynn.x[ i ][ j ][ km-1 ] ) / 2.;
			}
		}

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
	// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
				p_dynn.x[ i ][ 0 ][ k ] = c43 * p_dynn.x[ i ][ 1 ][ k ] - c13 * p_dynn.x[ i ][ 2 ][ k ];
				p_dynn.x[ i ][ jm-1 ][ k ] = c43 * p_dynn.x[ i ][ jm-2 ][ k ] - c13 * p_dynn.x[ i ][ jm-3 ][ k ];
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ];
					if ( h.x[ i ][ j ][ k ] == 1. )				p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ] = .0;
				}
			}
		}
	}
}









void Pressure_Hyd::computePressure_2D ( double r_0_water, Array_1D &rad, Array_1D &the, Array &p_dyn, Array &p_dynn, Array &h, Array &rhs_v, Array &rhs_w, Array &aux_v, Array &aux_w )
{
// Pressure using Euler equation ( 2. derivative of pressure added to the Poisson-right-hand-side )

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

// boundary conditions for the the-direction, loop index j
	for ( int k = 0; k < km; k++ )
	{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
		aux_w.x[ im-1 ][ 0 ][ k ] = c43 * aux_w.x[ im-1 ][ 1 ][ k ] - c13 * aux_w.x[ im-1 ][ 2 ][ k ];
		aux_w.x[ im-1 ][ jm-1 ][ k ] = c43 * aux_w.x[ im-1 ][ jm-2 ][ k ] - c13 * aux_w.x[ im-1 ][ jm-3 ][ k ];
	}

// boundary conditions for the phi-direction, loop index k
	for ( int j = 0; j < jm; j++ )
	{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
		aux_w.x[ im-1 ][ j ][ 0 ] = c43 * aux_w.x[ im-1 ][ j ][ 1 ] - c13 * aux_w.x[ im-1 ][ j ][ 2 ];
		aux_w.x[ im-1 ][ j ][ km-1 ] = c43 * aux_w.x[ im-1 ][ j ][ km-2 ] - c13 * aux_w.x[ im-1 ][ j ][ km-3 ];
		aux_w.x[ im-1 ][ j ][ 0 ] = aux_w.x[ im-1 ][ j ][ km-1 ] = ( aux_w.x[ im-1 ][ j ][ 0 ] + aux_w.x[ im-1 ][ j ][ km-1 ] ) / 2.;
	}


	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

	iter_prec = 0;
	while ( iter_prec <= 20 )																// iter_prec may be varied
	{
		iter_prec = iter_prec + 1;

		rm = rad.z[ 0 ];
		rm2 = rm * rm;

		for ( int j = 1; j < jm-1; j++ )
		{
			sinthe = sin( the.z[ j ] );
			rmsinthe = rm * sinthe;
			denom = 2. / ( rm * dthe2 ) + 2. / ( rmsinthe * dphi2 );
			num2 = 1. / ( rm * dthe2 );
			num3 = 1. / ( rmsinthe * dphi );

			for ( int k = 1; k < km-1; k++ )
			{
				if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j + 1 ][ k ] == 0. ) )			aux_v.x[ im-1 ][ j ][ k ] = aux_v.x[ im-1 ][ j + 1 ][ k ];
				if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j - 1 ][ k ] == 0. ) )			aux_v.x[ im-1 ][ j ][ k ] = aux_v.x[ im-1 ][ j - 1 ][ k ];

				if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k + 1 ] == 0. ) )			aux_w.x[ im-1 ][ j ][ k ] = aux_w.x[ im-1 ][ j ][ k + 1 ];
				if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h.x[ im-1 ][ j ][ k - 1 ] == 1. ) )			aux_w.x[ im-1 ][ j ][ k ] = aux_w.x[ im-1 ][ j ][ k - 1 ];

// gradients of RHS terms at mountain sides 2.order accurate in the-direction
				drhs_vdthe = ( aux_v.x[ im-1 ][ j+1 ][ k ] - aux_v.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );

				if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j + 1 ][ k ] == 0. ) )				drhs_vdthe = ( aux_v.x[ im-1 ][ j + 1 ][ k ] - aux_v.x[ im-1 ][ j ][ k ] ) / dthe;					// 1. order accurate
				if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j - 1 ][ k ] == 0. ) )				drhs_vdthe = ( aux_v.x[ im-1 ][ j ][ k ] - aux_v.x[ im-1 ][ j - 1 ][ k ] ) / dthe;					// 1. order accurate

// gradients of RHS terms at mountain sides 2.order accurate in phi-direction
				drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k+1 ] - aux_w.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );

				if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k + 1 ] == 0. ) )				drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k + 1 ] - aux_w.x[ im-1 ][ j ][ k ] ) / dphi;					// 1. order accurate
				if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h.x[ im-1 ][ j ][ k - 1 ] == 1. ) )				drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k ] - aux_w.x[ im-1 ][ j ][ k - 1 ] ) / dphi;					// 1. order accurate

				drhs_vdthe = ( aux_v.x[ im-1 ][ j+1 ][ k ] - aux_v.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );
				drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k+1 ] - aux_w.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );


				p_dynn.x[ im-1 ][ j ][ k ] = ( ( p_dyn.x[ im-1 ][ j+1 ][ k ] + p_dyn.x[ im-1 ][ j-1 ][ k ] ) * num2
													+ ( p_dyn.x[ im-1 ][ j ][ k+1 ] + p_dyn.x[ im-1 ][ j ][ k-1 ] ) * num3 
													+ r_0_water * ( drhs_vdthe + drhs_wdphi ) ) / denom;

				if ( h.x[ im-1 ][ j ][ k ] == 1. )				p_dynn.x[ im-1 ][ j ][ k ] = .0;
			}
		}


// boundary conditions for the the-direction, loop index j
		for ( int k = 0; k < km; k++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			p_dynn.x[ im-1 ][ 0 ][ k ] = 0.;
			p_dynn.x[ im-1 ][ jm-1 ][ k ] = 0.;
		}

// boundary conditions for the phi-direction, loop index k
		for ( int j = 1; j < jm - 1; j++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			p_dynn.x[ im-1 ][ j ][ 0 ] = c43 * p_dynn.x[ im-1 ][ j ][ 1 ] - c13 * p_dynn.x[ im-1 ][ j ][ 2 ];
			p_dynn.x[ im-1 ][ j ][ km-1 ] = c43 * p_dynn.x[ im-1 ][ j ][ km-2 ] - c13 * p_dynn.x[ im-1 ][ j ][ km-3 ];
			p_dynn.x[ im-1 ][ j ][ 0 ] = p_dynn.x[ im-1 ][ j ][ km-1 ] = ( p_dynn.x[ im-1 ][ j ][ 0 ] + p_dynn.x[ im-1 ][ j ][ km-1 ] ) / 2.;
		}
	}
}
