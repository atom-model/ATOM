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

#include "Pressure_Atm.h"
#include "Array.h"
#include "Array_2D.h"
#include "MinMax_Atm.h"

using namespace std;


Pressure_Atm::Pressure_Atm ( int im, int jm, int km, double dr, double dthe, double dphi )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> dr = dr;
	this-> dthe = dthe;
	this-> dphi = dphi;


	c43 = 4./3.;
	c13 = 1./3.;

// array "alfa" for Thomas algorithm
	alfa = 0L;

	alfa = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		alfa[ l ] = 0.;
	}


// array "beta" for Thomas algorithm
	beta = 0L;

	beta = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		beta[ l ] = 0.;
	}
}


Pressure_Atm::~Pressure_Atm ()
{
	delete [  ] alfa;
	delete [  ] beta;
}








void Pressure_Atm::computePressure_3D ( double pa, Array_1D &rad, Array_1D &the, Array &p_dyn, Array &p_dynn, Array &h, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &aux_u, Array &aux_v, Array &aux_w )
{
	cout.precision ( 6 );
	cout.setf ( ios::fixed );

// boundary conditions for the r-direction, loop index i
	for ( int j = 1; j < jm - 1; j++ )
	{
		for ( int k = 1; k < km - 1; k++ )
		{
//			aux_u.x[ 0 ][ j ][ k ] = aux_u.x[ 3 ][ j ][ k ] - 3. * aux_u.x[ 2 ][ j ][ k ] + 3. * aux_u.x[ 1 ][ j ][ k ];										// extrapolation
//			aux_u.x[ im-1 ][ j ][ k ] = aux_u.x[ im-4 ][ j ][ k ] - 3. * aux_u.x[ im-3 ][ j ][ k ] + 3. * aux_u.x[ im-2 ][ j ][ k ];					// extrapolation
//			aux_u.x[ 0 ][ j ][ k ] = c43 * aux_u.x[ 1 ][ j ][ k ] - c13 * aux_u.x[ 2 ][ j ][ k ];									// von Neumann
//			aux_u.x[ im-1 ][ j ][ k ] = c43 * aux_u.x[ im-2 ][ j ][ k ] - c13 * aux_u.x[ im-3 ][ j ][ k ];
			aux_u.x[ 0 ][ j ][ k ] = 0.;
			aux_u.x[ im-1 ][ j ][ k ] = 0.;
		}
	}


// boundary conditions for the the-direction, loop index j
	for ( int k = 0; k < km; k++ )
	{
		for ( int i = 0; i < im; i++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
//			aux_u.x[ i ][ 0 ][ k ] = c43 * aux_u.x[ i ][ 1 ][ k ] - c13 * aux_u.x[ i ][ 2 ][ k ];
//			aux_u.x[ i ][ jm-1 ][ k ] = c43 * aux_u.x[ i ][ jm-2 ][ k ] - c13 * aux_u.x[ i ][ jm-3 ][ k ];
			aux_u.x[ i ][ 0 ][ k ] = 0.;
			aux_u.x[ i ][ jm-1 ][ k ] = 0.;
//			aux_u.x[ i ][ 0 ][ k ] = aux_u.x[ i ][ 3 ][ k ] - 3. * aux_u.x[ i ][ 2 ][ k ] + 3. * aux_u.x[ i ][ 1 ][ k ];												// extrapolation
//			aux_u.x[ i ][ jm - 1 ][ k ] = aux_u.x[ i ][ jm - 4 ][ k ] - 3. * aux_u.x[ i ][ jm - 3 ][ k ] + 3. * aux_u.x[ i ][ jm - 2 ][ k ];					// extrapolation
		}
	}


// boundary conditions for the phi-direction, loop index k
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 1; j < jm - 1; j++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			aux_u.x[ i ][ j ][ 0 ] = c43 * aux_u.x[ i ][ j ][ 1 ] - c13 * aux_u.x[ i ][ j ][ 2 ];
			aux_u.x[ i ][ j ][ km-1 ] = c43 * aux_u.x[ i ][ j ][ km-2 ] - c13 * aux_u.x[ i ][ j ][ km-3 ];
			aux_u.x[ i ][ j ][ 0 ] = aux_u.x[ i ][ j ][ km-1 ] = ( aux_u.x[ i ][ j ][ 0 ] + aux_u.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}



// boundary conditions for the r-direction, loop index i
	for ( int j = 1; j < jm - 1; j++ )
	{
		for ( int k = 1; k < km - 1; k++ )
		{
//			aux_v.x[ 0 ][ j ][ k ] = aux_v.x[ 3 ][ j ][ k ] - 3. * aux_v.x[ 2 ][ j ][ k ] + 3. * aux_v.x[ 1 ][ j ][ k ];											// extrapolation
//			aux_v.x[ im-1 ][ j ][ k ] = aux_v.x[ im-4 ][ j ][ k ] - 3. * aux_v.x[ im-3 ][ j ][ k ] + 3. * aux_v.x[ im-2 ][ j ][ k ];					// extrapolation
			aux_v.x[ 0 ][ j ][ k ] = c43 * aux_v.x[ 1 ][ j ][ k ] - c13 * aux_v.x[ 2 ][ j ][ k ];									// von Neumann
//			aux_v.x[ im-1 ][ j ][ k ] = c43 * aux_v.x[ im-2 ][ j ][ k ] - c13 * aux_v.x[ im-3 ][ j ][ k ];
			aux_v.x[ im-1 ][ j ][ k ] = 0.;
		}
	}


// boundary conditions for the the-direction, loop index j
	for ( int k = 0; k < km; k++ )
	{
		for ( int i = 0; i < im; i++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
//			aux_v.x[ i ][ 0 ][ k ] = c43 * aux_v.x[ i ][ 1 ][ k ] - c13 * aux_v.x[ i ][ 2 ][ k ];
//			aux_v.x[ i ][ jm-1 ][ k ] = c43 * aux_v.x[ i ][ jm-2 ][ k ] - c13 * aux_v.x[ i ][ jm-3 ][ k ];
			aux_v.x[ i ][ 0 ][ k ] = 0.;
			aux_v.x[ i ][ jm-1 ][ k ] = 0.;
//			aux_v.x[ i ][ 0 ][ k ] = aux_v.x[ i ][ 3 ][ k ] - 3. * aux_v.x[ i ][ 2 ][ k ] + 3. * aux_v.x[ i ][ 1 ][ k ];												// extrapolation
//			aux_v.x[ i ][ jm - 1 ][ k ] = aux_v.x[ i ][ jm - 4 ][ k ] - 3. * aux_v.x[ i ][ jm - 3 ][ k ] + 3. * aux_v.x[ i ][ jm - 2 ][ k ];					// extrapolation
//			aux_v.x[ i ][ 0 ][ k ] = 0.;
//			aux_v.x[ i ][ jm-1 ][ k ] = 0.;
		}
	}


// boundary conditions for the phi-direction, loop index k
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 1; j < jm - 1; j++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			aux_v.x[ i ][ j ][ 0 ] = c43 * aux_v.x[ i ][ j ][ 1 ] - c13 * aux_v.x[ i ][ j ][ 2 ];
			aux_v.x[ i ][ j ][ km-1 ] = c43 * aux_v.x[ i ][ j ][ km-2 ] - c13 * aux_v.x[ i ][ j ][ km-3 ];
			aux_v.x[ i ][ j ][ 0 ] = aux_v.x[ i ][ j ][ km-1 ] = ( aux_v.x[ i ][ j ][ 0 ] + aux_v.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}



// boundary conditions for the r-direction, loop index i
	for ( int j = 1; j < jm - 1; j++ )
	{
		for ( int k = 1; k < km - 1; k++ )
		{
//			aux_w.x[ 0 ][ j ][ k ] = aux_w.x[ 3 ][ j ][ k ] - 3. * aux_w.x[ 2 ][ j ][ k ] + 3. * aux_w.x[ 1 ][ j ][ k ];											// extrapolation
//			aux_w.x[ im-1 ][ j ][ k ] = aux_w.x[ im-4 ][ j ][ k ] - 3. * aux_w.x[ im-3 ][ j ][ k ] + 3. * aux_w.x[ im-2 ][ j ][ k ];					// extrapolation
			aux_w.x[ 0 ][ j ][ k ] = c43 * aux_w.x[ 1 ][ j ][ k ] - c13 * aux_w.x[ 2 ][ j ][ k ];									// von Neumann
//			aux_w.x[ im-1 ][ j ][ k ] = c43 * aux_w.x[ im-2 ][ j ][ k ] - c13 * aux_w.x[ im-3 ][ j ][ k ];
			aux_w.x[ im-1 ][ j ][ k ] = 0.;
		}
	}


// boundary conditions for the the-direction, loop index j
	for ( int k = 0; k < km; k++ )
	{
		for ( int i = 0; i < im; i++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
//			aux_w.x[ i ][ 0 ][ k ] = c43 * aux_w.x[ i ][ 1 ][ k ] - c13 * aux_w.x[ i ][ 2 ][ k ];
//			aux_w.x[ i ][ jm-1 ][ k ] = c43 * aux_w.x[ i ][ jm-2 ][ k ] - c13 * aux_w.x[ i ][ jm-3 ][ k ];
			aux_w.x[ i ][ 0 ][ k ] = 0.;
			aux_w.x[ i ][ jm-1 ][ k ] = 0.;
//			aux_w.x[ i ][ 0 ][ k ] = aux_w.x[ i ][ 3 ][ k ] - 3. * aux_w.x[ i ][ 2 ][ k ] + 3. * aux_w.x[ i ][ 1 ][ k ];												// extrapolation
//			aux_w.x[ i ][ jm - 1 ][ k ] = aux_w.x[ i ][ jm - 4 ][ k ] - 3. * aux_w.x[ i ][ jm - 3 ][ k ] + 3. * aux_w.x[ i ][ jm - 2 ][ k ];					// extrapolation
//			aux_w.x[ i ][ 0 ][ k ] = 0.;
//			aux_w.x[ i ][ jm-1 ][ k ] = 0.;
		}
	}


// boundary conditions for the phi-direction, loop index k
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 1; j < jm - 1; j++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			aux_w.x[ i ][ j ][ 0 ] = c43 * aux_w.x[ i ][ j ][ 1 ] - c13 * aux_w.x[ i ][ j ][ 2 ];
			aux_w.x[ i ][ j ][ km-1 ] = c43 * aux_w.x[ i ][ j ][ km-2 ] - c13 * aux_w.x[ i ][ j ][ km-3 ];
			aux_w.x[ i ][ j ][ 0 ] = aux_w.x[ i ][ j ][ km-1 ] = ( aux_w.x[ i ][ j ][ 0 ] + aux_w.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}



	switch_pres = 0;

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
				rm2sinthe = rm2 * sinthe;
				rm2sinthe2 = rm2 * sinthe2;
				rm2dthe2 = rm2 * dthe2;
				denom = 2. / dr2 + 2. / ( rm2 * dthe2 ) + 2. / ( rm2sinthe2 * dphi2 );
				num1 = 1. / dr2;
				num2 = 1. / ( rm2 * dthe2 );
				num3 = 1. / ( rm2sinthe2 * dphi2 );

// gradients of RHS terms at mountain summits 2.order accurate in r-direction
				for ( int k = 1; k < km-1; k++ )
				{
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i + 1][ j ][ k ] == 0. ) )			aux_u.x[ i ][ j ][ k ] = c43 * aux_u.x[ i + 1 ][ j ][ k ] - c13 * aux_u.x[ i + 2 ][ j ][ k ];									// von Neumann
/*
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j + 1 ][ k ] == 0. ) )			aux_v.x[ i ][ j ][ k ] = c43 * aux_v.x[ i ][ j + 1 ][ k ] - c13 * aux_v.x[ i ][ j + 2 ][ k ];
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) )				aux_v.x[ i ][ j ][ k ] = c43 * aux_v.x[ i ][ j - 1 ][ k ] - c13 * aux_v.x[ i ][ j - 2 ][ k ];

					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) )			aux_w.x[ i ][ j ][ k ] = c43 * aux_w.x[ i ][ j ][ k + 1 ] - c13 * aux_w.x[ i ][ j ][ k + 2 ];
					if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h.x[ i ][ j ][ k - 1 ] == 1. ) )				aux_w.x[ i ][ j ][ k ] = c43 * aux_w.x[ i ][ j ][ k - 1 ] - c13 * aux_w.x[ i ][ j ][ k - 2 ];
*/
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j + 1 ][ k ] == 0. ) )			aux_v.x[ i ][ j ][ k ] = aux_v.x[ i ][ j + 1 ][ k ];
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) )				aux_v.x[ i ][ j ][ k ] = aux_v.x[ i ][ j - 1 ][ k ];

					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) )			aux_w.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k + 1 ];
					if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h.x[ i ][ j ][ k - 1 ] == 1. ) )				aux_w.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k - 1 ];


					drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
//					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i + 1 ][ j ][ k ] == 0. ) )			drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i ][ j ][ k ] ) / dr;					// 1. order accurate
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i + 1 ][ j ][ k ] == 0. ) )			drhs_udr = ( - 3. * aux_u.x[ i ][ j ][ k ] + 4. * aux_u.x[ i + 1 ][ j ][ k ] - aux_u.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );					// 2. order accurate



// gradients of RHS terms at mountain sides 2.order accurate in the-direction
					drhs_vdthe = ( aux_v.x[ i ][ j+1 ][ k ] - aux_v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe * rm );

					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j + 1 ][ k ] == 0. ) )			drhs_vdthe = ( aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j ][ k ] ) / ( rm * dthe );					// 1. order accurate
					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) )				drhs_vdthe = ( aux_v.x[ i ][ j ][ k ] - aux_v.x[ i ][ j - 1 ][ k ] ) / ( rm * dthe );					// 1. order accurate

/*
					if ( ( j >= 2 ) && ( j < jm - 2 ) )
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( ( h.x[ i ][ j + 1 ][ k ] == 0. ) && ( h.x[ i ][ j + 2 ][ k ] == 0. ) ) )
						{
							drhs_vdthe = ( - 3. * aux_v.x[ i ][ j ][ k ] + 4. * aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j + 2 ][ k ] ) / ( 2. * rm * dthe );					// 2. order accurate
						}
						else 			drhs_vdthe = ( aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j ][ k ] ) / ( rm * dthe );

						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) && ( h.x[ i ][ j - 2 ][ k ] == 0. ) )
						{
							drhs_vdthe = ( - 3. * aux_v.x[ i ][ j ][ k ] + 4. * aux_v.x[ i ][ j - 1 ][ k ] - aux_v.x[ i ][ j - 2 ][ k ] ) / ( 2. * rm * dthe );					// 2. order accurate
						}
						else			drhs_vdthe = ( aux_v.x[ i ][ j ][ k ] - aux_v.x[ i ][ j - 1 ][ k ] ) / ( rm * dthe );
					}
					else
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j + 1 ][ k ] == 0. ) )			drhs_vdthe = ( aux_v.x[ i ][ j + 1 ][ k ] - aux_v.x[ i ][ j ][ k ] ) / ( rm * dthe );
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) )				drhs_vdthe = ( aux_v.x[ i ][ j ][ k ] - aux_v.x[ i ][ j - 1 ][ k ] ) / ( rm * dthe );
					}
*/


// gradients of RHS terms at mountain sides 2.order accurate in phi-direction
					drhs_wdphi = ( aux_w.x[ i ][ j ][ k+1 ] - aux_w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe );

					if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) )			drhs_wdphi = ( aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k ] ) / ( rmsinthe * dphi );					// 1. order accurate
					if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h.x[ i ][ j ][ k - 1 ] == 1. ) )				drhs_wdphi = ( aux_w.x[ i ][ j ][ k ] - aux_w.x[ i ][ j ][ k - 1 ] ) / ( rmsinthe * dphi );					// 1. order accurate

/*
					if ( ( k >= 2 ) && ( k < km - 2 ) )
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) && ( h.x[ i ][ j ][ k + 2 ] == 0. ) )
						{
							drhs_wdphi = ( - 3. * aux_w.x[ i ][ j ][ k ] + 4. * aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k + 2 ] ) / ( 2. * rmsinthe * dphi );					// 2. order accurate
						}
						else 			drhs_wdphi = ( aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k ] ) / ( rmsinthe * dphi );

						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k - 1 ] == 0. ) && ( h.x[ i ][ j ][ k - 2 ] == 0. ) )
						{
							drhs_wdphi = ( - 3. * aux_w.x[ i ][ j ][ k ] + 4. * aux_w.x[ i ][ j ][ k - 1 ] - aux_w.x[ i ][ j ][ k - 2 ] ) / ( 2. * rmsinthe * dphi );					// 2. order accurate
						}
						else			drhs_wdphi = ( aux_w.x[ i ][ j ][ k ] - aux_w.x[ i ][ j ][ k - 1 ] ) / ( rmsinthe * dphi );
					}
					else
					{
						if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) )			drhs_wdphi = ( aux_w.x[ i ][ j ][ k + 1 ] - aux_w.x[ i ][ j ][ k ] ) / ( rmsinthe * dphi );
						if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h.x[ i ][ j ][ k - 1 ] == 1. ) )				drhs_wdphi = ( aux_w.x[ i ][ j ][ k ] - aux_w.x[ i ][ j ][ k - 1 ] ) / ( rmsinthe * dphi );
					}
*/


// explicit pressure computation based on the Poisson equation
					if ( switch_pres == 1 )
					{
						p_dynn.x[ i ][ j ][ k ] = ( ( p_dyn.x[ i+1 ][ j ][ k ] + p_dyn.x[ i-1 ][ j ][ k ] ) * num1
															+ ( p_dyn.x[ i ][ j+1 ][ k ] + p_dyn.x[ i ][ j-1 ][ k ] ) * num2
															+ ( p_dyn.x[ i ][ j ][ k+1 ] + p_dyn.x[ i ][ j ][ k-1 ] ) * num3 
															- ( drhs_udr + drhs_vdthe + drhs_wdphi ) ) / denom;
					}
					else 				p_dynn.x[ i ][ j ][ k ] = drhs_udr + drhs_vdthe + drhs_wdphi;
					if ( h.x[ i ][ j ][ k ] == 1. )				p_dynn.x[ i ][ j ][ k ] = .0;
				}
			}
		}



// boundary conditions for the r-direction, loop index i
	for ( int j = 1; j < jm - 1; j++ )
	{
		for ( int k = 1; k < km - 1; k++ )
		{
//			p_dynn.x[ 0 ][ j ][ k ] = 0.;
//			p_dynn.x[ 0 ][ j ][ k ] = c43 * p_dynn.x[ 1 ][ j ][ k ] - c13 * p_dynn.x[ 2 ][ j ][ k ];
			p_dynn.x[ im-1 ][ j ][ k ] = 0.;
		}
	}


// boundary conditions for the the-direction, loop index j
	for ( int k = 0; k < km; k++ )
	{
		for ( int i = 0; i < im; i++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			p_dynn.x[ i ][ 0 ][ k ] = 0.;
			p_dynn.x[ i ][ jm-1 ][ k ] = 0.;
//			p_dynn.x[ i ][ 0 ][ k ] = c43 * p_dynn.x[ i ][ 1 ][ k ] - c13 * p_dynn.x[ i ][ 2 ][ k ];
//			p_dynn.x[ i ][ jm-1 ][ k ] = c43 * p_dynn.x[ i ][ jm-2 ][ k ] - c13 * p_dynn.x[ i ][ jm-3 ][ k ];
		}
	}


// boundary conditions for the phi-direction, loop index k
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 1; j < jm - 1; j++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			p_dynn.x[ i ][ j ][ 0 ] = c43 * p_dynn.x[ i ][ j ][ 1 ] - c13 * p_dynn.x[ i ][ j ][ 2 ];
			p_dynn.x[ i ][ j ][ km-1 ] = c43 * p_dynn.x[ i ][ j ][ km-2 ] - c13 * p_dynn.x[ i ][ j ][ km-3 ];
			p_dynn.x[ i ][ j ][ 0 ] = p_dynn.x[ i ][ j ][ km-1 ] = ( p_dynn.x[ i ][ j ][ 0 ] + p_dynn.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}																								// end of while


		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				for ( int i = 0; i < im; i++ )
				{
//					if ( ( j == 62 ) && ( k == 87 ) )	cout << "  iter_prec = " << iter_prec << " recurrence:  i = "<< i << "  p_dyn = " << p_dyn.x[ i ][ j ][ k ] << "  p_dynn = " << p_dynn.x[ i ][ j ][ k ] << endl;

//					if ( switch_pres == 1 )					p_dyn.x[ i ][ j ][ k ] = .5 * ( p_dyn.x[ i ][ j ][ k ] + p_dynn.x[ i ][ j ][ k ] );
//					if ( switch_pres == 1 )					p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ];
					if ( switch_pres == 1 )					p_dyn.x[ i ][ j ][ k ] = 2. / 3. *p_dyn.x[ i ][ j ][ k ] + 1. / 3. *p_dynn.x[ i ][ j ][ k ];
//					if ( switch_pres == 1 )					p_dyn.x[ i ][ j ][ k ] = 1. / 4. *p_dyn.x[ i ][ j ][ k ] + 3. / 4. *p_dynn.x[ i ][ j ][ k ];
					else 												p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ];

					if ( h.x[ i ][ j ][ k ] == 1. )					p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ] = .0;
//					if ( p_dynn.x[ i ][ j ][ k ] >= .05 )		p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ] = .05;
//					if ( p_dynn.x[ i ][ j ][ k ] <= - .05 )		p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ] = - .05;

//					if ( ( j == 62 ) && ( k == 87 ) )	cout << "  iter_prec = " << iter_prec << " recurrence:  i = "<< i << "  p_dyn = " << p_dyn.x[ i ][ j ][ k ] << "  p_dynn = " << p_dynn.x[ i ][ j ][ k ] << endl;

				}
			}
		}
		switch_pres = 1;
	}
}








void Pressure_Atm::computePressure_2D ( double pa, Array_1D &rad, Array_1D &the, Array &p_dyn, Array &p_dynn, Array &h, Array &rhs_v, Array &rhs_w, Array &aux_v, Array &aux_w )
{
// Pressure using Euler equation ( 2. derivative of pressure added to the Poisson-right-hand-side )

// boundary conditions for the the-direction, loop index j
	for ( int k = 0; k < km; k++ )
	{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
		aux_v.x[ 0 ][ 0 ][ k ] = c43 * aux_v.x[ 0 ][ 1 ][ k ] - c13 * aux_v.x[ 0 ][ 2 ][ k ];
		aux_v.x[ 0 ][ jm-1 ][ k ] = c43 * aux_v.x[ 0 ][ jm-2 ][ k ] - c13 * aux_v.x[ 0 ][ jm-3 ][ k ];
//		aux_v.x[ 0 ][ 0 ][ k ] = 0.;
//		aux_v.x[ 0 ][ jm-1 ][ k ] = 0.;
	}


// boundary conditions for the phi-direction, loop index k
	for ( int j = 0; j < jm; j++ )
	{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
		aux_v.x[ 0 ][ j ][ 0 ] = c43 * aux_v.x[ 0 ][ j ][ 1 ] - c13 * aux_v.x[ 0 ][ j ][ 2 ];
		aux_v.x[ 0 ][ j ][ km-1 ] = c43 * aux_v.x[ 0 ][ j ][ km-2 ] - c13 * aux_v.x[ 0 ][ j ][ km-3 ];
		aux_v.x[ 0 ][ j ][ 0 ] = aux_v.x[ 0 ][ j ][ km-1 ] = ( aux_v.x[ 0 ][ j ][ 0 ] + aux_v.x[ 0 ][ j ][ km-1 ] ) / 2.;
	}



// boundary conditions for the the-direction, loop index j
	for ( int k = 0; k < km; k++ )
	{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
		aux_w.x[ 0 ][ 0 ][ k ] = c43 * aux_w.x[ 0 ][ 1 ][ k ] - c13 * aux_w.x[ 0 ][ 2 ][ k ];
		aux_w.x[ 0 ][ jm-1 ][ k ] = c43 * aux_w.x[ 0 ][ jm-2 ][ k ] - c13 * aux_w.x[ 0 ][ jm-3 ][ k ];
//		aux_w.x[ 0 ][ 0 ][ k ] = 0.;
//		aux_w.x[ 0 ][ jm-1 ][ k ] = 0.;
	}


// boundary conditions for the phi-direction, loop index k
	for ( int j = 0; j < jm; j++ )
	{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
		aux_w.x[ 0 ][ j ][ 0 ] = c43 * aux_w.x[ 0 ][ j ][ 1 ] - c13 * aux_w.x[ 0 ][ j ][ 2 ];
		aux_w.x[ 0 ][ j ][ km-1 ] = c43 * aux_w.x[ 0 ][ j ][ km-2 ] - c13 * aux_w.x[ 0 ][ j ][ km-3 ];
		aux_w.x[ 0 ][ j ][ 0 ] = aux_w.x[ 0 ][ j ][ km-1 ] = ( aux_w.x[ 0 ][ j ][ 0 ] + aux_w.x[ 0 ][ j ][ km-1 ] ) / 2.;
	}


	switch_pres = 0;


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
			sinthe2 = sinthe * sinthe;
			costhe = cos( the.z[ j ] );
			cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
			rmsinthe = rm * sinthe;
			rm2sinthe = rm2 * sinthe;
			rm2sinthe2 = rm2 * sinthe2;
			rm2dthe2 = rm2 * dthe2;
			denom = 2. / ( rm2 * dthe2 ) + 2. / ( rm2sinthe2 * dphi2 );
			num2 = 1. / ( rm2 * dthe2 );
			num3 = 1. / ( rm2sinthe2 * dphi2 );

			for ( int k = 1; k < km-1; k++ )
			{
/*
				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j + 1 ][ k ] == 0. ) )			aux_v.x[ 0 ][ j ][ k ] = c43 * aux_v.x[ 0 ][ j + 1 ][ k ] - c13 * aux_v.x[ 0 ][ j + 2 ][ k ];
				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j - 1 ][ k ] == 0. ) )			aux_v.x[ 0 ][ j ][ k ] = c43 * aux_v.x[ 0 ][ j - 1 ][ k ] - c13 * aux_v.x[ 0 ][ j - 2 ][ k ];

				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j ][ k + 1 ] == 0. ) )			aux_w.x[ 0 ][ j ][ k ] = c43 * aux_w.x[ 0 ][ j ][ k + 1 ] - c13 * aux_w.x[ 0 ][ j ][ k + 2 ];
				if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( h.x[ 0 ][ j ][ k - 1 ] == 1. ) )			aux_w.x[ 0 ][ j ][ k ] = c43 * aux_w.x[ 0 ][ j ][ k - 1 ] - c13 * aux_w.x[ 0 ][ j ][ k - 2 ];
*/
				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j + 1 ][ k ] == 0. ) )			aux_v.x[ 0 ][ j ][ k ] = aux_v.x[ 0 ][ j + 1 ][ k ];
				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j - 1 ][ k ] == 0. ) )			aux_v.x[ 0 ][ j ][ k ] = aux_v.x[ 0 ][ j - 1 ][ k ];

				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j ][ k + 1 ] == 0. ) )			aux_w.x[ 0 ][ j ][ k ] = aux_w.x[ 0 ][ j ][ k + 1 ];
				if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( h.x[ 0 ][ j ][ k - 1 ] == 1. ) )			aux_w.x[ 0 ][ j ][ k ] = aux_w.x[ 0 ][ j ][ k - 1 ];



// gradients of RHS terms at mountain sides 2.order accurate in the-direction
				drhs_vdthe = ( aux_v.x[ 0 ][ j+1 ][ k ] - aux_v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe * rm );

				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j + 1 ][ k ] == 0. ) )				drhs_vdthe = ( aux_v.x[ 0 ][ j + 1 ][ k ] - aux_v.x[ 0 ][ j ][ k ] ) / ( rm * dthe );					// 1. order accurate
				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j - 1 ][ k ] == 0. ) )				drhs_vdthe = ( aux_v.x[ 0 ][ j ][ k ] - aux_v.x[ 0 ][ j - 1 ][ k ] ) / ( rm * dthe );					// 1. order accurate


// gradients of RHS terms at mountain sides 2.order accurate in phi-direction
				drhs_wdphi = ( aux_w.x[ 0 ][ j ][ k+1 ] - aux_w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe );

				if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 0 ][ j ][ k + 1 ] == 0. ) )				drhs_wdphi = ( aux_w.x[ 0 ][ j ][ k + 1 ] - aux_w.x[ 0 ][ j ][ k ] ) / ( rmsinthe * dphi );					// 1. order accurate
				if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( h.x[ 0 ][ j ][ k - 1 ] == 1. ) )				drhs_wdphi = ( aux_w.x[ 0 ][ j ][ k ] - aux_w.x[ 0 ][ j ][ k - 1 ] ) / ( rmsinthe * dphi );					// 1. order accurate




				drhs_vdthe = ( aux_v.x[ 0 ][ j+1 ][ k ] - aux_v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe * rm );
				drhs_wdphi = ( aux_w.x[ 0 ][ j ][ k+1 ] - aux_w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe );

				if ( switch_pres == 1 )
				{
					p_dynn.x[ 0 ][ j ][ k ] = ( ( p_dyn.x[ 0 ][ j+1 ][ k ] + p_dyn.x[ 0 ][ j-1 ][ k ] ) * num2
														+ ( p_dyn.x[ 0 ][ j ][ k+1 ] + p_dyn.x[ 0 ][ j ][ k-1 ] ) * num3
														- ( drhs_vdthe + drhs_wdphi ) ) / denom;
				}
				else									p_dynn.x[ 0 ][ j ][ k ] = drhs_vdthe + drhs_wdphi;

				if ( h.x[ 0 ][ j ][ k ] == 1. )				p_dynn.x[ 0 ][ j ][ k ] = .0;


//	if ( ( j == 90 ) && ( k == 180 ) )	cout << "  RHS:  " << "  h = " << h.x[ 0 ][ j ][ k ] << "  aux_v = " << aux_v.x[ 0 ][ j ][ k ] << "  aux_w = " << aux_w.x[ 0 ][ j ][ k ] << "  drhs_vdthe = " << drhs_vdthe << "  drhs_wdphi = " << drhs_wdphi << "  d2pj = " << ( p_dyn.x[ 0 ][ j+1 ][ k ] - 2. * p_dyn.x[ 0 ][ j ][ k ] + p_dyn.x[ 0 ][ j-1 ][ k ] ) * num2 << "  d2pk = " << ( p_dyn.x[ 0 ][ j ][ k+1 ] - 2. * p_dyn.x[ 0 ][ j ][ k ] + p_dyn.x[ 0 ][ j ][ k-1 ] ) * num3 << " num2 = "<< num2 << " num3 = "<< num3 << endl;

			}
		}



// boundary conditions for the the-direction, loop index j
		for ( int k = 0; k < km; k++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
//			p_dynn.x[ 0 ][ 0 ][ k ] = c43 * p_dynn.x[ 0 ][ 1 ][ k ] - c13 * p_dynn.x[ 0 ][ 2 ][ k ];
//			p_dynn.x[ 0 ][ jm-1 ][ k ] = c43 * p_dynn.x[ 0 ][ jm-2 ][ k ] - c13 * p_dynn.x[ 0 ][ jm-3 ][ k ];
			p_dynn.x[ 0 ][ 0 ][ k ] = 0.;
			p_dynn.x[ 0 ][ jm-1 ][ k ] = 0.;
//			p_dynn.x[ 0 ][ 0 ][ k ] = p_dynn.x[ 0 ][ 3 ][ k ] - 3. * p_dynn.x[ 0 ][ 2 ][ k ] + 3. * p_dynn.x[ 0 ][ 1 ][ k ];												// extrapolation
//			p_dynn.x[ 0 ][ jm - 1 ][ k ] = p_dynn.x[ 0 ][ jm - 4 ][ k ] - 3. * p_dynn.x[ 0 ][ jm - 3 ][ k ] + 3. * p_dynn.x[ 0 ][ jm - 2 ][ k ];					// extrapolation
		}


// boundary conditions for the phi-direction, loop index k
		for ( int j = 1; j < jm - 1; j++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			p_dynn.x[ 0 ][ j ][ 0 ] = c43 * p_dynn.x[ 0 ][ j ][ 1 ] - c13 * p_dynn.x[ 0 ][ j ][ 2 ];
			p_dynn.x[ 0 ][ j ][ km-1 ] = c43 * p_dynn.x[ 0 ][ j ][ km-2 ] - c13 * p_dynn.x[ 0 ][ j ][ km-3 ];
			p_dynn.x[ 0 ][ j ][ 0 ] = p_dynn.x[ 0 ][ j ][ km-1 ] = ( p_dynn.x[ 0 ][ j ][ 0 ] + p_dynn.x[ 0 ][ j ][ km-1 ] ) / 2.;
		}



		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				if ( ( j == 90 ) && ( k == 180 ) )	cout << "  iter_prec = " << iter_prec << " recurrence: " << "  p_dyn = " << p_dyn.x[ 0 ][ j ][ k ] << "  p_dynn = " << p_dynn.x[ 0 ][ j ][ k ] << endl;

//				if ( switch_pres == 1 )						p_dyn.x[ 0 ][ j ][ k ] = .5 * ( p_dyn.x[ 0 ][ j ][ k ] + p_dynn.x[ 0 ][ j ][ k ] );
				if ( switch_pres == 1 )						p_dyn.x[ 0 ][ j ][ k ] = 2. / 3. *p_dyn.x[ 0 ][ j ][ k ] + 1. / 3. *p_dynn.x[ 0 ][ j ][ k ];
//				if ( switch_pres == 1 )						p_dyn.x[ 0 ][ j ][ k ] = 1. / 4. *p_dyn.x[ 0 ][ j ][ k ] + 3. / 4. *p_dynn.x[ 0 ][ j ][ k ];
//				if ( switch_pres == 1 )						p_dyn.x[ 0 ][ j ][ k ] = p_dynn.x[ 0 ][ j ][ k ];
				else 													p_dyn.x[ 0 ][ j ][ k ] = p_dynn.x[ 0 ][ j ][ k ];

				if ( h.x[ 0 ][ j ][ k ] == 1. )					p_dyn.x[ 0 ][ j ][ k ] = p_dynn.x[ 0 ][ j ][ k ] = .0;

//				if ( p_dynn.x[ 0 ][ j ][ k ] >= .05 )		p_dyn.x[ 0 ][ j ][ k ] = p_dynn.x[ 0 ][ j ][ k ] = .05;
//				if ( p_dynn.x[ 0 ][ j ][ k ] <= - .05 )	p_dyn.x[ 0 ][ j ][ k ] = p_dynn.x[ 0 ][ j ][ k ] = - .05;

			}
		}
		switch_pres = 1;
	}
}
