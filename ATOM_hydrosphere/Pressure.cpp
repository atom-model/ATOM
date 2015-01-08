/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to compute the pressure independent of the other variables
*/


#include <iostream>
#include <cmath>

#include "Pressure.h"

using namespace std;




Pressure::Pressure ( int im, int jm, int km, double dr, double dthe, double dphi )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> dr = dr;
	this-> dthe = dthe;
	this-> dphi = dphi;

	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;
}


Pressure::~Pressure () {}




void Pressure::computePressure ( double pa, Array_1D &rad, Array_1D &the, Array &p, Array &h, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &aux_u, Array &aux_v, Array &aux_w )
{
// Pressure using Euler equation ( 2. derivative of pressure added to the Poisson-right-hand-side )

	for ( int i = 1; i < im-1; i++ )
	{
		dr2 = dr * dr;
		dthe2 = dthe * dthe;
		dphi2 = dphi * dphi;
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
			denom = 2. / dr2 + 2. / ( rm2 * dthe2 ) + 2. / ( rmsinthe * dphi2 );
			num1 = 1. / dr2;
			num2 = 1. / ( rm2 * dthe2 );
			num3 = 1. / ( rmsinthe * dphi2 );

			for ( int k = 1; k < km-1; k++ )
			{
				drhs_udr = ( aux_u.x[ i+1 ][ j ][ k ] - aux_u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
				drhs_vdthe = ( aux_v.x[ i ][ j+1 ][ k ] - aux_v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe * rm );
				drhs_wdphi = ( aux_w.x[ i ][ j ][ k+1 ] - aux_w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe );

				if ( h.x[ i ][ j ][ k ] == 0. )
					{
						p.x[ i ][ j ][ k ] = ( ( p.x[ i+1 ][ j ][ k ] + p.x[ i-1 ][ j ][ k ] ) / num1 +
													( p.x[ i ][ j+1 ][ k ] + p.x[ i ][ j-1 ][ k ] ) / num2 +
													( p.x[ i ][ j ][ k+1 ] + p.x[ i ][ j ][ k-1 ] ) / num3 + 
													drhs_udr + drhs_vdthe + drhs_wdphi ) / denom;
					}
				else    p.x[ i ][ j ][ k ] = pa;
			}
		}
	}
}





void Pressure::computePressure_2D ( double pa, Array_1D &rad, Array_1D &the, Array &p, Array &h, Array &rhs_v, Array &rhs_w, Array &aux_v, Array &aux_w )
{
// Pressure using Euler equation ( 2. derivative of pressure added to the Poisson-right-hand-side )
	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

	rm = rad.z[ im-1 ];
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
		denom = 2. / ( rm2 * dthe2 ) + 2. / ( rmsinthe * dphi2 );
		num2 = 1. / rm2dthe2;
		num3 = 1. / ( rmsinthe * dphi2 );

		for ( int k = 1; k < km-1; k++ )
		{
			denom = 2. / ( rm2 * dthe2 ) + 2. / ( rmsinthe * dphi2 );
			drhs_vdthe = ( aux_v.x[ im-1 ][ j+1 ][ k ] - aux_v.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe * rm );
			drhs_wdphi = ( aux_w.x[ im-1 ][ j ][ k+1 ] - aux_w.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe );

			if ( h.x[ im-1 ][ j ][ k ] == 0. )
			{
				p.x[ im-1 ][ j ][ k ] = ( ( p.x[ im-1 ][ j+1 ][ k ] + p.x[ im-1 ][ j-1 ][ k ] ) * num2
										+ ( p.x[ im-1 ][ j ][ k+1 ] + p.x[ im-1 ][ j ][ k-1 ] ) * num3
										- ( drhs_vdthe + drhs_wdphi ) ) / denom;
			}
			else    p.x[ im-1 ][ j ][ k ] = pa;
		}
	}

}
