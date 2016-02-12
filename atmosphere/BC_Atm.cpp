/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the coordinate system
*/


#include <iostream>
#include <cmath>

#include "BC_Atm.h"

using namespace std;




BC_Atmosphere::BC_Atmosphere ( int im, int jm, int km )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;

	c43 = 4./3.;
	c13 = 1./3.;
}


BC_Atmosphere::~BC_Atmosphere() {}




void BC_Atmosphere::BC_radius ( double t_tropopause, double c_tropopause, double co2_tropopause, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &co2 )
{
// boundary conditions for the r-direction, loop index i
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			u.x[ 0 ][ j ][ k ] = 0.;
//			v.x[ 0 ][ j ][ k ] = 2. * dr / rad.z[ 1 ] * v.x[ 1 ][ j ][ k ] + v.x[ 2 ][ j ][ k ];				// sea surface
//			w.x[ 0 ][ j ][ k ] = 2. * dr / rad.z[ 1 ] * w.x[ 1 ][ j ][ k ] + w.x[ 2 ][ j ][ k ];			// sea surface
//			v.x[ 0 ][ j ][ k ] = c43 * v.x[ 1 ][ j ][ k ] - c13 * v.x[ 2 ][ j ][ k ];							// zero tangent
//			w.x[ 0 ][ j ][ k ] = c43 * w.x[ 1 ][ j ][ k ] - c13 * w.x[ 2 ][ j ][ k ];						// zero tangent
			v.x[ 0 ][ j ][ k ] = v.x[ 3 ][ j ][ k ] - 3. * v.x[ 2 ][ j ][ k ] + 3. * v.x[ 1 ][ j ][ k ];		// extrapolation
			w.x[ 0 ][ j ][ k ] = w.x[ 3 ][ j ][ k ] - 3. * w.x[ 2 ][ j ][ k ] + 3. * w.x[ 1 ][ j ][ k ];	// leads to round harmonic profiles

//			p_dyn.x[ 0 ][ j ][ k ] = c43 * p_dyn.x[ 1 ][ j ][ k ] - c13 * p_dyn.x[ 2 ][ j ][ k ];						// zero tangent
			p_dyn.x[ 0 ][ j ][ k ] = p_dyn.x[ 3 ][ j ][ k ] - 3. * p_dyn.x[ 2 ][ j ][ k ] + 3. * p_dyn.x[ 1 ][ j ][ k ];	// leads to round harmonic profiles

			u.x[ im-1 ][ j ][ k ] = 0.;
			v.x[ im-1 ][ j ][ k ] = 0.;																					// stratosphere
			w.x[ im-1 ][ j ][ k ] = 0.;

			t.x[ im-1 ][ j ][ k ] = t_tropopause;
			c.x[ im-1 ][ j ][ k ] = c_tropopause;
			co2.x[ im-1 ][ j ][ k ] = co2_tropopause;
			p_dyn.x[ im-1 ][ j ][ k ] = p_dyn.x[ im-4 ][ j ][ k ] - 3. * p_dyn.x[ im-3 ][ j ][ k ] + 3. * p_dyn.x[ im-2 ][ j ][ k ];

		}
	}
}






void BC_Atmosphere::BC_theta ( double t_tropopause, double c_tropopause, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &co2 )
{
// boundary conditions for the the-direction, loop index j

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )

			t.x[ i ][ 0 ][ k ] = c43 * t.x[ i ][ 1 ][ k ] - c13 * t.x[ i ][ 2 ][ k ];
			t.x[ i ][ jm-1 ][ k ] = c43 * t.x[ i ][ jm-2 ][ k ] - c13 * t.x[ i ][ jm-3 ][ k ];

			u.x[ i ][ 0 ][ k ] = c43 * u.x[ i ][ 1 ][ k ] - c13 * u.x[ i ][ 2 ][ k ];
			u.x[ i ][ jm-1 ][ k ] = c43 * u.x[ i ][ jm-2 ][ k ] - c13 * u.x[ i ][ jm-3 ][ k ];

//			v.x[ i ][ 0 ][ k ] = c43 * v.x[ i ][ 1 ][ k ] - c13 * v.x[ i ][ 2 ][ k ];
			v.x[ i ][ 0 ][ k ] = 0.;
//			v.x[ i ][ jm-1 ][ k ] = c43 * v.x[ i ][ jm-2 ][ k ] - c13 * v.x[ i ][ jm-3 ][ k ];
			v.x[ i ][ jm-1 ][ k ] = 0.;

//			w.x[ i ][ 0 ][ k ] = c43 * w.x[ i ][ 1 ][ k ] - c13 * w.x[ i ][ 2 ][ k ];
			w.x[ i ][ 0 ][ k ] = 0.;
//			w.x[ i ][ jm-1 ][ k ] = c43 * w.x[ i ][ jm-2 ][ k ] - c13 * w.x[ i ][ jm-3 ][ k ];
			w.x[ i ][ jm-1 ][ k ] = 0.;

			p_dyn.x[ i ][ 0 ][ k ] = c43 * p_dyn.x[ i ][ 1 ][ k ] - c13 * p_dyn.x[ i ][ 2 ][ k ];
			p_dyn.x[ i ][ jm-1 ][ k ] = c43 * p_dyn.x[ i ][ jm-2 ][ k ] - c13 * p_dyn.x[ i ][ jm-3 ][ k ];

			c.x[ i ][ 0 ][ k ] = c43 * c.x[ i ][ 1 ][ k ] - c13 * c.x[ i ][ 2 ][ k ];
			c.x[ i ][ jm-1 ][ k ] = c43 * c.x[ i ][ jm-2 ][ k ] - c13 * c.x[ i ][ jm-3 ][ k ];

			co2.x[ i ][ 0 ][ k ] = c43 * co2.x[ i ][ 1 ][ k ] - c13 * co2.x[ i ][ 2 ][ k ];
			co2.x[ i ][ jm-1 ][ k ] = c43 * co2.x[ i ][ jm-2 ][ k ] - c13 * co2.x[ i ][ jm-3 ][ k ];

		}
	}
}







void BC_Atmosphere::BC_phi ( Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &co2 )
{
// boundary conditions for the phi-direction, loop index k

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )

//			t.x[ i ][ j ][ 0 ] = c43 * t.x[ i ][ j ][ 1 ] - c13 * t.x[ i ][ j ][ 2 ];
//			t.x[ i ][ j ][ km-1 ] = c43 * t.x[ i ][ j ][ km-2 ] - c13 * t.x[ i ][ j ][ km-3 ];

			t.x[ i ][ j ][ 0 ] = t.x[ i ][ j ][ 3 ] - 3. * t.x[ i ][ j ][ 2 ] + 3. * t.x[ i ][ j ][ 1 ];		// extrapolation
			t.x[ i ][ j ][ km-1 ] = t.x[ i ][ j ][ km-4 ] - 3. * t.x[ i ][ j ][ km-3 ] + 3. * t.x[ i ][ j ][ km-2 ];		// extrapolation
//			t.x[ i ][ j ][ 0 ] = t.x[ i ][ j ][ km-1 ] = ( t.x[ i ][ j ][ 0 ] + t.x[ i ][ j ][ km-1 ] ) / 2.;

//			u.x[ i ][ j ][ 0 ] = c43 * u.x[ i ][ j ][ 1 ] - c13 * u.x[ i ][ j ][ 2 ];
//			u.x[ i ][ j ][ km-1 ] = c43 * u.x[ i ][ j ][ km-2 ] - c13 * u.x[ i ][ j ][ km-3 ];

			u.x[ i ][ j ][ 0 ] = u.x[ i ][ j ][ 3 ] - 3. * u.x[ i ][ j ][ 2 ] + 3. * u.x[ i ][ j ][ 1 ];		// extrapolation
			u.x[ i ][ j ][ km-1 ] = u.x[ i ][ j ][ km-4 ] - 3. * u.x[ i ][ j ][ km-3 ] + 3. * u.x[ i ][ j ][ km-2 ];		// extrapolation
//			u.x[ i ][ j ][ 0 ] = u.x[ i ][ j ][ km-1 ] = ( u.x[ i ][ j ][ 0 ] + u.x[ i ][ j ][ km-1 ] ) / 2.;

//			v.x[ i ][ j ][ 0 ] = c43 * v.x[ i ][ j ][ 1 ] - c13 * v.x[ i ][ j ][ 2 ];
//			v.x[ i ][ j ][ km-1 ] = c43 * v.x[ i ][ j ][ km-2 ] - c13 * v.x[ i ][ j ][ km-3 ];

			v.x[ i ][ j ][ 0 ] = v.x[ i ][ j ][ 3 ] - 3. * v.x[ i ][ j ][ 2 ] + 3. * v.x[ i ][ j ][ 1 ];		// extrapolation
			v.x[ i ][ j ][ km-1 ] = v.x[ i ][ j ][ km-4 ] - 3. * v.x[ i ][ j ][ km-3 ] + 3. * v.x[ i ][ j ][ km-2 ];		// extrapolation
//			v.x[ i ][ j ][ 0 ] = v.x[ i ][ j ][ km-1 ] = ( v.x[ i ][ j ][ 0 ] + v.x[ i ][ j ][ km-1 ] ) / 2.;

//			w.x[ i ][ j ][ 0 ] = c43 * w.x[ i ][ j ][ 1 ] - c13 * w.x[ i ][ j ][ 2 ];
//			w.x[ i ][ j ][ km-1 ] = c43 * w.x[ i ][ j ][ km-2 ] - c13 * w.x[ i ][ j ][ km-3 ];

			w.x[ i ][ j ][ 0 ] = w.x[ i ][ j ][ 3 ] - 3. * w.x[ i ][ j ][ 2 ] + 3. * w.x[ i ][ j ][ 1 ];		// extrapolation
			w.x[ i ][ j ][ km-1 ] = w.x[ i ][ j ][ km-4 ] - 3. * w.x[ i ][ j ][ km-3 ] + 3. * w.x[ i ][ j ][ km-2 ];		// extrapolation
//			w.x[ i ][ j ][ 0 ] = w.x[ i ][ j ][ km-1 ] = ( w.x[ i ][ j ][ 0 ] + w.x[ i ][ j ][ km-1 ] ) / 2.;

//			p_dyn.x[ i ][ j ][ 0 ] = c43 * p_dyn.x[ i ][ j ][ 1 ] - c13 * p_dyn.x[ i ][ j ][ 2 ];
//			p_dyn.x[ i ][ j ][ km-1 ] = c43 * p_dyn.x[ i ][ j ][ km-2 ] - c13 * p_dyn.x[ i ][ j ][ km-3 ];

			p_dyn.x[ i ][ j ][ 0 ] = p_dyn.x[ i ][ j ][ 3 ] - 3. * p_dyn.x[ i ][ j ][ 2 ] + 3. * p_dyn.x[ i ][ j ][ 1 ];		// extrapolation
			p_dyn.x[ i ][ j ][ km-1 ] = p_dyn.x[ i ][ j ][ km-4 ] - 3. * p_dyn.x[ i ][ j ][ km-3 ] + 3. * p_dyn.x[ i ][ j ][ km-2 ];		// extrapolation
//			p_dyn.x[ i ][ j ][ 0 ] = p_dyn.x[ i ][ j ][ km-1 ] = ( p_dyn.x[ i ][ j ][ 0 ] + p_dyn.x[ i ][ j ][ km-1 ] ) / 2.;

//			c.x[ i ][ j ][ 0 ] = c43 * c.x[ i ][ j ][ 1 ] - c13 * c.x[ i ][ j ][ 2 ];
//			c.x[ i ][ j ][ km-1 ] = c43 * c.x[ i ][ j ][ km-2 ] - c13 * c.x[ i ][ j ][ km-3 ];

			c.x[ i ][ j ][ 0 ] = c.x[ i ][ j ][ 3 ] - 3. * c.x[ i ][ j ][ 2 ] + 3. * c.x[ i ][ j ][ 1 ];		// extrapolation
			c.x[ i ][ j ][ km-1 ] = c.x[ i ][ j ][ km-4 ] - 3. * c.x[ i ][ j ][ km-3 ] + 3. * c.x[ i ][ j ][ km-2 ];		// extrapolation
//			c.x[ i ][ j ][ 0 ] = c.x[ i ][ j ][ km-1 ] = ( c.x[ i ][ j ][ 0 ] + c.x[ i ][ j ][ km-1 ] ) / 2.;

//			co2.x[ i ][ j ][ 0 ] = c43 * co2.x[ i ][ j ][ 1 ] - c13 * co2.x[ i ][ j ][ 2 ];
//			co2.x[ i ][ j ][ km-1 ] = c43 * co2.x[ i ][ j ][ km-2 ] - c13 * co2.x[ i ][ j ][ km-3 ];

			co2.x[ i ][ j ][ 0 ] = co2.x[ i ][ j ][ 3 ] - 3. * co2.x[ i ][ j ][ 2 ] + 3. * co2.x[ i ][ j ][ 1 ];		// extrapolation
			co2.x[ i ][ j ][ km-1 ] = co2.x[ i ][ j ][ km-4 ] - 3. * co2.x[ i ][ j ][ km-3 ] + 3. * co2.x[ i ][ j ][ km-2 ];		// extrapolation
//			co2.x[ i ][ j ][ 0 ] = co2.x[ i ][ j ][ km-1 ] = ( co2.x[ i ][ j ][ 0 ] + co2.x[ i ][ j ][ km-1 ] ) / 2.;

		}
	}
}
