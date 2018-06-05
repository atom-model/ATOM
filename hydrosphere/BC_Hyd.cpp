/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the coordinate system
*/


#include <iostream>
#include <cmath>

#include "BC_Hyd.h"

using namespace std;




BC_Hydrosphere::BC_Hydrosphere ( int im, int jm, int km )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;

	c43 = 4./3.;
	c13 = 1./3.;
}


BC_Hydrosphere::~BC_Hydrosphere() {}




void BC_Hydrosphere::RB_radius ( double ca, double ta, double pa, double dr, Array_1D &rad, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c )
{
// boundary conditions for the r-direction, loop index i

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
// temperature und salinity at sea level are constant values
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			u.x[ 0 ][ j ][ k ] = c43 * u.x[ 1 ][ j ][ k ] - c13 * u.x[ 2 ][ j ][ k ];		// Neumann
			v.x[ 0 ][ j ][ k ] = c43 * v.x[ 1 ][ j ][ k ] - c13 * v.x[ 2 ][ j ][ k ];		// Neumann
			w.x[ 0 ][ j ][ k ] = c43 * w.x[ 1 ][ j ][ k ] - c13 * w.x[ 2 ][ j ][ k ];		// Neumann
			p_dyn.x[ 0 ][ j ][ k ] = c43 * p_dyn.x[ 1 ][ j ][ k ] - c13 * p_dyn.x[ 2 ][ j ][ k ];		// Neumann

//			u.x[ 0 ][ j ][ k ] = u.x[ 3 ][ j ][ k ] - 3. * u.x[ 2 ][ j ][ k ] + 3. * u.x[ 1 ][ j ][ k ];		// extrapolation
//			v.x[ 0 ][ j ][ k ] = v.x[ 3 ][ j ][ k ] - 3. * v.x[ 2 ][ j ][ k ] + 3. * v.x[ 1 ][ j ][ k ];		// extrapolation
//			w.x[ 0 ][ j ][ k ] = w.x[ 3 ][ j ][ k ] - 3. * w.x[ 2 ][ j ][ k ] + 3. * w.x[ 1 ][ j ][ k ];		// extrapolation
//			t.x[ 0 ][ j ][ k ] = t.x[ 3 ][ j ][ k ] - 3. * t.x[ 2 ][ j ][ k ] + 3. * t.x[ 1 ][ j ][ k ];		// extrapolation
//			p_dyn.x[ 0 ][ j ][ k ] = p_dyn.x[ 3 ][ j ][ k ] - 3. * p_dyn.x[ 2 ][ j ][ k ] + 3. * p_dyn.x[ 1 ][ j ][ k ];		// extrapolation

			t.x[ 0 ][ j ][ k ] = ta;		// Dirichlet
			c.x[ 0 ][ j ][ k ] = ca;		// Dirichlet


			u.x[ im-1 ][ j ][ k ] = 0.;
//			u.x[ im-1 ][ j ][ k ] = c43 * u.x[ im-2 ][ j ][ k ] - c13 * u.x[ im-3 ][ j ][ k ];		// Neumann
			v.x[ im-1 ][ j ][ k ] = c43 * v.x[ im-2 ][ j ][ k ] - c13 * v.x[ im-3 ][ j ][ k ];		// Neumann
			w.x[ im-1 ][ j ][ k ] = c43 * w.x[ im-2 ][ j ][ k ] - c13 * w.x[ im-3 ][ j ][ k ];		// Neumann
			c.x[ im-1 ][ j ][ k ] = c43 * c.x[ im-2 ][ j ][ k ] - c13 * c.x[ im-3 ][ j ][ k ];		// Neumann
//			p_dyn.x[ im-1 ][ j ][ k ] = c43 * p_dyn.x[ im-2 ][ j ][ k ] - c13 * p_dyn.x[ im-3 ][ j ][ k ];		// Neumann

//			t.x[ im-1 ][ j ][ k ] = t.x[ im-4 ][ j ][ k ] - 3. * t.x[ im-3 ][ j ][ k ] + 3. * t.x[ im-2 ][ j ][ k ];		// extrapolation
//			u.x[ im-1 ][ j ][ k ] = u.x[ im-4 ][ j ][ k ] - 3. * u.x[ im-3 ][ j ][ k ] + 3. * u.x[ im-2 ][ j ][ k ];		// extrapolation
//			v.x[ im-1 ][ j ][ k ] = v.x[ im-4 ][ j ][ k ] - 3. * v.x[ im-3 ][ j ][ k ] + 3. * v.x[ im-2 ][ j ][ k ];		// extrapolation
//			w.x[ im-1 ][ j ][ k ] = w.x[ im-4 ][ j ][ k ] - 3. * w.x[ im-3 ][ j ][ k ] + 3. * w.x[ im-2 ][ j ][ k ];		// extrapolation
//			p_dyn.x[ im-1 ][ j ][ k ] = p_dyn.x[ im-4 ][ j ][ k ] - 3. * p_dyn.x[ im-3 ][ j ][ k ] + 3. * p_dyn.x[ im-2 ][ j ][ k ];		// extrapolation
//			c.x[ im-1 ][ j ][ k ] = c.x[ im-4 ][ j ][ k ] - 3. * c.x[ im-3 ][ j ][ k ] + 3. * c.x[ im-2 ][ j ][ k ];		// extrapolation

		}
	}
}






void BC_Hydrosphere::RB_theta ( double ca, double ta, double pa, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c )
{
// boundary conditions for the the-direction, loop index j

//	for ( int i = 0; i < im-1; i++ )
	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			t.x[ i ][ 0 ][ k ] = c43 * t.x[ i ][ 1 ][ k ] - c13 * t.x[ i ][ 2 ][ k ];
			t.x[ i ][ jm-1 ][ k ] = c43 * t.x[ i ][ jm-2 ][ k ] - c13 * t.x[ i ][ jm-3 ][ k ];

			u.x[ i ][ 0 ][ k ] = c43 * u.x[ i ][ 1 ][ k ] - c13 * u.x[ i ][ 2 ][ k ];
			u.x[ i ][ jm-1 ][ k ] = c43 * u.x[ i ][ jm-2 ][ k ] - c13 * u.x[ i ][ jm-3 ][ k ];

			v.x[ i ][ 0 ][ k ] = 0.;
			v.x[ i ][ jm-1 ][ k ] = 0.;

			w.x[ i ][ 0 ][ k ] = 0.;
			w.x[ i ][ jm-1 ][ k ] = 0.;

			p_dyn.x[ i ][ 0 ][ k ] = c43 * p_dyn.x[ i ][ 1 ][ k ] - c13 * p_dyn.x[ i ][ 2 ][ k ];
			p_dyn.x[ i ][ jm-1 ][ k ] = c43 * p_dyn.x[ i ][ jm-2 ][ k ] - c13 * p_dyn.x[ i ][ jm-3 ][ k ];

			c.x[ i ][ 0 ][ k ] = c43 * c.x[ i ][ 1 ][ k ] - c13 * c.x[ i ][ 2 ][ k ];
			c.x[ i ][ jm-1 ][ k ] = c43 * c.x[ i ][ jm-2 ][ k ] - c13 * c.x[ i ][ jm-3 ][ k ];
		}
	}
}






void BC_Hydrosphere::RB_phi ( Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c )
{
// boundary conditions for the phi-direction, loop index k
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )
			t.x[ i ][ j ][ 0 ] = c43 * t.x[ i ][ j ][ 1 ] - c13 * t.x[ i ][ j ][ 2 ];
			t.x[ i ][ j ][ km-1 ] = c43 * t.x[ i ][ j ][ km-2 ] - c13 * t.x[ i ][ j ][ km-3 ];
			t.x[ i ][ j ][ 0 ] = t.x[ i ][ j ][ km-1 ] = ( t.x[ i ][ j ][ 0 ] + t.x[ i ][ j ][ km-1 ] ) / 2.;

			u.x[ i ][ j ][ 0 ] = c43 * u.x[ i ][ j ][ 1 ] - c13 * u.x[ i ][ j ][ 2 ];
			u.x[ i ][ j ][ km-1 ] = c43 * u.x[ i ][ j ][ km-2 ] - c13 * u.x[ i ][ j ][ km-3 ];
//			u.x[ i ][ j ][ 0 ] = u.x[ i ][ j ][ 1 ];
//			u.x[ i ][ j ][ km-1 ] = u.x[ i ][ j ][ km-2 ];
			u.x[ i ][ j ][ 0 ] = u.x[ i ][ j ][ km-1 ] = ( u.x[ i ][ j ][ 0 ] + u.x[ i ][ j ][ km-1 ] ) / 2.;

			v.x[ i ][ j ][ 0 ] = c43 * v.x[ i ][ j ][ 1 ] - c13 * v.x[ i ][ j ][ 2 ];
			v.x[ i ][ j ][ km-1 ] = c43 * v.x[ i ][ j ][ km-2 ] - c13 * v.x[ i ][ j ][ km-3 ];
//			v.x[ i ][ j ][ 0 ] = v.x[ i ][ j ][ 1 ];
//			v.x[ i ][ j ][ km-1 ] = v.x[ i ][ j ][ km-2 ];
			v.x[ i ][ j ][ 0 ] = v.x[ i ][ j ][ km-1 ] = ( v.x[ i ][ j ][ 0 ] + v.x[ i ][ j ][ km-1 ] ) / 2.;

			w.x[ i ][ j ][ 0 ] = c43 * w.x[ i ][ j ][ 1 ] - c13 * w.x[ i ][ j ][ 2 ];
			w.x[ i ][ j ][ km-1 ] = c43 * w.x[ i ][ j ][ km-2 ] - c13 * w.x[ i ][ j ][ km-3 ];
//			w.x[ i ][ j ][ 0 ] = w.x[ i ][ j ][ 1 ];
//			w.x[ i ][ j ][ km-1 ] = w.x[ i ][ j ][ km-2 ];
			w.x[ i ][ j ][ 0 ] = w.x[ i ][ j ][ km-1 ] = ( w.x[ i ][ j ][ 0 ] + w.x[ i ][ j ][ km-1 ] ) / 2.;

			p_dyn.x[ i ][ j ][ 0 ] = c43 * p_dyn.x[ i ][ j ][ 1 ] - c13 * p_dyn.x[ i ][ j ][ 2 ];
			p_dyn.x[ i ][ j ][ km-1 ] = c43 * p_dyn.x[ i ][ j ][ km-2 ] - c13 * p_dyn.x[ i ][ j ][ km-3 ];
//			p_dyn.x[ i ][ j ][ 0 ] = p_dyn.x[ i ][ j ][ 1 ];
//			p_dyn.x[ i ][ j ][ km-1 ] = p_dyn.x[ i ][ j ][ km-2 ];
			p_dyn.x[ i ][ j ][ 0 ] = p_dyn.x[ i ][ j ][ km-1 ] = ( p_dyn.x[ i ][ j ][ 0 ] + p_dyn.x[ i ][ j ][ km-1 ] ) / 2.;

			c.x[ i ][ j ][ 0 ] = c43 * c.x[ i ][ j ][ 1 ] - c13 * c.x[ i ][ j ][ 2 ];
			c.x[ i ][ j ][ km-1 ] = c43 * c.x[ i ][ j ][ km-2 ] - c13 * c.x[ i ][ j ][ km-3 ];
			c.x[ i ][ j ][ 0 ] = c.x[ i ][ j ][ km-1 ] = ( c.x[ i ][ j ][ 0 ] + c.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}
}

