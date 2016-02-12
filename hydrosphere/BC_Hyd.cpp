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




void BC_Hydrosphere::RB_radius ( double ca, double ta, double pa, double dr, Array_1D &rad, Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_t, Array &rhs_c, Array &aux_u, Array &aux_v, Array &aux_w, Array &h, Array &Salt_Finger, Array &Salt_Diffusion, Array &Salt_Balance )
{
// boundary conditions for the r-direction, loop index i

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
// temperature und salinity at sea level are constant values
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )

			t.x[ 0 ][ j ][ k ] = ta;
			c.x[ 0 ][ j ][ k ] = ca;
			p.x[ 0 ][ j ][ k ] = pa;
			p.x[ im-1 ][ j ][ k ] = pa;

			u.x[ 0 ][ j ][ k ] = 0.;
			u.x[ im-1 ][ j ][ k ] = 0.;

			v.x[ 0 ][ j ][ k ] = 0.;
			w.x[ 0 ][ j ][ k ] = 0.;


// open development of the velocity components through von Neumann condition, dv/dr = dw/dr = 0 
//			v.x[ im-1 ][ j ][ k ] = c43 * v.x[ im-2 ][ j ][ k ] - c13 * v.x[ im-3 ][ j ][ k ];
//			w.x[ im-1 ][ j ][ k ] = c43 * w.x[ im-2 ][ j ][ k ] - c13 * w.x[ im-3 ][ j ][ k ];
//			v.x[ im-1 ][ j ][ k ] = 2. * dr / rad.z[ im-2 ] * v.x[ im-2 ][ j ][ k ] + v.x[ im-3 ][ j ][ k ];
//			w.x[ im-1 ][ j ][ k ] = 2. * dr / rad.z[ im-2 ] * w.x[ im-2 ][ j ][ k ] + w.x[ im-3 ][ j ][ k ];

//			v.x[ im-1 ][ j ][ k ] = c43 * v.x[ im-2 ][ j ][ k ] - c13 * v.x[ im-3 ][ j ][ k ];
//			w.x[ im-1 ][ j ][ k ] = c43 * w.x[ im-2 ][ j ][ k ] - c13 * w.x[ im-3 ][ j ][ k ];

// *********************  various boundary conditions at sea level ***************

//			v.x[ im-1 ][ j ][ k ] = 2. * dr / rad.z[ im-2 ] * v.x[ im-2 ][ j ][ k ] + v.x[ im-3 ][ j ][ k ];						// radial influence
//			w.x[ im-1 ][ j ][ k ] = 2. * dr / rad.z[ im-2 ] * w.x[ im-2 ][ j ][ k ] + w.x[ im-3 ][ j ][ k ];					// radial influence

//			v.x[ im-1 ][ j ][ k ] = c43 * v.x[ im-2 ][ j ][ k ] - c13 * v.x[ im-3 ][ j ][ k ];										// von Neumann ( not to be used )
//			w.x[ im-1 ][ j ][ k ] = c43 * w.x[ im-2 ][ j ][ k ] - c13 * w.x[ im-3 ][ j ][ k ];									// leads to a spreadig of north/south equatorial currents

//			v.x[ im-1 ][ j ][ k ] = v.x[ im-4 ][ j ][ k ] - 3. * v.x[ im-3 ][ j ][ k ] + 3. * v.x[ im-2 ][ j ][ k ];				// extrapolation
//			w.x[ im-1 ][ j ][ k ] = w.x[ im-4 ][ j ][ k ] - 3. * w.x[ im-3 ][ j ][ k ] + 3. * w.x[ im-2 ][ j ][ k ];			// leads to round harmonic profiles

//			p.x[ 0 ][ j ][ k ] = c43 * p.x[ 1 ][ j ][ k ] - c13 * p.x[ 2 ][ j ][ k ];
//			p.x[ im-1 ][ j ][ k ] = c43 * p.x[ im-2 ][ j ][ k ] - c13 * p.x[ im-3 ][ j ][ k ];

			rhs_u.x[ 0 ][ j ][ k ] = c43 * rhs_u.x[ 1 ][ j ][ k ] - c13 * rhs_u.x[ 2 ][ j ][ k ];
			rhs_u.x[ im-1 ][ j ][ k ] = c43 * rhs_u.x[ im-2 ][ j ][ k ] - c13 * rhs_u.x[ im-3 ][ j ][ k ];

			rhs_v.x[ 0 ][ j ][ k ] = c43 * rhs_v.x[ 1 ][ j ][ k ] - c13 * rhs_v.x[ 2 ][ j ][ k ];
			rhs_v.x[ im-1 ][ j ][ k ] = c43 * rhs_v.x[ im-2 ][ j ][ k ] - c13 * rhs_v.x[ im-3 ][ j ][ k ];

			rhs_w.x[ 0 ][ j ][ k ] = c43 * rhs_w.x[ 1 ][ j ][ k ] - c13 * rhs_w.x[ 2 ][ j ][ k ];
			rhs_w.x[ im-1 ][ j ][ k ] = c43 * rhs_w.x[ im-2 ][ j ][ k ] - c13 * rhs_w.x[ im-3 ][ j ][ k ];

			rhs_t.x[ 0 ][ j ][ k ] = c43 * rhs_t.x[ 1 ][ j ][ k ] - c13 * rhs_t.x[ 2 ][ j ][ k ];
			rhs_t.x[ im-1 ][ j ][ k ] = c43 * rhs_t.x[ im-2 ][ j ][ k ] - c13 * rhs_t.x[ im-3 ][ j ][ k ];

			rhs_c.x[ 0 ][ j ][ k ] = c43 * rhs_c.x[ 1 ][ j ][ k ] - c13 * rhs_c.x[ 2 ][ j ][ k ];
			rhs_c.x[ im-1 ][ j ][ k ] = c43 * rhs_c.x[ im-2 ][ j ][ k ] - c13 * rhs_c.x[ im-3 ][ j ][ k ];

			aux_u.x[ 0 ][ j ][ k ] = c43 * aux_u.x[ 1 ][ j ][ k ] - c13 * aux_u.x[ 2 ][ j ][ k ];
			aux_u.x[ im-1 ][ j ][ k ] = c43 * aux_u.x[ im-2 ][ j ][ k ] - c13 * aux_u.x[ im-3 ][ j ][ k ];

			aux_v.x[ 0 ][ j ][ k ] = c43 * aux_v.x[ 1 ][ j ][ k ] - c13 * aux_v.x[ 2 ][ j ][ k ];
			aux_v.x[ im-1 ][ j ][ k ] = c43 * aux_v.x[ im-2 ][ j ][ k ] - c13 * aux_v.x[ im-3 ][ j ][ k ];

			aux_w.x[ 0 ][ j ][ k ] = c43 * aux_w.x[ 1 ][ j ][ k ] - c13 * aux_w.x[ 2 ][ j ][ k ];
			aux_w.x[ im-1 ][ j ][ k ] = c43 * aux_w.x[ im-2 ][ j ][ k ] - c13 * aux_w.x[ im-3 ][ j ][ k ];

			Salt_Balance.x[ 0 ][ j ][ k ] = c43 * Salt_Balance.x[ 1 ][ j ][ k ] - c13 * Salt_Balance.x[ 2 ][ j ][ k ];
			Salt_Diffusion.x[ 0 ][ j ][ k ] = c43 * Salt_Diffusion.x[ 1 ][ j ][ k ] - c13 * Salt_Diffusion.x[ 2 ][ j ][ k ];
			Salt_Finger.x[ 0 ][ j ][ k ] = c43 * Salt_Finger.x[ 1 ][ j ][ k ] - c13 * Salt_Finger.x[ 2 ][ j ][ k ];

			Salt_Balance.x[ im-1 ][ j ][ k ] = c43 * Salt_Balance.x[ im-2 ][ j ][ k ] - c13 * Salt_Balance.x[ im-3 ][ j ][ k ];
			Salt_Diffusion.x[ im-1 ][ j ][ k ] = c43 * Salt_Diffusion.x[ im-2 ][ j ][ k ] - c13 * Salt_Diffusion.x[ im-3 ][ j ][ k ];
			Salt_Finger.x[ im-1 ][ j ][ k ] = c43 * Salt_Finger.x[ im-2 ][ j ][ k ] - c13 * Salt_Finger.x[ im-3 ][ j ][ k ];

		}
	}
}






void BC_Hydrosphere::RB_theta ( double ca, double ta, double pa, Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_t, Array &rhs_c, Array &aux_u, Array &aux_v, Array &aux_w, Array &h, Array &Salt_Finger, Array &Salt_Diffusion, Array &Salt_Balance )
{
// boundary conditions for the the-direction, loop index j

	for ( int i = 0; i < im-1; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
//			if ( h.x[ i ][ jm-1 ][ k ] != 1. )
			{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )

				t.x[ i ][ 0 ][ k ] = c43 * t.x[ i ][ 1 ][ k ] - c13 * t.x[ i ][ 2 ][ k ];
				t.x[ i ][ jm-1 ][ k ] = c43 * t.x[ i ][ jm-2 ][ k ] - c13 * t.x[ i ][ jm-3 ][ k ];

				u.x[ i ][ 0 ][ k ] = c43 * u.x[ i ][ 1 ][ k ] - c13 * u.x[ i ][ 2 ][ k ];
				u.x[ i ][ jm-1 ][ k ] = c43 * u.x[ i ][ jm-2 ][ k ] - c13 * u.x[ i ][ jm-3 ][ k ];

//				u.x[ i ][ 0 ][ k ] = 0.;
//				u.x[ i ][ jm-1 ][ k ] = 0.;

				v.x[ i ][ 0 ][ k ] = c43 * v.x[ i ][ 1 ][ k ] - c13 * v.x[ i ][ 2 ][ k ];
				v.x[ i ][ jm-1 ][ k ] = c43 * v.x[ i ][ jm-2 ][ k ] - c13 * v.x[ i ][ jm-3 ][ k ];

				w.x[ i ][ 0 ][ k ] = c43 * w.x[ i ][ 1 ][ k ] - c13 * w.x[ i ][ 2 ][ k ];
				w.x[ i ][ jm-1 ][ k ] = c43 * w.x[ i ][ jm-2 ][ k ] - c13 * w.x[ i ][ jm-3 ][ k ];

//				v.x[ i ][ 0 ][ k ] = 0.;
//				v.x[ i ][ jm-1 ][ k ] = 0.;

//				w.x[ i ][ 0 ][ k ] = 0.;
//				w.x[ i ][ jm-1 ][ k ] = 0.;

				p.x[ i ][ 0 ][ k ] = c43 * p.x[ i ][ 1 ][ k ] - c13 * p.x[ i ][ 2 ][ k ];
				p.x[ i ][ jm-1 ][ k ] = c43 * p.x[ i ][ jm-2 ][ k ] - c13 * p.x[ i ][ jm-3 ][ k ];

				c.x[ i ][ 0 ][ k ] = c43 * c.x[ i ][ 1 ][ k ] - c13 * c.x[ i ][ 2 ][ k ];
				c.x[ i ][ jm-1 ][ k ] = c43 * c.x[ i ][ jm-2 ][ k ] - c13 * c.x[ i ][ jm-3 ][ k ];

				rhs_u.x[ i ][ 0 ][ k ] = c43 * rhs_u.x[ i ][ 1 ][ k ] - c13 * rhs_u.x[ i ][ 2 ][ k ];
				rhs_u.x[ i ][ jm-1 ][ k ] = c43 * rhs_u.x[ i ][ jm-2 ][ k ] - c13 * rhs_u.x[ i ][ jm-3 ][ k ];

				rhs_v.x[ i ][ 0 ][ k ] = c43 * rhs_v.x[ i ][ 1 ][ k ] - c13 * rhs_v.x[ i ][ 2 ][ k ];
				rhs_v.x[ i ][ jm-1 ][ k ] = c43 * rhs_v.x[ i ][ jm-2 ][ k ] - c13 * rhs_v.x[ i ][ jm-3 ][ k ];

				rhs_w.x[ i ][ 0 ][ k ] = c43 * rhs_w.x[ i ][ 1 ][ k ] - c13 * rhs_w.x[ i ][ 2 ][ k ];
				rhs_w.x[ i ][ jm-1 ][ k ] = c43 * rhs_w.x[ i ][ jm-2 ][ k ] - c13 * rhs_w.x[ i ][ jm-3 ][ k ];

				rhs_t.x[ i ][ 0 ][ k ] = c43 * rhs_t.x[ i ][ 1 ][ k ] - c13 * rhs_t.x[ i ][ 2 ][ k ];
				rhs_t.x[ i ][ jm-1 ][ k ] = c43 * rhs_t.x[ i ][ jm-2 ][ k ] - c13 * rhs_t.x[ i ][ jm-3 ][ k ];

				rhs_c.x[ i ][ 0 ][ k ] = c43 * rhs_c.x[ i ][ 1 ][ k ] - c13 * rhs_c.x[ i ][ 2 ][ k ];
				rhs_c.x[ i ][ jm-1 ][ k ] = c43 * rhs_c.x[ i ][ jm-2 ][ k ] - c13 * rhs_c.x[ i ][ jm-3 ][ k ];

				aux_u.x[ i ][ 0 ][ k ] = c43 * aux_u.x[ i ][ 1 ][ k ] - c13 * aux_u.x[ i ][ 2 ][ k ];
				aux_u.x[ i ][ jm-1 ][ k ] = c43 * aux_u.x[ i ][ jm-2 ][ k ] - c13 * aux_u.x[ i ][ jm-3 ][ k ];

				aux_v.x[ i ][ 0 ][ k ] = c43 * aux_v.x[ i ][ 1 ][ k ] - c13 * aux_v.x[ i ][ 2 ][ k ];
				aux_v.x[ i ][ jm-1 ][ k ] = c43 * aux_v.x[ i ][ jm-2 ][ k ] - c13 * aux_v.x[ i ][ jm-3 ][ k ];

				aux_w.x[ i ][ 0 ][ k ] = c43 * aux_w.x[ i ][ 1 ][ k ] - c13 * aux_w.x[ i ][ 2 ][ k ];
				aux_w.x[ i ][ jm-1 ][ k ] = c43 * aux_w.x[ i ][ jm-2 ][ k ] - c13 * aux_w.x[ i ][ jm-3 ][ k ];


				Salt_Balance.x[ i ][ 0 ][ k ] = c43 * Salt_Balance.x[ i ][ 1 ][ k ] - c13 * Salt_Balance.x[ i ][ 2 ][ k ];
				Salt_Balance.x[ i ][ jm-1 ][ k ] = c43 * Salt_Balance.x[ i ][ jm-2 ][ k ] - c13 * Salt_Balance.x[ i ][ jm-3 ][ k ];

				Salt_Diffusion.x[ i ][ 0 ][ k ] = c43 * Salt_Diffusion.x[ i ][ 1 ][ k ] - c13 * Salt_Diffusion.x[ i ][ 2 ][ k ];
				Salt_Diffusion.x[ i ][ jm-1 ][ k ] = c43 * Salt_Diffusion.x[ i ][ jm-2 ][ k ] - c13 * Salt_Diffusion.x[ i ][ jm-3 ][ k ];

				Salt_Finger.x[ i ][ 0 ][ k ] = c43 * Salt_Finger.x[ i ][ 1 ][ k ] - c13 * Salt_Finger.x[ i ][ 2 ][ k ];
				Salt_Finger.x[ i ][ jm-1 ][ k ] = c43 * Salt_Finger.x[ i ][ jm-2 ][ k ] - c13 * Salt_Finger.x[ i ][ jm-3 ][ k ];


/*
				t.x[ i ][ 0 ][ k ] = ta;
				t.x[ i ][ jm-1 ][ k ] = ta;

				u.x[ i ][ 0 ][ k ] = 0.;
				u.x[ i ][ jm-1 ][ k ] = 0.;

				v.x[ i ][ 0 ][ k ] = 0.;
				v.x[ i ][ jm-1 ][ k ] = 0.;

				w.x[ i ][ 0 ][ k ] = 0.;
				w.x[ i ][ jm-1 ][ k ] = 0.;

				p.x[ i ][ 0 ][ k ] = pa;
				p.x[ i ][ jm-1 ][ k ] = pa;
*/
//				c.x[ i ][ 0 ][ k ] = ca;
//				c.x[ i ][ jm-1 ][ k ] = ca;
/*
				rhs_u.x[ i ][ 0 ][ k ] = 0.;
				rhs_u.x[ i ][ jm-1 ][ k ] = 0.;

				rhs_v.x[ i ][ 0 ][ k ] = 0.;
				rhs_v.x[ i ][ jm-1 ][ k ] = 0.;

				rhs_w.x[ i ][ 0 ][ k ] = 0.;
				rhs_w.x[ i ][ jm-1 ][ k ] = 0.;

				rhs_t.x[ i ][ 0 ][ k ] = 0.;
				rhs_t.x[ i ][ jm-1 ][ k ] = 0.;

				rhs_c.x[ i ][ 0 ][ k ] = 0.;
				rhs_c.x[ i ][ jm-1 ][ k ] = 0.;

				aux_u.x[ i ][ 0 ][ k ] = 0.;
				aux_u.x[ i ][ jm-1 ][ k ] = 0.;

				aux_v.x[ i ][ 0 ][ k ] = 0.;
				aux_v.x[ i ][ jm-1 ][ k ] = 0.;

				aux_w.x[ i ][ 0 ][ k ] = 0.;
				aux_w.x[ i ][ jm-1 ][ k ] = 0.;

				Salt_Balance.x[ i ][ 0 ][ k ] = 0.;
				Salt_Balance.x[ i ][ jm-1 ][ k ] = 0.;

				Salt_Diffusion.x[ i ][ 0 ][ k ] = 0.;
				Salt_Diffusion.x[ i ][ jm-1 ][ k ] = 0.;

				Salt_Finger.x[ i ][ 0 ][ k ] = 0.;
				Salt_Finger.x[ i ][ jm-1 ][ k ] = 0.;
*/
			}
		}
	}
}






void BC_Hydrosphere::RB_phi ( Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_t, Array &rhs_c, Array &aux_u, Array &aux_v, Array &aux_w, Array &h, Array &Salt_Finger, Array &Salt_Diffusion, Array &Salt_Balance )
{
// boundary conditions for the phi-direction, loop index k

	for ( int i = 0; i < im-1; i++ )
//	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
//			if ( h.x[ i ][ j ][ km-1 ] != 1. )
			{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )

				t.x[ i ][ j ][ 0 ] = c43 * t.x[ i ][ j ][ 1 ] - c13 * t.x[ i ][ j ][ 2 ];
				t.x[ i ][ j ][ km-1 ] = c43 * t.x[ i ][ j ][ km-2 ] - c13 * t.x[ i ][ j ][ km-3 ];
				t.x[ i ][ j ][ 0 ] = t.x[ i ][ j ][ km-1 ] = ( t.x[ i ][ j ][ 0 ] + t.x[ i ][ j ][ km-1 ] ) / 2.;

				u.x[ i ][ j ][ 0 ] = c43 * u.x[ i ][ j ][ 1 ] - c13 * u.x[ i ][ j ][ 2 ];
				u.x[ i ][ j ][ km-1 ] = c43 * u.x[ i ][ j ][ km-2 ] - c13 * u.x[ i ][ j ][ km-3 ];
				u.x[ i ][ j ][ 0 ] = u.x[ i ][ j ][ km-1 ] = ( u.x[ i ][ j ][ 0 ] + u.x[ i ][ j ][ km-1 ] ) / 2.;

				v.x[ i ][ j ][ 0 ] = c43 * v.x[ i ][ j ][ 1 ] - c13 * v.x[ i ][ j ][ 2 ];
				v.x[ i ][ j ][ km-1 ] = c43 * v.x[ i ][ j ][ km-2 ] - c13 * v.x[ i ][ j ][ km-3 ];
				v.x[ i ][ j ][ 0 ] = v.x[ i ][ j ][ km-1 ] = ( v.x[ i ][ j ][ 0 ] + v.x[ i ][ j ][ km-1 ] ) / 2.;

				w.x[ i ][ j ][ 0 ] = c43 * w.x[ i ][ j ][ 1 ] - c13 * w.x[ i ][ j ][ 2 ];
				w.x[ i ][ j ][ km-1 ] = c43 * w.x[ i ][ j ][ km-2 ] - c13 * w.x[ i ][ j ][ km-3 ];
				w.x[ i ][ j ][ 0 ] = w.x[ i ][ j ][ km-1 ] = ( w.x[ i ][ j ][ 0 ] + w.x[ i ][ j ][ km-1 ] ) / 2.;

				p.x[ i ][ j ][ 0 ] = c43 * p.x[ i ][ j ][ 1 ] - c13 * p.x[ i ][ j ][ 2 ];
				p.x[ i ][ j ][ km-1 ] = c43 * p.x[ i ][ j ][ km-2 ] - c13 * p.x[ i ][ j ][ km-3 ];
				p.x[ i ][ j ][ 0 ] = p.x[ i ][ j ][ km-1 ] = ( p.x[ i ][ j ][ 0 ] + p.x[ i ][ j ][ km-1 ] ) / 2.;

				c.x[ i ][ j ][ 0 ] = c43 * c.x[ i ][ j ][ 1 ] - c13 * c.x[ i ][ j ][ 2 ];
				c.x[ i ][ j ][ km-1 ] = c43 * c.x[ i ][ j ][ km-2 ] - c13 * c.x[ i ][ j ][ km-3 ];
				c.x[ i ][ j ][ 0 ] = c.x[ i ][ j ][ km-1 ] = ( c.x[ i ][ j ][ 0 ] + c.x[ i ][ j ][ km-1 ] ) / 2.;

				rhs_u.x[ i ][ j ][ 0 ] = c43 * rhs_u.x[ i ][ j ][ 1 ] - c13 * rhs_u.x[ i ][ j ][ 2 ];
				rhs_u.x[ i ][ j ][ km-1 ] = c43 * rhs_u.x[ i ][ j ][ km-2 ] - c13 * rhs_u.x[ i ][ j ][ km-3 ];
				rhs_u.x[ i ][ j ][ 0 ] = rhs_u.x[ i ][ j ][ km-1 ] = ( rhs_u.x[ i ][ j ][ 0 ] + rhs_u.x[ i ][ j ][ km-1 ] ) / 2.;

				rhs_v.x[ i ][ j ][ 0 ] = c43 * rhs_v.x[ i ][ j ][ 1 ] - c13 * rhs_v.x[ i ][ j ][ 2 ];
				rhs_v.x[ i ][ j ][ km-1 ] = c43 * rhs_v.x[ i ][ j ][ km-2 ] - c13 * rhs_v.x[ i ][ j ][ km-3 ];
				rhs_v.x[ i ][ j ][ 0 ] = rhs_v.x[ i ][ j ][ km-1 ] = ( rhs_v.x[ i ][ j ][ 0 ] + rhs_v.x[ i ][ j ][ km-1 ] ) / 2.;

				rhs_w.x[ i ][ j ][ 0 ] = c43 * rhs_w.x[ i ][ j ][ 1 ] - c13 * rhs_w.x[ i ][ j ][ 2 ];
				rhs_w.x[ i ][ j ][ km-1 ] = c43 * rhs_w.x[ i ][ j ][ km-2 ] - c13 * rhs_w.x[ i ][ j ][ km-3 ];
				rhs_w.x[ i ][ j ][ 0 ] = rhs_w.x[ i ][ j ][ km-1 ] = ( rhs_w.x[ i ][ j ][ 0 ] + rhs_w.x[ i ][ j ][ km-1 ] ) / 2.;

				rhs_t.x[ i ][ j ][ 0 ] = c43 * rhs_t.x[ i ][ j ][ 1 ] - c13 * rhs_t.x[ i ][ j ][ 2 ];
				rhs_t.x[ i ][ j ][ km-1 ] = c43 * rhs_t.x[ i ][ j ][ km-2 ] - c13 * rhs_t.x[ i ][ j ][ km-3 ];
				rhs_t.x[ i ][ j ][ 0 ] = rhs_t.x[ i ][ j ][ km-1 ] = ( rhs_t.x[ i ][ j ][ 0 ] + rhs_t.x[ i ][ j ][ km-1 ] ) / 2.;

				rhs_c.x[ i ][ j ][ 0 ] = c43 * rhs_c.x[ i ][ j ][ 1 ] - c13 * rhs_c.x[ i ][ j ][ 2 ];
				rhs_c.x[ i ][ j ][ km-1 ] = c43 * rhs_c.x[ i ][ j ][ km-2 ] - c13 * rhs_c.x[ i ][ j ][ km-3 ];
				rhs_c.x[ i ][ j ][ 0 ] = rhs_c.x[ i ][ j ][ km-1 ] = ( rhs_c.x[ i ][ j ][ 0 ] + rhs_c.x[ i ][ j ][ km-1 ] ) / 2.;

				aux_u.x[ i ][ j ][ 0 ] = c43 * aux_u.x[ i ][ j ][ 1 ] - c13 * aux_u.x[ i ][ j ][ 2 ];
				aux_u.x[ i ][ j ][ km-1 ] = c43 * aux_u.x[ i ][ j ][ km-2 ] - c13 * aux_u.x[ i ][ j ][ km-3 ];
				aux_u.x[ i ][ j ][ 0 ] = aux_u.x[ i ][ j ][ km-1 ] = ( aux_u.x[ i ][ j ][ 0 ] + aux_u.x[ i ][ j ][ km-1 ] ) / 2.;

				aux_v.x[ i ][ j ][ 0 ] = c43 * aux_v.x[ i ][ j ][ 1 ] - c13 * aux_v.x[ i ][ j ][ 2 ];
				aux_v.x[ i ][ j ][ km-1 ] = c43 * aux_v.x[ i ][ j ][ km-2 ] - c13 * aux_v.x[ i ][ j ][ km-3 ];
				aux_v.x[ i ][ j ][ 0 ] = aux_v.x[ i ][ j ][ km-1 ] = ( aux_v.x[ i ][ j ][ 0 ] + aux_v.x[ i ][ j ][ km-1 ] ) / 2.;

				aux_w.x[ i ][ j ][ 0 ] = c43 * aux_w.x[ i ][ j ][ 1 ] - c13 * aux_w.x[ i ][ j ][ 2 ];
				aux_w.x[ i ][ j ][ km-1 ] = c43 * aux_w.x[ i ][ j ][ km-2 ] - c13 * aux_w.x[ i ][ j ][ km-3 ];
				aux_w.x[ i ][ j ][ 0 ] = aux_w.x[ i ][ j ][ km-1 ] = ( aux_w.x[ i ][ j ][ 0 ] + aux_w.x[ i ][ j ][ km-1 ] ) / 2.;

				Salt_Balance.x[ i ][ j ][ 0 ] = c43 * Salt_Balance.x[ i ][ j ][ 1 ] - c13 * Salt_Balance.x[ i ][ j ][ 2 ];
				Salt_Balance.x[ i ][ j ][ km-1 ] = c43 * Salt_Balance.x[ i ][ j ][ km-2 ] - c13 * Salt_Balance.x[ i ][ j ][ km-3 ];
				Salt_Balance.x[ i ][ j ][ 0 ] = Salt_Balance.x[ i ][ j ][ km-1 ] = ( Salt_Balance.x[ i ][ j ][ 0 ] + Salt_Balance.x[ i ][ j ][ km-1 ] ) / 2.;

				Salt_Diffusion.x[ i ][ j ][ 0 ] = c43 * Salt_Diffusion.x[ i ][ j ][ 1 ] - c13 * Salt_Diffusion.x[ i ][ j ][ 2 ];
				Salt_Diffusion.x[ i ][ j ][ km-1 ] = c43 * Salt_Diffusion.x[ i ][ j ][ km-2 ] - c13 * Salt_Diffusion.x[ i ][ j ][ km-3 ];
				Salt_Diffusion.x[ i ][ j ][ 0 ] = Salt_Diffusion.x[ i ][ j ][ km-1 ] = ( Salt_Diffusion.x[ i ][ j ][ 0 ] + Salt_Diffusion.x[ i ][ j ][ km-1 ] ) / 2.;

				Salt_Finger.x[ i ][ j ][ 0 ] = c43 * Salt_Finger.x[ i ][ j ][ 1 ] - c13 * Salt_Finger.x[ i ][ j ][ 2 ];
				Salt_Finger.x[ i ][ j ][ km-1 ] = c43 * Salt_Finger.x[ i ][ j ][ km-2 ] - c13 * Salt_Finger.x[ i ][ j ][ km-3 ];
				Salt_Finger.x[ i ][ j ][ 0 ] = Salt_Finger.x[ i ][ j ][ km-1 ] = ( Salt_Finger.x[ i ][ j ][ 0 ] + Salt_Finger.x[ i ][ j ][ km-1 ] ) / 2.;

			}
		}
	}
}




void BC_Hydrosphere::BC_NST_control_2D ( double dr, double dthe, double dphi, double re, double mue_air, double mue_water, Array &h, Array &v, Array &w, Array &p, Array_2D &aux_2D_v, Array_2D &aux_2D_w, Array_1D &rad, Array_1D &the )
{
// boundary conditions for sea surface by an additional volume force

// 1. order derivatives for 3 spacial directions in Finite Difference Methods ( FDM )

//	mue = ( 1. + mue_water / mue_air ) / re;
	mue = 1. / re;

// collection of coefficients
	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

// collection of coefficients
	rm = rad.z[ 0 ];
	rm2 = rm * rm;

	for ( int j = 1; j < jm-1; j++ )
	{
// collection of coefficients
		sinthe = sin( the.z[ j ] );
		sinthe2 = sinthe * sinthe;
		costhe = cos( the.z[ j ] );
		rmsinthe = rm * sinthe;
		rm2sinthe = rm2 * sinthe;
		rm2sinthe2 = rm2 * sinthe2;

		for ( int k = 1; k < km-1; k++ )
		{
			if ( h.x[ 0 ][ j ][ k ] == 0. ) 
			{
				dpdthe = ( p.x[ 0 ][ j + 1 ][ k ] - p.x[ 0 ][ j - 1 ][ k ] ) / ( 2. * dthe );
				dvdthe = ( v.x[ 0 ][ j + 1 ][ k ] - v.x[ 0 ][ j - 1 ][ k ] ) / ( 2. * dthe );
				dwdthe = ( w.x[ 0 ][ j + 1 ][ k ] - w.x[ 0 ][ j - 1 ][ k ] ) / ( 2. * dthe );

				dpdphi = ( p.x[ 0 ][ j ][ k + 1 ] - p.x[ 0 ][ j ][ k - 1 ] ) / ( 2. * dphi );
				dvdphi = ( v.x[ 0 ][ j ][ k + 1 ] - v.x[ 0 ][ j ][ k - 1 ] ) / ( 2. * dphi );
				dwdphi = ( w.x[ 0 ][ j ][ k + 1 ] - w.x[ 0 ][ j ][ k - 1 ] ) / ( 2. * dphi );

				LHS_v = dvdthe / rm + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe2 ) / rm 
								+ 2. * mue / ( rm2sinthe2 * dphi2 ) + 6. * mue / rm / ( 2. * dr ) - mue / dr2;

				RHS_v = mue / ( rm2 * dthe2 ) * ( v.x[ 0 ][ j + 1 ][ k ] + v.x[ 0 ][ j - 1 ][ k ] ) 
								+ mue / ( rm2sinthe2 * dphi2 ) * ( v.x[ 0 ][ j ][ k + 1 ] + v.x[ 0 ][ j ][ k - 1 ] ) 
								+ mue * ( - 2. * v.x[ 1 ][ j ][ k ] + v.x[ 2 ][ j ][ k ] ) / dr2 - w.x[ 0 ][ j ][ k ] / rmsinthe * dvdphi 
								- dpdthe / rm + 2. * mue / rm * ( 4. * v.x[ 1 ][ j ][ k ] - v.x[ 2 ][ j ][ k ] ) / ( 2. * dr ) 
								+ mue * costhe / rm2sinthe * dvdthe - 2. * mue * costhe / rm2sinthe2 * dwdphi;

				aux_2D_v.y[ j ][ k ] = RHS_v / LHS_v;

/*
				LHS_v = dvdthe / rm + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe2 ) / rm 
								+ 2. * mue / ( rm2sinthe2 * dphi2 );

				RHS_v = mue / ( rm2 * dthe2 ) * ( v.x[ 0 ][ j + 1 ][ k ] + v.x[ 0 ][ j - 1 ][ k ] ) 
								+ mue / ( rm2sinthe2 * dphi2 ) * ( v.x[ 0 ][ j ][ k + 1 ] + v.x[ 0 ][ j ][ k - 1 ] ) 
								- w.x[ 0 ][ j ][ k ] / rmsinthe * dvdphi 
								- dpdthe / rm 
								+ mue * costhe / rm2sinthe * dvdthe - 2. * mue * costhe / rm2sinthe2 * dwdphi;

				aux_2D_v.y[ j ][ k ] = RHS_v / LHS_v;
*/

				LHS_w = dwdphi / rmsinthe + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe ) / rm
								+ 2. * mue / ( rm2sinthe2 * dphi2 ) + 6. * mue / rm / ( 2. * dr ) - mue / dr2;

				RHS_w = mue / ( rm2 * dthe2 ) * ( w.x[ 0 ][ j + 1 ][ k ] + w.x[ 0 ][ j - 1 ][ k ] ) 
								+ mue / ( rm2sinthe2 * dphi2 ) * ( w.x[ 0 ][ j ][ k + 1 ] + w.x[ 0 ][ j ][ k - 1 ] ) 
								+ mue * ( - 2. * w.x[ 1 ][ j ][ k ] + w.x[ 2 ][ j ][ k ] ) / dr2 - v.x[ 0 ][ j ][ k ] / rm * dwdthe 
								- dpdphi / rmsinthe + 2. * mue / rm * ( 4. * w.x[ 1 ][ j ][ k ] - w.x[ 2 ][ j ][ k ] ) / ( 2. * dr ) 
								+ mue * costhe / rm2sinthe * dwdthe + 2. * mue * costhe / rm2sinthe2 * dvdphi;

				aux_2D_w.y[ j ][ k ] = RHS_w / LHS_w;

 /*
				LHS_w = dwdphi / rmsinthe + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe ) / rm
								+ 2. * mue / ( rm2sinthe2 * dphi2 );

				RHS_w = mue / ( rm2 * dthe2 ) * ( w.x[ 0 ][ j + 1 ][ k ] + w.x[ 0 ][ j - 1 ][ k ] ) 
								+ mue / ( rm2sinthe2 * dphi2 ) * ( w.x[ 0 ][ j ][ k + 1 ] + w.x[ 0 ][ j ][ k - 1 ] ) 
								- v.x[ 0 ][ j ][ k ] / rm * dwdthe 
								- dpdphi / rmsinthe
								+ mue * costhe / rm2sinthe * dwdthe + 2. * mue * costhe / rm2sinthe2 * dvdphi;

				aux_2D_w.y[ j ][ k ] = RHS_w / LHS_w;
*/
			}
			else
			{
				aux_2D_v.y[ j ][ k ] = 0.;
				aux_2D_w.y[ j ][ k ] = 0.;
			}
		}
	}


	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			v.x[ 0 ][ j ][ k ] = aux_2D_v.y[ j ][ k ];
			w.x[ 0 ][ j ][ k ] = aux_2D_w.y[ j ][ k ];
		}
	}

}






void BC_Hydrosphere::BC_NST_control_3D ( double dr, double dthe, double dphi, double re, double mue_air, double mue_water, Array &h, Array &u, Array &v, Array &w, Array &t, Array &p, Array &c, Array &aux_u, Array &aux_v, Array &aux_w, Array_1D &rad, Array_1D &the )
{
// boundary conditions for sea surface by an additional volume force

// 1. order derivatives for 3 spacial directions in Finite Difference Methods ( FDM )

//	mue = ( 1. + mue_water / mue_air ) / re;
	mue = 1. / re;



// 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )

// collection of coefficients
	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

	for ( int i = 1; i < im-1; i++ )
	{
// collection of coefficients
		rm = rad.z[ i ];
		rm2 = rm * rm;

		for ( int j = 1; j < jm-1; j++ )
		{
// collection of coefficients
			sinthe = sin( the.z[ j ] );
			sinthe2 = sinthe * sinthe;
			costhe = cos( the.z[ j ] );
			cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
			rmsinthe = rm * sinthe;
			rm2sinthe = rm2 * sinthe;
			rm2sinthe2 = rm2 * sinthe2;

			for ( int k = 1; k < km-1; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. ) 
				{
// 1st order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
					dudr = ( u.x[ i+1 ][ j ][ k ] - u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dvdr = ( v.x[ i+1 ][ j ][ k ] - v.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dwdr = ( w.x[ i+1 ][ j ][ k ] - w.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dtdr = ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dpdr = ( p.x[ i+1 ][ j ][ k ] - p.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dcdr = ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );

					dudthe = ( u.x[ i ][ j+1 ][ k ] - u.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dvdthe = ( v.x[ i ][ j+1 ][ k ] - v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dwdthe = ( w.x[ i ][ j+1 ][ k ] - w.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dtdthe = ( t.x[ i ][ j+1 ][ k ] - t.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dpdthe = ( p.x[ i ][ j+1 ][ k ] - p.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dcdthe = ( c.x[ i ][ j+1 ][ k ] - c.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );


					dudphi = ( u.x[ i ][ j ][ k+1 ] - u.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dvdphi = ( v.x[ i ][ j ][ k+1 ] - v.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dwdphi = ( w.x[ i ][ j ][ k+1 ] - w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dtdphi = ( t.x[ i ][ j ][ k+1 ] - t.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dpdphi = ( p.x[ i ][ j ][ k+1 ] - p.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dcdphi = ( c.x[ i ][ j ][ k+1 ] - c.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );

				// 2nd order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components

					d2udr2 = ( u.x[ i+1 ][ j ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] ) / dr2;
					d2vdr2 = ( v.x[ i+1 ][ j ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] ) / dr2;
					d2wdr2 = ( w.x[ i+1 ][ j ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] ) / dr2;
					d2tdr2 = ( t.x[ i+1 ][ j ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i-1 ][ j ][ k ] ) / dr2;
					d2cdr2 = ( c.x[ i+1 ][ j ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i-1 ][ j ][ k ] ) / dr2;

					d2udthe2 = ( u.x[ i ][ j+1 ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2vdthe2 = ( v.x[ i ][ j+1 ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2wdthe2 = ( w.x[ i ][ j+1 ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2tdthe2 = ( t.x[ i ][ j+1 ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2cdthe2 = ( c.x[ i ][ j+1 ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j-1 ][ k ] ) / dthe2;

					d2udphi2 = ( u.x[ i ][ j ][ k+1 ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2vdphi2 = ( v.x[ i ][ j ][ k+1 ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2wdphi2 = ( w.x[ i ][ j ][ k+1 ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2tdphi2 = ( t.x[ i ][ j ][ k+1 ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2cdphi2 = ( c.x[ i ][ j ][ k+1 ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j ][ k-1 ] ) / dphi2;



					LHS_u = dudr + 2. * mue / dr2 - 2. * mue / rm2 + 2. * mue * ( rm2 * dthe2 )
									+ 2. * mue / ( rm2sinthe2 * dphi2 );

					RHS_u =  mue / dr2 * ( u.x[ i+1 ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] ) + mue / ( rm2 * dthe2 ) * ( u.x[ i ][ j + 1 ][ k ] + u.x[ i ][ j - 1 ][ k ] ) 
									+ mue / ( rm2sinthe2 * dphi2 ) * ( u.x[ i ][ j ][ k + 1 ] + u.x[ i ][ j ][ k - 1 ] ) 
									+ 4. * mue / rm * dudr 
									- dpdr + ( mue * costhe / rm2sinthe - v.x[ i ][ j ][ k ] / rm ) * dudthe
									- w.x[ i ][ j ][ k ] / rmsinthe * dudphi;

					aux_u.x[ i ][ j ][ k ] = RHS_u/ LHS_u;


					LHS_v = dvdthe / rm + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe2 ) / rm 
									+ 2. * mue / ( rm2sinthe2 * dphi2 ) + 2.* mue / dr2;

					RHS_v = mue / ( rm2 * dthe2 ) * ( v.x[ i ][ j + 1 ][ k ] + v.x[ i ][ j - 1 ][ k ] ) 
									+ mue / ( rm2sinthe2 * dphi2 ) * ( v.x[ i ][ j ][ k + 1 ] + v.x[ i ][ j ][ k - 1 ] ) 
									+ mue * ( v.x[ i+1 ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] ) / dr2 - w.x[ i ][ j ][ k ] / rmsinthe * dvdphi 
									- dpdthe / rm + 2. * mue / rm * dvdr 
									+ mue * costhe / rm2sinthe * dvdthe - 2. * mue * costhe / rm2sinthe2 * dwdphi
									- u.x[ i ][ j ][ k ] * dvdr + 2. * mue * dudthe;

					aux_v.x[ i ][ j ][ k ] = RHS_v / LHS_v;


					LHS_w = dwdphi / rmsinthe + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe ) / rm
									+ 2. * mue / ( rm2sinthe2 * dphi2 ) + 2.* mue / dr2;

					RHS_w = mue / ( rm2 * dthe2 ) * ( w.x[ i ][ j + 1 ][ k ] + w.x[ i ][ j - 1 ][ k ] ) 
									+ mue / ( rm2sinthe2 * dphi2 ) * ( w.x[ i ][ j ][ k + 1 ] + w.x[ i ][ j ][ k - 1 ] ) 
									+ mue * ( w.x[ i+1 ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] ) / dr2 - v.x[ i ][ j ][ k ] / rm * dwdthe 
									- dpdphi / rmsinthe + 2. * mue / rm * dwdr 
									+ mue * costhe / rm2sinthe * dwdthe + 2. * mue * costhe / rm2sinthe2 * dvdphi
									- u.x[ i ][ j ][ k ] * dwdr + 2. * mue / rmsinthe * dudphi;

					aux_w.x[ i ][ j ][ k ] = RHS_w / LHS_w;

				}
				else
				{
					dudr = 0.;
					dvdr = 0.;
					dwdr = 0.;
					dtdr = 0.;
					dpdr = 0.;
					dcdr = 0.;

					d2udr2 = 0.;
					d2vdr2 = 0.;
					d2wdr2 = 0.;
					d2tdr2 = 0.;
					d2cdr2 = 0.;

					dudthe = 0.;
					dvdthe = 0.;
					dwdthe = 0.;
					dtdthe = 0.;
					dpdthe = 0.;
					dcdthe = 0.;

					dudthe = 0.;
					d2udthe2 = 0.;
					d2vdthe2 = 0.;
					d2wdthe2 = 0.;
					d2tdthe2 = 0.;
					d2cdthe2 = 0.;

					dudphi = 0.;
					dvdphi = 0.;
					dwdphi = 0.;
					dtdphi = 0.;
					dpdphi = 0.;
					dcdphi = 0.;

					d2udphi2 = 0.;
					d2vdphi2 = 0.;
					d2wdphi2 = 0.;
					d2tdphi2 = 0.;
					d2cdphi2 = 0.;
				}
				
			}
		}
	}



	for ( int i = 1; i < im-1; i++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
			u.x[ i ][ j ][ k ] = aux_u.x[ i ][ j ][ k ];
			v.x[ i ][ j ][ k ] = aux_v.x[ i ][ j ][ k ];
			w.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k ];

			}
		}
	}


}
