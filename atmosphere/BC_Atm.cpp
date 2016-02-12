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




void BC_Atmosphere::BC_radius ( double tao, double tau, double pa, double ca, double co2a, double dr, Array_1D &rad, double co2_vegetation, double co2_ocean, double co2_land, Array_2D &Vegetation, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &co2, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_t, Array &rhs_c, Array &rhs_co2, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Rain, Array &Ice )
{
// boundary conditions for the r-direction, loop index i
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
//			if ( h.x[ 0 ][ j ][ k ] == 0. )																					// sea surface
			{
				u.x[ 0 ][ j ][ k ] = 0.;
//				v.x[ 0 ][ j ][ k ] = 0.;
//				w.x[ 0 ][ j ][ k ] = 0.;
//				v.x[ 0 ][ j ][ k ] = 2. * dr / rad.z[ 1 ] * v.x[ 1 ][ j ][ k ] + v.x[ 2 ][ j ][ k ];				// sea surface
//				w.x[ 0 ][ j ][ k ] = 2. * dr / rad.z[ 1 ] * w.x[ 1 ][ j ][ k ] + w.x[ 2 ][ j ][ k ];			// sea surface
//				v.x[ 0 ][ j ][ k ] = c43 * v.x[ 1 ][ j ][ k ] - c13 * v.x[ 2 ][ j ][ k ];
//				w.x[ 0 ][ j ][ k ] = c43 * w.x[ 1 ][ j ][ k ] - c13 * w.x[ 2 ][ j ][ k ];
				v.x[ 0 ][ j ][ k ] = v.x[ 3 ][ j ][ k ] - 3. * v.x[ 2 ][ j ][ k ] + 3. * v.x[ 1 ][ j ][ k ];				// extrapolation
				w.x[ 0 ][ j ][ k ] = w.x[ 3 ][ j ][ k ] - 3. * w.x[ 2 ][ j ][ k ] + 3. * w.x[ 1 ][ j ][ k ];			// leads to round harmonic profiles
//				p.x[ im-1 ][ j ][ k ] = pa;

				u.x[ im-1 ][ j ][ k ] = 0.;
				v.x[ im-1 ][ j ][ k ] = 0.;																					// stratosphere
				w.x[ im-1 ][ j ][ k ] = 0.;
				t.x[ im-1 ][ j ][ k ] = tao;

				p.x[ im-1 ][ j ][ k ] = c43 * p.x[ im-2 ][ j ][ k ] - c13 * p.x[ im-3 ][ j ][ k ];			// stratosphere
//				p.x[ 0 ][ j ][ k ] = c43 * p.x[ 1 ][ j ][ k ] - c13 * p.x[ 2 ][ j ][ k ];							// sea surface

				c.x[ im-1 ][ j ][ k ] = 0.;
				co2.x[ im-1 ][ j ][ k ] = c43 * co2.x[ im-2 ][ j ][ k ] - c13 * co2.x[ im-3 ][ j ][ k ];			// stratosphere
/*
				rhs_u.x[ 0 ][ j ][ k ] = c43 * rhs_u.x[ 1 ][ j ][ k ] - c13 * rhs_u.x[ 2 ][ j ][ k ];
				rhs_v.x[ 0 ][ j ][ k ] = c43 * rhs_v.x[ 1 ][ j ][ k ] - c13 * rhs_v.x[ 2 ][ j ][ k ];
				rhs_w.x[ 0 ][ j ][ k ] = c43 * rhs_w.x[ 1 ][ j ][ k ] - c13 * rhs_w.x[ 2 ][ j ][ k ];
				rhs_t.x[ 0 ][ j ][ k ] = c43 * rhs_t.x[ 1 ][ j ][ k ] - c13 * rhs_t.x[ 2 ][ j ][ k ];
				rhs_c.x[ 0 ][ j ][ k ] = c43 * rhs_c.x[ 1 ][ j ][ k ] - c13 * rhs_c.x[ 2 ][ j ][ k ];
				rhs_co2.x[ 0 ][ j ][ k ] = c43 * rhs_co2.x[ 1 ][ j ][ k ] - c13 * rhs_co2.x[ 2 ][ j ][ k ];
*/
/*
				rhs_u.x[ 0 ][ j ][ k ] = 0.;
				rhs_v.x[ 0 ][ j ][ k ] = 0.;
				rhs_w.x[ 0 ][ j ][ k ] = 0.;
				rhs_t.x[ 0 ][ j ][ k ] = 0.;
				rhs_c.x[ 0 ][ j ][ k ] = 0.;
				rhs_co2.x[ 0 ][ j ][ k ] = 0.;

				rhs_u.x[ im-1 ][ j ][ k ] = c43 * rhs_u.x[ im-2 ][ j ][ k ] - c13 * rhs_u.x[ im-3 ][ j ][ k ];			// stratosphere
				rhs_v.x[ im-1 ][ j ][ k ] = c43 * rhs_v.x[ im-2 ][ j ][ k ] - c13 * rhs_v.x[ im-3 ][ j ][ k ];			// stratosphere
				rhs_w.x[ im-1 ][ j ][ k ] = c43 * rhs_w.x[ im-2 ][ j ][ k ] - c13 * rhs_w.x[ im-3 ][ j ][ k ];		// stratosphere
				rhs_t.x[ im-1 ][ j ][ k ] = c43 * rhs_t.x[ im-2 ][ j ][ k ] - c13 * rhs_t.x[ im-3 ][ j ][ k ];			// stratosphere
				rhs_c.x[ im-1 ][ j ][ k ] = c43 * rhs_c.x[ im-2 ][ j ][ k ] - c13 * rhs_c.x[ im-3 ][ j ][ k ];
				rhs_co2.x[ im-1 ][ j ][ k ] = c43 * rhs_co2.x[ im-2 ][ j ][ k ] - c13 * rhs_co2.x[ im-3 ][ j ][ k ];
*/
/*
				aux_u.x[ 0 ][ j ][ k ] = c43 * aux_u.x[ 1 ][ j ][ k ] - c13 * aux_u.x[ 2 ][ j ][ k ];
				aux_v.x[ 0 ][ j ][ k ] = c43 * aux_v.x[ 1 ][ j ][ k ] - c13 * aux_v.x[ 2 ][ j ][ k ];
				aux_w.x[ 0 ][ j ][ k ] = c43 * aux_w.x[ 1 ][ j ][ k ] - c13 * aux_w.x[ 2 ][ j ][ k ];
*/
/*
				aux_u.x[ 0 ][ j ][ k ] = 0.;
				aux_v.x[ 0 ][ j ][ k ] = 0.;
				aux_w.x[ 0 ][ j ][ k ] = 0.;

				aux_u.x[ im-1 ][ j ][ k ] = c43 * aux_u.x[ im-2 ][ j ][ k ] - c13 * aux_u.x[ im-3 ][ j ][ k ];			// stratosphere
				aux_v.x[ im-1 ][ j ][ k ] = c43 * aux_v.x[ im-2 ][ j ][ k ] - c13 * aux_v.x[ im-3 ][ j ][ k ];			// stratosphere
				aux_w.x[ im-1 ][ j ][ k ] = c43 * aux_w.x[ im-2 ][ j ][ k ] - c13 * aux_w.x[ im-3 ][ j ][ k ];			// stratosphere

				Latency.x[ 0 ][ j ][ k ] = c43 * Latency.x[ 1 ][ j ][ k ] - c13 * Latency.x[ 2 ][ j ][ k ];
				Rain.x[ 0 ][ j ][ k ] = c43 * Rain.x[ 1 ][ j ][ k ] - c13 * Rain.x[ 2 ][ j ][ k ];
				Ice.x[ 0 ][ j ][ k ] = c43 * Ice.x[ 1 ][ j ][ k ] - c13 * Ice.x[ 2 ][ j ][ k ];
*/
/*
				Latency.x[ im-1 ][ j ][ k ] = c43 * Latency.x[ im-2 ][ j ][ k ] - c13 * Latency.x[ im-3 ][ j ][ k ];			// stratosphere
				Rain.x[ im-1 ][ j ][ k ] = c43 * Rain.x[ im-2 ][ j ][ k ] - c13 * Rain.x[ im-3 ][ j ][ k ];			// stratosphere
				Ice.x[ im-1 ][ j ][ k ] = c43 * Ice.x[ im-2 ][ j ][ k ] - c13 * Ice.x[ im-3 ][ j ][ k ];								// stratosphere
*/
/*
				Latency.x[ im-1 ][ j ][ k ] = 0.;										// stratosphere
				Rain.x[ im-1 ][ j ][ k ] = 0.;											// stratosphere
				Ice.x[ im-1 ][ j ][ k ] = 0.;												// stratosphere
*/
			}
/*
			else
			{
//  continental surface
				u.x[ 0 ][ j ][ k ] = 0.;
				v.x[ 0 ][ j ][ k ] = 0.;
				w.x[ 0 ][ j ][ k ] = 0.;

				u.x[ im-1 ][ j ][ k ] = 0.;																					// stratosphere
				v.x[ im-1 ][ j ][ k ] = 0.;
				w.x[ im-1 ][ j ][ k ] = 0.;

				t.x[ im-1 ][ j ][ k ] = tao;
				c.x[ im-1 ][ j ][ k ] = 0.;
				p.x[ im-1 ][ j ][ k ] = pa;

				p.x[ im-1 ][ j ][ k ] = c43 * p.x[ im-2 ][ j ][ k ] - c13 * p.x[ im-3 ][ j ][ k ];			// stratosphere
				p.x[ 0 ][ j ][ k ] = c43 * p.x[ 1 ][ j ][ k ] - c13 * p.x[ 2 ][ j ][ k ];							// continent surface

				c.x[ im-1 ][ j ][ k ] = 0.;
				co2.x[ im-1 ][ j ][ k ] = 1.;

				rhs_u.x[ 0 ][ j ][ k ] = 0.;
				rhs_v.x[ 0 ][ j ][ k ] = 0.;
				rhs_w.x[ 0 ][ j ][ k ] = 0.;
				rhs_t.x[ 0 ][ j ][ k ] = 0.;
				rhs_c.x[ 0 ][ j ][ k ] = 0.;
				rhs_co2.x[ 0 ][ j ][ k ] = 0.;

				aux_u.x[ 0 ][ j ][ k ] = 0.;
				aux_v.x[ 0 ][ j ][ k ] = 0.;
				aux_w.x[ 0 ][ j ][ k ] = 0.;

				Latency.x[ 0 ][ j ][ k ] = 0.;
				Rain.x[ 0 ][ j ][ k ] = 0.;
				Ice.x[ 0 ][ j ][ k ] = 0.;
			}
*/
		}
	}

}






void BC_Atmosphere::BC_theta ( Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &co2, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_t, Array &rhs_c, Array &rhs_co2, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Rain, Array &Ice )
{
// boundary conditions for the the-direction, loop index j

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
// zero tangent ( von Neumann condition ) or constant value ( Dirichlet condition )

			t.x[ i ][ 0 ][ k ] = c43 * t.x[ i ][ 1 ][ k ] - c13 * t.x[ i ][ 2 ][ k ];
			t.x[ i ][ jm-1 ][ k ] = c43 * t.x[ i ][ jm-2 ][ k ] - c13 * t.x[ i ][ jm-3 ][ k ];
/*
			u.x[ i ][ 0 ][ k ] = c43 * u.x[ i ][ 1 ][ k ] - c13 * u.x[ i ][ 2 ][ k ];
			u.x[ i ][ jm-1 ][ k ] = c43 * u.x[ i ][ jm-2 ][ k ] - c13 * u.x[ i ][ jm-3 ][ k ];

			v.x[ i ][ 0 ][ k ] = c43 * v.x[ i ][ 1 ][ k ] - c13 * v.x[ i ][ 2 ][ k ];
			v.x[ i ][ jm-1 ][ k ] = c43 * v.x[ i ][ jm-2 ][ k ] - c13 * v.x[ i ][ jm-3 ][ k ];

			w.x[ i ][ 0 ][ k ] = c43 * w.x[ i ][ 1 ][ k ] - c13 * w.x[ i ][ 2 ][ k ];
			w.x[ i ][ jm-1 ][ k ] = c43 * w.x[ i ][ jm-2 ][ k ] - c13 * w.x[ i ][ jm-3 ][ k ];
*/

			u.x[ i ][ 0 ][ k ] = 0.;
			u.x[ i ][ jm-1 ][ k ] = 0.;

			v.x[ i ][ 0 ][ k ] = 0.;
			v.x[ i ][ jm-1 ][ k ] = 0.;

			w.x[ i ][ 0 ][ k ] = 0.;
			w.x[ i ][ jm-1 ][ k ] = 0.;

			p.x[ i ][ 0 ][ k ] = c43 * p.x[ i ][ 1 ][ k ] - c13 * p.x[ i ][ 2 ][ k ];
			p.x[ i ][ jm-1 ][ k ] = c43 * p.x[ i ][ jm-2 ][ k ] - c13 * p.x[ i ][ jm-3 ][ k ];


			c.x[ i ][ 0 ][ k ] = c43 * c.x[ i ][ 1 ][ k ] - c13 * c.x[ i ][ 2 ][ k ];
			c.x[ i ][ jm-1 ][ k ] = c43 * c.x[ i ][ jm-2 ][ k ] - c13 * c.x[ i ][ jm-3 ][ k ];

			co2.x[ i ][ 0 ][ k ] = c43 * co2.x[ i ][ 1 ][ k ] - c13 * co2.x[ i ][ 2 ][ k ];
			co2.x[ i ][ jm-1 ][ k ] = c43 * co2.x[ i ][ jm-2 ][ k ] - c13 * co2.x[ i ][ jm-3 ][ k ];

			rhs_u.x[ i ][ 0 ][ k ] = c43 * rhs_u.x[ i ][ 1 ][ k ] - c13 * rhs_u.x[ i ][ 2 ][ k ];
			rhs_u.x[ i ][ jm-1 ][ k ] = c43 * rhs_u.x[ i ][ jm-2 ][ k ] - c13 * rhs_u.x[ i ][ jm-3 ][ k ];

			rhs_v.x[ i ][ 0 ][ k ] = c43 * rhs_v.x[ i ][ 1 ][ k ] - c13 * rhs_v.x[ i ][ 2 ][ k ];
			rhs_v.x[ i ][ jm-1 ][ k ] = c43 * rhs_v.x[ i ][ jm-2 ][ k ] - c13 * rhs_v.x[ i ][ jm-3 ][ k ];

			rhs_w.x[ i ][ 0 ][ k ] = c43 * rhs_w.x[ i ][ 1 ][ k ] - c13 * rhs_w.x[ i ][ 2 ][ k ];
			rhs_w.x[ i ][ jm-1 ][ k ] = c43 * rhs_w.x[ i ][ jm-2 ][ k ] - c13 * rhs_w.x[ i ][ jm-3 ][ k ];

			aux_u.x[ i ][ 0 ][ k ] = c43 * aux_u.x[ i ][ 1 ][ k ] - c13 * aux_u.x[ i ][ 2 ][ k ];
			aux_u.x[ i ][ jm-1 ][ k ] = c43 * aux_u.x[ i ][ jm-2 ][ k ] - c13 * aux_u.x[ i ][ jm-3 ][ k ];

			aux_v.x[ i ][ 0 ][ k ] = c43 * aux_v.x[ i ][ 1 ][ k ] - c13 * aux_v.x[ i ][ 2 ][ k ];
			aux_v.x[ i ][ jm-1 ][ k ] = c43 * aux_v.x[ i ][ jm-2 ][ k ] - c13 * aux_v.x[ i ][ jm-3 ][ k ];

			aux_w.x[ i ][ 0 ][ k ] = c43 * aux_w.x[ i ][ 1 ][ k ] - c13 * aux_w.x[ i ][ 2 ][ k ];
			aux_w.x[ i ][ jm-1 ][ k ] = c43 * aux_w.x[ i ][ jm-2 ][ k ] - c13 * aux_w.x[ i ][ jm-3 ][ k ];

			rhs_t.x[ i ][ 0 ][ k ] = c43 * rhs_t.x[ i ][ 1 ][ k ] - c13 * rhs_t.x[ i ][ 2 ][ k ];
			rhs_t.x[ i ][ jm-1 ][ k ] = c43 * rhs_t.x[ i ][ jm-2 ][ k ] - c13 * rhs_t.x[ i ][ jm-3 ][ k ];

			rhs_c.x[ i ][ 0 ][ k ] = c43 * rhs_c.x[ i ][ 1 ][ k ] - c13 * rhs_c.x[ i ][ 2 ][ k ];
			rhs_c.x[ i ][ jm-1 ][ k ] = c43 * rhs_c.x[ i ][ jm-2 ][ k ] - c13 * rhs_c.x[ i ][ jm-3 ][ k ];

			rhs_co2.x[ i ][ 0 ][ k ] = c43 * rhs_co2.x[ i ][ 1 ][ k ] - c13 * rhs_co2.x[ i ][ 2 ][ k ];
			rhs_co2.x[ i ][ jm-1 ][ k ] = c43 * rhs_co2.x[ i ][ jm-2 ][ k ] - c13 * rhs_co2.x[ i ][ jm-3 ][ k ];

			Latency.x[ i ][ 0 ][ k ] = c43 * Latency.x[ i ][ 1 ][ k ] - c13 * Latency.x[ i ][ 2 ][ k ];
			Latency.x[ i ][ jm-1 ][ k ] = c43 * Latency.x[ i ][ jm-2 ][ k ] - c13 * Latency.x[ i ][ jm-3 ][ k ];

			Rain.x[ i ][ 0 ][ k ] = c43 * Rain.x[ i ][ 1 ][ k ] - c13 * Rain.x[ i ][ 2 ][ k ];
			Rain.x[ i ][ jm-1 ][ k ] = c43 * Rain.x[ i ][ jm-2 ][ k ] - c13 * Rain.x[ i ][ jm-3 ][ k ];

			Ice.x[ i ][ 0 ][ k ] = c43 * Ice.x[ i ][ 1 ][ k ] - c13 * Ice.x[ i ][ 2 ][ k ];
			Ice.x[ i ][ jm-1 ][ k ] = c43 * Ice.x[ i ][ jm-2 ][ k ] - c13 * Ice.x[ i ][ jm-3 ][ k ];

		}
	}
}







void BC_Atmosphere::BC_phi ( Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &co2, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_t, Array &rhs_c, Array &rhs_co2, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Rain, Array &Ice )
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

			co2.x[ i ][ j ][ 0 ] = c43 * co2.x[ i ][ j ][ 1 ] - c13 * co2.x[ i ][ j ][ 2 ];
			co2.x[ i ][ j ][ km-1 ] = c43 * co2.x[ i ][ j ][ km-2 ] - c13 * co2.x[ i ][ j ][ km-3 ];
			co2.x[ i ][ j ][ 0 ] = co2.x[ i ][ j ][ km-1 ] = ( co2.x[ i ][ j ][ 0 ] + co2.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_u.x[ i ][ j ][ 0 ] = c43 * rhs_u.x[ i ][ j ][ 1 ] - c13 * rhs_u.x[ i ][ j ][ 2 ];
			rhs_u.x[ i ][ j ][ km-1 ] = c43 * rhs_u.x[ i ][ j ][ km-2 ] - c13 * rhs_u.x[ i ][ j ][ km-3 ];
			rhs_u.x[ i ][ j ][ 0 ] = rhs_u.x[ i ][ j ][ km-1 ] = ( rhs_u.x[ i ][ j ][ 0 ] + rhs_u.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_t.x[ i ][ j ][ 0 ] = c43 * rhs_t.x[ i ][ j ][ 1 ] - c13 * rhs_t.x[ i ][ j ][ 2 ];
			rhs_t.x[ i ][ j ][ km-1 ] = c43 * rhs_t.x[ i ][ j ][ km-2 ] - c13 * rhs_t.x[ i ][ j ][ km-3 ];
			rhs_t.x[ i ][ j ][ 0 ] = rhs_t.x[ i ][ j ][ km-1 ] = ( rhs_t.x[ i ][ j ][ 0 ] + rhs_t.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_v.x[ i ][ j ][ 0 ] = c43 * rhs_v.x[ i ][ j ][ 1 ] - c13 * rhs_v.x[ i ][ j ][ 2 ];
			rhs_v.x[ i ][ j ][ km-1 ] = c43 * rhs_v.x[ i ][ j ][ km-2 ] - c13 * rhs_v.x[ i ][ j ][ km-3 ];
			rhs_v.x[ i ][ j ][ 0 ] = rhs_v.x[ i ][ j ][ km-1 ] = ( rhs_v.x[ i ][ j ][ 0 ] + rhs_v.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_w.x[ i ][ j ][ 0 ] = c43 * rhs_w.x[ i ][ j ][ 1 ] - c13 * rhs_w.x[ i ][ j ][ 2 ];
			rhs_w.x[ i ][ j ][ km-1 ] = c43 * rhs_w.x[ i ][ j ][ km-2 ] - c13 * rhs_w.x[ i ][ j ][ km-3 ];
			rhs_w.x[ i ][ j ][ 0 ] = rhs_w.x[ i ][ j ][ km-1 ] = ( rhs_w.x[ i ][ j ][ 0 ] + rhs_w.x[ i ][ j ][ km-1 ] ) / 2.;

			aux_u.x[ i ][ j ][ 0 ] = c43 * aux_u.x[ i ][ j ][ 1 ] - c13 * aux_u.x[ i ][ j ][ 2 ];
			aux_u.x[ i ][ j ][ km-1 ] = c43 * aux_u.x[ i ][ j ][ km-2 ] - c13 * aux_u.x[ i ][ j ][ km-3 ];
			aux_u.x[ i ][ j ][ 0 ] = aux_u.x[ i ][ j ][ km-1 ] = ( aux_u.x[ i ][ j ][ 0 ] + aux_u.x[ i ][ j ][ km-1 ] ) / 2.;

			aux_v.x[ i ][ j ][ 0 ] = c43 * aux_v.x[ i ][ j ][ 1 ] - c13 * aux_v.x[ i ][ j ][ 2 ];
			aux_v.x[ i ][ j ][ km-1 ] = c43 * aux_v.x[ i ][ j ][ km-2 ] - c13 * aux_v.x[ i ][ j ][ km-3 ];
			aux_v.x[ i ][ j ][ 0 ] = aux_v.x[ i ][ j ][ km-1 ] = ( aux_v.x[ i ][ j ][ 0 ] + aux_v.x[ i ][ j ][ km-1 ] ) / 2.;

			aux_w.x[ i ][ j ][ 0 ] = c43 * aux_w.x[ i ][ j ][ 1 ] - c13 * aux_w.x[ i ][ j ][ 2 ];
			aux_w.x[ i ][ j ][ km-1 ] = c43 * aux_w.x[ i ][ j ][ km-2 ] - c13 * aux_w.x[ i ][ j ][ km-3 ];
			aux_w.x[ i ][ j ][ 0 ] = aux_w.x[ i ][ j ][ km-1 ] = ( aux_w.x[ i ][ j ][ 0 ] + aux_w.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_t.x[ i ][ j ][ 0 ] = c43 * rhs_t.x[ i ][ j ][ 1 ] - c13 * rhs_t.x[ i ][ j ][ 2 ];
			rhs_t.x[ i ][ j ][ km-1 ] = c43 * rhs_t.x[ i ][ j ][ km-2 ] - c13 * rhs_t.x[ i ][ j ][ km-3 ];
			rhs_t.x[ i ][ j ][ 0 ] = rhs_t.x[ i ][ j ][ km-1 ] = ( rhs_t.x[ i ][ j ][ 0 ] + rhs_t.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_c.x[ i ][ j ][ 0 ] = c43 * rhs_c.x[ i ][ j ][ 1 ] - c13 * rhs_c.x[ i ][ j ][ 2 ];
			rhs_c.x[ i ][ j ][ km-1 ] = c43 * rhs_c.x[ i ][ j ][ km-2 ] - c13 * rhs_c.x[ i ][ j ][ km-3 ];
			rhs_c.x[ i ][ j ][ 0 ] = rhs_c.x[ i ][ j ][ km-1 ] = ( rhs_c.x[ i ][ j ][ 0 ] + rhs_c.x[ i ][ j ][ km-1 ] ) / 2.;

			rhs_co2.x[ i ][ j ][ 0 ] = c43 * rhs_co2.x[ i ][ j ][ 1 ] - c13 * rhs_co2.x[ i ][ j ][ 2 ];
			rhs_co2.x[ i ][ j ][ km-1 ] = c43 * rhs_co2.x[ i ][ j ][ km-2 ] - c13 * rhs_co2.x[ i ][ j ][ km-3 ];
			rhs_co2.x[ i ][ j ][ 0 ] = rhs_co2.x[ i ][ j ][ km-1 ] = ( rhs_co2.x[ i ][ j ][ 0 ] + rhs_co2.x[ i ][ j ][ km-1 ] ) / 2.;

			Latency.x[ i ][ j ][ 0 ] = c43 * Latency.x[ i ][ j ][ 1 ] - c13 * Latency.x[ i ][ j ][ 2 ];
			Latency.x[ i ][ j ][ km-1 ] = c43 * Latency.x[ i ][ j ][ km-2 ] - c13 * Latency.x[ i ][ j ][ km-3 ];
			Latency.x[ i ][ j ][ 0 ] = Latency.x[ i ][ j ][ km-1 ] = ( Latency.x[ i ][ j ][ 0 ] + Latency.x[ i ][ j ][ km-1 ] ) / 2.;

			Rain.x[ i ][ j ][ 0 ] = c43 * Rain.x[ i ][ j ][ 1 ] - c13 * Rain.x[ i ][ j ][ 2 ];
			Rain.x[ i ][ j ][ km-1 ] = c43 * Rain.x[ i ][ j ][ km-2 ] - c13 * Rain.x[ i ][ j ][ km-3 ];
			Rain.x[ i ][ j ][ 0 ] = Rain.x[ i ][ j ][ km-1 ] = ( Rain.x[ i ][ j ][ 0 ] + Rain.x[ i ][ j ][ km-1 ] ) / 2.;

			Ice.x[ i ][ j ][ 0 ] = c43 * Ice.x[ i ][ j ][ 1 ] - c13 * Ice.x[ i ][ j ][ 2 ];
			Ice.x[ i ][ j ][ km-1 ] = c43 * Ice.x[ i ][ j ][ km-2 ] - c13 * Ice.x[ i ][ j ][ km-3 ];
			Ice.x[ i ][ j ][ 0 ] = Ice.x[ i ][ j ][ km-1 ] = ( Ice.x[ i ][ j ][ 0 ] + Ice.x[ i ][ j ][ km-1 ] ) / 2.;

		}
	}
}



void BC_Atmosphere::BC_NST_control_2D ( double dr, double dthe, double dphi, double re, double mue_air, double mue_water, Array &h, Array &v, Array &w, Array &p, Array_2D &aux_2D_v, Array_2D &aux_2D_w, Array_1D &rad, Array_1D &the )
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
//			if ( h.x[ 0 ][ j ][ k ] == 0. ) 
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
								- .5 * dpdthe / rm + 2. * mue / rm * ( 4. * v.x[ 1 ][ j ][ k ] - v.x[ 2 ][ j ][ k ] ) / ( 2. * dr ) 
								+ mue * costhe / rm2sinthe * dvdthe - 2. * mue * costhe / rm2sinthe2 * dwdphi;

				aux_2D_v.y[ j ][ k ] = RHS_v / LHS_v;

/*
				LHS_v = dvdthe / rm + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe2 ) / rm 
								+ 2. * mue / ( rm2sinthe2 * dphi2 );

				RHS_v = mue / ( rm2 * dthe2 ) * ( v.x[ 0 ][ j + 1 ][ k ] + v.x[ 0 ][ j - 1 ][ k ] ) 
								+ mue / ( rm2sinthe2 * dphi2 ) * ( v.x[ 0 ][ j ][ k + 1 ] + v.x[ 0 ][ j ][ k - 1 ] ) 
								- w.x[ 0 ][ j ][ k ] / rmsinthe * dvdphi 
								- .5 * dpdthe / rm 
								+ mue * costhe / rm2sinthe * dvdthe - 2. * mue * costhe / rm2sinthe2 * dwdphi;

				aux_2D_v.y[ j ][ k ] = RHS_v / LHS_v;
*/

				LHS_w = dwdphi / rmsinthe + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe ) / rm
								+ 2. * mue / ( rm2sinthe2 * dphi2 ) + 6. * mue / rm / ( 2. * dr ) - mue / dr2;

				RHS_w = mue / ( rm2 * dthe2 ) * ( w.x[ 0 ][ j + 1 ][ k ] + w.x[ 0 ][ j - 1 ][ k ] ) 
								+ mue / ( rm2sinthe2 * dphi2 ) * ( w.x[ 0 ][ j ][ k + 1 ] + w.x[ 0 ][ j ][ k - 1 ] ) 
								+ mue * ( - 2. * w.x[ 1 ][ j ][ k ] + w.x[ 2 ][ j ][ k ] ) / dr2 - v.x[ 0 ][ j ][ k ] / rm * dwdthe 
								- .5 * dpdphi / rmsinthe + 2. * mue / rm * ( 4. * w.x[ 1 ][ j ][ k ] - w.x[ 2 ][ j ][ k ] ) / ( 2. * dr ) 
								+ mue * costhe / rm2sinthe * dwdthe + 2. * mue * costhe / rm2sinthe2 * dvdphi;

				aux_2D_w.y[ j ][ k ] = RHS_w / LHS_w;

 /*
				LHS_w = dwdphi / rmsinthe + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe ) / rm
								+ 2. * mue / ( rm2sinthe2 * dphi2 );

				RHS_w = mue / ( rm2 * dthe2 ) * ( w.x[ 0 ][ j + 1 ][ k ] + w.x[ 0 ][ j - 1 ][ k ] ) 
								+ mue / ( rm2sinthe2 * dphi2 ) * ( w.x[ 0 ][ j ][ k + 1 ] + w.x[ 0 ][ j ][ k - 1 ] ) 
								- v.x[ 0 ][ j ][ k ] / rm * dwdthe 
								- .5 * dpdphi / rmsinthe
								+ mue * costhe / rm2sinthe * dwdthe + 2. * mue * costhe / rm2sinthe2 * dvdphi;

				aux_2D_w.y[ j ][ k ] = RHS_w / LHS_w;
*/
			}
//			else
			{
//				aux_2D_v.y[ j ][ k ] = 0.;
//				aux_2D_w.y[ j ][ k ] = 0.;
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






void BC_Atmosphere::BC_NST_control_3D ( double dr, double dthe, double dphi, double re, double mue_air, double mue_water, Array &h, Array &u, Array &v, Array &w, Array &t, Array &p, Array &c, Array &co2, Array &aux_u, Array &aux_v, Array &aux_w, Array_1D &rad, Array_1D &the )
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
//				if ( h.x[ i ][ j ][ k ] == 0. ) 
				{
// 1st order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
					dudr = ( u.x[ i+1 ][ j ][ k ] - u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dvdr = ( v.x[ i+1 ][ j ][ k ] - v.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dwdr = ( w.x[ i+1 ][ j ][ k ] - w.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dtdr = ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dpdr = ( p.x[ i+1 ][ j ][ k ] - p.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dcdr = ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
					dco2dr = ( co2.x[ i+1 ][ j ][ k ] - co2.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );

					dudthe = ( u.x[ i ][ j+1 ][ k ] - u.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dvdthe = ( v.x[ i ][ j+1 ][ k ] - v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dwdthe = ( w.x[ i ][ j+1 ][ k ] - w.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dtdthe = ( t.x[ i ][ j+1 ][ k ] - t.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dpdthe = ( p.x[ i ][ j+1 ][ k ] - p.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dcdthe = ( c.x[ i ][ j+1 ][ k ] - c.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
					dco2dthe = ( co2.x[ i ][ j+1 ][ k ] - co2.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );


					dudphi = ( u.x[ i ][ j ][ k+1 ] - u.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dvdphi = ( v.x[ i ][ j ][ k+1 ] - v.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dwdphi = ( w.x[ i ][ j ][ k+1 ] - w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dtdphi = ( t.x[ i ][ j ][ k+1 ] - t.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dpdphi = ( p.x[ i ][ j ][ k+1 ] - p.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dcdphi = ( c.x[ i ][ j ][ k+1 ] - c.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
					dco2dphi = ( co2.x[ i ][ j ][ k+1 ] - co2.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );

				// 2nd order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components

					d2udr2 = ( u.x[ i+1 ][ j ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] ) / dr2;
					d2vdr2 = ( v.x[ i+1 ][ j ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] ) / dr2;
					d2wdr2 = ( w.x[ i+1 ][ j ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] ) / dr2;
					d2tdr2 = ( t.x[ i+1 ][ j ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i-1 ][ j ][ k ] ) / dr2;
					d2cdr2 = ( c.x[ i+1 ][ j ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i-1 ][ j ][ k ] ) / dr2;
					d2co2dr2 = ( co2.x[ i+1 ][ j ][ k ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i-1 ][ j ][ k ] ) / dr2;

					d2udthe2 = ( u.x[ i ][ j+1 ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2vdthe2 = ( v.x[ i ][ j+1 ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2wdthe2 = ( w.x[ i ][ j+1 ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2tdthe2 = ( t.x[ i ][ j+1 ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2cdthe2 = ( c.x[ i ][ j+1 ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j-1 ][ k ] ) / dthe2;
					d2co2dthe2 = ( co2.x[ i ][ j+1 ][ k ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i ][ j-1 ][ k ] ) / dthe2;

					d2udphi2 = ( u.x[ i ][ j ][ k+1 ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2vdphi2 = ( v.x[ i ][ j ][ k+1 ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2wdphi2 = ( w.x[ i ][ j ][ k+1 ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2tdphi2 = ( t.x[ i ][ j ][ k+1 ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2cdphi2 = ( c.x[ i ][ j ][ k+1 ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j ][ k-1 ] ) / dphi2;
					d2co2dphi2 = ( co2.x[ i ][ j ][ k+1 ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i ][ j ][ k-1 ] ) / dphi2;



					LHS_u = dudr + 2. * mue / dr2 - 2. * mue / rm2 + 2. * mue * ( rm2 * dthe2 )
									+ 2. * mue / ( rm2sinthe2 * dphi2 );

					RHS_u =  mue / dr2 * ( u.x[ i+1 ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] ) + mue / ( rm2 * dthe2 ) * ( u.x[ i ][ j + 1 ][ k ] + u.x[ i ][ j - 1 ][ k ] ) 
									+ mue / ( rm2sinthe2 * dphi2 ) * ( u.x[ i ][ j ][ k + 1 ] + u.x[ i ][ j ][ k - 1 ] ) 
									+ 4. * mue / rm * dudr 
									- .5 * dpdr + ( mue * costhe / rm2sinthe - v.x[ i ][ j ][ k ] / rm ) * dudthe
									- w.x[ i ][ j ][ k ] / rmsinthe * dudphi;

					aux_u.x[ i ][ j ][ k ] = RHS_u/ LHS_u;


					LHS_v = dvdthe / rm + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe2 ) / rm 
									+ 2. * mue / ( rm2sinthe2 * dphi2 ) + 2.* mue / dr2;

					RHS_v = mue / ( rm2 * dthe2 ) * ( v.x[ i ][ j + 1 ][ k ] + v.x[ i ][ j - 1 ][ k ] ) 
									+ mue / ( rm2sinthe2 * dphi2 ) * ( v.x[ i ][ j ][ k + 1 ] + v.x[ i ][ j ][ k - 1 ] ) 
									+ mue * ( v.x[ i+1 ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] ) / dr2 - w.x[ i ][ j ][ k ] / rmsinthe * dvdphi 
									- .5 * dpdthe / rm + 2. * mue / rm * dvdr 
									+ mue * costhe / rm2sinthe * dvdthe - 2. * mue * costhe / rm2sinthe2 * dwdphi
									- u.x[ i ][ j ][ k ] * dvdr + 2. * mue * dudthe;

					aux_v.x[ i ][ j ][ k ] = RHS_v / LHS_v;


					LHS_w = dwdphi / rmsinthe + 2. * mue / ( rm2 * dthe2 ) + mue * ( 1. + costhe * costhe / sinthe ) / rm
									+ 2. * mue / ( rm2sinthe2 * dphi2 ) + 2.* mue / dr2;

					RHS_w = mue / ( rm2 * dthe2 ) * ( w.x[ i ][ j + 1 ][ k ] + w.x[ i ][ j - 1 ][ k ] ) 
									+ mue / ( rm2sinthe2 * dphi2 ) * ( w.x[ i ][ j ][ k + 1 ] + w.x[ i ][ j ][ k - 1 ] ) 
									+ mue * ( w.x[ i+1 ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] ) / dr2 - v.x[ i ][ j ][ k ] / rm * dwdthe 
									- .5 * dpdphi / rmsinthe + 2. * mue / rm * dwdr 
									+ mue * costhe / rm2sinthe * dwdthe + 2. * mue * costhe / rm2sinthe2 * dvdphi
									- u.x[ i ][ j ][ k ] * dwdr + 2. * mue / rmsinthe * dudphi;

					aux_w.x[ i ][ j ][ k ] = RHS_w / LHS_w;

				}
/*
				else
				{
					dudr = 0.;
					dvdr = 0.;
					dwdr = 0.;
					dtdr = 0.;
					dpdr = 0.;
					dcdr = 0.;
					dco2dr = 0.;

					d2udr2 = 0.;
					d2vdr2 = 0.;
					d2wdr2 = 0.;
					d2tdr2 = 0.;
					d2cdr2 = 0.;
					d2co2dr2 = 0.;

					dudthe = 0.;
					dvdthe = 0.;
					dwdthe = 0.;
					dtdthe = 0.;
					dpdthe = 0.;
					dcdthe = 0.;
					dco2dthe = 0.;

					dudthe = 0.;
					d2udthe2 = 0.;
					d2vdthe2 = 0.;
					d2wdthe2 = 0.;
					d2tdthe2 = 0.;
					d2cdthe2 = 0.;
					d2co2dthe2 = 0.;

					dudphi = 0.;
					dvdphi = 0.;
					dwdphi = 0.;
					dtdphi = 0.;
					dpdphi = 0.;
					dcdphi = 0.;
					dco2dphi = 0.;

					d2udphi2 = 0.;
					d2vdphi2 = 0.;
					d2wdphi2 = 0.;
					d2tdphi2 = 0.;
					d2cdphi2 = 0.;
					d2co2dphi2 = 0.;
				}
*/
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
