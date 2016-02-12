/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to restore the old by new values inside the iterational processes
*/

#include <iostream>
#include "Restore_Hyd.h"

using namespace std;


Restore::Restore ( int im, int jm, int km )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
}


Restore::~Restore () {}


void Restore::restoreOldNew_3D ( double coeff, Array &u, Array &v, Array &w, Array &t, Array &p_dyn, Array &c, Array &un, Array &vn, Array &wn, Array &tn, Array &pn_dyn, Array &cn )
{
// Restore from old to new values

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				tn.x[ i ][ j ][ k ] = coeff * t.x[ i ][ j ][ k ];
				cn.x[ i ][ j ][ k ] = coeff * c.x[ i ][ j ][ k ];
				un.x[ i ][ j ][ k ] = coeff * u.x[ i ][ j ][ k ];
				vn.x[ i ][ j ][ k ] = coeff * v.x[ i ][ j ][ k ];
				wn.x[ i ][ j ][ k ] = coeff * w.x[ i ][ j ][ k ];
				pn_dyn.x[ i ][ j ][ k ] = coeff * p_dyn.x[ i ][ j ][ k ];
			}
		}
	}
}


void Restore::restoreOldNew_2D ( double coeff, Array &v, Array &w, Array &p_dyn, Array &vn, Array &wn, Array &pn_dyn )
{
// Restore of velocity components and temperature at sea surface for the next time step

			for ( int j = 0; j < jm; j++ )
			{
				for ( int k = 0; k < km; k++ )
				{
					pn_dyn.x[ im-1 ][ j ][ k ] = coeff * p_dyn.x[ im-1 ][ j ][ k ];
					vn.x[ im-1 ][ j ][ k ] = coeff * v.x[ im-1 ][ j ][ k ];
					wn.x[ im-1 ][ j ][ k ] = coeff * w.x[ im-1 ][ j ][ k ];
				}
			}
}

