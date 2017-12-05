/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include <iostream>
#include <cmath>
#include "RungeKutta_Hyd.h"
#include "RHS_Hyd.h"

using namespace std;


RungeKutta_Hydrosphere::RungeKutta_Hydrosphere ( int n, int im, int jm, int km, double dt )
{
	this -> n = n;
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> dt = dt;
}

RungeKutta_Hydrosphere::~RungeKutta_Hydrosphere () {}



void RungeKutta_Hydrosphere::solveRungeKutta_3D_Hydrosphere ( RHS_Hydrosphere &prepare, int &n, int &n_pres, double L_hyd, double g, double cp_w, double u_0, double t_0, double c_0, double r_0_water, double ta, double pa, double ca, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_p, Array &rhs_c, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &p_dynn, Array &cn, Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion, Array &Buoyancy_Force, Array &Salt_Balance, Array &p_stat )
{
//  3D volume iterations
// Runge-Kutta 4. order for u, v and w component, temperature and salt concentration

	for ( int i = 1; i < im-1; i++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
// Runge-Kutta 4. order for k1 step ( dt )
				prepare.RK_RHS_3D_Hydrosphere ( i, j, k, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, rhs_t, rhs_u, rhs_v, rhs_w, rhs_p, rhs_c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Buoyancy_Force, Salt_Balance, p_stat );

				kt1 = dt * rhs_t.x[ i ][ j ][ k ];
				ku1 = dt * rhs_u.x[ i ][ j ][ k ];
				kv1 = dt * rhs_v.x[ i ][ j ][ k ];
				kw1 = dt * rhs_w.x[ i ][ j ][ k ];
				kp1 = dt * rhs_p.x[ i ][ j ][ k ];
				kc1 = dt * rhs_c.x[ i ][ j ][ k ];

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt1 * .5;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku1 * .5;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv1 * .5;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw1 * .5;
				p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ] + kp1 * .5;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc1 * .5;

// Runge-Kutta 4. order for k2 step ( dt )
				prepare.RK_RHS_3D_Hydrosphere ( i, j, k, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, rhs_t, rhs_u, rhs_v, rhs_w, rhs_p, rhs_c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Buoyancy_Force, Salt_Balance, p_stat );

				kt2 = dt * rhs_t.x[ i ][ j ][ k ] * .5;
				ku2 = dt * rhs_u.x[ i ][ j ][ k ] * .5;
				kv2 = dt * rhs_v.x[ i ][ j ][ k ] * .5;
				kw2 = dt * rhs_w.x[ i ][ j ][ k ] * .5;
				kp2 = dt * rhs_p.x[ i ][ j ][ k ] * .5;
				kc2 = dt * rhs_c.x[ i ][ j ][ k ] * .5;

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt2 * .5;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku2 * .5;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv2 * .5;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw2 * .5;
				p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ] + kp2 * .5;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc2 * .5;

	// Runge-Kutta 4. order for k3 step ( dt )
				prepare.RK_RHS_3D_Hydrosphere ( i, j, k, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, rhs_t, rhs_u, rhs_v, rhs_w, rhs_p, rhs_c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Buoyancy_Force, Salt_Balance, p_stat );

				kt3 = dt * rhs_t.x[ i ][ j ][ k ] * .5;
				ku3 = dt * rhs_u.x[ i ][ j ][ k ] * .5;
				kv3 = dt * rhs_v.x[ i ][ j ][ k ] * .5;
				kw3 = dt * rhs_w.x[ i ][ j ][ k ] * .5;
				kp3 = dt * rhs_p.x[ i ][ j ][ k ] * .5;
				kc3 = dt * rhs_c.x[ i ][ j ][ k ] * .5;

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt3;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku3;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv3;
				w.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] + kw3;
				p_dyn.x[ i ][ j ][ k ] = p_dyn.x[ i ][ j ][ k ] + kp3;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc3;

	// Runge-Kutta 4. order for k4 step ( dt )
				prepare.RK_RHS_3D_Hydrosphere ( i, j, k, L_hyd, g, cp_w, u_0, t_0, c_0, r_0_water, ta, pa, ca, rad, the, phi, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn, rhs_t, rhs_u, rhs_v, rhs_w, rhs_p, rhs_c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, Buoyancy_Force, Salt_Balance, p_stat );

				kt4 = dt * rhs_t.x[ i ][ j ][ k ];
				ku4 = dt * rhs_u.x[ i ][ j ][ k ];
				kv4 = dt * rhs_v.x[ i ][ j ][ k ];
				kw4 = dt * rhs_w.x[ i ][ j ][ k ];
				kp4 = dt * rhs_p.x[ i ][ j ][ k ];
				kc4 = dt * rhs_c.x[ i ][ j ][ k ];

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + ( kt1 + 2. * kt2 + 2. * kt3 + kt4 ) / 6.;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ( ku1 + 2. * ku2 + 2. * ku3 + ku4 ) / 6.;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
				p_dyn.x[ i ][ j ][ k ] = p_dynn.x[ i ][ j ][ k ] + ( kp1 + 2. * kp2 + 2. * kp3 + kp4 ) / 6.;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + ( kc1 + 2. * kc2 + 2. * kc3 + kc4 ) / 6.;
			}
		}
	}
}





void RungeKutta_Hydrosphere::solveRungeKutta_2D_Hydrosphere ( RHS_Hydrosphere &prepare_2D, int &n, int &n_pres, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &rhs_v, Array &rhs_w, Array &rhs_p, Array &h, Array &v, Array &w, Array &p_dyn, Array &vn, Array &wn, Array &p_dynn, Array &aux_v, Array &aux_w )
{
//  2D surface iterations
// Runge-Kutta 4. order for u, v and w component, temperature and salt concentration

	cout.precision ( 9 );
	cout.setf ( ios::fixed );

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
// Runge-Kutta 4. order for k1 step ( dt )
			prepare_2D.RK_RHS_2D_Hydrosphere ( j, k, rad, the, phi, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, rhs_p, aux_v, aux_w );

			kv1 = dt * rhs_v.x[ im-1 ][ j ][ k ];
			kw1 = dt * rhs_w.x[ im-1 ][ j ][ k ];
			kp1 = dt * rhs_p.x[ im-1 ][ j ][ k ];

			v.x[ im-1 ][ j ][ k ] = vn.x[ im-1 ][ j ][ k ] + kv1 * .5;
			w.x[ im-1 ][ j ][ k ] = wn.x[ im-1 ][ j ][ k ] + kw1 * .5;
			p_dyn.x[ im-1 ][ j ][ k ] = p_dynn.x[ im-1 ][ j ][ k ] + kp1 * .5;

	// Runge-Kutta 4. order for k2 step ( dt )
			prepare_2D.RK_RHS_2D_Hydrosphere ( j, k, rad, the, phi, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, rhs_p, aux_v, aux_w );

			kv2 = dt * rhs_v.x[ im-1 ][ j ][ k ] * .5;
			kw2 = dt * rhs_w.x[ im-1 ][ j ][ k ] * .5;
			kp2 = dt * rhs_p.x[ im-1 ][ j ][ k ] * .5;

			v.x[ im-1 ][ j ][ k ] = vn.x[ im-1 ][ j ][ k ] + kv2 * .5;
			w.x[ im-1 ][ j ][ k ] = wn.x[ im-1 ][ j ][ k ] + kw2 * .5;
			p_dyn.x[ im-1 ][ j ][ k ] = p_dynn.x[ im-1 ][ j ][ k ] + kp2 * .5;

		// Runge-Kutta 4. order for k3 step ( dt )
			prepare_2D.RK_RHS_2D_Hydrosphere ( j, k, rad, the, phi, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, rhs_p, aux_v, aux_w );

			kv3 = dt * rhs_v.x[ im-1 ][ j ][ k ] * .5;
			kw3 = dt * rhs_w.x[ im-1 ][ j ][ k ] * .5;
			kp3 = dt * rhs_p.x[ im-1 ][ j ][ k ] * .5;

			v.x[ im-1 ][ j ][ k ] = vn.x[ im-1 ][ j ][ k ] + kv3;
			w.x[ im-1 ][ j ][ k ] = wn.x[ im-1 ][ j ][ k ] + kw3;
			p_dyn.x[ im-1 ][ j ][ k ] = p_dynn.x[ im-1 ][ j ][ k ] + kp3;

		// Runge-Kutta 4. order for k4 step ( dt )
			prepare_2D.RK_RHS_2D_Hydrosphere ( j, k, rad, the, phi, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, rhs_p, aux_v, aux_w );

			kv4 = dt * rhs_v.x[ im-1 ][ j ][ k ];
			kw4 = dt * rhs_w.x[ im-1 ][ j ][ k ];
			kp4 = dt * rhs_p.x[ im-1 ][ j ][ k ];

			v.x[ im-1 ][ j ][ k ] = vn.x[ im-1 ][ j ][ k ] + ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
			w.x[ im-1 ][ j ][ k ] = wn.x[ im-1 ][ j ][ k ] + ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
			p_dyn.x[ im-1 ][ j ][ k ] = p_dynn.x[ im-1 ][ j ][ k ] + ( kp1 + 2. * kp2 + 2. * kp3 + kp4 ) / 6.;
		}
	}


	if ( n == n_pres )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				p_dyn.x[ im-1 ][ j ][ k ] = p_dynn.x[ im-1 ][ j ][ k ];
			}
		}
		n_pres = n_pres + 2;
	}
}
