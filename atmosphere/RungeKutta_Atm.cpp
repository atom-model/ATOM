/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include <iostream>
#include "RungeKutta_Atm.h"

using namespace std;


RungeKutta_Atmosphere::RungeKutta_Atmosphere ( int im, int jm, int km, double dt )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> dt = dt;
}

RungeKutta_Atmosphere::~RungeKutta_Atmosphere () {}


void RungeKutta_Atmosphere::solveRungeKutta_3D_Atmosphere ( RHS_Atmosphere &prepare, int *im_tropopause, double lv, double ls, double ep, double hp, double u_0, double t_0, double c_0, double co2_0, double p_0, double r_0_air, double r_0_water_vapour, double r_0_co2, double L_atm, double cp_l, double R_Air, double R_WaterVapour, double R_co2, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c, Array &rhs_co2, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &p_stat, Array &c, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &co2n, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Condensation_3D, Array &Evaporation_3D, Array &Rain, Array &Ice, Array &Rain_super, Array &IceLayer, Array &cloud, Array &BuoyancyForce, Array &Q_Sensible )
{
// Runge-Kutta 4. order for u, v and w component, temperature, water vapour and co2 content

	for ( int i = 1; i < im-1; i++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
// Runge-Kutta 4. order for k1 step ( dt )

				prepare.RK_RHS_3D_Atmosphere ( i, j, k, im_tropopause, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_0_air, r_0_water_vapour, r_0_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, h, t, u, v, w, p_dyn, p_stat, c, co2, tn, un, vn, wn, cn, co2n, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Condensation_3D, Evaporation_3D, Rain, Ice, Rain_super, IceLayer, cloud, BuoyancyForce, Q_Sensible );

				kt1 = dt * rhs_t.x[ i ][ j ][ k ];
				ku1 = dt * rhs_u.x[ i ][ j ][ k ];
				kv1 = dt * rhs_v.x[ i ][ j ][ k ];
				kw1 = dt * rhs_w.x[ i ][ j ][ k ];
				kc1 = dt * rhs_c.x[ i ][ j ][ k ];
				kco21 = dt * rhs_co2.x[ i ][ j ][ k ];

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt1 * .5;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku1 * .5;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv1 * .5;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw1 * .5;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc1 * .5;
				co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco21 * .5;

// Runge-Kutta 4. order for k2 step ( dt )

				prepare.RK_RHS_3D_Atmosphere ( i, j, k, im_tropopause, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_0_air, r_0_water_vapour, r_0_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, h, t, u, v, w, p_dyn, p_stat, c, co2, tn, un, vn, wn, cn, co2n, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Condensation_3D, Evaporation_3D, Rain, Ice, Rain_super, IceLayer, cloud, BuoyancyForce, Q_Sensible );

				kt2 = dt * rhs_t.x[ i ][ j ][ k ] * .5;
				ku2 = dt * rhs_u.x[ i ][ j ][ k ] * .5;
				kv2 = dt * rhs_v.x[ i ][ j ][ k ] * .5;
				kw2 = dt * rhs_w.x[ i ][ j ][ k ] * .5;
				kc2 = dt * rhs_c.x[ i ][ j ][ k ] * .5;
				kco22 = dt * rhs_co2.x[ i ][ j ][ k ] * .5;

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt2 * .5;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku2 * .5;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv2 * .5;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw2 * .5;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc2 * .5;
				co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco22 * .5;

// Runge-Kutta 4. order for k3 step ( dt )

				prepare.RK_RHS_3D_Atmosphere ( i, j, k, im_tropopause, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_0_air, r_0_water_vapour, r_0_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, h, t, u, v, w, p_dyn, p_stat, c, co2, tn, un, vn, wn, cn, co2n, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Condensation_3D, Evaporation_3D, Rain, Ice, Rain_super, IceLayer, cloud, BuoyancyForce, Q_Sensible );

				kt3 = dt * rhs_t.x[ i ][ j ][ k ] * .5;
				ku3 = dt * rhs_u.x[ i ][ j ][ k ] * .5;
				kv3 = dt * rhs_v.x[ i ][ j ][ k ] * .5;
				kw3 = dt * rhs_w.x[ i ][ j ][ k ] * .5;
				kc3 = dt * rhs_c.x[ i ][ j ][ k ] * .5;
				kco23 = dt * rhs_co2.x[ i ][ j ][ k ] * .5;

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt3;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku3;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv3;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw3;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc3;
				co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco23;

// Runge-Kutta 4. order for k4 step ( dt )

				prepare.RK_RHS_3D_Atmosphere ( i, j, k, im_tropopause, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_0_air, r_0_water_vapour, r_0_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, h, t, u, v, w, p_dyn, p_stat, c, co2, tn, un, vn, wn, cn, co2n, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_co2, aux_u, aux_v, aux_w, Latency, Condensation_3D, Evaporation_3D, Rain, Ice, Rain_super, IceLayer, cloud, BuoyancyForce, Q_Sensible );

				kt4 = dt * rhs_t.x[ i ][ j ][ k ];
				ku4 = dt * rhs_u.x[ i ][ j ][ k ];
				kv4 = dt * rhs_v.x[ i ][ j ][ k ];
				kw4 = dt * rhs_w.x[ i ][ j ][ k ];
				kc4 = dt * rhs_c.x[ i ][ j ][ k ];
				kco24 = dt * rhs_co2.x[ i ][ j ][ k ];

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + ( kt1 + 2. * kt2 + 2. * kt3 + kt4 ) / 6.;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ( ku1 + 2. * ku2 + 2. * ku3 + ku4 ) / 6.;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + ( kc1 + 2. * kc2 + 2. * kc3 + kc4 ) / 6.;
				co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + ( kco21 + 2. * kco22 + 2. * kco23 + kco24 ) / 6.;

			}
		}
	}
}






void RungeKutta_Atmosphere::solveRungeKutta_2D_Atmosphere ( RHS_Atmosphere &prepare, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &rhs_v, Array &rhs_w, Array &h, Array &v, Array &w, Array &p_dyn, Array &vn, Array &wn, Array &aux_v, Array &aux_w )
{
// Runge-Kutta 4. order for u, v and w component, temperature, water vapour and co2 content
//  2D surface iterations
// Runge-Kutta 4. order for u, v and w component, temperature and salt concentration

	cout.precision ( 9 );
	cout.setf ( ios::fixed );

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
// Runge-Kutta 4. order for k1 step ( dt )
			prepare.RK_RHS_2D_Atmosphere ( j, k, rad, the, phi, h, v, w, p_dyn, vn, wn, rhs_v, rhs_w, aux_v, aux_w );

			kv1 = dt * rhs_v.x[ 0 ][ j ][ k ];
			kw1 = dt * rhs_w.x[ 0 ][ j ][ k ];

			v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv1 * .5;
			w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw1 * .5;

	// Runge-Kutta 4. order for k2 step ( dt )
			prepare.RK_RHS_2D_Atmosphere ( j, k, rad, the, phi, h, v, w, p_dyn, vn, wn, rhs_v, rhs_w, aux_v, aux_w );

			kv2 = dt * rhs_v.x[ 0 ][ j ][ k ] * .5;
			kw2 = dt * rhs_w.x[ 0 ][ j ][ k ] * .5;

			v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv2 * .5;
			w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw2 * .5;

		// Runge-Kutta 4. order for k3 step ( dt )
			prepare.RK_RHS_2D_Atmosphere ( j, k, rad, the, phi, h, v, w, p_dyn, vn, wn, rhs_v, rhs_w, aux_v, aux_w );

			kv3 = dt * rhs_v.x[ 0 ][ j ][ k ] * .5;
			kw3 = dt * rhs_w.x[ 0 ][ j ][ k ] * .5;

			v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv3;
			w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw3;

		// Runge-Kutta 4. order for k4 step ( dt )
			prepare.RK_RHS_2D_Atmosphere ( j, k, rad, the, phi, h, v, w, p_dyn, vn, wn, rhs_v, rhs_w, aux_v, aux_w );

			kv4 = dt * rhs_v.x[ 0 ][ j ][ k ];
			kw4 = dt * rhs_w.x[ 0 ][ j ][ k ];

			v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
			w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
			}
		}
	}
