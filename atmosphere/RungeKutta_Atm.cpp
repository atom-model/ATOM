/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * class to produce results by the Runge-Kutta solution scheme
*/

#include <iostream>
#include <cmath>
#include "RungeKutta_Atm.h"
#include "RHS_Atm.h"

using namespace std;


RungeKutta_Atmosphere::RungeKutta_Atmosphere ( int im, int jm, int km, double dt, double dr, double dthe, double dphi )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> dt = dt;
	this -> dr = dr;
	this -> dphi = dphi;
	this -> dthe = dthe;
}

RungeKutta_Atmosphere::~RungeKutta_Atmosphere () {}


void RungeKutta_Atmosphere::solveRungeKutta_3D_Atmosphere ( RHS_Atmosphere &prepare, int &n, double lv, double ls, double ep, double hp, double u_0, double t_0, double c_0, double co2_0, double p_0, double r_air, double r_water, double r_water_vapour, double r_co2, double L_atm, double cp_l, double R_Air, double R_WaterVapour, double R_co2, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c, Array &rhs_cloud, Array &rhs_ice, Array &rhs_co2, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &p_stat, Array &c, Array &cloud, Array &ice, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &p_dynn, Array &cn, Array &cloudn, Array &icen, Array &co2n, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &BuoyancyForce, Array &Q_Sensible, Array &P_rain, Array &P_snow, Array &S_v, Array &S_c, Array &S_i, Array &S_r, Array &S_s, Array &S_c_c, Array_2D &Topography, Array_2D &Evaporation_Dalton, Array_2D &Precipitation )
{
// Runge-Kutta 4. order for u, v and w component, temperature, water vapour and co2 content

	for ( int i = 1; i < im-1; i++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
// Runge-Kutta 4. order for k1 step ( dt )

				tn.x[ i ][ j ][ k ] = t.x[ i ][ j ][ k ];
				un.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ];
				vn.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ];
				wn.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ];
				cn.x[ i ][ j ][ k ] =  c.x[ i ][ j ][ k ];
				cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ];
				icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ];
				co2n.x[ i ][ j ][ k ] = co2.x[ i ][ j ][ k ];


//	if ( ( i == 11 ) && ( j == 70 ) && ( k == 160 ) ) 				cout << "   vor prepare      i = " << i << "   j = " << j << "   k = " << k << "   u = " << u.x[ i ][ j ][ k ] << "   un = " << un.x[ i ][ j ][ k ] << "   v = " << v.x[ i ][ j ][ k ] << "   vn = " << vn.x[ i ][ j ][ k ] << "   w = " << w.x[ i ][ j ][ k ] << "   wn = " << wn.x[ i ][ j ][ k ] << "   p_dyn = " << p_dyn.x[ i ][ j ][ k ] << "   p_dynn = " << p_dynn.x[ i ][ j ][ k ] << endl;


				prepare.RK_RHS_3D_Atmosphere ( n, i, j, k, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, p_dynn, cn, cloudn, icen, co2n, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_cloud, rhs_ice, rhs_co2, aux_u, aux_v, aux_w, Latency, BuoyancyForce, Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c, Topography, Evaporation_Dalton, Precipitation );

//	if ( ( i == 11 ) && ( j == 70 ) && ( k == 160 ) ) 				cout << "   hinter prepare   i = " << i << "   j = " << j << "   k = " << k << "   u = " << u.x[ i ][ j ][ k ] << "   un = " << un.x[ i ][ j ][ k ] << "   v = " << v.x[ i ][ j ][ k ] << "   vn = " << vn.x[ i ][ j ][ k ] << "   w = " << w.x[ i ][ j ][ k ] << "   wn = " << wn.x[ i ][ j ][ k ] << "   p_dyn = " << p_dyn.x[ i ][ j ][ k ] << "   p_dynn = " << p_dynn.x[ i ][ j ][ k ] << endl;

				kt1 = rhs_t.x[ i ][ j ][ k ];
				ku1 = rhs_u.x[ i ][ j ][ k ];
				kv1 = rhs_v.x[ i ][ j ][ k ];
				kw1 = rhs_w.x[ i ][ j ][ k ];
				kc1 = rhs_c.x[ i ][ j ][ k ];
				kcloud1 = rhs_cloud.x[ i ][ j ][ k ];
				kice1 = rhs_ice.x[ i ][ j ][ k ];
				kco1 = rhs_co2.x[ i ][ j ][ k ];

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt1 * .5 * dt;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku1 * .5 * dt;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv1 * .5 * dt;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw1 * .5 * dt;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc1 * .5 * dt;
				cloud.x[ i ][ j ][ k ] = cloudn.x[ i ][ j ][ k ] + kcloud1 * .5 * dt;
				ice.x[ i ][ j ][ k ] = icen.x[ i ][ j ][ k ] + kice1 * .5 * dt;
				co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco1 * .5 * dt;

//	if ( ( i == 11 ) && ( j == 70 ) && ( k == 160 ) ) 				cout << "   i = " << i << "   j = " << j << "   k = " << k << "   ku1 = " << ku1 << "   kv1 = " << kv1 << "   kw1 = " << kw1 << "   kp1 = " << kp1 << "   u = " << u.x[ i ][ j ][ k ] << "   un = " << un.x[ i ][ j ][ k ] << "   v = " << v.x[ i ][ j ][ k ] << "   vn = " << vn.x[ i ][ j ][ k ] << "   w = " << w.x[ i ][ j ][ k ] << "   wn = " << wn.x[ i ][ j ][ k ] << "   p_dyn = " << p_dyn.x[ i ][ j ][ k ] << "   p_dynn = " << p_dynn.x[ i ][ j ][ k ] << endl;



// Runge-Kutta 4. order for k2 step ( dt )
				prepare.RK_RHS_3D_Atmosphere ( n, i, j, k, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, p_dynn, cn, cloudn, icen, co2n, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_cloud, rhs_ice, rhs_co2, aux_u, aux_v, aux_w, Latency, BuoyancyForce, Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c, Topography, Evaporation_Dalton, Precipitation );

				kt2 = rhs_t.x[ i ][ j ][ k ];
				ku2 = rhs_u.x[ i ][ j ][ k ];
				kv2 = rhs_v.x[ i ][ j ][ k ];
				kw2 = rhs_w.x[ i ][ j ][ k ];
				kc2 = rhs_c.x[ i ][ j ][ k ];
				kcloud2 = rhs_cloud.x[ i ][ j ][ k ];
				kice2 = rhs_ice.x[ i ][ j ][ k ];
				kco2 = rhs_co2.x[ i ][ j ][ k ];

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt2 * .5 * dt;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku2 * .5 * dt;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv2 * .5 * dt;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw2 * .5 * dt;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc2 * .5 * dt;
				cloud.x[ i ][ j ][ k ] = cloudn.x[ i ][ j ][ k ] + kcloud2 * .5 * dt;
				ice.x[ i ][ j ][ k ] = icen.x[ i ][ j ][ k ] + kice2 * .5 * dt;
				co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco2 * .5 * dt;

//	if ( ( i == 11 ) && ( j == 70 ) && ( k == 160 ) ) 				cout << "   i = " << i << "   j = " << j << "   k = " << k << "   ku2 = " << ku2 << "   kv2 = " << kv2 << "   kw2 = " << kw2 << "   kp2 = " << kp2 << "   u = " << u.x[ i ][ j ][ k ] << "   un = " << un.x[ i ][ j ][ k ] << "   v = " << v.x[ i ][ j ][ k ] << "   vn = " << vn.x[ i ][ j ][ k ] << "   w = " << w.x[ i ][ j ][ k ] << "   wn = " << wn.x[ i ][ j ][ k ] << "   p_dyn = " << p_dyn.x[ i ][ j ][ k ] << "   p_dynn = " << p_dynn.x[ i ][ j ][ k ] << endl;

// Runge-Kutta 4. order for k3 step ( dt )
				prepare.RK_RHS_3D_Atmosphere ( n, i, j, k, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, p_dynn, cn, cloudn, icen, co2n, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_cloud, rhs_ice, rhs_co2, aux_u, aux_v, aux_w, Latency, BuoyancyForce, Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c, Topography, Evaporation_Dalton, Precipitation );

				kt3 = rhs_t.x[ i ][ j ][ k ];
				ku3 = rhs_u.x[ i ][ j ][ k ];
				kv3 = rhs_v.x[ i ][ j ][ k ];
				kw3 = rhs_w.x[ i ][ j ][ k ];
				kc3 = rhs_c.x[ i ][ j ][ k ];
				kcloud3 = rhs_cloud.x[ i ][ j ][ k ];
				kice3 = rhs_ice.x[ i ][ j ][ k ];
				kco3 = rhs_co2.x[ i ][ j ][ k ];

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + kt3 * dt;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + ku3 * dt;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + kv3 * dt;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + kw3 * dt;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + kc3 * dt;
				cloud.x[ i ][ j ][ k ] = cloudn.x[ i ][ j ][ k ] + kcloud3 * dt;
				ice.x[ i ][ j ][ k ] = icen.x[ i ][ j ][ k ] + kice3 * dt;
				co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + kco3 * dt;

// Runge-Kutta 4. order for k4 step ( dt )
				prepare.RK_RHS_3D_Atmosphere ( n, i, j, k, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, p_dynn, cn, cloudn, icen, co2n, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_cloud, rhs_ice, rhs_co2, aux_u, aux_v, aux_w, Latency, BuoyancyForce, Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s, S_c_c, Topography, Evaporation_Dalton, Precipitation );

				kt4 = rhs_t.x[ i ][ j ][ k ];
				ku4 = rhs_u.x[ i ][ j ][ k ];
				kv4 = rhs_v.x[ i ][ j ][ k ];
				kw4 = rhs_w.x[ i ][ j ][ k ];
				kc4 = rhs_c.x[ i ][ j ][ k ];
				kcloud4 = rhs_cloud.x[ i ][ j ][ k ];
				kice4 = rhs_ice.x[ i ][ j ][ k ];
				kco4 = rhs_co2.x[ i ][ j ][ k ];

				t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] + dt * ( kt1 + 2. * kt2 + 2. * kt3 + kt4 ) / 6.;
				u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] + dt * ( ku1 + 2. * ku2 + 2. * ku3 + ku4 ) / 6.;
				v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] + dt * ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
				w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] + dt * ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;
				c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] + dt * ( kc1 + 2. * kc2 + 2. * kc3 + kc4 ) / 6.;
				cloud.x[ i ][ j ][ k ] = cloudn.x[ i ][ j ][ k ] + dt * ( kcloud1 + 2. * kcloud2 + 2. * kcloud3 + kcloud4 ) / 6.;
				ice.x[ i ][ j ][ k ] = icen.x[ i ][ j ][ k ] + dt * ( kice1 + 2. * kice2 + 2. * kice3 + kice4 ) / 6.;
				co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] + dt *( kco1 + 2. * kco2 + 2. * kco3 + kco4 ) / 6.;

				tn.x[ i ][ j ][ k ] = t.x[ i ][ j ][ k ];
				un.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ];
				vn.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ];
				wn.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ];
				cn.x[ i ][ j ][ k ] =  c.x[ i ][ j ][ k ];
				cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ];
				icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ];
				co2n.x[ i ][ j ][ k ] = co2.x[ i ][ j ][ k ];
			}
		}
	}
}






void RungeKutta_Atmosphere::solveRungeKutta_2D_Atmosphere ( RHS_Atmosphere &prepare_2D, int &n, double r_air, double u_0, double p_0, double L_atm, Array_1D &rad, Array_1D &the, Array &rhs_v, Array &rhs_w, Array &h, Array &v, Array &w, Array &p_dyn, Array &vn, Array &wn, Array &p_dynn, Array &aux_v, Array &aux_w )
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

			vn.x[ 0 ][ j ][ k ] = v.x[ 0 ][ j ][ k ];
			wn.x[ 0 ][ j ][ k ] = w.x[ 0 ][ j ][ k ];


			prepare_2D.RK_RHS_2D_Atmosphere ( j, k, r_air, u_0, p_0, L_atm, rad, the, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, aux_v, aux_w );

			kv1 = rhs_v.x[ 0 ][ j ][ k ];
			kw1 =rhs_w.x[ 0 ][ j ][ k ];

			v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv1 * .5 * dt;
			w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw1 * .5 * dt;

	// Runge-Kutta 4. order for k2 step ( dt )
			prepare_2D.RK_RHS_2D_Atmosphere ( j, k, r_air, u_0, p_0, L_atm, rad, the, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, aux_v, aux_w );

			kv2 = rhs_v.x[ 0 ][ j ][ k ];
			kw2 = rhs_w.x[ 0 ][ j ][ k ];

			v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv2 * .5 * dt;
			w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw2 * .5 * dt;

		// Runge-Kutta 4. order for k3 step ( dt )
			prepare_2D.RK_RHS_2D_Atmosphere ( j, k, r_air, u_0, p_0, L_atm, rad, the, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, aux_v, aux_w );

			kv3 = rhs_v.x[ 0 ][ j ][ k ];
			kw3 = rhs_w.x[ 0 ][ j ][ k ];

			v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + kv3 * dt;
			w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + kw3 * dt;

		// Runge-Kutta 4. order for k4 step ( dt )
			prepare_2D.RK_RHS_2D_Atmosphere ( j, k, r_air, u_0, p_0, L_atm, rad, the, h, v, w, p_dyn, vn, wn, p_dynn, rhs_v, rhs_w, aux_v, aux_w );

			kv4 = rhs_v.x[ 0 ][ j ][ k ];
			kw4 =rhs_w.x[ 0 ][ j ][ k ];

			v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] + dt * ( kv1 + 2. * kv2 + 2. * kv3 + kv4 ) / 6.;
			w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] + dt * ( kw1 + 2. * kw2 + 2. * kw3 + kw4 ) / 6.;

			vn.x[ 0 ][ j ][ k ] = v.x[ 0 ][ j ][ k ];
			wn.x[ 0 ][ j ][ k ] = w.x[ 0 ][ j ][ k ];


		}
	}
}
