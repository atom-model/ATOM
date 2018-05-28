/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to prepare the boundary and initial conditions for diverse variables
*/


#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>

#include "BC_Thermohalin.h"
#include "Array.h"

using namespace std;

BC_Thermohalin::BC_Thermohalin ( int im, int jm, int km, int i_beg, int i_max, int Ma, int Ma_max, int Ma_max_half, double dr, double g, double r_0_water, double ua, double va, double wa, double ta, double ca, double pa, double u_0, double p_0, double t_0, double c_0, double cp_w, double L_hyd, double t_average, double t_cretaceous_max, double t_equator, double t_pole, const string &input_path )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> i_max = i_max;
	this -> i_beg = i_beg;
	this -> Ma = Ma;
	this -> Ma_max = Ma_max;
	this -> Ma_max_half = Ma_max_half;
	this -> dr = dr;
	this -> g = g;
	this -> r_0_water = r_0_water;
	this -> ua = ua;
	this -> va = va;
	this -> wa = wa;
	this -> ta = ta;
	this -> ca = ca;
	this -> pa = pa;
	this -> u_0 = u_0;
	this -> p_0 = p_0;
	this -> t_0 = t_0;
	this -> c_0 = c_0;
	this -> cp_w = cp_w;
	this -> L_hyd = L_hyd;
	this -> t_average = t_average;
	this -> t_cretaceous_max = t_cretaceous_max;
	this -> t_equator = t_equator;
	this -> t_pole = t_pole;
	this -> input_path = input_path;


	i_middle = i_beg + 36;							// 0 + 36 = 36 for total depth 100m ( i_beg= 0 ), asymmetric with depth for stepsizes of 25m
	i_u_0 = i_beg + 20;								// 0 + 20 = 20 for total depth 500m ( i_beg= 0 ), asymmetric with depth for stepsizes of 25m

	j_half = ( jm - 1 ) / 2;
	j_half_1 = j_half + 1;

	d_i_beg = ( double ) i_beg;
	d_i_middle = ( double ) i_middle;
	d_i_max = ( double ) ( im - 1 );


// reduction or amplification of flow velocities along coasts
// for the artificial initial and boundary conditions
//	u_max = .01 / u_0;				// numerator radial velocity u-max in m/s          u_max = 0.01 m/s
	u_max = .1 / u_0;				// numerator radial velocity u-max in m/s          u_max = 0.01 m/s


// ocean surface velocity is about 3% of the wind velocity at the surface
	water_wind = .03;

	pi180 = 180./M_PI;

// Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
// a further turning downwards until the end of the shear layer such that finally 90° of turning are reached
	Ekman_angle = 45.0 / pi180;
}


BC_Thermohalin::~BC_Thermohalin(){}


void BC_Thermohalin::IC_v_w_EkmanSpiral ( Array_1D & rad, Array_1D & the, Array &h, Array &v, Array &w )
{
// initial conditions for v and w velocity components at the sea surface
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ im - 1 ][ j ][ k ] != 1. )
			{
				v.x[ im - 1 ][ j ][ k ] = water_wind * v.x[ im - 1 ][ j ][ k ] / u_0;
				w.x[ im - 1 ][ j ][ k ] = water_wind * w.x[ im - 1 ][ j ][ k ] / u_0;
			}
		}
	}


// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° )
	for ( int j = 0; j < 31; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ im-1 ][ j ][ k ] != 1. )
			{
				vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

				if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

				alfa = asin ( fabs ( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

				if ( w.x[ im-1 ][ j ][ k ] >= 0. )
				{
					angle = + alfa - Ekman_angle;
				}
				else
				{
					angle = - alfa - Ekman_angle;
				}

				if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}
			}
		}
	}


// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° )
	for ( int j = 31; j < 60; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ im-1 ][ j ][ k ] != 1. )
			{
				vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

				if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

				alfa = asin ( fabs ( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

				if ( v.x[ im-1 ][ j ][ k ] >= 0. )
				{
					angle = + alfa - Ekman_angle;
				}
				else
				{
					angle = - alfa - Ekman_angle;
				}

				if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = - vel_magnitude * sin ( angle );
				}

				if (  v.x[ im-1 ][ j ][ k ] <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = - vel_magnitude * sin ( angle );
				}

			}
		}
	}


// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° )
	for ( int j = 60; j < 91; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ im-1 ][ j ][ k ] != 1. )
			{
				vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

				if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

				alfa = asin ( fabs ( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

				if ( w.x[ im-1 ][ j ][ k ] >= 0. )
				{
					angle = + alfa - Ekman_angle;
				}
				else
				{
					angle = - alfa - Ekman_angle;
				}

				if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}
			}
		}
	}


// south equatorial Hadley cell ( from j=90 till j=120 compares to 0° till 30° )
	for ( int j = 91; j < 122; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ im-1 ][ j ][ k ] != 1. )
			{
				vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

				if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

				beta = asin ( fabs( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

				if ( w.x[ im-1 ][ j ][ k ] >= 0. )
				{
					angle = beta - Ekman_angle;
				}
				else
				{
					angle = - beta - Ekman_angle;
				}

				if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}
			}
		}
	}


// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° )
	for ( int j = 121; j < 151; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ im-1 ][ j ][ k ] != 1. )
			{
				vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

				if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

				beta = asin ( fabs( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

				if ( v.x[ im-1 ][ j ][ k ] <= 0. )
				{
					angle = beta - Ekman_angle;
				}
				else
				{
					angle = - beta - Ekman_angle;
				}

				if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = - vel_magnitude * sin ( angle );
				}

				if (  v.x[ im-1 ][ j ][ k ] >= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = - vel_magnitude * sin ( angle );
				}
			}
		}
	}


// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° )
	for ( int j = 151; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ im-1 ][ j ][ k ] != 1. )
			{
				vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

				if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
				if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

				beta = asin ( fabs( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

				if ( w.x[ im-1 ][ j ][ k ] >= 0. )
				{
					angle = beta - Ekman_angle;
				}
				else
				{
					angle = - beta - Ekman_angle;
				}

				if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}

				if (  w.x[ im-1 ][ j ][ k ] <= 0. )
				{
					v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
					w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
				}
			}
		}
	}


// surface wind vector driving the Ekman spiral in the Ekman layer
// northern hemisphere
	double Omega = 7.29e-5;
	double f = 0.;
	double K = 10.;
	double gam_f_K = 0.;
	double gam_z = 0.;
	double exp_gam_z = 0.;
	double sin_gam_z = 0;
	double cos_gam_z = 0;
	double v_g = 0.;
	double w_g = 0.;

	for ( int j = 0; j < jm; j++ )
	{
		f = 2. * Omega * abs ( sin ( the.z[ j ] ) );
		gam_f_K = sqrt ( abs ( f ) / ( 2. * K ) );

		for ( int i = im-2; i >= i_beg; i-- )
		{
//			gam_z = M_PI * ( double ) ( i - i_beg ) / ( double ) ( im - 1 - i_beg ) * gam_f_K;
			gam_z = M_PI * ( double ) ( i - i_beg ) / ( double ) ( im - 1 - i_beg );
			exp_gam_z = exp ( - gam_z );
			sin_gam_z = sin ( gam_z );
			cos_gam_z = cos ( gam_z );

			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v_g = v.x[ i + 1 ][ j ][ k ];
					w_g = w.x[ i + 1 ][ j ][ k ];

					if ( j <= j_half )
					{
						v.x[ i ][ j ][ k ] = w_g * exp_gam_z * sin_gam_z + v_g * ( 1. - exp_gam_z * cos_gam_z );
						w.x[ i ][ j ][ k ] = w_g * ( 1. - exp_gam_z * cos_gam_z ) - v_g * exp_gam_z * sin_gam_z;
					}
					else
					{
						sin_gam_z = - sin ( gam_z );
						v.x[ i ][ j ][ k ] = w_g * exp_gam_z * sin_gam_z + v_g * ( 1. - exp_gam_z * cos_gam_z );
						w.x[ i ][ j ][ k ] = w_g * ( 1. - exp_gam_z * cos_gam_z ) - v_g * exp_gam_z * sin_gam_z;
					}
				}
			}
		}
	}
}







void BC_Thermohalin::IC_v_w_WestEastCoast ( Array &h, Array &u, Array &v, Array &w, Array &c )
{
// initial conditions for v and w velocity components at the sea surface close to east or west coasts
// reversal of v velocity component between north and south equatorial current ommitted at respectively 10°
// w component unchanged

// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: east coast
	k_grad = 5;																			// extension of velocity change

	int l_add = 0;

	double d_k = 0.;
	double d_l = 0.;

	for ( int j = 0; j < jm; j++ )													// outer loop: latitude
	{
		for ( int k = 0; k < km - k_grad; k++ )											// inner loop: longitude
		{
			if ( ( h.x[ i_max ][ j ][ k ] == 1. ) && ( h.x[ i_max ][ j ][ k + 1 ] == 0. ) )			//  then water is closest to coast
			{
				l_add = 0;
				for ( int l = k; l < k + k_grad; l++ )								// extension of change, sign change in v-velocity and distribution of u-velocity with depth
				{
					d_l = ( double ) ( k_grad + l_add );
					d_k = ( double ) ( k_grad );

					for ( int i = i_middle; i <= i_max; i++ )
					{
						d_i = ( double ) i;
						u.x[ i ][ j ][ l ] = - ( d_i - d_i_max ) / ( d_i_middle - d_i_max ) * u_max * ( d_l * d_l / ( d_k * d_k ) );
					}

					for ( int i = i_u_0; i <= i_middle; i++ )
					{
						d_i = ( double ) i;
						u.x[ i ][ j ][ l ] = - ( d_i - d_i_beg ) / ( d_i_middle - d_i_beg ) * u_max * ( d_l * d_l / ( d_k * d_k ) );
					}
					l_add--;


					for ( int i = i_u_0; i <= i_max; i++ )
					{
						for ( int l = k; l < ( k + k_grad ); l++ )		// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension
						{
							v.x[ i ][ j ][ l ] = v.x[ i ][ j ][ k + k_grad ] / ( double )( ( k + k_grad ) - k ) * ( double )( l - k ); // extension of v-velocity
							w.x[ i ][ j ][ l ] = w.x[ i ][ j ][ k + k_grad ] / ( double )( ( k + k_grad ) - k ) * ( double )( l - k ); // extension of v-velocity
						}
					}
				}
			}
		}																						// end of longitudinal loop
	}																							// end of latitudinal loop



// search for west coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: west coast
	for ( int j = 0; j < jm; j++ )													// outer loop: latitude
	{
		for ( int k = km - 2; k > k_grad; k-- )											// inner loop: longitude
		{
			if ( ( h.x[ i_max ][ j ][ k ] == 1. ) && ( h.x[ i_max ][ j ][ k - 1 ] == 0. ) )			//  then water is closest to coast
			{
				l_add = 0;

				for ( int l = k; l > ( k - k_grad ); l-- )						// backward extention of velocity change: nothing changes
				{
					d_l = ( double ) ( k_grad - l_add );
					d_k = ( double ) ( k_grad );

					for ( int i = i_middle; i <= i_max; i++ )
					{
						d_i = ( double ) i;
						u.x[ i ][ j ][ l ] = + ( d_i - d_i_max ) / ( d_i_middle - d_i_max ) * u_max * ( d_l * d_l / ( d_k * d_k ) );
					}

					for ( int i = i_u_0; i <= i_middle; i++ )
					{
						d_i = ( double ) i;
						u.x[ i ][ j ][ l ] = + ( d_i - d_i_beg ) / ( d_i_middle - d_i_beg ) * u_max * ( d_l * d_l / ( d_k * d_k ) );
					}
					l_add++;


					for ( int i = i_u_0; i <= i_max; i++ )
					{
						for ( int l = k; l > ( k - k_grad ); l-- )			// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension
						{
							v.x[ i ][ j ][ l ] = v.x[ i ][ j ][ k - k_grad ] / ( double )( ( k - k_grad ) - k ) * ( double )( l - k ); // extension of v-velocity
							w.x[ i ][ j ][ l ] = w.x[ i ][ j ][ k - k_grad ] / ( double )( ( k - k_grad ) - k ) * ( double )( l - k ); // extension of v-velocity
						}
					}
				}
			}
		}
	}

}








void BC_Thermohalin::BC_Temperature_Salinity ( Array &h, Array &t, Array &c, Array &p_dyn )
{
// initial conditions for salt content and temperature and salinity decrease below the sea surface

// boundary condition of  temperature and salinity
// parabolic distribution from pole to pole accepted
// maximum temperature and salinity of earth's surface at equator t_max = 1.137 compares to 37° C compares to 307 K
// maximum temperature and salinity of earth's surface at equator t_max = 1.099 compares to 27° C compares to 297 K
// polar temperature and salinity of earth's surface at north and south pole t_pol = 1.0 compares to -3° C compares to 270 K
// minimum temperature and salinity at the tropopause t_min = .789 compares to -60° C compares to 213 K
// minimum temperature and salinity in the deep ocean t_deep_ocean = 1.026 compares to 4° C compares to 277 K

	j_half = ( jm -1 ) / 2;
	j_max = jm - 1;

	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

// temperature-distribution by Ruddiman approximated by a parabola
	t_coeff = t_pole - t_equator;
	t_cretaceous_coeff = t_cretaceous_max / ( ( double ) Ma_max_half - ( double ) ( Ma_max_half * Ma_max_half / Ma_max ) );   // in °C
	t_cretaceous = t_cretaceous_coeff * ( double ) ( - ( Ma * Ma ) / Ma_max + Ma );   // in °C
	if ( Ma == 0 ) 	t_cretaceous = 0.;


	cout.precision ( 3 );

	time_slice_comment = "      time slice of Cretaceous-AGCM:";
	time_slice_number = " Ma = ";
	time_slice_unit = " million years";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << time_slice_comment << resetiosflags ( ios::left ) << setw ( 6 ) << fixed << setfill ( ' ' ) << time_slice_number << setw ( 3 ) << Ma << setw ( 12 ) << time_slice_unit << endl << endl;


	temperature_comment = "      temperature increase at cretaceous times: ";
	temperature_gain = " t increase";
	temperature_modern = "      mean temperature at modern times: ";
	temperature_cretaceous = "      mean temperature at cretaceous times: ";
	temperature_average = " t modern";
	temperature_average_cret = " t cretaceous";
	temperature_unit =  "°C ";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << temperature_comment << resetiosflags ( ios::left ) << setw ( 12 ) << temperature_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << t_cretaceous << setw ( 5 ) << temperature_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_modern << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average << setw ( 5 ) << temperature_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_cretaceous << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average_cret  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average + t_cretaceous << setw ( 5 ) << temperature_unit << endl;


	c_average = ( t_average + 346. ) / 10.;									// in psu, relation taken from "Ocean Circulation, The Open University"
	c_cretaceous = ( t_average + t_cretaceous + 346. ) / 10.;		// in psu
	c_cretaceous = c_cretaceous - c_average;
	if ( Ma == 0 ) 	c_cretaceous = 0.;

	salinity_comment = "      salinity increase at cretaceous times: ";
	salinity_gain = " c increase";
	salinity_modern = "      mean salinity at modern times: ";
	salinity_cretaceous = "      mean salinity at cretaceous times: ";
	salinity_average = " c modern";
	salinity_average_cret = " c cretaceous";
	salinity_unit =  "psu ";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << salinity_comment << resetiosflags ( ios::left ) << setw ( 12 ) << salinity_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << c_cretaceous << setw ( 5 ) << salinity_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << salinity_modern << resetiosflags ( ios::left ) << setw ( 13 ) << salinity_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << c_average << setw ( 5 ) << salinity_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << salinity_cretaceous << resetiosflags ( ios::left ) << setw ( 13 ) << salinity_average_cret  << " = "  << setw ( 7 )  << setfill ( ' ' ) << c_average + c_cretaceous << setw ( 5 ) << salinity_unit << endl;
	cout << endl;


	t_cretaceous = ( t_cretaceous + t_average + t_0 ) / t_0 - ( ( t_average + t_0 ) / t_0 );    // non-dimensional

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( h.x[ im-1 ][ j ][ k ] == 0. )
			{
				d_j = ( double ) j;
				if ( Ma == 0 )			t.x[ im-1 ][ j ][ k ] = t.x[ im-1 ][ j ][ k ] + t_cretaceous; // paleo surface temperature
				else						t.x[ im-1 ][ j ][ k ] = t_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + t_pole + t_cretaceous; // paleo surface temperature

				t_Celsius = t.x[ im-1 ][ j ][ k ] * t_0 - t_0;

				if ( t_Celsius < 4. )
				{
					t.x[ im -1 ][ j ][ k ] = ta + t_cretaceous;
					t_Celsius = ( ta + t_cretaceous ) * t_0 - t_0;										// water temperature not below 4°C
				}

				if ( Ma == 0 )			c.x[ im-1 ][ j ][ k ] = c.x[ im-1 ][ j ][ k ] + c_cretaceous / c_0; // paleo surface temperature
				else 						c.x[ im-1 ][ j ][ k ] = c.x[ im-1 ][ j ][ k ];	// non-dimensional
			}
			else
			{
				t.x[ im-1 ][ j ][ k ] = ta + t_cretaceous;
				c.x[ im-1 ][ j ][ k ] = ca + c_cretaceous / c_0;
			}
		}
	}



	double tm_tbeg = 0.;
	double cm_cbeg = 0.;

// distribution of t and c with increasing depth till i_beg valid for all time slices including the modern world
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			tm_tbeg = ( t.x[ im-1 ][ j ][ k ] - ta ) / ( double ) ( i_max * i_max - i_beg * i_beg );
			cm_cbeg = ( c.x[ im-1 ][ j ][ k ] - ca ) / ( double ) ( i_max * i_max - i_beg * i_beg );
			for ( int i = i_beg; i < im-1; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
//					t.x[ i ][ j ][ k ] = ( t.x[ im-1 ][ j ][ k ] - ta ) * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg ); // linear approach
					t.x[ i ][ j ][ k ] = ta + tm_tbeg * ( double ) ( i * i - i_beg * i_beg );				// parabolic approach

					t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;
					c.x[ i ][ j ][ k ] = ca + cm_cbeg * ( double ) ( i * i - i_beg * i_beg );				// parabolic approach

					if ( t_Celsius < 4. )
					{
						t.x[ i ][ j ][ k ] = ta;
						c.x[ i ][ j ][ k ] = ca;
					}

				}
				else
				{
					t.x[ i ][ j ][ k ] = ta;
					c.x[ i ][ j ][ k ] = ca;
				}
			}
		}
	}


	for ( int i = 0; i < i_beg; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					t.x[ i ][ j ][ k ] = ta;
					c.x[ i ][ j ][ k ] = ca;
				}
			}
		}
	}

}






void BC_Thermohalin::BC_Pressure_Density ( Array &p_stat, Array &r_water, Array &r_salt_water, Array &t, Array &c, Array &h )
{
// hydrostatic pressure, equations of state for water and salt water density
// as functions of salinity, temperature and hydrostatic pressure

	double t_Celsius_0 = 0.;
	double t_Celsius_1 = 0.;
	double p_km = 0.;
	double C_p = 0.;
	double beta_p =  0.;
	double alfa_t_p =  0.;
	double gamma_t_p =  0.;

	double E_water = 2.15e9; 								// given in N/m²
	double beta_water = 8.8e-5;							// given in m³/( m³ * °C)
	double r_air = 1.2041;										// given in kg/m³
	double R_Air = 287.1;										// given in J/( kg*K )

// hydrostatic pressure, water and salt water density at the surface
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			p_stat.x[ im-1 ][ j ][ k ] =  .01 * ( r_air * R_Air * t.x[ im-1 ][ j ][ k ] * t_0 ) / 1000.;		// given in bar, isochoric approach, constant air density at the surface
			r_water.x[ im-1 ][ j ][ k ] = r_0_water;				// given in kg/m³

			t_Celsius_1 = t.x[ im-1 ][ j ][ k ] * t_0 - t_0;
			p_km = 0.;
			C_p = 999.83;
			beta_p = .808;
			alfa_t_p = .0708 * ( 1. + .068 * t_Celsius_1 );
			gamma_t_p = .003 * ( 1. - .012 * t_Celsius_1 );

			r_salt_water.x[ im-1 ][ j ][ k ] = C_p + beta_p * c.x[ im-1 ][ j ][ k ] * c_0 - alfa_t_p * t_Celsius_1 - gamma_t_p * ( c_0 - c.x[ im-1 ][ j ][ k ] * c_0 ) * t_Celsius_1;
		}
	}


// hydrostatic pressure, water and salt water density in the field
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = im-2; i >= 0; i-- )
			{
				d_i = ( double ) ( im - 1 - i );

				t_Celsius_1 = t.x[ i ][ j ][ k ] * t_0 - t_0;
				t_Celsius_0 = t.x[ i + 1 ][ j ][ k ] * t_0 - t_0;

				p_stat.x[ i ][ j ][ k ] = r_0_water * g * d_i * ( L_hyd / ( double ) ( im-1 ) ) / 100000. + p_0 / 1000.;				// hydrostatic pressure in bar
				r_water.x[ i ][ j ][ k ] = r_water.x[ i + 1 ][ j ][ k ] / ( 1. + beta_water * ( t_Celsius_1 - t_Celsius_0 ) ) / ( 1. - ( p_stat.x[ i ][ j ][ k ] - p_stat.x[ i+1 ][ j ][ k ] ) / E_water * 1e5 );				// given in kg/m³

				p_km  = ( double ) ( im - 1 - i ) * ( L_hyd / ( double ) ( im-1 ) ) / 1000.;
				C_p = 999.83 + 5.053 * p_km - .048 * p_km * p_km;
				beta_p = .808 - .0085* p_km;
				alfa_t_p = .0708 * ( 1. + .351 * p_km + .068 *( 1. - .0683 * p_km ) * t_Celsius_1 );
				gamma_t_p = .003 * ( 1. - .059 * p_km - .012 * ( 1. - .064 * p_km ) * t_Celsius_1 );

				r_salt_water.x[ i ][ j ][ k ] = C_p + beta_p * c.x[ i ][ j ][ k ] * c_0 - alfa_t_p * t_Celsius_1 - gamma_t_p * ( c_0 - c.x[ i ][ j ][ k ] * c_0 ) * t_Celsius_1;
			}
		}
	}
}






void BC_Thermohalin::BC_Surface_Temperature_NASA ( const string &Name_SurfaceTemperature_File, Array &t )
{
	// initial conditions for the temperature and salinity at the sea surface

	cout.precision ( 3 );
	cout.setf ( ios::fixed );

	ifstream Name_SurfaceTemperature_File_Read;
	Name_SurfaceTemperature_File_Read.open(Name_SurfaceTemperature_File);

	if (!Name_SurfaceTemperature_File_Read.is_open()) {
		cerr << "ERROR: could not open SurfaceTemperature_File file at " << Name_SurfaceTemperature_File << "\n";
		abort();
	}


	j = 0;
	k = 0;

	while ( ( k < km ) && ( !Name_SurfaceTemperature_File_Read.eof() ) ) {
		while ( j < jm ) {
			Name_SurfaceTemperature_File_Read >> dummy_1;
			Name_SurfaceTemperature_File_Read >> dummy_2;
			Name_SurfaceTemperature_File_Read >> dummy_3;

			t.x[ im-1 ][ j ][ k ] = ( dummy_3 + 273.15 ) / 273.15;

			j++;
		}
		j = 0;
		k++;
	}

	Name_SurfaceTemperature_File_Read.close();


	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			if ( k == 180 )	t.x[ im-1 ][ j ][ k ] = ( t.x[ im-1 ][ j ][ k + 1 ] + t.x[ im-1 ][ j ][ k - 1 ] ) * .5;
		}
	}

}






void BC_Thermohalin::BC_Surface_Salinity_NASA ( const string &Name_SurfaceSalinity_File, Array &c )
{
	// initial conditions for the salinity at the sea surface
	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 3 );
	cout.setf ( ios::fixed );

	ifstream Name_SurfaceSalinity_File_Read;
	Name_SurfaceSalinity_File_Read.open(Name_SurfaceSalinity_File);

	if (!Name_SurfaceSalinity_File_Read.is_open()) {
		cerr << "ERROR: could not open SurfaceSalinity_File file at " << Name_SurfaceSalinity_File << "\n";
		abort();
	}

	j = 0;
	k = 0;

	while ( ( k < km ) && ( !Name_SurfaceSalinity_File_Read.eof() ) ) {
		while ( j < jm ) {
			Name_SurfaceSalinity_File_Read >> dummy_1;
			Name_SurfaceSalinity_File_Read >> dummy_2;
			Name_SurfaceSalinity_File_Read >> dummy_3;

			if ( dummy_3 < 0. ) dummy_3 = ca;

			else 		c.x[ im-1 ][ j ][ k ] = dummy_3 / c_0;
			j++;
		}
		j = 0;
		k++;
	}

	Name_SurfaceSalinity_File_Read.close();
}





void BC_Thermohalin::IC_CircumPolar_Current ( Array &h, Array &u, Array &v, Array &w, Array &c )
{
// south polar sea
// antarctic circumpolar current ( -5000m deep ) ( from j=147 until j=152 compares to 57°S until 62°S,
//                                                                            from k=0 until k=km compares to 0° until 360° )


	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 147; j < 153; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
//					c.x[ i ][ j ][ k ] = ca;
					w.x[ i ][ j ][ k ] = .1 / u_0;					// 10 cm/s
				}
			}
		}
	}

}
