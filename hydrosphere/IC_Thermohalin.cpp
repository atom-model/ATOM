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

#include "IC_Thermohalin.h"

using namespace std;




IC_Thermohalin::IC_Thermohalin ( int im, int jm, int km )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;

// assumption of maximum depth of sea 6000 m compares to 40 steps times 150 m

	i_beg = 27;				// compares to an ocean depth of 1950 m
//	i_beg = 25;				// compares to an ocean depth of 2250 m

	i_max = im - 1;
	i_bottom = 8;
	i_deep = i_bottom + 17;
	i_half = i_beg;
	i_middle = i_beg - 5;
	j_half = ( jm - 1 ) / 2;

	d_i_half = ( double ) i_half;
	d_i_max = ( double ) ( im - 1 );


// reduction or amplification of flow velocities along coasts 
// for the artificial initial and boundary conditions
// u_0 = .45 m/s

//	IC_water = 1.;
//	IC_water = .1;
	IC_water = .5;					// no dimension, ( average velocity compares to          u_0 * IC_water = 0,225 m/s )

// ocean surface velocity is about 3% of the wind velocity at the surface

	water_wind = .03;

	pi180 = 180./M_PI;

// Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
// a further turning downwards until the end of the shear layer such that finally 90° of turning are reached

	Ekman_angle = 45.0 / pi180;
	Ekman_angle_add = 4.5 / pi180;

}


IC_Thermohalin::~IC_Thermohalin() {}





void IC_Thermohalin::IC_v_w_Atmosphere ( Array &h, Array &u, Array &v, Array &w )
{
// initial conditions for v and w velocity components at the sea surface

// surface velocity reduced to 3% of the wind velocity

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = water_wind * v.x[ i_max ][ j ][ k ] * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = water_wind * w.x[ i_max ][ j ][ k ] * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
				else
				{
					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
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
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
				else
				{
					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}



//  normal velocity component distribution as initial condition

//	i_middle = 2;
//	i_half = im - 2;
/*
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ 0 ][ j ][ k ] == 0. )
			{
				for ( int i = i_middle; i <= i_half; i++ )					// loop in radial direction, extension for u -velocity component, downwelling here
				{
					m = i_half - i;
					d_i = ( double ) i;
					u.x[ i ][ j ][ k ] = d_i / d_i_half * water_wind / ( double )( k + 1 );	// increase with depth, decrease with distance from coast
					u.x[ m ][ j ][ k ] = d_i / d_i_half * water_wind / ( double )( k + 1 );// decrease with depth, decrease with distance from coast
				}
			}
		}
	}
*/
}





void IC_Thermohalin::IC_v_w_Smoothing ( int iter, Array &h, Array &u, Array &v, Array &w, Array &t, Array &c )
{
// initial conditions for v and w velocity components at the sea surface
// after reading wind data, applying Ekman motion and formulation of west/east coast corretions 

	jmkm = ( double ) ( ( jm -1 ) * ( km -1 ) );

	if ( iter == 0 )
	{
		v_sum = 0.;
		w_sum = 0.;
		t_sum = 0.;
		c_sum = 0.;
		n_smooth = 0;

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v_sum += fabs ( v.x[ im-1 ][ j ][ k ] );
					w_sum += fabs ( w.x[ im-1 ][ j ][ k ] );
	//				t_sum += t.x[ im-1 ][ j ][ k ];
	//				c_sum += c.x[ im-1 ][ j ][ k ];
					n_smooth++;
				}
				else
				{
					v.x[ im-1 ][ j ][ k ] = 0.;
					w.x[ im-1 ][ j ][ k ] = 0.;
				}
			}
		}


		v_sum = v_sum / ( double ) n_smooth;
		w_sum = w_sum / ( double ) n_smooth;
		t_sum = t_sum / ( double ) n_smooth;
		c_sum = c_sum / ( double ) n_smooth;


		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					if ( fabs ( v.x[ im-1 ][ j ][ k ] ) > v_sum )			v.x[ im-1 ][ j ][ k ] = 0.;
					if ( fabs ( w.x[ im-1 ][ j ][ k ] ) > w_sum )			w.x[ im-1 ][ j ][ k ] = 0.;
	//				if ( fabs ( t.x[ im-1 ][ j ][ k ] ) > t_sum )				t.x[ im-1 ][ j ][ k ] = t_sum;
	//				if ( fabs ( c.x[ im-1 ][ j ][ k ] ) > c_sum )			c.x[ im-1 ][ j ][ k ] = c_sum;
				}
				else
				{
					v.x[ im-1 ][ j ][ k ] = 0.;
					w.x[ im-1 ][ j ][ k ] = 0.;
				}
			}
		}
	}



/*
	v_sum = 0.;
	w_sum = 0.;
	t_sum = 0.;
	c_sum = 0.;
	n_smooth = 0;


	for ( int i = 1; i < im-2; i++ )
	{
		for ( int j = 1; j < jm-2; j++ )
		{
			for ( int k = 1; k < km-2; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v_sum += fabs ( v.x[ i ][ j ][ k ] );
					w_sum += fabs ( w.x[ i ][ j ][ k ] );
					t_sum += t.x[ i ][ j ][ k ];
					c_sum += c.x[ i ][ j ][ k ];
					n_smooth++;
				}
			}
		}
	}


	u_sum = u_sum / ( double ) n_smooth;
	v_sum = v_sum / ( double ) n_smooth;
	w_sum = w_sum / ( double ) n_smooth;
	t_sum = t_sum / ( double ) n_smooth;
	c_sum = c_sum / ( double ) n_smooth;


	for ( int i = 1; i < im-2; i++ )
	{
		for ( int j = 1; j < jm-2; j++ )
		{
			for ( int k = 1; k < km-2; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					if ( fabs ( u.x[ i ][ j ][ k ] ) > u_sum )			u.x[ i ][ j ][ k ] = ( u.x[ i+1 ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] + u.x[ i ][ j+1 ][ k ] + u.x[ i ][ j-1 ][ k ] + u.x[ i ][ j ][ k+1 ] + u.x[ i ][ j ][ k-1 ] ) / 6.;
					if ( fabs ( v.x[ i ][ j ][ k ] ) > v_sum )			v.x[ i ][ j ][ k ] = ( v.x[ i+1 ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] + v.x[ i ][ j+1 ][ k ] + v.x[ i ][ j-1 ][ k ] + v.x[ i ][ j ][ k+1 ] + v.x[ i ][ j ][ k-1 ] ) / 6.;
					if ( fabs ( w.x[ i ][ j ][ k ] ) > w_sum )		w.x[ i ][ j ][ k ] = ( w.x[ i+1 ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] + w.x[ i ][ j+1 ][ k ] + w.x[ i ][ j-1 ][ k ] + w.x[ i ][ j ][ k+1 ] + w.x[ i ][ j ][ k-1 ] ) / 6.;
					if ( t.x[ i ][ j ][ k ] > t_sum )						t.x[ i ][ j ][ k ] = ( t.x[ i+1 ][ j ][ k ] + t.x[ i-1 ][ j ][ k ] + t.x[ i ][ j+1 ][ k ] + t.x[ i ][ j-1 ][ k ] + t.x[ i ][ j ][ k+1 ] + t.x[ i ][ j ][ k-1 ] ) / 6.;
					if ( c.x[ i ][ j ][ k ] > c_sum )						c.x[ i ][ j ][ k ] = ( c.x[ i+1 ][ j ][ k ] + c.x[ i-1 ][ j ][ k ] + c.x[ i ][ j+1 ][ k ] + c.x[ i ][ j-1 ][ k ] + c.x[ i ][ j ][ k+1 ] + c.x[ i ][ j ][ k-1 ] ) / 6.;
				}
			}
		}
	}
*/
}




void IC_Thermohalin::IC_v_w_WestEastCoast ( Array &h, Array &u, Array &v, Array &w, Array &c )
{
// initial conditions for v and w velocity components at the sea surface close to east or west coasts
// reversal of v velocity component between north and south equatorial current ommitted at respectively 10°
// w component unchanged

// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: east coast

	k_grad = 6;																			// extension of velocity change
	k_a = k_grad;																		// left distance
	k_b = 0;																				// right distance

	k_water = 0;																		// on water closest to coast
	k_sequel = 1;																		// on solid ground

	for ( int j = 0; j < 91; j++ )													// outer loop: latitude
	{
		for ( int k = 0; k < km - k_grad; k++ )											// inner loop: longitude
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. ) k_sequel = 0;				// if solid ground: k_sequel = 0

			if ( ( h.x[ i_max ][ j ][ k ] == 0. ) && ( k_sequel == 0 ) ) k_water = 0;	// if water and and k_sequel = 0 then is water closest to coast
			else k_water = 1;														// somewhere on water

			if ( ( h.x[ i_max ][ j ][ k ] == 0. ) && ( k_water == 0 ) )	// if water is closest to coast, change of velocity components begins
			{
				for ( int l = 0; l < k_grad; l++ )								// extension of change, sign change in v-velocity and distribution of u-velocity with depth
				{
					v.x[ i_max ][ j ][ k + l ] = - v.x[ i_max ][ j ][ k + l ];	// existing velocity changes sign

					for ( int i = i_middle; i <= i_half; i++ )					// loop in radial direction, extension for u -velocity component, downwelling here
					{
						m = i_half - i;
						d_i = ( double ) i;
						c.x[ i ][ j ][ k ] = 1.11;
						u.x[ i ][ j ][ k + l ] = - d_i / d_i_half * IC_water / ( ( double )( l + 1 ) );			// increase with depth, decrease with distance from coast
						u.x[ m ][ j ][ k + l ] = - d_i / d_i_half * IC_water / ( ( double )( l + 1 ) );			// decrease with depth, decrease with distance from coast
					}
				}

/*
				for ( int l = ( k + k_grad - k_a ); l < ( k + k_grad + k_b + 1 ); l++ ) // starting at local longitude + max extension - begin of smoothing k_a  until ending at  + k_b
				{
					v.x[ i_max ][ j ][ l ] = ( v.x[ i_max ][ j ][ k + k_grad + k_b ] - v.x[ i_max ][ j ][ k + k_grad - k_a ] ) / ( double )( ( k + k_grad + k_b ) -  ( k + k_grad - k_a ) ) * ( double )( l -  ( k + k_grad - k_a ) ) + v.x[ i_max ][ j ][ k + k_grad - k_a ]; // extension of v-velocity, smoothing algorithm by a linear equation 
				}

				for ( int l = k; l < ( k + k_grad + k_b + 1 ); l++ )		// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
				{
					w.x[ i_max ][ j ][ l ] = w.x[ i_max ][ j ][ k + k_grad + k_b ]  / ( double )( ( k + k_grad + k_b ) -  k ) * ( double )( l - k ); // extension of v-velocity
				}
*/
				k_sequel = 1;															// looking for another east coast
			}
		}																						// end of longitudinal loop
		k_water = 0;																	// starting at another latitude
	}																							// end of latitudinal loop




// southern hemisphere: east coast

	k_water = 0;
	k_sequel = 1;

	for ( int j = 91; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. ) k_sequel = 0;

			if ( ( h.x[ i_max ][ j ][ k ] == 0. ) && ( k_sequel == 0 ) ) k_water = 0;
			else k_water = 1;

			if ( ( h.x[ i_max ][ j ][ k ] == 0. ) && ( k_water == 0 ) )
			{
				for ( int l = 0; l < k_grad; l++ )
				{
					v.x[ i_max ][ j ][ k + l ] = - v.x[ i_max ][ j ][ k + l ];

					for ( int i = i_middle; i <= i_half; i++ )
					{
						m = i_half - i;
						d_i = ( double ) i;
						c.x[ i ][ j ][ k ] = 1.11;
						u.x[ i ][ j ][ k + l ] = - d_i / d_i_half * IC_water / ( ( double )( l + 1 ) );			// increase with depth, decrease with distance from coast
						u.x[ m ][ j ][ k + l ] = - d_i / d_i_half * IC_water / ( ( double )( l + 1 ) );			// decrease with depth, decrease with distance from coast
					}
				}
/*
				for ( int l = ( k + k_grad - k_a ); l < ( k + k_grad + k_b + 1 ); l++ )
				{
					v.x[ i_max ][ j ][ l ] = ( v.x[ i_max ][ j ][ k + k_grad + k_b ] - v.x[ i_max ][ j ][ k + k_grad - k_a ] ) / ( double )( ( k + k_grad + k_b ) -  ( k + k_grad - k_a ) ) * ( double )( l -  ( k + k_grad - k_a ) ) + v.x[ i_max ][ j ][ k + k_grad - k_a ];
				}

				for ( int l = k; l < ( k + k_grad + k_b + 1 ); l++ )
				{
					w.x[ i_max ][ j ][ l ] = w.x[ i_max ][ j ][ k + k_grad + k_b ]  / ( double )( ( k + k_grad + k_b ) -  k ) * ( double )( l - k );
				}
*/
				k_sequel = 1;
			}
		}
		k_water = 0;
	}




// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: west coast

	k_grad = 6;																			// extension of velocity change
	k_a = 0;																				// left distance

	k_water = 0;																		// somewhere on water
	flip = 0;																				// somewhere on water

	for ( int j = 0; j < 91; j++ )													// outer loop: latitude
	{
		for ( int k = k_grad; k < km; k++ )											// inner loop: longitude
		{
			if ( h.x[ i_max ][ j ][ k ] == 0. )										// if somewhere on water
			{
				k_water = 0;															// somewhere on water: k_water = 0
				flip = 0;																	// somewhere on water: flip = 0
			}
			else k_water = 1;														// first time on land

			if ( ( flip == 0 ) && ( k_water == 1 ) )							// on water closest to land
			{
				for ( int l = k; l > ( k - k_grad - 1 ); l-- )						// backward extention of velocity change: nothing changes
				{
					w.x[ i_max ][ j ][ l ] = - w.x[ i_max ][ j ][ l ];

					for ( int i = i_middle; i <= i_half; i++ )					// loop in radial direction, extension for u -velocity component, downwelling here
					{
						m = i_half - i;
						d_i = ( double ) i;
						c.x[ i ][ j ][ k ] = 1.11;
						u.x[ i ][ j ][ l ] = + d_i / d_i_half * IC_water / ( ( double )( k - l + 1 ) );			// increase with depth, decrease with distance from coast
						u.x[ m ][ j ][ l ] = + d_i / d_i_half * IC_water / ( ( double )( k - l + 1 ) );			// decrease with depth, decrease with distance from coast
					}
				}
/*
				for ( int l = k; l > ( k - k_grad - k_a + 1 ); l-- )			// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
				{
					v.x[ i_max ][ j ][ l ] = v.x[ i_max ][ j ][ k - k_grad - k_a ] / ( double )( ( k - k_grad - k_a ) - k ) * ( double )( l - k ); // extension of v-velocity
				}

				for ( int l = ( k - k_grad - 3 ); l < ( k - k_grad + 3 ); l++ )			// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
				{
					w.x[ i_max ][ j ][ l ] = (  - w.x[ i_max ][ j ][ k - k_grad - 3 ] + w.x[ i_max ][ j ][ k - k_grad + 3 ] ) * ( double ) ( l - ( k - k_grad - 3 ) ) / ( double ) ( ( k - k_grad + 3 ) - ( k - k_grad - 3 ) ) - w.x[ i_max ][ j ][ k - k_grad + 3 ];
				}
*/
				flip = 1;
			}
		}
		flip = 0;
	}






// southern hemisphere: west coast

	k_water = 0;
	flip = 0;

	for ( int j = 91; j < jm; j++ )
	{
		for ( int k = k_grad; k < km; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 0. )
			{
				k_water = 0;
				flip = 0;
			}
			else k_water = 1;

			if ( ( flip == 0 ) && ( k_water == 1 ) )
			{
				for ( int l = k; l > ( k - k_grad + 1 ); l-- )
				{
					w.x[ i_max ][ j ][ l ] = - w.x[ i_max ][ j ][ l ];

					for ( int i = i_middle; i <= i_half; i++ )
					{
						m = i_half - i;
						d_i = ( double ) i;
						c.x[ i ][ j ][ k ] = 1.11;
						u.x[ i ][ j ][ l ] = + d_i / d_i_half * IC_water / ( ( double )( k - l + 1 ) );
						u.x[ m ][ j ][ l ] = + d_i / d_i_half * IC_water / ( ( double )( k - l + 1 ) );
					}
				}
/*
				for ( int l = k; l > ( k - k_grad - k_a - 1 ); l-- )
				{
					v.x[ i_max ][ j ][ l ] = v.x[ i_max ][ j ][ k - k_grad - k_a ] / ( double )( ( k - k_grad - k_a ) - k ) * ( double )( l - k );
				}

				for ( int l = ( k - k_grad - 3 ); l < ( k - k_grad + 3 ); l++ )			// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
				{
					w.x[ i_max ][ j ][ l ] = ( - w.x[ i_max ][ j ][ k - k_grad - 3 ] + w.x[ i_max ][ j ][ k - k_grad + 3 ] ) * ( double ) ( l - ( k - k_grad - 3 ) ) / ( double ) ( ( k - k_grad + 3 ) - ( k - k_grad - 3 ) ) - w.x[ i_max ][ j ][ k - k_grad + 3 ];
				}
*/
				flip = 1;
			}
		}
		flip = 0;
	}

}






void IC_Thermohalin::IC_v_w_Ekman ( Array &h, Array &v, Array &w )
{
//	initial conditions for v and w velocity components at the sea surface
//	reduction of velocity with increasing depth for the purpose of the Ekman spiral

//	north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 0; j < 31; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					vel_magnitude = sqrt ( v.x[ i ][ j ][ k ] * v.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k ] * w.x[ i ][ j ][ k ] );

					if ( w.x[ i ][ j ][ k ] == 0. ) w.x[ i ][ j ][ k ] = 1.e-6;
					if ( v.x[ i ][ j ][ k ] == 0. ) v.x[ i ][ j ][ k ] = 1.e-6;
					if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

					alfa = asin ( fabs ( w.x[ i ][ j ][ k ] ) / vel_magnitude );
					Ekman = ( double ) ( i_max - i ) * Ekman_angle_add;

					if ( w.x[ im-1 ][ j ][ k ] > 0. )
					{
						angle = + alfa - Ekman_angle - Ekman;
					}
					else
					{
						angle = - alfa - Ekman_angle - Ekman;
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle > 0. ) )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle < 0. ) )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if (  w.x[ im-1 ][ j ][ k ] < 0. )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}
				}
				else
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}



// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 31; j < 60; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					vel_magnitude = sqrt ( v.x[ i ][ j ][ k ] * v.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k ] * w.x[ i ][ j ][ k ] );

					if ( w.x[ i ][ j ][ k ] == 0. ) w.x[ i ][ j ][ k ] = 1.e-6;
					if ( v.x[ i ][ j ][ k ] == 0. ) v.x[ i ][ j ][ k ] = 1.e-6;
					if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

					alfa = asin ( fabs ( w.x[ i ][ j ][ k ] ) / vel_magnitude );
					Ekman = ( double ) ( i_max - i ) * Ekman_angle_add;

					if ( v.x[ im-1 ][ j ][ k ] > 0. )
					{
						angle = + alfa - Ekman_angle - Ekman;
					}
					else
					{
						angle = - alfa - Ekman_angle - Ekman;
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle > 0. ) )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle < 0. ) )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = - vel_magnitude * sin ( angle );
					}

					if (  v.x[ i ][ j ][ k ] < 0. )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = - vel_magnitude * sin ( angle );
					}
				}
				else
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}




// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 60; j < 91; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					vel_magnitude = sqrt ( v.x[ i ][ j ][ k ] * v.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k ] * w.x[ i ][ j ][ k ] );

					if ( w.x[ i ][ j ][ k ] == 0. ) w.x[ i ][ j ][ k ] = 1.e-6;
					if ( v.x[ i ][ j ][ k ] == 0. ) v.x[ i ][ j ][ k ] = 1.e-6;
					if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

					alfa = asin ( fabs ( w.x[ i ][ j ][ k ] ) / vel_magnitude );
					Ekman = ( double ) ( i_max - i ) * Ekman_angle_add;

					if ( w.x[ im-1 ][ j ][ k ] > 0. )
					{
						angle = + alfa - Ekman_angle - Ekman;
					}
					else
					{
						angle = - alfa - Ekman_angle - Ekman;
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle > 0. ) )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle < 0. ) )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if (  w.x[ im-1 ][ j ][ k ] < 0. )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}
				}
				else
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}




// south equatorial Hadley cell ( from j=90 till j=120 compares to 0° till 30° )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 91; j < 122; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					vel_magnitude = sqrt ( v.x[ i ][ j ][ k ] * v.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k ] * w.x[ i ][ j ][ k ] );

					if ( w.x[ i ][ j ][ k ] == 0. ) w.x[ i ][ j ][ k ] = 1.e-6;
					if ( v.x[ i ][ j ][ k ] == 0. ) v.x[ i ][ j ][ k ] = 1.e-6;
					if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

					beta = asin ( fabs( w.x[ i ][ j ][ k ] ) / vel_magnitude );
					Ekman = ( double ) ( i_max - i ) * Ekman_angle_add;

					if ( w.x[ im-1 ][ j ][ k ] > 0. )
					{
						angle = beta - Ekman_angle - Ekman;
					}
					else
					{
						angle = - beta - Ekman_angle - Ekman;
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle > 0. ) )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle < 0. ) )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if (  w.x[ im-1 ][ j ][ k ] < 0. )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}
				}
				else
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}


// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 121; j < 151; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					vel_magnitude = sqrt ( v.x[ i ][ j ][ k ] * v.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k ] * w.x[ i ][ j ][ k ] );

					if ( w.x[ i ][ j ][ k ] == 0. ) w.x[ i ][ j ][ k ] = 1.e-6;
					if ( v.x[ i ][ j ][ k ] == 0. ) v.x[ i ][ j ][ k ] = 1.e-6;
					if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

					beta = asin ( fabs( w.x[ i ][ j ][ k ] ) / vel_magnitude );
					Ekman = ( double ) ( i_max - i ) * Ekman_angle_add;

					if ( v.x[ im-1 ][ j ][ k ] < 0. )
					{
						angle = beta - Ekman_angle - Ekman;
					}
					else
					{
						angle = - beta - Ekman_angle - Ekman;
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle > 0. ) )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle < 0. ) )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = - vel_magnitude * sin ( angle );
					}

					if (  v.x[ i ][ j ][ k ] > 0. )
					{
						v.x[ i ][ j ][ k ] = + vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = - vel_magnitude * sin ( angle );
					}
				}
				else
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}


// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 151; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					vel_magnitude = sqrt ( v.x[ i ][ j ][ k ] * v.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k ] * w.x[ i ][ j ][ k ] );

					if ( w.x[ i ][ j ][ k ] == 0. ) w.x[ i ][ j ][ k ] = 1.e-6;
					if ( v.x[ i ][ j ][ k ] == 0. ) v.x[ i ][ j ][ k ] = 1.e-6;
					if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

					beta = asin ( fabs( w.x[ i ][ j ][ k ] ) / vel_magnitude );
					Ekman = ( double ) ( i_max - i ) * Ekman_angle_add;

					if ( w.x[ im-1 ][ j ][ k ] > 0. )
					{
						angle = beta - Ekman_angle - Ekman;
					}
					else
					{
						angle = - beta - Ekman_angle - Ekman;
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle > 0. ) )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if ( ( w.x[ im-1 ][ j ][ k ] > 0. ) && ( angle < 0. ) )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}

					if (  w.x[ im-1 ][ j ][ k ] < 0. )
					{
						v.x[ i ][ j ][ k ] = - vel_magnitude * cos ( angle );
						w.x[ i ][ j ][ k ] = + vel_magnitude * sin ( angle );
					}
				}
				else
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}

}







void IC_Thermohalin::BC_Temperature_Salinity ( int Ma, int Ma_max, int Ma_max_half, double t_0, double p_0, double c_0, double t_Cretaceous_max, double t_Average, double t_equator, double t_pole, double ua, double va, double wa, double ta, double ca, double pa, Array_2D &t_j, Array_2D & c_j, Array_2D &p_j, Array &h, Array &t, Array &c, Array &p )
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
	k_half = ( km -1 ) / 2;


	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

	t_coeff = ( t_pole - t_equator );
	t_Cretaceous_coeff = t_Cretaceous_max / ( ( double ) Ma_max_half - ( double ) ( Ma_max_half * Ma_max_half / Ma_max ) );
	t_Cretaceous = t_Cretaceous_coeff * ( double ) ( - ( Ma * Ma ) / Ma_max + Ma );

	cout.precision ( 3 );

	time_slice_comment = "      time slice of Cretaceous-AGCM:";
	time_slice_number = " Ma = ";
	time_slice_unit = " million years";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << time_slice_comment << resetiosflags ( ios::left ) << setw ( 6 ) << fixed << setfill ( ' ' ) << time_slice_number << setw ( 3 ) << Ma << setw ( 12 ) << time_slice_unit << endl << endl;


	temperature_comment = "      temperature increase at cretaceous times: ";
	temperature_gain = " t cretaceous";
	temperature_modern = "      mean temperature at modern times: ";
	temperature_average = " t modern";
	temperature_unit =  "°C ";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << temperature_comment << resetiosflags ( ios::left ) << setw ( 12 ) << temperature_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << t_Cretaceous << setw ( 5 ) << temperature_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_modern << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_Average << setw ( 5 ) << temperature_unit << endl;


	c_Average = ( ( t_Average + 346. ) / 10. );
	c_Cretaceous = ( ( t_Average + t_Cretaceous + 346. ) / 10. );
	c_diff = c_Cretaceous - c_Average;

	salinity_comment = "      salinity increase at cretaceous times: ";
	salinity_gain = " c cretaceous";
	salinity_modern = "      mean salinity at modern times: ";
	salinity_average = " c modern";
	salinity_unit =  "psu ";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << salinity_comment << resetiosflags ( ios::left ) << setw ( 12 ) << salinity_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << c_diff << setw ( 5 ) << salinity_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << salinity_modern << resetiosflags ( ios::left ) << setw ( 13 ) << salinity_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << c_Average << setw ( 5 ) << salinity_unit << endl;
	cout << endl;


	t_Cretaceous = ( t_Cretaceous + t_Average + t_0 ) / t_0 - ( ( t_Average + t_0 ) / t_0 );    // non-dimensional


	if (  Ma > 0 )																			// when time slices are run, modern world excluded
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				if ( h.x[ im-1 ][ j ][ k ] == 0. )
				{
					d_j = ( double ) j;
					t.x[ im-1 ][ j ][ k ] = t_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + t_pole + t_Cretaceous;
					t_Celsius = t.x[ im-1 ][ j ][ k ] * t_0 - t_0;
					if ( t_Celsius <= 0. ) t_Celsius = 4.;						// water temperature not below 4°C
					c.x[ im-1 ][ j ][ k ] = ( ( t_Celsius + 346. ) / 10. ) / c_0;
				}
				else
				{
					t.x[ im-1 ][ j ][ k ] = ta;
					c.x[ im-1 ][ j ][ k ] = ca;
				}
			}
		}
	}                                      // end of while (  Ma > 0 ) when time slices are run, modern world excluded



// distribution of t and c with increasing depth till i_beg valid for all time slices including the modern world

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					t.x[ i ][ j ][ k ] = ( t.x[ im-1 ][ j ][ k ] - ta ) * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg ) + ta;
					t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;
					if ( t_Celsius <= 0. ) t_Celsius = 4.;						// water temperature not below 4°C
					c.x[ i ][ j ][ k ] = ( ( t_Celsius + 346. ) / 10. ) / c_0;
				}
				else
				{
					t.x[ i ][ j ][ k ] = ta;
					c.x[ i ][ j ][ k ] = ca;
				}
			}
		}
	}



// distribution of t and c from sea bottom to thermocline valid for all time slices including the modern world

	for ( int i = 0; i < i_beg; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					t.x[ i ][ j ][ k ] = ta;
					t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;
					if ( t_Celsius <= 0. ) t_Celsius = 4.;						// water temperature not below 4°C
					c.x[ i ][ j ][ k ] = ( ( t_Celsius + 346. ) / 10. ) / c_0;
				}
				else
				{
					t.x[ i ][ j ][ k ] = ta;
					c.x[ i ][ j ][ k ] = ca;
				}
			}
		}
	}
}







void IC_Thermohalin::IC_Atlantischer_Ozean ( Array &h, Array &u, Array &v, Array &w, Array &c )
{
// Ströme entlang der Küsten
// Schliessen der polaren, subpolaren und subtropischen atmosphärischen Wirbelsysteme


// Thermohaline Conveyor Belt

// Atlantic

// Südatlantik		diagonal from Kap Agulhas in Südafrika until Kap St. Roque in Südamerika

// Südäquatorial-Strom als Oberflächenströmung ( from j=96 until j=112 compares to 6°S until 22°S,
//                                               from k=325 until k=km compares to 35°W until 0° )

	j_beg = 96;
	j_end = 113;
	k_beg = 325;
	k_end = km;
	j_run = 0;
	k_run = 0;
	j_step = 2;
	k_step = 10;


	while ( ( j_beg + j_run ) <= j_end && ( k_beg + k_run + k_step ) <= ( k_end + k_step ) )
	{
		for ( int j = j_beg + j_run; j < ( j_beg + j_step + j_run ); j++ )
		{
			for ( int k = k_beg + k_run; k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		k_run++;
		}
	j_run++;
	}



// westlich des konstanten Geschwindigkeitsstreifens diagonal from Kap Agulhas in Südafrika until Kap St. Roque in Südamerika
// Südäquatorial-Strom als Oberflächenströmung ( from j=96 until j=112 compares to 6°S until 22°S,
//                                               from k=325 until k=km compares to 35°W until 0° )

	j_beg = 96;
	j_end = 113;
	k_beg = 325;
	k_end = km;
	j_run = 0;
	k_run = 0;
	j_step = 2;
	k_step = 10;
	k_exp = 10;

	while ( ( j_beg + j_run ) <= j_end && ( k_beg + k_run - k_exp ) <= ( k_end - k_exp ) )
	{
		for ( int j = j_beg + j_run; j < ( j_beg + j_step + j_run ); j++ )
		{
			for ( int k = ( k_beg + k_run - k_exp ); k < ( k_beg + k_run ); k++ )
			{
				k_z = k - ( k_beg + k_run - k_exp );
				k_n = ( k_beg + k_run ) - ( k_beg + k_run - k_exp );

				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_beg + k_run ] - v.x[ i ][ j ][ k_beg + k_run - k_exp ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg + k_run - k_exp ];
						w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_beg + k_run ] - w.x[ i ][ j ][ k_beg + k_run - k_exp ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg + k_run - k_exp ];
					}
				}
			}
			k_run++;
		}
	j_run++;
	}



// östlich des konstanten Geschwindigkeitsstreifens diagonal from Kap Agulhas in Südafrika until Kap St. Roque in Südamerika
// Südäquatorial-Strom als Oberflächenströmung ( from j=96 until j=112 compares to 6°S until 22°S,
//                                                                              from k=325 until k=km compares to 35°W until 0° )

	j_beg = 96;
	j_end = 113;
	k_beg = 325;
	k_end = km;
	j_run = 0;
	k_run = 0;
	j_step = 2;
	k_step = 10;
	k_exp = 3;

	while ( ( j_beg + j_run ) <= j_end && ( k_beg + k_run + k_step + k_exp ) <= ( k_end + k_exp ) )
	{
		for ( int j = j_beg + j_run; j < ( j_beg + j_step + j_run ); j++ )
		{
			for ( int k = ( k_beg + k_run + k_step ); k < ( k_beg + k_run + k_step + k_exp ); k++ )
			{
				if ( k >= k_end ) break;
				k_z = k - ( k_beg + k_run + k_step );
				k_n = ( k_beg + k_run + k_step + k_exp ) - ( k_beg + k_run + k_step ) ;

				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_beg + k_run + k_step + k_exp ] - v.x[ i ][ j ][ k_beg + k_run + k_step ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg + k_run + k_step ];
						w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_beg + k_run + k_step + k_exp ] - w.x[ i ][ j ][ k_beg + k_run + k_step ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg + k_run + k_step ];
					}
				}
			}
			k_run++;
		}
	j_run++;
	}







// Thermohalin Conveyor Belt

// Atlantic

// Südatlantik		diagonal from Kap Agulhas in Südafrika until Kap St. Roque in Südamerika

// Südäquatorial-Strom als Oberflächenströmung ( from j=117 until j=128 compares to 27°S until 38°S,
//                                                                               from k=0 until k=20 compares to 0° until 20°O )

	j_beg = 117;
	j_end = 129;
	k_beg = 0;
	k_end = 21;
	j_run = 0;
	k_run = 0;
	j_step = 2;
	k_step = 10;


	while ( ( j_beg + j_run ) <= j_end && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = j_beg + j_run; j < ( j_beg + j_step + j_run ); j++ )
		{
			for ( int k = k_beg + k_run; k < ( k_beg + k_step + k_run ); k++ )
			{
				if ( k >= k_end ) break;
				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		k_run++;
		}
	j_run++;
	}




// westlich des konstanten Geschwindigkeitsstreifens diagonal from Kap Agulhas in Südafrika until Kap St. Roque in Südamerika
// Südäquatorial-Strom als Oberflächenströmung ( from j=117 until j=134 compares to 27°S until 44°S,
//                                                                               from k=0 until k=20 compares to 0° until 20°O )

	j_beg = 117;
	j_end = 135;
	k_beg = 0;
	k_end = 21;
	j_run = 0;
	k_run = 0;
	j_step = 2;
	k_step = 10;
	k_exp = 10;

	while ( ( j_beg + j_run ) <= j_end && ( k_beg + k_run - k_exp ) <= ( k_end - k_exp ) )
	{
		for ( int j = j_beg + j_run; j < ( j_beg + j_step + j_run ); j++ )
		{
			for ( int k = ( k_beg + k_run - k_exp ); k < ( k_beg + k_run ); k++ )
			{
				if ( k >= k_end ) break;
				k_z = k - ( k_beg + k_run - k_exp );
				k_n = ( k_beg + k_run ) - ( k_beg + k_run - k_exp );

				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_beg + k_run ] - v.x[ i ][ j ][ k_beg + k_run - k_exp ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg + k_run - k_exp ];
						w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_beg + k_run ] - w.x[ i ][ j ][ k_beg + k_run - k_exp ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg + k_run - k_exp ];
					}
				}
			}
			k_run++;
		}
	j_run++;
	}



// östlich des konstanten Geschwindigkeitsstreifens diagonal from Kap Agulhas in Südafrika until Kap St. Roque in Südamerika
// Südäquatorial-Strom als Oberflächenströmung ( from j=117 until j=134 compares to 27°S until 44°S,
//                                                                               from k=0 until k=20 compares to 0° until 20°O )

	j_beg = 117;
	j_end = 113;
	k_beg = 0;
	k_end = 21;
	j_run = 0;
	k_run = - 1;
	j_step = 2;
	k_step = 0;
	k_exp = 6;

	while ( ( j_beg + j_run ) <= j_end && ( k_beg + k_run + k_step + k_exp ) <= ( k_end + k_exp ) )
	{
		for ( int j = j_beg + j_run; j < ( j_beg + j_step + j_run ); j++ )
		{
			for ( int k = ( k_beg + k_run + k_step ); k < ( k_beg + k_run + k_step + k_exp ); k++ )
			{
				if ( k >= k_end ) break;
				k_z = k - ( k_beg + k_run + k_step );
				k_n = ( k_beg + k_run + k_step + k_exp ) - ( k_beg + k_run + k_step ) ;

				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_beg + k_run + k_step + k_exp ] - v.x[ i ][ j ][ k_beg + k_run + k_step ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg + k_run + k_step ];
						w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_beg + k_run + k_step + k_exp ] - w.x[ i ][ j ][ k_beg + k_run + k_step ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg + k_run + k_step ];
					}
				}
			}
			k_run++;
		}
	j_run++;
	}






// Atlantic
// Küstenströmungen

// Südamerika		Nordküsten

// Guayana-Strom im Nord-Osten from Südamerika ( from j=71 until j=96 compares to 6°S until 19°N,
//                                                                                  from k=280 until k=325 compares to 80°W until 35°W )

	j_beg = 71;
	j_end = 97;
	k_beg = 280;
	k_end = 326;

	k_step = 30;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;

		for ( int k = k_b + 8; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 8 ] ) * ( double ) ( k - ( k_b + 8 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 8 ) ) + v.x[ i ][ j ][ k_b + 8 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 8 ] ) * ( double ) ( k - ( k_b + 8 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 8 ) ) + w.x[ i ][ j ][ k_b + 8 ];
			}
		}

	}





// Atlantischer Ozean

// Golf from Mexiko                                     ( from j=53 until j=58 compares to 23°N until 27°N,
// parallel zum Äquator ostwärts                 from k=234 until k=252 compares to 77°W until 96°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 53; j < 59; j++ )
		{
			for ( int k = 234; k < 253; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						v.x[ i ][ j ][ k ] = + 0.0;
						w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		}
	}




// Atlantischer Ozean

// Golf from Mexiko                                          ( from j=63 until j=67 compares to 23°N until 27°N,
// parallel zum Äquator westwärts                    from k=264 until k=283 compares to 77°W until 96°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 63; j < 68; j++ )
		{
			for ( int k = 264; k < 284; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						v.x[ i ][ j ][ k ] = + 0.0;
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		}
	}



// Atlantischer Ozean

// Golf from Mexiko                                          ( from j=63 until j=67 compares to 23°N until 27°N,
// parallel zum Äquator nordwärts                    from k=264 until k=283 compares to 77°W until 96°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 63; j < 68; j++ )
		{
			for ( int k = 264; k < 284; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = - 0.00;
					}
				}
			}
		}
	}





// Atlantic

// Verbindung Ostküste Golf from Mexiko
// Golf from Mexiko                                     ( from j=63 until j=67 compares to 23°N until 27°N,
// parallel zum Äquator Ostküste                 from k=264 until k=283 compares to 77°W until 96°W )

	j_beg = 63;
	j_end = 68;
	k_beg = 264;
	k_end = 284;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
//					u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_end ] - u.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + u.x[ i ][ j ][ k_beg ];
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}





// Atlantic

// Verbindung zwischen Guayana- und Golf from Mexiko-Strom

// Nordamerika		Ostküsten

// Floridastraße an der Südspitze from Florida ( from j=63 until j=67 compares to 23°N until 27°N,
// parallel zum Äquator                         from k=264 until k=283 compares to 77°W until 96°W )
/*
	j_beg = 63;
	j_end = 68
	k_beg = 264;
	k_end = 284;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}
*/





// Atlantic
// Nordamerika		Ostküsten

// Golf-Strom im Nord-Osten from Nordamerika ( from j=40 until j=70 compares to 20°N until 50°N,
// südlich from Neufundland                    from k=278 until k=308 compares to 82°W until 52°W )

	j_beg = 40;
	j_end = 71;
	k_beg = 278;
	k_end = 309;
	k_step = 20;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;

		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}

	}




// Atlantic

// Verbindung zwischen Guayana- und Golf-Strom

// Nordamerika		Ostküsten

// Golf-Strom im Nord-Osten from Nordamerika ( from j=44 until j=70 compares to 20°N until 46°N,
// südlich from Neufundland                    from k=280 until k=308 compares to 80°W until 52°W )

	j_beg = 44;
	j_end = 71;
	k_beg = 280;
	k_end = 309;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}






// Atlantic

// Sankt Lorenz-Golf		Geschwindigkeitsreduktion		Westküste

// Sankt Lorenz-Golf ( from j=38 until j=45 compares to 45°N until 52°N,
//                                 from k=295 until k=305 compares to 65°W until 55°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 38; j < 46; j++ )
		{
			for ( int k = 295; k < 306; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
//					u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}





// Atlantic
// Thermohalin Conveyor Belt

// unterhalb Neufundlands und dem Sankt Lorenz-Strom als Anschlussstück

// Golf-Strom auf dem Atlantik from Süd-West in Nord-Ost-Richtung ( from j=43 until j=51 compares to 39°N until 47°N,
//                                                                                                            from k=295 until k=320 compares to 40°W until 65°W )

	j_beg = 43;
	j_end = 52;
	k_beg = 295;
	k_end = 321;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 50;


	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run ); k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		k_run++;
		}
	j_run++;
	}






// Atlantic
// Thermohalin Conveyor Belt

// Verbindung unterhalb Neufundlands und dem Sankt Lorenz-Strom und dem Nordatlantik

// Golf-Strom auf dem Atlantik from Süd-West in Nord-Ost-Richtung ( from j=48 until j=54 compares to 36°N until 42°N,
//                                                                                                           from k=295 until k=340 compares to 20°W until 65°W )

	j_beg = 48;
	j_end = 55;
	k_beg = 295;
	k_end = 341;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}








// Atlantic
// Thermohalin Conveyor Belt

// diagonal from der Labradorspitze until nach Nordnorwegen

// Golf-Strom auf dem Atlantik from Süd-West in Nord-Ost-Richtung ( from j=37 until j=44 compares to 63°N until 46°N,
// Verlängerung                                                                                     from k=304 until k=km compares to 56°W until 0°W )

	j_beg = 37;
	j_end = 45;
	k_beg = 304;
	k_end = km;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 50;


	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run ); k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		k_run++;
		}
	j_run++;
	}




// Atlantic

// Übergang from Golfstrom zum Nordatlantik

// Grönland-Strom im Osten                                      ( from j=15 until j=27 compares to 75°N until 63°N,
// parallel zum Äquator                                                from k=350 until k=km compares to 10°W until 0° )

	j_beg = 15;
	j_end = 28;
	k_beg = 350;
	k_end = km;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}




// Atlantic

// Verbindung zwischen Atlantik und Nordmeer, Golfstrom Richtung Norwegen

// Nordatlantik
// diagonal from der Labradorspitze until nach Nordnorwegen

// Golf-Strom auf dem Atlantik from Süd-West in Nord-Ost-Richtung ( from j=27 until j=44 compares to 63°N until 46°N,
// Verlängerung                                                                                    from k=304 until k=km compares to 56°W until 0°W )

	j_beg = 27;
	j_end = 45;
	k_beg = 304;
	k_end = km;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}







// Atlantic
// Thermohalin Conveyor Belt

// diagonal als Verlängerung über 0° until nach Nordnorwegen

// Golf-Strom auf dem Atlantik from Süd-West in Nord-Ost-Richtung ( from j=18 until j=30 compares to 72°N until 60°N,
//                                                                                                           from k=0 until k=20 compares to 0°W until 20°O )

	j_beg = 18;
	j_end = 31;
	k_beg = 0;
	k_end = 21;
	k_step = 10;

	k_exp = 10;
	k_w = 4;

	v_grad = 0.0008;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b - k_step - k_exp; k < ( k_b - k_w ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b - k_w ] - v.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + v.x[ i ][ j ][ k_b - k_step - k_exp ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b - k_w ] - w.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + w.x[ i ][ j ][ k_b - k_step - k_exp ];
			}
		}
	}






// Atlantic

// Grönland		Westküsten

// Grönland-Strom im Westen until in die Baffinbucht ( from j=15 until j=31 compares to 75°N until 59°N,
//                                                   from k=300 until k=315 compares to 60°W until 45°W )

	j_beg = 15;
	j_end = 32;
	k_beg = 300;
	k_end = 316;
	k_step = 4;
	k_w = 2;
	k_exp = 4;

	v_grad = + 0.0010;

	k_a = k_b = 0;

	flip = 0;


	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b - k_step - k_exp; k < ( k_b - k_w ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b - k_w ] - v.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + v.x[ i ][ j ][ k_b - k_step - k_exp ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b - k_w ] - w.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + w.x[ i ][ j ][ k_b - k_step - k_exp ];
			}
		}
	}






// Atlantic

// Südspitze from Grönland		Westküsten

// Grönland-Strom im Westen until in die Baffinbucht ( from j=15 until j=31 compares to 75°N until 59°N,
// parallel zum Äquator                                                from k=310 until k=315 compares to 50°W until 45°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 15; j < 32; j++ )
		{
			for ( int k = 310; k < 316; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						v.x[ i ][ j ][ k ] = + 0.0;
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		}
	}






// Atlantic

// Übergang an der Südspitze from Grönland zum Atlantik		Westküsten

// Grönland-Strom im Westen until in die Baffinbucht ( from j=15 until j=31 compares to 75°N until 59°N,
// parallel zum Äquator                                                from k=300 until k=315 compares to 60°W until 45°W )

	j_beg = 15;
	j_end = 32;
	k_beg = 300;
	k_end = 316;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}






// Atlantic

// Grönland		Ostküsten
// Thermohalin Conveyor Belt

// Ost-Grönland-Strom im Nord-Osten from Grönland ( from j=10 until j=31 compares to 80°N until 59°N,
//                                                                                    from k=315 until k=342 compares to 45°W until 18°W )

	j_beg = 10;
	j_end = 32;
	k_beg = 315;
	k_end = 343;
	k_step = 8;
	k_w = 4;
	k_exp = 4;

	v_grad = + 0.0008;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg-3; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
//					u.x[ i ][ j ][ k ] = - IC_water;
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;

		for ( int k = k_b + k_w; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step ] - v.x[ i ][ j ][ k_b + k_w ] ) * ( double ) ( k - ( k_b + k_w ) ) / ( double ) ( ( k_b + k_step ) - ( k_b + k_w ) ) + v.x[ i ][ j ][ k_b + k_w ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step ] - w.x[ i ][ j ][ k_b + k_w ] ) * ( double ) ( k - ( k_b + k_w ) ) / ( double ) ( ( k_b + k_step ) - ( k_b + k_w ) ) + w.x[ i ][ j ][ k_b + k_w ];
			}
		}

	}





// Atlantic

// Island		Südküste
// Thermohalin Conveyor Belt

// Süd-Island-Strom ( from j=26 until j=29 compares to 64°N until 61°N,
//                                 from k=326 until k=347 compares to 34°W until 13°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 26; j < 30; j++ )
		{
			for ( int k = 326; k < 348; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}






// Atlantic

// südlicher Übergang bei Island		Südküste
// Thermohalin Conveyor Belt

// Süd-Island-Strom ( from j=24 until j=29 compares to 64°N until 61°N,
//                                 from k=326 until k=347 compares to 34°W until 13°W )

	j_beg = 24;
	j_end = 30;
	k_beg = 326;
	k_end = 348;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_end ][ k ] - u.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + u.x[ i ][ j_beg ][ k ];
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}






// Atlantic

// westlicher Übergang bei Island		Südküste
// Thermohalin Conveyor Belt

// Süd-Island-Strom ( from j=24 until j=29 compares to 64°N until 61°N,
//                                 from k=326 until k=347 compares to 34°W until 13°W )

	j_beg = 24;
	j_end = 30;
	k_beg = 326;
	k_end = 348;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_end ] - u.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + u.x[ i ][ j ][ k_beg ];
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}





// Atlantic

// östlicher Übergang bei Island		Südküste
// Thermohalin Conveyor Belt

// Süd-Island-Strom ( from j=24 until j=29 compares to 64°N until 61°N,
//                                from k=326 until k=347 compares to 34°W until 13°W )

	j_beg = 24;
	j_end = 30;
	k_beg = 326;
	k_end = 348;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_end ] - u.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + u.x[ i ][ j ][ k_beg ];
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}








// Atlantic
// Thermohalin Conveyor Belt

// Vorbereitung zur Abströmung des Golf-Stroms südlich from Island

// Golf-Strom auf dem Atlantik nach Süden ( from j=32 until j=40 compares to 58°N until 50°N,
//                                                                     from k=317 until k=327 compares to 43°W until 33°W )

	for ( int i = i_bottom; i < im-1; i++ )
	{
		for ( int j = 32; j < 41; j++ )
		{
			for ( int k = 317; k < 328; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = c.x[ im-1 ][ j ][ k ];
					u.x[ i ][ j ][ k ] = - IC_water;
//					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
//					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					v.x[ i ][ j ][ k ] = + IC_water;
					w.x[ i ][ j ][ k ] = + IC_water;
				}
			}
		}
	}






// Atlantic

// Südamerika		Ostküsten

// Brasil-Strom im Süd-Osten from Südamerika ( from j=98 until j=127 compares to 8°S until 37°S,
// Oberflächenströmung                                       from k=305 until k=325 compares to 55°W until 35°W )

	j_beg = 98;
	j_end = 128;
	k_beg = 305;
	k_end = 326;
	k_step = 20;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}

		for ( int k = k_b; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}

		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}


		k_a = k_b;
		flip = 0;
	}






// Atlantic

// Südamerika		Ostküsten

// Falkland-Strom im Süd-Osten from Südamerika ( from j=127 until j=146 compares to 37°S until 56°S,
//                                                                               from k=290 until k=330 compares to 70°W until 30°W )

	j_beg = 127;
	j_end = 147;
	k_beg = 290;
	k_end = 331;
	k_step = 20;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}

	}






// Verbindung zwischen Brasil- und Falklandstrom

// südlich vom Brasil-Strom im Süd-Osten from Südamerika ( from j=98 until j=127 compares to 8°S until 37°S,
// Oberflächenströmung  mit Falklandstrom                            from k=300 until k=325 compares to 60°W until 35°W )

	j_beg = 198;
	j_end = 128;
	k_beg = 300;
	k_end = 326;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}







// Atlantic

// Nordafrika		Westküsten

// Kanaren-Strom im Westen from Nordafrika ( from j=50 until j=77 compares to 13°N until 40°N,
//                                                                        from k=340 until k=352 compares to 20°W until 8°W )

	j_beg = 50;
	j_end = 78;
	k_beg = 340;
	k_end = 353;

	v_grad = + 0.001;
	k_step = 8;
	k_w = 4;
	k_exp = 4;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b - k_step - k_exp; k < ( k_b - k_w ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b - k_w ] - v.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + v.x[ i ][ j ][ k_b - k_step - k_exp ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b - k_w ] - w.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + w.x[ i ][ j ][ k_b - k_step - k_exp ];
			}
		}

	}





// Atlantic

// Nordafrika		Westküsten

// Guinea-Strom im Westen from Afrika ( from j=77 until j=89 compares to 13°N until 1°N,
//                                                              from k=340 until k=352 compares to 20°W until 8°W )

	j_beg = 77;
	j_end = 90;
	k_beg = 340;
	k_end = 353;

	v_grad = + 0.001;
	k_step = 8;
	k_w = 4;
	k_exp = 4;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b - k_step - k_exp; k < ( k_b - k_w ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b - k_w ] - v.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + v.x[ i ][ j ][ k_b - k_step - k_exp ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b - k_w ] - w.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + w.x[ i ][ j ][ k_b - k_step - k_exp ];
			}
		}

	}






// Atlantic

// Nordafrika		Westküsten
// Verbindung zwischen Kanaren- und Guineastrom

// Kanaren-Strom im Westen from Nordafrika ( from j=50 until j=77 compares to 13°N until 40°N,
//                                                                        from k=340 until k=352 compares to 20°W until 8°W )

	j_beg = 50;
	j_end = 78;
	k_beg = 340;
	k_end = 353;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}


// Atlantic

// Nordafrika		Westküsten
// Verbindung zwischen Kanarenstrom und Nordatlantik

// Kanaren-Strom im Westen from Nordafrika ( from j=50 until j=77 compares to 13°N until 40°N,
//                                                                        from k=340 until k=352 compares to 20°W until 8°W )

	j_beg = 50;
	j_end = 78;
	k_beg = 340;
	k_end = 353;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}






// Atlantic

// Nordafrika		Westküsten
// Verbindung zwischen Guineastrom und Südatlantik

// Guinea-Strom im Westen from Afrika ( from j=77 until j=89 compares to 13°N until 1°N,
//                                                              from k=340 until k=352 compares to 20°W until 8°W )

	j_beg = 77;
	j_end = 90;
	k_beg = 340;
	k_end = 353;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}
}













void IC_Thermohalin::IC_Indischer_Ozean ( Array &h, Array &u, Array &v, Array &w )
{
// Ströme entlang der Küsten
// Schliessen der polaren, subpolaren und subtropischen atmosphärischen Wirbelsysteme


// Thermohalin Conveyor Belt

// Indischer Ozean

// diagonal from Indonesien/Australien until Kap Agulhas in Südafrika

// Westlich laufender Strom als Oberflächenströmung ( from j=104 until j=135 compares to 14°S until 45°S,
//                                                                                     from k=10 until k=120 compares to 10°O until 120°O )

	j_beg = 104;
	j_end = 136;
	k_beg = 10;
	k_end = 121;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 20;


	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run ); k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		k_run = k_run + 3;
		}
	j_run++;
	}



// westlich des konstanten Geschwindigkeitsstreifens diagonal from Indonesien/Australien until Kap Agulhas in Südafrika
// Westlich laufender Strom als Oberflächenströmung ( from j=104 until j=135 compares to 14°S until 45°S,
//                                                                                      from k=10 until k=120 compares to 10°O until 120°O )

	j_beg = 104;
	j_end = 136;
	k_beg = 10;
	k_end = 121;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 20;
	k_exp = 30;

	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run + k_step - k_exp ) <= ( k_end - k_exp ) )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run + k_step - k_exp ); k < ( k_beg + k_run + k_step ); k++ )
			{
				k_z = k - ( k_beg + k_run + k_step - k_exp );
				k_n = ( k_beg + k_run + k_step ) - ( k_beg + k_run + k_step - k_exp ) ;

				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_beg + k_run + k_step - 1 ] - v.x[ i ][ j ][ k_beg + k_run + k_step - k_exp ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg + k_run + k_step - k_exp ];
						w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_beg + k_run + k_step - 1 ] - w.x[ i ][ j ][ k_beg + k_run + k_step - k_exp ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg + k_run + k_step - k_exp ];
					}
				}
			}
			k_run = k_run + 3;
		}
	j_run++;
	}



// östlich des konstanten Geschwindigkeitsstreifens diagonal from Indonesien/Australien until Kap Agulhas in Südafrika
// Westlich laufender Strom als Oberflächenströmung ( from j=104 until j=135 compares to 14°S until 45°S,
//                                                                                      from k=10 until k=120 compares to 10°O until 120°O )

	j_beg = 104;
	j_end = 136;
	k_beg = 10;
	k_end = 121;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 20;
	k_exp = 30;


	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run + k_step + k_exp ) <= ( k_end + k_exp ) )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run + k_step ); k < ( k_beg + k_run + k_step + k_exp ); k++ )
			{
				k_z = k - ( k_beg + k_run + k_step );
				k_n = ( k_beg + k_run + k_step + k_exp ) - ( k_beg + k_run + k_step ) ;

				for ( int i = i_beg; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_beg + k_run + k_step + k_exp ] - v.x[ i ][ j ][ k_beg + k_run + k_step  - 1 ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg + k_run + k_step - 1 ];
						w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_beg + k_run + k_step + k_exp ] - w.x[ i ][ j ][ k_beg + k_run + k_step  - 1 ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg + k_run + k_step - 1 ];
					}
				}
			}
			k_run = k_run + 3;
		}
	j_run++;
	}




// Küstenströmungen

// Indischer Ozean

// Südafrika/Madagaskar		Ostküsten

// Agulhas-Strom im Süd-Osten from Südafrika ( from j=100 until j=126 compares to 10°S until 36°S,
//                                                                          from k=20 until k=42 compares to 20°O until 42°O )

	j_beg = 100;
	j_end = 127;
	k_beg = 20;
	k_end = 43;
	k_step = 14;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}


	}





// Küstenströmungen

// Indischer Ozean

// Südafrika/Madagaskar		Ostküsten

// Madagaskar-Strom im Osten ( from j=103 until j=125 compares to 13°S until 35°S,
//                                                  from k=45 until k=54 compares to 45°O until 54°O )

	j_beg = 103;
	j_end = 126;
	k_beg = 45;
	k_end = 55;
	k_step = 20;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}


	}





// Indischer Ozean

// Agulhas-Stroms an der Südspitze from Südafrika ( from j=125 until j=128 compares to 35°S until 38°S,
// parallel zum Äquator                                               from k=18 until k=27 compares to 18°O until 27°O )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 125; j < 129; j++ )
		{
			for ( int k = 18; k < 28; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						v.x[ i ][ j ][ k ] = + 0.0;
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		}
	}



// Indischer Ozean

// Agulhas-Stroms an der Südspitze from Südafrika ( from j=125 until j=128 compares to 35°S until 38°S,
// westlicher Übergang parallel zum Äquator             from k=18 until k=27 compares to 18°O until 27°O )

	j_beg = 125;
	j_end = 129;
	k_beg = 18;
	k_end = 28;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}


// Indischer Ozean

// Agulhas-Stroms an der Südspitze from Südafrika ( from j=125 until j=128 compares to 35°S until 38°S,
// östlicher Übergang parallel zum Äquator                from k=18 until k=27 compares to 18°O until 27°O )

	j_beg = 125;
	j_end = 129;
	k_beg = 18;
	k_end = 28;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}



// Indischer Ozean

// Agulhas-Stroms an der Südspitze from Südafrika ( from j=125 until j=128 compares to 35°S until 38°S,
// südlicher Übergang parallel zum Äquator               from k=18 until k=27 compares to 18°O until 27°O )

	j_beg = 125;
	j_end = 129;
	k_beg = 18;
	k_end = 28;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}







// Indischer Ozean

// Nordafrika		Ostküsten

// Somali-Strom im Nord-Osten from Nordafrika ( from j=80 until j=100 compares to 10°N until 10°S,
//                                                                           from k=38 until k=58 compares to 38°O until 58°O )

	j_beg = 80;
	j_end = 101;
	k_beg = 38;
	k_end = 59;
	k_step = 20;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;

		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}

	}




// Indischer Ozean

// Nordafrika		Ostküsten

// Verbindung zwischen Somali- und Agulhas-Strom

// Somali-Strom im Nord-Osten from Nordafrika ( from j=80 until j=100 compares to 10°N until 10°S,
//                                                                           from k=38 until k=58 compares to 38°O until 58°O )

	j_beg = 81;
	j_end = 101;
	k_beg = 38;
	k_end = 59;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}








// Indischer Ozean

// Indonesien		Westküste

// Sumatra/Java-Strom im Westen from Indonesien ( from j=85 until j=99 compares to 5°N until 9°S,
//                                                                                 from k=95 until k=105 compares to 95°O until 105°O )

	j_beg = 85;
	j_end = 100;
	k_beg = 95;
	k_end = 106;

	k_exp = 5;
	k_w = 4;

	v_grad = + 0.0008;
	k_step = 8;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;

		for ( int k = k_b - k_step - k_exp; k < ( k_b - k_w ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b - k_w ] - v.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + v.x[ i ][ j ][ k_b - k_step - k_exp ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b - k_w ] - w.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + w.x[ i ][ j ][ k_b - k_step - k_exp ];
			}
		}

	}







// Indischer Ozean

// Australien		Westküste

// Westaustralien-Strom im Westen from Australien ( from j=112 until j=125 compares to 22°S until 35°S,
//                                                                                 from k=103 until k=115 compares to 103°O until 115°O )

	j_beg = 112;
	j_end = 126;
	k_beg = 103;
	k_end = 116;

	k_exp = 5;
	k_w = 4;

	v_grad = 0.0008;
	k_step = 8;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + 0.1 * v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b - k_step - k_exp; k < ( k_b - k_w ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b - k_w ] - v.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + v.x[ i ][ j ][ k_b - k_step - k_exp ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b - k_w ] - w.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + w.x[ i ][ j ][ k_b - k_step - k_exp ];
			}
		}
	}
}











void IC_Thermohalin::IC_Pazifischer_Ozean ( Array &h, Array &u, Array &v, Array &w )
{
// Thermohalin Conveyor Belt und Ströme entlang der Küsten
// Schliessen der hydrosphärischen Wirbelsysteme


// Pazifischer Ozean
// Thermohalin Conveyor Belt

// Oberflächenströmung westlich from Kalifornien ( from j=52 until j=80 compares to 10°N until 38°N
//                                                                             from k=220 until k=250 compares to 110°W until 140°W )
/*
	j_beg = 52;
	j_end = 81;
	k_beg = 220;
	k_end = 251;

	v_grad = 0.0010;
	k_step = 8;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_b - k_a;
				if ( k_grad <= - 3 ) k_grad = - 2;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;
	}
*/




// Pazifischer Ozean
// Thermohalin Conveyor Belt

// Oberflächenströmung südlich der Aleuten ( from j=50 until j=56 compares to 34°N until 40°N
//                                           from k=169 until k=230 compares to 169°O until 130°W )
/*
	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 50; j < 57; j++ )
		{
			for ( int k = 169; k < 231; k++ )
			{
				v.x[ i ][ j ][ k ] = + 0.;
				w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
			}
		}
	}
*/





// Pazifik

// Strömung from Japan until Kalifornien in Nordamerika
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=50 until j=55 compares to 35°N until 40°N,
//                             	                                                             from k=160 until k=220 compares to 160°O until 140°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 50; j < 56; j++ )
		{
			for ( int k = 160; k < 221; k++ )
			{
				v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
			}
		}
	}





// Pazifik

// Strömung from Japan until Kalifornien in Nordamerika
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=50 until j=55 compares to 35°N until 40°N,
// westlicher Übergang parallel zum Äquator                       from k=160 until k=220 compares to 160°O until 140°W )

	j_beg = 50;
	j_end = 56;
	k_beg = 160;
	k_end = 221;

	for ( int i = 5; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}




// Pazifik

// Strömung from Japan until Kalifornien in Nordamerika
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=50 until j=55 compares to 35°N until 40°N,
// östlicher Übergang                                                            from k=160 until k=220 compares to 160°O until 140°W )

	j_beg = 50;
	j_end = 56;
	k_beg = 160;
	k_end = 221;

	for ( int i = 5; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}




// Pazifik

// Strömung from Japan until Kalifornien in Nordamerika
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=50 until j=55 compares to 35°N until 40°N,
// nördlicher Übergang                                                         from k=160 until k=220 compares to 160°O until 140°W )

	j_beg = 50;
	j_end = 56;
	k_beg = 160;
	k_end = 221;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}



// Strömung from Japan until Kalifornien in Nordamerika
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=50 until j=55 compares to 35°N until 40°N,
// südlicher Übergang                                                           from k=160 until k=220 compares to 160°O until 140°W )

	j_beg = 50;
	j_end = 56;
	k_beg = 160;
	k_end = 221;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}





// Pazifik

// Strömung from den Aleuten until zur ostwärtslaufenden Strömung
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=27 until j=50 compares to 40°N until 63°N,
// Nord-Süd, Oststrang                                                         from k=195 until k=200 compares to 165°W until 160°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 27; j < 51; j++ )
		{
			for ( int k = 195; k < 201; k++ )
			{
				v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
			}
		}
	}



// Pazifik

// Strömung from der ostwärtslaufenden Strömung until zu den Aleuten
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=27 until j=50 compares to 40°N until 63°N,
// Süd-Nord, Weststrang                                                      from k=195 until k=200 compares to 170°W until 165°W )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 27; j < 51; j++ )
		{
			for ( int k = 195; k < 201; k++ )
			{
				v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
			}
		}
	}




// Pazifik

// Strömung from der ostwärtslaufenden Strömung until zu den Aleuten
// Strömung from den Aleuten until zur ostwärtslaufenden Strömung
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=40 until j=50 compares to 40°N until 50°N,
// Übergang                                                                          from k=194 until k=198 compares to 166°W until 162°W )

	j_beg = 40;
	j_end = 51;
	k_beg = 194;
	k_end = 199;

	for ( int i = 5; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}







// Pazifik

// Strömung from Japan until Kalifornien in Nordamerika
// Nord-Süd
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=27 until j=50 compares to 40°N until 63°N,
// Nord-Süd, Oststrang östlicher Übergang                         from k=198 until k=204 compares to 162°W until 156°W )

	j_beg = 27;
	j_end = 51;
	k_beg = 198;
	k_end = 205;

	for ( int i = 5; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}




// Pazifik

// Strömung from Japan until Kalifornien in Nordamerika
// Süd-Nord
// T-Stück
// Nach Osten laufender Strom als Oberflächenströmung ( from j=27 until j=50 compares to 40°N until 63°N,
// Süd-Nord, Weststrang  westlicher Übergang                   from k=183 until k=192 compares to 177°W until 168°W )

	j_beg = 27;
	j_end = 51;
	k_beg = 183;
	k_end = 193;

	for ( int i = 5; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}










// Pazifischer Ozean
// Thermohalin Conveyor Belt

// südlich, parallel zum Äquator from Indonesien/Australien until Kalifornien in Nordamerika

// Westlich laufender Strom als Oberflächenströmung ( from j=102 until j=107 compares to 12°S until 17°S,
//                                                                                     from k=170 until k=260 compares to 170°O until 90°W )
/*
	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 102; j < 108; j++ )
		{
			for ( int k = 701; k < 261; k++ )
			{
				v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
			}
		}
	}
*/





// Küstenströmungen

// Pazifik

// Südamerika		Westküsten

// Peru- oder Humboldt-Strom im Westen from Südamerika ( from j=95 until j=130 compares to 5°S until 40°S,
// Verlängerung als Oberflächenstrom                                    from k=270 until k=290 compares to 70°W until 90°W )
/*
	j_beg = 95;
	j_end = 131;
	k_beg = 270;
	k_end = 291;

	v_grad = + 0.0010;
	k_step = 8;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_b - k_a;
				if ( k_grad <= - 3 ) k_grad = - 2;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;
	}
*/




// Pazifik

// Mittelamerika		Westküsten

// "Mittelamerika"-Strom im Westen from Mittelamerika ( Westindien ) ( from j=60 until j=85 compares to 5°N until 30°N,
//                                                                                                              from k=245 until k=275 compares to 85°W until 115°W )
/*
	j_beg = 60;
	j_end = 86;
	k_beg = 245;
	k_end = 276;

	v_grad = 0.0010;
	k_step = 6;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_b - k_a;
				if ( k_grad <= - 3 ) k_grad = - 2;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;
	}
*/





// Pazifik

// Nordamerika		Westküsten

// Strom am Golf from Alaska ( from j=29 until j=50 compares to 40°N until 61°N,
//                                              from k=215 until k=240 compares to 145°W until 120°W )

	j_beg = 29;
	j_end = 51;
	k_beg = 215;
	k_end = 241;

	k_exp = 10;
	k_w = 4;

	v_grad = 0.0008;
	k_step = 20;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b - k_step - k_exp; k < ( k_b - k_w ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b - k_w ] - v.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + v.x[ i ][ j ][ k_b - k_step - k_exp ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b - k_w ] - w.x[ i ][ j ][ k_b - k_step - k_exp ] ) * ( double ) ( k - ( k_b - k_step - k_exp ) ) / ( double ) ( ( k_b - k_w ) - ( k_b - k_step - k_exp ) ) + w.x[ i ][ j ][ k_b - k_step - k_exp ];
			}
		}
	}





// Pazifik

// Verbindung zwischen Golf from Alaska- und NEC-Strom

// Nord-West-Amerika		Westküsten

// Strom am Golf from Alaska ( from j=29 until j=50 compares to 40°N until 61°N,
//                                              from k=215 until k=240 compares to 145°W until 120°W )

	j_beg = 29;
	j_end = 51;
	k_beg = 215;
	k_end = 241;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}







// Pazifik

// Nord-Ost-Asien		Ostküsten

// Kamtschatka-Strom im Nordosten from Japan ( from j=25 until j=40 compares to 65°N until 50°N,
//                                                                            from k=159 until k=180 compares to 159°O until 180°O )

	j_beg = 25;
	j_end = 41;
	k_beg =159;
	k_end = 181;
	k_step = 15;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}

		for ( int k = k_b; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + 0.5 * v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}

		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}


		k_a = k_b;
		flip = 0;
	}







// Pazifik

// Nord-Ost-Asien		Ostküsten

// Oyashio-Strom im Nordosten from Japan ( from j=39 until j=48 compares to 51°N until 42°N,
//                                                                    from k=142 until k=161 compares to 142°O until 161°O )

	j_beg = 39;
	j_end = 49;
	k_beg =142;
	k_end = 162;
	k_step = 30;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}

		for ( int k = k_b; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}

		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}


		k_a = k_b;
		flip = 0;
	}






// Pazifik

// Nord-Ost-Asien		Ostküsten

// Kuroshio-Strom im Südosten from Japan ( from j=48 until j=68 compares to 22°N until 42°N,
//                                                                    from k=120 until k=142 compares to 120°O until 142°O )

	j_beg = 48;
	j_end = 69;
	k_beg = 120;
	k_end = 143;
	k_step = 30;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}
	}




// Pazifik

// Verbindung zwischen Kamtschatka- und Oyashio-Strom

// Nord-Ost-Asien		Ostküsten

// Kamtschatka- und Oyashio-Strom im Norden from Japan ( from j=45 until j=75 compares to 15°N until 45°N,
//                                                                                              from k=120 until k=160 compares to 120°O until 160°O )

	j_beg = 45;
	j_end = 76;
	k_beg = 120;
	k_end = 161;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}





// Pazifik

// Verbindung zwischen Oyashio- und Kuroshio-Strom

// Nord-Ost-Asien		Ostküsten

// Oyashio- und Kuroshio-Strom im Osten from Japan ( from j=45 until j=75 compares to 15°N until 45°N,
//                                                                                     from k=120 until k=160 compares to 120°O until 160°O )

	j_beg = 45;
	j_end = 76;
	k_beg = 120;
	k_end = 161;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}





// Pazifik

// Verbindung zwischen Kuroshio- und Ostchinesisches Meer
// Kuroshio- und Ostchinesisches Meer		Ostküsten

// Kuroshio-Strom und Ostchinesisches Meer im Osten from Japan ( from j=45 until j=75 compares to 15°N until 45°N,
//                                                                                                          from k=120 until k=160 compares to 120°O until 160°O )

	j_beg = 45;
	j_end = 76;
	k_beg = 120;
	k_end = 161;

	for ( int i = 5; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}







// Verbindung zwischen Kuroshio- und NEC-Strom

// Nord-Ost-Asien		Ostküsten

// Oyashio- und Kuroshio-Strom im Osten from Japan ( from j=45 until j=75 compares to 15°N until 45°N,
//                                                                                     from k=120 until k=160 compares to 120°O until 160°O )

	j_beg = 45;
	j_end = 76;
	k_beg = 120;
	k_end = 161;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}





// Pazifik

// Nord-Ost-Asien		Ostküsten

// Südchinesisches Meer im Osten from Vietnam ( from j=60 until j=80 compares to 23°N until 0°,
//                                             from k=92 until k=107 compares to 104°O until 120°O )
/*
	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 60; j < 81; j++ )
		{
			for ( int k = 92; k < 108; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}
*/



// Pazifik

// Nord-Ost-Asien		Ostküsten

// Tsushima-Strom im Westen from Japan until Indonesien ( from j=40 until j=53 compares to 45°N until 30°N,
//                                                     from k=116 until k= 129 compares to 130°O until 145°O )
/*
	l = 1;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 40; j < 54; j++ )
		{
			for ( int k = 116; k < 130; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					if ( l <= 3 )
					{
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						l++;
					}
				}
			if ( h.x[ i ][ j ][ k ] == 1. ) l = 1;
			}
		}
	}
*/




// Pazifik

// Philippinen (Indonesien) 		Ostküsten
 
// Strom im Osten der Philippinen ( from j=60 until j=71 compares to 10°N until 22°N,
//                                  from k=111 until k=116 compares to 125°O until 130°O )
/*
	j_beg = 60;
	j_end = 72;
	k_beg = 111;
	k_end = 117;
	k_step = 20;

	v_grad = 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;
	}
*/



// Pazifik

// Philippinen (Indonesien) Neuguinea		Ostküsten
 
// Durchfluss Zwischen Indonesien und Neuguinea ( from j=80 until j=97 compares to 10°N until 7°S
// Fortsetzung Richtung Pazifik                                   from k=125 until k=135 compares to 125°O until 135°O )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 80; j < 98; j++ )
		{
			for ( int k = 125; k < 136; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}






// Pazifik

// Philippinen (Indonesien) Neuguinea		Ostküsten
 
// Durchfluss Zwischen Indonesien und Neuguinea ( from j=80 until j=97 compares to 10°N until 7°S
// westlicher Übergang parallel zum Äquator             from k=125 until k=135 compares to 125°O until 135°O )

	j_beg = 80;
	j_end = 98;
	k_beg = 125;
	k_end = 136;

	for ( int i = 5; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}




// Pazifik

// Philippinen (Indonesien) Neuguinea		Ostküsten
 
// Durchfluss Zwischen Indonesien und Neuguinea ( from j=80 until j=97 compares to 10°N until 7°S
// östlicher Übergang parallel zum Äquator                from k=125 until k=135 compares to 125°O until 135°O )

	j_beg = 80;
	j_end = 98;
	k_beg = 125;
	k_end = 136;

	for ( int i = 5; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}





// Pazifik

// Philippinen (Indonesien) Neuguinea		Ostküsten
 
// Durchfluss Zwischen Indonesien und Neuguinea ( from j=80 until j=97 compares to 10°N until 7°S
// nördlicher Übergang parallel zum Äquator              from k=120 until k=135 compares to 125°O until 135°O )

	j_beg = 80;
	j_end = 98;
	k_beg = 120;
	k_end = 136;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}


// Pazifik

// Philippinen (Indonesien) Neuguinea		Ostküsten
 
// Durchfluss Zwischen Indonesien und Neuguinea ( from j=80 until j=97 compares to 10°N until 7°S
// südlicher Übergang parallel zum Äquator               from k=125 until k=135 compares to 125°O until 135°O )

	j_beg = 80;
	j_end = 98;
	k_beg = 125;
	k_end = 136;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}





// Pazifik

// Papua-Neuguinea (Indonesien)		Nordküsten from Neuguinea
 
// Strom im Nord-Osten from Neuguinea ( from j=90 until j=100 compares to 0°S until 10°S,
//                                                                from k=135 until k=150 compares to 135°O until 150°O )

	j_beg = 90;
	j_end = 100;
	k_beg = 135;
	k_end = 151;
	k_step = 20;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = - v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}
	}






// Pazifik

// Australien		Ostküsten im Korallen- und Tasmanensee
 
// Strom im Osten Australien ( from j=101 until j=140 compares to 11°S until 50°S,
//                                               from k=142 until k=158 compares to 142°O until 158°O )

	j_beg = 101;
	j_end = 141;
	k_beg = 142;
	k_end = 159;
	k_step = 20;

	v_grad = + 0.001;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_end; k > k_beg; k-- )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_end; k > k_beg; k-- )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_a - k_b;
				if ( k_grad >= 2 ) k_grad = 1;
				if ( k_grad <= 0 ) k_grad = 1;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k < ( k_b + k_step ); k++ )
		{
			for ( int i = i_beg; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - 0.4 * v_grad * ( double ) ( k_grad ) * IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
		k_a = k_b;
		flip = 0;


		for ( int k = k_b + 4; k < ( k_b + k_step ); k++ ) 
		{
			for ( int i = i_beg; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_b + k_step +1 ] - v.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + v.x[ i ][ j ][ k_b + 4 ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_b + k_step +1 ] - w.x[ i ][ j ][ k_b + 4 ] ) * ( double ) ( k - ( k_b + 4 ) ) / ( double ) ( ( k_b + k_step -1 ) - ( k_b + 4 ) ) + w.x[ i ][ j ][ k_b + 4 ];
			}
		}
	}




// Pazifik

// Verbindung im Norden from Australien und Neuseeland

// Australien		Ostküsten im Korallen- und Tasmanensee
 
// Strom im Osten from Australien ( from j=110 until j=128 compares to 20°S until 38°S,
//                                                     from k=142 until k=158 compares to 142°O until 158°O )

	j_beg = 110;
	j_end = 129;
	k_beg = 142;
	k_end = 159;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}






// Pazifik

// Australien		Südküsten im Südaustralischen Becken
 
// Strom im Süden from Australien ( from j=123 until j=132 compares to 33°S until 42°S,
//                                                      from k=116 until k=148 compares to 116°O until 148°O )

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 123; j < 133; j++ )
		{
			for ( int k = 116; k < 149; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}





// Verbindung zwischen Australien und Zirkumpolarstrom

// Australien		Südküste

	j_beg = 128;
	j_end = 133;
	k_beg = 116;
	k_end = 149;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}
}






void IC_Thermohalin::IC_Nord_Polar_Meer ( Array &h, Array &u, Array &v, Array &w )
{
// Ströme entlang der Küsten
// Schliessen der polaren, subpolaren und subtropischen atmosphärischen Wirbelsysteme


// Arktische Strömungen

// Verlängerung zum Nordpol und Verbreiterung ( depth -45 m )
// Durchströmung der Bering-Straße ( from j=0 until j=30 compares to 60°N until 90°N,
//                                                           from k=180 until k=194 compares to 180°W until 166°W )


	for ( int i = im-4; i < im; i++ )
	{
		for ( int j = 0; j < 31; j++ )
		{
			for ( int k = 180; k < 195; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}






// Arktische Strömungen

// Verbindung zwischen Beringmeer und Nordpazifik

// Verlängerung until zum Nordpol und Verbreiterung ( depth -45 m )
// Durchströmung der Bering-Straße ( from j=0 until j=30 compares to 60°N until 90°N,
//                                                           from k=180 until k=194 compares to 180°W until 166°W )

	j_beg = 31;
	j_end = 36;
	k_beg = 180;
	k_end = 195;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}






// Arktische Strömungen

// westliche Verbindung zwischen Beringmeer und Nordpazifik

// Verlängerung until zum Nordpol und Verbreiterung ( depth -45 m )
// Durchströmung der Bering-Straße ( from j=0 until j=30 compares to 60°N until 90°N,
//                                                           from k=180 until k=194 compares to 180°W until 166°W )

	j_beg = 0;
	j_end = 36;
	k_beg = 180;
	k_end = 195;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
//					u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_end ] - u.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + u.x[ i ][ j ][ k_beg ];
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}






// Arktische Strömungen

// östliche Verbindung zwischen Beringmeer und Nordpazifik

// Verlängerung until zum Nordpol und Verbreiterung ( depth -45 m )
// Durchströmung der Bering-Straße ( from j=0 until j=30 compares to 60°N until 90°N,
//                                                           from k=180 until k=194 compares to 180°W until 166°W )

	j_beg = 0;
	j_end = 31;
	k_beg = 180;
	k_end = 195;

	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
//					u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_end ] - u.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + u.x[ i ][ j ][ k_beg ];
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}






// Fortsetzung der Beringstraßenstrom über den Nordpol hinaus ( depth -45 m )

// Verlängerung der Bering-Straße until zum Nordpol und Verbreiterung  ( from j=0 until j=20 compares to 90°N until 70°N,
//                                                                                                                  from k=320 until k=km compares to 20°W until 0° )
/*
	for ( int i = im-4; i < im; i++ )
	{
		for ( int j = 0; j < 21; j++ )
		{
			for ( int k = 320; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}
*/
}







void IC_Thermohalin::IC_South_Polar_Sea ( Array &h, Array &u, Array &v, Array &w, Array &c )
{
// flow along coasts
// closing the polar, subpolar and subtropic atmospheric circulation systems


// south polar sea

// setting the South Pole with zero velocity

// antarctic circumpolar current ( -5000m deep ) ( from j=147 until j=152 compares to 57°S until 62°S,
//                                                                            from k=0 until k=km compares to 0° until 360° )


	for ( int i = 10; i < im-3; i++ )
	{
		for ( int j = 147; j < 153; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
//					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - 5 ) / ( double ) ( im-2 - 5 );
//					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - 5 ) / ( double ) ( im-2 - 5 );
//					v.x[ i ][ j ][ k ] = - IC_water;
					w.x[ i ][ j ][ k ] = IC_water;
				}
			}
		}
	}



// antarctic circumpolar current ( -5000m deep )

// nördlich vom Antarktischen Zirkumpolarstrom
// antarctic circumpolar current ( from j=147 until j=152 compares to 57°S until 62°S,
//                                                        from k=0 until k=km compares to 0° until 360° )

	j_beg = 144;
	j_end = 150;
	k_beg = 0;
	k_end = km;

	for ( int i = 10; i < im - 3; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = ( c.x[ i ][ j_end ][ k ] - c.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + c.x[ i ][ j_beg ][ k ];
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}




// antarctic circumpolar current ( -5000m deep )

// south of antarctic circumpolar current
// antarctic circumpolar current ( from j=147 until j=152 compares to 57°S until 62°S,
//                                                 from k=0 until k=km compares to 0° until 360° )

	j_beg = 150;
	j_end = 155;
	k_beg = 0;
	k_end = km;

	for ( int i = 10; i < im - 3; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = ( c.x[ i ][ j_end ][ k ] - c.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + c.x[ i ][ j_beg ][ k ];
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}






// south polar sea

// subantarctic  mode water ( -1000m deep )

// subantarctcic front parallel to the circumpolar current ( from j=144 until j=146 compares to 54°S until 56°S,
//                                                                                        from k=0 until k=km compares to 0° until 360° )

	for ( int i = i_beg+10; i < im; i++ )
	{
		for ( int j = 144; j < 147; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
//					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg + 5 ) / ( double ) ( i_max - i_beg + 5 );
					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_beg + 5 ) / ( double ) ( i_max - i_beg + 5 );
				}
			}
		}
	}



// north of the constant velocity belt of the subantarctic mode water ( -1000m deep )
// subantarctcic front parallel to the circumpolar current ( from j=144 until j=146 compares to 54°S until 56°S,
//                                                                                        from k=0 until k=km compares to 0° until 360° )

	j_beg = 144;
	j_end = 147;
	k_beg = 0;
	k_end = km;

	for ( int i = i_beg+10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}




// north of the constant velocity belt of the subantarctic mode water ( -1000m deep )
// subantarctcic front parallel to the circumpolar current ( from j=144 until j=146 compares to 54°S until 56°S,
//                                                                                        from k=0 until k=km compares to 0° until 360° )

	j_beg = 144;
	j_end = 147;
	k_beg = 0;
	k_end = km;

	for ( int i = i_beg+10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}




// south polar sea

// Weddell sea

// diagonal from Lasarew sea until antarktic peninsula

// south-west directed flow at the surface ( from j=155 until j=165 compares to 65°S until 75°S,
//                                                                  from k=299 until k=km compares to 61°W until 0° )

	j_beg = 155;
	j_end = 166;
	k_beg = 299;
	k_end = km;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 50;


	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run ) <=  k_end )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run ); k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = 10; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
//						v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		k_run = k_run+1;
		}
	j_run++;
	}



// Weddell sea
// north of Lasarew sea until antarktic peninsula
// south-west directed flow at the surface ( from j=155 until j=165 compares to 65°S until 75°S,
//                                                                  from k=299 until k=km compares to 61°W until 0° )

	j_beg = 155;
	j_end = 166;
	k_beg = 299;
	k_end = km;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}







// Weddell sea
// east of of Lasarew sea until antarktic peninsula
// south-west directed flow at the surface ( from j=155 until j=165 compares to 65°S until 75°S,
//                                                                  from k=352 until k=km compares to 8°W until 0° )

	j_beg = 155;
	j_end = 166;
	k_beg = km - 8;
	k_end = km;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}


	j_beg = 155;
	j_end = 166;
	k_beg = km - 8;
	k_end = km;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}


	j_beg = 155;
	j_end = 166;
	k_beg = km - 43;
	k_end = km - 1;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}





// south polar sea

// Ross sea

// diagonally from Martin peninsula sea until Kap Adare

// south-west directed flow at the surface ( from j=155 until j=168 compares to 65°S until 78°S,
//                                                                  from k=180 until k=280 compares to 180°W until 80°W )

	j_beg = 155;
	j_end = 169;
	k_beg = 180;
	k_end = 281;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 90;


	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run ); k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = 10; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
//						v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		k_run = k_run+1;
		}
	j_run++;
	}





// Ross sea

// north of Martin peninsula sea until Kap Adare

// south-west directed flow at the surface ( from j=155 until j=168 compares to 65°S until 78°S,
//                                                                  from k=180 until k=245 compares to 180°W until 115°W )

	j_beg = 155;
	j_end = 169;
	k_beg = 180;
	k_end = 246;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}





// Ross sea

// east of Martin peninsula sea until Kap Adare

// south-west directed flow at the surface ( from j=155 until j=168 compares to 65°S until 78°S,
//                                                                  from k=180 until k=245 compares to 180°W until 115°W )


	j_beg = 155;
	j_end = 169;
	k_beg = 254;
	k_end = 262;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}


	j_beg = 155;
	j_end = 169;
	k_beg = 254;
	k_end = 262;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
				}
			}
		}
	}


	j_beg = 155;
	j_end = 169;
	k_beg = 190;
	k_end = 250;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}




	j_beg = 155;
	j_end = 169;
	k_beg = 246;
	k_end = 255;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}






// south polar sea

// antarctic Indic basin

// diagonally from Adelieland until Queen Mary Land

// south-west directed flow at the surface ( from j=152 until j=158 compares to 62°S until 68°S,
//                                                                  from k=0 until k=160 compares to 0°O until 160°O )

	j_beg = 152;
	j_end = 159;
	k_beg = 0;
	k_end = 161;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 110;


	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run ); k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = 10; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
//						v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg );
					}
				}
			}
		k_run = k_run+1;
		}
	j_run++;
	}


// antarctic Indic basin

// north of Adelieland until Queen Mary Land

// south-west directed flow at the surface ( from j=152 until j=158 compares to 62°S until 68°S,
//                                                                  from k=0 until k=160 compares to 0°O until 160°O )

	j_beg = 152;
	j_end = 159;
	k_beg = 0;
	k_end = 161;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end + 1; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				j_z = j - j_beg;
				j_n = j_end - j_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_end ][ k ] - v.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + v.x[ i ][ j_beg ][ k ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_end ][ k ] - w.x[ i ][ j_beg ][ k ] ) * ( double ) ( j_z ) / ( double ) ( j_n ) + w.x[ i ][ j_beg ][ k ];
				}
			}
		}
	}



	j_beg = 152;
	j_end = 159;
	k_beg = 132;
	k_end = 138;

	for ( int i = 10; i < im; i++ )
	{
		for ( int j = j_beg; j < j_end; j++ )
		{
			for ( int k = k_beg; k < k_end; k++ )
			{
				k_z = k - k_beg;
				k_n = k_end - k_beg;

				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_end ] - v.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + v.x[ i ][ j ][ k_beg ];
					w.x[ i ][ j ][ k ] = ( w.x[ i ][ j ][ k_end ] - w.x[ i ][ j ][ k_beg ] ) * ( double ) ( k_z ) / ( double ) ( k_n ) + w.x[ i ][ j ][ k_beg ];
				}
			}
		}
	}
}






void IC_Thermohalin::IC_EquatorialCurrents ( Array &h, Array &u, Array &v, Array &w )
{
// currents along the equator
// equatorial undercurrent - Cromwell flow, EUC
// upwelling at the end of the equatorial undercurrent - Cromwell flow, EUC
// equatorial intermediate current, EIC
// nothern and southern equatorial subsurface counter-currents, NSCC und SSCC
// nothern and southern equatorial counter-currents, NECC und SECC
// domes at the end of nothern and southern equatorial counter-currents, NECC und SECC

// for the depth of the sea compares i=1 to a depth of 150 m	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	i_EIC_u = im - 14;
	i_EIC_o = im - 6;
	i_SCC_u = im - 10;
	i_SCC_o = im - 6;
	i_ECC_u = im - 6;

	i_ECC_o = im;



// equatorial currents and counter-currents



//    §§§§§§§§§§§§§§§§§§§§§§§§§		Pacific ocean		§§§§§§§§§§§§§§§§§§§§§§§§§§



// Pacific ocean

// equatorial undercurrent - Cromwell current ( EUC, i=im-2 until i=im-1 compares to -100 until -200m depth )
// equatorial undercurrent - Cromwell current ( from j=87 until j=93 compares to 3°N until 3°S,
//                                                                         from k=145 until k=270 compares to 145°O until 90°W )


	for ( int i = im-4; i < im-2; i++ )
	{
		for ( int j = 87; j < 94; j++ )
		{
			for ( int k = 145; k < 271; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = + IC_water;
				}
			}
		}
	}




// Pacific ocean

// downward flow at the end of the equatorial undercurrent - Cromwell current ( EUC, -200m, 3°N - 3°S, 85°W - 90°W )

	for ( int i = im-14; i < im-2; i++ )
	{
		for ( int j = 87; j < 94; j++ )
		{
			for ( int k = 265; k < 271; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = + IC_water;
						w.x[ i ][ j ][ k ] = - IC_water;
					}
				}
			}
		}
	}


// Pacific ocean

// upward flow at the end of the equatorial undercurrent - Cromwell current ( EUC, -200m, 3°N - 3°S, 140°O - 145°O )

	for ( int i = im-14; i < im-2; i++ )
	{
		for ( int j = 87; j < 94; j++ )
		{
			for ( int k = 140; k < 146; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = - IC_water;
						w.x[ i ][ j ][ k ] = + IC_water;
					}
				}
			}
		}
	}



// Pacific ocean

// equatorial intermediate current ( EIC, i=im-4 until i=im-2 compares to -300 until -1000m depth )
// equatorial intermediate current ( from j=88 until j=92 compares to 2°N until 2°S,
//                                                      from k=145 until k=270 compares to 145°O until 90°W )

	for ( int i = i_EIC_u; i < i_EIC_o; i++ )
	{
		for ( int j = 88; j < 93; j++ )
		{
			for ( int k = 145; k < 271; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_EIC_u ) / ( double ) ( i_EIC_o - i_EIC_u );
				}
			}
		}
	}




// Pacific ocean

// equatorial northern and southern subsurface counter-current

// equatorial northern subsurface counter-current ( NSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial northern subsurface counter-current ( from j=86 until j=87 compares to 3°N until 4°N,
//                                                                               from k=135 until k=270 compares to 135°O until 90°W )

	for ( int i = i_SCC_u; i < i_SCC_o; i++ )
	{
		for ( int j = 86; j < 88; j++ )
		{
			for ( int k = 135; k < 271; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_SCC_o - i_SCC_u );
				}
			}
		}
	}


// Pacific ocean

// equatorial southern subsurface counter-current ( SSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial southern subsurface counter-current ( from j=93 until j=94 compares to 3°S until 4°S,
//                                                                               from k=165 until k=270 compares to 165°O until 90°W )

	for ( int i = i_SCC_u; i < i_SCC_o; i++ )
	{
		for ( int j = 93; j < 95; j++ )
		{
			for ( int k = 165; k < 271; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_SCC_o - i_SCC_u );
				}
			}
		}
	}



// Pacific ocean

// equatorial current at the surface


// equatorial northern counter-current ( NECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial northern counter-current ( from j=83 until j=87 compares to 3°N until 7°N,
//                                                            from k=135 until k=270 compares to 135°O until 90°W )


	for ( int i = i_ECC_u; i < i_ECC_o; i++ )
	{
		for ( int j =83; j < 88; j++ )
		{
			for ( int k = 135; k < 271; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_ECC_u ) / ( double ) ( i_ECC_o - i_ECC_u );
				}
			}
		}
	}




// equatorial current at the surface

// equatorial southern counter-current 

// Pacific ocean

// equatorial northern and southern counter-currents at the surface

// equatorial currents between NECC and SECC 

// equatorial northern counter-current ( from j=87 until j=93 compares to 3°N until 3°S,
// equatorial northern counter-current   from k=135 until k=270 compares to 135°O until 90°W )

	for ( int i = im-2; i < im; i++ )
	{
		for ( int j = 87; j < 94; j++ )
		{
			for ( int k = 135; k < 271; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  - IC_water;
				}
			}
		}
	}






// Pacific ocean

// equatorial northern and southern counter-currents at the surface

// equatorial currents between NECC and SECC 

// equatorial southern counter-current ( SECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial southern counter-current ( from j=93 until j=96 compares to 3°S until 6°S,
//                                                             from k=155 until k=270 compares to 155°O until 90°W )


	for ( int i = i_ECC_u; i < i_ECC_o; i++ )
	{
		for ( int j = 93; j < 97; j++ )
		{
			for ( int k = 155; k < 271; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_ECC_u ) / ( double ) ( i_ECC_o - i_ECC_u );
				}
			}
		}
	}




// Pacific ocean

// Guatemala dome at the end of the northern equatorial counter-current at the surface, Guinea current = NECC

// Guatemala-dome ( zum NECC gehörend, i=im-1 until i=im compares to 0 until -200m depth )
// Guatemala-dome ( from j=83 until j=87 compares to 3°N until 7°N,
//                                from k=265 until k=275 compares to 85°W until 95°W )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 83; j < 88; j++ )
		{
			for ( int k = 265; k < 276; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}


// Pacific ocean

// Equador dome at the end of the southern equatorial counter-current at the surface

// Equador-dome ( zum SECC gehörend, i=im-1 until i=im compares to 0 until -200m depth )
// Equador-dome ( from j=93 until j=96 compares to 3°S until 6°S,
//                           from k=265 until k=275 compares to 85°W until 95°W )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 93; j < 97; j++ )
		{
			for ( int k = 265; k < 276; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] =  - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}




// Pacific ocean

// downward flow opposit the Guatemala dome at the end of the northern equatorial counter-current at the surface, Guinea current = NECC

// Guatemala-dome ( from j=83 until j=87 compares to 3°N until 7°N,
//                                from k=135 until k=140 compares to 135°W until 140°W )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 83; j < 88; j++ )
		{
			for ( int k = 135; k < 141; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}


// Pacific ocean

// downward flow opposit the Equador-dome dome at the end of the northern equatorial counter-current at the surface, Guinea current = NECC

// Equador-dome ( from j=93 until j=96 compares to 3°S until 6°S,
//                            from k=155 until k=160 compares to 155°O until 160°O )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 93; j < 97; j++ )
		{
			for ( int k = 155; k < 161; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] =  - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] =  - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] =  - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}





//    §§§§§§§§§§§§§§§§§§§§§§§§§		Indic ocean		§§§§§§§§§§§§§§§§§§§§§§§§§§


// equator currents and counter-currents


// Indic ocean

// equatorial under current - Cromwell current ( EUC, i=im-2 until i=im-1 compares to -100 until -200m depth )
// equatorial under current - Cromwell current ( from j=89 until j=91 compares to 1°N until 1°S,
//                                                                             from k=55 until k=90 compares to 55°O until 90°O )

	for ( int i = im-2; i < im-1; i++ )
	{
		for ( int j = 89; j < 92; j++ )
		{
			for ( int k = 55; k < 91; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = IC_water;
				}
			}
		}
	}



// Indic ocean

// downward flow at the end of the equatorial under-current - Cromwell-current ( EUC, -200m, 3°N - 3°S, 85°W - 90°W )
// equatorial under current - Cromwell current ( from j=89 until j=91 compares to 1°N until 1°S,
//                                                                              from k=85 until k=90 compares to 85°O until 90°O )

	for ( int i = im-14; i < im-2; i++ )
	{
		for ( int j = 89; j < 92; j++ )
		{
			for ( int k = 84; k < 91; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = + IC_water;
						w.x[ i ][ j ][ k ] = - IC_water;
					}
				}
			}
		}
	}



// Indic ocean

// upward flow at the end of the equatorial under-current - Cromwell-current ( EUC, -200m, 3°N - 3°S, 140°W - 145°W )
// equatorial under current - Cromwell current ( from j=89 until j=91 compares to 1°N until 1°S,
//                                                                             from k=50 until k=55 compares to 50°O until 55°O )

	for ( int i = im-1; i < im-2; i++ )
	{
		for ( int j = 89; j < 92; j++ )
		{
			for ( int k = 50; k < 56; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = - IC_water;
						w.x[ i ][ j ][ k ] = + IC_water;
					}
				}
			}
		}
	}





// Indic ocean

// equatorial intermediat current ( EIC, i=im-4 until i=im-2 compares to -300 until -1000m depth )
// equatorial intermediat current ( from j=89 until j=91 compares to 1°N until 1°S,
//                                                      from k=55 until k=90 compares to 55°O until 90°O )

	for ( int i = i_EIC_u; i < i_EIC_o; i++ )
	{
		for ( int j = 89; j < 92; j++ )
		{
			for ( int k = 55; k < 91; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_EIC_u ) / ( double ) ( i_EIC_o - i_EIC_u );
				}
			}
		}
	}







// Indic ocean

// northern and southern equatorial subsurface counter-currents
// equatorial northern subsurface counter-current ( NSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial northern subsurface counter-current ( from j=85 until j=88 compares to 2°N until 5°N,
//                                                                         from k=55 until k=90 compares to 55°O until 90°O )

	for ( int i = i_SCC_u; i < i_SCC_o; i++ )
	{
		for ( int j = 85; j < 89; j++ )
		{
			for ( int k = 55; k < 91; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_SCC_o - i_SCC_u );
				}
			}
		}
	}



// Indic ocean

// equatorial northern subsurface counter-current ( SSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial northern subsurface counter-current ( from j=92 until j=95 compares to 2°S until 5°S,
//                                                                        from k=55 until k=90 compares to 55°O until 90°O )

	for ( int i = i_SCC_u; i < i_SCC_o; i++ )
	{
		for ( int j = 92; j < 96; j++ )
		{
			for ( int k = 55; k < 91; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_SCC_o - i_SCC_u );
				}
			}
		}
	}




// Indic ocean

// equatorial northern counter-current at the surface, SW Monsun-current = NECC

// equatorial northern counter-current ( NECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial northern counter-current ( from j=85 until j=88 compares to 2°N until 5°N,
//                                                                 from k=55 until k=90 compares to 55°O until 90°O )

	for ( int i = i_ECC_u; i < i_ECC_o; i++ )
	{
		for ( int j = 85; j < 99; j++ )
		{
			for ( int k = 55; k < 91; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_ECC_u ) / ( double ) ( i_ECC_o - i_ECC_u );
				}
			}
		}
	}






// Indic ocean

// Nördliche und südliche äquatoriale Gegenströmung an der Oberfläche

// Äquatoriale Strömung zwischen NECC und SECC 

// Äquatoriale Gegenströmungen ( from j=89 until j=91 compares to 1°N until 1°S,
//                                                     from k=55 until k=90 compares to 55°O until 90°O )

	for ( int i = im-2; i < im; i++ )
	{
		for ( int j = 89; j < 92; j++ )
		{
			for ( int k = 55; k < 91; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  - IC_water;
				}
			}
		}
	}







// Indic ocean

// Südliche äquatoriale Gegenströmung an der Oberfläche

// equatorial southern counter-current ( SECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial southern counter-current ( from j=92 until j=95 compares to 2°S until 5°S,
//                                                                from k=55 until k=90 compares to 55°O until 90°O )

	for ( int i = i_ECC_u; i < i_ECC_o; i++ )
	{
		for ( int j = 92; j < 95; j++ )
		{
			for ( int k = 55; k < 91; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_ECC_u ) / ( double ) ( i_ECC_o - i_ECC_u );
				}
			}
		}
	}



// Indic ocean

// Nordsumatra-Dom am Ende der Nördlichen äquatorialen Gegenströmung an der Oberfläche, Guinea-Strom = NECC

// Nordsumatra-Dom ( zum NECC gehörend, i=im-1 until i=im compares to 0 until -200m depth )
// Nordsumatra-Dom ( from j=85 until j=88 compares to 2°N until 5°N,
//                                  from k=88 until k=95 compares to 88°O until 95°O )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 85; j < 89; j++ )
		{
			for ( int k = 88; k < 96; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}


// Indic ocean

// Südsumatra-Dom am Ende der Südlichen äquatorialen Gegenströmung an der Oberfläche



// Südsumatra-Dom ( zum SECC gehörend, i=im-1 until i=im compares to 0 until -200m depth )
// Südsumatra-Dom ( from j=92 until j=95 compares to 2°S until 5°S,
//                                 from k=90 until k=98 compares to 90°O until 98°O )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 92; j < 96; j++ )
		{
			for ( int k = 90; k < 99; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}





// Indic ocean

// Abströmung gegenüber dem Nordsumatra-Dom am Ende der Nördlichen äquatorialen Gegenströmung an der Oberfläche, Guinea-Strom = NECC

// Abströmung gegenüber dem Nordsumatra-Dom ( zum NECC gehörend, i=im-1 until i=im compares to 0 until -200m depth )
// Abströmung gegenüber dem Nordsumatra-Dom ( from j=84 until j=87 compares to 3°N until 6°N,
//                                                                                 from k=55 until k=59 compares to 55°O until 59°O )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 84; j < 88; j++ )
		{
			for ( int k = 55; k < 60; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}


// Indic ocean

// Abströmung gegenüber dem Südsumatra-Dom am Ende der Südlichen äquatorialen Gegenströmung an der Oberfläche



// Abströmung gegenüber dem Südsumatra-Dom ( zum SECC gehörend, i=im-1 until i=im compares to 0 until -200m depth )
// Abströmung gegenüber dem Südsumatra-Dom ( from j=93 until j=96 compares to 3°S until 6°S,
//                                                                                from k=55 until k=59 compares to 55°O until 59°O )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 93; j < 97; j++ )
		{
			for ( int k = 55; k < 60; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] =  - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] =  - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}





//    §§§§§§§§§§§§§§§§§§§§§§§§§		Atlantic ocean		§§§§§§§§§§§§§§§§§§§§§§§§§§







// Atlantic ocean

// equator currents and counter-currents

// equatorial under current - Cromwell current ( EUC, i=im-2 until i=im-1 compares to -100 until -200m depth )
// equatorial under current - Cromwell current ( from j=89 until j=91 compares to 1°N until 1°S,
//                                                                             from k=326 until k=355 compares to 34°W until 5°W )

	for ( int i = im-4; i < im-2; i++ )
	{
		for ( int j = 89; j < 92; j++ )
		{
			for ( int k = 326; k < 356; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = IC_water;
				}
			}
		}
	}






// Atlantic ocean

// upward flow at the end of the equatorial under-current - Cromwell-current ( EUC, -200m, 3°N - 3°S, 30°W - 35°W )
// equatorial under current - Cromwell current ( from j=89 until j=91 compares to 1°N until 1°S,
//                                                                             from k=326 until k=330 compares to 34°W until 30°W )

	for ( int i = im-14; i < im-2; i++ )
	{
		for ( int j = 89; j < 92; j++ )
		{
			for ( int k = 326; k < 330; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = - IC_water;
						w.x[ i ][ j ][ k ] = + IC_water;
					}
				}
			}
		}
	}




// Atlantic ocean

// downward flow at the end of the equatorial under-current - Cromwell-current ( EUC, -200m, 2°N - 0°, 3°W - 7°W )
// equatorial under current - Cromwell current ( from j=89 until j=91 compares to 1°N until 1°S,
//                                                                             from k=350 until k=355 compares to 10°W until 5°W )

	for ( int i = im-14; i < im-2; i++ )
	{
		for ( int j = 89; j < 92; j++ )
		{
			for ( int k = 350; k < 356; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = + IC_water;
						w.x[ i ][ j ][ k ] = - IC_water;
					}
				}
			}
		}
	}




// Atlantic ocean

// equatorial intermediat current ( EIC, i=im-4 until i=im-2 compares to -300 until -1000m depth,  )
// equatorial intermediat current ( from j=88 until j=92 compares to 2°N until 2°S,
//                                                      from k=330 until k=355 compares to 30°W until 5°W )

	for ( int i = i_EIC_u; i < i_EIC_o; i++ )
	{
		for ( int j = 88; j < 93; j++ )
		{
			for ( int k = 330; k < 356; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_EIC_u ) / ( double ) ( i_EIC_o - i_EIC_u );
				}
			}
		}
	}




// Atlantic ocean

// northern and southern equatorial subsurface counter-currents
// equatorial northern subsurface counter-current ( NSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial northern subsurface counter-current ( from j=86 until j=87 compares to 3°N until 4°N,
//                                                                         from k=320 until k=340 compares to 40°W until 20°W )

	for ( int i = i_SCC_u; i < i_SCC_o; i++ )
	{
		for ( int j = 86; j < 88; j++ )
		{
			for ( int k = 320; k < 341; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_SCC_o - i_SCC_u );
				}
			}
		}
	}




// Atlantic ocean

// equatorial northern subsurface counter-current ( SSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial northern subsurface counter-current ( from j=93 until j=94 compares to 3°S until 4°S,
//                                                                        from k=329 until k=355 compares to 31°W until 5°W )

	for ( int i = i_SCC_u; i < i_SCC_o; i++ )
	{
		for ( int j = 93; j < 95; j++ )
		{
			for ( int k = 329; k < 355; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_SCC_o - i_SCC_u );
				}
			}
		}
	}






// Atlantic ocean

// equatorial currents at the surface

// equatorial northern counter-current ( NECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial northern counter-current ( from j=84 until j=87 compares to 3°N until 6°N,
//                                                                 from k=330 until k=340 compares to 30°W until 20°W )

	for ( int i = i_ECC_u; i < i_ECC_o; i++ )
	{
		for ( int j = 84; j < 88; j++ )
		{
			for ( int k = 330; k < 340; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_ECC_u ) / ( double ) ( i_ECC_o - i_ECC_u );
				}
			}
		}
	}






// Atlantic ocean

// northern and southern equatorial currents at the surface

// Äquatoriale Strömung zwischen NECC und SECC 

// equatorial northern counter-current ( from j=87 until j=93 compares to 3°N until 3°S,
//                                                                 from k=330 until k=355 compares to 30°W until 5°W )

	for ( int i = im-2; i < im; i++ )
	{
		for ( int j = 87; j < 94; j++ )
		{
			for ( int k = 330; k < 356; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  - IC_water;
				}
			}
		}
	}






// Atlantic ocean

// northern and southern equatorial currents at the surface

// equatorial southern counter-current ( SECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial southern counter-current ( from j=83 until j=85 compares to 3°S until 6°S,
//                                                               from k=293 until k=316 compares to 30°W until 5°W )
/*

	for ( int i = i_ECC_u; i < i_ECC_o; i++ )
	{
		for ( int j = 83; j < 86; j++ )
		{
			for ( int k = 293; k < 317; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_ECC_u ) / ( double ) ( i_ECC_o - i_ECC_u );
				}
			}
		}
	}
*/








// Atlantic ocean

// Guinea-dome at the end of the northern equatorial counter-current at the surface, Guinea-current = NECC

// Guinea-dome ( zum NECC gehörend, i=im-1 until i=im compares to 0 until -200m depth )
// Guinea-dome ( from j=84 until j=87 compares to 3°N until 6°N,
//                         from k=340 until k=345 compares to 15°W until 20°W )
/*
	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 84; j < 88; j++ )
		{
			for ( int k = 340; k < 345; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}






// Atlantic ocean

// Angola-dome at the end of the southern equatorial counter-current at the surface

// Angola-dome ( zum SECC gehörend, i=im-1 until i=im compares to 0 until -200m depth )
// Angola-dome ( from j=93 until j=96 compares to 3°S until 6°S,
//                        from k=330 until k=355 compares to 30°W until 5°W )

	for ( int i = i_SCC_u; i < i_ECC_o - 1; i++ )
	{
		for ( int j = 93; j < 97; j++ )
		{
			for ( int k = 330; k < 355; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					{
						u.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
						v.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
//						w.x[ i ][ j ][ k ] =  + IC_water * ( double ) ( i - i_SCC_u ) / ( double ) ( i_ECC_o - 1 - i_SCC_u );
					}
				}
			}
		}
	}
*/
}







void IC_Thermohalin::IC_DeepWater ( Array &h, Array &u, Array &v, Array &w, Array &c )
{
// initial conditions for v- and w-velocity components in deep flows
// Atlantic ocean
// Thermohaline Conveyor Belt

// downwards flow to deep flow     			  ( from j=32 until j=36 compares to 58°N until 54°N,
// Greenland until middel atlantic ridge       from k=317 until k=327 compares to 43°W until 33°W )


	for ( int i = i_bottom; i < im; i++ )
	{
		for ( int j = 32; j < 37; j++ )
		{
			for ( int k = 317; k < 328; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
					u.x[ i ][ j ][ k ] = - .1 * IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
//					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_deep - i_deep );
//					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_deep - i_deep );
				}
			}
		}
	}






// Atlantic ocean
// Thermohaline Conveyor Belt

// southerly directed current as deep flow ( from j=32 until j=50 compares to 58°N until 40°N,
// Greenland until middel atlantic ridge       from k=317 until k=327 compares to 43°W until 33°W )


	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 32; j < 51; j++ )
		{
			for ( int k = 317; k < 328; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
//					c.x[ i ][ j ][ k ] = 1.11;
//					u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}





// Atlantic ocean
// Thermohaline Conveyor Belt

// south of Newfoundland along the middle atlantic ridge

// south-west directed flow as deep flow ( from j=45 until j=60 compares to 30°N until 45°N,
//                                                                   from k=306 until k=320 compares to 54°W until 40°W )

	j_beg = 45;
	j_end = 61;
	k_beg = 306;
	k_end = 321;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 11;


	while ( ( j_end - j_run ) >= j_beg && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = ( j_end - j_run ); j > ( j_end - j_step - j_run ); j-- )
		{
			for ( int k = ( k_beg + k_run ); k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = i_bottom; i < i_deep + 1; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						c.x[ i ][ j ][ k ] = 1.11;
//						u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					}
				}
			}
		k_run++;
		}
	j_run++;
	}




// Atlantic

// parallel to the north atlantic ridge

// south of Newfoundland Guayana-fault ( from j=60 until j=88 compares to 2°N until 30°N,
//                                                                 from k=307 until k=316 compares to 53°W until 44°W )

	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 60; j < 89; j++ )
		{
			for ( int k = 307; k < 317; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
//					u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					w.x[ i ][ j ][ k ] = + 0.;
				}
			}
		}
	}





// Atlantic

// Southamerica		north coasts

// below the Guayana-current as deep flow ( from j=80 until j=100 compares to 10°N until 10°S,
//                                                                     from k=307 until k=327 compares to 53°W until 33°W )

	j_beg = 80;
	j_end = 100;
	k_beg = 307;
	k_end = 328;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 13;


	while ( ( j_beg + j_run ) <= j_end && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = j_beg + j_run; j < ( j_beg + j_step + j_run ); j++ )
		{
			for ( int k = k_beg + k_run; k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = i_bottom; i < i_deep + 1; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						c.x[ i ][ j ][ k ] = 1.11;
//						u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						v.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					}
				}
			}
		k_run++;
		}
	j_run++;
	}




// Atlantic

// Southamerica		east coasts

// Brasil-current bending east ( from j=107 until j=111 compares to 30°S until 35°S,
//                                               from k=273 until k=298 compares to 53°W until 25°W )
/*
	for ( int i = 0; i < i_beg; i++ )
	{
		for ( int j = 107; j < 112; j++ )
		{
			for ( int k = 273; k < 299; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					v.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}
*/



// Atlantic

// Southamerica		east coasts

// Parallel to the Brasil-current as deep flow ( from j=99 until j=149 compares to 9°S until 59°S,
//                                                                      from k=334 until k=339 compares to 26°W until 21°W )

	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 99; j < 150; j++ )
		{
			for ( int k = 334; k < 340; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
//					u.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}




// Pacific ocean

// Peru-current ( Humboldt-current ) coming from circumpolar current ( from j=107 until j=132 compares to 30°S until 58°S,
//                                                                                                               from k=246 until k=258 compares to 70°W until 83°W )
/*
	j_beg = 107;
	j_end = 133;
	k_beg = 246;
	k_end = 259;

	v_grad = + 0.001;
	k_step = 8;

	k_a = k_b = 0;

	flip = 0;

	for ( int k = k_beg; k < k_end; k++ )
	{
		if ( h.x[ i_max ][ j_beg ][ k ] == 1. )
		{
			k_a = k;
			flip = 1;
		}
		if ( flip == 1 ) break;
	}

	flip = 0;

	for ( int j = j_beg+1; j < j_end; j++ )
	{
		for ( int k = k_beg; k < k_end; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. )
			{
				k_b = k;
				k_grad = k_b - k_a;
				if ( k_grad <= - 3 ) k_grad = - 2;
				flip = 1;
			}
		if ( flip == 1 ) break;
		}


		for ( int k = k_b; k > ( k_b - k_step ); k-- )
		{
			for ( int i = 5; i < i_beg+2; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v.x[ i ][ j ][ k ] = + v_grad * IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					w.x[ i ][ j ][ k ] = + v_grad * ( double ) ( k_grad ) * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
		k_a = k_b;
		flip = 0;
	}
*/


// Pacific ocean

// Peru-current ( Humboldt-current ) coming from circumpolar current ( from j=120 until j=148 compares to 30°S until 58°S,
//                                                                                                              from k=277 until k=290 compares to 70°W until 83°W )
	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 120; j < 149; j++ )
		{
			for ( int k = 277; k < 291; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
//					u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}





// Pacific ocean

// deep flow to west out of Peru-current ( Humboldt-current ) ( from j=120 until j=125 compares to 30°S until 35°S,
// parallel to the equator                                                             from k=190 until k=285 compares to 75°W until 170°W )

	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 120; j < 126; j++ )
		{
			for ( int k = 190; k < 286; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
//					u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}




// Pacific ocean

// Newseeland in the east coming from circumpolar current ( from j=92 until j=152 compares to 2°S until 62°S,
// in direction to Japan                                                             from k=188 until k=197 compares to 172°W until 163°W )

	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 92; j < 153; j++ )
		{
			for ( int k = 188; k < 198; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
//					u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}





// Pacific ocean

// equator until Japan as deep flow ( from j=68 until j=93 compares to 3°S until 22°N,
//                                                        from k=162 until k=193 compares to 167°W until 162°O )

	j_beg = 68;
	j_end = 94;
	k_beg = 162;
	k_end = 194;
	j_run = 0;
	k_run = 0;
	j_step = 1;
	k_step = 10;


	while ( ( j_beg + j_run ) <= j_end && ( k_beg + k_run ) <= k_end )
	{
		for ( int j = j_beg + j_run; j < ( j_beg + j_step + j_run ); j++ )
		{
			for ( int k = k_beg + k_run; k < ( k_beg + k_step + k_run ); k++ )
			{
				for ( int i = i_bottom; i < i_deep + 1; i++ )
				{
					if ( h.x[ i ][ j ][ k ] != 1. )
					{
						c.x[ i ][ j ][ k ] = 1.11;
//						u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					}
				}
			}
		k_run++;
		}
	j_run++;
	}



// Pacific ocean

// east of Japan

// upward flow east of Japan ( from j=54 until j=67 compares to 23°N until 36°N
//                                              from k=162 until k=171 compares to 162°O until 171°O )

	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 54; j < 68; j++ )
		{
			for ( int k = 162; k < 172; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
					u.x[ i ][ j ][ k ] = + .1 * IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}





// Pacific ocean

// east of Japan

// upward flow east of Japan ( from j=50 until j=58 compares to 32°N until 40°N
//                                              from k=162 until k=210 compares to 162°O until 150°W )

	for ( int i = i_deep; i < im - 1; i++ )
	{
		for ( int j = 50; j < 59; j++ )
		{
			for ( int k = 162; k < 210; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
					u.x[ i ][ j ][ k ] = + .1 * IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
//					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
//					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
				}
			}
		}
	}



//	&&&&&&&&&&&&&&&&&&&&&&&&&&&		Indic ocean		&&&&&&&&&&&&&&&&&&&&&&&&&&


// Indic ocean
/*
// south east of Southafrica		east coasts

// south east of Southafrica coming from  the circumpolar current ( from j=83 until j=103 compares to 3°S until 26°S,
// prolongation of the circumpolar current                                         from k=55 until k=58 compares to 62°O until 65°O )

	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 83; j < 104; j++ )
		{
			for ( int k = 55; k < 59; k++ )
			{
			if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
//					u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}
*/



// Indic ocean

// south east of Southafrica		east coasts

// south east of Southafrica coming from  the circumpolar current ( from j=93 until j=152 compares to 3°S until 62°S,
//                                                                                                       from k=55 until k=60 compares to 55°O until 60°O )

	for ( int i = i_bottom; i < i_deep + 1; i++ )
	{
		for ( int j = 93; j < 153; j++ )
		{
			for ( int k = 55; k < 61; k++ )
			{
			if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
//					u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
//					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					//w.x[ i ][ j ][ k ] = IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
				}
			}
		}
	}




// Indic ocean

// south east of Southafrica		east coasts

// upwards flow north of Madagaskar ( from j=93 until j=97 compares to 3°S until 7°S,
//                                                           from k=55 until k=65 compares to 55°O until 65°O )

	for ( int i = i_deep; i < im - 2; i++ )
	{
		for ( int j = 93; j < 98; j++ )
		{
			for ( int k = 55; k < 66; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
					u.x[ i ][ j ][ k ] = + .1 * IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
//					v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
//					w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
				}
			}
		}
	}







// Indic ocean

// couth of India along the equator

// deep flow south of India ( from j=93 until j=97 compares to 3°S until 7°S,
// upwards flow                    from k=64 until k=95 compares to 64°O until 95°O )

	for ( int i = i_deep; i < im - 2; i++ )
	{
		for ( int j = 93; j < 98; j++ )
		{
			for ( int k = 64; k < 96; k++ )
			{
				if ( h.x[ i ][ j ][ k ] != 1. )
				{
					c.x[ i ][ j ][ k ] = 1.11;
					u.x[ i ][ j ][ k ] = + .1 * IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_deep ) / ( double ) ( i_max - i_deep );
				}
			}
		}
	}





// Indic ocean

// west of Corea along the equator

// upwards flow south of India ( from j=84 until j=89 compares to 5°S until 10°S,
// 90°W surface flow                  from k=78 until k=82 compares to 88°O until 92°O )
/*
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 84; j < 90; j++ )
		{
			for ( int k = 78; k < 83; k++ )
			{
				{
					{
						u.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						w.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					}
				}
			}
		}
	}
*/




// Indic ocean

// south of Corea

// 90°W surface flow south of India ( from j=84 until j=96 compares to 5°S until 18°S,
//                                                        from k=82 until k=85 compares to 92°O until 96°O )
/*
	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 84; j < 97; j++ )
		{
			for ( int k = 82; k < 86; k++ )
			{
				{
					{
						v.x[ i ][ j ][ k ] = + IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					}
				}
			}
		}
	}
*/



// Indic ocean

// south-east of Southafric		east coasts

// 90°W deep flow from circumpolar current ( from j=89 until j=132 compares to 10°S until 59°S,
//                                                                      from k=74 until k=78 compares to 83°O until 88°O )
/*
	for ( int i = 0; i < i_beg; i++ )
	{
		for ( int j = 89; j < 133; j++ )
		{
			for ( int k = 74; k < 79; k++ )
			{
				{
					{
						v.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
						w.x[ i ][ j ][ k ] = - IC_water * ( double ) ( i - i_bottom ) / ( double ) ( i_deep - i_bottom );
					}
				}
			}
		}
	}
*/
}




void IC_Thermohalin::BC_Surface_Temperature ( const string &Name_SurfaceTemperature_File, Array_2D &t_j, Array &t )
{
// initial conditions for the temperature and salinity at the sea surface

	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 3 );
	cout.setf ( ios::fixed );


// reading data from file Name_SurfaceTemperature_File_Read
	ifstream Name_SurfaceTemperature_File_Read;
	Name_SurfaceTemperature_File_Read.open ( Name_SurfaceTemperature_File.c_str(), ios_base::in );
	Name_SurfaceTemperature_File_Read.seekg ( 0L, ios::beg );
	anfangpos_1 = Name_SurfaceTemperature_File_Read.tellg ();


	if ( Name_SurfaceTemperature_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: could be opened" << endl;
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: begins at ::::::: " << anfangpos_1 << endl;
	}

	j = 0;
	k = 0;


	while ( ( k < km ) && ( !Name_SurfaceTemperature_File_Read.eof() ) )
	{
		while ( j < jm )
		{
			Name_SurfaceTemperature_File_Read >> dummy_1;
			Name_SurfaceTemperature_File_Read >> dummy_2;
			Name_SurfaceTemperature_File_Read >> dummy_3;

			t.x[ i_max ][ j ][ k ] = t_j.y[ j ][ k ] = ( dummy_3 + 273.15 ) / 273.15;

			j++;
		}
	j = 0;
	k++;
	}


//	cout << " ***** Ausdruck from 2D-Feldern ***** " << endl;
//	t_j.printArray_2D();
//	cout << " ***** printout of fields ***** " << endl;
//	t.printArray();

// Ende Lesen from Name_SurfaceTemperature_File

	Name_SurfaceTemperature_File_Read.seekg ( 0L, ios::end );
	endpos_1 = Name_SurfaceTemperature_File_Read.tellg ();

// Abschlussanweisungen für den Dateiabschluss (Dateiverwaltung)

	cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: ends at ::::::::: " << endpos_1 << endl;
	cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: has the length of ::::: " << endpos_1 - anfangpos_1 << " Bytes"<< endl;

// Im Falle eines Lesefehlers

	if ( Name_SurfaceTemperature_File_Read == NULL )
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: does not exist ::::::::: " << endl << endl << endl;
	}

	Name_SurfaceTemperature_File_Read.close();

	if ( Name_SurfaceTemperature_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: could be closed" << endl;
		cout << endl;
	}

	if ( Name_SurfaceTemperature_File_Read.fail() )
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: could not be closed" << endl;

// Ende Lesen from Name_SurfaceTemperature_File_Read
}






void IC_Thermohalin::BC_Surface_Salinity ( const string &Name_SurfaceSalinity_File, Array_2D &c_j, Array &c )
{
// initial conditions for the salinity at the sea surface

	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 3 );
	cout.setf ( ios::fixed );


// reading data from file Name_SurfaceSalinity_File_Read
	ifstream Name_SurfaceSalinity_File_Read;
	Name_SurfaceSalinity_File_Read.open ( Name_SurfaceSalinity_File.c_str(), ios_base::in );
	Name_SurfaceSalinity_File_Read.seekg ( 0L, ios::beg );
	anfangpos_2 = Name_SurfaceSalinity_File_Read.tellg ();


	if ( Name_SurfaceSalinity_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: could be opened" << endl;
		cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: begins at ::::::: " << anfangpos_2 << endl;
	}

	j = 0;
	k = 0;

	while ( ( k < km ) && ( !Name_SurfaceSalinity_File_Read.eof() ) )
	{
		while ( j < jm )
		{
			Name_SurfaceSalinity_File_Read >> dummy_1;
			Name_SurfaceSalinity_File_Read >> dummy_2;
			Name_SurfaceSalinity_File_Read >> dummy_3;

			if ( dummy_3 < 0. ) dummy_3 = 0.;

			c.x[ i_max ][ j ][ k ] = c_j.y[ j ][ k ] = dummy_3 / 38.8 / 1000.;
			j++;

//			cout << "\n***** Name_SurfaceSalinity_File_Read:   Länge = " << dummy_1 << "  Breite = " << dummy_2 << "  Tiefe = " << dummy_3 << endl;
//			cout << "***** Name_SurfaceSalinity_File_Read:   c = " << c.x[ i_max ][ j ][ k ] << "  j = " << j << "  k = " << k << endl;
		}
	j = 0;
	k++;
	}

//	cout << " ***** Ausdruck from 2D-Feldern ***** " << endl;
//	c_j.printArray_2D();
//	cout << " ***** printout of fields ***** " << endl;
//	c.printArray();

// Ende Lesen from Name_SurfaceSalinity_File

	Name_SurfaceSalinity_File_Read.seekg ( 0L, ios::end );
	endpos_2 = Name_SurfaceSalinity_File_Read.tellg ();

// Abschlussanweisungen für den Dateiabschluss (Dateiverwaltung)

	cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: ends at ::::::::: " << endpos_2 << endl;
	cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: has the length of ::::: " << endpos_2 - anfangpos_2 << " Bytes"<< endl;

// Im Falle eines Lesefehlers

	if ( Name_SurfaceSalinity_File_Read == NULL )
	{
		cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: does not exist ::::::::: " << endl << endl << endl;
	}

	Name_SurfaceSalinity_File_Read.close();

	if ( Name_SurfaceSalinity_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: could be closed" << endl;
		cout << endl;
	}

	if ( Name_SurfaceSalinity_File_Read.fail() )
		cout << "***** file ::::: " << Name_SurfaceSalinity_File << " ::::: could not be closed" << endl;

// Ende Lesen from Name_SurfaceSalinity_File_Read
}





void IC_Thermohalin::BC_Surface_Pressure ( double pa, double const gr, double const r_0_water, const double p_0, const double t_0, Array_2D &p_j, Array_2D &t_j, Array &h, Array &p, Array &t )
{
// boundary condition of surface pressure given by surface temperature through gas equation

//	rR = r_0_air * R_Air;
	rg = r_0_water * gr;
/*
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
//			p_j.y[ j ][ k ]  =  ( rR * t_j.y[ j ][ k ] * ( t_0 + 20. ) ) / p_0;		// pressure distribution on the surface in hPa / hPa
			p_j.y[ j ][ k ]  =  1.;		// pressure distribution on the surface in hPa / hPa
			p.x[ im-1 ][ j ][ k ] = p_j.y[ j ][ k ];

			if ( h.x[ im-1 ][ j ][ k ] == 1. )   
			{
				p_j.y[ j ][ k ]  =  0.;																				//  pressure distribution in solid surroundings
				p.x[ im-1 ][ j ][ k ] = p_j.y[ j ][ k ];
			}
		}
	}
*/

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = im-1; i >= 0; i-- )
			{
				if ( h.x[ i ][ j ][ k ] == 1. )   p.x[ i ][ j ][ k ] = pa;
				else
				{
					d_i = ( double ) i;
					p.x[ i ][ j ][ k ] = ( rg * ( d_i_max - d_i ) / 100. + p_0 ) / p_0;				// step size in 150 m   given in hPa / hPa
				}
			}
		}
	}
//	cout << " ***** Ausdruck from 2D-Feldern ***** " << endl;
//	p_j.printArray_2D();
//	cout << " ***** printout of fields ***** " << endl;
//	p.printArray();

}

