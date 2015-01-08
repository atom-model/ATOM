/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to surveil the accuracy of the iterations
*/


#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>

#include "Accuracy.h"

using namespace std;




Accuracy::Accuracy( int n, int im, int jm, int km, double dr, double dthe, double dphi )
{
	this-> n = n;
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> dr = dr;
	this-> dthe = dthe;
	this-> dphi = dphi;
}


Accuracy::Accuracy ( int im, int Ma, int n, int velocity_iter, int pressure_iter, double min, double L_atm )
{
	this-> im = im;
	this-> n = n;
	this-> velocity_iter = velocity_iter;
	this-> pressure_iter = pressure_iter;
	this-> min = min;
	this-> Ma = Ma;
	this-> L_atm = L_atm;

	i_u = j_u = k_u = i_v = j_v = k_v = i_w = j_w = k_w = i_t = j_t = k_t = i_c = j_c = k_c = i_co2 = j_co2 = k_co2 = i_p = j_p = k_p = 0;
	i_res = j_res = k_res = 0;
	residuum = max_u = max_v = max_w = max_t = max_c = max_co2 = max_p = 0.;
	residuum_alt = min = min_u = min_v = min_w = min_t = min_c = min_co2 = min_p = 0.;
}


Accuracy::~Accuracy () {}


double Accuracy::residuumQuery ( int &i_res, int &j_res, int &k_res, double &min, Array_1D &rad, Array_1D &the, Array &u, Array &v, Array &w )
{
// value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )

	min = 0.;

	for ( int i = 1; i < im-1; i++ )
	{
		rmsinthe = rad.z[ i ] * sinthe;

		for ( int j = 1; j < jm-1; j++ )
		{
			sinthe = sin( the.z[ j ] );
			costhe = cos( the.z[ j ] );

			for ( int k = 1; k < km-1; k++ )
			{
				dudr   = ( u.x[ i+1 ][ j ][ k ] - u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
				dvdthe = ( v.x[ i ][ j+1 ][ k ] - v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
				dwdphi = ( w.x[ i ][ j ][ k+1 ] - w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );

				residuum = dudr + 2. * u.x[ i ][ j ][ k ] / rad.z[ i ] + dvdthe / rad.z[ i ]
							+ costhe / rmsinthe * v.x[ i ][ j ][ k ] + dwdphi / rmsinthe;
				if ( fabs ( residuum ) >= min ) 
				{
					min = residuum;
					i_res = i;
					j_res = j;
					k_res = k;
				}
			}
		}
	}
	return 0;
}



double Accuracy::residuumQuery_2D ( int &j_res, int &k_res, double &min, Array_1D &rad, Array_1D &the, Array &v, Array &w )
{
// value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )

	min = 0.;

	rmsinthe = rad.z[ 0 ] * sinthe;

	for ( int j = 1; j < jm-1; j++ )
	{
		sinthe = sin( the.z[ j ] );
		costhe = cos( the.z[ j ] );

	for ( int k = 1; k < km-1; k++ )
		{
				dvdthe = ( v.x[ 0 ][ j+1 ][ k ] - v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );
				dwdphi = ( w.x[ 0 ][ j ][ k+1 ] - w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
				residuum = dvdthe / rad.z[ 0 ] + costhe / rmsinthe * v.x[ 0 ][ j ][ k ] + dwdphi / rmsinthe;
			if ( fabs ( residuum ) >= min ) 
			{
				min = residuum;
				j_res = j;
				k_res = k;
			}
		}
	}

	return 0;
}






double Accuracy::steadyQuery ( int &i_u, int &j_u, int &k_u, int &i_v, int &j_v, int &k_v, int &i_w, int &j_w, int &k_w, int &i_t, int &j_t, int &k_t, int &i_c, int &j_c, int &k_c, int &i_co2, int &j_co2, int &k_co2, int &i_p, int &j_p, int &k_p, double &min_u, double &min_v, double &min_w, double &min_t, double &min_c, double &min_co2, double &min_p, Array &u, Array &un, Array &v, Array &vn, Array &w, Array &wn, Array &t, Array &tn, Array &c, Array &cn, Array &co2, Array &co2n, Array &p, Array &pn )
{
// state of a steady solution ( min )
	min_u = 0.;
	min_v = 0.;
	min_w = 0.;
	min_t = 0.;
	min_c = 0.;
	min_p = 0.;

	for ( int i = 1; i < im-2; i++ )
	{
		for ( int j = 1; j < jm-2; j++ )
		{
			for ( int k = 1; k < km-2; k++ )
			{
				max_u = fabs ( ( u.x[ i ][ j ][ k ] - un.x[ i ][ j ][ k ] ) );
				if ( max_u >= min_u )
				{
					min_u = max_u;
					i_u = i;
					j_u = j;
					k_u = k;
				}

				max_v = fabs ( ( v.x[ i ][ j ][ k ] - vn.x[ i ][ j ][ k ] ) );
				if ( max_v >= min_v )
				{
					min_v = max_v;
					i_v = i;
					j_v= j;
					k_v = k;
				}

				max_w = fabs ( ( w.x[ i ][ j ][ k ] - wn.x[ i ][ j ][ k ] ) );
				if ( max_w >= min_w )
				{
					min_w = max_w;
					i_w = i;
					j_w = j;
					k_w = k;
				}

				max_t = fabs ( ( t.x[ i ][ j ][ k ] - tn.x[ i ][ j ][ k ] ) );
				if ( max_t >= min_t )
				{
					min_t = max_t;
					i_t = i;
					j_t = j;
					k_t = k;
				}

				max_c = fabs ( ( c.x[ i ][ j ][ k ] - cn.x[ i ][ j ][ k ] ) );
				if ( max_c >= min_c )
				{
					min_c = max_c;
					i_c = i;
					j_c = j;
					k_c = k;
				}

				max_co2 = fabs ( ( co2.x[ i ][ j ][ k ] - co2n.x[ i ][ j ][ k ] ) );
				if ( max_co2 >= min_co2 )
				{
					min_co2 = max_co2;
					i_co2 = i;
					j_co2 = j;
					k_co2 = k;
				}

				max_p = fabs ( ( p.x[ i ][ j ][ k ] - pn.x[ i ][ j ][ k ] ) );
				if ( max_p >= min_p )
				{
					min_p = max_p;
					i_p = i;
					j_p = j;
					k_p = k;
				}
			}
		}
	}
	return 0;
}





double Accuracy::steadyQuery_2D ( int &i_v, int &j_v, int &k_v, int &i_w, int &j_w, int &k_w, int &i_p, int &j_p, int &k_p, double &min_v, double &min_w, double &min_p, Array &v, Array &vn, Array &w, Array &wn, Array &p, Array &pn )
{
// state of a steady solution ( min )
	min_v = 0.;
	min_w = 0.;
	min_p = 0.;

	for ( int j = 2; j < jm-2; j++ )
	{
		for ( int k = 2; k < km-2; k++ )
		{
			max_v = fabs ( ( v.x[ 0 ][ j ][ k ] - vn.x[ 0 ][ j ][ k ] ) );
			if ( max_v >= min_v )
			{
				min_v = max_v;
				j_v= j;
				k_v = k;
			}

			max_w = fabs ( ( w.x[ 0 ][ j ][ k ] - wn.x[ 0 ][ j ][ k ] ) );
			if ( max_w >= min_w )
			{
				min_w = max_w;
				j_w = j;
				k_w = k;
			}

			max_p = fabs ( ( p.x[ 0 ][ j ][ k ] - pn.x[ 0 ][ j ][ k ] ) );
			if ( max_p >= min_p )
			{
				min_p = max_p;
				j_p = j;
				k_p = k;
			}
		}
	}
	return 0;
}




void Accuracy::iterationPrintout ( int &nm, int &velocity_iter_max, int &pressure_iter_max, int &i_res, int &j_res, int &k_res, int &i_u, int &j_u, int &k_u, int &i_v, int &j_v, int &k_v, int &i_w, int &j_w, int &k_w, int &i_t, int &j_t, int &k_t, int &i_c, int &j_c, int &k_c, int &i_co2, int &j_co2, int &k_co2, int &i_p, int &j_p, int &k_p, double &min_u, double &min_v, double &min_w, double &min_t, double &min_c, double &min_co2, double &min_p )
{
// statements on the convergence und iterational process
	cout.precision ( 6 );
	cout.setf ( ios::fixed );

	cout << endl << endl;
	cout << "      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	cout << "      3D iterational process" << endl;
	cout << "      max total iteration number nm = " << nm << endl;
	cout << "      outer pressure loop:  max iteration number pressure_iter_max = " << pressure_iter_max << endl;
	cout << "      inner velocity loop:  max iteration number velocity_iter_max = " << velocity_iter_max << endl << endl;

	cout << "      n = " << n << "     " << "velocity_iter = " << velocity_iter << "     " << "pressure_iter = " << pressure_iter<< "     " << "Ma = " << Ma << endl;
	cout << endl;


// printout of maximum and minimum absolute and relative errors of the computed values at their locations while iterating
	heading = " 3D iterational process for the surface boundary conditions\n printout of maximum and minimum absolute and relative errors of the computed values at their locations: level, latitude, longitude";

	cout << endl << endl << heading << endl << endl;

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	int choice = { 1 };

	preparation:

	switch ( choice )
	{
		case 1 :	name_Value = " residuum: continuity equation ";
						Value = min;
						i_loc = i_res;
						j_loc = j_res;
						k_loc = k_res;
						break;

		case 2 :	name_Value = " dp: pressure Poisson equation ";
						Value = min_p;
						i_loc = i_p;
						j_loc = j_p;
						k_loc = k_p;
						break;

		case 3 :	name_Value = " du: Navier Stokes equation ";
						Value = min_u;
						i_loc = i_u;
						j_loc = j_u;
						k_loc = k_u;
						break;

		case 4 :	name_Value = " dv: Navier Stokes equation ";
						Value = min_v;
						i_loc = i_v;
						j_loc = j_v;
						k_loc = k_v;
						break;

		case 5 :	name_Value = " dw: Navier Stokes equation ";
						Value = min_w;
						i_loc = i_w;
						j_loc = j_w;
						k_loc = k_w;
						break;

		case 6 :	name_Value = " dt: energy transport equation ";
						Value = min_t;
						i_loc = i_t;
						j_loc = j_t;
						k_loc = k_t;
						break;

		case 7 :	name_Value = " dc: vapour transport equation ";
						Value = min_c;
						i_loc = i_c;
						j_loc = j_c;
						k_loc = k_c;
						break;

		case 8 :	name_Value = " dco2: CO2 transport equation ";
						Value = min_co2;
						i_loc = i_co2;
						j_loc = j_co2;
						k_loc = k_co2;
						break;

		default : 	cout << choice << "error in iterationPrintout member function in class Accuracy" << endl;
	}
						i_loc_level = i_loc * int ( L_atm ) / ( im - 1 );

						if ( j_loc <= 90 )
						{
							j_loc_deg = 90 - j_loc;
							deg_lat = deg_north;
						}

						if ( j_loc > 90 )
						{
							j_loc_deg = j_loc - 90;
							deg_lat = deg_south;
						}


						if ( k_loc <= 180 )
						{
							k_loc_deg = k_loc;
							deg_lon = deg_east;
						}

						if ( k_loc > 180 )
						{
							k_loc_deg = 360 - k_loc;
							deg_lon = deg_east;
						}
	

	cout << setiosflags ( ios::left ) << setw ( 36 ) << setfill ( '.' ) << name_Value << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << Value << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << endl;

	choice++;
	if ( choice <= 8 ) goto preparation;

	cout << endl << endl;
}






void Accuracy::iterationPrintout_2D ( int &nm, int &velocity_iter_max_2D, int &pressure_iter_max_2D, int &j_res, int &k_res, int &i_v, int &j_v, int &k_v, int &i_w, int &j_w, int &k_w,  int &i_p, int &j_p, int &k_p, double &min_v, double &min_w, double &min_p )
{
// statements on the convergence und iterational process
	cout.precision ( 6 );
	cout.setf ( ios::fixed );

	cout << "      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	cout << "      2D iterational process" << endl;
	cout << "      max total iteration number nm = " << nm << endl;
	cout << "      outer pressure loop:  max iteration number pressure_iter_max_2D = " << pressure_iter_max_2D << endl;
	cout << "      inner velocity loop:  max iteration number velocity_iter_max_2D = " << velocity_iter_max_2D << endl << endl;

	cout << "      n = " << n << "     " << "velocity_iter = " << velocity_iter << "     " << "pressure_iter = " << pressure_iter<< "     " << "Ma = " << Ma << endl;
	cout << endl;


// printout of maximum and minimum absolute and relative errors of the computed values at their locations while iterating
	heading = " 2D iterational process for the surface boundary conditions\n printout of maximum and minimum absolute and relative errors of the computed values at their locations: level, latitude, longitude";

	cout << endl << endl << heading << endl << endl;

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	int choice = { 1 };

	preparation:

	switch ( choice )
	{
		case 1 :	name_Value = " residuum: continuity equation ";
						Value = min;
						j_loc = j_res;
						k_loc = k_res;
						break;

		case 2 :	name_Value = " dp: pressure Poisson equation ";
						Value = min_p;
						j_loc = j_p;
						k_loc = k_p;
						break;

		case 3 :	name_Value = " dv: Navier Stokes equation ";
						Value = min_v;
						j_loc = j_v;
						k_loc = k_v;
						break;

		case 4 :	name_Value = " dw: Navier Stokes equation ";
						Value = min_w;
						j_loc = j_w;
						k_loc = k_w;
						break;


		default : 	cout << choice << "error in iterationPrintout member function in class Accuracy" << endl;
	}

						if ( j_loc <= 90 )
						{
							j_loc_deg = 90 - j_loc;
							deg_lat = deg_north;
						}

						if ( j_loc > 90 )
						{
							j_loc_deg = j_loc - 90;
							deg_lat = deg_south;
						}


						if ( k_loc <= 180 )
						{
							k_loc_deg = k_loc;
							deg_lon = deg_east;
						}

						if ( k_loc > 180 )
						{
							k_loc_deg = 360 - k_loc;
							deg_lon = deg_west;
						}
	

	cout << setiosflags ( ios::left ) << setw ( 36 ) << setfill ( '.' ) << name_Value << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << Value << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon << endl;

	choice++;
	if ( choice <= 4 ) goto preparation;

	cout << endl << endl;
}


