/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to search min/max values of variables
*/


#include <iostream>
#include <iomanip>
#include <cstring>

#include "MinMax_Hyd.h"

using namespace std;




MinMax::MinMax ( int jm, int km, double c_0 )
{
	this-> jm = jm;
	this-> km = km;
	this-> c_0 = c_0;

	minValue = maxValue = 0.;
	jmax = kmax = jmin = kmin = 0;
}


MinMax::MinMax ( int im, int jm, int km, double c_0, double L_hyd )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> c_0 = c_0;
	this-> L_hyd = L_hyd;

	minValue = maxValue = 0.;
	imax = jmax = kmax = imin = jmin = kmin = 0;
}


MinMax::~MinMax () {}


void MinMax::searchMinMax_3D ( string &name_maxValue, string &name_minValue, string &name_unitValue, Array &value_3D )
{
// search for minimum and maximum values of the 3-dimensional data sets

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	heading = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";

	maxValue = 0.;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( value_3D.x[ i ][ j ][ k ] > maxValue ) 
				{
					maxValue = value_3D.x[ i ][ j ][ k ];
					imax = i;
					jmax = j;
					kmax = k;
				}
			}
		}
	}

	minValue = maxValue;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( value_3D.x[ i ][ j ][ k ] < minValue ) 
				{
					minValue = value_3D.x[ i ][ j ][ k ];
					imin = i;
					jmin = j;
					kmin = k;
				}
			}
		}
	}
	if ( minValue < 0. ) minValue =0.;


	imax_level = - imax * int ( L_hyd ) / ( im - 1 );
	imin_level = - imin * int ( L_hyd ) / ( im - 1 );


//	maximum latitude and longitude units recalculated
	if ( jmax <= 90 )
	{
		jmax_deg = 90 - jmax;
		deg_lat_max = deg_north;
	}

	if ( jmax > 90 )
	{
		jmax_deg = jmax - 90;
		deg_lat_max = deg_south;
	}


	if ( kmax <= 180 )
	{
		kmax_deg = 180 - kmax;
		deg_lon_max = deg_west;
	}

	if ( kmax > 180 )
	{
		kmax_deg = kmax - 180;
		deg_lon_max = deg_east;
	}

//	minimum latitude and longitude units recalculated
	if ( jmin <= 90 )
	{
		jmin_deg = 90 - jmin;
		deg_lat_min = deg_north;
	}

	if ( jmin > 90 )
	{
		jmin_deg = jmin - 90;
		deg_lat_min = deg_south;
	}


	if ( kmin <= 180 )
	{
		kmin_deg = 180 - kmin;
		deg_lon_min = deg_west;
	}

	if ( kmin > 180 )
	{
		kmin_deg = kmin - 180;
		deg_lon_min = deg_east;
	}

	cout.precision ( 6 );


	cout << endl << heading << endl << endl;

	if ( name_unitValue == "psu" )
	{
		maxValue = maxValue * c_0;
		minValue = minValue * c_0;
	}

	cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;

}




void MinMax::searchMinMax ( string &name_maxValue, string &name_minValue, string &name_unitValue, Array_2D &value )
{
// search for minimum and maximum values of the 3-dimensional data sets on the sea surface

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	maxValue = 0.;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			if ( value.y[ j ][ k ] > maxValue ) 
			{
				maxValue = value.y[ j ][ k ];
				jmax = j;
				kmax = k;
			}
		}
	}

	minValue = maxValue;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			if ( ( value.y[ j ][ k ] < minValue ) && ( value.y[ j ][ k ] != 0. ) )
			{
				minValue = value.y[ j ][ k ];
				jmin = j;
				kmin = k;
			}
		}
	}
	if ( minValue < 0. ) minValue =0.;




	imax_level = 0;
	imin_level = 0;

//	maximum latitude and longitude units recalculated
	if ( jmax <= 90 )
	{
		jmax_deg = 90 - jmax;
		deg_lat_max = deg_north;
	}

	if ( jmax > 90 )
	{
		jmax_deg = jmax - 90;
		deg_lat_max = deg_south;
	}


	if ( kmax <= 180 )
	{
		kmax_deg = 180 - kmax;
		deg_lon_max = deg_west;
	}

	if ( kmax > 180 )
	{
		kmax_deg = kmax - 180;
		deg_lon_max = deg_east;
	}

//	minimum latitude and longitude units recalculated
	if ( jmin <= 90 )
	{
		jmin_deg = 90 - jmin;
		deg_lat_min = deg_north;
	}

	if ( jmin > 90 )
	{
		jmin_deg = jmin - 90;
		deg_lat_min = deg_south;
	}


	if ( kmin <= 180 )
	{
		kmin_deg = 180 - kmin;
		deg_lon_min = deg_west;
	}

	if ( kmin > 180 )
	{
		kmin_deg = kmin - 180;
		deg_lon_min = deg_east;
	}


	cout.precision ( 6 );

		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;

}


double MinMax::out_maxValue (  ) const { return maxValue; }

double MinMax::out_minValue (  ) const { return minValue; }

int MinMax::out_imax (  ) const { return imax; }

int MinMax::out_jmax (  ) const { return jmax; }

int MinMax::out_kmax (  ) const { return kmax; }

int MinMax::out_imin (  ) const { return imin; }

int MinMax::out_jmin (  ) const { return jmin; }

int MinMax::out_kmin(  ) const { return kmin; }
