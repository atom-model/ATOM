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




MinMax_Hyd::MinMax_Hyd ( int jm, int km, double c_0 )
{
	this-> jm = jm;
	this-> km = km;
	this-> c_0 = c_0;

	minValue = maxValue = 0.;
	jmax = kmax = jmin = kmin = 0;
}


MinMax_Hyd::MinMax_Hyd ( int im, int jm, int km, double c_0, double L_hyd )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> c_0 = c_0;
	this-> L_hyd = L_hyd;

	minValue = maxValue = 0.;
	imax = jmax = kmax = imin = jmin = kmin = 0;
}


MinMax_Hyd::~MinMax_Hyd () {}


void MinMax_Hyd::searchMinMax_3D ( string &name_maxValue, string &name_minValue, string &name_unitValue, Array &val_D, Array &h )
{
// search for minimum and maximum values of the 3-dimensional data sets
	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	heading_1 = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";
	heading_2 = " results based on three dimensional considerations of the problem";

	maxValue = 0.;
	imax = 0;
	jmax = 0;
	kmax = 0;
	imin = 0;
	jmin = 0;
	kmin = 0;

	Array &value_3D = val_D;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
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

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
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


	imax_level = imax * int ( L_hyd ) / ( im - 1 ) - int ( L_hyd );
	imin_level = imin * int ( L_hyd ) / ( im - 1 ) - int ( L_hyd );

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
		kmax_deg = kmax;
		deg_lon_max = deg_east;
	}

	if ( kmax > 180 )
	{
		kmax_deg = 360 - kmax;
		deg_lon_max = deg_west;
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
		kmin_deg = kmin;
		deg_lon_min = deg_east;
	}

	if ( kmin > 180 )
	{
		kmin_deg = 360 - kmin;
		deg_lon_min = deg_west;
	}

	cout.precision ( 6 );


	if ( name_maxValue == " max temperature " )
	{
		cout << endl << heading_1 << endl << heading_2 << endl << endl;

		maxValue = maxValue * 273.15 - 273.15;
		minValue = minValue * 273.15 - 273.15;
	}


// coefficient for units in hPa for the dynamic pressure 		coeff = r_0_water * u_0 * u_0 = 1026 * 15 * 15 * 0.01 = 207.77 hPa
	if ( name_maxValue == " max 3D pressure dynamic " )
	{
		maxValue = maxValue * 207.77;
		minValue = minValue * 207.77;
	}


	if ( name_maxValue == " max salt concentration " )
	{
		maxValue = maxValue * c_0;
		minValue = minValue * c_0;
	}

		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;
}




void MinMax_Hyd::searchMinMax_2D ( string &name_maxValue, string &name_minValue, string &name_unitValue, Array_2D &val, Array &h )
{
// search for minimum and maximum values of the 2-dimensional data sets on the sea surface

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	heading_1 = " printout of maximum and minimum values of properties at their locations: latitude, longitude";
	heading_2 = " results based on two dimensional considerations of the problem";

	maxValue = 0.;

	Array_2D &value = val;

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
		kmax_deg = kmax;
		deg_lon_max = deg_east;
	}

	if ( kmax > 180 )
	{
		kmax_deg = 360 - kmax;
		deg_lon_max = deg_west;
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
		kmin_deg = kmin;
		deg_lon_min = deg_east;
	}

	if ( kmin > 180 )
	{
		kmin_deg = 360 - kmin;
		deg_lon_min = deg_west;
	}


	cout.precision ( 6 );

	if ( name_maxValue == " max salt total " )
	{
		cout << endl << endl << heading_1 << endl << heading_2 << endl << endl;
	}

		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;

}


double MinMax_Hyd::out_maxValue (  ) const { return maxValue; }

double MinMax_Hyd::out_minValue (  ) const { return minValue; }


