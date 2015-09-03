/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
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

#include "MinMax_Atm.h"

using namespace std;




MinMax::MinMax ( int jm, int km, double coeff_mmWS )
{
	this-> jm = jm;
	this-> km = km;
	this-> coeff_mmWS = coeff_mmWS;

	minValue = maxValue = 0.;
	jmax = kmax = jmin = kmin = 0;
}


MinMax::MinMax ( int im, int jm, int km )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;

	minValue = maxValue = 0.;
	imax = jmax = kmax = imin = jmin = kmin = 0;
}


MinMax::~MinMax () {}


void MinMax::searchMinMax_3D ( string &name_maxValue, string &name_minValue, string &name_unitValue, Array &value_3D, Array &hc )
{
// search for minimum and maximum values of the 3-dimensional data sets

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	heading = " printout of maximum and minimum values of properties at their locations: latitude, longitude, level";

	maxValue = 0.;
	imax = 0;
	jmax = 0;
	kmax = 0;
	imin = 0;
	jmin = 0;
	kmin = 0;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( hc.x[ i ][ j ][ k ] == 0. )
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
	}

	minValue = maxValue;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( hc.x[ i ][ j ][ k ] == 0. )
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
	}
	if ( minValue < 0. ) minValue =0.;

	imax_level = imax * 500;
	imin_level = imin * 500;

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


	if ( name_maxValue == " max water vapour " )
	{
		cout << endl << heading << endl << endl;

		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue * 1000. << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue * 1000. << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;
	}

	if ( name_maxValue == " max co2 " )
	{
		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue * 280. << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue * 280. << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;
	}
}



void MinMax::searchMinMax ( string &name_maxValue, string &name_minValue, string &name_unitValue, Array_2D &value, Array &hc )
{
// search for minimum and maximum values of the 2-dimensional data sets on the sea surface

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	maxValue = 0.;
	jmax = 0;
	kmax = 0;
	jmin = 0;
	kmin = 0;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			if ( hc.x[ 0 ][ j ][ k ] == 0. )
			{
				if ( value.y[ j ][ k ] > maxValue ) 
				{
					maxValue = value.y[ j ][ k ];
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
			if ( hc.x[ 0 ][ j ][ k ] == 0. )
			{
				if ( value.y[ j ][ k ] < minValue ) 
				{
					minValue = value.y[ j ][ k ];
					jmin = j;
					kmin = k;
				}
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

	if ( name_maxValue == " max precipitation " )
	{
		name_unitValue = "g/kg";
		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue / coeff_mmWS << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue / coeff_mmWS << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;
	}


	if ( name_maxValue == " max precipitable water " )
	{
		name_unitValue = "g/kg";
		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue / coeff_mmWS << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue / coeff_mmWS << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;
	}
}


double MinMax::out_maxValue (  ) const
{
	return maxValue;
}

double MinMax::out_minValue (  ) const
{
	return minValue;
}

int MinMax::out_imax (  ) const
{
	return imax;
}

int MinMax::out_jmax (  ) const
{
	return jmax;
}

int MinMax::out_kmax (  ) const
{
	return kmax;
}

int MinMax::out_imin (  ) const
{
	return imin;
}

int MinMax::out_jmin (  ) const
{
	return jmin;
}

int MinMax::out_kmin(  ) const
{
	return kmin;
}
