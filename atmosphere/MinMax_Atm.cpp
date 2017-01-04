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




MinMax_Atm::MinMax_Atm ( int jm, int km, double coeff_mmWS )
{
	this-> jm = jm;
	this-> km = km;
	this-> coeff_mmWS = coeff_mmWS;

	minValue = maxValue = 0.;
	jmax = kmax = jmin = kmin = 0;
}


MinMax_Atm::MinMax_Atm ( int im, int jm, int km )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;

	minValue = maxValue = 0.;
	imax = jmax = kmax = imin = jmin = kmin = 0;
}


MinMax_Atm::~MinMax_Atm () {}


void MinMax_Atm::searchMinMax_3D ( string &name_maxValue, string &name_minValue, string &name_unitValue, Array &val_D, Array &h )
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

	Array &value_D = val_D;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( value_D.x[ i ][ j ][ k ] > maxValue ) 
				{
					maxValue = value_D.x[ i ][ j ][ k ];
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
				if ( value_D.x[ i ][ j ][ k ] < minValue ) 
				{
					minValue = value_D.x[ i ][ j ][ k ];
					imin = i;
					jmin = j;
					kmin = k;
				}
			}
		}
	}


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


	if ( name_maxValue == " max 3D temperature " )
	{
		cout << endl << heading_1 << endl << heading_2 << endl << endl;

		maxValue = maxValue * 273.15 - 273.15;
		minValue = minValue * 273.15 - 273.15;
	}


	if ( name_maxValue == " max 3D water vapour " )
	{
		maxValue = maxValue * 1000.;
		minValue = minValue * 1000.;
	}

	if ( name_maxValue == " max 3D cloud water " )
	{
		maxValue = maxValue * 1000.;
		minValue = minValue * 1000.;
	}

	if ( name_maxValue == " max 3D cloud ice " )
	{
		maxValue = maxValue * 1000.;
		minValue = minValue * 1000.;
	}

	if ( name_maxValue == " max 3D rain " )
	{
		maxValue = maxValue * 1000.;
		minValue = minValue * 1000.;
	}

	if ( name_maxValue == " max 3D snow " )
	{
		maxValue = maxValue * 1000.;
		minValue = minValue * 1000.;
	}


// coefficient for units in hPa for the dynamic pressure 		coeff = 0.5 * r_Air * u_0 * u_0 = 0.5 * 1.2 * 15 * 15 = 135 * 500. * 0.01 = 675
	if ( name_maxValue == " max 3D pressure dynamic " )
	{
		maxValue = maxValue * 675.;
		minValue = minValue * 675.;
	}


	if ( name_maxValue == " max 3D co2 " )
	{
		maxValue = maxValue * 280.;
		minValue = minValue * 280.;
	}


		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;

}








void MinMax_Atm::searchMinMax_2D ( string &name_maxValue, string &name_minValue, string &name_unitValue, Array_2D &val, Array &h )
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
			if ( value.y[ j ][ k ] < minValue ) 
			{
				minValue = value.y[ j ][ k ];
				jmin = j;
				kmin = k;
			}
		}
	}


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

	if ( name_maxValue != " max co2_total " )
	{
		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;
	}


	if ( name_maxValue == " max co2_total " )
	{
		maxValue = 280. * maxValue;
		minValue = 280. * minValue;

		name_unitValue = "ppm";


		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;
	}

	if ( name_maxValue == " max precipitation " )
	{
		name_unitValue = "g/kg";

		cout << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_maxValue << " = " << resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << maxValue / coeff_mmWS << setw ( 6 ) << name_unitValue << setw ( 5 ) << jmax_deg << setw ( 3 ) << deg_lat_max << setw ( 4 ) << kmax_deg << setw ( 3 ) << deg_lon_max << setw ( 6 ) << imax_level << setw ( 2 ) << level << "          " << setiosflags ( ios::left ) << setw ( 26 ) << setfill ( '.' ) << name_minValue << " = "<< resetiosflags ( ios::left ) << setw ( 12 ) << fixed << setfill ( ' ' ) << minValue / coeff_mmWS << setw ( 6 ) << name_unitValue << setw ( 5 )  << jmin_deg << setw ( 3 ) << deg_lat_min << setw ( 4 ) << kmin_deg << setw ( 3 ) << deg_lon_min  << setw ( 6 ) << imin_level << setw ( 2 ) << level << endl;
	}


	if ( name_maxValue == " max precipitation_NASA " )
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






double MinMax_Atm::out_maxValue (  ) const
{
	return maxValue;
}

double MinMax_Atm::out_minValue (  ) const
{
	return minValue;
}

