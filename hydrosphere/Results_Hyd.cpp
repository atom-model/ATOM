/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to produce resulting data on mean sea level
*/


#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>

#include "Results_Hyd.h"

using namespace std;


Results_MSL_Hyd::Results_MSL_Hyd ( int im, int jm, int km )
: 	co2_vegetation ( 0.16667 ),										// value valid for 100/600Gt per year on the total area
	co2_ocean ( 0.001 ),													// value valid for 0.6/600Gt per year on the total area
	co2_land ( 0.00033 )													// value valid for 0.2/600Gt per year on the total area
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
}


Results_MSL_Hyd::~Results_MSL_Hyd () {}


void Results_MSL_Hyd::run_MSL_data ( double u_0, double c_0, Array &hc, Array &uc, Array &vc, Array &wc, Array &cc, Array &SF, Array &SD, Array_2D &Up, Array_2D &Do, Array_2D &Sf, Array_2D &Sd, Array_2D &St, Array_2D &Bw )
{
// total upwelling as sum on normal velocity component values in a virtual vertical column

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Up.y[ j ][ k ] = 0.;																// Upwelling
			Do.y[ j ][ k ] = 0.;																// Downwelling
			Bw.y[ j ][ k ] = 0.;																// Bottom water
			Sf.y[ j ][ k ] = 0.;																// SaltFinger
			Sd.y[ j ][ k ] = 0.;																// Saltdiffusion
			St.y[ j ][ k ] = 0.;																// total Salt
		}
	}



	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( hc.x[ i ][ j ][ k ] == 0. )
				{
					if ( uc.x[ i ][ j ][ k ] > 0. ) Up.y[ j ][ k ] += uc.x[ i ][ j ][ k ] * u_0;
					if ( uc.x[ i ][ j ][ k ] < 0. ) Do.y[ j ][ k ] += uc.x[ i ][ j ][ k ] * u_0;
					Sf.y[ j ][ k ] += SF.x[ i ][ j ][ k ] * c_0;
					Sd.y[ j ][ k ] += SD.x[ i ][ j ][ k ] * c_0;
					St.y[ j ][ k ] += cc.x[ i ][ j ][ k ] * c_0;
				}
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Do.y[ j ][ k ] = fabs ( Do.y[ j ][ k ] );
			Sd.y[ j ][ k ] = fabs ( Sd.y[ j ][ k ] );
		}
	}


//	ideal age		corresponds to deep currents

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < 30; i++ )
			{
				if ( hc.x[ i ][ j ][ k ] == 0. )
				{
					Bw.y[ j ][ k ] += sqrt ( uc.x[ i ][ j ][ k ] * uc.x[ i ][ j ][ k ] + vc.x[ i ][ j ][ k ] * vc.x[ i ][ j ][ k ] + wc.x[ i ][ j ][ k ] * wc.x[ i ][ j ][ k ] ) * u_0;
				}
			}
		}
	}
}




void Results_MSL_Hyd::land_oceanFraction ( Array &hc )
{
// calculation of the ratio ocean to land, also addition and substraction of CO2 of land, ocean and vegetation

	h_point_max =  ( jm - 1 ) * ( km - 1 );

	h_land = 0;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( hc.x[ im-1 ][ j ][ k ] == 1. )		h_land++;
		}
	}

	h_ocean = h_point_max - h_land;

 //																										for the purpose of testing
	co2_vegetation = co2_vegetation / ( double ) h_land;
	co2_ocean = co2_ocean / ( double ) h_ocean;
	co2_land = co2_land / ( double ) h_land;
	ozean_land = ( double ) h_ocean / ( double ) h_land;

	cout.precision ( 3 );

	cout << endl;
	cout << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      total number of points at constant hight " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_point_max << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the ocean surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_ocean << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the land surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_land << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      ocean/land ratio " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << ozean_land << endl << endl;

//	cout << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      addition of CO2 by ocean surface " << " = + " << resetiosflags ( ios::left ) << setw ( 7 ) << scientific << setfill ( ' ' ) << co2_ocean << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      addition of CO2 by land surface " << " = + " << resetiosflags ( ios::left ) << setw ( 7 ) << scientific << setfill ( ' ' ) << co2_land << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      substraction of CO2 by vegetation " << " = - " << resetiosflags ( ios::left ) << setw ( 7 ) << scientific << setfill ( ' ' ) << co2_vegetation << endl << setiosflags ( ios::left ) << setw ( 50 ) << "      valid for one single point on the surface"<< endl << endl;
	cout << endl;
}






void Results_MSL_Hyd::show_MSL_data ( double c_0, Array &hc, Array &cc, Array &tc, Array &pc, Array &uc, Array_2D &Up, Array_2D &Do, Array_2D &Bw, Array_2D &Sf, Array_2D &Sd, Array_2D &St  )
{
	cout.precision ( 2 );

// printout of surface data at one predefinded location

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	name_Value_1 = " downwelling ";
	name_Value_2 = " upwelling ";
	name_Value_3 = " bottom water ";
	name_Value_4 = " salt finger ";
	name_Value_5 = " salt diffusion ";
	name_Value_6 = " salt total ";

	name_unit_ms = " m/s";
	name_unit_psu = " psu";

	heading = " printout of surface data at predefinded locations: level, latitude, longitude";

	i_loc_level = 0;																		// only at sea level MSL, constant
	j_loc = 90;
	k_loc = 180;


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
		k_loc_deg = 180 - k_loc;
		deg_lon = deg_west;
	}

	if ( k_loc > 180 )
	{
		k_loc_deg = k_loc - 180;
		deg_lon = deg_east;
	}


	Value_1 = Do.y[ j_loc ][ k_loc ];
	Value_2 = Up.y[ j_loc ][ k_loc ];
	Value_3 = Bw.y[ j_loc ][ k_loc ];
	Value_4 = Sf.y[ j_loc ][ k_loc ];
	Value_5 = Sd.y[ j_loc ][ k_loc ];
	Value_6 = St.y[ j_loc ][ k_loc ];

	cout << endl << endl << heading << endl << endl;

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_1 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_1 << setw ( 6 ) << name_unit_ms << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_2 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_2 << setw ( 6 ) << name_unit_ms << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_3 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_3 << setw ( 6 ) << name_unit_ms << endl << "                       " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_4 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_4 << setw ( 6 ) << name_unit_psu << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_5 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_5 << setw ( 6 ) << name_unit_psu << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_6 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_6 << setw ( 6 ) << name_unit_psu << endl << endl << endl;
}
