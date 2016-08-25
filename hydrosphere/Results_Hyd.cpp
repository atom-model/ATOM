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
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;

	c43 = 4./3.;
	c13 = 1./3.;
}


Results_MSL_Hyd::~Results_MSL_Hyd () {}


void Results_MSL_Hyd::run_MSL_data ( double u_0, double c_0, Array &h, Array &u, Array &v, Array &w, Array &c, Array &Salt_Finger, Array &Salt_Diffusion, Array &Buoyancy_Force, Array_2D &Upwelling, Array_2D &Downwelling, Array_2D &SaltFinger, Array_2D &SaltDiffusion, Array_2D &BuoyancyForce, Array_2D &Salt_total, Array_2D &BottomWater )
{
// total upwelling as sum on normal velocity component values in a virtual vertical column

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Upwelling.y[ j ][ k ] = 0.;																// Upwelling
			Downwelling.y[ j ][ k ] = 0.;																// Downwelling
			BottomWater.y[ j ][ k ] = 0.;																// Bottom water
			SaltFinger.y[ j ][ k ] = 0.;																// SaltFinger
			SaltDiffusion.y[ j ][ k ] = 0.;																// SaltFinger
			BuoyancyForce.y[ j ][ k ] = 0.;																// Saltdiffusion
			Salt_total.y[ j ][ k ] = 0.;																// total Salt
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					if ( u.x[ i ][ j ][ k ] > 0. ) Upwelling.y[ j ][ k ] += u.x[ i ][ j ][ k ] * u_0;
					if ( u.x[ i ][ j ][ k ] < 0. ) Downwelling.y[ j ][ k ] += u.x[ i ][ j ][ k ] * u_0;
//					Buoyancy_Force.x[ i ][ j ][ k ] = fabs ( Buoyancy_Force.x[ i ][ j ][ k ] );
					SaltFinger.y[ j ][ k ] += Salt_Finger.x[ i ][ j ][ k ];
					SaltDiffusion.y[ j ][ k ] += Salt_Diffusion.x[ i ][ j ][ k ];
					BuoyancyForce.y[ j ][ k ] += Buoyancy_Force.x[ i ][ j ][ k ];
					Salt_total.y[ j ][ k ] += c.x[ i ][ j ][ k ] * c_0;
				}
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Downwelling.y[ j ][ k ] = fabs ( Downwelling.y[ j ][ k ] );
//			BuoyancyForce.y[ j ][ k ] = fabs ( BuoyancyForce.y[ j ][ k ] );
		}
	}


//	ideal age		corresponds to deep currents

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < 30; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					BottomWater.y[ j ][ k ] += sqrt ( u.x[ i ][ j ][ k ] * u.x[ i ][ j ][ k ] + v.x[ i ][ j ][ k ] * v.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k ] * w.x[ i ][ j ][ k ] ) * u_0;
				}
			}
		}
	}




// boundaries of buoyancy force

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
//			p_dyn.x[ 0 ][ j ][ k ] = c43 * p_dyn.x[ 1 ][ j ][ k ] - c13 * p_dyn.x[ 2 ][ j ][ k ];
//			p_dyn.x[ im-1 ][ j ][ k ] = c43 * p_dyn.x[ im-2 ][ j ][ k ] - c13 * p_dyn.x[ im-3 ][ j ][ k ];
//			p_dyn.x[ im-1 ][ j ][ k ] = 0.;

//			Buoyancy_Force.x[ 0 ][ j ][ k ] = c43 * Buoyancy_Force.x[ 1 ][ j ][ k ] - c13 * Buoyancy_Force.x[ 2 ][ j ][ k ];
			Buoyancy_Force.x[ 0 ][ j ][ k ] = 0.;
			Buoyancy_Force.x[ im-1 ][ j ][ k ] = c43 * Buoyancy_Force.x[ im-2 ][ j ][ k ] - c13 * Buoyancy_Force.x[ im-3 ][ j ][ k ];
//			Buoyancy_Force.x[ im-1 ][ j ][ k ] = 0.;
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int i = 0; i < im; i++ )
		{
//			Buoyancy_Force.x[ i ][ 0 ][ k ] = c43 * Buoyancy_Force.x[ i ][ 1 ][ k ] - c13 * Buoyancy_Force.x[ i ][ 2 ][ k ];
			Buoyancy_Force.x[ i ][ 0 ][ k ] = 0.;
			Buoyancy_Force.x[ i ][ jm-1 ][ k ] = c43 * Buoyancy_Force.x[ i ][ jm-2 ][ k ] - c13 * Buoyancy_Force.x[ i ][ jm-3 ][ k ];
		}
	}


	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Buoyancy_Force.x[ i ][ j ][ 0 ] = c43 * Buoyancy_Force.x[ i ][ j ][ 1 ] - c13 * Buoyancy_Force.x[ i ][ j ][ 2 ];
			Buoyancy_Force.x[ i ][ j ][ km-1 ] = c43 * Buoyancy_Force.x[ i ][ j ][ km-2 ] - c13 * Buoyancy_Force.x[ i ][ j ][ km-3 ];
			Buoyancy_Force.x[ i ][ j ][ 0 ] = Buoyancy_Force.x[ i ][ j ][ km-1 ] = ( Buoyancy_Force.x[ i ][ j ][ 0 ] + Buoyancy_Force.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}

}




void Results_MSL_Hyd::land_oceanFraction ( Array &h )
{
// calculation of the ratio ocean to land, also addition and substraction of CO2 of land, ocean and vegetation

	h_point_max =  ( jm - 1 ) * ( km - 1 );

	h_land = 0;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ im-1 ][ j ][ k ] == 1. )		h_land++;
		}
	}

	h_ocean = h_point_max - h_land;

	ozean_land = ( double ) h_ocean / ( double ) h_land;

	cout.precision ( 3 );

	cout << endl;
	cout << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      total number of points at constant hight " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_point_max << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the ocean surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_ocean << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the land surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_land << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      ocean/land ratio " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << ozean_land << endl << endl;
	cout << endl;
}






void Results_MSL_Hyd::show_MSL_data ( double c_0, Array &h, Array &c, Array &t, Array &p_dyn, Array &u, Array_2D &Upwelling, Array_2D &Downwelling, Array_2D &BottomWater, Array_2D &SaltFinger, Array_2D &BuoyancyForce, Array_2D &Salt_total  )
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


	Value_1 = Downwelling.y[ j_loc ][ k_loc ];
	Value_2 = Upwelling.y[ j_loc ][ k_loc ];
	Value_3 = BottomWater.y[ j_loc ][ k_loc ];
	Value_4 = SaltFinger.y[ j_loc ][ k_loc ];
	Value_5 = BuoyancyForce.y[ j_loc ][ k_loc ];
	Value_6 = Salt_total.y[ j_loc ][ k_loc ];

	cout << endl << endl << heading << endl << endl;

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_1 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_1 << setw ( 6 ) << name_unit_ms << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_2 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_2 << setw ( 6 ) << name_unit_ms << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_3 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_3 << setw ( 6 ) << name_unit_ms << endl << "                       " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_4 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_4 << setw ( 6 ) << name_unit_psu << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_5 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_5 << setw ( 6 ) << name_unit_psu << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_6 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_6 << setw ( 6 ) << name_unit_psu << endl << endl << endl;
}
