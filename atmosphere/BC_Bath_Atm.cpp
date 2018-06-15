/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to read and prepare the bathymetric and topografic data
*/


#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <sstream>

#include "BC_Bath_Atm.h"

using namespace std;




BC_Bathymetry_Atmosphere::BC_Bathymetry_Atmosphere ( int NASATemperature, int im, int jm, int km, double co2_vegetation, double co2_land, double co2_ocean )
{
	this-> NASATemperature = NASATemperature;
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> co2_vegetation = co2_vegetation;
	this -> co2_land = co2_land;
	this -> co2_ocean = co2_ocean;
}



BC_Bathymetry_Atmosphere::~BC_Bathymetry_Atmosphere(){}



void BC_Bathymetry_Atmosphere::BC_MountainSurface ( string &Name_Bathymetry_File, double L_atm, Array_2D &Topography, Array_2D &value_top, Array &h, Array &aux_w )
{
	cout.precision ( 8 );
	cout.setf ( ios::fixed );

// default adjustment, h must be 0 everywhere
	h.initArray(im, jm, km, 0.);

// reading data from file Name_Bathymetry_File_Read
	ifstream Name_Bathymetry_File_Read;
	Name_Bathymetry_File_Read.open(Name_Bathymetry_File);

	if ( !Name_Bathymetry_File_Read.is_open() ) {
		cerr << "ERROR: could not open bathymetry file at " << Name_Bathymetry_File << "\n";
		abort();
	}

	double i_h = 0.;

	for (int j = 0; j < jm && !Name_Bathymetry_File_Read.eof(); j++) {
		for (int k = 0; k < km; k++) {
			Name_Bathymetry_File_Read >> dummy_1;
			Name_Bathymetry_File_Read >> dummy_2;
			Name_Bathymetry_File_Read >> dummy_3;

			if ( dummy_3 < 0. )
			{
				h.x[ 0 ][ j ][ k ] = Topography.y[ j ][ k ] = 0.;
			}
			 else
			{
				Topography.y[ j ][ k ] = dummy_3;
				i_h = round ( dummy_3 / ( L_atm / ( double) ( im - 1 ) ) );
				hight = ( int ) i_h;

				for ( int i = 0; i <= hight; i++ )
				{
					h.x[ i ][ j ][ k ] = 1.;
				}
			}
		}
	}

	Name_Bathymetry_File_Read.close();

// rewriting bathymetrical data from -180° _ 0° _ +180° coordinate system to 0°- 360°
	l = 0;

	for ( int k = 180; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			value_top.y[ j ][ l ] = Topography.y[ j ][ k ];
			for ( int i = 0; i < im; i++ )
			{
				aux_w.x[ i ][ j ][ l ] = h.x[ i ][ j ][ k ];
			}
		}
		l++;
	}

	l = l - 1;

	for ( int k = 0; k < 181; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			value_top.y[ j ][ l ] = Topography.y[ j ][ k ];
			for ( int i = 0; i < im; i++ )
			{
				aux_w.x[ i ][ j ][ l ] = h.x[ i ][ j ][ k ];
			}
		}
		l++;
	}



	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Topography.y[ j ][ k ] = value_top.y[ j ][ k ];
			for ( int i = 0; i < im; i++ )
			{
				h.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k ];
				aux_w.x[ i ][ j ][ k ] = 0.;
			}
		}
	}

}







void BC_Bathymetry_Atmosphere::BC_SolidGround ( int RadiationModel, int Ma, double g, double hp, double ep, double r_air, double R_Air, double t_0, double t_land, double t_cretaceous, double t_equator, double t_pole, double t_tropopause, double c_land, double c_tropopause, double co2_0, double co2_equator, double co2_pole, double co2_tropopause, double co2_cretaceous, double pa, double gam, double sigma, Array &h, Array &u, Array &v, Array &w, Array &t, Array &p_dyn, Array &c, Array &cloud, Array &ice, Array &co2, Array &radiation_3D, Array_2D &Vegetation )
{

// boundary conditions for the total solid ground
	j_half = ( jm -1 ) / 2;
	j_max = jm - 1;

	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

	int i_mount = 0;

// boundary conditions for solid ground areas
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = im-2; i >= 0; i-- )													// i = 0 a must: to keep velocities at surfaces to zero
			{
				if ( h.x[ i ][ j ][ k ] == 1. )
				{
					if ( i_mount == 0 ) 		i_mount = i;

					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;

					cloud.x[ i ][ j ][ k ] = 0.;
					ice.x[ i ][ j ][ k ] = 0.;

					p_dyn.x[ i ][ j ][ k ] = 0.;

					if ( NASATemperature == 0 )
					{
						t.x[ i ][ j ][ k ] = t.x[ i_mount ][ j ][ k ];
						c.x[ i ][ j ][ k ] = c.x[ i_mount ][ j ][ k ];				// water vapour amount above mount surface repeated
						co2.x[ i ][ j ][ k ] = co2.x[ i_mount ][ j ][ k ];		// co2 amount above mount surface repeated
					}

					if ( NASATemperature == 1 )
					{
						t.x[ i ][ j ][ k ] = t.x[ i_mount ][ j ][ k ];
						c.x[ i ][ j ][ k ] = c.x[ i_mount ][ j ][ k ];				// water vapour amount above mount surface repeated
						co2.x[ i ][ j ][ k ] = co2.x[ i_mount ][ j ][ k ];		// co2 amount above mount surface repeated
					}
				}
			}
			i_mount = 0;
		}
	}
}





void BC_Bathymetry_Atmosphere::vegetationDistribution ( double max_Precipitation, Array_2D &Precipitation, Array_2D &Vegetation, Array &t, Array &h )
{
// description or vegetation areas following the local dimensionsles values of precipitation, maximum value is 1

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( t.x[ 0 ][ j ][ k ] >= 1. ) ) Vegetation.y[ j ][ k ] = Precipitation.y[ j ][ k ] / max_Precipitation;			// actual vegetation areas
			else Vegetation.y[ j ][ k ] = 0.;
			if ( max_Precipitation <= 0. ) Vegetation.y[ j ][ k ] = 0.;
		}
	}
}





void BC_Bathymetry_Atmosphere::land_oceanFraction ( Array &h )
{
// calculation of the ratio ocean to land, also addition and subtraction of CO2 of land, ocean and vegetation

	h_point_max =  ( jm - 1 ) * ( km - 1 );

	h_land = 0;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ 0 ][ j ][ k ] == 1. )		h_land = h_land + h.x[ 0 ][ j ][ k ];
		}
	}

	h_ocean = h_point_max - h_land;

	ozean_land = ( double ) h_ocean / ( double ) h_land;

	cout.precision ( 3 );

	cout << endl;
	cout << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      total number of points at constant hight " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_point_max << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the ocean surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_ocean << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the land surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_land << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      ocean/land ratio " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << ozean_land << endl << endl;


	cout << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      addition of CO2 by ocean surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << co2_ocean << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      addition of CO2 by land surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << co2_land << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      subtraction of CO2 by vegetation " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << co2_vegetation << endl << setiosflags ( ios::left ) << setw ( 50 ) << "      valid for one single point on the surface"<< endl << endl;
	cout << endl;
}







