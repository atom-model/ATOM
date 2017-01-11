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

#include "BC_Bath_Atm.h"

using namespace std;




BC_Bathymetry_Atmosphere::BC_Bathymetry_Atmosphere ( int im, int jm, int km, double co2_vegetation, double co2_land, double co2_ocean )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> co2_vegetation = co2_vegetation;
	this -> co2_land = co2_land;
	this -> co2_ocean = co2_ocean;
}



BC_Bathymetry_Atmosphere::~BC_Bathymetry_Atmosphere(){}



void BC_Bathymetry_Atmosphere::BC_MountainSurface ( string &Name_Bathymetry_File, double L_atm, Array &h, Array &aux_w )
{
	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 8 );
	cout.setf ( ios::fixed );

	// default adjustment, h must be 0 everywhere
	h.initArray(im, jm, km, 0.);

	// reading data from file Name_Bathymetry_File_Read
	ifstream Name_Bathymetry_File_Read;
	Name_Bathymetry_File_Read.open(Name_Bathymetry_File);

	if (!Name_Bathymetry_File_Read.is_open()) {
		cerr << "ERROR: could not open bathymetry file at " << Name_Bathymetry_File << "\n";
		abort();
	}

	for (int j = 0; j < jm && !Name_Bathymetry_File_Read.eof(); j++) {
		for (int k = 0; k < km; k++) {
			Name_Bathymetry_File_Read >> dummy_1;
			Name_Bathymetry_File_Read >> dummy_2;
			Name_Bathymetry_File_Read >> dummy_3;

			if ( dummy_3 < 0. ) {
				h.x[ 0 ][ j ][ k ] = 0.;
			} else {
				int i = ( im - 1 ) * ( dummy_3 / L_atm );
				i_SL = i;

				for ( i = 0; i <= i_SL; i++ ) {
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
			for ( int i = 0; i < im; i++ )
			{
				aux_w.x[ i ][ j ][ l ] = h.x[ i ][ j ][ k ];
			}
		}
		l++;
	}



	for ( int k = 0; k < 180; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
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
			for ( int i = 0; i < im; i++ )
			{
				h.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k ];
			}
		}
	}

	// end rewriting bathymetry
}





void BC_Bathymetry_Atmosphere::BC_IceShield ( int Ma, double t_0, Array &h, Array &t, Array &c, Array &IceLayer, Array_2D &Ice_Balance, Array_2D &Ice_Balance_add )
{
// computation of ice shield following the theorie by Milankowitsch

	min = 0.;
	max = 0.;

// prescribed reference level still to be fixed!

	Hoehe_equi_km = 3000;
	Hoehe_delta = 500;
	Hoehe_equi = Hoehe_equi_km / Hoehe_delta;

// limit of temperature for accumulation 1 and 2 and ablation

	t_equi_Celsius = - 8.;
	t_plus_Celsius = + 2.;
	t_pluss_Celsius = + 12.;
	t_minus_Celsius = - 15.;
	t_minuss_Celsius = - 30.;

	t_equi = ( t_equi_Celsius + t_0 ) / t_0;
	t_plus = ( t_plus_Celsius + t_0 ) / t_0;
	t_pluss = ( t_pluss_Celsius + t_0 ) / t_0;
	t_minus = ( t_minus_Celsius + t_0 ) / t_0;
	t_minuss = ( t_minuss_Celsius + t_0 ) / t_0;




	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( t.x[ 0 ][ j ][ k ] >= t_plus ) break;
			{
				for ( int i = 0; i < im; i++ )																		// ice balance performed over total hight
				{
					if ( ( t.x[ i ][ j ][ k ] >= t_minuss ) && ( t.x[ i ][ j ][ k ] <= t_plus ) )
					{
						Akkumulation_1 = 0.0125 * ( t.x[ i ][ j ][ k ] * t_0 ) - 2.9735;			//	formula valid from -30°C to +2°C by Ruddiman
					}
					else
					{
						Akkumulation_1 = 0.;
					}

					if ( ( t.x[ i ][ j ][ k ] > t_plus ) && ( t.x[ i ][ j ][ k ] <= t_pluss ) )
					{
						Akkumulation_2 = - 0.04 * ( t.x[ i ][ j ][ k ] * t_0 ) + 11.5;					//	formula valid from +2°C to +12°C by Ruddiman
					}
					else
					{
						Akkumulation_2 = 0.;
					}

					if ( ( t.x[ i ][ j ][ k ] >= t_minus ) && ( t.x[ i ][ j ][ k ] <= t_plus ) )
					{
						Ablation = 0.01384 * ( t.x[ i ][ j ][ k ] * t_0 ) * ( t.x[ i ][ j ][ k ] * t_0 ) - 7.14186 * ( t.x[ i ][ j ][ k ] * t_0 ) + 921.2458;	//	formula valid from -15°C to +2°C by Ruddiman
					}
					else
					{
						Ablation = 0.;
					}


					if ( ( t.x[ i ][ j ][ k ] >= t_minuss ) && ( t.x[ i ][ j ][ k ] <= t_plus ) )
					{
						Ice_Balance.y[ j ][ k ] = Akkumulation_1 - Ablation;						// ice balance for accumulation_1 minus ablation
					}


					if ( ( t.x[ i ][ j ][ k ] > t_plus ) && ( t.x[ i ][ j ][ k ] <= t_pluss ) )
					{
						Ice_Balance.y[ j ][ k ] = Akkumulation_2;										// ice balance for accumulation_2
					}

					Ice_Balance_add.y[ j ][ k ] = Ice_Balance_add.y[ j ][ k ] + Ice_Balance.y[ j ][ k ];		// average of ice balance out of the two parts across hight
				}
			}
		}
	}

// searching minimum and maximum values of mean ice balance

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( Ice_Balance_add.y[ j ][ k ] >= max ) 
			{
				max = Ice_Balance_add.y[ j ][ k ];
			}

			if ( Ice_Balance_add.y[ j ][ k ] <= min ) 
			{
				min = Ice_Balance_add.y[ j ][ k ];
			}
		}
	}


	Ice_Balance_add_diff = max - min;										// maximum crosswise extension of ice balance

// presentation of ice shield hight following the local dimensionless ice balance values, maximum value now 1

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Ice_Hoehe = ( Ice_Balance_add.y[ j ][ k ] - min ) / max * 2.;			// assumed actual ice shield hight
			i_Ice_lauf = ( int ) Ice_Hoehe;

			for ( int i = 0; i < i_Ice_lauf; i++ )
			{
				if ( ( t.x[ 0 ][ j ][ k ] >= t_plus ) || ( h.x[ 0 ][ j ][ k ] == 0. ) ) break;
				{
					IceLayer.x[ i ][ j ][ k ] = 1.;
				}
			}
		}
	}
}






void BC_Bathymetry_Atmosphere::BC_SolidGround ( int RadiationModel, int i_max, double g, double hp, double ep, double r_air, double R_Air, double t_0, double t_land, double t_cretaceous, double t_equator, double t_pole, double t_tropopause, double c_land, double c_tropopause, double co2_0, double co2_equator, double co2_pole, double co2_tropopause, double co2_cretaceous, double pa, double gam, Array &h, Array &u, Array &v, Array &w, Array &t, Array &p_dyn, Array &c, Array &cloud, Array &ice, Array &co2, Array_2D &Vegetation )
{

// boundary conditions for the total solid ground

	j_half = ( jm -1 ) / 2;
	j_max = jm - 1;

	d_i_max = ( double ) i_max;
	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

	t_co2_eff = t_pole - t_equator;
	co_co2_eff = co2_pole - co2_equator;


// velocity components, dynamioc pressure, temperature, water vapour and CO2-content as boundary condition at the sea surface and on ground

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( h.x[ 0 ][ j ][ k ] == 0. ) 
			{
				d_j = ( double ) j;
				co2.x[ 0 ][ j ][ k ] = ( co_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + co2_pole + co2_cretaceous + co2_ocean ) / co2_0; // non-dimensional
			}
			if ( h.x[ 0 ][ j ][ k ] == 1. ) 
			{
				d_j = ( double ) j;
				t.x[ 0 ][ j ][ k ] = t_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + t_pole + t_cretaceous + t_land;

//				c.x[ 0 ][ j ][ k ] = hp * ep * exp ( 17.0809 * ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) ) / ( ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01 );
//				c.x[ 0 ][ j ][ k ] = c_land * c.x[ 0 ][ j ][ k ];					// relativ water vapour contents on land reduced by factor

				co2.x[ 0 ][ j ][ k ] = ( co_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + co2_pole + co2_cretaceous + co2_land - co2_vegetation * Vegetation.y[ j ][ k ] ) / co2_0;	// parabolic distribution from pole to pole

				u.x[ 0 ][ j ][ k ] = 0.;
				v.x[ 0 ][ j ][ k ] = 0.;
				w.x[ 0 ][ j ][ k ] = 0.;
			}
		}
	}



// boundary conditions for solid ground areas
	for ( int i = 1; i < im-1; i++ )													// i = 0 a must: to keep velocities at surfaces to zero
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 1. )
				{
					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;

					p_dyn.x[ i ][ j ][ k ] = 0.;

					d_i = ( double ) i;

					t.x[ i ][ j ][ k ] = ( t_tropopause - t.x[ 0 ][ j ][ k ] ) / d_i_max * d_i + t.x[ 0 ][ j ][ k ];										// linear temperature decay up to tropopause

//					c.x[ i ][ j ][ k ] = c.x[ 0 ][ j ][ k ] - ( c_tropopause - c.x[ 0 ][ j ][ k ] ) * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) );	// radial distribution approximated by a parabola ( Weischet )

					co2.x[ i ][ j ][ k ] = co2.x[ 0 ][ j ][ k ] - ( ( co2_tropopause - co2.x[ 0 ][ j ][ k ] * co2_0 ) * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) ) ) / co2_0;
																																// radial distribution approximated by a parabola

				}
			}
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







