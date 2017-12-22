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

#include "BC_Thermohalin.h"

using namespace std;

BC_Thermohalin::BC_Thermohalin ( int im, int jm, int km, int i_beg, int i_max, int Ma, int Ma_max, int Ma_max_half, double dr, double g, double r_0_water, double ua, double va, double wa, double ta, double ca, double pa, double u_0, double p_0, double t_0, double c_0, double cp_w, double L_hyd, double t_average, double t_cretaceous_max, double t_equator, double t_pole, const string &input_path )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> i_max = i_max;
	this -> i_beg = i_beg;
	this -> Ma = Ma;
	this -> Ma_max = Ma_max;
	this -> Ma_max_half = Ma_max_half;
	this -> dr = dr;
	this -> g = g;
	this -> r_0_water = r_0_water;
	this -> ua = ua;
	this -> va = va;
	this -> wa = wa;
	this -> ta = ta;
	this -> ca = ca;
	this -> pa = pa;
	this -> u_0 = u_0;
	this -> p_0 = p_0;
	this -> t_0 = t_0;
	this -> c_0 = c_0;
	this -> cp_w = cp_w;
	this -> L_hyd = L_hyd;
	this -> t_average = t_average;
	this -> t_cretaceous_max = t_cretaceous_max;
	this -> t_equator = t_equator;
	this -> t_pole = t_pole;
	this->input_path = input_path;


	i_bottom = 8;
	i_deep = i_bottom + 17;
	i_half = i_beg + 6;
	i_middle = i_beg - 7;
	j_half = ( jm - 1 ) / 2;
	j_half_1 = j_half + 1;

	d_i_half = ( double ) i_half;
	d_i_max = ( double ) ( im - 1 );


// reduction or amplification of flow velocities along coasts
// for the artificial initial and boundary conditions
	IC_water = .0015;				// reduced radial velocity u ( average velocity compares to          u_0 * IC_water = 0,45 * IC_water = 0,000675 m/s )


// ocean surface velocity is about 3% of the wind velocity at the surface
	water_wind = .03;

	pi180 = 180./M_PI;

// Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
// a further turning downwards until the end of the shear layer such that finally 90° of turning are reached
	Ekman_angle = 45.0 / pi180;
	Ekman_angle_add = 4.5 / pi180;



// array "dr_var" for Thomas algorithm
	dr_var = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		dr_var[ l ] = 0.;
	}

}


BC_Thermohalin::~BC_Thermohalin()
{
	delete [  ] dr_var;
}





void BC_Thermohalin::IC_v_w_Atmosphere ( Array &h, Array &u, Array &v, Array &w )
{
// initial conditions for v and w velocity components at the sea surface
// surface velocity reduced to 3% of the wind velocity ( water_wind )
	for ( int i = i_beg; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
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
				if ( h.x[ i ][ j ][ k ] == 0. )
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
}






void BC_Thermohalin::IC_v_w_WestEastCoast ( Array &h, Array &u, Array &v, Array &w, Array &c )
{
// initial conditions for v and w velocity components at the sea surface close to east or west coasts
// reversal of v velocity component between north and south equatorial current ommitted at respectively 10°
// w component unchanged

// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: east coast
	k_grad = 8;																			// extension of velocity change

	k_water = 0;																		// on water closest to coast
	k_sequel = 1;																		// on solid ground

	for ( int j = 0; j < j_half_1; j++ )													// outer loop: latitude
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
					for ( int i = i_middle; i <= i_half; i++ )					// loop in radial direction, extension for u -velocity component, downwelling here
					{
						m = i_half - i;
						d_i = ( double ) i;
						u.x[ i ][ j ][ k + l ] = - d_i / d_i_half * IC_water / ( ( double )( l + 1 ) );			// increase with depth, decrease with distance from coast
						u.x[ m ][ j ][ k + l ] = - d_i / d_i_half * IC_water / ( ( double )( l + 1 ) );			// decrease with depth, decrease with distance from coast
					}
				}
				k_sequel = 1;															// looking for another east coast
			}
		}																						// end of longitudinal loop
		k_water = 0;																	// starting at another latitude
	}																							// end of latitudinal loop



// southern hemisphere: east coast
	k_water = 0;
	k_sequel = 1;

	for ( int j = j_half_1; j < jm; j++ )
	{
		for ( int k = 0; k < km - k_grad; k++ )
		{
			if ( h.x[ i_max ][ j ][ k ] == 1. ) k_sequel = 0;

			if ( ( h.x[ i_max ][ j ][ k ] == 0. ) && ( k_sequel == 0 ) ) k_water = 0;
			else k_water = 1;

			if ( ( h.x[ i_max ][ j ][ k ] == 0. ) && ( k_water == 0 ) )
			{
				for ( int l = 0; l < k_grad; l++ )
				{
					for ( int i = i_middle; i <= i_half; i++ )
					{
						m = i_half - i;
						d_i = ( double ) i;
						u.x[ i ][ j ][ k + l ] = - d_i / d_i_half * IC_water / ( ( double )( l + 1 ) );			// increase with depth, decrease with distance from coast
						u.x[ m ][ j ][ k + l ] = - d_i / d_i_half * IC_water / ( ( double )( l + 1 ) );			// decrease with depth, decrease with distance from coast
					}
				}
				k_sequel = 1;
			}
		}
		k_water = 0;
	}



// search for west coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: west coast
	k_grad = 8;																			// extension of velocity change

	k_water = 0;																		// somewhere on water
	flip = 0;																				// somewhere on water

	for ( int j = 0; j < j_half_1; j++ )													// outer loop: latitude
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
					for ( int i = i_middle; i <= i_half; i++ )					// loop in radial direction, extension for u -velocity component, upwelling here
					{
						m = i_half - i;
						d_i = ( double ) i;
						u.x[ i ][ j ][ l ] = + d_i / d_i_half * IC_water / ( ( double )( k - l + 1 ) );			// increase with depth, decrease with distance from coast
						u.x[ m ][ j ][ l ] = + d_i / d_i_half * IC_water / ( ( double )( k - l + 1 ) );			// decrease with depth, decrease with distance from coast
					}
				}
				flip = 1;
			}
		}
		flip = 0;
	}



// southern hemisphere: west coast
	k_water = 0;
	flip = 0;

	for ( int j = j_half_1; j < jm; j++ )
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
					for ( int i = i_middle; i <= i_half; i++ )
					{
						m = i_half - i;
						d_i = ( double ) i;
						u.x[ i ][ j ][ l ] = + d_i / d_i_half * IC_water / ( ( double )( k - l + 1 ) );
						u.x[ m ][ j ][ l ] = + d_i / d_i_half * IC_water / ( ( double )( k - l + 1 ) );
					}
				}
				flip = 1;
			}
		}
		flip = 0;
	}
}







void BC_Thermohalin::BC_Temperature_Salinity ( Array &h, Array &t, Array &c, Array &p_dyn )
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

	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

// temperature-distribution by Ruddiman approximated by a parabola
	t_coeff = t_pole - t_equator;
	t_cretaceous_coeff = t_cretaceous_max / ( ( double ) Ma_max_half - ( double ) ( Ma_max_half * Ma_max_half / Ma_max ) );   // in °C
	t_cretaceous = t_cretaceous_coeff * ( double ) ( - ( Ma * Ma ) / Ma_max + Ma );   // in °C
	if ( Ma == 0 ) 	t_cretaceous = 0.;


	cout.precision ( 3 );

	time_slice_comment = "      time slice of Cretaceous-AGCM:";
	time_slice_number = " Ma = ";
	time_slice_unit = " million years";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << time_slice_comment << resetiosflags ( ios::left ) << setw ( 6 ) << fixed << setfill ( ' ' ) << time_slice_number << setw ( 3 ) << Ma << setw ( 12 ) << time_slice_unit << endl << endl;


	temperature_comment = "      temperature increase at cretaceous times: ";
	temperature_gain = " t increase";
	temperature_modern = "      mean temperature at modern times: ";
	temperature_cretaceous = "      mean temperature at cretaceous times: ";
	temperature_average = " t modern";
	temperature_average_cret = " t cretaceous";
	temperature_unit =  "°C ";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << temperature_comment << resetiosflags ( ios::left ) << setw ( 12 ) << temperature_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << t_cretaceous << setw ( 5 ) << temperature_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_modern << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average << setw ( 5 ) << temperature_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_cretaceous << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average_cret  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average + t_cretaceous << setw ( 5 ) << temperature_unit << endl;


	c_average = ( t_average + 346. ) / 10.;									// in psu, relation taken from "Ocean Circulation, The Open University"
	c_cretaceous = ( t_average + t_cretaceous + 346. ) / 10.;		// in psu
	c_cretaceous = c_cretaceous - c_average;
	if ( Ma == 0 ) 	c_cretaceous = 0.;

	salinity_comment = "      salinity increase at cretaceous times: ";
	salinity_gain = " c increase";
	salinity_modern = "      mean salinity at modern times: ";
	salinity_cretaceous = "      mean salinity at cretaceous times: ";
	salinity_average = " c modern";
	salinity_average_cret = " c cretaceous";
	salinity_unit =  "psu ";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << salinity_comment << resetiosflags ( ios::left ) << setw ( 12 ) << salinity_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << c_cretaceous << setw ( 5 ) << salinity_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << salinity_modern << resetiosflags ( ios::left ) << setw ( 13 ) << salinity_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << c_average << setw ( 5 ) << salinity_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << salinity_cretaceous << resetiosflags ( ios::left ) << setw ( 13 ) << salinity_average_cret  << " = "  << setw ( 7 )  << setfill ( ' ' ) << c_average + c_cretaceous << setw ( 5 ) << salinity_unit << endl;
	cout << endl;


	t_cretaceous = ( t_cretaceous + t_average + t_0 ) / t_0 - ( ( t_average + t_0 ) / t_0 );    // non-dimensional

	if ( Ma >= 0 )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				if ( h.x[ im-1 ][ j ][ k ] == 0. )
				{
					d_j = ( double ) j;
					if ( Ma == 0 )			t.x[ im-1 ][ j ][ k ] = t.x[ im-1 ][ j ][ k ] + t_cretaceous; // paleo surface temperature
					else						t.x[ im-1 ][ j ][ k ] = t_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + t_pole + t_cretaceous; // paleo surface temperature

					t_Celsius = t.x[ im-1 ][ j ][ k ] * t_0 - t_0;

					if ( t_Celsius < 4. )
					{
						t.x[ im -1 ][ j ][ k ] = ta + t_cretaceous;
						t_Celsius = ( ta + t_cretaceous ) * t_0 - t_0;										// water temperature not below 4°C
					}

					c.x[ im-1 ][ j ][ k ] = ( ( t_Celsius + 346. ) / 10. + c_cretaceous ) / c_0;	// non-dimensional

				}
				else
				{
					t.x[ im-1 ][ j ][ k ] = ta + t_cretaceous;
					c.x[ im-1 ][ j ][ k ] = ca + c_cretaceous / c_0;
				}
			}
		}
	}



// distribution of t and c with increasing depth till i_beg valid for all time slices including the modern world
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
		for ( int i = i_beg; i < im-1; i++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					t.x[ i ][ j ][ k ] = ( t.x[ im-1 ][ j ][ k ] - ta ) * ( double ) ( i - i_beg ) / ( double ) ( i_max - i_beg ) + ta + t_cretaceous;
					t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;
					c.x[ i ][ j ][ k ] = ( ( t_Celsius + 346. ) / 10. + c_cretaceous ) / c_0;

					if ( t_Celsius < 4. )
					{
						t.x[ i ][ j ][ k ] = ta + t_cretaceous;
						c.x[ i ][ j ][ k ] = ca + c_cretaceous / c_0;
					}
				}
				else
				{
					t.x[ i ][ j ][ k ] = ta + t_cretaceous;
					c.x[ i ][ j ][ k ] = ca + c_cretaceous / c_0;
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
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					t.x[ i ][ j ][ k ] = ta + t_cretaceous;
					c.x[ i ][ j ][ k ] = ca + c_cretaceous / c_0;
				}
			}
		}
	}

}






void BC_Thermohalin::BC_Pressure ( Array &p_stat, Array &t, Array &h )
{
// hydrostatic pressure distribution at surface

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			p_stat.x[ im-1 ][ j ][ k ] = p_0 / 1000.;		// given in bar
		}
	}


// hydrostatic pressure distribution increasing with depth

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = im-2; i >= 0; i-- )
			{
				d_i = ( double ) ( im - 1 - i );
				p_stat.x[ i ][ j ][ k ] = r_0_water * g * d_i * ( L_hyd / ( double ) ( im-1 ) ) / 100000. + p_0 / 1000.;				// current water pressure in bar
			}
		}
	}

}






void BC_Thermohalin::BC_Surface_Temperature_NASA ( const string &Name_SurfaceTemperature_File, Array &t )
{
	// initial conditions for the temperature and salinity at the sea surface

	cout.precision ( 3 );
	cout.setf ( ios::fixed );

	ifstream Name_SurfaceTemperature_File_Read;
	Name_SurfaceTemperature_File_Read.open(Name_SurfaceTemperature_File);

	if (!Name_SurfaceTemperature_File_Read.is_open()) {
		cerr << "ERROR: could not open SurfaceTemperature_File file at " << Name_SurfaceTemperature_File << "\n";
		abort();
	}


	j = 0;
	k = 0;

	while ( ( k < km ) && ( !Name_SurfaceTemperature_File_Read.eof() ) ) {
		while ( j < jm ) {
			Name_SurfaceTemperature_File_Read >> dummy_1;
			Name_SurfaceTemperature_File_Read >> dummy_2;
			Name_SurfaceTemperature_File_Read >> dummy_3;

			t.x[ im-1 ][ j ][ k ] = ( dummy_3 + 273.15 ) / 273.15;

			j++;
		}
		j = 0;
		k++;
	}

	Name_SurfaceTemperature_File_Read.close();
}






void BC_Thermohalin::BC_Surface_Salinity_NASA ( const string &Name_SurfaceSalinity_File, Array &c )
{
	// initial conditions for the salinity at the sea surface
	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 3 );
	cout.setf ( ios::fixed );

	// reading data from file Name_SurfaceSalinity_File_Read
	ifstream Name_SurfaceSalinity_File_Read;
	string path = input_path + Name_SurfaceSalinity_File;
	Name_SurfaceSalinity_File_Read.open(path);

	if (!Name_SurfaceSalinity_File_Read.is_open()) {
		cout << "could not read " << path << " at " << __FILE__ << " line " << __LINE__ << endl;
		abort();
	}

	j = 0;
	k = 0;

	while ( ( k < km ) && ( !Name_SurfaceSalinity_File_Read.eof() ) ) {
		while ( j < jm ) {
			Name_SurfaceSalinity_File_Read >> dummy_1;
			Name_SurfaceSalinity_File_Read >> dummy_2;
			Name_SurfaceSalinity_File_Read >> dummy_3;

			if ( dummy_3 < 0. ) dummy_3 = 0.;

			c.x[ im-1 ][ j ][ k ] = dummy_3 / 38.8 / 1000.;
			j++;
		}
		j = 0;
		k++;
	}

	Name_SurfaceSalinity_File_Read.close();
}
