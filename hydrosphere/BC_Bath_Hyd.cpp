/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the salt concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to read and prepare the bathymetric and topografic data
*/


#include <iostream>
#include <fstream>
#include <cmath>


#include "BC_Bath_Hyd.h"

using namespace std;




BC_Bathymetry_Hydrosphere::BC_Bathymetry_Hydrosphere ( int im, int jm, int km )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
}


BC_Bathymetry_Hydrosphere::~BC_Bathymetry_Hydrosphere() {}



void BC_Bathymetry_Hydrosphere::BC_SeaGround ( const string &Name_Bathymetry_File, double L_hyd, Array &h, Array &aux_w )
{
	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 8 );
	cout.setf ( ios::fixed );


// default adjustment, h must be 0 everywhere
		for ( k = 0; k < km; k++ )
		{
			for ( j = 0; j < jm; j++ )
			{
				for ( i = 0; i < im; i++ )
				{
					h.x[ i ][ j ][ k ] = 0.;					// default
				}
			}
		}


// reading data from file Name_Bathymetry_File_Read
	ifstream Name_Bathymetry_File_Read;
	Name_Bathymetry_File_Read.open ( Name_Bathymetry_File.c_str(), ios_base::in );
	Name_Bathymetry_File_Read.seekg ( 0L, ios::beg );
	anfangpos_1 = Name_Bathymetry_File_Read.tellg ();


	if ( Name_Bathymetry_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: could be opened" << endl;
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: begins at ::::::: " << anfangpos_1 << endl;
	}


	j = 0;
	k = 0;

	L_hyd = - L_hyd;

	while ( ( j < jm ) && ( !Name_Bathymetry_File_Read.eof() ) )
	{
		while ( k < km )
		{
			Name_Bathymetry_File_Read >> dummy_1;
			Name_Bathymetry_File_Read >> dummy_2;
			Name_Bathymetry_File_Read >> dummy_3;

			if ( dummy_3 > 0. ) dummy_3 = 0.;

			i = ( im -1 ) - im * ( dummy_3 / L_hyd );
			i_boden = i;

			for ( i = 0; i <= i_boden; i++ ) 
			{
				h.x[ i ][ j ][ k ] = 1.;
//				cout << "\n***** Name_Bathymetry_File_Read:   Länge = " << dummy_1 << "  Breite = " << dummy_2 << "  Tiefe = " << dummy_3 << endl;
//				cout << "***** Name_Bathymetry_File_Read:   i = " << i << "  i_boden = " << i_boden << "  j = " << j << "  k = " << k << endl;
			}
			k++;
		}
	k = 0;
	j++;
	}


// end reading Name_Bathymetry_File

	Name_Bathymetry_File_Read.seekg ( 0L, ios::end );
	endpos_1 = Name_Bathymetry_File_Read.tellg ();

// final file administration

	cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: ends at ::::::::: " << endpos_1 << endl;
	cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: has the length of ::::: " << endpos_1 - anfangpos_1 << " bytes"<< endl;

// in case of reading error

	if ( Name_Bathymetry_File_Read == NULL )
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: not yet exists ::::::::: " << endl << endl << endl;
	}

	Name_Bathymetry_File_Read.close();

	if ( Name_Bathymetry_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: could be closed" << endl;
		cout << endl;
	}

	if ( Name_Bathymetry_File_Read.fail() )
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: could not be closed" << endl;

// end reading Name_Bathymetry_File_Read






// Umschreiben der bathymetrischen Daten von -180° - 0° - +180° Koordinatnesystem auf 0°- 360°

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

// Ende Umschreiben Bathymetrie
}





void BC_Bathymetry_Hydrosphere::BC_SolidGround ( double ca, double ta, double pa, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &pn_dyn, Array &cn )
{
// boundary conditions for the total solid ground

	for ( int i = 0; i < im-1; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 1. )
				{
					p_dyn.x[ i ][ j ][ k ] =  pn_dyn.x[ i ][ j ][ k ] = pa;
					t.x[ i ][ j ][ k ] = tn.x[ i ][ j ][ k ] = ta;
					c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] = ca;
					u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}

}
