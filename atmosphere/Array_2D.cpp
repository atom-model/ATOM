/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 2D arrays
*/


#include <iostream>

#include "Array_2D.h"

using namespace std;



Array_2D::Array_2D ( )
{
	jm = 181;
	km = 361;

	y = new double*[ jm ];

	for ( int j = 0; j < jm; j++ )
	{
		y[ j ] = new double[ km ];
	}

// arbitrary initialisation of the z-field
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			y[ j ][ k ] = 222.;
		}
	}

}




Array_2D::~Array_2D ( )
{
	for ( int j = 0; j < jm; j++ )
	{
		delete [  ] y[ j ];
	}
	delete [  ] y;
}





void Array_2D::initArray_2D ( int jm, int km, double bb )
{
// initialisation of the y-field
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			y[ j ][ k ] = bb;
		}
	}
}




void Array_2D::printArray_2D ( )
{
	cout.precision ( 3 );
	cout.setf ( ios::fixed );

	cout << endl;
	cout << "  phi = k-direction ======>  theta = j-direction downwards :::::::::: r-level " << endl;
	cout << endl;

	for ( int j = 0; j < jm; j+=4 )
//	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k+=20 )
//		for ( int k = 0; k < km; k++ )
		{
			cout.width ( 4 );
			cout.fill( ' ' );

// 			cout << y[ j ][ k ] << " (" << &y[ j ][ k ] << ")" << " ";
			cout << y[ j ][ k ] << " ";
		}
	cout << endl;
	}
	cout << endl;
}



