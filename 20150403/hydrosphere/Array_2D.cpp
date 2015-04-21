/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
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



Array_2D::Array_2D ( int jm, int km, double bb )
{
	this -> jm = jm;
	this -> km = km;
	this -> bb = bb;

	y = new double*[ jm ];

	for ( int j = 0; j < jm; j++ )
	{
		y[ j ] = new double[ km ];
	}



// initialisation of the y-field

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			y[ j ][ k ] = bb;
// 			cout << y[ j ][ k ] << " (" << &y[ j ][ k ] << ")" << "  ";
// 			cout << y[ j ][ k ] << "  ";
		}
// 		cout << endl;
	}
// 	cout << endl;
	}



Array_2D::~Array_2D()
{
	for ( int j = 0; j < jm; j++ )
	{
		delete y[ j ];
	}

	delete [  ] y;
}





void Array_2D::printArray_2D()
{
	cout.precision ( 4 );
	cout.setf ( ios::fixed );

	cout << "  phi = k-direction ======>  theta = j-direction downwards :::::::::: r-level " << endl;
	cout << endl;

	for ( int j = 0; j < jm; j+=4 )
//	for ( int j = 0; j < 20; j++ )
	{
		for ( int k = 0; k < km; k+=16 )
		{
			cout.width ( 6 );
			cout.fill( ' ' );

// 			cout << y[ j ][ k ] << " (" << &y[ j ][ k ] << ")" << " ";
			cout << y[ j ][ k ] << " ";
		}
	cout << endl;
	}
	cout << endl;
}



