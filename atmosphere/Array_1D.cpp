/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 1D arrays
*/


#include <iostream>

#include "Array_1D.h"

using namespace std;


Array_1D::Array_1D ( int i, double cc, double z0, double dz )
{
	im = i;
	this -> cc = cc;
	this -> z0 = z0;
	this -> dz = dz;

	z = new double[ im ];


// initialisation of the z-field

	for ( int l = 0; l < im; l++ )
	{
		z[ l ] = cc;
//		cout << z[ l ] << " (" << &z[ l ] << ")" << "  ";
//		cout << z[ l ] << "  ";
	}
//	cout << endl;
//	cout << endl;
}



Array_1D::~Array_1D()
{
	delete [  ] z;
}



void Array_1D::Coordinates (  )
{
	z[ 0 ] = z0;

	for ( int i = 1; i < im; i++ )
	{
		z[ i ] = z[ i - 1 ] + dz;
	}

//	cout << " ***** Printout of 1D-array in Array_1D ***** " << endl;
//	printArray_1D();

}



void Array_1D::printArray_1D()
{
	cout.precision ( 6 );
	cout.setf ( ios::fixed );

	cout << endl;
	cout << " coordinate-direction " << endl;
	cout << endl;

	for ( int i = 0; i < im; i++ )
//	for ( int i = 0; i < im; i++ )
	{
		cout.width ( 6 );
		cout.fill( ' ' );

// 		cout << z[ i ] << " (" << &i[ i ] << ")" << " ";
		cout << z[ i ] << " ";
	}
cout << endl;
}



