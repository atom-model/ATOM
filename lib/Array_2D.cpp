/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 2D arrays
*/

#include <cassert>
#include <iostream>

#include "Array_2D.h"

using namespace std;

#define MAXJ 181
#define MAXK 361


Array_2D::Array_2D(int jdim, int kdim, double val) : y(NULL) {
    initArray_2D(jdim, kdim, val);
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
    if(!y){//when y is null
	    assert(jm <= MAXJ);
	    assert(km <= MAXK);

        this->jm = jm;
        this->km = km;

        y = new double*[jm];

        for ( int j = 0; j < jm; j++ )
        {
            y[ j ] = new double[km];
            for ( int k = 0; k < km; k++ )
            {
                y[ j ][ k ] = bb;
            }
        }
    }else{
        assert(jm == this->jm);
        assert(km == this->km);
        for ( int j = 0; j < jm; j++ )
        {
            for ( int k = 0; k < km; k++ )
            {
                y[ j ][ k ] = bb;
            }
        }
    }  
}

void Array_2D::printArray_2D ( int jm, int km )
{
	assert(jm == this->jm);
	assert(km == this->km);

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
