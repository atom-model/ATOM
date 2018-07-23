/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 3D arrays
*/


#include <cassert>
#include <iostream>

#include "Array.h"

using namespace std;

#define MAXI 41
#define MAXJ 181
#define MAXK 361

Array::Array(int idim, int jdim, int kdim, double val):
    x(NULL) 
{
    initArray(idim, jdim, kdim, val);
}

Array::~Array ( )
{
    for ( int i = 0; i < im; i++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            delete [  ] x[ i ][ j ];
        }
    }

    for ( int i = 0; i < im; i++ )
    {
        delete [  ] x[ i ];
    }

    delete [  ] x;
}


void Array::initArray ( int im, int jm, int km, double aa )
{
    if(x){
        assert(im == this->im);
        assert(jm == this->jm);
        assert(km == this->km);

        for ( int i = 0; i < im; i++ )
        {
            for ( int j = 0; j < jm; j++ )
            {
                for ( int k = 0; k < km; k++ )
                {
                    x[ i ][ j ][ k ] = aa;
                }
            }
        }
    }else{
        assert(im <= MAXI);
        assert(jm <= MAXJ);
        assert(km <= MAXK);

        this->im = im;
        this->jm = jm;
        this->km = km;

        x = new double**[im];

        for ( int i = 0; i < im; i++ )
        {
            x[ i ] = new double*[jm];

            for ( int j = 0; j < jm; j++ )
            {
                x[ i ][ j ] = new double[km];
                for ( int k = 0; k < km; k++ )
                {
                    x[ i ][ j ][ k ] = aa;
                }
            }
        }    
    }
}


void Array::printArray ( int im, int jm, int km )
{
    assert(im == this->im);
    assert(jm == this->jm);
    assert(km == this->km);

    cout.precision ( 3 );
    cout.setf ( ios::fixed );

//      for ( int i = 0; i < im; i++ )
      for ( int i = 0; i <= 0; i++ )
    {
//        cout << "i = " << i << "   " << "im = " << im << "   " << "( transfer test of x in 'print_Array()' )" << endl;
//        cout << endl;
//        cout << "  phi = k-direction ======>  theta = j-direction downwards :::::::::: r-level = " << i << endl;
        cout << endl;

        for ( int j = 0; j < jm; j+=4 )
//    for ( int i = 0; i < im ; i++ )
//        for ( int j = 0; j < im; j++ )
        {
            for ( int k = 0; k < km; k+=20 )
//            for ( int k = 0; k < km; k++ )
            {
                cout.width ( 4 );
                cout.fill( ' ' );

//                 cout << x[ i ][ j ][ k ] << " (" << &x[ i ][ j ][ k ] << ")" << " ";
                cout << x[ i ][ j ][ k ] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}



