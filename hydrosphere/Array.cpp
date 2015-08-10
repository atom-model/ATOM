/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 3D arrays
*/


#include <iostream>

#include "Array.h"

using namespace std;



Array::Array ( int im, int jm, int km, double aa )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> aa = aa;
/*
	x = new double**[ im ]();

	for ( int i = 0; i < im; i++ )
	{
		x[ i ] = new double*[ jm ] ();
	}

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			x[ i ][ j ] = new double[ km ] ();
		}
	}
*/

	x = new double**[ im ];

	for ( int i = 0; i < im; i++ )
	{
		x[ i ] = new double*[ jm ];

		for ( int j = 0; j < jm; j++ )
		{
			x[ i ][ j ] = new double[ km ];
		}
	}





/*
#define HEIGHT 5
#define WIDTH 3
#define DEPTH 7

int main() {
  double ***p2DArray;

  // Allocate memory
  p2DArray = new double**[HEIGHT];
  for (int i = 0; i < HEIGHT; ++i) {
    p2DArray[i] = new double*[WIDTH];

    for (int j = 0; j < WIDTH; ++j)
      p2DArray[i][j] = new double[DEPTH];
  }

  // Assign values
  p2DArray[0][0][0] = 3.6;
  p2DArray[1][2][4] = 4.0;

  // De-Allocate memory to prevent memory leak
  for (int i = 0; i < HEIGHT; ++i) {
    for (int j = 0; j < WIDTH; ++j)
      delete [] p2DArray[i][j];

    delete [] p2DArray[i];
  }
  delete [] p2DArray;

  return 0;
}
*/








// initialisation of the x-field

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				x[ i ][ j ][ k ] = aa;
// 				cout << x[ i ][ j ][ k ] << " (" << &x[ i ][ j ][ k ] << ")" << "  ";
// 				cout << x[ i ][ j ][ k ] << "  ";
			}
// 			cout << endl;
		}
//		cout << endl;
	}
// 	cout << endl;
}



Array::~Array()
{
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			delete x[ i ][ j ];
		}
	}

	for ( int i = 0; i < im; i++ )
	{
		delete x[ i ];
	}

	delete [  ] x;
}







void Array::printArray()
{
	cout.precision ( 4 );
	cout.setf ( ios::fixed );

	for ( int i = im-1; i < im; i++ )
//	for ( int i = 1; i < im; i++ )
	{
		cout << "i = " << i << "   " << "im = " << im << "   " << "( transfer test of x in 'print_Array()' )" << endl;
		cout << endl;
		cout << "  phi = k-direction ======>  theta = j-direction downwards :::::::::: r-level = " << i << endl;
		cout << endl;

		for ( int j = 0; j < jm; j+=4 )
//		for ( int j = 0; j < 20; j++ )
		{
			for ( int k = 0; k < km; k+=16 )
//			for ( int k = 0; k < 24; k++ )
			{
				cout.width ( 7 );
//				cout.fill ( ' ' );
// 				cout << x[ i ][ j ][ k ] << " (" << &x[ i ][ j ][ k ] << ")" << "  ";
//				cout << x[ i ][ j ][ k ] << "  ";
				cout << x[ i ][ j ][ k ];
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;
}



