/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to print results on screen
*/


#include <iostream>
#include <cmath>

#include "Print_Atm.h"

using namespace std;




Print_Atmosphere::Print_Atmosphere ( int im, int jm, int km, int nm, int n, double time )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> nm = nm;
	this -> n = n;
	this -> time = time;
}


Print_Atmosphere::~Print_Atmosphere() {}


void Print_Atmosphere::printData ( Array &t, Array &u, Array &v, Array &w, Array &c, Array &p )
{
	cout.precision ( 4 );
	cout.setf ( ios::fixed );

	for ( int i = 0; i < im; i++ )
	{
		cout << "  phi = k-direction ======>  time = " << time << "   time block = " << n << "  r-level = " << i+1 << endl;
		cout << endl;

		for ( int j = 0; j < jm; j+=4 )
		{

//			cout << "SkalarFunktion c\n" << endl;

			for ( int k = 0; k < km; k+=16 )
			{
				cout.width ( 7 );
				cout.fill ( ' ' );
				cout << w.x[ i ][ j ][ k ] << " ";
			}
			cout << "\n";

		}
		cout << "\n\n";
	}
	cout << "\n";
}


void Print_Atmosphere::printIntro ( Array_1D &rad, Array_1D &the, Array_1D &phi )
{
	cout.precision ( 0 );
	cout.setf ( ios::fixed );

	cout << "\n\n";
	cout << "  radial pressure computation\n";
	cout << "\n\n";
	cout << "            ";
	cout << "radial position\n\n";
	cout << "       ";

	for ( int i = 1; i < im; i++ )
	{
		cout <<  i << "      ";
	}

	cout.precision ( 2 );
	cout.setf ( ios::fixed );

	cout << "\n\n";
	cout << "            ";
	cout << "radial steps in i-direction\n\n";
	cout << "        ";

	for ( int i = 0; i < im; i++ )
	{
		cout << rad.z[ i ] << "   ";
	}

	cout.precision ( 0 );
	cout.setf ( ios::fixed );

	cout << "\n\n";
	cout << "            ";
	cout << "theta steps in j-direction\n\n";
	cout << "        ";

	for ( int j = 0; j < jm; j++ )
	{
		cout << the.z[ j ] * 180. / M_PI << "     ";
	}

	cout << "\n\n";
	cout << "            ";
	cout << "phi steps in k-direction\n\n";
	cout << "        ";

	for ( int k = 0; k < km; k++ )
	{
		cout << phi.z[ k ] * 180. / M_PI << "     ";
	}
	cout << "\n\n\n\n" ;
	cout << "  time steps located on top of each other in " << nm;
	cout << "  theta-phi-r-time-blocks ( i-j-k-n ).\n\n";
	cout << "  r = i-direction consits of " << im;
	cout << "  theta-phi-blocks ( j-k ).\n\n";
	cout << "  theta = j-direction showing downwards\n\n";

	n = n + 1;
	cout << "  phi = k-direction showing to the right =========>---( time block = " << n << " ) \n\n\n";

	return;

}


