/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to search min/max values of variables
*/

#include <functional>
#include <iostream>
#include <cstring>

#include "Array.h"
#include "Array_2D.h"

#ifndef _MINMAX_
#define _MINMAX_

using namespace std;
namespace{
    std::function< double(double) > default_lambda=[](double i)->double{return i;};
}
class MinMax_Atm
{
    private:
        int im, jm, km, imax, jmax, kmax, imin, jmin, kmin;

        double maxValue, minValue;

    public:
        MinMax_Atm ( int, int );
        MinMax_Atm ( int, int, int );
        ~MinMax_Atm ();

        void searchMinMax_2D ( string , string , string , Array_2D &, Array &, double coeff=1.0);

        void searchMinMax_3D ( string , string , string , Array &, Array &, 
                               double coeff=1.0, 
                               std::function< double(double) > lambda = default_lambda,
                               bool print_heading=false );

        double out_maxValue (  ) const;

        double out_minValue (  ) const;
};
#endif
