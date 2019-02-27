/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 2D arrays
*/

#ifndef _ARRAY_2D_
#define _ARRAY_2D_

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class Array_2D
{
private:
    int jm, km;

public:
    double **y;

    Array_2D(int jdim, int kdim, double val);
    ~Array_2D ( );

    Array_2D(): jm(0), km(0), y(NULL){}

    void printArray_2D ( int, int );
    void initArray_2D ( int, int, double );

    void save(const string& fn){
        std::ofstream os(fn + ".bin", std::ios::binary | std::ios::out);
        for(int j=jm-1; j>=0; j--){
            os.write(reinterpret_cast<const char*>(y[j]+(km/2)), std::streamsize((km/2)*sizeof(double)));
            os.write(reinterpret_cast<const char*>(y[j]), std::streamsize((km/2+1)*sizeof(double)));
        }
        os.close();
    }
};
#endif
