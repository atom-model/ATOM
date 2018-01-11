/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to build 3D arrays
*/


#ifndef _ARRAY_
#define _ARRAY_

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Int3DArray{
    private:
        int im, jm, km;

    public:
        std::vector<std::vector<std::vector<int> > > x;

    Int3DArray(int idim, int jdim, int kdim, double val){
        im=idim; jm=jdim; km=kdim;
        x.resize(idim);
        for(int i=0; i<idim; i++){
            x[i].resize(jdim);
            for(int j=0; j<jdim; j++){
                 x[i][j].resize(kdim, val);
            }
        }
    }

    ~Int3DArray ( ){}
};

class Array
{
private:
    int im, jm, km;

public:
    double ***x;

    Array(int idim, int jdim, int kdim, double val);
    ~Array ( );

    void printArray ( int, int, int );
    void initArray ( int, int, int, double );

    Int3DArray to_Int3DArray(){
        Int3DArray a(im,jm,km,0);
        for(int i=0; i<im; i++){
            for(int j=0; j<jm; j++){
                for(int k=0; k<km; k++){
                    a.x[i][j][k]=int(round(x[i][j][k]));
                }
            }
        }
        return a;
    }

};
#endif
