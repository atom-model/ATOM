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
#include <cassert>

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

    Array(): im(0), jm(0), km(0), x(NULL)
    {}

    Array(const Array &a){
        im = a.im;
        jm = a.jm;
        km = a.km;

        x = new double**[im];

        for ( int i = 0; i < im; i++ )
        {
            x[ i ] = new double*[jm];

            for ( int j = 0; j < jm; j++ )
            {
                x[ i ][ j ] = new double[km];
                for ( int k = 0; k < km; k++ )
                {
                    x[ i ][ j ][ k ] = a.x[i][j][k];
                }
            }
        }
    }

    void printArray ( int, int, int );
    void initArray ( int, int, int, double );

    friend Array operator* (double coeff, const Array &a);

    void operator=(const Array &a){
        assert(this->im == a.im);
        assert(this->jm == a.jm);
        assert(this->km == a.km);
        for(int i=0; i<im; i++){
            for(int j=0; j<jm; j++){
                for(int k=0; k<km; k++){
                    this->x[i][j][k]=a.x[i][j][k];
                }
            }
        }
    }

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

    double max() const{
        assert(im && jm && km);
        double ret=x[0][0][0];
        for(int i=0; i<im; i++){
            for(int j=0; j<jm; j++){
                for(int k=0; k<km; k++){
                    ret=std::max(ret, x[i][j][k]);
                }
            }
        }
        return ret;
    }

    double min() const{
        assert(im && jm && km);
        double ret=x[0][0][0];
        for(int i=0; i<im; i++){
            for(int j=0; j<jm; j++){
                for(int k=0; k<km; k++){
                    ret=std::min(ret, x[i][j][k]);
                }
            }
        }
        return ret;
    }

    double mean() const {
        assert(im && jm && km);
        double ret=0;
        for(int i=0; i<im; i++){
            for(int j=0; j<jm; j++){
                for(int k=0; k<km; k++){
                    ret+=x[i][j][k];
                }
            }
        }
        return ret/(im*jm*km);
    }

    double max_2D() const{
        assert(jm && km);
        double ret=x[0][0][0];
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                ret=std::max(ret, x[0][j][k]);
            }
        }
        return ret;
    }

    double min_2D() const{
        assert(jm && km);
        double ret=x[0][0][0];
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                ret=std::min(ret, x[0][j][k]);
            }
        }
        return ret;
    }

    double mean_2D() const{
        assert(jm && km);
        double ret=0;
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                ret+=x[0][j][k];
            }
        }
        return ret/(jm*km);
    }

    bool has_nan() const{
        for ( int i = 0; i < im; i++ ){
            for ( int j = 0; j < jm; j++ ){
                for ( int k = 0; k < km; k++ ){
                    if(isnan(x[ i ][ j ][ k ])){
                        return true;
                    }
                }
            }
        }
        return false;
    }
};

inline Array operator* (double coeff, const Array &a){
    Array ret(a.im,a.jm,a.km, 0.);
    for(int i=0; i<a.im; i++){
        for(int j=0; j<a.jm; j++){
            for(int k=0; k<a.km; k++){
                ret.x[i][j][k]=coeff * a.x[i][j][k];
            }
        }
    }
    return ret;
}

#endif
