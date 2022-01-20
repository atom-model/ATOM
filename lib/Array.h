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
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <numeric>
#include <tuple>

using namespace std;

class Array
{
public:
    int im, jm, km;

public:
    double ***x;
    
    Array(int idim, int jdim, int kdim, double val);
    ~Array ( );

    Array(): im(0), jm(0), km(0), x(NULL)
    {}

    //copy constructor
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

    void printArray( string dir, int im, int jm, int km );
    void initArray( int im, int jm, int km, double value);

    int get_im() const {return im;}
    int get_jm() const {return jm;}
    int get_km() const {return km;}

    friend Array operator* (double coeff, const Array &a);

    //overload 
    void operator=(const Array &a){
        if(!x){
            initArray(a.im, a.jm, a.km, 0);
        }
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

    //overload
    Array operator+(double val){
        Array ret(im, jm, km, 0.);
        for(int i=0; i<im; i++){
            for(int j=0; j<jm; j++){
                for(int k=0; k<km; k++){
                    ret.x[i][j][k]=this->x[i][j][k]+val;
                }
            }
        }
        return ret;
    }

    //overload
    Array operator-(double val){
        return operator+(-val);
    }

    //overload
    Array operator*(double val){
        return val*(*this);
    }

    void save(const string& fn){
        for(int i=0; i<im; i++){
            save(fn, i);
        }
    }

    void save(const string& fn, int i){
        assert(i<im); 
        assert(i>=0);        
        std::ofstream os(fn + "_layer_" + to_string(i) + ".bin", std::ios::binary | std::ios::out);
        for(int j=jm-1; j>=0; j--){
            os.write(reinterpret_cast<const char*>(x[i][j]+(km/2)), std::streamsize((km/2)*sizeof(double)));
            os.write(reinterpret_cast<const char*>(x[i][j]), std::streamsize((km/2+1)*sizeof(double)));
        }
        os.close();
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
                    if(std::isnan(x[ i ][ j ][ k ])){
                        if ( std::isnan( x[ i ][ j ][ k ] ) )
                            cout << endl << "  ************* NaN detected ************  "
                            << endl << "   i = " << i << "   j = " << j << "   k = " << k << endl;
                        return true;
                    }
                }
            }
        }
        return false;
    }

    void inspect(const std::string& prefix="") const;
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

template <typename T = double>
class Vector3D {
public:
    Vector3D(size_t i, size_t j, size_t k, const T& t = T()) :
        m_i(i), m_j(j), m_k(k), m_data((!i?1:i)*j*k, t)
    {}

    //https://en.cppreference.com/w/cpp/language/operators
    //To provide multidimensional array access semantics, e.g. to implement a 3D array access a[i][j][k] = x;, 
    //operator[] has to return a reference to a 2D plane, which has to have its own operator[] which returns 
    //a reference to a 1D row, which has to have operator[] which returns a reference to the element. 
    //To avoid this complexity, some libraries opt for overloading operator() instead, 
    //so that 3D access expressions have the Fortran-like syntax a(i, j, k) = x;
    T& operator()(size_t i, size_t j, size_t k) {
        return m_data[i*m_j*m_k + j*m_k + k];
    }

    const T& operator()(size_t i, size_t j, size_t k) const {
        return m_data[i*m_j*m_k + j*m_k + k];
    }

    T mean() const{
        return std::accumulate(m_data.begin(), m_data.end(), 0.0) / m_data.size();
    }

    std::tuple<T, int, int, int>
    max() const
    {
        std::tuple<T, unsigned, unsigned, unsigned> ret(m_data[0], 0, 0, 0);
        for(std::size_t i=0; i<m_i; i++){
            for(std::size_t j=0; j<m_j; j++){
                for(std::size_t k=0; k<m_k; k++){
                    if(std::get<0>(ret) < m_data[i*m_j*m_k + j*m_k + k])
                    {
                        std::get<0>(ret) = m_data[i*m_j*m_k + j*m_k + k];
                        std::get<1>(ret) = i;
                        std::get<2>(ret) = j;
                        std::get<3>(ret) = k;
                    }
                }
            }
        }
        return ret;
    }

    std::tuple<T, int, int, int>
    min() const
    {
        std::tuple<T, unsigned, unsigned, unsigned> ret(m_data[0], 0, 0, 0);
        for(std::size_t i=0; i<m_i; i++){
            for(std::size_t j=0; j<m_j; j++){
                for(std::size_t k=0; k<m_k; k++){
                    if(std::get<0>(ret) > m_data[i*m_j*m_k + j*m_k + k])
                    {
                        std::get<0>(ret) = m_data[i*m_j*m_k + j*m_k + k];
                        std::get<1>(ret) = i;
                        std::get<2>(ret) = j;
                        std::get<3>(ret) = k;
                    }
                }
            }
        }
        return ret;
    }

private:
    size_t m_i, m_j, m_k;
    std::vector<T> m_data;
};

#endif
