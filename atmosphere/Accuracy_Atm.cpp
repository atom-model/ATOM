/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to surveil the accuracy of the iterations
*/


#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>

#include "Array_1D.h"
#include "Accuracy_Atm.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

Accuracy_Atm::Accuracy_Atm( int im, int jm, int km, double dthe, double dphi ):
    im(im),
    jm(jm),
    km(km),
    i_res(0), 
    j_res(0), 
    k_res(0),
    dr(0),
    dthe(dthe),
    dphi(dphi),
    min(0.),
    is_3d_flag(false)
{}


Accuracy_Atm::Accuracy_Atm( int im, int jm, int km, double dr, double dthe, double dphi ):
    im(im),
    jm(jm),
    km(km),
    i_res(0), 
    j_res(0), 
    k_res(0),
    dr(dr),
    dthe(dthe),
    dphi(dphi),
    min(0.),
    is_3d_flag(true)
{}

Accuracy_Atm::~Accuracy_Atm () {}

std::tuple<double, int, int> 
Accuracy_Atm::residuumQuery_2D ( Array_1D &rad, Array_1D &the, Array &v, Array &w, Vector3D<> &residuum_2d )
{
    // value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
    for ( int j = 1; j < jm-1; j++ )
    {
        double sinthe = sin( the.z[ j ] );
        double costhe = cos( the.z[ j ] );
        double rmsinthe = rad.z[ 0 ] * sinthe;

        for ( int k = 1; k < km-1; k++ )
        {
            double dvdthe = ( v.x[ 0 ][ j+1 ][ k ] - v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );
            double dwdphi = ( w.x[ 0 ][ j ][ k+1 ] - w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
            double &residuum = residuum_2d(0,j,k) = 
                dvdthe / rad.z[ 0 ] + costhe / rmsinthe * v.x[ 0 ][ j ][ k ] + dwdphi / rmsinthe;
                
            if ( fabs ( residuum ) > min )
            {
                min = residuum;
                j_res = j;
                k_res = k;
            }
        }
    }
    return std::make_tuple(min, j_res, k_res);
}

std::tuple<double, int, int, int>
Accuracy_Atm::residuumQuery_3D ( Array_1D &rad, Array_1D &the, Array &u, Array &v, Array &w, 
                                        Vector3D<> &residuum_3d )
{
    assert(is_3d_flag);
    // value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
    for ( int i = 1; i < im-1; i++ )
    {
        for ( int j = 1; j < jm-1; j++ )
        {
            double sinthe = sin( the.z[ j ] );
            double costhe = cos( the.z[ j ] );
            double rmsinthe = rad.z[ i ] * sinthe;

            for ( int k = 1; k < km-1; k++ )
            {
                double dudr = ( u.x[ i+1 ][ j ][ k ] - u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
                double dvdthe = ( v.x[ i ][ j+1 ][ k ] - v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
                double dwdphi = ( w.x[ i ][ j ][ k+1 ] - w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );

                double &residuum = residuum_3d(i,j,k) = dudr + 2. * u.x[ i ][ j ][ k ] / rad.z[ i ] + dvdthe / rad.z[ i ]
                            + costhe / rmsinthe * v.x[ i ][ j ][ k ] + dwdphi / rmsinthe;
                //logger() << residuum_3d(i,j,k) << "  " << residuum << "  residuum" << std::endl;
                if ( fabs ( residuum ) > min )
                {
                    min = residuum;
                    i_res = i;
                    j_res = j;
                    k_res = k;
                }
            }
        }
    }
    return std::tuple<double, int, int, int>(min, i_res, j_res, k_res);
}

void Accuracy_Atm::steadyQuery_2D ( Array &v, Array &vn, Array &w, Array &wn, Array &p_dyn, Array &p_dynn )
{
    // state of a steady solution ( min )
    double min_v = 0., min_w = 0., min_p = 0., tmp = 0.;
    int j_v, k_v, j_w, k_w, j_p, k_p;
    j_v = k_v = j_w = k_w = j_p = k_p = 0;    

    for ( int j = 0; j < jm; j++ )
    {
        for ( int k = 0; k < km; k++ )
        {
            tmp = fabs ( v.x[ 0 ][ j ][ k ] - vn.x[ 0 ][ j ][ k ] );
            if ( tmp > min_v )
            {
                min_v = tmp;
                j_v = j;
                k_v = k;
            }

            tmp = fabs ( w.x[ 0 ][ j ][ k ] - wn.x[ 0 ][ j ][ k ] );
            if ( tmp > min_w )
            {
                min_w = tmp;
                j_w = j;
                k_w = k;
            }

            tmp = fabs ( p_dyn.x[ 0 ][ j ][ k ] - p_dynn.x[ 0 ][ j ][ k ] );
            if ( tmp > min_p )
            {
                min_p = tmp;
                j_p = j;
                k_p = k;
            }
        }
    }


    // statements on the convergence und iterational process
    cout.precision ( 6 );
    cout.setf ( ios::fixed );

    // printout of maximum and minimum absolute and relative errors of the computed values at their locations while iterating

    cout << endl << endl << 
        " 2D iterational process for the surface boundary conditions\n printout of maximum and minimum absolute and " <<
        "relative errors of the computed values at their locations: level, latitude, longitude"
         << endl << endl;

    
    print( " residuum: continuity equation ", min, j_res, k_res);

    print( " dp: pressure Poisson equation ", min_p, j_p, k_p);

    print( " dv: Navier Stokes equation ", min_v, j_v, k_v);

    print( " dw: Navier Stokes equation ", min_w, j_w, k_w);

    return;
}


void Accuracy_Atm::print(const string& name, double value, int j, int k) const{

    AtomUtils::HemisphereCoords coords = AtomUtils::convert_coords(k, j);

    cout << setiosflags ( ios::left ) << setw ( 36 ) << setfill ( '.' ) << name << " = " << resetiosflags ( ios::left ) << 
        setw ( 12 ) << fixed << setfill ( ' ' ) << value << setw ( 5 ) << int(coords.lat) << setw ( 3 ) << coords.north_or_south 
        << setw ( 4 ) << int(coords.lon) << setw ( 3 ) << coords.east_or_west << endl;
}




void Accuracy_Atm::steadyQuery_3D ( Array &u, Array &un, Array &v, Array &vn, Array &w, Array &wn, Array &t, Array &tn, 
    Array &c, Array &cn, Array &cloud, Array &cloudn, Array &ice, Array &icen, Array &co2, Array &co2n, Array &p_dyn, 
    Array &p_dynn, double L_atm, Array_1D &rad )
{
    // state of a steady solution ( min )
    double min_u = 0.; 
    double min_v = 0.;
    double min_w = 0.;
    double min_t = 0.;
    double min_c = 0.;
    double min_p = 0.;
    double min_cloud = 0.;
    double min_ice = 0.;
    double min_co2 = 0., tmp = 0.;

    int i_u = 0, j_u = 0, k_u = 0; 
    int i_v = 0, j_v = 0, k_v = 0;
    int i_w = 0, j_w = 0, k_w = 0;
    int i_t = 0, j_t = 0, k_t = 0;
    int i_c = 0, j_c = 0, k_c = 0;
    int i_cloud = 0, j_cloud = 0, k_cloud = 0;
    int i_ice = 0, j_ice = 0, k_ice = 0;
    int i_co2 = 0, j_co2 = 0, k_co2 = 0; 
    int i_p = 0, j_p = 0, k_p = 0;

    for ( int i = 0; i < im; i++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int k = 0; k < km; k++ )
            {
                tmp = fabs ( u.x[ i ][ j ][ k ] - un.x[ i ][ j ][ k ] );
                if ( tmp > min_u )
                {
                    min_u = tmp;
                    i_u = i;
                    j_u = j;
                    k_u = k;
                }

                tmp = fabs ( v.x[ i ][ j ][ k ] - vn.x[ i ][ j ][ k ] );
                if ( tmp > min_v )
                {
                    min_v = tmp;
                    i_v = i;
                    j_v = j;
                    k_v = k;
                }

                tmp = fabs ( w.x[ i ][ j ][ k ] - wn.x[ i ][ j ][ k ] );
                if ( tmp > min_w )
                {
                    min_w = tmp;
                    i_w = i;
                    j_w = j;
                    k_w = k;
                }

                tmp = fabs ( t.x[ i ][ j ][ k ] - tn.x[ i ][ j ][ k ] );
                if ( tmp > min_t )
                {
                    min_t = tmp;
                    i_t = i;
                    j_t = j;
                    k_t = k;
                }

                tmp = fabs ( c.x[ i ][ j ][ k ] - cn.x[ i ][ j ][ k ] );
                if ( tmp > min_c )
                {
                    min_c = tmp;
                    i_c = i;
                    j_c = j;
                    k_c = k;
                }

                tmp = fabs ( cloud.x[ i ][ j ][ k ] - cloudn.x[ i ][ j ][ k ] );
                if ( tmp > min_cloud )
                {
                    min_cloud = tmp;
                    i_cloud = i;
                    j_cloud = j;
                    k_cloud = k;
                }

                tmp = fabs ( ice.x[ i ][ j ][ k ] - icen.x[ i ][ j ][ k ] );
                if ( tmp > min_ice )
                {
                    min_ice = tmp;
                    i_ice = i;
                    j_ice = j;
                    k_ice = k;
                }

                tmp = fabs ( co2.x[ i ][ j ][ k ] - co2n.x[ i ][ j ][ k ] );
                if ( tmp > min_co2 )
                {
                    min_co2 = tmp;
                    i_co2 = i;
                    j_co2 = j;
                    k_co2 = k;
                }

                tmp = fabs ( p_dyn.x[ i ][ j ][ k ] - p_dynn.x[ i ][ j ][ k ] );
                if ( tmp > min_p )
                {
                    min_p = tmp;
                    i_p = i;
                    j_p = j;
                    k_p = k;
                }
            }
        }
    }

    // statements on the convergence und iterational process
    cout.precision ( 6 );
    cout.setf ( ios::fixed );

    // printout of maximum and minimum absolute and relative errors of the computed values at their locations while iterating

    cout << endl << endl << " 3D iterational process for the surface boundary conditions\n printout of maximum and minimum " <<
        "absolute and relative errors of the computed values at their locations: level, latitude, longitude" << endl;
    
    print(" residuum: continuity equation ", min, i_res * int ( L_atm ) / ( im - 1 ), j_res, k_res, L_atm, rad );

    print(" dp: pressure Poisson equation ", min_p, i_p * int ( L_atm ) / ( im - 1 ), j_p, k_p, L_atm, rad );

    print(" du: Navier Stokes equation ", min_u, i_u * int ( L_atm ) / ( im - 1 ), j_u, k_u, L_atm, rad );

    print(" dv: Navier Stokes equation ", min_v, i_v * int ( L_atm ) / ( im - 1 ), j_v, k_v, L_atm, rad );

    print(" dw: Navier Stokes equation ", min_w, i_w * int ( L_atm ) / ( im - 1 ), j_w, k_w, L_atm, rad );

    print(" dt: energy transport equation ", min_t, i_t * int ( L_atm ) / ( im - 1 ), j_t, k_t, L_atm, rad );

    print(" dc: vapour transport equation ", min_c, i_c * int ( L_atm ) / ( im - 1 ), j_c, k_c, L_atm, rad );

    print(" dcloud: cloud transport equation ", min_cloud, i_cloud * int ( L_atm ) / ( im - 1 ), j_cloud, k_cloud, L_atm, rad );

    print(" dice: ice transport equation ", min_ice, i_ice * int ( L_atm ) / ( im - 1 ), j_ice, k_ice, L_atm, rad );

    print(" dco: CO2 transport equation ", min_co2, i_co2 * int ( L_atm ) / ( im - 1 ), j_co2, k_co2, L_atm, rad );

    return;
}

void Accuracy_Atm::print(const string& name, double value, int i, int j, int k, double L_atm, Array_1D &rad ) const{

    AtomUtils::HemisphereCoords coords = AtomUtils::convert_coords(k, j);

    double zeta = 3.715;

    i = ( int )( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching

    cout << setiosflags ( ios::left ) << setw ( 36 ) << setfill ( '.' ) << name << " = " << resetiosflags ( ios::left ) << 
        setw ( 12 ) << fixed << setfill ( ' ' ) << value << setw ( 5 ) << int(coords.lat) << setw ( 3 ) << coords.north_or_south 
        << setw ( 4 ) << int(coords.lon) << setw ( 3 ) << coords.east_or_west << setw ( 6 ) << i << setw ( 2 ) << "m" << endl;
}


double Accuracy_Atm::out_min (  ) const
{
    return min;
}

int Accuracy_Atm::out_i_res (  ) const
{
    return i_res;
}

int Accuracy_Atm::out_j_res (  ) const
{
    return j_res;
}

int Accuracy_Atm::out_k_res (  ) const
{
    return k_res;
}
