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

#include "Accuracy_Atm.h"
#include "Utils.h"

using namespace std;

Accuracy_Atm::Accuracy_Atm( int im, int jm, int km, double dthe, double dphi )
{
    this-> im = im;
    this-> jm = jm;
    this-> km = km;
    this-> dthe = dthe;
    this-> dphi = dphi;
}


Accuracy_Atm::Accuracy_Atm( int im, int jm, int km, double dr, double dthe, double dphi )
{
    this-> im = im;
    this-> jm = jm;
    this-> km = km;
    this-> dr = dr;
    this-> dthe = dthe;
    this-> dphi = dphi;
}


Accuracy_Atm::Accuracy_Atm ( int n, int nm, int Ma, int im, int jm, int km, double min, int j_res, int k_res )
{
    this-> n = n;
    this-> nm = nm;
    this-> Ma = Ma;
    this-> im = im;
    this-> jm = jm;
    this-> km = km;
    this-> min = min;
    this-> j_res = j_res;
    this-> k_res = k_res;
}


Accuracy_Atm::Accuracy_Atm ( int n, int nm, int Ma, int im, int jm, int km, double min, int i_res, int j_res, int k_res, double L_atm )
{
    this-> n = n;
    this-> nm = nm;
    this-> Ma = Ma;
    this-> im = im;
    this-> jm = jm;
    this-> km = km;
    this-> min = min;
    this-> i_res = i_res;
    this-> j_res = j_res;
    this-> k_res = k_res;
    this-> L_atm = L_atm;
}


Accuracy_Atm::~Accuracy_Atm () {}




double Accuracy_Atm::residuumQuery_2D ( Array_1D &rad, Array_1D &the, Array &v, Array &w )
{
    // value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
    min = residuum = 0.;

    for ( int j = 1; j < jm-1; j++ )
    {
        sinthe = sin( the.z[ j ] );
        costhe = cos( the.z[ j ] );
        rmsinthe = rad.z[ 0 ] * sinthe;

        for ( int k = 1; k < km-1; k++ )
        {
                dvdthe = ( v.x[ 0 ][ j+1 ][ k ] - v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );
                dwdphi = ( w.x[ 0 ][ j ][ k+1 ] - w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
                residuum = dvdthe / rad.z[ 0 ] + costhe / rmsinthe * v.x[ 0 ][ j ][ k ] + dwdphi / rmsinthe;
            if ( fabs ( residuum ) >= min )
            {
                min = residuum;
                j_res = j;
                k_res = k;
            }
        }
    }
    return 0;
}

double Accuracy_Atm::residuumQuery_3D ( Array_1D &rad, Array_1D &the, Array &u, Array &v, Array &w )
{
    // value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
    min = residuum = 0.;

    for ( int i = 1; i < im-1; i++ )
    {
        for ( int j = 1; j < jm-1; j++ )
        {
            sinthe = sin( the.z[ j ] );
            costhe = cos( the.z[ j ] );
            rmsinthe = rad.z[ i ] * sinthe;

            for ( int k = 1; k < km-1; k++ )
            {
                dudr = ( u.x[ i+1 ][ j ][ k ] - u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
                dvdthe = ( v.x[ i ][ j+1 ][ k ] - v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
                dwdphi = ( w.x[ i ][ j ][ k+1 ] - w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );

                residuum = dudr + 2. * u.x[ i ][ j ][ k ] / rad.z[ i ] + dvdthe / rad.z[ i ]
                            + costhe / rmsinthe * v.x[ i ][ j ][ k ] + dwdphi / rmsinthe;
                if ( fabs ( residuum ) >= min )
                {
                    min = residuum;
                    i_res = i;
                    j_res = j;
                    k_res = k;
                }
            }
        }
    }
    return 0;
}

double Accuracy_Atm::steadyQuery_2D ( Array &v, Array &vn, Array &w, Array &wn, Array &p_dyn, Array &p_dynn )
{
    // state of a steady solution ( min )
    min_v = max_v = 0.;
    min_w = max_w = 0.;
    min_p = max_p = 0.;

    for ( int j = 0; j < jm; j++ )
    {
        for ( int k = 0; k < km; k++ )
        {
            max_v = fabs ( v.x[ 0 ][ j ][ k ] - vn.x[ 0 ][ j ][ k ] );
            if ( max_v >= min_v )
            {
                min_v = max_v;
                j_v = j;
                k_v = k;
            }

            max_w = fabs ( w.x[ 0 ][ j ][ k ] - wn.x[ 0 ][ j ][ k ] );
            if ( max_w >= min_w )
            {
                min_w = max_w;
                j_w = j;
                k_w = k;
            }

            max_p = fabs ( p_dyn.x[ 0 ][ j ][ k ] - p_dynn.x[ 0 ][ j ][ k ] );
            if ( max_p >= min_p )
            {
                min_p = max_p;
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

    return 0;
}


void Accuracy_Atm::print(const string& name, double value, int j, int k) const{

    AtomUtils::HemisphereCoords coords = AtomUtils::convert_coords(k, j);

    cout << setiosflags ( ios::left ) << setw ( 36 ) << setfill ( '.' ) << name << " = " << resetiosflags ( ios::left ) << 
        setw ( 12 ) << fixed << setfill ( ' ' ) << value << setw ( 5 ) << coords.lat << setw ( 3 ) << coords.north_or_south 
        << setw ( 4 ) << coords.lon << setw ( 3 ) << coords.east_or_west << endl;
}




double Accuracy_Atm::steadyQuery_3D ( Array &u, Array &un, Array &v, Array &vn, Array &w, Array &wn, Array &t, Array &tn, 
    Array &c, Array &cn, Array &cloud, Array &cloudn, Array &ice, Array &icen, Array &co2, Array &co2n, Array &p_dyn, 
    Array &p_dynn )
{
    // state of a steady solution ( min )
    min_u = max_u = 0.;
    min_v = max_v = 0.;
    min_w = max_w = 0.;
    min_t = max_t = 0.;
    min_c = max_c = 0.;
    min_p = max_p = 0.;
    min_cloud = max_cloud = 0.;
    min_ice = max_ice = 0.;
    min_co2 = max_co2 = 0.;

    for ( int i = 0; i < im; i++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int k = 0; k < km; k++ )
            {
                max_u = fabs ( u.x[ i ][ j ][ k ] - un.x[ i ][ j ][ k ] );
                if ( max_u >= min_u )
                {
                    min_u = max_u;
                    i_u = i;
                    j_u = j;
                    k_u = k;
                }

                max_v = fabs ( v.x[ i ][ j ][ k ] - vn.x[ i ][ j ][ k ] );
                if ( max_v >= min_v )
                {
                    min_v = max_v;
                    i_v = i;
                    j_v = j;
                    k_v = k;
                }

                max_w = fabs ( w.x[ i ][ j ][ k ] - wn.x[ i ][ j ][ k ] );
                if ( max_w >= min_w )
                {
                    min_w = max_w;
                    i_w = i;
                    j_w = j;
                    k_w = k;
                }

                max_t = fabs ( t.x[ i ][ j ][ k ] - tn.x[ i ][ j ][ k ] );
                if ( max_t >= min_t )
                {
                    min_t = max_t;
                    i_t = i;
                    j_t = j;
                    k_t = k;
                }

                max_c = fabs ( c.x[ i ][ j ][ k ] - cn.x[ i ][ j ][ k ] );
                if ( max_c >= min_c )
                {
                    min_c = max_c;
                    i_c = i;
                    j_c = j;
                    k_c = k;
                }

                max_cloud = fabs ( cloud.x[ i ][ j ][ k ] - cloudn.x[ i ][ j ][ k ] );
                if ( max_cloud >= min_cloud )
                {
                    min_cloud = max_cloud;
                    i_cloud = i;
                    j_cloud = j;
                    k_cloud = k;
                }

                max_ice = fabs ( ice.x[ i ][ j ][ k ] - icen.x[ i ][ j ][ k ] );
                if ( max_ice >= min_ice )
                {
                    min_ice = max_ice;
                    i_ice = i;
                    j_ice = j;
                    k_ice = k;
                }

                max_co2 = fabs ( co2.x[ i ][ j ][ k ] - co2n.x[ i ][ j ][ k ] );
                if ( max_co2 >= min_co2 )
                {
                    min_co2= max_co2;
                    i_co2 = i;
                    j_co2 = j;
                    k_co2 = k;
                }

                max_p = fabs ( p_dyn.x[ i ][ j ][ k ] - p_dynn.x[ i ][ j ][ k ] );
                if ( max_p >= min_p )
                {
                    min_p = max_p;
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

    print(" residuum: continuity equation ", min, i_res * int ( L_atm ) / ( im - 1 ), j_res, k_res);

    print(" dp: pressure Poisson equation ", min_p, i_p * int ( L_atm ) / ( im - 1 ), j_p, k_p);

    print(" du: Navier Stokes equation ", min_u, i_u * int ( L_atm ) / ( im - 1 ), j_u, k_u);

    print(" dv: Navier Stokes equation ", min_v, i_v * int ( L_atm ) / ( im - 1 ), j_v, k_v);

    print(" dw: Navier Stokes equation ", min_w, i_w * int ( L_atm ) / ( im - 1 ), j_w, k_w);

    print(" dt: energy transport equation ", min_t, i_t * int ( L_atm ) / ( im - 1 ), j_t, k_t);

    print(" dc: vapour transport equation ", min_c, i_c * int ( L_atm ) / ( im - 1 ), j_c, k_c);

    print(" dcloud: cloud transport equation ", min_cloud, i_cloud * int ( L_atm ) / ( im - 1 ), j_cloud, k_cloud);

    print(" dice: ice transport equation ", min_ice, i_ice * int ( L_atm ) / ( im - 1 ), j_ice, k_ice);

    print(" dco: CO2 transport equation ", min_co2, i_co2 * int ( L_atm ) / ( im - 1 ), j_co2, k_co2);

    return 0;
}

void Accuracy_Atm::print(const string& name, double value, int i, int j, int k) const{

    AtomUtils::HemisphereCoords coords = AtomUtils::convert_coords(k, j);

    cout << setiosflags ( ios::left ) << setw ( 36 ) << setfill ( '.' ) << name << " = " << resetiosflags ( ios::left ) << 
        setw ( 12 ) << fixed << setfill ( ' ' ) << value << setw ( 5 ) << coords.lat << setw ( 3 ) << coords.north_or_south 
        << setw ( 4 ) << coords.lon << setw ( 3 ) << coords.east_or_west << setw ( 6 ) << i << setw ( 2 ) << "m" << endl;
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
