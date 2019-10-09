/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the salt concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to read and prepare the bathymetric and topografic data
*/


#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 

#include "BC_Bath_Hyd.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

BC_Bathymetry_Hydrosphere::BC_Bathymetry_Hydrosphere ( int im, int jm, int km ):
    im(im),
    jm(jm),
    km(km)
{ }

BC_Bathymetry_Hydrosphere::~BC_Bathymetry_Hydrosphere() {}

void BC_Bathymetry_Hydrosphere::BC_SeaGround(const string &bathymetry_file,
                                                          double L_hyd, Array &h, Array_2D &Bathymetry){
    // default adjustment, h must be 0 everywhere
    h.initArray(im, jm, km, 0.);

    ifstream ifile(bathymetry_file);

    if (!ifile.is_open()){
        cerr << "ERROR: could not open bathymetry file " << bathymetry_file << "\n";
        abort();
    }

    double lon, lat, depth;
    int j, k;
    for (j = 0; j < jm && !ifile.eof(); j++) {
        for (k = 0; k < km && !ifile.eof(); k++) {
            depth = 999; // in case the height is NaN
            ifile >> lon >> lat >> depth;

            if ( depth > 0. )
            {
                depth = 0;
            }
            
            int i_boden = (im-1) + int(floor( depth / L_hyd * ( im - 1 ) ) );
            Bathymetry.y[j][k] = -depth;
            for ( int i = 0; i <= i_boden; i++ ){
                h.x[i][j][k] = 1.;
            }

            if(ifile.fail())
            {
                ifile.clear();
                std::string tmp;
                std::getline(ifile, tmp);
                logger() << "bad data in topography at: " << lon << " " << lat << " " << tmp << std::endl;
            }
            //logger() << lon << " " << lat << " " << h.x[0][j][k] << std::endl;            
        }
    }

    if(j != jm || k != km ){
        std::cerr << "wrong topography file size! aborting..."<<std::endl;
        abort();
    }

    // rewrite bathymetric data from -180° - 0° - +180° to 0°- 360°
    for ( int j = 0; j < jm; j++ )
    {
        move_data(Bathymetry.y[j], km);
        for ( int i = 0; i < im; i++ )
        {
            move_data(h.x[i][j], km);
        }
    }

    // reduction and smoothing of peaks and needles in the bathymetry
    double h_center = 0.;

    for ( int k = 2; k < km-2; k++ ){
        for ( int j = 2; j < jm-2; j++ ){
            for ( int i = 2; i < im-2; i++ ){
                if ( ( is_land( h, i, j, k) ) && ( ( is_water( h, i, j+1, k) ) && ( is_water( h, i, j-1, k) ) )
                 && ( ( is_water( h, i, j, k+1) ) && ( is_water( h, i, j, k-1) ) ) ){
                     h_center = h.x[i][j][k];
                     h.x[i][j][k] = 0.;
                 }

                if ( ( h_center == 1. ) && ( ( is_water( h, i, j-2, k) ) && ( is_water( h, i, j, k+2) ) ) ){
                    h.x[i][j - 1][k + 1] = 0.;
                }
                if ( ( h_center == 1. ) && ( ( is_water( h, i, j-2, k) ) && ( is_water( h, i, j, k-2) ) ) ){ 
                    h.x[i][j - 1][k - 1] = 0.;
                }
                if ( ( h_center == 1. ) && ( ( is_water( h, i, j+2, k) ) && ( is_water( h, i, j, k-2) ) ) ){
                    h.x[i][j + 1][k - 1] = 0.;
                }
                if ( ( h_center == 1. ) && ( ( is_water( h, i, j+2, k) ) && ( is_water( h, i, j, k+2) ) ) ){ 
                    h.x[i][j + 1][k + 1] = 0.;
                }
            }
        }
    }
}





void BC_Bathymetry_Hydrosphere::BC_SolidGround ( Array &h, Array &t, Array &u, Array &v,
                                                    Array &w, Array &p_dyn, Array &c, Array &tn,
                                                    Array &un, Array &vn, Array &wn, Array &p_dynn, Array &cn){
    // boundary conditions for the total solid ground
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                if(is_land(h,  i,  j,  k )){
                    p_dyn.x[i][j][k] = p_dynn.x[i][j][k];
                    t.x[i][j][k] = tn.x[i][j][k];
                    c.x[i][j][k] = cn.x[i][j][k];
                    u.x[i][j][k] = un.x[i][j][k] = 0.;
                    v.x[i][j][k] = vn.x[i][j][k] = 0.;
                    w.x[i][j][k] = wn.x[i][j][k] = 0.;
                }
            }
        }
    }
}
