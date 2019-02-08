/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to produce resulting data on mean sea level
*/

#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include "Results_Hyd.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


Results_Hyd::Results_Hyd ( int im, int jm, int km ){
    this-> im = im;
    this-> jm = jm;
    this-> km = km;

    c43 = 4./3.;
    c13 = 1./3.;

// array "aux_grad_v" for for the computation of Ekman pumping
    aux_grad_v = 0L;

    aux_grad_v = new double[ im ];

    for ( int l = 0; l < im; l++ ){
        aux_grad_v[ l ] = 0.;
    }

// array "aux_grad_w" for for the computation of Ekman pumping
    aux_grad_w = 0L;

    aux_grad_w = new double[ im ];

    for ( int l = 0; l < im; l++ ){
        aux_grad_w[ l ] = 0.;
    }


// auxiliar 2D Array "aux_v" for the computation of Ekman pumping
    aux_v = 0L;

    aux_v = new double*[ jm ];

    for ( int l = 0; l < jm; l++ ){
        aux_v[ l ] = new double[ km ];
    }

// default values
    for ( int l = 0; l < jm; l++ ){
        for ( int n = 0; n < km; n++ )
        {
            aux_v[ l ][ n ] = 0.;
        }
    }


// auxiliar 2D Array "aux_w" for the computation of Ekman pumping
    aux_w = 0L;

    aux_w = new double*[ jm ];

    for ( int l = 0; l < jm; l++ ){
        aux_w[ l ] = new double[ km ];
    }

// default values
    for ( int l = 0; l < jm; l++ ){
        for ( int n = 0; n < km; n++ ){
            aux_w[ l ][ n ] = 0.;
        }
    }
}


Results_Hyd::~Results_Hyd (){
    for ( int j = 0; j < jm; j++ ){
        delete [  ] aux_v[ j ];
    }

    delete [  ] aux_v;

    for ( int j = 0; j < jm; j++ ){
        delete [  ] aux_w[ j ];
    }

    delete [  ] aux_w;
}




void Results_Hyd::run_data ( int i_beg, double dr, double dthe, double L_hyd, double u_0,
                            double c_0, Array_1D &rad, Array_1D &the, Array &h, Array &u, Array &v,
                            Array &w, Array &c, Array &Salt_Balance, Array &Salt_Finger,
                            Array &Salt_Diffusion, Array &Buoyancy_Force_3D, Array_2D &Upwelling,
                            Array_2D &Downwelling, Array_2D &SaltFinger, Array_2D &SaltDiffusion,
                            Array_2D &BuoyancyForce_2D, Array_2D &Salt_total,
                            Array_2D &BottomWater ){
    // total upwelling as sum on normal velocity component values in a virtual vertical column

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            Upwelling.y[ j ][ k ] = 0.; // Upwelling
            Downwelling.y[ j ][ k ] = 0.; // Downwelling
            aux_v[ j ][ k ] = 0.;  // auxiliar field for Upwelling
            aux_w[ j ][ k ] = 0.;  // auxiliar field for Downwelling
            BottomWater.y[ j ][ k ] = 0.;    // Bottom water
            SaltFinger.y[ j ][ k ] = 0.;     // SaltFinger
            SaltDiffusion.y[ j ][ k ] = 0.;    // SaltFinger
            BuoyancyForce_2D.y[ j ][ k ] = 0.;  // Saltdiffusion
            Salt_total.y[ j ][ k ] = 0.;   // total Salt
        }
    }

    double i_Ekman_layer = 500.;  // assumed Ekman-layer depth of 500m
    double coeff = i_Ekman_layer / L_hyd;
    double rmsinthe = 0.;

    int i_Ekman = ( im - 1 ) * ( 1. - coeff );
    int j_half = ( jm -1 ) / 2;
    int i_max = im - 1;
    int i_diff = i_max - i_Ekman;

    //Ekman pumping, upwelling, downwelling
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = i_Ekman; i < im-1; i++ ){
                if ( is_water( h, i, j, k) ){
                    aux_grad_v[ i ] = v.x[ i ][ j ][ k ];
                    aux_grad_w[ i ] = w.x[ i ][ j ][ k ];
                }else{
                    aux_v[ j ][ k ] = 0.;
                    aux_w[ j ][ k ] = 0.;

                    aux_grad_v[ i ] = 0.;
                    aux_grad_w[ i ] = 0.;
                }
            }
            if ( i_diff % 2 == 0 ){
                aux_v[ j ][ k ] = simpson ( i_diff, dr, aux_grad_v );
                aux_w[ j ][ k ] = simpson ( i_diff, dr, aux_grad_w );
//                aux_v[ j ][ k ] = trapezoidal ( i_diff, dr, aux_grad_v );
//                aux_w[ j ][ k ] = trapezoidal ( i_diff, dr, aux_grad_w );
//                aux_v[ j ][ k ] = rectangular ( i_diff, dr, aux_grad_v );
//                aux_w[ j ][ k ] = rectangular ( i_diff, dr, aux_grad_w );
            }
            else cout << "       i_diff = i_max - i_Ekman    must be an even number to use the Simpson integration method" << endl;
        }
    }

    for ( int k = 1; k < km-1; k++ ){
        for ( int j = 1; j < j_half-1; j++ ){
            rmsinthe = rad.z[ im-1 ] * sin( the.z[ j ] );

            BottomWater.y[ j ][ k ] = ( aux_v[ j + 1 ][ k ] - aux_v[ j - 1 ][ k ] ) /
                ( 2. * rad.z[ im-1 ] * dthe ) + ( aux_w[ j ][ k + 1 ] - aux_w[ j ][ k + 1 ] )
                / ( 2. * rmsinthe * dphi );
            BottomWater.y[ j ][ k ] = BottomWater.y[ j ][ k ] * u_0;
        }

        for ( int j = j_half+2; j < jm-1; j++ ){
            rmsinthe = rad.z[ im-1 ] * sin( the.z[ j ] );

            BottomWater.y[ j ][ k ] = ( aux_v[ j + 1 ][ k ] - aux_v[ j - 1 ][ k ] ) /
                ( 2. * rad.z[ im-1 ] * dthe )
                                                    + ( aux_w[ j ][ k + 1 ] - aux_w[ j ][ k + 1 ] ) /
                                                    ( 2. * rmsinthe * dphi );
            BottomWater.y[ j ][ k ] = BottomWater.y[ j ][ k ] * u_0;
        }

        for ( int j = 1; j < jm-1; j++ ){
            if ( BottomWater.y[ j ][ k ] >= 0. )  Upwelling.y[ j ][ k ] = BottomWater.y[ j ][ k ];
            else  Upwelling.y[ j ][ k ] = 0.;

            if ( BottomWater.y[ j ][ k ] < 0. )  Downwelling.y[ j ][ k ] = BottomWater.y[ j ][ k ];
            else  Downwelling.y[ j ][ k ] = 0.;
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            Downwelling.y[ j ][ k ] = fabs ( Downwelling.y[ j ][ k ] );
        }
    }


    for ( int k = 0; k < km; k++ ){
        for ( int i = 0; i < im; i++ ){
            Upwelling.y[ 0 ][ k ] = c43 * Upwelling.y[ 1 ][ k ] -
                c13 * Upwelling.y[ 2 ][ k ];
            Upwelling.y[ jm-1 ][ k ] = c43 * Upwelling.y[ jm-2 ][ k ] -
                c13 * Upwelling.y[ jm-3 ][ k ];

            Downwelling.y[ 0 ][ k ] = c43 * Downwelling.y[ 1 ][ k ] -
                c13 * Downwelling.y[ 2 ][ k ];
            Downwelling.y[ jm-1 ][ k ] = c43 * Downwelling.y[ jm-2 ][ k ] -
                c13 * Downwelling.y[ jm-3 ][ k ];
        }
    }

    for ( int i = 0; i < im; i++ ){
        for ( int j = 0; j < jm; j++ ){
            Upwelling.y[ j ][ 0 ] = c43 * Upwelling.y[ j ][ 1 ] - c13 * Upwelling.y[ j ][ 2 ];
            Upwelling.y[ j ][ km-1 ] = c43 * Upwelling.y[ j ][ km-2 ] -
                c13 * Upwelling.y[ j ][ km-3 ];
            Upwelling.y[ j ][ 0 ] = Upwelling.y[ j ][ km-1 ] = ( Upwelling.y[ j ][ 0 ] +
                Upwelling.y[ j ][ km-1 ] ) / 2.;

            Downwelling.y[ j ][ 0 ] = c43 * Downwelling.y[ j ][ 1 ] - c13 * Downwelling.y[ j ][ 2 ];
            Downwelling.y[ j ][ km-1 ] = c43 * Downwelling.y[ j ][ km-2 ] -
                c13 * Downwelling.y[ j ][ km-3 ];
            Downwelling.y[ j ][ 0 ] = Downwelling.y[ j ][ km-1 ] = ( Downwelling.y[ j ][ 0 ] +
                Downwelling.y[ j ][ km-1 ] ) / 2.;
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = im/2; i < im; i++ ){
                if ( is_water( h, i, j, k) ){
                    SaltFinger.y[ j ][ k ] += Salt_Finger.x[ i ][ j ][ k ];
                    SaltDiffusion.y[ j ][ k ] += Salt_Diffusion.x[ i ][ j ][ k ];
                    BuoyancyForce_2D.y[ j ][ k ] += Buoyancy_Force_3D.x[ i ][ j ][ k ];
                    Salt_total.y[ j ][ k ] += c.x[ i ][ j ][ k ] * c_0;
                }
            }
        }
    }

    // boundaries of buoyancy force
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            Buoyancy_Force_3D.x[ 0 ][ j ][ k ] = c43 * Buoyancy_Force_3D.x[ 1 ][ j ][ k ] -
                c13 * Buoyancy_Force_3D.x[ 2 ][ j ][ k ];
            Buoyancy_Force_3D.x[ im-1 ][ j ][ k ] = c43 * Buoyancy_Force_3D.x[ im-2 ][ j ][ k ] -
                c13 * Buoyancy_Force_3D.x[ im-3 ][ j ][ k ];

            Salt_Finger.x[ 0 ][ j ][ k ] = c43 * Salt_Finger.x[ 1 ][ j ][ k ] -
                c13 * Salt_Finger.x[ 2 ][ j ][ k ];
            Salt_Finger.x[ im-1 ][ j ][ k ] = c43 * Salt_Finger.x[ im-2 ][ j ][ k ] -
                c13 * Salt_Finger.x[ im-3 ][ j ][ k ];

            Salt_Diffusion.x[ 0 ][ j ][ k ] = c43 * Salt_Diffusion.x[ 1 ][ j ][ k ] -
                c13 * Salt_Diffusion.x[ 2 ][ j ][ k ];
            Salt_Diffusion.x[ im-1 ][ j ][ k ] = c43 * Salt_Diffusion.x[ im-2 ][ j ][ k ] -
                c13 * Salt_Diffusion.x[ im-3 ][ j ][ k ];

            Salt_Balance.x[ 0 ][ j ][ k ] = c43 * Salt_Balance.x[ 1 ][ j ][ k ] -
                c13 * Salt_Balance.x[ 2 ][ j ][ k ];
            Salt_Balance.x[ im-1 ][ j ][ k ] = c43 * Salt_Balance.x[ im-2 ][ j ][ k ] -
                c13 * Salt_Balance.x[ im-3 ][ j ][ k ];
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int i = 0; i < im; i++ ){
            Buoyancy_Force_3D.x[ i ][ 0 ][ k ] = c43 * Buoyancy_Force_3D.x[ i ][ 1 ][ k ] -
                c13 * Buoyancy_Force_3D.x[ i ][ 2 ][ k ];
            Buoyancy_Force_3D.x[ i ][ jm-1 ][ k ] = c43 * Buoyancy_Force_3D.x[ i ][ jm-2 ][ k ] -
                c13 * Buoyancy_Force_3D.x[ i ][ jm-3 ][ k ];

            Salt_Finger.x[ i ][ 0 ][ k ] = c43 * Salt_Finger.x[ i ][ 1 ][ k ] -
                c13 * Salt_Finger.x[ i ][ 2 ][ k ];
            Salt_Finger.x[ i ][ jm-1 ][ k ] = c43 * Salt_Finger.x[ i ][ jm-2 ][ k ] -
                c13 * Salt_Finger.x[ i ][ jm-3 ][ k ];

            Salt_Diffusion.x[ i ][ 0 ][ k ] = c43 * Salt_Diffusion.x[ i ][ 1 ][ k ] -
                 c13 * Salt_Diffusion.x[ i ][ 2 ][ k ];
            Salt_Diffusion.x[ i ][ jm-1 ][ k ] = c43 * Salt_Diffusion.x[ i ][ jm-2 ][ k ] -
                c13 * Salt_Diffusion.x[ i ][ jm-3 ][ k ];

            Salt_Balance.x[ i ][ 0 ][ k ] = c43 * Salt_Balance.x[ i ][ 1 ][ k ] -
                c13 * Salt_Balance.x[ i ][ 2 ][ k ];
            Salt_Balance.x[ i ][ jm-1 ][ k ] = c43 * Salt_Balance.x[ i ][ jm-2 ][ k ] -
                c13 * Salt_Balance.x[ i ][ jm-3 ][ k ];
        }
    }

    for ( int i = 0; i < im; i++ ){
        for ( int j = 0; j < jm; j++ ){
            Buoyancy_Force_3D.x[ i ][ j ][ 0 ] = c43 * Buoyancy_Force_3D.x[ i ][ j ][ 1 ] -
                c13 * Buoyancy_Force_3D.x[ i ][ j ][ 2 ];
            Buoyancy_Force_3D.x[ i ][ j ][ km-1 ] = c43 * Buoyancy_Force_3D.x[ i ][ j ][ km-2 ] -
                c13 * Buoyancy_Force_3D.x[ i ][ j ][ km-3 ];
            Buoyancy_Force_3D.x[ i ][ j ][ 0 ] = Buoyancy_Force_3D.x[ i ][ j ][ km-1 ] =
                ( Buoyancy_Force_3D.x[ i ][ j ][ 0 ] + Buoyancy_Force_3D.x[ i ][ j ][ km-1 ] ) / 2.;

            Salt_Finger.x[ i ][ j ][ 0 ] = c43 * Salt_Finger.x[ i ][ j ][ 1 ] -
                c13 * Salt_Finger.x[ i ][ j ][ 2 ];
            Salt_Finger.x[ i ][ j ][ km-1 ] = c43 * Salt_Finger.x[ i ][ j ][ km-2 ] -
                c13 * Salt_Finger.x[ i ][ j ][ km-3 ];
            Salt_Finger.x[ i ][ j ][ 0 ] = Salt_Finger.x[ i ][ j ][ km-1 ] = ( Salt_Finger.x[ i ][ j ][ 0 ] +
                 Salt_Finger.x[ i ][ j ][ km-1 ] ) / 2.;

            Salt_Diffusion.x[ i ][ j ][ 0 ] = c43 * Salt_Diffusion.x[ i ][ j ][ 1 ] -
                c13 * Salt_Diffusion.x[ i ][ j ][ 2 ];
            Salt_Diffusion.x[ i ][ j ][ km-1 ] = c43 * Salt_Diffusion.x[ i ][ j ][ km-2 ] -
                c13 * Salt_Diffusion.x[ i ][ j ][ km-3 ];
            Salt_Diffusion.x[ i ][ j ][ 0 ] = Salt_Diffusion.x[ i ][ j ][ km-1 ] =
                ( Salt_Diffusion.x[ i ][ j ][ 0 ] + Salt_Diffusion.x[ i ][ j ][ km-1 ] ) / 2.;

            Salt_Balance.x[ i ][ j ][ 0 ] = c43 * Salt_Balance.x[ i ][ j ][ 1 ] -
                c13 * Salt_Balance.x[ i ][ j ][ 2 ];
            Salt_Balance.x[ i ][ j ][ km-1 ] = c43 * Salt_Balance.x[ i ][ j ][ km-2 ] -
                c13 * Salt_Balance.x[ i ][ j ][ km-3 ];
            Salt_Balance.x[ i ][ j ][ 0 ] = Salt_Balance.x[ i ][ j ][ km-1 ] =
                ( Salt_Balance.x[ i ][ j ][ 0 ] + Salt_Balance.x[ i ][ j ][ km-1 ] ) / 2.;
        }
    }

    cout.precision ( 4 );

// printout of surface data at one predefinded location
    level = "m";
    deg_north = "°N";
    deg_south = "°S";
    deg_west = "°W";
    deg_east = "°E";

    name_Value_1 = " downwelling ";
    name_Value_2 = " upwelling ";
    name_Value_3 = " bottom water ";
    name_Value_4 = " salt finger ";
    name_Value_5 = " salt diffusion ";
    name_Value_6 = " salt total ";

    name_unit_ms = " m/s";
    name_unit_psu = " psu";

    heading = " printout of surface data at predefinded locations: level, latitude, longitude";

    i_loc_level = 0;                                                                        // only at sea level MSL, constant
    j_loc = 90;
    k_loc = 180;


    if ( j_loc <= 90 ){
        j_loc_deg = 90 - j_loc;
        deg_lat = deg_north;
    }

    if ( j_loc > 90 ){
        j_loc_deg = j_loc - 90;
        deg_lat = deg_south;
    }

    if ( k_loc <= 180 ){
        k_loc_deg = 180 - k_loc;
        deg_lon = deg_west;
    }

    if ( k_loc > 180 ){
        k_loc_deg = k_loc - 180;
        deg_lon = deg_east;
    }

    Value_1 = Downwelling.y[ j_loc ][ k_loc ];
    Value_2 = Upwelling.y[ j_loc ][ k_loc ];
    Value_3 = BottomWater.y[ j_loc ][ k_loc ];
    Value_4 = SaltFinger.y[ j_loc ][ k_loc ];
    Value_5 = BuoyancyForce_2D.y[ j_loc ][ k_loc ];
    Value_6 = Salt_total.y[ j_loc ][ k_loc ];

    cout << endl << endl << heading << endl << endl;

    cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 )
        << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg
        << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 20 )
        << setfill ( '.' ) << name_Value_1 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_1 << setw ( 6 ) << name_unit_ms
        << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_2
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_2 << setw ( 6 ) << name_unit_ms << "   " << setiosflags ( ios::left )
        << setw ( 20 ) << setfill ( '.' ) << name_Value_3 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_3 << setw ( 6 ) << name_unit_ms
        << endl << "                       " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' )
        << name_Value_4 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed
        << setfill ( ' ' ) << Value_4 << setw ( 6 ) << name_unit_psu << "   "
        << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_5
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_5 << setw ( 6 ) << name_unit_psu << "   " << setiosflags ( ios::left )
        << setw ( 20 ) << setfill ( '.' ) << name_Value_6 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_6 << setw ( 6 ) << name_unit_psu
        << endl << endl << endl;
}




void Results_Hyd::land_oceanFraction ( Array &h ){
// calculation of the ratio ocean to land, also addition and substraction of CO2 of land, ocean and vegetation
    h_point_max =  ( jm - 1 ) * ( km - 1 );

    h_land = 0;

    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_land( h, im-1, j, k) )        h_land++;
        }
    }

    h_ocean = h_point_max - h_land;

    ozean_land = ( double ) h_ocean / ( double ) h_land;

    cout.precision ( 3 );

    cout << endl;
    cout << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' )
        << "      total number of points at constant hight " << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << h_point_max << endl << setiosflags ( ios::left )
        << setw ( 50 ) << setfill ( '.' ) << "      number of points on the ocean surface " << " = "
        << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_ocean << endl
        << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the land surface "
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_land
        << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      ocean/land ratio "
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << ozean_land << endl << endl;
    cout << endl;
}


