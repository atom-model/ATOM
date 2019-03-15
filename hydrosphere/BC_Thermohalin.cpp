/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to prepare the boundary and initial conditions for diverse variables
*/


#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <Utils.h>

#include "BC_Thermohalin.h"
#include "Array.h"

using namespace std;
using namespace AtomUtils;

BC_Thermohalin::BC_Thermohalin ( int im, int jm, int km, int i_beg, int i_max,
                            int Ma, int Ma_max, int Ma_max_half, double dr, double g,
                            double r_0_water, double ua, double va, double wa, double ta,
                            double ca, double pa, double u_0, double p_0, double t_0,
                            double c_0, double cp_w, double L_hyd, double t_average,
                            double t_cretaceous_max, double t_equator, double t_pole,
                            const string &input_path ){
    this -> im = im;
    this -> jm = jm;
    this -> km = km;
    this -> i_max = i_max;
    this -> i_beg = i_beg;
    this -> Ma = Ma;
    this -> Ma_max = Ma_max;
    this -> Ma_max_half = Ma_max_half;
    this -> dr = dr;
    this -> g = g;
    this -> r_0_water = r_0_water;
    this -> ua = ua;
    this -> va = va;
    this -> wa = wa;
    this -> ta = ta;
    this -> ca = ca;
    this -> pa = pa;
    this -> u_0 = u_0;
    this -> p_0 = p_0;
    this -> t_0 = t_0;
    this -> c_0 = c_0;
    this -> cp_w = cp_w;
    this -> L_hyd = L_hyd;
    this -> t_average = t_average;
    this -> t_cretaceous_max = t_cretaceous_max;
    this -> t_equator = t_equator;
    this -> t_pole = t_pole;
    this -> input_path = input_path;

    j_half = ( jm - 1 ) / 2;
}


BC_Thermohalin::~BC_Thermohalin(){}


void BC_Thermohalin::IC_v_w_EkmanSpiral ( double Ma, Array_1D & rad, Array_1D & the,
                                    Array &h, Array &v, Array &w ){
    int j_30 = 30;
    int j_60 = 60;
    int j_90 = 90;
    int j_120 = 120;
    int j_150 = 150;

    pi180 = 180./M_PI;
    if ( Ma = 0 )  water_wind = .03;  // ocean surface velocity is about 3% of the wind velocity at the surface
    else           water_wind = 1.;

// Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
// a further turning downwards until the end of the shear layer such that finally 90° of turning are reached

//    Ekman_angle = 45.0 / pi180;
    Ekman_angle = 0.0;

// initial conditions for v and w velocity components at the sea surface
// ocean surface velocity is about 3% of the wind velocity at the surface
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_water( h, im-1, j, k ) ){
                v.x[ im - 1 ][ j ][ k ] = water_wind * v.x[ im - 1 ][ j ][ k ] / u_0;
                w.x[ im - 1 ][ j ][ k ] = water_wind * w.x[ im - 1 ][ j ][ k ] / u_0;
            }
        }
    }


// Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
// a further turning downwards until the end of the shear layer such that finally 90° of turning are reached
/*
// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° )
    for ( int j = 0; j < j_30; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_water( h, im-1, j, k ) ){
                vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] +
                    w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

                if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

                alfa = asin ( fabs ( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

                if ( w.x[ im-1 ][ j ][ k ] >= 0. ){
                    angle = + alfa - Ekman_angle;
                }else{
                    angle = - alfa - Ekman_angle;
                }

                if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }
            }
        }
    }

// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° )
    for ( int j = j_30; j < j_60; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_water( h, im-1, j, k ) ){
                vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] +
                    w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

                if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

                alfa = asin ( fabs ( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

                if ( v.x[ im-1 ][ j ][ k ] >= 0. ){
                    angle = + alfa - Ekman_angle;
                }else{
                    angle = - alfa - Ekman_angle;
                }

                if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = - vel_magnitude * sin ( angle );
                }

                if (  v.x[ im-1 ][ j ][ k ] <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = - vel_magnitude * sin ( angle );
                }
            }
        }
    }

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° )
    for ( int j = j_60; j < j_90; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_water( h, im-1, j, k ) ){
                vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] +
                     w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

                if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

                alfa = asin ( fabs ( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

                if ( w.x[ im-1 ][ j ][ k ] >= 0. ){
                    angle = + alfa - Ekman_angle;
                }else{
                    angle = - alfa - Ekman_angle;
                }

                if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }
            }
        }
    }

// south equatorial Hadley cell ( from j=90 till j=120 compares to 0° till 30° )
    for ( int j = j_90; j < j_120; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_water( h, im-1, j, k ) ){
                vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] +
                    w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

                if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

                beta = asin ( fabs( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

                if ( w.x[ im-1 ][ j ][ k ] >= 0. ){
                    angle = beta - Ekman_angle;
                }else{
                    angle = - beta - Ekman_angle;
                }

                if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }
            }
        }
    }

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° )
    for ( int j = j_120; j < j_150; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_water( h, im-1, j, k ) ){
                vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] +
                    w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

                if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

                beta = asin ( fabs( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

                if ( v.x[ im-1 ][ j ][ k ] <= 0. ){
                    angle = beta - Ekman_angle;
                }else{
                    angle = - beta - Ekman_angle;
                }

                if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = - vel_magnitude * sin ( angle );
                }

                if (  v.x[ im-1 ][ j ][ k ] >= 0. ){
                    v.x[ im-1 ][ j ][ k ] = + vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = - vel_magnitude * sin ( angle );
                }
            }
        }
    }

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° )
    for ( int j = j_150; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_water( h, im-1, j, k ) ){
                vel_magnitude = sqrt ( v.x[ im-1 ][ j ][ k ] * v.x[ im-1 ][ j ][ k ] 
                    + w.x[ im-1 ][ j ][ k ] * w.x[ im-1 ][ j ][ k ] );

                if ( w.x[ im-1 ][ j ][ k ] == 0. ) w.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( v.x[ im-1 ][ j ][ k ] == 0. ) v.x[ im-1 ][ j ][ k ] = 1.e-6;
                if ( vel_magnitude == 0. ) vel_magnitude = 1.e-6;

                beta = asin ( fabs( w.x[ im-1 ][ j ][ k ] ) / vel_magnitude );

                if ( w.x[ im-1 ][ j ][ k ] >= 0. ){
                    angle = beta - Ekman_angle;
                }else{
                    angle = - beta - Ekman_angle;
                }

                if ( w.x[ im-1 ][ j ][ k ] >= 0. && angle >= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] >= 0. && angle <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }

                if (  w.x[ im-1 ][ j ][ k ] <= 0. ){
                    v.x[ im-1 ][ j ][ k ] = - vel_magnitude * cos ( angle );
                    w.x[ im-1 ][ j ][ k ] = + vel_magnitude * sin ( angle );
                }
            }
        }
    }
*/

//    double i_Ekman_layer = 500.;                                    // assumed Ekman-layer depth of 500m
//    double i_Ekman_layer = 200.;                                    // assumed Ekman-layer depth of 200m
    double i_Ekman_layer = 100.;                                    // assumed Ekman-layer depth of 100m
    double coeff = i_Ekman_layer / L_hyd;

    int i_Ekman = ( im - 1 ) * ( 1. - coeff );

// surface wind vector driving the Ekman spiral in the Ekman layer
// northern hemisphere
    double gam_z = 0.;
    double exp_gam_z = 0.;
    double sin_gam_z = 0;
    double cos_gam_z = 0;
    double v_g = 0.;
    double w_g = 0.;

    for ( int j = 0; j < jm; j++ ){
        for ( int i = im-2; i >= i_Ekman; i-- ){
            gam_z = M_PI * ( double ) ( i - i_Ekman ) 
                / ( double ) ( im - 1 - i_Ekman );
            exp_gam_z = exp ( - gam_z );
            if ( j <= j_half ) sin_gam_z = sin ( gam_z );
            else               sin_gam_z = - sin ( gam_z );
            cos_gam_z = cos ( gam_z );
            for ( int k = 0; k < km; k++ ){
                if ( is_water( h, i, j, k) ){
                    v_g = v.x[ i + 1 ][ j ][ k ];
                    w_g = w.x[ i + 1 ][ j ][ k ];
                    v.x[ i ][ j ][ k ] = w_g * exp_gam_z * sin_gam_z 
                        + v_g * ( 1. - exp_gam_z * cos_gam_z );
                    w.x[ i ][ j ][ k ] = w_g * ( 1. - exp_gam_z * cos_gam_z ) 
                        - v_g * exp_gam_z * sin_gam_z;
                }
            }
        }
    }

    for ( int i = 0; i <= i_Ekman; i++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                v.x[ i ][ j ][ k ] = 0.;
                w.x[ i ][ j ][ k ] = 0.;
                }
            }
        }

    for ( int i = 0; i <= im-1; i++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                if ( is_land( h, i, j, k) ){
                    v.x[ i ][ j ][ k ] = 0.;
                    w.x[ i ][ j ][ k ] = 0.;
                }
            }
        }
    }
}






void BC_Thermohalin::BC_Temperature_Salinity ( Array_1D &rad, Array &h, Array &t, Array &c, Array &p_dyn ){
    // initial conditions for salt content and temperature and salinity decrease below the sea surface
    // boundary condition of  temperature and salinity
    // parabolic distribution from pole to pole accepted
    // maximum temperature and salinity of earth's surface at equator t_max = 1.137 compares to 37° C compares to 307 K
    // maximum temperature and salinity of earth's surface at equator t_max = 1.099 compares to 27° C compares to 297 K
    // polar temperature and salinity of earth's surface at north and south pole t_pol = 1.0 compares to -3° C compares to 270 K
    // minimum temperature and salinity at the tropopause t_min = .789 compares to -60° C compares to 213 K
    // minimum temperature and salinity in the deep ocean t_deep_ocean = 1.026 compares to 4° C compares to 277 K

    j_half = ( jm -1 ) / 2;
    j_max = jm - 1;

    d_j_half = ( double ) j_half;
    d_j_max = ( double ) j_max;

    // temperature-distribution by Ruddiman approximated by a parabola
    t_coeff = t_pole - t_equator;
    t_cretaceous_coeff = t_cretaceous_max / ( ( double ) Ma_max_half -
        ( double ) ( Ma_max_half * Ma_max_half / Ma_max ) );   // in °C
    t_cretaceous = t_cretaceous_coeff * ( double ) ( - ( Ma * Ma ) / Ma_max + Ma );   // in °C
    if ( Ma == 0 )  t_cretaceous = 0.;

    t_cretaceous = 0.;

    cout.precision ( 3 );

    time_slice_comment = "      time slice of Cretaceous-AGCM:";
    time_slice_number = " Ma = ";
    time_slice_unit = " million years";

    cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' )
        << time_slice_comment << setw ( 6 ) << std::fixed << setfill ( ' ' )
        << time_slice_number << setw ( 3 ) << Ma << setw ( 12 ) << time_slice_unit 
        << endl << endl;

    temperature_comment = "      temperature increase at cretaceous times: ";
    temperature_gain = " t increase";
    temperature_modern = "      mean temperature at modern times: ";
    temperature_cretaceous = "      mean temperature at cretaceous times: ";
    temperature_average = " t modern";
    temperature_average_cret = " t cretaceous";
    temperature_unit =  "°C ";

    cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << temperature_comment << 
        resetiosflags ( ios::left ) << setw ( 12 ) << temperature_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << 
        t_cretaceous << setw ( 5 ) << temperature_unit << endl << setw ( 50 ) << setfill ( '.' )  << 
        setiosflags ( ios::left ) << temperature_modern << resetiosflags ( ios::left ) << setw ( 13 ) << 
        temperature_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average << setw ( 5 ) << temperature_unit << 
        endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_cretaceous << 
        resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average_cret  << " = "  << setw ( 7 )  << setfill ( ' ' )
        << t_average + t_cretaceous << setw ( 5 ) << temperature_unit << endl;


    c_average = ( t_average + 346. ) / 10.;// in psu, relation taken from "Ocean Circulation, The Open University"
    c_cretaceous = ( t_average + t_cretaceous + 346. ) / 10.;// in psu
    c_cretaceous = c_cretaceous - c_average;
    if ( Ma == 0 )  c_cretaceous = 0.;

    salinity_comment = "      salinity increase at cretaceous times: ";
    salinity_gain = " c increase";
    salinity_modern = "      mean salinity at modern times: ";
    salinity_cretaceous = "      mean salinity at cretaceous times: ";
    salinity_average = " c modern";
    salinity_average_cret = " c cretaceous";
    salinity_unit =  "psu ";

    cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << salinity_comment << 
        resetiosflags ( ios::left ) << setw ( 12 ) << salinity_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << 
        c_cretaceous << setw ( 5 ) << salinity_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) 
        << salinity_modern << resetiosflags ( ios::left ) << setw ( 13 ) << salinity_average  << " = "  << setw ( 7 )  << 
        setfill ( ' ' ) << c_average << setw ( 5 ) << salinity_unit << endl << setw ( 50 ) << setfill ( '.' )  << 
        setiosflags ( ios::left ) << salinity_cretaceous << resetiosflags ( ios::left ) << setw ( 13 ) << 
        salinity_average_cret  << " = "  << setw ( 7 )  << setfill ( ' ' ) << c_average + c_cretaceous << setw ( 5 ) << 
        salinity_unit << endl;
    cout << endl;

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            if ( is_water( h, im-1, j, k ) ){
                t.x[ im-1 ][ j ][ k ] = t.x[ im-1 ][ j ][ k ] + t_cretaceous; // paleo surface temperature
                c.x[ im-1 ][ j ][ k ] = c.x[ im-1 ][ j ][ k ] + c_cretaceous / c_0;  // non-dimensional

                t_Celsius = t.x[ im-1 ][ j ][ k ] * t_0 - t_0;
/*
                if ( t_Celsius <= ( ta * t_0 - t_0 ) ){ // water temperature not below -1°C
                    t.x[ im-1 ][ j ][ k ] = ta;
                    c.x[ im-1 ][ j ][ k ] = ca;
                }
*/
            }
        }
    }

    double tm_tbeg = 0.;
    double cm_cbeg = 0.;
    double d_i_max = ( double ) i_max;
    double d_i_beg = ( double ) i_beg;
    double d_i = 0.;

// distribution of t and c with increasing depth till i_beg valid for all time slices including the modern world
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            tm_tbeg = ( t.x[ im-1 ][ j ][ k ] - ta ) 
                / ( d_i_max * d_i_max - d_i_beg * d_i_beg );
            cm_cbeg = ( c.x[ im-1 ][ j ][ k ] - ca ) 
                / ( d_i_max * d_i_max - d_i_beg * d_i_beg );
            for ( int i = i_beg; i < im-1; i++ ){
                if ( is_water( h, i, j, k ) ){
                    d_i = ( double ) i;
                    t.x[ i ][ j ][ k ] = ta + tm_tbeg 
                        * ( d_i * d_i - d_i_beg * d_i_beg );// parabolic approach
                    t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;
                    c.x[ i ][ j ][ k ] = ca + cm_cbeg 
                        * ( d_i * d_i - d_i_beg * d_i_beg );// parabolic approach
/*
                    if ( t_Celsius <= ( ta * t_0 - t_0 ) ){ // water temperature not below 4°C
                        t.x[ i ][ j ][ k ] = ta;
                        c.x[ i ][ j ][ k ] = ca;
                    }
*/
                }else{
                    t.x[ i ][ j ][ k ] = ta;
                    c.x[ i ][ j ][ k ] = ca;
                }
            }
        }
    }

    for ( int i = 0; i < i_beg; i++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                if ( is_water( h, i, j, k ) ){
                    t.x[ i ][ j ][ k ] = ta;
                    c.x[ i ][ j ][ k ] = ca;
                }
            }
        }
    }
}




void BC_Thermohalin::BC_Pressure_Density ( Array_1D &rad, Array &p_stat, Array &r_water,
                                    Array &r_salt_water, Array &t, Array &c, Array &h ){
// hydrostatic pressure, equations of state for water and salt water density
// as functions of salinity, temperature and hydrostatic pressure

    double t_Celsius_0 = 0.;
    double t_Celsius_1 = 0.;
    double p_km = 0.;
    double C_p = 0.;
    double beta_p =  0.;
    double alfa_t_p =  0.;
    double gamma_t_p =  0.;

    double E_water = 2.15e9;                            // given in N/m²
    double beta_water = 8.8e-5;                         // given in m³/( m³ * °C)
    double r_air = 1.2041;                              // given in kg/m³
    double R_Air = 287.1;                               // given in J/( kg*K )

// hydrostatic pressure, water and salt water density at the surface
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            p_stat.x[ im-1 ][ j ][ k ] =  .01 * ( r_air * R_Air 
                * t.x[ im-1 ][ j ][ k ] * t_0 ) / 1000.;
            // given in bar, isochoric approach, constant air density at the surface
            r_water.x[ im-1 ][ j ][ k ] = r_0_water;                // given in kg/m³

            t_Celsius_1 = t.x[ im-1 ][ j ][ k ] * t_0 - t_0;
            p_km = 0.;
            C_p = 999.83;
            beta_p = .808;
            alfa_t_p = .0708 * ( 1. + .068 * t_Celsius_1 );
            gamma_t_p = .003 * ( 1. - .012 * t_Celsius_1 );

            r_salt_water.x[ im-1 ][ j ][ k ] = C_p + beta_p 
                * c.x[ im-1 ][ j ][ k ] * c_0 - alfa_t_p * t_Celsius_1 
                - gamma_t_p * ( c_0 - c.x[ im-1 ][ j ][ k ] * c_0 ) 
                * t_Celsius_1;
        }
    }

// hydrostatic pressure, water and salt water density in the field
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = im-2; i >= 0; i-- ){
                d_i = ( double ) ( im - 1 - i );

                t_Celsius_1 = t.x[ i ][ j ][ k ] * t_0 - t_0;
                t_Celsius_0 = t.x[ i + 1 ][ j ][ k ] * t_0 - t_0;

                p_stat.x[ i ][ j ][ k ] = r_water.x[ i + 1 ][ j ][ k ] 
                    * g * d_i * ( L_hyd / ( double ) ( im-1 ) ) / 100000. 
                    + p_0 / 1000.;             // hydrostatic pressure in bar
                r_water.x[ i ][ j ][ k ] = r_water.x[ i + 1 ][ j ][ k ] 
                    / ( 1. + beta_water * ( t_Celsius_1 - t_Celsius_0 ) ) 
                    / ( 1. - ( p_stat.x[ i ][ j ][ k ] -
                    p_stat.x[ i+1 ][ j ][ k ] ) / E_water * 1e5 );  // given in kg/m³

                p_km  = ( double ) ( im - 1 - i ) * ( L_hyd 
                    / ( double ) ( im-1 ) ) / 1000.;
                C_p = 999.83 + 5.053 * p_km - .048 * p_km * p_km;
                beta_p = .808 - .0085* p_km;
                alfa_t_p = .0708 * ( 1. + .351 * p_km + .068 
                    * ( 1. - .0683 * p_km ) * t_Celsius_1 );
                gamma_t_p = .003 * ( 1. - .059 * p_km - .012 
                    * ( 1. - .064 * p_km ) * t_Celsius_1 );

                r_salt_water.x[ i ][ j ][ k ] = C_p + beta_p 
                    * c.x[ i ][ j ][ k ] * c_0 - alfa_t_p * t_Celsius_1 
                    - gamma_t_p * ( c_0 - c.x[ i ][ j ][ k ] * c_0 ) 
                    * t_Celsius_1;
            }
        }
    }
}






void BC_Thermohalin::BC_Surface_Temperature_NASA
                     ( const string &Name_SurfaceTemperature_File, Array &t ){
    // initial conditions for the temperature and salinity at the sea surface

    cout.precision ( 3 );
    cout.setf ( ios::fixed );

    ifstream Name_SurfaceTemperature_File_Read;
    Name_SurfaceTemperature_File_Read.open(Name_SurfaceTemperature_File);

    if (!Name_SurfaceTemperature_File_Read.is_open()) {
        cerr << "ERROR: could not open SurfaceTemperature_File file at " 
            << Name_SurfaceTemperature_File << "\n";
        abort();
    }

    j = 0;
    k = 0;

    while ( ( k < km ) && ( !Name_SurfaceTemperature_File_Read.eof() ) ){
        while ( j < jm ){
            Name_SurfaceTemperature_File_Read >> dummy_1;
            Name_SurfaceTemperature_File_Read >> dummy_2;
            Name_SurfaceTemperature_File_Read >> dummy_3;

            t.x[ im-1 ][ j ][ k ] = ( dummy_3 + 273.15 ) / 273.15;

            j++;
        }
        j = 0;
        k++;
    }

    Name_SurfaceTemperature_File_Read.close();

    for ( int j = 0; j < jm; j++ ){
        for ( int k = 1; k < km-1; k++ ){
            if ( k == 180 ) t.x[ im-1 ][ j ][ k ] = ( t.x[ im-1 ][ j ][ k + 1 ] 
                + t.x[ im-1 ][ j ][ k - 1 ] ) * .5;
        }
    }
}






void BC_Thermohalin::BC_Surface_Salinity_NASA 
                     ( const string &Name_SurfaceSalinity_File, Array &c ){
    // initial conditions for the salinity at the sea surface
    streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, 
        endpos_3, anfangpos_4, endpos_4;

    cout.precision ( 3 );
    cout.setf ( ios::fixed );

    ifstream Name_SurfaceSalinity_File_Read;
    Name_SurfaceSalinity_File_Read.open(Name_SurfaceSalinity_File);

    if (!Name_SurfaceSalinity_File_Read.is_open()){
        cerr << "ERROR: could not open SurfaceSalinity_File file at " 
        << Name_SurfaceSalinity_File << "\n";
        abort();
    }

    j = 0;
    k = 0;

    while ( ( k < km ) && ( !Name_SurfaceSalinity_File_Read.eof() ) ){
        while ( j < jm ){
            Name_SurfaceSalinity_File_Read >> dummy_1;
            Name_SurfaceSalinity_File_Read >> dummy_2;
            Name_SurfaceSalinity_File_Read >> dummy_3;

            if ( dummy_3 < 0. ) dummy_3 = ca;

            else        c.x[ im-1 ][ j ][ k ] = dummy_3 / c_0;
            j++;
        }
        j = 0;
        k++;
    }

    Name_SurfaceSalinity_File_Read.close();
}





void BC_Thermohalin::BC_Surface_v_Velocity_Ocean ( const string &Name_v_surface_ocean_File, Array &v ){
// initial conditions for the Name_v_surface_ocean_File at the sea surface
    cout.precision ( 3 );
    cout.setf ( ios::fixed );
    ifstream Name_v_surface_ocean_File_Read(Name_v_surface_ocean_File);
    if (!Name_v_surface_ocean_File_Read.is_open()){
        cerr << "ERROR: could not open Name_v_surface_ocean_File file at "
        << Name_v_surface_ocean_File << "\n";
        abort();
    }
    int k_half = ( km - 1 ) / 2;  // position at 180°E ( Greenwich )
    j = 0;
    k = 0;
    while ( ( k < km ) && !Name_v_surface_ocean_File_Read.eof() ){
        while ( j < jm ){
            double lat, lon, v_velocity;
            Name_v_surface_ocean_File_Read >> lat;
            Name_v_surface_ocean_File_Read >> lon;
            Name_v_surface_ocean_File_Read >> v_velocity;
            v.x[ im-1 ][ j ][ k ] = v_velocity;
            j++;
        }
    j = 0;
    k++;
    }
}




void BC_Thermohalin::BC_Surface_w_Velocity_Ocean ( const string &Name_w_surface_ocean_File, Array &w ){
// initial conditions for the Name_w_surface_ocean_File at the sea surface
    cout.precision ( 3 );
    cout.setf ( ios::fixed );
    ifstream Name_w_surface_ocean_File_Read(Name_w_surface_ocean_File);
    if (!Name_w_surface_ocean_File_Read.is_open()){
        cerr << "ERROR: could not open Name_w_surface_ocean_File file at "
        << Name_w_surface_ocean_File << "\n";
        abort();
    }
    int k_half = ( km - 1 ) / 2;  // position at 180°E ( Greenwich )
    j = 0;
    k = 0;
    while ( ( k < km ) && !Name_w_surface_ocean_File_Read.eof() ){
        while ( j < jm ){
            double lat, lon, w_velocity;
            Name_w_surface_ocean_File_Read >> lat;
            Name_w_surface_ocean_File_Read >> lon;
            Name_w_surface_ocean_File_Read >> w_velocity;
            w.x[ im-1 ][ j ][ k ] = w_velocity;
            j++;
        }
    j = 0;
    k++;
    }
}





void BC_Thermohalin::IC_Equatorial_Currents 
                     ( Array &h, Array &u, Array &v, Array &w ){
// currents along the equator
// equatorial undercurrent - Cromwell flow, EUC
// equatorial intermediate current, EIC
// nothern and southern equatorial subsurface counter-currents, NSCC und SSCC
// nothern and southern equatorial counter-currents, NECC und SECC

    double IC_water = 1.;  // no dimension, ( average velocity compares to   u_0 * IC_water = 0,25 )
/*
// one grid step compares to a depth of 25 m for L_hyd = 1000m   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    i_EIC_u = i_beg; // 1000m depth  // westward equatorial intermediate current, max 0.1 m/s
    i_EIC_o = 28; // 0m depth
    i_SCC_u = 8; // 800m depth  // eastward subsurface counter current, max 0.05 m/s
    i_SCC_o = 28; // 300m depth
    i_ECC_u = 28; // 400m depth  // eastward equatorial counter current, max 0.2 m/s
    i_ECC_o = im; // 0m depth
    i_EUC_u = 32; // 200m depth  // eastward equatorial under current ( Cromwell current ), max 0.8 m/s
    i_EUC_o = 36; // 100m depth
*/
// one grid step compares to a depth of 5 m for L_hyd = 200m   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    i_EIC_u = i_beg; // 200m depth  // westward equatorial intermediate current, max 0.1 m/s
    i_EIC_o = 0; // 0m depth
    i_SCC_u = 0; // 800m depth  // eastward subsurface counter current, max 0.05 m/s
    i_SCC_o = 0; // 300m depth
    i_ECC_u = 0; // 200m depth  // eastward equatorial counter current, max 0.2 m/s
    i_ECC_o = im; // 0m depth
    i_EUC_u = 0; // 200m depth  // eastward equatorial under current ( Cromwell current ), max 0.8 m/s
    i_EUC_o = 20; // 100m depth

// equatorial currents and counter-currents

//  §§§§§§§§§§§§§§§§§§§§§§§§§   valid for all paleo-ocean constellations along the equator   §§§§§§§§§§§§§§§§§§§§§§§§§§

    j_half = ( jm -1 ) / 2;
    k_beg = 0;
    k_end = 0;

// extention of land and ocean areas
    for ( int k = 1; k < km-1; k++ ){
        if ( k < km ){
            if ( ( is_water( h, im-1, j_half, k ) ) 
                && ( is_land( h, im-1, j_half, k+1 ) ) ){
                k_end = k;
            }
            if ( ( is_land( h, im-1, j_half, k ) ) 
                && ( is_water( h, im-1, j_half, k+1 ) ) ){
                k_beg = k;
            }
            if ( ( ( is_water( h, im-1, j_half, k ) ) 
                && ( is_water( h, im-1, j_half, k+1 ) ) ) 
                && ( k == km-2 ) ){
                k_end = k;
            }

// equatorial northern counter-current ( NECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial northern counter-current ( from j=83 until j=88 compares to 3°N until 8°N )
            for ( int i = i_ECC_u; i < i_ECC_o; i++ ){
                for ( int j =83; j < 88; j++ ){
                    for ( int k = k_beg; k <= k_end; k++ ){
                        if ( is_water( h, i, j, k ) ){
                            v.x[ i ][ j ][ k ] = 0.;
                            w.x[ i ][ j ][ k ] = + .01 * IC_water 
                                * ( double ) ( i - i_ECC_u ) 
                                / ( double ) ( i_ECC_o - i_ECC_u );
                        }
                    }
                }
            }

// equatorial southern counter-current ( SECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial southern counter-current ( from j=93 until j=96 compares to 3°S until 6°S )
            for ( int i = i_ECC_u; i < i_ECC_o; i++ ){
                for ( int j = 93; j < 98; j++ ){
                    for ( int k = k_beg; k <= k_end; k++ ){
                        if ( is_water( h, i, j, k ) ){
                            v.x[ i ][ j ][ k ] = 0.;
                            w.x[ i ][ j ][ k ] = + .01 * IC_water 
                                * ( double ) ( i - i_ECC_u ) 
                                / ( double ) ( i_ECC_o - i_ECC_u );
                        }
                    }
                }
            }

// equatorial undercurrent - Cromwell current ( EUC, i=im-2 until i=im-1 compares to -100 until -200m depth )
// equatorial undercurrent - Cromwell current ( from j=87 until j=93 compares to 3°N until 3°S )
            for ( int i = i_EUC_u; i < i_EUC_o; i++ ){
                for ( int j = 87; j < 94; j++ ){
                    for ( int k = k_beg; k <= k_end; k++ ){
                        if ( is_water( h, i, j, k ) ){
                            v.x[ i ][ j ][ k ] = 0.;
                            w.x[ i ][ j ][ k ] = + .4 * IC_water;
                        }
                    }
                }
            }

// equatorial intermediate current ( EIC, i=im-4 until i=im-2 compares to -300 until -1000m depth )
// equatorial intermediate current ( from j=88 until j=92 compares to 2°N until 2°S )
            for ( int i = i_EIC_u; i < i_EIC_o; i++ ){
                for ( int j = 88; j < 93; j++ ){
                    for ( int k = k_beg; k <= k_end; k++ ){
                        if ( is_water( h, i, j, k ) ){
                            v.x[ i ][ j ][ k ] = 0.;
                            w.x[ i ][ j ][ k ] = - .05 * IC_water 
                                * ( double ) ( i - i_EIC_u ) 
                                / ( double ) ( i_EIC_o - i_EIC_u );
                        }
                    }
                }
            }

// equatorial northern and southern subsurface counter-current
// equatorial northern subsurface counter-current ( NSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial northern subsurface counter-current ( from j=86 until j=87 compares to 3°N until 4°N )
            for ( int i = i_SCC_u; i < i_SCC_o; i++ ){
                for ( int j = 86; j < 88; j++ ){
                    for ( int k = k_beg; k <= k_end; k++ ){
                        if ( is_water( h, i, j, k ) ){
                            v.x[ i ][ j ][ k ] = 0.;
                            w.x[ i ][ j ][ k ] = .025 * IC_water 
                                * ( double ) ( i - i_SCC_u ) 
                                / ( double ) ( i_SCC_o - i_SCC_u );
                        }
                    }
                }
            }

// equatorial southern subsurface counter-current ( SSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial southern subsurface counter-current ( from j=93 until j=94 compares to 3°S until 4°S )
            for ( int i = i_SCC_u; i < i_SCC_o; i++ ){
                for ( int j = 93; j < 95; j++ ){
                    for ( int k = k_beg; k <= k_end; k++ ){
                        if ( is_water( h, i, j, k ) ){
                            v.x[ i ][ j ][ k ] = 0.;
                            w.x[ i ][ j ][ k ] = .025 * IC_water 
                                * ( double ) ( i - i_SCC_u ) 
                                / ( double ) ( i_SCC_o - i_SCC_u );
                        }
                    }
                }
            }
        }
    }
}






void BC_Thermohalin::IC_CircumPolar_Current 
                     ( Array &h, Array &u, Array &v, Array &w, Array &c ){
// south polar sea
// antarctic circumpolar current ( -5000m deep ) ( from j=147 until j=152 compares to 57°S until 62°S,
// from k=0 until k=km compares to 0° until 360° )
    for ( int i = i_beg; i < im; i++ ){
        for ( int j = 147; j < 153; j++ ){
            for ( int k = 0; k < km; k++ ){
                if ( is_water( h, i, j, k ) ){
//                  c.x[ i ][ j ][ k ] = ca;
                    w.x[ i ][ j ][ k ] = 2. * u_0;  // 0.5 m/s
                }
            }
        }
    }
}




void BC_Thermohalin::IC_u_WestEastCoast 
                     ( Array_1D &rad, Array &h, Array &u, Array &v, 
                     Array &w, Array &un, Array &vn, Array &wn ){
// initial conditions for v and w velocity components at the sea surface close to east or west coasts
// reversal of v velocity component between north and south equatorial current ommitted at respectively 10°
// w component unchanged

    int i_beg = 20;  // == 500m depth
    int i_middle = i_beg + 16;  // 0 + 36 = 36 for total depth 100m ( i_beg= 0 ), asymmetric with depth for stepsizes of 25m
    int i_half = 38;
    d_i_middle = ( double ) i_middle;
    d_i_half = ( double ) i_half;
    d_i_max = ( double ) i_max;

// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: east coast
    k_grad = 10;                                                                            // extension of velocity change
    k_a = k_grad;                                                                        // left distance
    k_b = 0;                                                                                // right distance

    k_water = 0;                                                                        // on water closest to coast
    k_sequel = 1;                                                                        // on solid ground k_sequel = 0

    for ( int j = 0; j < 91; j++ ){                                                    // outer loop: latitude
        for ( int k = 0; k < km; k++ ){                                            // inner loop: longitude
            if ( is_land( h, i_half, j, k ) ) k_sequel = 0;                // if solid ground: k_sequel = 0
            if ( ( is_water( h, i_half, j, k ) ) && ( k_sequel == 0 ) ) 
                k_water = 0;    // if water and and k_sequel = 0 then is water closest to coast
            else k_water = 1;                                                        // somewhere on water

            if ( ( is_water( h, i_half, j, k ) ) && ( k_water == 0 ) ){    // if water is closest to coast, change of velocity components begins
                for ( int l = 0; l < k_grad; l++ ){                                // extension of change, sign change in v-velocity and distribution of u-velocity with depth
                    for ( int i = i_beg; i < i_half; i++ ){                    // loop in radial direction, extension for u -velocity component, downwelling here
                        m = i + ( i_half - i_middle );
                        d_i = ( double ) i;
                        u.x[ i ][ j ][ k + l ] = - d_i / d_i_half * water_wind 
                            / ( ( double )( l + 1 ) );    // increase with depth, decrease with distance from coast
                        u.x[ m ][ j ][ k + l ] = - d_i / d_i_middle * water_wind 
                            / ( ( double )( l + 1 ) );// decrease with depth, decrease with distance from coast
                    }
                }
                k_sequel = 1;                                                            // looking for another east coast
            }
        }                                                                                        // end of longitudinal loop
        k_water = 0;                                                                    // starting at another latitude
    }                                                                                            // end of latitudinal loop


// southern hemisphere: east coast
    k_water = 0;
    k_sequel = 1;

    for ( int j = 91; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_land( h, i_half, j, k ) ) k_sequel = 0;
            if ( ( is_water( h, i_half, j, k ) ) && ( k_sequel == 0 ) ) 
                k_water = 0;
            else k_water = 1;
            if ( ( is_water( h, i_half, j, k ) ) && ( k_water == 0 ) ){
                for ( int l = 0; l < k_grad; l++ ){
                    for ( int i = i_beg; i < i_half; i++ ){
                        m = i + ( i_half - i_middle );
                        d_i = ( double ) i;
                        u.x[ i ][ j ][ k + l ] = - d_i / d_i_half * water_wind 
                            / ( ( double )( l + 1 ) );    // increase with depth, decrease with distance from coast
                        u.x[ m ][ j ][ k + l ] = - d_i / d_i_middle * water_wind 
                            / ( ( double )( l + 1 ) );// decrease with depth, decrease with distance from coast
                    }
                }
                k_sequel = 1;
            }
        }
        k_water = 0;
    }


// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: west coast
    k_a = 0;                                                                                // left distance
    k_water = 0;                                                                        // somewhere on water
    flip = 0;                                                                                // somewhere on water

    for ( int j = 0; j < 91; j++ ){                                                    // outer loop: latitude
        for ( int k = 0; k < km; k++ ){                                            // inner loop: longitude
            if ( is_water( h, i_half, j, k ) ){                                        // if somewhere on water
                k_water = 0;                                                            // somewhere on water: k_water = 0
                flip = 0;                                                                    // somewhere on water: flip = 0
            }
            else k_water = 1;                                                        // first time on land

            if ( ( flip == 0 ) && ( k_water == 1 ) ){                            // on water closest to land
                for ( int l = k; l > ( k - k_grad + 1 ); l-- ){                        // backward extention of velocity change: nothing changes
                    for ( int i = i_beg; i < i_half; i++ ){                    // loop in radial direction, extension for u -velocity component, downwelling here
                        m = i + ( i_half - i_middle );
                        d_i = ( double ) i;
                        u.x[ i ][ j ][ l ] = + d_i / d_i_half * water_wind 
                            / ( ( double )( k - l + 1 ) );    // increase with depth, decrease with distance from coast
                        u.x[ m ][ j ][ l ] = + d_i / d_i_middle * water_wind 
                            / ( ( double )( k - l + 1 ) );    // decrease with depth, decrease with distance from coast
                    }
                }
                flip = 1;
            }
        }
        flip = 0;
    }


// southern hemisphere: west coast
    k_water = 0;
    flip = 0;

    for ( int j = 91; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            if ( is_water( h, i_half, j, k ) ){
                k_water = 0;
                flip = 0;
            }
            else k_water = 1;

            if ( ( flip == 0 ) && ( k_water == 1 ) ){
                for ( int l = k; l > ( k - k_grad + 1 ); l-- ){
                    for ( int i = i_beg; i < i_half; i++ ){
                        m = i + ( i_half - i_middle );
                        d_i = ( double ) i;
                        u.x[ i ][ j ][ l ] = + d_i / d_i_half * water_wind 
                            / ( ( double )( k - l + 1 ) );
                        u.x[ m ][ j ][ l ] = + d_i / d_i_middle * water_wind 
                            / ( ( double )( k - l + 1 ) );
                    }
                }
                flip = 1;
            }
        }
        flip = 0;
    }

    for ( int i = 0; i < i_beg; i++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                u.x[ i ][ j ][ k ] = 0.;
                }
            }
        }
    }











void BC_Thermohalin::Value_Limitation_Hyd ( Array &h, Array &u, Array &v,
                                    Array &w, Array &p_dyn, Array &t, Array &c ){
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
                if ( u.x[ i ][ j ][ k ] >= 11.11 )  u.x[ i ][ j ][ k ] = 11.11;
                if ( u.x[ i ][ j ][ k ] <= - 11.11 )  u.x[ i ][ j ][ k ] = - 11.11;

                if ( v.x[ i ][ j ][ k ] >= .552 )  v.x[ i ][ j ][ k ] = .552;
                if ( v.x[ i ][ j ][ k ] <= - .552 )  v.x[ i ][ j ][ k ] = - .552;

                if ( w.x[ i ][ j ][ k ] >= .552 )  w.x[ i ][ j ][ k ] = .552;
                if ( w.x[ i ][ j ][ k ] <= - .552 )  w.x[ i ][ j ][ k ] = - .552;

                if ( t.x[ i ][ j ][ k ] >= 1.147 )  t.x[ i ][ j ][ k ] = 1.147;      // 40.15 °C
                if ( t.x[ i ][ j ][ k ] <= 0.9963 )  t.x[ i ][ j ][ k ] = 0.9963;    // -1.0 °C

                if ( c.x[ i ][ j ][ k ] >= 1.156 )  c.x[ i ][ j ][ k ] = 1.156;    // 40.0 psu
                if ( c.x[ i ][ j ][ k ] <= .95 )  c.x[ i ][ j ][ k ] = .95;           // 32.0 psu
            }
        }
    }
}



void BC_Thermohalin::Pressure_Limitation_Hyd ( Array &p_dyn, Array &p_dynn ){
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
                if ( p_dyn.x[ i ][ j ][ k ] >= .1 )  p_dyn.x[ i ][ j ][ k ] = .1;
                if ( p_dyn.x[ i ][ j ][ k ] <= - .1 )  p_dyn.x[ i ][ j ][ k ] = - .1;

                p_dynn.x[ i ][ j ][ k ] = p_dyn.x[ i ][ j ][ k ];
            }
        }
    }
}



void BC_Thermohalin::BC_Evaporation ( Array_2D &salinity_evaporation, Array_2D &Evaporation_Dalton, 
                     Array_2D &Precipitation, Array &h, Array &c, Array &r_water ){
// preparations for salinity increase due to evaporation and precipitation differences
// procedure given in Rui Xin Huang, Ocean Circulation, p. 165
    double salinity_surface = 0.;
    double coeff_salinity = 1.1574e-5 * L_hyd / u_0 * 1.e-3 * c_0;
                          // 1.1574-5 is the conversion from (Evap-Prec) in mm/d to m/s
                          // 1.e-3 * c_0 stands for the transformation 
                          // from non-dimensional salinity to salt mass fraction

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
    // the 1. order gradients are built at 1 level below sea surface, at sea level values diverge
    // this formula contains a 2. order accurate gradient of 1. order, needs 3 points, not always available along coasts
    // gradient formed 1 point below surface, causes problems when formed at the surface, which normally is correct
//            salinity_surface = - ( - 3. * c.x[ im - 2 ][ j ][ k ] +
//                               4. * c.x[ im - 3 ][ j ][ k ] - c.x[ im - 4 ][ j ][ k ] ) / ( 2. * dr ) *
//                               ( 1. - 2. * c.x[ im - 2 ][ j ][ k ] ) * 
//                               ( Evaporation_Dalton.y[ j ][ k ] - Precipitation.y[ j ][ k ] );

    // the 1. order gradients are built at 1 level below sea surface, at sea level values diverge
    // this formula contains a 1. order accurate gradient of 1. order, needs 2 points, in general available along coasts
    // gradient formed 1 point below surface, causes problems when formed at the surface, which normally is correct
            salinity_surface = - ( c.x[ im - 2 ][ j ][ k ] - c.x[ im - 3 ][ j ][ k ] ) / dr *
                       ( 1. - 2. * c.x[ im - 2 ][ j ][ k ] ) * 
                       ( Evaporation_Dalton.y[ j ][ k ] - Precipitation.y[ j ][ k ] );

            salinity_evaporation.y[ j ][ k ] = coeff_salinity * salinity_surface;

            if ( is_land( h, im-1, j, k) )    salinity_evaporation.y[ j ][ k ] = 0.;

//            salinity_evaporation.y[ j ][ k ] = 0.;    // test case

//            c.x[ im-1 ][ j ][ k ] = c.x[ im-1 ][ j ][ k ] + salinity_evaporation.y[ j ][ k ];

            if ( c.x[ im-1 ][ j ][ k ] >= 1.156 )  c.x[ im-1 ][ j ][ k ] = 1.156;  // 40.0 psu
            if ( c.x[ im-1 ][ j ][ k ] <= .95 )  c.x[ im-1 ][ j ][ k ] = .95;      // 32.0 psu


    cout.precision ( 8 );
    cout.setf ( ios::fixed );

//    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   j = " << j << "   k = " << k << "   salinity_evaporation = " << salinity_evaporation.y[ j ][ k ] << "   coeff_salinity = " << coeff_salinity << "   salinity_surface = " << salinity_surface << "   c = " << c.x[ im-1 ][ j ][ k ] << "   c_0 = " << c_0 << "   r_0_water = " << r_0_water << "   c * c_0 = " << c.x[ im-1 ][ j ][ k ] * c_0 << "   Evap-Prec = " << ( Evaporation_Dalton.y[ j ][ k ] - Precipitation.y[ j ][ k ] ) << "   Evap = " << Evaporation_Dalton.y[ j ][ k ] << "   Prec = " << Precipitation.y[ j ][ k ] << "   c_grad_1 = " << ( c.x[ im - 2 ][ j ][ k ] - c.x[ im - 3 ][ j ][ k ] ) / dr << "   c_grad_2 = " << - ( - 3. * c.x[ im - 2 ][ j ][ k ] + 4. * c.x[ im - 3 ][ j ][ k ] - c.x[ im - 4 ][ j ][ k ] ) / ( 2. * dr ) << endl;

        }
    }
}
