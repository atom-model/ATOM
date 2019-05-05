/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to produce resulting data on mean sea level
*/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "Results_Atm.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

Results_MSL_Atm::Results_MSL_Atm ( int im, int jm, int km, int sun, double g, double ep,
                            double hp, double u_0, double p_0, double t_0, double c_0, double co2_0,
                            double sigma, double albedo_equator, double lv, double ls, double cp_l,
                            double L_atm, double dt, double dr, double dthe, double dphi, double r_air,
                            double R_Air, double r_water_vapour, double R_WaterVapour,
                            double co2_vegetation, double co2_ocean, double co2_land,
                            double gam, double t_pole, double t_cretaceous, double t_average ):
f_Penman ( 2. ){
    this-> im = im;
    this-> jm = jm;
    this-> km = km;
    this-> L_atm = L_atm;
    this-> dt = dt;
    this-> dr = dr;
    this-> dthe = dthe;
    this-> dphi = dphi;
    this-> sun = sun;
    this-> g = g;
    this-> ep = ep;
    this-> hp = hp;
    this-> u_0 = u_0;
    this-> p_0 = p_0;
    this-> t_0 = t_0;
    this-> c_0 = c_0;
    this-> co2_0 = co2_0;
    this-> sigma = sigma;
    this-> albedo_equator = albedo_equator;
    this-> lv = lv;
    this-> ls = ls;
    this-> gam = gam;
    this-> cp_l = cp_l;
    this-> r_air = r_air;
    this-> R_Air = R_Air;
    this-> R_WaterVapour = R_WaterVapour;
    this-> r_water_vapour = r_water_vapour;
    this-> co2_vegetation = co2_vegetation;
    this-> co2_ocean = co2_ocean;
    this-> co2_land = co2_land;
    this-> t_pole = t_pole;
    this-> t_cretaceous = t_cretaceous;
    this-> t_average = t_average;

// array "vel_av" for velocity averaging
    vel_av = 0L;
    vel_av = new double[ im ];

    for ( int l = 0; l < im; l++ ){
        vel_av[ l ] = 0.;
    }

    coeff_mmWS = r_air / r_water_vapour;                                    // coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]
    coeff_lv = lv / ( cp_l * t_0 );                                                        // coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_lv = 9.1069 in [ / ]
    coeff_ls = ls / ( cp_l * t_0 );                                                        // coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_ls = 10.3091 in [ / ]

// array "temp_av" for air temperature averaging
    temp_av = 0L;
    temp_av = new double[ im ];

    for ( int l = 0; l < im; l++ ){
        temp_av[ l ] = 0.;
    }

// array "lat_" for air temperature averaging
    lat_av = 0L;
    lat_av = new double[ im ];

    for ( int l = 0; l < im; l++ ){
        lat_av[ l ] = 0.;
    }

// array "sen_av" for air temperature averaging
    sen_av = 0L;
    sen_av = new double[ im ];

    for ( int l = 0; l < im; l++ ){
        sen_av[ l ] = 0.;
    }

    coeff_mmWS = r_air / r_water_vapour;                                    // coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]
    coeff_lv = lv / ( cp_l * t_0 );                                                        // coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_lv = 9.1069 in [ / ]
    coeff_ls = ls / ( cp_l * t_0 );                                                        // coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_ls = 10.3091 in [ / ]

    c43 = 4./3.;
    c13 = 1./3.;
}


Results_MSL_Atm::~Results_MSL_Atm (){
        delete [  ] vel_av;
        delete [  ] temp_av;
        delete [  ] lat_av;
        delete [  ] sen_av;
}



void Results_MSL_Atm::run_MSL_data ( int n, int velocity_iter_max, int RadiationModel,
                                    double &t_cretaceous, Array_1D &rad, Array_1D &the, Array_1D &phi,
                                    Array &h, Array &c, Array &cn, Array &co2, Array &co2n, Array &t,
                                    Array &tn, Array &p_dyn, Array &p_stat, Array &BuoyancyForce,
                                    Array &u, Array &v, Array &w, Array &Q_Latent, Array &Q_Sensible,
                                    Array &radiation_3D, Array &cloud, Array &cloudn, Array &ice,
                                    Array &icen, Array &P_rain, Array &P_snow, Array &aux_u,
                                    Array &aux_v, Array &aux_w, Array_2D &temperature_NASA,
                                    Array_2D &precipitation_NASA, Array_2D &precipitable_water,
                                    Array_2D &Q_radiation, Array_2D &Q_Evaporation, Array_2D &Q_latent,
                                    Array_2D &Q_sensible, Array_2D &Q_bottom, Array_2D &Evaporation_Penman,
                                    Array_2D &Evaporation_Dalton, Array_2D &Vegetation, Array_2D &albedo,
                                    Array_2D &co2_total, Array_2D &Precipitation, Array &S_v, Array &S_c,
                                    Array &S_i, Array &S_r, Array &S_s, Array &S_c_c ){
// determination of temperature and pressure by the law of Clausius-Clapeyron for water vapour concentration
// reaching saturation of water vapour pressure leads to formation of rain or ice
// precipitation and cloud formation by formulas from Häckel
// dry adiabatic lapse rate and saturated adiabatic lapse rate = temperature decrease with hight
// SL stands for sea level

// calculation of a total quantity as sum on all values in a virtual column in r-direction
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            if ( RadiationModel >= 1 ) Q_radiation.y[ j ][ k ] = radiation_3D.x[ 0 ][ j ][ k ];         // two- and multi-layer radiation balance assumed

            precipitable_water.y[ j ][ k ] = 0.;                                // precipitable water
            Evaporation_Penman.y[ j ][ k ] = 0.;                                // Evaporation by Penman
            Evaporation_Dalton.y[ j ][ k ] = 0.;                                // Evaporation by Dalton
        }
    }

    double zeta = 3.715;

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im-1; i++ ){
// on the boundary between land and air searching for the top of mountains
                if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i+1, j, k ) ) ){
                    t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;

                    r_dry = 100. * p_stat.x[ i ][ j ][ k ] / ( R_Air * t.x[ i ][ j ][ k ] * t_0 );
                    r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );
                    // density of humid air, COSMO version without cloud and ice water, masses negligible

                    e = c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep;  // water vapour pressure in hPa

                    t_denom = t_Celsius + 234.175;
                    E = hp * exp ( 17.0809 * t_Celsius / t_denom );
                    // saturation vapour pressure in the water phase for t > 0°C in hPa

                    Delta = 4000. * E / ( t_denom * t_denom );
                    // gradient of the water vapour pressure curve in hPa/K, coef = 234.175 * 17.0809
                    sat_deficit = ( E - e );  // saturation deficit in hPa
                    gamma = p_stat.x[ 0 ][ j ][ k ] * cp_l / ( ep * lv );  // Psychrometer constant in hPa/K

                    E_a = .35 * ( 1. + .15 * sqrt ( ( v.x[ i + 1 ][ j ][ k ] * v.x[ i + 1 ][ j ][ k ] +
                        w.x[ i + 1 ][ j ][ k ] * w.x[ i + 1 ][ j ][ k ] ) / 2. ) * u_0 * 3.6 ) * sat_deficit;
                        // ventilation-humidity Penmans formula

                    if ( is_land ( h, i, j, k ) )  Q_Evaporation.y[ j ][ k ] = 2300.;  // minimum value used for printout

                    Q_latent.y[ j ][ k ] = Q_Latent.x[ i ][ j ][ k ];  // latente heat in [W/m2] from energy transport equation
                    Q_sensible.y[ j ][ k ] = Q_Sensible.x[ i ][ j ][ k ];  // sensible heat in [W/m2] from energy transport equation
                    Q_bottom.y[ j ][ k ] = - ( Q_radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] -
                        Q_sensible.y[ j ][ k ] );    // difference understood as heat into the ground

                    if ( t_Celsius >= 0. ){
/*
                        double vau = .5 * ( v.x[ i ][ j ][ k ] + v.x[ i + 1 ][ j ][ k ] );
                        double weh = .5 * ( w.x[ i ][ j ][ k ] + w.x[ i + 1 ][ j ][ k ] );
                        Evaporation_Dalton.y[ j ][ k ] = C_Dalton ( u_0, vau, weh ) * sat_deficit * 24.;
*/
                        Evaporation_Dalton.y[ j ][ k ] = C_Dalton ( u_0, v.x[ i + 1 ][ j ][ k ],
                            w.x[ i + 1 ][ j ][ k ] ) * sat_deficit * 24.;

//                        Evaporation_Dalton.y[ j ][ k ] = C_Dalton ( u_0, v.x[ i ][ j ][ k ],
//                            w.x[ i ][ j ][ k ] ) * sat_deficit * 24.;  // mm/h in mm/d

//                        Evaporation_Dalton.y[ j ][ k ] = 8.46e-4 * C_Dalton ( u_0, v.x[ i+1 ][ j ][ k ], w.x[ i+1 ][ j ][ k ] ) *
//                        sat_deficit * dt_dim / ( r_humid * dr_dim ) * 24.;  // mm/h in mm/d
//                        Evaporation_Dalton.y[ j ][ k ] = 8.46e-4 * C_Dalton ( u_0, v.x[ i+1 ][ j ][ k ], w.x[ i+1 ][ j ][ k ] ) *
//                        sat_deficit * dt_dim / ( r_humid * dr_dim ) * 24.;  // mm/h in mm/d
                        // simplified formula for Evaporation by Dalton law dependent on surface water velocity in kg/(m²*d) = mm/d
                        if ( Evaporation_Dalton.y[ j ][ k ] <= 0. )  Evaporation_Dalton.y[ j ][ k ] = 0.;

//    cout << "   i = " << i << "   j = " << j << "   k = " << k << "   Evaporation_Dalton = " << Evaporation_Dalton.y[ j ][ k ] << "   c = " << c.x[ i ][ j ][ k ] << "   v = " << v.x[ i ][ j ][ k ] << "   w = " << w.x[ i ][ j ][ k ] << "   v1 = " << v.x[ i+1 ][ j ][ k ] << "   w1 = " << w.x[ i+1 ][ j ][ k ] << "   p = " << p_stat.x[ i ][ j ][ k ] << "   p1 = " << p_stat.x[ i+1 ][ j ][ k ] << "   sat_deficit = " << sat_deficit << "   t_Celsius = " << t_Celsius << endl;
                    }else   Evaporation_Dalton.y[ j ][ k ] = 0.;


                    Evaporation_Penman.y[ j ][ k ] = f_Penman * .0346 * ( ( Q_radiation.y[ j ][ k ] +
                        Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma );
                        // .0346 coefficient W/m2 corresponds to mm/d (Kraus)
                    if ( Evaporation_Penman.y[ j ][ k ] <= 0. )  Evaporation_Penman.y[ j ][ k ] = 0.;
                }

// only on the sea surface
                if ( ( i == 0 ) && ( is_air ( h, 0, j, k ) ) ){
                    t_Celsius = t.x[ 0 ][ j ][ k ] * t_0 - t_0;
                    r_dry = 100. * p_stat.x[ 0 ][ j ][ k ] / ( R_Air * t.x[ 0 ][ j ][ k ] * t_0 );
                    r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );
                    // density of humid air, COSMO version withot cloud and ice water, masses negligible
                    e = c.x[ 0 ][ j ][ k ] * p_stat.x[ 0 ][ j ][ k ] / ep;  // water vapour pressure in hPa
                    t_denom = t_Celsius + 234.175;
                    E = hp * exp ( 17.0809 * t_Celsius / t_denom );  // saturation vapour pressure in the water phase for t > 0°C in hPa
                    Delta = 4000. * E / ( t_denom * t_denom );  // gradient of the water vapour pressure curve in hPa/K, coef = 234.175 * 17.0809
                    sat_deficit = ( E - e );  // saturation deficit in hPa/K
                    gamma = p_stat.x[ 0 ][ j ][ k ] * cp_l / ( ep * lv );  // Psychrometer constant in hPa/K
                    E_a = .35 * ( 1. + .15 * sqrt ( ( v.x[ 1 ][ j ][ k ] * v.x[ 1 ][ j ][ k ] +
                        w.x[ 1 ][ j ][ k ] * w.x[ 1 ][ j ][ k ] ) / 2. ) * u_0 * 3.6 ) * sat_deficit;  // ventilation-humidity Penmans formula

                    if ( t_Celsius >= - 2. )  Q_Evaporation.y[ j ][ k ] = ( 2500.8 - 2.372 *
                        ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) );    // heat of Evaporation of water in [kJ/kg] (Kuttler) => variable lv
                    else  Q_Evaporation.y[ j ][ k ] = ( 2500.8 - 2.372 *
                        ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) + 300.; // heat of Evaporation of ice + 300 [kJ/kg]

                    Q_latent.y[ j ][ k ] = Q_Latent.x[ 0 ][ j ][ k ];  // latente heat in [W/m2] from energy transport equation
                    Q_sensible.y[ j ][ k ] = Q_Sensible.x[ 0 ][ j ][ k ];  // sensible heat in [W/m2] from energy transport equation
                    Q_bottom.y[ j ][ k ] = - ( Q_radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] -
                        Q_sensible.y[ j ][ k ] );  // difference understood as heat of the ground

                    if ( t_Celsius >= 0. ){
                        Evaporation_Dalton.y[ j ][ k ] = C_Dalton ( u_0, v.x[ 0 ][ j ][ k ], w.x[ 0 ][ j ][ k ] ) *
                            sat_deficit * 24.;  // ocean surface evaporation, mm/h in mm/d
//                        Evaporation_Dalton.y[ j ][ k ] = C_Dalton ( u_0, v.x[ 0 ][ j ][ k ], w.x[ 0 ][ j ][ k ] ) *
//                        sat_deficit * dt_dim / ( r_humid * dr_dim ) * 24.;  // mm/h in mm/d
//                        Evaporation_Dalton.y[ j ][ k ] = 8.46e-4 * C_Dalton ( u_0, v.x[ 0 ][ j ][ k ], w.x[ 0 ][ j ][ k ] ) *
//                        sat_deficit * dt_dim / ( r_humid * dr_dim ) * 24.;  // mm/h in mm/d
                      // simplified formula for Evaporation by Dalton law dependent on surface water velocity in kg/(m²*d) = mm/d
                      // not air but ocean surface temperature
                      // should be involved in water vapour saturation difference, it is not the saturation deficit
                        if ( Evaporation_Dalton.y[ j ][ k ] <= 0. )  Evaporation_Dalton.y[ j ][ k ] = 0.;
                    }else   Evaporation_Dalton.y[ j ][ k ] = 0.;


                    Evaporation_Penman.y[ j ][ k ] = f_Penman * .0346 * ( ( Q_radiation.y[ j ][ k ] +
                        Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma );
                        // .0346 coefficient W/m2 corresponds to mm/d (Kraus)
                    if ( Evaporation_Penman.y[ j ][ k ] <= 0. )  Evaporation_Penman.y[ j ][ k ] = 0.;
                }
            }
        }
    }



// boundaries of various variables
    for ( int k = 1; k < km-1; k++ ){
        for ( int j = 1; j < jm-1; j++ ){
            BuoyancyForce.x[ 0 ][ j ][ k ] = c43 * BuoyancyForce.x[ 1 ][ j ][ k ] -
                c13 * BuoyancyForce.x[ 2 ][ j ][ k ];
            BuoyancyForce.x[ im-1 ][ j ][ k ] = c43 * BuoyancyForce.x[ im-2 ][ j ][ k ] -
                c13 * BuoyancyForce.x[ im-3 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) )  BuoyancyForce.x[ 0 ][ j ][ k ] = 0.;
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int i = 0; i < im; i++ ){
            BuoyancyForce.x[ i ][ 0 ][ k ] = c43 * BuoyancyForce.x[ i ][ 1 ][ k ] -
                 c13 * BuoyancyForce.x[ i ][ 2 ][ k ];
            BuoyancyForce.x[ i ][ jm-1 ][ k ] = c43 * BuoyancyForce.x[ i ][ jm-2 ][ k ] -
                 c13 * BuoyancyForce.x[ i ][ jm-3 ][ k ];
//            if ( h.x[ i ][ 0 ][ k ] == 1. )     BuoyancyForce.x[ i ][ 0 ][ k ] = 0.;
        }
    }

    for ( int i = 0; i < im; i++ ){
        for ( int j = 1; j < jm-1; j++ ){
            BuoyancyForce.x[ i ][ j ][ 0 ] = c43 * BuoyancyForce.x[ i ][ j ][ 1 ] -
                c13 * BuoyancyForce.x[ i ][ j ][ 2 ];
            BuoyancyForce.x[ i ][ j ][ km-1 ] = c43 * BuoyancyForce.x[ i ][ j ][ km-2 ] -
                c13 * BuoyancyForce.x[ i ][ j ][ km-3 ];
            BuoyancyForce.x[ i ][ j ][ 0 ] = BuoyancyForce.x[ i ][ j ][ km-1 ] =
                ( BuoyancyForce.x[ i ][ j ][ 0 ] + BuoyancyForce.x[ i ][ j ][ km-1 ] ) / 2.;
//            if ( h.x[ i ][ j ][ 0 ] == 1. )     BuoyancyForce.x[ i ][ j ][ 0 ] = 0.;
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
                if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
                if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
                if ( ice.x[ i ][ j ][ k ] < 0. )  ice.x[ i ][ j ][ k ] = 0.;
                if ( P_rain.x[ i ][ j ][ k ] < 0. )  P_rain.x[ i ][ j ][ k ] = 0.;
                if ( P_snow.x[ i ][ j ][ k ] < 0. )  P_snow.x[ i ][ j ][ k ] = 0.;
            }
        }
    }

    double step = 0.;

// surface values of Evaporation, Condensation, Water, Water_super, IceAir, precipitable_water only for radial printout
    precipitable_water.y[ 0 ][ 0 ] = 0.;
    precipitable_water.y[ 0 ][ 0 ] = 0.;
    co2_total.y[ 0 ][ 0 ] = 0.;

    precipitation_NASA_average = 0.;
    precipitablewater_average = 0.;
    precipitation_average = 0.;
    Evaporation_Penman_average = 0.;
    Evaporation_Dalton_average = 0.;
    co2_vegetation_average = 0.;
    temperature_surf_average = 0.;
    velocity_average = 0.;
    temperature_average = 0.;
    latent_heat_average = 0.;
    sensible_heat_average = 0.;
    lat_av_loc = 0.;
    sen_av_loc = 0.;

    for ( int i = 0; i < im; i++ ){
        vel_av[ i ] = 0.;
        temp_av[ i ] = 0.;
        lat_av[ i ] = 0.;
        sen_av[ i ] = 0.;
    }

// printout of surface data at one predefinded location
    int j_loc_Dresden = 39;
    int k_loc_Dresden = 346;

    int j_loc_Sydney = 123;
    int k_loc_Sydney = 151;

    int j_loc_Pacific = 90;
    int k_loc_Pacific = 180;

    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            co2_total.y[ j ][ k ] = co2.x[ 0 ][ j ][ k ];
            for ( int i = 0; i < im; i++ ){
                e = 100. * c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep;  // water vapour pressure in Pa
                a = e / ( R_WaterVapour * t.x[ i ][ j ][ k ] * t_0 );  // absolute humidity in kg/m³
//                precipitable_water.y[ j ][ k ] +=  a * L_atm / ( double ) ( im - 1 );  // mass of water in kg/m²

                step = ( ( exp( zeta * ( rad.z[ i + 1 ] - 1. ) ) - 1 ) 
                     - ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) ) 
                     * ( L_atm / ( double ) ( im-1 ) );
                precipitable_water.y[ j ][ k ] +=  a * step;
                    // mass of water in kg/m², coordinate stretching 
                    // precipitable_water mass in 1 kg/m² compares to 1 mm height, with water density kg/ ( m² * mm )
            }
        }
    }

    double coeff_prec = 8.64e+4;  // dimensions see below
//    double coeff_rain_snow = dr * ( L_atm / ( double ) ( im-1 ) );  // dimensionsional in kg/(m²*s) == mm/s

// surface values of precipitation and precipitable water
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
//            P_rain.x[ 0 ][ j ][ k ] = coeff_rain_snow * P_rain.x[ 0 ][ j ][ k ];
//            P_snow.x[ 0 ][ j ][ k ] = coeff_rain_snow * P_snow.x[ 0 ][ j ][ k ];
//            Precipitation.y[ j ][ k ] = coeff_prec * coeff_rain_snow * ( P_rain.x[ 0 ][ j ][ k ] + P_snow.x[ 0 ][ j ][ k ] );
            Precipitation.y[ j ][ k ] = coeff_prec * ( P_rain.x[ 0 ][ j ][ k ] + P_snow.x[ 0 ][ j ][ k ] );
            // 60 s * 60 min * 24 h = 86400 s == 1 d
            // Precipitation, P_rain and P_snow in kg/ ( m² * s ) = mm/s
            // Precipitation in 86400. * kg/ ( m² * d ) = 86400 mm/d
            // kg/ ( m² * s ) == mm/s ( Kraus, p. 94 )
            if ( Precipitation.y[ j ][ k ] >= 25. )  Precipitation.y[ j ][ k ] = 25.;
            if ( Precipitation.y[ j ][ k ] <= 0 )  Precipitation.y[ j ][ k ] = 0.;
        }
    }

// averaging of precipitation in mm/a and precipitable water in mm
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            precipitation_NASA_average += precipitation_NASA.y[ j ][ k ];
            precipitablewater_average += precipitable_water.y[ j ][ k ];
            precipitation_average += Precipitation.y[ j ][ k ];
            temperature_surf_average += t.x[ 0 ][ j ][ k ] * t_0 - t_0;
            Evaporation_Penman_average += Evaporation_Penman.y[ j ][ k ];
            Evaporation_Dalton_average += Evaporation_Dalton.y[ j ][ k ];
            co2_vegetation_average += co2_total.y[ j ][ k ];
        }
    }

    temperature_surf_average = ( GetMean_3D(jm, km, t) - 1. ) * t_0;
    precipitablewater_average = GetMean_2D (jm, km, precipitable_water);
    precipitation_average = 365. * GetMean_2D(jm, km, Precipitation);
    precipitation_NASA_average = 365. * GetMean_2D(jm, km, precipitation_NASA);
    co2_vegetation_average = GetMean_2D(jm, km, co2_total);
    Evaporation_Penman_average = 365. * GetMean_2D (jm, km, Evaporation_Penman);
    Evaporation_Dalton_average = 365. * GetMean_2D (jm, km, Evaporation_Dalton);




/*
    co2_vegetation_average = co2_vegetation_average / ( double ) ( ( jm - 1 ) * ( km - 1 ) );
    precipitation_NASA_average = precipitation_NASA_average  / ( double ) ( ( jm - 1 ) * ( km - 1 ) );
    precipitablewater_average = precipitablewater_average  / ( double ) ( ( jm - 1 ) * ( km - 1 ) );
    precipitation_average = precipitation_average  / ( double ) ( ( jm - 1 ) * ( km - 1 ) );
    Evaporation_Penman_average = Evaporation_Penman_average  / ( double ) ( ( jm - 1 ) * ( km - 1 ) );
    Evaporation_Dalton_average = Evaporation_Dalton_average  / ( double ) ( ( jm - 1 ) * ( km - 1 ) );
*/


/*
    precipitablewater_average = GetMean_2D (jm, km, precipitable_water);
    precipitation_average = GetMean_2D(jm, km, Precipitation) / coeff_mmWS * 86400.;  // 60 s * 60 min * 24 h = 86400 s == 1 d
//    precipitation_NASA_average = GetMean_2D(jm, km, precipitation_NASA) / coeff_mmWS * 86400.;  // 60 s * 60 min * 24 h = 86400 s == 1 d;
    precipitation_NASA_average = GetMean_2D(jm, km, precipitation_NASA);
*/

    cout.precision ( 2 );

    level = "m";
    deg_north = "°N";
    deg_south = "°S";
    deg_west = "°W";
    deg_east = "°E";

    name_Value_1 = " radiation emission";
    name_Value_2 = " latent heat ";
    name_Value_3 = " sensible heat ";
    name_Value_4 = " bottom heat ";
    name_Value_5 = " Evaporation Penman ";
    name_Value_6 = " Evaporation Dalton ";
    name_Value_7 = " precipitable water average ";
    name_Value_8 = " precipitation average per year ";
    name_Value_9 = " precipitation average per day ";
    name_Value_10 = " precipitation NASA average per year ";
    name_Value_11 = " precipitation NASA average per day ";
    name_Value_12 = " Evaporation_Penman_average per year ";
    name_Value_13 = " Evaporation_Penman_average per day ";
    name_Value_14 = " Evaporation_Dalton_average per year ";
    name_Value_15 = " Evaporation_Dalton_average per day ";
    name_Value_16 = " to fill ";
    name_Value_17 = " to fill ";
    name_Value_18 = " to fill ";
    name_Value_19 = " to fill ";
    name_Value_20 = " to fill ";
    name_Value_21 = " to fill ";
    name_Value_22 = " co2_average ";
    name_Value_23 = " precipitable water ";
    name_Value_24 = " precipitation ";
    name_Value_25 = " temperature_modern_average ";
    name_Value_26 = " temperature_surf_average ";
    name_Value_27 = " temperature_expected_average ";

    name_unit_wm2 = " W/m2";
    name_unit_mmd = " mm/d";
    name_unit_mm = " mm";
    name_unit_mma = " mm/a";
    name_unit_ppm = " ppm";
    name_unit_t = " C";

    heading = " printout of surface data at predefinded locations: level, latitude, longitude";
    heading_Dresden = " City of Dresden, Germany, Europe";
    heading_Sydney = " City of Sydney, New South Wales, Australia";
    heading_Pacific = " Equator in the central Pacific";

    cout << endl << endl << heading << endl << endl;

    int choice = { 1 };
    preparation:
    switch ( choice ){
        case 1 :    cout << heading_Dresden << endl;
                        i_loc_level = 0;  // sea level
                        j_loc = j_loc_Dresden;  // 51°N, Dresden Germany
                        k_loc = k_loc_Dresden;  // 14°W, Dresden Germany
                        break;
        case 2 :    cout << heading_Sydney << endl;
                        i_loc_level = 0;  // sea level
                        j_loc = j_loc_Sydney;  // 33°S, Dresden Germany
                        k_loc = k_loc_Sydney;  // 151°E, Dresden Germany
                        break;
        case 3 :    cout << heading_Pacific << endl;
                        i_loc_level = 0;  // sea level
                        j_loc = j_loc_Pacific;  // 0°N, Equator
                        k_loc = k_loc_Pacific;  // 180°E, central Pacific
                        break;
    default :     cout << choice << "error in iterationPrintout member function in class Accuracy" << endl;
    }

    if ( j_loc <= 90 ){
        j_loc_deg = 90 - j_loc;
        deg_lat = deg_north;
    }

    if ( j_loc > 90 ){
        j_loc_deg = j_loc - 90;
        deg_lat = deg_south;
    }

    if ( k_loc <= 180 ){
        k_loc_deg = k_loc;
        deg_lon = deg_east;
    }

    if ( k_loc > 180 ){
        k_loc_deg = 360 - k_loc;
        deg_lon = deg_east;
    }

    Value_1 = Q_radiation.y[ j_loc ][ k_loc ];
    Value_2 = Q_latent.y[ j_loc ][ k_loc ];
    Value_3 = Q_sensible.y[ j_loc ][ k_loc ];
    Value_4 = Q_bottom.y[ j_loc ][ k_loc ];

    cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg
        << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon
        << "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_1
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_1 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left )
        << setw ( 25 ) << setfill ( '.' ) << name_Value_2 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_2 << setw ( 6 ) << name_unit_wm2
        << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_3
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_3 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left )
        << setw ( 25 ) << setfill ( '.' ) << name_Value_4 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_4 << setw ( 6 )
        << name_unit_wm2 << endl;

    Value_17 = 0.;
    Value_18 = 0.;
    Value_19 = 0.;
    Value_23 = precipitable_water.y[ j_loc ][ k_loc ];

    cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg
        << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon
        << "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_23
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_23 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left )
        << setw ( 25 ) << setfill ( '.' ) << name_Value_19 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_17 << setw ( 6 ) << name_unit_wm2
        << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_20
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_18 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left )
        << setw ( 25 ) << setfill ( '.' ) << name_Value_21 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_19 << setw ( 6 )
        << name_unit_wm2 << endl;

    Value_14 = 0.;
    Value_15 = 0.;
    Value_16 = 0.;
    Value_24 = Precipitation.y[ j_loc ][ k_loc ];

    cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg
        << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon
        << "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_24
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_24 << setw ( 6 ) << name_unit_mmd << "   " << setiosflags ( ios::left )
        << setw ( 25 ) << setfill ( '.' ) << name_Value_16 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_14 << setw ( 6 ) << name_unit_wm2
        << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_17
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_15 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left )
        << setw ( 25 ) << setfill ( '.' ) << name_Value_18 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_16 << setw ( 6 )
        << name_unit_wm2 << endl;

    Value_5 = Evaporation_Penman.y[ j_loc ][ k_loc ];
    Value_6 = Evaporation_Dalton.y[ j_loc ][ k_loc ];

    cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg
        << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon
        << "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_5
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_5 << setw ( 6 ) << name_unit_mmd << "   " << setiosflags ( ios::left )
        << setw ( 25 ) << setfill ( '.' ) << name_Value_6 << " = " << resetiosflags ( ios::left )
        << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_6 << setw ( 6 )
        << name_unit_mmd << endl << endl;

    choice++;
    if ( choice <= 3 ) goto preparation;

    cout << endl;

    Value_7 = precipitablewater_average;
    Value_8 = precipitation_average;

    cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' )
        << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed
        << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left )
        << setw ( 40 ) << setfill ( '.' ) << name_Value_8 << " = " << resetiosflags ( ios::left )
       << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_8 << setw ( 6 ) << name_unit_mma
        << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_9
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_8 / 365. << setw ( 6 ) << name_unit_mmd << endl;

    Value_10 = precipitation_NASA_average;

    cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' )
        << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 )
        << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm
        << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_10
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_10 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left )
        << setw ( 40 ) << setfill ( '.' ) << name_Value_11 << " = "
        << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_10 / 365. << setw ( 6 ) << name_unit_mmd << endl;

    Value_13 = Evaporation_Dalton_average;

    cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' )
        << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 )
        << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm
        << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' )
        << name_Value_14 << " = " << resetiosflags ( ios::left ) << setw ( 7 )
        << fixed << setfill ( ' ' ) << Value_13 << setw ( 6 ) << name_unit_mma
        << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' )
        << name_Value_15 << " = " << resetiosflags ( ios::left ) << setw ( 7 )
        << fixed << setfill ( ' ' ) << Value_13 / 365. << setw ( 6 )
        << name_unit_mmd << endl;

    Value_9 = co2_vegetation_average * co2_0;
    Value_12 = Evaporation_Penman_average;

    cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' )
        << name_Value_22 << " = " << resetiosflags ( ios::left ) << setw ( 7 )
        << fixed << setfill ( ' ' ) << Value_9 << setw ( 6 ) << name_unit_ppm
        << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_12
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_12 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left )
        << setw ( 40 ) << setfill ( '.' ) << name_Value_13 << " = "
        << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_12 / 365. << setw ( 6 ) << name_unit_mmd
        << endl << endl << endl;

    Value_25 = ( GetMean_2D(jm, km, temperature_NASA) - 1. ) * t_0;
    Value_26 = temperature_surf_average;
    Value_27 = t_average + t_cretaceous * t_0;

    cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' )
        << name_Value_25 << " = " << resetiosflags ( ios::left ) << setw ( 7 )
        << fixed << setfill ( ' ' ) << Value_25 << setw ( 6 ) << name_unit_t << "   "
        << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_26
        << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_26 << setw ( 6 ) << name_unit_t << "   " << setiosflags ( ios::left )
        << setw ( 40 ) << setfill ( '.' ) << name_Value_27 << " = "
        << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' )
        << Value_27 << setw ( 6 ) << name_unit_t << endl << endl << endl;
}




float Results_MSL_Atm::GetMean_3D(int jm, int km, Array &val_3D){
    if(m_node_weights.size() != (unsigned)jm){
        CalculateNodeWeights(jm, km);
    }
    double ret=0., weight=0.;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (val_3D.x[0][j][k]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret+=val_3D.x[0][j][k]*m_node_weights[j][k];
            weight+=m_node_weights[j][k];
        }
    }
    return ret/weight;
}




float Results_MSL_Atm::GetMean_2D(int jm, int km, Array_2D &val_2D){
    if(m_node_weights.size() != (unsigned)jm){
        CalculateNodeWeights(jm, km);
    }
    double ret=0., weight=0.;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (val_2D.y[ j ][ k ]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret+=val_2D.y[ j ][ k ]*m_node_weights[j][k];
            weight+=m_node_weights[j][k];
        }
    }
    return ret/weight;
}




void Results_MSL_Atm::CalculateNodeWeights(int jm, int km){
    //use cosine of latitude as weights for now
    //longitudes: 0-360(km) latitudes: 90-(-90)(jm)
    double weight = 0.;
    m_node_weights.clear();
    for(int i=0; i<jm; i++){
        if(i<=90){
            weight = cos((90-i) * M_PI / 180.0 );
        }else{
            weight = cos((i-90) * M_PI / 180.0 );
        }
        m_node_weights.push_back(std::vector<double>());
        m_node_weights[i].resize(km, weight);
    }
    return;
}



double Results_MSL_Atm::C_Dalton ( double u_0, double v, double w ){
    // variation of the heat transfer coefficient in Dalton's evaporation law, parabola
    // air velocity measured 2m above sea level
    double C_max = - .053;  // for v_max = 10 m/s, but C is function of v, should be included
    // Geiger ( 1961 ) by > Zmarsly, Kuttler, Pethe in mm/( h * hPa ), p. 133
    double v_max = 10.;
    double fac = 10.;  // factor to adjust the ratio of NASA precipitation 
                       // to Dalton evaporation for the modern world, 
                       // relative difference between global precipitation and evaporation is 10%
    double vel_magnitude = sqrt ( v * v + w * w ) * u_0;
    return sqrt ( C_max * C_max / v_max * vel_magnitude ) / fac;  // result in mm/h
}

