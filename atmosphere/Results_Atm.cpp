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

#include <cAtmosphereModel.h>
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

void cAtmosphereModel::run_MSL_data(){
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            if(RadiationModel >= 1) 
                Q_radiation.y[ j ][ k ] = radiation_3D.x[ 0 ][ j ][ k ];
            precipitable_water.y[ j ][ k ] = 0.;
            Evaporation_Dalton.y[ j ][ k ] = 0.;
        }
    }
    float c43 = 4./3.;
    float c13 = 1./3.;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            int i = get_surface_layer(j,k);
            float t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;
            float e = c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep;  // water vapour pressure in hPa
            float t_denom = t_Celsius + 234.175;
            float E = hp * exp ( 17.0809 * t_Celsius / t_denom );
                // saturation vapour pressure in the water phase for t > 0°C in hPa
            float sat_deficit = ( E - e );  // saturation deficit in hPa
            if(is_land( h, i, j, k)){ //ocean surface
                if(t_Celsius >= - 2.)  Q_Evaporation.y[ j ][ k ] = 
                    ( 2500.8 - 2.372 * ( t.x[ i ][ j ][ k ] * t_0 - t_0 ) );
                // heat of Evaporation of water in [kJ/kg] (Kuttler) => variable lv
                else  Q_Evaporation.y[ j ][ k ] = ( 2500.8 - 2.372 * 
                          ( t.x[ i ][ j ][ k ] * t_0 - t_0 ) ) + 300.;
                // heat of Evaporation of ice + 300 [kJ/kg]
            }else  Q_Evaporation.y[ j ][ k ] = 2300.;  // minimum value used for printout
            Q_latent.y[ j ][ k ] = Q_Latent.x[ i ][ j ][ k ];  // latente heat in [W/m2] from energy transport equation
            Q_sensible.y[ j ][ k ] = Q_Sensible.x[ i ][ j ][ k ];  // sensible heat in [W/m2] from energy transport equation
            Q_bottom.y[ j ][ k ] = - Q_radiation.y[ j ][ k ] 
                - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ];    
                // difference understood as heat into the ground
            if(t_Celsius >= 0.){
                Evaporation_Dalton.y[ j ][ k ] = C_Dalton(u_0, v.x[ i ][ j ][ k ],
                    w.x[ i ][ j ][ k ]) * sat_deficit * 24.;
                if(is_land( h, i, j, k)){
                    t_Celsius = t.x[ i+1 ][ j ][ k ] * t_0 - t_0;
                    e = c.x[ i+1 ][ j ][ k ] * p_stat.x[ i+1 ][ j ][ k ] / ep;  // water vapour pressure in hPa
                    t_denom = t_Celsius + 234.175;
                    E = hp * exp ( 17.0809 * t_Celsius / t_denom );
                    sat_deficit = ( E - e );  // saturation deficit in hPa
                    Evaporation_Dalton.y[ j ][ k ] = 
                        C_Dalton(u_0, v.x[ i+1 ][ j ][ k ], w.x[ i+1 ][ j ][ k ]) 
                        * sat_deficit * 24.;
                }
                // simplified formula for Evaporation by Dalton law dependent on surface water velocity in kg/(m²*d) = mm/d
                if(Evaporation_Dalton.y[ j ][ k ] <= 0.)  
                    Evaporation_Dalton.y[ j ][ k ] = 0.;
            }else  Evaporation_Dalton.y[ j ][ k ] = 0.;
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            BuoyancyForce.x[ im-1 ][ j ][ k ] = c43 * BuoyancyForce.x[ im-2 ][ j ][ k ] -
                c13 * BuoyancyForce.x[ im-3 ][ j ][ k ];
            if(is_land ( h, 0, j, k))  
                BuoyancyForce.x[ 0 ][ j ][ k ] = 0.;
            else
                BuoyancyForce.x[ 0 ][ j ][ k ] = c43 * BuoyancyForce.x[ 1 ][ j ][ k ] -
                    c13 * BuoyancyForce.x[ 2 ][ j ][ k ];
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            BuoyancyForce.x[ i ][ 0 ][ k ] = c43 * BuoyancyForce.x[ i ][ 1 ][ k ] -
                 c13 * BuoyancyForce.x[ i ][ 2 ][ k ];
            BuoyancyForce.x[ i ][ jm-1 ][ k ] = c43 * BuoyancyForce.x[ i ][ jm-2 ][ k ] -
                 c13 * BuoyancyForce.x[ i ][ jm-3 ][ k ];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 1; j < jm-1; j++){
            float b1 = c43 * BuoyancyForce.x[ i ][ j ][ 1 ] - c13 * BuoyancyForce.x[ i ][ j ][ 2 ];
            float b2 = c43 * BuoyancyForce.x[ i ][ j ][ km-2 ] - c13 * BuoyancyForce.x[ i ][ j ][ km-3 ];
            BuoyancyForce.x[ i ][ j ][ 0 ] = BuoyancyForce.x[ i ][ j ][ km-1 ] = (b1+ b2) / 2.;
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[ i ][ j ][ k ] < 0.)  c.x[ i ][ j ][ k ] = 0.;
                if(cloud.x[ i ][ j ][ k ] < 0.)  cloud.x[ i ][ j ][ k ] = 0.;
                if(ice.x[ i ][ j ][ k ] < 0.)  ice.x[ i ][ j ][ k ] = 0.;
                if(P_rain.x[ i ][ j ][ k ] < 0.)  P_rain.x[ i ][ j ][ k ] = 0.;
                if(P_snow.x[ i ][ j ][ k ] < 0.)  P_snow.x[ i ][ j ][ k ] = 0.;
            }
        }
    }
    // surface values of Evaporation, Condensation, Water, Water_super, IceAir, precipitable_water only for radial printout
    precipitable_water.y[ 0 ][ 0 ] = 0.;
    precipitable_water.y[ 0 ][ 0 ] = 0.;
    co2_total.y[ 0 ][ 0 ] = 0.;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            co2_total.y[ j ][ k ] = co2.x[ 0 ][ j ][ k ];
            for(int i = 0; i < im; i++){
                float e = 100. * c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep;  // water vapour pressure in Pa
                float a = e / ( R_WaterVapour * t.x[ i ][ j ][ k ] * t_0 );  // absolute humidity in kg/m³
                float step = get_layer_height(i+1) - get_layer_height(i);
                precipitable_water.y[ j ][ k ] +=  a * step;
                 // mass of water in kg/m²
                // precipitable_water mass in 1 kg/m² compares to 1 mm hight, with water density kg/ ( m² * mm )
            }
        }
    }
    double coeff_prec = 86400.;  // dimensions see below
    // surface values of precipitation and precipitable water
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Precipitation.y[ j ][ k ] = coeff_prec * ( P_rain.x[ 0 ][ j ][ k ] + P_snow.x[ 0 ][ j ][ k ] );
            // 60 s * 60 min * 24 h = 86400 s == 1 d
            // Precipitation, P_rain and P_snow in kg/ ( m² * s ) = mm/s
            // Precipitation in 86400. * kg/ ( m² * d ) = 86400 mm/d
            // kg/ ( m² * s ) == mm/s ( Kraus, p. 94 )
        }
    }
}


