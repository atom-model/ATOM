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

using namespace std;

Results_MSL_Atm::Results_MSL_Atm ( int im, int jm, int km, int sun, double g, double ep, double hp, double u_0, double p_0, double t_0, double c_0, double co2_0, double sigma, double albedo_equator, double lv, double ls, double cp_l, double L_atm, double dt, double dr, double dthe, double dphi, double r_air, double R_Air, double r_water, double r_water_vapour, double R_WaterVapour, double co2_vegetation, double co2_ocean, double co2_land, double gam, double t_pole, double t_cretaceous )
:	coeff_Diffusion_latent ( .005 ),													// diffusion coefficient for latent heat in [m²/s]    0.015
	coeff_Diffusion_sensibel ( .2 ),													// diffusion coefficient for sensible heat in [m²/s]     0.03
	f_Haude ( .3 )																				// Haude factor for evapotranspiration 0.3 for low, dense vegetation as raw average value by Kuttler
{
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
	this-> r_water = r_water;
	this-> co2_vegetation = co2_vegetation;
	this-> co2_ocean = co2_ocean;
	this-> co2_land = co2_land;
	this-> t_pole = t_pole;
	this-> t_cretaceous = t_cretaceous;

	coeff_mmWS = r_air / r_water_vapour;									// coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]
	coeff_lv = lv / ( cp_l * t_0 );														// coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_lv = 9.1069 in [ / ]
	coeff_ls = ls / ( cp_l * t_0 );														// coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_ls = 10.3091 in [ / ]

	c43 = 4./3.;
	c13 = 1./3.;
}


Results_MSL_Atm::~Results_MSL_Atm (){}



void Results_MSL_Atm::run_MSL_data ( int n, int velocity_iter_max, int RadiationModel, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &c, Array &cn, Array &co2, Array &co2n, Array &t, Array &tn, Array &p_dyn, Array &p_stat, Array &BuoyancyForce, Array &u, Array &v, Array &w, Array &Latency, Array &Q_Sensible, Array &radiation_3D, Array &t_cond_3D, Array &t_evap_3D, Array &cloud, Array &cloudn, Array &ice, Array &icen, Array &P_rain, Array &P_snow, Array &aux_u, Array &aux_v, Array &aux_w, Array_2D &precipitation_NASA, Array_2D &Evaporation, Array_2D &Condensation, Array_2D &LatentHeat, Array_2D &precipitable_water, Array_2D &Q_Radiation, Array_2D &Q_Evaporation, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Q_bottom, Array_2D &Evaporation_Penman, Array_2D &Evaporation_Haude, Array_2D &Vegetation, Array_2D &Radiation_Balance, Array_2D &albedo, Array_2D &co2_total, Array_2D &Precipitation, Array &S_v, Array &S_c, Array &S_i, Array &S_r, Array &S_s, Array &S_c_c )
{
// determination of temperature and pressure by the law of Clausius-Clapeyron for water vapour concentration
// reaching saturation of water vapour pressure leads to formation of rain or ice
// precipitation and cloud formation by formulas from Häckel
// dry adiabatic lapse rate and saturated adiabatic lapse rate = temperature decrease with hight
// SL stands for sea level

// calculation of a total quantity as sum on all values in a virtual column in r-direction
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( RadiationModel >= 2 ) Q_Radiation.y[ j ][ k ] = radiation_3D.x[ 0 ][ j ][ k ];		 // two- and multi-layer radiation balance assumed

			Evaporation.y[ j ][ k ] = 0.;												// Evaporation
			Condensation.y[ j ][ k ] = 0.;											// Condensation
			LatentHeat.y[ j ][ k ] = 0.;												// Condensation
			precipitable_water.y[ j ][ k ] = 0.;									// precipitable water
			Evaporation_Penman.y[ j ][ k ] = 0.;								// Evaporation by Penman
			Evaporation_Haude.y[ j ][ k ] = 0.;								// Evaporation by Haude

		}
	}


	for ( int k = 1; k < km-1; k++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
// on the boundary between land and air searching for the top of mountains
				if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i + 1 ][ j ][ k ] == 0. ) )
				{
					if ( i == 0 ) 	p_stat.x[ 0 ][ j ][ k ] = ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01; // given in hPa
					else 	p_stat.x[ i ][ j ][ k ] = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t.x[ i ][ j ][ k ] * t_0 ) ) * p_stat.x[ 0 ][ j ][ k ];	// given in hPa
																																			// current air pressure, step size in 500 m, from a polytropic atmosphere in hPa

					if ( Latency.x[ i ][ j ][ k ] <= 0. ) 	t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0 + ( ( t_cond_3D.x[ i ][ j ][ k ] ) / t_0 );
					else  											t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0 + ( ( t_evap_3D.x[ i ][ j ][ k ] ) / t_0 );

					r_dry = 100. * p_stat.x[ i ][ j ][ k ] / ( R_Air * t.x[ i ][ j ][ k ] * t_0 );
					r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );				// density of humid air, COSMO version withot cloud and ice water, masses negligible

					e = c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep; 													// water vapour pressure in hPa
					a = 216.6 * e / ( t.x[ i ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3

					t_denom = t_Celsius + 234.175;
					E = hp * exp ( 17.0809 * t_Celsius / t_denom );													// saturation vapour pressure in the water phase for t > 0°C in hPa

					Delta = 4000. * E / ( t_denom * t_denom );															// gradient of the water vapour pressure curve in hPa/K, coef = 234.175 * 17.0809
					sat_deficit = ( E - e );																							// saturation deficit in hPa/K
					gamma = p_stat.x[ 0 ][ j ][ k ] * cp_l / ( ep * lv );												// Psychrometer constant in hPa/K

					c_grad = ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i + 1 ][ j ][ k ] - c.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );	// water vapour pressure gradient in g/(Kg m)
					t_grad = ( - 3. * t.x[ i ][ j ][ k ] + 4. * t.x[ i + 1 ][ j ][ k ] - t.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );	// temperature gradient in K/m

					E_a = .35 * ( 1. + .15 * sqrt ( v.x[ i + 1 ][ j ][ k ] * v.x[ i + 1 ][ j ][ k ] + w.x[ i + 1 ][ j ][ k ] * w.x[ i + 1 ][ j ][ k ] ) * u_0 ) * sat_deficit;	// ventilation-humidity Penmans formula

					Q_Evaporation.y[ j ][ k ] = ( 2500.8 - 2.372 * ( t.x[ i ][ j ][ k ] * t_0 - t_0 ) ) * 1.e-3;	// heat of Evaporation of water in [MJ/kg] (Kuttler) = variable lv
					Q_latent.y[ j ][ k ] = - a * lv * coeff_Diffusion_latent * c_grad / L_atm;					// latente heat in [W/m2] from energy transport equation
					Q_sensible.y[ j ][ k ] = - r_air * cp_l * coeff_Diffusion_sensibel * t_grad * t_0 / L_atm;		// sensible heat in [W/m2] from energy transport equation
					Q_bottom.y[ j ][ k ] = Q_Radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ];	// difference understood as heat into the ground

					Evaporation_Haude.y[ j ][ k ] = f_Haude * sat_deficit;											// simplified formula for Evaporation over day length of 12h by Haude, Häckel
					if ( Evaporation_Haude.y[ j ][ k ] <= 0. ) 		Evaporation_Haude.y[ j ][ k ] = 0.;
					Evaporation_Penman.y[ j ][ k ] = .0346 * ( ( Q_Radiation.y[ j ][ k ] - Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma );
					if ( Evaporation_Penman.y[ j ][ k ] <= 0. ) 		Evaporation_Penman.y[ j ][ k ] = 0.;	// .0346 coefficient W/m2 corresponds to mm/d (Kraus)

				}



// only on the sea surface
				if ( ( i == 0 ) && ( h.x[ 0 ][ j ][ k ] == 0. ) )
				{
					if ( i == 0 ) 	p_stat.x[ 0 ][ j ][ k ] = ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;		// given in hPa

					if ( Latency.x[ 0 ][ j ][ k ] <= 0. ) 	t_Celsius = t.x[ 0 ][ j ][ k ] * t_0 - t_0 + ( ( t_cond_3D.x[ 0 ][ j ][ k ] ) / t_0 );
					else  											t_Celsius = t.x[ 0 ][ j ][ k ] * t_0 - t_0 + ( ( t_evap_3D.x[ 0 ][ j ][ k ] ) / t_0 );

					r_dry = 100. * p_stat.x[ 0 ][ j ][ k ] / ( R_Air * t.x[ 0 ][ j ][ k ] * t_0 );
					r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );				// density of humid air, COSMO version withot cloud and ice water, masses negligible

					e = c.x[ 0 ][ j ][ k ] * p_stat.x[ 0 ][ j ][ k ] / ep; 													// water vapour pressure in hPa
					a = 216.6 * e / ( t.x[ 0 ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3

					t_denom = t_Celsius + 234.175;
					E = hp * exp ( 17.0809 * t_Celsius / t_denom );													// saturation vapour pressure in the water phase for t > 0°C in hPa
					Delta = 4000. * E / ( t_denom * t_denom );															// gradient of the water vapour pressure curve in hPa/K, coef = 234.175 * 17.0809

					sat_deficit = ( E - e );																							// saturation deficit in hPa/K
					gamma = p_stat.x[ 0 ][ j ][ k ] * cp_l / ( ep * lv );												// Psychrometer constant in hPa/K

					c_grad = ( - 3. * c.x[ 0 ][ j ][ k ] + 4. * c.x[ 1 ][ j ][ k ] - c.x[ 2 ][ j ][ k ] ) / ( 2. * dr );	// water vapour pressure gradient in g/(Kg m)
					t_grad = ( - 3. * t.x[ 0 ][ j ][ k ] + 4. * t.x[ 1 ][ j ][ k ] - t.x[ 2 ][ j ][ k ] ) / ( 2. * dr );	// temperature gradient in K/m

					E_a = .35 * ( 1. + .15 * sqrt ( v.x[ 1 ][ j ][ k ] * v.x[ 1 ][ j ][ k ] + w.x[ 1 ][ j ][ k ] * w.x[ 1 ][ j ][ k ] ) ) * sat_deficit;	// ventilation-humidity for Penman's formula

					Q_Evaporation.y[ j ][ k ] = ( 2500.8 - 2.372 * ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) * 1.e-3;	// heat of Evaporation of water in [MJ/kg] (Kuttler) = lv
					Q_latent.y[ j ][ k ] = - a * lv * coeff_Diffusion_latent * c_grad / L_atm;		// latente heat in [W/m2] from energy transport equation
					Q_sensible.y[ j ][ k ] = - r_air * cp_l * coeff_Diffusion_sensibel * t_grad * t_0 / L_atm;		// sensible heat in [W/m2] from energy transport equation
					Q_bottom.y[ j ][ k ] = Q_Radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ];	// difference understood as heat of the ground

					Evaporation_Haude.y[ j ][ k ] = 0.;
					if ( Evaporation_Haude.y[ j ][ k ] <= 0. ) 		Evaporation_Haude.y[ j ][ k ] = 0.;
//					Evaporation_Penman.y[ j ][ k ] = .0346 * ( ( Q_Radiation.y[ j ][ k ] - Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma );
					Evaporation_Penman.y[ j ][ k ] = 0.;
//					Evaporation_Penman.y[ j ][ k ] = .0346 * ( ( Q_Radiation.y[ j ][ k ] - Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma ) * 50.;
					if ( Evaporation_Penman.y[ j ][ k ] <= 0. ) 		Evaporation_Penman.y[ j ][ k ] = 0.;	// .0346 coefficient W/m2 corresponds to mm/d (Kraus)
				}
			}

		}
	}




// boundaries of various variables
	for ( int k = 1; k < km-1; k++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			Latency.x[ 0 ][ j ][ k ] = c43 * Latency.x[ 1 ][ j ][ k ] - c13 * Latency.x[ 2 ][ j ][ k ];
			Latency.x[ im-1 ][ j ][ k ] = c43 * Latency.x[ im-2 ][ j ][ k ] - c13 * Latency.x[ im-3 ][ j ][ k ];
			if ( h.x[ 0 ][ j ][ k ] == 1. ) 	Latency.x[ 0 ][ j ][ k ] = 0.;

			Q_Sensible.x[ 0 ][ j ][ k ] = c43 * Q_Sensible.x[ 1 ][ j ][ k ] - c13 * Q_Sensible.x[ 2 ][ j ][ k ];
			Q_Sensible.x[ im-1 ][ j ][ k ] = c43 * Q_Sensible.x[ im-2 ][ j ][ k ] - c13 * Q_Sensible.x[ im-3 ][ j ][ k ];
			if ( h.x[ 0 ][ j ][ k ] == 1. ) 	Q_Sensible.x[ 0 ][ j ][ k ] = 0.;

			t_cond_3D.x[ 0 ][ j ][ k ] = c43 * t_cond_3D.x[ 1 ][ j ][ k ] - c13 * t_cond_3D.x[ 2 ][ j ][ k ];
			t_cond_3D.x[ im-1 ][ j ][ k ] = c43 * t_cond_3D.x[ im-2 ][ j ][ k ] - c13 * t_cond_3D.x[ im-3 ][ j ][ k ];
			if ( h.x[ 0 ][ j ][ k ] == 1. ) 	t_cond_3D.x[ 0 ][ j ][ k ] = 0.;

			t_evap_3D.x[ 0 ][ j ][ k ] = c43 * t_evap_3D.x[ 1 ][ j ][ k ] - c13 * t_evap_3D.x[ 2 ][ j ][ k ];
			t_evap_3D.x[ im-1 ][ j ][ k ] = c43 * t_evap_3D.x[ im-2 ][ j ][ k ] - c13 * t_evap_3D.x[ im-3 ][ j ][ k ];
			if ( h.x[ 0 ][ j ][ k ] == 1. ) 	t_evap_3D.x[ 0 ][ j ][ k ] = 0.;

			BuoyancyForce.x[ 0 ][ j ][ k ] = c43 * BuoyancyForce.x[ 1 ][ j ][ k ] - c13 * BuoyancyForce.x[ 2 ][ j ][ k ];
			BuoyancyForce.x[ im-1 ][ j ][ k ] = c43 * BuoyancyForce.x[ im-2 ][ j ][ k ] - c13 * BuoyancyForce.x[ im-3 ][ j ][ k ];
			if ( h.x[ 0 ][ j ][ k ] == 1. ) 	BuoyancyForce.x[ 0 ][ j ][ k ] = 0.;
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int i = 0; i < im; i++ )
		{
			t.x[ i ][ 0 ][ k ] = c43 * t.x[ i ][ 1 ][ k ] - c13 * t.x[ i ][ 2 ][ k ];
			t.x[ i ][ jm-1 ][ k ] = c43 * t.x[ i ][ jm-2 ][ k ] - c13 * t.x[ i ][ jm-3 ][ k ];

			radiation_3D.x[ i ][ 0 ][ k ] = c43 * radiation_3D.x[ i ][ 1 ][ k ] - c13 * radiation_3D.x[ i ][ 2 ][ k ];
			radiation_3D.x[ i ][ jm-1 ][ k ] = c43 * radiation_3D.x[ i ][ jm-2 ][ k ] - c13 * radiation_3D.x[ i ][ jm-3 ][ k ];

			Latency.x[ i ][ 0 ][ k ] = c43 * Latency.x[ i ][ 1 ][ k ] - c13 * Latency.x[ i ][ 2 ][ k ];
			Latency.x[ i ][ jm-1 ][ k ] = c43 * Latency.x[ i ][ jm-2 ][ k ] - c13 * Latency.x[ i ][ jm-3 ][ k ];
			if ( h.x[ i ][ 0 ][ k ] == 1. ) 	Latency.x[ i ][ 0 ][ k ] = 0.;

			Q_Sensible.x[ i ][ 0 ][ k ] = c43 * Q_Sensible.x[ i ][ 1 ][ k ] - c13 * Q_Sensible.x[ i ][ 2 ][ k ];
			Q_Sensible.x[ i ][ jm-1 ][ k ] = c43 * Q_Sensible.x[ i ][ jm-2 ][ k ] - c13 * Q_Sensible.x[ i ][ jm-3 ][ k ];
			if ( h.x[ i ][ 0 ][ k ] == 1. ) 	Q_Sensible.x[ i ][ 0 ][ k ] = 0.;

			t_cond_3D.x[ i ][ 0 ][ k ] = c43 * t_cond_3D.x[ i ][ 1 ][ k ] - c13 * t_cond_3D.x[ i ][ 2 ][ k ];
			t_cond_3D.x[ i ][ jm-1 ][ k ] = c43 * t_cond_3D.x[ i ][ jm-2 ][ k ] - c13 * t_cond_3D.x[ i ][ jm-3 ][ k ];
			if ( h.x[ i ][ 0 ][ k ] == 1. ) 	t_cond_3D.x[ i ][ 0 ][ k ] = 0.;

			t_evap_3D.x[ i ][ 0 ][ k ] = c43 * t_evap_3D.x[ i ][ 1 ][ k ] - c13 * t_evap_3D.x[ i ][ 2 ][ k ];
			t_evap_3D.x[ i ][ jm-1 ][ k ] = c43 * t_evap_3D.x[ i ][ jm-2 ][ k ] - c13 * t_evap_3D.x[ i ][ jm-3 ][ k ];
			if ( h.x[ i ][ 0 ][ k ] == 1. ) 	t_evap_3D.x[ i ][ 0 ][ k ] = 0.;

			BuoyancyForce.x[ i ][ 0 ][ k ] = c43 * BuoyancyForce.x[ i ][ 1 ][ k ] - c13 * BuoyancyForce.x[ i ][ 2 ][ k ];
			BuoyancyForce.x[ i ][ jm-1 ][ k ] = c43 * BuoyancyForce.x[ i ][ jm-2 ][ k ] - c13 * BuoyancyForce.x[ i ][ jm-3 ][ k ];
			if ( h.x[ i ][ 0 ][ k ] == 1. ) 	BuoyancyForce.x[ i ][ 0 ][ k ] = 0.;


			c.x[ i ][ 0 ][ k ] = c43 * c.x[ i ][ 1 ][ k ] - c13 * c.x[ i ][ 2 ][ k ];
			c.x[ i ][ jm-1 ][ k ] = c43 * c.x[ i ][ jm-2 ][ k ] - c13 * c.x[ i ][ jm-3 ][ k ];

			cloud.x[ i ][ 0 ][ k ] = c43 * cloud.x[ i ][ 1 ][ k ] - c13 * cloud.x[ i ][ 2 ][ k ];
			cloud.x[ i ][ jm-1 ][ k ] = c43 * cloud.x[ i ][ jm-2 ][ k ] - c13 * cloud.x[ i ][ jm-3 ][ k ];

			ice.x[ i ][ 0 ][ k ] = c43 * ice.x[ i ][ 1 ][ k ] - c13 * ice.x[ i ][ 2 ][ k ];
			ice.x[ i ][ jm-1 ][ k ] = c43 * ice.x[ i ][ jm-2 ][ k ] - c13 * ice.x[ i ][ jm-3 ][ k ];

			P_rain.x[ i ][ 0 ][ k ] = c43 * P_rain.x[ i ][ 1 ][ k ] - c13 * P_rain.x[ i ][ 2 ][ k ];
			P_rain.x[ i ][ jm-1 ][ k ] = c43 * P_rain.x[ i ][ jm-2 ][ k ] - c13 * P_rain.x[ i ][ jm-3 ][ k ];

			P_snow.x[ i ][ 0 ][ k ] = c43 * P_snow.x[ i ][ 1 ][ k ] - c13 * P_snow.x[ i ][ 2 ][ k ];
			P_snow.x[ i ][ jm-1 ][ k ] = c43 * P_snow.x[ i ][ jm-2 ][ k ] - c13 * P_snow.x[ i ][ jm-3 ][ k ];
		}
	}


	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			t.x[ i ][ j ][ 0 ] = c43 * t.x[ i ][ j ][ 1 ] - c13 * t.x[ i ][ j ][ 2 ];
			t.x[ i ][ j ][ km-1 ] = c43 * t.x[ i ][ j ][ km-2 ] - c13 * t.x[ i ][ j ][ km-3 ];
			t.x[ i ][ j ][ 0 ] = t.x[ i ][ j ][ km-1 ] = ( t.x[ i ][ j ][ 0 ] + t.x[ i ][ j ][ km-1 ] ) / 2.;

			radiation_3D.x[ i ][ j ][ 0 ] = c43 * radiation_3D.x[ i ][ j ][ 1 ] - c13 * radiation_3D.x[ i ][ j ][ 2 ];
			radiation_3D.x[ i ][ j ][ km-1 ] = c43 * radiation_3D.x[ i ][ j ][ km-2 ] - c13 * radiation_3D.x[ i ][ j ][ km-3 ];
			radiation_3D.x[ i ][ j ][ 0 ] = radiation_3D.x[ i ][ j ][ km-1 ] = ( radiation_3D.x[ i ][ j ][ 0 ] + radiation_3D.x[ i ][ j ][ km-1 ] ) / 2.;

			Latency.x[ i ][ j ][ 0 ] = c43 * Latency.x[ i ][ j ][ 1 ] - c13 * Latency.x[ i ][ j ][ 2 ];
			Latency.x[ i ][ j ][ km-1 ] = c43 * Latency.x[ i ][ j ][ km-2 ] - c13 * Latency.x[ i ][ j ][ km-3 ];
			Latency.x[ i ][ j ][ 0 ] = Latency.x[ i ][ j ][ km-1 ] = ( Latency.x[ i ][ j ][ 0 ] + Latency.x[ i ][ j ][ km-1 ] ) / 2.;
			if ( h.x[ i ][ j ][ 0 ] == 1. ) 	Latency.x[ i ][ j ][ 0 ] = 0.;

			Q_Sensible.x[ i ][ j ][ 0 ] = c43 * Q_Sensible.x[ i ][ j ][ 1 ] - c13 * Q_Sensible.x[ i ][ j ][ 2 ];
			Q_Sensible.x[ i ][ j ][ km-1 ] = c43 * Q_Sensible.x[ i ][ j ][ km-2 ] - c13 * Q_Sensible.x[ i ][ j ][ km-3 ];
			Q_Sensible.x[ i ][ j ][ 0 ] = Q_Sensible.x[ i ][ j ][ km-1 ] = ( Q_Sensible.x[ i ][ j ][ 0 ] + Q_Sensible.x[ i ][ j ][ km-1 ] ) / 2.;
			if ( h.x[ i ][ j ][ 0 ] == 1. ) 	Q_Sensible.x[ i ][ j ][ 0 ] = 0.;

			t_cond_3D.x[ i ][ j ][ 0 ] = c43 * t_cond_3D.x[ i ][ j ][ 1 ] - c13 * t_cond_3D.x[ i ][ j ][ 2 ];
			t_cond_3D.x[ i ][ j ][ km-1 ] = c43 * t_cond_3D.x[ i ][ j ][ km-2 ] - c13 * t_cond_3D.x[ i ][ j ][ km-3 ];
			t_cond_3D.x[ i ][ j ][ 0 ] = t_cond_3D.x[ i ][ j ][ km-1 ] = ( t_cond_3D.x[ i ][ j ][ 0 ] + t_cond_3D.x[ i ][ j ][ km-1 ] ) / 2.;
			if ( h.x[ i ][ j ][ 0 ] == 1. ) 	t_cond_3D.x[ i ][ j ][ 0 ] = 0.;

			t_evap_3D.x[ i ][ j ][ 0 ] = c43 * t_evap_3D.x[ i ][ j ][ 1 ] - c13 * t_evap_3D.x[ i ][ j ][ 2 ];
			t_evap_3D.x[ i ][ j ][ km-1 ] = c43 * t_evap_3D.x[ i ][ j ][ km-2 ] - c13 * t_evap_3D.x[ i ][ j ][ km-3 ];
			t_evap_3D.x[ i ][ j ][ 0 ] = t_evap_3D.x[ i ][ j ][ km-1 ] = ( t_evap_3D.x[ i ][ j ][ 0 ] + t_evap_3D.x[ i ][ j ][ km-1 ] ) / 2.;
			if ( h.x[ i ][ j ][ 0 ] == 1. ) 	t_evap_3D.x[ i ][ j ][ 0 ] = 0.;

			BuoyancyForce.x[ i ][ j ][ 0 ] = c43 * BuoyancyForce.x[ i ][ j ][ 1 ] - c13 * BuoyancyForce.x[ i ][ j ][ 2 ];
			BuoyancyForce.x[ i ][ j ][ km-1 ] = c43 * BuoyancyForce.x[ i ][ j ][ km-2 ] - c13 * BuoyancyForce.x[ i ][ j ][ km-3 ];
			BuoyancyForce.x[ i ][ j ][ 0 ] = BuoyancyForce.x[ i ][ j ][ km-1 ] = ( BuoyancyForce.x[ i ][ j ][ 0 ] + BuoyancyForce.x[ i ][ j ][ km-1 ] ) / 2.;
			if ( h.x[ i ][ j ][ 0 ] == 1. ) 	BuoyancyForce.x[ i ][ j ][ 0 ] = 0.;

			cloud.x[ i ][ j ][ 0 ] = c43 * cloud.x[ i ][ j ][ 1 ] - c13 * cloud.x[ i ][ j ][ 2 ];
			cloud.x[ i ][ j ][ km-1 ] = c43 * cloud.x[ i ][ j ][ km-2 ] - c13 * cloud.x[ i ][ j ][ km-3 ];
			cloud.x[ i ][ j ][ 0 ] = cloud.x[ i ][ j ][ km-1 ] = ( cloud.x[ i ][ j ][ 0 ] + cloud.x[ i ][ j ][ km-1 ] ) / 2.;

			ice.x[ i ][ j ][ 0 ] = c43 * ice.x[ i ][ j ][ 1 ] - c13 * ice.x[ i ][ j ][ 2 ];
			ice.x[ i ][ j ][ km-1 ] = c43 * ice.x[ i ][ j ][ km-2 ] - c13 * ice.x[ i ][ j ][ km-3 ];
			ice.x[ i ][ j ][ 0 ] = ice.x[ i ][ j ][ km-1 ] = ( ice.x[ i ][ j ][ 0 ] + ice.x[ i ][ j ][ km-1 ] ) / 2.;

			P_rain.x[ i ][ j ][ 0 ] = c43 * P_rain.x[ i ][ j ][ 1 ] - c13 * P_rain.x[ i ][ j ][ 2 ];
			P_rain.x[ i ][ j ][ km-1 ] = c43 * P_rain.x[ i ][ j ][ km-2 ] - c13 * P_rain.x[ i ][ j ][ km-3 ];
			P_rain.x[ i ][ j ][ 0 ] = P_rain.x[ i ][ j ][ km-1 ] = ( P_rain.x[ i ][ j ][ 0 ] + P_rain.x[ i ][ j ][ km-1 ] ) / 2.;

			P_snow.x[ i ][ j ][ 0 ] = c43 * P_snow.x[ i ][ j ][ 1 ] - c13 * P_snow.x[ i ][ j ][ 2 ];
			P_snow.x[ i ][ j ][ km-1 ] = c43 * P_snow.x[ i ][ j ][ km-2 ] - c13 * P_snow.x[ i ][ j ][ km-3 ];
			P_snow.x[ i ][ j ][ 0 ] = P_snow.x[ i ][ j ][ km-1 ] = ( P_snow.x[ i ][ j ][ 0 ] + P_snow.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}



// surface values of Evaporation, Condensation, Water, Water_super, IceAir, precipitable_water only for radial printout
	precipitable_water.y[ 0 ][ 0 ] = 0.;
	precipitable_water.y[ 0 ][ 0 ] = 0.;
	co2_total.y[ 0 ][ 0 ] = 0.;

	precipitation_NASA_average = 0.;
	precipitablewater_average = 0.;
	precipitation_average = 0.;
	Evaporation_Penman_average = 0.;
	Evaporation_Haude_average = 0.;
	co2_vegetation_average = 0.;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			co2_total.y[ j ][ k ] = co2.x[ 0 ][ j ][ k ];

			if ( Latency.x[ 0 ][ j ][ k ] >= 0. )
			{
				Condensation.y[ j ][ k ] = Latency.x[ 0 ][ j ][ k ];
			}
			else
			{
				Evaporation.y[ j ][ k ] = Latency.x[ 0 ][ j ][ k ];
			}
			LatentHeat.y[ j ][ k ] = Latency.x[ 0 ][ j ][ k ];

			for ( int i = 0; i < im; i++ )
			{
				e = c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep; 													// water vapour pressure in hPa
				a = 216.6 * e / ( t.x[ i ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3

				precipitable_water.y[ j ][ k ] += a * L_atm / ( double ) ( im - 1 );					// kg/m³ * m
			}
		}
	}


	double coeff_prec = 86400. * 1000. / u_0;

// surface values of precipitation and precipitable water
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Precipitation.y[ j ][ k ] = coeff_prec * ( P_rain.x[ 0 ][ j ][ k ] + P_snow.x[ 0 ][ j ][ k ] );							// 60 s * 60 min * 24 h = 86400 s == 1 d
																																											// P_rain and P_snow in kg/ ( m² * s )
																																											// Precipitation in kg/ ( m² * d ) == mm / d
																																											// u_0 from a non-dimensionalisation process
																																											// 1000. from kg to g
			if ( Precipitation.y[ j ][ k ] <= 0 )			Precipitation.y[ j ][ k ] = 0.;

			precipitable_water.y[ j ][ k ] = precipitable_water.y[ j ][ k ] / 1000.;														// divided by water density ( 1000 kg/m³ ) results in m compares as well to mm ( absolute values identical )
		}
	}




// averaging of precipitation in mm/a and precipitable water in mm
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			precipitation_NASA_average += precipitation_NASA.y[ j ][ k ];
			precipitablewater_average += precipitable_water.y[ j ][ k ];
			precipitation_average += Precipitation.y[ j ][ k ];

//			Evaporation_Penman_average += Evaporation_Penman.y[ j ][ k ];
//			Evaporation_Haude_average += Evaporation_Haude.y[ j ][ k ];

			Evaporation_Penman_average += 0.;
			Evaporation_Haude_average += 0.;

			co2_vegetation_average += co2_total.y[ j ][ k ];
		}
	}

	co2_vegetation_average = co2_vegetation_average / ( double ) ( ( jm -1 ) * ( km - 1 ) );
	precipitablewater_average = precipitablewater_average / ( double ) ( ( jm -1 ) * ( km - 1 ) );
	precipitation_average = 365. * precipitation_average / ( double ) ( ( jm -1 ) * ( km - 1 ) );
	precipitation_NASA_average = 365. * precipitation_NASA_average / ( double ) ( ( jm -1 ) * ( km - 1 ) );
	Evaporation_Penman_average = 365. * Evaporation_Penman_average / ( double ) ( ( jm -1 ) * ( km - 1 ) );
	Evaporation_Haude_average = 365. * Evaporation_Haude_average / ( double ) ( ( jm -1 ) * ( km - 1 ) );




// boundaries of Q_bottom, Q_latent, Q_sensible
	for ( int j = 0; j < jm; j++ )
	{
		Q_bottom.y[ j ][ 0 ] = c43 * Q_bottom.y[ j ][ 1 ] - c13 * Q_bottom.y[ j ][ 2 ];
		Q_bottom.y[ j ][ km-1 ] = c43 * Q_bottom.y[ j ][ km-2 ] - c13 * Q_bottom.y[ j ][ km-3 ];
		Q_bottom.y[ j ][ 0 ] = Q_bottom.y[ j ][ km-1 ] = ( Q_bottom.y[ j ][ 0 ] + Q_bottom.y[ j ][ km-1 ] ) / 2.;

		Q_latent.y[ j ][ 0 ] = c43 * Q_latent.y[ j ][ 1 ] - c13 * Q_latent.y[ j ][ 2 ];
		Q_latent.y[ j ][ km-1 ] = c43 * Q_latent.y[ j ][ km-2 ] - c13 * Q_latent.y[ j ][ km-3 ];
		Q_latent.y[ j ][ 0 ] = Q_latent.y[ j ][ km-1 ] = ( Q_latent.y[ j ][ 0 ] + Q_latent.y[ j ][ km-1 ] ) / 2.;

		Q_sensible.y[ j ][ 0 ] = c43 * Q_sensible.y[ j ][ 1 ] - c13 * Q_sensible.y[ j ][ 2 ];
		Q_sensible.y[ j ][ km-1 ] = c43 * Q_sensible.y[ j ][ km-2 ] - c13 * Q_sensible.y[ j ][ km-3 ];
		Q_sensible.y[ j ][ 0 ] = Q_sensible.y[ j ][ km-1 ] = ( Q_sensible.y[ j ][ 0 ] + Q_sensible.y[ j ][ km-1 ] ) / 2.;

		Evaporation_Penman.y[ j ][ 0 ] = c43 * Evaporation_Penman.y[ j ][ 1 ] - c13 * Evaporation_Penman.y[ j ][ 2 ];
		Evaporation_Penman.y[ j ][ km-1 ] = c43 * Evaporation_Penman.y[ j ][ km-2 ] - c13 * Evaporation_Penman.y[ j ][ km-3 ];
		Evaporation_Penman.y[ j ][ 0 ] = Evaporation_Penman.y[ j ][ km-1 ] = ( Evaporation_Penman.y[ j ][ 0 ] + Evaporation_Penman.y[ j ][ km-1 ] ) / 2.;

		Evaporation_Haude.y[ j ][ 0 ] = c43 * Evaporation_Haude.y[ j ][ 1 ] - c13 * Evaporation_Haude.y[ j ][ 2 ];
		Evaporation_Haude.y[ j ][ km-1 ] = c43 * Evaporation_Haude.y[ j ][ km-2 ] - c13 * Evaporation_Haude.y[ j ][ km-3 ];
		Evaporation_Haude.y[ j ][ 0 ] = Evaporation_Haude.y[ j ][ km-1 ] = ( Evaporation_Haude.y[ j ][ 0 ] + Evaporation_Haude.y[ j ][ km-1 ] ) / 2.;

		Precipitation.y[ j ][ 0 ] = c43 * Precipitation.y[ j ][ 1 ] - c13 * Precipitation.y[ j ][ 2 ];
		Precipitation.y[ j ][ km-1 ] = c43 * Precipitation.y[ j ][ km-2 ] - c13 * Precipitation.y[ j ][ km-3 ];
		Precipitation.y[ j ][ 0 ] = Precipitation.y[ j ][ km-1 ] = ( Precipitation.y[ j ][ 0 ] + Precipitation.y[ j ][ km-1 ] ) / 2.;

		Evaporation_Haude.y[ j ][ 0 ] = c43 * Evaporation_Haude.y[ j ][ 1 ] - c13 * Evaporation_Haude.y[ j ][ 2 ];
		Evaporation_Haude.y[ j ][ km-1 ] = c43 * Evaporation_Haude.y[ j ][ km-2 ] - c13 * Evaporation_Haude.y[ j ][ km-3 ];
		Evaporation_Haude.y[ j ][ 0 ] = Evaporation_Haude.y[ j ][ km-1 ] = ( Evaporation_Haude.y[ j ][ 0 ] + Evaporation_Haude.y[ j ][ km-1 ] ) / 2.;

		Evaporation_Haude.y[ j ][ 0 ] = c43 * Evaporation_Haude.y[ j ][ 1 ] - c13 * Evaporation_Haude.y[ j ][ 2 ];
		Evaporation_Haude.y[ j ][ km-1 ] = c43 * Evaporation_Haude.y[ j ][ km-2 ] - c13 * Evaporation_Haude.y[ j ][ km-3 ];
		Evaporation_Haude.y[ j ][ 0 ] = Evaporation_Haude.y[ j ][ km-1 ] = ( Evaporation_Haude.y[ j ][ 0 ] + Evaporation_Haude.y[ j ][ km-1 ] ) / 2.;

	}



	for ( int k = 0; k < km; k++ )
	{
		Q_bottom.y[ 0 ][ k ] = c43 * Q_bottom.y[ 1 ][ k ] - c13 * Q_bottom.y[ 2 ][ k ];
		Q_bottom.y[ jm-1 ][ k ] = c43 * Q_bottom.y[ jm-2 ][ k ] - c13 * Q_bottom.y[ jm-3 ][ k ];

		Q_latent.y[ 0 ][ k ] = c43 * Q_latent.y[ 1 ][ k ] - c13 * Q_latent.y[ 2 ][ k ];
		Q_latent.y[ jm-1 ][ k ] = c43 * Q_latent.y[ jm-2 ][ k ] - c13 * Q_latent.y[ jm-3 ][ k ];

		Q_sensible.y[ 0 ][ k ] = c43 * Q_sensible.y[ 1 ][ k ] - c13 * Q_sensible.y[ 2 ][ k ];
		Q_sensible.y[ jm-1 ][ k ] = c43 * Q_sensible.y[ jm-2 ][ k ] - c13 * Q_sensible.y[ jm-3 ][ k ];

		Evaporation_Penman.y[ 0 ][ k ] = c43 * Evaporation_Penman.y[ 1 ][ k ] - c13 * Evaporation_Penman.y[ 2 ][ k ];
		Evaporation_Penman.y[ jm-1 ][ k ] = c43 * Evaporation_Penman.y[ jm-2 ][ k ] - c13 * Evaporation_Penman.y[ jm-3 ][ k ];

		Evaporation_Haude.y[ 0 ][ k ] = c43 * Evaporation_Haude.y[ 1 ][ k ] - c13 * Evaporation_Haude.y[ 2 ][ k ];
		Evaporation_Haude.y[ jm-1 ][ k ] = c43 * Evaporation_Haude.y[ jm-2 ][ k ] - c13 * Evaporation_Haude.y[ jm-3 ][ k ];

		Precipitation.y[ 0 ][ k ] = c43 * Precipitation.y[ 1 ][ k ] - c13 * Precipitation.y[ 2 ][ k ];
		Precipitation.y[ jm-1 ][ k ] = c43 * Precipitation.y[ jm-2 ][ k ] - c13 * Precipitation.y[ jm-3 ][ k ];

		Evaporation_Haude.y[ 0 ][ k ] = c43 * Evaporation_Haude.y[ 1 ][ k ] - c13 * Evaporation_Haude.y[ 2 ][ k ];
		Evaporation_Haude.y[ jm-1 ][ k ] = c43 * Evaporation_Haude.y[ jm-2 ][ k ] - c13 * Evaporation_Haude.y[ jm-3 ][ k ];

		Evaporation_Haude.y[ 0 ][ k ] = c43 * Evaporation_Haude.y[ 1 ][ k ] - c13 * Evaporation_Haude.y[ 2 ][ k ];
		Evaporation_Haude.y[ jm-1 ][ k ] = c43 * Evaporation_Haude.y[ jm-2 ][ k ] - c13 * Evaporation_Haude.y[ jm-3 ][ k ];

	}






	cout.precision ( 2 );

// printout of surface data at one predefinded location

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
	name_Value_6 = " Evaporation Haude ";
	name_Value_7 = " precipitable water average ";
	name_Value_8 = " precipitation average per year ";
	name_Value_9 = " precipitation average per day ";
	name_Value_10 = " precipitation NASA average per year ";
	name_Value_11 = " precipitation NASA average per day ";
	name_Value_12 = " Evaporation_Penman_average per year ";
	name_Value_13 = " Evaporation_Penman_average per day ";
	name_Value_14 = " Evaporation_Haude_average per year ";
	name_Value_15 = " Evaporation_Haude_average per day ";
	name_Value_16 = " latent heat 2D ";
	name_Value_17 = " Condensation heat 2D ";
	name_Value_18 = " Evaporation heat 2D ";
	name_Value_19 = " latent heat surf ";
	name_Value_20 = " Condensation heat surf ";
	name_Value_21 = " Evaporation heat surf ";
	name_Value_22 = " co2_average ";
	name_Value_23 = " precipitable water ";
	name_Value_24 = " precipitation ";

	name_unit_wm2 = " W/m2";
	name_unit_mmd = " mm/d";
	name_unit_mm = " mm";
	name_unit_mma = " mm/a";
	name_unit_ppm = " ppm";

	heading = " printout of surface data at predefinded locations: level, latitude, longitude";
	heading_Dresden = " City of Dresden, Germany, Europe";
	heading_Sydney = " City of Sydney, New South Wales, Australia";
	heading_Equator = " Equator in the central Pacific";

	cout << endl << endl << heading << endl << endl;

	int choice = { 1 };

	preparation:

	switch ( choice )
	{
		case 1 :	cout << heading_Dresden << endl;
						i_loc_level = 0;																		// sea level
						j_loc = 39;																				// 51°N, Dresden Germany
						k_loc = 346;																			// 14°W, Dresden Germany
						break;

		case 2 :	cout << heading_Sydney << endl;
						i_loc_level = 0;																		// sea level
						j_loc = 123;																			// 33°S, Dresden Germany
						k_loc = 151;																			// 151°E, Dresden Germany
						break;

		case 3 :	cout << heading_Equator << endl;
						i_loc_level = 0;																		// sea level
						j_loc = 90;																				// 0°N, Equator
						k_loc = 180;																			// 180°E, central Pacific
						break;

	default : 	cout << choice << "error in iterationPrintout member function in class Accuracy" << endl;
	}


	if ( j_loc <= 90 )
	{
		j_loc_deg = 90 - j_loc;
		deg_lat = deg_north;
	}

	if ( j_loc > 90 )
	{
		j_loc_deg = j_loc - 90;
		deg_lat = deg_south;
	}


	if ( k_loc <= 180 )
	{
		k_loc_deg = k_loc;
		deg_lon = deg_east;
	}

	if ( k_loc > 180 )
	{
		k_loc_deg = 360 - k_loc;
		deg_lon = deg_east;
	}


	Value_1 = Q_Radiation.y[ j_loc ][ k_loc ];
	Value_2 = Q_latent.y[ j_loc ][ k_loc ];
	Value_3 = Q_sensible.y[ j_loc ][ k_loc ];
	Value_4 = Q_bottom.y[ j_loc ][ k_loc ];

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_1 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_1 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_2 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_2 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_3 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_3 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_4 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_4 << setw ( 6 ) << name_unit_wm2 << endl;

	Value_17 = Latency.x[ 0 ][ j_loc ][ k_loc ];
	Value_18 = t_cond_3D.x[ 0 ][ j_loc ][ k_loc ];
	Value_19 = t_evap_3D.x[ 0 ][ j_loc ][ k_loc ];
	Value_23 = precipitable_water.y[ j_loc ][ k_loc ];

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_23 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_23 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_19 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_17 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_20 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_18 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_21 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_19 << setw ( 6 ) << name_unit_wm2 << endl;


	Value_14 = Condensation.y[ j_loc ][ k_loc ] + Evaporation.y[ j_loc ][ k_loc ];
	Value_15 = Condensation.y[ j_loc ][ k_loc ];
	Value_16 = Evaporation.y[ j_loc ][ k_loc ];
	Value_24 = Precipitation.y[ j_loc ][ k_loc ];

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_24 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_24 << setw ( 6 ) << name_unit_mmd << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_16 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_14 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_17 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_15 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_18 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_16 << setw ( 6 ) << name_unit_wm2 << endl;

	Value_5 = Evaporation_Penman.y[ j_loc ][ k_loc ];
	Value_6 = Evaporation_Haude.y[ j_loc ][ k_loc ];

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon << "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_5 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_5 << setw ( 6 ) << name_unit_mmd << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_6 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_6 << setw ( 6 ) << name_unit_mmd << endl << endl;


	choice++;
	if ( choice <= 3 ) goto preparation;

	cout << endl;




	Value_7 = precipitablewater_average;
	Value_8 = precipitation_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_8 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_8 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_9 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_8 / 365. << setw ( 6 ) << name_unit_mmd << endl;


	Value_10 = precipitation_NASA_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_10 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_10 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_11 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_10 / 365. << setw ( 6 ) << name_unit_mmd << endl;


	Value_13 = Evaporation_Haude_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_14 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_13 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_15 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_13 / 365. << setw ( 6 ) << name_unit_mmd << endl;


	Value_9 = co2_vegetation_average * co2_0;
	Value_12 = Evaporation_Penman_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_22 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_9 << setw ( 6 ) << name_unit_ppm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_12 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_12 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_13 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_12 / 365. << setw ( 6 ) << name_unit_mmd << endl << endl << endl;
}




