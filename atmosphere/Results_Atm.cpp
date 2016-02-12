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

#include "Results_Atm.h"

using namespace std;


Results_MSL_Atm::Results_MSL_Atm ( int im, int jm, int km, int sun, double g, double ep, double hp, double u_0, double p_0, double t_0, double c_0, double co2_0, double sigma, double albedo_extra, double lv, double cp_l, double L_atm, double dr, double dthe, double dphi, double r_0_air, double R_Air, double r_0_water_vapour, double R_WaterVapour, double co2_vegetation, double co2_ocean, double co2_land, double gam )
:	coeff_Diffusion_latent ( .005 ),													// diffusion coefficient for latent heat in [m²/s]    0.015
	coeff_Diffusion_sensibel ( .2 ),													// diffusion coefficient for sensible heat in [m²/s]     0.03
	f_Haude ( .3 )																			// Haude factor for evapotranspiration 0.3 for low, dense vegetation as raw average value by Kuttler
{
	coeff_mmWS = r_0_air / r_0_water_vapour;									// coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]
	coeff_lv = lv / ( cp_l * t_0 );														// coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_lv = 9.1069 in [ / ]

	c43 = 4./3.;
	c13 = 1./3.;

	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> L_atm = L_atm;
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
	this-> albedo_extra = albedo_extra;
	this-> lv = lv;
	this-> gam = gam;
	this-> cp_l = cp_l;
	this-> r_0_air = r_0_air;
	this-> R_Air = R_Air;
	this-> r_0_water_vapour = r_0_water_vapour;
	this-> R_WaterVapour = R_WaterVapour;
	this-> co2_vegetation = co2_vegetation;
	this-> co2_ocean = co2_ocean;
	this-> co2_land = co2_land;
}


Results_MSL_Atm::~Results_MSL_Atm () {}



void Results_MSL_Atm::run_MSL_data ( int RadiationModel, double max_Precipitation, double max_CO2_total, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &c, Array &co2, Array &t, Array &p_dyn, Array &p_stat, Array &BuoyancyForce, Array &u, Array &v, Array &w, Array &Rain, Array &Rain_super, Array &Ice, Array &Latency, Array &Q_Sensible, Array &radiation_3D, Array &t_cond_3D, Array &t_evap_3D, Array &aux_u, Array &aux_v, Array &aux_w, Array_2D &Precipitation, Array_2D &precipitation_j, Array_2D &Water, Array_2D &Water_super, Array_2D &IceAir, Array_2D &Evaporation, Array_2D &Condensation, Array_2D &LatentHeat, Array_2D &precipitable_water, Array_2D &Q_Radiation, Array_2D &Q_Evaporation, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Q_bottom, Array_2D &Evaporation_Penman, Array_2D &Evaporation_Haude, Array_2D &Vegetation, Array_2D &Radiation_Balance, Array_2D &Radiation_Balance_par, Array_2D &Radiation_Balance_bot, Array_2D &albedo, Array_2D &co2_total )
{
// determination of temperature and pressure by the law of Clausius-Clapeyron for water vapour concentration
// reaching saturation of water vapour pressure leads to formation of rain or ice
// precipitation and cloud formation by formulas from Häckel
// dry adiabatic lapse rate and saturated adiabatic lapse rate = temperature decrease with hight
// SL stands for sea level

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			t_Celsius_SL = t.x[ 0 ][ j ][ k ] * t_0 - t_0;																	// conversion from Kelvin to Celsius at sea surface = NN
			p_SL = p_stat.x[ 0 ][ j ][ k ] = ( r_0_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;												// from gas equation given in hPa

			e_SL = c.x[ 0 ][ j ][ k ] * p_SL / ep; 																						// water vapour pressure in hPa
			a_SL = 216.6 * e_SL / ( t.x[ 0 ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3 at sea level
			q_SL = c.x[ 0 ][ j ][ k ];																							// threshold value for water vapour at sea level in kg/kg

			E_Rain_SL = hp * exp ( 17.0809 * t_Celsius_SL / ( 234.175 + t_Celsius_SL ) );		// saturation water vapour pressure for the water phase at t > 0°C in hPa
			E_Rain_super_SL = hp * exp ( 17.8436 * t_Celsius_SL / ( 245.425 + t_Celsius_SL ) );	// saturation water vapour pressure for the water phase at t < 0°C, supercooled in hPa
			E_Ice_SL = hp * exp ( 22.4429 * t_Celsius_SL / ( 272.44 + t_Celsius_SL ) );			// saturation water vapour pressure for the ice phase in hPa

			q_Rain_SL = ep * E_Rain_SL / ( p_SL - E_Rain_SL );																				// water vapour amount at saturation with water formation in kg/kg
			q_Rain_super_SL = ep * E_Rain_super / ( p_SL - E_Rain_super );														// water vapour amount at saturation with water formation in kg/kg
			q_Ice_SL = ep * E_Ice_SL / ( p_SL - E_Ice_SL );																				// water vapour amount at saturation with ice formation in kg/kg

			Rain.x[ 0 ][ j ][ k ] = ( q_SL - q_Rain_SL );																	// liquid water as surplus from the local saturated water vapour  in g/kg
			Rain_super.x[ 0 ][ j ][ k ] = ( q_SL - q_Rain_super_SL );												// liquid water as surplus from the  
			Ice.x[ 0 ][ j ][ k ] = ( q_SL - q_Ice_SL );																			// ice formation as surplus above the supercooled saturated water vapour in g/kg

			t_dew_SL = ( 423.86 - 234.175 * log ( e_SL ) ) / ( log ( e_SL ) - 18.89 );					// dewpoint temperature on ground in °C 		by Häckel
			h_level = 122. * ( t_Celsius_SL - t_dew_SL );															// Condensation level in m		by Häckel + correction
			i_level = ( int ) ( h_level / 500. );																					// Condensation level in radial steps
			sat_deficit = E_Rain_SL - e_SL;																						// saturation deficit, if positive then saturation is less than 100%
			RF_e = e_SL / E_Rain_SL * 100.;																					// relative humidity at any point in %

//	if ( ( j == 67 ) && ( k == 278 ) )		cout << 0 << "   " << t_Celsius_SL << "   " << p_SL << "   " << p_SL << "   " << e_SL << "   " << E_Rain_SL << "   " << a_SL << "   " << q_SL << "   " << q_Rain_SL << "   " << t_dew_SL << "   " << h_level << "   " << sat_deficit << "   " << RF_e << endl;  // Havana comparison with the vertical water vapour distribution 

			for ( int i = 1; i < im; i++ )
			{
				t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;																		// conversion from Kelvin to Celsius
				p_h = p_stat.x[ i ][ j ][ k ] = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( r_0_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) ) * p_0;

				e_h = c.x[ i ][ j ][ k ] * p_h / ep; 																						// water vapour pressure in hPa
				a_h = 216.6 * e_h / ( t.x[ i ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3
				q_h = c.x[ i ][ j ][ k ];																							// threshold value for water vapour at local hight h in kg/kg

				E_Rain = hp * exp ( 17.0809 * t_Celsius / ( 234.175 + t_Celsius ) );						// saturation water vapour pressure for the water phase at t > 0°C in hPa
				E_Rain_super = hp * exp ( 17.8436 * t_Celsius / ( 245.425 + t_Celsius ) );			// saturation water vapour pressure for the water phase at t < 0°C, supercooled in hPa
				E_Ice = hp * exp ( 22.4429 * t_Celsius / ( 272.44 + t_Celsius ) );							// saturation water vapour pressure for the ice phase in hPa

				q_Rain  = ep * E_Rain / ( p_h - E_Rain );																					// water vapour amount at saturation with water formation in kg/kg
				q_Rain_super  = ep * E_Rain_super / ( p_h - E_Rain_super );															// water vapour amount at saturation with water formation in kg/kg
				q_Ice  = ep * E_Ice / ( p_h - E_Ice );																						// water vapour amount at saturation with ice formation in kg/kg

// precipitation and cloud formation from formulas by Häckel
// h stands for the local position

				h_h = - R_Air * t.x[ i ][ j ][ k ] * t_0 / g * log ( p_h / p_0 );										// barometric elevation formula solved for h corresponding to hight over ground in m
				t_dew = ( 423.86 - 234.175 * log ( e_h ) ) / ( log ( e_h ) - 18.89 );							// current dewpoint temperature in °C
				sat_deficit = E_Rain - e_h;																						// saturation deficit, if positive then saturation is less than 100%
				RF_e = e_h / E_Rain * 100.;																					// relative humidity at any point in %

//	if ( ( j == 67 ) && ( k == 278 ) )		cout << i << "   " << t_Celsius << "   " << p_h << "   " << e_h << "   " << E_Rain << "   " << a_h << "   " << q_h << "   " << q_Rain << "   " << t_dew << "   " << h_h << "   " << sat_deficit << "   " << RF_e << endl;  // Havana comparison with the vertical water vapour distribution 

// application of threshhold values for water vapour to compute rain, super cooled water and ice

				if ( ( q_h > q_Rain ) && ( Latency.x[ i ][ j ][ k ] <= 0. ) && ( t_Celsius >= 0. ) )			Rain.x[ i ][ j ][ k ] = ( q_h - q_Rain );// 		liquid water as surplus from the local saturated water vapour in kg/kg
				else	Rain.x[ i ][ j ][ k ] = 0.;

				if ( ( q_h > q_Rain_super ) && ( Latency.x[ i ][ j ][ k ] <= 0. ) && ( t_Celsius < 0. ) && ( t_Celsius >= - 20. ) )		Rain_super.x[ i ][ j ][ k ] = ( q_h - q_Rain_super );		 // liquid water as surplus from the local saturated supercooled water vapour in kg/kg
				else	Rain_super.x[ i ][ j ][ k ] = 0.;

				if ( ( q_h > q_Ice ) && ( Latency.x[ i ][ j ][ k ] <= 0. ) && ( t_Celsius < -20. ) )			Ice.x[ i ][ j ][ k ] = ( q_h - q_Ice );		// ice as surplus from the local saturated ice kg/kg
				else	Ice.x[ i ][ j ][ k ] = 0.;


/*
// printout for various thermodynamical quantities for the preticipation computations along the equator ( j = 90 )
				if ( ( i == 5 ) && ( j == 90 ) && ( k == 180 ) )
				{
					cout << endl;
					cout << " i = " << i << "   j = " << j << "   k = " << k << "   i_level = " << i_level  << "   h_level (m) = " << h_level << "   h_h (m) = " << h_h << endl << endl;

					cout << " t_h (°C) = " << t_Celsius << "   p_h (hPa) = " << p_h << "   a_h (g/m3) = " << a_h << "   c_h (g/Kg) = " << c.x[ i ][ j ][ k ] * 1000. << "   q_Rain (g/Kg) = " << q_Rain * 1000. << "   e_h (hPa) = " << e_h << "   E_Rain (hPa) = " << E_Rain  << endl << endl;

					cout << " Rain (g/kg) = " << Rain.x[ i ][ j ][ k ] * 1000. << "   Rain_super (g/kg) = " << Rain_super.x[ i ][ j ][ k ] * 1000. << "   Ice (g/kg) = " << Ice.x[ i ][ j ][ k ] * 1000. << "   E_Ice (hPa) = " << E_Ice << "   sat_deficit (hPa) = " << sat_deficit << "   t_dew (°C) = " << t_dew << "   E_Rain_SL (hPa) = " << E_Rain_SL << endl << endl;

					cout << " t_SL (°C) = " << t_Celsius_SL << "   p_SL (hPa) = " << p_SL << "   a_SL (g/m3) = " << a_SL << "   c_SL (g/Kg) = " << c.x[ 0 ][ j ][ k ] * 1000. << "   e_SL (hPa) = " << e_SL << "   q_Ice (g/Kg) = " << q_Ice * 1000. << "   t_dew_SL (°C) = " << t_dew_SL << endl << endl;

					cout << " Evap_Haude (mm/d) = " << Evap_Haude << "   RF_e (%) = " << RF_e  << "   E_Rain_super (hPa) = " << E_Rain_super << "   q_Rain_super (g/Kg) = " << q_Rain_super * 1000. << "   q_h (g/Kg) = " << q_h * 1000. << "   q_SL (g/Kg) = " << q_SL * 1000. << endl << endl;
				}
*/
			}
		}
	}


// calculation of a total quantity as sum on all values in a virtual column in r-direction

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( RadiationModel == 0 ) Q_Radiation.y[ j ][ k ] = Radiation_Balance_par.y[ j ][ k ]; // parabolic radiation balance assumed
			if ( RadiationModel >= 2 ) Q_Radiation.y[ j ][ k ] = radiation_3D.x[ 0 ][ j ][ k ];		 // two- and multi-layer radiation balance assumed

			Precipitation.y[ j ][ k ] = 0.;												// precipitation
			IceAir.y[ j ][ k ] = 0.; 														// ice
			Water.y[ j ][ k ] = 0.;															// rain water
			Water_super.y[ j ][ k ] = 0.;												// supercooled water
			Evaporation.y[ j ][ k ] = 0.;												// Evaporation
			Condensation.y[ j ][ k ] = 0.;												// Condensation
			LatentHeat.y[ j ][ k ] = 0.;													// Condensation
			precipitable_water.y[ j ][ k ] = 0.;										// precipitable water
			Evaporation_Penman.y[ j ][ k ] = 0.;									// Evaporation by Penman
			Evaporation_Haude.y[ j ][ k ] = 0.;										// Evaporation by Haude

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
					if ( i == 0 ) 	p_stat.x[ 0 ][ j ][ k ] = ( r_0_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;// given in hPa
					else 	p_stat.x[ i ][ j ][ k ] = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( r_0_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) ) * p_0;
																																			// current air pressure, step size in 500 m, from a polytropic atmosphere in hPa
					t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;																	// transforming Kelvin into Celsius

					e = c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep; 													// water vapour pressure in hPa
					a = 216.6 * e / ( t.x[ i ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3

					t_denom = t_Celsius + 234.175;
					E = hp * exp ( 17.0809 * t_Celsius / t_denom );													// saturation vapour pressure in the water phase for t > 0°C in hPa

					Delta = 4000. * E / ( t_denom * t_denom );															// gradient of the water vapour pressure curve in hPa/K, coef = 234.175 * 17.0809
					sat_deficit = ( E - e );																							// saturation deficit in hPa/K
					gamma = p_stat.x[ 0 ][ j ][ k ] * cp_l / ( ep * lv );												// Psychrometer constant in hPa/K

					c_grad = ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i + 1 ][ j ][ k ] - c.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );	// water vapour pressure gradient in g/(Kg m)
					t_grad = ( - 3. * t.x[ i ][ j ][ k ] + 4. * t.x[ i + 1 ][ j ][ k ] - t.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );	// temperature gradient in K/m

					E_a = .35 * ( 1. + .15 * sqrt ( v.x[ i + 1 ][ j ][ k ] * v.x[ i + 1 ][ j ][ k ] + w.x[ i + 1 ][ j ][ k ] * w.x[ i + 1 ][ j ][ k ] ) ) * sat_deficit;	// ventilation-humidity Penmans formula

					Q_Evaporation.y[ j ][ k ] = ( 2500.8 - 2.372 * ( t.x[ i ][ j ][ k ] * t_0 - t_0 ) ) * 1.e-3;	// heat of Evaporation of water in [MJ/kg] (Kuttler) = variable lv
					Q_latent.y[ j ][ k ] = - a * lv * coeff_Diffusion_latent * c_grad / L_atm;					// latente heat in [W/m2] from energy transport equation
					Q_sensible.y[ j ][ k ] = - r_0_air * cp_l * coeff_Diffusion_sensibel * t_grad * t_0 / L_atm;		// sensible heat in [W/m2] from energy transport equation
					Q_bottom.y[ j ][ k ] = Q_Radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ];	// difference understood as heat into the ground

					Evaporation_Haude.y[ j ][ k ] = f_Haude * sat_deficit;											// simplified formula for Evaporation over day length of 12h by Haude, Häckel
					if ( Evaporation_Haude.y[ j ][ k ] <= 0. ) 		Evaporation_Haude.y[ j ][ k ] = 0.;
					Evaporation_Penman.y[ j ][ k ] = .0346 * ( ( Q_Radiation.y[ j ][ k ] - Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma );
					if ( Evaporation_Penman.y[ j ][ k ] <= 0. ) 		Evaporation_Penman.y[ j ][ k ] = 0.;	// .0346 coefficient W/m2 corresponds to mm/d (Kraus)
				}

// only on the sea surface
				if ( ( i == 0 ) && ( h.x[ 0 ][ j ][ k ] == 0. ) )
				{
					if ( i == 0 ) 	p_stat.x[ 0 ][ j ][ k ] = ( r_0_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;		// given in hPa
					else 	p_stat.x[ i ][ j ][ k ] = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( r_0_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) ) * p_0;

					t_Celsius = t.x[ 0 ][ j ][ k ] * t_0 - t_0;																	// transforming Kelvin into Celsius

					e = c.x[ 0 ][ j ][ k ] * p_stat.x[ 0 ][ j ][ k ] / ep; 													// water vapour pressure in hPa
					a = 216.6 * e / ( t.x[ 0 ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3

					t_denom = t_Celsius + 234.175;
					E = hp * exp ( 17.0809 * t_Celsius / t_denom );													// saturation vapour pressure in the water phase for t > 0°C in hPa
					Delta = 4000. * E / ( t_denom * t_denom );															// gradient of the water vapour pressure curve in hPa/K, coef = 234.175 * 17.0809

					sat_deficit = ( E - e );																							// saturation deficit in hPa/K
					gamma = p_stat.x[ 0 ][ j ][ k ] * cp_l / ( ep * lv );												// Psychrometer constant in hPa/K

					c_grad = ( - 3. * c.x[ 0 ][ j ][ k ] + 4. * c.x[ 1 ][ j ][ k ] - c.x[ 2 ][ j ][ k ] ) / ( 2. * dr );	// water vapour pressure gradient in g/(Kg m)
					t_grad = ( - 3. * t.x[ 0 ][ j ][ k ] + 4. * t.x[ 1 ][ j ][ k ] - t.x[ 2 ][ j ][ k ] ) / ( 2. * dr );	// temperature gradient in K/m

					E_a = .35 * ( 1. + .15 * sqrt ( v.x[ 1 ][ j ][ k ] * v.x[ 1 ][ j ][ k ] + w.x[ 1 ][ j ][ k ] * w.x[ 1 ][ j ][ k ] ) ) * sat_deficit;	// ventilation-humidity for Penmans formula

					Q_Evaporation.y[ j ][ k ] = ( 2500.8 - 2.372 * ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) * 1.e-3;	// heat of Evaporation of water in [MJ/kg] (Kuttler) = lv
					Q_latent.y[ j ][ k ] = - a * lv * coeff_Diffusion_latent * c_grad / L_atm;		// latente heat in [W/m2] from energy transport equation
					Q_sensible.y[ j ][ k ] = - r_0_air * cp_l * coeff_Diffusion_sensibel * t_grad * t_0 / L_atm;		// sensible heat in [W/m2] from energy transport equation
					Q_bottom.y[ j ][ k ] = Q_Radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ];	// difference understood as heat of the ground

					Evaporation_Haude.y[ j ][ k ] = f_Haude * sat_deficit;										// simplified formula for Evaporation over day length of 12h by Haude, Häckel
					if ( Evaporation_Haude.y[ j ][ k ] <= 0. ) 		Evaporation_Haude.y[ j ][ k ] = 0.;
					Evaporation_Penman.y[ j ][ k ] = .0346 * ( ( Q_Radiation.y[ j ][ k ] - Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma );
					if ( Evaporation_Penman.y[ j ][ k ] <= 0. ) 		Evaporation_Penman.y[ j ][ k ] = 0.;	// .0346 coefficient W/m2 corresponds to mm/d (Kraus)

				}
			}

		}
	}




// boundaries of latent heat

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			p_dyn.x[ 0 ][ j ][ k ] = c43 * p_dyn.x[ 1 ][ j ][ k ] - c13 * p_dyn.x[ 2 ][ j ][ k ];
			p_dyn.x[ im-1 ][ j ][ k ] = c43 * p_dyn.x[ im-2 ][ j ][ k ] - c13 * p_dyn.x[ im-3 ][ j ][ k ];

			Latency.x[ 0 ][ j ][ k ] = c43 * Latency.x[ 1 ][ j ][ k ] - c13 * Latency.x[ 2 ][ j ][ k ];
			Latency.x[ im-1 ][ j ][ k ] = c43 * Latency.x[ im-2 ][ j ][ k ] - c13 * Latency.x[ im-3 ][ j ][ k ];

			Q_Sensible.x[ 0 ][ j ][ k ] = c43 * Q_Sensible.x[ 1 ][ j ][ k ] - c13 * Q_Sensible.x[ 2 ][ j ][ k ];
			Q_Sensible.x[ im-1 ][ j ][ k ] = c43 * Q_Sensible.x[ im-2 ][ j ][ k ] - c13 * Q_Sensible.x[ im-3 ][ j ][ k ];

			t_cond_3D.x[ 0 ][ j ][ k ] = c43 * t_cond_3D.x[ 1 ][ j ][ k ] - c13 * t_cond_3D.x[ 2 ][ j ][ k ];
			t_cond_3D.x[ im-1 ][ j ][ k ] = c43 * t_cond_3D.x[ im-2 ][ j ][ k ] - c13 * t_cond_3D.x[ im-3 ][ j ][ k ];
//			t_cond_3D.x[ im-1 ][ j ][ k ] = t_cond_3D.x[ im-4 ][ j ][ k ] - 3. * t_cond_3D.x[ im-3 ][ j ][ k ] + 3. * t_cond_3D.x[ im-2 ][ j ][ k ];		// extrapolation

			t_evap_3D.x[ 0 ][ j ][ k ] = c43 * t_evap_3D.x[ 1 ][ j ][ k ] - c13 * t_evap_3D.x[ 2 ][ j ][ k ];
			t_evap_3D.x[ im-1 ][ j ][ k ] = c43 * t_evap_3D.x[ im-2 ][ j ][ k ] - c13 * t_evap_3D.x[ im-3 ][ j ][ k ];
//			t_evap_3D.x[ im-1 ][ j ][ k ] = t_evap_3D.x[ im-4 ][ j ][ k ] - 3. * t_evap_3D.x[ im-3 ][ j ][ k ] + 3. * t_evap_3D.x[ im-2 ][ j ][ k ];		// extrapolation

			BuoyancyForce.x[ 0 ][ j ][ k ] = c43 * BuoyancyForce.x[ 1 ][ j ][ k ] - c13 * BuoyancyForce.x[ 2 ][ j ][ k ];
			BuoyancyForce.x[ im-1 ][ j ][ k ] = c43 * BuoyancyForce.x[ im-2 ][ j ][ k ] - c13 * BuoyancyForce.x[ im-3 ][ j ][ k ];
//			BuoyancyForce.x[ im-1 ][ j ][ k ] = BuoyancyForce.x[ im-4 ][ j ][ k ] - 3. * BuoyancyForce.x[ im-3 ][ j ][ k ] + 3. * BuoyancyForce.x[ im-2 ][ j ][ k ];		// extrapolation

			Rain.x[ 0 ][ j ][ k ] = c43 * Rain.x[ 1 ][ j ][ k ] - c13 * Rain.x[ 2 ][ j ][ k ];
			Rain.x[ im-1 ][ j ][ k ] = c43 * Rain.x[ im-2 ][ j ][ k ] - c13 * Rain.x[ im-3 ][ j ][ k ];

			Rain_super.x[ 0 ][ j ][ k ] = c43 * Rain_super.x[ 1 ][ j ][ k ] - c13 * Rain_super.x[ 2 ][ j ][ k ];
			Rain_super.x[ im-1 ][ j ][ k ] = c43 * Rain_super.x[ im-2 ][ j ][ k ] - c13 * Rain_super.x[ im-3 ][ j ][ k ];

			Ice.x[ 0 ][ j ][ k ] = c43 * Ice.x[ 1 ][ j ][ k ] - c13 * Ice.x[ 2 ][ j ][ k ];
			Ice.x[ im-1 ][ j ][ k ] = c43 * Ice.x[ im-2 ][ j ][ k ] - c13 * Ice.x[ im-3 ][ j ][ k ];

		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int i = 0; i < im; i++ )
		{
			Latency.x[ i ][ 0 ][ k ] = c43 * Latency.x[ i ][ 1 ][ k ] - c13 * Latency.x[ i ][ 2 ][ k ];
			Latency.x[ i ][ jm-1 ][ k ] = c43 * Latency.x[ i ][ jm-2 ][ k ] - c13 * Latency.x[ i ][ jm-3 ][ k ];

			Q_Sensible.x[ i ][ 0 ][ k ] = c43 * Q_Sensible.x[ i ][ 1 ][ k ] - c13 * Q_Sensible.x[ i ][ 2 ][ k ];
			Q_Sensible.x[ i ][ jm-1 ][ k ] = c43 * Q_Sensible.x[ i ][ jm-2 ][ k ] - c13 * Q_Sensible.x[ i ][ jm-3 ][ k ];

			t_cond_3D.x[ i ][ 0 ][ k ] = c43 * t_cond_3D.x[ i ][ 1 ][ k ] - c13 * t_cond_3D.x[ i ][ 2 ][ k ];
			t_cond_3D.x[ i ][ jm-1 ][ k ] = c43 * t_cond_3D.x[ i ][ jm-2 ][ k ] - c13 * t_cond_3D.x[ i ][ jm-3 ][ k ];

			t_evap_3D.x[ i ][ 0 ][ k ] = c43 * t_evap_3D.x[ i ][ 1 ][ k ] - c13 * t_evap_3D.x[ i ][ 2 ][ k ];
			t_evap_3D.x[ i ][ jm-1 ][ k ] = c43 * t_evap_3D.x[ i ][ jm-2 ][ k ] - c13 * t_evap_3D.x[ i ][ jm-3 ][ k ];

			BuoyancyForce.x[ i ][ 0 ][ k ] = c43 * BuoyancyForce.x[ i ][ 1 ][ k ] - c13 * BuoyancyForce.x[ i ][ 2 ][ k ];
			BuoyancyForce.x[ i ][ jm-1 ][ k ] = c43 * BuoyancyForce.x[ i ][ jm-2 ][ k ] - c13 * BuoyancyForce.x[ i ][ jm-3 ][ k ];

			Rain.x[ i ][ 0 ][ k ] = c43 * Rain.x[ i ][ 1 ][ k ] - c13 * Rain.x[ i ][ 2 ][ k ];
			Rain.x[ i ][ jm-1 ][ k ] = c43 * Rain.x[ i ][ jm-2 ][ k ] - c13 * Rain.x[ i ][ jm-3 ][ k ];

			Rain_super.x[ i ][ 0 ][ k ] = c43 * Rain_super.x[ i ][ 1 ][ k ] - c13 * Rain_super.x[ i ][ 2 ][ k ];
			Rain_super.x[ i ][ jm-1 ][ k ] = c43 * Rain_super.x[ i ][ jm-2 ][ k ] - c13 * Rain_super.x[ i ][ jm-3 ][ k ];

			Ice.x[ i ][ 0 ][ k ] = c43 * Ice.x[ i ][ 1 ][ k ] - c13 * Ice.x[ i ][ 2 ][ k ];
			Ice.x[ i ][ jm-1 ][ k ] = c43 * Ice.x[ i ][ jm-2 ][ k ] - c13 * Ice.x[ i ][ jm-3 ][ k ];
		}
	}


	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Latency.x[ i ][ j ][ 0 ] = c43 * Latency.x[ i ][ j ][ 1 ] - c13 * Latency.x[ i ][ j ][ 2 ];
			Latency.x[ i ][ j ][ km-1 ] = c43 * Latency.x[ i ][ j ][ km-2 ] - c13 * Latency.x[ i ][ j ][ km-3 ];
			Latency.x[ i ][ j ][ 0 ] = Latency.x[ i ][ j ][ km-1 ] = ( Latency.x[ i ][ j ][ 0 ] + Latency.x[ i ][ j ][ km-1 ] ) / 2.;

			Q_Sensible.x[ i ][ j ][ 0 ] = c43 * Q_Sensible.x[ i ][ j ][ 1 ] - c13 * Q_Sensible.x[ i ][ j ][ 2 ];
			Q_Sensible.x[ i ][ j ][ km-1 ] = c43 * Q_Sensible.x[ i ][ j ][ km-2 ] - c13 * Q_Sensible.x[ i ][ j ][ km-3 ];
			Q_Sensible.x[ i ][ j ][ 0 ] = Q_Sensible.x[ i ][ j ][ km-1 ] = ( Q_Sensible.x[ i ][ j ][ 0 ] + Q_Sensible.x[ i ][ j ][ km-1 ] ) / 2.;

			t_cond_3D.x[ i ][ j ][ 0 ] = c43 * t_cond_3D.x[ i ][ j ][ 1 ] - c13 * t_cond_3D.x[ i ][ j ][ 2 ];
			t_cond_3D.x[ i ][ j ][ km-1 ] = c43 * t_cond_3D.x[ i ][ j ][ km-2 ] - c13 * t_cond_3D.x[ i ][ j ][ km-3 ];
			t_cond_3D.x[ i ][ j ][ 0 ] = t_cond_3D.x[ i ][ j ][ km-1 ] = ( t_cond_3D.x[ i ][ j ][ 0 ] + t_cond_3D.x[ i ][ j ][ km-1 ] ) / 2.;

			t_evap_3D.x[ i ][ j ][ 0 ] = c43 * t_evap_3D.x[ i ][ j ][ 1 ] - c13 * t_evap_3D.x[ i ][ j ][ 2 ];
			t_evap_3D.x[ i ][ j ][ km-1 ] = c43 * t_evap_3D.x[ i ][ j ][ km-2 ] - c13 * t_evap_3D.x[ i ][ j ][ km-3 ];
			t_evap_3D.x[ i ][ j ][ 0 ] = t_evap_3D.x[ i ][ j ][ km-1 ] = ( t_evap_3D.x[ i ][ j ][ 0 ] + t_evap_3D.x[ i ][ j ][ km-1 ] ) / 2.;

			BuoyancyForce.x[ i ][ j ][ 0 ] = c43 * BuoyancyForce.x[ i ][ j ][ 1 ] - c13 * BuoyancyForce.x[ i ][ j ][ 2 ];
			BuoyancyForce.x[ i ][ j ][ km-1 ] = c43 * BuoyancyForce.x[ i ][ j ][ km-2 ] - c13 * BuoyancyForce.x[ i ][ j ][ km-3 ];
			BuoyancyForce.x[ i ][ j ][ 0 ] = BuoyancyForce.x[ i ][ j ][ km-1 ] = ( BuoyancyForce.x[ i ][ j ][ 0 ] + BuoyancyForce.x[ i ][ j ][ km-1 ] ) / 2.;

			Rain.x[ i ][ j ][ 0 ] = c43 * Rain.x[ i ][ j ][ 1 ] - c13 * Rain.x[ i ][ j ][ 2 ];
			Rain.x[ i ][ j ][ km-1 ] = c43 * Rain.x[ i ][ j ][ km-2 ] - c13 * Rain.x[ i ][ j ][ km-3 ];
			Rain.x[ i ][ j ][ 0 ] = Rain.x[ i ][ j ][ km-1 ] = ( Rain.x[ i ][ j ][ 0 ] + Rain.x[ i ][ j ][ km-1 ] ) / 2.;

			Rain_super.x[ i ][ j ][ 0 ] = c43 * Rain_super.x[ i ][ j ][ 1 ] - c13 * Rain_super.x[ i ][ j ][ 2 ];
			Rain_super.x[ i ][ j ][ km-1 ] = c43 * Rain_super.x[ i ][ j ][ km-2 ] - c13 * Rain_super.x[ i ][ j ][ km-3 ];
			Rain_super.x[ i ][ j ][ 0 ] = Rain_super.x[ i ][ j ][ km-1 ] = ( Rain_super.x[ i ][ j ][ 0 ] + Rain_super.x[ i ][ j ][ km-1 ] ) / 2.;

			Ice.x[ i ][ j ][ 0 ] = c43 * Ice.x[ i ][ j ][ 1 ] - c13 * Ice.x[ i ][ j ][ 2 ];
			Ice.x[ i ][ j ][ km-1 ] = c43 * Ice.x[ i ][ j ][ km-2 ] - c13 * Ice.x[ i ][ j ][ km-3 ];
			Ice.x[ i ][ j ][ 0 ] = Ice.x[ i ][ j ][ km-1 ] = ( Ice.x[ i ][ j ][ 0 ] + Ice.x[ i ][ j ][ km-1 ] ) / 2.;
		}
	}



// surface values of Evaporation, Condensation, Water, Water_super, IceAir, precipitable_water only for radial printout
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( Rain.x[ i ][ j ][ k ] >= 0. )		Water.y[ j ][ k ] += Rain.x[ i ][ j ][ k ];

				if ( Rain_super.x[ i ][ j ][ k ] >= 0. )		Water_super.y[ j ][ k ] += Rain_super.x[ i ][ j ][ k ];

				if ( Ice.x[ i ][ j ][ k ] >= 0. )		IceAir.y[ j ][ k ] += Ice.x[ i ][ j ][ k ];

				e = c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep; 													// water vapour pressure in hPa
				a = 216.6 * e / ( t.x[ i ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3

				precipitable_water.y[ j ][ k ] += a * ( double ) i * 500.;									//  kg/m³ * m

				co2_total.y[ j ][ k ] += co2.x[ i ][ j ][ k ];

				if ( Latency.x[ 0 ][ j ][ k ] >= 0. )
				{
					Condensation.y[ j ][ k ] = Latency.x[ 0 ][ j ][ k ];
				}
				else
				{
					Evaporation.y[ j ][ k ] = Latency.x[ 0 ][ j ][ k ];
				}
				LatentHeat.y[ j ][ k ] = Latency.x[ 0 ][ j ][ k ];
			}
		}
	}



// surface values of precipitation and precipitable water
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Precipitation.y[ j ][ k ] = coeff_mmWS * ( Water.y[ j ][ k ] + Water_super.y[ j ][ k ] + IceAir.y[ j ][ k ] );// precipitation consists of water vapour + supercooled water + ice
			precipitable_water.y[ j ][ k ] = precipitable_water.y[ j ][ k ] / 1000.;		// divided by water density ( 1000 kg/m³ ) results in m compares as well to mm ( absolute values identical )
		}
	}




// averaging of precipitation in mm/a and precipitable water in mm
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			precipitation_average += Precipitation.y[ j ][ k ];
			precipitation_NASA_average += precipitation_j.y[ j ][ k ];
			precipitablewater_average += precipitable_water.y[ j ][ k ];
			Evaporation_Penman_average += Evaporation_Penman.y[ j ][ k ];
			Evaporation_Haude_average += Evaporation_Haude.y[ j ][ k ];

			Water.y[ j ][ k ] = coeff_mmWS * Water.y[ j ][ k ];
			Water_super.y[ j ][ k ] = coeff_mmWS * Water_super.y[ j ][ k ];
			IceAir.y[ j ][ k ] = coeff_mmWS * IceAir.y[ j ][ k ];
			co2_vegetation += co2_total.y[ j ][ k ];
		}
	}

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
	}

}








void Results_MSL_Atm::show_MSL_data ( Array &h, Array &c, Array &t, Array &p_dyn, Array &u, Array &Rain, Array &Ice, Array &Latency, Array &Q_Sensible, Array &t_cond_3D, Array &t_evap_3D, Array_2D &Precipitation, Array_2D &precipitation_j, Array_2D &IceAir, Array_2D &Evaporation, Array_2D &Condensation, Array_2D &precipitable_water, Array_2D &Q_Radiation, Array_2D &Q_Evaporation, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Q_bottom, Array_2D &Evaporation_Penman, Array_2D &Evaporation_Haude )
{
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

	name_unit_wm2 = " W/m2";
	name_unit_mmd = " mm/d";
	name_unit_mm = " mm";
	name_unit_mma = " mm/a";

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

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_1 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_1 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_19 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_17 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_20 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_18 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_21 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_19 << setw ( 6 ) << name_unit_wm2 << endl;


	Value_14 = Condensation.y[ j_loc ][ k_loc ] + Evaporation.y[ j_loc ][ k_loc ];
	Value_15 = Condensation.y[ j_loc ][ k_loc ];
	Value_16 = Evaporation.y[ j_loc ][ k_loc ];

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_1 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_1 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_16 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_14 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_17 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_15 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 25 ) << setfill ( '.' ) << name_Value_18 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_16 << setw ( 6 ) << name_unit_wm2 << endl;

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


	Value_12 = Evaporation_Penman_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_12 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_12 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_13 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_12 / 365. << setw ( 6 ) << name_unit_mmd << endl;

/*
	Value_13 = Evaporation_Haude_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_14 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_13 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_15 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_13 / 365. << setw ( 6 ) << name_unit_mmd << endl << endl << endl;
*/
}


