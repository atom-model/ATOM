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
#include <cstring>

#include "Results_Atm.h"

using namespace std;


Results_MSL_Atm::Results_MSL_Atm ( int im, int jm, int km, int sun, double ep, double hp, double p_0, double t_0, double c_0, double sigma, double albedo, double lv, double cp_l, double L_atm, double dr, double r_0_air, double R_Air, double r_0_water_vapour, double R_WaterVapour, double co2_vegetation, double co2_ocean, double co2_land )
:	coeff_Diffusion_latent ( 1. ),													// diffusion coefficient for latent heat in [m2/s]
	coeff_Diffusion_sensibel ( .01 ),												// diffusion coefficient for sensible heat in [m2/s]
	f_Haude ( 12. * .30 )																// Haude factor for evapotranspiration, 12 for day, 0.3 for low, 
																									// dense vegetation as raw average value by Kuttler

{
	coeff_mmWS = r_0_air / r_0_water_vapour;
	f_Haude = f_Haude * .5;															// test factor

	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> L_atm = L_atm;
	this-> dr = dr;
	this-> sun = sun;
	this-> ep = ep;
	this-> hp = hp;
	this-> p_0 = p_0;
	this-> t_0 = t_0;
	this-> c_0 = c_0;
	this-> sigma = sigma;
	this-> albedo = albedo;
	this-> lv = lv;
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


void Results_MSL_Atm::run_MSL_data ( Array &h, Array &c, Array &t, Array &p, Array &u, Array &Rain, Array &Rain_super, Array &Ice, Array &Latency, Array &Condensation_3D, Array &Evaporation_3D, Array_2D &Precipitation, Array_2D &Water, Array_2D &Water_super, Array_2D &IceAir, Array_2D &Evaporation, Array_2D &Condensation, Array_2D &precipitable_water, Array_2D &Q_Balance_Radiation, Array_2D &Q_Evaporation, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Q_diff, Array_2D &Evaporation_Penman, Array_2D &Evaporation_Haude, Array_2D &t_j, Array_2D &c_j )
{
// calculation of a total quantity as sum on all values in a virtual column in r-direction

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Precipitation.y[ j ][ k ] = 0.;								// precipitation
			IceAir.y[ j ][ k ] = 0.; 										// ice
			Water.y[ j ][ k ] = 0.;											// rain water
			Water_super.y[ j ][ k ] = 0.;								// supercooled water
			Evaporation.y[ j ][ k ] = 0.;								// evaporation
			Condensation.y[ j ][ k ] = 0.;							// condensation
			precipitable_water.y[ j ][ k ] = 0.;						// precipitable water
			Q_Balance_Radiation.y[ j ][ k ] = 0.;				// radiation balance
			Q_latent.y[ j ][ k ] = 0.;										// latent heat
			Q_sensible.y[ j ][ k ] = 0.;									// sensible heat
			Q_diff.y[ j ][ k ] = 0.;											// bottom heat
			Evaporation_Penman.y[ j ][ k ] = 0.;				// evaporation by Penman
			Evaporation_Haude.y[ j ][ k ] = 0.;					// evaporation by Haude
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			p_baro = ( r_0_air * 287.1 * t.x[ 0 ][ j ][ k ] * ( t_0 + 20. ) ) / 100.;										// surface pressure from barometric formula in hPa

			Q_Evaporation.y[ j ][ k ] = ( 2500.8 - 2.372 * ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) * 1.e-3;				// heat of evaporation of water in [MJ/kg] (Kuttler)

			t_Celsius_SL = t.x[ 0 ][ j ][ k ] * t_0 - t_0;																			// transforming Kelvin into Celsius
			e_SL = ( r_0_water_vapour * R_WaterVapour * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;
			a_SL = 216.6 * e_SL / ( t.x[ 0 ][ j ][ k ] * t_0 );																		// absolute humidity in kg/m3

			c_grad = ( - 3. * c.x[ 0 ][ j ][ k ] + 4. * c.x[ 1 ][ j ][ k ] - c.x[ 2 ][ j ][ k ] ) / ( 2. * dr );				// water vapour pressure gradient in Kg/(Kg m)
			t_grad = ( - 3. * t.x[ 0 ][ j ][ k ] + 4. * t.x[ 1 ][ j ][ k ] - t.x[ 2 ][ j ][ k ] ) / ( 2. * dr );					// temperature gradient in K/m

			if ( h.x[ 0 ][ j ][ k ] == 0. )
			{
				if ( sun == 0. ) Q_Balance_Radiation.y[ j ][ k ] = sigma * ( 1. - albedo ) * pow ( t.x[ 0 ][ j ][ k ] * t_0, 4. );	// radiation balance by the surface temperature in W/m2
				Q_latent.y[ j ][ k ] = - a_SL * lv * c_0 * coeff_Diffusion_latent * c_grad / L_atm;			// latente heat in [W/m2] from energy transport equation
				Q_sensible.y[ j ][ k ] = - r_0_air * cp_l * coeff_Diffusion_sensibel * t.x[ 0 ][ j ][ k ] * t_0 * t_grad / L_atm;	// sensible heat in [W/m2] from energy transport equation
			}
			else 
			{
				Q_Balance_Radiation.y[ j ][ k ] = 0.;
				Q_latent.y[ j ][ k ] = 0.;
				Q_sensible.y[ j ][ k ] = 0.;
			}

			t_denom = t_Celsius_SL + 237.3;
			Delta = 4098. * ( .6108 * exp ( ( 17.27 * t_Celsius_SL )  / t_denom ) ) / ( t_denom * t_denom );// gradient of the water vapour pressure curve
			E_Rain = hp * exp ( 17.0809 * t_Celsius_SL / ( 234.175 + t_Celsius_SL ) );						// saturation vapour pressure in the water phase for t > 0°C in hPa
			sat_deficit = E_Rain - e_SL;																									// saturation deficit
			gamma = p_baro * cp_l / ( ep * lv );																						// Psychrometer constant
			E_a = .35 * ( 1. + .15 * u.x[ 1 ][ j ][ k ] * .01 ) * sat_deficit * .75;												// ventilation-humidity member for Penmans formula

			Evaporation_Penman.y[ j ][ k ] = ( Q_Balance_Radiation.y[ j ][ k ] * Delta * .75 + gamma * E_a ) / ( Delta * .75 + gamma );
			Evaporation_Haude.y[ j ][ k ] = f_Haude * sat_deficit;
			if ( Evaporation_Penman.y[ j ][ k ] <= 0. ) Evaporation_Penman.y[ j ][ k ] = 0.;
			if ( Evaporation_Haude.y[ j ][ k ] <= 0. ) Evaporation_Haude.y[ j ][ k ] = 0.;
																																						// simplified formula for evaporation over day length of 12h by Haude, Häckel
																																						// coefficient per day about 0.30 for gras-evapotranspiration given, Kuttler
			Q_diff.y[ j ][ k ] = Q_Balance_Radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ];										// difference understood as heat of the ground


			for ( int i = 0; i < im; i++ )
			{
				Water.y[ j ][ k ] += Rain.x[ i ][ j ][ k ];
				Water_super.y[ j ][ k ] += Rain_super.x[ i ][ j ][ k ];
				IceAir.y[ j ][ k ] += Ice.x[ i ][ j ][ k ];

				precipitable_water.y[ j ][ k ] += c.x[ i ][ j ][ k ];
			}

			for ( int i = 0; i < im; i++ )
			{
				if ( Latency.x[ i ][ j ][ k ] >= .0 )			Evaporation.y[ j ][ k ] += Latency.x[ i ][ j ][ k ];
				if ( Latency.x[ i ][ j ][ k ] <= .0 )			Condensation.y[ j ][ k ] += - Latency.x[ i ][ j ][ k ];
			}

			for ( int i = 0; i < im; i++ )
			{
				if ( Latency.x[ i ][ j ][ k ] >= .0 )			Evaporation_3D.x[ i ][ j ][ k ] = Latency.x[ i ][ j ][ k ];
				if ( Latency.x[ i ][ j ][ k ] <= .0 )			Condensation_3D.x[ i ][ j ][ k ] = - Latency.x[ i ][ j ][ k ];
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
//			Precipitation.y[ j ][ k ] = ( Water.y[ j ][ k ] + Water_super.y[ j ][ k ] + IceAir.y[ j ][ k ] ) * coeff_mmWS;	// precipitation consists of 100% saturated water vapour + 100% supercooled water + 100% ice
//			Precipitation.y[ j ][ k ] = ( Water.y[ j ][ k ] + .5 * Water_super.y[ j ][ k ] ) * coeff_mmWS;						// precipitation consists of 100% saturated water vapour + 50% supercooled water
			Precipitation.y[ j ][ k ] = ( Water.y[ j ][ k ] + Water_super.y[ j ][ k ] ) * coeff_mmWS;							// precipitation consists of 100% saturated water vapour + 100% supercooled water
//			Precipitation.y[ j ][ k ] = Water.y[ j ][ k ] * coeff_mmWS;																		// precipitation consists of 100% saturated water vapour
			precipitable_water.y[ j ][ k ] = precipitable_water.y[ j ][ k ] * coeff_mmWS;
		}
	}


// averaging of precipitation in mm/a and precipitable water in mm
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			precipitation_average += Precipitation.y[ j ][ k ];
			precipitablewater_average += precipitable_water.y[ j ][ k ];
		}
	}

	precipitablewater_average = precipitablewater_average / ( double ) ( ( jm -1 ) * ( km - 1 ) );
	precipitation_average = 365. * precipitation_average / ( double ) ( ( jm -1 ) * ( km - 1 ) );

}








void Results_MSL_Atm::show_MSL_data ( Array &h, Array &c, Array &t, Array &p, Array &u, Array &Rain, Array &Es, Array &Latency, Array_2D &Precipitation, Array_2D &IceAir, Array_2D &Evaporation, Array_2D &Condensation, Array_2D &precipitable_water, Array_2D &Q_Balance_Radiation, Array_2D &Q_Evaporation, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Q_diff, Array_2D &Evaporation_Penman, Array_2D &Evaporation_Haude )
{
	cout.precision ( 2 );

// printout of surface data at one predefinded location

	level = "m";
	deg_north = "°N";
	deg_south = "°S";
	deg_west = "°W";
	deg_east = "°E";

	name_Value_1 = " radiation ";
	name_Value_2 = " latent heat ";
	name_Value_3 = " sensible heat ";
	name_Value_4 = " bottom heat ";
	name_Value_5 = " evaporation Penman ";
	name_Value_6 = " evaporation Haude ";
	name_Value_7 = " precipitable water average ";
	name_Value_8 = " precipitation average per year ";
	name_Value_9 = " precipitation average per day ";

	name_unit_wm2 = " W/m2";
	name_unit_mmd = " mm/d";
	name_unit_mm = " mm";
	name_unit_mma = " mm/a";

	heading = " printout of surface data at predefinded locations: level, latitude, longitude";

	i_loc_level = 0;
	j_loc = 90;
	k_loc = 180;


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
		k_loc_deg = 180 - k_loc;
		deg_lon = deg_west;
	}

	if ( k_loc > 180 )
	{
		k_loc_deg = k_loc - 180;
		deg_lon = deg_east;
	}


	Value_1 = Q_Balance_Radiation.y[ j_loc ][ k_loc ];
	Value_2 = Q_latent.y[ j_loc ][ k_loc ];
	Value_3 = Q_sensible.y[ j_loc ][ k_loc ];
	Value_4 = Q_diff.y[ j_loc ][ k_loc ];

	cout << endl << endl << heading << endl << endl;

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_1 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_1 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_2 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_2 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_3 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_3 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_4 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_4 << setw ( 6 ) << name_unit_wm2 << endl;


	Value_5 = Evaporation_Penman.y[ j_loc ][ k_loc ];
	Value_6 = Evaporation_Haude.y[ j_loc ][ k_loc ];

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon << "  " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_5 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_5 << setw ( 6 ) << name_unit_mmd << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_6 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_6 << setw ( 6 ) << name_unit_mmd << endl << endl;


	Value_7 = precipitablewater_average;
	Value_8 = precipitation_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_8 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_8 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_9 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_8 / 365. << setw ( 6 ) << name_unit_mmd << endl << endl;





}






void Results_MSL_Atm::land_oceanFraction ( Array &h )
{
// calculation of the ratio ocean to land, also addition and substraction of CO2 of land, ocean and vegetation

	h_point_max =  ( jm - 1 ) * ( km - 1 );

	h_land = 0;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ 0 ][ j ][ k ] == 1. )		h_land++;
		}
	}

	h_ocean = h_point_max - h_land;

 //																										for the purpose of testing
	co2_vegetation = co2_vegetation / ( double ) h_land;
	co2_ocean = co2_ocean / ( double ) h_ocean;
	co2_land = co2_land / ( double ) h_land;
	ozean_land = ( double ) h_ocean / ( double ) h_land;

	cout.precision ( 3 );

	cout << endl;
	cout << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      total number of points at constant hight " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_point_max << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the ocean surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_ocean << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      number of points on the land surface " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << h_land << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      ocean/land ratio " << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << ozean_land << endl << endl;

//	cout << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      addition of CO2 by ocean surface " << " = + " << resetiosflags ( ios::left ) << setw ( 7 ) << scientific << setfill ( ' ' ) << co2_ocean << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      addition of CO2 by land surface " << " = + " << resetiosflags ( ios::left ) << setw ( 7 ) << scientific << setfill ( ' ' ) << co2_land << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << "      substraction of CO2 by vegetation " << " = - " << resetiosflags ( ios::left ) << setw ( 7 ) << scientific << setfill ( ' ' ) << co2_vegetation << endl << setiosflags ( ios::left ) << setw ( 50 ) << "      valid for one single point on the surface"<< endl << endl;
	cout << endl;
}




void Results_MSL_Atm::vegetationDistribution ( double max_Precipitation, Array_2D &Precipitation, Array_2D &veg, Array &t, Array &h )
{
// description or vegetation areas following the local dimensionsles values of precipitation, maximum value is 1

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( t.x[ 0 ][ j ][ k ] >= t_0 ) ) veg.y[ j ][ k ] = Precipitation.y[ j ][ k ] / max_Precipitation;			// actual vegetation areas
			else veg.y[ j ][ k ] = 0.;
		}
	}
}
