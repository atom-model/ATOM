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


void Results_MSL_Atm::run_MSL_data ( Array &hc, Array &cc, Array &tc, Array &pc, Array &uc, Array &Rai, Array &Rai_su, Array &Ic, Array &Lat, Array &cond_3D, Array &evap_3D, Array_2D &prec, Array_2D &rai, Array_2D &rai_su, Array_2D &ice, Array_2D &evap, Array_2D &cond, Array_2D &prec_wat, Array_2D &Q_Bal, Array_2D &Q_evap, Array_2D &Q_lat, Array_2D &Q_sen, Array_2D &Q_dif, Array_2D &evap_Pen, Array_2D &evap_Hau, Array_2D &t_j, Array_2D &c_j )
{
// calculation of a total quantity as sum on all values in a virtual column in r-direction

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			prec.y[ j ][ k ] = 0.;					// precipitation
			ice.y[ j ][ k ] = 0.; 						// ice
			rai.y[ j ][ k ] = 0.;						// rain water
			rai_su.y[ j ][ k ] = 0.;					// supercooled water
			evap.y[ j ][ k ] = 0.;					// evaporation
			cond.y[ j ][ k ] = 0.;					// condensation
			prec_wat.y[ j ][ k ] = 0.;			// precipitable water
			Q_Bal.y[ j ][ k ] = 0.;				// radiation balance
			Q_lat.y[ j ][ k ] = 0.;					// latent heat
			Q_sen.y[ j ][ k ] = 0.;				// sensible heat
			Q_dif.y[ j ][ k ] = 0.;					// bottom heat
			evap_Pen.y[ j ][ k ] = 0.;			// evaporation by Penman
			evap_Hau.y[ j ][ k ] = 0.;			// evaporation by Haude
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			p_baro = ( r_0_air * 287.1 * tc.x[ 0 ][ j ][ k ] * ( t_0 + 20. ) ) / 100.;											// surface pressure from barometric formula in hPa

			Q_evap.y[ j ][ k ] = ( 2500.8 - 2.372 * ( tc.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) * 1.e-3;							// heat of evaporation of water in [MJ/kg] (Kuttler)

			t_Celsius_SL = tc.x[ 0 ][ j ][ k ] * t_0 - t_0;																				// transforming Kelvin into Celsius
			e_SL = ( r_0_water_vapour * R_WaterVapour * tc.x[ 0 ][ j ][ k ] * t_0 ) * .01;
			a_SL = 216.6 * e_SL / ( tc.x[ 0 ][ j ][ k ] * t_0 );																		// absolute humidity in kg/m3

			c_grad = ( - 3. * cc.x[ 0 ][ j ][ k ] + 4. * cc.x[ 1 ][ j ][ k ] - cc.x[ 2 ][ j ][ k ] ) / ( 2. * dr );				// water vapour pressure gradient in Kg/(Kg m)
			t_grad = ( - 3. * tc.x[ 0 ][ j ][ k ] + 4. * tc.x[ 1 ][ j ][ k ] - tc.x[ 2 ][ j ][ k ] ) / ( 2. * dr );					// temperature gradient in K/m

			if ( hc.x[ 0 ][ j ][ k ] == 0. )
			{
				if ( sun == 0. ) Q_Bal.y[ j ][ k ] = sigma * ( 1. - albedo ) * pow ( tc.x[ 0 ][ j ][ k ] * t_0, 4. );	// radiation balance by prescription of the surface temperature in W/m2
				Q_lat.y[ j ][ k ] = - a_SL * lv * c_0 * coeff_Diffusion_latent * c_grad / L_atm;					// latente heat in [W/m2] from energy transport equation
				Q_sen.y[ j ][ k ] = - r_0_air * cp_l * coeff_Diffusion_sensibel * tc.x[ 0 ][ j ][ k ] * t_0 * t_grad / L_atm;	// sensible heat in [W/m2] from energy transport equation
			}
			else 
			{
				Q_Bal.y[ j ][ k ] = 0.;
				Q_lat.y[ j ][ k ] = 0.;
				Q_sen.y[ j ][ k ] = 0.;
			}

			t_denom = t_Celsius_SL + 237.3;
			Delta = 4098. * ( .6108 * exp ( ( 17.27 * t_Celsius_SL )  / t_denom ) ) / ( t_denom * t_denom );// gradient of the water vapour pressure curve
			E_Rain = hp * exp ( 17.0809 * t_Celsius_SL / ( 234.175 + t_Celsius_SL ) );					// saturation vapour pressure in the water phase for t > 0°C in hPa
			sat_deficit = E_Rain - e_SL;																								// saturation deficit
			gamma = p_baro * cp_l / ( ep * lv );																					// Psychrometer constant
			E_a = .35 * ( 1. + .15 * uc.x[ 1 ][ j ][ k ] * .01 ) * sat_deficit * .75;										// ventilation-humidity member for Penmans formula

			evap_Pen.y[ j ][ k ] = ( Q_Bal.y[ j ][ k ] * Delta * .75 + gamma * E_a ) / ( Delta * .75 + gamma );
			evap_Hau.y[ j ][ k ] = f_Haude * sat_deficit;
			if ( evap_Pen.y[ j ][ k ] <= 0. ) evap_Pen.y[ j ][ k ] = 0.;
			if ( evap_Hau.y[ j ][ k ] <= 0. ) evap_Hau.y[ j ][ k ] = 0.;
																																						// simplified formula for evaporation over day length of 12h by Haude, Häckel
																																						// coefficient per day about 0.30 for gras-evapotranspiration given, Kuttler
			Q_dif.y[ j ][ k ] = Q_Bal.y[ j ][ k ] - Q_lat.y[ j ][ k ] - Q_sen.y[ j ][ k ];										// difference understood as heat of the ground


			for ( int i = 0; i < im; i++ )
			{
//				p_h = exp ( - 9.8066 * (  double ) i * 500. / ( R_Air * tc.x[ i ][ j ][ k ] * t_0 ) ) * p_0;	// current air pressure, step size in 500 m, from barometric formula in hPa
//				e_h = ( r_0_water_vapour * R_WaterVapour * tc.x[ i ][ j ][ k ] * t_0 ) * .01;
//				a_h = 216.6 * e_h / ( tc.x[ i ][ j ][ k ] * t_0 );																	// absolute humidity in kg/m3 compares to density of water vapour

				rai.y[ j ][ k ] += Rai.x[ i ][ j ][ k ];
				rai_su.y[ j ][ k ] += Rai_su.x[ i ][ j ][ k ];
				ice.y[ j ][ k ] += Ic.x[ i ][ j ][ k ];

				prec_wat.y[ j ][ k ] += cc.x[ i ][ j ][ k ];
			}

			for ( int i = 0; i < im; i++ )
			{
				if ( Lat.x[ i ][ j ][ k ] >= .0 )			evap.y[ j ][ k ] += Lat.x[ i ][ j ][ k ];
				if ( Lat.x[ i ][ j ][ k ] <= .0 )			cond.y[ j ][ k ] += - Lat.x[ i ][ j ][ k ];
			}

			for ( int i = 0; i < im; i++ )
			{
				if ( Lat.x[ i ][ j ][ k ] >= .0 )			evap_3D.x[ i ][ j ][ k ] = Lat.x[ i ][ j ][ k ];
				if ( Lat.x[ i ][ j ][ k ] <= .0 )			cond_3D.x[ i ][ j ][ k ] = - Lat.x[ i ][ j ][ k ];
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
//			prec.y[ j ][ k ] = ( rai.y[ j ][ k ] + rai_su.y[ j ][ k ] + ice.y[ j ][ k ] ) * coeff_mmWS;	// precipitation consists of 100% saturated water vapour + 100% supercooled water + 100% ice
//			prec.y[ j ][ k ] = ( rai.y[ j ][ k ] + .5 * rai_su.y[ j ][ k ] ) * coeff_mmWS;						// precipitation consists of 100% saturated water vapour + 50% supercooled water
			prec.y[ j ][ k ] = ( rai.y[ j ][ k ] + rai_su.y[ j ][ k ] ) * coeff_mmWS;						// precipitation consists of 100% saturated water vapour + 100% supercooled water
//			prec.y[ j ][ k ] = rai.y[ j ][ k ] * coeff_mmWS;														// precipitation consists of 100% saturated water vapour
			prec_wat.y[ j ][ k ] = prec_wat.y[ j ][ k ] * coeff_mmWS;
		}
	}
}








void Results_MSL_Atm::show_MSL_data ( Array &hc, Array &cc, Array &tc, Array &pc, Array &uc, Array &Rai, Array &Es, Array &Lat, Array_2D &prec, Array_2D &ice, Array_2D &evap, Array_2D &cond, Array_2D &prec_wat, Array_2D &Q_Bal, Array_2D &Q_evap, Array_2D &Q_lat, Array_2D &Q_sen, Array_2D &Q_dif, Array_2D &evap_Pen, Array_2D &evap_Hau  )
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

	name_unit_wm2 = " W/m2";
	name_unit_mm = " mm/d";

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


	Value_1 = Q_Bal.y[ j_loc ][ k_loc ];
	Value_2 = Q_lat.y[ j_loc ][ k_loc ];
	Value_3 = Q_sen.y[ j_loc ][ k_loc ];
	Value_4 = Q_dif.y[ j_loc ][ k_loc ];

	cout << endl << endl << heading << endl << endl;

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_1 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_1 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_2 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_2 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_3 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_3 << setw ( 6 ) << name_unit_wm2 << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_4 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_4 << setw ( 6 ) << name_unit_wm2 << endl;


	Value_5 = evap_Pen.y[ j_loc ][ k_loc ];
	Value_6 = evap_Hau.y[ j_loc ][ k_loc ];

	cout << setw ( 6 ) << i_loc_level << setw ( 2 ) << level << setw ( 5 ) << j_loc_deg << setw ( 3 ) << deg_lat << setw ( 4 ) << k_loc_deg << setw ( 3 ) << deg_lon<< "  " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_5 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_5 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 20 ) << setfill ( '.' ) << name_Value_6 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_6 << setw ( 6 ) << name_unit_mm << endl << endl;
}






void Results_MSL_Atm::land_oceanFraction ( Array &hc )
{
// calculation of the ratio ocean to land, also addition and substraction of CO2 of land, ocean and vegetation

	h_point_max =  ( jm - 1 ) * ( km - 1 );

	h_land = 0;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( hc.x[ 0 ][ j ][ k ] == 1. )		h_land++;
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




void Results_MSL_Atm::vegetationDistribution ( double max_Precipitation, Array_2D &prec, Array_2D &veg, Array &tc, Array &hc )
{
// description or vegetation areas following the local dimensionsles values of precipitation, maximum value is 1

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( ( hc.x[ 0 ][ j ][ k ] == 1. ) && ( tc.x[ 0 ][ j ][ k ] >= t_0 ) ) veg.y[ j ][ k ] = prec.y[ j ][ k ] / max_Precipitation;			// actual vegetation areas
			else veg.y[ j ][ k ] = 0.;
		}
	}
}
