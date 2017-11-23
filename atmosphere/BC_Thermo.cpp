/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
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
#include <algorithm>

#include "BC_Thermo.h"
#include "Array.h"
#include "Array_2D.h"

using namespace std;



BC_Thermo::BC_Thermo ( string &output_path, int im, int jm, int km, int i_beg, int i_max, int RadiationModel, int NASATemperature, int sun, int declination, int sun_position_lat, int sun_position_lon, int Ma, int Ma_max, int Ma_max_half, double dt, double dr, double dthe, double dphi, double g, double ep, double hp, double u_0, double p_0, double t_0, double c_0, double sigma, double epsilon_extra, double lv, double ls, double cp_l, double L_atm, double r_air, double R_Air, double r_water_vapour, double R_WaterVapour, double co2_0, double co2_cretaceous, double co2_vegetation, double co2_ocean, double co2_land, double ik, double c_tropopause, double co2_tropopause, double c_ocean, double c_land, double t_average, double co2_average, double co2_pole, double t_cretaceous, double t_cret_cor, double t_cretaceous_max, double radiation_ocean, double radiation_pole, double radiation_equator, double t_land, double t_tropopause, double t_equator, double t_pole, double gam, double epsilon_equator, double epsilon_pole, double epsilon_tropopause, double albedo_equator, double albedo_pole, double ik_equator, double ik_pole )
{
    this ->output_path = output_path;
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
	this -> i_beg = i_beg;
	this -> i_max = i_max;
	this-> L_atm = L_atm;
	this-> dt = dt;
	this-> dr = dr;
	this-> dthe = dthe;
	this-> dphi = dphi;
	this-> RadiationModel = RadiationModel;
	this-> NASATemperature = NASATemperature;
	this-> sun = sun;
	this-> g = g;
	this-> ep = ep;
	this-> hp = hp;
	this-> u_0 = u_0;
	this-> p_0 = p_0;
	this-> t_0 = t_0;
	this-> c_0 = c_0;
	this-> sigma = sigma;
	this-> albedo_equator = albedo_equator;
	this-> albedo_pole = albedo_pole;
	this-> gam = gam;
	this-> lv = lv;
	this-> ls = ls;
	this-> cp_l = cp_l;
	this-> r_air = r_air;
	this-> R_Air = R_Air;
	this-> r_water_vapour = r_water_vapour;
	this-> R_WaterVapour = R_WaterVapour;
	this-> co2_0 = co2_0;
	this-> co2_cretaceous = co2_cretaceous;
	this-> co2_vegetation = co2_vegetation;
	this-> co2_ocean = co2_ocean;
	this-> co2_land = co2_land;
	this-> ik = ik;
	this-> ik_equator = ik_equator;
	this-> ik_pole = ik_pole;
	this-> epsilon_extra = epsilon_extra;
	this-> epsilon_pole = epsilon_pole;
	this-> epsilon_tropopause = epsilon_tropopause;
	this-> epsilon_equator = epsilon_equator;
	this-> c_tropopause = c_tropopause;
	this-> co2_tropopause = co2_tropopause;
	this-> c_ocean = c_ocean;
	this-> c_land = c_land;
	this-> t_average = t_average;
	this-> co2_average = co2_average;
	this-> co2_pole = co2_pole;
	this-> Ma = Ma;
	this-> Ma_max = Ma_max;
	this-> t_cretaceous_max = t_cretaceous_max;
	this-> t_cretaceous = t_cretaceous;
	this-> t_cret_cor = t_cret_cor;
	this-> radiation_ocean = radiation_ocean;
	this-> radiation_pole = radiation_pole;
	this-> radiation_equator = radiation_equator;
	this-> t_land = t_land;
	this-> t_tropopause = t_tropopause;
	this-> t_equator = t_equator;
	this-> t_pole = t_pole;
	this-> declination = declination;
	this-> sun_position_lat = sun_position_lat;
	this-> sun_position_lon = sun_position_lon;
	this-> Ma_max_half = Ma_max_half;


	coeff_mmWS = r_air / r_water_vapour;									// coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]
	coeff_lv = lv / ( cp_l * t_0 );														// coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_lv = 9.1069 in [ / ]
	coeff_ls = ls / ( cp_l * t_0 );														// coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_ls = 10.3091 in [ / ]

	c43 = 4./3.;
	c13 = 1./3.;

	pi180 = 180./M_PI;



																									// fall velocity for water droplets of 0.1 mm compares to 0.012 m/s
																									// fall velocity for water droplets of 0.01 mm compares to 0.8 m/s
																									// fall velocity for water droplets of 0.3 mm compares to 2.4 m/s
																									// fall velocity for water droplets of 1.0 mm compares to 6.3 m/s
																									// fall velocity for water droplets of 5.0 mm compares to 14.1 m/s
	dt_rain_dim = 250.;																// dt_rain_dim is the time  in 250 s to pass dr = 500 m, 500 m / 250 s = 2 m/s fallout velocity
//	dt_rain_dim = 50.;																	// dt_rain_dim is the time  in 250 s to pass dr = 500 m, 500 m / 250 s = 2 m/s fallout velocity
//	dt_rain_dim = L_atm / u_0 * dt;												// dt_rain_dim is the time  in 250 s to pass dr = 500 m, 500 m / 250 s = 2 m/s fallout velocity

																									// fall velocity for snow flakes of 1.0 mm compares to 0.9 m/s
																									// fall velocity for snow flakes of 2.0 mm compares to 1.2 m/s
																									// fall velocity for snow flakes of 4.0 mm compares to 1.4 m/s
	dt_snow_dim = 417.;																// dt_snow_dim is the time  in 417 s to pass dr = 500 m, 500 m / 417 s = 1.2 m/s fallout velocity
//	dt_snow_dim = 500.;																// dt_snow_dim is the time  in 417 s to pass dr = 500 m, 500 m / 500 s = 1.0 m/s fallout velocity
//	dt_snow_dim = L_atm / u_0 * dt;											// dt_snow_dim is the time  in 417 s to pass dr = 500 m, 500 m / 417 s = 1.2 m/s fallout velocity

	dt_dim = L_atm / u_0 * dt;														// dimensional time step of system
//	dt_dim = 500. / u_0 * dt;														// dimensional time step of system
//	dt_dim = 50.;																			// dimensional time step of system



	cout.precision ( 8 );
	cout.setf ( ios::fixed );

// array "jm_temp_asym" for configuring data due to latitude dependent tropopause
	jm_temp_asym = 0L;

	jm_temp_asym = new double[ jm ];

	for ( int l = 0; l < jm; l++ )
	{
		jm_temp_asym[ l ] = 0;
//		cout << jm_temp_asym[ l ] << endl;
	}


// array "alfa" for Thomas algorithm
	alfa = 0L;

	alfa = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		alfa[ l ] = 0.;
//		cout << alfa[ l ] << endl;
	}


// array "beta" for Thomas algorithm
	beta = 0L;

	beta = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		beta[ l ] = 0.;
//		cout << beta[ l ] << endl;
	}

// array "AA" for the multi-layer radiation computation
	AA = 0L;

	AA = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		AA[ l ] = 0.;
//		cout << AA[ l ] << endl;
	}


// array "cloud_max" for the multi-layer radiation computation
	cloud_max = 0L;

	cloud_max = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		cloud_max[ l ] = 0.;
//		cout << cloud_max[ l ] << endl;
	}




// Array "CC" for the multi-layer radiation computation
	CC = 0L;

	CC = new double*[ im ];

	for ( int l = 0; l < im; l++ )
	{
		CC[ l ] = new double[ im ];
	}

// default values
	for ( int l = 0; l < im; l++ )
	{
		for ( int n = 0; n < im; n++ )
		{
			CC[ l ][ n ] = 0.;
		}
	}

}



BC_Thermo::~BC_Thermo()
{
	for ( int i = 0; i < im; i++ )
	{
		delete [  ] CC[ i ];
	}

	delete [  ] CC;

	delete [  ] jm_temp_asym;
	delete [  ] alfa;
	delete [  ] beta;
	delete [  ] AA;
	delete [  ] cloud_max;
}





void BC_Thermo::BC_Radiation_2D_layer ( Array_2D &albedo, Array_2D &epsilon, Array_2D &precipitable_water, Array_2D &Ik, Array_2D &Q_Radiation, Array_2D &Radiation_Balance, Array_2D &Radiation_Balance_atm, Array_2D &Radiation_Balance_bot, Array_2D &temp_eff_atm, Array_2D &temp_eff_bot, Array_2D &temp_rad, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Q_bottom, Array_2D & co2_total, Array &t, Array &c, Array &h, Array &epsilon_3D, Array &radiation_3D )
{
// class element for the computation of the radiation balance
// computation of the local temperature based on the radiation balance
// two layer model

	cout.precision ( 6 );
	cout.setf ( ios::fixed );

	pi180 = 180./M_PI;

//	j_half = ( jm -1 ) / 2 + 30;															// position of the sun at 30°S ( 90 + 30 = 120 )
	j_half = ( jm -1 ) / 2;																	// position of the sun at 0°S
	j_max = jm - 1;

	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

	k_half = ( km -1 ) / 2;																// position of the sun at 0° oder 180° ( Greenwich )
	k_max = km - 1;

	d_k_half = ( double ) k_half;
	d_k_max = ( double ) k_max;


	t_eff_earth = 255.;																		// effectiv radiation temperature of 255 K compares to -18 °C for 30% albedo_equator
																									// comparable with measurements of temperature in 5000 m hight and pressure 550 hPa
																									// sigma * t-earth ** 4 = ( 1- albedo_equator ) Ik / 4 = 239 W/m2, valid as reference value

//	j_sun = 30;																				// lateral sun location, 30 = 60°N
	j_sun = 0;																					// equatorial sun location


	ik_equator = 341.3;
	ik_pole = 60.;
	ik_co2_eff = ik_pole - ik_equator;

	albedo_equator = albedo_equator;
	albedo_pole = .7;
	albedo_co2_eff = albedo_pole - albedo_equator;

// effective temperature, albedo and emissivity/absorptivity for the two layer model
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			d_j = ( double ) j;
			if ( h.x[ 0 ][ j ][ k ]  == 0. ) 	albedo.y[ j ][ k ] = albedo_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + albedo_pole;

			if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 1 ][ j ][ k ] == 0. ) ) 	albedo.y[ j ][ k ] = albedo_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + albedo_pole;

			for ( int i = 1; i < im-1; i++ )
			{
				if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i+1 ][ j ][ k ] == 0. ) ) 	albedo.y[ j ][ k ] = albedo_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + albedo_pole;

				if ( ( t.x[ i ][ j ][ k ] * t_0 - t_0 <= 0. ) && ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i+1 ][ j ][ k ] == 0. ) ) 	albedo.y[ j ][ k ] = albedo_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + albedo_pole;
			}


//				Ik.y[ j ][ k ] = ik * sin ( ( ( double ) j - j_sun ) / pi180 ) * sin ( ( ( double ) k - 90 ) / pi180 ) / 4.;	//	solar short wave radiation on a point location
//				Ik.y[ j ][ k ] = ik * sin ( ( ( double ) j - j_sun ) / pi180 ) / 4.;	// solar short wave radiation as zonal distribution
//				if ( Ik.y[ j ][ k ] < 60. ) Ik.y[ j ][ k ] = 60.;
				Ik.y[ j ][ k ] = ik_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + ik_pole;


				if ( ( Ma == 0 ) && ( h.x[ 0 ][ j ][ k ] == 1. ) )
				{
					t_Ik = pow ( Ik.y[ j ][ k ] / sigma, ( 1. / 4. ) );
					Ik.y[ j ][ k ] = sigma * pow ( ( t_Ik ), 4. );
				}

				t_eff_earth = pow ( ( 1. - albedo.y[ j ][ k ] ) * Ik.y[ j ][ k ] / sigma, ( 1. / 4. ) );

				epsilon.y[ j ][ k ] = 2. * ( 1. - 1. / pow ( ( t.x[ 0 ][ j ][ k ] * t_0 / t_eff_earth ), 4. ) ); // bottom temperature must be known in advance
//				if ( epsilon.y[ j ][ k ] >= 1. ) 		epsilon.y[ j ][ k ] = 1.;

//				e = ( r_water_vapour * R_WaterVapour * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;				// water vapour pressure in hPa
//				e = c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep; 												// water vapour pressure in hPa
//				epsilon.y[ j ][  k ] = .595 + .0416 * sqrt ( e );										// dependency given in Häckel


//				epsilon.y[ j ][  k ] = .71;

		}
	}


// radiation and temperature prediction for the two layer model
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Radiation_Balance_bot.y[ j ][ k ] = ( 1. - albedo.y[ j ][ k ] ) / ( 1. - epsilon.y[ j ][ k ] / 2. ) * Ik.y[ j ][ k ];
			Radiation_Balance_atm.y[ j ][ k ] = epsilon.y[ j ][ k ] * ( 1. - albedo.y[ j ][ k ] ) / ( 2. - epsilon.y[ j ][ k ] ) * Ik.y[ j ][ k ];

			Radiation_Balance.y[ j ][ k ] = ( 1. - albedo.y[ j ][ k ] ) * Ik.y[ j ][ k ] + Radiation_Balance_bot.y[ j ][ k ] - Radiation_Balance_atm.y[ j ][ k ] - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ] - Q_bottom.y[ j ][ k ];

			Q_Radiation.y[ j ][ k ] = Radiation_Balance_bot.y[ j ][ k ];

			temp_eff_bot.y[ j ][ k ] = pow ( Radiation_Balance_bot.y[ j ][ k ] / sigma, ( 1. / 4. ) );
			temp_eff_atm.y[ j ][ k ] = pow ( Radiation_Balance_atm.y[ j ][ k ] / ( epsilon.y[ j ][ k ] * sigma ), ( 1. / 4. ) );

			t.x[ 0 ][ j ][ k ] = pow ( Radiation_Balance_bot.y[ j ][ k ] / sigma, ( 1. / 4. ) ) / t_0;

			temp_rad.y[ j ][ k ] = t.x[ 0 ][ j ][ k ] * t_0 - t_0;


// zero layer model, not used
			D = Ik.y[ j ][ k ] * .26;															// direct sun radiation, short-wave (Häckel)
			H = Ik.y[ j ][ k ] * .29;															// diffusive sky radiation, short-wave (Häckel)
			G = D + H;																			// total solar radiation from above = global radiation, short-wave
			R_short = G * albedo_equator;												// total reflected solar radiation, short-wave (Häckel)
			Q_short = G - R_short;														// short-wave radiation balance

			AG = Ik.y[ j ][ k ] * 1.14;														// downward terrestrial radiation, long-wave (Häckel)
			A = Ik.y[ j ][ k ] * .95;															// outgoing long-wave radiation of the surface, long-wave (Häckel)
			R_long = Ik.y[ j ][ k ] * .05;													// total reflected terrestrial radiation, long-wave (Häckel)
			Q_long = AG - A - R_long;													// long-wave radiation balance

			Q_total = Q_short + Q_long;												// total radiation balance, short- and long-wave

			Q_lat = Ik.y[ j ][ k ] * .27;														// latent energy (Häckel)
			Q_sen = Ik.y[ j ][ k ] * .05;													// sensible energy (Häckel)
			Q_bot = Q_total - Q_lat - Q_sen;

			Q_rad = Q_short + AG - Q_lat - Q_sen - Q_bot;

			t_rad = pow ( Q_rad / ( epsilon_extra * sigma ), ( 1. / 4. ) );
		}
	}



/*
// printout for various radiation quantities for the terrestrial radiation balance

	cout.precision ( 3 );

			if ( ( j == 90 ) && ( k == 180 ) )
			{
				cout << endl;
				cout << "   j = " << j << "   k = " << k << endl << endl;

				cout << "__________     radiation data in W/m²    _______________" << endl << endl;

				cout << " t_SL (°C) = " << t.x[ 0 ][ j ][ k ] * t_0 - t_0 << "    Ik / 4 (W/m²) = " << Ik.y[ j ][ k ] << "   Radiation_Balance_bot = " << Radiation_Balance_bot.y[ j ][ k ] << "   Radiation_Balance_atm = " << Radiation_Balance_atm.y[ j ][ k ] << "   temp_eff_bot = " << temp_eff_bot.y[ j ][ k ] - 273.15 << "   temp_eff_atm = " << temp_eff_atm.y[ j ][ k ] - 273.15 << endl << endl;

				cout << " D  = " << D << "   H = " << H << "   G = " << G << "   R_short = " << R_short << "   Q_short = " << Q_short << "   epsilon_extra = " << epsilon_extra  << "   epsilon = " << epsilon.y[ j ][ k ] << "   albedo_equator = " << albedo_equator << "   albedo = " << albedo.y[ j ][ k ] << endl << endl;

				cout << " AG = " << AG << "   A = " << A << "   R_long  = " << R_long << "   Q_long  = " << Q_long << "   e = " << e << "   Q_total = " << Q_total << endl << endl;

				cout << " Q_lat = " << Q_lat << "   Q_sen = " << Q_sen << "   Q_bot  = " << Q_bot << "   Q_rad  = " << Q_rad << "   t_rad [ °C ] = " << t_rad << "   Q_total = " << Q_total << endl << endl;			Radiation_Balance_bot.y[ j ][ k ] = ( 1. - albedo.y[ j ][ k ] ) / ( 1. - epsilon.y[ j ][ k ] / 2. ) * Ik.y[ j ][ k ];
			Radiation_Balance_atm.y[ j ][ k ] = epsilon.y[ j ][ k ] * ( 1. - albedo.y[ j ][ k ] ) / ( 2. - epsilon.y[ j ][ k ] ) * Ik.y[ j ][ k ];


				cout << " Q_latent = " << Q_latent.y[ j ][ k ] << "   Q_sensible = " <<  + Q_sensible.y[ j ][ k ] << "   Q_bottom  = " << Q_bottom.y[ j ][ k ] << "   Radadiation_Balance  = " << Radiation_Balance.y[ j ][ k ] << "   t [ K ] = " << t.x[ 0 ][ j ][ k ] * t_0 << "   Q_total = " << Q_total << endl << endl << endl;

				cout << "__________     radiation data in % are based on the local extra terrestrial sun radiaton Ik in W/m²    _______________" << endl << endl;

				cout << " t_h (°C) = " << t.x[ 0 ][ j ][ k ] * t_0 - t_0 << "    Ik / 4 (%) = " << Ik.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   Radiation_Balance_bot = " << Radiation_Balance_bot.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   Radiation_Balance_atm = " << Radiation_Balance_atm.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   temp_eff_bot = " << temp_eff_bot.y[ j ][ k ] - 273.15 << "   temp_eff_atm = " << temp_eff_atm.y[ j ][ k ] - 273.15 << endl << endl;

				cout << " D  = " << D / Ik.y[ j ][ k ] * 100. << "   H = " << H / Ik.y[ j ][ k ] * 100. << "   G = " << G / Ik.y[ j ][ k ] * 100. << "   R_short = " << R_short / Ik.y[ j ][ k ] * 100. << "   Q_short = " << Q_short / Ik.y[ j ][ k ] * 100. << "   epsilon = " << epsilon.y[ j ][ k ] << "   albedo = " << albedo.y[ j ][ k ] << endl << endl;

				cout << " AG = " << AG / Ik.y[ j ][ k ] * 100. << "   A = " << A / Ik.y[ j ][ k ] * 100. << "   R_long  = " << R_long / Ik.y[ j ][ k ] * 100. << "   Q_long  = " << Q_long / Ik.y[ j ][ k ] * 100. << "   e (hPa) = " << e << "   Q_total = " << Q_total / Ik.y[ j ][ k ] * 100. << endl << endl;

				cout << " Q_lat = " << Q_lat / Ik.y[ j ][ k ] * 100. << "   Q_sen = " << Q_sen / Ik.y[ j ][ k ] * 100. << "   Q_bot  = " << Q_bot / Ik.y[ j ][ k ] * 100. << "   Q_rad  = " << Q_rad / Ik.y[ j ][ k ] * 100. << "   t_rad [ °C ] = " << t_rad << "   Q_total = " << Q_total / Ik.y[ j ][ k ] * 100. << endl << endl;

				cout << " Q_latent = " << Q_latent.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   Q_sensible = " <<  + Q_sensible.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   Q_bottom  = " << Q_bottom.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   Radadiation_Balance  = " << Radiation_Balance.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   t [ °C ] = " << t.x[ 0 ][ j ][ k ] * t_0 - t_0 << "   Q_total = " << Q_total << endl << endl << endl;

			}

	cout.precision ( 6 );
*/





/*
	cout << endl << " ***** printout of 2D-field solar radiation ***** " << endl << endl;
	Ik.printArray_2D();			t.x[ 0 ][ j ][ k ] = pow ( Radiation_Balance_bot.y[ j ][ k ] / sigma, ( 1. / 4. ) ) / t_0;


	cout << endl << " ***** printout of 2D-field albedo ***** " << endl << endl;
	albedo.printArray_2D();

	cout << endl << " ***** printout of 2D-field epsilon ***** " << endl << endl;
	epsilon.printArray_2D();

	cout << endl << " ***** printout of 2D-field Radiation_Balance ***** " << endl << endl;
	Radiation_Balance.printArray_2D();

	cout << endl << " ***** printout of 2D-field Radiation_Balance_bot ***** " << endl << endl;
	Radiation_Balance_bot.printArray_2D();

	cout << endl << " ***** printout of 2D-field Radiation_Balance_atm ***** " << endl << endl;
	Radiation_Balance_atm.printArray_2D();

	cout << endl << " ***** printout of 2D-field Q_Radiation ***** " << endl << endl;
	Q_Radiation.printArray_2D();

	cout << endl << " ***** printout of 2D-field Q_latent ***** " << endl << endl;
	Q_latent.printArray_2D();

	cout << endl << " ***** printout of 2D-field Q_sensible ***** " << endl << endl;
	Q_sensible.printArray_2D();

	cout << endl << " ***** printout of 2D-field Q_bottom ***** " << endl << endl;
	Q_bottom.printArray_2D();

	cout << endl << " ***** printout of 2D-field temperature_rad ***** " << endl << endl;
	temp_rad.printArray_2D();
*/
//	cout << endl << " ***** printout of 2D-field solar radiation ***** " << endl << endl;
//	Ik.printArray_2D();


}






void BC_Thermo::BC_Radiation_multi_layer ( int *im_tropopause, int n, Array_2D &t_j, Array_2D &albedo, Array_2D &epsilon, Array_2D &precipitable_water, Array_2D &Ik, Array_2D &Q_Radiation, Array_2D &Radiation_Balance, Array_2D &Radiation_Balance_atm, Array_2D &Radiation_Balance_bot, Array_2D &temp_eff_atm, Array_2D &temp_eff_bot, Array_2D &temp_rad, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Q_bottom, Array_2D & co2_total, Array &p_stat, Array &t, Array &c, Array &h, Array &epsilon_3D, Array &radiation_3D, Array &cloud, Array &ice )
{
// class element for the computation of the radiation and the temperature distribution
// computation of the local temperature based on short and long wave radiation
// multi layer radiation model

	cout.precision ( 4 );
	cout.setf ( ios::fixed );

	pi180 = 180./M_PI;

//	j_half = ( jm -1 ) / 2 + 30;															// position of the sun at 30°S ( 90 + 30 = 120 )
	j_half = ( jm -1 ) / 2;																	// position of the sun at 0°S
	j_max = jm - 1;

	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

	k_half = ( km -1 ) / 2;																// position of the sun at 0° oder 180° ( Greenwich )
	k_max = km - 1;

	d_k_half = ( double ) k_half;
	d_k_max = ( double ) k_max;

//	j_sun = 30;																				// lateral sun location, 30 = 60°N
	j_sun = 0;																					// equatorial sun location

	ik_co2_eff = ik_pole - ik_equator;

	albedo_co2_eff = albedo_pole - albedo_equator;



// effective temperature, albedo and emissivity/absorptivity for the two layer model
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			d_j = ( double ) j;
			if ( h.x[ 0 ][ j ][ k ]  == 0. ) 	albedo.y[ j ][ k ] = albedo_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + albedo_pole;

			if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ 1 ][ j ][ k ] == 0. ) ) 	albedo.y[ j ][ k ] = albedo_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + albedo_pole;

			for ( int i = 1; i < im-1; i++ )
			{
				if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i+1 ][ j ][ k ] == 0. ) ) 	albedo.y[ j ][ k ] = albedo_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + albedo_pole;

				if ( ( t.x[ i ][ j ][ k ] * t_0 - t_0 <= 0. ) && ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i+1 ][ j ][ k ] == 0. ) ) 	albedo.y[ j ][ k ] = albedo_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + albedo_pole;
			}

//				Ik.y[ j ][ k ] = ik * sin ( ( ( double ) j - j_sun ) / pi180 ) * sin ( ( ( double ) k - 90 ) / pi180 ) / 4.;	//	solar short wave radiation on a point location
//				Ik.y[ j ][ k ] = ik * sin ( ( ( double ) j - j_sun ) / pi180 ) / 4.;														// solar short wave radiation as zonal distribution
//				if ( Ik.y[ j ][ k ] < 60. ) Ik.y[ j ][ k ] = 60.;
				Ik.y[ j ][ k ] = ik_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + ik_pole;

				if ( ( Ma == 0 ) && ( h.x[ 0 ][ j ][ k ]  == 1. ) )
				{
					t_Ik = pow ( Ik.y[ j ][ k ] / sigma, ( 1. / 4. ) );
					Ik.y[ j ][ k ] = sigma * pow ( ( t_Ik ), 4. );
				}

			TK = .54;

			Ik.y[ j ][ k ] = TK * Ik.y[ j ][ k ];
		}
	}




// absorption/emissivity computation
//	epsilon_co2_eff_max = .594;															// constant  given by Häckel ( F. Baur and H. Philips, 1934 )

	epsilon_co2_eff_2D = epsilon_pole - epsilon_equator;

	d_i_max = ( double ) ( im - 1 );

	for ( int j = 0; j < jm; j++ )
	{
		d_j = ( double ) j;
		epsilon_co2_eff_max = epsilon_co2_eff_2D * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + epsilon_pole;

		for ( int k = 0; k < km; k++ )
		{


// above assumed tropopause constant temperature t_tropopause
/*
// temperature values for the tropopause
			for ( int i = im_tropopause[ j ] + 1; i < im; i++ )
			{
				t.x[ i ][ j ][ k ] = t_tropopause;
				epsilon_3D.x[ i ][ j ][ k ] = epsilon_tropopause;										// dependency given by Häckel ( F. Baur and H. Philips, 1934 )
			}


// temperature decreasing approaching the tropopause, above constant temperature following Standard Atmosphere
			for ( int i = 0; i <= im_tropopause[ j ]; i++ )
			{

				if ( c.x[ i ][ j ][ k ] <= 0. )	c.x[ i ][ j ][ k ] = 0.;
				if ( cloud.x[ i ][ j ][ k ] <= 0. )	cloud.x[ i ][ j ][ k ] = 0.;
				if ( ice.x[ i ][ j ][ k ] <= 0. )	ice.x[ i ][ j ][ k ] = 0.;

				d_i_max = ( double ) ( im_tropopause[ j ] );
				d_i = ( double ) i;
				e = ( c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k ] ) * p_stat.x[ i ][ j ][ k ] / ep; // COSMO water vapour pressure based on local water vapour, cloud water, cloud ice in hPa

				epsilon_co2_eff = epsilon_co2_eff_max - ( epsilon_tropopause - epsilon_co2_eff_max ) * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) ); // radial distribution
				epsilon_3D.x[ i ][ j ][ k ] = epsilon_co2_eff + .0416 * sqrt ( e );										// dependency given by Häckel ( F. Baur and H. Philips, 1934 )
			}
*/



			for ( int i = 0; i < im; i++ )
			{
				if ( c.x[ i ][ j ][ k ] <= 0. )											c.x[ i ][ j ][ k ] = 0.;
				if ( cloud.x[ i ][ j ][ k ] <= 0. )									cloud.x[ i ][ j ][ k ] = 0.;
				if ( ice.x[ i ][ j ][ k ] <= 0. )										ice.x[ i ][ j ][ k ] = 0.;

				e = ( c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k ] ) * p_stat.x[ i ][ j ][ k ] / ep; // COSMO water vapour pressure based on local water vapour, cloud water, cloud ice in hPa

				d_i = ( double ) i;

				epsilon_co2_eff = epsilon_co2_eff_max - ( epsilon_tropopause - epsilon_co2_eff_max ) * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) ); // radial distribution
				epsilon_3D.x[ i ][ j ][ k ] = epsilon_co2_eff + .0416 * sqrt ( e );										// dependency given by Häckel ( F. Baur and H. Philips, 1934 )


//				if ( ( j == 90 ) && ( k == 180 ) )	cout << " epsilon   i = "<< i << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  c = " << c.x[ i ][ j ][ k ] << "  cl = " << cloud.x[ i ][ j ][ k ] << "  ci = " << ice.x[ i ][ j ][ k ] << "  e = " << e << "  epsilon_co2_eff_max = " << epsilon_co2_eff_max << "  epsilon_co2_eff = " << epsilon_co2_eff << "  epsilon_co2_eff_2D = " << epsilon_co2_eff_2D << "  epsilon_sqrt = " << .0416 * sqrt ( e ) << "  epsilon_3D = " << epsilon_3D.x[ i ][ j ][ k ] << endl;
			}




		}
	}


//	iteration procedure for the computation of the temperature based on the multi-layer radiation model
// temperature needs an initial guess which must be corrected by the long and short wave radiation remaining in the atmosphere
	iter_rad = 0;
	while ( iter_rad <= 5 )																// iter_rad may be varied
	{
		iter_rad = iter_rad + 1;

// coefficient formed for the tridiogonal set of equations for the absorption/emission coefficient of the multi-layer radiation model
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
// radiation boundary conditions for the top ot the troposphere
				radiation_3D.x[ im - 1 ][ j ][ k ] = ( 1. - epsilon_3D.x[ im - 1 ][ j ][ k ] ) * sigma * pow ( t.x[ im - 1 ][ j ][ k ] * t_0, 4. ); // long wave radiation leaving the atmosphere above the tropopause, later needed for non-dimensionalisation

				rad_lon_terrestic = sigma * pow ( t.x[ 0 ][ j ][ k ] * t_0, 4. );									// long wave surface radiation based on local temperature, Stephan-Bolzmann law
				rad_lon_back = epsilon_3D.x[ 1 ][ j ][ k ] * sigma * pow ( t.x[ 1 ][ j ][ k ] * t_0, 4. );	// long wave back radiation absorbed from the first water vapour layer out of 40
 
				Ik_loss = .07 * Ik.y[ j ][ k ] / TK;																				// short wave radiation loss on the surface
				Ik_tot = rad_lon_back + Ik.y[ j ][ k ] - Ik_loss;															// total short and long wave radiation leaving the surface

				AA[ 0 ] = Ik_tot / radiation_3D.x[ im - 1 ][ j ][ k ];													// non-dimensional surface radiation
				CC[ 0 ][ 0 ] = 0.;																									// no absorption of radiation on the surface by water vapour

				radiation_3D.x[ 0 ][ j ][ k ] = ( 1. - epsilon_3D.x[ 0 ][ j ][ k ] ) * sigma * pow ( t.x[ 0 ][ j ][ k ] * t_0, 4. ) / radiation_3D.x[ im - 1 ][ j ][ k ]; // long wave radiation leaving the surface

				for ( int i = 1; i < im - 1; i++ )
				{
					AA[ i ] = AA[ i - 1 ] * ( 1. - epsilon_3D.x[ i ][ j ][ k ] );											// transmitted long wave radiation from each layer
					CC[ i ][ i ]= epsilon_3D.x[ i ][ j ][ k ] * sigma * pow ( t.x[ i ][ j ][ k ] * t_0, 4. ) / radiation_3D.x[ im - 1 ][ j ][ k ]; // absorbed long wave radiation in each layer

					radiation_3D.x[ i ][ j ][ k ] = ( 1. - epsilon_3D.x[ i ][ j ][ k ] ) * sigma * pow ( t.x[ i ][ j ][ k ] * t_0, 4. ) / radiation_3D.x[ im - 1 ][ j ][ k ]; // long wave radiation leaving each layer

					for ( int l = i + 1; l < im; l++ )
					{
						CC[ i ][ l ]= CC[ i ][ l - 1 ] * ( 1. - epsilon_3D.x[ l ][ j ][ k ] );								// additional transmitted radiation from layer to layer in radial direction
					}
				}


// Thomas algorithm to solve the tridiogonal equation system for the solution of the radiation with a recurrence formula
// additionally embedded in an iterational process
				for ( int i = 0; i < im-1; i++ )
				{
// values at the surface
					if ( i == 0 )
					{
						aa = 0.;
						bb = - radiation_3D.x[ 0 ][ j ][ k ];
						cc = radiation_3D.x[ 1 ][ j ][ k ];
						dd = - AA[ 0 ];
						CCC = 0.;
						DDD = 0.;
						alfa[ 0 ] = 0.;
						beta[ 0 ] = 0.;
					}

// values in the field
					if ( i == 1 )
					{
						aa = radiation_3D.x[ 0 ][ j ][ k ];
						bb = - 2. * radiation_3D.x[ 1 ][ j ][ k ];
						cc = radiation_3D.x[ 2 ][ j ][ k ];
						dd = - AA[ 0 ] + AA[ 1 ];
						CCC = 0.;
						DDD = 0.;
					}

					if ( i == 2 )
					{
						CCC = CC[ 1 ][ 2 ];
						aa = radiation_3D.x[ 1 ][ j ][ k ];
						bb = - 2. * radiation_3D.x[ 2 ][ j ][ k ];
						cc = radiation_3D.x[ 3 ][ j ][ k ];
						dd = - AA[ 1 ] + AA[ 2 ] + CCC;
					}

					if ( i > 2 )
					{
						CCC = 0.;
						DDD = 0.;

						for ( int l = 1; l <= i - 1; l++ )
						{
							CCC = CCC + CC[ l ][ i ];
						}

						for ( int l = 1; l <= i - 2; l++ )
						{
							DDD = DDD + CC[ l ][ i - 1 ];
						}

						aa = radiation_3D.x[ i - 1 ][ j ][ k ];
						bb = - 2. * radiation_3D.x[ i ][ j ][ k ];
						if ( i == im - 2 )  cc = radiation_3D.x[ i + 1 ][ j ][ k ] / radiation_3D.x[ im - 1 ][ j ][ k ];
						else cc = radiation_3D.x[ i + 1 ][ j ][ k ];
						dd = - AA[ i - 1 ] + AA[ i ] + CCC - DDD;
					}

					if ( i > 0 )
					{
						alfa[ i ] = cc / ( bb - aa * alfa[ i - 1 ] );
						beta[ i ] = ( dd - aa * beta[ i - 1 ] ) / ( bb - aa * alfa[ i - 1 ] );
					}
					else
					{
						alfa[ i ] = cc / bb;
						beta[ i ] = dd / bb;
					}

//					if ( ( j == 90 ) && ( k == 180 ) )	cout << " rad berechnung  i = "<< i << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  radiation_3D = " << radiation_3D.x[ i ][ j ][ k ] << "  p_stat = " << p_stat.x[ i ][ j ][ k ] << endl;

				}
 

				radiation_3D.x[ im - 1 ][ j ][ k ] = ( 1. - epsilon_3D.x[ im - 1 ][ j ][ k ] ) * sigma * pow ( t.x[ im - 1 ][ j ][ k ] * t_0, 4. ); // dimensional form of the radiation leaving the last layer




// recurrence formula for the radiation and temperature

				for ( int i = im - 2; i >= 0; i-- )
				{

// above assumed tropopause constant temperature t_tropopause
/*
					if ( i >= im_tropopause[ j ] )
					{
						radiation_3D.x[ i + 1 ][ j ][ k ] = radiation_3D.x[ im - 1 ][ j ][ k ];
					}
*/
					radiation_3D.x[ i ][ j ][ k ] = - alfa[ i ] * radiation_3D.x[ i + 1 ][ j ][ k ] + beta[ i ];							// Thomas algorithm, recurrence formula
					t.x[ i ][ j ][ k ] = .5 * ( t.x[ i ][ j ][ k ] + pow ( radiation_3D.x[ i ][ j ][ k ] / sigma, ( 1. / 4. ) ) / t_0 );	// averaging of temperature values to smooth the iterations

//					if ( ( j == 90 ) && ( k == 180 ) )	cout << " multi  i = "<< i << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  radiation_3D = " << radiation_3D.x[ i ][ j ][ k ] << "  p_stat = " << p_stat.x[ i ][ j ][ k ] << endl;

				}
			}
		}
	}

}






void BC_Thermo::BC_Temperature ( int *im_tropopause, Array_2D &temperature_NASA, Array &h, Array &t, Array &p_dyn, Array &p_stat )
{
// boundary condition of  temperature on land 
// parabolic distribution from pole to pole accepted
// temperature on land at equator t_max = 1.055 compares to 15° C compared to 288 K
// temperature at tropopause t_min = 0.77 compares to -62° C compares to 211 K
// temperature at tropopause t_min = 0.89 compares to -30° C compares to 243 K
// temperature difference from equator to pole   18°C compares to  t_delta = 0.0659  compares to  18 K
	j_half = ( jm -1 ) / 2;
	j_max = jm - 1;

	d_i_max = ( double ) i_max;
	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

// temperature-distribution by Ruddiman approximated by a parabola
	t_co2_eff = t_pole - t_equator;
	t_cretaceous_co2_eff = t_cretaceous_max / ( ( double ) Ma_max_half - ( double ) ( Ma_max_half * Ma_max_half / Ma_max ) );   // in °C
	t_cretaceous = t_cretaceous_co2_eff * ( double ) ( - ( Ma * Ma ) / Ma_max + Ma );   // in °C
	if ( Ma == 0 ) 	t_cretaceous = t_cret_cor = 0.;

	cout.precision ( 3 );

	time_slice_comment = "      time slice of Cretaceous-AGCM:";
	time_slice_number = " Ma = ";
	time_slice_unit = " million years";

	cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) << time_slice_comment << resetiosflags ( ios::left ) << setw ( 6 ) << fixed << setfill ( ' ' ) << time_slice_number << setw ( 3 ) << Ma << setw ( 12 ) << time_slice_unit << endl << endl;

	temperature_comment = "      temperature increase at cretaceous times: ";
	temperature_gain = " t increase";
	temperature_modern = "      mean temperature at modern times: ";
	temperature_cretaceous = "      mean temperature at cretaceous times: ";
	temperature_average = " t modern";
	temperature_average_cret = " t cretaceous";
	temperature_unit =  "°C ";

	cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) << temperature_comment << resetiosflags ( ios::left ) << setw ( 12 ) << temperature_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << t_cretaceous << setw ( 5 ) << temperature_unit << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_modern << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average << setw ( 5 ) << temperature_unit << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_cretaceous << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average_cret  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average + t_cretaceous << setw ( 5 ) << temperature_unit << endl;

	t_cretaceous = ( t_cretaceous + t_average + t_0 ) / t_0 - ( ( t_average + t_0 ) / t_0 );    // non-dimensional

	// temperatur distribution at aa prescribed sun position
	// sun_position_lat = 60,    position of sun j = 120 means 30°S, j = 60 means 30°N
	// sun_position_lon = 180, position of sun k = 180 means 0° or 180° E ( Greenwich, zero meridian )
	// asymmetric temperature distribution from pole to pole for  j_d  maximum temperature ( linear equation + parabola )

		if ( ( Ma > 0 ) && ( sun == 1 ) )
		{
			j_par = sun_position_lat;																	// position of maximum temperature, sun position
			j_par = j_par + declination;																// angle of sun axis, declination = 23,4°
			j_pol = jm - 1;

			j_par_f = ( double ) j_par;
			j_pol_f = ( double ) j_pol;

			aa = ( t_equator - t_pole ) / ( ( ( j_par_f * j_par_f ) - ( j_pol_f * j_pol_f ) ) - 2. * j_par_f * ( j_par_f - j_pol_f ) );
			bb = - 2. * aa * j_par_f;
			cc = t_equator + aa * j_par_f * j_par_f;
			j_d = sqrt ( ( cc - t_pole ) / aa );
			dd = 2. * aa * j_d + bb;
			t_dd = dd * j_d + t_pole;
			e = t_pole;


	// asymmetric temperature distribution from pole to pole for  j_d  maximum temperature ( linear equation + parabola )
			for ( int k = 0; k < km; k++ )
			{
				for ( int j = 0; j < jm; j++ )
				{
					d_j = ( double ) j;
					if ( d_j <= j_d )
					{
						t.x[ 0 ][ j ][ k ] = dd * d_j + e + t_cretaceous;
					}
					if ( d_j > j_d )
					{
						t.x[ 0 ][ j ][ k ] = aa * d_j * d_j + bb * d_j + cc + t_cretaceous;
					}
				}
			}


	// transfer of zonal constant temperature into aa 1D-temperature field
			for ( int j = 0; j < jm; j++ )
			{
				jm_temp_asym[ j ] = t.x[ 0 ][ j ][ 20 ];
			}


	// longitudinally variable temperature distribution from west to east in parabolic form
	// communicates the impression of local sun radiation on the southern hemisphere
			k_par = sun_position_lon;												// position of the sun at constant longitude
			k_pol = km - 1;

			t_360 = (  t_0 + 5. ) / t_0;

			for ( int j = 0; j < jm; j++ )
			{
				for ( int k = 0; k < km; k++ )
				{
					k_par_f = ( double ) k_par;
					k_pol_f = ( double ) k_pol;
					d_k = ( double ) k;

					aa = ( jm_temp_asym[ j ] - t_360 ) / ( ( ( k_par_f * k_par_f ) - ( k_pol_f * k_pol_f ) ) - 2. * k_par_f * ( k_par_f - k_pol_f ) );
					bb = - 2. * aa * k_par_f;
					cc = jm_temp_asym[ j ] + aa * k_par_f * k_par_f;

					t.x[ 0 ][ j ][ k ] = aa * d_k * d_k + bb * d_k + cc;
				}
			}
		}																								// temperatur distribution at aa prescribed sun position




//      cout << endl << "   NASATemperature = " << NASATemperature << "   t_pole = " << t_pole << "   t_cretaceous = " << t_cretaceous << "   t_land = " << t_land << "   t_co2_eff = " << t_co2_eff << "   d_j_half = " << d_j_half << "   sun = " << sun << endl;


	if ( RadiationModel == 3 )
	{		if ( sun == 0 )
		{
			for ( int k = 0; k < km; k++ )
			{
				for ( int j = 0; j < jm; j++ )
				{
					d_j = ( double ) j;
					if ( NASATemperature == 0 )									// if ( NASATemperature == 0 ) parabolic surface temperature is used
					{
						t.x[ 0 ][ j ][ k ] = t_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + t_pole + t_cretaceous;
						if ( h.x[ 0 ][ j ][ k ] == 1. ) 									t.x[ 0 ][ j ][ k ] = t_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + t_pole + t_cretaceous + t_land;	// end parabolic temperature distribution
					}

					if ( NASATemperature == 1 )									// if ( NASATemperature == 1 ) surface temperature is NASA based
					{
						t.x[ 0 ][ j ][ k ] = t.x[ 0 ][ j ][ k ] + ( t_cretaceous - t_cret_cor );
						if ( h.x[ 0 ][ j ][ k ] == 1. ) 									t.x[ 0 ][ j ][ k ] = t_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + t_pole + t_cretaceous + t_land;	// end parabolic temperature distribution
					}
				}
			}
		}
	}


// temperature approaching the tropopause, above constant temperature following Standard Atmosphere

	d_i_max = ( double ) ( im - 1 );

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 1; i <= im - 1; i++ )
			{
				hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
				d_i = ( double ) i;
				t.x[ i ][ j ][ k ] = ( t_tropopause - t.x[ 0 ][ j ][ k ] ) / d_i_max * d_i + t.x[ 0 ][ j ][ k ];				// linear temperature decay up to tropopause, privat approximation
//				t.x[ i ][ j ][ k ] = ( t.x[ 0 ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) / t_0;		// linear temperature decay up to tropopause, Standard Atmosphere with slope 0.0065°C/m, produces too low temperatures
			}
		}
	}




// above assumed tropopause constant temperature t_tropopause
/*
// temperature decreasing approaching the tropopause, above constant temperature following Standard Atmosphere
	for ( int j = 0; j < jm; j++ )
	{
		d_i_max = ( double ) im_tropopause[ j ];

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 1; i <= im_tropopause[ j ]; i++ )
			{
				hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
				d_i = ( double ) i;
//				t.x[ i ][ j ][ k ] = ( t_tropopause - t.x[ 0 ][ j ][ k ] ) / d_i_max * d_i + t.x[ 0 ][ j ][ k ];				// linear temperature decay up to tropopause, privat approximation
				t.x[ i ][ j ][ k ] = ( t.x[ 0 ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) / t_0;																// linear temperature decay up to tropopause, Standard Atmosphere with slope 0.0065°C/m
			}
		}
	}


// temperature values for the tropopause
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = im_tropopause[ j ] + 1; i < im; i++ )
			{
				t.x[ i ][ j ][ k ] = t_tropopause;
			}
		}
	}
*/

}







void BC_Thermo::BC_WaterVapour ( Array &h, Array &t, Array &c )
{
// initial and boundary conditions of water vapour on water and land surfaces
// parabolic water vapour distribution from pole to pole accepted

// maximum water vapour content on water surface at equator c_equator = 1.04 compares to 0.04 volume parts
// minimum water vapour at tropopause c_tropopause = 0.0 compares to 0.0 volume parts
// value 0.04 stands for the maximum value of 40 g/kg, g water vapour per kg dry air

	d_i_max = ( double ) i_max;
	j_half = ( jm -1 ) / 2;

	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

//	c_coeff = - 1.;

// water vapour contents computed by Clausius-Clapeyron-formula
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			d_j = ( double ) j;
			if ( h.x[ 0 ][ j ][ k ]  == 0. ) 
			{
				c.x[ 0 ][ j ][ k ]  = hp * ep *exp ( 17.0809 * ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) ) / ( ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01 );	// saturation of relative water vapour in kg/kg
				c.x[ 0 ][ j ][ k ] = c_ocean * c.x[ 0 ][ j ][ k ];														// relativ water vapour contents on ocean surface reduced by factor
			}

			if ( h.x[ 0 ][ j ][ k ]  == 1. ) 
			{
				c.x[ 0 ][ j ][ k ]  = hp * ep * exp ( 17.0809 * ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) ) / ( ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01 );
				c.x[ 0 ][ j ][ k ] = c_land * c.x[ 0 ][ j ][ k ];															// relativ water vapour contents on land reduced by factor
			}

//			c.x[ 0 ][ j ][ k ] = - c_coeff * c.x[ 0 ][ j ][ k ] * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half );

		}
	}



// water vapour distribution decreasing approaching tropopause
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 1; i < im; i++ )
			{
				d_i_max = ( double ) ( im - 1 );
				d_i = ( double ) i;

				c.x[ i ][ j ][ k ] = c.x[ 0 ][ j ][ k ] - ( c_tropopause - c.x[ 0 ][ j ][ k ] ) * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) ); 		// parabolic decrease
//				c.x[ i ][ j ][ k ] = ( c_tropopause - c.x[ 0 ][ j ][ k ] ) * d_i / d_i_max + c.x[ 0 ][ j ][ k ];												// linear decrease

			}																															// end i
		}																																// end j
	}																																	// end k
}







void BC_Thermo::BC_CloudWaterIce ( Array &h, Array &t, Array &cloud, Array &ice )
{
// initial and boundary conditions of cloud water and cloud ice in vertical direction

// initial vertical cloud ice distribution


	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				cloud.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ] = 0.;
			}
		}
	}



//	q_m = .002;
	q_m = .0005;
	r_m = 10;
	r_a = 5;
	r_b = 15;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = r_a; i <= r_m; i++ )
			{
				r_m_i = ( double ) i / ( double ) r_m ;
				cloud.x[ i ][ j ][ k ] = q_m * r_m_i * r_m_i * ( 3. - 2. * r_m_i );
			}																															// end i

			int l = 1;
			for ( int i = r_m + 1; i <= r_b; i++ )
			{
				cloud.x[ i ][ j ][ k ] = cloud.x[ r_m - l ][ j ][ k ];
				l++;
			}																															// end i

		}																																// end j
	}																																	// end k


// initial vertical cloud ice distribution
//	q_m = .001;
	q_m = .0001;
	r_m = 20;
	r_a = 15;
	r_b = 25;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = r_a; i <= r_m; i++ )
			{
				r_m_i = ( double ) i / ( double ) r_m ;
				ice.x[ i ][ j ][ k ] = q_m * r_m_i * r_m_i * ( 3. - 2. * r_m_i );
			}																															// end i

			int l = 1;
			for ( int i = r_m + 1; i <= r_b; i++ )
			{
				ice.x[ i ][ j ][ k ] = ice.x[ r_m - l ][ j ][ k ];
				l++;
			}																															// end i

		}																																// end j
	}																																	// end k

}







void BC_Thermo::BC_CO2 ( Array_2D &Vegetation, Array &h, Array &t, Array &p_dyn, Array &co2 )
{
// initial and boundary conditions of CO2 content on water and land surfaces
// parabolic CO2 content distribution from pole to pole accepted

	j_half = j_max / 2;

	t_cretaceous_co2_eff = t_cretaceous_max / ( ( double ) Ma_max_half - ( double ) ( Ma_max_half * Ma_max_half / Ma_max ) );   // in °C
	t_cretaceous = t_cretaceous_co2_eff * ( double ) ( - ( Ma * Ma ) / Ma_max + Ma );   // in °C
	if ( Ma == 0 ) 	t_cretaceous = 0.;

// CO2-distribution by Ruddiman approximated by a parabola
	co2_cretaceous = 3.2886 * pow ( ( t_cretaceous + t_average ), 2 ) - 32.8859 * ( t_cretaceous + t_average ) + 102.2148;  // in ppm
	co2_average = 3.2886 * pow ( t_average, 2 ) - 32.8859 * t_average + 102.2148;  // in ppm
	co2_cretaceous = co2_cretaceous - co2_average;
	if ( Ma == 0 ) 	co2_cretaceous = 0.;

	cout.precision ( 3 );

	co_comment = "      co2 increase at cretaceous times: ";
	co_gain = " co2 increase";
	co_modern = "      mean co2 at modern times: ";
	co_cretaceous_str = "      mean co2 at cretaceous times: ";
	co_average_str = " co2 modern";
	co_average_cret = " co2 cretaceous";
	co_unit =  "ppm ";

	cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) << co_comment << resetiosflags ( ios::left ) << setw ( 12 ) << co_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << co2_cretaceous << setw ( 5 ) << co_unit << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) << co_modern << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_str  << " = "  << setw ( 7 )  << setfill ( ' ' ) << co2_average << setw ( 5 ) << co_unit << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) << co_cretaceous_str << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_cret  << " = "  << setw ( 7 )  << setfill ( ' ' ) << co2_average + co2_cretaceous << setw ( 5 ) << co_unit << endl;
	cout << endl;

	d_i_max = ( double ) i_max;
	d_j_half = ( double ) j_half;

	co_co2_eff = co2_pole - co2_average;

// CO2-content as initial solution
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( h.x[ 0 ][ j ][ k ]  == 0. ) 
			{
				d_j = ( double ) j;
				co2.x[ 0 ][ j ][ k ] = ( co_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + co2_pole + co2_cretaceous + co2_ocean ) / co2_0; // non-dimensional
			}
			if ( h.x[ 0 ][ j ][ k ]  == 1. ) 
			{
				d_j = ( double ) j;
				co2.x[ 0 ][ j ][ k ] =  ( co_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + co2_pole + co2_cretaceous + co2_land - co2_vegetation * Vegetation.y[ j ][ k ] ) / co2_0;	// parabolic distribution from pole to pole
			}
		}
	}




// co2 distribution decreasing approaching tropopause, above no co2
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 1; i < im; i++ )
			{
					d_i = ( double ) i;
					co2.x[ i ][ j ][ k ] = co2.x[ 0 ][ j ][ k ] - ( ( co2_tropopause - co2.x[ 0 ][ j ][ k ] * co2_0 ) * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) ) ) / co2_0;
																																										// radial distribution approximated by a parabola
//					co2.x[ i ][ j ][ k ] = ( co2_tropopause - co2.x[ 0 ][ j ][ k ] ) / d_i_max * d_i + co2.x[ 0 ][ j ][ k ];			// linear co2 decay up to tropopause
			}
		}
	}

}





void BC_Thermo::TropopauseLocation ( int *im_tropopause )
{
// parabolic tropopause location distribution from pole to pole assumed

	j_half = ( jm -1 ) / 2;

	d_j_half = ( double ) j_half;

	trop_co2_eff = ( double ) ( i_beg - i_max );

// computation of the tropopause from pole to pole

	for ( int j = 0; j < jm; j++ )
	{
		d_j = ( double ) j;
		im_tropopause[ j ] = ( trop_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) ) + i_beg;
	}
}





void BC_Thermo::BC_Radiation_parabolic ( Array_2D &Radiation_Balance_par, Array &h )
{
// initial and boundary conditions for the radiation balance on water and land surfaces
// parabolic distribution from pole to pole accepted

	j_half = j_max / 2;

	rad_bal_minus = epsilon_extra * sigma * pow ( t_cretaceous * t_0, 4. );	// decrease of radiation during paleo climate

	radiation_land_co2_eff = radiation_pole - radiation_equator;
	radiation_ocean_co2_eff = radiation_pole - radiation_equator - radiation_ocean; // increase by 40 W/m² on ocean


// radiation-distribution as initial solution
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( h.x[ 0 ][ j ][ k ]  == 0. ) 
			{
				d_j = ( double ) j;
				Radiation_Balance_par.y[ j ][ k ] = radiation_ocean_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + radiation_pole;
				Radiation_Balance_par.y[ j ][ k ] = Radiation_Balance_par.y[ j ][ k ] - rad_bal_minus;
			}
			if ( h.x[ 0 ][ j ][ k ]  == 1. ) 
			{
				d_j = ( double ) j;
				Radiation_Balance_par.y[ j ][ k ] = radiation_land_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + radiation_pole; // parabolic distribution from pole to pole
				Radiation_Balance_par.y[ j ][ k ] = Radiation_Balance_par.y[ j ][ k ] - rad_bal_minus;
			}
		}
	}

}






void BC_Thermo::IC_CellStructure ( int *im_tropopause, Array &h, Array &u, Array &v, Array &w )
{
// boundary condition for the velocity components in the circulation cells

// latest version by Grotjahn ( Global Atmospheric Circulations, 1993 )
// default for the velocity components u, v, and w as initial conditions


// velocities given in m/s, 1 m/s compares to 3.6 km/h, non-dimensionalized by u_0 below
// do not change the velocity initial conditions !!

// velocity assumptions at the equator 0°
	ua_00 = 0.1;																				// in m/s compares to 1.08 km/h, non-dimensionalized by u_0 below

	va_equator_SL =  0.000;
	va_equator_Tropopause = 0.000;

	wa_equator_SL = - 5.;
	wa_equator_Tropopause = - 7.5;

// velocity assumptions for latitude at 15° and 30° in the Hadley cell
	ua_30 = - 0.1;

	va_Hadley_SL = 1.;
	va_Hadley_Tropopause = - 1.;

	va_Hadley_SL_15 = 1.;
	va_Hadley_Tropopause_15 = - 1.;

	wa_Hadley_SL = 5.;																	// at surface
	wa_Hadley_Tropopause = 30.;													// subtropic jet in m/s compares to 108 km/h

// velocity assumptions for latitude at 45° and 60° in the Ferrel cell
	ua_60 = 0.1;

	va_Ferrel_SL = - 0.75;
	va_Ferrel_Tropopause = 1.;

	va_Ferrel_SL_45 = - 1.;
	va_Ferrel_Tropopause_45 = 1.;

	wa_Ferrel_SL = 1.5;																	// subpolar jet
	wa_Ferrel_Tropopause = 10.;														// subpolar jet in m/s compares to 36 km/h

// velocity assumptions for latitude 90° in the Polar cell
	ua_90 = - 0.1;

	va_Polar_SL = 0.;
	va_Polar_Tropopause = 0.;

	va_Polar_SL_75 = 0.75;
	va_Polar_Tropopause_75 = - 0.75;

	va_Polar_SL_75 = 0.75;
	va_Polar_Tropopause_75 = - 0.75;

/*
// lower velocity assumtions for test case
	ua_00 = 0.1;																				// in m/s compares to 1.08 km/h, non-dimensionalized by u_0 below

	va_equator_SL =  0.000;
	va_equator_Tropopause = 0.000;

	wa_equator_SL = - .1;
	wa_equator_Tropopause = - .15;


	ua_30 = - 0.1;

	va_Hadley_SL = .1;
	va_Hadley_Tropopause = - .1;

	va_Hadley_SL_15 = .1;
	va_Hadley_Tropopause_15 = - .1;

	wa_Hadley_SL = 1.;																	// at surface
	wa_Hadley_Tropopause = .2;													// subtropic jet in m/s compares to 108 km/h


	ua_60 = 0.1;

	va_Ferrel_SL = - 0.075;
	va_Ferrel_Tropopause = .1;

	va_Ferrel_SL_45 = - .1;
	va_Ferrel_Tropopause_45 = .1;

	wa_Ferrel_SL = .15;																	// subpolar jet
	wa_Ferrel_Tropopause = 1.;														// subpolar jet in m/s compares to 36 km/h



	ua_90 = - 0.1;

	va_Polar_SL_75 = 0.075;
	va_Polar_Tropopause_75 = - 0.075;

	va_Polar_SL_75 = 0.075;
	va_Polar_Tropopause_75 = - 0.075;

	wa_Polar_SL = 0.;
	wa_Polar_Tropopause = 0.;
*/


// preparations for diagonal velocity value connections
	im_1 = im - 1;

	j_aeq = 90;

	j_pol_n = 0;
	j_pol_s = jm-1;
	j_pol_v_n = 15;
	j_pol_v_s = 165;

	j_fer_n = 30;
	j_fer_s = 150;
	j_fer_v_n = 45;
	j_fer_v_s = 135;

	j_had_n = 60;
	j_had_s = 120;
	j_had_v_n = 75;
	j_had_v_s = 105;



/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// initial velocity components in the northern and southern
// Pole, Ferrel and Hadley cells


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////// equator ///////////////////////////////////////


// equator ( at j=90 compares to 0° latitude )
// u-component up to tropopause and back on half distance (  i = 20 )
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_aeq; j < j_aeq + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			i_half = im_tropopause[ j ] / 2;
			d_i_half = ( double ) i_half;

			for ( int i = 0; i <= i_max; i++ )
			{
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = - ua_00 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


// equator ( at j=90 )
// v- and w-component up to tropopause and stratosphere above
	for ( int j = j_aeq; j < j_aeq + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) i_max;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_equator_Tropopause - va_equator_SL ) * d_i / d_i_max + va_equator_SL;
				w.x[ i ][ j ][ k ] = ( wa_equator_Tropopause - wa_equator_SL ) * d_i / d_i_max + wa_equator_SL;
			}
		}
	}


	for ( int j = j_aeq; j < j_aeq + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) i_max - ( double ) im_1;
		if ( i_max == 40 ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) i - ( double ) im_1;
				v.x[ i ][ j ][ k ] = va_equator_Tropopause * d_i / d_i_max;
				w.x[ i ][ j ][ k ] = wa_equator_Tropopause * d_i / d_i_max;
			}
		}
	}


/////////////////////////////////////// end equator ///////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////

// cell structure in northern hemisphere

////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////// northern polar cell /////////////////////////////////////////


// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// u-component up to tropopause and back on half distance
// extension around North Pole ( from j=0 till j=5 compares to 85° till 90° northern latitude )
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < j_pol_n +1; j++ )
		{
			i_max = im_tropopause[ j ];
			i_half = im_tropopause[ j ] / 2;
			d_i_half = ( double ) i_half;

			for ( int i = 0; i <= i_max; i++ )
			{
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = - ua_90 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// v- and w-component from Pole up to tropopause
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_n; j < j_fer_n + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			d_i_max = ( double ) i_max;

			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Polar_Tropopause - va_Polar_SL ) * d_i / d_i_max + va_Polar_SL;
				w.x[ i ][ j ][ k ] = ( wa_Polar_Tropopause - wa_Polar_SL ) * d_i / d_i_max + wa_Polar_SL;   // indifferent except at j_pol_n
			}
		}
	}

	for ( int j = j_pol_n; j < j_fer_n + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Polar_Tropopause * d_i / d_i_max;   // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = wa_Polar_Tropopause * d_i / d_i_max;   // replacement for forming diagonals
			}
		}
	}


// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// v-component up to tropopause and back on half distance
// extension around the North Pole ( from j=0 till j=5 compares to 85° till 90° northern latitude )
// v-component at 75°N
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_v_n; j <  j_pol_v_n + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			d_i_max = ( double ) i_max;

			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Polar_Tropopause_75 - va_Polar_SL_75 ) * d_i / d_i_max + va_Polar_SL_75;
			}
		}
	}

	for ( int j = j_pol_v_n; j < j_pol_v_n + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Polar_Tropopause_75 * d_i / d_i_max;
			}
		}
	}


/////////////////////////////////// end northern polar cell /////////////////////////////////////////



/////////////////////////////////// northern Ferrel cell /////////////////////////////////////////


// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° northern latitude )

// u-component up to tropopause and back on half distance
// extension around 60° northern latitude ( from j=29 till j=31 compares to 59° till 61° northern latitude )
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_n; j < j_fer_n + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			i_half = im_tropopause[ j ] / 2;
			d_i_half = ( double ) i_half;

			for ( int i = 0; i <= i_max; i++ )
			{
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = - ua_60 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


// north equatorial Ferrel cell
// v- and w-component up to tropopause and stratosphere above
// extension around 60° northern latitude ( from j=29 till j=31 compares to 59° till 61° northern latitude )

	for ( int j = j_fer_n; j < j_fer_n + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) i_max;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Ferrel_Tropopause - va_Ferrel_SL ) * d_i / d_i_max + va_Ferrel_SL;   // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = ( wa_Ferrel_Tropopause - wa_Ferrel_SL ) * d_i / d_i_max + wa_Ferrel_SL;   // replacement for forming diagonals
			}
		}
	}

	for ( int j = j_fer_n; j < j_fer_n + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Ferrel_Tropopause * d_i / d_i_max;   // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = wa_Ferrel_Tropopause * d_i / d_i_max;   // replacement for forming diagonals
			}
		}
	}


// north equatorial Ferrel cell
// v-component up to tropopause and stratosphere above
// extension around 60° northern latitude ( from j=29 till j=31 compares to 59° till 61° northern latitude )
// v-component at 45°N
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_v_n; j <  j_fer_v_n + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			d_i_max = ( double ) i_max;

			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Ferrel_Tropopause_45 - va_Ferrel_SL_45 ) * d_i / d_i_max + va_Ferrel_SL_45;
			}
		}
	}

	for ( int j = j_fer_v_n; j < j_fer_v_n + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Ferrel_Tropopause_45 * d_i / d_i_max;
			}
		}
	}


//////////////////////////////// end northern Ferrel cell ////////////////////////////////////////



//////////////////////////////// northern Hadley cell ////////////////////////////////////////


// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// u-component up to tropopause and back on half distance
// extension around 30° northern latitude ( from j=59 till j=61 compares to 29° till 31° northern latitude )
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_n; j < j_had_n + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			i_half = im_tropopause[ j ] / 2;
			d_i_half = ( double ) i_half;

			for ( int i = 0; i <= i_max; i++ )
			{
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = - ua_30 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


// north equatorial Hadley cell
// v- and w-component up to tropopause and stratosphere above
// extension around 30° northern latitude ( from j=59 till j=61 compares to 29° till 31° northern latitude )
	for ( int j = j_had_n; j < j_had_n + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) i_max;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Hadley_Tropopause - va_Hadley_SL ) * d_i / d_i_max + va_Hadley_SL;    // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = ( wa_Hadley_Tropopause - wa_Hadley_SL ) * d_i / d_i_max + wa_Hadley_SL;    // replacement for forming diagonals
			}
		}
	}


	for ( int j = j_had_n; j < j_had_n + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Hadley_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = wa_Hadley_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
			}
		}
	}


// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// v-component up to tropopause and stratosphere above
// v-component at 15°N
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_v_n; j <  j_had_v_n + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			d_i_max = ( double ) i_max;

			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Hadley_Tropopause_15 - va_Hadley_SL_15 ) * d_i / d_i_max + va_Hadley_SL_15;
			}
		}
	}

	for ( int j = j_had_v_n; j < j_had_v_n + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Hadley_Tropopause_15 * d_i / d_i_max;
			}
		}
	}


////////////////////////////////// end northern Hadley cell ////////////////////////////////////////



// ||||||||||||||||||||||||||||||||||||||| equator ||||||||||||||||||||||||||||||||||||||||||||||
// ||||||||||||||||||||||||||||||||||||||| equator ||||||||||||||||||||||||||||||||||||||||||||||



////////////////////////////////////////////////////////////////////////////////////////////////

// cell structure in southern hemisphere

////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////// southern Hadley cell ////////////////////////////////////////



// south equatorial Hadley cell ( from j=90 till j=120 compares to 0° till 30° southern latitude )
// u-component up to tropopause and back on half distance
// extension around 30° southern latitude ( from j=119 till j=121 compares to 29° till 31° southern latitude )
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_s; j < j_had_s + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			i_half = im_tropopause[ j ] / 2;
			d_i_half = ( double ) i_half;

			for ( int i = 0; i <= i_max; i++ )
			{
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = - ua_30 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


// south equatorial Hadley cell
// v- and w-component up to tropopause and stratosphere above
// extension around 30° southern latitude ( from j=119 till j=121 compares to 29° till 31° southern latitude )
	for ( int j = j_had_s; j < j_had_s + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) i_max;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Hadley_Tropopause - va_Hadley_SL ) * d_i / d_i_max + va_Hadley_SL;    // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = ( wa_Hadley_Tropopause - wa_Hadley_SL ) * d_i / d_i_max + wa_Hadley_SL;    // replacement for forming diagonals
			}
		}
	}

	for ( int j = j_had_s; j < j_had_s + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Hadley_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = wa_Hadley_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
			}
		}
	}

// extension around 30° southern latitude ( from j=119 till j=121 compares to 29° till 31° southern latitude )
// v-component up to tropopause and stratosphere above
// v-component at 15°S
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_v_s; j <  j_had_v_s + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			d_i_max = ( double ) i_max;

			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Hadley_Tropopause_15 - va_Hadley_SL_15 ) * d_i / d_i_max + va_Hadley_SL_15;
			}
		}
	}

	for ( int j = j_had_v_s; j < j_had_v_s + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Hadley_Tropopause_15 * d_i / d_i_max;
			}
		}
	}


////////////////////////////// end southern Hadley cell ////////////////////////////////////////////



////////////////////////////// southern Ferrel cell ////////////////////////////////////////////


// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// u-component up to tropopause and back on half distance
// extension around 60° southern latitude ( from j=139 till j=151 compares to 59° till 61° southern latitude )
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_s; j < j_fer_s + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			i_half = im_tropopause[ j ] / 2;
			d_i_half = ( double ) i_half;

			for ( int i = 0; i <= i_max; i++ )
			{
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = - ua_60 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


// south equatorial Ferrel cell
// v- and w-component up to tropopause and stratosphere above
	for ( int j = j_fer_s; j < j_fer_s + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) i_max;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Ferrel_Tropopause - va_Ferrel_SL ) * d_i / d_i_max + va_Ferrel_SL;    // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = ( wa_Ferrel_Tropopause - wa_Ferrel_SL ) * d_i / d_i_max + wa_Ferrel_SL;    // replacement for forming diagonals
			}
		}
	}

	for ( int j = j_fer_s; j < j_fer_s + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Ferrel_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = wa_Ferrel_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
			}
		}
	}


// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// v-component up to tropopause and stratosphere above
// v-component at 45°N
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_v_s; j <  j_fer_v_s + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			d_i_max = ( double ) i_max;

			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Ferrel_Tropopause_45 - va_Ferrel_SL_45 ) * d_i / d_i_max + va_Ferrel_SL_45;
			}
		}
	}

	for ( int j = j_fer_v_s; j < j_fer_v_s + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Ferrel_Tropopause_45 * d_i / d_i_max;
			}
		}
	}


///////////////////////////// end southern Ferrel cell /////////////////////////////////////////////



///////////////////////////// southern Polar cell /////////////////////////////////////////////


// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// u-component up to tropopause and back on half distance
// extension around South Pole ( from j=175 till j=180 compares to 85° till 90° southern latitude )
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_s; j < jm; j++ )
		{
			i_max = im_tropopause[ j ];
			i_half = im_tropopause[ j ] / 2;
			d_i_half = ( double ) i_half;

			for ( int i = 0; i <= i_max; i++ )
			{
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = - ua_90 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_s + 1; j < j_pol_s + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			d_i_max = ( double ) i_max;
			
			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( ( va_Polar_Tropopause - va_Polar_SL ) * d_i / d_i_max + va_Polar_SL );
				w.x[ i ][ j ][ k ] = ( wa_Polar_Tropopause - wa_Polar_SL ) * d_i / d_i_max + wa_Polar_SL;  // indifferent except at j_pol_s
			}
		}
	}

	for ( int j = j_fer_s + 1; j < j_pol_s + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Polar_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
				w.x[ i ][ j ][ k ] = wa_Polar_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
			}
		}
	}


// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// u-component up to tropopause and back on half distance
// extension around South Pole ( from j=175 till j=180 compares to 85° till 90° southern latitude )
// v-component at 75°S
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_v_s; j <  j_pol_v_s + 1; j++ )
		{
			i_max = im_tropopause[ j ];
			d_i_max = ( double ) i_max;

			for ( int i = 0; i < i_max; i++ )
			{
				d_i = ( double ) i;
				v.x[ i ][ j ][ k ] = ( va_Polar_Tropopause_75 - va_Polar_SL_75 ) * d_i / d_i_max + va_Polar_SL_75;
			}
		}
	}

	for ( int j = j_pol_v_s; j < j_pol_v_s + 1; j++ )
	{
		i_max = im_tropopause[ j ];
		d_i_max = ( double ) ( i_max - im_1 );
		if ( d_i_max == 0. ) d_i_max = 1.e-6;

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = i_max; i < im; i++ )
			{
				d_i = ( double ) ( i - im_1 );
				v.x[ i ][ j ][ k ] = va_Polar_Tropopause_75 * d_i / d_i_max;
			}
		}
	}


///////////////////////////////////////// end southern polar cell ///////////////////////////////////




/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////
/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////



////////////// meridional values of w-velocity component from Pol till Ferrel, from Ferrel till Hadley, from Hadley till equator ///////////////


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// w-component formed by the diagonal starting from subtropical jet and North Pole
	d_j_60n = ( double ) j_fer_n;
	d_j_90n = ( double ) j_pol_n;
	d_diff = d_j_60n - d_j_90n;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_n; j < j_fer_n + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_fer_n ][ k ] - u.x[ i ][ j_pol_n ][ k ] ) * ( d_j - d_j_90n ) / d_diff + u.x[ i ][ j_pol_n ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_fer_n ][ k ] - w.x[ i ][ j_pol_n ][ k ] ) * ( d_j - d_j_90n ) / d_diff + w.x[ i ][ j_pol_n ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// v-component formed by the diagonal starting from equator till 45°N
	d_j_90n = ( double ) j_pol_n;
	d_j_75n = ( double ) j_pol_v_n;
	d_diff = d_j_75n - d_j_90n;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_n; j < j_pol_v_n + 1; j++ )
		{
			d_j = ( double ) j;
			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_pol_v_n ][ k ] - v.x[ i ][ j_pol_n ][ k ] ) * ( d_j - d_j_90n ) / d_diff + v.x[ i ][ j_pol_n ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// v-component formed by the diagonal starting from 45°N till subtropical jet
	d_j_75n = ( double ) j_pol_v_n;
	d_j_60n = ( double ) j_fer_n;
	d_diff = d_j_60n - d_j_75n;
  
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_v_n; j < j_fer_n + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_fer_n ][ k ] - v.x[ i ][ j_pol_v_n ][ k ] ) * ( d_j - d_j_75n ) / d_diff + v.x[ i ][ j_pol_v_n ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° northern latitude )
// w-component formed by the diagonal starting from tropical jet and subtropical jet
	d_j_60n = ( double ) j_fer_n;
	d_j_30n = ( double ) j_had_n;
	d_diff = d_j_30n - d_j_60n;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_n + 1; j < j_had_n + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_had_n ][ k ] - u.x[ i ][ j_fer_n ][ k ] ) * ( d_j - d_j_60n ) / d_diff + u.x[ i ][ j_fer_n ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_had_n ][ k ] - w.x[ i ][ j_fer_n ][ k ] ) * ( d_j - d_j_60n ) / d_diff + w.x[ i ][ j_fer_n ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° northern latitude )
// v-component formed by the diagonal starting from equator till 45°N
	d_j_60n = ( double ) j_fer_n;
	d_j_45n = ( double ) j_fer_v_n;
	d_diff = d_j_45n - d_j_60n;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_n; j < j_fer_v_n + 1; j++ )
		{
			d_j = ( double ) j;
			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_fer_v_n ][ k ] - v.x[ i ][ j_fer_n ][ k ] ) * ( d_j - d_j_60n ) / d_diff + v.x[ i ][ j_fer_n ][ k ];
			}
		}
	}

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° northern latitude )
// v-component formed by the diagonal starting from 45°N till subtropical jet
	d_j_45n = ( double ) j_fer_v_n;
	d_j_30n = ( double ) j_had_n;
	d_diff = d_j_30n - d_j_45n;
  
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_v_n; j < j_had_n + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_had_n ][ k ] - v.x[ i ][ j_fer_v_n ][ k ] ) * ( d_j - d_j_45n ) / d_diff + v.x[ i ][ j_fer_v_n ][ k ];
			}
		}
	}



/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// v-component formed by the diagonal starting from equator till 45°S
	d_j_30s = ( double ) j_had_s;
	d_j_45s = ( double ) j_fer_v_s;
	d_diff = d_j_45s - d_j_30s;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_v_s; j > j_had_s  - 1; j-- )
		{
			d_j = ( double ) j;
			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_fer_v_s ][ k ] - v.x[ i ][ j_had_s ][ k ] ) * ( d_j - d_j_30s ) / d_diff + v.x[ i ][ j_had_s ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// u- and w-component formed by the diagonal starting from equator and tropical jet
	d_j_5n = ( double ) j_aeq;
	d_j_30n = ( double ) j_had_n;
	d_diff = d_j_5n - d_j_30n;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_n; j < j_aeq + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_aeq ][ k ] - u.x[ i ][ j_had_n ][ k ] ) * ( d_j - d_j_30n ) / d_diff + u.x[ i ][ j_had_n ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_aeq ][ k ] - w.x[ i ][ j_had_n ][ k ] ) * ( d_j - d_j_30n ) / d_diff + w.x[ i ][ j_had_n ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// v-component formed by the diagonal starting from equator till 15°N
	d_j_5n = ( double ) j_aeq;
	d_j_15n = ( double ) j_had_v_n;
	d_diff = d_j_5n - d_j_15n;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_v_n; j < j_aeq + 1; j++ )
		{
			d_j = ( double ) j;
			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_aeq ][ k ] - v.x[ i ][ j_had_v_n ][ k ] ) * ( d_j - d_j_15n ) / d_diff + v.x[ i ][ j_had_v_n ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// v-component formed by the diagonal starting from 15°N till subtropical jet
	d_j_15n = ( double ) j_had_v_n;
	d_j_30n = ( double ) j_had_n;
	d_diff = d_j_15n - d_j_30n;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_n; j < j_had_v_n + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_had_v_n ][ k ] - v.x[ i ][ j_had_n ][ k ] ) * ( d_j - d_j_30n ) / d_diff + v.x[ i ][ j_had_n ][ k ];
			}
		}
	}



/////////////////////////////////////////// change in j-direction in southern hemisphere //////////////////////////////////////////////




/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// w-component formed by the diagonal starting from tropical jet and subtropical jet
	d_j_60s = ( double ) j_fer_s;
	d_j_90s = ( double ) j_pol_s;
	d_diff = d_j_60s - d_j_90s;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_s; j < j_pol_s + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_fer_s ][ k ] - u.x[ i ][ j_pol_s ][ k ] ) * ( d_j - d_j_90s ) / d_diff + u.x[ i ][ j_pol_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_fer_s ][ k ] - w.x[ i ][ j_pol_s ][ k ] ) * ( d_j - d_j_90s ) / d_diff + w.x[ i ][ j_pol_s ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// v-component formed by the diagonal starting from 60°S till 75°S
	d_j_75s = ( double ) j_pol_v_s;
	d_j_90s = ( double ) j_pol_s;
	d_diff = d_j_90s - d_j_75s;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_s; j > j_pol_v_s - 1; j-- )
		{
			d_j = ( double ) j;
			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_pol_s ][ k ] - v.x[ i ][ j_pol_v_s ][ k ] ) * ( d_j - d_j_75s ) / d_diff + v.x[ i ][ j_pol_v_s ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// v-component formed by the diagonal starting from 60°S till 75°S
	d_j_60s = ( double ) j_fer_s;
	d_j_75s = ( double ) j_pol_v_s;
	d_diff = d_j_75s - d_j_60s;
  
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_pol_v_s; j > j_fer_s - 1; j-- )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_pol_v_s ][ k ] - v.x[ i ][ j_fer_s ][ k ] ) * ( d_j - d_j_60s ) / d_diff + v.x[ i ][ j_fer_s ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// w-component formed by the diagonal starting from polar jet and subtropical jet
	d_j_60s = ( double ) j_fer_s;
	d_j_30s = ( double ) j_had_s;
	d_diff = d_j_30s - d_j_60s;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_s; j < j_fer_s + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_had_s ][ k ] - u.x[ i ][ j_fer_s ][ k ] ) * ( d_j - d_j_60s ) / d_diff + u.x[ i ][ j_fer_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_had_s ][ k ] - w.x[ i ][ j_fer_s ][ k ] ) * ( d_j - d_j_60s ) / d_diff + w.x[ i ][ j_fer_s ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// v-component formed by the diagonal starting from 45°S till subtropical jet
	d_j_45s = ( double ) j_fer_v_s;
	d_j_60s = ( double ) j_fer_s;
	d_diff = d_j_60s - d_j_45s;
  
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_fer_s; j > j_fer_v_s - 1; j-- )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_fer_s ][ k ] - v.x[ i ][ j_fer_v_s ][ k ] ) * ( d_j - d_j_45s ) / d_diff + v.x[ i ][ j_fer_v_s ][ k ];
			}
		}
	}


///////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Hadley cell
// u- and w-component formed by the diagonal starting from equatorjet and subtropical jet
	d_j_5s = ( double ) j_aeq;
	d_j_30s = ( double ) j_had_s;
	d_diff = d_j_5s - d_j_30s;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_aeq; j < j_had_s + 1; j++ )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_aeq ][ k ] - u.x[ i ][ j_had_s ][ k ] ) * ( d_j - d_j_30s ) / d_diff + u.x[ i ][ j_had_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_aeq ][ k ] - w.x[ i ][ j_had_s ][ k ] ) * ( d_j - d_j_30s ) / d_diff + w.x[ i ][ j_had_s ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Hadley cell
// v-component formed by the diagonal starting from 15°S till subtropical jet
	d_j_15s = ( double ) j_had_v_s;
	d_j_30s = ( double ) j_had_s;
	d_diff = d_j_30s - d_j_15s;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_s; j > j_had_v_s - 1; j-- )
		{
			d_j = ( double ) j;

			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_had_s ][ k ] - v.x[ i ][ j_had_v_s ][ k ] ) * ( d_j - d_j_15s ) / d_diff + v.x[ i ][ j_had_v_s ][ k ];
			}
		}
	}


/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Hadley cell
// v-component formed by the diagonal starting from equator till 15°N
	d_j_5s = ( double ) j_aeq;
	d_j_15s = ( double ) j_had_v_s;
	d_diff = d_j_15s - d_j_5s;
 
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = j_had_v_s; j > j_aeq - 1; j-- )
		{
			d_j = ( double ) j;
			for ( int i = 0; i < im; i++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_had_v_s ][ k ] - v.x[ i ][ j_aeq ][ k ] ) * ( d_j - d_j_5s ) / d_diff + v.x[ i ][ j_aeq ][ k ];
			}
		}
	}


///////////////////////////////////////////////// change in sign of v-component /////////////////////////////////////////////////
///////////////////////////////////////////////// values identical with the northern hemisphere //////////////////////////////////////////////////
	for ( int i = 0; i < im; i++ )
	{
		for ( int j = j_aeq + 1; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				v.x[ i ][ j ][ k ] = - v.x[ i ][ j ][ k ];
			}
		}
	}


/////////////////////////////////////////////// end forming diagonals ///////////////////////////////////////////////////







/////////////////////////////////////////// forming diagonals to simulate thermic highs ///////////////////////////////////////////////////////
/////////////////////////////////////////// forming diagonals to simulate thermic highs ///////////////////////////////////////////////////////



///////////////////////////////////////////////// Pacific ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// North Pacific
///////////////////////////////////////////////// change in sign of v-component in the northern hemisphere //////////////////////////
//	ua_00 = .03;																		// in m/s compares to 1.08 km/h, non-dimensionalized by u_0 below
//	d_i_half = ( double ) i_half;

	j_had_n_end = j_had_n - 8;

	k_w = 120;
	k_w_end = k_w - 5;
	k_e = 240;

	d_j_w = ( double ) j_aeq;

	for ( int i = 0; i <= i_half; i++ )
	{
//		d_i = ( double ) i;

		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			v.x[ i ][ j ][ k_w ] = - .5 * v.x[ i ][ j ][ k_w ];
//			v.x[ i ][ j ][ k_w ] = - v.x[ i ][ j ][ k_w ];
//			d_j = ( double ) j;																// asymmetric extention of highs
//			v.x[ i ][ j ][ k_w ] = .5 * ( v.x[ i ][ j_had_n ][ k_w ] - v.x[ i ][ j_had_s ][ k_w ] ) * ( d_j * d_j / ( d_j_w * d_j_w ) - 2. * d_j / d_j_w );
		}
	}

/////////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_w; k <= k_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_e ] - v.x[ i ][ j ][ k_w ] ) / ( ( double ) ( k_e - k_w ) ) * ( double ) ( k - k_w ) + v.x[ i ][ j ][ k_w ];
//				u.x[ i ][ j ][ k ] = - ua_00 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}

///////////////////////////////////////////// smoothing transitions /////////////////////////////////////////////////

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_w_end + 1; k <= k_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_w ] - v.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + v.x[ i ][ j ][ k_w_end ];
//				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_w ] - u.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + u.x[ i ][ j ][ k_w_end ];
			}
		}
	}



///////////////////////////////////////////////// Pacific ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// South Pacific
///////////////////////////////////////////////// change in sign of v-component in the southern hemisphere /////////////////////////////////////////////////


	j_had_s_end = j_had_s + 8;

	k_w = 130;
	k_w_end = k_w - 5;
	k_e = 260;

	for ( int i = 0; i <= i_half; i++ )
	{
//		d_i = ( double ) i;

		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			v.x[ i ][ j ][ k_w ] = - .5 * v.x[ i ][ j ][ k_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			for ( int k = k_w; k <= k_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_e ] - v.x[ i ][ j ][ k_w ] ) / ( ( double ) ( k_e - k_w ) ) * ( double ) ( k - k_w ) + v.x[ i ][ j ][ k_w ];
//				u.x[ i ][ j ][ k ] = - ua_00 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			for ( int k = k_w_end+1; k <= k_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_w ] - v.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + v.x[ i ][ j ][ k_w_end ];
//				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_w ] - u.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + u.x[ i ][ j ][ k_w_end ];
			}
		}
	}





///////////////////////////////////////////////// Indic ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// North Indic .......................................................................... only on land
///////////////////////////////////////////////// change in sign of v-component in the northern hemisphere //////////////////////////

	j_had_n_end = j_had_n - 5;
	j_had_s_end = j_had_s + 5;

	k_w = 30;
	k_w_end = k_w - 10;
	k_e = 90;

	for ( int i = 0; i <= i_half; i++ )
	{
//		d_i = ( double ) i;

		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			v.x[ i ][ j ][ k_w ] = - .5 * v.x[ i ][ j ][ k_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_w; k <= k_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_e ] - v.x[ i ][ j ][ k_w ] ) / ( ( double ) ( k_e - k_w ) ) * ( double ) ( k - k_w ) + v.x[ i ][ j ][ k_w ];
//				u.x[ i ][ j ][ k ] = - ua_00 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_w_end+1; k <= k_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_w ] - v.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + v.x[ i ][ j ][ k_w_end ];
//				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_w ] - u.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + u.x[ i ][ j ][ k_w_end ];
			}
		}
	}



///////////////////////////////////////////////// Indicic ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// South Indic
///////////////////////////////////////////////// change in sign of v-component in the southern hemisphere /////////////////////////////////////////////////

	j_had_n_end = j_had_n - 5;
	j_had_s_end = j_had_s + 5;

	k_w = 35;
	k_w_end = k_w - 3;
	k_e = 90;

	for ( int i = 0; i <= i_half; i++ )
	{
//		d_i = ( double ) i;

		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			v.x[ i ][ j ][ k_w ] = - .5 * v.x[ i ][ j ][ k_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			for ( int k = k_w; k <= k_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_e ] - v.x[ i ][ j ][ k_w ] ) / ( ( double ) ( k_e - k_w ) ) * ( double ) ( k - k_w ) + v.x[ i ][ j ][ k_w ];
//				u.x[ i ][ j ][ k ] = - ua_00 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			for ( int k = k_w_end+1; k <= k_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_w ] - v.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + v.x[ i ][ j ][ k_w_end ];
//				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_w ] - u.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + u.x[ i ][ j ][ k_w_end ];
			}
		}
	}



///////////////////////////////////////////////// Atlantic ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// North Altlantic
///////////////////////////////////////////////// change in sign of v-component in the northern hemisphere //////////////////////////

	j_had_n_end = j_had_n - 5;
	j_had_s_end = j_had_s + 5;

	k_w = 280;
//	k_w_end = k_w - 10;
//	k_e = 340;
	k_w_end = k_w - 5;
	k_e = 330;

	for ( int i = 0; i <= i_half; i++ )
	{
//		d_i = ( double ) i;

		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			v.x[ i ][ j ][ k_w ] = - .5 * v.x[ i ][ j ][ k_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_w; k <= k_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_e ] - v.x[ i ][ j ][ k_w ] ) / ( ( double ) ( k_e - k_w ) ) * ( double ) ( k - k_w ) + v.x[ i ][ j ][ k_w ];
//				u.x[ i ][ j ][ k ] = - ua_00 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_w_end+1; k <= k_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_w ] - v.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + v.x[ i ][ j ][ k_w_end ];
//				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_w ] - u.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + u.x[ i ][ j ][ k_w_end ];
			}
		}
	}



///////////////////////////////////////////////// Atlantic ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// South Atlantic
///////////////////////////////////////////////// change in sign of v-component in the southern hemisphere /////////////////////////////////////////////////

	j_had_n_end = j_had_n - 5;
	j_had_s_end = j_had_s + 5;

	k_w = 320;
//	k_w_end = k_w - 10;
//	k_e = 360;
	k_w_end = k_w - 5;
	k_e = 360;

	for ( int i = 0; i <= i_half; i++ )
	{
//		d_i = ( double ) i;

		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			v.x[ i ][ j ][ k_w ] = - .5 * v.x[ i ][ j ][ k_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			for ( int k = k_w; k <= k_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_e ] - v.x[ i ][ j ][ k_w ] ) / ( ( double ) ( k_e - k_w ) ) * ( double ) ( k - k_w ) + v.x[ i ][ j ][ k_w ];
//				u.x[ i ][ j ][ k ] = - ua_00 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_had_s_end; j++ )
		{
			for ( int k = k_w_end+1; k <= k_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_w ] - v.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + v.x[ i ][ j ][ k_w_end ];
//				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j ][ k_w ] - u.x[ i ][ j ][ k_w_end ] ) / ( ( double ) ( k_w - k_w_end ) ) * ( double ) ( k - k_w_end ) + u.x[ i ][ j ][ k_w_end ];
			}
		}
	}


/////////////////////////////////////////// end of forming diagonals to simulate thermic highs ///////////////////////////////////////////////////////
/////////////////////////////////////////// end of forming diagonals to simulate thermic highs ///////////////////////////////////////////////////////





///////////////////////////////////////////////// smoothing transitions from cell to cell //////////////////////////
///////////////////////////////////////////////// smoothing transitions from cell to cell //////////////////////////


///////////////////////////////////////////////// Northern hemisphere ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Northern hemisphere
///////////////////////////////////////////////// smoothing transitions from Headley to Ferrel to Polar cells //////////////////////////

	j_s = j_had_n - 3;
	j_n = j_had_n + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}


	j_s = j_fer_n - 3;
	j_n = j_fer_n + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}



// Northern hemisphere
///////////////////////////////////////////////// smoothing transitions around jets /////////////////////////////////////////////////////

	j_s = j_had_v_n - 3;
	j_n = j_had_v_n + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}


	j_s = j_fer_v_n - 3;
	j_n = j_fer_v_n + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}



///////////////////////////////////////////////// smoothing transitions around equator //////////////////////////

	j_s = j_aeq - 3;
	j_n = j_aeq + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}



///////////////////////////////////////////////// Southern hemisphere ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Southern hemisphere
///////////////////////////////////////////////// smoothing transitions from Headley to Ferrel to Polar cells //////////////////////////

	j_s = j_had_s - 3;
	j_n = j_had_s + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}


	j_s = j_fer_s - 3;
	j_n = j_fer_s + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}


// Southern hemisphere
///////////////////////////////////////////////// smoothing transitions around jets /////////////////////////////////////////////////////

	j_s = j_had_v_s - 3;
	j_n = j_had_v_s + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}


	j_s = j_fer_v_s - 3;
	j_n = j_fer_v_s + 3;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = j_s; j <= j_n; j++ )
			{
				u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_n ][ k ] - u.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + u.x[ i ][ j_s ][ k ];
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_n ][ k ] - v.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + v.x[ i ][ j_s ][ k ];
				w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_n ][ k ] - w.x[ i ][ j_s ][ k ] ) / ( ( int ) ( j_n - j_s ) ) * ( int ) ( j - j_s ) + w.x[ i ][ j_s ][ k ];
			}
		}
	}



// non dimensionalization by u_0
	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				u.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ] / u_0;
				v.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] / u_0;
				w.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] / u_0;
//				if ( h.x[ i ][ j ][ k ] == 1. )		u.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] = 0.;
			}
		}
	}

///////////////////////////////////////////////// end of smoothing transitions from cell to cell //////////////////////////
///////////////////////////////////////////////// end of smoothing transitions from cell to cell //////////////////////////
}





void BC_Thermo::BC_Surface_Temperature_NASA ( const string &Name_SurfaceTemperature_File, Array_2D &temperature_NASA, Array &t )
{
// initial conditions for the Name_SurfaceTemperature_File at the sea surface

	cout.precision ( 3 );
	cout.setf ( ios::fixed );

	ifstream Name_SurfaceTemperature_File_Read;
	Name_SurfaceTemperature_File_Read.open(Name_SurfaceTemperature_File);

	if (!Name_SurfaceTemperature_File_Read.is_open()) {
		cerr << "ERROR: could not open SurfaceTemperature_File file at " << Name_SurfaceTemperature_File << "\n";
		abort();
	}

	k_half = ( km -1 ) / 2;																// position at 180°E ( Greenwich )

	j = 0;
	k = 0;


	while ( ( k < km ) && !Name_SurfaceTemperature_File_Read.eof() )
	{
		while ( j < jm )
		{
			Name_SurfaceTemperature_File_Read >> dummy_1;
			Name_SurfaceTemperature_File_Read >> dummy_2;
			Name_SurfaceTemperature_File_Read >> dummy_3;

			t.x[ 0 ][ j ][ k ] = temperature_NASA.y[ j ][ k ] = ( dummy_3 + 273.15 ) / 273.15;

			j++;
		}
	j = 0;
	k++;
	}

// correction of surface temperature around 180°E
	for ( int j = 0; j < jm; j++ )
	{
		t.x[ 0 ][ j ][ k_half ] = ( t.x[ 0 ][ j ][ k_half + 1 ] + t.x[ 0 ][ j ][ k_half - 1 ] ) / 2.;
		temperature_NASA.y[ j ][ k_half ] = ( temperature_NASA.y[ j ][ k_half + 1 ] + temperature_NASA.y[ j ][ k_half - 1 ] ) / 2.;
	}

	Name_SurfaceTemperature_File_Read.close();

}





void BC_Thermo::BC_NASAbasedSurfTempRead ( const string &Name_NASAbasedSurfaceTemperature_File, double &t_cretaceous, double &t_cret_cor, Array &t, Array &c, Array &cloud, Array &ice )
{
    //string NASAbasedSurfaceTemperature_File = output_path + "/NASAbasedSurfaceTemperature.xyz";
    ifstream Name_NASAbasedSurfaceTemperature_File_Read;
    Name_NASAbasedSurfaceTemperature_File_Read.precision(4);
    Name_NASAbasedSurfaceTemperature_File_Read.setf(ios::fixed);
    Name_NASAbasedSurfaceTemperature_File_Read.open(Name_NASAbasedSurfaceTemperature_File);

    if (!Name_NASAbasedSurfaceTemperature_File_Read.is_open()) {
        cerr << Name_NASAbasedSurfaceTemperature_File << endl;
        cerr << "ERROR: could not open NASAbased file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }


	Name_NASAbasedSurfaceTemperature_File_Read >> t_cretaceous >> t_cret_cor;

	for ( int k = 0; k < km; k++ ) 
	{
		for ( int j = 0; j < jm; j++ ) 
		{
//			for ( int i = 0; i < im; i++ ) 
			{
/*
				Name_NASAbasedSurfaceTemperature_File_Read >> t.x[ i ][ j ][ k ];
				Name_NASAbasedSurfaceTemperature_File_Read >> c.x[ i ][ j ][ k ];
				Name_NASAbasedSurfaceTemperature_File_Read >> cloud.x[ i ][ j ][ k ];
				Name_NASAbasedSurfaceTemperature_File_Read >> ice.x[ i ][ j ][ k ];
*/
				Name_NASAbasedSurfaceTemperature_File_Read >> t.x[ 0 ][ j ][ k ];
//				Name_NASAbasedSurfaceTemperature_File_Read >> c.x[ 0 ][ j ][ k ];
//				Name_NASAbasedSurfaceTemperature_File_Read >> cloud.x[ 0 ][ j ][ k ];
//				Name_NASAbasedSurfaceTemperature_File_Read >> ice.x[ 0 ][ j ][ k ];
			}
        }
    }
    Name_NASAbasedSurfaceTemperature_File_Read.close();
}




void BC_Thermo::BC_NASAbasedSurfTempWrite ( const string &Name_NASAbasedSurfaceTemperature_File, double &t_cretaceous, double &t_cret_cor, Array &t, Array &c, Array &cloud, Array &ice )
{
// initial conditions for the Name_NASAbasedSurfaceTemperature_File at the sea surface

	cout.precision ( 4 );
	cout.setf ( ios::fixed );

    string NASAbasedSurfaceTemperature_File = output_path + "/NASAbasedSurfaceTemperature.xyz";

	ofstream Name_NASAbasedSurfaceTemperature_File_Write;
	Name_NASAbasedSurfaceTemperature_File_Write.open(NASAbasedSurfaceTemperature_File);

	if (!Name_NASAbasedSurfaceTemperature_File_Write.is_open()) {
		cerr << "ERROR: could not open Name_NASAbasedSurfaceTemperature_File file at " << Name_NASAbasedSurfaceTemperature_File << "\n";
		abort();
	}

	Name_NASAbasedSurfaceTemperature_File_Write << t_cretaceous << " " << t_cret_cor << endl;

	for ( int k = 0; k < km; k++ ) 
	{
		for ( int j = 0; j < jm; j++ ) 
		{
//			for ( int i = 0; i < im; i++ ) 
			{
				Name_NASAbasedSurfaceTemperature_File_Write << t.x[ i ][ j ][ k ] << " " << c.x[ i ][ j ][ k ] << " " << cloud.x[ i ][ j ][ k ] << " " << ice.x[ i ][ j ][ k ] << endl;

//				Name_NASAbasedSurfaceTemperature_File_Write << t.x[ 0 ][ j ][ k ] << " " << c.x[ 0 ][ j ][ k ] << " " << cloud.x[ 0 ][ j ][ k ] << " " << ice.x[ 0 ][ j ][ k ] << endl;
				Name_NASAbasedSurfaceTemperature_File_Write << t.x[ 0 ][ j ][ k ] << endl;

//				c.x[ i ][ j ][ k ] = 0.;
//				cloud.x[ i ][ j ][ k ] = 0.;
//				ice.x[ i ][ j ][ k ] = 0.;
			}
		}
	}
	Name_NASAbasedSurfaceTemperature_File_Write.close();
}





void BC_Thermo::BC_Surface_Precipitation_NASA ( const string &Name_SurfacePrecipitation_File, Array_2D &precipitation_NASA )
{
// initial conditions for the Name_SurfacePrecipitation_File at the sea surface

	cout.precision ( 3 );
	cout.setf ( ios::fixed );

	ifstream Name_SurfacePrecipitation_File_Read;
	Name_SurfacePrecipitation_File_Read.open(Name_SurfacePrecipitation_File);

	if (!Name_SurfacePrecipitation_File_Read.is_open()) {
		cerr << "ERROR: could not open SurfacePrecipitation_File file at " << Name_SurfacePrecipitation_File << "\n";
		abort();
	}

	j = 0;
	k = 0;


	while ( ( k < km ) && !Name_SurfacePrecipitation_File_Read.eof() )
	{
		while ( j < jm )
		{
			Name_SurfacePrecipitation_File_Read >> dummy_1;
			Name_SurfacePrecipitation_File_Read >> dummy_2;
			Name_SurfacePrecipitation_File_Read >> dummy_3;

			precipitation_NASA.y[ j ][ k ] = dummy_3;
			j++;
		}
	j = 0;
	k++;
	}

	Name_SurfacePrecipitation_File_Read.close();

}




void BC_Thermo::BC_Pressure ( Array &p_stat, Array &p_dyn, Array &t, Array &h )
{
	exp_pressure = g / ( 1.e-2 * gam * R_Air );

// boundary condition of surface pressure given by surface temperature through gas equation
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			p_stat.x[ 0 ][ j ][ k ] =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );		// given in hPa
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 1; i < im; i++ )
			{
//				p_stat.x[ i ][ j ][ k ] = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t.x[ i ][ j ][ k ] * t_0 ) ) * p_stat.x[ 0 ][ j ][ k ];	// isothermal atmosphere T = const

				hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

				p_stat.x[ i ][ j ][ k ] = pow ( ( ( t.x[ 0 ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) / ( t.x[ 0 ][ j ][ k ] * t_0 ) ), exp_pressure ) * p_stat.x[ 0 ][ j ][ k ];	// linear temperature distribution T = T0 - gam * hight
																									// current air pressure, step size in 500 m, from politropic formula in hPa
//					if ( ( j == 90 ) && ( k == 180 ) )	cout << " pressure  i = "<< i << " exp_pressure = "<< exp_pressure << " hight = "<< hight << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  p_stat = " << p_stat.x[ i ][ j ][ k ] << endl;
			}
		}
	}


/*
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				p_dyn.x[ i ][ j ][ k ] = p_stat.x[ i ][ j ][ k ] / p_stat.x[ 0 ][ j ][ k ] / 10.;
			}
		}
	}
*/
}









void BC_Thermo::Latent_Heat ( Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &tn, Array &u, Array &v, Array &w, Array &p_dyn, Array &p_stat, Array &c, Array &Latency, Array &Q_Sensible, Array &t_cond_3D, Array &t_evap_3D, Array &radiation_3D )
{
	double Latency_Ice = 0.; 

// collection of coefficients for phase transformation
	coeff_L = 100. * r_air * lv * u_0 / L_atm;							// coefficient for Latency
	coeff_Q = 100. * cp_l * r_air * u_0 * t_0 / L_atm;				// coefficient for Q_Sensible

	c32 = 3. / 2.;
	c42 = 4. / 2.;
	c12 = 1. / 2.;

// 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )
// collection of coefficients
	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;
	rm = rad.z[ 0 ];

	for ( int j = 1; j < jm-1; j++ )
	{
// collection of coefficients
		sinthe = sin( the.z[ j ] );
		rmsinthe = rm * sinthe;

		for ( int k = 1; k < km-1; k++ )
		{
/*
			t_Celsius_0 = t.x[ 0 ][ j ][ k ] * t_0 - t_0;																// conversion from Kelvin to Celsius
			t_Celsius_nj = t.x[ 0 ][ j - 1 ][ k ] * t_0 - t_0;														// conversion from Kelvin to Celsius
			t_Celsius_pj = t.x[ 0 ][ j + 1 ][ k ] * t_0 - t_0;														// conversion from Kelvin to Celsius
			t_Celsius_nk = t.x[ 0 ][ j ][ k - 1 ] * t_0 - t_0;														// conversion from Kelvin to Celsius
			t_Celsius_pk= t.x[ 0 ][ j ][ k + 1 ] * t_0 - t_0;														// conversion from Kelvin to Celsius

			E_Rain = hp * exp_func ( t_Celsius_0, 17.0809, 234.175 );		// saturation water vapour pressure for the water phase at t > 0°C in hPa
			E_Ice = hp * exp_func ( t_Celsius_0, 22.4429, 272.44 );			// saturation water vapour pressure for the ice phase in hPa

			p_h = p_stat.x[ 0 ][ j ][ k ];

			q_Rain  = ep * E_Rain / ( p_h - E_Rain );															// water vapour amount at saturation with water formation in kg/kg
			q_Ice  = ep * E_Ice / ( p_h - E_Ice );																	// water vapour amount at saturation with ice formation in kg/kg

			dp_hdthe = ( p_stat.x[ 0 ][ j+1 ][ k ] - p_stat.x[ 0 ][ j-1 ][ k ] ) / ( 2. * rm * dthe );
			dp_hdphi = ( p_stat.x[ 0 ][ j ][ k+1 ] - p_stat.x[ 0 ][ j ][ k-1 ] ) / ( 2. * rmsinthe * dphi );

			E_dEdthe_Rain = hp * ( exp_func ( t_Celsius_pj, 17.0809, 234.175 ) - exp_func ( t_Celsius_nj, 17.0809, 234.175 ) ) / ( 2. * rm * dthe );
			E_dEdthe_Ice = hp * ( exp_func ( t_Celsius_pj, 22.4429, 272.44 ) - exp_func ( t_Celsius_nj, 22.4429, 272.44 ) ) / ( 2. * rm * dthe );

			E_dEdphi_Rain = hp * ( exp_func ( t_Celsius_pk, 17.0809, 234.175 ) - exp_func ( t_Celsius_nk, 17.0809, 234.175 ) ) / ( 2. * rmsinthe * dphi );
			E_dEdphi_Ice = hp * ( exp_func ( t_Celsius_pk, 22.4429, 272.44 ) - exp_func ( t_Celsius_nk, 22.4429, 272.44 ) ) / ( 2. * rmsinthe * dphi );

			v_av = .5 * ( v.x[ 0 ][ j+1 ][ k ] + v.x[ 0 ][ j-1 ][ k ] );
			w_av = .5 * ( w.x[ 0 ][ j ][ k+1 ] + w.x[ 0 ][ j ][ k-1 ] );
*/
//			Latency.x[ 0 ][ j ][ k ] = coeff_L * q_Rain * ( v_av * ( E_dEdthe_Rain / E_Rain - dp_hdthe / p_h ) + w_av * ( E_dEdphi_Rain / E_Rain - dp_hdphi / p_h ) );
			Latency.x[ 0 ][ j ][ k ] = 0.;

//			Latency.x[ 0 ][ j ][ k ] = Latency.x[ 0 ][ j ][ k ] + coeff_L * q_Ice * ( v_av * ( E_dEdthe_Ice / E_Ice - dp_hdthe / p_h ) + w_av * ( E_dEdphi_Ice / E_Ice - dp_hdphi / p_h ) );

//			Q_Sensible.x[ 0 ][ j ][ k ] = - coeff_Q * ( v_av * ( t.x[ 0 ][ j+1 ][ k ] - t.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe * rm ) + w_av * ( t.x[ 0 ][ j ][ k+1 ] - t.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe ) );							// sensible heat in [W/m2] from energy transport equation
			Q_Sensible.x[ 0 ][ j ][ k ] = 0.;							// sensible heat in [W/m2] from energy transport equation


			if ( Latency.x[ 0 ][ j ][ k ] <= 0. )			t_cond_3D.x[ 0 ][ j ][ k ] = pow ( ( radiation_3D.x[ 0 ][ j ][ k ] - Latency.x[ 0 ][ j ][ k ] - Q_Sensible.x[ 0 ][ j ][ k ] ) / sigma, .25 ) - pow ( radiation_3D.x[ 0 ][ j ][ k ] / sigma, .25 );						// temperature increase due to condensation
			else 													t_cond_3D.x[ 0 ][ j ][ k ] = 0.;

			if ( Latency.x[ 0 ][ j ][ k ] > 0. )			t_evap_3D.x[ 0 ][ j ][ k ] = pow ( ( radiation_3D.x[ 0 ][ j ][ k ] - Latency.x[ 0 ][ j ][ k ] - Q_Sensible.x[ 0 ][ j ][ k ] ) / sigma, .25 ) - pow ( radiation_3D.x[ 0 ][ j ][ k ] / sigma, .25 );						// temperature decrease due to evaporation
			else 													t_evap_3D.x[ 0 ][ j ][ k ] = 0.;


/*
			cout.precision ( 4 );
			if ( ( j == 90 ) && ( k == 180 ) )	cout << 0 << "   " << j << "   " << k << "   " << Latency.x[ 0 ][ j ][ k ] << "   " << E_dEdr_Rain << "   " << E_dEdthe_Rain << "   " << E_dEdphi_Rain << "   " << dp_hdr << "   " << dp_hdthe << "   " << dp_hdphi << "   " << p_h << "   " << u.x[ 0 ][ j ][ k ] << "   " << v.x[ 0 ][ j ][ k ] << "   " << w.x[ 0 ][ j ][ k ] << "   " << ( E_dEdr_Rain / E_Rain - dp_hdr / p_h ) << "   " << ( E_dEdthe_Rain / E_Rain - dp_hdthe / p_h ) << "   " << ( E_dEdphi_Rain / E_Rain - dp_hdphi / p_h ) << "   " << coeff_lv * u.x[ 0 ][ j ][ k ] * q_Rain * ( E_dEdr_Rain / E_Rain - dp_hdr / p_h ) << "   " << coeff_lv * v.x[ 0 ][ j ][ k ] * q_Rain * ( E_dEdthe_Rain / E_Rain - dp_hdthe / p_h ) << "   " << coeff_lv * w.x[ 0 ][ j ][ k ] * q_Rain * ( E_dEdphi_Rain / E_Rain - dp_hdphi / p_h ) << "   " << t.x[ 0 ][ j ][ k ] << endl;
			cout.precision ( 8 );
*/
		}
	}


	for ( int j = 1; j < jm-1; j++ )
	{
// collection of coefficients
		sinthe = sin( the.z[ j ] );
		sinthe2 = sinthe * sinthe;
		costhe = cos( the.z[ j ] );
		cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
		rmsinthe = rm * sinthe;
		rm2sinthe = rm2 * sinthe;
		rm2sinthe2 = rm2 * sinthe2;

// water vapour can condensate/evaporate und sublimate/vaporize
// water vapour turns to or developes from water or ice
// latent heat of water vapour

		for ( int k = 1; k < km-1; k++ )
		{
			int i_mount = 0;

			for ( int i = im-2; i >= 1; i-- )
			{
// collection of coefficients
				rm = rad.z[ i ];
				rm2 = rm * rm;

				t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;																		// conversion from Kelvin to Celsius
				t_Celsius_ni = t.x[ i - 1 ][ j ][ k ] * t_0 - t_0;																// conversion from Kelvin to Celsius
				t_Celsius_pi = t.x[ i + 1 ][ j ][ k ] * t_0 - t_0;															// conversion from Kelvin to Celsius
				T = t.x[ i ][ j ][ k ] * t_0;
//				t_Celsius_nj = t.x[ i ][ j - 1 ][ k ] * t_0 - t_0;																// conversion from Kelvin to Celsius
//				t_Celsius_pj = t.x[ i ][ j + 1 ][ k ] * t_0 - t_0;															// conversion from Kelvin to Celsius
//				t_Celsius_nk = t.x[ i ][ j ][ k - 1 ] * t_0 - t_0;															// conversion from Kelvin to Celsius
//				t_Celsius_pk= t.x[ i ][ j ][ k + 1 ] * t_0 - t_0;															// conversion from Kelvin to Celsius

				p_h = p_stat.x[ i ][ j ][ k ];
/*
				e_h = c.x[ i ][ j ][ k ] * p_h / ep; 																				// water vapour pressure in hPa
				a_h = 216.6 * e_h / ( t.x[ i ][ j ][ k ] * t_0 );																// absolute humidity in kg/m3
				q_h = c.x[ i ][ j ][ k ];																								// threshold value for water vapour at local hight h in kg/kg
*/
				E_Rain = hp * exp_func ( T, 17.2694, 35.86 );										// saturation water vapour pressure for the water phase at t > 0°C in hPa
				E_Ice = hp * exp_func ( T, 21.8746, 7.66 );											// saturation water vapour pressure for the ice phase in hPa

				q_Rain  = ep * E_Rain / ( p_h - E_Rain );																// water vapour amount at saturation with water formation in kg/kg
				q_Ice  = ep * E_Ice / ( p_h - E_Ice );																		// water vapour amount at saturation with ice formation in kg/kg

//				dp_hdr = ( p_stat.x[ i+1 ][ j ][ k ] - p_stat.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
//				dp_hdthe = ( p_stat.x[ i ][ j+1 ][ k ] - p_stat.x[ i ][ j-1 ][ k ] ) / ( 2. * rm * dthe );
//				dp_hdphi = ( p_stat.x[ i ][ j ][ k+1 ] - p_stat.x[ i ][ j ][ k-1 ] ) / ( 2. * rmsinthe * dphi );

//				E_dEdr_Rain = hp * ( exp_func ( t_Celsius_pi, 17.0809, 234.175 ) - exp_func ( t_Celsius_ni, 17.0809, 234.175 ) ) / ( 2. * dr );
//				E_dEdr_Ice = hp * ( exp_func ( t_Celsius_pi, 22.4429, 272.44 ) - exp_func ( t_Celsius_ni, 22.4429, 272.44 ) ) / ( 2. * dr );

//				E_dEdthe_Rain = hp * ( exp_func ( t_Celsius_pj, 17.0809, 234.175 ) - exp_func ( t_Celsius_nj, 17.0809, 234.175 ) ) / ( 2. * rm * dthe );
//				E_dEdthe_Ice = hp * ( exp_func ( t_Celsius_pj, 22.4429, 272.44 ) - exp_func ( t_Celsius_nj, 22.4429, 272.44 ) ) / ( 2. * rm * dthe );

//				E_dEdphi_Rain = hp * ( exp_func ( t_Celsius_pk, 17.0809, 234.175 ) - exp_func ( t_Celsius_nk, 17.0809, 234.175 ) ) / ( 2. * rmsinthe * dphi );
//				E_dEdphi_Ice = hp * ( exp_func ( t_Celsius_pk, 22.4429, 272.44 ) - exp_func ( t_Celsius_nk, 22.4429, 272.44 ) ) / ( 2. * rmsinthe * dphi );

				u_av = .5 * ( u.x[ i+1 ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] );
//				v_av = .5 * ( v.x[ i ][ j+1 ][ k ] + v.x[ i ][ j-1 ][ k ] );
//				w_av = .5 * ( w.x[ i ][ j ][ k+1 ] + w.x[ i ][ j ][ k-1 ] );

//				Latency.x[ i ][ j ][ k ] = coeff_L * ( u_av * q_Rain * ( E_dEdr_Rain / E_Rain - dp_hdr / p_h ) + v_av * q_Rain * ( E_dEdthe_Rain / E_Rain - dp_hdthe / p_h ) + w_av * q_Rain * ( E_dEdphi_Rain / E_Rain - dp_hdphi / p_h ) );
//				if ( c.x[ i ][ j ][ k ] >= q_Rain )			Latency.x[ i ][ j ][ k ] = coeff_L * ( u_av * q_Rain * ( E_dEdr_Rain / E_Rain - dp_hdr / p_h ) );
				if ( c.x[ i ][ j ][ k ] >= q_Rain )			Latency.x[ i ][ j ][ k ] = coeff_L * ( u_av * ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) );
				else 												Latency.x[ i ][ j ][ k ] = 0.;

//				Latency.x[ i ][ j ][ k ] = Latency.x[ i ][ j ][ k ] + coeff_L * ( u_av * q_Ice * ( E_dEdr_Ice / E_Ice - dp_hdr / p_h ) + v_av * q_Ice * ( E_dEdthe_Ice / E_Ice - dp_hdthe / p_h ) + w_av * q_Ice * ( E_dEdphi_Ice / E_Ice - dp_hdphi / p_h ) );
//				if ( c.x[ i ][ j ][ k ] >= q_Ice )				Latency_Ice = coeff_L * ( u_av * q_Ice * ( E_dEdr_Ice / E_Ice - dp_hdr / p_h ) );
				if ( c.x[ i ][ j ][ k ] >= q_Ice )			Latency_Ice = coeff_L * ( u_av * ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) );
				else 												Latency_Ice = 0.;
				Latency.x[ i ][ j ][ k ] = Latency.x[ i ][ j ][ k ] + Latency_Ice;

//				Q_Sensible.x[ i ][ j ][ k ] = - coeff_Q * ( u_av * ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) + v_av * ( t.x[ i ][ j+1 ][ k ] - t.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe * rm ) + w_av * ( t.x[ i ][ j ][ k+1 ] - t.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi * rmsinthe ) );	// sensible heat in [W/m2] from energy transport equation
				Q_Sensible.x[ i ][ j ][ k ] = - coeff_Q * ( u_av * ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) );	// sensible heat in [W/m2] from energy transport equation

				if ( Latency.x[ i ][ j ][ k ] <= 0. )		t_cond_3D.x[ i ][ j ][ k ] = pow ( ( radiation_3D.x[ i ][ j ][ k ] - Latency.x[ i ][ j ][ k ] - Q_Sensible.x[ i ][ j ][ k ] ) / sigma, .25 ) - pow ( radiation_3D.x[ i ][ j ][ k ] / sigma , .25 );																								// temperature increase due to condensation
				else 												t_cond_3D.x[ i ][ j ][ k ] = 0.;

				if ( Latency.x[ i ][ j ][ k ] > 0. )			t_evap_3D.x[ i ][ j ][ k ] = pow ( ( radiation_3D.x[ i ][ j ][ k ] - Latency.x[ i ][ j ][ k ] - Q_Sensible.x[ i ][ j ][ k ] ) / sigma, .25 ) - pow ( radiation_3D.x[ i ][ j ][ k ] / sigma, .25 );																									// temperature decrease due to evaporation
				else 												t_evap_3D.x[ i ][ j ][ k ] = 0.;


/*
				if ( h.x[ i ][ j ][ k ] == 1. )
				{
					if ( i_mount == 0 ) 		i_mount = i;

					Latency.x[ i ][ j ][ k ] = Latency.x[ i_mount + 1 ][ j ][ k ];							// water vapour amount above mount surface repeated
					Q_Sensible.x[ i ][ j ][ k ] = Q_Sensible.x[ i_mount + 1 ][ j ][ k ];				// water vapour amount above mount surface repeated
					t_cond_3D.x[ i ][ j ][ k ] = t_cond_3D.x[ i_mount + 1 ][ j ][ k ];				// water vapour amount above mount surface repeated
					t_evap_3D.x[ i ][ j ][ k ] = t_evap_3D.x[ i_mount + 1 ][ j ][ k ];				// water vapour amount above mount surface repeated
				}
*/

				if ( abs(h.x[ i ][ j ][ k ] - 1.) < std::numeric_limits<double>::epsilon() )
				{
					Latency.x[ i ][ j ][ k ] = 0.;
					Q_Sensible.x[ i ][ j ][ k ] = 0.;
					t_cond_3D.x[ i ][ j ][ k ] = 0.;
					t_evap_3D.x[ i ][ j ][ k ] = 0.;
				}

			}
			i_mount = 0;
		}
	}
}







void BC_Thermo::Ice_Water_Saturation_Adjustment ( int *im_tropopause, int n, int velocity_iter_max, int RadiationModel, Array &h, Array &c, Array &cn, Array &cloud, Array &cloudn, Array &ice, Array &icen, Array &t, Array &p_stat, Array &S_c_c )
{
	cout.precision ( 6 );
// Ice_Water_Saturation_Adjustment, distribution of cloud ice and cloud water dependent on water vapour amount and temperature

// constant coefficients for the adjustment of cloud water and cloud ice amount vice versa

	t_1 = 253.15;
	t_00 = 236.15;
	t_Celsius_1 = t_1 - t_0;																						// -20 °C
	t_Celsius_2 = t_00 - t_0;																						// -37 °C

// setting water vapour, cloud water and cloud ice into the proper thermodynamic ratio based on the local temperatures
// starting from a guessed parabolic temperature and water vapour distribution in north/south direction
//	for ( int k = 1; k < km-1; k++ )
//	{
//		for ( int j = 1; j < jm-1; j++ )
//		{
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			t_u = t.x[ 0 ][ j ][ k ] * t_0;

			t_Celsius_SL = t_u - t_0;
			p_SL = .01 * r_air * R_Air * t_u;																	// from gas equation given in hPa

			E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 );										// saturation water vapour pressure for the water phase at t > 0°C in hPa
			E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 );											// saturation water vapour pressure for the ice phase in hPa

			q_Rain  = ep * E_Rain / ( p_SL - E_Rain );													// water vapour amount at saturation with water formation in kg/kg
			q_Ice = ep * E_Ice / ( p_SL - E_Ice );															// water vapour amount at saturation with ice formation in kg/kg

			r_dry = 100. * p_SL / ( R_Air * t_u );
			r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ 0 ][ j ][ k ] );	// density of humid air, COSMO version withot cloud and ice water, masses negligible

			for ( int i = 1; i < im-1; i++ )
			{
				t_u = t.x[ i ][ j ][ k ] * t_0;																		// in K

				t_Celsius = t_u - t_0;																				// in C
	
				p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );		// given in hPa

				hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

//				if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_u ) ) * p_SL;
				if ( i != 0 ) 			p_h = pow ( ( ( t.x[ 0 ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) / ( t.x[ 0 ][ j ][ k ] * t_0 ) ), exp_pressure ) * p_SL;
				else 						p_h = p_SL;

				r_dry = 100. * p_h / ( R_Air * t_u );
				r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );	// density of humid air, COSMO version withot cloud and ice water, masses negligible

				e_h = .01 * r_humid * R_WaterVapour * t_u;											// delivers the same results

				a_h = 216.6 * e_h / t_u;																			// absolute humidity in kg/m3
				q_h = c.x[ i ][ j ][ k ];																				// threshold value for water vapour at local hight h in kg/kg

				E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 );										// saturation water vapour pressure for the water phase at t > 0°C in hPa
				E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 );											// saturation water vapour pressure for the ice phase in hPa

				q_Rain = ep * E_Rain / ( p_h - E_Rain );													// water vapour amount at saturation with water formation in kg/kg
				q_Ice = ep * E_Ice / ( p_h - E_Ice );														// water vapour amount at saturation with ice formation in kg/kg


// %%%%%%%%%%%%%%%%%%%%%%%%%%%     warm cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%     warm cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



// warm cloud phase in case water vapour is over-saturated
				if ( t_Celsius >= 0. )
				{
					q_T = c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ];										// total water content
					t_u = t.x[ i ][ j ][ k ] * t_0;																	// in K
					t_Celsius = t_u - t_0;

					p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );							// given in hPa

					hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

					if ( i != 0 ) 			p_h = pow ( ( ( t.x[ 0 ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) / ( t.x[ 0 ][ j ][ k ] * t_0 ) ), exp_pressure ) * p_SL;
					else 					p_h = p_SL;

					E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 );										// saturation water vapour pressure for the water phase at t > 0°C in hPa
					q_Rain = ep * E_Rain / ( p_h - E_Rain );												// water vapour amount at saturation with water formation in kg/kg
					q_Rain_n = q_Rain;

					if ( q_T <= q_Rain )																			// subsaturated
					{
						c.x[ i ][ j ][ k ] = q_T;																		// total water amount as water vapour
						cloud.x[ i ][ j ][ k ] = 0.;																	// no cloud water available
						ice.x[ i ][ j ][ k ] = 0.;																		// no cloud ice available above 0 °C

//	if ( ( j == 90 ) && ( k == 180 ) )	cout << "warm cloud q_T   " << " i = "<< i << "  p_h = " << p_h << "  r_humid = " << r_humid << "  t_Celsius = " << t_Celsius << "  t = " << t_u - t_0 << endl << "  q_T = " << q_T * 1000. << "  c = " << c.x[ i ][ j ][ k ] * 1000. << "  cl = " << cloud.x[ i ][ j ][ k ] * 1000. << "  q_Rain= " << q_Rain * 1000. << "  q_Rain_n = " << q_Rain_n * 1000. << endl << endl;

					}
					else 																									// oversaturated
					{
						iter_prec = 0;
						while ( iter_prec <= 20 )																// iter_prec may be varied
						{
							iter_prec = iter_prec + 1;

							T_it = ( t_u + lv / cp_l * c.x[ i ][ j ][ k ] - lv / cp_l * q_Rain );
							E_Rain = hp * exp_func ( T_it, 17.2694, 35.86 );										// saturation water vapour pressure for the water phase at t > 0°C in hPa
							q_Rain = ep * E_Rain / ( p_h - E_Rain );										// water vapour amount at saturation with water formation in kg/kg
							q_Rain = .5 * ( q_Rain_n + q_Rain );

							c.x[ i ][ j ][ k ] = q_Rain;																// water vapour restricted to saturated water vapour amount
							cloud.x[ i ][ j ][ k ] = q_T - c.x[ i ][ j ][ k ];									// cloud water amount
							ice.x[ i ][ j ][ k ] = 0.;																	// no cloud ice available

							q_T = c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ];

//	if ( ( j == 90 ) && ( k == 180 ) )	cout << "warm cloud iter   " << " i = "<< i << " iter_prec   "<< iter_prec << "  p_h = " << p_h << "  r_humid = " << r_humid << "  T_it = " << T_it - t_0 << "  t = " << t_u - t_0 << endl << "  q_T = " << q_T * 1000. << "  c = " << c.x[ i ][ j ][ k ] * 1000. << "  cl = " << cloud.x[ i ][ j ][ k ] * 1000. << "  q_Rain = " << q_Rain * 1000. << "  q_Rain_n = " << q_Rain_n * 1000. << endl << endl;

							if ( c.x[ i ][ j ][ k ] < 0. )									c.x[ i ][ j ][ k ] = 0.;
							if ( cloud.x[ i ][ j ][ k ] < 0. )								cloud.x[ i ][ j ][ k ] = 0.;

							if ( fabs ( q_Rain / q_Rain_n - 1. ) <= 1.e-5 ) 	break;
							q_Rain_n = q_Rain;
						}
					}

					cn.x[ i ][ j ][ k ] = c.x[ i ][ j ][ k ];
					cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ];
					icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ];

					t.x[ i ][ j ][ k ] = T_it / t_0;
				}																												// end ( t_Celsius > 0. )




// %%%%%%%%%%%%%%%%%%%%%%%%%%%     mixed cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%     mixed cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// mixed cloud phase, if 0°C > t > -38°C
//				if ( ( t_Celsius < 0. ) && ( t_Celsius > t_Celsius_2 ) )
				if ( t_Celsius < 0. )
				{
					if ( t_Celsius < t_Celsius_2 ) 										cloud.x[ i ][ j ][ k ] = 0.;
					if ( t_Celsius > 0. ) 													ice.x[ i ][ j ][ k ] = 0.;

					q_v_b = c.x[ i ][ j ][ k ];
					q_c_b = cloud.x[ i ][ j ][ k ];
					q_i_b = ice.x[ i ][ j ][ k ];

					q_T = q_v_b + q_c_b + q_i_b;															// total water content

					t_u = t.x[ i ][ j ][ k ] * t_0;																	// in K
					T = T_nue = t_u;																					// in K

					E_Rain = hp * exp_func ( T, 17.2694, 35.86 );										// saturation water vapour pressure for the water phase at t > 0°C in hPa
					E_Ice = hp * exp_func ( T, 21.8746, 7.66 );											// saturation water vapour pressure for the ice phase in hPa

					q_Rain = ep * E_Rain / ( p_h - E_Rain );											// water vapour amount at saturation with water formation in kg/kg
					q_Ice = ep * E_Ice / ( p_h - E_Ice );												// water vapour amount at saturation with ice formation in kg/kg

					q_v_hyp = q_Rain;

// §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

					iter_prec = 0;
					while ( iter_prec <= 20 )																	// iter_prec may be varied
					{
						iter_prec = iter_prec + 1;

						p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );						// given in hPa

						hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
//						if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * T ) ) * p_SL;						// given in hPa
						if ( i != 0 ) 			p_h = pow ( ( ( T - gam * hight * 1.e-2 ) / ( T ) ), exp_pressure ) * p_SL;						// given in hPa
						else 					p_h = p_SL;

						CND = ( T - t_00 ) / ( t_0 - t_00 );
						DEP = ( t_0 - T ) / ( t_0 - t_00 );

						d_q_v = - ( q_v_b - q_v_hyp );
						d_q_c = - d_q_v * CND;
						d_q_i = - d_q_v * DEP;

						d_t = ( lv * d_q_c + ls * d_q_i ) / cp_l / 1000.;											// in K

						T = T + d_t;																				// in K

						c.x[ i ][ j ][ k ] = q_v_b;
						cloud.x[ i ][ j ][ k ] = q_c_b;
						ice.x[ i ][ j ][ k ] = q_i_b;

						q_v_b = c.x[ i ][ j ][ k ] + d_q_v;
						q_c_b = cloud.x[ i ][ j ][ k ] + d_q_c;
						q_i_b = ice.x[ i ][ j ][ k ] + d_q_i;

						E_Rain = hp * exp_func ( T, 17.2694, 35.86 );										// saturation water vapour pressure for the water phase at t > 0°C in hPa
						E_Ice = hp * exp_func ( T, 21.8746, 7.66 );											// saturation water vapour pressure for the ice phase in hPa

						q_Rain = ep * E_Rain / ( p_h - E_Rain );											// water vapour amount at saturation with water formation in kg/kg
						q_Ice = ep * E_Ice / ( p_h - E_Ice );												// water vapour amount at saturation with ice formation in kg/kg

						if ( ( q_c_b <= 0. ) && ( q_i_b <= 0. ) )					q_v_hyp = 0.;
						else 																		q_v_hyp = ( q_c_b * q_Rain + q_i_b * q_Ice ) / ( q_c_b + q_i_b );

						S_c_c.x[ i ][ j ][ k ] = .5 * d_q_c / dt_rain_dim;								// rate of condensating or evaporating water vapour to form cloud water, 0.5 by COSMO

						if ( ( iter_prec >= 3 ) && ( fabs ( q_v_b / q_v_hyp - 1. ) <= 1.e-5 ) )		break;

						q_T = q_v_b + q_c_b + q_i_b;															// total water content

						q_v_b = .5 * ( c.x[ i ][ j ][ k ] + q_v_b );
						q_c_b = .5 * ( cloud.x[ i ][ j ][ k ] + q_c_b );
						q_i_b = .5 * ( ice.x[ i ][ j ][ k ] + q_i_b );

//	if ( ( j == 90 ) && ( k == 180 ) )	cout << " mixed cloud   " << " i = "<< i << " iter_prec   "<< iter_prec << "  p_h = " << p_h << "  r_humid = " << r_humid << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  T = " << T - t_0 << "  d_t = " << d_t << endl << "  q_T = " << q_T * 1000. << "  c = " << c.x[ i ][ j ][ k ] * 1000. << "  q_v_b = " << q_v_b * 1000. << "  cloud = " << cloud.x[ i ][ j ][ k ] * 1000. << "  q_c_b = " << q_c_b * 1000. << "  ice = " << ice.x[ i ][ j ][ k ] * 1000. << "  q_i_b = " << q_i_b * 1000. << endl << "  q_Rain = " << q_Rain * 1000. << "  q_Ice = " << q_Ice * 1000.  << "  q_v_hyp = " << q_v_hyp * 1000. << "  CND = " << CND << "  DEP = " << DEP << "  d_q_v = " << d_q_v * 1000. << "  d_q_c = " << d_q_c * 1000. << "  d_q_i = " << d_q_i * 1000. << endl << endl;

					}																												// iter_prec end

// above assumed tropopause constant temperature t_tropopause
/*
					if ( i >= im_tropopause[ j ] )
					{
						q_v_b = q_c_b = q_i_b = 0.;
					}
*/

					if ( q_v_b <= 0. )	q_v_b = 0.;
					if ( q_c_b <= 0. )	q_c_b = 0.;
					if ( q_i_b <= 0. )		q_i_b = 0.;

					cn.x[ i ][ j ][ k ] = c.x[ i ][ j ][ k ] = q_v_b;
					cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ] = q_c_b;
					icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ] = q_i_b;

					if ( t_Celsius < t_Celsius_2 ) 										cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ] = 0.;
					if ( t_Celsius > 0. ) 													icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ] = 0.;

					t.x[ i ][ j ][ k ] = T / t_0;


				}																												// end ( ( t_Celsius < 0. ) && ( t_Celsius >= t_Celsius_2 ) )





// %%%%%%%%%%%%%%%%%%%%%%%%%%%     ice cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%     ice cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



/*
// ice cloud phase, if t <= -37°C
				if ( t_Celsius <= t_Celsius_2 - 1 )
				{
					q_T = c.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k ];											// total water content
					t_u = t.x[ i ][ j ][ k ] * t_0;																	// in K
					t_Celsius = t_u - t_0;

					p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );							// given in hPa
					hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

//					if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_u ) ) * p_SL;
					if ( i != 0 ) 			p_h = pow ( ( ( t_u - gam * hight * 1.e-2 ) / ( t_u ) ), exp_pressure ) * p_SL;
					else 					p_h = p_SL;

					E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 );											// saturation water vapour pressure for the ice phase in hPa
					q_Ice = ep * E_Ice / ( p_h - E_Ice );													// water vapour amount at saturation with ice formation in kg/kg
					q_Ice_n = q_Ice;

					if ( q_T <= q_Ice )
					{
						c.x[ i ][ j ][ k ] = q_T;																		// total water amount as water vapour
						cloud.x[ i ][ j ][ k ] = 0.;																	// no cloud water available below -37 °C
						ice.x[ i ][ j ][ k ] = 0.;																		// no cloud ice available

//	if ( ( j == 90 ) && ( k == 180 ) )	cout << "ice cloud q_T   " << " i = "<< i << "  p_h = " << p_h << "  r_humid = " << r_humid << "  t_Celsius = " << t_Celsius << "  t = " << t_u - t_0 << "  q_T = " << q_T * 1000. << "  c = " << c.x[ i ][ j ][ k ] * 1000. << "  cl = " << cloud.x[ i ][ j ][ k ] * 1000. << "  ci = " << ice.x[ i ][ j ][ k ] * 1000. << "  q_Ice = " << q_Ice * 1000. << "  q_Ice_n = " << q_Ice_n * 1000. << endl << endl;

					}
					else
					{
						iter_prec = 0;
						while ( iter_prec <= 20 )																// iter_prec may be varied
						{
							iter_prec = iter_prec + 1;
							T_it = ( t_u + lv / cp_l * c.x[ i ][ j ][ k ] - ls / cp_l * q_Ice );
							E_Ice = hp * exp_func ( T_it, 21.8746, 7.66 );											// saturation water vapour pressure for the ice phase in hPa
							q_Ice  = ep * E_Ice / ( p_h - E_Ice );											// water vapour amount at saturation with water formation in kg/kg
							q_Ice = .5 * ( q_Ice_n + q_Ice );

							c.x[ i ][ j ][ k ] = q_Ice;																// water vapour restricted to saturated water vapour amount
							cloud.x[ i ][ j ][ k ] = 0.;																// cloud water amount
							ice.x[ i ][ j ][ k ] = q_T - c.x[ i ][ j ][ k ];										// cloud ice amount

							if ( c.x[ i ][ j ][ k ] < 0. )									c.x[ i ][ j ][ k ] = 0.;
							if ( cloud.x[ i ][ j ][ k ] < 0. )								cloud.x[ i ][ j ][ k ] = 0.;
							if ( ice.x[ i ][ j ][ k ] < 0. )									ice.x[ i ][ j ][ k ] = 0.;

							q_T = c.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k ];

//				if ( ( j == 90 ) && ( k == 180 ) )	cout << "ice cloud iter   " << " i = "<< i << " iter_prec   "<< iter_prec << "  p_h = " << p_h << "  r_humid = " << r_humid << "  T_it = " << T_it - t_0 << "  t = " << t_u - t_0 << "  q_T = " << q_T * 1000. << "  c = " << c.x[ i ][ j ][ k ] * 1000. << "  cl = " << cloud.x[ i ][ j ][ k ] * 1000. << "  ci = " << ice.x[ i ][ j ][ k ] * 1000. << "  q_Ice = " << q_Ice * 1000. << "  q_Ice_n = " << q_Ice_n * 1000. << endl << endl;

							if ( fabs ( q_Ice / q_Ice_n - 1. ) <= 1.e-5 ) 		break;
							q_Ice_n = q_Ice;

						}
					}

//					t.x[ i ][ j ][ k ] = ( t_Celsius + t_0 ) / t_0;
					t.x[ i ][ j ][ k ] = T_it / t_0;
				}																												// end ( t_Celsius < t_Celsius_2 )
*/
			}																													// end i
		}																														// end j
	}																															// end k


/*
	cout << endl << " ***** Results before Snow ***** " << endl << endl;
	cout << endl << " ***** printout of 3D-field water vapour ***** " << endl << endl;
	c.printArray();

	cout << endl << " ***** printout of 3D-field water cloud ***** " << endl << endl;
	cloud.printArray();

	cout << endl << " ***** printout of 3D-field water ice ***** " << endl << endl;
	ice.printArray();

	cout << endl << " ***** printout of 3D-field temperature ***** " << endl << endl;
	t.printArray();

//	cout << endl << " ***** printout of 3D-field epsilon ***** " << endl << endl;
//	epsilon_3D.printArray();

//	cout << endl << " ***** printout of 3D-field radiation_3D ***** " << endl << endl;
//	radiation_3D.printArray();

	cout << endl << " ***** printout of 3D-field p_stat ***** " << endl << endl;
	p_stat.printArray();

//	cout << endl << " ***** printout of 3D-field p_dyn ***** " << endl << endl;
//	p_dyn.printArray();

//	cout << endl << " ***** printout of 2D-field precipitable_water ***** " << endl << endl;
//	precipitable_water.printArray_2D();

//	cout << endl << " ***** printout of 2D-field Precipitation ***** " << endl << endl;
//	Precipitation.printArray_2D();

//	cout << endl << " ***** printout of 2D-field Evaporation_Penman ***** " << endl << endl;
//	Evaporation_Penman.printArray_2D();

//	cout << endl << " ***** printout of 3D-field P_rain ***** " << endl << endl;
//	P_rain.printArray();

//	cout << endl << " ***** printout of 3D-field P_snow ***** " << endl << endl;
//	P_snow.printArray();

//	cout << endl << " ***** printout of 3D-field co2 ***** " << endl << endl;
//	co2.printArray();
*/



/*
// printout for various thermodynamical quantities for the preticipation computations along the equator ( j = 90 )
	cout.precision ( 8 );
				if ( ( i == 5 ) && ( j == 90 ) && ( k == 180 ) )
				{
					cout << endl;
					cout << " i = " << i << "   j = " << j << "   k = " << k << endl << endl;

					cout << " t_h (°C) = " << t_Celsius << "   p_h (hPa) = " << p_h << "   a_h (g/m3) = " << a_h << "   c_h (g/Kg) = " << c.x[ i ][ j ][ k ] * 1000. << "   q_Rain (g/Kg) = " << q_Rain * 1000. << "   e_h (hPa) = " << e_h << "   E_Rain (hPa) = " << E_Rain  << endl << endl;

					cout << " Rain (g/kg) = " << Rain.x[ i ][ j ][ k ] * 1000. << "   Rain_super (g/kg) = " << Rain_super.x[ i ][ j ][ k ] * 1000. << "   Ice (g/kg) = " << Ice.x[ i ][ j ][ k ] * 1000. << "   E_Ice (hPa) = " << E_Ice << "   sat_deficit (hPa) = " << sat_deficit << "   t_dew (°C) = " << t_dew << "   E_Rain_SL (hPa) = " << E_Rain_SL << endl << endl;

					cout << " Evap_Haude (mm/d) = " << Evap_Haude << "   RF_e (%) = " << RF_e  << "   E_Rain_super (hPa) = " << E_Rain_super << "   q_Rain_super (g/Kg) = " << q_Rain_super * 1000. << "   q_h (g/Kg) = " << q_h * 1000. << "   q_SL (g/Kg) = " << q_SL * 1000. << endl << endl;
				}

			}																						// end i
		}																							// end j
	}																								// end k
*/


		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					if ( c.x[ i ][ j ][ k ] < 0. ) 											c.x[ i ][ j ][ k ] = 0.;
					if ( cloud.x[ i ][ j ][ k ] < 0. ) 									cloud.x[ i ][ j ][ k ] = 0.;
					if ( ice.x[ i ][ j ][ k ] < 0. ) 										ice.x[ i ][ j ][ k ] = 0.;
				}
			}
		}

}






void BC_Thermo::Two_Category_Ice_Scheme ( int n, int velocity_iter_max, int RadiationModel, Array &h, Array &c, Array &t, Array &p_stat, Array &cloud, Array &ice, Array &P_rain, Array &P_snow, Array &S_v, Array &S_c, Array &S_i, Array &S_r, Array &S_s, Array &S_c_c )
{
//	Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation distribution formed of rain and snow

// constant coefficients for the transport of cloud water and cloud ice amount vice versa, rain and snow in the parameterization procedures
	a_if = .66;
	c_ac = .24;
	c_rim = 18.6;
	bet_ev = 5.9;
	alf_melt = 7.2e-6;
	bet_melt = bet_dep = 13.;
	alf_if = 1.92e-6;
	alf_cf = 3.97e-5;
	E_cf = 5.0e-3;
	tau_r = 1.e4;
	tau_s = 1.e3;
	a_mc = .08;
	a_mv = .02;
	N_cf_0_surf = 2.e5;
	N_cf_0_6km = 1.e4;
	N_i_0 = 1.e2;																			// in m-3
	t_nuc = 267.15;																		// in K
	t_d = 248.15;																			// in K
	t_hn = 236.15;																		// in K
	t_r_frz = 271.15;																		// in K
	m_i_0 = 1.e-12;																		// in kg
	c_i_dep = 1.3e-5;
	m_i_max = 1.e-9;																	// in kg
	m_s_0 = 3.e-9;																		// in kg
	c_c_au = 4.e-4;																		// in 1/s
	c_i_au = 1.e-3;																		// in 1/s
	c_agg = 10.3;
	c_i_cri = .24;
	c_r_cri = 3.2e-5;
	alf_ev = 1.e-3;
	c_s_dep = 1.8e-2;
	bet_s_dep = 12.3;
//	c_s_dep = 3.2e-2;
//	bet_s_dep = 14.7;
	c_s_melt = 8.43e-5;
	b_s_melt = 12.05;
//	c_s_melt = 1.48e-4;
//	b_s_melt = 14.37;
	a_s_melt = 2.31e3;
	c_r_frz = 3.75e-2;
	a_i_m = 130.;																			// in kg/m3
	a_s_m = .038;																			// in kg/m3
	N_r_0 = 8.e6;																			// in 1/m4
	N_s_0 = 8.e5;																			// in 1/m4
	b_u = .3;
	alf_1 = 5.e-4;
	alf_2 = .011;
	p_ps = .05;
	bet_p = 2.e-3;																			// in s

	t_1 = 253.15;
	t_00 = 236.15;
	t_Celsius_1 = t_1 - t_0;
	t_Celsius_2 = t_00 - t_0;

// rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
// the choosen scheme is a Two Category Ice Scheme
// besides the transport equation for the water vapour exists two equations for the cloud water and the cloud ice transport
// since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation


	for ( int k = 1; k < km-1; k++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			P_rain.x[ im-1 ][ j ][ k ] = 0.;
			P_snow.x[ im-1 ][ j ][ k ] = 0.;
			S_r.x[ im-1 ][ j ][ k ] = 0.;
			S_s.x[ im-1 ][ j ][ k ] = 0.;

			for ( int i = im-2; i >= 0; i-- )
			{
				if ( c.x[ i ][ j ][ k ] < 0. ) 													c.x[ i ][ j ][ k ] = 0.;
				if ( cloud.x[ i ][ j ][ k ] < 0. ) 											cloud.x[ i ][ j ][ k ] = 0.;
				if ( ice.x[ i ][ j ][ k ] < 0. ) 												ice.x[ i ][ j ][ k ] = 0.;

				t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;
				t_u = t.x[ i ][ j ][ k ] * t_0;

				p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );	// given in hPa

				hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

//				if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_u ) ) * p_SL;
				if ( i != 0 ) 			p_h = pow ( ( ( t_u - gam * hight * 1.e-2 ) / ( t_u ) ), exp_pressure ) * p_SL;
				else 					p_h = p_SL;

				r_dry = 100. * p_h / ( R_Air * t_u );																// density of dry air in kg/m³
				r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );							// density of humid air, COSMO version without cloud and ice water, masses negligible in kg/m³

				q_h = c.x[ i ][ j ][ k ];																					// threshold value for water vapour at local hight h in kg/kg

				E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 );										// saturation water vapour pressure for the water phase at t > 0°C in hPa
				E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 );											// saturation water vapour pressure for the ice phase in hPa

				q_Rain  = ep * E_Rain / ( p_h - E_Rain );													// water vapour amount at saturation with water formation in kg/kg
				q_Ice  = ep * E_Ice / ( p_h - E_Ice );															// water vapour amount at saturation with ice formation in kg/kg

				if ( ( t_u >= t_0 ) && ( c_c_au * cloud.x[ i ][ j ][ k ] > 0. ) )					S_c_au = c_c_au * cloud.x[ i ][ j ][ k ];	// cloud water to rain, cloud droplet collection
				else 																								S_c_au = 0.;

				if ( ( t_u <= t_0 ) && ( c_i_au * ice.x[ i ][ j ][ k ] > 0. ) && ( ice.x[ i ][ j ][ k ] >= q_Ice ) )		S_i_au = c_i_au * ( ice.x[ i ][ j ][ k ] - q_Ice );		// cloud ice to snow, cloud ice crystal aggregation
				else 																																S_i_au = 0.;

				S_i_au = 0.;

				if ( t_u <= t_0 ) 																N_i = N_i_0 * exp ( .2 * ( t_0 - t_u ) );
				else 																				N_i = N_i_0;

				if ( t_u <= t_0 )
				{
					if ( ( r_humid * ice.x[ i ][ j ][ k ] / N_i <= m_i_max ) && ( ice.x[ i ][ j ][ k ] > 0. ) )			m_i = r_humid * ice.x[ i ][ j ][ k ] / N_i;
					else 																														m_i = m_i_max;
				}


				if ( ( t_Celsius < 0. ) && ( t_Celsius >= t_Celsius_2 ) )
				{
					if ( c.x[ i ][ j ][ k ] > q_Ice )											// supersaturation
					{
						S_i_dep = c_i_dep * N_i * pow ( m_i, ( 1. / 3. ) ) * ( c.x[ i ][ j ][ k ] - q_Ice );								// supersaturation
					}
					if ( ( c.x[ i ][ j ][ k ] < q_Ice ) && ( - ice.x[ i ][ j ][ k ] / dt_snow_dim > ( c.x[ i ][ j ][ k ] - q_Ice ) / dt_snow_dim ) ) 	S_i_dep = - ice.x[ i ][ j ][ k ] / dt_snow_dim;	// subsaturation
					else 																																								S_i_dep = ( c.x[ i ][ j ][ k ] - q_Ice ) / dt_snow_dim;	// subsaturation
				}
				else 			S_i_dep = 0.;

//				S_i_dep = 0.;																// 9



				S_r.x[ i ][ j ][ k ] = S_c_au;
				S_s.x[ i ][ j ][ k ] = S_i_au;

				if ( h.x[ i ][ j ][ k ] == 1. )
				{
					S_r.x[ i ][ j ][ k ] = 0.;
					S_s.x[ i ][ j ][ k ] = 0.;
				}

				P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ] + ( S_r.x[ i ][ j ][ k ] + S_r.x[ i + 1 ][ j ][ k ] ) * .5 * r_humid * dr;
				P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ] + ( S_s.x[ i ][ j ][ k ] + S_s.x[ i + 1 ][ j ][ k ] ) * .5 * r_humid * dr;

			}
		}
	}


	if ( n >= 2 )
	{
		iter_prec = 0;
		while ( iter_prec <= 5 )															// iter_prec may be varied
		{
			iter_prec = iter_prec + 1;

			for ( int k = 1; k < km-1; k++ )
			{
				for ( int j = 1; j < jm-1; j++ )
				{
					P_rain.x[ im-1 ][ j ][ k ] = 0.;
					P_snow.x[ im-1 ][ j ][ k ] = 0.;

					for ( int i = im-2; i >= 0; i-- )
					{
						t_u = t.x[ i ][ j ][ k ] * t_0;
						t_Celsius = t_u - t_0;

						p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );	// given in hPa

						hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

//						if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_u ) ) * p_SL;
						if ( i != 0 ) 			p_h = pow ( ( ( t_u - gam * hight * 1.e-2 ) / ( t_u ) ), exp_pressure ) * p_SL;
						else 						p_h = p_SL;

						r_dry = 100. * p_h / ( R_Air * t_u );																// density of dry air in kg/m³
						r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );		// density of humid air, COSMO version without cloud and ice water, masses negligible in kg/m³

						q_h = c.x[ i ][ j ][ k ];																					// threshold value for water vapour at local hight h in kg/kg

						E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 );											// saturation water vapour pressure for the water phase at t > 0°C in hPa
						E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 );												// saturation water vapour pressure for the ice phase in hPa

						q_Rain  = ep * E_Rain / ( p_h - E_Rain );														// water vapour amount at saturation with water formation in kg/kg
						q_Ice  = ep * E_Ice / ( p_h - E_Ice );															// water vapour amount at saturation with ice formation in kg/kg


// ice and snow average size
						if ( t_u <= t_0 ) 														N_i = N_i_0 * exp ( .2 * ( t_0 - t_u ) );
						else 																		N_i = N_i_0;

						if ( t_u <= t_0 )
						{
							if ( ( r_humid * ice.x[ i ][ j ][ k ] / N_i <= m_i_max ) && ( ice.x[ i ][ j ][ k ] > 0. ) )			m_i = r_humid * ice.x[ i ][ j ][ k ] / N_i;
							else 																																m_i = m_i_max;
						}


// nucleation and depositional growth of cloud ice
						if ( ice.x[ i ][ j ][ k ] == 0. )
						{
							if ( ( t_hn <= t_u ) && ( t_u < t_d ) && ( c.x[ i ][ j ][ k ] >= q_Ice ) ) 				S_nuc = m_i_0 / ( r_humid * dt_dim ) * N_i;		// nucleation of cloud ice
							if ( ( t_d <= t_u ) && ( t_u <= t_nuc ) && ( c.x[ i ][ j ][ k ] >= q_Rain ) )			S_nuc = m_i_0 / ( r_humid * dt_dim ) * N_i;		// nucleation of cloud ice
						}
						else 																																											S_nuc = 0.;

//						S_nuc = 0.;															// 14


						if ( ( t_u < t_hn ) && ( cloud.x[ i ][ j ][ k ] > 0. ) ) 		S_c_frz = cloud.x[ i ][ j ][ k ] / dt_rain_dim;		//nucleationof cloud ice due to freezing of cloud water
						else 																		S_c_frz = 0.;

//						S_c_frz = 0.;															// 16


						if ( ( t_Celsius < 0. ) && ( t_Celsius >= t_Celsius_2 ) )
						{	
							if ( c.x[ i ][ j ][ k ] > q_Ice )								// supersaturation
							{
								S_i_dep = c_i_dep * N_i * pow ( m_i, ( 1. / 3. ) ) * ( c.x[ i ][ j ][ k ] - q_Ice );								// supersaturation
							}

							if ( ( c.x[ i ][ j ][ k ] < q_Ice ) && ( - ice.x[ i ][ j ][ k ] / dt_snow_dim > ( c.x[ i ][ j ][ k ] - q_Ice ) / dt_snow_dim ) ) 	S_i_dep = - ice.x[ i ][ j ][ k ] / dt_snow_dim;	// subsaturation
							else 																																												S_i_dep = ( c.x[ i ][ j ][ k ] - q_Ice ) / dt_snow_dim;	// subsaturation
						}
						else 			S_i_dep = 0.;

//						S_i_dep = 0.;														// 9


// autoconversion processes
						if ( ( t_u >= t_0 ) && ( c_c_au * cloud.x[ i ][ j ][ k ] > 0. ) )					S_c_au = c_c_au * cloud.x[ i ][ j ][ k ];	// cloud water to rain, cloud droplet collection
						else 																								S_c_au = 0.;

						if ( ( t_u <= t_0 ) && ( c_i_au * ice.x[ i ][ j ][ k ] > 0. ) && ( ice.x[ i ][ j ][ k ] >= q_Ice ) )		S_i_au = c_i_au * ( ice.x[ i ][ j ][ k ] - q_Ice );		// cloud ice to snow, cloud ice crystal aggregation
						else 																																S_i_au = 0.;

						if ( ( t_u < t_0 ) && ( c.x[ i ][ j ][ k ] > q_Ice ) )									S_d_au = S_i_dep / ( 1.5 * ( pow ( ( m_s_0 / m_i ),  ( 2. / 3. ) ) - 1. ) );		// depositional growth of cloud ice
						else 																								S_d_au = 0.;

//						S_c_au = 0.;															// 1
//						S_i_au = 0.;															// 5
//						S_d_au = 0.;															// 17


// collection mechanism
						if ( t_u > t_0 ) 														S_ac = c_ac * cloud.x[ i ][ j ][ k ] * pow ( P_rain.x[ i ][ j ][ k ], ( 7. / 9. ) ); 	// accreation rate from depletion of cloud water due to collection by all rain drops
						else																		S_ac = 0.;										// accreation rate from depletion of cloud water due to collection by all rain drops

						if ( t_u < t_0 ) 														S_rim = c_rim * cloud.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];
						else																		S_rim = 0.;										// riming rate of snow mass due to collection of supercooled cloud droplets
																																								// by falling snow particles

						if ( t_u >= t_0 ) 														S_shed = c_rim * cloud.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];
						else 																		S_shed = 0.;									// rate of water shed by melting wet snow particles
																																						// collecting cloud droplets to produce rain

//						S_ac = 0.;																// 2
//						S_rim = 0.;															// 7
//						S_shed = 0.;															// 11



						if ( t_u < t_0 )
						{
							S_agg = c_agg * ice.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];														// collection of cloud ice by snow particles

							S_i_cri = c_i_cri * ice.x[ i ][ j ][ k ] * pow ( P_rain.x[ i ][ j ][ k ], ( 7. / 9. ) );								// decrease in cloud ice mass due to collision/coalescense interaction with raindrops

							S_r_cri = c_r_cri * ice.x[ i ][ j ][ k ] / m_i * pow ( P_rain.x[ i ][ j ][ k ], ( 13. / 9. ) );					// decrease of rainwater due to freezing resulting from collection of ice crystals
						}
						else
						{
							S_agg = 0.;
							S_i_cri = 0.;
							S_r_cri = 0.;
						}

//						S_agg = 0.;															// 6
//						S_i_cri = 0.;															// 13
//						S_r_cri = 0.;															// 12


// diffusional growth of rain and snow

						if ( t_u >= t_0 )			S_ev = alf_ev * ( 1. + bet_ev * pow ( P_rain.x[ i ][ j ][ k ], ( 1. / 6. ) ) ) * ( q_Rain - c.x[ i ][ j ][ k ] ) * pow ( P_rain.x[ i ][ j ][ k ], ( 4. / 9. ) );
																									// evaporation of rain due to water vapour diffusion
						else 							S_ev = 0.;

//						S_ev = 0.;																// 3


						if ( t_u < t_0 ) 			S_s_dep = c_s_dep * ( 1. + bet_s_dep * pow ( P_snow.x[ i ][ j ][ k ], ( 5. / 26. ) ) ) * ( c.x[ i ][ j ][ k ] - q_Ice ) * pow ( P_snow.x[ i ][ j ][ k ], ( 8. / 13. ) );
																									// deposition/sublimation of snow 
						else 							S_s_dep = 0.;

//						S_s_dep = 0.;														// 8


// melting and freezing
						if ( ( t_u > t_0 ) && ( ice.x[ i ][ j ][ k ] > 0. ) ) 			S_i_melt = ice.x[ i ][ j ][ k ] / dt_snow_dim; // cloud ice particles melting to cloud water
						else 																		S_i_melt = 0.;

//						S_i_melt = 0.;														// 15


						if ( t_u > t_0 )
						{
							p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );		// given in hPa

							hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

//							p_t_in = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_u ) ) * p_SL;	// given in hPa
							p_t_in = pow ( ( ( t_u - gam * hight * 1.e-2 ) / ( t_u ) ), exp_pressure ) * p_SL;	// given in hPa

							E_Rain_t_in = hp * exp_func ( t_0, 17.2694, 35.86 );								// saturation water vapour pressure for the water phase at t = 0°C in hPa
							q_Rain_t_in = ep * E_Rain_t_in / ( p_t_in - E_Rain_t_in );							// water vapour amount at saturation with water formation in kg/kg

//							S_s_melt = c_s_melt * ( 1. + b_s_melt * pow ( P_snow.x[ i ][ j ][ k ], ( 5. / 26. ) ) ) * ( ( t_u - t_0 ) + a_s_melt * ( c.x[ i ][ j ][ k ] - q_Rain_t_in ) ) * pow ( P_snow.x[ i ][ j ][ k ], ( 8. / 13. ) );																// melting rate of snow to form rain
							S_s_melt = 50. * c_s_melt * ( 1. + b_s_melt * pow ( P_snow.x[ i ][ j ][ k ], ( 5. / 26. ) ) ) * ( ( t_u - t_0 ) + a_s_melt * ( c.x[ i ][ j ][ k ] - q_Rain_t_in ) ) * pow ( P_snow.x[ i ][ j ][ k ], ( 8. / 13. ) );								// melting rate of snow to form rain
						}
						else 																	S_s_melt = 0.;

//					S_s_melt = 0.;														// 10


					if ( t_r_frz -  t_u > 0. ) 											S_r_frz = c_r_frz * pow ( ( t_r_frz -  t_u ), ( 3. / 2. ) ) * pow ( P_rain.x[ i ][ j ][ k ], ( 3. / 2. ) );
					else 																		S_r_frz = 0.;

//					S_r_frz = 0.;															// 4


// sinks and sources
					S_v.x[ i ][ j ][ k ] = - S_c_c.x[ i ][ j ][ k ] + S_ev - S_i_dep - S_s_dep - S_nuc;
					S_c.x[ i ][ j ][ k ] = S_c_c.x[ i ][ j ][ k ] - S_c_au - S_ac - S_c_frz + S_i_melt - S_rim - S_shed;
					S_i.x[ i ][ j ][ k ] = S_nuc + S_c_frz + S_i_dep - S_i_melt - S_i_au - S_d_au - S_agg - S_i_cri;
					S_r.x[ i ][ j ][ k ] = S_c_au + S_ac - S_ev + S_shed - S_r_cri - S_r_frz + S_s_melt;
					S_s.x[ i ][ j ][ k ] = S_i_au + S_d_au + S_agg + S_rim + S_s_dep + S_i_cri + S_r_cri + S_r_frz - S_s_melt;


// rain and snow integration
					if ( h.x[ i ][ j ][ k ] == 1. )
					{
						S_r.x[ i ][ j ][ k ] = 0.;
						S_s.x[ i ][ j ][ k ] = 0.;
					}

					P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ] + ( S_r.x[ i ][ j ][ k ] + S_r.x[ i + 1 ][ j ][ k ] ) * .5 * r_humid * dr;
					P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ] + ( S_s.x[ i ][ j ][ k ] + S_s.x[ i + 1 ][ j ][ k ] ) * .5 * r_humid * dr;

					if ( P_rain.x[ i ][ j ][ k ] < 0. )										P_rain.x[ i ][ j ][ k ] = 0.;
					if ( P_snow.x[ i ][ j ][ k ] < 0. )										P_snow.x[ i ][ j ][ k ] = 0.;

					if ( c.x[ i ][ j ][ k ] < 0. ) 											c.x[ i ][ j ][ k ] = 0.;
					if ( cloud.x[ i ][ j ][ k ] < 0. ) 										cloud.x[ i ][ j ][ k ] = 0.;
					if ( ice.x[ i ][ j ][ k ] < 0. ) 											ice.x[ i ][ j ][ k ] = 0.;


					cout.precision ( 8 );
					cout.setf ( ios::fixed );
//					cout.precision ( 2 );
//					cout.setf ( std::ios::scientific );

//					if ( ( j == 90 ) && ( k == 180 ) )	cout << " iter_prec = " << iter_prec << "  i = " << i << "  p_h = " << p_h << "  p_SL = " << p_SL << "  r_humid = " << r_humid << "  coeff_P = " << coeff_P << "  t_cretaceous = " << t_cretaceous << "  t_Celsius = " << t_Celsius << "  t = " << t_u - t_0 << "  q_T = " << ( c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k ] ) * 1000. << "  c = " << c.x[ i ][ j ][ k ] * 1000. << "  cl = " << cloud.x[ i ][ j ][ k ] * 1000. << "  ci = " << ice.x[ i ][ j ][ k ] * 1000. << "  q_Rain = " << q_Rain * 1000. << "  q_Ice = " << q_Ice * 1000. << "  S_i_dep = " << S_i_dep << "  N_i = " << N_i << "  m_i = " << m_i * 1000000. << "  S_nuc = " << S_nuc << "  S_c_frz = " << S_c_frz << "  S_c_au = " << S_c_au << "  S_i_au = " << S_i_au << "  S_d_au = " << S_d_au << "  S_ac = " << S_ac << "  S_rim = " << S_rim << "  S_shed = " << S_shed << "  S_agg = " << S_agg << "  S_i_cri = " << S_i_cri << "  S_r_cri = " << S_r_cri << "  S_ev = " << S_ev << "  S_s_dep = " << S_s_dep << "  S_i_melt = " << S_i_melt << "  S_s_melt = " << S_s_melt << "  S_r_frz = " << S_r_frz << "  S_v = " << S_v.x[ i ][ j ][ k ] << "  S_c = " << S_c.x[ i ][ j ][ k ] << "  S_i = " << S_i.x[ i ][ j ][ k ] << "  S_r = " << S_r.x[ i ][ j ][ k ] << "  S_s = " << S_s.x[ i ][ j ][ k ] << "  S_c_c = " << S_c_c.x[ i ][ j ][ k ] << "  P_rain = " << P_rain.x[ i ][ j ][ k ] << "  P_snow = " << P_snow.x[ i ][ j ][ k ] << endl;

					}																					// end i RainSnow
				}																						// end j
			}																							// end k
		}																								// end iter_prec
	}																								// end n

}







void BC_Thermo::IC_Temperature_WestEastCoast ( Array &h, Array &t )
{
// initial conditions for the temperature close to coast sides to damp out shades of preceeding timeslices

//	j_grad = 5;																				// extension for temperature change in zonal direction
	j_grad = 7;																				// extension for temperature change in zonal direction
//	k_grad = 5;																				// extension for temperature change in longitudinal direction
	k_grad = 7;																				// extension for temperature change in longitudinal direction



// search for north coasts to smooth the temperature

// northern and southern hemisphere: north coast

	j_water = 0;																			// somewhere on water
	flip = 0;																					// somewhere on water

	for ( int k = 1; k < km-1; k++ )											// outer loop: longitude
	{
		for ( int j = j_grad; j < jm-1; j++ )									// inner loop: latitude
		{
			if ( h.x[ 0 ][ j ][ k ] == 0. )												// if somewhere on water
			{
				j_water = 0;																// somewhere on water: j_water = 0
				flip = 0;																		// somewhere on water: flip = 0
			}
			else j_water = 1;															// first time on land

			if ( ( flip == 0 ) && ( j_water == 1 ) )							// on water closest to land
			{
				ll = j - j_grad;																// starting point of temperature smoothing

				for ( int l = ll; l < j; l++ )		t.x[ 0 ][ l ][ k ] = t.x[ 0 ][ ll ][ k ];		// replacement of temperature values

				flip = 1;																		// somewhere on land: flip = 1
			}
		}																							// end of latitudinal loop
		flip = 0;																				// somewhere on water: flip = 0
	}																								// end of longitudinal loop





// northern and southern hemisphere: south coast

	j_water = 0;																			// on water closest to coast
	j_sequel = 1;																			// on solid ground

	for ( int k = 1; k < km-1; k++ )											// outer loop: latitude
	{
		for ( int j = 0; j < jm - j_grad; j++ )									// inner loop: longitude
		{
			if ( h.x[ 0 ][ j ][ k ] == 1. ) j_sequel = 0;						// if solid ground: j_sequel = 0

			if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( j_sequel == 0 ) ) j_water = 0;	// if water and and j_sequel = 0 then is water closest to coast
			else j_water = 1;															// somewhere on water

			if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( j_water == 0 ) )		// if water is closest to coast, change of velocity components begins
			{
				ll = j + j_grad;															// starting point of temperature smoothing

				for ( int l = ll; l > j; l-- )		t.x[ 0 ][ l ][ k ] = t.x[ 0 ][ ll ][ k ];		// replacement of temperature values


				j_sequel = 1;																// looking for another south coast
			}
		}																							// end of longitudinal loop
		j_water = 0;																		// starting at another latitude
	}																								// end of latitudinal loop




// northern and southern hemisphere: east coast
	k_water = 0;																			// on water closest to coast
	k_sequel = 1;																			// on solid ground

	for ( int j = 1; j < jm-1; j++ )												// outer loop: latitude
	{
		for ( int k = k_grad; k < km-k_grad; k++ )						// inner loop: longitude
		{
			if ( h.x[ 0 ][ j ][ k ] == 1. ) k_sequel = 0;						// if solid ground: k_sequel = 0

			if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( k_sequel == 0 ) ) k_water = 0;		 // if water and and k_sequel = 0 then water lies closest to a coast
			else k_water = 1;															// somewhere on water

			if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( k_water == 0 ) )		// if water is closest to coast, change of velocity components begins
			{
				ll = k + k_grad;															// starting point of temperature smoothing

				for ( int l = ll; l > k; l-- )		t.x[ 0 ][ j ][ l ] = t.x[ 0 ][ j ][ ll ];		// replacement of temperature values

				k_sequel = 1;																// looking for another east coast
			}
		}																							// end of longitudinal loop
		k_water = 0;																		// starting at another longitude
	}																								// end of latitudinal loop





// northern and southern hemisphere: west coast

	k_water = 0;																			// somewhere on water
	flip = 0;																					// somewhere on water

	for ( int j = 1; j < jm-1; j++ )												// outer loop: latitude
	{
		for ( int k = k_grad; k < km-1; k++ )								// inner loop: longitude
		{
			if ( h.x[ 0 ][ j ][ k ] == 0. )												// if somewhere on water
			{
				k_water = 0;																// somewhere on water: k_water = 0
				flip = 0;																		// somewhere on water: flip = 0
			}
			else k_water = 1;															// first time on land

			if ( ( flip == 0 ) && ( k_water == 1 ) )							// on water closest to land
			{
				ll = k - k_grad;															// starting point of temperature smoothing

				for ( int l = ll; l < k; l++ )		t.x[ 0 ][ j ][ l ] = t.x[ 0 ][ j ][ ll ];		// replacement of temperature values

				flip = 1;																		// somewhere on land: flip = 1
			}
		}																							// end of longitudinal loop
		flip = 0;																				// somewhere on water: flip = 0
	}																								// end of latitudinal loop



}



double BC_Thermo::exp_func ( double &T_K, const double &co_1, const double &co_2 )
{
	return exp ( co_1 * ( T_K - 273.15 ) / ( T_K - co_2 ) );						// temperature in °K
}



 
double BC_Thermo::out_t_cretaceous (  ) const
{
	return t_cretaceous;
}

double BC_Thermo::out_t_cret_cor (  ) const
{
	return t_cret_cor;
}

double BC_Thermo::out_co2 (  ) const
{
	return co2_cretaceous;
}

