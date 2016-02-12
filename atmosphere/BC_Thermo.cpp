/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
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

#include "BC_Thermo.h"

using namespace std;


BC_Thermo::BC_Thermo ( int im, int jm, int km, Array &t, Array &c, Array &aux_v, Array &aux_w )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;

	c43 = 4./3.;
	c13 = 1./3.;

	pi180 = 180./M_PI;


//	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 8 );
	cout.setf ( ios::fixed );

// array "im_tropopause" for configuring data due to latitude dependent tropopause

	im_tropopause = new int[ jm ];

	for ( int l = 0; l < jm; l++ )
	{
		im_tropopause[ l ] = 0;
//		cout << im_tropopause[ l ] << endl;
	}


// array "jm_temp_asym" for configuring data due to latitude dependent tropopause

	jm_temp_asym = new double[ jm ];

	for ( int l = 0; l < jm; l++ )
	{
		jm_temp_asym[ l ] = 0;
//		cout << jm_temp_asym[ l ] << endl;
	}

}



BC_Thermo::~BC_Thermo()
{
	delete [  ]im_tropopause;
	delete [  ]jm_temp_asym;
}





void BC_Thermo::BC_RadiationBalance_comp ( double t_0, double ik, double sigma, double albedo, double epsilon, double r_0_water_vapour, double R_WaterVapour, Array_2D &Precipitation, Array_2D &Ik, Array_2D &Radiation_Balance, Array_2D &Radiation_Balance_atm, Array_2D &Radiation_Balance_bot, Array_2D &temp_eff_atm, Array_2D &temp_eff_bot, Array_2D &t_j, Array &t )
{
// class element for the computation of the radiation balance
// computation of the local temperature based on the radiation balance

	cout.precision ( 6 );
	cout.setf ( ios::fixed );

	pi180 = 180./M_PI;

//	j_half = ( jm -1 ) / 2 + 30;														// position of the sun at 30°S ( 90 + 30 = 120 )
	j_half = ( jm -1 ) / 2;																	// position of the sun at 0°S
	j_max = jm - 1;

	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

	k_half = ( km -1 ) / 2;																// position of the sun at 0° oder 180° ( Greenwich )
	k_max = km - 1;

	d_k_half = ( double ) k_half;
	d_k_max = ( double ) k_max;


	t_eff_earth = 255.;																	// effectiv radiation temperature of 255 K compares to -18 °C for 30% albedo
																									// comparable with measurements of temperature in 5000 m hight and pressure 550 hPa
																									// sigma * t-earth ** 4 = ( 1- albedo ) Ik / 4 = 239 W/m2, valid as reference value

//	j_sun = 30;																				// lateral sun location, 30 = 60°N
	j_sun = 0;																				// equatorial sun location

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
//			Ik.y[ j ][ k ] = ik * sin ( ( ( double ) j - j_sun ) / pi180 ) * sin ( ( ( double ) k - 90 ) / pi180 ) / 4.;	//	solar short wave radiation on a point location
			Ik.y[ j ][ k ] = ik * sin ( ( ( double ) j - j_sun ) / pi180 ) / 4.;	// solar short wave radiation as zonal distribution
			if ( Ik.y[ j ][ k ] < 40. ) Ik.y[ j ][ k ] = 40.;

			epsilon = .71;
//			epsilon = 2. * ( 1. - 1. / pow ( ( t_j.y[ j ][ k ] * t_0 / t_eff_earth ), 4. ) ); // bottom temperature must be known in advance

			Radiation_Balance_bot.y[ j ][ k ] = ( 1. - albedo ) / ( 1. - epsilon / 2. ) * Ik.y[ j ][ k ];
			Radiation_Balance_atm.y[ j ][ k ] = ( 1. - albedo ) / ( 2. - epsilon ) * Ik.y[ j ][ k ];

			temp_eff_bot.y[ j ][ k ] = pow ( Radiation_Balance_bot.y[ j ][ k ] / sigma, 1. / 4. );
			temp_eff_atm.y[ j ][ k ] = pow ( Radiation_Balance_atm.y[ j ][ k ] / sigma, 1. / 4. );

			epsilon = 2. * ( 1. - 1. / pow ( ( temp_eff_bot.y[ j ][ k ] / t_eff_earth ), 4. ) ); // bottom temperature must be known in advance

			e = ( r_0_water_vapour * R_WaterVapour * temp_eff_bot.y[ j ][ k ] ) * .01;	// water vapour pressure at the surface in hPa

			D = Ik.y[ j ][ k ] * .26;														// direct sun radiation, short-wave (Häckel)
			H = Ik.y[ j ][ k ] * .29;														// diffusive sky radiation, short-wave (Häckel)
			G = D + H;																		// total solar radiation from above = global radiation, short-wave
			R_short = G * albedo;														// total reflected solar radiation, short-wave
			Q_short = G - R_short;													// short-wave radiation balance

//			AG = ( .594 + .0416 * sqrt ( e ) ) * sigma * pow ( temp_eff_bot.y[ j ][ k ], 4. ); // downward terrestrial radiation, long-wave (Häckel)
//			A = epsilon * sigma * pow ( temp_eff_bot.y[ j ][ k ], 4. );	// outgoing long-wave radiation of the surface, long-wave
			AG = Ik.y[ j ][ k ] * 1.14;													// downward terrestrial radiation, long-wave (Häckel)
			A = Ik.y[ j ][ k ] * .95;														// outgoing long-wave radiation of the surface, long-wave (Häckel)
			R_long = AG * .05;															// total reflected terrestrial radiation, long-wave (Häckel)
			Q_long = AG - A - R_long;												// long-wave radiation balance

			Radiation_Balance.y[ j ][ k ] = Q_short + Q_long;			// total radiation balance, short- and long-wave



// printout for various radiation quantities for the terrestrial radiation balance
/*
	cout.precision ( 3 );

			cout << endl;
			cout << "   j = " << j << "   k = " << k << endl << endl;

			cout << "__________     radiation data in W/m²    _______________" << endl << endl;

			cout << " t_SL (°C) = " << t_j.y[ j ][ k ] * t_0 - t_0 << "    Ik / 4 (W/m²) = " << Ik.y[ j ][ k ] << "   Radiation_Balance_bot = " << Radiation_Balance_bot.y[ j ][ k ] << "   Radiation_Balance_atm = " << Radiation_Balance_atm.y[ j ][ k ] << "   temp_eff_bot = " << temp_eff_bot.y[ j ][ k ] << "   temp_eff_atm = " << temp_eff_atm.y[ j ][ k ] << endl << endl;

			cout << " D  = " << D << "   H = " << H << "   G = " << G << "   R_short = " << R_short << "   Q_short = " << Q_short << "   epsilon = " << epsilon << "   albedo = " << albedo << endl << endl;

			cout << " AG = " << AG << "   A = " << A << "               R_long  = " << R_long << "   Q_long  = " << Q_long << "   e = " << e << "   Radiation_Balance = " << Radiation_Balance.y[ j ][ k ] << endl << endl << endl;

			cout << "__________     radiation data in % are based on the local extra terrestrial sun radiaton Ik in W/m²    _______________" << endl << endl;

			cout << " t_h (°C) = " << t_j.y[ j ][ k ] * t_0 - t_0 << "    Ik / 4 (%) = " << Ik.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   Radiation_Balance_bot = " << Radiation_Balance_bot.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   Radiation_Balance_atm = " << Radiation_Balance_atm.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << "   temp_eff_bot = " << temp_eff_bot.y[ j ][ k ] << "   temp_eff_atm = " << temp_eff_atm.y[ j ][ k ] << endl << endl;

			cout << " D  = " << D / Ik.y[ j ][ k ] * 100. << "   H = " << H / Ik.y[ j ][ k ] * 100. << "   G = " << G / Ik.y[ j ][ k ] * 100. << "   R_short = " << R_short / Ik.y[ j ][ k ] * 100. << "   Q_short = " << Q_short / Ik.y[ j ][ k ] * 100. << "   epsilon = " << epsilon << "   albedo = " << albedo << endl << endl;

			cout << " AG = " << AG / Ik.y[ j ][ k ] * 100. << "   A = " << A / Ik.y[ j ][ k ] * 100. << "               R_long  = " << R_long / Ik.y[ j ][ k ] * 100. << "   Q_long  = " << Q_long / Ik.y[ j ][ k ] * 100. << "   e (hPa) = " << e << "   Radiation_Balance = " << Radiation_Balance.y[ j ][ k ] / Ik.y[ j ][ k ] * 100. << endl << endl << endl;

	cout.precision ( 6 );
*/

		}
	}
}






void BC_Thermo::BC_Temperature ( int i_max, int i_beg, int RadiationFluxDensity, int Ma, int Ma_max, int Ma_max_half, int sun_position_lat, int sun_position_lon, int declination, double sun, double ep, double hp, double t_0, double p_0, double t_land_plus, double t_cretaceous_max, double t_average, double t_equator, double t_pole, double t_tropopause, Array_2D &t_j, Array &h, Array &t, Array &p )
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

	t_coeff = t_pole - t_equator;

// temperature-distribution by Ruddiman approximated by a parabola
	t_Cretaceous_coeff = t_cretaceous_max / ( ( double ) Ma_max_half - ( double ) ( Ma_max_half * Ma_max_half / Ma_max ) );   // in °C
	if ( Ma == 0 ) 	t_Cretaceous_coeff = 0.;
	t_Cretaceous = t_Cretaceous_coeff * ( double ) ( - ( Ma * Ma ) / Ma_max + Ma );   // in °C


	cout.precision ( 3 );

	time_slice_comment = "      time slice of Cretaceous-AGCM:";
	time_slice_number = " Ma = ";
	time_slice_unit = " million years";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << time_slice_comment << resetiosflags ( ios::left ) << setw ( 6 ) << fixed << setfill ( ' ' ) << time_slice_number << setw ( 3 ) << Ma << setw ( 12 ) << time_slice_unit << endl << endl;


	temperature_comment = "      temperature increase at cretaceous times: ";
	temperature_gain = " t cretaceous";
	temperature_modern = "      mean temperature at modern times: ";
	temperature_average = " t modern";
	temperature_unit =  "°C ";

	cout << endl << setiosflags ( ios::left ) << setw ( 50 ) << setfill ( '.' ) << temperature_comment << resetiosflags ( ios::left ) << setw ( 12 ) << temperature_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << t_Cretaceous << setw ( 5 ) << temperature_unit << endl << setw ( 50 ) << setfill ( '.' )  << setiosflags ( ios::left ) << temperature_modern << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average << setw ( 5 ) << temperature_unit << endl;


	t_Cretaceous = ( t_Cretaceous + t_average + t_0 ) / t_0 - ( ( t_average + t_0 ) / t_0 );    // non-dimensional

// if Ma = 0 ( modern times ), NASA-temperature distribution is used

// mean zonal parabolic temperature distribution
	if ( Ma > 0 )
	{
		if ( sun == 0 )
		{
			for ( int k = 0; k < km; k++ )
			{
				for ( int j = 0; j < jm; j++ )
				{
					d_j = ( double ) j;
					t.x[ 0 ][ j ][ k ] = t_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + t_pole + t_Cretaceous;

					if ( h.x[ 0 ][ j ][ k ]  == 1. ) 
					{
						t.x[ 0 ][ j ][ k ] = t_j.y[ j ][ k ] = t.x[ 0 ][ j ][ k ] + t_land_plus;
					}
				}
			}                                                           
		}																							// end parabolic temperature distribution
	}																								// end of while (  Ma > 0 )



// temperatur distribution at a prescribed sun position
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

		a = ( t_equator - t_pole ) / ( ( ( j_par_f * j_par_f ) - ( j_pol_f * j_pol_f ) ) - 2. * j_par_f * ( j_par_f - j_pol_f ) );
		b = - 2. * a * j_par_f;
		cc = t_equator + a * j_par_f * j_par_f;
		j_d = sqrt ( ( cc - t_pole ) / a );
		d = 2. * a * j_d + b;
		t_d = d * j_d + t_pole;
		e = t_pole;


// asymmetric temperature distribution from pole to pole for  j_d  maximum temperature ( linear equation + parabola )
		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				d_j = ( double ) j;
				if ( d_j <= j_d )
				{
					t.x[ 0 ][ j ][ k ] = t_j.y[ j ][ k ] = d * d_j + e + t_Cretaceous;
				}
				if ( d_j > j_d )
				{
					t.x[ 0 ][ j ][ k ] = t_j.y[ j ][ k ] = a * d_j * d_j + b * d_j + cc + t_Cretaceous;
				}
			}
		}


// transfer of zonal constant temperature into a 1D-temperature field
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

				a = ( jm_temp_asym[ j ] - t_360 ) / ( ( ( k_par_f * k_par_f ) - ( k_pol_f * k_pol_f ) ) - 2. * k_par_f * ( k_par_f - k_pol_f ) );
				b = - 2. * a * k_par_f;
				cc = jm_temp_asym[ j ] + a * k_par_f * k_par_f;

				t.x[ 0 ][ j ][ k ] = t_j.y[ j ][ k ] = a * d_k * d_k + b * d_k + cc;
			}
		}
	}																								// temperatur distribution at a prescribed sun position




// temperature decreasing approaching the tropopause, above constant temperature following Standard Atmosphere
	for ( int j = 0; j < jm; j++ )
	{
		d_i_max = ( double ) im_tropopause[ j ];

		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i <= im_tropopause[ j ]; i++ )
			{
				{
					d_i = ( double ) i;
					t.x[ i ][ j ][ k ] = ( t_tropopause - t.x[ 0 ][ j ][ k ] ) / d_i_max * d_i + t.x[ 0 ][ j ][ k ];				// linear temperature decay up to tropopause
//					t.x[ i ][ j ][ k ] = pow ( p.x[ i ][ j ][ k ], .28596 );																// temperature from adiabatic formula
//					t.x[ i ][ j ][ k ] = ( t.x[ 0 ][ j ][ k ] * t_0 - .9433 * d_i * 5. ) / t_0;											// temperature from dry adiabatic formula ( correct ))
				}
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

}







void BC_Thermo::TropopauseLocation ( int i_max, int i_beg )
{
// parabolic tropopause location distribution from pole to pole assumed

	j_half = ( jm -1 ) / 2;

	d_j_half = ( double ) j_half;

	trop_coeff = ( double ) ( i_beg - i_max );

// computation of the tropopause from pole to pole

	for ( int j = 0; j < jm; j++ )
	{
		d_j = ( double ) j;
		im_tropopause[ j ] = ( int ) ( trop_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) ) + i_beg;
	}
}







void BC_Thermo::BC_WaterVapour ( int i_max, int i_beg, double ep, double hp, double t_0, double c_0, double p_0, double c_land_minus, double c_ocean_minus, double c_tropopause, Array_2D &c_j, Array &h, Array &t, Array &p, Array &c, Array_2D &t_j, Array_2D &precipitation_j )
{
// initial and boundary conditions of water vapour on water and land surfaces
// parabolic water vapour distribution from pole to pole accepted

// maximum water vapour content on water surface at equator c_equator = 1.04 compares to 0.04 volume parts
// polar water vapour contents on water surface at North and South Pole c_pol = 1.004 compares to 0.004 volume parts
// minimum water vapour at tropopause c_tropopause = 0.0 compares to 0.0 volume parts
// value 1.0 stands for the maximum value of 35 g/kg water vapour

	d_i_max = ( double ) i_max;

	d_j_half = ( double ) j_half;
	d_j_max = ( double ) j_max;

// water vapour contents computed by Clausius-Clapeyron-formula
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			d_j = ( double ) j;
			if ( h.x[ 0 ][ j ][ k ]  == 0. ) 
			{
//				e_h = ( r_0_water_vapour * R_WaterVapour * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;								// water vapour pressure at local hight h
//				c_j.y[ j ][ k ] = ep * e_h / p_0;
				c_j.y[ j ][ k ]  = hp * ep *exp ( 17.0809 * ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) ) / p_0; // saturation of relative water vapour 
				c_j.y[ j ][ k ] = c_ocean_minus * c_j.y[ j ][ k ];				// relativ water vapour contents on ocean surface reduced by factor
				if ( c_j.y[ j ][ k ] >= .04 )	c_j.y[ j ][ k ] = .04;				// 40 g is the maximum water vapour in dry air without outfall
				c.x[ 0 ][ j ][ k ] = c_j.y[ j ][ k ];										// relativ water vapour contents increased by by factor
			}

			if ( h.x[ 0 ][ j ][ k ]  == 1. ) 
			{
//				e_h = ( r_0_water_vapour * R_WaterVapour * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;								// water vapour pressure at local hight h
//				c_j.y[ j ][ k ] = ep * e_h / p_0;
				c_j.y[ j ][ k ]  = hp * ep * exp ( 17.0809 * ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 ) ) ) / p_0;
				c_j.y[ j ][ k ] = c_land_minus * c_j.y[ j ][ k ];					// relativ water vapour contents on land reduced by factor
				if ( c_j.y[ j ][ k ] >= .04 )	c_j.y[ j ][ k ] = .04;				// 40 g is the maximum water vapour in dry air without outfall
				c.x[ 0 ][ j ][ k ] = c_j.y[ j ][ k ];										// relativ water vapour contents increased by factor
			}
		}
	}





// water vapour distribution decreasing approaching tropopause, above no water vapour
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
					if ( i <= im_tropopause[ j ] - 3 )								// reduction by 5 because water vapour is finished at latest 12 km
					{
						d_i_max = ( double ) ( im_tropopause[ j ] - 3 );		// reduction by 5 because water vapour is finished at latest 12 km
						d_i = ( double ) i;

						c.x[ i ][ j ][ k ] = c_j.y[ j ][ k ] - ( c_tropopause - c_j.y[ j ][ k ] ) * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) );	// radial distribution approximated by a parabola ( Weischet )
					}
					else 			c.x[ i ][ j ][ k ] = c_tropopause;
			}
		}
	}



/*
// water vapour at the surface reduced by 20% to achieve evaporation at the surface for a positive water vapour gradient
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			c.x[ i ][ j ][ k ] = c_j.y[ j ][ k ] *.80;										// 80% saturation of relative water vapour 
		}
	}
*/


}





void BC_Thermo::BC_CO2 ( int i_max, int i_beg, int Ma, double co2_0, double co2_Average, double co2_equator, double co2_pole, double co2_tropopause, double co2_vegetation, double co2_ocean, double co2_land, double t_cretaceous_max, double t_average, Array_2D &co2_j, Array_2D &Vegetation, Array &h, Array &t, Array &p, Array &co2 )
{
// initial and boundary conditions of CO2 content on water and land surfaces
// parabolic CO2 content distribution from pole to pole accepted

	j_half = j_max / 2;


// CO2-distribution by Ruddiman approximated by a parabola
	co2_Cretaceous = 3.2886 * pow ( ( t_Cretaceous + t_average ), 2 ) - 32.8859 * ( t_Cretaceous + t_average ) + 102.2148;  // in ppm
	co2_Cretaceous = co2_Cretaceous - co2_Average;
	if ( Ma == 0 ) 	co2_Cretaceous = 0.;


	co2_comment = "      co2 increase at cretaceous times: ";
	co2_gain = " co2 cretaceous";
	co2_modern = "      mean co2 at modern times: ";
	co2_average = " co2 modern";
	co2_unit = "ppm ";

	cout << endl << setiosflags ( ios::left ) << setw ( 48 ) << setfill ( '.' ) << co2_comment << resetiosflags ( ios::left ) << setw ( 12 ) << co2_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << co2_Cretaceous << setw ( 5 ) << co2_unit << endl << setw ( 48 ) << setfill ( '.' )  << setiosflags ( ios::left ) << co2_modern << resetiosflags ( ios::left ) << setw ( 15 ) << co2_average << " = "  << setw ( 7 ) << setfill ( ' ' ) << co2_Average << setw ( 5 ) << co2_unit << endl << endl;


	co2_coeff = co2_pole - co2_equator;


// CO2-content as initial solution
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( h.x[ 0 ][ j ][ k ]  == 0. ) 
			{
				d_j = ( double ) j;
				co2_j.y[ j ][ k ] = co2_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + co2_pole;
				co2.x[ 0 ][ j ][ k ] = co2_j.y[ j ][ k ] + co2_ocean;															// 0.6/600Gt CO2-source on oceans per year 
			}
			if ( h.x[ 0 ][ j ][ k ]  == 1. ) 
			{
				d_j = ( double ) j;
				co2_j.y[ j ][ k ] = co2_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + co2_pole; // parabolic distribution from pole to pole
				co2_j.y[ j ][ k ] = co2_j.y[ j ][ k ] - co2_vegetation * Vegetation.y[ j ][ k ];						// 100/600Gt CO2-sink by vegetation on land per year
				co2.x[ 0 ][ j ][ k ] = co2_j.y[ j ][ k ] + co2_land;																// 0.2/600Gt CO2-source on land per year everywhere
			}
		}
	}




// co2 distribution decreasing approaching tropopause, above no co2
	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( i <= im_tropopause[ j ] )
				{
					d_i_max = ( double ) im_tropopause[ j ];
					d_i = ( double ) i;
					co2.x[ i ][ j ][ k ] = co2_j.y[ j ][ k ] - ( co2_tropopause - co2_j.y[ j ][ k ] ) * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) );	// radial distribution approximated by a parabola
//					co2.x[ i ][ j ][ k ] = ( co2_tropopause - co2.x[ 0 ][ j ][ k ] ) / d_i_max * d_i + co2.x[ 0 ][ j ][ k ];									// linear co2 decay up to tropopause
				}
				else 			co2.x[ i ][ j ][ k ] = co2_tropopause;
				if ( co2.x[ i ][ j ][ k ] < 0. ) co2.x[ i ][ j ][ k ] = 0.;
			}
		}
	}

}




void BC_Thermo::BC_RadiationBalance_parab ( double sigma, double epsilon, double radiation_ocean, double radiation_pole, double radiation_equator, Array_2D &Radiation_Balance_par, Array &h )
{
// initial and boundary conditions for the radiation balance on water and land surfaces
// parabolic distribution from pole to pole accepted

	j_half = j_max / 2;

	rad_bal_minus = epsilon * sigma * pow ( t_Cretaceous * t_0, 4. );	// decrease of radiation during paleo climate

	radiation_land_coeff = radiation_pole - radiation_equator;
	radiation_ocean_coeff = radiation_pole - radiation_equator - radiation_ocean; // increase by 40 W/m² on ocean


// radiation-distribution as initial solution
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( h.x[ 0 ][ j ][ k ]  == 0. ) 
			{
				d_j = ( double ) j;
				Radiation_Balance_par.y[ j ][ k ] = radiation_ocean_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + radiation_pole;
				Radiation_Balance_par.y[ j ][ k ] = Radiation_Balance_par.y[ j ][ k ] - rad_bal_minus;
			}
			if ( h.x[ 0 ][ j ][ k ]  == 1. ) 
			{
				d_j = ( double ) j;
				Radiation_Balance_par.y[ j ][ k ] = radiation_land_coeff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + radiation_pole; // parabolic distribution from pole to pole
				Radiation_Balance_par.y[ j ][ k ] = Radiation_Balance_par.y[ j ][ k ] - rad_bal_minus;
			}
		}
	}

}






void BC_Thermo::IC_CellStructure ( int i_max, int i_beg, double u_0, Array &u, Array &v, Array &w )
{
// boundary condition for the velocity components in the circulation cells

// latest version by Grotjahn ( Global Atmospheric Circulations, 1993 )
// default for the velocity components u, v, and w as initial conditions


// velocity assumptions at the equator 0°

	ua_00 = 1.;

	va_equator_SL =  0.000;
	va_equator_Tropopause = 0.000;

	wa_equator_SL = - 5.;
	wa_equator_Tropopause = - 7.5;

// velocity assumptions for latitude at 15° and 30° in the Hadley cell

	ua_30 = - 1.;																	// do not change the velocity initial conditions !!

//	va_Hadley_SL = 2.;
	va_Hadley_SL = 1.;
//	va_Hadley_Tropopause = - 2.;
	va_Hadley_Tropopause = - 1.;

//	va_Hadley_SL_15 = 5.;
	va_Hadley_SL_15 = 1.;
//	va_Hadley_Tropopause_15 = - 5.;
	va_Hadley_Tropopause_15 = - 1.;

	wa_Hadley_SL = 5.;									// at surface
	wa_Hadley_Tropopause = 30.;					// subtropic jet

// velocity assumptions for latitude at 45° and 60° in the Ferrel cell

	ua_60 = 0.75;

	va_Ferrel_SL = - 0.75;
//	va_Ferrel_Tropopause = 3.;
	va_Ferrel_Tropopause = 1.;

//	va_Ferrel_SL_45 = - 1.75;
	va_Ferrel_SL_45 = - 1.;
//	va_Ferrel_Tropopause_45 = 3.;
	va_Ferrel_Tropopause_45 = 1.;

	wa_Ferrel_SL = 1.5;								// subpolar jet
	wa_Ferrel_Tropopause = 10.;				// subpolar jet

// velocity assumptions for latitude 90° in the Polar cell

//	ua_90 = - 1.;
	ua_90 = - 0.75;

	va_Polar_SL = 0.75;
	va_Polar_Tropopause = - 0.75;

	va_Polar_SL_75 = 0.75;
	va_Polar_Tropopause_75 = - 0.75;

	va_Polar_SL_75 = 0.75;
	va_Polar_Tropopause_75 = - 0.75;

	wa_Polar_SL = - 0.75;
	wa_Polar_Tropopause = 0.75;





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

			for ( int i = 0; i <= i_half; i++ )
			{
				m = i_max - i;
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = d_i / d_i_half * ua_00;
				u.x[ m ][ j ][ k ] = d_i / d_i_half * ua_00;
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

			for ( int i = 0; i <= i_half; i++ )
			{
				m = i_max - i;
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = d_i / d_i_half * ua_90;
				u.x[ m ][ j ][ k ] = d_i / d_i_half * ua_90;
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

			for ( int i = 0; i <= i_half; i++ )
			{
				m = i_max - i;
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = d_i / d_i_half * ua_60;
				u.x[ m ][ j ][ k ] = d_i / d_i_half * ua_60;
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

			for ( int i = 0; i <= i_half; i++ )
			{
				m = i_max - i;
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = d_i / d_i_half * ua_30;
				u.x[ m ][ j ][ k ] = d_i / d_i_half * ua_30;
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

			for ( int i = 0; i <= i_half; i++ )
			{
				m = i_max - i;
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = d_i / d_i_half * ua_30;
				u.x[ m ][ j ][ k ] = d_i / d_i_half * ua_30;
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

			for ( int i = 0; i <= i_half; i++ )
			{
				m = i_max - i;
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = d_i / d_i_half * ua_60;
				u.x[ m ][ j ][ k ] = d_i / d_i_half * ua_60;
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

			for ( int i = 0; i <= i_half; i++ )
			{
				m = i_max - i;
				d_i = ( double ) i;
				u.x[ i ][ j ][ k ] = d_i / d_i_half * ua_90;
				u.x[ m ][ j ][ k ] = d_i / d_i_half * ua_90;
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

	j_pac_had_n_end = j_had_n - 8;

	k_pac_w = 120;
	k_pac_w_end = k_pac_w - 5;
	k_pac_e = 240;

	d_j_w = ( double ) j_aeq;

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			v.x[ i ][ j ][ k_pac_w ] = - .5 * v.x[ i ][ j ][ k_pac_w ];

//			d_j = ( double ) j;																// asymmetric extention of highs
//			v.x[ i ][ j ][ k_pac_w ] = .5 * ( v.x[ i ][ j_had_n ][ k_pac_w ] - v.x[ i ][ j_had_s ][ k_pac_w ] ) * ( d_j * d_j / ( d_j_w * d_j_w ) - 2. * d_j / d_j_w );

		}
	}

/////////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_pac_w; k <= k_pac_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_e ] - v.x[ i ][ j ][ k_pac_w ] ) / ( ( int ) ( k_pac_e - k_pac_w ) ) * ( int ) ( k - k_pac_w ) + v.x[ i ][ j ][ k_pac_w ];
			}
		}
	}

///////////////////////////////////////////// smoothing transitions /////////////////////////////////////////////////

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_pac_w_end+1; k <= k_pac_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_w ] - v.x[ i ][ j ][ k_pac_w_end ] ) / ( ( int ) ( k_pac_w - k_pac_w_end ) ) * ( int ) ( k - k_pac_w_end ) + v.x[ i ][ j ][ k_pac_w_end ];
			}
		}
	}



///////////////////////////////////////////////// Pacific ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// South Pacific
///////////////////////////////////////////////// change in sign of v-component in the southern hemisphere /////////////////////////////////////////////////


	j_pac_had_s_end = j_had_s + 8;

	k_pac_w = 130;
	k_pac_w_end = k_pac_w - 5;
	k_pac_e = 260;

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			v.x[ i ][ j ][ k_pac_w ] = - .5 * v.x[ i ][ j ][ k_pac_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			for ( int k = k_pac_w; k <= k_pac_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_e ] - v.x[ i ][ j ][ k_pac_w ] ) / ( ( int ) ( k_pac_e - k_pac_w ) ) * ( int ) ( k - k_pac_w ) + v.x[ i ][ j ][ k_pac_w ];
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			for ( int k = k_pac_w_end+1; k <= k_pac_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_w ] - v.x[ i ][ j ][ k_pac_w_end ] ) / ( ( int ) ( k_pac_w - k_pac_w_end ) ) * ( int ) ( k - k_pac_w_end ) + v.x[ i ][ j ][ k_pac_w_end ];
			}
		}
	}





///////////////////////////////////////////////// Indic ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// North Indic .......................................................................... only on land
///////////////////////////////////////////////// change in sign of v-component in the northern hemisphere //////////////////////////

	j_pac_had_n_end = j_had_n - 5;
	j_pac_had_s_end = j_had_s + 5;

	k_pac_w = 30;
	k_pac_w_end = k_pac_w - 10;
	k_pac_e = 90;

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			v.x[ i ][ j ][ k_pac_w ] = - .5 * v.x[ i ][ j ][ k_pac_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_pac_w; k <= k_pac_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_e ] - v.x[ i ][ j ][ k_pac_w ] ) / ( ( int ) ( k_pac_e - k_pac_w ) ) * ( int ) ( k - k_pac_w ) + v.x[ i ][ j ][ k_pac_w ];
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_pac_w_end+1; k <= k_pac_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_w ] - v.x[ i ][ j ][ k_pac_w_end ] ) / ( ( int ) ( k_pac_w - k_pac_w_end ) ) * ( int ) ( k - k_pac_w_end ) + v.x[ i ][ j ][ k_pac_w_end ];
			}
		}
	}



///////////////////////////////////////////////// Indicic ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// South Indic
///////////////////////////////////////////////// change in sign of v-component in the southern hemisphere /////////////////////////////////////////////////

	j_pac_had_n_end = j_had_n - 5;
	j_pac_had_s_end = j_had_s + 5;

	k_pac_w = 35;
	k_pac_w_end = k_pac_w - 3;
	k_pac_e = 90;

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			v.x[ i ][ j ][ k_pac_w ] = - .5 * v.x[ i ][ j ][ k_pac_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			for ( int k = k_pac_w; k <= k_pac_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_e ] - v.x[ i ][ j ][ k_pac_w ] ) / ( ( int ) ( k_pac_e - k_pac_w ) ) * ( int ) ( k - k_pac_w ) + v.x[ i ][ j ][ k_pac_w ];
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			for ( int k = k_pac_w_end+1; k <= k_pac_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_w ] - v.x[ i ][ j ][ k_pac_w_end ] ) / ( ( int ) ( k_pac_w - k_pac_w_end ) ) * ( int ) ( k - k_pac_w_end ) + v.x[ i ][ j ][ k_pac_w_end ];
			}
		}
	}



///////////////////////////////////////////////// Atlantic ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// North Altlantic
///////////////////////////////////////////////// change in sign of v-component in the northern hemisphere //////////////////////////

	j_pac_had_n_end = j_had_n - 5;
	j_pac_had_s_end = j_had_s + 5;

	k_pac_w = 280;
//	k_pac_w_end = k_pac_w - 10;
//	k_pac_e = 340;
	k_pac_w_end = k_pac_w - 5;
	k_pac_e = 330;

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			v.x[ i ][ j ][ k_pac_w ] = - .5 * v.x[ i ][ j ][ k_pac_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_pac_w; k <= k_pac_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_e ] - v.x[ i ][ j ][ k_pac_w ] ) / ( ( int ) ( k_pac_e - k_pac_w ) ) * ( int ) ( k - k_pac_w ) + v.x[ i ][ j ][ k_pac_w ];
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_pac_had_n_end; j < j_aeq; j++ )
		{
			for ( int k = k_pac_w_end+1; k <= k_pac_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_w ] - v.x[ i ][ j ][ k_pac_w_end ] ) / ( ( int ) ( k_pac_w - k_pac_w_end ) ) * ( int ) ( k - k_pac_w_end ) + v.x[ i ][ j ][ k_pac_w_end ];
			}
		}
	}



///////////////////////////////////////////////// Atlantic ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// South Atlantic
///////////////////////////////////////////////// change in sign of v-component in the southern hemisphere /////////////////////////////////////////////////

	j_pac_had_n_end = j_had_n - 5;
	j_pac_had_s_end = j_had_s + 5;

	k_pac_w = 320;
//	k_pac_w_end = k_pac_w - 10;
//	k_pac_e = 360;
	k_pac_w_end = k_pac_w - 5;
	k_pac_e = 360;

	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			v.x[ i ][ j ][ k_pac_w ] = - .5 * v.x[ i ][ j ][ k_pac_w ];
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			for ( int k = k_pac_w; k <= k_pac_e; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_e ] - v.x[ i ][ j ][ k_pac_w ] ) / ( ( int ) ( k_pac_e - k_pac_w ) ) * ( int ) ( k - k_pac_w ) + v.x[ i ][ j ][ k_pac_w ];
			}
		}
	}


	for ( int i = 0; i <= i_half; i++ )
	{
		for ( int j = j_aeq+1; j <= j_pac_had_s_end; j++ )
		{
			for ( int k = k_pac_w_end+1; k <= k_pac_w; k++ )
			{
				v.x[ i ][ j ][ k ] = ( v.x[ i ][ j ][ k_pac_w ] - v.x[ i ][ j ][ k_pac_w_end ] ) / ( ( int ) ( k_pac_w - k_pac_w_end ) ) * ( int ) ( k - k_pac_w_end ) + v.x[ i ][ j ][ k_pac_w_end ];
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
			}
		}
	}




///////////////////////////////////////////////// end of smoothing transitions from cell to cell //////////////////////////
///////////////////////////////////////////////// end of smoothing transitions from cell to cell //////////////////////////

}










void BC_Thermo::BC_Surface_Temperature ( const string &Name_SurfaceTemperature_File, Array_2D &t_j, Array &t )
{
// initial conditions for the Name_SurfaceTemperature_File at the sea surface

	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 3 );
	cout.setf ( ios::fixed );

// reading data from file Name_SurfaceTemperature_File_Read
	ifstream Name_SurfaceTemperature_File_Read;
	Name_SurfaceTemperature_File_Read.open ( Name_SurfaceTemperature_File.c_str(), ios_base::in );
	Name_SurfaceTemperature_File_Read.seekg ( 0L, ios::beg );
	anfangpos_1 = Name_SurfaceTemperature_File_Read.tellg ();


	if ( Name_SurfaceTemperature_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: could be opened" << endl;
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: begins at ::::::: " << anfangpos_1 << endl;
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

			t.x[ 0 ][ j ][ k ] = t_j.y[ j ][ k ] = ( dummy_3 + 273.15 ) / 273.15;
			j++;
		}
	j = 0;
	k++;
	}

// correction of surface temperature around 180°E
	for ( int j = 0; j < jm; j++ )
	{
		t.x[ 0 ][ j ][ k_half ] = t_j.y[ j ][ k_half ] = ( t_j.y[ j ][ k_half + 1 ] + t_j.y[ j ][ k_half - 1 ] ) / 2.;
	}


// Ende Lesen von Name_SurfaceTemperature_File

	Name_SurfaceTemperature_File_Read.seekg ( 0L, ios::end );
	endpos_1 = Name_SurfaceTemperature_File_Read.tellg ();

// Abschlussanweisungen für den Dateiabschluss (Dateiverwaltung)

	cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: ends at ::::::::: " << endpos_1 << endl;
	cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: has the length of ::::: " << endpos_1 - anfangpos_1 << " Bytes!"<< endl;

// Im Falle eines Lesefehlers

	if ( Name_SurfaceTemperature_File_Read == NULL )
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: not yet exists! ::::::::: " << endl << endl << endl;
	}

	Name_SurfaceTemperature_File_Read.close();

	if ( Name_SurfaceTemperature_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: could be closed." << endl;
		cout << endl;
	}

	if ( Name_SurfaceTemperature_File_Read.fail() )
		cout << "***** file ::::: " << Name_SurfaceTemperature_File << " ::::: could not be closed!" << endl;

// Ende Lesen von Name_SurfaceTemperature_File_Read
}







void BC_Thermo::BC_Surface_Precipitation ( const string &Name_SurfacePrecipitation_File, Array_2D &precipitation_j )
{
// initial conditions for the Name_SurfacePrecipitation_File at the sea surface

	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 3 );
	cout.setf ( ios::fixed );


// reading data from file Name_SurfacePrecipitation_File_Read
	ifstream Name_SurfacePrecipitation_File_Read;
	Name_SurfacePrecipitation_File_Read.open ( Name_SurfacePrecipitation_File.c_str(), ios_base::in );
	Name_SurfacePrecipitation_File_Read.seekg ( 0L, ios::beg );
	anfangpos_1 = Name_SurfacePrecipitation_File_Read.tellg ();


	if ( Name_SurfacePrecipitation_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: could be opened" << endl;
		cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: begins at ::::::: " << anfangpos_1 << endl;
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

			precipitation_j.y[ j ][ k ] = dummy_3;
			j++;
		}
	j = 0;
	k++;
	}



// Ende Lesen von Name_SurfacePrecipitation_File

	Name_SurfacePrecipitation_File_Read.seekg ( 0L, ios::end );
	endpos_1 = Name_SurfacePrecipitation_File_Read.tellg ();

// Abschlussanweisungen für den Dateiabschluss (Dateiverwaltung)

	cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: ends at ::::::::: " << endpos_1 << endl;
	cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: has the length of ::::: " << endpos_1 - anfangpos_1 << " Bytes!"<< endl;

// Im Falle eines Lesefehlers

	if ( Name_SurfacePrecipitation_File_Read == NULL )
	{
		cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: not yet exists! ::::::::: " << endl << endl << endl;
	}

	Name_SurfacePrecipitation_File_Read.close();

	if ( Name_SurfacePrecipitation_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: could be closed." << endl;
		cout << endl;
	}

	if ( Name_SurfacePrecipitation_File_Read.fail() )
		cout << "***** file ::::: " << Name_SurfacePrecipitation_File << " ::::: could not be closed!" << endl;

// Ende Lesen von Name_SurfacePrecipitation_File_Read
}





void BC_Thermo::BC_Pressure ( double const r_0_air, double const R_Air, const double p_0, const double t_0, Array_2D &p_j, Array_2D &t_j, Array &p, Array &t, Array &h )
{
// boundary condition of surface pressure given by surface temperature through gas equation

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			p_j.y[ j ][ k ]  =  ( r_0_air * R_Air * t_j.y[ j ][ k ] * ( t_0 + 20. ) ) / p_0 / 100.;		// given in hPa / hPa
			p.x[ 0 ][ j ][ k ] = p_j.y[ j ][ k ];
		}
	}

// only switch on, when you know what you do

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 1; i < im; i++ )
			{
				d_i = ( double ) i;
				if ( h.x[ i ][ j ][ k ] == 0. )   p.x[ i ][ j ][ k ] = exp ( - 9.8066 * (  double ) i * 500. / ( 287.1 * t.x[ i ][ j ][ k ] * t_0 ) );				// step size in 500 m   given in hPa / hPa
			}
		}
	}

}




void BC_Thermo::IC_v_w_WestEastCoast ( Array &h, Array &u, Array &v, Array &w )
{
// initial conditions for v and w velocity components at the sea surface close to east or west coasts
// reversal of v velocity component between north and south equatorial current ommitted at respectively 10°
// w component unchanged

// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: east coast

	d_i_max = ( double ) i_max;

	k_grad = 10;																			// extension of velocity change
	k_a = k_grad;																			// left distance
	k_b = 0;																					// right distance

	k_water = 0;																			// on water closest to coast
	k_sequel = 1;																			// on solid ground

	for ( int j = 0; j < 91; j++ )														// outer loop: latitude
	{
		for ( int k = 0; k < km; k++ )												// inner loop: longitude
		{
			if ( h.x[ 0 ][ j ][ k ] == 1. ) k_sequel = 0;							// if solid ground: k_sequel = 0

			if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( k_sequel == 0 ) ) k_water = 0;// if water and and k_sequel = 0 then is water closest to coast
			else k_water = 1;															// somewhere on water

			if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( k_water == 0 ) )			// if water is closest to coast, change of velocity components begins
			{
				for ( int l = 0; l < k_grad; l++ )									// extension of change, sign change in v-velocity and distribution of u-velocity with depth
				{
					v.x[ 0 ][ j ][ k + l ] = - v.x[ 0 ][ j ][ k + l ];					// existing velocity changes sign
				}


				for ( int l = ( k + k_grad - k_a ); l < ( k + k_grad + k_b + 1 ); l++ ) // starting at local longitude + max extension - begin of smoothing k_a  until ending at  + k_b
				{
					v.x[ 0 ][ j ][ l ] = ( v.x[ 0 ][ j ][ k + k_grad + k_b ] - v.x[ 0 ][ j ][ k + k_grad - k_a ] ) / ( double )( ( k + k_grad + k_b ) -  ( k + k_grad - k_a ) ) * ( double )( l -  ( k + k_grad - k_a ) ) + v.x[ 0 ][ j ][ k + k_grad - k_a ];	 // extension of v-velocity, smoothing algorithm by a linear equation 

					for ( int i = 1; i <= 11; i++ )										// loop in radial direction, extension for u -velocity component, downwelling here
					{
						d_i = ( double ) i;
						v.x[ i ][ j ][ k ] = - ( 0. - v.x[ 0 ][ j ][ l ] ) / d_i_max * d_i - v.x[ 0 ][ j ][ l ];	// radial distribution approximated by a straight line
					}

				}

				for ( int l = k; l < ( k + k_grad + k_b + 1 ); l++ )			// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
				{
					w.x[ 0 ][ j ][ l ] = w.x[ 0 ][ j ][ k + k_grad + k_b ]  / ( double )( ( k + k_grad + k_b ) -  k ) * ( double )( l - k ); // extension of v-velocity

					for ( int i = 1; i <= 11; i++ )										// loop in radial direction, extension for u -velocity component, downwelling here
					{
						d_i = ( double ) i;
						w.x[ i ][ j ][ k ] = - ( 0. - w.x[ 0 ][ j ][ l ] ) / d_i_max * d_i - w.x[ 0 ][ j ][ l ];					// radial distribution approximated by a straight line
					}

				}

				k_sequel = 1;																// looking for another east coast
			}
		}																							// end of longitudinal loop
		k_water = 0;																		// starting at another latitude
	}																								// end of latitudinal loop




// southern hemisphere: east coast

	k_water = 0;
	k_sequel = 1;

	for ( int j = 91; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ 0 ][ j ][ k ] == 1. ) k_sequel = 0;

			if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( k_sequel == 0 ) ) k_water = 0;
			else k_water = 1;

			if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( k_water == 0 ) )
			{
				for ( int l = 0; l < k_grad; l++ )
				{
					v.x[ 0 ][ j ][ k + l ] = - v.x[ 0 ][ j ][ k + l ];
				}

				for ( int l = ( k + k_grad - k_a ); l < ( k + k_grad + k_b + 1 ); l++ )
				{
					v.x[ 0 ][ j ][ l ] = ( v.x[ 0 ][ j ][ k + k_grad + k_b ] - v.x[ 0 ][ j ][ k + k_grad - k_a ] ) / ( double )( ( k + k_grad + k_b ) -  ( k + k_grad - k_a ) ) * ( double )( l -  ( k + k_grad - k_a ) ) + v.x[ 0 ][ j ][ k + k_grad - k_a ];
/*
					for ( int i = 1; i <= 11; i++ )										// loop in radial direction, extension for u -velocity component, downwelling here
					{
						d_i = ( double ) i;
						v.x[ i ][ j ][ k ] = - ( 0. - v.x[ 0 ][ j ][ l ] ) / d_i_max * d_i - v.x[ 0 ][ j ][ l ];					// radial distribution approximated by a straight line
					}
*/
				}

				for ( int l = k; l < ( k + k_grad + k_b + 1 ); l++ )
				{
					w.x[ 0 ][ j ][ l ] = w.x[ 0 ][ j ][ k + k_grad + k_b ]  / ( double )( ( k + k_grad + k_b ) -  k ) * ( double )( l - k );
/*
					for ( int i = 1; i <= 11; i++ )					// loop in radial direction, extension for u -velocity component, downwelling here
					{
						d_i = ( double ) i;
						w.x[ i ][ j ][ k ] = - ( 0. - w.x[ 0 ][ j ][ l ] ) / d_i_max * d_i - w.x[ 0 ][ j ][ l ];					// radial distribution approximated by a straight line
					}
*/
				}

				k_sequel = 1;
			}
		}
		k_water = 0;
	}




// search for west coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included

// northern hemisphere: west coast
/*
	k_grad = 6;																				// extension of velocity change
//	k_grad = 20;																			// extension of velocity change
	k_a = 0;																					// left distance

	k_water = 0;																			// somewhere on water
	flip = 0;																					// somewhere on water

	for ( int j = 0; j < 91; j++ )														// outer loop: latitude
	{
		for ( int k = 0; k < km; k++ )												// inner loop: longitude
		{
			if ( h.x[ 0 ][ j ][ k ] == 0. )													// if somewhere on water
			{
				k_water = 0;																// somewhere on water: k_water = 0
				flip = 0;																		// somewhere on water: flip = 0
			}
			else k_water = 1;															// first time on land

			if ( ( flip == 0 ) && ( k_water == 1 ) )								// on water closest to land
			{
				for ( int l = k; l > ( k - k_grad + 1 ); l-- )						// backward extention of velocity change: nothing changes
				{
					w.x[ 0 ][ j ][ l ] = - w.x[ 0 ][ j ][ l ];
				}

				for ( int l = k; l > ( k - k_grad - k_a - 1 ); l-- )				// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
				{
					v.x[ 0 ][ j ][ l ] = v.x[ 0 ][ j ][ k - k_grad - k_a ] / ( double )( ( k - k_grad - k_a ) - k ) * ( double )( l - k ); // extension of v-velocity

					for ( int i = 1; i <= 11; i++ )										// loop in radial direction, extension for u -velocity component, downwelling here
					{
						d_i = ( double ) i;
						v.x[ i ][ j ][ k ] = - ( 0. - v.x[ 0 ][ j ][ l ] ) / d_i_max * d_i - v.x[ 0 ][ j ][ l ];					// radial distribution approximated by a straight line
					}

				}

				for ( int l = ( k - k_grad - 3 ); l < ( k - k_grad + 3 ); l++ )			// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
				{
					w.x[ 0 ][ j ][ l ] = (  - w.x[ 0 ][ j ][ k - k_grad - 3 ] + w.x[ 0 ][ j ][ k - k_grad + 3 ] ) * ( double ) ( l - ( k - k_grad - 3 ) ) / ( double ) ( ( k - k_grad + 3 ) - ( k - k_grad - 3 ) ) - w.x[ 0 ][ j ][ k - k_grad + 3 ];

					for ( int i = 1; i <= 11; i++ )										// loop in radial direction, extension for u -velocity component, downwelling here
					{
						d_i = ( double ) i;
						w.x[ i ][ j ][ k ] = - ( 0. - w.x[ 0 ][ j ][ l ] ) / d_i_max * d_i - w.x[ 0 ][ j ][ l ];					// radial distribution approximated by a straight line
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

	for ( int j = 91; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ 0 ][ j ][ k ] == 0. )
			{
				k_water = 0;
				flip = 0;
			}
			else k_water = 1;

			if ( ( flip == 0 ) && ( k_water == 1 ) )
			{
				for ( int l = k; l > ( k - k_grad + 1 ); l-- )
				{
					w.x[ 0 ][ j ][ l ] = - w.x[ 0 ][ j ][ l ];
				}

				for ( int l = k; l > ( k - k_grad - k_a - 1 ); l-- )
				{
					v.x[ 0 ][ j ][ l ] = v.x[ 0 ][ j ][ k - k_grad - k_a ] / ( double )( ( k - k_grad - k_a ) - k ) * ( double )( l - k );

					for ( int i = 1; i <= 11; i++ )										// loop in radial direction, extension for u -velocity component, downwelling here
					{
						d_i = ( double ) i;
						v.x[ i ][ j ][ k ] = - ( 0. - v.x[ 0 ][ j ][ l ] ) / d_i_max * d_i - v.x[ 0 ][ j ][ l ];					// radial distribution approximated by a straight line
					}

				}

				for ( int l = ( k - k_grad - 3 ); l < ( k - k_grad + 3 ); l++ )			// smoothing algorithm by a linear equation, starting at local longitude until ending at max extension + k_b
				{
					w.x[ 0 ][ j ][ l ] = ( - w.x[ 0 ][ j ][ k - k_grad - 3 ] + w.x[ 0 ][ j ][ k - k_grad + 3 ] ) * ( double ) ( l - ( k - k_grad - 3 ) ) / ( double ) ( ( k - k_grad + 3 ) - ( k - k_grad - 3 ) ) - w.x[ 0 ][ j ][ k - k_grad + 3 ];

					for ( int i = 1; i <= 11; i++ )					// loop in radial direction, extension for u -velocity component, downwelling here
					{
						d_i = ( double ) i;
						w.x[ i ][ j ][ k ] = - ( 0. - w.x[ 0 ][ j ][ l ] ) / d_i_max * d_i - w.x[ 0 ][ j ][ l ];					// radial distribution approximated by a straight line
					}

				}

				flip = 1;
			}
		}
		flip = 0;
	}
*/
}







void BC_Thermo::IC_v_w_Smoothing ( int velocity_iter, Array &h, Array &u, Array &v, Array &w, Array &t, Array &c )
{
// initial conditions for v and w velocity components at the sea surface
// after reading wind data, applying Ekman motion and formulation of west/east coast corretions 

	jmkm = ( double ) ( ( jm -1 ) * ( km -1 ) );

	if ( velocity_iter == 0 )
	{
		v_sum = 0.;
		w_sum = 0.;
		t_sum = 0.;
		c_sum = 0.;
		n_smooth = 0;

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v_sum += fabs ( v.x[ im-1 ][ j ][ k ] );
					w_sum += fabs ( w.x[ im-1 ][ j ][ k ] );
					n_smooth++;
				}
				else
				{
					v.x[ im-1 ][ j ][ k ] = 0.;
					w.x[ im-1 ][ j ][ k ] = 0.;
				}
			}
		}


		v_sum = v_sum / ( double ) n_smooth;
		w_sum = w_sum / ( double ) n_smooth;
		t_sum = t_sum / ( double ) n_smooth;
		c_sum = c_sum / ( double ) n_smooth;


		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					if ( fabs ( v.x[ im-1 ][ j ][ k ] ) > v_sum )			v.x[ im-1 ][ j ][ k ] = 0.;
					if ( fabs ( w.x[ im-1 ][ j ][ k ] ) > w_sum )			w.x[ im-1 ][ j ][ k ] = 0.;
				}
				else
				{
					v.x[ im-1 ][ j ][ k ] = 0.;
					w.x[ im-1 ][ j ][ k ] = 0.;
				}
			}
		}
	}



/*
	v_sum = 0.;
	w_sum = 0.;
	t_sum = 0.;
	c_sum = 0.;
	n_smooth = 0;


	for ( int i = 1; i < im-2; i++ )
	{
		for ( int j = 1; j < jm-2; j++ )
		{
			for ( int k = 1; k < km-2; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					v_sum += fabs ( v.x[ i ][ j ][ k ] );
					w_sum += fabs ( w.x[ i ][ j ][ k ] );
					t_sum += t.x[ i ][ j ][ k ];
					c_sum += c.x[ i ][ j ][ k ];
					n_smooth++;
				}
			}
		}
	}


	u_sum = u_sum / ( double ) n_smooth;
	v_sum = v_sum / ( double ) n_smooth;
	w_sum = w_sum / ( double ) n_smooth;
	t_sum = t_sum / ( double ) n_smooth;
	c_sum = c_sum / ( double ) n_smooth;


	for ( int i = 1; i < im-2; i++ )
	{
		for ( int j = 1; j < jm-2; j++ )
		{
			for ( int k = 1; k < km-2; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 0. )
				{
					if ( fabs ( u.x[ i ][ j ][ k ] ) > u_sum )			u.x[ i ][ j ][ k ] = ( u.x[ i+1 ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] + u.x[ i ][ j+1 ][ k ] + u.x[ i ][ j-1 ][ k ] + u.x[ i ][ j ][ k+1 ] + u.x[ i ][ j ][ k-1 ] ) / 6.;
					if ( fabs ( v.x[ i ][ j ][ k ] ) > v_sum )			v.x[ i ][ j ][ k ] = ( v.x[ i+1 ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] + v.x[ i ][ j+1 ][ k ] + v.x[ i ][ j-1 ][ k ] + v.x[ i ][ j ][ k+1 ] + v.x[ i ][ j ][ k-1 ] ) / 6.;
					if ( fabs ( w.x[ i ][ j ][ k ] ) > w_sum )		w.x[ i ][ j ][ k ] = ( w.x[ i+1 ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] + w.x[ i ][ j+1 ][ k ] + w.x[ i ][ j-1 ][ k ] + w.x[ i ][ j ][ k+1 ] + w.x[ i ][ j ][ k-1 ] ) / 6.;
					if ( t.x[ i ][ j ][ k ] > t_sum )						t.x[ i ][ j ][ k ] = ( t.x[ i+1 ][ j ][ k ] + t.x[ i-1 ][ j ][ k ] + t.x[ i ][ j+1 ][ k ] + t.x[ i ][ j-1 ][ k ] + t.x[ i ][ j ][ k+1 ] + t.x[ i ][ j ][ k-1 ] ) / 6.;
					if ( c.x[ i ][ j ][ k ] > c_sum )						c.x[ i ][ j ][ k ] = ( c.x[ i+1 ][ j ][ k ] + c.x[ i-1 ][ j ][ k ] + c.x[ i ][ j+1 ][ k ] + c.x[ i ][ j-1 ][ k ] + c.x[ i ][ j ][ k+1 ] + c.x[ i ][ j ][ k-1 ] ) / 6.;
				}
			}
		}
	}
*/
}

