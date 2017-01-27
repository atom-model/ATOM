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


Results_MSL_Atm::Results_MSL_Atm ( int im, int jm, int km, double g, double ep, double hp, double u_0, double p_0, double t_0, double c_0, double co2_0, double sigma, double albedo_extra, double lv, double ls, double cp_l, double L_atm, double dt, double dr, double dthe, double dphi, double r_air, double R_Air, double r_water, double r_water_vapour, double R_WaterVapour, double co2_vegetation, double co2_ocean, double co2_land, double gam, double t_pole, double t_cretaceous )
:	coeff_Diffusion_latent ( .005 ),													// diffusion coefficient for latent heat in [m²/s]    0.015
	coeff_Diffusion_sensibel ( .2 ),													// diffusion coefficient for sensible heat in [m²/s]     0.03
	f_Haude ( .3 )																			// Haude factor for evapotranspiration 0.3 for low, dense vegetation as raw average value by Kuttler
{
	coeff_mmWS = r_air / r_water_vapour;									// coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]
	coeff_lv = lv / ( cp_l * t_0 );														// coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_lv = 9.1069 in [ / ]
	coeff_ls = ls / ( cp_l * t_0 );														// coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_ls = 10.3091 in [ / ]

	c43 = 4./3.;
	c13 = 1./3.;

	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> L_atm = L_atm;
	this-> dt = dt;
	this-> dr = dr;
	this-> dthe = dthe;
	this-> dphi = dphi;
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



// array "e_d" for the cumulus convection
	e_d = 0L;

	e_d = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		e_d[ l ] = 0.;
	}

// array "e_l" for the cumulus convection
	e_l = 0L;

	e_l = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		e_l[ l ] = 0.;
	}

// array "e_p" for the cumulus convection
	e_p = 0L;

	e_p = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		e_p[ l ] = 0.;
	}

// array "g_p" for the cumulus convection
	g_p = 0L;

	g_p = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		g_p[ l ] = 0.;
	}

// array "c_u" for the cumulus convection
	c_u = 0L;

	c_u = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		c_u[ l ] = 0.;
	}

// array "r_humid_u" for the cumulus convection
	r_humid_u = 0L;

	r_humid_u = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		r_humid_u[ l ] = 0.;
	}

// array "r_humid_u_parc" for the cumulus convection
	r_humid_u_parc = 0L;

	r_humid_u_parc = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		r_humid_u_parc[ l ] = 0.;
	}

// array "r_dry_u" for the cumulus convection
	r_dry_u = 0L;


	r_dry_u = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		r_dry_u[ l ] = 0.;
	}

// array "r_dry_u_parc" for the cumulus convection
	r_dry_u_parc = 0L;

	r_dry_u_parc = new double[ im ];

	for ( int l = 0; l < im; l++ )
	{
		r_dry_u_parc[ l ] = 0.;
	}

}


Results_MSL_Atm::~Results_MSL_Atm ()
{
	delete [  ] e_d;
	delete [  ] e_l;
	delete [  ] e_p;
	delete [  ] g_p;
	delete [  ] c_u;
	delete [  ] r_humid_u;
	delete [  ] r_humid_u_parc;
	delete [  ] r_dry_u;
	delete [  ] r_dry_u_parc;
}



void Results_MSL_Atm::run_MSL_data ( int n, int velocity_iter_max, int RadiationModel, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &c, Array &cn, Array &co2, Array &co2n, Array &t, Array &tn, Array &p_dyn, Array &p_stat, Array &BuoyancyForce, Array &u, Array &v, Array &w, Array &Latency, Array &Q_Sensible, Array &radiation_3D, Array &t_cond_3D, Array &t_evap_3D, Array &cloud, Array &cloudn, Array &ice, Array &icen, Array &P_rain, Array &P_snow, Array &aux_u, Array &aux_v, Array &aux_w, Array_2D &precipitation_NASA, Array_2D &Evaporation, Array_2D &Condensation, Array_2D &LatentHeat, Array_2D &precipitable_water, Array_2D &Q_Radiation, Array_2D &Q_Evaporation, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Q_bottom, Array_2D &Evaporation_Penman, Array_2D &Evaporation_Haude, Array_2D &Vegetation, Array_2D &Radiation_Balance, Array_2D &Radiation_Balance_par, Array_2D &Radiation_Balance_bot, Array_2D &albedo, Array_2D &co2_total, Array_2D &Precipitation, Array &S_v, Array &S_c, Array &S_i, Array &S_r, Array &S_s )
{
// determination of temperature and pressure by the law of Clausius-Clapeyron for water vapour concentration
// reaching saturation of water vapour pressure leads to formation of rain or ice
// precipitation and cloud formation by formulas from Häckel
// dry adiabatic lapse rate and saturated adiabatic lapse rate = temperature decrease with hight
// SL stands for sea level

// constant coefficients for the transport of cloud water and cloud ice amount vice versa, rain and snow in the parameterization procedures
	a_if = .66;
	c_ac = .24;
	c_rim = 18.6;
	bet_ev = 5.9;
	alf_melt = 7.2e-6;
	bet_melt = bet_dep = 13.;
	alf_if = 1.92e-6;
	alf_cf = 1.55e-3;
	E_cf = 5.0e-3;
	tau_r = 1.e4;
	tau_s = 1.e3;
	t_1 = 253.15;
	t_2 = 235.15;
	a_mc = .08;
	a_mv = .02;
	t_m1 = .5 * ( t_0 + t_1 );
	t_m2 = .5 * ( t_0 + t_2 );
	N_cf_0_surf = 2.e5;
	N_cf_0_6km = 1.e4;
	N_i_0 = 1.e2;																			// in m3
	t_nuc = 267.15;																		// in K
	t_d = 248.15;																		// in K
	t_hn = 236.15;																		// in K
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
	c_s_melt = 8.43e-5;
	b_s_melt = 12.05;
	a_s_melt = 2.31e3;
	a_i_m = 130.;																			// in kg/m3
	a_s_m = .038;																			// in kg/m3
	N_r_0 = 8.e6;																			// in 1/m4
	N_s_0 = 8.e5;																			// in 1/m4
	b_u = .3;
	alf_1 = 5.e-4;
	alf_2 = .011;
	p_ps = .05;
	bet_p = 2.e-3;																			// in s
	t_Celsius_1 = t_2 - t_0;
	T_t_in0 = t_0 - t_2;


// setting water vapour, cloud water and cloud ice into the proper thermodynamic ratio based on the local temperatures
// starting from a guessed parabolic temperature and water vapour distribution in north/south direction
//	for ( int k = 0; k < km; k++ )
//	{
//		for ( int j = 0; j < jm; j++ )
//		{

	for ( int k = 1; k < km-1; k++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			t_u = t.x[ 0 ][ j ][ k ] * t_0;

			t_Celsius_SL = t_u - t_0;
			p_SL = .01 * r_air * R_Air * t_u;																// from gas equation given in hPa

			E_Rain = hp * exp ( 17.0809 * t_Celsius_SL / ( 234.175 + t_Celsius_SL ) );		// saturation water vapour pressure for the water phase at t > 0°C in hPa
			q_Rain  = ep * E_Rain / ( p_SL - E_Rain );														// water vapour amount at saturation with water formation in kg/kg

			E_Ice = hp * exp ( 22.4429 * t_Celsius_SL / ( 272.44 + t_Celsius_SL ) );			// saturation water vapour pressure for the ice phase in hPa
			q_Ice = ep * E_Ice / ( p_SL - E_Ice );															// water vapour amount at saturation with ice formation in kg/kg

			r_dry = 100. * p_SL / ( R_Air * t_u );
			r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ 0 ][ j ][ k ] );		// density of humid air, COSMO version withot cloud and ice water, masses negligible

			e_SL = .01 * ( r_humid * R_WaterVapour * t_u );											// delivers the same results
			a_SL = 216.6 * e_SL / t_u;																			// absolute humidity in kg/m3 at sea level

			t_dew_SL = ( 423.86 - 234.175 * log ( e_SL ) ) / ( log ( e_SL ) - 18.89 );			// dewpoint temperature on ground in °C 		by Häckel
			h_level = 122. * ( ( t_u - t_0 ) - t_dew_SL );													// Condensation level in m		by Häckel + correction
			i_level = ( int ) ( h_level / 500. );																// Condensation level in radial steps

			for ( int i = 1; i < im-1; i++ )
			{
				q_v_b = c.x[ i ][ j ][ k ];
				q_c_b = cloud.x[ i ][ j ][ k ];
				q_i_b = ice.x[ i ][ j ][ k ];

				q_T = q_v_b + q_c_b + q_i_b;																	// total water content in the cumulus
				t_u = t.x[ i ][ j ][ k ] * t_0;																		// in K

				t_Celsius = t_u - t_0;																				// in C

				if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_u ) ) * p_SL;
				else 					p_h = p_SL;

				r_dry = 100. * p_h / ( R_Air * t_u );
				r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );	// density of humid air, COSMO version withot cloud and ice water, masses negligible

				e_h = .01 * r_humid * R_WaterVapour * t_u;											// delivers the same results

				a_h = 216.6 * e_h / t_u;																			// absolute humidity in kg/m3
				q_h = c.x[ i ][ j ][ k ];																				// threshold value for water vapour at local hight h in kg/kg

				E_Rain = hp * exp ( 17.0809 * t_Celsius / ( 234.175 + t_Celsius ) );			// saturation water vapour pressure for the water phase at t > 0°C in hPa
				E_Rain_super = hp * exp ( 17.8436 * t_Celsius / ( 245.425 + t_Celsius ) );	// saturation water vapour pressure for the water phase at t < 0°C, supercooled in hPa
				E_Ice = hp * exp ( 22.4429 * t_Celsius / ( 272.44 + t_Celsius ) );				// saturation water vapour pressure for the ice phase in hPa

				q_Rain = ep * E_Rain / ( p_h - E_Rain );													// water vapour amount at saturation with water formation in kg/kg
				q_Rain_super = ep * E_Rain_super / ( p_h - E_Rain_super );						// water vapour amount at saturation with water formation in kg/kg
				q_Ice = ep * E_Ice / ( p_h - E_Ice );															// water vapour amount at saturation with ice formation in kg/kg

// precipitation and cloud formation from formulas by Häckel
// h stands for the local position

				h_h = - R_Air * t_u / g * log ( p_h / p_0 );													// barometric elevation formula solved for h corresponding to hight over ground in m
				t_dew = ( 423.86 - 234.175 * log ( e_h ) ) / ( log ( e_h ) - 18.89 );				// current dewpoint temperature in °C 		by Häckel
				sat_deficit = E_Rain - e_h;																		// saturation deficit, if positive then saturation is less than 100%

				if ( h.x[ i ][ j ][ k ] == 0. )  RF_e = e_h / E_Rain * 100.;							// relative humidity at any point in %
				else  RF_e = 0.;

				if ( RF_e > 100. )  RF_e = 100.;																// remains at 100% relative humidity

// availability of cloud water and cloud ice
				if ( t_Celsius >= 0. )
				{
					ice.x[ i ][ j ][ k ] = 0.;																		// no cloud ice available
				}

				if ( t_Celsius <= t_Celsius_1 )
				{
					cloud.x[ i ][ j ][ k ] = 0.;																		// no cloud water available
				}




				if ( n == 1 )
				{
// warm cloud phase
					q_T = c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ];											// total water content
					t_Celsius = ( t.x[ i ][ j ][ k ] * t_0 - lv / cp_l * cloud.x[ i ][ j ][ k ] ) - t_0;

					if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_u ) ) * p_SL;
					else 					p_h = p_SL;

					E_Rain = hp * exp ( 17.0809 * t_Celsius / ( 234.175 + t_Celsius ) );		// saturation water vapour pressure for the water phase at t > 0°C in hPa
					q_Rain = ep * E_Rain / ( p_h - E_Rain );												// water vapour amount at saturation with water formation in kg/kg
					q_Rain_n = q_Rain;

					if ( q_T <= q_Rain )
					{
						c.x[ i ][ j ][ k ] = q_T;																		// total water amount as water vapour
						cloud.x[ i ][ j ][ k ] = 0.;																	// no cloud water available
					}
					else
					{
						iter_prec = 0;
						while ( iter_prec <= 20 )																// iter_rad may be varied
						{
							iter_prec = iter_prec + 1;

							t_Celsius = ( t_u + lv / cp_l * c.x[ i ][ j ][ k ] - lv / cp_l * q_Rain ) - t_0;
							E_Rain = hp * exp ( 17.0809 * t_Celsius / ( 234.175 + t_Celsius ) ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
							q_Rain = ep * E_Rain / ( p_h - E_Rain );										// water vapour amount at saturation with water formation in kg/kg
							q_Rain = .5 * ( q_Rain_n + q_Rain );

							c.x[ i ][ j ][ k ] = q_Rain;															// water vapour restricted to saturated water vapour amount
							cloud.x[ i ][ j ][ k ] = q_T - c.x[ i ][ j ][ k ];									// cloud water amount

							if ( cloud.x[ i ][ j ][ k ] >= ice.x[ i ][ j ][ k ] ) 			cloud.x[ i ][ j ][ k ] = q_Rain;
							if ( ice.x[ i ][ j ][ k ] >= cloud.x[ i ][ j ][ k ] ) 			ice.x[ i ][ j ][ k ] = q_Ice;

							if ( c.x[ i ][ j ][ k ] < 0. ) 											c.x[ i ][ j ][ k ] = 0.;
							if ( cloud.x[ i ][ j ][ k ] < 0. ) 									cloud.x[ i ][ j ][ k ] = 0.;

							if ( fabs ( q_Rain / q_Rain_n - 1. ) <= 1.e-5 ) 		break;
							q_Rain_n = q_Rain;
						}
					}

					t.x[ i ][ j ][ k ] = ( t_Celsius + t_0 ) / t_0;
					t_u = t.x[ i ][ j ][ k ] * t_0;																	// in K



// mixed cloud phase
					if ( ( t_Celsius < 0. ) && ( t_Celsius >= t_Celsius_1 ) )
					{
						q_v_hyp = ( q_c_b * q_Rain + q_i_b * q_Ice ) / ( q_c_b + q_i_b );
						q_v_hyp_n = q_v_hyp;

						if ( q_T >= q_Rain )
						{
							cloud.x[ i ][ j ][ k ] = q_Rain;
							ice.x[ i ][ j ][ k ] = q_Ice;
						}

						T = T_nue = t_u;																			// in K
						T_tilda_h = t_u + lv / cp_l * q_v_b;													// in K

						iter_prec = 0;
						while ( iter_prec <= 100 )																// iter_prec may be varied
						{
							iter_prec = iter_prec + 1;

							T_t_in = t_0 - T;
							T_t_in1 = T - t_2;

							CND = T_t_in1 / T_t_in0;
							DEP = T_t_in / T_t_in0;

							d_q_v = q_v_hyp_n - q_v_hyp;
							d_q_c = - d_q_v * CND;
							d_q_i = - d_q_v * DEP;

							d_t = ( lv * d_q_c + ls * d_q_i ) / cp_l;											// in K

							T_nue = T + d_t;																		// in K
							t_Celsius = T_nue - t_0;																// in C

							if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_u ) ) * p_SL;
							else 					p_h = p_SL;

							E_Rain = hp * exp ( 17.0809 * t_Celsius / ( 234.175 + t_Celsius ) );// saturation water vapour pressure for the water phase at t > 0°C in hPa
							q_Rain = ep * E_Rain / ( p_h - E_Rain );										// water vapour amount at saturation with water formation in kg/kg

							E_Ice = hp * exp ( 22.4429 * t_Celsius / ( 272.44 + t_Celsius ) );	// saturation water vapour pressure for the ice phase in hPa
							q_Ice = ep * E_Ice / ( p_h - E_Ice );												// water vapour amount at saturation with ice formation in kg/kg

							q_v_hyp_n = ( q_c_b * q_Rain + q_i_b * q_Ice ) / ( q_c_b + q_i_b );

							q_v_b = c.x[ i ][ j ][ k ] + d_q_v;
							q_c_b = cloud.x[ i ][ j ][ k ] + d_q_c;
							q_i_b = ice.x[ i ][ j ][ k ] + d_q_i;

							if ( cloud.x[ i ][ j ][ k ] >= ice.x[ i ][ j ][ k ] ) 				cloud.x[ i ][ j ][ k ] = q_Rain;
							if ( ice.x[ i ][ j ][ k ] >= cloud.x[ i ][ j ][ k ] ) 				ice.x[ i ][ j ][ k ] = q_Ice;

							if ( c.x[ i ][ j ][ k ] < 0. ) 												c.x[ i ][ j ][ k ] = 0.;
							if ( cloud.x[ i ][ j ][ k ] < 0. ) 										cloud.x[ i ][ j ][ k ] = 0.;
							if ( ice.x[ i ][ j ][ k ] < 0. ) 											ice.x[ i ][ j ][ k ] = 0.;

							if ( ( iter_prec >= 3 ) && ( fabs ( T - T_nue ) / T_nue <= 1.e-5 ) )		break;

							T = T_nue;
						} 																									// iter_prec end

						t.x[ i ][ j ][ k ] = ( t_Celsius + t_0 ) / t_0;
						c.x[ i ][ j ][ k ] = q_v_b;
						cloud.x[ i ][ j ][ k ] = q_c_b;
						ice.x[ i ][ j ][ k ] = q_i_b;
					}																										// ( ( t_Celsius < 0. ) && ( t_Celsius >= t_Celsius_1 ) )



// ice cloud phase
					if ( t_Celsius < t_Celsius_1 )
					{
						q_T = c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k ];			// total water content
						t_Celsius = ( t.x[ i ][ j ][ k ] * t_0 - ls / cp_l * ice.x[ i ][ j ][ k ] ) - t_0;

						if ( i != 0 ) 			p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t.x[ i ][ j ][ k ] * t_0 ) ) * p_SL;
						else 					p_h = p_SL;

						E_Ice = hp * exp ( 22.4429 * t_Celsius / ( 272.44 + t_Celsius ) );		// saturation water vapour pressure for the ice phase in hPa
						q_Ice = ep * E_Ice / ( p_h - E_Ice );													// water vapour amount at saturation with ice formation in kg/kg
						q_Ice_n = q_Ice;

						if ( q_T <= q_Ice )
						{
							c.x[ i ][ j ][ k ] = q_T;																	// total water amount as water vapour
							cloud.x[ i ][ j ][ k ] = 0.;																// no cloud water available
							ice.x[ i ][ j ][ k ] = 0.;																// no ice water available
						}
						else
						{
							iter_prec = 0;
							while ( iter_prec <= 20 )															// iter_rad may be varied
							{
								iter_prec = iter_prec + 1;
								t_Celsius = ( t_u + lv / cp_l * c.x[ i ][ j ][ k ] - ls / cp_l * q_Ice ) - t_0;
								E_Ice = hp * exp ( 17.0809 * t_Celsius / ( 234.175 + t_Celsius ) );	// saturation water vapour pressure for the water phase at t > 0°C in hPa
								q_Ice  = ep * E_Ice / ( p_h - E_Ice );											// water vapour amount at saturation with water formation in kg/kg
								q_Ice = .5 * ( q_Ice_n + q_Ice );

								c.x[ i ][ j ][ k ] = q_Ice;															// water vapour restricted to saturated water vapour amount
								cloud.x[ i ][ j ][ k ] = 0.;															// cloud water amount
								ice.x[ i ][ j ][ k ] = q_T - c.x[ i ][ j ][ k ];									// cloud water amount

								if ( c.x[ i ][ j ][ k ] < 0. ) c.x[ i ][ j ][ k ] = 0.;
								if ( cloud.x[ i ][ j ][ k ] < 0. ) cloud.x[ i ][ j ][ k ] = 0.;
								if ( ice.x[ i ][ j ][ k ] < 0. ) ice.x[ i ][ j ][ k ] = 0.;

								if ( fabs ( q_Ice / q_Ice_n - 1. ) <= 1.e-5 ) 		break;
								q_Ice_n = q_Ice;
							}
						}

						t.x[ i ][ j ][ k ] = ( t_Celsius + t_0 ) / t_0;
						t_u = t.x[ i ][ j ][ k ] * t_0;																// in K
					}																										// end ( t_Celsius < t_Celsius_1 )
				}																											// end n
			}
		}
	}




/*
// printout for various thermodynamical quantities for the preticipation computations along the equator ( j = 90 )
	cout.precision ( 8 );
				if ( ( i == 5 ) && ( j == 90 ) && ( k == 180 ) )
				{
					cout << endl;
					cout << " i = " << i << "   j = " << j << "   k = " << k << "   i_level = " << i_level  << "   h_level (m) = " << h_level << "   h_h (m) = " << h_h << endl << endl;

					cout << " t_h (°C) = " << t_Celsius << "   p_h (hPa) = " << p_h << "   a_h (g/m3) = " << a_h << "   c_h (g/Kg) = " << c.x[ i ][ j ][ k ] * 1000. << "   q_Rain (g/Kg) = " << q_Rain * 1000. << "   e_h (hPa) = " << e_h << "   E_Rain (hPa) = " << E_Rain  << endl << endl;

					cout << " Rain (g/kg) = " << Rain.x[ i ][ j ][ k ] * 1000. << "   Rain_super (g/kg) = " << Rain_super.x[ i ][ j ][ k ] * 1000. << "   Ice (g/kg) = " << Ice.x[ i ][ j ][ k ] * 1000. << "   E_Ice (hPa) = " << E_Ice << "   sat_deficit (hPa) = " << sat_deficit << "   t_dew (°C) = " << t_dew << "   E_Rain_SL (hPa) = " << E_Rain_SL << endl << endl;

					cout << " t_SL (°C) = " << t_Celsius_SL << "   p_SL (hPa) = " << p_SL << "   a_SL (g/m3) = " << a_SL << "   c_SL (g/Kg) = " << c.x[ 0 ][ j ][ k ] * 1000. << "   e_SL (hPa) = " << e_SL << "   q_Ice (g/Kg) = " << q_Ice * 1000. << "   t_dew_SL (°C) = " << t_dew_SL << endl << endl;

					cout << " Evap_Haude (mm/d) = " << Evap_Haude << "   RF_e (%) = " << RF_e  << "   E_Rain_super (hPa) = " << E_Rain_super << "   q_Rain_super (g/Kg) = " << q_Rain_super * 1000. << "   q_h (g/Kg) = " << q_h * 1000. << "   q_SL (g/Kg) = " << q_SL * 1000. << endl << endl;
				}

			}																						// end i
		}																							// end j
	}																								// end k
*/







// rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
// the choosen scheme is a Two Category Ice Scheme
// besides the transport equation for the water vapour exists two equations for the cloud water and the cloud ice transport
// since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation

	if ( n >= 2 )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					P_rain.x[ i ][ j ][ k ] = 0.;
					P_snow.x[ i ][ j ][ k ] = 0.;
				}

				P_rain.x[ im-1 ][ j ][ k ] = 0.;
				P_snow.x[ im-1 ][ j ][ k ] = 0.;

				for ( int i = im-2; i >= 0; i-- )
				{
					p_h = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t.x[ i ][ j ][ k ] * t_0 ) ) * p_SL;	// given in hPa

					r_dry = 100. * p_h / ( R_Air * t.x[ i ][ j ][ k ] * t_0 );
					r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );									// density of humid air, COSMO version withot cloud and ice water, masses negligible

					coeff_P = r_humid * ( L_atm / ( double ) ( im-1 ) );

					t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;

					q_h = c.x[ i ][ j ][ k ];																												// threshold value for water vapour at local hight h in kg/kg

					E_Rain = hp * exp ( 17.0809 * t_Celsius / ( 234.175 + t_Celsius ) );											// saturation water vapour pressure for the water phase at t > 0°C in hPa
					E_Ice = hp * exp ( 22.4429 * t_Celsius / ( 272.44 + t_Celsius ) );												// saturation water vapour pressure for the ice phase in hPa

					q_Rain  = ep * E_Rain / ( p_h - E_Rain );																					// water vapour amount at saturation with water formation in kg/kg
					q_Ice  = ep * E_Ice / ( p_h - E_Ice );																							// water vapour amount at saturation with ice formation in kg/kg


// ice and snow average size
					if ( t.x[ i ][ j ][ k ] * t_0 <= t_0 ) 									N_i = N_i_0 * exp ( .2 * ( t_0 - t.x[ i ][ j ][ k ] * t_0 ) );
					else 																			N_i = N_i_0;

					if ( ( r_humid * ice.x[ i ][ j ][ k ] / N_i <= m_i_max ) && ( ice.x[ i ][ j ][ k ] > 0. ) )			m_i = r_humid * ice.x[ i ][ j ][ k ] / N_i;
					else 																			m_i = m_i_max;
					if ( m_i <= 0. )																m_i = m_i_max;


// nucleation and depositional growth of cloud ice
					if ( ( t.x[ i ][ j ][ k ] * t_0 < t_d ) && ( ice.x[ i ][ j ][ k ] == 0. ) && ( c.x[ i ][ j ][ k ] >= q_Ice ) ) 															S_nuc = m_i_0 / ( r_humid * dt ) * N_i;
					else 																																																S_nuc = 0.;
					if ( ( t_d <= t.x[ i ][ j ][ k ] * t_0 ) && ( t.x[ i ][ j ][ k ] * t_0 <= t_nuc ) && ( ice.x[ i ][ j ][ k ] == 0. ) && ( c.x[ i ][ j ][ k ] >= q_Rain ) ) 	S_nuc = m_i_0 / ( r_humid * dt ) * N_i;
					else 																																																S_nuc = 0.;

					if ( ( t.x[ i ][ j ][ k ] * t_0 < t_hn ) && ( cloud.x[ i ][ j ][ k ] > 0. ) ) 	S_c_frz = cloud.x[ i ][ j ][ k ] / dt;
					else 																							S_c_frz = 0.;

//					S_c_frz = 0.;

					if ( t.x[ i ][ j ][ k ] * t_0 < t_0 )
					{
						S_i_dep = c_i_dep * N_i * pow ( m_i, 1. / 3. ) * ( c.x[ i ][ j ][ k ] - q_Ice );								// supersaturation
					}
					else  S_i_dep = 0.;

//					S_i_dep = 0.;


	// autoconversion processes
					if ( c_c_au * cloud.x[ i ][ j ][ k ] > 0. )								S_c_au = c_c_au * cloud.x[ i ][ j ][ k ];	// cloud water to rain, cloud droplet collection
					else 																			S_c_au = 0.;
					if ( c_i_au * ice.x[ i ][ j ][ k ] > 0. )									S_i_au = c_i_au * ice.x[ i ][ j ][ k ];		// cloud ice to snow, cloud ice crystal aggregation
					else 																			S_i_au = 0.;

					S_d_au = S_i_dep / ( 1.5 * ( pow ( m_s_0 / m_i, 2. / 3. ) - 1. ) );													// depositional growth of cloud ice

//					S_c_au = 0.;
//					S_i_au = 0.;
//					S_d_au = 0.;


// collection mechanism
					S_ac = c_ac * cloud.x[ i ][ j ][ k ] * pow ( P_rain.x[ i ][ j ][ k ], 7. / 9. ); 									// accreation rate from depletion of cloud water due to collection by all rain drops

					if ( t.x[ i ][ j ][ k ] * t_0 < t_0 ) 										S_rim = c_rim * cloud.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];
					if ( t.x[ i ][ j ][ k ] * t_0 >= t_0 ) 									S_rim = 0.;										// riming rate of snow mass due to collection of supercooled cloud droplets
																																								// by falling snow particles

					if ( t.x[ i ][ j ][ k ] * t_0 >= t_0 ) 									S_shed = c_rim * cloud.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];
					if ( t.x[ i ][ j ][ k ] * t_0 < t_0 ) 										S_shed = 0.;										// rate of water shed by melting wet snow particles
																																								// collecting cloud droplets to produce rain

//					S_ac = 0.;
//					S_rim = 0.;
//					S_shed = 0.;

					if ( t.x[ i ][ j ][ k ] * t_0 < t_0 )
					{
						S_agg = c_agg * ice.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];														// collection of cloud ice by snow particles

						S_i_cri = c_i_cri * ice.x[ i ][ j ][ k ] * pow ( P_rain.x[ i ][ j ][ k ], 7. / 9. );								// decrease in cloud ice mass due to collision/coalescense interaction with raindrops
						S_r_cri = c_r_cri * ice.x[ i ][ j ][ k ] / m_i * pow ( P_rain.x[ i ][ j ][ k ], 13. / 9. );						// decrease of rainwater due to freezing resulting from collection of ice crystals
					}
					else
					{
						S_agg = 0.;
						S_i_cri = 0.;
						S_r_cri = 0.;
					}

//					S_agg = 0.;
//					S_i_cri = 0.;
//					S_r_cri = 0.;


// diffusional growth of rain and snow
					if ( t.x[ i ][ j ][ k ] * t_0 >= t_0 )			S_ev = alf_ev * ( 1. + bet_ev * pow ( P_rain.x[ i ][ j ][ k ], 1. / 6. ) ) * ( q_Rain - c.x[ i ][ j ][ k ] ) * pow ( P_rain.x[ i ][ j ][ k ], 4. / 9. );
																																								// evaporation of rain due to water vapour diffusion
					else 												S_ev = 0.;

					if ( t.x[ i ][ j ][ k ] * t_0 < t_0 ) 			S_s_dep = c_s_dep * ( 1. + bet_s_dep * pow ( P_snow.x[ i ][ j ][ k ], 5. / 26. ) ) * ( c.x[ i ][ j ][ k ] - q_Ice ) * pow ( P_snow.x[ i ][ j ][ k ], 8. / 13. );
																																								// deposition/sublimation of snow 
					else 												S_s_dep = 0.;

//					S_ev = 0.;
//					S_s_dep = 0.;


// melting and freezing
					if ( ( t.x[ i ][ j ][ k ] * t_0 > t_0 ) && ( ice.x[ i ][ j ][ k ] > 0. ) ) 	S_i_melt = ice.x[ i ][ j ][ k ] / dt; // cloud ice particles melting to cloud water
					else 																						S_i_melt = 0.;

					p_t_in = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t_0 ) ) * p_SL;	// given in hPa
					t_Celsius_0 = 0.;
					E_Rain_t_in = hp * exp ( 17.0809 * t_Celsius_0 / ( 234.175 + t_Celsius_0 ) );								// saturation water vapour pressure for the water phase at t > 0°C in hPa
					q_Rain_t_in = ep * E_Rain_t_in / ( p_t_in - E_Rain_t_in );																// water vapour amount at saturation with water formation in kg/kg

					if ( t.x[ i ][ j ][ k ] * t_0 > t_0 ) 	S_s_melt = c_s_melt * ( 1. + b_s_melt * pow ( P_snow.x[ i ][ j ][ k ], 5. / 26. ) ) * ( ( t.x[ i ][ j ][ k ] * t_0 - t_0 ) + a_s_melt * ( c.x[ i ][ j ][ k ] - q_Rain_t_in ) ) * pow ( P_snow.x[ i ][ j ][ k ], 8. / 13. );																// melting rate of snow to form rain
					else 																							S_s_melt = 0.;

//					S_i_melt = 0.;
//					S_s_melt = 0.;

					if ( t.x[ i ][ j ][ k ] * t_0 <= t_0 )		S_if_frz = alf_if * ( exp ( a_if * ( t_0 - t.x[ i ][ j ][ k ] * t_0 ) ) - 1. ) * pow ( P_rain.x[ i ][ j ][ k ], 14. / 9. );	// freezing rate from immersion freezing
					else 											S_if_frz = 0.;

					if ( i <= 12 ) 																N_cf_0 = ( N_cf_0_6km - N_cf_0_surf ) / 12. * ( double ) i + N_cf_0_surf;	// linear distribution
					else 																			N_cf_0 = N_cf_0_6km;

					if ( t.x[ i ][ j ][ k ] * t_0 < 270.16 ) 									N_cf = N_cf_0 * pow ( 270.16 - t.x[ i ][ j ][ k ] * t_0, 1.3 );	// number density of contact nuclei
					if ( t.x[ i ][ j ][ k ] * t_0 >= 270.16 ) 								N_cf = 0.;

					S_cf_frz = alf_cf * E_cf * N_cf * pow ( P_rain.x[ i ][ j ][ k ], 13. / 9. );											// freezing rate from contact nucleation

//					S_if_frz = 0.;
//					S_cf_frz = 0.;


// sinks and sources
					S_v.x[ i ][ j ][ k ] = S_ev - S_i_dep - S_s_dep - S_nuc;
					S_c.x[ i ][ j ][ k ] = - S_c_au - S_ac - S_c_frz + S_i_melt - S_rim - S_shed - S_cf_frz;
					S_i.x[ i ][ j ][ k ] = S_nuc + S_c_frz + S_i_dep - S_i_melt - S_i_au - S_d_au - S_agg - S_i_cri + S_cf_frz;
					S_r.x[ i ][ j ][ k ] = S_c_au + S_ac - S_ev + S_shed - S_r_cri - S_if_frz + S_s_melt;
					S_s.x[ i ][ j ][ k ] = S_i_au + S_d_au + S_agg + S_rim + S_s_dep + S_i_cri + S_r_cri + S_if_frz - S_s_melt;


// rain and snow integration
					P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ] + S_r.x[ i ][ j ][ k ] * coeff_P;
					P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ] + S_s.x[ i ][ j ][ k ] * coeff_P;

					if ( t.x[ i ][ j ][ k ] * t_0 > t_0 ) 	P_snow.x[ i ][ j ][ k ] = 0.;

					if ( P_rain.x[ i ][ j ][ k ] < 0. )										P_rain.x[ i ][ j ][ k ] = 0.;
					if ( P_snow.x[ i ][ j ][ k ] < 0. )									P_snow.x[ i ][ j ][ k ] = 0.;

					if ( i <= im - 20 )
					{
						if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h.x[ i + 2 ][ j ][ k ] == 1. ) )
						{
							P_rain.x[ i ][ j ][ k ] = 0.;
							P_snow.x[ i ][ j ][ k ] = 0.;
						}
					}


					if ( c.x[ i ][ j ][ k ] < 0. ) 											c.x[ i ][ j ][ k ] = 0.;
					if ( cloud.x[ i ][ j ][ k ] < 0. ) 										cloud.x[ i ][ j ][ k ] = 0.;
					if ( ice.x[ i ][ j ][ k ] < 0. ) 											ice.x[ i ][ j ][ k ] = 0.;

				}																					// end i RainSnow
			}																						// end j
		}																							// end k

	}																								// end n









// calculation of a total quantity as sum on all values in a virtual column in r-direction
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( RadiationModel == 0 ) Q_Radiation.y[ j ][ k ] = Radiation_Balance_par.y[ j ][ k ]; // parabolic radiation balance assumed
			if ( RadiationModel >= 2 ) Q_Radiation.y[ j ][ k ] = radiation_3D.x[ 0 ][ j ][ k ];		 // two- and multi-layer radiation balance assumed

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
					if ( i == 0 ) 	p_stat.x[ 0 ][ j ][ k ] = ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01; // given in hPa
					else 	p_stat.x[ i ][ j ][ k ] = exp ( - g * ( double ) i * ( L_atm / ( double ) ( im-1 ) ) / ( R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) ) * p_stat.x[ 0 ][ j ][ k ];	// given in hPa
																																			// current air pressure, step size in 500 m, from a polytropic atmosphere in hPa

					if ( Latency.x[ i ][ j ][ k ] <= 0. ) 	t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0 + ( ( t_cond_3D.x[ i ][ j ][ k ] ) / t_0 );
					else  											t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0 + ( ( t_evap_3D.x[ i ][ j ][ k ] ) / t_0 );

					r_dry = 100. * p_stat.x[ i ][ j ][ k ] / ( R_Air * t.x[ i ][ j ][ k ] * t_0 );
					r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );				// density of humid air, COSMO version withot cloud and ice water, masses negligible

					e = c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep; 													// water vapour pressure in hPa
//					e = ( r_humid * R_WaterVapour * t.x[ i ][ j ][ k ] * t_0 ) * .01;								// delivers the same results
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
					Q_sensible.y[ j ][ k ] = - r_air * cp_l * coeff_Diffusion_sensibel * t_grad * t_0 / L_atm;		// sensible heat in [W/m2] from energy transport equation
					Q_bottom.y[ j ][ k ] = Q_Radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ];	// difference understood as heat into the ground

					Evaporation_Haude.y[ j ][ k ] = f_Haude * sat_deficit;											// simplified formula for Evaporation over day length of 12h by Haude, Häckel
					if ( Evaporation_Haude.y[ j ][ k ] <= 0. ) 		Evaporation_Haude.y[ j ][ k ] = 0.;
					Evaporation_Penman.y[ j ][ k ] = .0346 * ( ( Q_Radiation.y[ j ][ k ] - Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma );
					if ( Evaporation_Penman.y[ j ][ k ] <= 0. ) 		Evaporation_Penman.y[ j ][ k ] = 0.;	// .0346 coefficient W/m2 corresponds to mm/d (Kraus)

				}

//				Q_Sensible.x[ i ][ j ][ k ] = - cp_l * r_air * u_0 * t_0 / L_atm * ( u.x[ i ][ j ][ k ] * ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) + v.x[ i ][ j ][ k ] * ( t.x[ i ][ j+1 ][ k ] - t.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe ) / rm + w.x[ i ][ j ][ k ] * ( t.x[ i ][ j ][ k+1 ] - t.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi ) / rmsinthe );							// sensible heat in [W/m2] from energy transport equation
/*
				if ( Latency.x[ i ][ j ][ k ] <= 0. )		t_cond_3D.x[ i ][ j ][ k ] = pow ( ( fabs ( radiation_3D.x[ i ][ j ][ k ] ) ) / sigma, .25 );																								// temperature increase due to condensation
				if ( Latency.x[ i ][ j ][ k ] <= 0. )		t_cond_3D.x[ i ][ j ][ k ] = pow ( ( fabs ( Latency.x[ i ][ j ][ k ] ) ) / sigma, .25 );																								// temperature increase due to condensation
				if ( Latency.x[ i ][ j ][ k ] <= 0. )		t_cond_3D.x[ i ][ j ][ k ] = pow ( ( fabs ( Q_Sensible.x[ i ][ j ][ k ] ) ) / sigma, .25 );																								// temperature increase due to condensation
*/
/*
				if ( Latency.x[ i ][ j ][ k ] <= 0. )		t_cond_3D.x[ i ][ j ][ k ] = pow ( ( radiation_3D.x[ i ][ j ][ k ] - Latency.x[ i ][ j ][ k ] - Q_Sensible.x[ i ][ j ][ k ] ) / sigma, .25 ) - pow ( radiation_3D.x[ i ][ j ][ k ] / sigma, .25 );																							// temperature increase due to condensation
				else 												t_cond_3D.x[ i ][ j ][ k ] = 0.;

				if ( Latency.x[ i ][ j ][ k ] > 0. )			t_evap_3D.x[ i ][ j ][ k ] = pow ( ( radiation_3D.x[ i ][ j ][ k ] - Latency.x[ i ][ j ][ k ] - Q_Sensible.x[ i ][ j ][ k ] ) / sigma, .25 ) - pow ( radiation_3D.x[ i ][ j ][ k ] / sigma, .25 );																								// temperature decrease due to evaporation
				else 												t_evap_3D.x[ i ][ j ][ k ] = 0.;

				if ( h.x[ i ][ j ][ k ] == 1. )					t_cond_3D.x[ i ][ j ][ k ] = t_evap_3D.x[ i ][ j ][ k ] = 0.;
*/
//	cout << i << "   " << j << "   " << k << "   " << radiation_3D.x[ i ][ j ][ k ] << "   " << Latency.x[ i ][ j ][ k ] << "   " << Q_Sensible.x[ i ][ j ][ k ] << "   " << t_cond_3D.x[ i ][ j ][ k ] << "   " << t_evap_3D.x[ i ][ j ][ k ] << endl;


// only on the sea surface
				if ( ( i == 0 ) && ( h.x[ 0 ][ j ][ k ] == 0. ) )
				{
					if ( i == 0 ) 	p_stat.x[ 0 ][ j ][ k ] = ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;		// given in hPa

					if ( Latency.x[ 0 ][ j ][ k ] <= 0. ) 	t_Celsius = t.x[ 0 ][ j ][ k ] * t_0 - t_0 + ( ( t_cond_3D.x[ 0 ][ j ][ k ] ) / t_0 );
					else  											t_Celsius = t.x[ 0 ][ j ][ k ] * t_0 - t_0 + ( ( t_evap_3D.x[ 0 ][ j ][ k ] ) / t_0 );

					r_dry = 100. * p_stat.x[ 0 ][ j ][ k ] / ( R_Air * t.x[ 0 ][ j ][ k ] * t_0 );
					r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );				// density of humid air, COSMO version withot cloud and ice water, masses negligible

					e = c.x[ 0 ][ j ][ k ] * p_stat.x[ 0 ][ j ][ k ] / ep; 													// water vapour pressure in hPa
//					e = ( r_humid * R_WaterVapour * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;								// delivers the same results
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
					Q_sensible.y[ j ][ k ] = - r_air * cp_l * coeff_Diffusion_sensibel * t_grad * t_0 / L_atm;		// sensible heat in [W/m2] from energy transport equation
					Q_bottom.y[ j ][ k ] = Q_Radiation.y[ j ][ k ] - Q_latent.y[ j ][ k ] - Q_sensible.y[ j ][ k ];	// difference understood as heat of the ground

					Evaporation_Haude.y[ j ][ k ] = f_Haude * sat_deficit;											// simplified formula for Evaporation over day length of 12h by Haude, Häckel
					if ( Evaporation_Haude.y[ j ][ k ] <= 0. ) 		Evaporation_Haude.y[ j ][ k ] = 0.;
					Evaporation_Penman.y[ j ][ k ] = .0346 * ( ( Q_Radiation.y[ j ][ k ] - Q_bottom.y[ j ][ k ] ) * Delta + gamma * E_a ) / ( Delta + gamma );
					if ( Evaporation_Penman.y[ j ][ k ] <= 0. ) 		Evaporation_Penman.y[ j ][ k ] = 0.;	// .0346 coefficient W/m2 corresponds to mm/d (Kraus)

//					Q_Sensible.x[ 0 ][ j ][ k ] = - cp_l * r_air * u_0 * t_0 / L_atm * ( v.x[ 0 ][ j ][ k ] * ( t.x[ 0 ][ j+1 ][ k ] - t.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe ) / rm + w.x[ 0 ][ j ][ k ] * ( t.x[ 0 ][ j ][ k+1 ] - t.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi ) / rmsinthe );							// sensible heat in [W/m2] from energy transport equation
/*
					if ( Latency.x[ 0 ][ j ][ k ] <= 0. )		t_cond_3D.x[ 0 ][ j ][ k ] = pow ( ( radiation_3D.x[ 0 ][ j ][ k ] - Latency.x[ 0 ][ j ][ k ] - Q_Sensible.x[ 0 ][ j ][ k ] ) / sigma, .25 ) - pow ( radiation_3D.x[ 0 ][ j ][ k ] / sigma, .25 );																							// temperature increase due to condensation
					else 													t_cond_3D.x[ 0 ][ j ][ k ] = 0.;

					if ( Latency.x[ 0 ][ j ][ k ] > 0. )			t_evap_3D.x[ 0 ][ j ][ k ] = pow ( ( radiation_3D.x[ 0 ][ j ][ k ] - Latency.x[ 0 ][ j ][ k ] - Q_Sensible.x[ 0 ][ j ][ k ] ) / sigma, .25 ) - pow ( radiation_3D.x[ 0 ][ j ][ k ] / sigma, .25 );																								// temperature decrease due to evaporation
					else 													t_evap_3D.x[ 0 ][ j ][ k ] = 0.;

					if ( h.x[ 0 ][ j ][ k ] == 1. )					t_cond_3D.x[ 0 ][ j ][ k ] = t_evap_3D.x[ 0 ][ j ][ k ] = 0.;
*/
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

			Q_Sensible.x[ 0 ][ j ][ k ] = c43 * Q_Sensible.x[ 1 ][ j ][ k ] - c13 * Q_Sensible.x[ 2 ][ j ][ k ];
			Q_Sensible.x[ im-1 ][ j ][ k ] = c43 * Q_Sensible.x[ im-2 ][ j ][ k ] - c13 * Q_Sensible.x[ im-3 ][ j ][ k ];

			t_cond_3D.x[ 0 ][ j ][ k ] = c43 * t_cond_3D.x[ 1 ][ j ][ k ] - c13 * t_cond_3D.x[ 2 ][ j ][ k ];
			t_cond_3D.x[ im-1 ][ j ][ k ] = c43 * t_cond_3D.x[ im-2 ][ j ][ k ] - c13 * t_cond_3D.x[ im-3 ][ j ][ k ];

			t_evap_3D.x[ 0 ][ j ][ k ] = c43 * t_evap_3D.x[ 1 ][ j ][ k ] - c13 * t_evap_3D.x[ 2 ][ j ][ k ];
			t_evap_3D.x[ im-1 ][ j ][ k ] = c43 * t_evap_3D.x[ im-2 ][ j ][ k ] - c13 * t_evap_3D.x[ im-3 ][ j ][ k ];

			BuoyancyForce.x[ 0 ][ j ][ k ] = c43 * BuoyancyForce.x[ 1 ][ j ][ k ] - c13 * BuoyancyForce.x[ 2 ][ j ][ k ];
			BuoyancyForce.x[ im-1 ][ j ][ k ] = c43 * BuoyancyForce.x[ im-2 ][ j ][ k ] - c13 * BuoyancyForce.x[ im-3 ][ j ][ k ];

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

			Q_Sensible.x[ i ][ 0 ][ k ] = c43 * Q_Sensible.x[ i ][ 1 ][ k ] - c13 * Q_Sensible.x[ i ][ 2 ][ k ];
			Q_Sensible.x[ i ][ jm-1 ][ k ] = c43 * Q_Sensible.x[ i ][ jm-2 ][ k ] - c13 * Q_Sensible.x[ i ][ jm-3 ][ k ];

			t_cond_3D.x[ i ][ 0 ][ k ] = c43 * t_cond_3D.x[ i ][ 1 ][ k ] - c13 * t_cond_3D.x[ i ][ 2 ][ k ];
			t_cond_3D.x[ i ][ jm-1 ][ k ] = c43 * t_cond_3D.x[ i ][ jm-2 ][ k ] - c13 * t_cond_3D.x[ i ][ jm-3 ][ k ];

			t_evap_3D.x[ i ][ 0 ][ k ] = c43 * t_evap_3D.x[ i ][ 1 ][ k ] - c13 * t_evap_3D.x[ i ][ 2 ][ k ];
			t_evap_3D.x[ i ][ jm-1 ][ k ] = c43 * t_evap_3D.x[ i ][ jm-2 ][ k ] - c13 * t_evap_3D.x[ i ][ jm-3 ][ k ];

			BuoyancyForce.x[ i ][ 0 ][ k ] = c43 * BuoyancyForce.x[ i ][ 1 ][ k ] - c13 * BuoyancyForce.x[ i ][ 2 ][ k ];
			BuoyancyForce.x[ i ][ jm-1 ][ k ] = c43 * BuoyancyForce.x[ i ][ jm-2 ][ k ] - c13 * BuoyancyForce.x[ i ][ jm-3 ][ k ];


//			c.x[ i ][ 0 ][ k ] = c43 * c.x[ i ][ 1 ][ k ] - c13 * c.x[ i ][ 2 ][ k ];
//			c.x[ i ][ jm-1 ][ k ] = c43 * c.x[ i ][ jm-2 ][ k ] - c13 * c.x[ i ][ jm-3 ][ k ];
			if ( ( c.x[ i ][ 0 ][ k ] < 0. ) || ( c.x[ i ][ jm-1 ][ k ] < 0. ) ) c.x[ i ][ 0 ][ k ] = c.x[ i ][ jm-1 ][ k ] = 0.;

//			cloud.x[ i ][ 0 ][ k ] = c43 * cloud.x[ i ][ 1 ][ k ] - c13 * cloud.x[ i ][ 2 ][ k ];
//			cloud.x[ i ][ jm-1 ][ k ] = c43 * cloud.x[ i ][ jm-2 ][ k ] - c13 * cloud.x[ i ][ jm-3 ][ k ];
			if ( ( cloud.x[ i ][ 0 ][ k ] < 0. ) || ( cloud.x[ i ][ jm-1 ][ k ] < 0. ) ) cloud.x[ i ][ 0 ][ k ] = cloud.x[ i ][ jm-1 ][ k ] = 0.;

//			ice.x[ i ][ 0 ][ k ] = c43 * ice.x[ i ][ 1 ][ k ] - c13 * ice.x[ i ][ 2 ][ k ];
//			ice.x[ i ][ jm-1 ][ k ] = c43 * ice.x[ i ][ jm-2 ][ k ] - c13 * ice.x[ i ][ jm-3 ][ k ];
			if ( ( ice.x[ i ][ 0 ][ k ] < 0. ) || ( ice.x[ i ][ jm-1 ][ k ] < 0. ) ) ice.x[ i ][ 0 ][ k ] = ice.x[ i ][ jm-1 ][ k ] = 0.;

			P_rain.x[ i ][ 0 ][ k ] = c43 * P_rain.x[ i ][ 1 ][ k ] - c13 * P_rain.x[ i ][ 2 ][ k ];
			P_rain.x[ i ][ jm-1 ][ k ] = c43 * P_rain.x[ i ][ jm-2 ][ k ] - c13 * P_rain.x[ i ][ jm-3 ][ k ];
			if ( ( P_rain.x[ i ][ 0 ][ k ] < 0. ) || ( P_rain.x[ i ][ jm-1 ][ k ] < 0. ) ) P_rain.x[ i ][ 0 ][ k ] = P_rain.x[ i ][ jm-1 ][ k ] = 0.;

			P_snow.x[ i ][ 0 ][ k ] = c43 * P_snow.x[ i ][ 1 ][ k ] - c13 * P_snow.x[ i ][ 2 ][ k ];
			P_snow.x[ i ][ jm-1 ][ k ] = c43 * P_snow.x[ i ][ jm-2 ][ k ] - c13 * P_snow.x[ i ][ jm-3 ][ k ];
			if ( ( P_snow.x[ i ][ 0 ][ k ] < 0. ) || ( P_snow.x[ i ][ jm-1 ][ k ] < 0. ) ) P_snow.x[ i ][ 0 ][ k ] = P_snow.x[ i ][ jm-1 ][ k ] = 0.;
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

//			cloud.x[ i ][ j ][ 0 ] = c43 * cloud.x[ i ][ j ][ 1 ] - c13 * cloud.x[ i ][ j ][ 2 ];
//			cloud.x[ i ][ j ][ km-1 ] = c43 * cloud.x[ i ][ j ][ km-2 ] - c13 * cloud.x[ i ][ j ][ km-3 ];
//			if ( ( cloud.x[ i ][ j ][ 0 ] < 0. ) || ( cloud.x[ i ][ j ][ km-1 ] < 0. ) ) cloud.x[ i ][ j ][ 0 ] = cloud.x[ i ][ j ][ km-1 ] = 0.;
//			cloud.x[ i ][ j ][ 0 ] = cloud.x[ i ][ j ][ km-1 ] = ( cloud.x[ i ][ j ][ 0 ] + cloud.x[ i ][ j ][ km-1 ] ) / 2.;

//			ice.x[ i ][ j ][ 0 ] = c43 * ice.x[ i ][ j ][ 1 ] - c13 * ice.x[ i ][ j ][ 2 ];
//			ice.x[ i ][ j ][ km-1 ] = c43 * ice.x[ i ][ j ][ km-2 ] - c13 * ice.x[ i ][ j ][ km-3 ];
//			if ( ( ice.x[ i ][ j ][ 0 ] < 0. ) || ( ice.x[ i ][ j ][ km-1 ] < 0. ) ) ice.x[ i ][ j ][ 0 ] = ice.x[ i ][ j ][ km-1 ] = 0.;
//			ice.x[ i ][ j ][ 0 ] = ice.x[ i ][ j ][ km-1 ] = ( ice.x[ i ][ j ][ 0 ] + ice.x[ i ][ j ][ km-1 ] ) / 2.;

			P_rain.x[ i ][ j ][ 0 ] = c43 * P_rain.x[ i ][ j ][ 1 ] - c13 * P_rain.x[ i ][ j ][ 2 ];
			P_rain.x[ i ][ j ][ km-1 ] = c43 * P_rain.x[ i ][ j ][ km-2 ] - c13 * P_rain.x[ i ][ j ][ km-3 ];
			if ( ( P_rain.x[ i ][ j ][ 0 ] < 0. ) || ( P_rain.x[ i ][ j ][ km-1 ] < 0. ) ) P_rain.x[ i ][ j ][ 0 ] = P_rain.x[ i ][ j ][ km-1 ] = 0.;
			P_rain.x[ i ][ j ][ 0 ] = P_rain.x[ i ][ j ][ km-1 ] = ( P_rain.x[ i ][ j ][ 0 ] + P_rain.x[ i ][ j ][ km-1 ] ) / 2.;

			P_snow.x[ i ][ j ][ 0 ] = c43 * P_snow.x[ i ][ j ][ 1 ] - c13 * P_snow.x[ i ][ j ][ 2 ];
			P_snow.x[ i ][ j ][ km-1 ] = c43 * P_snow.x[ i ][ j ][ km-2 ] - c13 * P_snow.x[ i ][ j ][ km-3 ];
			if ( ( P_snow.x[ i ][ j ][ 0 ] < 0. ) || ( P_snow.x[ i ][ j ][ km-1 ] < 0. ) ) P_snow.x[ i ][ j ][ 0 ] = P_snow.x[ i ][ j ][ km-1 ] = 0.;
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

				precipitable_water.y[ j ][ k ] += a * L_atm / ( double ) ( im - 1 );						//  kg/m³ * m

				if ( h.x[ i ][ j ][ k ] == 1. )		cloud.x[ i ][ j ][ k ] = BuoyancyForce.x[ i ][ j ][ k ] = t_cond_3D.x[ i ][ j ][ k ] = t_evap_3D.x[ i ][ j ][ k ] = Latency.x[ i ][ j ][ k ] = 0.;

			}
		}
	}



// surface values of precipitation and precipitable water
	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Precipitation.y[ j ][ k ] = 1000. * ( P_rain.x[ 0 ][ j ][ k ] + P_snow.x[ 0 ][ j ][ k ] );

			precipitable_water.y[ j ][ k ] = precipitable_water.y[ j ][ k ] / 1000.;// divided by water density ( 1000 kg/m³ ) results in m compares as well to mm ( absolute values identical )
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

			Evaporation_Penman_average += Evaporation_Penman.y[ j ][ k ];
			Evaporation_Haude_average += Evaporation_Haude.y[ j ][ k ];

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
	}



	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				if ( cloud.x[ i ][ j ][ k ] < 0. ) 		cloud.x[ i ][ j ][ k ] = 0.;
				if ( ice.x[ i ][ j ][ k ] < 0. ) 		ice.x[ i ][ j ][ k ] = 0.;
				if ( P_rain.x[ i ][ j ][ k ] < 0. ) 		P_rain.x[ i ][ j ][ k ] = 0.;
				if ( P_snow.x[ i ][ j ][ k ] < 0. ) 		P_snow.x[ i ][ j ][ k ] = 0.;
			}
		}
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

/*
	Value_10 = precipitation_NASA_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_10 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_10 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_11 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_10 / 365. << setw ( 6 ) << name_unit_mmd << endl;
*/

	Value_9 = co2_vegetation_average * co2_0;
	Value_12 = Evaporation_Penman_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_22 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_9 << setw ( 6 ) << name_unit_ppm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_12 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_12 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_13 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_12 / 365. << setw ( 6 ) << name_unit_mmd << endl;

/*
	Value_13 = Evaporation_Haude_average;

	cout << setw ( 6 ) << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_7 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_7 << setw ( 6 ) << name_unit_mm << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_14 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_13 << setw ( 6 ) << name_unit_mma << "   " << setiosflags ( ios::left ) << setw ( 40 ) << setfill ( '.' ) << name_Value_15 << " = " << resetiosflags ( ios::left ) << setw ( 7 ) << fixed << setfill ( ' ' ) << Value_13 / 365. << setw ( 6 ) << name_unit_mmd << endl << endl << endl;
*/



}




