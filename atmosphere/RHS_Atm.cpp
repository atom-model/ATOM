/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/


#include <iostream>
#include <cmath>

#include "RHS_Atm.h"

using namespace std;





RHS_Atmosphere::RHS_Atmosphere ( int im, int jm, int km, double dt, double dr, double dthe, double dphi, double re, double ec, double sc_WaterVapour, double sc_CO2, double gr, double pr, double omega, double coriolis, double centrifugal, double WaterVapour, double buoyancy, double CO2 )
{
	this-> dt = dt;
	this-> dr = dr;
	this-> dthe = dthe;
	this-> dphi = dphi;
	this-> re = re;
	this-> ec = ec;
	this-> sc_WaterVapour = sc_WaterVapour;
	this-> sc_CO2 = sc_CO2;
	this-> gr = gr;
	this-> pr = pr;
	this-> omega = omega;
	this-> coriolis = coriolis;
	this-> centrifugal = centrifugal;
	this-> WaterVapour = WaterVapour;
	this-> buoyancy = buoyancy;
	this-> CO2 = CO2;
}


RHS_Atmosphere::~RHS_Atmosphere() {}



void RHS_Atmosphere::RK_RHS_Atmosphere ( int i, int j, int k, double lv, double ls, double ep, double hp, double u_0, double t_0, double t_Boussinesq, double c_0, double co2_0, double p_0, double r_0_air, double r_0_water_vapour, double r_0_co2, double L_atm, double cp_l, double R_Air, double R_WaterVapour, double R_co2, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &co2n, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c, Array &rhs_co2, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Rain, Array &Ice, Array &Rain_super, Array &IceLayer )
{
// collection of coefficients for phase transformation

	L_atm = L_atm / ( im - 1 );														// characteristic length for non-dimensionalisation

	coeff_lv = ( lv * c_0 ) / ( r_0_air * cp_l * t_0 * u_0 );				// coeff_lv = 0.009522
	coeff_ls = ( ls * c_0 ) / (  r_0_air * cp_l * t_0 * u_0 );				// coeff_ls = 0.008916518

	c43 = 4. / 3.;
	c13 = 1. / 3.;

//	k_Force = 1.;																			// factor for accelleration of convergence processes inside the immersed boundary conditions
	k_Force = 10.;																			// factor for accelleration of convergence processes inside the immersed boundary conditions

	cc = + 1.;

	h_check_i = h_check_j = h_check_k = 0;


// 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )

// collection of coefficients

	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

// collection of coefficients

	rm = rad.z[ i ];
	rm2 = rm * rm;

// collection of coefficients

	sinthe = sin( the.z[ j ] );
	sinthe2 = sinthe * sinthe;
	costhe = cos( the.z[ j ] );
	cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
	rmsinthe = rm * sinthe;
	rm2sinthe = rm2 * sinthe;
	rm2sinthe2 = rm2 * sinthe2;


//  3D volume iterations in case 1. and 2. order derivatives at walls are needed >>>>>>>>>>>>>>>>>>>>>>>>
// only in positive r-direction above ground 

	if ( ( h.x[ i - 1 ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k ] == 0. ) )
	{
		h_0_i = ( double ) ( i ) * dr - dr / 4.;

		if ( fabs ( ( ( double ) ( i + 1 ) * dr - h_0_i ) ) < dr )
		{
			h_c_i = cc * ( 1. - fabs ( ( ( double ) ( i + 1 ) * dr - h_0_i ) ) / dr ); 
//			h_c_i = cc * ( .5 * ( acos ( fabs ( ( double ) ( i + 1 ) * dr - h_0_i ) * 3.14 / dr ) + 1. ) ); 
			h_d_i = 1. - h_c_i;
			h_check_i = 1;
		}
	}


// 3D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in positive the-direction along northerly boundaries 

	if ( ( h.x[ i ][ j - 1 ][ k ] == 1. ) && ( h.x[ i ][ j ][ k ] == 0. ) )
	{
		h_0_j = ( double ) ( j ) * dthe + dthe / 4.;

		if ( fabs ( ( ( double ) ( j + 1 ) * dthe - h_0_j ) ) < dthe )
		{
			h_c_j = cc * ( 1. - fabs ( ( ( double ) ( j + 1 ) * dthe - h_0_j ) ) / dthe ); 
//			h_c_j = cc * ( .5 * ( acos ( fabs ( ( double ) ( j + 1 ) * dthe - h_0_j ) * 3.14 / dthe ) + 1. ) ); 
			h_d_j = 1. - h_c_j;
			h_check_j = 1;
		}
	}


// 3D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in negative the-direction along southerly boundaries 

	if ( ( h.x[ i ][ j + 1 ][ k ] == 1. ) && ( h.x[ i ][ j ][ k ] == 0. ) )
	{
		h_0_j = ( double ) ( j ) * dthe - dthe / 4.;

		if ( fabs ( ( ( double ) ( j - 1 ) * dthe - h_0_j ) ) < dthe )
		{
			h_c_j = cc * ( 1. - fabs ( ( ( double ) ( j - 1 ) * dthe - h_0_j ) ) / dthe ); 
//			h_c_j = cc * ( .5 * ( acos ( fabs ( ( double ) ( j - 1 ) * dthe - h_0_j ) * 3.14 / dthe ) + 1. ) ); 
			h_d_j = 1. - h_c_j;
			h_check_j = 1;
		}
	}


// 3D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in positive phi-direction on westerly boundaries 

	if ( ( h.x[ i ][ j ][ k - 1 ] == 1. ) && ( h.x[ i ][ j ][ k ] == 0. ) )
	{
		h_0_k = ( double ) ( k ) * dphi + dphi / 4.;

		if ( fabs ( ( ( double ) ( k + 1 ) * dphi - h_0_k ) ) < dphi )
		{
			h_c_k = cc * ( 1. - fabs ( ( ( double ) ( k + 1 ) * dphi - h_0_k ) ) / dphi ); 
//			h_c_k = cc * ( .5 * ( acos ( fabs ( ( double ) ( k + 1 ) * dphi - h_0_k ) * 3.14 / dphi ) + 1. ) ); 
			h_d_k = 1. - h_c_k;
			h_check_k = 1;
		}
	}


// 3D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in negative phi-direction along easterly boundaries 

	if ( ( h.x[ i ][ j ][ k + 1 ] == 1. ) && ( h.x[ i ][ j ][ k ] == 0. ) )
	{
		h_0_k = ( double ) ( k ) * dphi - dphi / 4.;

		if ( fabs ( ( ( double ) ( k - 1 ) * dphi - h_0_k ) ) < dphi )
		{
			h_c_k = cc * ( 1. - fabs ( ( ( double ) ( k - 1 ) * dphi - h_0_k ) ) / dphi ); 
//			h_c_k = cc * ( .5 * ( acos ( fabs ( ( double ) ( k - 1 ) * dphi - h_0_k ) * 3.14 / dphi ) + 1. ) ); 
			h_d_k = 1. - h_c_k;
			h_check_k = 1;
		}
	}


		if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h_check_i != 1 ) )
		{
			h_c_i = 0.; 
			h_d_i = 1. - h_c_i;
		}

		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h_check_i != 1 ) )
		{
			h_c_i = 1.; 
			h_d_i = 1. - h_c_i;
		}


		if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h_check_j != 1 ) )
		{
			h_c_j = 0.; 
			h_d_j = 1. - h_c_j;
		}

		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h_check_j != 1 ) )
		{
			h_c_j = 1.; 
			h_d_j = 1. - h_c_j;
		}


		if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h_check_k != 1 ) )
		{
			h_c_k = 0.; 
			h_d_k = 1. - h_c_k;
		}

		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h_check_k != 1 ) )
		{
			h_c_k = 1.; 
			h_d_k = 1. - h_c_k;
		}



// 1st order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components

	dudr = h_d_i * ( u.x[ i+1 ][ j ][ k ] - u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dvdr = h_d_i * ( v.x[ i+1 ][ j ][ k ] - v.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dwdr = h_d_i * ( w.x[ i+1 ][ j ][ k ] - w.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dtdr = h_d_i * ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dpdr = h_d_i * ( p.x[ i+1 ][ j ][ k ] - p.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dcdr = h_d_i * ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dco2dr = h_d_i * ( co2.x[ i+1 ][ j ][ k ] - co2.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );

	dudthe = h_d_j * ( u.x[ i ][ j+1 ][ k ] - u.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dvdthe = h_d_j * ( v.x[ i ][ j+1 ][ k ] - v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dwdthe = h_d_j * ( w.x[ i ][ j+1 ][ k ] - w.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dtdthe = h_d_j * ( t.x[ i ][ j+1 ][ k ] - t.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dpdthe = h_d_j * ( p.x[ i ][ j+1 ][ k ] - p.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dcdthe = h_d_j * ( c.x[ i ][ j+1 ][ k ] - c.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dco2dthe = h_d_j * ( co2.x[ i ][ j+1 ][ k ] - co2.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );


	dudphi = h_d_k * ( u.x[ i ][ j ][ k+1 ] - u.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dvdphi = h_d_k * ( v.x[ i ][ j ][ k+1 ] - v.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dwdphi = h_d_k * ( w.x[ i ][ j ][ k+1 ] - w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dtdphi = h_d_k * ( t.x[ i ][ j ][ k+1 ] - t.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dpdphi = h_d_k * ( p.x[ i ][ j ][ k+1 ] - p.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dcdphi = h_d_k * ( c.x[ i ][ j ][ k+1 ] - c.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dco2dphi = h_d_k * ( co2.x[ i ][ j ][ k+1 ] - co2.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );

// 2nd order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components

	d2udr2 = h_d_i * ( u.x[ i+1 ][ j ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] ) / dr2;
	d2vdr2 = h_d_i * ( v.x[ i+1 ][ j ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] ) / dr2;
	d2wdr2 = h_d_i * ( w.x[ i+1 ][ j ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] ) / dr2;
	d2tdr2 = h_d_i * ( t.x[ i+1 ][ j ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i-1 ][ j ][ k ] ) / dr2;
	d2cdr2 = h_d_i * ( c.x[ i+1 ][ j ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i-1 ][ j ][ k ] ) / dr2;
	d2co2dr2 = h_d_i * ( co2.x[ i+1 ][ j ][ k ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i-1 ][ j ][ k ] ) / dr2;

	d2udthe2 = h_d_j * ( u.x[ i ][ j+1 ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2vdthe2 = h_d_j * ( v.x[ i ][ j+1 ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2wdthe2 = h_d_j * ( w.x[ i ][ j+1 ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2tdthe2 = h_d_j * ( t.x[ i ][ j+1 ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2cdthe2 = h_d_j * ( c.x[ i ][ j+1 ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2co2dthe2 = h_d_j * ( co2.x[ i ][ j+1 ][ k ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i ][ j-1 ][ k ] ) / dthe2;

	d2udphi2 = h_d_k * ( u.x[ i ][ j ][ k+1 ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2vdphi2 = h_d_k * ( v.x[ i ][ j ][ k+1 ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2wdphi2 = h_d_k * ( w.x[ i ][ j ][ k+1 ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2tdphi2 = h_d_k * ( t.x[ i ][ j ][ k+1 ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2cdphi2 = h_d_k * ( c.x[ i ][ j ][ k+1 ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2co2dphi2 = h_d_k * ( co2.x[ i ][ j ][ k+1 ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i ][ j ][ k-1 ] ) / dphi2;



// determination of temperature and pressure by the law of Clausius-Clapeyron for water vapour concentration
// reaching saturation of water vapour pressure leads to formation of rain or ice
// precipitation and cloud formation by formulas from Häckel
// dry adiabatic lapse rate and saturated adiabatic lapse rate = temperature decrease with hight
// SL stands for sea level

	t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;																						// conversion from Kelvin to Celsius
	t_Celsius_SL = t.x[ 0 ][ j ][ k ] * t_0 - t_0;																				// conversion from Kelvin to Celsius at sea surface = NN

	p_SL = ( r_0_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;																// surface pressure from barometric formula in hPa
//	cout << r_0_air << "   " << R_Air << "   " << t_0 << "   " << t.x[ 0 ][ j ][ k ] << "   " << p_SL << endl;
	p_h = exp ( - 9.8066 * ( double ) i * 500. / ( R_Air * t.x[ i ][ j ][ k ] * t_0 ) ) * p_0;					// current air pressure, step size in 500 m, from barometric formula in hPa

	E_Rain_SL = hp * exp ( 17.0809 * t_Celsius_SL / ( 234.175 + t_Celsius_SL ) );					// saturation vapour pressure in the water phase for t > 0°C in hPa
	E_Rain = hp * exp ( 17.0809 * t_Celsius / ( 234.175 + t_Celsius ) );										// saturation water vapour pressure for the water phase at t > 0°C
	E_Rain_super = hp * exp ( 17.8436 * t_Celsius / ( 245.425 + t_Celsius ) );							// saturation water vapour pressure for the water phase at t < 0°C, supercooled
	E_Ice = hp * exp ( 22.4429 * t_Celsius / ( 272.44 + t_Celsius ) );											// saturation water vapour pressure for the ice phase

	q_Rain  = ep * E_Rain / p_h;																									// water vapour amount at saturation with water formation
	q_Rain_super  = ep * E_Rain_super / p_h;																			// water vapour amount at saturation with water formation
	q_Ice  = ep * E_Ice / p_h;																										// water vapour amount at saturation with ice formation

	e_SL = ( r_0_water_vapour * R_WaterVapour * t.x[ 0 ][ j ][ k ] * t_0 ) * .01;							// water vapour pressure  in hPa at sea level
	a_SL = 216.6 * e_SL / ( t.x[ 0 ][ j ][ k ] * t_0 );																			// absolute humidity in kg/m3 at sea level
	q_SL = ep * c.x[ 0 ][ j ][ k ];																										// threshold value for water vapour at sea level in kg/kg
//	q_SL  = ep * e_SL / p_SL;																										// water vapour amount at sea level

	t_tau_SL = ( 423.86 - 234.175 * log ( e_SL ) ) / ( log ( e_SL ) - 18.89 );								// dewpoint temperature on ground in °C 		by Häckel
	h_level = 122. * ( t_Celsius_SL - t_tau_SL );																			// condensation level in m		by Häckel + correction
	h_h = - R_Air * t.x[ i ][ j ][ k ] * t_0 / gr * log ( p_h / p_0 );														// barometric elevation formula solved for h corresponding to hight over ground
	i_level = ( int ) ( h_level / 500. );																								// condensation level in radial steps

// precipitation and cloud formation from formulas by Häckel
// h stands for the local position

	e_h = ( r_0_water_vapour * R_WaterVapour * t.x[ i ][ j ][ k ] * t_0 ) * .01;								// water vapour pressure  in hPa
	a_h = 216.6 * e_h / ( t.x[ i ][ j ][ k ] * t_0 );																				// absolute humidity in kg/m3
	q_h = ep * c.x[ i ][ j ][ k ];																											// threshold value for water vapour at local hight h in kg/kg
//	q_h  = ep * e_h / p_h;																												// water vapour amount at local hight h

	t_tau_h = ( 423.86 - 234.175 * log ( e_h ) ) / ( log ( e_h ) - 18.89 );										// current dewpoint temperature in °C
	sat_Deficit = E_Rain - e_h;																										// saturation deficit, if positive then saturation is less than 100%
	RF_e = e_h / E_Rain * 100.;																									// relative humidity at any point in %
	Evaporation_Haude = 3.6 * sat_Deficit;																					// simplified formula for evaporation over day length of 12h by Haude ( from Häckel )
																																					// coefficient per day about 0.36 for gras-evapotranspiration given in (Kuttler)

// application of threshhold values for water vapour to compute rain, super cooled water and ice
//	if ( ( q_h > q_Rain ) && ( t_Celsius >= 0. ) ) 	Rain.x[ i ][ j ][ k ] = ( q_h - q_Rain );				// liquid water as surplus of local and saturated water vapour
																				// only the latent heat of water vapour is used!!!!!!!
	if ( q_h > q_Rain )	Rain.x[ i ][ j ][ k ] = ( q_h - q_Rain );														// liquid water as surplus of local and saturated water vapour
	else 	Rain.x[ i ][ j ][ k ] = 0.;																									// only positiv saturation level for water vapour is used

	if ( ( q_h > q_Rain_super ) && ( t_Celsius <= 0. ) && ( t_Celsius >= - 12. ) )	Rain_super.x[ i ][ j ][ k ] = ( q_h - q_Rain_super ); // supercooled  saturated water vapour 
	else 	Rain_super.x[ i ][ j ][ k ] = 0.;																						// only positiv saturation level for water vapour is used

	if ( ( q_h > q_Ice ) && ( t_Celsius < -12. ) && ( c.x[ i ][ j ][ k ] > 0. ) )	Ice.x[ i ][ j ][ k ] = ( q_h - q_Ice );	// ice formation as surplus of local and saturated water vapour
	else  Ice.x[ i ][ j ][ k ] = 0.;																										// only positiv saturation level for water vapour is used



/*
// printout for various thermodynamical quantities for the preticipation computations along the equator ( j = 90 )
	if ( ( i == 5 ) && ( j == 90 ) && ( k == 180 ) )
	{
		cout << endl;
		cout << " i = " << i << "   j = " << j << "   k = " << k << "   i_level = " << i_level  << "   h_level (m) = " << h_level << "   h_h (m) = " << h_h << endl << endl;

		cout << " t_h (°C) = " << t_Celsius << "   p_h (hPa) = " << p_h << "   a_h (g/m3) = " << a_h << "   c_h (hPa/hPa) = " << c.x[ i ][ j ][ k ] << "   q_Rain (hPa/hPa) = " << q_Rain << "   e_h (hPa) = " << e_h << "   E_Rain (hPa) = " << E_Rain  << endl << endl;

		cout << " Rain (g/kg) = " << Rain.x[ i ][ j ][ k ] << "   Rain_super (g/kg) = " << Rain_super.x[ i ][ j ][ k ] << "   Ice (g/kg) = " << Ice.x[ i ][ j ][ k ] << "   E_Ice (hPa) = " << E_Ice << "   sat_Deficit (hPa) = " << sat_Deficit << "   t_tau_h (°C) = " << t_tau_h << "   E_Rain_SL (hPa) = " << E_Rain_SL << endl << endl;

		cout << " t_SL (°C) = " << t_Celsius_SL << "   p_SL (hPa) = " << p_SL << "   a_SL (g/m3) = " << a_SL << "   c_SL (hPa/hPa) = " << c.x[ 0 ][ j ][ k ] << "   e_SL (hPa) = " << e_SL << "   q_Ice (hPa/hPa) = " << q_Ice << "   t_tau_SL (°C) = " << t_tau_SL << endl << endl;

		cout << " Evap_Haude (mm/d) = " << Evaporation_Haude << "   RF_e (%) = " << RF_e  << "   E_Rain_super (hPa) = " << E_Rain_super << "   q_Rain_super (hPa/hPa) = " << q_Rain_super << "   q_h (hPa/hPa) = " << q_h << "   q_SL (hPa/hPa) = " << q_SL << endl << endl;
	}
*/


// 0°C limit separates water vapour from condensate/evaporate und sublimate/vaporize
// water vapour turns to or developes from water or ice

// force and source terms presented one by one:

	Rain_aux = + WaterVapour * coeff_lv * ep * ( dcdr + dcdthe / rm + dcdphi / rmsinthe ); 									// only the latent heat of water vapour is used!!!!!!!
	RS_LatentHeat_Energy_Rain = Rain_aux;

//	Rain_super_aux = + WaterVapour * coeff_lv * ep * ( dcdr + dcdthe / rm + dcdphi / rmsinthe );
//	Ice_aux = + WaterVapour * coeff_ls * ep * ( dcdr + dcdthe / rm + dcdphi / rmsinthe );
/*
	if ( ( Rain_aux < 0. ) && ( q_h > q_Rain ) )	RS_LatentHeat_Energy_Rain = Rain_aux;						// only the latent heat of water vapour is used!!!!!!!
	else     RS_LatentHeat_Energy_Rain = 0.;
*/

//	if ( ( Rain_aux > 0. ) && ( t_Celsius >= 0. ) )	RS_LatentHeat_Energy_Rain = Rain_aux;
//	if ( ( Rain_super_aux < 0. ) && ( q_h > q_Rain_super ) && ( t_Celsius <= 0. ) && ( t_Celsius >= - 12. ) )	RS_LatentHeat_Energy_Rain_super = Rain_super_aux;
//	if ( ( Rain_super_aux > 0. ) && ( t_Celsius <= 0. ) )	RS_LatentHeat_Energy_Rain_super = Rain_super_aux;

//	if ( ( Ice_aux < 0. ) && ( q_h > q_Ice ) && ( t_Celsius < -12. ) && ( c.x[ i ][ j ][ k ] > 0. ) )	RS_LatentHeat_Energy_Ice = Ice_aux;
//	if ( ( Ice_aux > 0. ) && ( t_Celsius <= - 12. ) )	RS_LatentHeat_Energy_Ice = Ice_aux;


//	Latency.x[ i ][ j ][ k ] = RS_LatentHeat_Energy_Rain + RS_LatentHeat_Energy_Rain_super + RS_LatentHeat_Energy_Ice;
//	Latency.x[ i ][ j ][ k ] = RS_LatentHeat_Energy_Rain + RS_LatentHeat_Energy_Rain_super;
	Latency.x[ i ][ j ][ k ] = RS_LatentHeat_Energy_Rain;																				// only the latent heat of water vapour is used!!!!!!!


// Coriolis and centrifugal terms in the energy and momentum equations
	if ( h.x[ i ][ j ][ k ] == 0. )
	{
		RS_Coriolis_Energy = ( + u.x[ i ][ j ][ k ] * coriolis * 2. * omega * sinthe * w.x[ i ][ j ][ k ]
												- w.x[ i ][ j ][ k ] * coriolis * ( 2. * omega * sinthe * u.x[ i ][ j ][ k ] + 2. * omega * costhe * v.x[ i ][ j ][ k ] )
												+ v.x[ i ][ j ][ k ] * coriolis * 2. * omega * costhe * w.x[ i ][ j ][ k ] ) * ec * pr;

		RS_centrifugal_Energy = + centrifugal * rad.z[ i ] * pow ( ( omega * sinthe ), 2 ) * ec * pr
													+ centrifugal * rad.z[ i ] * sinthe * costhe * pow ( ( omega ), 2 ) * ec * pr;

		RS_Coriolis_Momentum_rad = + h_d_i * coriolis * 2. * omega * sinthe * w.x[ i ][ j ][ k ];
		RS_Coriolis_Momentum_the = + h_d_j * coriolis * 2. * omega * costhe * w.x[ i ][ j ][ k ];
		RS_Coriolis_Momentum_phi = - h_d_k * coriolis * ( 2. * omega * sinthe * u.x[ i ][ j ][ k ] + 2. * omega * costhe * v.x[ i ][ j ][ k ] );

		RS_centrifugal_Momentum_rad = + centrifugal * rad.z[ i ] * pow ( ( omega * sinthe ), 2 );
		RS_centrifugal_Momentum_the = + centrifugal * rad.z[ i ] * sinthe * costhe * pow ( ( omega ), 2 );
	}
	else
	{
		RS_Coriolis_Energy = 0.;
		RS_centrifugal_Energy = 0.;
		RS_Coriolis_Momentum_rad = 0.;
		RS_Coriolis_Momentum_the = 0.;
		RS_Coriolis_Momentum_phi = 0.;
		RS_centrifugal_Momentum_rad = 0.;
		RS_centrifugal_Momentum_the = 0.;
	}


// Right Hand Side of the time derivative ot temperature, pressure, water vapour concentration and velocity components

	rhs_t.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dtdr + v.x[ i ][ j ][ k ] * dtdthe / rm + w.x[ i ][ j ][ k ] * dtdphi / rmsinthe )
			+ ec * ( u.x[ i ][ j ][ k ] * dpdr + v.x[ i ][ j ][ k ] / rm * dpdthe + w.x[ i ][ j ][ k ] / rmsinthe * dpdphi )
			+ h_d_i * gr * u.x[ i ][ j ][ k ] * ( t.x[ i ][ j ][ k ] - t_Boussinesq ) / t_Boussinesq
			+ ( d2tdr2 + dtdr * 2. / rm + d2tdthe2 / rm2 + dtdthe * costhe / rm2sinthe + d2tdphi2 / rm2sinthe2 ) / ( re * pr )
			+ 2. * ec / re * ( ( dudr * dudr) + pow ( ( dvdthe / rm + h_d_i * u.x[ i ][ j ][ k ] / rm ), 2. )
			+ pow ( ( dwdphi / rmsinthe + h_d_i * u.x[ i ][ j ][ k ] / rm + h_d_i * v.x[ i ][ j ][ k ] * cotthe / rm ), 2. ) )
			+ ec / re * ( pow ( ( dvdr - h_d_i * v.x[ i ][ j ][ k ] / rm + dudthe / rm ), 2. )
			+ pow ( ( dudphi / rmsinthe + dwdr - h_d_i * w.x[ i ][ j ][ k ] / rm ), 2. )
			+ pow ( ( dwdthe * sinthe / rm2 - h_d_i * w.x[ i ][ j ][ k ] * costhe / rmsinthe + dvdphi / rmsinthe ), 2. ) )
			+ RS_Coriolis_Energy + RS_centrifugal_Energy
			- RS_LatentHeat_Energy_Rain; 																						// only the latent heat of water vapour is used!!!!!!!
//			- RS_LatentHeat_Energy_Rain - RS_LatentHeat_Energy_Rain_super - RS_LatentHeat_Energy_Ice;


	rhs_u.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dudr + v.x[ i ][ j ][ k ] * dudthe / rm + w.x[ i ][ j ][ k ] * dudphi / rmsinthe )
			- dpdr + h_d_i * gr * ( t.x[ i ][ j ][ k ] - t_Boussinesq ) / t_Boussinesq
			 + ( d2udr2 + h_d_i * 2. * u.x[ i ][ j ][ k ] / rm2 + d2udthe2 / rm2 + 4. * dudr / rm + dudthe * costhe / rm2sinthe + d2udphi2 / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_rad + RS_centrifugal_Momentum_rad
			- h_c_i * u.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_v.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dvdr + v.x[ i ][ j ][ k ] * dvdthe / rm + w.x[ i ][ j ][ k ] * dvdphi / rmsinthe ) +
			- dpdthe / rm + ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
			- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * h_d_j * v.x[ i ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2
			+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_the + RS_centrifugal_Momentum_the
			- h_c_j * v.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_w.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dwdr + v.x[ i ][ j ][ k ] * dwdthe / rm + w.x[ i ][ j ][ k ] * dwdphi / rmsinthe ) +
			- dpdphi / rmsinthe + ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
			- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * h_d_k * w.x[ i ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2
			+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_phi
			- h_c_k * w.x[ i ][ j ][ k ] * k_Force / dphi2;					// immersed boundary condition as a negative force addition


	rhs_c.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dcdr + v.x[ i ][ j ][ k ] * dcdthe / rm + w.x[ i ][ j ][ k ] * dcdphi / rmsinthe )
			+ ( d2cdr2 + dcdr * 2. / rm + d2cdthe2 / rm2 + dcdthe * costhe / rm2sinthe + d2cdphi2 / rm2sinthe2 ) / ( sc_WaterVapour * re )
			- ( t.x[ i ][ j ][ k ] - t_Boussinesq ) / t_Boussinesq
			- h_c_i * c.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_co2.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dco2dr + v.x[ i ][ j ][ k ] * dco2dthe / rm + w.x[ i ][ j ][ k ] * dco2dphi / rmsinthe )
			+ ( d2co2dr2 + dco2dr * 2. / rm + d2co2dthe2 / rm2 + dco2dthe * costhe / rm2sinthe + d2co2dphi2 / rm2sinthe2 ) / ( sc_CO2 * re )
			- h_c_i * co2.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


// for the Poisson equation to solve for the pressure, pressure gradient substracted from the above RHS

	aux_u.x[ i ][ j ][ k ] = rhs_u.x[ i ][ j ][ k ] + dpdr;
	aux_v.x[ i ][ j ][ k ] = rhs_v.x[ i ][ j ][ k ] + dpdthe / rm;
	aux_w.x[ i ][ j ][ k ] = rhs_w.x[ i ][ j ][ k ] + dpdphi / rmsinthe;
}






void RHS_Atmosphere::RK_RHS_2D_Atmosphere ( int j, int k, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &v, Array &w, Array &p, Array &vn, Array &wn, Array &rhs_v, Array &rhs_w, Array &aux_v, Array &aux_w )
{
//  2D surface iterations

//	k_Force = 1.;																			// factor for accelleration of convergence processes inside the immersed boundary conditions
	k_Force = 10.;																			// factor for accelleration of convergence processes inside the immersed boundary conditions

	cc = + 1.;

	h_check_i = h_check_j = h_check_k = 0;


// collection of coefficients
	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

	rm = rad.z[ 0 ];
	rm2 = rm * rm;

// collection of coefficients
	sinthe = sin( the.z[ j ] );
	sinthe2 = sinthe * sinthe;
	costhe = cos( the.z[ j ] );
	rmsinthe = rm * sinthe;
	rm2sinthe = rm2 * sinthe;
	rm2sinthe2 = rm2 * sinthe2;

/*
// in case needed for the influence of the radial derivatives 
//  3D volume iterations in case 1. and 2. order derivatives at walls are needed >>>>>>>>>>>>>>>>>>>>>>>>
// only in positive r-direction above ground 

	if ( ( h.x[ i -1 ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k ] == 0. ) )
	{
		h_0_i = ( double ) ( i ) * dr + dr / 4.;

		if ( fabs ( ( ( double ) ( i + 1 ) * dr - h_0_i ) ) < dr )
		{
			h_c_i = cc * ( 1. - fabs ( ( ( double ) ( i + 1 ) * dr - h_0_i ) ) / dr ); 
//			h_c_i = cc * ( .5 * ( acos ( fabs ( ( double ) ( i + 1 ) * dr - h_0_i ) * 3.14 / dr ) + 1. ) ); 
			h_d_i = 1. - h_c_i;
			h_check_i = 1;
		}
	}
*/



// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in positive the-direction along northerly boundaries 

	if ( ( h.x[ 0 ][ j - 1 ][ k ] == 1. ) && ( h.x[ 0 ][ j ][ k ] == 0. ) )
	{
		h_0_j = ( double ) ( j ) * dthe + dthe / 4.;

		if ( fabs ( ( ( double ) ( j + 1 ) * dthe - h_0_j ) ) < dthe )
		{
			h_c_j = cc * ( 1. - fabs ( ( ( double ) ( j + 1 ) * dthe - h_0_j ) ) / dthe ); 
//			h_c_j = cc * ( .5 * ( acos ( fabs ( ( double ) ( j + 1 ) * dthe - h_0_j ) * 3.14 / dthe ) + 1. ) ); 
			h_d_j = 1. - h_c_j;
			h_check_j = 1;
		}
	}



// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in negative the-direction along southerly boundaries 

	if ( ( h.x[ 0 ][ j + 1 ][ k ] == 1. ) && ( h.x[ 0 ][ j ][ k ] == 0. ) )
	{
		h_0_j = ( double ) ( j ) * dthe - dthe / 4.;

		if ( fabs ( ( ( double ) ( j - 1 ) * dthe - h_0_j ) ) < dthe )
		{
			h_c_j = cc * ( 1. - fabs ( ( ( double ) ( j - 1 ) * dthe - h_0_j ) ) / dthe ); 
//			h_c_j = cc * ( .5 * ( acos ( fabs ( ( double ) ( j - 1 ) * dthe - h_0_j ) * 3.14 / dthe ) + 1. ) ); 
			h_d_j = 1. - h_c_j;
			h_check_j = 1;
		}
	}


// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in positive phi-direction on westerly boundaries 

	if ( ( h.x[ 0 ][ j ][ k - 1 ] == 1. ) && ( h.x[ 0 ][ j ][ k ] == 0. ) )
	{
		h_0_k = ( double ) ( k ) * dphi + dphi / 4.;

		if ( fabs ( ( ( double ) ( k + 1 ) * dphi - h_0_k ) ) < dphi )
		{
			h_c_k = cc * ( 1. - fabs ( ( ( double ) ( k + 1 ) * dphi - h_0_k ) ) / dphi ); 
//			h_c_k = cc * ( .5 * ( acos ( fabs ( ( double ) ( k + 1 ) * dphi - h_0_k ) * 3.14 / dphi ) + 1. ) ); 
			h_d_k = 1. - h_c_k;
			h_check_k = 1;
		}
	}



// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in negative phi-direction along easterly boundaries 

	if ( ( h.x[ 0 ][ j ][ k + 1 ] == 1. ) && ( h.x[ 0 ][ j ][ k ] == 0. ) )
	{
		h_0_k = ( double ) ( k ) * dphi - dphi / 4.;

		if ( fabs ( ( ( double ) ( k - 1 ) * dphi - h_0_k ) ) < dphi )
		{
			h_c_k = cc * ( 1. - fabs ( ( ( double ) ( k - 1 ) * dphi - h_0_k ) ) / dphi ); 
//			h_c_k = cc * ( .5 * ( acos ( fabs ( ( double ) ( k - 1 ) * dphi - h_0_k ) * 3.14 / dphi ) + 1. ) ); 
			h_d_k = 1. - h_c_k;
			h_check_k = 1;
		}
	}


		if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( h_check_j != 1 ) )
		{
			h_c_j = 0.; 
			h_d_j = 1. - h_c_j;
		}

		if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h_check_j != 1 ) )
		{
			h_c_j = 1.; 
			h_d_j = 1. - h_c_j;
		}


		if ( ( h.x[ 0 ][ j ][ k ] == 0. ) && ( h_check_k != 1 ) )
		{
			h_c_k = 0.; 
			h_d_k = 1. - h_c_k;
		}

		if ( ( h.x[ 0 ][ j ][ k ] == 1. ) && ( h_check_k != 1 ) )
		{
			h_c_k = 1.; 
			h_d_k = 1. - h_c_k;
		}

	if ( h.x[ 0 ][ j ][ k ] == 0. )
	{
		dvdthe = h_d_j * ( v.x[ 0 ][ j+1 ][ k ] - v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );
		dwdthe = h_d_j * ( w.x[ 0 ][ j+1 ][ k ] - w.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );
		dpdthe = h_d_j * ( p.x[ 0 ][ j+1 ][ k ] - p.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );

		dvdphi = h_d_k * ( v.x[ 0 ][ j ][ k+1 ] - v.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
		dwdphi = h_d_k * ( w.x[ 0 ][ j ][ k+1 ] - w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
		dpdphi = h_d_k * ( p.x[ 0 ][ j ][ k+1 ] - p.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );

		d2vdthe2 = h_d_j *  ( v.x[ 0 ][ j+1 ][ k ] - 2. * v.x[ 0 ][ j ][ k ] + v.x[ 0 ][ j-1 ][ k ] ) / dthe2;
		d2wdthe2 = h_d_j * ( w.x[ 0 ][ j+1 ][ k ] - 2. * w.x[ 0 ][ j ][ k ] + w.x[ 0 ][ j-1 ][ k ] ) / dthe2;

		d2vdphi2 = h_d_k * ( v.x[ 0 ][ j ][ k+1 ] - 2. * v.x[ 0 ][ j ][ k ] + v.x[ 0 ][ j ][ k-1 ] ) / dphi2;
		d2wdphi2 = h_d_k * ( w.x[ 0 ][ j ][ k+1 ] - 2. * w.x[ 0 ][ j ][ k ] + w.x[ 0 ][ j ][ k-1 ] ) / dphi2;


		RS_Coriolis_Momentum_the = + coriolis * 2. * omega * costhe * w.x[ 0 ][ j ][ k ] * h_d_j;
		RS_centrifugal_Momentum_the = + centrifugal * rad.z[ 0 ] * sinthe * costhe * pow ( ( omega ), 2 );
		RS_Coriolis_Momentum_phi = - coriolis * omega * costhe * v.x[ 0 ][ j ][ k ] * h_d_k;

		rhs_v.x[ 0 ][ j ][ k ] = - ( v.x[ 0 ][ j ][ k ] * dvdthe / rm + w.x[ 0 ][ j ][ k ] * dvdphi / rmsinthe ) +
				- dpdthe / rm + ( d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
				- ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[ 0 ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
				+ RS_Coriolis_Momentum_the + RS_centrifugal_Momentum_the
				- h_c_j * v.x[ 0 ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition

		rhs_w.x[ 0 ][ j ][ k ] = - ( v.x[ 0 ][ j ][ k ] * dwdthe / rm +  w.x[ 0 ][ j ][ k ] * dwdphi / rmsinthe ) +
				- dpdphi / rmsinthe + ( d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
				- ( 1. + costhe * costhe / sinthe2 ) * h_d_k * w.x[ 0 ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2 + dvdphi * 2. * costhe / rm2sinthe2 ) / re
				+ RS_Coriolis_Momentum_phi
				- h_c_k * w.x[ 0 ][ j ][ k ] * k_Force / dphi2;					// immersed boundary condition as a negative force addition

		aux_v.x[ 0 ][ j ][ k ] = rhs_v.x[ 0 ][ j ][ k ] + dpdthe / rm;
		aux_w.x[ 0 ][ j ][ k ] = rhs_w.x[ 0 ][ j ][ k ] + dpdphi / rmsinthe;

	}

	else
	{
		dvdthe = 0.;
		dwdthe = 0;
		dpdthe = 0.;

		dvdphi = 0.;
		dwdphi = 0.;
		dpdphi = 0.;

		d2vdthe2 = 0.;
		d2wdthe2 = 0.;

		d2vdphi2 = 0.;
		d2wdphi2 = 0.;

		rhs_v.x[ 0 ][ j ][ k ] = 0.;
		rhs_w.x[ 0 ][ j ][ k ] = 0.;

		aux_v.x[ 0 ][ j ][ k ] = 0.;
		aux_w.x[ 0 ][ j ][ k ] = 0.;

		v.x[ 0 ][ j ][ k ] = 0.;
		w.x[ 0 ][ j ][ k ] = 0.;
	}
}
