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


RHS_Atmosphere::RHS_Atmosphere ( int jm, int km, double dthe, double dphi, double re, double omega, double coriolis, double centrifugal )
{
	this-> jm = jm;
	this-> km = km;
	this-> dthe = dthe;
	this-> dphi = dphi;
	this-> re = re;
	this-> omega = omega;
	this-> coriolis = coriolis;
	this-> centrifugal = centrifugal;
}


RHS_Atmosphere::RHS_Atmosphere ( int im, int jm, int km, double dt, double dr, double dthe, double dphi, double re, double ec, double sc_WaterVapour, double sc_CO2, double g, double pr, double omega, double coriolis, double centrifugal, double WaterVapour, double buoyancy, double CO2, double gam, double sigma, double lambda )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> dt = dt;
	this-> dr = dr;
	this-> dthe = dthe;
	this-> dphi = dphi;
	this-> re = re;
	this-> ec = ec;
	this-> sc_WaterVapour = sc_WaterVapour;
	this-> sc_CO2 = sc_CO2;
	this-> g = g;
	this-> pr = pr;
	this-> gam = gam;
	this-> lambda = lambda;
	this-> omega = omega;
	this-> coriolis = coriolis;
	this-> centrifugal = centrifugal;
	this-> WaterVapour = WaterVapour;
	this-> buoyancy = buoyancy;
	this-> CO2 = CO2;
	this-> sigma = sigma;


// array "im_tropopause" for configuring data due to latitude dependent tropopause

	im_tropopause = new int*[ jm ];

	for ( int l = 0; l < jm; l++ )
	{
		im_tropopause[ l ] = 0;
//		cout << im_tropopause[ l ] << endl;
	}
}


RHS_Atmosphere::~RHS_Atmosphere() 
{
	for ( int j = 0; j < jm; j++ )
	{
		delete [  ] im_tropopause[ j ];
	}

	delete [  ] im_tropopause;
}



void RHS_Atmosphere::RK_RHS_3D_Atmosphere ( int n, int i, int j, int k, double lv, double ls, double ep, double hp, double u_0, double t_0, double c_0, double co2_0, double p_0, double r_air, double r_water, double r_water_vapour, double r_co2, double L_atm, double cp_l, double R_Air, double R_WaterVapour, double R_co2, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &p_stat, Array &c, Array &cloud, Array &ice, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &cloudn, Array &icen, Array &co2n, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c, Array &rhs_cloud, Array &rhs_ice, Array &rhs_co2, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &IceLayer, Array &BuoyancyForce, Array &Q_Sensible, Array &P_rain, Array &P_snow, Array &S_v, Array &S_c, Array &S_i, Array &S_r, Array &S_s, Array &S_c_c )
{
// collection of coefficients for phase transformation
	coeff_lv = lv / ( cp_l * t_0 );					// coefficient for the specific latent vapourisation heat ( condensation heat ), coeff_lv = 9.1069 in [ / ]
	coeff_ls = ls / ( cp_l * t_0 );					// coefficient for the specific latent sublimation heat ( sublimation heat ) coeff_ls = 10.9031 in [ / ]
	coeff_lf = ( lv + ls ) / ( cp_l * t_0 );			// coefficient for the specific latent fusion heat ( sublimation heat ) coeff_ls = 20.01 in [ / ]
	coeff_L_atm_u_0 = L_atm / u_0;					// coefficient for the source terms in the ice-water transport equations

	c43 = 4. / 3.;
	c13 = 1. / 3.;

	k_Force = 10.;																			// factor for accelleration of convergence processes inside the immersed boundary conditions

	cc = + 1.;

	h_check_i = h_check_j = h_check_k = 0;


// 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )
// collection of coefficients9
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
	dpdr = h_d_i * ( p_dyn.x[ i+1 ][ j ][ k ] - p_dyn.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dcdr = h_d_i * ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dclouddr = h_d_i * ( cloud.x[ i+1 ][ j ][ k ] - cloud.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dicedr = h_d_i * ( ice.x[ i+1 ][ j ][ k ] - ice.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dcodr = h_d_i * ( co2.x[ i+1 ][ j ][ k ] - co2.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );

	dudthe = h_d_j * ( u.x[ i ][ j+1 ][ k ] - u.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dvdthe = h_d_j * ( v.x[ i ][ j+1 ][ k ] - v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dwdthe = h_d_j * ( w.x[ i ][ j+1 ][ k ] - w.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dtdthe = h_d_j * ( t.x[ i ][ j+1 ][ k ] - t.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dpdthe = h_d_j * ( p_dyn.x[ i ][ j+1 ][ k ] - p_dyn.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dcdthe = h_d_j * ( c.x[ i ][ j+1 ][ k ] - c.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dclouddthe = h_d_j * ( cloud.x[ i ][ j+1 ][ k ] - cloud.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dicedthe = h_d_j * ( ice.x[ i ][ j+1 ][ k ] - ice.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dcodthe = h_d_j * ( co2.x[ i ][ j+1 ][ k ] - co2.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );

	dudphi = h_d_k * ( u.x[ i ][ j ][ k+1 ] - u.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dvdphi = h_d_k * ( v.x[ i ][ j ][ k+1 ] - v.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dwdphi = h_d_k * ( w.x[ i ][ j ][ k+1 ] - w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dtdphi = h_d_k * ( t.x[ i ][ j ][ k+1 ] - t.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dpdphi = h_d_k * ( p_dyn.x[ i ][ j ][ k+1 ] - p_dyn.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dcdphi = h_d_k * ( c.x[ i ][ j ][ k+1 ] - c.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dclouddphi = h_d_k * ( cloud.x[ i ][ j ][ k+1 ] - cloud.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dicedphi = h_d_k * ( ice.x[ i ][ j ][ k+1 ] - ice.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dcodphi = h_d_k * ( co2.x[ i ][ j ][ k+1 ] - co2.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );


// 2nd order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
	d2udr2 = h_d_i * ( u.x[ i+1 ][ j ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] ) / dr2;
	d2vdr2 = h_d_i * ( v.x[ i+1 ][ j ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] ) / dr2;
	d2wdr2 = h_d_i * ( w.x[ i+1 ][ j ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] ) / dr2;
	d2tdr2 = h_d_i * ( t.x[ i+1 ][ j ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i-1 ][ j ][ k ] ) / dr2;
	d2cdr2 = h_d_i * ( c.x[ i+1 ][ j ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i-1 ][ j ][ k ] ) / dr2;
	d2clouddr2 = h_d_i * ( cloud.x[ i+1 ][ j ][ k ] - 2. * cloud.x[ i ][ j ][ k ] + cloud.x[ i-1 ][ j ][ k ] ) / dr2;
	d2icedr2 = h_d_i * ( ice.x[ i+1 ][ j ][ k ] - 2. * ice.x[ i ][ j ][ k ] + ice.x[ i-1 ][ j ][ k ] ) / dr2;
	d2codr2 = h_d_i * ( co2.x[ i+1 ][ j ][ k ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i-1 ][ j ][ k ] ) / dr2;

	d2udthe2 = h_d_j * ( u.x[ i ][ j+1 ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2vdthe2 = h_d_j * ( v.x[ i ][ j+1 ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2wdthe2 = h_d_j * ( w.x[ i ][ j+1 ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2tdthe2 = h_d_j * ( t.x[ i ][ j+1 ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2cdthe2 = h_d_j * ( c.x[ i ][ j+1 ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2clouddthe2 = h_d_j * ( cloud.x[ i ][ j+1 ][ k ] - 2. * cloud.x[ i ][ j ][ k ] + cloud.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2icedthe2 = h_d_j * ( ice.x[ i ][ j+1 ][ k ] - 2. * ice.x[ i ][ j ][ k ] + ice.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2codthe2 = h_d_j * ( co2.x[ i ][ j+1 ][ k ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i ][ j-1 ][ k ] ) / dthe2;

	d2udphi2 = h_d_k * ( u.x[ i ][ j ][ k+1 ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2vdphi2 = h_d_k * ( v.x[ i ][ j ][ k+1 ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2wdphi2 = h_d_k * ( w.x[ i ][ j ][ k+1 ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2tdphi2 = h_d_k * ( t.x[ i ][ j ][ k+1 ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2cdphi2 = h_d_k * ( c.x[ i ][ j ][ k+1 ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2clouddphi2 = h_d_k * ( cloud.x[ i ][ j ][ k+1 ] - 2. * cloud.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2icedphi2 = h_d_k * ( ice.x[ i ][ j ][ k+1 ] - 2. * ice.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2codphi2 = h_d_k * ( co2.x[ i ][ j ][ k+1 ] - 2. * co2.x[ i ][ j ][ k ] + co2.x[ i ][ j ][ k-1 ] ) / dphi2;


// Coriolis and centrifugal terms in the energy and momentum equations
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


// latent heat of water vapour
//	RS_LatentHeat_Energy_Rain = + Latency.x[ i ][ j ][ k ] * coeff_lv / ( r_air * lv * u_0 / L_atm );
//	Q_Sensible.x[ i ][ j ][ k ] = - cp_l * r_air * u_0 * t_0 / L_atm * ( u.x[ i ][ j ][ k ] * ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) + v.x[ i ][ j ][ k ] * ( t.x[ i ][ j+1 ][ k ] - t.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe ) / rm + w.x[ i ][ j ][ k ] * ( t.x[ i ][ j ][ k+1 ] - t.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi ) / rmsinthe );							// sensible heat in [W/m2] from energy transport equation

// Boussineq-approximation for the buoyancy force caused by humid air lighter than dry air
	r_dry = 100. * p_stat.x[ i ][ j ][ k ] / ( R_Air * t.x[ i ][ j ][ k ] * t_0 );														// density of dry air
	r_humid = r_dry / ( 1. + ( R_WaterVapour / R_Air - 1. ) * c.x[ i ][ j ][ k ] );												// density of humid air, COSMO version withot cloud and ice water, masses negligible

	RS_buoyancy_Momentum = - L_atm / ( u_0 * u_0 ) * buoyancy * g * ( r_humid / r_dry - 1. ); 					// any humid air is less dense than dry air

	BuoyancyForce.x[ i ][ j ][ k ] = RS_buoyancy_Momentum * u_0 * u_0 / L_atm * r_air;							// in N/m² = Pa
	if ( h.x[ i ][ j ][ k ] == 1. ) 	BuoyancyForce.x[ i ][ j ][ k ] = 0.;

	RS_buoyancy_Water_Vapour = buoyancy * ( cloud.x[ i ][ j ][ k ] - cloudn.x[ i ][ j ][ k ] );



				cout.precision ( 8 );
				cout.setf ( ios::fixed );
//				cout.precision ( 2 );
//				cout.setf ( std::ios::scientific );

//				if ( ( j == 90 ) && ( k == 180 ) )	cout << " RHS   "<< i << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  c = " << c.x[ i ][ j ][ k ] * 1000. << "  cl = " << cloud.x[ i ][ j ][ k ] * 1000. << "  ci = " << ice.x[ i ][ j ][ k ] * 1000. << "  q_Rain = " << q_Rain * 1000. << "  q_Ice = " << q_Ice * 1000. << "  S_i_dep = " << S_i_dep * 1000. << "  N_i = " << N_i << "  m_i = " << m_i * 1000000. << "  S_nuc = " << S_nuc * 1000. << "  S_c_frz = " << S_c_frz * 1000. << "  S_c_au = " << S_c_au * 1000. << "  S_i_au = " << S_i_au * 1000. << "  S_d_au = " << S_d_au * 1000. << "  S_ac = " << S_ac * 1000. << "  S_rim = " << S_rim * 1000. << "  S_shed = " << S_shed * 1000. << "  S_agg = " << S_agg * 1000. << "  S_i_cri = " << S_i_cri * 1000. << "  S_r_cri = " << S_r_cri * 1000. << "  S_ev = " << S_ev * 1000. << "  S_s_dep = " << S_s_dep * 1000. << "  S_i_melt = " << S_i_melt * 1000. << "  S_s_melt = " << S_s_melt * 1000. << "  S_r_frz = " << S_r_frz * 1000. << "  S_v = " << S_v.x[ i ][ j ][ k ] * 1000. << "  S_c = " << S_c.x[ i ][ j ][ k ] * 1000. << "  S_i = " << S_i.x[ i ][ j ][ k ] * 1000. << "  S_r = " << S_r.x[ i ][ j ][ k ] * 1000. << "  S_s = " << S_s.x[ i ][ j ][ k ] * 1000. << "  P_rain = " << P_rain.x[ i ][ j ][ k ] << "  P_snow = " << P_snow.x[ i ][ j ][ k ] << endl;



	S_c.x[ i ][ j ][ k ] = S_i.x[ i ][ j ][ k ] = S_r.x[ i ][ j ][ k ] = S_s.x[ i ][ j ][ k ] = 0.;




// Right Hand Side of the time derivative ot temperature, pressure, water vapour concentration and velocity components
	rhs_t.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dtdr + v.x[ i ][ j ][ k ] * dtdthe / rm + w.x[ i ][ j ][ k ] * dtdphi / rmsinthe )
			+ ec * .5 * ( u.x[ i ][ j ][ k ] * dpdr + v.x[ i ][ j ][ k ] / rm * dpdthe + w.x[ i ][ j ][ k ] / rmsinthe * dpdphi )
			+ ( d2tdr2 + dtdr * 2. / rm + d2tdthe2 / rm2 + dtdthe * costhe / rm2sinthe + d2tdphi2 / rm2sinthe2 ) / ( re * pr )
			+ 2. * ec / re * ( ( dudr * dudr) + pow ( ( dvdthe / rm + h_d_i * u.x[ i ][ j ][ k ] / rm ), 2. )
			+ pow ( ( dwdphi / rmsinthe + h_d_i * u.x[ i ][ j ][ k ] / rm + h_d_j * v.x[ i ][ j ][ k ] * cotthe / rm ), 2. ) )
			+ ec / re * ( pow ( ( dvdr - h_d_j * v.x[ i ][ j ][ k ] / rm + dudthe / rm ), 2. )
			+ pow ( ( dudphi / rmsinthe + dwdr - h_d_k * w.x[ i ][ j ][ k ] / rm ), 2. )
			+ pow ( ( dwdthe * sinthe / rm2 - h_d_k * w.x[ i ][ j ][ k ] * costhe / rmsinthe + dvdphi / rmsinthe ), 2. ) )
			+ RS_Coriolis_Energy + RS_centrifugal_Energy
			+ coeff_lv * coeff_L_atm_u_0 * ( S_c.x[ i ][ j ][ k ] + S_r.x[ i ][ j ][ k ] )
			+ coeff_ls * coeff_L_atm_u_0 * ( S_i.x[ i ][ j ][ k ] + S_s.x[ i ][ j ][ k ] )
			- h_c_i * t.x[ i ][ j ][ k ] * k_Force / dthe2;																					// immersed boundary condition as a negative force addition

	rhs_u.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dudr + v.x[ i ][ j ][ k ] * dudthe / rm + w.x[ i ][ j ][ k ] * dudphi / rmsinthe )
			- .5 * dpdr + ( d2udr2 + h_d_i * 2. * u.x[ i ][ j ][ k ] / rm2 + d2udthe2 / rm2 + 4. * dudr / rm + dudthe * costhe / rm2sinthe + d2udphi2 / rm2sinthe2 ) / re
			+ RS_buoyancy_Momentum
			+ RS_Coriolis_Momentum_rad + RS_centrifugal_Momentum_rad
			- h_c_i * u.x[ i ][ j ][ k ] * k_Force / dthe2;																					// immersed boundary condition as a negative force addition

	rhs_v.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dvdr + v.x[ i ][ j ][ k ] * dvdthe / rm + w.x[ i ][ j ][ k ] * dvdphi / rmsinthe ) +
			- .5 * dpdthe / rm + ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
			- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * h_d_j * v.x[ i ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2
			+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_the + RS_centrifugal_Momentum_the
			- h_c_j * v.x[ i ][ j ][ k ] * k_Force / dthe2;																					// immersed boundary condition as a negative force addition

	rhs_w.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dwdr + v.x[ i ][ j ][ k ] * dwdthe / rm + w.x[ i ][ j ][ k ] * dwdphi / rmsinthe ) +
			- .5 * dpdphi / rmsinthe + ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
			- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * h_d_k * w.x[ i ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2
			+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_phi
			- h_c_k * w.x[ i ][ j ][ k ] * k_Force / dphi2;																					// immersed boundary condition as a negative force addition

	rhs_c.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dcdr + v.x[ i ][ j ][ k ] * dcdthe / rm + w.x[ i ][ j ][ k ] * dcdphi / rmsinthe )
			+ ( d2cdr2 + dcdr * 2. / rm + d2cdthe2 / rm2 + dcdthe * costhe / rm2sinthe + d2cdphi2 / rm2sinthe2 ) / ( sc_WaterVapour * re )
			+ S_v.x[ i ][ j ][ k ] * coeff_L_atm_u_0
			- h_c_i * c.x[ i ][ j ][ k ] * k_Force / dphi2;

	rhs_cloud.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dclouddr + v.x[ i ][ j ][ k ] * dclouddthe / rm + w.x[ i ][ j ][ k ] * dclouddphi / rmsinthe )
			+ ( d2clouddr2 + dclouddr * 2. / rm + d2clouddthe2 / rm2 + dclouddthe * costhe / rm2sinthe + d2clouddphi2 / rm2sinthe2 ) / ( sc_WaterVapour * re )
			+ S_c.x[ i ][ j ][ k ] * coeff_L_atm_u_0
			- h_c_i * cloud.x[ i ][ j ][ k ] * k_Force / dphi2;

	rhs_ice.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dicedr + v.x[ i ][ j ][ k ] * dicedthe / rm + w.x[ i ][ j ][ k ] * dicedphi / rmsinthe )
			+ ( d2icedr2 + dicedr * 2. / rm + d2icedthe2 / rm2 + dicedthe * costhe / rm2sinthe + d2icedphi2 / rm2sinthe2 ) / ( sc_WaterVapour * re )
			+ S_i.x[ i ][ j ][ k ] * coeff_L_atm_u_0
			- h_c_i * ice.x[ i ][ j ][ k ] * k_Force / dphi2;

	rhs_co2.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dcodr + v.x[ i ][ j ][ k ] * dcodthe / rm + w.x[ i ][ j ][ k ] * dcodphi / rmsinthe )
			+ ( d2codr2 + dcodr * 2. / rm + d2codthe2 / rm2 + dcodthe * costhe / rm2sinthe + d2codphi2 / rm2sinthe2 ) / ( sc_CO2 * re )
			- h_c_i * co2.x[ i ][ j ][ k ] * k_Force / dthe2;						// immersed boundary condition as a negative force addition


// for the Poisson equation to solve for the pressure, pressure gradient substracted from the above RHS
	aux_u.x[ i ][ j ][ k ] = rhs_u.x[ i ][ j ][ k ] + dpdr;
	aux_v.x[ i ][ j ][ k ] = rhs_v.x[ i ][ j ][ k ] + dpdthe / rm;
	aux_w.x[ i ][ j ][ k ] = rhs_w.x[ i ][ j ][ k ] + dpdphi / rmsinthe;

}
















void RHS_Atmosphere::RK_RHS_2D_Atmosphere ( int j, int k, Array_1D &rad, Array_1D &the, Array &h, Array &v, Array &w, Array &p_dyn, Array &vn, Array &wn, Array &rhs_v, Array &rhs_w, Array &aux_v, Array &aux_w )
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
		dpdthe = h_d_j * ( p_dyn.x[ 0 ][ j+1 ][ k ] - p_dyn.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );

		dvdphi = h_d_k * ( v.x[ 0 ][ j ][ k+1 ] - v.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
		dwdphi = h_d_k * ( w.x[ 0 ][ j ][ k+1 ] - w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
		dpdphi = h_d_k * ( p_dyn.x[ 0 ][ j ][ k+1 ] - p_dyn.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );

		d2vdthe2 = h_d_j *  ( v.x[ 0 ][ j+1 ][ k ] - 2. * v.x[ 0 ][ j ][ k ] + v.x[ 0 ][ j-1 ][ k ] ) / dthe2;
		d2wdthe2 = h_d_j * ( w.x[ 0 ][ j+1 ][ k ] - 2. * w.x[ 0 ][ j ][ k ] + w.x[ 0 ][ j-1 ][ k ] ) / dthe2;

		d2vdphi2 = h_d_k * ( v.x[ 0 ][ j ][ k+1 ] - 2. * v.x[ 0 ][ j ][ k ] + v.x[ 0 ][ j ][ k-1 ] ) / dphi2;
		d2wdphi2 = h_d_k * ( w.x[ 0 ][ j ][ k+1 ] - 2. * w.x[ 0 ][ j ][ k ] + w.x[ 0 ][ j ][ k-1 ] ) / dphi2;


		RS_Coriolis_Momentum_the = + coriolis * 2. * omega * costhe * w.x[ 0 ][ j ][ k ] * h_d_j;
		RS_centrifugal_Momentum_the = + centrifugal * rad.z[ 0 ] * sinthe * costhe * pow ( ( omega ), 2 );
		RS_Coriolis_Momentum_phi = - coriolis * omega * costhe * v.x[ 0 ][ j ][ k ] * h_d_k;

		rhs_v.x[ 0 ][ j ][ k ] = - ( v.x[ 0 ][ j ][ k ] * dvdthe / rm + w.x[ 0 ][ j ][ k ] * dvdphi / rmsinthe ) +
				-  .5 *dpdthe / rm + ( d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
				- ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[ 0 ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
				+ RS_Coriolis_Momentum_the + RS_centrifugal_Momentum_the
				- h_c_j * v.x[ 0 ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition

		rhs_w.x[ 0 ][ j ][ k ] = - ( v.x[ 0 ][ j ][ k ] * dwdthe / rm +  w.x[ 0 ][ j ][ k ] * dwdphi / rmsinthe ) +
				-  .5 *dpdphi / rmsinthe + ( d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
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
