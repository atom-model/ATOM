/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
*/

#include <iostream>
#include <cmath>

#include "RHS_Hyd.h"

using namespace std;



RHS_Hydrosphere::RHS_Hydrosphere ( int jm, int km, double dthe, double dphi, double re, double omega, double coriolis, double centrifugal )
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


RHS_Hydrosphere::RHS_Hydrosphere ( int im, int jm, int km, double r0, double dt, double dr, double dthe, double dphi, double re, double ec, double sc, double g, double pr, double omega, double coriolis, double centrifugal, double buoyancy )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;
	this-> r0 = r0;
	this-> dt = dt;
	this-> dr = dr;
	this-> dthe = dthe;
	this-> dphi = dphi;
	this-> re = re;
	this-> ec = ec;
	this-> sc = sc;
	this-> g = g;
	this-> pr = pr;
	this-> omega = omega;
	this-> coriolis = coriolis;
	this-> centrifugal = centrifugal;
	this-> buoyancy = buoyancy;

}


RHS_Hydrosphere::~RHS_Hydrosphere() {}


void RHS_Hydrosphere::RK_RHS_3D_Hydrosphere ( int i, int j, int k, double L_hyd, double g, double cp_w, double u_0, double t_0, double c_0, double r_0_water, double ta, double pa, double ca, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &p_dynn, Array &cn, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_p, Array &rhs_c, Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion, Array &BuoyancyForce_3D, Array &Salt_Balance, Array &p_stat )
{
// collection of coefficients for phase transformation
	L_hyd = L_hyd / ( im - 1 );														// characteristic length for non-dimensionalisation

	c43 = 4. / 3.;
	c13 = 1. / 3.;


	k_Force = 10.;																			// factor for accelleration of convergence processes inside the immersed boundary conditions

	cc = + 1.;

	h_check_i = h_check_j = h_check_k = 0;


// 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )

// collection of coefficients
	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

// collection of coefficients
	rm = rad.z[ i ] / r0;
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
			h_d_j = 1. - h_c_j;
			h_check_j = 1;
		}
	}


// 3D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in negative the-direction along southerly boundaries 
	if ( ( h.x[ i ][ j +1 ][ k ] == 1. ) && ( h.x[ i ][ j ][ k ] == 0. ) )
	{
		h_0_j = ( double ) ( j ) * dthe - dthe / 4.;

		if ( fabs ( ( ( double ) ( j - 1 ) * dthe - h_0_j ) ) < dthe )
		{
			h_c_j = cc * ( 1. - fabs ( ( ( double ) ( j - 1 ) * dthe - h_0_j ) ) / dthe ); 
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



// 1. order derivative for temperature, pressure, salt concentrations and velocity components

// computation of initial and boundary conditions for the v and w velocity component
// computation at the surface

//  3D volume iterations

// 1st order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
	dudr = h_d_i * ( u.x[ i+1 ][ j ][ k ] - u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dvdr = h_d_i * ( v.x[ i+1 ][ j ][ k ] - v.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dwdr = h_d_i * ( w.x[ i+1 ][ j ][ k ] - w.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dtdr = h_d_i * ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dpdr = h_d_i * ( p_dyn.x[ i+1 ][ j ][ k ] - p_dyn.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
	dcdr = h_d_i * ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );

	dudthe = h_d_j * ( u.x[ i ][ j+1 ][ k ] - u.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dvdthe = h_d_j * ( v.x[ i ][ j+1 ][ k ] - v.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dwdthe = h_d_j * ( w.x[ i ][ j+1 ][ k ] - w.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dtdthe = h_d_j * ( t.x[ i ][ j+1 ][ k ] - t.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dpdthe = h_d_j * ( p_dyn.x[ i ][ j+1 ][ k ] - p_dyn.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );
	dcdthe = h_d_j * ( c.x[ i ][ j+1 ][ k ] - c.x[ i ][ j-1 ][ k ] ) / ( 2. * dthe );

	dudphi = h_d_k * ( u.x[ i ][ j ][ k+1 ] - u.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dvdphi = h_d_k * ( v.x[ i ][ j ][ k+1 ] - v.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dwdphi = h_d_k * ( w.x[ i ][ j ][ k+1 ] - w.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dtdphi = h_d_k * ( t.x[ i ][ j ][ k+1 ] - t.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dpdphi = h_d_k * ( p_dyn.x[ i ][ j ][ k+1 ] - p_dyn.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );
	dcdphi = h_d_k * ( c.x[ i ][ j ][ k+1 ] - c.x[ i ][ j ][ k-1 ] ) / ( 2. * dphi );

// 2nd order derivative for temperature, pressure, water vapour and co2 concentrations and velocity components
	d2udr2 = h_d_i * ( u.x[ i+1 ][ j ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] ) / dr2;
	d2vdr2 = h_d_i * ( v.x[ i+1 ][ j ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] ) / dr2;
	d2wdr2 = h_d_i * ( w.x[ i+1 ][ j ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] ) / dr2;
	d2tdr2 = h_d_i * ( t.x[ i+1 ][ j ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i-1 ][ j ][ k ] ) / dr2;
	d2cdr2 = h_d_i * ( c.x[ i+1 ][ j ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i-1 ][ j ][ k ] ) / dr2;

	d2udthe2 = h_d_j * ( u.x[ i ][ j+1 ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2vdthe2 = h_d_j * ( v.x[ i ][ j+1 ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2wdthe2 = h_d_j * ( w.x[ i ][ j+1 ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2tdthe2 = h_d_j * ( t.x[ i ][ j+1 ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j-1 ][ k ] ) / dthe2;
	d2cdthe2 = h_d_j * ( c.x[ i ][ j+1 ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j-1 ][ k ] ) / dthe2;

	d2udphi2 = h_d_k * ( u.x[ i ][ j ][ k+1 ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2vdphi2 = h_d_k * ( v.x[ i ][ j ][ k+1 ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2wdphi2 = h_d_k * ( w.x[ i ][ j ][ k+1 ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2tdphi2 = h_d_k * ( t.x[ i ][ j ][ k+1 ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i ][ j ][ k-1 ] ) / dphi2;
	d2cdphi2 = h_d_k * ( c.x[ i ][ j ][ k+1 ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i ][ j ][ k-1 ] ) / dphi2;



	if ( i < im - 2 )
	{
		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i + 1 ][ j ][ k ] == 0. ) )
		{
			dudr = ( - 3. * u.x[ i ][ j ][ k ] + 4. * u.x[ i + 1 ][ j ][ k ] - u.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );			// 2. order accurate

			d2udr2 = ( 2 * u.x[ i ][ j ][ k ] - 2. * u.x[ i + 1 ][ j ][ k ] + u.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr );			// 2. order accurate
		}
	}
	else
	{
		dudr = ( u.x[ i+1 ][ j ][ k ] - u.x[ i ][ j ][ k ] ) / dr;

		d2udr2 = 0.;
	}



	if ( ( j >= 2 ) && ( j < jm - 3 ) )
	{
		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( ( h.x[ i ][ j + 1 ][ k ] == 0. ) && ( h.x[ i ][ j + 2 ][ k ] == 0. ) ) )
		{
			dudthe = h_d_j * ( - 3. * u.x[ i ][ j ][ k ] + 4. * u.x[ i ][ j + 1 ][ k ] - u.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dvdthe = h_d_j * ( - 3. * v.x[ i ][ j ][ k ] + 4. * v.x[ i ][ j + 1 ][ k ] - v.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dwdthe = h_d_j * ( - 3. * w.x[ i ][ j ][ k ] + 4. * w.x[ i ][ j + 1 ][ k ] - w.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dpdthe = h_d_j * ( - 3. * p_dyn.x[ i ][ j ][ k ] + 4. * p_dyn.x[ i ][ j + 1 ][ k ] - p_dyn.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dtdthe = h_d_j * ( - 3. * t.x[ i ][ j ][ k ] + 4. * t.x[ i ][ j + 1 ][ k ] - t.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dcdthe = h_d_j * ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i ][ j + 1 ][ k ] - c.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate

			d2udthe2 = h_d_j * ( 2 * u.x[ i ][ j ][ k ] - 2. * u.x[ i ][ j + 1 ][ k ] + u.x[ i ][ j + 2 ][ k ] ) / dthe2;										// 2. order accurate
			d2vdthe2 = h_d_j * ( 2 * v.x[ i ][ j ][ k ] - 2. * v.x[ i ][ j + 1 ][ k ] + v.x[ i ][ j + 2 ][ k ] ) / dthe2;										// 2. order accurate
			d2wdthe2 = h_d_j * ( 2 * w.x[ i ][ j ][ k ] - 2. * w.x[ i ][ j + 1 ][ k ] + w.x[ i ][ j + 2 ][ k ] ) / dthe2;										// 2. order accurate
			d2tdthe2 = h_d_j * ( 2 * t.x[ i ][ j ][ k ] - 2. * t.x[ i ][ j + 1 ][ k ] + t.x[ i ][ j + 2 ][ k ] ) / dthe2;										// 2. order accurate
			d2cdthe2 = h_d_j * ( 2 * c.x[ i ][ j ][ k ] - 2. * c.x[ i ][ j + 1 ][ k ] + c.x[ i ][ j + 2 ][ k ] ) / dthe2;										// 2. order accurate
		}
		else
		{
			dudthe = h_d_j * ( u.x[ i ][ j + 1 ][ k ] - u.x[ i ][ j ][ k ] ) / dthe;
			dvdthe = h_d_j * ( v.x[ i ][ j + 1 ][ k ] - v.x[ i ][ j ][ k ] ) / dthe;
			dwdthe = h_d_j * ( w.x[ i ][ j + 1 ][ k ] - w.x[ i ][ j ][ k ] ) / dthe;
			dpdthe = h_d_j * ( p_dyn.x[ i ][ j + 1 ][ k ] - p_dyn.x[ i ][ j ][ k ] ) / dthe;
			dtdthe = h_d_j * ( t.x[ i ][ j + 1 ][ k ] - t.x[ i ][ j ][ k ] ) / dthe;
			dcdthe = h_d_j * ( c.x[ i ][ j + 1 ][ k ] - c.x[ i ][ j ][ k ] ) / dthe;

			d2udthe2 = d2vdthe2 = d2wdthe2 = d2tdthe2 = d2cdthe2 = 0.;
		}


		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) && ( h.x[ i ][ j - 2 ][ k ] == 0. ) )
		{
			dudthe = h_d_j * ( - 3. * u.x[ i ][ j ][ k ] + 4. * u.x[ i ][ j - 1 ][ k ] - u.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dvdthe = h_d_j * ( - 3. * v.x[ i ][ j ][ k ] + 4. * v.x[ i ][ j - 1 ][ k ] - v.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dwdthe = h_d_j * ( - 3. * w.x[ i ][ j ][ k ] + 4. * w.x[ i ][ j - 1 ][ k ] - w.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dpdthe = h_d_j * ( - 3. * p_dyn.x[ i ][ j ][ k ] + 4. * p_dyn.x[ i ][ j - 1 ][ k ] - p_dyn.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dtdthe = h_d_j * ( - 3. * t.x[ i ][ j ][ k ] + 4. * t.x[ i ][ j - 1 ][ k ] - t.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dcdthe = h_d_j * ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i ][ j - 1 ][ k ] - c.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate

			d2udthe2 = h_d_j * ( 2 * u.x[ i ][ j ][ k ] - 2. * u.x[ i ][ j - 1 ][ k ] + u.x[ i ][ j - 2 ][ k ] ) / dthe2;										// 2. order accurate
			d2vdthe2 = h_d_j * ( 2 * v.x[ i ][ j ][ k ] - 2. * v.x[ i ][ j - 1 ][ k ] + v.x[ i ][ j - 2 ][ k ] ) / dthe2;										// 2. order accurate
			d2wdthe2 = h_d_j * ( 2 * w.x[ i ][ j ][ k ] - 2. * w.x[ i ][ j - 1 ][ k ] + w.x[ i ][ j - 2 ][ k ] ) / dthe2;										// 2. order accurate
			d2tdthe2 = h_d_j * ( 2 * t.x[ i ][ j ][ k ] - 2. * t.x[ i ][ j - 1 ][ k ] + t.x[ i ][ j - 2 ][ k ] ) / dthe2;										// 2. order accurate
			d2cdthe2 = h_d_j * ( 2 * c.x[ i ][ j ][ k ] - 2. * c.x[ i ][ j - 1 ][ k ] + c.x[ i ][ j - 2 ][ k ] ) / dthe2;										// 2. order accurate
		}
		else
		{
			dudthe = h_d_j * ( u.x[ i ][ j ][ k ] - u.x[ i ][ j - 1 ][ k ] ) / dthe;
			dvdthe = h_d_j * ( v.x[ i ][ j ][ k ] - v.x[ i ][ j - 1 ][ k ] ) / dthe;
			dwdthe = h_d_j * ( w.x[ i ][ j ][ k ] - w.x[ i ][ j - 1 ][ k ] ) / dthe;
			dpdthe = h_d_j * ( p_dyn.x[ i ][ j ][ k ] - p_dyn.x[ i ][ j - 1 ][ k ] ) / dthe;
			dtdthe = h_d_j * ( t.x[ i ][ j ][ k ] - t.x[ i ][ j - 1 ][ k ] ) / dthe;
			dcdthe = h_d_j * ( c.x[ i ][ j ][ k ] - c.x[ i ][ j - 1 ][ k ] ) / dthe;

			d2udthe2 = d2vdthe2 = d2wdthe2 = d2tdthe2 = d2cdthe2 = 0.;
		}
	}
	else
	{
		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j + 1 ][ k ] == 0. ) )
		{
			dudthe = h_d_j * ( u.x[ i ][ j + 1 ][ k ] - u.x[ i ][ j ][ k ] ) / dthe;
			dvdthe = h_d_j * ( v.x[ i ][ j + 1 ][ k ] - v.x[ i ][ j ][ k ] ) / dthe;
			dwdthe = h_d_j * ( w.x[ i ][ j + 1 ][ k ] - w.x[ i ][ j ][ k ] ) / dthe;
			dpdthe = h_d_j * ( p_dyn.x[ i ][ j + 1 ][ k ] - p_dyn.x[ i ][ j ][ k ] ) / dthe;
			dtdthe = h_d_j * ( t.x[ i ][ j + 1 ][ k ] - t.x[ i ][ j ][ k ] ) / dthe;
			dcdthe = h_d_j * ( c.x[ i ][ j + 1 ][ k ] - c.x[ i ][ j ][ k ] ) / dthe;
		}

		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j - 1 ][ k ] == 0. ) )
		{
			dudthe = h_d_j * ( u.x[ i ][ j ][ k ] - u.x[ i ][ j - 1 ][ k ] ) / dthe;
			dvdthe = h_d_j * ( v.x[ i ][ j ][ k ] - v.x[ i ][ j - 1 ][ k ] ) / dthe;
			dwdthe = h_d_j * ( w.x[ i ][ j ][ k ] - w.x[ i ][ j - 1 ][ k ] ) / dthe;
			dpdthe = h_d_j * ( p_dyn.x[ i ][ j ][ k ] - p_dyn.x[ i ][ j - 1 ][ k ] ) / dthe;
			dtdthe = h_d_j * ( t.x[ i ][ j ][ k ] - t.x[ i ][ j - 1 ][ k ] ) / dthe;
			dcdthe = h_d_j * ( c.x[ i ][ j ][ k ] - c.x[ i ][ j - 1 ][ k ] ) / dthe;
		}
		d2udthe2 = d2vdthe2 = d2wdthe2 = d2tdthe2 = d2cdthe2 = 0.;
	}


	if ( ( k >= 2 ) && ( k < km - 3 ) )
	{
		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) && ( h.x[ i ][ j ][ k + 2 ] == 0. ) )
		{
			dudphi = h_d_k * ( - 3. * u.x[ i ][ j ][ k ] + 4. * u.x[ i ][ j ][ k + 1 ] - u.x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dvdphi = h_d_k * ( - 3. * v.x[ i ][ j ][ k ] + 4. * v.x[ i ][ j ][ k + 1 ] - v.x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dwdphi = h_d_k * ( - 3. * w.x[ i ][ j ][ k ] + 4. * w.x[ i ][ j ][ k + 1 ] - w.x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dpdphi = h_d_k * ( - 3. * p_dyn.x[ i ][ j ][ k ] + 4. * p_dyn.x[ i ][ j ][ k + 1 ] - p_dyn.x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dtdphi = h_d_k * ( - 3. * t.x[ i ][ j ][ k ] + 4. * t.x[ i ][ j ][ k + 1 ] - t.x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dcdphi = h_d_k * ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i ][ j ][ k + 1 ] - c.x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate

			d2udthe2 = h_d_k * ( 2 * u.x[ i ][ j ][ k ] - 2. * u.x[ i ][ j ][ k + 1 ] + u.x[ i ][ j ][ k + 2 ] ) / dphi2;										// 2. order accurate
			d2vdthe2 = h_d_k * ( 2 * v.x[ i ][ j ][ k ] - 2. * v.x[ i ][ j ][ k + 1 ] + v.x[ i ][ j ][ k + 2 ] ) / dphi2;										// 2. order accurate
			d2wdthe2 = h_d_k * ( 2 * w.x[ i ][ j ][ k ] - 2. * w.x[ i ][ j ][ k + 1 ] + w.x[ i ][ j ][ k + 2 ] ) / dphi2;										// 2. order accurate
			d2tdthe2 = h_d_k * ( 2 * t.x[ i ][ j ][ k ] - 2. * t.x[ i ][ j ][ k + 1 ] + t.x[ i ][ j ][ k + 2 ] ) / dphi2;										// 2. order accurate
			d2cdthe2 = h_d_k * ( 2 * c.x[ i ][ j ][ k ] - 2. * c.x[ i ][ j ][ k + 1 ] + c.x[ i ][ j ][ k + 2 ] ) / dphi2;										// 2. order accurate
		}
		else
		{
			dudphi = h_d_k * ( u.x[ i ][ j ][ k + 1 ] - u.x[ i ][ j ][ k ] ) / dphi;
			dvdphi = h_d_k * ( v.x[ i ][ j ][ k + 1 ] - v.x[ i ][ j ][ k ] ) / dphi;
			dwdphi = h_d_k * ( w.x[ i ][ j ][ k + 1 ] - w.x[ i ][ j ][ k ] ) / dphi;
			dpdphi = h_d_k * ( p_dyn.x[ i ][ j ][ k + 1 ] - p_dyn.x[ i ][ j ][ k ] ) / dphi;
			dtdphi = h_d_k * ( t.x[ i ][ j ][ k + 1 ] - t.x[ i ][ j ][ k ] ) / dphi;
			dcdphi = h_d_k * ( c.x[ i ][ j ][ k + 1 ] - c.x[ i ][ j ][ k ] ) / dphi;

			d2udphi2 = d2vdphi2 = d2wdphi2 = d2tdphi2 = d2cdphi2 = 0.;
		}


		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k - 1 ] == 0. ) && ( h.x[ i ][ j ][ k - 2 ] == 0. ) )
		{
			dudphi = h_d_k * ( - 3. * u.x[ i ][ j ][ k ] + 4. * u.x[ i ][ j ][ k - 1 ] - u.x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dvdphi = h_d_k * ( - 3. * v.x[ i ][ j ][ k ] + 4. * v.x[ i ][ j ][ k - 1 ] - v.x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dwdphi = h_d_k * ( - 3. * w.x[ i ][ j ][ k ] + 4. * w.x[ i ][ j ][ k - 1 ] - w.x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dpdphi = h_d_k * ( - 3. * p_dyn.x[ i ][ j ][ k ] + 4. * p_dyn.x[ i ][ j ][ k - 1 ] - p_dyn.x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dtdphi = h_d_k * ( - 3. * t.x[ i ][ j ][ k ] + 4. * t.x[ i ][ j ][ k - 1 ] - t.x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dcdphi = h_d_k * ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i ][ j ][ k - 1 ] - c.x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate

			d2udthe2 = h_d_k * ( 2 * u.x[ i ][ j ][ k ] - 2. * u.x[ i ][ j ][ k - 1 ] + u.x[ i ][ j ][ k - 2 ] ) / dphi2;										// 2. order accurate
			d2vdthe2 = h_d_k * ( 2 * v.x[ i ][ j ][ k ] - 2. * v.x[ i ][ j ][ k - 1 ] + v.x[ i ][ j ][ k - 2 ] ) / dphi2;										// 2. order accurate
			d2wdthe2 = h_d_k * ( 2 * w.x[ i ][ j ][ k ] - 2. * w.x[ i ][ j ][ k - 1 ] + w.x[ i ][ j ][ k - 2 ] ) / dphi2;										// 2. order accurate
			d2tdthe2 = h_d_k * ( 2 * t.x[ i ][ j ][ k ] - 2. * t.x[ i ][ j ][ k - 1 ] + t.x[ i ][ j ][ k - 2 ] ) / dphi2;										// 2. order accurate
			d2cdthe2 = h_d_k * ( 2 * c.x[ i ][ j ][ k ] - 2. * c.x[ i ][ j ][ k - 1 ] + c.x[ i ][ j ][ k - 2 ] ) / dphi2;										// 2. order accurate
		}
		else
		{
			dudphi = h_d_k * ( u.x[ i ][ j ][ k ] - u.x[ i ][ j ][ k - 1 ] ) / dphi;
			dvdphi = h_d_k * ( v.x[ i ][ j ][ k ] - v.x[ i ][ j ][ k - 1 ] ) / dphi;
			dwdphi = h_d_k * ( w.x[ i ][ j ][ k ] - w.x[ i ][ j ][ k - 1 ] ) / dphi;
			dpdphi = h_d_k * ( p_dyn.x[ i ][ j ][ k ] - p_dyn.x[ i ][ j ][ k - 1 ] ) / dphi;
			dtdphi = h_d_k * ( t.x[ i ][ j ][ k ] - t.x[ i ][ j ][ k - 1 ] ) / dphi;
			dcdphi = h_d_k * ( c.x[ i ][ j ][ k ] - c.x[ i ][ j ][ k - 1 ] ) / dphi;

			d2udphi2 = d2vdphi2 = d2wdphi2 = d2tdphi2 = d2cdphi2 = 0.;
		}
	}
	else
	{
		if ( ( h.x[ i ][ j ][ k ] == 1. ) && ( h.x[ i ][ j ][ k + 1 ] == 0. ) )
		{
			dudphi = h_d_k * ( u.x[ i ][ j ][ k + 1 ] - u.x[ i ][ j ][ k ] ) / dphi;
			dvdphi = h_d_k * ( v.x[ i ][ j ][ k + 1 ] - v.x[ i ][ j ][ k ] ) / dphi;
			dwdphi = h_d_k * ( w.x[ i ][ j ][ k + 1 ] - w.x[ i ][ j ][ k ] ) / dphi;
			dpdphi = h_d_k * ( p_dyn.x[ i ][ j ][ k + 1 ] - p_dyn.x[ i ][ j ][ k ] ) / dphi;
			dtdphi = h_d_k * ( t.x[ i ][ j ][ k + 1 ] - t.x[ i ][ j ][ k ] ) / dphi;
			dcdphi = h_d_k * ( c.x[ i ][ j ][ k + 1 ] - c.x[ i ][ j ][ k ] ) / dphi;
		}

		if ( ( h.x[ i ][ j ][ k ] == 0. ) && ( h.x[ i ][ j ][ k - 1 ] == 1. ) )
		{
			dudphi = h_d_k * ( u.x[ i ][ j ][ k ] - u.x[ i ][ j ][ k - 1 ] ) / dphi;
			dvdphi = h_d_k * ( v.x[ i ][ j ][ k ] - v.x[ i ][ j ][ k - 1 ] ) / dphi;
			dwdphi = h_d_k * ( w.x[ i ][ j ][ k ] - w.x[ i ][ j ][ k - 1 ] ) / dphi;
			dpdphi = h_d_k * ( p_dyn.x[ i ][ j ][ k ] - p_dyn.x[ i ][ j ][ k - 1 ] ) / dphi;
			dtdphi = h_d_k * ( t.x[ i ][ j ][ k ] - t.x[ i ][ j ][ k - 1 ] ) / dphi;
			dcdphi = h_d_k * ( c.x[ i ][ j ][ k ] - c.x[ i ][ j ][ k - 1 ] ) / dphi;
		}
		d2udphi2 = d2vdphi2 = d2wdphi2 = d2tdphi2 = d2cdphi2 = 0.;
	}





// coriolis and centrifugal forces are due to the very small rotation number of the earth extremely small
// their effect will be imagined after a huge number of iterations
	RS_Coriolis_Energy = ( + u.x[ i ][ j ][ k ] * coriolis * 2. * omega * sinthe * w.x[ i ][ j ][ k ]
							- w.x[ i ][ j ][ k ] * coriolis * ( 2. * omega * sinthe * u.x[ i ][ j ][ k ] + 2. * omega * costhe * v.x[ i ][ j ][ k ] )
							+ v.x[ i ][ j ][ k ] * coriolis * 2. * omega * costhe * w.x[ i ][ j ][ k ] ) * ec * pr;

	RS_Centrifugal_Energy = + centrifugal * rad.z[ i ] * pow ( ( omega * sinthe ), 2 ) * ec * pr
							 + centrifugal * rad.z[ i ] * sinthe * costhe * pow ( ( omega ), 2 ) * ec * pr;

	RS_Coriolis_Momentum_rad = +  h_d_i * coriolis * 2. * omega * sinthe * w.x[ i ][ j ][ k ];
	RS_Coriolis_Momentum_the = +  h_d_j * coriolis * 2. * omega * costhe * w.x[ i ][ j ][ k ];
	RS_Coriolis_Momentum_phi = -  h_d_k * coriolis * ( 2. * omega * sinthe * u.x[ i ][ j ][ k ] + 2. * omega * costhe * v.x[ i ][ j ][ k ] );

	RS_Centrifugal_Momentum_rad = + centrifugal * rad.z[ i ] * pow ( ( omega * sinthe ), 2 );
	RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ i ] * sinthe * costhe * pow ( ( omega ), 2 );

	RS_buoyancy_Momentum = L_hyd / ( r_0_water * u_0 * u_0 ) * buoyancy * g * ( ( c.x[ i ][ j ][ k ] * ca - 7.3 ) / ( ca - 7.3 ) - 1. );		// buoyancy based on water density 

	RS_buoyancy_Energy = u_0 * ec * u.x[ i ][ j ][ k ] * RS_buoyancy_Momentum;																							// energy increase by buoyancy

	BuoyancyForce_3D.x[ i ][ j ][ k ] = RS_buoyancy_Momentum * ( r_0_water * u_0 * u_0 ) / L_hyd;																// dimension as pressure in N/m2

 	Salt_Balance.x[ i ][ j ][ k ] = ( c.x[ i ][ j ][ k ] - 1. ) * c_0;																																// difference of salinity compared to average


	if ( Salt_Balance.x[ i ][ j ][ k ] >= 0. )
	{
		Salt_Finger.x[ i ][ j ][ k ] = Salt_Balance.x[ i ][ j ][ k ];		// for positiv salinity balance
	}
	else
	{
		Salt_Finger.x[ i ][ j ][ k ] = 0.;
	}

	if ( Salt_Balance.x[ i ][ j ][ k ] < 0. )
	{
		Salt_Diffusion.x[ i ][ j ][ k ] = Salt_Balance.x[ i ][ j ][ k ];	// for negativ salinity balance
	}
	else
	{
		Salt_Diffusion.x[ i ][ j ][ k ] = 0.;
	}


	if ( h.x[ i ][ j ][ k ] == 1. )
	{
		Salt_Balance.x[ i ][ j ][ k ] = 0.;
		Salt_Finger.x[ i ][ j ][ k ] = 0.;
		BuoyancyForce_3D.x[ i ][ j ][ k ] = 0.;

		RS_Salt_Energy = 0.;

		RS_Coriolis_Energy = 0.;
		RS_Coriolis_Momentum_rad = 0.;
		RS_Coriolis_Momentum_the = 0.;
		RS_Coriolis_Momentum_phi = 0.;

		RS_Centrifugal_Energy = 0.;
		RS_Centrifugal_Momentum_rad = 0.;
		RS_Centrifugal_Momentum_the = 0.;
	}


// Right Hand Side of the time derivative ot temperature, pressure, salt concentration and velocity components

//  3D volume iterations
	rhs_t.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dtdr + v.x[ i ][ j ][ k ] * dtdthe / rm + w.x[ i ][ j ][ k ] * dtdphi / rmsinthe )
			- ec * ( u.x[ i ][ j ][ k ] * dpdr + v.x[ i ][ j ][ k ] / rm * dpdthe + w.x[ i ][ j ][ k ] / rmsinthe * dpdphi )
			+ ( d2tdr2 + dtdr * 2. / rm + d2tdthe2 / rm2 + dtdthe * costhe / rm2sinthe + d2tdphi2 / rm2sinthe2 ) / ( re * pr )
			+ 2. * ec / re * ( ( dudr * dudr) + pow ( ( dvdthe / rm + h_d_i * u.x[ i ][ j ][ k ] / rm ), 2. )
			+ pow ( ( dwdphi / rmsinthe + h_d_i * u.x[ i ][ j ][ k ] / rm + h_d_j * v.x[ i ][ j ][ k ] * cotthe / rm ), 2. ) )
			+ ec / re * ( pow ( ( dvdr - h_d_j * v.x[ i ][ j ][ k ] / rm + dudthe / rm ), 2. ) + pow ( ( dudphi / rmsinthe + dwdr - h_d_k * w.x[ i ][ j ][ k ] / rm ), 2. )
			+ pow ( ( dwdthe * sinthe / rm2 - h_d_k * w.x[ i ][ j ][ k ] * costhe / rmsinthe + dvdphi / rmsinthe ), 2. ) )
			+ RS_buoyancy_Energy
			+ RS_Coriolis_Energy + RS_Centrifugal_Energy;
//			- h_c_i * t.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_u.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dudr + v.x[ i ][ j ][ k ] * dudthe / rm + w.x[ i ][ j ][ k ] * dudphi / rmsinthe )
			- dpdr + ( d2udr2 + h_d_i * 2. * u.x[ i ][ j ][ k ] / rm2 + d2udthe2 / rm2 + 4. * dudr / rm + dudthe * costhe / rm2sinthe + d2udphi2 / rm2sinthe2 ) / re
			+ RS_buoyancy_Momentum
			+ RS_Coriolis_Momentum_rad + RS_Centrifugal_Momentum_rad
			- h_c_i * u.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_v.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dvdr + v.x[ i ][ j ][ k ] * dvdthe / rm + w.x[ i ][ j ][ k ] * dvdphi / rmsinthe )
			- dpdthe / rm + ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
			- ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[ i ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2
			+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_the + RS_Centrifugal_Momentum_the
			- h_c_j * v.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_w.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dwdr + v.x[ i ][ j ][ k ] * dwdthe / rm + w.x[ i ][ j ][ k ] * dwdphi / rmsinthe )
			- dpdphi / rmsinthe + ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
			- ( 1. + costhe * costhe / sinthe2 ) * h_d_k * w.x[ i ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2
			+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_phi
			- h_c_k * w.x[ i ][ j ][ k ] * k_Force / dphi2;					// immersed boundary condition as a negative force addition

	rhs_p.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] *  ( rhs_u.x[ i ][ j ][ k ] + dpdr ) + v.x[ i ][ j ][ k ] * ( rhs_v.x[ i ][ j ][ k ] + dpdthe / rm ) + w.x[ i ][ j ][ k ] * ( rhs_w.x[ i ][ j ][ k ] + dpdphi / rmsinthe ) )
				- h_c_k * p_dyn.x[ i ][ j ][ k ] * k_Force / dphi2;					// immersed boundary condition as a negative force addition


	rhs_c.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dcdr + v.x[ i ][ j ][ k ] * dcdthe / rm + w.x[ i ][ j ][ k ] * dcdphi / rmsinthe )
			+ ( d2cdr2 + dcdr * 2. / rm + d2cdthe2 / rm2 + dcdthe * costhe / rm2sinthe + d2cdphi2 / rm2sinthe2 ) / ( sc * re );
//			- h_c_i * c.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	// for the Poisson equation to solve for the pressure, pressure gradient sbstracted from the RHS

	aux_u.x[ i ][ j ][ k ] = rhs_u.x[ i ][ j ][ k ] + dpdr;
	aux_v.x[ i ][ j ][ k ] = rhs_v.x[ i ][ j ][ k ] + dpdthe / rm;
	aux_w.x[ i ][ j ][ k ] = rhs_w.x[ i ][ j ][ k ] + dpdphi / rmsinthe;

/*
	if ( u.x[ i ][ j ][ k ] >= .0002 )			u.x[ i ][ j ][ k ] = .0002;
	if ( u.x[ i ][ j ][ k ] <= - .0002 )			u.x[ i ][ j ][ k] = - .0002;

	if ( v.x[ i ][ j ][ k ] >= .0024 )					v.x[ i ][ j ][ k ] = .0024;
	if ( v.x[ i ][ j ][ k ] <= - .0024 )				v.x[ i ][ j ][ k ] = - .0024;

	if ( w.x[ i ][ j ][ k ] >= .054 )				w.x[ i ][ j ][ k ] = .054;
	if ( w.x[ i ][ j ][ k ] <= - 054 )				w.x[ i ][ j ][ k ] = - .054;
*/
}




void RHS_Hydrosphere::RK_RHS_2D_Hydrosphere ( int j, int k, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &v, Array &w, Array &p_dyn, Array &vn, Array &wn, Array &p_dynn, Array &rhs_v, Array &rhs_w, Array &rhs_p, Array &aux_v, Array &aux_w )
{
//  2D surface iterations
	im = 41;
	k_Force = 10.;																			// factor for accelleration of convergence processes inside the immersed boundary conditions

	cc = + 1.;

	h_check_i = h_check_j = h_check_k = 0;

// collection of coefficients
	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

	rm = rad.z[ im-1 ];
	rm2 = rm * rm;

// collection of coefficients
	sinthe = sin( the.z[ j ] );
	sinthe2 = sinthe * sinthe;
	costhe = cos( the.z[ j ] );
	rmsinthe = rm * sinthe;
	rm2sinthe = rm2 * sinthe;
	rm2sinthe2 = rm2 * sinthe2;

// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in positive the-direction along northerly boundaries 
	if ( ( h.x[ im-1 ][ j - 1 ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k ] == 0. ) )
	{
		h_0_j = ( double ) ( j ) * dthe + dthe / 4.;

		if ( fabs ( ( ( double ) ( j + 1 ) * dthe - h_0_j ) ) < dthe )
		{
			h_c_j = 0.; 
			h_d_j = 1. - h_c_j;
			h_check_j = 1;
		}
	}


// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in negative the-direction along southerly boundaries 
	if ( ( h.x[ im-1 ][ j + 1 ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k ] == 0. ) )
	{
		h_0_j = ( double ) ( j ) * dthe - dthe / 4.;

		if ( fabs ( ( ( double ) ( j - 1 ) * dthe - h_0_j ) ) < dthe )
		{
			h_c_j = 0.; 
			h_d_j = 1. - h_c_j;
			h_check_j = 1;
		}
	}


// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in positive phi-direction on westerly boundaries 
	if ( ( h.x[ im-1 ][ j ][ k - 1 ] == 1. ) && ( h.x[ im-1 ][ j ][ k ] == 0. ) )
	{
		h_0_k = ( double ) ( k ) * dphi + dphi / 4.;

		if ( fabs ( ( ( double ) ( k + 1 ) * dphi - h_0_k ) ) < dphi )
		{
			h_c_k = 0.; 
			h_d_k = 1. - h_c_k;
			h_check_k = 1;
		}
	}


// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in negative phi-direction along easterly boundaries 
	if ( ( h.x[ im-1 ][ j ][ k + 1 ] == 1. ) && ( h.x[ im-1 ][ j ][ k ] == 0. ) )
	{
		h_0_k = ( double ) ( k ) * dphi - dphi / 4.;

		if ( fabs ( ( ( double ) ( k - 1 ) * dphi - h_0_k ) ) < dphi )
		{
			h_c_k = 0.; 
			h_d_k = 1. - h_c_k;
			h_check_k = 1;
		}
	}


	if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h_check_j != 1 ) )
	{
		h_c_j = 0.; 
		h_d_j = 1. - h_c_j;
		}

	if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h_check_j != 1 ) )
	{
		h_c_j = 1.; 
		h_d_j = 1. - h_c_j;
	}


	if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h_check_k != 1 ) )
	{
		h_c_k = 0.; 
		h_d_k = 1. - h_c_k;
	}

	if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h_check_k != 1 ) )
	{
		h_c_k = 1.; 
		h_d_k = 1. - h_c_k;
	}


	if ( p_dyn.x[ im-1 ][ j ][ k ] >= .01 )			p_dyn.x[ im-1 ][ j ][ k ] = .01;
	if ( p_dyn.x[ im-1 ][ j ][ k ] <= - .01 )			p_dyn.x[ im-1 ][ j ][ k ] = - .01;

	if ( v.x[ im-1 ][ j ][ k ] >= .003 )			v.x[ im-1 ][ j ][ k ] = .003;
	if ( v.x[ im-1 ][ j ][ k ] <= - .003 )			v.x[ im-1 ][ j ][ k ] = - .003;

	if ( w.x[ im-1 ][ j ][ k ] >= .01 )				w.x[ im-1 ][ j ][ k ] = .01;
	if ( w.x[ im-1 ][ j ][ k ] <= - .01 )			w.x[ im-1 ][ j ][ k ] = .01;


	dvdthe = h_d_j * ( v.x[ im-1 ][ j+1 ][ k ] - v.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );
	dwdthe = h_d_j * ( w.x[ im-1 ][ j+1 ][ k ] - w.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );
	dpdthe = h_d_j * ( p_dyn.x[ im-1 ][ j+1 ][ k ] - p_dyn.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );

	dvdphi = h_d_k * ( v.x[ im-1 ][ j ][ k+1 ] - v.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );
	dwdphi = h_d_k * ( w.x[ im-1 ][ j ][ k+1 ] - w.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );
	dpdphi = h_d_k * ( p_dyn.x[ im-1 ][ j ][ k+1 ] - p_dyn.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );

	d2vdthe2 = h_d_j *  ( v.x[ im-1 ][ j+1 ][ k ] - 2. * v.x[ im-1 ][ j ][ k ] + v.x[ im-1 ][ j-1 ][ k ] ) / dthe2;
	d2wdthe2 = h_d_j * ( w.x[ im-1 ][ j+1 ][ k ] - 2. * w.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j-1 ][ k ] ) / dthe2;

	d2vdphi2 = h_d_k * ( v.x[ im-1 ][ j ][ k+1 ] - 2. * v.x[ im-1 ][ j ][ k ] + v.x[ im-1 ][ j ][ k-1 ] ) / dphi2;
	d2wdphi2 = h_d_k * ( w.x[ im-1 ][ j ][ k+1 ] - 2. * w.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k-1 ] ) / dphi2;


	if ( ( j >= 2 ) && ( j < jm - 3 ) )
	{
		if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( ( h.x[ im-1 ][ j + 1 ][ k ] == 0. ) && ( h.x[ im-1 ][ j + 2 ][ k ] == 0. ) ) )
		{
			dvdthe = h_d_j * ( - 3. * v.x[ im-1 ][ j ][ k ] + 4. * v.x[ im-1 ][ j + 1 ][ k ] - v.x[ im-1 ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dwdthe = h_d_j * ( - 3. * w.x[ im-1 ][ j ][ k ] + 4. * w.x[ im-1 ][ j + 1 ][ k ] - w.x[ im-1 ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dpdthe = h_d_j * ( - 3. * p_dyn.x[ im-1 ][ j ][ k ] + 4. * p_dyn.x[ im-1 ][ j + 1 ][ k ] - p_dyn.x[ im-1 ][ j + 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
		}
		else
		{
			dvdthe = h_d_j * ( v.x[ im-1 ][ j + 1 ][ k ] - v.x[ im-1 ][ j ][ k ] ) / dthe;
			dwdthe = h_d_j * ( w.x[ im-1 ][ j + 1 ][ k ] - w.x[ im-1 ][ j ][ k ] ) / dthe;
			dpdthe = h_d_j * ( p_dyn.x[ im-1 ][ j + 1 ][ k ] - p_dyn.x[ im-1 ][ j ][ k ] ) / dthe;
		}


		if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j - 1 ][ k ] == 0. ) && ( h.x[ im-1 ][ j - 2 ][ k ] == 0. ) )
		{
			dvdthe = h_d_j * ( - 3. * v.x[ im-1 ][ j ][ k ] + 4. * v.x[ im-1 ][ j - 1 ][ k ] - v.x[ im-1 ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dwdthe = h_d_j * ( - 3. * w.x[ im-1 ][ j ][ k ] + 4. * w.x[ im-1 ][ j - 1 ][ k ] - w.x[ im-1 ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
			dpdthe = h_d_j * ( - 3. * p_dyn.x[ im-1 ][ j ][ k ] + 4. * p_dyn.x[ im-1 ][ j - 1 ][ k ] - p_dyn.x[ im-1 ][ j - 2 ][ k ] ) / ( 2. * dthe );					// 2. order accurate
		}
		else
		{
			dvdthe = h_d_j * ( v.x[ im-1 ][ j ][ k ] - v.x[ im-1 ][ j - 1 ][ k ] ) / dthe;
			dwdthe = h_d_j * ( w.x[ im-1 ][ j ][ k ] - w.x[ im-1 ][ j - 1 ][ k ] ) / dthe;
			dpdthe = h_d_j * ( p_dyn.x[ im-1 ][ j ][ k ] - p_dyn.x[ im-1 ][ j - 1 ][ k ] ) / dthe;
		}

		d2vdthe2 = d2wdthe2 = 0.;
	}
	else
	{
		if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j + 1 ][ k ] == 0. ) )
		{
			dvdthe = h_d_j * ( v.x[ im-1 ][ j + 1 ][ k ] - v.x[ im-1 ][ j ][ k ] ) / dthe;
			dwdthe = h_d_j * ( w.x[ im-1 ][ j + 1 ][ k ] - w.x[ im-1 ][ j ][ k ] ) / dthe;
			dpdthe = h_d_j * ( p_dyn.x[ im-1 ][ j + 1 ][ k ] - p_dyn.x[ im-1 ][ j ][ k ] ) / dthe;
		}

		if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j - 1 ][ k ] == 0. ) )
		{
			dvdthe = h_d_j * ( v.x[ im-1 ][ j ][ k ] - v.x[ im-1 ][ j - 1 ][ k ] ) / dthe;
			dwdthe = h_d_j * ( w.x[ im-1 ][ j ][ k ] - w.x[ im-1 ][ j - 1 ][ k ] ) / dthe;
			dpdthe = h_d_j * ( p_dyn.x[ im-1 ][ j ][ k ] - p_dyn.x[ im-1 ][ j - 1 ][ k ] ) / dthe;
		}

		d2vdthe2 = d2wdthe2 = 0.;
	}


	if ( ( k >= 2 ) && ( k < km - 3 ) )
	{
		if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k + 1 ] == 0. ) && ( h.x[ im-1 ][ j ][ k + 2 ] == 0. ) )
		{
			dvdphi = h_d_k * ( - 3. * v.x[ im-1 ][ j ][ k ] + 4. * v.x[ im-1 ][ j ][ k + 1 ] - v.x[ im-1 ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dwdphi = h_d_k * ( - 3. * w.x[ im-1 ][ j ][ k ] + 4. * w.x[ im-1 ][ j ][ k + 1 ] - w.x[ im-1 ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dpdphi = h_d_k * ( - 3. * p_dyn.x[ im-1 ][ j ][ k ] + 4. * p_dyn.x[ im-1 ][ j ][ k + 1 ] - p_dyn.x[ im-1 ][ j ][ k + 2 ] ) / ( 2. * dphi );					// 2. order accurate
		}
		else
		{
			dvdphi = h_d_k * ( v.x[ im-1 ][ j ][ k + 1 ] - v.x[ im-1 ][ j ][ k ] ) / dphi;
			dwdphi = h_d_k * ( w.x[ im-1 ][ j ][ k + 1 ] - w.x[ im-1 ][ j ][ k ] ) / dphi;
			dpdphi = h_d_k * ( p_dyn.x[ im-1 ][ j ][ k + 1 ] - p_dyn.x[ im-1 ][ j ][ k ] ) / dphi;
		}


		if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k - 1 ] == 0. ) && ( h.x[ im-1 ][ j ][ k - 2 ] == 0. ) )
		{
			dvdphi = h_d_k * ( - 3. * v.x[ im-1 ][ j ][ k ] + 4. * v.x[ im-1 ][ j ][ k - 1 ] - v.x[ im-1 ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dwdphi = h_d_k * ( - 3. * w.x[ im-1 ][ j ][ k ] + 4. * w.x[ im-1 ][ j ][ k - 1 ] - w.x[ im-1 ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
			dpdphi = h_d_k * ( - 3. * p_dyn.x[ im-1 ][ j ][ k ] + 4. * p_dyn.x[ im-1 ][ j ][ k - 1 ] - p_dyn.x[ im-1 ][ j ][ k - 2 ] ) / ( 2. * dphi );					// 2. order accurate
		}
		else
		{
			dvdphi = h_d_k * ( v.x[ im-1 ][ j ][ k ] - v.x[ im-1 ][ j ][ k - 1 ] ) / dphi;
			dwdphi = h_d_k * ( w.x[ im-1 ][ j ][ k ] - w.x[ im-1 ][ j ][ k - 1 ] ) / dphi;
			dpdphi = h_d_k * ( p_dyn.x[ im-1 ][ j ][ k ] - p_dyn.x[ im-1 ][ j ][ k - 1 ] ) / dphi;
		}

		d2vdphi2 = d2wdphi2 = 0.;
	}
	else
	{
		if ( ( h.x[ im-1 ][ j ][ k ] == 1. ) && ( h.x[ im-1 ][ j ][ k + 1 ] == 0. ) )
		{
			dvdphi = h_d_k * ( v.x[ im-1 ][ j ][ k + 1 ] - v.x[ im-1 ][ j ][ k ] ) / dphi;
			dwdphi = h_d_k * ( w.x[ im-1 ][ j ][ k + 1 ] - w.x[ im-1 ][ j ][ k ] ) / dphi;
			dpdphi = h_d_k * ( p_dyn.x[ im-1 ][ j ][ k + 1 ] - p_dyn.x[ im-1 ][ j ][ k ] ) / dphi;
		}

		if ( ( h.x[ im-1 ][ j ][ k ] == 0. ) && ( h.x[ im-1 ][ j ][ k - 1 ] == 1. ) )
		{
			dvdphi = h_d_k * ( v.x[ im-1 ][ j ][ k ] - v.x[ im-1 ][ j ][ k - 1 ] ) / dphi;
			dwdphi = h_d_k * ( w.x[ im-1 ][ j ][ k ] - w.x[ im-1 ][ j ][ k - 1 ] ) / dphi;
			dpdphi = h_d_k * ( p_dyn.x[ im-1 ][ j ][ k ] - p_dyn.x[ im-1 ][ j ][ k - 1 ] ) / dphi;
		}

		d2vdphi2 = d2wdphi2 = 0.;
	}





	RS_Coriolis_Momentum_the = + coriolis * 2. * omega * costhe * w.x[ im-1 ][ j ][ k ] * h_d_j;
	RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ 0 ] * sinthe * costhe * pow ( ( omega ), 2 );
	RS_Coriolis_Momentum_phi = - coriolis * omega * costhe * v.x[ im-1 ][ j ][ k ] * h_d_k;

	rhs_v.x[ im-1 ][ j ][ k ] = - ( v.x[ im-1 ][ j ][ k ] * dvdthe / rm + w.x[ im-1 ][ j ][ k ] * dvdphi / rmsinthe ) +
				- dpdthe / rm - ( d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
				- ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[ im-1 ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
				+ RS_Coriolis_Momentum_the + RS_Centrifugal_Momentum_the
				- h_c_j * v.x[ im-1 ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition

	rhs_w.x[ im-1 ][ j ][ k ] = - ( v.x[ im-1 ][ j ][ k ] * dwdthe / rm +  w.x[ im-1 ][ j ][ k ] * dwdphi / rmsinthe ) +
				- dpdphi / rmsinthe + ( d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
				- ( 1. + costhe * costhe / sinthe2 ) * h_d_k * w.x[ im-1 ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2 + dvdphi * 2. * costhe / rm2sinthe2 ) / re
				+ RS_Coriolis_Momentum_phi
				- h_c_k * w.x[ im-1 ][ j ][ k ] * k_Force / dphi2;					// immersed boundary condition as a negative force addition

	rhs_p.x[ im-1 ][ j ][ k ] = - ( v.x[ im-1 ][ j ][ k ] * ( rhs_v.x[ im-1 ][ j ][ k ] + dpdthe / rm ) + w.x[ im-1 ][ j ][ k ] * ( rhs_w.x[ im-1 ][ j ][ k ] + dpdphi / rmsinthe ) )
				- h_c_k * p_dyn.x[ im-1 ][ j ][ k ] * k_Force / dphi2;					// immersed boundary condition as a negative force addition

	aux_v.x[ im-1 ][ j ][ k ] = rhs_v.x[ im-1 ][ j ][ k ] + dpdthe / rm;
	aux_w.x[ im-1 ][ j ][ k ] = rhs_w.x[ im-1 ][ j ][ k ] + dpdphi / rmsinthe;
 

	if ( h.x[ im-1 ][ j ][ k ] == 1. )					aux_v.x[ im-1 ][ j ][ k ] = aux_w.x[ im-1 ][ j ][ k ] = 0.;
	if ( h.x[ im-1 ][ j ][ k ] == 1. )					v.x[ im-1 ][ j ][ k ] = w.x[ im-1 ][ j ][ k ] = 0.;
}





