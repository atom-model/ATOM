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


void RHS_Hydrosphere::RK_RHS_3D_Hydrosphere ( int i, int j, int k, double L_hyd, double g, double cp_w, double u_0, double t_0, double c_0, double r_0_water, double ta, double pa, double ca, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c, Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion, Array &BuoyancyForce_3D, Array &Salt_Balance )
{

// collection of coefficients for phase transformation

	L_hyd = L_hyd / ( im - 1 );														// characteristic length for non-dimensionalisation

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
//	rm = rad.z[ i ];
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

	if ( ( h.x[ i ][ j +1 ][ k ] == 1. ) && ( h.x[ i ][ j ][ k ] == 0. ) )
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

/*
	dudr = h_d_i * ( u.x[ i+1 ][ j ][ k ] - u.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * rm_1;
	dvdr = h_d_i * ( v.x[ i+1 ][ j ][ k ] - v.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * rm_1;
	dwdr = h_d_i * ( w.x[ i+1 ][ j ][ k ] - w.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * rm_1;
	dtdr = h_d_i * ( t.x[ i+1 ][ j ][ k ] - t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * rm_1;
	dpdr = h_d_i * ( p_dyn.x[ i+1 ][ j ][ k ] - p_dyn.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * rm_1;
	dcdr = h_d_i * ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr ) * rm_1;
*/
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

/*
	d2udr2 = h_d_i * ( u.x[ i+1 ][ j ][ k ] - 2. * u.x[ i ][ j ][ k ] + u.x[ i-1 ][ j ][ k ] ) / dr2 * rm_2 - 2. / rm_1 * dudr;
	d2vdr2 = h_d_i * ( v.x[ i+1 ][ j ][ k ] - 2. * v.x[ i ][ j ][ k ] + v.x[ i-1 ][ j ][ k ] ) / dr2 * rm_2 - 2. / rm_1 * dvdr;
	d2wdr2 = h_d_i * ( w.x[ i+1 ][ j ][ k ] - 2. * w.x[ i ][ j ][ k ] + w.x[ i-1 ][ j ][ k ] ) / dr2 * rm_2 - 2. / rm_1 * dwdr;
	d2tdr2 = h_d_i * ( t.x[ i+1 ][ j ][ k ] - 2. * t.x[ i ][ j ][ k ] + t.x[ i-1 ][ j ][ k ] ) / dr2 * rm_2 - 2. / rm_1 * dtdr;
	d2cdr2 = h_d_i * ( c.x[ i+1 ][ j ][ k ] - 2. * c.x[ i ][ j ][ k ] + c.x[ i-1 ][ j ][ k ] ) / dr2 * rm_2 - 2. / rm_1 * dcdr;
*/
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

/*
	if ( t.x[ i ][ j ][ k ] < ta ) 
	{
		tn.x[ i ][ j ][ k ] = ta;
		cn.x[ i ][ j ][ k ] = ca;
		c_Boussinesq = ( ( ( ta * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
	}
	else
	{
		c_Boussinesq = ( ( ( tn.x[ i ][ j ][ k ] * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
	}
*/
//		c_Boussinesq = .9682;															// for c = 0.9682 compares to a salinity of 33.5 psu
//		c_Boussinesq = .7225;															// for c = 0.7225 compares to a salinity of 25.0 psu taken from Burchard
//		c_Boussinesq = 1.01156;														// for c = 1.01156 compares to a salinity of 35.0 psu, ca corresponds to ta = 1.01464  ( = 4°C )
		c_Boussinesq = 1.0571;															// for c = 1.0571 compares to a salinity of 36.58 psu

//		c_Boussinesq = r_0_water;														// compares to 1025 kg/m³ of salt water, for c = 0.7225 compares to a salinity of 25.0 psu, taken from Burchard


	RS_buoyancy_Momentum = - L_hyd / ( u_0 * u_0 ) * buoyancy * g * ( ( c.x[ i ][ j ][ k ] * c_0 + r_0_water ) - ( c_Boussinesq * c_0 + r_0_water ) ) / ( c_Boussinesq * c_0 + r_0_water );																								// bouyancy based on salt water density 
//	RS_buoyancy_Momentum = - L_hyd / ( u_0 * u_0 ) * buoyancy * g * ( ( c.x[ i ][ j ][ k ] * c_0 + 1000. ) - c_Boussinesq ) / c_Boussinesq;																								// bouyancy based on salt water density 
//	RS_buoyancy_Momentum = 0.;
	RS_buoyancy_Energy = u_0 * ec * u.x[ i ][ j ][ k ] * RS_buoyancy_Momentum;		// bouyancy based on salt water density
	BuoyancyForce_3D.x[ i ][ j ][ k ] = RS_buoyancy_Momentum;
 	Salt_Balance.x[ i ][ j ][ k ] = RS_Salt_Balance = ( c.x[ i ][ j ][ k ] - c_Boussinesq ) * c_0; 

//	cout << i << "   " << j << "   " << k << "   " << RS_Salt_Balance << "   " << c.x[ i ][ j ][ k ] << "   " << c_Boussinesq << "   " << RS_buoyancy_Momentum << "   " << g << "   " << buoyancy << endl;

	if ( Salt_Balance.x[ i ][ j ][ k ] >= 0. )
	{
		Salt_Finger.x[ i ][ j ][ k ] = Salt_Balance.x[ i ][ j ][ k ];
	}
	else
	{
		Salt_Finger.x[ i ][ j ][ k ] = 0.;
	}

	if ( Salt_Balance.x[ i ][ j ][ k ] < 0. )
	{
		Salt_Diffusion.x[ i ][ j ][ k ] = Salt_Balance.x[ i ][ j ][ k ];
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



/*
			cout << endl;
			cout << " i = " << i << "   j = " << j << "   k = " << k << endl;
			cout << " RS_Coriolis_Momentum_rad = " << RS_Coriolis_Momentum_rad * 1000. << "     RS_Coriolis_Momentum_the = " << RS_Coriolis_Momentum_the * 1000. << "     RS_Coriolis_Momentum_phi = " << RS_Coriolis_Momentum_phi * 1000. << endl << endl;
			cout << " RS_Centrifugal_Momentum_rad = " << RS_Centrifugal_Momentum_rad * 1000. << "     RS_Centrifugal_Momentum_the = " << RS_Centrifugal_Momentum_the * 1000. << endl;
			cout << " RS_Coriolis_Energy = " << RS_Coriolis_Energy * 1000. << "     RS_Centrifugal_Energy = " << RS_Centrifugal_Energy * 1000. << endl;
			cout << endl;
*/


// Right Hand Side of the time derivative ot temperature, pressure, salt concentration and velocity components

//  3D volume iterations

	rhs_t.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dtdr + v.x[ i ][ j ][ k ] * dtdthe / rm + w.x[ i ][ j ][ k ] * dtdphi / rmsinthe )
			- .5 * ec * ( u.x[ i ][ j ][ k ] * dpdr + v.x[ i ][ j ][ k ] / rm * dpdthe + w.x[ i ][ j ][ k ] / rmsinthe * dpdphi )
			+ ( d2tdr2 + dtdr * 2. / rm + d2tdthe2 / rm2 + dtdthe * costhe / rm2sinthe + d2tdphi2 / rm2sinthe2 ) / ( re * pr )
			+ 2. * ec / re * ( ( dudr * dudr) + pow ( ( dvdthe / rm + h_d_i * u.x[ i ][ j ][ k ] / rm ), 2. )
			+ pow ( ( dwdphi / rmsinthe + h_d_i * u.x[ i ][ j ][ k ] / rm + h_d_j * v.x[ i ][ j ][ k ] * cotthe / rm ), 2. ) )
			+ ec / re * ( pow ( ( dvdr - h_d_j * v.x[ i ][ j ][ k ] / rm + dudthe / rm ), 2. ) + pow ( ( dudphi / rmsinthe + dwdr - h_d_k * w.x[ i ][ j ][ k ] / rm ), 2. )
			+ pow ( ( dwdthe * sinthe / rm2 - h_d_k * w.x[ i ][ j ][ k ] * costhe / rmsinthe + dvdphi / rmsinthe ), 2. ) )
			+ RS_buoyancy_Energy
			+ RS_Coriolis_Energy + RS_Centrifugal_Energy;
//			- h_c_i * t.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_u.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dudr + v.x[ i ][ j ][ k ] * dudthe / rm + w.x[ i ][ j ][ k ] * dudphi / rmsinthe )
			- .5 * dpdr + ( d2udr2 + h_d_i * 2. * u.x[ i ][ j ][ k ] / rm2 + d2udthe2 / rm2 + 4. * dudr / rm + dudthe * costhe / rm2sinthe + d2udphi2 / rm2sinthe2 ) / re
			+ RS_buoyancy_Momentum
			+ RS_Coriolis_Momentum_rad + RS_Centrifugal_Momentum_rad
			- h_c_i * u.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_v.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dvdr + v.x[ i ][ j ][ k ] * dvdthe / rm + w.x[ i ][ j ][ k ] * dvdphi / rmsinthe )
			- .5 * dpdthe / rm + ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
			- ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[ i ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2
			+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_the + RS_Centrifugal_Momentum_the
			- h_c_j * v.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	rhs_w.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dwdr + v.x[ i ][ j ][ k ] * dwdthe / rm + w.x[ i ][ j ][ k ] * dwdphi / rmsinthe )
			- .5 * dpdphi / rmsinthe + ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
			- ( 1. + costhe * costhe / sinthe2 ) * h_d_k * w.x[ i ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2
			+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
			+ RS_Coriolis_Momentum_phi
			- h_c_k * w.x[ i ][ j ][ k ] * k_Force / dphi2;					// immersed boundary condition as a negative force addition


	rhs_c.x[ i ][ j ][ k ] = - ( u.x[ i ][ j ][ k ] * dcdr + v.x[ i ][ j ][ k ] * dcdthe / rm + w.x[ i ][ j ][ k ] * dcdphi / rmsinthe )
			+ ( d2cdr2 + dcdr * 2. / rm + d2cdthe2 / rm2 + dcdthe * costhe / rm2sinthe + d2cdphi2 / rm2sinthe2 ) / ( sc * re );
//			- h_c_i * c.x[ i ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition


	// for the Poisson equation to solve for the pressure, pressure gradient sbstracted from the RHS

	aux_u.x[ i ][ j ][ k ] = rhs_u.x[ i ][ j ][ k ] + dpdr;
	aux_v.x[ i ][ j ][ k ] = rhs_v.x[ i ][ j ][ k ] + dpdthe / rm;
	aux_w.x[ i ][ j ][ k ] = rhs_w.x[ i ][ j ][ k ] + dpdphi / rmsinthe;

}




void RHS_Hydrosphere::RK_RHS_2D_Hydrosphere ( int j, int k, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &v, Array &w, Array &p_dyn, Array &vn, Array &wn, Array &rhs_v, Array &rhs_w, Array &aux_v, Array &aux_w )
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

	rm = rad.z[ im-1 ] / r0;
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

	if ( ( h.x[ im - 1 ][ j ][ k ] == 1. ) && ( h.x[ im - 2 ][ j ][ k ] == 0. ) )
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
*/



// 2D adapted immersed boundary method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// only in positive the-direction along northerly boundaries 

	if ( ( h.x[ im - 1 ][ j - 1 ][ k ] == 1. ) && ( h.x[ im - 1 ][ j ][ k ] == 0. ) )
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

	if ( ( h.x[ im - 1 ][ j + 1 ][ k ] == 1. ) && ( h.x[ im - 1 ][ j ][ k ] == 0. ) )
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

	if ( ( h.x[ im - 1 ][ j ][ k - 1 ] == 1. ) && ( h.x[ im - 1 ][ j ][ k ] == 0. ) )
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

	if ( ( h.x[ im - 1 ][ j ][ k + 1 ] == 1. ) && ( h.x[ im - 1 ][ j ][ k ] == 0. ) )
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


		if ( ( h.x[ im - 1 ][ j ][ k ] == 0. ) && ( h_check_j != 1 ) )
		{
			h_c_j = 0.; 
			h_d_j = 1. - h_c_j;
		}

		if ( ( h.x[ im - 1 ][ j ][ k ] == 1. ) && ( h_check_j != 1 ) )
		{
			h_c_j = 1.; 
			h_d_j = 1. - h_c_j;
		}


		if ( ( h.x[ im - 1 ][ j ][ k ] == 0. ) && ( h_check_k != 1 ) )
		{
			h_c_k = 0.; 
			h_d_k = 1. - h_c_k;
		}

		if ( ( h.x[ im - 1 ][ j ][ k ] == 1. ) && ( h_check_k != 1 ) )
		{
			h_c_k = 1.; 
			h_d_k = 1. - h_c_k;
		}


		dvdthe = h_d_j * ( v.x[ im-1 ][ j+1 ][ k ] - v.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );
		dwdthe = h_d_j * ( w.x[ im-1 ][ j+1 ][ k ] - w.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );
		dpdthe = h_d_j * ( p_dyn.x[ im-1 ][ j+1 ][ k ] - p_dyn.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );

		dvdphi = h_d_k * ( v.x[ im-1 ][ j ][ k+1 ] - v.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );
		dwdphi = h_d_k * ( w.x[ im-1 ][ j ][ k+1 ] - w.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );
		dpdphi = h_d_k * ( p_dyn.x[ im-1 ][ j ][ k+1 ] - p_dyn.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );

		d2vdthe2 = h_d_j * ( v.x[ im-1 ][ j+1 ][ k ] - 2. * v.x[ im-1 ][ j ][ k ] + v.x[ im-1 ][ j-1 ][ k ] ) / dthe2;
		d2wdthe2 = h_d_j * ( w.x[ im-1 ][ j+1 ][ k ] - 2. * w.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j-1 ][ k ] ) / dthe2;

		d2vdphi2 = h_d_k * ( v.x[ im-1 ][ j ][ k+1 ] - 2. * v.x[ im-1 ][ j ][ k ] + v.x[ im-1 ][ j ][ k-1 ] ) / dphi2;
		d2wdphi2 = h_d_k * ( w.x[ im-1 ][ j ][ k+1 ] - 2. * w.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k-1 ] ) / dphi2;

		RS_Coriolis_Momentum_the = + coriolis * 2. * omega * costhe * w.x[ im-1 ][ j ][ k ] * h_d_j;
		RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ im-1 ] * sinthe * costhe * pow ( ( omega ), 2 );
		RS_Coriolis_Momentum_phi = - coriolis * omega * costhe * v.x[ im-1 ][ j ][ k ] * h_d_k;

		rhs_v.x[ im-1 ][ j ][ k ] = - ( v.x[ im-1 ][ j ][ k ] * dvdthe / rm + w.x[ im-1 ][ j ][ k ] * dvdphi / rmsinthe ) +
				- .5 * dpdthe / rm + ( d2vdr2 + 2. * dvdr / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
				- ( 1. + costhe * costhe / sinthe2 ) * h_d_j * v.x[ im-1 ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
				+ RS_Coriolis_Momentum_the + RS_Centrifugal_Momentum_the
				- h_c_j * v.x[ im - 1 ][ j ][ k ] * k_Force / dthe2;					// immersed boundary condition as a negative force addition

		rhs_w.x[ im-1 ][ j ][ k ] = - ( v.x[ im-1 ][ j ][ k ] * dwdthe / rm + w.x[ im-1 ][ j ][ k ] * dwdphi / rmsinthe ) +
				- .5 * dpdphi / rmsinthe + ( d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
				- ( 1. + costhe * costhe / sinthe2 ) * h_d_k * w.x[ im-1 ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2	+ dvdphi * 2. * costhe / rm2sinthe2 ) / re
				+ RS_Coriolis_Momentum_phi
				- h_c_k * w.x[ im - 1 ][ j ][ k ] * k_Force / dphi2;					// immersed boundary condition as a negative force addition

 
		aux_v.x[ im-1 ][ j ][ k ] = rhs_v.x[ im-1 ][ j ][ k ] + dpdthe / rm;
		aux_w.x[ im-1 ][ j ][ k ] = rhs_w.x[ im-1 ][ j ][ k ] + dpdphi / rmsinthe;

}







void RHS_Hydrosphere::Pressure_RHS_Hydrosphere ( double ta, double ca, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c, Array &aux_u, Array &aux_v, Array &aux_w )
{
// collection of coefficients for phase transformation
	c43 = 4. / 3.;
	c13 = 1. / 3.;

	k_Force = 10.;																			// factor for accelleration of convergence processes inside the immersed boundary conditions


// 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )

// collection of coefficients
	dr2 = dr * dr;
	dthe2 = dthe * dthe;
	dphi2 = dphi * dphi;

//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    level i = 0    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// collection of coefficients
	rm = rad.z[ 0 ];
	rm2 = rm * rm;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
// collection of coefficients
			sinthe = sin( the.z[ j ] );
			sinthe2 = sinthe * sinthe;
			costhe = cos( the.z[ j ] );
			cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
			rmsinthe = rm * sinthe;
			rm2sinthe = rm2 * sinthe;
			rm2sinthe2 = rm2 * sinthe2;

// 1st order derivative for pressure and velocity components
			dudr = ( 4. * u.x[ 1 ][ j ][ k ] - u.x[ 2 ][ j ][ k ] ) / ( 2. * dr );
			dvdr = ( - 3. * v.x[ 0 ][ j ][ k ] + 4. * v.x[ 1 ][ j ][ k ] - v.x[ 2 ][ j ][ k ] ) / ( 2. * dr );
			dwdr = ( - 3. * w.x[ 0 ][ j ][ k ] + 4. * w.x[ 1 ][ j ][ k ] - w.x[ 2 ][ j ][ k ] ) / ( 2. * dr );

			dvdthe = ( v.x[ 0 ][ j+1 ][ k ] - v.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );
			dwdthe = ( w.x[ 0 ][ j+1 ][ k ] - w.x[ 0 ][ j-1 ][ k ] ) / ( 2. * dthe );

			dvdphi = ( v.x[ 0 ][ j ][ k+1 ] - v.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );
			dwdphi = ( w.x[ 0 ][ j ][ k+1 ] - w.x[ 0 ][ j ][ k-1 ] ) / ( 2. * dphi );

// 2nd order derivative for pressure and velocity components
			d2udr2 = ( - 2. * u.x[ 1 ][ j ][ k ] + u.x[ 2 ][ j ][ k ] ) / dr2;
			d2vdr2 = ( v.x[ 0 ][ j ][ k ] - 2. * v.x[ 1 ][ j ][ k ] + v.x[ 2 ][ j ][ k ] ) / dr2;
			d2wdr2 = ( w.x[ 0 ][ j ][ k ] - 2. * w.x[ 1 ][ j ][ k ] + w.x[ 2 ][ j ][ k ] ) / dr2;

			d2vdthe2 = ( v.x[ 0 ][ j+1 ][ k ] - 2. * v.x[ 0 ][ j ][ k ] + v.x[ 0 ][ j-1 ][ k ] ) / dthe2;
			d2wdthe2 = ( w.x[ 0 ][ j+1 ][ k ] - 2. * w.x[ 0 ][ j ][ k ] + w.x[ 0 ][ j-1 ][ k ] ) / dthe2;

			d2vdphi2 = ( v.x[ 0 ][ j ][ k+1 ] - 2. * v.x[ 0 ][ j ][ k ] + v.x[ 0 ][ j ][ k-1 ] ) / dphi2;
			d2wdphi2 = ( w.x[ 0 ][ j ][ k+1 ] - 2. * w.x[ 0 ][ j ][ k ] + w.x[ 0 ][ j ][ k-1 ] ) / dphi2;

// Coriolis and centrifugal terms in the energy and momentum equations
			RS_Coriolis_Momentum_rad = + coriolis * 2. * omega * sinthe * w.x[ 0 ][ j ][ k ];
			RS_Coriolis_Momentum_the = + coriolis * 2. * omega * costhe * w.x[ 0 ][ j ][ k ];
			RS_Coriolis_Momentum_phi = - coriolis * ( 2. * omega * sinthe * u.x[ 0 ][ j ][ k ] + 2. * omega * costhe * v.x[ 0 ][ j ][ k ] );

			RS_Centrifugal_Momentum_rad = + centrifugal * rad.z[ 0 ] * pow ( ( omega * sinthe ), 2 );
			RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ 0 ] * sinthe * costhe * pow ( ( omega ), 2 );


			if ( t.x[ 0 ][ j ][ k ] < ta ) 
			{
				tn.x[ 0 ][ j ][ k ] = ta;
				cn.x[ 0 ][ j ][ k ] = ca;
				c_Boussinesq = ( ( ( ta * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}
			else
			{
				c_Boussinesq = ( ( ( tn.x[ 0 ][ j ][ k ] * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}

			RS_buoyancy_Momentum = - buoyancy * .5 * g * ( c.x[ 0 ][ j ][ k ] - c_Boussinesq ) / c_Boussinesq;		// bouyancy based on water vapour 


// Right Hand Side of the Navier-Stokes equations at the boundaries for the pressure Poisson equation
			aux_u.x[ 0 ][ j ][ k ] = ( d2udr2 + d2udthe2 / rm2 + 4. * dudr / rm + d2udphi2 / rm2sinthe2 ) / re
												+ RS_buoyancy_Momentum
												+ RS_Coriolis_Momentum_rad + RS_Centrifugal_Momentum_rad;

			aux_v.x[ 0 ][ j ][ k ] = - ( v.x[ 0 ][ j ][ k ] * dvdthe / rm + w.x[ 0 ][ j ][ k ] * dvdphi / rmsinthe ) +
												+ ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
												- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * v.x[ 0 ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2
												+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
												+ RS_Coriolis_Momentum_the + RS_Centrifugal_Momentum_the
												- v.x[ 0 ][ j ][ k ] * k_Force / dthe2;

			aux_w.x[ 0 ][ j ][ k ] = - ( v.x[ 0 ][ j ][ k ] * dwdthe / rm + w.x[ 0 ][ j ][ k ] * dwdphi / rmsinthe ) +
												+ ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
												- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * w.x[ 0 ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2
												+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
												+ RS_Coriolis_Momentum_phi
												- w.x[ 0 ][ j ][ k ] * k_Force / dphi2;



//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    level i = im-1    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// collection of coefficients
			rm = rad.z[ im-1 ];
			rm2 = rm * rm;

// collection of coefficients
			sinthe = sin( the.z[ j ] );
			sinthe2 = sinthe * sinthe;
			costhe = cos( the.z[ j ] );
			cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
			rmsinthe = rm * sinthe;
			rm2sinthe = rm2 * sinthe;
			rm2sinthe2 = rm2 * sinthe2;

// 1st order derivative for pressure and velocity components
			dudr = ( 4. * u.x[ im-2 ][ j ][ k ] - u.x[ im-3 ][ j ][ k ] ) / ( 2. * dr );
			dvdr = ( - 3. * v.x[ im-1 ][ j ][ k ] + 4. * v.x[ im-2 ][ j ][ k ] - v.x[ im-3 ][ j ][ k ] ) / ( 2. * dr );
			dwdr = ( - 3. * w.x[ im-1 ][ j ][ k ] + 4. * w.x[ im-2 ][ j ][ k ] - w.x[ im-3 ][ j ][ k ] ) / ( 2. * dr );

			dvdthe = ( v.x[ im-1 ][ j+1 ][ k ] - v.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );
			dwdthe = ( w.x[ im-1 ][ j+1 ][ k ] - w.x[ im-1 ][ j-1 ][ k ] ) / ( 2. * dthe );

			dvdphi = ( v.x[ im-1 ][ j ][ k+1 ] - v.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );
			dwdphi = ( w.x[ im-1 ][ j ][ k+1 ] - w.x[ im-1 ][ j ][ k-1 ] ) / ( 2. * dphi );

// 2nd order derivative for pressure and velocity components
			d2udr2 = ( - 2. * u.x[ im-2 ][ j ][ k ] + u.x[ im-3 ][ j ][ k ] ) / dr2;
			d2vdr2 = ( v.x[ im-1 ][ j ][ k ] - 2. * v.x[ im-2 ][ j ][ k ] + v.x[ im-3 ][ j ][ k ] ) / dr2;
			d2wdr2 = ( w.x[ im-1 ][ j ][ k ] - 2. * w.x[ im-2 ][ j ][ k ] + w.x[ im-3 ][ j ][ k ] ) / dr2;

			d2vdthe2 = ( v.x[ im-1 ][ j+1 ][ k ] - 2. * v.x[ im-1 ][ j ][ k ] + v.x[ im-1 ][ j-1 ][ k ] ) / dthe2;
			d2wdthe2 = ( w.x[ im-1 ][ j+1 ][ k ] - 2. * w.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j-1 ][ k ] ) / dthe2;

			d2vdphi2 = ( v.x[ im-1 ][ j ][ k+1 ] - 2. * v.x[ im-1 ][ j ][ k ] + v.x[ im-1 ][ j ][ k-1 ] ) / dphi2;
			d2wdphi2 = ( w.x[ im-1 ][ j ][ k+1 ] - 2. * w.x[ im-1 ][ j ][ k ] + w.x[ im-1 ][ j ][ k-1 ] ) / dphi2;

// Coriolis and centrifugal terms in the energy and momentum equations
			RS_Coriolis_Momentum_rad = + coriolis * 2. * omega * sinthe * w.x[ im-1 ][ j ][ k ];
			RS_Coriolis_Momentum_the = + coriolis * 2. * omega * costhe * w.x[ im-1 ][ j ][ k ];
			RS_Coriolis_Momentum_phi = - coriolis * ( 2. * omega * sinthe * u.x[ im-1 ][ j ][ k ] + 2. * omega * costhe * v.x[ im-1 ][ j ][ k ] );

			RS_Centrifugal_Momentum_rad = + centrifugal * rad.z[ im-1 ] * pow ( ( omega * sinthe ), 2 );
			RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ im-1 ] * sinthe * costhe * pow ( ( omega ), 2 );

			if ( t.x[ im-1 ][ j ][ k ] < ta ) 
			{
				tn.x[ im-1 ][ j ][ k ] = ta;
				cn.x[ im-1 ][ j ][ k ] = ca;
				c_Boussinesq = ( ( ( ta * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}
			else
			{
				c_Boussinesq = ( ( ( tn.x[ im-1 ][ j ][ k ] * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}

			RS_buoyancy_Momentum = - buoyancy * .5 * g * ( c.x[ im-1 ][ j ][ k ] - c_Boussinesq ) / c_Boussinesq;		// bouyancy based on water vapour 

// Right Hand Side of the Navier-Stokes equations at the boundaries for the pressure Poisson equation
			aux_u.x[ im-1 ][ j ][ k ] = + ( d2udr2 + d2udthe2 / rm2 + 4. * dudr / rm + d2udphi2 / rm2sinthe2 ) / re
													+ RS_buoyancy_Momentum
													+ RS_Coriolis_Momentum_rad + RS_Centrifugal_Momentum_rad;

			aux_v.x[ im-1 ][ j ][ k ] =  - ( v.x[ im-1 ][ j ][ k ] * dvdthe / rm + w.x[ im-1 ][ j ][ k ] * dvdphi / rmsinthe ) +
													+ ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
													- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * v.x[ im-1 ][ j ][ k ] / rm + d2vdphi2 / rm2sinthe2
													+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
													+ RS_Coriolis_Momentum_the + RS_Centrifugal_Momentum_the
													- v.x[ im-1 ][ j ][ k ] * k_Force / dthe2;

			aux_w.x[ im-1 ][ j ][ k ] =  - ( v.x[ im-1 ][ j ][ k ] * dwdthe / rm + w.x[ im-1 ][ j ][ k ] * dwdphi / rmsinthe ) +
														+ ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
														- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * w.x[ im-1 ][ j ][ k ] / rm + d2wdphi2 / rm2sinthe2
														+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
														+ RS_Coriolis_Momentum_phi
														- w.x[ im-1 ][ j ][ k ] * k_Force / dphi2;
		}
	}



//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    level k = 0    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	for ( int i = 1; i < im-1; i++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
// collection of coefficients
			rm = rad.z[ i ];
			rm2 = rm * rm;
			sinthe = sin( the.z[ j ] );
			sinthe2 = sinthe * sinthe;
			costhe = cos( the.z[ j ] );
			cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
			rmsinthe = rm * sinthe;
			rm2sinthe = rm2 * sinthe;
			rm2sinthe2 = rm2 * sinthe2;

// 1st order derivative for pressure and velocity components
			dudr = ( u.x[ i+1 ][ j ][ 0 ] - u.x[ i-1 ][ j ][ 0 ] ) / ( 2. * dr );
			dvdr = ( v.x[ i+1 ][ j ][ 0 ] - v.x[ i-1 ][ j ][ 0 ] ) / ( 2. * dr );
			dwdr = ( w.x[ i+1 ][ j ][ 0 ] - w.x[ i-1 ][ j ][ 0 ] ) / ( 2. * dr );

			dudthe = ( u.x[ i ][ j+1 ][ 0 ] - u.x[ i ][ j-1 ][ 0 ] ) / ( 2. * dthe );
			dvdthe = ( v.x[ i ][ j+1 ][ 0 ] - v.x[ i ][ j-1 ][ 0 ] ) / ( 2. * dthe );
			dwdthe = ( w.x[ i ][ j+1 ][ 0 ] - w.x[ i ][ j-1 ][ 0 ] ) / ( 2. * dthe );

			dudphi = ( - 3. * u.x[ i ][ j ][ 0 ] + 4. * u.x[ i ][ j ][ 1 ] - u.x[ i ][ j ][ 2 ] ) / ( 2. * dphi );
			dvdphi = ( - 3. * v.x[ i ][ j ][ 0 ] + 4. * v.x[ i ][ j ][ 1 ] - v.x[ i ][ j ][ 2 ] ) / ( 2. * dphi );
			dwdphi = ( - 3. * w.x[ i ][ j ][ 0 ] + 4. * w.x[ i ][ j ][ 1 ] - w.x[ i ][ j ][ 2 ] ) / ( 2. * dphi );

// 2nd order derivative for pressure and velocity components
			d2udr2 = ( u.x[ i ][ j ][ 0 ] - 2. * u.x[ i ][ j ][ 0 ] + u.x[ i ][ j ][ 0 ] ) / dr2;
			d2vdr2 = ( v.x[ i ][ j ][ 0 ] - 2. * v.x[ i ][ j ][ 0 ] + v.x[ i ][ j ][ 0 ] ) / dr2;
			d2wdr2 = ( w.x[ i ][ j ][ 0 ] - 2. * w.x[ i ][ j ][ 0 ] + w.x[ i ][ j ][ 0 ] ) / dr2;

			d2udthe2 = ( u.x[ i ][ j+1 ][ 0 ] - 2. * u.x[ i ][ j ][ 0 ] + u.x[ i ][ j-1 ][ 0 ] ) / dthe2;
			d2vdthe2 = ( v.x[ i ][ j+1 ][ 0 ] - 2. * v.x[ i ][ j ][ 0 ] + v.x[ i ][ j-1 ][ 0 ] ) / dthe2;
			d2wdthe2 = ( w.x[ i ][ j+1 ][ 0 ] - 2. * w.x[ i ][ j ][ 0 ] + w.x[ i ][ j-1 ][ 0 ] ) / dthe2;

			d2udphi2 = ( u.x[ i ][ j ][ 0 ] - 2. * u.x[ i ][ j ][ 1 ] + u.x[ i ][ j ][ 2 ] ) / dphi2;
			d2vdphi2 = ( v.x[ i ][ j ][ 0 ] - 2. * v.x[ i ][ j ][ 1 ] + v.x[ i ][ j ][ 2 ] ) / dphi2;
			d2wdphi2 = ( w.x[ i ][ j ][ 0 ] - 2. * w.x[ i ][ j ][ 1 ] + w.x[ i ][ j ][ 2 ] ) / dphi2;

// Coriolis and centrifugal terms in the energy and momentum equations
			RS_Coriolis_Momentum_rad = + coriolis * 2. * omega * sinthe * w.x[ i ][ j ][ 0 ];
			RS_Coriolis_Momentum_the = + coriolis * 2. * omega * costhe * w.x[ i ][ j ][ 0 ];
			RS_Coriolis_Momentum_phi = - coriolis * ( 2. * omega * sinthe * u.x[ i ][ j ][ 0 ] + 2. * omega * costhe * v.x[ i ][ j ][ 0 ] );

			RS_Centrifugal_Momentum_rad = + centrifugal * rad.z[ i ] * pow ( ( omega * sinthe ), 2 );
			RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ i ] * sinthe * costhe * pow ( ( omega ), 2 );

			if ( t.x[ i ][ j ][ 0 ] < ta ) 
			{
				tn.x[ i ][ j ][ 0 ] = ta;
				cn.x[ i ][ j ][ 0 ] = ca;
				c_Boussinesq = ( ( ( ta * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}
			else
			{
				c_Boussinesq = ( ( ( tn.x[ i ][ j ][ 0 ] * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}

			RS_buoyancy_Momentum = - buoyancy * .5 * g * ( c.x[ i ][ j ][ 0 ] - c_Boussinesq ) / c_Boussinesq;		// bouyancy based on water vapour 

// Right Hand Side of the Navier-Stokes equations at the boundaries for the pressure Poisson equation
			aux_u.x[ i ][ j ][ 0 ] = - ( u.x[ i ][ j ][ 0 ] * dudr + v.x[ i ][ j ][ 0 ] * dudthe / rm + w.x[ i ][ j ][ 0 ] * dudphi / rmsinthe )
												- ( d2udr2 + h_d_i * 2. * u.x[ i ][ j ][ 0 ] / rm2 + d2udthe2 / rm2 + 4. * dudr / rm + dudthe * costhe / rm2sinthe + d2udphi2 / rm2sinthe2 ) / re
												+ RS_buoyancy_Momentum
												+ RS_Coriolis_Momentum_rad + RS_Centrifugal_Momentum_rad;

			aux_v.x[ i ][ j ][ 0 ] = - ( u.x[ i ][ j ][ 0 ] * dvdr + v.x[ i ][ j ][ 0 ] * dvdthe / rm + w.x[ i ][ j ][ 0 ] * dvdphi / rmsinthe ) +
												+ ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
												- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * v.x[ i ][ j ][ 0 ] / rm + d2vdphi2 / rm2sinthe2
												+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
												+ RS_Coriolis_Momentum_the + RS_Centrifugal_Momentum_the
												- v.x[ i ][ j ][ 0 ] * k_Force / dthe2;

			aux_w.x[ i ][ j ][ 0 ] = - ( u.x[ i ][ j ][ 0 ] * dwdr + v.x[ i ][ j ][ 0 ] * dwdthe / rm + w.x[ i ][ j ][ 0 ] * dwdphi / rmsinthe ) +
												+ ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
												- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * w.x[ i ][ j ][ 0 ] / rm + d2wdphi2 / rm2sinthe2
												+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
												+ RS_Coriolis_Momentum_phi
												- w.x[ i ][ j ][ 0 ] * k_Force / dphi2;



//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    level k = km-1    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// collection of coefficients
			rm = rad.z[ i ];
			rm2 = rm * rm;
			sinthe = sin( the.z[ j ] );
			sinthe2 = sinthe * sinthe;
			costhe = cos( the.z[ j ] );
			cotthe = cos( the.z[ j ] ) / sin( the.z[ j ] );
			rmsinthe = rm * sinthe;
			rm2sinthe = rm2 * sinthe;
			rm2sinthe2 = rm2 * sinthe2;

// 1st order derivative for pressure and velocity components
			dudr = ( u.x[ i+1 ][ j ][ km-1 ] - u.x[ i-1 ][ j ][ km-1 ] ) / ( 2. * dr );
			dvdr = ( v.x[ i+1 ][ j ][ km-1 ] - v.x[ i-1 ][ j ][ km-1 ] ) / ( 2. * dr );
			dwdr = ( w.x[ i+1 ][ j ][ km-1 ] - w.x[ i-1 ][ j ][ km-1 ] ) / ( 2. * dr );

			dudthe = ( u.x[ i ][ j+1 ][ km-1 ] - u.x[ i ][ j-1 ][ km-1 ] ) / ( 2. * dthe );
			dvdthe = ( v.x[ i ][ j+1 ][ km-1 ] - v.x[ i ][ j-1 ][ km-1 ] ) / ( 2. * dthe );
			dwdthe = ( w.x[ i ][ j+1 ][ km-1 ] - w.x[ i ][ j-1 ][ km-1 ] ) / ( 2. * dthe );

			dudphi = ( - 3. * u.x[ i ][ j ][ km-1 ] + 4. * u.x[ i ][ j ][ km-2 ] - u.x[ i ][ j ][ km-3 ] ) / ( 2. * dphi );
			dvdphi = ( - 3. * v.x[ i ][ j ][ km-1 ] + 4. * v.x[ i ][ j ][ km-2 ] - v.x[ i ][ j ][ km-3 ] ) / ( 2. * dphi );
			dwdphi = ( - 3. * w.x[ i ][ j ][ km-1 ] + 4. * w.x[ i ][ j ][ km-2 ] - w.x[ i ][ j ][ km-3 ] ) / ( 2. * dphi );

// 2nd order derivative for pressure and velocity components
			d2udr2 = ( u.x[ i ][ j ][ km-1 ] - 2. * u.x[ i ][ j ][ km-1 ] + u.x[ i ][ j ][ km-1 ] ) / dr2;
			d2vdr2 = ( v.x[ i ][ j ][ km-1 ] - 2. * v.x[ i ][ j ][ km-1 ] + v.x[ i ][ j ][ km-1 ] ) / dr2;
			d2wdr2 = ( w.x[ i ][ j ][ km-1 ] - 2. * w.x[ i ][ j ][ km-1 ] + w.x[ i ][ j ][ km-1 ] ) / dr2;

			d2udthe2 = ( u.x[ i ][ j+1 ][ km-1 ] - 2. * u.x[ i ][ j ][ km-1 ] + u.x[ i ][ j-1 ][ km-1 ] ) / dthe2;
			d2vdthe2 = ( v.x[ i ][ j+1 ][ km-1 ] - 2. * v.x[ i ][ j ][ km-1 ] + v.x[ i ][ j-1 ][ km-1 ] ) / dthe2;
			d2wdthe2 = ( w.x[ i ][ j+1 ][ km-1 ] - 2. * w.x[ i ][ j ][ km-1 ] + w.x[ i ][ j-1 ][ km-1 ] ) / dthe2;

			d2udphi2 = ( u.x[ i ][ j ][ km-1 ] - 2. * u.x[ i ][ j ][ km-2 ] + u.x[ i ][ j ][ km-3 ] ) / dphi2;
			d2vdphi2 = ( v.x[ i ][ j ][ km-1 ] - 2. * v.x[ i ][ j ][ km-2 ] + v.x[ i ][ j ][ km-3 ] ) / dphi2;
			d2wdphi2 = ( w.x[ i ][ j ][ km-1 ] - 2. * w.x[ i ][ j ][ 2 ] + w.x[ i ][ j ][ km-3 ] ) / dphi2;

// Coriolis and centrifugal terms in the energy and momentum equations
			RS_Coriolis_Momentum_rad = + coriolis * 2. * omega * sinthe * w.x[ i ][ j ][ km-1 ];
			RS_Coriolis_Momentum_the = + coriolis * 2. * omega * costhe * w.x[ i ][ j ][ km-1 ];
			RS_Coriolis_Momentum_phi = - coriolis * ( 2. * omega * sinthe * u.x[ i ][ j ][ km-1 ] + 2. * omega * costhe * v.x[ i ][ j ][ km-1 ] );

			RS_Centrifugal_Momentum_rad = + centrifugal * rad.z[ i ] * pow ( ( omega * sinthe ), 2 );
			RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ i ] * sinthe * costhe * pow ( ( omega ), 2 );

			if ( t.x[ i ][ j ][ km-1 ] < ta ) 
			{
				tn.x[ i ][ j ][ km-1 ] = ta;
				cn.x[ i ][ j ][ km-1 ] = ca;
				c_Boussinesq = ( ( ( ta * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}
			else
			{
				c_Boussinesq = ( ( ( tn.x[ i ][ j ][ km-1 ] * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}

			RS_buoyancy_Momentum = - buoyancy * .5 * g * ( c.x[ i ][ j ][ km-1 ] - c_Boussinesq ) / c_Boussinesq;		// bouyancy based on water vapour 


// Right Hand Side of the Navier-Stokes equations at the boundaries for the pressure Poisson equation
			aux_u.x[ i ][ j ][ km-1 ] = - ( u.x[ i ][ j ][ km-1 ] * dudr + v.x[ i ][ j ][ km-1 ] * dudthe / rm + w.x[ i ][ j ][ km-1 ] * dudphi / rmsinthe )
												+ ( d2udr2 + d2udthe2 / rm2 + 4. * dudr / rm + d2udphi2 / rm2sinthe2 ) / re
												+ RS_buoyancy_Momentum
												+ RS_Coriolis_Momentum_rad + RS_Centrifugal_Momentum_rad;

			aux_v.x[ i ][ j ][ km-1 ] = - ( u.x[ i ][ j ][ km-1 ] * dvdr + v.x[ i ][ j ][ km-1 ] * dvdthe / rm + w.x[ i ][ j ][ km-1 ] * dvdphi / rmsinthe ) +
												+ ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
												- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * v.x[ i ][ j ][ km-1 ] / rm + d2vdphi2 / rm2sinthe2
												+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
												+ RS_Coriolis_Momentum_the + RS_Centrifugal_Momentum_the
												- v.x[ i ][ j ][ km-1 ] * k_Force / dthe2;

			aux_w.x[ i ][ j ][ km-1 ] = - ( u.x[ i ][ j ][ km-1 ] * dwdr + v.x[ i ][ j ][ km-1 ] * dwdthe / rm + w.x[ i ][ j ][ km-1 ] * dwdphi / rmsinthe ) +
												+ ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
												- ( 1. + costhe * costhe / ( rm * sinthe2 ) ) * w.x[ i ][ j ][ km-1 ] / rm + d2wdphi2 / rm2sinthe2
												+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
												+ RS_Coriolis_Momentum_phi
												- w.x[ i ][ j ][ km-1 ] * k_Force / dphi2;
			}
	}



//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    level j = 0    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	for ( int i = 1; i < im-1; i++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
// collection of coefficients
			rm = rad.z[ i ];
			rm2 = rm * rm;
			sinthe = sin( the.z[ 0 ] );
			sinthe2 = sinthe * sinthe;
			costhe = cos( the.z[ 0 ] );
			cotthe = cos( the.z[ 0 ] ) / sin( the.z[ 0 ] );
			rmsinthe = rm * sinthe;
			rm2sinthe = rm2 * sinthe;
			rm2sinthe2 = rm2 * sinthe2;

// 1st order derivative for pressure and velocity components
			dudr = ( u.x[ i+1 ][ 0 ][ k ] - u.x[ i-1 ][ 0 ][ k ] ) / ( 2. * dr );
			dvdr = ( v.x[ i+1 ][ 0 ][ k ] - v.x[ i-1 ][ 0 ][ k ] ) / ( 2. * dr );
			dwdr = ( w.x[ i+1 ][ 0 ][ k ] - w.x[ i-1 ][ 0 ][ k ] ) / ( 2. * dr );

			dudthe = ( - 3. * u.x[ i ][ 0 ][ k ] + 4. * u.x[ i ][ 1 ][ k ] - u.x[ i ][ 2 ][ k ] ) / ( 2. * dthe );
			dvdthe = ( - 3. * v.x[ i ][ 0 ][ k ] + 4. * v.x[ i ][ 1 ][ k ] - v.x[ i ][ 2 ][ k ] ) / ( 2. * dthe );
			dwdthe = ( - 3. * w.x[ i ][ 0 ][ k ] + 4. * w.x[ i ][ 1 ][ k ] - w.x[ i ][ 2 ][ k ] ) / ( 2. * dthe );

			dudphi = ( u.x[ i ][ 0 ][ k+1 ] - u.x[ i ][ 0 ][ k-1 ] ) / ( 2. * dphi );
			dvdphi = ( v.x[ i ][ 0 ][ k+1 ] - v.x[ i ][ 0 ][ k-1 ] ) / ( 2. * dphi );
			dwdphi = ( w.x[ i ][ 0 ][ k+1 ] - w.x[ i ][ 0 ][ k-1 ] ) / ( 2. * dphi );

// 2nd order derivative for pressure and velocity components
			d2udr2 = ( u.x[ i+1 ][ 0 ][ k ] - 2. * u.x[ i ][ 0 ][ k ] + u.x[ i-1 ][ 0 ][ k ] ) / dr2;
			d2vdr2 = ( v.x[ i+1 ][ 0 ][ k ] - 2. * v.x[ i ][ 0 ][ k ] + v.x[ i-1 ][ 0 ][ k ] ) / dr2;
			d2wdr2 = ( w.x[ i+1 ][ 0 ][ k ] - 2. * w.x[ i ][ 0 ][ k ] + w.x[ i-1 ][ 0 ][ k ] ) / dr2;

			d2udthe2 = ( u.x[ i ][ 0 ][ k ] - 2. * u.x[ i ][ 1 ][ k ] + u.x[ i ][ 2 ][ k ] ) / dthe2;
			d2vdthe2 = ( v.x[ i ][ 0 ][ k ] - 2. * v.x[ i ][ 1 ][ k ] + v.x[ i ][ 2 ][ k ] ) / dthe2;
			d2wdthe2 = ( w.x[ i ][ 0 ][ k ] - 2. * w.x[ i ][ 1 ][ k ] + w.x[ i ][ 2 ][ k ] ) / dthe2;

			d2udphi2 = ( u.x[ i ][ 0 ][ k+1 ] - 2. * u.x[ i ][ 0 ][ k ] + u.x[ i ][ 0 ][ k-1 ] ) / dphi2;
			d2vdphi2 = ( v.x[ i ][ 0 ][ k+1 ] - 2. * v.x[ i ][ 0 ][ k ] + v.x[ i ][ 0 ][ k-1 ] ) / dphi2;
			d2wdphi2 = ( w.x[ i ][ 0 ][ k+1 ] - 2. * w.x[ i ][ 0 ][ k ] + w.x[ i ][ 0 ][ k-1 ] ) / dphi2;

// Coriolis and centrifugal terms in the energy and momentum equations
			RS_Centrifugal_Momentum_rad = + centrifugal * rad.z[ i ] * pow ( ( omega * sinthe ), 2 );
			RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ i ] * sinthe * costhe * pow ( ( omega ), 2 );

			if ( t.x[ i ][ 0 ][ k ] < ta ) 
			{
				tn.x[ i ][ 0 ][ k ] = ta;
				cn.x[ i ][ 0 ][ k ] = ca;
				c_Boussinesq = ( ( ( ta * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}
			else
			{
				c_Boussinesq = ( ( ( tn.x[ i ][ 0 ][ k ] * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}

			RS_buoyancy_Momentum = - buoyancy * .5 * g * ( c.x[ i ][ 0 ][ k ] - c_Boussinesq ) / c_Boussinesq;		// bouyancy based on water vapour 


// Right Hand Side of the Navier-Stokes equations at the boundaries for the pressure Poisson equation
			aux_u.x[ i ][ 0 ][ k ] = - u.x[ i ][ 0 ][ k ] * dudr
												+ ( d2udr2 + d2udthe2 / rm2 + 4. * dudr / rm + d2udphi2 / rm2sinthe2 ) / re
												+ RS_buoyancy_Momentum
												+ RS_Centrifugal_Momentum_rad;

			aux_v.x[ i ][ 0 ][ k ] = 0.;

/*
			aux_v.x[ i ][ 0 ][ k ] = - u.x[ i ][ 0 ][ k ] * dvdr
												+ ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
												+ d2vdphi2 / rm2sinthe2
												+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
												+ RS_Centrifugal_Momentum_the
												- v.x[ i ][ 0 ][ k ] * k_Force / dthe2;
*/
			aux_w.x[ i ][ 0 ][ k ] = 0.;
/*
			aux_w.x[ i ][ 0 ][ k ] = - u.x[ i ][ 0 ][ k ] * dwdr
												+ ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
												+ d2wdphi2 / rm2sinthe2
												+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
												- w.x[ i ][ 0 ][ k ] * k_Force / dphi2;
*/


//	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    level j = jm-1    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// collection of coefficients
			rm = rad.z[ i ];
			rm2 = rm * rm;
			sinthe = sin( the.z[ jm-1 ] );
			sinthe2 = sinthe * sinthe;
			costhe = cos( the.z[ jm-1 ] );
			cotthe = cos( the.z[ jm-1 ] ) / sin( the.z[ jm-1 ] );
			rmsinthe = rm * sinthe;
			rm2sinthe = rm2 * sinthe;
			rm2sinthe2 = rm2 * sinthe2;

// 1st order derivative for pressure and velocity components
			dudr = ( u.x[ i+1 ][ jm-1 ][ k ] - u.x[ i-1 ][ jm-1 ][ k ] ) / ( 2. * dr );
			dvdr = ( v.x[ i+1 ][ jm-1 ][ k ] - v.x[ i-1 ][ jm-1 ][ k ] ) / ( 2. * dr );
			dwdr = ( w.x[ i+1 ][ jm-1 ][ k ] - w.x[ i-1 ][ jm-1 ][ k ] ) / ( 2. * dr );

			dudthe = ( - 3. * u.x[ i ][ jm-1 ][ k ] + 4. * u.x[ i ][ 1 ][ k ] - u.x[ i ][ 2 ][ k ] ) / ( 2. * dthe );
			dvdthe = ( - 3. * v.x[ i ][ jm-1 ][ k ] + 4. * v.x[ i ][ 1 ][ k ] - v.x[ i ][ 2 ][ k ] ) / ( 2. * dthe );
			dwdthe = ( - 3. * w.x[ i ][ jm-1 ][ k ] + 4. * w.x[ i ][ 1 ][ k ] - w.x[ i ][ 2 ][ k ] ) / ( 2. * dthe );

			dudphi = ( u.x[ i ][ jm-1 ][ k+1 ] - u.x[ i ][ jm-1 ][ k-1 ] ) / ( 2. * dphi );
			dvdphi = ( v.x[ i ][ jm-1 ][ k+1 ] - v.x[ i ][ jm-1 ][ k-1 ] ) / ( 2. * dphi );
			dwdphi = ( w.x[ i ][ jm-1 ][ k+1 ] - w.x[ i ][ jm-1 ][ k-1 ] ) / ( 2. * dphi );

// 2nd order derivative for pressure and velocity components
			d2udr2 = ( u.x[ i+1 ][ jm-1 ][ k ] - 2. * u.x[ i ][ jm-1 ][ k ] + u.x[ i-1 ][ jm-1 ][ k ] ) / dr2;
			d2vdr2 = ( v.x[ i+1 ][ jm-1 ][ k ] - 2. * v.x[ i ][ jm-1 ][ k ] + v.x[ i-1 ][ jm-1 ][ k ] ) / dr2;
			d2wdr2 = ( w.x[ i+1 ][ jm-1 ][ k ] - 2. * w.x[ i ][ jm-1 ][ k ] + w.x[ i-1 ][ jm-1 ][ k ] ) / dr2;

			d2udthe2 = ( u.x[ i ][ jm-1 ][ k ] - 2. * u.x[ i ][ 1 ][ k ] + u.x[ i ][ 2 ][ k ] ) / dthe2;
			d2vdthe2 = ( v.x[ i ][ jm-1 ][ k ] - 2. * v.x[ i ][ 1 ][ k ] + v.x[ i ][ 2 ][ k ] ) / dthe2;
			d2wdthe2 = ( w.x[ i ][ jm-1 ][ k ] - 2. * w.x[ i ][ 1 ][ k ] + w.x[ i ][ 2 ][ k ] ) / dthe2;

			d2udphi2 = ( u.x[ i ][ jm-1 ][ k+1 ] - 2. * u.x[ i ][ jm-1 ][ k ] + u.x[ i ][ jm-1 ][ k-1 ] ) / dphi2;
			d2vdphi2 = ( v.x[ i ][ jm-1 ][ k+1 ] - 2. * v.x[ i ][ jm-1 ][ k ] + v.x[ i ][ jm-1 ][ k-1 ] ) / dphi2;
			d2wdphi2 = ( w.x[ i ][ jm-1 ][ k+1 ] - 2. * w.x[ i ][ jm-1 ][ k ] + w.x[ i ][ jm-1 ][ k-1 ] ) / dphi2;

// Coriolis and centrifugal terms in the energy and momentum equations
			RS_Centrifugal_Momentum_rad = + centrifugal * rad.z[ i ] * pow ( ( omega * sinthe ), 2 );
			RS_Centrifugal_Momentum_the = + centrifugal * rad.z[ i ] * sinthe * costhe * pow ( ( omega ), 2 );

			if ( t.x[ i ][ jm-1 ][ k ] < ta ) 
			{
				tn.x[ i ][ jm-1 ][ k ] = ta;
				cn.x[ i ][ jm-1 ][ k ] = ca;
				c_Boussinesq = ( ( ( ta * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}
			else
			{
				c_Boussinesq = ( ( ( tn.x[ i ][ jm-1 ][ k ] * t_0 - t_0 ) + 346. ) / 10. ) / c_0;
			}

			RS_buoyancy_Momentum = - buoyancy * .5 * g * ( c.x[ i ][ jm-1 ][ k ] - c_Boussinesq ) / c_Boussinesq;		// bouyancy based on water vapour 


// Right Hand Side of the Navier-Stokes equations at the boundaries for the pressure Poisson equation
			aux_u.x[ i ][ jm-1 ][ k ] = - u.x[ i ][ jm-1 ][ k ] * dudr
												+ ( d2udr2 + d2udthe2 / rm2 + 4. * dudr / rm + d2udphi2 / rm2sinthe2 ) / re
												+ RS_buoyancy_Momentum
												+ RS_Centrifugal_Momentum_rad;

			aux_v.x[ i ][ jm-1 ][ k ] = - u.x[ i ][ jm-1 ][ k ] * dvdr
												+ ( d2vdr2 + dvdr * 2. / rm + d2vdthe2 / rm2 + dvdthe / rm2sinthe * costhe
												+ d2vdphi2 / rm2sinthe2
												+ 2. * dudthe / rm2 - dwdphi * 2. * costhe / rm2sinthe2 ) / re
												+ RS_Centrifugal_Momentum_the
												- v.x[ i ][ jm-1 ][ k ] * k_Force / dthe2;

			aux_w.x[ i ][ jm-1 ][ k ] = - u.x[ i ][ jm-1 ][ k ] * dwdr
												+ ( d2wdr2 + dwdr * 2. / rm + d2wdthe2 / rm2 + dwdthe / rm2sinthe  * costhe
												+ d2wdphi2 / rm2sinthe2
												+ 2. * dudphi / rm2sinthe + dvdphi * 2. * costhe / rm2sinthe2 ) / re
												- w.x[ i ][ jm-1 ][ k ] * k_Force / dphi2;
		}
	}

}

