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
#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"

#ifndef _RHS_HYDROSPHERE_
#define _RHS_HYDROSPHERE_

using namespace std;

class RHS_Hydrosphere
{
	private:
		int im, jm, km;
		
        double dt, dr, dthe, dphi, r0;
		double re, pr, ec, sc, m_g;
		
		double Buoyancy;

	public:
        RHS_Hydrosphere ( int jm, int km, double dthe, double dphi, double re );

	    RHS_Hydrosphere ( int im, int jm, int km, double r0, double dt, double dr, double dthe, double dphi, double re, 
            double ec, double sc, double g, double pr, double buoyancy );

        ~RHS_Hydrosphere ();

        void RK_RHS_3D_Hydrosphere ( int i, int j, int k, double L_hyd, double g, double cp_w, double u_0, double t_0, 
            double c_0, 
            double r_0_water, double ta, double pa, double ca, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, 
            Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &c, Array &tn, Array &un, Array &vn, Array &wn, 
            Array &p_dynn, Array &cn, Array &rhs_t, Array &rhs_u, Array &rhs_v, Array &rhs_w, Array &rhs_c, 
            Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion, Array &BuoyancyForce_3D, 
            Array &Salt_Balance, Array &p_stat, Array &ro_water, Array &ro_salt_water );

        void RK_RHS_2D_Hydrosphere ( int j, int k, double r_0_water, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &v, Array &w, 
            Array &p_dyn, Array &vn, Array &wn, Array &p_dynn, Array &rhs_v, Array &rhs_w, Array &aux_v, 
            Array &aux_w );

};
#endif
