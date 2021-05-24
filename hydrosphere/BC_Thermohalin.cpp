/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
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
#include <Utils.h>
#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "cHydrosphereModel.h"
#include "AtomMath.h"

using namespace std;
using namespace AtomUtils;

void cHydrosphereModel::PresStat_SaltWaterDens(){
cout << endl << "      PresStat_SaltWaterDens" << endl;
// hydrostatic pressure, equations of state for water and salt water density
// as functions of salinity, temperature and hydrostatic pressure
    double t_Celsius_0 = 0.;
    double t_Celsius_1 = 0.;
    double p_km = 0.;
    double C_p = 0.;
    double beta_p =  0.;
    double alfa_t_p =  0.;
    double gamma_t_p =  0.;
    double E_water = 2.15e9;  // given in N/m²
    double beta_water = 8.8e-5;  // given in m³/(m³ * °C)
    double r_air = 1.2041;  // given in kg/m³
    double R_Air = 287.1;  // given in J/(kg*K)
    double r_0_saltwater = 1027.;  // in kg/m³
// hydrostatic pressure, water and salt water density at the surface
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            p_stat.x[im-1][j][k] = .01 *(r_air * R_Air * t.x[im-1][j][k] * t_0)/1000.;
            // given in bar, isochoric approach, constant air density at the surface
            r_water.x[im-1][j][k] = r_0_water;  // given in kg/m³
            t_Celsius_1 = t.x[im-1][j][k] * t_0 - t_0;
            p_km = 0.;
            C_p = 999.83;
            beta_p = .808;
            alfa_t_p = .0708 *(1. + .068 * t_Celsius_1);
            gamma_t_p = .003 *(1. - .012 * t_Celsius_1);
            r_salt_water.x[im-1][j][k] = C_p + beta_p * c.x[im-1][j][k] * c_0  //in kg/m³ approximation by Gill (saltwater.pdf)
                - alfa_t_p * t_Celsius_1 - gamma_t_p
                * ( 1. - c.x[im-1][j][k]) * c_0 * t_Celsius_1;
        }
    }
// hydrostatic pressure, water and salt water density in the flow field
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = im-2; i >= 0; i--){
                double d_i = (double)(im-1-i);
                t_Celsius_1 = t.x[i][j][k] * t_0 - t_0;
                t_Celsius_0 = t.x[i+1][j][k] * t_0 - t_0;
//              p_stat.x[i][j][k] = r_0_water * g * d_i *(L_hyd /(double)(im-1))
//                   100000. + p_0/1000.;                // hydrostatic pressure in bar
                p_stat.x[i][j][k] = r_water.x[i+1][j][k] * g * d_i 
                    * (L_hyd/(double)(im-1))/100000. + p_0/1000.;           // hydrostatic pressure in bar
                r_water.x[i][j][k] = r_water.x[i+1][j][k]/(1. + beta_water  //in kg/m³
                    * (t_Celsius_1 - t_Celsius_0))/(1. - (p_stat.x[i][j][k]
                    - p_stat.x[i+1][j][k])/E_water * 1e5);
                p_km  = (double)(im-1-i) * (L_hyd/(double)(im-1))/1000.;  // depth in km
                C_p = 999.83 + 5.053 * p_km - .048 * p_km * p_km;
                beta_p = .808 - .0085* p_km;
                alfa_t_p = .0708 *(1. + .351 * p_km 
                    + .068 * (1. - .0683 * p_km) * t_Celsius_1);
                gamma_t_p = .003 *(1. - .059 * p_km 
                    - .012 * (1. - .064 * p_km) * t_Celsius_1);
                r_salt_water.x[i][j][k] = C_p + beta_p * c.x[i][j][k] * c_0  //in kg/m³ approximation by Gill (saltwater.pdf)
                    - alfa_t_p * t_Celsius_1 - gamma_t_p
                    * ( 1. - c.x[i][j][k]) * c_0 * t_Celsius_1;
            }
            for(int i = im-1; i >= 0; i--){
                if(is_land(h, i, j, k))  r_salt_water.x[i][j][k] = r_0_saltwater;  // arbitrary choosen for paraview plotting
            }
        }
    }
cout << "      PresStat_SaltWaterDens ended" << endl;
}
/*
*
*/
void cHydrosphereModel::Value_Limitation_Hyd(){
cout << endl << "      Value_Limitation_Hyd" << endl;
// the limiting values depend on local singular behaviour
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
/*
            if(EkmanPumping.y[j][k] >= .01)  EkmanPumping.y[j][k] = .01; // in m/s
            if(EkmanPumping.y[j][k] < .01)  EkmanPumping.y[j][k] = -.01; // in m/s
            if(is_land(h, im-1, j, k))  EkmanPumping.y[j][k] = 0.;
            if(Upwelling.y[j][k] >= .01)  Upwelling.y[j][k] = .01; // in m/s
            else  Upwelling.y[j][k] = 0.;
            if(Downwelling.y[j][k] < -.01)  Downwelling.y[j][k] = -.01; // in m/s
            else  Downwelling.y[j][k] = 0.;
*/
            for(int i = 0; i < im; i++){
/*
                if(u.x[i][j][k] >= 0.01/u_0)  u.x[i][j][k] = 0.01/u_0; // non-dimensional
                if(u.x[i][j][k] <= - 0.01/u_0)  u.x[i][j][k] = - 0.01/u_0; // non-dimensional
*/
                if(v.x[i][j][k] >= .125/u_0)  v.x[i][j][k] = .125/u_0;
                if(v.x[i][j][k] <= - .125/u_0)  v.x[i][j][k] = - .125/u_0;
                if(w.x[i][j][k] >= .125/u_0)  w.x[i][j][k] = .125/u_0;
                if(w.x[i][j][k] <= - .125/u_0)  w.x[i][j][k] = - .125/u_0;

                if(t.x[i][j][k] >= 1.147)  t.x[i][j][k] = 1.147; //40.15 °C
                if(t.x[i][j][k] <= 0.9927)  t.x[i][j][k] = 0.9927;// -1.0 °C
                if(c.x[i][j][k] * c_0 >= 50.)  c.x[i][j][k] = 1.4451;  // 50.0 psu
                if(c.x[i][j][k] * c_0 <= 32.)  c.x[i][j][k] = 0.9249;      // 32.0 psu
            }
        }
    }
cout << "      Value_Limitation_Hyd ended" << endl;
}
/*
*
*/
void cHydrosphereModel::SalinityEvaporation(){
cout << endl << "      SalinityEvaporation" << endl;
// preparations for salinity increase due to evaporation minus precipitation
// procedure given in Rui Xin Huang, Ocean Circulation, p. 165
    double coeff_salinity = 1.1574e-8 * L_hyd/(r_0_water * u_0 * c_0);  // 1.1574-8 is the conversion from (Evap-Prec) in mm/d to m/s
    double evap_precip = 0.;
    double sal_flux = 0.;
    double step = L_hyd/(double)(im-1);  // in m, 200/40 = 5m
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            c_fix.y[j][k] = c.x[im-1][j][k];
                evap_precip = coeff_salinity 
//                    * (Evaporation_Dalton.y[j][k] - Precipitation.y[j][k]);  // in m/s
                    * (Evaporation_Penman.y[j][k] - Precipitation.y[j][k]);  // in m/s
//            sal_flux = - (- 3. * c.x[im-1][j][k] + 4. * c.x[im-2][j][k]  // 1. order derivative, 2. order accurate
//                - c.x[im-3][j][k])/(2. * step) 
//                * (1. - 2. * c.x[im-1][j][k]);
            sal_flux = - (c.x[im-1][j][k] - c.x[im-2][j][k])/step  // 1. order derivative, 1. order accurate
                * (1. - 2. * c.x[im-1][j][k]);
            salinity_evaporation.y[j][k] = r_salt_water.x[im-1][j][k] 
                * sal_flux * evap_precip;
            if(is_land(h, im-1, j, k))  
                salinity_evaporation.y[j][k] = 0.;
            c.x[im-1][j][k] = c_fix.y[j][k] 
                + salinity_evaporation.y[j][k];
/*
            cout.precision(8);
            cout.setf(ios::fixed);
            if((j == 90)&&(k == 180)) cout << endl
                << "  SalinityEvaporation" << endl
                << "   j = " << j << "   k = " << k << endl
                << "   coeff_salinity = " << coeff_salinity << endl
                << "   r_0_water = " << r_0_water 
                << "   r_salt_water = " << r_salt_water.x[im-1][j][k] << endl
                << "   evap_precip = " << evap_precip/coeff_salinity
//                << "   evap = " << Evaporation_Dalton.y[j][k] 
                << "   evap = " << Evaporation_Penman.y[j][k] 
                << "   prec = " << Precipitation.y[j][k] << endl
                << "   sal_flux = " << sal_flux * c_0  << endl
                << "   salinity_evaporation = " << salinity_evaporation.y[j][k] * c_0 << endl
                << "   c_fix = " << c_fix.y[j][k] * c_0 
                << "   c = " << c.x[im-1][j][k] * c_0 
                << "   c_0 = " << c_0 << endl;
*/
/*
                int i_trans = 37; // 15 m vertical extension of ocean surface evaporation, arbitrary 
//                int i_trans = 20; // 15 m vertical extension of ocean surface evaporation, arbitrary 
                for(int i = i_trans; i <= im-1; i++){
//                    if(i <= im-1){
//                        if(i > 0){
//                            double x = get_layer_height(i) 
//                               /get_layer_height(i_trans); 
                            double x = (double)i/(double)i_trans; 
//                            c.x[i][j][k] = parabola_interp(c.x[i_trans][j][k], 
//                                c.x[im-1][j][k], x); // parabolic transition
                            c.x[i][j][k] = parabola_interp(c.x[i_trans][j][k], 
                                c.x[im-1][j][k], x); // parabolic transition

                            c.x[i][j][k] = (c.x[i_trans][j][k] - c.x[0][j][k]) 
                                * (get_layer_height(i)/get_layer_height(i_trans)) 
                                + c.x[im-1][j][k]; // linear transition

//                        }
//                    }
                } // end i
*/
        }
    }
cout << "      SalinityEvaporation ended" << endl;
}
/*
*
*/
