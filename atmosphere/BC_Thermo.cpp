/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
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

#include "Array.h"
#include "Array_2D.h"
#include "cAtmosphereModel.h"
#include "Utils.h"
#include "AtomMath.h"

using namespace std;
using namespace AtomUtils;

void cAtmosphereModel::RadiationMultiLayer(){
    cout << endl << "      RadiationMultiLayer" << endl;
    // class element for the computation of the radiation and the temperature distribution
    // computation of the local temperature based on short and long wave radiation
    // multi layer radiation model
//    temp_tropopause = std::vector<double>(jm, t_tropopause_pole);
    std::vector<double> step(im, 0);
    std::vector<double> alfa(im, 0.0);
    std::vector<double> beta(im, 0.0);
    std::vector<double> AA(im, 0.0);
    std::vector<double> CA(im, 0.0);
    std::vector<std::vector<double> > CC(im, std::vector<double>(im, 0.0));
    short_wave_radiation = std::vector<double>(jm, rad_pole_short);
    radiation_original = std::vector<double>(im, 0.0);
    int j_max = jm-1;
    int j_half = j_max/2;
    double albedo_eff = albedo_pole - albedo_equator;
    // effective temperature, albedo and emissivity/absorptivity for the two layer model
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im-1; i++){
                if(is_ocean_surface(h, i, j, k)
                    ||is_land_surface(h, i, j, k)){
                    albedo.y[j][k] = albedo_eff * parabola((double)j
                        /(double)j_half) + albedo_pole;
                }
            }
        }
    }
// short wave radiation
    double rad_short_eff = rad_pole_short - rad_equator_short;
    for(int j=j_half; j>=0; j--){
        short_wave_radiation[j] = rad_short_eff * parabola((double)j
            /(double)j_half) + rad_pole_short;
    }
    for(int j=j_max; j>j_half; j--){
        short_wave_radiation[j] = short_wave_radiation[j_max-j];
    }
    int i_trop = 0;
    int i_mount = 0;
    for(int j = 0; j < jm; j++){
//        i_trop = m_model-> get_tropopause_layer(j);
        i_trop = im-1;
        for(int k = 0; k < km; k++){
//            i_mount = i_topography[j][k];
            i_mount = 0;
            for(int i = i_mount; i <= i_trop; i++){
                radiation.x[i][j][k] = sigma * pow(t.x[i][j][k] * t_0, 4.0);
                radiation_original[i] = radiation.x[i][j][k];
            }
            for(int i = i_mount; i <= i_trop; i++){
//                double t_u = t.x[i][j][k] * t_0;
//                step[i] = log((get_layer_height(i+1) + 1.0)
//                    /(get_layer_height(i) + 1.0));
                step[i] = get_layer_height(i+1) - get_layer_height(i);
//                step[i] = log(((double)(i+1) + 1.0)/((double)(i) + 1.0));
                // dependency given by Häckel, Meteorologie, p. 205 (law by F. Baur and H. Philips, 1934)
                // epsilon_eff describe the effect in the emissivity computation of other gases like CO2
                // in the original formula this value is 0.594, for reasons of adjustment to the modern atmosphere,
                // this constant became a variable in zonal direction
                // this variable reacts very sensitive and changes the temperature field extremely
                // the second term describes the influence of water vapour only 
                // applicable models described in K. H. Byun and L.-D. Chen: 
                // Total emissivity of CO°2 near earth condition, Journal of Mechanical Science and Technology 27 (10)(2013)3183-3189)
                double P_c = (1e-6 * p_hydro.x[i][j][k] 
                             * co2.x[i][j][k] * co2_0/p_0);  // in atm, partial pressure of CO2
                double u_c = P_c * step[i] * 100.0; // in atm*cm, model by Atwater and Ball
//                u_c = 300.0/0.9869 * P_c * step[i] * 100.0/(t.x[i][j][k] * t_0); 
                        // modell by Yamamoto and Sasamori
                 // influence on the emissivity by carbon dioxcide by the law by Bliss
//                x_c = 10.0 * P_c/(R_co2 * t.x[i][j][k] * t_0) * step[i] * p_0;
                  // coefficient in the Bliss law, 10.0 => conversion from kg/m² to g/cm²
//                eps_co2 = 0.185 * (1.0 - exp (- 50.0 * x_c)); // modell by Bliss
                double eps_co2 = 0.185 * (1.0 - exp(- 0.3919 * pow(u_c, 0.4))); // model by Atwater and Ball
                eps_co2 = 0.5 * eps_co2; // model by Atwater and Ball, 
                eps_co2 = 0.0;  // no influence of CO2 on the emissivity
                // factor .5 fits the Atwater and Ball results to Yamamoto and Sasamori,
                // which fit best to HITRAN band data
                // COSMO water vapour pressure based on local water vapour in hPa
                double e = c.x[i][j][k] * p_hydro.x[i][j][k]/ep;  // in hPa, simplified approach
//                double e = c.x[i][j][k] * p_hydro.x[i][j][k]/(c.x[i][j][k] + ep);  // in hPa, complete equation

//                epsilon.x[i][j][k] = eps_co2 + 0.594 + 0.0416 * sqrt(e); // constant  given in Häckel, p. 205 (F. Baur and H. Philips, 1934), ok
//                epsilon.x[i][j][k] = eps_co2 + 0.52 + 0.065 * sqrt(e); // (Brunt, 1932), leads to a little lower radiation values
                epsilon.x[i][j][k] = eps_co2 + 0.684 + 0.0056 * e; // (Bignami, 1995, -40 <=> 45°C), best choice

//                epsilon.x[i][j][k] = eps_co2 + (1.0 - 0.26 * exp(- 7.77e-4 
//                    * (t_0 - t_u))); // (Idso and Jackson, 1969, -29 <=> 37°C), nearly no difference to the e based laws
//                epsilon.x[i][j][k] = eps_co2 + 9.2e-6 * t_u * t_u; // (Swinbank, 1963, 2 <=> 29°C), tends to lower radiation results tha tht e based laws

//                epsilon.x[i][j][k] = eps_co2 + (0.7 + 5.95e-5 * e
//                    * exp(1500.0/t_u)); // (Idso, 1981, -40 <=> 45°C), nearly no difference to the e based laws
//                epsilon.x[i][j][k] = eps_co2 + 1.24 * pow((e/t_u),1/7); // (Brutsaert, 1975, -40 <=> 45°C), nearly no difference to the e based laws
//                epsilon.x[i][j][k] = 0.98 * pow((e/t_u), 0.0687); // atmospheric emissivity by Vogel/Bliss, produces NANs

//                epsilon.x[i][j][k] = eps_co2 + cloudiness.x[i][j][k]; // introduction of the cloudiness by Xu and Randall, not applicable
/*
                cout.precision(6);
//                if((i == 20)&&(j == 90)&&(k == 180))  cout << endl 
                if((j == 90)&&(k == 180))  cout << endl 
                    << "  Multi Layer Radiation Model preparation   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
                    << "   i = " << i << "  j = " << j << "  k = " << k << endl 
                    << "  step = " << step[i]
                    << "  e = " << e
                    << "  P_c = " << P_c
                    << "  u_c = " << u_c << endl
                    << "  eps_co2 = " << eps_co2
                    << "  cloudiness = " << cloudiness.x[i][j][k]
                    << "  epsilon = " << epsilon.x[i][j][k] << endl
                    << "  c = " << c.x[i][j][k]
                    << "  cloud = " << cloud.x[i][j][k]
                    << "  ice = " << ice.x[i][j][k] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0
                    << "  t_tropopause_equator = " << t_tropopause_equator * t_0 - t_0 << endl
                    << "  radiation_original[i] = " << radiation_original[i]
                    << "  short_wave_radiation[j] = " << short_wave_radiation[j] << endl
                    << "  radiation_mount = " << radiation.x[i_mount][j][k]
                    << "  radiation = " << radiation.x[i][j][k]
                    << "  radiation_trop = " << sigma
                        * pow(t.x[i_trop][j][k] * t_0, 4.0) << endl;
*/
                if(epsilon.x[i][j][k] > 1.)  epsilon.x[i][j][k] = 1.;
            }  // end i
            epsilon_2D.y[j][k] = epsilon.x[0][j][k];
/*
            // inside mountains
            for(int i = i_mount-1; i >= 0; i--){
                if(is_land(h, i, j, k)){
                    epsilon.x[i][j][k] = epsilon.x[i_mount][j][k];
                    radiation.x[i][j][k] = radiation.x[i_mount][j][k];
                    radiation_original[i] = radiation.x[i_mount][j][k];
                }
            } // end i
*/
/*
            // above tropopause
            for(int i = i_trop; i < im; i++){
                epsilon.x[i][j][k] = epsilon.x[i_trop][j][k];
                radiation_original[i_trop] = sigma * pow(t.x[i_trop][j][k] * t_0, 4.0);
                aux_t.x[i][j][k] = t.x[i_trop][j][k];
            } // end i
*/
            AA[i_mount] = radiation.x[i_mount][j][k];// non-dimensional surface radiation
            CC[i_mount][i_mount] = epsilon.x[i_mount][j][k] 
                * radiation.x[i_mount][j][k]; // no absorption of radiation on the surface by water vapour
            for(int i = i_mount+1; i <= i_trop; i++){
                AA[i] = AA[i-1] * (1.0 - epsilon.x[i][j][k]); // transmitted radiation from each layer
                CC[i][i]= epsilon.x[i][j][k] * radiation.x[i][j][k]; // absorbed radiation in each layer
            }
            for(int i = i_mount+2; i <= i_trop; i++){
                CA[i] = 0.0;
                for(int l = 1; l <= i-1; l++){
                    CC[l][i] = CC[l][i-1] * (1.0 - epsilon.x[i][j][k]); // transmitted radiation leaving the actual layer i
                    CA[i] += CC[l][i]; // sum on all transmitted radiations with index l passing the actual layer i
/*
                    cout.precision(6);
                    if((j == 90)&&(k == 180))  cout << endl 
                        << "  Multi Layer Radiation Model AA[i]  BA[i]     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
                        << "   iter_rad = " << iter_rad << endl
                        << "   l = " << l << "  i = " << i << "  j = " << j << "  k = " << k << endl 
                        << "  CC[i][i] = " << CC[i][i]
                        << "  CC[l][i] = " << CC[l][i]
                        << "  CA[i] = " << CA[i] << endl;
*/
                } // end l
            } // end i
/*
            cout.precision(6);
             if((j == 90)&&(k == 180))  cout << endl
                << "  Multi Layer Radiation Model absoption = emmission    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl
                << "   iter_rad = " << iter_rad << endl
                << "   j = " << j << "  k = " << k << endl 
                << "  albedo = " << albedo.y[j][k]
                << "  epsilon[i_mount] = " << epsilon.x[i_mount][j][k]
                << "  epsilon[i_mount+1] = " << epsilon.x[i_mount+1][j][k] << endl
                << "  short_wave_radiation = " << short_wave_radiation[j]
                << "  short_wave_radiation_albedo = " << (1.0 - albedo.y[j][k]) * short_wave_radiation[j] << endl
                << "  radiation_original[i_mount] = " << radiation_original[i_mount]  << endl
                << "  radiation[i_mount] = " << radiation.x[i_mount][j][k]
                << "  radiation[i_mount+1] = " << radiation.x[i_mount+1][j][k] << endl
                << "  AA[i_mount] = " << AA[i_mount]
                << "  AA[i_mount+1] = " << AA[i_mount+1]
                << "  CC[i_mount][i_mount] = " << CC[i_mount][i_mount]
                << "  CC[i_mount+1][i_mount+1] = " << CC[i_mount+1][i_mount+1]
                << "  CA[i_mount+1] = " << CA[i_mount+1] << endl << endl
                << "  t[i_mount] = " << t.x[i_mount][j][k] * t_0 - t_0
                << "  t[i_mount+1] = " << t.x[i_mount+1][j][k] * t_0 - t_0 << endl << endl << endl;
*/
// Thomas algorithm to solve the tridiogonal equation system for the solution of the radiation with a recurrence formula
            double aa, bb, cc, dd;
            for(int i = i_mount; i < i_trop; i++){
                if(i == i_mount){
                    aa = 0.0;
                    bb = - 2.0 * radiation.x[i][j][k];
                    cc = epsilon.x[i+1][j][k] * radiation.x[i+1][j][k];
                    dd = - (1.0 - albedo.y[j][k]) 
                        * short_wave_radiation[j];
                    alfa[i] = - cc/bb;
                    beta[i] = + dd/bb;
                    }
                    if(i == i_mount+1){
                        aa = radiation.x[i-1][j][k];
                        bb = - 2.0 * epsilon.x[i][j][k] * radiation.x[i][j][k];
                        cc = epsilon.x[i+1][j][k] * radiation.x[i+1][j][k];
                        dd = AA[i];
                        alfa[i] = - cc/(bb + aa * alfa[i-1]);
                        beta[i] = + (dd - aa * beta[i-1]) 
                            /(bb + aa * alfa[i-1]);
                    }
                    if(i == i_mount+2){
                        aa = epsilon.x[i-1][j][k] * radiation.x[i-1][j][k];
                        bb = - 2.0 * epsilon.x[i][j][k] * radiation.x[i][j][k];
                        cc = epsilon.x[i+1][j][k] * radiation.x[i+1][j][k];
                        dd = - AA[i-1] + AA[i] + CC[i-1][i];
                        alfa[i] = - cc/(bb + aa * alfa[i-1]);
                        beta[i] = + (dd - aa * beta[i-1]) 
                            /(bb + aa * alfa[i-1]);
                    }
                    if(i > i_mount+2){
                        aa = epsilon.x[i-1][j][k] * radiation.x[i-1][j][k];
                        bb = - 2.0 * epsilon.x[i][j][k] * radiation.x[i][j][k];
                        cc = epsilon.x[i+1][j][k] * radiation.x[i+1][j][k];
                        dd = - AA[i-1] + AA[i] - CA[i-1] + CA[i];
                        alfa[i] = - cc/(bb + aa * alfa[i-1]);
                        beta[i] = + (dd - alfa[i] * beta[i-1]) 
                            /(bb + aa * alfa[i-1]);
                    }
/*
                    cout.precision(6);
                    if((j == 90)&&(k == 180))  cout << endl 
                        << "  Multi Layer Radiation Model    /////////////////////////////////////////////////////" << endl
                        << "   i = " << i << "  j = " << j << "  k = " << k << endl 
                        << "  epsilon[i] = " << epsilon.x[i][j][k]
                        << "  epsilon[i+1] = " << epsilon.x[i+1][j][k] << endl
                        << "  radiation_original[i] = " << radiation_original[i]
                        << "  radiation_original[i+1] = " << radiation_original[i+1] << endl
                        << "  radiation[i] = " << radiation.x[i][j][k]
                        << "  radiation[i+1] = " << radiation.x[i+1][j][k] << endl
                        << "  t[i] = " << t.x[i][j][k] * t_0 - t_0
                        << "  t[i+1] = " << t.x[i+1][j][k] * t_0 - t_0 << endl
                        << "  AA[i] = " << AA[i]
                        << "  AA[i+1] = " << AA[i+1]
                        << "  CC[i][i] = " << CC[i][i]
                        << "  CC[i+1][i+1] = " << CC[i+1][i+1]
                        << "  CA[i] = " << CA[i]
                        << "  CA[i+1] = " << CA[i+1]
                        << "  CA[i+1] - CA[i] = " << CA[i+1] - CA[i] << endl
                        << "  aa = " << aa
                        << "  bb = " << bb
                        << "  cc = " << cc
                        << "  dd = " << dd << endl
                        << "  alfa[i] = " << alfa[i]
                        << "  beta[i] = " << beta[i] << endl << endl;
*/
            } // end i
            radiation.x[i_trop][j][k] = epsilon.x[i_trop][j][k] 
                * radiation.x[i_trop][j][k]/( - CA[i_trop-1] + CA[i_trop]);
//            for(int i = i_trop-1; i >= i_mount; i--){
            for(int i = i_trop-1; i >= 0; i--){
            // Thomas algorithm, recurrence formula
                radiation.x[i][j][k] = alfa[i] * radiation.x[i+1][j][k] 
                    + beta[i];
/*
                cout.precision(6);
                if((j == 90)&&(k == 180))  cout << endl 
                    << "  Multi Layer Radiation Model recurrence" << endl
                    << "   i = " << i << "  j = " << j << "  k = " << k << endl 
                    << "  c = " << c.x[i][j][k]
                    << "  cloud = " << cloud.x[i][j][k]
                    << "  ice = " << ice.x[i][j][k] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl
                    << "  epsilon = " << epsilon.x[i][j][k] << endl
                    << "  albedo = " << albedo.y[j][k]
                    << "  short_wave_radiation = " << short_wave_radiation[j]
                    << "  short_wave_radiation_albedo = " << (1.0 - albedo.y[j][k]) * short_wave_radiation[j] << endl
                    << "  alfa[i] = " << alfa[i]
                    << "  beta[i] = " << beta[i] << endl
                    << "  radiation_mount = " << radiation.x[i_mount][j][k]
                    << "  radiation = " << radiation.x[i][j][k]
                    << "  radiation_trop = " << radiation.x[i_trop][j][k] << endl
                    << "  rad_original[i_mount] = " << radiation_original[i_mount]
                    << "  rad_original[i] = " << radiation_original[i]
                    << "  rad_original[i_trop] = " << radiation_original[i_trop] << endl
                    << "  AA[i] = " << AA[i]
                    << "  AA[i+1] = " << AA[i+1]
                    << "  CC[i][i] = " << CC[i][i]
                    << "  CC[i+1][i+1] = " << CC[i+1][i+1]
                    << "  CA[i] = " << CA[i]
                    << "  CA[i+1] = " << CA[i+1]
                    << "  CA[i+1] - CA[i] = " << CA[i+1] - CA[i] << endl << endl;
*/
            } // end i
//            for(int i = i_mount; i <= i_trop; i++){
            for(int i = 0; i <= i_trop; i++){
                radiation.x[i][j][k] = radiation_original[i] + radiation.x[i][j][k];
                t.x[i][j][k] = pow(radiation.x[i][j][k]/sigma, 0.25)/t_0;
/*
                cout.precision(6);
                 if((j == 90)&&(k == 180))  cout << endl 
                    << "  Multi Layer Radiation Model temperature iteration    §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl
                    << "   i = " << i << "  j = " << j << "  k = " << k << endl 
                    << "  epsilon = " << epsilon.x[i][j][k] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl
                    << "  radiation = " << radiation.x[i][j][k]
                    << "  radiation_original[i] = " << radiation_original[i] << endl << endl;
*/
            } // end i
/*
            for(int i = i_trop; i < im; i++){ // above tropopause
                t.x[i][j][k] = temp_tropopause[j];
                radiation.x[i][j][k] = sigma
                    * pow(t.x[i][j][k] * t_0, 4.0);
            } // end i
*/
        } // end k
    } // end j
/*
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i <= i_mount)){
                    radiation.x[i][j][k] = radiation.x[i_mount][j][k];
                    t.x[i][j][k] = t.x[i_mount][j][k];
                }
            }
        }
    }
*/
    cout << "      RadiationMultiLayer ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::PressureDensity(){
    cout << endl << "      PressureDensity" << endl;
// hydrostatic pressure is understood by a zero velocity field
    double R_W_R_A = R_WaterVapour/R_Air;
    double beta = 42.0; // in K, COSMO
    double t_u = 0.0;
    double height = 0.0;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
//            int i_mount = i_topography[j][k];
            int i_mount = 0;
            height = get_layer_height(i_mount);
            t_u = temp_reconst.y[j][k] + t_0;
            p_hydro.x[0][j][k] = 1e-2 * (r_air * R_Air * t_u);  // reference static pressure given in hPa, by gas equation
            r_dry.x[0][j][k] = 1e2 * p_hydro.x[0][j][k]
                /(R_Air * t_u); // in kg/m3, COSMO
            for(int i = i_mount; i < im; i++){
                height = get_layer_height(i);
                p_hydro.x[i][j][k] = p_hydro.x[0][j][k] // potential static pressure expected at height by barometric height formula,
                    * exp(- t_u/beta * (1.0 - sqrt(1.0 
                    - (2.0 * beta * g * height)/(R_Air * t_u * t_u)))); // COSMO
                r_dry.x[i][j][k] = 1e2 * p_hydro.x[i][j][k]/(R_Air * t_u);
                r_humid.x[i][j][k] = 1e2 * p_hydro.x[i][j][k]
                    /(R_Air * (1. + (R_W_R_A - 1.0) * c.x[i][j][k] 
                    - cloud.x[i][j][k] - ice.x[i][j][k]) * t_u); 

//                p_hydro.x[i][j][k] = 1e-2 * r_dry.x[i][j][k] * R_Air 
//                    * t.x[i][j][k] * t_0 * (1.0 + 0.608 * c.x[i][j][k]);
/*
                cout.precision(6);
                cout.setf(ios::fixed);
                if((j == 90) && (k == 180)) cout << endl 
                    << "  PressureDensity" << endl
                    << "  i = " << i 
                    << "  height = " << get_layer_height(i) << endl
                    << "  p_stat_sl = " << p_hydro.x[0][j][k]
                    << "  p_hydro = " << p_hydro.x[i][j][k]
                    << "  r_dry = " << r_dry.x[i][j][k]
                    << "  r_humid = " << r_humid.x[i][j][k]
                    << "  t_u = " << t_u - t_0
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl << endl;
*/
            }
        }
    }
/*
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            p_stat_landscape.y[j][k] = p_hydro.x[i_mount][j][k]; // in hPa
            r_dry_landscape.y[j][k] = r_dry.x[i_mount][j][k]; // in kg/m³
            r_humid_landscape.y[j][k] = r_humid.x[i_mount][j][k]; // in kg/m³

            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i <= i_mount)){
                    p_hydro.x[i][j][k] = p_hydro.x[i_mount][j][k];
                    r_dry.x[i][j][k] = r_dry.x[i_mount][j][k];
                    r_humid.x[i][j][k] = r_humid.x[i_mount][j][k];
                }
            }

        }
    }
*/
    cout << "      PressureDensity ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::LatentSensibleHeat(){
    cout << endl << "      LatentSensibleHeat" << endl;
    float Q_Latent_Ice = 0.0; 
    float coeff_S = cp_l * r_air * t_0; // coefficient for Q_Sensible
    float coeff_L = r_water_vapour * c_0; // coefficient for Q_Latent
    float coeff_lat = 5.4e4; // diffusion coefficient for water vapour in air, coefficient meets the ratio 27 % to 5 % of latent to sensible heat
    float coeff_sen = 2.7; // diffusion coefficient for sensible heat transfer in air
    float step, step_p, step_m;
    for(int j = 1; j < jm-1; j++){
    // water vapour can condensate/evaporate und sublimate/vaporize
    // water vapour turns to or developes from water or ice
        for(int k = 1; k < km-1; k++){
            for(int i = 1; i < im-1; i++){
                step_p = (get_layer_height(i+1) - get_layer_height(i));
                step_m = (get_layer_height(i) - get_layer_height(i-1));
                step = step_p + step_m;
                Q_Latent.x[i][j][k] = - lv
                    * (c.x[i+1][j][k] - c.x[i-1][j][k])/step;
                Q_Latent_Ice = - ls 
                    * (ice.x[i+1][j][k] - ice.x[i-1][j][k])/step;
                Q_Sensible.x[i][j][k] = + coeff_sen * coeff_S 
                    * (t.x[i+1][j][k] - t.x[i-1][j][k])/step;  // sensible heat in [W/m2] from energy transport equation
                Q_Latent.x[i][j][k] = 
                    coeff_lat * coeff_L 
                        * (Q_Latent.x[i][j][k] + Q_Latent_Ice);
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Q_Latent.x[im-1][j][k] = Q_Latent.x[im-4][j][k] 
                - 3.0 * Q_Latent.x[im-3][j][k] 
                + 3.0 * Q_Latent.x[im-2][j][k];
            Q_Sensible.x[im-1][j][k] = Q_Sensible.x[im-4][j][k] 
                - 3.0 * Q_Sensible.x[im-3][j][k] 
                + 3.0 * Q_Sensible.x[im-2][j][k];
                }
            }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            Q_Latent.x[i][0][k] = c43 * Q_Latent.x[i][1][k] 
                - c13 * Q_Latent.x[i][2][k];
            Q_Latent.x[i][jm-1][k] = c43 * Q_Latent.x[i][jm-2][k] 
                - c13 * Q_Latent.x[i][jm-3][k];
            Q_Sensible.x[i][0][k] = c43 * Q_Sensible.x[i][1][k] 
                - c13 * Q_Sensible.x[i][2][k];
            Q_Sensible.x[i][jm-1][k] = c43 * Q_Sensible.x[i][jm-2][k] 
                - c13 * Q_Sensible.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Q_Latent.x[i][j][0] = c43 * Q_Latent.x[i][j][1] 
                - c13 * Q_Latent.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            Q_Latent.x[i][j][km-1] = c43 * Q_Latent.x[i][j][km-2] 
                - c13 * Q_Latent.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
/*
            Q_Latent.x[i][j][0] = Q_Latent.x[i][j][1] 
                - 3.0 * Q_Latent.x[i][j][2] 
                + 3.0 * Q_Latent.x[i][j][3];  // extrapolation
            Q_Latent.x[i][j][km-1] = Q_Latent.x[im-4][j][km-4] 
                - 3.0 * Q_Latent.x[im-3][j][km-3] 
                + 3.0 * Q_Latent.x[im-2][j][km-2];  // extrapolation
*/
            Q_Latent.x[i][j][0] = Q_Latent.x[i][j][km-1] 
                = (Q_Latent.x[i][j][0] + Q_Latent.x[i][j][km-1])/2.0;
            Q_Sensible.x[i][j][0] = c43 * Q_Sensible.x[i][j][1] 
                - c13 * Q_Sensible.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            Q_Sensible.x[i][j][km-1] = c43 * Q_Sensible.x[i][j][km-2] 
                - c13 * Q_Sensible.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
/*
            Q_Sensible.x[i][j][0] = Q_Sensible.x[i][j][1] 
                - 3.0 * Q_Sensible.x[i][j][2] 
                + 3.0 * Q_Sensible.x[i][j][3];  // extrapolation
            Q_Sensible.x[i][j][km-1] = Q_Sensible.x[im-4][j][km-4] 
                - 3.0 * Q_Sensible.x[im-3][j][km-3] 
                + 3.0 * Q_Sensible.x[im-2][j][km-2];  // extrapolation
            Q_Sensible.x[i][j][0] = Q_Sensible.x[i][j][km-1] 
                = (Q_Sensible.x[i][j][0] + Q_Sensible.x[i][j][km-1])/2.0;
*/
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i < i_mount)){
                    Q_Latent.x[i][j][k] = Q_Latent.x[i_mount][j][k];
                    Q_Sensible.x[i][j][k] = Q_Sensible.x[i_mount][j][k];
                }
            }
        }
    }
    cout << "      LatentSensibleHeat ended" << endl;
    return;
}
/*
*
*/
//Tao, W.-K., Simpson, J., and McCumber, M.: 
//An Ice-Water Saturation Adjustment, American Meteorological Society, Notes and 
//Correspondence, Volume 1, 321–235, 1988, 1988.
//Ice_Water_SaturationAdjustment, distribution of cloud ice and cloud water 
//dependent on water vapour amount and temperature
void cAtmosphereModel::SaturationAdjustment(){ 
    cout << endl << "      SaturationAdjustment" << endl;
    cout.precision(9);
    bool satadjust = false;
    double saturation = 0.0;
    int i_sat = 0;
    int j_sat = 0;
    int k_sat = 0;
    double height_sat = 0.0;
    double t_latent = 0.0;
    int iter_prec_end = 20;
//    int iter_prec_end = 10;
//    int iter_prec_end = 3;
    int iter_prec = 0;
    double q_v_hyp = 0.0;
    double dt_dim = (dr + dthe + dphi)/3.0 * L_atm/u_0; // dimensional time step = 39.933 s
//    double dt_dim = dt * L_atm/u_0; // dimensional time step = 0.1 s
    double coeff_evap_cond = 0.5; // original COSMO
    double t_u = 0.0;
    double T = 0.0;
    double d_t = 0.0;
    double E_Rain = 0.0;
    double E_Ice = 0.0;
    double q_Rain = 0.0;
    double q_Ice = 0.0;
    double q_v_b = 0.0;
    double q_c_b = 0.0;
    double q_i_b = 0.0;
    double CND = 0.0;
    double DEP = 0.0;
    double d_q_v = 0.0;
    double d_q_c = 0.0;
    double d_q_i = 0.0;
// setting water vapour, cloud water and cloud ice into the proper thermodynamic ratio based on the local temperatures
// starting from a guessed parabolic temperature and water vapour distribution in north/south direction
    double cloud_loc_equator = 22.0;
    double cloud_loc_pole = 18.0;
    for(int k = 1; k < km-1; k++){
        cloud_loc = std::vector<double>(jm, cloud_loc_pole); // radial location of cloud water maximum
        double cloud_loc_eff = cloud_loc_pole - cloud_loc_equator;  // coefficient for the zonal parabolic cloudwater extention
        double d_j_half = (double)(jm-1)/2.0;
        for(int j = 1; j < jm-1; j++){
            double d_j = (double)j;
            cloud_loc[j] = cloud_loc_eff 
                * parabola((double)d_j/(double)d_j_half) + cloud_loc_pole;
            for(int i = 0; i < im; i++){
                double height = get_layer_height(i);
                t_u = t.x[i][j][k] * t_0; // in K
                if(t_u <= t_00){
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k] = 0.0;
                }
                if(t_u > t_0)  ice.x[i][j][k] = 0.0;
                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] < 0.0)  ice.x[i][j][k] = 0.0;
                T = t_u; // in K
                E_Rain = hp * exp_func(T, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Ice = hp * exp_func(T, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_v_b = c.x[i][j][k];
                q_c_b = cloud.x[i][j][k];
                q_i_b = ice.x[i][j][k];
                q_v_hyp = q_v_b;

                if((c.x[i][j][k] >= 0.85 * q_Rain)&&(T > t_00)){ // condition for cloud water and ice formation, available water vapor greater than at saturation
                    // §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§
                    satadjust = true;
//                    saturation = c.x[i][j][k] - 0.85 * q_Rain;
                    for(iter_prec = 1; iter_prec <= iter_prec_end; iter_prec++){ // iter_prec = 2 given by COSMO
                    // condensation ==> water vapor saturation for cloud water formation, deposition ==> ice crystal for cloud ice formation
                        CND = (T - t_00)/(t_0 - t_00); // t_00 = 236.15°C, t_0 = 273.15°C
                        DEP = (t_0 - T)/(t_0 - t_00);
                        if(T <= t_00){
                            CND = 0.0;
                            DEP = 1.0;
                        }
                        if(T >= t_0){
                            CND = 1.0;
                            DEP = 0.0;
                        }
                        d_q_v = q_v_hyp - q_v_b;  // changes in water vapour causing cloud water and cloud ice
                        d_q_c = - d_q_v * CND;
                        d_q_i = - d_q_v * DEP;
                        d_t = (lv * d_q_c + ls * d_q_i)/cp_l; // in K, temperature changes
                        T = T + d_t; // in K
                        t_latent = T - t_u;
                        q_v_b = q_v_b + d_q_v;  // new values
                        q_c_b = q_c_b + d_q_c;
                        q_i_b = q_i_b + d_q_i;
                        if(q_c_b <= 0.0)  q_c_b = 0.0;
                        if(q_i_b <= 0.0)  q_i_b = 0.0;
                        E_Rain = hp * exp_func(T, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        E_Ice = hp * exp_func(T, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
                        q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                        q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                        if((q_c_b > 0.0) && (q_i_b > 0.0))
                            q_v_hyp = (q_c_b * q_Rain + q_i_b * q_Ice) 
                                /(q_c_b + q_i_b);
                        if((q_c_b >= 0.0) && (q_i_b == 0.0))  q_v_hyp = q_Rain;
                        if((q_c_b == 0.0) && (q_i_b > 0.0))  q_v_hyp = q_Ice;
                        if(T >= t_0) q_i_b = 0.0;
                        if(T <= t_00){
                            q_v_b = 0.0;
                            q_c_b = 0.0;
                            q_i_b = 0.0;
                        }
                        saturation = q_v_b - 0.85 * q_Rain;
                        t_latent = T - t_u;
                        S_c_c.x[i][j][k] = coeff_evap_cond  // negative values == condensation, positive values == evaporation
                            * d_q_c/dt_dim;
/*
                        cout.precision(10);
                        cout.setf(ios::fixed);
                        if((j == 90)&&(k == 180))  cout << endl
                            << "  SaturationAdjustment" << endl 
                            << "  Ma = " << (int)*get_current_time() << endl
                            << "  height = "<< height << endl
                            << "  saturation = "<< saturation * 1e3 << endl
                            << "  i = " << i << "  j = " << j << "  k = " << k << endl
                            << "  iprec = "<< iter_prec << endl
                            << "  CND = " << CND
                            << "  DEP = " << DEP << endl
                            << "  p_hydro = " << p_hydro.x[i][j][k] 
                            << "  p_stat = " << p_stat.x[i][j][k] << endl
                            << "  r_dry = " << r_dry.x[i][j][k] 
                            << "  r_humid = " << r_humid.x[i][j][k] << endl
                            << "  dt = " << d_t 
                            << "  t_latent = " << t_latent
                            << "  T = " << T - t_0 
                            << "  t = " << t_u - t_0 << endl
                            << "  q_Rain = " << q_Rain * 1e3 
                            << "  q_Ice = " << q_Ice * 1e3 << endl
                            << "  q_dif = " << fabs(q_v_b/q_v_hyp - 1.) 
                            << "  d_q_v = " << d_q_v * 1e3 
                            << "  d_q_c = " << d_q_c * 1e3 
                            << "  d_q_i = " << d_q_i * 1e3 << endl
                            << "  q_v_hyp = " << q_v_hyp * 1e3
                            << "  q_v_b = " << q_v_b * 1e3 
                            << "  q_c_b = " << q_c_b * 1e3 
                            << "  q_i_b = " << q_i_b * 1e3 << endl
                            << "  c = " << c.x[i][j][k] * 1e3 
                            << "  cloud = " << q_v_b * 1e3 
                            << "  ice = " << q_c_b * 1e3 << endl
                            << "  Scc = " << S_c_c.x[i][j][k] * 1e3 << endl;
*/
                        if(fabs(q_v_b/q_v_hyp - 1.0) <= 1.e-4){
                            saturation = q_v_b - 0.85 * q_Rain;
                            t_latent = T - t_u;
                            i_sat = i;
                            j_sat = j;
                            k_sat = k;
                            break;
                        }
                        else q_v_hyp = 0.5 * (q_v_hyp + q_v_b);  // has smoothing effect
                    } // iter_prec end

                    c.x[i][j][k] = q_v_b;  // new values achieved after converged iterations
                    cloud.x[i][j][k] = q_c_b;
                    ice.x[i][j][k] = q_i_b;
                    t.x[i][j][k] = t.x[i][j][k] - (T - t_u)/t_0;
                } // iterations for mixed cloud phase
            } // end i
        } // end j
    } // end k
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
                if(is_land(h, 0, j, k)){
                    c.x[0][j][k] = c.x[i_mount][j][k];
                    cloud.x[0][j][k] = cloud.x[i_mount][j][k];
                    ice.x[0][j][k] = ice.x[i_mount][j][k];
                }
            for(int i = i_mount; i >= 0; i--){

                if(is_land(h, i, j, k)){
                    c.x[i][j][k] = c.x[i_mount][j][k];
                    cloud.x[i][j][k] = cloud.x[i_mount][j][k];
                    ice.x[i][j][k] = ice.x[i_mount][j][k];
                }

                if(t.x[i][j][k] * t_0 <= t_00){
                    c.x[i][j][k] = 0.0;
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k] = 0.0;
                }
            }
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(cloud.x[i][j][k] <= 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] <= 0.0) ice.x[i][j][k] = 0.0;
            }
        }
    }
    if(satadjust == false)  
        cout << "      no saturation of water vapour in SaturationAdjustment found" 
        << endl;
    else
        cout << "      saturation of water vapour in SaturationAdjustment found" 
        << endl
        << "      iter_prec = " << iter_prec << endl
        << "      i_sat = " << i_sat
        << "   j_sat = " << j_sat
        << "   k_sat = " << k_sat
        << "   height_sat[m] = " << height_sat
        << "   t_latent = " << t_latent
        << "   saturation[g/kg] = " << saturation * 1e3 << endl;
    cout << "      SaturationAdjustment ended" << endl;
    return;
}
/*
*
*/
// One-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
// resulting the precipitation distribution formed of rain and snow
void cAtmosphereModel::ZeroCategoryIceScheme(){   
    cout << endl << "      ZeroCategoryIceScheme" << endl;
    // constant coefficients for the transport of cloud water and cloud ice amount, 
    // rain and snow in the parameterization procedures
    bool rain = false;
    double P_Rain = 0.0;
    double Rain = 0.0;
    int i_rain = 0;
    int j_rain = 0;
    int k_rain = 0;
    double height_rain = 0.0;
    double maxValue_rain = 0.0;
    double N_r_0 = 8.0e+6,  // 1/m4
           v_r_0 = 4.9,  // m/s
           c_r_t = 12.63,  // 
           tau_r = 3.3e3,  // s          can be adjusted to fit best the average NASA-precipitation of 2.68 mm/d
           c_ac = 0.24,  // m2/kg
           b_ev = 8.05;  // m2*s/kg
    double N_cf_0 = 0.0;
    double a_ev = 0.0;
    double N_cf = 0.0;
    double E_Rain = 0.0;
    double q_Rain = 0.0;
    double A_r = r_0_water * M_PI * N_r_0,
           B_r = r_0_water * M_PI * N_r_0 * v_r_0 
                * tgamma(4.5)/tgamma(4.0);
    double S_au, S_ac, S_ev;
    double r_q_r, v_r_t;
    std::vector<double> step(im, 0.0);
    // rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
    // the choosen scheme is a Two Category Ice Scheme
    // besides the transport equation for the water vapour exists two equations for the cloud water and the cloud ice transport
    // since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                P_rain.x[i][j][k] = 0.0;
            }
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            P_rain.x[im-1][j][k] = 0.0;
            P_rain.x[im-2][j][k] = 0.0;
            S_r.x[im-1][j][k] = 0.0;
            S_r.x[im-2][j][k] = 0.0;
            for(int i = im-2; i >= 0; i--){
                Rain = P_rain.x[i][j][k];
                r_q_r = A_r * pow(B_r,-(8.0/9.0)) 
                    * pow(P_rain.x[i+1][j][k],(8.0/9.0));
                v_r_t = c_r_t * pow(r_q_r, (1.0/8.0));
                double t_u = t.x[i][j][k] * t_0;
                step[i] = get_layer_height(i) - get_layer_height(i-1);  // local atmospheric shell thickness
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
    // collection mechanisms accretion
    // accretion of cloud water by raindrops
                S_ac = c_ac * cloud.x[i][j][k]  // c_ac = 0.24, in m2/kg
                    * pow(Rain,(7.0/9.0));
    // autoconversion processes
    // autoconversion of cloud water to form rain
                S_au = max(cloud.x[i][j][k], 0.0)/tau_r;
    // diffusional growth of rain
    // evaporation of rain water
                if(t_u >= t_0){  // temperature below zero
                    a_ev = 2.76e-3 * exp(0.055 * (t_0 - t_u));
                    S_ev = a_ev * (1.0 + b_ev * pow(Rain,(1.0/6.0)))  // evaporation of rain due to water vapour diffusion
                        * (q_Rain - c.x[i][j][k]) 
                        * pow(Rain, (4.0/9.0));
                }else  S_ev = 0.0;
    // sinks and sources
                S_v.x[i][j][k] = - S_c_c.x[i][j][k] + S_ev;
                S_c.x[i][j][k] = S_c_c.x[i][j][k] - S_au - S_ac; 
                S_r.x[i][j][k] = S_au + S_ac - S_ev;

                if(is_land(h,i,j,k)){
                    S_c_c.x[i][j][k] = 0.0;
                    S_v.x[i][j][k] = 0.0;
                    S_c.x[i][j][k] = 0.0;
                    S_r.x[i][j][k] = 0.0;
                }
    // rain and snow integration
                if(t_u >= t_0)
/*
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * 0.5 * (S_r.x[i+1][j][k] 
                        + S_r.x[i][j][k]) * step[i];  // in kg/(m2 * s) == mm/s 
*/
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_r.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] >= 20.0/8.64e4)  P_rain.x[i][j][k] = 20.0/8.64e4;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] > 0.0){
                    rain = true;
                    if(P_rain.x[i][j][k] > maxValue_rain){
                        maxValue_rain = P_rain.x[i][j][k];
                        P_Rain = maxValue_rain * 8.64e4; // in mm/d
                        i_rain = i;
                        j_rain = j;
                        k_rain = k;
                        height_rain = get_layer_height(i);
/*
                        cout.precision(9);
                        cout << "      max values of P_rain " 
                            << endl
                            << "      i__rain = " << i_rain
                            << "      j_rain = " << j_rain
                            << "   k_rain = " << k_rain
                            << "   height_rain = " << height_rain
                            << "   P_Rain = " << P_Rain << endl;
*/
                    }
                }
/*
                cout.precision(10);
                cout.setf(ios::fixed);
                if((j == 124) && (k == 87)) cout << endl 
                    << "  ZeroCategoryIceScheme ----- SINKS and SOURCES ----- " << endl 
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  height = " << get_layer_height(i)
                    << "  step = " << step[i] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl << endl

                    << "  b)   autoconversion of rain" << endl
                    << "  Sau = " << S_au * 1e6 << endl << endl

                    << "  b)   collection mechanism" << endl
                    << "  Sac = " << S_ac * 1e6 << endl << endl

                    << "  d)   evaporation of rain" << endl
                    << "  Sev = " << S_ev * 1e6 << endl << endl;
*/
/*
                cout.precision(10);
                cout.setf(ios::fixed);
                if((j == 124) && (k == 87)) cout << endl 
                    << "  ZeroCategoryIceScheme ----- RAIN and SNOW PRODUCTION -----" << endl
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  step = " << step[i] << endl
                    << "  N_cf = " << N_cf
                    << "  N_cf_0 = " << N_cf_0 << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl
                    << "  c = " << c.x[i][j][k] * 1e3
                    << "  cl = " << cloud.x[i][j][k] * 1e3
                    << "  ice = " << ice.x[i][j][k] * 1e3 << endl
                    << "  r_dry = " << r_dry.x[i][j][k] 
                    << "  r_humid = " << r_humid.x[i][j][k] << endl
                    << "  p_hydro = " << p_hydro.x[i][j][k]
                    << "  p_stat = " << p_stat.x[i][j][k] << endl << endl

                    << "  a)   source and sink term due to condensation and evaporation" << endl
                    << "  S_c_c = " << S_c_c.x[i][j][k] * 1e3 << endl << endl

                    << "  b)   source and sink terms for the water categories: water vapour, cloud water and cloud ice" << endl
                    << "  S_v = - S_c_c + S_ev" << endl
                    << "  S_c = S_c_c - S_au - S_ac" << endl
                    << "  S_v = " << S_v.x[i][j][k] * 1e6 
                    << "  S_c = " << S_c.x[i][j][k] * 1e6  << endl << endl

                    << "  c)   source and sink terms for the water categories: rain" << endl
                    << "  S_r = S_au + S_ac - S_ev" << endl
                    << "  S_s = S_nuc + S_rim + S_dep + S_frz - S_melt" << endl
                    << "  S_r = " << S_r.x[i][j][k] * 1e6 
                    << "  S_s = " << S_s.x[i][j][k] * 1e6 << endl << endl

                    << "  d)   rain and snow formed by the sources and sinks in mm/d" << endl
                    << "  r_q_r = " << r_q_r * 8.64e4
                    << "  v_r_t = " << v_r_t << endl
                    << "  A_r = " << A_r
                    << "  B_r = " << B_r << endl << endl

                    << "  P_rain = " << P_rain.x[i][j][k] * 8.64e4
                    << "  P_conv_midl = " << P_conv_midl.x[i][j][k] * 8.64e4
                    << "  P_tot = " << P_rain.x[i][j][k] * 8.64e4 + P_snow.x[i][j][k] * 8.64e4 + P_conv_midl.x[i][j][k] * 8.64e4
                    << "  P_NASA = " << precipitation_NASA.y[j][k]  << endl << endl << endl;
*/
            }  // end i RainSnow
        }  // end j
    }  // end k

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            P_rain.x[im-1][j][k] = c43 * P_rain.x[im-2][j][k] 
                - c13 * P_rain.x[im-3][j][k];
            P_snow.x[im-1][j][k] = c43 * P_snow.x[im-2][j][k] 
                - c13 * P_snow.x[im-3][j][k];
                }
            }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            P_rain.x[i][0][k] = c43 * P_rain.x[i][1][k] 
                - c13 * P_rain.x[i][2][k];
            P_rain.x[i][jm-1][k] = c43 * P_rain.x[i][jm-2][k] 
                - c13 * P_rain.x[i][jm-3][k];
            P_snow.x[i][0][k] = c43 * P_snow.x[i][1][k] 
                - c13 * P_snow.x[i][2][k];
            P_snow.x[i][jm-1][k] = c43 * P_snow.x[i][jm-2][k] 
                - c13 * P_snow.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] 
                - c13 * P_rain.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] 
                - c13 * P_rain.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] 
                = (P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/2.0;
            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] 
                - c13 * P_snow.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] 
                - c13 * P_snow.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
        }
    }
/*
*
*/
    AtomUtils::fft_gaussian_filter_3d(P_rain,1);
    AtomUtils::fft_gaussian_filter_3d(P_snow,1);
/*
*
*/
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i < i_mount)){
                    P_rain.x[i][j][k] = P_rain.x[i_mount][j][k];
                    P_snow.x[i][j][k] = P_snow.x[i_mount][j][k];
                }
            }
        }
    }
    if(rain == false)  
        cout << "      no rain fall in ZeroCategoryIceScheme found" 
        << endl;
    else
        cout << "      rain fall in ZeroCategoryIceScheme found" 
        << endl
        << "      i_rain = " << i_rain
        << "      j_rain = " << j_rain
        << "   k_rain = " << k_rain
        << "   height_rain[m] = " << height_rain
        << "   P_Rain[mm/d] = " << P_Rain << endl;
    cout << "      ZeroCategoryIceScheme ended" << endl;
    return;
}
/*
*
*/
// One-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
// resulting the precipitation distribution formed of rain and snow
void cAtmosphereModel::OneCategoryIceScheme(){   
    cout << endl << "      OneCategoryIceScheme" << endl;
    // constant coefficients for the transport of cloud water and cloud ice amount, 
    // rain and snow in the parameterization procedures
    bool rain = false;
    bool snow = false;
    double P_Rain = 0.0;
    double P_Snow = 0.0;
    double Rain = 0.0;
    double Snow = 0.0;
    int i_rain = 0;
    int i_snow = 0;
    int j_rain = 0;
    int j_snow = 0;
    int k_rain = 0;
    int k_snow = 0;
    double height_rain = 0.0;
    double height_snow = 0.0;
    double maxValue_rain = 0.0;
    double maxValue_snow = 0.0;
    double N_r_0 = 8.0e+6,  // in 1/m-4
           N_s_0 = 4.0e+5,  // in 1/m-4
           v_s_0 = 130.0,  // in m/s
           v_r_0 = 4.9,  // in m/s
           c_r_t = 12.63,  // in 
           c_s_t = 2.87,  // in 
           a_s_m = 0.038,  // in kg/m2
           c_ac = 0.24,  // m2/kg
           c_rim = 0.69,  // m2/kg
           b_ev = 5.98,  // m2*s/kg
           a_melt = 3.90e-6, // K/(kg/kg)
           b_melt = 10.50,  // m2*s/kg
           b_dep = 10.50,  // m2*s/kg
           a_if = 1.92e-6,
           a_cf = 3.97e-5,
           E_cf = 5.0e-3,
           N_cf_0_surf = 2.0e5,  // 1/m3
           N_cf_0_top = 1.0e4,  // 1/m3
//           tau_r = 1.0e4,  // s
           tau_r = 3.3e3,  // s          can be adjusted to fit best the average NASA-precipitation of 2.68 mm/d
           tau_s = 1.0e3,  // s
           a_mc = 0.08,  // kg/m2
           a_mv = 0.02;  // kg/m2
    double N_cf_0 = 0.0;
    double eps_t = 0.0;
    double a_m = 0.0;
    double a_ev = 0.0;
    double a_dep = 0.0;
    double N_cf = 0.0;
    double t_m1 = 0.5 * (t_0 + t_000);
    double t_m2 = 0.5 * (t_0 + t_00);
    double q_Rain = 0.0,
           E_Rain = 0.0;
    double q_Ice = 0.0,
           E_Ice = 0.0;
    double A_r = r_0_water * M_PI * N_r_0,
           A_s = 2.0 * a_s_m * N_s_0, 
           B_r = r_0_water * M_PI * N_r_0 * v_r_0 * tgamma(4.5)/tgamma(4.0), 
           B_s = N_s_0 * v_s_0 * a_s_m * tgamma(3.25);
    double S_nuc, S_frz, S_cf_frz, S_if_frz, S_dep, S_au, S_ac, S_rim, S_shed, S_ev, S_melt;
    double r_q_r, r_q_s, v_r_t, v_s_t;
    std::vector<double> step(im, 0.0);
    // rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
    // the choosen scheme is a Two Category Ice Scheme
    // besides the transport equation for the water vapour exists two equations for the cloud water and the cloud ice transport
    // since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] < 0.0)  ice.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                P_rain.x[i][j][k] = 0.0;
                P_snow.x[i][j][k] = 0.0;
            }
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            P_rain.x[im-1][j][k] = 0.0;
            P_snow.x[im-1][j][k] = 0.0;
            P_rain.x[im-2][j][k] = 0.0;
            P_snow.x[im-2][j][k] = 0.0;
            S_r.x[im-1][j][k] = 0.0;
            S_s.x[im-1][j][k] = 0.0;
            S_r.x[im-2][j][k] = 0.0;
            S_s.x[im-2][j][k] = 0.0;
            for(int i = im-2; i >= 0; i--){
                Rain = P_rain.x[i][j][k];
                Snow = P_snow.x[i][j][k];
                r_q_r = A_r * pow(B_r,-(8.0/9.0)) 
                    * pow(P_rain.x[i][j][k],(8.0/9.0));
                r_q_s = A_s * pow(B_s,-(12.0/13.0)) 
                    * pow(P_snow.x[i][j][k],(12.0/13.0));
                v_r_t = c_r_t * pow(r_q_r, (1.0/8.0));
                v_s_t = c_s_t * pow(r_q_s, (1.0/12.0));
                double t_u = t.x[i][j][k] * t_0;
                step[i] = get_layer_height(i) - get_layer_height(i-1);  // local atmospheric shell thickness
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
    // mass size relation of circular plates
                if((t_u < t_0)&&(t_u >= t_m1))
                    a_m = a_mc - a_mv * (1.0 + cos(2.0 * M_PI 
                        * (t_u - t_m1)/(t_0 - t_000)));
                else  a_m = a_mc;
    // autoconversion and nucleation process, epsilon(T)
                if((t_u > t_00)&&(t_u <= t_0))
                    eps_t = 0.5 * (1.0 + sin(M_PI*(t_m2 - t_u)/(t_0 - t_00)));
                else
                    eps_t = 0.0;
//                if (cloud.x[i][j][k] >= 0.0){
//                    S_au = (1.0 - eps_t)/tau_r * cloud.x[i][j][k];
//                    S_nuc = eps_t/tau_s * cloud.x[i][j][k];
                    S_au = (1.0 - eps_t)/tau_r * max(0.0, cloud.x[i][j][k]);
                    S_nuc = eps_t/tau_s * max(0.0, cloud.x[i][j][k]);
//                }
    // collection mechanisms accretion, riming and shedding
    // accretion of cloud water by raindrops
                S_ac = (1.0 - eps_t) * c_ac * cloud.x[i][j][k]  // c_ac = 0.24, in m2/kg
                    * pow(Rain,(7.0/9.0));
                if(t_u < t_0)
                    S_rim = c_rim/a_m * cloud.x[i][j][k] * Snow;  // c_rim = 18.6,  m2/kg
                else  S_rim = 0.0;  // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >
                                   // by falling snow particles
                if(t_u >= t_0)  
                    S_shed = c_rim/a_m * cloud.x[i][j][k] * Snow;
                else  S_shed = 0.0;  // rate of water shed by melting wet snow particles, < IX >
                                    // collecting cloud droplets to produce rain
    // diffusional growth of rain and snow
    // evaporation of rain water
                a_ev = 2.76e-3 * exp(0.055 * (t_0 - t_u));
                S_ev = a_ev * (1.0 + b_ev * pow(Rain,(1.0/6.0)))  // evaporation of rain due to water vapour diffusion, < XIII >
                    * (q_Rain - c.x[i][j][k]) 
                    * pow(Rain, (4.0/9.0));
    // deposition growth and sublimation of cloud ice
                a_dep = 1.13e-3 * exp(0.073 * (t_0 - t_u));
                S_dep = a_dep/pow(a_m, - 0.5) * (1.0 + b_dep   // melting rate of snow to form rain, < XVI >
                    * pow(a_m, - 0.25) * pow(Snow, (0.9/4.3)))  // c_s_melt = 8.43e-5, (m2*s)/(K*kg)
                    * ((c.x[i][j][k] - q_Ice) * pow(Snow, (5.0/8.6)));
    // melting of snow to form cloud water
                S_melt = a_melt/pow(a_mc, - 0.5) * (1.0 + b_melt   // melting rate of snow to form rain, < XVI >
                    * pow(a_mc, - 0.25) * pow(Snow, (0.9/4.3)))  // c_s_melt = 8.43e-5, (m2*s)/(K*kg)
                    * ((t_u - t_0) * pow(Snow, (5.0/8.6)));
    // freezing of rain to form snow
                N_cf_0 = N_cf_0_surf + (N_cf_0_surf - N_cf_0_top)
                    /(p_hydro.x[0][j][k] - 500.0) 
                    * (p_hydro.x[i][j][k] - p_hydro.x[0][j][k]);
                if(p_hydro.x[i][j][k] <= 500.0)  N_cf_0 = N_cf_0_top;
                if(t_u < 270.16){
                    N_cf = N_cf_0 * pow((270.16 - t_u), 1.3);
                }
                else  N_cf = 0.0;
                if((t_u >= t_0)&&(t_u <= t_00)){
                    S_if_frz = a_if * (exp(a_if * (t_u - t_0)) - 1.0)  // immersion freezing
                        * pow(Rain, (14.0/9.0));
                    S_cf_frz = a_cf * E_cf * N_cf                       // contact freezing nucleation
                        * pow(Rain, (13.0/9.0));
                    S_frz = S_if_frz + S_cf_frz;                        // mixing ratio of snow
                }
                else  S_frz = 0.0;
    // sinks and sources
                S_v.x[i][j][k] = - S_c_c.x[i][j][k] + S_ev - S_dep;
                S_c.x[i][j][k] = S_c_c.x[i][j][k] - S_au - S_ac 
                    - S_nuc - S_rim - S_shed;
                S_r.x[i][j][k] = S_au + S_ac - S_ev + S_shed 
                    - S_frz + S_melt;
                S_s.x[i][j][k] = S_nuc + S_rim + S_dep + S_frz - S_melt;

                if(is_land(h,i,j,k)){
                    S_c_c.x[i][j][k] = 0.0;
                    S_v.x[i][j][k] = 0.0;
                    S_c.x[i][j][k] = 0.0;
                    S_r.x[i][j][k] = 0.0;
                    S_s.x[i][j][k] = 0.0;
                }
    // rain and snow integration
                if(t_u >= t_0)
/*
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * 0.5 * (S_r.x[i+1][j][k] 
                        + S_r.x[i][j][k]) * step[i];  // in kg/(m2 * s) == mm/s 
*/
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_r.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] >= 20.0/8.64e4)  P_rain.x[i][j][k] = 20.0/8.64e4;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] > 0.0){
                    rain = true;
                    if(P_rain.x[i][j][k] > maxValue_rain){
                        maxValue_rain = P_rain.x[i][j][k];
                        P_Rain = maxValue_rain * 8.64e4; // in mm/d
                        i_rain = i;
                        j_rain = j;
                        k_rain = k;
                        height_rain = get_layer_height(i);
/*
                        cout.precision(9);
                        cout << "      max values of P_rain " 
                            << endl
                            << "      i__rain = " << i_rain
                            << "      j_rain = " << j_rain
                            << "   k_rain = " << k_rain
                            << "   height_rain = " << height_rain
                            << "   P_Rain = " << P_Rain << endl;
*/
                    }
                }
                if((t_u < t_0)&&(t_u >= t_00))
/*
                    P_snow.x[i][j][k] = P_snow.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * 0.5 * (S_s.x[i+1][j][k] 
                        + S_s.x[i][j][k]) * step[i];  // in kg/(m2 * s) == mm/s 
*/
                    P_snow.x[i][j][k] = P_snow.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_s.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_snow.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] > 0.0){
                    snow = true;
                    if(P_snow.x[i][j][k] > maxValue_snow){
                        maxValue_snow = P_snow.x[i][j][k];
                        P_Snow = maxValue_snow * 8.64e4; // in mm/d
                        i_snow = i;
                        j_snow = j;
                        k_snow = k;
                        height_snow = get_layer_height(i);
/*
                        cout.precision(9);
                        cout << "      max values of P_snow " 
                            << endl
                            << "      i__snow = " << i_snow
                            << "      j_snow = " << j_snow
                            << "   k_snow = " << k_snow
                            << "   height_snow = " << height_snow
                            << "   P_Snow = " << P_Snow << endl;
*/
                    }
                }
/*
                cout.precision(10);
                cout.setf(ios::fixed);
                if((j == 124) && (k == 87)) cout << endl 
                    << "  OneCategoryIceScheme ----- SINKS and SOURCES ----- " << endl 
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  height = " << get_layer_height(i)
                    << "  step = " << step[i] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl << endl

                    << "  a)   nucleation and depositional growth of snow" << endl
                    << "  Snuc = " << S_nuc * 1e3 
                    << "  Sdep = " << S_dep * 1e3 << endl << endl

                    << "  b)   autoconversion of rain" << endl
                    << "  Sau = " << S_au * 1e3 << endl << endl

                    << "  c)   collection mechanism" << endl
                    << "  Sac = " << S_ac * 1e3
                    << "  Srim = " << S_rim * 1e3 
                    << "  Sshed = " << S_shed * 1e3 << endl << endl

                    << "  d)   evaporation of rain" << endl
                    << "  Sev = " << S_ev * 1e3 << endl << endl

                    << "  e)   melting and freezing of snow" << endl
                    << "  Smelt = " << S_melt * 1e3 
                    << "  Sfrz = " << S_frz * 1e3 << endl << endl;
*/
/*
                cout.precision(10);
                cout.setf(ios::fixed);
                if((j == 124) && (k == 87)) cout << endl 
                    << "  OneCategoryIceScheme ----- RAIN and SNOW PRODUCTION -----" << endl
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  step = " << step[i] << endl
                    << "  N_cf = " << N_cf
                    << "  N_cf_0 = " << N_cf_0 << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl
                    << "  c = " << c.x[i][j][k] * 1.0e6
                    << "  cl = " << cloud.x[i][j][k] * 1.0e6
                    << "  ice = " << ice.x[i][j][k] * 1.0e6 << endl
                    << "  r_dry = " << r_dry.x[i][j][k] 
                    << "  r_humid = " << r_humid.x[i][j][k] << endl
                    << "  p_hydro = " << p_hydro.x[i][j][k]
                    << "  p_stat = " << p_stat.x[i][j][k] << endl << endl

                    << "  a)   source and sink term due to condensation and evaporation" << endl
                    << "  S_c_c = " << S_c_c.x[i][j][k] * 1.0e6 << endl << endl

                    << "  b)   source and sink terms for the water categories: water vapour, cloud water and cloud ice" << endl
                    << "  S_v = - S_c_c + S_ev - S_dep" << endl
                    << "  S_c = S_c_c - S_au - S_ac - S_nuc - S_rim - S_shed" << endl
                    << "  S_v = " << S_v.x[i][j][k] * 1.0e6 
                    << "  S_c = " << S_c.x[i][j][k] * 1.0e6  << endl << endl

                    << "  c)   source and sink terms for the water categories: rain and snow" << endl
                    << "  S_r = S_au + S_ac - S_ev + S_shed - S_frz + S_melt" << endl
                    << "  S_s = S_nuc + S_rim + S_dep + S_frz - S_melt" << endl
                    << "  S_r = " << S_r.x[i][j][k] * 1.0e6 
                    << "  S_s = " << S_s.x[i][j][k] * 1.0e6 << endl << endl

                    << "  d)   rain and snow formed by the sources and sinks in mm/d" << endl
                    << "  r_q_r = " << r_q_r * 8.64e4
                    << "  r_q_s = " << r_q_s * 8.64e4 << endl
                    << "  v_r_t = " << v_r_t
                    << "  v_s_t = " << v_s_t << endl
                    << "  A_r = " << A_r
                    << "  A_s = " << A_s
                    << "  B_r = " << B_r
                    << "  B_s = " << B_s << endl << endl

                    << "  P_rain = " << P_rain.x[i][j][k] * 8.64e4
                    << "  P_snow = " << P_snow.x[i][j][k] * 8.64e4
                    << "  P_conv_midl = " << P_conv_midl.x[i][j][k] * 8.64e4
                    << "  P_tot = " << P_rain.x[i][j][k] * 8.64e4 + P_snow.x[i][j][k] * 8.64e4 + P_conv_midl.x[i][j][k] * 8.64e4
                    << "  P_NASA = " << precipitation_NASA.y[j][k]  << endl << endl << endl;
*/
            }  // end i RainSnow
        }  // end j
    }  // end k

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            P_rain.x[im-1][j][k] = c43 * P_rain.x[im-2][j][k] 
                - c13 * P_rain.x[im-3][j][k];
            P_snow.x[im-1][j][k] = c43 * P_snow.x[im-2][j][k] 
                - c13 * P_snow.x[im-3][j][k];
                }
            }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            P_rain.x[i][0][k] = c43 * P_rain.x[i][1][k] 
                - c13 * P_rain.x[i][2][k];
            P_rain.x[i][jm-1][k] = c43 * P_rain.x[i][jm-2][k] 
                - c13 * P_rain.x[i][jm-3][k];
            P_snow.x[i][0][k] = c43 * P_snow.x[i][1][k] 
                - c13 * P_snow.x[i][2][k];
            P_snow.x[i][jm-1][k] = c43 * P_snow.x[i][jm-2][k] 
                - c13 * P_snow.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] 
                - c13 * P_rain.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] 
                - c13 * P_rain.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] 
                = (P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/2.0;
            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] 
                - c13 * P_snow.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] 
                - c13 * P_snow.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
        }
    }
/*
*
*/
    AtomUtils::fft_gaussian_filter_3d(P_rain,1);
    AtomUtils::fft_gaussian_filter_3d(P_snow,1);
/*
*
*/
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i < i_mount)){
                    P_rain.x[i][j][k] = P_rain.x[i_mount][j][k];
                    P_snow.x[i][j][k] = P_snow.x[i_mount][j][k];
                }
            }
        }
    }
    if(rain == false)  
        cout << "      no rain fall in OneCategoryIceScheme found" 
        << endl;
    else
        cout << "      rain fall in OneCategoryIceScheme found" 
        << endl
        << "      i_rain = " << i_rain
        << "      j_rain = " << j_rain
        << "   k_rain = " << k_rain
        << "   height_rain[m] = " << height_rain
        << "   P_Rain[mm/d] = " << P_Rain << endl;
    if(snow == false)  
        cout << "      no snow fall in OneCategoryIceScheme found" 
        << endl;
    else
        cout << "      snow fall in OneCategoryIceScheme found" 
        << endl
        << "      i_snow = " << i_snow
        << "      j_snow = " << j_snow
        << "   k_snow = " << k_snow
        << "   height_snow[m] = " << height_snow
        << "   P_Snow[mm/d] = " << P_Snow << endl;
    cout << "      OneCategoryIceScheme ended" << endl;
    return;
}
// Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
// resulting the precipitation distribution formed of rain and snow
void cAtmosphereModel::TwoCategoryIceScheme(){   
    cout << endl << "      TwoCategoryIceScheme" << endl;
    // constant coefficients for the transport of cloud water and cloud ice amount, 
    // rain and snow in the parameterization procedures
    bool rain = false;
    bool snow = false;
    double P_Rain = 0.0;
    double P_Snow = 0.0;
    double Rain = 0.0;
    double Snow = 0.0;
    int i_rain = 0;
    int i_snow = 0;
    int j_rain = 0;
    int j_snow = 0;
    int k_rain = 0;
    int k_snow = 0;
    double height_rain = 0.0;
    double height_snow = 0.0;
    double maxValue_rain = 0.0;
    double maxValue_snow = 0.0;
    double N_i_0 = 1.0e2,  // in m-3
           m_i_0 = 1.0e-12,  // in kg
           m_i_max = 1.0e-9,  // in kg
           m_s_0 = 3.0e-9,  // in kg
           N_r_0 = 8.0e+6,  // in 1/m-4
           N_s_0 = 8.0e+5,  // in 1/m-4
           v_s_0 = 9.356,  // in m/s
           v_r_0 = 130.0,  // in m/s
           c_r_t = 12.63,  // in 
           c_s_t = 2.49,  // in 
           a_s_m = 0.038,  // in kg/m2
           c_i_dep = 1.3e-5,  // in m3/s
           c_c_au = 4.0e-4,  // in 1/s COSMO
           c_i_au = 1.0e-3,  // in 1/s COSMO
           c_ac = 0.24,  // m2/kg
           c_rim = 18.6,  // m2/kg
           c_agg = 10.3,  // m2/kg
           c_i_cri = 0.24,  // m2
           c_r_cri = 3.2e-5,  // m2
           a_ev = 1.0e-3,  // m2/kg
           b_ev = 5.9,  // m2*s/kg
           c_s_dep = 1.8e-2,  // m2/kg
           b_s_dep = 12.3,  // m2*s/kg
           c_s_melt = 8.43e-5,  // (m2*s)/(K*kg)
           b_s_melt = 12.05,  // m2*s/kg
           a_s_melt = 2.31e3, // K/(kg/kg)
           c_r_frz = 3.75e-2,  // m2/(K*kg)
           t_nuc = 267.15,  // in K    -6 °C
           t_d = 248.15,  // in K    -25 °C
           t_hn = 236.15,  // in K    -37 °C
           t_r_frz = 271.15;  // in K    -2 °C
    double q_Rain = 0.0,
           E_Rain = 0.0;
    double q_Ice = 0.0,
           E_Ice = 0.0;
    double dt_snow_dim = 0.0,  // dt_snow_dim is the fallout velocity
           dt_rain_dim = 0.0;  // dt_rain_dim is the fallout velocity
    double m_i = m_i_max;  
    double A_r = r_0_water * M_PI * N_r_0,
           A_s = 2.0 * a_s_m * N_s_0, 
           B_r = r_0_water * M_PI * N_r_0 * v_r_0 * tgamma(4.5)/tgamma(4.0), 
           B_s = N_s_0 * v_s_0 * a_s_m * tgamma(3.25);
    double N_i, S_nuc, S_c_frz, S_i_dep, S_c_au, S_i_au, S_d_au, 
           S_ac, S_rim, S_shed;
    double S_agg, S_i_cri, S_r_cri, S_ev, S_s_dep, S_i_melt, S_s_melt, S_r_frz;
    double r_q_r, r_q_s, v_r_t, v_s_t;
    std::vector<double> step(im, 0.0);
    // rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
    // the choosen scheme is a Two Category Ice Scheme
    // besides the transport equation for the water vapour exists two equations for the cloud water and the cloud ice transport
    // since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] < 0.0)  ice.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                P_rain.x[i][j][k] = 0.0;
                P_snow.x[i][j][k] = 0.0;
            }
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            P_rain.x[im-1][j][k] = 0.0;
            P_snow.x[im-1][j][k] = 0.0;
            P_rain.x[im-2][j][k] = 0.0;
            P_snow.x[im-2][j][k] = 0.0;
            S_r.x[im-1][j][k] = 0.0;
            S_s.x[im-1][j][k] = 0.0;
            S_r.x[im-2][j][k] = 0.0;
            S_s.x[im-2][j][k] = 0.0;
            for(int i = im-2; i >= 0; i--){
                Rain = P_rain.x[i][j][k];
                Snow = P_snow.x[i][j][k];
                r_q_r = A_r * pow(B_r,-(8.0/9.0)) 
                    * pow(P_rain.x[i][j][k],(8.0/9.0));
                r_q_s = A_s * pow(B_s,-(12.0/13.0)) 
                    * pow(P_snow.x[i][j][k],(12.0/13.0));
                v_r_t = c_r_t * pow(r_q_r, (1.0/8.0));
                v_s_t = c_s_t * pow(r_q_s, (1.0/12.0));
                double t_u = t.x[i][j][k] * t_0;
                step[i] = get_layer_height(i) - get_layer_height(i-1);  // local atmospheric shell thickness
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                dt_rain_dim = step[i]/1.6; // adjusted rain fall time step by fixed velocities == 1.6 m/s by variable local step size
                dt_snow_dim = step[i]/0.96; // adjusted snow fall time step by fixed velocities == 0.96 m/s by variable local step size
        // number density of ice particles and mass of cloud ice
                if((!(t_u > t_0)&&(!(t_u <= t_hn)))){  
                    N_i = N_i_0 * exp(0.2 * (t_0 - t_u)); // in 1/m³
                    m_i = min(r_humid.x[i][j][k] * ice.x[i][j][k]/N_i, m_i_max); // in kg
                    if(m_i > m_i_max) m_i = m_i_max;  // m_i_max = 1.0e-9, in kg
                    if(m_i < m_i_0) m_i = m_i_0;  // m_i_0 = 1.0e-12, in kg
/*
                cout.precision(15);
                cout.setf(ios::fixed);
                if((j == 124) && (k == 87)) cout << endl 
                    << "  TwoCategoryIceScheme ----- N_i and m_i ----- number density of ice particles and mass of cloud ice" << endl 
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  step = " << step[i] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl
                    << "  dt_rain_dim = " << dt_rain_dim
                    << "  dt_snow_dim = " << dt_snow_dim << endl
                    << "  number density of ice particles and mass of cloud ice" << endl
                    << "  N_i_0 = " << N_i_0 << "  N_i = " << N_i << endl
                    << "  m_i_0 = " << m_i_0 << "  m_i_max = " << m_i_max << "  m_s_0 = " << m_s_0 << "  m_i = " << m_i << endl << endl;
*/
                }
        // nucleation and depositional growth of cloud ice
                if(ice.x[i][j][k] == 0.0){
                    if(((t_u < t_d)&&(c.x[i][j][k] >= q_Ice))
                        ||(((t_d <= t_u)&&(t_u <= t_nuc)) 
                        &&(c.x[i][j][k] >= q_Rain)))
                        S_nuc = m_i_0/(r_humid.x[i][j][k] * dt_snow_dim) * N_i;  // nucleation of cloud ice, < I > in kg/(kg*s)
                }
                else  S_nuc = 0.0;
                if((t_u < t_hn)&&(cloud.x[i][j][k] > 0.0))  // happens only below -37°C
                    S_c_frz = cloud.x[i][j][k]/dt_rain_dim;  //nucleation of cloud ice due to freezing of cloud water, < II > in kg/(kg*s)
                else  S_c_frz = 0.0;
                if(t_u <= t_0){
                    if(c.x[i][j][k] > q_Ice)  // supersaturation
                        S_i_dep = c_i_dep * N_i * pow(m_i, (1.0/3.0)) // c_i_dep = 1.3e-5, in m3/s
                            * (c.x[i][j][k] - q_Ice);  // supersaturation, < III >
                    else  S_i_dep = 0.0;
                    if(- ice.x[i][j][k]/dt_snow_dim > (c.x[i][j][k] - q_Ice)/dt_snow_dim) 
                        S_i_dep = - ice.x[i][j][k]/dt_snow_dim; // subsaturation, < III >
                    else  S_i_dep = (c.x[i][j][k] - q_Ice)/dt_snow_dim; // subsaturation, < III >
                }
                else  S_i_dep = 0.0; // temperature > 0
    // autoconversion processes
                if((t_u >= t_0)&&(cloud.x[i][j][k] > 0.0))  // c_c_au = 4.0e-4, in 1/s
                    S_c_au = max(c_c_au * cloud.x[i][j][k], 0.0);  // cloud water to rain, cloud droplet collection, < IV >
                else  S_c_au = 0.0;
                if((t_u <= t_0)&&(ice.x[i][j][k] > 0.0))  // c_i_au = 1.0e-3, in 1/s
                    S_i_au = max(c_i_au * ice.x[i][j][k], 0.0);  // cloud ice to snow, cloud ice crystal aggregation, < V >
                else  S_i_au = 0.0;
                if(t_u <= t_0)
                    S_d_au = S_i_dep/(1.5 * (pow((m_s_0/m_i),  // m_s_0 = 3.0e-9, in kg
                        (2./3.)) - 1.0));  // autoconversion due to depositional growth of cloud ice, < VI >
                else  S_d_au = 0.0;
    // collection mechanism
                S_ac = c_ac * cloud.x[i][j][k]  // c_ac = 0.24, in m2/kg
//                    * pow(r_q_r,(7.0/9.0));
                    * pow(Rain,(7.0/9.0));
                    // accretion rate from depletion of cloud water due to collection by all rain drops, < VII >
                if(t_u < t_0)  
                    S_rim = c_rim * cloud.x[i][j][k] * Snow;  // c_rim = 18.6,  m2/kg
//                    S_rim = c_rim * cloud.x[i][j][k] 
//                        * pow(r_q_s,(13.0/12.0));  // c_rim = 18.6,  m2/kg
                else  S_rim = 0.0;  // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >
                                   // by falling snow particles
                if(t_u >= t_0)  
                    S_shed = c_rim * cloud.x[i][j][k] * Snow;
//                    S_shed = c_rim * cloud.x[i][j][k] 
//                        * pow(r_q_s,(13.0/12.0));
                else  S_shed = 0.0;  // rate of water shed by melting wet snow particles, < IX >
                                    // collecting cloud droplets to produce rain
                if(t_u <= t_0){
                    S_agg = c_agg * ice.x[i][j][k] * Snow; // collection of cloud ice by snow particles, < X >, c_agg = 10.3, m2/kg
//                    S_agg = c_agg * ice.x[i][j][k] 
//                        * pow(r_q_s,(13.0/12.0)); // collection of cloud ice by snow particles, < X >, c_agg = 10.3, m2/kg
                    S_i_cri = c_i_cri * ice.x[i][j][k] // c_i_cri = 0.24, m2
                        * pow(Rain,(7.0/9.0));
//                        * pow(r_q_r,(7.0/8.0));
                        // decrease in cloud ice mass due to collision/coalescense interaction with raindrops, < XI >
                    S_r_cri = c_r_cri * ice.x[i][j][k]/m_i // c_r_cri = 3.2e-5, m2
                        * pow(Rain,(13.0/9.0));
//                        * pow(r_q_r,(13.0/8.0));
                        // decrease of rainwater due to freezing resulting from collection of ice crystals, < XII >
                }else{
                    S_agg = 0.0;
                    S_i_cri = 0.0;
                    S_r_cri = 0.0;
                }
    // diffusional growth of rain and snow
                if(t_u < t_0){  // temperature below zero
                    S_s_dep = c_s_dep * (1.0 + b_s_dep // deposition/sublimation of snow, < XIV >
                        * pow(Snow,(5.0/26.0))) * (c.x[i][j][k] - q_Ice) // c_s_dep = 1.8e-2, m2/kg,
                        * pow(Snow,(8.0/13.0));
//                        * pow(r_q_s,(5.0/24.0))) * (c.x[i][j][k] - q_Ice) // c_s_dep = 1.8e-2, m2/kg,
//                        * pow(r_q_s,(2.0/3.0));
                    S_ev = 0.0;
                }else{  
                    S_ev = a_ev * (1.0 + b_ev * pow(Rain,(1.0/6.0)))  // evaporation of rain due to water vapour diffusion, < XIII >
                        * (q_Rain - c.x[i][j][k]) 
                        * pow(Rain,(4.0/9.0));
//                    S_ev = a_ev * (1.0 + b_ev * pow(r_q_r,(5.0/24.0)))  // evaporation of rain due to water vapour diffusion, < XIII >
//                        * (q_Rain - c.x[i][j][k]) 
//                        * pow(r_q_r,(4.0/9.0));
                    S_s_dep = 0.0;
                }
    // melting and freezing
                if(t_u > t_0){ // temperature above zero
                    if(ice.x[i][j][k] > 0.0)
                        S_i_melt = ice.x[i][j][k]/dt_snow_dim; // cloud ice particles melting to cloud water, < XV >
                    double p_t_in = p_hydro.x[i][j][k];
                    double E_Rain_t_in = hp * exp_func(t_0, 17.2694, 35.86);
                        // saturation water vapour pressure for the water phase at t = 0°C in hPa
                    double q_Rain_t_in = ep * E_Rain_t_in/(p_t_in 
                        - E_Rain_t_in);
                        // water vapour amount at saturation with water formation in kg/kg
                    S_s_melt = c_s_melt * (1.0 + b_s_melt   // melting rate of snow to form rain, < XVI >
                        * pow(Snow,(5.0/26.0))) *  // c_s_melt = 8.43e-5, (m2*s)/(K*kg)
//                        * pow(r_q_s,(5.0/24.0))) *  // c_s_melt = 8.43e-5, (m2*s)/(K*kg)
                        ((t_u - t_0) + a_s_melt * (c.x[i][j][k] 
                        - q_Rain_t_in)) * pow(Snow,(8.0/13.0));
//                        - q_Rain_t_in)) * pow(r_q_s,(2.0/3.0));
                }
                if(t_r_frz - t_u > 0.0){
                    double t_frz = max((t_r_frz - t_u), 0.0);
                    S_r_frz = c_r_frz * pow(t_frz, (3.0/2.0)) 
                        * pow(Rain,(3.0/2.0));
//                        * pow(r_q_r,(27.0/16.0));
                        // immersion freezing and contact nucleation, < XVII >
                }
                else  S_r_frz = 0.0;
    // sinks and sources
                S_v.x[i][j][k] = - S_c_c.x[i][j][k] + S_ev - S_i_dep 
                    - S_s_dep - S_nuc;
                S_c.x[i][j][k] = S_c_c.x[i][j][k] - S_c_au - S_ac 
                    - S_c_frz + S_i_melt - S_rim - S_shed;
                S_i.x[i][j][k] = S_nuc + S_c_frz + S_i_dep - S_i_melt 
                    - S_i_au - S_d_au - S_agg - S_i_cri;
                S_r.x[i][j][k] = S_c_au + S_ac - S_ev + S_shed 
                    - S_r_cri - S_r_frz + S_s_melt;
                S_s.x[i][j][k] = S_i_au + S_d_au + S_agg + S_rim 
                    + S_s_dep + S_i_cri + S_r_cri + S_r_frz - S_s_melt;

                if(is_land(h,i,j,k)){
                    S_c_c.x[i][j][k] = 0.0;
                    S_v.x[i][j][k] = 0.0;
                    S_c.x[i][j][k] = 0.0;
                    S_i.x[i][j][k] = 0.0;
                    S_r.x[i][j][k] = 0.0;
                    S_s.x[i][j][k] = 0.0;
                }
    // rain and snow integration
                if(t_u >= t_0)
/*
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * 0.5 * (S_r.x[i+1][j][k] 
                        + S_r.x[i][j][k]) * step[i];  // in kg/(m2 * s) == mm/s 
*/
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_r.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] >= 20.0/8.64e4)  P_rain.x[i][j][k] = 20.0/8.64e4;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] > 0.0){
                    rain = true;
                    if(P_rain.x[i][j][k] > maxValue_rain){
                        maxValue_rain = P_rain.x[i][j][k];
                        P_Rain = maxValue_rain * 8.64e4; // in mm/d
                        i_rain = i;
                        j_rain = j;
                        k_rain = k;
                        height_rain = get_layer_height(i);
                    }
                }
                if((t_u < t_0)&&(t_u >= t_00))
                    P_snow.x[i][j][k] = P_snow.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_s.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_snow.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] >= 1.0/8.64e4)  P_snow.x[i][j][k] = 1.0/8.64e4;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] > 0.0){
                    snow = true;
                    if(P_snow.x[i][j][k] > maxValue_snow){
                        maxValue_snow = P_snow.x[i][j][k];
                        P_Snow = maxValue_snow * 8.64e4; // in mm/d
                        i_snow = i;
                        j_snow = j;
                        k_snow = k;
                        height_snow = get_layer_height(i);
                    }
                }
/*
                cout.precision(10);
                cout.setf(ios::fixed);
                if((j == 124) && (k == 87)) cout << endl 
                    << "  TwoCategoryIceScheme ----- SINKS and SOURCES ----- " << endl 
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  step = " << step[i] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl << endl

                    << "  a)   nucleation and depositional growth of cloud ice" << endl
                    << "  Snuc = " << S_nuc * 1.0e6 
                    << "  Scfrz = " << S_c_frz * 1.0e6 
                    << "  Sidep = " << S_i_dep * 1.0e6 << endl << endl

                    << "  b)   autoconversion processes" << endl
                    << "  Siau = " << S_i_au * 1.0e6
                    << "  Scau = " << S_c_au * 1.0e6 
                    << "  Sdau = " << S_d_au * 1.0e6 << endl << endl

                    << "  c)   collection mechanism" << endl
                    << "  Sac = " << S_ac * 1.0e6
                    << "  Sagg = " << S_agg * 1.0e6
                    << "  Srim = " << S_rim * 1.0e6 
                    << "  Sshed = " << S_shed * 1.0e6
                    << "  Sagg = " << S_agg * 1.0e6
                    << "  Sicri = " << S_i_cri * 1.0e6
                    << "  Srcri = " << S_r_cri * 1.0e6 << endl << endl

                    << "  d)   diffusional growth of rain and snow" << endl
                    << "  Sev = " << S_ev * 1.0e6
                    << "  Ssdep = " << S_s_dep * 1.0e6 << endl << endl

                    << "  e)   melting and freezing" << endl
                    << "  Simelt = " << S_i_melt * 1.0e6 
                    << "  Ssmelt = " << S_s_melt * 1.0e6
                    << "  Srfrz = " << S_r_frz * 1.0e6 << endl << endl;
*/
/*
                cout.precision(10);
                cout.setf(ios::fixed);
                if((j == 124) && (k == 87)) cout << endl 
                    << "  TwoCategoryIceScheme ----- RAIN and SNOW PRODUCTION -----" << endl
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  step = " << step[i] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl
                    << "  c = " << c.x[i][j][k] * 1e3
                    << "  cl = " << cloud.x[i][j][k] * 1e3
                    << "  ice = " << ice.x[i][j][k] * 1e3 << endl
                    << "  r_dry = " << r_dry.x[i][j][k] 
                    << "  r_humid = " << r_humid.x[i][j][k] << endl
                    << "  p_hydro = " << p_hydro.x[i][j][k]
                    << "  p_stat = " << p_stat.x[i][j][k] << endl << endl

                    << "  a)   source and sink term due to condensation and evaporation" << endl
                    << "  S_c_c = " << S_c_c.x[i][j][k] * 1e3 << endl << endl

                    << "  b)   source and sink terms for the water categories: water vapour, cloud water and cloud ice" << endl
                    << "  S_v = - S_c_c + S_ev - S_i_dep - S_s_dep - S_nuc" << endl
                    << "  S_c = S_c_c - S_c_au - S_ac - S_c_frz + S_i_melt - S_rim - S_shed" << endl
                    << "  S_i = S_nuc + S_c_frz + S_i_dep - S_i_melt - S_i_au - S_d_au - S_agg - S_i_cri" << endl
                    << "  S_v = " << S_v.x[i][j][k] * 1.0e6 
                    << "  S_c = " << S_c.x[i][j][k] * 1.0e6 
                    << "  S_i = " << S_i.x[i][j][k] * 1.0e6 << endl << endl

                    << "  c)   source and sink terms for the water categories: rain and snow" << endl
                    << "  S_r = S_c_au + S_ac - S_ev + S_shed - S_r_cri - S_r_frz + S_s_melt" << endl
                    << "  S_s = S_i_au + S_d_au + S_agg + S_rim + S_s_dep + S_i_cri + S_r_cri + S_r_frz - S_s_melt" << endl
                    << "  S_r = " << S_r.x[i][j][k] * 1.0e6 
                    << "  S_s = " << S_s.x[i][j][k] * 1.0e6 << endl << endl

                    << "  d)   rain and snow formed by the sources and sinks in mm/d" << endl
                    << "  r_q_r = " << r_q_r * 8.64e4
                    << "  r_q_s = " << r_q_s * 8.64e4 << endl
                    << "  v_r_t = " << v_r_t
                    << "  v_s_t = " << v_s_t << endl
                    << "  A_r = " << A_r
                    << "  A_s = " << A_s
                    << "  B_r = " << B_r
                    << "  B_s = " << B_s << endl
                    << "  B_s = " << B_s << endl << endl

                    << "  P_rain = " << P_rain.x[i][j][k] * 8.64e4
                    << "  P_snow = " << P_snow.x[i][j][k] * 8.64e4
                    << "  P_conv_midl = " << P_conv_midl.x[i][j][k] * 8.64e4
                    << "  P_tot = " << P_rain.x[i][j][k] * 8.64e4 + P_snow.x[i][j][k] * 8.64e4 + P_conv_midl.x[i][j][k] * 8.64e4
                    << "  P_NASA = " << precipitation_NASA.y[j][k]  << endl << endl << endl;
*/
            }  // end i RainSnow
        }  // end j
    }  // end k

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            P_rain.x[im-1][j][k] = c43 * P_rain.x[im-2][j][k] 
                - c13 * P_rain.x[im-3][j][k];
            P_snow.x[im-1][j][k] = c43 * P_snow.x[im-2][j][k] 
                - c13 * P_snow.x[im-3][j][k];
                }
            }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            P_rain.x[i][0][k] = c43 * P_rain.x[i][1][k] 
                - c13 * P_rain.x[i][2][k];
            P_rain.x[i][jm-1][k] = c43 * P_rain.x[i][jm-2][k] 
                - c13 * P_rain.x[i][jm-3][k];
            P_snow.x[i][0][k] = c43 * P_snow.x[i][1][k] 
                - c13 * P_snow.x[i][2][k];
            P_snow.x[i][jm-1][k] = c43 * P_snow.x[i][jm-2][k] 
                - c13 * P_snow.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] 
                - c13 * P_rain.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] 
                - c13 * P_rain.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] 
                = (P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/2.0;
            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] 
                - c13 * P_snow.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] 
                - c13 * P_snow.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
        }
    }
/*
*
*/
    AtomUtils::fft_gaussian_filter_3d(P_rain,1);
    AtomUtils::fft_gaussian_filter_3d(P_snow,1);
/*
*
*/
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i < i_mount)){
                    P_rain.x[i][j][k] = P_rain.x[i_mount][j][k];
                    P_snow.x[i][j][k] = P_snow.x[i_mount][j][k];
                }
            }
        }
    }
    if(rain == false)  
        cout << "      no rain fall in TwoCategoryIceScheme found" 
        << endl;
    else
        cout << "      rain fall in TwoCategoryIceScheme found" 
        << endl
        << "      i_rain = " << i_rain
        << "      j_rain = " << j_rain
        << "   k_rain = " << k_rain
        << "   height_rain[m] = " << height_rain
        << "   P_Rain[mm/d] = " << P_Rain << endl;
    if(snow == false)  
        cout << "      no snow fall in TwoCategoryIceScheme found" 
        << endl;
    else
        cout << "      snow fall in TwoCategoryIceScheme found" 
        << endl
        << "      i_snow = " << i_snow
        << "      j_snow = " << j_snow
        << "   k_snow = " << k_snow
        << "   height_snow[m] = " << height_snow
        << "   P_Snow[mm/d] = " << P_Snow << endl;
    cout << "      TwoCategoryIceScheme ended" << endl;
    return;
}
/*
*
*/
// Three-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
// resulting the precipitation distribution formed of rain and snow
void cAtmosphereModel::ThreeCategoryIceScheme(){   
    cout << endl << "      ThreeCategoryIceScheme" << endl;
    // constant coefficients for the transport of cloud water, cloud ice and cloud graupel amount, 
    // rain snow and graupel in the parameterization procedures
    bool rain = false;
    bool snow = false;
    bool graupel = false;
    double P_Rain = 0.0;
    double P_Snow = 0.0;
    double P_Graupel = 0.0;
    double Rain = 0.0;
    double Snow = 0.0;
    int i_rain = 0;
    int i_snow = 0;
    int i_graupel = 0;
    int j_rain = 0;
    int j_snow = 0;
    int j_graupel = 0;
    int k_rain = 0;
    int k_snow = 0;
    int k_graupel = 0;
    double height_rain = 0.0;
    double height_snow = 0.0;
    double height_graupel = 0.0;
    double maxValue_rain = 0.0;
    double maxValue_snow = 0.0;
    double maxValue_graupel = 0.0;
    double N_i_0 = 1.0e2,  // in m-3
           N_g_0 = 4.0e6,  // in m-4
           m_i_0 = 1.0e-12,  // in kg
           m_i_max = 1.0e-9,  // in kg
           m_s_0 = 3.0e-9,  // in kg
           N_r_0 = 8.0e+6,  // in 1/m-4
           N_s_0 = 8.0e+5,  // in 1/m-4
           v_s_0 = 130.0,  // in m/s
           v_r_0 = 4.9,  // in m/s
           v_g_0 = 442.0,  // in m/s
           c_r_t = 12.63,  // in 
           c_s_t = 2.49,  // in 
           a_s_m = 0.038,  // in kg/m2
           a_g_m = 169.6,  // in kg/m2
           c_i_dep = 1.3e-5,  // in m3/s
           c_c_au = 4.0e-4,  // in 1/s
           c_i_au = 1.0e-3,  // in 1/s
           c_ac = 0.24,  // m2/kg
           c_rim = 18.6,  // m2/kg
           c_agg = 10.3,  // m2/kg
           c_i_cri = 0.24,  // m2
           c_r_cri = 3.2e-5,  // m2
           a_ev = 1.0e-3,  // m2/kg
           b_ev = 5.9,  // m2*s/kg
           c_r_frz = 3.75e-2,  // m2/(K*kg)
           t_nuc = 267.15,  // in K    -6 °C
           t_d = 248.15,  // in K    -25 °C
           t_hn = 236.15,  // in K    -37 °C
           t_r_frz = 271.15,  // in K    -2 °C
           z_csg = 0.5,
           c_g_rim = 4.43,
           c_g_agg = 2.46,
           c_s_1 = 0.0,
           c_s_2 = 0.0,
           c_s_3 = 0.0,
           c_s_4 = 0.0,
           c_g_1 = 0.0,
           c_g_2 = 0.0,
           c_g_3 = 0.0,
           c_g_4 = 0.0,
           t_crit = 0.0,
           d_v = 0.0,
           l_h = 0.0;
    double dt_snow_dim = 0.0,  // dt_snow_dim is the fallout velocity
           dt_rain_dim = 0.0;  // dt_rain_dim is the fallout velocity
    double t_u, p_u;
    double height;
    double q_Rain = 0.0,
           E_Rain = 0.0;
    double q_Ice = 0.0,
           E_Ice = 0.0;
    double m_i = m_i_max;  
    double A_r = r_0_water * M_PI * N_r_0,
           A_s = 2.0 * a_s_m * N_s_0, 
           A_g = 2.0 * a_g_m * N_g_0, 
           B_r = r_0_water * M_PI * N_r_0 * v_r_0 * tgamma(4.5)/tgamma(4.0), 
           B_s = N_s_0 * v_s_0 * a_s_m * tgamma(3.25),
           B_g = N_g_0 * v_g_0 * a_g_m * tgamma(3.25);
    double N_i, S_nuc, S_c_frz, S_i_dep, S_c_au, S_i_au, S_d_au, 
           S_ac, S_s_rim, S_g_rim, S_s_shed, S_g_shed;
    double S_s_agg, S_g_agg, S_i_cri, S_r_cri, S_ev, S_s_dep, S_i_melt, S_s_melt, S_r_frz;
    double S_csg, S_g_dep, S_g_melt;
    double r_q_r, r_q_s, r_q_g, v_r_t, v_s_t;
    std::vector<double> step(im, 0.0);
    // rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
    // the choosen scheme is a Two Category Ice Scheme
    // besides the transport equation for the water vapour exists two equations for the cloud water and the cloud ice transport
    // since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] < 0.0)  ice.x[i][j][k] = 0.0;
                if(gr.x[i][j][k] < 0.0)  gr.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                if(P_graupel.x[i][j][k] < 0.0)  P_graupel.x[i][j][k] = 0.0;
                P_rain.x[i][j][k] = 0.0;
                P_snow.x[i][j][k] = 0.0;
                P_graupel.x[i][j][k] = 0.0;
            }
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            P_rain.x[im-1][j][k] = 0.0;
            P_snow.x[im-1][j][k] = 0.0;
            P_graupel.x[im-1][j][k] = 0.0;
            P_rain.x[im-2][j][k] = 0.0;
            P_snow.x[im-2][j][k] = 0.0;
            P_graupel.x[im-2][j][k] = 0.0;
            S_r.x[im-1][j][k] = 0.0;
            S_s.x[im-1][j][k] = 0.0;
            S_g.x[im-1][j][k] = 0.0;
            S_r.x[im-2][j][k] = 0.0;
            S_s.x[im-2][j][k] = 0.0;
            S_g.x[im-2][j][k] = 0.0;
            for(int i = im-2; i >= 0; i--){
                height = get_layer_height(i);
                Rain = P_rain.x[i][j][k];
                Snow = P_snow.x[i][j][k];
                r_q_r = A_r * pow(B_r,-(8.0/9.0)) 
                    * pow(P_rain.x[i][j][k],(8.0/9.0));
                r_q_s = A_s * pow(B_s,-(12.0/13.0)) 
                    * pow(P_snow.x[i][j][k],(12.0/13.0));
                r_q_g = A_g * pow(B_g,-(12.0/13.0)) 
                    * pow(P_graupel.x[i][j][k],(12.0/13.0));
                v_r_t = c_r_t * pow(r_q_r, (1.0/8.0));
                v_s_t = c_s_t * pow(r_q_s, (1.0/12.0));
                t_u = t.x[i][j][k] * t_0;
                p_u = 1.0e2 * p_hydro.x[i][j][k];  // Pa
                step[i] = get_layer_height(i) - get_layer_height(i-1);  // local atmospheric shell thickness
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                dt_rain_dim = step[i]/1.6; // adjusted rain fall time step by fixed velocities == 1.6 m/s by variable local step size
                dt_snow_dim = step[i]/0.96; // adjusted snow fall time step by fixed velocities == 0.96 m/s by variable local step size
    // number density of ice particles and mass of cloud ice
                if((!(t_u > t_0)&&(!(t_u <= t_hn)))){  
                    N_i = N_i_0 * exp(0.2 * (t_0 - t_u)); // in 1/m³
                    m_i = min(r_humid.x[i][j][k] * ice.x[i][j][k]/N_i, m_i_max); // in kg
                    if(m_i > m_i_max) m_i = m_i_max;  // m_i_max = 1.0e-9, in kg
                    if(m_i < m_i_0) m_i = m_i_0;  // m_i_0 = 1.0e-12, in kg
                }
    // heterogeneous nucleation and depositional growth of cloud ice
                if(ice.x[i][j][k] == 0.0){
                    if(((t_u < t_d)&&(c.x[i][j][k] >= q_Ice))
                        ||(((t_d <= t_u)&&(t_u <= t_nuc)) 
                        &&(c.x[i][j][k] >= q_Rain)))
                        S_nuc = m_i_0/(r_humid.x[i][j][k] * dt_snow_dim) * N_i;  // nucleation of cloud ice, < I > in kg/(kg*s)
                }
                else  S_nuc = 0.0;
                if((t_u < t_hn)&&(cloud.x[i][j][k] > 0.0))  // happens only below -37°C
                    S_c_frz = cloud.x[i][j][k]/dt_rain_dim;  //nucleation of cloud ice due to freezing of cloud water, < II > in kg/(kg*s)
                else  S_c_frz = 0.0;
    // deposition growth and sublimation of cloud ice
                if(t_u <= t_0){
                    if(c.x[i][j][k] > q_Ice)  // supersaturation
                        S_i_dep = c_i_dep * N_i * pow(m_i, (1.0/3.0))
                            * (c.x[i][j][k] - q_Ice);
                    else  S_i_dep = 0.0;
                    if(- ice.x[i][j][k]/dt_snow_dim > (c.x[i][j][k] - q_Ice)/dt_snow_dim) 
                        S_i_dep = - ice.x[i][j][k]/dt_snow_dim; // subsaturation, < III >
                    else  S_i_dep = (c.x[i][j][k] - q_Ice)/dt_snow_dim; // subsaturation, < III >
                }
                else  S_i_dep = 0.0; // temperature > 0
    // autoconversion processes
    // autoconversion of cloud water to form rain
                if((t_u >= t_0)&&(cloud.x[i][j][k] > 0.0))
//                    S_c_au = max(c_c_au * (cloud.x[i][j][k] - 0.0002), 0.0);  // cloud water to rain, cloud droplet collection
                    S_c_au = max(c_c_au * cloud.x[i][j][k], 0.0);  // cloud water to rain, cloud droplet collection
                else  S_c_au = 0.0;
    // autoconversion of cloud ice to form snow due to aggregation
                if((t_u <= t_0)&&(ice.x[i][j][k] > 0.0))
                    S_i_au = max(c_i_au * ice.x[i][j][k], 0.0);  // cloud ice to snow, cloud ice crystal aggregation
                else  S_i_au = 0.0;
    // autoconversion of cloud ice to form snow due to deposition
                if(t_u <= t_0)
                    S_d_au = S_i_dep/(1.5 * (pow((m_s_0/m_i),
                        (2./3.)) - 1.0));  // autoconversion due to depositional growth of cloud ice
                else  S_d_au = 0.0;
    // collection mechanism
    // accretion of cloud water by raindrops
                if(t_u >= t_0)  
                    S_ac = c_ac * cloud.x[i][j][k]  // accretion rate from depletion of cloud water due to collection by all rain drops
                        * pow(Rain,(7.0/9.0));
                else  S_ac = 0.0;
    // collection of cloud water by snow or graupel (riming)
                if(t_u < t_0){  
                    S_s_rim = c_rim * cloud.x[i][j][k] * Snow;  // c_rim = 18.6,  m2/kg
                    S_g_rim = c_rim * cloud.x[i][j][k] 
                        * pow((r_humid.x[i][j][k] * r_q_g), 0.94878);  // c_rim = 18.6,  m2/kg
                }else{
                    S_s_rim = 0.0;
                    S_g_rim = 0.0;
                }
    // collection of cloud water by wet snow or graupel to form rain (shedding)
                if(t_u >= t_0){
                    S_s_shed = c_rim * cloud.x[i][j][k] * Snow;
                    S_g_shed = c_g_rim * cloud.x[i][j][k] 
                        * pow((r_humid.x[i][j][k] * r_q_g), 0.94878);
                }else{
                    S_s_shed = 0.0;
                    S_g_shed = 0.0;
                }
    // collection of cloud ice by snow or graupel
                if(t_u <= t_0){
                    S_s_agg = c_agg * ice.x[i][j][k] * Snow; // collection of cloud ice by snow particles, < X >, c_agg = 10.3, m2/kg
                    S_g_agg = c_g_agg * ice.x[i][j][k] // collection of cloud ice by graupel particles
                        * pow((r_humid.x[i][j][k] * r_q_g), 0.94878);
    // collection of cloud ice by rain to form graupel
                    S_i_cri = c_i_cri * ice.x[i][j][k]
                        * pow(Rain,(7.0/9.0));
    // freezing of rain due to collection of cloud ice to form graupel (aggregation)
                    S_r_cri = c_r_cri * ice.x[i][j][k]/m_i
                        * pow(Rain,(13.0/9.0));
                }else{
                    S_s_agg = 0.0;
                    S_g_agg = 0.0;
                    S_i_cri = 0.0;
                    S_r_cri = 0.0;
                }
    // evaporation of rain water
                if(t_u < t_0)  // temperature below zero
                    S_ev = 0.0;
                else
                    S_ev = a_ev * (1.0 + b_ev * pow(Rain,(1.0/6.0)))  // evaporation of rain due to water vapour diffusion, < XIII >
                        * (q_Rain - c.x[i][j][k]) 
                        * pow(Rain,(4.0/9.0));
    // deposition growth and sublimation of snow or graupel
                if(t_u >= t_0)
                    d_v = 101325.0/p_u 
                        * (2.22e-5 + 1.46e-7 * (t_u - t_0));
                else
                    d_v = 101325.0/p_u 
                        * (2.22e-5 + 1.25e-7 * (t_u - t_0));
                l_h = 0.024 + 8.0e-5 * (t_u - t_0);
                double a_melt = r_humid.x[i][j][k] * lv * d_v/l_h;
                double E_Rain_t_0 = hp * exp_func(t_0, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                double q_Rain_t_0 = ep * E_Rain_t_0
                           /(1.0e-2 * p_u - E_Rain_t_0); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                if(c.x[i][j][k] >= q_Rain_t_0)
                    t_crit = t_0 - 1.0/l_h * ls * d_v * r_humid.x[i][j][k] 
                        * (c.x[i][j][k] - q_Rain_t_0);
                else  t_crit = t_0;
                if(t_u < t_0){  // temperature below zero
                    c_s_1 = 2.91955;
                    c_s_2 = 0.0109928;
                    c_s_3 = 15871.3;
                    c_s_4 = 1.74744e-6;
                    c_g_1 = 0.398561;
                    c_g_2 = 0.00152398;
                    c_g_3 = 2554.99;
                    c_g_4 = 2.6531e-7;
                    S_s_dep = 
                        (c_s_1 - c_s_2 * t_u + c_s_3/p_u + c_s_4 * p_u) 
                        * (c.x[i][j][k] - q_Ice) 
                        * pow((r_humid.x[i][j][k] * r_q_s), 0.8); // deposition/sublimation of snow
                    S_g_dep = 
                        (c_g_1 - c_g_2 * t_u + c_g_3/p_u + c_g_4 * p_u) 
                        * (c.x[i][j][k] - q_Ice) 
                        * pow((r_humid.x[i][j][k] * r_q_g), 0.6); // deposition/sublimation of graupel
                }
                if((t_0 < t_u)&&(t_u < t_crit)){  // temperature greater zero and less the critical temperature
                    c_s_1 = 0.28003;
                    c_s_4 = 0.146293e-6;
                    c_g_1 = 0.0418521;
                    c_g_4 = 4.7524e-8;
                    S_s_dep = 
                        (c_s_1 - c_s_4 * p_u) 
                        * (c.x[i][j][k] - q_Ice) 
                        * pow((r_humid.x[i][j][k] * r_q_s), 0.8); // deposition/sublimation of snow, < XIV >
                    S_g_dep = 
                        (c_g_1 - c_g_4 * p_u) 
                        * (c.x[i][j][k] - q_Ice) 
                        * pow((r_humid.x[i][j][k] * r_q_g), 0.6); // deposition/sublimation of snow, < XIV >
                }
                if((t_0 < t_crit)&&(t_crit < t_u)){  // critical temperature greater zero and less temperature
                    c_s_1 = 2.41897,
                    c_s_3 = 31282.3,
                    c_g_1 = 0.153907,
                    c_g_4 = 7.86703e-7;
                    S_s_dep = 
                        (c_s_1 + c_s_3/p_u) 
                        * (c.x[i][j][k] - q_Ice) 
                        * pow((r_humid.x[i][j][k] * r_q_s), 0.8); // deposition/sublimation of snow, < XIV >
                    S_g_dep = 
                        (c_g_1 - c_g_4 * p_u) 
                        * (c.x[i][j][k] - q_Ice) 
                        * pow((r_humid.x[i][j][k] * r_q_g), 0.6); // deposition/sublimation of snow, < XIV >
                }
    // melting of snow or graupel to form cloud water
                if(t_u >= t_0){  // temperature above zero
                    c_s_1 = 0.612654e-3;
                    c_s_3 = 79.6863,
                    c_g_1 = 7.39441e-5,
                    c_g_3 = 12.31698,
                    a_melt = 2.95e+3;
                    S_s_melt = ((t_u - t_0) + a_melt 
                        * (c.x[i][j][k] - q_Rain_t_0)) 
                        * (c_s_1 + c_s_3/p_u)
                        * pow((r_humid.x[i][j][k] * r_q_s), 0.8);
                    S_g_melt = ((t_u - t_0) + a_melt 
                        * (c.x[i][j][k] - q_Rain_t_0)) 
                        * (c_g_1 + c_g_3/p_u)
                        * pow((r_humid.x[i][j][k] * r_q_g), 0.6);
                }
    // melting of cloud ice to form cloud water
                S_i_melt = 0.0;
                if(t_u > t_0){ // temperature above zero
                    if(ice.x[i][j][k] > 0.0)
                        S_i_melt = ice.x[i][j][k]/dt_snow_dim;
                }
    // freezing of rain to form graupel
                if(t_r_frz - t_u > 0.0){
//                    S_r_frz = c_r_frz * pow((t_r_frz - t_u), (3.0/2.0)) 
//                        * pow(r_q_r, (27.0/16.0));
                    double t_frz = max((t_r_frz - t_u), 0.0);
                    S_r_frz = c_r_frz * pow(t_frz, (3.0/2.0)) 
                        * pow(r_q_r, (27.0/16.0));
//                        * pow(Rain, (27.0/16.0));
                }else  S_r_frz = Rain/dt_rain_dim;
    // conversion of snow to graupel due to riming
                S_csg = 0.0;
                if(cloud.x[i][j][k] > 0.0002){ // in kg/kg
                    S_csg = z_csg * cloud.x[i][j][k] 
                    * pow((r_humid.x[i][j][k] * r_q_s), 0.75);
//                    * pow(Snow, 0.75);
                }
    // collection mechanism
                if(t_u < t_0)  
                    S_g_rim = c_rim * cloud.x[i][j][k] 
                        * pow((r_humid.x[i][j][k] * gr.x[i][j][k]), 0.94878);  // c_rim = 18.6,  m2/kg
                else  S_g_rim = 0.0;  // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >
                                   // by falling snow particles
                if(t_u >= t_0)  
                    S_g_shed = c_g_rim * cloud.x[i][j][k] 
                    * pow((r_humid.x[i][j][k] * gr.x[i][j][k]), 0.94878);
                else  S_g_shed = 0.0;  // rate of water shed by melting wet snow particles, < IX >
                                    // collecting cloud droplets to produce rain
                if(t_u <= t_0){
                    S_g_agg = c_g_agg * ice.x[i][j][k] 
                    * pow((r_humid.x[i][j][k] * gr.x[i][j][k]), 0.94878); // collection of cloud ice by snow particles, < X >, c_agg = 10.3, m2/kg
                }else{
                    S_g_agg = 0.0;
                }
    // sinks and sources
                S_v.x[i][j][k] = - S_c_c.x[i][j][k] + S_ev - S_i_dep 
                    - S_s_dep - S_g_dep - S_nuc;
                S_c.x[i][j][k] = S_c_c.x[i][j][k] - S_c_au - S_ac 
                    - S_c_frz + S_i_melt - S_s_rim - S_g_rim - S_s_shed 
                    - S_g_shed;
                S_i.x[i][j][k] = S_nuc + S_c_frz + S_i_dep - S_i_melt 
                    - S_i_au - S_d_au - S_s_agg - S_g_agg - S_i_cri;
                S_r.x[i][j][k] = S_c_au + S_ac - S_ev + S_s_shed 
                    + S_g_shed - S_r_cri - S_r_frz + S_s_melt + S_g_melt;
                S_s.x[i][j][k] = S_i_au + S_d_au + S_s_agg + S_s_rim 
                    + S_s_dep + S_i_cri + S_r_cri 
                    - S_s_melt - S_csg;
                S_g.x[i][j][k] = S_g_agg + S_g_rim + S_g_dep + S_i_cri 
                    + S_r_cri + S_r_frz - S_g_melt + S_csg;
                if(is_land(h,i,j,k)){
                    S_c_c.x[i][j][k] = 0.0;
                    S_v.x[i][j][k] = 0.0;
                    S_c.x[i][j][k] = 0.0;
                    S_i.x[i][j][k] = 0.0;
                    S_g.x[i][j][k] = 0.0;
                    S_r.x[i][j][k] = 0.0;
                    S_s.x[i][j][k] = 0.0;
                }
    // rain and snow integration
                if(t_u >= t_0)
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_r.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] >= 10.0/8.64e4)  P_rain.x[i][j][k] = 10.0/8.64e4;  // (mm/d)/(24h * 60min * 60 *s) == mm/s
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] > 0.0){
                    rain = true;
                    if(P_rain.x[i][j][k] > maxValue_rain){
                        maxValue_rain = P_rain.x[i][j][k];
                        P_Rain = maxValue_rain * 8.64e4; // in mm/d
                        i_rain = i;
                        j_rain = j;
                        k_rain = k;
                        height_rain = get_layer_height(i);
                    }
                }
                if((t_u < t_0)&&(t_u >= t_00))
                    P_snow.x[i][j][k] = P_snow.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_s.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_snow.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] >= 1.0/8.64e4)  P_snow.x[i][j][k] = 1.0/8.64e4;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] > 0.0){
                    snow = true;
                    if(P_snow.x[i][j][k] > maxValue_snow){
                        maxValue_snow = P_snow.x[i][j][k];
                        P_Snow = maxValue_snow * 8.64e4; // in mm/d
                        i_snow = i;
                        j_snow = j;
                        k_snow = k;
                        height_snow = get_layer_height(i);
                    }
                }
                if(t_u >= t_00)
                    P_graupel.x[i][j][k] = P_graupel.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_g.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_graupel.x[i][j][k] = 0.0;
                if(P_graupel.x[i][j][k] >= 1.0/8.64e4)  P_graupel.x[i][j][k] = 1.0/8.64e4;
                if(P_graupel.x[i][j][k] < 0.0)  P_graupel.x[i][j][k] = 0.0;
                if(P_graupel.x[i][j][k] > 0.0){
                    graupel = true;
                    if(P_graupel.x[i][j][k] > maxValue_graupel){
                        maxValue_graupel = P_graupel.x[i][j][k];
                        P_Graupel = maxValue_graupel * 8.64e4; // in mm/d
                        i_graupel = i;
                        j_graupel = j;
                        k_graupel = k;
                        height_graupel = get_layer_height(i);
                    }
                }
/*
                cout.precision(10);
                cout.setf(ios::fixed);
                if((j == 35) && (k == 87)) cout << endl 
                    << "  ThreeCategoryIceScheme ----- SINKS and SOURCES ----- " << endl 
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  step = " << step[i] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl << endl

                    << "  a)   nucleation and depositional growth of cloud ice" << endl
                    << "  Snuc = " << S_nuc * 1.0e6 
                    << "  Scfrz = " << S_c_frz * 1.0e6 
                    << "  Sidep = " << S_i_dep * 1.0e6 << endl << endl

                    << "  b)   autoconversion processes" << endl
                    << "  Siau = " << S_i_au * 1.0e6
                    << "  Scau = " << S_c_au * 1.0e6 
                    << "  Sdau = " << S_d_au * 1.0e6 
                    << "  Scsg = " << S_csg * 1.0e6 << endl << endl

                    << "  c)   collection mechanism" << endl
                    << "  Sac = " << S_ac * 1.0e6
                    << "  Ssagg = " << S_s_agg * 1.0e6
                    << "  Sgagg = " << S_g_agg * 1.0e6
                    << "  Ssrim = " << S_s_rim * 1.0e6 
                    << "  Sgrim = " << S_g_rim * 1.0e6 << endl
                    << "  Ssshed = " << S_s_shed * 1.0e6
                    << "  Sgshed = " << S_g_shed * 1.0e6
                    << "  Ssagg = " << S_s_agg * 1.0e6
                    << "  Sgagg = " << S_g_agg * 1.0e6 << endl
                    << "  Sicri = " << S_i_cri * 1.0e6
                    << "  Srcri = " << S_r_cri * 1.0e6 << endl << endl

                    << "  d)   diffusional growth of rain, snow and graupel" << endl
                    << "  Sev = " << S_ev * 1.0e6
                    << "  Ssdep = " << S_s_dep * 1.0e6 << endl
                    << "  Sgdep = " << S_g_dep * 1.0e6 << endl << endl

                    << "  e)   melting and freezing" << endl
                    << "  Simelt = " << S_i_melt * 1.0e6 
                    << "  Ssmelt = " << S_s_melt * 1.0e6
                    << "  Sgmelt = " << S_g_melt * 1.0e6
                    << "  Srfrz = " << S_r_frz * 1.0e6 << endl << endl;
*/
/*
                cout.precision(10);
                cout.setf(ios::fixed);
                if((j == 35) && (k == 87)) cout << endl 
                    << "  ThreeCategoryIceScheme ----- RAIN, SNOW and GRAUPEL PRODUCTION -----" << endl
                    << "  Ma = " << (int)*get_current_time() << endl
                    << "  i = " << i << "  j = " << j << "  k = " << k << endl
                    << "  step = " << step[i]
                    << "  height = " << height << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl
                    << "  d_v = " << d_v
                    << "  l_h = " << l_h
                    << "  t_crit = " << t_crit - t_0 << endl
                    << "  c = " << c.x[i][j][k] * 1e3
                    << "  cl = " << cloud.x[i][j][k] * 1e3
                    << "  ice = " << ice.x[i][j][k] * 1e3
                    << "  q_Rain_t_0 = " << q_Rain_t_0 * 1e3 << endl
                    << "  r_dry = " << r_dry.x[i][j][k] 
                    << "  r_humid = " << r_humid.x[i][j][k] << endl
                    << "  p_hydro = " << p_hydro.x[i][j][k]
                    << "  p_stat = " << p_stat.x[i][j][k] << endl << endl

                    << "  a)   source and sink term due to condensation and evaporation" << endl
                    << "  S_c_c = " << S_c_c.x[i][j][k] * 1.0e6 << endl << endl

                    << "  b)   source and sink terms for the water categories: water vapour, cloud water, cloud ice and cloud graupel" << endl
                    << "  S_v = - S_c_c + S_ev - S_i_dep - S_s_dep - S_g_dep - S_nuc" << endl
                    << "  S_c = S_c_c - S_c_au - S_ac - S_c_frz + S_i_melt - S_s_rim - S_g_rim - S_s_shed - S_g_shed" << endl
                    << "  S_i = S_nuc + S_c_frz + S_i_dep - S_i_melt - S_i_au - S_d_au - S_s_agg - S_g_agg - S_i_cri" << endl
                    << "  S_s = S_i_au + S_d_au + S_s_agg + S_s_rim + S_s_dep + S_i_cri + S_r_cri - S_s_melt - S_csg" << endl
                    << "  S_v = " << S_v.x[i][j][k] * 1.0e6 
                    << "  S_c = " << S_c.x[i][j][k] * 1.0e6 
                    << "  S_i = " << S_i.x[i][j][k] * 1.0e6
                    << "  S_g = " << S_g.x[i][j][k] * 1.0e6 << endl << endl

                    << "  c)   source and sink terms for the water categories: rain, snow and graupel" << endl
                    << "  S_r = S_c_au + S_ac - S_ev + S_s_shed + S_g_shed - S_r_cri - S_r_frz + S_s_melt + S_g_melt" << endl
                    << "  S_s = S_i_au + S_d_au + S_agg + S_rim + S_s_dep + S_i_cri + S_r_cri + S_r_frz - S_s_melt" << endl
                    << "  S_g = S_g_agg + S_g_rim + S_g_dep + S_i_cri + S_r_cri + S_r_frz - S_g_melt + S_csg" << endl
                    << "  S_r = " << S_r.x[i][j][k] * 1.0e6 
                    << "  S_s = " << S_s.x[i][j][k] * 1.0e6 
                    << "  S_g = " << S_g.x[i][j][k] * 1.0e6 << endl << endl

                    << "  d)   rain and snow formed by the sources and sinks in mm/d" << endl
                    << "  r_q_r = " << r_q_r * 8.64e4
                    << "  r_q_s = " << r_q_s * 8.64e4
                    << "  r_q_g = " << r_q_s * 8.64e4 << endl
                    << "  v_r_t = " << v_r_t
                    << "  v_s_t = " << v_s_t << endl
                    << "  A_r = " << A_r
                    << "  A_s = " << A_s
                    << "  B_r = " << B_r
                    << "  B_s = " << B_s << endl << endl

                    << "  P_rain = " << P_rain.x[i][j][k] * 8.64e4
                    << "  P_snow = " << P_snow.x[i][j][k] * 8.64e4
                    << "  P_graupel = " << P_graupel.x[i][j][k] * 8.64e4 << endl
                    << "  P_conv_midl = " << P_conv_midl.x[i][j][k] * 8.64e4
                    << "  P_conv_shal = " << P_conv_shall.x[i][j][k] * 8.64e4 << endl
                    << "  P_tot = " << P_rain.x[i][j][k] * 8.64e4 + P_snow.x[i][j][k] * 8.64e4 + P_conv_midl.x[i][j][k] * 8.64e4
                    << "  P_NASA = " << precipitation_NASA.y[j][k]  << endl << endl << endl;
*/
            }  // end i RainSnowGraupel
        }  // end j
    }  // end k
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            P_rain.x[im-1][j][k] = c43 * P_rain.x[im-2][j][k] 
                - c13 * P_rain.x[im-3][j][k];
            P_snow.x[im-1][j][k] = c43 * P_snow.x[im-2][j][k] 
                - c13 * P_snow.x[im-3][j][k];
            P_graupel.x[im-1][j][k] = c43 * P_graupel.x[im-2][j][k] 
                - c13 * P_graupel.x[im-3][j][k];
                }
            }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            P_rain.x[i][0][k] = c43 * P_rain.x[i][1][k] 
                - c13 * P_rain.x[i][2][k];
            P_rain.x[i][jm-1][k] = c43 * P_rain.x[i][jm-2][k] 
                - c13 * P_rain.x[i][jm-3][k];
            P_snow.x[i][0][k] = c43 * P_snow.x[i][1][k] 
                - c13 * P_snow.x[i][2][k];
            P_snow.x[i][jm-1][k] = c43 * P_snow.x[i][jm-2][k] 
                - c13 * P_snow.x[i][jm-3][k];
            P_graupel.x[i][0][k] = c43 * P_graupel.x[i][1][k] 
                - c13 * P_graupel.x[i][2][k];
            P_graupel.x[i][jm-1][k] = c43 * P_graupel.x[i][jm-2][k] 
                - c13 * P_graupel.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] 
                - c13 * P_rain.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] 
                - c13 * P_rain.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] 
                = (P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/2.0;
            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] 
                - c13 * P_snow.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] 
                - c13 * P_snow.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
            P_graupel.x[i][j][0] = c43 * P_graupel.x[i][j][1] 
                - c13 * P_graupel.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_graupel.x[i][j][km-1] = c43 * P_graupel.x[i][j][km-2] 
                - c13 * P_graupel.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
        }
    }
/*
*
*/
    AtomUtils::fft_gaussian_filter_3d(P_rain,1);
    AtomUtils::fft_gaussian_filter_3d(P_snow,1);
    AtomUtils::fft_gaussian_filter_3d(P_graupel,1);
/*
*
*/
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i < i_mount)){
                    P_rain.x[i][j][k] = P_rain.x[i_mount][j][k];
                    P_snow.x[i][j][k] = P_snow.x[i_mount][j][k];
                    P_graupel.x[i][j][k] = P_graupel.x[i_mount][j][k];
                }
            }
        }
    }
    if(rain == false)  
        cout << "      no rain fall in ThreeCategoryIceScheme found" 
        << endl;
    else
        cout << "      rain fall in ThreeCategoryIceScheme found" 
        << endl
        << "      i_rain = " << i_rain
        << "      j_rain = " << j_rain
        << "   k_rain = " << k_rain
        << "   height_rain[m] = " << height_rain
        << "   P_Rain[mm/d] = " << P_Rain << endl;
    if(snow == false)  
        cout << "      no snow fall in ThreeCategoryIceScheme found" 
        << endl;
    else
        cout << "      snow fall in ThreeCategoryIceScheme found" 
        << endl
        << "      i_snow = " << i_snow
        << "      j_snow = " << j_snow
        << "   k_snow = " << k_snow
        << "   height_snow[m] = " << height_snow
        << "   P_Snow[mm/d] = " << P_Snow << endl;
    if(graupel == false)  
        cout << "      no graupel fall in ThreeCategoryIceScheme found" 
        << endl;
    else
        cout << "      graupel fall in ThreeCategoryIceScheme found" 
        << endl
        << "      i_graupel = " << i_graupel
        << "      j_graupel = " << j_graupel
        << "   k_graupel = " << k_graupel
        << "   height_graupel[m] = " << height_graupel
        << "   P_Graupel[mm/d] = " << P_Graupel << endl;
    cout << "      ThreeCategoryIceScheme ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::MoistConvectionMidL(){
    cout << endl << "      MoistConvectionMidL" << endl;
// collection of coefficients for phase transformation
    int i_base = 0;
    int ii_base_beg = 1;
    int i_lfs = 0;
    double a_ev = 1.0e-3;  // m2/kg
    double b_ev = 5.9;  // m2*s/kg
    double R_W_R_A = R_WaterVapour/R_Air;
    double eps = 1.0/R_W_R_A;
    double alf = (1.0 - eps)/eps;
    double t_00 = 236.15; // = -37°C
    bool moistconv = false;
    double P_moistconv = 0.0;
    int i_moistconv = 0;
    int j_moistconv = 0;
    int k_moistconv = 0;
    double height_moistconv = 0.0;
    double maxValue = 0.0;
/*
// constants used for the COSMO parameterization
    double bet = 42.0; // in K given by COSMO
    double b_u = 0.3;
    double alf_1 = 5.0e-4; // in 1/s
    double alf_2 = 0.011; // in 1/s
    double C_p = 0.05; // in ./.
    double bet_p = 2.0e-3;  // in 1/s, from the COSMO paper
    double eps_u = 1.0e-4;  // in 1/m
    double eps_d = 2.0e-4;
    double del_u = eps_u;
    double del_d = eps_d;
    double delta_i_c = 0.0;
    double K_p = 0.0;  // in 1/s
*/
    double gam_d = - 0.3;  // from original paper by Tiedtke, gam_d = - 0.2, COSMO approach
    double a_u = 1.0e0;  // area fraction in updraft to define the vertical velocity
    double a_d = 1.0e0;  // area fraction in downdraft to define the vertical velocit
    double HumidityRelative = 0.0;
    double L_latent = 0.0;
    double q_buoy_base = 0.0;
    double q_buoy_lfs = 0.0;
    double E_Rain = 0.0;
    double q_Rain = 0.0;
    double e = 0.0;
//    double E_Ice = 0.0;
    double E_Rain_add = 0.0;
    double E_Rain_base = 0.0;
    double q_Rain_add = 0.0;
    double q_Rain_base = 0.0;
//    double E_Ice_add = 0.0;
//    double q_Ice = 0.0;
//    double q_Ice_add = 0.0;
//    double r_humid_add = 0.0;
//    double r_dry_add = 0.0;
//    double p_stat_add = 0.0;
    double r_humid_add = 0.0;
//    double r_dry_add = 0.0;
    double t_u = 0.0;
    double t_u_add = 0.0;
    double t_u_base = 0.0;
    double height = 0.0;
//    double height_base = 0.0;
//    double height_lfs = 0.0;
    double rm = 0.0;
    double exp_rm = 0.0;
    double M_u_denom = 0.0;
    double M_d_denom = 0.0;
    double M_u_num = 0.0;
    double M_d_num = 0.0;
    double t_vir = 0.0;
    double t_vir_env = 0.0;
    double t_add_u = 0.000732; //  in K == + 0.2 °C temperature increase for a lifting air parcel, given by ECMWF = European Center for Medium-Range Weather Forecast
//    double t_add_u = 0.0; //  in K == + 0.2 °C temperature increase for a lifting air parcel, given by ECMWF = European Center for Medium-Range Weather Forecast
    double q_v_u_add = 1.0e-4; // == 0.0001 kg/kg increase for a lifting air parcel, given by ECMWF
//    double q_v_u_add = 0.5e-4; // == 0.00005 kg/kg increase for a lifting air parcel, given by ECMWF
//    double q_v_u_add = 0.0; // == 0.001 kg/kg increase for a lifting air parcel, testing
    std::vector<double> step(im, 0);
    std::vector<std::vector<double> > M_u_Base(jm, std::vector<double> (km, 0));
    std::vector<std::vector<double> > M_d_LFS(jm, std::vector<double> (km, 0));
    std::vector<std::vector<double> > u_cloud_base(jm, std::vector<double> (km, 0));
    rad.Coordinates(im, r0, dr);
/*
    double dcdr = 0.0, dcdthe = 0.0, dcdphi = 0.0; // COSMO data
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            int i_mount = i_topography[j][k];
            sinthe = sin(the.z[j]);
            for(int i = im-1; i >= i_mount ; i--){
                rm = rad.z[i] * get_layer_height(i);
                rmsinthe = rm * sinthe;
                exp_rm = 1.0/(rm + 1.0);
                dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k])/(2. * dr) * exp_rm;
                if(i < im-2){
                    if((is_land(h, i, j, k))&&((is_air(h, i+1, j, k))
                        &&(is_air(h, i+2, j, k)))){
                        dcdr = (- 3. * c.x[i][j][k] + 4. * c.x[i+1][j][k] 
                            - c.x[i+2][j][k])/(2. * dr ) * exp_rm;
                    }
                }
                if(c.x[i][j][k] == 0.0) c.x[i][j][k] = c_tropopause;
                dcdr = 0.0; 
                dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k])
                    /(2.0 * dthe * rm);
                dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1])
                    /(2.0 * dphi * rmsinthe);

                if((j >= 2)&&(j < jm-3)){
                    if((is_land(h, i, j, k))&&((is_air(h, i, j+1, k))
                        &&(is_air(h, i, j+2, k))))
                        dcdthe = (- 3. * c.x[i][j][k] + 4. * c.x[i][j+1][k] 
                            - c.x[i][j+2][k])/(2. * dthe * rm);
                    if((is_land(h, i, j, k))&&(is_air(h, i, j-1, k))
                        &&(is_air(h, i, j-2, k)))
                        dcdthe = (- 3. * c.x[i][j][k] + 4. * c.x[i][j-1][k] 
                            - c.x[i][j-2][k])/(2. * dthe * rm);
                }
                if((k >= 2)&&(k < km-3)){
                    if((is_land(h, i, j, k))&&(is_air(h, i, j, k+1))
                        &&(is_air(h, i, j, k+2)))
                        dcdphi = (- 3. * c.x[i][j][k] + 4. * c.x[i][j][k+1] 
                            - c.x[i][j][k+2])/(2. * dphi * rmsinthe);
                    if((is_land(h, i, j, k))&&(is_air(h, i, j, k-1))
                        &&(is_air(h, i, j, k-2)))
                        dcdphi = (- 3. * c.x[i][j][k] + 4. * c.x[i][j][k-1] 
                            - c.x[i][j][k-2])/(2. * dphi * rmsinthe);
                }
                if((j >= 2)&&(j < jm-3)){
                    if((is_land(h, i, j, k))&&((is_air(h, i, j+1, k))
                        &&(is_air(h, i, j+2, k))))
                        dcdthe = (c.x[i][j+1][k] - c.x[i][j][k])/(dthe * rm);
                    if((is_land(h, i, j, k))&&(is_air(h, i, j-1, k))
                        &&(is_air(h, i, j-2, k)))
                        dcdthe = (c.x[i][j][k] - c.x[i][j-1][k])/(dthe * rm);
                }
                if((k >= 2)&&(k < km-3)){
                    if((is_land(h, i, j, k))&&(is_air(h, i, j, k+1))
                        &&(is_air(h, i, j, k+2)))
                        dcdphi = (c.x[i][j][k+1] - c.x[i][j][k])/(dphi * rmsinthe);
                    if((is_land(h, i, j, k))&&(is_air(h, i, j, k-1))
                        &&(is_air(h, i, j, k-2)))
                        dcdphi = (c.x[i][j][k] - c.x[i][j][k-1])/(dphi * rmsinthe);
                }
// specification of entrainment in the updraft
                if(c.x[i][j][k] == 0.0) c.x[i][j][k] = 1e+6;
                E_u.x[i][j][k] = - r_humid.x[i][j][k]/c.x[i][j][k] 
                    * (u.x[i][j][k] * dcdr 
                    + v.x[i][j][k] * dcdthe 
                    + w.x[i][j][k] * dcdphi) * u_0; // in kg/(m³s)
                if(E_u.x[i][j][k] >= 0.2) E_u.x[i][j][k] = 0.2;
                if(E_u.x[i][j][k] <= - 0.2) E_u.x[i][j][k] = - 0.2;
            }
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            E_u.x[0][j][k] = c43 * E_u.x[1][j][k] -
                c13 * E_u.x[2][j][k];
            E_u.x[im-1][j][k] = c43 * E_u.x[im-2][j][k] -
                c13 * E_u.x[im-3][j][k];
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            E_u.x[i][0][k] = c43 * E_u.x[i][1][k] -
                 c13 * E_u.x[i][2][k];
            E_u.x[i][jm-1][k] = c43 * E_u.x[i][jm-2][k] -
                 c13 * E_u.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            E_u.x[i][j][0] = c43 * E_u.x[i][j][1] -
                c13 * E_u.x[i][j][2];
            E_u.x[i][j][km-1] = c43 * E_u.x[i][j][km-2] -
                c13 * E_u.x[i][j][km-3];
            E_u.x[i][j][0] = E_u.x[i][j][km-1] =
               (E_u.x[i][j][0] + E_u.x[i][j][km-1])/2.;
        }
    }
    AtomUtils::fft_gaussian_filter_3d(E_u,1);
*/
// parameterisation of moist convection
    CAPE = std::vector<double>(im, 0.0);
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
// boundary conditions on the base of clouds (level of free convection)
            int i_mount = i_topography[j][k];
            if(i_mount >= 0) 
                ii_base_beg = i_mount;
            else ii_base_beg = 0;
            for(int ii_base = ii_base_beg; ii_base < im-1; ii_base++){
                t_u = t.x[ii_base][j][k] * t_0; // in K
                t_u_add = t_add_u * t_0 + t_u; // in K
                q_v_u.x[ii_base][j][k] = c.x[ii_base][j][k] + q_v_u_add;
                q_c_u.x[ii_base][j][k] = cloud.x[ii_base][j][k];

                r_humid_add = 1e2 * p_hydro.x[ii_base][j][k]
                    /(R_Air * (1. + (R_W_R_A - 1.) * q_v_u.x[ii_base][j][k] 
                    - q_c_u.x[ii_base][j][k]) * t_u_add); 

                t_vir = t_u_add * (1.0 + alf * q_v_u.x[ii_base][j][k] 
                    - q_c_u.x[ii_base][j][k]);
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Rain_add = hp * exp_func(t_vir, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_hydro.x[ii_base][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Rain_add = ep * E_Rain_add/(p_hydro.x[ii_base][j][k] - E_Rain_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                e = q_v_u.x[ii_base][j][k] * p_hydro.x[ii_base][j][k]/ep;  // water vapour pressure in hPa
                HumidityRelative = e/E_Rain_add * 100.0;
                q_buoy_base = r_humid_add - r_humid.x[ii_base][j][k];
/*
                cout.precision(6);
                if((j == 90) &&(k == 180))  cout << endl 
                    << "  MoistConvectionMidL   Cloud Base" << endl 
                    << "  Ma = " << Ma << endl
                    << "  i_base UPDRAFT" << endl 
                    << "  j = " << j << "  k = " << k << endl 
                    << "  ii_base = " << ii_base 
                    << "  r_dry.x[ii_base] = " << r_dry.x[ii_base][j][k] 
                    << "  r_humid.x[ii_base] = " << r_humid.x[ii_base][j][k] 
                    << "  r_humid_add = " << r_humid_add << endl
                    << "  E_Rain = " << E_Rain 
                    << "  E_Rain_add = " << E_Rain_add 
                    << "  q_Rain = " << q_Rain 
                    << "  q_Rain_add = " << q_Rain_add << endl
                    << "  HumidityRelative = " << HumidityRelative
                    << "  q_buoy_base = " << q_buoy_base  << endl << endl << endl;
*/
                if((HumidityRelative >= 80.0)&&(q_buoy_base >= 0.0)  // lifting convection level (LCL)
                    &&(height >= 500.0)){ // search for cloud base, height threshold of 200m prescribed by ECMWF for midlevel convection
                    i_Base.y[j][k] = (double)ii_base;
                    i_base = ii_base;
                    height = get_layer_height(i_base);
                    if(is_land(h, i_base, j, k)){
                       int i_mount = i_topography[j][k];
                       u.x[i_base][j][k] = 0.0;
                       v.x[i_base][j][k] = 0.0;
                       w.x[i_base][j][k] = 0.0;
                       c.x[i_base][j][k] = c.x[i_mount][j][k];
                       cloud.x[i_base][j][k] = cloud.x[i_mount][j][k];
                       ice.x[i_base][j][k] = ice.x[i_mount][j][k];
                    }
                    t_u = t.x[i_base][j][k] * t_0; // in K
                    t_u_add = t_add_u * t_0 + t_u; // in K
                    r_humid_add = 1e2 * p_hydro.x[i_base][j][k]
                        /(R_Air * (1. + (R_W_R_A - 1.) * q_v_u.x[i_base][j][k] 
                        - q_c_u.x[i_base][j][k]) * t_u_add); 
//                    if(iter_cnt == -1)  u_cloud_base[j][k] = u.x[i_base][j][k];
                    u.x[i_base][j][k] = u.x[i_base][j][k];
                    if(iter_cnt == -1)  u_cloud_base[j][k] = u.x[i_base][j][k];
                    M_u.x[i_base][j][k] = r_humid_add 
                        * u.x[i_base][j][k] * u_0;  // in kg/(m²s) [== mm/s] updraft at cloud base 
//                        * u_cloud_base[j][k] * u_0;  // in kg/(m²s) [== mm/s] updraft at cloud base 
//                    if(M_u.x[i_base][j][k] < 0.0) M_u.x[i_base][j][k] = 0.0;
                    M_d_LFS[j][k] = gam_d * M_u.x[i_base][j][k];  // in kg/(m²s) [== mm/s] downdraft at cloud top 
                    if(is_land(h, i_base, j, k))  
                        M_u.x[i_base][j][k] = M_d_LFS[j][k] = 0.0;
                    s.x[i_base][j][k] = (cp_l * t_u + g * height)/s_0;
                    q_v_u.x[i_base][j][k] = c.x[i_base][j][k] + q_v_u_add;
                    q_c_u.x[i_base][j][k] = cloud.x[i_base][j][k];
                    u_u.x[i_base][j][k] = u.x[i_base][j][k];
                    v_u.x[i_base][j][k] = v.x[i_base][j][k];
                    w_u.x[i_base][j][k] = w.x[i_base][j][k];
                    s_u.x[i_base][j][k] = (cp_l * t_vir + g * height)/s_0;
                    g_p.x[i_base][j][k] = 0.0;
                    c_u.x[i_base][j][k] = 0.0;
                    e_l.x[i_base][j][k] = 0.0;
                    e_p.x[i_base][j][k] = 0.0;
                    rm = rad.z[i_base];
                    exp_rm = 1./(rm + 1.0);
                    t_vir = t_u_add * (1.0 + alf * q_v_u.x[i_base][j][k] 
                        - q_c_u.x[i_base][j][k]);
                    t_vir_env = t_u * (1.0 + alf * c.x[i_base][j][k] 
                        - cloud.x[i_base][j][k]);
                    step[i_base] = get_layer_height(i_base+1) 
                        - get_layer_height(i_base);  // in m local atmospheric shell thickness
                    CAPE[i_base] = g * step[i_base]
                        /exp_rm * (t_vir - t_vir_env)/t_vir_env;
/*
                    cout.precision(6);
                    if((j == 90) &&(k == 180))  cout << endl 
                        << "  MoistConvectionMidL   Cloud Base    +++++++++++++++++++++++++++++++++++++++++++++++++++" << endl 
                        << "  Ma = " << Ma << endl
                        << "  i_base UPDRAFT" << endl 
                        << "  j = " << j << "  k = " << k << endl 
                        << "  height_base = " << height << endl
                        << "  ii_base = " << ii_base 
                        << "  i_base = " << (int)i_Base.y[j][k]
                        << "  i_lfs = " << (int)i_LFS.y[j][k] << endl
                        << "  i_tropo = " << i_topography[j][k] << endl
                        << "  s = " << s.x[i_base][j][k] 
                        << "  s_u = " << s_u.x[i_base][j][k] << endl
                        << "  E_Rain = " << E_Rain 
                        << "  E_Rain_add = " << E_Rain_add 
                        << "  q_Rain = " << q_Rain 
                        << "  q_Rain_add = " << q_Rain_add << endl
                        << "  q_buoy_base = " << q_buoy_base 
                        << "  q_buoy_lfs = " << q_buoy_lfs 
                        << "  c = " << c.x[i_base][j][k] 
                        << "  q_v_u = " << q_v_u.x[i_base][j][k] 
                        << "  cloud = " << cloud.x[i_base][j][k] 
                        << "  q_c_u = " << q_c_u.x[i_base][j][k] << endl
                        << "  u = " << u.x[i_base][j][k] 
                        << "  u_u = " << u_u.x[i_base][j][k] 
                        << "  v = " << v.x[i_base][j][k] 
                        << "  v_u = " << v_u.x[i_base][j][k] 
                        << "  w = " << w.x[i_base][j][k] 
                        << "  w_u = " << w_u.x[i_base][j][k] << endl
                        << "  t_u = " << t_u - t_0 
                        << "  t_u_add = " << t_add_u * t_u - t_0 << endl
                        << "  HumidityRelative = " << HumidityRelative 
                        << "  r_dry = " << r_dry.x[i_base][j][k] 
                        << "  r_dry_add = " << r_dry_add 
                        << "  r_humid = " << r_humid.x[i_base][j][k] 
                        << "  r_humid_add = " << r_humid_add << endl
                        << "  p_stat = " << p_stat.x[i_base][j][k] 
                        << "  p_hydro = " << p_hydro.x[i_base][j][k] 
                        << "  M_u = " << M_u.x[i_base][j][k] 
                        << "  M_d_LFS = " << M_d_LFS[j][k] << endl
                        << "  rad = " << rad.z[i_base]
                        << "  step_non_dim = " << g * (rad.z[i_base] - rad.z[i_base-1])/exp_rm
                        << "  temp_non_dim = " << (t_u_add - t_u)/t_u
                        << "  CAPE = " << CAPE[i_base] << endl << endl;
*/
                    break;
                } // i_base-loop
            } // ii_base-loop
// level of free sinking
/*
*
*/
// boundary conditions on top of the cloud
            for(int ii_lfs = im-2; ii_lfs > 0; ii_lfs--){
                height = get_layer_height(ii_lfs); // in m
                t_u = t.x[ii_lfs][j][k] * t_0; // in K
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_hydro.x[ii_lfs][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                t_u_add = t_add_u * t_0 + t_u; // in K
                r_humid_add = 1e2 * p_hydro.x[ii_lfs][j][k]/(R_Air 
                    * (1. + (R_W_R_A - 1.) * c.x[ii_lfs][j][k] 
                    - cloud.x[ii_lfs][j][k] - ice.x[ii_lfs][j][k]) 
                    * t_u_add); 
                q_v_d.x[ii_lfs][j][k] = c.x[ii_lfs][j][k];
                q_buoy_lfs = 0.5 * (q_c_u.x[ii_lfs][j][k] + q_Rain) 
                    - c.x[ii_lfs][j][k];
//                q_buoy_lfs = r_humid_add - r_humid.x[ii_lfs][j][k];
                if((q_buoy_lfs <= 0.0)&&(t_u >= t_00)){ // search for cloud top
                    i_LFS.y[j][k] = (double)ii_lfs;
                    i_lfs = ii_lfs;
                    height = get_layer_height(i_lfs); // in m
                    if((t.x[i_lfs][j][k] * t_0) >= t_0) 
                        L_latent = lv;
                    else L_latent = ls;
                    s.x[i_lfs][j][k] = (cp_l * t_u + g * height)/s_0;
                    q_v_d.x[i_lfs][j][k] = c.x[i_lfs][j][k];
                    M_d.x[i_lfs][j][k] = M_d_LFS[j][k];  // in kg/(m²s) [== mm/s] downdraft at cloud top
                    u_d.x[i_lfs][j][k] = u.x[i_lfs][j][k] 
                        + M_d.x[i_lfs][j][k]/(a_d * r_humid.x[i_lfs][j][k]); 
                    v_d.x[i_lfs][j][k] = v.x[i_lfs][j][k];
                    w_d.x[i_lfs][j][k] = w.x[i_lfs][j][k];
                    s_d.x[i_lfs][j][k] = s.x[i_lfs][j][k];
                    e_p.x[i_lfs][j][k] = 0.0;
                    e_d.x[i_lfs][j][k] = 0.0;
/*
                    cout.precision(6);
                    if((j == 90) &&(k == 180))  cout << endl 
                        << "  MoistConvectionMidL   Cloud Top     ---------------------------------------------"  << endl 
                        << "  Ma = " << Ma << endl
                        << "  i_lfs DOWNDRAFT" << endl
                        << "  j = " << j<< "  k = " << k << endl 
                        << "  height_lfs = " << height << endl
                        << "  ii_lfs = " << ii_lfs 
                        << "  i_base = " << (int)i_Base.y[j][k] 
                        << "  i_lfs = " << (int)i_LFS.y[j][k] << endl
                        << "  s = " << s.x[i_lfs][j][k] 
                        << "  s_d = " << s_d.x[i_lfs][j][k] << endl
                        << "  E_Rain = " << E_Rain 
                        << "  q_Rain = " << q_Rain << endl
                        << "  q_buoy_base = " << q_buoy_base 
                        << "  c = " << c.x[i_lfs][j][k] 
                        << "  q_v_d = " << q_v_d.x[i_lfs][j][k] 
                        << "  cloud = " << cloud.x[i_lfs][j][k] << endl
                        << "  u = " << u.x[i_lfs][j][k] 
                        << "  u_d = " << u_d.x[i_lfs][j][k] 
                        << "  v = " << v.x[i_lfs][j][k] 
                        << "  v_d = " << v_d.x[i_lfs][j][k] 
                        << "  w = " << w.x[i_lfs][j][k] 
                        << "  w_d = " << w_d.x[i_lfs][j][k] << endl
                        << "  t = " << t.x[i_lfs][j][k] * t_0 - t_0 << endl
                        << "  HumidityRelative = " << HumidityRelative 
                        << "  r_dry = " << r_dry.x[i_lfs][j][k] 
                        << "  r_humid = " << r_humid.x[i_lfs][j][k] << endl
                        << "  p_stat = " << p_stat.x[i_lfs][j][k] 
                        << "  p_hydro = " << p_hydro.x[i_lfs][j][k] 
                        << "  M_d = " << M_d.x[i_lfs][j][k] 
                        << "  M_d_LFS = " << M_d_LFS[j][k] << endl << endl;
*/
                    break;
                }
            } // ii_lfs-loop
/*
*
*/
// zero moist values below i_base
            for(int i = 0; i < i_base; i++){
                t_u = t.x[i][j][k] * t_0; // in K
                if(t_u >= t_0) L_latent = lv;
                else L_latent = ls;
                height = get_layer_height(i);
                s.x[i][j][k] = (cp_l * t_u + g * height)/s_0;
                s_u.x[i][j][k] = 0.0;
                s_d.x[i][j][k] = 0.0;
                q_v_u.x[i][j][k] = 0.0;
                q_c_u.x[i][j][k] = 0.0;
                u_u.x[i][j][k] = 0.0;
                v_u.x[i][j][k] = 0.0;
                w_u.x[i][j][k] = 0.0;
                g_p.x[i][j][k] = 0.0;
                c_u.x[i][j][k] = 0.0;
                e_l.x[i][j][k] = 0.0;
            } // zero moist values below i_base
// zero moist values above i_lfs
            for(int i = i_lfs+1; i < im; i++){
                t_u = t.x[i][j][k] * t_0; // in K
                height = get_layer_height(i);
                s.x[i][j][k] = (cp_l * t_u + g * height)/s_0;
                s_u.x[i][j][k] = 0.0;
                s_d.x[i][j][k] = 0.0;
                q_v_d.x[i][j][k] = 0.0;
                u_d.x[i][j][k] = 0.0;
                v_d.x[i][j][k] = 0.0;
                w_d.x[i][j][k] = 0.0;
            } // zero moist values above i_lfs
        }  // end j for cloud base and cloud top
    }  // end k for cloud base and cloud top
/*
*
*/
// ECMWF data (European Center of Mid-Range Weather Forecast)
    double gam_d_K = 0.5;
    double c_d = 0.506;
    double bet_d = 1.875;
    double f_t = 2.0;
    double del_u_1 = 0.75e-4; // 1/m
    double K_u_sqrt = 0.0;
    double K_d_sqrt = 0.0;
    double D_u_2 = 0.0;
    double del_M_u = 0.0;
    double D_d_2 = 0.0;
    double del_M_d = 0.0;
    double eps_u_deep = 1.75e-3; // 1/m
    double alf_11 = 5.44e-4; // s
    double alf_22 = 5.09e-3;
    double alf_33 = 0.577;
    double RH_crit_water = 0.9; // over water
    double RH_crit_land = 0.7; // over land
    double sig = 0.05;
//
// CONVECTIVE UPDRAFT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            i_base = (int)i_Base.y[j][k];
            i_lfs = (int)i_LFS.y[j][k];
            K_u = std::vector<double>(im, 0.0);
            K_u[i_base] = 0.5 * u_u.x[i_base][j][k] * u_u.x[i_base][j][k]; // kinetic energy
            for(int i = i_base; i < im-1; i++){
                height = get_layer_height(i);
                step[i] = get_layer_height(i+1) - get_layer_height(i);  // in m local atmospheric shell thickness
                t_u = t.x[i][j][k] * t_0; // in K
                t_u_add = t_add_u * t_0 + t_u; // in K
                t_vir = t_u_add * (1.0 + alf * q_v_u.x[i][j][k] 
                    - q_c_u.x[i][j][k]);
                t_vir_env = t_u * (1.0 + alf * c.x[i][j][k] 
                    - cloud.x[i][j][k]);
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
//                E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                E_Rain_add = hp * exp_func(t_vir, 17.2694, 35.86);
//                E_Ice_add = hp * exp_func(t_u, 21.8746, 7.66);
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Rain_add = ep * E_Rain_add/(p_hydro.x[i][j][k] - E_Rain_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
//                q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
//                q_Ice_add = ep * E_Ice_add/(p_hydro.x[i][j][k] - E_Ice_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                t_u_base = t.x[i_base][j][k] * t_0;
                E_Rain_base = hp * exp_func(t_u_base, 17.2694, 35.86);
                q_Rain_base = ep * E_Rain_base
                    /(p_hydro.x[i_base][j][k] - E_Rain_base); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                e = c.x[i][j][k] * p_hydro.x[i][j][k]/ep;  // water vapour pressure in hPa
                HumidityRelative = e/E_Rain * 100.0;
                r_humid_add = 1e2 * p_hydro.x[i][j][k]
                    /(R_Air * (1.0 + (R_W_R_A - 1.) * q_v_u.x[i][j][k] 
                    - q_c_u.x[i][j][k]) * t_u_add); 
                q_buoy_base = r_humid_add - r_humid.x[i][j][k];
                if(step[i] >= t_0) L_latent = lv;
                else L_latent = ls;
                s.x[i][j][k] = (cp_l * t_u + g * height)/s_0;
/*
// COSMO data
                eps_u = 1.e-4;
                del_u = 1.e-4;
                if(t_u_add < t_00){
                    eps_u = 0.0;
                    del_u = 0.0;
                }
*/
// CAPE (Convective Available Potential Energy) == difference between the moist adiabatic and the environmental temperature
                rm = rad.z[i];
                exp_rm = 1./(rm + 1.0);
                CAPE[i+1] = CAPE[i] + g * step[i]
                    /exp_rm * (t_vir - t_vir_env)/t_vir_env;
// water vapour in the updraft
                if(q_v_u.x[i][j][k] >= q_Rain_add){
                    q_c_u.x[i][j][k] = q_v_u.x[i][j][k] - q_Rain_add;
                    if(q_c_u.x[i][j][k] >= 0.0) 
                        q_v_u.x[i][j][k] = q_Rain_add;
                }
// condensation within the updraft
                if((c.x[i][j][k] - q_Rain_add) >= 0.0)
                    c_u.x[i][j][k] = (c.x[i][j][k] - q_Rain_add) 
                        * M_u.x[i][j][k]/(r_humid.x[i][j][k] * step[i]);  // in kg/(kg*s)
// specification of entrainment in the updraft
//                E_u.x[i][j][k] = E_u.x[i][j][k] // COSMO data
//                    + eps_u * M_u.x[i][j][k]; // in kg/(m³s)
                E_u.x[i][j][k] = eps_u_deep * M_u.x[i][j][k]
                    /r_humid.x[i][j][k] * (1.3 - HumidityRelative/100.0) 
                    * pow(q_Rain/q_Rain_base, 3.0);

                if(E_u.x[i][j][k] > 0.01) E_u.x[i][j][k] = 0.01;
                if(E_u.x[i][j][k] < - 0.01) E_u.x[i][j][k] = - 0.01;

// specification of detrainment in the updraft
/*
// COSMO data
                if(i == i_lfs)
                    D_u.x[i][j][k] = (1. - b_u) 
                        * M_u.x[i][j][k]/step[i]
                        + del_u * M_u.x[i][j][k];
                if(i == i_lfs+1)
                    D_u.x[i][j][k] = b_u 
                        * M_u.x[i][j][k]/step[i]
                        + del_u * M_u.x[i][j][k];  // in kg/(m³s)
*/
// ECMWF approach
                if(q_buoy_base < 0.0){ // search for non-buoyant situation
                    if(M_u.x[i][j][k] == 0.0)  M_u.x[i][j][k] = 1e-6;
                    K_u[i+1] = K_u[i] + step[i]  // kinetic energy
                        * (- E_u.x[i][j][k]/M_u.x[i][j][k] 
                        * (1.0 + bet_d * c_d) 
                        * 2.0 * K_u[i] + g/(f_t * (1.0 + gam_d_K)) 
                        * (t_u_add - t_u)/t_u);
                    K_u_sqrt = sqrt(K_u[i]/K_u[i+1]);
                    del_M_u = (1.6 - HumidityRelative/100.0) * K_u_sqrt;
                    D_u_2 = M_u.x[i][j][k]/(r_humid.x[i][j][k] 
                        * step[i]) * (1.0/del_M_u - 1.0);
                }else{
                    K_u[i] = 0.0;
                    K_u[i+1] = 0.0;
                    K_u_sqrt = 0.0;
                    del_M_u = 0.0;
                    D_u_2 = 0.0;
                }
                D_u.x[i][j][k] = del_u_1 * M_u.x[i][j][k]  // in kg/(m³s)
                    /r_humid.x[i][j][k] + D_u_2;
                if(t_u < t_00)  D_u.x[i][j][k] = 0.0;

                if(D_u.x[i][j][k] > 0.01) D_u.x[i][j][k] = 0.01;
                if(D_u.x[i][j][k] < - 0.01) D_u.x[i][j][k] = - 0.01;

// evaporation of cloud water in the updraft
                e_l.x[i][j][k] = D_u.x[i][j][k]/r_humid.x[i][j][k]
                    * q_c_u.x[i][j][k];  // in kg/(kg*s)
// formation of precipitation within the updraft
/*
// COSMO data
                if(is_land(h, 0, j, k)) delta_i_c = 3000.0; // over land in m
                if(is_air(h, 0, j, k))  delta_i_c = 1500.0; // over ocean in m
                height = get_layer_height(i); // in m
                height_base = get_layer_height(i_base); // in m
                if(height > height_base + delta_i_c)  K_p = bet_p;  // in 1/s
                else  K_p = 0.0;
 */
                if(t_u_add >= t_00){
//                    g_p.x[i][j][k] = K_p * q_c_u.x[i][j][k];  // in kg/(kg*s)// COSMO approach
                    double c_00 = 1.5e-3; // 1/s
                    double q_krit = 5.0e-4; // kg/kg
                    double w_u_u = 10.0; // m/s

                    double w_u_g = 0.0;
                    w_u_g = fabs(u_u.x[i][j][k]) * u_0;
                    if(w_u_g >= w_u_u)  w_u_g = w_u_u;
                    if(w_u_g == 0.0)  w_u_g = 1.0e-6;
                    if(w_u_g <= 1.0)  w_u_g = 1.0;
//                    if(w_u_g <= 0.16)  w_u_g = 0.16;

                    g_p.x[i][j][k] = M_u.x[i][j][k]/r_humid.x[i][j][k]  // original ECMWF formula
                        * (c_00/(0.75 * w_u_u)) 
//                        * (c_00/(0.75 * w_u_g))
                        * q_c_u.x[i][j][k] 
                        * (1.0 - exp(- pow(q_c_u.x[i][j][k]/q_krit, 2.0)));  // in kg/(kg*s)
                }
                else  g_p.x[i][j][k] = 0.0;
//                if(g_p.x[i][j][k] < 0.0)  g_p.x[i][j][k] = 0.0;
// mass flux within the updraft
                M_u.x[i+1][j][k] = M_u.x[i][j][k] 
                    - step[i] * (E_u.x[i][j][k] - D_u.x[i][j][k]);    // in kg/(m²s) integration from cloud base upwards
                M_u_denom = M_u.x[i+1][j][k];
                M_u_num = M_u.x[i][j][k] - step[i] * D_u.x[i][j][k];
                if(M_u_denom == 0.0)  M_u_denom = 1e-6;

                if(M_u.x[i+1][j][k] > 5.0) M_u.x[i+1][j][k] = 5.0;
                if(M_u.x[i+1][j][k] < - 5.0) M_u.x[i+1][j][k] = - 5.0;

// dry entropy within the updraft
                s_u.x[i+1][j][k] = (M_u_num * s_u.x[i][j][k] 
                    - step[i] * (E_u.x[i][j][k] * s.x[i][j][k] 
                    - D_u.x[i][j][k] * s_u.x[i][j][k]
                    + L_latent * r_humid.x[i][j][k] * c_u.x[i][j][k])/s_0)
                    /M_u_denom;
// cloud water within the updraft
                q_v_u.x[i+1][j][k] = (M_u_num * q_v_u.x[i][j][k]  // in kg/kg
                    - step[i] * (E_u.x[i][j][k] * c.x[i][j][k] 
                    - r_humid.x[i][j][k] * c_u.x[i][j][k]))/M_u_denom;
                if(q_v_u.x[i+1][j][k] >= q_Rain_add)
                    q_v_u.x[i+1][j][k] = q_Rain_add;
// cloud water within the updraft
                q_c_u.x[i+1][j][k] = (M_u_num * q_c_u.x[i][j][k]   // in kg/kg
                    - step[i] * r_humid.x[i][j][k] 
                    * (c_u.x[i][j][k] - g_p.x[i][j][k]))/M_u_denom;
// velocity components within the updraft
                u_u.x[i+1][j][k] = u.x[i][j][k] 
                    + M_u.x[i][j][k]/(a_u * r_humid.x[i][j][k]); 
                if(u_u.x[i+1][j][k] > 2.5) u_u.x[i+1][j][k] = 2.5;
                if(u_u.x[i+1][j][k] < - 2.5) u_u.x[i+1][j][k] = - 2.5;
                v_u.x[i+1][j][k] = (M_u_num * v_u.x[i][j][k] 
                    - step[i] * E_u.x[i][j][k] * v.x[i][j][k])/M_u_denom;  // in m/s
                if(v_u.x[i+1][j][k] > 0.1) v_u.x[i+1][j][k] = 0.1;
                if(v_u.x[i+1][j][k] < - 0.1) v_u.x[i+1][j][k] = - 0.1;
                w_u.x[i+1][j][k] = (M_u_num * w_u.x[i][j][k] 
                    - step[i] * E_u.x[i][j][k] * w.x[i][j][k])/M_u_denom;  // in m/s
                if(w_u.x[i+1][j][k] > 10.0) w_u.x[i+1][j][k] = 10.0;
                if(w_u.x[i+1][j][k] < -10.0) w_u.x[i+1][j][k] = - 10.0;
//                if(t_u <= t_00){
                if(is_land(h, i, j, k)||(t_u <= t_00)){
                    E_u.x[i+1][j][k] = 0.0;
                    M_u.x[i+1][j][k] = 0.0;
                    s_u.x[i+1][j][k] = 0.0;
                    q_v_u.x[i+1][j][k] = 0.0;
                    q_c_u.x[i+1][j][k] = 0.0;
                    u_u.x[i+1][j][k] = 0.0;
                    v_u.x[i+1][j][k] = 0.0;
                    w_u.x[i+1][j][k] = 0.0;
                    g_p.x[i][j][k] = 0.0;
                    c_u.x[i][j][k] = 0.0;
                    e_l.x[i][j][k] = 0.0;
                }
/*    cout.precision(9);
            if((j == 90) &&(k == 180))  cout << endl 
                << "  MoistConvectionMidL   UPDRAFT     +++++++++++++++++++++++++++++++++++++++++++++++" << endl
                << "  Ma = " << (int)*get_current_time() << endl
                << "   i = " << i 
                << "  j = " << j << "  k = " << k << endl 
                << "   i_base = " << (int)i_Base.y[j][k] 
                << "   i_lfs = " << (int)i_LFS.y[j][k] << endl
                << "  s = " << s.x[i][j][k] 
                << "  s_u = " << s_u.x[i][j][k] << endl
                << "  E_u = " << E_u.x[i][j][k] 
                << "  D_u = " << D_u.x[i][j][k] << endl
                << "  e_l = " << e_l.x[i][j][k] 
                << "  c_u = " << c_u.x[i][j][k] 
                << "  g_p = " << g_p.x[i][j][k] << endl
                << "  E_Rain = " << E_Rain 
                << "  E_Rain_add = " << E_Rain_add 
                << "  q_Rain = " << q_Rain 
                << "  q_Rain_add = " << q_Rain_add << endl
                << "  q_buoy_base = " << q_buoy_base 
                << "  c = " << c.x[i][j][k] 
                << "  q_v_u = " << q_v_u.x[i][j][k] 
                << "  cloud = " << cloud.x[i][j][k] 
                << "  q_c_u = " << q_c_u.x[i][j][k] << endl
                << "  u = " << u.x[i][j][k] 
                << "  u_u = " << u_u.x[i][j][k] 
                << "  v = " << v.x[i][j][k] 
                << "  v_u = " << v_u.x[i][j][k] 
                << "  w = " << w.x[i][j][k] 
                << "  w_u = " << w_u.x[i][j][k] << endl
                << "  t_u = " << t_u - t_0 
                << "  t_u_add = " << t_u_add - t_0
                << "  t_vir = " << t_vir - t_0 << endl
                << "  HumidityRelative = " << HumidityRelative 
                << "  r_dry = " << r_dry.x[i][j][k] 
                << "  r_humid = " << r_humid.x[i][j][k] 
                << "  r_humid_add = " << r_humid_add << endl
                << "  p_stat = " << p_stat.x[i][j][k] 
                << "  p_hydro = " << p_hydro.x[i][j][k] 
                << "  height = " << height 
//                << "  height_base = " << height_base 
                << "  step = " << step[i] << endl
                << "  M_u = " << M_u.x[i][j][k] 
                << "  M_d_LFS = " << M_d_LFS[j][k] << endl
                << "  rad = " << rad.z[i]
                << "  step_non_dim = " << g * (rad.z[i+1] - rad.z[i])/exp_rm
                << "  temp_non_dim = " << (t_u_add - t_u)/t_u
                << "  CAPE = " << CAPE[i] << endl
                << "  del_M_u = " << del_M_u
                << "  D_u_2 = " << D_u_2
                << "  D_u = " << D_u.x[i][j][k]
                << "  K_u = " << K_u[i] << endl << endl;
*/
            } // end i convective updraft
        }  // end j convective updraft
    }  // end k convective updraft
//
// CONVECTIVE DOWNDRAFT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            i_base = (int)i_Base.y[j][k];
            i_lfs = (int)i_LFS.y[j][k];
            M_d.x[i_lfs][j][k] = M_d_LFS[j][k];
//            height_lfs = get_layer_height(i_lfs);
            P_conv_midl.x[i_lfs][j][k] = 0.0;
            t_u_add = t_add_u * t.x[i_base][j][k] * t_0; // in K
            E_Rain_add = hp * exp_func(t_u_add, 17.2694, 35.86);
            K_d = std::vector<double>(im, 0.0);
            K_d[i_lfs] = 0.5 * u_d.x[i_lfs][j][k] * u_d.x[i_lfs][j][k];
            for(int i = i_lfs; i > 0; i--){
                height = get_layer_height(i);
                step[i] = get_layer_height(i) - get_layer_height(i-1);  // local atmospheric shell thickness
                t_u = t.x[i][j][k] * t_0; // in K
                e = c.x[i][j][k] * p_hydro.x[i][j][k]/ep;  // water vapour pressure in hPa
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                t_u_add = t_add_u * t_0 + t_u; // in K
                t_vir = t_u_add * (1.0 + alf * q_v_d.x[i][j][k]
                    - cloud.x[i][j][k]);
                E_Rain_add = hp * exp_func(t_vir, 17.2694, 35.86);
//                E_Ice_add = hp * exp_func(t_vir, 21.8746, 7.66);
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Rain_add = ep * E_Rain_add/(p_hydro.x[i][j][k] - E_Rain_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
//              q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
//                q_Ice_add = ep * E_Ice_add/(p_hydro.x[i][j][k] - E_Ice_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                r_humid_add = 1e2 * p_hydro.x[i][j][k]
                    /(R_Air * (1.0 + (R_W_R_A - 1.) * q_v_d.x[i][j][k]) 
                    * t_u_add); 
                q_buoy_base = r_humid_add - r_humid.x[i][j][k];
                HumidityRelative = e/E_Rain * 100.0;
                if(t_u >= t_0) L_latent = lv;
                else L_latent = ls;
                s.x[i][j][k] = (cp_l * t_u + g * height)/s_0;
// evaporation of precipitation below the cloud base
                if(i <= i_base){
//                    e_p.x[i][j][k] = C_p * alf_1   // in kg/(kg*s) // COSMO data
//                        * (q_Rain - q_v_d.x[i][j][k])
//                        * sqrt(p_ps/alf_2 * P_conv_midl.x[i][j][k]/C_p); // original COSMO formula
// ECMWF approach
                    double p_ps = p_hydro.x[i][j][k]/p_hydro.x[0][j][k];
                    if(is_land(h, 0, j, k))
                        e_p.x[i][j][k] = sig * alf_11 
                            * (RH_crit_land * q_Rain - c.x[i][j][k])
                            * pow((sqrt(p_ps)/alf_22 
                            * P_conv_midl.x[i][j][k]/sig), alf_33); // original ECMWF formula, P_conv_midl in kg/(m²*s) == mm/s
                    else 
                        e_p.x[i][j][k] = sig * alf_11 
                            * (RH_crit_water * q_Rain - c.x[i][j][k])
                            * pow((sqrt(p_ps)/alf_22 
                            * P_conv_midl.x[i][j][k]/sig), alf_33); // original ECMWF formula
                }
                if(e_p.x[i][j][k] <= 0.0) e_p.x[i][j][k] = 0.0;
// evaporation of rain precipitation within the downdraft
                if((t_u_add >= t_0)&&(q_v_d.x[i][j][k] <= q_Rain_add))
                    e_d.x[i][j][k] = a_ev 
                        * (1.0 + b_ev * pow(P_conv_midl.x[i][j][k],(1.0/6.0)))  // evaporation of rain due to water vapour diffusion, < XIII >
                        * (q_Rain_add - q_v_d.x[i][j][k]) 
                        * pow(P_conv_midl.x[i][j][k],(4.0/9.0));
                if(e_d.x[i][j][k] < 0.0) e_d.x[i][j][k] = 0.0;
// specification of entrainment and detrainment in the up/downdraft
// COSMO approach
//                E_d.x[i][j][k] = eps_d * fabs(M_d.x[i][j][k]);
//                D_d.x[i][j][k] = del_d * fabs(M_d.x[i][j][k]);

// ECMWF approach
                t_u_base = t.x[i_base][j][k] * t_0;
                E_Rain_base = hp * exp_func(t_u_base, 17.2694, 35.86);
                q_Rain_base = ep * E_Rain_base
                    /(p_hydro.x[i_base][j][k] - E_Rain_base); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                E_d.x[i][j][k] = eps_u_deep * M_d.x[i][j][k]
                    /r_humid.x[i][j][k] * (1.3 - HumidityRelative/100.0) 
                    * pow(q_Rain/q_Rain_base, 3.0);

                if(E_d.x[i][j][k] > 0.01) E_d.x[i][j][k] = 0.01;
                if(E_d.x[i][j][k] < - 0.01) E_d.x[i][j][k] = - 0.01;

// ECMWF kinetic energy (vertical velocity of the downdraft)
                if(q_buoy_base < 0.0){ // search for non-buoyant situation
                    if(M_d.x[i][j][k] == 0.0)  M_d.x[i][j][k] = 1e-6;
                    K_d[i-1] = K_d[i] + step[i]  // kinetic energy
                        * ((- E_d.x[i][j][k]/M_d.x[i][j][k]) 
                        * (1.0 + bet_d * c_d) 
                        * 2.0 * K_d[i] + g/(f_t * (1.0 + gam_d_K)) 
                        * ((t_u_add - t_u)/t_u));
                    K_d_sqrt = sqrt(K_d[i-1]/K_d[i]);
                    del_M_d = (1.6 - HumidityRelative/100.0) * K_d_sqrt;
                    D_u_2 = M_d.x[i][j][k]/(r_humid.x[i][j][k] 
                        * step[i]) * (1.0/del_M_d - 1.0);
                }else{
                    K_d[i] = 0.0;
                    K_d[i-1] = 0.0;
                    K_d_sqrt = 0.0;
                    del_M_d = 0.0;
                    D_d_2 = 0.0;
                }
                D_d.x[i][j][k] = del_u_1 * M_d.x[i][j][k]  // in kg/(m³s)
                    /r_humid.x[i][j][k] + D_d_2;

                if(D_d.x[i][j][k] > 0.01) D_d.x[i][j][k] = 0.01;
                if(D_d.x[i][j][k] < - 0.01) D_d.x[i][j][k] = - 0.01;

// mass flux within the downdraft
                M_d.x[i-1][j][k] = M_d.x[i][j][k] 
                    + step[i] * (E_d.x[i][j][k] - D_d.x[i][j][k]);  // integration from LFS downwards
                M_d_denom = M_d.x[i-1][j][k];
                M_d_num = M_d.x[i][j][k] + step[i] * D_d.x[i][j][k];
                if(M_d_denom == 0.0)  M_d_denom = 1e-6;

                if(M_d.x[i-1][j][k] > 5.0) M_d.x[i-1][j][k] = 5.0;
                if(M_d.x[i-1][j][k] < - 5.0) M_d.x[i-1][j][k] = - 5.0;

// dry entropy within the downdraft
                if(t_u >= t_0) L_latent = lv;
                else L_latent = ls;
                s_d.x[i-1][j][k] = (M_d_num * s_d.x[i][j][k] 
                    - step[i] * (E_d.x[i][j][k] * s.x[i][j][k] 
                    - D_d.x[i][j][k] * s_d.x[i][j][k]
                   + L_latent * r_humid.x[i][j][k] * e_d.x[i][j][k])/s_0)
                    /M_d_denom;
// cloud water within the downdraft
                q_v_d.x[i-1][j][k] = (M_d_num * q_v_d.x[i][j][k] 
                    + step[i] * (E_d.x[i][j][k] * c.x[i][j][k] 
                    + r_humid.x[i][j][k] * e_d.x[i][j][k]))/M_d_denom;
                if(q_v_d.x[i-1][j][k] >= q_Rain) 
                    q_v_d.x[i-1][j][k] = q_Rain;
                if(q_v_d.x[i-1][j][k] <= 0.0) q_v_d.x[i-1][j][k] = q_Rain;
// velocity components within the downdraft
                u_d.x[i-1][j][k] = u.x[i][j][k] 
                    + M_d.x[i-1][j][k]/(a_d * r_humid.x[i][j][k]); 
                if(u_d.x[i-1][j][k] > 2.5) u_d.x[i-1][j][k] = 2.5;
                if(u_d.x[i-1][j][k] < - 2.5) u_d.x[i-1][j][k] = - 2.5;
                v_d.x[i-1][j][k] = (M_d_num * v_d.x[i][j][k] 
                    + step[i] * E_d.x[i][j][k] * v.x[i][j][k])/M_d_denom;
                if(v_d.x[i+1][j][k] > .1) v_d.x[i+1][j][k] = 0.1;
                if(v_d.x[i+1][j][k] < -.1) v_d.x[i+1][j][k] = - 0.1;
                w_d.x[i-1][j][k] = (M_d_num * w_d.x[i][j][k] 
                    + step[i] * E_d.x[i][j][k] * w.x[i][j][k])/M_d_denom;
                if(w_d.x[i+1][j][k] > 10.0) w_d.x[i+1][j][k] = 10.0;
                if(w_d.x[i+1][j][k] < -10.0) w_d.x[i+1][j][k] = - 10.0;
// rain water formed by cloud convection
                P_conv_midl.x[i-1][j][k] = P_conv_midl.x[i][j][k]   // in kg/(m²*s) == mm/s
                    + step[i] * r_humid.x[i][j][k]
                    * (g_p.x[i][j][k]   // in kg/(kg*s)
                    - e_d.x[i][j][k]   // in kg/(kg*s)
                    - e_p.x[i][j][k]);  // in kg/(kg*s)
                if(P_conv_midl.x[i-1][j][k] <= 0.0) P_conv_midl.x[i-1][j][k] = 0.0;
                if(P_conv_midl.x[i-1][j][k] > 8.0/8.64e4)  P_conv_midl.x[i-1][j][k] = 8.0/8.64e4; //assumed max value
                if(P_conv_midl.x[i-1][j][k] > 0.0){
                    moistconv = true;
                    if(P_conv_midl.x[i-1][j][k] > maxValue){
                        maxValue = P_conv_midl.x[i-1][j][k];
                        P_moistconv = maxValue * 8.64e4; // in mm/d
                        i_moistconv = i;
                        j_moistconv = j;
                        k_moistconv = k;
                        height_moistconv = get_layer_height(i-1);
/*
                        cout.precision(9);
                        cout << "      max values of MoistConvectionMidL " 
                        << endl
                        << "      i_moistconv = " << i_moistconv
                        << "      j_moistconv = " << j_moistconv
                        << "   k_moistconv = " << k_moistconv
                        << "   height_moistconv = " << height_moistconv
                        << "   P_moistconv = " << P_moistconv << endl;
*/
                    }
                }
                if(is_land(h, i, j, k)||(t_u <= t_00)){
                    M_d.x[i-1][j][k] = 0.0;
                    s_d.x[i-1][j][k] = 1.0;
                    q_v_d.x[i-1][j][k] = 0.0;
                    u_d.x[i-1][j][k] = 0.0;
                    v_d.x[i-1][j][k] = 0.0;
                    w_d.x[i-1][j][k] = 0.0;
                    e_d.x[i][j][k] = 0.0;
                    e_p.x[i][j][k] = 0.0;
                }
/*
            cout.precision(9);
            if((j == 90) &&(k == 180))  cout << endl 
                << "  MoistConvectionMidL   DOWNDRAFT     -------------------------------------------" << endl
                << "  Ma = " << (int)*get_current_time() << endl
                << "  i = " << i 
                << "  j = " << j<< "  k = " << k << endl 
                << "  i_base = " << int(i_Base.y[j][k]) 
                << "  i_lfs = " << int(i_LFS.y[j][k]) << endl
                << "  s = " << s.x[i][j][k] 
                << "  s_d = " << s_d.x[i][j][k] << endl
                << "  E_d = " << E_d.x[i][j][k] 
                << "  D_d = " << D_d.x[i][j][k] << endl
                << "  e_d = " << e_d.x[i][j][k] 
                << "  e_l = " << e_l.x[i][j][k] 
                << "  e_p = " << e_p.x[i][j][k] 
                << "  c_u = " << c_u.x[i][j][k] 
                << "  g_p = " << g_p.x[i][j][k] << endl
                << "  P_convshall_mm/s = " << P_conv_shall.x[i][j][k]
                << "  P_convshall_mm/d = " << P_conv_shall.x[i][j][k] * 8.64e4 << endl
                << "  P_conv_mm/s = " << P_conv_midl.x[i][j][k]
                << "  P_conv_mm/d = " << P_conv_midl.x[i][j][k] * 8.64e4 << endl
                << "  P_rain_mm/s = " << P_rain.x[i][j][k]
                << "  P_rain_mm/d = " << P_rain.x[i][j][k] * 8.64e4 << endl
                << "  P_snow_mm/s = " << P_snow.x[i][j][k]
                << "  P_snow_mm/d = " << P_snow.x[i][j][k] * 8.64e4 << endl
                << "  E_Rain = " << E_Rain 
                << "  q_Rain_add = " << q_Rain_add
                << "  q_Rain_base = " << q_Rain_base
                << "  q_Rain = " << q_Rain << endl
                << "  q_buoy_base = " << q_buoy_base 
                << "  c = " << c.x[i][j][k] 
                << "  q_v_d = " << q_v_d.x[i][j][k] 
                << "  cloud = " << cloud.x[i][j][k] << endl
                << "  u = " << u.x[i][j][k] 
                << "  u_d = " << u_d.x[i][j][k] 
                << "  v = " << v.x[i][j][k] 
                << "  v_d = " << v_d.x[i][j][k] 
                << "  w = " << w.x[i][j][k] 
                << "  w_d = " << w_d.x[i][j][k] << endl
                << "  t_u = " << t_u - t_0 << endl
                << "  HumidityRelative = " << HumidityRelative 
                << "  r_dry = " << r_dry.x[i][j][k] 
                << "  r_humid = " << r_humid.x[i][j][k] 
                << "  r_humid_add = " << r_humid_add << endl
                << "  p_stat = " << p_stat.x[i][j][k] 
                << "  p_hydro = " << p_hydro.x[i][j][k] 
                << "  height = " << height 
//                << "  height_lfs = " << height_lfs 
                << "  step = " << step[i] << endl
                << "  M_d = " << M_d.x[i][j][k] 
                << "  M_d_LFS = " << M_d_LFS[j][k] << endl
//                << "  del_M_d = " << del_M_d
//                << "  D_d_2 = " << D_d_2
                << "  D_d = " << D_d.x[i][j][k]
                << "  K_d = " << K_d[i] << endl << endl;
*/
            } // end i convective downdraft
            e_d.x[0][j][k] = c43 * e_d.x[1][j][k] - c13 * e_d.x[2][j][k];
            e_p.x[0][j][k] = c43 * e_p.x[1][j][k] - c13 * e_p.x[2][j][k];
            E_d.x[0][j][k] = c43 * E_d.x[1][j][k] - c13 * E_d.x[2][j][k];
            D_d.x[0][j][k] = c43 * D_d.x[1][j][k] - c13 * D_d.x[2][j][k];
        }  // end j convective downdraft
    }  // end k convective downdraft
/*
*
*/
    AtomUtils::fft_gaussian_filter_3d(P_conv_midl,1);
/*
*
*/
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
// RHS for thermodynamic forcing due to moist convection
            for(int i = 0; i < im-1; i++){
                step[i] = get_layer_height(i+1) - get_layer_height(i);  // local atmospheric shell thickness
                if((t.x[i][j][k] * t_0) >= t_0) L_latent = lv;
                else L_latent = ls;
                MC_t.x[i][j][k] = 
                    - ((M_u.x[i][j][k] * (s_u.x[i][j][k] - s.x[i][j][k]) 
                    + M_d.x[i][j][k] * (s_d.x[i][j][k] - s.x[i][j][k])) 
                    - (M_u.x[i][j][k] * (s_u.x[i][j][k] - s.x[i][j][k]) 
                    + M_d.x[i][j][k] * (s_d.x[i][j][k] - s.x[i][j][k]))) 
                    /(r_humid.x[i][j][k] * step[i])/cp_l * s_0
                    + L_latent/cp_l * (c_u.x[i][j][k] - e_d.x[i][j][k] 
                    - e_l.x[i][j][k] - e_p.x[i][j][k]);  // in (K)(1/s)
                MC_q.x[i][j][k] = 
                    - ((M_u.x[i][j][k] * (q_v_u.x[i][j][k] - c.x[i][j][k]) 
                    + M_d.x[i][j][k] * (q_v_d.x[i][j][k] - c.x[i][j][k])) 
                    - (M_u.x[i][j][k] * (q_v_u.x[i][j][k] - c.x[i][j][k]) 
                    + M_d.x[i][j][k] * (q_v_d.x[i][j][k] - c.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k])
                    - (c_u.x[i][j][k] - e_d.x[i][j][k] - e_l.x[i][j][k] 
                    - e_p.x[i][j][k]);  // in (kg/kg)(1/s)
                MC_v.x[i][j][k] = 
                    - ((M_u.x[i][j][k] * (v_u.x[i][j][k] - v.x[i][j][k]) 
                    + M_d.x[i][j][k] * (v_d.x[i][j][k] - v.x[i][j][k])) 
                    - (M_u.x[i][j][k] * (v_u.x[i][j][k] - v.x[i][j][k]) 
                    + M_d.x[i][j][k] * (v_d.x[i][j][k] - v.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k]) * u_0;  // in (m/s)(1/s)
                MC_w.x[i][j][k] = 
                    - ((M_u.x[i][j][k] * (w_u.x[i][j][k] - w.x[i][j][k]) 
                    + M_d.x[i][j][k] * (w_d.x[i][j][k] - w.x[i][j][k])) 
                    - (M_u.x[i][j][k] * (w_u.x[i][j][k] - w.x[i][j][k]) 
                    + M_d.x[i][j][k] * (w_d.x[i][j][k] - w.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k]) * u_0;  // in (m/s)(1/s)
            }  // MC-loop
        } // j-loop
    } // k-loop
/*
// boundaries of various variables
    for(int k = 1; k < km-1; k++){
        for(int j = 2; j < jm-2; j++){
            M_u.x[0][j][k] = c43 * M_u.x[1][j][k] -
                c13 * M_u.x[2][j][k];
            M_u.x[im-1][j][k] = c43 * M_u.x[im-2][j][k] -
                c13 * M_u.x[im-3][j][k];
            M_d.x[0][j][k] = c43 * M_d.x[1][j][k] -
                c13 * M_d.x[2][j][k];
            M_d.x[im-1][j][k] = c43 * M_d.x[im-2][j][k] -
                c13 * M_d.x[im-3][j][k];
            s.x[0][j][k] = c43 * s.x[1][j][k] -
                c13 * s.x[2][j][k];
            s.x[im-1][j][k] = c43 * s.x[im-2][j][k] -
                c13 * s.x[im-3][j][k];
            s_u.x[0][j][k] = c43 * s_u.x[1][j][k] -
                c13 * s_u.x[2][j][k];
            s_u.x[im-1][j][k] = c43 * s_u.x[im-2][j][k] -
                c13 * s_u.x[im-3][j][k];
            s_d.x[0][j][k] = c43 * s_d.x[1][j][k] -
                c13 * s_d.x[2][j][k];
            s_d.x[im-1][j][k] = c43 * s_d.x[im-2][j][k] -
                c13 * s_d.x[im-3][j][k];
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int i = 0; i < im; i++){
            M_u.x[i][0][k] = c43 * M_u.x[i][1][k] -
                 c13 * M_u.x[i][2][k];
            M_u.x[i][jm-1][k] = c43 * M_u.x[i][jm-2][k] -
                 c13 * M_u.x[i][jm-3][k];
            M_d.x[i][0][k] = c43 * M_d.x[i][1][k] -
                 c13 * M_d.x[i][2][k];
            M_d.x[i][jm-1][k] = c43 * M_d.x[i][jm-2][k] -
                 c13 * M_d.x[i][jm-3][k];
            s.x[i][0][k] = c43 * s.x[i][1][k] -
                 c13 * s.x[i][2][k];
            s.x[i][jm-1][k] = c43 * s.x[i][jm-2][k] -
                 c13 * s.x[i][jm-3][k];
            s_u.x[i][0][k] = c43 * s_u.x[i][1][k] -
                 c13 * s_u.x[i][2][k];
            s_u.x[i][jm-1][k] = c43 * s_u.x[i][jm-2][k] -
                 c13 * s_u.x[i][jm-3][k];
            s_d.x[i][0][k] = c43 * s_d.x[i][1][k] -
                 c13 * s_d.x[i][2][k];
            s_d.x[i][jm-1][k] = c43 * s_d.x[i][jm-2][k] -
                 c13 * s_d.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 1; j < jm-1; j++){
            M_u.x[i][j][0] = c43 * M_u.x[i][j][1] -
                c13 * M_u.x[i][j][2];
            M_u.x[i][j][km-1] = c43 * M_u.x[i][j][km-2] -
                c13 * M_u.x[i][j][km-3];
            M_u.x[i][j][0] = M_u.x[i][j][km-1] =
               (M_u.x[i][j][0] + M_u.x[i][j][km-1])/2.;
            M_d.x[i][j][0] = c43 * M_d.x[i][j][1] -
                c13 * M_d.x[i][j][2];
            M_d.x[i][j][km-1] = c43 * M_d.x[i][j][km-2] -
                c13 * M_d.x[i][j][km-3];
            M_d.x[i][j][0] = M_d.x[i][j][km-1] =
               (M_d.x[i][j][0] + M_d.x[i][j][km-1])/2.;
            s.x[i][j][0] = c43 * s.x[i][j][1] -
                c13 * s.x[i][j][2];
            s.x[i][j][km-1] = c43 * s.x[i][j][km-2] -
                c13 * s.x[i][j][km-3];
            s.x[i][j][0] = s.x[i][j][km-1] =
               (s.x[i][j][0] + s.x[i][j][km-1])/2.;
            s_u.x[i][j][0] = c43 * s_u.x[i][j][1] -
                c13 * s_u.x[i][j][2];
            s_u.x[i][j][km-1] = c43 * s_u.x[i][j][km-2] -
                c13 * s_u.x[i][j][km-3];
            s_u.x[i][j][0] = s_u.x[i][j][km-1] =
               (s_u.x[i][j][0] + s_u.x[i][j][km-1])/2.;
            s_d.x[i][j][0] = c43 * s_d.x[i][j][1] -
                c13 * s_d.x[i][j][2];
            s_d.x[i][j][km-1] = c43 * s_d.x[i][j][km-2] -
                c13 * s_d.x[i][j][km-3];
            s_d.x[i][j][0] = s_d.x[i][j][km-1] =
               (s_d.x[i][j][0] + s_d.x[i][j][km-1])/2.;
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if(is_land(h, i, j, k)){
                    P_conv_midl.x[i][j][k] = P_conv_midl.x[i_mount][j][k];
                    M_u.x[i][j][k] = M_u.x[i_mount][j][k];
                    M_d.x[i][j][k] = M_d.x[i_mount][j][k];
                    MC_t.x[i][j][k] = MC_t.x[i_mount][j][k];
                    MC_q.x[i][j][k] = MC_q.x[i_mount][j][k];
                    MC_v.x[i][j][k] = MC_v.x[i_mount][j][k];
                    MC_w.x[i][j][k] = MC_w.x[i_mount][j][k];
                    s.x[i][j][k] = s.x[i_mount][j][k];
                    s_u.x[i][j][k] = s_u.x[i_mount][j][k];
                    s_d.x[i][j][k] = s_d.x[i_mount][j][k];
                    q_v_u.x[i][j][k] = q_v_u.x[i_mount][j][k];
                    q_v_d.x[i][j][k] = q_v_d.x[i_mount][j][k];
                    q_c_u.x[i][j][k] = q_c_u.x[i_mount][j][k];
                    E_u.x[i][j][k] = E_u.x[i_mount][j][k];
                    E_d.x[i][j][k] = E_d.x[i_mount][j][k];
                    D_u.x[i][j][k] = D_u.x[i_mount][j][k];
                    D_d.x[i][j][k] = D_d.x[i_mount][j][k];
                    c_u.x[i][j][k] = c_u.x[i_mount][j][k];
                    g_p.x[i][j][k] = g_p.x[i_mount][j][k];
                    e_l.x[i][j][k] = e_l.x[i_mount][j][k];
                    e_p.x[i][j][k] = e_p.x[i_mount][j][k];
                    e_d.x[i][j][k] = e_d.x[i_mount][j][k];
                    u_u.x[i][j][k] = u_u.x[i_mount][j][k];
                    u_d.x[i][j][k] = u_d.x[i_mount][j][k];
                    v_u.x[i][j][k] = v_u.x[i_mount][j][k];
                    v_d.x[i][j][k] = v_d.x[i_mount][j][k];
                    w_u.x[i][j][k] = w_u.x[i_mount][j][k];
                    w_d.x[i][j][k] = w_d.x[i_mount][j][k];
                }
            }
        }
    }
*/
    if(moistconv == false)  
        cout << "      no saturation of water vapour in MoistConvectionMidL found" 
        << endl;
    else
        cout << "      saturation of water vapour in MoistConvectionMidL found" 
        << endl
        << "      i_moistconv = " << i_moistconv
        << "      j_moistconv = " << j_moistconv
        << "   k_moistconv = " << k_moistconv
        << "   height_moistconv[m] = " << height_moistconv
        << "   P_moistconv[mm/d] = " << P_moistconv << endl;
    cout << "      MoistConvectionMidL ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::MoistConvectionShall(){
    cout << endl << "      MoistConvectionShall" << endl;
// collection of coefficients for phase transformation
    int i_base = 0;
    int ii_base_beg = 1;
    int i_lfs = 0;
    double a_ev = 1.0e-3;  // m2/kg
    double b_ev = 5.9;  // m2*s/kg
    double R_W_R_A = R_WaterVapour/R_Air;
    double eps = 1.0/R_W_R_A;
    double alf = (1.0 - eps)/eps;
    double t_00 = 236.15; // = -37°C
    bool moistconv = false;
    double P_moistconvshall = 0.0;
    int i_moistconv = 0;
    int j_moistconv = 0;
    int k_moistconv = 0;
    double height_moistconv = 0.0;
    double maxValue = 0.0;
/*
// constants used for the COSMO parameterization
    double bet = 42.0; // in K given by COSMO
    double b_u = 0.3;
    double alf_1 = 5.0e-4; // in 1/s
    double alf_2 = 0.011; // in 1/s
    double C_p = 0.05; // in ./.
    double bet_p = 2.0e-3;  // in 1/s, from the COSMO paper
    double eps_u = 1.0e-4;  // in 1/m
    double eps_d = 2.0e-4;
    double del_u = eps_u;
    double del_d = eps_d;
    double delta_i_c = 0.0;
    double K_p = 0.0;  // in 1/s
*/
    double gam_d = - 0.3;  // from original paper by Tiedtke, gam_d = - 0.2, COSMO approach
    double a_u = 1.0e0;  // area fraction in updraft to define the vertical velocity
    double a_d = 1.0e0;  // area fraction in downdraft to define the vertical velocit
    double HumidityRelative = 0.0;
    double L_latent = 0.0;
    double q_buoy_base = 0.0;
    double E_Rain = 0.0;
    double q_Rain = 0.0;
    double e = 0.0;
//    double E_Ice = 0.0;
    double E_Rain_add = 0.0;
    double E_Rain_base = 0.0;
    double q_Rain_add = 0.0;
    double q_Rain_base = 0.0;
//    double E_Ice_add = 0.0;
//    double q_Ice = 0.0;
//    double q_Ice_add = 0.0;
//    double r_humid_add = 0.0;
//    double r_dry_add = 0.0;
//    double p_stat_add = 0.0;
    double r_humid_add = 0.0;
//    double r_dry_add = 0.0;
    double t_u = 0.0;
    double t_u_add = 0.0;
    double t_u_base = 0.0;
    double height = 0.0;
//    double height_base = 0.0;
//    double height_lfs = 0.0;
    double rm = 0.0;
    double exp_rm = 0.0;
    double M_u_denom = 0.0;
    double M_d_denom = 0.0;
    double M_u_num = 0.0;
    double M_d_num = 0.0;
    double t_vir = 0.0;
    double t_vir_env = 0.0;
//    double t_add_u = 0.0; //  in K == + 0.0 °C temperature increase for a lifting air parcel, given by ECMWF = European Center for Medium-Range Weather Forecast
    double t_add_u = 0.00732; //  in K == + 2.0 °C temperature increase for a lifting air parcel, given by ECMWF = European Center for Medium-Range Weather Forecast
//    double q_v_u_add = 0.0; // == 0.000 kg/kg increase for a lifting air parcel, testing
    double q_v_u_add = 5.0e-3; // == 0.005 kg/kg increase for a lifting air parcel, given by ECMWF
    double p_stat_check = 0.0;
    std::vector<double> step(im, 0);
    std::vector<std::vector<double> > M_u_Base(jm, std::vector<double> (km, 0));
    std::vector<std::vector<double> > M_d_LFS(jm, std::vector<double> (km, 0));
    std::vector<std::vector<double> > u_cloud_base(jm, std::vector<double> (km, 0));
    rad.Coordinates(im, r0, dr);
/*
    double dcdr = 0.0, dcdthe = 0.0, dcdphi = 0.0; // COSMO data
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            int i_mount = i_topography[j][k];
            sinthe = sin(the.z[j]);
            for(int i = im-1; i >= i_mount ; i--){
                rm = rad.z[i] * get_layer_height(i);
                rmsinthe = rm * sinthe;
                exp_rm = 1.0/(rm + 1.0);
                dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k])/(2. * dr) * exp_rm;
                if(i < im-2){
                    if((is_land(h, i, j, k))&&((is_air(h, i+1, j, k))
                        &&(is_air(h, i+2, j, k)))){
                        dcdr = (- 3. * c.x[i][j][k] + 4. * c.x[i+1][j][k] 
                            - c.x[i+2][j][k])/(2. * dr ) * exp_rm;
                    }
                }
                if(c.x[i][j][k] == 0.0) c.x[i][j][k] = c_tropopause;
                dcdr = 0.0; 
                dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k])
                    /(2.0 * dthe * rm);
                dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1])
                    /(2.0 * dphi * rmsinthe);

                if((j >= 2)&&(j < jm-3)){
                    if((is_land(h, i, j, k))&&((is_air(h, i, j+1, k))
                        &&(is_air(h, i, j+2, k))))
                        dcdthe = (- 3. * c.x[i][j][k] + 4. * c.x[i][j+1][k] 
                            - c.x[i][j+2][k])/(2. * dthe * rm);
                    if((is_land(h, i, j, k))&&(is_air(h, i, j-1, k))
                        &&(is_air(h, i, j-2, k)))
                        dcdthe = (- 3. * c.x[i][j][k] + 4. * c.x[i][j-1][k] 
                            - c.x[i][j-2][k])/(2. * dthe * rm);
                }
                if((k >= 2)&&(k < km-3)){
                    if((is_land(h, i, j, k))&&(is_air(h, i, j, k+1))
                        &&(is_air(h, i, j, k+2)))
                        dcdphi = (- 3. * c.x[i][j][k] + 4. * c.x[i][j][k+1] 
                            - c.x[i][j][k+2])/(2. * dphi * rmsinthe);
                    if((is_land(h, i, j, k))&&(is_air(h, i, j, k-1))
                        &&(is_air(h, i, j, k-2)))
                        dcdphi = (- 3. * c.x[i][j][k] + 4. * c.x[i][j][k-1] 
                            - c.x[i][j][k-2])/(2. * dphi * rmsinthe);
                }
                if((j >= 2)&&(j < jm-3)){
                    if((is_land(h, i, j, k))&&((is_air(h, i, j+1, k))
                        &&(is_air(h, i, j+2, k))))
                        dcdthe = (c.x[i][j+1][k] - c.x[i][j][k])/(dthe * rm);
                    if((is_land(h, i, j, k))&&(is_air(h, i, j-1, k))
                        &&(is_air(h, i, j-2, k)))
                        dcdthe = (c.x[i][j][k] - c.x[i][j-1][k])/(dthe * rm);
                }
                if((k >= 2)&&(k < km-3)){
                    if((is_land(h, i, j, k))&&(is_air(h, i, j, k+1))
                        &&(is_air(h, i, j, k+2)))
                        dcdphi = (c.x[i][j][k+1] - c.x[i][j][k])/(dphi * rmsinthe);
                    if((is_land(h, i, j, k))&&(is_air(h, i, j, k-1))
                        &&(is_air(h, i, j, k-2)))
                        dcdphi = (c.x[i][j][k] - c.x[i][j][k-1])/(dphi * rmsinthe);
                }
// specification of entrainment in the updraft
                if(c.x[i][j][k] == 0.0) c.x[i][j][k] = 1e+6;
                E_u.x[i][j][k] = - r_humid.x[i][j][k]/c.x[i][j][k] 
                    * (u.x[i][j][k] * dcdr 
                    + v.x[i][j][k] * dcdthe 
                    + w.x[i][j][k] * dcdphi) * u_0; // in kg/(m³s)
                if(E_u.x[i][j][k] >= 0.2) E_u.x[i][j][k] = 0.2;
                if(E_u.x[i][j][k] <= - 0.2) E_u.x[i][j][k] = - 0.2;
            }
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            E_u.x[0][j][k] = c43 * E_u.x[1][j][k] -
                c13 * E_u.x[2][j][k];
            E_u.x[im-1][j][k] = c43 * E_u.x[im-2][j][k] -
                c13 * E_u.x[im-3][j][k];
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            E_u.x[i][0][k] = c43 * E_u.x[i][1][k] -
                 c13 * E_u.x[i][2][k];
            E_u.x[i][jm-1][k] = c43 * E_u.x[i][jm-2][k] -
                 c13 * E_u.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            E_u.x[i][j][0] = c43 * E_u.x[i][j][1] -
                c13 * E_u.x[i][j][2];
            E_u.x[i][j][km-1] = c43 * E_u.x[i][j][km-2] -
                c13 * E_u.x[i][j][km-3];
            E_u.x[i][j][0] = E_u.x[i][j][km-1] =
               (E_u.x[i][j][0] + E_u.x[i][j][km-1])/2.;
        }
    }
    AtomUtils::fft_gaussian_filter_3d(E_u,1);
*/
// parameterisation of moist convection
    CAPE = std::vector<double>(im, 0.0);
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
// boundary conditions on the base of clouds (level of free convection)
            int i_mount = i_topography[j][k];
            if(i_mount >= 0) 
                ii_base_beg = i_mount;
            else ii_base_beg = 0;
            for(int ii_base = ii_base_beg; ii_base < im-1; ii_base++){
                t_u = t.x[ii_base][j][k] * t_0; // in K
                t_u_add = t_add_u * t_0 + t_u; // in K
                q_v_u.x[ii_base][j][k] = c.x[ii_base][j][k] + q_v_u_add;
                q_c_u.x[ii_base][j][k] = cloud.x[ii_base][j][k];
                r_humid_add = 1e2 * p_hydro.x[ii_base][j][k]
                    /(R_Air * (1. + (R_W_R_A - 1.) * q_v_u.x[ii_base][j][k] 
                    - q_c_u.x[ii_base][j][k]) * t_u_add); 
                t_vir = t_u_add * (1.0 + alf * q_v_u.x[ii_base][j][k] 
                    - q_c_u.x[ii_base][j][k]);
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Rain_add = hp * exp_func(t_vir, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_hydro.x[ii_base][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Rain_add = ep * E_Rain_add/(p_hydro.x[ii_base][j][k] - E_Rain_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                e = q_v_u.x[ii_base][j][k] * p_hydro.x[ii_base][j][k]/ep;  // water vapour pressure in hPa
                HumidityRelative = e/E_Rain_add * 100.0;
                q_buoy_base = r_humid_add - r_humid.x[ii_base][j][k];
/*
                cout.precision(6);
                if((j == 90) &&(k == 180))  cout << endl 
                    << "  MoistConvectionMidL   Cloud Base" << endl 
                    << "  Ma = " << Ma << endl
                    << "  i_base UPDRAFT" << endl 
                    << "  j = " << j << "  k = " << k << endl 
                    << "  ii_base = " << ii_base 
                    << "  r_dry.x[ii_base] = " << r_dry.x[ii_base][j][k] 
                    << "  r_humid.x[ii_base] = " << r_humid.x[ii_base][j][k] 
                    << "  r_humid_add = " << r_humid_add << endl
                    << "  E_Rain = " << E_Rain 
                    << "  E_Rain_add = " << E_Rain_add 
                    << "  q_Rain = " << q_Rain 
                    << "  q_Rain_add = " << q_Rain_add << endl
                    << "  HumidityRelative = " << HumidityRelative
                    << "  q_buoy_base = " << q_buoy_base  << endl << endl << endl;
*/
                if(q_buoy_base >= 0.0){ // search for cloud base
                    i_Base.y[j][k] = (double)ii_base;
                    i_base = ii_base;
                    height = get_layer_height(i_base);
                    if(is_land(h, i_base, j, k)){
                       int i_mount = i_topography[j][k];
                       u.x[i_base][j][k] = 0.0;
                       v.x[i_base][j][k] = 0.0;
                       w.x[i_base][j][k] = 0.0;
                       c.x[i_base][j][k] = c.x[i_mount][j][k];
                       cloud.x[i_base][j][k] = cloud.x[i_mount][j][k];
                       ice.x[i_base][j][k] = ice.x[i_mount][j][k];
                    }
                    t_u = t.x[i_base][j][k] * t_0; // in K
                    t_u_add = t_add_u * t_0 + t_u; // in K
                    r_humid_add = 1e2 * p_hydro.x[i_base][j][k]
                        /(R_Air * (1. + (R_W_R_A - 1.) * q_v_u.x[i_base][j][k] 
                        - q_c_u.x[i_base][j][k]) * t_u_add); 
                    if(iter_cnt == -1)  u_cloud_base[j][k] = u.x[i_base][j][k];
                    s.x[i_base][j][k] = (cp_l * t_u + g * height
                        + lv * c.x[i_base][j][k])/s_0;
                    q_v_u.x[i_base][j][k] = c.x[i_base][j][k] + q_v_u_add;
                    q_c_u.x[i_base][j][k] = cloud.x[i_base][j][k];
                    u_u.x[i_base][j][k] = u.x[i_base][j][k];
                    v_u.x[i_base][j][k] = v.x[i_base][j][k];
                    w_u.x[i_base][j][k] = w.x[i_base][j][k];
                    s_u.x[i_base][j][k] = (cp_l * t_vir + g * height 
                        + lv * q_v_u.x[i_base][j][k])/s_0;
                    g_p.x[i_base][j][k] = 0.0;
                    c_u.x[i_base][j][k] = 0.0;
                    e_l.x[i_base][j][k] = 0.0;
                    e_p.x[i_base][j][k] = 0.0;
                    M_u.x[i_base][j][k] = r_humid.x[i_base][j][k]
                        * u.x[i_base][j][k] * u_0 * s_u.x[i_base][j][k]
                        /(s_u.x[i_base][j][k] - s.x[i_base][j][k]);  // in kg/(m²s) [== mm/s] updraft at cloud base 
                    M_d_LFS[j][k] = gam_d * M_u.x[i_base][j][k];  // in kg/(m²s) [== mm/s] downdraft at cloud top 
                    if(is_land(h, i_base, j, k))  
                        M_u.x[i_base][j][k] = M_d_LFS[j][k] = 0.0;
//                    u_u.x[i_base][j][k] = u.x[i_base][j][k] 
//                        + M_u.x[i_base][j][k]/(a_d * r_humid_add); 
                    rm = rad.z[i_base];
                    exp_rm = 1./(rm + 1.0);
                    t_vir = t_u_add * (1.0 + alf * q_v_u.x[i_base][j][k] 
                        - q_c_u.x[i_base][j][k]);
                    t_vir_env = t_u * (1.0 + alf * c.x[i_base][j][k] 
                        - cloud.x[i_base][j][k]);
                    step[i_base] = get_layer_height(i_base+1) 
                        - get_layer_height(i_base);  // in m local atmospheric shell thickness
                    CAPE[i_base] = g * step[i_base] // convective available potential energy
                        /exp_rm * (t_vir - t_vir_env)/t_vir_env;
/*
                    cout.precision(6);
                    if((j == 90) &&(k == 180))  cout << endl 
                        << "  MoistConvectionMidL   Cloud Base    +++++++++++++++++++++++++++++++++++++++++++++++++++" << endl 
                        << "  Ma = " << Ma << endl
                        << "  i_base UPDRAFT" << endl 
                        << "  j = " << j << "  k = " << k << endl 
                        << "  height_base = " << height << endl
                        << "  ii_base = " << ii_base 
                        << "  i_base = " << (int)i_Base.y[j][k]
                        << "  i_lfs = " << (int)i_LFS.y[j][k] << endl
                        << "  i_tropo = " << i_topography[j][k] << endl
                        << "  s = " << s.x[i_base][j][k] 
                        << "  s_u = " << s_u.x[i_base][j][k] << endl
                        << "  E_Rain = " << E_Rain 
                        << "  E_Rain_add = " << E_Rain_add 
                        << "  q_Rain = " << q_Rain 
                        << "  q_Rain_add = " << q_Rain_add << endl
                        << "  q_buoy_base = " << q_buoy_base 
                        << "  c = " << c.x[i_base][j][k] 
                        << "  q_v_u = " << q_v_u.x[i_base][j][k] 
                        << "  cloud = " << cloud.x[i_base][j][k] 
                        << "  q_c_u = " << q_c_u.x[i_base][j][k] << endl
                        << "  u = " << u.x[i_base][j][k] 
                        << "  u_u = " << u_u.x[i_base][j][k] 
                        << "  v = " << v.x[i_base][j][k] 
                        << "  v_u = " << v_u.x[i_base][j][k] 
                        << "  w = " << w.x[i_base][j][k] 
                        << "  w_u = " << w_u.x[i_base][j][k] << endl
                        << "  t_u = " << t_u - t_0 
                        << "  t_u_add = " << t_add_u * t_u - t_0 << endl
                        << "  HumidityRelative = " << HumidityRelative 
                        << "  r_dry = " << r_dry.x[i_base][j][k] 
                        << "  r_dry_add = " << r_dry_add 
                        << "  r_humid = " << r_humid.x[i_base][j][k] 
                        << "  r_humid_add = " << r_humid_add << endl
                        << "  p_stat = " << p_stat.x[i_base][j][k] 
                        << "  p_hydro = " << p_hydro.x[i_base][j][k] 
                        << "  M_u = " << M_u.x[i_base][j][k] 
                        << "  M_d_LFS = " << M_d_LFS[j][k] << endl
                        << "  rad = " << rad.z[i_base]
                        << "  step_non_dim = " << g * (rad.z[i_base] - rad.z[i_base-1])/exp_rm
                        << "  temp_non_dim = " << (t_u_add - t_u)/t_u
                        << "  CAPE = " << CAPE[i_base] << endl << endl;
*/
                    break;
                } // i_base-loop
            } // ii_base-loop
// level of free sinking
/*
*
*/
// boundary conditions on top of the cloud
            for(int ii_lfs = im-2; ii_lfs > 0; ii_lfs--){
                t_u = t.x[ii_lfs][j][k] * t_0; // in K
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_hydro.x[ii_lfs][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                t_u_add = t_add_u * t_0 + t_u; // in K
                r_humid_add = 1e2 * p_hydro.x[ii_lfs][j][k]/(R_Air 
                    * (1. + (R_W_R_A - 1.) * c.x[ii_lfs][j][k] 
                    - cloud.x[ii_lfs][j][k] - ice.x[ii_lfs][j][k]) 
                    * t_u_add); 
                q_v_d.x[ii_lfs][j][k] = c.x[ii_lfs][j][k];
                p_stat_check = p_hydro.x[i_base][j][k] 
                    - p_hydro.x[ii_lfs][j][k];


//                if((u_u.x[ii_lfs][j][k] <= 0.0)&&(t_u >= t_0)){ // search for cloud top
                if((u_u.x[ii_lfs][j][k] <= 0.0)
                    &&(t_u >= t_0)
                    &&(p_stat_check <= 200.0)){ // search for cloud top
                    i_LFS.y[j][k] = (double)ii_lfs;
                    i_lfs = ii_lfs;
                    height = get_layer_height(i_lfs); // in m
                    if((t.x[i_lfs][j][k] * t_0) >= t_0) 
                        L_latent = lv;
                    else L_latent = ls;
                    s.x[i_lfs][j][k] = (cp_l * t_u + g * height
                        + lv * c.x[i_lfs][j][k])/s_0;
                    q_v_d.x[i_lfs][j][k] = c.x[i_lfs][j][k];
                    M_d.x[i_lfs][j][k] = M_d_LFS[j][k];  // in kg/(m²s) [== mm/s] downdraft at cloud top
//                    u_d.x[i_lfs][j][k] = u.x[i_lfs][j][k] 
//                        + M_d.x[i_lfs][j][k]/(a_d * r_humid.x[i_lfs][j][k]); 
                    u_d.x[i_lfs][j][k] = u.x[i_lfs][j][k]; 
                    v_d.x[i_lfs][j][k] = v.x[i_lfs][j][k];
                    w_d.x[i_lfs][j][k] = w.x[i_lfs][j][k];
                    s_d.x[i_lfs][j][k] = s.x[i_lfs][j][k];
                    e_p.x[i_lfs][j][k] = 0.0;
                    e_d.x[i_lfs][j][k] = 0.0;
/*
                    cout.precision(6);
                    if((j == 90) &&(k == 180))  cout << endl 
                        << "  MoistConvectionMidL   Cloud Top     ---------------------------------------------"  << endl 
                        << "  Ma = " << Ma << endl
                        << "  i_lfs DOWNDRAFT" << endl
                        << "  j = " << j<< "  k = " << k << endl 
                        << "  height_lfs = " << height << endl
                        << "  ii_lfs = " << ii_lfs 
                        << "  i_base = " << (int)i_Base.y[j][k] 
                        << "  i_lfs = " << (int)i_LFS.y[j][k] << endl
                        << "  s = " << s.x[i_lfs][j][k] 
                        << "  s_d = " << s_d.x[i_lfs][j][k] << endl
                        << "  E_Rain = " << E_Rain 
                        << "  q_Rain = " << q_Rain << endl
                        << "  q_buoy_base = " << q_buoy_base 
                        << "  c = " << c.x[i_lfs][j][k] 
                        << "  q_v_d = " << q_v_d.x[i_lfs][j][k] 
                        << "  cloud = " << cloud.x[i_lfs][j][k] << endl
                        << "  u = " << u.x[i_lfs][j][k] 
                        << "  u_d = " << u_d.x[i_lfs][j][k] 
                        << "  v = " << v.x[i_lfs][j][k] 
                        << "  v_d = " << v_d.x[i_lfs][j][k] 
                        << "  w = " << w.x[i_lfs][j][k] 
                        << "  w_d = " << w_d.x[i_lfs][j][k] << endl
                        << "  t = " << t.x[i_lfs][j][k] * t_0 - t_0 << endl
                        << "  HumidityRelative = " << HumidityRelative 
                        << "  r_dry = " << r_dry.x[i_lfs][j][k] 
                        << "  r_humid = " << r_humid.x[i_lfs][j][k] << endl
                        << "  p_stat = " << p_stat.x[i_lfs][j][k] 
                        << "  p_hydro = " << p_hydro.x[i_lfs][j][k] 
                        << "  M_d = " << M_d.x[i_lfs][j][k] 
                        << "  M_d_LFS = " << M_d_LFS[j][k] << endl << endl;
*/
                    break;
                }
            } // ii_lfs-loop
/*
*
*/
// zero moist values below i_base
            for(int i = 0; i < i_base; i++){
                t_u = t.x[i][j][k] * t_0; // in K
                if(t_u >= t_0) L_latent = lv;
                else L_latent = ls;
                height = get_layer_height(i);
                s.x[i][j][k] = (cp_l * t_u + g * height 
                    + lv * c.x[i][j][k])/s_0;
                s_u.x[i][j][k] = 0.0;
                s_d.x[i][j][k] = 0.0;
                q_v_u.x[i][j][k] = 0.0;
                q_c_u.x[i][j][k] = 0.0;
                u_u.x[i][j][k] = 0.0;
                v_u.x[i][j][k] = 0.0;
                w_u.x[i][j][k] = 0.0;
                g_p.x[i][j][k] = 0.0;
                c_u.x[i][j][k] = 0.0;
                e_l.x[i][j][k] = 0.0;
            } // zero moist values below i_base
// zero moist values above i_lfs
            for(int i = i_lfs+1; i < im; i++){
                t_u = t.x[i][j][k] * t_0; // in K
                height = get_layer_height(i);
                s.x[i][j][k] = (cp_l * t_u + g * height 
                    + lv * c.x[i][j][k])/s_0;
                s_u.x[i][j][k] = 0.0;
                s_d.x[i][j][k] = 0.0;
                q_v_d.x[i][j][k] = 0.0;
                u_d.x[i][j][k] = 0.0;
                v_d.x[i][j][k] = 0.0;
                w_d.x[i][j][k] = 0.0;
            } // zero moist values above i_lfs
        }  // end j for cloud base and cloud top
    }  // end k for cloud base and cloud top
/*
*
*/
// ECMWF data (European Center of Mid-Range Weather Forecast)
    double gam_d_K = 0.5;
    double c_d = 0.506;
    double bet_d = 1.875;
    double f_t = 2.0;
    double del_u_1 = 0.75e-4; // 1/m
    double K_u_sqrt = 0.0;
    double K_d_sqrt = 0.0;
    double D_u_2 = 0.0;
    double del_M_u = 0.0;
    double D_d_2 = 0.0;
    double del_M_d = 0.0;
    double eps_u_deep = 1.75e-3; // 1/m
    double eps_u_shallow = 2.0 * eps_u_deep; // 1/m
    double alf_11 = 5.44e-4; // s
    double alf_22 = 5.09e-3;
    double alf_33 = 0.577;
    double RH_crit_water = 0.9; // over water
    double RH_crit_land = 0.7; // over land
    double sig = 0.05;
//
// CONVECTIVE UPDRAFT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            i_base = (int)i_Base.y[j][k];
            i_lfs = (int)i_LFS.y[j][k];
            K_u = std::vector<double>(im, 0.0);
            K_u[i_base] = 0.5 * u_u.x[i_base][j][k] * u_u.x[i_base][j][k];
            for(int i = i_base; i < im-1; i++){
                p_stat_check = p_hydro.x[i_base][j][k] 
                    - p_hydro.x[i_lfs][j][k];
//            if((u_u.x[i][j][k] >= 0.0)
//                &&(p_stat_check <= 200.0) // hPa

                if(
                    (p_stat_check <= 200.0) // hPa
                    &&(t.x[i][j][k] * t_0 > t_0)){

                height = get_layer_height(i);
                step[i] = get_layer_height(i+1) - get_layer_height(i);  // in m local atmospheric shell thickness
                t_u = t.x[i][j][k] * t_0; // in K
                t_u_add = t_add_u * t_0 + t_u; // in K
                t_vir = t_u_add * (1.0 + alf * q_v_u.x[i][j][k] 
                    - q_c_u.x[i][j][k]);
                t_vir_env = t_u * (1.0 + alf * c.x[i][j][k] 
                    - cloud.x[i][j][k]);
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
//                E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                E_Rain_add = hp * exp_func(t_vir, 17.2694, 35.86);
//                E_Ice_add = hp * exp_func(t_u, 21.8746, 7.66);
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Rain_add = ep * E_Rain_add/(p_hydro.x[i][j][k] - E_Rain_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
//                q_Ice = ep * E_Ice/(p_hydro.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
//                q_Ice_add = ep * E_Ice_add/(p_hydro.x[i][j][k] - E_Ice_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                t_u_base = t.x[i_base][j][k] * t_0;
                E_Rain_base = hp * exp_func(t_u_base, 17.2694, 35.86);
                q_Rain_base = ep * E_Rain_base
                    /(p_hydro.x[i_base][j][k] - E_Rain_base); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                e = c.x[i][j][k] * p_hydro.x[i][j][k]/ep;  // water vapour pressure in hPa
                HumidityRelative = e/E_Rain * 100.0;
                r_humid_add = 1e2 * p_hydro.x[i][j][k]
                    /(R_Air * (1.0 + (R_W_R_A - 1.) * q_v_u.x[i][j][k] 
                    - q_c_u.x[i][j][k]) * t_u_add); 
                q_buoy_base = r_humid_add - r_humid.x[i][j][k];
                if(step[i] >= t_0) L_latent = lv;
                else L_latent = ls;
                s.x[i][j][k] = (cp_l * t_u + g * height 
                    + lv * c.x[i][j][k])/s_0;
/*
// COSMO data
                eps_u = 1.e-4;
                del_u = 1.e-4;
                if(t_u_add < t_00){
                    eps_u = 0.0;
                    del_u = 0.0;
                }
*/
// CAPE (Convective Available Potential Energy) == difference between the moist adiabatic and the environmental temperature
                rm = rad.z[i];
                exp_rm = 1./(rm + 1.0);
                CAPE[i+1] = CAPE[i] + g * step[i]
                    /exp_rm * (t_vir - t_vir_env)/t_vir_env;
// water vapour in the updraft
                if(q_v_u.x[i][j][k] >= q_Rain_add){
                    q_c_u.x[i][j][k] = q_v_u.x[i][j][k] - q_Rain_add;
                    if(q_c_u.x[i][j][k] >= 0.0) 
                        q_v_u.x[i][j][k] = q_Rain_add;
                }
// condensation within the updraft
                if((c.x[i][j][k] - q_Rain_add) >= 0.0)
                    c_u.x[i][j][k] = (c.x[i][j][k] - q_Rain_add) 
                        * M_u.x[i][j][k]/(r_humid.x[i][j][k] * step[i]);  // in kg/(kg*s)
// specification of entrainment in the updraft
//                E_u.x[i][j][k] = E_u.x[i][j][k] // COSMO data
//                    + eps_u * M_u.x[i][j][k]; // in kg/(m³s)
                E_u.x[i][j][k] = eps_u_shallow * M_u.x[i][j][k]
                    /r_humid.x[i][j][k] * (1.3 - HumidityRelative/100.0) 
                    * pow(q_Rain/q_Rain_base, 3.0);

                if(E_u.x[i][j][k] > 0.01) E_u.x[i][j][k] = 0.01;
                if(E_u.x[i][j][k] < - 0.01) E_u.x[i][j][k] = - 0.01;

// specification of detrainment in the updraft
/*
// COSMO data
                if(i == i_lfs)
                    D_u.x[i][j][k] = (1. - b_u) 
                        * M_u.x[i][j][k]/step[i]
                        + del_u * M_u.x[i][j][k];
                if(i == i_lfs+1)
                    D_u.x[i][j][k] = b_u 
                        * M_u.x[i][j][k]/step[i]
                        + del_u * M_u.x[i][j][k];  // in kg/(m³s)
*/
// ECMWF approach
                if(q_buoy_base < 0.0){ // search for non-buoyant situation
                    if(M_u.x[i][j][k] == 0.0)  M_u.x[i][j][k] = 1e-6;
                    K_u[i+1] = K_u[i] + step[i]  // kinetic energy
                        * (- E_u.x[i][j][k]/M_u.x[i][j][k] 
                        * (1.0 + bet_d * c_d) 
                        * 2.0 * K_u[i] + g/(f_t * (1.0 + gam_d_K)) 
                        * (t_u_add - t_u)/t_u);
                    K_u_sqrt = sqrt(K_u[i]/K_u[i+1]);
                    del_M_u = (1.6 - HumidityRelative/100.0) * K_u_sqrt;
                    D_u_2 = M_u.x[i][j][k]/(r_humid.x[i][j][k] 
                        * step[i]) * (1.0/del_M_u - 1.0);
                }else{
                    K_u[i] = 0.0;
                    K_u[i+1] = 0.0;
                    K_u_sqrt = 0.0;
                    del_M_u = 0.0;
                    D_u_2 = 0.0;
                }
                D_u.x[i][j][k] = del_u_1 * M_u.x[i][j][k]  // in kg/(m³s)
                    /r_humid.x[i][j][k] + D_u_2;
                if(t_u < t_00)  D_u.x[i][j][k] = 0.0;

                if(D_u.x[i][j][k] > 0.01) D_u.x[i][j][k] = 0.01;
                if(D_u.x[i][j][k] < - 0.01) D_u.x[i][j][k] = - 0.01;

// evaporation of cloud water in the updraft
                e_l.x[i][j][k] = D_u.x[i][j][k]/r_humid.x[i][j][k]
                    * q_c_u.x[i][j][k];  // in kg/(kg*s)
// formation of precipitation within the updraft
/*
// COSMO data
                if(is_land(h, 0, j, k)) delta_i_c = 3000.0; // over land in m
                if(is_air(h, 0, j, k))  delta_i_c = 1500.0; // over ocean in m
                height = get_layer_height(i); // in m
                height_base = get_layer_height(i_base); // in m
                if(height > height_base + delta_i_c)  K_p = bet_p;  // in 1/s
                else  K_p = 0.0;
 */
                if(t_u_add >= t_00){
//                    g_p.x[i][j][k] = K_p * q_c_u.x[i][j][k];  // in kg/(kg*s)// COSMO approach
                    double c_00 = 1.5e-3; // 1/s
                    double q_krit = 5.0e-4; // kg/kg
                    double w_u_u = 10.0; // m/s
                    double w_u_g = 0.0;
                    w_u_g = fabs(u_u.x[i][j][k]) * u_0;
                    if(w_u_g >= w_u_u)  w_u_g = w_u_u;
                    if(w_u_g <= 0.0)  w_u_g = 1.0e-6;
                    if(w_u_g <= 1.0)  w_u_g = 1.0;
                    g_p.x[i][j][k] = M_u.x[i][j][k]/r_humid.x[i][j][k]  // original ECMWF formula
//                        * (c_00/(0.75 * w_u_u)) 
                        * (c_00/(0.75 * w_u_g))
                        * q_c_u.x[i][j][k] 
                        * (1.0 - exp(- pow(q_c_u.x[i][j][k]/q_krit, 2.0)));  // in kg/(kg*s)
                }
                else  g_p.x[i][j][k] = 0.0;
//                if(g_p.x[i][j][k] < 0.0)  g_p.x[i][j][k] = 0.0;
// mass flux within the updraft
                M_u.x[i+1][j][k] = M_u.x[i][j][k] 
                    - step[i] * (E_u.x[i][j][k] - D_u.x[i][j][k]);    // in kg/(m²s) integration from cloud base upwards
                M_u_denom = M_u.x[i+1][j][k];
                M_u_num = M_u.x[i][j][k] - step[i] * D_u.x[i][j][k];
                if(M_u_denom == 0.0)  M_u_denom = 1e-6;

                if(M_u.x[i+1][j][k] > 5.0) M_u.x[i+1][j][k] = 5.0;
                if(M_u.x[i+1][j][k] < - 5.0) M_u.x[i+1][j][k] = - 5.0;

// dry entropy within the updraft
                s_u.x[i+1][j][k] = (M_u_num * s_u.x[i][j][k] 
                    - step[i] * (E_u.x[i][j][k] * s.x[i][j][k] 
                    - D_u.x[i][j][k] * s_u.x[i][j][k]
                    + L_latent * r_humid.x[i][j][k] * c_u.x[i][j][k])/s_0)
                    /M_u_denom;
// cloud water within the updraft
                q_v_u.x[i+1][j][k] = (M_u_num * q_v_u.x[i][j][k]  // in kg/kg
                    - step[i] * (E_u.x[i][j][k] * c.x[i][j][k] 
                    - r_humid.x[i][j][k] * c_u.x[i][j][k]))/M_u_denom;
                if(q_v_u.x[i+1][j][k] >= q_Rain_add)
                    q_v_u.x[i+1][j][k] = q_Rain_add;
// cloud water within the updraft
                q_c_u.x[i+1][j][k] = (M_u_num * q_c_u.x[i][j][k]   // in kg/kg
                    - step[i] * r_humid.x[i][j][k] 
                    * (c_u.x[i][j][k] - g_p.x[i][j][k]))/M_u_denom;
// velocity components within the updraft
                u_u.x[i+1][j][k] = u.x[i][j][k] 
                    + M_u.x[i][j][k]/(a_u * r_humid.x[i][j][k]); 
                if(u_u.x[i+1][j][k] > 2.5) u_u.x[i+1][j][k] = 2.5;
                if(u_u.x[i+1][j][k] < - 2.5) u_u.x[i+1][j][k] = - 2.5;
                v_u.x[i+1][j][k] = (M_u_num * v_u.x[i][j][k] 
                    - step[i] * E_u.x[i][j][k] * v.x[i][j][k])/M_u_denom;  // in m/s
                if(v_u.x[i+1][j][k] > 0.1) v_u.x[i+1][j][k] = 0.1;
                if(v_u.x[i+1][j][k] < - 0.1) v_u.x[i+1][j][k] = - 0.1;
                w_u.x[i+1][j][k] = (M_u_num * w_u.x[i][j][k] 
                    - step[i] * E_u.x[i][j][k] * w.x[i][j][k])/M_u_denom;  // in m/s
                if(w_u.x[i+1][j][k] > 10.0) w_u.x[i+1][j][k] = 10.0;
                if(w_u.x[i+1][j][k] < -10.0) w_u.x[i+1][j][k] = - 10.0;

//                if(is_land(h, i, j, k)||(t_u <= t_00)){
                if(is_land(h, i+1, j, k)||(t_u <= t_00)){
                    E_u.x[i+1][j][k] = 0.0;
                    D_u.x[i+1][j][k] = 0.0;
                    M_u.x[i+1][j][k] = 0.0;
                    s_u.x[i+1][j][k] = 0.0;
                    q_v_u.x[i+1][j][k] = 0.0;
                    q_c_u.x[i+1][j][k] = 0.0;
                    u_u.x[i+1][j][k] = 0.0;
                    v_u.x[i+1][j][k] = 0.0;
                    w_u.x[i+1][j][k] = 0.0;
                    g_p.x[i][j][k] = 0.0;
                    c_u.x[i][j][k] = 0.0;
                    e_l.x[i][j][k] = 0.0;

                    E_u.x[i][j][k] = 0.0;
                    D_u.x[i][j][k] = 0.0;
                    M_u.x[i][j][k] = 0.0;
                    s_u.x[i][j][k] = 0.0;
                    q_v_u.x[i][j][k] = 0.0;
                    q_c_u.x[i][j][k] = 0.0;
                    u_u.x[i][j][k] = 0.0;
                    v_u.x[i][j][k] = 0.0;
                    w_u.x[i][j][k] = 0.0;
                }
/*
            cout.precision(9);
            if((j == 90) &&(k == 180))  cout << endl 
                << "  MoistConvectionMidL   UPDRAFT     +++++++++++++++++++++++++++++++++++++++++++++++" << endl
                << "  Ma = " << (int)*get_current_time() << endl
                << "   i = " << i 
                << "  j = " << j << "  k = " << k << endl 
                << "   i_base = " << (int)i_Base.y[j][k] 
                << "   i_lfs = " << (int)i_LFS.y[j][k] << endl
                << "  s = " << s.x[i][j][k] 
                << "  s_u = " << s_u.x[i][j][k] << endl
                << "  E_u = " << E_u.x[i][j][k] 
                << "  D_u = " << D_u.x[i][j][k] << endl
                << "  e_l = " << e_l.x[i][j][k] 
                << "  c_u = " << c_u.x[i][j][k] 
                << "  g_p = " << g_p.x[i][j][k] << endl
                << "  E_Rain = " << E_Rain 
                << "  E_Rain_add = " << E_Rain_add 
                << "  q_Rain = " << q_Rain 
                << "  q_Rain_add = " << q_Rain_add << endl
                << "  q_buoy_base = " << q_buoy_base 
                << "  c = " << c.x[i][j][k] 
                << "  q_v_u = " << q_v_u.x[i][j][k] 
                << "  cloud = " << cloud.x[i][j][k] 
                << "  q_c_u = " << q_c_u.x[i][j][k] << endl
                << "  u = " << u.x[i][j][k] 
                << "  u_u = " << u_u.x[i][j][k] 
                << "  v = " << v.x[i][j][k] 
                << "  v_u = " << v_u.x[i][j][k] 
                << "  w = " << w.x[i][j][k] 
                << "  w_u = " << w_u.x[i][j][k] << endl
                << "  t_u = " << t_u - t_0 
                << "  t_u_add = " << t_u_add - t_0
                << "  t_vir = " << t_vir - t_0 << endl
                << "  HumidityRelative = " << HumidityRelative 
                << "  r_dry = " << r_dry.x[i][j][k] 
                << "  r_humid = " << r_humid.x[i][j][k] 
                << "  r_humid_add = " << r_humid_add << endl
                << "  p_stat = " << p_stat.x[i][j][k] 
                << "  p_hydro = " << p_hydro.x[i][j][k] 
                << "  height = " << height 
//                << "  height_base = " << height_base 
                << "  step = " << step[i] << endl
                << "  M_u = " << M_u.x[i][j][k] 
                << "  M_d_LFS = " << M_d_LFS[j][k] << endl
                << "  rad = " << rad.z[i]
                << "  step_non_dim = " << g * (rad.z[i+1] - rad.z[i])/exp_rm
                << "  temp_non_dim = " << (t_u_add - t_u)/t_u
                << "  CAPE = " << CAPE[i] << endl
                << "  del_M_u = " << del_M_u
                << "  D_u_2 = " << D_u_2
                << "  D_u = " << D_u.x[i][j][k]
                << "  K_u = " << K_u[i] << endl << endl;
*/
                }
            } // end i convective updraft
        }  // end j convective updraft
    }  // end k convective updraft
//
// CONVECTIVE DOWNDRAFT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            i_base = (int)i_Base.y[j][k];
            i_lfs = (int)i_LFS.y[j][k];
            M_d.x[i_lfs][j][k] = M_d_LFS[j][k];
//            height_lfs = get_layer_height(i_lfs);
            P_conv_shall.x[i_lfs][j][k] = 0.0;
            t_u_add = t_add_u * t.x[i_base][j][k] * t_0; // in K
            E_Rain_add = hp * exp_func(t_u_add, 17.2694, 35.86);
            K_d = std::vector<double>(im, 0.0);
            K_d[i_lfs] = 0.5 * u_d.x[i_lfs][j][k] * u_d.x[i_lfs][j][k];
            for(int i = i_lfs; i > 0; i--){
                height = get_layer_height(i);
                step[i] = get_layer_height(i) - get_layer_height(i-1);  // local atmospheric shell thickness
                t_u = t.x[i][j][k] * t_0; // in K
                e = c.x[i][j][k] * p_hydro.x[i][j][k]/ep;  // water vapour pressure in hPa
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                t_u_add = t_add_u * t_0 + t_u; // in K
                t_vir = t_u_add * (1.0 + alf * q_v_d.x[i][j][k]
                    - cloud.x[i][j][k]);
                E_Rain_add = hp * exp_func(t_vir, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_hydro.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Rain_add = ep * E_Rain_add/(p_hydro.x[i][j][k] - E_Rain_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                r_humid_add = 1e2 * p_hydro.x[i][j][k]
                    /(R_Air * (1.0 + (R_W_R_A - 1.) * q_v_d.x[i][j][k]) 
                    * t_u_add); 
                q_buoy_base = r_humid_add - r_humid.x[i][j][k];
                HumidityRelative = e/E_Rain * 100.0;
                if(t_u >= t_0) L_latent = lv;
                else L_latent = ls;
                s.x[i][j][k] = (cp_l * t_u + g * height 
                    + lv * c.x[i][j][k])/s_0;
// evaporation of precipitation below the cloud base
                if(i <= i_base){
//                    e_p.x[i][j][k] = C_p * alf_1   // in kg/(kg*s) // COSMO data
//                        * (q_Rain - q_v_d.x[i][j][k])
//                        * sqrt(p_ps/alf_2 * P_conv_shall.x[i][j][k]/C_p); // original COSMO formula
// ECMWF approach
                    double p_ps = p_hydro.x[i][j][k]/p_hydro.x[0][j][k];
                    if(is_land(h, 0, j, k))
                        e_p.x[i][j][k] = sig * alf_11 
                            * (RH_crit_land * q_Rain - c.x[i][j][k])
                            * pow((sqrt(p_ps)/alf_22 * P_conv_shall.x[i][j][k]/sig), 
                            alf_33); // original ECMWF formula
                    else 
                        e_p.x[i][j][k] = sig * alf_11 
                            * (RH_crit_water * q_Rain - c.x[i][j][k])
                            * pow((sqrt(p_ps)/alf_22 * P_conv_shall.x[i][j][k]/sig), 
                            alf_33); // original ECMWF formula
                }
                if(e_p.x[i][j][k] <= 0.0) e_p.x[i][j][k] = 0.0;
// evaporation of rain precipitation within the downdraft
                if((t_u_add >= t_0)&&(q_v_d.x[i][j][k] <= q_Rain_add))
                    e_d.x[i][j][k] = a_ev 
                        * (1.0 + b_ev * pow(P_conv_shall.x[i][j][k],(1.0/6.0)))  // evaporation of rain due to water vapour diffusion, < XIII >
                        * (q_Rain_add - q_v_d.x[i][j][k]) 
                        * pow(P_conv_shall.x[i][j][k],(4.0/9.0));
                if(e_d.x[i][j][k] < 0.0) e_d.x[i][j][k] = 0.0;
// specification of entrainment and detrainment in the up/downdraft
// COSMO approach
//                E_d.x[i][j][k] = eps_d * fabs(M_d.x[i][j][k]);
//                D_d.x[i][j][k] = del_d * fabs(M_d.x[i][j][k]);

// ECMWF approach
                t_u_base = t.x[i_base][j][k] * t_0;
                E_Rain_base = hp * exp_func(t_u_base, 17.2694, 35.86);
                q_Rain_base = ep * E_Rain_base
                    /(p_hydro.x[i_base][j][k] - E_Rain_base); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                E_d.x[i][j][k] = eps_u_shallow * M_d.x[i][j][k]
                    /r_humid.x[i][j][k] * (1.3 - HumidityRelative/100.0) 
                    * pow(q_Rain/q_Rain_base, 3.0);

                if(E_d.x[i][j][k] > 0.01) E_d.x[i][j][k] = 0.01;
                if(E_d.x[i][j][k] < - 0.01) E_d.x[i][j][k] = - 0.01;

// ECMWF kinetic energy (vertical velocity of the downdraft)
                if(q_buoy_base < 0.0){ // search for non-buoyant situation
                    if(M_d.x[i][j][k] == 0.0)  M_d.x[i][j][k] = 1e-6;
                    K_d[i-1] = K_d[i] + step[i]  // kinetic energy
                        * (- E_d.x[i][j][k]/M_d.x[i][j][k] 
                        * (1.0 + bet_d * c_d) 
                        * 2.0 * K_d[i] + g/(f_t * (1.0 + gam_d_K)) 
                        * (t_u_add - t_u)/t_u);
                    K_d_sqrt = sqrt(K_d[i-1]/K_d[i]);
                    del_M_d = (1.6 - HumidityRelative/100.0) * K_d_sqrt;
                    D_u_2 = M_d.x[i][j][k]/(r_humid.x[i][j][k] 
                        * step[i]) * (1.0/del_M_d - 1.0);
                }else{
                    K_d[i] = 0.0;
                    K_d[i-1] = 0.0;
                    K_d_sqrt = 0.0;
                    del_M_d = 0.0;
                    D_d_2 = 0.0;
                }
                D_d.x[i][j][k] = del_u_1 * M_d.x[i][j][k]  // in kg/(m³s)
                    /r_humid.x[i][j][k] + D_d_2;

                if(D_d.x[i][j][k] > 0.01) D_d.x[i][j][k] = 0.01;
                if(D_d.x[i][j][k] < - 0.01) D_d.x[i][j][k] = - 0.01;

// mass flux within the downdraft
                M_d.x[i-1][j][k] = M_d.x[i][j][k] 
                    + step[i] * (E_d.x[i][j][k] - D_d.x[i][j][k]);  // integration from LFS downwards
                M_d_denom = M_d.x[i-1][j][k];
                M_d_num = M_d.x[i][j][k] + step[i] * D_d.x[i][j][k];
                if(M_d_denom == 0.0)  M_d_denom = 1e-6;

                if(M_d.x[i-1][j][k] > 5.0) M_d.x[i-1][j][k] = 5.0;
                if(M_d.x[i-1][j][k] < - 5.0) M_d.x[i-1][j][k] = - 5.0;

// dry entropy within the downdraft
                if(t_u >= t_0) L_latent = lv;
                else L_latent = ls;
                s_d.x[i-1][j][k] = (M_d_num * s_d.x[i][j][k] 
                    - step[i] * (E_d.x[i][j][k] * s.x[i][j][k] 
                    - D_d.x[i][j][k] * s_d.x[i][j][k]
                   + L_latent * r_humid.x[i][j][k] * e_d.x[i][j][k])/s_0)
                    /M_d_denom;
// cloud water within the downdraft
                q_v_d.x[i-1][j][k] = (M_d_num * q_v_d.x[i][j][k] 
                    + step[i] * (E_d.x[i][j][k] * c.x[i][j][k] 
                    + r_humid.x[i][j][k] * e_d.x[i][j][k]))/M_d_denom;
                if(q_v_d.x[i-1][j][k] >= q_Rain) 
                    q_v_d.x[i-1][j][k] = q_Rain;
                if(q_v_d.x[i-1][j][k] <= 0.0) q_v_d.x[i-1][j][k] = q_Rain;
// velocity components within the downdraft
                u_d.x[i-1][j][k] = u.x[i][j][k] 
                    + M_d.x[i-1][j][k]/(a_d * r_humid.x[i][j][k]); 
                if(u_d.x[i-1][j][k] > 2.5) u_d.x[i-1][j][k] = 2.5;
                if(u_d.x[i-1][j][k] < - 2.5) u_d.x[i-1][j][k] = - 2.5;
                v_d.x[i-1][j][k] = (M_d_num * v_d.x[i][j][k] 
                    + step[i] * E_d.x[i][j][k] * v.x[i][j][k])/M_d_denom;
                if(v_d.x[i+1][j][k] > 0.1) v_d.x[i+1][j][k] = 0.1;
                if(v_d.x[i+1][j][k] < -0.1) v_d.x[i+1][j][k] = - 0.1;
                w_d.x[i-1][j][k] = (M_d_num * w_d.x[i][j][k] 
                    + step[i] * E_d.x[i][j][k] * w.x[i][j][k])/M_d_denom;
                if(w_d.x[i+1][j][k] > 10.0) w_d.x[i+1][j][k] = 10.0;
                if(w_d.x[i+1][j][k] < -10.0) w_d.x[i+1][j][k] = - 10.0;
// rain water formed by cloud convection
                P_conv_shall.x[i-1][j][k] = P_conv_shall.x[i][j][k]   // in kg/(m²*s) == mm/s
                    + step[i] * r_humid.x[i][j][k]
                    * (g_p.x[i][j][k]   // in kg/(kg*s)
                    - e_d.x[i][j][k]   // in kg/(kg*s)
                    - e_p.x[i][j][k]);  // in kg/(kg*s)
                if(P_conv_shall.x[i-1][j][k] <= 0.0) P_conv_shall.x[i-1][j][k] = 0.0;
                if(t.x[i-1][j][k] * t_0 <= t_0) P_conv_shall.x[i-1][j][k] = 0.0;
                if(P_conv_shall.x[i-1][j][k] > 8.0/8.64e4)  P_conv_shall.x[i-1][j][k] = 8.0/8.64e4; //assumed max value
                if(P_conv_shall.x[i-1][j][k] > 0.0){
                    moistconv = true;
                    if(P_conv_shall.x[i-1][j][k] > maxValue){
                        maxValue = P_conv_shall.x[i-1][j][k];
                        P_moistconvshall = maxValue * 8.64e4; // in mm/d
                        i_moistconv = i;
                        j_moistconv = j;
                        k_moistconv = k;
                        height_moistconv = get_layer_height(i-1);
/*
                        cout.precision(9);
                        cout << "      max values of MoistConvectionShall " 
                        << endl
                        << "      i_moistconvshall = " << i_moistconv
                        << "      j_moistconvshall = " << j_moistconv
                        << "   k_moistconvshall = " << k_moistconv
                        << "   height_moistconvshall = " << height_moistconv
                        << "   P_moistconvshall = " << P_moistconvshall << endl;
*/
                    }
                }

//                if(is_land(h, i, j, k)||(t_u <= t_00)){
                if(is_land(h, i-1, j, k)||(t_u <= t_00)){
                    E_d.x[i-1][j][k] = 0.0;
                    D_d.x[i-1][j][k] = 0.0;
                    M_d.x[i-1][j][k] = 0.0;
                    s_d.x[i-1][j][k] = 1.0;
                    q_v_d.x[i-1][j][k] = 0.0;
                    u_d.x[i-1][j][k] = 0.0;
                    v_d.x[i-1][j][k] = 0.0;
                    w_d.x[i-1][j][k] = 0.0;
                    e_d.x[i][j][k] = 0.0;
                    e_p.x[i][j][k] = 0.0;

                    E_d.x[i][j][k] = 0.0;
                    D_d.x[i][j][k] = 0.0;
                    M_d.x[i][j][k] = 0.0;
                    s_d.x[i][j][k] = 1.0;
                    q_v_d.x[i][j][k] = 0.0;
                    u_d.x[i][j][k] = 0.0;
                    v_d.x[i][j][k] = 0.0;
                    w_d.x[i][j][k] = 0.0;
                }
/*
            cout.precision(9);
            if((j == 90) &&(k == 180))  cout << endl 
                << "  MoistConvectionMidL   DOWNDRAFT     -------------------------------------------" << endl
                << "  Ma = " << (int)*get_current_time() << endl
                << "  i = " << i 
                << "  j = " << j<< "  k = " << k << endl 
                << "  i_base = " << int(i_Base.y[j][k]) 
                << "  i_lfs = " << int(i_LFS.y[j][k]) << endl
                << "  s = " << s.x[i][j][k] 
                << "  s_d = " << s_d.x[i][j][k] << endl
                << "  E_d = " << E_d.x[i][j][k] 
                << "  D_d = " << D_d.x[i][j][k] << endl
                << "  e_d = " << e_d.x[i][j][k] 
                << "  e_l = " << e_l.x[i][j][k] 
                << "  e_p = " << e_p.x[i][j][k] 
                << "  c_u = " << c_u.x[i][j][k] 
                << "  g_p = " << g_p.x[i][j][k] << endl
                << "  P_conv_mm/s = " << P_conv_midl.x[i][j][k]
                << "  P_conv_mm/d = " << P_conv_midl.x[i][j][k] * 8.64e4 << endl
                << "  P_conv_shall_mm/s = " << P_conv_shall.x[i][j][k]
                << "  P_conv_shall_mm/d = " << P_conv_shall.x[i][j][k] * 8.64e4 << endl
                << "  P_rain_mm/s = " << P_rain.x[i][j][k]
                << "  P_rain_mm/d = " << P_rain.x[i][j][k] * 8.64e4 << endl
                << "  P_snow_mm/s = " << P_snow.x[i][j][k]
                << "  P_snow_mm/d = " << P_snow.x[i][j][k] * 8.64e4 << endl
                << "  E_Rain = " << E_Rain 
                << "  q_Rain_add = " << q_Rain_add
                << "  q_Rain_base = " << q_Rain_base
                << "  q_Rain = " << q_Rain << endl
                << "  q_buoy_base = " << q_buoy_base 
                << "  c = " << c.x[i][j][k] 
                << "  q_v_d = " << q_v_d.x[i][j][k] 
                << "  cloud = " << cloud.x[i][j][k] << endl
                << "  u = " << u.x[i][j][k] 
                << "  u_d = " << u_d.x[i][j][k] 
                << "  v = " << v.x[i][j][k] 
                << "  v_d = " << v_d.x[i][j][k] 
                << "  w = " << w.x[i][j][k] 
                << "  w_d = " << w_d.x[i][j][k] << endl
                << "  t_u = " << t_u - t_0 << endl
                << "  HumidityRelative = " << HumidityRelative 
                << "  r_dry = " << r_dry.x[i][j][k] 
                << "  r_humid = " << r_humid.x[i][j][k] 
                << "  r_humid_add = " << r_humid_add << endl
                << "  p_stat = " << p_stat.x[i][j][k] 
                << "  p_hydro = " << p_hydro.x[i][j][k] 
                << "  height = " << height 
//                << "  height_lfs = " << height_lfs 
                << "  step = " << step[i] << endl
                << "  M_d = " << M_d.x[i][j][k] 
                << "  M_d_LFS = " << M_d_LFS[j][k] << endl
//                << "  del_M_d = " << del_M_d
//                << "  D_d_2 = " << D_d_2
                << "  D_d = " << D_d.x[i][j][k]
                << "  K_d = " << K_d[i] << endl << endl;
*/
            } // end i convective downdraft
            e_d.x[0][j][k] = c43 * e_d.x[1][j][k] - c13 * e_d.x[2][j][k];
            e_p.x[0][j][k] = c43 * e_p.x[1][j][k] - c13 * e_p.x[2][j][k];
            E_d.x[0][j][k] = c43 * E_d.x[1][j][k] - c13 * E_d.x[2][j][k];
            D_d.x[0][j][k] = c43 * D_d.x[1][j][k] - c13 * D_d.x[2][j][k];
        }  // end j convective downdraft
    }  // end k convective downdraft
/*
*
*/
    AtomUtils::fft_gaussian_filter_3d(P_conv_shall,1);
/*
*
*/

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
// RHS for thermodynamic forcing due to moist convection
            for(int i = 0; i < im-1; i++){
                step[i] = get_layer_height(i+1) - get_layer_height(i);  // local atmospheric shell thickness
                if((t.x[i][j][k] * t_0) >= t_0) L_latent = lv;
                else L_latent = ls;
                MC_t.x[i][j][k] = 
                    - ((M_u.x[i][j][k] * (s_u.x[i][j][k] - s.x[i][j][k]) 
                    + M_d.x[i][j][k] * (s_d.x[i][j][k] - s.x[i][j][k])) 
                    - (M_u.x[i][j][k] * (s_u.x[i][j][k] - s.x[i][j][k]) 
                    + M_d.x[i][j][k] * (s_d.x[i][j][k] - s.x[i][j][k]))) 
                    /(r_humid.x[i][j][k] * step[i])/cp_l * s_0
                    + L_latent/cp_l * (c_u.x[i][j][k] - e_d.x[i][j][k] 
                    - e_l.x[i][j][k] - e_p.x[i][j][k]);  // in (K)(1/s)
                MC_q.x[i][j][k] = 
                    - ((M_u.x[i][j][k] * (q_v_u.x[i][j][k] - c.x[i][j][k]) 
                    + M_d.x[i][j][k] * (q_v_d.x[i][j][k] - c.x[i][j][k])) 
                    - (M_u.x[i][j][k] * (q_v_u.x[i][j][k] - c.x[i][j][k]) 
                    + M_d.x[i][j][k] * (q_v_d.x[i][j][k] - c.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k])
                    - (c_u.x[i][j][k] - e_d.x[i][j][k] - e_l.x[i][j][k] 
                    - e_p.x[i][j][k]);  // in (kg/kg)(1/s)
                MC_v.x[i][j][k] = 
                    - ((M_u.x[i][j][k] * (v_u.x[i][j][k] - v.x[i][j][k]) 
                    + M_d.x[i][j][k] * (v_d.x[i][j][k] - v.x[i][j][k])) 
                    - (M_u.x[i][j][k] * (v_u.x[i][j][k] - v.x[i][j][k]) 
                    + M_d.x[i][j][k] * (v_d.x[i][j][k] - v.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k]) * u_0;  // in (m/s)(1/s)
                MC_w.x[i][j][k] = 
                    - ((M_u.x[i][j][k] * (w_u.x[i][j][k] - w.x[i][j][k]) 
                    + M_d.x[i][j][k] * (w_d.x[i][j][k] - w.x[i][j][k])) 
                    - (M_u.x[i][j][k] * (w_u.x[i][j][k] - w.x[i][j][k]) 
                    + M_d.x[i][j][k] * (w_d.x[i][j][k] - w.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k]) * u_0;  // in (m/s)(1/s)
            }  // MC-loop
        } // j-loop
    } // k-loop
// boundaries of various variables
    for(int k = 1; k < km-1; k++){
        for(int j = 2; j < jm-2; j++){
            M_u.x[0][j][k] = c43 * M_u.x[1][j][k] -
                c13 * M_u.x[2][j][k];
            M_u.x[im-1][j][k] = c43 * M_u.x[im-2][j][k] -
                c13 * M_u.x[im-3][j][k];
            M_d.x[0][j][k] = c43 * M_d.x[1][j][k] -
                c13 * M_d.x[2][j][k];
            M_d.x[im-1][j][k] = c43 * M_d.x[im-2][j][k] -
                c13 * M_d.x[im-3][j][k];
            s.x[0][j][k] = c43 * s.x[1][j][k] -
                c13 * s.x[2][j][k];
            s.x[im-1][j][k] = c43 * s.x[im-2][j][k] -
                c13 * s.x[im-3][j][k];
            s_u.x[0][j][k] = c43 * s_u.x[1][j][k] -
                c13 * s_u.x[2][j][k];
            s_u.x[im-1][j][k] = c43 * s_u.x[im-2][j][k] -
                c13 * s_u.x[im-3][j][k];
            s_d.x[0][j][k] = c43 * s_d.x[1][j][k] -
                c13 * s_d.x[2][j][k];
            s_d.x[im-1][j][k] = c43 * s_d.x[im-2][j][k] -
                c13 * s_d.x[im-3][j][k];
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int i = 0; i < im; i++){
            M_u.x[i][0][k] = c43 * M_u.x[i][1][k] -
                 c13 * M_u.x[i][2][k];
            M_u.x[i][jm-1][k] = c43 * M_u.x[i][jm-2][k] -
                 c13 * M_u.x[i][jm-3][k];
            M_d.x[i][0][k] = c43 * M_d.x[i][1][k] -
                 c13 * M_d.x[i][2][k];
            M_d.x[i][jm-1][k] = c43 * M_d.x[i][jm-2][k] -
                 c13 * M_d.x[i][jm-3][k];
            s.x[i][0][k] = c43 * s.x[i][1][k] -
                 c13 * s.x[i][2][k];
            s.x[i][jm-1][k] = c43 * s.x[i][jm-2][k] -
                 c13 * s.x[i][jm-3][k];
            s_u.x[i][0][k] = c43 * s_u.x[i][1][k] -
                 c13 * s_u.x[i][2][k];
            s_u.x[i][jm-1][k] = c43 * s_u.x[i][jm-2][k] -
                 c13 * s_u.x[i][jm-3][k];
            s_d.x[i][0][k] = c43 * s_d.x[i][1][k] -
                 c13 * s_d.x[i][2][k];
            s_d.x[i][jm-1][k] = c43 * s_d.x[i][jm-2][k] -
                 c13 * s_d.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 1; j < jm-1; j++){
            M_u.x[i][j][0] = c43 * M_u.x[i][j][1] -
                c13 * M_u.x[i][j][2];
            M_u.x[i][j][km-1] = c43 * M_u.x[i][j][km-2] -
                c13 * M_u.x[i][j][km-3];
            M_u.x[i][j][0] = M_u.x[i][j][km-1] =
               (M_u.x[i][j][0] + M_u.x[i][j][km-1])/2.;
            M_d.x[i][j][0] = c43 * M_d.x[i][j][1] -
                c13 * M_d.x[i][j][2];
            M_d.x[i][j][km-1] = c43 * M_d.x[i][j][km-2] -
                c13 * M_d.x[i][j][km-3];
            M_d.x[i][j][0] = M_d.x[i][j][km-1] =
               (M_d.x[i][j][0] + M_d.x[i][j][km-1])/2.;
            s.x[i][j][0] = c43 * s.x[i][j][1] -
                c13 * s.x[i][j][2];
            s.x[i][j][km-1] = c43 * s.x[i][j][km-2] -
                c13 * s.x[i][j][km-3];
            s.x[i][j][0] = s.x[i][j][km-1] =
               (s.x[i][j][0] + s.x[i][j][km-1])/2.;
            s_u.x[i][j][0] = c43 * s_u.x[i][j][1] -
                c13 * s_u.x[i][j][2];
            s_u.x[i][j][km-1] = c43 * s_u.x[i][j][km-2] -
                c13 * s_u.x[i][j][km-3];
            s_u.x[i][j][0] = s_u.x[i][j][km-1] =
               (s_u.x[i][j][0] + s_u.x[i][j][km-1])/2.;
            s_d.x[i][j][0] = c43 * s_d.x[i][j][1] -
                c13 * s_d.x[i][j][2];
            s_d.x[i][j][km-1] = c43 * s_d.x[i][j][km-2] -
                c13 * s_d.x[i][j][km-3];
            s_d.x[i][j][0] = s_d.x[i][j][km-1] =
               (s_d.x[i][j][0] + s_d.x[i][j][km-1])/2.;
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if(is_land(h, i, j, k)){
                    P_conv_shall.x[i][j][k] = P_conv_shall.x[i_mount][j][k];
                    M_u.x[i][j][k] = M_u.x[i_mount][j][k];
                    M_d.x[i][j][k] = M_d.x[i_mount][j][k];
                    MC_t.x[i][j][k] = MC_t.x[i_mount][j][k];
                    MC_q.x[i][j][k] = MC_q.x[i_mount][j][k];
                    MC_v.x[i][j][k] = MC_v.x[i_mount][j][k];
                    MC_w.x[i][j][k] = MC_w.x[i_mount][j][k];
                    s.x[i][j][k] = s.x[i_mount][j][k];
                    s_u.x[i][j][k] = s_u.x[i_mount][j][k];
                    s_d.x[i][j][k] = s_d.x[i_mount][j][k];
                    q_v_u.x[i][j][k] = q_v_u.x[i_mount][j][k];
                    q_v_d.x[i][j][k] = q_v_d.x[i_mount][j][k];
                    q_c_u.x[i][j][k] = q_c_u.x[i_mount][j][k];
                    E_u.x[i][j][k] = E_u.x[i_mount][j][k];
                    E_d.x[i][j][k] = E_d.x[i_mount][j][k];
                    D_u.x[i][j][k] = D_u.x[i_mount][j][k];
                    D_d.x[i][j][k] = D_d.x[i_mount][j][k];
                    c_u.x[i][j][k] = c_u.x[i_mount][j][k];
                    g_p.x[i][j][k] = g_p.x[i_mount][j][k];
                    e_l.x[i][j][k] = e_l.x[i_mount][j][k];
                    e_p.x[i][j][k] = e_p.x[i_mount][j][k];
                    e_d.x[i][j][k] = e_d.x[i_mount][j][k];
                    u_u.x[i][j][k] = u_u.x[i_mount][j][k];
                    u_d.x[i][j][k] = u_d.x[i_mount][j][k];
                    v_u.x[i][j][k] = v_u.x[i_mount][j][k];
                    v_d.x[i][j][k] = v_d.x[i_mount][j][k];
                    w_u.x[i][j][k] = w_u.x[i_mount][j][k];
                    w_d.x[i][j][k] = w_d.x[i_mount][j][k];
                }
            }
        }
    }

    if(moistconv == false)  
        cout << "      no saturation of water vapour in MoistConvectionShall found" 
        << endl;
    else
        cout << "      saturation of water vapour in MoistConvectionShall found" 
        << endl
        << "      i_moistconvshall = " << i_moistconv
        << "      j_moistconvshall = " << j_moistconv
        << "   k_moistconvshall = " << k_moistconv
        << "   height_moistconvshall[m] = " << height_moistconv
        << "   P_moistconvshall[mm/d] = " << P_moistconvshall << endl;
    cout << "      MoistConvectionShall ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::WaterVapourEvaporation(){ 
    cout << endl << "      WaterVapourEvaporation" << endl;
// preparations for water vapour increase due to the differences between evaporation and precipitation 
// procedure given in Rui Xin Huang, Ocean Circulation, p. 165
// amount of additional water vapour by the difference of evaporation and precipitation is negligible but functional
    rad.Coordinates(im, r0, dr);
    std::vector<double> step(im, 0);
//    double c_surf_ocean = 0.8; // 80% humidity
//    double c_surf_ocean = 0.9; // 90% humidity
//    double c_surf_ocean = 1.0; // 100% humidity
//    double rm = log(get_layer_height(0)/L_atm + 1.0);
//    double rm = log(rad.z[0] + 1.0);
//    double rm = rad.z[0];
//    double exp_rm = 1./(rm + 1.0);
    double evap_precip = 0.0;
    double vapour_surface_n = 0.0;
    double vapour_surface = 0.0;
    double coeff_evap = 1.1574e-8;  // 1.1574-8 is the conversion of (Evap-Prec) from in mm/d to m/s
    double e = 0.0;  // water vapour pressure in Pa
    double e_sa = 0.0;  // water vapour pressure in Pa
    double r = 0.0;
    double u_bar = 0.0;
    double R_net = 0.0;
    double t_u = 0.0;
    double E_Rain = 0.0;
    double E_Ice = 0.0;
    double sat_deficit = 0.0;
    double del_gam = 0.0;
    double gam_del = 0.0;

// short wave radiation
    short_wave_radiation = std::vector<double>(jm, rad_pole_short);
    int j_max = jm-1;
    int j_half = j_max/2;
    double rad_short_eff = rad_pole_short - rad_equator_short;
    for(int j=j_half; j>=0; j--){
        short_wave_radiation[j] = rad_short_eff * parabola((double)j
            /(double)j_half) + rad_pole_short;
    }
    for(int j=j_max; j>j_half; j--){
        short_wave_radiation[j] = short_wave_radiation[j_max-j];
    }

  // saturation vapour pressure in the water phase for t > 0°C in hPa
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            step[0] = get_layer_height(1) - get_layer_height(0);
/*
            c.x[0][j][k] = c_surf_ocean * hp * ep * exp(17.0809 
                * (t.x[0][j][k] * t_0 - t_0)/(234.175 
                + (t.x[0][j][k] * t_0 - t_0)))/((r_air * R_Air 
                * t.x[0][j][k] * t_0) * .01);
*/
            c_fix.y[j][k] = c.x[0][j][k];
            t_u = t.x[0][j][k] * t_0;
            del_gam = 0.439 + 0.0112 * t.x[0][j][k];
            gam_del = 0.5495 + 0.01119 * t.x[0][j][k];
//            double u_bar = sqrt((v.x[1][j][k] 
//                * v.x[1][j][k] + w.x[1][j][k] 
//                * w.x[1][j][k])/2.0) * u_0;
            u_bar = sqrt((v.x[0][j][k] 
                * v.x[0][j][k] + w.x[0][j][k] 
                * w.x[0][j][k])/2.0) * u_0;
            R_net = short_wave_radiation[j]/radiation.x[0][j][k];
            if((t.x[0][j][k] * t_0) >= t_0){
                if(is_land(h, 0, j, k)){
                    t_u = t.x[0][j][k] * t_0;
                    e = c.x[0][j][k] * p_hydro.x[0][j][k]/ep;  // water vapour pressure in hPa
                    E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                    sat_deficit = E_Rain - e;  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = 
                        C_Dalton(0, j, k, coeff_Dalton, u_0, v, w) 
                        * sat_deficit * 24.0;  // mm/h in mm/d
                    e_sa = E_Rain * 100.0;
                    r = e/e_sa;
                    Evaporation_Penman.y[j][k] = del_gam * R_net/lv // in mm/d
                        + gam_del * 0.0026 * (1.0 + 0.54 * u_bar) 
                        * (1.0 - r) * e_sa;
                } // simplified formula for Evapotranspiration by Dalton or Penman law  in kg/(m²*s) = mm/s
                if(is_water(h, 0, j, k)){
                    e = c.x[0][j][k] * p_hydro.x[0][j][k]/ep;  // water vapour pressure in hPa
                    t_u = t.x[0][j][k] * t_0;
                    E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                    sat_deficit = E_Rain - e;  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = 
                        C_Dalton(0, j, k, coeff_Dalton, u_0, v, w) 
                        * sat_deficit * 24.0;  // mm/h in mm/d
                    e_sa = E_Rain * 100.0;
                    r = e/e_sa;
                    Evaporation_Penman.y[j][k] = del_gam * R_net/lv  // in mm/d
                        + gam_del * 0.0026 * (1.0 + 0.54 * u_bar) 
                        * (1.0 - r) * e_sa;
                }
            }
            if((t.x[0][j][k] * t_0) < t_0){
                if(is_land(h, 0, j, k)){
                    t_u = t.x[0][j][k] * t_0;
                    e = c.x[0][j][k] * p_hydro.x[0][j][k]/ep;  // water vapour pressure in hPa
                    E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                    sat_deficit = (E_Ice - e);  // saturation deficit in hPa
//                    Evaporation_Dalton.y[j][k] = 
//                        C_Dalton(1, j, k, coeff_Dalton, u_0, v, w) 
//                        * sat_deficit * 24.;  // mm/h in mm/d
                    Evaporation_Dalton.y[j][k] = 0.0;
                    e_sa = E_Ice * 100.0;
                    r = e/e_sa;
                    Evaporation_Penman.y[j][k] = del_gam * R_net/lv  // in mm/d
                        + gam_del * 0.0026 * (1.0 + 0.54 * u_bar) 
                        * (1.0 - r) * e_sa;
                }
                if(is_water(h, 0, j, k)){
                    e = c.x[0][j][k] * p_hydro.x[0][j][k]/ep;  // water vapour pressure in hPa
                    t_u = t.x[0][j][k] * t_0;
                    E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                    sat_deficit = (E_Ice - e);  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = 
                        C_Dalton(0, j, k, coeff_Dalton, u_0, v, w) 
                        * sat_deficit * 24.0;  // mm/h in mm/d
                    e_sa = E_Ice * 100.0;
                    r = e/e_sa;
//                    Evaporation_Penman.y[j][k] = del_gam * R_net/lv // in mm/d evapo-transpiration
//                        + gam_del * 0.0026 * (1.0 + 0.54 * u_bar) 
//                        * (1.0 - r) * e_sa;
                    Evaporation_Penman.y[j][k] = 0.0;
                }
            }
            if(Evaporation_Dalton.y[j][k] <= 0.0)  
                Evaporation_Dalton.y[j][k] = 0.0;
            if(Evaporation_Penman.y[j][k] <= 0.0)  
                Evaporation_Penman.y[j][k] = 0.0;
            double Evaporation = Evaporation_Dalton.y[j][k] 
                + Evaporation_Penman.y[j][k];
//            evap_precip = coeff_evap * (Evaporation_Dalton.y[j][k] 
//                - 8.64e4 * Precipitation.y[j][k]); // in mm/d
//            evap_precip = coeff_evap * (Evaporation_Penman.y[j][k] 
            evap_precip = coeff_evap * (Evaporation 
                - 8.64e4 * Precipitation.y[j][k]); // in mm/s
            for(int iter_prec = 1; iter_prec <= 10; iter_prec++){ // iter_prec may be varied
                vapour_surface = + ((- 3.0 * c.x[0][j][k] 
                    + 4.0 * c.x[1][j][k] - c.x[2][j][k])/(2.0 * step[0]) 
                    * (1.0 - 2. * c.x[0][j][k]) * evap_precip);  // 1. order derivative, 2. order accurate
                if(iter_prec == 1)  vapour_surface = 0.9 * vapour_surface;
                vapour_evaporation.y[j][k] = vapour_surface;
                if(is_land(h, 0, j, k)) vapour_evaporation.y[j][k] = 0.0;
                c.x[0][j][k] = c_fix.y[j][k] + vapour_evaporation.y[j][k];
/*
                cout.precision(8);
                cout.setf(ios::fixed);
                if((j == 120)&&(k == 180)) cout << endl
                    << "  WaterVapourEvaporation" << endl
                    << "  it = " << iter_prec << endl
                    << "  j = " << j << "  k = " << k << endl
                    << "  vap_evap = " << vapour_evaporation.y[j][k] * 1e3 << endl
                    << "  vap_diff = " << fabs(vapour_surface/vapour_surface_n - 1.)  * 1e3
                    << "  vap_surf = " << vapour_surface * 1e3
                    << "  vap_surf_n = " << vapour_surface_n * 1e3 << endl
                    << "  c_fix = " << c_fix.y[j][k] * 1e3 
                    << "  c = " << c.x[0][j][k] * 1e3 << endl
                    << "  v = " << v.x[0][j][k] * u_0
                    << "  w = " << w.x[0][j][k] * u_0
                    << "  vel_magnitude = " << sqrt(v.x[0][j][k] * v.x[0][j][k] 
                        + w.x[0][j][k] * w.x[0][j][k]) * u_0 
                    << "  c_grad_1 = " << (c.x[1][j][k] - c.x[0][j][k])/step[0]
                    << "  c_grad_2 = " << (- 3. * c.x[0][j][k] + 4. * c.x[1][j][k] 
                        - c.x[2][j][k]) /(2. * step[0]) << endl
                    << "  sat_deficit = " << sat_deficit 
                    << "  e = " << e 
                    << "  E_Rain = " << E_Rain << endl
                    << "  C_Dalton = " << C_Dalton(1, j, k, coeff_Dalton, u_0, v, w)
                        * sat_deficit * 2.778e-4
                    << "  Evap-Prec = " << evap_precip 
                    << "  Evap_Dalton = " << Evaporation_Dalton.y[j][k] 
                    << "  Evap_Penman = " << Evaporation_Penman.y[j][k] 
                    << "  Evap_global = " << Evaporation 
                    << "  Prec = " << 8.64e4 * Precipitation.y[j][k] << endl;
*/
/*
                int i_trans = 3; // 2 ~ 129 m vertical extension of ocean surface evaporation, arbitrary 
                for(int i = 0; i < i_trans; i++){
                    if(i <= im-1){
//                        if(i > 0){
//                            double x = get_layer_height(i) 
//                               /get_layer_height(i_trans); 
                            double x = (double)i/(double)i_trans; 
                            c.x[i][j][k] = parabola_interp(c.x[i_trans][j][k], 
                                c.x[0][j][k], x); // parabolic transition

                            c.x[i][j][k] = (c.x[i_trans][j][k] - c.x[0][j][k]) 
                                * (get_layer_height(i)/get_layer_height(i_trans)) 
                                + c.x[0][j][k]; // linear transition

//                        }
                    }
                } // end i
*/
                if(fabs(vapour_surface/vapour_surface_n - 1.0) 
                    <= 1.e-5)  break;  
                vapour_surface_n = vapour_surface;
            } // iter_prec 
        } // end j
    } // end k
    cout << "      WaterVapourEvaporation ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::StandAtm_DewPoint_HumidRel(){
    cout << endl << "      StandAtm_DewPoint_HumidRel" << endl;
    double e, E, t_u, height;
//    double beta = 42.0; // in K, COSMO
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            height = get_layer_height(i_mount);
            TempStand.x[i_mount][j][k] = t.x[0][j][k] * t_0 - t_0;
            TempStand.x[0][j][k] = TempStand.x[i_mount][j][k] 
                + 6.5 * height/1000.0; // International Standard Atmosphere (ISA), higher temperatures than COSMO
            for(int i = 0; i < im; i++){
                height = get_layer_height(i);
                t_u = t.x[i][j][k] * t_0; // in K
                if(t_u >= t_0)
                    E = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                else
                    E = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the ice phase at t < 0°C in hPa
                e = c.x[i][j][k] * p_hydro.x[i][j][k]/ep;  // water vapour pressure in hPa
                if(e <= 0) e = 1e-3;
                // linear temperature decay up to tropopause
                // compares to the International Standard Atmosphere (ISA) 
                TempStand.x[i][j][k] = TempStand.x[0][j][k] 
                        - 6.5 * height/1000.0;
//                TempStand.x[i][j][k] = TempStand.x[0][j][k] * sqrt(1.0 
//                    - (2.0 * beta * g * height)
//                    /(R_Air * TempStand.x[0][j][k] * TempStand.x[0][j][k])); // COSMO
                TempDewPoint.x[i][j][k] = (423.86 - 234.175 * log(e))
                    /(log(e) - 18.89);  // in C°, Häckel, Meteorologie, p. 82
                HumidityRel.x[i][j][k] = e/E * 100.0;
                if(HumidityRel.x[i][j][k] > 100.0)
                    HumidityRel.x[i][j][k] = 100.0;
/*
                cout.precision(6);
                cout.setf(ios::fixed);
                if((j == 90) && (k == 180)) cout << endl 
//                if((j == 62) && (k == 87)) cout << endl 
                    << "  StandAtm_DewPoint_HumidRel" << endl
                    << "  i = " << i 
                    << "  height = " << get_layer_height(i) << endl
                    << "  E = " << E
                    << "  e = " << e << endl
                    << "  TempStand = " << TempStand.x[i][j][k]
                    << "  TempDewPoint = " << TempDewPoint.x[i][j][k]
                    << "  HumidityRel = " << HumidityRel.x[i][j][k] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl << endl;
*/
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if(is_land(h, i, j, k)){
                    TempStand.x[i][j][k] = TempStand.x[i_mount][j][k];
                    TempDewPoint.x[i][j][k] = TempDewPoint.x[i_mount][j][k];
                    HumidityRel.x[i][j][k] = HumidityRel.x[i_mount][j][k];
                }
            }
        }
    }
    cout << "      StandAtm_DewPoint_HumidRel ended" << endl;
    return;
}

