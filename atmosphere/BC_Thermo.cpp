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
    // class element for the computation of the radiation and the temperature distribution
    // computation of the local temperature based on short and long wave radiation
    // multi layer radiation model
    temp_tropopause = std::vector<double>(jm, t_tropopause_pole);
    std::vector<double> step(im, 0);
    int j_max = jm-1;
    int j_half = j_max/2;
//    double temp_eff = t_tropopause_pole - t_tropopause_equator;
    for(int j=j_half; j>=0; j--){
        temp_tropopause[j] = t_tropopause_equator;
    }
    for(int j=j_max; j>j_half; j--){
        temp_tropopause[j] = temp_tropopause[j_max-j];
    }
//    AtomUtils::smooth_tropopause(jm, temp_tropopause);
    std::map<float, float> pole_temp_map;  // Stein/Rüdiger/Parish linear pole temperature (Ma) distribution
    load_map_from_file(pole_temperature_file, pole_temp_map); 
    double rad_eff = rad_pole - rad_equator;
    double albedo_co2_eff = albedo_pole - albedo_equator;
    // effective temperature, albedo and emissivity/absorptivity for the two layer model
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im-1; i++){
                if(is_ocean_surface(h, i, j, k)
                    ||is_land_surface(h, i, j, k)){
                    albedo.y[j][k] = albedo_co2_eff * parabola((double)j
                        /(double)j_half) + albedo_pole;
                }
            }
        }
    }
    // absorption/emissivity computation
//    double epsilon_eff_2D = epsilon_pole - epsilon_equator;
    // iteration procedure for the computation of the temperature based on the multi-layer radiation model
    // temperature needs an initial guess which must be corrected by the long wave radiation remaining in the atmosphere
    int i_trop = 0;
    int i_mount = 0;
    for(int j = 0; j < jm; j++){
//        i_trop = m_model->get_tropopause_layer(j);
        i_trop = im-1;
        // on zero level, lateral parabolic distribution
        for(int k = 0; k < km; k++){
            i_mount = i_topography[j][k];
            // shortwave radiation from the sun absobed at the surface 47%
            // in W/m², assumption of parabolic surface radiation at zero level
            radiation_surface.y[j][k] = rad_eff * parabola((double)j
                /(double)j_half) + rad_pole;
            t.x[i_trop][j][k] = t_tropopause_equator;
            for(int i = 0; i < i_trop; i++){
                step[i] = get_layer_height(i+1)
                    - get_layer_height(i);  // in m local atmospheric shell thickness
                // radial parabolic distribution, start on zero level
                // dependency given by Häckel, Meteorologie, p. 205 (law by F. Baur and H. Philips, 1934)
                // epsilon_eff describe the effect in the emissivity computation of other gases like CO2
                // in the original formula this value is 0.594, for reasons of adjustment to the modern atmosphere,
                // this constant became a variable in zonal direction
                // this variable reacts very sensitive and changes the temperature field extremely
                // the second term describes the influence of water vapour only 
                // applicable models described in K. H. Byun and L.-D. Chen: 
                // Total emissivity of CO°2 near earth condition, Journal of Mechanical Science and Technology 27 (10)(2013)3183-3189)
                double P_c = (1e-6 * p_stat.x[i][j][k] 
                             * co2.x[i][j][k] * co2_0/p_0);  // in atm partial pressure of CO2
                double u_c = P_c * step[i] * 100.0; // in atm*cm model by Atwater and Ball
//                u_c = 300.0/0.9869 * P_c * step[i] * 100.0/(t.x[i][j][k] * t_0); 
                        // modell by Yamamoto and Sasamori
                 // influence on the emissivity by carbon dioxcide by the law by Bliss
//                x_c = 10.0 * P_c/(R_co2 * t.x[i][j][k] * t_0) * step[i] * p_0;
                  // coefficient in the Bliss law, 10.0 => conversion from kg/m² to g/cm²
//                eps_co2 = 0.185 * (1.0 - exp (- 50.0 * x_c)); // modell by Bliss
                double eps_co2 = 0.185 * (1.0 - exp(- 0.3919 * pow(u_c, 0.4))); // model by Atwater and Ball
                eps_co2 = 0.5 * eps_co2; // model by Atwater and Ball, 
                eps_co2 = 0.0; 
                // factor .5 fits the Atwater and Ball results to Yamamoto and Sasamori,
                // which fit best to HITRAN band data
//                eps_co2 = 0.0; // no influence of CO2 on the emissivity
//                eps_co2 = 1.5 * eps_co2; // Test, 1.5 is maximum, otherwise nans develop
//                double t_u = t.x[i][j][k] * t_0;
                // COSMO water vapour pressure based on local water vapour in hPa
                double e = c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // in hPa, simplified approach
//                double e = c.x[i][j][k] * p_stat.x[i][j][k]/(c.x[i][j][k] + ep);  // in hPa, complete approach
//                epsilon.x[i][j][k] = eps_co2 + 0.594 + 0.0416 * sqrt(e); // constant  given in Häckel, p. 205 (F. Baur and H. Philips, 1934)
//                epsilon.x[i][j][k] = eps_co2 + 0.52 + 0.065 * sqrt(e); // (Brunt, 1932)
                epsilon.x[i][j][k] = eps_co2 + 0.684 + 0.0056 * e; // (Bignami, 1995, -40 <=> 45°C), best choice
//                epsilon.x[i][j][k] = eps_co2 + 0.8; // (Bignami, 1995, -40 <=> 45°C), best choice
//                epsilon.x[i][j][k] = eps_co2 + (1.0 - 0.26 * exp(- 7.77e-4 
//                    * (t_0 - t_u)); // (Idso and Jackson, 1969, -29 <=> 37°C)
//                epsilon.x[i][j][k] = eps_co2 + 9.2e-6 * t_u * t_u; // (Swinbank, 1963, 2 <=> 29°C)
//                epsilon.x[i][j][k] = eps_co2 + (0.7 + 5.95e-5 * e
//                    * exp(1500.0/t_u)); // (Idso, 1981, -40 <=> 45°C) %%%%%%%%%% good performance, values too high
//                epsilon.x[i][j][k] = eps_co2 + 1.24 * pow((e/t_u),1/7); // (Brutsaert, 1975, -40 <=> 45°C) %%%%%%%%%% good performance, out of function, nans
//                epsilon.x[i][j][k] = 0.98 * pow((e/(t.x[i][j][k] * t_0)), 0.0687); // atmospheric emissivity by Vogel/Bliss, emissivity seems too high
                radiation.x[i][j][k] = sigma * pow(t.x[i][j][k] * t_0, 4.0);
                aux_t.x[i][j][k] = t.x[i][j][k];
/*
                cout.precision(6);
                if((i == 20)&&(j == 90)&&(k == 180))  cout << endl 
                    << "  Multi Layer Radiation Model preparation   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
                    << "   i = " << i << "  j = " << j << "  k = " << k << endl 
                    << "  step = " << step[i]
                    << "  e = " << e
                    << "  P_c = " << P_c
                    << "  u_c = " << u_c << endl
                    << "  eps_co2 = " << eps_co2
                    << "  epsilon = " << epsilon.x[i][j][k] << endl
                    << "  c = " << c.x[i][j][k]
                    << "  cloud = " << cloud.x[i][j][k]
                    << "  ice = " << ice.x[i][j][k] << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0
                    << "  t_tropopause_equator = " << t_tropopause_equator * t_0 - t_0 << endl
                    << "  radiation_mount = " << radiation.x[i_mount][j][k]
                    << "  radiation = " << radiation.x[i][j][k]
                    << "  radiation_trop = " << sigma
                        * pow(t.x[i_trop][j][k] * t_0, 4.0) << endl;
*/
                if(epsilon.x[i][j][k] > 1.)  epsilon.x[i][j][k] = 1.;
            }  // end i
            epsilon.x[i_trop][j][k] = epsilon.x[i_trop-3][j][k] 
                - 3.0 * epsilon.x[i_trop-2][j][k] 
                + 3.0 * epsilon.x[i_trop-1][j][k];
            epsilon_2D.y[j][k] = epsilon.x[0][j][k];
            // inside mountains
            for(int i = i_mount-1; i >= 0; i--){
                if(is_land(h, i, j, k)){
                    epsilon.x[i][j][k] = epsilon.x[i_mount][j][k];
                    radiation.x[i][j][k] = radiation.x[i_mount][j][k];
                }
            }
            // above tropopause
            for(int i = i_trop; i < im; i++){
                epsilon.x[i][j][k] = epsilon.x[i][j][k];
                radiation.x[i][j][k] = sigma * pow(t.x[i][j][k] * t_0, 4.0);
                aux_t.x[i][j][k] = t.x[i][j][k];
            }
        }  // end k
    }  // end j
    // iteration procedure for the computation of the temperature based on the multi-layer radiation model
    // temperature needs an initial guess which must be corrected by the long wave radiation remaining in the atmosphere
        // coefficient formed for the tridiogonal set of equations for the absorption/emission coefficient of the multi-layer radiation model
    for(int j = 0; j < jm; j++){
//        i_trop = get_tropopause_layer(j);
        i_trop = im-1;
        for(int k = 0; k < km; k++){
//            i_mount = i_topography[j][k];
            i_mount = 0;
            std::vector<double> alfa(im, 0.0);
            std::vector<double> beta(im, 0.0);
            std::vector<double> AA(im, 0.0);
            std::vector<std::vector<double> > CC(im, std::vector<double>(im, 0.0));
            double CCC = 0.0, DDD = 0.0;
            // radiation leaving the atmosphere above the tropopause, 
            // emitted by atmophere 49%, Clouds 9% and 12% atmospheric window
            for(int iter_rad = 1;  iter_rad <= 10; iter_rad++){ // iter_rad may be varied
//                t.x[i_trop][j][k] = temp_tropopause[j];
                t.x[i_trop][j][k] = t_tropopause_equator;
                radiation.x[i_trop][j][k] = sigma
                    * pow(t.x[i_trop][j][k] * t_0, 4.0); // radiation at the tropopause
                radiation.x[i_mount][j][k] = sigma
                    * pow(t.x[i_mount][j][k] * t_0, 4.0); // radiation leaving the surface

//                AA[i_mount] = radiation_surface.y[j][k]/radiation.x[i_trop][j][k];// non-dimensional surface radiation
                AA[i_mount] = radiation.x[i_mount][j][k]/radiation.x[i_trop][j][k];// non-dimensional surface radiation
                CC[i_mount][i_mount] = 0.; // no absorption of radiation on the surface by water vapour
                for(int i = i_mount+1; i <= i_trop; i++){
                    aux_t.x[i][j][k] = t.x[i][j][k];
                    AA[i] = AA[i-1] * (1.0 - epsilon.x[i][j][k]); // transmitted radiation from each layer
                    double rad_tmp = sigma * pow(t.x[i][j][k] * t_0, 4.0) 
                        /radiation.x[i_trop][j][k];
                    CC[i][i]= epsilon.x[i][j][k] * rad_tmp; // absorbed radiation in each layer
                    radiation.x[i][j][k] = 
                        (1.0 - epsilon.x[i][j][k]) * rad_tmp; // radiation leaving each layer
                    for(int l = i_mount+1; l <= i_trop; l++){
                        // additional transmitted radiation from layer to layer in radial direction
                        CC[i][l] = CC[i][l-1] * (1.0 - epsilon.x[l][j][k]);
                    }
                }
                // Thomas algorithm to solve the tridiogonal equation system for the solution of the radiation with a recurrence formula
                double aa, bb, cc, dd;
                for(int i = i_mount; i < i_trop; i++){
                    if(i == i_mount){
                        bb = - radiation.x[i][j][k];
                        cc = radiation.x[i+1][j][k];
                        dd = - AA[i];
                        alfa[i] = cc/bb;
                        beta[i] = dd/bb;
                    }else{
                        for(int l = i_mount+1; l <= i-1; l++){
                            CCC = CCC + CC[l][i];
                        }
                        for(int l = i_mount+1; l <= i-2; l++){
                            DDD = DDD + CC[l][i-1];
                        }
                        aa = radiation.x[i-1][j][k];
                        bb = - 2.0 * radiation.x[i][j][k];
                        cc = radiation.x[i+1][j][k];
                        dd = - AA[i-1] + AA[i] + CCC - DDD;
                        alfa[i] = cc/(bb - aa * alfa[i-1]);
                        beta[i] = (dd - aa * beta[i-1]) 
                            /(bb - aa * alfa[i-1]);
                    }
                }
                t.x[i_trop][j][k] = t_tropopause_equator;
                radiation.x[i_trop][j][k] = sigma
                    * pow(t.x[i_trop][j][k] * t_0, 4.0);

                // recurrence formula for the radiation and temperature
                for(int i = i_trop-1; i >= i_mount; i--){
                    // above assumed tropopause constant temperature t_tropopause_equator
                    // Thomas algorithm, recurrence formula
                    radiation.x[i][j][k] = - alfa[i] 
                        * radiation.x[i+1][j][k] + beta[i];
                    t.x[i][j][k] = pow(radiation.x[i][j][k] 
                        /sigma, 0.25)/t_0;
/*
                cout.precision(6);
                 if((i == 20)&&(j == 90)&&(k == 180))  cout << endl 
                    << "  Multi Layer Radiation Model final" << endl
                    << "   iter_rad = " << iter_rad << endl
                    << "   i = " << i << "  j = " << j << "  k = " << k << endl 
                    << "  epsilon = " << epsilon.x[i][j][k] << endl
                    << "  c = " << c.x[i][j][k]
                    << "  cloud = " << cloud.x[i][j][k]
                    << "  ice = " << ice.x[i][j][k] << endl
                    << "  t_tropopause_equator = " << t_tropopause_equator * t_0 - t_0 << endl
                    << "  t = " << t.x[i][j][k] * t_0 - t_0
                    << "  aux_t = " << aux_t.x[i][j][k] * t_0 - t_0
                    << "  t_diff_abs = " << (aux_t.x[i][j][k] - t.x[i][j][k]) * t_0
                    << "  t_diff_rel = " << 1.0 - aux_t.x[i][j][k]/t.x[i][j][k] << endl
                    << "  radiation_mount = " << radiation.x[i_mount][j][k]
                    << "  radiation = " << radiation.x[i][j][k]
                    << "  radiation_trop = " << radiation.x[i_trop][j][k] << endl;
*/
                    t.x[i][j][k] = 0.5 * (t.x[i][j][k] + aux_t.x[i][j][k]); // averaging of temperature values to smooth the iterations
                } // end i
                for(int i = i_trop; i < im; i++){ // above tropopause
//                    t.x[i][j][k] = temp_tropopause[j];
                    t.x[i][j][k] = t_tropopause_equator;
                    radiation.x[i][j][k] = sigma
                        * pow(t.x[i][j][k] * t_0, 4.0);
                } // end i
                radiation.x[i_mount][j][k] = radiation.x[i_mount+3][j][k] 
                    - 3.0 * radiation.x[i_mount+2][j][k] 
                    + 3.0 * radiation.x[i_mount+1][j][k];
                t.x[i_mount][j][k] = t.x[i_mount+3][j][k] 
                    - 3.0 * t.x[i_mount+2][j][k] 
                    + 3.0 * t.x[i_mount+1][j][k];
//              if(fabs(t.x[i][j][k] - aux_t.x[i][j][k]) <= 1.e-3)  break;
            } // end iter_rad
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
}
/*
*
*/
void cAtmosphereModel::PressureDensity(){
// static pressure is understood by a zero velocity field
//    double bet = 42.; // in K given by COSMO
    double t_u = 0.;
    double R_W_R_A = R_WaterVapour/R_Air;
// boundary condition of surface pressure given by gas equation
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            t_u = t.x[0][j][k] * t_0; // in K
            p_stat.x[0][j][k] = 1e-2 * (r_air * R_Air * t_u);  // given in hPa
            r_dry.x[0][j][k] = 1e2 * p_stat.x[0][j][k]/(R_Air * t_u);
            r_humid.x[0][j][k] = 1e2 * p_stat.x[0][j][k]/(R_Air 
                * (1. + (R_W_R_A - 1.) * c.x[0][j][k] 
                - cloud.x[0][j][k] - ice.x[0][j][k]) * t_u); 
/*
    cout.precision(6);
    cout.setf(ios::fixed);
    if((j == 90) && (k == 180)) cout << endl 
        << "  PressureDensity" << endl
        << "  i = " << 0 
        << "  height = " << get_layer_height(0) 
        << "  p_stat_sl = " << p_stat.x[0][j][k]
        << "  p_stat = " << p_stat.x[0][j][k]
        << "  r_dry = " << r_dry.x[0][j][k]
        << "  r_humid = " << r_humid.x[0][j][k]
        << "  t = " << t.x[0][j][k] * t_0 - t_0 << endl;
*/
        }
    }
    double exp_pressure = g/(1.e-2 * gam * R_Air);
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 1; i < im; i++){
                double height = get_layer_height(i);
                double height_0 = get_layer_height(0);
                t_u = t.x[0][j][k] * t_0;
                p_stat.x[i][j][k] = p_stat.x[0][j][k] 
                    * pow(((t_u - gam * height * 1.e-2)/(t_u - gam * height_0 * 1.e-2)), exp_pressure);
                // current air pressure, politropic formula in hPa, Kraus p. 28
                r_dry.x[i][j][k] = 1e2 * p_stat.x[i][j][k]/(R_Air * t_u);
                r_humid.x[i][j][k] = 1e2 * p_stat.x[i][j][k]/(R_Air 
                    * (1. + (R_W_R_A - 1.) * c.x[i][j][k] 
                    - cloud.x[i][j][k] - ice.x[i][j][k]) * t_u); 
                if(is_land(h, i, j, k)){
                    r_dry.x[i][j][k] = r_air;
                    r_humid.x[i][j][k] = r_air;
                }
/*
    cout.precision(6);
    cout.setf(ios::fixed);
    if((j == 90) && (k == 180)) cout << endl 
        << "  PressureDensity" << endl
        << "  i = " << i 
        << "  height = " << get_layer_height(i) 
        << "  p_stat_sl = " << p_stat.x[0][j][k]
        << "  p_stat = " << p_stat.x[i][j][k]
        << "  r_dry = " << r_dry.x[i][j][k]
        << "  r_humid = " << r_humid.x[i][j][k]
        << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl;
*/
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i <= i_mount)){
                    r_dry.x[i][j][k] = r_dry.x[i_mount][j][k];
                    r_humid.x[i][j][k] = r_humid.x[i_mount][j][k];
                }
            }
        }
    }
}
/*
*
*/
void cAtmosphereModel::LatentHeat(){
    float Q_Latent_Ice = 0.; 
    float coeff_Q = cp_l * r_air * t_0; // coefficient for Q_Sensible
    float coeff_lat = 1.e6 * 50.; // diffusion coefficient for water vapour in air (50)
    float coeff_sen = 2.7; // diffusion coefficient for sensible heat transfer in air
    float a, e, step, step_p, step_m;
    for(int j = 0; j < jm; j++){
    // water vapour can condensate/evaporate und sublimate/vaporize
    // water vapour turns to or developes from water or ice
    // latent heat of water vapour
        for(int k = 0; k < km; k++){
            int i_mount = get_surface_layer(j, k);
            for(int i = i_mount+1; i < im-2; i++){
                step_p = (get_layer_height(i+1) - get_layer_height(i));
                step_m = (get_layer_height(i) - get_layer_height(i-1));
/*
                step_p = log((get_layer_height(i+1)/L_atm + 1.)
                    /(get_layer_height(i)/L_atm + 1.)) * L_atm;  // local atmospheric shell thickness
                step_m = log((get_layer_height(i)/L_atm + 1.)
                    /(get_layer_height(i-1)/L_atm + 1.)) * L_atm;  // local atmospheric shell thickness
*/
                step = step_p + step_m;

                e = .01 * c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // water vapour pressure in Pa
                if(i > i_mount){
                    a = e/(R_WaterVapour * t.x[i][j][k] * t_0);  // absolute humidity in kg/m³
                    Q_Latent.x[i][j][k] = - lv * a 
                        * (c.x[i+1][j][k] - c.x[i-1][j][k])/step;
                    Q_Latent_Ice = - ls * a * (ice.x[i+1][j][k] 
                        - ice.x[i-1][j][k])/step;
                    Q_Sensible.x[i][j][k] = - coeff_sen * coeff_Q 
                        * (t.x[i+1][j][k] - t.x[i-1][j][k])/step;  // sensible heat in [W/m2] from energy transport equation
                }
                if(i == i_mount){
                    a = e/(R_WaterVapour * t.x[i][j][k] * t_0);  // absolute humidity in kg/m³
                    Q_Latent.x[i][j][k] = - lv * a 
                        * (- 3. * c.x[i][j][k] + 4. * c.x[i+1][j][k] 
                        - c.x[i+2][j][k])/step;
                    Q_Latent_Ice = - ls * a * (- 3. * ice.x[i][j][k] +
                        4. * ice.x[i+1][j][k] 
                        - ice.x[i+2][j][k])/step;
                    // sensible heat in [W/m2] from energy transport equation
                    Q_Sensible.x[i][j][k] = - coeff_sen * coeff_Q 
                        * (- 3. * t.x[i][j][k] +
                        4. * t.x[i+1][j][k] - t.x[i+2][j][k])/step;
                }
                Q_Latent.x[i][j][k] = 
                    coeff_lat * (Q_Latent.x[i][j][k] + Q_Latent_Ice);
            }
        }
    }
    float c43 = 4./3.;
    float c13 = 1./3.;
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            int i_mount = get_surface_layer(j, k);
            Q_Latent.x[im-1][j][k] = c43 * Q_Latent.x[im-2][j][k] -
                c13 * Q_Latent.x[im-3][j][k];
            Q_Sensible.x[im-1][j][k] = c43 * Q_Sensible.x[im-2][j][k] -
                c13 * Q_Sensible.x[im-3][j][k];
            if(is_land (h, 0, j, k)){  
                Q_Latent.x[0][j][k] = Q_Latent.x[i_mount][j][k];
                Q_Sensible.x[0][j][k] = Q_Sensible.x[i_mount][j][k];
            }else{
                Q_Latent.x[0][j][k] = c43 * Q_Latent.x[1][j][k] -
                    c13 * Q_Latent.x[2][j][k];
                Q_Sensible.x[0][j][k] = c43 * Q_Sensible.x[1][j][k] -
                    c13 * Q_Sensible.x[2][j][k];
            }
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            Q_Latent.x[i][0][k] = c43 * Q_Latent.x[i][1][k] -
                 c13 * Q_Latent.x[i][2][k];
            Q_Latent.x[i][jm-1][k] = c43 * Q_Latent.x[i][jm-2][k] -
                 c13 * Q_Latent.x[i][jm-3][k];
            Q_Sensible.x[i][0][k] = c43 * Q_Sensible.x[i][1][k] -
                 c13 * Q_Sensible.x[i][2][k];
            Q_Sensible.x[i][jm-1][k] = c43 * Q_Sensible.x[i][jm-2][k] -
                 c13 * Q_Sensible.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 1; j < jm-1; j++){
            float b1 = c43 * Q_Latent.x[i][j][1] 
                - c13 * Q_Latent.x[i][j][2];
            float b2 = c43 * Q_Latent.x[i][j][km-2] 
                - c13 * Q_Latent.x[i][j][km-3];
            Q_Latent.x[i][j][0] = Q_Latent.x[i][j][km-1] = (b1+ b2)/2.;
            b1 = c43 * Q_Sensible.x[i][j][1] 
                - c13 * Q_Sensible.x[i][j][2];
            b2 = c43 * Q_Sensible.x[i][j][km-2] 
                - c13 * Q_Sensible.x[i][j][km-3];
            Q_Sensible.x[i][j][0] = Q_Sensible.x[i][j][km-1] = (b1+ b2)/2.;
        }
    }
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
// constant coefficients for the adjustment of cloud water and cloud ice amount vice versa
    int iter_prec_end = 10;
    int iter_prec = 0;
//    double R_A_R_W = R_Air/R_WaterVapour;
    float t_00 = 236.15;
    float t_Celsius_0 = 0.0; // in Celsius = 0 °C
    float q_v_hyp = 0.0;
    float dt_dim = (dr + dthe + dphi)/3.0 * L_atm/u_0; // dimensional time step = 56.6 s
    float coeff_evap_cond = 0.5; // original COSMO
    float t_u = 0.0;
    float T = 0.0;
    float d_t = 0.0;
    float t_Celsius = 0.0;
    float E_Rain = 0.0;
    float E_Ice = 0.0;
    float q_Rain = 0.0;
    float q_Ice = 0.0;
    float q_v_b = 0.0;
    float q_c_b = 0.0;
    float q_i_b = 0.0;
    float CND = 0.0;
    float DEP = 0.0;
    float d_q_v = 0.0;
    float d_q_c = 0.0;
    float d_q_i = 0.0;
// setting water vapour, cloud water and cloud ice into the proper thermodynamic ratio based on the local temperatures
// starting from a guessed parabolic temperature and water vapour distribution in north/south direction
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            for(int i = 0; i < im; i++){
                iter_prec = 1;
                double height = get_layer_height(i);
                t_u = t.x[i][j][k] * t_0; // in K
                t_Celsius = t_u - t_0; // in °C
                if(t_Celsius > t_Celsius_0)  ice.x[i][j][k] = 0.0;
                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] < 0.0)  ice.x[i][j][k] = 0.0;
                q_v_b = c.x[i][j][k];
                q_c_b = cloud.x[i][j][k];
                q_i_b = ice.x[i][j][k];
                T = t_u; // in K
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Ice = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
                q_Rain = ep * E_Rain/(p_stat.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Ice = ep * E_Ice/(p_stat.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_v_hyp = q_v_b;
                if(c.x[i][j][k] >= q_Rain){
                    // §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§
                    for(iter_prec = 1; iter_prec <= iter_prec_end; iter_prec++){ // iter_prec = 2 given by COSMO
                    // condensation == water vapor saturation for cloud water formation, deposition == ice crystal for cloud ice formation
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
                        q_v_b = q_v_b + d_q_v;  // new values
                        q_c_b = q_c_b + d_q_c;
                        q_i_b = q_i_b + d_q_i;
                        if(q_c_b <= 0.0)  q_c_b = 0.0;
                        if(q_i_b <= 0.0)  q_i_b = 0.0;
                        E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        E_Ice = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
                        q_Rain = ep * E_Rain/(p_stat.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                        q_Ice = ep * E_Ice/(p_stat.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                        if((q_c_b > 0.0) && (q_i_b > 0.0))
                            q_v_hyp = (q_c_b * q_Rain + q_i_b * q_Ice) 
                                /(q_c_b + q_i_b);
                        if((q_c_b >= 0.0) && (q_i_b == 0.0))  q_v_hyp = q_Rain;
                        if((q_c_b == 0.0) && (q_i_b > 0.0))  q_v_hyp = q_Ice;
                        if(T >= t_0) q_i_b = 0.0;
                        if(T <= t_00){
                            q_c_b = 0.0;
                            q_i_b = 0.0;
                        }
/*
                        cout.precision(5);
                        cout.setf(ios::fixed);
                        if((j == 60)&&(k == 180))  cout << endl
                            << "  SaturationAdjustment" << endl 
                            << "  i = " << i << "  j = " << j << "  k = " << k << endl
                            << "  iprec = "<< iter_prec << endl
                            << "  height = "<< height << endl
                            << "  CND = " << CND
                            << "  DEP = " << DEP << endl
                            << "  p_stat = " << p_stat.x[i][j][k] 
                            << "  p_dyn = " << p_dyn.x[i][j][k] << endl
                            << "  r_dry = " << r_dry.x[i][j][k] 
                            << "  r_humid = " << r_humid.x[i][j][k] << endl
                            << "  dt = " << d_t 
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
                        if(fabs(q_v_b/q_v_hyp - 1.) <= 1.e-3)  break;
                        else q_v_hyp = 0.5 * (q_v_hyp + q_v_b);  // has smoothing effect
                    } // iter_prec end
                }
                S_c_c.x[i][j][k] = coeff_evap_cond  // negative values == condensation, positive values == evaporation
                    * (cloud.x[i][j][k] - q_c_b)/dt_dim;
                c.x[i][j][k] = q_v_b;  // new values achieved after converged iterations
                cloud.x[i][j][k] = q_c_b;
                ice.x[i][j][k] = q_i_b;
                t.x[i][j][k] = T/t_0;
            } // end i
        } // end j
    } // end k
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            t.x[i][0][k] = c43 * t.x[i][1][k] -
                 c13 * t.x[i][2][k];
            t.x[i][jm-1][k] = c43 * t.x[i][jm-2][k] -
                 c13 * t.x[i][jm-3][k];
            c.x[i][0][k] = c43 * c.x[i][1][k] -
                 c13 * c.x[i][2][k];  
            c.x[i][jm-1][k] = c43 * c.x[i][jm-2][k] -
                 c13 * c.x[i][jm-3][k];
            cloud.x[i][0][k] = c43 * cloud.x[i][1][k] -
                 c13 * cloud.x[i][2][k];
            cloud.x[i][jm-1][k] = c43 * cloud.x[i][jm-2][k] -
                 c13 * cloud.x[i][jm-3][k];
            ice.x[i][0][k] = c43 * ice.x[i][1][k] -
                 c13 * ice.x[i][2][k];
            ice.x[i][jm-1][k] = c43 * ice.x[i][jm-2][k] -
                 c13 * ice.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            t.x[i][j][0] = c43 * t.x[i][j][1] -
                c13 * t.x[i][j][2];
            t.x[i][j][km-1] = c43 * t.x[i][j][km-2] -
                c13 * t.x[i][j][km-3];
            t.x[i][j][0] = t.x[i][j][km-1] =
               (t.x[i][j][0] + t.x[i][j][km-1])/2.;
            c.x[i][j][0] = c43 * c.x[i][j][1] -
                c13 * c.x[i][j][2];
            c.x[i][j][km-1] = c43 * c.x[i][j][km-2] -
                c13 * c.x[i][j][km-3];
            c.x[i][j][0] = c.x[i][j][km-1] =
               (c.x[i][j][0] + c.x[i][j][km-1])/2.;
            cloud.x[i][j][0] = cloud.x[i][j][km-1] =
               (cloud.x[i][j][0] + cloud.x[i][j][km-1])/2.;
            cloud.x[i][j][0] = c43 * cloud.x[i][j][1] -
                c13 * cloud.x[i][j][2];
            cloud.x[i][j][0] = cloud.x[i][j][km-1] =
               (cloud.x[i][j][0] + cloud.x[i][j][km-1])/2.;
            ice.x[i][j][km-1] = c43 * ice.x[i][j][km-2] -
                c13 * ice.x[i][j][km-3];
            ice.x[i][j][0] = ice.x[i][j][km-1] =
               (ice.x[i][j][0] + ice.x[i][j][km-1])/2.;
            ice.x[i][j][0] = ice.x[i][j][km-1] =
               (ice.x[i][j][0] + ice.x[i][j][km-1])/2.;
        }
    }
/*
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.0)      c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] < 0.0)    ice.x[i][j][k] = 0.0;
                if(is_land(h, i, j, k)){
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k] = 0.0;
                    S_c_c.x[i][j][k] = 0.0;
                }
            }
        }
    }
*/
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
/*
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i <= i_mount)){
                    c.x[i][j][k] = c.x[i_mount][j][k];
                    cloud.x[i][j][k] = cloud.x[i_mount][j][k];
                    ice.x[i][j][k] = ice.x[i_mount][j][k];
//                    t.x[i][j][k] = t.x[i_mount][j][k];
                }
            }
*/
            if(is_land(h, 0, j, k)){
                c.x[0][j][k] = c.x[i_mount][j][k];
                cloud.x[0][j][k] = cloud.x[i_mount][j][k];
                ice.x[0][j][k] = ice.x[i_mount][j][k];
                t.x[0][j][k] = t.x[i_mount][j][k];
            }
        }
    }
    return;
}
/*
*
*/
// Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
// resulting the precipitation distribution formed of rain and snow
void cAtmosphereModel::TwoCategoryIceScheme(){   
    // constant coefficients for the transport of cloud water and cloud ice amount vice versa, 
    // rain and snow in the parameterization procedures
//    double R_A_R_W = R_Air/R_WaterVapour;
    float N_i_0 = 1.0e2,  // in m-3
          m_i_0 = 1.0e-12,  // in kg
          m_i_max = 1.0e-9,  // in kg
          m_s_0 = 3.0e-9,  // in kg
          c_i_dep = 1.3e-5,  // in m3/(kg*s)
          c_c_au = 4.0e-4,  // in 1/s
          c_i_au = 1.0e-3,  // in 1/s
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
          t_hn = 236.15,  // in K    -40 °C
          t_r_frz = 271.15;  // in K    -2 °C
    float dt_snow_dim = 417.,  // dt_snow_dim is the time  in 417 s to pass dr = 400 m, 400 m/417 s = .96 m/s fallout velocity
          dt_rain_dim = 250.0;  // dt_rain_dim is the time  in 250 s to pass dr = 400 m, 400 m/250 s = 1.6 m/s fallout velocity
    float q_Rain = 0.0,
          E_Rain = 0.0;
    float q_Ice = 0.0,
          E_Ice = 0.0;
    double m_i = m_i_max;  
    float N_i, S_nuc, S_c_frz, S_i_dep=0, S_c_au, S_i_au, S_d_au, 
        S_ac, S_rim, S_shed;
    float S_agg, S_i_cri, S_r_cri, S_ev, S_s_dep, S_i_melt, S_s_melt, S_r_frz;
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
            }
        }
    }
    if(true){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                P_rain.x[im-1][j][k] = 0.0;
                P_snow.x[im-1][j][k] = 0.0;
                S_r.x[im-1][j][k] = 0.0;
                S_s.x[im-1][j][k] = 0.0;
                for(int i = im-2; i >= 0; i--){
                    float t_u = t.x[i][j][k] * t_0;
                    float t_Celsius = t_u - t_0;
                    step[i] = get_layer_height(i+1) - get_layer_height(i);  // local atmospheric shell thickness
                    E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                    q_Rain = ep * E_Rain/(p_stat.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                    E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                    q_Ice = ep * E_Ice/(p_stat.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                    dt_rain_dim = step[i]/1.6; 
                    // adjusted rain fall time step by fixed velocities == 1.6 m/s by variable local step size
                    dt_snow_dim = step[i]/0.96; 
                    // ice and snow average size 
                    if(!(t_u > t_0)){  
                        N_i = N_i_0 * exp(0.2 * (t_0 - t_u));
                        m_i = r_humid.x[i][j][k] * ice.x[i][j][k]/N_i;
                        if(m_i > m_i_max) { m_i = m_i_max; }
                        if(m_i < m_i_0) { m_i = m_i_0; }
                    }
                    // nucleation and depositional growth of cloud ice
                    if(ice.x[i][j][k] == 0.0){
                        if(((t_u < t_d) && (c.x[i][j][k] >= q_Ice))
                            || (((t_d <= t_u) && (t_u <= t_nuc)) 
                            && (c.x[i][j][k] >= q_Rain)))
                            S_nuc = m_i_0/(r_humid.x[i][j][k] * dt_snow_dim) * N_i;  // nucleation of cloud ice, < I >
                    }
                    else  S_nuc = 0.0;
                    if((t_u < t_hn) && (cloud.x[i][j][k] > 0.0))
                        S_c_frz = cloud.x[i][j][k]/dt_rain_dim;  //nucleation of cloud ice due to freezing of cloud water, < II >
                    else  S_c_frz = 0.0;
                    if(!(t_Celsius > 0.0)){  //temperature <= 0
                        if(c.x[i][j][k] > q_Ice){  // supersaturation
                            S_i_dep = c_i_dep * N_i * pow(m_i, (1./3.)) 
                                * (c.x[i][j][k] - q_Ice);  // supersaturation, < III >
                        }else if(- ice.x[i][j][k]  > c.x[i][j][k] - q_Ice) 
                            S_i_dep = - ice.x[i][j][k]/dt_snow_dim; // subsaturation, < III >
                    }
                    else  S_i_dep = 0.0; //temperature > 0
                    // autoconversion processes
                    if((t_u >= t_0) && (cloud.x[i][j][k] > 0.0))
                        S_c_au = c_c_au * cloud.x[i][j][k];  // cloud water to rain, cloud droplet collection, < IV >
                    else  S_c_au = 0.0;
                    if((t_u <= t_0) && (ice.x[i][j][k] > 0.0))
                        S_i_au = c_i_au * ice.x[i][j][k];  // cloud ice to snow, cloud ice crystal aggregation, < V >
                    else  S_i_au = 0.0;
                    if(t_u <= t_0)
                        S_d_au = S_i_dep/(1.5 * (pow((m_s_0/m_i), 
                            (2./3.)) - 1.0));  // autoconversion due to depositional growth of cloud ice, < VI >
                    else  S_d_au = 0.0;
                    // collection mechanism
                    if(t_u > t_0)  
                        S_ac = c_ac * cloud.x[i][j][k] 
                            * pow(P_rain.x[i][j][k],(7.0/9.0));
                        // accreation rate from depletion of cloud water due to collection by all rain drops, < VII >
                    else  S_ac = 0.0;  // accreation rate from depletion of cloud water due to collection by all rain drops
                    if(t_u < t_0)  
                        S_rim = c_rim * cloud.x[i][j][k] * P_snow.x[i][j][k];
                    else  S_rim = 0.0;  // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >
                                       // by falling snow particles
                    if(t_u >= t_0)  
                        S_shed = c_rim * cloud.x[i][j][k] * P_snow.x[i][j][k];
                    else  S_shed = 0.0;  // rate of water shed by melting wet snow particles, < IX >
                                        // collecting cloud droplets to produce rain
                    if(t_u <= t_0){
                        S_agg = c_agg * ice.x[i][j][k] * P_snow.x[i][j][k];  // collection of cloud ice by snow particles, < X >
                        S_i_cri = c_i_cri * ice.x[i][j][k] 
                            * pow(P_rain.x[i][j][k],(7.0/9.0));
                            // decrease in cloud ice mass due to collision/coalescense interaction with raindrops, < XI >
                        S_r_cri = c_r_cri * ice.x[i][j][k]/m_i 
                            * pow(P_rain.x[i][j][k],(13.0/9.0));
                            // decrease of rainwater due to freezing resulting from collection of ice crystals, < XII >
                    }else{
                        S_agg = 0.0;
                        S_i_cri = 0.0;
                        S_r_cri = 0.0;
                    }
                    // diffusional growth of rain and snow
                    if(t_u < t_0){  // temperature below zero
                        S_s_dep = c_s_dep * (1. + b_s_dep 
                            * pow(P_snow.x[i][j][k],(5.0/26.0))) * 
                            (c.x[i][j][k] - q_Ice) 
                            * pow(P_snow.x[i][j][k],(8.0/13.0));
                         // deposition/sublimation of snow, < XIV >
                        S_ev = 0.0;
                    }else{  
                        S_ev = a_ev * (1. + b_ev 
                            * pow(P_rain.x[i][j][k],(1.0/6.0))) *
                            (q_Rain - c.x[i][j][k]) 
                            * pow(P_rain.x[i][j][k],(4.0/9.0));
                        // evaporation of rain due to water vapour diffusion, < XIII >
                        S_s_dep = 0.0;
                    }
                    // melting and freezing
                    S_i_melt = 0.0;
                    S_s_melt = 0.0;
                    if(t_u > t_0){ // temperature above zero
                        if(ice.x[i][j][k] > 0.0)
                            S_i_melt = ice.x[i][j][k]/dt_snow_dim; // cloud ice particles melting to cloud water, < XV >
//                            float p_t_in = p_stat.x[i][j][k];
                        float p_t_in = p_stat.x[i][j][k];
                        float E_Rain_t_in = hp * exp_func(t_0, 17.2694, 35.86);
                            // saturation water vapour pressure for the water phase at t = 0°C in hPa
                        float q_Rain_t_in = ep * E_Rain_t_in/(p_t_in 
                            - E_Rain_t_in);
                            // water vapour amount at saturation with water formation in kg/kg
                        S_s_melt = c_s_melt * (1. + b_s_melt 
                            * pow(P_snow.x[i][j][k],(5.0/26.0))) *
                            ((t_u - t_0) + a_s_melt * (c.x[i][j][k] 
                            - q_Rain_t_in)) * pow(P_snow.x[i][j][k], 
                            (8./13.));  // melting rate of snow to form rain, < XVI >
                    }
                    if(t_r_frz - t_u > 0.0)
                        S_r_frz = c_r_frz * pow((t_r_frz - t_u), 
                            (3./2.)) * pow(P_rain.x[i][j][k],(3.0/2.0));
                            // immersion freezing and contact nucleation, < XVII >
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
                    if((is_land(h,i,j,k)) && (is_land(h,i+1,j,k))){
                        S_c_c.x[i][j][k] = 0.0;
                        S_v.x[i][j][k] = 0.0;
                        S_c.x[i][j][k] = 0.0;
                        S_i.x[i][j][k] = 0.0;
                        S_r.x[i][j][k] = 0.0;
                        S_s.x[i][j][k] = 0.0;
                    }
                    // rain and snow integration
//                    int i_mount = i_topography[j][k];
                    if(P_rain.x[i+1][j][k] < 0.0)  P_rain.x[i+1][j][k] = 0.0;
                    if(P_snow.x[i+1][j][k] < 0.0)  P_snow.x[i+1][j][k] = 0.0;
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * 0.5 * (S_r.x[i+1][j][k] 
                        + S_r.x[i][j][k]) * step[i];  // in kg/(m2 * s) == mm/s 
                    if(t_Celsius < 0.0)  P_snow.x[i][j][k] = P_snow.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * 0.5 * (S_s.x[i+1][j][k] 
                        + S_s.x[i][j][k]) * step[i];  // in kg/(m2 * s) == mm/s 
                    else  P_snow.x[i][j][k] = 0.0;
/*
                    cout.precision(6);
                    cout.setf(ios::fixed);
                    if((j == 124) && (k == 87)) cout << endl 
                        << "  TwoCategoryIceScheme ----- SINKS and SOURCES" << endl 
                        << "  i = " << i << "  j = " << j << "  k = " << k << endl
                        << "  step = " << step[i] << endl
                        << "  Sev = " << S_ev * 1e3 
                        << "  Sidep = " << S_i_dep * 1e3 
                        << "  Ssdep = " << S_s_dep * 1e3 
                        << "  Snuc = " << S_nuc * 1e3 
                        << "  Scau = " << S_c_au * 1e3 
                        << "  Sagg = " << S_agg * 1e3 << endl 
                        << "  Scfrz = " << S_c_frz * 1e3 
                        << "  Simelt = " << S_i_melt * 1e3 
                        << "  Srim = " << S_rim * 1e3 
                        << "  Sshed = " << S_shed * 1e3 
                        << "  Siau = " << S_i_au * 1e3
                        << "  Sdau = " << S_d_au * 1e3
                        << "  Sac = " << S_ac * 1e3 << endl;
*/
/*
                    cout.precision(6);
                    cout.setf(ios::fixed);
                    if((j == 124) && (k == 87)) cout << endl 
                        << "  TwoCategoryIceScheme ----- RAIN and SNOW PRODUCTION" << endl
                        << "  i = " << i << "  j = " << j << "  k = " << k << endl
                        << "  step = " << step[i] << endl
                        << "  S_c_c = " << S_c_c.x[i][j][k] * 1e3 
                        << "  S_v = " << S_v.x[i][j][k] * 1e3 
                        << "  S_c = " << S_c.x[i][j][k] * 1e3 
                        << "  S_i = " << S_i.x[i][j][k] * 1e3 
                        << "  S_r = " << S_r.x[i][j][k] * 1e3 
                        << "  S_s = " << S_s.x[i][j][k] * 1e3 << endl
                        << "  P_rain = " << P_rain.x[i][j][k] * 1.0 
                        << "  P_snow = " << P_snow.x[i][j][k] * 1.0
                        << "  P_conv = " << P_conv.x[i][j][k] * 1.0 << endl
                        << "  c = " << c.x[i][j][k] * 1e3
                        << "  cl = " << cloud.x[i][j][k] * 1e3
                        << "  ice = " << ice.x[i][j][k] * 1e3 << endl
                        << "  t = " << t.x[i][j][k] * t_0 - t_0 << endl
                        << "  r_dry = " << r_dry.x[i][j][k] 
                        << "  r_humid = " << r_humid.x[i][j][k] << endl
                        << "  p_stat = " << p_stat.x[i][j][k]
                        << "  p_dyn = " << p_dyn.x[i][j][k] << endl;
*/
                }  // end i RainSnow
            }  // end j
        }  // end k
    }  // end if true
    fft_gaussian_filter_3d(P_rain,1);
    fft_gaussian_filter_3d(P_snow,1);
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i <= i_mount)){
                    P_rain.x[i][j][k] = P_rain.x[i_mount][j][k];
                    P_snow.x[i][j][k] = P_snow.x[i_mount][j][k];
                }
            }
        }
    }
    return;
}
/*
*
*/
void cAtmosphereModel::ValueLimitationAtm(){
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            if(Precipitation.y[j][k] >= 25.0)  Precipitation.y[j][k] = 25.0;
            if(Precipitation.y[j][k] <= 0.0)  Precipitation.y[j][k] = 0.0;
            for(int i = 0; i < im; i++){
                if(u.x[i][j][k] >= 0.106)  u.x[i][j][k] = 0.106;
                if(u.x[i][j][k] <= - 0.106)  u.x[i][j][k] = - 0.106;
                if(v.x[i][j][k] >= 0.5)  v.x[i][j][k] = 0.5;
                if(v.x[i][j][k] <= - 0.5)  v.x[i][j][k] = - 0.5;
                if(w.x[i][j][k] >= 6.25)  w.x[i][j][k] = 6.25;
                if(w.x[i][j][k] <= - 2.0)  w.x[i][j][k] = - 2.0;
                if(t.x[i][j][k] >= 1.165)  t.x[i][j][k] = 1.165;  // == 45.0 °C
                if(t.x[i][j][k] <= - 0.78)  t.x[i][j][k] = - 0.78;  // == 59.82 °C
                if(c.x[i][j][k] >= 0.04)  c.x[i][j][k] = 0.04;
                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] >= 0.02)  cloud.x[i][j][k] = 0.02;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] >= 0.01)  ice.x[i][j][k] = 0.01;
                if(ice.x[i][j][k] < 0.0)  ice.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] >= 10.0)  P_rain.x[i][j][k] = 10.0;
                if(P_snow.x[i][j][k] >= 1.0)  P_snow.x[i][j][k] = 1.0;
                if(P_conv.x[i][j][k] >= 10.0)  P_conv.x[i][j][k] = 10.0;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                if(P_conv.x[i][j][k] < 0.0)  P_conv.x[i][j][k] = 0.0;
                if(co2.x[i][j][k] >= 5.36)  co2.x[i][j][k] = 5.36;
                if(co2.x[i][j][k] <= 1.0)  co2.x[i][j][k] = 1.0;
            }
        }
    }
}
/*
*
*/
void cAtmosphereModel::MoistConvection(){
// collection of coefficients for phase transformation
    int i_base = 0;
    int ii_base_beg = 1;
    int i_lfs = 0;
//    int i_tropopause = 0;
    double R_W_R_A = R_WaterVapour/R_Air;
    double t_00 = 236.15; // = -37°C
    double bet = 42.; // in K given by COSMO
    double b_u = .3;
    double alf_1 = 5.e-4; // in 1/s
    double alf_2 = .011; // in 1/s
    double p_ps = 50.;
    double C_p = .05; // in ./.
    double bet_p = 2e-3;  // in 1/s, from the COSMO paper
    double gam_d = 0.3;  // from original paper by Tiedtke, gam_d = - 0.2
    double eps_u = 1.e-4;  // in 1/m
    double eps_d = 2.e-4;
    double del_u = eps_u;
    double del_d = eps_d;
    double delta_i_c = 0.0;
    double K_p = 0.0;  // in 1/s
    double HumidityRelative = 0.0, L_latent = 0.0;
    double q_buoy_base = 0.0, q_buoy_lfs = 0.0, cloud_buoy = 0.0;
    double E_Rain = 0.0, q_Rain = 0.0;
    double E_Ice = 0.0;
    double E_Rain_add = 0.0, q_Rain_add = 0.0;
    double E_Ice_add = 0.0, q_Ice = 0.0, q_Ice_add = 0.0;
    double r_humid_add = 0.0, r_dry_add = 0.0, p_stat_add = 0.0;
    double t_u = 0.0, t_u_add = 0.0, height = 0.0, height_base = 0.0;
    double rm = 0.0, exp_rm = 0.0, sinthe = 0.0, rmsinthe = 0.0;
    double M_u_denom = 0.0, M_d_denom = 0.0, M_u_num = 0.0, M_d_num = 0.0;
    double dcdr = 0.0, dcdthe = 0.0, dcdphi = 0.0;
//    double t_add_u = 1.0; // == + 0.0 °C temperature increase for a lifting air parcel, given by COSMO for Tiedtke-Scheme 1989
    double t_add_u = 1.001; // == + 0.3 °C temperature increase for a lifting air parcel, given by COSMO for Tiedtke-Scheme 1989
//    double t_add_u = 1.01; // == + 2.73 °C temperature increase for a lifting air parcel, given by COSMO for Tiedtke-Scheme 1989
    double area = 0.0;
    double exp_pressure = g/(1.0e-2 * gam * R_Air);
    std::vector<double> step(im, 0);
    std::vector<std::vector<double> > M_u_Base(jm, std::vector<double> (km, 0));
    std::vector<std::vector<double> > M_d_LFS(jm, std::vector<double> (km, 0));
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            int i_mount = i_topography[j][k];
            sinthe = sin(the.z[j]);
            for(int i = im-1; i >= i_mount ; i--){
                rm = rad.z[i] * get_layer_height(i);
                rmsinthe = rm * sinthe;
                exp_rm = 1.0/(rm + 1.0);
/*
                dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k])/(2. * dr) * exp_rm;
                if(i < im-2){
                    if((is_land(h, i, j, k))&&((is_air(h, i+1, j, k))
                        &&(is_air(h, i+2, j, k)))){
                        dcdr = (- 3. * c.x[i][j][k] + 4. * c.x[i+1][j][k] 
                            - c.x[i+2][j][k])/(2. * dr ) * exp_rm;
                    }
                }
*/
                if(c.x[i][j][k] == 0.0) c.x[i][j][k] = c_tropopause;
                dcdr = 0.0; 
                dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k])
                    /(2.0 * dthe * rm);
                dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1])
                    /(2.0 * dphi * rmsinthe);
/*
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
*/
/*
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
*/
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
    fft_gaussian_filter_3d(E_u,1);
// parameterisation of moist convection
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
//            costhe = cos(the.z[j]);
//            if(j <= 90) area = 111.11 * 111.11 * (1.- costhe) * 1e6; // in m², simplified formula
//            else area = - 111.11 * 111.11 * (1.+ costhe) * 1e6;
            area = 1.; // in m²
// boundary conditions on the base of clouds (level of free convection)
            if(i_topography[j][k] > 0) 
                ii_base_beg = i_topography[j][k] + 1;
            else ii_base_beg = 1;
            for(int ii_base = ii_base_beg; ii_base < im-1; ii_base++){
//                int i_mount = i_topography[j][k];
                height = get_layer_height(ii_base);
                t_u = t.x[ii_base][j][k] * t_0; // in K
                t_u_add = t_add_u * t.x[ii_base][j][k] * t_0; // in K
                p_stat_add = p_stat.x[0][j][k] 
                    * pow(((t_u_add - gam * height * 1.e-2)
                    /(t_u_add - gam * height * 1.e-2)), exp_pressure);
                // current air pressure, politropic formula in hPa, Kraus p. 28
                r_dry_add = 1e2 * p_stat_add/(R_Air * t_u_add);
                r_humid_add = 1e2 * p_stat_add/(R_Air 
                    * (1. + (R_W_R_A - 1.) * c.x[ii_base][j][k] 
                    - cloud.x[ii_base][j][k] - ice.x[ii_base][j][k]) * t_u_add); 
                q_v_u.x[ii_base][j][k] = c.x[ii_base][j][k];
                q_c_u.x[ii_base][j][k] = cloud.x[ii_base][j][k];
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Rain_add = hp * exp_func(t_u_add, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_stat.x[ii_base][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Rain_add = ep * E_Rain_add/(p_stat.x[ii_base][j][k] - E_Rain_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                HumidityRelative = r_dry.x[ii_base][j][k]
                    /r_humid.x[ii_base][j][k] * 100.0;
                cloud_buoy = r_dry.x[ii_base][j][k] 
                    - r_humid.x[ii_base][j][k];
                q_buoy_base = r_humid.x[ii_base][j][k] - r_humid_add;
                if((HumidityRelative >= 90.)&&(q_buoy_base <= 0.0)){ // search for cloud base
                    i_Base.y[j][k] = (double)ii_base;
                    i_base = ii_base;
                    M_u.x[i_base][j][k] = r_humid.x[i_base][j][k] 
                        * u.x[i_base][j][k] * u_0;  // in kg/(m²s) [== mm/s] updraft at cloud base 
                    M_d_LFS[j][k] = gam_d * M_u.x[i_base][j][k];  // in kg/(m²s) [== mm/s] downdraft at cloud top 
//                    if((i_base >= i_mount)&&(i_mount != 0))
//                        M_d_LFS[j][k] = 1e-1 * gam_d * M_u.x[i_base][j][k];
                    s.x[i_base][j][k] = (cp_l * t_u + g * height)/s_0;
                    q_v_u.x[i_base][j][k] = c.x[i_base][j][k];
                    q_c_u.x[i_base][j][k] = cloud.x[i_base][j][k];
                    u_u.x[i_base][j][k] = u.x[i_base][j][k] 
                        + M_u.x[i_base][j][k]/(r_humid.x[i_base][j][k] * area);  // in m/s
                    v_u.x[i_base][j][k] = v.x[i_base][j][k];
                    w_u.x[i_base][j][k] = w.x[i_base][j][k];
                    s_u.x[i_base][j][k] = s.x[i_base][j][k];
                    g_p.x[i_base][j][k] = 0.0;
                    c_u.x[i_base][j][k] = 0.0;
                    e_l.x[i_base][j][k] = 0.0;
                    e_p.x[i_base][j][k] = 0.0;
/*
                    cout.precision(6);
                    if((j == 90) &&(k == 180))  cout << endl 
                        << "  MoistConvection   Cloud Base" << endl 
                        << "  i_base UPDRAFT" << endl 
                        << "  j = " << j << "  k = " << k << endl 
                        << "  ii_base = " << ii_base 
                        << "  i_base = " << (int)i_Base.y[j][k]
                        << "  i_lfs = " << (int)i_LFS.y[j][k]
                        << "  area = " << area << endl
                        << "  i_tropo = " << i_topography[j][k] << endl
                        << "  s = " << s.x[i_base][j][k] 
                        << "  s_u = " << s_u.x[i_base][j][k] << endl
                        << "  E_Rain = " << E_Rain 
                        << "  E_Rain_add = " << E_Rain_add 
                        << "  q_Rain = " << q_Rain 
                        << "  q_Rain_add = " << q_Rain_add << endl
                        << "  cloud_buoy = " << cloud_buoy 
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
                        << "  p_dyn = " << p_dyn.x[i_base][j][k] 
                        << "  p_stat = " << p_stat.x[i_base][j][k] 
                        << "  p_stat_add = " << p_stat_add 
                        << "  M_u = " << M_u.x[i_base][j][k] 
                        << "  M_d_LFS = " << M_d_LFS[j][k] << endl;
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
                q_Rain = ep * E_Rain/(p_stat.x[ii_lfs][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_v_d.x[ii_lfs][j][k] = c.x[ii_lfs][j][k];
                double cloud_buoy = r_dry.x[i_lfs][j][k] - r_humid.x[i_lfs][j][k];
                q_buoy_lfs = .5 * (q_c_u.x[ii_lfs][j][k] + q_Rain);
                if((q_buoy_lfs <= c.x[ii_lfs][j][k])&&(t_u >= t_00)){ // search for cloud top
                    i_LFS.y[j][k] = (double)ii_lfs;
                    i_lfs = ii_lfs;
                    if((t.x[i_lfs][j][k] * t_0) >= t_0) 
                        L_latent = lv;
                    else L_latent = ls;
                    s.x[i_lfs][j][k] = (cp_l * t_u + g * height)/s_0;
                    q_v_d.x[i_lfs][j][k] = c.x[i_lfs][j][k];
                    M_d.x[i_lfs][j][k] = M_d_LFS[j][k];  // in kg/(m²s) [== mm/s] downdraft at cloud top
                    u_d.x[i_lfs][j][k] = u.x[i_lfs][j][k] 
                        + M_d.x[i_lfs][j][k]/(r_humid.x[i_lfs][j][k] * area); 
                    v_d.x[i_lfs][j][k] = v.x[i_lfs][j][k];
                    w_d.x[i_lfs][j][k] = w.x[i_lfs][j][k];
                    s_d.x[i_lfs][j][k] = s.x[i_lfs][j][k];
                    e_p.x[i_lfs][j][k] = 0.0;
                    e_d.x[i_lfs][j][k] = q_v_d.x[i_lfs][j][k] - q_Rain;
/*
                    cout.precision(6);
                    if((j == 90) &&(k == 180))  cout << endl 
                        << "  MoistConvection   Cloud Top"  << endl 
                        << "  i_lfs DOWNDRAFT" << endl
                        << "  j = " << j<< "  k = " << k << endl 
                        << "  ii_lfs = " << ii_lfs 
                        << "  i_base = " << (int)i_Base.y[j][k] 
                        << "  i_lfs = " << (int)i_LFS.y[j][k] << endl
                        << "  s = " << s.x[i_lfs][j][k] 
                        << "  s_d = " << s_d.x[i_lfs][j][k] << endl
                        << "  E_Rain = " << E_Rain 
                        << "  q_Rain = " << q_Rain << endl
                        << "  cloud_buoy = " << cloud_buoy 
                        << "  q_buoy_base = " << q_buoy_base 
                        << "  q_buoy_lfs = " << q_buoy_lfs << endl
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
                        << "  p_dyn = " << p_dyn.x[i_lfs][j][k] 
                        << "  p_stat = " << p_stat.x[i_lfs][j][k] 
                        << "  M_d = " << M_d.x[i_lfs][j][k] 
                        << "  M_d_LFS = " << M_d_LFS[j][k] << endl;
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
                s_u.x[i][j][k] = s.x[i][j][k];
                s_d.x[i][j][k] = s.x[i][j][k];
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
                s_u.x[i][j][k] = s.x[i][j][k];
                s_d.x[i][j][k] = s.x[i][j][k];
                q_v_d.x[i][j][k] = 0.0;
                u_d.x[i][j][k] = 0.0;
                v_d.x[i][j][k] = 0.0;
                w_d.x[i][j][k] = 0.0;
            } // zero moist values above i_lfs
        }  // end j for cloud base and cloud top
    }  // end k for cloud base and cloud top
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
// CONVECTIVE UPDRAFT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            i_base = (int)i_Base.y[j][k];
            i_lfs = (int)i_LFS.y[j][k];
//            i_tropopause = get_tropopause_layer(j);
            double i_tropopause = im-1;
            for(int i = i_base; i < im-1; i++){
//            for(int i = i_base; i <= i_lfs; i++){
                step[i] = get_layer_height(i+1) - get_layer_height(i);  // in m local atmospheric shell thickness
                double cloud_buoy = r_dry.x[i][j][k] - r_humid.x[i][j][k];
                HumidityRelative = r_dry.x[i][j][k]/r_humid.x[i][j][k] * 100.;
                t_u = t.x[i][j][k] * t_0; // in K
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                t_u_add = t_add_u * t.x[i][j][k] * t_0; // in K
                p_stat_add = p_stat.x[0][j][k] * exp(- t_u_add/bet 
                    * (1. - sqrt(1. - 2. * bet * g * height
                    /(R_Air * t_u_add * t_u_add))));
                E_Rain_add = hp * exp_func(t_u_add, 17.2694, 35.86);
                E_Ice_add = hp * exp_func(t_u, 21.8746, 7.66);
                q_Rain = ep * E_Rain/(p_stat.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Rain_add = ep * E_Rain_add/(p_stat.x[i][j][k] - E_Rain_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Ice = ep * E_Ice/(p_stat.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Ice_add = ep * E_Ice_add/(p_stat.x[i][j][k] - E_Ice_add); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                r_dry_add = 1e2 * p_stat_add/(R_Air * t_u_add);
                r_humid_add = 1e2 * p_stat_add/(R_Air * (1. 
                    + (R_W_R_A - 1.) * q_v_u.x[i][j][k] 
                    - q_c_u.x[i][j][k]) * t_u_add); 
                q_buoy_base = r_humid.x[i][j][k] - r_humid_add;
                if(t_u_add >= t_0) L_latent = lv;
                else L_latent = ls;
                s.x[i][j][k] = (cp_l * t_u + g * height)/s_0;
                eps_u = 1.e-4;
                del_u = 1.e-4;
                if(t_u_add < t_00){
                    eps_u = 0.0;
                    del_u = 0.0;
                }
// water vapour in the updraft
                if(q_v_u.x[i][j][k] >= q_Rain_add){
                    q_c_u.x[i][j][k] = q_v_u.x[i][j][k] - q_Rain_add;
                    if(q_c_u.x[i][j][k] <= 0.0) q_c_u.x[i][j][k] = 0.0;
                    if(q_c_u.x[i][j][k] >= 0.0) 
                        q_v_u.x[i][j][k] = q_Rain_add;
                    if(q_v_u.x[i][j][k] < 0.0) q_v_u.x[i][j][k] = q_Rain_add;
                }else  q_c_u.x[i][j][k] = 0.0;
// condensation within the updraft
                if((t_u_add >= t_0)&&(q_c_u.x[i][j][k] >= 0.0))
                    c_u.x[i][j][k] = q_v_u.x[i][j][k] - q_Rain_add;
                else  c_u.x[i][j][k] = q_v_u.x[i][j][k] - q_Ice_add;
                if(c_u.x[i][j][k] < 0.0) c_u.x[i][j][k] = 0.0;
// specification of entrainment in the updraft
                E_u.x[i][j][k] = E_u.x[i][j][k] 
                    + eps_u * M_u.x[i][j][k]; // in kg/(m³s)
// specification of detrainment in the updraft
                if(i == i_lfs)
                    D_u.x[i][j][k] = (1. - b_u) 
                        * M_u.x[i][j][k]/step[i]
                        + del_u * M_u.x[i][j][k];
                if(i == i_lfs+1)
                    D_u.x[i][j][k] = b_u 
                        * M_u.x[i][j][k]/step[i]
                        + del_u * M_u.x[i][j][k];  // in kg/(m³s)
// evaporation of cloud water in the updraft
                e_l.x[i][j][k] = D_u.x[i][j][k]/r_humid.x[i][j][k]
                    * q_c_u.x[i][j][k];  // in kg/(kg*s)
// formation of precipitation within the updraft
                if(is_land(h, 0, j, k)) delta_i_c = 3000.0; // over land in m
                if(is_air(h, 0, j, k))  delta_i_c = 1500.0; // over ocean in m
                height = get_layer_height(i); // in m
                height_base = get_layer_height(i_base); // in m
                if(height > height_base + delta_i_c)  K_p = bet_p;  // in 1/s
                else  K_p = 0.0;
                if(t_u_add >= t_00)
                    g_p.x[i][j][k] = K_p * q_c_u.x[i][j][k];  // in kg/(kg*s)
                else  g_p.x[i][j][k] = 0.0;
                if(g_p.x[i][j][k] < 0.0)  g_p.x[i][j][k] = 0.0;
// mass flux within the updraft
                if(is_land(h, i, j, k))  M_u.x[i][j][k] = 0.0;
                M_u.x[i+1][j][k] = M_u.x[i][j][k] 
                    + step[i] * (E_u.x[i][j][k] - D_u.x[i][j][k]);    // in kg/(m²s) integration from cloud base upwards
                M_u_denom = M_u.x[i+1][j][k];
                M_u_num = M_u.x[i][j][k] - step[i] * D_u.x[i][j][k];
                if(M_u_denom == 0.0)  M_u_denom = 1e+6;
// conditions within the updraft
                u_u.x[i][j][k] = u.x[i][j][k] 
                    + M_u.x[i][j][k]/(r_humid.x[i][j][k] * area);  // in m/s
// dry entropy within the updraft
                s_u.x[i+1][j][k] = (M_u_num * s_u.x[i][j][k] 
                    + step[i] * (E_u.x[i][j][k] * s.x[i][j][k] 
                    + L_latent * r_humid.x[i][j][k] * c_u.x[i][j][k]))
                    /M_u_denom;
                if(s_u.x[i+1][j][k] > 1.5) s_u.x[i+1][j][k] = s.x[i+1][j][k];
                if(s_u.x[i+1][j][k] < 0.0) s_u.x[i+1][j][k] = s.x[i+1][j][k];
// cloud water within the updraft
                q_v_u.x[i+1][j][k] = (M_u_num * q_v_u.x[i][j][k]  // in kg/kg
                    + step[i] * (E_u.x[i][j][k] * c.x[i][j][k] 
                    - r_humid.x[i][j][k] * c_u.x[i][j][k]))/M_u_denom;
                if(q_v_u.x[i+1][j][k] >= q_Rain_add)
                    q_v_u.x[i+1][j][k] = q_Rain_add;
                if(q_v_u.x[i+1][j][k] < 0.0) q_v_u.x[i+1][j][k] = q_Rain_add;
// cloud water within the updraft
                q_c_u.x[i+1][j][k] = (M_u_num * q_c_u.x[i][j][k]   // in kg/kg
                    + step[i] * r_humid.x[i][j][k] 
                    * (c_u.x[i][j][k] - g_p.x[i][j][k]))/M_u_denom;
                if(q_c_u.x[i+1][j][k] < 0.0) q_c_u.x[i+1][j][k] = 0.0;
// velocity components within the updraft
                u_u.x[i+1][j][k] = u.x[i+1][j][k] 
                    + M_u.x[i+1][j][k]/(r_humid.x[i+1][j][k] * area); 
                v_u.x[i+1][j][k] = (M_u_num * v_u.x[i][j][k] 
                    + step[i] * E_u.x[i][j][k] * v.x[i][j][k])/M_u_denom;  // in m/s
                if(v_u.x[i+1][j][k] > .1) v_u.x[i+1][j][k] = .1;
                if(v_u.x[i+1][j][k] < -.1) v_u.x[i+1][j][k] = -.1;
                w_u.x[i+1][j][k] = (M_u_num * w_u.x[i][j][k] 
                    + step[i] * E_u.x[i][j][k] * w.x[i][j][k])/M_u_denom;  // in m/s
                if(w_u.x[i+1][j][k] > 10.0) w_u.x[i+1][j][k] = 10.0;
                if(w_u.x[i+1][j][k] < -10.0) w_u.x[i+1][j][k] = 10.0;
                if(t_u <= t_00){
                    M_u.x[i+1][j][k] = 0.0;
                    s_u.x[i+1][j][k] = 1.0;
                    q_v_u.x[i+1][j][k] = 0.0;
                    q_c_u.x[i+1][j][k] = 0.0;
                    u_u.x[i+1][j][k] = 0.0;
                    v_u.x[i+1][j][k] = 0.0;
                    w_u.x[i+1][j][k] = 0.0;
                    g_p.x[i][j][k] = 0.0;
                    c_u.x[i][j][k] = 0.0;
                    e_l.x[i][j][k] = 0.0;
                }
/*
            cout.precision(6);
            if((j == 90) &&(k == 180))  cout << endl 
                << "  MoistConvection   UPDRAFT" << endl
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
                << "  cloud_buoy = " << cloud_buoy 
                << "  q_buoy_base = " << q_buoy_base 
                << "  q_buoy_lfs = " << q_buoy_lfs << endl
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
                << "  t_u_add = " << t_add_u * t_u - t_0 << endl
                << "  HumidityRelative = " << HumidityRelative 
                << "  r_dry = " << r_dry.x[i][j][k] 
                << "  r_dry_add = " << r_dry_add 
                << "  r_humid = " << r_humid.x[i][j][k] 
                << "  r_humid_add = " << r_humid_add << endl
                << "  p_dyn = " << p_dyn.x[i][j][k] 
                << "  p_stat = " << p_stat.x[i][j][k] 
                << "  p_stat_add = " << p_stat_add
                << "  h = " << height 
                << "  h_base = " << height_base 
                << "  step = " << step[i] << endl
                << "  M_u = " << M_u.x[i][j][k] 
                << "  M_d_LFS = " << M_d_LFS[j][k] << endl;
*/
            } // end i convective updraft
        }  // end j
    }  // end k
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
// CONVECTIVE DOWNDRAFT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//            i_tropopause = get_tropopause_layer(j);
            double i_tropopause = im-1;
            i_base = (int)i_Base.y[j][k];
            i_lfs = (int)i_LFS.y[j][k];
            P_conv.x[i_lfs][j][k] = 0.0;
            for(int i = i_lfs; i > 0; i--){
                double cloud_buoy = r_dry.x[i][j][k] - r_humid.x[i][j][k];
                HumidityRelative = r_dry.x[i][j][k]/r_humid.x[i][j][k] * 100.;
                step[i] = get_layer_height(i) - get_layer_height(i-1);;  // local atmospheric shell thickness
                t_u = t.x[i][j][k] * t_0; // in K
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                t_u = t.x[i][j][k] * t_0; // in K
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                q_Rain = ep * E_Rain/(p_stat.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                if(t_u >= t_0) L_latent = lv;
                else L_latent = ls;
                s.x[i][j][k] = (cp_l * t_u + g * height)/s_0;
// evaporation of precipitation below the cloud base
                if(i <= i_base)
                    e_p.x[i][j][k] = C_p * alf_1 
                        * (q_Rain - q_v_d.x[i][j][k])
                        * sqrt(p_ps/alf_2 * P_conv.x[i][j][k]/C_p); // original COSMO formula
                if(e_p.x[i][j][k] <= 0.0) e_p.x[i][j][k] = 0.0;
// evaporation of rain precipitation within the downdraft
                e_d.x[i][j][k] = q_v_d.x[i][j][k] - q_Rain;
                if(e_d.x[i][j][k] <= 0.0) e_d.x[i][j][k] = 0.0;
// specification of entrainment and detrainment in the up/downdraft
                E_d.x[i][j][k] = eps_d * fabs(M_d.x[i][j][k]);
                D_d.x[i][j][k] = del_d * fabs(M_d.x[i][j][k]);
// mass flux within the downdraft
                if(is_land(h, i, j, k))  M_d.x[i][j][k] = 0.0;
                M_d.x[i-1][j][k] = M_d.x[i][j][k] 
                    + step[i] * (E_d.x[i][j][k] - D_d.x[i][j][k]);  // integration from LFS downwards
                M_d_denom = M_d.x[i-1][j][k];
                M_d_num = M_d.x[i][j][k] + step[i] * D_d.x[i][j][k];
                if(M_d_denom == 0.0)  M_d_denom = 1e+6;
// dry entropy within the downdraft
                if(t_u >= t_0) L_latent = lv;
                else L_latent = ls;
                s_d.x[i-1][j][k] = (M_d_num * s_d.x[i][j][k] 
                    + step[i] * (E_d.x[i][j][k] * s.x[i][j][k] 
                    - L_latent * r_humid.x[i][j][k] * e_d.x[i][j][k]))
                    /M_d_denom;
                if(s_d.x[i-1][j][k] > 1.5) s_d.x[i-1][j][k] = s.x[i-1][j][k];
                if(s_d.x[i-1][j][k] < 0.0) s_d.x[i-1][j][k] = s.x[i-1][j][k];
// cloud water within the downdraft
                q_v_d.x[i-1][j][k] = (M_d_num * q_v_d.x[i][j][k] 
                    + step[i] * (E_d.x[i][j][k] * c.x[i][j][k] 
                    + r_humid.x[i][j][k] * e_d.x[i][j][k]))/M_d_denom;
                if(q_v_d.x[i-1][j][k] >= q_Rain) 
                    q_v_d.x[i-1][j][k] = q_Rain;
                if(q_v_d.x[i-1][j][k] <= 0.0) q_v_d.x[i-1][j][k] = q_Rain;
// velocity components within the downdraft
                u_d.x[i-1][j][k] = u.x[i-1][j][k] 
                    + M_d.x[i-1][j][k]/(r_humid.x[i-1][j][k] * area); 
                v_d.x[i-1][j][k] = (M_d_num * v_d.x[i][j][k] 
                    + step[i] * E_d.x[i][j][k] * v.x[i][j][k])/M_d_denom;
                w_d.x[i-1][j][k] = (M_d_num * w_d.x[i][j][k] 
                    + step[i] * E_d.x[i][j][k] * w.x[i][j][k])/M_d_denom;
// rain water formed by cloud convection
                    if(u_d.x[i-1][j][k] <= 0.0){
                    P_conv.x[i-1][j][k] = P_conv.x[i][j][k] 
                        + step[i] * r_humid.x[i][j][k]
                        * (g_p.x[i][j][k] 
                        - e_d.x[i][j][k] 
                        - e_p.x[i][j][k]);
                    if(P_conv.x[i-1][j][k] <= 0.0) P_conv.x[i-1][j][k] = 0.0;
                    }
                if(t_u <= t_00){
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
            cout.precision(6);
            if((j == 90) &&(k == 180))  cout << endl 
                << "  MoistConvection   DOWNDRAFT" << endl
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
                << "  P_conv = " << P_conv.x[i][j][k] << endl
                << "  E_Rain = " << E_Rain 
                << "  q_Rain = " << q_Rain << endl
                << "  cloud_buoy = " << cloud_buoy 
                << "  q_buoy_base = " << q_buoy_base 
                << "  q_buoy_lfs = " << q_buoy_lfs << endl
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
                << "  r_humid = " << r_humid.x[i][j][k] << endl
                << "  p_dyn = " << p_dyn.x[i][j][k] 
                << "  p_stat = " << p_stat.x[i][j][k] 
                << "  p_stat = " << p_stat
                << "  h = " << height 
                << "  h_base = " << height_base 
                << "  step = " << step[i] << endl
                << "  M_d = " << M_d.x[i][j][k] 
                << "  M_d_LFS = " << M_d_LFS[j][k] << endl;
*/
            } // end i convective downdraft
            e_d.x[0][j][k] = c43 * e_d.x[1][j][k] - c13 * e_d.x[2][j][k];
            e_p.x[0][j][k] = c43 * e_p.x[1][j][k] - c13 * e_p.x[2][j][k];
            E_d.x[0][j][k] = c43 * E_d.x[1][j][k] - c13 * E_d.x[2][j][k];
            D_d.x[0][j][k] = c43 * D_d.x[1][j][k] - c13 * D_d.x[2][j][k];
        }  // end j
    }  // end k
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
// RHS for thermodynamic forcing due to moist convection
            for(int i = 0; i < im-1; i++){
//                int i_mount = i_topography[j][k];
                step[i] = get_layer_height(i+1) - get_layer_height(i);  // local atmospheric shell thickness
                if((t.x[i][j][k] * t_0) >= t_0) L_latent = lv;
                else L_latent = ls;
                MC_t.x[i][j][k] = 
                    - ((M_u.x[i+1][j][k] * (s_u.x[i+1][j][k] - s.x[i+1][j][k]) 
                    + M_d.x[i+1][j][k] * (s_d.x[i+1][j][k] - s.x[i+1][j][k])) 
                    - (M_u.x[i][j][k] * (s_u.x[i][j][k] - s.x[i][j][k]) 
                    + M_d.x[i][j][k] * (s_d.x[i][j][k] - s.x[i][j][k]))) 
                    /(r_humid.x[i][j][k] * step[i])/cp_l * s_0
                    + L_latent/cp_l * (c_u.x[i][j][k] - e_d.x[i][j][k] 
                    - e_l.x[i][j][k] - e_p.x[i][j][k]);  // in (K)(1/s)
                MC_q.x[i][j][k] = 
                    - ((M_u.x[i+1][j][k] * (q_v_u.x[i+1][j][k] - c.x[i+1][j][k]) 
                    + M_d.x[i+1][j][k] * (q_v_d.x[i+1][j][k] - c.x[i+1][j][k])) 
                    - (M_u.x[i][j][k] * (q_v_u.x[i][j][k] - c.x[i][j][k]) 
                    + M_d.x[i][j][k] * (q_v_d.x[i][j][k] - c.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k])
                    - (c_u.x[i][j][k] - e_d.x[i][j][k] - e_l.x[i][j][k] 
                    - e_p.x[i][j][k]);  // in (kg/kg)(1/s)
                MC_v.x[i][j][k] = 
                    - ((M_u.x[i+1][j][k] * (v_u.x[i+1][j][k] - v.x[i+1][j][k]) 
                    + M_d.x[i+1][j][k] * (v_d.x[i+1][j][k] - v.x[i+1][j][k])) 
                    - (M_u.x[i][j][k] * (v_u.x[i][j][k] - v.x[i][j][k]) 
                    + M_d.x[i][j][k] * (v_d.x[i][j][k] - v.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k]) * u_0;  // in (m/s)(1/s)
                MC_w.x[i][j][k] = 
                    - ((M_u.x[i+1][j][k] * (w_u.x[i+1][j][k] - w.x[i+1][j][k]) 
                    + M_d.x[i+1][j][k] * (w_d.x[i+1][j][k] - w.x[i+1][j][k])) 
                    - (M_u.x[i][j][k] * (w_u.x[i][j][k] - w.x[i][j][k]) 
                    + M_d.x[i][j][k] * (w_d.x[i][j][k] - w.x[i][j][k]))) 
                    /(step[i] * r_humid.x[i][j][k]) * u_0;  // in (m/s)(1/s)
            }  // MC-loop
        } // j-loop
    } // k-loop
/*
*
*/
    fft_gaussian_filter_3d(P_conv,1);
    fft_gaussian_filter_3d(MC_t,1);
    fft_gaussian_filter_3d(MC_q,1);
    fft_gaussian_filter_3d(MC_v,1);
    fft_gaussian_filter_3d(MC_w,1);
/*
*
*/
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
            P_conv.x[i][0][k] = c43 * P_conv.x[i][1][k] -
                 c13 * P_conv.x[i][2][k];
            P_conv.x[i][jm-1][k] = c43 * P_conv.x[i][jm-2][k] -
                 c13 * P_conv.x[i][jm-3][k];
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
            P_conv.x[i][j][0] = c43 * P_conv.x[i][j][1] -
                c13 * P_conv.x[i][j][2];
            P_conv.x[i][j][km-1] = c43 * P_conv.x[i][j][km-2] -
                c13 * P_conv.x[i][j][km-3];
            P_conv.x[i][j][0] = P_conv.x[i][j][km-1] =
               (P_conv.x[i][j][0] + P_conv.x[i][j][km-1])/2.;
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i <= i_mount)){
//                    t.x[i][j][k] = t.x[i_mount][j][k];
                    P_conv.x[i][j][k] = P_conv.x[i_mount][j][k];
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
//                    c.x[i][j][k] = c.x[i_mount][j][k];
//                    cloud.x[i][j][k] = cloud.x[i_mount][j][k];
//                    ice.x[i][j][k] = ice.x[i_mount][j][k];
                }
            }
        }
    }
}
/*
*
*/
void cAtmosphereModel::WaterVapourEvaporation(){ 
// preparations for water vapour increase due to the differences between evaporation and precipitation 
// procedure given in Rui Xin Huang, Ocean Circulation, p. 165
// amount of additional water vapour by the difference of evaporation and precipitation is negligible but functional
//    double c_surf_ocean = .8; // 80% humidity
//    double c_surf_ocean = .9; // 90% humidity
    double c_surf_ocean = 1.; // 100% humidity
//    double rm = log(get_layer_height(0)/L_atm + 1.);
    double rm = log(rad.z[0] + 1.);
    double exp_rm = 1./(rm + 1.);
    double evap_precip = 0.0;
    double vapour_surface_n = 0.0;
    double vapour_surface = 0.0;
    double coeff_evap = 1.1574e-8;  // 1.1574-8 is the conversion from (Evap-Prec) in mm/d to m/s
    float e, E_Rain, E_Ice, sat_deficit;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            int i = get_surface_layer(j,k);
            c.x[0][j][k] = c_surf_ocean * hp * ep * exp(17.0809 
                * (t.x[0][j][k] * t_0 - t_0)/(234.175 
                + (t.x[0][j][k] * t_0 - t_0)))/((r_air * R_Air 
                * t.x[0][j][k] * t_0) * .01);
            c_fix.y[j][k] = c.x[0][j][k];
            float t_u = t.x[i][j][k] * t_0;
            float t_Celsius = t_u - t_0;
            if(t_Celsius >= 0.0){
                if(is_land(h, i, j, k)){
                    t_u = (t.x[0][j][k] + t_land) * t_0;
//                    e = c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // water vapour pressure in hPa
                    e = c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // water vapour pressure in hPa
                    E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                    sat_deficit = E_Rain - e;  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = C_Dalton(i+1, j, k, 
                        coeff_Dalton, u_0, v, w) * sat_deficit * 24.;  // mm/h in mm/d
                } // simplified formula for Evaporation by Dalton law dependent on surface water velocity in kg/(m²*s) = mm/s
                if(is_water(h, 0, j, k)){
//                    e = c.x[0][j][k] * p_stat.x[0][j][k]/ep;  // water vapour pressure in hPa
                    e = c.x[0][j][k] * p_stat.x[0][j][k]/ep;  // water vapour pressure in hPa
                    t_u = t.x[0][j][k] * t_0;
                    E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                    sat_deficit = E_Rain - e;  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = C_Dalton(0, j, k, 
                        coeff_Dalton, u_0, v, w) * sat_deficit * 24.;  // mm/h in mm/d
                }
            }
            if(t_Celsius <= 0.0){
                if(is_land(h, i, j, k)){
                    t_u = (t.x[0][j][k] + t_land) * t_0;
//                    e = c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // water vapour pressure in hPa
                    e = c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // water vapour pressure in hPa
                    E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                    sat_deficit = (E_Ice - e);  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = C_Dalton(i+1, j, k, coeff_Dalton, u_0, v,
                        w) * sat_deficit * 24.;  // mm/h in mm/d
                }
                if(is_water(h, 0, j, k)){
                    e = c.x[0][j][k] * p_stat.x[0][j][k]/ep;  // water vapour pressure in hPa
                    t_u = (t.x[0][j][k] + t_land) * t_0;
                    E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                    sat_deficit = (E_Ice - e);  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = C_Dalton(0, j, k, coeff_Dalton, u_0, v,
                        w) * sat_deficit * 24.;  // mm/h in mm/d
                }
            }
            if(Evaporation_Dalton.y[j][k] <= 0.0)  
                Evaporation_Dalton.y[j][k] = 0.0;
            for(int iter_prec = 1; iter_prec <= 10; iter_prec++){ // iter_prec may be varied
                evap_precip = coeff_evap * (Evaporation_Dalton.y[j][k] - Precipitation.y[j][k]); // in mm/s
                vapour_surface = - ((- 3. * c.x[0][j][k] + 4. * c.x[1][j][k] 
                    - c.x[2][j][k])/(2. * dr) * exp_rm // 1. order derivative, 2. order accurate
                    * (1. - 2. * c.x[0][j][k]) * evap_precip) * r_humid.x[0][j][k]; 
//                vapour_surface = - (c.x[1][j][k] - c.x[0][j][k]) // 1. order derivative, 1. order accurate
//                   /dr * exp_rm * (1. + 2. * c.x[0][j][k]) * evap_precip;
//                if(iter_prec == 1)  vapour_surface_n = .9 * vapour_surface;
                if(iter_prec == 1)  vapour_surface = .9 * vapour_surface_n;
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
                    << "  vap_diff = " << fabs(vapour_surface/vapour_surface_n - 1.) 
                    << "  vap_surf = " << vapour_surface 
                    << "  vap_surf_n = " << vapour_surface_n << endl
                    << "  c_fix = " << c_fix.y[j][k] * 1e3 
                    << "  c = " << c.x[0][j][k] * 1e3 << endl
                    << "  v = " << v.x[0][j][k] * u_0
                    << "  w = " << w.x[0][j][k] * u_0
                    << "  vel_magnitude = " << sqrt(v.x[0][j][k] * v.x[0][j][k] 
                        + w.x[0][j][k] * w.x[0][j][k]) * u_0 
                    << "  c_grad_1 = " << (c.x[1][j][k] - c.x[0][j][k])/dr * exp_rm
                    << "  c_grad_2 = " << (- 3. * c.x[0][j][k] + 4. * c.x[1][j][k] 
                        - c.x[2][j][k]) /(2. * dr) * exp_rm << endl
                    << "  sat_deficit = " << sat_deficit 
                    << "  e = " << e 
                    << "  E_Rain = " << E_Rain << endl
                    << "  C_Dalton = " << C_Dalton(i+1, j, k, coeff_Dalton, u_0, v, w)
                        * sat_deficit * 2.778e-4
                    << "  Evap-Prec = " << evap_precip 
                    << "  Evap = " << Evaporation_Dalton.y[j][k] 
                    << "  Prec = " << Precipitation.y[j][k] << endl;
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
                if(fabs(vapour_surface/vapour_surface_n - 1.) 
                    <= 1.e-5)  break;  
                vapour_surface_n = vapour_surface;
            } // iter_prec 
        } // end j
    } // end k
}
/*
*
*/
void cAtmosphereModel::USStand_DewPoint_HumidRel(){
    temp_tropopause = std::vector<double>(jm, t_tropopause_pole);
    int j_max = jm-1;
    int j_half = j_max/2;
    double e, E, t_u;
    double temp_eff = t_tropopause_pole - t_tropopause_equator;
    for(int j=j_half; j>=0; j--){
//        double x = 1e-3 * (double)(j_half-j)/(double)j_half;
//        double x = (double)(j_half-j)/(double)j_half;
//        temp_tropopause[j] = Agnesi(t_tropopause_equator, x); 
        temp_tropopause[j] = temp_eff * parabola((double)j
            /(double)j_half) + t_tropopause_pole;
    }
    for(int j=j_max; j>j_half; j--){
        temp_tropopause[j] = temp_tropopause[j_max-j];
    }
    AtomUtils::smooth_tropopause(jm, temp_tropopause);
    for(int j = 0; j < jm; j++){
//        int i_trop = get_tropopause_layer(j);
        int i_trop = im-1;
        for(int k = 0; k < km; k++){
//            int i_mount = i_topography[j][k];
            int i_mount = 0;
            for(int i = i_mount; i < im; i++){
                t_u = t.x[i][j][k] * t_0; // in K
                TempStand.x[0][j][k] = t.x[0][j][k] * t_0 - t_0;
                e = c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // water vapour pressure in hPa
                if(t_u >= t_0)
                    E = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                else
//                    E = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the water phase at t < 0°C in hPa
                    E = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                if(i <= i_trop){
                        // linear temperature decay up to tropopause
                        // compares to the US Standard Atmosphere 
                        // with variable tropopause locations
                        // US Standard Atmosphere
                        TempStand.x[i][j][k] = TempStand.x[0][j][k] 
                            - 6.5 * (double)(get_layer_height(i))/1000.0;
                        TempDewPoint.x[i][j][k] = (423.86 - 234.175 * log(e))
                            /(log(e) - 18.89);  // in C°, Häckel, Meteorologie, p. 82
                        HumidityRel.x[i][j][k] = e/E * 100.0;
                        if(HumidityRel.x[i][j][k] > 100.0)
                            HumidityRel.x[i][j][k] = 100.0;
                }else{ // above tropopause
                    TempStand.x[i][j][k] = TempStand.x[i_trop][j][k];
                    TempDewPoint.x[i][j][k] = TempDewPoint.x[i_trop][j][k];
                    HumidityRel.x[i][j][k] = HumidityRel.x[i_trop][j][k];
                }
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i >= 0; i--){
                if((is_land(h, i, j, k))&&(i <= i_mount)){
                    TempStand.x[i][j][k] = TempStand.x[i_mount][j][k];
                    TempDewPoint.x[i][j][k] = TempDewPoint.x[i_mount][j][k];
                    HumidityRel.x[i][j][k] = HumidityRel.x[i_mount][j][k];
                }
            }
        }
    }
}
/*
*
*/
void cAtmosphereModel::PressureLimitationAtm(){
    // class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(is_land(h, i, j, k)){
                    p_dyn.x[i][j][k] = 0.0;
                }
                p_dynn.x[i][j][k] = p_dyn.x[i][j][k];
            }
        }
    }
}
