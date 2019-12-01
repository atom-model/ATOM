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

double get_pole_temperature(int Ma, const std::map<float, float> &pole_temp_map);

void cAtmosphereModel::BC_Radiation_multi_layer(){
    // class element for the computation of the radiation and the temperature distribution
    // computation of the local temperature based on short and long wave radiation
    // multi layer radiation model
    if(debug){
        Array tmp = (t-1)*t_0;
        logger()<<"20180912: Enter RML ... "<<std::endl;
        tmp.inspect("20180912: ");
    }
    temp_tropopause = std::vector<double>(jm, t_tropopause_pole);
    int j_max = jm-1;
    int j_half = j_max/2;
    for(int j=j_half; j>=0; j--){
        double x = sqrt(pow(t_tropopause,3) / t_tropopause_pole 
            - pow(t_tropopause,2)) * (double)(j_half-j) 
            / (double)j_half;
        temp_tropopause[j] = Agnesi(x, t_tropopause); 
    }
    for(int j=j_max; j>j_half; j--){
        temp_tropopause[j] = temp_tropopause[j_max-j];
    }
    AtomUtils::smooth_tropopause(jm, temp_tropopause);
    std::map<float, float> pole_temp_map;  // Stein/Rüdiger/Parish linear pole temperature (Ma) distribution
    load_map_from_file(pole_temperature_file, pole_temp_map); 
    double rad_eff = rad_pole - rad_equator;
    double albedo_co2_eff = albedo_pole - albedo_equator;
    double j_max_half = (jm-1)/2;
    // effective temperature, albedo and emissivity/absorptivity for the two layer model
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im-1; i++){
                if(is_ocean_surface(h, i, j, k) || is_land_surface(h, i, j, k)){
                    albedo.y[j][k] = albedo_co2_eff * parabola((double)j
                        /(double)j_max_half) + albedo_pole;
                }
            }
        }
    }
    // absorption/emissivity computation
    double epsilon_eff_max = .594; // constant  given by Häckel (F. Baur and H. Philips, 1934)
    // constant value stands for other non-condensable gases than water vapour in the equation for epsilon
    double epsilon_eff_2D = epsilon_pole - epsilon_equator;
    for(int j = 0; j < jm; j++){
        int i_trop = m_model->get_tropopause_layer(j);
        // on zero level, lateral parabolic distribution
        epsilon_eff_max = epsilon_eff_2D * parabola((double)j
            /(double)j_max_half) + epsilon_pole;
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            // shortwave radiation from the sun absobed at the surface 47%
            // in W/m², assumption of parabolic surface radiation at zero level
            radiation_surface.y[j][k] = rad_eff * parabola((double)j
                /(double)j_max_half) + rad_pole;
            for(int i = 0; i <= i_trop; i++){
                if(c.x[i][j][k] < 0.)      c.x[i][j][k] = 0.;
                if(cloud.x[i][j][k] < 0.)  cloud.x[i][j][k] = 0.;
                if(ice.x[i][j][k] < 0.)    ice.x[i][j][k] = 0.;
                // COSMO water vapour pressure based on local water vapour in hPa
                double e = c.x[i][j][k] * p_stat.x[i][j][k] / ep;
                // radial parabolic distribution, start on zero level
                double epsilon_eff = epsilon_eff_max 
                    * (1. - get_layer_height(i) / get_layer_height(i_trop));
                // dependency given by Häckel, Meteorologie, p. 205 (law by F. Baur and H. Philips, 1934)
                // epsilon_eff describe the effect in the emissivity computation of other gases like CO2
                // in the original formula this value is 0.594, for reasons of adjustment to the modern atmosphere,
                // this constant became a variable in zonal direction
                // this variable reacts very sensitive and changes the temperature field extremely
                // the second term describes the influence of water vapour only 
                if(i >= i_mount){ //start from the mountain top
                    double step = get_layer_height(i+1) - get_layer_height(i); // local atmospheric shell thickness
                    double P_c = (1.e-6 * p_stat.x[i][j][k] 
                                 * co2.x[i][j][k] * co2_0 / p_0);
                    double u_c = P_c * step * 100. * (L_atm / (double)(im-1)); // model by Atwater and Ball
//                    u_c = 300. / 0.9869 * P_c * step * 100. / (t.x[i][j][k] * t_0); 
                            // modell by Yamamoto and Sasamori
                 // influence on the emissivity by carbon dioxcide by the law by Bliss
//                    x_c = 10. * P_c / (R_co2 * t.x[i][j][k] * t_0) * step * p_0;
                  // coefficient in the Bliss law, 10. => conversion from kg/m² to g/cm²
//                    eps_co2 = .185 * (1. - exp (- 50. * x_c)); // modell by Bliss
                    double eps_co2 = .185 * (1. - exp(- 0.3919 * pow(u_c, 0.4))); // model by Atwater and Ball
                    eps_co2 = .5 * eps_co2; // model by Atwater and Ball, 
                                            // factor .5 fits the Atwater and Ball results to Yamamoto and Sasamori,
                                            // which fit best to HITRAN band data
//                    eps_co2 = 1.5 * eps_co2; // Test, 1.5 is maximum, otherwise nans develop

                    epsilon_3D.x[i][j][k] = eps_co2 + epsilon_eff 
                        + .0416 * sqrt(e);
                    // atmospheric emissivity by Vogel/Bliss, emissivity seems too high
//                    epsilon_3D.x[i][j][k] = 0.98 * pow((e/(t.x[i][j][k] * t_0)), 0.0687);

                    radiation_3D.x[i][j][k] = (1.-epsilon_3D.x[i][j][k]) 
                        * sigma * pow(t.x[i][j][k] * t_0, 4.);
                }
                if(epsilon_3D.x[i][j][k] > 1.)  epsilon_3D.x[i][j][k] = 1.;
            }
            epsilon.y[j][k] = epsilon_3D.x[0][j][k];
            // inside mountains
            for(int i = i_mount - 1; i >= 0; i--){
                epsilon_3D.x[i][j][k] = epsilon_3D.x[i_mount][j][k];
                radiation_3D.x[i][j][k] = radiation_3D.x[i_mount][j][k];
            }
            //above tropopause
            for(int i = i_trop; i < im; i++){
                epsilon_3D.x[i][j][k] = epsilon_3D.x[i_trop][j][k];
                t.x[i][j][k] = temp_tropopause[j];
                radiation_3D.x[i][j][k] = (1. - epsilon_3D.x[i][j][k]) 
                    * sigma * pow(t.x[i][j][k] * t_0, 4.);
            }
        }
    }
    // iteration procedure for the computation of the temperature based on the multi-layer radiation model
    // temperature needs an initial guess which must be corrected by the long wave radiation remaining in the atmosphere
    for(int iter_rad = 1;  iter_rad <= 4; iter_rad++){ // iter_rad may be varied
        // coefficient formed for the tridiogonal set of equations for the absorption/emission coefficient of the multi-layer radiation model
        for(int j = 0; j < jm; j++){
            int i_trop = get_tropopause_layer(j);
            for(int k = 0; k < km; k++){
                int i_mount = i_topography[j][k];
                std::vector<double> alfa(im, 0);
                std::vector<double> beta(im, 0);
                std::vector<double> AA(im, 0);
                std::vector<std::vector<double> > CC(im, std::vector<double>(im, 0));
                double CCC = 0, DDD = 0;
                // radiation leaving the atmosphere above the tropopause, 
                // emitted of atmophere 49%, Clouds 9% and 12% atmospheric window
                radiation_3D.x[i_trop][j][k] = .70 * radiation_surface.y[j][k]; 
                // longwave back radiation absorbed by the surface 116% of the shortwave incoming sun radiation
                double radiation_back = radiation_surface.y[j][k];
                // longwave back radiation loss of 12% through the atmospheric window
                double atmospheric_window = .12 * radiation_surface.y[j][k]; 
                // latent heat leaving the surface for cloud formation 24% of the shortwave incoming sun radiation
                double Latent_Heat = .24 * radiation_surface.y[j][k];
                // sensible heat leaving the surface by convection 5% of the shortwave incoming sun radiation
                double Sensible_Heat = .05 * radiation_surface.y[j][k];
                double rad_surf_diff = radiation_back + Latent_Heat
                    + Sensible_Heat + atmospheric_window; // radiation leaving the surface
                AA[i_mount] = rad_surf_diff / radiation_3D.x[i_trop][j][k];// non-dimensional surface radiation
                CC[i_mount][i_mount] = 0.; // no absorption of radiation on the surface by water vapour
                radiation_3D.x[i_mount][j][k] = 
                    (1. - epsilon_3D.x[i_mount][j][k]) * sigma
                    * pow(t.x[i_mount][j][k] * t_0, 4.) 
                    / radiation_3D.x[i_trop][j][k]; // radiation leaving the surface
                for(int i = i_mount+1; i <= i_trop; i++){
                    AA[i] = AA[i-1] * (1. - epsilon_3D.x[i][j][k]); // transmitted radiation from each layer
                    double rad_tmp = sigma * pow(t.x[i][j][k] * t_0, 4.) 
                        / radiation_3D.x[i_trop][j][k];
                    CC[i][i]= epsilon_3D.x[i][j][k] * rad_tmp; // absorbed radiation in each layer
                    radiation_3D.x[i][j][k] = 
                        (1. - epsilon_3D.x[i][j][k]) * rad_tmp; // radiation leaving each layer
                    for(int l = i_mount+1; l <= i_trop; l++){
                        // additional transmitted radiation from layer to layer in radial direction
                        CC[i][l] = CC[i][l-1] * (1. - epsilon_3D.x[l][j][k]);
                    }
                }
                // Thomas algorithm to solve the tridiogonal equation system for the solution of the radiation with a recurrence formula
                // additionally embedded in an iterational process
                double aa, bb, cc, dd;
                for(int i = i_mount; i < i_trop; i++){ // values at the surface
                    if(i == i_mount){
                        bb = - radiation_3D.x[i][j][k];
                        cc = radiation_3D.x[i+1][j][k];
                        dd = - AA[i];
                        alfa[i] = cc / bb;
                        beta[i] = dd / bb;
                    }else{
                        for(int l = i_mount+1; l <= i-1; l++){
                            CCC = CCC + CC[l][i];
                        }
                        for(int l = i_mount+1; l <= i-2; l++){
                            DDD = DDD + CC[l][i-1];
                        }
                        aa = radiation_3D.x[i-1][j][k];
                        bb = - 2. * radiation_3D.x[i][j][k];
                        cc = radiation_3D.x[i+1][j][k];
                        dd = - AA[i-1] + AA[i] + CCC - DDD;
                        alfa[i] = cc / (bb - aa * alfa[i-1]);
                        beta[i] = (dd - aa * beta[i-1]) 
                            / (bb - aa * alfa[i-1]);
                    }
                }
                // radiation leaving the atmosphere above the tropopause, later needed for non-dimensionalisation
                t.x[i_trop][j][k] = temp_tropopause[j];
                radiation_3D.x[i_trop][j][k] = 
                    (1. - epsilon_3D.x[i_trop][j][k]) * sigma
                    * pow(t.x[i_trop][j][k] * t_0, 4.);
                // recurrence formula for the radiation and temperature
                for(int i = i_trop - 1; i >= i_mount; i--){
                    // above assumed tropopause constant temperature t_tropopause
                    // Thomas algorithm, recurrence formula
                    radiation_3D.x[i][j][k] = - alfa[i] 
                        * radiation_3D.x[i+1][j][k] + beta[i];
                    t.x[i][j][k] = pow(radiation_3D.x[i][j][k] / sigma, 
                        (1. / 4.)) / t_0;    // averaging of temperature values to smooth the iterations
                } // end i
                for(int i = i_trop; i < im; i++){ // above tropopause
                    t.x[i][j][k] = temp_tropopause[j];
                    radiation_3D.x[i][j][k] = 
                        (1. - epsilon_3D.x[i_trop][j][k]) * sigma
                        * pow(t.x[i][j][k] * t_0, 4.);
                } // end i
            } // end k
        } // end j
    }
    logger() << "exit BC_Radiation_multi_layer: temperature max: " 
        << (t.max() - 1)*t_0 << std::endl << std::endl;
    if(debug){
        Array tmp = (t-1)*t_0;
        logger()<<"20180912: Exit RML ... "<<std::endl;
        tmp.inspect("20180912: ");
    }
}


void cAtmosphereModel::init_temperature(){
    if(debug){
        Array tmp = (t-1)*t_0;
        logger()<<"20180912: Enter BCT ... "<<std::endl;
        tmp.inspect("20180912: ");
    }
    // logger() << std::endl << "enter BC_Temperature: temperature max: " << (t.max()-1)*t_0 << std::endl;
    // logger() << "enter BC_Temperature: temperature min: " << (t.min()-1)*t_0 << std::endl << std::endl;

    // Lenton_etal_COPSE_time_temp, constant paleo mean temperature, added to the surface initial temperature
    // difference between mean temperature (Ma) and mean temperature (previous Ma) == t_paleo_add
    double t_paleo_add = 0; 
    if(!is_first_time_slice()){
        if((NASATemperature != 0)&&(*get_current_time() > 0))  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - get_mean_temperature_from_curve(*get_previous_time());
        if(NASATemperature == 0)  
            t_paleo_add = get_mean_temperature_from_curve(*get_current_time())
                - t_average;
        t_paleo_add /= t_0; // non-dimensional 
    }
    cout.precision(3);
    int Ma = *get_current_time();
    const char* time_slice_comment = "      time slice of Paleo-AGCM:";
    const char* time_slice_number = " Ma = ";
    const char* time_slice_unit = " million years";
    cout << endl << setiosflags(ios::left) << setw(55) << setfill('.') 
        << time_slice_comment << resetiosflags(ios::left) << setw (6) 
        << fixed << setfill(' ') << time_slice_number << setw (3) << Ma 
        << setw(12) << time_slice_unit << endl << endl;
    const char* temperature_comment = "      temperature increase at paleo times: ";
    const char* temperature_gain = " t increase";
    const char* temperature_modern = "      mean temperature at modern times: ";
    const char* temperature_paleo = "      mean temperature at paleo times: ";
    const char* temperature_average = " t modern";
    const char* temperature_average_pal = " t paleo";
    const char* temperature_unit =  "°C ";
    cout << endl << setiosflags(ios::left) << setw(55) << setfill('.') 
        << temperature_comment << resetiosflags(ios::left) << setw(13) 
        << temperature_gain << " = " << setw(7) << setfill(' ') 
        << t_paleo << setw(5) << temperature_unit << endl << setw(55) 
        << setfill('.') << setiosflags(ios::left) << temperature_modern 
        << resetiosflags(ios::left) << setw(13) << temperature_average 
        << " = " << setw(7) << setfill(' ') << t_average << setw(5) 
        << temperature_unit << endl << setw(55) << setfill('.') 
        << setiosflags(ios::left) << temperature_paleo << resetiosflags(ios::left) 
        << setw(13) << temperature_average_pal << " = "  << setw(7) 
        << setfill(' ') << t_average + t_paleo << setw(5) 
        << temperature_unit << endl;
    // temperatur distribution at a prescribed sun position
    // sun_position_lat = 60,    position of sun j = 120 means 30°S, j = 60 means 30°N
    // sun_position_lon = 180, position of sun k = 180 means 0° or 180° E (Greenwich, zero meridian)
    // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature (linear equation + parabola)
    if((*get_current_time() > 0) && (sun == 1)){
        double j_par = sun_position_lat; // position of maximum temperature, sun position
        j_par = j_par + declination; // angle of sun axis, declination = 23,4°
        double j_pol = jm-1;
        double j_par_f = (double)j_par;
        double j_pol_f = (double)j_pol;
        double aa = (t_equator - t_pole) / (((j_par_f * j_par_f) 
            - (j_pol_f * j_pol_f)) - 2. * j_par_f * (j_par_f - j_pol_f));
        double bb = - 2. * aa * j_par_f;
        double cc = t_equator + aa * j_par_f * j_par_f;
        double j_d = sqrt ((cc - t_pole) / aa);
        double dd = 2. * aa * j_d + bb;
        double e = t_pole;
        // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature (linear equation + parabola)
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                double d_j = (double)j;
                if(d_j <= j_d){
                    t.x[0][j][k] = dd * d_j + e + t_paleo_add;
                }
                if(d_j > j_d){
                    t.x[0][j][k] = aa * d_j * d_j + bb * d_j 
                        + cc + t_paleo_add;
                }
            }
        }
        // longitudinally variable temperature distribution from west to east in parabolic form
        // communicates the impression of local sun radiation on the southern hemisphere
        double k_par = sun_position_lon;  // position of the sun at constant longitude
        double k_pol = km - 1;
        double t_360 = (t_0 + 5.) / t_0;
        for(int j = 0; j < jm; j++){
            double jm_temp_asym = t.x[0][j][20];//transfer of zonal constant temperature into aa 1D-temperature field
            for(int k = 0; k < km; k++){
                double k_par_f = (double)k_par;
                double k_pol_f = (double)k_pol;
                double d_k = (double) k;
                aa = (jm_temp_asym - t_360) / (((k_par_f * k_par_f) 
                    - (k_pol_f * k_pol_f)) - 2. * k_par_f * 
                    (k_par_f - k_pol_f));
                bb = - 2. * aa * k_par_f;
                cc = jm_temp_asym + aa * k_par_f * k_par_f;
                t.x[0][j][k] = aa * d_k * d_k + bb * d_k + cc;
            }
        }
    }// temperatur distribution at aa prescribed sun position
    // pole temperature adjustment, combination of linear time dependent functions 
    // Stein/Rüdiger/Parish locally constant pole temperature
    // difference between pole temperature (Ma) and pole temperature (previous Ma)
    std::map<float, float> pole_temp_map;  // Stein/Rüdiger/Parish linear pole temperature (Ma) distribution
    load_map_from_file(pole_temperature_file, pole_temp_map); 
    double d_j_half = (double)(jm-1)/2.;
    float t_pole_diff_ocean = 0.;
    float t_pole_diff_land = 0.;
    //the t_pole_diff_ocean should be the difference between this time slice and the previous one
    if(!is_first_time_slice()){
        t_pole_diff_ocean = get_pole_temperature(*get_current_time(), pole_temp_map)
            - get_pole_temperature(*get_previous_time(), pole_temp_map);
        t_pole_diff_land = get_pole_temperature(*get_current_time(), pole_temp_map)
            - get_pole_temperature(*get_previous_time(), pole_temp_map);
    }
    // in °C, constant local pole temperature as function of Ma for hothouse climates 
    float pole_temperature = 1 + get_pole_temperature(*get_current_time(), 
        pole_temp_map) / t_0;
    float t_eff = pole_temperature - t_equator;  // coefficient for the zonal parabolic temperature distribution
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double d_j = (double)j;
            if(NASATemperature == 0){  // parabolic ocean surface temperature assumed
                t.x[0][j][k] = t_eff * parabola((double)d_j/(double)d_j_half) 
                    + pole_temperature + t_paleo_add;
//                srand(time(NULL));
//                t.x[0][j][k] += (rand() % 10 - 5) / 50. / t_0;
                if(is_land(h, 0, j, k)){
                    t.x[0][j][k] += m_model->t_land;
                }
            }else{  // if(NASATemperature == 1) ocean surface temperature based on NASA temperature distribution
                // transported for later time slices Ma by use_earthbyte_reconstruction
                if(is_land (h, 0, j, k)){  // on land a parabolic distribution assumed, no NASA based data transportable
                    if(*get_current_time() > 0){
//                        t.x[0][j][k] = (t_eff * parabola(d_j/d_j_half) 
//                            + pole_temperature) + t_paleo_add + m_model->t_land;
                        t.x[0][j][k] += t_paleo_add + m_model->t_land
                            + t_pole_diff_land * fabs(parabola((double)d_j
                            /(double)d_j_half) + 1.) / t_0;
                        // land surface temperature increased by mean t_paleo_add
                        // and by a zonally equator wards decreasing temperature difference is added
                        // Stein/Rüdiger/Parish pole temperature decreasing equator wards
                    }
                    if(*get_current_time() == 0)
                        t.x[0][j][k] = temperature_NASA.y[j][k];  // initial temperature by NASA for Ma=0
                }else{ // if the location is ocean
                    if(*get_current_time() > 0){        
                        // ocean surface temperature increased by mean t_paleo_add
                        // and by a zonally equator wards decreasing temperature difference is added
                        // Stein/Rüdiger/Parish pole temperature decreasing equator wards
                        t.x[0][j][k] += t_paleo_add 
                            + t_pole_diff_ocean * fabs(parabola((double)d_j 
                            / (double)d_j_half) + 1.) / t_0;
                    }
                    if(*get_current_time() == 0)
                        t.x[0][j][k] = temperature_NASA.y[j][k];  // initial temperature by NASA for Ma=0
                }
            }// else(NASATemperature == 1)
        }// for j
    }// for k
    // zonal temperature along tropopause
//    double t_eff_tropo = t_tropopause_pole - t_tropopause;
    //use "linear temperature decay" to generate temperature data for layers between mountain top and tropopause
    //use "mountain top temperature" for the layers below mountain top
    //use "tropopause tempeature" for the layers above tropopause
    // temperature approaching the tropopause, above constant temperature following Standard Atmosphere
    temp_tropopause = std::vector<double>(jm, t_tropopause_pole);
    int j_max = jm-1;
    int j_half = j_max/2;
    for(int j=j_half; j>=0; j--){
        double x = sqrt(pow(t_tropopause,3) / t_tropopause_pole 
            - pow(t_tropopause,2)) * (double)(j_half-j) 
            / (double)j_half;
        temp_tropopause[j] = Agnesi(x, t_tropopause); 
/*
    cout << "   j = " << j << "   x = " << x 
        << "   Agnesi = " << Agnesi(x, t_tropopause) 
        << "   t_tropopause = " << t_tropopause 
        << "   temp_tropopause = " << temp_tropopause[j] << endl;
*/
    }
    for(int j=j_max; j>j_half; j--){
        temp_tropopause[j] = temp_tropopause[j_max-j];
    }
    AtomUtils::smooth_tropopause(jm, temp_tropopause);
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            int i_trop = get_tropopause_layer(j);
//            for(int i = i_mount; i < im; i++){
            for(int i = 0; i < im; i++){
                if(i < i_trop){
                    if(i>i_mount){
                        // linear temperature decay up to tropopause, privat  approximation
                        t.x[i][j][k] = (temp_tropopause[j] - t.x[0][j][k]) * 
                            (get_layer_height(i) / get_layer_height(i_trop)) 
                            + t.x[0][j][k]; 
                        // US Standard Atmosphere
                    }else{
                        t.x[i][j][k] = (temp_tropopause[j] - t.x[0][j][k]) * 
                            (get_layer_height(i) / get_layer_height(i_trop)) 
                            + t.x[0][j][k]; 
                    }
                }else{ // above tropopause
                    t.x[i][j][k] = temp_tropopause[j];
                }
            }
        }
    }
//    fft_gaussian_filter(t, 5);
    logger() << "exit BC_Temperature: temperature max: " << (t.max()-1)*t_0 
        << std::endl << std::endl;
    if(debug){
        Array tmp = (t-1)*t_0;
        logger()<<"20180912: Exit BCT ... "<<std::endl;
        tmp.inspect("20180912: ");
    }
}


void cAtmosphereModel::read_NASA_temperature(const string &fn){
    // initial conditions for the Name_SurfaceTemperature_File at the sea surface
    ifstream ifs(fn);
    if(!ifs.is_open()) {
        cerr << "ERROR: unable to open SurfaceTemperature_File file: "<< fn << "\n";
        abort();
    }
    float lat, lon, temperature;
    for(int k=0; k < km && !ifs.eof(); k++){
        for(int j=0; j < jm; j++){
            ifs >> lat >> lon >> temperature;
            t.x[0][j][k] = temperature_NASA.y[j][k] = 
                (temperature + t_0) / t_0;
        }
    }
    // correction of surface temperature around 180°E
    int k_half = (km -1)/2;
    for(int j = 0; j < jm; j++){
        t.x[0][j][k_half] = (t.x[0][j][k_half + 1] 
            + t.x[0][j][k_half - 1]) / 2.;
        temperature_NASA.y[j][k_half] = (temperature_NASA.y[j][k_half + 1] +
            temperature_NASA.y[j][k_half - 1]) / 2.;
    }
}


void cAtmosphereModel::read_NASA_precipitation(const string &fn){
    // initial conditions for the Name_SurfacePrecipitation_File at the sea surface
    ifstream ifs(fn);
    if (!ifs.is_open()) {
        cerr << "ERROR: unable to open SurfacePrecipitation_File file: " << fn << "\n";
        abort();
    }
    double lat, lon, precipitation;
    for(int k=0; k < km && !ifs.eof(); k++){
        for(int j=0; j < jm; j++){
            ifs >> lat >> lon >> precipitation;
            precipitation_NASA.y[j][k] = precipitation;
        }
    }
}

void cAtmosphereModel::BC_Pressure(){
    float exp_pressure = g / (1.e-2 * gam * R_Air);
    // boundary condition of surface pressure given by surface temperature through gas equation
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            p_stat.x[0][j][k] = .01 * (r_air * R_Air 
                * t.x[0][j][k] * t_0);      // given in hPa
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 1; i < im; i++){
                double height = get_layer_height(i);
                p_stat.x[i][j][k] = pow(((t.x[0][j][k] 
                    * t_0 - gam * height * 1.e-2) /
                    (t.x[0][j][k] * t_0)), exp_pressure) 
                    * p_stat.x[0][j][k];
                // linear temperature distribution T = T0 - gam * height
                // current air pressure, step size in 500 m, from politropic formula in hPa
            }
        }
    }
}

void cAtmosphereModel::Latent_Heat(){
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
                step_p = ( get_layer_height(i+1) - get_layer_height(i) );
                step_m = ( get_layer_height(i) - get_layer_height(i-1) );
                step = step_p + step_m;
                e = .01 * c.x[i][j][k] * p_stat.x[i][j][k] / ep;  // water vapour pressure in Pa
                if(i > i_mount){
                    a = e / (R_WaterVapour * t.x[i][j][k] * t_0);  // absolute humidity in kg/m³
                    Q_Latent.x[i][j][k] = - lv * a 
                        * (c.x[i+1][j][k] - c.x[i-1][j][k]) / step;
                    Q_Latent_Ice = - ls * a * (ice.x[i+1][j][k] 
                        - ice.x[i-1][j][k]) / step;
                    Q_Sensible.x[i][j][k] = - coeff_sen * coeff_Q 
                        * (t.x[i+1][j][k] - t.x[i-1][j][k]) / step;  // sensible heat in [W/m2] from energy transport equation
                }
                if(i == i_mount){
                    a = e / (R_WaterVapour * t.x[i][j][k] * t_0);  // absolute humidity in kg/m³
                    Q_Latent.x[i][j][k] = - lv * a 
                        * (- 3. * c.x[i][j][k] + 4. * c.x[i+1][j][k] 
                        - c.x[i+2][j][k]) / step;
                    Q_Latent_Ice = - ls * a * (- 3. * ice.x[i][j][k] +
                        4. * ice.x[i+1][j][k] 
                        - ice.x[i+2][j][k]) / step;
                    // sensible heat in [W/m2] from energy transport equation
                    Q_Sensible.x[i][j][k] = - coeff_sen * coeff_Q 
                        * (- 3. * t.x[i][j][k] +
                        4. * t.x[i+1][j][k] - t.x[i+2][j][k]) / step;
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
            Q_Latent.x[i][j][0] = Q_Latent.x[i][j][km-1] = (b1+ b2) / 2.;
            b1 = c43 * Q_Sensible.x[i][j][1] 
                - c13 * Q_Sensible.x[i][j][2];
            b2 = c43 * Q_Sensible.x[i][j][km-2] 
                - c13 * Q_Sensible.x[i][j][km-3];
            Q_Sensible.x[i][j][0] = Q_Sensible.x[i][j][km-1] = (b1+ b2) / 2.;
        }
    }
}





//Tao, W.-K., Simpson, J., and McCumber, M.: 
//An Ice-Water Saturation Adjustment, American Meteorological Society, Notes and 
//Correspon-dence, Volume 1, 321–235, 1988, 1988.
//Ice_Water_Saturation_Adjustment, distribution of cloud ice and cloud water 
//dependent on water vapour amount and temperature
void cAtmosphereModel::Ice_Water_Saturation_Adjustment(){ 
    if(debug){
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
        assert(!c.has_nan());
        assert(!t.has_nan());
    }
    // constant coefficients for the adjustment of cloud water and cloud ice amount vice versa
    float t_00 = 236.15;
    float t_Celsius_0 = 0.; // in Celsius = 0 °C
    float exp_pressure = g / (1.e-2 * gam * R_Air);
    float q_v_hyp = 0., q_T = 0.;
    float dt_dim = L_atm / u_0 * dt;// dimensional time step of the system in s == 0.02 s
    float t_u = 0.;
    float T_it = 0.;
    float T = 0.;
    float d_t = 0.;
    float t_Celsius = 0.;
    float p_SL = 0.;
    float p_h = 0.;
    float height = 0.;
    float E_Rain = 0.;
    float E_Ice = 0.;
    float q_Rain = 0.;
    float q_Ice = 0.;
    float q_Rain_n = 0.;
    float q_v_b = 0.;
    float q_c_b = 0.;
    float q_i_b = 0.;
    float CND = 0.;
    float DEP = 0.;
    float d_q_v = 0.;
    float d_q_c = 0.;
    float d_q_i = 0.;
    cout.precision(5);
    cout.setf(ios::fixed);
    // setting water vapour, cloud water and cloud ice into the proper thermodynamic ratio based on the local temperatures
    // starting from a guessed parabolic temperature and water vapour distribution in north/south direction
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                t_u = t.x[i][j][k] * t_0; // in K
                t_Celsius = t_u - t_0; // in C
                p_SL = .01 * ( r_air * R_Air * t.x[0][j][k] * t_0 ); // given in hPa
                height = get_layer_height(i);
                if(i != 0)
                      p_h = pow((t_u - gam * height * 1.e-2) / t_u, 
                            exp_pressure) * p_SL;
                else  p_h = p_SL;
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Ice = hp * exp_func(t_u, 21.8746, 7.66); // saturation water vapour pressure for the ice phase in hPa
                q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg
                /** %%%%%%%%%%%%%%%%%%%%%%%%%%%     warm cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% **/
                // warm cloud phase in case water vapour is over-saturated
                if(t_Celsius >= t_Celsius_0){ //temperature above 0 Celsius
                    q_T = c.x[i][j][k] + cloud.x[i][j][k]; // total water content
                    if(!(q_T > q_Rain)){ //     subsaturated
                        c.x[i][j][k] = q_T; // total water amount as water vapour
                        cloud.x[i][j][k] = 0.; // no cloud water available
                        ice.x[i][j][k] = 0.; // no cloud ice available above 0 °C
                        T_it = t_u;
                    }else{ //     oversaturated
                        for(int iter_prec = 1; iter_prec <= 20; iter_prec++){ // iter_prec may be varied
                            if(i != 0)
                                p_h = pow((T_it - gam * height * 1.e-2) / t_u, 
                            exp_pressure) * p_SL;
                            else  p_h = p_SL;
                            T_it = t_u + lv / cp_l * ( c.x[i][j][k] - q_Rain );
                            E_Rain = hp * exp_func(T_it, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                            q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                            q_Rain = .5 * ( q_Rain_n + q_Rain );  // smoothing the iteration process
                            c.x[i][j][k] = q_Rain; // water vapour restricted to saturated water vapour amount
                            cloud.x[i][j][k] = q_T - c.x[i][j][k]; // cloud water amount
                            ice.x[i][j][k] = 0.; // no cloud ice available
                            q_T = c.x[i][j][k] + cloud.x[i][j][k];
                            if(c.x[i][j][k] < 0.)  c.x[i][j][k] = 0.;
                            if(cloud.x[i][j][k] < 0.)  cloud.x[i][j][k] = 0.;
                            if(fabs(q_Rain / q_Rain_n - 1.) < 1.e-3)  break;  
                            q_Rain_n = q_Rain;
                        }//end of for 
                    }//oversaturated
                    cn.x[i][j][k] = c.x[i][j][k];
                    cloudn.x[i][j][k] = cloud.x[i][j][k];
                    icen.x[i][j][k] = ice.x[i][j][k];
                    tn.x[i][j][k] = t.x[i][j][k] = T_it / t_0;
                } // end (t_Celsius > 0.)
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%     end  warm cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                // %%%%%%%%%%%%%%%%%%%%%%%%%%%     begin mixed cloud phase         %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else{ //temperature below 0 Celsius
                    if(t_Celsius > t_Celsius_0)  ice.x[i][j][k] = 0.;
                    q_v_b = c.x[i][j][k];
                    q_c_b = cloud.x[i][j][k];
                    q_i_b = ice.x[i][j][k];
                    q_T = q_v_b + q_c_b + q_i_b; // total water content
                    T = t_u; // in K
                    E_Rain = hp * exp_func(T, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                    E_Ice = hp * exp_func(T, 21.8746, 7.66); // saturation water vapour pressure for the ice phase in hPa
                    q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                    q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg
                    if(q_c_b > 0. && q_i_b > 0.)
                        q_v_hyp = ( q_c_b * q_Rain + q_i_b * q_Ice ) 
                                  / ( q_c_b + q_i_b );
                    if((q_c_b >= 0.) && (q_i_b == 0.))  q_v_hyp = q_Rain;
                    if((q_c_b == 0.) && (q_i_b >= 0.))  q_v_hyp = q_Ice;
                    // §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§
                    for(int iter_prec = 1; iter_prec <= 10; iter_prec++){ // iter_prec = 2 given by COSMO
                    // condensation == water vapor saturation for cloud water formation, deposition == ice crystal for cloud ice formation
                        CND = ( T - t_00 ) / ( t_0 - t_00 ); // t_00 = 236.15°C, t_0 = 273.15°C
                        DEP = ( t_0 - T ) / ( t_0 - t_00 );
                        if(T <= t_00){
                            CND = 0.;
                            DEP = 1.;
                        }
                        if(T >= t_0){
                            CND = 1.;
                            DEP = 0.;
                        }
                        d_q_v = q_v_hyp - q_v_b;  // changes in water vapour causing cloud water and cloud ice
                        d_q_v = .5 * d_q_v; // but why?
                        d_q_c = - d_q_v * CND;
                        d_q_i = - d_q_v * DEP;
                        d_t = ( lv * d_q_c + ls * d_q_i ) / cp_l; // in K, temperature changes
                        T = t_u + d_t; // in K
                        q_v_b = c.x[i][j][k] - d_q_v;  // new values
                        q_c_b = cloud.x[i][j][k] + d_q_c;
                        q_i_b = ice.x[i][j][k] + d_q_i;
                        if(q_v_b < 0.)  q_v_b = 0.;  // negative values excluded, when iteration starts
                        if(q_c_b < 0.)  q_c_b = 0.;
                        if(q_i_b < 0.)  q_i_b = 0.;
                        if(i != 0) p_h = pow(((T - gam * height * 1.e-2) / T), 
                                       exp_pressure) * p_SL; // given in hPa
                        else  p_h = p_SL;
                        E_Rain = hp * exp_func(T, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        E_Ice = hp * exp_func(T, 21.8746, 7.66); // saturation water vapour pressure for the ice phase in hPa
                        q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                        q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg
                        if((q_c_b > 0.) && (q_i_b > 0.))
                            q_v_hyp = ( q_c_b * q_Rain + q_i_b * q_Ice ) 
                                      / ( q_c_b + q_i_b );
                        if((q_c_b >= 0.) && (q_i_b == 0.))  q_v_hyp = q_Rain;
                        if((q_c_b == 0.) && (q_i_b >= 0.))  q_v_hyp = q_Ice;
                        if((q_c_b > 0.) && (q_i_b > 0.))
                        q_T = q_v_b + q_c_b + q_i_b; // total water content
                        if(iter_prec >= 3 && fabs(q_v_b / q_v_hyp - 1.) 
                            <= 1.e-3)  break;
                        q_v_b = .5 * ( q_v_hyp + q_v_b );  // has smoothing effect
                    } // iter_prec end
                    //** §§§§§§§§§§§§§   end iterations for mixed cloud phase     §§§§§§§**
                    // rate of condensating or evaporating water vapour to form cloud water in kg/(kg*s), 0.5 given by COSMO
                    S_c_c.x[i][j][k] = 0.5 * ( cloud.x[i][j][k] 
                        - q_c_b ) / dt_dim;
                    if(is_land(h, i, j, k))  S_c_c.x[i][j][k] = 0.;
                    cn.x[i][j][k] = c.x[i][j][k] = q_v_b;  // new values achieved after converged iterations
                    cloudn.x[i][j][k] = cloud.x[i][j][k] = q_c_b;
                    icen.x[i][j][k] = ice.x[i][j][k] = q_i_b;
                    tn.x[i][j][k] = t.x[i][j][k] = T / t_0;
                    if(t_Celsius > t_Celsius_0)  icen.x[i][j][k] = 
                                            ice.x[i][j][k] = 0.;
                } // temperature below 0 Celsius
            } // end i
        } // end j
    } // end k
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.)      c.x[i][j][k] = 0.;
                if(cloud.x[i][j][k] < 0.)  cloud.x[i][j][k] = 0.;
                if(ice.x[i][j][k] < 0.)    ice.x[i][j][k] = 0.;
//                if(is_land (h, i, j, k)){
//                    cloudn.x[i][j][k] = cloud.x[i][j][k] = 0.;
//                    icen.x[i][j][k] = ice.x[i][j][k] = 0.;
//                }
            }
        }
    }
    if(debug){
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
        assert(!c.has_nan());
        assert(!t.has_nan());
    }
//    fft_gaussian_filter(cloud, 5);
//    fft_gaussian_filter(ice, 5);
//    fft_gaussian_filter(c, 5);
}



// Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
// resulting the precipitation distribution formed of rain and snow
void cAtmosphereModel::Two_Category_Ice_Scheme(){   
    // constant coefficients for the transport of cloud water and cloud ice amount vice versa, 
    // rain and snow in the parameterization procedures
    float N_i_0 = 1.e2,  // in m-3
          m_i_0 = 1.e-12,  // in kg
          m_i_max = 1.e-9,  // in kg
          m_s_0 = 3.e-9,  // in kg
          c_i_dep = 1.3e-5,  // in m3/(kg*s)
          c_c_au = 4.e-4,  // in 1/s
          c_i_au = 1.e-3,  // in 1/s
          c_ac = .24,  // m2/kg
          c_rim = 18.6,  // m2/kg
          c_agg = 10.3,  // m2/kg
          c_i_cri = .24,  // m2
          c_r_cri = 3.2e-5,  // m2
          a_ev = 1.e-3,  // m2/kg
          b_ev = 5.9,  // m2*s/kg
          c_s_dep = 1.8e-2,  // m2/kg
          b_s_dep = 12.3,  // m2*s/kg
          c_s_melt = 8.43e-5,  // (m2*s)/(K*kg)
          b_s_melt = 12.05,  // m2*s/kg
          a_s_melt = 2.31e3, // K/(kg/kg)
          c_r_frz = 3.75e-2,  // (m2*s)/(K*kg)
          t_nuc = 267.15,  // in K    -6 °C
          t_d = 248.15,  // in K    -25 °C
          t_hn = 236.15,  // in K    -40 °C
          t_r_frz = 271.15;  // in K    -2 °C
    double exp_pressure = g / (1.e-2 * gam * R_Air);
    float dt_snow_dim = 417.,  // dt_snow_dim is the time  in 417 s to pass dr = 400 m, 400 m / 417 s = .96 m/s fallout velocity
          dt_rain_dim = 250.;  // dt_rain_dim is the time  in 250 s to pass dr = 400 m, 400 m / 250 s = 1.6 m/s fallout velocity
    double coeff_snow = .01;
    double m_i = m_i_max;  
    float p_h, N_i, S_nuc, S_c_frz, S_i_dep=0, S_c_au, S_i_au, S_d_au, 
        S_ac, S_rim, S_shed;
    float S_agg, S_i_cri, S_r_cri, S_ev, S_s_dep, S_i_melt, S_s_melt, S_r_frz;
    // rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
    // the choosen scheme is a Two Category Ice Scheme
    // besides the transport equation for the water vapour exists two equations for the cloud water and the cloud ice transport
    // since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.)  c.x[i][j][k] = 0.;
                if(cloud.x[i][j][k] < 0.)  cloud.x[i][j][k] = 0.;
                if(ice.x[i][j][k] < 0.)  ice.x[i][j][k] = 0.;
                if(P_rain.x[i][j][k] < 0.)  P_rain.x[i][j][k] = 0.;
                if(P_snow.x[i][j][k] < 0.)  P_snow.x[i][j][k] = 0.;
            }
        }
    }
    /******************* initial values for rain and snow calculation *********************/
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            P_rain.x[im-1][j][k] = 0.;
            P_snow.x[im-1][j][k] = 0.;
            S_r.x[im-1][j][k] = 0.;
            S_s.x[im-1][j][k] = 0.;
            for(int i = im-2; i >= 0; i--){
                float t_u = t.x[i][j][k] * t_0;
                if((is_land(h, i, j, k)) && (is_land(h, i+1, j, k))){
                    S_r.x[i][j][k] = 0.;
                    S_s.x[i][j][k] = 0.;
                }else{
                    if(t_u >= t_0){
                        if(cloud.x[i][j][k] > 0.)
                            S_c_au = c_c_au * cloud.x[i][j][k];// cloud water to rain, cloud droplet collection in kg/(kg*s)
                        S_i_au = 0.;
                    }else{ // temperature < 0 
                        S_c_au = 0.;
                        if(ice.x[i][j][k] > 0.)
                            S_i_au = c_i_au * ice.x[i][j][k];  // cloud ice to snow, cloud ice crystal aggregation
                    }
                    S_r.x[i][j][k] = S_c_au;  // in kg / (kg * s)
                    S_s.x[i][j][k] = S_i_au;
                }
                if(P_rain.x[i+1][j][k] < 0.)  P_rain.x[i+1][j][k] = 0.;
                if(P_snow.x[i+1][j][k] < 0.)  P_snow.x[i+1][j][k] = 0.;
                float p_SL = .01 * (r_air * R_Air * t.x[0][j][k] * t_0); // given in hPa
                double height = get_layer_height(i);
                float p_h;
                if(i != 0)  
                    p_h = pow(((t_u - gam * height * 1.e-2) / t_u), 
                        exp_pressure) * p_SL;
                else  p_h = p_SL;
                float r_dry = 100. * p_h / (R_Air * t_u);  // density of dry air in kg/m³
                float r_humid = r_dry 
                    / (1. + ( R_WaterVapour / R_Air - 1. ) * c.x[i][j][k] 
                    - cloud.x[i][j][k] - ice.x[i][j][k]);                
                double step = get_layer_height(i+1) - get_layer_height(i);
                float t_Celsius = t_u - t_0;
                P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                    + coeff_Precipitation * r_humid * S_r.x[i+1][j][k] 
                    * step;  // in kg / (m2 * s) == mm/s
                if(t_Celsius < 0.)  P_snow.x[i][j][k] = P_snow.x[i+1][j][k]
                    + coeff_Precipitation * r_humid * coeff_snow
                    * S_s.x[i+1][j][k] * step; 
                if(P_rain.x[i][j][k] < 0.)  P_rain.x[i][j][k] = 0.;
                if(P_snow.x[i][j][k] < 0.)  P_snow.x[i][j][k] = 0.;
            }
        }
    }
    /******************* main part for rain and snow calculation *********************/
    if(true){
        for(int iter_prec = 1; iter_prec < 2; iter_prec++){  // no iterations made
            for(int k = 0; k < km; k++){
                for(int j = 0; j < jm; j++){
                    P_rain.x[im-1][j][k] = 0.;
                    P_snow.x[im-1][j][k] = 0.;
                    S_r.x[im-1][j][k] = 0.;
                    S_s.x[im-1][j][k] = 0.;
                    for(int i = im-2; i >= 0; i--){
                        float t_u = t.x[i][j][k] * t_0;
                        float t_Celsius = t_u - t_0;
                        float p_SL =  .01 * (r_air * R_Air * t.x[0][j][k] * t_0); // given in hPa
                        float height = get_layer_height(i);
                        if(i != 0)  
                            p_h = pow(((t_u - gam * height * 1.e-2) / t_u), 
                                exp_pressure) * p_SL;
                        else  p_h = p_SL;
                        float r_dry = 100. * p_h / (R_Air * t_u);  // density of dry air in kg/m³
                        float r_humid = r_dry 
                            / (1. + ( R_WaterVapour / R_Air - 1. ) * c.x[i][j][k] 
                            - cloud.x[i][j][k] - ice.x[i][j][k]);                
                        float E_Rain = hp * exp_func(t_u, 17.2694, 35.86);  // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        float E_Ice = hp * exp_func(t_u, 21.8746, 7.66);  // saturation water vapour pressure for the ice phase in hPa
                        float q_Rain = ep * E_Rain / (p_h - E_Rain);  // water vapour amount at saturation with water formation in kg/kg
                        float q_Ice  = ep * E_Ice / (p_h - E_Ice);  // water vapour amount at saturation with ice formation in kg/kg
                        dt_rain_dim = ( get_layer_height(i+1) 
                            - height ) / 1.6; 
                        // adjusted rain fall time step by fixed velocities == 1.6 m/s by variable local step size
                        dt_snow_dim = ( get_layer_height(i+1) 
                            - height ) / 0.96; 
                        // adjusted snow fall time step by fixed velocities == 0.96 m/s by variable local step size
                        // ice and snow average siz 
                        if(!(t_u > t_0)){  
                            N_i = N_i_0 * exp(.2 * (t_0 - t_u));
                            m_i = r_humid * ice.x[i][j][k] / N_i;
                            if(m_i > m_i_max) { m_i = m_i_max; }
                            if(m_i < m_i_0) { m_i = m_i_0; }
                        }
                        // nucleation and depositional growth of cloud ice
                        if(ice.x[i][j][k] == 0.){
                            if(((t_u < t_d) && (c.x[i][j][k] >= q_Ice))
                                || (((t_d <= t_u) && (t_u <= t_nuc)) 
                                && (c.x[i][j][k] >= q_Rain)))
                                S_nuc = m_i_0 / (r_humid * dt_snow_dim) * N_i;  // nucleation of cloud ice, < I >
                        }
                        else  S_nuc = 0.;
                        if((t_u < t_hn) && (cloud.x[i][j][k] > 0.))
                            S_c_frz = cloud.x[i][j][k] / dt_rain_dim;  //nucleation of cloud ice due to freezing of cloud water, < II >
                        else  S_c_frz = 0.;
                        if(!(t_Celsius > 0.)){  //temperature <= 0
                            if(c.x[i][j][k] > q_Ice){  // supersaturation
                                S_i_dep = c_i_dep * N_i * pow(m_i, (1. / 3.)) 
                                    * (c.x[i][j][k] - q_Ice);  // supersaturation, < III >
                            }else if(- ice.x[i][j][k]  > c.x[i][j][k] - q_Ice) 
                                S_i_dep = - ice.x[i][j][k] / dt_snow_dim; // subsaturation, < III >
                        }
                        else  S_i_dep = 0.; //temperature > 0
                        // autoconversion processes
                        if((t_u >= t_0) && (cloud.x[i][j][k] > 0.))
                            S_c_au = c_c_au * cloud.x[i][j][k];  // cloud water to rain, cloud droplet collection, < IV >
                        else  S_c_au = 0.;
                        if((t_u <= t_0) && (ice.x[i][j][k] > 0.))
                            S_i_au = c_i_au * ice.x[i][j][k];  // cloud ice to snow, cloud ice crystal aggregation, < V >
                        else  S_i_au = 0.;
                        if(t_u <= t_0)
                            S_d_au = S_i_dep / (1.5 * (pow((m_s_0 / m_i), 
                                (2. / 3.)) - 1.));  // autoconversion due to depositional growth of cloud ice, < VI >
                        else  S_d_au = 0.;
                        // collection mechanism
                        if(t_u > t_0)  
                            S_ac = c_ac * cloud.x[i][j][k] 
                                * pow(P_rain.x[i][j][k], (7. / 9.));
                            // accreation rate from depletion of cloud water due to collection by all rain drops, < VII >
                        else  S_ac = 0.;  // accreation rate from depletion of cloud water due to collection by all rain drops
                        if(t_u < t_0)  
                            S_rim = c_rim * cloud.x[i][j][k] * P_snow.x[i][j][k];
                        else  S_rim = 0.;  // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >
                                           // by falling snow particles
                        if(t_u >= t_0)  
                            S_shed = c_rim * cloud.x[i][j][k] * P_snow.x[i][j][k];
                        else  S_shed = 0.;  // rate of water shed by melting wet snow particles, < IX >
                                            // collecting cloud droplets to produce rain
                        if(t_u <= t_0){
                            S_agg = c_agg * ice.x[i][j][k] * P_snow.x[i][j][k];  // collection of cloud ice by snow particles, < X >
                            S_i_cri = c_i_cri * ice.x[i][j][k] 
                                * pow(P_rain.x[i][j][k], (7. / 9.));
                                // decrease in cloud ice mass due to collision/coalescense interaction with raindrops, < XI >
                            S_r_cri = c_r_cri * ice.x[i][j][k] / m_i 
                                * pow(P_rain.x[i][j][k], (13. / 9.));
                                // decrease of rainwater due to freezing resulting from collection of ice crystals, < XII >
                        }else{
                            S_agg = 0.;
                            S_i_cri = 0.;
                            S_r_cri = 0.;
                        }
                        // diffusional growth of rain and snow
                        if(t_u < t_0){  // temperature below zero
                            S_s_dep = c_s_dep * (1. + b_s_dep 
                                * pow(P_snow.x[i][j][k], (5. / 26.))) * 
                                (c.x[i][j][k] - q_Ice) 
                                * pow(P_snow.x[i][j][k], (8. / 13.));
                             // deposition/sublimation of snow, < XIV >
                            S_ev = 0.;
                        }else{  
                            S_ev = a_ev * (1. + b_ev 
                                * pow(P_rain.x[i][j][k], (1. / 6.))) *
                                (q_Rain - c.x[i][j][k]) 
                                * pow(P_rain.x[i][j][k], (4. / 9.));
                            // evaporation of rain due to water vapour diffusion, < XIII >
                            S_s_dep = 0.;
                        }
                        // melting and freezing
                        S_i_melt = 0.;
                        S_s_melt = 0.;
                        if(t_u > t_0){ // temperature above zero
                            if(ice.x[i][j][k] > 0.)
                                S_i_melt = ice.x[i][j][k] / dt_snow_dim; // cloud ice particles melting to cloud water, < XV >
                            float p_t_in = pow(((t_0 - gam * height * 1.e-2) 
                                / t_0), exp_pressure) * p_SL;  // given in hPa
                            float E_Rain_t_in = hp * exp_func(t_0, 17.2694, 35.86);
                                // saturation water vapour pressure for the water phase at t = 0°C in hPa
                            float q_Rain_t_in = ep * E_Rain_t_in / (p_t_in 
                                - E_Rain_t_in);
                                // water vapour amount at saturation with water formation in kg/kg
                            S_s_melt = c_s_melt * (1. + b_s_melt 
                                * pow(P_snow.x[i][j][k], (5. / 26.))) *
                                ((t_u - t_0) + a_s_melt * (c.x[i][j][k] 
                                - q_Rain_t_in)) * pow(P_snow.x[i][j][k], 
                                (8. / 13.));  // melting rate of snow to form rain, < XVI >
                                // arbitrary factor 50 adjusts the average precipitation in mm/d
                        }
                        if(t_r_frz - t_u > 0.)
                            S_r_frz = c_r_frz * pow((t_r_frz - t_u), 
                                (3. / 2.)) * pow(P_rain.x[i][j][k], (3. / 2.));
                                // immersion freezing and contact nucleation, < XVII >
                        else  S_r_frz = 0.;

// switching on, causes: 
//                        S_s_melt=0; // causes 
//                        S_r_frz=0; // causes 
//                        S_r_cri=0; // causes 
//                        S_i_cri=0;// causes 

//                        S_agg=0;// causes no snow around pole regions
//                        S_rim=0; // causes less snow around pole regions

//                        S_i_au=0; // causes less snow around pole regions
//                        S_s_dep=0; // causes less snow in near pole regions
//                        S_d_au=0; //  causes considerably more snow, TODO
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
                        if((is_land(h, i, j, k)) && (is_land(h, i+1, j, k))){
                            S_c_c.x[i][j][k] = 0.;
                            S_v.x[i][j][k] = 0.;
                            S_c.x[i][j][k] = 0.;
                            S_i.x[i][j][k] = 0.;
                            S_r.x[i][j][k] = 0.;
                            S_s.x[i][j][k] = 0.;
                        }
                        // rain and snow integration
                        if(P_rain.x[i+1][j][k] < 0.)  P_rain.x[i+1][j][k] = 0.;
                        if(P_snow.x[i+1][j][k] < 0.)  P_snow.x[i+1][j][k] = 0.;
                        double step = get_layer_height(i+1) - get_layer_height(i);
                        P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                             + coeff_Precipitation * r_humid 
                             * S_r.x[i+1][j][k] * step;  // in kg / (m2 * s) == mm/s 
                        if(t_Celsius < 0.)  P_snow.x[i][j][k] = P_snow.x[i+1][j][k]
                             + coeff_Precipitation * r_humid * coeff_snow 
                             * S_s.x[i+1][j][k] * step; 
                        else  P_snow.x[i][j][k] = 0.;
                        if(P_rain.x[i][j][k] < 0.)  P_rain.x[i][j][k] = 0.;
                        if(P_snow.x[i][j][k] < 0.)  P_snow.x[i][j][k] = 0.;
                    }  // end i RainSnow
                }  // end j
            }  // end k
        }  // end iter_prec
    }  // end if true
}





void cAtmosphereModel::Value_Limitation_Atm(){
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            if(Precipitation.y[j][k] >= 25.)  Precipitation.y[j][k] = 25.;
            if(Precipitation.y[j][k] <= 0)  Precipitation.y[j][k] = 0.;
            for(int i = 0; i < im; i++){
/*
                if(u.x[i][j][k] >= .106)  u.x[i][j][k] = .106;
                if(u.x[i][j][k] <= - .106)  u.x[i][j][k] = - .106;
                if(v.x[i][j][k] >= 1.125)  v.x[i][j][k] = 1.125;
                if(v.x[i][j][k] <= - 1.125)  v.x[i][j][k] = - 1.125;
                if(w.x[i][j][k] >= 3.5)  w.x[i][j][k] = 3.5;
                if(w.x[i][j][k] <= - 1.469)  w.x[i][j][k] = - 1.469;
                if(t.x[i][j][k] >= 1.165)  t.x[i][j][k] = 1.165;  // == 45 °C
                if(t.x[i][j][k] <= - .78)  t.x[i][j][k] = - .78;  // == 59.82 °C
                if(c.x[i][j][k] >= .035)  c.x[i][j][k] = .035;
                if(c.x[i][j][k] < 0.)  c.x[i][j][k] = 0.;
                if(cloud.x[i][j][k] >= .02)  cloud.x[i][j][k] = .02;
                if(cloud.x[i][j][k] < 0.)  cloud.x[i][j][k] = 0.;
                if(ice.x[i][j][k] >= .01)  ice.x[i][j][k] = .01;
                if(ice.x[i][j][k] < 0.)  ice.x[i][j][k] = 0.;
                if(P_rain.x[i][j][k] >= 10.)  P_rain.x[i][j][k] = 10.;
                if(P_rain.x[i][j][k] < 0.)  P_rain.x[i][j][k] = 0.;
                if(P_snow.x[i][j][k] >= 1.)  P_snow.x[i][j][k] = 1.;
                if(P_snow.x[i][j][k] < 0.)  P_snow.x[i][j][k] = 0.;
                if(co2.x[i][j][k] >= 5.36)  co2.x[i][j][k] = 5.36;
                if(co2.x[i][j][k] <= 1.)  co2.x[i][j][k] = 1.;
*/
                if(is_land(h, i, j, k)){
                    u.x[i][j][k] = 0.;
                    v.x[i][j][k] = 0.;
                    w.x[i][j][k] = 0.;
//                    t.x[i][j][k] = 1.;  // = 273.15 K
//                    c.x[i][j][k] = 0.;
                    cloud.x[i][j][k] = 0.;
                    ice.x[i][j][k] = 0.;
//                    co2.x[i][j][k] = 1.;  // = 280 ppm
                    p_dyn.x[i][j][k] = 0.;
                }
            }
        }
    }
}


void cAtmosphereModel::IC_v_w_WestEastCoast(){
// initial conditions for v and w velocity components at the sea surface close to east or west coasts
// reversal of v velocity component between north and south equatorial current ommitted at respectively 10°
// w component not changed in sign
// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included
// east coast
    int i_max = 20;
    double d_i_max = get_layer_height(i_max);
    double d_i = 0.;
    double v_neg = 0.;
    int k_mid = 30;
    int k_beg = 1;
    for(int j = 0; j < jm; j++){
        for(int k = 1; k < km-1; k++){
            if((is_air(h, 0, j, k))&&(is_land(h, 0, j, k-1))){
                while(k >= k_beg){
                    if(is_land(h, 0, j, k+1))  break;
                    if(is_air(h, 0, j, k-2))  break;
                    if(is_air(h, 0, j, k-3))  break;
                    v_neg = - v.x[0][j][k+k_mid];
                    for(int l = 0; l <= k_mid; l++){
                        v.x[0][j][k+l] = ( v.x[0][j][k+k_mid] - v_neg ) 
                            / (double)k_mid * (double)l + v_neg;
                        for(int i = 1; i <= i_max; i++){
                            d_i = get_layer_height(i);
                            v.x[i][j][k+l] = ( v.x[i_max][j][k+l] - v.x[0][j][k+l] ) 
                                / d_i_max * d_i + v.x[0][j][k+l];
                        }
                        if(is_land(h, 0, j, k+l))  v.x[0][j][k+l] = 0.;
                    }
                    for(int l = 0; l <= k_mid; l++){
                        w.x[0][j][k+l] = w.x[0][j][k+k_mid] 
                            / (double)k_mid * (double)l;
                        for(int i = 1; i <= i_max; i++){
                            d_i = get_layer_height(i);
                            w.x[i][j][k+l] = ( w.x[i_max][j][k+l] - w.x[0][j][k+l] ) 
                                / d_i_max * d_i + w.x[0][j][k+l];
                        }
                        if(is_land(h, 0, j, k+l))  w.x[0][j][k+l] = 0.;
                    }
                    k_beg = k + k_mid;
                    break;
                } // end while
            } // end if
        } // end k
        k_beg = 1;
    } // end j
// west coast
    k_mid = 10;
    k_beg = k_mid;
    for(int j = 0; j < jm; j++){
        for(int k = k_mid; k < km-1; k++){
            if((is_air(h, 0, j, k-1))&&(is_land(h, 0, j, k))){
                while(k >= k_beg){
                    if(is_air(h, 0, j, k+1))  break;
/*
//                    v_neg = - v.x[0][j][k+k_mid];
                    v_neg = v.x[0][j][k+k_mid];
                    for(int l = k_mid; l >= 0; l--){
                        v.x[0][j][k-l] = ( v.x[0][j][k-k_mid] - v_neg ) 
                            / (double)k_mid * (double)l + v_neg;
                        for(int i = 1; i <= i_max; i++){
                            d_i = get_layer_height(i);
                            v.x[i][j][k-l] = ( v.x[i_max][j][k-l] - v.x[0][j][k-l] ) 
                                / d_i_max * d_i + v.x[0][j][k-l];
                        }
                        if(is_land(h, 0, j, k-l))  v.x[0][j][k-l] = 0.;
                    }
*/
                    for(int l = k_mid; l >= 0; l--){
                        w.x[0][j][k-l] = w.x[0][j][k-k_mid] 
                            / (double)k_mid * (double)l;
                        for(int i = 1; i <= i_max; i++){
                            d_i = get_layer_height(i);
                            w.x[i][j][k-l] = ( w.x[i_max][j][k-l] - w.x[0][j][k-l] ) 
                                / d_i_max * d_i + w.x[0][j][k-l];
                        }
                        if(is_land(h, 0, j, k-l))  w.x[0][j][k-l] = 0.;
                    }
                    k_beg = k - k_mid;
                    break;
                } // end while
            } // end if
        } // end k
        k_beg = k_mid;
    } // end j
}








void cAtmosphereModel::WaterVapourEvaporation(){ 
// preparations for water vapour increase due to the differences between evaporation and precipitation 
// procedure given in Rui Xin Huang, Ocean Circulation, p. 165
// amount of additional water vapour by the difference of evaporation and precipitation is negligible but functional
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.)  c.x[i][j][k] = 0.;
                if(P_rain.x[i][j][k] < 0.)  P_rain.x[i][j][k] = 0.;
                if(P_snow.x[i][j][k] < 0.)  P_snow.x[i][j][k] = 0.;
            }
        }
    }
    // surface values of Evaporation, Condensation, Water, Water_super, IceAir, precipitable_water only for radial printout
    precipitable_water.y[0][0] = 0.;
    precipitable_water.y[0][0] = 0.;
    co2_total.y[0][0] = 0.;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            co2_total.y[j][k] = co2.x[0][j][k];
            for(int i = 0; i < im; i++){
                float e = 100. * c.x[i][j][k] * p_stat.x[i][j][k] / ep;  // water vapour pressure in Pa
                float a = e / (R_WaterVapour * t.x[i][j][k] * t_0);  // absolute humidity in kg/m³
                double step = get_layer_height(i+1) 
                    - get_layer_height(i);
                precipitable_water.y[j][k] +=  a * step;
                 // mass of water in kg/m²
                // precipitable_water mass in 1 kg/m² compares to 1 mm height, with water density kg/ (m² * mm)
            }
        }
    }
    double coeff_prec = 86400.;  // dimensions see below
    // surface values of precipitation and precipitable water
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            if((j == 0) && (k == 0)) P_snow.x[0][0][0] = 0.;
            Precipitation.y[j][k] = coeff_prec * (P_rain.x[0][j][k] 
                + P_snow.x[0][j][k]); // in mm/d
            // 60 s * 60 min * 24 h = 86400 s == 1 d
            // Precipitation, P_rain and P_snow in kg/ (m² * s) == mm/s
            // Precipitation in 86400. * kg/ (m² * d) = 86400 * mm/d
            // kg/ (m² * s) == mm/s (Kraus, p. 94)
        }
    }
    float e, E_Rain, E_Ice, sat_deficit;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            int i = get_surface_layer(j,k);
            float t_u = t.x[i][j][k] * t_0;
            float t_Celsius = t_u - t_0;
            if(t_Celsius >= 0.){
                if(is_land(h, i, j, k)){
                    t_u = t.x[i][j][k] * t_0;
                    t_Celsius = t_u - t_0;
                    e = c.x[i][j][k] * p_stat.x[i][j][k] / ep;  // water vapour pressure in hPa
                    E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                    sat_deficit = ( E_Rain - e );  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = 
                        C_Dalton(coeff_Dalton, u_0, v.x[i+1][j][k], w.x[i+1][j][k]) 
                        * sat_deficit * 24.; // since at the suface velocity is 0, one grid point added
                } // simplified formula for Evaporation by Dalton law dependent on surface water velocity in kg/(m²*d) = mm/d
                if(is_water(h, 0, j, k)){
                    t_Celsius = ( t.x[0][j][k] + 0.007322 ) * t_0 - t_0; // assumption that ocean surface temperature is 2°C higher
                    e = c.x[0][j][k] * p_stat.x[0][j][k] / ep;  // water vapour pressure in hPa
                    t_u = ( t.x[0][j][k] + 0.007322 ) * t_0;
                    E_Rain = hp * exp_func(t_u, 17.2694, 35.86);
                    sat_deficit = ( E_Rain - e );  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = 
                        C_Dalton(coeff_Dalton, u_0, v.x[0][j][k], w.x[0][j][k]) 
                        * sat_deficit * 24.;
                }
            }
            if(t_Celsius <= 0.){
                if(is_land(h, i, j, k)){
                    t_u = t.x[i][j][k] * t_0;
                    t_Celsius = t_u - t_0;
                    e = c.x[i][j][k] * p_stat.x[i][j][k] / ep;  // water vapour pressure in hPa
                    E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                    sat_deficit = ( E_Ice - e );  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = 
                        C_Dalton(coeff_Dalton, u_0, v.x[i+1][j][k], w.x[i+1][j][k]) 
                        * sat_deficit * 24.; // since at the suface velocity is 0, one grid point added
                }
                if(is_water(h, 0, j, k)){
                    e = c.x[0][j][k] * p_stat.x[0][j][k] / ep;  // water vapour pressure in hPa
                    t_u = ( t.x[0][j][k] + 0.007322 ) * t_0;
                    t_Celsius = t_u - t_0; // assumption that ocean surface temperature is 2°C higher
                    E_Ice = hp * exp_func(t_u, 21.8746, 7.66);
                    sat_deficit = ( E_Ice - e );  // saturation deficit in hPa
                    Evaporation_Dalton.y[j][k] = 
                        C_Dalton(coeff_Dalton, u_0, v.x[0][j][k], w.x[0][j][k]) 
                        * sat_deficit * 24.;
                }
            }
            if(Evaporation_Dalton.y[j][k] <= 0.)  
                Evaporation_Dalton.y[j][k] = 0.;
        }
    }
    double coeff_vapour = 1.1574e-8 * L_atm / ( c_0 * u_0 );  // 1.1574-8 is the conversion from (Evap-Prec) in mm/d to mm/s
    double zeta = 3.715;
    double rm = rad.z[0];
    double exp_rm = 1. / exp(zeta * rm);
    double evap_precip = 0.;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            if(iter_cnt_3d == 1)  c_fix.y[j][k] = c.x[0][j][k];
            double vapour_surface_n = 0.;
            double vapour_surface = 0.;
            for(int iter_prec = 1; iter_prec <= 20; iter_prec++){ // iter_prec may be varied
                evap_precip = Evaporation_Dalton.y[j][k] - Precipitation.y[j][k];
//                vapour_surface = - ( - 3. * c.x[0][j][k] + 4. * c.x[1][j][k] 
//                    - c.x[2][j][k] ) / ( 2. * dr * exp_rm ) // 1. order derivative, 2. order accurate
//                    * ( 1. + 2. * c.x[0][j][k] ) * evap_precip; 
                vapour_surface = - ( c.x[1][j][k] - c.x[0][j][k] ) // 1. order derivative, 1. order accurate
                    / ( dr * exp_rm ) * ( 1. + 2. * c.x[0][j][k] ) 
                    * evap_precip;
                if(iter_prec == 1)  vapour_surface_n = .9 * vapour_surface;
                vapour_evaporation.y[j][k] = coeff_vapour * vapour_surface;
                if(is_land(h, 0, j, k))
                    vapour_evaporation.y[j][k] = 0.;
                c.x[0][j][k] = c_fix.y[j][k] + vapour_evaporation.y[j][k];
/*
    cout.precision(5);
    cout.setf(ios::fixed);
    if((j == 90) && (k == 180)) cout << "  it = " << iter_prec 
        << "  vap_evap = " << vapour_evaporation.y[j][k] * 1000. 
        << "  coeff_vap = " << coeff_vapour  << "  vap_surf = " << vapour_surface 
        << "  vap_surf_n = " << vapour_surface_n << "  c_fix = " << c_fix.y[j][k] * 1000. 
        << "  c = " << c.x[0][j][k] * 1000. << "  Evap-Prec = " << evap_precip 
        << "  Evap = " << Evaporation_Dalton.y[j][k] 
        << "  Prec = " << Precipitation.y[j][k] << "  c_grad_1 = " 
        << ( c.x[1][j][k] - c.x[0][j][k] ) / ( dr * exp_rm ) << "  c_grad_2 = " 
        << ( - 3. * c.x[0][j][k] + 4. * c.x[1][j][k] - c.x[2][j][k] ) 
        / ( 2. * dr * exp_rm ) << endl;
*/
                for(int i = 0; i < im; i++){
                    if(i <= im-1){
                        if(i > 0){
                            double x = get_layer_height(i) 
                                / get_layer_height(im-1); 
                            c.x[i][j][k] = parabola_interp(c_tropopause, 
                                c.x[0][j][k], x); 
                        }
                    }
                } // end i
                if(iter_prec >= 5 && fabs(vapour_surface / vapour_surface_n - 1.) 
                    <= 1.e-3)  break;  
                vapour_surface_n = vapour_surface;
            } // iter_prec 
        } // end j
    } // end k
}


void cAtmosphereModel::Pressure_Limitation_Atm(){
    // class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(p_dyn.x[i][j][k] >= .25)  p_dyn.x[i][j][k] = .25;
                if(p_dyn.x[i][j][k] <= - .25)  p_dyn.x[i][j][k] = - .25;

                if(is_land(h, i, j, k)){
                    p_dyn.x[i][j][k] = 0.;
                }
                p_dynn.x[i][j][k] = p_dyn.x[i][j][k];
            }
        }
    }
}

double get_pole_temperature(int Ma, int Ma_1, int Ma_2, double t_1, double t_2){
    return (t_2 - t_1) / (double) (Ma_2 - Ma_1) * (double) (Ma - Ma_1) + t_1;
}

double get_pole_temperature(int Ma, const std::map<float, float> &pole_temp_map){
    assert(pole_temp_map.size()>1);
    std::pair<int, double> up = *pole_temp_map.begin(), bottom = *++pole_temp_map.begin();
    // when Ma out of boundary
    if(Ma <= pole_temp_map.begin()->first){
        return pole_temp_map.begin()->second; 
    }else if(Ma > (--pole_temp_map.end())->first){
        return (--pole_temp_map.end())->second;
    }
    for(const auto& pair : pole_temp_map){
        if(pair.first>=Ma){
            bottom = pair;
            break;
        }else{
            up = pair;
        }
    }
    return get_pole_temperature(Ma, up.first, bottom.first, up.second, bottom.second);
}

