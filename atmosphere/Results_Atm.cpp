/*
 * Atmosphere General Circulation Modell(AGCM)applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "Utils.h"
#include "cAtmosphereModel.h"

using namespace std;
using namespace AtomUtils;

void cAtmosphereModel::print_min_max_atm(){
    cout << endl << " flow properties: " << endl << endl;
    searchMinMax_3D(" max temperature ", " min temperature ", 
        " deg", t, 273.15, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max u-component ", " min u-component ", 
        "m/s", u, u_0);
    searchMinMax_3D(" max v-component ", " min v-component ", 
        "m/s", v, u_0);
    searchMinMax_3D(" max w-component ", " min w-component ", 
        "m/s", w, u_0);
    searchMinMax_3D(" max pressure dynamic ", " min pressure dynamic ", 
        "hPa", p_dyn, r_air * u_0 * u_0 * 1e-2);
    searchMinMax_3D(" max pressure static ", " min pressure static ", 
        "hPa", p_stat, 1.0);
    searchMinMax_3D(" max density dry air ", " min density dry air ", 
        "kg/m3", r_dry, 1.0);
    searchMinMax_3D(" max density wet air ", " min density wet air ", 
        "kg/m3", r_humid, 1.0);
    cout << endl << " energies: " << endl << endl;
    searchMinMax_3D(" max radiation ",  " min radiation ", 
        "W/m2", radiation, 1.0);
    searchMinMax_3D(" max sensible heat ", " min sensible heat ", 
        "W/m2", Q_Sensible, 1.0);
    searchMinMax_3D(" max latent heat ", " min latent heat ", 
        "W/m2", Q_Latent, 1.0);
    cout << endl << " cloud thermodynamics: " << endl << endl;
    searchMinMax_3D(" max water vapour ",  " min water vapour ", 
        "g/kg", c, 1000.0);
    searchMinMax_3D(" max cloud water ", " min cloud water ", 
        "g/kg", cloud, 1000.0);
    searchMinMax_3D(" max cloud ice ", " min cloud ice ", 
        "g/kg", ice, 1000.0);
    cout << endl << " precipitation: " << endl << endl;
    searchMinMax_3D(" max rain ", " min rain ", "mm/d", P_rain, 8.64e4);
    searchMinMax_3D(" max snow ", " min snow ", "mm/d", P_snow, 8.64e4);
    searchMinMax_3D(" max conv ", " min conv ", "mm/d", P_conv, 8.64e4);
    cout << endl << " moist convection: " << endl;
    cout << endl << " --- up/downdraft: " << endl << endl;
    searchMinMax_3D(" max E_u ", " min E_u ", "kg/m3s", E_u, 1.0);
    searchMinMax_3D(" max D_u ", " min D_u ", "kg/m3s", D_u, 1.0);
    searchMinMax_3D(" max M_u ", " min M_u ", "kg/m2s", M_u, 1.0);
    searchMinMax_3D(" max E_d ", " min E_d ", "kg/m3s", E_d, 1.0);
    searchMinMax_3D(" max D_d ", " min D_d ", "kg/m3s", D_d, 1.0);
    searchMinMax_3D(" max M_d ", " min M_d ", "kg/m2s", M_d, 1.0);
    cout << endl << " --- source terms: " << endl << endl;
    searchMinMax_3D(" max MC_t ", " min MC_t ", "K/s", MC_t, 1.0);
    searchMinMax_3D(" max MC_q ", " min MC_q ", "g/kgs", MC_q, 1000.0);
    searchMinMax_3D(" max MC_v ", " min MC_v ", "m/ss", MC_v, u_0);
    searchMinMax_3D(" max MC_w ", " min MC_w ", "m/ss", MC_w, u_0);
    cout << endl << " --- velocity components in the up- and downdraft: " 
        << endl << endl;
    searchMinMax_3D(" max u_u ", " min u_u ", "m/s", u_u, 1.0);
    searchMinMax_3D(" max v_u ", " min v_u ", "m/s", v_u, 1.0);
    searchMinMax_3D(" max w_u ", " min w_u ", "m/s", w_u, 1.0);
    searchMinMax_3D(" max u_d ", " min u_d ", "m/s", u_d, 1.0);
    searchMinMax_3D(" max v_d ", " min v_d ", "m/s", v_d, 1.0);
    searchMinMax_3D(" max w_d ", " min w_u ", "m/s", w_d, 1.0);
    cout << endl 
    << " --- condensation and evaporation terms in the up- and downdraft: " 
        << endl << endl;
    searchMinMax_3D(" max c_u ", " min c_u ", "g/kgs", c_u, 1000.0);
    searchMinMax_3D(" max e_d ", " min e_d ", "g/kgs", e_d, 1000.0);
    searchMinMax_3D(" max e_l ", " min e_l ", "g/kgs", e_l, 1000.0);
    searchMinMax_3D(" max e_p ", " min e_p ", "g/kgs", e_p, 1000.0);
    searchMinMax_3D(" max g_p ", " min g_p ", "g/kgs", g_p, 1000.0);
    cout << endl 
    << " --- water vapour, cloud water and entropies in the up-and downdraft: " 
        << endl << endl;
    searchMinMax_3D(" max q_v_u ", " min q_v_u ", "g/kg", q_v_u, 1000.0);
    searchMinMax_3D(" max q_c_u ", " min q_c_u ", "g/kg", q_c_u, 1000.0);
    searchMinMax_3D(" max q_v_d ", " min q_v_d ", "g/kg", q_v_d, 1000.0);
    searchMinMax_3D(" max s ", " min s ", "./.", s, 1.0);
    searchMinMax_3D(" max s_u ", " min s_u ", "./.", s_u, 1.0);
    searchMinMax_3D(" max s_d ", " min s_d ", "./.", s_d, 1.0);

    cout << endl << " greenhouse gas: " << endl << endl;
    searchMinMax_3D(" max co2 ", " min co2 ", "ppm", co2, co2_0);
    searchMinMax_3D(" max epsilon ",  " min epsilon ", "%", epsilon, 1.0);

    cout << endl << " forces per unit volume: " << endl << endl;
    searchMinMax_3D(" max pressure force ", " min pressure force ", 
        "N", PressureGradientForce, 1.0);
    searchMinMax_3D(" max buoyancy force ", " min buoyancy force ", 
        "N", BuoyancyForce, 1.0);
    searchMinMax_3D(" max Coriolis force ", " min Coriolis force ", 
        "N", CoriolisForce, 1.0);

    cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" 
        << endl <<
        " results based on two dimensional considerations of the problem" 
        << endl;
    cout << endl << " co2 distribution: " << endl << endl;
    searchMinMax_2D(" max co2_total ", " min co2_total ", 
         " ppm ", co2_total, co2_0);
    cout << endl << " precipitation: " << endl << endl;
    searchMinMax_2D(" max precipitation ", " min precipitation ", 
        "mm/d", Precipitation, 8.64e4);
    searchMinMax_2D(" max precipitable NASA ", " min precipitable NASA ", 
        "mm/d", precipitation_NASA, 1.0);
    searchMinMax_2D(" max precipitable water ", " min precipitable water ", 
        "mm", precipitable_water, 1.0);

    cout << endl << " radiation, latent and sensible energies: " 
        << endl << endl;
    searchMinMax_3D(" max 3D Q radiation ", " min 3D Q radiation ", 
        "W/m2", radiation, 1.0);
    searchMinMax_3D(" max 3D Q latent ", " min 3D Q latent ", 
        "W/m2", Q_Latent, 1.0);
    searchMinMax_3D(" max 3D Q sensible ", " min 3D Q sensible ", 
        "W/m2", Q_Sensible, 1.0);

    cout << endl << " energies at see level without convection influence: " 
        << endl << endl;
    searchMinMax_2D(" max 2D Q bottom ", " min 2D Q bottom heat ", 
        "W/m2", Q_bottom, 1.0);
    searchMinMax_2D(" max 2D Q latent ", " min 2D Q latent heat ", 
        "W/m2", Q_latent, 1.0);
    searchMinMax_2D(" max 2D Q sensible ", " min 2D Q sensible heat ", 
        "W/m2", Q_sensible, 1.0);

    cout << endl << " secondary data: " << endl << endl;
    searchMinMax_2D(" max Vegetation ", " min Vegetion ", 
        " ./.", Vegetation, 1.0);
    searchMinMax_2D(" max heat Evaporation ", " min heat Evaporation ", 
        " kJ/kg", Q_Evaporation, 1.0);
    searchMinMax_2D(" max Evaporation Dalton ", " min Evaporation Dalton ", 
        "mm/d", Evaporation_Dalton, 1.0);
    searchMinMax_2D(" max Evaporation Penman ", " min Evaporation Penman ", 
        "mm/d", Evaporation_Penman, 1.0);

    cout << endl << " properties of the atmosphere at the surface: " 
        << endl << endl;
    searchMinMax_2D(" max 2D albedo ", " min 2D albedo ", 
        "%", albedo, 1.0);
    searchMinMax_2D(" max 2D epsilon_2D ", " min 2D epsilon_2D ", 
        "%", epsilon_2D, 1.0);

    cout << endl << " topography and the vertical cloud extension: " 
        << endl << endl;
    searchMinMax_2D(" max 2D topography ", " min 2D topography ", 
        "m", Topography, 1.0);
    searchMinMax_2D(" max cloud base ",  " min cloud base ", 
        "./.", i_Base, 1);
    searchMinMax_2D(" max cloud top ",  " min cloud top ", 
        "./.", i_LFS, 1);
    return;
}
/*
*
*/
void cAtmosphereModel::run_data_atm(){
// determination of temperature and pressure by the law of Clausius-Clapeyron for water vapour concentration
// reaching saturation of water vapour pressure leads to formation of rain or ice
// precipitation and cloud formation by formulas from Häckel
// dry adiabatic lapse rate and saturated adiabatic lapse rate = temperature decrease with hight
// SL stands for sea level
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
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            precipitable_water.y[j][k] = 0.;
            Evaporation_Penman.y[j][k] = 0.;
            Evaporation_Dalton.y[j][k] = 0.;
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            // only on the sea surface
            if((is_air(h, 0, j, k))){
                p_stat.x[0][j][k] = (r_air * R_Air 
                    * t.x[0][j][k] * t_0) * 0.01;  // given in hPa
                double t_Celsius = t.x[0][j][k] * t_0 - t_0;
                double e = c.x[0][j][k] * p_stat.x[0][j][k]/ep;  // water vapour pressure in Pa
                double t_denom = t_Celsius + 234.175;
                double E = hp * exp(17.0809 * t_Celsius/t_denom);  // saturation vapour pressure in the water phase for t > 0°C in hPa
                double sat_deficit = E - e;  // saturation deficit in hPa/K
                double del_gam = 0.439 + 0.0112 * t.x[0][j][k];
                double gam_del = 0.5495 + 0.01119 * t.x[0][j][k];
                double e_sa = E * 100.0;
                double r = e/e_sa;
                double u_bar = sqrt((v.x[0][j][k] 
                    * v.x[0][j][k] + w.x[0][j][k] 
                    * w.x[0][j][k])/2.0) * u_0;
                double R_net = short_wave_radiation[j]/radiation.x[0][j][k];
                if(t_Celsius >= - 2.0) Q_Evaporation.y[j][k] =
                    (2500.8 - 2.372 * (t.x[0][j][k] * t_0 - t_0));    // heat of Evaporation of water in [kJ/kg] (Kuttler) => variable lv
                else  Q_Evaporation.y[j][k] = (2500.8 - 2.372 *
                    (t.x[0][j][k] * t_0 - t_0)) + 300.0; // heat of Evaporation of ice + 300 [kJ/kg]
                Q_radiation.y[j][k] = radiation.x[0][j][k];  // long wave radiation in [W/m2]
                Q_latent.y[j][k] = Q_Latent.x[0][j][k];  // latente heat in [W/m2]
                Q_sensible.y[j][k] = Q_Sensible.x[0][j][k];  // sensible heat in [W/m2]
                Q_bottom.y[j][k] = - (radiation.x[0][j][k] 
                    - Q_Latent.x[0][j][k] - Q_Sensible.x[0][j][k]);  // difference understood as heat of the ground
                if(is_air(h, 0, j, k))
                    Evaporation_Dalton.y[j][k] = 
                        C_Dalton(1, j, k, coeff_Dalton, u_0, v, w) // for the ocean surface
                        * sat_deficit * 24.0;  // mm/h in mm/d
                else  Evaporation_Dalton.y[j][k] = 0.0; 
                if(Evaporation_Dalton.y[j][k] <= 0.0) 
                    Evaporation_Dalton.y[j][k] = 0.0;
                if(is_land(h, 0, j, k))
                    Evaporation_Penman.y[j][k] = del_gam * R_net/lv  // in mm/d
                        + gam_del * 0.0026 * (1.0 + 0.54 * u_bar) 
                        * (1.0 - r) * e_sa;                          // for the land surface with vegetation
                else  Evaporation_Penman.y[j][k] = 0.0;
                if(Evaporation_Penman.y[j][k] <= 0.0) Evaporation_Penman.y[j][k] = 0.0;
            }
            for(int i = 0; i <= i_topography[j][k]+1; i++){
                if((is_land(h, i, j, k))
                    &&(is_air(h, i+1, j, k))){
                    double t_Celsius = t.x[i][j][k] * t_0 - t_0;
                    double e = c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // water vapour pressure in hPa
                    double t_denom = t_Celsius + 234.175;
                    double E = hp * exp(17.0809 * t_Celsius/t_denom);
                    double sat_deficit = E - e;  // saturation deficit in hPa
                    double del_gam = 0.439 + 0.0112 * t.x[i][j][k];
                    double gam_del = 0.5495 + 0.01119 * t.x[i][j][k];
                    double e_sa = E * 100.0;
                    double r = e/e_sa;
                    double u_bar = sqrt((v.x[i][j][k] 
                        * v.x[i][j][k] + w.x[i][j][k] 
                        * w.x[i][j][k])/2.0) * u_0;
                    double R_net = short_wave_radiation[j]/radiation.x[0][j][k];
                    if(is_land(h, i, j, k))  Q_Evaporation.y[j][k] = 2300.0;  // minimum value used for printout
                    if(is_water(h, 0, j, k))
                        Evaporation_Dalton.y[j][k] = 
                            C_Dalton(1, j, k, coeff_Dalton, u_0, v, w) 
                            * sat_deficit * 24.0;  // mm/h in mm/d
                    else  Evaporation_Dalton.y[j][k] = 0.0; 
                    if(Evaporation_Dalton.y[j][k] <= 0.0) 
                        Evaporation_Dalton.y[j][k] = 0.0;
                    if(is_land(h, 0, j, k))
                        Evaporation_Penman.y[j][k] = del_gam * R_net/lv  // in mm/d
                            + gam_del * 0.0026 * (1.0 + 0.54 * u_bar) 
                            * (1.0 - r) * e_sa;
                    else  Evaporation_Penman.y[j][k] = 0.0;
                    if(Evaporation_Penman.y[j][k] <= 0.0) Evaporation_Penman.y[j][k] = 0.0;
                }
            }
        }
    }
    double coriolis = 1.0;
    double coeff_Coriolis = r_air * u_0; // coefficient for Coriolis term = 9.6328
    double coeff_buoy = r_air * (u_0 * u_0)/L_atm; // coefficient for bouancy term = 0.2871
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            BuoyancyForce.x[im-1][j][k] = 
                BuoyancyForce.x[im-4][j][k] 
                - 3.0 * BuoyancyForce.x[im-3][j][k] 
                + 3.0 * BuoyancyForce.x[im-2][j][k];  // extrapolation
            CoriolisForce.x[im-1][j][k] = 
                CoriolisForce.x[im-4][j][k] 
                - 3.0 * CoriolisForce.x[im-3][j][k] 
                + 3.0 * CoriolisForce.x[im-2][j][k];  // extrapolation
            PressureGradientForce.x[im-1][j][k] = 
                PressureGradientForce.x[im-4][j][k] 
                - 3.0 * PressureGradientForce.x[im-3][j][k] 
                + 3.0 * PressureGradientForce.x[im-2][j][k];  // extrapolation

            double dpdr = p_dyn.x[0][j][k] = p_dyn.x[3][j][k] 
                - 3. * p_dyn.x[2][j][k] + 3. * p_dyn.x[1][j][k];  // extrapolation
            double sinthe = sin(the.z[j]);
            double costhe = cos(the.z[j]);
            double coriolis_rad = coriolis * 2.0 * omega
                * costhe * w.x[0][j][k];
            double coriolis_the = - coriolis * 2.0 * omega
                * sinthe * w.x[0][j][k];
            double coriolis_phi = coriolis * 2.0 * omega
                * (sinthe * v.x[0][j][k] 
                - costhe * u.x[0][j][k]);
            CoriolisForce.x[0][j][k] = coeff_Coriolis 
                * sqrt((pow (coriolis_rad,2) 
                + pow (coriolis_the,2) 
                + pow (coriolis_phi,2))/3.0);
            BuoyancyForce.x[0][j][k] = buoyancy 
                * r_air * (t.x[0][j][k] - 1.0) * g;
            PressureGradientForce.x[0][j][k] = - coeff_buoy * dpdr;
            if(is_land(h, 0, j, k)){
                BuoyancyForce.x[0][j][k] = 0.0;
                PressureGradientForce.x[0][j][k] = 0.0;
                CoriolisForce.x[0][j][k] = 0.0;
            }
        }
    }
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            BuoyancyForce.x[i][0][k] = c43 * BuoyancyForce.x[i][1][k] 
                - c13 * BuoyancyForce.x[i][2][k];
            BuoyancyForce.x[i][jm-1][k] = c43 * BuoyancyForce.x[i][jm-2][k]       
                - c13 * BuoyancyForce.x[i][jm-3][k];
            CoriolisForce.x[i][0][k] = c43 * CoriolisForce.x[i][1][k] 
                - c13 * CoriolisForce.x[i][2][k];
            CoriolisForce.x[i][jm-1][k] = c43 * CoriolisForce.x[i][jm-2][k]       
                - c13 * CoriolisForce.x[i][jm-3][k];
            PressureGradientForce.x[i][0][k] = c43 * PressureGradientForce.x[i][1][k] 
                - c13 * PressureGradientForce.x[i][2][k];
            PressureGradientForce.x[i][jm-1][k] = c43 * PressureGradientForce.x[i][jm-2][k]       
                - c13 * PressureGradientForce.x[i][jm-3][k];
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 1; j < jm-1; j++){
            BuoyancyForce.x[i][j][0] = c43 * BuoyancyForce.x[i][j][1] 
                - c13 * BuoyancyForce.x[i][j][2];
            BuoyancyForce.x[i][j][km-1] = c43 * BuoyancyForce.x[i][j][km-2] 
                - c13 * BuoyancyForce.x[i][j][km-3];
            BuoyancyForce.x[i][j][0] = BuoyancyForce.x[i][j][km-1] =
               (BuoyancyForce.x[i][j][0] + BuoyancyForce.x[i][j][km-1])/ 2.;
            CoriolisForce.x[i][j][0] = c43 * CoriolisForce.x[i][j][1] 
                - c13 * CoriolisForce.x[i][j][2];
            CoriolisForce.x[i][j][km-1] = c43 * CoriolisForce.x[i][j][km-2] 
                - c13 * CoriolisForce.x[i][j][km-3];
            CoriolisForce.x[i][j][0] = CoriolisForce.x[i][j][km-1] =
               (CoriolisForce.x[i][j][0] + CoriolisForce.x[i][j][km-1])/ 2.;
            PressureGradientForce.x[i][j][0] = c43 * PressureGradientForce.x[i][j][1] 
                - c13 * PressureGradientForce.x[i][j][2];
            PressureGradientForce.x[i][j][km-1] = c43 * PressureGradientForce.x[i][j][km-2] 
                - c13 * PressureGradientForce.x[i][j][km-3];
            PressureGradientForce.x[i][j][0] = PressureGradientForce.x[i][j][km-1] =
               (PressureGradientForce.x[i][j][0] + PressureGradientForce.x[i][j][km-1])/ 2.;
        }
    }

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] < 0.) c.x[i][j][k] = 0.;
                if(cloud.x[i][j][k] < 0.) cloud.x[i][j][k] = 0.;
                if(ice.x[i][j][k] < 0.) ice.x[i][j][k] = 0.;
            }
        }
    }
    precipitable_water.y[0][0] = 0.;
    precipitable_water.y[0][0] = 0.;
    co2_total.y[0][0] = 0.;
//    int j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;
    string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon;
    double precipitation_NASA_average = 0.;
    double precipitablewater_average = 0.;
    double precipitation_average = 0.;
    double Evaporation_Penman_average = 0.;
    double Evaporation_Dalton_average = 0.;
    double Evaporation_average = 0.;
    double co2_average = 0.;
    double temperature_NASA_average = 0.;
    double temperature_average = 0.;
    double temperature_expected_average = 0.;
/*
    int j_loc_Dresden = 39;
    int k_loc_Dresden = 346;
    int j_loc_Sydney = 123;
    int k_loc_Sydney = 151;
    int j_loc_Pacific = 90;
    int k_loc_Pacific = 180;
*/
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            co2_total.y[j][k] = co2.x[0][j][k];
            for(int i = 0; i < im-1; i++){
                double e = 100.0 * c.x[i][j][k] * p_stat.x[i][j][k]/ep;  // water vapour pressure in Pa
                double a = e/(R_WaterVapour * t.x[i][j][k] * t_0);  // absolute humidity in kg/m³
                double step = get_layer_height(i+1) - get_layer_height(i);
                precipitable_water.y[j][k] +=  a * step;
                // precipitable_water mass in 1 kg/m² compares to 1 mm height, with water density kg/(m² * mm)
/*
                if((j==90)&&(k==180))  cout << endl
                    << "  precipitable water" << endl
                    << "  i = " << i << endl
                    << "  a = " << a 
                    << "  e = " << e << endl 
                    << "   " << height << "   " << precipitable_water.y[j][k] <<  endl;
*/
            }
        }
    }
//    double coeff_prec = 8.64e4;  // dimensions see below
    double coeff_prec = 1.0;  // dimensions see below    no convertion from mm/s to mm/a
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Precipitation.y[j][k] = coeff_prec * (P_rain.x[0][j][k] 
                + P_snow.x[0][j][k] + P_conv.x[0][j][k]);
            // 60 s * 60 m * 24 h = 86400 s == 1 d
            // Precipitation, P_conv, P_rain and P_snow in kg/(m²*s) = mm/s
            // Precipitation in 86400 * kg/(m²*d) = 86400 mm/d
            // kg/(m² * s) == mm/s(Kraus, p. 94)
        }
    }


    // vegetation as function of precipitation, timberline and temperature
    double max_Precipitation = 0.0;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(Precipitation.y[j][k] > max_Precipitation){
                max_Precipitation = Precipitation.y[j][k];
            }
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            if(max_Precipitation > 0 && is_land(h, 0, j, k) 
                && !(get_layer_height(i_mount) > 4400.0) // above the vegetation line of 4.4 km no vegetation
                && !((t.x[0][j][k] * t_0 - t_0) < - 40.0)){ //  vegetation >= -40°C
                Vegetation.y[j][k] = Precipitation.y[j][k] 
                    /max_Precipitation; // actual vegetation areas
            }else{
                Vegetation.y[j][k] = 0.;
            }
        }
    }
    temperature_NASA_average = AtomUtils::GetMean_2D(jm, km, temperature_NASA);
    temperature_average = (AtomUtils::GetMean_3D(jm, km, t) - 1.0) * t_0;
//    temperature_expected_average = get_global_temperature_from_curve(*get_current_time());
    temperature_expected_average = 
        get_temperatures_from_curve(*get_current_time(), 
        m_global_temperature_curve);
    precipitablewater_average = AtomUtils::GetMean_2D(jm, km, precipitable_water);
    precipitation_average = 365. * AtomUtils::GetMean_2D(jm, km, Precipitation);
    precipitation_NASA_average = 365. * AtomUtils::GetMean_2D(jm, km, precipitation_NASA);
    co2_average = AtomUtils::GetMean_2D(jm, km, co2_total);
    Evaporation_Penman_average = 365. * AtomUtils::GetMean_2D(jm, km, Evaporation_Penman);
    Evaporation_Dalton_average = 365. * AtomUtils::GetMean_2D(jm, km, Evaporation_Dalton);
    Evaporation_average = Evaporation_Penman_average + Evaporation_Dalton_average;
    cout.precision(2);
    level = "m";
    deg_north = "°N";
    deg_south = "°S";
    deg_west = "°W";
    deg_east = "°E";
    string name_Value_1 = " radiation emission";
    string name_Value_2 = " latent heat ";
    string name_Value_3 = " sensible heat ";
    string name_Value_4 = " bottom heat ";
    string name_Value_5 = " Evaporation Penman ";
    string name_Value_6 = " Evaporation Dalton ";
    string name_Value_7 = " precipitable water average ";
    string name_Value_8 = " precipitation average per year ";
    string name_Value_9 = " precipitation average per day ";
    string name_Value_10 = " precipitation NASA average per year ";
    string name_Value_11 = " precipitation NASA average per day ";
    string name_Value_12 = " Evaporation_Penman_average per year ";
    string name_Value_13 = " Evaporation_Penman_average per day ";
    string name_Value_14 = " Evaporation_Dalton_average per year ";
    string name_Value_15 = " Evaporation_Dalton_average per day ";
    string name_Value_16 = " Evaporation_average per year ";
    string name_Value_17 = " Evaporation_average per day ";
    string name_Value_18 = " to fill ";
    string name_Value_19 = " to fill ";
    string name_Value_20 = " to fill ";
    string name_Value_21 = " to fill ";
    string name_Value_22 = " co2_average ";
    string name_Value_23 = " precipitable water ";
    string name_Value_24 = " precipitation ";
    string name_Value_25 = " temperature_NASA_average ";
    string name_Value_26 = " temperature_average ";
    string name_Value_27 = " temperature_expected_average ";
    string name_unit_wm2 = " W/m2";
    string name_unit_mmd = " mm/d";
    string name_unit_mm = " mm";
    string name_unit_mma = " mm/a";
    string name_unit_ppm = " ppm";
    string name_unit_t = " deg";
    string heading = " printout of surface data at predefinded locations: level, latitude, longitude";
    string heading_Dresden = " City of Dresden, Germany, Europe";
    string heading_Sydney = " City of Sydney, New South Wales, Australia";
    string heading_Pacific = " Equator in the central Pacific";
    cout << endl << endl;
/*
    cout << endl << endl << heading << endl << endl;
    int choice = {1};
    preparation:
    switch(choice){
        case 1 :    cout << heading_Dresden << endl;
                        i_loc_level = 0;  // sea level
                        j_loc = j_loc_Dresden;  // 51°N, Dresden Germany
                        k_loc = k_loc_Dresden;  // 14°W, Dresden Germany
                        break;
        case 2 :    cout << heading_Sydney << endl;
                        i_loc_level = 0;  // sea level
                        j_loc = j_loc_Sydney;  // 33°S, Dresden Germany
                        k_loc = k_loc_Sydney;  // 151°E, Dresden Germany
                        break;
        case 3 :    cout << heading_Pacific << endl;
                        i_loc_level = 0;  // sea level
                        j_loc = j_loc_Pacific;  // 0°N, Equator
                        k_loc = k_loc_Pacific;  // 180°E, central Pacific
                        break;
    default :     cout << choice << "error in iterationPrintout member function in class Accuracy" << endl;
    }
    if(j_loc <= 90){
        j_loc_deg = 90 - j_loc;
        deg_lat = deg_north;
    }
    if(j_loc > 90){
        j_loc_deg = j_loc - 90;
        deg_lat = deg_south;
    }
    if(k_loc <= 180){
        k_loc_deg = k_loc;
        deg_lon = deg_east;
    }
    if(k_loc > 180){
        k_loc_deg = 360 - k_loc;
        deg_lon = deg_east;
    }
    double Value_1 = Q_radiation.y[j_loc][k_loc];
    double Value_2 = Q_latent.y[j_loc][k_loc];
    double Value_3 = Q_sensible.y[j_loc][k_loc];
    double Value_4 = Q_bottom.y[j_loc][k_loc];
    cout << setw(6)<< i_loc_level << setw(2)<< level << setw(5)<< j_loc_deg
        << setw(3)<< deg_lat << setw(4)<< k_loc_deg << setw(3)<< deg_lon
        << "  " << setiosflags(ios::left)<< setw(25)<< setfill('.')<< name_Value_1
        << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')
        << Value_1 << setw(6)<< name_unit_wm2 << "   " << setiosflags(ios::left)
        << setw(25)<< setfill('.')<< name_Value_2 << " = " << resetiosflags(ios::left)
        << setw(7)<< fixed << setfill(' ')<< Value_2 << setw(6)<< name_unit_wm2
        << "   " << setiosflags(ios::left)<< setw(25)<< setfill('.')<< name_Value_3
        << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')
        << Value_3 << setw(6)<< name_unit_wm2 << "   " << setiosflags(ios::left)
        << setw(25)<< setfill('.')<< name_Value_4 << " = " << resetiosflags(ios::left)
        << setw(7)<< fixed << setfill(' ')<< Value_4 << setw(6)
        << name_unit_wm2 << endl;
    double Value_17 = 0.;
    double Value_18 = 0.;
    double Value_19 = 0.;
    double Value_23 = precipitable_water.y[j_loc][k_loc];
    cout << setw(6)<< i_loc_level << setw(2)<< level << setw(5)<< j_loc_deg
        << setw(3)<< deg_lat << setw(4)<< k_loc_deg << setw(3)<< deg_lon
        << "  " << setiosflags(ios::left)<< setw(25)<< setfill('.')<< name_Value_23
        << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')
        << Value_23 << setw(6)<< name_unit_mm << "   " << setiosflags(ios::left)
        << setw(25)<< setfill('.')<< name_Value_19 << " = " << resetiosflags(ios::left)
        << setw(7)<< fixed << setfill(' ')<< Value_17 << setw(6)<< name_unit_wm2
        << "   " << setiosflags(ios::left)<< setw(25)<< setfill('.')<< name_Value_20
        << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')
        << Value_18 << setw(6)<< name_unit_wm2 << "   " << setiosflags(ios::left)
        << setw(25)<< setfill('.')<< name_Value_21 << " = " << resetiosflags(ios::left)
        << setw(7)<< fixed << setfill(' ')<< Value_19 << setw(6)
        << name_unit_wm2 << endl;
    double Value_14 = 0.;
    double Value_15 = 0.;
    double Value_16 = 0.;
    double Value_24 = Precipitation.y[j_loc][k_loc];
    cout << setw(6)<< i_loc_level << setw(2)<< level << setw(5)<< j_loc_deg
        << setw(3)<< deg_lat << setw(4)<< k_loc_deg << setw(3)<< deg_lon
        << "  " << setiosflags(ios::left)<< setw(25)<< setfill('.')<< name_Value_24
        << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')
        << Value_24 << setw(6)<< name_unit_mmd << "   " << setiosflags(ios::left)
        << setw(25)<< setfill('.')<< name_Value_16 << " = " << resetiosflags(ios::left)
        << setw(7)<< fixed << setfill(' ')<< Value_14 << setw(6)<< name_unit_wm2
        << "   " << setiosflags(ios::left)<< setw(25)<< setfill('.')<< name_Value_17
        << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')
        << Value_15 << setw(6)<< name_unit_wm2 << "   " << setiosflags(ios::left)
        << setw(25)<< setfill('.')<< name_Value_18 << " = " << resetiosflags(ios::left)
        << setw(7)<< fixed << setfill(' ')<< Value_16 << setw(6)
        << name_unit_wm2 << endl;
    double Value_5 = Evaporation_Penman.y[j_loc][k_loc];
    double Value_6 = Evaporation_Dalton.y[j_loc][k_loc];
    cout << setw(6)<< i_loc_level << setw(2)<< level << setw(5)<< j_loc_deg
        << setw(3)<< deg_lat << setw(4)<< k_loc_deg << setw(3)<< deg_lon
        << "  " << setiosflags(ios::left)<< setw(25)<< setfill('.')<< name_Value_5
        << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')
        << Value_5 << setw(6)<< name_unit_mmd << "   " << setiosflags(ios::left)
        << setw(25)<< setfill('.')<< name_Value_6 << " = " << resetiosflags(ios::left)
        << setw(7)<< fixed << setfill(' ')<< Value_6 << setw(6)
        << name_unit_mmd << endl << endl;
    choice++;
    if(choice <= 3) goto preparation;
    cout << endl;
*/
    double Value_7 = precipitablewater_average;
    double Value_8 = precipitation_average * 8.64e4;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_7 << " = " << resetiosflags(ios::left) << setw(7) << fixed
        << setfill(' ') << Value_7 << setw(6) << name_unit_mm << "   " << setiosflags(ios::left)
        << setw(40) << setfill('.') << name_Value_8 << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << Value_8 << setw(6) << name_unit_mma
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.') << name_Value_9
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_8/365. << setw(6)<< name_unit_mmd << endl;
    double Value_10 = precipitation_NASA_average;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_7 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_7 << setw(6) << name_unit_mm
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.') << name_Value_10
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_10 << setw(6) << name_unit_mma << "   " << setiosflags(ios::left)
        << setw(40) << setfill('.') << name_Value_11 << " = "
        << resetiosflags(ios::left)<< setw(7) << fixed << setfill(' ')
        << Value_10/365. << setw(6) << name_unit_mmd << endl;
    double Value_13 = Evaporation_Dalton_average;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_7 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_7 << setw(6) << name_unit_mm
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_14 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_13 << setw(6) << name_unit_mma
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_15 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_13/365. << setw(6)
        << name_unit_mmd << endl;
    double Value_9 = co2_average * co2_0;
    double Value_12 = Evaporation_Penman_average;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_22 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_9 << setw(6) << name_unit_ppm
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.') << name_Value_12
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_12 << setw(6) << name_unit_mma << "   " << setiosflags(ios::left)
        << setw(40) << setfill('.') << name_Value_13 << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_12/365. << setw(6) << name_unit_mmd << endl;
    double Value_20 = Evaporation_average;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_22 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_9 << setw(6) << name_unit_ppm
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.') << name_Value_16
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_20 << setw(6) << name_unit_mma << "   " << setiosflags(ios::left)
        << setw(40) << setfill('.') << name_Value_17 << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_20/365. << setw(6) << name_unit_mmd << endl << endl;
    double Value_25 =temperature_NASA_average;
    double Value_26 = temperature_average;
    double Value_27 = temperature_expected_average;
    cout << setw(6)<< setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_25 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_25 << setw(6) << name_unit_t << "   "
        << setiosflags(ios::left) << setw(40) << setfill('.') << name_Value_26
        << " = " << resetiosflags(ios::left)<< setw(7) << fixed << setfill(' ')
        << Value_26 << setw(6) << name_unit_t << "   " << setiosflags(ios::left)
        << setw(40) << setfill('.') << name_Value_27 << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_27 << setw(6) << name_unit_t << endl << endl << endl;
    return;
}
/*
*
*/
