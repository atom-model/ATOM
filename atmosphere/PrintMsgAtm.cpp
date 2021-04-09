#include <iomanip>
#include "cAtmosphereModel.h"

void cAtmosphereModel::print_welcome_msg(){
    if(verbose){
        cout << endl << endl << endl;
        cout << "***** Atmosphere General Circulation Model (AGCM) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 4 additional transport equations to describe the water vapour, cloud water, cloud ice and co2 concentration" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** water vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code (German Weather Forecast)" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** co2 concentration appears in the absorptivity of the radiation models" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz (roger.grundmann@web.de)" << endl << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
        has_welcome_msg_printed = true;
    }
}
/*
*
*/
void cAtmosphereModel::print_final_remarks(){
    cout << endl << "***** end of the Atmosphere General Circulation Modell (AGCM) *****" << endl << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
    cout << "\n\n\n\n";
}
/*
*
*/
void cAtmosphereModel::print_min_max_atm(){
    cout << endl << " flow properties: " << endl << endl;
    searchMinMax_3D(" max temperature ", " min temperature ", 
        " deg", t, 273.15, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max u-component ", " min u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max v-component ", " min v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max w-component ", " min w-component ", "m/s", w, u_0);
    searchMinMax_3D(" max pressure dynamic ", " min pressure dynamic ", "hPa", p_dyn, r_air * u_0 * u_0 * 1e-2);
    searchMinMax_3D(" max pressure static ", " min pressure static ", "hPa", p_stat, 1.);
    searchMinMax_3D(" max density dry air ", " min density dry air ", "kg/m3", r_dry, 1.);
    searchMinMax_3D(" max density wet air ", " min density wet air ", "kg/m3", r_humid, 1.);
    cout << endl << " energies: " << endl << endl;
    searchMinMax_3D(" max radiation ",  " min radiation ",  "W/m2", radiation, 1.);
    searchMinMax_3D(" max sensible heat ", " min sensible heat ", "W/m2", Q_Sensible, 1.);
    searchMinMax_3D(" max latent heat ", " min latent heat ", "W/m2", Q_Latent, 1.);
    cout << endl << " cloud thermodynamics: " << endl << endl;
//    searchMinMax_3D(" max cloud base ",  " min cloud base ", "./.", i_Base, 1);
//    searchMinMax_3D(" max cloud top ",  " min cloud top ", "./.", i_LFS, 1);
    searchMinMax_3D(" max water vapour ",  " min water vapour ", "g/kg", c, 1000.);
    searchMinMax_3D(" max cloud water ", " min cloud water ", "g/kg", cloud, 1000.);
    searchMinMax_3D(" max cloud ice ", " min cloud ice ", "g/kg", ice, 1000.);
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
    searchMinMax_3D(" max MC_t ", " min MC_t ", "K/s", MC_t, 1.);
    searchMinMax_3D(" max MC_q ", " min MC_q ", "g/kgs", MC_q, 1000.);
    searchMinMax_3D(" max MC_v ", " min MC_v ", "m/ss", MC_v, u_0);
    searchMinMax_3D(" max MC_w ", " min MC_w ", "m/ss", MC_w, u_0);
    cout << endl << " --- velocity components in the up- and downdraft: " << endl << endl;
    searchMinMax_3D(" max u_u ", " min u_u ", "m/s", u_u, 1.);
    searchMinMax_3D(" max v_u ", " min v_u ", "m/s", v_u, 1.);
    searchMinMax_3D(" max w_u ", " min w_u ", "m/s", w_u, 1.);
    searchMinMax_3D(" max u_d ", " min u_d ", "m/s", u_d, 1.);
    searchMinMax_3D(" max v_d ", " min v_d ", "m/s", v_d, 1.);
    searchMinMax_3D(" max w_d ", " min w_u ", "m/s", w_d, 1.);
    cout << endl << " --- condensation and evaporation terms in the up- and downdraft: " << endl << endl;
    searchMinMax_3D(" max c_u ", " min c_u ", "g/kgs", c_u, 1000.);
    searchMinMax_3D(" max e_d ", " min e_d ", "g/kgs", e_d, 1000.);
    searchMinMax_3D(" max e_l ", " min e_l ", "g/kgs", e_l, 1000.);
    searchMinMax_3D(" max e_p ", " min e_p ", "g/kgs", e_p, 1000.);
    searchMinMax_3D(" max g_p ", " min g_p ", "g/kgs", g_p, 1000.);
    cout << endl << " --- water vapour, cloud water and entropies in the up-and downdraft: " << endl << endl;
    searchMinMax_3D(" max q_v_u ", " min q_v_u ", "g/kg", q_v_u, 1000.);
    searchMinMax_3D(" max q_c_u ", " min q_c_u ", "g/kg", q_c_u, 1000.);
    searchMinMax_3D(" max q_v_d ", " min q_v_d ", "g/kg", q_v_d, 1000.);
    searchMinMax_3D(" max s ", " min s ", "./.", s, 1.);
    searchMinMax_3D(" max s_u ", " min s_u ", "./.", s_u, 1.);
    searchMinMax_3D(" max s_d ", " min s_d ", "./.", s_d, 1.);
    cout << endl << " greenhouse gas: " << endl << endl;
    searchMinMax_3D(" max co2 ", " min co2 ", "ppm", co2, co2_0);
    searchMinMax_3D(" max epsilon ",  " min epsilon ", "%", epsilon, 1.);
    cout << endl << " forces per unit volume: " << endl << endl;
    searchMinMax_3D(" max pressure force ", " min pressure force ", "N/m3", PressureGradientForce, 1.);
    searchMinMax_3D(" max buoyancy force ", " min buoyancy force ", "N/m3", BuoyancyForce, 1.);
    searchMinMax_3D(" max Coriolis force ", " min Coriolis force ", "N/m3", CoriolisForce, 1.);
    cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" 
        << endl <<
        " results based on two dimensional considerations of the problem" << endl;
    cout << endl << " co2 distribution: " << endl << endl;
   searchMinMax_2D(" max co2_total ", " min co2_total ", " ppm ", co2_total, co2_0);
    cout << endl << " precipitation: " << endl << endl;
    searchMinMax_2D(" max precipitation ", " min precipitation ", "mm/d", Precipitation, 8.64e4);
    max_Precipitation = out_maxValue();
    searchMinMax_2D(" max precipitable NASA ", " min precipitable NASA ", "mm/d", precipitation_NASA, 1.);
    searchMinMax_2D(" max precipitable water ", " min precipitable water ", "mm", precipitable_water, 1.);
    cout << endl << " energies at see level without convection influence: " 
        << endl << endl;
    searchMinMax_2D(" max 2D Q radiation ", " min 2D Q radiation ", "W/m2", Q_radiation, 1.);
    searchMinMax_2D(" max 2D Q latent ", " min 2D Q latent ", "W/m2", Q_latent, 1.);
    searchMinMax_2D(" max 2D Q sensible ", " min 2D Q sensible ", "W/m2", Q_sensible, 1.);
    searchMinMax_2D(" max 2D Q bottom ", " min 2D Q bottom heat ", "W/m2", Q_bottom, 1.);
    cout << endl << " secondary data: " << endl << endl;
    searchMinMax_2D(" max heat Evaporation ", " min heat Evaporation ", " kJ/kg", Q_Evaporation, 1.);
    searchMinMax_2D(" max Evaporation Dalton ", " min Evaporation Dalton ", "mm/d", Evaporation_Dalton, 1.);
    cout << endl << " properties of the atmosphere at the surface: " 
        << endl << endl;
    searchMinMax_2D(" max 2D albedo ", " min 2D albedo ", "%", albedo, 1.);
    searchMinMax_2D(" max 2D epsilon_2D ", " min 2D epsilon_2D ", "%", epsilon_2D, 1.);
    searchMinMax_2D(" max 2D topography ", " min 2D topography ", "m", Topography, 1.);
}
/*
*
*/
float cAtmosphereModel::GetMean_3D(int jm, int km, Array &val_3D){
    if(m_node_weights.size() != (unsigned)jm){
        CalculateNodeWeights(jm, km);
    }
    double ret=0., weight=0.;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (val_3D.x[0][j][k]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret+=val_3D.x[0][j][k]*m_node_weights[j][k];
            weight+=m_node_weights[j][k];
        }
    }
    return ret/weight;
}
/*
*
*/
float cAtmosphereModel::GetMean_2D(int jm, int km, Array_2D &val_2D){
    if(m_node_weights.size() != (unsigned)jm){
        CalculateNodeWeights(jm, km);
    }
    double ret=0., weight=0.;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (val_2D.y[j][k]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret+=val_2D.y[j][k]*m_node_weights[j][k];
            weight+=m_node_weights[j][k];
        }
    }
    return ret/weight;
}
/*
*
*/
void cAtmosphereModel::CalculateNodeWeights(int jm, int km){
    //use cosine of latitude as weights for now
    //longitudes: 0-360(km) latitudes: 90-(-90)(jm)
    double weight = 0.;
    m_node_weights.clear();
    for(int i=0; i<jm; i++){
        if(i<=90){
            weight = cos((90-i) * M_PI/180.0);
        }else{
            weight = cos((i-90) * M_PI/180.0);
        }
        m_node_weights.push_back(std::vector<double>());
        m_node_weights[i].resize(km, weight);
    }
    return;
}


