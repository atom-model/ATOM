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
    double coeff_p = 1e-2 * r_air * u_0 * u_0; // in hPa = 0.77056
    searchMinMax_3D(" max 3D temperature ", " min 3D temperature ", 
        " deg", t, 273.15, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max 3D u-component ", " min 3D u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max 3D v-component ", " min 3D v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max 3D w-component ", " min 3D w-component ", "m/s", w, u_0);
    searchMinMax_3D(" max 3D pressure dynamic ", " min 3D pressure dynamic ", "hPa", p_dyn, coeff_p);
    searchMinMax_3D(" max 3D pressure static ", " min 3D pressure static ", "hPa", p_stat, 1.);
    cout << endl << " energies in the three dimensional space: " << endl << endl;
    searchMinMax_3D(" max 3D radiation ",  " min 3D radiation ",  "W/m2", radiation_3D, 1.);
    searchMinMax_3D(" max 3D sensible heat ", " min 3D sensible heat ", "W/m2", Q_Sensible, 1.);
    searchMinMax_3D(" max 3D latent heat ", " min 3D latent heat ", "W/m2", Q_Latent, 1.);
    cout << endl << " greenhouse gases: " << endl << endl;
    searchMinMax_3D(" max 3D water vapour ",  " min 3D water vapour ", "g/kg", c, 1000.);
    searchMinMax_3D(" max 3D cloud water ", " min 3D cloud water ", "g/kg", cloud, 1000.);
    searchMinMax_3D(" max 3D cloud ice ", " min 3D cloud ice ", "g/kg", ice, 1000.);
//    searchMinMax_3D(" max 3D rain ", " min 3D rain ", "mm/d", P_rain, 8.64e4);
//    searchMinMax_3D(" max 3D snow ", " min 3D snow ", "mm/d", P_snow, 8.64e4);
    searchMinMax_3D(" max 3D rain ", " min 3D rain ", "mm/d", P_rain, 1.);
    searchMinMax_3D(" max 3D snow ", " min 3D snow ", "mm/d", P_snow, 1.);
    searchMinMax_3D(" max 3D co2 ", " min 3D co2 ", "ppm", co2, co2_0);
    searchMinMax_3D(" max 3D epsilon ",  " min 3D epsilon ", "%", epsilon_3D, 1.);
    searchMinMax_3D(" max 3D buoyancy force ", " min 3D buoyancy force ", "kN/m3", BuoyancyForce, 1.);
    cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" 
        << endl <<
        " results based on two dimensional considerations of the problem" << endl;
    cout << endl << " co2 distribution row-wise: " << endl << endl;
   searchMinMax_2D(" max co2_total ", " min co2_total ", " ppm ", co2_total, co2_0);
    cout << endl << " precipitation: " << endl << endl;
    searchMinMax_2D(" max precipitation ", " min precipitation ", "mm/d", Precipitation, 1.);
    max_Precipitation = out_maxValue();
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
    searchMinMax_2D(" max 2D epsilon ", " min 2D epsilon ", "%", epsilon, 1.);
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


