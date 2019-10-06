#include <iomanip>
#include "cAtmosphereModel.h"
#include "MinMax_Atm.h"

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

void cAtmosphereModel::print_final_remarks(){
    cout << endl << "***** end of the Atmosphere General Circulation Modell (AGCM) *****" << endl << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
    cout << "\n\n\n\n";
}

void cAtmosphereModel::print_min_max_values(){
    MinMax_Atm min_max_3d(im, jm, km);
    min_max_3d.searchMinMax_3D(" max 3D temperature ", " min 3D temperature ", 
        "°C", t, h, 273.15, [](double i)->double{return i - 273.15;}, true);
    min_max_3d.searchMinMax_3D(" max 3D u-component ", " min 3D u-component ", 
        "m/s", u, h, u_0);
    min_max_3d.searchMinMax_3D(" max 3D v-component ", " min 3D v-component ", 
        "m/s", v, h, u_0);
    min_max_3d.searchMinMax_3D(" max 3D w-component ", " min 3D w-component ", 
        "m/s", w, h, u_0);
    min_max_3d.searchMinMax_3D(" max 3D pressure dynamic ", " min 3D pressure dynamic ", 
        "hPa", p_dyn, h, 0.768); // 0.768 = 0.01 * r_air *u_0*u_0 in hPa
    min_max_3d.searchMinMax_3D(" max 3D pressure static ", " min 3D pressure static ", 
        "hPa", p_stat, h);
    cout << endl << " energies in the three dimensional space: " << endl << endl;
    min_max_3d.searchMinMax_3D(" max 3D radiation ",  " min 3D radiation ",  "W/m2", 
        radiation_3D, h);
    min_max_3d.searchMinMax_3D(" max 3D sensible heat ", " min 3D sensible heat ", 
        "W/m2", Q_Sensible, h);
    min_max_3d.searchMinMax_3D(" max 3D latent heat ", " min 3D latent heat ", 
        "W/m2", Q_Latent, h);
    cout << endl << " greenhouse gases: " << endl << endl;
    min_max_3d.searchMinMax_3D(" max 3D water vapour ",  " min 3D water vapour ", 
        "g/kg", c, h, 1000.);
    min_max_3d.searchMinMax_3D(" max 3D cloud water ", " min 3D cloud water ", 
        "g/kg", cloud, h, 1000.);
    min_max_3d.searchMinMax_3D(" max 3D cloud ice ", " min 3D cloud ice ", 
        "g/kg", ice, h, 1000.);
    min_max_3d.searchMinMax_3D(" max 3D rain ", " min 3D rain ", "mm/d", 
        P_rain, h, 8.64e4);
    min_max_3d.searchMinMax_3D(" max 3D snow ", " min 3D snow ", "mm/d", 
        P_snow, h, 8.64e4);
    min_max_3d.searchMinMax_3D(" max 3D co2 ", " min 3D co2 ", "ppm", 
        co2, h, 280.);
    min_max_3d.searchMinMax_3D(" max 3D epsilon ",  " min 3D epsilon ", "%", 
        epsilon_3D, h);
    min_max_3d.searchMinMax_3D(" max 3D buoyancy force ", " min 3D buoyancy force ", 
        "kN/m2", BuoyancyForce, h);
    cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" 
        << endl <<
        " results based on two dimensional considerations of the problem" << endl;
    cout << endl << " co2 distribution row-wise: " << endl << endl;
    MinMax_Atm  min_max_2d(jm, km);
    min_max_2d.searchMinMax_2D(" max co2_total ", " min co2_total ", " ppm ", 
        co2_total, h, 280.);
    cout << endl << " precipitation: " << endl << endl;
    min_max_2d.searchMinMax_2D(" max precipitation ", " min precipitation ", 
        "mm/d", Precipitation, h, 1.);
    max_Precipitation = min_max_2d.out_maxValue();
    min_max_2d.searchMinMax_2D(" max precipitable water ", " min precipitable water ", 
        "mm", precipitable_water, h, 1.);
    cout << endl << " energies at see level without convection influence: " 
        << endl << endl;
    min_max_2d.searchMinMax_2D(" max 2D Q radiation ", " min 2D Q radiation ", 
        "W/m2", Q_radiation, h);
    min_max_2d.searchMinMax_2D(" max 2D Q latent ", " min 2D Q latent ", 
        "W/m2", Q_latent, h);
    min_max_2d.searchMinMax_2D(" max 2D Q sensible ", " min 2D Q sensible ", 
        "W/m2", Q_sensible, h);
    min_max_2d.searchMinMax_2D(" max 2D Q bottom ", " min 2D Q bottom heat ", 
        "W/m2", Q_bottom, h);
    cout << endl << " secondary data: " << endl << endl;
    min_max_2d.searchMinMax_2D(" max heat Evaporation ", " min heat Evaporation ", 
        " kJ/kg", Q_Evaporation, h);
    min_max_2d.searchMinMax_2D(" max Evaporation Dalton ", " min Evaporation Dalton ", 
        "mm/d", Evaporation_Dalton, h);
    cout << endl << " properties of the atmosphere at the surface: " 
        << endl << endl;
    min_max_2d.searchMinMax_2D(" max 2D albedo ", " min 2D albedo ", "%", 
        albedo, h);
    min_max_2d.searchMinMax_2D(" max 2D epsilon ", " min 2D epsilon ", "%", 
        epsilon, h);
    min_max_2d.searchMinMax_2D(" max 2D topography ", " min 2D topography ", 
        "m", Topography, h);
    int j_loc_Dresden = 39;
    int k_loc_Dresden = 346;
    int j_loc_Sydney = 123;
    int k_loc_Sydney = 151;
    int j_loc_Pacific = 90;
    int k_loc_Pacific = 180;
    double precipitation_NASA_average = 0.;
    double precipitablewater_average = 0.;
    double precipitation_average = 0.;
    double temperature_surf_average = 0.;
    double Evaporation_Dalton_average = 0.;
    double co2_vegetation_average = 0.;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            precipitation_NASA_average += precipitation_NASA.y[j][k];
            precipitablewater_average += precipitable_water.y[j][k];
            precipitation_average += Precipitation.y[j][k];
            temperature_surf_average += t.x[0][j][k] * t_0 - t_0;
            Evaporation_Dalton_average += Evaporation_Dalton.y[j][k];
            co2_vegetation_average += co2_total.y[j][k];
        }
    }
    temperature_surf_average = (GetMean_3D(jm, km, t) - 1.) * t_0;
    precipitablewater_average = GetMean_2D (jm, km, precipitable_water);
    precipitation_average = 365. * GetMean_2D(jm, km, Precipitation);
    precipitation_NASA_average = 365. * GetMean_2D(jm, km, precipitation_NASA);
    co2_vegetation_average = GetMean_2D(jm, km, co2_total);
    Evaporation_Dalton_average = 365. * GetMean_2D (jm, km, Evaporation_Dalton);
    cout.precision(2);
    const char* level = "m";
    const char* deg_north = "°N";
    const char* deg_south = "°S";
//    const char* deg_west = "°W";
    const char* deg_east = "°E";
    const char* name_Value_1 = " radiation emission";
    const char* name_Value_2 = " latent heat ";
    const char* name_Value_3 = " sensible heat ";
    const char* name_Value_4 = " bottom heat ";
    const char* name_Value_6 = " Evaporation Dalton ";
    const char* name_Value_7 = " precipitable water average ";
    const char* name_Value_8 = " precipitation average per year ";
    const char* name_Value_9 = " precipitation average per day ";
    const char* name_Value_10 = " precipitation NASA average per year ";
    const char* name_Value_11 = " precipitation NASA average per day ";
    const char* name_Value_14 = " Evaporation_Dalton_average per year ";
    const char* name_Value_15 = " Evaporation_Dalton_average per day ";
    const char* name_Value_16 = " to fill ";
    const char* name_Value_17 = " to fill ";
    const char* name_Value_18 = " to fill ";
    const char* name_Value_19 = " to fill ";
    const char* name_Value_20 = " to fill ";
    const char* name_Value_21 = " to fill ";
    const char* name_Value_22 = " co2_average ";
    const char* name_Value_23 = " precipitable water ";
    const char* name_Value_24 = " precipitation ";
    const char* name_Value_25 = " temperature_modern_average ";
    const char* name_Value_26 = " temperature_surf_average ";
    const char* name_Value_27 = " temperature_expected_average ";
    const char* name_unit_wm2 = " W/m2";
    const char* name_unit_mmd = " mm/d";
    const char* name_unit_mm = " mm";
    const char* name_unit_mma = " mm/a";
    const char* name_unit_ppm = " ppm";
    const char* name_unit_t = " C";
    const char* heading = " printout of surface data at predefinded locations: level, latitude, longitude";
    const char* heading_Dresden = " City of Dresden, Germany, Europe";
    const char* heading_Sydney = " City of Sydney, New South Wales, Australia";
    const char* heading_Pacific = " Equator in the central Pacific";
    int i_loc_level = 0;
    int j_loc = 0;
    int k_loc = 0;
    int j_loc_deg = 0;
    int k_loc_deg = 0;
    const char* deg_lat = 0;
    const char* deg_lon = 0;
/*
    int deg_north = 0;
    int deg_south = 0;
    int deg_west = 0;
    int deg_east = 0;

*/
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
    default :    cout << choice << "error in iterationPrintout member function in class Accuracy" << endl;
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
    cout << setw(6) << i_loc_level << setw(2) << level << setw(5) << j_loc_deg
        << setw(3) << deg_lat << setw(4) << k_loc_deg << setw(3) << deg_lon
        << "  " << setiosflags(ios::left) << setw(25) << setfill('.') << name_Value_1
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_1 << setw(6) << name_unit_wm2 << "   " << setiosflags(ios::left)
        << setw(25) << setfill('.') << name_Value_2 << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << Value_2 << setw(6) << name_unit_wm2
        << "   " << setiosflags(ios::left) << setw(25) << setfill('.') << name_Value_3
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_3 << setw(6) << name_unit_wm2 << "   " << setiosflags(ios::left)
        << setw(25) << setfill('.') << name_Value_4 << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << Value_4 << setw(6)
        << name_unit_wm2 << endl;
    double Value_17 = 0.;
    double Value_18 = 0.;
    double Value_19 = 0.;
    double Value_23 = precipitable_water.y[j_loc][k_loc];
    cout << setw(6) << i_loc_level << setw(2) << level << setw(5) << j_loc_deg
        << setw(3) << deg_lat << setw(4) << k_loc_deg << setw(3) << deg_lon
        << "  " << setiosflags(ios::left) << setw(25) << setfill('.') << name_Value_23
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_23 << setw(6) << name_unit_mm << "   " << setiosflags(ios::left)
        << setw(25) << setfill('.') << name_Value_19 << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << Value_17 << setw(6) << name_unit_wm2
        << "   " << setiosflags(ios::left) << setw(25) << setfill('.') << name_Value_20
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_18 << setw(6) << name_unit_wm2 << "   " << setiosflags (ios::left)
        << setw (25) << setfill ('.') << name_Value_21 << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << Value_19 << setw(6)
        << name_unit_wm2 << endl;
    double Value_14 = 0.;
    double Value_15 = 0.;
    double Value_16 = 0.;
    double Value_24 = Precipitation.y[j_loc][k_loc];
    cout << setw(6) << i_loc_level << setw(2) << level << setw(5) << j_loc_deg
        << setw(3) << deg_lat << setw(4) << k_loc_deg << setw(3) << deg_lon
        << "  " << setiosflags(ios::left) << setw(25) << setfill('.') << name_Value_24
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_24 << setw(6) << name_unit_mmd << "   " << setiosflags(ios::left)
        << setw(25) << setfill('.') << name_Value_16 << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << Value_14 << setw(6) << name_unit_wm2
        << "   " << setiosflags(ios::left) << setw(25) << setfill('.') << name_Value_17
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_15 << setw(6) << name_unit_wm2 << "   " << setiosflags(ios::left)
        << setw(25) << setfill('.') << name_Value_18 << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << Value_16 << setw(6)
        << name_unit_wm2 << endl;
    double Value_6 = Evaporation_Dalton.y[j_loc][k_loc];
    cout << setw(6) << i_loc_level << setw(2) << level << setw(5) << j_loc_deg
        << setw(3) << deg_lat << setw(4) << k_loc_deg << setw(3) << deg_lon
        << "  " << setiosflags(ios::left)
        << setw(25) << setfill('.') << name_Value_6 << " = " << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << Value_6 << setw(6)
        << name_unit_mmd << endl << endl;
    choice++;
    if(choice <= 3) goto preparation;
    cout << endl;
    double Value_7 = precipitablewater_average;
    double Value_8 = precipitation_average;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_7 << " = " << resetiosflags(ios::left) << setw(7) << fixed
        << setfill(' ') << Value_7 << setw(6) << name_unit_mm << "   " << setiosflags(ios::left)
        << setw(40) << setfill('.') << name_Value_8 << " = " << resetiosflags(ios::left)
       << setw(7) << fixed << setfill(' ') << Value_8 << setw(6) << name_unit_mma
        << "   " << setiosflags(ios::left) << setw(40) << setfill ('.') << name_Value_9
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_8 / 365. << setw(6) << name_unit_mmd << endl;
    double Value_10 = precipitation_NASA_average;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_7 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_7 << setw(6) << name_unit_mm
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.') << name_Value_10
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_10 << setw(6) << name_unit_mma << "   " << setiosflags(ios::left)
        << setw(40) << setfill('.') << name_Value_11 << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_10 / 365. << setw(6) << name_unit_mmd << endl;
    double Value_13 = Evaporation_Dalton_average;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_7 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_7 << setw(6) << name_unit_mm
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_14 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_13 << setw(6) << name_unit_mma
        << "   " << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_15 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_13 / 365. << setw(6)
        << name_unit_mmd << endl;
    double Value_9 = co2_vegetation_average * co2_0;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_22 << " = " << resetiosflags(ios::left) << setw (7)
        << fixed << setfill(' ') << Value_9 << setw(6) << name_unit_ppm
        << endl << endl << endl;
//    double Value_25 = (GetMean_2D(jm, km, temperature_NASA) - 1.) * t_0;
    double Value_25 = (GetMean_2D(jm, km, temperature_NASA) - 1.);
    double Value_26 = temperature_surf_average;
    double Value_27 = t_average + t_paleo * t_0;
    cout << setw(6) << setiosflags(ios::left) << setw(40) << setfill('.')
        << name_Value_25 << " = " << resetiosflags(ios::left) << setw(7)
        << fixed << setfill(' ') << Value_25 << setw(6) << name_unit_t << "   "
        << setiosflags(ios::left) << setw(40) << setfill('.') << name_Value_26
        << " = " << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_26 << setw(6) << name_unit_t << "   " << setiosflags(ios::left)
        << setw(40) << setfill('.') << name_Value_27 << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << Value_27 << setw(6) << name_unit_t << endl << endl << endl;
}



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


void cAtmosphereModel::CalculateNodeWeights(int jm, int km){
    //use cosine of latitude as weights for now
    //longitudes: 0-360(km) latitudes: 90-(-90)(jm)
    double weight = 0.;
    m_node_weights.clear();
    for(int i=0; i<jm; i++){
        if(i<=90){
            weight = cos((90-i) * M_PI / 180.0);
        }else{
            weight = cos((i-90) * M_PI / 180.0);
        }
        m_node_weights.push_back(std::vector<double>());
        m_node_weights[i].resize(km, weight);
    }
    return;
}


