/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
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

#include "BC_Thermohalin.h"
#include "cHydrosphereModel.h"
#include "Array.h"

using namespace std;
using namespace AtomUtils;

BC_Thermohalin::BC_Thermohalin(int im, int jm, int km, int i_beg, int i_max,
                            int Ma, int Ma_max, int Ma_max_half, double dr, double g,
                            double r_0_water, double u_0, double p_0, double t_0,
                            double c_0, double cp_w, double L_hyd, double t_average,
                            double t_cretaceous_max, double t_equator, double t_pole,
                            const string &input_path){
    this -> im = im;
    this -> jm = jm;
    this -> km = km;
    this -> i_max = i_max;
    this -> i_beg = i_beg;
    this -> Ma = Ma;
    this -> Ma_max = Ma_max;
    this -> Ma_max_half = Ma_max_half;
    this -> dr = dr;
    this -> g = g;
    this -> r_0_water = r_0_water;
    this -> u_0 = u_0;
    this -> p_0 = p_0;
    this -> t_0 = t_0;
    this -> c_0 = c_0;
    this -> cp_w = cp_w;
    this -> L_hyd = L_hyd;
    this -> t_average = t_average;
    this -> t_cretaceous_max = t_cretaceous_max;
    this -> t_equator = t_equator;
    this -> t_pole = t_pole;
    this -> input_path = input_path;

    j_half = ( jm - 1 ) / 2;
}

BC_Thermohalin::~BC_Thermohalin(){}


void cHydrosphereModel::IC_v_w_EkmanSpiral(){
    float water_wind = .03;  // ocean surface velocity is about 3% of the wind velocity at the surface
    //Ekman spiral demands 45° turning of the water flow compared to the air flow at contact surface
    //a further turning downwards until the end of the shear layer such that finally 90° of turning are reached
//    float Ekman_angle = 45.0 / pi180;
    float Ekman_angle = 0.0;
    // initial conditions for v and w velocity components at the sea surface
    // ocean surface velocity is about 3% of the wind velocity at the surface
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_water(h, im-1, j, k)){
                v.x[im - 1][j][k] = water_wind * v.x[im - 1][j][k] / u_0;
                w.x[im - 1][j][k] = water_wind * w.x[im - 1][j][k] / u_0;
            }
        }
    }
    // surface wind vector driving the Ekman spiral in the Ekman layer
    // northern and southern hemisphere
    int i_Ekman = 0;
    double sinthe = 0.;
    double vel_magnit = 0.;
    double vel_mag = 0.;
    double i_Ekman_layer = 0.;
    double coeff = 0.;
    double gam_z = 0.;
    double exp_gam_z = 0.;
    double sin_gam_z = 0;
    double cos_gam_z = 0;
    double alfa = 0;
    float angle = 0;
    for(int j = 1; j < jm-1; j++){
        for(int k = 1; k < km-1; k++){
            sinthe = sin( fabs( the.z[j] - M_PI / 2. ) );
            vel_magnit = sqrt( v.x[im-1][j][k] * v.x[im-1][j][k] 
                         + w.x[im-1][j][k] * w.x[im-1][j][k] ) 
                         / water_wind * u_0; // dimensional surface wind velocity U_10
//          original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 139, eq. 9.16
            i_Ekman_layer = 3. * 7.6 / sqrt( sinthe ) * vel_magnit;
            // assumed depth at i_Ekman = 3 x pi / gam, velocity opposit plus pi to surface velocity
            coeff = i_Ekman_layer / L_hyd;
            if(coeff >= 1.) coeff = 1.;
            i_Ekman = ( im - 1 ) * ( 1. - coeff );

//          prevention of singular values
            if( w.x[im-1][j][k] == 0. ) w.x[im-1][j][k] = 1.e-6;
            if( v.x[im-1][j][k] == 0. ) v.x[im-1][j][k] = 1.e-6;
            if( vel_magnit == 0. ) vel_magnit = 1.e-6;

//          turning the wind velocities to water velocities along the surface, local angle +/- Ekman_angle
//            alfa = atan( fabs( v.x[im-1][j][k] / w.x[im-1][j][k] ));
            if(j <= (jm-1)/2){
                alfa = atan( fabs( v.x[im-1][j][k] / w.x[im-1][j][k] ));
                if(( w.x[im-1][j][k] >= 0. ) && ( v.x[im-1][j][k] >= 0. )){
                    angle = alfa - Ekman_angle;
                }
                if(( w.x[im-1][j][k] <= 0. ) && ( v.x[im-1][j][k] >= 0. )){
                    angle = .5 * M_PI - ( alfa + Ekman_angle ) + .5 * M_PI;
                }
                if(( w.x[im-1][j][k] <= 0. ) && ( v.x[im-1][j][k] <= 0. )){
                    angle = alfa - Ekman_angle + M_PI;
                }
                if(( w.x[im-1][j][k] >= 0. ) && ( v.x[im-1][j][k] <= 0. )){
                    angle = .5 * M_PI - ( alfa + Ekman_angle ) + 1.5 * M_PI;
                }
            }else{
                alfa = - atan( fabs( v.x[im-1][j][k] / w.x[im-1][j][k] ));
                if(( w.x[im-1][j][k] >= 0. ) && ( v.x[im-1][j][k] >= 0. )){
                    angle = alfa + Ekman_angle;
                }
                if(( w.x[im-1][j][k] <= 0. ) && ( v.x[im-1][j][k] >= 0. )){
                    angle = .5 * M_PI - ( alfa - Ekman_angle ) + .5 * M_PI;
                }
                if(( w.x[im-1][j][k] <= 0. ) && ( v.x[im-1][j][k] <= 0. )){
                    angle = alfa + Ekman_angle + M_PI;
                }
                if(( w.x[im-1][j][k] >= 0. ) && ( v.x[im-1][j][k] <= 0. )){
                    angle = .5 * M_PI - ( alfa - Ekman_angle ) + 1.5 * M_PI;
                }
            }
//          original law in Robert H. Stewart, Introduction to Physical Oceanography, p. 138, eq. 9.11a/b
            v.x[im-1][j][k] = v.x[im-1][j][k] * water_wind / u_0; // non-dimensional surface water velocity V_0
            w.x[im-1][j][k] = w.x[im-1][j][k] * water_wind / u_0; // non-dimensional surface water velocity V_0
            vel_mag = sqrt( v.x[im-1][j][k] * v.x[im-1][j][k] 
                          + w.x[im-1][j][k] * w.x[im-1][j][k] );
            v.x[im-1][j][k] = vel_mag * cos( angle );
            w.x[im-1][j][k] = - vel_mag * sin( angle );
            for(int i = im-1; i >= i_Ekman; i--){
//            for(int i = im-2; i >= i_Ekman; i--){
//          original laws in Robert H. Stewart, Introduction to Physical Oceanography, p. 137, eq. 9.9a/b
                gam_z = - 3. * M_PI * (double)( i - im-1 ) 
                    / (double)( i_Ekman - im-1 );
                exp_gam_z = exp ( gam_z );
                sin_gam_z = - sin( angle + gam_z );
                cos_gam_z = cos( angle + gam_z );
                if(j <= (jm-1)/2){
                    v.x[i][j][k] = vel_mag * exp_gam_z * sin_gam_z; 
                    w.x[i][j][k] = - vel_mag * exp_gam_z * cos_gam_z;
                }else{
                    v.x[i][j][k] = - vel_mag * exp_gam_z * sin_gam_z; 
                    w.x[i][j][k] = - vel_mag * exp_gam_z * cos_gam_z;
                }
//    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   i = " << i << "   j = " << j << "   k = " << k << "   gam_z = " << gam_z << "   exp_gam_z = " << exp_gam_z << "   sin_gam_z = " << sin_gam_z << "   cos_gam_z = " << cos_gam_z << "   i_Ekman = " << i_Ekman << "   v = " << v.x[i][j][k] << "   w = " << w.x[i][j][k] << "   vel_mag = " << vel_mag << "   vel_magnit = " << vel_magnit << endl;
            }
        }
    }
    for(int i = 0; i <= im-1; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                if(is_land( h, i, j, k)){
                    v.x[i][j][k] = 0.;
                    w.x[i][j][k] = 0.;
                }
            }
        }
    }
}


void BC_Thermohalin::BC_Temperature_Salinity(Array &h, Array &t, Array &c, 
    Array &p_dyn, Array_2D &Evaporation_Dalton, Array_2D &Precipitation, 
    Array_2D &salinity_evaporation){
    // initial conditions for salt content and temperature and salinity decrease below the sea surface
    j_half = ( jm -1 ) / 2;
    j_max = jm - 1;
    d_j_half = (double)j_half;
    d_j_max = (double)j_max;
    // temperature-distribution by Ruddiman approximated by a parabola
    t_coeff = t_pole - t_equator;
    t_cretaceous_coeff = t_cretaceous_max / ( (double)Ma_max_half -
        (double)( Ma_max_half * Ma_max_half / Ma_max ) );   // in °C
    t_cretaceous = t_cretaceous_coeff * (double)( - ( Ma * Ma ) / Ma_max + Ma );   // in °C
    if(Ma == 0)  t_cretaceous = 0.;
    t_cretaceous = 0.;
    cout.precision(3);
    time_slice_comment = "      time slice of Cretaceous-AGCM:";
    time_slice_number = " Ma = ";
    time_slice_unit = " million years";
    cout << endl << setiosflags(ios::left) << setw (50) << setfill('.')
        << time_slice_comment << setw(6) << std::fixed << setfill(' ')
        << time_slice_number << setw(3) << Ma << setw(12) << time_slice_unit 
        << endl << endl;
    temperature_comment = "      temperature increase at cretaceous times: ";
    temperature_gain = " t increase";
    temperature_modern = "      mean temperature at modern times: ";
    temperature_cretaceous = "      mean temperature at cretaceous times: ";
    temperature_average = " t modern";
    temperature_average_cret = " t cretaceous";
    temperature_unit =  "°C ";
    cout << endl << setiosflags(ios::left) << setw(50) << setfill('.') 
        << temperature_comment << resetiosflags(ios::left) << setw(12) 
        << temperature_gain << " = " << setw(7) << setfill(' ') 
        << t_cretaceous << setw(5) << temperature_unit << endl << setw(50) 
        << setfill('.') << setiosflags(ios::left) << temperature_modern 
        << resetiosflags(ios::left) << setw(13) << temperature_average 
        << " = " << setw(7) << setfill(' ') << t_average << setw(5) 
        << temperature_unit << endl << setw(50) << setfill('.') 
        << setiosflags(ios::left) << temperature_cretaceous 
        << resetiosflags(ios::left) << setw(13) << temperature_average_cret 
        << " = " << setw(7) << setfill(' ') << t_average + t_cretaceous 
        << setw(5) << temperature_unit << endl;
    c_average = ( t_average + 346. ) / 10.;// in psu, relation taken from "Ocean Circulation, The Open University"
    c_cretaceous = ( t_average + t_cretaceous + 346. ) / 10.;// in psu
    c_cretaceous = c_cretaceous - c_average;
    if(Ma == 0)  c_cretaceous = 0.;
    salinity_comment = "      salinity increase at cretaceous times: ";
    salinity_gain = " c increase";
    salinity_modern = "      mean salinity at modern times: ";
    salinity_cretaceous = "      mean salinity at cretaceous times: ";
    salinity_average = " c modern";
    salinity_average_cret = " c cretaceous";
    salinity_unit =  "psu ";
    cout << endl << setiosflags(ios::left) << setw(50) << setfill('.') 
        << salinity_comment << resetiosflags(ios::left) << setw(12) 
        << salinity_gain << " = " << setw(7) << setfill(' ') << c_cretaceous 
        << setw(5) << salinity_unit << endl << setw(50) << setfill('.') 
        << setiosflags(ios::left) << salinity_modern << resetiosflags(ios::left) 
        << setw(13) << salinity_average  << " = " << setw(7) << setfill(' ') 
        << c_average << setw(5) << salinity_unit << endl << setw(50) 
        << setfill('.') << setiosflags(ios::left) << salinity_cretaceous 
        << resetiosflags(ios::left) << setw(13) << salinity_average_cret 
        << " = " << setw(7) << setfill(' ') << c_average + c_cretaceous 
        << setw(5) << salinity_unit << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            if(is_water( h, im-1, j, k)){
                d_j = ( double ) j;
                t_Celsius = t.x[im-1][j][k] * t_0 - t_0;
                if(t_Celsius <= t_pole * t_0 - t_0)
                    t_Celsius = t_pole * t_0 - t_0;
                c.x[im-1][j][k] = ( ( t_Celsius + 346. ) / 10. ) / c_0;
                if ( c.x[im-1][j][k] <= 1. ) c.x[im-1][j][k] = 1.;
            }else{
                c.x[im-1][j][k] = 1.;
            }
        }
    }
    double tm_tbeg = 0.;
    double cm_cbeg = 0.;
    double d_i_max = (double)i_max;
    double d_i_beg = (double)i_beg;
    double d_i = 0.;
// distribution of t and c with increasing depth
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            tm_tbeg = ( t.x[im-1][j][k] - 1. ) 
                / ( d_i_max * d_i_max - d_i_beg * d_i_beg );
            if(t.x[im-1][j][k] <= t_pole)
                t.x[im-1][j][k] = t_pole;
            cm_cbeg = ( c.x[im-1][j][k] - 1. ) 
                / ( d_i_max * d_i_max - d_i_beg * d_i_beg );
            for(int i = i_beg; i < im-1; i++){
                    d_i = (double)i;
                    t.x[i][j][k] = 1. + tm_tbeg 
                        * ( d_i * d_i - d_i_beg * d_i_beg );// parabolic approach
                    if(t.x[i][j][k] <= t_pole)
                        t.x[i][j][k] = t_pole;
                    c.x[i][j][k] = 1. + cm_cbeg 
                        * ( d_i * d_i - d_i_beg * d_i_beg );// parabolic approach
            }
        }
    }
// preparations for salinity increase due to the differences between evaporation and precipitation 
// procedure given in Rui Xin Huang, Ocean Circulation, p. 165
// amount of additional salinity by the difference of evaporation and precipitation is negligible but functional
    double coeff_salinity = 1.1574e-8 * L_hyd / ( c_0 * u_0 );  // 1.1574-8 is the conversion from (Evap-Prec) in mm/d to mm/s
    double evap_precip = 0.;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            std::vector<std::vector<double> > c_fix(jm, std::vector<double>(km, 0));
            c_fix[j][k] = c.x[im-1][j][k];
            double salinity_surface_n = 0.;
            double salinity_surface = 0.;
                for(int iter_prec = 1; iter_prec <= 20; iter_prec++){// iter_prec may be varied
                    evap_precip = Evaporation_Dalton.y[j][k] 
                        - Precipitation.y[j][k];
//                salinity_surface = - ( - 3. * c.x[im - 1][j][k]
//                               + 4. * c.x[im - 2][j][k] 
//                               - c.x[im - 3][j][k] ) / ( 2. * dr ) // 1. order derivative, 2. order accurate
//                               * ( 1. - 2. * c.x[im - 1][j][k] ) 
//                               * evap_precip; 
                    salinity_surface = ( c.x[im-1][j][k] 
                        - c.x[im-2][j][k] ) / dr // 1. order derivative, 1. order accurate
                        * ( 1. - 2. * c.x[im-1][j][k] ) * evap_precip;
                    if(iter_prec == 1) 
                        salinity_surface_n = .9 * salinity_surface;
                    salinity_evaporation.y[j][k] = coeff_salinity 
                        * salinity_surface;
                    if(is_land(h, im-1, j, k))
                        salinity_evaporation.y[j][k] = 0.;
                    c.x[im-1][j][k] = c_fix[j][k] 
                        + salinity_evaporation.y[j][k];
                    double cm_cbeg = 0.;
                    double d_i_max = (double)i_max;
                    double d_i_beg = (double)i_beg;
                    double d_i = 0.;
                    cm_cbeg = ( c.x[im-1][j][k] - 1. ) 
                        / ( d_i_max * d_i_max - d_i_beg * d_i_beg );
                    for(int i = i_beg; i < im-1; i++){
                        d_i = (double)i;
                        c.x[i][j][k] = 1. + cm_cbeg 
                            * ( d_i * d_i - d_i_beg * d_i_beg );// parabolic approach
                    }
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   iter_prec = " << iter_prec << "   salinity_evaporation = " << salinity_evaporation.y[j][k] * c_0 << "   coeff_salinity = " << coeff_salinity << "   salinity_surface = " << salinity_surface * c_0 << "   salinity_surface_n = " << salinity_surface_n * c_0 << "   c_fix = " << c_fix[j][k] * c_0 << "   c = " << c.x[im-1][j][k] * c_0 << "   Evap-Prec = " << evap_precip << "   Evap = " << Evaporation_Dalton.y[j][k] << "   Prec = " << Precipitation.y[j][k] << "   c_grad_1 = " << ( c.x[im-1][j][k] - c.x[im-2][j][k] ) / dr << "   c_grad_2 = " << - ( - 3. * c.x[im-1][j][k] + 4. * c.x[im-2][j][k] - c.x[im-3][j][k] ) / ( 2. * dr ) << endl;
*/
                    if(fabs(salinity_surface - salinity_surface_n) 
                        < 1.e-5 * salinity_surface_n)  break;  
                    salinity_surface_n = salinity_surface;
            } // iter_prec 
        } // end j
    } // end k
    for(int i = 0; i < i_beg; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                if(is_water(h, i, j, k)){
                    t.x[i][j][k] = 1.;
                    c.x[i][j][k] = 1.;
                }
            }
        }
    }
}




void BC_Thermohalin::BC_Pressure_Density(Array &p_stat, Array &r_water,
                     Array &r_salt_water, Array &t, Array &c, Array &h){
// hydrostatic pressure, equations of state for water and salt water density
// as functions of salinity, temperature and hydrostatic pressure
    double t_Celsius_0 = 0.;
    double t_Celsius_1 = 0.;
    double p_km = 0.;
    double C_p = 0.;
    double beta_p =  0.;
    double alfa_t_p =  0.;
    double gamma_t_p =  0.;
    double E_water = 2.15e9;                                // given in N/m²
    double beta_water = 8.8e-5;                         // given in m³/( m³ * °C)
    double r_air = 1.2041;                                      // given in kg/m³
    double R_Air = 287.1;                                       // given in J/( kg*K )
// hydrostatic pressure, water and salt water density at the surface
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            p_stat.x[im-1][j][k] =  .01 * ( r_air * R_Air * t.x[im-1][j][k] * t_0 ) / 1000.;
            // given in bar, isochoric approach, constant air density at the surface
            r_water.x[im-1][j][k] = r_0_water;                // given in kg/m³
            t_Celsius_1 = t.x[im-1][j][k] * t_0 - t_0;
            p_km = 0.;
            C_p = 999.83;
            beta_p = .808;
            alfa_t_p = .0708 * ( 1. + .068 * t_Celsius_1 );
            gamma_t_p = .003 * ( 1. - .012 * t_Celsius_1 );
            r_salt_water.x[im-1][j][k] = C_p + beta_p * c.x[im-1][j][k] * c_0 -
                alfa_t_p * t_Celsius_1 - gamma_t_p * ( c_0 - c.x[im-1][j][k] * c_0 ) * t_Celsius_1;
        }
    }
// hydrostatic pressure, water and salt water density in the field
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = im-2; i >= 0; i--){
                d_i = (double)( im - 1 - i );
                t_Celsius_1 = t.x[i][j][k] * t_0 - t_0;
                t_Celsius_0 = t.x[i + 1][j][k] * t_0 - t_0;
//              p_stat.x[i][j][k] = r_0_water * g * d_i * ( L_hyd / ( double ) ( im-1 ) )
//                   100000. + p_0 / 1000.;                // hydrostatic pressure in bar
                p_stat.x[i][j][k] = r_water.x[i + 1][j][k] * g * d_i * ( L_hyd /
                    (double)( im-1 ) ) / 100000. + p_0 / 1000.;             // hydrostatic pressure in bar
                r_water.x[i][j][k] = r_water.x[i + 1][j][k] / ( 1. + beta_water *
                    ( t_Celsius_1 - t_Celsius_0 ) ) / ( 1. - ( p_stat.x[i][j][k] -
                    p_stat.x[i+1][j][k] ) / E_water * 1e5 );               // given in kg/m³
                p_km  = (double)( im - 1 - i ) * ( L_hyd / (double)( im-1 ) ) / 1000.;
                C_p = 999.83 + 5.053 * p_km - .048 * p_km * p_km;
                beta_p = .808 - .0085* p_km;
                alfa_t_p = .0708 * ( 1. + .351 * p_km + .068 *( 1. - .0683 * p_km ) * t_Celsius_1 );
                gamma_t_p = .003 * ( 1. - .059 * p_km - .012 * ( 1. - .064 * p_km ) * t_Celsius_1 );
                r_salt_water.x[i][j][k] = C_p + beta_p * c.x[i][j][k] * c_0 -
                    alfa_t_p * t_Celsius_1 - gamma_t_p * ( c_0 - c.x[i][j][k] * c_0 ) * t_Celsius_1;
            }
        }
    }
}






void BC_Thermohalin::BC_Surface_Temperature_NASA
    (const string &Name_SurfaceTemperature_File, Array &t){
    cout.precision(3);
    cout.setf(ios::fixed);
    ifstream Name_SurfaceTemperature_File_Read;
    Name_SurfaceTemperature_File_Read.open(Name_SurfaceTemperature_File);
    if(!Name_SurfaceTemperature_File_Read.is_open()){
        cerr << "ERROR: could not open SurfaceTemperature_File file at " << Name_SurfaceTemperature_File << "\n";
        abort();
    }
    j = 0;
    k = 0;
    while((k < km) && (!Name_SurfaceTemperature_File_Read.eof())){
        while(j < jm){
            Name_SurfaceTemperature_File_Read >> dummy_1;
            Name_SurfaceTemperature_File_Read >> dummy_2;
            Name_SurfaceTemperature_File_Read >> dummy_3;
            t.x[im-1][j][k] = ( dummy_3 + 273.15 ) / 273.15;
            j++;
        }
        j = 0;
        k++;
    }
    Name_SurfaceTemperature_File_Read.close();
    for(int j = 0; j < jm; j++){
        for(int k = 1; k < km-1; k++){
            if(k == 180) t.x[im-1][j][k] = ( t.x[im-1][j][k + 1] 
                + t.x[im-1][j][k - 1] ) * .5;
        }
    }
}






void BC_Thermohalin::BC_Surface_Salinity_NASA(const string &Name_SurfaceSalinity_File, Array &c){
    // initial conditions for the salinity at the sea surface
    streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;
    cout.precision(3);
    cout.setf(ios::fixed);
    ifstream Name_SurfaceSalinity_File_Read;
    Name_SurfaceSalinity_File_Read.open(Name_SurfaceSalinity_File);
    if(!Name_SurfaceSalinity_File_Read.is_open()){
        cerr << "ERROR: could not open SurfaceSalinity_File file at " << Name_SurfaceSalinity_File << "\n";
        abort();
    }
    j = 0;
    k = 0;
    while((k < km) && (!Name_SurfaceSalinity_File_Read.eof()) ){
        while(j < jm){
            Name_SurfaceSalinity_File_Read >> dummy_1;
            Name_SurfaceSalinity_File_Read >> dummy_2;
            Name_SurfaceSalinity_File_Read >> dummy_3;
            if(dummy_3 < 0.) dummy_3 = 1.;
            else  c.x[im-1][j][k] = dummy_3 / c_0;
            j++;
        }
        j = 0;
        k++;
    }
    Name_SurfaceSalinity_File_Read.close();
}

void BC_Thermohalin::IC_Equatorial_Currents 
                     (Array &h, Array &u, Array &v, Array &w){
// currents along the equator
// equatorial undercurrent - Cromwell flow, EUC
// equatorial intermediate current, EIC
// nothern and southern equatorial subsurface counter-currents, NSCC und SSCC
// nothern and southern equatorial counter-currents, NECC und SECC
    double IC_water = 1.;  // no dimension, ( average velocity compares to   u_0 * IC_water = 0,25 )
/*
// one grid step compares to a depth of 25 m for L_hyd = 1000m   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    i_EIC_u = i_beg; // 1000m depth  // westward equatorial intermediate current, max 0.1 m/s
    i_EIC_o = 28; // 0m depth
    i_SCC_u = 8; // 800m depth  // eastward subsurface counter current, max 0.05 m/s
    i_SCC_o = 28; // 300m depth
    i_ECC_u = 28; // 400m depth  // eastward equatorial counter current, max 0.2 m/s
    i_ECC_o = im; // 0m depth
    i_EUC_u = 32; // 200m depth  // eastward equatorial under current ( Cromwell current ), max 0.8 m/s
    i_EUC_o = 36; // 100m depth
*/
// one grid step compares to a depth of 5 m for L_hyd = 200m   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    i_EIC_u = i_beg; // 200m depth  // westward equatorial intermediate current, max 0.1 m/s
    i_EIC_o = 0; // 0m depth
    i_SCC_u = 0; // 800m depth  // eastward subsurface counter current, max 0.05 m/s
    i_SCC_o = 0; // 300m depth
    i_ECC_u = 0; // 200m depth  // eastward equatorial counter current, max 0.2 m/s
    i_ECC_o = im; // 0m depth
    i_EUC_u = 0; // 200m depth  // eastward equatorial under current ( Cromwell current ), max 0.8 m/s
    i_EUC_o = 20; // 100m depth
// equatorial currents and counter-currents
//  §§§§§§§§§§§§§§§§§§§§§§§§§   valid for all paleo-ocean constellations along the equator   §§§§§§§§§§§§§§§§§§§§§§§§§§
    j_half = ( jm -1 ) / 2;
    k_beg = 0;
    k_end = 0;
// extention of land and ocean areas
    for(int k = 1; k < km-1; k++){
        if( k < km ){
            if((is_water(h, im-1, j_half, k)) 
                && (is_land(h, im-1, j_half, k+1)) ){
                k_end = k;
            }
            if((is_land(h, im-1, j_half, k)) 
                && (is_water(h, im-1, j_half, k+1)) ){
                k_beg = k;
            }
            if(((is_water(h, im-1, j_half, k)) 
                && (is_water(h, im-1, j_half, k+1)) ) 
                && ( k == km-2 ) ){
                k_end = k;
            }
// equatorial northern counter-current ( NECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial northern counter-current ( from j=83 until j=88 compares to 3°N until 8°N )
            for(int i = i_ECC_u; i < i_ECC_o; i++){
                for(int j =83; j < 88; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = + .01 * IC_water 
                                * (double)( i - i_ECC_u ) 
                                / (double)( i_ECC_o - i_ECC_u );
                        }
                    }
                }
            }
// equatorial southern counter-current ( SECC, i=im-1 until i=im compares to 0 until -200m depth )
// equatorial southern counter-current ( from j=93 until j=96 compares to 3°S until 6°S )
            for(int i = i_ECC_u; i < i_ECC_o; i++){
                for(int j = 93; j < 98; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = + .01 * IC_water 
                                * (double)( i - i_ECC_u ) 
                                / (double)( i_ECC_o - i_ECC_u );
                        }
                    }
                }
            }
// equatorial undercurrent - Cromwell current ( EUC, i=im-2 until i=im-1 compares to -100 until -200m depth )
// equatorial undercurrent - Cromwell current ( from j=87 until j=93 compares to 3°N until 3°S )
            for(int i = i_EUC_u; i < i_EUC_o; i++){
                for(int j = 87; j < 94; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = + .4 * IC_water;
                        }
                    }
                }
            }
// equatorial intermediate current ( EIC, i=im-4 until i=im-2 compares to -300 until -1000m depth )
// equatorial intermediate current ( from j=88 until j=92 compares to 2°N until 2°S )
            for(int i = i_EIC_u; i < i_EIC_o; i++){
                for(int j = 88; j < 93; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = - .05 * IC_water 
                                * (double)( i - i_EIC_u ) 
                                / (double)( i_EIC_o - i_EIC_u );
                        }
                    }
                }
            }
// equatorial northern and southern subsurface counter-current
// equatorial northern subsurface counter-current ( NSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial northern subsurface counter-current ( from j=86 until j=87 compares to 3°N until 4°N )
            for(int i = i_SCC_u; i < i_SCC_o; i++){
                for(int j = 86; j < 88; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = .025 * IC_water 
                                * (double)( i - i_SCC_u ) 
                                / (double)( i_SCC_o - i_SCC_u );
                        }
                    }
                }
            }
// equatorial southern subsurface counter-current ( SSCC, i=im-3 until i=im-2 compares to -300 until -800m depth )
// equatorial southern subsurface counter-current ( from j=93 until j=94 compares to 3°S until 4°S )
            for(int i = i_SCC_u; i < i_SCC_o; i++){
                for(int j = 93; j < 95; j++){
                    for(int k = k_beg; k <= k_end; k++){
                        if(is_water(h, i, j, k)){
                            v.x[i][j][k] = 0.;
                            w.x[i][j][k] = .025 * IC_water 
                                * (double)( i - i_SCC_u ) 
                                / (double)( i_SCC_o - i_SCC_u );
                        }
                    }
                }
            }
        }
    }
}






void BC_Thermohalin::IC_CircumPolar_Current 
                     (Array &h, Array &u, Array &v, Array &w, Array &c){
// south polar sea
// antarctic circumpolar current ( -5000m deep ) ( from j=147 until j=152 compares to 57°S until 62°S,
// from k=0 until k=km compares to 0° until 360° )
    for(int i = i_beg; i < im; i++){
        for(int j = 147; j < 153; j++){
            for(int k = 0; k < km; k++){
                if(is_water(h, i, j, k)){
//                  c.x[i][j][k] = 1.;
                    w.x[i][j][k] = 2. * u_0;  // 0.5 m/s
                }
            }
        }
    }
}




void BC_Thermohalin::IC_u_WestEastCoast
                     (Array_1D &rad, Array &h, Array &u, Array &v, 
                     Array &w, Array &un, Array &vn, Array &wn){
// initial conditions for v and w velocity components at the sea surface close to east or west coasts
// reversal of v velocity component between north and south equatorial current ommitted at respectively 10°
// w component unchanged
    int i_beg = 20;  // == 500m depth
    int i_middle = i_beg + 16;  // 0 + 36 = 36 for total depth 100m ( i_beg= 0 ), asymmetric with depth for stepsizes of 25m
    int i_half = 38;
    d_i_middle = (double)i_middle;
    d_i_half = (double)i_half;
    d_i_max =  (double)i_max;
// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included
// northern hemisphere: east coast
//    k_grad = 6;                                                       // extension of velocity change
    k_grad = 10;                                                        // extension of velocity change
    k_a = k_grad;                                                       // left distance
    k_b = 0;                                                            // right distance
    k_water = 0;                                                        // on water closest to coast
    k_sequel = 1;                                                       // on solid ground
    for(int j = 0; j < 91; j++){                                        // outer loop: latitude
        for(int k = 0; k < km; k++){                                    // inner loop: longitude
            if(is_land(h, i_half, j, k)) k_sequel = 0;                 // if solid ground: k_sequel = 0
            if((is_water(h, i_half, j, k)) && (k_sequel == 0)) 
                k_water = 0;    // if water and and k_sequel = 0 then is water closest to coast
            else k_water = 1;                                           // somewhere on water
            if((is_water(h, i_half, j, k)) && (k_water == 0)){     // if water is closest to coast, change of velocity components begins
                for(int l = 0; l < k_grad; l++){// extension of change, sign change in v-velocity and distribution of u-velocity with depth
                    if(k+l > km-1) break;
                    for(int i = i_beg; i < i_half; i++){                // loop in radial direction, extension for u -velocity component, downwelling here
                        m = i + ( i_half - i_middle );
                        d_i = ( double ) i;
                        u.x[i][j][k + l] = - d_i / d_i_half * water_wind 
                            / ( ( double )( l + 1 ) );    // increase with depth, decrease with distance from coast
                        u.x[m][j][k + l] = - d_i / d_i_middle * water_wind 
                            / ( ( double )( l + 1 ) );// decrease with depth, decrease with distance from coast
                    }
                }
                k_sequel = 1;                                           // looking for another east coast
            }
        }                                                               // end of longitudinal loop
        k_water = 0;                                                    // starting at another latitude
    }                                                                   // end of latitudinal loop
// southern hemisphere: east coast
    k_water = 0;
    k_sequel = 1;
    for(int j = 91; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_land(h, i_half, j, k)) k_sequel = 0;
            if((is_water(h, i_half, j, k)) && (k_sequel == 0)) 
                k_water = 0;
            else k_water = 1;
            if((is_water(h, i_half, j, k)) && (k_water == 0)){
                for(int l = 0; l < k_grad; l++){
                    if(k+l > km-1) break; 
                    for(int i = i_beg; i < i_half; i++){
                        m = i + ( i_half - i_middle );
                        d_i = (double)i;
                        u.x[i][j][k + l] = - d_i / d_i_half * water_wind 
                            / ( (double)( l + 1 ) );    // increase with depth, decrease with distance from coast
                        u.x[m][j][k + l] = - d_i / d_i_middle * water_wind 
                            / ( (double)( l + 1 ) );// decrease with depth, decrease with distance from coast
                    }
                }
                k_sequel = 1;
            }
        }
        k_water = 0;
    }
// search for east coasts and associated velocity components to close the circulations
// transition between coast flows and open sea flows included
// northern hemisphere: west coast
//    k_grad = 6;                                                       // extension of velocity change
    k_a = 0;                                                            // left distance
    k_water = 0;                                                        // somewhere on water
    flip = 0;                                                           // somewhere on water
    for(int j = 0; j < 91; j++){                                        // outer loop: latitude
        for(int k = 0; k < km; k++){                                    // inner loop: longitude
            if(is_water(h, i_half, j, k)){                              // if somewhere on water
                k_water = 0;                                            // somewhere on water: k_water = 0
                flip = 0;                                               // somewhere on water: flip = 0
            }
            else k_water = 1;                                           // first time on land
            if((flip == 0) && (k_water == 1)){                      // on water closest to land
                for(int l = k; l > ( k - k_grad + 1 ); l--){            // backward extention of velocity change: nothing changes
                    if(l < 0) break;
                    for(int i = i_beg; i < i_half; i++){                // loop in radial direction, extension for u -velocity component, downwelling here
                        m = i + ( i_half - i_middle );
                        d_i = (double)i;
                        u.x[i][j][l] = + d_i / d_i_half * water_wind 
                            / ( (double)( k - l + 1 ) );    // increase with depth, decrease with distance from coast
                        u.x[m][j][l] = + d_i / d_i_middle * water_wind 
                            / ( (double)( k - l + 1 ) );    // decrease with depth, decrease with distance from coast
                    }
                }
                flip = 1;
            }
        }
        flip = 0;
    }
// southern hemisphere: west coast
    k_water = 0;
    flip = 0;
    for(int j = 91; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_water(h, i_half, j, k)){
                k_water = 0;
                flip = 0;
            }
            else k_water = 1;
            if((flip == 0) && (k_water == 1)){
                for(int l = k; l > ( k - k_grad + 1 ); l--){
                    if(l < 0) break;    
                    for(int i = i_beg; i < i_half; i++){
                        m = i + ( i_half - i_middle );
                        d_i = (double)i;
                        u.x[i][j][l] = + d_i / d_i_half * water_wind 
                            / ( (double)( k - l + 1 ) );
                        u.x[m][j][l] = + d_i / d_i_middle * water_wind 
                            / ( (double)( k - l + 1 ) );
                    }
                }
                flip = 1;
            }
        }
        flip = 0;
    }
    for(int i = 0; i < i_beg; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                u.x[i][j][k] = 0.;
            }
        }
    }
}


void BC_Thermohalin::Value_Limitation_Hyd(Array &h, Array &u, Array &v,
                                    Array &w, Array &p_dyn, Array &t, Array &c){
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if( u.x[i][j][k] >= 11.11 )  u.x[i][j][k] = 11.11;
                if( u.x[i][j][k] <= - 11.11 )  u.x[i][j][k] = - 11.11;
                if( v.x[i][j][k] >= .552 )  v.x[i][j][k] = .552;
                if( v.x[i][j][k] <= - .552 )  v.x[i][j][k] = - .552;
                if( w.x[i][j][k] >= .552 )  w.x[i][j][k] = .552;
                if( w.x[i][j][k] <= - .552 )  w.x[i][j][k] = - .552;
                if( t.x[i][j][k] >= 1.147 )  t.x[i][j][k] = 1.147; //40.15 °C
                if( t.x[i][j][k] <= 0.9963 )  t.x[i][j][k] = 0.9963;// -1.0 °C
                if( c.x[i][j][k] >= 1.156 )  c.x[i][j][k] = 1.156;  // 40.0 psu
                if( c.x[i][j][k] <= .95 )  c.x[i][j][k] = .95;      // 32.0 psu
            }
        }
    }
}



void BC_Thermohalin::Pressure_Limitation_Hyd(Array &p_dyn, Array &p_dynn){
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(p_dyn.x[i][j][k] >= .1)  p_dyn.x[i][j][k] = .1;
                if(p_dyn.x[i][j][k] <= - .1)  p_dyn.x[i][j][k] = - .1;
                p_dynn.x[i][j][k] = p_dyn.x[i][j][k];
            }
        }
    }
}

