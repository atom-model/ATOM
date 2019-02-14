/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
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

#include "BC_Thermo.h"
#include "Array.h"
#include "Array_2D.h"
#include "cAtmosphereModel.h"
#include "Utils.h"
#include "AtomMath.h"

using namespace std;
using namespace AtomUtils;

double get_pole_temperature(int Ma, const std::map<int, double> &pole_temp_map);

BC_Thermo::BC_Thermo (cAtmosphereModel* model, int im, int jm, int km, Array& h): 
        m_model(model),
        im(im),
        jm(jm),
        km(km),
        h(h),
        i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0)))
{
    this-> tropopause_equator = model->tropopause_equator;
    this-> tropopause_pole = model->tropopause_pole;
    this-> L_atm = model->L_atm;
    this-> dt = model->dt;
    this-> dr = model->dr;
    this-> dthe = model->dthe;
    this-> dphi = model->dphi;
    this-> RadiationModel = model->RadiationModel;
    this-> NASATemperature = model->NASATemperature;
    this-> sun = model->sun;
    this-> g = model->g;
    this-> ep = model->ep;
    this-> hp = model->hp;
    this-> u_0 = model->u_0;
    this-> p_0 = model->p_0;
    this-> t_0 = model->t_0;
    this-> sigma = model->sigma;
    this-> albedo_equator = model->albedo_equator;
    this-> albedo_pole = model->albedo_pole;
    this-> gam = model->gam;
    this-> lv = model->lv;
    this-> ls = model->ls;
    this-> cp_l = model->cp_l;
    this-> r_air = model->r_air;
    this-> R_Air = model->R_Air;
    this-> r_water_vapour = model->r_water_vapour;
    this-> R_WaterVapour = model->R_WaterVapour;
    this-> co2_cretaceous = model->co2_cretaceous;
    this-> co2_vegetation = model->co2_vegetation;
    this-> co2_ocean = model->co2_ocean;
    this-> co2_land = model->co2_land;
    this-> co2_factor = model->co2_factor;
    this-> rad_equator = model->rad_equator;
    this-> rad_pole = model->rad_pole;
    this-> epsilon_pole = model->epsilon_pole;
    this-> epsilon_tropopause = model->epsilon_tropopause;
    this-> epsilon_equator = model->epsilon_equator;
    this-> c_tropopause = model->c_tropopause;
    this-> co2_tropopause = model->co2_tropopause;
    this-> c_ocean = model->c_ocean;
    this-> t_average = model->t_average;
    this-> co2_average = model->co2_average;
    this-> co2_pole = model->co2_pole;
    this-> co2_equator = model->co2_equator;
    this-> t_tropopause = model->t_tropopause;
    this-> t_equator = model->t_equator;
    this-> t_pole = model->t_pole;
    this-> declination = model->declination;
    this-> sun_position_lat = model->sun_position_lat;
    this-> sun_position_lon = model->sun_position_lon;

    im_tropopause = model->get_tropopause();

    Ma = int(round(*m_model->get_current_time()));

    t_cretaceous = m_model->get_mean_temperature_from_curve(Ma) -
        m_model->get_mean_temperature_from_curve(0);

    coeff_mmWS = r_air / r_water_vapour;// coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]
    // coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_lv = 9.1069 in [ / ]
    coeff_lv = lv / ( cp_l * t_0 );
    // coefficient for the specific latent Evaporation heat ( Condensation heat ), coeff_ls = 10.3091 in [ / ]
    coeff_ls = ls / ( cp_l * t_0 );

    c43 = 4./3.;
    c13 = 1./3.;

    pi180 = 180./M_PI;

    // fall velocity for water droplets of 0.1 mm compares to 0.012 m/s
    // fall velocity for water droplets of 0.01 mm compares to 0.8 m/s
    // fall velocity for water droplets of 0.3 mm compares to 2.4 m/s
    // fall velocity for water droplets of 1.0 mm compares to 6.3 m/s
    // fall velocity for water droplets of 5.0 mm compares to 14.1 m/s

    dt_rain_dim = 250.;// dt_rain_dim is the time  in 250 s to pass dr = 400 m, 400 m / 250 s = 1.6 m/s fallout velocity

    // fall velocity for snow flakes of 1.0 mm compares to 0.9 m/s
    // fall velocity for snow flakes of 2.0 mm compares to 1.2 m/s
    // fall velocity for snow flakes of 4.0 mm compares to 1.4 m/s

    dt_snow_dim = 417.;// dt_snow_dim is the time  in 417 s to pass dr = 400 m, 400 m / 417 s = .96 m/s fallout velocity

    dt_dim = L_atm / u_0 * dt;// dimensional time step of the system in s == 0.02 s

    dr_dim = dr * L_atm;  // = 0.025 * 16000 = 400

    cout.precision ( 8 );
    cout.setf ( ios::fixed );

    // Array "i_topography" integer field for mapping the topography
    // land surface in a 2D field
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            for ( int i = im-2; i >= 0; i-- ){
                if ( is_land ( h, i, j, k ) ){
                    i_topography[ j ][ k ] = i;
                    break;
                }
            }
        }
    }

    i_half = ( im -1 ) / 2;
    i_max = im - 1;

    d_i_half = ( double ) i_half ;
    d_i_max = ( double ) i_max;

    j_half = ( jm -1 ) / 2;  // position of the sun at 0°S
    j_max = jm - 1;

    d_j_half = ( double ) j_half;
    d_j_max = ( double ) j_max;

    k_half = ( km -1 ) / 2;  // position of the sun at 0° oder 180° ( Greenwich )
    k_max = km - 1;

    d_k_half = ( double ) k_half;
    d_k_max = ( double ) k_max;
}

BC_Thermo::~BC_Thermo(){}

void cAtmosphereModel::BC_Radiation_multi_layer(){
    // class element for the computation of the radiation and the temperature distribution
    // computation of the local temperature based on short and long wave radiation
    // multi layer radiation model
    if(debug){
        Array tmp = (t-1)*t_0;
        logger()<<"20180912: Enter RML ... "<<std::endl;
        tmp.inspect("20180912: ");
    }

    double rad_eff = rad_pole - rad_equator;
    double albedo_co2_eff = albedo_pole - albedo_equator;

    double j_max_half = ( jm -1 ) / 2;
    // effective temperature, albedo and emissivity/absorptivity for the two layer model
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            for ( int i = 0; i < im-1; i++ ){
                if ( is_ocean_surface(h, i, j, k) || is_land_surface(h, i, j, k) ){    
                    albedo.y[ j ][ k ] = albedo_co2_eff * parabola( j / j_max_half ) + albedo_pole;
                }
            }
        }
    }

    // absorption/emissivity computation
    double epsilon_eff_max = .594; // constant  given by Häckel ( F. Baur and H. Philips, 1934 )
    // constant value stands for other non-condensable gases than water vapour in the equation for epsilon
    double epsilon_eff_2D = epsilon_pole - epsilon_equator;

    for ( int j = 0; j < jm; j++ ){
        int i_trop = m_model->get_tropopause_layer(j);

        // on zero level, lateral parabolic distribution
        epsilon_eff_max = epsilon_eff_2D * parabola( j / j_max_half ) + epsilon_pole;

        for ( int k = 0; k < km; k++ ){
            int i_mount = i_topography[ j ][ k ];

            // in W/m², assumption of parabolic surface radiation at zero level
            radiation_surface.y[ j ][ k ] = rad_eff * parabola( j / j_max_half ) + rad_pole;

            for ( int i = 0; i <= i_trop; i++ ){
                if ( c.x[ i ][ j ][ k ] < 0. )      c.x[ i ][ j ][ k ] = 0.;
                if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
                if ( ice.x[ i ][ j ][ k ] < 0. )    ice.x[ i ][ j ][ k ] = 0.;

                // COSMO water vapour pressure based on local water vapour, cloud water, cloud ice in hPa
                double e = ( c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k ] ) * p_stat.x[ i ][ j ][ k ] / ep;
                
                // radial parabolic distribution, start on zero level
                double epsilon_eff = epsilon_eff_max - ( epsilon_tropopause - epsilon_eff_max ) *
                    parabola( m_model->get_layer_height(i) / m_model->get_layer_height(i_trop) );
                
                double co2_coeff = 1.;
                if ( fabs(m_model->CO2 - 1) < std::numeric_limits<double>::epsilon() ){
                    // influence of co2 in the atmosphere, co2_coeff = 1. means no influence
                    co2_coeff = co2_factor * ( co2_equator / co2_tropopause );
                }

                // dependency given by Häckel ( F. Baur and H. Philips, 1934 )
                if( i >= i_mount ){ //start from the mountain top
                    epsilon_3D.x[ i ][ j ][ k ] = co2_coeff * epsilon_eff + .0416 * sqrt ( e );
                    radiation_3D.x[ i ][ j ][ k ] = ( 1. - epsilon_3D.x[ i ][ j ][ k ] ) * sigma * 
                                    pow ( t.x[ i ][ j ][ k ] * t_0, 4. );
                }
                if ( epsilon_3D.x[ i ][ j ][ k ] > 1. )  epsilon_3D.x[ i ][ j ][ k ] = 1.;
            }
            epsilon.y[ j ][ k ] = epsilon_3D.x[ i_trop ][ j ][ k ];

            // inside mountains
            for ( int i = i_mount - 1; i >= 0; i-- ){
                epsilon_3D.x[ i ][ j ][ k ] = epsilon_3D.x[ i_mount ][ j ][ k ];
                radiation_3D.x[ i ][ j ][ k ] = radiation_3D.x[ i_mount ][ j ][ k ];
            }

            //above tropopause
            for ( int i = i_trop; i < im; i++ ){
                epsilon_3D.x[ i ][ j ][ k ] = epsilon_3D.x[ i_trop ][ j ][ k ];
                t.x[ i ][ j ][ k ] = t_tropopause;
                radiation_3D.x[ i ][ j ][ k ] = ( 1. - epsilon_3D.x[ i ][ j ][ k ] ) * sigma *
                                 pow ( t.x[ i ][ j ][ k ] * t_0, 4. );
            }
        }
    }

    // iteration procedure for the computation of the temperature based on the multi-layer radiation model
    // temperature needs an initial guess which must be corrected by the long wave radiation remaining in the atmosphere

    for ( int iter_rad = 1;  iter_rad <= 4; iter_rad++ ){ // iter_rad may be varied
        // coefficient formed for the tridiogonal set of equations for the absorption/emission coefficient of the multi-layer radiation model
        for ( int j = 0; j < jm; j++ ){
            int i_trop = get_tropopause_layer(j);

            for ( int k = 0; k < km; k++ ){
                int i_mount = i_topography[ j ][ k ];

                std::vector<double> alfa(im, 0);
                std::vector<double> beta(im, 0);
                std::vector<double> AA(im, 0);
                std::vector<std::vector<double> > CC(im, std::vector<double>(im, 0));
                double CCC = 0, DDD = 0;

                // radiation leaving the atmosphere above the tropopause, later needed for non-dimensionalisation
                radiation_3D.x[ i_trop ][ j ][ k ] = ( 1. - epsilon_3D.x[ i_trop ][ j ][ k ] ) * sigma * 
                    pow ( t.x[ i_trop ][ j ][ k ] * t_0, 4. ); 

                // back radiation absorbed by the first water vapour layer out of 40
                double radiation_back = epsilon_3D.x[ i_mount + 1 ][ j ][ k ] * sigma * 
                    pow ( t.x[ i_mount + 1 ][ j ][ k ] * t_0, 4. );
                double atmospheric_window = .1007 * radiation_surface.y[ j ][ k ]; // radiation loss through the atmospheric window
                double rad_surf_diff = radiation_back + radiation_surface.y[ j ][ k ] - atmospheric_window; // radiation leaving the surface

                double fac_rad = ( double ) i_mount * .07 + 1.;  // linear increase with hight, best choice for Ma>0
                // compensation of the missing water vapour at the place of mountain areas to result in a higher emissivity 
                //for higher back radiation
                rad_surf_diff = fac_rad * rad_surf_diff;

                AA[ i_mount ] = rad_surf_diff / radiation_3D.x[ i_trop ][ j ][ k ];// non-dimensional surface radiation
                CC[ i_mount ][ i_mount ] = 0.; // no absorption of radiation on the surface by water vapour

                radiation_3D.x[ i_mount ][ j ][ k ] = ( 1. - epsilon_3D.x[ i_mount ][ j ][ k ] ) * sigma * 
                    pow ( t.x[ i_mount ][ j ][ k ] * t_0, 4. ) / radiation_3D.x[ i_trop ][ j ][ k ]; // radiation leaving the surface

                for ( int i = i_mount + 1; i <= i_trop; i++ ){
                    AA[ i ] = AA[ i - 1 ] * ( 1. - epsilon_3D.x[ i ][ j ][ k ] ); // transmitted radiation from each layer
                    double tmp = sigma * pow ( t.x[ i ][ j ][ k ] * t_0, 4. ) / radiation_3D.x[ i_trop ][ j ][ k ];
                    CC[ i ][ i ]= epsilon_3D.x[ i ][ j ][ k ] * tmp; // absorbed radiation in each layer
                    radiation_3D.x[ i ][ j ][ k ] = ( 1. - epsilon_3D.x[ i ][ j ][ k ] ) * tmp; // radiation leaving each layer

                    for ( int l = i_mount + 1; l <= i_trop; l++ ){
                        // additional transmitted radiation from layer to layer in radial direction
                        CC[ i ][ l ] = CC[ i ][ l - 1 ] * ( 1. - epsilon_3D.x[ l ][ j ][ k ] );
                    }
                }

                // Thomas algorithm to solve the tridiogonal equation system for the solution of the radiation with a recurrence formula
                // additionally embedded in an iterational process
                double aa, bb, cc, dd;
                for ( int i = i_mount; i < i_trop; i++ ){ // values at the surface
                    if ( i == i_mount ){
                        bb = - radiation_3D.x[ i ][ j ][ k ];
                        cc = radiation_3D.x[ i + 1 ][ j ][ k ];
                        dd = - AA[ i ];
                        alfa[ i ] = cc / bb;
                        beta[ i ] = dd / bb;
                    }else{
                        for ( int l = i_mount + 1; l <= i - 1; l++ ){
                            CCC = CCC + CC[ l ][ i ];
                        }
                        for ( int l = i_mount + 1; l <= i - 2; l++ ){
                            DDD = DDD + CC[ l ][ i - 1 ];
                        }
                        aa = radiation_3D.x[ i - 1 ][ j ][ k ];
                        bb = - 2. * radiation_3D.x[ i ][ j ][ k ];
                        cc = radiation_3D.x[ i + 1 ][ j ][ k ];
                        dd = - AA[ i - 1 ] + AA[ i ] + CCC - DDD;
                        alfa[ i ] = cc / ( bb - aa * alfa[ i - 1 ] );
                        beta[ i ] = ( dd - aa * beta[ i - 1 ] ) / ( bb - aa * alfa[ i - 1 ] );
                    }
                }

                // radiation leaving the atmosphere above the tropopause, later needed for non-dimensionalisation
                t.x[ i_trop ][ j ][ k ] = t_tropopause;
                radiation_3D.x[ i_trop ][ j ][ k ] = ( 1. - epsilon_3D.x[ i_trop ][ j ][ k ] ) * sigma * 
                        pow ( t.x[ i_trop ][ j ][ k ] * t_0, 4. );

                // recurrence formula for the radiation and temperature
                for ( int i = i_trop - 1; i >= i_mount; i-- ){
                    // above assumed tropopause constant temperature t_tropopause
                    // Thomas algorithm, recurrence formula
                    radiation_3D.x[ i ][ j ][ k ] = - alfa[ i ] * radiation_3D.x[ i + 1 ][ j ][ k ] + beta[ i ];
                    t.x[ i ][ j ][ k ] = .5 * ( t.x[ i ][ j ][ k ] + pow ( radiation_3D.x[ i ][ j ][ k ] / sigma, 
                        ( 1. / 4. ) ) / t_0 );    // averaging of temperature values to smooth the iterations
                }

                for ( int i = i_trop; i < im; i++ ){ // above tropopause
                    t.x[ i ][ j ][ k ] = t_tropopause;
                    radiation_3D.x[ i ][ j ][ k ] = ( 1. - epsilon_3D.x[ i_trop ][ j ][ k ] ) * sigma * 
                        pow ( t.x[ i ][ j ][ k ] * t_0, 4. );
                }
            }
        }
    }
    logger() << "exit BC_Radiation_multi_layer: temperature max: " << (t.max() - 1)*t_0 << std::endl << std::endl;

    if(debug){
        Array tmp = (t-1)*t_0;
        logger()<<"20180912: Exit RML ... "<<std::endl;
        tmp.inspect("20180912: ");
    }
    
}


void cAtmosphereModel::init_temperature()
{
    if(debug){
        Array tmp = (t-1)*t_0;
        logger()<<"20180912: Enter BCT ... "<<std::endl;
        tmp.inspect("20180912: ");
    }
    // boundary condition of  temperature on land 
    // parabolic distribution from pole to pole accepted
    // temperature on land at equator t_max = 1.055 compares to 15° C compared to 288 K
    // temperature at tropopause t_min = 0.77 compares to -62° C compares to 211 K
    // temperature at tropopause t_min = 0.89 compares to -30° C compares to 243 K
    // temperature difference from equator to pole   18°C compares to  t_delta = 0.0659  compares to  18 K

    // logger() << std::endl << "enter BC_Temperature: temperature max: " << (t.max()-1)*t_0 << std::endl;
    // logger() << "enter BC_Temperature: temperature min: " << (t.min()-1)*t_0 << std::endl << std::endl;

    // Lenton_etal_COPSE_time_temp, constant cretaceous mean temperature, added to the surface initial temperature
    // difference between mean temperature ( Ma ) and mean temperature ( previous Ma ) == t_cretaceous_add
    double t_cretaceous_add = 0; 
    if(!is_first_time_slice()){
        t_cretaceous_add = get_mean_temperature_from_curve(*get_current_time()) - 
            get_mean_temperature_from_curve(*get_previous_time());
        t_cretaceous_add /= t_0; 
    }

    // temperatur distribution at aa prescribed sun position
    // sun_position_lat = 60,    position of sun j = 120 means 30°S, j = 60 means 30°N
    // sun_position_lon = 180, position of sun k = 180 means 0° or 180° E ( Greenwich, zero meridian )
    // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature ( linear equation + parabola )

    if ( ( *get_current_time() > 0 ) && ( sun == 1 ) ){
        double j_par = sun_position_lat; // position of maximum temperature, sun position
        j_par = j_par + declination; // angle of sun axis, declination = 23,4°
        double j_pol = jm - 1;
        double j_par_f = ( double ) j_par;
        double j_pol_f = ( double ) j_pol;

        double aa = ( t_equator - t_pole ) / ( ( ( j_par_f * j_par_f ) - ( j_pol_f * j_pol_f ) ) - 2. * j_par_f * 
            ( j_par_f - j_pol_f ) );
        double bb = - 2. * aa * j_par_f;
        double cc = t_equator + aa * j_par_f * j_par_f;
        double j_d = sqrt ( ( cc - t_pole ) / aa );
        double dd = 2. * aa * j_d + bb;
        double e = t_pole;

        // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature ( linear equation + parabola )
        for ( int k = 0; k < km; k++ ){
            for ( int j = 0; j < jm; j++ ){
                double d_j = ( double ) j;
                if ( d_j <= j_d ){
                    t.x[ 0 ][ j ][ k ] = dd * d_j + e + t_cretaceous_add;
                }
                if ( d_j > j_d ){
                    t.x[ 0 ][ j ][ k ] = aa * d_j * d_j + bb * d_j + cc + t_cretaceous_add;
                }
            }
        }

        // longitudinally variable temperature distribution from west to east in parabolic form
        // communicates the impression of local sun radiation on the southern hemisphere
        double k_par = sun_position_lon;  // position of the sun at constant longitude
        double k_pol = km - 1;

        double t_360 = (  t_0 + 5. ) / t_0;

        for ( int j = 0; j < jm; j++ ){
            double jm_temp_asym = t.x[ 0 ][ j ][ 20 ];//transfer of zonal constant temperature into aa 1D-temperature field
            for ( int k = 0; k < km; k++ ){
                double k_par_f = ( double ) k_par;
                double k_pol_f = ( double ) k_pol;
                double d_k = ( double ) k;

                aa = ( jm_temp_asym - t_360 ) / ( ( ( k_par_f * k_par_f ) - ( k_pol_f * k_pol_f ) ) - 2. * k_par_f * 
                    ( k_par_f - k_pol_f ) );
                bb = - 2. * aa * k_par_f;
                cc = jm_temp_asym + aa * k_par_f * k_par_f;

                t.x[ 0 ][ j ][ k ] = aa * d_k * d_k + bb * d_k + cc;
            }
        }
    }// temperatur distribution at aa prescribed sun position


    // pole temperature adjustment, combination of linear time dependent functions 
    // Stein/Rüdiger/Parish locally constant pole temperature
    // difference between pole temperature ( Ma ) and pole temperature ( previous Ma )
    double t_pole_diff_ocean = 0., t_pole_diff_land;

    std::map<int, double> pole_temp_map{  // Stein/Rüdiger/Parish linear pole temperature ( Ma ) distribution
        {0, 0.},
        {40, 22. },
        {45, 23.5 },
        {50, 24.1 },
        {55, 24.3 },
        {60, 22.4 },
        {70, 24.2 },
        {80, 23.7 },
        {90, 22.8 },
        {100, 21.8},
        {120, 19.},
        {130, 17.8},
        {140, 16.9},
        {150, 16.4},
        {160, 16.},
        {340, 16.}
    }; 

    double d_j_half = ( double ) ( jm -1 ) / 2.0;
    int i_mount = 0;

    // temperature initial conditions along the surface

    if ( RadiationModel == 1 ){

//    logger() << "RadiationModel BC_Temperature: temperature max: " << (t.max()-1)*t_0 << std::endl;
//    logger() << "   RadiationModel = " << RadiationModel << endl;
//    logger() << "   temperature in RadiationModel == 1 " << endl << "RadiationModel BC_Temperature: temperature = " <<  t.x[ 0 ][ 0 ][ 0 ] << endl;

        //the t_pole_diff_ocean should be the difference between this time slice and the previous one, right? -- mchin
        if(!is_first_time_slice()){
            t_pole_diff_ocean = get_pole_temperature(*get_current_time(), pole_temp_map) - 
                get_pole_temperature(*get_previous_time(), pole_temp_map);
        }
        //on land, the difference is between this time slice and present day because 
        //the temperature data from previous time on land is not used.
        t_pole_diff_land = get_pole_temperature(*get_current_time(), pole_temp_map );
        // in °C, constant local pole temperature as function of Ma for hothouse climates 

        double t_eff = t_pole - t_equator;  // coefficient for the zonal parabolic temperature distribution
//        t_eff = ( t_pole_ma + t_0 ) / t_0 - t_equator;  // coefficient for the zonal parabolic temperature distribution

        for ( int k = 0; k < km; k++ ){
            for ( int j = 0; j < jm; j++ ){
                double d_j = ( double ) j;
                if ( NASATemperature == 0 ){  // parabolic ocean surface temperature assumed
                    t.x[ i_mount ][ j ][ k ] = t_eff * parabola( d_j / d_j_half ) + t_pole + t_cretaceous_add;
                                                                      // increasing pole and mean temperature ( Ma ) incorporated

                    if ( is_land ( h, 0, j, k ) ){  // parabolic land surface temperature assumed
                        t.x[ i_mount ][ j ][ k ] = t_eff * parabola( d_j / d_j_half ) + t_pole
                            + t_cretaceous_add + m_model->t_land;
                                                                      // increasing pole and mean temperature ( Ma ) incorporated
                                                                      // in case land temperature is assumed to be
                                                                      // globally higher than ocean temperature, t_land is added too
                    }
                }else{  // if ( NASATemperature == 1 ) ocean surface temperature based on NASA temperature distribution
                    // transported for later time slices Ma by use_earthbyte_reconstruction
                    if ( is_land (h, 0, j, k ) ){  // on land a parabolic distribution assumed, no NASA based data transportable
                        t.x[ 0 ][ j ][ k ] = t_eff * parabola( d_j / d_j_half ) + t_pole
//                        t.x[ i_mount ][ j ][ k ] = t_eff * parabola( d_j / d_j_half ) + ( t_pole_ma + t_0 ) / t_0
                            + t_cretaceous_add + m_model->t_land;

                            // Stein/Rüdiger/Parish pole temperature decreasing equator wards
                            t.x[ 0 ][ j ][ k ] += t_pole_diff_land * fabs ( parabola( d_j / d_j_half ) + 1. ) / t_0;
                    }else{ // if the this location is ocean
                        if(*get_current_time() > 0){//when the current time is 0, 
                                                    //the temperature data has already been read into t.x[ 0 ][ j ][ k ]        
                            // ocean surface temperature increased by mean t_cretaceous_add
                            // and by a zonally equator wards decreasing temperature difference is added
                            // Stein/Rüdiger/Parish pole temperature decreasing equator wards
                            t.x[ 0 ][ j ][ k ] += t_cretaceous_add + 
                                t_pole_diff_ocean * fabs ( parabola( d_j / d_j_half ) + 1. ) / t_0;
                        }
                    }
                }// else ( NASATemperature == 1 )
            }// for j
        }// for k
    }// if ( RadiationModel == 1 )

    // zonal temperature along tropopause
    double t_eff_tropo = t_tropopause_pole - t_tropopause;

    //use "linear temperature decay" to generate temperature data for layers between mountain top and tropopause
    //use "mountain top temperature" for the layers below mountain top
    //use "tropopause tempeature" for the layers above tropopause
    // temperature approaching the tropopause, above constant temperature following Standard Atmosphere
    for ( int j = 0; j < jm; j++ ){
        double temp_tropopause =  t_eff_tropo * parabola( j / ((jm-1)/2.0) ) +
                t_tropopause_pole + t_cretaceous_add;   //temperature at tropopause     

        for ( int k = 0; k < km; k++ ){
            int i_mount = i_topography[ j ][ k ];
            int i_trop = get_tropopause_layer(j);

            double t_mount_top = ( temp_tropopause - t.x[ 0 ][ j ][ k ] ) *
                (get_layer_height(i_mount) / get_layer_height(i_trop)) + 
                t.x[ 0 ][ j ][ k ]; //temperature at mountain top

            for ( int i = 1; i < im; i++ ){
                if ( i < i_trop+1 ){
                    if(i>i_mount){
                        // linear temperature decay up to tropopause, privat  approximation
                        t.x[ i ][ j ][ k ] = ( temp_tropopause - t.x[ 0 ][ j ][ k ] ) * 
                            (get_layer_height(i) / get_layer_height(i_trop)) + 
                            t.x[ 0 ][ j ][ k ]; 
                    }else{
                        t.x[ i ][ j ][ k ] = t_mount_top; //inside mountain
                    }
                }else{ // above tropopause
                    t.x[ i ][ j ][ k ] = temp_tropopause;
                }
            }
            t.x[ 0 ][ j ][ k ] = t_mount_top;
        }
    }

    logger() << "exit BC_Temperature: temperature max: " << (t.max()-1)*t_0 << std::endl << std::endl;
    if(debug){
        Array tmp = (t-1)*t_0;
        logger()<<"20180912: Exit BCT ... "<<std::endl;
        tmp.inspect("20180912: ");
    }
}


void BC_Thermo::BC_CO2( Array_2D &Vegetation, Array &h, Array &t, Array &p_dyn, Array &co2 ){
    // initial and boundary conditions of CO2 content on water and land surfaces
    // parabolic CO2 content distribution from pole to pole accepted

    j_half = j_max / 2;

    // temperature-distribution by Ruddiman approximated by a parabola
    //t_cretaceous_eff = t_cretaceous_max / ( ( double ) Ma_max_half - ( double ) ( Ma_max_half * Ma_max_half / ( double ) Ma_max ) );   // in °C
    //t_cretaceous = t_cretaceous_eff * ( double ) ( - ( Ma * Ma ) / ( double ) Ma_max + Ma );   // in °C
    //if ( Ma == 0 )  t_cretaceous = t_cretaceous_prev = 0.;

    // CO2-distribution by Ruddiman approximated by a parabola
    co2_cretaceous = 3.2886 * pow ( ( t_cretaceous + t_average ), 2 ) - 32.8859 *
        ( t_cretaceous + t_average ) + 102.2148;  // in ppm
    co2_average = 3.2886 * pow ( t_average, 2 ) - 32.8859 * t_average + 102.2148;  // in ppm
    co2_cretaceous = co2_cretaceous - co2_average;

    cout.precision ( 3 );

    string co_comment = "      co2 increase at cretaceous times: ";
    string co_gain = " co2 increase";
    string co_modern = "      mean co2 at modern times: ";
    string co_cretaceous_str = "      mean co2 at cretaceous times: ";
    string co_average_str = " co2 modern";
    string co_average_cret = " co2 cretaceous";
    string co_unit =  "ppm ";

    cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) <<
        co_comment << resetiosflags ( ios::left )         << setw ( 12 ) << co_gain << " = "
        << setw ( 7 ) << setfill ( ' ' ) << co2_cretaceous << setw ( 5 ) << co_unit << 
        endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) << co_modern
        << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_str  << " = "
        << setw ( 7 )  << setfill ( ' ' ) << co2_average << setw ( 5 ) << co_unit 
        << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left )
        << co_cretaceous_str << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_cret
        << " = "  << setw ( 7 )  << setfill ( ' ' ) << co2_average + co2_cretaceous
        << setw ( 5 ) << co_unit << endl;
    cout << endl;

    d_i_max = ( double ) i_max;
    d_j_half = ( double ) j_half;

    double co2_0 = m_model->co2_0;
    co2_equator = co2_equator / co2_0;
    co2_pole = co2_pole / co2_0;
    co2_cretaceous = co2_cretaceous / co2_0;
    co2_land = co2_land / co2_0;
    co2_ocean = co2_ocean / co2_0;
    co2_vegetation = co2_vegetation / co2_0;
    co2_tropopause = co2_tropopause / co2_0;

    co2_eff = co2_pole - co2_equator;

    // CO2-content as initial solution
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            int i_mount = i_topography[ j ][ k ];

            if ( is_air ( h, i_mount, j, k ) ){
                d_j = ( double ) j;
                co2.x[ i_mount ][ j ][ k ] = co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + 
                    co2_pole + co2_cretaceous + co2_ocean; // non-dimensional
            }
            if ( is_land ( h, i_mount, j, k ) ){
                d_j = ( double ) j;
                co2.x[ i_mount ][ j ][ k ] = co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + 
                    co2_pole + co2_cretaceous + co2_land - co2_vegetation * Vegetation.y[ j ][ k ];  // parabolic distribution from pole to pole
            }
        }
    }

    // co2 distribution decreasing approaching tropopause, above no co2
    for ( int j = 0; j < jm; j++ ){
        // i_trop = im_tropopause[ j ];
        int i_trop = im_tropopause[ j ] + GetTropopauseHightAdd ( t_cretaceous / t_0 );
        d_i_max = ( double ) i_trop;
        for ( int k = 0; k < km; k++ ){
            int i_mount = i_topography[ j ][ k ];
            for ( int i = 1; i <= im - 1; i++ ){
                if ( i <= i_trop ){
                    d_i = ( double ) i;
                    co2.x[ i ][ j ][ k ] = co2.x[ i_mount ][ j ][ k ] - ( co2_tropopause - co2.x[ i_mount ][ j ][ k ] ) * 
                        ( d_i / d_i_max * ( d_i / d_i_max - 2. ) );// radial distribution approximated by a parabola
                }
                //else co2.x[ i ][ j ][ k ] = co2.x[ i_trop ][ j ][ k ];
                else  co2.x[ i ][ j ][ k ] = co2_tropopause;
            }
            for ( int i = i_trop - 1; i >= 0; i-- ){
                if ( ( is_land ( h, i, j, k ) ) != ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i+1, j, k ) ) ) ){
                    co2.x[ i ][ j ][ k ] = co2.x[ i_mount ][ j ][ k ];
                }
            }
        }
    }
}

void BC_Thermo::TropopauseLocation(){
// parabolic tropopause location distribution from pole to pole assumed
    j_max = jm - 1;
    j_half = ( jm -1 ) / 2;
//  j_infl = 45;                                                        // fixed value due to the cubic function for the location of the tropopause
    j_infl = 43;
    // flattens the equator peak
    d_j_half = ( double ) j_half;
    d_j_infl = 3. * ( double ) j_infl;
    trop_co2_eff = ( double ) ( tropopause_pole - tropopause_equator );

// no stripes in longitudinal direction
// computation of the tropopause from pole to pole, constant approach
    for ( int j = 0; j < jm; j++ ){ // parabolic approach
        im_tropopause[ j ] = tropopause_equator; // constant approach
    }

/*
// minor stripes in longitudinal direction
// computation of the tropopause from pole to pole, parabolic approach
    for ( int j = 0; j < jm; j++ ) // parabolic approach
    {
        d_j = ( double ) j;
        im_tropopause[ j ] = ( trop_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) -
            2. * d_j / d_j_half ) ) + tropopause_pole; // parabolic approach
    }
*/

/*
// visibal stripes in longitudinal direction
    double trop = 0;  // cubic approach

// computation of the tropopause from pole to pole, cubic approach
    for ( int j = 0; j <= j_half; j++ )  // cubic approach
    {
        d_j = ( double ) j;
        trop = ( - trop_co2_eff * ( d_j * d_j * d_j - d_j_infl *d_j * d_j ) /
            ( d_j_half * d_j_half * d_j_half - d_j_infl * d_j_half * d_j_half ) +
            ( double ) tropopause_pole );  // cubic approach

        im_tropopause[ j ] = ( int ) trop;  // cubic approach
    }

    for ( int j = j_half + 1; j < jm; j++ )  // cubic approach
    {
        im_tropopause[ j ] = im_tropopause[ j_max - j ];  // cubic approach
    }
*/
}

void BC_Thermo::smooth_transition(Array &u, Array &v, Array &w, int lat){
    int start = lat-3, end = lat+3;
    for ( int i = 0; i < im; i++ ){
        for ( int k = 0; k < km; k++ ){
            for ( int j = start; j <= end; j++ ){
                u.x[ i ][ j ][ k ] = ( u.x[ i ][ end ][ k ] - u.x[ i ][ start ][ k ] ) *
                    ( ( double ) ( j - start )  / ( end - start ) ) + u.x[ i ][ start ][ k ];
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ end ][ k ] - v.x[ i ][ start ][ k ] ) *
                    ( ( double ) ( j - start )  / ( end - start ) ) + v.x[ i ][ start ][ k ];
                w.x[ i ][ j ][ k ] = ( w.x[ i ][ end ][ k ] - w.x[ i ][ start ][ k ] ) *
                    ( ( double ) ( j - start )  / ( end - start ) ) + w.x[ i ][ start ][ k ];
            }
        }
    }
}


void BC_Thermo::form_diagonals(Array &a, int start, int end){
    for ( int k = 0; k < km; k++ ){
        for ( int j = start; j < end; j++ ){
            for ( int i = 0; i < im; i++ ){
                a.x[ i ][ j ][ k ] = ( a.x[ i ][ end ][ k ] - a.x[ i ][ start ][ k ] ) *
                    ( j - start ) / (double)(end - start) + a.x[ i ][ start ][ k ];
            }
        }
    }
}

void  BC_Thermo::init_u(Array &u, int lat_1, int lat_2, double coefficient){
    for(int j = lat_1; j < lat_2+1; j++ ){
        int tropopause_layer = m_model->get_tropopause_layer(j);             
        double tropopause_height = m_model->get_layer_height(tropopause_layer);
        for( int k = 0; k < km; k++ ){
            for( int i = 0; i < tropopause_layer; i++ ){
                u.x[ i ][ j ][ k ] = -coefficient * 
                    parabola_interp(-1, 0, m_model->get_layer_height(i)*2/tropopause_height);
            }
        }
    }
}

void  BC_Thermo::init_v_or_w(Array &v_or_w, int lat_1, int lat_2, double coeff_trop, double coeff_sl){
    for(int j = lat_1; j < lat_2+1; j++ ){
        int tropopause_layer = m_model->get_tropopause_layer(j);
        double tropopause_height = m_model->get_layer_height(tropopause_layer);
        for( int k = 0; k < km; k++ ){
            for( int i = 0; i < tropopause_layer; i++ ){
                v_or_w.x[ i ][ j ][ k ] = ( coeff_trop - coeff_sl ) *
                    m_model->get_layer_height(i)/tropopause_height + coeff_sl;
            }
        }
    }
    init_v_or_w_above_tropopause(v_or_w, lat_1, lat_2, coeff_trop);
}

void  BC_Thermo::init_v_or_w_above_tropopause(Array &v_or_w, int lat_1, int lat_2, double coeff){
    for(int j = lat_1; j < lat_2+1; j++ ){
        int tropopause_layer = m_model->get_tropopause_layer(j);
        if(tropopause_layer >= im-1) return;
        double tropopause_height = m_model->get_layer_height(tropopause_layer);
        for( int k = 0; k < km; k++ ){
            for( int i = tropopause_layer; i < im; i++ ){
                v_or_w.x[ i ][ j ][ k ] = coeff * (m_model->get_layer_height(im-1) - m_model->get_layer_height(i)) / 
                    (m_model->get_layer_height(im-1) - tropopause_height);
            }
        }
    }
}


void BC_Thermo::IC_CellStructure ( Array &h, Array &u, Array &v, Array &w ){
// boundary condition for the velocity components in the circulation cells

// latest version by Grotjahn ( Global Atmospheric Circulations, 1993 )
// default for the velocity components u, v, and w as initial conditions

// velocities given in m/s, 1 m/s compares to 3.6 km/h, non-dimensionalized by u_0 at the end of this class element
// do not change the velocity initial conditions !!

// preparations for diagonal velocity value connections
    j_aeq = 90;
    j_pol_n = 0;
    j_pol_s = jm-1;
    j_pol_v_n = 15;
    j_pol_v_s = 165;
    j_fer_n = 30;
    j_fer_s = 150;
    j_fer_v_n = 45;
    j_fer_v_s = 135;
    j_had_n = 60;
    j_had_s = 120;
    j_had_v_n = 75;
    j_had_v_s = 105;

    // initial velocity components in the northern and southern
    // Pole, Ferrel and Hadley cells

    // equator ( at j=90 compares to 0° latitude )
    // u-component up to tropopause and back on half distance (  i = 20 )

    double ua_00 = 1.;  // in m/s compares to 3.6 km/h, non-dimensionalized by u_0 at the end of this function
    double va_equator_SL =  0.000;
    double va_equator_Tropopause = 0.000;
    double wa_equator_SL = - 1.;
    double wa_equator_Tropopause = - 7.5;

    //equator
    init_u(u,90,90,ua_00);//lat: 0
    
    init_v_or_w(v,90,90,va_equator_Tropopause,va_equator_SL);
    init_v_or_w(w,90,90,wa_equator_Tropopause,wa_equator_SL);
    
    //polar cell
    double ua_90 = - 0.5;
    double va_Polar_SL = 0.;
    double va_Polar_Tropopause = 0.;
    double va_Polar_SL_75 = .5;
    double va_Polar_Tropopause_75 = - 1.;
    double wa_Polar_SL = - 0.01;
    double wa_Polar_Tropopause = 0.;

    //northern polar cell
    init_u(u, 0, 0, ua_90); //lat: 90

    init_v_or_w(v,0,30,va_Polar_Tropopause,va_Polar_SL); //lat: 90-60
    init_v_or_w(w,0,30,wa_Polar_Tropopause,wa_Polar_SL); //lat: 90-60
    init_v_or_w(v,15,15,va_Polar_Tropopause_75,va_Polar_SL_75); //lat: 75

    //southern polar cell
    init_u(u, 180, 180, ua_90); //lat: -90

    init_v_or_w(v,150,180,va_Polar_Tropopause,va_Polar_SL); //lat: -90-(-60)
    init_v_or_w(w,150,180,wa_Polar_Tropopause,wa_Polar_SL); //lat: -90-(-60)
    init_v_or_w(v,165,165,va_Polar_Tropopause_75,va_Polar_SL_75); //lat: -75

    //Ferrel cell
    double ua_60 = 0.5;
    double va_Ferrel_SL = 0.5;
    double va_Ferrel_Tropopause = 1.;
    double va_Ferrel_SL_45 = - 0.1;
    double va_Ferrel_Tropopause_45 = 1.;
    double wa_Ferrel_SL = -0.2;    // subpolar jet
    double wa_Ferrel_Tropopause = 10.;           

    //northern Ferrel cell
    init_u(u, 30, 30, ua_60); //lat: 60 

    init_v_or_w(v,30,30,va_Ferrel_Tropopause,va_Ferrel_SL); //lat: 60
    init_v_or_w(w,30,30,wa_Ferrel_Tropopause,wa_Ferrel_SL); //lat: 60
    init_v_or_w(v,45,45,va_Ferrel_Tropopause_45,va_Ferrel_SL_45); //lat: 45   

    //southern Ferrel cell
    init_u(u, 150, 150, ua_60); //lat: -60 

    init_v_or_w(v,150,150,va_Ferrel_Tropopause,va_Ferrel_SL); //lat: -60
    init_v_or_w(w,150,150,wa_Ferrel_Tropopause,wa_Ferrel_SL); //lat: -60
    init_v_or_w(v,135,135,va_Ferrel_Tropopause_45,va_Ferrel_SL_45); //lat: -45   
 
    // Hadley cell
    double ua_30 = - 1.;
    double va_Hadley_SL = .25;
    double va_Hadley_Tropopause = - 1.;
    double va_Hadley_SL_15 = 1.;
    double va_Hadley_Tropopause_15 = - 1.;
    double wa_Hadley_SL = 1.;            // at surface
    double wa_Hadley_Tropopause = 30.;  // subtropic jet in m/s compares to 108 km/h

    //northern Hadley cell
    init_u(u, 60, 60, ua_30); //lat: 30 

    init_v_or_w(v,60,60,va_Hadley_Tropopause,va_Hadley_SL); //lat: 30
    init_v_or_w(w,60,60,wa_Hadley_Tropopause,wa_Hadley_SL); //lat: 30
    init_v_or_w(v,75,75,va_Hadley_Tropopause_15,va_Hadley_SL_15); //lat: 15   

    //southern Hadley cell
    init_u(u, 120, 120, ua_30); //lat: -30 

    init_v_or_w(v,120,120,va_Hadley_Tropopause,va_Hadley_SL); //lat: -30
    init_v_or_w(w,120,120,wa_Hadley_Tropopause,wa_Hadley_SL); //lat: -30
    init_v_or_w(v,105,105,va_Hadley_Tropopause_15,va_Hadley_SL_15); //lat: -15 

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////
/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

////////////// meridional values of w-velocity component from Pol till Ferrel, from Ferrel till Hadley, from Hadley till equator ///////////////

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// w-component formed by the diagonal starting from subtropical jet and North Pole
    //northen hemisphere
    form_diagonals(u, 0, 30);
    form_diagonals(w, 0, 30);
    form_diagonals(v, 0, 15);
    form_diagonals(v, 15, 30);

    form_diagonals(u, 30, 60);
    form_diagonals(w, 30, 60);
    form_diagonals(v, 30, 45);
    form_diagonals(v, 45, 60);

    form_diagonals(u, 60, 90);
    form_diagonals(w, 60, 90);
    form_diagonals(v, 60, 75);
    form_diagonals(v, 75, 90);

    //southen hemisphere
    form_diagonals(u, 90, 120);
    form_diagonals(w, 90, 120);
    form_diagonals(v, 90, 105);
    form_diagonals(v, 105, 120);

    form_diagonals(u, 120, 150);
    form_diagonals(w, 120, 150);
    form_diagonals(v, 120, 135);
    form_diagonals(v, 135, 150);

    form_diagonals(u, 150, 180);
    form_diagonals(w, 150, 180);
    form_diagonals(v, 150, 165);
    form_diagonals(v, 165, 180);

    for ( int i = 0; i < im; i++ ){
        for ( int j = 91; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                v.x[ i ][ j ][ k ] = - v.x[ i ][ j ][ k ];
            }
        }
    }


///////////////////////////////////////////////// smoothing transitions from cell to cell //////////////////////////
///////////////////////////////////////////////// smoothing transitions from cell to cell //////////////////////////

///////////////////////////////////////////////// Northern hemisphere ////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Northern hemisphere
///////////////////////////////////////////////// smoothing transitions from Headley to Ferrel to Polar cells //////////////////////////
    smooth_transition(u,v,w,60); 
    smooth_transition(u,v,w,30);    
    smooth_transition(u,v,w,75);
    smooth_transition(u,v,w,45);

    smooth_transition(u,v,w,90); 

    smooth_transition(u,v,w,120);
    smooth_transition(u,v,w,150);
    smooth_transition(u,v,w,105);
    smooth_transition(u,v,w,135);

    // non dimensionalization by u_0
    for ( int i = 0; i < im; i++ ){
        for ( int k = 0; k < km; k++ ){
            for ( int j = 0; j < jm; j++ ){
                u.x[ i ][ j ][ k ] = u.x[ i ][ j ][ k ] / u_0;
                v.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] / u_0;
                w.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] / u_0;
                if ( is_land ( h, i, j, k ) )     u.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] = 0.;
            }
        }
    }

///////////////////////////////////////////////// end of smoothing transitions from cell to cell //////////////////////////
///////////////////////////////////////////////// end of smoothing transitions from cell to cell //////////////////////////
}





void BC_Thermo::BC_Surface_Temperature_NASA ( const string &Name_SurfaceTemperature_File,
                             Array_2D &temperature_NASA, Array &t ){
// initial conditions for the Name_SurfaceTemperature_File at the sea surface

    cout.precision ( 3 );
    cout.setf ( ios::fixed );

    ifstream Name_SurfaceTemperature_File_Read(Name_SurfaceTemperature_File);

    if (!Name_SurfaceTemperature_File_Read.is_open()) {
        cerr << "ERROR: could not open SurfaceTemperature_File file at "
        << Name_SurfaceTemperature_File << "\n";
        abort();
    }

    k_half = ( km -1 ) / 2;                                                             // position at 180°E ( Greenwich )
    int j = 0;
    int k = 0;

    while ( ( k < km ) && !Name_SurfaceTemperature_File_Read.eof() ){
        while ( j < jm ){
            double lat, lon, temperature;
            Name_SurfaceTemperature_File_Read >> lat;
            Name_SurfaceTemperature_File_Read >> lon;
            Name_SurfaceTemperature_File_Read >> temperature;

            t.x[ 0 ][ j ][ k ] = temperature_NASA.y[ j ][ k ] = ( temperature + t_0 ) / t_0;
            j++;
        }
    j = 0;
    k++;
    }

    // correction of surface temperature around 180°E
    for ( int j = 0; j < jm; j++ ){
        t.x[ 0 ][ j ][ k_half ] = ( t.x[ 0 ][ j ][ k_half + 1 ] + t.x[ 0 ][ j ][ k_half - 1 ] ) / 2.;
        temperature_NASA.y[ j ][ k_half ] = ( temperature_NASA.y[ j ][ k_half + 1 ] +
            temperature_NASA.y[ j ][ k_half - 1 ] ) / 2.;
    }
}


void BC_Thermo::BC_Surface_Precipitation_NASA ( const string &Name_SurfacePrecipitation_File,
                             Array_2D &precipitation_NASA ){
    // initial conditions for the Name_SurfacePrecipitation_File at the sea surface
    cout.precision ( 3 );
    cout.setf ( ios::fixed );

    ifstream Name_SurfacePrecipitation_File_Read(Name_SurfacePrecipitation_File);

    if (!Name_SurfacePrecipitation_File_Read.is_open()) {
        cerr << "ERROR: could not open SurfacePrecipitation_File file at "
        << Name_SurfacePrecipitation_File << "\n";
        abort();
    }

    int j = 0;
    int k = 0;

    while ( ( k < km ) && !Name_SurfacePrecipitation_File_Read.eof() ){
        while ( j < jm ){
            double lat, lon, precipitation;
            Name_SurfacePrecipitation_File_Read >> lat;
            Name_SurfacePrecipitation_File_Read >> lon;
            Name_SurfacePrecipitation_File_Read >> precipitation;

            precipitation_NASA.y[ j ][ k ] = precipitation;
            j++;
        }
    j = 0;
    k++;
    }
}




void BC_Thermo::BC_Pressure ( Array &p_stat, Array &p_dyn, Array &t, Array &h ){
    exp_pressure = g / ( 1.e-2 * gam * R_Air );

// boundary condition of surface pressure given by surface temperature through gas equation
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            p_stat.x[ 0 ][ j ][ k ] =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );      // given in hPa
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 1; i < im; i++ ){
                hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
                p_stat.x[ i ][ j ][ k ] = pow ( ( ( t.x[ 0 ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) /
                    ( t.x[ 0 ][ j ][ k ] * t_0 ) ), exp_pressure ) * p_stat.x[ 0 ][ j ][ k ];
                // linear temperature distribution T = T0 - gam * hight
                // current air pressure, step size in 500 m, from politropic formula in hPa
            }
        }
    }
}









void BC_Thermo::Latent_Heat ( Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h,
                             Array &t, Array &tn, Array &u, Array &v, Array &w, Array &p_dyn,
                             Array &p_stat, Array &c, Array &ice, Array &Q_Latent, Array &Q_Sensible,
                             Array &radiation_3D, Array_2D &Q_radiation, Array_2D &Q_latent,
                             Array_2D &Q_sensible, Array_2D &Q_bottom ){
    double Q_Latent_Ice = 0.; 

// collection of coefficients for phase transformation
    coeff_Lv = lv / ( L_atm / ( double ) ( im-1 ) );                                        // coefficient for Q_latent generated by cloud water
    coeff_Ls = ls / ( L_atm / ( double ) ( im-1 ) );                                        // coefficient for Q_latent generated by cloud ice
    coeff_Q = cp_l * r_air * t_0 / ( L_atm / ( double ) ( im-1 ) );             // coefficient for Q_Sensible
    double coeff_lat = .079;
    double coeff_sen = .15;
    c32 = 3. / 2.;
    c42 = 4. / 2.;
    c12 = 1. / 2.;

// 1. and 2. derivatives for 3 spacial directions and and time in Finite Difference Methods ( FDM )
// collection of coefficients
    dr2 = dr * dr;
    dthe2 = dthe * dthe;
    dphi2 = dphi * dphi;
    rm = rad.z[ 0 ];
    for ( int j = 0; j < jm; j++ ){
// collection of coefficients
        sinthe = sin( the.z[ j ] );
        rmsinthe = rm * sinthe;
        for ( int k = 0; k < km; k++ ){
            int i_mount = i_topography[ j ][ k ];
            t_Celsius = t.x[ i_mount ][ j ][ k ] * t_0 - t_0;  // conversion from Kelvin to Celsius
            T = t.x[ i_mount ][ j ][ k ] * t_0;
            p_h = p_stat.x[ i_mount ][ j ][ k ];
            E_Rain = hp * exp_func ( T, 17.2694, 35.86 );  // saturation water vapour pressure for the water phase at t > 0°C in hPa
            E_Ice = hp * exp_func ( T, 21.8746, 7.66 );  // saturation water vapour pressure for the ice phase in hPa
            q_Rain  = ep * E_Rain / ( p_h - E_Rain );  // water vapour amount at saturation with water formation in kg/kg
            q_Ice  = ep * E_Ice / ( p_h - E_Ice );  // water vapour amount at saturation with ice formation in kg/kg
            e = .01 * c.x[ i_mount ][ j ][ k ] * p_stat.x[ i_mount ][ j ][ k ] / ep;  // water vapour pressure in Pa
            a = e / ( R_WaterVapour * t.x[ i_mount ][ j ][ k ] * t_0 );  // absolute humidity in kg/m³
            Q_Latent.x[ i_mount ][ j ][ k ] = - coeff_Lv * a * ( - 3. * c.x[ i_mount ][ j ][ k ] +
                4. * c.x[ i_mount + 1 ][ j ][ k ] - c.x[ i_mount + 2 ][ j ][ k ] ) / ( 2. * dr );
            Q_Latent_Ice = - coeff_Ls * a * ( - 3. * ice.x[ i_mount ][ j ][ k ] +
                4. * ice.x[ i_mount + 1 ][ j ][ k ] - ice.x[ i_mount + 2 ][ j ][ k ] ) / ( 2. * dr );
            Q_Latent.x[ i_mount ][ j ][ k ] = coeff_lat * ( Q_Latent.x[ i_mount ][ j ][ k ] + Q_Latent_Ice );
            Q_Sensible.x[ i_mount ][ j ][ k ] = - coeff_sen * coeff_Q * ( - 3. * t.x[ i_mount ][ j ][ k ] +
                4. * t.x[ i_mount + 1 ][ j ][ k ] - t.x[ i_mount + 2 ][ j ][ k ] ) / ( 2. * dr );   // sensible heat in [W/m2] from energy transport equation
        }
    }

    for ( int j = 0; j < jm; j++ ){
// collection of coefficients + 
        sinthe = sin( the.z[ j ] );
        sinthe2 = sinthe * sinthe;
        costhe = cos( the.z[ j ] );
        rmsinthe = rm * sinthe;
        rm2sinthe = rm2 * sinthe;
        rm2sinthe2 = rm2 * sinthe2;

// water vapour can condensate/evaporate und sublimate/vaporize
// water vapour turns to or developes from water or ice
// latent heat of water vapour

        for ( int k = 0; k < km; k++ ){
            int i_mount = i_topography[ j ][ k ];
            for ( int i = i_mount + 1; i < im-2; i++ ){
// collection of coefficients
                rm = rad.z[ i ];
                rm2 = rm * rm;
                t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;  // conversion from Kelvin to Celsius
                T = t.x[ i ][ j ][ k ] * t_0; // in K
                p_h = p_stat.x[ i ][ j ][ k ];
                E_Rain = hp * exp_func ( T, 17.2694, 35.86 );  // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Ice = hp * exp_func ( T, 21.8746, 7.66 );  // saturation water vapour pressure for the ice phase in hPa
                q_Rain  = ep * E_Rain / ( p_h - E_Rain );  // water vapour amount at saturation with water formation in kg/kg
                q_Ice  = ep * E_Ice / ( p_h - E_Ice );  // water vapour amount at saturation with ice formation in kg/kg
                e = .01 * c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep;  // water vapour pressure in Pa
                a = e / ( r_water_vapour * t.x[ i ][ j ][ k ] * t_0 );  // absolute humidity in kg/m³
                Q_Latent.x[ i ][ j ][ k ] = - coeff_Lv * a * ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
                Q_Latent_Ice = - coeff_Ls * a * ( ice.x[ i+1 ][ j ][ k ] - ice.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );
                Q_Latent.x[ i ][ j ][ k ] = coeff_lat * ( Q_Latent.x[ i ][ j ][ k ] + Q_Latent_Ice );
                Q_Sensible.x[ i ][ j ][ k ] = - coeff_sen * coeff_Q * ( t.x[ i+1 ][ j ][ k ] -
                    t.x[ i-1 ][ j ][ k ] ) / ( 2. * dr );  // sensible heat in [W/m2] from energy transport equation
            }
        }
    }
}






void cAtmosphereModel::Ice_Water_Saturation_Adjustment()
{ 
    if(debug){
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
        assert(!c.has_nan());
        assert(!t.has_nan());
    }
    cout.precision ( 6 );
// Ice_Water_Saturation_Adjustment, distribution of cloud ice and cloud water dependent on water vapour amount and temperature
/*
    logger() << std::endl << std::endl << "enter %%%%%%%%%%%% Ice_Water_Saturation_Adjustment: temperature max: "
    << (t.max() - 1)*t_0 << std::endl << std::endl;

    logger() << "enter Ice_Water_Saturation_Adjustment: water vapour max: "
        << c.max() * 1000. << std::endl;
    logger() << "enter Ice_Water_Saturation_Adjustment: cloud water max: "
        << cloud.max() * 1000. << std::endl;
    logger() << "enter Ice_Water_Saturation_Adjustment: cloud ice max: "
        << ice.max() * 1000. << std::endl << std::endl;
*/
// constant coefficients for the adjustment of cloud water and cloud ice amount vice versa
    float t_00 = 236.15;
    float t_Celsius_2 = t_00 - t_0; // in Kelvin = -37 °C
    float exp_pressure = g / ( 1.e-2 * gam * R_Air );
    float q_v_hyp, q_T;
    float dt_dim = L_atm / u_0 * dt;// dimensional time step of the system in s == 0.02 s
    // setting water vapour, cloud water and cloud ice into the proper thermodynamic ratio based on the local temperatures
    // starting from a guessed parabolic temperature and water vapour distribution in north/south direction
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
/** %%%%%%%%%%%%%%%%%%%%%%%%%%     saturation pressure     %%%%%%%%%%%%%%%%%%%%%%%%%%%% **/
                float t_u = t.x[ i ][ j ][ k ] * t_0; // in K
                float t_Celsius = t_u - t_0; // in C

                float p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                float hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
                float p_h;
                if ( i != 0 )           p_h = pow ( ( ( t.x[ i ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) /
                                                      ( t.x[ i ][ j ][ k ] * t_0 ) ), exp_pressure ) * p_SL;
                else                     p_h = p_SL;

                float E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                float E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 ); // saturation water vapour pressure for the ice phase in hPa
                float q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                float q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg

/** %%%%%%%%%%%%%%%%%%%%%%%%%%%     warm cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% **/

// warm cloud phase in case water vapour is over-saturated
                if ( t_Celsius >= 0. ){
                    q_T = c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ]; // total water content
                    t_u = t.x[ i ][ j ][ k ] * t_0; // in K
                    t_Celsius = t_u - t_0;

                    p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                    hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
                    if ( i != 0 )           p_h = pow ( ( ( t.x[ i ][ j ][ k ] * t_0 - gam * hight * 1.e-2 ) /
                                                          ( t.x[ i ][ j ][ k ] * t_0 ) ), exp_pressure ) * p_SL;
                    else                    p_h = p_SL;

                    E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                    q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                    float q_Rain_n = q_Rain;
                    float T_it;
                    if ( q_T <= q_Rain ){ /**     subsaturated     **/
                        c.x[ i ][ j ][ k ] = q_T; // total water amount as water vapour
                        cloud.x[ i ][ j ][ k ] = 0.; // no cloud water available
                        ice.x[ i ][ j ][ k ] = 0.; // no cloud ice available above 0 °C
/*
if ( ( i == 5 ) && ( j == 90 ) && ( k == 180 ) ){
    logger() << "no cloud               Ice_Water_Saturation_Adjustment: temperature max: " << (t.max() - 1)*t_0 <<"          iter_prec: " << iter_prec << std::endl;
    logger() << "no cloud               Ice_Water_Saturation_Adjustment: water vapour max: " << c.max() * 1000. << std::endl;
    logger() << "no cloud               Ice_Water_Saturation_Adjustment: cloud water max: " << cloud.max() * 1000. << std::endl;
    logger() << "no cloud               Ice_Water_Saturation_Adjustment: cloud ice max: " << ice.max() * 1000. << std::endl << std::endl;
}
*/
                    }else{ /**     oversaturated     **/
                        for(int iter_prec = 1; iter_prec <= 20; iter_prec++ ){ // iter_prec may be varied
/*
if ( ( i == 5 ) && ( j == 90 ) && ( k == 180 ) ){
    logger() << "warm cloud               Ice_Water_Saturation_Adjustment: temperature max: " << (t.max() - 1)*t_0 <<"          iter_prec: " << iter_prec << std::endl;
    logger() << "warm cloud               Ice_Water_Saturation_Adjustment: water vapour max: " << c.max() * 1000. << std::endl;
    logger() << "warm cloud               Ice_Water_Saturation_Adjustment: cloud water max: " << cloud.max() * 1000. << std::endl;
    logger() << "warm cloud               Ice_Water_Saturation_Adjustment: cloud ice max: " << ice.max() * 1000. << std::endl << std::endl;
}
*/
                            T_it = ( t_u + lv / cp_l * c.x[ i ][ j ][ k ] - lv / cp_l * q_Rain );
                            E_Rain = hp * exp_func ( T_it, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                            q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                            q_Rain = .5 * ( q_Rain_n + q_Rain );  // smoothing the iteration process

                            c.x[ i ][ j ][ k ] = q_Rain; // water vapour restricted to saturated water vapour amount
                            cloud.x[ i ][ j ][ k ] = q_T - c.x[ i ][ j ][ k ]; // cloud water amount
                            ice.x[ i ][ j ][ k ] = 0.; // no cloud ice available
                            q_T = c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ];

                            if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
                            if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;

//                            if ( fabs ( q_Rain / q_Rain_n - 1. ) <= 1.e-5 )     break;

                            if( ( q_Rain_n ) > std::numeric_limits<double>::epsilon() &&
                                fabs ( q_Rain / q_Rain_n - 1. ) <= 1.e-5 )    break;  // make sure q_Rain_n is not 0 divisor

                            q_Rain_n = q_Rain;
                        }
                    }
                    cn.x[ i ][ j ][ k ] = c.x[ i ][ j ][ k ];
                    cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ];
                    icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ];
                    t.x[ i ][ j ][ k ] = T_it / t_0;
                } // end ( t_Celsius > 0. )
/** %%%%%%%%%%%%%%%%%%%%%%%%%%%     end          warm cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% **/


/** %%%%%%%%%%%%%%%%%%%%%%%%%%%     mixed cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% **/

/** %%%%%%%%%%%%%%%%%%%%%%%%%%     saturation pressure     %%%%%%%%%%%%%%%%%%%%%%%%%%%% **/
// mixed cloud phase, if 0°C > t > -37°C
                if ( t_Celsius < 0. ){
                    if ( t_Celsius < t_Celsius_2 )  cloud.x[ i ][ j ][ k ] = 0.;
                    if ( t_Celsius > 0. )  ice.x[ i ][ j ][ k ] = 0.;
                    float q_v_b = c.x[ i ][ j ][ k ];
                    float q_c_b = cloud.x[ i ][ j ][ k ];
                    float q_i_b = ice.x[ i ][ j ][ k ];
                    q_T = q_v_b + q_c_b + q_i_b; // total water content
/*
if ( ( i == 13 ) && ( j == 90 ) && ( k == 180 ) ){
    logger() << "   q_v_b = " << q_v_b * 1000. << "   q_c_b = " << q_c_b * 1000.
        << "   q_i_b = " << q_i_b * 1000. << "   q_T = " << q_T * 1000. << endl << endl;
}
*/
                    t_u = t.x[ i ][ j ][ k ] * t_0; // in K
                    float T = t_u; // in K

                    E_Rain = hp * exp_func ( T, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                    E_Ice = hp * exp_func ( T, 21.8746, 7.66 ); // saturation water vapour pressure for the ice phase in hPa
                    q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                    q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg

                    if ( ( q_c_b > 0. ) && ( q_i_b > 0. ) )
                        q_v_hyp = ( q_c_b * q_Rain + q_i_b * q_Ice ) / ( q_c_b + q_i_b );
                    if ( ( q_c_b >= 0. ) && ( q_i_b == 0. ) )  q_v_hyp = q_Rain;
                    if ( ( q_c_b == 0. ) && ( q_i_b > 0. ) )  q_v_hyp = q_Ice;

/** §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ **/

                    for(int iter_prec = 1; iter_prec <= 20; iter_prec++ ){ // iter_prec may be varied

/*
if ( ( i == 13 ) && ( j == 90 ) && ( k == 180 ) ){
    logger() << "mixed cloud               Ice_Water_Saturation_Adjustment: temperature max: "
        << (t.max() - 1)*t_0 <<"          iter_prec: " << iter_prec << std::endl;
    logger() << "mixed cloud               Ice_Water_Saturation_Adjustment: water vapour max: "
        << c.max() * 1000. << std::endl;
    logger() << "mixed cloud               Ice_Water_Saturation_Adjustment: cloud water max: "
        << cloud.max() * 1000. << std::endl;
    logger() << "mixed cloud               Ice_Water_Saturation_Adjustment: cloud ice max: "
        << ice.max() * 1000. << std::endl << std::endl;
}
*/
/** condensation == water vapor saturation for cloud water formation, deposition == ice crystal for cloud ice formation **/
                        // t_0 = 273.15 K == 0 °C, t_00 = 236.15 K == -37 °C
                        float CND = ( T - t_00 ) / ( t_0 - t_00 );
                        if ( T < t_00 ) CND = 0.;
                         // T = t_00 => CND = 0 ( no condensation == no cloud water ),
                         // T = t_0 => CND = 1 ( max condensation == max cloud water )

                        float DEP = ( t_0 - T ) / ( t_0 - t_00 );
                        if ( T > t_0 ) DEP = 0.;
                        // T = t_0 => DEP = 0 ( no deposition == no cloud ice ),
                        // T = t_00 => DEP = 1 ( max deposition == max cloud ice )

                        float d_q_v = q_v_hyp - q_v_b;  // changes in water vapour causing cloud water and cloud ice
                        float d_q_c = - d_q_v * CND;
                        float d_q_i = - d_q_v * DEP;

                        float d_t = ( lv * d_q_c + ls * d_q_i ) / cp_l; // in K, temperature changes
                        T = T + d_t; // in K

                        q_v_b = c.x[ i ][ j ][ k ] + d_q_v;  // new values
                        q_c_b = cloud.x[ i ][ j ][ k ] + d_q_c;
                        q_i_b = ice.x[ i ][ j ][ k ] + d_q_i;

                        if ( q_v_b < 0. )  q_v_b = 0.;  // negative values excluded, when iteration starts
                        if ( q_c_b < 0. )  q_c_b = 0.;
                        if ( q_i_b < 0. )  q_i_b = 0.;

                        p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa, needed for saturation water vapour
                        hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
                        if ( i != 0 ){
                            if( T > gam * hight * 1.e-2){
                                p_h = pow ( ( ( T - gam * hight * 1.e-2 ) / ( T ) ), exp_pressure ) * p_SL; // given in hPa
                            }else{
                                logger()<<"WARNING: T is less than gam * hight * 1.e-2. "<< __LINE__<<" "
                                    << __FILE__<<std::endl; 
                                p_h = p_SL;
                            }
                        }
                        else  p_h = p_SL;

                        E_Rain = hp * exp_func ( T, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        E_Ice = hp * exp_func ( T, 21.8746, 7.66 ); // saturation water vapour pressure for the ice phase in hPa
                        q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                        q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg

                        if ( ( q_c_b > 0. ) && ( q_i_b > 0. ) )
                            q_v_hyp = ( q_c_b * q_Rain + q_i_b * q_Ice ) / ( q_c_b + q_i_b );
                        // average amount of ater vapour based on temperature changes
                        if ( ( q_c_b >= 0. ) && ( q_i_b == 0. ) )  q_v_hyp = q_Rain;
                        if ( ( q_c_b == 0. ) && ( q_i_b > 0. ) )  q_v_hyp = q_Ice;

                        // rate of condensating or evaporating water vapour to form cloud water, 0.5 given by COSMO
//                        S_c_c.x[ i ][ j ][ k ] = .5 * d_q_c / dt_dim;
                        S_c_c.x[ i ][ j ][ k ] = .5 * ( cn.x[ i ][ j ][ k ] - c.x[ i ][ j ][ k ] ) / dt_dim;
                        if ( is_land ( h, i, j, k ) )  S_c_c.x[ i ][ j ][ k ] = 0.;

                        q_T = q_v_b + q_c_b + q_i_b; // total water content, not used except for print out for mass conservation test
/*
if ( ( i == 13 ) && ( j == 90 ) && ( k == 180 ) ){
    logger() << "   iter_prec = " << iter_prec << endl;
    logger() << "   CND = " << CND << "   DEP = " << DEP << "   d_t = " << d_t << endl;
    logger() << "   d_q_v = " << d_q_v * 1000. << "   d_q_c = " << d_q_c * 1000.
        << "   d_q_i = " << d_q_i * 1000. << "   q_v_hyp = " << q_v_hyp * 1000.
        << "   T = " << T << endl << endl;
    logger() << "   q_v_b = " << q_v_b * 1000. << "   q_c_b = " << q_c_b * 1000.
        << "   q_i_b = " << q_i_b * 1000. << "   q_T = " << q_T * 1000. << endl
        << "   ( d_q_v + d_q_c + d_q_i ) = " << ( d_q_v + d_q_c + d_q_i ) * 1000.
        << endl << "   fabs ( q_v_b / q_v_hyp - 1. ) = " << fabs ( q_v_b / q_v_hyp - 1. )
        << endl << endl;
}
*/
                        if( iter_prec >= 3 && (q_v_hyp) > std::numeric_limits<double>::epsilon() &&
                            fabs ( q_v_b / q_v_hyp - 1. ) <= 1.e-5 )    break;  // make sure q_v_hyp is not 0 divisor

                        q_v_b = .5 * ( q_v_hyp + q_v_b );  // has smoothing effect
                    } // iter_prec end

/** §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     end          iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ **/

                    cn.x[ i ][ j ][ k ] = c.x[ i ][ j ][ k ] = q_v_b;  // new values achieved after converged iterations
                    cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ] = q_c_b;
                    icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ] = q_i_b;

                    if ( t_Celsius < t_Celsius_2 )     cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ] = 0.;
                    if ( t_Celsius > 0. )                     icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ] = 0.;
                    t.x[ i ][ j ][ k ] = T / t_0;
                } // end ( ( t_Celsius < 0. ) && ( t_Celsius >= t_Celsius_2 ) )
            } // end i
        } // end j
    } // end k

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
                if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;  // in case negative values appear
                if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
                if ( ice.x[ i ][ j ][ k ] < 0. )  ice.x[ i ][ j ][ k ] = 0.;
            }
        }
    }
/*
    logger() << "end Ice_Water_Saturation_Adjustment: water vapour max: "
        << c.max() * 1000. << std::endl;
    logger() << "end Ice_Water_Saturation_Adjustment: cloud water max: "
        << cloud.max() * 1000. << std::endl;
    logger() << "end Ice_Water_Saturation_Adjustment: cloud ice max: "
        << ice.max() * 1000. << std::endl << std::endl;

    logger() << "end %%%%%%%%%%%% Ice_Water_Saturation_Adjustment: temperature max: "
    << (t.max() - 1)*t_0 << std::endl << std::endl << std::endl;
*/
    if(m_model->debug){
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
        assert(!c.has_nan());
        assert(!t.has_nan());
    }
}






void BC_Thermo::Two_Category_Ice_Scheme ( Array &h, Array &c, Array &t, Array &p_stat, 
        Array &cloud, Array &ice, Array &P_rain, Array &P_snow, Array &S_v, Array &S_c, Array &S_i, Array &S_r, 
        Array &S_s, Array &S_c_c ){
    if(m_model->debug){
        assert(!c.has_nan());
        assert(!t.has_nan());
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
    }
    //  Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation distribution formed of rain and snow
    // constant coefficients for the transport of cloud water and cloud ice amount vice versa, rain and snow in the parameterization procedures
    double N_i_0 = 1.e2;  // in m-3
    double m_i_0 = 1.e-12;  // in kg
    double m_i_max = 1.e-9;  // in kg
    double m_s_0 = 3.e-9;  // in kg

    double c_i_dep = 1.3e-5;  // in m3/(kg*s)
    double c_c_au = 4.e-4;  // in 1/s
    double c_i_au = 1.e-3;  // in 1/s
    double c_ac = .24;  // m2/kg
    double c_rim = 18.6;  // m2/kg
    double c_agg = 10.3;  // m2/kg
    double c_i_cri = .24;  // m2
    double c_r_cri = 3.2e-5;  // m2
    double a_ev = 1.e-3;  // m2/kg
    double b_ev = 5.9;  // m2*s/kg
    double c_s_dep = 1.8e-2;  // m2/kg
    double b_s_dep = 12.3;  // m2*s/kg
    double c_s_melt = 8.43e-5;  // (m2*s)/(K*kg)
    double b_s_melt = 12.05;  // m2*s/kg
    double a_s_melt = 2.31e3; // K/(kg/kg)
    double c_r_frz = 3.75e-2;  // (m2*s)/(K*kg)

    double t_nuc = 267.15;  // in K    -6 °C
    double t_d = 248.15;  // in K    -25 °C
    double t_hn = 236.15;  // in K    -40 °C
    double t_r_frz = 271.15;  // in K    -2 °C

    double p_t_in = 0.;
    double E_Rain_t_in = 0.;
    double q_Rain_t_in = 0.;

//    double t_1 = 253.15;  // in K    -20 °C
//    double t_00 = 236.15;  // in K    -40 °C
//    double t_Celsius_1 = t_1 - t_0;  // -20 °C
//    double t_Celsius_2 = t_00 - t_0;  // -37 °C

    double exp_pressure = g / ( 1.e-2 * gam * R_Air );

    //The m_i is only used in this function and I see no reason it should be a member of this class.
    //So, I move it out of class.
    //I also don't see any reason the above horde of variables should stay as class members. -- mchin
    double m_i = m_i_max; //initialize m_i local variable. If uninitialized, the value in m_i can be anything(bad, very bad). 

// rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
// the choosen scheme is a Two Category Ice Scheme
// besides the transport equation for the water vapour exists two equations for the cloud water and the cloud ice transport
// since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
                if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
                if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
                if ( ice.x[ i ][ j ][ k ] < 0. )  ice.x[ i ][ j ][ k ] = 0.;
                if ( P_rain.x[ i ][ j ][ k ] < 0. )  P_rain.x[ i ][ j ][ k ] = 0.;
                if ( P_snow.x[ i ][ j ][ k ] < 0. )  P_snow.x[ i ][ j ][ k ] = 0.;
            }
        }
    }

/******************* initial values for rain and snow calculation *********************/

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            P_rain.x[ im-1 ][ j ][ k ] = 0.;
            P_snow.x[ im-1 ][ j ][ k ] = 0.;
            S_r.x[ im-1 ][ j ][ k ] = 0.;
            S_s.x[ im-1 ][ j ][ k ] = 0.;

            for ( int i = im-2; i >= 0; i-- ){
                if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
                if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
                if ( ice.x[ i ][ j ][ k ] < 0. )  ice.x[ i ][ j ][ k ] = 0.;

                t_Celsius = t.x[ i ][ j ][ k ] * t_0 - t_0;
                t_u = t.x[ i ][ j ][ k ] * t_0;

                p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

                if ( i != 0 )  p_h = pow ( ( ( t_u - gam * hight * 1.e-2 ) / ( t_u ) ), exp_pressure ) * p_SL;
                else  p_h = p_SL;

                r_dry = 100. * p_h / ( R_Air * t_u );  // density of dry air in kg/m³
                r_humid = r_dry * ( 1. + c.x[ i ][ j ][ k ] ) / ( 1. + R_WaterVapour / R_Air * c.x[ i ][ j ][ k ] );
                q_h = c.x[ i ][ j ][ k ];  // threshold value for water vapour at local hight h in kg/kg
                E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 );  // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 );  // saturation water vapour pressure for the ice phase in hPa
                q_Rain  = ep * E_Rain / ( p_h - E_Rain );  // water vapour amount at saturation with water formation in kg/kg
                q_Ice  = ep * E_Ice / ( p_h - E_Ice );  // water vapour amount at saturation with ice formation in kg/kg

                if ( ( t_u >= t_0 ) && ( c_c_au * cloud.x[ i ][ j ][ k ] > 0. ) )
                    S_c_au = c_c_au * cloud.x[ i ][ j ][ k ];
                   // cloud water to rain, cloud droplet collection in kg / ( kg * s )
                else  S_c_au = 0.;

                if ( ( t_u < t_0 ) && ( c_i_au * ice.x[ i ][ j ][ k ] > 0. ) )
                    S_i_au = c_i_au * ice.x[ i ][ j ][ k ];     // cloud ice to snow, cloud ice crystal aggregation
                else  S_i_au = 0.;

// ice and snow average size
                if ( t_u <= t_0 )  N_i = N_i_0 * exp ( .2 * ( t_0 - t_u ) );

                if ( t_u <= t_0 ){
                    m_i = r_humid * ice.x[ i ][ j ][ k ] / N_i;
                    if ( m_i > m_i_max ) { m_i = m_i_max; }
                    else  m_i = r_humid * ice.x[ i ][ j ][ k ] / N_i;
                    if ( m_i < m_i_0 ) { m_i = m_i_0; }
                }

                if ( t_Celsius <= 0. ){
                    if ( c.x[ i ][ j ][ k ] > q_Ice ){  // supersaturation
                        S_i_dep = c_i_dep * N_i * pow ( m_i, ( 1. / 3. ) ) * ( c.x[ i ][ j ][ k ] - q_Ice );  // supersaturation, < III >
                    }
                    if ( ( c.x[ i ][ j ][ k ] < q_Ice ) && ( - ice.x[ i ][ j ][ k ] / dt_snow_dim >
                        ( c.x[ i ][ j ][ k ] - q_Ice ) / dt_snow_dim ) )
                        S_i_dep = - ice.x[ i ][ j ][ k ] / dt_snow_dim; // subsaturation, < III >
                }
                else  S_i_dep = 0.;


                S_r.x[ i ][ j ][ k ] = S_c_au;  // in kg / ( kg * s )
                S_s.x[ i ][ j ][ k ] = S_i_au;
                if ( ( is_land ( h, i, j, k ) ) && ( is_land ( h, i+1, j, k ) ) ){
                    S_r.x[ i ][ j ][ k ] = 0.;
                    S_s.x[ i ][ j ][ k ] = 0.;
                }

                if ( P_rain.x[ i + 1 ][ j ][ k ] < 0. )  P_rain.x[ i + 1 ][ j ][ k ] = 0.;
                if ( P_snow.x[ i + 1 ][ j ][ k ] < 0. )  P_snow.x[ i + 1 ][ j ][ k ] = 0.;

                /*
                float step = m_model->get_layer_height(i+1) - m_model->get_layer_height(i);
                P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ]
                     + r_humid * ( ( S_r.x[ i ][ j ][ k ] - S_r.x[ i + 1 ][ j ][ k ] )
                     / 2. * step ) / 2.;  // in kg / ( m2 * s ) == mm/s
                P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ]
                     + r_humid * ( ( S_s.x[ i ][ j ][ k ] - S_s.x[ i + 1 ][ j ][ k ] )
                     / 2. * step ) / 2.;
                */
                P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ] +
                    ( S_r.x[ i ][ j ][ k ] + S_r.x[ i + 1 ][ j ][ k ] ) * .5 * r_humid * dr * 200.;  // in kg / ( m2 * s ) == mm/s
                P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ] +
                    ( S_s.x[ i ][ j ][ k ] + S_s.x[ i + 1 ][ j ][ k ] ) * .5 * r_humid * dr * 200.;
                      // arbitrary factor 200 adjusts precipitation for the modern world.

                if ( P_rain.x[ i ][ j ][ k ] < 0. )  P_rain.x[ i ][ j ][ k ] = 0.;
                if ( P_snow.x[ i ][ j ][ k ] < 0. )  P_snow.x[ i ][ j ][ k ] = 0.;
            }
        }
    }


/******************* main part for rain and snow calculation *********************/

    if ( true ){
        for(int iter_prec = 1; iter_prec <= 5; iter_prec++ ){ // iter_prec may be varied, but is sufficient
            for ( int k = 0; k < km; k++ ){
                for ( int j = 0; j < jm; j++ ){
                    P_rain.x[ im-1 ][ j ][ k ] = 0.;
                    P_snow.x[ im-1 ][ j ][ k ] = 0.;

                    for ( int i = im-2; i >= 0; i-- ){
                        t_u = t.x[ i ][ j ][ k ] * t_0;
                        t_Celsius = t_u - t_0;

                        p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                        hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );

                        if ( i != 0 )  p_h = pow ( ( ( t_u - gam * hight * 1.e-2 ) / ( t_u ) ), exp_pressure ) * p_SL;
                        else  p_h = p_SL;

                        r_dry = 100. * p_h / ( R_Air * t_u );  // density of dry air in kg/m³
                        r_humid = r_dry * ( 1. + c.x[ i ][ j ][ k ] ) / ( 1. + R_WaterVapour / R_Air * c.x[ i ][ j ][ k ] );
                        E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 );  // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 );  // saturation water vapour pressure for the ice phase in hPa
                        q_Rain = ep * E_Rain / ( p_h - E_Rain );  // water vapour amount at saturation with water formation in kg/kg
                        q_Ice  = ep * E_Ice / ( p_h - E_Ice );  // water vapour amount at saturation with ice formation in kg/kg

// ice and snow average size
                        if ( t_u <= t_0 )  N_i = N_i_0 * exp ( .2 * ( t_0 - t_u ) );

                        if ( t_u <= t_0 ){
                            m_i = r_humid * ice.x[ i ][ j ][ k ] / N_i;
                            if ( m_i > m_i_max ) { m_i = m_i_max; }
                            else  m_i = r_humid * ice.x[ i ][ j ][ k ] / N_i;
                            if ( m_i < m_i_0 ) { m_i = m_i_0; }
                        }

// nucleation and depositional growth of cloud ice
                        if ( ice.x[ i ][ j ][ k ] == 0. ){
                            if ( ( t_u < t_d ) && ( c.x[ i ][ j ][ k ] >= q_Ice ) )
                                S_nuc = m_i_0 / ( r_humid * dt_snow_dim ) * N_i;  // nucleation of cloud ice, < I >

                            if ( ( ( t_d <= t_u ) && ( t_u <= t_nuc ) ) && ( c.x[ i ][ j ][ k ] >= q_Rain ) )
                                S_nuc = m_i_0 / ( r_humid * dt_snow_dim ) * N_i;  // nucleation of cloud ice, < I >
                        }
                        else  S_nuc = 0.;

                        if ( ( t_u < t_hn ) && ( cloud.x[ i ][ j ][ k ] > 0. ) )
                            S_c_frz = cloud.x[ i ][ j ][ k ] / dt_rain_dim;  //nucleation of cloud ice due to freezing of cloud water, < II >
                        else  S_c_frz = 0.;

                        if ( t_Celsius <= 0. ){
                            if ( c.x[ i ][ j ][ k ] > q_Ice ){  // supersaturation
                                S_i_dep = c_i_dep * N_i * pow ( m_i, ( 1. / 3. ) ) * ( c.x[ i ][ j ][ k ] - q_Ice );  // supersaturation, < III >
                            }
                            if ( ( c.x[ i ][ j ][ k ] < q_Ice ) && ( - ice.x[ i ][ j ][ k ] / dt_snow_dim >
                                ( c.x[ i ][ j ][ k ] - q_Ice ) / dt_snow_dim ) )
                                S_i_dep = - ice.x[ i ][ j ][ k ] / dt_snow_dim; // subsaturation, < III >
                        }
                        else  S_i_dep = 0.;

// autoconversion processes
                        if ( ( t_u >= t_0 ) && ( cloud.x[ i ][ j ][ k ] > 0. ) )
                            S_c_au = c_c_au * cloud.x[ i ][ j ][ k ];  // cloud water to rain, cloud droplet collection, < IV >
                        else  S_c_au = 0.;

                        if ( ( t_u <= t_0 ) && ( ice.x[ i ][ j ][ k ] > 0. ) )
                            S_i_au = c_i_au * ice.x[ i ][ j ][ k ];  // cloud ice to snow, cloud ice crystal aggregation, < V >
                        else  S_i_au = 0.;

                        if ( t_u <= t_0 )
                            S_d_au = S_i_dep / ( 1.5 * ( pow ( ( m_s_0 / m_i ),  ( 2. / 3. ) ) - 1. ) );  // depositional growth of cloud ice, < VI >
                        else  S_d_au = 0.;

// collection mechanism
                        if ( t_u > t_0 )  S_ac = c_ac * cloud.x[ i ][ j ][ k ] * pow ( P_rain.x[ i ][ j ][ k ], ( 7. / 9. ) );
                            // accreation rate from depletion of cloud water due to collection by all rain drops, < VII >
                        else  S_ac = 0.;  // accreation rate from depletion of cloud water due to collection by all rain drops

                        if ( t_u < t_0 )  S_rim = c_rim * cloud.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];
                        else  S_rim = 0.;  // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >
                                                     // by falling snow particles

                        if ( t_u >= t_0 )  S_shed = c_rim * cloud.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];
                        else  S_shed = 0.;  // rate of water shed by melting wet snow particles, < IX >
                                                       // collecting cloud droplets to produce rain

                        if ( t_u <= t_0 ){
                            S_agg = c_agg * ice.x[ i ][ j ][ k ] * P_snow.x[ i ][ j ][ k ];  // collection of cloud ice by snow particles, < X >

                            S_i_cri = c_i_cri * ice.x[ i ][ j ][ k ] * pow ( P_rain.x[ i ][ j ][ k ], ( 7. / 9. ) );
                                // decrease in cloud ice mass due to collision/coalescense interaction with raindrops, < XI >
                            S_r_cri = c_r_cri * ice.x[ i ][ j ][ k ] / m_i * pow ( P_rain.x[ i ][ j ][ k ], ( 13. / 9. ) );
                                // decrease of rainwater due to freezing resulting from collection of ice crystals, < XII >
                        }else{
                            S_agg = 0.;
                            S_i_cri = 0.;
                            S_r_cri = 0.;
                        }

// diffusional growth of rain and snow
                        if ( t_u >= t_0 )  S_ev = a_ev * ( 1. + b_ev * pow ( P_rain.x[ i ][ j ][ k ],
                            ( 1. / 6. ) ) ) * ( q_Rain - c.x[ i ][ j ][ k ] ) * pow ( P_rain.x[ i ][ j ][ k ], ( 4. / 9. ) );
                            // evaporation of rain due to water vapour diffusion, < XIII >
                        else  S_ev = 0.;

                        if ( t_u < t_0 )  S_s_dep = c_s_dep * ( 1. + b_s_dep * pow ( P_snow.x[ i ][ j ][ k ],
                            ( 5. / 26. ) ) ) * ( c.x[ i ][ j ][ k ] - q_Ice ) * pow ( P_snow.x[ i ][ j ][ k ], ( 8. / 13. ) );
                             // deposition/sublimation of snow, < XIV >
                        else  S_s_dep = 0.;

// melting and freezing
                        if ( ( t_u > t_0 ) && ( ice.x[ i ][ j ][ k ] > 0. ) )
                            S_i_melt = ice.x[ i ][ j ][ k ] / dt_snow_dim; // cloud ice particles melting to cloud water, < XV >
                        else  S_i_melt = 0.;

                        if ( t_u > t_0 ){
                            p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );  // given in hPa
                            hight = ( double ) i * ( L_atm / ( double ) ( im-1 ) );
                            p_t_in = pow ( ( ( t_0 - gam * hight * 1.e-2 ) / t_0 ), exp_pressure ) * p_SL;  // given in hPa
                            E_Rain_t_in = hp * exp_func ( t_0, 17.2694, 35.86 );
                                // saturation water vapour pressure for the water phase at t = 0°C in hPa
                            q_Rain_t_in = ep * E_Rain_t_in / ( p_t_in - E_Rain_t_in );
                                // water vapour amount at saturation with water formation in kg/kg

                            S_s_melt = c_s_melt * ( 1. + b_s_melt * pow ( P_snow.x[ i ][ j ][ k ], ( 5. / 26. ) ) ) *
                                ( ( t_u - t_0 ) + a_s_melt * ( c.x[ i ][ j ][ k ] - q_Rain_t_in ) ) *
                                pow ( P_snow.x[ i ][ j ][ k ], ( 8. / 13. ) );  // melting rate of snow to form rain, < XVI >
                                // arbitrary factor 50 adjusts the average precipitation in mm/d
                        }
                        else  S_s_melt = 0.;

                        if ( t_r_frz -  t_u > 0. )  S_r_frz = c_r_frz * pow ( ( t_r_frz -  t_u ), ( 3. / 2. ) ) *
                            pow ( P_rain.x[ i ][ j ][ k ], ( 3. / 2. ) );
                            // immersion freezing and contact nucleation, < XVII >
                        else  S_r_frz = 0.;

// sinks and sources
                        S_v.x[ i ][ j ][ k ] = - S_c_c.x[ i ][ j ][ k ] + S_ev - S_i_dep -
                            S_s_dep - S_nuc;
                        S_c.x[ i ][ j ][ k ] = S_c_c.x[ i ][ j ][ k ] - S_c_au - S_ac - S_c_frz +
                            S_i_melt - S_rim - S_shed;
                        S_i.x[ i ][ j ][ k ] = S_nuc + S_c_frz + S_i_dep - S_i_melt -
                            S_i_au - S_d_au - S_agg - S_i_cri;
                        S_r.x[ i ][ j ][ k ] = S_c_au + S_ac - S_ev + S_shed - S_r_cri -
                            S_r_frz + S_s_melt;
                        S_s.x[ i ][ j ][ k ] = S_d_au + S_s_dep + S_i_au + S_rim +
                            S_agg + S_i_cri + S_r_cri + S_r_frz - S_s_melt;

                        if ( ( is_land ( h, i, j, k ) ) && ( is_land ( h, i+1, j, k ) ) ){
                            S_c_c.x[ i ][ j ][ k ] = 0.;
                            S_v.x[ i ][ j ][ k ] = 0.;
                            S_c.x[ i ][ j ][ k ] = 0.;
                            S_i.x[ i ][ j ][ k ] = 0.;
                            S_r.x[ i ][ j ][ k ] = 0.;
                            S_s.x[ i ][ j ][ k ] = 0.;
                        }

// rain and snow integration
                        if ( P_rain.x[ i + 1 ][ j ][ k ] < 0. )  P_rain.x[ i + 1 ][ j ][ k ] = 0.;
                        if ( P_snow.x[ i + 1 ][ j ][ k ] < 0. )  P_snow.x[ i + 1 ][ j ][ k ] = 0.;

                        /*
                        float step = m_model->get_layer_height(i+1) - m_model->get_layer_height(i-1);
                        P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ]
                             + r_humid * ( ( S_r.x[ i ][ j ][ k ] - S_r.x[ i + 1 ][ j ][ k ] )
                             / 2. * step ) / 2.;  // in kg / ( m2 * s ) == mm/s
                        P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ]
                             + r_humid * ( ( S_s.x[ i ][ j ][ k ] - S_s.x[ i + 1 ][ j ][ k ] )
                             / 2. * step ) / 2.;
                        */
                        P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ] +
                            ( S_r.x[ i ][ j ][ k ] + S_r.x[ i + 1 ][ j ][ k ] ) * .5 * r_humid * dr * 200.;  // in kg / ( m2 * s ) == mm/s
                        P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ] +
                            ( S_s.x[ i ][ j ][ k ] + S_s.x[ i + 1 ][ j ][ k ] ) * .5 * r_humid * dr * 200.;
                      // arbitrary factor 200 adjusts precipitation for the modern world.
                        if ( P_rain.x[ i ][ j ][ k ] < 0. )  P_rain.x[ i ][ j ][ k ] = 0.;
                        if ( P_snow.x[ i ][ j ][ k ] < 0. )  P_snow.x[ i ][ j ][ k ] = 0.;

                        if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
                        if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
                        if ( ice.x[ i ][ j ][ k ] < 0. )  ice.x[ i ][ j ][ k ] = 0.;
                    }  // end i RainSnow
                }  // end j
            }  // end k
        }  // end iter_prec
    }  // end n
    if(m_model->debug){
        assert(!c.has_nan());
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
        assert(!P_snow.has_nan());
        assert(!P_rain.has_nan());
    }
}







void BC_Thermo::IC_Temperature_WestEastCoast ( Array &h, Array &t ){
// initial conditions for the temperature close to coast sides to damp out shades of preceeding timeslices
    j_grad = 7;                                                                             // extension for temperature change in zonal direction
    k_grad = 7;                                                                             // extension for temperature change in longitudinal direction

// search for north coasts to smooth the temperature

// northern and southern hemisphere: north coast
    j_air = 0;                                                                            // somewhere on air
    flip = 0;                                                                                   // somewhere on air
    for ( int k = 1; k < km-1; k++ ){                                            // outer loop: longitude
        for ( int j = j_grad; j < jm-1; j++ ){                                   // inner loop: latitude
            if ( is_air ( h, 0, j, k ) ){                                             // if somewhere on air
                j_air = 0;                                                                // somewhere on air: j_air = 0
                flip = 0;                                                                       // somewhere on air: flip = 0
            }
            else j_air = 1;                                                           // first time on land
            if ( ( flip == 0 ) && ( j_air == 1 ) ){                            // on air closest to land
                ll = j - j_grad;                                                                // starting point of temperature smoothing
                for ( int l = ll; l < j; l++ )      t.x[ 0 ][ l ][ k ] = t.x[ 0 ][ ll ][ k ];       // replacement of temperature values
                flip = 1;                                                                       // somewhere on land: flip = 1
            }
        }                                                                                           // end of latitudinal loop
        flip = 0;                                                                               // somewhere on air: flip = 0
    }                                                                                               // end of longitudinal loop

// northern and southern hemisphere: south coast
    j_air = 0;                                                                            // on air closest to coast
    j_sequel = 1;                                                                           // on solid ground
    for ( int k = 1; k < km-1; k++ ){                                           // outer loop: latitude
        for ( int j = 0; j < jm - j_grad; j++ ){                                 // inner loop: longitude
            if ( is_land ( h, 0, j, k ) ) j_sequel = 0;                       // if solid ground: j_sequel = 0
            if ( ( is_air ( h, 0, j, k ) ) && ( j_sequel == 0 ) ) j_air = 0;   // if air and and j_sequel = 0 then is air closest to coast
            else j_air = 1;                                                           // somewhere on air
            if ( ( is_air ( h, 0, j, k ) ) && ( j_air == 0 ) ){     // if air is closest to coast, change of velocity components begins
                ll = j + j_grad;                                                            // starting point of temperature smoothing
                for ( int l = ll; l > j; l-- )      t.x[ 0 ][ l ][ k ] = t.x[ 0 ][ ll ][ k ];       // replacement of temperature values
                j_sequel = 1;                                                               // looking for another south coast
            }
        }                                                                                           // end of longitudinal loop
        j_air = 0;                                                                        // starting at another latitude
    }                                                                                               // end of latitudinal loop

// northern and southern hemisphere: east coast
    k_air = 0;                                                                            // on air closest to coast
    k_sequel = 1;                                                                           // on solid ground
    for ( int j = 1; j < jm-1; j++ ){                                                // outer loop: latitude
        for ( int k = k_grad; k < km-k_grad; k++ ){                      // inner loop: longitude
            if ( is_land ( h, 0, j, k ) ) k_sequel = 0;                       // if solid ground: k_sequel = 0
            if ( ( is_air ( h, 0, j, k ) ) && ( k_sequel == 0 ) ) k_air = 0;        // if air and and k_sequel = 0 then air lies closest to a coast
            else k_air = 1;                                                           // somewhere on air
            if ( ( is_air ( h, 0, j, k ) ) && ( k_air == 0 ) ){     // if air is closest to coast, change of velocity components begins
                ll = k + k_grad;                                                            // starting point of temperature smoothing
                for ( int l = ll; l > k; l-- )      t.x[ 0 ][ j ][ l ] = t.x[ 0 ][ j ][ ll ];       // replacement of temperature values
                k_sequel = 1;                                                               // looking for another east coast
            }
        }                                                                                           // end of longitudinal loop
        k_air = 0;                                                                        // starting at another longitude
    }                                                                                               // end of latitudinal loop

// northern and southern hemisphere: west coast
    k_air = 0;                                                                            // somewhere on air
    flip = 0;                                                                                   // somewhere on air
    for ( int j = 1; j < jm-1; j++ ){                                                // outer loop: latitude
        for ( int k = k_grad; k < km-1; k++ ){                               // inner loop: longitude
            if ( is_air ( h, 0, j, k ) ){                                             // if somewhere on air
                k_air = 0;                                                                // somewhere on air: k_air = 0
                flip = 0;                                                                       // somewhere on air: flip = 0
            }
            else k_air = 1;                                                           // first time on land
            if ( ( flip == 0 ) && ( k_air == 1 ) ){                            // on air closest to land
                ll = k - k_grad;                                                            // starting point of temperature smoothing
                for ( int l = ll; l < k; l++ )      t.x[ 0 ][ j ][ l ] = t.x[ 0 ][ j ][ ll ];       // replacement of temperature values
                flip = 1;                                                                       // somewhere on land: flip = 1
            }
        }                                                                                           // end of longitudinal loop
        flip = 0;                                                                               // somewhere on air: flip = 0
    }                                                                                               // end of latitudinal loop
}





void BC_Thermo::Value_Limitation_Atm ( Array &h, Array &u, Array &v, Array &w,
                            Array &p_dyn, Array &t, Array &c, Array &cloud, Array &ice, Array &co2 ){
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
                if ( u.x[ i ][ j ][ k ] >= .106 )  u.x[ i ][ j ][ k ] = .106;
                if ( u.x[ i ][ j ][ k ] <= - .106 )  u.x[ i ][ j ][ k ] = - .106;
                if ( v.x[ i ][ j ][ k ] >= .125 )  v.x[ i ][ j ][ k ] = .125;
                if ( v.x[ i ][ j ][ k ] <= - .125 )  v.x[ i ][ j ][ k ] = - .125;
                if ( w.x[ i ][ j ][ k ] >= 3.5 )  w.x[ i ][ j ][ k ] = 3.5;
                if ( w.x[ i ][ j ][ k ] <= - .469 )  w.x[ i ][ j ][ k ] = - .469;
                if ( t.x[ i ][ j ][ k ] >= 1.165 )  t.x[ i ][ j ][ k ] = 1.165;  // == 45 °C
                if ( t.x[ i ][ j ][ k ] <= - .78 )  t.x[ i ][ j ][ k ] = - .78;  // == 59.82 °C
//                if ( c.x[ i ][ j ][ k ] >= .022 )  c.x[ i ][ j ][ k ] = .022;
                if ( c.x[ i ][ j ][ k ] >= .03 )  c.x[ i ][ j ][ k ] = .03;
                if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
//                if ( cloud.x[ i ][ j ][ k ] >= .008 )  cloud.x[ i ][ j ][ k ] = .008;
                if ( cloud.x[ i ][ j ][ k ] >= .01 )  cloud.x[ i ][ j ][ k ] = .01;
                if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
//                if ( ice.x[ i ][ j ][ k ] >= .0025 )  ice.x[ i ][ j ][ k ] = .0025;
                if ( ice.x[ i ][ j ][ k ] >= .005 )  ice.x[ i ][ j ][ k ] = .005;
                if ( ice.x[ i ][ j ][ k ] < 0. )  ice.x[ i ][ j ][ k ] = 0.;
                if ( co2.x[ i ][ j ][ k ] >= 5.36 )  co2.x[ i ][ j ][ k ] = 5.36;
                if ( co2.x[ i ][ j ][ k ] <= 1. )  co2.x[ i ][ j ][ k ] = 1.;

                if ( is_land ( h, i, j, k ) ){
                    u.x[ i ][ j ][ k ] = 0.;
                    v.x[ i ][ j ][ k ] = 0.;
                    w.x[ i ][ j ][ k ] = 0.;
//                    t.x[ i ][ j ][ k ] = 1.;  // = 273.15 K
//                    c.x[ i ][ j ][ k ] = 0.;
                    cloud.x[ i ][ j ][ k ] = 0.;
                    ice.x[ i ][ j ][ k ] = 0.;
                    co2.x[ i ][ j ][ k ] = 1.;  // = 280 ppm
                    p_dyn.x[ i ][ j ][ k ] = 0.;
                }
            }
        }
    }
}




void BC_Thermo::Pressure_Limitation_Atm ( Array &p_dyn, Array &p_dynn ){
// class element for the limitation of flow properties, to avoid unwanted growth around geometrical singularities
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
                if ( p_dyn.x[ i ][ j ][ k ] >= .25 )  p_dyn.x[ i ][ j ][ k ] = .25;
                if ( p_dyn.x[ i ][ j ][ k ] <= - .25 )  p_dyn.x[ i ][ j ][ k ] = - .25;

                if ( is_land ( h, i, j, k ) ){
                    p_dyn.x[ i ][ j ][ k ] = 0.;
                }
                p_dynn.x[ i ][ j ][ k ] = p_dyn.x[ i ][ j ][ k ];
            }
        }
    }
}



int BC_Thermo::GetTropopauseHightAdd(double t_cret){
    double d_i_h_round = round((t_cret * t_0) / 2.6);
        // adiabatic slope of radial temperature 0.65/100m, stepsize 400m => 2.6/400m
    return ( int ) d_i_h_round;
}


double get_pole_temperature(int Ma, int Ma_1, int Ma_2, double t_1, double t_2){
    return (t_2 - t_1) / (double) (Ma_2 - Ma_1) * (double) (Ma - Ma_1) + t_1;
}

double get_pole_temperature(int Ma, const std::map<int, double> &pole_temp_map){
    assert(pole_temp_map.size()>1);
    
    std::pair<int, double> up = *pole_temp_map.begin(), bottom = *++pole_temp_map.begin();
    
    // when Ma out of boundary
    if(Ma <= pole_temp_map.begin()->first){
        return pole_temp_map.begin()->second; 
    }else if(Ma > (--pole_temp_map.end())->first){
        return (--pole_temp_map.end())->second;
    }

    for( const auto& pair : pole_temp_map ){
        if(pair.first>=Ma){
            bottom = pair;
            break;
        }else{
            up = pair;
        }
    }
    return get_pole_temperature(Ma, up.first, bottom.first, up.second, bottom.second);
}

