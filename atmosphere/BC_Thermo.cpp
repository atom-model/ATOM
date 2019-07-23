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

using namespace std;
using namespace AtomUtils;

bool debug = false;

BC_Thermo::BC_Thermo (cAtmosphereModel* model, int im, int jm, int km, double c_0, double c_land, double t_land, double co2_0, Array& h): 
        m_model(model),
        im(im),
        jm(jm),
        km(km),
        h(h),
        i_topography(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
        t_pal(std::vector<double>(im, 0))
{
    this -> tropopause_equator = model->tropopause_equator;
    this -> tropopause_pole = model->tropopause_pole;
    this-> L_atm = model->L_atm;
    this-> dt = model->dt;
    this-> dr = model->dr;
    this-> dthe = model->dthe;
    this-> dphi = model->dphi;
    this-> RadiationModel = model->RadiationModel;
    this-> NASATemperature = model->NASATemperature;
    this-> sun = model->sun;
    this-> g =  model->g;
    this-> ep =  model->ep;
    this-> hp =  model->hp;
    this-> u_0 =  model->u_0;
    this-> p_0 =  model->p_0;
    this-> t_0 =  model->t_0;
    this-> c_0 =  model->c_0;
    this-> co2_0 =  model->co2_0;
    this-> sigma =  model->sigma;
    this-> albedo_equator =  model->albedo_equator;
    this-> albedo_pole =  model->albedo_pole;
    this-> gam =  model->gam;
    this-> lv =  model->lv;
    this-> ls =  model->ls;
    this-> cp_l =  model->cp_l;
    this-> r_air =  model->r_air;
    this-> R_Air =  model->R_Air;
    this-> r_water_vapour =  model->r_water_vapour;
    this-> R_WaterVapour =  model->R_WaterVapour;
    this-> R_co2 =  model->R_co2;
    this-> co2_paleo =  model->co2_paleo;
    this-> co2_vegetation =  model->co2_vegetation;
    this-> co2_ocean =  model->co2_ocean;
    this-> co2_land =  model->co2_land;
    this-> emissivity_add =  model->emissivity_add;
    this-> rad_equator =  model->rad_equator;
    this-> rad_pole =  model->rad_pole;
    this-> irr =  model->irr;
    this-> epsilon_pole =  model->epsilon_pole;
    this-> epsilon_tropopause =  model->epsilon_tropopause;
    this-> epsilon_equator =  model->epsilon_equator;
    this-> c_tropopause =  model->c_tropopause;
    this-> co2_tropopause =  model->co2_tropopause;
    this-> c_ocean =  model->c_ocean;
    this-> c_land =  model->c_land;
    this-> t_average =  model->t_average;
    this-> co2_average =  model->co2_average;
    this-> co2_pole =  model->co2_pole;
    this-> co2_equator =  model->co2_equator;
    this-> t_land =  model->t_land;
    this-> t_tropopause =  model->t_tropopause;
    this-> t_equator =  model->t_equator;
    this-> t_pole =  model->t_pole;
    this-> declination =  model->declination;
    this-> sun_position_lat =  model->sun_position_lat;
    this-> sun_position_lon =  model->sun_position_lon;
    this-> zeta = model->zeta;

    im_tropopause = model->get_tropopause();

    Ma = int(round(*m_model->get_current_time()));

    t_paleo = m_model->get_mean_temperature_from_curve(Ma) -
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




void BC_Thermo::BC_Radiation_multi_layer ( Array_1D &rad, Array_2D &albedo, Array_2D &epsilon,
                             Array_2D &radiation_surface, Array &p_stat, Array &t, Array &c,
                             Array &h, Array &epsilon_3D, Array &radiation_3D, Array &cloud,
                             Array &ice, Array &co2 ){
    // class element for the computation of the radiation and the temperature distribution
    // computation of the local temperature based on short and long wave radiation
    // multi layer radiation model
    cout.precision ( 4 );
    cout.setf ( ios::fixed );

    pi180 = 180./M_PI;

    j_half = ( jm -1 ) / 2;  // position of the sun at 0°S
    j_max = jm - 1;
    d_j_half = ( double ) j_half;
    d_j_max = ( double ) j_max;
    k_half = ( km -1 ) / 2;  // position of the sun at 0° oder 180° ( Greenwich )
    k_max = km - 1;
    d_k_half = ( double ) k_half;
    d_k_max = ( double ) k_max;
    double j_max_half = ( jm -1 ) / 2;

    rad_eff = rad_pole - rad_equator;
    albedo_co2_eff = albedo_pole - albedo_equator;

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
    epsilon_eff_max = .594; // constant  given by Häckel ( F. Baur and H. Philips, 1934 )
    // constant value stands for other non-condensable gases than water vapour in the equation for epsilon
    epsilon_eff_2D = epsilon_pole - epsilon_equator;
//    double x_c = 0.;
    double u_c = 0.;
    double eps_co2 = 0.;
    double P_c = 0.;
    double step = 0.;

    for ( int j = 0; j < jm; j++ ){
        i_trop = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0);

        height_tropo = ( exp( zeta * ( rad.z[ i_trop ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
        d_i_max = height_tropo;

        // on zero level, lateral parabolic distribution
        epsilon_eff_max = epsilon_eff_2D * parabola( j / j_max_half ) + epsilon_pole;

        for ( int k = 0; k < km; k++ ){
//            i_mount = i_topography[ j ][ k ];
            i_mount = 0;
/*
    // albedo computation
            if ( is_land( h, 0, j, k ) )
                albedo.y[ j ][ k ] = albedo_co2_eff * parabola( j / j_max_half ) + albedo_pole;

            if ( ( is_land( h, 0, j, k ) ) && ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 <= 0. ) )  albedo.y[ j ][ k ] = .95;  // assumption of 75% reflection
            if ( ( is_land( h, 0, j, k ) ) && ( t.x[ 0 ][ j ][ k ] * t_0 - t_0 >= 24. ) )  albedo.y[ j ][ k ] = .3;  // assumption for vegetation
            if ( is_land( h, 0, j, k ) ){
                if ( ( ( j >= 60 ) && ( j <= 80 ) ) && ( ( k == 0 ) || ( k <= 60 ) ) )  albedo.y[ j ][ k ] = .40;  // assumption for desert land
                if ( ( ( j >= 60 ) && ( j <= 80 ) ) && ( ( k >= 340 ) && ( k < km ) ) )  albedo.y[ j ][ k ] = .40;  // assumption for desert land
            }
            if ( is_water( h, 0, j, k ) ){
                double albedo_pole_water = 0.03;
                double albedo_equator_water = .14;
                double albedo_co2_eff_water = albedo_pole_water - albedo_equator_water;
                albedo.y[ j ][ k ] = albedo_co2_eff_water * parabola( j / j_max_half ) + albedo_pole_water;
                albedo.y[ j ][ k ] = .14;  // assumption for reflecting water surface 0.03 - 0.1
            }
            radiation_surface.y[ j ][ k ] = irr * albedo.y[ j ][ k ];
*/
            // in W/m², assumption of parabolic surface radiation at zero level
            radiation_surface.y[ j ][ k ] = rad_eff * parabola( j / j_max_half ) + rad_pole;  // without albedo computation

            for ( int i = 0; i <= i_trop; i++ ){
                if ( c.x[ i ][ j ][ k ] < 0. )      c.x[ i ][ j ][ k ] = 0.;
                if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
                if ( ice.x[ i ][ j ][ k ] < 0. )    ice.x[ i ][ j ][ k ] = 0.;

                // COSMO water vapour pressure based on local water vapour, cloud water, cloud ice in hPa
                e = ( c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ] + ice.x[ i ][ j ][ k ] ) * p_stat.x[ i ][ j ][ k ] / ep;

                // radial parabolic distribution, start on zero level
                    height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                    d_i = height;
                    epsilon_eff = epsilon_eff_max - ( epsilon_tropopause - epsilon_eff_max ) *
                        ( d_i / d_i_max * ( d_i / d_i_max - 2. ) );// radial distribution approximated by a parabola

                // dependency given by Häckel, Meteorologie, p. 205 ( law by F. Baur and H. Philips, 1934 )
                // co2_coeff * epsilon_eff describe the effect in the emissivity computation of other gases like CO2
                // in the original formula this value is 0.594, for reasons of adjustment to the modern atmosphere,
                // this constant became a variable in zonal direction
                // this variable reacts very sensitive and changes the temperature field extremely
                // the second term describes the influence of water vapour only
                if( i >= i_mount ){ //start from the mountain top
                    step = ( ( exp( zeta * ( rad.z[ i + 1 ] - 1. ) ) - 1 ) 
                             - ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) ) 
                             * ( L_atm / ( double ) ( im-1 ) ); // local atmospheric shell thickness
                    P_c = ( 1.e-6 * p_stat.x[ i ][ j ][ k ] * co2.x[ i ][ j ][ k ] * co2_0 / p_0 );
                    u_c = P_c * step * 100.; // model by Atwater and Ball
//                    u_c = 300. / 0.9869 * P_c * step * 100. / ( t.x[ i ][ j ][ k ] * t_0 ); 
                            // modell by Yamamoto and Sasamori

                 // influence on the emissivity by carbon dioxcide by the law by Bliss
//                    x_c = 10. * P_c / ( R_co2 * t.x[ i ][ j ][ k ] * t_0 ) * step * p_0;
                  // coefficient in the Bliss law, 10. => conversion from kg/m² to g/cm²
//                    eps_co2 = .185 * ( 1. - exp ( - 50. * x_c ) ); // modell by Bliss
                    eps_co2 = .185 * ( 1. - exp ( - 0.3919 * pow ( u_c, 0.4 ) ) ); // model by Atwater and Ball
                    eps_co2 = .5 * eps_co2; // model by Atwater and Ball, 
                                            // factor .5 fits the Atwater and Ball results to Yamamoto and Sasamori,
                                            // which fit best to HITRAN band data
//                    eps_co2 = 0.; // Test

                    epsilon_3D.x[ i ][ j ][ k ] = emissivity_add * epsilon_eff + .0416 * sqrt ( e ) + eps_co2;

                    // atmospheric emissivity by Vogel/Bliss, emissivity seems too high
//                    epsilon_3D.x[ i ][ j ][ k ] = 0.98 * pow( ( e / ( t.x[ i ][ j ][ k ] * t_0 ) ), 0.0687 );

                    radiation_3D.x[ i ][ j ][ k ] = ( 1. - epsilon_3D.x[ i ][ j ][ k ] ) * sigma * 
                                    pow ( t.x[ i ][ j ][ k ] * t_0, 4. );
                }
                if ( epsilon_3D.x[ i ][ j ][ k ] > 1. )  epsilon_3D.x[ i ][ j ][ k ] = 1.;

                epsilon.y[ j ][ k ] = epsilon_3D.x[ 0 ][ j ][ k ];

/*
if ( ( j == 90 ) && ( k == 180 ) ) cout << "   i = " << i << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_average = " << t_average << "   t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "   co2_coeff = " << co2_coeff << "   co2_eff = " << ( co2_pole - co2_equator ) << "   co2_pole = " << co2_pole << "   co2_equator = " << co2_equator << "   co2_paleo = " << co2_paleo << "   co2_average = " << co2_average << "   co2 = " << co2.x[ i ][ j ][ k ] * co2_0 << "   epsilon_3D = " << epsilon_3D.x[ i ][ j ][ k ] << "   epsilon_eff = " << epsilon_eff << "   co2max = " << co2max * co2_0 << "   x_c = " << x_c << "   u_c = " << u_c << "   eps_co2 = " << eps_co2 << "   p_stat_co2 = " << 1.e-6 * p_stat.x[ i ][ j ][ k ] * co2.x[ i ][ j ][ k ] * co2_0 << "   p_stat = " << p_stat.x[ i ][ j ][ k ] << "   step = " << step << "   PcL = " << P_c * step * 100. << endl;
*/
            }

    cout.precision ( 8 );
    cout.setf ( ios::fixed );
/*
    if ( k == 180 ) cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_average = " << t_average << "   t = " << t.x[ i_mount ][ j ][ k ] * t_0 - t_0 << "   co2_coeff = " << co2_coeff << "   co2_eff = " << ( co2_pole - co2_equator ) << "   co2_pole = " << co2_pole << "   co2_equator = " << co2_equator << "   co2_paleo = " << co2_paleo << "   co2_average = " << co2_average << "   co2 = " << co2.x[ i_mount ][ j ][ k ] * co2_0 << "   epsilon_3D = " << epsilon_3D.x[ i_mount ][ j ][ k ] << "   epsilon_eff = " << epsilon_eff << "   co2max = " << co2max * co2_0 << "   x_c = " << x_c << "   u_c = " << u_c << "   eps_co2 = " << eps_co2 << "   p_stat_co2 = " << 1.e-6 * p_stat.x[ i_mount ][ j ][ k ] * co2.x[ i ][ j ][ k ] * co2_0 << "   p_stat = " << p_stat.x[ i_mount ][ j ][ k ] << "   step = " << step << "   PcL = " << P_c * step * 100. << endl;
*/
/*
    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_average = " << t_average << "   t = " << t.x[ i_mount ][ j ][ k ] * t_0 - t_0 << "   co2_coeff = " << co2_coeff << "   co2_eff = " << ( co2_pole - co2_equator ) << "   co2_pole = " << co2_pole << "   co2_equator = " << co2_equator << "   co2_paleo = " << co2_paleo << "   co2_average = " << co2_average << "   co2 = " << co2.x[ i_mount ][ j ][ k ] * co2_0 << "   epsilon_3D = " << epsilon_3D.x[ i_mount ][ j ][ k ] << "   epsilon_eff = " << epsilon_eff << "   co2max = " << co2max * co2_0 << "   x_c = " << x_c << "   u_c = " << u_c << "   eps_co2 = " << eps_co2 << "   p_stat_co2 = " << 1.e-6 * p_stat.x[ i_mount ][ j ][ k ] * co2.x[ i_mount ][ j ][ k ] * co2_0 << "   p_stat = " << p_stat.x[ i_mount ][ j ][ k ] << "   step = " << step << "   PcL = " << P_c * step * 100. << endl;
*/
/*
            // inside mountains
            for ( int i = i_mount - 1; i >= 0; i-- ){
                epsilon_3D.x[ i ][ j ][ k ] = epsilon_3D.x[ i_mount ][ j ][ k ];
                radiation_3D.x[ i ][ j ][ k ] = radiation_3D.x[ i_mount ][ j ][ k ];
            }
*/
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
            i_trop = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );

            for ( int k = 0; k < km; k++ ){
//                i_mount = i_topography[ j ][ k ];
                i_mount = 0;

                std::vector<double> alfa(im, 0);
                std::vector<double> beta(im, 0);
                std::vector<double> AA(im, 0);
                std::vector<std::vector<double> > CC(im, std::vector<double>(im, 0));
                double CCC = 0, DDD = 0;

                // radiation leaving the atmosphere above the tropopause, later needed for non-dimensionalisation
                radiation_3D.x[ i_trop ][ j ][ k ] = ( 1. - epsilon_3D.x[ i_trop ][ j ][ k ] ) * sigma * 
                    pow ( t.x[ i_trop ][ j ][ k ] * t_0, 4. ); 

                // back radiation absorbed by the first water vapour layer out of 40
                radiation_back = epsilon_3D.x[ i_mount + 1 ][ j ][ k ] * sigma * 
                    pow ( t.x[ i_mount + 1 ][ j ][ k ] * t_0, 4. );
                atmospheric_window = .1007 * radiation_surface.y[ j ][ k ]; // radiation loss through the atmospheric window
                rad_surf_diff = radiation_back + radiation_surface.y[ j ][ k ] - atmospheric_window; // radiation leaving the surface

//                fac_rad = ( double ) i_mount * .07 + 1.;  // linear increase with height, best choice for Ma>0
                fac_rad = 1.;  // linear increase with height, best choice for Ma>0
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
            }  // k-loop
        }  // j-loop
    }  // iter_rad-loop
}



void BC_Thermo::BC_Temperature( Array_1D &rad, Array_2D &temperature_NASA, Array &h,
                             Array &t, Array &tn, Array &p_dyn, Array &p_stat ){
//    logger() << std::endl << "enter BC_Temperature: temperature max: " << (t.max()-1)*t_0 << std::endl;
//    logger() << "enter BC_Temperature: temperature min: " << (t.min()-1)*t_0 << std::endl << std::endl;
    double t_paleo_add = 0;
// Lenton_etal_COPSE_time_temp, constant paleo mean temperature, added to the surface initial temperature
// difference between mean temperature ( Ma ) and mean temperature ( previous Ma ) == t_paleo_add
    if(!m_model->is_first_time_slice()){
        t_paleo_add = m_model->get_mean_temperature_from_curve(Ma) -  
            m_model->get_mean_temperature_from_curve(*m_model->get_previous_time());
        t_paleo_add /= t_0;  // non-dimensional
    }

//    i_mount = i_topography[ j ][ k ];
    i_mount = 0;
    t_pal[ i_mount ] = t_paleo_add;
//    t_pal[ i_mount ] = t_paleo;
    t_paleo_add = t_paleo;
    i_trop = im_tropopause[ j ];
    height_tropo = ( exp( zeta * ( rad.z[ i_trop ] - 1. ) ) - 1 ) 
        * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
    d_i_max = height_tropo;
    for ( int i = i_mount; i < im; i++ ){
        if ( i <= i_trop ){
            height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) 
                * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
            d_i = height;
            t_pal[ i ] = - t_paleo_add * ( d_i / d_i_max - 1. );  // linear t_paleo_add decay up to tropopause
//            t_pal[ i ] = - t_paleo * ( d_i / d_i_max - 1. );  // linear t_paleo_add decay up to tropopause

//    cout.precision ( 8 );
//    cout.setf ( ios::fixed );
//    cout << "   i_mount = " << i_mount << "   i = " << i << "   t_paleo = " << t_paleo * t_0 << "   t_paleo_add = " << t_paleo_add * t_0 << "   t_pal = " << t_pal[ i ] * t_0 << endl;

        }else{ // above the tropopause
            t_pal[ i ] = 0.;
        }
    }


    cout.precision ( 3 );

    time_slice_comment = "      time slice of Paleo-AGCM:";
    time_slice_number = " Ma = ";
    time_slice_unit = " million years";

    cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) << time_slice_comment << 
        resetiosflags ( ios::left ) << setw ( 6 ) << fixed << setfill ( ' ' ) << time_slice_number << 
        setw ( 3 ) << Ma << setw ( 12 ) << time_slice_unit << endl << endl;

    temperature_comment = "      temperature increase at paleo times: ";
    temperature_gain = " t increase";
    temperature_modern = "      mean temperature at modern times: ";
    temperature_paleo = "      mean temperature at paleo times: ";
    temperature_average = " t modern";
    temperature_average_pal = " t paleo";
    temperature_unit =  "°C ";

    cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) << temperature_comment << 
        resetiosflags ( ios::left ) << setw ( 12 ) << temperature_gain << " = " << setw ( 7 ) << setfill ( ' ' ) << 
        t_paleo << setw ( 5 ) << temperature_unit << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) 
        << temperature_modern << resetiosflags ( ios::left ) << setw ( 13 ) << temperature_average  << " = "  << setw ( 7 )  
        << setfill ( ' ' ) << t_average << setw ( 5 ) << temperature_unit << endl << setw ( 55 ) << setfill ( '.' )  << 
        setiosflags ( ios::left ) << temperature_paleo << resetiosflags ( ios::left ) << setw ( 13 ) << 
        temperature_average_pal  << " = "  << setw ( 7 )  << setfill ( ' ' ) << t_average + t_paleo << setw ( 5 ) << 
        temperature_unit << endl;

    // temperatur distribution at a prescribed sun position
    // sun_position_lat = 60,    position of sun j = 120 means 30°S, j = 60 means 30°N
    // sun_position_lon = 180, position of sun k = 180 means 0° or 180° E ( Greenwich, zero meridian )
    // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature ( linear equation + parabola )

    if ( ( Ma > 0 ) && ( sun == 1 ) ){
        j_par = sun_position_lat; // position of maximum temperature, sun position
        j_par = j_par + declination; // angle of sun axis, declination = 23,4°
        j_pol = jm - 1;
        j_par_f = ( double ) j_par;
        j_pol_f = ( double ) j_pol;

        aa = ( t_equator - t_pole ) / ( ( ( j_par_f * j_par_f ) - ( j_pol_f * j_pol_f ) ) - 2. * j_par_f * 
            ( j_par_f - j_pol_f ) );
        bb = - 2. * aa * j_par_f;
        cc = t_equator + aa * j_par_f * j_par_f;
        j_d = sqrt ( ( cc - t_pole ) / aa );
        dd = 2. * aa * j_d + bb;
        t_dd = dd * j_d + t_pole;
        e = t_pole;

        // asymmetric temperature distribution from pole to pole for  j_d  maximum temperature ( linear equation + parabola )
        for ( int k = 0; k < km; k++ ){
            for ( int j = 0; j < jm; j++ ){
                d_j = ( double ) j;
                if ( d_j <= j_d ){
                    t.x[ 0 ][ j ][ k ] = dd * d_j + e + t_paleo_add;
                }
                if ( d_j > j_d ){
                    t.x[ 0 ][ j ][ k ] = aa * d_j * d_j + bb * d_j + cc + t_paleo_add;
                }
            }
        }

        // longitudinally variable temperature distribution from west to east in parabolic form
        // communicates the impression of local sun radiation on the southern hemisphere
        k_par = sun_position_lon;  // position of the sun at constant longitude
        k_pol = km - 1;

        double t_360 = (  t_0 + 5. ) / t_0;

        for ( int j = 0; j < jm; j++ ){
            double jm_temp_asym = t.x[ 0 ][ j ][ 20 ];//transfer of zonal constant temperature into aa 1D-temperature field
            for ( int k = 0; k < km; k++ ){
                k_par_f = ( double ) k_par;
                k_pol_f = ( double ) k_pol;
                d_k = ( double ) k;

                aa = ( jm_temp_asym - t_360 ) / ( ( ( k_par_f * k_par_f ) - ( k_pol_f * k_pol_f ) ) - 2. * k_par_f * 
                    ( k_par_f - k_pol_f ) );
                bb = - 2. * aa * k_par_f;
                cc = jm_temp_asym + aa * k_par_f * k_par_f;

                t.x[ 0 ][ j ][ k ] = aa * d_k * d_k + bb * d_k + cc;
            }
        }
    }// temperatur distribution at aa prescribed sun position


// pole temperature adjustment, combination of linear time dependent functions 
    double t_pole_ma = 0.;  // Stein/Rüdiger/Parish locally constant pole temperature
                            // difference between pole temperature ( Ma ) and pole temperature ( previous Ma ) == t_pole_ma

    std::map<int, double> pole_temp_map{  // Stein/Rüdiger/Parish linear pole temperature ( Ma ) distribution

//        {0, 0.},
        {0, -15.4},
        {40, 22.},
        {45, 23.5},
        {50, 24.1},
        {55, 24.3},
        {60, 22.4},
        {70, 24.2},
        {80, 23.7},
        {90, 22.8},
        {100, 21.8},
        {120, 19.},
        {130, 17.8},
        {140, 16.9},
        {150, 16.4},
        {160, 16.},
        {340, 16.}
    }; 

    double d_j_half = ( double ) ( jm-1 ) / 2.;
    // temperature initial conditions along the surface
    if ( RadiationModel == 1 ){  // Multi-Layer-Radiation Model is active
        t_pole_ma = ( GetPoleTemperature ( Ma, pole_temp_map ) + t_0 ) / t_0;
    // non-dimensional, constant local pole temperature as function of Ma for hothouse climates, Stein/Rüdiger/Parish pole temperature 
        t_eff = t_pole - t_equator;
          // non-dimensional, coefficient for the zonal parabolic temperature distribution
        for ( int k = 0; k < km; k++ ){
            for ( int j = 1; j < jm-1; j++ ){
//                i_mount = i_topography[ j ][ k ];  // data along the topography
                i_mount = 0;  // data along the suface
                d_j = ( double ) j;
                if ( NASATemperature == 0 ){  // parabolic ocean surface temperature assumed
                    t.x[ i_mount ][ j ][ k ] = t_eff * parabola( d_j / d_j_half ) + t_pole + t_paleo_add;
                              //  constant pole and increasing mean temperature( Ma ) incorporated
                    if ( is_land ( h, 0, j, k ) ){  // parabolic land surface temperature assumed
                        t.x[ i_mount ][ j ][ k ] = t_eff * parabola( d_j / d_j_half ) + t_pole
                            + t_paleo_add + t_land;
                              //  constant pole and increasing mean temperature( Ma ) incorporated
                              // t_land, if mean land temperature is higher
                    }
                }else{  // for ( NASATemperature == 1 ) surface temperature by NASA applied
                        // transported for later time slices Ma by use_earthbyte_reconstruction
                    if ( is_land (h, 0, j, k ) ){  // on land a parabolic distribution is assumed, no NASA based data transportable
                        t.x[ i_mount ][ j ][ k ] = t_eff * parabola( d_j / d_j_half ) + t_pole
                            + t_paleo_add + t_land;
                              //  constant pole and increasing mean temperature( Ma ) incorporated
                              // t_land, if mean land temperature is higher
                        if ( t.x[ i_mount ][ 0 ][ k ] < t_pole_ma )
                            t.x[ i_mount ][ j ][ k ] = t.x[ i_mount ][ j ][ k ]
                                + ( t_pole_ma - t.x[ i_mount ][ 0 ][ k ] ) * fabs( parabola( d_j / d_j_half ) + 1. );
                                  // Stein/Rüdiger/Parish pole temperature decreasing parabolically equator wards
                        if ( Ma == 0 )  t.x[ i_mount ][ j ][ k ] = temperature_NASA.y[ j ][ k ];  // initial temperature by NASA for Ma=0
                    }
                    if ( is_air ( h, 0, j, k ) ){  // NASA based ocean surface temperature by use_earthbyte_reconstruction
                        if( Ma == 0 ){
                            t.x[ i_mount ][ j ][ k ] = temperature_NASA.y[ j ][ k ];  // initial temperature by NASA for Ma=0
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_pole_ma = " << t_pole_ma * t_0 - t_0 << "   t_paleo = " << t_paleo * t_0 << "   t_paleo_add = " << t_paleo_add * t_0 << "   fabs_parabola = " << parabola( d_j / d_j_half ) + 1. << "   t_pole_ma-fabs_parabola = " << ( t_pole_ma - temperature_NASA.y[ 0 ][ k ] ) * fabs( parabola( d_j / d_j_half ) + 1. ) * t_0 << "   temperature_NASA = " << temperature_NASA.y[ 0 ][ k ] * t_0 - t_0 << "   t = " << t.x[ i_mount ][ 0 ][ k ] * t_0 - t_0 << endl;
*/
                        }else{
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    if ( ( j == 0 ) && ( k == 180 ) ) cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_pole_ma = " << t_pole_ma * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_paleo_add = " << t_paleo_add * t_0 << "   fabs_parabola = " << parabola( d_j / d_j_half ) + 1. << "   t_pole_ma-fabs_parabola = " << ( t_pole_ma - t.x[ i_mount ][ 0 ][ k ] ) * fabs( parabola( d_j / d_j_half ) + 1. ) * t_0 << "   temperature_NASA = " << temperature_NASA.y[ 0 ][ k ] * t_0 - t_0 << "   t = " << t.x[ i_mount ][ 0 ][ k ] * t_0 - t_0 << endl;
*/
                            t.x[ i_mount ][ j ][ k ] = t.x[ i_mount ][ j ][ k ] + t_paleo_add;
                            if ( t.x[ i_mount ][ 0 ][ k ] < t_pole_ma )
                                t.x[ i_mount ][ j ][ k ] = t.x[ i_mount ][ j ][ k ]
                                + ( t_pole_ma - t.x[ i_mount ][ 0 ][ k ] ) * fabs( parabola( d_j / d_j_half ) + 1. );
                                  // Stein/Rüdiger/Parish pole temperature decreasing parabolically equator wards
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    if ( ( j == 0 ) && ( k == 180 ) ) cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_pole_ma = " << t_pole_ma * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_paleo_add = " << t_paleo_add * t_0 << "   fabs_parabola = " << parabola( d_j / d_j_half ) + 1. << "   t_pole_ma-fabs_parabola = " << ( t_pole_ma - t.x[ i_mount ][ 0 ][ k ] ) * fabs( parabola( d_j / d_j_half ) + 1. ) * t_0 << "   temperature_NASA = " << temperature_NASA.y[ 0 ][ k ] * t_0 - t_0 << "   t = " << t.x[ i_mount ][ 0 ][ k ] * t_0 - t_0 << endl;
*/
                        }// Ma > 0
                    }// is_air
                }// else ( NASATemperature == 1 )
            }// for j
        }// for k
    }// if ( RadiationModel == 1 )

    // zonal temperature along tropopause
    t_tropopause_pole = - 4.; // temperature reduction at poles in°C
    t_tropopause_pole = t_tropopause + t_tropopause_pole / t_0;
    t_eff_tropo = t_tropopause_pole - t_tropopause;
    double temp_tropopause = 0.;
    // temperature approaching the tropopause given by the Standard Atmosphere, above temperature is constant
    for ( int j = 0; j < jm; j++ ){
//        i_trop = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo_add );
        i_trop = im_tropopause[ j ] + GetTropopauseHightAdd ( t_pal[ i_trop ] );
        d_j = ( double ) j;
        height_tropo = ( exp( zeta * ( rad.z[ i_trop ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
        d_i_max = height_tropo;
//        double temp_tropopause = t_eff_tropo * parabola( d_j / d_j_half ) +
//                t_tropopause_pole + t_paleo_add;        

        for ( int k = 0; k < km; k++ ){
//            i_mount = i_topography[ j ][ k ];
            i_mount = 0;
            for ( int i = 0; i < im; i++ ){
                if ( i <= i_trop ){
                    temp_tropopause = t_eff_tropo * parabola( d_j / d_j_half ) +
                        t_tropopause_pole + t_pal[ i ];        
                    height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                    d_i = height;
                    t.x[ i ][ j ][ k ] = ( temp_tropopause - t.x[ i_mount ][ j ][ k ] ) * ( d_i / d_i_max ) + 
                        t.x[ i_mount ][ j ][ k ];  // linear temperature decay up to tropopause, privat  approximation
                }else{ // above the tropopause
                    t.x[ i ][ j ][ k ] = temp_tropopause;
                }
            }
//            for ( int i = 0; i < i_trop; i++ ){
//                if ( is_land ( h, i, j, k ) )  t.x[ i ][ j ][ k ] = 1.; // inside mountains
//                if ( is_land ( h, i, j, k ) )  t.x[ i ][ j ][ k ] = t.x[ i_mount ][ j ][ k ]; // inside mountains
//            }
        }
    }
}









void BC_Thermo::BC_WaterVapour ( Array_1D &rad, Array &h, Array &p_stat, Array &t, Array &c,
                            Array &v, Array &w ){
    // initial and boundary conditions of water vapour on water and land surfaces
    // parabolic water vapour distribution from pole to pole accepted

    // maximum water vapour content on water surface at equator c_equator = 1.04 compares to 0.04 volume parts
    // minimum water vapour at tropopause c_tropopause = 0.0 compares to 0.0 volume parts
    // value 0.04 stands for the maximum value of 40 g/kg, g water vapour per kg dry air
    // water vapour contents computed by Clausius-Clapeyron-formula
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            i_mount = i_topography[ j ][ k ];

            if ( is_air ( h, 0, j, k ) ){
                i_mount = 0;
                c.x[ i_mount ][ j ][ k ] = hp * ep * exp ( 17.0809 * ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + 
                    ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) ) ) / ( ( r_air * R_Air * t.x[ i_mount ][ j ][ k ] * t_0 ) * .01 );
                // saturation of relative water vapour in kg/kg
                c.x[ i_mount ][ j ][ k ] = c_ocean * c.x[ i_mount ][ j ][ k ];
                // relativ water vapour contents on ocean surface reduced by factor
            }
            if ( is_land ( h, 0, j, k ) ){
                i_mount = 0;
                c.x[ i_mount ][ j ][ k ] = hp * ep * exp ( 17.0809 * ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) / ( 234.175 + 
                    ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) ) ) / ( ( r_air * R_Air * t.x[ i_mount ][ j ][ k ] * t_0 ) * .01 );
                c.x[ i_mount ][ j ][ k ] = c_land * c.x[ i_mount ][ j ][ k ];
                // relativ water vapour contents on land reduced by factor
            }
        }
    }


    // water vapour distribution decreasing approaching tropopause
    for ( int j = 0; j < jm; j++ ){
//        i_trop = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
//        height_tropo = ( exp( zeta * ( rad.z[ i_trop ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
//        d_i_max = height_tropo;
        for ( int k = 0; k < km; k++ ){
            for ( int i = 1; i < im; i++ ){  // if temperature is below -37°C, there will exist c_tropopause
                if ( ( t.x[ i ][ j ][ k ] * t_0 - t_0 ) <= -37. ){
                i_trop = i;
                height_tropo = ( exp( zeta * ( rad.z[ i_trop ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i_max = height_tropo;
                break;
                }
            }
//            i_mount = i_topography[ j ][ k ];
            i_mount = 0;
            for ( int i = 1; i < im; i++ ){
                if ( i < i_trop ){
                    height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                    d_i = height;
                    c.x[ i ][ j ][ k ] = c.x[ i_mount ][ j ][ k ] - ( c_tropopause - c.x[ i_mount ][ j ][ k ] ) * 
                        ( d_i / d_i_max * ( d_i / d_i_max - 2. ) );  // radial parabolic decrease
                }else{
                    c.x[ i ][ j ][ k ] = c_tropopause;
                }
//                if ( is_land ( h, i, j, k ) ){
//                    c.x[ i ][ j ][ k ] = c.x[ i_mount ][ j ][ k ];
//                }
            } // end i
        }// end k
    }// end j
}






void BC_Thermo::BC_CO2( int Ma, double L_atm, Array_1D &rad, Array_2D &Vegetation, Array &h, Array &t, Array &p_dyn, Array &co2 ){
    // initial and boundary conditions of CO2 content on water and land surfaces
    // parabolic CO2 content distribution from pole to pole accepted

    j_half = j_max / 2;

    // CO2-distribution by Ruddiman approximated by a parabola
    co2_paleo = 3.2886 * pow ( ( t_paleo + t_average ), 2 ) - 32.8859 *
        ( t_paleo + t_average ) + 102.2148;  // in ppm
    co2_average = 3.2886 * pow ( t_average, 2 ) - 32.8859 * t_average + 102.2148;  // in ppm
    co2_paleo = co2_paleo - co2_average;

    cout.precision ( 3 );

    co_comment = "      co2 increase at paleo times: ";
    co_gain = " co2 increase";
    co_modern = "      mean co2 at modern times: ";
    co_paleo_str = "      mean co2 at paleo times: ";
    co_average_str = " co2 modern";
    co_average_pal = " co2 paleo";
    co_unit =  "ppm ";

    cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) <<
        co_comment << resetiosflags ( ios::left )         << setw ( 12 ) << co_gain << " = "
        << setw ( 7 ) << setfill ( ' ' ) << co2_paleo << setw ( 5 ) << co_unit << 
        endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) << co_modern
        << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_str  << " = "
        << setw ( 7 )  << setfill ( ' ' ) << co2_average << setw ( 5 ) << co_unit 
        << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left )
        << co_paleo_str << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_pal
        << " = "  << setw ( 7 )  << setfill ( ' ' ) << co2_average + co2_paleo
        << setw ( 5 ) << co_unit << endl;
    cout << endl;

    d_i_max = ( double ) i_max;
    d_j_half = ( double ) j_half;

    co2_equator = co2_equator / co2_0;
    co2_pole = co2_pole / co2_0;
    co2_paleo = co2_paleo / co2_0;
    co2_land = co2_land / co2_0;
    co2_ocean = co2_ocean / co2_0;
    co2_tropopause = co2_tropopause / co2_0;
    co2_eff = co2_pole - co2_equator;

    double emittancy_total = 0.423; // in W/m²
    double coeff_em = 5.6697e-8; // in W/(m² K)
    double delta_T = 0.02; // in K

    // CO2-content as initial solution
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
//            i_mount = i_topography[ j ][ k ];
            i_mount = 0;
            if ( is_air ( h, i_mount, j, k ) ){
                d_j = ( double ) j;
/*
                co2.x[ i_mount ][ j ][ k ] = co2_eff * ( d_j * d_j 
                    / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + 
                    co2_pole + co2_paleo + co2_ocean; // non-dimensional, parabolic assumption
*/
                co2.x[ i_mount ][ j ][ k ] = exp ( 4. * delta_T * coeff_em 
                     * pow( ( t.x[ i_mount ][ j ][ k ] * t_0 ), 3 ) / emittancy_total )
                     + co2_paleo + co2_ocean;
                     // reciprocal formula for the temperature increase by co2, 
                     // Temp_co2_add by Nasif Nahle Sabag in PostProcess_Atm.cpp

                // taken over from Ruddiman, p 86, effect of co2 on global temperature
//                co2.x[ i_mount ][ j ][ k ] = ( 3.2886 * pow ( ( t.x[ i_mount ][ j ][ k ] 
//                    * t_0 - t_0 ), 2 ) - 32.8859 * ( t.x[ i_mount ][ j ][ k ] 
//                    * t_0 - t_0 ) + 102.2148 + co2_ocean ) / co2_0;  // non-dimensional
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_average = " << t_average << "   t = " << t.x[ i_mount ][ 0 ][ k ] * t_0 - t_0 << "   co2_eff = " << ( co2_pole - co2_equator ) * co2_0 << "   co2_pole = " << co2_pole * co2_0 << "   co2_equator = " << co2_equator * co2_0 << "   co2_paleo = " << co2_paleo * co2_0 << "   co2_average = " << co2_average << "   co2 = " << co2.x[ i_mount ][ 0 ][ k ] * co2_0 << endl;
*/
            }
            if ( is_land ( h, i_mount, j, k ) ){
                d_j = ( double ) j;
/*
                co2.x[ i_mount ][ j ][ k ] = co2_eff * ( d_j * d_j 
                    / ( d_j_half * d_j_half ) - 2. * d_j / d_j_half ) + 
                    co2_pole + co2_paleo + co2_land + co2_vegetation 
                    * Vegetation.y[ j ][ k ] / co2_0;  // parabolic distribution from pole to pole
*/
                co2.x[ i_mount ][ j ][ k ] = exp ( 4. * delta_T * coeff_em 
                     * pow( ( t.x[ i_mount ][ j ][ k ] * t_0 ), 3 ) / emittancy_total )
                     + co2_paleo + co2_land + co2_vegetation
                    * Vegetation.y[ j ][ k ] / co2_0;
                     // reciprocal formula for the temperature increase by co2, 
                     // Temp_co2_add by Nasif Nahle Sabag in PostProcess_Atm.cpp




                // taken over from Ruddiman, p 86, effect of co2 on global temperature
//                co2.x[ i_mount ][ j ][ k ] = ( 3.2886 * pow ( ( t.x[ i_mount ][ j ][ k ] 
//                    * t_0 - t_0 ), 2 ) - 32.8859 * ( t.x[ i_mount ][ j ][ k ] * t_0 - t_0 ) 
//                    + 102.2148 + co2_land - co2_vegetation * Vegetation.y[ j ][ k ] ) / co2_0;  // non-dimensional
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    if ( ( j == 90 ) && ( k == 30 ) ) cout << "   i_mount = " << i_mount << "   j = " << j << "   k = " << k << "   t_eff = " << ( t_pole - t_equator ) * t_0 << "   t_pole = " << t_pole * t_0 - t_0 << "   t_equator = " << t_equator * t_0 - t_0 << "   t_paleo = " << t_paleo << "   t_average = " << t_average << "   t = " << t.x[ i_mount ][ 0 ][ k ] * t_0 - t_0 << "   co2_eff = " << ( co2_pole - co2_equator ) * co2_0 << "   co2_pole = " << co2_pole * co2_0 << "   co2_equator = " << co2_equator * co2_0 << "   co2_paleo = " << co2_paleo * co2_0 << "   co2_average = " << co2_average << "   co2 = " << co2.x[ i_mount ][ 0 ][ k ] * co2_0 << "   co2_vegetation = " << co2_vegetation << "   Vegetation = " << Vegetation.y[ j ][ k ] << endl;
*/
            }
        }
    }

    // co2 distribution decreasing approaching tropopause, above no co2
    for ( int j = 0; j < jm; j++ ){
        i_trop = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
//        d_i_max = ( double ) i_trop;
        height_tropo = ( exp( zeta * ( rad.z[ i_trop ] - 1. ) ) - 1 ) 
            * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
        d_i_max = height_tropo;

        for ( int k = 0; k < km; k++ ){
//            i_mount = i_topography[ j ][ k ];
            i_mount = 0;
            for ( int i = 0; i <= im - 1; i++ ){
                if ( i <= i_trop ){
                    height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) 
                        * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                    d_i = height;
                    co2.x[ i ][ j ][ k ] = co2.x[ i_mount ][ j ][ k ] 
                        - ( co2_tropopause - co2.x[ i_mount ][ j ][ k ] ) 
                        * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) );
                        // radial distribution approximated by a parabola
                }
                else  co2.x[ i ][ j ][ k ] = co2_tropopause;
            }
        }
    }
}





void BC_Thermo::BC_CO2_Iter( int Ma, double L_atm, Array_1D &rad, Array_2D &Vegetation,
                Array_2D &Topography, Array &h, Array &t, Array &p_dyn, Array &co2 ){
    // initial and boundary conditions of CO2 content on water and land surfaces
    // parabolic CO2 content distribution from pole to pole accepted
    j_half = j_max / 2;

    // CO2-distribution by Ruddiman approximated by a parabola
    co2_paleo = 3.2886 * pow ( ( t_paleo + t_average ), 2 ) - 32.8859 *
        ( t_paleo + t_average ) + 102.2148;  // in ppm
    co2_average = 3.2886 * pow ( t_average, 2 ) - 32.8859 * t_average + 102.2148;  // in ppm
    co2_paleo = co2_paleo - co2_average;

    cout.precision ( 3 );

    co_comment = "      co2 increase at paleo times: ";
    co_gain = " co2 increase";
    co_modern = "      mean co2 at modern times: ";
    co_paleo_str = "      mean co2 at paleo times: ";
    co_average_str = " co2 modern";
    co_average_pal = " co2 paleo";
    co_unit =  "ppm ";

    cout << endl << setiosflags ( ios::left ) << setw ( 55 ) << setfill ( '.' ) <<
        co_comment << resetiosflags ( ios::left )         << setw ( 12 ) << co_gain << " = "
        << setw ( 7 ) << setfill ( ' ' ) << co2_paleo << setw ( 5 ) << co_unit << 
        endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left ) << co_modern
        << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_str  << " = "
        << setw ( 7 )  << setfill ( ' ' ) << co2_average << setw ( 5 ) << co_unit 
        << endl << setw ( 55 ) << setfill ( '.' )  << setiosflags ( ios::left )
        << co_paleo_str << resetiosflags ( ios::left ) << setw ( 13 ) << co_average_pal
        << " = "  << setw ( 7 )  << setfill ( ' ' ) << co2_average + co2_paleo
        << setw ( 5 ) << co_unit << endl;
    cout << endl;

    d_i_max = ( double ) i_max;
    d_j_half = ( double ) j_half;

    co2_equator = co2_equator / co2_0;
    co2_pole = co2_pole / co2_0;
    co2_paleo = co2_paleo / co2_0;
    co2_land = co2_land / co2_0;
    co2_ocean = co2_ocean / co2_0;
    co2_tropopause = co2_tropopause / co2_0;

    co2_eff = co2_pole - co2_equator;

    // CO2-content as initial solution
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
//            i_mount = i_topography[ j ][ k ];
            i_mount = 0;
            if ( ( Ma > 0 ) && ( is_land ( h, i_mount, j, k ) ) ) 
                co2.x[ i_mount ][ j ][ k ] = co2.x[ i_mount ][ j ][ k ]
                    + co2_vegetation * Vegetation.y[ j ][ k ] / co2_0; 
        }
    }

    // co2 distribution decreasing approaching tropopause, above no co2
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            height_tropo = Topography.y[ j ][ k ];
            d_i_max = height_tropo;

            for ( int i = 0; i <= im - 1; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) 
                    * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                if ( height <= height_tropo ) i_trop = i;
            }

//            i_mount = i_topography[ j ][ k ];
            i_mount = 0;
            if ( Ma > 0 ){
                for ( int i = 0; i <= im - 1; i++ ){
                    if ( is_land ( h, i, j, k ) ){
                        height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) 
                            * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                        d_i = height;
                        co2.x[ i ][ j ][ k ] = co2.x[ i_mount ][ j ][ k ] 
                            - ( co2.x[ i_trop ][ j ][ k ] - co2.x[ i_mount ][ j ][ k ] ) 
                            * ( d_i / d_i_max * ( d_i / d_i_max - 2. ) );
                            // radial distribution approximated by a parabola
/*
    cout.precision ( 2 );
    cout.setf ( ios::fixed );
    if ( ( j == 60 ) && ( k == 87 ) ) cout << "   i = " << i << "   j = " << j << "   k = " << k << "   i_trop = " << i_trop << "   Topography = " << Topography.y[ j ][ k ] << "   height_tropo = " << height_tropo << "   height = " << height << "   t = " << t.x[ i ][ 0 ][ k ] * t_0 - t_0 << "   co2_paleo = " << co2_paleo * co2_0 << "   co2_tropopause = " << co2_tropopause << "   co2_mount = " << co2.x[ i_mount ][ j ][ k ] * co2_0 << "   co2 = " << co2.x[ i ][ j ][ k ] * co2_0 << endl;
*/
                    }
                }
            }
        }
    }
}





void BC_Thermo::TropopauseLocation(){
// parabolic tropopause location distribution from pole to pole assumed
    j_max = jm - 1;
    j_half = ( jm -1 ) / 2;
    j_infl = 43;  // fixed value due to the cubic function for the location of the tropopause


    // flattens the equator peak
    d_j_half = ( double ) j_half;
    d_j_infl = 3. * ( double ) j_infl;
    trop_pole = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
        * log( ( double ) ( im - 1 ) * tropopause_pole / L_atm  + 1. ) + 1. ) );  // coordinate stretching
    trop_equator = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
        * log( ( double ) ( im - 1 ) * tropopause_equator / L_atm  + 1. ) + 1. ) );  // coordinate stretching
    trop_co2_eff = ( double ) ( trop_pole - trop_equator );


// no stripes in longitudinal direction
// computation of the tropopause from pole to pole, constant approach
    for ( int j = 0; j < jm; j++ ){ // parabolic approach
        im_tropopause[ j ] = ( int ) trop_equator; // constant approach
    }

/*
// minor stripes in longitudinal direction
// computation of the tropopause from pole to pole, parabolic approach
    for ( int j = 0; j < jm; j++ ) // parabolic approach
    {
        d_j = ( double ) j;
        im_tropopause[ j ] = ( trop_co2_eff * ( d_j * d_j / ( d_j_half * d_j_half ) -
            2. * d_j / d_j_half ) ) + trop_pole; // parabolic approach
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
            ( double ) trop_pole );  // cubic approach

        im_tropopause[ j ] = ( int ) trop;  // cubic approach
    }

    for ( int j = j_half + 1; j < jm; j++ )  // cubic approach
    {
        im_tropopause[ j ] = im_tropopause[ j_max - j ];  // cubic approach
    }
*/
}

void BC_Thermo::IC_CellStructure ( Array_1D &rad, Array &h, Array &u, Array &v, Array &w ){
// boundary condition for the velocity components in the circulation cells

// latest version by Grotjahn ( Global Atmospheric Circulations, 1993 )
// default for the velocity components u, v, and w as initial conditions

// velocities given in m/s, 1 m/s compares to 3.6 km/h, non-dimensionalized by u_0 at the end of this class element
// do not change the velocity initial conditions !!

// velocity assumptions at the equator 0°
    ua_00 = 1.;  // in m/s compares to 3.6 km/h, non-dimensionalized by u_0 at the end of this class elemen
//    va_equator_SL =  0.000;
    va_equator_Tropopause = 0.000;
//    wa_equator_SL = - 1.;
    wa_equator_Tropopause = - 7.5;
    double va_equator_SL_land =  0.000;
    double wa_equator_SL_land = - 1.;

// velocity assumptions for latitude at 15° and 30° in the Hadley cell
    ua_30 = + 1.;
//    va_Hadley_SL = .25;
    va_Hadley_Tropopause = - 1.;
//    va_Hadley_SL_15 = 1.;
    va_Hadley_Tropopause_15 = - 1.;
//    wa_Hadley_SL = 1.;                                                // at surface
    wa_Hadley_Tropopause = 30.;                                         // subtropic jet in m/s compares to 108 km/h
    double va_Hadley_SL_land = .25;
    double va_Hadley_SL_15_land = 1.;
    double wa_Hadley_SL_land = 1.;                                                // at surface

// velocity assumptions for latitude at 45° and 60° in the Ferrel cell
    ua_60 = .5;
//    va_Ferrel_SL = 0.5;
    va_Ferrel_Tropopause = 1.;
//    va_Ferrel_SL_45 = -.1;
    va_Ferrel_Tropopause_45 = 1.;
//    wa_Ferrel_SL = -.2;                                               // subpolar jet
    wa_Ferrel_Tropopause = 10.;                                         // subpolar jet in m/s compares to 36 km/h
    double va_Ferrel_SL_land = 0.5;
    double va_Ferrel_SL_45_land = -.1;
    double wa_Ferrel_SL_land = -.2;                                               // subpolar jet

// velocity assumptions for latitude 90° in the Polar cell
    ua_90 = .5;
//    va_Polar_SL = 0.;
    va_Polar_Tropopause = 0.;
//    va_Polar_SL_75 = .5;
    va_Polar_Tropopause_75 = - 1.;
//    wa_Polar_SL = - 0.01;
    wa_Polar_Tropopause = 0.;
    double va_Polar_SL_land = 0.;
    double va_Polar_SL_75_land = .5;
    double wa_Polar_SL_land = - 0.01;

// preparations for diagonal velocity value connections
    im_1 = im - 1;
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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// initial velocity components in the northern and southern
// Pole, Ferrel and Hadley cells

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////// equator ///////////////////////////////////////
// equator ( at j=90 compares to 0° latitude )
// u-component up to tropopause and back on half distance (  i = 20 )
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_aeq; j < j_aeq + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            i_half = ( im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 ) ) / 2;
            height_tropo = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
                * log( ( double ) ( im - 1 ) * ( .5 * tropopause_equator ) 
                / L_atm  + 1. ) + 1. ) );  // coordinate stretching
            i_half = ( int ) height_tropo;
            height_tropo = ( exp( zeta * ( rad.z[ i_half ] - 1. ) ) - 1 ) 
                * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
            d_i_half = height_tropo;
            for ( int i = 0; i < i_half; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = ua_00 * ( d_i / d_i_half );
            }
            for ( int i = i_half; i < i_max; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = 1. - ( ua_00 * ( ( d_i - i_max ) / ( d_i_half - i_max ) ) - 1. );
            }
        }
    }

// equator ( at j=90 )
// v- and w-component up to tropopause and stratosphere above
    for ( int j = j_aeq; j < j_aeq + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) i_max;
        for ( int k = 0; k < km; k++ ){
            va_equator_SL = v.x[ 0 ][ j ][ k ];
            wa_equator_SL = w.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_equator_SL = va_equator_SL_land;
                wa_equator_SL = wa_equator_SL_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_equator_Tropopause - va_equator_SL ) *
                    d_i / d_i_max + va_equator_SL;
                w.x[ i ][ j ][ k ] = ( wa_equator_Tropopause - wa_equator_SL ) *
                    d_i / d_i_max + wa_equator_SL;
            }
        }
    }
    for ( int j = j_aeq; j < j_aeq + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) i_max - ( double ) im_1;
        if ( i_max == 40 ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) i - ( double ) im_1;
                v.x[ i ][ j ][ k ] = va_equator_Tropopause * d_i / d_i_max;
                w.x[ i ][ j ][ k ] = wa_equator_Tropopause * d_i / d_i_max;
            }
        }
    }

//    goto Printout;

/////////////////////////////////////// end equator ///////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////

// cell structure in northern hemisphere

////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////// northern polar cell /////////////////////////////////////////


// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// u-component up to tropopause and back on half distance
// extension around North Pole ( from j=0 till j=5 compares to 85° till 90° northern latitude )
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < j_pol_n +1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            i_half = ( im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 ) ) / 2;
            height_tropo = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
                * log( ( double ) ( im - 1 ) * ( .5 * tropopause_equator ) 
                / L_atm  + 1. ) + 1. ) );  // coordinate stretching
            i_half = ( int ) height_tropo;
            height_tropo = ( exp( zeta * ( rad.z[ i_half ] - 1. ) ) - 1 ) 
                * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
            d_i_half = height_tropo;

            for ( int i = 0; i < i_half; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = ua_90 * ( - ( d_i / d_i_half ) );
            }
            for ( int i = i_half; i < i_max; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = ua_90 * ( - ( 1. - ( ( ( d_i - i_max ) / ( d_i_half - i_max ) ) - 1. ) ) );
            }

        }
    }


// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// v- and w-component from Pole up to tropopause
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_n; j < j_fer_n + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            d_i_max = ( double ) i_max;
            va_Polar_SL = v.x[ 0 ][ j ][ k ];
            wa_Polar_SL = w.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Polar_SL = va_Polar_SL_land;
                wa_Polar_SL = wa_Polar_SL_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Polar_Tropopause - va_Polar_SL ) *
                    d_i / d_i_max + va_Polar_SL;
                w.x[ i ][ j ][ k ] = ( wa_Polar_Tropopause - wa_Polar_SL ) *
                    d_i / d_i_max + wa_Polar_SL;   // indifferent except at j_pol_n
            }
        }
    }
    for ( int j = j_pol_n; j < j_fer_n + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Polar_Tropopause * d_i / d_i_max;   // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = wa_Polar_Tropopause * d_i / d_i_max;   // replacement for forming diagonals
            }
        }
    }

// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// v-component up to tropopause and back on half distance
// extension around the North Pole ( from j=0 till j=5 compares to 85° till 90° northern latitude )
// v-component at 75°N
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_v_n; j <  j_pol_v_n + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            d_i_max = ( double ) i_max;
            va_Polar_SL_75 = v.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Polar_SL_75 = va_Polar_SL_75_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Polar_Tropopause_75 - va_Polar_SL_75 ) *
                    d_i / d_i_max + va_Polar_SL_75;
            }
        }
    }
    for ( int j = j_pol_v_n; j < j_pol_v_n + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Polar_Tropopause_75 * d_i / d_i_max;
            }
        }
    }

//    goto Printout;


/////////////////////////////////// end northern polar cell /////////////////////////////////////////

/////////////////////////////////// northern Ferrel cell /////////////////////////////////////////

// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° northern latitude )

// u-component up to tropopause and back on half distance
// extension around 60° northern latitude ( from j=29 till j=31 compares to 59° till 61° northern latitude )
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_n; j < j_fer_n + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            i_half = ( im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 ) ) / 2;
            height_tropo = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
                * log( ( double ) ( im - 1 ) * ( .5 * tropopause_equator ) 
                / L_atm  + 1. ) + 1. ) );  // coordinate stretching
            i_half = ( int ) height_tropo;
            height_tropo = ( exp( zeta * ( rad.z[ i_half ] - 1. ) ) - 1 ) 
                * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
            d_i_half = height_tropo;

            for ( int i = 0; i < i_half; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = ua_60 * ( d_i / d_i_half );
            }
            for ( int i = i_half; i < i_max; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = ua_60 * ( 1. - ( ( ( d_i - i_max ) / ( d_i_half - i_max ) ) - 1. ) );
            }

        }
    }

// north equatorial Ferrel cell
// v- and w-component up to tropopause and stratosphere above
// extension around 60° northern latitude ( from j=29 till j=31 compares to 59° till 61° northern latitude )
    for ( int j = j_fer_n; j < j_fer_n + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) i_max;
        for ( int k = 0; k < km; k++ ){
            va_Ferrel_SL = v.x[ 0 ][ j ][ k ];
            wa_Ferrel_SL = w.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Ferrel_SL = va_Ferrel_SL;
                wa_Ferrel_SL = wa_Ferrel_SL_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Ferrel_Tropopause - va_Ferrel_SL ) *
                    d_i / d_i_max + va_Ferrel_SL;   // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = ( wa_Ferrel_Tropopause - wa_Ferrel_SL ) *
                    d_i / d_i_max + wa_Ferrel_SL;   // replacement for forming diagonals
            }
        }
    }
    for ( int j = j_fer_n; j < j_fer_n + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Ferrel_Tropopause * d_i / d_i_max;   // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = wa_Ferrel_Tropopause * d_i / d_i_max;   // replacement for forming diagonals
            }
        }
    }

// north equatorial Ferrel cell
// v-component up to tropopause and stratosphere above
// extension around 60° northern latitude ( from j=29 till j=31 compares to 59° till 61° northern latitude )
// v-component at 45°N
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_v_n; j <  j_fer_v_n + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            d_i_max = ( double ) i_max;
            va_Ferrel_SL_45 = v.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Ferrel_SL_45 = va_Ferrel_SL_45_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Ferrel_Tropopause_45 - va_Ferrel_SL_45 ) *
                    d_i / d_i_max + va_Ferrel_SL_45;
            }
        }
    }
    for ( int j = j_fer_v_n; j < j_fer_v_n + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Ferrel_Tropopause_45 * d_i / d_i_max;
            }
        }
    }

//////////////////////////////// end northern Ferrel cell ////////////////////////////////////////

//////////////////////////////// northern Hadley cell ////////////////////////////////////////

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// u-component up to tropopause and back on half distance
// extension around 30° northern latitude ( from j=59 till j=61 compares to 29° till 31° northern latitude )
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_n; j < j_had_n + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            i_half = ( im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 ) ) / 2;
            height_tropo = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
                * log( ( double ) ( im - 1 ) * ( .5 * tropopause_equator ) 
                / L_atm  + 1. ) + 1. ) );  // coordinate stretching
            i_half = ( int ) height_tropo;
            height_tropo = ( exp( zeta * ( rad.z[ i_half ] - 1. ) ) - 1 ) 
                * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
            d_i_half = height_tropo;

            for ( int i = 0; i < i_half; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = - ua_30 * ( d_i / d_i_half );
            }
            for ( int i = i_half; i < i_max; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = - ( 1. - ( ua_30 * ( ( d_i - i_max ) / ( d_i_half - i_max ) ) - 1. ) );
            }

        }
    }

// north equatorial Hadley cell
// v- and w-component up to tropopause and stratosphere above
// extension around 30° northern latitude ( from j=59 till j=61 compares to 29° till 31° northern latitude )
    for ( int j = j_had_n; j < j_had_n + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) i_max;
        for ( int k = 0; k < km; k++ ){
            va_Hadley_SL = v.x[ 0 ][ j ][ k ];
            wa_Hadley_SL = w.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Hadley_SL = va_Hadley_SL;
                wa_Hadley_SL = wa_Hadley_SL_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Hadley_Tropopause - va_Hadley_SL ) *
                    d_i / d_i_max + va_Hadley_SL;    // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = ( wa_Hadley_Tropopause - wa_Hadley_SL ) *
                    d_i / d_i_max + wa_Hadley_SL;    // replacement for forming diagonals
            }
        }
    }
    for ( int j = j_had_n; j < j_had_n + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Hadley_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = wa_Hadley_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
            }
        }
    }

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// v-component up to tropopause and stratosphere above
// v-component at 15°N
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_v_n; j <  j_had_v_n + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            d_i_max = ( double ) i_max;
            va_Hadley_SL_15 = v.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Hadley_SL_15 = va_Hadley_SL_15_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Hadley_Tropopause_15 - va_Hadley_SL_15 ) *
                    d_i / d_i_max + va_Hadley_SL_15;
            }
        }
    }
    for ( int j = j_had_v_n; j < j_had_v_n + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Hadley_Tropopause_15 * d_i / d_i_max;
            }
        }
    }

//    goto Printout;

////////////////////////////////// end northern Hadley cell ////////////////////////////////////////

// ||||||||||||||||||||||||||||||||||||||| equator ||||||||||||||||||||||||||||||||||||||||||||||
// ||||||||||||||||||||||||||||||||||||||| equator ||||||||||||||||||||||||||||||||||||||||||||||

////////////////////////////////////////////////////////////////////////////////////////////////

// cell structure in southern hemisphere

////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////// southern Hadley cell ////////////////////////////////////////


// south equatorial Hadley cell ( from j=90 till j=120 compares to 0° till 30° southern latitude )
// u-component up to tropopause and back on half distance
// extension around 30° southern latitude ( from j=119 till j=121 compares to 29° till 31° southern latitude )
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_s; j < j_had_s + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            i_half = ( im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 ) ) / 2;
            height_tropo = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
                * log( ( double ) ( im - 1 ) * ( .5 * tropopause_equator ) 
                / L_atm  + 1. ) + 1. ) );  // coordinate stretching
            i_half = ( int ) height_tropo;
            height_tropo = ( exp( zeta * ( rad.z[ i_half ] - 1. ) ) - 1 ) 
                * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
            d_i_half = height_tropo;

            for ( int i = 0; i < i_half; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = - ua_30 * ( d_i / d_i_half );
            }
            for ( int i = i_half; i < i_max; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = - ( 1. - ( ua_30 * ( ( d_i - i_max ) / ( d_i_half - i_max ) ) - 1. ) );
            }

        }
    }

// south equatorial Hadley cell
// v- and w-component up to tropopause and stratosphere above
// extension around 30° southern latitude ( from j=119 till j=121 compares to 29° till 31° southern latitude )
    for ( int j = j_had_s; j < j_had_s + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) i_max;
        for ( int k = 0; k < km; k++ ){
            va_Hadley_SL = v.x[ 0 ][ j ][ k ];
            wa_Hadley_SL = w.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Hadley_SL = va_Hadley_SL_land;
                wa_Hadley_SL = wa_Hadley_SL_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Hadley_Tropopause - va_Hadley_SL ) *
                    d_i / d_i_max + va_Hadley_SL;    // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = ( wa_Hadley_Tropopause - wa_Hadley_SL ) *
                    d_i / d_i_max + wa_Hadley_SL;    // replacement for forming diagonals
            }
        }
    }
    for ( int j = j_had_s; j < j_had_s + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Hadley_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = wa_Hadley_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
            }
        }
    }

//     goto Printout;

// extension around 30° southern latitude ( from j=119 till j=121 compares to 29° till 31° southern latitude )
// v-component up to tropopause and stratosphere above
// v-component at 15°S
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_v_s; j <  j_had_v_s + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            d_i_max = ( double ) i_max;
            va_Hadley_SL_15 = v.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Hadley_SL_15 = va_Hadley_SL_15_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Hadley_Tropopause_15 - va_Hadley_SL_15 ) *
                    d_i / d_i_max + va_Hadley_SL_15;
            }
        }
    }
    for ( int j = j_had_v_s; j < j_had_v_s + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Hadley_Tropopause_15 * d_i / d_i_max;
            }
        }
    }

////////////////////////////// end southern Hadley cell ////////////////////////////////////////////

//    goto Printout;

////////////////////////////// southern Ferrel cell ////////////////////////////////////////////

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// u-component up to tropopause and back on half distance
// extension around 60° southern latitude ( from j=139 till j=151 compares to 59° till 61° southern latitude )
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_s; j < j_fer_s + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            i_half = ( im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 ) ) / 2;
            height_tropo = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
                * log( ( double ) ( im - 1 ) * ( .5 * tropopause_equator ) 
                / L_atm  + 1. ) + 1. ) );  // coordinate stretching
            i_half = ( int ) height_tropo;
            height_tropo = ( exp( zeta * ( rad.z[ i_half ] - 1. ) ) - 1 ) 
                * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
            d_i_half = height_tropo;

            for ( int i = 0; i < i_half; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = ua_60 * ( d_i / d_i_half );
            }
            for ( int i = i_half; i < i_max; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = ua_60 * ( 1. - ( ( ( d_i - i_max ) / ( d_i_half - i_max ) ) - 1. ) );
            }

        }
    }

// south equatorial Ferrel cell
// v- and w-component up to tropopause and stratosphere above
    for ( int j = j_fer_s; j < j_fer_s + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) i_max;
        for ( int k = 0; k < km; k++ ){
            va_Ferrel_SL = v.x[ 0 ][ j ][ k ];
            wa_Ferrel_SL = w.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Ferrel_SL = va_Ferrel_SL_land;
                wa_Ferrel_SL = wa_Ferrel_SL_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Ferrel_Tropopause - va_Ferrel_SL ) *
                    d_i / d_i_max + va_Ferrel_SL;    // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = ( wa_Ferrel_Tropopause - wa_Ferrel_SL ) *
                    d_i / d_i_max + wa_Ferrel_SL;    // replacement for forming diagonals
            }
        }
    }
    for ( int j = j_fer_s; j < j_fer_s + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Ferrel_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = wa_Ferrel_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
            }
        }
    }

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// v-component up to tropopause and stratosphere above
// v-component at 45°N
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_v_s; j <  j_fer_v_s + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            d_i_max = ( double ) i_max;
            va_Ferrel_SL_45 = v.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Ferrel_SL_45 = va_Ferrel_SL_45_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Ferrel_Tropopause_45 - va_Ferrel_SL_45 ) *
                    d_i / d_i_max + va_Ferrel_SL_45;
            }
        }
    }
    for ( int j = j_fer_v_s; j < j_fer_v_s + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Ferrel_Tropopause_45 * d_i / d_i_max;
            }
        }
    }

//    goto Printout;

///////////////////////////// end southern Ferrel cell /////////////////////////////////////////////

///////////////////////////// southern Polar cell /////////////////////////////////////////////

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// u-component up to tropopause and back on half distance
// extension around South Pole ( from j=175 till j=180 compares to 85° till 90° southern latitude )
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_s; j < jm; j++ ){
/*
            i_max = im_tropopause[ j ];
            i_half = im_tropopause[ j ] / 2;
            d_i_half = ( double ) i_half;
            for ( int i = 0; i <= i_max; i++ ){
                d_i = ( double ) i;
                u.x[ i ][ j ][ k ] = - ua_90 * ( d_i * d_i / ( d_i_half * d_i_half ) - 2. * d_i / d_i_half );
            }
*/

            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            i_half = ( im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 ) ) / 2;
            height_tropo = round ( ( double ) ( im - 1 ) / 2. * ( 1. / zeta 
                * log( ( double ) ( im - 1 ) * ( .5 * tropopause_equator ) 
                / L_atm  + 1. ) + 1. ) );  // coordinate stretching
            i_half = ( int ) height_tropo;
            height_tropo = ( exp( zeta * ( rad.z[ i_half ] - 1. ) ) - 1 ) 
                * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
            d_i_half = height_tropo;


            for ( int i = 0; i < i_half; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = - ua_90 * ( d_i / d_i_half );
//     cout << "    i = " << i << "    j = " << j << "    k = " << k << "    zeta = " << zeta << "    i_max = " << i_max << "    d_i_half = " << d_i_half << "    i_half = " << i_half << "    d_i = " << d_i << "    trop_equ = " << tropopause_equator << "    u = " << u.x[ i ][ j ][ k ] << endl;
            }
            for ( int i = i_half; i < i_max; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                d_i = height;
                u.x[ i ][ j ][ k ] = - ua_90 * ( 1. - ( ( ( d_i - i_max ) / ( d_i_half - i_max ) ) - 1. ) );
//     cout << "    i = " << i << "    j = " << j << "    k = " << k << "    zeta = " << zeta << "    i_max = " << i_max << "    d_i_half = " << d_i_half << "    i_half = " << i_half << "    d_i = " << d_i << "    trop_equ = " << tropopause_equator << "    u = " << u.x[ i ][ j ][ k ] << endl;
            }

        }
    }

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_s + 1; j < j_pol_s + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            d_i_max = ( double ) i_max;
            va_Polar_SL = v.x[ 0 ][ j ][ k ];
            wa_Polar_SL = w.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Polar_SL = va_Polar_SL_land;
                wa_Polar_SL = wa_Polar_SL_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( ( va_Polar_Tropopause - va_Polar_SL ) *
                    d_i / d_i_max + va_Polar_SL );
                w.x[ i ][ j ][ k ] = ( wa_Polar_Tropopause - wa_Polar_SL ) *
                    d_i / d_i_max + wa_Polar_SL;  // indifferent except at j_pol_s
            }
        }
    }
    for ( int j = j_fer_s + 1; j < j_pol_s + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Polar_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
                w.x[ i ][ j ][ k ] = wa_Polar_Tropopause * d_i / d_i_max;    // replacement for forming diagonals
            }
        }
    }

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// u-component up to tropopause and back on half distance
// extension around South Pole ( from j=175 till j=180 compares to 85° till 90° southern latitude )
// v-component at 75°S
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_v_s; j <  j_pol_v_s + 1; j++ ){
            i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
            d_i_max = ( double ) i_max;
            va_Polar_SL_75 = v.x[ 0 ][ j ][ k ];
            if ( is_land ( h, 0, j, k ) ){
                va_Polar_SL_75 = va_Polar_SL_75_land;
                }
            for ( int i = 0; i < i_max; i++ ){
                d_i = ( double ) i;
                v.x[ i ][ j ][ k ] = ( va_Polar_Tropopause_75 - va_Polar_SL_75 ) *
                    d_i / d_i_max + va_Polar_SL_75;
            }
        }
    }
    for ( int j = j_pol_v_s; j < j_pol_v_s + 1; j++ ){
        i_max = im_tropopause[ j ] + GetTropopauseHightAdd ( t_paleo / t_0 );
        d_i_max = ( double ) ( i_max - im_1 );
        if ( d_i_max == 0. ) d_i_max = 1.e-6;
        for ( int k = 0; k < km; k++ ){
            for ( int i = i_max; i < im; i++ ){
                d_i = ( double ) ( i - im_1 );
                v.x[ i ][ j ][ k ] = va_Polar_Tropopause_75 * d_i / d_i_max;
            }
        }
    }

//    goto Printout;

///////////////////////////////////////// end southern polar cell ///////////////////////////////////

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////
/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

////////////// meridional values of w-velocity component from Pol till Ferrel, from Ferrel till Hadley, from Hadley till equator ///////////////

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// w-component formed by the diagonal starting from subtropical jet and North Pole
    d_j_60n = ( double ) j_fer_n;
    d_j_90n = ( double ) j_pol_n;
    d_diff = d_j_60n - d_j_90n;
 
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_n; j < j_fer_n + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_fer_n ][ k ] - u.x[ i ][ j_pol_n ][ k ] ) *
                    ( d_j - d_j_90n ) / d_diff + u.x[ i ][ j_pol_n ][ k ];
                w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_fer_n ][ k ] - w.x[ i ][ j_pol_n ][ k ] ) *
                    ( d_j - d_j_90n ) / d_diff + w.x[ i ][ j_pol_n ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// v-component formed by the diagonal starting from equator till 45°N
    d_j_90n = ( double ) j_pol_n;
    d_j_75n = ( double ) j_pol_v_n;
    d_diff = d_j_75n - d_j_90n;
 
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_n; j < j_pol_v_n + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_pol_v_n ][ k ] - v.x[ i ][ j_pol_n ][ k ] ) *
                    ( d_j - d_j_90n ) / d_diff + v.x[ i ][ j_pol_n ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial polar cell ( from j=0 till j=30 compares to 60° till 90° northern latitude )
// v-component formed by the diagonal starting from 45°N till subtropical jet
    d_j_75n = ( double ) j_pol_v_n;
    d_j_60n = ( double ) j_fer_n;
    d_diff = d_j_60n - d_j_75n;
  
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_v_n; j < j_fer_n + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_fer_n ][ k ] - v.x[ i ][ j_pol_v_n ][ k ] ) *
                    ( d_j - d_j_75n ) / d_diff + v.x[ i ][ j_pol_v_n ][ k ];
            }
        }
    }

//    goto Printout;

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° northern latitude )
// w-component formed by the diagonal starting from tropical jet and subtropical jet
    d_j_60n = ( double ) j_fer_n;
    d_j_30n = ( double ) j_had_n;
    d_diff = d_j_30n - d_j_60n;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_n + 1; j < j_had_n + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_had_n ][ k ] - u.x[ i ][ j_fer_n ][ k ] ) *
                    ( d_j - d_j_60n ) / d_diff + u.x[ i ][ j_fer_n ][ k ];
                w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_had_n ][ k ] - w.x[ i ][ j_fer_n ][ k ] ) *
                    ( d_j - d_j_60n ) / d_diff + w.x[ i ][ j_fer_n ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° northern latitude )
// v-component formed by the diagonal starting from equator till 45°N
    d_j_60n = ( double ) j_fer_n;
    d_j_45n = ( double ) j_fer_v_n;
    d_diff = d_j_45n - d_j_60n;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_n; j < j_fer_v_n + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_fer_v_n ][ k ] - v.x[ i ][ j_fer_n ][ k ] ) *
                    ( d_j - d_j_60n ) / d_diff + v.x[ i ][ j_fer_n ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Ferrel cell ( from j=30 till j=60 compares to 30° till 60° northern latitude )
// v-component formed by the diagonal starting from 45°N till subtropical jet
    d_j_45n = ( double ) j_fer_v_n;
    d_j_30n = ( double ) j_had_n;
    d_diff = d_j_30n - d_j_45n;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_v_n; j < j_had_n + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_had_n ][ k ] - v.x[ i ][ j_fer_v_n ][ k ] ) *
                    ( d_j - d_j_45n ) / d_diff + v.x[ i ][ j_fer_v_n ][ k ];
            }
        }
    }

//    goto Printout;

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// v-component formed by the diagonal starting from equator till 45°S
    d_j_30s = ( double ) j_had_s;
    d_j_45s = ( double ) j_fer_v_s;
    d_diff = d_j_45s - d_j_30s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_v_s; j > j_had_s  - 1; j-- ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_fer_v_s ][ k ] - v.x[ i ][ j_had_s ][ k ] ) *
                    ( d_j - d_j_30s ) / d_diff + v.x[ i ][ j_had_s ][ k ];
            }
        }
    }

//    goto Printout;

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// u- and w-component formed by the diagonal starting from equator and tropical jet
    d_j_5n = ( double ) j_aeq;
    d_j_30n = ( double ) j_had_n;
    d_diff = d_j_5n - d_j_30n;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_n; j < j_aeq + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_aeq ][ k ] - u.x[ i ][ j_had_n ][ k ] ) *
                    ( d_j - d_j_30n ) / d_diff + u.x[ i ][ j_had_n ][ k ];
                w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_aeq ][ k ] - w.x[ i ][ j_had_n ][ k ] ) *
                    ( d_j - d_j_30n ) / d_diff + w.x[ i ][ j_had_n ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// v-component formed by the diagonal starting from equator till 15°N
    d_j_5n = ( double ) j_aeq;
    d_j_15n = ( double ) j_had_v_n;
    d_diff = d_j_5n - d_j_15n;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_v_n; j < j_aeq + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_aeq ][ k ] - v.x[ i ][ j_had_v_n ][ k ] ) *
                    ( d_j - d_j_15n ) / d_diff + v.x[ i ][ j_had_v_n ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// north equatorial Hadley cell ( from j=60 till j=90 compares to 0° till 30° northern latitude )
// v-component formed by the diagonal starting from 15°N till subtropical jet
    d_j_15n = ( double ) j_had_v_n;
    d_j_30n = ( double ) j_had_n;
    d_diff = d_j_15n - d_j_30n;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_n; j < j_had_v_n + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_had_v_n ][ k ] - v.x[ i ][ j_had_n ][ k ] ) *
                    ( d_j - d_j_30n ) / d_diff + v.x[ i ][ j_had_n ][ k ];
            }
        }
    }

/////////////////////////////////////////// change in j-direction in southern hemisphere //////////////////////////////////////////////

//    goto Printout;

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// w-component formed by the diagonal starting from tropical jet and subtropical jet
    d_j_60s = ( double ) j_fer_s;
    d_j_90s = ( double ) j_pol_s;
    d_diff = d_j_60s - d_j_90s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_s; j < j_pol_s + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_fer_s ][ k ] - u.x[ i ][ j_pol_s ][ k ] ) *
                    ( d_j - d_j_90s ) / d_diff + u.x[ i ][ j_pol_s ][ k ];
                w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_fer_s ][ k ] - w.x[ i ][ j_pol_s ][ k ] ) *
                    ( d_j - d_j_90s ) / d_diff + w.x[ i ][ j_pol_s ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// v-component formed by the diagonal starting from 60°S till 75°S
    d_j_75s = ( double ) j_pol_v_s;
    d_j_90s = ( double ) j_pol_s;
    d_diff = d_j_90s - d_j_75s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_s; j > j_pol_v_s - 1; j-- ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_pol_s ][ k ] - v.x[ i ][ j_pol_v_s ][ k ] ) *
                    ( d_j - d_j_75s ) / d_diff + v.x[ i ][ j_pol_v_s ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial polar cell ( from j=150 till j=180 compares to 60° till 90° southern latitude )
// v-component formed by the diagonal starting from 60°S till 75°S
    d_j_60s = ( double ) j_fer_s;
    d_j_75s = ( double ) j_pol_v_s;
    d_diff = d_j_75s - d_j_60s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_pol_v_s; j > j_fer_s - 1; j-- ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_pol_v_s ][ k ] - v.x[ i ][ j_fer_s ][ k ] ) *
                    ( d_j - d_j_60s ) / d_diff + v.x[ i ][ j_fer_s ][ k ];
            }
        }
    }

//    goto Printout;

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// w-component formed by the diagonal starting from polar jet and subtropical jet
    d_j_60s = ( double ) j_fer_s;
    d_j_30s = ( double ) j_had_s;
    d_diff = d_j_30s - d_j_60s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_s; j < j_fer_s + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_had_s ][ k ] - u.x[ i ][ j_fer_s ][ k ] ) *
                    ( d_j - d_j_60s ) / d_diff + u.x[ i ][ j_fer_s ][ k ];
                w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_had_s ][ k ] - w.x[ i ][ j_fer_s ][ k ] ) *
                    ( d_j - d_j_60s ) / d_diff + w.x[ i ][ j_fer_s ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Ferrel cell ( from j=120 till j=150 compares to 30° till 60° southern latitude )
// v-component formed by the diagonal starting from 45°S till subtropical jet
    d_j_45s = ( double ) j_fer_v_s;
    d_j_60s = ( double ) j_fer_s;
    d_diff = d_j_60s - d_j_45s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_fer_s; j > j_fer_v_s - 1; j-- ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_fer_s ][ k ] - v.x[ i ][ j_fer_v_s ][ k ] ) *
                    ( d_j - d_j_45s ) / d_diff + v.x[ i ][ j_fer_v_s ][ k ];
            }
        }
    }

//    goto Printout;

///////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Hadley cell
// u- and w-component formed by the diagonal starting from equatorjet and subtropical jet
    d_j_5s = ( double ) j_aeq;
    d_j_30s = ( double ) j_had_s;
    d_diff = d_j_5s - d_j_30s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_aeq; j < j_had_s + 1; j++ ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                u.x[ i ][ j ][ k ] = ( u.x[ i ][ j_aeq ][ k ] - u.x[ i ][ j_had_s ][ k ] ) *
                    ( d_j - d_j_30s ) / d_diff + u.x[ i ][ j_had_s ][ k ];
                w.x[ i ][ j ][ k ] = ( w.x[ i ][ j_aeq ][ k ] - w.x[ i ][ j_had_s ][ k ] ) *
                    ( d_j - d_j_30s ) / d_diff + w.x[ i ][ j_had_s ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Hadley cell
// v-component formed by the diagonal starting from 15°S till subtropical jet
    d_j_15s = ( double ) j_had_v_s;
    d_j_30s = ( double ) j_had_s;
    d_diff = d_j_30s - d_j_15s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_s; j > j_had_v_s - 1; j-- ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_had_s ][ k ] - v.x[ i ][ j_had_v_s ][ k ] ) *
                    ( d_j - d_j_15s ) / d_diff + v.x[ i ][ j_had_v_s ][ k ];
            }
        }
    }

/////////////////////////////////////////// forming diagonals ///////////////////////////////////////////////////////

// south equatorial Hadley cell
// v-component formed by the diagonal starting from equator till 15°N
    d_j_5s = ( double ) j_aeq;
    d_j_15s = ( double ) j_had_v_s;
    d_diff = d_j_15s - d_j_5s;
    for ( int k = 0; k < km; k++ ){
        for ( int j = j_had_v_s; j > j_aeq - 1; j-- ){
            d_j = ( double ) j;
            for ( int i = 1; i < im; i++ ){
                v.x[ i ][ j ][ k ] = ( v.x[ i ][ j_had_v_s ][ k ] - v.x[ i ][ j_aeq ][ k ] ) *
                    ( d_j - d_j_5s ) / d_diff + v.x[ i ][ j_aeq ][ k ];
            }
        }
    }

//    goto Printout;
/*
///////////////////////////////////////////////// change in sign of v-component /////////////////////////////////////////////////
///////////////////////////////////////////////// values identical with the northern hemisphere //////////////////////////////////////////////////
    for ( int i = 1; i < im; i++ ){
        for ( int j = j_aeq + 1; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                v.x[ i ][ j ][ k ] = - v.x[ i ][ j ][ k ];
            }
        }
    }
*/
/////////////////////////////////////////////// end forming diagonals ///////////////////////////////////////////////////

//    goto Printout;

//    Printout:

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

//    Printout:
                    va_equator_SL = va_equator_SL;
    
///////////////////////////////////////////////// end of smoothing transitions from cell to cell //////////////////////////
///////////////////////////////////////////////// end of smoothing transitions from cell to cell //////////////////////////
}





void BC_Thermo::BC_Surface_v_Velocity ( const string &Name_v_surface_File, Array &v ){
// initial conditions for the Name_v_surface_File at the sea surface
    cout.precision ( 3 );
    cout.setf ( ios::fixed );
    ifstream Name_v_surface_File_Read(Name_v_surface_File);
    if (!Name_v_surface_File_Read.is_open()){
        cerr << "ERROR: could not open Name_v_surface_File file at "
        << Name_v_surface_File << "\n";
        abort();
    }
    k_half = ( km - 1 ) / 2;  // position at 180°E ( Greenwich )
    j = 0;
    k = 0;
    while ( ( k < km ) && !Name_v_surface_File_Read.eof() ){
        while ( j < jm ){
            double lat, lon, v_velocity;
            Name_v_surface_File_Read >> lat;
            Name_v_surface_File_Read >> lon;
            Name_v_surface_File_Read >> v_velocity;
            v.x[ 0 ][ j ][ k ] = v_velocity;
            j++;
        }
    j = 0;
    k++;
    }
}




void BC_Thermo::BC_Surface_w_Velocity ( const string &Name_w_surface_File, Array &w ){
// initial conditions for the Name_w_surface_File at the sea surface
    cout.precision ( 3 );
    cout.setf ( ios::fixed );
    ifstream Name_w_surface_File_Read(Name_w_surface_File);
    if (!Name_w_surface_File_Read.is_open()){
        cerr << "ERROR: could not open Name_w_surface_File file at "
        << Name_w_surface_File << "\n";
        abort();
    }
    k_half = ( km - 1 ) / 2;  // position at 180°E ( Greenwich )
    j = 0;
    k = 0;
    while ( ( k < km ) && !Name_w_surface_File_Read.eof() ){
        while ( j < jm ){
            double lat, lon, w_velocity;
            Name_w_surface_File_Read >> lat;
            Name_w_surface_File_Read >> lon;
            Name_w_surface_File_Read >> w_velocity;
            w.x[ 0 ][ j ][ k ] = w_velocity;
            j++;
        }
    j = 0;
    k++;
    }
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
    k_half = ( km -1 ) / 2;                                             // position at 180°E ( Greenwich )
    j = 0;
    k = 0;
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

    if (!Name_SurfacePrecipitation_File_Read.is_open()){
        cerr << "ERROR: could not open SurfacePrecipitation_File file at "
        << Name_SurfacePrecipitation_File << "\n";
        abort();
    }

    j = 0;
    k = 0;

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




void BC_Thermo::BC_Pressure ( double L_atm, Array_1D &rad, Array &p_stat, Array &p_dyn, Array &t, Array &h ){
    exp_pressure = g / ( 1.e-2 * gam * R_Air );
// boundary condition of surface pressure given by surface temperature through gas equation
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            p_stat.x[ 0 ][ j ][ k ] = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 );  // given in hPa
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 1; i < im; i++ ){
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                p_stat.x[ i ][ j ][ k ] = pow ( ( ( t.x[ 0 ][ j ][ k ] * t_0 - gam * height * 1.e-2 ) /
                    ( t.x[ 0 ][ j ][ k ] * t_0 ) ), exp_pressure ) * p_stat.x[ 0 ][ j ][ k ];
                // international standard athmosphere formula in hPa
                // linear temperature distribution T = T0 - gam * height
            }
        }
    }
}









void BC_Thermo::Latent_Heat ( Array_1D &rad, Array_1D &the, Array_1D &phi, 
                              Array &h, Array &t, Array &tn, Array &u, 
                              Array &v, Array &w, Array &p_dyn, Array &p_stat, 
                              Array &c, Array &ice, Array &Q_Latent, 
                              Array &Q_Sensible, Array &radiation_3D, 
                              Array_2D &Q_radiation, Array_2D &Q_latent, 
                              Array_2D &Q_sensible, Array_2D &Q_bottom ){
    double Q_Latent_Ice = 0.; 
    double step = 0.;
    int i_mount = 0;
    coeff_Lv = lv * r_air * u_0;         // coefficient for Q_latent generated by cloud water
    coeff_Ls = ls * r_air * u_0;         // coefficient for Q_latent generated by cloud ice
    coeff_Q = cp_l * r_air * t_0 * u_0;  // coefficient for Q_Sensible
    double coeff_lat = 10.;  // arbitrary adjustment to fit Q_latent/ Q_sensible == 27%/5% of irradiance
    double coeff_sen = .12;
/*
    step = ( ( exp( zeta * ( rad.z[ 2 ] - 1. ) ) - 1 ) 
                  - ( exp( zeta * ( rad.z[ 0 ] - 1. ) ) - 1 ) ) 
                  * ( L_atm / ( double ) ( im-1 ) ); // local atmospheric shell thickness
    // surface values of latent and sensible heat built from water vapour and temperature gradients
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            i_mount = i_topography[ j ][ k ];
            e = .01 * c.x[ i_mount ][ j ][ k ] * p_stat.x[ i_mount ][ j ][ k ] / ep;  // water vapour pressure in Pa
            a = e / ( R_WaterVapour * t.x[ i_mount ][ j ][ k ] * t_0 );  // absolute humidity in kg/m³
//            Q_Latent.x[ i_mount ][ j ][ k ] = - coeff_Lv * a * u.x[ i_mount ][ j ][ k ] // depends on the definition
            Q_Latent.x[ i_mount ][ j ][ k ] = + coeff_Lv * a  // depends on the definition
                * ( - 3. * c.x[ i_mount ][ j ][ k ] + 4. * c.x[ i_mount+1 ][ j ][ k ] 
                - c.x[ i_mount+2 ][ j ][ k ] ) / step;
//            Q_Latent_Ice = - coeff_Ls * a * u.x[ i_mount ][ j ][ k ] 
            Q_Latent_Ice = + coeff_Ls * a
                * ( - 3. * ice.x[ i_mount ][ j ][ k ] + 4. * ice.x[ i_mount+1 ][ j ][ k ] 
                - ice.x[ i_mount+2 ][ j ][ k ] ) / step;
            Q_Latent.x[ i_mount ][ j ][ k ] = coeff_lat 
                * ( Q_Latent.x[ i_mount ][ j ][ k ] + Q_Latent_Ice );
            Q_Sensible.x[ i_mount ][ j ][ k ] = + coeff_sen * coeff_Q 
//                * u.x[ i_mount ][ j ][ k ] * ( - 3. * t.x[ i_mount ][ j ][ k ] 
                * ( - 3. * t.x[ i_mount ][ j ][ k ] 
                + 4. * t.x[ i_mount+1 ][ j ][ k ] - t.x[ i_mount+2 ][ j ][ k ] ) / step;
        }
    }
*/
    // values of latent and sensible heat above surface built from water vapour and temperature gradients
    for ( int j = 0; j < jm; j++ ){
        for ( int k = 0; k < km; k++ ){
            i_mount = i_topography[ j ][ k ];
//            i_mount = 0.;
            for ( int i = i_mount; i < im-2; i++ ){
                step = ( ( exp( zeta * ( rad.z[ i+1 ] - 1. ) ) - 1 ) 
                         - ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) ) 
                         * ( L_atm / ( double ) ( im-1 ) ); // local atmospheric shell thickness
                e = .01 * c.x[ i ][ j ][ k ] * p_stat.x[ i ][ j ][ k ] / ep;  // water vapour pressure in Pa
                a = e / ( r_water_vapour * t.x[ i ][ j ][ k ] * t_0 );  // absolute humidity in kg/m³
//                Q_Latent.x[ i ][ j ][ k ] = - coeff_Lv * a * u.x[ i ][ j ][ k ] 
                Q_Latent.x[ i ][ j ][ k ] = - coeff_Lv * a 
                    * ( c.x[ i+1 ][ j ][ k ] - c.x[ i ][ j ][ k ] ) / step;
//                Q_Latent_Ice = - coeff_Ls * a * u.x[ i ][ j ][ k ] 
                Q_Latent_Ice = - coeff_Ls * a 
                    * ( ice.x[ i+1 ][ j ][ k ] - ice.x[ i ][ j ][ k ] ) / step;
                Q_Latent.x[ i ][ j ][ k ] = coeff_lat * ( Q_Latent.x[ i ][ j ][ k ] 
                    + Q_Latent_Ice );
                Q_Sensible.x[ i ][ j ][ k ] = - coeff_sen * coeff_Q 
//                    * u.x[ i ][ j ][ k ] * ( t.x[ i+1 ][ j ][ k ] -
                    * ( t.x[ i+1 ][ j ][ k ] -
                    t.x[ i ][ j ][ k ] ) / step;
            }
            for ( int i = 0; i <= i_mount; i++ ){
//                if ( is_land( h, 0, j, k ) ){
// small wiggles in temperature cause strange gradients building Q_Latent and Q_Sensible
                    Q_Latent.x[ i_mount+1 ][ j ][ k ] = Q_Latent.x[ i_mount+4 ][ j ][ k ] 
                        - 3. * Q_Latent.x[ i_mount+3 ][ j ][ k ] 
                        + 3. * Q_Latent.x[ i_mount+2 ][ j ][ k ];
                    Q_Sensible.x[ i_mount+1 ][ j ][ k ] = Q_Sensible.x[ i_mount+4 ][ j ][ k ] 
                        - 3. * Q_Sensible.x[ i_mount+3 ][ j ][ k ] 
                        + 3. * Q_Sensible.x[ i_mount+2 ][ j ][ k ];
                    Q_Latent.x[ i ][ j ][ k ] = Q_Latent.x[ i_mount+4 ][ j ][ k ] 
                        - 3. * Q_Latent.x[ i_mount+3 ][ j ][ k ] 
                        + 3. * Q_Latent.x[ i_mount+2 ][ j ][ k ];
                    Q_Sensible.x[ i ][ j ][ k ] = Q_Sensible.x[ i_mount+4 ][ j ][ k ] 
                        - 3. * Q_Sensible.x[ i_mount+3 ][ j ][ k ] 
                        + 3. * Q_Sensible.x[ i_mount+2 ][ j ][ k ];
//                }
            }
        }
    }
}






void BC_Thermo::Ice_Water_Saturation_Adjustment ( Array_1D &rad, Array &h, 
                Array &c, Array &cn, Array &cloud, Array &cloudn, Array &ice, 
                Array &icen, Array &t, Array &p_stat, Array &S_c_c ){
    if(debug){
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
        assert(!c.has_nan());
        assert(!t.has_nan());
    }
    cout.precision ( 6 );
// Ice_Water_Saturation_Adjustment, distribution of cloud ice and cloud water dependent on water vapour amount and temperature
// constant coefficients for the adjustment of cloud water and cloud ice amount vice versa
    t_00 = 236.15;
    t_Celsius_2 = t_00 - t_0; // in Kelvin = -37 °C

    exp_pressure = g / ( 1.e-2 * gam * R_Air );

    // setting water vapour, cloud water and cloud ice into the proper thermodynamic ratio based on the local temperatures
    // starting from a guessed parabolic temperature and water vapour distribution in north/south direction
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int i = 0; i < im; i++ ){
/** %%%%%%%%%%%%%%%%%%%%%%%%%%     saturation pressure     %%%%%%%%%%%%%%%%%%%%%%%%%%%% **/
                t_u = t.x[ i ][ j ][ k ] * t_0; // in K
                t_Celsius = t_u - t_0; // in C

                p_SL =  .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) 
                         * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                if ( i != 0 )  p_h = pow ( ( ( t.x[ i ][ j ][ k ] * t_0 
                                     - gam * height * 1.e-2 ) /
                                     ( t.x[ i ][ j ][ k ] * t_0 ) ), 
                                     exp_pressure ) * p_SL;
                else  p_h = p_SL;

                r_dry = 100. * p_h / ( R_Air * t_u );
                r_humid = r_dry * ( 1. + c.x[ i ][ j ][ k ] ) 
                          / ( 1. + R_WaterVapour / R_Air * c.x[ i ][ j ][ k ] );
                e_h = .01 * r_humid * R_WaterVapour * t_u;
                a_h = 216.6 * e_h / t_u; // absolute humidity in kg/m3
                E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                E_Ice = hp * exp_func ( t_u, 21.8746, 7.66 ); // saturation water vapour pressure for the ice phase in hPa
                q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg

/** %%%%%%%%%%%%%%%%%%%%%%%%%%%     warm cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% **/

// warm cloud phase in case water vapour is over-saturated
                if ( t_Celsius >= 0. ){
                    q_T = c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ]; // total water content
                    t_u = t.x[ i ][ j ][ k ] * t_0; // in K
                    t_Celsius = t_u - t_0;

                    p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                    height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) 
                             * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                    if ( i != 0 )  p_h = pow ( ( ( t.x[ i ][ j ][ k ] 
                                         * t_0 - gam * height * 1.e-2 ) 
                                         / ( t.x[ i ][ j ][ k ] * t_0 ) ), 
                                         exp_pressure ) * p_SL;
                    else  p_h = p_SL;

                    E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                    q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                    q_Rain_n = q_Rain;

                    if ( q_T <= q_Rain ){ /**     subsaturated     **/
                        c.x[ i ][ j ][ k ] = q_T; // total water amount as water vapour
                        cloud.x[ i ][ j ][ k ] = 0.; // no cloud water available
                        ice.x[ i ][ j ][ k ] = 0.; // no cloud ice available above 0 °C
                        T_it = t_u;

                    }else{ /**     oversaturated     **/
                        iter_prec = 0;
                        while ( iter_prec <= 20 ){ // iter_prec may be varied
                            ++iter_prec;

                            T_it = ( t_u + lv / cp_l * c.x[ i ][ j ][ k ] 
                                   - lv / cp_l * q_Rain );
                            E_Rain = hp * exp_func ( T_it, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                            q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                            q_Rain = .5 * ( q_Rain_n + q_Rain );  // smoothing the iteration process

                            c.x[ i ][ j ][ k ] = q_Rain; // water vapour restricted to saturated water vapour amount
                            cloud.x[ i ][ j ][ k ] = q_T - c.x[ i ][ j ][ k ]; // cloud water amount
                            ice.x[ i ][ j ][ k ] = 0.; // no cloud ice available
                            q_T = c.x[ i ][ j ][ k ] + cloud.x[ i ][ j ][ k ];

                            if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
                            if ( cloud.x[ i ][ j ][ k ] < 0. )  
                                 cloud.x[ i ][ j ][ k ] = 0.;

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
// mixed cloud phase, if 0°C > t > -37°C (t_Celsius_2)
                if ( t_Celsius < 0. ){
                    if ( t_Celsius < t_Celsius_2 )  cloud.x[ i ][ j ][ k ] = 0.;
                    if ( t_Celsius < t_Celsius_2 )  ice.x[ i ][ j ][ k ] = 0.;
                    if ( t_Celsius > 0. )  ice.x[ i ][ j ][ k ] = 0.;
                    q_v_b = c.x[ i ][ j ][ k ];
                    q_c_b = cloud.x[ i ][ j ][ k ];
                    q_i_b = ice.x[ i ][ j ][ k ];
                    q_T = q_v_b + q_c_b + q_i_b; // total water content
                    t_u = t.x[ i ][ j ][ k ] * t_0; // in K
                    T = T_nue = t_u; // in K

                    E_Rain = hp * exp_func ( T, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                    E_Ice = hp * exp_func ( T, 21.8746, 7.66 ); // saturation water vapour pressure for the ice phase in hPa
                    q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                    q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg

                    if ( ( q_c_b > 0. ) && ( q_i_b > 0. ) )
                           q_v_hyp = ( q_c_b * q_Rain + q_i_b * q_Ice ) 
                           / ( q_c_b + q_i_b );
                    if ( ( q_c_b >= 0. ) && ( q_i_b == 0. ) )  q_v_hyp = q_Rain;
                    if ( ( q_c_b == 0. ) && ( q_i_b > 0. ) )  q_v_hyp = q_Ice;

/** §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ **/

                    iter_prec = 0;
                    while ( iter_prec <= 20 ){ // iter_prec may be varied
                        ++iter_prec;

/** condensation == water vapor saturation for cloud water formation, deposition == ice crystal for cloud ice formation **/
                        CND = ( T - t_00 ) / ( t_0 - t_00 );
                        if ( T < t_00 ) CND = 0.;

                        DEP = ( t_0 - T ) / ( t_0 - t_00 );
                        if ( T > t_0 ) DEP = 0.;

                        d_q_v = q_v_hyp - q_v_b;  // changes in water vapour causing cloud water and cloud ice
                        d_q_c = - d_q_v * CND;
                        d_q_i = - d_q_v * DEP;

                        d_t = ( lv * d_q_c + ls * d_q_i ) / cp_l; // in K, temperature changes
                        T = T + d_t; // in K

                        q_v_b = c.x[ i ][ j ][ k ] + d_q_v;  // new values
                        q_c_b = cloud.x[ i ][ j ][ k ] + d_q_c;
                        q_i_b = ice.x[ i ][ j ][ k ] + d_q_i;

                        if ( q_v_b < 0. )  q_v_b = 0.;  // negative values excluded, when iteration starts
                        if ( q_c_b < 0. )  q_c_b = 0.;
                        if ( q_i_b < 0. )  q_i_b = 0.;

                        p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa, needed for saturation water vapour
                        height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) 
                                 * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                        if ( i != 0 ){
                            if( T > gam * height * 1.e-2){
                                p_h = pow ( ( ( T - gam * height * 1.e-2 ) 
                                / ( T ) ), exp_pressure ) * p_SL; // given in hPa
                            }else{
                                logger()<<"WARNING: T is less than gam * height * 1.e-2. "
                                        << __LINE__<<" " << __FILE__<<std::endl; 
                                p_h = p_SL;
                            }
                        }
                        else  p_h = p_SL;

                        E_Rain = hp * exp_func ( T, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        E_Ice = hp * exp_func ( T, 21.8746, 7.66 ); // saturation water vapour pressure for the ice phase in hPa
                        q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                        q_Ice = ep * E_Ice / ( p_h - E_Ice ); // water vapour amount at saturation with ice formation in kg/kg

                        if ( ( q_c_b > 0. ) && ( q_i_b > 0. ) )
                            q_v_hyp = ( q_c_b * q_Rain + q_i_b * q_Ice ) 
                            / ( q_c_b + q_i_b );
                        // average amount of ater vapour based on temperature changes
                        if ( ( q_c_b >= 0. ) && ( q_i_b == 0. ) )  q_v_hyp = q_Rain;
                        if ( ( q_c_b == 0. ) && ( q_i_b > 0. ) )  q_v_hyp = q_Ice;

                        // rate of condensating or evaporating water vapour to form cloud water, 0.5 given by COSMO
//                        S_c_c.x[ i ][ j ][ k ] = .5 * d_q_c / dt_dim;
                        S_c_c.x[ i ][ j ][ k ] = .5 * ( cn.x[ i ][ j ][ k ] 
                                                 - c.x[ i ][ j ][ k ] ) / dt_dim;
                        if ( is_land ( h, i, j, k ) )  S_c_c.x[ i ][ j ][ k ] = 0.;

                        q_T = q_v_b + q_c_b + q_i_b; // total water content, not used except for print out for mass conservation test

                        if( iter_prec >= 3 && (q_v_hyp) > 
                            std::numeric_limits<double>::epsilon() &&
                            fabs ( q_v_b / q_v_hyp - 1. ) <= 1.e-5 )  break;  // make sure q_v_hyp is not 0 divisor

                        q_v_b = .5 * ( q_v_hyp + q_v_b );  // has smoothing effect
                    } // iter_prec end

/** §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     end          iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ **/

                    cn.x[ i ][ j ][ k ] = c.x[ i ][ j ][ k ] = q_v_b;  // new values achieved after converged iterations
                    cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ] = q_c_b;
                    icen.x[ i ][ j ][ k ] = ice.x[ i ][ j ][ k ] = q_i_b;

                    if ( t_Celsius < t_Celsius_2 )  cloudn.x[ i ][ j ][ k ] = 
                                                    cloud.x[ i ][ j ][ k ] = 0.;
                    if ( t_Celsius < t_Celsius_2 )  icen.x[ i ][ j ][ k ] = 
                                                    ice.x[ i ][ j ][ k ] = 0.;
                    if ( t_Celsius > 0. )           icen.x[ i ][ j ][ k ] = 
                                                    ice.x[ i ][ j ][ k ] = 0.;
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
    if(debug){
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
        assert(!c.has_nan());
        assert(!t.has_nan());
    }
}






void BC_Thermo::Two_Category_Ice_Scheme ( Array_1D &rad, Array &h, Array &c, Array &t, Array &p_stat, 
        Array &cloud, Array &ice, Array &P_rain, Array &P_snow, Array &S_v, Array &S_c, Array &S_i, Array &S_r, 
        Array &S_s, Array &S_c_c ){
    if(debug){
        assert(!c.has_nan());
        assert(!t.has_nan());
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
    }
    //  Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, resulting the precipitation distribution formed of rain and snow
    // constant coefficients for the transport of cloud water and cloud ice amount vice versa, rain and snow in the parameterization procedures
    N_i_0 = 1.e2;  // in m-3
    m_i_0 = 1.e-12;  // in kg
    m_i_max = 1.e-9;  // in kg
    m_s_0 = 3.e-9;  // in kg

    c_i_dep = 1.3e-5;  // in m3/(kg*s)
    c_c_au = 4.e-4;  // in 1/s
    c_i_au = 1.e-3;  // in 1/s
    c_ac = .24;  // m2/kg
    c_rim = 18.6;  // m2/kg
    c_agg = 10.3;  // m2/kg
    c_i_cri = .24;  // m2
    c_r_cri = 3.2e-5;  // m2
    a_ev = 1.e-3;  // m2/kg
    b_ev = 5.9;  // m2*s/kg
    c_s_dep = 1.8e-2;  // m2/kg
    b_s_dep = 12.3;  // m2*s/kg
    c_s_melt = 8.43e-5;  // (m2*s)/(K*kg)
    b_s_melt = 12.05;  // m2*s/kg
    a_s_melt = 2.31e3; // K/(kg/kg)
    c_r_frz = 3.75e-2;  // (m2*s)/(K*kg)

    t_nuc = 267.15;  // in K    -6 °C
    t_d = 248.15;  // in K    -25 °C
    t_hn = 236.15;  // in K    -40 °C
    t_r_frz = 271.15;  // in K    -2 °C

    t_1 = 253.15;  // in K    -20 °C
    t_00 = 236.15;  // in K    -40 °C
    t_Celsius_1 = t_1 - t_0;  // -20 °C
    t_Celsius_2 = t_00 - t_0;  // -37 °C

    exp_pressure = g / ( 1.e-2 * gam * R_Air );

    //The m_i is only used in this function and I see no reason it should be a member of this class.
    //So, I move it out of class.
    //I also don't see any reason the above horde of variables should stay as class members. -- mchin
    double m_i = m_i_max; //initialize m_i local variable. If uninitialized, the value in m_i can be anything(bad, very bad). 
    double step = 0.;

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

                p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching

                if ( i != 0 )  p_h = pow ( ( ( t_u - gam * height * 1.e-2 ) / ( t_u ) ), exp_pressure ) * p_SL;
                else  p_h = p_SL;

                r_dry = 100. * p_h / ( R_Air * t_u );  // density of dry air in kg/m³
                r_humid = r_dry * ( 1. + c.x[ i ][ j ][ k ] ) / ( 1. + R_WaterVapour / R_Air * c.x[ i ][ j ][ k ] );
                q_h = c.x[ i ][ j ][ k ];  // threshold value for water vapour at local height h in kg/kg
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

                step = ( ( exp( zeta * ( rad.z[ i + 1 ] - 1. ) ) - 1 ) 
                     - ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) ) 
                     * ( L_atm / ( double ) ( im-1 ) );
                P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ]
                     + r_humid * ( ( S_r.x[ i ][ j ][ k ] - S_r.x[ i + 1 ][ j ][ k ] ) 
                     / 2. * step ) / 2.;  // in kg / ( m2 * s ) == mm/s
                P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ]
                     + r_humid * ( ( S_s.x[ i ][ j ][ k ] - S_s.x[ i + 1 ][ j ][ k ] ) 
                     / 2. * step ) / 2.;

//    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   i = " << i << "   j = " << j << "   k = " << k << "   zeta = " << zeta << "   rad = " << rad.z[ i ] << "   step = " << step << "   height = " << ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) ) << "   P_rain = " << P_rain.x[ i ][ j ][ k ] << "   P_snow = " << P_snow.x[ i ][ j ][ k ] << "   S_r = " << S_r.x[ i ][ j ][ k ] << "   S_s = " << S_s.x[ i ][ j ][ k ] << "   r_humid = " << r_humid << endl;

                if ( P_rain.x[ i ][ j ][ k ] < 0. )  P_rain.x[ i ][ j ][ k ] = 0.;
                if ( P_snow.x[ i ][ j ][ k ] < 0. )  P_snow.x[ i ][ j ][ k ] = 0.;
            }
        }
    }


/******************* main part for rain and snow calculation *********************/

    if ( true ){
        iter_prec = 0;
//        while ( iter_prec <= 5 ){  // iter_prec may be varied, but is sufficient
        while ( iter_prec <= 0 ){  // iter_prec may be varied, but is sufficient
            iter_prec = iter_prec + 1;
            for ( int k = 0; k < km; k++ ){
                for ( int j = 0; j < jm; j++ ){
                    P_rain.x[ im-1 ][ j ][ k ] = 0.;
                    P_snow.x[ im-1 ][ j ][ k ] = 0.;

                    for ( int i = im-2; i >= 0; i-- ){
                        t_u = t.x[ i ][ j ][ k ] * t_0;
                        t_Celsius = t_u - t_0;

                        p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                        height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching

                        if ( i != 0 )  p_h = pow ( ( ( t_u - gam * height * 1.e-2 ) / ( t_u ) ), exp_pressure ) * p_SL;
                        else  p_h = p_SL;

                        r_dry = 100. * p_h / ( R_Air * t_u );  // density of dry air in kg/m³
                        r_humid = r_dry * ( 1. + c.x[ i ][ j ][ k ] ) / ( 1. + R_WaterVapour / R_Air * c.x[ i ][ j ][ k ] );
                        q_h = c.x[ i ][ j ][ k ];  // threshold value for water vapour at local height h in kg/kg
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
                            height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
                            p_t_in = pow ( ( ( t_0 - gam * height * 1.e-2 ) / t_0 ), exp_pressure ) * p_SL;  // given in hPa
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

                        step = ( ( exp( zeta * ( rad.z[ i + 1 ] - 1. ) ) - 1 ) 
                             - ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) ) 
                             * ( L_atm / ( double ) ( im-1 ) );
                        P_rain.x[ i ][ j ][ k ] = P_rain.x[ i + 1 ][ j ][ k ]
                             + r_humid * ( ( S_r.x[ i ][ j ][ k ] - S_r.x[ i + 1 ][ j ][ k ] ) 
                             / 2. * step ) / 2.;  // in kg / ( m2 * s ) == mm/s
                        P_snow.x[ i ][ j ][ k ] = P_snow.x[ i + 1 ][ j ][ k ]
                             + r_humid * ( ( S_s.x[ i ][ j ][ k ] - S_s.x[ i + 1 ][ j ][ k ] ) 
                             / 2. * step ) / 2.;

//    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   i = " << i << "   j = " << j << "   k = " << k << "   zeta = " << zeta << "   rad = " << rad.z[ i ] << "   step = " << step << "   height = " << ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) * ( L_atm / ( double ) ( im-1 ) ) << "   P_rain = " << P_rain.x[ i ][ j ][ k ] << "   P_snow = " << P_snow.x[ i ][ j ][ k ] << "   S_r = " << S_r.x[ i ][ j ][ k ] << "   S_s = " << S_s.x[ i ][ j ][ k ] << "   r_humid = " << r_humid << endl;
/*
                        if ( ( P_rain.x[ 0 ][ j ][ k ] ) >= 7.e-5 )  
                               P_rain.x[ 0 ][ j ][ k ] = 7.e-5;
                        if ( ( P_snow.x[ 0 ][ j ][ k ] ) >= 1.e-5 )  
                               P_snow.x[ 0 ][ j ][ k ] = 1.e-5;
*/
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
    if(debug){
        assert(!c.has_nan());
        assert(!cloud.has_nan());
        assert(!ice.has_nan());
        assert(!P_snow.has_nan());
        assert(!P_rain.has_nan());
    }
}




void BC_Thermo::BC_Evaporation ( Array_1D &rad, Array_2D &vapour_evaporation, Array_2D &Evaporation_Dalton, 
                     Array_2D &Precipitation, Array &h, Array &c, Array &cn ){
// mass flux of water vapour follows the same rules as given for the salinity flux in oceans
// preparations for salinity increase due to evaporation and precipitation differences
// procedure given in Rui Xin Huang, Ocean Circulation, p. 165

    double vapour_surface = 0.;
    double evap_precip = 0.;
    double coeff_vapour = 1.1574e-8;  // 1.1574-8 is the conversion from (Evap-Prec) in mm/d to m/s
//    double coeff_vapour = 1.1574e-8 * 2000.;  
//                                      2000. is fantasy, but it produces a small increase 
//                                      of surface water vapour due to evaporation, TODO
    double rm = rad.z[ 0 ];
    double exp_rm = 1. / exp( zeta * rm );

// additional water vapour as a source term due to evaporation at ocean surface ( i = 0 )
    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            evap_precip = Evaporation_Dalton.y[ j ][ k ] - Precipitation.y[ j ][ k ];

    // this formula contains a 2. order accurate gradient of 1. order, needs 3 points
            vapour_surface = ( - 3. * c.x[ 0 ][ j ][ k ] + 4. * c.x[ 1 ][ j ][ k ] - c.x[ 2 ][ j ][ k ] ) / 
                             ( 2. * dr * exp_rm ) * ( 1. - 2. * c.x[ 0 ][ j ][ k ] ) * evap_precip;     // 2. ord.
    // this formula contains a 1. order accurate gradient of 1. order, needs 2 points
//          vapour_surface = ( c.x[ 0 ][ j ][ k ] - c.x[ 1 ][ j ][ k ] ) / ( dr * exp_rm ) * ( 1. - 2. * c.x[ 0 ][ j ][ k ] ) * evap_precip;

            vapour_evaporation.y[ j ][ k ] = - coeff_vapour * vapour_surface;
            c.x[ 0 ][ j ][ k ] = c.x[ 0 ][ j ][ k ] + vapour_evaporation.y[ j ][ k ];
/*
    cout.precision ( 8 );
    cout.setf ( ios::fixed );
    if ( ( j == 90 ) && ( k == 180 ) ) cout << "   j = " << j << "   k = " << k << "   vapour_evaporation = " << vapour_evaporation.y[ j ][ k ] << "   coeff_vapour = " << coeff_vapour << "   vapour_surface = " << vapour_surface << "   c [g/kg] = " << c.x[ 0 ][ j ][ k ] << "   c_0 = " << c_0 << "   c [kg/kg] = " << c.x[ 0 ][ j ][ k ] * 1000. << "   Evap-Prec = " << evap_precip << "   Evap = " << Evaporation_Dalton.y[ j ][ k ] << "   Prec = " << Precipitation.y[ j ][ k ] << "   c_grad_1 = " << ( c.x[ 0 ][ j ][ k ] - c.x[ 1 ][ j ][ k ] ) / dr << "   c_grad_2 = " << - ( - 3. * c.x[ 0 ][ j ][ k ] + 4. * c.x[ 1 ][ j ][ k ] - c.x[ 2 ][ j ][ k ] ) / ( 2. * dr ) << endl;
*/
        }
     }
 }




void BC_Thermo::Moist_Convection ( Array_1D &rad, Array_1D &the, Array_1D &phi, 
    Array &h, Array &t, Array &u, Array &v, Array &w, Array &p_dyn, Array &p_stat, 
    Array &c, Array &cloud, Array &ice, Array &co2, Array &tn, Array &un, 
    Array &vn, Array &wn, Array &cn, Array &cloudn, Array &icen, Array &co2n, 
    Array &P_rain, Array &P_snow, Array &P_conv, Array &M_u, Array &M_d, 
    Array &MC_s, Array &MC_q, Array &MC_v, Array &MC_w ){
// collection of coefficients for phase transformation
    int i_LFS = 0;
    int i_b = 0;
    double b_u = .3;
    double alf_1 = 5.e-4;
    double alf_2 = .011;
    double p_ps = .05;
    double C_p = p_ps;
    double bet_p = 2.e-3;  // in 1/s
//    double eps_u = 1.e-4;  // in 1/m
    double eps_u = 0.;  // in 1/m
    double gam_d = - 0.3;  // negative from original paper by Tiedtke
    double eps_d = eps_u;
    double del_u = eps_u;
    double del_d = eps_u;
    double delta_i_c = 0.;
    double K_p = 0.;
    double E_u = 0.;
    double E_d = 0.;
    double D_u = 0.;
    double D_d = 0.;
    double humidity_rel = 0.;
    double cloud_buoy = 0.;
    double M_u_denom = 0.;
    double M_d_denom = 0.;
    std::vector<double> s(im, 0);
    std::vector<double> c_u(im, 0);
    std::vector<double> e_d(im, 0);
    std::vector<double> e_l(im, 0);
    std::vector<double> e_p(im, 0);
    std::vector<double> g_p(im, 0);
    std::vector<double> s_u(im, 0);
    std::vector<double> s_d(im, 0);
    std::vector<double> u_u(im, 0);
    std::vector<double> u_d(im, 0);
    std::vector<double> v_u(im, 0);
    std::vector<double> v_d(im, 0);
    std::vector<double> w_u(im, 0);
    std::vector<double> w_d(im, 0);
    std::vector<double> q_v_u(im, 0);
    std::vector<double> q_v_d(im, 0);
    std::vector<double> q_c_u(im, 0);
    std::vector<double> r_humid(im, 0);
    std::vector<double> step(im, 0);

// cumulus cloud model: taken from COSMO based on Tiedtke mass-flux scheme
// parameterisation of moist convection
    for ( int k = 1; k < km-1; k++ ){
        for ( int j = 1; j < jm-1; j++ ){
            i_LFS = 0;
            i_b = 0;
            i_mount = i_topography[ j ][ k ];
            for ( int i = 1; i < im-1; i++ ){
                rm = rad.z[ i ];
                double exp_rm = 1. / exp( zeta * rm );
                sinthe = sin( the.z[ j ] );
                rmsinthe = rm * sinthe;
                q_v_u[ i ] = q_v_d[ i ] = c.x[ i ][ j ][ k ]; 
                q_c_u[ i ] = cloud.x[ i ][ j ][ k ]; 
                t_u = t.x[ i ][ j ][ k ] * t_0; // in K
//                if( i == 1 )  t_u = t.x[ i ][ j ][ k ] * t_0 + 1.; // in K
//                if( i == 9 )  t_u = t.x[ i ][ j ][ k ] * t_0 + 10.; // in K
                t_Celsius = t_u - t_0; // in C
// warm cloud phase in case water vapour is over-saturated
//                if ( t_Celsius >= 0. ){
                q_T = q_v_u[ i ] + q_c_u[ i ]; // total water content
                t_Celsius = t_u - t_0;
                p_SL = .01 * ( r_air * R_Air * t.x[ 0 ][ j ][ k ] * t_0 ); // given in hPa
                height = ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) 
                         * ( L_atm / ( double ) ( im-1 ) );  // coordinate stretching
//                    step[ i ] = ( ( exp( zeta * ( rad.z[ i+1 ] - 1. ) ) - 1 ) 
//                             - ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) ) 
//                             * ( L_atm / ( double ) ( im-1 ) ); // local atmospheric shell thickness
                step[ i ] = ( ( exp( zeta * ( rad.z[ i+1 ] - 1. ) ) - 1 ) 
                         - ( exp( zeta * ( rad.z[ i ] - 1. ) ) - 1 ) );  // local atmospheric shell thickness
                if ( i != 0 )  p_h = pow ( ( ( t.x[ i ][ j ][ k ] 
                                     * t_0 - gam * height * 1.e-2 ) 
                                     / ( t.x[ i ][ j ][ k ] * t_0 ) ), 
                                     exp_pressure ) * p_SL;
                else  p_h = p_SL;
                r_dry = 100. * p_h / ( R_Air * t_u );
                r_humid[ i ] = r_dry * ( 1. + c.x[ i ][ j ][ k ] ) 
                    / ( 1. + R_WaterVapour / R_Air * c.x[ i ][ j ][ k ] );
                E_Rain = hp * exp_func ( t_u, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                q_Rain_n = q_Rain;
                if ( q_T <= q_Rain ){ /**     subsaturated     **/
                    q_v_u[ i ] = q_T; // total water amount as water vapour
                    q_c_u[ i ] = 0.; // no cloud water available
                    T_it = t_u;
                }else{ /**     oversaturated     **/
                    iter_prec = 0;
                    while ( iter_prec <= 20 ){ // iter_prec may be varied
                        ++iter_prec;
                        T_it = ( t_u + lv / cp_l * q_v_u[ i ] 
                               - lv / cp_l * q_Rain );
                        e_h = q_v_u[ i ] * p_stat.x[ i ][ j ][ k ] / ep;  // water vapour pressure in hPa
                        a_h = 216.6 * e_h / t_u; // absolute humidity in kg/m3
                        E_Rain = hp * exp_func ( T_it, 17.2694, 35.86 ); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                        q_Rain = .5 * ( q_Rain_n + q_Rain );  // smoothing the iteration process
                        q_v_u[ i ] = q_Rain; // water vapour restricted to saturated water vapour amount
                        q_c_u[ i ] = q_T - q_v_u[ i ]; // cloud water amount
                        q_T = q_v_u[ i ] + q_c_u[ i ];
                        if ( q_v_u[ i ] < 0. )  q_v_u[ i ] = 0.;
                        if ( q_c_u[ i ] < 0. )  q_c_u[ i ] = 0.;
                        if( ( q_Rain_n ) > std::numeric_limits<double>::epsilon() &&
                            fabs ( q_Rain / q_Rain_n - 1. ) <= 1.e-5 )    break;  // make sure q_Rain_n is not 0 divisor
                        q_Rain_n = q_Rain;

//    cout << "    in warm cloud   iter_prec = " << iter_prec << "   i = " << i << "   q_v_u = " << q_v_u[ i ] << "   q_c_u = " << q_c_u[ i ] << "   q_T = " << q_T << "   q_Rain = " << q_Rain << "   q_Rain_n = " << q_Rain_n << endl; 

                    }  // end iter_prec
                }  // end oversaturated

//                cn.x[ i ][ j ][ k ] = c.x[ i ][ j ][ k ];
//                cloudn.x[ i ][ j ][ k ] = cloud.x[ i ][ j ][ k ];
                t.x[ i ][ j ][ k ] = T_it / t_0;

                if ( q_v_u[ i ] < 0. )  q_v_u[ i ] = 0.;
                if ( q_c_u[ i ] < 0. )  q_c_u[ i ] = 0.;
                humidity_rel = e_h / E_Rain * 100.;
/** %%%%%%%%%%%%%%%%%%%%%%%%%%%     end          warm cloud phase     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% **/

                double dcdr = ( c.x[ i+1 ][ j ][ k ] - c.x[ i-1 ][ j ][ k ] ) 
                              / ( 2. * dr ) * exp_rm;  // Tiedtke included this term
                double dcdthe = ( c.x[ i ][ j+1 ][ k ] - c.x[ i ][ j-1 ][ k ] ) 
                                / ( 2. * dthe ) / rm;
                double dcdphi = ( c.x[ i ][ j ][ k+1 ] - c.x[ i ][ j ][ k-1 ] ) 
                                / ( 2. * dphi ) / rmsinthe;
                dcdr = 0.;
//                    dcdthe = 0.;
//                    dcdphi = 0.;
// dry static energy
                s[ i ] = cp_l * t.x[ i ][ j ][ k ] * t_0 + g * height;
// level of the cloud base ( level of free convection )
                if( ( humidity_rel >= 90. ) 
                    && ( u.x[ i ][ j ][ k ] >= 0. ) 
                    && ( i_b == 0 ) )  i_b = i;
//   if ( ( j == 90 ) && ( k == 180 ) )  cout << "   LCB    i = " << i << "   i_mount = " << i_mount << "   i_b = " << i_b << "   i_LFS = " << i_LFS << endl;
// specification of entrainment in the updraft
                if( ( i_b != 0 ) && ( i >= i_b )
                    && ( u_u[ i ] >= 0. ) ){
                    M_u.x[ i_b ][ j ][ k ] = r_humid[ i ] * u.x[ i_b ][ j ][ k ];  // updraft at cloud base
                    E_u = - r_humid[ i ] / c.x[ i ][ j ][ k ] 
                        * ( u.x[ i ][ j ][ k ] * dcdr + v.x[ i ][ j ][ k ] * dcdthe 
                        + w.x[ i ][ j ][ k ] * dcdphi )
                        + eps_u * M_u.x[ i ][ j ][ k ];
                    u_u[ i ] = u.x[ i ][ j ][ k ] + M_u.x[ i ][ j ][ k ] 
                               / ( dthe / rm * dphi / rmsinthe );
                }else  E_u = 0.;
// level of free sinking
                cloud_buoy = ( q_v_u[ i ] + q_Rain );
                if ( ( cloud_buoy < c.x[ i ][ j ][ k ] ) 
                    && ( i_LFS == 0 ) )  i_LFS = i;
                // updraft at cloud base == downdraft at level of free sinking (LFS)
                M_d.x[ i_LFS ][ j ][ k ] = gam_d * M_u.x[ i_b ][ j ][ k ];  // downdraft at cloud top
                u_d[ i ] = u.x[ i ][ j ][ k ] + M_d.x[ i ][ j ][ k ] 
                               / ( dthe / rm * dphi / rmsinthe );

//   if ( ( j == 90 ) && ( k == 180 ) )  cout << "   LFS    i = " << i << "   i = " << i-1 << "   i_mount = " << i_mount << "   i_b = " << i_b << "   i_LFS = " << i_LFS << "   M_u = " << M_u.x[ i ][ j ][ k ] << "   E_u = " << E_u << "   M_d = " << M_d.x[ i ][ j ][ k ] << endl;

// specification of detrainment in the updraft
                if ( i == i_LFS )  D_u = ( 1. - b_u ) * M_u.x[ i ][ j ][ k ] / step[ i ]
                                        + del_u * M_u.x[ i ][ j ][ k ];
                else               D_u = 0.;
                if ( i == i_LFS + 1 )  D_u = b_u * M_u.x[ i ][ j ][ k ] / step[ i ]
                                           + del_u * M_u.x[ i ][ j ][ k ];
                if ( i > i_LFS + 1 )  D_u = 0.; 
// specification of entrainment and detrainment in the udowndraft
                eps_d = eps_u;
                E_d = eps_d * fabs( M_d.x[ i ][ j ][ k ] );
                D_d = del_d * fabs( M_d.x[ i ][ j ][ k ] );
// evaporation of cloud water in the environment
                e_l[ i ] = D_u / r_humid[ i ] * q_c_u[ i ];
// evaporation of precipitation below cloud base
                if ( i <= i_b )
                      e_p[ i ] = C_p * alf_1 * ( q_Rain - q_v_u[ i ] ) 
                      * sqrt ( p_ps / alf_2 * P_conv.x[ i ][ j ][ k ] / C_p );
                else  e_p[ i ] = 0.;
// evaporation of precipitation within the downdraft
                e_d[ i ] = q_v_d[ i ] - q_Rain;
                if ( e_d[ i ] > 0. )  q_v_d[ i ] = q_Rain;
                if ( e_d[ i ] < 0. )  e_d[ i ] = 0.;
// condensation/deposition within the updraft
                c_u[ i ] = q_c_u[ i ];
                if ( c_u[ i ] >= q_Rain )  q_v_u[ i ] = q_Rain;
                if ( c_u[ i ] < 0. )  c_u[ i ] = 0.;
// formation of precipitation within the updraft
                if ( is_land( h, 0, j, k ) )  delta_i_c = 3000.;
                else                          delta_i_c = 1500.;  // in m
                double height_it = ( exp( zeta * ( rad.z[ i_LFS ] - 1. ) ) - 1 ) 
                    * ( L_atm / ( double ) ( im-1 ) );  // in m
                if ( height > height_it + delta_i_c )  K_p = 0.;
                else  K_p = bet_p;
                g_p[ i ] = K_p * q_c_u[ i ];
// cloud model
// convective updraft
                if ( ( i > i_b-1 ) && ( i_b != 0 ) && ( i_LFS == 0 )
                    && ( u.x[ i ][ j ][ k ] >= 0. ) ){
                    M_u.x[ i+1 ][ j ][ k ] = M_u.x[ i ][ j ][ k ] 
                                             + step[ i ] * ( E_u - D_u );  // integration from cloud base upwards
                    M_u_denom = M_u.x[ i ][ j ][ k ] + step[ i ] * D_u;
                    if ( M_u_denom == 0. )  M_u_denom = 1.e+6;
                    s_u[ i ] = step[ i ] * ( E_u * s[ i ] + lv * r_humid[ i ] * c_u[ i ] ) / M_u_denom;
                    q_v_u[ i ] = step[ i ] * ( E_u * c.x[ i ][ j ][ k ] - r_humid[ i ] * c_u[ i ] ) 
                        / M_u_denom;
                    q_c_u[ i ] = step[ i ] * r_humid[ i ] *  ( c_u[ i ] - g_p[ i ] ) / M_u_denom;
                    v_u[ i ] = step[ i ] * ( E_u * v.x[ i ][ j ][ k ] ) / M_u_denom;
                    w_u[ i ] = step[ i ] * ( E_u * w.x[ i ][ j ][ k ] ) / M_u_denom;
                    if ( M_u.x[ i ][ j ][ k ] < 0. )  M_u.x[ i ][ j ][ k ] = 0.;
                    if ( M_u.x[ i+1 ][ j ][ k ] < 0. )  M_u.x[ i+1 ][ j ][ k ] = 0.;

    cout.precision ( 8 );
    if ( ( j == 90 ) && ( k == 180 ) )  cout << "   updraft   i = " << i << "   i_mount = " << i_mount << "   i_b = " << i_b << "   i_LFS = " << i_LFS << "  s = " << s[ i ] << "  s_u = " << s_u[ i ] << "  s_d = " << s_d[ i ] << "  E_u = " << E_u << "  D_u = " << D_u << "  E_d = " << E_d << "  D_d = " << D_d << "  e_l = " << e_l[ i ] << "  e_p = " << e_p[ i ] << "  e_d = " << e_d[ i ] << "  c_u = " << c_u[ i ] << "  g_p = " << g_p[ i ] << "  M_u = " << M_u.x[ i ][ j ][ k ] << "  M_u_1 = " << M_u.x[ i+1 ][ j ][ k ] << "  M_d = " << M_d.x[ i ][ j ][ k ] << "  MC_s = " << MC_s.x[ i ][ j ][ k ] << "  MC_q = " << MC_q.x[ i ][ j ][ k ] << "  MC_v = " << MC_v.x[ i ][ j ][ k ] << "  MC_w = " << MC_w.x[ i ][ j ][ k ] << "  c = " << c.x[ i ][ j ][ k ] << "  q_Rain = " << q_Rain << "  cloud_buoy = " << cloud_buoy - c.x[ i ][ j ][ k ] << "  q_v_u = " << q_v_u[ i ] << "  q_v_d = " << q_v_d[ i ] << "  cloud = " << cloud.x[ i ][ j ][ k ] << "  q_c_u = " << q_c_u[ i ] << "  u = " << u.x[ i ][ j ][ k ] << "  u_u = " << u_u[ i ] << "  u_d = " << u_d[ i ] << "  v = " << v.x[ i ][ j ][ k ] << "  v_u = " << v_u[ i ] << "  v_d = " << v_d[ i ] << "  w = " << w.x[ i ][ j ][ k ] << "  w_u = " << w_u[ i ] << "  w_d = " << w_d[ i ] << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  T_it = " << T_it - t_0 << "  P_conv = " << P_conv.x[ i ][ j ][ k ] << "  e_h = " << e_h << "  E_Rain = " << E_Rain << "  humidity_rel = " << humidity_rel << "  r_dry = " << r_dry << "  r_humid = " << r_humid[ i ] << "  a_h = " << a_h << "  step = " << step[ i ] << "  height = " << height << "  dcdthe = " << dcdthe << "  dcdphi = " << dcdphi << endl;

                    if( ( t.x[ i ][ j ][ k ] * t_0 - t_0 ) <= -37. )  M_u.x[ i ][ j ][ k ] = 0.;
                }  // convective updraft
            }  // end i-loop

    if ( ( j == 90 ) && ( k == 180 ) )  cout << endl << endl;

// convective downdraft
            for ( int i = im-1; i >= 1; i-- ){
                if( i <= i_LFS ){
                   if( ( t.x[ i ][ j ][ k ] * t_0 - t_0 ) <= -37. )  M_d.x[ i ][ j ][ k ] = 0.;
                   M_d.x[ i-1 ][ j ][ k ] = M_d.x[ i ][ j ][ k ] + step[ i ] * ( E_d - D_d );  // integration from LFS downwards
                   M_d_denom = M_d.x[ i ][ j ][ k ] + step[ i ] * D_d;
                   if ( M_d_denom == 0. )  M_d_denom = 1.e+6;
                   s_d[ i ] = step[ i ] * ( E_d * s[ i ] - lv * r_humid[ i ] * e_d[ i ] ) / M_d_denom;
                   q_v_d[ i ] = step[ i ] * ( E_d * c.x[ i ][ j ][ k ] 
                       + r_humid[ i ] * e_d[ i ] ) / M_d_denom;
                   v_d[ i ] = step[ i ] * ( E_d * v.x[ i ][ j ][ k ] ) / M_d_denom;
                   w_d[ i ] = step[ i ] * ( E_d * w.x[ i ][ j ][ k ] ) / M_d_denom;

    cout.precision ( 8 );
    if ( ( j == 90 ) && ( k == 180 ) )  cout << "   downdraft   i = " << i << "   i_mount = " << i_mount << "   i_b = " << i_b << "   i_LFS = " << i_LFS << "  s = " << s[ i ] << "  s_u = " << s_u[ i ] << "  s_d = " << s_d[ i ] << "  E_u = " << E_u << "  D_u = " << D_u << "  E_d = " << E_d << "  D_d = " << D_d << "  e_l = " << e_l[ i ] << "  e_p = " << e_p[ i ] << "  e_d = " << e_d[ i ] << "  c_u = " << c_u[ i ] << "  g_p = " << g_p[ i ] << "  M_u = " << M_u.x[ i ][ j ][ k ] << "  M_u_1 = " << M_u.x[ i+1 ][ j ][ k ] << "  M_d = " << M_d.x[ i ][ j ][ k ] << "  M_d_1 = " << M_d.x[ i-1 ][ j ][ k ] << "  MC_s = " << MC_s.x[ i ][ j ][ k ] << "  MC_q = " << MC_q.x[ i ][ j ][ k ] << "  MC_v = " << MC_v.x[ i ][ j ][ k ] << "  MC_w = " << MC_w.x[ i ][ j ][ k ] << "  c = " << c.x[ i ][ j ][ k ] << "  q_Rain = " << q_Rain << "  cloud_buoy = " << cloud_buoy - c.x[ i ][ j ][ k ] << "  q_v_u = " << q_v_u[ i ] << "  q_v_d = " << q_v_d[ i ] << "  cloud = " << cloud.x[ i ][ j ][ k ] << "  q_c_u = " << q_c_u[ i ] << "  u = " << u.x[ i ][ j ][ k ] << "  v = " << v.x[ i ][ j ][ k ] << "  v_u = " << v_u[ i ] << "  v_d = " << v_d[ i ] << "  w = " << w.x[ i ][ j ][ k ] << "  w_u = " << w_u[ i ] << "  w_d = " << w_d[ i ] << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  T_it = " << T_it - t_0 << "  P_conv = " << P_conv.x[ i ][ j ][ k ] << "  e_h = " << e_h << "  E_Rain = " << E_Rain << "  humidity_rel = " << humidity_rel << "  r_dry = " << r_dry << "  r_humid = " << r_humid[ i ] << "  a_h = " << a_h << "  step = " << step[ i ] << "  height = " << height << endl;

// rain water formed by cloud convection
                   P_conv.x[ i-1 ][ j ][ k ] = P_conv.x[ i ][ j ][ k ] 
                       - step[ i ] * r_humid[ i ] * ( g_p[ i ] - e_d[ i ] - e_p[ i ] );
                   if ( P_conv.x[ i-1 ][ j ][ k ] < 0. )
                        P_conv.x[ i-1 ][ j ][ k ] = 0.;
                   if ( is_land( h, i, j, k ) )  P_conv.x[ i-1 ][ j ][ k ] = 0.;
               }  // convective downdraft
           }  // end i-loop



// RHS for thermodynamic forcing due to moist convection
            for ( int i = 1; i < im-1; i++ ){
//                if( ( T_it - t_0 ) <= -37. )  M_u.x[ i ][ j ][ k ] = 
//                    M_d.x[ i ][ j ][ k ] = M_u.x[ i+1 ][ j ][ k ] = 
//                    M_d.x[ i+1 ][ j ][ k ] = 0.;
                MC_s.x[ i ][ j ][ k ] = - ( ( M_u.x[ i + 1 ][ j ][ k ] 
                    * ( s_u[ i ] - s[ i ] ) + M_d.x[ i + 1 ][ j ][ k ] * ( s_d[ i ] - s[ i ] ) ) 
                    - ( M_u.x[ i ][ j ][ k ] * ( s_u[ i ] - s[ i ] ) 
                    + M_d.x[ i ][ j ][ k ] * ( s_d[ i ] - s[ i ] ) ) ) 
                    / ( cp_l * r_humid[ i ] * step[ i ] ) 
                    + lv / cp_l * ( c_u[ i ] - e_d[ i ] - e_l[ i ] - e_p[ i ] );
                MC_q.x[ i ][ j ][ k ] = - ( ( M_u.x[ i + 1 ][ j ][ k ] 
                    * ( q_v_u[ i ] - c.x[ i + 1 ][ j ][ k ] ) 
                    + M_d.x[ i + 1 ][ j ][ k ] * ( q_v_d[ i ] - c.x[ i + 1 ][ j ][ k ] ) ) 
                    - ( M_u.x[ i ][ j ][ k ] * ( q_v_u[ i ] - c.x[ i ][ j ][ k ] ) 
                    + M_d.x[ i ][ j ][ k ] * ( q_v_d[ i ] - c.x[ i ][ j ][ k ] ) ) ) 
                    / ( r_humid[ i ] * step[ i ] ) - ( c_u[ i ] - e_d[ i ] - e_l[ i ] - e_p[ i ] );
                MC_v.x[ i ][ j ][ k ] = - ( ( M_u.x[ i + 1 ][ j ][ k ] 
                    * ( v_u[ i ] - v.x[ i + 1 ][ j ][ k ] ) 
                    + M_d.x[ i + 1 ][ j ][ k ] * ( v_d[ i ] - v.x[ i + 1 ][ j ][ k ] ) ) 
                    - ( M_u.x[ i ][ j ][ k ] * ( v_u[ i ] - v.x[ i ][ j ][ k ] ) 
                    + M_d.x[ i ][ j ][ k ] * ( v_d[ i ] - v.x[ i ][ j ][ k ] ) ) ) 
                    / ( r_humid[ i ] * step[ i ] );
                MC_w.x[ i ][ j ][ k ] = - ( ( M_u.x[ i + 1 ][ j ][ k ] 
                    * ( w_u[ i ] - w.x[ i + 1 ][ j ][ k ] ) 
                    + M_d.x[ i + 1 ][ j ][ k ] * ( w_d[ i ] - w.x[ i + 1 ][ j ][ k ] ) ) 
                    - ( M_u.x[ i ][ j ][ k ] * ( w_u[ i ] - w.x[ i ][ j ][ k ] ) 
                    + M_d.x[ i ][ j ][ k ] * ( w_d[ i ] - w.x[ i ][ j ][ k ] ) ) ) 
                    / ( r_humid[ i ] * step[ i ] );

//    cout.precision ( 8 );
//    if ( ( j == 90 ) && ( k == 180 ) )  cout << "   i = " << i << "   i_mount = " << i_mount << "   i_b = " << i_b << "   i_LFS = " << i_LFS << "  s = " << s[ i ] << "  s_u = " << s_u[ i ] << "  s_d = " << s_d[ i ] << "  E_u = " << E_u << "  D_u = " << D_u << "  E_d = " << E_d << "  D_d = " << D_d << "  e_l = " << e_l << "  e_p = " << e_p << "  e_d = " << e_d << "  c_u = " << c_u << "  g_p = " << g_p << "  M_u = " << M_u.x[ i ][ j ][ k ] << "  M_u_add = " << M_u_add << "  M_d = " << M_d.x[ ii ][ j ][ k ] << "  MC_s = " << MC_s.x[ i ][ j ][ k ] << "  MC_q = " << MC_q.x[ i ][ j ][ k ] << "  MC_v = " << MC_v.x[ i ][ j ][ k ] << "  MC_w = " << MC_w.x[ i ][ j ][ k ] << "  c = " << c.x[ i ][ j ][ k ] << "  q_Rain = " << q_Rain << "  cloud_buoy = " << cloud_buoy - c.x[ i ][ j ][ k ] << "  q_v_u = " << q_v_u[ i ] << "  q_v_d = " << q_v_d[ i ] << "  cloud = " << cloud.x[ i ][ j ][ k ] << "  q_c_u = " << q_c_u[ i ] << "  u = " << u.x[ i ][ j ][ k ] << "  v = " << v.x[ i ][ j ][ k ] << "  v_u = " << v_u[ i ] << "  v_d = " << v_d[ i ] << "  w = " << w.x[ i ][ j ][ k ] << "  w_u = " << w_u[ i ] << "  w_d = " << w_d[ i ] << "  t = " << t.x[ i ][ j ][ k ] * t_0 - t_0 << "  T_it = " << T_it - t_0 << "  P_conv = " << P_conv.x[ i ][ j ][ k ] << "  e_h = " << e_h << "  E_Rain = " << E_Rain << "  humidity_rel = " << humidity_rel << "  r_dry = " << r_dry << "  r_humid = " << r_humid[ i ] << "  a_h = " << a_h << "  step = " << step[ i ] << "  height = " << height << "  dcdthe = " << dcdthe << "  dcdphi = " << dcdphi << endl;

            }  // end i-loop
        }
    }


// boundaries of various variables
    for ( int k = 1; k < km-1; k++ ){
        for ( int j = 1; j < jm-1; j++ ){
            M_u.x[ 0 ][ j ][ k ] = c43 * M_u.x[ 1 ][ j ][ k ] -
                c13 * M_u.x[ 2 ][ j ][ k ];
            M_u.x[ im-1 ][ j ][ k ] = c43 * M_u.x[ im-2 ][ j ][ k ] -
                c13 * M_u.x[ im-3 ][ j ][ k ];
            M_d.x[ 0 ][ j ][ k ] = c43 * M_d.x[ 1 ][ j ][ k ] -
                c13 * M_d.x[ 2 ][ j ][ k ];
            M_d.x[ im-1 ][ j ][ k ] = c43 * M_d.x[ im-2 ][ j ][ k ] -
                c13 * M_d.x[ im-3 ][ j ][ k ];
        }
    }

    for ( int k = 0; k < km; k++ ){
        for ( int i = 0; i < im; i++ ){
            M_u.x[ i ][ 0 ][ k ] = c43 * M_u.x[ i ][ 1 ][ k ] -
                 c13 * M_u.x[ i ][ 2 ][ k ];
            M_u.x[ i ][ jm-1 ][ k ] = c43 * M_u.x[ i ][ jm-2 ][ k ] -
                 c13 * M_u.x[ i ][ jm-3 ][ k ];
            M_d.x[ i ][ 0 ][ k ] = c43 * M_d.x[ i ][ 1 ][ k ] -
                 c13 * M_d.x[ i ][ 2 ][ k ];
            M_d.x[ i ][ jm-1 ][ k ] = c43 * M_d.x[ i ][ jm-2 ][ k ] -
                 c13 * M_d.x[ i ][ jm-3 ][ k ];
        }
    }

    for ( int i = 0; i < im; i++ ){
        for ( int j = 1; j < jm-1; j++ ){
            M_u.x[ i ][ j ][ 0 ] = c43 * M_u.x[ i ][ j ][ 1 ] -
                c13 * M_u.x[ i ][ j ][ 2 ];
            M_u.x[ i ][ j ][ km-1 ] = c43 * M_u.x[ i ][ j ][ km-2 ] -
                c13 * M_u.x[ i ][ j ][ km-3 ];
            M_u.x[ i ][ j ][ 0 ] = M_u.x[ i ][ j ][ km-1 ] =
                ( M_u.x[ i ][ j ][ 0 ] + M_u.x[ i ][ j ][ km-1 ] ) / 2.;
            M_d.x[ i ][ j ][ 0 ] = c43 * M_d.x[ i ][ j ][ 1 ] -
                c13 * M_d.x[ i ][ j ][ 2 ];
            M_d.x[ i ][ j ][ km-1 ] = c43 * M_d.x[ i ][ j ][ km-2 ] -
                c13 * M_d.x[ i ][ j ][ km-3 ];
            M_d.x[ i ][ j ][ 0 ] = M_d.x[ i ][ j ][ km-1 ] =
                ( M_d.x[ i ][ j ][ 0 ] + M_d.x[ i ][ j ][ km-1 ] ) / 2.;
        }
    }

/*
// formation of water vapour gradients along mountain sides
            // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in radial direction
                if ( i < im - 2 ){
                    if ( ( is_land ( h, i, j, k ) ) && ( ( is_air ( h, i+1, j, k ) ) && ( is_air ( h, i+2, j, k ) ) ) ){
                        dcdr = ( - 3. * u.x[ i ][ j ][ k ] + 4. * u.x[ i + 1 ][ j ][ k ] - u.x[ i + 2 ][ j ][ k ] ) / ( 2. * dr ) * exp_rm;
                    }
                }

            // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral southward direction
                if ( ( j >= 2 ) && ( j <= jm - 3 ) ){
                    if ( ( is_land ( h, i, j, k ) ) && ( ( is_air ( h, i, j+1, k ) ) && ( is_air ( h, i, j+2, k ) ) ) ){
                        dcdthe = ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i ][ j + 1 ][ k ] - c.x[ i ][ j + 2 ][ k ] ) / ( 2. * dthe );
                    }

            // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral northward direction
                    if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j-1, k ) ) && ( is_air ( h, i, j-2, k ) ) ){
                        dcdthe = ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i ][ j - 1 ][ k ] - c.x[ i ][ j - 2 ][ k ] ) / ( 2. * dthe );
                    }

            // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral southward direction
                    if ( ( ( is_land ( h, i, j, k ) ) 
                        && ( ( is_air ( h, i, j+1, k ) ) && ( is_land ( h, i, j+2, k ) ) ) ) 
                        || ( ( j == jm - 2 )
                        && ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j+1, k ) ) ) ) ){
                        dcdthe = ( c.x[ i ][ j + 1 ][ k ] - c.x[ i ][ j ][ k ] ) / dthe;
                    }

            // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral northward direction
                    if ( ( ( is_land ( h, i, j, k ) ) 
                        && ( ( is_air ( h, i, j-1, k ) ) && ( is_land ( h, i, j-2, k ) ) ) ) 
                        || ( ( j == 1 )
                        && ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j-1, k ) ) ) ) ){
                        dcdthe = ( c.x[ i ][ j - 1 ][ k ] - c.x[ i ][ j ][ k ] ) / dthe;
                    }
                }

            // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in longitudial eastward direction
                if ( ( k >= 2 ) && ( k <= km - 3 ) ){
                    if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k+1 ) ) && ( is_air ( h, i, j, k+2 ) ) ){
                        dcdphi = ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i ][ j ][ k + 1 ] - c.x[ i ][ j ][ k + 2 ] ) / ( 2. * dphi );
                    }

            // 2. order accurate finite differences 1. and 2. order starting from solid surfaces in longitudial westward direction
                    if ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k-1 ) ) && ( is_air ( h, i, j, k-2 ) ) ){
                        dcdphi = ( - 3. * c.x[ i ][ j ][ k ] + 4. * c.x[ i ][ j ][ k - 1 ] - c.x[ i ][ j ][ k - 2 ] ) / ( 2. * dphi );
                    }

            // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral eastward direction
                    if ( ( ( is_land ( h, i, j, k ) ) 
                        && ( ( is_air ( h, i, j, k+1 ) ) && ( is_land ( h, i, j, k+2 ) ) ) ) 
                        || ( ( k == km - 2 )
                        && ( ( is_air ( h, i, j, k ) ) && ( is_land ( h, i, j, k+1 ) ) ) ) ){
                        dcdphi = ( c.x[ i ][ j ][ k + 1 ] - c.x[ i ][ j ][ k ] ) / dphi;
                    }

            // 1. order accurate finite differences 1. and 2. order starting from solid surfaces in lateral westward direction
                    if ( ( ( is_land ( h, i, j, k ) ) 
                        && ( ( is_air ( h, i, j, k-1 ) ) && ( is_land ( h, i, j, k-2 ) ) ) ) 
                        || ( ( k == 1 )
                        && ( ( is_land ( h, i, j, k ) ) && ( is_air ( h, i, j, k-1 ) ) ) ) ){
                        dcdphi = ( c.x[ i ][ j ][ k - 1 ] - c.x[ i ][ j ][ k ] ) / dphi;
                    }
                }
*/
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
                if ( c.x[ i ][ j ][ k ] >= .02 )  c.x[ i ][ j ][ k ] = .02;
//                if ( c.x[ i ][ j ][ k ] >= .03 )  c.x[ i ][ j ][ k ] = .03;
                if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
//                if ( cloud.x[ i ][ j ][ k ] >= .008 )  cloud.x[ i ][ j ][ k ] = .008;
                if ( cloud.x[ i ][ j ][ k ] >= .01 )  cloud.x[ i ][ j ][ k ] = .01;
                if ( cloud.x[ i ][ j ][ k ] < 0. )  cloud.x[ i ][ j ][ k ] = 0.;
//                if ( ice.x[ i ][ j ][ k ] >= .0025 )  ice.x[ i ][ j ][ k ] = .0025;
                if ( ice.x[ i ][ j ][ k ] >= .0005 )  ice.x[ i ][ j ][ k ] = .0005;
                if ( ice.x[ i ][ j ][ k ] < 0. )  ice.x[ i ][ j ][ k ] = 0.;
                if ( co2.x[ i ][ j ][ k ] >= 8.92 )  co2.x[ i ][ j ][ k ] = 8.92; // == 2500 ppm
                if ( co2.x[ i ][ j ][ k ] <= 1. )  co2.x[ i ][ j ][ k ] = 1.; // == 280 ppm

                if ( is_land ( h, i, j, k ) ){
                    u.x[ i ][ j ][ k ] = 0.;
                    v.x[ i ][ j ][ k ] = 0.;
                    w.x[ i ][ j ][ k ] = 0.;
//                    t.x[ i ][ j ][ k ] = 1.;  // = 273.15 K
//                    c.x[ i ][ j ][ k ] = 0.;
/*
                    if ( t.x[ i ][ j ][ k ] >= 1.165 )  t.x[ i ][ j ][ k ] = 1.165;  // == 45 °C
                    if ( t.x[ i ][ j ][ k ] <= - .78 )  t.x[ i ][ j ][ k ] = - .78;  // == 59.82 °C
                    if ( c.x[ i ][ j ][ k ] >= .03 )  c.x[ i ][ j ][ k ] = .03;
                    if ( c.x[ i ][ j ][ k ] < 0. )  c.x[ i ][ j ][ k ] = 0.;
*/
                    cloud.x[ i ][ j ][ k ] = 0.;
                    ice.x[ i ][ j ][ k ] = 0.;
//                    co2.x[ i ][ j ][ k ] = 1.;  // = 280 ppm
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



int BC_Thermo::GetTropopauseHightAdd(double t_pal){
    double d_i_h_round = round((t_pal * t_0) / 2.6);
        // adiabatic slope of radial temperature 0.65/100m, stepsize 400m => 2.6/400m
    return ( int ) d_i_h_round;
}


double BC_Thermo::GetPoleTemperature(int Ma, int Ma_1, int Ma_2, double t_1, double t_2){
    return (t_2 - t_1) / (double) (Ma_2 - Ma_1) * (double) (Ma - Ma_1) + t_1;
}
/*
double BC_Thermo::GetPoleTemperature(int Ma, const std::map<int, double> &pole_temp_map){
    assert(pole_temp_map.size()>1);
    std::pair<int, double> up = *pole_temp_map.begin(), bottom = *++pole_temp_map.begin();
    if(Ma <= pole_temp_map.begin()->first || Ma > (--pole_temp_map.end())->first){
        return t_pole; // when Ma out of boundary
    }

    for( const auto& pair : pole_temp_map ){
        if(pair.first>=Ma){
            bottom = pair;
            break;
        }else{
            up = pair;
        }
    }
    return GetPoleTemperature(Ma, up.first, bottom.first, up.second, bottom.second);
}
*/

double BC_Thermo::GetPoleTemperature(int Ma, const std::map<int, double> &pole_temp_map){
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
    return GetPoleTemperature(Ma, up.first, bottom.first, up.second, bottom.second);
}

