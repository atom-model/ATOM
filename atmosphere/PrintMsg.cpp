#include "cAtmosphereModel.h"
#include "MinMax_Atm.h"

void cAtmosphereModel::print_welcome_msg()
{
    if (verbose) {
        cout << endl << endl << endl;
        cout << "***** Atmosphere General Circulation Model ( AGCM ) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 4 additional transport equations to describe the water vapour, cloud water, cloud ice and co2 concentration" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** water vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code ( German Weather Forecast )" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** co2 concentration appears in the absorptivity of the radiation models" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )" << endl << endl;

        //cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
        has_welcome_msg_printed = true;
    }
}

void cAtmosphereModel::print_final_remarks()
{
    //  final remarks
    cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
    cout << "\n\n\n\n";
}

void cAtmosphereModel::print_min_max_values()
{
    MinMax_Atm min_max_3d( im, jm, km );

    //  searching of maximum and minimum values of temperature
    min_max_3d.searchMinMax_3D( " max 3D temperature ", " min 3D temperature ", "°C", t, h, 273.15,
                                [](double i)->double{return i - 273.15;},
                                true );

    //  searching of maximum and minimum values of u-component
    min_max_3d.searchMinMax_3D ( " max 3D u-component ", " min 3D u-component ", "m/s", u, h, u_0);

    //  searching of maximum and minimum values of v-component
    min_max_3d.searchMinMax_3D ( " max 3D v-component ", " min 3D v-component ", "m/s", v, h, u_0 );

    //  searching of maximum and minimum values of w-component
    min_max_3d.searchMinMax_3D ( " max 3D w-component ", " min 3D w-component ", "m/s", w, h, u_0 );

    //  searching of maximum and minimum values of dynamic pressure
    min_max_3d.searchMinMax_3D ( " max 3D pressure dynamic ", " min 3D pressure dynamic ", "hPa", p_dyn, h, 0.768 ); // 0.768 = 0.01 * r_air *u_0*u_0 in hPa

    //  searching of maximum and minimum values of static pressure
    min_max_3d.searchMinMax_3D ( " max 3D pressure static ", " min 3D pressure static ", "hPa", p_stat, h );

    cout << endl << " energies in the three dimensional space: " << endl << endl;

    //  searching of maximum and minimum values of radiation_3D
    min_max_3d.searchMinMax_3D ( " max 3D radiation ",  " min 3D radiation ",  "W/m2", radiation_3D, h );

//  searching of maximum and minimum values of sensible heat
    min_max_3d.searchMinMax_3D ( " max 3D sensible heat ", " min 3D sensible heat ", "W/m2", Q_Sensible, h );

    //  searching of maximum and minimum values of latency
    min_max_3d.searchMinMax_3D ( " max 3D latent heat ", " min 3D latent heat ", "W/m2", Q_Latent, h );

    cout << endl << " greenhouse gases: " << endl << endl;

    //  searching of maximum and minimum values of water vapour
    min_max_3d.searchMinMax_3D ( " max 3D water vapour ",  " min 3D water vapour ", "g/kg", c, h, 1000. );

    //  searching of maximum and minimum values of cloud water
    min_max_3d.searchMinMax_3D ( " max 3D cloud water ", " min 3D cloud water ", "g/kg", cloud, h, 1000. );

    //  searching of maximum and minimum values of cloud ice
    min_max_3d.searchMinMax_3D ( " max 3D cloud ice ", " min 3D cloud ice ", "g/kg", ice, h, 1000. );

    //  searching of maximum and minimum values of rain precipitation
    min_max_3d.searchMinMax_3D ( " max 3D rain ", " min 3D rain ", "mm/d", P_rain, h, 8.46e4 );

    //  searching of maximum and minimum values of snow precipitation
    min_max_3d.searchMinMax_3D (  " max 3D snow ", " min 3D snow ", "mm/d", P_snow, h, 8.46e4 );

    //  searching of maximum and minimum values of co2
    min_max_3d.searchMinMax_3D ( " max 3D co2 ", " min 3D co2 ", "ppm", co2, h, 280. );

    //  searching of maximum and minimum values of epsilon
    min_max_3d.searchMinMax_3D ( " max 3D epsilon ",  " min 3D epsilon ", "%", epsilon_3D, h );

    //  searching of maximum and minimum values of buoyancy force
    min_max_3d.searchMinMax_3D (  " max 3D buoyancy force ", " min 3D buoyancy force ", "kN/m2", BuoyancyForce, h );

    // 2D-fields

    cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" << endl <<
        " results based on two dimensional considerations of the problem" << endl;

    cout << endl << " co2 distribution row-wise: " << endl << endl;

    MinMax_Atm  min_max_2d ( jm, km );

    //  searching of maximum and minimum values of co2 total
    min_max_2d.searchMinMax_2D ( " max co2_total ", " min co2_total ", " ppm ", co2_total, h, 280. );

    cout << endl << " precipitation: " << endl << endl;

    //  searching of maximum and minimum values of precipitation
    min_max_2d.searchMinMax_2D (  " max precipitation ", " min precipitation ", "mm/d", Precipitation, h, 1. );
    max_Precipitation = min_max_2d.out_maxValue (  );

    //  searching of maximum and minimum values of precipitable water
    min_max_2d.searchMinMax_2D ( " max precipitable water ", " min precipitable water ", "mm",
                                 precipitable_water, h, 1. );

    cout << endl << " energies at see level without convection influence: " << endl << endl;

    //  searching of maximum and minimum values of radiation
    min_max_2d.searchMinMax_2D ( " max 2D Q radiation ", " min 2D Q radiation ",  "W/m2", Q_radiation, h );

    //  searching of maximum and minimum values of latent energy
    min_max_2d.searchMinMax_2D ( " max 2D Q latent ", " min 2D Q latent ", "W/m2", Q_latent, h );

    //  searching of maximum and minimum values of sensible energy
    min_max_2d.searchMinMax_2D ( " max 2D Q sensible ", " min 2D Q sensible ", "W/m2", Q_sensible, h );
    //  searching of maximum and minimum values of bottom heat
    min_max_2d.searchMinMax_2D ( " max 2D Q bottom ", " min 2D Q bottom heat ", "W/m2", Q_bottom, h );

    cout << endl << " secondary data: " << endl << endl;

    //  searching of maximum and minimum values of Evaporation
    min_max_2d.searchMinMax_2D ( " max heat Evaporation ", " min heat Evaporation ", " kJ/kg", Q_Evaporation, h );

    //  searching of maximum and minimum values of Evaporation by Dalton
    min_max_2d.searchMinMax_2D ( " max Evaporation Dalton ", " min Evaporation Dalton ", "mm/d", Evaporation_Dalton, h );

    //  searching of maximum and minimum values of Evaporation by Penman
    min_max_2d.searchMinMax_2D ( " max Evaporation Penman ", " min Evaporation Penman ", "mm/d", Evaporation_Penman, h );

    cout << endl << " properties of the atmosphere at the surface: " << endl << endl;

    //  searching of maximum and minimum values of albedo
    min_max_2d.searchMinMax_2D (  " max 2D albedo ", " min 2D albedo ", "%", albedo, h );

    //  searching of maximum and minimum values of epsilon
    min_max_2d.searchMinMax_2D ( " max 2D epsilon ", " min 2D epsilon ", "%", epsilon, h );

    //  searching of maximum and minimum values of topography
    min_max_2d.searchMinMax_2D ( " max 2D topography ", " min 2D topography ", "m", Topography, h );
}

