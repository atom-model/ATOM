#include "cHydrosphereModel.h"
#include "MinMax_Hyd.h"

void cHydrosphereModel::print_welcome_msg(){
   if (verbose) {
        cout << endl << endl << endl;
        cout << "***** Hydrosphere General Circulation Model ( OGCM ) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 1 additional transport equations to describe the salinity" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** salinity is part of the Boussinesq approximation" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )" << endl << endl;

        //cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
        has_printed_welcome_msg = true;    
    }
}

void cHydrosphereModel::print_min_max()
{
    // 3D_fields
    // searching of maximum and minimum values of temperature
    string str_max_temperature = " max temperature ", str_min_temperature = " min temperature ", str_unit_temperature = "C";
    MinMax_Hyd      minmaxTemperature ( im, jm, km, u_0, c_0, L_hyd );
    minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

    //  searching of maximum and minimum values of u-component
    string str_max_u = " max 3D u-component ", str_min_u = " min 3D u-component ", str_unit_u = "mm/s";
    MinMax_Hyd      minmax_u ( im, jm, km, u_0, c_0, L_hyd );
    minmax_u.searchMinMax_3D ( str_max_u, str_min_u, str_unit_u, u, h );

    //  searching of maximum and minimum values of v-component
    string str_max_v = " max 3D v-component ", str_min_v = " min 3D v-component ", str_unit_v = "m/s";
    MinMax_Hyd      minmax_v ( im, jm, km, u_0, c_0, L_hyd );
    minmax_v.searchMinMax_3D ( str_max_v, str_min_v, str_unit_v, v, h );

    //  searching of maximum and minimum values of w-component
    string str_max_w = " max 3D w-component ", str_min_w = " min 3D w-component ", str_unit_w = "m/s";
    MinMax_Hyd      minmax_w ( im, jm, km, u_0, c_0, L_hyd );
    minmax_w.searchMinMax_3D ( str_max_w, str_min_w, str_unit_w, w, h );

    //      searching of maximum and minimum values of pressure
    string str_max_pressure = " max pressure dynamic ", str_min_pressure = " min pressure dynamic ", str_unit_pressure = "hPa";
    MinMax_Hyd      minmaxPressure ( im, jm, km, u_0, c_0, L_hyd );
    minmaxPressure.searchMinMax_3D ( str_max_pressure, str_min_pressure, str_unit_pressure, p_dyn, h );

    //      searching of maximum and minimum values of static pressure
    string str_max_pressure_stat = " max pressure static ", str_min_pressure_stat = " min pressure static ", str_unit_pressure_stat = "bar";
    MinMax_Hyd      minmaxPressure_stat ( im, jm, km, u_0, c_0, L_hyd );
    minmaxPressure_stat.searchMinMax_3D ( str_max_pressure_stat, str_min_pressure_stat, str_unit_pressure_stat, p_stat, h );

    cout << endl << " salinity based results in the three dimensional space: " << endl << endl;
    //  searching of maximum and minimum values of salt concentration
    string str_max_salt_concentration = " max salt concentration ", str_min_salt_concentration = " min salt concentration ", str_unit_salt_concentration = "psu";
    MinMax_Hyd      minmaxSalt ( im, jm, km, u_0, c_0, L_hyd );
    minmaxSalt.searchMinMax_3D ( str_max_salt_concentration, str_min_salt_concentration, str_unit_salt_concentration, c, h );

    //  searching of maximum and minimum values of salt balance
    string str_max_salt_balance = " max salt balance ", str_min_salt_balance = " min salt balance ", str_unit_salt_balance = "psu";
    MinMax_Hyd      minmaxSaltBalance ( im, jm, km, u_0, c_0, L_hyd );
    minmaxSaltBalance.searchMinMax_3D ( str_max_salt_balance, str_min_salt_balance, str_unit_salt_balance, Salt_Balance, h );

    //  searching of maximum and minimum values of salt finger
    string str_max_salt_finger = " max salt finger ", str_min_salt_finger = " min salt finger ", str_unit_salt_finger = "psu";
    MinMax_Hyd      minmaxSaltFinger ( im, jm, km, u_0, c_0, L_hyd );
    minmaxSaltFinger.searchMinMax_3D ( str_max_salt_finger, str_min_salt_finger, str_unit_salt_finger, Salt_Finger, h );

    //  searching of maximum and minimum values of salt diffusion
    string str_max_salt_diffusion = " max salt diffusion ", str_min_salt_diffusion = " min salt diffusion ", str_unit_salt_diffusion = "psu";
    MinMax_Hyd      minmaxSaltDiffusion ( im, jm, km, u_0, c_0, L_hyd );
    minmaxSaltDiffusion.searchMinMax_3D ( str_max_salt_diffusion, str_min_salt_diffusion, str_unit_salt_diffusion, Salt_Diffusion, h );

    //  searching of maximum and minimum values of buoyancy force
    string str_max_BuoyancyForce_3D = " max buoyancy force ", str_min_BuoyancyForce_3D = " min buoyancy force ", str_unit_BuoyancyForce_3D = "kN/m2";
    MinMax_Hyd      minmaxBuoyancyForce_3D ( im, jm, km, u_0, c_0, L_hyd );
    minmaxBuoyancyForce_3D.searchMinMax_3D ( str_max_BuoyancyForce_3D, str_min_BuoyancyForce_3D, str_unit_BuoyancyForce_3D, BuoyancyForce_3D, h );
    // 2D_fields

    //  searching of maximum and minimum values of total salt volume in a column
    string str_max_salt_total = " max salt total ", str_min_salt_total = " min salt total ", str_unit_salt_total = "psu";
    MinMax_Hyd      minmaxSalt_total ( jm, km, c_0 );
    minmaxSalt_total.searchMinMax_2D ( str_max_salt_total, str_min_salt_total, str_unit_salt_total, Salt_total, h );

    //  searching of maximum and minimum values of salt finger volume in a column
    string str_max_Salt_Finger = " max Salt_Finger ", str_min_Salt_Finger = " min Salt_Finger ", str_unit_Salt_Finger = "psu";
    MinMax_Hyd      minmaxSalt_finger ( jm, km, c_0 );
    minmaxSalt_finger.searchMinMax_2D ( str_max_Salt_Finger, str_min_Salt_Finger, str_unit_Salt_Finger, SaltFinger, h );

    //  searching of maximum and minimum values of salt diffusion volume in a column
    string str_max_Salt_Diffusion = " max Salt_Diffusion ", str_min_Salt_Diffusion = " min Salt_Diffusion ", str_unit_Salt_Diffusion = "psu";
    MinMax_Hyd      minmaxSalt_diffusion ( jm, km, c_0 );
    minmaxSalt_diffusion.searchMinMax_2D ( str_max_Salt_Diffusion, str_min_Salt_Diffusion, str_unit_Salt_Diffusion, SaltDiffusion, h );

    //  searching of maximum and minimum values of salt diffusion volume in a column
    string str_max_BuoyancyForce_2D = " max BuoyancyForce_2D ", str_min_BuoyancyForce_2D = " min BuoyancyForce_2D ", str_unit_BuoyancyForce_2D = "N";
    MinMax_Hyd      minmaxBuoyancyForce_2D ( jm, km, c_0 );
    minmaxBuoyancyForce_2D.searchMinMax_2D ( str_max_BuoyancyForce_2D, str_min_BuoyancyForce_2D, str_unit_BuoyancyForce_2D, BuoyancyForce_2D, h );

    cout << endl << " deep currents averaged for a two dimensional plane: " << endl << endl;

    //  searching of maximum and minimum values of upwelling volume in a column
    string str_max_upwelling = " max upwelling ", str_min_upwelling = " min upwelling ", str_unit_upwelling = "m/s";
    MinMax_Hyd      minmaxUpwelling ( jm, km, c_0 );
    minmaxUpwelling.searchMinMax_2D ( str_max_upwelling, str_min_upwelling, str_unit_upwelling, Upwelling, h );

    //  searching of maximum and minimum values of downwelling volume in a column
    string str_max_downwelling = " max downwelling ", str_min_downwelling = " min downwelling ", str_unit_downwelling = "m/s";
    MinMax_Hyd      minmaxDownwelling ( jm, km, c_0 );
    minmaxDownwelling.searchMinMax_2D ( str_max_downwelling, str_min_downwelling, str_unit_downwelling, Downwelling, h );

    //  searching of maximum and minimum values of bottom water volume in a column
    string str_max_bottom_water = " max bottom water ", str_min_bottom_water = " min bottom water ", str_unit_bottom_water = "m/s";
    MinMax_Hyd      minmaxBottom_water ( jm, km, c_0 );
    minmaxBottom_water.searchMinMax_2D ( str_max_bottom_water, str_min_bottom_water, str_unit_bottom_water, EkmanPumping, h );

    //  searching of maximum and minimum values of the bathymetry
    string str_max_bathymetry = " max bathymetry ", str_min_bathymetry = " min bathymetry ", str_unit_bathymetry = "m";
    MinMax_Hyd      minmaxBathymetry ( jm, km, c_0 );
    minmaxBathymetry.searchMinMax_2D ( str_max_bathymetry, str_min_bathymetry, str_unit_bathymetry, Bathymetry, h );
}
