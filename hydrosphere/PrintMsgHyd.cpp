#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "cHydrosphereModel.h"

void cHydrosphereModel::print_welcome_msg(){
   if(verbose){
        cout << endl << endl << endl;
        cout << "***** Hydrosphere General Circulation Model(OGCM) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 1 additional transport equations to describe the salinity" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** salinity is part of the Boussinesq approximation" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)" << endl << endl;

        //cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
        has_printed_welcome_msg = true;    
    }
}
/*
 * 
*/
void cHydrosphereModel::print_min_max_hyd(){
    cout << endl << " flow properties: " << endl << endl;
    searchMinMax_3D(" max temperature ", " min temperature ", 
        " deg", t, 273.15, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max u-component ", " min u-component ", "m/s", u, u_0);
    searchMinMax_3D(" max v-component ", " min v-component ", "m/s", v, u_0);
    searchMinMax_3D(" max w-component ", " min w-component ", "m/s", w, u_0);
    searchMinMax_3D(" max pressure dynamic ", " min pressure dynamic ", "hPa", p_dyn, r_0_water * u_0 * u_0 * 1e-2);
    searchMinMax_3D(" max pressure static ", " min pressure static ", "bar", p_stat, 1.);
    searchMinMax_3D(" max water density ", " min water density ", "kg/m3", r_water, 1.);
    searchMinMax_3D(" max salt water density ", " min salt water density ", "kg/m3", r_salt_water, 1.);
    cout << endl << " salinity based results in the three dimensional space: " << endl << endl;
    searchMinMax_3D(" max salt concentration ", " min salt concentration ", "psu", c, c_0);
    searchMinMax_3D(" max salt balance ", " min salt balance ", "kg/m3", Salt_Balance, 1.);
    searchMinMax_3D(" max salt finger ", " min salt finger ", "kg/m3", Salt_Finger, 1.);
    searchMinMax_3D(" max salt diffusion ", " min salt diffusion ", "kg/m3", Salt_Diffusion, 1.);
    cout << endl << " forces per unit volume: " << endl << endl;
    searchMinMax_3D(" max pressure force ", " min pressure force ", "N/m3", PressureGradientForce, 1.);
    searchMinMax_3D(" max buoyancy force ", " min buoyancy force ", "N/m3", BuoyancyForce, 1.);
    searchMinMax_3D(" max Coriolis force ", " min Coriolis force ", "N/m3", CoriolisForce, 1.);
    cout << endl << " salt concentration averaged for the two dimensional surface plane: " << endl << endl;
    searchMinMax_2D(" max salt total ", " min salt total ", "psu", Salt_total, c_0);
    searchMinMax_2D(" max Salt_Finger ", " min Salt_Finger ", "kg/m3", SaltFinger, 1.);
    searchMinMax_2D(" max Salt_Diffusion ", " min Salt_Diffusion ", "kg/m3", SaltDiffusion, 1.);
    searchMinMax_2D(" max BuoyancyForce_2D ", " min BuoyancyForce_2D ", "N", BuoyancyForce_2D, 1.);
    cout << endl << " deep currents averaged for a two dimensional plane: " << endl << endl;
    searchMinMax_2D(" max EkmanPumping ", " min EkmanPumping ", "m/a", EkmanPumping, 1.);
    searchMinMax_2D(" max upwelling ", " min upwelling ", "m/a", Upwelling, 1.);
    searchMinMax_2D(" max downwelling ", " min downwelling ", "m/a", Downwelling, 1.);
    searchMinMax_2D(" max bathymetry ", " min bathymetry ", "m", Bathymetry, 1.);
}
