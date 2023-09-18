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
        has_printed_welcome_msg = true;
    }
    return;
}
/*
*
*/
void cAtmosphereModel::print_final_remarks(){
    cout << endl 
    << "***** end of the Atmosphere General Circulation Modell (AGCM) *****" 
        << endl << endl;
    cout << "\n\n";
    return;
}
/*
*
*/
void cAtmosphereModel::print_loop_3D_headings(){
    cout.precision(7);
    cout << endl << endl;
    cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    cout << " 3D Euler flow for a AGCM iterational process solved by 4. order Runge-Kutta scheme" << endl;
    if(use_stretched_coordinate_system)
        cout << endl << "      coordinate stretching is in use" << endl << endl;
    cout << " max total iteration number nm = " << nm << endl << endl;
    cout << " present state of the computation: " << endl 
        << endl << " current time slice" 
        << endl << "      Ma = " << (int)*get_current_time() << endl;

    cout << endl << " number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations, control when to write output files (every how many pressure iterations) " << endl        << "      n = " << iter_cnt << "    velocity_iter_max = " 
        << velocity_iter_max 
        << "     velocity_iter = " << velocity_iter 
        << "    pressure_iter_max = " << pressure_iter_max 
        << "    pressure_iter = " << pressure_iter
        << "    checkpoint = " << checkpoint << endl;
    return;
}
/*
*
*/

