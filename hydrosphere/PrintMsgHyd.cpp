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
void cHydrosphereModel::print_final_remarks(){
    cout << endl 
    << "***** end of the Hydrosphere General Circulation Modell (OGCM) *****" 
        << endl << endl;
    cout << "\n\n";
    return;
}
/*
*
*/
void cHydrosphereModel::print_loop_3D_headings(){
    cout << endl << endl;
    cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    cout << " 3D OGCM iterational process" << endl;
    cout << " max total iteration number nm = " << nm << endl << endl;
    cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum \
        and current number of velocity iterations, maximum and current number of pressure iterations " << endl
        << endl << " Ma = " << (int)*get_current_time() 
        << "     n = " << iter_cnt << "    velocity_iter_max = " 
        << velocity_iter_max 
        << "     velocity_iter = " << velocity_iter 
        << "    pressure_iter_max = " << pressure_iter_max 
        << "    pressure_iter = " << pressure_iter << endl;
    return;
}
/*
*
*/
/*
*
*/
