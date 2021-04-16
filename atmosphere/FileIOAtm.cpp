/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include <iostream>
#include <fstream>

#include "cAtmosphereModel.h"

using namespace std;

void cAtmosphereModel::Atmosphere_v_w_Transfer(const string &Name_Bathymetry_File){
    string Name_v_w_Transfer_File = output_path + "/[" + Name_Bathymetry_File + "]_Transfer_Atm.vw";
    ofstream v_w_Transfer_File;
    v_w_Transfer_File.precision(4);
    v_w_Transfer_File.setf(ios::fixed);
    v_w_Transfer_File.open(Name_v_w_Transfer_File);
    if(!v_w_Transfer_File.is_open()){
        cout << "ERROR: transfer file name in atmosphere: " << Name_v_w_Transfer_File << "\n";
        cerr << "ERROR: could not open transfer file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
// velocity components in m/s, p_dyn in hPa, Precipitation and Evaporation_Dalton in mm/d
//            v_w_Transfer_File << v.x[0][j][k] * u_0 << " " << w.x[0][j][k] * u_0 << " " << 
//                t.x[0][j][k] << " " << p_dyn.x[0][j][k] << " " << Evaporation_Dalton.y[j][k] << 
//                " " << Precipitation.y[j][k] << endl;
            v_w_Transfer_File << v.x[0][j][k] * u_0 << " " << w.x[0][j][k] * u_0 << " " << 
                t.x[0][j][k] << " " << p_dyn.x[0][j][k] << " " << Evaporation_Dalton.y[j][k] << 
                " " << Precipitation.y[j][k] << endl;
        }
    }
    v_w_Transfer_File.close();
}


void cAtmosphereModel::Atmosphere_PlotData (string &Name_Bathymetry_File, int iter_cnt){
    string Name_PlotData_File = output_path + "/[" + Name_Bathymetry_File + "]_PlotData_Atm"+
        (iter_cnt > 0 ? "_"+to_string(iter_cnt) : "") + ".xyz";
    ofstream PlotData_File(Name_PlotData_File);
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);
    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    PlotData_File << "lons(deg)" << ", " << "lats(deg)" << ", " << "topography(m)" << ", " << "v-velocity(m/s)" << ", " 
        << "w-velocity(m/s)" << ", " << "velocity-mag(m/s)" << ", " << "temperature(Celsius)" << ", " 
        << "water_vapour(g/kg)" << ", " << "precipitation(mm/d)" << ", " <<  "precipitable water(mm)" << ", " 
        << "Evaporation_Dalton (mm/d) " <<endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double vel_mag = sqrt(pow(v.x[0][j][k] * u_0, 2) + pow(w.x[0][j][k] * u_0, 2));
            PlotData_File << k << " " << 90-j << " " << h.x[0][j][k] << " " << v.x[0][j][k] * u_0 << " " 
                << w.x[0][j][k] * u_0 << " " << vel_mag << " " << t.x[0][j][k] * t_0 - t_0 << " " 
                << c.x[0][j][k] * 1000. << " " << Precipitation.y[j][k] << " " << precipitable_water.y[j][k] 
                << " " <<  Evaporation_Dalton.y[j][k] << " " << endl;
        }
    }
}

