#include <iostream>
#include <fstream>
#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "cHydrosphereModel.h"

using namespace std;


void cHydrosphereModel::Atmosphere_v_w_Transfer(const string &Name_Bathymetry_File){
    ifstream v_w_Transfer_File;
    string Name_v_w_Transfer_File = output_path + "/[" + Name_Bathymetry_File + "]_Transfer_Atm.vw";
    v_w_Transfer_File.precision(4);
    v_w_Transfer_File.setf(ios::fixed);
    v_w_Transfer_File.open(Name_v_w_Transfer_File);
    if(!v_w_Transfer_File.is_open()){
        cout << "ERROR: transfer file name in hydrosphere: " << Name_v_w_Transfer_File << "\n";
        cerr << "ERROR: could not open transfer file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            v_w_Transfer_File >> v.x[im-1][j][k];  //dimensional in m/s
            v_w_Transfer_File >> w.x[im-1][j][k];  //dimensional in m/s
            v_w_Transfer_File >> t.x[im-1][j][k];  //dimensional in m/s
            v_w_Transfer_File >> p_dyn.x[im-1][j][k];  //dimensional in hPa
            v_w_Transfer_File >> Evaporation_Dalton.y[j][ k];  //dimensional in mm/d
            v_w_Transfer_File >> Precipitation.y[j][k];  //dimensional in mm/d
        }
    }
    v_w_Transfer_File.close();
}


void cHydrosphereModel::Hydrosphere_PlotData(const string &Name_Bathymetry_File, int iter_cnt){
    string path = output_path + "/[" + Name_Bathymetry_File + "]_PlotData_Hyd"+
        (iter_cnt > 0 ? "_" + (iter_cnt) : "") + ".xyz";
    ofstream PlotData_File(path);
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);
    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    PlotData_File << "lons(deg)" << ", " << "lats(deg)" << ", " << "topography" << ", " << "v-velocity(m/s)" << ", " 
        << "w-velocity(m/s)" << ", " << "velocity-mag(m/s)" << ", " << "temperature(Celsius)" << ", " << "salinity(psu)" 
        << ", " << "Ekman_pumping(m/s)" << ", " <<  "upwelling(m/s)" << ", " <<  "downwelling(m/s)" << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double vel_mag = sqrt(pow(v.x[im-1][j][k] * u_0 , 2) + pow(w.x[im-1][j][k] * u_0, 2));
            PlotData_File << k << " " << 90-j << " " << h.x[im-1][j][k] << " " << v.x[im-1][j][k] * u_0 
            << " " << w.x[im-1][j][k] * u_0 << " " << vel_mag << " " << t.x[im-1][j][k] * t_0 - t_0 
            << " " << c.x[im-1][j][k] << " " << EkmanPumping.y[j][k] << " " << Upwelling.y[j][k] 
            << "   " << -Downwelling.y[j][k] << " " <<  endl;
        }
    }
}

