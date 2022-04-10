#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cstdlib>

#include "cHydrosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

void cHydrosphereModel::read_Hydrosphere_Surface_Data(int Ma){
cout << endl << "      OGCM: read_Hydrosphere_Surface_Data ......................." << endl;
    if(!is_global_temperature_curve_loaded()) 
        load_global_temperature_curve();
    if(!is_equat_temperature_curve_loaded()) 
        load_equat_temperature_curve();
    if(!is_pole_temperature_curve_loaded()) 
        load_pole_temperature_curve();
    //Prepare the temperature, precipitation and salinity, Name_Transfer_File
//    cout.precision(6);
//    cout.setf(ios::fixed);
    string Name_Transfer_File;
    stringstream ssName_v_w_Transfer_File;
    string Name_SurfaceTemperature_File = temperature_file;
    string Name_SurfaceNASATemperature_File  = temperature_file;
    string Name_SurfacePrecipitation_File = precipitation_file;
    string Name_SurfaceSalinity_File = salinity_file;
    if(Ma != 0 && use_earthbyte_reconstruction){
        Name_SurfaceTemperature_File = output_path + "/" 
            + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";
//        Name_SurfacePrecipitation_File = output_path + "/" 
//            + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";    
        Name_SurfaceSalinity_File = output_path + "/" + std::to_string(Ma) 
            + "Ma_Reconstructed_Salinity.xyz";
        struct stat info;
        if(stat(Name_SurfaceSalinity_File.c_str(), &info) != 0){
//           stat(Name_SurfacePrecipitation_File.c_str(), &info) != 0 ||
            std::string cmd_str = "python " + reconstruction_script_path 
                + " " + std::to_string(Ma - time_step) + " " + std::to_string(Ma) 
                + " " + output_path + " " + BathymetrySuffix +" hyd";
            int ret = system(cmd_str.c_str());
            std::cout << " reconstruction script returned: " << ret << std::endl;
        }
    }
    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    if(!has_printed_welcome_msg)  print_welcome_msg();
    HydrosphereDataTransfer(bathymetry_name);
    cout << endl << "      bathymetry given by the x-y-z data set:    " 
        << bathymetry_name.c_str();
    init_bathymetry(bathymetry_path + "/" + bathymetry_name);
    if(use_NASA_velocity){
        read_IC(velocity_v_file, v.x[im-1], jm, km);
        read_IC(velocity_w_file, w.x[im-1], jm, km);
    }
    if((Ma != 0)&&(use_earthbyte_reconstruction))
        read_IC(Name_SurfaceTemperature_File, t.x[im-1], jm, km);  // reconstructed temperature in °C
    read_IC(Name_SurfaceNASATemperature_File, temperature_NASA.y, jm, km);
    read_IC(Name_SurfacePrecipitation_File, precipitation_NASA.y, jm, km);
    cout << endl << "      OGCM: read_Hydrosphere_Surface_Data ended ................." << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::HydrosphereDataTransfer(const string &Name_Bathymetry_File){
cout << endl << "      OGCM: AtmosphereDataTransfer" << endl;
    ifstream Transfer_File;
    string Name_v_w_Transfer_File = output_path + "/[" + Name_Bathymetry_File + "]_Transfer_Atm.vw";
    Transfer_File.precision(4);
    Transfer_File.setf(ios::fixed);
    Transfer_File.open(Name_v_w_Transfer_File);
    if(!Transfer_File.is_open()){
        cout << "ERROR: transfer file name in hydrosphere: " << Name_v_w_Transfer_File << "\n";
        cerr << "ERROR: could not open transfer file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Transfer_File >> v.x[im-1][j][k];  //non-dimensional
            Transfer_File >> w.x[im-1][j][k];  //non-dimensional
            Transfer_File >> t.x[im-1][j][k];  //non-dimensional
            Transfer_File >> p_stat.x[im-1][j][k];  //dimensional in hPa
            Transfer_File >> Evaporation_Dalton.y[j][ k];  //dimensional in mm/d
            Transfer_File >> Precipitation.y[j][k];  //dimensional in mm/d
        }
    }
    Transfer_File.close();
    cout << "      OGCM: AtmosphereDataTransfer ended" << endl;
    return;
}
/*
*
*/
void cHydrosphereModel::HydrospherePlotData(const string &Name_Bathymetry_File, int iter_cnt){
cout << endl << "      OGCM: HydrospherePlotData" << endl;
    string Name_PlotData_File = output_path + "/[" + Name_Bathymetry_File + "]_PlotData_Hyd.xyz";
    ofstream PlotData_File(Name_PlotData_File);
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);
    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    PlotData_File << "lons(deg)" 
        << ", " << "lats(deg)" 
        << ", " << "topography(m)" 
        << ", " << "v-velocity(m/s)" 
        << ", " << "w-velocity(m/s)" 
        << ", " << "velocity-mag(m/s)" 
        << ", " << "temperature(Celsius)" 
        << ", " << "salinity(g/kg)" 
        << ", " << "Ekman_pumping(cm/d)" 
        << ", " << "upwelling(cm/d)" 
        << ", " << "downwelling(cm/d)"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double vel_mag = sqrt(pow(v.x[im-1][j][k] * u_0, 2) + pow(w.x[im-1][j][k] * u_0, 2));
            PlotData_File << k << " " << 90-j 
                << " " << h.x[im-1][j][k] 
                << " " << v.x[im-1][j][k] * u_0
                << " " << w.x[im-1][j][k] * u_0 
                << " " << vel_mag
                << " " << t.x[im-1][j][k] * t_0 - t_0
                << " " << c.x[im-1][j][k] * c_0 
                << " " << EkmanPumping.y[j][k]
                << " " << Upwelling.y[j][k]
                << " " << Downwelling.y[j][k] << " " << endl;
        }
    }
    cout << "      OGCM: HydrospherePlotData ended" << endl;
    return;
}
