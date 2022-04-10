#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cstdlib>

#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;
 
void cAtmosphereModel::read_Atmosphere_Surface_Data(int Ma){
cout << endl << "      AGCM: read_Atmosphere_Surface_Data ......................." << endl;
    if(!is_global_temperature_curve_loaded()) 
        load_global_temperature_curve();
    if(!is_equat_temperature_curve_loaded()) 
        load_equat_temperature_curve();
    if(!is_pole_temperature_curve_loaded()) 
        load_pole_temperature_curve();
    //Prepare the temperature, precipitation and velocity data file
    string Name_SurfaceTemperature_File  = temperature_file;
    string Name_SurfaceNASATemperature_File  = temperature_file;
    string Name_SurfaceNASAPrecipitation_File = precipitation_file;
    if(Ma != 0 && use_earthbyte_reconstruction){
        Name_SurfaceTemperature_File = output_path 
            + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";  // reconstructed temperature in °C
//        Name_SurfaceNASAPrecipitation_File = output_path 
//            + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";    
        velocity_v_file = output_path + std::to_string(Ma) 
            + "Ma_Reconstructed_wind_v.xyz";
        velocity_w_file = output_path + std::to_string(Ma) 
            + "Ma_Reconstructed_wind_w.xyz";
        struct stat info;
//        char direction_i = 'i', direction_j = 'j', direction_k = 'k';
        if(stat(output_path.c_str(), &info) != 0){
             mkdir(output_path.c_str(), 0777);
        }
        if(stat(Name_SurfaceTemperature_File.c_str(), &info) != 0 || 
//           stat(Name_SurfaceNASAPrecipitation_File.c_str(), &info) != 0 ||
           stat(velocity_v_file.c_str(), &info) != 0 ||
           stat(velocity_w_file.c_str(), &info) != 0){
               std::string cmd_str = "python " + reconstruction_script_path 
                   + " " + std::to_string(Ma - time_step) + " " 
                   + std::to_string(Ma) + " " + output_path + " " 
                   + BathymetrySuffix + " atm";
               int ret = system(cmd_str.c_str());
               std::cout << " reconstruction script returned: " 
                   << ret << std::endl;
        } 
    }
    if(!has_printed_welcome_msg)  print_welcome_msg();
    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    cout << endl << "      topography given by the x-y-z data set:    " 
        << bathymetry_name.c_str();
    init_topography(bathymetry_path + "/" + bathymetry_name);
    if(use_NASA_velocity){
        read_IC(velocity_v_file, v.x[0], jm, km);
        read_IC(velocity_w_file, w.x[0], jm, km);    
    }
    if((Ma != 0)&&(use_earthbyte_reconstruction))
        read_IC(Name_SurfaceTemperature_File, t.x[0], jm, km);  // reconstructed temperature in °C
    read_IC(Name_SurfaceNASATemperature_File, temperature_NASA.y, jm, km);
    read_IC(Name_SurfaceNASAPrecipitation_File, precipitation_NASA.y, jm, km);
    cout << endl << "      AGCM: read_Atmosphere_Surface_Data ended ................." << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::AtmosphereDataTransfer(const string &Name_Bathymetry_File){
cout << endl << "      AGCM: AtmosphereDataTransfer" << endl;
    ofstream Transfer_File;
    string Name_v_w_Transfer_File = output_path + "/[" + Name_Bathymetry_File + "]_Transfer_Atm.vw";
    Transfer_File.precision(4);
    Transfer_File.setf(ios::fixed);
    Transfer_File.open(Name_v_w_Transfer_File);
    if(!Transfer_File.is_open()){
        cout << "ERROR: transfer file name in atmosphere: " << Name_v_w_Transfer_File << "\n";
        cerr << "ERROR: could not open transfer file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Transfer_File << v.x[0][j][k] << " " << w.x[0][j][k] << " " <<   //v and w non-dimensional, t non-dim,
                t.x[0][j][k] << " " << p_stat.x[0][j][k] << " " << Evaporation_Dalton.y[j][k] << 
                " " << Precipitation.y[j][k] << endl;
        }
    }
    Transfer_File.close();
    cout << "      AGCM: AtmosphereDataTransfer ended" << endl;
    return;
}
/*
*
*/
void cAtmosphereModel::AtmospherePlotData(const string &Name_Bathymetry_File, int iter_cnt){
cout << endl << "      AGCM: AtmospherePlotData" << endl;
    string Name_PlotData_File = output_path + "/[" + Name_Bathymetry_File + "]_PlotData_Atm.xyz";
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
        << ", " << "water_vapour(g/kg)" 
        << ", " << "precipitation(mm/d)" 
        << ", " << "precipitable water(mm)" 
        << ", " << "Evaporation_Penman (mm/d) " 
        << ", " << "temp_reconst (Celsius) " 
        << ", " << "temp_landscape (Celsius) " 
        << ", " << "p_stat_landscape (hPa) " 
        << ", " << "r_dry_landscape (kg/m3) " 
        << ", " << "r_humid_landscape (kg/m3) " << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            double vel_mag = sqrt(pow(v.x[0][j][k] * u_0, 2) + pow(w.x[0][j][k] * u_0, 2));
            PlotData_File << k << " " << 90-j 
                << " " << h.x[0][j][k] 
                << " " << v.x[0][j][k] * u_0
                << " " << w.x[0][j][k] * u_0 
                << " " << vel_mag
                << " " << t.x[0][j][k] * t_0 - t_0
                << " " << c.x[0][j][k] * 1000. 
                << " " << Precipitation.y[j][k] * 86400.0
                << " " << precipitable_water.y[j][k] 
                << " " << Evaporation_Penman.y[j][k]
                << " " << temp_reconst.y[j][k]
                << " " << temp_landscape.y[j][k]
                << " " << p_stat_landscape.y[j][k]
                << " " << r_dry_landscape.y[j][k]
                << " " << r_humid_landscape.y[j][k] << " " << endl;
        }
    }
//    v.printArray("AGCM", im, jm, km);
//    w.printArray("AGCM", im, jm, km);
    cout << "      AGCM: AtmospherePlotData ended" << endl;
    return;
}
/*
*
*/
