/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include <string>
#include <fstream>

#include "Array.h"
#include "Array_2D.h"
#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

namespace ParaViewAtm{
    void dump_array(const string &name, Array &a, double multiplier, ofstream &f){
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for(int k = 0; k < a.km; k++){
            for(int j = 0; j < a.jm; j++){
                for(int i = 0; i < a.im; i++){
                    f << (a.x[i][j][k] * multiplier) << endl;
                }
                f << "\n";
            }
            f << "\n";
        }
        f << "\n";
        f << "    </DataArray>\n";
    }
/*
 * 
*/
    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << (a.y[j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f){
        f <<  "SCALARS " << desc << " float " << 1 << endl;
        f <<  "LOOKUP_TABLE default" << endl;
        for(int i = 0; i < a.im; i++){
            for(int j = 0; j < a.jm; j++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for(int i = 0; i < a.im; i++){
            for(int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
}
/*
 * 
*/
void cAtmosphereModel::paraview_panorama_vts(string &Name_Bathymetry_File, int n){
    using namespace ParaViewAtm;
    double x, y, z, dx, dy, dz;
    string Atmosphere_panorama_vts_File_Name = output_path + "/[" 
        + Name_Bathymetry_File + "]_Atm_panorama_" + std::to_string(n) + ".vts";
    ofstream Atmosphere_panorama_vts_File;
    Atmosphere_panorama_vts_File.precision(4);
    Atmosphere_panorama_vts_File.setf(ios::fixed);
    Atmosphere_panorama_vts_File.open(Atmosphere_panorama_vts_File_Name);
    if(!Atmosphere_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open panorama_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Atmosphere_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Atmosphere_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography u-component v-component w-component Temperature CondensationTemp EvaporationTemp Epsilon_3D PressureDynamic PressureStatic WaterVapour CloudWater CloudIce CO2-Concentration Q_Latent Rain RainSuper Ice PrecipitationRain PrecipitationSnow PrecipitationConv Updraft Downdraft\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Atmosphere_panorama_vts_File << u.x[i][j][k] << " " 
                    << v.x[i][j][k] << " " << w.x[i][j][k] << endl;
            }
            Atmosphere_panorama_vts_File <<  "\n"  << endl;
        }
        Atmosphere_panorama_vts_File <<  "\n"  << endl;
    }
    Atmosphere_panorama_vts_File <<  "\n"  << endl;
    Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;
    dump_array("Topography", h, 1.0, Atmosphere_panorama_vts_File);
    dump_array("u-component", u, u_0, Atmosphere_panorama_vts_File);
    dump_array("v-component", v, u_0, Atmosphere_panorama_vts_File);
    dump_array("w-component", w, u_0, Atmosphere_panorama_vts_File);
    Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Atmosphere_panorama_vts_File << t.x[i][j][k] * t_0 - t_0 << endl;
            }
            Atmosphere_panorama_vts_File <<  "\n"  << endl;
        }
        Atmosphere_panorama_vts_File <<  "\n"  << endl;
    }
    Atmosphere_panorama_vts_File <<  "\n"  << endl;
    Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;
    dump_array("Epsilon_3D", epsilon, 1.0, Atmosphere_panorama_vts_File);
    dump_array("WaterVapour", c, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("CloudWater", cloud, 1e6, Atmosphere_panorama_vts_File);
    dump_array("CloudIce", ice, 1e6, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationRain", P_rain, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationSnow", P_snow, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationConv", P_conv, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("PressureDynamic", p_dyn, 1.0, Atmosphere_panorama_vts_File);
    dump_array("PressureStatic", p_stat, 1.0, Atmosphere_panorama_vts_File);
    dump_array("BuoyancyForce", BuoyancyForce, 1.0, Atmosphere_panorama_vts_File);
    dump_array("CoriolisForce", CoriolisForce, 1.0, Atmosphere_panorama_vts_File);
    dump_array("PressureGradientForce", PressureGradientForce, 1., Atmosphere_panorama_vts_File);
    dump_array("CO2-Concentration", co2, co2_0, Atmosphere_panorama_vts_File);
    dump_array("Q_Latent", Q_Latent, 1.0, Atmosphere_panorama_vts_File);
    dump_array("Q_Radiation", radiation, 1.0, Atmosphere_panorama_vts_File);
    dump_array("Q_Sensible", Q_Sensible, 1.0, Atmosphere_panorama_vts_File);
//    dump_array("CloudBase", i_Base, 1, Atmosphere_panorama_vts_File);
//    dump_array("CloudLFS", i_LFS, 1, Atmosphere_panorama_vts_File);
    dump_array("u_u", u_u, 1.0, Atmosphere_panorama_vts_File);
    dump_array("u_d", u_d, 1.0, Atmosphere_panorama_vts_File);
    dump_array("M_u", M_u, 1.0, Atmosphere_panorama_vts_File);
    dump_array("M_d", M_d, 1.0, Atmosphere_panorama_vts_File);
//    dump_array("q_v_u", q_v_u, 1000.0, Atmosphere_panorama_vts_File);
//    dump_array("q_v_d", q_v_d, 1000.0, Atmosphere_panorama_vts_File);
//    dump_array("q_c_u", q_c_u, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("g_p", g_p, 1000.0, Atmosphere_panorama_vts_File);
//    dump_array("c_u", c_u, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("e_d", e_d, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("e_p", e_p, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationConv", P_conv, 8.64e4, Atmosphere_panorama_vts_File);
    dump_array("MassStreamfunction", stream, 1.0, Atmosphere_panorama_vts_File);

    Atmosphere_panorama_vts_File <<  "   </PointData>\n" << endl;
    Atmosphere_panorama_vts_File <<  "   <Points>\n"  << endl;
    Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    x = 0.;
    y = 0.;
    z = 0.;
    dx = .1;
    dy = .1;
    dz = .1;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if (k == 0 || j == 0) x = 0.;
                else x = x + dx;
                Atmosphere_panorama_vts_File << x << " " << y << " " << z  << endl;
            }
            x = 0;
            y = y + dy;
            Atmosphere_panorama_vts_File <<  "\n"  << endl;
        }
        y = 0;
        z = z + dz;
        Atmosphere_panorama_vts_File <<  "\n"  << endl;
    }
    Atmosphere_panorama_vts_File <<  "    </DataArray>\n"  << endl;
    Atmosphere_panorama_vts_File <<  "   </Points>\n"  << endl;
    Atmosphere_panorama_vts_File <<  "  </Piece>\n"  << endl;
    Atmosphere_panorama_vts_File <<  " </StructuredGrid>\n"  << endl;
    Atmosphere_panorama_vts_File <<  "</VTKFile>\n"  << endl;
    Atmosphere_panorama_vts_File.close();
    cout << "   File:  " << "[" 
        + Name_Bathymetry_File + "]_Atm_panorama_" + std::to_string(n) + ".vts" 
        << "  has been written to Directory:  " 
        << output_path << endl;
    return;
}
/*
 * 
*/
void cAtmosphereModel::paraview_vtk_radial(string &Name_Bathymetry_File, 
    int Ma, int i_radial, int n){
    using namespace ParaViewAtm;
    double x, y, z, dx, dy;
    string Atmosphere_radial_File_Name = output_path + "/[" 
        + Name_Bathymetry_File + "]_Atm_radial_" + std::to_string(i_radial) 
        + "_" + std::to_string(n) + ".vtk";
    ofstream Atmosphere_vtk_radial_File;
    Atmosphere_vtk_radial_File.precision(4);
    Atmosphere_vtk_radial_File.setf(ios::fixed);
    Atmosphere_vtk_radial_File.open(Atmosphere_radial_File_Name);
    if(!Atmosphere_vtk_radial_File.is_open()){
        cerr << "ERROR: could not open paraview_vtk file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Atmosphere_vtk_radial_File <<  "# vtk DataFile Version 3.0" << endl;
    Atmosphere_vtk_radial_File <<  "Radial_Data_Atmosphere_Circulation\n";
    Atmosphere_vtk_radial_File <<  "ASCII" << endl;
    Atmosphere_vtk_radial_File <<  "DATASET STRUCTURED_GRID" << endl;
    Atmosphere_vtk_radial_File <<  "DIMENSIONS " << km << " "<< jm << " " << 1 << endl;
    Atmosphere_vtk_radial_File <<  "POINTS " << jm * km << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.1;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(k == 0) y = 0.;
            else y = y + dy;
            Atmosphere_vtk_radial_File << x << " " << y << " "<< z << endl;
        }
        y = 0.;
        x = x + dx;
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
//            temp_reconst.y[j][k] = temp_reconst.y[j][k] * t_0 - t_0;  // in °C
            double Evaporation = Evaporation_Dalton.y[j][k] 
                + Evaporation_Penman.y[j][k];
            aux_w.x[0][j][k] = Evaporation
                - 8.64e4 * Precipitation.y[j][k];
        }
    }
    Atmosphere_vtk_radial_File <<  "POINT_DATA " << jm * km << endl;
    dump_radial("u-Component", u, u_0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("v-Component", v, u_0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("w-Component", w, u_0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("rhs_u", u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("rhs_v", v, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("rhs_w", w, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("u_u", u_u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("v_u", v_u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("w_u", w_u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("u_d", u_d, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("v_d", v_d, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("w_d", w_d, 1.0, i_radial, Atmosphere_vtk_radial_File);

    Atmosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_radial_File << t.x[i_radial][j][k] * t_0 
                - t_0 << endl;
        }
    }
    dump_radial_2d("Temperature_NASA", temperature_NASA, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Temperature_Reconst", temp_reconst, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Temperature_Landscape", temp_landscape, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("PressureStatic_Landscape", p_stat_landscape, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("r_dry_Landscape", r_dry_landscape, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("r_humid_Landscape", r_humid_landscape, 1.0, Atmosphere_vtk_radial_File);
    dump_radial("Topography", h, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Topography_m", Topography, 1.0, Atmosphere_vtk_radial_File);
    dump_radial("WaterVapour", c, 1000.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudWater", cloud, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudIce", ice, 1e3, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationRain", P_rain, 8.64e4, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationSnow", P_snow, 8.64e4, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationConv", P_conv, 8.64e4, i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Precipitation", Precipitation, 8.64e4, Atmosphere_vtk_radial_File);
    dump_radial_2d("PrecipitableWater", precipitable_water, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Precipitation_NASA", precipitation_NASA, 1.0, Atmosphere_vtk_radial_File);
    dump_radial("PressureDynamic", p_dyn, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PressureStatic", p_stat, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("r_dry", r_dry, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("r_humid", r_humid, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Epsilon", epsilon, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("albedo_2D", albedo, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("epsilon_2D", epsilon_2D, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_radiation_2D", Q_radiation, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_bottom_2D", Q_bottom, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_latent_2D", Q_latent, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_sensible_2D", Q_sensible, 1.0, Atmosphere_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CoriolisForce", CoriolisForce, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PressureGradientForce", PressureGradientForce, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Q_Radiation", radiation, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Q_Latent", Q_Latent, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Q_Sensible", Q_Sensible, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("TempStandard", TempStand, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("TempDewPoint", TempDewPoint, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Humidity_rel", HumidityRel, 0.1, i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Evaporation_Dalton", Evaporation_Dalton, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Evaporation_Penman", Evaporation_Penman, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("Heat_Evaporation", Q_Evaporation, 1.0, Atmosphere_vtk_radial_File);
    dump_radial("Evap-Precip", aux_w, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Vegetation", Vegetation, 1.0, Atmosphere_vtk_radial_File);
    dump_radial("CO2-Concentration", co2, co2_0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("M_updraft", M_u, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("M_downdraft", M_d, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("CloudBase", i_Base, 1.0, Atmosphere_vtk_radial_File);
    dump_radial_2d("CloudLFS", i_LFS, 1, Atmosphere_vtk_radial_File);
//    dump_radial("c_u", c_u, 1000.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("e_d", e_d, 1000.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("q_v_d", q_v_d, 1000.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("MassStreamfunction", stream, 1.0, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("u_Streamfunction", u_stream, 1.0, i_radial, Atmosphere_vtk_radial_File);

    Atmosphere_vtk_radial_File <<  "VECTORS v-w-Cell float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_radial_File << v.x[i_radial][j][k] 
                << " " << w.x[i_radial][j][k] << " " << z << endl;
        }
    }

    Atmosphere_vtk_radial_File <<  "VECTORS v-w-Updraft float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_radial_File << v_u.x[i_radial][j][k] 
                << " " << w_u.x[i_radial][j][k] << " " << z << endl;
        }
    }

    Atmosphere_vtk_radial_File <<  "VECTORS v-w-Downdraft float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_radial_File << v_d.x[i_radial][j][k] 
                << " " << w_d.x[i_radial][j][k] << " " << z << endl;
        }
    }

    Atmosphere_vtk_radial_File.close();
    cout << "   File:  " << "[" 
        + Name_Bathymetry_File + "]_Atm_radial_" + std::to_string(i_radial) 
        + "_" + std::to_string(n) + ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
/*
 * 
*/
void cAtmosphereModel::paraview_vtk_zonal(string &Name_Bathymetry_File, 
    int k_zonal, int n){
    using namespace ParaViewAtm;
    double x, y, z, dx, dy;
    string Atmosphere_zonal_File_Name = output_path + "/[" + Name_Bathymetry_File 
        + "]_Atm_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk";
    ofstream Atmosphere_vtk_zonal_File;
    Atmosphere_vtk_zonal_File.precision(4);
    Atmosphere_vtk_zonal_File.setf(ios::fixed);
    Atmosphere_vtk_zonal_File.open(Atmosphere_zonal_File_Name);
    if(!Atmosphere_vtk_zonal_File.is_open()){
        cerr << "ERROR: could not open vtk_zonal file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Atmosphere_vtk_zonal_File <<  "# vtk DataFile Version 3.0" << endl;
    Atmosphere_vtk_zonal_File <<  "Zonal_Data_Atmosphere_Circulation\n";
    Atmosphere_vtk_zonal_File <<  "ASCII" << endl;
    Atmosphere_vtk_zonal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Atmosphere_vtk_zonal_File <<  "DIMENSIONS " << jm << " "<< im << " " << 1 << endl;
    Atmosphere_vtk_zonal_File <<  "POINTS " << im * jm << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.05;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            if (j == 0) y = 0.;
            else y = y + dy;
            Atmosphere_vtk_zonal_File << x << " " << y << " "<< z << endl;
        }
        y = 0.0;
        x = x + dx;
    }
    Atmosphere_vtk_zonal_File <<  "POINT_DATA " << im * jm << endl;
    dump_zonal("u-Component", u, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("rhs_u", u, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("rhs_v", v, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("rhs_w", w, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    Atmosphere_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
    Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Atmosphere_vtk_zonal_File << t.x[i][j][k_zonal] * t_0 
                - t_0 << endl;
        }
    }
    for(int i = 0; i <= im-1; i++){
        for(int j = 0; j < jm; j++){
            double height = get_layer_height(i);
            aux_w.x[i][j][k_zonal] = c.x[i][j][k_zonal] 
                + cloud.x[i][j][k_zonal] + ice.x[i][j][k_zonal]
                - MC_q.x[i][j][k_zonal];
            aux_u.x[i][j][k_zonal] = i_Base.y[j][k_zonal];
            aux_v.x[i][j][k_zonal] = i_LFS.y[j][k_zonal];
            aux_t.x[i][j][k_zonal] = height;
        }
    }
    dump_zonal("Topography", h, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("i_Base", aux_u, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("i_LFS", aux_v, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("WaterVapour", c, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("TempStandard", TempStand, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("TempDewPoint", TempDewPoint, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Humidity_rel", HumidityRel, 0.1, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudWater", cloud, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudIce", ice, 1e3, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Cloudiness", cloudiness, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Total_Water", aux_w, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationRain", P_rain, 8.64e4, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationSnow", P_snow, 8.64e4, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationConv", P_conv, 8.64e4, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_WaterVapour", S_v, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_CloudWater", S_c, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_CloudIce", S_i, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_Rain", S_r, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_Snow", S_s, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_CloudWater_CondEvap", S_c_c, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("r_dry", r_dry, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("r_humid", r_humid, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("M_u", M_u, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("M_d", M_d, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("MC_t", MC_t, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("MC_q", MC_q, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("MC_v", MC_v, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("MC_w", MC_w, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("g_p", g_p, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("e_d", e_d, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("e_p", e_p, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("e_l", e_l, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("c_u", c_u, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("q_c_u", q_c_u, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("q_v_u", q_v_u, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("q_v_d", q_v_d, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("s", s, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("s_u", s_u, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("s_d", s_d, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("u_u", u_u, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("u_d", u_d, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("v_u", v_u, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("v_d", v_d, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("w_u", w_u, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("w_d", w_d, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("E_u", E_u, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("D_u", D_u, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("E_d", E_d, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("D_d", D_d, 1000.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CoriolisForce", CoriolisForce, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PressureGradientForce", PressureGradientForce, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Q_Radiation", radiation, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Epsilon", epsilon, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Q_Latent", Q_Latent, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Q_Sensible", Q_Sensible, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PressureDynamic", p_dyn, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PressureStatic", p_stat, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CO2-Concentration", co2, co2_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("height", aux_t, 0.001, k_zonal, Atmosphere_vtk_zonal_File);
//    dump_zonal("CloudBase", i_Base, 1, k_zonal, Atmosphere_vtk_zonal_File);
//    dump_zonal("CloudLFS", i_LFS, 1, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("MassStreamfunction", stream, 1.0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("u_Streamfunction", u_stream, 1.0, k_zonal, Atmosphere_vtk_zonal_File);

    Atmosphere_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Atmosphere_vtk_zonal_File << u.x[i][j][k_zonal]/u_0 << " " 
                << v.x[i][j][k_zonal]/u_0 << " " << z << endl;
        }
    }

    Atmosphere_vtk_zonal_File <<  "VECTORS u-v-Stream float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Atmosphere_vtk_zonal_File << u_stream.x[i][j][k_zonal]/u_0 << " " 
                << v.x[i][j][k_zonal]/u_0 << " " << z << endl;
        }
    }

    Atmosphere_vtk_zonal_File <<  "VECTORS u-v-Updraft float " << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Atmosphere_vtk_zonal_File << u_u.x[i][j][k_zonal]/u_0 
                << " " << v_u.x[i][j][k_zonal]/u_0 << " " << z << endl;
        }
    }

    Atmosphere_vtk_zonal_File <<  "VECTORS u_v-Downdraft float " << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Atmosphere_vtk_zonal_File << u_d.x[i][j][k_zonal]/u_0 
                << " " << v_d.x[i][j][k_zonal]/u_0 << " " << z << endl;
        }
    }

    Atmosphere_vtk_zonal_File.close();
    cout << "   File:  " << "[" + Name_Bathymetry_File 
        + "]_Atm_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
/*
 * 
*/
void cAtmosphereModel::paraview_vtk_longal(string &Name_Bathymetry_File, 
    int j_longal, int n){
    using namespace ParaViewAtm;
    double x, y, z, dx, dz;
    string Atmosphere_longal_File_Name = output_path + "/[" + Name_Bathymetry_File 
        + "]_Atm_longal_" + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk";
    ofstream Atmosphere_vtk_longal_File;
    Atmosphere_vtk_longal_File.precision(4);
    Atmosphere_vtk_longal_File.setf(ios::fixed);
    Atmosphere_vtk_longal_File.open(Atmosphere_longal_File_Name);
    if(!Atmosphere_vtk_longal_File.is_open()){
        cerr << "ERROR: could not open vtk_longal file " << __FILE__ 
            << " at line " << __LINE__ << "\n";
        abort();
    }
    Atmosphere_vtk_longal_File <<  "# vtk DataFile Version 3.0" << endl;
    Atmosphere_vtk_longal_File <<  "Longitudinal_Data_Atmosphere_Circulation\n";
    Atmosphere_vtk_longal_File <<  "ASCII" << endl;
    Atmosphere_vtk_longal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Atmosphere_vtk_longal_File <<  "DIMENSIONS " << km << " "<< im << " " << 1 << endl;
    Atmosphere_vtk_longal_File <<  "POINTS " << im * km << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dz = 0.025;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            if(k == 0)  z = 0.0;
            else  z = z + dz;
            Atmosphere_vtk_longal_File << x << " " << y << " "<< z << endl;
        }
        z = 0.0;
        x = x + dx;
    }
    Atmosphere_vtk_longal_File <<  "POINT_DATA " << im * km << endl;
    dump_longal("u-Component", u, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("v-Component", v, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("w-Component", w, u_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("rhs_u", u, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("rhs_v", v, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("rhs_w", w, 1.0, j_longal, Atmosphere_vtk_longal_File);
    Atmosphere_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
    Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_longal_File << t.x[i][j_longal][k] 
                * t_0 - t_0 << endl;
        }
    }
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            aux_w.x[i][j_longal][k] = c.x[i][j_longal][k] 
                + cloud.x[i][j_longal][k] + ice.x[i][j_longal][k]
                - MC_q.x[i][j_longal][k];
            double height = get_layer_height(i);
            aux_t.x[i][j_longal][k] = height;
        }
    }
    dump_longal("Topography", h, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("WaterVapour", c, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudWater", cloud, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudIce", ice, 1e3, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Total_Water", aux_w, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationRain", P_rain, 8.64e4, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationSnow", P_snow, 8.64e4, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationConv", P_conv, 8.64e4, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PressureDynamic", p_dyn, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PressureStatic", p_stat, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("r_dry", r_dry, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("r_humid", r_humid, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CoriolisForce", CoriolisForce, 1.0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PressureGradientForce", PressureGradientForce, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Q_Latent", Q_Latent, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Q_Sensible", Q_Sensible, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Q_Radiation", radiation, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("TempStandard", TempStand, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("TempDewPoint", TempDewPoint, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Humidity_rel", HumidityRel, .1, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Epsilon", epsilon, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CO2-Concentration", co2, co2_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("height", aux_t, .001, j_longal, Atmosphere_vtk_longal_File);
//    dump_longal("CloudBase", i_Base, 1, j_longal, Atmosphere_vtk_longal_File);
//    dump_longal("CloudLFS", i_LFS, 1, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("M_u", M_u, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("M_d", M_d, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MC_t", MC_t, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MC_q", MC_q, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MC_v", MC_v, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MC_w", MC_w, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("g_p", g_p, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("e_d", e_d, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("e_p", e_p, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("e_l", e_l, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("c_u", c_u, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("q_c_u", q_c_u, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("q_v_u", q_v_u, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("q_v_d", q_v_d, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("s", s, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("s_u", s_u, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("s_d", s_d, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("u_u", u_u, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("u_d", u_d, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("v_u", v_u, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("v_d", v_d, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("w_u", w_u, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("w_d", w_d, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("E_u", E_u, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("D_u", D_u, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("E_d", E_d, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("D_d", D_d, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("MassStreamfunction", stream, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("u_Streamfunction", u_stream, 1., j_longal, Atmosphere_vtk_longal_File);

    Atmosphere_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_longal_File << u.x[i][j_longal][k] 
                << " " << y << " " << w.x[i][j_longal][k] << endl;
        }
    }

    Atmosphere_vtk_longal_File <<  "VECTORS u-w-Updraft float " << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_longal_File << u_u.x[i][j_longal][k] 
                << " " << w_u.x[i][j_longal][k] << " " << z << endl;
        }
    }

    Atmosphere_vtk_longal_File <<  "VECTORS u-w-Downdraft float " << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_longal_File << u_d.x[i][j_longal][k] 
                << " " << w_d.x[i][j_longal][k] << " " << z << endl;
        }
    }

    Atmosphere_vtk_longal_File.close();
    cout << "   File:  " << "[" + Name_Bathymetry_File 
        + "]_Atm_longal_" + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
