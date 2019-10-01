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

    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }

    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << (a.y[j][k] * multiplier) << endl;
            }
        }
    }

    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f){
        f <<  "SCALARS " << desc << " float " << 1 << endl;
        f <<  "LOOKUP_TABLE default" << endl;
        for(int i = 0; i < a.im; i++){
            for(int j = 0; j < a.jm; j++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }

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
    dump_array("u-component", u, 1.0, Atmosphere_panorama_vts_File);
    dump_array("v-component", v, 1.0, Atmosphere_panorama_vts_File);
    dump_array("w-component", w, 1.0, Atmosphere_panorama_vts_File);
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
    dump_array("Epsilon_3D", epsilon_3D, 1.0, Atmosphere_panorama_vts_File);
    dump_array("WaterVapour", c, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("CloudWater", cloud, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("CloudIce", ice, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationRain", P_rain, 86400., Atmosphere_panorama_vts_File);
    dump_array("PrecipitationSnow", P_snow, 86400., Atmosphere_panorama_vts_File);
    dump_array("PressureDynamic", p_dyn, u_0 * u_0 * r_air *.01, Atmosphere_panorama_vts_File);
    dump_array("PressureStatic", p_stat, 1.0, Atmosphere_panorama_vts_File);
    dump_array("BuoyancyForce", BuoyancyForce, 1.0, Atmosphere_panorama_vts_File);
    dump_array("CO2-Concentration", co2, 1.0, Atmosphere_panorama_vts_File);
    dump_array("Q_Latent", Q_Latent, 1.0, Atmosphere_panorama_vts_File);
    dump_array("Q_Sensible", Q_Sensible, 1.0, Atmosphere_panorama_vts_File);
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
}

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
    x = 0.;
    y = 0.;
    z = 0.;
    dx = .1;
    dy = .1;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(k == 0) y = 0.;
            else y = y + dy;
            Atmosphere_vtk_radial_File << x << " " << y << " "<< z << endl;
        }
        y = 0.;
        x = x + dx;
    }
    if(Ma == 0){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                temp_NASA.y[j][k] = temperature_NASA.y[j][k] * t_0 - t_0;
                for(int i = im-2; i >= 0; i--){
                    if(Ma == 0){
                        if((is_air(h, i+1, j, k)) && (is_land(h, i, j, k)))
                            aux_v.x[i_radial][j][k] = (t.x[0][j][k] 
                            - temperature_NASA.y[j][k]) * t_0;
                        if(is_air(h, 0, j, k))
                            aux_v.x[i_radial][j][k] = (t.x[0][j][k] 
                            - temperature_NASA.y[j][k]) * t_0;
                    }else  aux_v.x[i_radial][j][k] = 0.;
                    aux_w.x[i_radial][j][k] = Evaporation_Dalton.y[j][k] 
                        - Precipitation.y[j][k];
                }
            }
        }
    }
    float exp_pressure = g / (1.e-2 * gam * R_Air);
    float p_h, p_SL, t_u, height, E_Rain, q_Rain;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                t_u = t.x[i][j][k] * t_0; // in K
                p_SL = .01 * ( r_air * R_Air * t.x[0][j][k] * t_0 ); // given in hPa
                height = get_layer_height(i);
                if(i != 0)
                      p_h = pow((t_u - gam * height * 1.e-2) / t_u, 
                            exp_pressure) * p_SL;
                else  p_h = p_SL;
                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                q_Rain = ep * E_Rain / ( p_h - E_Rain ); // water vapour amount at saturation with water formation in kg/kg
                if(c.x[i][j][k] >= q_Rain){
                    dew_point_temperature.y[j][k] = t_u - t_0;
                    condensation_level.y[j][k] = height;
                    break;
                }
            }
        }
    }
    Atmosphere_vtk_radial_File <<  "POINT_DATA " << jm * km << endl;
    dump_radial("u-Component", u, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("v-Component", v, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("w-Component", w, 1., i_radial, Atmosphere_vtk_radial_File);
    Atmosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_radial_File << t.x[i_radial][j][k] * t_0 
                - t_0 << endl;
        }
    }
    dump_radial_2d("Temp_NASA", temp_NASA, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("DewPointTemp", dew_point_temperature, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("CondensLevel", condensation_level, 1.e-3, Atmosphere_vtk_radial_File);
    dump_radial("Temp_NASA_diff", aux_v, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Topography", h, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Topography_m", Topography, 1., Atmosphere_vtk_radial_File);
    dump_radial("WaterVapour", c, 1000., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudWater", cloud, 1000., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudIce", ice, 1000., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationRain", P_rain, 86400., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationSnow", P_snow, 86400., i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Precipitation", Precipitation, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("PrecipitableWater_2D", precipitable_water, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Precipitation_NASA", precipitation_NASA, 1., Atmosphere_vtk_radial_File);
    dump_radial("PressureDynamic", p_dyn, u_0 * u_0 * r_air *.01, i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PressureStatic", p_stat, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Epsilon_3D", epsilon_3D, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("albedo_2D", albedo, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("epsilon_2D", epsilon, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_radiation_2D", Q_radiation, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_bottom_2D", Q_bottom, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_latent_2D", Q_latent, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_sensible_2D", Q_sensible, 1., Atmosphere_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Q_Radiation", radiation_3D, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Q_Latent", Q_Latent, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Q_Sensible", Q_Sensible, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Evaporation_Dalton", Evaporation_Dalton, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Heat_Evaporation", Q_Evaporation, 1., Atmosphere_vtk_radial_File);
//    dump_radial("Evap-Precip", aux_w, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Evap-Precip", aux_w, -1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Vegetation", Vegetation, 1., Atmosphere_vtk_radial_File);
    dump_radial("CO2-Concentration", co2, co2_0, i_radial, Atmosphere_vtk_radial_File);
    Atmosphere_vtk_radial_File <<  "VECTORS v-w-Cell float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_radial_File << v.x[i_radial][j][k] 
                << " " << w.x[i_radial][j][k] << " " << z << endl;
        }
    }
    Atmosphere_vtk_radial_File.close();
}

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
    x = 0.;
    y = 0.;
    z = 0.;
    dx = .1;
    dy = .05;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            if (j == 0) y = 0.;
            else y = y + dy;
            Atmosphere_vtk_zonal_File << x << " " << y << " "<< z << endl;
        }
        y = 0.;
        x = x + dx;
    }
    Atmosphere_vtk_zonal_File <<  "POINT_DATA " << im * jm << endl;
    dump_zonal("u-Component", u, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("v-Component", v, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("w-Component", w, 1., k_zonal, Atmosphere_vtk_zonal_File);
    Atmosphere_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
    Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Atmosphere_vtk_zonal_File << t.x[i][j][k_zonal] * t_0 
                - t_0 << endl;
        }
    }
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            aux_w.x[i][j][k_zonal] = c.x[i][j][k_zonal] 
                + cloud.x[i][j][k_zonal] + ice.x[i][j][k_zonal];
            float t_u = t.x[i][j][k_zonal] * t_0;
            float T = t_u;
            float p_SL = .01 * (r_air * R_Air * t.x[0][j][k_zonal] * t_0);
            float p_h;
            if(i != 0)   p_h = exp (- g * (double)i * (L_atm 
                / (double)(im-1)) / (R_Air * t_u)) * p_SL;
            else         p_h = p_SL;
            float E_Rain = hp * exp_func(T, 17.2694, 35.86);
            float E_Ice = hp * exp_func(T, 21.8746, 7.66);
            float q_Rain = ep * E_Rain / (p_h - E_Rain);
            float q_Ice = ep * E_Ice / (p_h - E_Ice);
            aux_u.x[i][j][k_zonal] = c.x[i][j][k_zonal] - q_Rain;
            if(T <= t_0){
                aux_v.x[i][j][k_zonal] = c.x[i][j][k_zonal] - q_Ice;
            }
            else  aux_v.x[i][j][k_zonal] = 0.;
            float height = get_layer_height(i);
            BuoyancyForce.x[i][j][k_zonal] = height;
        }
    }
    for(int j = 0; j < jm; j++){
        double t_eff_tropo = t_tropopause_pole - t_tropopause;
        double temp_tropopause = t_eff_tropo * parabola(j/((jm-1)/2.0)) +
            t_tropopause_pole;   //temperature at tropopause     
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            int i_trop = get_tropopause_layer(j);
            M_u.x[0][j][k] = t.x[0][j][k];
            double t_mount_top = (temp_tropopause - M_u.x[0][j][k]) *
                (get_layer_height(i_mount) / get_layer_height(i_trop)) + 
                M_u.x[0][j][k]; //temperature at mountain top
            for(int i = i_mount; i < im; i++){
                if(i < i_trop+1){
                    if(i>i_mount){
                        // linear temperature decay up to tropopause, privat  approximation
//                        M_u.x[i][j][k] = (temp_tropopause - M_u.x[0][j][k]) * 
//                            (get_layer_height(i) / get_layer_height(i_trop)) + 
//                            M_u.x[0][j][k]; 
                        // US Standard Atmosphere
                        M_u.x[i][j][k] = M_u.x[0][j][k] - 6.5 * 
                            (get_layer_height(i) / 1000.) / t_0;
                    }else{
                        M_u.x[i][j][k] = t_mount_top; //inside mountain
                    }
                }else{ // above tropopause
                    M_u.x[i][j][k] = temp_tropopause;
                }
            }
            M_u.x[0][j][k] = t_mount_top;
        }
    }
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                M_u.x[i][j][k] = M_u.x[i][j][k] * t_0 - t_0;
            }
        }
    }
    dump_zonal("Topography", h, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("WaterVapour", c, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Temp_Standard", M_u, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudWater", cloud, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudIce", ice, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Saturation_Water", aux_u, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Saturation_Ice", aux_v, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Total_Water", aux_w, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationRain", P_rain, 86400., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationSnow", P_snow, 86400., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_WaterVapour", S_v, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_CloudWater", S_c, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_CloudIce", S_i, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_Rain", S_r, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_Snow", S_s, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_CloudWater_CondEvap", S_c_c, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Q_Radiation", radiation_3D, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Epsilon_3D", epsilon_3D, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Q_Latent", Q_Latent, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Q_Sensible", Q_Sensible, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PressureDynamic", p_dyn, u_0 * u_0 * r_air *.01, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PressureStatic", p_stat, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CO2-Concentration", co2, co2_0, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("height", BuoyancyForce, .001, k_zonal, Atmosphere_vtk_zonal_File);
    Atmosphere_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Atmosphere_vtk_zonal_File << u.x[i][j][k_zonal] << " " 
                << v.x[i][j][k_zonal] << " " << z << endl;
        }
    }
    Atmosphere_vtk_zonal_File.close();
}

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
    x = 0.;
    y = 0.;
    z = 0.;
    dx = .1;
    dz = .025;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            if(k == 0)  z = 0.;
            else  z = z + dz;
            Atmosphere_vtk_longal_File << x << " " << y << " "<< z << endl;
        }
        z = 0.;
        x = x + dx;
    }
    Atmosphere_vtk_longal_File <<  "POINT_DATA " << im * km << endl;
    dump_longal("u-Component", u, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("v-Component", v, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("w-Component", w, 1., j_longal, Atmosphere_vtk_longal_File);
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
                + cloud.x[i][j_longal][k] + ice.x[i][j_longal][k];
            float height = get_layer_height(i);
            aux_v.x[i][j_longal][k] = height;
        }
    }
    dump_longal("Topography", h, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("WaterVapour", c, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudWater", cloud, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudIce", ice, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Total_Water", aux_w, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationRain", P_rain, 86400., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationSnow", P_snow, 86400., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PressureDynamic", p_dyn, u_0 * u_0 * r_air *.01, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PressureStatic", p_stat, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Q_Latent", Q_Latent, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Q_Sensible", Q_Sensible, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Epsilon_3D", epsilon_3D, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CO2-Concentration", co2, co2_0, j_longal, Atmosphere_vtk_longal_File);
    dump_longal("height", aux_v, .001, j_longal, Atmosphere_vtk_longal_File);
    Atmosphere_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Atmosphere_vtk_longal_File << u.x[i][j_longal][k] 
                << " " << y << " " << w.x[i][j_longal][k] << endl;
        }
    }
    Atmosphere_vtk_longal_File.close();
}
