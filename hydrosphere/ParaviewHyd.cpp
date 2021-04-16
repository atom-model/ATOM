/*
 * Ocean General Circulation Modell(OGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/
#include <string>
#include <fstream>
#include "Array.h"
#include "Array_2D.h"
#include "cHydrosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

namespace ParaViewHyd{
    void dump_array(const string &name, Array &a, double multiplier, ofstream &f){
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for(int k = 0; k < a.km; k++){
            for(int j = 0; j < a.jm; j++){
                for(int i = 0; i < a.im; i++){
                    f <<(a.x[i][j][k] * multiplier) << endl;
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
                f <<(a.x[i][j][k] * multiplier) << endl;
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
                f <<(a.y[j][k] * multiplier) << endl;
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
                f <<(a.x[i][j][k] * multiplier) << endl;
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
                f <<(a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
}
/*
*
*/
void cHydrosphereModel::paraview_vts(const string &Name_Bathymetry_File, int n){
    using namespace ParaViewHyd;
    double x, y, z, sinthe, sinphi, costhe, cosphi;
    string Hydrosphere_vts_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Hyd" + std::to_string(n) + ".vts";
    ofstream Hydrosphere_vts_File;
    Hydrosphere_vts_File.precision(4);
    Hydrosphere_vts_File.setf(ios::fixed);
    string path = output_path + Name_Bathymetry_File;
    Hydrosphere_vts_File.open(path);
    if(!Hydrosphere_vts_File.is_open()){
        cerr << "ERROR: could not open paraview_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Hydrosphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Hydrosphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Hydrosphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Hydrosphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Hydrosphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Bathymetry Temperature PressureDyn Salinity\">\n"  << endl;
    Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        sinphi = sin(phi.z[k]);
        cosphi = cos(phi.z[k]);
        for(int j = 0; j < jm; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            for(int i = 0; i < im; i++){
                aux_u.x[i][j][k] = sinthe * cosphi * u.x[i][j][k] + costhe * cosphi *
                    v.x[i][j][k] - sinphi * w.x[i][j][k];
                aux_v.x[i][j][k] = sinthe * sinphi * u.x[i][j][k] + sinphi * costhe *
                    v.x[i][j][k] + cosphi * w.x[i][j][k];
                aux_w.x[i][j][k] = costhe * u.x[i][j][k] - sinthe * v.x[i][j][k];

                Hydrosphere_vts_File << aux_u.x[i][j][k] << " " << aux_v.x[i][j][k]
                    << " " << aux_w.x[i][j][k]  << endl;
            }
            Hydrosphere_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_vts_File <<  "\n"  << endl;
    }
    Hydrosphere_vts_File <<  "\n"  << endl;
    Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;
    Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Bathymetry\" format=\"ascii\">\n"  << endl;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Hydrosphere_vts_File << h.x[i][j][k]  << endl;
            }
            Hydrosphere_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_vts_File <<  "\n"  << endl;
    }
    Hydrosphere_vts_File <<  "\n"  << endl;
    Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;
    Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Hydrosphere_vts_File << t.x[i][j][k]  << endl;
            }
            Hydrosphere_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_vts_File <<  "\n"  << endl;
    }
    Hydrosphere_vts_File <<  "\n"  << endl;
    Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;
    Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureDyn\" format=\"ascii\">\n"  << endl;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Hydrosphere_vts_File << p_dyn.x[i][j][k]  << endl;
            }
            Hydrosphere_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_vts_File <<  "\n"  << endl;
    }
    Hydrosphere_vts_File <<  "\n"  << endl;
    Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;
    Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Salinity\" format=\"ascii\">\n"  << endl;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Hydrosphere_vts_File << c.x[i][j][k]  << endl;
            }
            Hydrosphere_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_vts_File <<  "\n"  << endl;
    }
    Hydrosphere_vts_File <<  "\n"  << endl;
    Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;
    Hydrosphere_vts_File <<  "   </PointData>\n" << endl;
    Hydrosphere_vts_File <<  "   <Points>\n"  << endl;
    Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                x = rad.z[i] * sin(the.z[j]) * cos(phi.z[k]);
                y = rad.z[i] * sin(the.z[j]) * sin(phi.z[k]);
                z = rad.z[i] * cos(the.z[j]);

                Hydrosphere_vts_File << x << " " << y << " " << z  << endl;
            }
            Hydrosphere_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_vts_File <<  "\n"  << endl;
    }
    Hydrosphere_vts_File <<  "    </DataArray>\n"  << endl;
    Hydrosphere_vts_File <<  "   </Points>\n"  << endl;
    Hydrosphere_vts_File <<  "  </Piece>\n"  << endl;
    Hydrosphere_vts_File <<  " </StructuredGrid>\n"  << endl;
    Hydrosphere_vts_File <<  "</VTKFile>\n"  << endl;
    Hydrosphere_vts_File.close();
}
/*
*
*/
void cHydrosphereModel::paraview_panorama_vts(const string &Name_Bathymetry_File, int n){
    using namespace ParaViewHyd;
    double x, y, z, dx, dy, dz;
    string Atmosphere_panorama_vts_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Hyd_panorama_" + std::to_string(n) + ".vts";
    ofstream Hydrosphere_panorama_vts_File;
    Hydrosphere_panorama_vts_File.precision(4);
    Hydrosphere_panorama_vts_File.setf(ios::fixed);
    Hydrosphere_panorama_vts_File.open(Atmosphere_panorama_vts_File_Name);
    if(!Hydrosphere_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open panorama_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Hydrosphere_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Hydrosphere_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Bathymetry Temperature u-velocity v-velocity w-velocity PressureDynamic PressureStatic Salinity DensityWater DensitySaltWater SaltFinger SaltDiffusion SaltBalance  BuoyancyForce\">\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Hydrosphere_panorama_vts_File << u.x[i][j][k] << " " << v.x[i][j][k] << " " << w.x[i][j][k]  << endl;
            }
            Hydrosphere_panorama_vts_File <<  "\n"  << endl;
        }
        Hydrosphere_panorama_vts_File <<  "\n"  << endl;
    }
    Hydrosphere_panorama_vts_File <<  "\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;
    dump_array("Bathymetry", h, 1., Hydrosphere_panorama_vts_File);
    dump_array("u-velocity", u, u_0, Hydrosphere_panorama_vts_File);
    dump_array("v-velocity", v, u_0, Hydrosphere_panorama_vts_File);
    dump_array("w-velocity", w, u_0, Hydrosphere_panorama_vts_File);
    dump_array("Temperature", t, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("PressureDynamic", p_dyn, r_0_water * u_0 * u_0 * 1e-2, Hydrosphere_panorama_vts_File);
    dump_array("PressureStatic", p_stat, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("Salinity", c, c_0, Hydrosphere_panorama_vts_File);
    dump_array("DensityWater", r_water, 1., Hydrosphere_panorama_vts_File);
    dump_array("DensitySaltWater", r_salt_water, 1., Hydrosphere_panorama_vts_File);
    dump_array("SaltFinger", Salt_Finger, 1., Hydrosphere_panorama_vts_File);
    dump_array("SaltDiffusion", Salt_Diffusion, 1., Hydrosphere_panorama_vts_File);
    dump_array("SaltBalance", Salt_Balance, 1., Hydrosphere_panorama_vts_File);
    dump_array("BuoyancyForce", BuoyancyForce, 1., Hydrosphere_panorama_vts_File);
    dump_array("CoriolisForce", CoriolisForce, 7.292e5, Hydrosphere_panorama_vts_File);
    dump_array("PressureGradientForce", PressureGradientForce, 1., Hydrosphere_panorama_vts_File);
    Hydrosphere_panorama_vts_File <<  "   </PointData>\n" << endl;
    Hydrosphere_panorama_vts_File <<  "   <Points>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    x = 0.;
    y = 0.;
    z = 0.;
    dx = .025;
    dy = .1;
    dz = .1;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(k == 0 || j == 0) x = 0.;
                else x = x + dx;
                Hydrosphere_panorama_vts_File << x << " " << y << " " << z  << endl;
            }
            x = 0;
            y = y + dy;
            Hydrosphere_panorama_vts_File <<  "\n"  << endl;
        }
        y = 0;
        z = z + dz;
        Hydrosphere_panorama_vts_File <<  "\n"  << endl;
    }
    Hydrosphere_panorama_vts_File <<  "    </DataArray>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "   </Points>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "  </Piece>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  " </StructuredGrid>\n"  << endl;
    Hydrosphere_panorama_vts_File <<  "</VTKFile>\n"  << endl;
    Hydrosphere_panorama_vts_File.close();
    cout << "   File:  " << "[" 
        + Name_Bathymetry_File + "]_Hyd_panorama_" + std::to_string(n) + ".vts" 
        << "  has been written to Directory:  " 
        << output_path << endl;
}
/*
*
*/
void cHydrosphereModel::paraview_vtk_longal(const string &Name_Bathymetry_File, 
    int j_longal, int n){
    using namespace ParaViewHyd;
    double x, y, z, dx, dz;
    string Hydrosphere_longal_File_Name = output_path + "/[" + Name_Bathymetry_File 
        + "]_Hyd_longal_" + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk";
    ofstream Hydrosphere_vtk_longal_File;
    Hydrosphere_vtk_longal_File.precision(4);
    Hydrosphere_vtk_longal_File.setf(ios::fixed);
    Hydrosphere_vtk_longal_File.open(Hydrosphere_longal_File_Name);
    if(!Hydrosphere_vtk_longal_File.is_open()){
        cerr << "ERROR: could not open vtk_longal file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Hydrosphere_vtk_longal_File <<  "# vtk DataFile Version 3.0" << endl;
    Hydrosphere_vtk_longal_File <<  "Longitudinal_Data_Hydrosphere_Circulation\n";
    Hydrosphere_vtk_longal_File <<  "ASCII" << endl;
    Hydrosphere_vtk_longal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Hydrosphere_vtk_longal_File <<  "DIMENSIONS " << km << " "<< im << " " << 1 << endl;
    Hydrosphere_vtk_longal_File <<  "POINTS " << im * km << " float" << endl;
    x = 0.;
    y = 0.;
    z = 0.;
    dx = .1;
    dz = .025;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            if(k == 0) z = 0.;
            else z = z + dz;

            Hydrosphere_vtk_longal_File << x << " " << y << " "<< z << endl;
        }
        z = 0.;
        x = x + dx;
    }
    Hydrosphere_vtk_longal_File <<  "POINT_DATA " << im * km << endl;
    dump_longal("Bathymetry", h, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("u-Component", u, u_0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("v-Component", v, u_0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("w-Component", w, u_0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("Temperature", t, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("PressureDynamic", p_dyn, r_0_water * u_0 * u_0 * 1e-2, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("PressureStatic", p_stat, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("Salinity", c, c_0, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("DensityWater", r_water, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("DensitySaltWater", r_salt_water, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("SaltFinger", Salt_Finger, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("SaltDiffusion", Salt_Diffusion, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("SaltBalance", Salt_Balance, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("CoriolisForce", CoriolisForce, 7.292e5, j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("PressureGradientForce", PressureGradientForce, 1., j_longal, Hydrosphere_vtk_longal_File);
    Hydrosphere_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_longal_File << u.x[i][j_longal][k] << " " << y << " " << w.x[i][j_longal][k] << endl;
        }
    }
    Hydrosphere_vtk_longal_File.close();
    cout << "   File:  " << "[" + Name_Bathymetry_File 
        + "]_Hyd_longal_" + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
}
/*
*
*/
void cHydrosphereModel::paraview_vtk_radial(const string &Name_Bathymetry_File, 
    int i_radial, int n){
    using namespace ParaViewHyd;
    double x, y, z, dx, dy;
    string Hydrosphere_radial_File_Name = output_path + "/[" + Name_Bathymetry_File 
        + "]_Hyd_radial_" + std::to_string(i_radial) + "_" + std::to_string(n) + ".vtk";
    ofstream Hydrosphere_vtk_radial_File;
    Hydrosphere_vtk_radial_File.precision(4);
    Hydrosphere_vtk_radial_File.setf(ios::fixed);
    Hydrosphere_vtk_radial_File.open(Hydrosphere_radial_File_Name);
    if(!Hydrosphere_vtk_radial_File.is_open()){
        cerr << "ERROR: could not open paraview_vtk file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Hydrosphere_vtk_radial_File <<  "# vtk DataFile Version 3.0" << endl;
    Hydrosphere_vtk_radial_File <<  "Radial_Data_Hydrosphere_Circulation\n";
    Hydrosphere_vtk_radial_File <<  "ASCII" << endl;
    Hydrosphere_vtk_radial_File <<  "DATASET STRUCTURED_GRID" << endl;
    Hydrosphere_vtk_radial_File <<  "DIMENSIONS " << km << " "<< jm << " " << 1 << endl;
    Hydrosphere_vtk_radial_File <<  "POINTS " << jm * km << " float" << endl;
    x = 0.;
    y = 0.;
    z = 0.;
    dx = .1;
    dy = .1;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(k == 0) y = 0.;
            else y = y + dy;

            Hydrosphere_vtk_radial_File << x << " " << y << " "<< z << endl;
        }
        y = 0.;
        x = x + dx;
    }
    Hydrosphere_vtk_radial_File <<  "POINT_DATA " << jm * km << endl;
    Hydrosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_radial_File << t.x[i_radial][j][k] * t_0 - t_0 << endl;
            aux_v.x[i_radial][j][k] = Evaporation_Dalton.y[j][k] - Precipitation.y[j][k];
            if(is_land(h, i_radial, j, k))  aux_v.x[i_radial][j][k] = 0.;
        }
    }
    dump_radial("Bathymetry", h, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial_2d("Bathymetry_m", Bathymetry, 1., Hydrosphere_vtk_radial_File);
    dump_radial("u-Component", u, u_0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("v-Component", v, u_0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("w-Component", w, u_0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("PressureDynamic", p_dyn, r_0_water * u_0 * u_0 * 1e-2, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("PressureStatic", p_stat, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("Salinity", c, c_0, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("DensityWater", r_water, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("DensitySaltWater", r_salt_water, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("SaltFinger", Salt_Finger, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("SaltDiffusion", Salt_Diffusion, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("SaltBalance", Salt_Balance, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("CoriolisForce", CoriolisForce, 7.292e5, i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("PressureGradientForce", PressureGradientForce, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial_2d("EkmanPumping", EkmanPumping, 1., Hydrosphere_vtk_radial_File);
    dump_radial_2d("Upwelling", Upwelling, 1., Hydrosphere_vtk_radial_File);
    dump_radial_2d("Downwelling", Downwelling, 1., Hydrosphere_vtk_radial_File);
    dump_radial_2d("Evaporation_Dalton", Evaporation_Dalton, 1., Hydrosphere_vtk_radial_File);
    dump_radial_2d("Precipitation", Precipitation, 1., Hydrosphere_vtk_radial_File);
    dump_radial_2d("SalinityEvapPrec", salinity_evaporation, 1e8 * c_0, Hydrosphere_vtk_radial_File);
    dump_radial("Evap-Precip", aux_v, 1., i_radial, Hydrosphere_vtk_radial_File);
    Hydrosphere_vtk_radial_File <<  "VECTORS v-w-Cell float" << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Hydrosphere_vtk_radial_File << v.x[i_radial][j][k] << " " << w.x[i_radial][j][k] << " " << z << endl;
        }
    }
    Hydrosphere_vtk_radial_File.close();
    cout << "   File:  " << "[" 
        + Name_Bathymetry_File + "]_Hyd_radial_" + std::to_string(i_radial) 
        + "_" + std::to_string(n) + ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
}
/*
*
*/
void cHydrosphereModel::paraview_vtk_zonal(const string &Name_Bathymetry_File, 
    int k_zonal, int n){
    using namespace ParaViewHyd;
    double x, y, z, dx, dy;
    string Hydrosphere_zongal_File_Name = output_path + "/[" + Name_Bathymetry_File 
        + "]_Hyd_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk";
    ofstream Hydrosphere_vtk_zonal_File;
    Hydrosphere_vtk_zonal_File.precision(4);
    Hydrosphere_vtk_zonal_File.setf(ios::fixed);
    Hydrosphere_vtk_zonal_File.open(Hydrosphere_zongal_File_Name);
    if(!Hydrosphere_vtk_zonal_File.is_open()){
        cerr << "ERROR: could not open vtk_zonal file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Hydrosphere_vtk_zonal_File <<  "# vtk DataFile Version 3.0" << endl;
    Hydrosphere_vtk_zonal_File <<  "Zonal_Data_Hydrosphere_Circulation\n";
    Hydrosphere_vtk_zonal_File <<  "ASCII" << endl;
    Hydrosphere_vtk_zonal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Hydrosphere_vtk_zonal_File <<  "DIMENSIONS " << jm << " "<< im << " " << 1 << endl;
    Hydrosphere_vtk_zonal_File <<  "POINTS " << im * jm << " float" << endl;
    x = 0.;
    y = 0.;
    z = 0.;
    dx = .1;
    dy = .05;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            if(j == 0) y = 0.;
            else y = y + dy;
            Hydrosphere_vtk_zonal_File << x << " " << y << " "<< z << endl;
        }
        y = 0.;
        x = x + dx;
    }
    Hydrosphere_vtk_zonal_File <<  "POINT_DATA " << im * jm << endl;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                aux_u.x[i][j][k] = t.x[i][j][k] * t_0 - t_0;
            }
        }
    }
    dump_zonal("Bathymetry", h, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("u-Component", u, u_0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("Temperature", aux_u, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("PressureDynamic", p_dyn, r_0_water * u_0 * u_0 * 1e-2, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("PressureStatic", p_stat, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("Salinity", c, c_0, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("DensityWater", r_water, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("DensitySaltWater", r_salt_water, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("SaltFinger", Salt_Finger, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("SaltDiffusion", Salt_Diffusion, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("SaltBalance", Salt_Balance, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("CoriolisForce", CoriolisForce, 7.292e5, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("PressureGradientForce", PressureGradientForce, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    Hydrosphere_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Hydrosphere_vtk_zonal_File << u.x[i][j][k_zonal] << " " << v.x[i][j][k_zonal] << " " << z << endl;
        }
    }
    Hydrosphere_vtk_zonal_File.close();
    cout << "   File:  " << "[" + Name_Bathymetry_File 
        + "]_Hyd_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
}
