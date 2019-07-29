/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
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

using namespace std;

namespace ParaViewAtm{

    void dump_array( const string &name, Array &a, double multiplier, ofstream &f )
    {
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";

        for (int k = 0; k < a.km; k++) {
            for (int j = 0; j < a.jm; j++) {
                for (int i = 0; i < a.im; i++) {
                    f << (a.x[i][j][k] * multiplier) << endl;
                }
                f << "\n";
            }
            f << "\n";
        }
        f << "\n";
        f << "    </DataArray>\n";
    }

    void dump_radial( const string &desc, Array &a, double multiplier, int i, ofstream &f )
    {
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;

        for (int j = 0; j < a.jm; j++) {
            for (int k = 0; k < a.km; k++) {
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }

    void dump_radial_2d( const string &desc, Array_2D &a, double multiplier, ofstream &f )
    {
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;

        for (int j = 0; j < a.jm; j++) {
            for (int k = 0; k < a.km; k++) {
                f << (a.y[j][k] * multiplier) << endl;
            }
        }
    }

    void dump_zonal( const string &desc, Array &a, double multiplier, int k, ofstream &f )
    {
        f <<  "SCALARS " << desc << " float " << 1 << endl;
        f <<  "LOOKUP_TABLE default" << endl;

        for ( int i = 0; i < a.im; i++ ) {
            for ( int j = 0; j < a.jm; j++ ) {
                f << (a.x[ i ][ j ][ k ] * multiplier) << endl;
            }
        }
    }

    void dump_longal( const string &desc, Array &a, double multiplier, int j, ofstream &f )
    {
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;

        for (int i = 0; i < a.im; i++) {
            for (int k = 0; k < a.km; k++) {
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }

}


void cAtmosphereModel::paraview_panorama_vts (string &Name_Bathymetry_File, int n)
{
using namespace ParaViewAtm;

    double x, y, z, dx, dy, dz;

    string Atmosphere_panorama_vts_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Atm_panorama_" + std::to_string(n) + ".vts";
    ofstream Atmosphere_panorama_vts_File;
    Atmosphere_panorama_vts_File.precision(4);
    Atmosphere_panorama_vts_File.setf(ios::fixed);
    Atmosphere_panorama_vts_File.open(Atmosphere_panorama_vts_File_Name);

    if (!Atmosphere_panorama_vts_File.is_open())
    {
        cerr << "ERROR: could not open panorama_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    Atmosphere_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Atmosphere_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

    Atmosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography u-component v-component w-component Temperature CondensationTemp EvaporationTemp Epsilon_3D PressureDynamic PressureStatic WaterVapour CloudWater CloudIce CO2-Concentration Q_Latent Rain RainSuper Ice PrecipitationRain PrecipitationSnow PrecipitationConv Updraft Downdraft\">\n"  << endl;

// writing u, v und w velocity components in cartesian coordinates
    Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;


    for ( int k = 0; k < km; k++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int i = 0; i < im; i++ )
            {
// transformtion from spherical to cartesian coordinates for representation in ParaView
                Atmosphere_panorama_vts_File << u.x[ i ][ j ][ k ] << " " << v.x[ i ][ j ][ k ] << " " << w.x[ i ][ j ][ k ] << endl;
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

// writing of temperature
    Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for ( int k = 0; k < km; k++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int i = 0; i < im; i++ )
            {
                Atmosphere_panorama_vts_File << t.x[ i ][ j ][ k ] * t_0 - t_0 << endl;
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

// writing cartesian coordinates
    x = 0.;
    y = 0.;
    z = 0.;

    dx = .1;
    dy = .1;
    dz = .1;

    for ( int k = 0; k < km; k++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int i = 0; i < im; i++ )
            {
                if ( k == 0 || j == 0 ) x = 0.;
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

