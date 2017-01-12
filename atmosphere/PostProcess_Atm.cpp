/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "PostProcess_Atm.h"

using namespace std;

PostProcess_Atmosphere::PostProcess_Atmosphere(int im, int jm, int km, string &output_path) {
    this->im = im;
    this->jm = jm;
    this->km = km;
    this->output_path = output_path;
}

PostProcess_Atmosphere::~PostProcess_Atmosphere() {}

void PostProcess_Atmosphere::Atmosphere_SequelFile_write(string &Name_Bathymetry_File, int &n, double &time, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &co2n) {
    stringstream Name_Sequel_File;

    streampos anfangpos, endpos;

    // file administration
    // FIXME: where are we supposed to read/write these sequel files?
    Name_Sequel_File << output_path << "/[" << Name_Bathymetry_File << "]_Sequel_Atm.seq";
    ofstream Sequel_File;
    Sequel_File.precision(4);
    Sequel_File.setf(ios::fixed);
    Sequel_File.open(Name_Sequel_File.str());

    if (!Sequel_File.is_open()) {
        cerr << "ERROR: could not open sequel file " << Name_Sequel_File.str() << " for writing at " << __FILE__ << " line " << __LINE__ << "\n";
        abort();
    }

    // begin writing
    cout << "***** Atmosphere_SequelFile_write:   n = " << n << "  time = " << time << endl;
    Sequel_File << n << " " << time << endl;

    for ( int i = 0; i < im; i++ ) {
        Sequel_File << rad.z[ i ] << endl;
    }

    for ( int j = 0; j < jm; j++ ) {
        Sequel_File << the.z[ j ] << endl;
    }

    for ( int k = 0; k < km; k++ ) {
        Sequel_File << phi.z[ k ] << endl;
    }

    for ( int i = 0; i < im; i++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int k = 0; k < km; k++ ) {
                Sequel_File << u.x[ i ][ j ][ k ] << " " << v.x[ i ][ j ][ k ] << " " << w.x[ i ][ j ][ k ]  << endl;
            }
        }
    }

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Sequel_File << un.x[ i ][ j ][ k ] << " " << vn.x[ i ][ j ][ k ] << " " << wn.x[ i ][ j ][ k ]  << endl;
            }
        }
    }

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Sequel_File << t.x[ i ][ j ][ k ] << " " << tn.x[ i ][ j ][ k ]  << endl;
            }
        }
    }

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Sequel_File << c.x[ i ][ j ][ k ] << " " << cn.x[ i ][ j ][ k ]  << endl;
            }
        }
    }

    // end writing
    Sequel_File.close();
}

void PostProcess_Atmosphere::Atmosphere_SequelFile_read(string &Name_Bathymetry_File, int &n, double &time, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &co2n) {
    stringstream Name_Sequel_File;

    streampos anfangpos, endpos;

    // file administration
    // FIXME: where are we supposed to read/write these sequel files?
    Name_Sequel_File << output_path << "/[" << Name_Bathymetry_File << "]_Sequel_Atm.seq";
    ifstream Sequel_File;
    Sequel_File.open(Name_Sequel_File.str());

    if (!Sequel_File.is_open()) {
        cerr << "WARNING: could not open sequel file " << Name_Sequel_File.str() << " for reading at " << __FILE__ << " line " << __LINE__ << "\n";
        return; // we tolerate it for now as we don't know what the sequel files are
    }

    // begin reading
    Sequel_File >> n;
    Sequel_File >> time;

    cout << "***** Atmosphere_SequelFile_read:   n = " << n << "  time = " << time << endl;

    for ( int i = 0; i < im; i++ ) {
        Sequel_File >> rad.z[ i ];
    }

    for ( int j = 0; j < jm; j++ ) {
        Sequel_File >> the.z[ j ];
    }

    for ( int k = 0; k < km; k++ ) {
        Sequel_File >> phi.z[ k ];
    }

    for ( int i = 0; i < im; i++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int   k = 0; k < km; k++ ) {
                Sequel_File >> u.x[ i ][ j ][ k ];
                Sequel_File >> v.x[ i ][ j ][ k ];
                Sequel_File >> w.x[ i ][ j ][ k ];
            }
        }
    }

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Sequel_File >> un.x[ i ][ j ][ k ];
                Sequel_File >> vn.x[ i ][ j ][ k ];
                Sequel_File >> wn.x[ i ][ j ][ k ];
            }
        }
    }

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Sequel_File >> t.x[ i ][ j ][ k ];
                Sequel_File >> tn.x[ i ][ j ][ k ];
            }
        }
    }

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Sequel_File >> c.x[ i ][ j ][ k ];
                Sequel_File >> cn.x[ i ][ j ][ k ];
            }
        }
    }

    // end reading
    Sequel_File.close();
}

void PostProcess_Atmosphere::Atmosphere_v_w_Transfer ( string &Name_Bathymetry_File, Array &v, Array &w, Array &p_dyn )
{
    // file administration
    string Name_v_w_Transfer_File = output_path + "/[" + Name_Bathymetry_File + "]_Transfer_Atm.vw";
    ofstream v_w_Transfer_File;
    v_w_Transfer_File.precision(4);
    v_w_Transfer_File.setf(ios::fixed);
    v_w_Transfer_File.open(Name_v_w_Transfer_File);

    if (!v_w_Transfer_File.is_open()) {
        cerr << "ERROR: could not open transfer file at " << __FILE__ << " line " << __LINE__ << "\n";
        abort();
    }

    for ( int j = 0; j < jm; j++ ) {
        for ( int k = 0; k < km; k++ ) {
            v_w_Transfer_File << v.x[ 0 ][ j ][ k ] << " " << w.x[ 0 ][ j ][ k ] << " " << p_dyn.x[ 0 ][ j ][ k ]  << endl;
        }
    }

    v_w_Transfer_File.close();
}

void PostProcess_Atmosphere::dump_array(const string &name, Array &a, double multiplier, ofstream &f) {
    f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";

    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            for (int i = 0; i < im; i++) {
                f << (a.x[i][j][k] * multiplier) << endl;
            }
            f << "\n";
        }
        f << "\n";
    }
    f << "\n";
    f << "    </DataArray>\n";
}

void PostProcess_Atmosphere::paraview_vts ( string &Name_Bathymetry_File, int &n, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &p_dyn, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Rain, Array &Ice, Array &Rain_super, Array &IceLayer ) {
    double x, y, z, sinthe, sinphi, costhe, cosphi;

    string Atmosphere_vts_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Atm" + std::to_string(n) + ".vts";
    ofstream Atmosphere_vts_File;
    Atmosphere_vts_File.precision(4);
    Atmosphere_vts_File.setf(ios::fixed);
    Atmosphere_vts_File.open(Atmosphere_vts_File_Name);
    if (!Atmosphere_vts_File.is_open()) {
        cerr << "ERROR: could not open paraview_vts file at " << __FILE__ << " line " << __LINE__ << "\n";
        abort();
    }

    Atmosphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Atmosphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Atmosphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Atmosphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

    //        Atmosphere_vts_File <<  "   <PointData Vectors=\"Velocity Rotation\" Scalars=\"Topography Temperature Pressure WaterVapour u-Component v-Component w-Component\">\n"  << endl;
    Atmosphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography Temperature Pressure WaterVapour\">\n"  << endl;

    // writing u, v und w velocity components in cartesian coordinates
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;

    for ( int k = 0; k < km; k++ ) {
        sinphi = sin( phi.z[ k ] );
        cosphi = cos( phi.z[ k ] );

        for ( int j = 0; j < jm; j++ ) {
            sinthe = sin( the.z[ j ] );
            costhe = cos( the.z[ j ] );

            for ( int i = 0; i < im; i++ ) {
                aux_u.x[ i ][ j ][ k ] = sinthe * cosphi * u.x[ i ][ j ][ k ] + costhe * cosphi * v.x[ i ][ j ][ k ] - sinphi * w.x[ i ][ j ][ k ];
                aux_v.x[ i ][ j ][ k ] = sinthe * sinphi * u.x[ i ][ j ][ k ] + sinphi * costhe * v.x[ i ][ j ][ k ] + cosphi * w.x[ i ][ j ][ k ];
                aux_w.x[ i ][ j ][ k ] = costhe * u.x[ i ][ j ][ k ] - sinthe * v.x[ i ][ j ][ k ];

                Atmosphere_vts_File << aux_u.x[ i ][ j ][ k ] << " " << aux_v.x[ i ][ j ][ k ] << " " << aux_w.x[ i ][ j ][ k ]  << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }
    Atmosphere_vts_File <<  "\n"  << endl;
    Atmosphere_vts_File <<  "    </DataArray>\n" << endl;

    // writing topography
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Topography\" format=\"ascii\">\n"  << endl;

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Atmosphere_vts_File << h.x[ i ][ j ][ k ]  << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }
    Atmosphere_vts_File <<  "\n"  << endl;
    Atmosphere_vts_File <<  "    </DataArray>\n" << endl;

    // writing of temperature
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Atmosphere_vts_File << t.x[ i ][ j ][ k ]  << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }
    Atmosphere_vts_File <<  "\n"  << endl;
    Atmosphere_vts_File <<  "    </DataArray>\n" << endl;

    // writing pressure
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">\n"  << endl;

    for ( int k = 0; k < km; k++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int i = 0; i < im; i++ )
            {
                Atmosphere_vts_File << p_dyn.x[ i ][ j ][ k ] * 100. << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }
    Atmosphere_vts_File <<  "\n"  << endl;
    Atmosphere_vts_File <<  "    </DataArray>\n" << endl;

    // writing scalar function c
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"WaterVapour\" format=\"ascii\">\n"  << endl;

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Atmosphere_vts_File << c.x[ i ][ j ][ k ]  << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }
    Atmosphere_vts_File <<  "\n"  << endl;
    Atmosphere_vts_File <<  "    </DataArray>\n" << endl;

    // writing u-component
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"u-Component\" format=\"ascii\">\n"  << endl;

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                Atmosphere_vts_File << aux_u.x[ i ][ j ][ k ]  << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }
    Atmosphere_vts_File <<  "\n"  << endl;
    Atmosphere_vts_File <<  "    </DataArray>\n" << endl;

    // writing v-component
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"v-Component\" format=\"ascii\">\n"  << endl;

    for ( int k = 0; k < km; k++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int i = 0; i < im; i++ )
            {
                Atmosphere_vts_File << aux_v.x[ i ][ j ][ k ]  << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }
    Atmosphere_vts_File <<  "\n"  << endl;
    Atmosphere_vts_File <<  "    </DataArray>\n" << endl;

    // writing w-component
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"w-Component\" format=\"ascii\">\n"  << endl;

    for ( int k = 0; k < km; k++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int i = 0; i < im; i++ )
            {
                Atmosphere_vts_File << aux_w.x[ i ][ j ][ k ]  << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }
    Atmosphere_vts_File <<  "\n"  << endl;
    Atmosphere_vts_File <<  "    </DataArray>\n"  << endl;
    Atmosphere_vts_File <<  "   </PointData>\n" << endl;
    Atmosphere_vts_File <<  "   <Points>\n"  << endl;
    Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;

    // transformation from spherical to cartesian coordinates
    for ( int k = 0; k < km; k++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int i = 0; i < im; i++ )
            {
                x = rad.z[ i ] * sin( the.z[ j ] ) * cos ( phi.z[ k ] );
                y = rad.z[ i ] * sin( the.z[ j ] ) * sin ( phi.z[ k ] );
                z = rad.z[ i ] * cos( the.z[ j ] );

                Atmosphere_vts_File << x << " " << y << " " << z  << endl;
            }
            Atmosphere_vts_File <<  "\n"  << endl;
        }
        Atmosphere_vts_File <<  "\n"  << endl;
    }

    Atmosphere_vts_File <<  "    </DataArray>\n"  << endl;
    Atmosphere_vts_File <<  "   </Points>\n"  << endl;
    Atmosphere_vts_File <<  "  </Piece>\n"  << endl;
    Atmosphere_vts_File <<  " </StructuredGrid>\n"  << endl;
    Atmosphere_vts_File <<  "</VTKFile>\n"  << endl;

    Atmosphere_vts_File.close();
}

void PostProcess_Atmosphere::paraview_panorama_vts (string &Name_Bathymetry_File, int &pressure_iter, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, double &co2_0, Array &h, Array &t, Array &p_dyn, Array &p_stat, Array &BuoyancyForce, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Q_Sensible, Array &IceLayer, Array &epsilon_3D, Array &P_rain, Array &P_snow ) {
    double x, y, z, dx, dy, dz;

    // max_u = max_v = max_w = max_t = max_c = max_co2 = max_cloud = max_ice = max_P_rain = max_P_snow = max_P_conv = max_M_u = max_M_d = max_p_dyn = max_p_stat = 0.;
    // max_Rain = max_Rain_super = max_Ice = max_Latency = max_Q_Sensible = max_Precipitation = 0.;
    // max_t_Evaporation = max_t_Condensation = max_t_evap_3D = max_t_cond_3D = 0.;
    // max_precipitable_water = max_IceAir = max_Q_bottom = max_Q_latent = max_Q_sensible = 0.;
    // max_t_Evaporation_Penman = max_t_Evaporation_Haude = max_Q_Radiation = max_buoyancy_force = 0.;
    // max_Q_t_Evaporation = max_precipitation_NASA = max_Water = max_Water_super = max_Vegetation = max_IceLayer = 0.;

    string Atmosphere_panorama_vts_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Atm_panorama_" + std::to_string(pressure_iter) + ".vts";
    ofstream Atmosphere_panorama_vts_File;
    Atmosphere_panorama_vts_File.precision(4);
    Atmosphere_panorama_vts_File.setf(ios::fixed);
    Atmosphere_panorama_vts_File.open(Atmosphere_panorama_vts_File_Name);

    if (!Atmosphere_panorama_vts_File.is_open()) {
        cerr << "ERROR: could not open panorama_vts file at " << __FILE__ << " line " << __LINE__ << "\n";
        abort();
    }

    // begin writing
    Atmosphere_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Atmosphere_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

    //        Atmosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography Temperature Pressure WaterVapour Latency Rain Ice RainSuper IceLayer\">\n"  << endl;
    Atmosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography Temperature CondensationTemp EvaporationTemp Epsilon_3D PressureDynamic PressureStatic WaterVapour CloudWater CloudIce CO2-Concentration Latency Rain RainSuper Ice PrecipitationRain PrecipitationSnow PrecipitationConv Updraft Downdraft\">\n"  << endl;

    // writing u, v und w velocity components in cartesian coordinates
    Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    /*
        for ( int k = 1; k < km-1; k++ )
        {
            for ( int j = 1; j < jm-1; j++ )
            {
                for ( int i = 1; i < im-1; i++ )
                {
                    if ( fabs ( u.x[ i ][ j ][ k ] ) > max_u ) 
                    {
                        max_u = fabs ( u.x[ i ][ j ][ k ] );
                        if ( max_u == 0. ) max_u = 1.e-6;
                    }
                    if ( fabs ( v.x[ i ][ j ][ k ] ) > max_v ) 
                    {
                        max_v = fabs ( v.x[ i ][ j ][ k ] );
                        if ( max_v == 0. ) max_v = 1.e-6;
                    }
                    if ( fabs ( w.x[ i ][ j ][ k ] ) > max_w ) 
                    {
                        max_w = fabs ( w.x[ i ][ j ][ k ] );
                        if ( max_w == 0. ) max_w = 1.e-6;
                    }
                }
            }
        }
    */

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                // transformtion from spherical to cartesian coordinates for representation in ParaView
                //                    Atmosphere_panorama_vts_File << u.x[ i ][ j ][ k ] / max_u << " " << v.x[ i ][ j ][ k ] / max_v << " " << w.x[ i ][ j ][ k ] / max_w << endl;
                Atmosphere_panorama_vts_File << u.x[ i ][ j ][ k ] << " " << v.x[ i ][ j ][ k ] << " " << w.x[ i ][ j ][ k ] << endl;
            }
            Atmosphere_panorama_vts_File <<  "\n"  << endl;
        }
        Atmosphere_panorama_vts_File <<  "\n"  << endl;
    }
    Atmosphere_panorama_vts_File <<  "\n"  << endl;
    Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

    dump_array("Topography", h, 1.0, Atmosphere_panorama_vts_File);

    // writing of temperature
    Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
                //                    Atmosphere_panorama_vts_File << t.x[ i ][ j ][ k ] / max_t << endl;
                Atmosphere_panorama_vts_File << t.x[ i ][ j ][ k ] * t_0 - t_0 << endl;
            }
            Atmosphere_panorama_vts_File <<  "\n"  << endl;
        }
        Atmosphere_panorama_vts_File <<  "\n"  << endl;
    }
    Atmosphere_panorama_vts_File <<  "\n"  << endl;
    Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

    // dump_array("CondensationTemp", t_cond_3D, 1.0, Atmosphere_panorama_vts_File);
    // dump_array("EvaporationTemp", t_evap_3D, 1.0, Atmosphere_panorama_vts_File);
    dump_array("Epsilon_3D", epsilon_3D, 1.0, Atmosphere_panorama_vts_File);
    dump_array("PressureDynamic", p_dyn, u_0 * u_0 * r_air, Atmosphere_panorama_vts_File);
    // dump_array("PressureStatic", p_stat, 1.0, Atmosphere_panorama_vts_File);
    dump_array("BuoyancyForce", BuoyancyForce, 1.0, Atmosphere_panorama_vts_File);
    dump_array("WaterVapour", c, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("CloudWater", cloud, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("CloudIce", ice, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationRain", P_rain, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("PrecipitationSnow", P_snow, 1000.0, Atmosphere_panorama_vts_File);
    // dump_array("Updraft", M_u, 1.0, Atmosphere_panorama_vts_File);
    // dump_array("Downdraft", M_d, 1.0, Atmosphere_panorama_vts_File);
    // dump_array("PrecipitationConv", P_co2_nv, 1000.0, Atmosphere_panorama_vts_File);
    dump_array("CO2-Concentration", co2, 1.0, Atmosphere_panorama_vts_File);
    dump_array("Latency", Latency, 1.0, Atmosphere_panorama_vts_File);
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

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            for ( int i = 0; i < im; i++ ) {
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

void PostProcess_Atmosphere::dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f) {
    f << "SCALARS " << desc << " float " << 1 << endl;
    f << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            f << (a.x[i][j][k] * multiplier) << endl;
        }
    }
}

void PostProcess_Atmosphere::dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f) {
    f << "SCALARS " << desc << " float " << 1 << endl;
    f << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            f << (a.y[j][k] * multiplier) << endl;
        }
    }
}

void PostProcess_Atmosphere::paraview_vtk_radial( string &Name_Bathymetry_File, int &i_radial, int &pressure_iter, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, double &co2_0, double &radiation_equator, Array &h, Array &p_dyn, Array &p_stat, Array &t_cond_3D, Array &t_evap_3D, Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Q_Sensible, Array &IceLayer, Array &epsilon_3D, Array &P_rain, Array &P_snow, Array_2D &Evaporation, Array_2D &Condensation, Array_2D &precipitable_water, Array_2D &Q_bottom, Array_2D &Radiation_Balance, Array_2D &Q_Radiation, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Evaporation_Penman, Array_2D &Evaporation_Haude, Array_2D &Q_Evaporation, Array_2D &precipitation_NASA, Array_2D &Vegetation, Array_2D &albedo, Array_2D &epsilon, Array_2D &Precipitation) {
    double x, y, z, dx, dy;

    // max_u = max_v = max_w = max_t = max_c = max_co2 = max_cloud = max_ice = max_P_rain = max_P_snow = max_P_conv = max_M_u = max_M_d = max_p_dyn = max_p_stat = 0.;
    // max_Rain = max_Rain_super = max_Ice = max_Latency = max_Q_Sensible = max_Precipitation = max_albedo = max_epsilon = 0.;
    // max_t_Evaporation = max_t_Condensation = max_t_evap_3D = max_t_cond_3D = 0.;
    // max_precipitable_water = max_IceAir = max_Q_bottom = max_Q_latent = max_Q_sensible = 0.;
    // max_t_Evaporation_Penman = max_t_Evaporation_Haude = max_Q_Radiation = max_Radiation_Balance = max_buoyancy_force = 0.;
    // max_Q_t_Evaporation = max_precipitation_NASA = max_Water = max_Water_super = max_Vegetation = max_IceLayer = max_radiation_3D = 0.;

    string Atmosphere_radial_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Atm_radial_" + std::to_string(i_radial) + "_" + std::to_string(pressure_iter) + ".vtk";
    ofstream Atmosphere_vtk_radial_File;
    Atmosphere_vtk_radial_File.precision (4);
    Atmosphere_vtk_radial_File.setf(ios::fixed);
    Atmosphere_vtk_radial_File.open(Atmosphere_radial_File_Name);

    if (!Atmosphere_vtk_radial_File.is_open()) {
        cerr << "ERROR: could not open paraview_vtk file '" << Atmosphere_radial_File_Name << "' at " << __FILE__ << " line " << __LINE__ << "\n";
        abort();
    }

    // begin writing
    Atmosphere_vtk_radial_File <<  "# vtk DataFile Version 3.0" << endl;
    Atmosphere_vtk_radial_File <<  "Radial_Data_Atmosphere_Circulation\n";
    Atmosphere_vtk_radial_File <<  "ASCII" << endl;
    Atmosphere_vtk_radial_File <<  "DATASET STRUCTURED_GRID" << endl;
    Atmosphere_vtk_radial_File <<  "DIMENSIONS " << km << " "<< jm << " " << 1 << endl;
    Atmosphere_vtk_radial_File <<  "POINTS " << jm * km << " float" << endl;

    // transformation from spherical to cartesian coordinates
    x = 0.;
    y = 0.;
    z = 0.;

    dx = .1;
    dy = .1;

    for ( int j = 0; j < jm; j++ ) {
        for ( int k = 0; k < km; k++ ) {
            if ( k == 0 ) y = 0.;
            else y = y + dy;

            Atmosphere_vtk_radial_File << x << " " << y << " "<< z << endl;
        }
        y = 0.;
        x = x + dx;
    }

    Atmosphere_vtk_radial_File <<  "POINT_DATA " << jm * km << endl;

    if ( i_radial != 0 ) {
        dump_radial("u-Component", u, 1., i_radial, Atmosphere_vtk_radial_File);
    }

    dump_radial("v-Component", v, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("w-Component", w, 1., i_radial, Atmosphere_vtk_radial_File);

    // writing temperature
    Atmosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for ( int j = 0; j < jm; j++ ) {
        for ( int k = 0; k < km; k++ ) {
//                Atmosphere_vtk_radial_File << t.x[ i_radial ][ j ][ k ] / max_t << endl;
            Atmosphere_vtk_radial_File << t.x[ i_radial ][ j ][ k ] * t_0 - t_0 << endl;
        }
    }

    dump_radial("CondensationTemp", t_cond_3D, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("EvaporationTemp", t_evap_3D, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Epsilon_3D", epsilon_3D, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("WaterVapour", c, 1000., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudWater", cloud, 1000., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CloudIce", ice, 1000., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("CO2-Concentration", co2, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PressureDynamic", p_dyn, u_0 * u_0 * r_air, i_radial, Atmosphere_vtk_radial_File);
    // dump_radial("PressureStatic", p_stat, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Topography", h, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Latency", Latency, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("Q_Sensible", Q_Sensible, 1., i_radial, Atmosphere_vtk_radial_File);
    // dump_radial("Evaporation", Evaporation, 1., i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Condensation", Condensation, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Precipitation_NASA", precipitation_NASA, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("albedo", albedo, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("epsilon", epsilon, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("PrecipitableWater", precipitable_water, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Radiation_Balance", Radiation_Balance, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_Radiation", Q_Radiation, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_bottom", Q_bottom, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_latent", Q_latent, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Q_sensible", Q_sensible, 1., Atmosphere_vtk_radial_File);
    dump_radial_2d("Evaporation_Penman", Evaporation_Penman, 1., Atmosphere_vtk_radial_File);
    // dump_radial_2d("t_EvaporationHaude", Evaporation_Haude, 1., Atmosphere_vtk_radial_File); // normally this would be normalized
    // dump_radial_2d("Q_Evaporation", Q_Evaporation, 1., Atmosphere_vtk_radial_File); // normally this would be normalized
    dump_radial_2d("Vegetation", Vegetation, 1., Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationRain", P_rain, 1000., i_radial, Atmosphere_vtk_radial_File);
    dump_radial("PrecipitationSnow", P_snow, 1000., i_radial, Atmosphere_vtk_radial_File);
    dump_radial_2d("Precipitation", Precipitation, 1., Atmosphere_vtk_radial_File);

    // writing zonal u-v cell structure
    Atmosphere_vtk_radial_File <<  "VECTORS v-w-Cell float " << endl;
    for ( int j = 0; j < jm; j++ ) {
        for ( int k = 0; k < km; k++ ) {
//                Atmosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] / max_v << " " << w.x[ i_radial ][ j ][ k ] / max_w << " " << z << endl;
            Atmosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] << " " << w.x[ i_radial ][ j ][ k ] << " " << z << endl;
        }
    }

    Atmosphere_vtk_radial_File.close();
}

void PostProcess_Atmosphere::dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f) {
    f <<  "SCALARS " << desc << " float " << 1 << endl;
    f <<  "LOOKUP_TABLE default" << endl;

    for ( int i = 0; i < im; i++ ) {
        for ( int j = 0; j < jm; j++ ) {
            // Atmosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] / max_u << endl;
            f << (a.x[ i ][ j ][ k ] * multiplier) << endl;
        }
    }
}

void PostProcess_Atmosphere::paraview_vtk_zonal ( string &Name_Bathymetry_File, int &k_zonal, int &pressure_iter, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, double &co2_0, double &radiation_equator, Array &h, Array &p_dyn, Array &p_stat, Array &t_cond_3D, Array &t_evap_3D, Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Q_Sensible, Array &radiation_3D, Array &epsilon_3D, Array &P_rain, Array &P_snow, Array &S_v, Array &S_c, Array &S_i, Array &S_r, Array &S_s ) {
    double x, y, z, dx, dy;

    // stringstream Atmosphere_zonal_File_Name;

    // max_u = max_v = max_w = max_t = max_c = max_co2 = max_cloud = max_ice = max_P_rain = max_P_snow = max_P_conv = max_M_u = max_M_d = max_p_dyn = max_p_stat = 0.;
    // max_Rain = max_Rain_super = max_Ice = max_Latency = max_Q_Sensible = max_Precipitation = 0.;
    // max_t_Evaporation = max_t_Condensation = max_t_evap_3D = max_t_cond_3D = 0.;
    // max_precipitable_water = max_IceAir = max_Q_bottom = max_Q_latent = max_Q_sensible = 0.;
    // max_t_Evaporation_Penman = max_t_Evaporation_Haude = max_Q_Radiation = max_buoyancy_force = 0.;
    // max_Q_t_Evaporation = max_precipitation_NASA = max_Water = max_Water_super = max_Vegetation = max_IceLayer = 0.;

    string Atmosphere_zonal_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Atm_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(pressure_iter) + ".vtk";
    ofstream Atmosphere_vtk_zonal_File;
    Atmosphere_vtk_zonal_File.precision ( 4 );
    Atmosphere_vtk_zonal_File.setf ( ios::fixed );
    Atmosphere_vtk_zonal_File.open ( Atmosphere_zonal_File_Name);

    if (!Atmosphere_vtk_zonal_File.is_open()) {
        cerr << "ERROR: could not open vtk_zonal file at " << __FILE__ << " line " << __LINE__ << "\n";
        abort();
    }

    Atmosphere_vtk_zonal_File <<  "# vtk DataFile Version 3.0" << endl;
    Atmosphere_vtk_zonal_File <<  "Zonal_Data_Atmosphere_Circulation\n";
    Atmosphere_vtk_zonal_File <<  "ASCII" << endl;
    Atmosphere_vtk_zonal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Atmosphere_vtk_zonal_File <<  "DIMENSIONS " << jm << " "<< im << " " << 1 << endl;
    Atmosphere_vtk_zonal_File <<  "POINTS " << im * jm << " float" << endl;

    // transformation from spherical to cartesian coordinates
    x = 0.;
    y = 0.;
    z = 0.;

    dx = .1;
    dy = .05;

    for ( int i = 0; i < im; i++ ) {
        for ( int j = 0; j < jm; j++ ) {
            if ( j == 0 ) y = 0.;
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

    // writing of temperature
    Atmosphere_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
    Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;
    for ( int i = 0; i < im; i++ ) {
        for ( int j = 0; j < jm; j++ ) {
//                Atmosphere_vtk_zonal_File << t.x[ i ][ j ][ k_zonal ] / max_t << endl;
            Atmosphere_vtk_zonal_File << t.x[ i ][ j ][ k_zonal ] * t_0 - t_0 << endl;
        }
    }

    dump_zonal("CondensationTemp", t_cond_3D, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("EvaporationTemp", t_evap_3D, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Epsilon_3D", epsilon_3D, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("WaterVapour", c, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudWater", cloud, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CloudIce", ice, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationRain", P_rain, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PrecipitationSnow", P_snow, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_WaterVapour", S_v, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_CloudWater", S_c, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_CloudIce", S_i, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_Rain", S_r, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Source_Snow", S_s, 1000., k_zonal, Atmosphere_vtk_zonal_File);
    // dump_zonal("Updraft", M_u, 1., k_zonal, Atmosphere_vtk_zonal_File);
    // dump_zonal("Downdraft", M_d, 1., k_zonal, Atmosphere_vtk_zonal_File);
    // dump_zonal("PrecipitationConv", P_co2_nv, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("CO2-Concentration", co2, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PressureDynamic", p_dyn, u_0 * u_0 * r_air, k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("PressureStatic", p_stat, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Topography", h, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Latency", Latency, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Q_Sensible", Q_Sensible, 1., k_zonal, Atmosphere_vtk_zonal_File);
    // dump_zonal("t_cond_3D", t_cond_3D, 1., k_zonal, Atmosphere_vtk_zonal_File);
    // dump_zonal("t_evap_3D", t_Evap_3D, 1., k_zonal, Atmosphere_vtk_zonal_File);
    dump_zonal("Radiation", radiation_3D, 1., k_zonal, Atmosphere_vtk_zonal_File);

    // writing zonal u-v cell structure
    Atmosphere_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for ( int i = 0; i < im; i++ ) {
        for ( int j = 0; j < jm; j++ ) {
            // Atmosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] / max_u << " " << v.x[ i ][ j ][ k_zonal ] / max_v << " " << z << endl;
            Atmosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] << " " << v.x[ i ][ j ][ k_zonal ] << " " << z << endl;
        }
    }

    Atmosphere_vtk_zonal_File.close();
}

void PostProcess_Atmosphere::dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f) {
    f << "SCALARS " << desc << " float " << 1 << endl;
    f << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < im; i++) {
        for (int k = 0; k < km; k++) {
            // Atmosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] / max_u << endl;
            f << (a.x[i][j][k] * multiplier) << endl;
        }
    }
}

void PostProcess_Atmosphere::paraview_vtk_longal (string &Name_Bathymetry_File, int &j_longal, int &pressure_iter, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, double &co2_0, double &radiation_equator, Array &h, Array &p_dyn, Array &p_stat, Array &t_cond_3D, Array &t_evap_3D, Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Q_Sensible, Array &IceLayer, Array &epsilon_3D, Array &P_rain, Array &P_snow )
{
    double x, y, z, dx, dz;

    string Atmosphere_longal_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Atm_longal_" + std::to_string(j_longal) + "_" + std::to_string(pressure_iter) + ".vtk";
    ofstream Atmosphere_vtk_longal_File;
    Atmosphere_vtk_longal_File.precision(4);
    Atmosphere_vtk_longal_File.setf(ios::fixed);
    Atmosphere_vtk_longal_File.open(Atmosphere_longal_File_Name);

    if (!Atmosphere_vtk_longal_File.is_open()) {
        cerr << "ERROR: could not open vtk_longal file at " << __FILE__ << " line " << __LINE__ << "\n";
        abort();
    }

    Atmosphere_vtk_longal_File <<  "# vtk DataFile Version 3.0" << endl;
    Atmosphere_vtk_longal_File <<  "Longitudinal_Data_Atmosphere_Circulation\n";
    Atmosphere_vtk_longal_File <<  "ASCII" << endl;
    Atmosphere_vtk_longal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Atmosphere_vtk_longal_File <<  "DIMENSIONS " << km << " "<< im << " " << 1 << endl;
    Atmosphere_vtk_longal_File <<  "POINTS " << im * km << " float" << endl;

    // transformation from spherical to cartesian coordinates
    x = 0.;
    y = 0.;
    z = 0.;

    dx = .1;
    dz = .025;

    for ( int i = 0; i < im; i++ ) {
        for ( int k = 0; k < km; k++ ) {
            if ( k == 0 ) {
                z = 0.;
            } else {
                z = z + dz;
            }

            Atmosphere_vtk_longal_File << x << " " << y << " "<< z << endl;
        }
        z = 0.;
        x = x + dx;
    }

    Atmosphere_vtk_longal_File <<  "POINT_DATA " << im * km << endl;

    dump_longal("u-Component", u, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("v-Component", v, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("w-Component", v, 1., j_longal, Atmosphere_vtk_longal_File);

    // writing of temperature
    Atmosphere_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
    Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;
    for ( int i = 0; i < im; i++ ) {
        for ( int k = 0; k < km; k++ ) {
//                Atmosphere_vtk_longal_File << t.x[ i ][ j_longal ][ k ] / max_t << endl;
            Atmosphere_vtk_longal_File << t.x[ i ][ j_longal ][ k ] * t_0 - t_0 << endl;
        }
    }

    dump_longal("CondensationTemp", t_cond_3D, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("EvaporationTemp", t_evap_3D, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Epsilon_3D", epsilon_3D, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("WaterVapour", c, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudWater", cloud, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CloudIce", ice, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationRain", P_rain, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PrecipitationSnow", P_snow, 1000., j_longal, Atmosphere_vtk_longal_File);
    // dump_longal("Updraft", M_u, 1., j_longal, Atmosphere_vtk_longal_File);
    // dump_longal("Downdraft", M_d, 1., j_longal, Atmosphere_vtk_longal_File);
    // dump_longal("PrecipitationConv", P_co2_nv, 1000., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("CO2-Concentration", co2, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("PressureDynamic", p_dyn, u_0 * u_0 * r_air, j_longal, Atmosphere_vtk_longal_File);
    // dump_longal("PressureStatic", p_stat, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Topography", h, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Latency", Latency, 1., j_longal, Atmosphere_vtk_longal_File);
    dump_longal("Q_Sensible", Q_Sensible, 1., j_longal, Atmosphere_vtk_longal_File);
    // dump_longal("t_cond_3D", t_cond_3D, 1., j_longal, Atmosphere_vtk_longal_File);
    // dump_longal("t_evap_3D", t_evap_3D, 1., j_longal, Atmosphere_vtk_longal_File);
    // dump_longal("IceLayer", IceLayer, 1., j_longal, Atmosphere_vtk_longal_File);

    // writing longitudinal u-v cell structure
    Atmosphere_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for ( int i = 0; i < im; i++ ) {
        for ( int k = 0; k < km; k++ ) {
            // Atmosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] / max_u << " " << y << " " << w.x[ i ][ j_longal ][ k ] / max_w << endl;
            Atmosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] << " " << y << " " << w.x[ i ][ j_longal ][ k ] << endl;
        }
    }

    Atmosphere_vtk_longal_File.close();
}

void PostProcess_Atmosphere::Atmosphere_PlotData ( string &Name_Bathymetry_File, double u_0, double t_0, Array &v, Array &w, Array &t, Array &c, Array_2D &Precipitation, Array_2D &precipitable_water ) {
    string Name_PlotData_File = output_path + "/[" + Name_Bathymetry_File + "]_PlotData_Atm.xyz";
    ofstream PlotData_File;
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);
    PlotData_File.open(Name_PlotData_File);

    if (!PlotData_File.is_open()) {
        cerr << "ERROR: could not open PlotData file at " << __FILE__ << " line " << __LINE__ << "\n";
        abort();
    }

    PlotData_File << " latitude ( ° )" << "  , " << "longitude ( ° )" << "  ,    " << "v-velocity ( m/s )" << "   ,   " << "w-velocity ( m/s )" << "   ,   " << "temperature ( °C )" << "   ,  " << "water_vapour ( g/kg )" << "   ,   " << "precipitation ( mm )" << "   ,   " <<  "precipitable water ( mm )" << endl;

    for ( int k = 0; k < km; k++ ) {
        for ( int j = 0; j < jm; j++ ) {
            PlotData_File << k << " " << j << " " << v.x[ 0 ][ j ][ k ] * u_0 << " " << w.x[ 0 ][ j ][ k ] * u_0 << " " << t.x[ 0 ][ j ][ k ] * t_0 - t_0 << " " << c.x[ 0 ][ j ][ k ] * 1000. << " " << Precipitation.y[ j ][ k ] << " " << precipitable_water.y[ j ][ k ] << " " <<  endl;
        }
    }

    PlotData_File.close();
}
