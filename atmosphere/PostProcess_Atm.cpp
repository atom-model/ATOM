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
#include <iomanip>

#include "PostProcess_Atm.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

PostProcess_Atmosphere::PostProcess_Atmosphere( int im, int jm, int km, string &output_path )
{
    this->im = im;
    this->jm = jm;
    this->km = km;
    this->output_path = output_path;
}

void PostProcess_Atmosphere::Atmosphere_v_w_Transfer ( string &Name_Bathymetry_File, double u_0, Array &v, Array &w, Array &t, Array &p_dyn, Array_2D &Evaporation_Dalton, Array_2D &Precipitation )
{
    string Name_v_w_Transfer_File = output_path + "/[" + Name_Bathymetry_File + "]_Transfer_Atm.vw";
    ofstream v_w_Transfer_File;
    v_w_Transfer_File.precision(4);
    v_w_Transfer_File.setf(ios::fixed);
    v_w_Transfer_File.open(Name_v_w_Transfer_File);

    if (!v_w_Transfer_File.is_open()) {
        cout << "ERROR: transfer file name in atmosphere: " << Name_v_w_Transfer_File << "\n";
        cerr << "ERROR: could not open transfer file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    for ( int j = 0; j < jm; j++ ) {
        for ( int k = 0; k < km; k++ ) {
            v_w_Transfer_File << v.x[ 0 ][ j ][ k ] * u_0 << " " << w.x[ 0 ][ j ][ k ] * u_0 << " " << t.x[ 0 ][ j ][ k ] << " " << p_dyn.x[ 0 ][ j ][ k ] << " " << Evaporation_Dalton.y[ j ][ k ] << " " << Precipitation.y[ j ][ k ] << endl;    // dimensional v and w values
        }
    }
    v_w_Transfer_File.close();
}


void PostProcess_Atmosphere::Atmosphere_PlotData ( string &Name_Bathymetry_File, int iter_cnt, double u_0, double t_0, 
        Array &h, Array &v, Array &w, Array &t, Array &c, Array_2D &Precipitation, Array_2D &precipitable_water, 
        Array_2D &Evaporation_Dalton )
{
    string Name_PlotData_File = output_path + "/[" + Name_Bathymetry_File + "]_PlotData_Atm"+
        (iter_cnt > 0 ? "_"+to_string(iter_cnt) : "") + ".xyz";
    ofstream PlotData_File(Name_PlotData_File);
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);

    if (!PlotData_File.is_open())
    {
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    PlotData_File << "lons(deg)" << ", " << "lats(deg)" << ", " << "topography" << ", " << "v-velocity(m/s)" << ", " 
        << "w-velocity(m/s)" << ", " << "velocity-mag(m/s)" << ", " << "temperature(Celsius)" << ", " 
        << "water_vapour(g/kg)" << ", " << "precipitation(mm)" << ", " <<  "precipitable water(mm)" << ", " 
        << "Evaporation_Dalton (mm/day) " <<endl;

    for ( int k = 0; k < km; k++ ){
        for ( int j = 0; j < jm; j++ ){
            double vel_mag = sqrt ( pow ( v.x[ 0 ][ j ][ k ] * u_0, 2 ) + pow ( w.x[ 0 ][ j ][ k ] * u_0, 2 ) );
            PlotData_File << k << " " << 90-j << " " << h.x[ 0 ][ j ][ k ] << " " << v.x[ 0 ][ j ][ k ] * u_0 << " " 
                << w.x[ 0 ][ j ][ k ] * u_0 << " " << vel_mag << " " << t.x[ 0 ][ j ][ k ] * t_0 - t_0 << " " 
                << c.x[ 0 ][ j ][ k ] * 1000. << " " << Precipitation.y[ j ][ k ] << " " << precipitable_water.y[ j ][ k ] 
                << " " <<  Evaporation_Dalton.y[ j ][ k ] << " " << endl;
        }
    }
}


void PostProcess_Atmosphere::save(  const string &filename, const std::vector<string> &field_names, 
                                    const std::vector<Vector3D<>* > &data, unsigned layer )
{
    assert(field_names.size() == data.size());
    ofstream ofs(filename);

    if (!ofs.is_open())
    {
        cerr << "ERROR: unable to open file " << filename << "\n";
        abort();
    }

    for(std::size_t i = 0; i < field_names.size(); i++)
    {
        ofs << "lons(deg)" << " ";
        ofs << "lats(deg)" << " ";
        ofs << field_names[i] << " ";
    }
    ofs << std::endl;

    for ( int k = 0; k < km; k++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            ofs << k << " " << j << " ";
            for(std::size_t i = 0; i < data.size(); i++)
            {   
                ofs << (*data[i])(layer,j,k) << " ";
            }
            ofs << std::endl;
        }
    }
}

