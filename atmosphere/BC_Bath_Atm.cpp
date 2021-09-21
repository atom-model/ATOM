/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <algorithm> 
#include "BC_Bath_Atm.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

BC_Bathymetry_Atmosphere::BC_Bathymetry_Atmosphere(int NASATemperature,
    int im, int jm, int km, double co2_vegetation, 
    double co2_land, double co2_ocean):
    im(im),
    jm(jm),
    km(km),
    NASATemperature(NASATemperature),
    co2_vegetation(co2_vegetation),
    co2_ocean(co2_ocean),
    co2_land(co2_land)
{}

BC_Bathymetry_Atmosphere::~BC_Bathymetry_Atmosphere(){}

void BC_Bathymetry_Atmosphere::BC_MountainSurface(string &topo_filename,
    double L_atm, Array_1D &rad, Array_2D &Topography, Array &h){
    cout.precision(8);
    cout.setf(ios::fixed);
    double zeta = 3.715;
    // default adjustment, h must be 0 everywhere
    h.initArray(im, jm, km, 0.);
    // reading data from file Name_Bathymetry_File_Read
    ifstream ifile(topo_filename);
    if(! ifile.is_open()) {
        std::cerr << "ERROR: could not open Name_Bathymetry_File file: " <<  topo_filename << std::endl;
        abort();
    }
    double lon, lat, height;
    int j, k;
    for(j = 0; j < jm && !ifile.eof(); j++){
        for(k = 0; k < km && !ifile.eof(); k++){
            height = -999; // in case the height is NaN
            ifile >> lon >> lat >> height;
            if(height < 0.){
                height = Topography.y[j][k] = 0.;
            }else{
                Topography.y[j][k] = height;
                for(int i = 0; i < im; i++){
                    if(height <= (exp(zeta * (rad.z[i] - 1.)) - 1) 
                        * (L_atm /(double)(im-1))){
                        i_h = i;
                        break;
                    }
                }
                for(int i = 0; i <= i_h; i++){
                    h.x[i][j][k] = 1.;
                }
            }
            if(ifile.fail()){
                ifile.clear();
                std::string tmp;
                std::getline(ifile, tmp);
                logger() << "bad data in topography at: " << lon << " " 
                    << lat << " " << tmp << std::endl;
            }
            //logger() << lon << " " << lat << " " << h.x[0][j][k] << std::endl;            
        }
    }
    if(j != jm || k != km){
        std::cerr << "wrong topography file size! aborting..."<<std::endl;
        abort();
    }
    // rewriting bathymetrical data from -180° _ 0° _ +180° coordinate system to 0°- 360°
    for(int j = 0; j < jm; j++){
        move_data(Topography.y[j], km);
        for(int i = 0; i < im; i++){
            move_data(h.x[i][j], km);
        }
    }
//  reduction and smoothing of peaks and needles in the bathymetry
    double h_center = 0.;
    for(int k = 2; k < km-2; k++){
        for(int j = 2; j < jm-2; j++){
            for(int i = 0; i <= im-2; i++){
                if((is_land(h, i, j, k)) && ((is_water(h, i, j+1, k)) 
                    && (is_water(h, i, j-1, k))) && ((is_water(h, i, j, k+1)) 
                    && (is_water(h, i, j, k-1)))){
                    h_center = h.x[i][j][k];
                    h.x[i][j][k] = 0.;
                 }
/*
                if((h_center == 1.) && ((is_water(h, i, j-2, k)) 
                    && (is_water(h, i, j, k+2)))){
                    h.x[i][j-1][k+1] = 0.;
                }
                if((h_center == 1.) && ((is_water(h, i, j-2, k)) 
                    && (is_water(h, i, j, k-2)))){ 
                    h.x[i][j-1][k-1] = 0.;
                }
                if((h_center == 1.) && ((is_water(h, i, j+2, k)) 
                    && (is_water(h, i, j, k-2)))){
                    h.x[i][j+1][k-1] = 0.;
                }
                if((h_center == 1.) && ((is_water(h, i, j+2, k)) 
                    && (is_water(h, i, j, k+2)))){ 
                    h.x[i][j+1][k+1] = 0.;
                }
*/
/*
                if((h_center == 1.) && ((is_land(h, i, j, k)) 
                    && (is_water(h, i, j, k+2)))){
                    h.x[i][j][k] = h.x[i][j][k+1] = 0.;
                }
                if((h_center == 1.) && ((is_land(h, i, j, k)) 
                    && (is_water(h, i, j+2, k)))){ 
                    h.x[i][j][k] = h.x[i][j+1][k] = 0.;
                }
*/
                if((i >= 1) && (h_center == 1.) && (is_water(h, 2, j, k))){
                    h.x[0][j][k] = 0.;
                }
/*
                if((h_center == 1.) && ((is_water(h, i, j+2, k)) 
                    && (is_water(h, i, j, k+2)))){ 
                    h.x[i][j+1][k+1] = 0.;
                }
*/
            }
        }
    }
    return;
}
/*
*
*/
void BC_Bathymetry_Atmosphere::land_oceanFraction(Array &h){
    // calculation of the ratio ocean to land, also addition and subtraction of CO2 of land, ocean and vegetation
    h_point_max = (jm-1) * (km-1);
    h_land = 0;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_land(h, 0, j, k))  h_land = h_land + h.x[0][j][k];
        }
    }
    h_ocean = h_point_max - h_land;
    ozean_land = (double) h_ocean/(double) h_land;
    cout.precision(3);
    cout << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      total number of points at constant height " << " = " 
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ') 
        << h_point_max << endl << setiosflags(ios::left) << setw(50) 
        << setfill('.') << "      number of points on the ocean surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << h_ocean << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      number of points on the land surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << h_land << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      ocean/land ratio " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << ozean_land 
        << endl << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      addition of CO2 by ocean surface " << " = " 
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ') 
        << co2_ocean << endl << setiosflags(ios::left) << setw(50) 
        << setfill('.') << "      addition of CO2 by land surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << co2_land << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      subtraction of CO2 by vegetation " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << co2_vegetation << endl << setiosflags(ios::left) 
        << setw(50) << "      valid for one single point on the surface"<< endl << endl;
    cout << endl;
    return;
}

