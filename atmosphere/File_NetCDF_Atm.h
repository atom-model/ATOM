/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write computational results in netCDF-file format
*/


#include <iostream>
#include <fstream>
#include <cstring>

#include "Array.h"
#include "Array_2D.h"

#include <netcdf.h>


#ifndef _FILE_NETCDF_
#define _FILE_NETCDF_

// Definitions to build netCDF-files

#define ERR(e) { printf ( "Error in %s line %d: %s\n" , __FILE__, __LINE__, nc_strerror ( e ) ); return 2; }

using namespace std;

class File_NetCDF {
	private:
		int ncid, lon_dimid, lat_dimid, lvl_dimid, rec_dimid, v_varid, w_varid, h_varid, prec_varid, precwat_varid;
		int lat_varid, lon_varid;
		int dimids[ 4 ];
		int rec;
		int retval;
		int im, jm, km;
		int NDIMS, NLAT, NLON, START_LAT, START_LON, NREC, NLVL;

		size_t start[ 4 ], count[ 4 ];					// store locations for netCDF-file

		double lats[ 182 ], lons[ 362 ];			// coordinate arrays for netCDF
		double v_air[ 3 ][ 182 ][ 362 ];			// auxilliar arrays for v-componente for the netCDF-file
		double w_air[ 3 ][ 182 ][ 362 ];			// auxilliar arrays for w-componente for the netCDF-file
		double h_air[ 3 ][ 182 ][ 362 ];			// auxilliar arrays for depth h for the netCDF-file

	public:
		File_NetCDF ( int, int, int );
		~File_NetCDF ();
		int out_NetCDF (const string &, const string &, Array &, Array &, Array &, Array_2D &, Array_2D &);
};
#endif

