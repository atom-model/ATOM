/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
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


#ifndef _FILE_NETCDF_
#define _FILE_NETCDF_

// Definitions to build netCDF-files

#define ERR(e) { printf ( "Error : %s\n" , nc_strerror ( e ) ); return 2; }

using namespace std;



class File_NetCDF
{
	private:
		int ncid, lon_dimid, lat_dimid, lvl_dimid, rec_dimid, v_varid, w_varid, h_varid, upwelling_varid, downwelling_varid, bottomwater_varid;
		int lat_varid, lon_varid;
		int dimids[ 4 ];
		int lvl, lat, lon, rec;
		int retval;
		int im, jm, km;
		int NDIMS, NLAT, NLON, START_LAT, START_LON, NREC, NLVL;

		size_t start[ 4 ], count[ 4 ];					// store locations for netCDF-file

		double lats[ 182 ], lons[ 362 ];			// coordinate arrays for netCDF
		double v_water[ 3 ][ 182 ][ 362 ];		// auxilliar arrays for v-velocity component for netCDF-file
		double w_water[ 3 ][ 182 ][ 362 ];		// auxilliar arrays for w-velocity component for netCDF-file
		double h_water[ 3 ][ 182 ][ 362 ];		// auxilliar arrays for depth h for netCDF-file

		static const char LAT_NAME [ ], LON_NAME [ ], LVL_NAME [ ], REC_NAME [ ], UNITS [ ], DEGREES_EAST [ ], DEGREES_NORTH [ ], V_NAME [ ], W_NAME [ ], H_NAME [ ], UPWELLING_NAME [ ], DOWNWELLING_NAME [ ], BOTTOMWATER_NAME [ ], v_units [ ], w_units [ ], h_units [ ], upwelling_units [ ], downwelling_units [ ], bottomwater_units [ ];

		Array_2D		up ( int, int, double ), down ( int, int, double ), bot ( int, int, double );
		Array				v_w ( int, int, int, double ), w_w ( int, int, int, double ), h_w ( int, int, int, double );


	public:
		File_NetCDF (  int, int, int );

		~File_NetCDF ();

		double out_NetCDF ( const string &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D & );

};
#endif

