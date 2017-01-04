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
#include <cstring>
#include <fstream>
#include <netcdf.h>

#include "File_NetCDF_Hyd.h"

using namespace std;


const char File_NetCDF_Hyd::LAT_NAME [ ] = "latitude";
const char File_NetCDF_Hyd::LON_NAME [ ] = "longitude";
const char File_NetCDF_Hyd::LVL_NAME [ ] = "level";
const char File_NetCDF_Hyd::REC_NAME [ ] = "time";
const char File_NetCDF_Hyd::UNITS [ ] = "units";
const char File_NetCDF_Hyd::DEGREES_EAST [ ] = "degrees_east";
const char File_NetCDF_Hyd::DEGREES_NORTH [ ] = "degrees_north";
const char File_NetCDF_Hyd::V_NAME [ ] = "v_velocity_component";
const char File_NetCDF_Hyd::W_NAME [ ] = "w_velocity_component";
const char File_NetCDF_Hyd::H_NAME [ ] = "level_over_NN";
const char File_NetCDF_Hyd::UPWELLING_NAME [ ] = "upwelling";
const char File_NetCDF_Hyd::DOWNWELLING_NAME [ ] = "downwelling";
const char File_NetCDF_Hyd::BOTTOMWATER_NAME [ ] = "bottomwater";
const char File_NetCDF_Hyd::v_units [ ] = "m/s";
const char File_NetCDF_Hyd::w_units [ ] = "m/s";
const char File_NetCDF_Hyd::h_units [ ] = "m";
const char File_NetCDF_Hyd::upwelling_units [ ] = "m/s";
const char File_NetCDF_Hyd::downwelling_units [ ] = "m/s";
const char File_NetCDF_Hyd::bottomwater_units [ ] = "m/s";




File_NetCDF_Hyd::File_NetCDF_Hyd ( int im, int jm, int km )
{
	this-> im = im;
	this-> jm = jm;
	this-> km = km;

	NDIMS = 4;
	NLAT = jm +1;
	NLON = km + 1;
	START_LAT = 90;
	START_LON = -180;
	NREC = 1;
	NLVL = 3;
}


File_NetCDF_Hyd::~File_NetCDF_Hyd (){}




double File_NetCDF_Hyd::out_NetCDF ( const string &F_N, Array &v_w, Array &w_w, Array &h_w, Array_2D &up, Array_2D &down, Array_2D &bot )
{
//	results written in netCDF format

	for ( lat = 0; lat < NLAT; lat++ )		lats[ lat ] = START_LAT - ( double ) lat;
	for ( lon = 0; lon < NLON; lon++ )	lons[ lon ] = START_LON + ( double ) lon;

	if ( ( retval = nc_create ( F_N.c_str(), NC_CLOBBER, &ncid ) ) ) ERR ( retval );

	printf ( "***** successful creation of a netCDF-file through class NetCDF -> netCDF-file: %s\n\n", F_N.c_str() );

	if ( ( retval = nc_def_dim ( ncid, LVL_NAME, NLVL, &lvl_dimid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_dim ( ncid, LAT_NAME, NLAT, &lat_dimid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_dim ( ncid, LON_NAME, NLON, &lon_dimid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_dim ( ncid, REC_NAME, NC_UNLIMITED, &rec_dimid ) ) ) ERR ( retval );

	if ( ( retval = nc_def_var (ncid, LAT_NAME, NC_DOUBLE, 1, &lat_dimid, &lat_varid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_var ( ncid, LON_NAME, NC_DOUBLE, 1, &lon_dimid, &lon_varid ) ) ) ERR ( retval );

	if ( ( retval = nc_put_att_text (ncid, lat_varid, UNITS, strlen ( DEGREES_NORTH ), DEGREES_NORTH ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, lon_varid, UNITS, strlen ( DEGREES_EAST ), DEGREES_EAST ) ) ) ERR ( retval );

	dimids[ 0 ] = rec_dimid;
	dimids[ 1 ] = lvl_dimid;
	dimids[ 2 ] = lat_dimid;
	dimids[ 3 ] = lon_dimid;

	if ( ( retval = nc_def_var ( ncid, V_NAME, NC_DOUBLE, NDIMS, dimids, &v_varid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_var ( ncid, W_NAME, NC_DOUBLE, NDIMS, dimids, &w_varid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_var ( ncid, H_NAME, NC_DOUBLE, NDIMS, dimids, &h_varid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_var ( ncid, UPWELLING_NAME, NC_DOUBLE, NDIMS, dimids, &upwelling_varid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_var ( ncid, DOWNWELLING_NAME, NC_DOUBLE, NDIMS, dimids, &downwelling_varid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_var ( ncid, BOTTOMWATER_NAME, NC_DOUBLE, NDIMS, dimids, &bottomwater_varid ) ) ) ERR ( retval );

	if ( ( retval = nc_put_att_text ( ncid, v_varid, UNITS, strlen ( v_units ), v_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, w_varid, UNITS, strlen ( w_units ), w_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, h_varid, UNITS, strlen ( h_units ), h_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, upwelling_varid, UNITS, strlen ( upwelling_units ), upwelling_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, downwelling_varid, UNITS, strlen ( downwelling_units ), downwelling_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, bottomwater_varid, UNITS, strlen ( bottomwater_units ), bottomwater_units ) ) ) ERR ( retval );

	if ( ( retval = nc_enddef ( ncid ) ) ) ERR ( retval );

	if ( ( retval = nc_put_var_double ( ncid, lat_varid, &lats[ 0 ] ) ) ) ERR ( retval );
	if ( ( retval = nc_put_var_double ( ncid, lon_varid, &lons[ 0 ] ) ) ) ERR ( retval );

	count[ 0 ] = 1;
	count[ 1 ] = NLVL;
	count[ 2 ] = NLAT;
	count[ 3 ] = NLON;

	start[ 1 ] = 0;
	start[ 2 ] = 0;
	start[ 3 ] = 0;


	for ( lvl = 0; lvl < NLVL; lvl++ )
	{
		for ( lat = 0; lat < NLAT-1; lat++ )
		{
			for ( lon = 0; lon < NLON-1; lon++ )
			{
				v_water[ lvl ][ lat ][ lon ] = v_w.x[ im-1-lvl ][ lat ][ lon ];
				w_water[ lvl ][ lat ][ lon ] = w_w.x[ im-1-lvl ][ lat ][ lon ];
				h_water[ lvl ][ lat ][ lon ] = h_w.x[ im-1-lvl ][ lat ][ lon ];
			}
		}
	}


	for ( rec = 0; rec < NREC; rec++ )
	{
		start[ 0 ] = rec;
		if ( ( retval = nc_put_vara_double ( ncid, v_varid, start, count, &v_water[ 0 ][ 0 ][ 0 ] ) ) ) ERR ( retval );
		if ( ( retval = nc_put_vara_double ( ncid, w_varid, start, count, &w_water[ 0 ][ 0 ][ 0 ] ) ) ) ERR ( retval );
		if ( ( retval = nc_put_vara_double ( ncid, h_varid, start, count, &h_water[ 0 ][ 0 ][ 0 ] ) ) ) ERR ( retval );
	}

	if ( ( retval = nc_put_var_double ( ncid, upwelling_varid, &up.y[ 0 ][ 0 ] ) ) ) ERR ( retval );
	if ( ( retval = nc_put_var_double ( ncid, downwelling_varid, &down.y[ 0 ][ 0 ] ) ) ) ERR ( retval );
	if ( ( retval = nc_put_var_double ( ncid, bottomwater_varid, &bot.y[ 0 ][ 0 ] ) ) ) ERR ( retval );

	if ( ( retval = nc_close ( ncid ) ) ) ERR ( retval );

	printf ( "***** successful writing of results through class NetCDF -> netCDF-file: %s\n\n", F_N.c_str() );

	return 0;
}
