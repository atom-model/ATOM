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
#include <netcdf.h>

#include "File_NetCDF.h"

using namespace std;


File_NetCDF::File_NetCDF ( int im, int jm, int km )
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


File_NetCDF::~File_NetCDF () {}



double File_NetCDF::out_NetCDF (const string &Name_netCDF_File, Array &v, Array &w, Array &h, Array_2D &Precipitation, Array_2D &precipitable_water )
{
	const char LAT_NAME [ ] = "latitude";
	const char LON_NAME [ ] = "longitude";
	const char LVL_NAME [ ] = "level";
	const char REC_NAME [ ] = "time";
	const char UNITS [ ] = "units";
	const char DEGREES_EAST [ ] = "degrees_east";
	const char DEGREES_NORTH [ ] = "degrees_north";
	const char V_NAME [ ] = "v_velocity_component";
	const char W_NAME [ ] = "w_velocity_component";
	const char H_NAME [ ] = "level_over_NN";
	const char PREC_NAME [ ] = "precipitation";
	const char PRECWAT_NAME [ ] = "precipitable_water";
	const char v_units [ ] = "m/s";
	const char w_units [ ] = "m/s";
	const char h_units [ ] = "m";
	const char prec_units [ ] = "mm/d";
	const char precwat_units [ ] = "mm/d";


//	results written in netCDF format

	for ( int lat = 0; lat < NLAT; lat++ )		lats[ lat ] = START_LAT - ( double ) lat;
	for ( int lon = 0; lon < NLON; lon++ )	lons[ lon ] = START_LON + ( double ) lon;

	if ( ( retval = nc_create ( Name_netCDF_File.c_str(), NC_CLOBBER, &ncid ) ) ) ERR ( retval );

	printf ( "***** successful creation of a netCDF-file through class NetCDF -> netCDF-file: %s\n\n", Name_netCDF_File.c_str() );

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
	if ( ( retval = nc_def_var ( ncid, PREC_NAME, NC_DOUBLE, NDIMS, dimids, &prec_varid ) ) ) ERR ( retval );
	if ( ( retval = nc_def_var ( ncid, PRECWAT_NAME, NC_DOUBLE, NDIMS, dimids, &precwat_varid ) ) ) ERR ( retval );

	if ( ( retval = nc_put_att_text ( ncid, v_varid, UNITS, strlen ( v_units ), v_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, w_varid, UNITS, strlen ( w_units ), w_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, h_varid, UNITS, strlen ( h_units ), h_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, prec_varid, UNITS, strlen ( prec_units ), prec_units ) ) ) ERR ( retval );
	if ( ( retval = nc_put_att_text ( ncid, precwat_varid, UNITS, strlen ( precwat_units ), precwat_units ) ) ) ERR ( retval );

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

	for ( int lvl = 0; lvl < NLVL; lvl++ )
	{
		for ( int lat = 0; lat < NLAT-1; lat++ )
		{
			for ( int lon = 0; lon < NLON-1; lon++ )
			{
				v_air[ lvl ][ lat ][ lon ] = v.x[ im-1-lvl ][ lat ][ lon ];
				w_air[ lvl ][ lat ][ lon ] = w.x[ im-1-lvl ][ lat ][ lon ];
				h_air[ lvl ][ lat ][ lon ] = h.x[ im-1-lvl ][ lat ][ lon ];
			}
		}
	}

	for ( rec = 0; rec < NREC; rec++ )
	{
		start[ 0 ] = rec;
		if ( ( retval = nc_put_vara_double ( ncid, v_varid, start, count, &v_air[ 0 ][ 0 ][ 0 ] ) ) ) ERR ( retval );
		if ( ( retval = nc_put_vara_double ( ncid, w_varid, start, count, &w_air[ 0 ][ 0 ][ 0 ] ) ) ) ERR ( retval );
		if ( ( retval = nc_put_vara_double ( ncid, h_varid, start, count, &h_air[ 0 ][ 0 ][ 0 ] ) ) ) ERR ( retval );
	}

	if ( ( retval = nc_put_var_double ( ncid, prec_varid, &Precipitation.y[ 0 ][ 0 ] ) ) ) ERR ( retval );
	if ( ( retval = nc_put_var_double ( ncid, precwat_varid, &precipitable_water.y[ 0 ][ 0 ] ) ) ) ERR ( retval );

	if ( ( retval = nc_close ( ncid ) ) ) ERR ( retval );

	printf ( "***** successful writing of results through class NetCDF -> netCDF-file: %s\n\n", Name_netCDF_File.c_str() );

	return 0;
}
