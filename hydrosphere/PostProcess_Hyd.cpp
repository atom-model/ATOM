/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
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

#include "PostProcess_Hyd.h"

using namespace std;



PostProcess_Hydrosphere::PostProcess_Hydrosphere(int im, int jm, int km, const string &input_path, const string &output_path) {
	this->im = im;
	this->jm = jm;
	this->km = km;

	this->input_path = output_path;
	this->output_path = output_path;

}

PostProcess_Hydrosphere::~PostProcess_Hydrosphere() {}







void PostProcess_Hydrosphere::dump_array(const string &name, Array &a, double multiplier, ofstream &f) {
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




void PostProcess_Hydrosphere::dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f) {
    f << "SCALARS " << desc << " float " << 1 << endl;
    f << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            f << (a.x[i][j][k] * multiplier) << endl;
        }
    }
}



void PostProcess_Hydrosphere::dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f) {
    f << "SCALARS " << desc << " float " << 1 << endl;
    f << "LOOKUP_TABLE default" << endl;

    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            f << (a.y[j][k] * multiplier) << endl;
        }
    }
}




void PostProcess_Hydrosphere::dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f) {
    f <<  "SCALARS " << desc << " float " << 1 << endl;
    f <<  "LOOKUP_TABLE default" << endl;

    for ( int i = 0; i < im; i++ ) {
        for ( int j = 0; j < jm; j++ ) {
            f << (a.x[ i ][ j ][ k ] * multiplier) << endl;
        }
    }
}




void PostProcess_Hydrosphere::dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f) {
    f << "SCALARS " << desc << " float " << 1 << endl;
    f << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < im; i++) {
        for (int k = 0; k < km; k++) {
            f << (a.x[i][j][k] * multiplier) << endl;
        }
    }
}








void PostProcess_Hydrosphere::paraview_vts ( const string &Name_Bathymetry_File, int &n, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &p, Array &u, Array &v, Array &w, Array &c, Array &fup, Array &fvp, Array &fwp, Array &fcp, Array &fpp, Array &ftp, Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Buoyancy_Force, Array &Salt_Balance )
{
	double x, y, z, sinthe, sinphi, costhe, cosphi;

	// file administration
    string Hydrosphere_vts_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Hyd" + std::to_string(n) + ".vts";
//	string Hydrosphere_vts_File_Name << "/[" << Name_Bathymetry_File << "]_Hyd_Kreide_" << n << ".vts";
	ofstream Hydrosphere_vts_File;
	Hydrosphere_vts_File.precision ( 4 );
	Hydrosphere_vts_File.setf ( ios::fixed );
	string path = output_path + Name_Bathymetry_File;
	Hydrosphere_vts_File.open (path);

	if (!Hydrosphere_vts_File.is_open()) {
        cerr << "ERROR: could not open paraview_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
		abort();
	}

	Hydrosphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
	Hydrosphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
	Hydrosphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
	Hydrosphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

	//		Hydrosphere_vts_File <<  "   <PointData Vectors=\"Velocity Rotation\" Scalars=\"Topography Temperature Pressure Salinity u-Component v-Component w-Component\">\n"  << endl;
	Hydrosphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography Temperature Pressure Salinity\">\n"  << endl;

	// writing u, v und w velocity components in cartesian coordinates

	Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;


	for ( int k = 0; k < km; k++ ) {
		sinphi = sin( phi.z[ k ] );
		cosphi = cos( phi.z[ k ] );

		for ( int j = 0; j < jm; j++ ) {
			sinthe = sin( the.z[ j ] );
			costhe = cos( the.z[ j ] );

			for ( int i = 0; i < im; i++ ) {
// transformation from spherical to cartesian coordinates for presentation in ParaView
/*
				if ( h.x[ i ][ j ][ k ] == 1. )
				{
					u.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = 0.;
				}
*/
				fup.x[ i ][ j ][ k ] = sinthe * cosphi * u.x[ i ][ j ][ k ] + costhe * cosphi * v.x[ i ][ j ][ k ] - sinphi * w.x[ i ][ j ][ k ];
				fvp.x[ i ][ j ][ k ] = sinthe * sinphi * u.x[ i ][ j ][ k ] + sinphi * costhe * v.x[ i ][ j ][ k ] + cosphi * w.x[ i ][ j ][ k ];
				fwp.x[ i ][ j ][ k ] = costhe * u.x[ i ][ j ][ k ] - sinthe * v.x[ i ][ j ][ k ];

				Hydrosphere_vts_File << fup.x[ i ][ j ][ k ] << " " << fvp.x[ i ][ j ][ k ] << " " << fwp.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_vts_File <<  "\n"  << endl;
	Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;


/*
// writing aux_u, aux_v and aux_w components of velocity in cartesian coordinates

		Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Rotation\" format=\"ascii\">\n"  << endl;

		for ( int k = 0; k < km; k++ )
		{
			sinphi = sin( phi.z[ k ] );
			cosphi = cos( phi.z[ k ] );

			for ( int j = 0; j < jm; j++ )
			{
				sinthe = sin( the.z[ j ] );
				costhe = cos( the.z[ j ] );

				for ( int i = 0; i < im; i++ )
				{
// transformation from spherical to cartesian coordinates for presentation in ParaView

					rotu = sinthe * cosphi * aux_u.x[ i ][ j ][ k ] + costhe * cosphi * aux_v.x[ i ][ j ][ k ] - sinphi * aux_w.x[ i ][ j ][ k ];
					rotv = sinthe * sinphi * aux_u.x[ i ][ j ][ k ] + sinphi * costhe * aux_v.x[ i ][ j ][ k ] + cosphi * aux_w.x[ i ][ j ][ k ];
					rotw = costhe * aux_u.x[ i ][ j ][ k ] - sinthe * aux_v.x[ i ][ j ][ k ];

					Hydrosphere_vts_File << rotu << " " << rotv << " " << rotw  << endl;
				}
				Hydrosphere_vts_File <<  "\n"  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
		Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;
*/


	// writing of sea ground
	Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Topography\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_vts_File << h.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_vts_File <<  "\n"  << endl;
	Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;

	// writing of temperature
	Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_vts_File << t.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_vts_File <<  "\n"  << endl;
	Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;

	// writing temperature
	Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_vts_File << p.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_vts_File <<  "\n"  << endl;
	Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;

	// writing scalar function c
	Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Salinity\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_vts_File << c.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_vts_File <<  "\n"  << endl;
	Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;

/*
// writing u-component

		Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"u-Component\" format=\"ascii\">\n"  << endl;

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] == 1. ) u.x[ i ][ j ][ k ] = 0.;
//					Hydrosphere_vts_File << fup.x[ i ][ j ][ k ]  << endl;
					Hydrosphere_vts_File << u.x[ i ][ j ][ k ]  << endl;
				}
				Hydrosphere_vts_File <<  "\n"  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
		Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;



// writing v-component

		Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"v-Component\" format=\"ascii\">\n"  << endl;

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] == 1. ) v.x[ i ][ j ][ k ] = 0.;
//					Hydrosphere_vts_File << fvp.x[ i ][ j ][ k ]  << endl;
					Hydrosphere_vts_File << v.x[ i ][ j ][ k ]  << endl;
				}
				Hydrosphere_vts_File <<  "\n"  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
		Hydrosphere_vts_File <<  "    </DataArray>\n" << endl;



// writing w-component

		Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"w-Component\" format=\"ascii\">\n"  << endl;

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					if ( h.x[ i ][ j ][ k ] == 1. ) w.x[ i ][ j ][ k ] = 0.;
//					Hydrosphere_vts_File << fwp.x[ i ][ j ][ k ]  << endl;
					Hydrosphere_vts_File << w.x[ i ][ j ][ k ]  << endl;
				}
				Hydrosphere_vts_File <<  "\n"  << endl;
			}
			Hydrosphere_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_vts_File <<  "\n"  << endl;
		Hydrosphere_vts_File <<  "    </DataArray>\n"  << endl;
*/

	Hydrosphere_vts_File <<  "   </PointData>\n" << endl;
	Hydrosphere_vts_File <<  "   <Points>\n"  << endl;
	Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;

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

	// end writing

	Hydrosphere_vts_File.close();
}






void PostProcess_Hydrosphere::paraview_panorama_vts ( const string &Name_Bathymetry_File, int &pressure_iter, double &u_0, double &r_0_water, Array &h, Array &t, Array &p_dyn, Array &p_stat, Array &u, Array &v, Array &w, Array &c, Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion, Array &Buoyancy_Force, Array &Salt_Balance )
{
	double x, y, z, dx, dy, dz;


    string Atmosphere_panorama_vts_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Hyd_panorama_" + std::to_string(pressure_iter) + ".vts";
	ofstream Hydrosphere_panorama_vts_File;
	Hydrosphere_panorama_vts_File.precision ( 4 );
	Hydrosphere_panorama_vts_File.setf ( ios::fixed );
    Hydrosphere_panorama_vts_File.open(Atmosphere_panorama_vts_File_Name);

	if (!Hydrosphere_panorama_vts_File.is_open()) {
        cerr << "ERROR: could not open panorama_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
		abort();
	}

	// begin writing
	Hydrosphere_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
	Hydrosphere_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

	//		Hydrosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography Temperature Pressure Salinity SaltFinger BuoyancyForce SaltBalance\">\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography Temperature PressureDynamic PressureStatic Salinity\">\n"  << endl;

	// writing u, v und w velocity components in cartesian coordinates
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ ) {
		for ( int j = 0; j < jm; j++ ) {
			for ( int i = 0; i < im; i++ ) {
				// transformtion from spherical to cartesian coordinates for representation in ParaView
				Hydrosphere_panorama_vts_File << u.x[ i ][ j ][ k ] << " " << v.x[ i ][ j ][ k ] << " " << w.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


    dump_array("Topography", h, 1.0, Hydrosphere_panorama_vts_File);

    dump_array("u-velocity", u, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("v-velocity", v, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("w-velocity", w, 1.0, Hydrosphere_panorama_vts_File);

    dump_array("Temperature", t, 1.0, Hydrosphere_panorama_vts_File);

    dump_array("PressureDynamic", p_dyn, u_0 * u_0 * r_0_water, Hydrosphere_panorama_vts_File);
//    dump_array("PressureStatic", p_stat, 1.0, Hydrosphere_panorama_vts_File);

    dump_array("Salinity", c, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("Salt_Finger", Salt_Finger, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("SaltDiffusion", Salt_Diffusion, 1.0, Hydrosphere_panorama_vts_File);
    dump_array("SaltBalance", Salt_Balance, 1.0, Hydrosphere_panorama_vts_File);

    dump_array("BuoyancyForce", Buoyancy_Force, 1.0, Hydrosphere_panorama_vts_File);



/*
	// writing u-velocity
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"u-velocity\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ ) {
		for ( int j = 0; j < jm; j++ ) {
			for ( int i = 0; i < im; i++ ) {
				Hydrosphere_panorama_vts_File << u.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

	// writing v-velocity
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"v-velocity\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << v.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

	// writing w-velocity
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"w-velocity\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << w.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

	// writing of sea ground
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Topography\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << h.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

	// writing of temperature
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << t.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

	// writing dynamic pressure
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureDynamic\" format=\"ascii\">\n"  << endl;

	for ( int k = 1; k < km-1; k++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int i = 1; i < im-1; i++ )
			{
				if ( fabs ( p_dyn.x[ i ][ j ][ k ] ) > max_p_dyn )
				{
					max_p_dyn = fabs ( p_dyn.x[ i ][ j ][ k ] );
					if ( max_p_dyn == 0. ) max_p_dyn = 1.e-6;
				}
			}
		}
	}

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
//					Hydrosphere_panorama_vts_File << p_dyn.x[ i ][ j ][ k ] / max_p_dyn << endl;
				Hydrosphere_panorama_vts_File << p_dyn.x[ i ][ j ][ k ] << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing static pressure
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureStatic\" format=\"ascii\">\n"  << endl;

	for ( int k = 1; k < km-1; k++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int i = 1; i < im-1; i++ )
			{
				if ( fabs ( p_stat.x[ i ][ j ][ k ] ) > max_p_stat )
				{
					max_p_stat = fabs ( p_stat.x[ i ][ j ][ k ] );
					if ( max_p_stat == 0. ) max_p_stat = 1.e-6;
				}
			}
		}
	}

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
//					Hydrosphere_panorama_vts_File << p_stat.x[ i ][ j ][ k ] / max_p_stat << endl;
				Hydrosphere_panorama_vts_File << p_stat.x[ i ][ j ][ k ] << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


	// writing scalar function c
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Salinity\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << c.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

	// writing salt fingers
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SaltFinger\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << Salt_Finger.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


	// writing salt diffusion
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SaltDiffusion\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << Salt_Diffusion.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


	// writing salt balance
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SaltBalance\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << Salt_Balance.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


	// writing buoyancy force
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"BuoyancyForce\" format=\"ascii\">\n"  << endl;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Hydrosphere_panorama_vts_File << Buoyancy_Force.x[ i ][ j ][ k ]  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	}
	Hydrosphere_panorama_vts_File <<  "\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;
*/


	Hydrosphere_panorama_vts_File <<  "   </PointData>\n" << endl;
	Hydrosphere_panorama_vts_File <<  "   <Points>\n"  << endl;
	Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;


	// writing cartesian coordinates
	x = 0.;
	y = 0.;
	z = 0.;

	dx = .025;
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

	// end writing

	Hydrosphere_panorama_vts_File.close();
}







void PostProcess_Hydrosphere::paraview_vtk_longal ( const string &Name_Bathymetry_File, int &j_longal, int &pressure_iter, double &u_0, double &r_0_water, Array &h, Array &p_dyn, Array &p_stat, Array &t, Array &u, Array &v, Array &w, Array &c, Array &aux_u, Array &aux_v, Array &Salt_Finger, Array &Salt_Diffusion, Array &Buoyancy_Force, Array &Salt_Balance )
{
	double x, y, z, dx, dz;

	i_max = im;
	k_max = km;

// file administration

    string Hydrosphere_longal_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Hyd_longal_" + std::to_string(j_longal) + "_" + std::to_string(pressure_iter) + ".vtk";
	ofstream Hydrosphere_vtk_longal_File;
	Hydrosphere_vtk_longal_File.precision ( 4 );
	Hydrosphere_vtk_longal_File.setf ( ios::fixed );
    Hydrosphere_vtk_longal_File.open(Hydrosphere_longal_File_Name);


	if (!Hydrosphere_vtk_longal_File.is_open()) {
        cerr << "ERROR: could not open vtk_longal file " << __FILE__ << " at line " << __LINE__ << "\n";
		abort();
	}

	// begin writing
	Hydrosphere_vtk_longal_File <<  "# vtk DataFile Version 3.0" << endl;
	Hydrosphere_vtk_longal_File <<  "Longitudinal_Data_Hydrosphere_Circulation\n";
	Hydrosphere_vtk_longal_File <<  "ASCII" << endl;
	Hydrosphere_vtk_longal_File <<  "DATASET STRUCTURED_GRID" << endl;
	Hydrosphere_vtk_longal_File <<  "DIMENSIONS " << k_max << " "<< i_max << " " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "POINTS " << i_max * k_max << " float" << endl;

	// transformation from spherical to cartesian coordinates
	x = 0.;
	y = 0.;
	z = 0.;

	dx = .1;
	dz = .025;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( k == 0 ) z = 0.;
			else z = z + dz;

			Hydrosphere_vtk_longal_File << x << " " << y << " "<< z << endl;
		}
		z = 0.;
		x = x + dx;
	}



	Hydrosphere_vtk_longal_File <<  "POINT_DATA " << i_max * k_max << endl;




    dump_longal("Topography", h, 1., j_longal, Hydrosphere_vtk_longal_File);

    dump_longal("u-Component", u, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("v-Component", v, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("w-Component", w, 1., j_longal, Hydrosphere_vtk_longal_File);

    dump_longal("Temperature", t, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("Salinity", c, 1., j_longal, Hydrosphere_vtk_longal_File);

    dump_longal("PressureDynamic", p_dyn, u_0 * u_0 * r_0_water, j_longal, Hydrosphere_vtk_longal_File);
    // dump_longal("PressureStatic", p_stat, 1., j_longal, Hydrosphere_vtk_longal_File);

    dump_longal("SaltFinger", Salt_Finger, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("SaltDiffusion", Salt_Diffusion, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("SaltBalance", Salt_Balance, 1., j_longal, Hydrosphere_vtk_longal_File);
    dump_longal("BuoyancyForce", Buoyancy_Force, 1., j_longal, Hydrosphere_vtk_longal_File);







/*
	// writing u-component
	Hydrosphere_vtk_longal_File <<  "SCALARS u-Component float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing v-component
	Hydrosphere_vtk_longal_File <<  "SCALARS v-Component float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )

		{
			Hydrosphere_vtk_longal_File << v.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing w-component
	Hydrosphere_vtk_longal_File <<  "SCALARS w-Component float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" <<endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )

		{
			Hydrosphere_vtk_longal_File << w.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing of temperature
	Hydrosphere_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )

		{
			Hydrosphere_vtk_longal_File << t.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing salinity
	Hydrosphere_vtk_longal_File <<  "SCALARS Salinity float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )

		{
			Hydrosphere_vtk_longal_File << c.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing pressure
	Hydrosphere_vtk_longal_File <<  "SCALARS Pressure_dynamic float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ ) {
		for ( int k = 1; k < km-1; k++ ) {
			if ( fabs ( p_dyn.x[ i ][ j_longal ][ k ] ) > max_p_dyn )
			{
				max_p_dyn = fabs ( p_dyn.x[ i ][ j_longal ][ k ] );
				if ( max_p_dyn == 0. ) max_p_dyn = 1.e-6;
			}
		}
	}

	for ( int i = 0; i < im; i++ ) {
		for ( int k = 0; k < km; k++ ) {
//				Hydrosphere_vtk_longal_File << p_dyn.x[ i ][ j_longal ][ k ] / max_p_dyn << endl;
			Hydrosphere_vtk_longal_File << p_dyn.x[ i ][ j_longal ][ k ] << endl;
		}
	}


// writing static pressure
		Hydrosphere_vtk_longal_File <<  "SCALARS Pressure_static float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( p_stat.x[ i ][ j_longal ][ k ] ) > max_p_stat )
				{
					max_p_stat = fabs ( p_stat.x[ i ][ j_longal ][ k ] );
					if ( max_p_stat == 0. ) max_p_stat = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Hydrosphere_vtk_longal_File << p_stat.x[ i ][ j_longal ][ k ] / max_p_stat << endl;
				Hydrosphere_vtk_longal_File << p_stat.x[ i ][ j_longal ][ k ] << endl;
			}
		}


	// writing salt finger
	Hydrosphere_vtk_longal_File <<  "SCALARS SaltFinger float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_longal_File << Salt_Finger.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing salt diffusion
	Hydrosphere_vtk_longal_File <<  "SCALARS SaltDiffusion float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_longal_File << Salt_Diffusion.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing buoyancy force
	Hydrosphere_vtk_longal_File <<  "SCALARS BuoyancyForce float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_longal_File << Buoyancy_Force.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing salt balance
	Hydrosphere_vtk_longal_File <<  "SCALARS SaltBalance float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_longal_File << Salt_Balance.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// writing sea ground
	Hydrosphere_vtk_longal_File <<  "SCALARS Topography float " << 1 << endl;
	Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_longal_File << h.x[ i ][ j_longal ][ k ] << endl;
		}
	}
*/


	// writing longitudinal u-w-cell structure
	Hydrosphere_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] << " " << y << " " << w.x[ i ][ j_longal ][ k ] << endl;
		}
	}

	// end writing

	Hydrosphere_vtk_longal_File.close();
}







void PostProcess_Hydrosphere::paraview_vtk_radial ( const string &Name_Bathymetry_File, int &i_radial, int &pressure_iter, double &u_0, double &r_0_water, Array &h, Array &p_dyn, Array &p_stat, Array &t, Array &u, Array &v, Array &w, Array &c, Array &aux_u, Array &aux_v, Array &Salt_Finger, Array &Salt_Diffusion, Array &Buoyancy_Force, Array &Salt_Balance, Array_2D &Upwelling, Array_2D &Downwelling, Array_2D &SaltFinger, Array_2D &SaltDiffusion, Array_2D &BuoyancyForce, Array_2D &BottomWater )
{
	double x, y, z, dx, dy;

	j_max = jm;
	k_max = km;

	// file administration
    string Hydrosphere_radial_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Hyd_radial_" + std::to_string(i_radial) + "_" + std::to_string(pressure_iter) + ".vtk";
	ofstream Hydrosphere_vtk_radial_File;
	Hydrosphere_vtk_radial_File.precision ( 4 );
	Hydrosphere_vtk_radial_File.setf ( ios::fixed );
    Hydrosphere_vtk_radial_File.open(Hydrosphere_radial_File_Name);

	if (!Hydrosphere_vtk_radial_File.is_open()) {
        cerr << "ERROR: could not open paraview_vtk file " << __FILE__ << " at line " << __LINE__ << "\n";
		abort();
	}

	// begin writing
	Hydrosphere_vtk_radial_File <<  "# vtk DataFile Version 3.0" << endl;
	Hydrosphere_vtk_radial_File <<  "Radial_Data_Hydrosphere_Circulation\n";
	Hydrosphere_vtk_radial_File <<  "ASCII" << endl;
	Hydrosphere_vtk_radial_File <<  "DATASET STRUCTURED_GRID" << endl;
	Hydrosphere_vtk_radial_File <<  "DIMENSIONS " << k_max << " "<< j_max << " " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "POINTS " << j_max * k_max << " float" << endl;

	// transformation from spherical to cartesian coordinates
	x = 0.;
	y = 0.;
	z = 0.;

	dx = .1;
	dy = .1;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( k == 0 ) y = 0.;
			else y = y + dy;

			Hydrosphere_vtk_radial_File << x << " " << y << " "<< z << endl;
		}
		y = 0.;
		x = x + dx;
	}

	Hydrosphere_vtk_radial_File <<  "POINT_DATA " << j_max * k_max << endl;



    dump_radial("Topography", h, 1., i_radial, Hydrosphere_vtk_radial_File);

    dump_radial("u-Component", u, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("v-Component", v, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("w-Component", w, 1., i_radial, Hydrosphere_vtk_radial_File);

    dump_radial("Temperature", t, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("Salinity", c, 1., i_radial, Hydrosphere_vtk_radial_File);

    dump_radial("SaltFinger", Salt_Finger, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("SaltDiffusion", Salt_Diffusion, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("SaltBalance", Salt_Balance, 1., i_radial, Hydrosphere_vtk_radial_File);
    dump_radial("BuoyancyForce", Buoyancy_Force, 1., i_radial, Hydrosphere_vtk_radial_File);

    dump_radial("PressureDynamic", p_dyn, u_0 * u_0 * r_0_water, i_radial, Hydrosphere_vtk_radial_File);
    // dump_radial("PressureStatic", p_stat, 1., i_radial, Hydrosphere_vtk_radial_File);

    dump_radial_2d("Upwelling", Upwelling, 1., Hydrosphere_vtk_radial_File);
    dump_radial_2d("Downwelling", Downwelling, 1., Hydrosphere_vtk_radial_File);
    dump_radial_2d("BottomWater", BottomWater, 1., Hydrosphere_vtk_radial_File);






/*
	// writing u-component
	Hydrosphere_vtk_radial_File <<  "SCALARS u-Component float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << u.x[ i_radial ][ j ][ k ] << endl;
		}
	}

	// writing v-component
	Hydrosphere_vtk_radial_File <<  "SCALARS v-Component float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] << endl;
		}
	}

	// writing w-component
	Hydrosphere_vtk_radial_File <<  "SCALARS w-Component float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" <<endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << w.x[ i_radial ][ j ][ k ] << endl;
		}
	}

	// writing of temperature
	Hydrosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << t.x[ i_radial ][ j ][ k ] << endl;
		}
	}

	// writing salinity
	Hydrosphere_vtk_radial_File <<  "SCALARS Salinity float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << c.x[ i_radial ][ j ][ k ] << endl;
		}
	}


	// writing salt finger
	Hydrosphere_vtk_radial_File <<  "SCALARS SaltFinger float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << Salt_Finger.x[ i_radial ][ j ][ k ] << endl;
		}
	}

	// writing salt diffusion
	Hydrosphere_vtk_radial_File <<  "SCALARS SaltDiffusion float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << Salt_Diffusion.x[ i_radial ][ j ][ k ] << endl;
		}
	}

	// writing balance
	Hydrosphere_vtk_radial_File <<  "SCALARS SaltBalance float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << Salt_Balance.x[ i_radial ][ j ][ k ] << endl;
		}
	}

	// writing pressure
	Hydrosphere_vtk_radial_File <<  "SCALARS Pressure_dynamic float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

	for ( int j = 1; j < jm-1; j++ )
	{
		for ( int k = 1; k < km-1; k++ )
		{
			if ( fabs ( p_dyn.x[ i_radial ][ j ][ k ] ) > max_p_dyn )
			{
				max_p_dyn = fabs ( p_dyn.x[ i_radial ][ j ][ k ] );
				if ( max_p_dyn == 0. ) max_p_dyn = 1.e-6;
			}
		}
	}

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
//				Hydrosphere_vtk_radial_File << p_dyn.x[ i_radial ][ j ][ k ] / max_p_dyn << endl;
			Hydrosphere_vtk_radial_File << p_dyn.x[ i_radial ][ j ][ k ] << endl;
		}
	}


// writing static pressure
		Hydrosphere_vtk_radial_File <<  "SCALARS Pressure_static float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( p_stat.x[ i_radial ][ j ][ k ] ) > max_p_stat )
				{
					max_p_stat = fabs ( p_stat.x[ i_radial ][ j ][ k ] );
					if ( max_p_stat == 0. ) max_p_stat = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Hydrosphere_vtk_radial_File << p_stat.x[ i_radial ][ j ][ k ] / max_p_stat << endl;
				Hydrosphere_vtk_radial_File << p_stat.x[ i_radial ][ j ][ k ] << endl;
			}
		}


	// writing upwelling
	Hydrosphere_vtk_radial_File <<  "SCALARS Upwelling float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << Upwelling.y[ j ][ k ] << endl;
		}
	}

	// writing downwelling
	Hydrosphere_vtk_radial_File <<  "SCALARS Downwelling float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << Downwelling.y[ j ][ k ] << endl;
		}
	}

	// writing BottomWater
	Hydrosphere_vtk_radial_File <<  "SCALARS BottomWater float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << BottomWater.y[ j ][ k ] << endl;
		}
	}

	// writing salt finger
	Hydrosphere_vtk_radial_File <<  "SCALARS SaltFinger float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << SaltFinger.y[ j ][ k ] << endl;
		}
	}

	// writing salt diffusion
	Hydrosphere_vtk_radial_File <<  "SCALARS SaltDiffusion float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << SaltDiffusion.y[ j ][ k ] << endl;
		}
	}

	// writing bouyancy force
	Hydrosphere_vtk_radial_File <<  "SCALARS BuoyancyForce float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << BuoyancyForce.y[ j ][ k ] << endl;
		}
	}

	// writing sea ground
	Hydrosphere_vtk_radial_File <<  "SCALARS Topography float " << 1 << endl;
	Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << h.x[ i_radial ][ j ][ k ] << endl;
		}
	}
*/
	// writing zonal v-w cell structure
	Hydrosphere_vtk_radial_File <<  "VECTORS v-w-Cell float" << endl;

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Hydrosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] << " " << w.x[ i_radial ][ j ][ k ] << " " << z << endl;
		}
	}

	// end writing

	Hydrosphere_vtk_radial_File.close();
}







void PostProcess_Hydrosphere::paraview_vtk_zonal ( const string &Name_Bathymetry_File, int &k_zonal, int &pressure_iter, double &u_0, double &r_0_water, Array &h, Array &p_dyn, Array &p_stat, Array &t, Array &u, Array &v, Array &w, Array &c, Array &Salt_Finger, Array &Salt_Diffusion, Array &Buoyancy_Force, Array &Salt_Balance )
{
	double x, y, z, dx, dy;

	i_max = im;
	j_max = jm;

	// file administration
	streampos anfangpos, endpos;

    string Hydrosphere_zongal_File_Name = output_path + "/[" + Name_Bathymetry_File + "]_Hyd_zonal_" + std::to_string(k_zonal) + "_" + std::to_string(pressure_iter) + ".vtk";
//	string Hydrosphere_zonal_File_Name <<  "/[" << Name_Bathymetry_File << "]_Hyd_zonal_" << k_zonal << "_" << pressure_iter << ".vtk";
	ofstream Hydrosphere_vtk_zonal_File;
	Hydrosphere_vtk_zonal_File.precision ( 4 );
	Hydrosphere_vtk_zonal_File.setf ( ios::fixed );
    Hydrosphere_vtk_zonal_File.open ( Hydrosphere_zongal_File_Name);

	if (!Hydrosphere_vtk_zonal_File.is_open()) {
        cerr << "ERROR: could not open vtk_zonal file " << __FILE__ << " at line " << __LINE__ << "\n";
		abort();
	}

	// begin writing
	Hydrosphere_vtk_zonal_File <<  "# vtk DataFile Version 3.0" << endl;
	Hydrosphere_vtk_zonal_File <<  "Zonal_Data_Hydrosphere_Circulation\n";
	Hydrosphere_vtk_zonal_File <<  "ASCII" << endl;
	Hydrosphere_vtk_zonal_File <<  "DATASET STRUCTURED_GRID" << endl;
	Hydrosphere_vtk_zonal_File <<  "DIMENSIONS " << j_max << " "<< i_max << " " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "POINTS " << i_max * j_max << " float" << endl;

	// transformation from spherical to cartesian coordinates
	x = 0.;
	y = 0.;
	z = 0.;

	dx = .1;
	dy = .05;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			if ( j == 0 ) y = 0.;
			else y = y + dy;

			Hydrosphere_vtk_zonal_File << x << " " << y << " "<< z << endl;
		}
		y = 0.;
		x = x + dx;
	}


	Hydrosphere_vtk_zonal_File <<  "POINT_DATA " << i_max * j_max << endl;




    dump_zonal("Topography", h, 1., k_zonal, Hydrosphere_vtk_zonal_File);

    dump_zonal("u-Component", u, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("v-Component", v, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("w-Component", w, 1., k_zonal, Hydrosphere_vtk_zonal_File);

    dump_zonal("Temperature", t, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("Salinity", c, 1., k_zonal, Hydrosphere_vtk_zonal_File);

    dump_zonal("Salinity", c, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("SaltFinger", Salt_Finger, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("SaltDiffusion", Salt_Diffusion, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("SaltBalance", Salt_Balance, 1., k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("BuoyancyForce", Buoyancy_Force, 1., k_zonal, Hydrosphere_vtk_zonal_File);

    dump_zonal("PressureDynamic", p_dyn, u_0 * u_0 * r_0_water, k_zonal, Hydrosphere_vtk_zonal_File);
    dump_zonal("PressureStatic", p_stat, 1., k_zonal, Hydrosphere_vtk_zonal_File);





/*
	// writing u-component
	Hydrosphere_vtk_zonal_File <<  "SCALARS u-Component float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing v-component
	Hydrosphere_vtk_zonal_File <<  "SCALARS v-Component float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << v.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing w-component
	Hydrosphere_vtk_zonal_File <<  "SCALARS w-Component float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" <<endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << w.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing of temperature
	Hydrosphere_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << t.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing salinity
	Hydrosphere_vtk_zonal_File <<  "SCALARS Salinity float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << c.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing salt finger
	Hydrosphere_vtk_zonal_File <<  "SCALARS SaltFinger float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << Salt_Finger.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing salt diffusion
	Hydrosphere_vtk_zonal_File <<  "SCALARS SaltDiffusion float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << Salt_Diffusion.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing salt balance
	Hydrosphere_vtk_zonal_File <<  "SCALARS SaltBalance float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << Salt_Balance.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing buoancy force
	Hydrosphere_vtk_zonal_File <<  "SCALARS BuoyancyForce float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << Buoyancy_Force.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing salt balance
	Hydrosphere_vtk_zonal_File <<  "SCALARS SaltBalance float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << Salt_Balance.x[ i ][ j ][ k_zonal ] << endl;
		}
	}

	// writing dynamic pressure
	Hydrosphere_vtk_zonal_File <<  "SCALARS Pressure_dynamic float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
	{
		for ( int j = 1; j < jm-1; j++ )
		{
			if ( fabs ( p_dyn.x[ i ][ j ][ k_zonal ] ) > max_p_dyn )
			{
				max_p_dyn = fabs ( p_dyn.x[ i ][ j ][ k_zonal ] );
				if ( max_p_dyn == 0. ) max_p_dyn = 1.e-6;
			}
		}
	}

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
//				Hydrosphere_vtk_zonal_File << p_dyn.x[ i ][ j ][ k_zonal ] / max_p_dyn << endl;
			Hydrosphere_vtk_zonal_File << p_dyn.x[ i ][ j ][ k_zonal ] << endl;
		}
	}


// writing static pressure
		Hydrosphere_vtk_zonal_File <<  "SCALARS Pressure_static float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( p_stat.x[ i ][ j ][ k_zonal ] ) > max_p_stat )
				{
					max_p_stat = fabs ( p_stat.x[ i ][ j ][ k_zonal ] );
					if ( max_p_stat == 0. ) max_p_stat = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Hydrosphere_vtk_zonal_File << p_stat.x[ i ][ j ][ k_zonal ] / max_p_stat << endl;
				Hydrosphere_vtk_zonal_File << p_stat.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


	// writing sea ground
	Hydrosphere_vtk_zonal_File <<  "SCALARS Topography float " << 1 << endl;
	Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << h.x[ i ][ j ][ k_zonal ] << endl;
		}
	}
*/
	// writing zonal u-v cell structure
	Hydrosphere_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			Hydrosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] << " " << v.x[ i ][ j ][ k_zonal ] << " " << z << endl;
		}
	}

	// end writing

	Hydrosphere_vtk_zonal_File.close();
}






void PostProcess_Hydrosphere::Atmosphere_TransferFile_read ( const string &Name_Bathymetry_File, Array &v, Array &w, Array &p )
{

	// file administration
	ifstream v_w_Transfer_File;
    string Name_v_w_Transfer_File = input_path + "/[" + Name_Bathymetry_File + "]_Transfer_Atm.vw";
    v_w_Transfer_File.precision(4);
    v_w_Transfer_File.setf(ios::fixed);
	v_w_Transfer_File.open(Name_v_w_Transfer_File);

	if (!v_w_Transfer_File.is_open()) {
        cout << "ERROR: transfer file name in hydrosphere: " << Name_v_w_Transfer_File << "\n";
        cerr << "ERROR: could not open transfer file " << __FILE__ << " at line " << __LINE__ << "\n";
		abort();
	}

	// begin reading
	for ( int j = 0; j < jm; j++ ) {
		for ( int   k = 0; k < km; k++ ) {
			v_w_Transfer_File >> v.x[ im-1 ][ j ][ k ];
			v_w_Transfer_File >> w.x[ im-1 ][ j ][ k ];
			v_w_Transfer_File >> p.x[ im-1 ][ j ][ k ];
		}
	}

	// final file administration

	v_w_Transfer_File.close();
}





void PostProcess_Hydrosphere::Hydrosphere_PlotData ( const string &Name_Bathymetry_File, double &u_0, Array &h, Array &v, Array &w, Array &t, Array &c, Array_2D &BottomWater, Array_2D & Upwelling, Array_2D & Downwelling )
{

	// file administration
	ofstream PlotData_File;
	PlotData_File.precision ( 4 );
	PlotData_File.setf ( ios::fixed );
	string path = output_path + "/[" + Name_Bathymetry_File + "]_PlotData_Hyd.xyz";
	PlotData_File.open (path);

	if (!PlotData_File.is_open()) {
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
		abort();
	}

	PlotData_File << "lons (deg)" << ", " << "lats (deg)" << ", " << "topography" << ", " << "v-velocity (m/s)" << ", " << "w-velocity (m/s)" << ", " << "velocity-mag (m/s)" << ", " << "temperature (Celsius)" << ", " << "salinity (psu)" << ", " << "bottom_water (m/s)" << ", " <<  "upwelling (m/s)" << ", " <<  "downwelling (m/s)" << endl;

	double vel_mag;

	for ( int k = 0; k < km; k++ ) {
		for ( int j = 0; j < jm; j++ ) {
			vel_mag = sqrt ( pow ( v.x[ im-1 ][ j ][ k ] * u_0 , 2 ) + pow ( w.x[ im-1 ][ j ][ k ] * u_0, 2 ) );
			PlotData_File << k << " " << j << " " << h.x[ im-1 ][ j ][ k ] << " " << v.x[ im-1 ][ j ][ k ] * u_0 << " " << w.x[ im-1 ][ j ][ k ] * u_0 << " " << vel_mag << " " << t.x[ im-1 ][ j ][ k ] * 273.15 - 273.15 << " " << c.x[ im-1 ][ j ][ k ] * 35. << " " << BottomWater.y[ j ][ k ] << " " << Upwelling.y[ j ][ k ] << "   " << Downwelling.y[ j ][ k ] << " " <<  endl;
		}
	}

	PlotData_File.close();
}
