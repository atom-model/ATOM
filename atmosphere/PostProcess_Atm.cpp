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



PostProcess_Atmosphere::PostProcess_Atmosphere(int im, int jm, int km) {
	this->im = im;
	this->jm = jm;
	this->km = km;
}

PostProcess_Atmosphere::~PostProcess_Atmosphere() {}

void PostProcess_Atmosphere::Atmosphere_SequelFile_write(string &output_path, string &Name_Bathymetry_File, int &n, double &time, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &co2n) {
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

	for ( int i = 0; i < im; i++ )
	{
		Sequel_File << rad.z[ i ] << endl;
	}


	for ( int j = 0; j < jm; j++ )
	{
		Sequel_File << the.z[ j ] << endl;
	}


	for ( int k = 0; k < km; k++ )
	{
		Sequel_File << phi.z[ k ] << endl;
	}


	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Sequel_File << u.x[ i ][ j ][ k ] << " " << v.x[ i ][ j ][ k ] << " " << w.x[ i ][ j ][ k ]  << endl;
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Sequel_File << un.x[ i ][ j ][ k ] << " " << vn.x[ i ][ j ][ k ] << " " << wn.x[ i ][ j ][ k ]  << endl;
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Sequel_File << t.x[ i ][ j ][ k ] << " " << tn.x[ i ][ j ][ k ]  << endl;
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Sequel_File << c.x[ i ][ j ][ k ] << " " << cn.x[ i ][ j ][ k ]  << endl;
			}
		}
	}

	// end writing
	Sequel_File.close();
}

void PostProcess_Atmosphere::Atmosphere_SequelFile_read(string &output_path, string &Name_Bathymetry_File, int &n, double &time, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &co2n) {
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

	for ( int i = 0; i < im; i++ )
	{
		Sequel_File >> rad.z[ i ];
	}


	for ( int j = 0; j < jm; j++ )
	{
		Sequel_File >> the.z[ j ];
	}


	for ( int k = 0; k < km; k++ )
	{
		Sequel_File >> phi.z[ k ];
	}


	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int   k = 0; k < km; k++ )
			{
				Sequel_File >> u.x[ i ][ j ][ k ];
				Sequel_File >> v.x[ i ][ j ][ k ];
				Sequel_File >> w.x[ i ][ j ][ k ];
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Sequel_File >> un.x[ i ][ j ][ k ];
				Sequel_File >> vn.x[ i ][ j ][ k ];
				Sequel_File >> wn.x[ i ][ j ][ k ];
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
				Sequel_File >> t.x[ i ][ j ][ k ];
				Sequel_File >> tn.x[ i ][ j ][ k ];
			}
		}
	}


	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int i = 0; i < im; i++ )
			{
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
	stringstream Name_v_w_Transfer_File;

	streampos anfangpos, endpos;

// file administration
	Name_v_w_Transfer_File << "[" << Name_Bathymetry_File << "]_Transfer_Atm.vw";
	ofstream v_w_Transfer_File;
	v_w_Transfer_File.precision ( 4 );
	v_w_Transfer_File.setf ( ios::fixed );
	v_w_Transfer_File.open ( Name_v_w_Transfer_File.str().c_str(), ios_base::out );

	if ( v_w_Transfer_File.good() )
	{
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: could be opened" << endl;
		cout << endl;
		anfangpos = v_w_Transfer_File.tellp();
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: starts at ::::::: " << anfangpos << endl;


	// begin writing

		cout << "\n\n”***** Atmosphere_v_w_Transfer_File_write:   begin of writing!" << endl << endl;



		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				v_w_Transfer_File << v.x[ 0 ][ j ][ k ] << " " << w.x[ 0 ][ j ][ k ] << " " << p_dyn.x[ 0 ][ j ][ k ]  << endl;
			}
		}



// end writing


// final file administration

		endpos = v_w_Transfer_File.tellp();
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
		cout << endl;
	}


// in case of failing

	if ( v_w_Transfer_File.fail() )
	{
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: is lost!" << endl;
		return;
	}

	v_w_Transfer_File.close();

	if ( v_w_Transfer_File.good() )
	{
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: could be closed after writing" << endl;
		cout << endl;
	}

	if ( v_w_Transfer_File.fail() )
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: could not be closed properly!" << endl << endl << endl;


return;
}






void PostProcess_Atmosphere::paraview_vts ( string &Name_Bathymetry_File, int &n, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &p_dyn, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Rain, Array &Ice, Array &Rain_super, Array &IceLayer )
{
	double x, y, z, sinthe, sinphi, costhe, cosphi;

	stringstream Atmosphere_vts_File_Name;




// file administration
	streampos anfangpos, endpos;

	Atmosphere_vts_File_Name << "[" << Name_Bathymetry_File << "]_Atm" << n << ".vts";
	ofstream Atmosphere_vts_File;
	Atmosphere_vts_File.precision ( 4 );
	Atmosphere_vts_File.setf ( ios::fixed );
	Atmosphere_vts_File.open ( Atmosphere_vts_File_Name.str().c_str(), ios_base::out );

	if ( Atmosphere_vts_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_vts_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Atmosphere_vts_File.tellp();
		cout << "***** file ::::: " << Atmosphere_vts_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;




// begin writing

		Atmosphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
		Atmosphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
		Atmosphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
		Atmosphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

//		Atmosphere_vts_File <<  "   <PointData Vectors=\"Velocity Rotation\" Scalars=\"Topography Temperature Pressure WaterVapour u-Component v-Component w-Component\">\n"  << endl;
		Atmosphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography Temperature Pressure WaterVapour\">\n"  << endl;




// writing u, v und w velocity components in cartesian coordinates

		Atmosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;


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

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
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

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
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

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
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

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
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

// end writing


// final file administration

		endpos = Atmosphere_vts_File.tellp();
		cout << "***** file ::::: " << Atmosphere_vts_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Atmosphere_vts_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
	}


// in case of failing

	if ( Atmosphere_vts_File.fail() )
	{
		cout << "***** file ::::: " << Atmosphere_vts_File_Name.str() << " ::::: is lost!" << endl;
		return;
	}

	Atmosphere_vts_File.close();

	if ( Atmosphere_vts_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_vts_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Atmosphere_vts_File.fail() )
		cout << "***** file ::::: " << Atmosphere_vts_File_Name.str() << " ::::: could not be closed properly!" << endl;
return;
}





void PostProcess_Atmosphere::paraview_panorama_vts (string &Name_Bathymetry_File, int &pressure_iter, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, double &co2_0, Array &h, Array &t, Array &p_dyn, Array &p_stat, Array &BuoyancyForce, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Q_Sensible, Array &IceLayer, Array &epsilon_3D, Array &P_rain, Array &P_snow )
{
	double x, y, z, dx, dy, dz;

	stringstream Atmosphere_panorama_vts_File_Name;

	max_u = max_v = max_w = max_t = max_c = max_co2 = max_cloud = max_ice = max_P_rain = max_P_snow = max_P_conv = max_M_u = max_M_d = max_p_dyn = max_p_stat = 0.;
	max_Rain = max_Rain_super = max_Ice = max_Latency = max_Q_Sensible = max_Precipitation = 0.;
	max_t_Evaporation = max_t_Condensation = max_t_evap_3D = max_t_cond_3D = 0.;
	max_precipitable_water = max_IceAir = max_Q_bottom = max_Q_latent = max_Q_sensible = 0.;
	max_t_Evaporation_Penman = max_t_Evaporation_Haude = max_Q_Radiation = max_buoyancy_force = 0.;
	max_Q_t_Evaporation = max_precipitation_NASA = max_Water = max_Water_super = max_Vegetation = max_IceLayer = 0.;


// file administration
	streampos anfangpos, endpos;

	Atmosphere_panorama_vts_File_Name <<  "[" << Name_Bathymetry_File << "]_Atm_panorama_" << pressure_iter << ".vts";
	ofstream Atmosphere_panorama_vts_File;
	Atmosphere_panorama_vts_File.precision ( 4 );
	Atmosphere_panorama_vts_File.setf ( ios::fixed );
	Atmosphere_panorama_vts_File.open ( Atmosphere_panorama_vts_File_Name.str().c_str(), ios_base::out );

	if ( Atmosphere_panorama_vts_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_panorama_vts_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Atmosphere_panorama_vts_File.tellp();
		cout << "***** file ::::: " << Atmosphere_panorama_vts_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;


// begin writing
		Atmosphere_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
		Atmosphere_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
		Atmosphere_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
		Atmosphere_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

//		Atmosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Topography Temperature Pressure WaterVapour Latency Rain Ice RainSuper IceLayer\">\n"  << endl;
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

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
// transformtion from spherical to cartesian coordinates for representation in ParaView
//					Atmosphere_panorama_vts_File << u.x[ i ][ j ][ k ] / max_u << " " << v.x[ i ][ j ][ k ] / max_v << " " << w.x[ i ][ j ][ k ] / max_w << endl;
					Atmosphere_panorama_vts_File << u.x[ i ][ j ][ k ] << " " << v.x[ i ][ j ][ k ] << " " << w.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing of mountain surface
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Topography\" format=\"ascii\">\n"  << endl;

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					Atmosphere_panorama_vts_File << h.x[ i ][ j ][ k ]  << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing of temperature
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( t.x[ i ][ j ][ k ] ) > max_t ) 
					{
						max_t = fabs ( t.x[ i ][ j ][ k ] );
						if ( max_t == 0. ) max_t = 1.e-6;
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
//					Atmosphere_panorama_vts_File << t.x[ i ][ j ][ k ] / max_t << endl;
					Atmosphere_panorama_vts_File << t.x[ i ][ j ][ k ] * t_0 - t_0 << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


/*
// writing of condensation temperature
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CondensationTemp\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( t_cond_3D.x[ i ][ j ][ k ] ) > max_t ) 
					{
						max_t = fabs ( t_cond_3D.x[ i ][ j ][ k ] );
						if ( max_t == 0. ) max_t = 1.e-6;
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
//					Atmosphere_panorama_vts_File << t_cond_3D.x[ i ][ j ][ k ] / max_t << endl;
					Atmosphere_panorama_vts_File << t_cond_3D.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;



// writing of evaporation temperature
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"EvaporationTemp\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( t_evap_3D.x[ i ][ j ][ k ] ) > max_t ) 
					{
						max_t = fabs ( t_evap_3D.x[ i ][ j ][ k ] );
						if ( max_t == 0. ) max_t = 1.e-6;
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
//					Atmosphere_panorama_vts_File << t_evap_3D.x[ i ][ j ][ k ] / max_t << endl;
					Atmosphere_panorama_vts_File << t_evap_3D.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;
*/




// writing of epsilon_3D
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Epsilon_3D\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( epsilon_3D.x[ i ][ j ][ k ] ) > max_t ) 
					{
						max_t = fabs ( epsilon_3D.x[ i ][ j ][ k ] );
						if ( max_t == 0. ) max_t = 1.e-6;
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
//					Atmosphere_panorama_vts_File << epsilon_3D.x[ i ][ j ][ k ] / max_t << endl;
					Atmosphere_panorama_vts_File << epsilon_3D.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;





// writing dynamic pressure
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureDynamic\" format=\"ascii\">\n"  << endl;

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
//					Atmosphere_panorama_vts_File << p_dyn.x[ i ][ j ][ k ] / max_p_dyn << endl;
					Atmosphere_panorama_vts_File << p_dyn.x[ i ][ j ][ k ] * u_0 * u_0 * r_air << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;

/*
// writing static pressure
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureStatic\" format=\"ascii\">\n"  << endl;

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
//					Atmosphere_panorama_vts_File << p_stat.x[ i ][ j ][ k ] / max_p_stat << endl;
					Atmosphere_panorama_vts_File << p_stat.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;
*/

// writing buoyancy force
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"BuoyancyForce\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( BuoyancyForce.x[ i ][ j ][ k ] ) > max_buoyancy_force ) 
					{
						max_buoyancy_force = fabs ( BuoyancyForce.x[ i ][ j ][ k ] );
						if ( max_buoyancy_force == 0. ) max_buoyancy_force = 1.e-6;
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
//					Atmosphere_panorama_vts_File << BuoyancyForce.x[ i ][ j ][ k ] / max_p_stat << endl;
					Atmosphere_panorama_vts_File << BuoyancyForce.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing water vapour
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"WaterVapour\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( c.x[ i ][ j ][ k ] ) > max_c ) 
					{
						max_c = fabs ( c.x[ i ][ j ][ k ] );
						if ( max_c == 0. ) max_c = 1.e-6;
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
//					Atmosphere_panorama_vts_File << c.x[ i ][ j ][ k ] / max_c << endl;
					Atmosphere_panorama_vts_File << c.x[ i ][ j ][ k ] * 1000. << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing cloud water
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CloudWater\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( cloud.x[ i ][ j ][ k ] ) > max_cloud ) 
					{
						max_cloud = fabs ( cloud.x[ i ][ j ][ k ] );
						if ( max_cloud == 0. ) max_cloud = 1.e-6;
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
//					Atmosphere_panorama_vts_File << cloud.x[ i ][ j ][ k ] / max_cloud << endl;
					Atmosphere_panorama_vts_File << cloud.x[ i ][ j ][ k ] * 1000. << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing cloud ice
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CloudIce\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( ice.x[ i ][ j ][ k ] ) > max_ice ) 
					{
						max_ice = fabs ( ice.x[ i ][ j ][ k ] );
						if ( max_ice == 0. ) max_ice = 1.e-6;
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
//					Atmosphere_panorama_vts_File << ice.x[ i ][ j ][ k ] / max_ice << endl;
					Atmosphere_panorama_vts_File << ice.x[ i ][ j ][ k ] * 1000. << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing P_rain
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PrecipitationRain\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( P_rain.x[ i ][ j ][ k ] ) > max_P_rain ) 
					{
						max_P_rain = fabs ( P_rain.x[ i ][ j ][ k ] );
						if ( max_P_rain == 0. ) max_P_rain = 1.e-6;
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
//					Atmosphere_panorama_vts_File << P_rain.x[ i ][ j ][ k ] / max_P_rain << endl;
					Atmosphere_panorama_vts_File << P_rain.x[ i ][ j ][ k ] * 1000. << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing P_snow
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PrecipitationSnow\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					max_P_snow = fabs ( P_snow.x[ i ][ j ][ k ] );
					if ( max_P_snow == 0. ) max_P_snow = 1.e-6;
				}
			}
		}

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
//					Atmosphere_panorama_vts_File << P_snow.x[ i ][ j ][ k ] / max_P_snow << endl;
					Atmosphere_panorama_vts_File << P_snow.x[ i ][ j ][ k ] * 1000. << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;




/*
// writing M_u
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Updraft\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					max_M_u = fabs ( M_u.x[ i ][ j ][ k ] );
					if ( max_M_u == 0. ) max_M_u = 1.e-6;
				}
			}
		}

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
//					Atmosphere_panorama_vts_File << M_u.x[ i ][ j ][ k ] / max_M_u << endl;
					Atmosphere_panorama_vts_File << M_u.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing M_d
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Downdraft\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					max_M_d = fabs ( M_d.x[ i ][ j ][ k ] );
					if ( max_M_d == 0. ) max_M_d = 1.e-6;
				}
			}
		}

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
//					Atmosphere_panorama_vts_File << M_d.x[ i ][ j ][ k ] / max_M_d << endl;
					Atmosphere_panorama_vts_File << M_d.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;



// writing P_co2_nv
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PrecipitationConv\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					max_P_conv = fabs ( P_co2_nv.x[ i ][ j ][ k ] );
					if ( max_P_conv == 0. ) max_P_conv = 1.e-6;
				}
			}
		}

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
//					Atmosphere_panorama_vts_File << P_co2_nv.x[ i ][ j ][ k ] / max_P_conv << endl;
					Atmosphere_panorama_vts_File << P_co2_nv.x[ i ][ j ][ k ] * 1000. << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;
*/



// writing co2
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CO2-Concentration\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( co2.x[ i ][ j ][ k ] ) > max_co2 ) 
					{
						max_co2 = fabs ( co2.x[ i ][ j ][ k ] );
						if ( max_co2 == 0. ) max_co2 = 1.e-6;
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
//					Atmosphere_panorama_vts_File << co2.x[ i ][ j ][ k ] / max_co2 << endl;
					Atmosphere_panorama_vts_File << co2.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing Latency
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Latency\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( Latency.x[ i ][ j ][ k ] ) > max_Latency ) 
					{
						max_Latency = fabs ( Latency.x[ i ][ j ][ k ] );
						if ( max_Latency == 0. ) max_Latency = 1.e-6;
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
//					Atmosphere_panorama_vts_File << Latency.x[ i ][ j ][ k ] / max_Latency << endl;
					Atmosphere_panorama_vts_File << Latency.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing Q_Sensible
		Atmosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Q_Sensible\" format=\"ascii\">\n"  << endl;

		for ( int k = 1; k < km-1; k++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int i = 1; i < im-1; i++ )
				{
					if ( fabs ( Q_Sensible.x[ i ][ j ][ k ] ) > max_Q_Sensible ) 
					{
						max_Q_Sensible = fabs ( Q_Sensible.x[ i ][ j ][ k ] );
						if ( max_Q_Sensible == 0. ) max_Q_Sensible = 1.e-6;
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
//					Atmosphere_panorama_vts_File << Q_Sensible.x[ i ][ j ][ k ] / max_Q_Sensible << endl;
					Atmosphere_panorama_vts_File << Q_Sensible.x[ i ][ j ][ k ] << endl;
				}
				Atmosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Atmosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Atmosphere_panorama_vts_File <<  "\n"  << endl;
		Atmosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;





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

// end writing


// final file administration

		endpos = Atmosphere_panorama_vts_File.tellp();
		cout << "***** file ::::: " << Atmosphere_panorama_vts_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Atmosphere_panorama_vts_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
	}


// in case of failing

	if ( Atmosphere_panorama_vts_File.fail() )
	{
		cout << "***** file ::::: " << Atmosphere_panorama_vts_File_Name.str() << " ::::: is lost!" << endl;
		return;
	}

	Atmosphere_panorama_vts_File.close();

	if ( Atmosphere_panorama_vts_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_panorama_vts_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Atmosphere_panorama_vts_File.fail() )
		cout << "***** file ::::: " << Atmosphere_panorama_vts_File_Name.str() << " ::::: could not be closed properly!" << endl;
return;
}







void PostProcess_Atmosphere::paraview_vtk_radial ( string &Name_Bathymetry_File, int &i_radial, int &pressure_iter, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, double &co2_0, double &radiation_equator, Array &h, Array &p_dyn, Array &p_stat, Array &t_cond_3D, Array &t_evap_3D, Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Q_Sensible, Array &IceLayer, Array &epsilon_3D, Array &P_rain, Array &P_snow, Array_2D &Evaporation, Array_2D &Condensation, Array_2D &precipitable_water, Array_2D &Q_bottom, Array_2D &Radiation_Balance, Array_2D &Q_Radiation, Array_2D &Q_latent, Array_2D &Q_sensible, Array_2D &Evaporation_Penman, Array_2D &Evaporation_Haude, Array_2D &Q_Evaporation, Array_2D &precipitation_NASA, Array_2D &Vegetation, Array_2D &albedo, Array_2D &epsilon, Array_2D &Precipitation )
{
	double x, y, z, dx, dy;

	stringstream Atmosphere_radial_File_Name;

	j_max = jm;
	k_max = km;

	max_u = max_v = max_w = max_t = max_c = max_co2 = max_cloud = max_ice = max_P_rain = max_P_snow = max_P_conv = max_M_u = max_M_d = max_p_dyn = max_p_stat = 0.;
	max_Rain = max_Rain_super = max_Ice = max_Latency = max_Q_Sensible = max_Precipitation = max_albedo = max_epsilon = 0.;
	max_t_Evaporation = max_t_Condensation = max_t_evap_3D = max_t_cond_3D = 0.;
	max_precipitable_water = max_IceAir = max_Q_bottom = max_Q_latent = max_Q_sensible = 0.;
	max_t_Evaporation_Penman = max_t_Evaporation_Haude = max_Q_Radiation = max_Radiation_Balance = max_buoyancy_force = 0.;
	max_Q_t_Evaporation = max_precipitation_NASA = max_Water = max_Water_super = max_Vegetation = max_IceLayer = max_radiation_3D = 0.;


 
// file administration
	streampos anfangpos, endpos;

	Atmosphere_radial_File_Name << "[" <<  Name_Bathymetry_File << "]_Atm_radial_" << i_radial << "_" << pressure_iter << ".vtk";
	ofstream Atmosphere_vtk_radial_File;
	Atmosphere_vtk_radial_File.precision ( 4 );
	Atmosphere_vtk_radial_File.setf ( ios::fixed );
	Atmosphere_vtk_radial_File.open ( Atmosphere_radial_File_Name.str().c_str(), ios_base::out );

	if ( Atmosphere_vtk_radial_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_radial_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Atmosphere_vtk_radial_File.tellp();
		cout << "***** file ::::: " << Atmosphere_radial_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;


// begin writing
		Atmosphere_vtk_radial_File <<  "# vtk DataFile Version 3.0" << endl;
		Atmosphere_vtk_radial_File <<  "Radial_Data_Atmosphere_Circulation\n";
		Atmosphere_vtk_radial_File <<  "ASCII" << endl;
		Atmosphere_vtk_radial_File <<  "DATASET STRUCTURED_GRID" << endl;
		Atmosphere_vtk_radial_File <<  "DIMENSIONS " << k_max << " "<< j_max << " " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "POINTS " << j_max * k_max << " float" << endl;


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

				Atmosphere_vtk_radial_File << x << " " << y << " "<< z << endl;
			}
			y = 0.;
			x = x + dx;
		}


		Atmosphere_vtk_radial_File <<  "POINT_DATA " << j_max * k_max << endl;

		if ( i_radial != 0 )
		{
// writing u-component
			Atmosphere_vtk_radial_File <<  "SCALARS u-Component float " << 1 << endl;
			Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

			for ( int j = 1; j < jm-1; j++ )
			{
				for ( int k = 1; k < km-1; k++ )
				{
					if ( fabs ( u.x[ i_radial ][ j ][ k ] ) > max_u ) 
					{
						max_u = fabs ( u.x[ i_radial ][ j ][ k ] );
						if ( max_u == 0. ) max_u = 1.e-6;
					}
				}
			}

			for ( int j = 0; j < jm; j++ )
			{
				for ( int k = 0; k < km; k++ )
				{
//					Atmosphere_vtk_radial_File << u.x[ i_radial ][ j ][ k ] / max_u << endl;
					Atmosphere_vtk_radial_File << u.x[ i_radial ][ j ][ k ] << endl;
				}
			}
		}


// writing v-component
		Atmosphere_vtk_radial_File <<  "SCALARS v-Component float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( v.x[ i_radial ][ j ][ k ] ) > max_v ) 
				{
					max_v = fabs ( v.x[ i_radial ][ j ][ k ] );
					if ( max_v == 0. ) max_v = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] / max_v << endl;
				Atmosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing of temperature
		Atmosphere_vtk_radial_File <<  "SCALARS w-Component float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" <<endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( w.x[ i_radial ][ j ][ k ] ) > max_w ) 
				{
					max_w = fabs ( w.x[ i_radial ][ j ][ k ] );
					if ( max_w == 0. ) max_w = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << w.x[ i_radial ][ j ][ k ] / max_w << endl;
				Atmosphere_vtk_radial_File << w.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing temperature
		Atmosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( t.x[ i_radial ][ j ][ k ] ) > max_t ) 
				{
					max_t = fabs ( t.x[ i_radial ][ j ][ k ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << t.x[ i_radial ][ j ][ k ] / max_t << endl;
				Atmosphere_vtk_radial_File << t.x[ i_radial ][ j ][ k ] * t_0 - t_0 << endl;
			}
		}


// writing condensation temperature
		Atmosphere_vtk_radial_File <<  "SCALARS CondensationTemp float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( t_cond_3D.x[ i_radial ][ j ][ k ] ) > max_t ) 
				{
					max_t = fabs ( t_cond_3D.x[ i_radial ][ j ][ k ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << t_cond_3D.x[ i_radial ][ j ][ k ] / max_t << endl;
				Atmosphere_vtk_radial_File << t_cond_3D.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing evaporation temperature
		Atmosphere_vtk_radial_File <<  "SCALARS EvaporationTemp float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( t_evap_3D.x[ i_radial ][ j ][ k ] ) > max_t ) 
				{
					max_t = fabs ( t_evap_3D.x[ i_radial ][ j ][ k ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << t_evap_3D.x[ i_radial ][ j ][ k ] / max_t << endl;
				Atmosphere_vtk_radial_File << t_evap_3D.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing epsilon_3D
		Atmosphere_vtk_radial_File <<  "SCALARS Epsilon_3D float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( epsilon_3D.x[ i_radial ][ j ][ k ] ) > max_t ) 
				{
					max_t = fabs ( epsilon_3D.x[ i_radial ][ j ][ k ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << epsilon_3D.x[ i_radial ][ j ][ k ] / max_t << endl;
				Atmosphere_vtk_radial_File << epsilon_3D.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing water vapour
		Atmosphere_vtk_radial_File <<  "SCALARS WaterVapour float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( c.x[ i_radial ][ j ][ k ] ) > max_c ) 
				{
					max_c = fabs ( c.x[ i_radial ][ j ][ k ] );
					if ( max_c == 0. ) max_c = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << c.x[ i_radial ][ j ][ k ] / max_c << endl;
				Atmosphere_vtk_radial_File << c.x[ i_radial ][ j ][ k ] * 1000. << endl;
			}
		}


// writing cloud water
		Atmosphere_vtk_radial_File <<  "SCALARS CloudWater float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( cloud.x[ i_radial ][ j ][ k ] ) > max_cloud ) 
				{
					max_cloud = fabs ( cloud.x[ i_radial ][ j ][ k ] );
					if ( max_cloud == 0. ) max_cloud = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << cloud.x[ i_radial ][ j ][ k ] / max_cloud << endl;
				Atmosphere_vtk_radial_File << cloud.x[ i_radial ][ j ][ k ] * 1000. << endl;
			}
		}


// writing cloud ice
		Atmosphere_vtk_radial_File <<  "SCALARS CloudIce float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( ice.x[ i_radial ][ j ][ k ] ) > max_ice ) 
				{
					max_ice = fabs ( ice.x[ i_radial ][ j ][ k ] );
					if ( max_ice == 0. ) max_ice = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << ice.x[ i_radial ][ j ][ k ] / max_ice << endl;
				Atmosphere_vtk_radial_File << ice.x[ i_radial ][ j ][ k ] * 1000. << endl;
			}
		}

// writing CO2
		Atmosphere_vtk_radial_File <<  "SCALARS CO2-Concentration float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( co2.x[ i_radial ][ j ][ k ] ) > max_co2 ) 
				{
					max_co2 = fabs ( co2.x[ i_radial ][ j ][ k ] );
					if ( max_co2 == 0. ) max_co2 = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << co2.x[ i_radial ][ j ][ k ] / max_co2 << endl;
				Atmosphere_vtk_radial_File << co2.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing dynamic pressure
		Atmosphere_vtk_radial_File <<  "SCALARS PressureDynamic float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

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
//				Atmosphere_vtk_radial_File << p_dyn.x[ i_radial ][ j ][ k ] / max_p_dyn << endl;
				Atmosphere_vtk_radial_File << p_dyn.x[ i_radial ][ j ][ k ] * u_0 * u_0 * r_air << endl;
			}
		}

/*
// writing static pressure
		Atmosphere_vtk_radial_File <<  "SCALARS PressureStatic float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

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
//				Atmosphere_vtk_radial_File << p_stat.x[ i_radial ][ j ][ k ] / max_p_stat << endl;
				Atmosphere_vtk_radial_File << p_stat.x[ i_radial ][ j ][ k ] << endl;
			}
		}
*/

// writing buoyancy force
		Atmosphere_vtk_radial_File <<  "SCALARS BuoyancyForce float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( BuoyancyForce.x[ i_radial ][ j ][ k ] ) > max_buoyancy_force ) 
				{
					max_buoyancy_force = fabs ( BuoyancyForce.x[ i_radial ][ j ][ k ] );
					if ( max_buoyancy_force == 0. ) max_buoyancy_force = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << BuoyancyForce.x[ i_radial ][ j ][ k ] / max_buoyancy_force << endl;
				Atmosphere_vtk_radial_File << BuoyancyForce.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing topography
		Atmosphere_vtk_radial_File <<  "SCALARS Topography float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Atmosphere_vtk_radial_File << h.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing latent heat
		Atmosphere_vtk_radial_File <<  "SCALARS Latency float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Latency.x[ i_radial ][ j ][ k ] ) >= max_Latency ) 
				{
					max_Latency = fabs ( Latency.x[ i_radial ][ j ][ k ] );
					if ( max_Latency == 0. ) max_Latency = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Latency.x[ i_radial ][ j ][ k ] / max_Latency << endl;
				Atmosphere_vtk_radial_File << Latency.x[ i_radial ][ j ][ k ] << endl;
			}
		}


// writing sensible heat
		Atmosphere_vtk_radial_File <<  "SCALARS Q_Sensible float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Q_Sensible.x[ i_radial ][ j ][ k ] ) >= max_Q_Sensible ) 
				{
					max_Q_Sensible = fabs ( Q_Sensible.x[ i_radial ][ j ][ k ] );
					if ( max_Q_Sensible == 0. ) max_Q_Sensible = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Q_Sensible.x[ i_radial ][ j ][ k ] / max_Q_Sensible << endl;
				Atmosphere_vtk_radial_File << Q_Sensible.x[ i_radial ][ j ][ k ] << endl;
			}
		}


/*
// writing Evaporation
		Atmosphere_vtk_radial_File <<  "SCALARS Evaporation float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Evaporation.y[ j ][ k ] ) > max_t_Evaporation ) 
				{
					max_t_Evaporation = fabs ( Evaporation.y[ j ][ k ] );
					if ( max_t_Evaporation == 0. ) max_t_Evaporation = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Evaporation.y[ j ][ k ] / max_t_Evaporation << endl;
				Atmosphere_vtk_radial_File << Evaporation.y[ j ][ k ] << endl;
			}
		}



// writing Condensation
		Atmosphere_vtk_radial_File <<  "SCALARS Condensation float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Condensation.y[ j ][ k ] ) > max_t_Condensation ) 
				{
					max_t_Condensation = fabs ( Condensation.y[ j ][ k ] );
					if ( max_t_Condensation == 0. ) max_t_Condensation = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Condensation.y[ j ][ k ] / max_t_Condensation << endl;
				Atmosphere_vtk_radial_File << Condensation.y[ j ][ k ] << endl;
			}
		}
*/


// writing precipitation_NASA
		Atmosphere_vtk_radial_File <<  "SCALARS Precipitation_NASA float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( precipitation_NASA.y[ j ][ k ] ) > max_precipitation_NASA ) 
				{
					max_precipitation_NASA = fabs ( precipitation_NASA.y[ j ][ k ] );
					if ( max_precipitation_NASA == 0. ) max_precipitation_NASA = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << precipitation_NASA.y[ j ][ k ] / max_precipitation_NASA << endl;
				Atmosphere_vtk_radial_File << precipitation_NASA.y[ j ][ k ] << endl;
			}
		}



// writing albedo
		Atmosphere_vtk_radial_File <<  "SCALARS albedo float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( albedo.y[ j ][ k ] ) > max_albedo ) 
				{
					max_albedo = fabs ( albedo.y[ j ][ k ] );
					if ( max_albedo == 0. ) max_albedo = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << albedo.y[ j ][ k ] / max_albedo << endl;
				Atmosphere_vtk_radial_File << albedo.y[ j ][ k ] << endl;
			}
		}


// writing epsilon
		Atmosphere_vtk_radial_File <<  "SCALARS epsilon float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( epsilon.y[ j ][ k ] ) > max_epsilon ) 
				{
					max_epsilon = fabs ( epsilon.y[ j ][ k ] );
					if ( max_epsilon == 0. ) max_epsilon = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << epsilon.y[ j ][ k ] / max_epsilon << endl;
				Atmosphere_vtk_radial_File << epsilon.y[ j ][ k ] << endl;
			}
		}


// writing precipitable water
		Atmosphere_vtk_radial_File <<  "SCALARS PrecipitableWater float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( precipitable_water.y[ j ][ k ] ) > max_precipitable_water ) 
				{
					max_precipitable_water = fabs ( precipitable_water.y[ j ][ k ] );
					if ( max_precipitable_water == 0. ) max_precipitable_water = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << precipitable_water.y[ j ][ k ] / max_precipitable_water << endl;
				Atmosphere_vtk_radial_File << precipitable_water.y[ j ][ k ] << endl;
			}
		}


// writing radiation balance
		Atmosphere_vtk_radial_File <<  "SCALARS Radiation_Balance float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Radiation_Balance.y[ j ][ k ] ) > max_Radiation_Balance ) 
				{
					max_Radiation_Balance = fabs ( Radiation_Balance.y[ j ][ k ] );
					if ( max_Radiation_Balance == 0. ) max_Radiation_Balance = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Radiation_Balance.y[ j ][ k ] / max_Radiation_Balance << endl;
				Atmosphere_vtk_radial_File << Radiation_Balance.y[ j ][ k ] << endl;
			}
		}


// writing radiation Radiation_Balance
		Atmosphere_vtk_radial_File <<  "SCALARS Q_Radiation float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Q_Radiation.y[ j ][ k ] ) > max_Q_Radiation ) 
				{
					max_Q_Radiation = fabs ( Q_Radiation.y[ j ][ k ] );
					if ( max_Q_Radiation == 0. ) max_Q_Radiation = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Q_Radiation.y[ j ][ k ] / max_Q_Radiation << endl;
				Atmosphere_vtk_radial_File << Q_Radiation.y[ j ][ k ] << endl;
			}
		}


// writing energy difference
		Atmosphere_vtk_radial_File <<  "SCALARS Q_bottom float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Q_bottom.y[ j ][ k ] ) > max_Q_bottom ) 
				{
					max_Q_bottom = fabs ( Q_bottom.y[ j ][ k ] );
					if ( max_Q_bottom == 0. ) max_Q_bottom = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Q_bottom.y[ j ][ k ] / max_Q_bottom << endl;
				Atmosphere_vtk_radial_File << Q_bottom.y[ j ][ k ] << endl;
			}
		}


// writing Q_latent
		Atmosphere_vtk_radial_File <<  "SCALARS Q_latent  float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Q_latent.y[ j ][ k ] ) > max_Q_latent ) 
				{
					max_Q_latent = fabs ( Q_latent.y[ j ][ k ] );
					if ( max_Q_latent == 0. ) max_Q_latent = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Q_latent.y[ j ][ k ] / max_Q_latent << endl;
				Atmosphere_vtk_radial_File << Q_latent.y[ j ][ k ] << endl;
			}
		}


// writing Q_sensible
		Atmosphere_vtk_radial_File <<  "SCALARS Q_sensible float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Q_sensible.y[ j ][ k ] ) > max_Q_sensible ) 
				{
					max_Q_sensible = fabs ( Q_sensible.y[ j ][ k ] );
					if ( max_Q_sensible == 0. ) max_Q_sensible = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Q_sensible.y[ j ][ k ] / max_Q_sensible << endl;
				Atmosphere_vtk_radial_File << Q_sensible.y[ j ][ k ] << endl;
			}
		}


// writing Evaporation by Penman
		Atmosphere_vtk_radial_File <<  "SCALARS EvaporationPenman float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Evaporation_Penman.y[ j ][ k ] ) > max_t_Evaporation_Penman ) 
				{
					max_t_Evaporation_Penman = fabs ( Evaporation_Penman.y[ j ][ k ] );
					if ( max_t_Evaporation_Penman == 0. ) max_t_Evaporation_Penman = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Evaporation_Penman.y[ j ][ k ] / max_t_Evaporation_Penman << endl;
				Atmosphere_vtk_radial_File << Evaporation_Penman.y[ j ][ k ] << endl;
			}
		}


// writing Evaporation by Haude
/*
		Atmosphere_vtk_radial_File <<  "SCALARS t_EvaporationHaude float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Evaporation_Haude.y[ j ][ k ] ) > max_t_Evaporation_Haude ) 
				{
					max_t_Evaporation_Haude = fabs ( Evaporation_Haude.y[ j ][ k ] );
					if ( max_t_Evaporation_Haude == 0. ) max_t_Evaporation_Haude = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Atmosphere_vtk_radial_File << Evaporation_Haude.y[ j ][ k ] / max_t_Evaporation_Haude << endl;
			}
		}
*/

/*
// writing Evaporation by Q_Evaporation
		Atmosphere_vtk_radial_File <<  "SCALARS Q_Evaporation float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Q_Evaporation.y[ j ][ k ] ) > max_Q_t_Evaporation ) 
				{
					max_Q_t_Evaporation = fabs ( Q_Evaporation.y[ j ][ k ] );
					if ( max_Q_t_Evaporation == 0. ) max_Q_t_Evaporation = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Atmosphere_vtk_radial_File << Q_Evaporation.y[ j ][ k ] / max_Q_t_Evaporation<< endl;
			}
		}
*/





// writing Vegetation
		Atmosphere_vtk_radial_File <<  "SCALARS Vegetation float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Vegetation.y[ j ][ k ] ) > max_Vegetation ) 
				{
					max_Vegetation = fabs ( Vegetation.y[ j ][ k ] );
					if ( max_Vegetation == 0. ) max_Vegetation = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << Vegetation.y[ j ][ k ] / max_Vegetation << endl;
				Atmosphere_vtk_radial_File << Vegetation.y[ j ][ k ] << endl;
			}
		}



// writing P_rain
		Atmosphere_vtk_radial_File <<  "SCALARS PrecipitationRain float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( P_rain.x[ i_radial ][ j ][ k ] ) > max_P_rain ) 
				{
					max_P_rain = fabs ( P_rain.x[ i_radial ][ j ][ k ] );
					if ( max_P_rain == 0. ) max_P_rain = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << P_rain.x[ i_radial ][ j ][ k ] / max_P_rain << endl;
				Atmosphere_vtk_radial_File << P_rain.x[ i_radial ][ j ][ k ] * 1000. << endl;
			}
		}


// writing P_snow
		Atmosphere_vtk_radial_File <<  "SCALARS PrecipitationSnow float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( P_snow.x[ i_radial ][ j ][ k ] ) > max_P_snow ) 
				{
					max_P_snow = fabs ( P_snow.x[ i_radial ][ j ][ k ] );
					if ( max_P_snow == 0. ) max_P_snow = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				if ( t.x[ i_radial ][ j ][ k ] * t_0 - t_0 >= 0. )		P_snow.x[ i_radial ][ j ][ k ] = 0.;
//				Atmosphere_vtk_radial_File << P_snow.x[ i_radial ][ j ][ k ] / max_P_snow << endl;
				Atmosphere_vtk_radial_File << P_snow.x[ i_radial ][ j ][ k ] * 1000. << endl;
			}
		}


// writing precipitation_total

		Atmosphere_vtk_radial_File <<  "SCALARS Precipitation float " << 1 << endl;
		Atmosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Atmosphere_vtk_radial_File << Precipitation.y[ j ][ k ] << endl;
			}
		}


// writing zonal u-v cell structure
		Atmosphere_vtk_radial_File <<  "VECTORS v-w-Cell float " << endl;

		for ( int j = 1; j < jm-1; j++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( v.x[ i_radial ][ j ][ k ] ) > max_v ) 
				{
					max_v = fabs ( v.x[ i_radial ][ j ][ k ] );
					if ( max_v == 0. ) max_v = 1.e-6;
				}
				if ( fabs ( w.x[ i_radial ][ j ][ k ] ) > max_w ) 
				{
					max_w = fabs ( w.x[ i_radial ][ j ][ k ] );
					if ( max_w == 0. ) max_w = 1.e-6;
				}
			}
		}

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] / max_v << " " << w.x[ i_radial ][ j ][ k ] / max_w << " " << z << endl;
				Atmosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] << " " << w.x[ i_radial ][ j ][ k ] << " " << z << endl;
			}
		}






// end writing



// final file administration

		endpos = Atmosphere_vtk_radial_File.tellp();
		cout << "***** file ::::: " << Atmosphere_radial_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Atmosphere_radial_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
	}


// in case of failing

	if ( Atmosphere_vtk_radial_File.fail() )
	{
	  cout << "***** file ::::: " << Atmosphere_radial_File_Name.str() << " ::::: is lost!" << endl;
		return;
	}

	Atmosphere_vtk_radial_File.close();

	if ( Atmosphere_vtk_radial_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_radial_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Atmosphere_vtk_radial_File.fail() )
		cout << "***** file ::::: " << Atmosphere_radial_File_Name.str() << " ::::: could not be closed properly!" << endl;

return;
}




void PostProcess_Atmosphere::paraview_vtk_zonal ( string &Name_Bathymetry_File, int &k_zonal, int &pressure_iter, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, double &co2_0, double &radiation_equator, Array &h, Array &p_dyn, Array &p_stat, Array &t_cond_3D, Array &t_evap_3D, Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Q_Sensible, Array &radiation_3D, Array &epsilon_3D, Array &P_rain, Array &P_snow, Array &S_v, Array &S_c, Array &S_i, Array &S_r, Array &S_s )
{
	double x, y, z, dx, dy;

	stringstream Atmosphere_zonal_File_Name;

	max_u = max_v = max_w = max_t = max_c = max_co2 = max_cloud = max_ice = max_P_rain = max_P_snow = max_P_conv = max_M_u = max_M_d = max_p_dyn = max_p_stat = 0.;
	max_Rain = max_Rain_super = max_Ice = max_Latency = max_Q_Sensible = max_Precipitation = 0.;
	max_t_Evaporation = max_t_Condensation = max_t_evap_3D = max_t_cond_3D = 0.;
	max_precipitable_water = max_IceAir = max_Q_bottom = max_Q_latent = max_Q_sensible = 0.;
	max_t_Evaporation_Penman = max_t_Evaporation_Haude = max_Q_Radiation = max_buoyancy_force = 0.;
	max_Q_t_Evaporation = max_precipitation_NASA = max_Water = max_Water_super = max_Vegetation = max_IceLayer = 0.;

	i_max = im;
	j_max = jm;

// file administration
	streampos anfangpos, endpos;

	Atmosphere_zonal_File_Name <<  "[" << Name_Bathymetry_File << "]_Atm_zonal_" << k_zonal << "_" << pressure_iter << ".vtk";
	ofstream Atmosphere_vtk_zonal_File;
	Atmosphere_vtk_zonal_File.precision ( 4 );
	Atmosphere_vtk_zonal_File.setf ( ios::fixed );
	Atmosphere_vtk_zonal_File.open ( Atmosphere_zonal_File_Name.str().c_str(), ios_base::out );

	if ( Atmosphere_vtk_zonal_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_zonal_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Atmosphere_vtk_zonal_File.tellp();
		cout << "***** file ::::: " << Atmosphere_zonal_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;




// begin writing
		Atmosphere_vtk_zonal_File <<  "# vtk DataFile Version 3.0" << endl;
		Atmosphere_vtk_zonal_File <<  "Zonal_Data_Atmosphere_Circulation\n";
		Atmosphere_vtk_zonal_File <<  "ASCII" << endl;
		Atmosphere_vtk_zonal_File <<  "DATASET STRUCTURED_GRID" << endl;
		Atmosphere_vtk_zonal_File <<  "DIMENSIONS " << j_max << " "<< i_max << " " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "POINTS " << i_max * j_max << " float" << endl;

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

				Atmosphere_vtk_zonal_File << x << " " << y << " "<< z << endl;
			}
			y = 0.;
			x = x + dx;
		}



		Atmosphere_vtk_zonal_File <<  "POINT_DATA " << i_max * j_max << endl;

// writing u-component
		Atmosphere_vtk_zonal_File <<  "SCALARS u-Component float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( u.x[ i ][ j ][ k_zonal ] ) > max_u ) 
				{
					max_u = fabs ( u.x[ i ][ j ][ k_zonal ] );
					if ( max_u == 0. ) max_u = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] / max_u << endl;
				Atmosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing v-component
		Atmosphere_vtk_zonal_File <<  "SCALARS v-Component float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( v.x[ i ][ j ][ k_zonal ] ) > max_v ) 
				{
					max_v = fabs ( v.x[ i ][ j ][ k_zonal ] );
					if ( max_v == 0. ) max_v = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << v.x[ i ][ j ][ k_zonal ] / max_v << endl;
				Atmosphere_vtk_zonal_File << v.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing w-component
		Atmosphere_vtk_zonal_File <<  "SCALARS w-Component float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" <<endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( w.x[ i ][ j ][ k_zonal ] ) > max_w ) 
				{
					max_w = fabs ( w.x[ i ][ j ][ k_zonal ] );
					if ( max_w == 0. ) max_w = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << w.x[ i ][ j ][ k_zonal ] / max_w << endl;
				Atmosphere_vtk_zonal_File << w.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing of temperature
		Atmosphere_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( t.x[ i ][ j ][ k_zonal ] ) > max_t ) 
				{
					max_t = fabs ( t.x[ i ][ j ][ k_zonal ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << t.x[ i ][ j ][ k_zonal ] / max_t << endl;
				Atmosphere_vtk_zonal_File << t.x[ i ][ j ][ k_zonal ] * t_0 - t_0 << endl;
			}
		}


// writing of condensation temperature
		Atmosphere_vtk_zonal_File <<  "SCALARS CondensationTemp float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( t_cond_3D.x[ i ][ j ][ k_zonal ] ) > max_t ) 
				{
					max_t = fabs ( t_cond_3D.x[ i ][ j ][ k_zonal ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << t_cond_3D.x[ i ][ j ][ k_zonal ] / max_t << endl;
				Atmosphere_vtk_zonal_File << t_cond_3D.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing of evaporation temperature
		Atmosphere_vtk_zonal_File <<  "SCALARS EvaporationTemp float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( t_evap_3D.x[ i ][ j ][ k_zonal ] ) > max_t ) 
				{
					max_t = fabs ( t_evap_3D.x[ i ][ j ][ k_zonal ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << t_evap_3D.x[ i ][ j ][ k_zonal ] / max_t << endl;
				Atmosphere_vtk_zonal_File << t_evap_3D.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing of epsilon_3D
		Atmosphere_vtk_zonal_File <<  "SCALARS Epsilon_3D float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( epsilon_3D.x[ i ][ j ][ k_zonal ] ) > max_t ) 
				{
					max_t = fabs ( epsilon_3D.x[ i ][ j ][ k_zonal ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << epsilon_3D.x[ i ][ j ][ k_zonal ] / max_t << endl;
				Atmosphere_vtk_zonal_File << epsilon_3D.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing water vapour
		Atmosphere_vtk_zonal_File <<  "SCALARS WaterVapour float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( c.x[ i ][ j ][ k_zonal ] ) > max_c ) 
				{
					max_c = fabs ( c.x[ i ][ j ][ k_zonal ] );
					if ( max_c == 0. ) max_c = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << c.x[ i ][ j ][ k_zonal ] / max_c << endl;
				Atmosphere_vtk_zonal_File << c.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing cloud water
		Atmosphere_vtk_zonal_File <<  "SCALARS CloudWater float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( cloud.x[ i ][ j ][ k_zonal ] ) > max_cloud ) 
				{
					max_cloud = fabs ( cloud.x[ i ][ j ][ k_zonal ] );
					if ( max_cloud == 0. ) max_cloud = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << cloud.x[ i ][ j ][ k_zonal ] / max_cloud << endl;
				Atmosphere_vtk_zonal_File << cloud.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing cloud ice
		Atmosphere_vtk_zonal_File <<  "SCALARS CloudIce float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( ice.x[ i ][ j ][ k_zonal ] ) > max_ice ) 
				{
					max_ice = fabs ( ice.x[ i ][ j ][ k_zonal ] );
					if ( max_ice == 0. ) max_ice = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << ice.x[ i ][ j ][ k_zonal ] / max_ice << endl;
				Atmosphere_vtk_zonal_File << ice.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing P_rain
		Atmosphere_vtk_zonal_File <<  "SCALARS PrecipitationRain float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( P_rain.x[ i ][ j ][ k_zonal ] ) > max_P_rain ) 
				{
					max_P_rain = fabs ( P_rain.x[ i ][ j ][ k_zonal ] );
					if ( max_P_rain == 0. ) max_P_rain = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << P_rain.x[ i ][ j ][ k_zonal ] / max_P_rain << endl;
				Atmosphere_vtk_zonal_File << P_rain.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing P_snow
		Atmosphere_vtk_zonal_File <<  "SCALARS PrecipitationSnow float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				max_P_snow = fabs ( P_snow.x[ i ][ j ][ k_zonal ] );
				if ( max_P_snow == 0. ) max_P_snow = 1.e-6;
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				if ( t.x[ i ][ j ][ k_zonal ] * t_0 - t_0 >= 0. )		P_snow.x[ i ][ j ][ k_zonal ] = 0.;
//				Atmosphere_vtk_zonal_File << P_snow.x[ i ][ j ][ k_zonal ] / max_P_snow << endl;
				Atmosphere_vtk_zonal_File << P_snow.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing S_v
		Atmosphere_vtk_zonal_File <<  "SCALARS Source_WaterVapour float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Atmosphere_vtk_zonal_File << S_v.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing S_c
		Atmosphere_vtk_zonal_File <<  "SCALARS Source_CloudWater float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Atmosphere_vtk_zonal_File << S_c.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing S_i
		Atmosphere_vtk_zonal_File <<  "SCALARS Source_CloudIce float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Atmosphere_vtk_zonal_File << S_i.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing S_r
		Atmosphere_vtk_zonal_File <<  "SCALARS Source_Rain float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Atmosphere_vtk_zonal_File << S_r.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


// writing S_s
		Atmosphere_vtk_zonal_File <<  "SCALARS Source_Snow float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Atmosphere_vtk_zonal_File << S_s.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}


/*
// writing M_u
		Atmosphere_vtk_zonal_File <<  "SCALARS Updraft float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				max_M_u = fabs ( M_u.x[ i ][ j ][ k_zonal ] );
				if ( max_M_u == 0. ) max_M_u = 1.e-6;
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				if ( t.x[ i ][ j ][ k_zonal ] * t_0 - t_0 >= 0. )		M_u.x[ i ][ j ][ k_zonal ] = 0.;
//				Atmosphere_vtk_zonal_File << M_u.x[ i ][ j ][ k_zonal ] / max_M_u << endl;
				Atmosphere_vtk_zonal_File << M_u.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing M_d
		Atmosphere_vtk_zonal_File <<  "SCALARS Downdraft float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				max_M_d = fabs ( M_d.x[ i ][ j ][ k_zonal ] );
				if ( max_M_d == 0. ) max_M_d = 1.e-6;
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				if ( t.x[ i ][ j ][ k_zonal ] * t_0 - t_0 >= 0. )		M_d.x[ i ][ j ][ k_zonal ] = 0.;
//				Atmosphere_vtk_zonal_File << M_d.x[ i ][ j ][ k_zonal ] / max_M_d << endl;
				Atmosphere_vtk_zonal_File << M_d.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing P_co2_nv
		Atmosphere_vtk_zonal_File <<  "SCALARS PrecipitationConv float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				max_P_conv = fabs ( P_co2_nv.x[ i ][ j ][ k_zonal ] );
				if ( max_P_conv == 0. ) max_P_conv = 1.e-6;
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				if ( t.x[ i ][ j ][ k_zonal ] * t_0 - t_0 >= 0. )		P_co2_nv.x[ i ][ j ][ k_zonal ] = 0.;
//				Atmosphere_vtk_zonal_File << P_co2_nv.x[ i ][ j ][ k_zonal ] / max_P_conv << endl;
				Atmosphere_vtk_zonal_File << P_co2_nv.x[ i ][ j ][ k_zonal ] * 1000. << endl;
			}
		}
*/

// writing CO2
		Atmosphere_vtk_zonal_File <<  "SCALARS CO2-Concentration float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( co2.x[ i ][ j ][ k_zonal ] ) > max_co2 ) 
				{
					max_co2 = fabs ( co2.x[ i ][ j ][ k_zonal ] );
					if ( max_co2 == 0. ) max_co2 = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << co2.x[ i ][ j ][ k_zonal ] / max_co2 << endl;
				Atmosphere_vtk_zonal_File << co2.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing dynamic pressure
		Atmosphere_vtk_zonal_File <<  "SCALARS PressureDynamic float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

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
//				Atmosphere_vtk_zonal_File << p_dyn.x[ i ][ j ][ k_zonal ] / max_p_dyn << endl;
				Atmosphere_vtk_zonal_File << p_dyn.x[ i ][ j ][ k_zonal ] * u_0 * u_0 * r_air << endl;
			}
		}


// writing static pressure
		Atmosphere_vtk_zonal_File <<  "SCALARS PressureStatic float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

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
//				Atmosphere_vtk_zonal_File << p_stat.x[ i ][ j ][ k_zonal ] / max_p_stat << endl;
				Atmosphere_vtk_zonal_File << p_stat.x[ i ][ j ][ k_zonal ] << endl;
			}
		}



// writing buoyancy force
		Atmosphere_vtk_zonal_File <<  "SCALARS BuoyancyForce float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( BuoyancyForce.x[ i ][ j ][ k_zonal ] ) > max_buoyancy_force ) 
				{
					max_buoyancy_force = fabs ( BuoyancyForce.x[ i ][ j ][ k_zonal ] );
					if ( max_buoyancy_force == 0. ) max_buoyancy_force = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << BuoyancyForce.x[ i ][ j ][ k_zonal ] / max_buoyancy_force << endl;
				Atmosphere_vtk_zonal_File << BuoyancyForce.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing topography
		Atmosphere_vtk_zonal_File <<  "SCALARS Topography float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Atmosphere_vtk_zonal_File << h.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing Latency
		Atmosphere_vtk_zonal_File <<  "SCALARS Latency float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( Latency.x[ i ][ j ][ k_zonal ] ) > max_Latency ) 
				{
					max_Latency = fabs ( Latency.x[ i ][ j ][ k_zonal ] );
					if ( max_Latency == 0. ) max_Latency = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << Latency.x[ i ][ j ][ k_zonal ] / max_Latency << endl;
				Atmosphere_vtk_zonal_File << Latency.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing Q_Sensible
		Atmosphere_vtk_zonal_File <<  "SCALARS Q_Sensible float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( Q_Sensible.x[ i ][ j ][ k_zonal ] ) > max_Q_Sensible ) 
				{
					max_Q_Sensible = fabs ( Q_Sensible.x[ i ][ j ][ k_zonal ] );
					if ( max_Q_Sensible == 0. ) max_Q_Sensible = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << Q_Sensible.x[ i ][ j ][ k_zonal ] / max_Q_Sensible << endl;
				Atmosphere_vtk_zonal_File << Q_Sensible.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


/*
// writing t_cond_3D
		Atmosphere_vtk_zonal_File <<  "SCALARS t_cond_3D float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( t_cond_3D.x[ i ][ j ][ k_zonal ] ) > max_t_t_Condensation_3D ) 
				{
					max_t_t_Condensation_3D = fabs ( t_cond_3D.x[ i ][ j ][ k_zonal ] );
					if ( max_t_t_Condensation_3D == 0. ) max_t_t_Condensation_3D = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << t_cond_3D.x[ i ][ j ][ k_zonal ] / max_t_t_Condensation_3D << endl;
				Atmosphere_vtk_zonal_File << t_cond_3D.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing t_evap_3D
		Atmosphere_vtk_zonal_File <<  "SCALARS t_evap_3D float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( t_evap_3D.x[ i ][ j ][ k_zonal ] ) > max_t_t_Evaporation_3D ) 
				{
					max_t_t_Evaporation_3D = fabs ( t_evap_3D.x[ i ][ j ][ k_zonal ] );
					if ( max_t_t_Evaporation_3D == 0. ) max_t_t_Evaporation_3D = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << t_evap_3D.x[ i ][ j ][ k_zonal ] / max_t_t_Evaporation_3D << endl;
				Atmosphere_vtk_zonal_File << t_evap_3D.x[ i ][ j ][ k_zonal ] << endl;
			}
		}
*/



// writing radiation_3D
		Atmosphere_vtk_zonal_File <<  "SCALARS Radiation_3D float " << 1 << endl;
		Atmosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( radiation_3D.x[ i ][ j ][ k_zonal ] ) > max_radiation_3D ) 
				{
					max_radiation_3D = fabs ( radiation_3D.x[ i ][ j ][ k_zonal ] );
					if ( max_radiation_3D == 0. ) max_radiation_3D = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << radiation_3D.x[ i ][ j ][ k_zonal ] / max_radiation_3D << endl;
				Atmosphere_vtk_zonal_File << radiation_3D.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


// writing zonal u-v cell structure
		Atmosphere_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int j = 1; j < jm-1; j++ )
			{
				if ( fabs ( u.x[ i ][ j ][ k_zonal ] ) > max_u ) 
				{
					max_u = fabs ( u.x[ i ][ j ][ k_zonal ] );
					if ( max_u == 0. ) max_u = 1.e-6;
				}
				if ( fabs ( v.x[ i ][ j ][ k_zonal ] ) > max_v ) 
				{
					max_v = fabs ( v.x[ i ][ j ][ k_zonal ] );
					if ( max_v == 0. ) max_v = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
//				Atmosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] / max_u << " " << v.x[ i ][ j ][ k_zonal ] / max_v << " " << z << endl;
				Atmosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] << " " << v.x[ i ][ j ][ k_zonal ] << " " << z << endl;
			}
		}


// end writing



// final file administration
		endpos = Atmosphere_vtk_zonal_File.tellp();
		cout << "***** file ::::: " << Atmosphere_zonal_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Atmosphere_zonal_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
	}


// in case of failing

	if ( Atmosphere_vtk_zonal_File.fail() )
	{
		cout << "***** file ::::: " << Atmosphere_zonal_File_Name.str() << " ::::: is lost!" << endl;
		return;
	}

	Atmosphere_vtk_zonal_File.close();

	if ( Atmosphere_vtk_zonal_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_zonal_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Atmosphere_vtk_zonal_File.fail() )
		cout << "***** file ::::: " << Atmosphere_zonal_File_Name.str() << " ::::: could not be closed properly!" << endl;

return;
}





void PostProcess_Atmosphere::paraview_vtk_longal (string &Name_Bathymetry_File, int &j_longal, int &pressure_iter, double &u_0, double &t_0, double &p_0, double &r_air, double &c_0, double &co2_0, double &radiation_equator, Array &h, Array &p_dyn, Array &p_stat, Array &t_cond_3D, Array &t_evap_3D, Array &BuoyancyForce, Array &t, Array &u, Array &v, Array &w, Array &c, Array &co2, Array &cloud, Array &ice, Array &aux_u, Array &aux_v, Array &aux_w, Array &Latency, Array &Q_Sensible, Array &IceLayer, Array &epsilon_3D, Array &P_rain, Array &P_snow )
{
	double x, y, z, dx, dz;

	stringstream Atmosphere_longal_File_Name;

	max_u = max_v = max_w = max_t = max_c = max_co2 = max_cloud = max_ice = max_P_rain = max_P_snow = max_P_conv = max_M_u = max_M_d = max_p_dyn = max_p_stat = 0.;
	max_Rain = max_Rain_super = max_Ice = max_Latency = max_Q_Sensible = max_Precipitation = 0.;
	max_t_Evaporation = max_t_Condensation = max_t_evap_3D = max_t_cond_3D = 0.;
	max_precipitable_water = max_IceAir = max_Q_bottom = max_Q_latent = max_Q_sensible = 0.;
	max_t_Evaporation_Penman = max_t_Evaporation_Haude = max_Q_Radiation = max_buoyancy_force = 0.;
	max_Q_t_Evaporation = max_precipitation_NASA = max_Water = max_Water_super = max_Vegetation = max_IceLayer = 0.;

	i_max = im;
	k_max = km;

// file administration
	streampos anfangpos, endpos;

	Atmosphere_longal_File_Name <<  "[" << Name_Bathymetry_File << "]_Atm_longal_" << j_longal << "_" << pressure_iter << ".vtk";
	ofstream Atmosphere_vtk_longal_File;
	Atmosphere_vtk_longal_File.precision ( 4 );
	Atmosphere_vtk_longal_File.setf ( ios::fixed );
	Atmosphere_vtk_longal_File.open ( Atmosphere_longal_File_Name.str().c_str(), ios_base::out );

	if ( Atmosphere_vtk_longal_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_longal_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Atmosphere_vtk_longal_File.tellp();
		cout << "***** file ::::: " << Atmosphere_longal_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;


// begin writing
		Atmosphere_vtk_longal_File <<  "# vtk DataFile Version 3.0" << endl;
		Atmosphere_vtk_longal_File <<  "Longitudinal_Data_Atmosphere_Circulation\n";
		Atmosphere_vtk_longal_File <<  "ASCII" << endl;
		Atmosphere_vtk_longal_File <<  "DATASET STRUCTURED_GRID" << endl;
		Atmosphere_vtk_longal_File <<  "DIMENSIONS " << k_max << " "<< i_max << " " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "POINTS " << i_max * k_max << " float" << endl;


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

				Atmosphere_vtk_longal_File << x << " " << y << " "<< z << endl;
			}
			z = 0.;
			x = x + dx;
		}



		Atmosphere_vtk_longal_File <<  "POINT_DATA " << i_max * k_max << endl;

// writing u-component
		Atmosphere_vtk_longal_File <<  "SCALARS u-Component float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( u.x[ i ][ j_longal ][ k ] ) > max_u ) 
				{
					max_u = fabs ( u.x[ i ][ j_longal ][ k ] );
					if ( max_u == 0. ) max_u = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] / max_u << endl;
				Atmosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing v-component
		Atmosphere_vtk_longal_File <<  "SCALARS v-Component float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( v.x[ i ][ j_longal ][ k ] ) > max_v ) 
				{
					max_v = fabs ( v.x[ i ][ j_longal ][ k ] );
					if ( max_v == 0. ) max_v = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << v.x[ i ][ j_longal ][ k ] / max_v << endl;
				Atmosphere_vtk_longal_File << v.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing w-component
		Atmosphere_vtk_longal_File <<  "SCALARS w-Component float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" <<endl;

		for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( w.x[ i ][ j_longal ][ k ] ) > max_w ) 
				{
					max_w = fabs ( w.x[ i ][ j_longal ][ k ] );
					if ( max_w == 0. ) max_w = 1.e-6;
				}
			}
		}


		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << w.x[ i ][ j_longal ][ k ] / max_w << endl;
				Atmosphere_vtk_longal_File << w.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing of temperature
		Atmosphere_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( t.x[ i ][ j_longal ][ k ] ) > max_t ) 
				{
					max_t = fabs ( t.x[ i ][ j_longal ][ k ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << t.x[ i ][ j_longal ][ k ] / max_t << endl;
				Atmosphere_vtk_longal_File << t.x[ i ][ j_longal ][ k ] * t_0 - t_0 << endl;
			}
		}


// writing of condensation temperature
		Atmosphere_vtk_longal_File <<  "SCALARS CondensationTemp float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( t_cond_3D.x[ i ][ j_longal ][ k ] ) > max_t ) 
				{
					max_t = fabs ( t_cond_3D.x[ i ][ j_longal ][ k ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << t_cond_3D.x[ i ][ j_longal ][ k ] / max_t << endl;
				Atmosphere_vtk_longal_File << t_cond_3D.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing of evaporation temperature
		Atmosphere_vtk_longal_File <<  "SCALARS EvaporationTemp float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( t_evap_3D.x[ i ][ j_longal ][ k ] ) > max_t ) 
				{
					max_t = fabs ( t_evap_3D.x[ i ][ j_longal ][ k ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << t_evap_3D.x[ i ][ j_longal ][ k ] / max_t << endl;
				Atmosphere_vtk_longal_File << t_evap_3D.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing of epsilon_3D
		Atmosphere_vtk_longal_File <<  "SCALARS Epsilon_3D float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( epsilon_3D.x[ i ][ j_longal ][ k ] ) > max_t ) 
				{
					max_t = fabs ( epsilon_3D.x[ i ][ j_longal ][ k ] );
					if ( max_t == 0. ) max_t = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << epsilon_3D.x[ i ][ j_longal ][ k ] / max_t << endl;
				Atmosphere_vtk_longal_File << epsilon_3D.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing water vapour
		Atmosphere_vtk_longal_File <<  "SCALARS WaterVapour float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( c.x[ i ][ j_longal ][ k ] ) > max_c ) 
				{
					max_c = fabs ( c.x[ i ][ j_longal ][ k ] );
					if ( max_c == 0. ) max_c = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << c.x[ i ][ j_longal ][ k ] / max_c << endl;
				Atmosphere_vtk_longal_File << c.x[ i ][ j_longal ][ k ] * 1000. << endl;
			}
		}


// writing cloud water
		Atmosphere_vtk_longal_File <<  "SCALARS CloudWater float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( cloud.x[ i ][ j_longal ][ k ] ) > max_cloud ) 
				{
					max_cloud = fabs ( cloud.x[ i ][ j_longal ][ k ] );
					if ( max_cloud == 0. ) max_cloud = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << cloud.x[ i ][ j_longal ][ k ] / max_cloud << endl;
				Atmosphere_vtk_longal_File << cloud.x[ i ][ j_longal ][ k ] * 1000. << endl;
			}
		}


// writing cloud ice
		Atmosphere_vtk_longal_File <<  "SCALARS CloudIce float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( ice.x[ i ][ j_longal ][ k ] ) > max_ice ) 
				{
					max_ice = fabs ( ice.x[ i ][ j_longal ][ k ] );
					if ( max_ice == 0. ) max_ice = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << ice.x[ i ][ j_longal ][ k ] / max_ice << endl;
				Atmosphere_vtk_longal_File << ice.x[ i ][ j_longal ][ k ] * 1000. << endl;
			}
		}


// writing P_rain
		Atmosphere_vtk_longal_File <<  "SCALARS PrecipitationRain float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( P_rain.x[ i ][ j_longal ][ k ] ) > max_P_rain ) 
				{
					max_P_rain = fabs ( P_rain.x[ i ][ j_longal ][ k ] );
					if ( max_P_rain == 0. ) max_P_rain = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << P_rain.x[ i ][ j_longal ][ k ] / max_P_rain << endl;
				Atmosphere_vtk_longal_File << P_rain.x[ i ][ j_longal ][ k ] * 1000. << endl;
			}
		}


// writing P_snow
		Atmosphere_vtk_longal_File <<  "SCALARS PrecipitationSnow float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				max_P_snow = fabs ( P_snow.x[ i ][ j_longal ][ k ] );
				if ( max_P_snow == 0. ) max_P_snow = 1.e-6;
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				if ( t.x[ i ][ j_longal ][ k ] * t_0 - t_0 >= 0. )		P_snow.x[ i ][ j_longal ][ k ] = 0.;
//				Atmosphere_vtk_longal_File << P_snow.x[ i ][ j_longal ][ k ] / max_P_snow << endl;
				Atmosphere_vtk_longal_File << P_snow.x[ i ][ j_longal ][ k ] * 1000. << endl;
			}
		}

/*
// writing M_u
		Atmosphere_vtk_longal_File <<  "SCALARS Updraft float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				max_M_u = fabs ( M_u.x[ i ][ j_longal ][ k ] );
				if ( max_M_u == 0. ) max_M_u = 1.e-6;
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				if ( t.x[ i ][ j_longal ][ k ] * t_0 - t_0 >= 0. )		M_u.x[ i ][ j_longal ][ k ] = 0.;
//				Atmosphere_vtk_longal_File << M_u.x[ i ][ j_longal ][ k ] / max_M_u << endl;
				Atmosphere_vtk_longal_File << M_u.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing M_d
		Atmosphere_vtk_longal_File <<  "SCALARS Downdraft float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				max_M_d = fabs ( M_d.x[ i ][ j_longal ][ k ] );
				if ( max_M_d == 0. ) max_M_d = 1.e-6;
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				if ( t.x[ i ][ j_longal ][ k ] * t_0 - t_0 >= 0. )		M_d.x[ i ][ j_longal ][ k ] = 0.;
//				Atmosphere_vtk_longal_File << M_d.x[ i ][ j_longal ][ k ] / max_M_d << endl;
				Atmosphere_vtk_longal_File << M_d.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing P_co2_nv
		Atmosphere_vtk_longal_File <<  "SCALARS PrecipitationConv float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				max_P_conv = fabs ( P_co2_nv.x[ i ][ j_longal ][ k ] );
				if ( max_P_conv == 0. ) max_P_conv = 1.e-6;
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				if ( t.x[ i ][ j_longal ][ k ] * t_0 - t_0 >= 0. )		P_co2_nv.x[ i ][ j_longal ][ k ] = 0.;
//				Atmosphere_vtk_longal_File << P_co2_nv.x[ i ][ j_longal ][ k ] / max_P_conv << endl;
				Atmosphere_vtk_longal_File << P_co2_nv.x[ i ][ j_longal ][ k ] * 1000. << endl;
			}
		}
*/

// writing CO2
		Atmosphere_vtk_longal_File <<  "SCALARS CO2-Concentration float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( co2.x[ i ][ j_longal ][ k ] ) > max_co2 ) 
				{
					max_co2 = fabs ( co2.x[ i ][ j_longal ][ k ] );
					if ( max_co2 == 0. ) max_co2 = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << co2.x[ i ][ j_longal ][ k ] / max_co2 << endl;
				Atmosphere_vtk_longal_File << co2.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing dynamic pressure
		Atmosphere_vtk_longal_File <<  "SCALARS PressureDynamic float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( p_dyn.x[ i ][ j_longal ][ k ] ) > max_p_dyn ) 
				{
					max_p_dyn = fabs ( p_dyn.x[ i ][ j_longal ][ k ] );
					if ( max_p_dyn == 0. ) max_p_dyn = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << p_dyn.x[ i ][ j_longal ][ k ] / max_p_dyn << endl;
				Atmosphere_vtk_longal_File << p_dyn.x[ i ][ j_longal ][ k ] * u_0 * u_0 * r_air << endl;
			}
		}

/*
// writing static pressure
		Atmosphere_vtk_longal_File <<  "SCALARS PressureStatic float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

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
//				Atmosphere_vtk_longal_File << p_stat.x[ i ][ j_longal ][ k ] / max_p_stat << endl;
				Atmosphere_vtk_longal_File << p_stat.x[ i ][ j_longal ][ k ] << endl;
			}
		}
*/

// writing buoyancy force
		Atmosphere_vtk_longal_File <<  "SCALARS BuoyancyForce float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( BuoyancyForce.x[ i ][ j_longal ][ k ] ) > max_buoyancy_force ) 
				{
					max_buoyancy_force = fabs ( BuoyancyForce.x[ i ][ j_longal ][ k ] );
					if ( max_buoyancy_force == 0. ) max_buoyancy_force = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << BuoyancyForce[ i ][ j_longal ][ k ] / max_buoyancy_force << endl;
				Atmosphere_vtk_longal_File << BuoyancyForce.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing topography
		Atmosphere_vtk_longal_File <<  "SCALARS Topography float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Atmosphere_vtk_longal_File << h.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing Latency
		Atmosphere_vtk_longal_File <<  "SCALARS Latency float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Latency.x[ i ][ j_longal ][ k ] ) > max_Latency ) 
				{
					max_Latency = fabs ( Latency.x[ i ][ j_longal ][ k ] );
					if ( max_Latency == 0. ) max_Latency = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << Latency.x[ i ][ j_longal ][ k ] / max_Latency << endl;
				Atmosphere_vtk_longal_File << Latency.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing Q_Sensible
		Atmosphere_vtk_longal_File <<  "SCALARS Q_Sensible float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( Q_Sensible.x[ i ][ j_longal ][ k ] ) > max_Q_Sensible ) 
				{
					max_Q_Sensible = fabs ( Q_Sensible.x[ i ][ j_longal ][ k ] );
					if ( max_Q_Sensible == 0. ) max_Q_Sensible = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << Q_Sensible.x[ i ][ j_longal ][ k ] / max_Q_Sensible << endl;
				Atmosphere_vtk_longal_File << Q_Sensible.x[ i ][ j_longal ][ k ] << endl;
			}
		}


/*
// writing t_cond_3D
		Atmosphere_vtk_longal_File <<  "SCALARS t_cond_3D float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( t_cond_3D.x[ i ][ j_longal ][ k ] ) > max_t_t_Condensation_3D ) 
				{
					max_t_t_Condensation_3D = fabs ( t_cond_3D.x[ i ][ j_longal ][ k ] );
					if ( max_t_t_Condensation_3D == 0. ) max_t_t_Condensation_3D = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << t_cond_3D.x[ i ][ j_longal ][ k ] / max_t_t_Condensation_3D << endl;
				Atmosphere_vtk_longal_File << t_cond_3D.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// writing t_evap_3D
		Atmosphere_vtk_longal_File <<  "SCALARS t_evap_3D float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( t_evap_3D.x[ i ][ j_longal ][ k ] ) > max_t_t_Evaporation_3D ) 
				{
					max_t_t_Evaporation_3D = fabs ( t_evap_3D.x[ i ][ j_longal ][ k ] );
					if ( max_t_t_Evaporation_3D == 0. ) max_t_t_Evaporation_3D = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
//				Atmosphere_vtk_longal_File << t_evap_3D.x[ i ][ j_longal ][ k ] / max_t_t_Evaporation_3D << endl;
				Atmosphere_vtk_longal_File << t_evap_3D.x[ i ][ j_longal ][ k ] << endl;
			}
		}
*/



/*
// writing ice layer
		Atmosphere_vtk_longal_File <<  "SCALARS IceLayer float " << 1 << endl;
		Atmosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( IceLayer.x[ i ][ j_longal ][ k ] ) > max_IceLayer ) 
				{
					max_IceLayer = fabs ( IceLayer.x[ i ][ j_longal ][ k ] );
					if ( max_IceLayer == 0. ) max_IceLayer = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Atmosphere_vtk_longal_File << IceLayer.x[ i ][ j_longal ][ k ] / max_IceLayer<< endl;
			}
		}
*/



// writing longitudinal u-v cell structure
		Atmosphere_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;

	for ( int i = 1; i < im-1; i++ )
		{
			for ( int k = 1; k < km-1; k++ )
			{
				if ( fabs ( u.x[ i ][ j_longal ][ k ] ) > max_u ) 
				{
					max_u = fabs ( u.x[ i ][ j_longal ][ k ] );
					if ( max_u == 0. ) max_u = 1.e-6;
				}
				if ( fabs ( w.x[ i ][ j_longal ][ k ] ) > max_w ) 
				{
					max_w = fabs ( w.x[ i ][ j_longal ][ k ] );
					if ( max_w == 0. ) max_w = 1.e-6;
				}
			}
		}

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
//				Atmosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] / max_u << " " << y << " " << w.x[ i ][ j_longal ][ k ] / max_w << endl;
				Atmosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] << " " << y << " " << w.x[ i ][ j_longal ][ k ] << endl;
			}
		}


// end writing



// final file administration

		endpos = Atmosphere_vtk_longal_File.tellp();
		cout << "***** file ::::: " << Atmosphere_longal_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Atmosphere_longal_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
	}


// in case of failing

	if ( Atmosphere_vtk_longal_File.fail() )
	{
		cout << "***** file ::::: " << Atmosphere_longal_File_Name.str() << " ::::: is lost!" << endl;
		return;
	}

	Atmosphere_vtk_longal_File.close();

	if ( Atmosphere_vtk_longal_File.good() )
	{
		cout << "***** file ::::: " << Atmosphere_longal_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Atmosphere_vtk_longal_File.fail() )
	  cout << "***** file ::::: " << Atmosphere_longal_File_Name.str() << " ::::: could not be closed properly!" << endl;

return;
}






void PostProcess_Atmosphere::Atmosphere_PlotData ( string &Name_Bathymetry_File, double u_0, double t_0, Array &v, Array &w, Array &t, Array &c, Array_2D &Precipitation, Array_2D &precipitable_water )
{
	stringstream Name_PlotData_File;

	streampos anfangpos, endpos;

// file administration
	Name_PlotData_File << "[" << Name_Bathymetry_File << "]_PlotData_Atm.xyz";
	ofstream PlotData_File;
	PlotData_File.precision ( 4 );
	PlotData_File.setf ( ios::fixed );
	PlotData_File.open ( Name_PlotData_File.str().c_str(), ios_base::out );

	if ( PlotData_File.good() )
	{
		cout << "***** file ::::: " << Name_PlotData_File.str() << " ::::: could be opened" << endl;
		cout << endl;
		anfangpos = PlotData_File.tellp();
		cout << "***** file ::::: " << Name_PlotData_File.str() << " ::::: starts at ::::::: " << anfangpos << endl;


// begin writing

		cout << "\n\n”***** Atmosphere_PlotData_File_write:   begin of writing!" << endl << endl;


		PlotData_File << " latitude ( ° )" << "  , " << "longitude ( ° )" << "  ,    " << "v-velocity ( m/s )" << "   ,   " << "w-velocity ( m/s )" << "   ,   " << "temperature ( °C )" << "   ,  " << "water_vapour ( g/kg )" << "   ,   " << "precipitation ( mm )" << "   ,   " <<  "precipitable water ( mm )" << endl;

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				PlotData_File << k << " " << j << " " << v.x[ 0 ][ j ][ k ] * u_0 << " " << w.x[ 0 ][ j ][ k ] * u_0 << " " << t.x[ 0 ][ j ][ k ] * t_0 - t_0 << " " << c.x[ 0 ][ j ][ k ] * 1000. << " " << Precipitation.y[ j ][ k ] << " " << precipitable_water.y[ j ][ k ] << " " <<  endl;
			}
		}



// end writing


// final file administration

		endpos = PlotData_File.tellp();
		cout << "***** file ::::: " << Name_PlotData_File.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Name_PlotData_File.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
		cout << endl;
	}


// in case of failing

	if ( PlotData_File.fail() )
	{
		cout << "***** file ::::: " << Name_PlotData_File.str() << " ::::: is lost!" << endl;
		return;
	}

	PlotData_File.close();

	if ( PlotData_File.good() )
	{
		cout << "***** file ::::: " << Name_PlotData_File.str() << " ::::: could be closed after writing" << endl;
		cout << endl;
	}

	if ( PlotData_File.fail() )
		cout << "***** file ::::: " << Name_PlotData_File.str() << " ::::: could not be closed properly!" << endl << endl << endl;


return;
}






