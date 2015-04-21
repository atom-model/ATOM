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

#include "PostProcess_Hydrosphere.h"

using namespace std;



PostProcess_Hydrosphere::PostProcess_Hydrosphere ( int im, int jm, int km )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
}


PostProcess_Hydrosphere::~PostProcess_Hydrosphere() {}




void PostProcess_Hydrosphere::paraview_vts ( const string &Name_Bathymetry_File, int &n, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &p, Array &u, Array &v, Array &w, Array &c, Array &fup, Array &fvp, Array &fwp, Array &fcp, Array &fpp, Array &ftp, Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion, Array &Salt_Balance )
{
	double x, y, z, sinthe, sinphi, costhe, cosphi;

	stringstream Hydrosphere_vts_File_Name;

	streampos anfangpos, endpos;

// file administration
	Hydrosphere_vts_File_Name << "[" << Name_Bathymetry_File << "]_Hyd_Kreide_" << n << ".vts";
	ofstream Hydrosphere_vts_File;
	Hydrosphere_vts_File.precision ( 4 );
	Hydrosphere_vts_File.setf ( ios::fixed );
	Hydrosphere_vts_File.open ( Hydrosphere_vts_File_Name.str().c_str(), ios_base::out );

	if ( Hydrosphere_vts_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_vts_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Hydrosphere_vts_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_vts_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;


// file administration
// begin writing

		Hydrosphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
		Hydrosphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
		Hydrosphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
		Hydrosphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

//		Hydrosphere_vts_File <<  "   <PointData Vectors=\"Velocity Rotation\" Scalars=\"SeaGround Temperature Pressure SaltConcentration u-Component v-Component w-Component\">\n"  << endl;
		Hydrosphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"SeaGround Temperature Pressure SaltConcentration\">\n"  << endl;



// writing u, v und w velocity components in cartesian coordinates

		Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;


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

		Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SeaGround\" format=\"ascii\">\n"  << endl;

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

		Hydrosphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SaltConcentration\" format=\"ascii\">\n"  << endl;

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


// final file administration
		endpos = Hydrosphere_vts_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_vts_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Hydrosphere_vts_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes"<< endl;
	}


// in case of failing

	if ( Hydrosphere_vts_File.fail() )
	{
		cout << "***** file ::::: " << Hydrosphere_vts_File_Name.str() << " ::::: is lost" << endl;
		return;
	}

	Hydrosphere_vts_File.close();

	if ( Hydrosphere_vts_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_vts_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Hydrosphere_vts_File.fail() )
		cout << "***** file ::::: " << Hydrosphere_vts_File_Name.str() << " ::::: could not be closed properly" << endl;
return;
}





void PostProcess_Hydrosphere::paraview_panorama_vts ( const string &Name_Bathymetry_File, int &pressure_iter, Array &h, Array &t, Array &p, Array &u, Array &v, Array &w, Array &c, Array &aux_u, Array &aux_v, Array &aux_w, Array &Salt_Finger, Array &Salt_Diffusion, Array &Salt_Balance )
{
	double x, y, z, dx, dy, dz;

	stringstream Hydrosphere_panorama_vts_File_Name;


// file administration
	streampos anfangpos, endpos;

//	sprintf ( Hydrosphere_panorama_vts_File_Name, "[%s]_Hyd_panorama_%i.vts", Name_Bathymetry_File, n );
	Hydrosphere_panorama_vts_File_Name <<  "[" << Name_Bathymetry_File << "]_Hyd_panorama_" << pressure_iter << ".vts";
	ofstream Hydrosphere_panorama_vts_File;
	Hydrosphere_panorama_vts_File.precision ( 4 );
	Hydrosphere_panorama_vts_File.setf ( ios::fixed );
	Hydrosphere_panorama_vts_File.open ( Hydrosphere_panorama_vts_File_Name.str().c_str(), ios_base::out );

	if ( Hydrosphere_panorama_vts_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_panorama_vts_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Hydrosphere_panorama_vts_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_panorama_vts_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;



// begin writing

		Hydrosphere_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
		Hydrosphere_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
		Hydrosphere_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
		Hydrosphere_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;

//		Hydrosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"SeaGround Temperature Pressure SaltConcentration SaltFinger.y SaltDiffusion.y SaltBalance\">\n"  << endl;
		Hydrosphere_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"SeaGround Temperature Pressure SaltConcentration\">\n"  << endl;



// writing u, v und w velocity components in cartesian coordinates

		Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;




		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
// transformtion from spherical to cartesian coordinates for representation in ParaView
					Hydrosphere_panorama_vts_File << u.x[ i ][ j ][ k ] << " " << v.x[ i ][ j ][ k ] << " " << w.x[ i ][ j ][ k ]  << endl;
				}
				Hydrosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;




// writing of sea ground

		Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SeaGround\" format=\"ascii\">\n"  << endl;

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



// writing temperature

		Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">\n"  << endl;

		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					Hydrosphere_panorama_vts_File << p.x[ i ][ j ][ k ] << endl;
				}
				Hydrosphere_panorama_vts_File <<  "\n"  << endl;
			}
			Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		}
		Hydrosphere_panorama_vts_File <<  "\n"  << endl;
		Hydrosphere_panorama_vts_File <<  "    </DataArray>\n" << endl;


// writing scalar function c

		Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SaltConcentration\" format=\"ascii\">\n"  << endl;

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

/*
		Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SaltFinger.y\" format=\"ascii\">\n"  << endl;

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


// writing salt diffusion

		Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"SaltDiffusion.y\" format=\"ascii\">\n"  << endl;

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
*/




		Hydrosphere_panorama_vts_File <<  "   </PointData>\n" << endl;
		Hydrosphere_panorama_vts_File <<  "   <Points>\n"  << endl;
		Hydrosphere_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;






// Kartesische Coordinates ausschreiben


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


// final file administration
		endpos = Hydrosphere_panorama_vts_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_panorama_vts_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Hydrosphere_panorama_vts_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes"<< endl;
	}


// in case of failing

	if ( Hydrosphere_panorama_vts_File.fail() )
	{
		cout << "***** file ::::: " << Hydrosphere_panorama_vts_File_Name.str() << " ::::: is lost" << endl;
		return;
	}

	Hydrosphere_panorama_vts_File.close();

	if ( Hydrosphere_panorama_vts_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_panorama_vts_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Hydrosphere_panorama_vts_File.fail() )
		cout << "***** file ::::: " << Hydrosphere_panorama_vts_File_Name.str() << " ::::: could not be closed properly" << endl;
}









void PostProcess_Hydrosphere::paraview_vtk_longal ( const string &Name_Bathymetry_File, int &j_longal, int &pressure_iter, Array &h, Array &p, Array &t, Array &u, Array &v, Array &w, Array &c, Array &aux_u, Array &aux_v, Array &Salt_Finger, Array &Salt_Diffusion, Array &Salt_Balance )
{
	double x, y, z, dx, dz;

	stringstream Hydrosphere_longal_File_Name;

	i_max = im;
	k_max = km;

// file administration
	streampos anfangpos, endpos;

//	sprintf ( Hydrosphere_longal_File_Name, "[%s]_Hyd_longal_%i_%i.vtk", Name_Bathymetry_File, j_longal, n );
	Hydrosphere_longal_File_Name <<  "[" << Name_Bathymetry_File << "]_Hyd_longal_" << j_longal << "_" << pressure_iter << ".vtk";
	ofstream Hydrosphere_vtk_longal_File;
	Hydrosphere_vtk_longal_File.precision ( 4 );
	Hydrosphere_vtk_longal_File.setf ( ios::fixed );
	Hydrosphere_vtk_longal_File.open ( Hydrosphere_longal_File_Name.str().c_str(), ios_base::out );

	if ( Hydrosphere_vtk_longal_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_longal_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Hydrosphere_vtk_longal_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_longal_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;




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
//		dx = .025;
//		dz = .05;

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


		Hydrosphere_vtk_longal_File <<  "SCALARS u-Component float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;


// writing u-component

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_longal_File <<  "SCALARS v-Component float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;


// writing v-component

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << v.x[ i ][ j_longal ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_longal_File <<  "SCALARS w-Component float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" <<endl;


// writing w-component

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << w.x[ i ][ j_longal ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;


// writing of temperature

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << t.x[ i ][ j_longal ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_longal_File <<  "SCALARS Salinity float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;


// writing salinity

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << t.x[ i ][ j_longal ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_longal_File <<  "SCALARS SaltFinger float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;


// writing salt finger

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << Salt_Finger.x[ i ][ j_longal ][ k ] << endl;
			}
		}






		Hydrosphere_vtk_longal_File <<  "SCALARS SaltDiffusion float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;


// writing salt diffusion

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << Salt_Diffusion.x[ i ][ j_longal ][ k ] << endl;
			}
		}


		Hydrosphere_vtk_longal_File <<  "SCALARS SaltBalance float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;


// writing salt balance

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << Salt_Balance.x[ i ][ j_longal ][ k ] << endl;
			}
		}





		Hydrosphere_vtk_longal_File <<  "SCALARS Pressure float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;


// writing pressure

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << p.x[ i ][ j_longal ][ k ] << endl;
			}
		}


		Hydrosphere_vtk_longal_File <<  "SCALARS SeaGround float " << 1 << endl;
		Hydrosphere_vtk_longal_File <<  "LOOKUP_TABLE default" << endl;


// writing sea ground

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << h.x[ i ][ j_longal ][ k ] << endl;
			}
		}




		Hydrosphere_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;


// Schreiben der longitudinalen u-w-Zellenstruktur

		for ( int i = 0; i < im; i++ )
		{
			for ( int k = 0; k < km; k++ )

			{
				Hydrosphere_vtk_longal_File << u.x[ i ][ j_longal ][ k ] << " " << y << " " << w.x[ i ][ j_longal ][ k ] << endl;
			}
		}



// end writing



// final file administration
		endpos = Hydrosphere_vtk_longal_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_longal_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Hydrosphere_longal_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes"<< endl;
	}


// in case of failing

	if ( Hydrosphere_vtk_longal_File.fail() )
	{
		cout << "***** file ::::: " << Hydrosphere_longal_File_Name.str() << " ::::: is lost" << endl;
		return;
	}

	Hydrosphere_vtk_longal_File.close();

	if ( Hydrosphere_vtk_longal_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_longal_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Hydrosphere_vtk_longal_File.fail() )
		cout << "***** file ::::: " << Hydrosphere_longal_File_Name.str() << " ::::: could not be closed properly!" << endl;
}




void PostProcess_Hydrosphere::paraview_vtk_radial ( const string &Name_Bathymetry_File, int &i_radial, int &pressure_iter, Array &h, Array &p, Array &t, Array &u, Array &v, Array &w, Array &c, Array &aux_u, Array &aux_v, Array &Salt_Finger, Array &Salt_Diffusion, Array &Salt_Balance, Array_2D &Upwelling, Array_2D &Downwelling, Array_2D &SaltFinger, Array_2D &SaltDiffusion, Array_2D &BottomWater )
{
	double x, y, z, dx, dy;

	stringstream Hydrosphere_radial_File_Name;

	j_max = jm;
	k_max = km;

// file administration
	streampos anfangpos, endpos;

	Hydrosphere_radial_File_Name << "[" <<  Name_Bathymetry_File << "]_Hyd_radial_" << i_radial << "_" << pressure_iter << ".vtk";
	ofstream Hydrosphere_vtk_radial_File;
	Hydrosphere_vtk_radial_File.precision ( 4 );
	Hydrosphere_vtk_radial_File.setf ( ios::fixed );
	Hydrosphere_vtk_radial_File.open ( Hydrosphere_radial_File_Name.str().c_str(), ios_base::out );

	if ( Hydrosphere_vtk_radial_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_radial_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Hydrosphere_vtk_radial_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_radial_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;




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


		Hydrosphere_vtk_radial_File <<  "SCALARS u-Component float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;


// writing u-component

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << u.x[ i_radial ][ j ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_radial_File <<  "SCALARS v-Component float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;


// writing v-component

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_radial_File <<  "SCALARS w-Component float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" <<endl;


// writing w-component

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << w.x[ i_radial ][ j ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;


// writing of temperature

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << t.x[ i_radial ][ j ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_radial_File <<  "SCALARS Salinity float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;


// writing salinity

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << c.x[ i_radial ][ j ][ k ] << endl;
			}
		}



		Hydrosphere_vtk_radial_File <<  "SCALARS Upwelling float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;


// writing upwelling

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << Upwelling.y[ j ][ k ] << endl;
			}
		}




		Hydrosphere_vtk_radial_File <<  "SCALARS Downwelling float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;


// writing downwelling

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << Downwelling.y[ j ][ k ] << endl;
			}
		}


		Hydrosphere_vtk_radial_File <<  "SCALARS BottomWater float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;


// writing BottomWater

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << BottomWater.y[ j ][ k ] << endl;
			}
		}






		Hydrosphere_vtk_radial_File <<  "SCALARS SaltFinger float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;


// writing salt finger

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << SaltFinger.y[ j ][ k ] << endl;
			}
		}




		Hydrosphere_vtk_radial_File <<  "SCALARS SaltDiffusion float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;


// writing salt diffusion

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << SaltDiffusion.y[ j ][ k ] << endl;
			}
		}




		Hydrosphere_vtk_radial_File <<  "SCALARS Pressure float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;


// writing pressure

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << p.x[ i_radial ][ j ][ k ] << endl;
			}
		}


		Hydrosphere_vtk_radial_File <<  "SCALARS SeaGround float " << 1 << endl;
		Hydrosphere_vtk_radial_File <<  "LOOKUP_TABLE default" << endl;


// writing sea ground

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << h.x[ i_radial ][ j ][ k ] << endl;
			}
		}




		Hydrosphere_vtk_radial_File <<  "VECTORS v-w-Cell float" << endl;


// writing zonal v-w cell structure

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				Hydrosphere_vtk_radial_File << v.x[ i_radial ][ j ][ k ] << " " << w.x[ i_radial ][ j ][ k ] << " " << z << endl;
			}
		}



// end writing



// final file administration
		endpos = Hydrosphere_vtk_radial_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_radial_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Hydrosphere_radial_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes"<< endl;
	}


// in case of failing

	if ( Hydrosphere_vtk_radial_File.fail() )
	{
		cout << "***** file ::::: " << Hydrosphere_radial_File_Name.str() << " ::::: is lost" << endl;
		return;
	}

	Hydrosphere_vtk_radial_File.close();

	if ( Hydrosphere_vtk_radial_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_radial_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Hydrosphere_vtk_radial_File.fail() )
		cout << "***** file ::::: " << Hydrosphere_radial_File_Name.str() << " ::::: could not be closed properly!" << endl;
}




void PostProcess_Hydrosphere::paraview_vtk_zonal ( const string &Name_Bathymetry_File, int &k_zonal, int &pressure_iter, Array &h, Array &p, Array &t, Array &u, Array &v, Array &w, Array &c, Array &Salt_Finger, Array &Salt_Diffusion, Array &Salt_Balance )
{
	double x, y, z, dx, dy;

	stringstream Hydrosphere_zonal_File_Name;

	i_max = im;
	j_max = jm;

// file administration
	streampos anfangpos, endpos;

	Hydrosphere_zonal_File_Name <<  "[" << Name_Bathymetry_File << "]_Hyd_zonal_" << k_zonal << "_" << pressure_iter << ".vtk";
	ofstream Hydrosphere_vtk_zonal_File;
	Hydrosphere_vtk_zonal_File.precision ( 4 );
	Hydrosphere_vtk_zonal_File.setf ( ios::fixed );
	Hydrosphere_vtk_zonal_File.open ( Hydrosphere_zonal_File_Name.str().c_str(), ios_base::out );

	if ( Hydrosphere_vtk_zonal_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_zonal_File_Name.str() << " ::::: could be opened" << endl;
		anfangpos = Hydrosphere_vtk_zonal_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_zonal_File_Name.str() << " ::::: starts at ::::::: " << anfangpos << endl;




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


		Hydrosphere_vtk_zonal_File <<  "SCALARS u-Component float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;


// writing u-component

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] << endl;
			}
		}



		Hydrosphere_vtk_zonal_File <<  "SCALARS v-Component float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;


// writing v-component

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << v.x[ i ][ j ][ k_zonal ] << endl;
			}
		}



		Hydrosphere_vtk_zonal_File <<  "SCALARS w-Component float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" <<endl;


// writing w-component

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << w.x[ i ][ j ][ k_zonal ] << endl;
			}
		}



		Hydrosphere_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;


// writing of temperature

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << t.x[ i ][ j ][ k_zonal ] << endl;
			}
		}



		Hydrosphere_vtk_zonal_File <<  "SCALARS Salinity float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;


// writing salinity

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << c.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


		Hydrosphere_vtk_zonal_File <<  "SCALARS SaltFinger float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;


// writing salt finger

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << Salt_Finger.x[ i ][ j ][ k_zonal ] << endl;
			}
		}




		Hydrosphere_vtk_zonal_File <<  "SCALARS SaltDiffusion float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;


// writing salt diffusion

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << Salt_Diffusion.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


		Hydrosphere_vtk_zonal_File <<  "SCALARS SaltBalance float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;


// writing salt balance

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << Salt_Balance.x[ i ][ j ][ k_zonal ] << endl;
			}
		}







		Hydrosphere_vtk_zonal_File <<  "SCALARS Pressure float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;


// writing pressure

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << p.x[ i ][ j ][ k_zonal ] << endl;
			}
		}


		Hydrosphere_vtk_zonal_File <<  "SCALARS SeaGround float " << 1 << endl;
		Hydrosphere_vtk_zonal_File <<  "LOOKUP_TABLE default" << endl;


// writing sea ground

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << h.x[ i ][ j ][ k_zonal ] << endl;
			}
		}




		Hydrosphere_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;


// writing zonal u-v cell structure

		for ( int i = 0; i < im; i++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Hydrosphere_vtk_zonal_File << u.x[ i ][ j ][ k_zonal ] << " " << v.x[ i ][ j ][ k_zonal ] << " " << z << endl;
			}
		}



// end writing



// final file administration
		endpos = Hydrosphere_vtk_zonal_File.tellp();
		cout << "***** file ::::: " << Hydrosphere_zonal_File_Name.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Hydrosphere_zonal_File_Name.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes"<< endl;
	}


// in case of failing

	if ( Hydrosphere_vtk_zonal_File.fail() )
	{
		cout << "***** file ::::: " << Hydrosphere_zonal_File_Name.str() << " ::::: is lost!" << endl;
		return;
	}

	Hydrosphere_vtk_zonal_File.close();

	if ( Hydrosphere_vtk_zonal_File.good() )
	{
		cout << "***** file ::::: " << Hydrosphere_zonal_File_Name.str() << " ::::: could be closed after writing" << endl;
		cout << endl << endl << endl;
	}

	if ( Hydrosphere_vtk_zonal_File.fail() )
		cout << "***** file ::::: " << Hydrosphere_zonal_File_Name.str() << " ::::: could not be closed properly" << endl;
}



void PostProcess_Hydrosphere::Hydrosphere_SequelFile_write ( const string &Name_Bathymetry_File, int &n, int &iter_RB, double &time, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &aux_u, Array &aux_v, Array &aux_w, Array_2D &t_j, Array_2D &c_j )
{
	stringstream Name_Sequel_File;

	streampos anfangpos, endpos;

// file administration
	Name_Sequel_File << "[" << Name_Bathymetry_File << "]_Sequel_Hyd.seq";
	ofstream Sequel_File;
	Sequel_File.precision ( 4 );
	Sequel_File.setf ( ios::fixed );
	Sequel_File.open ( Name_Sequel_File.str().c_str(), ios_base::out );


	if ( Sequel_File.good() )
	{
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: could be opened" << endl;
		cout << endl;
		anfangpos = Sequel_File.tellp();
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: starts at ::::::: " << anfangpos << endl;


// begin writing

		cout << "\n\n***** Hydrosphere_SequelFile_write:   n = " << n << "  iter_RB = " << iter_RB << "  time = " << time << endl << endl << endl;
		Sequel_File << n << " " << iter_RB << " " << time << endl;

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


		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Sequel_File << t_j.y[ j ][ k ] << " " << c_j.y[ j ][ k ] << endl;
			}
		}


// end writing


// final file administration
		endpos = Sequel_File.tellp();
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
		cout << endl;
	}


// in case of failing

	if ( Sequel_File.fail() )
	{
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: is lost!" << endl;
		return;
	}

	Sequel_File.close();

	if ( Sequel_File.good() )
	{
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: could be closed after writing" << endl;
		cout << endl;
	}

	if ( Sequel_File.fail() )
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: could not be closed properly!" << endl << endl << endl;
}










void PostProcess_Hydrosphere::Hydrosphere_SequelFile_read ( const string &Name_Bathymetry_File, int &n, int &iter_BC, double &time, Array_1D &rad, Array_1D &the, Array_1D &phi, Array &h, Array &t, Array &u, Array &v, Array &w, Array &c, Array &tn, Array &un, Array &vn, Array &wn, Array &cn, Array &aux_u, Array &aux_v, Array &aux_w, Array_2D &t_j, Array_2D &c_j )
{
	stringstream Name_Sequel_File;

	streampos anfangpos, endpos;

// file administration
	Name_Sequel_File << "[" << Name_Bathymetry_File << "]_Sequel_Hyd.seq";
	ifstream Sequel_File;
	Sequel_File.open ( Name_Sequel_File.str().c_str(), ios_base::in );
	Sequel_File.seekg ( 0L, ios::beg );
	anfangpos = Sequel_File.tellg ();


	if ( Sequel_File.good() )
	{
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: could be opened." << endl;
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: starts at ::::::: " << anfangpos << endl;


// begin reading

		Sequel_File >> n;
		Sequel_File >> iter_BC;
		Sequel_File >> time;

		cout << "\n\n***** Hydrosphere_SequelFile_read:   n = " << n << "  iter_BC = " << iter_BC << "  time = " << time << endl << endl << endl;

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


		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				Sequel_File >> t_j.y[ j ][ k ];
				Sequel_File >> c_j.y[ j ][ k ];
			}
		}



// end reading

	Sequel_File.seekg ( 0L, ios::end );
	endpos = Sequel_File.tellg ();

// final file administration
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes!"<< endl;
	}


// in case of failing

	if ( Sequel_File == NULL )
	{
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: does not exist! ::::::::: " << endl << endl << endl;
		return;
	}

	Sequel_File.close();

	if ( Sequel_File.good() )
	{
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: could be closed." << endl;
		cout << endl;
	}

	if ( Sequel_File.fail() )
		cout << "***** file ::::: " << Name_Sequel_File.str() << " ::::: could not be closed!" << endl;
}





void PostProcess_Hydrosphere::Atmosphere_TransferFile_read ( const string &Name_Bathymetry_File, Array &v, Array &w, Array &p )
{
	stringstream Name_v_w_Transfer_File;

	streampos anfangpos, endpos;

// file administration
	Name_v_w_Transfer_File << "[" << Name_Bathymetry_File << "]_Transfer_Atm.vw";
	ifstream v_w_Transfer_File;
	v_w_Transfer_File.open ( Name_v_w_Transfer_File.str().c_str(), ios_base::in );
	v_w_Transfer_File.seekg ( 0L, ios::beg );
	anfangpos = v_w_Transfer_File.tellg ();


	if ( v_w_Transfer_File.good() )
	{
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: could be opened" << endl;
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: starts at ::::::: " << anfangpos << endl;


// begin reading

		for ( int j = 0; j < jm; j++ )
		{
			for ( int   k = 0; k < km; k++ )
			{
				v_w_Transfer_File >> v.x[ im-1 ][ j ][ k ];
				v_w_Transfer_File >> w.x[ im-1 ][ j ][ k ];
				v_w_Transfer_File >> p.x[ im-1 ][ j ][ k ];
			}
		}


// end reading

	v_w_Transfer_File.seekg ( 0L, ios::end );
	endpos = v_w_Transfer_File.tellg ();

// final file administration
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: ends at ::::::::: " << endpos << endl;
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: has the length of ::::: " << endpos - anfangpos << " bytes"<< endl;
	}


// in case of failing

	if ( v_w_Transfer_File == NULL )
	{
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: does not exist ::::::::: " << endl << endl << endl;
		return;
	}

	v_w_Transfer_File.close();

	if ( v_w_Transfer_File.good() )
	{
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: could be closed" << endl;
	}

	if ( v_w_Transfer_File.fail() )
		cout << "***** file ::::: " << Name_v_w_Transfer_File.str() << " ::::: could not be closed" << endl;
}






void PostProcess_Hydrosphere::Hydrosphere_PlotData ( const string &Name_Bathymetry_File, Array &v, Array &w, Array &t, Array &c, Array_2D &BottomWater, Array_2D & Upwelling, Array_2D & Downwelling )
{
	stringstream Name_PlotData_File;

	streampos anfangpos, endpos;

// file administration
	Name_PlotData_File << "[" << Name_Bathymetry_File << "]_PlotData_Hyd.xyz";
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





		PlotData_File << " latitude ( ° )" << "  ,  " << "longitude ( ° )" << "  ,    " << "v-velocity ( 0.724 * x * m/s )" << "  ,    " << "w-velocity ( 0.724 * x * m/s )" << "   ,   " << "temperature ( 273.15 * x - 273.15 )" << "  ,   " << "salinity ( 35 * x * psu )" << "   ,   " << "bottom_water" << "   ,   " <<  "upwelling" << "   ,   " <<  "downwelling" << endl;

		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				PlotData_File << j << " " << k << " " << v.x[ im-1 ][ j ][ k ] * .724 << " " << w.x[ im-1 ][ j ][ k ] * .724 << " " << t.x[ im-1 ][ j ][ k ] * 273.15 - 273.15 << " " << c.x[ im-1 ][ j ][ k ] * 35. << " " << BottomWater.y[ j ][ k ] << " " << Upwelling.y[ j ][ k ] << "   " << Downwelling.y[ j ][ k ] << " " <<  endl;

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

