/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to read and prepare the bathymetric and topografic data
*/


#include <iostream>
#include <cmath>
#include <fstream>

#include "BC_Bathymetry_Atmosphere.h"

using namespace std;




BC_Bathymetry_Atmosphere::BC_Bathymetry_Atmosphere ( int im, int jm, int km )
{
	this -> im = im;
	this -> jm = jm;
	this -> km = km;
}



BC_Bathymetry_Atmosphere::~BC_Bathymetry_Atmosphere(){}



void BC_Bathymetry_Atmosphere::BC_MountainSurface ( const string &Name_Bathymetry_File, double L_atm, Array_2D &aux_2D_h, Array &h, Array &aux_w )
{
	streampos anfangpos_1, endpos_1, anfangpos_2, endpos_2, anfangpos_3, endpos_3, anfangpos_4, endpos_4;

	cout.precision ( 8 );
	cout.setf ( ios::fixed );

// default adjustment, h must be 0 everywhere
		for ( k = 0; k < km; k++ )
		{
			for ( j = 0; j < jm; j++ )
			{
				for ( i = 0; i < im; i++ )
				{
					h.x[ i ][ j ][ k ] = 0.;					// default
				}
			}
		}



// reading data from file Name_Bathymetry_File_Read
	ifstream Name_Bathymetry_File_Read;
	Name_Bathymetry_File_Read.open ( Name_Bathymetry_File.c_str(), ios_base::in );
	Name_Bathymetry_File_Read.seekg ( 0L, ios::beg );
	anfangpos_1 = Name_Bathymetry_File_Read.tellg ();


	if ( Name_Bathymetry_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: could be opened" << endl;
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: begins at ::::::: " << anfangpos_1 << endl;
	}


	j = 0;
	k = 0;

	while ( ( j < jm ) && ( !Name_Bathymetry_File_Read.eof() ) )
	{
		while ( k < km )
		{
			Name_Bathymetry_File_Read >> dummy_1;
			Name_Bathymetry_File_Read >> dummy_2;
			Name_Bathymetry_File_Read >> dummy_3;

			if ( dummy_3 < 0. )
			{
				h.x[ 0 ][ j ][ k ] = 0.;
			}

			if ( dummy_3 >= 0. )
			{
				i = ( im - 1 ) * ( dummy_3 / L_atm );
				i_SL = i;

				for ( i = 0; i <= i_SL; i++ )			h.x[ i ][ j ][ k ] = 1.;
			}
			k++;
		}
	k = 0;
	j++;
	}


// end reading Name_Bathymetry_File

	Name_Bathymetry_File_Read.seekg ( 0L, ios::end );
	endpos_1 = Name_Bathymetry_File_Read.tellg ();

// final file administration

	cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: ends at ::::::::: " << endpos_1 << endl;
	cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: has the length of ::::: " << endpos_1 - anfangpos_1 << " bytes!"<< endl;

// in case of reading error

	if ( Name_Bathymetry_File_Read == NULL )
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: not yet exists! ::::::::: " << endl << endl << endl;
	}

	Name_Bathymetry_File_Read.close();

	if ( Name_Bathymetry_File_Read.good() )
	{
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: could be closed." << endl;
		cout << endl;
	}

	if ( Name_Bathymetry_File_Read.fail() )
		cout << "***** file ::::: " << Name_Bathymetry_File << " ::::: could not be closed!" << endl;

// end reading Name_Bathymetry_File_Read




// Umschreiben der bathymetrischen Daten von -180° - 0° - +180° Koordinatnesystem auf 0°- 360°

		l = 0;

		for ( int k = 180; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					aux_w.x[ i ][ j ][ l ] = h.x[ i ][ j ][ k ];
				}
			}
			l++;
		}



		for ( int k = 0; k < 180; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					aux_w.x[ i ][ j ][ l ] = h.x[ i ][ j ][ k ];
				}
			}
			l++;
		}



		for ( int k = 0; k < km; k++ )
		{
			for ( int j = 0; j < jm; j++ )
			{
				for ( int i = 0; i < im; i++ )
				{
					h.x[ i ][ j ][ k ] = aux_w.x[ i ][ j ][ k ];
				}
			}
		}

// Ende Umschreiben Bathymetrie


// gap control on the surface for the fixed temperature boundary conditions
// preparation for the next time slice
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( ( h.x[ 0 ][ j ][ k ]  == 1. ) || ( aux_2D_h.y[ j ][ k ]  == 1. ) )		// forming the gap, ( h.x ) new - ( aux_2D_h.y ) old bathymetry
				{
					aux_2D_h.y[ j ][ k ] = fabs ( aux_2D_h.y[ j ][ k ] - h.x[ 0 ][ j ][ k ] );							// land in the gap
				}
				else
				{
					aux_2D_h.y[ j ][ k ] = 0.;															// no gap, ocean in both bathymetries
				}
			}
		}
}





void BC_Bathymetry_Atmosphere::BC_IceShield ( int Ma, double t_0, Array &h, Array &t, Array &c, Array &IceLayer, Array_2D &Ice_Balance, Array_2D &Ice_Balance_add )
{
// computation of ice shield following the theorie by Milankowitsch

	min = 0.;
	max = 0.;

// prescribed reference level still to be fixed!

	Hoehe_equi_km = 3000;
	Hoehe_delta = 500;
	Hoehe_equi = Hoehe_equi_km / Hoehe_delta;

// limit of temperature for accumulation 1 and 2 and ablation

	t_equi_Celsius = - 8.;
	t_plus_Celsius = + 2.;
	t_pluss_Celsius = + 12.;
	t_minus_Celsius = - 15.;
	t_minuss_Celsius = - 30.;

	t_equi = ( t_equi_Celsius + t_0 ) / t_0;
	t_plus = ( t_plus_Celsius + t_0 ) / t_0;
	t_pluss = ( t_pluss_Celsius + t_0 ) / t_0;
	t_minus = ( t_minus_Celsius + t_0 ) / t_0;
	t_minuss = ( t_minuss_Celsius + t_0 ) / t_0;




	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( t.x[ 0 ][ j ][ k ] >= t_plus ) break;
			{
				for ( int i = 0; i < im; i++ )																		// ice balance performed over total hight
				{
					if ( ( t.x[ i ][ j ][ k ] >= t_minuss ) && ( t.x[ i ][ j ][ k ] <= t_plus ) )
					{
						Akkumulation_1 = 0.0125 * ( t.x[ i ][ j ][ k ] * t_0 ) - 2.9735;			//	formula valid from -30°C to +2°C by Ruddiman
					}
					else
					{
						Akkumulation_1 = 0.;
					}

					if ( ( t.x[ i ][ j ][ k ] > t_plus ) && ( t.x[ i ][ j ][ k ] <= t_pluss ) )
					{
						Akkumulation_2 = - 0.04 * ( t.x[ i ][ j ][ k ] * t_0 ) + 11.5;					//	formula valid from +2°C to +12°C by Ruddiman
					}
					else
					{
						Akkumulation_2 = 0.;
					}

					if ( ( t.x[ i ][ j ][ k ] >= t_minus ) && ( t.x[ i ][ j ][ k ] <= t_plus ) )
					{
						Ablation = 0.01384 * ( t.x[ i ][ j ][ k ] * t_0 ) * ( t.x[ i ][ j ][ k ] * t_0 ) - 7.14186 * ( t.x[ i ][ j ][ k ] * t_0 ) + 921.2458;	//	formula valid from -15°C to +2°C by Ruddiman
					}
					else
					{
						Ablation = 0.;
					}


					if ( ( t.x[ i ][ j ][ k ] >= t_minuss ) && ( t.x[ i ][ j ][ k ] <= t_plus ) )
					{
						Ice_Balance.y[ j ][ k ] = Akkumulation_1 - Ablation;						// ice balance for accumulation_1 minus ablation
					}


					if ( ( t.x[ i ][ j ][ k ] > t_plus ) && ( t.x[ i ][ j ][ k ] <= t_pluss ) )
					{
						Ice_Balance.y[ j ][ k ] = Akkumulation_2;										// ice balance for accumulation_2
					}

					Ice_Balance_add.y[ j ][ k ] = Ice_Balance_add.y[ j ][ k ] + Ice_Balance.y[ j ][ k ];		// average of ice balance out of the two parts across hight
				}
			}
		}
	}

// searching minimum and maximum values of mean ice balance

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( Ice_Balance_add.y[ j ][ k ] >= max ) 
			{
				max = Ice_Balance_add.y[ j ][ k ];
			}

			if ( Ice_Balance_add.y[ j ][ k ] <= min ) 
			{
				min = Ice_Balance_add.y[ j ][ k ];
			}
		}
	}


	Ice_Balance_add_diff = max - min;														// maximum crosswise extension of ice balance

// presentation of ice shield hight following the local dimensionless ice balance values, maximum value now 1

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			Ice_Hoehe = ( Ice_Balance_add.y[ j ][ k ] - min ) / max * 2.;			// assumed actual ice shield hight
			i_Ice_lauf = ( int ) Ice_Hoehe;

			for ( int i = 0; i < i_Ice_lauf; i++ )
			{
				if ( ( t.x[ 0 ][ j ][ k ] >= t_plus ) || ( h.x[ 0 ][ j ][ k ] == 0. ) ) break;
				{
					IceLayer.x[ i ][ j ][ k ] = 1.;
				}
			}
		}
	}
}




void BC_Bathymetry_Atmosphere::BC_Radiation ( double t_0, double ik, double sigma, double albedo_extra, double epsilon_atmos, Array_2D &Precipitation, Array_2D &Ik, Array_2D &Radiation_Balance_atm, Array_2D &Radiation_Balance_bot, Array_2D &temp_eff_atm, Array_2D &temp_eff_bot, Array_2D &t_j, Array &t )
{
// class element for the computation of the radiation balance
// computation of the local temperature of the radiation balance

	cout.precision ( 6 );
	cout.setf ( ios::fixed );

	pi180 = 180./M_PI;

	j_halb = ( jm -1 ) / 2 + 30;													// position of the sun at 30°S ( 90 + 30 = 120 )
	j_max = jm - 1;

	d_j_halb = ( double ) j_halb;
	d_j_max = ( double ) j_max;

	k_halb = ( km -1 ) / 2;															// position of the sun at 0° oder 180° ( Greenwich )
	k_max = km - 1;

	d_k_halb = ( double ) k_halb;
	d_k_max = ( double ) k_max;


	t_eff_Erde = 255.;														// effectiv radiation temperature of 255 K compares to -18 °C for 30% albedo
																						// comparable with measurements of temperature in 5000 m hight and pressure 550 hPa
																						// sigma * T-Erde 4 = ( 1- albedo ) Ik / 4 = 239 W/m2, valid as reference value

	j_sun = 30;

	for ( int k = 0; k < km; k++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			epsilon_atmos = 2. * ( 1. - 1. / pow ( ( t_j.y[ j ][ k ] * t_0 / t_eff_Erde ), 4. ) );
//			epsilon_atmos = 1.;
 
			Ik.y[ j ][ k ] = ik * sin ( ( ( double ) j - j_sun ) / pi180 ) * sin ( ( ( double ) k - 90 ) / pi180 );

			r = sqrt ( ( ( double ) j - 30 ) * ( ( double ) j - 30 ) + ( ( double ) k - 180 ) * ( ( double ) k - 180 ) );

			if ( Ik.y[ j ][ k ] < 400. )
//			if ( ( j >= ( int ) r ) || ( k >= ( int ) r ) )
			{
				j_r = j;
				k_r = k;
//				Ik.y[ j ][ k ] = Ik.y[ j ][ k_r ];
				Ik.y[ j ][ k ] = 400.;
			}

			Radiation_Balance_bot.y[ j ][ k ] = ( 1. - albedo_extra ) / ( 1. - epsilon_atmos / 2. ) * .25 * Ik.y[ j ][ k ];
			Radiation_Balance_atm.y[ j ][ k ] = ( 1. - albedo_extra ) / ( 2. - epsilon_atmos ) * .25 * Ik.y[ j ][ k ];
			temp_eff_bot.y[ j ][ k ] = pow ( Radiation_Balance_bot.y[ j ][ k ] / sigma, 1. / 4. );
			temp_eff_atm.y[ j ][ k ] = pow ( Radiation_Balance_atm.y[ j ][ k ] / sigma, 1. / 4. );
			temp_eff_bot.y[ j ][ k ] = t.x[ 0 ][ j ][ k ] = temp_eff_bot.y[ j ][ k ] / t_0;
		}
	}
}





void BC_Bathymetry_Atmosphere::BC_SolidGround ( Array &h, Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &pn, Array &cn, Array &co2n, Array &fup, Array &fvp, Array &fwp, Array &ftp, Array &fcp, Array &fco2p )
{
// boundary conditions for the total solid ground

	for ( int i = 0; i < im; i++ )
	{
		for ( int j = 0; j < jm; j++ )
		{
			for ( int k = 0; k < km; k++ )
			{
				if ( h.x[ i ][ j ][ k ] == 1. )
				{
//					p.x[ i ][ j ][ k ] = pn.x[ i ][ j ][ k ] = pa;

					u.x[ i ][ j ][ k ] = un.x[ i ][ j ][ k ] = 0.;
					v.x[ i ][ j ][ k ] = vn.x[ i ][ j ][ k ] = 0.;
					w.x[ i ][ j ][ k ] = wn.x[ i ][ j ][ k ] = 0.;

					c.x[ i ][ j ][ k ] = cn.x[ i ][ j ][ k ] = 0.;
					co2.x[ i ][ j ][ k ] = co2n.x[ i ][ j ][ k ] = 0.;

					fup.x[ i ][ j ][ k ] = 0.;
					fvp.x[ i ][ j ][ k ] = 0.;
					fwp.x[ i ][ j ][ k ] = 0.;
					ftp.x[ i ][ j ][ k ] = 0.;
					fcp.x[ i ][ j ][ k ] = 0.;
					fco2p.x[ i ][ j ][ k ] = 0.;
				}
			}
		}
	}
}







void BC_Bathymetry_Atmosphere::BC_SolidGround_2D ( Array &h, Array &t, Array &u, Array &v, Array &w, Array &p, Array &c, Array &co2, Array &tn, Array &un, Array &vn, Array &wn, Array &pn, Array &cn, Array &co2n, Array &fup, Array &fvp, Array &fwp, Array &ftp, Array &fcp, Array &fco2p )
{
// boundary conditions for the total solid ground

	for ( int j = 0; j < jm; j++ )
	{
		for ( int k = 0; k < km; k++ )
		{
			if ( h.x[ 0 ][ j ][ k ] == 1. )
			{
//				p.x[ 0 ][ j ][ k ] = pn.x[ 0 ][ j ][ k ] = pa;

				u.x[ 0 ][ j ][ k ] = un.x[ 0 ][ j ][ k ] = 0.;
				v.x[ 0 ][ j ][ k ] = vn.x[ 0 ][ j ][ k ] = 0.;
				w.x[ 0 ][ j ][ k ] = wn.x[ 0 ][ j ][ k ] = 0.;

				c.x[ 0 ][ j ][ k ] = cn.x[ 0 ][ j ][ k ] = 0.;
				co2.x[ 0 ][ j ][ k ] = co2n.x[ 0 ][ j ][ k ] = 0.;

				fup.x[ 0 ][ j ][ k ] = 0.;
				fvp.x[ 0 ][ j ][ k ] = 0.;
				fwp.x[ 0 ][ j ][ k ] = 0.;
				ftp.x[ 0 ][ j ][ k ] = 0.;
				fcp.x[ 0 ][ j ][ k ] = 0.;
				fco2p.x[ 0 ][ j ][ k ] = 0.;
			}
		}
	}
}

