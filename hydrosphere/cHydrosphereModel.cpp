/*
 * Ocean General Circulation Modell ( OGCM ) applied to laminar flow
 * program for the computation of geo-atmospherical circulating flows in a spherical shell
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )
*/
#include "cHydrosphereModel.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "BC_Hyd.h"
#include "BC_Bath_Hyd.h"
#include "BC_Thermohalin.h"
#include "Accuracy_Hyd.h"
#include "PostProcess_Hyd.h"
#include "Pressure_Hyd.h"
#include "MinMax_Hyd.h"
#include "Results_Hyd.h"
#include "Utils.h"

#include "Config.h"
#include "tinyxml2.h"
#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

// Earth's radius is r_earth = 6731 km compares to 6.731 [ / ]
// for 6 km expansion of the area of circulation compares to 0.02 [ / ] with 40 steps of size 0.0005 

// Definition of meridional and longitudinal step sizes 
// for example: dthe = the_Grad / pi180 = 1.125 / 57.3 = 0.01963

// maximum velocity on the sea surface  w_max = 0.29 [ / ] compares to 0.21 m/s = 0.78 km/h as annual mean 
// mean velocity for sea level is 0.5 to 1 km/h compares to 0.14 to 0.28 m/s

// maximum temperature of earth's surface at equator t_max = 1.1355 compares to 37° C compares to 310 K
// maximum temperature of earth's surface at equator t_max = 1.0974 compares to 27° C compares to 300 K
// minimum temperature at the poles t_pol = .7803 compares to -60° C compares to 213.15 K
// minimum temperature in the deep ocean t_deep_ocean = 1.0146 compares to 4° C compares to 277.15 K
// temperature t_0 = 1.000 compares to 0° C compares to 273,15 K
// temperature t_0 = 0.003661 compares to 1° C compares to 1 K

// 1 PSU ( Practical Salt Unit ) = 1 g/kg, means g of salt per kg sweet water
// mass of water compares to 1.0, rate of salt compares to 0.0346
// c_0 compares to the total mass for mean salinity of 34.6 psu or dimensionsless 1.
// for c = 0.9249 compares to a salinity of 32.0 psu
// for c = 0.9682 compares to a salinity of 33.5 psu
// for c = 1.0000 compares to a salinity of 34.6 psu
// for c = 1.0983 compares to a salinity of 38.0 psu

const float cHydrosphereModel::dr = 0.025;  // 0.025 x 40 = 1.0 compares to 16 km : 40 = 150 m for 1 radial step
const float cHydrosphereModel::dt = 0.000001; // time step satisfies the CFL condition

const float cHydrosphereModel::pi180 = 180./ M_PI;  // pi180 = 57.3
const float cHydrosphereModel::the_degree = 1.;   // compares to 1° step size laterally
const float cHydrosphereModel::phi_degree = 1.;  // compares to 1° step size longitudinally

const float cHydrosphereModel::dthe = the_degree / pi180; // dthe = the_degree / pi180 = 1.0 / 57.3 = 0.01745, 180 * .01745 = 3.141
const float cHydrosphereModel::dphi = phi_degree / pi180; // dphi = phi_degree / pi180 = 1.0 / 57.3 = 0.01745, 360 * .01745 = 6.282

cHydrosphereModel::cHydrosphereModel() :
    has_printed_welcome_msg(false)
{
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.
    PythonStream::OverrideCout();

    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

    // set default configuration
    SetDefaultConfig();
}

cHydrosphereModel::~cHydrosphereModel() { }

#include "cHydrosphereDefaults.cpp.inc"

void cHydrosphereModel::LoadConfig(const char *filename) {
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    try{
        if (err) {
            doc.PrintError();
            throw std::invalid_argument(std::string("unable to load config file:  ") + filename);
        }

        XMLElement *atom = doc.FirstChildElement("atom"), *elem_common = NULL, *elem_hydrosphere = NULL;
        if (!atom) {
            throw std::invalid_argument(std::string("Failed to find the 'atom' element in config file: ") + filename);
        }else{
            elem_common = atom->FirstChildElement( "common" );
            if(!elem_common){
                throw std::invalid_argument(std::string(
                    "Failed to find the 'common' element in 'atom' element in config file: ") + filename);
            }
            elem_hydrosphere = atom->FirstChildElement( "hydrosphere" );
            if (!elem_hydrosphere) {
                throw std::invalid_argument(std::string(
                    "Failed to find the 'hydrosphere' element in 'atom' element in config file: ") + filename);
            }
        }

        #include "HydrosphereLoadConfig.cpp.inc"

    }catch(const std::exception &exc){
        std::cerr << exc.what() << std::endl;
        abort();
    }
}

void cHydrosphereModel::RunTimeSlice(int Ma)
{
    // maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km ), 
    // maximum number of overall iterations ( n ),
    // maximum number of inner velocity loop iterations ( velocity_iter_max ),
    // maximum number of outer pressure loop iterations ( pressure_iter_max )
    m_current_time = Ma;

    reset_arrays();

    mkdir(output_path.c_str(), 0777);

    int Ma_max = 300;   // parabolic temperature distribution 300 Ma back
    int Ma_max_half = 150;  // half of time scale

    cout.precision ( 6 );
    cout.setf ( ios::fixed );

    //  Coordinate system in form of a spherical shell
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.Coordinates ( im, r0, dr );
    the.Coordinates ( jm, the0, dthe );
    phi.Coordinates ( km, phi0, dphi );

    int switch_2D = 0;

    // radial expansion of the computational field for the computation of initial values
    int i_max = im - 1;   // corresponds to sea level
    int i_beg = 0;  // compares to an ocean depth of 1000 m with step size 25m, location of the thermocline

    // naming a transfer file for the v- and w-vlocity component of the atmosphere at sea level
    string Name_v_w_Transfer_File;
    stringstream ssName_v_w_Transfer_File;

    //Prepare the temperature and precipitation data file
    string Name_SurfaceTemperature_File  = temperature_file;
    string Name_SurfaceSalinity_File = salinity_file;

    if(Ma != 0 && use_earthbyte_reconstruction){
        Name_SurfaceTemperature_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";
        Name_SurfaceSalinity_File = output_path + "/" + std::to_string(Ma) + "Ma_Reconstructed_Salinity.xyz";

        struct stat info;
        if( stat( Name_SurfaceSalinity_File.c_str(), &info ) != 0){
            std::string cmd_str = "python " + reconstruction_script_path + " " + std::to_string(Ma - time_step) + " " +
                    std::to_string(Ma) + " " + output_path + " " + BathymetrySuffix +" hyd";
            int ret = system(cmd_str.c_str());
            std::cout << " reconstruction script returned: " << ret << std::endl;
        }
    }

    string bathymetry_name = std::to_string(Ma) + BathymetrySuffix;

    if(!has_printed_welcome_msg)
        print_welcome_msg();


    //  class PostProcess for data transport, read and write
    PostProcess_Hydrosphere     read_Transfer ( im, jm, km, output_path );
    read_Transfer.Atmosphere_TransferFile_read ( bathymetry_name, v, w, t, p_dyn, Evaporation_Dalton, Precipitation );


    cout << "***** time slice for the Oceanic Global Circulation Modell ( OGCM ) is:    Ma = " << Ma << " million years" 
        << endl << endl;
    cout << "***** bathymetry/topography given y the x-y-z data set:    " << bathymetry_name.c_str() << endl << endl;


    // class BC_Bathymetry_Hydrosphere for the geometrical boundary condition of the computational area
    BC_Bathymetry_Hydrosphere       depth ( im, jm, km );

    //  class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents 
    //  and the ocean ground
    depth.BC_SeaGround(bathymetry_path + "/" + bathymetry_name, L_hyd, h, Bathymetry);

    // class BC_Hydrosphere for the boundary conditions for the variables at the spherical shell surfaces and the 
    // meridional interface
    BC_Hydrosphere      boundary ( im, jm, km );

    // class Pressure for the subsequent computation of the pressure by a separat Euler equation
    Pressure_Hyd        startPressure ( im, jm, km, dr, dthe, dphi );

    // class Results_MSL_Hyd to compute and show results on the mean sea level, MSL
    Results_Hyd     calculate_MSL ( im, jm, km );

    // class BC_Thermohalin for the initial and boundary conditions of the flow properties
    BC_Thermohalin      oceanflow ( im, jm, km, i_beg, i_max, Ma, Ma_max, Ma_max_half, dr, g, r_0_water, ua, va, wa, ta, 
                                    ca, pa, u_0, p_0, t_0, c_0, cp_w, L_hyd, t_average, t_cretaceous_max, t_equator, 
                                    t_pole, output_path );

    //  surface temperature from World Ocean Atlas 2009 given as boundary condition
    if ( Ma == 0 || use_earthbyte_reconstruction) 
    {
        //oceanflow.BC_Surface_Temperature_NASA ( Name_SurfaceTemperature_File, t );
    }
    //  surface salinity from World Ocean Atlas 2009 given as boundary condition
    if ( Ma == 0 || use_earthbyte_reconstruction) 
    {
        oceanflow.BC_Surface_Salinity_NASA ( Name_SurfaceSalinity_File, c );
    }

    if(Ma == 0) read_IC(velocity_v_file, v.x[im-1], jm, km);
    if(Ma == 0) read_IC(velocity_w_file, w.x[im-1], jm, km);

    iter_cnt_3d = -1;
    if(debug) save_data();
    iter_cnt_3d++;

    //  initial conditions for u-v-w-velocity components following the Ekman spiral
    oceanflow.IC_v_w_EkmanSpiral ( rad, the, h, v, w );

    //  initial conditions for u-velocity component
    oceanflow.IC_u_WestEastCoast ( rad, h, u, v, w, un, vn, wn );

    //  initial conditions for the equatorial currents
    //oceanflow.IC_Equatorial_Currents ( h, u, v, w );

    //  initial conditions for u, v and w velocity components in the circumpolar current
    //if ( Ma <= 41 )     oceanflow.IC_CircumPolar_Current ( h, u, v, w, c ); // Drake passage closed 41 Ma ago

    //  salinity distribution as initial condition in 3 dimensions
    oceanflow.BC_Temperature_Salinity ( h, t, c, p_dyn );

    //  surface pressure computed by surface temperature with gas equation
    oceanflow.BC_Pressure_Density ( p_stat, r_water, r_salt_water, t, c, h );


    //  storing of velocity components, pressure and temperature for iteration start
    store_intermediate_data_2D();
    store_intermediate_data_3D();

    // computation of the ratio ocean to land areas
    calculate_MSL.land_oceanFraction ( h );


    /*******************************************   start of pressure and velocity iterations ******************************/

    double emin = epsres * 100.;
    iter_cnt = 1;
    
    // ::::::::::::::::::::::::::::::::::::::   begin of 2D loop for initial surface conditions: if ( switch_2D == 0 )   ::::::
    if ( switch_2D != 1 )
    {
        // ******   iteration of initial conditions on the surface for the correction of flows close to coasts   ************
        // ******   start of pressure and velocity iterations for the 2D iterational process   *********************************
        //:::::::::::::   begin of pressure loop_2D : if ( pressure_iter_2D > pressure_iter_max_2D )   ::::::::::::::::
        for ( int pressure_iter_2D = 1; pressure_iter_2D <= pressure_iter_max_2D; pressure_iter_2D++)
        {
            // :::::  begin of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::
            for( int velocity_iter_2D = 1; velocity_iter_2D <= velocity_iter_max_2D; velocity_iter_2D++)
            {
                cout << endl << endl;
                cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    2D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
                cout << " 2D OGCM iterational process" << endl;
                cout << " max total iteration number nm = " << nm << endl << endl;

                cout << " present state of the 2D computation " << endl << "  current time slice, number of iterations, \
                    maximum and current number of velocity iterations, maximum and current number of pressure iterations " 
                    << endl << endl << " Ma = " << Ma << "     n = " << iter_cnt << "    velocity_iter_max_2D = " << 
                    velocity_iter_max_2D << "     velocity_iter_2D = " << velocity_iter_2D << "    pressure_iter_max_2D = " 
                    << pressure_iter_max_2D << "    pressure_iter_2D = " << pressure_iter_2D << endl;

                // class BC_Atmosphaere for the geometry of a shell of a sphere
                boundary.RB_theta ( t, u, v, w, p_dyn, c );
                boundary.RB_phi ( t, u, v, w, p_dyn, c );

                oceanflow.Value_Limitation_Hyd ( h, u, v, w, p_dyn, t, c );
                // composition of results
                calculate_MSL.run_data ( i_beg, dr, dthe, L_hyd, u_0, c_0, rad, the, h, u, v, w, c, Salt_Balance, Salt_Finger, 
                    Salt_Diffusion, BuoyancyForce_3D, Upwelling, Downwelling, SaltFinger, SaltDiffusion, BuoyancyForce_2D, 
                    Salt_total, BottomWater );

                // class RungeKutta for the solution of the differential equations describing the flow properties
                solveRungeKutta_2D_Hydrosphere();

                store_intermediate_data_2D();

                iter_cnt++;
            }
            //  :::::::::::::   end of velocity loop_2D: if ( velocity_iter_2D > velocity_iter_max_2D )   ::::::::::::::


            //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
            startPressure.computePressure_2D ( u_0, r_0_water, rad, the, p_dyn, p_dynn, h, aux_v, aux_w );

            // limit of the computation in the sense of time steps
            if ( iter_cnt > nm )
            {
                cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
                break;
            }
        } //end of pressure loop_2D: if ( pressure_iter_2D > pressure_iter_max_2D ) 
        iter_cnt = 1;
        switch_2D = 1;                      // 2D calculations finished
        emin = epsres * 100.;
    }// end of 2D loop for initial surface conditions: if ( switch_2D == 0 )   ::::::::::

    cout << endl << endl;

    iter_cnt_3d = 0;
    if(debug) save_data();
    // ::::   begin of 3D pressure loop : if ( pressure_iter > pressure_iter_max )   ::::::::::::::::::::::::
    for ( int pressure_iter = 1; pressure_iter <= pressure_iter_max; pressure_iter++ )
    {
        //   begin of 3D velocity loop : if ( velocity_iter > velocity_iter_max )   :::::::::::
        for( int velocity_iter = 1; velocity_iter <= velocity_iter_max; velocity_iter++ )
        {
            cout << endl << endl;
            cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
            cout << " 3D OGCM iterational process" << endl;
            cout << " max total iteration number nm = " << nm << endl << endl;
            cout << " present state of the computation " << endl << " current time slice, number of iterations, maximum \
                and current number of velocity iterations, maximum and current number of pressure iterations " << endl 
                << endl << " Ma = " << Ma << "     n = " << iter_cnt << "    velocity_iter_max = " << velocity_iter_max << 
                "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << 
                "    pressure_iter = " << pressure_iter << endl;

            // class RB_Hydrosphaere for the geometry of a shell of a sphere
            boundary.RB_radius ( ca, ta, pa, dr, rad, t, u, v, w, p_dyn, c );
            boundary.RB_theta ( t, u, v, w, p_dyn, c );
            boundary.RB_phi ( t, u, v, w, p_dyn, c );

            // class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of 
            // the continents and the ocean ground
            depth.BC_SolidGround ( ca, ta, pa, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn );

            // surface pressure computed by surface temperature with gas equation
            oceanflow.BC_Pressure_Density ( p_stat, r_water, r_salt_water, t, c, h );

            // preparations for salinity increase due to evaporation and precipitation differences
            oceanflow.BC_Evaporation ( Evaporation_Dalton, Precipitation, h, c, r_water );

            // limiting the increase of flow properties around geometrical peaks and corners
            oceanflow.Value_Limitation_Hyd ( h, u, v, w, p_dyn, t, c );

            // class RungeKutta for the solution of the differential equations describing the flow properties
            solveRungeKutta_3D_Hydrosphere(); 

            // class RB_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of 
            // the continents and the ocean ground
            depth.BC_SolidGround ( ca, ta, pa, h, t, u, v, w, p_dyn, c, tn, un, vn, wn, p_dynn, cn );

            print_min_max();

            // composition of results
            calculate_MSL.run_data ( i_beg, dr, dthe, L_hyd, u_0, c_0, rad, the, h, u, v, w, c, Salt_Balance, Salt_Finger, 
                    Salt_Diffusion, BuoyancyForce_3D, Upwelling, Downwelling, SaltFinger, SaltDiffusion, BuoyancyForce_2D, 
                    Salt_total, BottomWater );

            //  restoring the velocity component and the temperature for the new time step
            store_intermediate_data_3D();

            iter_cnt++;
            iter_cnt_3d++;
            if(debug) save_data();
        }
        //  ::::::  end of velocity loop_3D: if ( velocity_iter > velocity_iter_max )   :::::::::::::::::::::::


        //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
        startPressure.computePressure_3D ( u_0, r_0_water, rad, the, p_dyn, p_dynn, h, aux_u, aux_v, aux_w );

        if( debug && pressure_iter % checkpoint == 0 ){
            write_file(bathymetry_name, output_path);
        }

        //  limit of the computation in the sense of time steps
        if ( iter_cnt > nm )
        {
            cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!" << endl;
            break;
        }
    }// end of pressure loop_3D: if ( pressure_iter > pressure_iter_max )   :::::::::::

    cout << endl << endl;

    write_file(bathymetry_name, output_path, true);
    
    iter_cnt_3d++;
    save_data();
    
    //  final remarks
    cout << endl << "***** end of the Hydrosphere General Circulation Modell ( OGCM ) *****" << endl << endl;

    if ( emin <= epsres )        cout << "***** steady solution reached! *****" << endl;
}


void cHydrosphereModel::reset_arrays()
{
    // 1D arrays
    rad.initArray_1D(im, 1.); // radial coordinate direction
    the.initArray_1D(jm, 2.); // lateral coordinate direction
    phi.initArray_1D(km, 3.); // longitudinal coordinate direction

    // 2D arrays
    Bathymetry.initArray_2D(jm, km, 0.); // Bathymetry in m

    Upwelling.initArray_2D(jm, km, 0.); // upwelling
    Downwelling.initArray_2D(jm, km, 0.); // downwelling
    BottomWater.initArray_2D(jm, km, 0.); // 2D bottom water summed up in a vertical column

    SaltFinger.initArray_2D(jm, km, 0.);      // salt bulge of higher density
    SaltDiffusion.initArray_2D(jm, km, 0.);   // salt bulge of lower density
    Salt_total.initArray_2D(jm, km, 0.);     // rate of salt summed up in a vertical column

    BuoyancyForce_2D.initArray_2D(jm, km, 0.); // radiation balance at the surface

    Evaporation_Dalton.initArray_2D(jm, km, 0.); // evaporation by Penman in [mm/d]
    Precipitation.initArray_2D(jm, km, 0.); // areas of higher precipitation

    // 3D arrays
    h.initArray(im, jm, km, 0.); // bathymetry, depth from sea level

    t.initArray(im, jm, km, ta); // temperature
    u.initArray(im, jm, km, ua); // u-component velocity component in r-direction
    v.initArray(im, jm, km, va); // v-component velocity component in theta-direction
    w.initArray(im, jm, km, wa); // w-component velocity component in phi-direction
    c.initArray(im, jm, km, ca); // water vapour

    tn.initArray(im, jm, km, ta); // temperature new
    un.initArray(im, jm, km, ua); // u-velocity component in r-direction new
    vn.initArray(im, jm, km, va); // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, wa); // w-velocity component in phi-direction new
    cn.initArray(im, jm, km, ca); // water vapour new

    p_dyn.initArray(im, jm, km, pa); // dynamic pressure
    p_dynn.initArray(im, jm, km, pa); // dynamic pressure new
    p_stat.initArray(im, jm, km, pa); // static pressure

    rhs_t.initArray(im, jm, km, 0.); // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.); // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.); // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.); // auxilliar field RHS w-velocity component
    rhs_c.initArray(im, jm, km, 0.); // auxilliar field RHS water vapour

    aux_u.initArray(im, jm, km, 0.); // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.); // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.); // auxilliar field w-velocity component

    Salt_Finger.initArray(im, jm, km, 0.); // salt bulge of higher density
    Salt_Diffusion.initArray(im, jm, km, 0.); // salt bulge of lowerer density and temperature
    Salt_Balance.initArray(im, jm, km, 0.); // +/- salt balance

    r_water.initArray(im, jm, km, 0.); // water density as function of pressure
    r_salt_water.initArray(im, jm, km, 0.); // salt water density as function of pressure and temperature
    BuoyancyForce_3D.initArray(im, jm, km, 0.); // 3D buoyancy force
}

void cHydrosphereModel::Run() 
{
    mkdir(output_path.c_str(), 0777);

    cout << "Output is being written to " << output_path << "\n";

    // write out the config for reproducibility
    // disabled for now
    // std::stringstream output_config_path;
    // output_config_path << output_path << "/config_hyd.xml";
    // WriteConfig(output_config_path.str().c_str());


    if (verbose) {
        cout << endl << endl << endl;
        cout << "***** Hydrosphere General Circulation Model ( OGCM ) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 1 additional transport equations to describe the salinity" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** salinity is part of the Boussinesq approximation" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )" << endl << endl;

        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
    }


    for(int i = time_start; i <= time_end; i+=time_step)
    {
        RunTimeSlice(i);
    }

    //  final remarks
    cout << endl << "***** end of the Hydrosphere General Circulation Modell ( OGCM ) *****" << endl << endl;
    cout << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-hydrospheric circulation *****";
    cout << "\n\n\n\n";
}

void cHydrosphereModel::write_file( std::string &bathymetry_name, string& filepath, bool is_final_result)
{
/*
    //  printout in ParaView and plot files
    //  class PostProcess_Hydrosphaere for the printing of results
    PostProcess_Hydrosphere     write_File ( im, jm, km, filepath );

    int j_longal = 75;
    //  int j_longal = 90;
    write_File.paraview_vtk_longal ( bathymetry_name, j_longal, iter_cnt-1, u_0, r_0_water, h, p_dyn, p_stat, r_water, 
        r_salt_water, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

    //  zonal data along constant longitudes
    int k_zonal = 185;
    //  int k_zonal = 140;
    write_File.paraview_vtk_zonal ( bathymetry_name, k_zonal, iter_cnt-1, u_0, r_0_water, h, p_dyn, p_stat, r_water, r_salt_water, 
            t, u, v, w, c, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );

    //  radial data along constant hight above ground
    int i_radial = 40;
    //  int i_radial = 39;
    write_File.paraview_vtk_radial ( bathymetry_name, i_radial, iter_cnt-1, u_0, t_0, r_0_water, h, p_dyn, p_stat, r_water, 
        r_salt_water, t, u, v, w, c, aux_u, aux_v, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance, 
        Upwelling, Downwelling, SaltFinger, SaltDiffusion, BuoyancyForce_2D, BottomWater, Evaporation_Dalton, 
        Precipitation, Bathymetry );

    //  3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
    if(paraview_panorama_vts){
        write_File.paraview_panorama_vts ( bathymetry_name, iter_cnt-1, u_0, r_0_water, h, t, p_dyn, p_stat, r_water, r_salt_water, 
            u, v, w, c, aux_u, aux_v, aux_w, Salt_Finger, Salt_Diffusion, BuoyancyForce_3D, Salt_Balance );
    }
*/
    //  writing of plot data in the PlotData file
    PostProcess_Hydrosphere     ppa ( im, jm, km, output_path );
    Hydrosphere_PlotData ( bathymetry_name, (is_final_result ? -1 : iter_cnt-1));

}


void  cHydrosphereModel::save_data(){
    struct stat info;
    string path = output_path + "/bin_data/";
    if( stat( path.c_str(), &info ) != 0 ){
        mkdir(path.c_str(), 0777);
    }
    std::ostringstream ss;
    if(iter_cnt_3d == pressure_iter_max * velocity_iter_max + 1)
        ss << "_time_" << m_current_time << "_iter_n";
    else
        ss << "_time_" << m_current_time << "_iter_" << iter_cnt_3d;
    std::string postfix_str = ss.str();

    Array  v_t(im, jm, km, 0), w_t(im, jm, km, 0), c_t(im, jm, km, 0), t_t(im, jm, km, 0);
    for(int i=0; i<im; i++){
        for(int j=0; j<jm; j++){
            for(int k=0; k<km; k++){
                c_t.x[ i ][ j ][ k ] = c.x[ i ][ j ][ k ] * 35.;
                v_t.x[ i ][ j ][ k ] = v.x[ i ][ j ][ k ] * u_0;
                w_t.x[ i ][ j ][ k ] = w.x[ i ][ j ][ k ] * u_0;
                t_t.x[ i ][ j ][ k ] = t.x[ i ][ j ][ k ] * 273.15 - 273.15;
            }
        }
    }
    v_t.save(path + std::string("hyd_v")+postfix_str, im-1);
    w_t.save(path + std::string("hyd_w")+postfix_str, im-1);
    v_t.save(path + std::string("hyd_v")+postfix_str, im-11);
    w_t.save(path + std::string("hyd_w")+postfix_str, im-11);
    v_t.save(path + std::string("hyd_v")+postfix_str, im-21);
    w_t.save(path + std::string("hyd_w")+postfix_str, im-21);
    c_t.save(path + std::string("hyd_s")+postfix_str, im-1);
    t_t.save(path + std::string("hyd_t")+postfix_str, im-1);
}

void cHydrosphereModel::store_intermediate_data_3D(float coeff)
{
    for ( int i = 0; i < im; i++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int k = 0; k < km; k++ )
            {
                un.x[ i ][ j ][ k ] = coeff * u.x[ i ][ j ][ k ];
                vn.x[ i ][ j ][ k ] = coeff * v.x[ i ][ j ][ k ];
                wn.x[ i ][ j ][ k ] = coeff * w.x[ i ][ j ][ k ];
                p_dynn.x[ i ][ j ][ k ] = coeff * p_dyn.x[ i ][ j ][ k ];
                tn.x[ i ][ j ][ k ] = coeff * t.x[ i ][ j ][ k ];
                cn.x[ i ][ j ][ k ] = coeff * c.x[ i ][ j ][ k ];
            }
        }
    }
}

void cHydrosphereModel::store_intermediate_data_2D(float coeff)
{   
    for ( int j = 0; j < jm; j++ )
    {   
        for ( int k = 0; k < km; k++ )
        {   
            vn.x[ im - 1 ][ j ][ k ] = coeff * v.x[ im - 1 ][ j ][ k ];
            wn.x[ im - 1 ][ j ][ k ] = coeff * w.x[ im - 1 ][ j ][ k ]; 
            p_dynn.x[ im - 1 ][ j ][ k ] = coeff * p_dyn.x[ im - 1 ][ j ][ k ];
        }
    }
}

