#include "cAtmosphereModel.h"

#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "BC_Atm.h"
#include "BC_Bath_Atm.h"
#include "BC_Thermo.h"
#include "Accuracy_Atm.h"
#include "RHS_Atm.h"
#include "RungeKutta_Atm.h"
#include "Print_Atm.h"
#include "PostProcess_Atm.h"
#include "Pressure_Atm.h"
#include "Restore_Atm.h"
#include "Results_Atm.h"
#include "MinMax_Atm.h"

#include "tinyxml2.h"

#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;

cAtmosphereModel::cAtmosphereModel() {
    // Python and Notebooks can't capture stdout from this module. We override
    // cout's streambuf with a class that redirects stdout out to Python.
    PythonStream::OverrideCout();

    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);

    // set default configuration
    SetDefaultConfig();
}

cAtmosphereModel::~cAtmosphereModel() { }

#include "cAtmosphereDefaults.cpp.inc"

void cAtmosphereModel::RunTimeSlice(int Ma) {
    // maximum numbers of grid points in r-, theta- and phi-direction ( im, jm, km )
    // maximum number of overall iterations ( n )
    // maximum number of inner velocity loop iterations ( velocity_iter_max )
    // maximum number of outer pressure loop iterations ( pressure_iter_max )

    mkdir(output_path.c_str(), 0777);

    const int im = 41, jm = 181, km = 361, nm = 200;

    int j_res = 0, k_res = 0;


    // TODO: lots of scoping issues here and alias issue with min(); too much to safely sort out right now
    double residuum = 0, residuum_old = 0, min = 0;

    constexpr double pi180 = 180./M_PI;                       // pi180 = 57.3
    constexpr double the_degree = 1.;                         // compares to 1° step size laterally
    constexpr double phi_degree = 1.;                         // compares to 1° step size longitudinally

    double dthe = the_degree / pi180;           //dthe = the_degree / pi180 = 1.0 / 57.3 = 0.01745, 180 * .01745 = 3.141
    double dphi = phi_degree / pi180;               //dphi = phi_degree / pi180 = 1.0 / 57.3 = 0.01745, 360 * .01745 = 6.282

    const double the0 = 0.;                                       // North Pole
    const double phi0 = 0.;                                       // zero meridian in Greenwich
    const double r0 = 6.731;                                  // earth's radius is r_earth = 6731 km compares to 6.731 [ / ] * 1000 km, circumference of the earth 40074 km

    const double coeff_mmWS = r_air / r_water_vapour; // coeff_mmWS = 1.2041 / 0.0094 [ kg/m³ / kg/m³ ] = 128,0827 [ / ]

    int *im_tropopause = new int [ jm ];            // location of the tropopause

    //  class Array for 1-D, 2-D and 3-D field declarations

    // 1D arrays
    Array_1D rad(im, 1.); // radial coordinate direction
    Array_1D the(jm, 2.); // lateral coordinate direction
    Array_1D phi(km, 3.); // longitudinal coordinate direction

    // 2D arrays
    Array_2D Vegetation(jm, km, 0.); // vegetation via precipitation

    Array_2D LatentHeat(jm, km, 0.); // areas of higher latent heat
    Array_2D Condensation(jm, km, 0.); // areas of higher condensation
    Array_2D Evaporation(jm, km, 0.); // areas of higher evaporation

    Array_2D Precipitation(jm, km, 0.); // areas of higher precipitation
    Array_2D precipitable_water(jm, km, 0.); // areas of precipitable water in the air
    Array_2D precipitation_NASA(jm, km, 0.); // surface precipitation from NASA

    Array_2D Ice_Balance(jm, km, 0.); // rate of the ice shield
    Array_2D Ice_Balance_add(jm, km, 0.); // addition of all ice layers

    Array_2D Ik(jm, km, 0.); // direct sun radiation, short wave
    Array_2D Radiation_Balance(jm, km, 0.); // radiation balance at the surface
    Array_2D Radiation_Balance_par(jm, km, 0.); // radiation balance at the surface parabolic distribution
    Array_2D Radiation_Balance_atm(jm, km, 0.); // radiation balance of the atmosphere
    Array_2D Radiation_Balance_bot(jm, km, 0.); // radiation balance of the ground

    Array_2D t_j(jm, km, 0.); // 2D temperature at the surface
    Array_2D temp_eff_atm(jm, km, 0.); // effektive temperature in the atmosphere
    Array_2D temp_eff_bot(jm, km, 0.); // effektive temperature on the ground
    Array_2D temp_rad(jm, km, 0.); // effektive temperature based on radiation

    Array_2D albedo(jm, km, 0.); // albedo = reflectivity
    Array_2D epsilon(jm, km, 0.); // epsilon = absorptivity

    Array_2D Q_Radiation(jm, km, 0.); // heat from the radiation balance in [W/m2]
    Array_2D Q_Evaporation(jm, km, 0.); // evaporation heat of water by Kuttler
    Array_2D Q_latent(jm, km, 0.); // latent heat from bottom values by the energy transport equation
    Array_2D Q_sensible(jm, km, 0.); // sensible heat from bottom values by the energy transport equation
    Array_2D Q_bottom(jm, km, 0.); // difference by Q_Radiation - Q_latent - Q_sensible

    Array_2D Evaporation_Haude(jm, km, 0.); // evaporation by Haude in [mm/d]
    Array_2D Evaporation_Penman(jm, km, 0.); // evaporation by Penman in [mm/d]

    Array_2D MaxCloud(jm, km, 0.); // maximum cloud water mass in a vertical column
    Array_2D MaxIce(jm, km, 0.); // maximum cloud ice mass in a vertical column

    Array_2D co2_total(jm, km, 0.); // areas of higher co2 concentration

    Array_2D aux_2D_v(jm, km, 0.); // auxilliar field v
    Array_2D aux_2D_w(jm, km, 0.); // auxilliar field w

    // 3D arrays
    Array h(im, jm, km, 0.); // bathymetry, depth from sea level

    Array t(im, jm, km, ta); // temperature
    Array u(im, jm, km, ua); // u-component velocity component in r-direction
    Array v(im, jm, km, va); // v-component velocity component in theta-direction
    Array w(im, jm, km, wa); // w-component velocity component in phi-direction
    Array c(im, jm, km, ca); // water vapour
    Array cloud(im, jm, km, 0.); // cloud water
    Array ice(im, jm, km, 0.); // cloud ice
    Array co2(im, jm, km, coa); // CO2

    Array tn(im, jm, km, ta); // temperature new
    Array un(im, jm, km, ua); // u-velocity component in r-direction new
    Array vn(im, jm, km, va); // v-velocity component in theta-direction new
    Array wn(im, jm, km, wa); // w-velocity component in phi-direction new
    Array cn(im, jm, km, ca); // water vapour new
    Array cloudn(im, jm, km, 0.); // cloud water new
    Array icen(im, jm, km, 0.); // cloud ice new
    Array co2n(im, jm, km, coa); // CO2 new

    Array p_dyn(im, jm, km, pa); // dynamic pressure
    Array p_stat(im, jm, km, pa); // static pressure

    Array rhs_t(im, jm, km, 0.); // auxilliar field RHS temperature
    Array rhs_u(im, jm, km, 0.); // auxilliar field RHS u-velocity component
    Array rhs_v(im, jm, km, 0.); // auxilliar field RHS v-velocity component
    Array rhs_w(im, jm, km, 0.); // auxilliar field RHS w-velocity component
    Array rhs_c(im, jm, km, 0.); // auxilliar field RHS water vapour
    Array rhs_cloud(im, jm, km, 0.); // auxilliar field RHS cloud water
    Array rhs_ice(im, jm, km, 0.); // auxilliar field RHS cloud ice
    Array rhs_co2(im, jm, km, 0.); // auxilliar field RHS CO2

    Array aux_u(im, jm, km, 0.); // auxilliar field u-velocity component
    Array aux_v(im, jm, km, 0.); // auxilliar field v-velocity component
    Array aux_w(im, jm, km, 0.); // auxilliar field w-velocity component
    Array aux_p(im, jm, km, pa); // auxilliar field p

    Array Latency(im, jm, km, 0.); // latent heat
    Array Q_Sensible(im, jm, km, 0.); // sensible heat
    Array IceLayer(im, jm, km, 0.); // ice shield
    Array t_cond_3D(im, jm, km, 0.); // condensation temperature
    Array t_evap_3D(im, jm, km, 0.); // evaporation temperature
    Array BuoyancyForce(im, jm, km, 0.); // buoyancy force, Boussinesque approximation
    Array epsilon_3D(im, jm, km, 0.); // emissivity/ absorptivity
    Array radiation_3D(im, jm, km, 0.); // radiation

    Array P_rain(im, jm, km, 0.); // rain precipitation mass rate
    Array P_snow(im, jm, km, 0.); // snow precipitation mass rate
    Array S_v(im, jm, km, 0.); // water vapour mass rate due to category two ice scheme
    Array S_c(im, jm, km, 0.); // cloud water mass rate due to category two ice scheme
    Array S_i(im, jm, km, 0.); // cloud ice mass rate due to category two ice scheme
    Array S_r(im, jm, km, 0.); // rain mass rate due to category two ice scheme
    Array S_s(im, jm, km, 0.); // snow mass rate due to category two ice scheme

    //  cout << endl << " ***** printout of 3D-field temperature ***** " << endl << endl;
    //  t.printArray();

    //  cout << endl << " ***** printout of 2D-field vegetation ***** " << endl << endl;
    //  Vegetation.printArray_2D();

    //  cout << endl << " ***** printout of 1D-field radius ***** " << endl << endl;
    //  rad.printArray_1D();

    cout.precision ( 6 );
    cout.setf ( ios::fixed );

    //  Coordinate system in form of a spherical shell
    //  rad for r-direction normal to the surface of the earth, the for lateral and phi for longitudinal direction
    rad.Coordinates ( im, r0, dr );
    the.Coordinates ( jm, the0, dthe );
    phi.Coordinates ( km, phi0, dphi );

    //  initial values for the number of computed steps and the time
    int n = 0;
    double time = dt;
    int pressure_iter = 1;
    int pressure_iter_2D = 1;
    int switch_2D = 0;
    residuum_old = 0.;

    // radial expansion of the computational field for the computation of initial values
    int i_max = 32;         // corresponds to about 16 km above sea level, maximum hight of the tropopause at equator
    int i_beg = 16;         // corresponds to about 8 km above sea level, maximum hight of the tropopause at poles

    // naming a file to read the surface temperature of the modern world
    string Name_SurfaceTemperature_File;
    stringstream ssNameSurfaceTemperature;
    ssNameSurfaceTemperature << output_path << "SurfaceTemperature.xyz";
    Name_SurfaceTemperature_File = ssNameSurfaceTemperature.str();

    // naming a file to read the surface precipitation by NASA
    string Name_SurfacePrecipitation_File;
    stringstream ssNameSurfacePrecipitation;
    ssNameSurfacePrecipitation << "SurfacePrecipitation_NASA.xyz";
    Name_SurfacePrecipitation_File = ssNameSurfacePrecipitation.str();

    // TODO we're assuming that our output dir has been set
    cout << "Output is being written to " << output_path << "\n";
    cout << "Ma = " << Ma << "\n";

    // stringstream My;
    string bathymetry_name = std::to_string(Ma) + BathymetrySuffix;
    string bathymetry_filepath = bathymetry_path + "/" + bathymetry_name;
    // string Name_netCDF_File = std::to_string(Ma) + "Ma_atmosphere.nc";

    PostProcess_Atmosphere read_File(im, jm, km, output_path);
    read_File.Atmosphere_SequelFile_read(bathymetry_name, n, time, rad, the, phi, h, t, u, v, w, c, co2, tn, un, vn, wn, cn, co2n);
    // n++;

    // initialization of the bathymetry/topography

    // class BC_Bathymetry_Atmosphere for the geometrical boundary condition of the computational area
    BC_Bathymetry_Atmosphere LandArea(im, jm, km, co2_vegetation, co2_land, co2_ocean);

    // topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
    LandArea.BC_MountainSurface(bathymetry_filepath, L_atm, h, aux_w);

    if (IceShield) {
        LandArea.BC_IceShield(Ma, t_0, h, t, c, IceLayer, Ice_Balance, Ice_Balance_add);
    }

    // class calls for the solution of the flow properties

    // class BC_Atmosphere for the boundary conditions for the variables at the spherical shell surfaces and the meridional interface
    BC_Atmosphere boundary(im, jm, km, t_tropopause);

    // class RHS_Atmosphere for the preparation of the time independent right hand sides of the Navier-Stokes equations
    RHS_Atmosphere prepare(im, jm, km, dt, dr, dthe, dphi, re, ec, sc_WaterVapour, sc_CO2, g, pr, omega, coriolis, centrifugal, WaterVapour, buoyancy, CO2, gam, sigma, Lambda);
    RHS_Atmosphere prepare_2D (jm, km, dthe, dphi, re, omega, coriolis, centrifugal);

    // class RungeKutta_Atmosphere for the explicit solution of the Navier-Stokes equations
    RungeKutta_Atmosphere result(im, jm, km, dt, dr, dphi, dthe);

    // class Pressure for the subsequent computation of the pressure by a separat Euler equation
    Pressure_Atm startPressure(im, jm, km, dr, dthe, dphi);

    // class Restore to restore the iterational values from new to old
    Restore oldnew(im, jm, km);

    // class Results_MSL_Atm to compute and show results on the mean sea level, MSL
    Results_MSL_Atm calculate_MSL(im, jm, km, g, ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_extra, lv, ls, cp_l, L_atm, dt, dr, dthe, dphi, r_air, R_Air, r_water, r_water_vapour, R_WaterVapour, co2_vegetation, co2_ocean, co2_land, gam, t_pole, t_cretaceous);

    // configuration of the initial and boundary conditions for the temperature, CO2 und water vapour on land and ocean surfaces

    // class BC_Thermo for the initial and boundary conditions of the flow properties
    BC_Thermo circulation(im, jm, km, i_beg, i_max, RadiationModel, sun, declination, sun_position_lat, sun_position_lon, Ma, Ma_max, Ma_max_half, dr, dthe, dphi, g, ep, hp, u_0, p_0, t_0, c_0, sigma, albedo_extra, epsilon_extra, lv, cp_l, L_atm, r_air, R_Air, r_water_vapour, R_WaterVapour, co2_0, co2_cretaceous, co2_vegetation, co2_ocean, co2_land, ik, c_tropopause, co2_tropopause, c_ocean, c_land, t_average, co2_average, co2_pole, t_cretaceous, t_cretaceous_max, radiation_ocean, radiation_pole, radiation_equator, t_land, t_tropopause, t_equator, t_pole, gam, set_sun_position, verbose);

    // class element for the tropopause location as a parabolic distribution from pole to pole
    circulation.TropopauseLocation(im_tropopause);

    //  class element for the surface temperature from World Ocean Atlas 2009 given as boundary condition
    //  if ( Ma == 0 ) circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );
    //  circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );

    //  class element for the surface precipitation from NASA for comparison
    //  if ( Ma == 0 ) circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );
    //  circulation.BC_Surface_Temperature ( Name_SurfaceTemperature_File, t );

    // class element for the parabolic temperature distribution from pol to pol, maximum temperature at equator
    circulation.BC_Temperature(h, t, p_dyn, p_stat);
    t_cretaceous = circulation.out_temperature();

    // class element for the surface pressure computed by surface temperature with gas equation
    circulation.BC_Pressure(p_stat, t, h);

    // parabolic water vapour distribution from pol to pol, maximum water vapour volume at equator
    circulation.BC_WaterVapour(h, t, p_stat, c, cloud, ice, P_rain, P_snow, precipitation_NASA);

    //  circulation.IC_Cloud_Ice ( h, c, t, p_stat, cloud, ice );

    // class element for the parabolic CO2 distribution from pol to pol, maximum CO2 volume at equator
    circulation.BC_CO2(Vegetation, h, t, p_dyn, co2);
    co2_cretaceous = circulation.out_co2();

    // class element for the initial conditions for u-v-w-velocity components
    circulation.IC_CellStructure(im_tropopause, u, v, w);

    // class element for the surface temperature computation by radiation flux density
    if (RadiationModel == 3) {
        circulation.BC_Radiation_multi_layer(n, t_j, albedo, epsilon, precipitable_water, Ik, Q_Radiation, Radiation_Balance, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, temp_rad, Q_latent, Q_sensible, Q_bottom, co2_total, p_stat, t, c, h, epsilon_3D, radiation_3D, cloud, ice);
    }

    // class element for the surface temperature computation by radiation flux density
    //  if ( RadiationModel >= 2 ) circulation.BC_Radiation_2D_layer ( im_tropopause, albedo, epsilon, precipitable_water, Ik, Q_Radiation, Radiation_Balance, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, temp_rad, Q_latent, Q_sensible, Q_bottom, co2_total, t, c, h, epsilon_3D, radiation_3D );

    // class element for the parabolic radiation balance distribution from pole to pole, maximum radiation balance amount at equator
    if (RadiationModel == 1) {
        circulation.BC_Radiation_parabolic(Radiation_Balance_par, h);
    }

    //  class element for the initial conditions the latent heat
    //  circulation.Latent_Heat ( rad, the, phi, h, t, tn, u, v, w, p_dyn, p_stat, c, Latency, Q_Sensible, t_cond_3D, t_evap_3D, radiation_3D );

    // class element for the storing of velocity components, pressure and temperature for iteration start
    oldnew.restoreOldNew_3D(.9, u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, aux_p, cn, cloudn, icen, co2n);
    oldnew.restoreOldNew_2D(.9, v, w, p_dyn, aux_p, vn, wn);

    // class element for the computation of the ratio ocean to land areas, also supply and removal of CO2 on land, ocean and by vegetation
    LandArea.land_oceanFraction(h);

    // ******************************************   start of pressure and velocity iterations ************************************************************************

    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of pressure loop : if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::

Pressure_loop:

    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of velocity loop: while ( min >= epsres )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    // min = min_u = min_v = min_w = min_t = min_c = min_p = epsres * 3.;
    min = epsres * 3.;
    int velocity_iter = 0;
    int velocity_iter_2D = 0;

    // query to realize zero divergence of the continuity equation ( div c = 0 )
    while (min >= epsres) {
        // limit of the computation in the sense of time stepscoeff_P
        n++;
        if (n > nm) {
            cout << "       nm = " << nm << "     .....     maximum number of iterations   nm   reached!\n" << endl;
            break;
        }

        // limit of maximum number of iterations ( velocity_iter_max )
        velocity_iter++;

        cout << "  present state of the computation " << endl << "  current time slice, number of iterations, maximum and current number of velocity iterations, maximum and current number of pressure iterations " << endl << endl << "  Ma = " << Ma << "     n = " << n << "    velocity_iter_max = " << velocity_iter_max << "     velocity_iter = " << velocity_iter << "    pressure_iter_max = " << pressure_iter_max << "    pressure_iter = " << pressure_iter << endl;

        if (velocity_iter > velocity_iter_max)
        {
            n--;
            velocity_iter--;
            break;
        }

        if (switch_2D != 1) {

            // **********************************   start of pressure and velocity iterations for the 2D iterational process   *********************************

            // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of pressure loop_2D : if ( pressure_iter_2D > pressure_iter_max_2D )   :::::::::::::::::::::::::::::::::::::::::::

Pressure_loop_2D:

    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of velocity loop_2D: while ( min >= epsres )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            while ( velocity_iter_2D <= velocity_iter_max_2D )
            {

    //      limit of maximum number of iterations ( velocity_iter_max_2D )
                velocity_iter_2D++;
                if ( velocity_iter_2D > velocity_iter_max_2D )
                {
                    goto Pressure_iteration_2D;
                }

    //      class BC_Atmosphaere for the geometry of a shell of a sphere
                boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
                boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

    //      old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
                Accuracy_Atm        min_Residuum_old_2D ( im, jm, km, dthe, dphi );
                min_Residuum_old_2D.residuumQuery_2D ( rad, the, v, w );
                min = min_Residuum_old_2D.out_min (  );

                residuum_old = min;

    //      class RungeKutta for the solution of the differential equations describing the flow properties
                result.solveRungeKutta_2D_Atmosphere ( prepare_2D, rad, the, rhs_v, rhs_w, h, v, w, p_dyn, vn, wn, aux_v, aux_w );

    //      new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
                Accuracy_Atm        min_Residuum_2D ( im, jm, km, dthe, dphi );
                min_Residuum_2D.residuumQuery_2D ( rad, the, v, w );
                min = min_Residuum_2D.out_min (  );
                j_res = min_Residuum_2D.out_j_res (  );
                k_res = min_Residuum_2D.out_k_res (  );

                residuum = min;
                min = fabs ( ( residuum - residuum_old ) / residuum_old );

    //      state of a steady solution resulting from the pressure equation ( min_p ) for pn from the actual solution step
                Accuracy_Atm        min_Stationary_2D ( n, nm, Ma, im, jm, km, min, j_res, k_res, velocity_iter_2D, pressure_iter_2D, velocity_iter_max_2D, pressure_iter_max_2D );
                min_Stationary_2D.steadyQuery_2D ( v, vn, w, wn, p_dyn, aux_p );

                oldnew.restoreOldNew_2D ( 1., v, w, p_dyn, aux_p, vn, wn );


            }



    //  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop_2D: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Pressure_iteration_2D:

    //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
        startPressure.computePressure_2D ( pa, rad, the, p_dyn, h, rhs_v, rhs_w, aux_v, aux_w, aux_p );     // 2D pressure computation causes a pressure jump in radial direction along coast lines, 3D treatment needed later, 2D velocities are though corrected

    //  statements on the convergence und iterational process
        pressure_iter_2D++;
        velocity_iter_2D = 0;

        if ( pressure_iter_2D < pressure_iter_max_2D + 1 ) goto Pressure_loop_2D;

        switch_2D = 1;
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of pressure loop_2D: if ( pressure_iter_2D > pressure_iter_max_2D )   :::::::::::::::::::::::::::::::::::::::::::
    }

        if ( min >= epsres ) {
            time = time + dt;
        }

//      old value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
        Accuracy_Atm        min_Residuum_old ( im, jm, km, dr, dthe, dphi );
        min_Residuum_old.residuumQuery_3D ( rad, the, u, v, w );
        min = min_Residuum_old.out_min (  );

        residuum_old = min;

//      class BC_Atmosphaere for the geometry of a shell of a sphere
        boundary.BC_radius ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
        boundary.BC_theta ( t, u, v, w, p_dyn, c, cloud, ice, co2 );
        boundary.BC_phi ( t, u, v, w, p_dyn, c, cloud, ice, co2 );

//      class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
        LandArea.BC_SolidGround ( RadiationModel, i_max, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, co2_cretaceous, pa, gam, h, u, v, w, t, p_dyn, c, cloud, ice, co2, Vegetation );

//      class RungeKutta for the solution of the differential equations describing the flow properties
        result.solveRungeKutta_3D_Atmosphere ( prepare, n, lv, ls, ep, hp, u_0, t_0, c_0, co2_0, p_0, r_air, r_water, r_water_vapour, r_co2, L_atm, cp_l, R_Air, R_WaterVapour, R_co2, rad, the, phi, rhs_t, rhs_u, rhs_v, rhs_w, rhs_c, rhs_cloud, rhs_ice, rhs_co2, h, t, u, v, w, p_dyn, p_stat, c, cloud, ice, co2, tn, un, vn, wn, cn, cloudn, icen, co2n, aux_u, aux_v, aux_w, Latency, t_cond_3D, t_evap_3D, IceLayer, BuoyancyForce, Q_Sensible, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s );

//      class BC_Bathymetrie for the topography and bathymetry as boundary conditions for the structures of the continents and the ocean ground
        LandArea.BC_SolidGround ( RadiationModel, i_max, g, hp, ep, r_air, R_Air, t_0, t_land, t_cretaceous, t_equator, t_pole, t_tropopause, c_land, c_tropopause, co2_0, co2_equator, co2_pole, co2_tropopause, co2_cretaceous, pa, gam, h, u, v, w, t, p_dyn, c, cloud,ice, co2, Vegetation );

// class element for the surface temperature computation by radiation flux density
        if ( RadiationModel == 3 )                  circulation.BC_Radiation_multi_layer ( n, t_j, albedo, epsilon, precipitable_water, Ik, Q_Radiation, Radiation_Balance, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, temp_rad, Q_latent, Q_sensible, Q_bottom, co2_total, p_stat, t, c, h, epsilon_3D, radiation_3D, cloud, ice );

// class element for the surface temperature computation by radiation flux density
//      if ( RadiationModel >= 2 )                  circulation.BC_Radiation_2D_layer ( im_tropopause, albedo, epsilon, precipitable_water, Ik, Q_Radiation, Radiation_Balance, Radiation_Balance_atm, Radiation_Balance_bot, temp_eff_atm, temp_eff_bot, temp_rad, Q_latent, Q_sensible, Q_bottom, co2_total, t, c, h, epsilon_3D, radiation_3D );

//      new value of the residuum ( div c = 0 ) for the computation of the continuity equation ( min )
        Accuracy_Atm        min_Residuum ( im, jm, km, dr, dthe, dphi );
        min_Residuum.residuumQuery_3D ( rad, the, u, v, w );
        min = min_Residuum.out_min (  );
        int i_res = min_Residuum.out_i_res (  );
        j_res = min_Residuum.out_j_res (  );
        k_res = min_Residuum.out_k_res (  );

        residuum = min;
        min = fabs ( ( residuum - residuum_old ) / residuum_old );

//      statements on the convergence und iterational process
        Accuracy_Atm        min_Stationary ( n, nm, Ma, im, jm, km, min, i_res, j_res, k_res, velocity_iter, pressure_iter, velocity_iter_max, pressure_iter_max, L_atm );
        min_Stationary.steadyQuery_3D ( u, un, v, vn, w, wn, t, tn, c, cn, cloud, cloudn, ice, icen, co2, co2n, p_dyn, aux_p );

// 3D_fields

//      searching of maximum and minimum values of temperature
        string str_max_temperature = " max 3D temperature ", str_min_temperature = " min 3D temperature ", str_unit_temperature = "C";
        MinMax_Atm      minmaxTemperature ( im, jm, km );
        minmaxTemperature.searchMinMax_3D ( str_max_temperature, str_min_temperature, str_unit_temperature, t, h );

//      searching of maximum and minimum values of dynamic pressure
        string str_max_pressure = " max 3D pressure dynamic ", str_min_pressure = " min 3D pressure dynamic ", str_unit_pressure = "hPa";
        MinMax_Atm      minmaxPressure ( im, jm, km );
        minmaxPressure.searchMinMax_3D ( str_max_pressure, str_min_pressure, str_unit_pressure, p_dyn, h );

//      searching of maximum and minimum values of static pressure
        string str_max_pressure_stat = " max 3D pressure static ", str_min_pressure_stat = " min 3D pressure static ", str_unit_pressure_stat = "hPa";
        MinMax_Atm      minmaxPressure_stat ( im, jm, km );
        minmaxPressure_stat.searchMinMax_3D ( str_max_pressure_stat, str_min_pressure_stat, str_unit_pressure_stat, p_stat, h );

        cout << endl << " energies in the three dimensional space: " << endl << endl;

//      searching of maximum and minimum values of radiation_3D
        string str_max_radiation_3D = " max 3D radiation ", str_min_radiation_3D = " min 3D radiation ", str_unit_radiation_3D = "W/m2";
        MinMax_Atm      minmaxLatency ( im, jm, km );
        minmaxLatency.searchMinMax_3D ( str_max_radiation_3D, str_min_radiation_3D, str_unit_radiation_3D, radiation_3D, h );

//      searching of maximum and minimum values of sensible heat
        string str_max_Q_Sensible = " max 3D sensible heat ", str_min_Q_Sensible = " min 3D sensible heat ", str_unit_Q_Sensible = "W/m2";
        MinMax_Atm      minmaxQ_Sensible ( im, jm, km );
        minmaxQ_Sensible.searchMinMax_3D ( str_max_Q_Sensible, str_min_Q_Sensible, str_unit_Q_Sensible, Q_Sensible, h );

//      searching of maximum and minimum values of latency
        string str_max_latency = " max 3D latent heat ", str_min_latency = " min 3D latent heat ", str_unit_latency = "W/m2";
        MinMax_Atm      minmaxRadiation ( im, jm, km );
        minmaxRadiation.searchMinMax_3D ( str_max_latency, str_min_latency, str_unit_latency, Latency, h );

//      searching of maximum and minimum values of t_cond_3D
        string str_max_t_cond_3D = " max 3D condensation temp ", str_min_t_cond_3D = " min 3D condensation temp ", str_unit_t_cond_3D = "C";
        MinMax_Atm      minmaxt_cond_3D ( im, jm, km );
        minmaxt_cond_3D.searchMinMax_3D ( str_max_t_cond_3D, str_min_t_cond_3D, str_unit_t_cond_3D, t_cond_3D, h );

//      searching of maximum and minimum values of t_evap_3D
        string str_max_t_evap_3D = " max 3D evaporation temp ", str_min_t_evap_3D = " min 3D evaporation temp ", str_unit_t_evap_3D = "C";
        MinMax_Atm      minmaxt_evap_3D ( im, jm, km );
        minmaxt_evap_3D.searchMinMax_3D ( str_max_t_evap_3D, str_min_t_evap_3D, str_unit_t_evap_3D, t_evap_3D, h );

        cout << endl << " greenhouse gases: " << endl << endl;

//      searching of maximum and minimum values of water vapour
        string str_max_water_vapour = " max 3D water vapour ", str_min_water_vapour = " min 3D water vapour ", str_unit_water_vapour = "g/kg";
        MinMax_Atm      minmaxWaterVapour ( im, jm, km );
        minmaxWaterVapour.searchMinMax_3D ( str_max_water_vapour, str_min_water_vapour, str_unit_water_vapour, c, h );

//      searching of maximum and minimum values of cloud water
        string str_max_cloud_water = " max 3D cloud water ", str_min_cloud_water = " min 3D cloud water ", str_unit_cloud_water = "g/kg";
        MinMax_Atm      minmaxCloudWater ( im, jm, km );
        minmaxCloudWater.searchMinMax_3D ( str_max_cloud_water, str_min_cloud_water, str_unit_cloud_water, cloud, h );

//      searching of maximum and minimum values of cloud ice
        string str_max_cloud_ice = " max 3D cloud ice ", str_min_cloud_ice = " min 3D cloud ice ", str_unit_cloud_ice = "g/kg";
        MinMax_Atm      minmaxCloudIce ( im, jm, km );
        minmaxCloudIce.searchMinMax_3D ( str_max_cloud_ice, str_min_cloud_ice, str_unit_cloud_ice, ice, h );

//      searching of maximum and minimum values of rain precipitation
        string str_max_P_rain = " max 3D rain ", str_min_P_rain = " min 3D rain ", str_unit_P_rain = "g/kg";
        MinMax_Atm      minmaxPRain ( im, jm, km );
        minmaxPRain.searchMinMax_3D ( str_max_P_rain, str_min_P_rain, str_unit_P_rain, P_rain, h );

//      searching of maximum and minimum values of snow precipitation
        string str_max_P_snow = " max 3D snow ", str_min_P_snow = " min 3D snow ", str_unit_P_snow = "g/kg";
        MinMax_Atm      minmaxPSnow ( im, jm, km );
        minmaxPSnow.searchMinMax_3D ( str_max_P_snow, str_min_P_snow, str_unit_P_snow, P_snow, h );

//      searching of maximum and minimum values of co2
        string str_max_co2 = " max 3D co2 ", str_min_co2 = " min 3D co2 ", str_unit_co2 = "ppm";
        MinMax_Atm      minmaxCO2 ( im, jm, km );
        minmaxCO2.searchMinMax_3D ( str_max_co2, str_min_co2, str_unit_co2, co2, h );

//      searching of maximum and minimum values of epsilon
        string str_max_epsilon = " max 3D epsilon ", str_min_epsilon = " min 3D epsilon ", str_unit_epsilon = "%";
        MinMax_Atm      minmaxEpsilon_3D ( im, jm, km );
        minmaxEpsilon_3D.searchMinMax_3D ( str_max_epsilon, str_min_epsilon, str_unit_epsilon, epsilon_3D, h );

//      searching of maximum and minimum values of buoyancy force
        string str_max_buoyancy_force = " max 3D buoyancy force ", str_min_buoyancy_force = " min 3D buoyancy force ", str_unit_buoyancy_force = "N/m2";
        MinMax_Atm      minmaxBuoyancyForce ( im, jm, km );
        minmaxBuoyancyForce.searchMinMax_3D ( str_max_buoyancy_force, str_min_buoyancy_force, str_unit_buoyancy_force, BuoyancyForce, h );



// 2D-fields

//      searching of maximum and minimum values of co2 total
        cout << endl << " printout of maximum and minimum values of properties at their locations: latitude, longitude" << endl << " results based on two dimensional considerations of the problem" << endl;

        cout << endl << " co2 distribution columnwise: " << endl << endl;

        string str_max_co2_total = " max co2_total ", str_min_co2_total = " min co2_total ", str_unit_co2_total = " / ";
        MinMax_Atm      minmaxCO2_total ( jm, km, coeff_mmWS );
        minmaxCO2_total.searchMinMax_2D ( str_max_co2_total, str_min_co2_total, str_unit_co2_total, co2_total, h );
//      max_CO2_total = minmaxCO2_total.out_maxValue (  );

        cout << endl << " precipitation: " << endl << endl;

//      searching of maximum and minimum values of precipitation
        string str_max_precipitation = " max precipitation ", str_min_precipitation = " min precipitation ", str_unit_precipitation = "mm";
        MinMax_Atm      minmaxPrecipitation ( jm, km, coeff_mmWS );
        minmaxPrecipitation.searchMinMax_2D ( str_max_precipitation, str_min_precipitation, str_unit_precipitation, Precipitation, h );
        double max_Precipitation = minmaxPrecipitation.out_maxValue (  );

/*
//      searching of maximum and minimum values of NASA precipitation
        if ( Ma == 0 )
        {
            string str_max_precipitation = " max precipitation_NASA ", str_min_precipitation = " min precipitation_NASA ", str_unit_precipitation = "mm";
            MinMax      minmaxPrecipitation_NASA ( jm, km, coeff_mmWS );
            minmaxPrecipitation_NASA.searchMinMax_2D ( str_max_precipitation, str_min_precipitation, str_unit_precipitation, precipitation_NASA, h );
        }
*/
//      searching of maximum and minimum values of precipitable water
        string str_max_precipitable_water = " max precipitable water ", str_min_precipitable_water = " min precipitable water ", str_unit_precipitable_water = "mm";
        MinMax_Atm      minmaxPrecipitable_water ( jm, km, coeff_mmWS );
        minmaxPrecipitable_water.searchMinMax_2D ( str_max_precipitable_water, str_min_precipitable_water, str_unit_precipitable_water, precipitable_water, h );
//      max_precipitable_water = minmaxPrecipitable_water.out_minValue (  );

        cout << endl << " energies at see level without convection influence: " << endl << endl;

//      searching of maximum and minimum values of radiation balance
//      string str_max_Radiation_Balance = " max radiation balance ", str_min_Radiation_Balance = " min radiation balance ", str_unit_Radiation_Balance = "W/m2";
//      MinMax      minmaxRadiation_Balance ( jm, km, coeff_mmWS );
//      minmaxRadiation_Balance.searchMinMax_2D ( str_max_Radiation_Balance, str_min_Radiation_Balance, str_unit_Radiation_Balance, Radiation_Balance, h );
//      min_Radiation_Balance = minmaxRadiation_Balance.out_minValue (  );

//      searching of maximum and minimum values of radiation
        string str_max_Q_Radiation = " max 2D Q radiation ", str_min_Q_Radiation = " min 2D Q radiation ", str_unit_Q_Radiation = "W/m2";
        MinMax_Atm      minmaxQ_Radiation ( jm, km, coeff_mmWS );
        minmaxQ_Radiation.searchMinMax_2D ( str_max_Q_Radiation, str_min_Q_Radiation, str_unit_Q_Radiation, Q_Radiation, h );
//      min_Q_Radiation = minmaxQ_Radiation.out_minValue (  );

//      searching of maximum and minimum values of latent energy
        string str_max_Q_latent = " max 2D Q latent ", str_min_Q_latent = " min 2D Q latent ", str_unit_Q_latent = "W/m2";
        MinMax_Atm      minmaxQ_latent ( jm, km, coeff_mmWS );
        minmaxQ_latent.searchMinMax_2D ( str_max_Q_latent, str_min_Q_latent, str_unit_Q_latent, Q_latent, h );

//      searching of maximum and minimum values of sensible energy
        string str_max_Q_sensible = " max 2D Q sensible ", str_min_Q_sensible = " min 2D Q sensible ", str_unit_Q_sensible = "W/m2";
        MinMax_Atm      minmaxQ_sensible ( jm, km, coeff_mmWS );
        minmaxQ_sensible.searchMinMax_2D ( str_max_Q_sensible, str_min_Q_sensible, str_unit_Q_sensible, Q_sensible, h );

//      searching of maximum and minimum values of bottom heat

        string str_max_Q_bottom = " max 2D Q bottom ", str_min_Q_bottom = " min 2D Q bottom heat ", str_unit_Q_bottom = "W/m2";
        MinMax_Atm      minmaxQ_bottom ( jm, km, coeff_mmWS );
        minmaxQ_bottom.searchMinMax_2D ( str_max_Q_bottom, str_min_Q_bottom, str_unit_Q_bottom, Q_bottom, h );


//      searching of maximum and minimum values of latent heat
        string str_max_LatentHeat = " max 2D latent heat ", str_min_LatentHeat = " min 2D latent heat ", str_unit_LatentHeat = "W/m2";
        MinMax_Atm      minmaxLatentHeat_2D ( jm, km, coeff_mmWS );
        minmaxLatentHeat_2D.searchMinMax_2D ( str_max_LatentHeat, str_min_LatentHeat, str_unit_LatentHeat, LatentHeat, h );

//      searching of maximum and minimum values of t_cond_2D
        string str_max_t_cond = " max 2D Condensation ", str_min_t_cond = " min 2D Condensation ", str_unit_t_cond = "W/m2";
        MinMax_Atm      minmaxt_cond_2D ( jm, km, coeff_mmWS );
        minmaxt_cond_2D.searchMinMax_2D ( str_max_t_cond, str_min_t_cond, str_unit_t_cond, Condensation, h );


//      searching of maximum and minimum values of t_evap_2D
        string str_max_t_evap = " max 2D Evaporation ", str_min_t_evap = " min 2D Evaporation ", str_unit_t_evap = "W/m2";
        MinMax_Atm      minmaxt_evap_2D ( jm, km, coeff_mmWS );
        minmaxt_evap_2D.searchMinMax_2D ( str_max_t_evap, str_min_t_evap, str_unit_t_evap, Evaporation, h );

        cout << endl << " secondary data: " << endl << endl;


/*
//      searching of maximum and minimum values of Evaporation
        string str_max_heat_t_Evaporation = " max heat Evaporation ", str_min_heat_t_Evaporation = " min heat Evaporation ", str_unit_heat_t_Evaporation = " W/m2";
        MinMax      minmaxQ_t_Evaporation ( jm, km, coeff_mmWS );
        minmaxQ_t_Evaporation.searchMinMax_2D ( str_max_heat_t_Evaporation, str_min_heat_t_Evaporation, str_unit_heat_t_Evaporation, Q_Evaporation, h );
*/
/*
//      searching of maximum and minimum values of Evaporation by Haude
        string str_max_t_Evaporation_Haude = " max Evaporation Haude ", str_min_t_Evaporation_Haude = " min Evaporation Haude ", str_unit_t_Evaporation_Haude = "mm/d";
        MinMax      minmaxt_Evaporation_Haude ( jm, km, coeff_mmWS );
        minmaxt_Evaporation_Haude.searchMinMax_2D ( str_max_t_Evaporation_Haude, str_min_t_Evaporation_Haude, str_unit_t_Evaporation_Haude, Evaporation_Haude, h );
*/
//      searching of maximum and minimum values of Evaporation by Penman
        string str_max_t_Evaporation_Penman = " max Evaporation Penman ", str_min_t_Evaporation_Penman = " min Evaporation Penman ", str_unit_t_Evaporation_Penman = "mm/d";
        MinMax_Atm      minmaxt_Evaporation_Penman ( jm, km, coeff_mmWS );
        minmaxt_Evaporation_Penman.searchMinMax_2D ( str_max_t_Evaporation_Penman, str_min_t_Evaporation_Penman, str_unit_t_Evaporation_Penman, Evaporation_Penman, h );

        cout << endl << " properties of the atmosphere at the surface: " << endl << endl;

//      searching of maximum and minimum values of albedo
        string str_max_albedo = " max 2D albedo ", str_min_albedo = " min 2D albedo ", str_unit_albedo = "%";
        MinMax_Atm      minmaxAlbedo ( jm, km, coeff_mmWS );
        minmaxAlbedo.searchMinMax_2D ( str_max_albedo, str_min_albedo, str_unit_albedo, albedo, h );
/*
//      searching of maximum and minimum values of epsilon
        string str_max_epsilon = " max 2D epsilon ", str_min_epsilon = " min 2D epsilon ", str_unit_epsilon = "%";
        MinMax      minmaxEpsilon ( jm, km, coeff_mmWS );
        minmaxEpsilon.searchMinMax_2D ( str_max_epsilon, str_min_epsilon, str_unit_epsilon, epsilon, h );
*/

//      computation of vegetation areas
        LandArea.vegetationDistribution ( max_Precipitation, Precipitation, Vegetation, t, h );


//  class element for the initial conditions the latent heat
        circulation.Latent_Heat ( rad, the, phi, h, t, tn, u, v, w, p_dyn, p_stat, c, Latency, Q_Sensible, t_cond_3D, t_evap_3D, radiation_3D );


//      composition of results
        calculate_MSL.run_MSL_data ( n, velocity_iter_max, RadiationModel, rad, the, phi, h, c, cn, co2, co2n, t, tn, p_dyn, p_stat, BuoyancyForce, u, v, w, Latency, Q_Sensible, radiation_3D, t_cond_3D, t_evap_3D, cloud, cloudn, ice, icen, P_rain, P_snow, aux_u, aux_v, aux_w, precipitation_NASA, Evaporation, Condensation, LatentHeat, precipitable_water, Q_Radiation, Q_Evaporation, Q_latent, Q_sensible, Q_bottom, Evaporation_Penman, Evaporation_Haude, Vegetation, Radiation_Balance, Radiation_Balance_par, Radiation_Balance_bot, albedo, co2_total, Precipitation, S_v, S_c, S_i, S_r, S_s );



// printout of results at certain positions
//      calculate_MSL.show_MSL_data ( h, c, t, p_dyn, u, Latency, Q_Sensible, t_cond_3D, t_evap_3D, precipitation_NASA, Evaporation, Condensation, precipitable_water, Q_Radiation, Q_Evaporation, Q_latent, Q_sensible, Q_bottom, Evaporation_Penman, Evaporation_Haude );

//  restoring the velocity component and the temperature for the new time step
        oldnew.restoreOldNew_3D ( 1., u, v, w, t, p_dyn, c, cloud, ice, co2, un, vn, wn, tn, aux_p, cn, cloudn, icen, co2n );

    }


//  ::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of loop: while ( min >= epsres )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



    //  pressure from the Euler equation ( 2. order derivatives of the pressure by adding the Poisson right hand sides )
    startPressure.computePressure_3D ( pa, rad, the, p_dyn, h, rhs_u, rhs_v, rhs_w, aux_u, aux_v, aux_w, aux_p );

    //  statements on the convergence und iterational process
    pressure_iter++;
    velocity_iter = 0;

    if ( pressure_iter < pressure_iter_max + 1 ) {
        velocity_iter_2D = velocity_iter_max_2D;
        pressure_iter_2D = pressure_iter_max_2D;
        switch_2D = 1;
        goto Pressure_loop;
    }

    //  printout in ParaView files, netCDF files and sequel files
    int pressure_iter_aux = pressure_iter - 1;

    //  results written in netCDF format
    // class File_NetCDF to write results in the format of a netCDF-file
    // File_NetCDF printoutNetCDF(im, jm, km);
    // printoutNetCDF.out_NetCDF(output_path, Name_netCDF_File, v, w, h, Precipitation, precipitable_water);

    //  class PostProcess_Atmosphaere for the printing of results
    PostProcess_Atmosphere write_File(im, jm, km, output_path);

    //  writing of data in ParaView files
    //  radial data along constant hight above ground
    int i_radial = 0;
    write_File.paraview_vtk_radial ( bathymetry_name, i_radial, pressure_iter_aux, u_0, t_0, p_0, r_air, c_0, co2_0, radiation_equator, h, p_dyn, p_stat, t_cond_3D, t_evap_3D , BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain, P_snow, Evaporation, Condensation, precipitable_water, Q_bottom, Radiation_Balance, Q_Radiation, Q_latent, Q_sensible, Evaporation_Penman, Evaporation_Haude, Q_Evaporation, precipitation_NASA, Vegetation, albedo, epsilon, Precipitation );

    //  londitudinal data along constant latitudes
    int j_longal = 75;
    write_File.paraview_vtk_longal ( bathymetry_name, j_longal, pressure_iter_aux, u_0, t_0, p_0, r_air, c_0, co2_0, radiation_equator, h, p_dyn, p_stat, t_cond_3D, t_evap_3D, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain, P_snow );

    int k_zonal = 145;
    write_File.paraview_vtk_zonal ( bathymetry_name, k_zonal, pressure_iter_aux, u_0, t_0, p_0, r_air, c_0, co2_0, radiation_equator, h, p_dyn, p_stat, t_cond_3D, t_evap_3D, BuoyancyForce, t, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, radiation_3D, epsilon_3D, P_rain, P_snow, S_v, S_c, S_i, S_r, S_s );

    //  3-dimensional data in cartesian coordinate system for a streamline pattern in panorama view
    //  write_File.paraview_panorama_vts ( bathymetry_name, pressure_iter_aux, u_0, t_0, p_0, r_air, c_0, co2_0, h, t, p_dyn, p_stat, BuoyancyForce, u, v, w, c, co2, cloud, ice, aux_u, aux_v, aux_w, Latency, Q_Sensible, IceLayer, epsilon_3D, P_rain );

    //  3-dimensional data in spherical coordinate system for a streamline pattern in a shell of a sphere
    //  write_File.paraview_vts ( bathymetry_name, n, rad, the, phi, h, t, p_dyn, u, v, w, c, co2, aux_u, aux_v, aux_w, Latency, Rain, Ice, Rain_super, IceLayer );

    //  writing of v-w-data in the v_w_transfer file
    PostProcess_Atmosphere ppa(im, jm, km, output_path);
    ppa.Atmosphere_v_w_Transfer(bathymetry_name, v, w, p_dyn);
    ppa.Atmosphere_PlotData(bathymetry_name, u_0, t_0, v, w, t, c, Precipitation, precipitable_water);

    //  statements on the convergence und iterational process
    velocity_iter = 0;
    n--;

    if ( pressure_iter <= pressure_iter_max ) goto Pressure_loop;

    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of pressure loop: if ( pressure_iter > pressure_iter_max )   :::::::::::::::::::::::::::::::::::::::::::

    //   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    //  final remarks
    cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;

    if ( velocity_iter == velocity_iter_max )   cout << "***** number of time steps      n = " << n << ", end of program reached because of limit of maximum time steps ***** \n\n" << endl;

    if ( min <= epsres )        cout << "***** steady solution reached! *****" << endl;
}

void cAtmosphereModel::Run() {
    // create the output dir
    mkdir(output_path.c_str(), 0777);

    cout << "Output is being written to " << output_path << "\n";

    // write out the config for reproducibility
    // disabled for now
    // std::stringstream output_config_path;
    // output_config_path << output_path << "/config.xml";
    // WriteConfig(output_config_path.str().c_str());


    if (verbose) {
        cout << endl << endl << endl;
        cout << "***** Atmosphere General Circulation Model ( AGCM ) applied to laminar flow" << endl;
        cout << "***** program for the computation of geo-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 4 additional transport equations to describe the water vapour, cloud water, cloud ice and co2 concentration" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** multi-layer and two-layer radiation model for the computation of the surface temperature" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonaly constant" << endl;
        cout << "***** water vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code ( German Weather Forecast )" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** co2 concentration appears in the absorptivity of the radiation models" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz ( roger.grundmann@web.de )" << endl << endl;

        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
    }

    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   begin of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    // time slices to be run after actualizing
    int i_time_slice_max = 15;
    int *time_slice = new int [ i_time_slice_max ];     // time slices in Ma

    time_slice [ 0 ] = 0;                                   // Golonka Bathymetry and Topography
    time_slice [ 1 ] = 10;
    time_slice [ 2 ] = 20;
    time_slice [ 3 ] = 30;
    time_slice [ 4 ] = 40;
    time_slice [ 5 ] = 50;
    time_slice [ 6 ] = 60;
    time_slice [ 7 ] = 70;
    time_slice [ 8 ] = 80;
    time_slice [ 9 ] = 90;
    time_slice [ 10 ] = 100;
    time_slice [ 11 ] = 110;
    time_slice [ 12 ] = 120;
    time_slice [ 13 ] = 130;
    time_slice [ 14 ] = 140;

    //  choice of the time slice to be computed
    for (int i_time_slice = 0; i_time_slice < i_time_slice_max; i_time_slice++) {
        int Ma = time_slice[i_time_slice];
        cout << "Ma = " << Ma << "\n";

        RunTimeSlice(Ma);
    }
//   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   end of time slice loop: if ( i_time_slice >= i_time_slice_max )   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//  final remarks
    cout << endl << "***** end of the Atmosphere General Circulation Modell ( AGCM ) *****" << endl << endl;

    cout << endl;
    cout << "***** end of object oriented C++ program for the computation of 3D-atmospheric circulation *****";
    cout << "\n\n\n\n";
}
