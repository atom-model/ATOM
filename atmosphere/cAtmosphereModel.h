#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

#include <string>
#include <set>

#include "Array.h"
#include "Array_2D.h"
#include "Array_1D.h"
#include "tinyxml2.h"
#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;

class BC_Atmosphere;
class RungeKutta_Atmosphere;
class BC_Bathymetry_Atmosphere;
class RHS_Atmosphere;
class Pressure_Atm;
class Restore_Atm;
class Results_MSL_Atm;
class BC_Thermo;

class cAtmosphereModel {
public:
    const char *filename;

    cAtmosphereModel();
    ~cAtmosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);
    void PrintMaxMinValues();

    std::set<int>::const_iterator get_current_time() const{
        return m_current_time;
    }

    std::set<int>::const_iterator get_previous_time() const{
        if(m_current_time != m_time_list.begin()){
            std::set<int>::const_iterator ret = m_current_time;
            ret--;
            return ret;
        }
        else{
            return m_current_time;
        }
    }

    #include "AtmosphereParams.h.inc"

private:
    void SetDefaultConfig();
    void Reset();
    void WriteFile(int n, std::string &bathymetry_name, string& filepath);
    void Run2DLoop(int Ma, int n, int nm, int i_max, int pressure_iter_max_2D, 
                   int velocity_iter_max_2D, BC_Atmosphere &boundary, RungeKutta_Atmosphere &result,
                   BC_Bathymetry_Atmosphere &LandArea, RHS_Atmosphere &prepare_2D,
                   Pressure_Atm &startPressure, Restore_Atm &oldnew);

    void Run3DLoop(int Ma, int n, int nm, int i_max, int pressure_iter_max_3D, 
                   int velocity_iter_max_3D, BC_Atmosphere &boundary, RungeKutta_Atmosphere &result,
                   BC_Bathymetry_Atmosphere &LandArea, RHS_Atmosphere &prepare,
                   Pressure_Atm &startPressure, Restore_Atm &oldnew, Results_MSL_Atm &calculate_MSL, BC_Thermo &circulation);

    PythonStream ps;
    std::streambuf *backup;

    std::set<int> m_time_list;
    std::set<int>::const_iterator m_current_time;

    static const int im=41, jm=181, km=361, nm=200;
    int j_res, k_res;    

    double residuum_min;

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr, dt;
    static const double the0, phi0, r0;

    double coeff_mmWS, d_max_Precipitation;
    int *im_tropopause;  // location of the tropopause

//  class Array for 1-D, 2-D and 3-D field declarations

// 1D arrays
    Array_1D rad; // radial coordinate direction
    Array_1D the; // lateral coordinate direction
    Array_1D phi; // longitudinal coordinate direction

// 2D arrays
    Array_2D Topography; // topography
    Array_2D value_top; // auxiliar topography

    Array_2D Vegetation; // vegetation via precipitation

    Array_2D LatentHeat; // areas of higher latent heat
    Array_2D Condensation; // areas of higher condensation
    Array_2D Evaporation; // areas of higher evaporation

    Array_2D Precipitation; // areas of higher precipitation
    Array_2D precipitable_water; // areas of precipitable water in the air
    Array_2D precipitation_NASA; // surface precipitation from NASA

    Array_2D Ik; // direct sun radiation, short wave
    Array_2D Radiation_Balance; // radiation balance at the surface

    Array_2D temperature_NASA; // surface temperature from NASA

    Array_2D albedo; // albedo = reflectivity
    Array_2D epsilon; // epsilon = absorptivity

    Array_2D Q_Radiation; // heat from the radiation balance in [W/m2]
    Array_2D Q_Evaporation; // evaporation heat of water by Kuttler
    Array_2D Q_latent; // latent heat from bottom values by the energy transport equation
    Array_2D Q_sensible; // sensible heat from bottom values by the energy transport equation
    Array_2D Q_bottom; // difference by Q_Radiation - Q_latent - Q_sensible

    Array_2D Evaporation_Haude; // evaporation by Haude in [mm/d]
    Array_2D Evaporation_Penman; // evaporation by Penman in [mm/d]

    Array_2D co2_total; // areas of higher co2 concentration

// 3D arrays
    Array h; // bathymetry, depth from sea level
    Array t; // temperature
    Array u; // u-component velocity component in r-direction
    Array v; // v-component velocity component in theta-direction
    Array w; // w-component velocity component in phi-direction
    Array c; // water vapour
    Array cloud; // cloud water
    Array ice; // cloud ice
    Array co2; // CO2

    Array tn; // temperature new
    Array un; // u-velocity component in r-direction new
    Array vn; // v-velocity component in theta-direction new
    Array wn; // w-velocity component in phi-direction new
    Array cn; // water vapour new
    Array cloudn; // cloud water new
    Array icen; // cloud ice new
    Array co2n; // CO2 new

    Array p_dyn; // dynamic pressure
    Array p_dynn; // dynamic pressure
    Array p_stat; // static pressure

    Array rhs_t; // auxilliar field RHS temperature
    Array rhs_u; // auxilliar field RHS u-velocity component
    Array rhs_v; // auxilliar field RHS v-velocity component
    Array rhs_w; // auxilliar field RHS w-velocity component
    Array rhs_p; // auxilliar field RHS w-velocity component
    Array rhs_c; // auxilliar field RHS water vapour
    Array rhs_cloud; // auxilliar field RHS cloud water
    Array rhs_ice; // auxilliar field RHS cloud ice
    Array rhs_co2; // auxilliar field RHS CO2

    Array aux_u; // auxilliar field u-velocity component
    Array aux_v; // auxilliar field v-velocity component
    Array aux_w; // auxilliar field w-velocity component

    Array Latency; // latent heat
    Array Q_Sensible; // sensible heat
    Array IceLayer; // ice shield
    Array t_cond_3D; // condensation temperature
    Array t_evap_3D; // evaporation temperature
    Array BuoyancyForce; // buoyancy force, Boussinesque approximation
    Array epsilon_3D; // emissivity/ absorptivity
    Array radiation_3D; // radiation

    Array P_rain; // rain precipitation mass rate
    Array P_snow; // snow precipitation mass rate
    Array S_v; // water vapour mass rate due to category two ice scheme
    Array S_c; // cloud water mass rate due to category two ice scheme
    Array S_i; // cloud ice mass rate due to category two ice scheme
    Array S_r; // rain mass rate due to category two ice scheme
    Array S_s; // snow mass rate due to category two ice scheme
    Array S_c_c; // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique
};

#endif
