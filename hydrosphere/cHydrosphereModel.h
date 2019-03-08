#ifndef CHYDROSPHEREMODEL_H
#define CHYDROSPHEREMODEL_H

#include <string>
#include <vector>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

class cHydrosphereModel{
public:

    cHydrosphereModel();
    ~cHydrosphereModel();

    // FUNCTIONS
    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

    #include "HydrosphereParams.h.inc"

private:
    void SetDefaultConfig();
    void reset_arrays();
    void write_file( std::string &bathymetry_name, string& filepath, bool is_final_result = false);
    void save_data();

    void solveRungeKutta_2D_Hydrosphere();
    void solveRungeKutta_3D_Hydrosphere();

    void RK_RHS_3D_Hydrosphere(int i, int j, int k);
    void RK_RHS_2D_Hydrosphere(int j, int k);

    void store_intermediate_data_2D(float coeff=1);
    void store_intermediate_data_3D(float coeff=1);

    void print_welcome_msg();
    
    void print_min_max();
    
    const int im = 41, jm = 181, km = 361, nm = 200;

    static const float dr, dt, pi180, the_degree, phi_degree, dthe, dphi;

    int iter_cnt, iter_cnt_3d;
    int m_current_time;

    bool has_printed_welcome_msg;
    
    // 1D arrays
    Array_1D rad; // radial coordinate direction
    Array_1D the; // lateral coordinate direction
    Array_1D phi; // longitudinal coordinate direction

    // 2D arrays
    Array_2D Bathymetry; // Bathymetry in m
    Array_2D value_top; // auxiliar field for bathymetzry

    Array_2D Upwelling; // upwelling
    Array_2D Downwelling; // downwelling
    Array_2D BottomWater; // 2D bottom water summed up in a vertical column

    Array_2D SaltFinger;   // salt bulge of higher density
    Array_2D SaltDiffusion; // salt bulge of lower density
    Array_2D Salt_total;   // rate of salt summed up in a vertical column

    Array_2D BuoyancyForce_2D; // radiation balance at the surface

    Array_2D Evaporation_Dalton; // evaporation by Penman in [mm/d]
    Array_2D Precipitation; // areas of higher precipitation

    // 3D arrays
    Array h; // bathymetry, depth from sea level

    Array t; // temperature
    Array u; // u-component velocity component in r-direction
    Array v; // v-component velocity component in theta-direction
    Array w; // w-component velocity component in phi-direction
    Array c; // water vapour

    Array tn; // temperature new
    Array un; // u-velocity component in r-direction new
    Array vn; // v-velocity component in theta-direction new
    Array wn; // w-velocity component in phi-direction new
    Array cn; // water vapour new

    Array p_dyn; // dynamic pressure
    Array p_dynn; // dynamic pressure new
    Array p_stat; // static pressure

    Array rhs_t; // auxilliar field RHS temperature
    Array rhs_u; // auxilliar field RHS u-velocity component
    Array rhs_v; // auxilliar field RHS v-velocity component
    Array rhs_w; // auxilliar field RHS w-velocity component
    Array rhs_c; // auxilliar field RHS water vapour

    Array aux_u; // auxilliar field u-velocity component
    Array aux_v; // auxilliar field v-velocity component
    Array aux_w; // auxilliar field w-velocity component

    Array Salt_Finger; // salt bulge of higher density
    Array Salt_Diffusion; // salt bulge of lowerer density and temperature
    Array Salt_Balance; // +/- salt balance

    Array r_water; // water density as function of pressure
    Array r_salt_water; // salt water density as function of pressure and temperature
    Array BuoyancyForce_3D; // 3D buoyancy force
};
#endif
