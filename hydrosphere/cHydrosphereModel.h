#ifndef CHYDROSPHEREMODEL_H
#define CHYDROSPHEREMODEL_H

#include <string>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "tinyxml2.h"
#include "PythonStream.h"

using namespace std;
using namespace tinyxml2;

class cHydrosphereModel{

public:

    #include "HydrosphereParams.h.inc"

    cHydrosphereModel();
    ~cHydrosphereModel();
 
    cHydrosphereModel(const cHydrosphereModel&) = delete;
    cHydrosphereModel& operator=(const cHydrosphereModel&) = delete;

    static cHydrosphereModel* get_model(){
        if(!m_model){
            m_model = new cHydrosphereModel();
        }
        return m_model;
    }

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr, dt;
    static const double the0, phi0, r0, residuum_ref_hyd;

    double residuum_old, residuum_loop;

    int Ma;

    bool is_final_result = false;
    bool is_print_mode = false;

    bool is_first_time_slice() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        return (m_current_time == m_time_list.begin());
    }

    std::vector<std::vector<int> > i_bathymetry;

    std::map<float,float> m;
    float get_temperatures_from_curve(float time, std::map<float, float>& m) const;

    std::set<float>::const_iterator get_current_time() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }else{
            return m_current_time;
        }
    }

    std::set<float>::const_iterator get_previous_time() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        if(m_current_time != m_time_list.begin()){
            std::set<float>::const_iterator ret = m_current_time;
            ret--;
            return ret;
        }
        else{
            throw("The current time is the only time slice for now. There is no previous time yet.");
        }
    }

    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

private:

    static cHydrosphereModel* m_model;

    PythonStream ps;
    std::streambuf *backup;

    string bathymetry_name;
    string hy = "OGCM";

    bool has_printed_welcome_msg;
    bool is_global_temperature_curve_loaded(){
        return !m_global_temperature_curve.empty();
    }
    bool is_equat_temperature_curve_loaded(){
        return !m_equat_temperature_curve.empty();
    }
    bool is_pole_temperature_curve_loaded(){
        return !m_pole_temperature_curve.empty();
    }

    int iter_cnt, iter_cnt_3d;
    int velocity_iter, pressure_iter;
    int velocity_iter_2D, pressure_iter_2D;

    const int c43 = 4.0/3.0, c13 = 1.0/3.0;

    static const int im = 41, jm = 181, km = 361;

    int iter_n = 0;

    double t_paleo_total = 0.0;
    double t_pole_total = 0.0;
    double t_global_mean = 0.0;

    std::set<float> m_time_list;
    std::set<float>::const_iterator m_current_time;

    std::map<float,float> m_global_temperature_curve;
    std::map<float,float> m_equat_temperature_curve;
    std::map<float,float> m_pole_temperature_curve;

    double find_residuum_hyd();

    void SetDefaultConfig();
    void reset_arrays();
    void write_file(std::string &bathymetry_name, 
        std::string &output_path, bool is_final_result);
    void load_temperature_curve();
    void save_data();
    void run_data_hyd();
    void IC_t_WestEastCoast();
    void IC_u_WestEastCoast();
    void IC_Equatorial_Currents();
    void IC_CircumPolar_Current();
    void PresStat_SaltWaterDens();
    void solveRungeKutta_3D_Hydrosphere();
    void RK_RHS_3D_Hydrosphere(int i, int j, int k);
    void computePressure_3D();
    void store_intermediate_data_3D(float coeff=1);
    void print_welcome_msg();
    void print_final_remarks();
    void print_loop_3D_headings();
    void print_min_max_hyd();
    void BC_radius();
    void BC_theta();
    void BC_phi();
    void BC_SolidGround();
    void init_bathymetry(const string &bathymetry_file);
    void BC_Surface_Salinity_NASA(const string &Name_SurfaceSalinity_File);
    void BC_Surface_Temperature_NASA
       (const string &Name_SurfaceTemperature_File);
    void init_EkmanSpiral();
    void init_EkmanSpiral_surface();
    void init_temperature(int Ma);
    void init_salinity();
    void SalinityEvaporation();
    void Pressure_Limitation_Hyd();
    void ValueLimitationHyd();
    void paraview_vts(const string &Name_Bathymetry_File, int n);
    void paraview_panorama_vts(const string &Name_Bathymetry_File, int n);
    void paraview_vtk_zonal(const string &Name_Bathymetry_File, int k_zonal, int n);
    void paraview_vtk_radial(const string &Name_Bathymetry_File, int i_radial, int n);
    void paraview_vtk_longal(const string &Name_Bathymetry_File, int j_longal, int n);
    void HydrospherePlotData(const string &Name_Bathymetry_File, int iter_cnt);
    void HydrosphereDataTransfer(const string &Name_Bathymetry_File);
    void run_2D_loop();
    void run_3D_loop();
    void searchMinMax_2D(string, string, 
        string, Array_2D &, double coeff=1.);
    void searchMinMax_3D(string, string, 
        string, Array &, double coeff=1.,
        std::function< double(double) > lambda = [](double i)->double{return i;},
        bool print_heading=false);
    void load_global_temperature_curve();
    void load_equat_temperature_curve();
    void load_pole_temperature_curve();
    void read_Hydrosphere_Surface_Data(int Ma);

    Array_1D rad; // radial coordinate direction
    Array_1D the; // lateral coordinate direction
    Array_1D phi; // longitudinal coordinate direction
    Array_1D aux_grad_v; // auxilliar array
    Array_1D aux_grad_w; // auxilliar array
    
    Array_2D Bathymetry; // Bathymetry in m
    Array_2D value_top; // auxiliar field for bathymetzry
    Array_2D Upwelling; // upwelling
    Array_2D Downwelling; // downwelling
    Array_2D EkmanPumping; // 2D bottom water summed up in a vertical column
    Array_2D BuoyancyForce_2D; // radiation balance at the surface
    Array_2D salinity_evaporation; // additional salinity by evaporation
    Array_2D Evaporation_Dalton; // evaporation by Dalton in [mm/d]
    Array_2D Evaporation_Penman; // evaporation by Penman in [mm/d]
    Array_2D Precipitation; // areas of higher precipitation
    Array_2D precipitation_NASA; // surface precipitation from NASA
    Array_2D temperature_NASA; // surface temperature from NASA
    Array_2D temp_landscape; // landscape temperature
    Array_2D temp_reconst; // surface temperature from reconstuction tool
    Array_2D c_fix; // local surface salinity fixed for iterations
    Array_2D v_wind; // v-component of surface wind
    Array_2D w_wind; // w-component of surface wind

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
    Array p_hydro; // static pressure

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
    Array BuoyancyForce; // 3D buoyancy force
    Array CoriolisForce; // Coriolis force terms
    Array PressureGradientForce; // Force caused by normal pressure gradient
};
#endif
