#ifndef CATMOSPHEREMODEL_H
#define CATMOSPHEREMODEL_H

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

class cAtmosphereModel{

public:

    #include "AtmosphereParams.h.inc"

    cAtmosphereModel();
    ~cAtmosphereModel();

    cAtmosphereModel(const cAtmosphereModel&) = delete;
    cAtmosphereModel& operator=(const cAtmosphereModel&) = delete;

    static cAtmosphereModel* get_model(){
        if(!m_model){
            m_model = new cAtmosphereModel();
        }
        return m_model;
    }

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr, dt;
    static const double the0, phi0, r0;

    double residuum_old, residuum_loop;
    double tke_surf, dis_surf, tke_source_surf, dis_source_surf, nue_surf, prod_surf;

    int Ma;

    bool is_final_result = false;
    bool is_print_mode = false;

    bool is_first_time_slice() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        return (m_current_time == m_time_list.begin());
    }

    std::vector<std::vector<int> > i_topography;

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

    /*
     * Given a latitude, return the layer index of tropopause
    */
    int get_tropopause_layer(int j){
        assert(j>=0);
        assert(j<jm);
        //refer to  BC_Thermo::TropopauseLocation and BC_Thermo::GetTropopauseHightAdd
        //tropopause height is proportional to the mean tropospheric temperature.
        //higher near the equator - warm troposphere
        //lower at the poles - cold troposphere
        return tropopause_layers[j];
    }
    /*
     *
    */
    int get_surface_layer(int j, int k){
        return i_topography[j][k];
    }
    /*
     * This function must be called after init_layer_heights()
     * Given a layer index i, return the height of this layer
    */
    float get_layer_height(int i){
        if(0>i || i>im-1){
            return -1;
        }
        return m_layer_heights[i];
    }
    std::vector<float> get_layer_heights(){
        return m_layer_heights;
    }   
    /*
    * Given a altitude, return the layer index
    */
    int get_layer_index(float height){
        std::size_t i = 0;
        for(; i<m_layer_heights.size(); i++){
            if(height<m_layer_heights[i])
                return i-1;
        }
        return i;
    }

    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

private:

    static cAtmosphereModel* m_model;

    PythonStream ps;
    std::streambuf *backup;

    string at = "AGCM";
    string bathymetry_name;

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

    double t_paleo_total = 0.0;
    double t_pole_total = 0.0;
    double t_global_mean = 0.0;

    std::set<float> m_time_list;
    std::set<float>::const_iterator m_current_time;

    std::map<float,float> m_global_temperature_curve;
    std::map<float,float> m_equat_temperature_curve;
    std::map<float,float> m_pole_temperature_curve;

    std::vector<std::vector<double> > M_u_Base;
    std::vector<std::vector<double> > M_d_LFS;
    std::vector<double> step;
    std::vector<double> short_wave_radiation; // lateral short wave radiation
    std::vector<double> radiation_original; // original radiation
    std::vector<double> tropopause_layers; // keep the tropopause layer index
//    std::vector<double> temp_tropopause; // lateral temperature distribution along the tropopause
    std::vector<double> cloud_loc; // lateral cloudwater distribution
    std::vector<double> cloud_max; // lateral cloud_max distribution
    std::vector<double> ice_loc; // lateral ice distribution
    std::vector<double> ice_max; // lateral ice_max distribution
    std::vector<double> c_land_red;
    std::vector<double> c_ocean_red;
    std::vector<double> CAPE;
    std::vector<double> K_u;
    std::vector<double> K_d;
    std::vector<double> lapse_rate;
    std::vector<float> m_layer_heights;

    double find_residuum_atm();

    void SetDefaultConfig();
    void reset_arrays();
    void print_min_max_atm();
    void write_file(std::string &bathymetry_name, 
        std::string &output_path, bool is_final_result);
    void run_3D_loop();
    void load_global_temperature_curve();
    void load_equat_temperature_curve();
    void load_pole_temperature_curve();
    void calculate_node_weights();
    void init_steps();
    void init_tropopause_layers();
    void RK_RHS_3D_Atmosphere(int i, int j, int k);
    void solveRungeKutta_3D_Atmosphere();
    void fft(Array &);
    void ValueLimitationAtm();
    void BC_SolidGround(); 
    void init_co2(int Ma);
    void init_temperature(int Ma);
    void init_temperature_adiabatic(int Ma);
    void init_water_vapour();
    void init_velocities();
    void init_layer_heights(){
        double h = L_atm/(double)(im-1);
        for(int i=0; i<im; i++){
            if(use_stretched_coordinate_system)
                m_layer_heights.push_back((exp(zeta * (rad.z[i] - 1.0)) - 1.0) * h); // stretched coordinate system
            else
                m_layer_heights.push_back((rad.z[i] - 1.0) * L_atm);  // unstretched coordinate system
//            std::cout << i << "     " << h << "     " << rad.z[i] - 1.0 << "     " << m_layer_heights.back() << std::endl;
        } 
        return;
    }
    void init_topography(const string &topo_filename);
    void init_u(Array &u, int j);
    void init_v_or_w(Array &v_or_w, int lat, double coeff_trop, double coeff_sl);
    void init_v_or_w_above_tropopause(Array &v_or_w, int lat, double coeff);
    void init_v_or_w_below_base(Array &v_or_w, int i_base, int lat, double coeff);
    void form_diagonals(Array &a, int start, int end);
    void smooth_transition(Array &u, Array &v, Array &w, int lat);
    void RadiationMultiLayer();
    void save_data();
    void save_array(const string& fn, const Array& a);
    void SaturationAdjustment();
    void ZeroCategoryIceScheme();
    void OneCategoryIceScheme();
    void TwoCategoryIceScheme();
    void ThreeCategoryIceScheme();
    void PressureDensity();
    void LatentSensibleHeat();
    void MassStreamfunction();
    void StandAtm_DewPoint_HumidRel();
    void IC_vwt_WestEastCoast();   
    void IC_t_WestEastCoast();
    void BC_radius();
    void BC_theta();
    void BC_phi();
    void BC_pole();
    void computePressure_3D();
    void print_welcome_msg();
    void print_final_remarks();
    void print_loop_3D_headings();
    void WaterVapourEvaporation();
    void MoistConvectionMidL();
    void MoistConvectionShall();
    void store_intermediate_data_3D(double coeff=1);
    void adjust_temperature_IC(double** t, int jm, int km);
    void check_data();
    void check_data(Array& a, Array&an, const std::string& name);
    void paraview_panorama_vts(string &Name_Bathymetry_File, int n);
    void paraview_vtk_radial(string &Name_Bathymetry_File, int Ma, int i_radial, int n);
    void paraview_vtk_zonal(string &Name_Bathymetry_File, int k_zonal, int n);
    void paraview_vtk_longal(string &Name_Bathymetry_File, int j_longal, int n); 
    void AtmosphereDataTransfer(const string &Name_Bathymetry_File);
    void AtmospherePlotData(const string &Name_Bathymetry_File, int iter_cnt);
    void read_Atmosphere_Surface_Data(int Ma);
    void run_data_atm();
    void searchMinMax_2D(string, string, 
        string, Array_2D &, double coeff=1.0);
    void searchMinMax_3D(string, string, 
        string, Array &, double coeff=1.0,
        std::function< double(double) > lambda = [](double i)->double{return i;},
        bool print_heading=false);

    Array_1D rad; // radial coordinate direction
    Array_1D the; // lateral coordinate direction
    Array_1D phi; // longitudinal coordinate direction
    Array_1D aux_grad_v; // auxilliar array

    Array_2D Topography; // topography
    Array_2D value_top; // auxiliar topography
    Array_2D Vegetation; // vegetation via precipitation
    Array_2D Precipitation; // areas of higher precipitation
    Array_2D precipitable_water; // areas of precipitable water in the air
    Array_2D precipitation_NASA; // surface precipitation from NASA
    Array_2D temperature_NASA; // surface temperature from NASA
    Array_2D velocity_v_NASA; // surface v-velocity from NASA
    Array_2D velocity_w_NASA; // surface w-velocity from NASA
    Array_2D temp_reconst; // surface temperature from reconstuction tool
    Array_2D temp_landscape; // landscape temperature
    Array_2D p_hydro_landscape; // landscape static pressure
    Array_2D r_dry_landscape; // landscape dry air density
    Array_2D r_humid_landscape; // landscape humid air density
    Array_2D albedo; // albedo = reflectivity
    Array_2D epsilon_2D; // epsilon = absorptivity
    Array_2D Q_radiation; // heat from the radiation balance in [W/m2]
    Array_2D Q_Evaporation; // evaporation heat of water by Kuttler
    Array_2D Q_latent; // latent heat from bottom values by the energy transport equation
    Array_2D Q_sensible; // sensible heat from bottom values by the energy transport equation
    Array_2D Q_bottom; // difference by Q_radiation - Q_latent - Q_sensible
    Array_2D vapour_evaporation; // water vapour by evaporation in [mm/d]
    Array_2D Evaporation_Dalton; // evaporation by Dalton in [mm/d]
    Array_2D Evaporation_Penman; // evaporation by Penman in [mm/d]
    Array_2D co2_total; // areas of higher co2 concentration
    Array_2D dew_point_temperature; // dew point temperature
    Array_2D condensation_level; // local condensation level
    Array_2D c_fix; // local surface water vapour fixed for iterations
    Array_2D tropopause_height; // local height of the tropopause

    Array h; // bathymetry, depth from sea level
    Array t; // temperature
    Array u; // u-component velocity component in r-direction
    Array v; // v-component velocity component in theta-direction
    Array w; // w-component velocity component in phi-direction
    Array c; // water vapour
    Array cloud; // cloud water
    Array ice; // cloud ice
    Array gr; // cloud graupel
    Array cloudiness; // cloudiness, N in literature
    Array co2; // CO2

    Array tn; // temperature new
    Array un; // u-velocity component in r-direction new
    Array vn; // v-velocity component in theta-direction new
    Array wn; // w-velocity component in phi-direction new
    Array cn; // water vapour new
    Array cloudn; // cloud water new
    Array icen; // cloud ice new
    Array grn; // cloud ice new
    Array co2n; // CO2 new

    Array p_dyn; // dynamic pressure
    Array stream; // mass stream function
    Array u_stream; // u-velocity by mass stream function
    Array p_hydro; // static pressure
    Array TempStand; // US Standard Atmosphere Temperature
    Array TempDewPoint; // Dew Point Temperature
    Array HumidityRel; // relative humidity

    Array rhs_t; // auxilliar field RHS temperature
    Array rhs_u; // auxilliar field RHS u-velocity component
    Array rhs_v; // auxilliar field RHS v-velocity component
    Array rhs_w; // auxilliar field RHS w-velocity component
    Array rhs_c; // auxilliar field RHS water vapour
    Array rhs_cloud; // auxilliar field RHS cloud water
    Array rhs_ice; // auxilliar field RHS cloud ice
    Array rhs_g; // auxilliar field RHS cloud graupel
    Array rhs_co2; // auxilliar field RHS CO2

    Array aux_u; // auxilliar field u-velocity component
    Array aux_v; // auxilliar field v-velocity component
    Array aux_w; // auxilliar field w-velocity component
    Array aux_t; // auxilliar field t
    Array Q_Latent; // latent heat
    Array Q_Sensible; // sensible heat
    Array BuoyancyForce; // buoyancy force, Boussinesque approximation
    Array CoriolisForce; // Coriolis force terms
    Array PressureGradientForce; // Force caused by normal pressure gradient
    Array epsilon; // emissivity/ absorptivity
    Array radiation; // radiation
    Array P_rain; // rain precipitation mass rate
    Array P_snow; // snow precipitation mass rate
    Array P_graupel; // graupel precipitation mass rate
    Array P_conv_midl; // rain formation by mid-level cloud convection
    Array P_conv_shall; // rain formation by shallow cloud convection
    Array S_v; // water vapour mass rate
    Array S_c; // cloud water mass rate
    Array S_i; // cloud ice mass rate
    Array S_r; // rain mass rate
    Array S_s; // snow mass rate
    Array S_g; // graupel mass rate
    Array S_c_c; // cloud water mass rate due to condensation and evaporation in the saturation adjustment technique
    Array M_u; // moist convection within the updraft
    Array M_d; // moist convection within the downdraft
    Array MC_t; // moist convection acting on dry static energy
    Array MC_q; // moist convection acting on water vapour development
    Array MC_v; // moist convection acting on v-velocity component
    Array MC_w; // moist convection acting on w-velocity component
    Array r_dry; // density of dry air
    Array r_humid; // density of humid air
    Array g_p; // conversion cloud droplets to raindrops
    Array c_u; // condensation in the updraft
    Array e_d; // evaporation of precipitation in the downdraft
    Array e_l; // evaporation of cloud water in the environment
    Array e_p; // evaporation of cloud water in the environment
    Array s; // dry static energy
    Array s_u; // dry static energy in the updraft
    Array s_d; // dry static energy in the downdraft
    Array u_u; // u-velocity component in the updraft
    Array u_d; // u-velocity component in the downdraft
    Array v_u; // u-velocity component in the updraft
    Array v_d; // u-velocity component in the downdraft
    Array w_u; // u-velocity component in the updraft
    Array w_d; // u-velocity component in the downdraft
    Array q_v_u; // water vapour in the updraft
    Array q_v_d; // water vapour in the downdraft
    Array q_c_u; // cloud water in the updraft
    Array E_u; // moist entrainment in the updraft
    Array D_u; // moist detrainment in the updraft
    Array E_d; // moist entrainment in the downdraft
    Array D_d; // moist detrainment in the downdraft
    Array_2D i_Base; // locations of the cloud base
    Array_2D i_LFS; // locations of the cloud top
};
#endif
