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

namespace{
    std::function<double(double)> default_lambda=[](double i)->double{return i;};
}

class cAtmosphereModel{
public:
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

    void LoadConfig(const char *filename);
    void Run();
    void RunTimeSlice(int time_slice);

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

    bool is_first_time_slice() const{
        if(m_time_list.empty()){
            throw("The time list is empty. It is likely the model has not started yet.");
        }
        return (m_current_time == m_time_list.begin());
    }

    #include "AtmosphereParams.h.inc"

    float get_mean_temperature_from_curve(float time) const;
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

    float calculate_mean_temperature(const Array& t);

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr, dt;
    static const double the0, phi0, r0;
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
private:
    void SetDefaultConfig();
    void reset_arrays();
    void print_min_max_atm();
    void write_file(std::string &bathymetry_name, string& filepath, 
        bool is_final_result = false);

    float GetMean_3D(int jm, int km, Array &val_3D);
    float GetMean_2D(int jm, int km, Array_2D &val_2D);

    void run_2D_loop();
    void run_3D_loop();
public:
    void RK_RHS_2D_Atmosphere(int j, int k);
    void RK_RHS_3D_Atmosphere(int i, int j, int k);
    void solveRungeKutta_2D_Atmosphere();
    void solveRungeKutta_3D_Atmosphere();
    void fft(Array &);
private:
    void load_temperature_curve();

    std::map<float,float> m_temperature_curve;

    bool is_temperature_curve_loaded(){
        return !m_temperature_curve.empty();
    }

    void calculate_node_weights();
    void CalculateNodeWeights(int jm, int km);

    float calculate_mean_temperature(){
        return calculate_mean_temperature(t);
    }

    const  int c43 = 4./3., c13 = 1./3.;

    static const int im=41, jm=181, km=361, nm=200;

    std::vector<std::vector<int> > i_topography;
//    std::vector<std::vector<int> > i_Base;
//    std::vector<std::vector<int> > i_LFS;
    std::vector<std::vector<double> > M_u_Base;
    std::vector<std::vector<double> > M_d_LFS;
    std::vector<double> tropopause_layers; // keep the tropopause layer index
    std::vector<double> temp_tropopause; // lateral temperature distribution along the tropopause
    std::vector<double> cloud_loc; // lateral cloudwater distribution
    std::vector<double> cloud_max; // lateral cloud_max distribution
    std::vector<double> ice_loc; // lateral ice distribution
    std::vector<double> ice_max; // lateral ice_max distribution
    std::vector<double> c_land_red;
    std::vector<double> c_ocean_red;

    void init_tropopause_layers();
    void restrain_temperature();
    void PressureLimitationAtm();
    void ValueLimitationAtm();
    void vegetation_distribution();
    void BC_SolidGround(); 
    void init_co2();
    void init_temperature();
    void init_dynamic_pressure();
    void init_water_vapour();
    void init_velocities();
    void init_layer_heights(){
        const float zeta = 3.715;
        float h = L_atm/(im-1);
        for(int i=0; i<im; i++){
            m_layer_heights.push_back((exp(zeta * (rad.z[ i ] - 1.)) - 1) * h);
//            m_layer_heights.push_back((exp(i) - 1.) * h);
//            std::cout << m_layer_heights.back() << std::endl;
        } 
        return;
    }
    void init_topography(const string &topo_filename);
    void init_u(Array &u, int j);
    void init_v_or_w(Array &v_or_w, int lat, double coeff_trop, double coeff_sl);
    void init_v_or_w_above_tropopause(Array &v_or_w, int lat, double coeff);
//    void init_v_or_w(Array &v_or_w, int lat_1, int lat_2, double coeff_trop, double coeff_sl);
//    void init_v_or_w_above_tropopause(Array &v_or_w, int lat_1, int lat_2, double coeff);
    void form_diagonals(Array &a, int start, int end);
    void smooth_transition(Array &u, Array &v, Array &w, int lat);
    void RadiationMultiLayer();
    void save_data();
    void save_array(const string& fn, const Array& a);
    void SaturationAdjustment();
    void TwoCategoryIceScheme();
    void PressureDensity();
    void LatentHeat();
    void MassStreamfunction();
    void USStand_DewPoint_HumidRel();
    void IC_vwt_WestEastCoast();   
    void IC_t_WestEastCoast();
    void read_NASA_temperature(const string &fn);
    void read_NASA_precipitation(const string&);
    void BC_radius();
    void BC_theta();
    void BC_phi();
    void BC_pole();
    void computePressure_3D();
    void computePressure_2D();
    void print_welcome_msg();
    void print_final_remarks();
    void WaterVapourEvaporation();
    void MoistConvection();
    void store_intermediate_data_2D(float coeff=1);
    void store_intermediate_data_3D(float coeff=1);
    void adjust_temperature_IC(double** t, int jm, int km);
    void check_data();
    void check_data(Array& a, Array&an, const std::string& name);
    void paraview_panorama_vts(string &Name_Bathymetry_File, int n);
    void paraview_vtk_radial(string &Name_Bathymetry_File, int Ma, int i_radial, int n);
    void paraview_vtk_zonal(string &Name_Bathymetry_File, int k_zonal, int n);
    void paraview_vtk_longal(string &Name_Bathymetry_File, int j_longal, int n); 
    void Atmosphere_v_w_Transfer(const string &Name_Bathymetry_File);
    void Atmosphere_PlotData(string &Name_Bathymetry_File, int iter_cnt);
    void run_data_atm();
    void searchMinMax_2D(string, string, 
        string, Array_2D &, double coeff=1.);
    void searchMinMax_3D(string, string, 
        string, Array &, double coeff=1., 
        std::function< double(double) > lambda = default_lambda,
        bool print_heading=false);

    double out_maxValue() const;
    double out_minValue() const;

    static cAtmosphereModel* m_model;

    PythonStream ps;
    std::streambuf *backup;

    std::set<float> m_time_list;
    std::set<float>::const_iterator m_current_time;

    int iter_cnt, iter_cnt_3d, iter_cnt_2d; // iteration count

    string bathymetry_name;

    double coeff_mmWS;    // coeff_mmWS = 1.2041/0.0094 [ kg/m³/kg/m³ ] = 128,0827 [/]
    double max_Precipitation;
    double maxValue;
    double minValue;

    bool is_node_weights_initialised;
    bool has_welcome_msg_printed;

    std::vector<std::vector<double> > m_node_weights;
    std::vector<float> m_layer_heights;

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
    Array_2D radiation_surface; // direct sun radiation, short wave
    Array_2D temperature_NASA; // surface temperature from NASA
    Array_2D temp_NASA; // surface temperature from NASA for print function
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
//    Array_2D aux_2D; // local surface water vapour fixed for iterations
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
    Array stream; // mass stream function
    Array u_stream; // u-velocity by mass stream function
    Array p_stat; // static pressure
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
    Array P_conv; // rain formation by cloud convection
    Array S_v; // water vapour mass rate due to category two ice scheme
    Array S_c; // cloud water mass rate due to category two ice scheme
    Array S_i; // cloud ice mass rate due to category two ice scheme
    Array S_r; // rain mass rate due to category two ice scheme
    Array S_s; // snow mass rate due to category two ice scheme
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
