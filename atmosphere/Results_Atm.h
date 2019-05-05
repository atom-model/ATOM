/*
 * Atmosphere General Circulation Modell ( AGCM ) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to produce resulting data on mean sea level
*/

#include <iostream>
#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"

#ifndef Results_MSL_ATM
#define Results_MSL_ATM

using namespace std;

class Results_MSL_Atm
{
    private:
        int im, jm, km;
        int j_loc, k_loc, i_loc_level, j_loc_deg, k_loc_deg;
        int iter_prec;

        double E;
        double E_Rain, E_Rain_super, E_Ice, q_Rain, q_Rain_n, q_Rain_super, q_Ice, q_Ice_n;
        double e, a, e_SL, a_SL, p_SL;
        double t_dew_SL, t_Celsius_SL;
        double sun, albedo_equator, q_T;
        double e_h, a_h, p_h, q_h, t_dew, t_Celsius, t_Celsius_1, t_denom, Delta, E_a, gamma, g, gam;
        double i_level, h_level, h_h, sat_deficit, RF_e;
        double Evaporation_Penman_average, Evaporation_Dalton_average;
        double ep, hp, u_0, p_0, t_0, c_0, co2_0, sigma, albedo_extra, lv, ls, cp_l, r_air, dt, dr;
        double L_atm, c13, c43;
        double R_Air, height;
        double r_water_vapour, R_WaterVapour, precipitablewater_average, precipitation_average, precipitation_NASA_average, co2_vegetation_average;
        double coeff_mmWS, coeff_lv, coeff_ls, coeff_Diffusion_latent, coeff_Diffusion_sensibel, f_Haude, f_Penman, co2_vegetation, co2_ocean, co2_land;
        double r_humid, r_dry, lat_av_loc, sen_av_loc;
        double Value_1, Value_2, Value_3, Value_4, Value_5, Value_6, Value_7, Value_8, Value_9, Value_10, Value_12, Value_13;
        double Value_14, Value_15, Value_16, Value_17, Value_18, Value_19, Value_23, Value_24, Value_25, Value_26, Value_27;
        double dthe, dphi;
        double alf_ev;
        double S_nuc, S_ac, S_rim, S_shed, S_ev;
        double S_i_dep;
        double S_if_frz, S_cf_frz;
        double S_c_frz, S_c_au, S_i_au, S_d_au, S_agg, S_i_cri, S_r_cri, S_s_dep, S_i_melt, S_s_melt;
        double a_if, c_ac, c_rim, bet_ev, alf_melt, bet_melt, bet_dep, bet_s_dep, alf_if, alf_cf, a_s_melt, b_s_melt, E_cf, N_cf, N_cf_0, N_cf_0_surf, N_cf_0_6km;
        double tau_r, tau_s, t_1, t_2, t_m1, t_m2;
        double a_mc, a_mv;
        double coeff_P, a_i_m, a_s_m, N_r_0, N_s_0;
        double N_i_0, N_i, t_nuc, t_d, t_hn, m_i, m_i_0, m_i_max, m_s_0, c_i_dep, c_c_au, c_i_au, c_agg, c_i_cri, c_r_cri, c_s_dep, c_s_melt;
        double q_v_b, q_i_b;
        double q_c_b;
        double b_u, alf_1, alf_2, p_ps, bet_p;
        double t_u;
        double T, T_nue, T_tilda_h;
        double q_v_hyp, q_v_hyp_n;
        double T_t_in;
        double T_t_in0, T_t_in1, DEP, CND,d_q_v, d_q_c, d_q_i, d_t, p_t_in, t_Celsius_0, E_Rain_t_in, q_Rain_t_in, t_pole, t_cretaceous, temperature_surf_average, t_average, velocity_average, temperature_average, latent_heat_average, sensible_heat_average;
        double *e_d, *e_l, *e_p, *g_p, *c_u, *r_humid_u, *r_humid_u_parc, *r_dry_u, *r_dry_u_parc, *vel_av, *temp_av, *lat_av, *sen_av;

        string name_Value_1, name_Value_2, name_Value_3, name_Value_4, name_Value_5, name_Value_6, name_Value_7, name_Value_8, name_Value_9, name_Value_10, name_Value_11, name_Value_12, name_Value_13, name_Value_14, name_Value_15, name_Value_16, name_Value_17, name_Value_18, name_Value_19, name_Value_20, name_Value_21, name_Value_22, name_Value_23, name_Value_24, name_Value_25, name_Value_26, name_Value_27, name_unit_wm2, name_unit_mm, name_unit_mmd, name_unit_mma, name_unit_ppm, name_unit_t;
        string level, deg_north, deg_south, deg_west, deg_east, deg_lat, deg_lon, heading, heading_Dresden, heading_Sydney, heading_Pacific;

    public:
        Results_MSL_Atm ( int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double ); 


        ~Results_MSL_Atm (  );

        void run_MSL_data ( int, int, int, double &, Array_1D &, Array_1D &, Array_1D &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array_2D &, Array &, Array &, Array &, Array &, Array &, Array & );

        std::vector<std::vector<double> > m_node_weights;

        float GetMean_3D(int jm, int km, Array &val_3D);
        float GetMean_2D(int jm, int km, Array_2D &val_2D);

        void CalculateNodeWeights(int jm, int km);

        double C_Dalton ( double u_0, double v, double w );
};
#endif
