#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;
/*
*
*/
void cAtmosphereModel::init_velocities(){
    // boundary condition for the velocity components in the circulation cells

    // latest version by Grotjahn (Global Atmospheric Circulations, 1993)
    // default for the velocity components u, v, and w as initial conditions

    // velocities given in m/s, 1 m/s compares to 3.6 km/h, non-dimensionalized by u_0 at the end of this class element
    // do not change the velocity initial conditions !!

    // initial velocity components in the northern and southern
    // Pole, Ferrel and Hadley cells

    // equator (at j=90 compares to 0° latitude)
    // u-component up to tropopause and back on half distance (i = 20)

    init_u(u,0);
    init_u(u,30);
    init_u(u,60);
    init_u(u,90);
    init_u(u,120);
    init_u(u,150);
    init_u(u,180);

    //initialise v, tropopause, surface
    //equator
    init_v_or_w(v, 90, 0, 0); // lat:0
    
    //polar cell
    //northern polar cell
    init_v_or_w(v, 0, 0, 0);  //lat: 90
//    init_v_or_w(v, 15, -1, 0.5); //lat: 75
    init_v_or_w(v, 15, -.1, 0.05); //lat: 75

    //southern polar cell
    init_v_or_w(v, 180, 0, 0); //lat: -90
//    init_v_or_w(v, 165, -1, 0.5); //lat: -75
    init_v_or_w(v, 165, -.1, 0.05); //lat: -75

    //Ferrel cell
    //northern Ferrel cell
//    init_v_or_w(v, 30, 1, 0.5); //lat: 60
    init_v_or_w(v, 30, .2, 0.1); //lat: 60
    init_v_or_w(v, 45, 1, -0.1); //lat: 45   

    //southern Ferrel cell
//    init_v_or_w(v, 150, 1, 0.5); //lat: -60
    init_v_or_w(v, 150, .2, 0.1); //lat: -60
    init_v_or_w(v, 135, 1, -0.1); //lat: -45   

    // Hadley cell
    //northern Hadley cell
    init_v_or_w(v, 60, -1, 0.25); //lat: 30
    init_v_or_w(v, 75, -1, 1); //lat: 15   

    //southern Hadley cell
    init_v_or_w(v, 120, -1, 0.25); //lat: -30
    init_v_or_w(v, 105, -1, 1); //lat: -15 

    //initialise w, tropopause, surface
    //equator
    init_v_or_w(w, 90, -7.5, -1); //lat: 0
//    init_v_or_w(w, 90, -7.5, -.5); //lat: 0
    
    //polar cell
    //northern polar cell
    init_v_or_w(w, 0, 0, -0.01); //lat: 90

    //southern polar cell
    init_v_or_w(w, 180, 0, -0.01); //lat: -90

    //Ferrel cell
    //northern Ferrel cell
    init_v_or_w(w, 30, 10, -0.2); //lat: 60

    //southern Ferrel cell
    init_v_or_w(w, 150, 10, -0.2); //lat: -60
 
    // Hadley cell
    //northern Hadley cell
    init_v_or_w(w, 60, 30, 1); //lat: 30

    //southern Hadley cell
    init_v_or_w(w, 120, 30, 1); //lat: -30

    // forming diagonals 
    //northen hemisphere
    form_diagonals(u, 0, 30);
    form_diagonals(w, 0, 30);
    form_diagonals(v, 0, 15);
    form_diagonals(v, 15, 30);

    form_diagonals(u, 30, 60);
    form_diagonals(w, 30, 60);
    form_diagonals(v, 30, 45);
    form_diagonals(v, 45, 60);

    form_diagonals(u, 60, 90);
    form_diagonals(w, 60, 90);
    form_diagonals(v, 60, 75);
    form_diagonals(v, 75, 90);

    //southen hemisphere
    form_diagonals(u, 90, 120);
    form_diagonals(w, 90, 120);
    form_diagonals(v, 90, 105);
    form_diagonals(v, 105, 120);

    form_diagonals(u, 120, 150);
    form_diagonals(w, 120, 150);
    form_diagonals(v, 120, 135);
    form_diagonals(v, 135, 150);

    form_diagonals(u, 150, 180);
    form_diagonals(w, 150, 180);
    form_diagonals(v, 150, 165);
    form_diagonals(v, 165, 180);

    //change the direction for southern hemisphere
    if(!use_NASA_velocity){
        for (int i = 0; i < im; i++){
            for (int j = 91; j < jm; j++){
                for (int k = 0; k < km; k++){
                    v.x[i][j][k] = - v.x[i][j][k];
                }
            }
        }
    }

    //smoothing transitions from cell to cell
    /*smooth_transition(u,v,w,60); 
    smooth_transition(u,v,w,30);    
    smooth_transition(u,v,w,75);
    smooth_transition(u,v,w,45);
    smooth_transition(u,v,w,90); 
    smooth_transition(u,v,w,120);
    smooth_transition(u,v,w,150);
    smooth_transition(u,v,w,105);
    smooth_transition(u,v,w,135);
    */
    // non dimensionalization by u_0
    for (int i = 0; i < im; i++){
        for (int k = 0; k < km; k++){
            for (int j = 0; j < jm; j++){
                if (is_land (h, i, j, k))     
                    u.x[i][j][k] = v.x[i][j][k] = w.x[i][j][k] = 0.;
                else{
                    u.x[i][j][k] = u.x[i][j][k] / u_0;
                    v.x[i][j][k] = v.x[i][j][k] / u_0;
                    w.x[i][j][k] = w.x[i][j][k] / u_0;
                }
            }
        }
    }
}
/*
*
*/
void cAtmosphereModel::smooth_transition(Array &u, Array &v, Array &w, int lat){
    int start = lat-3, end = lat+3;
    for (int i = 0; i < im; i++){
        for (int k = 0; k < km; k++){
            for (int j = start; j <= end; j++){
                u.x[i][j][k] = (u.x[i][end][k] - u.x[i][start][k]) *
                    ((double) (j - start)  / (end - start)) + u.x[i][start][k];
                v.x[i][j][k] = (v.x[i][end][k] - v.x[i][start][k]) *
                    ((double) (j - start)  / (end - start)) + v.x[i][start][k];
                w.x[i][j][k] = (w.x[i][end][k] - w.x[i][start][k]) *
                    ((double) (j - start)  / (end - start)) + w.x[i][start][k];
            }
        }
    }
}

/*
*
*/
void cAtmosphereModel::form_diagonals(Array &a, int start, int end){
    for (int k = 0; k < km; k++){
        for (int j = start; j < end; j++){
            for (int i = 0; i < im; i++){
                a.x[i][j][k] = (a.x[i][end][k] - a.x[i][start][k]) *
                    (j - start) / (double)(end - start) + a.x[i][start][k];
            }
        }
    }
}

/*
*
*/
void  cAtmosphereModel::init_u(Array &u, int j){
    float ua_00 = 1,
//          ua_30 = 1,
          ua_30 = .8,
//          ua_60 = 0.5,
          ua_60 = 0.2,
//          ua_90 = 0.5;
          //ua_90 = 0.2;
          ua_90 = 0.05;
    int tropopause_layer = get_tropopause_layer(j);
    float tropopause_height = get_layer_height(tropopause_layer);
    for(int k = 0; k < km; k++){
        for(int i = 0; i < tropopause_layer; i++){
            float layer_height = get_layer_height(i), half_tropopause_height = tropopause_height/2.;
            float ratio;
            if(layer_height < half_tropopause_height){    
                ratio = layer_height / half_tropopause_height; 
                switch(j){
                    case 90:  u.x[i][90][k] = ua_00 * ratio; break;
                    case 60:  u.x[i][60][k] = - ua_30 * ratio; break;
                    case 120: u.x[i][120][k] = - ua_30 * ratio; break;
                    case 30:  u.x[i][30][k] = ua_60 * ratio; break;
                    case 150: u.x[i][150][k] = ua_60 * ratio; break;                
                    case 0:   u.x[i][0][k] = - ua_90 *  ratio ; break;
                    case 180: u.x[i][180][k] = - ua_90 * ratio; break;
                }
            }else{
                ratio = (tropopause_height-layer_height) / half_tropopause_height;
                switch(j){
                    case 90:  u.x[i][90][k] = ua_00 * ratio; break;
                    case 60:  u.x[i][60][k] = -ua_30 * ratio; break;
                    case 120: u.x[i][120][k] = -ua_30 * ratio; break;
                    case 150: u.x[i][150][k] = ua_60 * ratio; break;
                    case 30:  u.x[i][30][k] = ua_60 * ratio; break;
                    case 0:   u.x[i][0][k] = - ua_90 * ratio; break;
                    case 180: u.x[i][180][k] = - ua_90 * ratio; break;
                }
            }
        }
    }
}

/*
*
*/
void  cAtmosphereModel::init_v_or_w(Array &v_or_w, int j, double coeff_trop, double coeff_sl){
    int tropopause_layer = get_tropopause_layer(j);
    double tropopause_height = get_layer_height(tropopause_layer);
    for(int k = 0; k < km; k++){
        if(use_NASA_velocity && is_ocean_surface(h, 0, j, k))
        {
            coeff_sl = v_or_w.x[0][j][k];
        }
        if(!use_NASA_velocity && !is_ocean_surface(h, 0, j, k))
        {
            //coeff_sl = 0;
        }
        for(int i = 0; i < tropopause_layer; i++){
            v_or_w.x[i][j][k] = (coeff_trop - coeff_sl) *
                get_layer_height(i)/tropopause_height + coeff_sl;
        }
    }
    init_v_or_w_above_tropopause(v_or_w, j, coeff_trop);
}

/*
*
*/
void  cAtmosphereModel::init_v_or_w_above_tropopause(Array &v_or_w, int j, double coeff){
    int tropopause_layer = get_tropopause_layer(j);
    if(tropopause_layer >= im-1) return;
    double tropopause_height = get_layer_height(tropopause_layer);
    for(int k = 0; k < km; k++){
        for(int i = tropopause_layer; i < im; i++){
            v_or_w.x[i][j][k] = coeff * (get_layer_height(im-1) - get_layer_height(i)) / 
                (get_layer_height(im-1) - tropopause_height);
        }
    }
}
