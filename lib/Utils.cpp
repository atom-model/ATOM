#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <Utils.h>
#include "cAtmosphereModel.h"
#include "cHydrosphereModel.h"

using namespace AtomUtils;
using namespace std;

std::vector<std::vector<double> > m_node_weights;

std::ofstream& AtomUtils::get_logger(){
    static std::ofstream o("atom_log.txt", std::ofstream::out);
    return o;
}
/*
*
*/
HemisphereCoords AtomUtils::convert_coords(double lon, double lat){
    HemisphereCoords ret;
    if ( lat > 90 ){
        ret.lat = lat - 90;
        ret.north_or_south = "°S";
    }else{
        ret.lat = 90 - lat;
        ret.north_or_south = "°N";
    }
    if ( lon > 180 ){
        ret.lon = 360 - lon;
        ret.east_or_west = "°W";
    }else{
        ret.lon = lon;
        ret.east_or_west = "°E";
    }
    return ret;
}
/*
* [-180, 180] to [0, 360]
*/
int AtomUtils::lon_2_index(double lon){
    if(lon<0)
        return int(round(360+lon));
    else
        return int(round(lon));
}
/*
* [-90, 90] to [0, 180]
*/
int AtomUtils::lat_2_index(double lat){
    return int(round(90-lat));
}
/*
*
*/
//change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
void AtomUtils::move_data(double* data, int len){
    for(int i=0; i<len/2; i++){
        std::iter_swap(data+i, data+len/2+i);
    }
    data[len-1] = data[0];
    return;
}
/*
*
*/
void AtomUtils::move_data(std::vector<int>& data, int len){   
    for(int i=0; i<len/2; i++){
        std::iter_swap(data.begin()+i, data.begin()+len/2+i);
    }
    data[len-1] = data[0];
    return;
}
/*
*
*/
//the size of arrays must be the same with the size of new_arrays
//the values in the new_arrays will be assigned to coeff * old_values
void AtomUtils::move_data_to_new_arrays(int im, int jm, int km, double coeff, 
    std::vector<Array*>& arrays, std::vector<Array*>& new_arrays){
    assert(arrays.size() == new_arrays.size());
    for ( int i = 0; i < im; i++ ){
        move_data_to_new_arrays(jm, km, coeff, arrays, new_arrays, i);
    }
    return;
}
/*
*
*/
void AtomUtils::move_data_to_new_arrays(int jm, int km, double coeff, 
    std::vector<Array*>& arrays, std::vector<Array*>& new_arrays, int i){
    assert(arrays.size() == new_arrays.size());
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++ ){
            for(std::size_t n=0; n<arrays.size(); n++){
                new_arrays[n]->x[ i ][ j ][ k ] = coeff * arrays[n]->x[ i ][ j ][ k ];
            }
        }
    }
    return;
}
/*
*
*/
std::tuple<double, int, int, int>
AtomUtils::max_diff(int im, int jm, int km, const Array &a1, const Array &a2){
    std::tuple<double, int, int, int> ret;
    for ( int i = 0; i < im; i++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                double tmp = fabs ( a1.x[ i ][ j ][ k ] - a2.x[ i ][ j ][ k ] );
                if(tmp > std::get<0>(ret)){
                    std::get<0>(ret) = tmp;
                    std::get<1>(ret) = i;
                    std::get<2>(ret) = j;
                    std::get<3>(ret) = k;
                }
            }
        }
    }
    return ret;
}
/*
*
*/
double AtomUtils::C_Dalton(int i, int j, int k, double coeff_Dalton, 
    double u_0, Array &v, Array &w){
    // variation of the evaporation coefficient in Dalton's evaporation law, parabola
    // air velocity measured 2m above sea level
    // for v_max = 10 m/s, but C_Dalton is function of v, should be included
    // Geiger ( 1961 ) by > Zmarsly, Kuttler, Pethe in mm/( h * hPa ), p. 133
    double v_max = 10.;  // in m/s ..... Geiger ( 1961 ) by Zmarsly, Kuttler, Pethe in m/s, p. 133
                       // for the modern world the global precipitation is 10% higher than evaporation
    double c_Dalton_max = 0.053; // in mm/(h * hPa)
    double vel_magnitude = sqrt(v.x[i][j][k] * v.x[i][j][k] 
        + w.x[i][j][k] * w.x[i][j][k]) * u_0; // in m/s
    return coeff_Dalton * sqrt(c_Dalton_max * c_Dalton_max * vel_magnitude/v_max);  // result in mm/(h * hPa)
}
/*
*
*/
/*
double AtomUtils::Humility_critical(double x, double Hu_cr_max, double Hu_cr_mid){
    // critical humidity: assumption of a parabolic dependency in vertical direction 
    // by Jean Coiffier, Fundamentals of Numerical Weather Prediction,p 210
    // decreases the amount of cloud water, at a static pressure value of 550 hPa exists only 80%
        double Hu_cr_max;
        double Hu_cr_mid;
        double Hu_cr = (Hu_cr_max - Hu_cr_mid) * (x * x - 2.0 * x) 
            + Hu_cr_max;
    return Hu_cr;
}
*/
/*
*
*/
void AtomUtils::read_IC(const string& fn, double** a, int jm, int km){
    cout << endl << "      read_IC" << endl;
    cout << "      " << fn << endl;
    ifstream ifs(fn);
    if(!ifs.is_open()){
        cerr << "ERROR: unable to open " << fn << "\n";
        abort();
    }
    double lat, lon, d;
    for(int k=0; k < km && !ifs.eof(); k++){
        for(int j=0; j < jm; j++){
            ifs >> lat >> lon >> d;
            a[j][k] = d;
        }
    }
    cout << "      read_IC ended" << endl;
    return;
}
/*
*
*/
void AtomUtils::smooth_steps(int k, int im, int jm, Array &value){
    double coeff = 0.5;
    std::vector<std::vector<double> > inter(im, std::vector<double>(jm, 0));
    for(int i = 0; i < im ; i++){
        for(int j = 0; j < jm; j++){
            inter[i][j] = value.x[i][j][k];
        }
    }
    for(int i = 1; i < im-1 ; i++){
        for(int j = 1; j < jm-1; j++){
            value.x[i][j][k] = (coeff * inter[i][j] 
                + (1. - coeff) * (inter[i+1][j] 
                + inter[i-1][j] + inter[i][j+1] 
                + inter[i][j-1])/4.0);
        } 
    }
    return;
}
/*
*
*/
void AtomUtils::smooth_tropopause(int jm, std::vector<double> &value){
    double coeff = 0.5;
    std::vector<double>inter(jm, 0);
    for(int j = 0; j < jm; j++){
        inter[j] = value[j];
    }
//    for(int j = 1; j < jm-1; j++){
//        value[j] = (coeff * inter[j] + (1.0 - coeff) 
//            * (inter[j+1] + inter[j-1])/2.0);
    for(int j = 2; j < jm-2; j++){
        value[j] = (coeff * inter[j] + (1. - coeff) 
            * (inter[j+1] + inter[j-1] + inter[j+2] + inter[j-2])/4.0);
    } 
    return;
}
/*
*
*/
double AtomUtils::simpson(int n1, int n2, double dstep, Array_1D &value){
        double sum_even=0, sum_odd=0;
        if(n2 % 2 == 0){
            for(int i = n1+1; i < n2; i+=2){sum_odd += 4 * value.z[i];}
            for(int i = n1+2; i < n2; i+=2){sum_even += 2 * value.z[i];}
        }else cout << "       n2    must be an even number to use the \
            Simpson integration method" << endl;
    return dstep/3 * (value.z[n1] + sum_odd + sum_even + value.z[n2]); // Simpson Rule integration
}
/*
*
*/
double AtomUtils::trapezoidal(int n1, int n2, double dstep, Array_1D &value){
        double sum = 0;
        for(int i = n1+1; i < n2; i++){sum += 2 * value.z[i];}
        return dstep/2 * (value.z[n1] + sum + value.z[n2]);  // Trapezoidal Rule integration
}
/*
*
*/
double AtomUtils::rectangular(int n1, int n2, double dstep, Array_1D &value){
        double sum = 0;
        for(int i = n1; i <= n2; i++){sum += value.z[i];}
        return dstep * sum;  // Rectangular Rule integration
}
/*
*
*/
void AtomUtils::load_map_from_file(const std::string& fn, 
    std::map<float, float>& m){
    std::string line;
    std::ifstream f(fn);
    if(!f.is_open()){
        std::cout << "error occurred while opening file: "<< fn << std::endl;
    }
    float key=0, val=0;
    while(getline(f, line)){
        if(line.find("#")==0) continue;
        std::stringstream(line) >> key >> val;
        m.insert(std::pair<float,float>(key, val));
        //std::cout << key <<"  " << val << std::endl;
    }
    return;
}
/*
*
*/
double AtomUtils::GetMean_3D(int jm, int km, Array &val_3D){
    if(m_node_weights.size() != (unsigned)jm){
        CalculateNodeWeights(jm, km);
    }
    double ret=0.0, weight=0.0;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (val_3D.x[0][j][k]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret += val_3D.x[0][j][k] * m_node_weights[j][k];
            weight += m_node_weights[j][k];
        }
    }
    return ret/weight;
}
/*
*
*/
double AtomUtils::GetMean_2D(int jm, int km, Array_2D &val_2D){
    if(m_node_weights.size() != (unsigned)jm){
        CalculateNodeWeights(jm, km);
    }
    double ret=0.0, weight=0.0;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            //std::cout << (val_2D.y[j][k]-1)*t_0 << "  " << m_node_weights[j][k] << std::endl;
            ret += val_2D.y[j][k] * m_node_weights[j][k];
            weight += m_node_weights[j][k];
        }
    }
    return ret/weight;
}
/*
*
*/
void AtomUtils::CalculateNodeWeights(int jm, int km){
    //use cosine of latitude as weights for now
    //longitudes: 0-360(km) latitudes: 90-(-90)(jm)
    double weight = 0.0;
    m_node_weights.clear();
    for(int i=0; i<jm; i++){
        if(i<=90){
            weight = cos((90-i) * M_PI/180.0);
        }else{
            weight = cos((i-90) * M_PI/180.0);
        }
        m_node_weights.push_back(std::vector<double>());
        m_node_weights[i].resize(km, weight);
    }
    return;
}
/*
*
*/
void AtomUtils::land_oceanFraction(int i_local, int jm, int km, Array &val_3D){
    int h_point_max = (jm-1) * (km-1);
    int h_land = 0;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_land( val_3D, i_local, j, k))  h_land++;
        }
    }
    int h_ocean = h_point_max - h_land;
    float land_percent = 100.0 * (float)h_land/(float)h_point_max;
    float ocean_percent = 100.0 * (float)h_ocean/(float)h_point_max;

    cout.precision(3);
    cout << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      total number of points at constant hight " << " = " 
        << resetiosflags(ios::left)
        << setw(7) << fixed << setfill(' ') << h_point_max << endl 
        << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      number of points on the ocean surface " << " = "
        << resetiosflags(ios::left) 
        << setw(7) << fixed << setfill(' ') << h_ocean << endl
        << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      number of points on the land surface " << " = " 
        << resetiosflags(ios::left) 
        << setw(7) << fixed << setfill(' ') << h_land << endl 
        << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      ocean in % " << " = " << resetiosflags(ios::left) 
        << setw(7) << fixed << setfill(' ')
        << ocean_percent << endl
        << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      land in % " << " = " << resetiosflags(ios::left) 
        << setw(7) << fixed << setfill(' ')
        << land_percent << endl << endl;
    cout << endl;
    return;
}
/*
*
*/
int AtomUtils::RunStart(string comment){
    string at = "AGCM";
    string hy = "OGCM";
    string comment_1 = "";
    string comment_2 = "";
    if(comment.compare(at) == 0){
        comment_1 = " ... AGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
             Atmosphere General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    if(comment.compare(hy) == 0){
        comment_1 = " ... OGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
            Ocean General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    std::time_t Run_start;
    struct tm * timeinfo_begin;
    std::time(&Run_start);
    timeinfo_begin = std::localtime(&Run_start);
    std::cout << std::endl << std::endl;
    std::cout << comment_2 
        << std::endl << std::endl;
    std::cout << comment_1 
        << std::asctime(timeinfo_begin);
    return Run_start;
}
/*
*
*/
int AtomUtils::RunEnd(string comment, int Ma, int Run_start){
    string at = "AGCM";
    string hy = "OGCM";
    string comment_1 = "";
    string comment_2 = "";
        if(comment.compare(at) == 0){
        comment_1 = " ... AGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
             Atmosphere General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    if(comment.compare(hy) == 0){
        comment_1 = " ... OGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
            Ocean General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    std::time_t Run_end;
    std::time(&Run_end);
    struct tm * timeinfo_end;
    timeinfo_end = std::localtime(&Run_end);
    std::cout << std::endl << std::endl;
    std::cout << comment_2 
        << std::endl << std::endl;
    std::cout << comment_1 
        << std::asctime(timeinfo_end);
    int Run_total = Run_end - Run_start;
    int Run_total_minutes = Run_total/60;
    int Run_total_hours = Run_total_minutes/60;
    std::cout << std::endl << " ... computed time slice:"
        << "  Ma = " << Ma << std::endl << std::endl
        << " ... computer time needed:" << std::endl << std::endl
        << setw(20) << setfill(' ') << Run_total 
        << " seconds" << std::endl
        << " ... compares to:" << std::endl << std::endl
        << setw(20) << setfill(' ') << Run_total_hours 
        << " hours" << std::endl
        << setw(20) << setfill(' ') << Run_total_minutes 
        << " minutes" << std::endl
        << setw(20) << setfill(' ') << Run_total%60 
        << " seconds" << std::endl << std::endl
        << std::endl;
    return 0;
}
/*
*
*/
void AtomUtils::use_presets(bool preset_1, bool preset_2, bool preset_3){
    if(preset_1)
        cout << endl << "      use_earthbyte_reconstruction is in use" << endl;
    else
        cout << endl << "      use_earthbyte_reconstruction is not in use" << endl;
    if(preset_2)
        cout << "      use_NASA_temperature is in use" << endl;
    else
        cout << "      use_NASA_temperature is not in use" << endl;
    if(preset_3)
        cout << "      use_NASA_velocity is in use" << endl;
    else
        cout << "      use_NASA_velocity is not in use" << endl;
    return;
}
