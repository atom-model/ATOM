#include <sstream>

#include <Utils.h>

using namespace AtomUtils;

std::ofstream& AtomUtils::get_logger(){
    static std::ofstream o("atom_log.txt", std::ofstream::out);
    return o;
}

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

//change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
void AtomUtils::move_data(double* data, int len){
    for(int i=0; i<len/2; i++){
        std::iter_swap(data+i, data+len/2+i);
    }
    data[len-1] = data[0];
}

void AtomUtils::move_data(std::vector<int>& data, int len){   
    for(int i=0; i<len/2; i++){
        std::iter_swap(data.begin()+i, data.begin()+len/2+i);
    }
    data[len-1] = data[0];
}

//the size of arrays must be the same with the size of new_arrays
//the values in the new_arrays will be assigned to coeff * old_values
void AtomUtils::move_data_to_new_arrays(int im, int jm, int km, double coeff, 
    std::vector<Array*>& arrays, std::vector<Array*>& new_arrays){
    assert(arrays.size() == new_arrays.size());
    for ( int i = 0; i < im; i++ ){
        move_data_to_new_arrays(jm, km, coeff, arrays, new_arrays, i);
    }
}


void AtomUtils::move_data_to_new_arrays(int jm, int km, double coeff, std::vector<Array*>& arrays, 
                               std::vector<Array*>& new_arrays, int i){
    assert(arrays.size() == new_arrays.size());
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++ ){
            for(std::size_t n=0; n<arrays.size(); n++){
                new_arrays[n]->x[ i ][ j ][ k ] = coeff * arrays[n]->x[ i ][ j ][ k ];
            }
        }
    }
}


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


double AtomUtils::C_Dalton(double u_0, double v, double w){
    // variation of the heat transfer coefficient in Dalton's evaporation law, parabola
    // air velocity measured 2m above sea level
    double C_max = .053;  // for v_max = 10 m/s, but C is function of v, should be included
    // Geiger ( 1961 ) by > Zmarsly, Kuttler, Pethe in mm/( h * hPa ), p. 133
    double v_max = 10.;  // Geiger ( 1961 ) by Zmarsly, Kuttler, Pethe in m/s, p. 133
    double fac = 1.;  // factor to adjust the ratio of NASA precipitation 
                       // to Dalton evaporation for the modern world, 
                       // for the modern world the global precipitation is 10% higher than evaporation
    double vel_magnitude = sqrt(v * v + w * w) * u_0;
    return fac * sqrt(C_max * C_max / v_max * vel_magnitude);  // result in mm/h
}


void AtomUtils::read_IC(const string& fn, double** a, int jm, int km){
    ifstream ifs(fn);
    if(!ifs.is_open()){
        cerr << "ERROR: unable to open " << fn << "\n";
        abort();
    }
    double lat, lon, d;
    for(int k=0; k < km && !ifs.eof(); k++){
        for(int j=0; j < jm; j++){
            ifs >> lat >> lon >> d;
            a[ j ][ k ] = d;
        }
    }
}


double AtomUtils::smooth_steps(int k, int im, int jm, Array &value_in, 
    Array &value_out){
//    double coeff = .1;
    double coeff = .5;
//    double coeff = .0;
    std::vector<std::vector<double> > inter(im, std::vector<double>(jm, 0));
    for(int i = 0; i < im ; i++){
        for(int j = 0; j < jm; j++){
            inter[i][j] = value_in.x[i][j][k];
/*
            if(i==35 && j==90)  
                cout << "   smooth_steps 1   " 
                    << "   i = " << i << "   j = " << j << "   k = " << k 
                    << "   inter = " << inter[i][j] 
                    << "   value_in = " << value_in.x[i][j][k] 
                    << "   value_out = " << value_out.x[i][j][k] << endl;
*/
        }
    }
    for(int i = 1; i < im-1 ; i++){
        for(int j = 1; j < jm-1; j++){
/*
            value_out.x[i][j][k] = ( coeff * value_in.x[i][j][k] 
                + ( 1. - coeff ) * ( value_in.x[i+1][j][k] 
                + value_in.x[i-1][j][k] + value_in.x[i][j+1][k] 
                + value_in.x[i][j-1][k] ) / 4. );
*/
            value_out.x[i][j][k] = ( coeff * inter[i][j] 
                + ( 1. - coeff ) * ( inter[i+1][j] 
                + inter[i-1][j] + inter[i][j+1] 
                + inter[i][j-1] ) / 4. );
/*
        if(i==35 && j==90)  
            cout << "   smooth_steps 2   " 
                << "   i = " << i << "   j = " << j << "   k = " << k 
                << "   inter = " << inter[i][j] 
                << "   value_in = " << value_in.x[i][j][k] 
                << "   value_out = " << value_out.x[i][j][k] << endl;
*/
        } 
    }
    return value_out.x[0][0][0];
}


double AtomUtils::smooth_tropopause(int jm, double* value_in, double* value_out){
//    double coeff = .1;
    double coeff = .5;
//    double coeff = .0;
    std::vector<double>inter(jm, 0);
    for(int j = 0; j < jm; j++){
        inter[j] = value_in[j];
        cout << "   smooth_tropopause 1   " 
            << "   j = " << j << "   inter = " << inter[j] 
            << "   value_in = " << value_in[j] 
            << "   value_out = " << value_out[j] << endl;
    }
//    for(int j = 1; j < jm-1; j++){
    for(int j = 2; j < jm-2; j++){
//        value_out[j] = ( coeff * value_in[j] + ( 1. - coeff ) 
//            * ( value_in[j+1] + value_in[j-1] ) / 2. );
        value_out[j] = ( coeff * inter[j] + ( 1. - coeff ) 
//            * ( inter[j+1] + inter[j-1] ) / 2. );
            * ( inter[j+1] + inter[j-1] + inter[j+2] + inter[j-2] ) / 4. );

    cout << "   smooth_tropopause 2   " 
        << "   j = " << j << "   inter = " << inter[j] 
        << "   value_in = " << value_in[j] 
        << "   value_out = " << value_out[j] << endl;

    } 
    return value_out[0];
}


double AtomUtils::simpson(int n1, int n2, double dstep, Array_1D &value){
        double sum_even=0, sum_odd=0;
        if(n2 % 2 == 0){
            for(int i = n1 + 1; i < n2; i+=2){sum_odd += 4*value.z[i];}
            for(int i = n1 + 2; i < n2; i+=2){sum_even += 2*value.z[i];}
        }else cout << "       n2    must be an even number to use the Simpson integration method" << endl;
    return dstep/3 * (value.z[n1] + sum_odd + sum_even + value.z[n2]); // Simpson Rule integration
}


double AtomUtils::trapezoidal(int n1, int n2, double dstep, Array_1D &value){
        double sum=0;
        for(int i = n1+1; i < n2; i++){sum += 2*value.z[i];}
        return dstep/2 * (value.z[n1] + sum + value.z[n2]);  // Trapezoidal Rule integration
}


double AtomUtils::rectangular(int n1, int n2, double dstep, Array_1D &value){
        double sum = 0;
        for(int i = n1; i <= n2; i++){sum += value.z[i];}
        return dstep * sum;  // Rectangular Rule integration
}


void AtomUtils::load_map_from_file(const std::string& fn, std::map<float, float>& m){
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
}
