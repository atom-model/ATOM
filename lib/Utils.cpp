#include <vector>
#include <Utils.h>

using namespace AtomUtils;

std::ofstream& AtomUtils::logger(){
    static std::ofstream logger("atom_log.txt", std::ofstream::out);
    return logger;
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

//change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
void AtomUtils::move_data(double* data, int len)
{
    for(int i=0; i<len/2; i++){
        std::iter_swap(data+i, data+len/2+i);
    }
    data[len-1] = data[0];
}

//the size of arrays must be the same with the size of new_arrays
//the values in the new_arrays will be assigned to coeff * old_values
void AtomUtils::move_data_to_new_arrays( int im, int jm, int km, double coeff, std::vector<Array*>& arrays, 
                               std::vector<Array*>& new_arrays)
{
    assert(arrays.size() == new_arrays.size());

    for ( int i = 0; i < im; i++ )
    {
        move_data_to_new_arrays( jm, km, coeff, arrays, new_arrays, i);
    }
}


void AtomUtils::move_data_to_new_arrays( int jm, int km, double coeff, std::vector<Array*>& arrays, 
                               std::vector<Array*>& new_arrays, int i)
{
    assert(arrays.size() == new_arrays.size());

    for ( int j = 0; j < jm; j++ )
    {
        for ( int k = 0; k < km; k++ )
        {
            for(std::size_t n=0; n<arrays.size(); n++)
            {
                new_arrays[n]->x[ i ][ j ][ k ] = coeff * arrays[n]->x[ i ][ j ][ k ];
            }
        }
    }
}


std::tuple<double, int, int, int>
AtomUtils::max_diff(int im, int jm, int km, const Array &a1, const Array &a2)
{
    std::tuple<double, int, int, int> ret;
    for ( int i = 0; i < im; i++ )
    {
        for ( int j = 0; j < jm; j++ )
        {
            for ( int k = 0; k < km; k++ )
            {
                double tmp = fabs ( a1.x[ i ][ j ][ k ] - a2.x[ i ][ j ][ k ] );
                if ( tmp > std::get<0>(ret) )
                {
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


double AtomUtils::simpson(int &n1, int &n2, double &dstep, Array_1D &value){
        double sum_even=0, sum_odd=0;
        if (n2 % 2 == 0){
            for (int i = n1 + 1; i < n2; i+=2){sum_odd += 4*value.z[i];}
            for (int i = n1 + 2; i < n2; i+=2){sum_even += 2*value.z[i];}
        }else cout << "       n2    must be an even number to use the Simpson integration method" << endl;
    return dstep/3 * (value.z[n1] + sum_odd + sum_even + value.z[n2]); // Simpson Rule integration
}


double AtomUtils::trapezoidal(int &n1, int &n2, double &dstep, Array_1D &value){
        double sum=0;
        for (int i = n1+1; i < n2; i++){sum += 2*value.z[i];}
        return dstep/2 * (value.z[n1] + sum + value.z[n2]);  // Trapezoidal Rule integration
    }

double AtomUtils::rectangular(int &n1, int &n2, double &dstep, Array_1D &value){
        double sum = 0;
        for (int i = n1; i <= n2; i++){sum += value.z[i];}
        return dstep * sum;  // Rectangular Rule integration
    }
