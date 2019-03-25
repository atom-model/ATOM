#include <iostream>
#include <fstream>
#include <limits>

#include "Array.h"
#include "Array_1D.h"

namespace AtomUtils{
    using namespace std;
    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    }; 

    HemisphereCoords convert_coords(double lon, double lat);

    std::ofstream& logger();

    inline bool is_land(const Array& h, int i, int j, int k){
        return fabs(h.x[i][j][k] - 1) < std::numeric_limits<double>::epsilon();
    }

    inline bool is_air(const Array& h, int i, int j, int k){
        return !is_land(h,i,j,k);
    }

    inline bool is_water(const Array& h, int i, int j, int k){
        return !is_land(h,i,j,k);
    }

    inline bool is_ocean_surface(const Array& h, int i, int j, int k){
        return i==0 && !is_land(h, i, j, k);
    }

    inline bool is_land_surface(const Array& h, int i, int j, int k){
        return is_land(h, i, j, k) && !is_land(h, i+1, j, k);
    }

    inline double parabola(double x){
        return x*x - 2*x;
    }

    //change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
    void move_data(double* data, int len);

    void move_data_to_new_arrays(int im, int jm, int km, double coeff, std::vector<Array*>& arrays, 
                        std::vector<Array*>& new_arrays);

    void move_data_to_new_arrays( int jm, int km, double coeff, std::vector<Array*>& arrays, 
                        std::vector<Array*>& new_arrays, int i = 0);

    std::tuple<double, int, int, int>
    max_diff(int i, int j, int k, const Array &a1, const Array &a2);

    // integration tools
    double simpson(int &n1, int &n2, double &dstep, Array_1D &value);
    double trapezoidal(int &n1, int &n2, double &dstep, Array_1D &value);
    double rectangular(int &n1, int &n2, double &dstep, Array_1D &value);

    inline double exp_func(double &T_K, const double &co_1, const double &co_2){
        return exp(co_1 * (T_K - 273.15) / (T_K - co_2));                        // temperature in °K
    }
}
