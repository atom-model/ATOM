#ifndef _UTILS_
#define _UTILS_

#include <iostream>
#include <fstream>
#include <limits>

#include "Array.h"

#define logger() \
if (false) ; \
else get_logger()

namespace AtomUtils{
    using namespace std;
    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    }; 

    HemisphereCoords convert_coords(double lon, double lat);

    int lon_2_index(double lon);
    int lat_2_index(double lat);

    std::ofstream& get_logger();

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
    void move_data(std::vector<int>& data, int len);

    void move_data_to_new_arrays(int im, int jm, int km, double coeff, std::vector<Array*>& arrays, 
                        std::vector<Array*>& new_arrays);

    void move_data_to_new_arrays( int jm, int km, double coeff, std::vector<Array*>& arrays, 
                        std::vector<Array*>& new_arrays, int i = 0);

    std::tuple<double, int, int, int>
    max_diff(int i, int j, int k, const Array &a1, const Array &a2);

    inline double simpson(int & n, double &dstep, double *value){
        double sum_even=0, sum_odd=0;
        if (n % 2 == 0){
            for (int i = 1; i < n; i+=2){sum_odd += 4*value[i];}
            for (int i = 2; i < n; i+=2){sum_even += 2*value[i];}
        }else cout << "       n    must be an even number to use the Simpson integration method" << endl;
        return dstep/3 * (value[0] + sum_odd + sum_even + value[n]);             // Simpson Rule integration
    }

    inline double trapezoidal(int n, double dstep, double *value){
        double sum=0;
        for (int i = 1; i < n; i++){sum += 2*value[i];}
		return dstep/2 * (value[0] + sum + value[n]);                 // Trapezoidal Rule integration
    }

    inline double rectangular(int n, double dstep, double *value){
        double sum = 0;
        for (int i = 0; i <= n; i++){sum += value[i];}
        return dstep * sum;                // Rectangular Rule integration
    }

    inline double exp_func(double T_K, const double co_1, const double co_2){
        return exp(co_1 * (T_K - 273.15) / (T_K - co_2));                        // temperature in °K
    }

    double C_Dalton ( double u_0, double v, double w );

}

#endif
