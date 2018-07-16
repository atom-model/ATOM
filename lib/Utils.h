#include <iostream>
#include <fstream>
#include <limits>

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

    inline bool is_ocean_surface(const Array& h, int i, int j, int k){
        return i==0 && !is_land(h, i, j, k);
    }

    inline bool is_land_surface(const Array& h, int i, int j, int k){
        return is_land(h, i, j, k) && !is_land(h, i+1, j, k);
    }
}
