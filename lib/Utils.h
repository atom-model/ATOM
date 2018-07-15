#include <iostream>
#include <fstream>

namespace AtomUtils{
    using namespace std;
    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    }; 

    HemisphereCoords convert_coords(double lon, double lat);

    std::ofstream& logger();
}
