#include"cAtmosphereModel.h"
#include"Utils.h"
/*
 *
*/
double cAtmosphereModel::get_temperature(int time, int height, int lon, int lat){
    double nan = std::numeric_limits<double>::quiet_NaN();
    if(m_t_s.empty() || m_t_s.find(time) == m_t_s.end() || height < 0)
        return nan;

    int i = get_layer_index(height);
    int j = AtomUtils::lat_2_index(lat);
    int k = AtomUtils::lon_2_index(lon);

    if(!(0<=i&&i<im) || !(0<=j&&j<jm) || !(0<=k&&k<km))
        return nan;

    return (m_t_s[time].back().x[i][j][k]-1)*t_0;//273.15
}

