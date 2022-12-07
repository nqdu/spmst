#include <cmath>
/**
 * @brief compute great circle distance between two points
 * 
 * @param lon0,lat0 coordinates for the first point, in deg 
 * @param lon1,lat1 coordinates for the second point, in deg 
 * @param radius earth radius 
 * @return float great circle distance 
 */
float 
gps2dist(float lon0, float lon1, float lat0,float lat1,float radius)
{
    const float deg2rad = M_PI / 180.0;
    using std::sin; using std::atan2;

    float dlat,dlon;
    dlat = (lat0 - lat1) * deg2rad; dlon = (lon1 - lon0) * deg2rad;
    float a = sin(dlat * 0.5) * sin(dlat * 0.5) + sin(dlon * 0.5) *
                sin(dlon * 0.5) * cos(lat0 * deg2rad) * cos(lat1 * deg2rad);
    float deg = 2.0 * atan2(std::sqrt(a),std::sqrt(1-a));

    return radius * deg;
}
