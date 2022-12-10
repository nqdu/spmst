#include "spm2d.hpp"

static float 
__cal_dist_spherical(float lon1,float lat1,float h1,
                        float lon2,float lat2,float h2)
{
    const float deg2rad = M_PI / 180.;
    float r1 = h1 + earth, r2 = h2 + earth;
    float x1 = r1 * std::cos(lat1 * deg2rad) * std::sin(lon1 * deg2rad);
    float x2 = r2 * std::cos(lat2 * deg2rad) * std::sin(lon2 * deg2rad);
    float y1 = r1 * std::cos(lat1 * deg2rad) * std::cos(lon1 * deg2rad);
    float y2 = r2 * std::cos(lat2 * deg2rad) * std::cos(lon2 * deg2rad);
    float z1 = r1 * std::sin(lat1 * deg2rad), z2 = r2 * std::sin(lat2 * deg2rad);

    float dist = (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2);
    dist = std::sqrt(dist);

    return dist;
}

/**
 * @brief compute the length of one ray segment
 * 
 * @param x1,y1,z1 first point 
 * @param x2,y2,z3 the second point 
 * @return float 
 */
float SPM2D::
compute_length(float x1,float y1,float z1,float x2,float y2,float z2) const
{
    float dist;
    if(is_spherical){
        dist =  __cal_dist_spherical(x1,y1,z1,x2,y2,z2);
    }
    else{
        dist = (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2);
        dist = std::sqrt(dist);
    }
   
   return dist;
}