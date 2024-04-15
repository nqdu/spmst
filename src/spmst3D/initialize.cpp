#include <iostream>
#include "numerical.hpp"
#include "spmst3D/spmst3D.hpp"
#include "shared/IO.hpp"
#include "shared/bilinear.hpp"
#include "shared/gps2dist.hpp"

/**
 * @brief compute two point distance according to different coordinate system
 * 
 * @param x0,y0 (x,y)/(lon,lat) for the first point
 * @param x1,y1  (x,y)/(lon,lat) for the second point
 * @return float 
 */
float SPMST3D::
compute_distance(float x0,float y0,float x1,float y1) const
{
    float dist{};
    if(is_spherical){
        dist = gps2dist(x0,x1,y0,y1,earth);
    }
    else{
        dist = std::hypot(x0-x1,y0-y1);
    }

    return dist;
}

/**
 * @brief read parameter file 
 * 
 * @param filename 
 */
void SPMST3D::
read_params(const char *filename)
{
    param.read_file(filename);
}

/**
 * @brief read topography 
 * 
 * @param filename filename of topography
 */
void SPMST3D:: 
read_topography(const char *filename)
{
    printf("\nreading topography from %s ...\n",filename);
    mesh.read_topography(filename);
    printf("min and max topography in km: %f %f\n",mesh.topo_z.minCoeff(),mesh.topo_z.maxCoeff());

    // change depth according to the topography
    int nx = mesh.topo_x.size(), ny = mesh.topo_y.size();
    for(int ilat = 0; ilat < nlat; ilat ++){
    for(int ilon = 0; ilon < nlon; ilon ++){
        float altitude;
        float lat0 = model_lat[ilat];
        float lon0 = model_lon[ilon];
        altitude = interp2d(mesh.topo_x.data(),mesh.topo_y.data(),
                            mesh.topo_z.data(),nx,ny,lon0,lat0);

        if(shift_depth) {
            // shift depth according to models 
            float zshift = -depth(0,ilat,ilon) - altitude;
            for(int iz = 0; iz < nz; iz ++){
                depth(iz,ilat,ilon) += zshift;
            }
        }
        else{
            depth(0,ilat,ilon) = -altitude;
        }

        
    }}
}