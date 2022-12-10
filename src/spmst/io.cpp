#include "spmst.hpp"
#include "shared/gps2dist.hpp"

/**
 * @brief compute two point distance according to different coordinate system
 * 
 * @param x0,y0 (x,y)/(lon,lat) for the first point
 * @param x1,y1  (x,y)/(lon,lat) for the second point
 * @return float 
 */
float SPMST::
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
 * @brief write travel time to file
 * 
 * @param tsyn travel time data vector
 * @param filename output file
 */
void SPMST:: 
write_traveltime(const fvec &tsyn,const char *filename) const 
{
    FILE *fp = fopen(filename,"w");
    if(fp == NULL){
        printf("cannot open %s\n",filename);
        exit(1);
    }
    int inum = 0;
    for(int ievt = 0; ievt < nevents; ievt ++){
        int nsta = nrecvs_per_event[ievt];
        for(int ir = 0; ir < nsta; ir ++){
            float dist = compute_distance(evlon[ievt],evlat[ievt],
                                    stalon(ievt,ir),stalat(ievt,ir));
            fprintf(fp,"%g %g %g %g %g %g\n",evlon[ievt],evlat[ievt],
                    stalon(ievt,ir),stalat(ievt,ir),
                    dist,dist/tsyn[inum+ir]);
        }
        inum += nsta;
    }
    fclose(fp);
}