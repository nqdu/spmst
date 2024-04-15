#ifndef _SPMST2D_CLASS_H
#define _SPMST2D_CLASS_H

#include "spm2d/spm2d.hpp"
#include "invparam.hpp"

class SPMST2D {
public:
    // source and receiver
    int nevents;
    fvec evlon,evlat; // shape(nevents)
    std::vector<int> nrecvs_per_event; // shape(nevents)
    fmat2 stalon,stalat; // shape(nevents,# of receivers)
    fvec tobs; // travel time data 
    SPM2DMesh mesh; // SPM mesh

    // velocity model
    fmat2 velinit, veltrue; // shape (nlon,nlat)

    // inverse params
    InverseParamsBase param;

private: 
    std::vector<float> model_lon,model_lat;
    bool is_spherical; // if it's spherical coordinates

public:
    SPMST2D() {};
    // IO 
    void read_obsdata(const char *filename);
    void read_model(const char *filename,fmat2 &veloc_in, bool init_mesh = false);
    void set_velocity(const fmat2 &veloc_in);
    void read_topography(const char *filename);
    void read_params(const char *filename);
    void write_traveltime(const fvec &tsyn,const char *filename) const;
    void write_model(const char *filename,const fmat2 &vel) const;

    // functions
    void synthetic(const fmat2 &vel,fvec &tsyn) const;
    int compute_frechet(const fmat2 &vel,fvec &tsyn,const char *outfile) const;
    void tomography_iter(int iter,fmat2 &vel,fvec &tsyn) const;

private:
    // useful distance function
    float compute_distance(float x0,float y0,float x1,float y1) const;
};

#endif