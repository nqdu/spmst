#ifndef _SPMST_CLASS_H
#define _SPMST_CLASS_H
#include "spm2d.hpp"

class SPMST {
public:
    // source and receiver
    int nevents;
    fvec evlon,evlat; // shape(nevents)
    std::vector<int> nrecvs_per_event; // shape(nevents)
    fmat2 stalon,stalat; // shape(nevents,# of receivers)
    fvec tobs; // travel time data 
    SPM2D spm2dbase; // forward solver

    bool do_synthetic; // if do synthetic test

public:
    int nlon,nlat;
    float lonmin,lonmax,latmin,latmax;
    float smooth,damp;
    float vmin,vmax;
    int niters;
    fmat2 velinit,veltrue;

public:
    SPMST() {};
    // IO 
    void read_obsdata(const char *filename);
    void read_velocity(const char *filename,fmat2 &veloc_in);
    void set_velocity(const fmat2 &veloc_in,SPM2D &spm2d);
    void read_topography(const char *filename);
    void read_spmst_params(const char *filename);

    // functions
    void synthetic(const fmat2 &vel,fvec &tsyn);
    std::vector<int> compute_frechet(const fmat2 &vel,fvec &tsyn,const char *outfile);
    void tomography_iter(int iter,fmat2 &vel,fvec &tsyn);
};

#endif