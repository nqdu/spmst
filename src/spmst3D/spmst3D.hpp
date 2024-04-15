#ifndef _SPMST3D_CLASS_H
#define _SPMST3D_CLASS_H

#include "spm2d/spm2d.hpp"
#include "invparam.hpp"

struct StationPair
{
    int nreceivers;
    float evlon,evlat; // event coordinates
    fvec stlon,stlat; // receiver coodinates
    fvec ttime; // observed travel time

    // data type
    std::string swdtp;
    int period_id; // period index
    int counter; // the global data index
};

class SPMST3D {
public:
    std::vector<StationPair> stapairs; //  event-stations pairs
    fvec tobs; // travel time data in total
    SPM2DMesh mesh; // spm mesh
    InverseParamsBase param; // parameters for 
    bool is_spherical; // if use spherical coordinates
    bool shift_depth; // shift depth if required

    // period vector
    dvec tRc,tRg,tLc,tLg;

public:
    int nlon,nlat,nz;
    std::vector<float> model_lon,model_lat;
    fmat3 depth; // shape(nz,nlat,nlon)
    fmat3 velinit,veltrue; // shape (nz,nlat,nlon)

public:
    SPMST3D() {};
    // IO 
    void read_obsdata(const char *filename);
    void read_model(const char *filename,fmat3 &veloc_in,bool init_mesh);
    void read_topography(const char *filename);
    void read_params(const char *filename);
    void write_traveltime(const fvec &tsyn,const char *filename) const;
    void save_model(const char *filename,const fmat3 &veloc) const;

    // functions
    void synthetic(const fmat3 &vel,fvec &tsyn) const;
    int compute_frechet(const fmat3 &vel,fvec &tsyn,const char *outfile) const;
    void tomography_iter(int iter,fmat3 &vel,fvec &tsyn) const;

private:
    // useful distance function
    float compute_distance(float x0,float y0,float x1,float y1) const;
    int get_period_index(int period_idx,const std::string &swdtp) const;
    void empirical_relation(const float &vs,float &vp, float &rho) const;
    void empirical_deriv(float vp,float vs,float &drda,float &dadb) const;

};

#endif