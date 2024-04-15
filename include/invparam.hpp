#ifndef _SPMST_INV_PARAM_H
#define _SPMST_INV_PARAM_H

class InverseParamsBase {
public:
    float minvel,maxvel; // velocity prior information
    int maxiter; // maxiteration
    int ifsyn; // if or not do checkerboard test
    int iter_cur; // current iteration number (start from 0)

    // lsmr
    float smooth,damp; // parameters for lsmr
    int nthreads;

    // noise level
    float noiselevel; 

    // read parameter files
    void read_file(const char *filename);
};

#endif 