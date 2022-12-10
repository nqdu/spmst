#include "spmst.hpp"
#include <dirent.h>
#include <sys/stat.h>
#include <string>

static void 
create_directory(const char *dirname);

static void 
save_model(const char *filename, const SPMST &tomo,const fmat2 &vel);

int main(){
    // read parameters and init tomo
    SPMST tomo;
    tomo.read_spmst_params("spmst.in");
    tomo.read_obsdata("surfdata.txt");
    tomo.read_velocity("velocinit.in",tomo.velinit);
    if(tomo.do_synthetic){
        tomo.read_velocity("veloctrue.in",tomo.veltrue);
    }
    tomo.read_topography("topo.dat");

    // travel time
    int nt = tomo.tobs.size();
    fvec tsyn(nt);

    // checker board test if required
    if(tomo.do_synthetic){
        printf("\nCheckerboard Resolution Test Begin ...\n");
        tomo.synthetic(tomo.veltrue,tsyn);
        tomo.tobs = tsyn;
    }

    //tomo.read_topography("topo1.dat");

    // init model
    fmat2 vel = tomo.velinit;
    create_directory("models");
    std::string outfile = "models/mod_iter"+ std::to_string(0) +  ".dat";
    save_model(outfile.data(),tomo,vel);
    tomo.write_traveltime(tomo.tobs,"models/disper_obs.dat");

    // now do inversion
    for(int iter = 0; iter < tomo.niters; iter ++){
        printf("\nIteration %d ...\n",iter + 1);
        printf("computing frechet kernel ...\n");
        tomo.tomography_iter(iter,vel,tsyn);

        // rms
        float rms = (tsyn - tomo.tobs).square().sum() / nt;
        printf("rms for model %d: %g\n",iter,rms);

        // save current model
        outfile = "models/mod_iter"+ std::to_string(iter+1) +  ".dat";
        save_model(outfile.data(),tomo,vel);
        outfile = "models/disper_iter" + std::to_string(iter) + ".dat";
        tomo.write_traveltime(tsyn,outfile.data());
    }
    
    // compute travel time for the last iteration
    printf("\nsynthetic traveltime for the final model ...\n");
    tomo.synthetic(vel,tsyn);
    float rms = (tsyn - tomo.tobs).square().sum() / nt;
    printf("rms for the final model: %g\n",rms);
    outfile = "models/disper_iter" + std::to_string(tomo.niters) + ".dat";
    tomo.write_traveltime(tsyn,outfile.data());

    return 0;
}

static void 
create_directory(const char *dirname)
{
    if(opendir(dirname) == NULL){
        mkdir(dirname,0777);
    } 
}

static void 
save_model(const char *filename, const SPMST &tomo,const fmat2 &vel)
{
    FILE *fp = fopen(filename,"w");
    int nlat = tomo.nlat,nlon = tomo.nlon;
    for(int ilat = 0; ilat < nlat; ilat ++){
    for(int ilon = 0; ilon < nlon; ilon ++){
        float lat = tomo.latmin + (tomo.latmax - tomo.latmin) / (nlat-1) * ilat;
        float lon = tomo.lonmin + (tomo.lonmax - tomo.lonmin) / (nlon-1) * ilon;
        fprintf(fp,"%f %f %f\n",lon,lat,vel(ilat,ilon));
    }}
    fclose(fp);
}