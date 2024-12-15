#include "spmst2D/spmst2D.hpp"
#include <string>
#include "shared/IO.hpp"
#include "shared/gaussian.hpp"

int main(int argc, char* argv[]){
   // check input parameters
    std::string paramfile,modfile,datafile,modtrue,topofile;
    printf("\n**************************************\n");
    printf("*** 2-D SPM Surface Wave Tomography *****\n");
    printf("**************************************\n");
    if(argc == 5 || argc == 6){
        paramfile = argv[1];
        modfile = argv[4];
        datafile = argv[2];
        topofile = argv[3];
        modtrue = "None";
        if(argc == 6) modtrue = argv[5];
    }
    else {
        printf("Please run this executable file by:\n");
        printf("./this paramfile datafile topofile initmod (truemod)\n");
        exit(1);
    }
    printf("\n");

    // read parameters and init tomo
    SPMST2D tomo;
    tomo.read_params(paramfile.c_str());
    tomo.read_obsdata(datafile.c_str());
    tomo.read_model(modfile.c_str(),tomo.velinit,true);
    if(tomo.param.ifsyn){
        if(argc !=6) {
            printf("you should input a truemod when enabling SYN_TEST!\n");
            exit(1);
        }
        tomo.read_model(modtrue.c_str(),tomo.veltrue);
    }
    tomo.read_topography(topofile.c_str());

    // parameters used
    const auto &param = tomo.param;

    // travel time
    int nt = tomo.tobs.size();
    fvec tsyn(nt);

    // checker board test if required
    if(tomo.param.ifsyn){
        printf("\nCheckerboard Resolution Test Begin ...\n");
        tomo.synthetic(tomo.veltrue,tsyn);

        // add noise 
        for(int it = 0; it < nt; it ++) {
            tsyn[it] += tomo.param.noiselevel * gaussian(0.0,1.0);
        }
        tomo.tobs = tsyn;
    }

    // init model
    fmat2 vel = tomo.velinit;
    create_directory("models");
    std::string outfile = "models/mod_iter"+ std::to_string(param.iter_cur) +  ".dat";
    tomo.write_model(outfile.data(),vel);
    tomo.write_traveltime(tomo.tobs,"models/disper_obs.dat");

    // now do inversion
    for(int iter = param.iter_cur; iter < param.iter_cur + param.maxiter; iter ++){
        printf("\nIteration %d ...\n",iter + 1);
        printf("computing frechet kernel ...\n");
        tomo.tomography_iter(iter,vel,tsyn);

        // misfit
        float rms = (tsyn - tomo.tobs).square().sum() / nt;
        printf("L2 misfit for model %d: %g\n",iter,rms);

        // save current model
        outfile = "models/mod_iter"+ std::to_string(iter + 1) +  ".dat";
        tomo.write_model(outfile.data(),vel);
        outfile = "models/disper_iter" + std::to_string(iter) + ".dat";
        tomo.write_traveltime(tsyn,outfile.data());
    }
    
    // compute travel time for the last iteration
    printf("\nsynthetic traveltime for the final model ...\n");
    tomo.synthetic(vel,tsyn);
    float rms = (tsyn - tomo.tobs).square().sum() / nt;
    printf("L2 misfit for the final model: %g\n",rms);
    outfile = "models/disper_iter" + std::to_string(param.maxiter + param.iter_cur) + ".dat";
    tomo.write_traveltime(tsyn,outfile.data());

    return 0;
}