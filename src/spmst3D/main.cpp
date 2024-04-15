#include "spmst3D/spmst3D.hpp"
#include <string>
#include "shared/IO.hpp"
#include "shared/gaussian.hpp"

int main(int argc,char** argv) {
    std::string paramfile,modfile,datafile,modtrue,topofile;
    printf("\n**********************************************\n");
    printf("*** 3-D SPM Direct Surface Wave Tomography *****\n");
    printf("*************************************************\n");
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

    // create directory to store results
    const std::string outdir = "models/";
    std::string filename;
    create_directory(outdir.c_str());

    // read parameters and init tomo
    SPMST3D tomo;
    tomo.read_params(paramfile.c_str());
    tomo.read_obsdata(datafile.c_str());
    tomo.read_model(modfile.c_str(),tomo.velinit,true);
    if(tomo.param.ifsyn){
        tomo.read_model(modtrue.c_str(),tomo.veltrue,false);
    }
    tomo.read_topography(topofile.c_str());

    // travel time
    int nt = tomo.tobs.size();
    fvec tsyn(nt);

    // checker board test if required
    const auto &param = tomo.param;
    if(param.ifsyn){
        printf("\nCheckerboard Resolution Test, noise level = %f ...\n",tomo.param.noiselevel);
        tomo.synthetic(tomo.veltrue,tsyn);
        filename = outdir + "mod_true.dat";
        tomo.save_model(filename.c_str(),tomo.veltrue);

        // add gaussian noise
        for(int it = 0; it < nt; it ++){
            tsyn[it] += gaussian(0.,1.) * tomo.param.noiselevel;
        }
        tomo.tobs = tsyn * 1.0f;
    }

    // save obs data
    filename = outdir + "disper_obs.dat";
    tomo.write_traveltime(tomo.tobs,filename.c_str());

    // init model
    fmat3 vel = tomo.velinit;
    filename = outdir + "mod_iter"+ std::to_string(param.iter_cur) +  ".dat";
    tomo.save_model(filename.c_str(),vel);

    // now do inversion
    int max_iter = param.maxiter;
    for(int ii = 0; ii < max_iter; ii ++){
        int iter = param.iter_cur + ii;
        printf("\nIteration %d ...\n",iter + 1);
        tomo.tomography_iter(iter,vel,tsyn);

        // rms
        float rms = (tsyn - tomo.tobs).square().sum() / nt;
        printf("L2 misfit for model %d: %g\n",iter,rms);

        // save current model
        filename = "models/mod_iter"+ std::to_string(iter+1) +  ".dat";
        tomo.save_model(filename.c_str(),vel);
        filename = "models/disper_iter" + std::to_string(iter) + ".dat";
        tomo.write_traveltime(tsyn,filename.data());
    }
    
    // compute travel time for the last iteration
    printf("\nsynthetic traveltime for the final model ...\n");
    tomo.synthetic(vel,tsyn);
    float rms = (tsyn - tomo.tobs).square().sum() / nt;
    printf("L2 misfit for final model: %g\n",rms);
    filename = "models/disper_iter" + std::to_string(max_iter + param.iter_cur) + ".dat";
    tomo.write_traveltime(tsyn,filename.data());

    return 0;
}