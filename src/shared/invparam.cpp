#include "invparam.hpp"
#include "shared/IO.hpp"
#include <sstream>

void InverseParamsBase :: 
read_file(const char* paramfile) {
    // open file
    std::ifstream infile; infile.open(paramfile);

    // read parameters from file
    read_par_regex("NITERS",maxiter,infile);
    int ierr = read_par_regex("ITER_CURRENT",iter_cur,infile);
    if(ierr == 1) {
        iter_cur = 0;
    }

    // constraints
    read_par_regex("MIN_VELOC",minvel,infile);
    read_par_regex("MAX_VELOC",maxvel,infile);

    // read inv params based on inv_method
    read_par_regex("SMOOTH",smooth,infile);
    read_par_regex("DAMP",damp,infile);
    read_par_regex("NTHREADS",nthreads,infile);

    // synthetic test 
    read_par_regex("SYN_TEST",ifsyn,infile);
    if(ifsyn) {
        read_par_regex("NOISE_LEVEL",noiselevel,infile);
    }

    // print on the screen 
    printf("Inversion Parameters:\n");
    printf("===================================\n");
    printf("Min amd max velocity(km/s) = %f, %f\n",minvel,maxvel);
    printf("Max iterations = %d\n",maxiter);
    printf("current model = %d\n",iter_cur);

    printf("use LSQR solver: ");
    printf("Number of Threads Used = %d\n",nthreads);
    printf("smooth = %f, damp = %f\n",smooth,damp);

    // close file
    infile.close();
    
}