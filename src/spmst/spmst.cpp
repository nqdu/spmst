#include "spmst.hpp"
#include "shared/bilinear.hpp"
#include "shared/csr_matrix.hpp"
#include <fstream>
#include <omp.h>

const int USER_OMP_NUM_THREADS = 8;

/**
 * @brief compute synthetic data for current model
 * 
 * @param vel current model, shape(nlat,nlon)
 * @param tsyn synthetic travel time
 */
void SPMST:: 
synthetic(const fmat2 &vel,fvec &tsyn)
{
    // interpolate velocity to spm2d mesh
    set_velocity(vel,spm2dbase);

    // set data index
    std::vector<int> data_idx(nevents);
    data_idx[0] = 0;
    for(int ievt = 1; ievt < nevents; ievt++){
        data_idx[ievt] = data_idx[ievt-1] + nrecvs_per_event[ievt-1];
    }

    omp_set_num_threads(USER_OMP_NUM_THREADS);
    #pragma omp parallel for shared(data_idx,spm2dbase)
    for(int ievt = 0; ievt < nevents; ievt ++){
        SPM2D spm2d = spm2dbase;
        spm2d.locate_source_stations(evlon[ievt],evlat[ievt],&stalon(ievt,0),
                                     &stalat(ievt,0),nrecvs_per_event[ievt]);
        spm2d.compute_traveltime();
        int num = data_idx[ievt];
        for(int ir = 0; ir < nrecvs_per_event[ievt]; ir++){
            tsyn[num + ir] = spm2d.ttime_recv[ir];
        }
    }
}

/**
 * @brief merge frechet kernel for all procs into a single file
 * 
 * @param outfile 
 */
static void 
merge_frechet_kernel(const char *outfile)
{
    char filename[256];
    const int nprocs = USER_OMP_NUM_THREADS;
    if(nprocs != 1){
        FILE *fp,*fpout;
        char line[128];
        fpout = fopen(outfile,"w");
        for(int irank = 0; irank < nprocs; irank ++){
            sprintf(filename,"%s.%d.txt",outfile,irank);
            fp = fopen(filename,"r");
            while(fgets(line,sizeof(line),fp) != NULL){
                fputs(line,fpout);
            }
            fclose(fp);
            
            // remove this file
            std::remove(filename);
        }
        fclose(fpout);
    }
    else{
        sprintf(filename,"%s.%d.txt",outfile,0);
        std::rename(filename,outfile);
    }
}

/**
 * @brief compute travel time and frechet kernel
 * 
 * @param vel current model, shape(nlat,nlon)
 * @param tsyn synthetic data, shape(nt)
 * @param outfile file to write frechet kernel
 * @return int # of nonzeros for frechet matrix
 */
std::vector<int> SPMST::
compute_frechet(const fmat2 &vel,fvec &tsyn,const char *outfile)
{    
    // interpolate velocity to spm2d mesh
    set_velocity(vel,spm2dbase);

    // set data index
    std::vector<int> data_idx(nevents);
    data_idx[0] = 0;
    for(int ievt = 1; ievt < nevents; ievt++){
        data_idx[ievt] = data_idx[ievt-1] + nrecvs_per_event[ievt-1];
    }

    // for bilinear interpolation
    float lon[nlon],lat[nlat];
    for(int i = 0; i < nlon; i++){
        lon[i] = lonmin + (lonmax - lonmin) / (nlon - 1) * i;
    }
    for(int i = 0; i < nlat; i++){
        lat[i] = latmin + (latmax - latmin) / (nlat - 1) * i;
    }

    // # nonzeros in 
    int nt = tobs.size();
    std::vector<int> nonzeros(nt);

    omp_set_num_threads(USER_OMP_NUM_THREADS);
    #pragma omp parallel for shared(data_idx,spm2dbase,lon,lat,nonzeros)
    for(int myrank = 0; myrank < USER_OMP_NUM_THREADS; myrank ++){
        fvec fdm0(spm2dbase.nptstot);
        fmat2 fdm(nlat,nlon);
        SPM2D spm2d = spm2dbase;

        // open outfile to write out frechet kernel
        char filename[100];
        sprintf(filename,"%s.%d.txt",outfile,myrank);
        FILE *fp = fopen(filename,"w"); 

        // allocate jobs to each rank
        int sub_n = nevents / USER_OMP_NUM_THREADS;
        int num_larger_procs = nevents - USER_OMP_NUM_THREADS * sub_n;
        int start,end;
        if (myrank < num_larger_procs){ sub_n = sub_n + 1;
            start = 0 + myrank * sub_n;
        }
        else if (sub_n > 0){ 
            start = 0 + num_larger_procs + myrank * sub_n;
        }
        else { // this process has only zero elements
            start = -1;
            sub_n = 0;
        }
        end = start + sub_n - 1;

        // loop every event in this proc
        for(int ievt = start; ievt <= end; ievt ++){
            int num = data_idx[ievt];

            // compute travel time
            spm2d.locate_source_stations(evlon[ievt],evlat[ievt],&stalon(ievt,0),
                                        &stalat(ievt,0),nrecvs_per_event[ievt]);
            spm2d.compute_traveltime();
            for(int ir = 0; ir < nrecvs_per_event[ievt]; ir++){
                tsyn[num + ir] = spm2d.ttime_recv[ir];
            }

            // frechet kernel
            for(int ir = 0; ir < nrecvs_per_event[ievt]; ir++){
                spm2d.frechet_kernel(ir,fdm0);
                //printf("%f\n",fdm0.sum());
                fdm.setConstant(0);
                for(int inode = 0; inode < spm2d.nptstot; inode ++){
                    float x = spm2d.xstore[inode], y = spm2d.ystore[inode];
                    float term = fdm0[inode];
                    if(std::abs(term) > 0.0 ){
                        int iy,ix;
                        float coef[4];
                        bilinear(lon,lat,nlon,nlat,x,y,ix,iy,coef);
                        for(int i = 0; i < 2; i++){
                        for(int j = 0; j < 2; j++){
                            fdm(iy+i,ix+j) += term * coef[i*2+j];
                        }}
                    }
                }

                // write out
                int nar = (fdm.abs() > 0.0).cast<int>().sum();
                fprintf(fp,"# %d %d\n",num+ir,nar);
                for(int i = 0; i < nlat; i++){
                for(int j = 0; j< nlon; j++){
                    if(std::abs(fdm(i,j)) > 0.0){
                        fprintf(fp,"%d %g\n",i*nlon+j,fdm(i,j));
                    }
                }}  
                nonzeros[num+ir] = nar;
            }
        }

        // close file;
        fclose(fp);
    }

    // merge all files into one file
    merge_frechet_kernel(outfile);

    return nonzeros;
}

static void 
add_regularization_terms(csr_matrix &smat,float weight,int nlat,int nlon)
{
    // get parameters required
    int n = nlat * nlon; // model dimension
    int nar = smat.nonzeros - n * 5; // nonzeros excluding smooth term
    int m = smat.rows() - n; // data dimension

    int count = 0;
    for(int j = 0; j< nlat; j++){
    for(int i = 0; i< nlon; i++){
        if( i==0 || i == nlon-1 || j==0 || j == nlat-1){
            
            // and more restrictions to boundary points
            if(nar + 1 > smat.nonzeros){
                printf("please increase sparse ratio!\n");
                exit(1);
            }
            int rwc = count + m;
            smat.data[nar] = 2.0 * weight;
            smat.indices[nar] = j * nlon + i;
            smat.indptr[rwc + 1] = smat.indptr[rwc] + 1;
            nar += 1;
            count += 1;
            
           continue;
        }
        else{
            if(nar  + 5 > smat.nonzeros){
                printf("please increase sparse ratio!\n");
                exit(1);
            }
            int rwc = count + m;  // current row
            smat.indptr[rwc +1] = smat.indptr[rwc] + 5;
            int clc = j * nlon + i;// current column
            smat.data[nar] = 4.0 * weight;
            smat.indices[nar] = clc;

            // x direction
            smat.data[nar + 1] = -weight;
            smat.indices[nar + 1] = clc - 1;
            smat.data[nar + 2] = -weight;
            smat.indices[nar + 2] = clc + 1;

            // y direction
            smat.data[nar + 3] = -weight;
            smat.indices[nar + 3] = clc - nlon;
            smat.data[nar + 4] = -weight;
            smat.indices[nar + 4] = clc + nlon;

            nar += 5;
            count += 1;
        }
    }}
}

void SPMST:: 
tomography_iter(int iter,fmat2 &vel,fvec &tsyn)
{
    // compute travel time and save frechet kernel
    const char *frechet_file = "frechet.out";
    std::vector<int> nonzeros =  compute_frechet(vel,tsyn,frechet_file);

    // create csr matrix
    int nar = 0, nt = tobs.size();
    for(int i = 0; i < nt; i++) nar += nonzeros[i];

    csr_matrix smat(nt+nlat*nlon,nlat*nlon,nar + nlat*nlon*5);
    smat.indptr[0] = 0;
    for(int i=0;i<nt;i++){
        smat.indptr[i+1] = smat.indptr[i] +  nonzeros[i];
    }
    for(int i=nt;i<nt+nlat*nlon;i++) smat.indptr[i+1] = smat.indptr[nt];

    // read csr matrix from frechet.out
    smat.read(frechet_file);
    add_regularization_terms(smat,smooth,nlat,nlon);

    // lsqr solver
    printf("solving linear systems by LSQR ...\n");
    fvec res(smat.rows()); res.setConstant(0.0); 
    res.segment(0,nt) = tobs - tsyn;
    fmat2 dvel(nlat,nlon);
    LSQRDict dict(nlon*nlat*2,damp,smooth);
    smat.lsqr_solver(res.data(),dvel.data(),dict);
    printf("min and max variations: %f %f\n",dvel.minCoeff(),dvel.maxCoeff());

    // update model
    for(int ilat = 0; ilat < nlat; ilat ++){
    for(int ilon = 0; ilon < nlon; ilon ++){
        if(dvel(ilat,ilon) > 0.5) dvel(ilat,ilon) = 0.5;
        if(dvel(ilat,ilon) < -0.5) dvel(ilat,ilon) = -0.5;
        float v0 = vel(ilat,ilon) + dvel(ilat,ilon);
        if(v0 > vmax) v0 = vmax;
        if(v0 < vmin) v0 = vmin;
        vel(ilat,ilon) = v0;
    }}
}
