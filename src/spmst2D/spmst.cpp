#include "spmst2D/spmst2D.hpp"
#include "shared/bilinear.hpp"
#include "shared/csr_matrix.hpp"
#include "shared/parallel_tools.hpp"
#include <fstream>
#include <omp.h>

/**
 * @brief compute synthetic data for current model
 * 
 * @param vel current model, shape(nlat,nlon)
 * @param tsyn synthetic travel time
 */
void SPMST2D:: 
synthetic(const fmat2 &vel,fvec &tsyn) const
{
    // interpolate velocity to spm2d mesh
    fmat2 veloc(mesh.nelmnts,NPT2);
    mesh.interp_velocity(model_lon.data(),model_lat.data(),vel.data(),
                         model_lon.size(),model_lat.size(),veloc);

    // set data index
    std::vector<int> data_idx(nevents);
    data_idx[0] = 0;
    for(int ievt = 1; ievt < nevents; ievt++){
        data_idx[ievt] = data_idx[ievt-1] + nrecvs_per_event[ievt-1];
    }

    // parallel 
    #pragma omp parallel for shared(data_idx,veloc)
    for(int ievt = 0; ievt < nevents; ievt ++){
        SPM2DSolver sol(mesh,veloc);
        sol.locate_source_stations(evlon[ievt],evlat[ievt],&stalon(ievt,0),
                                     &stalat(ievt,0),nrecvs_per_event[ievt],mesh);
        sol.compute_traveltime(mesh,veloc);
        int num = data_idx[ievt];
        for(int ir = 0; ir < nrecvs_per_event[ievt]; ir++){
            tsyn[num + ir] = sol.ttime_recv[ir];
        }
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
int SPMST2D::
compute_frechet(const fmat2 &vel,fvec &tsyn,const char *outfile) const
{    
    // interpolate velocity to spm2d mesh
    fmat2 veloc(mesh.nelmnts,NPT2);
    mesh.interp_velocity(model_lon.data(),model_lat.data(),vel.data(),
                         model_lon.size(),model_lat.size(),veloc);

    // set data index
    std::vector<int> data_idx(nevents);
    data_idx[0] = 0;
    for(int ievt = 1; ievt < nevents; ievt++){
        data_idx[ievt] = data_idx[ievt-1] + nrecvs_per_event[ievt-1];
    }

    // get # of threads
    int nprocs = 1;
    #pragma omp parallel 
    {
        nprocs = omp_get_num_threads();
    }
    ivec nonzeros(nprocs); nonzeros.setZero();

    #pragma omp parallel for shared(data_idx,nonzeros,veloc)
    for(int myrank = 0; myrank < nprocs; myrank ++){
        fvec fdm0(mesh.nptstot);
        int nlat = model_lat.size(), nlon = model_lon.size();
        fmat2 fdm(nlat,nlon);

        // open outfile to write out frechet kernel
        std::string filename = std::string(outfile) + "." + std::to_string(myrank);
        std::ofstream fp; fp.open(filename,std::ios::binary);
        int model_dim = nlat * nlon;

        // write dimension
        int nt = tobs.size();
        fp.write((char*)&model_dim,sizeof(int));
        fp.write((char*)&nt,sizeof(int));

        // allocate jobs to each rank
        int start,end;
        allocate_tasks(nevents,nprocs,myrank,start,end);

        // loop every event in this proc
        for(int ievt = start; ievt <= end; ievt ++){
            int num = data_idx[ievt];

            // compute travel time
            SPM2DSolver sol(mesh,veloc);
            sol.locate_source_stations(evlon[ievt],evlat[ievt],&stalon(ievt,0),
                                        &stalat(ievt,0),nrecvs_per_event[ievt],mesh);
            sol.compute_traveltime(mesh,veloc);
            for(int ir = 0; ir < nrecvs_per_event[ievt]; ir++){
                tsyn[num + ir] = sol.ttime_recv[ir];
            }

            // frechet kernel
            for(int ir = 0; ir < nrecvs_per_event[ievt]; ir++){
                sol.frechet_kernel(ir,mesh,veloc,fdm0);
                fdm.setConstant(0);
                for(int inode = 0; inode < mesh.nptstot; inode ++){
                    float x = mesh.xstore[inode], y = mesh.ystore[inode];
                    float term = fdm0[inode];
                    if(std::abs(term) > 0.0 ){
                        int iy,ix;
                        float coef[4];
                        bilinear(model_lon.data(),model_lat.data(),nlon,nlat,x,y,ix,iy,coef);
                        for(int i = 0; i < 2; i++){
                        for(int j = 0; j < 2; j++){
                            fdm(iy+i,ix+j) += term * coef[i*2+j];
                        }}
                    }
                }

                // allocate temp arrays for kernel writing
                int nar = (fdm.abs() > 0.0).cast<int>().sum();
                std::vector<int> myidx(nar);
                std::vector<float> value(nar);

                // write rowidx/nonzeros for this row
                myidx[0] = num + ir;
                fp.write((char*)&myidx[0],sizeof(int));
                fp.write((char*)&nar,sizeof(int));
                nonzeros[myrank] += nar;

                // save nonzeros to myidx/value
                int c = 0;
                for(int i = 0; i < nlat; i++){
                for(int j = 0; j < nlon; j++){
                    if(std::abs(fdm(i,j)) > 0.0){
                        myidx[c] = i * nlon + j;
                        value[c] = fdm(i,j);
                        c += 1;
                    }
                }}

                // write to binary file
                fp.write((char*)myidx.data(),nar * sizeof(int));
                fp.write((char*)value.data(),sizeof(float) * nar);
            }
        }

        // close file;
        fp.close();
    }

    // merge all files into one file
    csr_matrix::merge_csr_files(nprocs,outfile);

    return nonzeros.sum();
}

/**
 * @brief add Tihkonov 2-nd regularization terms for 2-D problem
 * 
 * @param smat csr sparse matrix
 * @param weight regularization term
 * @param nlat # of nodes in x
 * @param nlon  # of nodes in y
 */
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

/**
 * @brief iteration for one step
 * 
 * @param iter step index
 * @param vel velocity model, shape(nlat,nlon)
 * @param tsyn synthetic data
 */
void SPMST2D:: 
tomography_iter(int iter,fmat2 &vel,fvec &tsyn) const
{
    // compute travel time and save frechet kernel
    const char *frechet_file = "frechet.bin";
    int nar =  compute_frechet(vel,tsyn,frechet_file);

    // create csr matrix
    int nlat = model_lat.size(), nlon = model_lon.size();
    int nt = tobs.size();
    csr_matrix smat(nt+nlat*nlon,nlat*nlon,nar + nlat*nlon*5);

    // read csr matrix from frechet.bin
    std::ifstream fp(frechet_file,std::ios::binary);
    fp.read((char*)&smat.indptr[0],sizeof(int));
    fp.read((char*)&smat.indptr[0],sizeof(int));
    smat.indptr[0] = 0;
    for(int i = 0; i < nt; i ++) {
        int col,nar1;
        fp.read((char*)&col,sizeof(int));
        fp.read((char*)&nar1,sizeof(int));
        int start = smat.indptr[col];
        smat.indptr[col+1] = start +  nar1;
        fp.read((char*)(smat.indices + start),sizeof(int)*nar1);
        fp.read((char*)(smat.data + start),sizeof(float)*nar1);
    }
    fp.close();
    for(int i=nt;i<nt+nlat*nlon;i++) smat.indptr[i+1] = smat.indptr[nt];

    // delete kernel file
    std::remove("frechet.bin");

    // add regularization term
    add_regularization_terms(smat,param.smooth,nlat,nlon);

    // lsqr solver
    printf("solving linear systems by LSQR ...\n");
    fvec res(smat.rows()); res.setConstant(0.0); 
    res.segment(0,nt) = tobs - tsyn;
    fmat2 dvel(nlat,nlon);
    LSQRDict dict(nlon*nlat*2,param.damp,param.smooth);
    smat.lsqr_solver(res.data(),dvel.data(),dict);
    printf("min and max variations: %f %f\n",dvel.minCoeff(),dvel.maxCoeff());

    // update model
    for(int ilat = 0; ilat < nlat; ilat ++){
    for(int ilon = 0; ilon < nlon; ilon ++){
        float v0 = vel(ilat,ilon) + dvel(ilat,ilon);
        if(v0 > param.maxvel) v0 = param.maxvel;
        if(v0 < param.minvel) v0 = param.minvel;
        vel(ilat,ilon) = v0;
    }}
}
