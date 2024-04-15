#include "spmst3D/spmst3D.hpp"
#include "shared/csr_matrix.hpp"
#include <fstream>

static void 
add_regularization_terms(csr_matrix &smat,float weight,int nz,int nlat,int nlon)
{
    // get parameters required
    int n = nlat * nlon * nz; // model dimension
    int nar = smat.nonzeros - n * 7; // nonzeros excluding smooth term
    int m = smat.rows() - n; // data dimension

    int count = 0;
    for(int k = 0; k < nz; k ++) {
    for(int i = 0; i < nlat; i++){
    for(int j = 0; j < nlon; j++){
    
        int idx = k * nlon * nlat + i * nlon + j;
        if( j ==0 || j == nlon-1 || i==0 || i == nlat-1 || k ==0 || k == nz-1){
            
            // and more restrictions to boundary points
            if(nar + 1 > smat.nonzeros){
                printf("please increase sparse ratio!\n");
                exit(1);
            }
            int rwc = count + m;
            smat.data[nar] = 4.0 * weight;
            //if( k == nz-1) smat.data[nar] = 8.0 * weight;
            smat.indices[nar] = idx;
            smat.indptr[rwc + 1] = smat.indptr[rwc] + 1;
            nar += 1;
            count += 1;
            
           continue;
        }
        else{
            if(nar  + 7 > smat.nonzeros){
                printf("please increase sparse ratio!\n");
                exit(1);
            }
            int rwc = count + m;  // current row
            smat.indptr[rwc +1] = smat.indptr[rwc] + 7;
            int clc = idx;// current column
            smat.data[nar] = 6.0 * weight;
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

            // z direction
            smat.data[nar + 5] = -weight;
            smat.indices[nar + 5] = clc - nlon * nlat;
            smat.data[nar + 6] = -weight;
            smat.indices[nar + 6] = clc + nlon * nlat;

            nar += 7;
            count += 1;
        }
    }}}
    smat.nonzeros = nar;
}

void SPMST3D:: 
tomography_iter(int iter,fmat3 &vel,fvec &tsyn) const
{
    // compute travel time and save frechet kernel
    const char *frechet_file = "frechet.bin";
    int nar = this -> compute_frechet(vel,tsyn,frechet_file);

    // create csr matrix
    int nt = tobs.size();
    int n = nlat * nlon * (nz-1);
    csr_matrix smat(nt+n,n,nar + n * 7);

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
    for(int i=nt;i<nt+n;i++) smat.indptr[i+1] = smat.indptr[nt];

    // read csr matrix from frechet.out and add regularization terms
    add_regularization_terms(smat,param.smooth,nz-1,nlat,nlon);
    std::remove(frechet_file);

    // lsqr solver
    printf("solving linear systems by LSQR ...\n");
    fvec res(smat.rows()); res.setConstant(0.0); 
    res.segment(0,nt) = tobs - tsyn;
    fvec dvel(nlat*nlon*(nz-1));
    LSQRDict dict(dvel.size()*2,param.damp,param.smooth,param.nthreads);
    smat.lsqr_solver(res.data(),dvel.data(),dict);
    printf("min and max variations: %f %f\n",dvel.minCoeff(),dvel.maxCoeff());

    // update model
    for(int iz = 0; iz < nz - 1; iz++){
    for(int ilat = 0; ilat < nlat; ilat ++){
    for(int ilon = 0; ilon < nlon; ilon ++){
        int idx = iz * nlon * nlat + ilat * nlon + ilon;
        float v0 = vel(iz,ilat,ilon) + dvel[idx];
        if(v0 > param.maxvel) v0 = param.maxvel;
        if(v0 < param.minvel) v0 = param.minvel;
        vel(iz,ilat,ilon) = v0;
    }}}
}