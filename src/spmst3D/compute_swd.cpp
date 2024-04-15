#include "spmst3D/spmst3D.hpp"
#include "SWD/swd.hpp"
#include "numerical.hpp"
#include "shared/bilinear.hpp"
#include "shared/csr_matrix.hpp"
#include "shared/parallel_tools.hpp"
#include <fstream>
#include <omp.h>

/**
 * @brief Get the period index  for current swdtp and 
 * 
 * @param period_idx 
 * @param swdtp 
 * @return int 
 */
int SPMST3D :: 
get_period_index(int period_idx,const std::string &swdtp) const
{
    int idx;
    int kmaxRc = tRc.size(), kmaxRg = tRg.size();
    int kmaxLc = tLc.size();
    if(swdtp == "Rc"){ 
        idx = period_idx;
    }
    else if(swdtp == "Rg") {
        idx = period_idx + kmaxRc;
    }
    else if(swdtp == "Lc") {
        idx = period_idx + kmaxRc + kmaxRg;
    }
    else {
        idx = period_idx + kmaxRc + kmaxRg + kmaxLc;
    }

    return idx;
}

/**
 * @brief Brocher's relation 
 * 
 * @param vs 
 * @param vp 
 * @param rho 
 */
void SPMST3D:: 
empirical_relation(const float &vs,float &vp, float &rho) const
{
    vp = 0.9409 + 2.0947*vs - 0.8206*std::pow(vs,2)+ 
            0.2683*std::pow(vs,3) - 0.0251*std::pow(vs,4);
    rho = 1.6612 * vp - 0.4721 * std::pow(vp,2) + 
            0.0671 * std::pow(vp,3) - 0.0043 * std::pow(vp,4) + 
            0.000106 * std::pow(vp,5);
}

void SPMST3D::
empirical_deriv(float vp,float vs,float &drda,float &dadb) const 
{
    drda = 1.6612 - 0.4721*2*vp + 0.0671*3*std::pow(vp,2) - 
           0.0043*4*std::pow(vp,3) + 0.000106*5*std::pow(vp,4);
    dadb = 2.0947 - 0.8206*2*vs + 0.2683*3 * std::pow(vs,2)
           - 0.0251*4*std::pow(vs,3);
}

/**
 * Convert kernels of layer-based model to that of grid-based model
 * @param nz no. of grid points for grid-based model
 * @param rmax no. of layers for layer-based model
 * @param dep,vp,vs,rho grid-based model
 * @param rvp,rvs,rrho rthk layer-based model
 * @param rdcdm parameters sensitivity kernel for layer-based model 
 * @param dcdm parameters sensitivity kernel for grid-based model
 */
static void
convert_param(const float *dep,const float *vp,const float *vs,
              const float *rho,int nz,int sublayer,float *rthk,
              float *rvp,float *rvs,float *rrho,
              int rmax,double *rdcdm,double *dcdm)
{
    // check input parameters
    int nr = nz + (nz-1) * sublayer;
    if(nr !=rmax && vs[0]!= 0.0){
        printf("please check input parameters\n");
        exit(1);
    }

    // convert 
    int n=sublayer + 1;
    if(vs[0]!=0.0)for(int i=0;i<nz;i++){ // no water layers
        double tup = 0.0, tdwn = 0.0;
        int idup = (i-1)*n,iddw = i*n;
        dcdm[i] = 0.0;
        for(int j=0;j<n;j++){
            double tmp = (2.*j+1.0)/(2.*n);
            if(idup +j>=0) tup += rdcdm[idup+j] * tmp;
            if(iddw+j<rmax) tdwn += rdcdm[iddw+j]*(1.0 - tmp);
            if(iddw == rmax-1) tdwn += rdcdm[rmax-1];
        }
        dcdm[i] += tup + tdwn;
    }
    else{ // water layer in the top
        dcdm[0] = rdcdm[0];
        for(int i=1;i<nz;i++){
            dcdm[i] = 0.0;
            double tup = 0.0, tdwn = 0.0;
            int idup = (i-2)*n+1,iddw = (i-1)*n+1;
            for(int j=0;j<n;j++){
                double tmp = (2.*j+1.0)/(2.*n);
                if(idup +j>=1) tup += rdcdm[idup+j] * tmp;
                if(iddw+j<rmax) tdwn += rdcdm[iddw+j]*(1.0 - tmp);
                if(iddw == rmax-1) tdwn += rdcdm[rmax-1];
            }
            dcdm[i] += tup + tdwn;
        }
    }
}

/**
 * refine grid based model to layerd based model, water layers considered
 * 
 * @param dep,vp,vs,rho grid-based model
 * @param nz no. of points in grid-based model
 * @param sublayer insert sublayer points to construct layered model
 * @param rdep, rvp, rvs, rrho, rthk layered velocity model 
 * @return rmax no. of layers of refined model
 */
static int 
grid2LayerModel(const float *dep,const float *vp,const float *vs,
                const float *rho,int nz, int sublayers,float *rthk,
                float *rvp,float *rvs, float *rrho)
{
    float thk;
    int n = sublayers + 1;
    int rmax,istart;
    if(vs[0] == 0.0){ // tackle water layer
        rmax = nz-1 + (nz-2) * sublayers + 1;
        istart = 1;
        rthk[0] = dep[1] - dep[0];
        rrho[0] = rho[0];
        rvs[0] = 0.0;
        rvp[0] = vp[0];
    }
    else{
        rmax = nz + (nz-1) * sublayers;
        istart = 0;
    }

    int k = istart;
    for(int i=istart;i<nz-1;i++){
        thk = (dep[i+1] - dep[i]) / n;
        // parameters in each layer are midpoint values in grid-based model
        for(int j=0;j<n;j++){ 
            rthk[k] = thk;
            rvp[k] = vp[i] + (2*j+1)*(vp[i+1]-vp[i])/(2.0*n);
            rvs[k] = vs[i] + (2*j+1)*(vs[i+1]-vs[i])/(2.0*n);
            rrho[k] = rho[i] + (2*j+1)*(rho[i+1]-rho[i])/(2.0*n);
            k++ ;
        }
    }

    // half space
    k = rmax-1;
    rthk[k] = 0.0;
    rvp[k] = vp[nz-1];
    rvs[k] = vs[nz-1];
    rrho[k] = rho[nz-1];

    return rmax;
}


/**
 * @brief internal function to compute dispersion
 * 
 * @param vs,vp,rho elastic parameters, shape(nz)
 * @param dep depth, shape(nz), in km
 * @param nz 
 * @param sublayer # how many layers to insert into the original layer
 * @param tRc,tRg,tLc,tLg period vector for each wave type
 * @param vc,vout dispersion curve for phase/used velocity, some of elements vg = vc 
 */
static void 
disp1D(const float *vs,const float *vp,const float *rho,
      const float *dep,int nz,int sublayer,const dvec &tRc,const dvec &tRg,
      const dvec &tLc,const dvec &tLg,double *vc,double *vout,bool is_sph_model)
{
    // convert grid-based model to layer-based model
    int rmax = nz + (nz-1) * sublayer;
    float *rvp = new float [rmax],*rvs = new float [rmax];
    float *rthk = new float [rmax],*rrho = new float [rmax];
    rmax = grid2LayerModel(dep,vp,vs,rho,nz,sublayer,rthk,rvp,rvs,rrho);

    // prepare parameters
    int kmaxRc = tRc.size(),kmaxRg = tRg.size();
    int kmaxLc = tLc.size(), kmaxLg = tLg.size();
    bool sphere = is_sph_model, keep_flat = false;
    int mode = 0;

    // compute dispersion
    if(kmaxRc > 0){
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRc.data(),vc,
                kmaxRc,"Rc",mode,sphere,keep_flat);
        memcpy(vout,vc,kmaxRc*sizeof(double));
    }
    if(kmaxRg > 0){
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRg.data(),vc + kmaxRc,
                kmaxRg,"Rc",mode,sphere,keep_flat);
        groupvel_r(rthk,rvp,rvs,rrho,rmax,tRg.data(),
                  vout + kmaxRc,kmaxRg,mode,sphere);
    }
    if(kmaxLc > 0){
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLc.data(),vc + kmaxRc+kmaxRg,
                 kmaxLc,"Lc",mode,sphere,keep_flat);
        memcpy(vout+ kmaxRc+kmaxRg,vc+ kmaxRc+kmaxRg,kmaxLc*sizeof(double));
    }
    if(kmaxLg > 0){
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLg.data(),vc +kmaxRc+kmaxRg+kmaxLc,
                 kmaxLg,"Lc",mode,sphere,keep_flat);
        groupvel_l(rthk,rvs,rrho,rmax,tLg.data(),vout+kmaxRc+kmaxRg+kmaxLc,
                    kmaxLg,mode,sphere);
    }

    // free space
    delete[] rrho;
    delete[] rvp;
    delete[] rvs;
    delete[] rthk;
}

/**
 * @brief internal function to compute dispersion and sensitivity kernel
 * 
 * @param vs,vp,rho elastic parameters, shape(nz)
 * @param dep depth, shape(nz), in km
 * @param nz 
 * @param sublayer # how many layers to insert into the original layer
 * @param tRc,tRg,tLc,tLg period vector for each wave type
 * @param vdis,vdiso  phase/used dispersion curve 
 * @param kvs,kvp,krho sensitivity kernels
 */
static void 
Kernel1D(const float *vs,const float *vp,const float *rho,
        const float *dep,int nz,int sublayer,const dvec &tRc,const dvec &tRg,
        const dvec &tLc, const dvec &tLg,double *vdis,double *vdiso,double *kvs,double *kvp,
        double *krho,bool is_sph_model)
{
    // convert grid-based model to layer-based model
    int rmax = nz + (nz-1) * sublayer;
    float *rvp = new float [rmax],*rvs = new float [rmax];
    float *rthk = new float [rmax],*rrho = new float [rmax];
    rmax = grid2LayerModel(dep,vp,vs,rho,nz,sublayer,rthk,rvp,rvs,rrho);

    // allocate space for kernel
    double *ur = new double [rmax], *uz = new double[rmax]; 
    double *tr = new double [rmax], *tz = new double[rmax];
    double *dudar = new double[rmax], *dudbr = new double[rmax];
    double *dudrr = new double [rmax],*dudhr = new double [rmax];
    double *dcdar = new double[rmax], *dcdbr = new double[rmax];
    double *dcdrr = new double [rmax],*dcdhr = new double [rmax];

    // prepare parameters
    int kmaxRc = tRc.size(),kmaxRg = tRg.size();
    int kmaxLc= tLc.size(), kmaxLg = tLg.size();
    bool sphere = is_sph_model;
    int mode = 0, iflsph = is_sph_model == true;

    // compute dispersion and sensitivity kernel
    if(kmaxRc > 0){ // rayleigh group
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRc.data(),vdis,
                kmaxRc,"Rc",mode,sphere);
        for(int i=0;i<kmaxRc;i++){
            double cg;
            sregn96_(rthk,rvp,rvs,rrho,rmax,&tRc[i],vdis+i,&cg,
                    ur,uz,tr,tz,dcdar,dcdbr,dcdhr,dcdrr,iflsph);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdar,kvp+i*nz);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdbr,kvs+i*nz);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdrr,krho+i*nz);
            vdiso[i] = vdis[i];
        }   
    }
    if(kmaxRg > 0){
        double cp[kmaxRg],t1[kmaxRg],t2[kmaxRg],c1[kmaxRg],c2[kmaxRg];
        double dt = 0.01;
        for(int i=0;i<kmaxRg;i++){
            t1[i] = tRg[i] * (1.0 + 0.5 * dt);
            t2[i] = tRg[i] * (1.0 - 0.5 * dt);
        }
        surfdisp(rthk,rvp,rvs,rrho,rmax,tRg.data(),cp,kmaxRg,"Rc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t1,c1,kmaxRg,"Rc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t2,c2,kmaxRg,"Rc",mode,sphere);
        for(int i=0;i<kmaxRg;i++){
            int k = i + kmaxRc;
            sregnpu_(rthk,rvp,rvs,rrho,rmax,&tRg[i],cp+i,vdis+k,
                     ur,uz,tr,tz,t1+i,c1+i,t2+i,c2+i,dcdar,
                     dcdbr,dcdhr,dcdrr,dudar,dudbr,dudhr,dudrr,iflsph);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudar,kvp+k*nz);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudbr,kvs+k*nz);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudrr,krho+k*nz);
            vdiso[k] = vdis[k];
            vdis[k] = cp[i];        
        }
    }
    if(kmaxLc > 0){
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLc.data(),vdis+kmaxRc+kmaxRg,
                kmaxLc,"Lc",mode,sphere);
        for(int i=0;i<kmaxLc;i++){
            double cg;
            int k = i + kmaxRc + kmaxRg;
            slegn96_(rthk,rvs,rrho,rmax,&tLc[i],vdis+k,
                    &cg,ur,tr,dcdbr,dcdhr,dcdrr,iflsph);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdbr,kvs+k*nz);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dcdrr,krho+k*nz);

            // copy
            vdiso[k] = vdis[k];          
        }  
    }
    if(kmaxLg > 0){
        double cp[kmaxLg],t1[kmaxLg],t2[kmaxLg],c1[kmaxLg],c2[kmaxLg];
        double dt = 0.01;
        for(int i=0;i<kmaxLg;i++){
            t1[i] = tLg[i] * (1.0 + 0.5 * dt);
            t2[i] = tLg[i] * (1.0 - 0.5 * dt);
        }
        surfdisp(rthk,rvp,rvs,rrho,rmax,tLg.data(),cp,
                kmaxLg,"Lc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t1,c1,
                kmaxLg,"Lc",mode,sphere);
        surfdisp(rthk,rvp,rvs,rrho,rmax,t2,c2,
                kmaxLg,"Lc",mode,sphere);
        for(int i=0;i<kmaxLg;i++){
            int k = i + kmaxRc + kmaxRg + kmaxLc;
            slegnpu_(rthk,rvs,rrho,rmax,&tLg[i],cp+i,vdis+k,
                    ur,tr,t1+i,c1+i,t2+i,c2+i,dcdbr,dcdhr,dcdrr,
                    dudbr,dudhr,dudrr,iflsph);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudbr,kvs+k*nz);
            convert_param(dep,vp,vs,rho,nz,sublayer,
                        rthk,rvp,rvs,rrho,rmax,dudrr,krho+k*nz);

            // copy
            vdiso[k] = vdis[k];
            vdis[k] = cp[i];        
        } 
    }

    // free space
    delete[] rrho;delete[] rvp;delete[] rvs;delete[] rthk;
    delete[] ur; delete [] uz; delete [] tr; delete[] tz;
    delete [] dudrr; delete[] dudar; delete[] dudbr; delete[] dudhr;
    delete [] dcdrr; delete[] dcdar; delete[] dcdbr; delete[] dcdhr;
}

/**
 * @brief synthetic travel time for a given model
 * 
 * @param vel 3D model
 * @param tsyn synthetic traveltime 
 */
void SPMST3D:: 
synthetic(const fmat3 &vel,fvec &tsyn)  const
{
    // allocate space for dispersion map
    int nt = tRc.size() + tRg.size() + tLc.size() + tLg.size();
    fmat2 vc(nt,nlat*nlon),vcout(nt,nlat*nlon); // shape(nt,nlat*nlon)

    // step 1: compute swd map
    printf("computing swd ...\n");
    #pragma omp parallel for shared(vc,vcout) collapse(2)
    for(int ilat = 0; ilat < nlat; ilat ++){
    for(int ilon = 0; ilon < nlon; ilon ++){
        int ii = ilat * nlon + ilon;

        // get vp and rho through brocher's relation
        std::vector<float> vs(nz), vp(nz),rho(nz),z(nz);
        dvec cg(nt),cgout(nt);
        for(int iz = 0; iz < nz; iz ++){
            vs[iz] = vel(iz,ilat,ilon);
            z[iz] = depth(iz,ilat,ilon);
            this-> empirical_relation(vs[iz],vp[iz],rho[iz]);
        }

        // compute dispersion map
        const int sublayer = 3;
        disp1D(vs.data(),vp.data(),rho.data(),z.data(),nz,
                sublayer,tRc,tRg,tLc,tLg,cg.data(),cgout.data(),
                is_spherical);
        
        // copy disp curve to vel2d
        vc.col(ii) = cg.cast<float>();
        vcout.col(ii) = cgout.cast<float>();
    }}

    // get num threads
    int nprocs{};
    #pragma omp parallel
    {
        nprocs = omp_get_num_threads();
    }

    // step2 compute travel time
    printf("computing traveltime ...\n");
    int nevents = stapairs.size();
    #pragma omp parallel for
    for(int irank = 0; irank < nprocs; irank ++) {
        int startid,endid;
        allocate_tasks(nevents,nprocs,irank,startid,endid);
        for(int ievt = startid; ievt <= endid; ievt ++){
            const auto &pair = stapairs[ievt];
            int pid = this-> get_period_index(pair.period_id,pair.swdtp);

            // print progress bar
            if(irank == 0) {
                int size = endid - startid + 1;
                float per = (ievt-startid+1.) / size * 100;
                print_progressbar(per);
            }

            // get velocity 
            fmat2 veloc(mesh.nelmnts,NPT2);
            mesh.interp_velocity(model_lon.data(),model_lat.data(),&vc(pid,0),
                                model_lon.size(),model_lat.size(),veloc);
            
            // initialize solver
            SPM2DSolver sol(mesh,veloc);
            sol.locate_source_stations(pair.evlon,pair.evlat,pair.stlon.data(),
                                        pair.stlat.data(),pair.nreceivers,mesh);
            sol.compute_traveltime(mesh,veloc);

            if(pair.swdtp == "Rg" || pair.swdtp == "Lg") { // for group velocity, the ray path is determined by phase veloc
                mesh.interp_velocity(model_lon.data(),model_lat.data(),&vcout(pid,0),
                                    model_lon.size(),model_lat.size(),veloc);
            
                sol.recompute_time(mesh,veloc);
            }
            for(int ir = 0; ir < pair.nreceivers; ir ++){
                tsyn[ir + pair.counter] = sol.ttime_recv[ir];
            }
        }
    }
    printf("\n");
}

/**
 * @brief 
 * 
 * @param vel current 3D model
 * @param tsyn synthetic data
 * @param outfile output file to store the frechet kernel
 * @return std::vector<int> 
 */
int SPMST3D:: 
compute_frechet(const fmat3 &vel,fvec &tsyn,const char *outfile) const
{
    // allocate space for dispersion map
    int nt = tRc.size() + tRg.size() + tLc.size() + tLg.size();
    fmat2 vc(nt,nlat*nlon),vcout(nt,nlat*nlon);
    dmat3 ker(nt,nz,nlat*nlon);

    // first compute dispmap and 1-D sensitivity kernel
    printf("computing Surface Wave 1-D Frechet Kernel ...\n");
    #pragma omp parallel for shared(vc,vcout,ker) collapse(2)
    for(int ilat = 0; ilat < nlat; ilat ++) {
    for(int ilon = 0; ilon < nlon; ilon ++) {
        int ii = ilat * nlon + ilon;

        // get vp and rho through brocher's relation
        std::vector<float> vs(nz), vp(nz),rho(nz),z(nz);
        dvec cg(nt),cgout(nt);
        std::vector<double> krho(nt*nz),kvp(nt*nz),kvs(nt*nz);
        for(int iz = 0; iz < nz; iz ++){
            vs[iz] = vel(iz,ilat,ilon);
            z[iz] = depth(iz,ilat,ilon);
            this-> empirical_relation(vs[iz],vp[iz],rho[iz]);
        }

        // compute dispersion map and sensitivity kerel
        const int sublayer = 3;
        Kernel1D(vs.data(),vp.data(),rho.data(),z.data(),nz,
                 sublayer,tRc,tRg,tLc,tLg,cg.data(),cgout.data(),
                 kvs.data(),kvp.data(),krho.data(),this->is_spherical);
        
        // copy disp curve to vel2d
        for(int it = 0; it < nt; it ++){
            vc(it,ii) = cg[it];
            vcout(it,ii) = cgout[it];
            for(int iz = 0; iz < nz; iz ++){
                float drda,dadb;
                this -> empirical_deriv(vp[iz],vs[iz],drda,dadb);
                ker(it,iz,ii) = kvs[it*nz+iz]+ kvp[it*nz+iz] * dadb 
                                    + krho[it*nz+iz] * drda * dadb;
            }
        }
    }}

    // get num threads
    int nprocs{};
    #pragma omp parallel
    {
        nprocs = omp_get_num_threads();
    }
    ivec nonzeros(nprocs); nonzeros.setZero();

    // frechet kernel for 2D 
    printf("computing traveltimes and 2-D Frechet Kernel ...\n");
    int nevents = stapairs.size();
    #pragma omp parallel for shared(vc,vcout,nonzeros,ker)
    for(int myrank = 0; myrank < nprocs; myrank ++){
        fvec fdm0(mesh.nptstot);
        fmat2 fdm(nlat,nlon);
        fvec frechet(nlat*nlon*(nz-1));

        // open outfile to write out frechet kernel
        std::string filename = std::string(outfile) + "." + std::to_string(myrank);
        std::ofstream fp; fp.open(filename,std::ios::binary);
        int model_dim = nlat * nlon * (nz-1);

        // write dimension
        int m = tobs.size();
        fp.write((char*)&m,sizeof(int));
        fp.write((char*)&model_dim,sizeof(int));

        // allocate jobs to each rank
        int startid,endid;
        allocate_tasks(nevents,nprocs,myrank,startid,endid);

        // loop every event in this proc
        for(int ievt = startid; ievt <= endid; ievt ++){
            auto const &pair = stapairs[ievt];

           // print progress bar
            if(myrank == 0) {
                int size = endid - startid + 1;
                float per = (ievt-startid+1.) / size * 100;
                print_progressbar(per);
            }
            
            // get velocity
            int pid = this-> get_period_index(pair.period_id,pair.swdtp);
            fmat2 veloc(mesh.nelmnts,NPT2);
            mesh.interp_velocity(model_lon.data(),model_lat.data(),&vc(pid,0),
                                model_lon.size(),model_lat.size(),veloc);

            // initialize solver 
            SPM2DSolver sol(mesh,veloc);
            sol.locate_source_stations(pair.evlon,pair.evlat,pair.stlon.data(),
                                        pair.stlat.data(),pair.nreceivers,mesh);
            sol.compute_traveltime(mesh,veloc);

            if(pair.swdtp == "Rg" || pair.swdtp == "Lg") { // for group velocity, the ray path is determined by phase veloc
                mesh.interp_velocity(model_lon.data(),model_lat.data(),&vcout(pid,0),
                                    model_lon.size(),model_lat.size(),veloc);
                sol.recompute_time(mesh,veloc);
            }

            for(int ir = 0; ir < pair.nreceivers; ir ++){
                tsyn[ir + pair.counter] = sol.ttime_recv[ir];
            }
 
            // frechet kernel
            for(int ir = 0; ir < pair.nreceivers; ir++){
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

                // compute 3D frechet kerel
                frechet.setZero();
                for(int iz = 0; iz < nz-1; iz ++){
                for(int ilat = 0; ilat < nlat; ilat ++){
                for(int ilon = 0; ilon < nlon; ilon ++){
                    int idx = iz * nlon * nlat + ilat * nlon + ilon;
                    frechet[idx] = ker(pid,iz,ilat*nlon+ilon) * fdm(ilat,ilon);
                }}}

                // statistics
                int nar = (frechet.abs() > 0.0).cast<int>().sum();
                int counter = pair.counter + ir;
                std::vector<int> myindx(nar);
                std::vector<float> value(nar);
                nonzeros[myrank] += nar;

                // write out flags
                fp.write((char*)&counter,sizeof(int));
                fp.write((char*)&nar,sizeof(int));

                // write data
                counter = 0;
                for(int iz = 0; iz < nz-1; iz ++){
                for(int ilat = 0; ilat < nlat; ilat ++){
                for(int ilon = 0; ilon < nlon; ilon ++){
                    int idx = iz * nlon * nlat + ilat * nlon + ilon;
                    float f = frechet[idx];
                    if(std::abs(f) > 0.0){
                        myindx[counter] = idx;
                        value[counter] = f;
                        counter += 1;
                    }
                }}}

                // write to binary file
                fp.write((char*)myindx.data(),nar * sizeof(int));
                fp.write((char*)value.data(),sizeof(float) * nar);
            }
        }

        // close file;
        fp.close();
    }
    printf("\n");

    // merge all frechet kernels in a file
    csr_matrix::merge_csr_files(nprocs,outfile);

    return nonzeros.sum();

}