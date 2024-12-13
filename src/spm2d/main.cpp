#include "spm2d/spm2d.hpp"
#include "shared/IO.hpp"
#include "shared/gps2dist.hpp"
#include "shared/bilinear.hpp"
#include <iostream>
#include <chrono>

int main(int argc, char **argv){
    // check input args
    if(argc != 4 && argc != 5){
        printf("Usage: ./syn spm2d.in surfdata.txt topo.txt [veloc.txt]\n");
        exit(1);
    }
    const char *headfile = argv[1], *datafile=argv[2], *topofile = argv[3];
    char * velocfile = NULL;
    if(argc == 5) velocfile = argv[4];

    // read header file
    int nx,nz;
    float zmin,xmin,zmax,xmax;
    float dx,dz;
    int use_sph;
    printf("\nreading mesh file %s...\n",headfile);
    std::ifstream fpin; fpin.open(headfile);
    if(!fpin.is_open()){
        printf("cannot open %s\n",headfile);
        exit(1);
    }
    read_file_param(fpin,xmin,xmax);
    read_file_param(fpin,zmin,zmax);
    read_file_param(fpin,nx,nz);
    read_file_param(fpin,use_sph);
    fpin.close();
    dx = (xmax - xmin) / nx; dz = (zmax - zmin) / nz;
    if(use_sph == 1){
        printf("using spherical coordinates ...\n");
        printf("lonmin = %f lonmax = %f\n",xmin,xmax);
        printf("latmin = %f latmax = %f\n",zmin,zmax);
        printf("nlat = %d nlon = %d\n",nz,nx);
    }
    else{
        printf("using cartesian coordinates ...\n");
        printf("xmin = %f xmax = %f\n",xmin,xmax);
        printf("ymin = %f ymax = %f\n",zmin,zmax);
        printf("ny = %d nx = %d\n",nz,nx);
    }
	printf("NPTX = %d NPTZ = %d\n",NPTX,NPTZ);

    // source and receiver
    float xs,zs; 
    int nr; 
    char dummy;
    fpin.open(datafile);
    if(!fpin.is_open()){
        printf("cannot open %s\n",datafile);
        exit(1);
    }
    std::string line;
    std::getline(fpin,line);
    sscanf(line.data(),"%c%f%f%d",&dummy,&xs,&zs,&nr);
    std::vector<float> xr(nr),zr(nr);
    for(int i = 0; i < nr; i++){
        read_file_param(fpin,xr[i],zr[i]);
    }
    fpin.close();

    // initialize
    SPM2DMesh mesh(xmin,zmin,dx,dz,nx,nz,use_sph==1);
    mesh.create_graph();

    // set topography
    printf("\nreading topography from %s ...\n",topofile);
    mesh.read_topography(topofile);

    // write vtk
    mesh.write_vtk("topo.vtk");
    
    // set velocity
    printf("\nsetting velocity ...\n");
    if(velocfile == NULL){
        printf("veloc file is not given, set velocity to 1.5\n");
        fmat2 veloc(nz,nx);
        std::vector<float> xcord(nx),zcord(nz);
        for(int iz = 0; iz < nz; iz ++) zcord[iz] = zmin + dz * iz;
        for(int ix = 0; ix < nx; ix ++) xcord[ix] = xmin + dx * ix;
        veloc.setConstant(1.5);
        // float radius = std::max(xmax-xmin,zmax-zmin) / 5.0;
        // for(int iz = 0; iz < nz; iz ++){
        //     zcord[iz] = zmin + dz * iz;
        //     for(int ix = 0; ix < nx; ix ++){
        //         xcord[ix] = xmin + dx * ix;
        //         float dist = std::hypot(xcord[ix] - xmin - nx/2. * dx,
        //                                 zcord[iz] - zmin - dz * nz/4.);
                
        //         dist = std::pow(dist /radius,2);
        //         veloc(iz,ix) += 0.5 * std::exp(-dist);
                
        //         dist = std::hypot(xcord[ix] - xmin - nx/2. * dx,
        //                             zcord[iz] - zmin - dz * nz/4. * 3);
        //         dist = std::pow(dist /radius,2);
        //         veloc(iz,ix) -= 0.5 * std::exp(-dist);
        //     }
        // }
        mesh.set_velocity(xcord.data(),zcord.data(),veloc.data(),nx,nz);
    }
    else{
        printf("reading velocity from %s\n",velocfile);
        fpin.open(velocfile);
        if(!fpin.is_open()){
            printf("cannot open %s\n",velocfile);
            exit(1);
        }
        int nnz,nnx;
        float zzmin,zzmax,xxmin,xxmax;
        int use_sph_tmp;
        read_file_param(fpin,xxmin,xxmax);
        read_file_param(fpin,zzmin,zzmax);
        read_file_param(fpin,nnx,nnz);
        read_file_param(fpin,use_sph_tmp);
        fmat2 veloc(nnz,nnx);
        std::vector<float> xcord(nnx),zcord(nnz);
        for(int iz = 0; iz < nnz; iz ++){
            zcord[iz] = zzmin + (zzmax - zzmin) / (nnz-1) * iz;
            for(int ix = 0; ix < nnx; ix ++){
                xcord[ix] = xxmin + (xxmax - xxmin) / (nnx-1) * ix;
                read_file_param(fpin,veloc(iz,ix));
            }
        }
        fpin.close();
        mesh.set_velocity(xcord.data(),zcord.data(),veloc.data(),nnx,nnz);
    }

    // create solver
    SPM2DSolver sol(mesh,mesh.veloc);

    // locate source and receivers
    sol.locate_source_stations(xs,zs,xr.data(),zr.data(),nr,mesh);

    // compute
    printf("\ncompute travel time ...\n");
    auto start = std::chrono::system_clock::now();
    sol.compute_traveltime(mesh,mesh.veloc);
    auto end = std::chrono::system_clock::now();
    auto elapse = std::chrono::duration<float>(end- start);
    printf("time elapsed = %g\n",elapse.count());

    //travel time
    for(int i = 0; i < nr; i ++){
        printf("%f %f %f\n",xr[i],zr[i],sol.ttime_recv[i]);
    }

    // compute homogeneous time
    fmat2 temp(mesh.nptstot,4);
    for(int i = 0; i < mesh.nptstot; i ++) {
        temp(i,0) = mesh.xstore[i];
        temp(i,1) = mesh.ystore[i];
        temp(i,2) = mesh.zstore[i];
        float dist;
        if(mesh.is_spherical) {
            dist = gps2dist(xs,mesh.xstore[i],zs,mesh.ystore[i],6371);
        }
        else {
            dist = std::hypot(mesh.xstore[i]-xs,mesh.ystore[i]-zs);
        }

        if(velocfile == NULL) {
            temp(i,3) = dist / 1.5;
        }
        else {
            temp(i,3) = dist / mesh.veloc.sum() * (mesh.nelmnts * NPT2);
        }
        
    }
    std::ofstream fpout("time.homo.out",std::ios::binary);
    fpout.write((char*)temp.data(),sizeof(float) * temp.size());
    fpout.close();

    // compute frechet kernel
    // fvec fdm(mesh.nptstot);
    // float lon[nx],lat[nz];
    // for(int i = 0; i < nx; i ++) {
    //     lon[i] = xmin + (xmax - xmin) / (nx - 1) * i;
    // }
    // for(int i = 0; i < nz; i ++) {
    //     lat[i] = zmin + (zmax - zmin) / (nz-1) * i;
    // }
    // fmat3 fre(nz,nx,4); fre.setZero();
    // for(int ir = 0; ir < nr; ir ++) {
    //     fdm.setZero();
    //     sol.frechet_kernel(ir,mesh,mesh.veloc,fdm);
    //     for(int inode = 0; inode < mesh.nptstot; inode ++) {
    //         float x = mesh.xstore[inode], y = mesh.ystore[inode];
    //         int iy,ix;
    //         float coef[4];
    //         bilinear(lon,lat,nx,nz,x,y,ix,iy,coef);
    //         for(int i = 0; i < 2; i++){
    //         for(int j = 0; j < 2; j++){
    //             fre(iy+i,ix+j,3) += fdm[inode] * coef[i*2+j];
    //         }}
    //     }
    // }

    // for(int iy = 0; iy < nz; iy ++) {
    // for(int ix = 0; ix < nx; ix ++) {
    //     fre(iy,ix,0) = lon[ix];
    //     fre(iy,ix,1) = lat[iy];
    //     fre(iy,ix,2) = 0.;
    // }}
    
    // fpout.open("kernel.bin",std::ios::binary);
    // fpout.write((char*)fre.data(),sizeof(float) * fre.size());
    // fpout.close();

    // save results
    sol.write_timefield(mesh,"time.out");
    sol.write_raypath(mesh,"ray.dat");

    return 0;
}
