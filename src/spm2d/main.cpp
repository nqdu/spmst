#include "spm2d/spm2d.hpp"
#include "shared/IO.hpp"
#include "shared/gps2dist.hpp"
#include "shared/bilinear.hpp"
#include <iostream>
#include <chrono>

int main(int argc, char **argv){
    // check input args
    if(argc != 4){
        printf("Usage: ./syn veloc.txt surfdata.txt topo.txt\n");
        exit(1);
    }
    const char *velocfile = argv[1], *datafile=argv[2], *topofile = argv[3];

    // read header file
    int nx,nz;
    float zmin,xmin,zmax,xmax;
    float dx,dz;
    int use_sph;
    printf("\nreading velocity file %s...\n",velocfile);
    std::ifstream fpin; fpin.open(velocfile);
    if(!fpin.is_open()){
        printf("cannot open %s\n",velocfile);
        exit(1);
    }
    read_file_param(fpin,nx,nz);
    read_file_param(fpin,xmin,xmax);
    read_file_param(fpin,zmin,zmax);
    read_file_param(fpin,use_sph);
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

    // read velocity
    fmat2 veloc(nz,nx);
    std::vector<float> xcord(nx),zcord(nz);
    for(int iz = 0; iz < nz; iz ++){
        zcord[iz] = zmin + (zmax - zmin) / (nz-1) * iz;
        for(int ix = 0; ix < nx; ix ++){
            xcord[ix] = xmin + (xmax - xmin) / (nx-1) * ix;
            read_file_param(fpin,veloc(iz,ix));
        }
    }
    fpin.close();

    // initialize mesh
    SPM2DMesh mesh(xmin,zmin,dx,dz,nx*1.2,nz*1.2,use_sph==1);
    mesh.create_graph();

    // set velocity
    mesh.set_velocity(xcord.data(),zcord.data(),veloc.data(),nx,nz);

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

    // set topography
    printf("\nreading topography from %s ...\n",topofile);
    mesh.read_topography(topofile);

    // write vtk
    mesh.write_vtk("topo.vtk");
    
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

        // use average velocity
        temp(i,3) = dist / mesh.veloc.sum() * (mesh.nelmnts * NPT2);
        
    }
    std::ofstream fpout("time.homo.out",std::ios::binary);
    fpout.write((char*)temp.data(),sizeof(float) * temp.size());
    fpout.close();

    // save results
    sol.write_timefield(mesh,"time.out");
    sol.write_raypath(mesh,"ray.dat");

    return 0;
}
