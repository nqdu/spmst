#include "spm2d.hpp"
#include <iostream>
#include <chrono>

static void 
read_line(char* restrict __s, int __n, FILE* restrict __stream)
{
    __s[0] = '#';
    while (__s[0] == '#' || __s[0] == '\n'){
        assert(fgets(__s,__n,__stream) != NULL);
    }
    
}

int main(){
    int nx,nz;
    float zmin,xmin,zmax,xmax;
    float dx,dz;

    char line[256];
    FILE *fp = fopen("spmst.in","r");
    if(fp == NULL){
        printf("cannot open spmst.in\n");
        exit(1);
    }
    read_line(line,sizeof(line),fp);
    sscanf(line,"%f%f",&xmin,&zmin);
    read_line(line,sizeof(line),fp);
    sscanf(line,"%f%f",&xmax,&zmax);
    read_line(line,sizeof(line),fp);
    sscanf(line,"%d%d",&nx,&nz);
    dx = (xmax - xmin) / nx; dz = (zmax - zmin) / nz;
    fclose(fp);

    // source and receiver
    float xs,zs,dummyf; 
    int nr,ierr; 
    char dummy;
    fp = fopen("surfdata.txt","r");
    if(fp == NULL){
        printf("cannot open surfdata.txt\n");
        exit(1);
    }
    ierr = fscanf(fp,"%c%f%f%d",&dummy,&xs,&zs,&nr); assert(ierr == 4);
    float xr[nr],zr[nr];
    for(int i = 0; i < nr; i++){
        ierr = fscanf(fp,"%f%f%f",xr+i,zr+i,&dummyf);
        assert(ierr == 3);
    }
    fclose(fp);

    // set initial
    SPM2D spm2d(xmin,zmin,dx,dz,nx,nz);
    spm2d.veloc.setConstant(3.0);

    // set topography
    printf("reading topography ...\n");
    spm2d.set_topology("topo.dat");

    // locate source and receivers
    spm2d.locate_source_stations(xs,zs,xr,zr,nr);

    auto start = std::chrono::system_clock::now();
    spm2d.compute_traveltime();
    auto end = std::chrono::system_clock::now();
    auto elapse = std::chrono::duration<float>(end- start);
    printf("time elapsed = %g\n",elapse.count());

    // save results
    spm2d.write_timefield("time.out");
    spm2d.write_raypath("ray.dat");
    spm2d.write_vtk("veloc.vtk");

    return 0;
}