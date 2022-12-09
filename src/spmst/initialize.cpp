#include "spmst.hpp"
#include "shared/bilinear.hpp"
#include "shared/gps2dist.hpp"
#include <iostream>

/**
 * @brief read observed data for SPMST  
 * 
 * @param filename observed data file
 */
void SPMST:: 
read_obsdata(const char *filename)
{
    printf("\nreading observed data ...\n");
    FILE *fp; 
    if((fp = fopen(filename,"r")) == NULL){
        printf("cannot open %s\n",filename);
        exit(1);
    }

    // scan the file for the first time to read source-station pairs
    char line[256];
    char dummy;
    nevents = 0;
    int max_sta = 0, ndata = 0;
    while(!feof(fp)){
        if(fgets(line,sizeof(line),fp) == NULL) break;
        nevents += 1;
        int nsta;
        float x,y;
        sscanf(line,"%c%f%f%d",&dummy,&x,&y,&nsta);
        max_sta = std::max(max_sta,nsta);
        ndata += nsta;

        // skip all receivers
        for(int i = 0; i < nsta; i++){
            assert(fgets(line,sizeof(line),fp) != NULL);
        }
    }
    fclose(fp);

    // allocate space
    evlon.resize(nevents); evlat.resize(nevents);
    nrecvs_per_event.resize(nevents);
    stalon.resize(nevents,max_sta); stalat.resize(nevents,max_sta);
    tobs.resize(ndata);

    // read again to get all coordinates and data
    ndata = 0;
    fp = fopen(filename,"r");
    for(int ievt = 0; ievt < nevents; ievt ++){
        assert(fgets(line,sizeof(line),fp) != NULL);
        sscanf(line,"%c%f%f%d",&dummy,&evlon[ievt],&evlat[ievt],&nrecvs_per_event[ievt]);
        for(int ir = 0; ir < nrecvs_per_event[ievt]; ir++){
            assert(fgets(line,sizeof(line),fp) != NULL);
            float v0;
            sscanf(line,"%f%f%f",&stalon(ievt,ir),&stalat(ievt,ir),&v0);
            tobs[ndata] = gps2dist(evlon[ievt],stalon(ievt,ir),evlat[ievt],stalat(ievt,ir),earth) / v0;
            ndata += 1;
        }
    }
    fclose(fp);
}

/**
 * @brief read topography file 
 * 
 * @param filename topography file
 * @param lon/lat longitude and latitude for the topo
 * @param topo 2D topography, shape(nlat,nlon)
 */
static void 
read_topo_file(const char *filename,fvec &lon,fvec &lat,fmat2 &topo)
{
    FILE *fp = fopen(filename,"r");
    if(fp == NULL) {
        printf("cannot find %s\n",filename);
        exit(1);
    }

    // read meta data 
    int ierr;
    int nlon,nlat; 
    float dlon,dlat,lonmin,latmin,lonmax,latmax;
    ierr = fscanf(fp,"%d%d",&nlon,&nlat); assert(ierr == 2);
    ierr = fscanf(fp,"%f%f%f%f",&lonmin,&lonmax,&latmin,&latmax); assert(ierr == 4);
    dlon = (lonmax - lonmin) / (nlon - 1); dlat = (latmax - latmin) / (nlat - 1);
    
    // set coordinates
    lon.resize(nlon),lat.resize(nlat);
    topo.resize(nlat,nlon);
    for(int i = 0; i < nlon; i++) lon[i] = lonmin + dlon * i;
    for(int i = 0; i < nlat; i++) lat[i] = latmin + dlat * i;

    // read topo
    for(int iy = 0; iy < nlat; iy++){
    for(int ix = 0; ix < nlon; ix++){
        ierr = fscanf(fp,"%f",&topo(iy,ix)); assert(ierr == 1);
    }}
    topo = topo * 0.001;

    fclose(fp);
}

void SPMST:: 
read_topography(const char *filename)
{
    printf("\nreading topography ...\n");
    fvec lon,lat; 
    fmat2 topo;
    read_topo_file(filename,lon,lat,topo);
    spm2dbase.set_topology(lon,lat,topo);
}

static void 
read_line(char* restrict __s, int __n, FILE* restrict __stream)
{
    __s[0] = '#';
    while (__s[0] == '#' || __s[0] == '\n'){
        assert(fgets(__s,__n,__stream) != NULL);
    }
    
}

/**
 * @brief read spmst parameter file
 * 
 * @param filename 
 */
void SPMST:: 
read_spmst_params(const char *filename)
{
    FILE *fp;
    if((fp = fopen(filename,"r")) == NULL){
        printf("cannot open %s\n",filename);
        exit(1);
    }
    char line[256];

    read_line(line,sizeof(line),fp);
    sscanf(line,"%f%f",&lonmin,&latmin);
    read_line(line,sizeof(line),fp);
    sscanf(line,"%f%f",&lonmax,&latmax);
    read_line(line,sizeof(line),fp);
    sscanf(line,"%d%d",&nlon,&nlat);

    // print information
    printf("\nMesh Information:\n");
    printf("lonmin = %f lonmax = %f\n",lonmin,lonmax);
    printf("latmin = %f latmax = %f\n",latmin,latmax);
    printf("nlat = %d nlon = %d\n",nlat,nlon);

    // read lsqr params
    read_line(line,sizeof(line),fp); sscanf(line,"%f%f",&smooth,&damp);
    read_line(line,sizeof(line),fp); sscanf(line,"%f%f",&vmin,&vmax);
    read_line(line,sizeof(line),fp); sscanf(line,"%d",&niters);

    // if synthetic data
    int is_syn;
    read_line(line,sizeof(line),fp); sscanf(line,"%d",&is_syn);
    do_synthetic = is_syn == 1;
    velinit.resize(nlat,nlon);
    if(do_synthetic) veltrue.resize(nlat,nlon);

    fclose(fp);

    // allocate space for solver (with more refined grid)
    float dlon = (lonmax - lonmin) / (nlon * 4 - 4);
    float dlat = (latmax - latmin) / (nlat * 4 - 4); 
    spm2dbase.initialize(lonmin,latmin,dlon,dlat,nlon*4 - 4,nlat*4 - 4);
}

void SPMST:: 
read_velocity(const char *filename,fmat2 &veloc_in)
{
    FILE *fp;
    if((fp = fopen(filename,"r")) == NULL){
        printf("cannot open %s\n",filename);
        exit(1);
    }

    for(int ilat = 0; ilat < nlat; ilat ++){
    for(int ilon = 0; ilon < nlon; ilon ++){
        assert(fscanf(fp,"%f",&veloc_in(ilat,ilon)) == 1);
    }}
    fclose(fp);
}

void SPMST:: 
set_velocity(const fmat2 &veloc_in,SPM2D &spm2d)
{
    float lon[nlon],lat[nlat];
    for(int i = 0; i < nlon; i++){
        lon[i] = lonmin + (lonmax - lonmin) / (nlon - 1) * i;
    }
    for(int i = 0; i < nlat; i++){
        lat[i] = latmin + (latmax - latmin) / (nlat - 1) * i;
    }

    // now interp velocity on the spm2d grid
    for(int ielem = 0; ielem < spm2d.nelmnts; ielem ++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = spm2d.ibool(ielem,ipt);
        float x = spm2d.xstore[inode], y = spm2d.ystore[inode];
        spm2d.veloc(ielem,ipt) = interp2d(lon,lat,veloc_in.data(),nlon,nlat,x,y);
    }}
}