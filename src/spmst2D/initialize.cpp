#include "spmst2D/spmst2D.hpp"
#include "shared/bilinear.hpp"
#include "shared/IO.hpp"
#include <iostream>
#include <fstream>

/**
 * @brief read observed data for SPMST2D  
 * 
 * @param filename observed data file
 */
void SPMST2D:: 
read_obsdata(const char *filename)
{
    // open file
    printf("\nreading observed data ...\n");
    std::ifstream fp; fp.open(filename);
    if(!fp.is_open()){
        printf("cannot open %s\n",filename);
        exit(1);
    }

    // scan the file for the first time to read source-station pairs
    std::string line;
    char dummy;
    nevents = 0;
    int max_sta = 0, ndata = 0;
    while(std::getline(fp,line)) {
        nevents += 1;
        int nsta;
        float x,y;
        sscanf(line.data(),"%c%f%f%d",&dummy,&x,&y,&nsta);
        max_sta = std::max(max_sta,nsta);
        ndata += nsta;

        // skip all receivers
        for(int i = 0; i < nsta; i++){
            getline(fp,line);
        }
    }
    fp.close();

    // allocate space
    evlon.resize(nevents); evlat.resize(nevents);
    nrecvs_per_event.resize(nevents);
    stalon.resize(nevents,max_sta); stalat.resize(nevents,max_sta);
    tobs.resize(ndata);

    // read again to get all coordinates and data
    ndata = 0;
    fp.open(filename);
    for(int ievt = 0; ievt < nevents; ievt ++){
        std::getline(fp,line);
        sscanf(line.data(),"%c%f%f%d",&dummy,&evlon[ievt],&evlat[ievt],&nrecvs_per_event[ievt]);
        for(int ir = 0; ir < nrecvs_per_event[ievt]; ir++){
            std::getline(fp,line);
            float v0,dist;
            sscanf(line.data(),"%f%f%f",&stalon(ievt,ir),&stalat(ievt,ir),&v0);
            dist = compute_distance(evlon[ievt],evlat[ievt],
                                    stalon(ievt,ir),stalat(ievt,ir));
            tobs[ndata] = dist / v0;
            
            ndata += 1;
        }
    }
    fp.close();
    printf("# of obeservations = %d\n",ndata);
}

/**
 * @brief read topography file and assign it to mesh
 * 
 * @param filename 
 */
void SPMST2D:: 
read_topography(const char *filename)
{
    printf("\nreading topography from %s ...\n",filename);
    mesh.read_topography(filename);
}

/**
 * @brief read spmst2D inverse parameters file
 * 
 * @param filename 
 */
void SPMST2D:: 
read_params(const char *filename)
{
    param.read_file(filename);
}

/**
 * @brief read velocity model 
 * 
 * @param filename velocity model file
 * @param veloc_in velocity model, shape(nlat,nlon)
 * @param init_mesh create spm mesh and print information
 */
void SPMST2D:: 
read_model(const char *filename,fmat2 &veloc_in,bool init_mesh)
{
    // open model file
    std::ifstream fp; fp.open(filename);
    if(!fp.is_open()){
        printf("cannot open %s\n",filename);
        exit(1);
    }

    // read model dimension
    int nlat,nlon;
    float lonmin,lonmax,latmin,latmax;
    read_file_param(fp,nlon,nlat);
    read_file_param(fp,lonmin,lonmax);
    read_file_param(fp,latmin,latmax);
    read_file_param(fp,is_spherical);

    // create lon/lat vector
    model_lat.resize(nlat); model_lon.resize(nlon);
    for(int i = 0; i < nlon; i ++) {
        model_lon[i] = lonmin + (lonmax - lonmin) / (nlon-1.) * i;
    }
    for(int i = 0; i < nlat; i ++) {
        model_lat[i] = latmin + (latmax - latmin) / (nlat-1.) * i;
    }

    // print information
    if(init_mesh) {
        printf("\nVelocity Model Information:\n");
        if(is_spherical){
            printf("using spherical coordinates ...\n");
            printf("lonmin = %f lonmax = %f\n",lonmin,lonmax);
            printf("latmin = %f latmax = %f\n",latmin,latmax);
            printf("nlat = %d nlon = %d\n",nlat,nlon);
        }
        else{
            printf("using cartesian coordinates ...\n");
            printf("xmin = %f xmax = %f\n",lonmin,lonmax);
            printf("ymin = %f ymax = %f\n",latmin,latmax);
            printf("ny = %d nx = %d\n",nlat,nlon);
        }
    }

    // resize veloc and read it
    veloc_in.resize(nlat,nlon);
    for(int ilat = 0; ilat < nlat; ilat ++) {
    for(int ilon = 0; ilon < nlon; ilon ++) {
        read_file_param(fp,veloc_in(ilat,ilon));
    }}

    // close file
    fp.close();

    // initialize mesh if required
    if(init_mesh) {
        const int nrefine = 2;
        int nlonr = (nlon - 1) * nrefine;
        int nlatr = (nlat - 1) * nrefine;
        float dlonr = (lonmax - lonmin) / nlonr;
        float dlatr = (latmax - latmin) / nlatr;
        mesh.initialize(lonmin,latmin,dlonr,dlatr,nlonr,nlatr,is_spherical);
        mesh.create_graph(false);
    }
}

/**
 * @brief Set the velocity in SPM2DMesh
 * 
 * @param veloc_in node-based velocity model, shape(nlat,nlon)
 */
void SPMST2D:: 
set_velocity(const fmat2 &veloc_in)
{
    int nlat = veloc_in.rows(), nlon = veloc_in.cols();

    // now interp velocity on the spm2d grid
    for(int ielem = 0; ielem < mesh.nelmnts; ielem ++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = mesh.ibool(ielem,ipt);
        float x = mesh.xstore[inode], y = mesh.ystore[inode];
        mesh.veloc(ielem,ipt) = interp2d(model_lon.data(),model_lat.data(),veloc_in.data(),nlon,nlat,x,y);
    }}
}
