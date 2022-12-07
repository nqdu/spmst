#include "spm2d.hpp"
#include "shared/gps2dist.hpp"
#include <iostream> 

void SPM2D:: 
write_timefield(const char *filename) const
{
    FILE *fp = fopen(filename,"w");
    for(int ielem = 0; ielem < nelmnts; ielem++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = ibool(ielem,ipt);
        float x = xstore[inode], y = ystore[inode], z = zstore[inode];
        //float t = gps2dist(x,xsource,y,ysource,earth) / 3.0;
        fprintf(fp,"%f %f %f %f\n",x,y,z,ttime[inode]);
    }
    }
    fclose(fp);
}

void SPM2D:: 
write_raypath(const char *filename) const
{
    FILE *fp = fopen(filename,"w");
    for(int ir = 0; ir < nreceivers; ir++){
        fprintf(fp,"%f %f %f\n",stalocs(ir,0),stalocs(ir,1),stalocs(ir,2));
        int ielem = comming_node_recv(ir,0), ipt = comming_node_recv(ir,1);
        while(ipt !=-1){
            int inode = ibool(ielem,ipt);
            fprintf(fp,"%f %f %f\n",xstore[inode],ystore[inode],zstore[inode]);
            ielem = comming_node(inode,0); 
            ipt = comming_node(inode,1);
        }
        fprintf(fp,"%f %f %f\n",xsource,ysource,zsource);
        fprintf(fp,">\n");
    }
    fclose(fp);
}

void SPM2D::
write_velocity(const char *filename) const
{
    FILE *fp = fopen(filename,"w");
    for(int ielem = 0; ielem < nelmnts; ielem++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = ibool(ielem,ipt);
        fprintf(fp,"%f %f %f\n",xstore[inode],ystore[inode],veloc(ielem,ipt));
    }
    }
    fclose(fp);
}

void SPM2D::
set_topology(const char *filename)
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
    fvec lon(nlon),lat(nlat);
    fmat2 topo(nlat,nlon);
    for(int i = 0; i < nlon; i++) lon[i] = lonmin + dlon * i;
    for(int i = 0; i < nlat; i++) lat[i] = latmin + dlat * i;

    // read topo
    for(int iy = 0; iy < nlat; iy++){
    for(int ix = 0; ix < nlon; ix++){
        ierr = fscanf(fp,"%f",&topo(iy,ix)); assert(ierr == 1);
    }}
    topo = topo * 0.001;

    fclose(fp);

    // copy arrays to private data
    topo_x = lon; topo_y = lat;
    topo_z = topo;

    // now set topology to the SPM class
    set_topology();
}

void SPM2D::
set_topology(fvec &lon, fvec &lat, fmat2 &topo)
{
    topo_x.resize(lon.size()); topo_y.resize(lat.size());
    topo_z.resize(topo.rows(),topo.cols());
    topo_x = lon; topo_y = lat; topo_z = topo;
    set_topology();
}