#include "spm2d.hpp"
#include "shared/gps2dist.hpp"
#include <iostream> 
#include "shared/bilinear.hpp"

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
    this -> set_topology();
}

void SPM2D::
set_topology(fvec &lon, fvec &lat, fmat2 &topo)
{
    topo_x.resize(lon.size()); topo_y.resize(lat.size());
    topo_z.resize(topo.rows(),topo.cols());
    topo_x = lon; topo_y = lat; topo_z = topo;
    this -> set_topology();
}

void SPM2D::
set_topology()
{
    int nx = topo_z.cols(), ny = topo_z.rows();
    for(int ielem = 0; ielem < nelmnts; ielem++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = ibool(ielem,ipt);
        float x0 = xstore[inode], y0 = ystore[inode];
        zstore[inode] = interp2d(topo_x.data(),topo_y.data(),topo_z.data(),
                                 nx,ny,x0,y0);
    }}
}

/**
 * @brief write vtk file for current velocity
 * 
 * @param vtkfile 
 */
void SPM2D:: 
write_vtk(const char *vtkfile) const
{
    FILE *fp = fopen(vtkfile,"w");
    if (fp == NULL) {
        printf("cannot open %s\n",vtkfile);
        exit(1);
    }

    // writing flags
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "material model VTK file\n");
    fprintf(fp, "ASCII\n");

    // writing grid
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp,"POINTS %d float\n",nptstot);
    for(int inode = 0; inode < nptstot; inode ++){
        fprintf(fp,"%f %f %f\n",xstore[inode],ystore[inode],zstore[inode]);
    }
    fprintf(fp,"\n");

    // writing cells
    fprintf(fp,"CELLS %d %d\n",nelmnts,nelmnts*5);
    for(int ielem = 0; ielem < nelmnts; ielem ++){
        fprintf(fp,"4 ");
        fprintf(fp,"%d %d %d %d\n",ibool(ielem,0),ibool(ielem,1),ibool(ielem,3),ibool(ielem,2));
    }
    fprintf(fp,"\n");

    // writing pixel
    fprintf(fp,"CELL_TYPES %d\n",nelmnts);
    for(int ielem = 0; ielem < nelmnts; ielem ++){
        fprintf(fp,"9\n");
    }
    fprintf(fp,"\n");

    // allocate flag 
    std::vector<int> mask_ibool(nptstot);
    memset(mask_ibool.data(),0,sizeof(int)*nptstot);
    std::vector<float> flag_val(nptstot);
    for(int ielem = 0; ielem < nelmnts; ielem ++){
    for(int ipts = 0; ipts < NPT2; ipts ++ ){
        int inode = ibool(ielem,ipts);
        if(mask_ibool[inode] == 0){
            flag_val[inode] = veloc(ielem,ipts);
            mask_ibool[inode] = 1;
        }
    }}

    // write point data
    fprintf(fp,"POINT_DATA %d\n",nptstot);
    fprintf(fp,"SCALARS Velocity float\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for(int inode = 0; inode < nptstot; inode ++){
        fprintf(fp,"%f\n",flag_val[inode]);
    }
    fprintf(fp,"\n");

    fclose(fp);
}