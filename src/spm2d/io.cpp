#include "spm2d.hpp"
#include "shared/gps2dist.hpp"
#include "shared/IO.hpp"
#include <iostream> 

void SPM2DSolver:: 
write_timefield(const SPM2DMesh &mesh,const char *filename) const
{
    int nelmnts = mesh.nelmnts;
    FILE *fp = fopen(filename,"wb");
    fmat2 temp(nelmnts*NPT2,4);

    int c = 0;
    for(int ielem = 0; ielem < nelmnts; ielem++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = mesh.ibool(ielem,ipt);
        float x = mesh.xstore[inode], y = mesh.ystore[inode], z = mesh.zstore[inode];
        //float t = gps2dist(x,xsource,y,ysource,earth) / 3.0;
        temp(c,0) = x; temp(c,1) = y; temp(c,2) = z;
        temp(c,3) = ttime[inode];
        c += 1;
    }}

    fwrite(temp.data(),sizeof(float),temp.size(),fp);
    fclose(fp);
}

void SPM2DSolver:: 
write_raypath(const SPM2DMesh &mesh,const char *filename) const
{
    FILE *fp = fopen(filename,"w");
    for(int ir = 0; ir < nreceivers; ir++){
        fprintf(fp,"%f %f %f\n",stalocs(ir,0),stalocs(ir,1),stalocs(ir,2));
        int ielem = comming_node_recv(ir,0), ipt = comming_node_recv(ir,1);
        while(ipt !=-1){
            int inode = mesh.ibool(ielem,ipt);
            fprintf(fp,"%f %f %f\n",mesh.xstore[inode],mesh.ystore[inode],mesh.zstore[inode]);
            ielem = comming_node(inode,0); 
            ipt = comming_node(inode,1);
        }
        fprintf(fp,"%f %f %f\n",xsource,ysource,zsource);
        fprintf(fp,">\n");
    }
    fclose(fp);
}

void SPM2DMesh::
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

/**
 * @brief read topography from file
 * 
 * @param filename 
 */
void SPM2DMesh::
read_topography(const char *filename)
{
    std::ifstream fp; fp.open(filename);
    if(!fp.is_open()) {
        printf("cannot find %s\n",filename);
        exit(1);
    }

    // read meta data 
    int nlon,nlat; 
    float dlon,dlat,lonmin,latmin,lonmax,latmax;
    read_file_param(fp,nlon,nlat);
    read_file_param(fp,lonmin,lonmax);
    read_file_param(fp,latmin,latmax);
    dlon = (lonmax - lonmin) / (nlon - 1); dlat = (latmax - latmin) / (nlat - 1);
    
    // set coordinates
    fvec lon(nlon),lat(nlat);
    fmat2 topo(nlat,nlon);
    for(int i = 0; i < nlon; i++) lon[i] = lonmin + dlon * i;
    for(int i = 0; i < nlat; i++) lat[i] = latmin + dlat * i;

    // read topo
    for(int iy = 0; iy < nlat; iy++){
    for(int ix = 0; ix < nlon; ix++){
        read_file_param(fp,topo(iy,ix));
    }}
    topo = topo * 0.001;

    fp.close();

    // copy arrays to private data
    topo_x = lon; topo_y = lat;
    topo_z = topo;

    // now set topology to the SPM class
    this -> set_topography(); 
}

/**
 * @brief write vtk file for current velocity
 * 
 * @param vtkfile 
 */
void SPM2DMesh:: 
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
        fprintf(fp,"%g\n",flag_val[inode]);
    }
    fprintf(fp,"\n");

    fclose(fp);
}