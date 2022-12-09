#ifndef _SPM2D_CLASS_H
#define _SPM2D_CLASS_H

#include "numerical.hpp"

const int NPTX = 10,NPTZ = 10, NPT2 = NPTX*2 + NPTZ*2 -4; // # of auxiliary nodes in x and z direction
const float earth = 6371.0;

class SPM2D{
public:
    int nelemx,nelemy; // # of elements in lon and lat direction
    int nelmnts; // # of all elements
    int nptstot; // # of unique points 
    imat2 ibool; // connectivity matrix shape(nelemx*nelemz,NPT2)
    fvec xstore,ystore,zstore; // x,y and z coordinates, shape(nptstot)
    fvec ttime; // traveltime for all nodes, shape(nptstot)
    fmat2 veloc; // velocity of each node shape(nelmnts,NPT2)
    imat2 comming_node; // which node this node is updated from, shape(nptstot,2) (ielem,ipt) 

    // source and receivers
    float xsource,ysource,zsource; // source location lon/lat/r
    int nreceivers; // # of receivers
    imat2 comming_node_recv; // which node this receiver is updated from, shape(nreceivers,2)
    fvec ttime_recv; // travel tiem at receivers, shape(nreceivers)
    fmat2 stalocs; // coordinates of receivers, shape(nreceivers,3) lon/lat/r

private: 
    Eigen::Array<bool,-1,1> is_visited; // if this element is visited, shape(nelemx*nelemz),
    int ielem_src; // source element
    std::vector<int> ielem_recv; // receivers element ,shape(nreceivers)
    float XMIN,YMIN,DX,DY; // mesh info lon and lat min and the size of each element
    std::vector<int> btree; // binary tree to do heap sort
    int treesize;// size of btree

    // topography 
    fvec topo_x,topo_y;
    fmat2 topo_z;
    

public:

    SPM2D() {};
    SPM2D(float lonmin,float latmin, float dlon, float dlat,int nlon,int nlat);
    void initialize(float lonmin,float latmin, float dlon, float dlat,int nlon,int nlat);
    void set_topology(fvec &lon, fvec &lat, fmat2 &topo);
    void set_topology(const char *filename);
    void locate_source_stations(float evlo,float evla,float* restrict stlo,
                                float* restrict stla,int nr);

    float compute_distance(float x1,float y1,float z1,float x2,float y2,float z2) const;
    void compute_traveltime();
    void frechet_kernel(int ir, fvec &fdm) const;

    // IO 
    void write_timefield(const char *filename) const;
    void write_raypath(const char *filename) const;
    void write_velocity(const char *filename) const;

private:
    float get_velocity(float x,float y) const;
    void set_topology();
};

#endif