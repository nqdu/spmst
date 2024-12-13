#ifndef _SPM2D_CLASS_H
#define _SPM2D_CLASS_H

#include "numerical.hpp"

const int NPTX = 13,NPTZ = 13;
const int NPT2 = NPTX*2 + NPTZ*2 -4; // # of auxiliary nodes in x and z direction
const float earth = 6371.0;

class SPM2DMesh {
public:
    int nelemx,nelemy; // # of elements in lon and lat direction
    int nelmnts; // # of all elements
    int nptstot; // # of unique points 
    imat2 ibool; // connectivity matrix shape(nelemx*nelemz,NPT2)
    fvec xstore,ystore,zstore; // x,y and z coordinates, shape(nptstot)
    fmat2 veloc; // velocity of each node shape(nelmnts,NPT2)
    float XMIN,YMIN,DX,DY; // mesh info lon and lat min and the size of each element

    // csr graph 
    std::vector<int> xadj; // shape (nptstot)
    std::vector<int> adjncy_iel,adjncy_ipt; // (ielem,ipt) pair for adjacent vertices, shape(nadj)
    std::vector<int> adjncy_self_ipt; // ipt = adjncy_self_ipt[k] for this node, in the same element of adjncy_iel[k]

    // topography 
    fvec topo_x,topo_y;
    fmat2 topo_z;

    // spherical coordinates?
    bool is_spherical;

private:
    Eigen::Array<bool,NPT2,NPT2,1> is_connected;
public:

    SPM2DMesh() {};
    SPM2DMesh(float lonmin,float latmin, float dlon, float dlat,int nlon,int nlat,bool sph=true);
    void create_graph(bool has_discon = false);
    void initialize(float lonmin,float latmin, float dlon, float dlat,int nlon,int nlat,bool sph=true);
    void read_topography(const char *filename);
    void set_velocity(const float* lon, const float* lat,
                     const float* veloc_in,int nlon,int nlat);
    void interp_velocity(const float* lon, const float* lat,
                        const float* veloc_in,int nlon,int nlat,
                        fmat2 &veloc_out) const;

    int locate_point(float x,float y) const;
    float compute_length(float x1,float y1,float z1,float x2,float y2,float z2) const;
    
    // IO 
    void write_velocity(const char *filename) const;
    void write_vtk(const char *vtkfile) const;

    float get_velocity(float lon,float lat) const;
    float get_velocity(float lon,float lat,const fmat2 &veloc_ex) const;

private:
    void get_connected_nodes();
    void set_topography();
};

/**
 * @brief SPM2D Solver
 * 
 */
class SPM2DSolver {
public:
    fvec ttime; // travel time, shape (nptstot)
    Eigen::Array<bool,-1,1> is_fixed; // shape(nptstot)

public:    
    // source and receiver
    float xsource,ysource,zsource;
    int nreceivers;
    fmat2 stalocs; // coordinates of receivers, shape(nreceivers,3) lon/lat/r
    int ielem_src; // source element
    std::vector<int> ielem_recv; // receivers element ,shape(nreceivers)
    fvec ttime_recv; // travel time, shape (nr)

private:
    imat2 comming_node; // which node this node is updated from, shape(nptstot,2) (ielem,ipt) 
    imat2 comming_node_recv; // which node this receiver is updated from, shape(nreceivers,2)
    std::vector<float> weights; // weights for each edge, for forward simulation

    // binary tree
    int treesize;
    std::vector<std::array<float,2>> btree;
    Eigen::Array<bool,-1,1> is_visited; // shape(nptstot)

public:
    //functions
    void locate_source_stations(float evlo,float evla,const float* __restrict stlo,
                                const float* __restrict stla,int nr,
                                const SPM2DMesh &mesh);
    SPM2DSolver(const SPM2DMesh &mesh,const fmat2 &veloc);
    void initialize(const SPM2DMesh &mesh,const fmat2 &veloc);
    void write_timefield(const SPM2DMesh &mesh,const char *filename) const;
    void write_raypath(const SPM2DMesh &mesh,const char *filename) const;

    // forward simulation
    void compute_traveltime(const SPM2DMesh &mesh,const fmat2 &veloc);
    void frechet_kernel(int ir,const SPM2DMesh &mesh,const fmat2 &veloc, fvec &fdm) const;
    void recompute_time(const SPM2DMesh &mesh,const fmat2 &veloc);

private:
    void dijkstra_cpu(const SPM2DMesh &mesh);
};

#endif
