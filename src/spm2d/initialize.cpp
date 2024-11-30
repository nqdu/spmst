#include "spm2d.hpp"
#include "shared/bilinear.hpp"
#include <iostream>
#include <tuple>
#include <unordered_map>

/**
 * @brief Get the unique points of a mesh (element by element) by returning its index
 * 
 * @param[in] x,z flattened mesh coordinates
 * @param[out] index indexes for unique points 
 * @param[in] n size(x)
 * @return  # of unique points
 */
static int
get_unique(const float* __restrict x, const float* __restrict z, 
           int* __restrict index, int n)
{
    // get min and max value
    float xmin = x[0],xmax = x[0];
    float zmin = z[0],zmax = z[0];
    for(int i = 1 ; i < n; i++){
        xmin = std::min(xmin,x[i]); xmax = std::max(xmax,x[i]);
        zmin = std::min(zmin,z[i]); zmax = std::max(zmax,z[i]);
    }

    // tolerance
    double xtol = 1.0e-4 * (xmax - xmin);
    double ztol = 1.0e-4 * (zmax - zmin);

    int nc = 0;
    for(int i = 0; i < n; i ++)index[i] = -1;

    typedef std::tuple<double,double,int> fpair;
    std::vector<fpair> cords(n);
    for(int i = 0 ; i < n; i++){
        cords[i] = std::make_tuple(x[i],z[i],i);
    }
    using std::get;
    std::sort(cords.begin(), cords.end(),[xtol,ztol](fpair & r1, fpair& r2){
        if(std::abs(get<0>(r1) -get<0>(r2)) <= xtol){
            return get<1>(r1) < get<1>(r2) - ztol;
        }
        else{
            return get<0>(r1) < get<0>(r2) - xtol;
        }
    });

    int prev = 0;
    index[get<2>(cords[0])] = nc;
    for(int i = 1; i < n; i++){
        if ( (std::abs(get<0>(cords[i]) - get<0>(cords[prev])) > xtol) ||
            (std::abs(get<1>(cords[i]) - get<1>(cords[prev])) > ztol) ){
            prev = i;
            nc += 1; 
        }

        index[get<2>(cords[i])] = nc;
    }

    return nc + 1;
}

/**
 * @brief Get the connectivity matrix in one cell
 * 
 */
void SPM2DMesh::
get_connected_nodes()
{
    // get coordinates for 4 anchor points
    fmat2 index(NPT2,2);
    int ipt = 0;
    for(int iy = 0; iy < 2; iy++){
    for(int ix = 0; ix < 2; ix++){
        index(ipt,0) = (NPTX - 1) * ix;
        index(ipt,1) = (NPTZ - 1) * iy; 
        ipt += 1;
    }}
    for(int ix = 0; ix < NPTX-2; ix++){
        index(ipt,0) = ix + 1.;
        index(ipt + NPTX-2,0) = index(ipt,0);
        index(ipt,1) = 0;
        index(ipt+NPTX-2,1) = NPTZ - 1;
        ipt += 1;
    }
    ipt += NPTX - 2;
    for(int iy = 0; iy < NPTZ - 2; iy ++){
        index(ipt,0) = 0;
        index(ipt + NPTZ-2,0) = NPTZ-1;
        index(ipt,1) = iy + 1;
        index(ipt+NPTZ-2,1) = iy + 1;
        ipt += 1;
    }

    Eigen::Matrix<float,4,2,1> side_vec; 
    for(int i = 0; i < 2; i ++){
        side_vec(i,0) = index(0,0) - index(i+1,0);
        side_vec(i,1) = index(0,1) - index(i+1,1);
        side_vec(i+2,0) = index(3,0) - index(2-i,0);
        side_vec(i+2,1) = index(3,1) - index(2-i,1);
    }

    // now loop around all nodes
    const float tol = 1.0e-5;
    is_connected.setConstant(true);
    for(int ipt0 = 0; ipt0 < NPT2; ipt0 ++){
        Eigen::Vector2f p0,vec;
        p0 = index.row(ipt0);
        for(int ipt1 = 0; ipt1 < NPT2; ipt1 ++){
            if(ipt0 == ipt1){
                is_connected(ipt0,ipt1) = false;
                continue;
            }

            // vec
            Eigen::Vector2f p1;
            p1 = index.row(ipt1);
            vec = p0 - p1;

            // now check if it's parallel to side vector
            for(int i = 0; i < 4; i ++){
                float s = vec.dot(side_vec.row(i));
                if(std::abs(s) < tol){ // parallel
                    float l = vec.norm();
                    bool flagx = p0[0] == p1[0];
                    bool flagy = p0[1] == p1[1];
                    bool flag0 = flagx && !flagy && (p0[0] == 0 || p0[0] == NPTX-1);
                    bool flag1 = !flagx && flagy && (p0[1] == 0 || p0[1] == NPTZ-1);
                    bool flag = std::abs(l - 1) > tol;
                    if((flag0 && flag) ||(flag && flag1)){
                        is_connected(ipt0,ipt1) = false;
                    }
                }
            }
        }
    }
}

/**
 * @brief Get the optimal global numbering to reduce cache miss
 * 
 * @param ibool connectivity matrix, shape(nelmnts,NPT2)
 */
static void 
get_indirect_connectivity(imat2 &ibool)
{
    int nelmnts = ibool.rows();
    int nglob = ibool.maxCoeff() + 1;
    Eigen::Array<int,-1,1> mask_ibool; mask_ibool.resize(nglob); mask_ibool.setConstant(-1);
    imat2 ibool_copy = ibool;

    // access each node to reduce cache miss
    int inum = 0;
    for(int ielem = 0; ielem < nelmnts; ielem ++){
    for(int ipt = 0; ipt < NPT2; ipt ++){
        if (mask_ibool[ibool_copy(ielem,ipt)] == -1){
            ibool(ielem,ipt) = inum;
            mask_ibool[ibool_copy(ielem,ipt)] = inum;
            inum += 1;
        }
        else{
            ibool(ielem,ipt) = mask_ibool[ibool_copy(ielem,ipt)];
        }
    }}
}

/**
 * @brief Get the velocity for a given point by bilinear interpolation
 * 
 * @param lon,lat  x and y coordinates 
 * @return float 
 */
float SPM2DMesh:: 
get_velocity(float lon,float lat) const
{
    int iy = ( lat - YMIN) / DY, ix = (lon - XMIN) / DX;
    int ielem = iy * nelemx + ix;
    float x1 = XMIN + DX * ix, x2  = x1 + DX;
    float y1 = YMIN + DY * iy, y2 = y1 + DY;
    float f1 = (x2 - lon)/DX * veloc(ielem,0) + (lon - x1) / DX * veloc(ielem,1);
    float f2 = (x2 - lon)/DX * veloc(ielem,2) + (lon - x1) / DX * veloc(ielem,3);
    float out = (y2 - lat) / DY * f1 + (lat - y1) / DY * f2;

    return out;
}

/**
 * @brief Get the velocity for a given point by bilinear interpolation of external velocity
 * 
 * @param lon,lat  x and y coordinates 
 * @param veloc_ex external velocity,shape(nelmtns,NPT2)
 * @return float 
 */
float SPM2DMesh:: 
get_velocity(float lon,float lat,const fmat2 &veloc_ex) const
{
    int iy = ( lat - YMIN) / DY, ix = (lon - XMIN) / DX;
    int ielem = iy * nelemx + ix;
    float x1 = XMIN + DX * ix, x2  = x1 + DX;
    float y1 = YMIN + DY * iy, y2 = y1 + DY;
    float f1 = (x2 - lon)/DX * veloc_ex(ielem,0) + (lon - x1) / DX * veloc_ex(ielem,1);
    float f2 = (x2 - lon)/DX * veloc_ex(ielem,2) + (lon - x1) / DX * veloc_ex(ielem,3);
    float out = (y2 - lat) / DY * f1 + (lat - y1) / DY * f2;

    return out;
}

SPM2DMesh:: 
SPM2DMesh(float lonmin,float latmin, float dlon, float dlat,int nlon,int nlat,bool sph)
{
    this-> initialize(lonmin,latmin,dlon,dlat,nlon,nlat,sph);
}


/**
 * @brief initialze SPM2DMesh class
 * 
 * @param lonmin,latmin mesh lower left point 
 * @param dlon,dlat element size in lon and lat direction
 * @param nlon,nlat # of elements in lon and lat direction 
 * @param sph = true if it is spherical coordinates
 */
void SPM2DMesh::
initialize(float lonmin,float latmin, float dlon, float dlat,int nlon,int nlat,bool sph)
{
    // allocate space
    nelemx = nlon; nelemy = nlat;
    nelmnts = nelemx * nelemy;
    ibool.resize(nelmnts,NPT2); veloc.resize(nelmnts,NPT2);
    XMIN = lonmin; YMIN = latmin; DX = dlon; DY = dlat;

    // allocate space for lon/lat and compute coordinates
    fmat2 xmesh(nelmnts,NPT2), ymesh(nelmnts,NPT2);
    for(int ielemy = 0; ielemy < nelemy; ielemy++){
    for(int ielemx = 0; ielemx < nelemx; ielemx++){
        // get lower left coordinates 
        float x = XMIN + DX * ielemx;
        float y = YMIN + DY * ielemy;
        int ielem = ielemy * nelemx + ielemx;

        // get coordinates for 4 anchor points
        int ipt = 0;
        for(int iy = 0; iy < 2; iy++){
        for(int ix = 0; ix < 2; ix++){
            xmesh(ielem,ipt) = x + DX * ix;
            ymesh(ielem,ipt) = y + DY * iy; 
            ipt += 1;
        }}
        for(int ix = 0; ix < NPTX-2; ix++){
            xmesh(ielem,ipt) = x + DX / (NPTX-1) * (ix + 1);
            xmesh(ielem,ipt + NPTX-2) = xmesh(ielem,ipt);
            ymesh(ielem,ipt) = y;
            ymesh(ielem,ipt + NPTX -2) = y + DY;
            ipt += 1;
        }
        ipt += NPTX - 2;
        for(int iy = 0; iy < NPTZ - 2; iy ++){
            xmesh(ielem,ipt) = x;
            xmesh(ielem,ipt + NPTZ-2) = x + DX;
            ymesh(ielem,ipt) = y + DY / (NPTZ - 1) * (iy + 1);
            ymesh(ielem,ipt + NPTZ -2) = ymesh(ielem,ipt);
            ipt += 1;
        }
    }}
    int npts = (NPTX - 2) * nelemx * (nelemy + 1) + 
                (NPTZ-2) * nelemy * (nelemx +1) + 
                (nelemx + 1) * (nelemy + 1);

    // get connectivity matrix and reduce cache missing
    nptstot = get_unique(xmesh.data(),ymesh.data(),ibool.data(),ibool.size());
    if(npts != nptstot){
        printf("The threshold in get_unique is not adequate\n");
        printf("%d %d\n",npts,nptstot);
        exit(1);
    }
    get_indirect_connectivity(ibool);

    // find connectivity matrix in one cell
    this -> get_connected_nodes();

    // get coordinates
    xstore.resize(nptstot); ystore.resize(nptstot); zstore.resize(nptstot);
    for(int ielem = 0; ielem < nelmnts; ielem ++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = ibool(ielem,ipt);
        xstore[inode] = xmesh(ielem,ipt);
        ystore[inode] = ymesh(ielem,ipt);
    }}

    // set is_spherical
    is_spherical = sph;
}

int SPM2DMesh::
locate_point(float x,float y) const
{
    // locate source/receiver locations
    int iy,ix;
    iy = (y - YMIN) / DY;
    ix = (x - XMIN) / DX;
    if(ix < 0 || ix >= nelemx || iy < 0 || iy >= nelemy){
        printf("point is outside the domain lon=%f lat=%f\n",x,y);
        exit(1);
    }

    return iy * nelemx + ix;
}

void SPM2DMesh::
set_topography()
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
 * @brief interp mesh-velocity based on continuous input velocity
 * 
 * @param lon x coordinates for topo, shape (nlon) 
 * @param lat y coordinates for topo, shape (nlat) 
 * @param veloc_in velocity, shape(nlat,nlon)
 * @param veloc_out velocity on mesh, shape(nelmnts,NPT2)
 */
void SPM2DMesh:: 
interp_velocity(const float* lon, const float* lat,
                const float* veloc_in,int nlon,int nlat,
                fmat2 &veloc_out) const
{
    for(int ielem = 0; ielem < nelmnts; ielem++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = ibool(ielem,ipt);
        float x0 = xstore[inode], y0 = ystore[inode];
        veloc_out(ielem,ipt) = interp2d(lon,lat,veloc_in,
                                    nlon,nlat,x0,y0);
    }}
}

/**
 * @brief Set velocity for SPM mesh
 * 
 * @param lon x coordinates for topo, shape (nlon) 
 * @param lat y coordinates for topo, shape (nlat) 
 * @param veloc_in velocity, shape(nlat,nlon)
 */
void SPM2DMesh::
set_velocity(const float* lon, const float* lat,
             const float* veloc_in,int nlon,int nlat)
{
    this -> interp_velocity(lon,lat,veloc_in,nlon,nlat,veloc);
}

/**
 * @brief create dual graph based on is_connected nodes
 * @param has_discon if the velocity model contains discontinuities
 * 
 */
void SPM2DMesh::create_graph(bool has_discon)
{
    // create adjacent vertices
    if(has_discon) {
        printf("\ncreating graph, with discontinuities ...\n");
    }
    else {
        printf("\ncreating graph, continuous model ...\n");
    }
    
    std::vector<std::vector<int64_t>> adjlist;
    adjlist.resize(nptstot);

    // loop to add adjacent nodes
    for(int ielem = 0; ielem < nelmnts; ielem ++) {
    for(int ipt = 0; ipt < NPT2; ipt ++) {
        int inode = ibool(ielem,ipt);
        for(int ipt1 = 0; ipt1 < NPT2; ipt1 ++) {
            if(is_connected(ipt,ipt1)) {
                int64_t idx1 = ielem * NPT2 + ipt1;
                adjlist[inode].push_back(idx1);
            }
        }
    }}

    // handle continuous model
    if(!has_discon) {
        for(int inode = 0; inode < nptstot; inode ++) {
            auto &list = adjlist[inode];
            int n = list.size();
            std::unordered_map<int,int64_t> vmap;
            for(int i = 0; i < n; i ++) {
                int64_t idx = list[i];
                int ielem = idx / NPT2, ipt = idx % NPT2;
                vmap[ibool(ielem,ipt)] = idx;
            }

            // get unique v instead
            int n_uniq = vmap.size();
            list.resize(0); list.reserve(n_uniq);
            for(auto it : vmap) {
                list.push_back(it.second);
            }
        }
    }

    // count edges
    int64_t n = 0;
    for(int64_t ipt = 0; ipt < nptstot; ipt ++) {
        n += adjlist[ipt].size();
    }

    // allcoate csr graph
    xadj.resize(nptstot + 1);
    adjncy_iel.resize(n); adjncy_ipt.resize(n);
    adjncy_self_ipt.resize(n);
    printf("# of nodes/edges = %d %ld\n",nptstot,(long)(n));

    // set value
    xadj[0] = 0;
    for(int inode = 0; inode < nptstot; inode ++) {
        int istart = xadj[inode];
        int iend = istart + adjlist[inode].size();
        xadj[inode + 1] = iend;

        for(int i = istart; i < iend; i ++) {
            int64_t idx = adjlist[inode][i-istart];
            int ielem = idx / NPT2, ipt = idx % NPT2;
            adjncy_iel[i] = ielem;
            adjncy_ipt[i] = ipt;
            int ipt0 = 0;
            for(;ipt0 < NPT2; ipt0 ++) {
                if(ibool(ielem,ipt0) == inode) {
                    break;
                }
            }
            adjncy_self_ipt[i] = ipt0;
        }
    }
}

/**
 * @brief locate source and receivers in SPM mesh
 * 
 * @param evlo,evla  x/y coordinates of source 
 * @param stlo,stlo  x/y coordinates of receivers, shape(nr) 
 * @param nr  # of stations
 */
void SPM2DSolver::
locate_source_stations(float evlo,float evla,const float* __restrict stlo,
                        const float* __restrict stla,int nr,
                        const SPM2DMesh &mesh)
{
    // copy source and receiver information 
    xsource = evlo; ysource = evla;
    nreceivers = nr;
    stalocs.resize(nreceivers,3); ielem_recv.resize(nreceivers);
    ttime_recv.resize(nreceivers);
    comming_node_recv.resize(nreceivers,2);
    for(int ir = 0; ir < nr; ir++){
        stalocs(ir,0) = stlo[ir];
        stalocs(ir,1) = stla[ir];
    }

    // locate source 
    ielem_src = mesh.locate_point(xsource,ysource);
    zsource = interp2d(mesh.topo_x.data(),mesh.topo_y.data(),
                       mesh.topo_z.data(),mesh.topo_x.size(),
                       mesh.topo_y.size(),xsource,ysource);

    // locate receiver elements
    for(int ir = 0; ir < nreceivers; ir ++){
        ielem_recv[ir] = mesh.locate_point(stalocs(ir,0),stalocs(ir,1));
        stalocs(ir,2) = interp2d(mesh.topo_x.data(),mesh.topo_y.data(),
                        mesh.topo_z.data(),mesh.topo_x.size(),
                        mesh.topo_y.size(),stalocs(ir,0),stalocs(ir,1));
    }
}

SPM2DSolver::
SPM2DSolver(const SPM2DMesh &mesh,const fmat2 &veloc)
{
    this -> initialize(mesh,veloc);
}

/**
 * @brief initialize solver, allocating arrays 
 * 
 * @param mesh SPM2DMesh class
 * @param veloc velocity model, shape(nelmnts,NPT2)
 */
void SPM2DSolver::
initialize(const SPM2DMesh &mesh,const fmat2 &veloc)
{
    int nptstot = mesh.nptstot;
    ttime.resize(nptstot);
    comming_node.resize(nptstot,2);
    is_fixed.resize(nptstot);
    is_visited.resize(nptstot);
    btree.resize(nptstot);

    // compute weights on the graph
    weights.resize(mesh.adjncy_iel.size());
    for(int i = 0; i < nptstot; i ++) {
        int istart = mesh.xadj[i];
        int iend = mesh.xadj[i + 1];

        for(int j = istart; j < iend; j ++) {
            int ielem = mesh.adjncy_iel[j], ipt = mesh.adjncy_ipt[j];
            int inode = mesh.ibool(ielem,ipt);

            // compute distance
            float dist = mesh.compute_length(
                mesh.xstore[i],mesh.ystore[i],mesh.zstore[i],
                mesh.xstore[inode],mesh.ystore[inode],mesh.zstore[inode]
            );

            // check which node it concide in mesh (ielem)
            int ipt0 = mesh.adjncy_self_ipt[j];

            // get velocity and weights
            float v0 = veloc(ielem,ipt0);
            float v1 = veloc(ielem,ipt);
            weights[j] = dist * 0.5 * (1. / v1 + 1. / v0);
        }
    }
}