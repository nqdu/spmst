#include "spm2d.hpp"
#include "shared/bilinear.hpp"
#include <iostream>
#include <tuple>

/**
 * @brief Get the unique points of a mesh (element by element) by returning its index
 * 
 * @param[in] x,z flattened mesh coordinates
 * @param[out] index indexes for unique points 
 * @param[in] n size(x)
 * @return  # of unique points
 */

static int
get_unique(const float* restrict x, const float* restrict z, 
           int* restrict index, int n)
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

float SPM2D:: 
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


SPM2D:: 
SPM2D(float lonmin,float latmin, float dlon, float dlat,int nlon,int nlat,bool sph)
{
    this-> initialize(lonmin,latmin,dlon,dlat,nlon,nlat,sph);
}


/**
 * @brief initialze SPM2D class
 * 
 * @param lonmin,latmin mesh lower left point 
 * @param dlon,dlat element size in lon and lat direction
 * @param nlon,nlat # of elements in lon and lat direction 
 * @param sph = true if it is spherical coordinates
 */
void SPM2D::
initialize(float lonmin,float latmin, float dlon, float dlat,int nlon,int nlat,bool sph)
{
    // allocate space
    nelemx = nlon; nelemy = nlat;
    nelmnts = nelemx * nelemy;
    ibool.resize(nelmnts,NPT2); veloc.resize(nelmnts,NPT2);
    XMIN = lonmin; YMIN = latmin; DX = dlon; DY = dlat;
    is_visited.resize(nelmnts); is_visited.setConstant(false);
    btree.resize(nelmnts);

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

    // get coordinates
    xstore.resize(nptstot); ystore.resize(nptstot); zstore.resize(nptstot);
    ttime.resize(nptstot);
    for(int ielem = 0; ielem < nelmnts; ielem ++){
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = ibool(ielem,ipt);
        xstore[inode] = xmesh(ielem,ipt);
        ystore[inode] = ymesh(ielem,ipt);
    }}
    ttime.setConstant(inf);

    // resize raypath nodes
    comming_node.resize(nptstot,2);

    // set is_spherical
    is_spherical = sph;
}

void SPM2D::
locate_source_stations(float evlo,float evla,float* restrict stlo,
                        float* restrict stla,int nr)
{
    // copy source and receiver information 
    xsource = evlo; ysource = evla;
    nreceivers = nr;
    stalocs.resize(nreceivers,3); ielem_recv.resize(nreceivers);
    ttime_recv.resize(nreceivers);
    for(int ir = 0; ir < nr; ir++){
        stalocs(ir,0) = stlo[ir];
        stalocs(ir,1) = stla[ir];
    }

    // locate source/receiver locations
    int iy,ix;
    iy = (ysource - YMIN) / DY;
    ix = (xsource - XMIN) / DX;
    if(ix < 0 || ix >= nelemx || iy < 0 || iy >= nelemy){
        printf("source is outside the domain lon=%f lat=%f\n",xsource,ysource);
        exit(1);
    }
    ielem_src = iy * nelemx + ix;
    for(int ir = 0; ir < nreceivers; ir ++){
        iy = (stalocs(ir,1) - YMIN) / DY;
        ix = (stalocs(ir,0) - XMIN) / DX;
        if(ix < 0 || ix >= nelemx || iy < 0 || iy >= nelemy){
            printf("receiver is outside the domain lon=%f lat=%f\n",stalocs(ir,0),stalocs(ir,1));
            exit(1);
        }
        ielem_recv[ir] = iy * nelemx + ix;
    }

    // resize raypath
    comming_node_recv.resize(nreceivers,2);

    // compute source and receiver 
    zsource = interp2d(topo_x.data(),topo_y.data(),topo_z.data(),
                        topo_x.size(),topo_y.size(),xsource,ysource);
    for(int ir = 0; ir < nreceivers; ir++){
        stalocs(ir,2) = interp2d(topo_x.data(),topo_y.data(),topo_z.data(),
                                topo_x.size(),topo_y.size(),stalocs(ir,0),
                                stalocs(ir,1));
    }
}