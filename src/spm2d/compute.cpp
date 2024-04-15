#include "spm2d.hpp"
#include "shared/bilinear.hpp"

/**
 * @brief compute strait line distance by spherical coordinates
 * 
 * @param lon1,lat1,h1: (lon,lat,height) for point 1 
 * @param lon2,lat2,h2: (lon,lat,height) for point 2 
 * @return float 
 */
static float 
cal_dist_spherical(float lon1,float lat1,float h1,
                   float lon2,float lat2,float h2)
{
    const double deg2rad = M_PI / 180.;
    double r1 = h1 + earth, r2 = h2 + earth;
    double x1 = r1 * std::cos(lat1 * deg2rad) * std::sin(lon1 * deg2rad);
    double x2 = r2 * std::cos(lat2 * deg2rad) * std::sin(lon2 * deg2rad);
    double y1 = r1 * std::cos(lat1 * deg2rad) * std::cos(lon1 * deg2rad);
    double y2 = r2 * std::cos(lat2 * deg2rad) * std::cos(lon2 * deg2rad);
    double z1 = r1 * std::sin(lat1 * deg2rad), z2 = r2 * std::sin(lat2 * deg2rad);

    double dist = (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2);
    dist = std::sqrt(dist);

    return (float)dist;
}

/**
 * @brief compute the length of one ray segment
 * 
 * @param x1,y1,z1 first point 
 * @param x2,y2,z3 the second point 
 * @return float 
 */
float SPM2DMesh::
compute_length(float x1,float y1,float z1,float x2,float y2,float z2) const
{
    float dist;
    if(is_spherical){
        dist =  cal_dist_spherical(x1,y1,z1,x2,y2,z2);
    }
    else{
        dist = (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2);
        dist = std::sqrt(dist);
    }
   
   return dist;
}

/**
 * @brief compute travel time for a given velocity model
 * 
 * @param mesh SPM2DMesh class
 * @param veloc input velcocity model, shape(nelmnts,NPT2)
 */
void SPM2DSolver::
compute_traveltime(const SPM2DMesh &mesh,const fmat2 &veloc)
{
    // zero out 
    is_fixed.setConstant(false);
    is_visited.setConstant(false);
    ttime.setConstant(inf);

    // set initial time by using source location
    treesize = 0;
    float vs;
    std::array<float,2> cordx,cordy;
    cordx = {mesh.xstore[mesh.ibool(ielem_src,0)],mesh.xstore[mesh.ibool(ielem_src,1)]};
    cordy = {mesh.ystore[mesh.ibool(ielem_src,0)],mesh.ystore[mesh.ibool(ielem_src,3)]};
    vs = interp2d(cordx.data(),cordy.data(),&veloc(ielem_src,0),2,2,xsource,ysource);
    for(int ipt = 0; ipt < NPT2; ipt ++) {
        int inode = mesh.ibool(ielem_src,ipt);
        float x = mesh.xstore[inode],y = mesh.ystore[inode], z = mesh.zstore[inode];
        float dist = mesh.compute_length(xsource,ysource,zsource,x,y,z);
        float t =  dist * 0.5 * (1. / vs + 1. / veloc(ielem_src,ipt));
        ttime[inode] = t;
        is_fixed[inode] = true;
        is_visited[inode] = true;

        // update btree
        btree[treesize] = inode;
        treesize += 1;

        // update comming node
        comming_node(inode,0) = ielem_src;
        comming_node(inode,1) = -1;
    }

    // now loop to compute SPM path
    this -> dijkstra_cpu(mesh);

    // get travel time
    for(int ir = 0; ir < nreceivers; ir ++) {
        float vr = mesh.get_velocity(stalocs(ir,0),stalocs(ir,1),veloc);
        int ielem = ielem_recv[ir];

        // initialize  
        float t = inf;
        int ipt_opt = -1;

        if(ielem == ielem_src) { // receiver is in source cell
            comming_node_recv(ir,0) = ielem_src;
            comming_node_recv(ir,1) = -1;
            float vs = mesh.get_velocity(xsource,ysource,veloc);
            float dist = mesh.compute_length(xsource,ysource,zsource,
                                            stalocs(ir,0),stalocs(ir,1),
                                            stalocs(ir,2));
            ttime_recv[ir] = dist * 0.5 * (1. / vr + 1. / vs);
            break;
        }

        // find minimum travel time
        for(int ipt = 0; ipt < NPT2; ipt ++) {
            int inode = mesh.ibool(ielem,ipt);
            float x = mesh.xstore[inode],y = mesh.ystore[inode];
            float z = mesh.zstore[inode];
            float dist = mesh.compute_length(x,y,z,stalocs(ir,0),
                                            stalocs(ir,1),stalocs(ir,2));
            float t0 = ttime[inode] + dist * 0.5 * (1. / vr + 1. / veloc(ielem,ipt));
            if (t > t0) {
                t = t0;
                ipt_opt = ipt;
            }
        }

        // get comming nodes
        ttime_recv[ir] = t;
        comming_node_recv(ir,0) = ielem;
        comming_node_recv(ir,1) = ipt_opt;
    }
}

/**
 * @brief Make a Min heap  
 * 
 * @param ttime travel timetime matrix, shape(nptstot)
 * @param ibool connectivity matrix, shape(nelmnts,NPT2)
 * @param btree vector of FMMIndex, shape
 * @param TreeSize current treesize
 */
static void 
create_minheap(const fvec &ttime,std::vector<int> &btree,int treesize)
{
    for(int i = treesize / 2 -1; i >= 0; i --){
        int childid = 2 * i + 1;
        int parentid = i;
        int inode_p =  btree[parentid];
        float timepar = ttime[inode_p];
        int inode_c  = btree[childid];
        float timechi = ttime[inode_c];

        while(childid <= treesize-1){
            if(childid < treesize - 1){
                int inode_c1 = btree[childid+1];
                float time_c1 = ttime[inode_c1];
                if(timechi > time_c1){
                    childid += 1;
                    inode_c = inode_c1;
                    timechi = time_c1;
                }
            }
            if(timepar > timechi){
                std::swap(btree[parentid],btree[childid]);
                parentid = childid;
                childid = parentid * 2 + 1;
                if(childid < treesize){
                    inode_c  = btree[childid];
                    timechi = ttime[inode_c];
                }
            }
            else{
                break;
            }
        }
    }
}

/**
 * @brief apply Dijkstra's Algorithm to compute travel time
 * 
 * @param mesh SPM2DMesh class
 */
void SPM2DSolver::
dijkstra_cpu(const SPM2DMesh &mesh)
{
    while(treesize > 0) {
        create_minheap(ttime,btree,treesize);
        int inode = btree[0];
        std::swap(btree[treesize-1],btree[0]);
        treesize -= 1;
        is_fixed[inode] = true;

        // compute adjacent nodes to update traveltime
        int istart = mesh.xadj[inode];
        int iend = mesh.xadj[inode + 1];
        for(int i = istart; i < iend; i ++) {
            int ielem = mesh.adjncy_iel[i],ipt = mesh.adjncy_ipt[i];
            int inode_a = mesh.ibool(ielem,ipt);

            // skip fixed node
            if(is_fixed[inode_a]) continue;

            // update traveltime
            float t = ttime[inode] + weights[i];
            if(is_visited[inode_a]) {
                if(t < ttime[inode_a]) {
                    ttime[inode_a] = t;
                    comming_node(inode_a,0) = ielem;
                    comming_node(inode_a,1) = mesh.adjncy_self_ipt[i];
                }
            }
            else {
                ttime[inode_a] = t;
                comming_node(inode_a,0) = ielem;
                comming_node(inode_a,1) = mesh.adjncy_self_ipt[i];
                btree[treesize] = inode_a;
                is_visited[inode_a] = true;
                treesize += 1;
            }
        }
    }
}

/**
 * @brief compute frechet kernel for current velocity model
 * 
 * @param ir receiver index
 * @param Mesh SPM2DMesh class
 * @param veloc velocity model
 * @param fdm frechet kernel, shape(nptstot)
 */
void SPM2DSolver::
frechet_kernel(int ir,const SPM2DMesh &mesh,const fmat2 &veloc, fvec &fdm) const 
{
    // initialize fdm
    fdm.setConstant(0.0);

    // some variables
    float x0,y0,z0,x1,y1,z1; // coordinates for two end points of each segment
    float v0,v1; //velocity at end of this ray path segment
    float vs; // velocity at source
    std::array<float,4> coefr{},coefs{};
    std::array<float,2> cordx,cordy;
    int ix,iy;

    // get source velocity
    const auto &xstore = mesh.xstore, &ystore = mesh.ystore;
    const auto &zstore = mesh.zstore;
    const auto &ibool = mesh.ibool;
    cordx = {xstore[ibool(ielem_src,0)],xstore[ibool(ielem_src,1)]};
    cordy = {ystore[ibool(ielem_src,0)],ystore[ibool(ielem_src,3)]};
    bilinear(cordx.data(),cordy.data(),2,2,xsource,ysource,ix,iy,coefs.data());
    vs = mesh.get_velocity(xsource,ysource,veloc);

    // get the frechet kernel for receiver ir
    // get velocity and coefr on this station
    cordx = {xstore[ibool(ielem_recv[ir],0)],xstore[ibool(ielem_recv[ir],1)]};
    cordy = {ystore[ibool(ielem_recv[ir],0)],ystore[ibool(ielem_recv[ir],3)]};
    bilinear(cordx.data(),cordy.data(),2,2,stalocs(ir,0),stalocs(ir,1),ix,iy,coefr.data());
    v0 = mesh.get_velocity(stalocs(ir,0),stalocs(ir,1),veloc); // set v0 to veloc at stations

    // now loop around all connected points to find frechet kernel
    x0 = stalocs(ir,0); y0 = stalocs(ir,1); z0 = stalocs(ir,2);
    int ielem0 = ielem_recv[ir], ipt0 = -1;
    int ielem1 = comming_node_recv(ir,0), ipt1 = comming_node_recv(ir,1);
    while(ipt1 != -1){
        int inode1 = ibool(ielem1,ipt1);
        x1 = xstore[inode1]; y1 = ystore[inode1]; z1 = zstore[inode1];
        v1 = veloc(ielem1,ipt1);
        float dist = mesh.compute_length(x0,y0,z0,x1,y1,z1);
        float term0 = -0.5 * dist / (v0 * v0);
        float term1 = -0.5 * dist / (v1 * v1);
        //printf("%f %f\n",term0,term1);
        
        // contribution from inode0/1
        fdm[inode1] += term1;
        if(ipt0 == -1){ // at station
            for(int i = 0; i< 4; i++){
                int inode0 = ibool(ielem0,i);
                fdm[inode0] += coefr[i] * term0;
            }
        }
        else{
            int inode0 = ibool(ielem0,ipt0);
            fdm[inode0] += term0;
        }

        // swap 
        v0 = v1;
        x0 = x1; y0 = y1; z0 = z1;
        ielem0 = ielem1; ipt0 = ipt1;
        ielem1 = comming_node(inode1,0), ipt1 = comming_node(inode1,1);
    }

    // contribution from source segment
    float dist = mesh.compute_length(x0,y0,z0,xsource,ysource,zsource);
    float term0 = -0.5 * dist / (v0 * v0);
    float term1 = -0.5 * dist / (vs * vs);
    int inode0 = ibool(ielem0,ipt0);
    fdm[inode0] += term0;
    for(int i = 0; i < 4 ;i++){
        int inode1 = ibool(ielem_src,i);
        fdm[inode1] += coefs[i] * term1;
    }
}

/**
 * @brief recompute travel time based on existing ray-path and new velocity
 * 
 * @param mesh Mesh class
 * @param veloc input velocity model
 */
void SPM2DSolver:: 
recompute_time(const SPM2DMesh &mesh,const fmat2 &veloc)
{
    float x0,y0,z0;
    float x1,y1,z1;
    float vs = mesh.get_velocity(xsource,ysource,veloc); // velocity at source
    int nr = ttime_recv.size();
    for(int ir = 0; ir < nr; ir ++) {
        float v0 = mesh.get_velocity(stalocs(ir,0),stalocs(ir,1),veloc); // set v0 to veloc at stations;
        float v1,t = 0.;

        // now loop around all connected points to find frechet kernel
        x0 = stalocs(ir,0); y0 = stalocs(ir,1); z0 = stalocs(ir,2);
        int ielem = comming_node_recv(ir,0), ipt = comming_node_recv(ir,1);
        while(ipt != -1){
            int inode1 = mesh.ibool(ielem,ipt);
            x1 = mesh.xstore[inode1]; y1 = mesh.ystore[inode1]; z1 = mesh.zstore[inode1];
            v1 = veloc(ielem,ipt);
            float dist = mesh.compute_length(x0,y0,z0,x1,y1,z1);
            t += dist * 0.5 * (1. / v0 + 1. / v1);

            // swap 
            v0 = v1;
            x0 = x1; y0 = y1; z0 = z1;
            ielem = comming_node(inode1,0), ipt = comming_node(inode1,1);
        }

        // contribution from source segment
        v1 = vs;
        float dist = mesh.compute_length(x0,y0,z0,xsource,ysource,zsource);
        t += dist * 0.5 * (1. / v0 + 1. / v1);

        ttime_recv[ir] = t;
    }
}