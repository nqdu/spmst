#include "spm2d.hpp"
#include "shared/bilinear.hpp"

/**
 * @brief find minimum traveltime in this cell with index ielem
 * 
 * @param ielem element index
 * @param ibool connectivity matrix, shape(nelmnts,NPT2)
 * @param ttime travel time on each node, shape(nptstot)
 * @return float 
 */
static float 
find_min_time(int ielem,const imat2 &ibool,const fvec &ttime)
{
    float t0 = inf;
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = ibool(ielem,ipt);
        t0 = std::min(t0,ttime[inode]);
    }

    return t0;
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
create_minheap(const fvec &ttime,const imat2 &ibool,std::vector<int> &btree,int treesize)
{
   
    for(int i=1;i<treesize;i++){ // loop over all child nodes
        int childid = i, parentid = (i-1) / 2;
        int ielem_child  = btree[childid];
        int ielem_parent =  btree[parentid];
        float timepar = find_min_time(ielem_parent,ibool,ttime);
        float timechi = find_min_time(ielem_child,ibool,ttime);
        
        while(parentid >=0 && timechi < timepar){
            std::swap(btree[parentid],btree[childid]);

            // update child and parent id
            childid = parentid;
            parentid = (childid - 1) / 2;
            ielem_parent =  btree[parentid];
            timepar = find_min_time(ielem_parent,ibool,ttime);
            // childid = parentid;
            // parentid = (childid - 1) / 2;
            // ielem_child  = btree[childid];
            // ielem_parent =  btree[parentid];
            // timepar = find_min_time(ielem_parent,ibool,ttime);
            // timechi = find_min_time(ielem_child,ibool,ttime);
        }
    }
}

/**
 * @brief compute travel time for current source-stations settings
 * 
 */
void SPM2D::
compute_traveltime()
{
    // initialize 
    is_visited.setConstant(false);
    ttime.setConstant(inf);

    // first compute travel time in source cell and add it to binary tree
    float vs = get_velocity(xsource,ysource);
    for(int ipt = 0; ipt < NPT2; ipt++){
        int inode = ibool(ielem_src,ipt);
        float dist = compute_length(xstore[inode],ystore[inode],zstore[inode],
                                     xsource,ysource,zsource);
        ttime[inode] = dist * 2.0 / (vs + veloc(ielem_src,ipt));
        
        // the minimum travel time is comming from source 
        comming_node(inode,0) = ielem_src;
        comming_node(inode,1) = -1;
    }
    
    // compute receiver travel time if it's in the source cell
    for(int ir = 0; ir < nreceivers; ir++){
        if(ielem_src != ielem_recv[ir]) continue;
        float vr = get_velocity(stalocs(ir,0),stalocs(ir,1));
        float dist = compute_length(xsource,ysource,zsource,stalocs(ir,0),
                                      stalocs(ir,1),stalocs(ir,2));
        ttime_recv[ir] = dist * 2.0 / (vs + vr);
        comming_node_recv(ir,0) = ielem_src;
        comming_node_recv(ir,1) = -1;
    }
    treesize = 1;
    btree[treesize-1] = ielem_src;
    is_visited[ielem_src] = true;

    while(treesize > 0){
        // heap sort to get cell with minimum travel time
        create_minheap(ttime,ibool,btree,treesize);
        std::swap(btree[treesize-1],btree[0]);
        int ielem = btree[treesize-1];
        //printf("%d %f\n",ielem,find_min_time(ielem,ibool,ttime));
        treesize -= 1; // remove this point

        // add neighbor cells into btree
        int iz0 = ielem / nelemx, ix0 = ielem % nelemx;
        for(int i = -1; i<=1; i++){
        for(int j = -1; j<=1; j++){
            int iz = iz0 + i, ix = ix0 + j;
            if(iz >=0 && iz < nelemy && ix >=0 && ix < nelemx){
                int im = iz * nelemx + ix;
                if(!is_visited[im]){
                    btree[treesize] = im;
                    treesize += 1;
                    is_visited[im] = true;
                 }
            }
        }}

        // skip source cell because it has already done
        if(ielem == ielem_src) continue;

        // get all touched nodes in this cell
        std::vector<int> ifront; ifront.reserve(NPT2);
        for(int ipt = 0; ipt < NPT2; ipt++){
            int inode = ibool(ielem,ipt);
            if(ttime[inode] != inf) {
                ifront.push_back(ipt);
            }
        }

        // loop over all nodes and compute minimum travel time from all possible nodes
        for(int ipt = 0; ipt < NPT2; ipt++){ 
            int inode = ibool(ielem,ipt);
            if(ttime[inode] != inf) continue;
            float t0 = inf, t{};
            int ipt0 = -1;
            for(int ipt1 :ifront){
                int inode1 = ibool(ielem,ipt1);
                float dist = compute_length(xstore[inode],ystore[inode],zstore[inode],
                                              xstore[inode1],ystore[inode1],zstore[inode1]);
                t = ttime[inode1] + dist * 2.0 / (veloc(ielem,ipt1) + veloc(ielem,ipt));
                if(t < t0){
                    t0 = t;
                    ipt0 = ipt1;
                }
            }

            // now cache minimum traveltime and which node it is comming from
            ttime[inode] = t0;
            comming_node(inode,0) = ielem;
            comming_node(inode,1) = ipt0;
        }

        // compute receiver traveltime if required
        for(int ir = 0; ir < nreceivers; ir++){
            if(ielem != ielem_recv[ir]) continue;
            float v0 = get_velocity(stalocs(ir,0),stalocs(ir,1));
            float t0 = inf;
            int ipt0 = -1;;
            for(int ipt1 :ifront){
                int inode1 = ibool(ielem,ipt1);
                float dist = compute_length(stalocs(ir,0),stalocs(ir,1),stalocs(ir,2),
                                            xstore[inode1],ystore[inode1],zstore[inode1]);
                float t = ttime[inode1] + dist * 2.0 / (v0 + veloc(ielem,ipt1));
                if(t < t0){
                    t0 = t;
                    ipt0 = ipt1;
                }
            }
            ttime_recv[ir] = t0;
            comming_node_recv(ir,0) = ielem;
            comming_node_recv(ir,1) = ipt0;
        }
    }
}

/**
 * @brief compute frechet kernel for current velocity model
 * 
 * @param ir receiver index
 * @param fdm frechet kernel, shape(nptstot)
 */
void SPM2D::
frechet_kernel(int ir,fvec &fdm) const 
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
    cordx = {xstore[ibool(ielem_src,0)],xstore[ibool(ielem_src,1)]};
    cordy = {ystore[ibool(ielem_src,0)],ystore[ibool(ielem_src,3)]};
    bilinear(cordx.data(),cordy.data(),2,2,xsource,ysource,ix,iy,coefs.data());
    vs = get_velocity(xsource,ysource);

    // get the frechet kernel for receiver ir
    // get velocity and coefr on this station
    cordx = {xstore[ibool(ielem_recv[ir],0)],xstore[ibool(ielem_recv[ir],1)]};
    cordy = {ystore[ibool(ielem_recv[ir],0)],ystore[ibool(ielem_recv[ir],3)]};
    bilinear(cordx.data(),cordy.data(),2,2,stalocs(ir,0),stalocs(ir,1),ix,iy,coefr.data());
    v0 = get_velocity(stalocs(ir,0),stalocs(ir,1)); // set v0 to veloc at stations

    // now loop around all connected points to find frechet kernel
    x0 = stalocs(ir,0); y0 = stalocs(ir,1); z0 = stalocs(ir,2);
    int ielem0 = ielem_recv[ir], ipt0 = -1;
    int ielem1 = comming_node_recv(ir,0), ipt1 = comming_node_recv(ir,1);
    while(ipt1 != -1){
        int inode1 = ibool(ielem1,ipt1);
        x1 = xstore[inode1]; y1 = ystore[inode1]; z1 = zstore[inode1];
        v1 = veloc(ielem1,ipt1);
        float dist = compute_length(x0,y0,z0,x1,y1,z1);
        float term = -2.0 * dist / std::pow(v0 + v1,2);
        
        // contribution from inode0/1
        fdm[inode1] += term;
        if(ipt0 == -1){ // at station
            for(int i = 0; i< 4; i++){
                int inode0 = ibool(ielem0,i);
                fdm[inode0] += coefr[i] * term;
            }
        }
        else{
            int inode0 = ibool(ielem0,ipt0);
            fdm[inode0] += term;
        }

        // swap 
        v0 = v1;
        x0 = x1; y0 = y1; z0 = z1;
        ielem0 = ielem1; ipt0 = ipt1;
        ielem1 = comming_node(inode1,0), ipt1 = comming_node(inode1,1);
    }

    // contribution from source segment
    float dist = compute_length(x0,y0,z0,xsource,ysource,zsource);
    float term = -2.0 * dist / std::pow(v0 + vs,2);
    int inode0 = ibool(ielem0,ipt0);
    fdm[inode0] += term;
    for(int i = 0; i < 4 ;i++){
        int inode1 = ibool(ielem_src,i);
        fdm[inode1] += coefs[i] * term;
    }
}