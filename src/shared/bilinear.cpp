#include "spm2d.hpp"

/**
 * @brief bilinear interpolation to get interp coefs 
 * 
 * @param x/y x,y coordinates, both in ascending order and even space
 * @param n size of x
 * @param x0,y0 interp point
 * @return float 
 */
void bilinear(const float* restrict x, const float* restrict y,int nx,int ny,
                float x0,float y0,int &ix,int &iy,float* restrict coef)
{
    float dx = x[1] - x[0], dy = y[1] - y[0];
    ix = (x0 - x[0]) / dx, iy = (y0 - y[0]) / dy;
    if(ix < 0 ) ix = 0;  
    if(ix > nx-2) ix = nx-2;
    if(iy < 0) iy = 0; 
    if(iy > ny -2) iy = ny - 2;

    float x1 = x[ix], x2  = x[ix+1];
    float y1 = y[iy], y2 = y[iy+1];
    float f1 = (x2 - x0)/ dx, f2 = (x0 - x1) / dx;
    float f3 = (x2 - x0) / dx, f4 = (x0 - x1) / dx;
    float f5 = (y2 - y0) / dy, f6 = (y0 - y1) / dy;
    coef[0] = f1 * f5; coef[1] = f2 * f5;
    coef[2] = f3 * f6; coef[3] = f4 * f6;
}

/**
 * @brief 
 * 
 * @param x,y coordinates, shape(nx) (ny)
 * @param z value for this point, shape(ny,nx)
 * @param nx 
 * @param ny 
 * @param x0,y0 interpolated points 
 * @return float 
 */
float 
interp2d(const float* restrict x, const float* restrict y,
        const float* restrict z,int nx,int ny,float x0,float y0)
{
    float coef[4];
    int ix,iy;
    bilinear(x,y,nx,ny,x0,y0,ix,iy,coef);
    float out = z[iy*nx+ix] * coef[0] + z[iy*nx+ix+1] * coef[1] + 
                z[(iy+1)*nx+ix] * coef[2] + z[(iy+1)*nx+ix+1] * coef[3];
    
    return out ;
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