#ifndef _SPMST_INTERP_H
#define _SPMST_INTERP_H
void bilinear(const float* x, const float* y,int nx,int ny,
                float x0,float y0,int &ix,int &iy,float* __restrict coef);

float interp2d(const float* x, const float* y,
        const float* z,int nx,int ny,float x0,float y0);

#endif