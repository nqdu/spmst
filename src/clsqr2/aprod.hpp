#ifndef _LSQR_APROD_H
#define _LSQR_APROD_H

#include "clsqr_const.hpp"

void aprod1(int m,int n,const real_t *x,real_t* __restrict y,
           const real_t* val, const int* indices,
           const int* indptr,int nproc);

void aprod2(int m,int n,real_t* __restrict x,const real_t* y,
           const real_t* val, const int* indices,
           const int* indptr,int nproc);
#endif