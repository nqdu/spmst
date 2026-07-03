//                If mode = 1, compute  y = y + A*x.
//                If mode = 2, compute  x = x + A(transpose)*y.
#include "clsqr_const.hpp"
#include <omp.h>

/**
 * @brief  perfrom y = y + A * x
 * 
 * @param m,n shape of the sparse matrix
 * @param x vector x, shape(n)
 * @param y vector y shape (m)
 * @param val  value for sparse matrix, shape(# of non-zero elements)
 * @param indices shape(# of non-zero elements)
 * @param indptr shape(m+1)
 * @param nproc # of procs used
 */
void aprod1(int m,int n,const real_t *x,real_t* __restrict y,
           const real_t* val, const int* indices,
           const int* indptr,int nproc)
{   
    // backup global omp
    int nproc_bak = 1;
    #pragma omp parallel 
    {
        nproc_bak = omp_get_num_threads();
    }
    omp_set_num_threads(nproc);

    #pragma omp parallel for shared(indptr,indices,x,y,val)
    for(int i = 0; i < m;i++){
        for(int j = indptr[i]; j < indptr[i+1]; j++){
            y[i] += val[j] * x[indices[j]];
        }
    }

    // set global nprocs back
    omp_set_num_threads(nproc_bak);
}

/**
 * @brief  perfrom x = x + A.T * y
 * 
 * @param m,n shape of the sparse matrix
 * @param x vector x, shape(n)
 * @param y vector y shape (m)
 * @param val  value for sparse matrix, shape(# of non-zero elements)
 * @param indices shape(# of non-zero elements)
 * @param indptr shape(m+1)
 * @param nproc # of procs used
 */
void aprod2(int m,int n,real_t* __restrict x,const real_t* y,
           const real_t* val, const int* indices,
           const int* indptr,int nproc)
{
    if(nproc == 1) {
        for(int i = 0; i < m;i++){
        for(int j = indptr[i]; j < indptr[i+1]; j++){
            x[indices[j]] += val[j] * y[i];
        }}

        return;
    }

    // backup global omp
    int nproc_bak = 1;
    #pragma omp parallel 
    {
        nproc_bak = omp_get_num_threads();
    }
    omp_set_num_threads(nproc);

    // atomic add 
    #pragma omp parallel for shared(indptr,indices,x,y,val)
    for(int i = 0; i < m;i++) {
    for(int j = indptr[i]; j < indptr[i+1]; j++) {
        #pragma omp atomic
        x[indices[j]] += val[j] * y[i];
    }}

    // set global nprocs back
    omp_set_num_threads(nproc_bak);
}