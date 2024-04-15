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
void aprod1(int m,int n,const real_t *x,real_t* restrict y,
           const real_t* val, const int* indices,
           const int* indptr,int nproc)
{   
    // serial code
    if(nproc == 1) {
        for(int i = 0; i < m;i++){
        for(int j = indptr[i]; j < indptr[i+1]; j++){
            y[i] += val[j] * x[indices[j]];
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

    // new ytmp
    real_t *ytmp = new real_t[nproc * m]();
    #pragma omp parallel for shared(ytmp)
    for(int irank = 0; irank < nproc; irank ++) {
        for(int i = irank; i < m; i += nproc) {
        for(int j = indptr[i]; j < indptr[i+1]; j ++ ) {
            ytmp[irank * m + i] += val[j] * x[indices[j]]; 
        }}
    }

    // copy to global
    for(int i = 0; i < m; i ++) {
        double s = 0.;
        for(int irank = 0; irank < nproc; irank ++) {
            s += ytmp[irank * m + i];
        }
        y[i] += s;
    }
    delete [] ytmp;

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
void aprod2(int m,int n,real_t* restrict x,const real_t* y,
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

    // allocate space 
    real_t *xtmp = new real_t[nproc * n]();

    // compute
    #pragma omp parallel for shared(xtmp)
    for(int irank = 0; irank < nproc; irank ++) {
        for(int i = irank; i < n; i += nproc) {
        for(int j = indptr[i]; j < indptr[i+1]; j ++ ) {
            xtmp[irank * n + indices[j]] += val[j] * y[i];
        }}
    }

    // copy to local
    for(int i = 0; i < n; i ++) {
        double s = 0.;
        for(int irank = 0; irank < nproc; irank ++) {
            s += xtmp[irank * n + i];
        }
        x[i] += s;
    }
    delete[] xtmp;

    // set global nprocs back
    omp_set_num_threads(nproc_bak);
}