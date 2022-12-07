//                If mode = 1, compute  y = y + A*x.
//                If mode = 2, compute  x = x + A(transpose)*y.
#include "clsqr2/clsqr2.hpp"

/**
 * @brief  perfrom matrix-vector multiplication for CSR matrix
 * 
 * @param mode If mode = 1, compute  y = y + A*x, else compute  x = x + A(transpose)*y.
 * @param m,n shape of the sparse matrix
 * @param x vectorx, shape(n)
 * @param y vector y shape (m)
 * @param val  value for sparse matrix, shape(# of non-zero elements)
 * @param indices shape(# of non-zero elements)
 * @param indptr shape(m+1)
 */
void 
aprod(int mode, int m, int n, real_t* restrict x, real_t* restrict y,
           const real_t* restrict val, const int* restrict indices,
           const int *restrict indptr)
{
    if(mode == 1){
        for(int i = 0; i < m;i++){
        for(int j = indptr[i]; j < indptr[i+1]; j++){
            y[i] += val[j] * x[indices[j]];
        }}
    }
    else{
        for(int i = 0; i < m;i++){
        for(int j = indptr[i]; j < indptr[i+1]; j++){
            x[indices[j]] += val[j] * y[i];
        }}
    }
}