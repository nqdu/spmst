#include "clsqr2.hpp"

void aprod(int mode, int m, int n, real_t* restrict x, real_t* restrict y,
           const real_t* restrict val, const int* restrict indices,
           const int *restrict indptr);