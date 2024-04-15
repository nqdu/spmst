/* bccblas.c
   $Revision: 231 $ $Date: 2006-04-15 18:47:05 -0700 (Sat, 15 Apr 2006) $

   ----------------------------------------------------------------------
   This file is part of BCLS (Bound-Constrained Least Squares).

   Copyright (C) 2006 Michael P. Friedlander, Department of Computer
   Science, University of British Columbia, Canada. All rights
   reserved. E-mail: <mpf@cs.ubc.ca>.
   
   BCLS is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of the
   License, or (at your option) any later version.
   
   BCLS is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
   Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with BCLS; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
   USA
   ----------------------------------------------------------------------
*/
/*!
   \file

   This file contains C-wrappers to the BLAS (Basic Linear Algebra
   Subprograms) routines that are used by BCLS.  Whenever possible,
   they should be replaced by corresponding BLAS routines that have
   been optimized to the machine being used.

   Included BLAS routines:

   - cblas_daxpy
   - cblas_dcopy
   - cblas_ddot
   - cblas_dnrm2
   - cblas_dscal
*/

#include "cblas.hpp"

/*!
  \param[in]     N
  \param[in]     alpha
  \param[in]     X      
  \param[in]     incX
  \param[in,out] Y
  \param[in]     incY
*/
void
cblas_daxpy( const int N, const real_t alpha, const real_t *X,
             const int incX, real_t* restrict Y, const int incY)
{
  int i;

  if (N     <= 0  ) return;
  if (alpha == 0.0) return;

  if (incX == 1 && incY == 1) {
      const int m = N % 4;

      for (i = 0; i < m; i++)
          Y[i] += alpha * X[i];
      
      for (i = m; i + 3 < N; i += 4) {
          Y[i    ] += alpha * X[i    ];
          Y[i + 1] += alpha * X[i + 1];
          Y[i + 2] += alpha * X[i + 2];
          Y[i + 3] += alpha * X[i + 3];
      }
  } else {
      int ix = OFFSET(N, incX);
      int iy = OFFSET(N, incY);

      for (i = 0; i < N; i++) {
          Y[iy] += alpha * X[ix];
          ix    += incX;
          iy    += incY;
      }
  }
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX
  \param[out]    Y
  \param[in]     incY
*/
void
cblas_dcopy( const int N, const real_t* X,
             const int incX, real_t* restrict Y, const int incY)
{
  int i;
  int ix = OFFSET(N, incX);
  int iy = OFFSET(N, incY);

  for (i = 0; i < N; i++) {
      Y[iy]  = X[ix];
      ix    += incX;
      iy    += incY;
  }
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX
  \param[in]     Y
  \param[in]     incY
  
  \return  Dot product of X and Y.

*/
real_t
cblas_ddot( const int N, const real_t* X,
            const int incX, const real_t* Y, const int incY)
{
  real_t r  = 0.0;
  int    i;
  int    ix = OFFSET(N, incX);
  int    iy = OFFSET(N, incY);

  for (i = 0; i < N; i++) {
      r  += X[ix] * Y[iy];
      ix += incX;
      iy += incY;
  }
  
  return r;
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX

  \return Two-norm of X.
*/
real_t
cblas_dnrm2( const int N, const real_t* X, const int incX) 
{
  real_t
      scale = 0.0,
      ssq   = 1.0;
  int
      i,
      ix    = 0;

  if (N <= 0 || incX <= 0) return 0;
  else if (N == 1)         return fabs(X[0]);

  for (i = 0; i < N; i++) {
      const real_t x = X[ix];

      if (x != 0.0) {
          const real_t ax = fabs(x);

          if (scale < ax) {
              ssq   = 1.0 + ssq * (scale / ax) * (scale / ax);
              scale = ax;
          } else {
              ssq += (ax / scale) * (ax / scale);
          }
      }

      ix += incX;
  }
  
  return scale * sqrt(ssq);
}

/*!
  \param[in]     N
  \param[in]     alpha
  \param[in,out] X
  \param[in]     incX
*/
void
cblas_dscal(const int N, const real_t alpha, real_t* restrict X, const int incX)
{
    int i, ix;

    if (incX <= 0) return;

    ix = OFFSET(N, incX);
    
    for (i = 0; i < N; i++) {
        X[ix] *= alpha;
        ix    += incX;
    }
}
