#ifndef _LSQR_MAIN_H
#define _LSQR_MAIN_H
#include "clsqr_const.hpp"
/* lsqr.h
   $Revision: 229 $ $Date: 2006-04-15 18:40:08 -0700 (Sat, 15 Apr 2006) $
   $Revision: 229 $ $Date: 2022-12-07 by nqdu, C++ version $
*/
/*!
   \file
   Header file for ISO C version of LSQR.
*/
void
lsqr( int m,
      int n,
      const real_t* value,
      const int* indices,
      const int* indptr,
      real_t damp,
      const real_t* b,     // len = m
      real_t* __restrict x,
      real_t atol,
      real_t btol,
      real_t conlim,
      int    itnlim,
      void   *filestream,
      // The remaining variables are output only.
      int    *istop_out,
      int    *itn_out,
      real_t *anorm_out,
      real_t *acond_out,
      real_t *rnorm_out,
      real_t *arnorm_out,
      real_t *xnorm_out,
      int nprocs_used
);
#endif