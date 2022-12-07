#ifndef _LSQR_H
#define _LSQR_H
#include "clsqr2.hpp"
/* lsqr.h
   $Revision: 229 $ $Date: 2006-04-15 18:40:08 -0700 (Sat, 15 Apr 2006) $
*/
/*!
   \file
   Header file for ISO C version of LSQR.
*/
void
lsqr( int m,
      int n,
      const real_t* restrict value,
      const int* restrict indices,
      const int* restrict indptr,
      real_t damp,
      const real_t* restrict b,     // len = m
      real_t* restrict x,
      real_t atol,
      real_t btol,
      real_t conlim,
      int    itnlim,
      void   *nout,
      // The remaining variables are output only.
      int    *istop_out,
      int    *itn_out,
      real_t *anorm_out,
      real_t *acond_out,
      real_t *rnorm_out,
      real_t *arnorm_out,
      real_t *xnorm_out
);
#endif