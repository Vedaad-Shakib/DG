/*------------------------------------------------------------------*/
/* XFLOW: A discontinuous Galerkin finite element software library. */
/*                                                                  */
/*                    Copyright  2007-2008                          */
/*           Krzysztof J. Fidkowski, kfid@alum.mit.edu              */
/*                                                                  */
/*                    Copyright  2008-2012                          */
/*                 The University of Michigan                       */
/*                    All rights reserved                           */
/*                                                                  */
/* This library is intended to be useful but is distributed without */
/* any warranty, not even merchantability or fitness for a          */
/* particular purpose.  It is free software: you can redistribute   */
/* it and/or modify it under the terms of the GNU Lesser General    */
/* Public License (LGPLv3).                                         */
/*                                                                  */
/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free       */
/* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.        */
/*------------------------------------------------------------------*/

/*
  FILE:  xf_MathBlas.c

  This file contains top-level wrappers for math functions that are
  implemented in BLAS.

*/

#include "xf.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"


// critical matrix size (row*col) for using BLAS
#define xf_USEBLAS 1000

/* BLAS prototypes */
// SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
void *dgemm_(char *, char *, int *, int *, int *, const real *, const real *, int *, 
	     const real *, int *, const real *, real *, int *);


#include "xf_MathBlasCommon.c"

/******************************************************************/
//   FUNCTION Definition: xf_MxM_Add
void
xf_MxM_Add(const real *A, const real *B, int rA, int n, int cB, real *C)
{
  // In c, we want:    C   = +A*B + C
  // in fortran   :    C^T = +B^T*A^T + C^T
  char transA[] = "N"; // no transpose for A (recall Fortran is column-major)
  char transB[] = "N"; // no transpose for B
  real alpha =  1.0;
  real beta  =  1.0;

  if (max(rA*n, cB*n) > xf_USEBLAS){
    //DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    //  C := alpha*op( A )*op( B ) + beta*C,
    dgemm_(transB, transA, &cB, &rA, &n, &alpha, B, &cB, A, &n, &beta, C, &cB);
  }
  else{
    xfc_MxM_Add(A, B, rA, n, cB, C);
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_MxM_Sub
void
xf_MxM_Sub(const real *A, const real *B, int rA, int n, int cB, real *C)
{
  // In c, we want:    C   = -A*B + C
  // in fortran   :    C^T = -B^T*A^T + C^T
  char transA[] = "N"; // no transpose for A (recall Fortran is column-major)
  char transB[] = "N"; // no transpose for B
  real alpha = -1.0;
  real beta  =  1.0;

  if (max(rA*n, cB*n) > xf_USEBLAS){
    //DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    //  C := alpha*op( A )*op( B ) + beta*C,
    dgemm_(transB, transA, &cB, &rA, &n, &alpha, B, &cB, A, &n, &beta, C, &cB);
  }
  else{
    xfc_MxM_Sub(A, B, rA, n, cB, C);
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_MxM_Set
void
xf_MxM_Set(const real *A, const real *B, int rA, int n, int cB, real *C)
{
  // In c, we want:    C   = A*B
  // in fortran   :    C^T = B^T*A^T
  char transA[] = "N"; // no transpose for A (recall Fortran is column-major)
  char transB[] = "N"; // no transpose for B
  real alpha =  1.0;
  real beta  =  0.0;

  if (max(rA*n, cB*n) > xf_USEBLAS){
    //DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    //  C := alpha*op( A )*op( B ) + beta*C,
    dgemm_(transB, transA, &cB, &rA, &n, &alpha, B, &cB, A, &n, &beta, C, &cB);
  }
  else{
    xfc_MxM_Set(A, B, rA, n, cB, C);
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_MxM_Neg
void
xf_MxM_Neg(const real *A, const real *B, int rA, int n, int cB, real *C)
{
  // In c, we want:    C   = -A*B
  // in fortran   :    C^T = -B^T*A^T
  char transA[] = "N"; // no transpose for A (recall Fortran is column-major)
  char transB[] = "N"; // no transpose for B
  real alpha = -1.0;
  real beta  =  0.0;

  if (max(rA*n, cB*n) > xf_USEBLAS){
    //DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    //  C := alpha*op( A )*op( B ) + beta*C,
    dgemm_(transB, transA, &cB, &rA, &n, &alpha, B, &cB, A, &n, &beta, C, &cB);
  }
  else{
    xfc_MxM_Neg(A, B, rA, n, cB, C);
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_MTxM_Add
void
xf_MTxM_Add(const real * RSTRCT A, const real * RSTRCT B, 
	    int cA, int n, int cB, real * RSTRCT C)
{
  // In c, we want:    C   = +A^T*B + C
  // in fortran   :    C^T = +B^T*A + C^T

  char transA[] = "T"; // transpose for A (recall Fortran is column-major)
  char transB[] = "N"; // no transpose for B
  real alpha =  1.0;
  real beta  =  1.0;

  if (max(cA*n, cB*n) > xf_USEBLAS){
    //DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    //  C := alpha*op( A )*op( B ) + beta*C,
    dgemm_(transB, transA, &cB, &cA, &n, &alpha, B, &cB, A, &cA, &beta, C, &cB);
  }
  else{
    xfc_MTxM_Add(A, B, cA, n, cB, C);
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_MTxM_Sub
void
xf_MTxM_Sub(const real * RSTRCT A, const real * RSTRCT B, 
	    int cA, int n, int cB, real * RSTRCT C)
{
  // In c, we want:    C   = -A^T*B + C
  // in fortran   :    C^T = -B^T*A + C^T

  char transA[] = "T"; // transpose for A (recall Fortran is column-major)
  char transB[] = "N"; // no transpose for B
  real alpha = -1.0;
  real beta  =  1.0;

  if (max(cA*n, cB*n) > xf_USEBLAS){
    //DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    //  C := alpha*op( A )*op( B ) + beta*C,
    dgemm_(transB, transA, &cB, &cA, &n, &alpha, B, &cB, A, &cA, &beta, C, &cB);
  }
  else{
    xfc_MTxM_Sub(A, B, cA, n, cB, C);
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_MTxM_Set
void
xf_MTxM_Set(const real *A, const real *B, int cA, int n, int cB, real *C)
{
  // In c, we want:    C   =  A^T*B
  // in fortran   :    C^T =  B^T*A

  char transA[] = "T"; // transpose for A (recall Fortran is column-major)
  char transB[] = "N"; // no transpose for B
  real alpha =  1.0;
  real beta  =  0.0;

  if (max(cA*n, cB*n) > xf_USEBLAS){
    //DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    //  C := alpha*op( A )*op( B ) + beta*C,
    dgemm_(transB, transA, &cB, &cA, &n, &alpha, B, &cB, A, &cA, &beta, C, &cB);
  }
  else{
    xfc_MTxM_Set(A, B, cA, n, cB, C);
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_MTxM_Neg
void
xf_MTxM_Neg(const real *A, const real *B, int cA, int n, int cB, real *C)
{
  // In c, we want:    C   = -A^T*B
  // in fortran   :    C^T = -B^T*A

  char transA[] = "T"; // transpose for A (recall Fortran is column-major)
  char transB[] = "N"; // no transpose for B
  real alpha = -1.0;
  real beta  =  0.0;
  
  if (max(cA*n, cB*n) > xf_USEBLAS){
    //DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    //  C := alpha*op( A )*op( B ) + beta*C,
    dgemm_(transB, transA, &cB, &cA, &n, &alpha, B, &cB, A, &cA, &beta, C, &cB);
  }
  else{
    xfc_MTxM_Neg(A, B, cA, n, cB, C);
  }
}


#if( UNIT_TEST==1 )
#include "xf_MathBlas.test.in"
#endif
