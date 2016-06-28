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
  FILE:  xf_MathBlasCommon.c

  This file contains home-grown versions of BLAS functions, called
  either when BLAS is not available or when matrices are small enough
  not to warrant BLAS.

*/

#include "xf.h"

/******************************************************************/
//   FUNCTION Definition: xfc_MxM_Set
void
xfc_MxM_Set(const real *A, const real *B, int rA, int n, int cB, real *C)
{
  
  int k, i, col;
  int ik2, ik3, ik;
  real temp, t0, t1, t2, t3;

  for(k=0; k<rA; k++){
    ik3 = cB*k;
    ik2 = n*k;
    i = 0;
    while ((i+3)<cB){
      t0 = 0; t1 = 0; t2 = 0; t3 = 0;
      
      for (col=0; col<n; col++){
	temp = A[ik2+col];
	ik = col*cB+i;
	t0 += temp*B[ik+0];
	t1 += temp*B[ik+1];
	t2 += temp*B[ik+2];
	t3 += temp*B[ik+3]; 
      }
      C[ik3+i+0] = t0;
      C[ik3+i+1] = t1;
      C[ik3+i+2] = t2;
      C[ik3+i+3] = t3;
      i += 4;
    }
    while (i<cB){
      t0 = 0;
      for (col=0; col<n; col++)
	t0 += A[ik2+col]*B[col*cB+i];
      
      C[ik3+i] = t0;
      i += 1;
    }
  }
}

/******************************************************************/
//   FUNCTION Definition: xfc_MxM_Neg
void
xfc_MxM_Neg(const real *A, const real *B, int rA, int n, int cB, real *C)
{
  
  int k, i, col;
  int ik2, ik3, ik;
  real temp, t0, t1, t2, t3;

  for(k=0; k<rA; k++){
    ik3 = cB*k;
    ik2 = n*k;
    i = 0;
    while ((i+3)<cB){
      t0 = 0; t1 = 0; t2 = 0; t3 = 0;
      
      for (col=0; col<n; col++){
	temp = A[ik2+col];
	ik = col*cB+i;
	t0 += temp*B[ik+0];
	t1 += temp*B[ik+1];
	t2 += temp*B[ik+2];
	t3 += temp*B[ik+3]; 
      }
      C[ik3+i+0] = -t0;
      C[ik3+i+1] = -t1;
      C[ik3+i+2] = -t2;
      C[ik3+i+3] = -t3;
      i += 4;
    }
    while (i<cB){
      t0 = 0;
      for (col=0; col<n; col++)
	t0 += A[ik2+col]*B[col*cB+i];
      
      C[ik3+i] = -t0;
      i += 1;
    }
  }
}

/******************************************************************/
//   FUNCTION Definition: xfc_MxM_Add
void
xfc_MxM_Add(const real *A, const real *B, int rA, int n, int cB, real *C)
{
  
  int k, i, col;
  int ik2, ik3, ik;
  real temp, t0, t1, t2, t3;

  for(k=0; k<rA; k++){
    ik3 = cB*k;
    ik2 = n*k;
    i = 0;
    while ((i+3)<cB){
      t0 = 0; t1 = 0; t2 = 0; t3 = 0;
      
      for (col=0; col<n; col++){
	temp = A[ik2+col];
	ik = col*cB+i;
	t0 += temp*B[ik+0];
	t1 += temp*B[ik+1];
	t2 += temp*B[ik+2];
	t3 += temp*B[ik+3]; 
      }
      C[ik3+i+0] += t0;
      C[ik3+i+1] += t1;
      C[ik3+i+2] += t2;
      C[ik3+i+3] += t3;
      i += 4;
    }
    while (i<cB){
      t0 = 0;
      for (col=0; col<n; col++)
	t0 += A[ik2+col]*B[col*cB+i];
      
      C[ik3+i] += t0;
      i += 1;
    }
  }
}


/******************************************************************/
//   FUNCTION Definition: xfc_MxM_Sub
void
xfc_MxM_Sub(const real *A, const real *B, int rA, int n, int cB, real *C)
{
  
  int k, i, col;
  int ik2, ik3, ik;
  real temp, t0, t1, t2, t3;

  for(k=0; k<rA; k++){
    ik3 = cB*k;
    ik2 = n*k;
    i = 0;
    while ((i+3)<cB){
      t0 = 0; t1 = 0; t2 = 0; t3 = 0;
      
      for (col=0; col<n; col++){
	temp = A[ik2+col];
	ik = col*cB+i;
	t0 += temp*B[ik+0];
	t1 += temp*B[ik+1];
	t2 += temp*B[ik+2];
	t3 += temp*B[ik+3]; 
      }
      C[ik3+i+0] -= t0;
      C[ik3+i+1] -= t1;
      C[ik3+i+2] -= t2;
      C[ik3+i+3] -= t3;
      i += 4;
    }
    while (i<cB){
      t0 = 0;
      for (col=0; col<n; col++)
	t0 += A[ik2+col]*B[col*cB+i];
      
      C[ik3+i] -= t0;
      i += 1;
    }
  }
}



/******************************************************************/
//   FUNCTION Definition: xfc_MTxM_Add
static void
xfc_MTxM_Add(const real * RSTRCT A, const real * RSTRCT B, 
	     int cA, int n, int cB, real * RSTRCT C)
{  
  int k, i, col, ik;
  real t0, t1, t2, t3, tB;

  for (i=0; i<cB; i++){
    
    k=0;
    while (k+3<cA){
      t0 = t1 = t2 = t3 = 0.;
      for (col=0; col<n; col++){
	ik = col*cA+k;
	tB = B[col*cB+i];
	t0 += A[ik]*tB;
	t1 += A[ik+1]*tB;
	t2 += A[ik+2]*tB;
	t3 += A[ik+3]*tB;
      }
      ik = cB*k+i;
      C[ik]    += t0;
      C[ik+cB] += t1;
      C[ik+2*cB] += t2;
      C[ik+3*cB] += t3;
      k+=4;
    }
    
    while (k<cA){
      for (col=0, t0=0.; col<n; col++)
	t0 += A[col*cA+k]*B[col*cB+i];
      C[cB*k+i] += t0;
      k++;
    }
  }

}


/******************************************************************/
//   FUNCTION Definition: xfc_MTxM_Sub
static void
xfc_MTxM_Sub(const real * RSTRCT A, const real * RSTRCT B, 
	     int cA, int n, int cB, real * RSTRCT C)
{
  int k, i, col, ik;
  real t0, t1, t2, t3, tB;
  
  //xf_printf("cA = %d, n = %d, cB = %d\n", cA, n, cB);
  for (i=0; i<cB; i++){
    
    k=0;
    while (k+3<cA){
      t0 = t1 = t2 = t3 = 0.;
      for (col=0; col<n; col++){
	ik = col*cA+k;
	tB = B[col*cB+i];
	t0 += A[ik]*tB;
	t1 += A[ik+1]*tB;
	t2 += A[ik+2]*tB;
	t3 += A[ik+3]*tB;
      }
      ik = cB*k+i;
      C[ik]    -= t0;
      C[ik+cB] -= t1;
      C[ik+2*cB] -= t2;
      C[ik+3*cB] -= t3;
      k+=4;
    }
    
    while (k<cA){
      for (col=0, t0=0.; col<n; col++)
	t0 += A[col*cA+k]*B[col*cB+i];
      C[cB*k+i] -= t0;
      k++;
    }
  }
  
}


/******************************************************************/
//   FUNCTION Definition: xfc_MTxM_Set
void
xfc_MTxM_Set(const real *A, const real *B, int cA, int n, int cB, real *C)
{
  int k, i, col;
  real t0;

  for(k=0; k<cA; k++){
    for (i=0; i<cB; i++){
      t0 = 0;
      for (col=0; col<n; col++){
	t0 += A[col*cA+k]*B[col*cB+i];
      }
      C[cB*k+i] = t0;
    }
  }
}



/******************************************************************/
//   FUNCTION Definition: xfc_MTxM_Neg
void
xfc_MTxM_Neg(const real *A, const real *B, int cA, int n, int cB, real *C)
{
  int k, i, col;
  real t0;

  for(k=0; k<cA; k++){
    for (i=0; i<cB; i++){
      t0 = 0;
      for (col=0; col<n; col++){
	t0 += A[col*cA+k]*B[col*cB+i];
      }
      C[cB*k+i] = -t0;
    }
  }
}

/* The functions in this file are unit-tested from the files that
   include it: either xf_MathBlas.c or xf_MathNoBlas.c */
