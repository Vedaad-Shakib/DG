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
  FILE:  xf_Math.c

  This file contains math functions.

*/

#include "xf.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "xf_Dynamic.h"
#include "xf_MathBlas.h"


/******************************************************************/
//   FUNCTION Definition: xf_GetAddFlag2
enum xfe_AddType
xf_GetAddFlag2(enum xfe_AddType AddFlag)
{ 
  // always return either add or sub
  if ((AddFlag == xfe_Add) || (AddFlag == xfe_Set)) 
    return xfe_Add;
  else
    return xfe_Sub;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetAddFlagNeg
enum xfe_AddType
xf_GetAddFlagNeg(enum xfe_AddType AddFlag)
{  
  switch(AddFlag){
  case xfe_Set: return xfe_Neg; break;
  case xfe_Neg: return xfe_Set; break;
  case xfe_Add: return xfe_Sub; break;
  case xfe_Sub: return xfe_Add; break;
  default:      return xfe_Add; break;
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_V_Add
void 
xf_V_Add(const real *u, int n, enum xfe_AddType AddFlag, real *v)
{
  int k;

  switch(AddFlag){
  case xfe_Set: for (k=0; k<n; k++) v[k]  =  u[k]; break;
  case xfe_Neg: for (k=0; k<n; k++) v[k]  = -u[k]; break;
  case xfe_Add: for (k=0; k<n; k++) v[k] +=  u[k]; break;
  case xfe_Sub: for (k=0; k<n; k++) v[k] -=  u[k]; break;
  default: xf_Error(xf_INPUT_ERROR); break;
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_cV_Add
void 
xf_cV_Add(const real *u, real c, int n, enum xfe_AddType AddFlag, real *v)
{
  int k;

  switch(AddFlag){
  case xfe_Set: for (k=0; k<n; k++) v[k]  =  c*u[k]; break;
  case xfe_Neg: for (k=0; k<n; k++) v[k]  = -c*u[k]; break;
  case xfe_Add: for (k=0; k<n; k++) v[k] +=  c*u[k]; break;
  case xfe_Sub: for (k=0; k<n; k++) v[k] -=  c*u[k]; break;
  default: xf_Error(xf_INPUT_ERROR); break;
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_PowInt
real
xf_PowInt(real x, int n)
{
  int k;
  real t;
  if (n < 0)
    t = 0.;
  else
    for (k=0, t=1.; k<n; k++) t *= x;
  return t;
}


/******************************************************************/
//   FUNCTION Definition: xf_Distance
real
xf_Distance(real *x0, real *x1, int dim)
{
  int i;
  real d2;

  for (i=0, d2=0.; i<dim; i++)
    d2 += (x1[i]-x0[i])*(x1[i]-x0[i]);
  
  return sqrt(d2);

}

/******************************************************************/
//   FUNCTION Definition: xf_MxV_Add
void 
xf_MxV_Add(const real *A, const real *u, int rA, int cA, real *v)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rA, cA, 1.e0, A, cA, u, 1, 1.e0, v, 1);
    /*int k, ik, j;
  real t;

  ik = 0;
  for (k=0; k<rA; k++){
    t = 0;
    for (j=0; j<cA; j++)
      t += A[ik+j]*u[j];
    v[k] += t;
    ik += cA;
    }*/
}


/******************************************************************/
//   FUNCTION Definition: xf_MxV_Sub
void 
xf_MxV_Sub(const real *A, const real *u, int rA, int cA, real *v)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rA, cA, -1.e0, A, cA, u, 1, 1.e0, v, 1);
    /*
  int k, ik, j;
  real t;

  ik = 0;
  for (k=0; k<rA; k++){
    t = 0;
    for (j=0; j<cA; j++)
      t += A[ik+j]*u[j];
    v[k] -= t;
    ik += cA;
  }
    */
}


/******************************************************************/
//   FUNCTION Definition: xf_MxV_Set
void 
xf_MxV_Set(const real *A, const real *u, int rA, int cA, real *v)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rA, cA, 1.e0, A, cA, u, 1, 0, v, 1);
    /*int k, ik, j;
  real t;

  ik = 0;
  for (k=0; k<rA; k++){
    t = 0;
    for (j=0; j<cA; j++)
      t += A[ik+j]*u[j];
    v[k] = t;
    ik += cA;
    }*/
}


/******************************************************************/
//   FUNCTION Definition: xf_MxV_Neg
void 
xf_MxV_Neg(const real *A, const real *u, int rA, int cA, real *v)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rA, cA, -1.e0, A, cA, u, 1, 0, v, 1);
    /*int k, ik, j;
  real t;

  ik = 0;
  for (k=0; k<rA; k++){
    t = 0;
    for (j=0; j<cA; j++)
      t -= A[ik+j]*u[j];
    v[k] = t;
    ik += cA;
    }*/
}



/******************************************************************/
//   FUNCTION Definition: xf_MxV
void
xf_MxV(const real *A, const real *u, int rA, int cA, 
       enum xfe_AddType AddFlag, real *v)
{
  switch (AddFlag){
  case xfe_Set: xf_MxV_Set(A, u, rA, cA, v);  break;
  case xfe_Add: xf_MxV_Add(A, u, rA, cA, v);  break;
  case xfe_Sub: xf_MxV_Sub(A, u, rA, cA, v);  break;
  case xfe_Neg: xf_MxV_Neg(A, u, rA, cA, v);  break;
  default: xf_Error(xf_INPUT_ERROR); break;
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_cMxV_Add
void 
xf_cMxV_Add(real c, real *A, real *u, int rA, int cA, real *v)
{
  int k, ik, j;
  real t;

  ik = 0;
  for (k=0; k<rA; k++){
    t = 0;
    for (j=0; j<cA; j++)
      t += A[ik+j]*u[j];
    v[k] += c*t;
    ik += cA;
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_MTxV_Set
static void 
xf_MTxV_Set(const real *A, const real *u, int cA, int rA, real *v)
{
  int k, j;
  real t;

  for (k=0; k<cA; k++){
    t = 0;
    for (j=0; j<rA; j++)
      t += A[j*cA+k]*u[j];
    v[k] = t;
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_MTxV_Add
static void 
xf_MTxV_Add(const real *A, const real *u, int cA, int rA, real *v)
{
  int k, j;
  real t;

  for (k=0; k<cA; k++){
    t = 0;
    for (j=0; j<rA; j++)
      t += A[j*cA+k]*u[j];
    v[k] += t;
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_MTxV_Sub
static void 
xf_MTxV_Sub(const real *A, const real *u, int cA, int rA, real *v)
{
  int k, j;
  real t;

  for (k=0; k<cA; k++){
    t = 0;
    for (j=0; j<rA; j++)
      t += A[j*cA+k]*u[j];
    v[k] -= t;
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_MTxV
void
xf_MTxV(const real *A, const real *u, int cA, int rA, 
	enum xfe_AddType AddFlag, real *v)
{
  switch (AddFlag){
  case xfe_Set: xf_MTxV_Set(A, u, cA, rA, v);  break;
  case xfe_Add: xf_MTxV_Add(A, u, cA, rA, v);  break;
  case xfe_Sub: xf_MTxV_Sub(A, u, cA, rA, v);  break;
  default: xf_Error(xf_INPUT_ERROR); break;
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_ndMTxVc
void
xf_ndMTxVc(int n, const real *A, const real *f, int d, real *u)
{
  // u{i,:} = A{i,:,:}*f{i}*u{i,:},  where (:) = [0,...,d-1]
  int i, j;
  real t, t1[3], t2[3];
  
  if (d>3) xf_printf("Warning, d=%d but only d<=3 is possible in xf_ndMxMT\n", d);
  for (i=0; i<n; i++){
    for (j=0; j<d; j++) t1[j] = u[i*d+j];
    xf_MTxV_Set(A+i*d*d, t1, d, d, t2);
    for (j=0, t=f[i]; j<d; j++) u[i*d+j] = t*t2[j];
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_ndMTxVic
void
xf_ndMTxVic(int n, const real *A, const real *f, int d, real *u)
{
  // u{i,:} = A{i,:,:}*(1/f{i})*u{i,:},  where (:) = [0,...,d-1]
  int i, j;
  real t, t1[3], t2[3];
  
  if (d>3) xf_printf("Warning, d=%d but only d<=3 is possible in xf_ndMxMT\n", d);
  for (i=0; i<n; i++){
    for (j=0; j<d; j++) t1[j] = u[i*d+j];
    xf_MTxV_Set(A+i*d*d, t1, d, d, t2);
    for (j=0, t=1./f[i]; j<d; j++) u[i*d+j] = t*t2[j];
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_MxM
static void
xf_MxM(const real *A, const real *B, int rA, int n, int cB, 
       enum xfe_AddType AddFlag, real *C)
{
  switch (AddFlag){
  case xfe_Set: xf_MxM_Set(A, B, rA, n, cB, C);  break;
  case xfe_Neg: xf_MxM_Neg(A, B, rA, n, cB, C);  break;
  case xfe_Add: xf_MxM_Add(A, B, rA, n, cB, C);  break;
  case xfe_Sub: xf_MxM_Sub(A, B, rA, n, cB, C);  break;
  default: xf_Error(xf_INPUT_ERROR);
  }

}

/******************************************************************/
//   FUNCTION Definition: xf_nMxM_Set
void
xf_nMxM_Set(int n, const real *A, const real *B, int rA, int cA, 
	    int cB, real *C)
{
  int i, Asize, Bsize, Csize;

  Asize = rA*cA;  Bsize = cA*cB;  Csize = rA*cB;

  for (i=0; i<n; i++)
    xf_MxM_Set(A+i*Asize, B+i*Bsize, rA, cA, cB, C+i*Csize);  
}

/******************************************************************/
//   FUNCTION Definition: xf_nMxM_Add
void
xf_nMxM_Add(int n, const real *A, const real *B, int rA, int cA, 
	    int cB, real *C)
{
  int i, Asize, Bsize, Csize;

  Asize = rA*cA;  Bsize = cA*cB;  Csize = rA*cB;

  for (i=0; i<n; i++)
    xf_MxM_Add(A+i*Asize, B+i*Asize, rA, cA, cB, C+i*Csize);  
}



/******************************************************************/
//   FUNCTION Definition: xf_cMxM_Add
void
xf_cMxM_Add(real c, real *A, real *B, int rA, int n, int cB, real *C)
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
        ik = cB*col;
	ik += i;
        t0 += temp*B[ik+0];
        t1 += temp*B[ik+1];
        t2 += temp*B[ik+2];
        t3 += temp*B[ik+3];
      }
      C[ik3+i+0] += c*t0;
      C[ik3+i+1] += c*t1;
      C[ik3+i+2] += c*t2;
      C[ik3+i+3] += c*t3;
      i += 4;
    }
    while (i<cB){
      t0 = 0;
      for (col=0; col<n; col++)
	t0 += A[ik2+col]*B[col*cB+i];
      
      C[ik3+i] += c*t0;
      i += 1;
    }
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_MTxM
void
xf_MTxM(const real *A, const real *B, int cA, int n, int cB, 
	enum xfe_AddType AddFlag, real *C)
{
  switch (AddFlag){
  case xfe_Set: xf_MTxM_Set(A, B, cA, n, cB, C);  break;
  case xfe_Neg: xf_MTxM_Neg(A, B, cA, n, cB, C);  break;
  case xfe_Add: xf_MTxM_Add(A, B, cA, n, cB, C);  break;
  case xfe_Sub: xf_MTxM_Sub(A, B, cA, n, cB, C);  break;
  default: xf_Error(xf_INPUT_ERROR);
  }
}




/******************************************************************/
//   FUNCTION Definition: xf_MTxwM_Set
void
xf_MTxwM_Set(const real *A, const real *w, const real *B, int cA, int n, int cB, real *C)
{
  int k, i, col;
  real t0;

  for(k=0; k<cA; k++){
    for (i=0; i<cB; i++){
      t0 = 0;
      for (col=0; col<n; col++){
	t0 += w[col]*A[col*cA+k]*B[col*cB+i];
      }
      C[cB*k+i] = t0;
    }
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_MTxwM_Add
void
xf_MTxwM_Add(const real *A, const real *w, const real *B, int cA, int n, int cB, real *C)
{
  int k, i, col;
  real t0;

  for(k=0; k<cA; k++){
    for (i=0; i<cB; i++){
      t0 = 0;
      for (col=0; col<n; col++){
	t0 += A[col*cA+k]*B[col*cB+i]*w[col];
      }
      C[cB*k+i] += t0;
    }
  }
}




/******************************************************************/
//   FUNCTION Definition: xf_MxMT_Set
void
xf_MxMT_Set(const real *A, const real *B, int rA, int n, int rB, real *C)
{
  int k, i, col, kn, in;
  real t0;

  for(k=0; k<rA; k++){
    kn = k*n;
    for (i=0; i<rB; i++){
      in = i*n;
      t0 = 0.;
      for (col=0; col<n; col++){
	t0 += A[kn+col]*B[in+col];
      }
      C[rB*k+i] = t0;
    }
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_ColMult
void
xf_ColMult(real *A, const real *v, int rA, int cA, int dv)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv];
    for (c=0; c<cA; c++){
      A[k+c] *= t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_ColDiv
void
xf_ColDiv(real *A, const real *v, int rA, int cA, int dv)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv];
    for (c=0; c<cA; c++){
      A[k+c] /= t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_ColMult_Set
void
xf_ColMult_Set(const real *A, const  real *v, int rA, int cA, int dv, 
	       real *B)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv];
    for (c=0; c<cA; c++){
      B[k+c] = A[k+c]*t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_ColMult_Add
void
xf_ColMult_Add(const real *A, const real *v, int rA, int cA, int dv, 
	       real *B)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv];
    for (c=0; c<cA; c++){
      B[k+c] += A[k+c]*t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_ColMult_Sub
void
xf_ColMult_Sub(const real *A, const real *v, int rA, int cA, int dv, 
	       real *B)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv];
    for (c=0; c<cA; c++){
      B[k+c] -= A[k+c]*t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_ColcMult
void
xf_ColcMult(real *A, const real *v, int rA, int cA, int dv, real f)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv]*f;
    for (c=0; c<cA; c++){
      A[k+c] *= t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_ColcMult_Set
void
xf_ColcMult_Set(const real *A, const real *v, int rA, int cA, int dv, 
		real f, real *B)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv]*f;
    for (c=0; c<cA; c++){
      B[k+c] = A[k+c]*t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_ColcMult_Add
void
xf_ColcMult_Add(const real *A, const real *v, int rA, int cA, int dv, 
		real f, real *B)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv]*f;
    for (c=0; c<cA; c++){
      B[k+c] += A[k+c]*t;
    } // c
  } // r
}


/******************************************************************/
//   FUNCTION Definition: xf_2ColMult_Set
void
xf_2ColMult_Set(const real *A, const real *v, const real *w, int rA, 
		int cA, int dv, int dw, real *B)
{
  int r, c, k;
  real t;
  
  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv]*w[r*dw];
    for (c=0; c<cA; c++){
      B[k+c] = A[k+c]*t;
    } // c
  } // r
}


/******************************************************************/
//   FUNCTION Definition: xf_2ColMult_Add
void
xf_2ColMult_Add(const real *A, const real *v, const real *w, int rA, 
		int cA, int dv, int dw, real *B)
{
  int r, c, k;
  real t;
  
  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv]*w[r*dw];
    for (c=0; c<cA; c++){
      B[k+c] += A[k+c]*t;
    } // c
  } // r
}


/******************************************************************/
//   FUNCTION Definition: xf_2ColcMult_Set
void
xf_2ColcMult_Set(const real *A, const real *v, const real *w, int rA, 
		 int cA, int dv, int dw, real f, real *B)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv]*w[r*dw]*f;
    for (c=0; c<cA; c++){
      B[k+c] = A[k+c]*t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_2ColcMult_Add
void
xf_2ColcMult_Add(const real *A, const real *v, const real *w, int rA, 
		 int cA, int dv, int dw, real f, real *B)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv]*w[r*dw]*f;
    for (c=0; c<cA; c++){
      B[k+c] += A[k+c]*t;
    } // c
  } // r
}

/******************************************************************/
//   FUNCTION Definition: xf_3ColcMult_Add
void
xf_3ColcMult_Add(const real *A, const real *v, const real *w, 
		 const real *u, int rA, int cA, int dv, int dw, 
		 int du, real f, real *B)
{
  int r, c, k;
  real t;

  for (r=0; r<rA; r++){
    k = r*cA;
    t = v[r*dv]*w[r*dw]*u[r*du]*f;
    for (c=0; c<cA; c++){
      B[k+c] += A[k+c]*t;
    } // c
  } // r
}



/******************************************************************/
//   FUNCTION Definition: xf_MatDetInv
int
xf_MatDetInv(real *A, int d, real *detA, real *iA)
{
  real det;

  if (d == 1){
    det = A[0];
    if (detA != NULL) (*detA) = det;
    if (iA != NULL){
      if (det == 0.) return xf_Error(xf_SINGULAR);
      iA[0] =  1./det;
    }
  }
  else if (d == 2){
    det = A[0]*A[3] - A[2]*A[1];
    if (detA != NULL) (*detA) = det;
    if (iA != NULL){
      if (det == 0.) return xf_Error(xf_SINGULAR);
      iA[0] =  A[3]/det;  iA[1] = -A[1]/det;
      iA[2] = -A[2]/det;  iA[3] =  A[0]/det;
    }
  }
  else if (d == 3){
    det = A[0]*A[4]*A[8] + A[1]*A[5]*A[6] + A[2]*A[3]*A[7]
      -A[6]*A[4]*A[2] - A[7]*A[5]*A[0] - A[8]*A[3]*A[1];
    
    if (detA != NULL) (*detA) = det;

    if (iA != NULL){
      if (det == 0.) return xf_Error(xf_SINGULAR);
      iA[0] = (A[4]*A[8]-A[5]*A[7])/det;
      iA[1] = (A[7]*A[2]-A[8]*A[1])/det;
      iA[2] = (A[1]*A[5]-A[2]*A[4])/det;
      iA[3] = (A[5]*A[6]-A[3]*A[8])/det;
      iA[4] = (A[8]*A[0]-A[6]*A[2])/det;
      iA[5] = (A[2]*A[3]-A[0]*A[5])/det;
      iA[6] = (A[3]*A[7]-A[4]*A[6])/det;
      iA[7] = (A[6]*A[1]-A[7]*A[0])/det;
      iA[8] = (A[0]*A[4]-A[1]*A[3])/det;
    }
  }
  else
    return xf_Error(xf_OUT_OF_BOUNDS); // d must be 2 or 3
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_dMxM
void
xf_dMxM(const real *D, real f, int d, int n, int offset, real *A)
{
  int i, j;
  real t1[3], t2[3];

  if (d>3) xf_printf("Warning, d=%d but only d<=3 is possible in xf_dMxM\n", d);

  for (i=0; i<n; i++){
    for (j=0; j<d; j++) t1[j] = A[i+j*offset];
    xf_MxV_Set(D, t1, d, d, t2);
    for (j=0; j<d; j++) A[i+j*offset] = f*t2[j];
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_dMxMT
void
xf_dMxMT(const real *D, real f, int d, int n, int offset, real *A)
{
  int i, j;
  real t1[3], t2[3];

  if (d>3) xf_printf("Warning, d=%d but only d<=3 is possible in xf_dMxMT\n", d);

  for (i=0; i<n; i++){
    for (j=0; j<d; j++) t1[j] = A[i+j*offset];
    xf_MTxV_Set(D, t1, d, d, t2);
    for (j=0; j<d; j++) A[i+j*offset] = f*t2[j];
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_DotProduct
void
xf_DotProduct(const real *a, const real *b, int n, real *dp)
{
  int i;
  real t;

  for (i=0, t=0.; i<n; i++) t += a[i]*b[i];
  (*dp) = t;
}


/******************************************************************/
//   FUNCTION Definition: xf_CrossProduct
void
xf_CrossProduct(real *a, real *b, real *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}


/******************************************************************/
//   FUNCTION Definition: xf_ComputePLU
int 
xf_ComputePLU(real * RSTRCT A, int r, int * RSTRCT P)
{

  int k, j, i, kmax, tk;
  int kr, jr, ik;
  real l, Akk, Amax, Amax0, tmp, t0, t1, t2, t3;

  if (r == 1){
    P[0] = 0;
    if (fabs(A[0]) < 0.0) return xf_Error(xf_SINGULAR);
    return xf_OK;
  }


  /* initialize P vector */
  for (k=0; k<r; k++) P[k] = k;

  if ((r == 1) && (fabs(A[0]) <= 0.0)) return xf_Error(xf_SINGULAR);

  /* Scale A by maximum entry (for better conditioning) */
  Amax0 = 0.0;
  for (k=0; k<(r*r); k++)
    if (fabs(A[k]) > Amax0) Amax0 = fabs(A[k]);

  if (Amax0 < 0.0) return xf_Error(xf_SINGULAR);

  if (Amax0 > sqrt(MEPS)) Amax0 = 1.0; // no need to scale if entries reasonable

  for (k=0; k<(r*r); k++) A[k] /= Amax0;
  
  for (k=0; k<(r-1); k++){     /*loop through columns*/
    kr = k*r;

    /* find maximum (absolute value) entry below pivot */
    Amax = fabs(A[kr + k]);
    kmax = k;
    for (i=(k+1); i<r; i++){
      if (fabs(A[i*r + k]) > Amax){
        Amax = fabs(A[i*r + k]);
        kmax = i;
      }
    }

    /* switch rows k and kmax if necessary */
    if (kmax != k){
      for (i=0; i<r; i++){
        tmp = A[kr + i]; A[kr + i] = A[kmax*r + i]; A[kmax*r + i] = tmp;
      }
      tk = P[k]; P[k] = P[kmax]; P[kmax] = tk;
    }

    Akk = A[kr + k];
    if (fabs(Akk) < 0.0) return xf_Error(xf_SINGULAR);

    for (j=(k+1); j<r; j++){
      jr = j*r;
      A[jr+k] = A[jr + k]/Akk;
    }

    i = k+1;
    while ((i+3)<r){
      ik = kr+i;
      t0 = A[ik+0]; t1 = A[ik+1]; t2 = A[ik+2]; t3 = A[ik+3];
      for (j=(k+1); j<r; j++){
	jr = j*r+i;
	l = A[j*r+k];
	A[jr+0] -= l*t0; A[jr+1] -= l*t1; A[jr+2] -= l*t2; A[jr+3] -= l*t3;
      }
      i += 4;
    }
    while (i<r){
      ik = kr+i;
      t0 = A[ik];
      for (j=(k+1); j<r; j++){
	jr = j*r;
	A[jr+i] -= A[jr+k]*t0;
      }
      i += 1;
    }

  }

  /* Scale A back; note, since A = PLU, only scale U */
  for (i=0; i<r; i++)
    for (j=i; j<r; j++)
      A[i*r+j] *= Amax0;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SolveU
int 
xf_SolveU(real *U, int r, real *u)
{
  int j, k, ik;
  real t, Ukk;

  if (r == 1){
    u[0] /= U[0];
    return xf_OK;
  }

  /* Solve U*u=u */
  ik = (r-1)*r;
  for (k=(r-1); k>=0; k--){
    t = u[k];
    for (j=(k+1); j<r; j++)
      t -= U[ik + j]*u[j];
    Ukk = U[ik + k];
    if (fabs(Ukk) < 0.0) return xf_Error(xf_SINGULAR);
    u[k] = t/U[ik + k];
    ik -= r;
  }
  return xf_OK;
} 


/******************************************************************/
//   FUNCTION Definition: xf_SolvePLU
int 
xf_SolvePLU(real * RSTRCT A, int * RSTRCT P, real * RSTRCT b, int r, 
	    real * RSTRCT u, real * RSTRCT yin)
{
  // yin is optional
  int ierr;
  int j, k, ik;
  real t, *y;

  if (r == 1){
    u[0] = b[0]/A[0];
    return xf_OK;
  }

  if (yin != NULL)
    y = yin;
  else{
    ierr = xf_Error(xf_Alloc((void **) &y, r, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  /* Solve Ly=Pb */
  ik = 0;
  for (k=0; k<r; k++){
    t = ((P == NULL) ? b[k] :  b[P[k]]);
    for (j=0; j<k; j++)
      t -= A[ik+j]*y[j];
    y[k] = t;
    ik += r;
  }

  /* Solve Uu=y */
  ik = (r-1)*r;
  for (k=(r-1); k>=0; k--){
    t = y[k];
    for (j=(k+1); j<r; j++){
      t -= A[ik + j]*u[j];
    }
    u[k] = t/A[ik + k];
    ik -= r;
  }

  if (yin == NULL)
    xf_Release((void *) y);

  return xf_OK;
} 


/******************************************************************/
//   FUNCTION Definition: xf_SolvePLUT
static int 
xf_SolvePLUT(real *A, int *P, real *b, int r, real *u, real *yin)
{
  // yin is optional
  int ierr;
  int j, k, ik;
  real t, *y;

  if (r == 1){
    u[0] = b[0]/A[0];
    return xf_OK;
  }

  if (yin != NULL)
    y = yin;
  else{
    ierr = xf_Error(xf_Alloc((void **) &y, r, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  /* Solve U^T y = b */
  for (k=0; k<r; k++){
    t = b[k];
    for (j=0, ik=0; j<k; j++, ik+=r)
      t -= A[ik+k]*y[j];
    y[k] = t/A[k*r+k];
  }

  /* Solve L^T (P^T u) = y */
  for (k=(r-1); k>=0; k--){
    t = y[k];
    if (P == NULL){
      for (j=(k+1), ik=r*(k+1); j<r; j++, ik+=r)
	t -= A[ik + k]*u[j];
      u[k] = t;
    }
    else{
      for (j=(k+1), ik=r*(k+1); j<r; j++, ik+=r)
	t -= A[ik + k]*u[P[j]];
      u[P[k]] = t;
    }
  }

  if (yin == NULL)
    xf_Release((void *) y);

  return xf_OK;
} 


/******************************************************************/
//   FUNCTION Definition: xf_SolvePLU_Matrix
int
xf_SolvePLU_Matrix(real *LU, int *P, int r, int cB, real *B)
{
  int ierr, i, j, k;
  real *C, t;

  ierr = xf_Error(xf_Alloc((void **) &C, r*cB, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // Solve LC = PB
  for (i=0; i<cB; i++){ // loop over columns of B
    for (k=0; k<r; k++){
      t = B[P[k]*cB+i];
      for (j=0; j<k; j++)
	t -= LU[k*r+j]*C[j*cB+i];
      C[k*cB+i] = t;
    }
  }

  // Solve UB = C
  for (i=0; i<cB; i++){ // loop over columns of B
    for (k=(r-1); k>=0; k--){   
      t = C[k*cB+i];
      for (j=k+1; j<r; j++)
	t -= LU[k*r+j]*B[j*cB+i];
      B[k*cB+i] = t/LU[k*r+k];
    }
  }

  xf_Release( (void *) C);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SolvePLU_MatrixR
static int
xf_SolvePLU_MatrixR(real *LU, int *P, int r, int rB, real *B)
{
  // B = B*A^{-1}, B is rB x r, A is r x r and PLU factored
  int ierr, i, j, k;
  real *C, *pB, *pC, t;

  ierr = xf_Error(xf_Alloc((void **) &C, r*rB, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // Set B = B*U^{-1}
  for (i=0; i<r; i++){ // loop over columns of B
    if (i > 0)
      for (k=0; k<rB; k++) {
	pB = B + k*r;
	for (j=0; j<i; j++)  // loop over B columns before i
	  pB[i] -= pB[j]*LU[j*r+i];
      }
    t = LU[i*r+i];
    for (k=0; k<rB; k++)
      B[k*r+i] /= t;

  } // i

  // Set C = B*L^{-1}
  for (i=r-1; i>=0; i--){
    for (k=0; k<rB; k++)
      C[k*r+i] = B[k*r+i];
    if (i < (r-1))
      for (k=0; k<rB; k++){
	pB = B + k*r;
	pC = C + k*r;
	for (j=i+1; j<r; j++){
	  pC[i] -= pC[j]*LU[j*r+i];
	}
      } // k
  } // i

  // Set B = C*P
  for (i=0; i<r; i++)
    for (j=0; j<rB; j++)
      B[j*r+P[i]] = C[j*r+i];
  
  xf_Release( (void *) C);
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ComputeBlockPLU
int
xf_ComputeBlockPLU(real *A, int n, int s, int *P)
{
  int ierr, i, j, k, s2, *pP;
  real *pA, *pL;

  s2 = s*s;

  // set P = 0
  for (i=0; i<n*s; i++) P[i] = 0;

  // test trivial case
  if ((n==1) && (s==1) && (fabs(A[0]) <= 0.0)) return xf_SINGULAR;

  // loop over row blocks
  for (i=0; i<n; i++){

    pA = A+(i*n+i)*s2;
    pP = P+i*s;

    // compute inplace PLU of on-diagonal block (A_{i,i})
    ierr = xf_Error(xf_ComputePLU(pA, s, pP));
    if (ierr != xf_OK) return ierr;

    // loop over row blocks below
    for (j=i+1; j<n; j++){
      // set L_{j,i} = A_{j,i}*A_{i,i}^{-1}
      pL = A + (j*n+i)*s2;
      ierr = xf_Error(xf_SolvePLU_MatrixR(pA, pP, s, s, pL));
      if (ierr != xf_OK) return ierr;

      // loop over column blocks
      for (k=i+1; k<n; k++){
	// set U{j,k} = A{j,k}-L_{j,i}*A{i,k}
	xf_MxM_Sub(pL, pA+(k-i)*s2, s, s, s, pL+(k-i)*s2);
      } // k
    } // j
  } // i

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SolveBlockPLU
int
xf_SolveBlockPLU(real *A, int n, int s, int *P, real *b, 
		 enum xfe_AddType AddFlag, real *u)
{
  int ierr, i, j, k, is, s2;
  real *y, *ty, *pL, *pU;

  ierr = xf_Error(xf_Alloc((void **) &y, n*s+s, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ty = y+n*s;

  if (s == 1){ // for speedup
    ierr = xf_Error(xf_SolvePLU(A, NULL, b, n, y, NULL));
    if (ierr != xf_OK) return ierr;
    xf_V_Add(y, n, AddFlag, u);
    xf_Release( (void *) y);
    return xf_OK;
  }


  s2 = s*s;

  // solve Ly = b (block form)
  for (i=0; i<n; i++){
    is = i*s;
    for (k=0; k<s; k++) y[is+k] = b[is+k];
    pL = A + i*n*s2;
    for (j=0; j<i; j++)
      xf_MxV_Sub(pL+j*s2, y+j*s, s, s, y+is);
  }

  // solve Uy = y
  for (i=n-1; i>=0; i--){
    is = i*s;
    pU = A + i*n*s2;
    for (j=i+1; j<n; j++)
      xf_MxV_Sub(pU+j*s2, y+j*s, s, s, y+is);
    ierr = xf_Error(xf_SolvePLU(pU+i*s2, P+is, y+is, s, y+is, ty));
    if (ierr != xf_OK) return ierr; 
  }

  // set/add/sub
  xf_V_Add(y, n*s, AddFlag, u);
  
  xf_Release( (void *) y);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SolveBlockPLUT
int
xf_SolveBlockPLUT(real *A, int n, int s, int *P, real *b, 
		  enum xfe_AddType AddFlag, real *u)
{
  int ierr, i, j, k, is, s2;
  real *y, *ty, *pL, *pU;

  ierr = xf_Error(xf_Alloc((void **) &y, n*s+s, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ty = y+n*s;

  if (s == 1){ // for speedup
    ierr = xf_Error(xf_SolvePLUT(A, NULL, b, n, y, NULL));
    if (ierr != xf_OK) return ierr;
    xf_V_Add(y, n, AddFlag, u);
    xf_Release( (void *) y);
    return xf_OK;
  }


  s2 = s*s;

  // solve U^T y = b (block form)
  for (i=0; i<n; i++){
    is = i*s;
    for (k=0; k<s; k++) y[is+k] = b[is+k];
    pU = A + i*s2;
    for (j=0; j<i; j++)
      xf_MTxV_Sub(pU+j*n*s2, y+j*s, s, s, y+is);
    ierr = xf_Error(xf_SolvePLUT(pU+i*n*s2, P+is, y+is, s, y+is, ty));
    if (ierr != xf_OK) return ierr; 
  }

  // solve L^T y = y
  for (i=n-1; i>=0; i--){
    is = i*s;
    pL = A + i*s2;
    for (j=i+1; j<n; j++)
      xf_MTxV_Sub(pL+j*n*s2, y+j*s, s, s, y+is);
  }

  // set/add/sub 
  xf_V_Add(y, n*s, AddFlag, u);

  xf_Release( (void *) y);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SolveBlockPLU_Matrix
int
xf_SolveBlockPLU_Matrix(real *A, int n, int s, int *P, real *B, 
			int m, enum xfe_AddType AddFlag, real *U)
{
  int ierr, i, j, k, is, s2, ims2;
  real *Y, *pL, *pU;

  s2 = s*s;

  ierr = xf_Error(xf_Alloc((void **) &Y, n*s2*m, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // solve LY = B (block form)
  for (i=0; i<n; i++){
    ims2 = i*m*s2;
    for (k=0; k<m*s2; k++) Y[ims2+k] = B[ims2+k];
    for (j=0; j<i; j++){
      pL = A + (i*n+j)*s2;
      for (k=0; k<m; k++)
	xf_MxM_Sub(pL, Y+(j*m+k)*s2, s, s, s, Y+(i*m+k)*s2);
    }
  }

  // solve UY = Y
  for (i=n-1; i>=0; i--){
    is = i*s;
    pU = A + i*n*s2;
    for (j=i+1; j<n; j++)
      for (k=0; k<m; k++)
	xf_MxM_Sub(pU+j*s2, Y+(j*m+k)*s2, s, s, s, Y+(i*m+k)*s2);
    for (k=0; k<m; k++){
      ierr = xf_Error(xf_SolvePLU_Matrix(pU+i*s2, P+is, s, s, Y+(i*m+k)*s2));
      if (ierr != xf_OK) return ierr; 
    }
  }

  // set/add/sub 
  xf_V_Add(Y, n*s2*m, AddFlag, U);

  xf_Release( (void *) Y);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PLUMxV_Set
static int 
xf_PLUMxV_Set(real *A, int *P, real *u, int r, real *b, real *yin)
{
  int ierr;
  int k, ik, j;
  real t, *y;

  if (yin != NULL)
    y = yin;
  else{
    ierr = xf_Error(xf_Alloc((void **) &y, r, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  /* y = U * u  */
  ik = 0;
  for(k=0; k<r; k++){
    t = 0.0;
    for(j=k; j<r; j++)
      t += A[ik+j] * u[j];
    y[k] = t;
    ik += r;
  }
  
  /* y = L * y  */
  ik = (r-1)*r;
  for(k=r-1; k>=0; k--){
    t = y[k];
    for(j=0; j<k; j++)
      t += A[ik+j] * y[j];
    y[k] = t;
    ik -= r;
  }
  
  /* b = inv(P) * y */
  for(k=0; k<r; k++)
    b[P[k]] = y[k];

  if (yin == NULL)
    xf_Release((void *)y);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PLUMTxV_Set
static int 
xf_PLUMTxV_Set(real *A, int *P, real *u, int r, real *b, real *yin)
{
  int ierr;
  int k, ik, j;
  real t, *y;

  if (yin != NULL)
    y = yin;
  else{
    ierr = xf_Error(xf_Alloc((void **) &y, r, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  /* y = P u */
  for(k=0; k<r; k++)
    y[k] = u[P[k]];

  /* y = L^T y */
  for(k=0; k<r; k++){
    t = y[k];
    for(j=k+1, ik=(k+1)*r; j<r; j++, ik+=r)
      t += A[ik+k] * y[j];
    y[k] = t;
  }

  /* b = U^T y */
  for(k=r-1; k>=0; k--){
    t = 0.0;
    for(j=0, ik=0; j<=k; j++, ik+=r)
      t += A[ik+k] * y[j];
    b[k] = t;
  }

  if (yin == NULL)
    xf_Release((void *) y);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BlockPLUMxV
int
xf_BlockPLUMxV(real *A, int n, int s, int *P, real *u, 
	       enum xfe_AddType AddFlag, real *b)
{
  int ierr, i, j, k, is, s2;
  real *y, *ty, *pL, *pU;

  if ((n==1) && (s==1)){
    switch(AddFlag){
    case xfe_Set:  b[0]  =  A[0]*u[0]; break;
    case xfe_Neg:  b[0]  = -A[0]*u[0]; break;
    case xfe_Add:  b[0] +=  A[0]*u[0]; break;
    case xfe_Sub:  b[0] -=  A[0]*u[0]; break;
    default: xf_Error(xf_INPUT_ERROR); break;
    }
    return xf_OK;
  }

  ierr = xf_Error(xf_Alloc((void **) &y, n*s+s, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ty = y + n*s;

  s2 = s*s;

  // y = U*u
  for (i=0; i<n; i++){
    is = i*s;
    pU = A + i*n*s2;
    ierr = xf_Error(xf_PLUMxV_Set(pU+i*s2, P+is, u+is, s, y+is, ty));
    if (ierr != xf_OK) return ierr;
    for (j=i+1; j<n; j++)
      xf_MxV_Add(pU+j*s2, u+j*s, s, s, y+is);
  }

  // y = L*y
  for (i=n-1; i>0; i--){
    is = i*s;
    pL = A + i*n*s2;
    for (j=0; j<i; j++)
      xf_MxV_Add(pL+j*s2, y+j*s, s, s, y+is);
  }
  
  // b @= y
  xf_V_Add(y, n*s, AddFlag, b);
 
  xf_Release( (void *) y);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BlockPLUMTxV
int
xf_BlockPLUMTxV(real *A, int n, int s, int *P, real *u, 
		enum xfe_AddType AddFlag, real *b)
{
  int ierr, i, j, k, is, s2;
  real *y, *ty, *pL, *pU;

  ierr = xf_Error(xf_Alloc((void **) &y, n*s+s, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ty = y + n*s;
  
  s2 = s*s;

  // y = L^T*u
  for (i=0; i<n; i++){
    is = i*s;
    pL = A + i*s2;
    for (k=0; k<s; k++) y[is+k] = u[is+k];
    for (j=i+1; j<n; j++)
      xf_MTxV_Add(pL+j*n*s2, u+j*s, s, s, y+is);
  }

  // y = U^T*y
  for (i=n-1; i>=0; i--){
    is = i*s;
    pU = A + i*s2;
    ierr = xf_Error(xf_PLUMTxV_Set(pU+i*n*s2, P+is, y+is, s, y+is, ty));
    if (ierr != xf_OK) return ierr;
    for (j=0; j<i; j++)
      xf_MTxV_Add(pU+j*n*s2, y+j*s, s, s, y+is);
  }
  
  // b @= y
  xf_V_Add(y, n*s, AddFlag, b);

  xf_Release( (void *) y);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BlockMxV
int
xf_BlockMxV(real *A, int n, int s, int m, real *u, 
	    enum xfe_AddType AddFlag, real *b)
{
  int i, j, k, l, ks, s2;
  real t, *pb, *pA, *pu, sg;

  if (s == 1){ // for speedup
    if ((n==1) && (m==1)){
      switch(AddFlag){
      case xfe_Set:  b[0]  =  A[0]*u[0]; break;
      case xfe_Neg:  b[0]  = -A[0]*u[0]; break;
      case xfe_Add:  b[0] +=  A[0]*u[0]; break;
      case xfe_Sub:  b[0] -=  A[0]*u[0]; break;
      default: xf_Error(xf_INPUT_ERROR); break;
      }
      return xf_OK;
    }
    xf_MxV(A, u, n, m, AddFlag, b);
    return xf_OK;
  }

  if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg))
    for (k=0; k<n*s; k++) b[k] = 0.0;

  switch(AddFlag){
  case xfe_Set: 
  case xfe_Add: 
    sg = 1.; break;
  case xfe_Neg: 
  case xfe_Sub: 
    sg = -1.0; break;
  default: return xf_Error(xf_OUT_OF_BOUNDS); break;
  }

  s2 = s*s;

  for (i=0; i<n; i++){
    pb = b + i*s;
    for (j=0; j<m; j++){
      pA = A + (i*m + j)*s2;
      pu = u + j*s;
      for (k=0; k<s; k++){
	ks = k*s;
	for (l=0, t=0.; l<s; l++) t += pA[ks+l]*pu[l];
	pb[k] += t*sg;
      } // k
    } // j
  } // i
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BlockMTxV
int
xf_BlockMTxV(real *A, int n, int s, int m, real *u, 
	     enum xfe_AddType AddFlag, real *b)
{
  // A^T is n x m blocks (each s x s)
  int i, j, k, l;
  real t, *pb, *pA, *pu;

  if (s == 1){ // for speedup
    if ((n==1) && (m==1)){
      switch(AddFlag){
      case xfe_Set:  b[0]  =  A[0]*u[0]; break;
      case xfe_Neg:  b[0]  = -A[0]*u[0]; break;
      case xfe_Add:  b[0] +=  A[0]*u[0]; break;
      case xfe_Sub:  b[0] -=  A[0]*u[0]; break;
      default: xf_Error(xf_INPUT_ERROR); break;
      }
      return xf_OK;
    }
    xf_MTxV(A, u, n, m, AddFlag, b);
    return xf_OK;
  }


  for (i=0; i<n; i++){
    pb = b + i*s;
    if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg))
      for (k=0; k<s; k++) pb[k] = 0.0;
    for (j=0; j<m; j++){
      pA = A + (j*n + i)*s*s;
      pu = u + j*s;
      for (k=0; k<s; k++){
	t = 0.0;
	for (l=0; l<s; l++) t += pA[l*s+k]*pu[l];
	switch(AddFlag){
	case xfe_Set: 
	case xfe_Add: 
	  pb[k] += t; break;
	case xfe_Neg: 
	case xfe_Sub: 
	  pb[k] -= t; break;
	default: return xf_Error(xf_OUT_OF_BOUNDS); break;
	}
      } // k
    } // j
  } // i
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BlockMxM
void
xf_BlockMxM(const real *A, int n, int s, int m, const real *U, 
	    int l, enum xfe_AddType AddFlag, real *B)
{
  // U is not a block matrix
  int in, s2;

  s2 = s*s;
  for (in=0; in<n; in++)
    xf_MTxM(U, A + in*m*s2, l, m, s2, AddFlag, B+in*l*s2);
}



/******************************************************************/
//   FUNCTION Definition: xf_MxBlockM
void
xf_MxBlockM(const real *U, int l, int n, const real *A, int s,
	    int m, enum xfe_AddType AddFlag, real *B)
{
  // U is not a block matrix
  int in, s2;

  s2 = s*s;
  
  xf_MxM(U, A, l, n, m*s2, AddFlag, B);
}


/******************************************************************/
//   FUNCTION Definition: xf_MTxBlockM
void
xf_MTxBlockM(const real *U, int l, int n, const real *A, int s,
	     int m, enum xfe_AddType AddFlag, real *B)
{
  // U is not a block matrix
  int in, s2;

  s2 = s*s;
  
  xf_MTxM(U, A, l, n, m*s2, AddFlag, B);
}



/******************************************************************/
//   FUNCTION Definition: xf_BlockMxBlockM
int
xf_BlockMxBlockM(const real *A, int n, int s, int m, const real *U, 
		 int l, enum xfe_AddType AddFlag, real *B)
{
  int in, im, il, k, s2;
  const real *pA, *pU;
  real *pB;
  enum xfe_AddType AddFlag2;

  s2 = s*s;
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  for (in=0; in<n; in++){
    pB = B + in*l*s2;
    if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg))
      for (k=0; k<l*s2; k++) pB[k] = 0.0;
    for (im=0; im<m; im++){
      pA = A + (in*m+im)*s2;
      for (il=0; il<l; il++){
	pU = U + (im*l+il)*s2;
	pB = B + (in*l+il)*s2;
	xf_MxM(pA, pU, s, s, s, AddFlag2, pB);
      } // il
    } // im
  } // in

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_BlockIdentity
void
xf_BlockIdentity(int n, int s, real * RSTRCT A)
{
  int in, k, s2, off;

  s2 = s*s;
  
  // zero A out
  for (k=0; k<n*n*s2; k++) A[k] = 0.;

  for (in=0; in<n; in++){
    off = in*(n+1)*s2;
    for (k=0; k<s2; k+=(s+1)) A[off+k] = 1.;
  } // in
}



/******************************************************************/
//   FUNCTION Definition: xf_OutProd_Add
static void
xf_OutProd_Add(const real *u, const real *v, int n, int m, real *A)
{
  int in, im, inm;

  for (in=0; in<n; in++)
    for (im=0, inm=in*m; im<m; im++)
      A[inm+im] += u[in]*v[im];
}


/******************************************************************/
//   FUNCTION Definition: xf_OutProd_Sub
static void
xf_OutProd_Sub(const real *u, const real *v, int n, int m, real *A)
{
  int in, im, inm;

  for (in=0; in<n; in++)
    for (im=0, inm=in*m; im<m; im++)
      A[inm+im] -= u[in]*v[im];
}

/******************************************************************/
//   FUNCTION Definition: xf_BlockOutProd_Add
void
xf_BlockOutProd_Add(const real *u, const real *v, int n, int s, 
		    int m, real *A)
{
  int in, im, s2;
  real *pA;

  s2 = s*s;

  for (in=0; in<n; in++){
    pA = A + in*m*s2;
    for (im=0; im<m; im++)
      xf_OutProd_Add(u+in*s, v+im*s, s, s, pA+im*s2);
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_BlockOutProd_Sub
void
xf_BlockOutProd_Sub(const real *u, const real *v, int n, int s, 
		    int m, real *A)
{
  int in, im, s2;
  real *pA;

  s2 = s*s;

  for (in=0; in<n; in++){
    pA = A + in*m*s2;
    for (im=0; im<m; im++)
      xf_OutProd_Sub(u+in*s, v+im*s, s, s, pA+im*s2);
  }
}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateHouseholder
static int 
xf_CalculateHouseholder(real *R, int rA, int cA, int ir, real *H)
{
  int ierr, i, j, len;
  real *u, *v, sign, norm;

  // initialize H to be the identity
  for (i=0; i<rA*rA; i++) H[i] = 0.;
  for (i=0; i<rA*rA; i+=(rA+1)) H[i] = 1.;

  if ((len = rA-ir) == 0) return xf_OK;

  // allocate memory
  ierr = xf_Error(xf_Alloc( (void **) &u, len, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  v = R + ir*cA + ir;
  
  sign = ((v[0] < 0.) ? -1.0 : 1.0);
  
  for (i=0, norm=0.; i<len; i++) norm += v[i*cA]*v[i*cA];
  norm = sqrt(norm);

  for (i=0; i<len; i++) u[i] = v[i*cA];
  u[0] += sign*norm;

  for (i=0, norm=0.; i<len; i++) norm += u[i]*u[i];

  for (i=ir; i<rA; i++)
    for (j=ir; j<rA; j++)
      H[i*rA+j] -= 2*u[i-ir]*u[j-ir]/norm;

  xf_Release( (void *) u);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_QRFactorHouseholder
int 
xf_QRFactorHouseholder(real *A, int rA, int cA, real *Q, real *R)
{
  int ierr, i, j; 
  real *H, *T;

  if (cA > rA) return xf_Error(xf_INPUT_ERROR);

  // allocate matrices
  ierr = xf_Error(xf_Alloc( (void **) &T, rA*rA, sizeof(real)));
  if (ierr != xf_OK) return ierr;
   
  ierr = xf_Error(xf_Alloc( (void **) &H, rA*rA, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // R = A 
  for (i=0; i<rA*cA; i++) R[i] = A[i]; 

  for(i=0; i<cA; i++){ // loop over columns

    // calculate householder matrix H using column of i'th minor of R
    ierr = xf_Error(xf_CalculateHouseholder(R, rA, cA, i, H));
    if (ierr != xf_OK) return ierr;

    // R = H * R
    for (j=0; j<rA*cA; j++) T[j] = R[j];
    xf_MxM_Set(H, T, rA, rA, cA, R);

    // Q = Q * H
    if(i == 0)
      for(j=0; j<rA*rA; j++) Q[j] = H[j];
    else{
      for (j=0; j<rA*rA; j++) T[j] = Q[j];
      xf_MxM_Set(T, H, rA, rA, rA, Q);
    }
  }
  
  xf_Release( (void *) H);
  xf_Release( (void *) T);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_EigSymmetricQR
int 
xf_EigSymmetricQR(real *A, int r, real tol, int maxiter, 
		  real *EG, real *EV)
{
  int ierr, iter, i, j, jmax;
  int *eflag;
  enum xfe_Bool converged;
  real emax;
  real *B, *S, *Stemp;
  real *Q, *R;


  // allocate temporary matrices
  ierr = xf_Error(xf_Alloc((void **) &S, r*r, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) &Stemp, r*r, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) &Q, r*r, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) &R, r*r, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // allocate a vector B for the iteration
  ierr = xf_Error(xf_Alloc((void **) &B, r*r, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<r*r; i++) B[i] = A[i]; // B = A

  iter = 0;
  converged = xfe_False;
  while(!converged){

    // factor B = QR
    ierr = xf_Error(xf_QRFactorHouseholder(B, r, r, Q, R));
    if (ierr != xf_OK) return ierr;
    
    if (iter == 0)
      for (i=0; i<r*r; i++) S[i] = Q[i];
    else      
      {
	// S = S*Q
	for (i=0; i<r*r; i++) Stemp[i] = S[i];
	xf_MxM_Set(Stemp, Q, r, r, r, S);
      }
    
    // B = R * Q
    xf_MxM_Set(R, Q, r, r, r, B);
    
    // check tolerances on D
    converged = xfe_True;
    for(i=0; i<r; i++){
      for(j=i+1; j<r; j++)
	if(fabs(B[i*r + j]) > tol){
	  converged = xfe_False;
	  break;
	}
      if (!converged) break;
    }// end for all below diagonal

    if(iter++ > maxiter) break;   
  }// end while    

  // for subsequent insert sort
  ierr = xf_Error(xf_Alloc((void **) &eflag, r, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<r; i++) eflag[i] = 0;

  // sort eigs + vectors via insert sort
  for (i=0; i<r; i++){
    emax = -1.0;
    jmax = -1;
    for (j=0; j<r; j++){
      if (eflag[j]) continue;
      if (fabs(B[j*r+j]) > emax){
	emax = fabs(B[j*r+j]);
	jmax = j;
      }
    } // j
    if (jmax < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
    EG[i] = B[jmax*r+jmax];
    for (j=0; j<r; j++)
      EV[i*r+j] = S[j*r+jmax]; // note, storing eigvs in rows of EV
    eflag[jmax] = 1;
  } // i

  xf_Release( (void *) S);
  xf_Release( (void *) Stemp);
  xf_Release( (void *) Q);
  xf_Release( (void *) R);
  xf_Release( (void *) B);
  xf_Release( (void *) eflag);
    
  if (iter > maxiter) return xf_NOT_CONVERGED;

  return xf_OK;

}

/******************************************************************/
//   FUNCTION Definition: xf_SVDGolubReinsch
int 
xf_SVDGolubReinsch(real *U, int m, int n, enum xfe_Bool NeedU, real *W, real *V)
{
  int ierr, i, j, k, l;
  int i1, l1, k1, its, its_max = 1000;
  enum xfe_Bool needsplit, converged;
  real *E;
  real f, g, h, scale, x, y, z, s, c, temp, eps;

  // check for input error
  if (m < n) return xf_Error(xf_INPUT_ERROR);

  // allocate temporary storage
  ierr = xf_Error(xf_Alloc((void **) &E, n, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  /*------------------------------------------*/
  /* Householder reduction to bidiagonal form */
  /*------------------------------------------*/
  g = 0.; scale = 0.; x = 0.;
  for (i=0; i<n; i++){
    l = i+1;
    E[i] = scale*g;
    g = 0.; s = 0.;
    
    /* left transformations that zero out subdiag elems of ith col */
    for (k=i, scale=0.; k<m; k++) scale += fabs(U[k*n+i]);
    if (scale != 0.){
      for (k=i; k<m; k++){
	temp = U[k*n+i] = U[k*n+i]/scale;
	s += temp*temp;
      }
      f = U[i*n+i];
      g = ((f < 0) ? sqrt(s) : -sqrt(s));
      h = f*g - s;
      U[i*n+i] = f - g;
      // apply left transformations to remaining columns of A
      for (j=l; j<n; j++){
	for (k=i, s=0.; k<m; k++) s += U[k*n+i]*U[k*n+j];
	f = s/h;
	for (k=i; k<m; k++) U[k*n+j] += f*U[k*n+i];
      } // j
      
      for (k=i; k<m; k++) U[k*n+i] *= scale;
    }
    W[i] = scale*g;

    /* compute the right transformations */
    g = 0.; s = 0.;
    for (k=l, scale=0.; k<n; k++) scale += fabs(U[i*n+k]);
    if ((l < n) && (scale != 0.)){
      for (k=l; k<n; k++){
	temp = U[i*n+k] = U[i*n+k]/scale;
	s += temp*temp;
      }
      f = U[i*n+l];
      g = ((f < 0) ? sqrt(s) : -sqrt(s));
      h = f*g - s;
      U[i*n+l] = f - g;
      for (k=l; k<n; k++) E[k] = U[i*n+k]/h;
      for (j=l; j<m; j++){
	for (k=l, s=0.; k<n; k++) s += U[j*n+k]*U[i*n+k];
	for (k=l; k<n; k++) U[j*n+k] += s*E[k];
      } // j
      for (k=l; k<n; k++) U[i*n+k] *= scale;
    }
    
    x = max(x, fabs(W[i]) + fabs(E[i]));

  } // i

  /* Accumulation of right-hand transformations into V */
  if (V != NULL){
    for (i=n-1; i>=0; i--){
      l = i+1;
      if ((l<n) && (g != 0.)){
	for (j=l; j<n; j++) V[j*n+i] = (U[i*n+j]/U[i*n+l])/g;
	for (j=l; j<n; j++){
	  for (k=l, s=0.; k<n; k++) s += U[i*n+k]*V[k*n+j];
	  for (k=l; k<n; k++) V[k*n+j] += s*V[k*n+i];
	}
      }
      for (j=l; j<n; j++) V[i*n+j] = V[j*n+i] = 0.;
      V[i*n+i] = 1.;
      g = E[i];
    } // i
  }

  /* Accumulation of left-hand transformations in U */
  if (NeedU){
    for (i=n-1; i>=0; i--){
      l = i+1;
      g = W[i];
      for (j=l; j<n; j++) U[i*n+j] = 0.;
      if (g != 0.){
	for (j=l; j<n; j++){
	  for (k=l, s=0.; k<m; k++) s += U[k*n+i]*U[k*n+j];
	  f = (s/U[i*n+i])/g;
	  for (k=i; k<m; k++) U[k*n+j] += f*U[k*n+i];
	}
	for (j=i; j<m; j++) U[j*n+i] /= g;
      }
      else
	for (j=i; j<m; j++) U[j*n+i] = 0.;
      U[i*n+i] += 1.0;
    } // i
  } // NeedU

  

  /*------------------------------------*/
  /* Diagonalization of bidiagonal form */
  /*------------------------------------*/

  eps = MEPS*x;
  for (k=n-1; k>=0; k--){
    k1 = k-1;
    converged = xfe_False;
    its = 0;
    while ((!converged) && (its < its_max)){
      
      // test for splitting
      needsplit = xfe_False;
      for (l=k; l>=0; l--){
	if (fabs(E[l]) <= eps) break;
	// E[0] is always 0, so should never get here with l==0
	if ((l1 = l-1) < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
	if (fabs(W[l1]) <= eps){
	  needsplit = xfe_True;
	  break;
	}
      }
      
      if (needsplit){
	// cancellation
	c = 0.; s = 1.;
	for (i=l; i<k; i++){
	  f = s*E[i]; E[i] = c*E[i];
	  if (fabs(f) <= eps) break;
	  g = W[i]; h = sqrt(f*f+g*g); c = g/h; s = -f/h;

	  if (NeedU){
	    for (j=0; j<m; j++){
	      y = U[j*n+l1]; z = U[j*n+i];
	      U[j*n+l1] =  y*c + z*s;
	      U[j*n+i ] = -y*s + z*c;
	    } // j
	  }
	} // i
      }
      
      z = W[k];
      if (l == k){// test for convergence
	converged = xfe_True;
	break; // out of while loop
      }
      // shift from bottom 2 x 2 minor
      if (k1 < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
      x = W[l]; y = W[k1]; g = E[k1]; h = E[k];
      f = ((y-z)*(y+z) + (g-h)*(g+h))/(2.*h*y);
      g = sqrt(f*f+1.);
      temp = ((f<0) ? f-g : f+g);
      f = ((x-z)*(x+z) + h*(y/temp-h))/x;
      
      // next QR transformation
      c = 1.; s = 1.;
      for (i1=l; i1<k; i1++){
	i = i1+1;
	g = E[i]; y = W[i]; h = s*g; g = c*g; z = sqrt(f*f+h*h); 
	E[i1] = z; c = f/z; s = h/z;
	f = x*c + g*s; g = -x*s + g*c; h = y*s; y = y*c;
	if (V != NULL){
	  for (j=0; j<n; j++){
	    x = V[j*n+i1]; z = V[j*n+i];
	    V[j*n+i1] = x*c + z*s;
	    V[j*n+i] = -x*s + z*c;
	  } // j
	}
	z = sqrt(f*f+h*h);
	W[i1] = z;
	if (z != 0.){
	  c = f/z; 
	  s = h/z;
	}
	f =  c*g + s*y;
	x = -s*g + c*y;

	if (NeedU){
	  for (j=0; j<m; j++){
	    y = U[j*n+i1];
	    z = U[j*n+i ];
	    U[j*n+i1] =  y*c + z*s;
	    U[j*n+i ] = -y*s + z*c;
	  } // j
	}

      } // i1
      E[l] = 0.;
      E[k] = f;
      W[k] = x;
    
      its++;
    } // while !converged

    if (!converged){
      xf_printf("GR SVD: inner loop not converged after %d iterations.\n", its);
      return xf_NOT_CONVERGED;
    }

    if (z < 0){ 
      // make W[k] non-negative
      W[k] = -W[k];
      if (V != NULL)
	for (j=0; j<n; j++) V[j*n+k] = -V[j*n+k];
    }

  } // k

  xf_Release( (void *) E);

  return xf_OK;

}


/******************************************************************/
//   FUNCTION Definition: xf_GivensRotation
static void 
xf_GivensRotation(real f, real g, real *c, real *s, real *r)
{
  real rr;
  // Calculates a Givens rotation from f and g

  if (g == 0.0){
    (*c) = 1.0;
    (*s) = 0.0;
    (*r) = f;
  }
  else if (f == 0.0){
    (*c) = 0.0;
    (*s) = 1.0;
    (*r) = g;
  }
  else{
    rr = sqrt(f*f + g*g);
    (*c) = f/rr;
    (*s) = g/rr;
    (*r) = rr;
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_QRBulgeChase
int 
xf_QRBulgeChase(real *shift, int nk, int np, real *alpha,
		real *beta, real *Q)
{
  int nm, i, j, k, jj;
  int itop, istart, iend;
  enum xfe_Bool done_blocks;
  real big, f, g, a1, a2, a3, a4, c, s, r;

  // check for input error
  if (np == 0) return xf_OK;
  if ((nk <= 0) || (np < 0)) return xf_Error(xf_INPUT_ERROR);

  nm = np+nk;

  // set Q = identity (Q must be already allocated)
  for (i=0; i<nm*nm; i++) Q[i] = 0.0;
  for (i=0; i<nm*nm; i+=nm+1) Q[i] = 1.0;

  itop = 0;

  // loop over shifts
  for (jj=0; jj<np; jj++){
    istart = itop;
    
    done_blocks = xfe_False;
    while (!done_blocks){

      // set iend based on non-negligible off-diagonal elements
      iend = nm-1;
      for (i=istart; i<nm-1; i++){
	big = fabs(alpha[i]) + fabs(alpha[i+1]);
	if (beta[i+1] <= MEPS*big){
	  beta[i+1] = 0.0;
	  iend = i;
	}
      } // i

      if (istart < iend) {
	
	// Construct plane rotation to drive h(istart+1,1) to zero
	f = alpha[istart] - shift[jj];
	g = beta[istart+1];
	xf_GivensRotation(f, g, &c, &s, &r);
	
	// Apply rotation to left and right of H
	a1 = c*alpha[istart]   + s* beta[istart+1];
	a2 = c* beta[istart+1] + s*alpha[istart+1];
	a4 = c*alpha[istart+1] - s* beta[istart+1];
	a3 = c* beta[istart+1] - s*alpha[istart];
	alpha[istart]   = c*a1 + s*a2;
	alpha[istart+1] = c*a4 - s*a3;
	beta[istart+1]  = c*a3 + s*a4;
    
	// Accumulate rotation into Q
	for (j=0; j<=min(istart+jj+1, nm-1); j++){
	  a1               =  c*Q[j*nm+istart] + s*Q[j*nm+istart+1];
	  Q[j*nm+istart+1] = -s*Q[j*nm+istart] + c*Q[j*nm+istart+1];
	  Q[j*nm+istart] = a1;
	}

	// Now chase the bulge
	for (i=istart+1; i<iend; i++){

	  // Construct plane rotation to zero out bulge created by previous rotation
	  f = beta[i];
	  g = s*beta[i+1];

	  // final update with previous rotation
	  beta[i+1] *= c;
      
	  xf_GivensRotation(f, g, &c, &s, &r);
      
	  // ensure off-diagonal elements of H remain nonnegative
	  if (r < 0){
	    r = -r; c = -c; s = -s;
	  }

	  // Apply rotation to left and right of H
	  beta[i] = r;

	  a1 = c*alpha[i]   + s* beta[i+1];
	  a2 = c* beta[i+1] + s*alpha[i+1];
	  a3 = c* beta[i+1] - s*alpha[i];
	  a4 = c*alpha[i+1] - s* beta[i+1];
	  alpha[i]   = c*a1 + s*a2;
	  alpha[i+1] = c*a4 - s*a3;
	  beta[i+1]  = c*a3 + s*a4;

	  // Accumulate rotation into Q
	  for (j=0; j<=min(i+jj+1, nm-1); j++){
	    a1          =  c*Q[j*nm+i] + s*Q[j*nm+i+1];
	    Q[j*nm+i+1] = -s*Q[j*nm+i] + c*Q[j*nm+i+1];
	    Q[j*nm+i] = a1;
	  }
	} // i
      } // end if (istart < iend)
    
      istart = iend + 1;

      // Ensure beta(iend) is nonnegative
      if (beta[iend] < 0.0){
	beta[iend] = -beta[iend];
	for (k=0; k<nm; k++) Q[k*nm+iend] *= -1.0;
      }

      // Apply the same shift to the next block if there is any
      done_blocks = (iend >= nm-1);
    }
    
    // Increase the start block if possible
    for (i=itop; i<nm; i++){
      if (beta[i+1] > 0.0) break;
      itop++;
    }

    // Finished applying jj-th shift

  } // j

  // Check for more possible deflation after the last shift
  for (i=itop; i<nm-1; i++){
    big = fabs(alpha[i]) + fabs(alpha[i+1]);
    if (beta[i+1] <= MEPS*big) beta[i+1] = 0.0;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SortIntPos
int 
xf_SortIntPos(int *u, int n, int *pos, enum xfe_Bool InverseFlag)
{
  int ierr, beg[300], end[300], L, R, s, i=0, k, ipiv;
  int *p;
  int piv;

  ierr = xf_Error(xf_Alloc( (void **) &p, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  for (k=0; k<n; k++) p[k] = k; // initial positions

  beg[0]=0; end[0]=n;
  while (i>=0) {
    if (i >= 300) return xf_Error(xf_OUT_OF_BOUNDS);
    L = beg[i]; R = end[i]-1;
    if (L<R) {
      piv=u[L]; ipiv = p[L];
      while (L<R) {
        while ((u[R]>=piv) && (L<R)) R--; 
	if (L<R){
	  u[L]=u[R]; p[L] = p[R];
	  L++;
	}
	while ((u[L]<=piv) && (L<R)) L++; 
	if (L<R){
	  u[R]=u[L]; p[R] = p[L];
	  R--;
	}
      }
      u[L]=piv; p[L] = ipiv;
      beg[i+1]=L+1; 
      end[i+1]=end[i]; 
      end[i++]=L;
      if ((end[i]-beg[i]) > (end[i-1]-beg[i-1])) {
        s=beg[i]; beg[i]=beg[i-1]; beg[i-1]=s;
        s=end[i]; end[i]=end[i-1]; end[i-1]=s; 
      }
    }
    else {
      i--; 
    }
  }


  if (InverseFlag)   // pos is inverse of the mapping we calculated
    for (k=0; k<n; k++) pos[p[k]] = k; 
  else               // pos is the mapping we calculated
    for (k=0; k<n; k++) pos[k] = p[k]; 

  xf_Release( (void *) p);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SortReal
static int 
xf_SortReal(real *u, int n, int *pos, enum xfe_Bool InverseFlag)
{
  /* 
     PURPOSE

     Sorts n entries in u in ascending order and returns ranks in pos.

     Function modified from public-domain C implementation written by
     Darel Rex Finley

     INPUTS:

     u : vector of real entries to be sorted
     n : number of entries to be sorted
     InverseFlag : whether forward or inverse mapping is desired (see below)
  
     OUTPUTS:


     pos : must be pre-allocated.  pos[i] stores the post-sorting rank of
     the ith (original) entry in u if InverseFlag == True.
     Otherwise, pos[i] stores the index of the rank i entry.

     RETURNS:  Error code

  */
  int ierr, beg[300], end[300], L, R, s, i=0, k, ipiv;
  int *p;
  real piv;

  ierr = xf_Error(xf_Alloc( (void **) &p, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  for (k=0; k<n; k++) p[k] = k; // initial positions

  beg[0]=0; end[0]=n;
  while (i>=0) {
    if (i >= 300) return xf_Error(xf_OUT_OF_BOUNDS);
    L = beg[i]; R = end[i]-1;
    if (L<R) {
      piv=u[L]; ipiv = p[L];
      while (L<R) {
        while ((u[R]>=piv) && (L<R)) R--; 
	if (L<R){
	  u[L]=u[R]; p[L] = p[R];
	  L++;
	}
	while ((u[L]<=piv) && (L<R)) L++; 
	if (L<R){
	  u[R]=u[L]; p[R] = p[L];
	  R--;
	}
      }
      u[L]=piv; p[L] = ipiv;
      beg[i+1]=L+1; 
      end[i+1]=end[i]; 
      end[i++]=L;
      if ((end[i]-beg[i]) > (end[i-1]-beg[i-1])) {
        s=beg[i]; beg[i]=beg[i-1]; beg[i-1]=s;
        s=end[i]; end[i]=end[i-1]; end[i-1]=s; 
      }
    }
    else {
      i--; 
    }
  }


  if (InverseFlag)   // pos is inverse of the mapping we calculated
    for (k=0; k<n; k++) pos[p[k]] = k; 
  else               // pos is the mapping we calculated
    for (k=0; k<n; k++) pos[k] = p[k]; 

  xf_Release( (void *) p);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_MergeSort
static int 
xf_MergeSort(int nlist, real **uList, int *N, int **iList)
{
  /* 
     PURPOSE

     Sorts locally-sorted data in nlist lists, stored in uList[ilist],
     ilist<nlist, and returns global rank of each data element via the
     integer array iList[ilist].

     Uses an ascending sort.

     INPUTS:

     nlist : number of lists to merge
     uList : N[ilist] real data values in each list
     N : N[ilist] = number of data values in ilist
  
     OUTPUTS:

     iList : N[ilist] global ranks

     RETURNS:  Error code

  */
  int ierr;
  int i, itot, ntot;
  int imin;
  int *count;
  real rmin, uval;
  
  // ntot = total number of elements
  for (i=0, ntot=0; i<nlist; i++) ntot += N[i];
  if (ntot <= 0) return xf_Error(xf_INPUT_ERROR);

  // count[i] will keep track of how far we are in each list
  ierr = xf_Error(xf_Alloc( (void **) &count, nlist, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nlist; i++) count[i] = 0;


  // loop over total number, and pull off next largest data one by one
  for (itot=0; itot<ntot; itot++){
    
    imin = -1;
    rmin = 0.0;
    for (i=0; i<nlist; i++)
      if (count[i] < N[i]){
	uval = uList[i][count[i]];
	if ((imin < 0) || (uval < rmin)){
	  imin = i;
	  rmin = uval;
	}
      }
    if (imin < 0) return xf_Error(xf_INPUT_ERROR);

    iList[imin][count[imin]++] = itot;
  } // itot
  
  xf_Release( (void *) count);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SortRealParallel
int 
xf_SortRealParallel(real *u, int n, enum xfe_Bool SerialFlag, int *pos)
{
  int ierr, k;
  int myRank, nProc, iProc;
  int *N = NULL;
  int *loc2elem=NULL;
  int *loc2glob=NULL;
  int *ibuf = NULL;
  int **iList = NULL;
  real *rbuf = NULL;
  real **uList = NULL;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if ((nProc == 1) || (SerialFlag)){ // serial case is easy
    ierr = xf_Error(xf_SortReal(u, n, pos, xfe_True));
    if (ierr != xf_OK) return ierr;
    return xf_OK;
  }

  // allocate temporary arrays
  ierr = xf_Error(xf_Alloc( (void **) &loc2elem, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &loc2glob, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // sort local vectors: loc2elem[i] = index of locally ith smallest elem
  ierr = xf_Error(xf_SortReal(u, n, loc2elem, xfe_False));
  if (ierr != xf_OK) return ierr;

  // alloc memory specific to proc 0
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &N, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }

  // gather number of elems to inform proc 0 of # on each proc
  ierr = xf_Error(xf_MPI_Gather((void *) &n, (void *) N, 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  // proc 0 allocates uList=values and iList=ranks
  if (myRank == 0){
    ierr = xf_Error(xf_VAlloc2( (void ***) &uList, nProc, N, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_VAlloc2( (void ***) &iList, nProc, N, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }

  // variable gather u data onto proc 0
  rbuf = ((myRank == 0) ? (void *) uList[0] : NULL);
  ierr = xf_Error(xf_MPI_Gatherv((void *) u, n, rbuf, N, sizeof(real), 0));
  if (ierr != xf_OK) return ierr;

  // proc 0: merge-sort data in uList, store ranks in iList
  if (myRank == 0){
    ierr = xf_Error(xf_MergeSort(nProc, uList, N, iList));
    if (ierr != xf_OK) return ierr;
  }

  // variable scatter global ranks to each proc
  ibuf = ((myRank == 0) ? (void *) iList[0] : NULL);
  ierr = xf_Error(xf_MPI_Scatterv((void *) ibuf, N, (void *) loc2glob, n, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  // set pos[i] = glob rank of elem i
  for (k=0; k<n; k++) pos[loc2elem[k]] = loc2glob[k];

  xf_Release ( (void  *) loc2elem);
  xf_Release ( (void  *) loc2glob);
  xf_Release ( (void  *) N);
  xf_Release2( (void **) uList);
  xf_Release2( (void **) iList);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_QuadraticRoots
int 
xf_QuadraticRoots(real *coeff, int *nroots, real *roots)
{
  real a, b, c;
  real disc, sqdisc;

  // pull off coefficients
  a = coeff[0];
  b = coeff[1];
  c = coeff[2];

  // discriminant
  disc = b*b-4.*a*c;
  
  if (disc > 0.){
    // two real roots
    (*nroots) = 2;
    sqdisc = sqrt(disc);
    roots[0] = (-b - sqdisc)/(2.*a);
    roots[1] = (-b + sqdisc)/(2.*a);
  }
  else if (disc == 0.){
    // one real root
    (*nroots) = 1;
    roots[0] = -b/(2.*a);
  }
  else{
    // no real roots
    (*nroots) = 0.;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_IntComp
static int 
xf_IntComp(const void *a, const void *b)
{
  const int *ia = (const int *)a; // casting pointer types
  const int *ib = (const int *)b;
  return *ia  - *ib; 
  /* integer comparison: returns negative if b > a 
   and positive if a > b */
}


/******************************************************************/
//   FUNCTION Definition: xf_SortInt
int
xf_SortInt(int n, int *a)
{
  qsort(a, n, sizeof(int), xf_IntComp);	
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RandUniform
int 
xf_RandUniform(int n, real *pval)
{
  int ierr, myRank;
  int i;

  // check parallel
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  // root calculates random values
  if (myRank == 0){
    for (i=0; i<n; i++){
      pval[i] = ((real) rand())/((real) RAND_MAX);
    } // i
  }

  // broadcast to other procs
  ierr = xf_Error(xf_MPI_Bcast((void *) pval, n*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RandNormal
int 
xf_RandNormal(int n, real *pval)
{
  int ierr, myRank;
  int i;
  real x0, x1;

  // check parallel
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  // root calculates random values
  if (myRank == 0){
    for (i=0; i<n; i+=2){
      // two uniform random values
      x0 = ((real) rand())/((real) RAND_MAX);
      x1 = ((real) rand())/((real) RAND_MAX);
      // use Box-Muller transform
      pval[i] = sqrt(-2.0*log(x0))*cos(2.0*PI*x1);
      if ((i+1)< n)
	pval[i+1] = sqrt(-2.0*log(x0))*sin(2.0*PI*x1);
    }
  }

  // broadcast to other procs
  ierr = xf_Error(xf_MPI_Bcast((void *) pval, n*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BinSearch
int xf_BinSearch(const int target, const int *set, 
                 int begin, int end, int *rank)
{
  int center;
  int size;
  
  if (end < begin)
    return xf_INPUT_ERROR;
  
  size = end-begin+1;
  if (size == 1)
    if (target == set[begin]){
      if (rank != NULL)
        (*rank) = begin;
      return xf_OK;
    }
    else{
      if (rank != NULL){
        if (target < set[begin])
          (*rank) = begin-1;
        else
          (*rank) = end+1;
      }
      return xf_NOT_FOUND;
    }
  
  if (target == set[begin]){
    if (rank != NULL)
      (*rank) = begin;
    return xf_OK;
  }
  
  if (target == set[end]){
    if (rank != NULL)
      (*rank) = end;
    return xf_OK;
  }
  
  //if within range
  if (target > set[end]){
    if (rank != NULL)
      (*rank) = end+1;
    return xf_NOT_FOUND;
  }
  else if (target < set[begin]){
    if (rank != NULL)
      (*rank) = begin-1;
    return xf_NOT_FOUND;
  }
  else{
    while (size > 1){
      center = begin + size/2-1;
      if (target > set[center]){
        //in top portion of set
        begin = center+1;
      }
      else {
        end = center;
      }
      size = end-begin+1;
    }
    if (target == set[begin]){
      if (rank != NULL)
        (*rank) = begin;
      return xf_OK;
    }
    else{
      if (target > set[begin]){
        if (rank != NULL)
          (*rank) = begin+1;
      }
      else{
        if (rank != NULL)
          (*rank) = begin;
      }
      return xf_NOT_FOUND;
    }
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_RealNorm
real 
xf_RealNorm(const real * RSTRCT A, int r)
{
  int i;
  real val;
  for (i=0,val=0.; i<r; i++) val += A[i]*A[i];
  return sqrt(val);
}

/******************************************************************/
//   FUNCTION Definition: xf_Add2OrderedSet
int xf_Add2OrderedSet(const int entry, int *set_size, int **set,
                      int **orig_rank)
{
  int ierr, rank, movesize, src, dest;
  
  rank = 0;
  if ((*set_size) == 0)
    ierr = xf_NOT_FOUND;
  else 
    ierr = xf_BinSearch(entry, (*set), 0, (*set_size)-1, &rank);
  
  if (ierr == xf_NOT_FOUND || ierr == xf_OK){
    ierr = xf_Error(xf_ReAlloc((void**)&(*set), (*set_size)+1, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    if (orig_rank != NULL){
      ierr = xf_Error(xf_ReAlloc((void**)&(*orig_rank), (*set_size)+1, 
                                 sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    //make sure to keep the crescent order
    if (rank == (*set_size)){
      movesize = 0;
    }
    else if (rank == -1) {
      //move all one position forward
      rank = 0;
      src = 0;
      dest = 1;
      movesize = (*set_size);
    }
    else {
      src = rank;
      dest = rank+1;
      movesize = (*set_size)-rank;
    }
    if (movesize > 0) {
      if (memmove((*set)+dest,(*set)+src,movesize*sizeof(int)) == NULL)
        return xf_Error(xf_MEMORY_ERROR);
      if (orig_rank != NULL)
        if (memmove((*orig_rank)+dest,(*orig_rank)+src,
                    movesize*sizeof(int)) == NULL)
          return xf_Error(xf_MEMORY_ERROR);
    }
    if (orig_rank != NULL)
      orig_rank[0][rank] = (*set_size);
    set[0][rank] = entry;
    (*set_size)++;
  }
  else if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

#if( UNIT_TEST==1 )
#include "xf_Math.test.in"
#endif
