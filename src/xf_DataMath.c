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
  FILE:  xf_DataMath.c

  This file contains math functions operating on Data structures.

*/
#include <stdlib.h>

#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Data.h"
#include "xf_Basis.h"
#include "xf_MeshTools.h"
#include "xf_Math.h"
#include "xf_MathLapack.h"
#include "xf_Solver.h"
#include "xf_Line.h"


/******************************************************************/
//   FUNCTION Definition: xf_SetVector
int 
xf_SetVector( const xf_Vector *B, enum xfe_AddType AddFlag, xf_Vector *A)
{
  // A @= B, @ = AddFlag
  int i;
  xf_long k, ntot;
  real *rvalA, *rvalB;
  xf_GenArray *gA, *gB;

  if (!xf_CompatibleVectors(A,B)) return xf_Error(xf_INCOMPATIBLE);

  for (i = 0; i<A->nArray; i++){
    gA = A->GenArray+i;
    gB = B->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);
    
    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rvalA = gA->rValue[0];
    rvalB = gB->rValue[0];
    switch (AddFlag){
    case xfe_Set:
      for (k=0; k<ntot; k++) rvalA[k]  =  rvalB[k];
      break;
    case xfe_Neg:
      for (k=0; k<ntot; k++) rvalA[k]  = -rvalB[k];
      break;
    case xfe_Add:
      for (k=0; k<ntot; k++) rvalA[k] +=  rvalB[k];
      break;
    case xfe_Sub:
      for (k=0; k<ntot; k++) rvalA[k] -=  rvalB[k];
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    }
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetZeroVector
int 
xf_SetZeroVector( xf_Vector *A){
  // A = 0
  int i;
  xf_long k, ntot;
  int *ival;
  real *rval;
  xf_GenArray *gA;

  for (i = 0; i<A->nArray; i++){
    gA = A->GenArray+i;
    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    if (gA->Size == xfe_SizeInt){
      ival = gA->iValue[0];
      for (k=0; k<ntot; k++) ival[k] = 0;
    }
    else if (gA->Size == xfe_SizeReal){
      rval = gA->rValue[0];
      for (k=0; k<ntot; k++) rval[k] = 0.0;
    }
    else return xf_Error(xf_NOT_SUPPORTED);	
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SetZeroHaloVector
int 
xf_SetZeroHaloVector( xf_Vector *A){
  // A = 0 on halo elements only
  int i;
  xf_long k, ntot;
  int *ival;
  real *rval;
  xf_GenArray *gA;

  for (i = A->nArraySelf; i<A->nArray; i++){
    gA = A->GenArray+i;
    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    if (gA->Size == xfe_SizeInt){
      ival = gA->iValue[0];
      for (k=0; k<ntot; k++) ival[k] = 0;
    }
    else if (gA->Size == xfe_SizeReal){
      rval = gA->rValue[0];
      for (k=0; k<ntot; k++) rval[k] = 0.0;
    }
    else return xf_Error(xf_NOT_SUPPORTED);	
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetConstVector
int 
xf_SetConstVector( xf_Vector *A, int iv, real rv){
  // A = ival (for integer) or rval (for real)
  int i;
  xf_long k, ntot;
  int *ival;
  real *rval;
  xf_GenArray *gA;

  for (i = 0; i<A->nArray; i++){
    gA = A->GenArray+i;
    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    if (gA->Size == xfe_SizeInt){
      ival = gA->iValue[0];
      for (k=0; k<ntot; k++) ival[k] = iv;
    }
    else if (gA->Size == xfe_SizeReal){
      rval = gA->rValue[0];
      for (k=0; k<ntot; k++) rval[k] = rv;
    }
    else return xf_Error(xf_NOT_SUPPORTED);	
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_VectorAbs
int 
xf_VectorAbs( xf_Vector *A){
  // A = abs(A), component-wise
  int i;
  xf_long k, ntot;
  int *ival;
  real *rval;
  xf_GenArray *gA;

  for (i = 0; i<A->nArray; i++){
    gA = A->GenArray+i;
    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    if (gA->Size == xfe_SizeInt){
      ival = gA->iValue[0];
      for (k=0; k<ntot; k++) ival[k] = abs(ival[k]);
    }
    else if (gA->Size == xfe_SizeReal){
      rval = gA->rValue[0];
      for (k=0; k<ntot; k++) rval[k] = fabs(rval[k]);
    }
    else return xf_Error(xf_NOT_SUPPORTED);	
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorMin
int 
xf_VectorMin( xf_Vector *A, int *imin, real *rmin){
  // min stored in ival (for integer) or rval (for real)
  int ierr;
  int i;
  xf_long k, ntot;
  int *ival;
  real *rval;
  xf_GenArray *gA;

  for (i = 0; i<A->nArraySelf; i++){
    gA = A->GenArray+i;
    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    if (gA->Size == xfe_SizeInt){
      if (imin == NULL) return xf_Error(xf_INPUT_ERROR);
      ival = gA->iValue[0];
      for (k=1, (*imin)=ival[0]; k<ntot; k++) (*imin) = min((*imin), ival[k]);
      // reduce-min
      ierr = xf_Error(xf_MPI_Allreduce(imin, 1, xfe_SizeInt, xfe_MPI_MIN));
      if (ierr != xf_OK) return ierr;
    }
    else if (gA->Size == xfe_SizeReal){
      if (rmin == NULL) return xf_Error(xf_INPUT_ERROR);
      rval = gA->rValue[0];
      for (k=1, (*rmin)=rval[0]; k<ntot; k++) (*rmin) = min((*rmin), rval[k]);
      // reduce-min
      ierr = xf_Error(xf_MPI_Allreduce(rmin, 1, xfe_SizeReal, xfe_MPI_MIN));
      if (ierr != xf_OK) return ierr;
    }
    else return xf_Error(xf_NOT_SUPPORTED);	
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorMask
int 
xf_VectorMask( xf_Vector *A, xf_Vector *Mask, int MaskVal){
  // Components of A for which Mask is not equal MaskVal
  // are set to zero
  int i;
  xf_long k, ntot;
  int *ival;
  real *rval;
  xf_GenArray *gA;

  for (i = 0; i<A->nArray; i++){
    gA = A->GenArray+i;
    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    if (gA->Size == xfe_SizeInt){
      ival = gA->iValue[0];
      for (k=0; k<ntot; k++) 
	if (Mask->GenArray[i].iValue[k][0] != MaskVal) ival[k] = 0;
    }
    else if (gA->Size == xfe_SizeReal){
      rval = gA->rValue[0];
      for (k=0; k<ntot; k++) 
	if (Mask->GenArray[i].iValue[k][0] != MaskVal) rval[k] = 0.;
    }
    else return xf_Error(xf_NOT_SUPPORTED);	
  }
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_VectorDot
int 
xf_VectorDot( const xf_Vector *A, const xf_Vector *B, real *pdp)
{
  // dp = dot(A,B)
  int ierr, i;
  xf_long k, j, ntot;
  real *rvalA, *rvalB, dp, *rA, *rB;
  xf_GenArray *gA, *gB;

  if (!xf_CompatibleVectors(A,B)) return xf_Error(xf_INCOMPATIBLE);

  dp = 0.0;
  for (i = 0; i<A->nArraySelf; i++){
    gA = A->GenArray+i;
    gB = B->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);

    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rvalA = gA->rValue[0];
    rvalB = gB->rValue[0];
      
    j = 0;
    for (k=0; k<ntot/4; k++, j+=4) {
      rA = rvalA+j;
      rB = rvalB+j;
      dp += rA[0]*rB[0] + rA[1]*rB[1] + rA[2]*rB[2] + rA[3]*rB[3];
    }
    for (k=j; k<ntot; k++)
      dp += rvalA[k] * rvalB[k];

  }

  // reduce-sum
  ierr = xf_Error(xf_MPI_Allreduce(&dp, 1, xfe_SizeReal, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;


  (*pdp) = dp;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorMult
int 
xf_VectorMult( xf_Vector *A, real c)
{
  // A *= c
  int i;
  xf_long k, ntot;
  real *rvalA;
  xf_GenArray *gA;

  if (c == 1.0) return xf_OK; // nothing to be done

  for (i = 0; i<A->nArraySelf; i++){
    gA = A->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);

    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rvalA = gA->rValue[0];
    
    for (k=0; k<ntot; k++) rvalA[k] *= c;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorVectorMult
int 
xf_VectorVectorMult( const xf_Vector *B, xf_Vector *A){
  // A *= B,
  int i;
  xf_long k, ntot;
  real *rvalA, *rvalB;
  xf_GenArray *gA, *gB;

  if (!xf_CompatibleVectors(A,B)) return xf_Error(xf_INCOMPATIBLE);

  for (i = 0; i<A->nArraySelf; i++){
    gA = A->GenArray+i;
    gB = B->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);
    
    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rvalA = gA->rValue[0];
    rvalB = gB->rValue[0];
      
    for (k=0; k<ntot; k++) rvalA[k]  *=  rvalB[k];
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorInv
int 
xf_VectorInv( xf_Vector *A)
{
  // A = 1/A
  int i;
  xf_long k, ntot;
  real *rvalA;
  xf_GenArray *gA;

  for (i = 0; i<A->nArraySelf; i++){
    gA = A->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);

    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rvalA = gA->rValue[0];
    
    for (k=0; k<ntot; k++){
      if (rvalA[k] == 0.) return xf_Error(xf_SINGULAR);
      rvalA[k] = 1./rvalA[k];
    } // k
  } // i
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorAdd
int 
xf_VectorAdd( xf_Vector *A, real c)
{
  // A += c
  int i;
  xf_long k, ntot;
  real *rvalA;
  xf_GenArray *gA;

  if (c == 0.0) return xf_OK; // nothing to be done

  for (i = 0; i<A->nArraySelf; i++){
    gA = A->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);

    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rvalA = gA->rValue[0];
    
    for (k=0; k<ntot; k++) rvalA[k] += c;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_VectorMultSet
int 
xf_VectorMultSet( const xf_Vector *B, real c, enum xfe_AddType AddFlag, xf_Vector *A){
  // A &= c*B, where &= is one of {=, +=, -=, =-}
  int i;
  xf_long k, ntot;
  real *rvalA, *rvalB;
  xf_GenArray *gA, *gB;

  if (!xf_CompatibleVectors(A,B)) return xf_Error(xf_INCOMPATIBLE);

  for (i = 0; i<A->nArraySelf; i++){
    gA = A->GenArray+i;
    gB = B->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);

    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rvalA = gA->rValue[0];
    rvalB = gB->rValue[0];
      
    switch (AddFlag){
    case xfe_Set:
      for (k=0; k<ntot; k++) rvalA[k]  =  c*rvalB[k];
      break;
    case xfe_Neg:
      for (k=0; k<ntot; k++) rvalA[k]  = -c*rvalB[k];
      break;
    case xfe_Add:
      for (k=0; k<ntot; k++) rvalA[k] +=  c*rvalB[k];
      break;
    case xfe_Sub:
      for (k=0; k<ntot; k++) rvalA[k] -=  c*rvalB[k];
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    }
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorRand
int 
xf_VectorRand(xf_Vector *A, int PseudoStep)
{
  // A = random between [0,1]; or equal increments (PseudoStep > 0)
  int i, j=0;
  xf_long k, ntot;
  real *rvalA;
  xf_GenArray *gA;
  
  for (i = 0; i<A->nArraySelf; i++){
    gA = A->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);

    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rvalA = gA->rValue[0];
    
    if (PseudoStep <= 0)
      for (k=0; k<ntot; k++)
	rvalA[k] = ((real) rand())/((real) RAND_MAX);
    else
      for (k=0; k<ntot; k++, j = (++j)%PseudoStep) 
	rvalA[k] = ((real) j) / ((real) PseudoStep);
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MatrixInvert
int 
xf_MatrixInvert( xf_Matrix *M, int nn, xf_Matrix *iM)
{
  // iM = inv(M), 
  int ierr, i, k;
  int *P;
  real *T;
  xf_GenArray *gM, *giM;

  gM  =  M->GenArray;
  giM = iM->GenArray;

  if (gM->Size != xfe_SizeReal) return xf_Error(xf_NOT_SUPPORTED);
  if (gM->n != giM->n) xf_Error(xf_INPUT_ERROR);
  if (gM->r != giM->r) xf_Error(xf_INPUT_ERROR);
  if (gM->r != nn*nn) xf_Error(xf_INPUT_ERROR);

  // P is a permutation vector for PLU
  ierr = xf_Error(xf_Alloc((void **) &P, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // M should not be modified, so use T
  ierr = xf_Error(xf_Alloc((void **) &T, nn*nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<gM->n; i++){
    for (k=0; k<nn*nn; k++) T[k] = gM->rValue[i][k];
    ierr = xf_Error(xf_ComputePLU(T, nn, P));
    if (ierr != xf_OK) return ierr;
    for (k=0; k<nn*nn; k++      ) giM->rValue[i][k] = 0.0;
    for (k=0; k<nn*nn; k+=(nn+1)) giM->rValue[i][k] = 1.0;
    ierr = xf_Error(xf_SolvePLU_Matrix(T, P, nn, nn, giM->rValue[i]));
    if (ierr != xf_OK) return ierr;
  }

  xf_Release( (void *) T);
  xf_Release( (void *) P);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_VectorNorm
int 
xf_VectorNorm( const xf_Vector *A, const int p, real *pnorm)
{
  // (*norm) = p-norm of A
  int i, j, ierr, rr;
  xf_long k, ntot;
  real norm;
  real *rval;
  xf_GenArray *gA;

  if ((p < 0) || (p > 2)) return xf_Error(xf_NOT_SUPPORTED);
  
  norm = 0.0;
  for (i = 0; i<A->nArraySelf; i++){ // only loop over arrays on own proc
    gA = A->GenArray+i;
    if (gA->Size != xfe_SizeReal) return xf_Error(xf_INPUT_ERROR);

    if (gA->vr == NULL)
      ntot = (xf_long) gA->n*gA->r;
      //if only consider one state variable
    //{
    //   ntot = (xf_long) gA->n;
    //   rr = gA->r / 5;
    //}
    else
      for (k=0, ntot=0; k<gA->n; k++) ntot += gA->vr[k];
    if (ntot == 0) continue;

    rval = gA->rValue[0];
    if (p == 0)
      for (k=0; k<ntot; k++) norm += rval[k];
      //for (k=0; k<ntot; k++) 
      //   for (j=0; j<rr; j++)
      //   norm += rval[k*gA->r + j*5];
    else if (p == 1)
      for (k=0; k<ntot; k++) norm += fabs(rval[k]);
      //for (k=0; k<ntot; k++) 
      //   for (j=0; j<rr; j++)
      //   norm += fabs(rval[k*gA->r + j*5]);
    else
      for (k=0; k<ntot; k++) norm += rval[k]*rval[k];
      //for (k=0; k<ntot; k++) 
      //   for (j=0; j<rr; j++)
      //   norm += rval[k*gA->r + j*5]*rval[k*gA->r + j*5];
  }

  // reduce-sum
  ierr = xf_Error(xf_MPI_Allreduce(&norm, 1, xfe_SizeReal, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;

  if (p == 2) norm = sqrt(norm);

  (*pnorm) = norm;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_VectorStateNorms
int 
xf_VectorStateNorms(xf_All *All, xf_Vector *A, real *StateNorms,
                    int p)
{
  int ierr, nn, sr, n, r;
  int elem, egrp;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;
  
  if (p != 2) return xf_Error(xf_NOT_SUPPORTED);
  
  for (r = 0; r < sr; r++)
    StateNorms[r] = 0.0;
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    if (A->GenArray[egrp].Size != xfe_SizeReal)
      return xf_Error(xf_INPUT_ERROR);
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      if (A->GenArray[egrp].vr == NULL)
        nn = A->GenArray[egrp].r/sr;
      else 
        nn = A->GenArray[egrp].vr[elem]/sr;
      
      for (n = 0; n < nn; n++) 
        for (r = 0; r < sr; r++) 
          StateNorms[r] += pow(A->GenArray[egrp].rValue[elem][n*sr+r],2.0);
    }
  }
  ierr = xf_Error(xf_MPI_Allreduce(StateNorms, sr, xfe_SizeReal, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;

  for (r = 0; r < sr; r++)
    StateNorms[r] = sqrt(StateNorms[r]);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelReAllocElemArray
static int 
xf_ParallelReAllocElemArray( xf_GenArray *ga)
{
  int ierr, j;
  int myRank, nProc, iProc;
  int *nSend;
  xf_ArrayParallelInfo *PInfo;
  
  if ((PInfo = ga->ParallelInfo) == NULL) return xf_OK; // nothing to do

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // if not parallel, return immediately
  if (nProc == 1) return xf_OK;

  // ReAllocate send buffers if not on halo
  if (!PInfo->HaloFlag){
    ierr = xf_Error(xf_Alloc( (void **) &nSend, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    for (iProc=0; iProc<nProc; iProc++){
      if (ga->vr == NULL)
	nSend[iProc] = PInfo->nSendElem[iProc]*ga->r;
      else // account for variable orders
	for (j=0, nSend[iProc]=0; j<PInfo->nSendElem[iProc]; j++)
	  nSend[iProc] += ga->vr[PInfo->SendElem[iProc][j]];
    }
    
    switch (ga->Size){
    case xfe_SizeInt:
      ierr = xf_Error(xf_VReAlloc2((void ***) &PInfo->iSendBuf, nProc, nSend, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_SizeReal:
      ierr = xf_Error(xf_VReAlloc2((void ***) &PInfo->rSendBuf, nProc, nSend, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    }
    
    xf_Release( (void *) nSend);
  }
 
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ProjectVector
int 
xf_ProjectVector( xf_All *All, const xf_Vector *A, 
		  enum xfe_Bool TransposeFlag, xf_Vector *B)
{
  int ierr, i, j, k, sr;
  int n1, n2;
  int  OrderA,  OrderB;
  int pOrderA, pOrderB;
  enum xfe_Bool VariableOrderA, VariableOrderB, ConstantOrder, ReCalc;
  enum xfe_BasisType BasisA, BasisB;
  real *TT;
  xf_Matrix *T = NULL;

  if (A->nArray <= 0){
    xf_printf("Warning, no data to project.\n");
    return xf_OK;
  }

  if ((A->Basis == NULL) || (A->Order == NULL) ||
      (B->Basis == NULL) || (B->Order == NULL)){
    xf_printf("Error, vectors to project are not interpolated.\n");
    return xf_Error(xf_INPUT_ERROR);
  }

  if ( ((A->Linkage != xfe_LinkageGlobElem) && (A->Linkage != xfe_LinkageElem)) || 
       ((B->Linkage != xfe_LinkageGlobElem) && (B->Linkage != xfe_LinkageElem)) ){
    xf_printf("Only element-linked vectors can be projected.\n");
    return xf_Error(xf_INPUT_ERROR);
  }

  // make sure StateRank matches
  if ((sr = A->StateRank) != B->StateRank) return xf_Error(xf_INPUT_ERROR);
  
  // flags denoting variable order
  VariableOrderA = ((A->nComp != NULL) && (A->vOrder != NULL));
  VariableOrderB = ((B->nComp != NULL) && (B->vOrder != NULL));
  ConstantOrder  = ((!VariableOrderA)  && (!VariableOrderB));

  // loop over and project each array
  for (i=0; i<A->nArray; i++){

    BasisA = A->Basis[i];  OrderA = A->Order[i];
    BasisB = B->Basis[i];  OrderB = B->Order[i];

    if ((ConstantOrder) && (BasisA == BasisB) && (OrderA == OrderB)){
      // set B[i] = A[i] and continue
      for (j=0; j<A->GenArray[i].n; j++)
	for (k=0; k<A->GenArray[i].r; k++)
	  B->GenArray[i].rValue[j][k] = A->GenArray[i].rValue[j][k];
      continue;
    }
    
    pOrderA = -1;
    pOrderB = -1;

    // begin loop over elements
    for (j=0; j<A->GenArray[i].n; j++){

      // should we recalculate transfer matrix?
      OrderA = ((VariableOrderA) ? A->vOrder[i][j] : A->Order[i]);
      OrderB = ((VariableOrderB) ? B->vOrder[i][j] : B->Order[i]);
      ReCalc = ((OrderA != pOrderA) || (OrderB != pOrderB));

      if (ReCalc){ // yes, recalculate
	// obtain n1 = current # unknowns
	ierr = xf_Error(xf_Order2nNode(BasisA, OrderA, &n1));
	if (ierr != xf_OK) return ierr;
	
	// obtain n2 = desired # unknowns
	ierr = xf_Error(xf_Order2nNode(BasisB, OrderB, &n2));
	if (ierr != xf_OK) return ierr;

	// locate transfer matrix from current basis/order to desired basis/order
	if (TransposeFlag){
	  // transpose of reverse operator is requested
	  ierr = xf_Error(xf_FindTransferMatrix(NULL, BasisB, OrderB, BasisA, OrderA, &T));
	  if (ierr != xf_OK) return ierr;
	}
	else{
	  ierr = xf_Error(xf_FindTransferMatrix(NULL, BasisA, OrderA, BasisB, OrderB, &T));
	  if (ierr != xf_OK) return ierr;
	}

	pOrderA = OrderA;
	pOrderB = OrderB;
      }
      
      // pull off real array of the matrix
      TT = T->GenArray->rValue[0];
      
      //  apply transformation to set rValue
      if (TransposeFlag)
	xf_MTxM_Set(TT, A->GenArray[i].rValue[j], n2, n1, sr, B->GenArray[i].rValue[j]);
      else
	xf_MxM_Set(TT, A->GenArray[i].rValue[j], n2, n1, sr, B->GenArray[i].rValue[j]);
    
    } // j (over elements)
    
  } // i (over arrays)

  // T was created stand-alone and must be destroyed
  ierr = xf_Error(xf_DestroyMatrix(T));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetDesiredOrder
static int
xf_GetDesiredOrder(int InOrder, int OrderScal, int OrderIncrement)
{
  // used by xf_ProjectVectorInPlace
  int Order;

  Order = ((OrderScal >= 0) ? OrderScal : InOrder + OrderIncrement);
  if (Order < 0) Order = InOrder;

  return Order;
}

/******************************************************************/
//   FUNCTION Definition: xf_PrepVectorVariableOrder
static int 
xf_PrepVectorVariableOrder( xf_Vector *V)
{
  // prepares V for variable order storage
  int ierr, i, j;

  if ((V->nComp != NULL) && (V->vOrder != NULL)) return xf_OK; // already variable order

  // do not support "half-way" structure allocations
  if (V->nComp  != NULL) return xf_Error(xf_INPUT_ERROR);
  if (V->vOrder != NULL) return xf_Error(xf_INPUT_ERROR);

  // allocate and set nComp
  ierr = xf_Error(xf_Alloc((void **) &V->nComp, V->nArray, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<V->nArray; i++) V->nComp[i] = V->GenArray[i].n;

  // allocate and set vOrder
  ierr = xf_Error(xf_VAlloc2((void ***) &V->vOrder, V->nArray, V->nComp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<V->nArray; i++)
    for (j=0; j<V->nComp[i]; j++)
      V->vOrder[i][j] = V->Order[i];

  // loop over arrays
  for (i=0; i<V->nArray; i++){
    // allocate and set vr
    ierr = xf_Error(xf_Alloc( (void **) &V->GenArray[i].vr, V->nComp[i], sizeof(int))); 
    if (ierr != xf_OK) return ierr;
    for (j=0; j<V->nComp[i]; j++) V->GenArray[i].vr[j] = V->GenArray[i].r;
    // no need to adjust .rValue or .iValue (already allocated adequately)
  } // i
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ProjectVectorInPlace_General
static int 
xf_ProjectVectorInPlace_General( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
				 enum xfe_BasisType *BasisVec, enum xfe_BasisType BasisScal, 
				 xf_Vector *VOrder, int **vOrder, int *OrderVec, int OrderScal,
				 int OrderIncrement)
{
  // description under specific wrapper functions in xf_DataMath.h
  int ierr, i, j, sr, k, negrp;
  int n1, n2, r1, r2, r, Order;
  int pOrder, pCurrOrder, CurrOrder, ReqOrder;
  int *vr = NULL;
  enum xfe_Bool VariableOrder;
  enum xfe_BasisType Basis;
  real *TT, **rValue, **rtemp;
  xf_Matrix *T = NULL;

  if (V->nArray <= 0){
    xf_printf("Warning, no data to project.\n");
    return xf_OK;
  }

  if ((V->Basis == NULL) || (V->Order == NULL)){
    xf_printf("Error, vector to project is not interpolated.\n");
    return xf_Error(xf_INPUT_ERROR);
  }

  if ((V->Linkage != xfe_LinkageGlobElem) && (V->Linkage != xfe_LinkageElem)){
    xf_printf("Only element-linked vectors can be projected.\n");
    return xf_Error(xf_INPUT_ERROR);
  }

  sr = V->StateRank;

  // if VOrder or vorder are specified, make sure V is ready for variable order
  if ((VOrder != NULL) || (vOrder != NULL)){
    ierr = xf_Error(xf_PrepVectorVariableOrder(V));
    if (ierr != xf_OK) return ierr;
  }

  // are we dealing with variable orders?
  VariableOrder = ((V->nComp != NULL) && (V->vOrder != NULL));

  // if mesh is parallel,input basis and order may not be defined separately for halo
  negrp = ((Mesh != NULL) ? Mesh->nElemGroup : V->nArray);

  // loop over and project each array
  for (i=0; i<V->nArray; i++){

    // requested order
    ReqOrder = ((OrderVec==NULL) ? OrderScal : OrderVec[i%negrp]);

    // determine basis
    Basis = ((BasisVec == NULL) ? BasisScal : BasisVec[i%negrp]);
    if (Basis == xfe_BasisLast) Basis = V->Basis[i];

    if (!VariableOrder){ // constant order vector
      // determine order
      Order = xf_GetDesiredOrder(V->Order[i], ReqOrder, OrderIncrement);
      // continue if nothing to do
      if ((V->Basis[i] == Basis) && (V->Order[i] == Order)) continue;
    }

    // verify that we're dealing with a real vector
    if (V->GenArray[i].Size != xfe_SizeReal) return xf_Error(xf_OUT_OF_BOUNDS);


    if (VariableOrder){ // in case of variable order ...

      // allocate a new vr
      ierr = xf_Error(xf_Alloc((void **) &vr, V->GenArray[i].n, sizeof(int)));
      if (ierr != xf_OK) return ierr;

      // fill in the new vr
      for (j=0; j<V->GenArray[i].n; j++){

	// requested order is element-specific when VOrder is specified 
	// or when vOrder is specified
	if (VOrder != NULL) ReqOrder = VOrder->GenArray[i].iValue[j][0];
	else if (vOrder != NULL) ReqOrder = vOrder[i][j];

	// current order
	CurrOrder = ((VariableOrder) ? V->vOrder[i][j] : V->Order[i]);
	
	// desired order
	Order = xf_GetDesiredOrder(CurrOrder, ReqOrder, OrderIncrement);
	
	// r2 = desired # unknowns * sr per element
	ierr = xf_Error(xf_Order2nNode(Basis, Order, &n2));
	if (ierr != xf_OK) return ierr;
	r2 = n2*sr;
	
	// set value in vr
	vr[j] = r2;

      } // j
      
      // set V->GenArray[i].r = max(vr)
      V->GenArray[i].r = 0;
      for (j=0; j<V->GenArray[i].n; j++)
	V->GenArray[i].r = max(V->GenArray[i].r, vr[j]);
      
      // allocate new rValue with vr
      ierr = xf_Error(xf_VAlloc2((void ***) &rValue, V->GenArray[i].n, vr, sizeof(real)));
      if (ierr != xf_OK) return ierr;

      // no longer need vr
      xf_Release( (void *) vr);

    }
    else{ // for constant order ...

      // r2 = desired # unknowns * sr per element
      ierr = xf_Error(xf_Order2nNode(Basis, Order, &n2));
      if (ierr != xf_OK) return ierr;
      r2 = n2*sr;

      // allocate an rValue with r2
      ierr = xf_Error(xf_Alloc2((void ***) &rValue, V->GenArray[i].n, r2, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }

    pOrder = -1;
    pCurrOrder = -1;

    // loop over elements and project data
    for (j=0; j<V->GenArray[i].n; j++){
      
      // requested order is element-specific when VOrder is specified 
      // or when vorder is specified
      if (VOrder != NULL) ReqOrder = VOrder->GenArray[i].iValue[j][0];
      else if (vOrder != NULL) ReqOrder = vOrder[i][j];
   
      // current order
      CurrOrder = ((VariableOrder) ? V->vOrder[i][j] : V->Order[i]);

      // desired order
      Order = xf_GetDesiredOrder(CurrOrder, ReqOrder, OrderIncrement);
      
      // obtain r1 = current # unknowns * sr per element
      ierr = xf_Error(xf_Order2nNode(V->Basis[i], CurrOrder, &n1));
      if (ierr != xf_OK) return ierr;
      r1 = n1*sr;

      // check r1 against existing r (should be same)
      r = ((VariableOrder) ? V->GenArray[i].vr[j] : V->GenArray[i].r);
      if (r != r1) return xf_Error(xf_INPUT_ERROR);
      
      // obtain r2 = desired # unknowns * sr per element
      ierr = xf_Error(xf_Order2nNode(Basis, Order, &n2));
      if (ierr != xf_OK) return ierr;
      r2 = n2*sr;

      // recalculate transfer matrix if order changes
      if ((pOrder != Order) || (pCurrOrder != CurrOrder)){
	pOrder = Order;
	pCurrOrder = CurrOrder;

	// locate transfer matrix from current basis/order to desired basis/order
	ierr = xf_Error(xf_FindTransferMatrix(DataSet, V->Basis[i], CurrOrder, Basis, Order, &T));
	if (ierr != xf_OK) return ierr;
      }

      // pull off real array of the matrix
      TT = T->GenArray->rValue[0];

      if ((T->GenArray->n != 1) || (T->GenArray->r != n1*n2)){
	xf_printf("T->GenArray->r = %d\n", T->GenArray->r);
	xf_printf("n1=%d, n2=%d, n1*n2 = %d\n", n1, n2, n1*n2);
	return xf_Error(xf_OUT_OF_BOUNDS);
      }

      // apply transformation to set rValue
      xf_MxM_Set(TT, V->GenArray[i].rValue[j], n2, n1, sr, rValue[j]);

      // set .vr[j] and Order if variable order
      if (VariableOrder){
	V->GenArray[i].vr[j] = r2;
	V->vOrder[i][j] = Order;
      }

    } // j over elements

    // set .Order to max(.vOrder)
    if ((VariableOrder) && (V->GenArray[i].n > 0)){
      V->Order[i]=0;
      for (j=0; j<V->GenArray[i].n; j++)
	V->Order[i] = max(V->Order[i], V->vOrder[i][j]);
    }

    // special consistency enforcement in case of zero elements
    if (V->GenArray[i].n == 0){
      ierr = xf_Error(xf_Order2nNode(V->Basis[i], V->Order[i], &V->GenArray[i].r));
      if (ierr != xf_OK) return ierr;
      V->GenArray[i].r *= sr;
    }
	  
    // swap rValue with V->GenArray[i].rValue
    swap(rValue, V->GenArray[i].rValue, rtemp);

    // release rValue
    xf_Release2( (void **) rValue);

    // set V->GenArray[i].r = r2, new Basis, and new Order
    V->Basis[i] = Basis;
    if (!VariableOrder){
      V->GenArray[i].r = r2;
      V->Order[i] = Order;
    }

    // reallocate parallel buffers
    if (V->GenArray[i].ParallelInfo != NULL){
      ierr = xf_Error(xf_ParallelReAllocElemArray(V->GenArray+i));
      if (ierr != xf_OK) return ierr;
    }

    
  } // i

  // Null DataSet means T was created stand-alone and must be destroyed
  if (DataSet == NULL){ 
    ierr = xf_Error(xf_DestroyMatrix(T));
    if (ierr != xf_OK) return ierr;
  }
  

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ProjectVectorInPlace
int 
xf_ProjectVectorInPlace( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
			 enum xfe_BasisType *BasisVec, enum xfe_BasisType BasisScal, 
			 xf_Vector *VOrder, int *OrderVec, int OrderScal,
			 int OrderIncrement)
{
  // wrapper for general function
  return xf_Error(xf_ProjectVectorInPlace_General(Mesh, DataSet, V, BasisVec, BasisScal,
						  VOrder, NULL, OrderVec, OrderScal, 
						  OrderIncrement));
}


/******************************************************************/
//   FUNCTION Definition: xf_ProjectVectorInPlace_Basis
int 
xf_ProjectVectorInPlace_Basis( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
			       enum xfe_BasisType *BasisVec, enum xfe_BasisType BasisScal)
{
  // wrapper for general function
  return xf_Error(xf_ProjectVectorInPlace_General(Mesh, DataSet, V, BasisVec, BasisScal,
						  NULL, NULL, NULL, -1, 
						  0));
}


/******************************************************************/
//   FUNCTION Definition: xf_ProjectVectorInPlace_OrderSet
int 
xf_ProjectVectorInPlace_OrderSet( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
				  enum xfe_BasisType *BasisVec, enum xfe_BasisType BasisScal, 
				  int **vOrder, int *OrderVec, int OrderScal)
{
  // wrapper for general function
  return xf_Error(xf_ProjectVectorInPlace_General(Mesh, DataSet, V, BasisVec, BasisScal,
						  NULL, vOrder, OrderVec, OrderScal, 
						  0));
}


/******************************************************************/
//   FUNCTION Definition: xf_ProjectVectorInPlace_OrderIncrement
int 
xf_ProjectVectorInPlace_OrderIncrement( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
					enum xfe_BasisType *BasisVec, 
					enum xfe_BasisType BasisScal, int OrderIncrement)
{
  // wrapper for general function
  return xf_Error(xf_ProjectVectorInPlace_General(Mesh, DataSet, V, BasisVec, BasisScal,
						  NULL, NULL, NULL, -1, OrderIncrement));
}

/******************************************************************/
//   FUNCTION Definition: xf_ProjectVectorInPlace_VOrder
int 
xf_ProjectVectorInPlace_VOrder( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
				enum xfe_BasisType *BasisVec, 
				enum xfe_BasisType BasisScal, xf_Vector *VOrder)
{
  // wrapper for general function
  return xf_Error(xf_ProjectVectorInPlace_General(Mesh, DataSet, V, BasisVec, BasisScal,
						  VOrder, NULL, NULL, -1, 0));
}


/******************************************************************/
//   FUNCTION Definition: xf_ProjectVectorInPlace_Vector
int 
xf_ProjectVectorInPlace_Vector( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
				xf_Vector *U)
{
  // wrapper for general function
  return xf_Error(xf_ProjectVectorInPlace_General(Mesh, DataSet, V, U->Basis, xfe_BasisLast,
						  NULL, U->vOrder, U->Order, -1, 0));
}

/******************************************************************/
//   FUNCTION Definition: xf_ProjectVectors_VOrderFile
int 
xf_ProjectVectors_VOrderFile( xf_Mesh *Mesh, xf_DataSet *DataSet, int nVector,
			      xf_Vector **Vi, const char *fname)
{
  int ierr, terr, i;
  xf_DataSet *DS = NULL;
  xf_Vector *VOrder = NULL;
  xf_Data *D;

  // does dynamic refinement exist on disk for this time slab?
  ierr = xf_Error(xf_CreateDataSet(&DS));
  if (ierr != xf_OK) return ierr;
  ierr = xf_ReadDataSetBinary(Mesh, NULL, fname, DS);
  if (ierr == xf_NOT_FOUND){
    terr = xf_Error(xf_DestroyDataSet(DS));
    if (terr != xf_OK) return terr;
  }
  if (ierr != xf_OK) return ierr;

  // use first piece of data as VOrder
  D = DS->Head;
  VOrder =  (xf_Vector *) D->Data;
  
  // project vectors in place
  for (i=0; i<nVector; i++){
    ierr = xf_Error(xf_ProjectVectorInPlace(Mesh, DataSet, Vi[i], Vi[i]->Basis, xfe_BasisLast,
					    VOrder, NULL, -1, 0));
    if (ierr != xf_OK) return ierr;
  }

  // destroy dataset
  ierr = xf_Error(xf_DestroyDataSet(DS));
  if (ierr != xf_OK) return ierr;
}




/******************************************************************/
//   FUNCTION Definition: xf_ProjectJacobian
int 
xf_ProjectJacobian( xf_All *All, const xf_JacobianMatrix *R_UA,
		    xf_JacobianMatrix *R_UB)
{
  int ierr, egrp, elem, face, sr;
  int egN, eN, negrp, negrphalo;
  int n1A, n1B, n2A, n2B, Tsize = -1;
  int Order1A, Order1B, Order2A, Order2B;
  int pOrder1A, pOrder1B, pOrder2A, pOrder2B;
  enum xfe_Bool VariableOrderA = xfe_False, VariableOrderB = xfe_False;
  enum xfe_BasisType Basis1A, Basis1B, Basis2A, Basis2B;
  enum xfe_BasisType pBasis1A, pBasis1B, pBasis2A, pBasis2B;
  real *TT1, *TT2, *TTemp = NULL;
  xf_Matrix *T1 = NULL, *T2 = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  negrp = Mesh->nElemGroup;

  if ((sr = R_UA->StateRank) != R_UB->StateRank) return xf_Error(xf_INPUT_ERROR);
  
  // check for variable order in either R_UA or R_UB
  if (R_UA->vnvec != NULL){
    VariableOrderA = xfe_True;
    if ((R_UA->U == NULL) || (R_UA->U->nComp == NULL) || (R_UA->U->vOrder == NULL))
      return xf_Error(xf_INPUT_ERROR);
  }
  if (R_UB->vnvec != NULL){
    VariableOrderB = xfe_True;
    if ((R_UB->U == NULL) || (R_UB->U->nComp == NULL) || (R_UB->U->vOrder == NULL))
      return xf_Error(xf_INPUT_ERROR);
  }

  /* Halo element groups are present in parallel runs */  
  negrphalo = ((Mesh->ParallelInfo == NULL) ? Mesh->nElemGroup : 2*Mesh->nElemGroup);

  // initialize previous values so that TransferMatrix is called first time through
  pBasis1A = pBasis1B = pBasis2A = pBasis2B = xfe_BasisLast;
  pOrder1A = pOrder1B = pOrder2A = pOrder2B = -1;

  for (egrp = 0; egrp<negrphalo; egrp++){
    Basis1A = R_UA->Basis[egrp%negrp];  Order1A = R_UA->Order[egrp%negrp];
    Basis1B = R_UB->Basis[egrp%negrp];  Order1B = R_UB->Order[egrp%negrp];

    n1A = R_UA->nvec[egrp];
    n1B = R_UB->nvec[egrp];

    for (elem = 0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // pull off variable orders
      if (VariableOrderA) Order1A = R_UA->U->vOrder[egrp][elem];
      if (VariableOrderB) Order1B = R_UB->U->vOrder[egrp][elem];

      // set n1A and n1B for variable order
      if (VariableOrderA) n1A = R_UA->vnvec[egrp][elem];
      if (VariableOrderB) n1B = R_UB->vnvec[egrp][elem];

      for (face = -1; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
	if (face == -1){
	  egN = egrp;
	  eN  = elem;
	}
	else{
	  egN = R_UA->egrpN[egrp][elem][face];
	  eN  = R_UA->elemN[egrp][elem][face];
	}
	if (egN < 0) continue;

	Basis2A = R_UA->Basis[egN%negrp];  Order2A = R_UA->Order[egN%negrp];
	Basis2B = R_UB->Basis[egN%negrp];  Order2B = R_UB->Order[egN%negrp];

	n2A = R_UA->nvec[egN];
	n2B = R_UB->nvec[egN];
	
	// pull off variable orders on neighbor
	if (VariableOrderA) Order2A = R_UA->U->vOrder[egN][eN];
	if (VariableOrderB) Order2B = R_UB->U->vOrder[egN][eN];

	// set n2A and n2B for variable order
	if (VariableOrderA) n2A = R_UA->vnvec[egN][eN];
	if (VariableOrderB) n2B = R_UB->vnvec[egN][eN];
	
	// make sure have big enough TTemp 
	if (n1A*n2B*sr*sr > Tsize){
	  Tsize = n1A*n2B*sr*sr;
	  ierr = xf_Error(xf_ReAlloc((void **) &TTemp, Tsize, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}


	// locate prolongation transfer matrix from 1B to 1A
	if ((Basis1A != pBasis1A) || (Basis1B != pBasis1B) ||
	    (Order1A != pOrder1A) || (Order1B != pOrder1B)){
	  ierr = xf_Error(xf_FindTransferMatrix(NULL, Basis1B, Order1B, Basis1A, Order1A, &T1));
	  if (ierr != xf_OK) return ierr;
	  pBasis1A = Basis1A; pBasis1B = Basis1B; pOrder1A = Order1A; pOrder1B = Order1B;
	}

	// pull off real array of the matrix
	TT1 = T1->GenArray->rValue[0];
	
	// locate prolongation transfer matrix from 2B to 2A
	if ((Basis2A != pBasis2A) || (Basis2B != pBasis2B) ||
	    (Order2A != pOrder2A) || (Order2B != pOrder2B)){
	  ierr = xf_Error(xf_FindTransferMatrix(NULL, Basis2B, Order2B, Basis2A, Order2A, &T2));
	  if (ierr != xf_OK) return ierr;
	  pBasis2A = Basis2A; pBasis2B = Basis2B; pOrder2A = Order2A; pOrder2B = Order2B;
	}

	// pull off real array of the matrix
	TT2 = T2->GenArray->rValue[0];


	// TTemp = R_UA*TT2;  R_UA is block n1A x n2A; TT2 is n2A x n2B; TTemp is block n1A x n2B
	xf_BlockMxM(R_UA->Value[egrp][elem][1+face], n1A, sr, n2A, TT2, 
		    n2B, xfe_Set, TTemp);


	// R_UB = TT1^T*TTemp;  TT1^T is n1B x n1A; R_UB is block n1B x n2B
	xf_MTxBlockM(TT1, n1B, n1A, TTemp, sr, n2B, xfe_Set, 
		     R_UB->Value[egrp][elem][1+face]);

      } // face
    } // elem
  } // egrp


  // T1 and T2 were created stand-alone and must be destroyed
  ierr = xf_Error(xf_DestroyMatrix(T1));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyMatrix(T2));
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) TTemp);

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_VectorSetUpperHouse
static int 
xf_VectorSetUpperHouse( xf_VectorSet *VS, int n, real *W, real *R)
{
/*
  Calculation performed in parallel.  VS is overwritten.  Total number
  of vector entries should be greater than n.  The vector entries are
  ordered along processors first, then arrays, then elements, then
  ranks.
*/
  int ierr, i, j, k, r, is;
  int myRank, nProc, iProc;
  int *vtemp, *vproc, *vegrp, *velem, *vr;
  int iproc, iegrp, ielem, ir;
  enum xfe_Bool done, ParallelFlag;
  xf_long nentriesLong, *nEntriesLong, ntot;
  int nentries, *nEntries, istart, *iStart;
  real g, s, scale, f, h, *E, *EV;
  xf_Vector *V, *Vi, *Vj;
  
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  ParallelFlag = (nProc > 1);

  // no support yet for variable orders
  if (VS->Vector[0].vOrder != NULL) return xf_Error(xf_NOT_SUPPORTED);

  /* First, construct vectors identifying ith row with processor
     number (iproc), element group number (iegrp), element number
     (ielem), and rank number within the element (ir). */

  // each proc sums up # entries it has
  V = VS->Vector+0;
  for (i=0, nentriesLong=0; i<V->nArraySelf; i++)
    nentriesLong += V->GenArray[i].n*V->GenArray[i].r;

  // allocate nEntriesLong on 0
  nEntriesLong = NULL;
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc( (void **) &nEntriesLong, nProc, sizeof(xf_long)));
    if (ierr != xf_OK) return ierr;
  }

  // gather number of entries to proc 0
  if (ParallelFlag){
    ierr = xf_Error(xf_MPI_Gather((void *) &nentriesLong, (void *) nEntriesLong, 
				  1*sizeof(xf_long), 0));
    if (ierr != xf_OK) return ierr;
  }
  else{
    nEntriesLong[0] = nentriesLong;
  }

  // each proc allocates vproc: vproc[i] = proc # for row i
  ierr = xf_Error(xf_Alloc( (void **) &vproc, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // proc 0 counts first n entries: stores # for each proc in nEntries
  iStart = NULL;
  nEntries = NULL;
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc( (void **) &iStart, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &nEntries, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    iProc = 0; iStart[0] = 0;  nEntries[0] = (int) nEntriesLong[0];
    ntot = nEntriesLong[iProc];
    while (ntot < (xf_long) n){
      iProc++;
      if (iProc >= nProc){
	xf_printf("Number of vectors exceeds # vector entries.\n");
	return xf_Error(xf_OUT_OF_BOUNDS);
      }
      iStart[iProc] = (int) ntot;
      nEntries[iProc] = (int) nEntriesLong[iProc];
      ntot += nEntriesLong[iProc];
    }
    nEntries[iProc] -= ((int) ntot - n);
    for (i=iProc+1; i<nProc; i++) nEntries[i] = 0;
    for (i=iProc+1; i<nProc; i++) iStart[i] = 0;

    for (i=0; i<nProc; i++)
      for (k=0; k<nEntries[i]; k++)
	vproc[iStart[i]+k] = i;
  }

  // proc 0 scatters iStart and nEntries
  if (ParallelFlag){
    ierr = xf_Error(xf_MPI_Scatter((void *) iStart, (void *) &istart, 
				   1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_MPI_Scatter((void *) nEntries, (void *) &nentries, 
				   1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
  }
  else{
    istart = iStart[0];
    nentries = nEntries[0];
  }

  // proc 0 bcasts vproc
  if (ParallelFlag){
    ierr = xf_Error(xf_MPI_Bcast((void *) vproc, n*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  // each proc builds vegrp, velem, vr
  ierr = xf_Error(xf_Alloc( (void **) &vtemp, 3*n, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<3*n; k++) vtemp[k] = -1;
  vegrp = vtemp; velem = vegrp+n; vr = velem+n;

  if (nentries > 0){
    k = istart;
    done = xfe_False;
    for (i=0; (i<V->nArraySelf) && (!done); i++)
      for (j=0; (j<V->GenArray[i].n) && (!done); j++)
	for (r=0; (r<V->GenArray[i].r) && (!done); r++){
	  vegrp[k] = i;
	  velem[k] = j;
	  vr[k]    = r;
	  k++;
	  if (k >= (istart+nentries)) done = xfe_True;
	} // r
  }

  // Allocate memory
  ierr = xf_Error(xf_Alloc( (void **) &E, n, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  /* for (i=0; i<n; i++){ */
  /*     xf_pprintf("i=%d: iproc = %d, iegrp=%d, ielem=%d, ir = %d\n", */
  /* 	       i, vproc[i], vegrp[i], velem[i], vr[i]); */
  /*   } */
 
  
  // Start main loop
  for (i=0; i<n; i++){
    /* xf_printf("A = \n"); */
    /*     V = VS->Vector+0; */
    /*     for (iegrp=0; iegrp<V->nArraySelf; iegrp++){ */
    /*       for (ielem=0; ielem<V->GenArray[iegrp].n; ielem++){ */
    /* 	for (ir=0; ir<V->GenArray[iegrp].r; ir++){ */
    /* 	  for (j=0; j<n; j++){ */
    /* 	    EV = VS->Vector[j].GenArray[iegrp].rValue[ielem]; */
    /* 	    xf_printf("%.10E ", EV[ir]); */
    /* 	  } */
    /* 	  xf_printf("\n"); */
    /* 	} */
    /*       } */
    /*     } */
    /*     xf_printf("R = \n"); */
    /*     for (j=0; j<i*n; j++){ */
    /*       xf_printf("%.10E ", R[j]); fflush(stdout); */
    /*       if ((j+1)%n == 0) xf_printf("\n"); */
    /*     } */
    
    Vi = VS->Vector + i; // ith vector in set
    iproc = vproc[i]; iegrp = vegrp[i]; ielem = velem[i]; ir = vr[i];
  

    g = s = 0.;
    ierr = xf_Error(xf_VectorNorm(Vi, 1, &scale));
    if (ierr != xf_OK) return ierr;
    
    if (scale > 0){
      // Vi *= 1/scale
      ierr = xf_Error(xf_VectorMult(Vi, 1./scale));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_VectorNorm(Vi, 2, &s));
      if (ierr != xf_OK) return ierr;

      if (myRank == iproc)
	f = Vi->GenArray[iegrp].rValue[ielem][ir];
      if (ParallelFlag){
	ierr = xf_Error(xf_MPI_Bcast((void *) &f, sizeof(real), iproc));
	if (ierr != xf_OK) return ierr;
      }
      g = ((f > 0) ? -s : s);
      h = f*g - s*s;
      if (myRank == iproc)
	Vi->GenArray[iegrp].rValue[ielem][ir] = f-g;
      E[i] = h;
      W[i] = scale*g;

      // apply transformations to remaining columns
      for (j=i+1; j<n; j++){
	Vj = VS->Vector + j;
	ierr = xf_Error(xf_VectorDot(Vi, Vj, &s));
	if (ierr != xf_OK) return ierr;
	f = s/h;
	ierr = xf_Error(xf_VectorMultSet(Vi, f, xfe_Add, Vj));
	if (ierr != xf_OK) return ierr;
      }
    }

    if (iegrp >= 0){ 
      // copy ith row of VS into R (starting at column i)
      // zero out ith row of VS (starting at column i+1)
      is = i-istart;
      if ((is < 0) || (is > nentries)) return xf_Error(xf_CODE_LOGIC_ERROR);
      for (j=0; j<i; j++) R[is*n+j] = 0.;
      for (j=i; j<n; j++){
	EV = VS->Vector[j].GenArray[iegrp].rValue[ielem];
	R[is*n+j] = EV[ir];
	if (j > i) EV[ir] = 0.;
      }
      R[is*n+i] = W[i];
    }

  } // i

  /* Accumulate Householder transformations in VS */
  for (i=n-1; i>=0; i--){
    iproc = vproc[i]; iegrp = vegrp[i]; ielem = velem[i]; ir = vr[i];
    Vi = VS->Vector + i;
    g = W[i];
    h = E[i];
    if (h == 0.0){ // set Vi to 0
      ierr = xf_Error(xf_SetZeroVector(Vi));
      if (ierr != xf_OK) return ierr;
      continue;
    }
    for (j=i+1; j<n; j++){
      Vj = VS->Vector + j;
      ierr = xf_Error(xf_VectorDot(Vi, Vj, &s));
      if (ierr != xf_OK) return ierr;
      f = s/h;
      ierr = xf_Error(xf_VectorMultSet(Vi, f, xfe_Add, Vj));
      if (ierr != xf_OK) return ierr;
    }
    if (myRank == iproc)
      s = Vi->GenArray[iegrp].rValue[ielem][ir]/h;
    if (ParallelFlag){
      ierr = xf_Error(xf_MPI_Bcast((void *) &s, sizeof(real), iproc));
      if (ierr != xf_OK) return ierr;
    }

    // Vi *= s
    ierr = xf_Error(xf_VectorMult(Vi, s));
    if (ierr != xf_OK) return ierr;

    if (myRank == iproc)
      Vi->GenArray[iegrp].rValue[ielem][ir] += 1.0;
  } //i

  /* Gather Ri onto root -> R */
  if (ParallelFlag){
    i = 0;
    if (myRank == 0){
      i = nentries*n;
      nEntries[0] = nentries = 0;
    }
    ierr = xf_Error(xf_MPI_Gatherv((void *) R, nentries, (void *) R+i, nEntries, 
				   n*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
  }

  // release memory
  xf_Release( (void *) vproc);
  xf_Release( (void *) vtemp);
  xf_Release( (void *) nEntriesLong);
  xf_Release( (void *) nEntries);
  xf_Release( (void *) iStart);
  xf_Release( (void *) E);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_VectorSetSVD
int 
xf_VectorSetSVD( xf_VectorSet *VS, int n, real *W)
{
  int ierr;
  int myRank, nProc;
  int egrp, elem, nr, r, i, j;
  enum xfe_Bool ParallelFlag;
  real *R, **EV, *T, *A, s;
  xf_Vector *V;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  ParallelFlag = (nProc > 1);
  
  // Allocate memory
  ierr = xf_Error(xf_Alloc( (void **) &R, n*n, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Create upper-triangular matrix R [n*n] corresponding to VS
  ierr = xf_Error(xf_VectorSetUpperHouse(VS, n, W, R));
  if (ierr != xf_OK) return ierr;

  // no support yet for variable orders
  if (VS->Vector[0].vOrder != NULL) return xf_Error(xf_NOT_SUPPORTED);

  /*  xf_printf("R = [\n"); */
  /*   for (j=0; j<n*n; j++){ */
  /*     xf_printf("%.10E ", R[j]); fflush(stdout); */
  /*     if ((j+1)%n == 0) xf_printf("\n"); */
  /*   } */

  // Call SVD calculation on R
  if (myRank == 0){
    ierr = xf_Error(xf_SVDGolubReinsch(R, n, n, xfe_True, W, NULL));
    if (ierr != xf_OK) return ierr;
    // note, R now contains the left singular vectors of the nxn system
  }
  
  if (ParallelFlag){ // Broadcast W, R if in parallel
    ierr = xf_Error(xf_MPI_Bcast((void *) W, n*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_MPI_Bcast((void *) R, n*n*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
  }

  // Set VS = VS * R  [VS was altered by UpperHouse] 
  ierr = xf_Error(xf_Alloc( (void **) &EV, n, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  A = NULL;
  V = VS->Vector + 0;
  for (egrp=0; egrp<V->nArraySelf; egrp++){
    nr    = V->GenArray[egrp].r;
    // allocate temporary matrices (for speedup)
    ierr = xf_Error(xf_ReAlloc( (void **) &A, 2*nr*n, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    T = A + nr*n;

    for (elem=0; elem<V->GenArray[egrp].n; elem++){
      // set A = VS((elem,all r), :)
      for (j=0; j<n; j++) EV[j] = VS->Vector[j].GenArray[egrp].rValue[elem];
      for (j=0; j<n; j++)
	for (r=0; r<nr; r++)
	  A[r*n+j] = EV[j][r];
      // set T = A*R
      xf_MxM_Set(A, R, nr, n, n, T);
      // put T back into VS
      for (j=0; j<n; j++)
	for (r=0; r<nr; r++)
	  EV[j][r] = T[r*n+j];
    } // elem
  } // egrp


  // Release memory
  xf_Release( (void *) R);
  xf_Release( (void *) EV);
  xf_Release( (void *) A);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorSetPOD
int 
xf_VectorSetPOD(xf_VectorSet **pSnapSet, int nSnap, int nBasis,
		enum xfe_Verbosity Verbosity)
{
  int ierr, i, j, jmax;
  int *used;
  real *W, wmax;
  xf_Data *D;
  xf_Vector *VList, *Vtemp;
  xf_VectorSet *SnapSet;

  SnapSet = (*pSnapSet);
  
  // allocate memory for singular values
  ierr = xf_Error(xf_Alloc( (void **) &W, nSnap, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // compute SVD of SnapSet; left singular vectors returned in SnapSet
  if (Verbosity != xfe_VerbosityLow)
    xf_printf("  Calculating SVD ... "); fflush(stdout);
  ierr = xf_Error(xf_VectorSetSVD(SnapSet, nSnap, W)); // W = singular values
  if (ierr != xf_OK) return ierr;
  if (Verbosity != xfe_VerbosityLow) xf_printf("done.\n");
  
  if (Verbosity == xfe_VerbosityHigh){
    xf_printf("W = [\n");
    for (i=0; i<nSnap; i++)
      xf_printf("%.12E\n", W[i]);
    xf_printf("]\n");
  }

  // allocate vector list
  ierr = xf_Error(xf_Alloc( (void **) &VList, nSnap, sizeof(xf_Vector)));
  if (ierr != xf_OK) return ierr;


  // temporary memory for sorting
  ierr = xf_Error(xf_Alloc( (void **) &used, nSnap, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nSnap; i++) used[i] = 0;

  // sort out the first nBasis largest vectors
  for (i=0; i<nBasis; i++){
    wmax = -1.;
    jmax = -1;
    for (j=0; j<nSnap; j++)
      if ((!used[j]) && (W[j] > wmax)){
	jmax = j;
	wmax = W[j];
      }
    if (wmax < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
    used[jmax] = 1;
    VList[i] = SnapSet->Vector[jmax];
  } // i
  
  // tag remaining vectors along to prevent memory problems
  i = nBasis;
  for (j=0; j<nSnap; j++) 
    if (!used[j])
      VList[i++] = SnapSet->Vector[j];
  if (i != nSnap) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  // swap SnapSet->Vector and VList
  swap(SnapSet->Vector, VList, Vtemp);

  // set pointer to SnapSet
  pSnapSet = &SnapSet;


  /* Release memory */
  xf_Release( (void *) VList);
  xf_Release( (void *) W);
  xf_Release( (void *) used);


  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_OrthonormalizeBasisSetMass
static int 
xf_OrthonormalizeBasisSetMass(xf_All *All, int nBasis, xf_VectorSet *BasisSet)
{
  /* Does what the name suggests.  If All is NULL, Mass matrix will
     not be used */
  int ierr, i, j;
  real dp, rnorm;
  xf_Vector *Utemp;

  xf_printf(" Orthonormalizing basis vectors... ");

  // allocate memory
  if (All != NULL){
    ierr = xf_Error(xf_FindSimilarVector(All, BasisSet->Vector+0, "Utemp",
					 xfe_True, xfe_False, NULL, &Utemp, NULL));
    if (ierr != xf_OK) return ierr;
  }

  for (i=0; i<nBasis; i++){
    if (All == NULL)
      Utemp = BasisSet->Vector+i;
    else{
      //Utemp = Vectori
      ierr = xf_Error(xf_SetVector(BasisSet->Vector+i, xfe_Set, Utemp));
      if (ierr != xf_OK) return ierr;
      // Utemp = M*Utemp
      ierr = xf_Error(xf_MultMassMatrix(All, 1.0, Utemp));
      if (ierr != xf_OK) return ierr;
    }

    // dp = Vectori^T * M * Vectori
    ierr = xf_Error(xf_VectorDot(Utemp, BasisSet->Vector+i, &dp));
    if (ierr != xf_OK) return ierr;

    rnorm = sqrt(dp);

    //xf_printf("i = %d, rnorm = %.10E, dp = %.10E\n", i, rnorm, dp);

    // Vectori /= rnorm
    ierr = xf_Error(xf_VectorMult(BasisSet->Vector+i, 1.0/rnorm));
    if (ierr != xf_OK) return ierr;

    // Utemp = Vectori
    if (All == NULL)
      Utemp = BasisSet->Vector+i;
    else{
      ierr = xf_Error(xf_SetVector(BasisSet->Vector+i, xfe_Set, Utemp));
      if (ierr != xf_OK) return ierr;
      // Utemp = M*Utemp
      ierr = xf_Error(xf_MultMassMatrix(All, 1.0, Utemp));
      if (ierr != xf_OK) return ierr;
    }

    for (j=i+1; j<nBasis; j++){
      // dp = Vectori^T * M * Vectorj = Utemp^T * Vectorj
      ierr = xf_Error(xf_VectorDot(Utemp, BasisSet->Vector+j, &dp));
      if (ierr != xf_OK) return ierr;

      if (dp > 1e-10){ // The incoming vectors should be mostly orthogonal
	xf_printf("Warning, Vi and Vj do not appear orthogonal. Vi^T * M * Vj = %.10E\n", dp);
      }

      // Vectorj -= dp*Vectori
      ierr = xf_Error(xf_VectorMultSet(BasisSet->Vector+i, dp, xfe_Sub, BasisSet->Vector+j));
      if (ierr != xf_OK) return ierr;

    } // j

  } // i

  // destroy Utemp
  if (All != NULL){
    ierr = xf_Error(xf_DestroyVector(Utemp, xfe_True));
    if (ierr != xf_OK) return ierr;
  }    

  xf_printf(" done.\n");

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_VectorSetPODMass
int 
xf_VectorSetPODMass(xf_All *All, const xf_VectorSet *SnapSet, 
		    int nSnap, int nBasis, xf_VectorSet *BasisSet)
{
  // If All is NULL, Mass matrix will not be used
  int ierr, i, j;
  real dp;
  real *K, *EG;
  xf_Data *D;
  xf_Vector *Utemp;

  // allocate memory
  ierr = xf_Error(xf_Alloc( (void **) &K, nSnap*nSnap, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  if (All != NULL){
    ierr = xf_Error(xf_FindSimilarVector(All, SnapSet->Vector+0, "Utemp",
					 xfe_True, xfe_False, NULL, &Utemp, NULL));
    if (ierr != xf_OK) return ierr;
  }

  xf_printf("  Calculating K ... ");
  for (i=0; i<nSnap; i++){
    if (All == NULL)
      Utemp = SnapSet->Vector+i;
    else{
      // Utemp = Vector+i
      ierr = xf_Error(xf_SetVector(SnapSet->Vector+i, xfe_Set, Utemp));
      if (ierr != xf_OK) return ierr;
      // Utemp = M*Utemp
      ierr = xf_Error(xf_MultMassMatrix(All, 1.0, Utemp));
      if (ierr != xf_OK) return ierr;
    }

    for (j=i; j<nSnap; j++){
      ierr = xf_Error(xf_VectorDot(Utemp, SnapSet->Vector+j, &dp));
      if (ierr != xf_OK) return ierr;
      dp = dp/( (real) nSnap);
      K[i*nSnap+j] = dp;
      K[j*nSnap+i] = dp; // K is symmetric
    } // j
  } // i
  xf_printf("done.\n");

  // destroy Utemp
  if (All != NULL){
    ierr = xf_Error(xf_DestroyVector(Utemp, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  
  /*  xf_printf("K = ["); */
  /*   for (i=0; i<nSnap; i++) */
  /*     for (j=0, xf_printf("\n"); j<nSnap; j++) */
  /*       xf_printf("%.6E ", K[i*nSnap+j]); */
  /*   xf_printf("\n]\n"); */


  // allocate memory for eigenvalues and eigenvectors
  ierr = xf_Error(xf_Alloc( (void **) &EG, nSnap, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  /* Calculate nBasis max eigs and eigvs of K.*/
  xf_printf("  Calculating eigenvalues ... "); fflush(stdout);
  ierr = xf_Error(xf_EigSym(nSnap, K, EG));  // V^T overwrites K
  if (ierr != xf_OK) return ierr;  // note, eigs are increasing order
  for (i=0; i<nSnap; i++) xf_printf("\n %.10E", EG[i]);
  xf_printf("\ndone.\n");

  
  /* Construct basis vectors, using eigvectors corresponding to
     largest eigs (last ones) first */
  xf_printf("  Constructing Basis vectors ... ");
  for (i=0; i<nBasis; i++){
    ierr = xf_Error(xf_SetZeroVector(BasisSet->Vector+i));
    if (ierr != xf_OK) return ierr;
    for (j=0; j<nSnap; j++){
      ierr = xf_Error(xf_VectorMultSet(SnapSet->Vector+j, K[(nSnap-1-i)*nSnap+j], 
				       xfe_Add, BasisSet->Vector+i));
      if (ierr != xf_OK) return ierr;
    }
  }
  xf_printf("done.\n");

  
  // Orthonormalize BasisSet
  ierr = xf_Error(xf_OrthonormalizeBasisSetMass(All, nBasis, BasisSet));
  if (ierr != xf_OK) return ierr;


  /* Release memory */
  xf_Release( (void *) K);
  xf_Release( (void *) EG);


  return xf_OK;
}

/******************************************************************/
//   FUNCTION Prototype: xf_MatrixFrobNorm
int 
xf_MatrixFrobNorm(real *M, int n, int m, real *Fnorm)
{
  int i, j, pos;
  
  if (M == NULL) return xf_Error(xf_INPUT_ERROR);
  
  (*Fnorm) = 0.0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      pos = i*m + j;
      (*Fnorm) += pow(M[pos],2.0);
    }
  }
  
  (*Fnorm) = sqrt((*Fnorm));
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Prototype: xf_ConvertVectorFromLagrange
int 
xf_ConvertVectorFromLagrange(xf_Vector *U)
{
  int ierr, i;
  enum xfe_BasisType *BasisNew = NULL;

  // U must have a basis
  if (U->Basis == NULL) return xf_Error(xf_INPUT_ERROR);

  // allocate vector for storing  bases
  ierr = xf_Error(xf_Alloc((void **) &BasisNew, U->nArray, sizeof(enum xfe_BasisType)));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<U->nArray; i++){
    // get an appropriate Lagrange basis
    BasisNew[i] = U->Basis[i];
    ierr = xf_Error(xf_Basis2Lagrange(BasisNew[i], U->Basis + i));
    if (ierr != xf_OK) return ierr;
  }
  
  // convert from Lagrange basis to actual Basis
  ierr = xf_Error(xf_ProjectVectorInPlace_Basis(NULL, NULL, U, BasisNew, xfe_BasisLast));
  if (ierr != xf_OK) return ierr;
  
  xf_Release((void *) BasisNew);
  
  return xf_OK;
}

  
#if( UNIT_TEST==1 )
#include "xf_DataMath.test.in"
#endif