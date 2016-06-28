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
  FILE:  xf_LinQueue.c

  This file contains functions that assist residual linearization.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_Math.h"
#include "xf_LinQueueStruct.h"


/******************************************************************/
//   FUNCTION Definition: xf_InitLinQueue
int
xf_InitLinQueue(xf_LinQueueData *LinQ)
{
  int k;
  
  for (k=0; k<xfe_LinQTerm_Last; k++){
    LinQ->F_u[k]    = NULL;
    LinQ->Active[k] = xfe_False;
    LinQ->Size[k]   = 0;
  }
  LinQ->T     = NULL;
  LinQ->Tsize = 0;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyLinQueue
void
xf_DestroyLinQueue(xf_LinQueueData *LinQ)
{
  int k;
  
  for (k=0; k<xfe_LinQTerm_Last; k++)
    xf_Release( (void *) LinQ->F_u[k]);

  xf_Release( (void *) LinQ->T);
}


/******************************************************************/
//   FUNCTION Definition: xf_ClearLinQueue
void
xf_ClearLinQueue(xf_LinQueueData *LinQ)
{
  int k;
  
  for (k=0; k<xfe_LinQTerm_Last; k++)
    LinQ->Active[k] = xfe_False;

}

/******************************************************************/
//   FUNCTION Definition: xf_GetTermFromLinQueue
int
xf_GetTermFromLinQueue(xf_LinQueueData *LinQ, enum xfe_LinQTermType Term,
		       int nq, int dim, int sr2, real **pF_u)
{
  int ierr, k;
  int Fsize;

  // determine required size
  Fsize = nq*sr2;
  if (Term != xfe_LinQTerm_PhiPhi) Fsize *= dim;
  if (Term == xfe_LinQTerm_GPhiGPhi) Fsize *= dim;
  
  // reallocate term
  if (Fsize > LinQ->Size[Term]){
    if (LinQ->Active[Term]) return xf_Error(xf_OUT_OF_BOUNDS); // term should not be active
    ierr = xf_Error(xf_ReAlloc((void **) &LinQ->F_u[Term], Fsize, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    LinQ->Size[Term] = Fsize;
  }
  
  // zero out and activate if not active
  if (!LinQ->Active[Term]){
    for (k=0; k<Fsize; k++) LinQ->F_u[Term][k] = 0.;
    LinQ->Active[Term] = xfe_True;
  }

  // return term
  (*pF_u) = LinQ->F_u[Term];

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AddToLinQueue
int
xf_AddToLinQueue(real *F_u, enum xfe_LinQTermType Term,
		 int nq, int dim, int sr2, int idim, real *w,
		 int dw, real fac, xf_LinQueueData *LinQ)
{
  int ierr, iq, k;
  int Fsize;
  int TermDim;
  int off, ndim, iqsr2;
  real val;

  // dimension of F_u associated with term
  TermDim = (Term != xfe_LinQTerm_PhiPhi) ? dim : 1;
  if (Term == xfe_LinQTerm_GPhiGPhi) TermDim = dim*dim;

  // determine required size of matrix
  Fsize = nq*sr2*TermDim;
  
  // reallocate term
  if (Fsize > LinQ->Size[Term]){
    if (LinQ->Active[Term]) return xf_Error(xf_OUT_OF_BOUNDS); // term should not be active
    ierr = xf_Error(xf_ReAlloc((void **) &LinQ->F_u[Term], Fsize, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    LinQ->Size[Term] = Fsize;
  }
  
  // zero out and activate if not active
  if (!LinQ->Active[Term]){
    for (k=0; k<Fsize; k++) LinQ->F_u[Term][k] = 0.;
    LinQ->Active[Term] = xfe_True;
  }

  // determine offset and number of dimensions to add
  off  = (idim < 0) ? 0 : idim*nq*sr2; // offset for adding to LinQ->F_u
  ndim = (idim < 0) ? TermDim: 1;      // how many dimensions to consider

  // add F_u to LinQ->F_u
  for (iq=0; iq<nq*ndim; iq++){
    val = (w == NULL) ? fac : fac*w[(iq%nq)*dw];
    iqsr2 = iq*sr2;
    for (k=0; k<sr2; k++) LinQ->F_u[Term][off+iqsr2+k] += val*F_u[iqsr2+k];
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyLinQueue
int
xf_ApplyLinQueue(xf_LinQueueData *LinQ, xf_BasisData *RowPhiData,
		 xf_BasisData *ColPhiData, int sr2, real *R_U)
{
  int ierr;
  int nq, n, i, j;
  int nnRow, nnCol;
  int sizeReq, dim;
  enum xfe_LinQTermType Term;
  real *T;

  // number of basis functions
  nnRow = RowPhiData->nn;
  nnCol = ColPhiData->nn;

  // number of quad points and dimension
  nq  = RowPhiData->nq; 
  dim = RowPhiData->dim;

  // reallocate temporary storage if necessary
  sizeReq = nq*max(nnRow,nnCol);
  if (sizeReq > LinQ->Tsize){
    LinQ->Tsize = sizeReq;
    ierr = xf_Error(xf_ReAlloc( (void **) &LinQ->T, LinQ->Tsize, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  T = LinQ->T;
 
  for (n=0; n<nnRow; n++){
    // Phi^T  * F_u *  Phi
    Term = xfe_LinQTerm_PhiPhi;
    if (LinQ->Active[Term]){
      xf_ColMult_Set(ColPhiData->Phi, RowPhiData->Phi+n, nq, nnCol, nnRow, T);
      xf_MTxM_Add(T, LinQ->F_u[Term], nnCol, nq, sr2, R_U+n*nnCol*sr2);
    }
    // GPhi^T * F_u *  Phi
    Term = xfe_LinQTerm_GPhiPhi;
    if (LinQ->Active[Term]){
      for (i=0; i<dim; i++){
	xf_ColMult_Set(ColPhiData->Phi, RowPhiData->gPhi+nnRow*nq*i+n, nq, nnCol, nnRow, T);
	xf_MTxM_Add(T, LinQ->F_u[Term]+i*nq*sr2, nnCol, nq, sr2, R_U+n*nnCol*sr2);
      }
    }
    // Phi^T *  F_u * GPhi
    Term = xfe_LinQTerm_PhiGPhi;
    if (LinQ->Active[Term]){
      for (i=0; i<dim; i++){
	xf_ColMult_Set(ColPhiData->gPhi+nnCol*nq*i, RowPhiData->Phi+n, nq, nnCol, nnRow, T);
	xf_MTxM_Add(T, LinQ->F_u[Term]+i*nq*sr2, nnCol, nq, sr2, R_U+n*nnCol*sr2);
      }
    }
    // GPhi^T * F_u * GPhi
    Term = xfe_LinQTerm_GPhiGPhi;
    if (LinQ->Active[Term]){
      for (i=0; i<dim; i++){ // row
	for (j=0; j<dim; j++){ // col
	  xf_ColMult_Set(ColPhiData->gPhi+nnCol*nq*j, RowPhiData->gPhi+nnRow*nq*i+n, 
			 nq, nnCol, nnRow, T);
	  xf_MTxM_Add(T, LinQ->F_u[Term]+nq*sr2*(i*dim+j), nnCol, nq, sr2, R_U+n*nnCol*sr2);
	}
      }
    }
  } // n

  return xf_OK;
}





