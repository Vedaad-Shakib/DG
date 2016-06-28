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
  FILE:  xf_LeanSolver.c

  This file contains the memory-lean solvers.

*/

#include "xf_AllStruct.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_SolverTools.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Residual.h"
#include "xf_LinearSolver.h"
#include "xf_Line.h"
#include "xf_Basis.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Memory.h"
#include "xf_Quad.h"
#include "xf_EqnSetHook.h"
#include "xf_MeshTools.h"
#include "xf_Log.h"
#include "xf_Output.h"




/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_BlockJacobiLean
int
xf_Jacobian_SolveM_BlockJacobiLean(xf_All *All, xf_JacobianMatrix *R_U,
				   xf_Vector *X, enum xfe_AddType AddFlag,
				   enum xfe_Bool TransposeFlag, 
				   xf_SolverData *SolverData)
{  
  int ierr, i, k, sr, sr2, nn, nnmax;
  int egrp, elem, pOrder, Order;
  real c, fac, val, *MM;
  real *ER = NULL, *ER_EU = NULL;
  xf_Vector *U;
  xf_Vector *dt;
  xf_Matrix *M;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  // pull off state
  if ((U = R_U->U) == NULL) return xf_Error(xf_INPUT_ERROR);

  sr = U->StateRank;
  sr2 = sr*sr;

  // pull off real/artificial time step info
  c  = SolverData->c;  // real time step constant on M
  dt = SolverData->dt; // artificial time step vector

  nnmax = -1;

  // loop over elements, PLU solve diagonal blocks, apply to X
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    pOrder = -1;

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);

      // determine nn = number of basis functions per element
      nn = xf_Jacobian_n(R_U,egrp,elem);

      if (nn > nnmax){
	nnmax = nn;
	// allocate space for ER, and ER_EU
	ierr = xf_Error(xf_ReAlloc( (void **) &ER, nn*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;

	ierr = xf_Error(xf_ReAlloc( (void **) &ER_EU, nn*nn*sr*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      if (Order != pOrder){
	pOrder = Order;
	// find mass matrix data if necessary
	if ((c != 0.) || (dt != NULL)){
	  ierr = xf_Error(xf_FindMassMatrixData(All, egrp, U->Basis[egrp], Order, &M));
	  if (ierr != xf_OK) return ierr;
	}
      }

      // calculate residual and self Jacobian, ER_EU
      ierr = xf_Error(xf_CalculateResidualLeanElem(All, egrp, elem, U, ER, ER_EU, 
						   NULL, NULL, R_U, SolverData));
      if (ierr != xf_OK) return ierr;

      // add (c+1/dt)*M to Jacobian.  Note, dt = artificial time step
      if ((c != 0.) || (dt != NULL)){
	ierr = xf_Error(xf_ElemMassMatrix(All, egrp, elem, U->Basis[egrp], 
					  Order, M, NULL, &MM, &fac));
	if (ierr != xf_OK) return ierr;

	if (dt == NULL)
	  fac *= c;
	else
	  fac *= (c + 1.0/dt->GenArray[egrp].rValue[elem][0]);
	
	for (i=0; i<nn*nn; i++){
	  val = MM[i]*fac;
	  for (k=0; k<sr2; k+=(sr+1))
	    ER_EU[i*sr2+k] += val;
	} // i
      }

      // compute block PLU of ER_EU
      ierr = xf_Error(xf_ComputeBlockPLU(ER_EU, nn, R_U->StateRank, R_U->P[egrp][elem]));
      if (ierr != xf_OK) return ierr;


      if (!TransposeFlag){ // X(elem) @= (ER_EU)^{-1} * X(elem)
	ierr = xf_Error(xf_SolveBlockPLU(ER_EU, nn,
					 R_U->StateRank, R_U->P[egrp][elem],
					 X->GenArray[egrp].rValue[elem], AddFlag,
					 X->GenArray[egrp].rValue[elem]));
	if (ierr != xf_OK) return ierr;
      }
      else{ // X(elem) @= (ER_EU)^{-T} * X(elem)
	ierr = xf_Error(xf_SolveBlockPLUT(ER_EU, nn,
					  R_U->StateRank, R_U->P[egrp][elem],
					  X->GenArray[egrp].rValue[elem], AddFlag,
					  X->GenArray[egrp].rValue[elem]));
	if (ierr != xf_OK) return ierr;
      }
	
    } // elem
  } // egrp

  xf_Release( (void *) ER);
  xf_Release( (void *) ER_EU);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultM_BlockJacobiLean
int
xf_Jacobian_MultM_BlockJacobiLean(xf_All *All, xf_JacobianMatrix *R_U,
				  xf_Vector *X, enum xfe_AddType AddFlag,
				  enum xfe_Bool TransposeFlag, 
				  xf_SolverData *SolverData, xf_Vector *Y)
{  
  int ierr, i, k, sr, sr2, nn, nnmax;
  int egrp, elem, pOrder, Order;
  real c, fac, val, *MM;
  real *ER = NULL, *ER_EU = NULL;
  xf_Vector *U;
  xf_Vector *dt;
  xf_Matrix *M;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  // pull off state
  if ((U = R_U->U) == NULL) return xf_Error(xf_INPUT_ERROR);

  sr = U->StateRank;
  sr2 = sr*sr;

  // pull off real/artificial time step info
  c  = SolverData->c;  // real time step constant on M
  dt = SolverData->dt; // artificial time step vector

  nnmax = -1;

  // loop over elements, PLU solve diagonal blocks, apply to X
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    pOrder = -1;

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);

      // determine nn = number of basis functions per element
      nn = xf_Jacobian_n(R_U,egrp,elem);

      if (nn > nnmax){
	nnmax = nn;
	// allocate space for ER, and ER_EU
	ierr = xf_Error(xf_ReAlloc( (void **) &ER, nn*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_ReAlloc( (void **) &ER_EU, nn*nn*sr*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }
      
      if (Order != pOrder){
	pOrder = Order;
	// find mass matrix data if necessary
	if ((c != 0.) || (dt != NULL)){
	  ierr = xf_Error(xf_FindMassMatrixData(All, egrp, U->Basis[egrp], Order, &M));
	  if (ierr != xf_OK) return ierr;
	}
      }

      // calculate residual
      ierr = xf_Error(xf_CalculateResidualLeanElem(All, egrp, elem, U, ER, ER_EU, 
						   NULL, NULL, R_U, SolverData));
      if (ierr != xf_OK) return ierr;

      // add (c+1/dt)*M to Jacobian.  Note, dt = artificial time step
      if ((c != 0.) || (dt != NULL)){
	ierr = xf_Error(xf_ElemMassMatrix(All, egrp, elem, U->Basis[egrp], 
					  Order, M, NULL, &MM, &fac));
	if (ierr != xf_OK) return ierr;

	if (dt == NULL)
	  fac *= c;
	else
	  fac *= (c + 1.0/dt->GenArray[egrp].rValue[elem][0]);
	
	for (i=0; i<nn*nn; i++){
	  val = MM[i]*fac;
	  for (k=0; k<sr2; k+=(sr+1))
	    ER_EU[i*sr2+k] += val;
	} // i
      }

      if (!TransposeFlag){  // Y(elem) @= (ER_EU) * X(elem)
	ierr = xf_Error(xf_BlockMxV(ER_EU, nn, R_U->StateRank, nn,
				    X->GenArray[egrp].rValue[elem], AddFlag,
				    Y->GenArray[egrp].rValue[elem]));
	if (ierr != xf_OK) return ierr;
      }
      else{ // Y(elem) @= (ER_EU)^T * X(elem)
	ierr = xf_Error(xf_BlockMTxV(ER_EU, nn, R_U->StateRank, nn,
				     X->GenArray[egrp].rValue[elem], AddFlag,
				     Y->GenArray[egrp].rValue[elem]));
	if (ierr != xf_OK) return ierr;
      }
	
    } // elem
  } // egrp

  xf_Release( (void *) ER);
  xf_Release( (void *) ER_EU);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultN_BlockJacobiLean
int
xf_Jacobian_MultN_BlockJacobiLean(xf_All *All, xf_JacobianMatrix *R_U,
				  xf_Vector *X, enum xfe_AddType AddFlag,
				  enum xfe_Bool TransposeFlag, 
				  xf_SolverData *SolverData,
				  enum xfe_Bool HaloFlag, xf_Vector *Y)
{ 
  int ierr, k, sr, nn, r, nnN, nnmax;
  int egrp, elem, face, negrp;
  int egrpN, elemN, faceN;
  int nface, pnface, egN;
  int rrN, *rrNvec = NULL;
  enum xfe_Bool ReAllocFlag, NearHalo;
  enum xfe_AddType AddFlag2; 
  real *rY;
  real *ER = NULL;
  real **ER_NU = NULL, **NR_EU = NULL;
  xf_Vector *U;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  // pull off state
  if ((U = R_U->U) == NULL) return xf_Error(xf_INPUT_ERROR);

  sr = U->StateRank;
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);

  pnface = -1;
  nnmax = -1;
  
  // loop over elements, PLU solve diagonal blocks, apply to X
  negrp = Mesh->nElemGroup;
  for (egrp=0; egrp<negrp; egrp++){
  
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // determine nn = number of basis functions per element
      nn = xf_Jacobian_n(R_U,egrp,elem);
      r  = sr*nn;
      
      if (nn > nnmax){
	nnmax = nn;
	// allocate space for ER
	ierr = xf_Error(xf_ReAlloc( (void **) &ER, nn*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      // number of faces for this element
      nface = Mesh->ElemGroup[egrp].nFace[elem];
      
      // is elem adjacent to halo?
      NearHalo = xfe_False;
      for (face=0; face<nface; face++)
	NearHalo = (NearHalo || (R_U->egrpN[egrp][elem][face] >= negrp));

      // continue if only doing halo and not near a halo elem
      if ((!NearHalo) && HaloFlag) continue;

      // do we need to reallocate Jacobians?
      ReAllocFlag = xfe_False;
      if (pnface < nface){
	ierr = xf_Error(xf_ReAlloc( (void **) &rrNvec, nface, sizeof(int)));
	if (ierr != xf_OK) return ierr;
	ReAllocFlag = xfe_True;
	for (face=max(pnface,0); face<nface; face++) rrNvec[face] = 0;
	pnface = nface;
      }

      // check sizes of neighbor blocks
      for (face=0; face<nface; face++){
	egrpN = R_U->egrpN[egrp][elem][face];
	elemN = R_U->elemN[egrp][elem][face];
	nnN = xf_Jacobian_n(R_U,egrpN,elemN);
	rrN = ((egrpN < 0) ? 0 : nnN*sr*r);
	if (rrNvec[face] < rrN){
	  ReAllocFlag = xfe_True;
	  rrNvec[face] = rrN;
	}
      }
	
      // rrNvec[face] = # size of rectangular off-diagonal Jacobian
      if (ReAllocFlag){
	ierr = xf_Error(xf_VReAlloc2( (void ***) &ER_NU, pnface, rrNvec, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VReAlloc2( (void ***) &NR_EU, pnface, rrNvec, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }


      // calculate residual
      ierr = xf_Error(xf_CalculateResidualLeanElem(All, egrp, elem, U, ER, NULL, 
						   ER_NU, NR_EU, R_U, SolverData));
      if (ierr != xf_OK) return ierr;


      if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){ // Set to 0 if necessary
	rY = Y->GenArray[egrp].rValue[elem];
	for (k=0; k<nn*sr; k++) rY[k] = 0.0;
      }
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
	egrpN = R_U->egrpN[egrp][elem][face];
	elemN = R_U->elemN[egrp][elem][face];
	nnN = xf_Jacobian_n(R_U,egrpN,elemN);

	if (  HaloFlag  &&  (egrpN < negrp)) continue; // only considering halo
	if ((!HaloFlag) && ((egrpN < 0) || (egrpN >= negrp))) continue; // X on halo is not here yet

	if (!TransposeFlag){ // Y(elem) @= (ER_NU) * X(elemN)
	  ierr = xf_Error(xf_BlockMxV(ER_NU[face], nn, R_U->StateRank, 
				      nnN, X->GenArray[egrpN].rValue[elemN], 
				      AddFlag2, Y->GenArray[egrp].rValue[elem]));
	  if (ierr != xf_OK) return ierr;
	}
	else{ // Y(elem) @= (NR_EU)^T * X(elemN)
	  faceN = R_U->faceN[egrp][elem][face];
	  ierr = xf_Error(xf_BlockMTxV(NR_EU[face], nn, R_U->StateRank, 
				       nnN, X->GenArray[egrpN].rValue[elemN], 
				       AddFlag2, Y->GenArray[egrp].rValue[elem]));
	  if (ierr != xf_OK) return ierr;
	}
      } // face
    } // elem
  } // egrp

  xf_Release(  (void  *) rrNvec);
  xf_Release(  (void  *) ER);
  xf_Release2( (void **) ER_NU);
  xf_Release2( (void **) NR_EU);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_LeanSolverStepBlockElem
int 
xf_LeanSolverStepBlockElem( xf_All *All, int egrp, int elem, 
			    xf_SolverData *SolverData, xf_Vector *U, 
			    xf_JacobianMatrix *R_U, xf_Vector *dU)
{
  int ierr, i, k, sr, sr2, nn, Order;
  int *P = NULL;
  real c, fac, val, *MM;
  real *ER = NULL, *ER_EU = NULL;
  xf_Mesh *Mesh;
  xf_Vector *dt;

  Mesh = All->Mesh;
  sr = U->StateRank;
  sr2 = sr*sr;

  // determine nn = number of basis functions per element
  nn = xf_Jacobian_n(R_U,egrp,elem);
  
  // get interpolation order
  Order = xf_InterpOrder(U, egrp, elem);
  
  // allocate space for P, ER, and ER_EU
  ierr = xf_Error(xf_ReAlloc( (void **) &P, nn*sr, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_ReAlloc( (void **) &ER, nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_ReAlloc( (void **) &ER_EU, nn*nn*sr*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // do not need lines
  if (SolverData != NULL) SolverData->CRequired = xfe_False;

    
  ierr = xf_Error(xf_CalculateResidualLeanElem(All, egrp, elem, U, ER, ER_EU, 
					       NULL, NULL, R_U, SolverData));
  if (ierr != xf_OK) return ierr;
  
  // add (c+1/dt)*M to Jacobian.  Note, dt = artificial time step
  c  = SolverData->c;  // real time step constant on M
  dt = SolverData->dt; // artificial time step vector
  if ((c != 0.) || (dt != NULL)){
    ierr = xf_Error(xf_ElemMassMatrix(All, egrp, elem, U->Basis[egrp], 
				      Order, NULL, NULL, &MM, &fac));
    if (ierr != xf_OK) return ierr;
    
    if (dt == NULL)
      fac *= c;
    else
      fac *= (c + 1.0/dt->GenArray[egrp].rValue[elem][0]);
    
    for (i=0; i<nn*nn; i++){
      val = MM[i]*fac;
      for (k=0; k<sr2; k+=(sr+1))
	ER_EU[i*sr2+k] += val;
    } // i
  }


  /* dU(elem) = - (ER_EU)^{-1} * ER   (carefully as ER_EU is in block format) */
  
  ierr = xf_Error(xf_ComputeBlockPLU(ER_EU, nn, sr, P));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SolveBlockPLU(ER_EU, nn, sr, P, ER, xfe_Neg, 
				   dU->GenArray[egrp].rValue[elem]));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) P);
  xf_Release( (void *) ER);
  xf_Release( (void *) ER_EU);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_LeanSolverStepBlock
static int 
xf_LeanSolverStepBlock( xf_All *All, xf_SolverData *SolverData, 
			xf_Vector *U, xf_JacobianMatrix *R_U, xf_Vector *dU)
{
/*
PURPOSE:

  Applies one step of a memory-lean block solver to obtain a state
  update, dU.

INPUTS:

  All : All structure
  SolverData : structure of useful solver data + params
  U : state vector
  R_U : Jacobian matrix without values, but with connectivity

OUTPUTS: 

  dU : state update

RETURN:

  Error Code
*/
  int ierr, k, sr, nn, nnmax;
  int egrp, elem, *P = NULL;
  real *ER = NULL, *ER_EU = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  sr = U->StateRank;
  
  nnmax = -1;

  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      // determine nn = number of basis functions per element
      nn = R_U->nvec[egrp];
      
      if (nn > nnmax){
	nnmax = nn;
	// allocate space for P, ER, and ER_EU
	ierr = xf_Error(xf_ReAlloc( (void **) &P, nn*sr, sizeof(int)));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_ReAlloc( (void **) &ER, nn*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_ReAlloc( (void **) &ER_EU, nn*nn*sr*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }	

      ierr = xf_Error(xf_CalculateResidualLeanElem(All, egrp, elem, U, ER, ER_EU, 
						   NULL, NULL, R_U, SolverData));
      if (ierr != xf_OK) return ierr;

      /* dU(elem) = - (ER_EU)^{-1} * ER   (carefully as ER_EU is in block format) */

      ierr = xf_Error(xf_ComputeBlockPLU(ER_EU, nn, sr, P));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_SolveBlockPLU(ER_EU, nn, sr, P, ER, xfe_Neg, 
				       dU->GenArray[egrp].rValue[elem]));
      if (ierr != xf_OK) return ierr;

    } // elem

  } // egrp

  xf_Release( (void *) P);
  xf_Release( (void *) ER);
  xf_Release( (void *) ER_EU);

  return xf_OK;
}


/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_LeanSolverStep */
/* int  */
/* xf_LeanSolverStep( xf_All *All, enum xfe_NonlinearSolverType NonlinearSolver,  */
/* 		   xf_SolverData *SolverData, xf_Vector *U, */
/* 		   xf_JacobianMatrix *R_U, xf_Vector *dU) */
/* { */
/*   int ierr; */
  
/*   switch(NonlinearSolver){ */
/*   case xfe_NonlinearSolverLeanBlock:  */
/*     ierr = xf_Error(xf_LeanSolverStepBlock(All, SolverData, U, R_U, dU)); */
/*     if (ierr != xf_OK) return ierr; */
/*     break; */
/*   case xfe_NonlinearSolverNewton:  */
/*   case xfe_NonlinearSolverNone: */
/*   default: */
/*     return xf_Error(xf_INPUT_ERROR); */
/*     break; */
/*   } */

/*   return xf_OK; */
/* } */
