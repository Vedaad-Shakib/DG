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
 FILE:  xf_Penalty.c
 
 This file contains interior penalty functions.
 
 */

#include "xf.h"
#include "xf_AllStruct.h"
#include "xf_Solver.h"
#include "xf_SolverStruct.h"
#include "xf_Residual.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Basis.h"
#include "xf_Math.h"
#include "xf_Quad.h"
#include "xf_Memory.h"
#include "xf_Penalty.h"
#include "xf_LinearSolver.h"
#include "xf_EqnSetHook.h"
#include "xf_MeshTools.h"
#include "xf_MPI.h"

/******************************************************************/
//   FUNCTION Definition: xf_FindPenaltyHessian
extern int
xf_FindPenaltyHessian(xf_All *All, int egrp, xf_Vector *U, xf_Matrix **pHp)
{
  int ierr, Order, nElem, Hrank;
  enum xfe_BasisType Basis;
  enum xfe_Bool found;
  
  Order = U->Order[egrp];
  Basis = U->Basis[egrp];
  
  nElem = All->Mesh->ElemGroup[egrp].nElem;
  Hrank = U->GenArray[egrp].r;//Hrank = nn*sr
  
  ierr = xf_Error(xf_FindMatrix(All->DataSet, "PenaltyHessian", xfe_LinkageElem, egrp,
                                Order, Order, Basis, Basis, xfe_SizeReal,
                                nElem, Hrank*Hrank, xfe_True, NULL, &(*pHp), &found));
  if (ierr != xf_OK) return ierr;
  
  //Permutation vector
  if (!found){
    ierr = xf_Error(xf_Alloc2((void ***) &((*pHp)->P), nElem, Hrank, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculatePenalty
int
xf_CalculatePenalty(xf_All *All, xf_Vector *U, xf_Vector *P, real *pP)
{
  int ierr, egrp, elem, Order, iq, nq, sr, pnq, nn, *IParam, QuadOrder;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged;
  real *xq, *Pq, *u, *EU, *RParam, Pelem;
  xf_BasisData *PhiData;
  xf_QuadData *QuadData;
  
  PhiData = NULL;
  QuadData = NULL;
  RParam = NULL;
  xq = NULL;
  u = NULL;
  sr = All->EqnSet->StateRank;
  Pq = NULL;
  pnq = -1;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  if (pP != NULL)
    (*pP) = 0.0;
  
  for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++){
    Basis = U->Basis[egrp];

    for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
      Order = xf_InterpOrder(U, egrp, elem);
      // determine required integration order
      ierr = xf_Error(xf_GetQuadOrderElem(All->Mesh, All->EqnSet, egrp, Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      /* Pull off quad points for the element; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(All->Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions (and grads) if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **)  &u,  nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ReAlloc( (void **)  &Pq,  nq, sizeof(real)));
        if (ierr != xf_OK) return ierr;
      }
      
      nn = PhiData->nn;
      EU  =  U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
      // interpolate state at quad points
      xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u);
      
      //computing the penalty term at the element
      ierr = xf_Error(xf_EqnSetPenaltyTerm(All->EqnSet, RParam, IParam, nq, u, NULL, Pq, NULL, NULL));
      if (ierr != xf_OK) return ierr;
      
      Pelem = 0.0;
      for (iq = 0; iq < nq; iq++) {
        Pelem += QuadData->wquad[iq]*Pq[iq];
      }//iq
      if (pP != NULL){
        (*pP) = max(Pelem,(*pP));
      }
      P->GenArray[egrp].rValue[elem][0] = Pelem;
      
      pnq = nq;
    }//elem
  }//egrp
  
  if (pP != NULL){
    //sum among all processors
    ierr = xf_Error(xf_MPI_Allreduce(&(*pP), 1, xfe_SizeReal, xfe_MPI_MAX));
    if (ierr != xf_OK) return ierr;
  }
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  // release memory
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) u);
  xf_Release( (void *) Pq);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculatePenaltyGradient
int
xf_CalculatePenaltyGradient(xf_All *All, xf_Vector *U, xf_Vector *GP)
{
  int ierr, egrp, elem, Order, nq, sr, pnq, q, i, j, nn, *IParam, QuadOrder;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged;
  real *xq, *u, *Pe_u, *EU, *RParam;
  xf_BasisData *PhiData;
  xf_QuadData *QuadData;
  
  PhiData = NULL;
  QuadData = NULL;
  xq = NULL;
  u = Pe_u = RParam = NULL;
  sr = All->EqnSet->StateRank;
  pnq = -1;
  
  //Zeroing the penalty gradient
  ierr = xf_Error(xf_SetZeroVector(GP));
  if (ierr != xf_OK) return ierr;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++){
    Basis = U->Basis[egrp];
    
    for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
      Order = xf_InterpOrder(U, egrp, elem);
      // determine required integration order
      ierr = xf_Error(xf_GetQuadOrderElem(All->Mesh, All->EqnSet, egrp, Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      /* Pull off quad points for the element; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(All->Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions (and grads) if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **)  &u,  nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ReAlloc( (void **)  &Pe_u,  nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        
        pnq = nq;
      }
      
      nn = PhiData->nn;
      EU  =  U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
      //interpolating the state at the quadrature points
      xf_MxM_Set(PhiData->Phi, EU,  nq, nn, sr, u);
      
      //computing the penalty gradient at each element
      ierr = xf_Error(xf_EqnSetPenaltyTerm(All->EqnSet, RParam, IParam, nq, u, NULL, NULL, Pe_u, NULL));
      if (ierr != xf_OK) return ierr;
      
      for (j = 0; j < nn; j++){//loop over basis functions
        for (i = 0; i < sr; i++){
          for (q = 0; q < nq; q++){
            GP->GenArray[egrp].rValue[elem][j*sr+i] += QuadData->wquad[q]*PhiData->Phi[q*nn+j]*Pe_u[q*sr+i];
          }//nq
        }//sr
      }//nn
      
    }//elem
  }//egrp
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  // release memory
  xf_Release( (void *) RParam);
  xf_Release( (void *) IParam);
  xf_Release( (void *) u);
  xf_Release( (void *) Pe_u);
  
  return xf_OK;
}

/****************************************************************************/
//   FUNCTION Definition: xf_CalculatePenaltyHessian
int
xf_CalculatePenaltyHessian(xf_All *All, xf_Vector *U, xf_Matrix *Hp, int egrp, int elem)
{
  int ierr, nq, nn, i, j, si, sj, q, sr, *IParam;
  int shift_i, shift_j, shift_si, shift_sj, shift_q, pos, Order, pnq, QuadOrder;
  real *Pe_uu, *Pe_UU, *xq, *u, *EU, Phi_i_at_q, Phi_j_at_q, *RParam;
  real ElemVol;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged;
  xf_BasisData *PhiData;
  xf_QuadData *QuadData;
  xf_Mesh *Mesh;
  
  PhiData = NULL;
  QuadData = NULL;
  xq = NULL;
  u = Pe_uu = Pe_UU = RParam = NULL;
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;
  pnq = nn = -1;
  
  Basis = U->Basis[egrp];
  Order = xf_InterpOrder(U, egrp, elem);
  // determine required integration order
  ierr = xf_Error(xf_GetQuadOrderElem(All->Mesh, All->EqnSet, egrp, Order, &QuadOrder));
  if (ierr != xf_OK) return ierr;
  
  /* Pull off quad points for the element; will not recalculate if
   Basis/Order have not changed. */
  ierr = xf_Error(xf_QuadElem(All->Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
  if (ierr != xf_OK) return ierr;
  
  nq = QuadData->nquad;
  xq = QuadData->xquad;
  
  // compute basis functions (and grads) if quad or basis or order changed
  ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
  if (ierr != xf_OK) return ierr;
  
  nn = PhiData->nn;
  
  // re-allocate data if quad points increased
  if (nq > pnq){
    ierr = xf_Error(xf_ReAlloc((void **)&u, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReAlloc((void **)&Pe_uu, nq*sr*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    pnq = nq;
  }
  
  //pointing to penalty Hessian
  Pe_UU = Hp->GenArray[0].rValue[elem];
  
  for (i = 0; i < nn*sr*nn*sr; i++)
    Pe_UU[i] = 0.0;
  
  EU  =  U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
  //interpolating the state at the limit points
  xf_MxM_Set(PhiData->Phi, EU,  nq, nn, sr, u);
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  //calculate penalty function Hessian
  ierr = xf_Error(xf_EqnSetPenaltyTerm(All->EqnSet, RParam, IParam, nq, u, NULL, NULL, NULL, Pe_uu));
  if (ierr != xf_OK) return ierr;
  
  //please look at documentation for reference about this calculation
  //this is a tensor product
  for (i = 0; i < nn; i++){//first basis function
    shift_i = i*nn*sr*sr;
    for (j = 0; j < nn; j++){//second basis function
      shift_j = j*sr*sr;
      for (si = 0; si < sr; si++){//first state variable
        shift_si = si*sr;
        for (sj = 0; sj < sr; sj++){//second state variable
          shift_sj = sj;
          pos = shift_i + shift_j + shift_si + shift_sj;
          Pe_UU[pos] = 0.0;
          for (q = 0; q < nq; q++){//quad points
            shift_q = q*sr*sr;
            Phi_i_at_q = PhiData->Phi[q*nn+i];
            Phi_j_at_q = PhiData->Phi[q*nn+j];
            Pe_UU[pos] += QuadData->wquad[q]*Phi_i_at_q*Phi_j_at_q*
                          Pe_uu[shift_q+shift_si+shift_sj];
          }//q
        }//sj
      }//si
    }//j
  }//i
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  // release memory
  Pe_UU = NULL;
  xf_Release( (void *) u);
  xf_Release( (void *) Pe_uu);
  xf_Release( (void *) RParam);
  xf_Release( (void *) IParam);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyPenaltyHessian
int
xf_ApplyPenaltyHessian(xf_All *All, xf_Vector *U, xf_Vector *Y, 
                       xf_Vector *W, enum xfe_Bool CalcHessian, 
                       enum xfe_Bool InverseFlag, enum xfe_AddType AddFlag)
{
  int ierr, egrp, elem, sr, nn;
  xf_Matrix *Hp;
  
  sr = All->EqnSet->StateRank;
  
  /* printf("CalcHessian: %d InverseFlag %d AddFlag: %d\n",CalcHessian,InverseFlag,AddFlag); */
  
  for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++){
    // find/allocate
    ierr = xf_Error(xf_FindPenaltyHessian(All, egrp, U, &Hp));
    if (ierr != xf_OK) return ierr;
    
    nn = U->GenArray[egrp].r/sr;
    
    for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
      if (CalcHessian){
        //the Hessian is stored in PLU form
        ierr = xf_Error(xf_CalculatePenaltyHessian(All, U, Hp, egrp, elem));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ComputeBlockPLU(Hp->GenArray[0].rValue[elem], nn, sr, Hp->P[elem]));
        if (ierr != xf_OK) return ierr;
      }
      
      if (InverseFlag){
        //W @= Hp^-1*Y where @ is AddFlag
        ierr = xf_Error(xf_SolveBlockPLU(Hp->GenArray[0].rValue[elem], nn, sr, Hp->P[elem], 
                                         Y->GenArray[egrp].rValue[elem], AddFlag, 
                                         W->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr;
      }
      else {
        //W @= Hp*Y where @ is AddFlag
        ierr = xf_Error(xf_BlockPLUMxV(Hp->GenArray[0].rValue[elem], nn, sr, Hp->P[elem], 
                                               Y->GenArray[egrp].rValue[elem], AddFlag, 
                                               W->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr; 
        
      }
    }
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyResidualHessian
int
xf_ApplyResidualHessian(xf_All *All, xf_Vector *U, xf_Vector *Y, 
                        xf_Vector *W, xf_SolverData *SolverData, 
                        enum xfe_AddType AddFlag)
{
  int ierr;
  xf_JacobianMatrix *R_U;
  xf_Vector *Y1, *Y2;
  enum xfe_Bool found;
  
  //Set pointer to Jacobian Matrix
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_False, NULL, &R_U, &found));
  if (ierr != xf_OK) return ierr;
  
  if (!found){
    xf_printf("Jacobian Matrix should exist at this point.\n");
    return xf_NOT_FOUND;
  }
  
  //Create the temporary vector Y
  ierr = xf_Error(xf_FindSimilarVector(All, U, "TempVector1", xfe_False, 
                                       xfe_False, NULL, &Y1, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "TempVector2", xfe_False, 
                                       xfe_False, NULL, &Y2, NULL));
  if (ierr != xf_OK) return ierr;
  
  //applying the Hessian without the penalty
  //First step: Y1 = dR_dU*Y
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Y, xfe_Set, xfe_False, SolverData, Y1));
  if (ierr != xf_OK) return ierr;
  
  //Second step: Y2 = 2*(dR_dU)^T*Y1
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Y1, xfe_Set, xfe_True, SolverData, Y2));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VectorMult(Y2,2.0));
  if (ierr != xf_OK) return ierr;
  
  //Set W @= Y2
  ierr = xf_Error(xf_SetVector(Y2, AddFlag, W));
  if (ierr != xf_OK) return ierr;
  
  //cleaning-up
  ierr = xf_Error(xf_DestroyVector(Y1,xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(Y2,xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyHessian
int
xf_ApplyHessian(xf_All *All, xf_Vector *U, xf_Vector *Y, xf_Vector *W,
                xf_SolverData *SolverData, enum xfe_Bool CalcHessian)
{
  //This function computes:
  //W = 2*(dR_dU)^T*dR_dU*Y + d2Pe_dUidUj*Y
  int ierr;
 
  
  //W = 2*(dR_dU)^T*dR_dU*Y
  ierr = xf_Error(xf_ApplyResidualHessian(All, U, Y, W, SolverData, xfe_Set));
  if (ierr != xf_OK) return ierr;
    
  //applying the penalty
  //W @= d2Pe_dUidUj*Y
  ierr = xf_Error(xf_ApplyPenaltyHessian(All, U, Y, W, CalcHessian, xfe_False, xfe_Add));
  if (ierr != xf_OK) return ierr;

  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyHessianPrecond
int
xf_ApplyHessianPrecond(xf_All *All, xf_Vector *U, xf_Vector *Y, xf_Vector *W,
                       xf_SolverData *SolverData, enum xfe_Bool CalcHessian)
{
  //This function computes:
  //W = (I - Hp^-1*Hr + (Hp^-1*Hr)*(Hp^-1*Hr))*Hp^-1*Y
  int ierr;
  xf_JacobianMatrix *R_U;
  xf_Vector *V, *Temp1, *Temp2;
  
  //Create the working vectors
  ierr = xf_Error(xf_FindSimilarVector(All, U, "VectorV", xfe_False, 
                                       xfe_False, NULL, &V, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "TempVector1", xfe_False, 
                                       xfe_False, NULL, &Temp1, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "TempVector2", xfe_False, 
                                       xfe_False, NULL, &Temp2, NULL));
  if (ierr != xf_OK) return ierr;
  
  //Set pointer to Jacobian Matrix
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_True, NULL, &R_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  //V = Hp^-1*Y
  ierr = xf_Error(xf_ApplyPenaltyHessian(All, U, Y, V, CalcHessian, xfe_True, xfe_Set));
  if (ierr != xf_OK) return ierr;
  
  //W = V
  ierr = xf_Error(xf_SetVector(V, xfe_Set, W));
  if (ierr != xf_OK) return ierr;
  
  //Temp1 = Hr*V (Temp1 = Hr*Hp^-1*Y)
  ierr = xf_Error(xf_ApplyResidualHessian(All, U, V, Temp1, SolverData, xfe_Set));
  if (ierr != xf_OK) return ierr;
  
  //Temp2 =  Hp^-1*Temp1
  ierr = xf_Error(xf_ApplyPenaltyHessian(All, U, Temp1, Temp2, CalcHessian, xfe_True, xfe_Set));
  if (ierr != xf_OK) return ierr;
  
  //W -= Temp2
  ierr = xf_Error(xf_SetVector(Temp2, xfe_Sub, W));
  if (ierr != xf_OK) return ierr;
  
  //Temp1 = Hr*Temp2
  ierr = xf_Error(xf_ApplyResidualHessian(All, U, Temp2, Temp1, SolverData, xfe_Set));
  if (ierr != xf_OK) return ierr;
  
  //W += Hp^-1*Temp1
  ierr = xf_Error(xf_ApplyPenaltyHessian(All, U, Temp1, W, CalcHessian, xfe_True, xfe_Add));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(V,xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(Temp1,xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(Temp2,xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_PenalizeResidual
int
xf_PenalizeResidual(xf_All *All, xf_Vector *U, xf_Vector *R, 
                    xf_SolverData *SolverData)
{
  int ierr, terr, egrp, elem, i, ntot;
  xf_Vector *P;
  
  if (SolverData == NULL) return xf_Error(xf_INPUT_ERROR);
  
  ierr = xf_CalculateResidual(All, U, R, NULL, SolverData);
  terr = xf_Error(xf_MPI_Allreduce(&ierr, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (terr != xf_OK) return terr;
  if (ierr != xf_OK) return ierr;
  
  if (SolverData->PenalizeResidual){
    P = SolverData->Pvec;//set the pointer to the penalty vector
    
    ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
    if (ierr != xf_OK) return ierr;
    
    for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++) {
      for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++) {
        if (U->GenArray[egrp].vr == NULL)
          ntot = U->GenArray[egrp].r;
        else 
          ntot = U->GenArray[egrp].vr[elem];
        for (i = 0; i < ntot; i++) {
          R->GenArray[egrp].rValue[elem][i] *= (1.0+P->GenArray[egrp].rValue[elem][0]);
        }
      }
    }
  
  }
  return xf_OK; 
}
