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
 FILE:  xf_Solver.c
 
 This file contains the high level solver driver functions.
 
 */

#include "xf_AllStruct.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_LinearSolverStruct.h"
#include "xf_SolverTools.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Residual.h"
#include "xf_ResidualStab.h"
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
#include "xf_Adapt.h"
#include "xf_AdaptStruct.h"
#include "xf_ErrEst.h"
#include "xf_All.h"
#include "xf_Penalty.h"
#include "xf_MeshMotion.h"
#include "xf_MeshMotionGCL.h"
#include "xf_Partition.h"




/* For communicating solver errors or flags within this file */
enum xfe_SolverFlag{
  xfe_SolverFlag_None,
  xfe_SolverFlag_Recoverable,
  xfe_SolverFlag_NotRecoverable,
  xfe_SolverFlag_Last,
};


/******************************************************************/
//   FUNCTION Definition: xf_CreateSolverData
int 
xf_CreateSolverData( xf_All *All, xf_SolverData **pSolverData)
{
  int ierr, i, nOutput;
  enum xfe_PreconditionerType Preconditioner;
  
  ierr = xf_Error(xf_Alloc((void **) pSolverData, 1, sizeof(xf_SolverData)));
  if (ierr != xf_OK) return ierr;
  
  (*pSolverData)->iIter          = 0;
  (*pSolverData)->CFL            = 0.0;
  (*pSolverData)->UpdateFrac     = 0.0;
  (*pSolverData)->ResNorm        = 0.0;
  (*pSolverData)->RnormPrev      = 0.0;
  (*pSolverData)->Rnorm          = 0.0;
  (*pSolverData)->ResPenaltyPrev = 1.0;
  (*pSolverData)->ResPenalty     = 1.0;
  (*pSolverData)->muPrev         = 1.0;
  (*pSolverData)->mu             = 1.0;
  (*pSolverData)->gfnormPrev     = 1.0;
  (*pSolverData)->gfnorm         = 1.0;
  (*pSolverData)->MaxCFLAchieved = 0.0;
  
  if (All->EqnSet->Outputs != NULL)
    nOutput = All->EqnSet->Outputs->nOutput;
  else
    nOutput = 0;
  
  ierr = xf_Error(xf_Alloc((void **) &(*pSolverData)->AdjResNorm, nOutput, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &(*pSolverData)->Output, nOutput, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<nOutput; i++){
    (*pSolverData)->AdjResNorm[i] = 0.0;
    (*pSolverData)->Output[i] = 0.0;
  }
  
  (*pSolverData)->normdU_prev = -1.0;
  (*pSolverData)->LinResTol = -1.0;
  (*pSolverData)->CRequired = xfe_False;
  (*pSolverData)->SortLines = xfe_False;
  (*pSolverData)->C = NULL; 
  
  // locate preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
                                     xfe_PreconditionerName, 
                                     (int ) xfe_PreconditionerLast, 
                                     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;
  
  // Check if penalization is on
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "PenalizeResidual", 
                            &(*pSolverData)->PenalizeResidual);
  if (ierr != xf_OK) return ierr;
  
  // Check if we are reusing dU
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "ReusedU", 
                            &(*pSolverData)->ReusedU);
  if (ierr != xf_OK) return ierr;
  
  // locate line connectivity matrix if preconditioner requires it
  xf_PreconditionerLineCheck(Preconditioner, &(*pSolverData)->CRequired, 
                             &(*pSolverData)->SortLines);
  if ((*pSolverData)->CRequired){
    ierr = xf_Error(xf_FindLineConnectivity(All, &(*pSolverData)->C));
    if (ierr != xf_OK) return ierr;
  }
  
  
  xf_InitStabData(&((*pSolverData)->StabData));
  (*pSolverData)->StabRequired = xfe_False;
  
  (*pSolverData)->c  = 0.; 
  (*pSolverData)->dt = NULL;
  (*pSolverData)->Pvec = NULL;
  
  (*pSolverData)->SkipParallelExchange = xfe_False;

  (*pSolverData)->ResidualOrderIncrement = 0;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroySolverData
int 
xf_DestroySolverData( xf_SolverData *SolverData)
{
  if (SolverData == NULL) return xf_OK;
  
  xf_Release( (void *) SolverData->AdjResNorm);
  xf_Release( (void *) SolverData->Output);
  xf_Release( (void *) SolverData);
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CheckSolverFlag
static int
xf_CheckSolverFlag(int inputerr, enum xfe_Bool Notify,
                   enum xfe_Bool ParallelFlag,
                   enum xfe_SolverFlag *ErrorFlagOut)
{
  int ierr;
  int myRank;
  enum xfe_SolverFlag ErrorFlag;
  /* Note, increasing value = increasing severity:
   0 : no error, return with xfe_False
   1 : possibly-recoverable error, return with xfe_False but (*Rewind) == True
   2 : unrecoverable error; return with xfe_True
   */
  
  ErrorFlag = (*ErrorFlagOut) = xfe_SolverFlag_None;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  switch(inputerr){
    case xf_OK:
      ErrorFlag = xfe_SolverFlag_None;
      break;
    case xf_NON_PHYSICAL:
      if (Notify) if (myRank == 0){ xf_printf("xf_NON_PHYSICAL occured.  Attempting to recover.\n"); }
      ErrorFlag = xfe_SolverFlag_Recoverable;
      break;
    case xf_SINGULAR:
      if (Notify) if (myRank == 0){ xf_printf("xf_SINGULAR occured.  Attempting to recover.\n"); }
      ErrorFlag = xfe_SolverFlag_Recoverable;
      break;
    case xf_NOT_CONVERGED:
      if (Notify) if (myRank == 0){ xf_printf("xf_NOT_CONVERGED occured.  Attempting to recover.\n"); }
      ErrorFlag = xfe_SolverFlag_Recoverable;
      break;
    case xf_NO_UPDATE:
      if (Notify) if (myRank == 0){ xf_printf("xf_NO_UPDATE occured.  Attempting to recover.\n"); }
      ErrorFlag = xfe_SolverFlag_Recoverable;
      break;
    default:
      if (Notify) if (myRank == 0){ xf_printf("Error = %d is not Recoverable.\n", inputerr); }
      ErrorFlag = xfe_SolverFlag_NotRecoverable;
      break;
  }
  
  // reduce-max ErrorFlag if in Parallel
  if (ParallelFlag){
    ierr = xf_Error(xf_MPI_Allreduce((int *) &ErrorFlag, 1, xfe_SizeInt, xfe_MPI_MAX));
    if (ierr != xf_OK) return ierr;
  }
  
  // store flag as output
  (*ErrorFlagOut) = ErrorFlag;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CheckSolverError
static enum xfe_Bool
xf_CheckSolverError(int inputerr, enum xfe_Bool Notify, xf_SolverData *SolverData,
                    xf_Vector *USafe, xf_Vector *U, enum xfe_Bool *Rewind)
{
  int ierr;
  enum xfe_SolverFlag ErrorFlag;
  
  // determine error flag
  ierr = xf_Error(xf_CheckSolverFlag(inputerr, Notify, 
                                     ((U != NULL) ? U->ParallelFlag : xfe_True), 
                                     &ErrorFlag));
  if (ierr != xf_OK) return ierr;
  
  if (ErrorFlag == xfe_SolverFlag_None){
    (*Rewind) = xfe_False; // all is ok
    return xfe_False;
  }
  else if (ErrorFlag == xfe_SolverFlag_NotRecoverable){
    (*Rewind) = xfe_True; // need to set rewind flag to catch solver error
    return xfe_True; // show's over
  }
  
  /* at this point we need to attempt to recover */
  (*Rewind) = xfe_True;
  
  // set U = USafe if data was provided
  if ((USafe != NULL) && (U != NULL)){
    ierr = xf_Error(xf_SetVector(USafe, xfe_Set, U));
    if (ierr != xf_OK) return ierr;
  }
  
  // Set CFL = CFL/CFLDecreaseFactor
  if (SolverData != NULL){
    SolverData->CFL = min(SolverData->CFLSafe, SolverData->CFL/SolverData->CFLDecreaseFactor);
    xf_printf("Decreasing CFL, now  = %.10E, min = %.10E\n", SolverData->CFL, SolverData->CFLMin);
    if (SolverData->CFL < SolverData->CFLMin){
      xf_printf("CFL is below CFLMin.  Exiting.\n");
      return xfe_True; // show's over
    }
  }
  
  return xfe_False;
}


/******************************************************************/
//   FUNCTION Definition: xf_CheckSolverErrorUnsteady
static enum xfe_Bool
xf_CheckSolverErrorUnsteady(int inputerr, enum xfe_Bool Notify, xf_SolverData *SolverData,
                            int nUSafe, xf_Vector **USafe, int nU, xf_Vector **U, 
                            real *pTimeStep, enum xfe_Bool *Rewind)
{
  int ierr, k;
  enum xfe_SolverFlag ErrorFlag;
  real TimeStepDecrease;
  
  // determine error flag
  ierr = xf_Error(xf_CheckSolverFlag(inputerr, Notify, 
                                     ((U != NULL) ? (*U)->ParallelFlag : xfe_True), 
                                     &ErrorFlag));
  if (ierr != xf_OK) return ierr;
  
  if (ErrorFlag == xfe_SolverFlag_None){
    (*Rewind) = xfe_False; // all is ok
    return xfe_False;
  }
  else if (ErrorFlag == xfe_SolverFlag_NotRecoverable){
    return xfe_True; // show's over
  }
  
  /* at this point we need to attempt to recover */
  (*Rewind) = xfe_True;
  
  // set U = USafe if data was provided
  if ((USafe != NULL) && (U != NULL)){
    if ((nUSafe != 1) && (nUSafe != nU)) return xf_Error(xf_INPUT_ERROR);
    for (k=0; k<nU; k++){
      ierr = xf_Error(xf_SetVector(USafe[k%nUSafe], xfe_Set, U[k]));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // Set TimeStep = TimeStep/TimeStepDecreaseFactor
  if (pTimeStep != NULL){
    TimeStepDecrease = ((SolverData != NULL) ? SolverData->TimeStepDecreaseFactor : 2.0);
    (*pTimeStep) = (*pTimeStep) / TimeStepDecrease;
    if (Notify)
      xf_printf("Decreasing TimeStep, now  = %.10E\n", (*pTimeStep));
  }

  
  // Set CFL = CFL/CFLDecreaseFactor
  if (SolverData != NULL){
    SolverData->CFL = min(SolverData->CFLSafe, SolverData->CFL/SolverData->CFLDecreaseFactor);
    if (Notify) xf_printf("Decreasing CFL, now  = %.10E, min = %.10E\n", 
			  SolverData->CFL, SolverData->CFLMin);
    if (SolverData->CFL < SolverData->CFLMin){
      xf_printf("CFL is below CFLMin.  Exiting.\n");
      return xfe_True; // show's over
    }
  }
  
  return xfe_False;
}




/******************************************************************/
//   FUNCTION Definition: xf_CalculateArtificialTimeStep
int
xf_CalculateArtificialTimeStep(xf_All *All, xf_Vector *U, real CFL, xf_Vector *dt)
{
  int ierr, dim, egrp, elem, Order, QuadOrder;
  int nq, pnq, nn, sr, *IParam, iAux, nAuxU;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged;
  real *xq, *EU, *u, *xglob, *RParam, **Auxu;
  real h, vmax, dtmax, dte;
  xf_BasisData *PhiData, *GeomPhiData, **AuxPhiData;
  xf_QuadData *QuadData;
  xf_Vector *EG, *V, **AuxU;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  sr   = All->EqnSet->StateRank;
  
  // obtain elem geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  if (ierr != xf_OK) return ierr;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, &nAuxU, &AuxU));
  if (ierr != xf_OK) return ierr;
  
  GeomPhiData = NULL;
  QuadData    = NULL;
  PhiData     = NULL;
  u           = NULL;
  xglob       = NULL;
  pnq         = -1;
  
  // Allocate vector of basis data for interpolating auxiliary vectors
  ierr = xf_Error(xf_Alloc( (void **) &AuxPhiData, nAuxU, sizeof(xf_BasisData *)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &Auxu, nAuxU, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  for (iAux=0; iAux<nAuxU; iAux++){
    AuxPhiData[iAux] = NULL;
    Auxu[iAux]       = NULL;
  }
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    Basis = U->Basis[egrp];
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // pull off order
      Order = xf_InterpOrder(U, egrp, elem);
      
      // determine eqn-set-based quadrature order
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, All->EqnSet, egrp, Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      /* Pull off quad points for the element; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **) &u, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
        if (ierr != xf_OK) return ierr;
      }
      
      // obtain global coords of quad points
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
                                      nq, xq, xglob));
      if (ierr != xf_OK) return ierr;
      
      nn = PhiData->nn;
      EU = U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
      
      // interpolate state
      xf_MxM_Set(PhiData->Phi, EU,  nq, nn, sr, u);
      
      // interpolate any auxiliary vectors, using AuxPhiData
      for (iAux=0; iAux<nAuxU; iAux++){
        V  = AuxU[iAux];
        if (nq > pnq){ // reallocate Auxu[iAux] if necessary
          ierr = xf_Error(xf_ReAlloc( (void **) Auxu+iAux, nq*V->StateRank, sizeof(real)));
          if (ierr != xf_OK) return ierr;
        }
        ierr = xf_Error(xf_EvalBasis(V->Basis[egrp], xf_InterpOrder(V,egrp,elem), QuadChanged, 
                                     nq, xq, xfb_Phi, AuxPhiData+iAux));
        if (ierr != xf_OK) return ierr;
        xf_MxM_Set(AuxPhiData[iAux]->Phi, V->GenArray[egrp].rValue[elem], nq, 
                   AuxPhiData[iAux]->nn, V->StateRank, Auxu[iAux]);
      }
      
      dtmax = -1.;
      
      ierr = xf_EqnSetMaxCharSpeed(All->EqnSet, nq, u, Auxu, xglob, IParam,
                                   RParam, &vmax, &dtmax, NULL, NULL);
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ElemSize(All, egrp, elem, EG, &h)); // A/perim
      if (ierr != xf_OK) return ierr;
      
      // set elemental dt
      dte = CFL*h/vmax;
      
      // if eqnset imposed a constraint on dt, use it
      if (dtmax > 0){
        dte = min(dte, CFL*dtmax);
      }
      
      dt->GenArray[egrp].rValue[elem][0] = dte;
      
      pnq = nq;
    } // elem
  } // egrp
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Auxiliary Vector Data */
  for (iAux=0; iAux<nAuxU; iAux++){
    ierr = xf_Error(xf_DestroyBasisData(AuxPhiData[iAux], xfe_True));
    if (ierr != xf_OK) return ierr;
    xf_Release( (void *) Auxu[iAux]);
  }
  xf_Release( (void *) AuxPhiData);
  xf_Release( (void *) Auxu);
  xf_Release( (void *) AuxU);
  
  // Release memory
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) u);
  xf_Release( (void *) xglob);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_LimitPoints
int 
xf_LimitPoints( xf_Mesh *Mesh, int egrp, int elem, xf_Vector *U,
               xf_EqnSet *EqnSet, int *pnq, real **pxq)
{
  int ierr, dim, k, nqtot;
  int face, egrpL, egrpR, faceorient;
  int Order, OrderL, OrderR, QuadOrder;
  enum xfe_Bool QuadChanged;
  xf_QuadData *QuadData;
  xf_IFace IFace;
  xf_BFace BFace;
  xf_Face Face;
  
  dim = Mesh->Dim;
  
  if (U->Order == NULL) return xf_Error(xf_INPUT_ERROR);
  if (U->Basis == NULL) return xf_Error(xf_INPUT_ERROR);
  
  (*pnq) = 0;
  
  // interpolation order
  Order = xf_InterpOrder(U, egrp, elem);
  
  // quadrature order on elem
  ierr = xf_Error(xf_GetQuadOrderElem(Mesh, EqnSet, egrp, Order, &QuadOrder));
  if (ierr != xf_OK) return ierr;
  
  // quad points on elem
  QuadData = NULL;
  ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
  if (ierr != xf_OK) return ierr;
  
  // reallocate memory
  nqtot = (*pnq) + QuadData->nquad;
  ierr = xf_Error(xf_ReAlloc((void **) pxq, nqtot*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  for (k=0; k<dim*QuadData->nquad; k++) (*pxq)[(*pnq)+k] = QuadData->xquad[k];
  
  (*pnq) = nqtot;
  
  
  // loop over faces
  for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
    Face = Mesh->ElemGroup[egrp].Face[elem][face];
    if (Face.Group == xf_INTERIORFACE){ // interior face
      IFace = Mesh->IFace[Face.Number];
      egrpL = IFace.ElemGroupL;
      egrpR = IFace.ElemGroupR;
      
      OrderL = xf_InterpOrder(U, egrpL, IFace.ElemL);
      OrderR = xf_InterpOrder(U, egrpR, IFace.ElemR);
      
      ierr = xf_Error(xf_GetQuadOrderIFace(Mesh, EqnSet, IFace, max(OrderL, OrderR),
                                           &QuadOrder));
      if (ierr != xf_OK) return ierr;
      if (egrpL == egrp)
        faceorient = IFace.OrientL;
      else if (egrpR == egrp)
        faceorient = IFace.OrientR;
      else
        return xf_Error(xf_CODE_LOGIC_ERROR);
    }
    else if (Face.Group >= 0){ // boundary face
      BFace = Mesh->BFaceGroup[Face.Group].BFace[Face.Number];
      ierr = xf_Error(xf_GetQuadOrderBFace(Mesh, EqnSet, BFace, Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      faceorient = BFace.Orient;
    }
    else continue;
    
    // quad points on face 
    ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face, 
                                QuadOrder, &QuadData, &QuadChanged));
    if (ierr != xf_OK) return ierr;
    
    // reallocate memory
    nqtot = (*pnq) + QuadData->nquad;
    ierr = xf_Error(xf_ReAlloc((void **) pxq, nqtot*dim, sizeof(real)));
    
    // convert xq (face coords) to elem ref coords
    ierr = xf_Error(xf_RefFace2Interpol(Mesh, egrp, elem, face, faceorient,
                                        QuadData->nquad, QuadData->xquad, 
                                        (*pxq) + dim*(*pnq)));
    if (ierr != xf_OK) return ierr;
    
    (*pnq) = nqtot;
    
  } // face
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UpdateStateMaxPrimChange
static int
xf_UpdateStateMaxPrimChange(xf_All *All, xf_Vector *U, xf_Vector *dU,
                            enum xfe_Bool *LimitFlag, 
                            enum xfe_Bool *UpdateFlag,
                            xf_SolverData *SolverData)
{
  /* General idea:
   set update frac = 1
   loop over elements; set frac = min(frac, current)
   if below minfrac, add to list, do not decrease frac
   if list not empty and not allowing local updates, return no update
   if list empty, apply frac to all elements
   if list not empty, apply frac to all elements except those in list
   */  
  int ierr, k, egrp, elem, Order;
  int nq, pnq, nn, sr, r, ib, nb, nbmax, nb0;
  int *IParam = NULL;
  enum xfe_BasisType Basis;
  enum xfe_Bool LocUpdate;
  real MinUpdateFrac, UpdateFrac, frac;
  real *xq, *EU, *EdU, *u, *du;
  real *RParam = NULL;
  xf_ElemReal *blist;
  xf_BasisData *PhiData;
  xf_Mesh *Mesh;
  
  
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;
  
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "MinUpdateFraction", 
                            &MinUpdateFrac);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "LocalStateUpdate", 
                            &LocUpdate);
  if (ierr != xf_OK) return ierr;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, 
                                       NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  if (SolverData != NULL) SolverData->UpdateFrac = 0.0;
  UpdateFrac = 1.0;
  nb = 0; nbmax = 0; nb0 = 10;
  blist  = NULL;
  
  PhiData  = NULL;
  u = du = NULL;
  xq = NULL;
  pnq = -1;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    Basis = U->Basis[egrp];
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // pull off order
      Order = xf_InterpOrder(U, egrp, elem);
      
      /* Pull off limit points -- effectively the quad points,
       including those on faces */
      ierr = xf_Error(xf_LimitPoints(Mesh, egrp, elem, U, All->EqnSet, &nq, &xq));
      if (ierr != xf_OK) return ierr;
      
      // compute basis functions 
      ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **)  &u,  nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **)  &du, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
      }
      
      nn = PhiData->nn;
      EU  =  U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
      EdU = dU->GenArray[egrp].rValue[elem]; // dU on elem [nn*sr]
      
      xf_MxM_Set(PhiData->Phi, EU,  nq, nn, sr, u);
      xf_MxM_Set(PhiData->Phi, EdU, nq, nn, sr, du);
      
      ierr = xf_Error(xf_EqnSetUpdateFraction(All->EqnSet, IParam, RParam, nq, u, du, &frac));
      if (ierr != xf_OK) return ierr;
      
      if (frac < MinUpdateFrac){
        if (!LocUpdate){ 
          UpdateFrac = -1.0;  // means MinUpdateFrac was violated
        }
        else{ // allow local updates; save update info
          nb++;
          if (nb > nbmax){
            nbmax += nb0;
            nb0 *= 2;
            ierr = xf_Error(xf_ReAlloc((void **) blist, nbmax, sizeof(xf_ElemReal)));
            if (ierr != xf_OK) return ierr;
          }
          blist[nb-1].egrp  = egrp;
          blist[nb-1].elem  = elem;
          blist[nb-1].value = frac;
        }
      }
      else{
        UpdateFrac = min(UpdateFrac, frac);
      }
    } // elem
  } // egrp
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  // release memory
  xf_Release( (void *) u);
  xf_Release( (void *) du);
  xf_Release( (void *) xq);
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  
  
  // reduce-min UpdateFrac
  ierr = xf_Error(xf_MPI_Allreduce(&UpdateFrac, 1, xfe_SizeReal, xfe_MPI_MIN));
  if (ierr != xf_OK) return ierr;
  
  // negative UpdateFrac means cannot go on (MinUpdateFrac violated)
  if (UpdateFrac < 0){
    if (SolverData != NULL) SolverData->UpdateFrac = 0.0;
    (*UpdateFlag) = xfe_False;
    (*LimitFlag)  = xfe_False;
    return xf_OK; // exit gracefully with UpdateFlag = False;
  }    
  
  if (SolverData != NULL) SolverData->UpdateFrac = UpdateFrac;
  
  if (UpdateFrac < 1.0) 
    (*LimitFlag) = xfe_True;
  else
    (*LimitFlag) = xfe_False;
  (*UpdateFlag) = xfe_True;
  
  // apply update to U
  ib = 0;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      EU  =  U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
      EdU = dU->GenArray[egrp].rValue[elem]; // dU on elem [nn*sr]
      r = ((U->GenArray[egrp].vr==NULL) ? U->GenArray[egrp].r : U->GenArray[egrp].vr[elem]);
      frac = UpdateFrac;
      if ((nb > ib) && (blist[ib].egrp == egrp) && (blist[ib].elem == elem)){
        frac = blist[ib].value;
        ib++;
      }
      for (k=0; k<r; k++) EU[k] += EdU[k]*frac;
    } // elem
  } // egrp
  
  // Release blist
  xf_Release( (void *) blist);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UpdateStateLineSearch
static int
xf_UpdateStateLineSearch(xf_All *All, xf_Vector *U, xf_Vector *dU,
               enum xfe_Bool *LimitFlag, enum xfe_Bool *UpdateFlag,
               xf_SolverData *SolverData)
{
  
  int ierr, k, egrp, elem, Order, ntot;
  int nq, pnq, nn, sr, r, terr;
  int *IParam = NULL;
  enum xfe_BasisType Basis;
  enum xfe_Bool Rewind, Greedy;
  real MinUpdateFrac, frac, pFrac, proj, ampfac;
  real *xq, *EU, *EdU, *u, *du, ReductionFactor;
  real *RParam = NULL, *RStateNorms, *RtStateNorms;
  xf_BasisData *PhiData;
  xf_Mesh *Mesh;
  xf_Vector *Utrial, *Rt, *R, *dt, *Rs, *Dt, *Gp;
  xf_Data *D;
  
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;

  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "MinUpdateFraction",
                            &MinUpdateFrac);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "UpdateFracReduction",
                            &ReductionFactor);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "LineSearchGreedy", 
                            &Greedy);
  if (ierr != xf_OK) return ierr;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, 
                                       NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  //SolverData stores essential information
  if (SolverData == NULL) return xf_Error(xf_INPUT_ERROR);
  
  SolverData->UpdateFrac = 0.0;
  
  PhiData  = NULL;
  u = du = NULL;
  xq = NULL;
  pnq = -1;
  
  pFrac = 1.0;
  
  ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "TimeStep", xfe_Vector, &D));
  if (ierr != xf_OK) return ierr;
  dt = (xf_Vector *)D->Data;
  
  /* create a timestep vector compatible with U */
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Dt", xfe_False, 
                                       xfe_False, NULL, &Dt, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "Residual", xfe_Vector, &D));
  if (ierr != xf_OK) return ierr;
  R = (xf_Vector *)D->Data;
  
  ierr = xf_Error(xf_FindSimilarVector(All, R, "Rs", xfe_False, 
                                       xfe_False, NULL, &Rs, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **)&RStateNorms, sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  /* if penalization is off, residual only gets updated */
  ierr = xf_Error(xf_PenalizeResidual(All, U, Rs, SolverData));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VectorStateNorms(All, Rs, RStateNorms, 2));
  if (ierr != xf_OK) return ierr;
  
  //find the penalty gradient vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "GradPenalty", xfe_True, 
                                       xfe_True, NULL, &Gp, NULL));
  if (ierr != xf_OK) return ierr;
  
  //calculate the penalty function gradient
  ierr = xf_Error(xf_CalculatePenaltyGradient(All, U, Gp));
  if (ierr != xf_OK) return ierr;
  
  //get max update that keeps the state physical
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    Basis = U->Basis[egrp];
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      if (U->GenArray[egrp].vr == NULL)
        ntot = U->GenArray[egrp].r;
      else 
        ntot = U->GenArray[egrp].vr[elem];
      
      //verify if dU is in the safe direction
      proj = 0.0;
      for (k = 0; k < ntot; k++){
        Dt->GenArray[egrp].rValue[elem][k] = dt->GenArray[egrp].rValue[elem][0];
        proj = max(dU->GenArray[egrp].rValue[elem][k]*
                           Gp->GenArray[egrp].rValue[elem][k],0.0);
      }
      if (proj > 0.0){//not in the safe direction
        // pull off order
        Order = xf_InterpOrder(U, egrp, elem);
        
        /* Pull off limit points -- effectively the quad points,
         including those on faces */
        ierr = xf_Error(xf_LimitPoints(Mesh, egrp, elem, U, All->EqnSet, &nq, &xq));
        if (ierr != xf_OK) return ierr;
        
        // compute basis functions 
        ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xq, xfb_Phi, &PhiData));
        if (ierr != xf_OK) return ierr;
        
        // re-allocate data if quad points increased
        if (nq > pnq){
          ierr = xf_Error(xf_ReAlloc( (void **)  &u,  nq*sr, sizeof(real)));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_ReAlloc( (void **)  &du, nq*sr, sizeof(real)));
          if (ierr != xf_OK) return ierr;
        }
        
        nn = PhiData->nn;
        EU  =  U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
        EdU = dU->GenArray[egrp].rValue[elem]; // dU on elem [nn*sr]
        
        xf_MxM_Set(PhiData->Phi, EU,  nq, nn, sr, u);
        xf_MxM_Set(PhiData->Phi, EdU, nq, nn, sr, du);
        
        ierr = xf_Error(xf_EqnSetUpdateFraction(All->EqnSet, IParam, RParam, 
                                                nq, u, du, &frac));
        if (ierr != xf_OK) return ierr;
        pFrac = min(pFrac,frac);
      }
    }
  }
  
  ierr = xf_Error(xf_MPI_Allreduce(&pFrac, 1, xfe_SizeReal, xfe_MPI_MIN));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  // release memory
  xf_Release( (void *) u);
  xf_Release( (void *) du);
  xf_Release( (void *) xq);
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  
  ierr = xf_Error(xf_FindSimilarVector(All, R, "Rt", xfe_False, 
                                       xfe_False, NULL, &Rt, NULL));
  if (ierr != xf_OK) return ierr;
  /* create trial state vector */
  ierr = xf_Error(xf_CreateVector(&Utrial));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_CopyVector(Mesh, U, Utrial));
  if (ierr != xf_OK) return ierr;
  
  frac = 1.0;
  (*UpdateFlag) = xfe_False;
  (*LimitFlag) = xfe_True;
  
  ierr = xf_Error(xf_Alloc((void **)&RtStateNorms, sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  
  while (!(*UpdateFlag) && frac*pFrac >= MinUpdateFrac) {
  /* set Utrial = U + frac*dU */
    ierr = xf_Error(xf_SetVector(U, xfe_Set, Utrial));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_VectorMultSet(dU, frac*pFrac, xfe_Add, Utrial));
    if (ierr != xf_OK) return ierr;
    
    /* calculate "unsteady" term  */
    ierr = xf_Error(xf_SetVector(Dt, xfe_Set, Rt));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_VectorInv(Rt));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_VectorVectorMult(dU, Rt));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_MultMassMatrix(All, frac*pFrac, Rt));
    if (ierr != xf_OK) return ierr;
        
    /* add steady part */
    ierr = xf_Error(xf_PenalizeResidual(All, Utrial, Rs, SolverData));
    terr = xf_Error(xf_MPI_Allreduce(&ierr, 1, xfe_SizeInt, xfe_MPI_MAX));
    if (terr != xf_OK) return terr;

    if (ierr == xf_NON_PHYSICAL){
      frac *= ReductionFactor;
      continue;
    }
    else if (ierr != xf_OK) return ierr;
        
    ierr = xf_Error(xf_SetVector(Rs, xfe_Add, Rt));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_VectorStateNorms(All, Rt, RtStateNorms, 2));
    if (ierr != xf_OK) return ierr;
    
    Rewind = xfe_False;
    for (r = 0; r < sr; r++){ 
      if ((RtStateNorms[r] > RStateNorms[r])) {
        Rewind = xfe_True;
        break;
      }
    }
    
    if (Rewind){
      frac *= ReductionFactor;
      continue;
    }
    else{
      (*UpdateFlag) = xfe_True;
      break;
    }
    SolverData->UpdateFrac = 0.0;
  }
  
  ampfac = 1.1;
  //Can we get greedy?
  ierr=xf_Error(xf_VectorDot(dU, Gp, &proj));
  if (ierr != xf_OK) return ierr;
  
  if (fabs(frac-1.)<MEPS && Greedy && proj < 0.0)//no back-tracking occurred
    //xf_printf("\n\n\n\nGreedy\n\n\n");
    while ((frac*pFrac < 1.) && (pFrac < 1.)) {
      frac *= ampfac;//amplify frac
      /* set Utrial = U + frac*dU */
      ierr = xf_Error(xf_SetVector(U, xfe_Set, Utrial));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_VectorMultSet(dU, frac*pFrac, xfe_Add, Utrial));
      if (ierr != xf_OK) return ierr;
      
      /* calculate "unsteady" term  */
      ierr = xf_Error(xf_SetVector(Dt, xfe_Set, Rt));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_VectorInv(Rt));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_VectorVectorMult(dU, Rt));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_MultMassMatrix(All, frac*pFrac, Rt));
      if (ierr != xf_OK) return ierr;
      
      /* add steady part */
      ierr = xf_Error(xf_PenalizeResidual(All, Utrial, Rs, SolverData));
      terr = xf_Error(xf_MPI_Allreduce(&ierr, 1, xfe_SizeInt, xfe_MPI_MAX));
      if (terr != xf_OK) return terr;
      
      if (ierr == xf_NON_PHYSICAL){
        frac *= ReductionFactor;//back-track and exit
        break;
      }
      else if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_SetVector(Rs, xfe_Add, Rt));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_VectorStateNorms(All, Rt, RtStateNorms, 2));
      if (ierr != xf_OK) return ierr;
      
      Rewind = xfe_False;
      for (r = 0; r < sr; r++){ 
        if ((RtStateNorms[r] > RStateNorms[r])) {
          Rewind = xfe_True;
          break;
        }
      }
      
      if (Rewind){
        frac /= ampfac;//back-track and exit
        break;
      }
      
    }
  
  if (frac*pFrac < MinUpdateFrac){
    SolverData->UpdateFrac = 0.0;
    (*UpdateFlag) = xfe_False;
    (*LimitFlag) = xfe_True;
  }
  else{ 
    (*UpdateFlag) = xfe_True;
    if(frac*pFrac > 1.)
      SolverData->UpdateFrac = 1.;
    else
      SolverData->UpdateFrac = frac*pFrac;
  }
  
  if ((*UpdateFlag)) {
    if (frac*pFrac < 1.0) (*LimitFlag) = xfe_True;
    else (*LimitFlag) = xfe_False;
    ierr = xf_Error(xf_VectorMultSet(dU, SolverData->UpdateFrac,
                                     xfe_Add, U));
    if (ierr != xf_OK) return ierr;
  }
  
  ierr = xf_Error(xf_DestroyVector(Rt, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(Rs, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(Utrial, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(Dt, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  xf_Release((void*)RStateNorms);
  xf_Release((void*)RtStateNorms);

  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UpdateState
static int
xf_UpdateState(xf_All *All, xf_Vector *U, xf_Vector *dU,
               enum xfe_Bool *LimitFlag, enum xfe_Bool *UpdateFlag,
               xf_SolverData *SolverData)
{
  int ierr;
  enum xfe_StateUpdateMethodType UpdateMethod;
  
  ierr = xf_GetKeyValueEnum(All->Param->KeyValue, "StateUpdateMethod", 
                            xfe_StateUpdateMethodName, 
                            (int)xfe_StateUpdateMethodLast, 
                            (int*)&UpdateMethod);
  if (ierr != xf_OK) return ierr;
  
  switch (UpdateMethod) {
    case xfe_StateUpdateMaxPrimChange:
      /* limited primitive variables change */
      ierr = xf_Error(xf_UpdateStateMaxPrimChange(All, U, dU, LimitFlag, 
                                                  UpdateFlag, SolverData));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_StateUpdateLineSearch:
      /* linesearch */
      ierr = xf_Error(xf_UpdateStateLineSearch(All, U, dU, LimitFlag, 
                                               UpdateFlag, SolverData));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FindCFLData
static int
xf_FindCFLData(xf_KeyValue KeyValue, xf_SolverData *SolverData)
{
  int ierr;
  
  ierr = xf_GetKeyValueReal(KeyValue, "CFL", &SolverData->CFL);
  if (ierr != xf_OK) return ierr;
  SolverData->CFLSafe = SolverData->CFL;
  
  ierr = xf_GetKeyValueReal(KeyValue, "CFLIncreaseFactor", 
                            &SolverData->CFLIncreaseFactor);
  if (ierr != xf_OK) return ierr;
  ierr = xf_GetKeyValueReal(KeyValue, "CFLDecreaseFactor", 
                            &SolverData->CFLDecreaseFactor);
  if (ierr != xf_OK) return ierr;
  ierr = xf_GetKeyValueReal(KeyValue, "CFLMax", &SolverData->CFLMax);
  if (ierr != xf_OK) return ierr;
  ierr = xf_GetKeyValueReal(KeyValue, "CFLMin", &SolverData->CFLMin);
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValueEnum(KeyValue, "CFLEvolution", 
                                     xfe_CFLEvolutionName, (int )xfe_CFL_Last, 
                                     (int *)&SolverData->CFLEvolution));
  if (ierr != xf_OK) return ierr;
  ierr = xf_GetKeyValueReal(KeyValue, "TimeStepDecreaseFactor", 
                            &SolverData->TimeStepDecreaseFactor);
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_EvolveCFL
static int
xf_EvolveCFL(xf_SolverData *SolverData, enum xfe_Bool LimitFlag)
{
  real ratio1, CFLSERA, CFLSERB;
  
  SolverData->CFLPrev = SolverData->CFL;
  
  switch (SolverData->CFLEvolution) {
    case xfe_CFL_Exp:
      if (!LimitFlag){
        SolverData->CFL *= SolverData->CFLIncreaseFactor;
        SolverData->CFL = min(SolverData->CFL, SolverData->CFLMax);
      }
      break;
    case xfe_CFL_SERA:
      ratio1 = SolverData->RnormPrev/SolverData->Rnorm;
      SolverData->CFL = min(SolverData->CFLPrev*ratio1,
                            SolverData->CFLMax);
      break;
    case xfe_CFL_SERB:
      ratio1 = SolverData->CFL/SolverData->normdU_prev;
      SolverData->CFL = min(ratio1,SolverData->CFLMax);
      break;
    case xfe_CFL_SERAB:
      //SER-A
      ratio1 = SolverData->RnormPrev/SolverData->Rnorm;
      CFLSERA = min(SolverData->CFLPrev*ratio1,
                    SolverData->CFLMax);
      
      //SER-B
      ratio1 = SolverData->CFL/SolverData->normdU_prev;
      CFLSERB = min(ratio1,SolverData->CFLMax);
      SolverData->CFL = 0.5*(CFLSERA+CFLSERB);
      break;
    case xfe_CFL_RDM:
        ratio1 = (SolverData->RnormPrev-SolverData->Rnorm)/SolverData->RnormPrev;
        SolverData->CFL *= min(pow(SolverData->CFLIncreaseFactor,ratio1),
                               SolverData->CFLMax);
      break;
    case xfe_CFL_mRDM:
      ratio1 = (SolverData->RnormPrev-SolverData->Rnorm)/SolverData->RnormPrev;
      SolverData->CFL *= min(pow(SolverData->CFLIncreaseFactor,max(ratio1,0.0)),
                             SolverData->CFLMax);
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_PreconditionerLeanCheck
int
xf_PreconditionerLeanCheck(enum xfe_PreconditionerType Preconditioner,
                           enum xfe_Bool *pLeanFlag)
{
  /* Checks if the preconditioner is memory-lean */
  
  (*pLeanFlag) = (Preconditioner == xfe_PreconditionerBlockJacobiLean);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_InitPenFactor
static int
xf_InitPenFactor(xf_All *All, xf_Vector *U, xf_SolverData *SolverData)
{
  int ierr, egrp, elem, sr, nn, n, r, i, j, is, js, nelem_glob;
  real P, MdtFNorm, PFNorm, temp, *Mdata, ratio;
  
  
  //set penalty factor to 1 and calculate penalty
  ierr = xf_Error(xf_SetKeyValueReal(All->EqnSet->KeyValue, 
                                     "PenaltyFcnFactor", 
                                     1.0));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_CalculatePenalty(All, U, SolverData->Pvec, &P));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetnElem(All->Mesh, NULL, &nelem_glob));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Allreduce(&nelem_glob, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  //P = P/nelem_glob;
  
  SolverData->mu = (pow(10.0,0.25)-1.0)/P;//1.0/ratio;;
    
  ierr = xf_Error(xf_SetKeyValueReal(All->EqnSet->KeyValue, 
                                     "PenaltyFcnFactor", 
                                     SolverData->mu));
  if (ierr != xf_OK) return ierr;
  
  SolverData->ResPenalty = 1.0+P;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AddPenaltyMatrix
int
xf_AddPenaltyMatrix(xf_All *All, xf_Vector *U, xf_Vector *R, 
                    xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
  int ierr, elem, egrp, i, j, is, js, nn, sr, sr2, nelem_glob;
  real P, aux1, *A, fac;
  xf_Vector *Gp, *Pvec;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "GradPenalty", xfe_False, 
                                       xfe_True, NULL, &Gp, NULL));
  if (ierr != xf_OK) return ierr;
  
  Pvec = SolverData->Pvec;
  
  //compute penalty quantities
  ierr = xf_Error(xf_CalculatePenalty(All, U, SolverData->Pvec, &P));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetnElem(All->Mesh, NULL, &nelem_glob));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Allreduce(&nelem_glob, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  //P = P/nelem_glob;
  //store scalars for logging and evolving mu
  SolverData->ResPenalty = 1.0+P;
  
  ierr = xf_Error(xf_CalculatePenaltyGradient(All, U, Gp));
  if (ierr != xf_OK) return ierr;
  
  //compute augumented Jacobian
  sr = U->StateRank;
  sr2 = sr*sr;
  for (egrp=0; egrp<All->Mesh->nElemGroup; egrp++){
    for (elem=0; elem<All->Mesh->ElemGroup[egrp].nElem; elem++){
      if (U->GenArray[egrp].vr == NULL)
        nn = U->GenArray[egrp].r/sr;
      else 
        nn = U->GenArray[egrp].vr[elem]/sr;
      fac = 1.0/(1.0+SolverData->Pvec->GenArray[egrp].rValue[elem][0]);
      //point to self block of the Jacobian
      A = R_U->Value[egrp][elem][0];
      for (j = 0; j < nn; j++) {
        for (i = 0; i < nn; i++) {
          for (js = 0; js < sr; js++) {
            for (is = 0; is < sr; is++) {
              aux1 = Gp->GenArray[egrp].rValue[elem][i*sr+is]*R->GenArray[egrp].rValue[elem][j*sr+js];
              A[j*nn*sr2 + i*sr2 + js*sr + is] += fac*aux1;
            }//is
          }//js
        }//i
      }//j
    }//elem
  }//egrp
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SolveNonlinearSystem_Newton
static int
xf_SolveNonlinearSystem_Newton(xf_All *All, real c, enum xfe_Bool LinearFlag,
                               enum xfe_Bool ErrorSolve, xf_Vector *S, xf_Vector *U,
                               xf_SolverData **pSolverData)
{
  /*
   PURPOSE:
   
   Solves:
   
   c*M*U + R(U) + S = 0
   
   using a Newton method.  For standard steady-state runs, pass c=0.
   For unsteady runs, c should be the appropriate coefficient from the
   time discretization.
   
   For robustness, artificial time-stepping is used in the course of
   the nonlinear solve, so that at each nonlinear iteration, the
   system looks like
   
   (c + 1/dta)*M*U + R(U) + S = 0,
   
   where dta is the (possibly element-local) artificial time step.
   During the course of the nonlinear solve, 1/dta -> 0.
   
   if ErrorSolve == True, and S != NULL, S is modified on the first
   iteration according to:
   
   S = S - (c*M*U0 + R(U0))
   
   This behavior is used in coarse level solves in p-Multigrid, with FAS.
   
   INPUTS:
   
   All: All structure
   c : constant in front of M*U product
   LinearFlag : If True, the Jacobian matrix will not be recalculated if
   it already exists.
   ErrorSolve : If True and S != NULL, source is modified on first iteration
   S : source vector (or NULL)
   U : state vector
   
   OUTPUTS:
   
   U : modified state vector
   
   RETURN:
   
   Error code
   
   */
  
  int ierr, iIter0, nIter, nRewindMax, nRewind;
  int WriteInterval;
  enum xfe_Bool UpdateFlag, LimitFlag, Halted, Rewind, Converged, RestartFlag;
  enum xfe_Bool WriteResidual, WriteUpdate, Found, LeanFlag, AdaptRobust;
  enum xfe_Bool SolverError = xfe_False;
  enum xfe_PreconditionerType Preconditioner;
  enum xfe_Verbosity Verbosity;
  char OutputFile[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN];
  real ResidualTolerance, ratio1, ratio2, ratio3, ratioMax;
  xf_SolverData *SolverData;
  xf_Vector *R, *dU, *USafe, *dt, *Pvec, *UInitial, *gf;
  xf_JacobianMatrix *R_U, *R_U0;
  xf_Data *D;
  
  // locate preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
                                     xfe_PreconditionerName, (int ) xfe_PreconditionerLast, 
                                     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;
  
  // Is the preconditioner memory-lean?
  ierr = xf_Error(xf_PreconditionerLeanCheck(Preconditioner, &LeanFlag));
  if (ierr != xf_OK) return ierr;
  
  
  // should we make the residual vector writable?
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "WriteResidual", &WriteResidual);
  if (ierr != xf_OK) return ierr;
  
  // should we make the update vector writable?
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "WriteUpdate", &WriteUpdate);
  if (ierr != xf_OK) return ierr;
  
  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  // Steady write interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "SteadyWriteInterval", 
                                    &WriteInterval));
  if (ierr != xf_OK) return ierr;
  
  // SavePrefix
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;
  
  
  // number of nonlinear iterations
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", &iIter0);
  if (ierr != xf_OK) return ierr;
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", &nIter);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "Restart", &RestartFlag));
  if (ierr != xf_OK) return ierr;
  
  // create/allocate SolverData
  if (pSolverData == NULL){
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    if (ierr != xf_OK) return ierr;
  }
  else {
    SolverData = (*pSolverData);
  }
  
  // Set time-stepping constant
  SolverData->c = c;
  
  // locate Jacobian vector (do not allocate values if lean)
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        !LeanFlag, NULL, &R_U0, &Found));
  if (ierr != xf_OK) return ierr;
  
  // If Jacobian exists and system is linear, do not recalculate the Jacobian
  R_U = ((LinearFlag && Found) ? NULL : R_U0);
  
  // locate Residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, xfe_True, &D, &R, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = WriteResidual;
  if (WriteResidual && (nIter==0)){ // avoid writing uninitialized data
    ierr = xf_Error(xf_SetZeroVector(R));
    if (ierr != xf_OK) return ierr;
  }
  
  // locate update vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_True, xfe_True, &D, &dU, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = WriteUpdate;
  if (WriteUpdate){ // avoid writing uninitialized data in case never compute dU
    ierr = xf_Error(xf_SetZeroVector(dU));
    if (ierr != xf_OK) return ierr;
  }
  
  // locate "safe" U vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "USafe", xfe_False, xfe_True, NULL, &USafe, NULL));
  if (ierr != xf_OK) return ierr;
  
  // locate time step vector for artificial time stepping 
  if (LinearFlag)
    dt = NULL; // (do not need if linear)
  else{
    ierr = xf_Error(xf_FindVector(All, "TimeStep", xfe_LinkageGlobElem, 1, NULL, 0, 0, NULL, 
                                  NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  xfe_True, NULL, 
                                  &dt, NULL));
    if (ierr != xf_OK) return ierr;
  }
  if (SolverData->PenalizeResidual){
    ierr = xf_Error(xf_FindVector(All, "PenaltyIntegral", xfe_LinkageGlobElem, 1, NULL, 0, 0, NULL, 
                                  NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  xfe_True, NULL, 
                                  &Pvec, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  // pull off CFL, CFLIncreaseFactor, CFLDecreaseFactor, etc.
  ierr = xf_Error(xf_FindCFLData(All->Param->KeyValue, SolverData));
  if (ierr != xf_OK) return ierr;
  
  
  // pull off residual tolerance
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "ResidualTolerance", &ResidualTolerance);
  if (ierr != xf_OK) return ierr;
  
  // set Safe quantities
  SolverData->CFLSafe = SolverData->CFL;
  ierr = xf_Error(xf_SetVector(U, xfe_Set, USafe));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "AdaptRobust", &AdaptRobust));
  if (ierr != xf_OK) return ierr;
  
  if (AdaptRobust){
    // locate UInitial vector (for robustness adaptation)
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "UInitial", xfe_Vector, &D));
    if (ierr != xf_OK) return ierr;
    UInitial = (xf_Vector *) D->Data;
  }
  
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nCFLReducedMax", &nRewindMax));
  if (ierr != xf_OK) return ierr;
  nRewind = 0;
  
  // start nonlinear outer steps
  Converged = xfe_False;
  for (SolverData->iIter=iIter0; SolverData->iIter < (iIter0 + nIter); SolverData->iIter++){
    
    Halted     = xfe_False;
    UpdateFlag = xfe_False;
    Rewind     = xfe_False;
    while (UpdateFlag == xfe_False){
      if (Halted = xf_CheckUserHalt(NULL)) break;
      if (Rewind) nRewind++;
      if (AdaptRobust && ((nRewind >= nRewindMax) || (SolverError))){
        ierr = xf_Error(xf_SetVector(UInitial, xfe_Set, U));
        if (ierr != xf_OK) return ierr;
        
        return xf_NEED_ADAPT;
      }
      if (SolverError) return xf_Error(xf_SOLVER_ERROR);
      
      // convert CFL to (local) artificial time step if need one
      if (dt != NULL){
        ierr = xf_CalculateArtificialTimeStep(All, U, SolverData->CFL, dt);
        SolverError = xf_CheckSolverError(ierr, xfe_True, SolverData, USafe, U, &Rewind);
        if (Rewind){
          continue;
        }
        if ((SolverData->CFLSafe > SolverData->MaxCFLAchieved)&&
            (SolverData->iIter > iIter0)){
          SolverData->MaxCFLAchieved = SolverData->CFLSafe;
          if (AdaptRobust){
            ierr = xf_Error(xf_SetVector(U, xfe_Set, UInitial));
            if (ierr != xf_OK) return ierr;
            
          }
        }
      }
      
      SolverData->dt = dt; // store pointer to dt in SolverData
      
      // calculate/add spatial residual + Jacobian
      ierr = xf_CalculateResidual(All, U, R, R_U, SolverData);
      SolverError = xf_CheckSolverError(ierr, xfe_True, SolverData, USafe, U, &Rewind);
      if (Rewind) continue;
      
      if (SolverData->PenalizeResidual && SolverData->iIter == iIter0) {
        //initialize pvec to 1
        ierr = xf_Error(xf_SetConstVector(Pvec, 0,0.0));
        if (ierr != xf_OK) return ierr;
        
        SolverData->Pvec = Pvec;//store pointer to penalty vector
      }
      
      
      // add c*M*U to residual and (c+1/dt)*M to Jacobian.  Note, dt = artificial time step
      ierr = xf_Error(xf_AddMassMatrix(All, c, dt, U, R, R_U, SolverData));
      if (ierr != xf_OK) return ierr;
      
      SolverData->muPrev = SolverData->mu;
      if (SolverData->PenalizeResidual){
        if (SolverData->iIter == iIter0) {
          if (!RestartFlag){
            //calculate initial mu
            ierr = xf_Error(xf_InitPenFactor(All, U, SolverData));
            if (ierr != xf_OK) return ierr;
            Pvec = SolverData->Pvec;
          }
          else{
            ierr = xf_Error(xf_GetKeyValueReal(All->EqnSet->KeyValue, 
                                               "PenaltyFcnFactor", 
                                               &SolverData->mu));
            if (ierr != xf_OK) return ierr;
          }
          SolverData->muPrev = SolverData->mu;
        }
        else {
          //SER for mu
          ratio1 = SolverData->ResPenalty/SolverData->ResPenaltyPrev;
          //maximum mu hard coded as 1e30
          SolverData->mu = min(SolverData->muPrev*ratio1,1e30);
          ierr = xf_Error(xf_SetKeyValueReal(All->EqnSet->KeyValue, 
                                             "PenaltyFcnFactor", 
                                             SolverData->mu));
          if (ierr != xf_OK) return ierr;
        }
        SolverData->ResPenaltyPrev = SolverData->ResPenalty;
        //add penalization to Jacobian
        ierr = xf_Error(xf_AddPenaltyMatrix(All, U, R, R_U, SolverData));
        if (ierr != xf_OK) return ierr;
        
        //store initial penalization 
        if (SolverData->iIter == iIter0)
          SolverData->ResPenaltyPrev = SolverData->ResPenalty;
      }
      
      // modify source for error solve if requested
      if ((SolverData->iIter==iIter0) && (S!=NULL) && (ErrorSolve==xfe_True)){
        ierr = xf_Error(xf_SetVector(R, xfe_Sub, S));
        if (ierr != xf_OK) return ierr;
      }
      
      // add source, S, if present
      if (S != NULL){
        ierr = xf_Error(xf_SetVector(S, xfe_Add, R));
        if (ierr != xf_OK) return ierr;
      }
      
      // compute residual norms
      ierr = xf_Error(xf_VectorNorm(R, 1, &SolverData->ResNorm));
      if (ierr != xf_OK) return ierr;
      
      SolverData->RnormPrev = SolverData->Rnorm;
      ierr = xf_Error(xf_VectorNorm(R, 2, &SolverData->Rnorm));
      if (ierr != xf_OK) return ierr;
      //store residual l2-norm for SER
      if (SolverData->iIter == iIter0){ 
        SolverData->RnormPrev = SolverData->Rnorm;
        SolverData->CFLPrev = SolverData->CFL;
      }
      
      // On first iteration, make ResidualTolerance relative if LinearFlag is True
      if ((LinearFlag) && (SolverData->iIter==iIter0))
        ResidualTolerance = max(SolverData->ResNorm*ResidualTolerance, MEPS);
      
      // write log entry (also print to stdout)
      ierr = xf_Error(xf_WriteLogEntry(All, SolverData, U));
      if (ierr != xf_OK) return ierr;
      
      // convergence check (note this is a relative check if LinearFlag is True)
      if (SolverData->ResNorm < ResidualTolerance){
        Converged = xfe_True;
        break;
      }
      
      if (SolverData->iIter == iIter0){
        ierr = xf_Error(xf_SetZeroVector(dU));
        if (ierr != xf_OK) return ierr;
      }
      
      // solve linear system: R_U*dU = -R with appropriate linear solver
      ierr = xf_SolveLinearSystem(All, R_U0, R, xfe_False, -1, SolverData, dU);
      SolverError = xf_CheckSolverError(ierr, xfe_True, SolverData, USafe, U, &Rewind);
      if (Rewind) continue;
      
      // (U + dU) may be non-physical; check what fraction of dU we can add
      ierr = xf_Error(xf_UpdateState(All, U, dU, &LimitFlag, &UpdateFlag, SolverData));
      if (ierr != xf_OK) return ierr;
      
      SolverError = xf_CheckSolverError( ((UpdateFlag) ? xf_OK : xf_NO_UPDATE), xfe_True, 
                                        SolverData, USafe, U, &Rewind);
      if (Rewind) continue;
      
      // adjust CFL and store USafe if update was not limited
      if (!LimitFlag){
        SolverData->CFLSafe = SolverData->CFL;
        ierr = xf_Error(xf_SetVector(U, xfe_Set, USafe));
        if (ierr != xf_OK) return ierr;
      }
      //evolve the CFL number
      ierr = xf_Error(xf_EvolveCFL(SolverData, LimitFlag));
      if (ierr != xf_OK) return ierr;
    } // end while UpdateFlag == False
    
    if (Converged) break; // case converged
    if (Halted) break; // indicates user halt
    
    /* Write out U to hard disk if at requested interval */
    if ((SolverData->iIter > 0) && (WriteInterval > 0) &&
        ((SolverData->iIter % WriteInterval) == 0)){
      sprintf(OutputFile, "%s_Iter%d.data\0", SavePrefix, SolverData->iIter);
      ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "State", U, OutputFile));
      if (ierr != xf_OK) return ierr;
    }
    
  } // iiter
  
  if ((Converged) && (Verbosity != xfe_VerbosityLow))
    xf_printf("Nonlinear solver converged to tolerance.\n");
  
  
  // store parameters
  ierr = xf_SetKeyValueReal(All->Param->KeyValue, "CFL", SolverData->CFL);
  if (ierr != xf_OK) return ierr;
  ierr = xf_SetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", SolverData->iIter);
  if (ierr != xf_OK) return ierr;
  
  // do not release temporary vectors
  
  // destroy SolverData
  if (pSolverData == NULL){
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    if (ierr != xf_OK) return ierr;
  }
  
  
  return xf_OK;
}


#include "xf_SolverPMG.c"

/******************************************************************/
//   FUNCTION Definition: xf_SolveNonlinearSystem
int
xf_SolveNonlinearSystem(xf_All *All, real c, enum xfe_Bool LinearFlag,
                        xf_Vector *S, xf_Vector **pU)
{
  int ierr, terr, iAdaptRobust;
  enum xfe_NonlinearSolverType NonlinearSolver;
  enum xfe_Bool Adapted;
  xf_Data *D;
  xf_Vector *UInitial;
  xf_SolverData *SolverData;
  
  //Make copy of initial condition
  ierr = xf_Error(xf_DuplicateVector(All->Mesh, (*pU), &UInitial));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DataSetRemove(All->DataSet, "UInitial", xfe_False));
  if (ierr != xf_OK) return ierr;
  
  //Adding to the data structure so it gets transferred whe adaptation happens.
  ierr = xf_Error(xf_DataSetAdd(All->DataSet, "UInitial", xfe_Vector, xfe_True, 
                                (void *) UInitial, NULL));
  if (ierr != xf_OK) return ierr;
  
  // get solver type
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "NonlinearSolver", 
                                     xfe_NonlinearSolverName, (int ) xfe_NonlinearSolverLast, 
                                     (int *) &NonlinearSolver));
  if (ierr != xf_OK) return ierr;
  
  Adapted = xfe_True;
  
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "iAdaptRobust", &iAdaptRobust));
  if (ierr != xf_OK) return ierr;
  
  while (Adapted) {
    // create SolverData
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    if (ierr != xf_OK) return ierr;
    
    Adapted = xfe_False;
    switch (NonlinearSolver){
      case xfe_NonlinearSolverNewton:
        terr = xf_Error(xf_SolveNonlinearSystem_Newton(All, c, LinearFlag, xfe_False, S, 
                                                       (*pU), &SolverData));
        break;
      case xfe_NonlinearSolverpMultigrid:
        terr = xf_Error(xf_SolveNonlinearSystem_pMultigrid(All, c, LinearFlag, S, (*pU)));
        break;
      default:
        return xf_Error(xf_NOT_SUPPORTED);
        break;
    }
    if (terr == xf_NEED_ADAPT){
      ierr = xf_Error(xf_AdaptRobust(All, (*pU), &iAdaptRobust, &Adapted, SolverData));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "UInitial", xfe_Vector, &D));
      if (ierr != xf_OK) return ierr;
      UInitial = (xf_Vector *) D->Data;
      
      ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "State", xfe_Vector, &D));
      if (ierr != xf_OK) return ierr;
      (*pU) = (xf_Vector *) D->Data;
      
      ierr = xf_Error(xf_SetVector(UInitial, xfe_Set, (*pU)));
      if (ierr != xf_OK) return ierr;
      
    }
    else if (terr != xf_OK)return terr;
    
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    if (ierr != xf_OK) return ierr;
  }
  
  ierr = xf_Error(xf_DataSetRemove(All->DataSet, "UInitial", xfe_True));
  if (ierr != xf_OK) return ierr;
  
  UInitial = NULL;
  
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "iAdaptRobust", iAdaptRobust));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SolveAdjoints
int
xf_SolveAdjoints(xf_All *All, real c, real d, enum xfe_Bool ReuseJacobian,
                 xf_Vector *U, int nPsi, xf_Vector **S, xf_Vector **Psi,
                 xf_Vector **LinR, enum xfe_Bool CalcSolError, 
                 enum xfe_Bool SetZero)
{
  int ierr, nIter, iAdjoint;
  enum xfe_Bool WriteAdjoint, Found, LeanFlag;
  enum xfe_PreconditionerType Preconditioner;
  enum xfe_Verbosity Verbosity;
  char *OutputName;
  real LinResTol, J;
  xf_SolverData *SolverData;
  xf_Vector *Adj, *AdjR, *dAdj;
  xf_Vector *GCL=NULL, *AdjGCL=NULL, *AdjRGCL=NULL;
  xf_JacobianMatrix *R_U;
  xf_Data *D;
  xf_Output *Output;
  
  // check input args
  if (nPsi <= 0) return xf_Error(xf_INPUT_ERROR);
  
  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  /****************************************************************************/
  // Loop over adjoints to check if there is a combined output
  if (CalcSolError){
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      
      // current adjoint vector
      OutputName = Psi[iAdjoint]->OutputName;     
      
      ierr = xf_Error(xf_FindOutput(All->EqnSet, OutputName, &Output));
      if (ierr != xf_OK) return ierr;
      
      if (Output->nSumOutput > 1){
        if (Verbosity != xfe_VerbosityLow)
          xf_printf("Solving for the signs in weights in %s Output:\n", OutputName);
        
        ierr = xf_Error(xf_ErrEstSolution(All, U, OutputName, xfe_True, xfe_False));
        if (ierr != xf_OK) return ierr;
      }
    }
  }
  /****************************************************************************/
  
  // should we make the adjoint vector(s) writable?
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "WriteAdjoint", &WriteAdjoint);
  if (ierr != xf_OK) return ierr;
  
  // pull off (max) number of adjoint iterations
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "nIterAdjoint", &nIter);
  if (ierr != xf_OK) return ierr;
  
  // locate Residual vector for adjoint (same name as primal residual; should never conflict)
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, xfe_True, NULL, &AdjR, &Found));
  if (ierr != xf_OK) return ierr;  

  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;  
  
  
  // locate Adjoint Update vector, dAdj -- do not add to All
  ierr = xf_Error(xf_FindSimilarVector(All, U, "dAdj", xfe_True, xfe_False, NULL, &dAdj, NULL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetZeroVector(dAdj));
  if (ierr != xf_OK) return ierr;
  
  // locate preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
                                     xfe_PreconditionerName, (int ) xfe_PreconditionerLast, 
                                     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;
  
  // Is the preconditioner memory-lean?
  ierr = xf_Error(xf_PreconditionerLeanCheck(Preconditioner, &LeanFlag));
  if (ierr != xf_OK) return ierr;
  
  // pull off adjoint residual tolerance (relative)
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "AdjointResidualTolerance", 
                            &LinResTol);
  if (ierr != xf_OK) return ierr;
  SolverData->LinResTol = LinResTol;
  
  
  // locate Jacobian and aux vectors; check size; create if necessary
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        !LeanFlag, NULL, &R_U, &Found));
  if (ierr != xf_OK) return ierr;
  
  // Re-compute Jacobian R_U if does not exist or if not working with a linear primal system
  if (((!Found) || (!ReuseJacobian)) && (!LeanFlag)){
    ierr = xf_Error(xf_CalculateResidual(All, U, AdjR, R_U, SolverData));
    if (ierr != xf_OK){
      xf_printf("Error during calculation of Jacobian, R_U, for adjoint solve(s).\n");
      return ierr;
    }
    
    // Add c*(Mass matrix) to R_U to accommodate unsteady-mode
    ierr = xf_Error(xf_AddMassMatrix(All, c, NULL, U, NULL, R_U, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  
  // Loop over adjoints
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    
    // current adjoint vector
    Adj = Psi[iAdjoint];
    if (SetZero) {
      ierr = xf_Error(xf_SetZeroVector(Adj));
      if (ierr != xf_OK) return ierr;
    }
    
    // associated output name (make sure not NULL)
    if ((OutputName = Adj->OutputName) == NULL) return xf_Error(xf_INPUT_ERROR);
    
    if (Verbosity != xfe_VerbosityLow)
      xf_printf("Solving for %s Adjoint:\n", OutputName);
    
    // calculate output, J, and linearization, J_U (x d, added to AdjR)
    if (d != 0.0){
      J = 0.0;
      ierr = xf_Error(xf_CalculateOutput(All, OutputName, U, &J, AdjR, xfe_Set));
      if (ierr != xf_OK) return ierr;
      if (d != 1.0){
        ierr = xf_Error(xf_VectorMult(AdjR, d));
        if (ierr != xf_OK) return ierr;
      }
    }
    else{
      ierr = xf_Error(xf_SetZeroVector(AdjR));
      if (ierr != xf_OK) return ierr;
    }
    
    // add source to AdjR
    if (S != NULL){
      ierr = xf_Error(xf_SetVector(S[iAdjoint], xfe_Add, AdjR));
      if (ierr != xf_OK) return ierr;
    }
    
    // compute adjoint residual and store it in AdjR (AdjR gets incremented by R_U^T*Adj)
    ierr = xf_Error(xf_CalculateLinearResidual(All, R_U, Adj, NULL, xfe_True,
                                               SolverData, AdjR));
    if (ierr != xf_OK) return ierr;
    
    // store pointer to AdjR if desired
    if (LinR != NULL) LinR[iAdjoint] = AdjR;
    
    // solve linear system: R_U^T*dAdj + AdjR = 0 with appropriate linear solver
    SolverData->ReusedU = xfe_False;
    ierr = xf_Error(xf_SolveLinearSystem(All, R_U, AdjR, xfe_True, nIter, SolverData, dAdj));
    if (ierr != xf_OK)
      xf_printf("Warning: error during adjoint linear solve for %s adjoint. Continuing\n",
                OutputName);
    
    // set Adj = Adj + dAdj
    ierr = xf_Error(xf_SetVector(dAdj, xfe_Add, Adj));
    if (ierr != xf_OK) return ierr;
    
  } // iAdjoint
  
  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;
  
  // destroy dAdj
  ierr = xf_Error(xf_DestroyVector(dAdj, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_TakeExplicitStep
static int
xf_TakeExplicitStep(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                    real Time, real EndTime, real *pTimeStep, 
		    xf_SolverData *SolverData, xf_Vector **Ui, xf_Vector *R, 
		    xf_Vector *Unp1, enum xfe_Bool *pAtEndTime)
{
  /*
   PURPOSE:
   
   Takes a step of an explicit solver.  Stores result in the vector
   Unp1 = "U at n+1".  For example, forward Euler solves
   
   M/dt*(U^{n+1} - U^n) + R(U^n) = 0
   
   INPUTS:
   
   All         : all structure
   TimeScheme  : what time scheme to use
   Time        : current time
   pTimeStep   : pointer to delta t; this function is allowed to change delta t
   SolverData  : solver data structure (contains CFL)
   Ui          : set of state vectors at [n, n-1, ...]
   R           : residual vector (so don't have to locate it again)
   
   OUTPUTS: 
   
   Unp1       : updated state
   *pTimeStep : calculated time step (if stepping by CFL, for example)
   *pAtEndTime : used only when marching by CFL.  Flag to indicate that end of time
                 was reached.
   
   RETURN: Error code
   
   */
  int ierr, iRewind;
  enum xfe_Bool LocalTimeStepping, UseConstDt, Rewind, UseGCL = xfe_False;
  enum xfe_Verbosity Verbosity;
  real ResNorm, TimeStep, *pDt;
  xf_Vector *U, **F = NULL;
  xf_Vector *dtVec = NULL, *dtVec0 = NULL;
  xf_Vector *GCL = NULL, *GCL0 = NULL, **FGCL = NULL;

  char title[xf_MAXSTRLEN];

  int i, j, nstage;
  real *tstage;
  real *astage, a;

  // Forward Euler (FE) coefficients
  int nstage_FE = 1;
  real tstage_FE[1] = {0};
  real astage_FE[1] = {1};

  // Second-order Runge-Kutta (RK2) coefficients
  int nstage_RK2 = 2;
  real tstage_RK2[2] = {0, 0.5};
  real astage_RK2[4] = {0.5, 0,
		          0, 1};

  // Fourth-order Runge-Kutta (RK4) coefficients
  int nstage_RK4 = 4;
  real tstage_RK4[4] = {0, 0.5, 0.5, 1}; 
  real astage_RK4[16] = {  0.5,     0,     0,     0,
			     0,   0.5,     0,     0,
			     0,     0,   1.0,     0,
			 1./6., 1./3., 1./3., 1./6.};


  // set coefficients depending on the time scheme
  switch (TimeScheme){
  case xfe_TimeSchemeFE:
    nstage = nstage_FE ; tstage = tstage_FE ; astage = astage_FE ; break;
  case xfe_TimeSchemeRK2:
    nstage = nstage_RK2; tstage = tstage_RK2; astage = astage_RK2; break;
  case xfe_TimeSchemeRK4:
    nstage = nstage_RK4; tstage = tstage_RK4; astage = astage_RK4; break;
  default:
    return xf_Error(xf_NOT_SUPPORTED); break;
  }
   
  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    ierr = xf_Error(xf_InitMeshMotionGCLVector(All, Ui[0], -1, xfe_False, &GCL));
    if (ierr != xf_OK) return ierr;
  }

  // set TimeStep to the input value
  TimeStep = (*pTimeStep);


  // determine if using a constant time step
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseConstDt", 
                                     &UseConstDt));
  if (ierr != xf_OK) return ierr;
  
  // determine if stepping locally
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "LocalTimeStepping", 
                                     &LocalTimeStepping));
  if (ierr != xf_OK) return ierr;
  
  // can only step locally if not using a constant time step
  LocalTimeStepping = ((LocalTimeStepping) && (!UseConstDt));

  if ((LocalTimeStepping) || (!UseConstDt)){    
    // locate time step vector
    ierr = xf_Error(xf_FindVector(All, "TimeStep", xfe_LinkageGlobElem, 1, NULL, 0, 0, NULL, 
                                  NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  xfe_True, NULL, 
                                  &dtVec0, NULL));
    if (ierr != xf_OK) return ierr;
    dtVec = dtVec0;
  }
  
  // locate temporary vectors
  ierr = xf_Error(xf_Alloc( (void **) &F, nstage, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nstage; i++){
    sprintf(title, "ExplicitTimeStepF_%d", i);
    ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], title, xfe_True, xfe_True, NULL, F+i, NULL));
    if (ierr != xf_OK) return ierr;
  }

  // locate temporary vectors for GCL
  if (UseGCL){
    sprintf(title, "ExplicitTimeStepGCL0_%d", i);
    ierr = xf_Error(xf_FindSimilarVector(All, GCL, title, xfe_True, xfe_True, NULL, &GCL0, NULL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &FGCL, nstage, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nstage; i++){
      sprintf(title, "ExplicitTimeStepFGCL_%d", i);
      ierr = xf_Error(xf_FindSimilarVector(All, GCL, title, xfe_True, xfe_True, NULL, FGCL+i, NULL));
      if (ierr != xf_OK) return ierr;
    }
  }
    

  // Rewinds are a safety net for failed residual calculations; TimeStep or CFL is lowered.
  Rewind  = xfe_True;
  iRewind = 0;
  while (Rewind){
    Rewind = xfe_False;
    if (iRewind > 0) xf_printf("Warning, explicit step rewinding.\n");

    // calculate time step vector
    if ((LocalTimeStepping) || (!UseConstDt)){
      ierr = xf_CalculateArtificialTimeStep(All, Ui[0], SolverData->CFL, dtVec0);
      if (ierr != xf_OK){
	xf_printf("Calculation of time step failed.\n");
	if (LocalTimeStepping) xf_printf("Not stepping locally.\n");
	if (!UseConstDt) xf_printf("Using constant dt.\n");
	LocalTimeStepping = xfe_False;
	UseConstDt = xfe_True;
	dtVec = NULL;
      }
      else{
	if (LocalTimeStepping){
	  TimeStep = 1.0; // this is necessary for M^{-1} multiplications
	  dtVec = dtVec0;
	}
	else{ // using a non-constant dt
	  // compute minimum dt, store this as TimeStep
	  ierr = xf_Error(xf_VectorMin(dtVec0, NULL, &TimeStep));
	  if (ierr != xf_OK) return ierr;
	  // we don't want to go past EndTime
	  if (Time+TimeStep > EndTime){
	    TimeStep = EndTime - Time;
	    (*pAtEndTime) = xfe_True;
	  }
	  dtVec = NULL;
	}
      }
    }
    else SolverData->CFL = 1.0;

    iRewind++;

    // notify user what the currrent CFL and time step are
    if (Verbosity != xfe_VerbosityLow){
      if (LocalTimeStepping)
	xf_printf("%s Step: CFL = %.6E, Local time stepping\n", 
		  xfe_TimeSchemeName[TimeScheme], SolverData->CFL);
      else if (!UseConstDt)
	xf_printf("%s Step: CFL = %.6E, Global TimeStep = %.6E, Time = %.6E\n", 
		  xfe_TimeSchemeName[TimeScheme], SolverData->CFL, TimeStep, Time);
      else
	xf_printf("%s Step: Fixed TimeStep = %.6E, Time = %.6E\n", 
		  xfe_TimeSchemeName[TimeScheme], TimeStep, Time);
    }

    // make sure Time is set correctly
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
    if (ierr != xf_OK) return ierr;

    // pointer to local time step in case of rewinds
    pDt = ((dtVec == NULL) ? &TimeStep : NULL);

    // set U (pointer to Unp1) = Ui[0]
    U = Unp1;
    ierr = xf_Error(xf_SetVector(Ui[0], xfe_Set, U));
    if (ierr != xf_OK) return ierr;

    // if using GCL, set GCL0 = GCL
    if (UseGCL){
      ierr = xf_Error(xf_SetVector(GCL, xfe_Set, GCL0));
      if (ierr != xf_OK) return ierr;
    }

    // begin loop over stages
    for (i=0; i<nstage; i++){

      // set time using tstage
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time+TimeStep*tstage[i]));
      if (ierr != xf_OK) return ierr;
     
      // calculate the GCL residual if using a GCL
      // do this before state res calc because halo exchange of GCL goes on here
      if (UseGCL){
	ierr = xf_MeshMotionGCLResidual(All, GCL, xfe_True, FGCL[i]);
	if (ierr != xf_OK) return ierr;
      }

      // F[i] = residual(U)
      ierr = xf_CalculateResidual(All, U, F[i], NULL, SolverData);
      if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 0, NULL, 1, &Ui[0], pDt, &Rewind)) 
	return xf_Error(xf_SOLVER_ERROR);
      if (Rewind) break;
      // compute residual norm if on first stage
      if (i == 0){
	ierr = xf_Error(xf_VectorNorm(F[i], 1, &ResNorm));
	if (ierr != xf_OK) return ierr;
      }
            
      // F[i] = -dt*M^{-1}*F[i]
      ierr = xf_Error(xf_MultInvMassMatrix(All, -TimeStep, dtVec, F[i]));
      if (ierr != xf_OK) return ierr;
      if (UseGCL){
	// FGCL[i] = -dt*M^{-1}*FGCL[i]
	ierr = xf_Error(xf_MultInvMassMatrix(All, -TimeStep, dtVec, FGCL[i]));
	if (ierr != xf_OK) return ierr;
      }

      // U = U0 + sum_{j=0}^i astage(i,j) F(j) ... this sets U for the next stage
      // same steps for GCL
      ierr = xf_Error(xf_SetVector(Ui[0], xfe_Set, U));
      if (ierr != xf_OK) return ierr;
      if (UseGCL){
	ierr = xf_Error(xf_SetVector(GCL0, xfe_Set, GCL));
	if (ierr != xf_OK) return ierr;
      }
      for (j=0; j<=i; j++){
	if ((a=astage[i*nstage+j]) == 0.) continue;
	ierr = xf_Error(xf_VectorMultSet(F[j], a, xfe_Add, U));
	if(ierr != xf_OK) return ierr;
	if (UseGCL){
	  ierr = xf_Error(xf_VectorMultSet(FGCL[j], a, xfe_Add, GCL));
	  if(ierr != xf_OK) return ierr;
	}
      } // j

    } // end for i over stages

    // note, since U points to Unp1, after last stage, the answer is in the right place

    if (Rewind) continue;

  } // while Rewind
      
  // store residual norm for logging
  SolverData->ResNorm = ResNorm;

  // Reset Time
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
  if (ierr != xf_OK) return ierr;

  // evolve the CFL number
  if ((LocalTimeStepping) || (!UseConstDt)){
    ierr = xf_Error(xf_EvolveCFL(SolverData, xfe_False));
    if (ierr != xf_OK) return ierr;
  }

  // store TimeStep as a return value in case it was changed
  (*pTimeStep) = TimeStep;

  // release memory
  xf_Release( (void *) F);
  if (UseGCL) xf_Release( (void *) FGCL);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_StoreTimeHistData
static int
xf_StoreTimeHistData(xf_All *All, xf_TimeHistData *TimeHistData,
                     int iTime, xf_Vector *U)
{
  /*
   PURPOSE:
   
   Stores outputs into time history data at time index iTime
   
   INPUTS:
   
   All          : All data structure
   TimeHistData : time history data structure
   iTime        : time index
   U            : state used for calculating outputs
   
   OUTPUTS:
   
   None, outputs are computed and stored
   
   RETURN: Error code
   
   */
  int ierr, i;
  
  if (TimeHistData == NULL) return xf_OK; // nothing to do
  
  if (iTime >= TimeHistData->nTime){
    xf_printf("Error: TimeHistData index exceeded (iTime = %d vs. nTime = %d)\n",
              iTime, TimeHistData->nTime);
    return xf_Error(xf_OUT_OF_BOUNDS);
  }
  for (i=0; i<TimeHistData->nOutput; i++){
    ierr = xf_Error(xf_CalculateOutput(All, TimeHistData->OutputNames[i], U,
                                       &TimeHistData->OutputValues[i][iTime],
                                       NULL, xfe_Set));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCION Definition: xf_WriteTimeHist
int 
xf_WriteTimeHist(xf_TimeHistData *TimeHistData, const char *fname)
{
  int ierr, myRank;
  int iTime, i, nOutput;
  FILE *fid;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  // open file for writing (error handled in parallel)
  ierr = xf_Error(xf_fopen(fname, "w", &fid));
  if (ierr != xf_OK) return ierr;
  
  // only root writes
  if (myRank == 0){
    // write out header
    fprintf(fid, "%% Time history data\n");
    fprintf(fid, "%% Time schemes: \n");
    for (i=0; i<xfe_TimeSchemeLast; i++)
      fprintf(fid, "%%   %d = %s\n", i, xfe_TimeSchemeName[i]);
    fprintf(fid, "%% %16s %16s %16s %16s", "iTime", "Time", "TimeStep", "TimeScheme");
    nOutput = TimeHistData->nOutput;
    for (i=0; i<nOutput; i++)
      fprintf(fid, " %16s", TimeHistData->OutputNames[i]);
    fprintf(fid, "\n");
    
    for (iTime=0; iTime<TimeHistData->nTime; iTime++){
      fprintf(fid, "  %16d %16.8E %16.8E %16d", iTime, TimeHistData->Time[iTime],
              TimeHistData->TimeStep[iTime], TimeHistData->TimeScheme[iTime]);
      for (i=0; i<nOutput; i++)
        fprintf(fid, " %16.8E", TimeHistData->OutputValues[i][iTime]);
      fprintf(fid, "\n");
    }
  }
  
  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

#include "xf_SolverDGTime.c"
#include "xf_SolverIRK.c"

/******************************************************************/
//   FUNCTION Definition: xf_ApplyTimeScheme_MultiStepStage
static int
xf_ApplyTimeScheme_MultiStepStage(xf_All *All, const char *SavePrefix,
                                  enum xfe_Bool RestartFlag, xf_Vector **pU0,
                                  xf_TimeHistData *TimeHistData,
                                  xf_SolverData **pSolverData)
{
  int ierr, i, iStep, nTimeStep, iTime, myRank, nProc;
  int WriteInterval, RepartEvery;
  enum xfe_Bool LinearFlag, ConstTimeStep, ReuseJacobian;
  enum xfe_Bool Explicit, MultiStage, AtEndTime, UseGCL, ProcViz;
  enum xfe_TimeSchemeType TimeScheme;
  enum xfe_TimeSchemeType CurrentTimeScheme, PreviousTimeScheme;
  enum xfe_Verbosity Verbosity;
  char StateName[xf_MAXSTRLEN];
  char PreHeader[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  real Time, TimeStep, EndTime, CFLStart, c, J;
  xf_MultiStepData MSData;
  xf_Vector **Ui, *U0; // vector of all state vectors (at each time index)
  xf_Vector *S; // source vector
  xf_Vector *Rn; // residual vector
  xf_Vector *GCL  =NULL; // GCL vector
  xf_Vector **GCLi=NULL; // vector of all GCL vectors (at each time index)
  xf_Vector *SGCL =NULL; // GCL source vector
  xf_Vector *RnGCL=NULL; // GCL residual vector

  xf_GenArray *gA;
  xf_Data *D;
  xf_DataSet *DataSet = NULL;
  xf_SolverData *SolverData = (*pSolverData);
  U0 = (*pU0);
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
    
  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
                                     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
                                     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;
  
  // need TimeHistData
  if (TimeHistData == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // TimeHistData must be filled in
  nTimeStep = TimeHistData->nTime-1;
  if (nTimeStep <  0) return xf_Error(xf_INPUT_ERROR);
  if (nTimeStep == 0) return xf_OK; // nothing to do
  
  // is time step constant?
  ConstTimeStep = TimeHistData->ConstTimeStep;
  
  // Unsteady write interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 
                                    &WriteInterval));
  if (ierr != xf_OK) return ierr;
  
  // Repartition interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "DD_RepartEveryNdt", 
                                    &RepartEvery));
  if (ierr != xf_OK) return ierr;
  if (RepartEvery <= 0) //no repartitioning
    RepartEvery = nTimeStep+2;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ProcViz", 
                                     &ProcViz));
  if (ierr != xf_OK) return ierr;
  
  // Is the system linear (as prescribed by the user)?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "LinearFlag", &LinearFlag));
  if (ierr != xf_OK) return ierr;
  
  
  // Initial CFL for artificial time stepping
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFLStart));
  if (ierr != xf_OK) return ierr;

  // End Time (used if marching by CFL in explicit solver)
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "EndTime", &EndTime));
  if (ierr != xf_OK) return ierr;
  
  // determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U0, -1, xfe_False, &GCL));
    if (ierr != xf_OK) return ierr;
  }

  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  
  // allocate space for state vector pointers
  ierr = xf_Error(xf_Alloc( (void **) &Ui, MultiStepData[TimeScheme].nStep+1, 
                           sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    ierr = xf_Error(xf_Alloc( (void **) &GCLi, MultiStepData[TimeScheme].nStep+1, 
			      sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
  }
  
  
  // the first state vector is U0 (TimeIndex == 0)
  Ui[0] = U0;
  if (UseGCL) GCLi[0] = GCL;
  
  // locate required state vectors; create if necessary
  for (iStep = 1; iStep <= MultiStepData[TimeScheme].nStep; iStep++){
    ierr = xf_FindPrimalState(All->DataSet, iStep, &D, NULL);
    if (ierr == xf_OK){
      Ui[iStep] = (xf_Vector *) D->Data;
    }
    else if (ierr == xf_NOT_FOUND){
      sprintf(StateName, "State_%d", iStep);
      ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], StateName, xfe_True, xfe_True,  
                                           NULL, Ui + iStep, NULL));
      if (ierr != xf_OK) return ierr;
    }
    else return xf_Error(ierr);
    Ui[iStep]->SolverRole = xfe_SolverRolePrimalState;
    Ui[iStep]->TimeIndex = iStep;
  }
  // same for GCL
  if (UseGCL){
    for (iStep = 1; iStep <= MultiStepData[TimeScheme].nStep; iStep++){
      sprintf(StateName, "%s_%d", GCLVectorTitle, iStep);
      ierr = xf_Error(xf_FindSimilarVector(All, GCLi[0], StateName, xfe_True, xfe_True,  
                                           NULL, GCLi + iStep, NULL));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // locate source vector
  ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], "Source", xfe_True, xfe_True,  
                                       NULL, &S, NULL));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    ierr = xf_Error(xf_FindSimilarVector(All, GCLi[0], "SourceGCL", xfe_True, xfe_True,  
					 NULL, &SGCL, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  // locate Residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], "Residual", xfe_False, xfe_True, 
                                       NULL, &Rn, NULL));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    ierr = xf_Error(xf_FindSimilarVector(All, GCLi[0], "ResidualGCL", xfe_False, xfe_True, 
					 NULL, &RnGCL, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  // calculate outputs at t=0, store in TimeHistData
  ierr = xf_Error(xf_StoreTimeHistData(All, TimeHistData, 0, Ui[0])); 
  if (ierr != xf_OK) return ierr;

  
  /* Create a dataset for writing the state at each WriteInterval */
  if (SavePrefix != NULL){
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;
    // add state to dataset
    ierr = xf_Error(xf_DataSetAdd(DataSet, "State", xfe_Vector,
				  xfe_True, (void *) Ui[0], NULL));
    if (ierr != xf_OK) return ierr;
    // add GCL vector to dataset
    if (UseGCL){
      ierr = xf_Error(xf_DataSetAdd(DataSet, GCLVectorTitle, xfe_Vector,
				    xfe_True, (void *) GCL, NULL));
      if (ierr != xf_OK) return ierr;
      GCL->SolverRole = xfe_SolverRoleOther;
    }
  }

  /* Write out iTime == 0 vector (initial condition) */
  if (SavePrefix != NULL){
    sprintf(OutputFile, "%s_U%d.data\0", SavePrefix, 0);
    ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, OutputFile));
    if (ierr != xf_OK) return ierr;
  }


  // flag indicating that the end of the simulation time has been reached
  AtEndTime = xfe_False;
  
  // Will hold the time scheme used at previous step to detect change
  PreviousTimeScheme = xfe_TimeSchemeLast;
  
  // Begin loop over time steps
  for (iTime=1; iTime<=nTimeStep; iTime++){
    if (iTime%RepartEvery == 0 && nProc >1 && iTime > 1){
      xf_printf("Rebalancing CPU-work load...\n");fflush(stdout);
      //rebalance the cpu work load
      ierr = xf_Error(xf_CPULoadBalance(All, Ui, MultiStepData[TimeScheme].nStep+1));
      if (ierr != xf_OK) return ierr;
      
      //reset pointers in DataSet since the data was reparallelized in xf_CPULoadBalance
      DataSet->Head = NULL;
      DataSet->Tail = NULL;
      
      //find primal state and GCL vectors
      ierr = xf_Error(xf_DataSetAdd(DataSet, "State", xfe_Vector,
                                    xfe_True, (void *) Ui[0], NULL));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_InitMeshMotionGCLVector(All, Ui[0], -1, xfe_False, &GCL));
      if (ierr != xf_OK) return ierr;
      
      if (UseGCL){
        ierr = xf_Error(xf_DataSetAdd(DataSet, GCLVectorTitle, xfe_Vector,
                                      xfe_True, (void *) GCL, NULL));
        if (ierr != xf_OK) return ierr;
      }
      
      if (ProcViz) {
        ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "ProcID", xfe_Vector, &D));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_DataSetAdd(DataSet, D->Title, xfe_Vector, 
                                      xfe_True, D->Data, NULL));
        if (ierr != xf_OK) return ierr;
      }
      //reset the pointer
      ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D, NULL));
      if (ierr != xf_OK) return ierr;
      U0 = (xf_Vector *)D->Data;
      (*pU0) = (xf_Vector *)D->Data;
      
      
      xf_printf("done.\n");fflush(stdout);
    }
    
    // current time scheme, time, and time step
    CurrentTimeScheme = TimeHistData->TimeScheme[iTime];
    TimeStep          = TimeHistData->TimeStep[iTime];
    Time              = TimeHistData->Time[iTime];
    
    // Obtain multi-step data
    MSData = MultiStepData[CurrentTimeScheme];

    
    // Is Time Scheme explicit?
    ierr = xf_Error(xf_GetTimeSchemeInfo(CurrentTimeScheme, &Explicit, 
                                         NULL, &MultiStage));
    if (ierr != xf_OK) return ierr;
    
    
    if (Explicit){ // explicit multistep or multistage
      // If explicit, call separate function
      ierr = xf_Error(xf_TakeExplicitStep(All, CurrentTimeScheme,
                                          TimeHistData->Time[iTime-1], EndTime,
                                          &TimeStep, SolverData, Ui, Rn, S, &AtEndTime));
      if (ierr != xf_OK) return ierr;
      
      // Shift data in Ui vectors by a series of swaps
      for (iStep=MSData.nStep; iStep>0; iStep--)
        swap(Ui[iStep]->GenArray, Ui[iStep-1]->GenArray, gA);
      
      // set Ui[0] = S (the calculated state)
      ierr = xf_Error(xf_SetVector(S, xfe_Set, Ui[0]));
      if (ierr != xf_OK) return ierr;

      // The explicit step function is allowed to adjust the TimeStep
      Time = TimeHistData->Time[iTime] = TimeHistData->Time[iTime-1] + TimeStep;
      TimeHistData->TimeStep[iTime] = TimeStep;
      
      // Set Time (end of time step)
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
      if (ierr != xf_OK) return ierr;
      
    }
    else if (MultiStage){ // implicit multistage
      ierr = xf_Error(xf_TakeImplicitIRKStep(All, SavePrefix, RestartFlag,
					     TimeHistData->Time[iTime-1], TimeStep, 
					     SolverData, Ui, TimeScheme));

      if (ierr != xf_OK) return ierr;

    }
    else { // implicit multistep
      
      // construct S = c1*u^n + c2*u^{n-1} + ... 
      ierr = xf_Error(xf_VectorMultSet(Ui[0], MSData.alpha[1], xfe_Set, S));
      if (ierr != xf_OK) return ierr;
      for (iStep = 2; iStep <= MSData.nStep; iStep++){
        ierr = xf_Error(xf_VectorMultSet(Ui[iStep-1], MSData.alpha[iStep], xfe_Add, S));
        if (ierr != xf_OK) return ierr;
      }
      if (UseGCL){ // same with GCL
	ierr = xf_Error(xf_VectorMultSet(GCLi[0], MSData.alpha[1], xfe_Set, SGCL));
	if (ierr != xf_OK) return ierr;
	for (iStep = 2; iStep <= MSData.nStep; iStep++){
	  ierr = xf_Error(xf_VectorMultSet(GCLi[iStep-1], MSData.alpha[iStep], xfe_Add, SGCL));
	  if (ierr != xf_OK) return ierr;
	}
      }
      
      // Set S = M/dt * S
      ierr = xf_Error(xf_MultMassMatrix(All, 1.0/TimeStep, S));
      if (ierr != xf_OK) return ierr;
      if (UseGCL){
	ierr = xf_Error(xf_MultMassMatrix(All, 1.0/TimeStep, SGCL));
	if (ierr != xf_OK) return ierr;
      }
      
      // Add current residual (e.g. trapezoidal): S += beta * Rn
      if (MSData.beta != 0.0){      
        /* Calculate residual at (beginning of) current time step */
        ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", 
                                           TimeHistData->Time[iTime-1]));
        if (ierr != xf_OK) return ierr;
	// GCL residual calculation goes first
	if (UseGCL){
	  ierr = xf_Error(xf_MeshMotionGCLResidual(All, GCLi[0], xfe_True, RnGCL));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_VectorMultSet(RnGCL, MSData.beta, xfe_Add, SGCL));
	  if (ierr != xf_OK) return ierr;
	}
        ierr = xf_CalculateResidual(All, Ui[0], Rn, NULL, SolverData);
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_VectorMultSet(Rn, MSData.beta, xfe_Add, S));
        if (ierr != xf_OK) return ierr;
      }
      
      // Shift data in Ui vectors by a series of swaps
      for (iStep=MSData.nStep; iStep>0; iStep--)
        swap(Ui[iStep]->GenArray, Ui[iStep-1]->GenArray, gA);
      if (UseGCL)
	for (iStep=MSData.nStep; iStep>0; iStep--)
	  swap(GCLi[iStep]->GenArray, GCLi[iStep-1]->GenArray, gA);

      // keep Ui[0] = Ui[1] (the current state)    
      ierr = xf_Error(xf_SetVector(Ui[1], xfe_Set, Ui[0]));
      if (ierr != xf_OK) return ierr;
      if (UseGCL){
	ierr = xf_Error(xf_SetVector(GCLi[1], xfe_Set, GCLi[0]));
	if (ierr != xf_OK) return ierr;
      }
      
      // Set initial CFL
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart));
      if (ierr != xf_OK) return ierr;
      
      // Set Time (end of time step)
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
      if (ierr != xf_OK) return ierr;
      
      // Write out header for this time step
      sprintf(PreHeader, "%s %.15E", "%% Time = ", Time);
      ierr = xf_Error(xf_WriteLogHeader(All, PreHeader));
      if (ierr != xf_OK) return ierr;
      
      c = MSData.alpha[0]/TimeStep; // coefficient on M*u^{n+1}
      
      /* Can reuse R_U for linear problems provided dt is constant and
       the time scheme has not changed from previous iteration. */
      ReuseJacobian = (LinearFlag && TimeHistData->ConstTimeStep &&
                       (CurrentTimeScheme == PreviousTimeScheme) && (!UseGCL));
      
      /* Call non-linear solver.  Ui[0] contains most recent state,
       which will be used as the starting point for the nonlinear
       solve. */
      if (UseGCL){
	// GCL first; this is a trivial linear "solve" with the mass matrix
	ierr = xf_Error(xf_MeshMotionGCLSolveSystem(All, c, SGCL, GCLi[0]));
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_SolveNonlinearSystem(All, c, ReuseJacobian, S, Ui+0));
      if (ierr != xf_OK){
        xf_printf("Error: nonlinear solve failed at iTime = %d, Time = %.6E\n", iTime, Time);
        return ierr;
      }
    } // end else implicit
    
    
    /* Calculate any desired outputs for time history */
    ierr = xf_Error(xf_StoreTimeHistData(All, TimeHistData, iTime, Ui[0]));
    if (ierr != xf_OK) return ierr;
    
    /* Calculate unsteady outputs (that are logged).  Note that the
       time passed in is from the *beginning* of the time step. */
    ierr = xf_Error(xf_IncrementUnsteadyOutputs(All, NULL, CurrentTimeScheme, 1, 
                                                Ui, S, Time-TimeStep, TimeStep, (iTime==1), 
                                                (iTime==nTimeStep), &J, NULL));
    if (ierr != xf_OK) return ierr; 
    
    /* Write out Ui[0] to hard disk if at requested interval */
    if ((SavePrefix != NULL) && (WriteInterval > 0) && ((iTime % WriteInterval) == 0)){
      sprintf(OutputFile, "%s_U%d.data\0", SavePrefix, iTime);
      ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, OutputFile));
      if (ierr != xf_OK) return ierr;
    }

    // write log entry (also print to stdout)
    SolverData->iIter = iTime;
    ierr = xf_Error(xf_WriteLogEntry(All, SolverData, Ui[0]));
    if (ierr != xf_OK) return ierr;

    // break if at or past EndTime
    if ((Time >= EndTime) || (AtEndTime)){
      if (Verbosity != xfe_VerbosityLow)
	xf_printf("Time = %.4f >= EndTime = %.4f.  Stopping unsteady simulation.\n", Time, EndTime);
      TimeHistData->nTime = iTime+1;
      break;
    }

    // break out if user requests a halt
    if (xf_CheckUserHalt(NULL)) break;
    
    PreviousTimeScheme = CurrentTimeScheme;
    
  } // iTime

  // destroy DataSet
  if (SavePrefix != NULL){
    for (D=DataSet->Head; D != NULL; D=D->Next) D->Data = NULL; // do not touch the vectors
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;
  }

  // release memory
  xf_Release( (void *) Ui);
  xf_Release( (void *) GCLi);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ApplyTimeScheme
int
xf_ApplyTimeScheme(xf_All *All, const char *SavePrefix,
                   enum xfe_Bool RestartFlag, xf_Vector **pU0,
                   xf_TimeHistData *TimeHistData)
{
  int ierr;
  enum xfe_TimeSchemeType TimeScheme;
  char OutputFile[xf_MAXSTRLEN];  
  xf_SolverData *SolverData;
  
  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
                                     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
                                     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;
  
  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;
  
  // pull off CFL, CFLIncreaseFactor, CFLDecreaseFactor, etc.
  ierr = xf_Error(xf_FindCFLData(All->Param->KeyValue, SolverData));
  if (ierr != xf_OK) return ierr;
  
  // TimeScheme should not be steady
  if (TimeScheme == xfe_TimeSchemeSteady) return xf_Error(xf_INPUT_ERROR);
  
  // need a TimeHistData structure
  if (TimeHistData == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // Call appropriate time-marching scheme
  switch (TimeScheme){
    case xfe_TimeSchemeDG1:
    case xfe_TimeSchemeDG2:
      ierr = xf_Error(xf_ApplyTimeScheme_DG(All, SavePrefix, RestartFlag, pU0,
                                            TimeHistData, &SolverData));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_TimeSchemeBDF1:
    case xfe_TimeSchemeBDF2:
    case xfe_TimeSchemeBDF3:
    case xfe_TimeSchemeTrapezoidal:
    case xfe_TimeSchemeFE:
    case xfe_TimeSchemeRK4:
    case xfe_TimeSchemeRK2:
    case xfe_TimeSchemeESDIRK3:
    case xfe_TimeSchemeESDIRK4:
    case xfe_TimeSchemeESDIRK5:
    case xfe_TimeSchemeDIRK3:
    case xfe_TimeSchemeDIRK4:
      ierr = xf_Error(xf_ApplyTimeScheme_MultiStepStage(All, SavePrefix, RestartFlag, pU0, 
                                                        TimeHistData, &SolverData));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  // End of unsteady computation: print out unsteady outputs
  ierr = xf_Error(xf_LogUnsteadyOutputs(All));
  if (ierr != xf_OK) return ierr;
  
  // print out time history
  if (SavePrefix != NULL){
    sprintf(OutputFile, "%s_TimeHist.txt\0", SavePrefix);
    ierr = xf_Error(xf_WriteTimeHist(TimeHistData, OutputFile));
    if (ierr != xf_OK) return ierr;
  }
  
  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xfYu_TakeExplicitStep
static int
xfYu_TakeExplicitStep(xf_All *All, Yu_Model *Model, enum xfe_TimeSchemeType TimeScheme, 
                    real Time, real EndTime, real *pTimeStep, 
                    xf_SolverData *SolverData, xf_Vector **Ui, xf_Vector *R, 
                    xf_Vector *Unp1, enum xfe_Bool *pAtEndTime)
{
    /*
     PURPOSE:
     
     Takes a step of an explicit solver.  Stores result in the vector
     Unp1 = "U at n+1".  For example, forward Euler solves
     
     M/dt*(U^{n+1} - U^n) + R(U^n) = 0
     
     INPUTS:
     
     All         : all structure
     TimeScheme  : what time scheme to use
     Time        : current time
     pTimeStep   : pointer to delta t; this function is allowed to change delta t
     SolverData  : solver data structure (contains CFL)
     Ui          : set of state vectors at [n, n-1, ...]
     R           : residual vector (so don't have to locate it again)
     
     OUTPUTS: 
     
     Unp1       : updated state
     *pTimeStep : calculated time step (if stepping by CFL, for example)
     *pAtEndTime : used only when marching by CFL.  Flag to indicate that end of time
     was reached.
     
     RETURN: Error code
     
     */
    int ierr, iRewind;
    enum xfe_Bool LocalTimeStepping, UseConstDt, Rewind, UseGCL = xfe_False;
    enum xfe_Verbosity Verbosity;
    real ResNorm, TimeStep, *pDt;
    xf_Vector *U, **F = NULL;
    xf_Vector *dtVec = NULL, *dtVec0 = NULL;
    xf_Vector *GCL = NULL, *GCL0 = NULL, **FGCL = NULL;
    
    char title[xf_MAXSTRLEN];
    
    int i, j, nstage;
    real *tstage;
    real *astage, a;
    
    // Forward Euler (FE) coefficients
    int nstage_FE = 1;
    real tstage_FE[1] = {0};
    real astage_FE[1] = {1};
    
    // Second-order Runge-Kutta (RK2) coefficients
    int nstage_RK2 = 2;
    real tstage_RK2[2] = {0, 0.5};
    real astage_RK2[4] = {0.5, 0,
        0, 1};
    
    // Fourth-order Runge-Kutta (RK4) coefficients
    int nstage_RK4 = 4;
    real tstage_RK4[4] = {0, 0.5, 0.5, 1}; 
    real astage_RK4[16] = {  0.5,     0,     0,     0,
        0,   0.5,     0,     0,
        0,     0,   1.0,     0,
        1./6., 1./3., 1./3., 1./6.};
    
    
    // set coefficients depending on the time scheme
    switch (TimeScheme){
        case xfe_TimeSchemeFE:
            nstage = nstage_FE ; tstage = tstage_FE ; astage = astage_FE ; break;
        case xfe_TimeSchemeRK2:
            nstage = nstage_RK2; tstage = tstage_RK2; astage = astage_RK2; break;
        case xfe_TimeSchemeRK4:
            nstage = nstage_RK4; tstage = tstage_RK4; astage = astage_RK4; break;
        default:
            return xf_Error(xf_NOT_SUPPORTED); break;
    }
    
    // determine verbosity level
    ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                       xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                       (int *) &Verbosity));
    if (ierr != xf_OK) return ierr;
    
    // determine if using a Geometric Conservation Law
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
    if (ierr != xf_OK) return ierr;
    if (UseGCL){
        ierr = xf_Error(xf_InitMeshMotionGCLVector(All, Ui[0], -1, xfe_False, &GCL));
        if (ierr != xf_OK) return ierr;
    }
    
    // set TimeStep to the input value
    TimeStep = (*pTimeStep);
    
    
    // determine if using a constant time step
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseConstDt", 
                                       &UseConstDt));
    if (ierr != xf_OK) return ierr;
    
    // determine if stepping locally
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "LocalTimeStepping", 
                                       &LocalTimeStepping));
    if (ierr != xf_OK) return ierr;
    
    // can only step locally if not using a constant time step
    LocalTimeStepping = ((LocalTimeStepping) && (!UseConstDt));
    
    if ((LocalTimeStepping) || (!UseConstDt)){    
        // locate time step vector
        ierr = xf_Error(xf_FindVector(All, "TimeStep", xfe_LinkageGlobElem, 1, NULL, 0, 0, NULL, 
                                      NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  xfe_True, NULL, 
                                      &dtVec0, NULL));
        if (ierr != xf_OK) return ierr;
        dtVec = dtVec0;
    }
    
    // locate temporary vectors
    ierr = xf_Error(xf_Alloc( (void **) &F, nstage, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nstage; i++){
        sprintf(title, "ExplicitTimeStepF_%d", i);
        ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], title, xfe_True, xfe_True, NULL, F+i, NULL));
        if (ierr != xf_OK) return ierr;
    }
    
    // locate temporary vectors for GCL
    if (UseGCL){
        sprintf(title, "ExplicitTimeStepGCL0_%d", i);
        ierr = xf_Error(xf_FindSimilarVector(All, GCL, title, xfe_True, xfe_True, NULL, &GCL0, NULL));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_Alloc( (void **) &FGCL, nstage, sizeof(xf_Vector *)));
        if (ierr != xf_OK) return ierr;
        for (i=0; i<nstage; i++){
            sprintf(title, "ExplicitTimeStepFGCL_%d", i);
            ierr = xf_Error(xf_FindSimilarVector(All, GCL, title, xfe_True, xfe_True, NULL, FGCL+i, NULL));
            if (ierr != xf_OK) return ierr;
        }
    }
    
    
    // Rewinds are a safety net for failed residual calculations; TimeStep or CFL is lowered.
    Rewind  = xfe_True;
    iRewind = 0;
    while (Rewind){
        Rewind = xfe_False;
        if (iRewind > 0) xf_printf("Warning, explicit step rewinding.\n");
        
        // calculate time step vector
        if ((LocalTimeStepping) || (!UseConstDt)){
            ierr = xf_CalculateArtificialTimeStep(All, Ui[0], SolverData->CFL, dtVec0);
            if (ierr != xf_OK){
                xf_printf("Calculation of time step failed.\n");
                if (LocalTimeStepping) xf_printf("Not stepping locally.\n");
                if (!UseConstDt) xf_printf("Using constant dt.\n");
                LocalTimeStepping = xfe_False;
                UseConstDt = xfe_True;
                dtVec = NULL;
            }
            else{
                if (LocalTimeStepping){
                    TimeStep = 1.0; // this is necessary for M^{-1} multiplications
                    dtVec = dtVec0;
                }
                else{ // using a non-constant dt
                    // compute minimum dt, store this as TimeStep
                    ierr = xf_Error(xf_VectorMin(dtVec0, NULL, &TimeStep));
                    if (ierr != xf_OK) return ierr;
                    // we don't want to go past EndTime
                    if (Time+TimeStep > EndTime){
                        TimeStep = EndTime - Time;
                        (*pAtEndTime) = xfe_True;
                    }
                    dtVec = NULL;
                }
            }
        }
        else SolverData->CFL = 1.0;
        
        iRewind++;
        
        // notify user what the currrent CFL and time step are
        if (Verbosity != xfe_VerbosityLow){
            if (LocalTimeStepping)
                xf_printf("%s Step: CFL = %.6E, Local time stepping\n", 
                          xfe_TimeSchemeName[TimeScheme], SolverData->CFL);
            else if (!UseConstDt)
                xf_printf("%s Step: CFL = %.6E, Global TimeStep = %.6E, Time = %.6E\n", 
                          xfe_TimeSchemeName[TimeScheme], SolverData->CFL, TimeStep, Time);
            else
                xf_printf("%s Step: Fixed TimeStep = %.6E, Time = %.6E\n", 
                          xfe_TimeSchemeName[TimeScheme], TimeStep, Time);
        }
        
        // make sure Time is set correctly
        ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
        if (ierr != xf_OK) return ierr;
        
        // pointer to local time step in case of rewinds
        pDt = ((dtVec == NULL) ? &TimeStep : NULL);
        
        // set U (pointer to Unp1) = Ui[0]
        U = Unp1;
        ierr = xf_Error(xf_SetVector(Ui[0], xfe_Set, U));
        if (ierr != xf_OK) return ierr;
        
        // if using GCL, set GCL0 = GCL
        if (UseGCL){
            ierr = xf_Error(xf_SetVector(GCL, xfe_Set, GCL0));
            if (ierr != xf_OK) return ierr;
        }
        
        // begin loop over stages
        for (i=0; i<nstage; i++){
            
            // set time using tstage
            ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time+TimeStep*tstage[i]));
            if (ierr != xf_OK) return ierr;
            
            // calculate the GCL residual if using a GCL
            // do this before state res calc because halo exchange of GCL goes on here
            if (UseGCL){
                ierr = xf_MeshMotionGCLResidual(All, GCL, xfe_True, FGCL[i]);
                if (ierr != xf_OK) return ierr;
            }
            
            // F[i] = residual(U)
            ierr = xfYu_CalculateResidual(All, Model, U, F[i], NULL, SolverData);
            if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 0, NULL, 1, &Ui[0], pDt, &Rewind)) 
                return xf_Error(xf_SOLVER_ERROR);
            if (Rewind) break;
            // compute residual norm if on first stage
            if (i == 0){
                ierr = xf_Error(xf_VectorNorm(F[i], 1, &ResNorm));
                if (ierr != xf_OK) return ierr;
            }
            
            // F[i] = -dt*M^{-1}*F[i]
            ierr = xf_Error(xf_MultInvMassMatrix(All, -TimeStep, dtVec, F[i]));
            if (ierr != xf_OK) return ierr;
            if (UseGCL){
                // FGCL[i] = -dt*M^{-1}*FGCL[i]
                ierr = xf_Error(xf_MultInvMassMatrix(All, -TimeStep, dtVec, FGCL[i]));
                if (ierr != xf_OK) return ierr;
            }
            
            // U = U0 + sum_{j=0}^i astage(i,j) F(j) ... this sets U for the next stage
            // same steps for GCL
            ierr = xf_Error(xf_SetVector(Ui[0], xfe_Set, U));
            if (ierr != xf_OK) return ierr;
            if (UseGCL){
                ierr = xf_Error(xf_SetVector(GCL0, xfe_Set, GCL));
                if (ierr != xf_OK) return ierr;
            }
            for (j=0; j<=i; j++){
                if ((a=astage[i*nstage+j]) == 0.) continue;
                ierr = xf_Error(xf_VectorMultSet(F[j], a, xfe_Add, U));
                if(ierr != xf_OK) return ierr;
                if (UseGCL){
                    ierr = xf_Error(xf_VectorMultSet(FGCL[j], a, xfe_Add, GCL));
                    if(ierr != xf_OK) return ierr;
                }
            } // j
            
        } // end for i over stages
        
        // note, since U points to Unp1, after last stage, the answer is in the right place
        
        if (Rewind) continue;
        
    } // while Rewind
    
    // store residual norm for logging
    SolverData->ResNorm = ResNorm;
    
    // Reset Time
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
    if (ierr != xf_OK) return ierr;
    
    // evolve the CFL number
    if ((LocalTimeStepping) || (!UseConstDt)){
        ierr = xf_Error(xf_EvolveCFL(SolverData, xfe_False));
        if (ierr != xf_OK) return ierr;
    }
    
    // store TimeStep as a return value in case it was changed
    (*pTimeStep) = TimeStep;
    
    // release memory
    xf_Release( (void *) F);
    if (UseGCL) xf_Release( (void *) FGCL);
    
    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xfYu_ApplyTimeScheme_MultiStepStage
static int
xfYu_ApplyTimeScheme_MultiStepStage(xf_All *All, Yu_Model *Model, const char *SavePrefix,
                                    enum xfe_Bool RestartFlag, xf_Vector **pU0,
                                    xf_TimeHistData *TimeHistData,
                                    xf_SolverData **pSolverData)
{
    int ierr, i, iStep, nTimeStep, iTime, myRank, nProc;
    int WriteInterval, RepartEvery;
    enum xfe_Bool LinearFlag, ConstTimeStep, ReuseJacobian;
    enum xfe_Bool Explicit, MultiStage, AtEndTime, UseGCL, ProcViz;
    enum xfe_TimeSchemeType TimeScheme;
    enum xfe_TimeSchemeType CurrentTimeScheme, PreviousTimeScheme;
    enum xfe_Verbosity Verbosity;
    char StateName[xf_MAXSTRLEN];
    char PreHeader[xf_MAXSTRLEN];
    char OutputFile[xf_MAXSTRLEN];
    real Time, TimeStep, EndTime, CFLStart, c, J;
    xf_MultiStepData MSData;
    xf_Vector **Ui, *U0; // vector of all state vectors (at each time index)
    xf_Vector *S; // source vector
    xf_Vector *Rn; // residual vector
    xf_Vector *GCL  =NULL; // GCL vector
    xf_Vector **GCLi=NULL; // vector of all GCL vectors (at each time index)
    xf_Vector *SGCL =NULL; // GCL source vector
    xf_Vector *RnGCL=NULL; // GCL residual vector
    
    xf_GenArray *gA;
    xf_Data *D;
    xf_DataSet *DataSet = NULL;
    xf_SolverData *SolverData = (*pSolverData);
    U0 = (*pU0);
    ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
    if (ierr != xf_OK) return ierr;
    
    // Time scheme
    ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
                                       xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
                                       (int *) &TimeScheme));
    if (ierr != xf_OK) return ierr;
    
    // need TimeHistData
    if (TimeHistData == NULL) return xf_Error(xf_INPUT_ERROR);
    
    // TimeHistData must be filled in
    nTimeStep = TimeHistData->nTime-1;
    if (nTimeStep <  0) return xf_Error(xf_INPUT_ERROR);
    if (nTimeStep == 0) return xf_OK; // nothing to do
    
    // is time step constant?
    ConstTimeStep = TimeHistData->ConstTimeStep;
    
    // Unsteady write interval
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 
                                      &WriteInterval));
    if (ierr != xf_OK) return ierr;
    
    // Repartition interval
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "DD_RepartEveryNdt", 
                                      &RepartEvery));
    if (ierr != xf_OK) return ierr;
    if (RepartEvery <= 0) //no repartitioning
        RepartEvery = nTimeStep+2;
    
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ProcViz", 
                                       &ProcViz));
    if (ierr != xf_OK) return ierr;
    
    // Is the system linear (as prescribed by the user)?
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "LinearFlag", &LinearFlag));
    if (ierr != xf_OK) return ierr;
    
    
    // Initial CFL for artificial time stepping
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFLStart));
    if (ierr != xf_OK) return ierr;
    
    // End Time (used if marching by CFL in explicit solver)
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "EndTime", &EndTime));
    if (ierr != xf_OK) return ierr;
    
    // determine if using a Geometric Conservation Law
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
    if (ierr != xf_OK) return ierr;
    if (UseGCL){
        ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U0, -1, xfe_False, &GCL));
        if (ierr != xf_OK) return ierr;
    }
    
    // determine verbosity level
    ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                       xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                       (int *) &Verbosity));
    if (ierr != xf_OK) return ierr;
    
    
    // allocate space for state vector pointers
    ierr = xf_Error(xf_Alloc( (void **) &Ui, MultiStepData[TimeScheme].nStep+1, 
                             sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    if (UseGCL){
        ierr = xf_Error(xf_Alloc( (void **) &GCLi, MultiStepData[TimeScheme].nStep+1, 
                                 sizeof(xf_Vector *)));
        if (ierr != xf_OK) return ierr;
    }
    
    
    // the first state vector is U0 (TimeIndex == 0)
    Ui[0] = U0;
    if (UseGCL) GCLi[0] = GCL;
    
    // locate required state vectors; create if necessary
    for (iStep = 1; iStep <= MultiStepData[TimeScheme].nStep; iStep++){
        ierr = xf_FindPrimalState(All->DataSet, iStep, &D, NULL);
        if (ierr == xf_OK){
            Ui[iStep] = (xf_Vector *) D->Data;
        }
        else if (ierr == xf_NOT_FOUND){
            sprintf(StateName, "State_%d", iStep);
            ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], StateName, xfe_True, xfe_True,  
                                                 NULL, Ui + iStep, NULL));
            if (ierr != xf_OK) return ierr;
        }
        else return xf_Error(ierr);
        Ui[iStep]->SolverRole = xfe_SolverRolePrimalState;
        Ui[iStep]->TimeIndex = iStep;
    }
    // same for GCL
    if (UseGCL){
        for (iStep = 1; iStep <= MultiStepData[TimeScheme].nStep; iStep++){
            sprintf(StateName, "%s_%d", GCLVectorTitle, iStep);
            ierr = xf_Error(xf_FindSimilarVector(All, GCLi[0], StateName, xfe_True, xfe_True,  
                                                 NULL, GCLi + iStep, NULL));
            if (ierr != xf_OK) return ierr;
        }
    }
    
    // locate source vector
    ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], "Source", xfe_True, xfe_True,  
                                         NULL, &S, NULL));
    if (ierr != xf_OK) return ierr;
    if (UseGCL){
        ierr = xf_Error(xf_FindSimilarVector(All, GCLi[0], "SourceGCL", xfe_True, xfe_True,  
                                             NULL, &SGCL, NULL));
        if (ierr != xf_OK) return ierr;
    }
    
    // locate Residual vector
    ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], "Residual", xfe_False, xfe_True, 
                                         NULL, &Rn, NULL));
    if (ierr != xf_OK) return ierr;
    if (UseGCL){
        ierr = xf_Error(xf_FindSimilarVector(All, GCLi[0], "ResidualGCL", xfe_False, xfe_True, 
                                             NULL, &RnGCL, NULL));
        if (ierr != xf_OK) return ierr;
    }
    
    // calculate outputs at t=0, store in TimeHistData
    ierr = xf_Error(xf_StoreTimeHistData(All, TimeHistData, 0, Ui[0])); 
    if (ierr != xf_OK) return ierr;
    
    
    /* Create a dataset for writing the state at each WriteInterval */
    if (SavePrefix != NULL){
        ierr = xf_Error(xf_CreateDataSet(&DataSet));
        if (ierr != xf_OK) return ierr;
        // add state to dataset
        ierr = xf_Error(xf_DataSetAdd(DataSet, "State", xfe_Vector,
                                      xfe_True, (void *) Ui[0], NULL));
        if (ierr != xf_OK) return ierr;
        // add GCL vector to dataset
        if (UseGCL){
            ierr = xf_Error(xf_DataSetAdd(DataSet, GCLVectorTitle, xfe_Vector,
                                          xfe_True, (void *) GCL, NULL));
            if (ierr != xf_OK) return ierr;
            GCL->SolverRole = xfe_SolverRoleOther;
        }
    }
    
    /* Write out iTime == 0 vector (initial condition) */
    if (SavePrefix != NULL){
        sprintf(OutputFile, "%s_U%d.data\0", SavePrefix, 0);
        ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, OutputFile));
        if (ierr != xf_OK) return ierr;
    }
    
    
    // flag indicating that the end of the simulation time has been reached
    AtEndTime = xfe_False;
    
    // Will hold the time scheme used at previous step to detect change
    PreviousTimeScheme = xfe_TimeSchemeLast;
    
    // Begin loop over time steps
    for (iTime=1; iTime<=nTimeStep; iTime++){
        if (iTime%RepartEvery == 0 && nProc >1 && iTime > 1){
            xf_printf("Rebalancing CPU-work load...\n");fflush(stdout);
            //rebalance the cpu work load
            ierr = xf_Error(xf_CPULoadBalance(All, Ui, MultiStepData[TimeScheme].nStep+1));
            if (ierr != xf_OK) return ierr;
            
            //reset pointers in DataSet since the data was reparallelized in xf_CPULoadBalance
            DataSet->Head = NULL;
            DataSet->Tail = NULL;
            
            //find primal state and GCL vectors
            ierr = xf_Error(xf_DataSetAdd(DataSet, "State", xfe_Vector,
                                          xfe_True, (void *) Ui[0], NULL));
            if (ierr != xf_OK) return ierr;
            
            ierr = xf_Error(xf_InitMeshMotionGCLVector(All, Ui[0], -1, xfe_False, &GCL));
            if (ierr != xf_OK) return ierr;
            
            if (UseGCL){
                ierr = xf_Error(xf_DataSetAdd(DataSet, GCLVectorTitle, xfe_Vector,
                                              xfe_True, (void *) GCL, NULL));
                if (ierr != xf_OK) return ierr;
            }
            
            if (ProcViz) {
                ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "ProcID", xfe_Vector, &D));
                if (ierr != xf_OK) return ierr;
                ierr = xf_Error(xf_DataSetAdd(DataSet, D->Title, xfe_Vector, 
                                              xfe_True, D->Data, NULL));
                if (ierr != xf_OK) return ierr;
            }
            //reset the pointer
            ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D, NULL));
            if (ierr != xf_OK) return ierr;
            U0 = (xf_Vector *)D->Data;
            (*pU0) = (xf_Vector *)D->Data;
            
            
            xf_printf("done.\n");fflush(stdout);
        }
        
        // current time scheme, time, and time step
        CurrentTimeScheme = TimeHistData->TimeScheme[iTime];
        TimeStep          = TimeHistData->TimeStep[iTime];
        Time              = TimeHistData->Time[iTime];
        
        // Obtain multi-step data
        MSData = MultiStepData[CurrentTimeScheme];
        
        
        // Is Time Scheme explicit?
        ierr = xf_Error(xf_GetTimeSchemeInfo(CurrentTimeScheme, &Explicit, 
                                             NULL, &MultiStage));
        if (ierr != xf_OK) return ierr;
        
        
        if (Explicit){ // explicit multistep or multistage
            // If explicit, call separate function
            ierr = xf_Error(xfYu_TakeExplicitStep(All, Model, CurrentTimeScheme,
                                                TimeHistData->Time[iTime-1], EndTime,
                                                &TimeStep, SolverData, Ui, Rn, S, &AtEndTime));
            if (ierr != xf_OK) return ierr;
            
            // Shift data in Ui vectors by a series of swaps
            for (iStep=MSData.nStep; iStep>0; iStep--)
                swap(Ui[iStep]->GenArray, Ui[iStep-1]->GenArray, gA);
            
            // set Ui[0] = S (the calculated state)
            ierr = xf_Error(xf_SetVector(S, xfe_Set, Ui[0]));
            if (ierr != xf_OK) return ierr;
            
            // The explicit step function is allowed to adjust the TimeStep
            Time = TimeHistData->Time[iTime] = TimeHistData->Time[iTime-1] + TimeStep;
            TimeHistData->TimeStep[iTime] = TimeStep;
            
            // Set Time (end of time step)
            ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
            if (ierr != xf_OK) return ierr;
            
        }
        else if (MultiStage){ // implicit multistage
            ierr = xf_Error(xf_TakeImplicitIRKStep(All, SavePrefix, RestartFlag,
                                                   TimeHistData->Time[iTime-1], TimeStep, 
                                                   SolverData, Ui, TimeScheme));
            
            if (ierr != xf_OK) return ierr;
            
        }
        else { // implicit multistep
            
            // construct S = c1*u^n + c2*u^{n-1} + ... 
            ierr = xf_Error(xf_VectorMultSet(Ui[0], MSData.alpha[1], xfe_Set, S));
            if (ierr != xf_OK) return ierr;
            for (iStep = 2; iStep <= MSData.nStep; iStep++){
                ierr = xf_Error(xf_VectorMultSet(Ui[iStep-1], MSData.alpha[iStep], xfe_Add, S));
                if (ierr != xf_OK) return ierr;
            }
            if (UseGCL){ // same with GCL
                ierr = xf_Error(xf_VectorMultSet(GCLi[0], MSData.alpha[1], xfe_Set, SGCL));
                if (ierr != xf_OK) return ierr;
                for (iStep = 2; iStep <= MSData.nStep; iStep++){
                    ierr = xf_Error(xf_VectorMultSet(GCLi[iStep-1], MSData.alpha[iStep], xfe_Add, SGCL));
                    if (ierr != xf_OK) return ierr;
                }
            }
            
            // Set S = M/dt * S
            ierr = xf_Error(xf_MultMassMatrix(All, 1.0/TimeStep, S));
            if (ierr != xf_OK) return ierr;
            if (UseGCL){
                ierr = xf_Error(xf_MultMassMatrix(All, 1.0/TimeStep, SGCL));
                if (ierr != xf_OK) return ierr;
            }
            
            // Add current residual (e.g. trapezoidal): S += beta * Rn
            if (MSData.beta != 0.0){      
                /* Calculate residual at (beginning of) current time step */
                ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", 
                                                   TimeHistData->Time[iTime-1]));
                if (ierr != xf_OK) return ierr;
                // GCL residual calculation goes first
                if (UseGCL){
                    ierr = xf_Error(xf_MeshMotionGCLResidual(All, GCLi[0], xfe_True, RnGCL));
                    if (ierr != xf_OK) return ierr;
                    ierr = xf_Error(xf_VectorMultSet(RnGCL, MSData.beta, xfe_Add, SGCL));
                    if (ierr != xf_OK) return ierr;
                }
                //ierr = xf_CalculateResidual(All, Ui[0], Rn, NULL, SolverData);
                ierr = xfYu_CalculateResidual(All, &Model, Ui[0], Rn, NULL, SolverData); 
                if (ierr != xf_OK) return ierr;
                ierr = xf_Error(xf_VectorMultSet(Rn, MSData.beta, xfe_Add, S));
                if (ierr != xf_OK) return ierr;
            }
            
            // Shift data in Ui vectors by a series of swaps
            for (iStep=MSData.nStep; iStep>0; iStep--)
                swap(Ui[iStep]->GenArray, Ui[iStep-1]->GenArray, gA);
            if (UseGCL)
                for (iStep=MSData.nStep; iStep>0; iStep--)
                    swap(GCLi[iStep]->GenArray, GCLi[iStep-1]->GenArray, gA);
            
            // keep Ui[0] = Ui[1] (the current state)    
            ierr = xf_Error(xf_SetVector(Ui[1], xfe_Set, Ui[0]));
            if (ierr != xf_OK) return ierr;
            if (UseGCL){
                ierr = xf_Error(xf_SetVector(GCLi[1], xfe_Set, GCLi[0]));
                if (ierr != xf_OK) return ierr;
            }
            
            // Set initial CFL
            ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart));
            if (ierr != xf_OK) return ierr;
            
            // Set Time (end of time step)
            ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
            if (ierr != xf_OK) return ierr;
            
            // Write out header for this time step
            sprintf(PreHeader, "%s %.15E", "%% Time = ", Time);
            ierr = xf_Error(xf_WriteLogHeader(All, PreHeader));
            if (ierr != xf_OK) return ierr;
            
            c = MSData.alpha[0]/TimeStep; // coefficient on M*u^{n+1}
            
            /* Can reuse R_U for linear problems provided dt is constant and
             the time scheme has not changed from previous iteration. */
            ReuseJacobian = (LinearFlag && TimeHistData->ConstTimeStep &&
                             (CurrentTimeScheme == PreviousTimeScheme) && (!UseGCL));
            
            /* Call non-linear solver.  Ui[0] contains most recent state,
             which will be used as the starting point for the nonlinear
             solve. */
            if (UseGCL){
                // GCL first; this is a trivial linear "solve" with the mass matrix
                ierr = xf_Error(xf_MeshMotionGCLSolveSystem(All, c, SGCL, GCLi[0]));
                if (ierr != xf_OK) return ierr;
            }
            ierr = xf_Error(xf_SolveNonlinearSystem(All, c, ReuseJacobian, S, Ui+0));
            if (ierr != xf_OK){
                xf_printf("Error: nonlinear solve failed at iTime = %d, Time = %.6E\n", iTime, Time);
                return ierr;
            }
        } // end else implicit
        
        
        /* Calculate any desired outputs for time history */
        ierr = xf_Error(xf_StoreTimeHistData(All, TimeHistData, iTime, Ui[0]));
        if (ierr != xf_OK) return ierr;
        
        /* Calculate unsteady outputs (that are logged).  Note that the
         time passed in is from the *beginning* of the time step. */
        ierr = xf_Error(xf_IncrementUnsteadyOutputs(All, NULL, CurrentTimeScheme, 1, 
                                                    Ui, S, Time-TimeStep, TimeStep, (iTime==1), 
                                                    (iTime==nTimeStep), &J, NULL));
        if (ierr != xf_OK) return ierr; 
        
        /* Write out Ui[0] to hard disk if at requested interval */
        if ((SavePrefix != NULL) && (WriteInterval > 0) && ((iTime % WriteInterval) == 0)){
            sprintf(OutputFile, "%s_U%d.data\0", SavePrefix, iTime);
            ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, OutputFile));
            if (ierr != xf_OK) return ierr;
        }
        
        // write log entry (also print to stdout)
        SolverData->iIter = iTime;
        ierr = xf_Error(xf_WriteLogEntry(All, SolverData, Ui[0]));
        if (ierr != xf_OK) return ierr;
        
        // break if at or past EndTime
        if ((Time >= EndTime) || (AtEndTime)){
            if (Verbosity != xfe_VerbosityLow)
                xf_printf("Time = %.4f >= EndTime = %.4f.  Stopping unsteady simulation.\n", Time, EndTime);
            TimeHistData->nTime = iTime+1;
            break;
        }
        
        // break out if user requests a halt
        if (xf_CheckUserHalt(NULL)) break;
        
        PreviousTimeScheme = CurrentTimeScheme;
        
    } // iTime
    
    // destroy DataSet
    if (SavePrefix != NULL){
        for (D=DataSet->Head; D != NULL; D=D->Next) D->Data = NULL; // do not touch the vectors
        ierr = xf_Error(xf_DestroyDataSet(DataSet));
        if (ierr != xf_OK) return ierr;
    }
    
    // release memory
    xf_Release( (void *) Ui);
    xf_Release( (void *) GCLi);
    
    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xfYu_ApplyTimeScheme
static int
xfYu_ApplyTimeScheme(xf_All *All, Yu_Model *Model, const char *SavePrefix,
                   enum xfe_Bool RestartFlag, xf_Vector **pU0,
                   xf_TimeHistData *TimeHistData)
{
    int ierr;
    enum xfe_TimeSchemeType TimeScheme;
    char OutputFile[xf_MAXSTRLEN];  
    xf_SolverData *SolverData;
    
    // Time scheme
    ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
                                       xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
                                       (int *) &TimeScheme));
    if (ierr != xf_OK) return ierr;
    
    // create/allocate SolverData
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    if (ierr != xf_OK) return ierr;
    
    // pull off CFL, CFLIncreaseFactor, CFLDecreaseFactor, etc.
    ierr = xf_Error(xf_FindCFLData(All->Param->KeyValue, SolverData));
    if (ierr != xf_OK) return ierr;
    
    // TimeScheme should not be steady
    if (TimeScheme == xfe_TimeSchemeSteady) return xf_Error(xf_INPUT_ERROR);
    
    // need a TimeHistData structure
    if (TimeHistData == NULL) return xf_Error(xf_INPUT_ERROR);
    
    // Call appropriate time-marching scheme
    ierr = xf_Error(xfYu_ApplyTimeScheme_MultiStepStage(All, Model, SavePrefix, RestartFlag, pU0, 
                                                        TimeHistData, &SolverData));
    if (ierr != xf_OK) return ierr;

    
    // End of unsteady computation: print out unsteady outputs
    ierr = xf_Error(xf_LogUnsteadyOutputs(All));
    if (ierr != xf_OK) return ierr;
    
    // print out time history
    if (SavePrefix != NULL){
        sprintf(OutputFile, "%s_TimeHist.txt\0", SavePrefix);
        ierr = xf_Error(xf_WriteTimeHist(TimeHistData, OutputFile));
        if (ierr != xf_OK) return ierr;
    }
    
    // destroy SolverData
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    if (ierr != xf_OK) return ierr;
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyTimeSchemeAdjoint_MultiStep
static int
xf_ApplyTimeSchemeAdjoint_MultiStep(xf_All *All, const char *SavePrefix,
                                    xf_Vector *U0, int nPsi, xf_Vector **Psi, 
                                    xf_TimeHistData *TimeHistData)
{
  int ierr, i, j, iStep, MaxStep, nTime, iTime, jTime;
  int WriteInterval, iAdjoint, nSumWeight;
  enum xfe_Bool LinearFlag, ReuseJacobian, TimeSchemeChanged;
  enum xfe_Bool MotionOn, UseGCL;
  enum xfe_TimeSchemeType TimeScheme, jTimeScheme;
  enum xfe_TimeSchemeType CurrentTimeScheme;
  char Title[xf_MAXSTRLEN];
  char PreHeader[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  real *TimeWeights, **SumOutputWeights;
  real c, d, Time;
  real TimeStep[xf_MAXMULTISTEP];
  xf_MultiStepData MSData;
  xf_Vector ***Psii    = NULL; // vector of all adjoint vectors (at each time index)
  xf_Vector ***PsiiGCL = NULL; // similarly for GCL
  xf_Vector **PsiGCL   = NULL; // for storing current-time GCL adjoint vectors
  xf_Vector **S = NULL, **SGCL = NULL;
  xf_Vector *U = NULL, *GCL = NULL;
  xf_Output *Output;
  xf_Data *D, **DataPsi;
  xf_DataSet *DataSetPsi;
  xf_DataSet *DataSet = NULL;
  xf_GenArray *gA;
  
  // set convenient variables
  nTime = TimeHistData->nTime; // # of time points
  TimeWeights = TimeHistData->TimeWeights;
  nSumWeight = TimeHistData->nSumWeight;
  SumOutputWeights = TimeHistData->SumOutputWeights;
  
  // Is the system linear (as prescribed by the user)?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "LinearFlag", &LinearFlag));
  if (ierr != xf_OK) return ierr;
  
  // Default Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme",
                                     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast,
                                     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;
  
  // maximum number of steps at any given time point
  MaxStep = MultiStepData[TimeScheme].nStep;
  
  
  // Unsteady write interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 
                                    &WriteInterval));
  if (ierr != xf_OK) return ierr;

  // Are we doing mesh motion?
  MotionOn = ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));

  // determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    // need a GCL vector for evaluting conservative residuals, Jacobian, outputs
    ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U0, -1, xfe_False, &GCL));
    if (ierr != xf_OK) return ierr;
  }

  // cannot assume linear system if motion on or if GCL is on
  LinearFlag = (LinearFlag && (!UseGCL) & (!MotionOn));
  
  // Allocate vector of adjoint vectors
  ierr = xf_Error(xf_Alloc( (void **) &Psii, nPsi, sizeof(xf_Vector **)));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){ // same with GCL if using it
    ierr = xf_Error(xf_Alloc( (void **) &PsiiGCL, nPsi, sizeof(xf_Vector **)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &PsiGCL, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
  }
  
  
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    // allocate space for adjoint vector pointers
    ierr = xf_Error(xf_Alloc( (void **) (Psii+iAdjoint), MultiStepData[TimeScheme].nStep+1, 
                             sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    // the "initial" adjoint vector is in Psi
    Psii[iAdjoint][0] = Psi[iAdjoint];
    // locate required adjoint vectors; create if necessary
    for (iStep = 1; iStep <= MultiStepData[TimeScheme].nStep; iStep++){
      sprintf(Title, "%s_Adjoint_%d", Psi[iAdjoint]->OutputName, iStep);
      ierr = xf_Error(xf_FindSimilarVector(All, Psii[iAdjoint][0], Title, xfe_True, 
                                           xfe_True,  NULL, Psii[iAdjoint] + iStep, NULL));
      if (ierr != xf_OK) return ierr;
    }
    // same with GCL
    if (UseGCL){
      ierr = xf_Error(xf_Alloc( (void **) (PsiiGCL+iAdjoint), MultiStepData[TimeScheme].nStep+1, 
				sizeof(xf_Vector *)));
      if (ierr != xf_OK) return ierr;
      // find a GCL adjoint vector based on regular adjoint, Psi
      ierr = xf_Error(xf_InitMeshMotionGCLAdjoint(All, Psi[iAdjoint], -1, xfe_True, PsiiGCL[iAdjoint]));
      if (ierr != xf_OK) return ierr;
      PsiGCL[iAdjoint] = PsiiGCL[iAdjoint][0];
      // locate GCL adjoint vectors for the multistep/stage scheme
      for (iStep = 1; iStep <= MultiStepData[TimeScheme].nStep; iStep++){
	sprintf(Title, "%s_%s_%d", GCLAdjointTitle, Psi[iAdjoint]->OutputName, iStep);
	ierr = xf_Error(xf_FindSimilarVector(All, PsiiGCL[iAdjoint][0], Title, xfe_True, xfe_True,
					     NULL, PsiiGCL[iAdjoint] + iStep, NULL));
	if (ierr != xf_OK) return ierr;
      }
    }
  }
  
  // Allocate vector of source vectors
  ierr = xf_Error(xf_Alloc( (void **) &S, nPsi, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    ierr = xf_Error(xf_Alloc( (void **) &SGCL, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
  }
  
  // locate source vectors
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    sprintf(Title, "%s_AdjointSource", Psi[iAdjoint]->OutputName);
    ierr = xf_Error(xf_FindSimilarVector(All, Psii[iAdjoint][0], Title, xfe_False, 
                                         xfe_True, NULL, S + iAdjoint, NULL));
    if (ierr != xf_OK) return ierr;
    // initialize source to zero
    ierr = xf_Error(xf_SetZeroVector(S[iAdjoint]));
    if (ierr != xf_OK) return ierr;
  }
  if (UseGCL){
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      sprintf(Title, "%s_GCLAdjointSource", Psi[iAdjoint]->OutputName);
      ierr = xf_Error(xf_FindSimilarVector(All, PsiiGCL[iAdjoint][0], Title, xfe_True, 
					   xfe_True, NULL, SGCL + iAdjoint, NULL));
      if (ierr != xf_OK) return ierr;
      // initialize source to zero
      ierr = xf_Error(xf_SetZeroVector(SGCL[iAdjoint]));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  /* Create a dataset for writing the adjoint at each WriteInterval */
  ierr = xf_Error(xf_CreateDataSet(&DataSetPsi));
  if (ierr != xf_OK) return ierr;
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    sprintf(Title, "%s_Adjoint", Psi[iAdjoint]->OutputName);
    ierr = xf_Error(xf_DataSetAdd(DataSetPsi, Title, xfe_Vector, xfe_True, 
				  (void *) Psii[iAdjoint][0], NULL));
    if (ierr != xf_OK) return ierr;
    // add GCL vector to dataset
    if (UseGCL){
      sprintf(Title, "%s_GCLAdjoint", Psi[iAdjoint]->OutputName);
      ierr = xf_Error(xf_DataSetAdd(DataSetPsi, Title, xfe_Vector,
				    xfe_True, (void *) PsiiGCL[iAdjoint][0], NULL));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // Begin loop over time steps (backwards in time)
  for (iTime=nTime-1; iTime>0; iTime--){
    
    // break out if user requests a halt
    if (xf_CheckUserHalt(NULL)) break;
    
    CurrentTimeScheme = TimeHistData->TimeScheme[iTime];
    
    TimeSchemeChanged = ((iTime == nTime-1) || 
                         (CurrentTimeScheme != TimeHistData->TimeScheme[iTime+1]));
    
    MSData = MultiStepData[CurrentTimeScheme];
    
    // Set Time
    Time = TimeHistData->Time[iTime];
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
    if (ierr != xf_OK) return ierr;
    
    // Calculate time steps
    for (iStep=MSData.nStep; iStep>0; iStep--)
      TimeStep[iStep] = TimeStep[iStep-1];
    TimeStep[0] = TimeHistData->Time[iTime] - TimeHistData->Time[iTime-1];
    
    
    // Set weights on any SumOutputs; this allows for time-varying weights
    if ((nSumWeight > 0) && (SumOutputWeights != NULL)){
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
        ierr = xf_Error(xf_FindOutput(All->EqnSet, Psi[iAdjoint]->OutputName, &Output));
        if (ierr != xf_OK) return ierr;
        if (Output->Type == xfe_SumOutput){
          if (Output->nSumOutput != nSumWeight) return xf_Error(xf_OUT_OF_BOUNDS);
          for (i=0; i<nSumWeight; i++){
            Output->SumOutputWeights[i] = SumOutputWeights[iTime][i];
          }
        }
      } // iAdjoint
    }
    
    // Shift data in Psii vectors
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      for (iStep=MSData.nStep; iStep>0; iStep--)
        swap(Psii[iAdjoint][iStep]->GenArray, Psii[iAdjoint][iStep-1]->GenArray, gA);
      // keep Psii[iAdjoint][0] = Psii[iAdjoint][1] (the current state)    
      ierr = xf_Error(xf_SetVector(Psii[iAdjoint][1], xfe_Set, Psii[iAdjoint][0]));
      if (ierr != xf_OK) return ierr;
      if (UseGCL){
	for (iStep=MSData.nStep; iStep>0; iStep--)
	  swap(PsiiGCL[iAdjoint][iStep]->GenArray, PsiiGCL[iAdjoint][iStep-1]->GenArray, gA);
	// keep PsiiGCL[iAdjoint][0] = PsiiGCL[iAdjoint][1] (the current state)    
	ierr = xf_Error(xf_SetVector(PsiiGCL[iAdjoint][1], xfe_Set, PsiiGCL[iAdjoint][0]));
	if (ierr != xf_OK) return ierr;
      }
    }      

    // Read in U from hard-disk if problem is not linear
    if (LinearFlag)
      U = U0; // use the default input state
    else{
      // need SavePrefix
      if (SavePrefix == NULL) return xf_Error(xf_INPUT_ERROR);
      
      // create data set for reading states 
      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;
      // read .data from file
      sprintf(Title, "%s_U%d.data\0", SavePrefix, iTime);
      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, Title, DataSet));
      if (ierr != xf_OK) return ierr;
      
      // use first (usually only) piece of data
      D = DataSet->Head;
      U = (xf_Vector *) D->Data;
    }

    // also read in GCL vector for evaluating conservative state residuals
    if (UseGCL){
      // Read in GCL vectors from DataSet as well
      D = D->Next; // GCL comes after state
      ierr = xf_Error(xf_SetVector((xf_Vector *) D->Data, xfe_Set, GCL));
      if (ierr != xf_OK) return ierr;
    }    
    
    // The linear system to be inverted will be [c*M + R_U]
    c = MSData.alpha[0]/TimeStep[0];
        
    /* ReuseJacobian only if LinearFlag && ConstTimeStep && TimeScheme
       same as in previous iteration.  Need a constant time
       step/scheme because the inverted matrix is [alpha*M/dt + R_U],
       which changes if dt or alpha changes.
     */
    ReuseJacobian = (LinearFlag && TimeHistData->ConstTimeStep &&
                     (!TimeSchemeChanged));
    
    // Write out header for this time step
    sprintf(PreHeader, "%s %.15E", "%% Adjoint Solve: Time = ", Time);
    ierr = xf_Error(xf_WriteLogHeader(All, PreHeader));
    if (ierr != xf_OK) return ierr;
    
    // Set time-weighting of output
    if (TimeWeights == NULL) // NULL means final-time weight is 1.0; all other are zero
      d = ((iTime == (nTime-1)) ? 1.0 : 0.0);
    else d = TimeWeights[iTime];
    
    
    // Solve linear system with existing S (initialized to 0 before loop)
    ierr = xf_Error(xf_SolveAdjoints(All, c, d, ReuseJacobian, U, nPsi, S, Psi, NULL,xfe_False,xfe_False));
    if (ierr != xf_OK){
      xf_printf("Error: adjoint solve failed at iTime = %d, Time = %.6E\n", iTime, Time);
      return ierr;
    }
    if (UseGCL){
      // Add contribution of "U" adjoint to source for GCL adjoint solve 
      ierr = xf_Error(xf_MeshMotionGCL_PsiTxR_G(All, U, nPsi, Psi, SGCL));
      if (ierr != xf_OK) return ierr;
      // GCL adjoint solve 
      ierr = xf_Error(xf_MeshMotionGCLSolveAdjoints(All, c, d, nPsi, SGCL, PsiGCL));
      if (ierr != xf_OK) return ierr;
    }    
    
    /* write out Psi to hard disk if at requested interval (all
     adjoints are written to one dataset) */
    if ((SavePrefix != NULL) && (iTime % WriteInterval) == 0){
      sprintf(OutputFile, "%s_Psi%d.data\0", SavePrefix, iTime);
      ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSetPsi, NULL, OutputFile));
      if (ierr != xf_OK) return ierr;
    }
    
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      /*
       Carefully set source, S, using transpose of unsteady Jacobian:
       
       S = alpha1[TimeScheme(iTime  )]/dt(iTime  ) * Psi(iTime  )
         + alpha2[TimeScheme(iTime+1)]/dt(iTime+1) * Psi(iTime+1)
         + ...
       
       where the alphaj are the j'th coefficients in the respective
       time schemes.  In the sum, expect at most MaxStep terms;
       however, do not go beyond nTime-1.
      */
      
      ierr = xf_Error(xf_VectorMultSet(Psii[iAdjoint][0], MSData.alpha[1]/TimeStep[0], 
                                       xfe_Set, S[iAdjoint]));
      if (ierr != xf_OK) return ierr;
      if (UseGCL){
	ierr = xf_Error(xf_VectorMultSet(PsiiGCL[iAdjoint][0], MSData.alpha[1]/TimeStep[0], 
					 xfe_Set, SGCL[iAdjoint]));
	if (ierr != xf_OK) return ierr;
      }
      
      for (j=2; j <= min(MaxStep, nTime-iTime); j++){
        jTime = iTime + j-1;
        jTimeScheme = TimeHistData->TimeScheme[jTime];
        MSData = MultiStepData[jTimeScheme];
        if (MSData.nStep < j) continue;
        ierr = xf_Error(xf_VectorMultSet(Psii[iAdjoint][j-1], 
                                         MSData.alpha[j]/TimeStep[j-1], 
                                         xfe_Add, S[iAdjoint]));
        if (ierr != xf_OK) return ierr;
	if (UseGCL){
	  ierr = xf_Error(xf_VectorMultSet(PsiiGCL[iAdjoint][j-1], 
					   MSData.alpha[j]/TimeStep[j-1], 
					   xfe_Add, SGCL[iAdjoint]));
	  if (ierr != xf_OK) return ierr;
	}
      }
      
      // Set S = M^T * S  (using fact that the mass matrix, M, is symmetric)
      ierr = xf_Error(xf_MultMassMatrix(All, 1.0, S[iAdjoint]));
      if (ierr != xf_OK) return ierr;
      if (UseGCL){
	ierr = xf_Error(xf_MultMassMatrix(All, 1.0, SGCL[iAdjoint]));
	if (ierr != xf_OK) return ierr;
      }
    }
    
    // Destroy data set (and contained states) in the nonlinear case
    if (!LinearFlag){
      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;
    }
    
  } // iTime
  
  // Set Time to first value in history
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", TimeHistData->Time[0]));
  if (ierr != xf_OK) return ierr;
  
  
  /* Following loop, set Psi(t=0) = S.  Note that after the last step
     in the above iTime loop, S contains [R(i)_U(0)]^T*Psi(i), which
     is exactly what we want to return for Psi(t=0) -- see description
     in xf_Solver.h */
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    ierr = xf_Error(xf_SetVector(S[iAdjoint], xfe_Set, Psii[iAdjoint][0]));
    if (ierr != xf_OK) return ierr;
  }
  if (UseGCL){
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      ierr = xf_Error(xf_SetVector(SGCL[iAdjoint], xfe_Set, PsiiGCL[iAdjoint][0]));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // Destroy DataSetPsi
  for (D=DataSetPsi->Head; D != NULL; D=D->Next) D->Data = NULL; // do not touch the vectors
  ierr = xf_Error(xf_DestroyDataSet(DataSetPsi));
  if (ierr != xf_OK) return ierr;
  
  // release memory
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) xf_Release( (void *) Psii[iAdjoint]);
  if (UseGCL)
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) xf_Release( (void *) PsiiGCL[iAdjoint]);
  
  xf_Release( (void *) Psii);
  xf_Release( (void *) PsiiGCL);
  xf_Release( (void *) PsiGCL);
  xf_Release( (void *) S);
  xf_Release( (void *) SGCL);

  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ApplyTimeSchemeAdjoint
int
xf_ApplyTimeSchemeAdjoint(xf_All *All, const char *SavePrefix,
                          xf_Vector *U0, int nPsi, xf_Vector **Psi, 
                          xf_TimeHistData *TimeHistData)
{
  int ierr;
  enum xfe_TimeSchemeType TimeScheme;
  
  
  // need time history data
  if ((TimeHistData == NULL) || (TimeHistData->nTime < 1))
    return xf_Error(xf_INPUT_ERROR);
  
  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
                                     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
                                     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;
  
  
  // TimeScheme should not be steady
  if (TimeScheme == xfe_TimeSchemeSteady) return xf_Error(xf_INPUT_ERROR);
  
  // Call appropriate adjoint time-marching scheme
  switch (TimeScheme){
    case xfe_TimeSchemeDG1:
    case xfe_TimeSchemeDG2:
      ierr = xf_Error(xf_ApplyTimeSchemeAdjoint_DG(All, SavePrefix, U0,
                                                   nPsi, Psi, TimeHistData));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_TimeSchemeBDF1:
    case xfe_TimeSchemeBDF2:
    case xfe_TimeSchemeBDF3:
    case xfe_TimeSchemeTrapezoidal:
    case xfe_TimeSchemeFE:
    case xfe_TimeSchemeRK4:
    case xfe_TimeSchemeRK2:
      ierr = xf_Error(xf_ApplyTimeSchemeAdjoint_MultiStep(All, SavePrefix, 
                                                          U0, nPsi, Psi, 
                                                          TimeHistData));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ApplyTimeSchemeAdapt
int
xf_ApplyTimeSchemeAdapt(xf_All *All, const char *SavePrefix,
                        enum xfe_AdaptOnType AdaptOn,
                        xf_TimeHistData *TimeHistData)
{
  int ierr;
  enum xfe_TimeSchemeType TimeScheme;
  
  // need time history data
  if ((TimeHistData == NULL) || (TimeHistData->nTime < 1))
    return xf_Error(xf_INPUT_ERROR);
  
  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
                                     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
                                     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;
  
  
  // TimeScheme should not be steady
  if (TimeScheme == xfe_TimeSchemeSteady) return xf_Error(xf_INPUT_ERROR);
  
  // Call appropriate adaptation routine for each time-marching scheme
  switch (TimeScheme){
    case xfe_TimeSchemeDG1:
    case xfe_TimeSchemeDG2:
      ierr = xf_Error(xf_ApplyTimeSchemeAdapt_DG(All, SavePrefix, AdaptOn, TimeHistData));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_TimeSchemeBDF1:
    case xfe_TimeSchemeBDF2:
    case xfe_TimeSchemeBDF3:
    case xfe_TimeSchemeTrapezoidal:
    case xfe_TimeSchemeFE:
    case xfe_TimeSchemeRK4:
    case xfe_TimeSchemeRK2:
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_PingResidual
int
xf_PingResidual(xf_All *All, xf_Vector *U, real eps, real tol)
{
  int ierr, ieps, i, j, k, egrp, elem, face;
  int in, is, jn, js, sr;
  int egN, eN;
  int r, rM, rN, nFace;
  enum xfe_Verbosity Verbosity;
  real ep, **veps = NULL, rate, fd, fa, *R0, **R_U0;
  xf_Mesh *Mesh;
  xf_Vector *R;
  xf_JacobianMatrix *R_U;
  xf_SolverData *SolverData;
  
  
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;
  
  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  
  R0   = NULL;
  R_U0 = NULL;
  
  // locate Jacobian and aux vectors; check size; create if necessary
  // do not resize Jacobian if larger than necessary
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_True, NULL, &R_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, xfe_True, NULL, &R, NULL));
  if (ierr != xf_OK) return ierr;
  
  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;
  
  
  ierr = xf_CalculateResidual(All, U, R, R_U, SolverData);
  if (ierr != xf_OK) return ierr;
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      if (Verbosity != xfe_VerbosityLow)
        xf_printf("egrp=%d, elem=%d\n", egrp, elem);
      
      nFace = Mesh->ElemGroup[egrp].nFace[elem];
      r = ((U->GenArray[egrp].vr==NULL) ? U->GenArray[egrp].r : U->GenArray[egrp].vr[elem]);
      rM = r;
      for (face=0; face<nFace; face++){
        egN = R_U->egrpN[egrp][elem][face];
        eN  = R_U->elemN[egrp][elem][face];
        if (egN < 0) continue;
        rN = ((U->GenArray[egN].vr==NULL) ? U->GenArray[egN].r : U->GenArray[egN].vr[eN]);
        rM = max(rM, rN);
      }
      
      // calculate R and R_U
      ierr = xf_CalculateResidual(All, U, R, R_U, SolverData);
      if (ierr != xf_OK) return ierr;
      
      
      // store residual in R0 and R_U0
      ierr = xf_Error(xf_ReAlloc((void **) &R0, r, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      
      for (k=0; k<r; k++) R0[k] = R->GenArray[egrp].rValue[elem][k];
      
      ierr = xf_Error(xf_ReAlloc2((void ***) &R_U0, 1+nFace, r*rM, sizeof(real)));
      if (ierr != xf_OK) return ierr;

      // reallocate veps
      ierr = xf_Error(xf_ReAlloc2((void ***) &veps, r, 2, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      
      for (face = -1; face<nFace; face++){
        if (face == -1){
          egN = egrp;
          eN  = elem;
        }
        else{
          egN = R_U->egrpN[egrp][elem][face];
          eN  = R_U->elemN[egrp][elem][face];
        }
        if (egN < 0) {
          if (Verbosity != xfe_VerbosityLow) xf_printf("Boundary.\n");
          continue;
        }
        
        rN = ((U->GenArray[egN].vr==NULL) ? U->GenArray[egN].r : U->GenArray[egN].vr[eN]);
  
        for (k=0; k<r*rN; k++) R_U0[1+face][k] = R_U->Value[egrp][elem][1+face][k];
      } // face
      
      // loop over faces
      for (face = -1; face<nFace; face++){
        if (Verbosity != xfe_VerbosityLow) xf_printf("face = %d\n", face);
        if (face == -1){
          egN = egrp;
          eN  = elem;
        }
        else{
          egN = R_U->egrpN[egrp][elem][face];
          eN  = R_U->elemN[egrp][elem][face];
        }
        if (egN < 0) continue;
        
        rN = ((U->GenArray[egN].vr==NULL) ? U->GenArray[egN].r : U->GenArray[egN].vr[eN]);
        
        for (j=0; j<rN; j++){
          jn = j/sr;
          js = j%sr;
	  for (ieps=0; ieps<2; ieps++){ // loop over two epsilons
	    ep = eps / ((real) ieps + 1.0);
	    // increment U
	    U->GenArray[egN].rValue[eN][j] += ep;
          
	    // calculate R and R_U
	    ierr = xf_CalculateResidual(All, U, R, R_U, SolverData);
	    if (ierr != xf_OK) return ierr;
          
	    // check Delta R vs R_U
	    for (i=0; i<r; i++){
	      in = i/sr;
	      is = i%sr;	      
	      fd = (R->GenArray[egrp].rValue[elem][i] - R0[i])/ep;
	      fa = 0.5*(R_U->Value[egrp][elem][1+face][in*sr*rN + jn*sr*sr + is*sr+js] 
			+ R_U0[1+face][in*sr*rN + jn*sr*sr + is*sr+js]);
              if (Verbosity != xfe_VerbosityLow){
                xf_printf("[ieps=%d: %3d, %3d] fd=%20.8e, fa=%20.8e, diff=%20.8e\n", 
                          ieps, i, j, fd, fa, fabs(fd-fa));
              }
	      veps[i][ieps] = fabs(fd-fa);
	    } // i
          
	    // decrement U
	    U->GenArray[egN].rValue[eN][j] -= ep;
	  } // ieps
	  for (i=0; i<r; i++){
	    if (veps[i][0] <= 5e-8) continue; // too small to check rate
	    rate = log(veps[i][0]/veps[i][1])/log(2.0);
	    //xf_printf("[%d, %d] veps[0] = %.6e, veps[1] = %.6e, rate = %.6f\n",
	    //	      i, j, veps[i][0], veps[i][1], rate);
	    if ((Verbosity != xfe_VerbosityLow) ||
                ((tol > 0.) && (rate < 2.0-tol))){
              xf_printf("[%d, %d] rate=%.5f\n", i, j, rate);
            }
            if ((tol > 0.) && (rate < 2.0-tol)){
              xf_printf("[%d, %d] rate=%.5f (expected at least %.5f)\n", i, j, rate, 2.0-tol);
              return xf_Error(xf_PING_FAILED);
            }
	  } // i
        } // j
        
      } // face
      
    } // elem
  } // egrp
  
  
  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void  *) R0);
  xf_Release2((void **) R_U0);
  xf_Release2((void **) veps);
  
  return xf_OK;
  
}


#if( UNIT_TEST==1 )
#include "xf_Solver.test.in"
#endif
