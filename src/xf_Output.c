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
 FILE:  xf_Output.c
 
 This file contains functions for computing outputs.
 
 */

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_Memory.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Basis.h"
#include "xf_Quad.h"
#include "xf_QuadRule.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_String.h"
#include "xf_EqnSetHook.h"
#include "xf_MeshTools.h"
#include "xf_Residual.h"
#include "xf_SolverTools.h"
#include "xf_Solver.h"
#include "xf_MeshMotion.h"
#include "xf_MeshMotionGCL.h"
#include <time.h>



/******************************************************************/
//   FUNCTION Definition: xf_FindOutputGCLLinearization
static int
xf_FindOutputGCLLinearization(xf_All *All, const char *OutputName, 
                              enum xfe_AddType AddFlag, xf_Vector **pJ_GCL)
{
  int ierr;
  enum xfe_Bool UseGCL;
  
  // determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  
  // find J_GCL if it exists (e.g. for GCL adjoint solve)
  if (UseGCL){ 
    ierr = xf_FindMeshMotionGCLLinearization(All, OutputName, -1, pJ_GCL);
    if (ierr == xf_OK){
      if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
        ierr = xf_Error(xf_SetZeroVector(*pJ_GCL));
        if (ierr != xf_OK) return ierr;
      }
    }
    else (*pJ_GCL) = NULL;
  }
  else (*pJ_GCL) = NULL;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateCutPlaneIntersect
static int
xf_CreateCutPlaneIntersect(xf_CutPlaneIntersect **pCutPlaneIntersect)
{
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pCutPlaneIntersect, 1, 
                           sizeof(xf_CutPlaneIntersect)));
  if (ierr != xf_OK) return ierr;
  
  (*pCutPlaneIntersect)->nelem = 0;
  (*pCutPlaneIntersect)->egrp  = NULL;
  (*pCutPlaneIntersect)->elem  = NULL;
  (*pCutPlaneIntersect)->nquad = NULL;
  (*pCutPlaneIntersect)->xquad = NULL;
  (*pCutPlaneIntersect)->wquad = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyCutPlaneIntersect
int
xf_DestroyCutPlaneIntersect(xf_CutPlaneIntersect *CutPlaneIntersect)
{
  
  if (CutPlaneIntersect == NULL) return xf_OK;
  
  xf_Release( (void  *) CutPlaneIntersect->egrp);
  xf_Release( (void  *) CutPlaneIntersect->elem);
  xf_Release( (void  *) CutPlaneIntersect->nquad);
  xf_Release2((void **) CutPlaneIntersect->xquad);
  xf_Release2((void **) CutPlaneIntersect->wquad);
  
  xf_Release( (void  *) CutPlaneIntersect);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_OutputExactErrorNorm
int
xf_OutputExactErrorNorm(xf_All *All, xf_Vector *U, xf_Vector *Ue, real *Value)
{
  int ierr, k, sr, sr2, iq, nq, pnq, dim;
  int Order, Ordere, QuadOrder;
  int egrp, elem;
  enum xfe_BasisType Basis, Basise;
  enum xfe_Bool found, QuadChanged;
  real *xq, *wq, *u, *ue, temp, Volume;
  real *EU, *EUe;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData, *PhiDatae;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  sr   = All->EqnSet->StateRank;
  
  // set Value to 0
  (*Value) = 0.0;
  
  
  // initialize vars to NULL
  QuadData    = NULL;
  PhiData     = NULL;
  PhiDatae    = NULL;
  JData       = NULL;
  u           = NULL;
  wq          = NULL;
  ue          = NULL;
  Volume      = 0.;
  pnq         = -1;  // previous number of quad points
  
  // consistency check
  if ((U->nArraySelf  != Mesh->nElemGroup) ||
      (Ue->nArraySelf != Mesh->nElemGroup)) return xf_Error(xf_INPUT_ERROR);
  
  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    // consistency check
    if (U->GenArray[egrp].n != Ue->GenArray[egrp].n) return xf_Error(xf_INPUT_ERROR);
    
    // Determine Basis and Order for U and Ue
    Basis  =  U->Basis[egrp];
    Basise = Ue->Basis[egrp];
    
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // Determine interpolation order for U and Ue
      Order  = xf_InterpOrder(U , egrp, elem);
      Ordere = xf_InterpOrder(Ue, egrp, elem);
      
      // determine required integration order
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*max(Order, Ordere), &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      
      /* Pull off quad points for the element; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      
      // compute basis functions
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_EvalBasis(Basise, Ordere, QuadChanged, nq, xq, xfb_Phi, &PhiDatae));
      if (ierr != xf_OK) return ierr;
      
      
      /* Compute geometry Jacobian; if not constant, compute at quad
       points.  Note if jacobian is constant, only one Jacobian will
       be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &ue, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
        if (ierr != xf_OK) return ierr;
      }
      
      EU  =  U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
      EUe = Ue->GenArray[egrp].rValue[elem]; // Ue on elem [nn*sr]
      
      // interpolate state at quad points
      xf_MxM_Set(PhiData->Phi , EU , nq, PhiData->nn , sr, u); 
      xf_MxM_Set(PhiDatae->Phi, EUe, nq, PhiDatae->nn, sr, ue);     
      
      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
        wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      // sum (u-ue)*wq over quad points, add to Value
      for (iq=0; iq<nq; iq++){
        for (k=0; k<sr; k++){
          temp = u[iq*sr+k] - ue[iq*sr+k];
          //xf_printf("[%d,%d,%d] %.15E\n", egrp, elem, iq, temp);
          (*Value) += wq[iq]*temp*temp;
        }
        Volume += wq[iq];
      }
      
      pnq = nq;
      
    } // elem
    
  } // egrp
  
  // divide out the Volume
  (*Value) /= Volume;
  
  // take sqrt of Value
  (*Value) = sqrt((*Value));
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(PhiDatae, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) u);
  xf_Release( (void *) wq);
  xf_Release( (void *) ue);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_OutputDomainIntegral
static int 
xf_OutputDomainIntegral( xf_All *All, xf_Output *Output, const xf_Vector *U,
                        real weight, real *Value, xf_Vector *Value_U, 
                        enum xfe_AddType AddFlag)
{
  int ierr, i, k, sr, sr2, iq, nq, pnq, dim;
  int Order, QuadOrder, QuadOrder0, IntOrder, nn;
  int egrp, elem, nelemtot;
  int *IParam;
  enum xfe_BasisType Basis;
  enum xfe_Bool found, QuadChanged, Need_WD, DoOffset = xfe_False;
  enum xfe_Bool MotionOn = xfe_False;
  real *RParam, *xq, *wq, *u, *gu, *s, *s_u;
  real *EU, *EV_U, *EV_G, LocValue, Time;
  real *xglob, *f, *wd, wdf, offset;
  real DiscreteL2Sum, vol, Volume, ElemLocValue, dp;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData, *GeomPhiData, *WDPhiData;
  xf_Vector *Vtemp, *WD;
  xf_Vector *J_GCL = NULL;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  xf_MotionData *MD = NULL;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  
  EqnSet = All->EqnSet;
  sr     = EqnSet->StateRank;
  sr2    = sr*sr;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // set LocValue to 0
  LocValue = 0.0;
  
  // used for discrete L2 error norms
  DiscreteL2Sum = 0;
  
  // volume for PerVol L2 norm
  Volume = 0.;
  
  // need a temporary vector if doing an L2 norm
  if ( ((Output->DomainNorm == xfe_DomainNormL2) ||
        (Output->DomainNorm == xfe_DomainNormL2PerVol))
      && (Value_U != NULL)){
    ierr = xf_Error(xf_FindSimilarVector(All, Value_U, "Vtemp", xfe_False, 
                                         xfe_True, NULL, &Vtemp, NULL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_SetZeroVector(Vtemp));
    if (ierr != xf_OK) return ierr;
  }
  
  // linearization should not be requested on an error norm output ...
  // not impossible to implement, but cannot think of when would need it
  if ((Output->DomainNorm == xfe_DomainNormL2Error) && (Value_U != NULL))
    return xf_Error(xf_NOT_SUPPORTED);
  
  // linearization not implemented for a discrete L2 norm
  if ((Output->DomainNorm == xfe_DomainNormL2Discrete) && (Value_U != NULL))
    return xf_Error(xf_NOT_SUPPORTED);
  
  // check if need a wall distance or offset
  Need_WD  = xfe_False;
  DoOffset = xfe_False;
  if ((Output->Function != NULL) && (strncmp(Output->Function, "WallDistance", 12) == 0)){
    /*
     The wall distance is used to introduce a multiplicative factor
     into an interior domain integral.  This factor is exp(-wdf*wd^2),
     where wdf is the wall distance decay factor, prescribed in
     Output->Data.
     */
    if ((Output->Data != NULL) && (sscanf(Output->Data, "%lf", &wdf) == 1)){ 
      ierr = xf_Error(xf_FindSupportedVector(All, "WallDistance", &WD));
      if (ierr != xf_OK) return ierr;
      Need_WD = xfe_True;
    }
    else
      xf_printf("Warning: wall distance requested but decay factor not specified in Data.\n");
  }
  else{
    // a value in Output->Data indicates a desired offset in scalar value
    if ((Output->Data != NULL) && (sscanf(Output->Data, "%lf", &offset) == 1)) 
      DoOffset = xfe_True;
  }
  
  // determine if we need mesh motion
  MotionOn = ((Mesh->Motion != NULL) && (Mesh->Motion->Active));
  if (MotionOn){
    ierr = xf_Error(xf_CreateMotionData(All, &MD));
    if (ierr != xf_OK) return ierr;
  }
  
  // locate J_GCL (zero out if AddFlag is set/neg)
  if (Value_U != NULL){ // only if also calculating Value_U
    ierr = xf_Error(xf_FindOutputGCLLinearization(All, Output->Name, AddFlag, &J_GCL));
    if (ierr != xf_OK) return ierr;
  }
  else J_GCL = NULL;
  
  // determine Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
  if (ierr != xf_OK) return ierr;
  
  // initialize vars to NULL
  QuadData    = NULL;
  PhiData     = NULL;
  JData       = NULL;
  u           = NULL;
  gu          = NULL;
  wq          = NULL;
  s           = NULL;
  s_u         = NULL;
  f           = NULL;
  xglob       = NULL;
  GeomPhiData = NULL;
  WDPhiData   = NULL; // for wall distance basis
  wd          = NULL; // for wall distance values
  
  pnq         = -1;  // previous number of quad points
  
  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    // Determine Basis and Order from the state, U
    Basis = U->Basis[egrp];
    
    ElemLocValue = LocValue;
    
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);
      
      // determine required integration order using eqn-set specific rule
      IntOrder = Order;
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, EqnSet, egrp, IntOrder, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      /* need higher order integration if doing L2, since integrand is
       squared but do not account for eqn-set here; could lose
       accuracy if calculating the L2 norm of a highly-nonlinear
       function of the state.  If we did account for eqn-set, we
       would run into max quad order problems in some cases. */
      if ((Output->DomainNorm == xfe_DomainNormL2) ||
          (Output->DomainNorm == xfe_DomainNormL2PerVol) ||
          (Output->DomainNorm == xfe_DomainNormL2Error)){
        ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*IntOrder, &QuadOrder0));
        if (ierr != xf_OK) return ierr;
        QuadOrder = max(QuadOrder0, QuadOrder);
      }
      
      /* Pull off quad points for the element; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions (and grads) if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi | xfb_GPhi | xfb_gPhi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      /* Compute geometry Jacobian; if not constant, compute at quad
       points.  Note if jacobian is constant, only one Jacobian will
       be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;
      
      nn = PhiData->nn; // number of interpolation nodes
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **)  &gu, dim*nq*sr, sizeof(real))); 
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **)  &s, nq, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        if (Value_U != NULL){
          ierr = xf_Error(xf_ReAlloc( (void **)  &s_u, nq*sr, sizeof(real)));
          if (ierr != xf_OK) return ierr;
        }
        if (Output->DomainNorm == xfe_DomainNormL2Error){
          ierr = xf_Error(xf_ReAlloc( (void **)  &f, nq*sr, sizeof(real)));
          if (ierr != xf_OK) return ierr;
        }
      }
      
      // obtain global coords of quad points
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
                                      nq, xq, xglob));
      if (ierr != xf_OK) return ierr;
      
      // obtain transformation map if doing mesh motion
      if (MotionOn){
        ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, Mesh->Motion, 
                                         nq, dim, Time, xglob, MD));
        if (ierr != xf_OK) return ierr;
      }
      
      EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
      
      // interpolate state at quad points
      xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u);     
      
      // transform state to physical
      if (MotionOn) xf_ModMotionPreEqnCall(nq, dim, sr, MD, u, NULL);
      
      /* Compute geometry Jacobian */
      ierr = xf_Error(xf_ElemJacobian(All->Mesh, egrp, elem, nq, xq, 
                                      xfb_detJ | xfb_iJ | xfb_J, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;
      
      /* convert reference basis grads (GPhi) to physical grads, gPhi */
      ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData)); 
      if (ierr != xf_OK) return ierr;
      
      // compute state gradient
      for (i=0; i<dim; i++)
        xf_MxM_Set(PhiData->gPhi+nn*nq*i, EU, nq, nn, sr, gu + nq*sr*i);
      
      // transform state gradient to physical
      if (MotionOn) xf_ModMotionPhysGrad(nq, dim, sr, MD, u, gu);
      
      // form detJ-multiplied quad weight vector, wq, multiplied by weight
      for (iq=0; iq<nq; iq++) 
        wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)]*weight;
      
      // account for mesh motion in quad weights (elements are deformed in physical space)
      if (MotionOn) xf_ColMult(wq, MD->g, nq, 1, 1);
      
      // interpolate wall distance if necessary and adjust quad weights
      if (Need_WD){
        if (nq > pnq){
          ierr = xf_Error(xf_ReAlloc( (void **) &wd, nq*WD->StateRank, sizeof(real)));
          if (ierr != xf_OK) return ierr;
        }
        ierr = xf_Error(xf_EvalBasis(WD->Basis[egrp], WD->Order[egrp], QuadChanged, 
                                     nq, xq, xfb_Phi, &WDPhiData));
        if (ierr != xf_OK) return ierr;
        xf_MxM_Set(WDPhiData->Phi, WD->GenArray[egrp].rValue[elem], nq, 
                   WDPhiData->nn, WD->StateRank, wd);
        
        // quad weights are adjusted here to incorporate the wall distance
        for (iq=0; iq<nq; iq++) wq[iq] *= exp(-wdf*wd[iq]*wd[iq]);
      }
      
      
      /** Distinguish between state error and scalar norms **/
      
      if (Output->DomainNorm == xfe_DomainNormL2Error){
        
        // call eqnset specific function for state values
        ierr = xf_Error(xf_EqnSetFcnState(EqnSet, Output->Function, Output->Data, IParam, 
                                          RParam, nq, (MotionOn) ? MD->x : xglob, &Time, f));
        if (ierr != xf_OK) return ierr;
        
        // set s = square error at each quad point
        for (iq=0; iq<nq; iq++){
          for (k=0, s[iq]=0.; k<sr; k++){
            i = iq*sr+k;
            s[iq] += (u[i]-f[i])*(u[i]-f[i]);
          }
        }
        
        // linearization is possible, but currently see no use for it
        if (Value_U != NULL) return xf_Error(xf_NOT_SUPPORTED);
        
      }
      else{
        
        // call eqnset specific function for scalar
        ierr = xf_Error(xf_EqnSetScalar(EqnSet, Output->ScalarName, IParam, 
                                        RParam, nq, u, gu, s, s_u, NULL, NULL, 0.0));
        if (ierr != xf_OK) return ierr;
        
        // Divide output linearization by gbar to make it wrt reference state
        if ((MotionOn) && (Value_U != NULL) && (s_u != NULL))
          xf_ColDiv(s_u, MD->gb, nq, sr, 1);
        
        // subtract an offset if specified
        if (DoOffset) for (iq=0; iq<nq; iq++) s[iq] -= offset;
        
        // subtract linearization from adjoint residual if not null
        if (Value_U != NULL){
          if ((Output->DomainNorm == xfe_DomainNormL2) || (Output->DomainNorm == xfe_DomainNormL2PerVol))
            EV_U = Vtemp->GenArray[egrp].rValue[elem];
          else
            EV_U = Value_U->GenArray[egrp].rValue[elem]; // Value_U on elem [nn*sr]
          
          // multiply s_u by wq
          xf_ColMult(s_u, wq, nq, sr, 1);
          
          /* 
           Need to set, add-to, subtract-from, etc. Value_U:
           
           No norm: Value   = int(s   dx)
           Value_U = int(s_U dx)
           
           L2 Norm: Value   = sqrt(int(s^2 dx))
           Value_U = int(s*s_U dx) * 1/Value
           
           On each elem, s_U{n,k} = s_u{q,k} * Phi{q,n}
           */
          
          if ((Output->DomainNorm == xfe_DomainNormL2) || (Output->DomainNorm == xfe_DomainNormL2PerVol))
            xf_ColMult(s_u, s, nq, sr, 1); // multiply s_u by s
          else if (Output->DomainNorm != xfe_DomainNormNone)
            return xf_Error(xf_NOT_SUPPORTED);
          
          // Value_U{n,k} @= int(s_u{q,k} * Phi{q,n}) 
          xf_MTxM(PhiData->Phi, s_u, nn, nq, sr, AddFlag, EV_U);
          
          // Value_G{n} @= int(s_u{q,k} * (-u{q,k}) * Phi{q,n}) -- u is physical here
          if (J_GCL != NULL){
            EV_G = J_GCL->GenArray[egrp].rValue[elem]; // J_GCL on elem [nn*sr]
            for (iq=0; iq<nq; iq++){
              xf_DotProduct(s_u+iq*sr, u+iq*sr, sr, &dp);
              s_u[iq] = -dp; // s_u is overwritten here
            }
            xf_MTxM(PhiData->Phi, s_u, nn, nq, 1, AddFlag, EV_G);
          }
          
          
        }
        
        // modify scalar depending on norm
        if (Output->DomainNorm == xfe_DomainNormL1)
          for (iq=0; iq<nq; iq++) s[iq] = fabs(s[iq]);
        else if ((Output->DomainNorm == xfe_DomainNormL2) || (Output->DomainNorm == xfe_DomainNormL2PerVol))
          for (iq=0; iq<nq; iq++) s[iq] *= s[iq];
        // else do nothing (L2Discrete norm just uses s[iq])
        
      }
      
      // sum scalar*wq over quad points, add to Value
      for (iq=0; iq<nq; iq++) LocValue += wq[iq]*s[iq];
      
      // sum quad weights to get volume
      for (iq=0, vol=0; iq<nq; iq++) vol += wq[iq];
      Volume += vol; // increment running total	
      
      // discrete L2 norm requires element averages
      if (Output->DomainNorm == xfe_DomainNormL2Discrete){
        ElemLocValue = LocValue - ElemLocValue;
        ElemLocValue /= vol;
        DiscreteL2Sum += ElemLocValue*ElemLocValue; // sum of square averages
      }
      
      pnq = nq;
    } // elem
    
  } // egrp
  
  if (Output->DomainNorm == xfe_DomainNormL2Discrete){
    // special sum/reduce for discrete L2 norm
    ierr = xf_Error(xf_MPI_Allreduce(&DiscreteL2Sum, 1, xfe_SizeReal, xfe_MPI_SUM));
    if (ierr != xf_OK) return ierr;
    
    nelemtot = Mesh->ElemGroup[egrp].nElem;
    ierr = xf_Error(xf_MPI_Allreduce(&nelemtot, 1, xfe_SizeInt, xfe_MPI_SUM));
    if (ierr != xf_OK) return ierr;
    
    LocValue = DiscreteL2Sum / ( (real) nelemtot);
  }
  else{
    // sum reduce LocValue
    ierr = xf_Error(xf_MPI_Allreduce(&LocValue, 1, xfe_SizeReal, xfe_MPI_SUM));
    if (ierr != xf_OK) return ierr;
  }
  
  // take sqrt of Value if L2 norm
  if ((Output->DomainNorm == xfe_DomainNormL2) ||
      (Output->DomainNorm == xfe_DomainNormL2Error) ||
      (Output->DomainNorm == xfe_DomainNormL2PerVol) ||
      (Output->DomainNorm == xfe_DomainNormL2Discrete)){
    
    if (Volume == 0.) return xf_Error(xf_OUT_OF_BOUNDS);
    if (Output->DomainNorm == xfe_DomainNormL2PerVol) weight *= 1./Volume; // per volume comes in here
    
    LocValue = sqrt(fabs(LocValue*weight));
    if (weight < 0) LocValue = -LocValue;
    
    // need to multiply linearization by 1/Value (see above description)
    if ((weight != 0.) && (Value_U != NULL)){
      ierr = xf_Error(xf_VectorMultSet(Vtemp, fabs(weight/LocValue), AddFlag, Value_U));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  if (Value != NULL) (*Value) = LocValue;
  
  
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Wall Distance Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(WDPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy mesh motion data */
  xf_DestroyMotionData(MD);
  
  xf_Release( (void  *) IParam);
  xf_Release( (void  *) RParam);
  xf_Release( (void *) u);
  xf_Release( (void *) gu);
  xf_Release( (void *) wq);
  xf_Release( (void *) s);
  xf_Release( (void *) s_u);
  xf_Release( (void *) f);
  xf_Release( (void *) xglob);
  xf_Release( (void *) wd);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_OutputBoundaryIntegral_VectorScalar
static int 
xf_OutputBoundaryIntegral_VectorScalar( xf_All *All, xf_Output *Output, enum xfe_Bool VectorFlag,
                                       const xf_Vector *U, real weight, real *Value, 
                                       xf_Vector *Value_U, enum xfe_AddType AddFlag)
{
  int ierr, k, i, d;
  int sr, sr2, dim, nq, iq;
  int nset, nBFG, *BFGs;
  int egrp, elem, face;
  int ibfg, ibfgrp, ibface;
  int Order, QuadOrder, pnq, nn;
  int *IParam;
  enum xfe_Bool QuadChanged;
  enum xfe_Bool MotionOn;
  enum xfe_AddType AddFlag2;
  enum xfe_BasisType Basis;
  real LocValue, dp, motionfac;
  real *wn, *xglob, *xelem;
  real *uI, *uB, *uB_uI;
  real *vn_uI, *v_uB, *v;
  real *xq, *wq;
  real *EU, *EV_U, *EV_G, *RParam;
  real NN, Time;
  xf_Vector *J_GCL = NULL;
  xf_BC *BC;
  xf_BFace BFace;
  xf_QuadData *QuadData;
  xf_BasisTable *PhiTable;
  xf_BasisData *PhiData, *GeomPhiData;
  xf_NormalData *NData;
  xf_MotionData *MD = NULL;
  xf_EqnSet *EqnSet;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  EqnSet = All->EqnSet;
  dim = Mesh->Dim;
  sr = EqnSet->StateRank;
  sr2 = sr*sr;
  
  LocValue = 0.0;
  
  AddFlag2 = xf_GetAddFlag2(AddFlag); // used for additional opers; add orsub
  
  if (Value_U != NULL){
    // zero out Value_U if requesting a set or neg
    if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
      ierr = xf_Error(xf_SetZeroVector(Value_U));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, 
                                       NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // Sort EqnSet->BCs to match boundary face groups
  ierr = xf_Error(xf_SortEqnSetBCs(Mesh, EqnSet->BCs+0));
  if (ierr != xf_OK) return ierr;
  BC = EqnSet->BCs[0].BC;
  
  // determine which boundary face groups to integrate over
  nBFG = Output->nBFG;
  ierr = xf_Error(xf_Alloc((void **) &BFGs, nBFG, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  nset = 0;
  for (i=0; i<nBFG; i++){
    for (k=0; k<All->Mesh->nBFaceGroup; k++)
      if (strcmp(Output->BFGTitles[i], All->Mesh->BFaceGroup[k].Title) == 0){
        BFGs[i] = k;
        nset++;
      }
  }
  if ((nset != nBFG) || (nBFG > All->Mesh->nBFaceGroup))
    return xf_Error(xf_OUT_OF_BOUNDS);
  
  // determine Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
  if (ierr != xf_OK) return ierr;
  
  // Are we doing mesh motion?
  MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));
  if (MotionOn){
    ierr = xf_Error(xf_CreateMotionData(All, &MD));
    if (ierr != xf_OK) return ierr;
  }
  
  // locate J_GCL (zero out if AddFlag is set/neg)
  if (Value_U != NULL){ // only if also calculating Value_U
    ierr = xf_Error(xf_FindOutputGCLLinearization(All, Output->Name, AddFlag, &J_GCL));
    if (ierr != xf_OK) return ierr;
  }
  else J_GCL = NULL;
  
  
  
  /* Create a basis table, PhiTable, that will store computed basis
   functions specific to each [element shape, face in element,
   orientation of face] combination, for quick lookup. */
  ierr = xf_Error(xf_CreateBasisTable(&PhiTable));
  if (ierr != xf_OK) return ierr;
  
  // set variables pre-loop
  QuadData    =  NULL;
  PhiData     =  NULL; 
  NData       =  NULL;
  GeomPhiData =  NULL;
  wn	      =  NULL;
  xglob       =  NULL;
  xelem       =  NULL; 
  uI          =  NULL; 
  uB          =  NULL;
  uB_uI       =  NULL;
  vn_uI       =  NULL;
  v_uB        =  NULL;
  v           =  NULL;
  pnq         =  -1;
  
  // loop over desired boundary groups
  for (ibfg=0; ibfg<nBFG; ibfg++){
    ibfgrp = BFGs[ibfg];
    
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      
      BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
      egrp = BFace.ElemGroup;
      elem = BFace.Elem;
      face = BFace.Face;
      
      Basis = U->Basis[egrp];
      Order = xf_InterpOrder(U, egrp, elem);
      
      // determine required integration order
      ierr = xf_Error(xf_GetQuadOrderBFace(Mesh, EqnSet, BFace, Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      /* Pull off quad points for the bface; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face, QuadOrder, 
                                  &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      wq = QuadData->wquad;
      
      // compute basis functions if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrp, elem, face, BFace.Orient,
                                                   Basis, Order, QuadChanged, nq, xq, 
                                                   xfb_Phi, &PhiData, PhiTable, &xelem));
      if (ierr != xf_OK) return ierr;
      
      /* Compute normal(s) at quad points.  If face is straight, only
       one normal will be computed/returned. */
      ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, nq, xq, &NData, NULL));
      if (ierr != xf_OK) return ierr;
      
      nn = PhiData->nn;
      
      EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **) &wn, nq*dim, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &uI, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &uB, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **)  &v, nq*dim, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        if (Value_U != NULL){
          ierr = xf_Error(xf_ReAlloc( (void **) &uB_uI, nq*sr2, sizeof(real)));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_ReAlloc( (void **)  &v_uB, nq*dim*sr, sizeof(real)));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_ReAlloc( (void **)  &vn_uI, nq*dim*sr, sizeof(real)));
          if (ierr != xf_OK) return ierr;
        }
        
      }
      
      // obtain global coords of quad points
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, xfe_True, 
                                      nq, xelem, xglob));
      if (ierr != xf_OK) return ierr;      
      
      // interpolate state at quad points
      xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, uI);
      
      // obtain transformation map if doing mesh motion
      if (MotionOn){
        ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, All->Mesh->Motion, 
                                         nq, dim, Time, xglob, MD));
        if (ierr != xf_OK) return ierr;
      }
      
      /* Initially, set wn = just the normals.  Will include quad
       weights after call to EqnSetBCState, since some weights may
       be negative while EqnSetBCState needs to have the true
       outward-pointing normal at each point.
       */
      for (d=0; d<dim; d++)
        for (iq=0;iq<nq; iq++) 
          wn[iq*dim+d] = NData->n[iq*dim*(NData->nq!=1)+d];
      
      // transform state + normal to physical
      if (MotionOn) xf_ModMotionPreEqnCall(nq, dim, sr, MD, uI, wn);
      
      // obtain boundary state and derivative at quad points
      ierr = xf_Error(xf_EqnSetBCState(EqnSet, BC+ibfgrp, IParam, RParam,
                                       nq, wn, (MotionOn) ? MD->x : xglob, 
                                       &Time, (MotionOn) ? MD->vg : NULL,
                                       uI, uB, uB_uI));
      if (ierr != xf_OK) return ierr;
      
      if (VectorFlag){
        
        // multiply wn by quad weights for further use in integration
        for (iq=0;iq<nq; iq++) 
          for (d=0; d<dim; d++)
            wn[iq*dim+d] *= wq[iq]*weight; // weight included here
        
        // call eqnset specific function for vector
        ierr = xf_Error(xf_EqnSetVector(EqnSet, Output->VectorName, IParam, 
                                        RParam, nq, uB, v, v_uB));
        if (ierr != xf_OK) return ierr;
        
        // dot v with normal, add to value
        xf_DotProduct(wn, v, nq*dim, &dp);
        LocValue += dp;
        
        // add linearization
        if (Value_U != NULL){
          
          EV_U = Value_U->GenArray[egrp].rValue[elem]; // Value_U on elem [nn*sr]
          
          // multiply v_uB by wn
          for (iq=0; iq<nq; iq++){
            //if mesh motion active, divide linearization by gbar to make wrt reference state
            motionfac = (MotionOn) ? 1./(MD->gb[iq]) : 1.; 
            for (d=0; d<dim; d++)
              for (k=0; k<sr; k++) v_uB[iq*dim*sr+d*sr+k] *= wn[iq*dim+d]*motionfac;
          }
          // form vn_uI = uB_uI^T*v_uB
          for (k=0; k<nq*sr; k++) vn_uI[k] = 0;
          for (iq=0; iq<nq; iq++)
            for (d=0; d<dim; d++)
              xf_MTxM_Add(uB_uI+iq*sr*sr, v_uB+iq*dim*sr+d*sr, sr, sr, 1, vn_uI+iq*sr);
          
          // EV_U{n,k} @= Phi{q,n}^T * vn_uI{q,k} 
          xf_MTxM(PhiData->Phi, vn_uI, nn, nq, sr, AddFlag2, EV_U);
          
          // EV_G{n} @= Phi{q,n}^T * vn_uI{q,k}*(-uI{q,k})  -- note uI is physical here
          if (J_GCL != NULL){
            EV_G = J_GCL->GenArray[egrp].rValue[elem]; // J_GCL on elem [nn*sr]
            for (iq=0; iq<nq; iq++){
              xf_DotProduct(vn_uI+iq*sr, uI+iq*sr, sr, v_uB+iq); // v_uB is temp storage here
              v_uB[iq] /= -1.0;
            }
            xf_MTxM(PhiData->Phi, v_uB, nn, nq, 1, AddFlag2, EV_G);
          }
        }
      }
      else{
        
        // multiply wn by quad weights for further use in integration
        for (iq=0;iq<nq; iq++) {
          for (NN=0., d=0; d<dim; d++)
            NN += wn[iq*dim+d]*wn[iq*dim+d]; // NN is face area Jacobian
          wn[iq] = sqrt(NN)*wq[iq]*weight; // weight included here (overwriting wn)
        }
        
        // call eqnset specific function for scalar
        ierr = xf_Error(xf_EqnSetScalar(EqnSet, Output->ScalarName, IParam, 
                                        RParam, nq, uB, NULL, v, v_uB, NULL, NULL, 0.0));
        if (ierr != xf_OK) return ierr;
        
        //Divide output linearization by gbar to make it wrt reference state
        if ((MotionOn) && (Value_U != NULL) && (v_uB != NULL)){
          xf_ColDiv(v_uB, MD->gb, nq, sr, 1);
        }
        
        // dot v with quad weights, add to value
        xf_DotProduct(wn, v, nq, &dp);
        LocValue += dp;
        
        // add linearization
        if (Value_U != NULL){
          
          EV_U = Value_U->GenArray[egrp].rValue[elem]; // Value_U on elem [nn*sr]
          
          // multiply v_uB by wn
          for (iq=0; iq<nq; iq++)
            for (k=0; k<sr; k++) v_uB[iq*sr+k] *= wn[iq];
          
          // form vn_uI = uB_uI^T*v_uB
          for (k=0; k<nq*sr; k++) vn_uI[k] = 0;
          for (iq=0; iq<nq; iq++)
            xf_MTxM_Add(uB_uI+iq*sr*sr, v_uB+iq*sr, sr, sr, 1, vn_uI+iq*sr);
          
          // EV_U{n,k} @= vn_uI{q,k} Phi{q,n}) 
          xf_MTxM(PhiData->Phi, vn_uI, nn, nq, sr, AddFlag2, EV_U);
          
          // EV_G{n} @= Phi{q,n}^T * vn_uI{q,k}*(-uI{q,k})  -- note uI is physical here
          if (J_GCL != NULL){
            EV_G = J_GCL->GenArray[egrp].rValue[elem]; // J_GCL on elem [nn*sr]
            for (iq=0; iq<nq; iq++){
              xf_DotProduct(vn_uI+iq*sr, uI+iq*sr, sr, v_uB+iq); // v_uB is temp storage here
              v_uB[iq] *= -1.0;
            }
            xf_MTxM(PhiData->Phi, v_uB, nn, nq, 1, AddFlag2, EV_G);
          }
          
        } // end if Value_U != NULL
        
      }
      
    } // ibface
  } // ibfg
  
  
  if (Value != NULL) (*Value) = LocValue;
  
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Table */
  ierr = xf_Error(xf_DestroyBasisTable(PhiTable));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Normal Data */
  ierr = xf_Error(xf_DestroyNormalData(NData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy mesh motion data */
  xf_DestroyMotionData(MD);
  
  xf_Release( (void *) BFGs);
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) wn);
  xf_Release( (void *) xglob);
  xf_Release( (void *) xelem);
  xf_Release( (void *) uI);
  xf_Release( (void *) uB);
  xf_Release( (void *) uB_uI);
  xf_Release( (void *) vn_uI);
  xf_Release( (void *) v_uB);
  xf_Release( (void *) v);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_OutputBoundaryIntegral
static int 
xf_OutputBoundaryIntegral( xf_All *All, xf_Output *Output, const xf_Vector *U, 
                          real weight, real *Value, xf_Vector *Value_U, 
                          enum xfe_AddType AddFlag)
{
  int ierr, k, i, sr, nset, nBFG, *BFGs;
  int myRank, nProc;
  enum xfe_Bool NegativeFlag = xfe_False;
  char DumpFile[xf_MAXSTRLEN];
  int  *FluxMoments;
  real *FluxWeights;
  xf_Vector *J_GCL = NULL;
  xf_EqnSet *EqnSet;
  xf_OutputEvalData OutputEval; 
  
  EqnSet = All->EqnSet;
  sr = EqnSet->StateRank;
  
  // If not using a flux, can use a vector dot n, or a scalar
  if (!Output->UsesFlux){
    return xf_Error(xf_OutputBoundaryIntegral_VectorScalar(All, Output, Output->VectorName != NULL, 
                                                           U, weight, Value, Value_U, AddFlag));
  }
  
  if (Output->nFluxComponent <= 0) return xf_Error(xf_INPUT_ERROR);
  if (Output->FluxComponentWeights == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // allocate a flux weight vector
  ierr = xf_Error(xf_Alloc((void **) &FluxWeights, sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &FluxMoments, sr, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (k=0; k<sr; k++) FluxWeights[k] =  0.;
  for (k=0; k<sr; k++) FluxMoments[k] = -1;
  
  // set weights according to request; include output weight here
  nset = 0;
  for (i=0; i<Output->nFluxComponent; i++){
    for (k=0; k<sr; k++)
      if (strcmp(Output->FluxComponentNames[i], EqnSet->StateName[k]) == 0){
        FluxWeights[k] = Output->FluxComponentWeights[i]*weight;
        if (Output->FluxComponentMoments != NULL)
          FluxMoments[k] = Output->FluxComponentMoments[i];
        nset++;
      }
  }
  if (nset != Output->nFluxComponent) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // determine which boundary face groups to integrate over
  nBFG = Output->nBFG;
  ierr = xf_Error(xf_Alloc((void **) &BFGs, nBFG, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  nset = 0;
  for (i=0; i<nBFG; i++){
    for (k=0; k<All->Mesh->nBFaceGroup; k++)
      if (strcmp(Output->BFGTitles[i], All->Mesh->BFaceGroup[k].Title) == 0){
        BFGs[i] = k;
        nset++;
      }
  }
  if ((nset != nBFG) || (nBFG > All->Mesh->nBFaceGroup))
    return xf_Error(xf_OUT_OF_BOUNDS);
  
  if (Value_U != NULL){
    // zero out Value_U if requesting a set or neg
    if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
      ierr = xf_Error(xf_SetZeroVector(Value_U));
      if (ierr != xf_OK) return ierr;
    }
    // we will use a function that adds to Value_U; so set weights if want neg
    if ((AddFlag == xfe_Sub) || (AddFlag == xfe_Neg)){
      NegativeFlag = xfe_True;
      for (k=0; k<sr; k++) FluxWeights[k] *= -1.;
    }
  }
  
  // locate J_GCL (zero out if AddFlag is set/neg)
  if (Value_U != NULL){ // only if also calculating Value_U
    ierr = xf_Error(xf_FindOutputGCLLinearization(All, Output->Name, AddFlag, &J_GCL));
    if (ierr != xf_OK) return ierr;
  }
  else J_GCL = NULL;
  
  // roll output evaluation inputs into a structure
  OutputEval.Value       = Value;
  OutputEval.Value_U     = Value_U;
  OutputEval.Value_G     = J_GCL;
  OutputEval.FluxWeights = FluxWeights;
  OutputEval.FluxMoments = FluxMoments;
  OutputEval.fidDump     = NULL;
  // DumpFile
  if (Output->DumpFile != NULL){
    // check if parallel
    ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
    if (ierr != xf_OK) return ierr;
    
    // use proc number to write to different files in parallel
    strcpy(DumpFile, Output->DumpFile);
    if (nProc > 1) sprintf(DumpFile, "%s.%d\0", Output->DumpFile, myRank);
    
    OutputEval.fidDump = fopen(DumpFile, "w");
    if (OutputEval.fidDump == NULL) 
      xf_printf("Warning, could not open the file %s for output.\n", DumpFile);
  }
  
  // call boundary residual routine to evaluate the output
  ierr = xf_Error(xf_CalculateResidualBFaces(All, U, NULL, NULL, &OutputEval, nBFG, BFGs, NULL));
  if (ierr != xf_OK) return ierr;
  
  // close DumpFile
  if (Output->DumpFile != NULL){
    fclose(OutputEval.fidDump);
  }
  
  // correct Value back to positive if was using negative weights
  if ((NegativeFlag) && (Value != NULL))
    (*Value) *= -1.;
  
  // include weight in output
  if (Value != NULL) (*Value) *= weight;
  
  xf_Release( (void *) FluxWeights);
  xf_Release( (void *) FluxMoments);
  xf_Release( (void *) BFGs);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateCutPlaneIntersect
static int 
xf_CalculateCutPlaneIntersect( xf_All *All, const real *CutPlane, const xf_Vector *U,
                              xf_CutPlaneIntersect **pCPI)
{
  // only for linear elements for now
  // 2D and 3D both supported
  int ierr, egrp, elem, k, d, dim;
  int ielem, QuadOrder, nq0, iq;
  int Order, pOrder, ntri, itri, nqtot;
  int s, nnode;
  int *Node, *vtri, n[3];
  enum xfe_Bool hit;
  real val, jac, *xq0, *wq0, *xi;
  real *xnode, *xref, *xglob, *xtemp, *wtemp;
  real vg1[3], vg2[3], vr1[3], vr2[3], v[3];
  xf_BasisData *PhiData;
  xf_CutPlaneIntersect *CPI;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim = Mesh->Dim;
  
  // no support for mesh motion in cut planes
  if ((Mesh->Motion != NULL) && (Mesh->Motion->Active)) return xf_Error(xf_NOT_SUPPORTED);
  
  // create (*pCPI)
  ierr = xf_Error(xf_CreateCutPlaneIntersect(pCPI));
  if (ierr != xf_OK) return ierr;
  CPI = (*pCPI);
  
  // (over)allocate CPI->egrp, CPI->elem
  for (egrp=0, CPI->nelem=0; egrp < Mesh->nElemGroup; egrp++)
    CPI->nelem += Mesh->ElemGroup[egrp].nElem;  
  
  ierr = xf_Error(xf_ReAlloc((void **) &CPI->egrp, CPI->nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_ReAlloc((void **) &CPI->elem, CPI->nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // count number of possibly-intersecting elems (fill CPI->egrp, CPI->elem)
  for (egrp=0, ielem=0; egrp < Mesh->nElemGroup; egrp++){
    
    // skip non Q1 groups
    if (Mesh->ElemGroup[egrp].QOrder != 1) continue;
    
    for (elem=0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      
      Node = Mesh->ElemGroup[egrp].Node[elem];
      
      // test for intersection using element nodes
      hit = xfe_False;
      for (k=0; k<Mesh->ElemGroup[egrp].nNode; k++){
        xnode = Mesh->Coord[Node[k]];
        for (d=0, val=CutPlane[dim]; d<dim; d++) val += CutPlane[d]*xnode[d];
        if (fabs(val) < MEPS){
          hit = xfe_True;
          break;
        }
        if (k == 0) s = sign(val);
        else{
          if (sign(val) != s){
            hit = xfe_True;
            break;
          }
        }
      } // k
      if (hit){
        CPI->egrp[ielem] = egrp;
        CPI->elem[ielem] = elem;
        ielem++;
      }
    } // elem
    
  } // egrp
  
  // reallocate CPI->egrp, CPI->elem
  CPI->nelem = ielem;
  
  ierr = xf_Error(xf_ReAlloc((void **) &CPI->egrp, CPI->nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_ReAlloc((void **) &CPI->elem, CPI->nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // allocate nquad, xquad, wquad
  ierr = xf_Error(xf_Alloc((void **) &CPI->nquad, CPI->nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &CPI->xquad, CPI->nelem, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &CPI->wquad, CPI->nelem, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  
  // determine max required integration order (for allocation)
  QuadOrder = 0;
  for (egrp=0; egrp < Mesh->nElemGroup; egrp++){
    ierr = xf_Error(xf_GetQuadOrderElem(Mesh, All->EqnSet, egrp, U->Order[egrp], &Order));
    if (ierr != xf_OK) return ierr;
    QuadOrder = max(Order, QuadOrder);
  }
  
  // obtain quad points on reference triangle
  xq0 = NULL;
  wq0 = NULL;
  if (dim == 3){
    ierr = xf_Error(xf_QuadTriangle(QuadOrder, &nq0, &xq0, &wq0));
    if (ierr != xf_OK) return ierr;
  }
  else{
    ierr = xf_Error(xf_QuadLine(QuadOrder, &nq0, &xq0, &wq0));
    if (ierr != xf_OK) return ierr;
  }
  
  // overallocate xtemp, wtemp (12 = 6 tets per hex * 2 tris per tet)
  nqtot = ((dim == 3) ? 12*nq0*CPI->nelem : nq0*CPI->nelem);
  ierr = xf_Error(xf_Alloc((void **) &xtemp, nqtot*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &wtemp, nqtot, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // loop over intersected elements
  pOrder = -1;
  PhiData = NULL;
  nqtot = 0;
  xref = NULL;
  vtri = NULL;
  xglob = NULL;
  for (ielem=0; ielem<CPI->nelem; ielem++){
    
    // elem info
    egrp = CPI->egrp[ielem];
    elem = CPI->elem[ielem];
    
    // get interpolation order
    Order = xf_InterpOrder(U, egrp, elem);
    
    // quadrature rule
    if (Order != pOrder){
      pOrder = Order;
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, All->EqnSet, egrp, Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      if (dim == 3){
        ierr = xf_Error(xf_QuadTriangle(QuadOrder, &nq0, &xq0, &wq0));
        if (ierr != xf_OK) return ierr;
      }
      else{
        ierr = xf_Error(xf_QuadLine(QuadOrder, &nq0, &xq0, &wq0));
        if (ierr != xf_OK) return ierr;
      }
    }
    
    // call element intersection function
    ierr = xf_Error(xf_IntersectElemWithPlane(Mesh, egrp, elem, 1, CutPlane, NULL,
                                              xfe_True, &PhiData, &nnode, &xref, 
                                              &ntri, &vtri, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    
    // realloc xglob
    ierr = xf_Error(xf_ReAlloc((void **) &xglob, dim*nnode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // determine global coords corresponding to xref
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &PhiData, xfe_True, 
                                    nnode, xref, xglob));
    if (ierr != xf_OK) return ierr;
    
    // set pointers
    CPI->xquad[ielem] = xtemp + dim*nqtot;
    CPI->wquad[ielem] = wtemp + nqtot;
    CPI->nquad[ielem] = 0;
    
    if (dim == 2){
      for (k=0; k<2; k++) v[k] = xglob[2+k] - xglob[k];
      for (k=0, jac=0.; k<2; k++) jac += v[k]*v[k];
      jac = sqrt(jac);
      for (iq=0; iq<nq0; iq++){
        xi = xq0 + iq;
        for (k=0; k<2; k++)
          CPI->xquad[ielem][2*iq+k] =  xref[k] + (xref[2+k]-xref[k])*xi[0];
        CPI->wquad[ielem][iq] = wq0[iq]*jac;
      } // iq
      CPI->nquad[ielem] += nq0;
    }
    else{
      // loop over number of intersecting triangles
      
      for (itri=0; itri<ntri; itri++){
        
        for (k=0; k<3; k++) n[k] = vtri[3*itri+k];
        
        // calculate jac = twice area of triangle
        for (k=0; k<3; k++) vg1[k] = xglob[3*n[1]+k] - xglob[3*n[0]+k];
        for (k=0; k<3; k++) vg2[k] = xglob[3*n[2]+k] - xglob[3*n[0]+k];
        xf_CrossProduct(vg1, vg2, v);
        for (k=0, jac=0.; k<3; k++) jac += v[k]*v[k];
        jac = sqrt(jac);
        
        //xf_printf("ielem = %d, itri = %d, jac = %.10E\n", ielem, itri, jac);
        
        // vectors in ref space
        for (k=0; k<3; k++) vr1[k] = xref[3*n[1]+k] - xref[3*n[0]+k];
        for (k=0; k<3; k++) vr2[k] = xref[3*n[2]+k] - xref[3*n[0]+k];
        
        for (iq=0; iq<nq0; iq++){
          xi = xq0 + 2*iq;
          
          // map ref points into global space; store in xquad
          for (k=0; k<3; k++)
            CPI->xquad[ielem][3*(itri*nq0+iq)+k] = 
            xref[3*n[0]+k] + vr1[k]*xi[0] + vr2[k]*xi[1];
          
          // scale weights by triangle area; store in wquad
          CPI->wquad[ielem][itri*nq0+iq] = wq0[iq]*jac;
        }
        CPI->nquad[ielem] += nq0;
      }
    }
    nqtot += CPI->nquad[ielem];
    
  } // ielem
  
  // reallocate xtemp, wtemp 
  ierr = xf_Error(xf_ReAlloc((void **) &xtemp, dim*nqtot, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc((void **) &wtemp, nqtot, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) xq0);
  xf_Release( (void *) wq0);
  xf_Release( (void *) xref);
  xf_Release( (void *) xglob);
  xf_Release( (void *) vtri);
  
  
  return xf_OK;
}



/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_CalculateLineIntersect */
/* static int  */
/* xf_CalculateLineIntersect( xf_All *All, const real *LineCoord, const xf_Vector *U, */
/* 			   xf_CutPlaneIntersect **pCPI) */
/* { */
/*   // 2D and 3D both supported, for linear and curved elements */
/*   int ierr, egrp, elem, k, d, dim; */
/*   int ielem, QuadOrder, nq0, iq; */
/*   int Order, pOrder, ntri, itri, nqtot; */
/*   int s, nnode; */
/*   int *Node, *vtri, n[3]; */
/*   enum xfe_Bool hit; */
/*   real val, jac, *xq0, *wq0, *xi; */
/*   real *xnode, *xref, *xglob, *xtemp, *wtemp; */
/*   real vg1[3], vg2[3], vr1[3], vr2[3], v[3]; */
/*   xf_BasisData *PhiData; */
/*   xf_CutPlaneIntersect *CPI; */
/*   xf_Mesh *Mesh; */

/*   Mesh = All->Mesh; */
/*   dim = Mesh->Dim; */

/*   // create (*pCPI) */
/*   ierr = xf_Error(xf_CreateCutPlaneIntersect(pCPI)); */
/*   if (ierr != xf_OK) return ierr; */
/*   CPI = (*pCPI); */

/*   // Initialize to NULL */
/*   CPI->nelem = 0; */
/*   CPI->egrp  = NULL; */
/*   CPI->elem  = NULL; */
/*   CPI->nquad = NULL; */
/*   CPI->xquad = NULL; */
/*   CPI->wquad = NULL; */
/*   pOrder     = -1; */
/*   xq0        = NULL; */
/*   wq0        = NULL; */
/*   PhiData1   = NULL; */
/*   PhiData2   = NULL; */
/*   xqglob     = NULL; */


/*   // Loop over element groups */
/*   for (egrp=0, ielem=0; egrp < Mesh->nElemGroup; egrp++){ */
/*     // quadrature rule and points on a [-1,1] line */
/*     if (U->Order[egrp] != pOrder){ */
/*       pOrder = U->Order[egrp]; */
/*       ierr = xf_Error(xf_GetQuadOrderElem(Mesh, All->EqnSet, egrp, U->Order[egrp],  */
/* 					  &QuadOrder)); */
/*       if (ierr != xf_OK) return ierr; */
/*       ierr = xf_Error(xf_QuadLine(QuadOrder, &nq0, &xq0, &wq0)); */
/*       if (ierr != xf_OK) return ierr; */
/*       ierr = xf_Error(xf_ReAlloc( (void **) &xqglob, nq0, sizeof(int))); */
/*       if (ierr != xf_OK) return ierr; */
/*     } */

/*     // loop over elements */
/*     for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){ */
/*       // call element intersection function */
/*       ierr = xf_Error(xf_IntersectElemWithLine(Mesh, egrp, elem, LineCoord, &PhiData1,  */
/* 					       &nseg, xref, xglob)); */
/*       if (ierr != xf_OK) return ierr; */

/*       if (nseg <= 0) continue; // no intersections for this element */

/*       ielem = CPI->nelem++; */

/*       // reallocate CPI */
/*       ierr = xf_Error(xf_ReAlloc((void **) &CPI->egrp, CPI->nelem, sizeof(int))); */
/*       if (ierr != xf_OK) return ierr; */
/*       ierr = xf_Error(xf_ReAlloc((void **) &CPI->elem, CPI->nelem, sizeof(int))); */
/*       if (ierr != xf_OK) return ierr; */
/*       ierr = xf_Error(xf_ReAlloc((void **) &CPI->nquad, CPI->nelem, sizeof(int))); */
/*       if (ierr != xf_OK) return ierr; */
/*       ierr = xf_Error(xf_ReAlloc((void **) &CPI->xquad, CPI->nelem, sizeof(real *))); */
/*       if (ierr != xf_OK) return ierr; */
/*       ierr = xf_Error(xf_ReAlloc((void **) &CPI->wquad, CPI->nelem, sizeof(real *))); */
/*       if (ierr != xf_OK) return ierr; */

/*       // set info for current elem */
/*       CPI->egrp[ielem]  = egrp; */
/*       CPI->elem[ielem]  = elem; */
/*       CPI->nquad[ielem] = nqtot = nq0*nseg; */
/*       ierr = xf_Error(xf_Alloc((void **) &CPI->xquad[ielem], nqtot*dim, sizeof(real))); */
/*       if (ierr != xf_OK) return ierr; */
/*       ierr = xf_Error(xf_Alloc((void **) &CPI->wquad[ielem], nqtot, sizeof(real))); */
/*       if (ierr != xf_OK) return ierr; */

/*       // loop over segments */
/*       for (iseg=0; iseg<nseg; iseg++){ */
/* 	// quad points in global space */
/* 	for (iq=0; iq<nq0; iq++) */
/* 	  for (d=0, s=0.5*(1.0+xq0[iq]); d<dim; d++) */
/* 	    xqglob[iq*dim+d] = xglob[2*dim*iseg+d]*(1.0-s)+xglob[2*dim*iseg+dim+d]*(1.0+s); */

/* 	// convert quad points in reference space and store them */
/* 	ierr = xf_Error(xf_Glob2RefElem(Mesh, egrp, elem, &PhiData2, xfe_True,  */
/* 					nnode, xqglob, CPI->xquad[ielem]+iseg*nq0)); */
/* 	if (ierr != xf_OK) return ierr; */

/* 	// segment length and quad weights */
/* 	slen = 0.; */
/* 	for (d=0; d<dim; d++){ */
/* 	  s = xglob[2*dim*iseg+dim+d]-xglob[2*dim*iseg+d]; */
/* 	  slen += s*s; */
/* 	} */
/* 	slen = sqrt(slen); */
/* 	for (iq=0; iq<nq0; iq++) CPI->wquad[ielem][iseg*nq0+iq]=0.5*slen*wq[iq]; */

/*       } // iseg */
/*     } // elem */
/*   } // egrp */

/*   /\* Destroy Basis Data *\/ */
/*   ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */

/*   ierr = xf_Error(xf_DestroyBasisData(PhiData2, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */


/*   xf_Release( (void *) xq0); */
/*   xf_Release( (void *) wq0); */
/*   xf_Release( (void *) xqglob); */


/*   return xf_OK; */
/* } */



/******************************************************************/
//   FUNCTION Definition: xf_OutputLineCutPlaneIntegral
static int 
xf_OutputLineCutPlaneIntegral( xf_All *All, xf_Output *Output, const xf_Vector *U, 
                              real weight, real *Value, xf_Vector *Value_U, 
                              enum xfe_AddType AddFlag)
{
  int ierr, k, i, sr, egrp, elem, dim, d, ifac;
  int Order, nn, iq, nq, pnq, ielem;
  int *IParam;
  enum xfe_Bool VectorFlag;
  enum xfe_BasisType Basis;
  real *RParam, *xq, *wq, *u, *s, *s_u, *gu;
  real *EU, *EV_U, LocValue, N[3], dp;
  xf_BasisData *PhiData;
  xf_CutPlaneIntersect *CPI;
  xf_EqnSet *EqnSet;
  xf_JacobianData *JData; 
  
  if (Output->UsesFlux) return xf_Error(xf_NOT_SUPPORTED);
  
  dim = All->Mesh->Dim;
  
  // does this output integrate a vector?
  VectorFlag = (Output->VectorName != NULL);
  
  EqnSet = All->EqnSet;
  sr = EqnSet->StateRank;
  
  if (Value_U != NULL){
    // zero out Value_U if requesting a set or neg
    if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
      ierr = xf_Error(xf_SetZeroVector(Value_U));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // initialize LocValue
  LocValue = 0.0;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // calculate line or cut-plane intersection if do not have it
  if (Output->CutPlaneIntersect == NULL){
    if (Output->Type == xfe_CutPlaneIntegral){
      ierr = xf_Error(xf_CalculateCutPlaneIntersect(All, Output->CutPlane, U,
                                                    &Output->CutPlaneIntersect));
      if (ierr != xf_OK) return ierr;
    }
    else if (Output->Type == xfe_LineIntegral){
      return xf_Error(xf_NOT_SUPPORTED);
      /*       ierr = xf_Error(xf_CalculateCutLineIntersect(All, Output->LineCoord, U, */
      /* 						   &Output->CutPlaneIntersect)); */
      /*       if (ierr != xf_OK) return ierr;       */
    }
  }
  
  // for easy typing
  CPI = Output->CutPlaneIntersect;
  
  // loop over intersected elements
  PhiData = NULL;
  u       = NULL;
  gu      = NULL;
  s       = NULL;
  s_u     = NULL;
  pnq     = -1;
  for (ielem=0; ielem<CPI->nelem; ielem++){
    
    JData = NULL;
    
    // elem info
    egrp = CPI->egrp[ielem];
    elem = CPI->elem[ielem];
    
    // basis and order of state vector on egrp
    Basis = U->Basis[egrp]; 
    Order = xf_InterpOrder(U, egrp, elem);
    
    // quad point information
    nq = CPI->nquad[ielem];
    xq = CPI->xquad[ielem];
    wq = CPI->wquad[ielem];
    
    // compute basis functions (and grads) if quad or basis or order changed
    ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xq, xfb_Phi | xfb_GPhi | xfb_gPhi, &PhiData));
    
    if (ierr != xf_OK) return ierr;
    
    nn = PhiData->nn; // number of interpolation nodes
    
    // re-allocate data if quad points increased
    if (nq > pnq){
      ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc( (void **)  &gu, dim*nq*sr, sizeof(real))); 
      if (ierr != xf_OK) return ierr;
      ifac = ((VectorFlag) ? dim : 1);
      ierr = xf_Error(xf_ReAlloc( (void **)  &s, nq*ifac, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      if (Value_U != NULL){
        ierr = xf_Error(xf_ReAlloc( (void **)  &s_u, nq*sr*ifac, sizeof(real)));
        if (ierr != xf_OK) return ierr;
      }
    }
    
    EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
    
    // interpolate state at quad points
    xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u);
    
    /* Compute geometry Jacobian */
    ierr = xf_Error(xf_ElemJacobian(All->Mesh, egrp, elem, nq, xq, 
                                    xfb_detJ | xfb_iJ | xfb_J, xfe_True, &JData));
    if (ierr != xf_OK) return ierr;
    
    /* convert reference basis grads (GPhi) to physical grads, gPhi */
    ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData)); 
    if (ierr != xf_OK) return ierr;
    
    for (d=0; d<dim; d++)
      xf_MxM_Set(PhiData->gPhi+nn*nq*d, EU, nq, nn, sr, gu + nq*sr*d);
    
    if (VectorFlag){
      // normal
      for (d=0; d<dim; d++) N[d] = Output->CutPlane[d];
      for (d=0, dp=0.; d<dim; d++) dp += N[d]*N[d];
      dp = sqrt(dp);
      for (d=0.; d<dim; d++) N[d] /= dp;
      
      // call eqnset specific function for vector
      ierr = xf_Error(xf_EqnSetVector(EqnSet, Output->VectorName, IParam, 
                                      RParam, nq, u, s, s_u));
      if (ierr != xf_OK) return ierr;
      for (iq=0; iq<nq; iq++){
        for (d=0, dp=0.; d<dim; d++) dp += N[d]*s[iq*dim+d];
        LocValue += wq[iq] * dp;
        if (Value_U != NULL){ // dot with normal, do not include wq yet
          for (k=0; k<sr; k++) s_u[iq*sr+k] = s_u[iq*dim*sr + 0*sr+k]*N[0];
          for (d=1; d<dim; d++)
            for (k=0; k<sr; k++) s_u[iq*sr+k] += s_u[iq*dim*sr + d*sr + k]*N[d];
        }
      }
      // Value_U{n,k} @= int(s_u{q,k} * Phi{q,n}) :: weight included here
      if (Value_U != NULL){
        xf_ColcMult(s_u, wq, nq, sr, 1, weight);// multiply s_u by wq*weight
        EV_U = Value_U->GenArray[egrp].rValue[elem]; // Value_U on elem [nn*sr]
        xf_MTxM(PhiData->Phi, s_u, nn, nq, sr, AddFlag, EV_U);
      }
      
    }
    else{
      // call eqnset specific function for scalar
      ierr = xf_Error(xf_EqnSetScalar(EqnSet, Output->ScalarName, IParam, 
                                      RParam, nq, u, gu, s, s_u, NULL, NULL,0.0));
      if (ierr != xf_OK) return ierr;
      
      if (Output->DomainNorm == xfe_DomainNormL2){
        for (iq=0; iq<nq; iq++)
          LocValue += wq[iq]*s[iq]*s[iq]; //integrate square of scalar
      }
      else{
        // sum scalar*wq over quad points, add to LocValue
        for (iq=0; iq<nq; iq++)
          LocValue += wq[iq]*s[iq];
      }
      
      // Value_U{n,k} @= int(s_u{q,k} * Phi{q,n}) :: weight included here
      if (Value_U != NULL){
        xf_ColcMult(s_u, wq, nq, sr, 1, weight);// multiply s_u by wq*weight
        EV_U = Value_U->GenArray[egrp].rValue[elem]; // Value_U on elem [nn*sr]
        xf_MTxM(PhiData->Phi, s_u, nn, nq, sr, AddFlag, EV_U);
      }
    }
    
    /* Destroy geometry Jacobian Data */
    ierr = xf_Error(xf_DestroyJacobianData(JData));
    if (ierr != xf_OK) return ierr;
    
    pnq = nq;
  } // ielem
  
  /* Set Value :: weight included here */
  if (Value != NULL) (*Value) = LocValue * weight;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void  *) IParam);
  xf_Release( (void  *) RParam);
  xf_Release( (void *) u);
  xf_Release( (void *) s);
  xf_Release( (void *) s_u);
  xf_Release( (void *) gu);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_GetLocalElem
int 
xf_GetLocalElem(xf_Mesh *Mesh, int egrpglob, int elemglob, 
                int *egrploc, int *elemloc)
{
  /* Calculates local element corresponding to global element
   (egrp,elem).  Only difference will be in parallel. -1 is returned
   in egrploc and elemloc if the global element does not exist on
   this processor.  This function uses a brute-force search and
   should not be called often. */
  
  int ierr, elem;
  int myRank, nProc;
  int **EL2G;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc == 1){ // return immediately if serial
    (*egrploc) = egrpglob;
    (*elemloc) = elemglob;
    return xf_OK;
  }
  
  if ((Mesh->ParallelInfo == NULL) ||
      (EL2G = Mesh->ParallelInfo->ElemLoc2Glob) == NULL)
    return xf_Error(xf_PARALLEL_ERROR);
  
  (*egrploc) = -1;
  (*elemloc) = -1;
  for (elem=0; elem<Mesh->ElemGroup[egrpglob].nElem; elem++)
    if (EL2G[egrpglob][elem] == elemglob){
      (*egrploc) = egrpglob;
      (*elemloc) = elem;
      break;
    }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_OutputPointValue
static int 
xf_OutputPointValue( xf_All *All, xf_Output *Output, const xf_Vector *U, 
                    real weight, real *Value, xf_Vector *Value_U, 
                    enum xfe_AddType AddFlag)
{
  int ierr, k, i, n, sr, egrp, elem;
  int Order, nn, dim;
  int *IParam;
  real Gamma;
  enum xfe_BasisType Basis;
  enum xfe_Bool MotionOn = xfe_False;
  real *RParam, *xq, *wq, *u, *s, *s_u;
  real *xref, *EU, *EV_U;
  real Time, xglob[3], dp;
  xf_BasisData *PhiData;
  xf_Vector *J_GCL, *GammaVec;
  xf_Data *GammaDat;
  xf_EqnSet *EqnSet;
  xf_Mesh *Mesh;
  xf_MotionData *MD = NULL;
  
  if (Output->UsesFlux) return xf_Error(xf_INPUT_ERROR);
 
  ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
  if(ierr == xf_NOT_FOUND)
  {
     xf_printf("Cannot find heat capacity ratio...\n");
     return ierr;
  }
  else
     GammaVec = (xf_Vector *) GammaDat->Data;


  EqnSet = All->EqnSet;
  sr     = EqnSet->StateRank;
  Mesh   = All->Mesh;
  dim    = Mesh->Dim;
  
  if (Value_U != NULL){
    // zero out Value_U if requesting a set or neg
    if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
      ierr = xf_Error(xf_SetZeroVector(Value_U));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // determine if we need mesh motion
  MotionOn = ((Mesh->Motion != NULL) && (Mesh->Motion->Active));
  if (MotionOn){
    ierr = xf_Error(xf_CreateMotionData(All, &MD));
    if (ierr != xf_OK) return ierr;
  }
  
  // locate J_GCL (zero out if AddFlag is set/neg)
  if (Value_U != NULL){ // only if also calculating Value_U
    ierr = xf_Error(xf_FindOutputGCLLinearization(All, Output->Name, AddFlag, &J_GCL));
    if (ierr != xf_OK) return ierr;
  }
  else J_GCL = NULL;
  
  // determine Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
  if (ierr != xf_OK) return ierr;
  
  // initialize vars to NULL
  PhiData = NULL;
  u       = NULL;
  s       = NULL;
  s_u     = NULL;
  
  /* Calculate element info.  May require a global-to-local map, but
   this is done just once. */
  if (Output->elemLocal == NULL){
    ierr = xf_Error(xf_Alloc( (void **) &Output->elemLocal, 2, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_GetLocalElem(All->Mesh, Output->egrp, Output->elem, 
                                    Output->elemLocal + 0, Output->elemLocal + 1));
    if (ierr != xf_OK) return ierr;
  }
  
  // element information
  egrp = Output->elemLocal[0];
  elem = Output->elemLocal[1];
  xref = Output->xref;
  
  if (egrp < 0){ // this is the case in parallel if loc is not on proc
    if (Value != NULL) (*Value) = 0.0;
    return xf_OK;
  }
  
  // need global coordinate of point if motion is on
  if (MotionOn){
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True, 1, xref, xglob));
    if (ierr != xf_OK) return ierr;
  }
  
  // basis and order of state vector on egrp
  Basis = U->Basis[egrp]; 
  Order = xf_InterpOrder(U, egrp, elem);
  
  // compute basis functions at xref
  ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, 1, xref, xfb_Phi, &PhiData));
  if (ierr != xf_OK) return ierr;
  
  nn = PhiData->nn; // number of interpolation nodes
  
  // obtain transformation map if doing mesh motion
  if (MotionOn){
    ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, Mesh->Motion, 
                                     1, dim, Time, xglob, MD));
    if (ierr != xf_OK) return ierr;
  }
  
  // allocate data
  ierr = xf_Error(xf_ReAlloc( (void **)  &u, sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **)  &s, 1, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  if (Value_U != NULL){
    ierr = xf_Error(xf_ReAlloc( (void **)  &s_u, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  
  EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
  Gamma = GammaVec->GenArray[egrp].rValue[elem][0];

  // interpolate state at quad points
  xf_MxM_Set(PhiData->Phi, EU, 1, nn, sr, u);
  
  // transform state to physical
  if (MotionOn) xf_ModMotionPreEqnCall(1, dim, sr, MD, u, NULL);
  
  // call eqnset specific function for scalar
  ierr = xf_Error(xf_EqnSetScalar(EqnSet, Output->ScalarName, IParam,
                                  RParam, 1, u, NULL, s, s_u, NULL, NULL, Gamma));
  if (ierr != xf_OK) return ierr;
  
  // Divide output linearization by gbar to make it wrt reference state
  if ((MotionOn) && (Value_U != NULL) && (s_u != NULL))
    xf_ColDiv(s_u, MD->gb, 1, sr, 1);
  
  if (Value != NULL){
    if (Output->DomainNorm == xfe_DomainNormL2){
      (*Value) = s[0]*s[0]*weight; //square of scalar
    }
    else if (Output->DomainNorm == xfe_DomainNormNone) (*Value) = s[0]*weight;
    // set Value :: include weight here
  }
  
  // Value_U{n,k} @= s_u{1,k} * Phi{1,n}
  if (Value_U != NULL){
    EV_U = Value_U->GenArray[egrp].rValue[elem]; // Value_U on elem [nn*sr]
    for (k=0; k<sr; k++) s_u[k] *= weight; // :: include weight here
    xf_MTxM(PhiData->Phi, s_u, nn, 1, sr, AddFlag, EV_U);
    if (J_GCL != NULL){
      // J_GCL{n} = J_ubar{k} dot (-u{k}) * Phi{1,n}
      // note, u{k} is the physical state
      xf_DotProduct(s_u, u, sr, &dp);
      xf_cV_Add(PhiData->Phi, -dp, nn, AddFlag, J_GCL->GenArray[egrp].rValue[elem]);
    }
  }
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy mesh motion data */
  xf_DestroyMotionData(MD);
  
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) u);
  xf_Release( (void *) s);
  xf_Release( (void *) s_u);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FindSensitivity
int 
xf_FindSensitivity(xf_Output *Output, xf_Sensitivity **pSensitivity, 
                   enum xfe_SensitivityType SensType, char *ParamName, 
                   int nvar, xf_KeyValue KeyValue, 
                   enum xfe_Bool *Found)
{
  int ierr, i, nmatch;
  enum xfe_Bool found = xfe_False, PMatch, SMatch, VMatch;
  
  /* for (i = 0; i < Output->nSensitivity; i++) {
   <#statements#>
   } */
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FindOutput
int 
xf_FindOutput( const xf_EqnSet *EqnSet, const char *Name, xf_Output **pOutput)
{
  int i;
  enum xfe_Bool found;
  xf_Output *Output;
  
  // no outputs -> return error
  if (EqnSet->Outputs == NULL) return xf_NOT_FOUND;
  
  Output = EqnSet->Outputs->Output;
  
  // find output with correct Name
  for (i=0, found = xfe_False; i<EqnSet->Outputs->nOutput; i++){
    Output = EqnSet->Outputs->Output + i;
    if (strcmp(Name, Output->Name) == 0){
      found = xfe_True;
      break;
    }
  } // i
  if (!found) return xf_NOT_FOUND;
  
  (*pOutput) = Output;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_IsOutputUnsteady
int 
xf_IsOutputUnsteady( const xf_EqnSet *EqnSet, const char *Name, 
                    enum xfe_Bool *IsUnsteady)
{
  int ierr;
  xf_Output *Output;
  
  ierr = xf_Error(xf_FindOutput(EqnSet, Name, &Output));
  if (ierr != xf_OK) return ierr;
  
  (*IsUnsteady) = (Output->TimeNorm != xfe_TimeNormNone);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadStoredOutput
int 
xf_ReadStoredOutput( const xf_EqnSet *EqnSet, const char *Name, 
                    real *Value, real *ErrEst)
{
  int ierr;
  xf_Output *Output;
  
  ierr = xf_Error(xf_FindOutput(EqnSet, Name, &Output));
  if (ierr != xf_OK) return ierr;
  
  if (Value  != NULL) (*Value ) = Output->Value;
  if (ErrEst != NULL) (*ErrEst) = Output->ErrEst;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateOutput
int 
xf_CalculateOutput( xf_All *All, const char *Name, const xf_Vector *U, 
                   real *Value, xf_Vector *Value_U, enum xfe_AddType AddFlag1)
{
  int ierr, nOutput, iOutput;
  xf_Output *Output0, *Output;
  enum xfe_AddType AddFlag, AddFlag2;
  real weight, LocValue=0.;
//  Yu_Model Model;
//  xf_Vector *GammaVec;
//  xf_Data   *GammaDat;

//  ierr = xf_Error(PullinModel(&Model));
//  if (ierr != xf_OK) return ierr;

//  ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, NULL));
//  if(ierr != xf_OK) return ierr;

//  ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
//  if(ierr == xf_NOT_FOUND)
//  {
//     xf_printf("Cannot find heat capacity ratio...\n");
//     return ierr;
//  }
//  else
//     GammaVec = (xf_Vector *) GammaDat->Data;
//  ierr = xf_Error(Yu_GammaVectorUpdate(All, &Model, GammaVec, U));
//  if (ierr != xf_OK) return ierr;

  if (Value != NULL) (*Value) = 0.0;
  
  ierr = xf_Error(xf_FindOutput(All->EqnSet, Name, &Output0));
  if (ierr != xf_OK) return ierr;
  
  nOutput = ((Output0->Type == xfe_SumOutput) ? Output0->nSumOutput : 1);
  
  // zero out Value
  if (Value != NULL) (*Value) = 0.0;
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag1);
  
  for (iOutput=0; iOutput<nOutput; iOutput++){
    
    // determine whether we're setting or adding to Value_U
    AddFlag = ((iOutput == 0) ? AddFlag1 : AddFlag2);
    
    // pull off desired output (only different if summing outputs)
    if (Output0->Type == xfe_SumOutput){
      ierr = xf_Error(xf_FindOutput(All->EqnSet, Output0->SumOutputNames[iOutput], 
                                    &Output));
      if (ierr != xf_OK) return ierr;
      weight = Output0->SumOutputWeights[iOutput];
    }
    else{
      Output = Output0;
      weight = 1.0;
    }
    
    switch (Output->Type){
      case xfe_DomainIntegral:
        ierr = xf_Error(xf_OutputDomainIntegral(All, Output, U, weight, &LocValue, 
                                                Value_U, AddFlag));
        if (ierr != xf_OK) return ierr;
        break;
        
      case xfe_BoundaryIntegral:
        ierr = xf_Error(xf_OutputBoundaryIntegral(All, Output, U, weight, &LocValue, 
                                                  Value_U, AddFlag));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_MPI_Allreduce(&LocValue, 1, xfe_SizeReal, xfe_MPI_SUM));
        if (ierr != xf_OK) return ierr;
        break;
        
      case xfe_CutPlaneIntegral:
      case xfe_LineIntegral:
        ierr = xf_Error(xf_OutputLineCutPlaneIntegral(All, Output, U, weight, &LocValue, 
                                                      Value_U, AddFlag));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_MPI_Allreduce(&LocValue, 1, xfe_SizeReal, xfe_MPI_SUM));
        if (ierr != xf_OK) return ierr;
        break;
        
      case xfe_PointValue:
        ierr = xf_Error(xf_OutputPointValue(All, Output, U, weight, &LocValue, 
                                            Value_U, AddFlag));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_MPI_Allreduce(&LocValue, 1, xfe_SizeReal, xfe_MPI_SUM));
        if (ierr != xf_OK) return ierr;
        break;
        
      default:
        return xf_Error(xf_NOT_SUPPORTED);
        break;
    }
    
    
    // add to requested value (weight was included in individual functions)
    if (Value != NULL) (*Value) += LocValue;
    
  } // iOutput
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_IncrementUnsteadyOutputs
int 
xf_IncrementUnsteadyOutputs( xf_All *All, const char *Name, 
                            enum xfe_TimeSchemeType TimeScheme, int nU,
                            xf_Vector **Ui, xf_Vector *Utemp, real Time, 
                            real TimeStep, enum xfe_Bool AtStart, 
                            enum xfe_Bool AtEnd, real *Value, 
                            xf_Vector **Value_Ui)
{
  int ierr, nOutput, iOutput, iU;
  int iq, nq, OrderTime, QuadOrder;
  int nSumOutput, iSum;
  enum xfe_Bool MultiStep;
  enum xfe_Bool UseGCL;
  enum xfe_TimeNormType TimeNorm;
  char LogOutput[xf_MAXLINELEN];
  const char *OutputStrings;
  char **OutputList = NULL;
  real LocValue = 0., rtemp, wtemp, weight;
  real OutputTime, xt, xtstart, xtend;
  real TimeOrig;
  real *tq = NULL, *wq = NULL;
  real phi[xf_MAXDGTIMENODE];
  xf_Output *Output, *Output0;
  xf_Vector *Value_Utemp = NULL;
  xf_Vector *J_GCL = NULL;
  xf_Vector **J_GCLi = NULL;
  
  if (Name == NULL){
    // use LogOutput
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "LogOutput", LogOutput));
    if (ierr != xf_OK) return ierr;
    OutputStrings = LogOutput;
  }
  else OutputStrings = Name;
  
  if (xf_NotNull(OutputStrings)){
    OutputList = NULL;
    ierr = xf_Error(xf_ScanXStringAlloc(OutputStrings, xf_MAXSTRLEN, 
                                        &nOutput, &OutputList));
    if (ierr != xf_OK) return ierr;
  }
  else return xf_OK; // nothing to do
  
  // Is TimeScheme MultiStep?
  ierr = xf_Error(xf_GetTimeSchemeInfo(TimeScheme, NULL, &MultiStep, NULL));
  if (ierr != xf_OK) return ierr;
  
  // check number of U vectors passed in
  if (MultiStep && (nU != 1)) return xf_Error(xf_INPUT_ERROR);
  if ((TimeScheme == xfe_TimeSchemeDG1) && (nU != 2)) return xf_Error(xf_INPUT_ERROR);
  
  // This function does not ultimately alter "Time" so as not to disrupt calling function
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &TimeOrig));
  if (ierr != xf_OK) return ierr;
  
  // Determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  
  
  // loop over outputs
  for (iOutput=0; iOutput<nOutput; iOutput++){
    
    // find output structure
    ierr = xf_Error(xf_FindOutput(All->EqnSet, OutputList[iOutput], &Output0));
    if (ierr != xf_OK) return ierr;
    
    // initialize Output0 if at start of time integration (only if Value != NULL)
    if ((Value != NULL) && AtStart) Output0->Value = 0.;
    
    // number of sub-outputs to deal with sumoutputs
    nSumOutput = ((Output0->Type == xfe_SumOutput) ? Output0->nSumOutput : 1);
    
    // allocate/pull off J_GCLi if using GCL and need derivatives
    if ((UseGCL) && (Value_Ui != NULL)){
      
      // Generic J_GCL needs to exist
      ierr = xf_Error(xf_FindMeshMotionGCLLinearization(All, Output0->Name, -1, &J_GCL));
      if (ierr != xf_OK) return ierr;
      
      // allocate J_GCLi
      ierr = xf_Error(xf_Alloc( (void **) &J_GCLi, nU, sizeof(xf_Vector *)));
      if (ierr != xf_OK) return ierr;
      
      for (iU=0; iU<nU; iU++){
        // J_GCLi already needs to exist
        ierr = xf_Error(xf_FindMeshMotionGCLLinearization(All, Output0->Name, iU, J_GCLi+iU));
        if (ierr != xf_OK) return ierr;
        
        // set J_GCLi to zero
        ierr = xf_Error(xf_SetZeroVector(J_GCLi[iU]));
        if (ierr != xf_OK) return ierr;
      } // iU
      
    }
    else{
      J_GCL = NULL;
      J_GCLi = NULL;
    }
    
    // initialize Value_U to zero
    if (Value_Ui != NULL){
      for (iU=0; iU<nU; iU++){
        ierr = xf_Error(xf_SetZeroVector(Value_Ui[iU]));
        if (ierr != xf_OK) return ierr;
      }
    }
    
    // loop over sub-outputs
    for (iSum=0; iSum<nSumOutput; iSum++){
      
      // pull off desired output (only different if summing outputs)
      if (Output0->Type == xfe_SumOutput){
        ierr = xf_Error(xf_FindOutput(All->EqnSet, Output0->SumOutputNames[iSum],
                                      &Output));
        if (ierr != xf_OK) return ierr;
        weight = Output0->SumOutputWeights[iSum];
      }
      else{
        Output = Output0;
        weight = 1.0;
      }
      
      // Only considering unsteady outputs
      if ((TimeNorm = Output->TimeNorm) == xfe_TimeNormNone) continue;
      
      LocValue = 0.;
      
      switch (TimeNorm){
        case xfe_TimeNormNone:
          continue; // only considering unsteady outputs
          break;
          
        case xfe_TimeNormFinal:
          if (AtEnd){ // only calculate when at end of simulation
            
            // set Time to end time
            ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time + TimeStep));
            if (ierr != xf_OK) return ierr;
            
            // If DG with GCL, interpolate GCL to end of time slab
            if( (UseGCL) && ((TimeScheme == xfe_TimeSchemeDG1)||(TimeScheme == xfe_TimeSchemeDG2))){
              ierr = xf_Error(xf_DGTimeInterpolateGCL(All, TimeScheme, 1., NULL));
              if (ierr != xf_OK) return ierr;
            }
            
            // calculate output (note: for DG with GCL, GCL should already be interpolated to end of time slab)
            ierr = xf_Error(xf_CalculateOutput(All, Output->Name, Ui[nU-1], &LocValue, 
                                               ((Value_Ui == NULL) ? NULL : Utemp), xfe_Set));
            if (ierr != xf_OK) return ierr;
            // include weight in linearization
            if (Value_Ui != NULL){
              ierr = xf_Error(xf_VectorMultSet(Utemp, weight, xfe_Add, Value_Ui[nU-1]));
              if (ierr != xf_OK) return ierr;
            }
            if (J_GCLi != NULL){ // same with GCL
              ierr = xf_Error(xf_VectorMultSet(J_GCL, weight, xfe_Add, J_GCLi[nU-1]));
              if (ierr != xf_OK) return ierr;
            }
          }
          break;
          
        case xfe_TimeNormPoint:
          if (MultiStep) return xf_Error(xf_NOT_SUPPORTED);
          
          // Determine if output time falls into time slab
          OutputTime = Output->StartTime;
          if ((OutputTime >= Time) && (OutputTime < Time+TimeStep)){
            
            // set Time to OutputTime
            ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", OutputTime));
            if (ierr != xf_OK) return ierr;
            
            xt = (OutputTime-Time)/TimeStep; // 1 position
            // Interpolate state
            ierr = xf_Error(xf_DGTimeInterpolateState(All, TimeScheme, Ui, &xt, &Utemp, phi));
            if (ierr != xf_OK) return ierr;
            // need a temporary vector if Value_Ui is requested
            if (Value_Ui != NULL){
              ierr = xf_Error(xf_FindSimilarVector(All, Value_Ui[0], "Value_Utemp", 
                                                   xfe_True, xfe_True, 
                                                   NULL, &Value_Utemp, NULL));
              if (ierr != xf_OK) return ierr;
            }
            // Calculate the output
            ierr = xf_Error(xf_CalculateOutput(All, Output->Name, Utemp, &LocValue, 
                                               ((Value_Ui == NULL) ? NULL : Value_Utemp), xfe_Set));
            if (ierr != xf_OK) return ierr;
            // distribute Value_Utemp into Value_Ui appropriately
            if (Value_Ui != NULL){
              // chain rule: Value_Ui[i] = Value_Utemp * Utemp_U[i]; include weight
              for (iU=0; iU<nU; iU++){
                ierr = xf_Error(xf_VectorMultSet(Value_Utemp, phi[iU]*weight, xfe_Add, Value_Ui[iU]));
                if (ierr != xf_OK) return ierr;
              }
            }
            // same with GCL
            if (J_GCLi != NULL){
              for (iU=0; iU<nU; iU++){
                ierr = xf_Error(xf_VectorMultSet(J_GCL, phi[iU]*weight, xfe_Add, J_GCLi[iU]));
                if (ierr != xf_OK) return ierr;
              }
            }
          }
          break;
          
        case xfe_TimeNormIntegral:
        case xfe_TimeNormSquareIntegral:
          if (MultiStep) return xf_Error(xf_NOT_SUPPORTED);
          
          // Determine if output time interval falls into time slab
          if ((Output->StartTime <  Time+TimeStep) && 
              (Output->EndTime   >  Time         )){
            
            // time interval = [xtstart, xtend] in [0,1]
            xtstart = ((Output->StartTime <  Time         ) ? 0. : 
                       (Output->StartTime-Time)/TimeStep);
            xtend   = ((Output->EndTime   >  Time+TimeStep) ? 1. : 
                       (Output->EndTime-Time)/TimeStep);
            // determine order in time
            ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
            if (ierr != xf_OK) return ierr;
            // quad order is semi-hardcoded
            QuadOrder = OrderTime+1;
            if (TimeNorm == xfe_TimeNormSquareIntegral) QuadOrder *= 2;
            // quad points
            tq = NULL;
            wq = NULL;
            ierr = xf_Error(xf_QuadLine(QuadOrder, &nq, &tq, &wq));
            if (ierr != xf_OK) return ierr;
            
            // need a temporary vector if Value_Ui is requested
            if (Value_Ui != NULL){
              ierr = xf_Error(xf_FindSimilarVector(All, Value_Ui[0], "Value_Utemp", 
                                                   xfe_True, xfe_True, 
                                                   NULL, &Value_Utemp, NULL));
              if (ierr != xf_OK) return ierr;
            }
            
            // loop over quad points
            for (iq=0; iq<nq; iq++){
              // quadrature weight at this point
              wtemp =  wq[iq]*TimeStep*(xtend-xtstart);
              xt    = xtstart + tq[iq]*(xtend-xtstart);
              // set Time to current time at quad point
              ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time+xt*TimeStep));
              if (ierr != xf_OK) return ierr;
              // Interpolate state
              ierr = xf_Error(xf_DGTimeInterpolateState(All, TimeScheme, Ui, &xt, &Utemp, phi));
              if (ierr != xf_OK) return ierr;
              // Calculate the output at Utemp 
              ierr = xf_Error(xf_CalculateOutput(All, Output->Name, Utemp, &rtemp, 
                                                 ((Value_Ui == NULL) ? NULL : Value_Utemp),
                                                 xfe_Set));
              if (ierr != xf_OK) return ierr;
              // need output squared if doing square integral
              if (TimeNorm == xfe_TimeNormSquareIntegral) wtemp *= rtemp; 
              // Add to running total:
              LocValue += wtemp*rtemp;
              // distribute Value_Ui appropriately
              if (Value_Ui != NULL){
                // chain rule: Value_Ui[i] += wtemp*Value_Utemp * Utemp_U[i]
                //  (square) : Value_Ui[i] += 2*wtemp*Value_Utemp * Utemp_U[i]
                //             (note wtemp contains rtemp by above)
                if (TimeNorm == xfe_TimeNormSquareIntegral) wtemp *= 2.; 
                for (iU=0; iU<nU; iU++){
                  ierr = xf_Error(xf_VectorMultSet(Value_Utemp, weight*wtemp*phi[iU], 
                                                   xfe_Add, Value_Ui[iU]));
                  if (ierr != xf_OK) return ierr;
                }
                // same with GCL
                if (J_GCLi != NULL){
                  for (iU=0; iU<nU; iU++){
                    ierr = xf_Error(xf_VectorMultSet(J_GCL, weight*wtemp*phi[iU], 
                                                     xfe_Add, J_GCLi[iU]));
                    if (ierr != xf_OK) return ierr;
                  }
                }
              }
            } // iq
            
            xf_Release( (void *) tq);
            xf_Release( (void *) wq);
            
          }
          
          break;
          
        default:
          return xf_Error(xf_NOT_SUPPORTED);
          break;
      }
      
      if (Value != NULL){
        Output0->Value += LocValue*weight;
        (*Value) = Output0->Value;
      }
      
    } // iSum
    
    xf_Release( (void *) J_GCLi);
    
  } // iOutput
  
  // Reset time to original value when this function was called
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", TimeOrig));
  if (ierr != xf_OK) return ierr;
  
  
  xf_Release2( (void **) OutputList);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PingCompareFDtoAnalytical
static int
xf_PingCompareFDtoAnalytical(real fd, real fa, real ep, 
                             enum xfe_Verbosity Verbosity)
{
  if (Verbosity != xfe_VerbosityLow){
    if (fabs(fd-fa) > 1000.0*ep*ep)
      xf_printf("<-----***\n");	  
    else if (fabs(fd-fa) > 100.0*ep*ep)
      xf_printf("<-----**\n");
    else if (fabs(fd-fa) > 10.0*ep*ep)
      xf_printf("<-----*\n");
    else
      xf_printf("\n");
  }
  if (fabs(fd-fa) > 10.0*ep*ep) return xf_Error(xf_PING_FAILED);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PingOutput
int
xf_PingOutput(xf_All *All, const char *Name, xf_Vector *U)
{
  int ierr, j, egrp, elem;
  int sr, r;
  enum xfe_Verbosity Verbosity;
  real ep, fd, fa, s0, *s_U0, s, *s_U;
  xf_Mesh *Mesh;
  xf_Vector *Value_U0, *Value_U;
  xf_GenArray *ga;
  
  ep = 1.0e-5;
  
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;
  
  // Determine verbosity for printing purposes
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", xfe_VerbosityName, 
                                     (int ) xfe_VerbosityLast, (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  // locate vectors
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Value_U0", xfe_False, xfe_False, 
                                       NULL, &Value_U0, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Value_U", xfe_False, xfe_False, 
                                       NULL, &Value_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  // calculate s0 and s_U0
  ierr = xf_Error(xf_CalculateOutput(All, Name, U, &s0, Value_U0, xfe_Set));
  if (ierr != xf_OK) return ierr;
  
  if (Verbosity != xfe_VerbosityLow)
    xf_printf("s0 = %.15E\n", s0);
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      ga = U->GenArray+egrp;
      r = ((ga->vr==NULL) ? ga->r : ga->vr[elem]);
      for (j=0; j<r; j++){
        // increment U
        U->GenArray[egrp].rValue[elem][j] += ep;
        
        // calculate s and s_U
        ierr = xf_Error(xf_CalculateOutput(All, Name, U, &s, Value_U, xfe_Set));
        if (ierr != xf_OK) return ierr;
        
        if (Verbosity != xfe_VerbosityLow)
          xf_printf("s = %.15E, s-s0 = %.15E\n", s, s-s0);
        
        // check Delta s vs s_U
        fd = (s - s0)/ep;
        fa = 0.5*(Value_U->GenArray[egrp].rValue[elem][j] +
                  Value_U0->GenArray[egrp].rValue[elem][j]);
        if (Verbosity != xfe_VerbosityLow)
          xf_printf("[%d, %d, %d] fd = %22.15E, fa = %22.15E", egrp, elem, j, fd, fa);
        if (Verbosity != xfe_VerbosityLow){
          if (fabs(fd-fa) > 1000.0*ep*ep)
            xf_printf("<-----***\n");	  
          else if (fabs(fd-fa) > 100.0*ep*ep)
            xf_printf("<-----**\n");
          else if (fabs(fd-fa) > 10.0*ep*ep)
            xf_printf("<-----*\n");
          else
            xf_printf("\n");
        }
        if (fabs(fd-fa) > 10.0*ep*ep) return xf_Error(xf_PING_FAILED);
        
        // decrement U
        U->GenArray[egrp].rValue[elem][j] -= ep;
      } // j
      
    } // elem
  } // egrp
  
  ierr = xf_Error(xf_DestroyVector(Value_U, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(Value_U0, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
  
}


/******************************************************************/
//   FUNCTION Definition: xf_PingOutputGCL
static int
xf_PingOutputGCL(xf_All *All, const char *Name, xf_Vector *U, xf_Vector *GCL)
{
  int ierr, j, egrp, elem;
  int sr, r;
  enum xfe_Verbosity Verbosity;
  real ep, fd, fa, s0, *s_U0, s, *s_U;
  xf_Mesh *Mesh;
  xf_Vector *J_GCL, *J_GCL0, *Value_U;
  xf_GenArray *ga;
  
  ep = 1.0e-5;
  
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;
  
  // Determine verbosity for printing purposes
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", xfe_VerbosityName, 
                                     (int ) xfe_VerbosityLast, (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  // locate vectors
  ierr = xf_Error(xf_FindMeshMotionGCLLinearization(All, Name, -1, &J_GCL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, J_GCL, "J_GCL0", xfe_False, xfe_False, 
                                       NULL, &J_GCL0, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Value_U", xfe_False, xfe_False, 
                                       NULL, &Value_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  // calculate s0 and s_U0
  ierr = xf_Error(xf_CalculateOutput(All, Name, U, &s0, Value_U, xfe_Set));
  if (ierr != xf_OK) return ierr;
  
  // copy over J_GCL -> J_GCL0
  ierr = xf_Error(xf_SetVector(J_GCL, xfe_Set, J_GCL0));
  if (ierr != xf_OK) return ierr;
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      ga = GCL->GenArray+egrp;
      r = ((ga->vr==NULL) ? ga->r : ga->vr[elem]);
      for (j=0; j<r; j++){
        // increment GCL
        GCL->GenArray[egrp].rValue[elem][j] += ep;
        
        // calculate s and s_U
        ierr = xf_Error(xf_CalculateOutput(All, Name, U, &s, Value_U, xfe_Set));
        if (ierr != xf_OK) return ierr;
        
        // check Delta s vs s_U
        fd = (s - s0)/ep;
        fa = 0.5*(J_GCL->GenArray[egrp].rValue[elem][j] +
                  J_GCL0->GenArray[egrp].rValue[elem][j]);
        if (Verbosity != xfe_VerbosityLow)
          xf_printf("[%d, %d, %d] fd = %22.15E, fa = %22.15E, fd-fa = %22.15E", 
                    egrp, elem, j, fd, fa, fabs(fd-fa));
        
        ierr = xf_Error(xf_PingCompareFDtoAnalytical(fd, fa, ep, Verbosity));
        if (ierr != xf_OK) return ierr;
        
        // decrement U
        GCL->GenArray[egrp].rValue[elem][j] -= ep;
      } // j
      
    } // elem
  } // egrp
  
  ierr = xf_Error(xf_DestroyVector(Value_U, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(J_GCL0, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
  
}



/******************************************************************/
//   FUNCTION Definition: xf_ComputeBCLinkedOutputs
static int
xf_ComputeBCLinkedOutputs(xf_All *All, xf_Vector *U, xf_KeyValue *pOutputsForBCs)
{
  int ierr, i, ibfgrp;
  int iOutput, nOutput;
  enum xfe_Bool found;
  char *OutputLinkage = NULL;
  char *OutputName = NULL;
  char **OutputNames = NULL;
  real Value;
  xf_BC *BC;
  
  // look for output-dependent BCs in All->EqnSet: add output names/values to OutputsForBCs
  for (ibfgrp=0; ibfgrp<All->Mesh->nBFaceGroup; ibfgrp++){
    // loop over mesh BCs because these do not change in possible sort of eqnset bcs
    found = xfe_False;
    for (i=0; i<All->EqnSet->BCs->nBC; i++)
      if (strcmp(All->Mesh->BFaceGroup[ibfgrp].Title, All->EqnSet->BCs->BC[i].BFGTitle) == 0){
        found = xfe_True;
        break;
      }
    if (!found) return xf_Error(xf_NOT_FOUND);
    BC = All->EqnSet->BCs->BC+i;
    if ((OutputLinkage=BC->OutputLinkage) != NULL){
      // account for there possibly being more than one output in line
      ierr = xf_Error(xf_ScanXStringAlloc(OutputLinkage, xf_MAXSTRLEN, &nOutput, &OutputNames));
      if (ierr != xf_OK) return ierr;
      
      for (iOutput=0; iOutput<nOutput; iOutput++){
        OutputName = OutputNames[iOutput];
        // add this output name to the list, if it does not exist
        ierr = xf_Error(xf_AddKeyValue(pOutputsForBCs, OutputName, "0.0", xfe_False));
        if (ierr != xf_OK) return ierr;
        // evaluate this output
        ierr = xf_Error(xf_CalculateOutput(All, OutputName, U, &Value, NULL, xfe_Set));
        if (ierr != xf_OK) return ierr;	
        // store the output in the given key-value structure
        ierr = xf_Error(xf_SetKeyValueReal(*pOutputsForBCs, OutputName, Value));
        if (ierr != xf_OK) return ierr;
      } // iOutput
      
      xf_Release2( (void **) OutputNames);
    }
  } // ibfgrp
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ComputeBCSensitivity
static int
xf_ComputeSensitivity(xf_All *All, xf_Vector *U, int nPsi, xf_Vector **Psi,
                      real *epsilon, real *Sensitivity)
{
  int ierr, i;
  xf_Vector *R;
  xf_SolverData *SolverData = NULL;
  
  // locate a residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, 
                                       xfe_True, NULL, &R, NULL));
  if (ierr != xf_OK) return ierr;  
  // create/allocate SolverData on All
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;
  // Calculate residual on All (nonzero because we adjusted/perturbed the BCs)
  ierr = xf_Error(xf_CalculateResidual(All, U, R, NULL, SolverData));
  if (ierr != xf_OK) return ierr;
  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;
  
  // loop over the adjoints and compute the error estimates
  for (i=0; i<nPsi; i++){
    // dJ = Psi^T * Rh
    ierr = xf_Error(xf_VectorDot(R, Psi[i], &Sensitivity[i]));
    if (ierr != xf_OK) return ierr;
    //here we assume only one parameter so we have only one epsilon
    //epsilon would be NULL for error estimation
    if (epsilon != NULL)
      Sensitivity[i] /= (*epsilon);//Note the sign convention in Psi
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetOutputDependentBCs
int
xf_SetOutputDependentBCs(xf_All *All, xf_Vector *U, xf_KeyValue *pOutputsForBCs, 
                         enum xfe_Bool PostSolve, enum xfe_Bool *Converged, 
                         enum xfe_Bool *Trim)
{
  int ierr, i, k;
  real epsilon, Sensitivity;
  enum xfe_Bool found = xfe_False;
  xf_Vector **Adj;
  clock_t clock_Start, clock_End;
  
  if (pOutputsForBCs->nKey > 1)
    return xf_Error(xf_NOT_SUPPORTED);
  //do nothing if user halted
  if (xf_CheckUserHalt(NULL)) return xf_OK;
  
  // compute outputs if there are no outputs given
  ierr = xf_Error(xf_ComputeBCLinkedOutputs(All, U, pOutputsForBCs));
  if (ierr != xf_OK) return ierr;
  
  (*Converged) = xfe_True;
  (*Trim) = xfe_False;
  // nothing to do if no bc-linked outputs
  if (pOutputsForBCs->nKey == 0) return xf_OK;
  
  if (PostSolve) {
    //Calculate epsilon and perturb bc
    ierr = xf_Error(xf_EqnSetOutputDependentBCs(All->EqnSet, pOutputsForBCs,
                                                NULL, &epsilon, Converged));
    if (ierr != xf_OK) return ierr;
    
    if (epsilon > 0.){//if epsilon < 0, bc is not of trim type
      //Compute sensitivity
      //dJ/dParam ~= Psi^T * Rh/(epsilon)
      for (k = 0; k < pOutputsForBCs->nKey; k++){
        ierr = xf_Error(xf_FindAdjointVectors(All, U, pOutputsForBCs->Key[k], 
                                              xfe_False, xfe_False, &i, &Adj, &found));
        if (ierr != xf_OK) return ierr;
        if (!(*Converged)){
          xf_printf("\nTrimming BC:\nCalling steady adjoint solver for BC-linked %s output.\n",
                    pOutputsForBCs->Key[k]);
          clock_Start = clock(); // start timer
          if (found)
            xf_printf("%s adjoint found. Re-using it.\n",pOutputsForBCs->Key[k]);
          
          ierr = xf_Error(xf_SolveAdjoints(All, 0.0, 1.0, xfe_False, U, i, 
                                           NULL, Adj, NULL,xfe_False, !found));
          if (ierr != xf_OK) return ierr;
          
          clock_End = clock(); // end timer
          xf_printf("\nSteady adjoint solve CPU time = %.10E for BC-linked %s output.\n", 
                    ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC,
                    pOutputsForBCs->Key[k]);
          
          ierr = xf_Error(xf_ComputeSensitivity(All, U, 1, Adj, 
                                                &epsilon, &Sensitivity));
          if (ierr != xf_OK) return ierr;
          xf_printf("Sensitivity BC: %1.12e\n",Sensitivity);
        }
        //Correct bc and exit to re-solve.
        ierr = xf_Error(xf_EqnSetOutputDependentBCs(All->EqnSet, pOutputsForBCs,
                                                    &Sensitivity, &epsilon, Converged));
        if (ierr != xf_OK) return ierr;
        
        //just destroy the pointer
        xf_Release( (void *) Adj);
        Adj = NULL;
      }
      (*Trim) = xfe_True;
    }
    else {
      (*Trim) = xfe_False;
    }
    
  }
  else {
    // set output-dependent BCs via an equation-set-specific function
    ierr = xf_Error(xf_EqnSetOutputDependentBCs(All->EqnSet, pOutputsForBCs,
                                                NULL, NULL, Converged));
    if (ierr != xf_OK) return ierr;
  }
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CompOutputSensitivity
int
xf_CompOutputSensitivities(xf_All *All, xf_Vector *U, xf_Output *Output)
{
  int ierr, iSens, nPsi;
  real epsilon;
  xf_Vector **Psi;
  enum xfe_Bool found;
  
  //do nothing if user halted
  if (xf_CheckUserHalt(NULL)) return xf_OK;
  
  ierr = xf_Error(xf_CalculateOutput(All, Output->Name, U, 
                                     &Output->Value, NULL, xfe_Set));
  if (ierr != xf_OK) return ierr;
      
  for (iSens = 0; iSens < Output->nSensitivity; iSens++) {
    //Find adjoint vector for Output
    found = xfe_False;
    ierr = xf_Error(xf_FindAdjointVectors(All, U, Output->Name, xfe_False, 
                                          xfe_False, &nPsi, &Psi, &found));
    if (ierr != xf_OK) return ierr;
    if (nPsi != 1) return xf_Error(xf_MULTIPLE_MATCHES);
    
    //Solve adjoint equation
    ierr = xf_Error(xf_SolveAdjoints(All, 0.0, 1.0, xfe_False, U, nPsi, 
                                     NULL, Psi, NULL, xfe_False, !found));
    if (ierr != xf_OK) return ierr;

    switch (Output->Sensitivity[iSens].Type) {
      case xfe_SensitivityEqnset:
        //Perturb eqnset parameter
        ierr = xf_Error(xf_EqnSetPerturbParam(All->EqnSet, 
                                              Output->Sensitivity+iSens, 
                                              &epsilon, xfe_False));
        if (ierr != xf_OK) return ierr;
        if (Output->Sensitivity[iSens].nVar < 0){//not yet initialized
          Output->Sensitivity[iSens].nVar = 1;
          ierr = xf_Error(xf_Alloc((void **)&Output->Sensitivity[iSens].value, 
                                   Output->Sensitivity[iSens].nVar, sizeof(real)));
          if (ierr != xf_OK) return ierr;
        }
        //calculate residual
        //dJ/dPar = Psi*R^T/dPar
        ierr = xf_Error(xf_ComputeSensitivity(All, U, nPsi, Psi, &epsilon, 
                                              &Output->Sensitivity[iSens].value[0]));
        if (ierr != xf_OK) return ierr;
        xf_printf("%s sensitivity w.r.t. %s: %1.12e\n",Output->Name,
                  Output->Sensitivity[iSens].ParamName,
                  Output->Sensitivity[iSens].value[0]);
        //return eqnset parameter back to its original value
        epsilon *= -1.0; //reverse sign and perturb parameter back
        ierr = xf_Error(xf_EqnSetPerturbParam(All->EqnSet, 
                                              Output->Sensitivity+iSens, 
                                              &epsilon, xfe_True));
        if (ierr != xf_OK) return ierr;
        break;
      case xfe_SensitivityMesh:
      case xfe_SensitivitySolver:
        return xf_Error(xf_NOT_SUPPORTED);
      default:
        xf_printf("Unknown sensitivity type.\n");
        return xf_Error(xf_INPUT_ERROR);
        break;
    }
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ErrEstBC
int
xf_ErrEstBC(xf_All *All, xf_Vector *U, int nPsi, xf_Vector **Psi,
            xf_KeyValue *pOutputsForBCs)
{
  int ierr, i;
  real *OutputError = NULL;
  enum xfe_Bool Converged, Trim;
  xf_KeyValue NewOutputsForBCs;
  
  // nothing to do if no bc-linked outputs
  if (pOutputsForBCs->nKey == 0) return xf_OK;
  
  // initialize key-value structure for output-dependent BCs
  ierr = xf_Error(xf_InitKeyValue(&NewOutputsForBCs));
  if (ierr != xf_OK) return ierr;
  
  // set BCs using current state outputs
  ierr = xf_Error(xf_SetOutputDependentBCs(All, U, &NewOutputsForBCs, 
                                           xfe_True, &Converged, &Trim));
  if (ierr != xf_OK) return ierr;
  if (!Trim){//BC is not a trim-type
    ierr = xf_Error(xf_Alloc((void **)&OutputError, nPsi, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // *** Compute errors via adjoint-weighted residuals ***
    // dJ = Psi^T * Rh
    ierr = xf_Error(xf_ComputeSensitivity(All, U, nPsi, Psi, NULL, 
                                          OutputError));
    if (ierr != xf_OK) return ierr;
    
    // loop over the adjoints and print the error estimates(sensitivities)
    for (i=0; i<nPsi; i++)
      xf_printf("BC Error estimation: AWR says %s should be incremented by %.15f\n", 
                Psi[i]->OutputName, OutputError[i]);
    xf_Release((void *)OutputError);
    
    // (re)set BCs using old state outputs
    ierr = xf_Error(xf_SetOutputDependentBCs(All, U, pOutputsForBCs, 
                                             xfe_True, &Converged, &Trim));
    if (ierr != xf_OK) return ierr;
    
    // store new outputs in OutputsForBCs
    ierr = xf_Error(xf_CopyKeyValue(&NewOutputsForBCs, pOutputsForBCs));
    if (ierr != xf_OK) return ierr;
  }
  // destroy key-value structure for output-dependent BCs
  ierr = xf_Error(xf_DestroyKeyValue(&NewOutputsForBCs));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
  
}



#if( UNIT_TEST==1 )
#include "xf_Output.test.in"
#endif

