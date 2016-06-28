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
 FILE:  xf_Adapt.c
 
 This file contains top-level functions for adaptation.
 
 */

#include "xf.h"
#include "xf_AllStruct.h"
#include "xf_All.h"
#include "xf_Memory.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Data.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_MeshTools.h"
#include "xf_Math.h"
#include "xf_AdaptStruct.h"
#include "xf_AdaptHang.h"
#include "xf_ErrEst.h"
#include "xf_MPI.h"
#include "xf_Mesh.h"
#include "xf_DataMath.h"
#include "xf_Quad.h"
#include "xf_EqnSetHook.h"
#include "xf_Residual.h"
#include "xf_AllPull.h"
#include "xf_Solver.h"
#include "xf_LeanSolver.h"
#include "xf_Output.h"
#include "xf_LinearSolver.h"
#include "xf_Penalty.h"
#include "xf_EqnSet.h"
#include "xf_SolverTools.h"
#include <time.h>

//Yu's own source code
#include "./xfYu_Adaption.c"

/******************************************************************/
//   FUNCTION Definition: xf_GetNumRefOpt
int 
xf_GetNumRefOpt(enum xfe_ShapeType Shape, int *n)
{
  switch (Shape){
    case xfe_Segment:
      (*n) = xfe_SegRefLast;
      break;
    case xfe_Quadrilateral:
      (*n) = xfe_QuadRefLast;
      break; 
    case xfe_Triangle:
      (*n) = 5;  // considering only single edge splitting for now
      break;
    case xfe_Tetrahedron:
      (*n) = 8; // also only considering single edge splitting
      break;
    case xfe_Hexahedron: 
      (*n) = xfe_HexRefLast;
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  return xf_OK;
  
}

/******************************************************************/
//   FUNCTION Definition: xf_FindAdaptIndicator
static int 
xf_FindAdaptIndicator(xf_All *All, enum xfe_Bool WithRefOpt, 
                      xf_Vector **pAdaptIndicator)
{
  int ierr, egrp, nopt;
  int *rvec = NULL;
  enum xfe_ShapeType Shape;
  xf_Data *D;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  
  if (WithRefOpt){
    ierr = xf_Error(xf_Alloc( (void **) &rvec, Mesh->nElemGroup, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_GetNumRefOpt(Shape, &nopt));
      if (ierr != xf_OK) return ierr;
      
      rvec[egrp] = nopt; // includes no refinement as an option
    } // egrp
  }
  
  ierr = xf_Error(xf_FindVector(All, "AdaptIndicator", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                NULL, NULL, NULL, NULL, rvec, xfe_SizeReal, xfe_False,  xfe_True, &D, 
                                pAdaptIndicator, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True; // make indicator writeable
  
  xf_Release( (void *) rvec);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AdaptIndicatorScalar
static int 
xf_AdaptIndicatorScalar(xf_All *All, char *ScalarName, 
                        xf_Vector *AdaptIndicator)
{
  int ierr, sr, iq, nq, pnq, dim;
  int Order, QuadOrder, nn;
  int egrp, elem;
  int *IParam;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged;
  real *RParam, *xq, *u, *s;
  real *EU, LocValue;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData, *GeomPhiData;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  xf_Data *D;
  xf_Vector *U;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  
  EqnSet = All->EqnSet;
  sr     = EqnSet->StateRank;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // zero out adaptive indicator
  ierr = xf_Error(xf_SetZeroVector(AdaptIndicator));
  if (ierr != xf_OK) return ierr;
  
  // pull off primal state
  ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D, NULL));
  if (ierr != xf_OK) return ierr;
  U = (xf_Vector *) D->Data;
  
  // initialize vars to NULL
  QuadData    = NULL;
  PhiData     = NULL;
  JData       = NULL;
  u           = NULL;
  s           = NULL;
  GeomPhiData = NULL;
  
  pnq         = -1;  // previous number of quad points
  
  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    // Determine Basis and Order from the state, U
    Basis = U->Basis[egrp];
    
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      Order = xf_InterpOrder(U, egrp, elem);
      
      // determine required integration order
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, EqnSet, egrp, Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      /* Pull off quad points for the element; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      
      // compute basis functions (and grads) if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
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
        ierr = xf_Error(xf_ReAlloc( (void **)  &s, nq, sizeof(real)));
        if (ierr != xf_OK) return ierr;
      }
      
      
      EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
      
      // interpolate state at quad points
      xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u);      
      
      // call eqnset specific function for scalar
      ierr = xf_Error(xf_EqnSetScalar(EqnSet, ScalarName, IParam, 
                                      RParam, nq, u, NULL, s, NULL, NULL, NULL, 0.0));
      if (ierr != xf_OK) return ierr;
      
      // modify scalar (absolute value, quad weights, Jacobian data)
      for (iq=0; iq<nq; iq++) s[iq] = fabs(s[iq])*QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      // sum scalar over quad points, add to LocValue
      for (iq=0, LocValue = 0.; iq<nq; iq++) LocValue += s[iq];
      
      AdaptIndicator->GenArray[egrp].rValue[elem][0] = LocValue;
      
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
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void  *) IParam);
  xf_Release( (void  *) RParam);
  xf_Release( (void *) u);
  xf_Release( (void *) s);
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AdaptIndicatorTest
static int 
xf_AdaptIndicatorTest(xf_All *All, xf_Vector *AdaptIndicator)
{
  // A hard-coded adaptive indicator for testing
  int ierr, egrp, elem, k;
  int ibfgrp, ibface;
  real err;
  xf_Mesh *Mesh;
  xf_BFace BFace;
  
  Mesh = All->Mesh;
  err = 0.0;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      for (k=0; k<AdaptIndicator->GenArray[egrp].r; k++)
        AdaptIndicator->GenArray[egrp].rValue[elem][k] = 0.0; //err;
      //AdaptIndicator->GenArray[egrp].rValue[elem][0] *= 2.0;
      err += .001;
    } // elem
  } // egrp
  
  ibfgrp = 0;
  for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
    BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
    egrp = BFace.ElemGroup;
    elem = BFace.Elem;
    AdaptIndicator->GenArray[egrp].rValue[elem][0] = 1.0;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AdaptSmoothRef
static int 
xf_AdaptSmoothRef(xf_All *All, xf_Vector *RefIndicator)
{
  /*
   PURPOSE:
   
   Smoothes adaptation to eliminate islands, fill voids, etc.
   
   INPUTS:
   
   All : All structure
   RefIndicator : coming in, elements to be refined have this set to 1;
   all other elements have this set to 0
   
   OUTPUTS: 
   
   RefIndicator : modified refinement indicator, based on smoothing operations
   
   RETURN:
   
   Error Code
   */
  int ierr, egrp, elem, face;
  int count, nneigh, nface, negrp;
  int egN, eN, faceN, iiface, hang;
  int *RI;
  real frac;
  const real frac_island = 0.7;
  const real frac_void   = 0.7;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  
  xf_printf("Smoothing refinement indicator.\n");
  
  negrp = Mesh->nElemGroup;
  
  // remove islands
  for (egrp=0; egrp<negrp; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      RI = RefIndicator->GenArray[egrp].iValue[elem];
      if (RI[0] == 0) continue;
      nface = Mesh->ElemGroup[egrp].nFace[elem];
      count  = 0;
      nneigh = 0;
      for (face=0; face<nface; face++){
        ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, &egN, &eN, &faceN));
        if (ierr != xf_OK) return ierr;
        if (egN >= 0){
          nneigh++;
          // increment count if face is not hanging and other side is not refined
          iiface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
          if ((Mesh->IFace[iiface].HangNumber            == 0) &&
              (RefIndicator->GenArray[egN].iValue[eN][0] == 0)) count++;
        }
      }
      frac = ((real) count) / ((real) nneigh);
      if (frac > frac_island) RI[0] = 0;
    }
  
  
  // fill voids
  for (egrp=0; egrp<negrp; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      RI = RefIndicator->GenArray[egrp].iValue[elem];
      if (RI[0] > 0) continue;
      nface = Mesh->ElemGroup[egrp].nFace[elem];
      count  = 0;
      nneigh = 0;
      for (face=0; face<nface; face++){
        ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, &egN, &eN, &faceN));
        if (ierr != xf_OK) return ierr;
        if (egN >= 0){
          nneigh++;
          // check if on coarse side of a hanging face
          ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, NULL, NULL, NULL));
          if (ierr != xf_OK) return ierr;
          if ((hang) || (RefIndicator->GenArray[egN].iValue[eN][0]==1)) count++;
        }
      }
      if (nneigh == 0) continue;
      frac = ((real) count) / ((real) nneigh);
      if (frac > frac_void) RI[0] = 1;
    }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AdaptOrderRef
static int 
xf_AdaptOrderRef(xf_All *All, const char *SavePrefix, xf_Vector *RefIndicator, 
                 xf_TimeHistData *TimeHistData, xf_TimeHistData *OldTimeHistData,
                 xf_Vector *BasisOrderElem)
{
/*

PURPOSE:
   
  Adapts spatial mesh in order, p.  Called in steady state, unsteady
  static-mesh, and unsteady dynamic-mesh cases.  Time histories and
  basis/order info in time are needed for unsteady cases.
      
INPUTS:
   
   All : All structure
   SavePrefix : prefix to use for saving files (e.g. order information for unsteady)
   RefIndicator : spatial refinement/coarsening indicator (for all time slabs if dynamic)
   TimeHistData : new (adapted) time history
   OldTimeHistData : old (un-adapted) time history
   BasisOrderElem : basis/order information (for all time slabs if dynamic)
   
OUTPUTS: 
   
   None: order is changed in all existing vectors and order vectors
   are written for use in next adaptive iteration for dynamic
   refinement cases.
   
RETURN:
   
   Error Code
 */

  // order refinement
  int ierr, i, j;
  int Order, ref;
  int iTime, nTime;
  int iTimeNew, nTimeNew;
  int AdaptOrderMin, AdaptOrderMax;
  int NeedRefine = 0;
  enum xfe_Bool Interpolated, ParallelFlag, WriteVOrder;
  enum xfe_Bool StaticMesh, DynamicMesh;
  enum xfe_Bool HitMin, HitMax;
  enum xfe_Verbosity Verbosity;
  real TargetTime;
  char OutputFile[xf_MAXSTRLEN];
  xf_Vector *V, *VOrder;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Mesh *Mesh;

  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  // verify if we are writing the VOrder vector. 
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "AdaptWriteVOrder", 
                                     &WriteVOrder));
  if (ierr != xf_OK) return ierr;
  
  // minimum order for p-adaptation
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "AdaptOrderMin", &AdaptOrderMin));
  if (ierr != xf_OK) return ierr;

  // maximum order for p-adaptation
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "AdaptOrderMax", &AdaptOrderMax));
  if (ierr != xf_OK) return ierr;

  // locate a vector for specifying desired order
  ierr = xf_Error(xf_FindVector(All, "VOrder", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_False,
                                xfe_False, NULL, &VOrder, NULL));
  if (ierr != xf_OK) return ierr;
  
  Mesh = All->Mesh;
  DataSet = All->DataSet;
  D = DataSet->Head;
  
  /*when doing hp-adaptation, each processor will 
    independently call this function with serialized meshes 
    and this function should have a serial behavior*/
  if (Mesh->ParallelInfo != NULL) ParallelFlag = xfe_True;
  else ParallelFlag = xfe_False;
  
  if (ParallelFlag){
    // make sure spatial ref info gets to halos
    // begin communication of halo data
    ierr = xf_Error(xf_HaloExchangeVectorBegin(RefIndicator));
    if (ierr != xf_OK) return ierr;
    // end communication of halo data
    ierr = xf_Error(xf_HaloExchangeVectorEnd(RefIndicator));
    if (ierr != xf_OK) return ierr;
  }

  // is this a dynamic mesh refinement?
  DynamicMesh = ((RefIndicator->StateRank != 1) && (nTime > 1));
  
  // do not project IC vectors for dynamic mesh refinement
  if (!DynamicMesh){
    // project existing vectors based on RefIndicator->GenArray[*].iValue[*][0]
    while (D != NULL){
    
      if (D->Type == xfe_Vector){
        V = (xf_Vector *) D->Data;
        Interpolated = ((V->Basis != NULL) && (V->Order != NULL));
      
        /** Only adapt interpolated state and adjoint data **/
        if ((V->Linkage == xfe_LinkageGlobElem) && (Interpolated) &&
            ((V->SolverRole == xfe_SolverRolePrimalState) ||
             (V->SolverRole == xfe_SolverRoleAdjointState))){
            
          // only adapt real vectors with data
          if ((V->nArray <= 0) || (V->GenArray[0].Size != xfe_SizeReal)) continue;
            
          NeedRefine = 0;
            
          // set VOrder
          for (i=0; i<VOrder->nArray; i++){
            for (j=0; j<VOrder->GenArray[i].n; j++){
              // current order
              Order = xf_InterpOrder(V, i, j);
              // desired order (refindicator is treated like an order increment)
              VOrder->GenArray[i].iValue[j][0] = max(0,Order + RefIndicator->GenArray[i].iValue[j][0]);
              if (RefIndicator->GenArray[i].iValue[j][0] != 0){
                NeedRefine = 1;
                //xf_printf("Setting order of (%d, %d) = %d\n", i,j, VOrder->GenArray[i].iValue[j][0]);
              }
            }
          } // i
          if (ParallelFlag){
            // all-reduce NeedRefine
            ierr = xf_Error(xf_MPI_Allreduce(&NeedRefine, 1, xfe_SizeInt, xfe_MPI_MAX));
            if (ierr != xf_OK) return ierr;
          }
          if (NeedRefine == 1){
              
            // project V in place
            if (Verbosity != xfe_VerbosityLow)
              xf_printf("Projecting vector %s: to p-refined orders\n", D->Title);

            ierr = xf_Error(xf_ProjectVectorInPlace_VOrder(All->Mesh, All->DataSet, V, V->Basis, 
                                                           xfe_BasisLast, VOrder));
            if (ierr != xf_OK) return ierr;
              
            // if steady, write out VOrder file, and set key value for later use
            if ((RefIndicator->StateRank == 1) && (WriteVOrder)){
              sprintf(OutputFile, "%s_VOrder0.data\0", SavePrefix);
              ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "VOrder", VOrder, OutputFile));
              if (ierr != xf_OK) return ierr;
              ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "VOrderFile", OutputFile));
              if (ierr != xf_OK) return ierr;
            }
          }
            
        }
      }
      D = D->Next;
    } // while D != NULL
  }  
  
  // we might be doing unsteady ... (even static meshes execute this)
  
  if ((BasisOrderElem != NULL) && (TimeHistData != NULL)) { // yes we are doing unsteady
    
    // consistency checks
    if (OldTimeHistData == NULL) return xf_Error(xf_INPUT_ERROR);
    nTime = OldTimeHistData->nTime; // # old time slabs

    // is this a static mesh refinement?
    StaticMesh = ((RefIndicator->StateRank == 1) && (nTime > 1));

    // if dynamic, we need time steps to match
    if ((!StaticMesh) && (RefIndicator->StateRank != nTime)) 
      return xf_Error(xf_INPUT_ERROR);
    
    // # new time slabs
    nTimeNew = TimeHistData->nTime;

    // starting index of old time slab
    iTime = 0;

    // special case for static mesh: just write one order file
    if (StaticMesh){
      nTime = nTimeNew = 1;
    }

    if (BasisOrderElem->StateRank != 2*nTime) return xf_Error(xf_INPUT_ERROR);

    // begin loop over new time slabs
    HitMax = HitMin = xfe_False;
    for (iTimeNew=0; iTimeNew<nTimeNew; iTimeNew++){

      // Determine index of old time slab based on time histories
      if (!StaticMesh){
        // use midpoint of current time slab as target
        TargetTime = TimeHistData->Time[iTimeNew] + 0.5*TimeHistData->TimeStep[iTimeNew];
        
        // Bring old time history up to target time
        while ( ((OldTimeHistData->Time[iTime] + OldTimeHistData->TimeStep[iTime]) < TargetTime)
                && (iTime < nTime))
          iTime++;
        
        // should never get past the end of the old time history (mismatch of times)
        if (iTime >= nTime) return xf_Error(xf_OUT_OF_BOUNDS);
      }

      // set orders (in VOrder) via old time slab index, iTime, and BasisOrderElem
      for (i=0; i<VOrder->nArray; i++){
        for (j=0; j<VOrder->GenArray[i].n; j++){
          // current Order at target time
          Order = BasisOrderElem->GenArray[i].iValue[j][1*nTime+iTime];
          ref = RefIndicator->GenArray[i].iValue[j][iTime];
          // new order: coarsened if ref<0, refined if ref>0
          Order = Order+ref;
          if (Order < AdaptOrderMin) HitMin = xfe_True;
          Order = max(AdaptOrderMin, Order); // bound from below
          if (Order > AdaptOrderMax)  HitMax = xfe_True;
          Order = min(AdaptOrderMax, Order); // bound from above
          VOrder->GenArray[i].iValue[j][0] = Order;
          /* 	  xf_printf("SavePrefixNext = %s, iTime = %d, iTimeNew = %d, i=%d, j=%d, Order=%d\n", */
          /* 		    SavePrefix, iTime, iTimeNew, i, j, Order); */
        }
      } // i
      // write out VOrder for every new time slab (can consolidate later)
      sprintf(OutputFile, "%s_VOrder%d.data\0", SavePrefix, iTimeNew); // use new time index
      /* 	  xf_printf("printing: %s\n", OutputFile); */
      /* 	  for (i=0; i<VOrder->nArray; i++){ */
      /* 	    for (j=0; j<VOrder->GenArray[i].n; j++) */
      /* 	      xf_printf("VOrder(%d,%d) = %d\n",  i, j, VOrder->GenArray[i].iValue[j][0]); */
      /* 	  } */
      ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "VOrder", VOrder, OutputFile));
      if (ierr != xf_OK) return ierr;      

    } // iTimeNew
      
    if (HitMin) xf_printf("Hit minimum order limit of %d at least once.\n", AdaptOrderMin);
    if (HitMax) xf_printf("Hit maximum order limit of %d at least once.\n", AdaptOrderMax);
    
    // set key-value that gives VOrderFile for IC (assume same as first(zeroth) time slab)
    // NOTE: leaving IC alone until we have a better system for adapting the IC order
    /* sprintf(OutputFile, "%s_VOrder%d.data\0", SavePrefix, 0); */
    /* ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "VOrderFile", OutputFile)); */
    /* if (ierr != xf_OK) return ierr; */
    
  } // end if unsteady
  
  
  // destroy VOrder
  ierr = xf_Error(xf_DestroyVector(VOrder, xfe_True));
  if (ierr != xf_OK) return ierr;
    
  return xf_OK;  
}

/******************************************************************/
//   FUNCTION Definition: xf_OptimizeHangAnisotropic
static int 
xf_OptimizeHangAnisotropic(xf_All *All, xf_Vector *AdaptIndicator, 
                           xf_Vector *RefIndicator, int **OrderToRefine, 
                           int nelemref, enum xfe_Bool AdaptIsotropic)
{
  /*
   PURPOSE:
   
   Determines whether anisotropic hanging node refinement may be more
   efficient.  Stores optimal refinement choices in RefIndicator
   
   INPUTS:
   
   All : All structure
   AdaptIndicator: error indicator for each element of All (real numbers)
   RefIndicator : coming in, elements to be refined have this set to 1;
   all other elements have this set to 0
   AdaptIsotropic: if true, cells are marked fo isotropic refinement
   
   OUTPUTS: 
   
   RefIndicator : on return, this integer vector will be set to the
   optimal refinement option for that element (according
   to GetNumRefOpt and if p-ref is an option).
   
   RETURN:
   
   Error Code
   */
  int ierr, *nAddElem, i, j, k, OptimumRef, nLimitedP, Pmax;
  int nAddElemTot, myRank, nProc, iref, nRefineTot, OrderIncrement;
  int nref_proc, index, nAdd, nPsi, nn, sr, n, s, nz_old, nz_new, Order, OrderL, OrderR;
  int RefOpt, egrp, elem, *nRefOpt, *nElem, nElemTot, egL, eL, egR, eR;
  int **RefOptStat, nNonPhysical, *nHRefOpt, nPRefOpt, dof_new, dof_old;
  char AdaptOutput[xf_MAXSTRLEN], AdaptVariableSet[xf_MAXSTRLEN];
  char statsfname[xf_MAXSTRLEN], SavePrefix[xf_MAXSTRLEN];
  real *CostRefOpt, Rnorm, Psinorm, projection;
  real temp, growth, cost, benefit;
  enum xfe_ShapeType Shape;
  enum xfe_Bool found, Problem, LimitedP;
  enum xfe_Bool ReachedGrowth;
  enum xfe_Bool ElemPicking, ConsiderPref;
  enum xfe_AdaptOnType AdaptOn;
  enum xfe_AdaptCostType CostMetric;
  xf_Mesh *Mesh, *Mesh_Small;
  xf_All *All_Small;
  xf_Vector *Tracer, *RefIndSmall;
  xf_Vector *USmall, *RSmall, *Psi_H, **Psi;
  xf_Vector *AdaptIndicatorSmall;
  xf_Data *DSmall, *D;
  xf_SolverData *SolverData;
  FILE *statsfile;
  //temp
  clock_t start, end;
  double elapsed;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptOn", 
                                     xfe_AdaptOnName, (int ) xfe_AdaptOnLast, 
                                     (int *) &AdaptOn));
  if (ierr != xf_OK) return ierr;
  //return immediately if uniform refinement
  if (AdaptOn == xfe_AdaptOnUniform)
    return xf_OK;
  
  //count the time on processor 0
  if (myRank == 0)
    start = clock();
  
  //Broadcast mesh
  ierr = xf_Error(xf_BcastMesh(&All->Mesh));
  if (ierr != xf_OK) return ierr;
  
  //Broadcast dataset
  ierr = xf_Error(xf_BcastDataSet(&All->DataSet));
  if (ierr != xf_OK) return ierr;
  
  Mesh = All->Mesh;
  
  //reset pointers locally
  if (nProc > 1){
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, 
                                       "AdaptIndicator", 
                                       xfe_Vector, &D));
    if (ierr != xf_OK) return ierr;
    AdaptIndicator = (xf_Vector *) D->Data;
    
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, 
                                       "RefIndicator", 
                                       xfe_Vector, &D));
    if (ierr != xf_OK) return ierr;
    RefIndicator = (xf_Vector *) D->Data;
  }
  
  //should we consider P-refinement
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                     "AdaptIncludeP", 
                                     &ConsiderPref));
  if (ierr != xf_OK) return ierr;
  if (ConsiderPref)
    nPRefOpt = 1;
  else
    nPRefOpt = 0;
  OrderIncrement = 1;
  
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, 
                                    "AdaptPmax", &Pmax));
  if (ierr != xf_OK) return ierr;
  
  // determine the element flagging method
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                     "AdaptFixedGrowth", 
                                     &ElemPicking));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, 
                                     "AdaptFixedGrowthFactor", &growth));
  if (ierr != xf_OK) return ierr;
  
  ReachedGrowth = xfe_False;
  //Cost metric
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptCostMetric", 
                                     xfe_AdaptCostName, (int ) xfe_AdaptCostLast, 
                                     (int *) &CostMetric));
  if (ierr != xf_OK) return ierr;
  
  // SavePrefix
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetnElem(Mesh, &nElem, &nElemTot));
  if (ierr != xf_OK) return ierr;
  //no need to reduce because each processor has a full copy of the mesh
  
  ierr = xf_Error(xf_Alloc((void **)&nRefOpt, Mesh->nElemGroup, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **)&nHRefOpt, Mesh->nElemGroup, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++)
    nRefOpt[egrp] = 0;
  
  nNonPhysical = 0;
  nLimitedP = 0;
  nAddElemTot = 0;
  sr = All->EqnSet->StateRank;
  nref_proc = nelemref/nProc + 1;
  
  /****************************************************************************/
  for (iref = 0; iref < nref_proc; iref++){
    if (!ReachedGrowth){
      index = iref*nProc + myRank;
      
      if (index < nelemref){
        egrp = OrderToRefine[index][0];
        elem = OrderToRefine[index][1];
        nz_old = -1;
        //getting the element shape
        ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
        if (ierr != xf_OK) return ierr;
        
        //getting number of H-refinement options
        if (AdaptIsotropic)
          nHRefOpt[egrp] = 2;//none+uniform
        else {
          ierr = xf_Error(xf_GetNumRefOpt(Shape, &nHRefOpt[egrp]));
          if (ierr != xf_OK) return ierr;
        }
        nRefOpt[egrp] = nHRefOpt[egrp]+nPRefOpt;
        
        //allocating array to store the cost of each refinement option
        ierr = xf_Error(xf_Alloc((void **)&CostRefOpt, nRefOpt[egrp], sizeof(real)));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_Alloc((void **)&nAddElem, nRefOpt[egrp], sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        CostRefOpt[0] = AdaptIndicator->GenArray[egrp].rValue[elem][0];
        LimitedP = xfe_False;
        Problem = xfe_False;
        //evaluating all the refinement options (excluding none)
        if (!AdaptIsotropic || ConsiderPref){
          for (RefOpt = 1; RefOpt < nRefOpt[egrp]; RefOpt++){
            //pulling element
            ierr = xf_Error(xf_CreateAll(&All_Small, xfe_False));
            if (ierr != xf_OK) return ierr;
            
            ierr = xf_Error(xf_AllPull(All, All_Small, egrp, elem));
            if (ierr != xf_OK) return ierr;
            
            Mesh_Small = All_Small->Mesh;
            
            ierr = xf_Error(xf_FindPrimalState(All_Small->DataSet, 
                                               0, &DSmall,
                                               NULL));
            if (ierr != xf_OK) return ierr;
            
            USmall = (xf_Vector *)DSmall->Data;
            
            //count nz
            if (nz_old < 0){
              nz_old = 0;
              Order = xf_InterpOrder(USmall, egrp, 0);//central element
              nz_old += pow(Order+1,2*Mesh_Small->Dim);
              
              if (USmall->GenArray[egrp].vr == NULL)
                dof_old = USmall->GenArray[egrp].r/sr;
              else 
                dof_old = USmall->GenArray[egrp].vr[0]/sr;
              
              for (k = 0; k < Mesh_Small->ElemGroup[egrp].nFace[0]; k++) {
                if (Mesh_Small->ElemGroup[egrp].Face[0][k].Group == xf_INTERIORFACE) {
                  n = Mesh_Small->ElemGroup[egrp].Face[0][k].Number;
                  egL = Mesh_Small->IFace[n].ElemGroupL;
                  eL = Mesh_Small->IFace[n].ElemL;
                  egR = Mesh_Small->IFace[n].ElemGroupR;
                  eR = Mesh_Small->IFace[n].ElemR;
                  OrderL = xf_InterpOrder(USmall, egL, eL);
                  OrderR = xf_InterpOrder(USmall, egR, eR);
                  nz_old += pow(OrderL+1,Mesh_Small->Dim)*pow(OrderR+1,Mesh_Small->Dim);
                }
              }//k
            }
            
            //Tracing central element
            ierr = xf_Error(xf_FindVector(All_Small, "Tracer", xfe_LinkageGlobElem, 
                                          1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, 
                                          xfe_False, xfe_True, &DSmall, &Tracer, NULL));
            if (ierr != xf_OK) return ierr;
            
            //Making Tracer writable so it will be transferred to the refined mesh
            DSmall->ReadWrite = xfe_True;
            
            //Setting tracer to zero
            ierr = xf_Error(xf_SetZeroVector(Tracer));
            if (ierr != xf_OK) return ierr;
            
            //Marking central element
            Tracer->GenArray[egrp].iValue[0][0] = 1;
            
            //Refining central element
            ierr = xf_Error(xf_FindVector(All_Small, "RefIndicator", xfe_LinkageGlobElem, 
                                          1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, 
                                          xfe_False, xfe_False, NULL, &RefIndSmall, &found));
            if (ierr != xf_OK) return ierr;
            //if (!found && nProc > 1) return xf_Error(xf_CODE_LOGIC_ERROR);
            
            ierr = xf_Error(xf_SetZeroVector(RefIndSmall));
            if (ierr != xf_OK) return ierr;
            
            //turn off verbosity
            ierr = xf_Error(xf_SetKeyValue(All_Small->Param->KeyValue, "Verbosity", 
                                           "Low"));
            if (ierr != xf_OK) return ierr;
            //do not write VOrder vector file
            ierr = xf_Error(xf_SetKeyValueBool(All_Small->Param->KeyValue, "AdaptWriteVOrder", 
                                               xfe_False));
            if (ierr != xf_OK) return ierr;
            
            if (RefOpt < nHRefOpt[egrp]){
              //Central element is the 0th element of group: egrp
              RefIndSmall->GenArray[egrp].iValue[0][0] = RefOpt;
              
              ierr = xf_Error(xf_AdaptHang(All_Small, RefIndSmall));
              if (ierr != xf_OK) return ierr;
              
              ierr = xf_Error(xf_DestroyVector(RefIndSmall, xfe_True));
              if (ierr != xf_OK) return ierr;
            }
            else if (ConsiderPref) {
              if (xf_InterpOrder(USmall, egrp, 0) + OrderIncrement > Pmax) {
                LimitedP = xfe_True;
                nLimitedP++;
                break;
              }
              else {
                RefIndSmall->GenArray[egrp].iValue[0][0] = OrderIncrement;//order increment of 1
                
                ierr = xf_Error(xf_AdaptOrderRef(All_Small, SavePrefix, RefIndSmall, NULL, NULL, NULL));
                if (ierr != xf_OK) return ierr;
              }
            }
            
            Mesh_Small = All_Small->Mesh;
            
            //Updating State in the refined mesh (U_H)
            ierr = xf_Error(xf_FindPrimalState(All_Small->DataSet, 
                                               0, &DSmall,
                                               NULL));
            if (ierr != xf_OK) return ierr;
            
            USmall = (xf_Vector *)DSmall->Data;
            
            ierr = xf_Error(xf_CreateSolverData(All_Small, &SolverData));
            if (ierr != xf_OK) return ierr;
            
            //compute residual
            ierr = xf_Error(xf_FindSimilarVector(All_Small, USmall, "Residual", xfe_False, 
                                                 xfe_True, NULL, &RSmall, NULL));
            if (ierr != xf_OK) return ierr;
            
            ierr = xf_Error(xf_CalculateResidual(All_Small, USmall, RSmall, NULL, SolverData));
            //Checking for non physical error
            if (ierr != xf_OK){
              //Proceeding to the next element
              Problem = xfe_True;
              break;
            }
            else if (ierr != xf_OK) return ierr;
          
            /* Determine which output we are adapting on*/
            ierr = xf_Error(xf_GetKeyValue(All_Small->Param->KeyValue, "AdaptOutput", 
                                           AdaptOutput));
            if (ierr != xf_OK) return ierr;
            
            ierr = xf_Error(xf_FindAdjointVectors(All_Small, USmall, AdaptOutput, xfe_False,
                                                  xfe_False, &nPsi, &Psi, &found));
            if (ierr != xf_OK) return ierr;
            if (nPsi != 1) return xf_Error(xf_MULTIPLE_MATCHES);
            if (!found) return xf_Error(xf_CODE_LOGIC_ERROR);
            Psi_H = Psi[0];
            
            //Projecting the residual in the adjoint direction
            ierr = xf_Error(xf_FindAdaptIndicator(All_Small, xfe_False, &AdaptIndicatorSmall));
            if (ierr != xf_OK) return ierr;
            
            //Setting Tracer pointer to the refined mesh
            found = xfe_False;
            ierr = xf_Error(xf_FindVector(All_Small, "Tracer", xfe_LinkageGlobElem, 
                                          1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, 
                                          xfe_False, xfe_True, NULL, &Tracer, &found));
            if (ierr != xf_OK) return ierr;
            if (found != xfe_True) return xf_Error(xf_CODE_LOGIC_ERROR);
            
            sr = All_Small->EqnSet->StateRank;

            nz_new = 0;
            dof_new = 0;
            for (i = 0; i < Mesh_Small->nElemGroup; i++){
              for (j = 0; j < Mesh_Small->ElemGroup[i].nElem; j++){
                if (Tracer->GenArray[i].iValue[j][0] == 1){
                  Order = xf_InterpOrder(USmall, i, j);
                  nz_new += pow(Order+1,2*Mesh_Small->Dim);
                  for (k = 0; k < Mesh_Small->ElemGroup[i].nFace[j]; k++) {
                    if (Mesh_Small->ElemGroup[i].Face[j][k].Group == xf_INTERIORFACE) {
                      n = Mesh_Small->ElemGroup[i].Face[j][k].Number;
                      egL = Mesh_Small->IFace[n].ElemGroupL;
                      eL = Mesh_Small->IFace[n].ElemL;
                      egR = Mesh_Small->IFace[n].ElemGroupR;
                      eR = Mesh_Small->IFace[n].ElemR;
                      OrderL = xf_InterpOrder(USmall, egL, eL);
                      OrderR = xf_InterpOrder(USmall, egR, eR);
                      nz_new += pow(OrderL+1,Mesh_Small->Dim)*pow(OrderR+1,Mesh_Small->Dim);
                    }
                  }//k 
                  projection = 0.0;
                  AdaptIndicatorSmall->GenArray[i].rValue[j][0] = 0.0;
                  if (USmall->GenArray[i].vr == NULL)
                    nn = USmall->GenArray[i].r/sr;
                  else 
                    nn = USmall->GenArray[i].vr[j]/sr;
                  dof_new += nn;
                  for (n = 0; n < nn; n++) {
                    for (s = 0; s < sr; s++) {
                      projection += fabs(Psi_H->GenArray[i].rValue[j][n*sr+s]*
                                    RSmall->GenArray[i].rValue[j][n*sr+s]);
                    }//s
                  }//n
                  AdaptIndicatorSmall->GenArray[i].rValue[j][0] = projection;
                }//Tracer[elem] == 1
              }//j (elements)
            }//i (element groups)
            if (CostMetric == xfe_AdaptCostDeltaDOF)
              cost = dof_new-dof_old;
            else if (CostMetric == xfe_AdaptCostDOF)
	      cost = dof_new;
            else if (CostMetric == xfe_AdaptCostDeltaNonZeros)
	      cost = nz_new-nz_old;
            else if (CostMetric == xfe_AdaptCostNonZeros)
              cost = nz_new;
            else
              return xf_Error(xf_INPUT_ERROR);
            
            //Considering only the central element's children
            benefit = 0.0;
            for (i = 0; i < Mesh_Small->nElemGroup; i++){
              for (j = 0; j < Mesh_Small->ElemGroup[i].nElem; j++){
                if (Tracer->GenArray[i].iValue[j][0] == 1){
                  benefit += AdaptIndicatorSmall->GenArray[i].rValue[j][0];
                }
              }
            }
            
            //evaluating the cost of the current refinement option
            CostRefOpt[RefOpt] = cost/benefit;
            //xf_pprintf("egrp: %d elem: %d RefOpt: %d Cost: %1.5g Benefit: %1.5g\n",egrp,elem,RefOpt,cost,benefit);
            //Cleaning-up
            ierr = xf_Error(xf_DestroySolverData(SolverData));
            if (ierr != xf_OK) return ierr;
            
            All_Small->EqnSet = NULL;
            
            ierr = xf_Error(xf_DestroyAll(All_Small));
            if (ierr != xf_OK) return ierr;
            
          }//RefOpt
        }
        //Choosing the optimal refinement option
        OptimumRef = 1;//start with isotropic H
        if (!Problem){
          temp = CostRefOpt[1];
          for (i = 2; i < nRefOpt[egrp]- ((LimitedP) ? nPRefOpt : 0); i++){
            if (CostRefOpt[i] < temp){
              temp = CostRefOpt[i];
              OptimumRef = i;
            }
          }
        }
        else {
          nNonPhysical++;
          //Setting the problematic cell to isotropic refinement
          xf_pprintf("xf_NON_PHYSICAL occurred on elem:%d egrp:%d for iRefOpt: %d\n",elem, egrp, RefOpt);
          xf_pprintf("Setting it to isotropic refinement\n");
          fflush(stdout);
          
          RefIndicator->GenArray[egrp].iValue[elem][0] = 1;
          //Cleaning up
          ierr = xf_Error(xf_DestroySolverData(SolverData));
          if (ierr != xf_OK) return ierr;
          
          All_Small->EqnSet = NULL;
          
          ierr = xf_Error(xf_DestroyAll(All_Small));
          if (ierr != xf_OK) return ierr;
        }
        RefIndicator->GenArray[egrp].iValue[elem][0] = OptimumRef;
        
        ierr = xf_Error(xf_Ref2nElem(Shape, OptimumRef, &nAdd));
        if (ierr != xf_OK) return ierr;
        nAdd--;//not considering the previously existing element
        
        xf_Release((void *)CostRefOpt);
        xf_Release((void *)nAddElem);
        
      }//index < nelemref
      else {
        nAdd = 0;
      }
      //only need to keep track of this when running with fixed growth
      if (ElemPicking == xfe_True){
        ierr = xf_Error(xf_MPI_Allreduce(&nAdd, 1, xfe_SizeInt, xfe_MPI_SUM));
        if (ierr != xf_OK) return ierr;
      }
      
      nAddElemTot += nAdd;
      
      if ((nAddElemTot+nElemTot)/nElemTot >= growth)
        if (ElemPicking == xfe_True)//Fixed growth method
          ReachedGrowth = xfe_True;
      
    }//not reached growth
    else
      RefIndicator->GenArray[egrp].iValue[elem][0] = 0;
    
  }//iref
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    ierr = xf_Error(xf_MPI_Allreduce(RefIndicator->GenArray[egrp].iValue[0], 
                                     RefIndicator->GenArray[egrp].n, 
                                     xfe_SizeInt, xfe_MPI_MAX));
    if (ierr != xf_OK) return ierr;
  }
  xf_MPI_Barrier();
  ierr = xf_Error(xf_MPI_Allreduce(&nNonPhysical, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Allreduce(&nLimitedP, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  //Refinement statistics
  //Do it on processor 0
  //only print statistics if there are more than one option
  if (!AdaptIsotropic || ConsiderPref){
    if (myRank == 0){
      /* sprintf(statsfname,"%s.stats");
       statsfile = fopen(statsfname, "a"); */
      xf_printf("xf_NON_PHYSICAL occurrences %d\n", nNonPhysical);
      ierr = xf_Error(xf_VAlloc2((void ***)&RefOptStat, Mesh->nElemGroup, nRefOpt, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      
      //zero-out the stats array
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++)
        for (RefOpt = 0; RefOpt < nRefOpt[egrp]; RefOpt++)
          RefOptStat[egrp][RefOpt] = 0;
      
      nRefineTot = 0;
      for (iref = 0; iref < nelemref; iref++){
        egrp = OrderToRefine[iref][0];
        elem = OrderToRefine[iref][1];
        
        RefOpt = RefIndicator->GenArray[egrp].iValue[elem][0];
        RefOptStat[egrp][RefOpt]++;
        
        if (RefOpt != 0)
          nRefineTot++;
      }
      if (ConsiderPref)
        xf_printf("%d elements were limited for P-ref\n", nLimitedP);
      //printing stats on the screen
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
        xf_printf("Refinement choices statistics for egrp %d:\n",egrp);
        for (RefOpt = 1; RefOpt < nHRefOpt[egrp]; RefOpt++){
          xf_printf("H-refinement option %d: %d/%d %3.1f%%\n", RefOpt, 
                    RefOptStat[egrp][RefOpt], nRefineTot, 100*(double)RefOptStat[egrp][RefOpt]/(double)nRefineTot);
        }
        for (RefOpt = 0; RefOpt < nPRefOpt; RefOpt++){
          xf_printf("P-refinement option %d: %d/%d %3.1f%%\n", RefOpt+1, 
                    RefOptStat[egrp][RefOpt+nHRefOpt[egrp]], nRefineTot, 100.0*(double)RefOptStat[egrp][RefOpt+nHRefOpt[egrp]]/(double)nRefineTot);
          fflush(stdout);
        }
        //fclose(statsfile);
        xf_Release2( (void **) RefOptStat);
        
        for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++)
          if (RefIndicator->GenArray[egrp].iValue[elem][0] >= nHRefOpt[egrp])
            //negative indicates p-refinement
            RefIndicator->GenArray[egrp].iValue[elem][0] = -OrderIncrement;
      }//egrp
      end = clock();
      elapsed = ((double)(end - start))/CLOCKS_PER_SEC;
    }//myRank == 0
    
    xf_printf("Anisotropic adaptation time statistics:\n",elapsed,elapsed/nRefineTot);
    xf_printf("Time: %1.3es ; %1.3es/elem\n",elapsed,elapsed/nRefineTot);
  }//!AdaptIsotropic or ConsiderPref
  //clean up on all processors
  xf_Release( (void *) nElem);
  xf_Release( (void *) nRefOpt);
  xf_Release( (void *) nHRefOpt);
  
  if (myRank != 0){
    //destroy mesh and dataset on non-root processors
    ierr = xf_Error(xf_DestroyDataSet(All->DataSet));
    if (ierr != xf_OK) return ierr;
    
    RefIndicator = NULL;
    AdaptIndicator = NULL;
    
    ierr = xf_Error(xf_DestroyMesh(All->Mesh));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ErrInd2ElemPos
static int 
xf_ErrInd2ElemPos(xf_All *All, xf_Vector *ErrIndicator, int ***pElemPos)
{
  /*
   PURPOSE:
   
   Determines the global position numbers of each element when sorted
   according to ErrIndicator.  Handles parallel.
   
   INPUTS:
   
   All : All structure
   ErrIndicator : real-valued vector of abs-value error estimates per elem
   
   OUTPUTS: 
   
   ElemPos[egrp][elem] = global position number of egrp,elem
   
   RETURN: Error code
   
   */
  int ierr, egrp, elem, negrp, sr, k;
  int nelemref, **ElemPos = NULL;
  int nelemtot, *nElem = NULL, *RI;
  int nelemtot_glob;
  real **Indicator = NULL;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  // used for sorting more than one value per element
  sr = ErrIndicator->StateRank;
  
  // create indicator vector over all elems (for sorting -> fixed fraction)
  ierr = xf_Error(xf_GetnElem(Mesh, &nElem, &nelemtot));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VAlloc2((void ***) &Indicator, negrp, nElem, sr*sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // copy over indicator data
  for (egrp=0; egrp<negrp; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      for (k=0; k<sr; k++)
        Indicator[egrp][sr*elem+k] = ErrIndicator->GenArray[egrp].rValue[elem][k];
  
  // create vector for storing element position post sorting
  ierr = xf_Error(xf_VAlloc2((void ***) pElemPos, negrp, nElem, sr*sizeof(int)));
  if (ierr != xf_OK) return ierr; 
  
  // sort indicator (ascending)
  ierr = xf_Error(xf_SortRealParallel(Indicator[0], sr*nelemtot, xfe_False, (*pElemPos)[0]));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) nElem);
  xf_Release2( (void **) Indicator);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ErrInd2RefInd
static int 
xf_ErrInd2RefInd(xf_All *All, xf_Vector *ErrIndicator, real frac, real fracCoarsen,
                 xf_Vector *RefIndicator, int ***pOrderToRefine, int *pnelemref,
                 int *pnelemcoarse)
{
  /*
   PURPOSE:
   
   Converts a real-valued error indicator, defined on each element of
   the mesh, to an integer-valued refinement indicator, also defined
   on each element.  Elements with the top fraction (frac) error
   indicators are flagged (RefIndicator = 1), while all other elements
   have RefIndicator = 0.  If fracCoarsen > 0, then elements with
   lowest fraction (fracCoarsen) error indicators are flagged for
   coarsening (RefIndicator = -1).
   
   INPUTS:
   
   All : All structure
   ErrIndicator : real-valued vector of error estimates per elem
   frac : fraction of elements to flag for refinement
   fracCoarsen : fraction of elements to flag for coarsening
   
   OUTPUTS: 
   
   RefIndicator : filled in (must be preallocated as an integer vector)
   OrderToRefine : Table that makes the elements to be refined in 
   decreasing order of error
   
   RETURN: Error code
   
   */
  int ierr, egrp, elem, elem_glob, negrp, i;
  int nelemref, **ElemPos = NULL;
  int nelemcoarse;
  int nelemtot, *RI;
  int nelemtot_glob;
  enum xfe_Bool AdaptRobust;
  enum xfe_AdaptOnType AdaptOn;
  xf_Mesh *Mesh;
  
  
  Mesh = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  //Check if we are using robustness adaptation
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                     "AdaptRobust", &AdaptRobust));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptOn", 
                                     xfe_AdaptOnName, (int ) xfe_AdaptOnLast, 
                                     (int *) &AdaptOn));
  if (ierr != xf_OK) return ierr;
  
  // create indicator vector over all elems (for sorting -> fixed fraction)
  ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;
  
  // Obtain ElemPos[egrp][elem] = global position number
  ierr = xf_Error(xf_ErrInd2ElemPos(All, ErrIndicator, &ElemPos));
  if (ierr != xf_OK) return ierr;
  
  // nelemtot_glob = global number of elements
  nelemtot_glob = nelemtot;
  ierr = xf_Error(xf_MPI_Allreduce(&nelemtot_glob, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  // number of elements to refine
  nelemref = nelemtot_glob*frac;
  
  (*pnelemref) = nelemref;
  
  if (pOrderToRefine != NULL){
    ierr = xf_Error(xf_Alloc2((void ***) pOrderToRefine, nelemref, 2, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    for (i = 0; i < nelemref; i++){//i is going to indicate the priority for refining
      (*pOrderToRefine)[i][0] = -1;//egrp
      (*pOrderToRefine)[i][1] = -1;//elem
    }
  }
  
  
  // flag all elements with ElemPos >= nelemtot_glob-nelemref; this is still isotropic
  for (egrp=0; egrp<negrp; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      RI = RefIndicator->GenArray[egrp].iValue[elem];
      RI[0] = 0;
      if (ElemPos[egrp][elem] >= (nelemtot_glob-nelemref)) {
        /*If we are doing robustness adaptation with penalty projection, 
          only the positive indicators should be considered*/
        if ((ErrIndicator->GenArray[egrp].rValue[elem][0] >= 0) || 
            ((ErrIndicator->GenArray[egrp].rValue[elem][0] < 0) && 
             (AdaptRobust != xfe_True) && (AdaptOn != xfe_AdaptOnPenalty))){
              
          RI[0] = 1;
          i = nelemtot_glob - ElemPos[egrp][elem]-1;
          if (Mesh->ParallelInfo != NULL)
            elem_glob = Mesh->ParallelInfo->ElemLoc2Glob[egrp][elem];
          else
            elem_glob = elem;
          if (pOrderToRefine != NULL){
            (*pOrderToRefine)[i][0] = egrp;
            (*pOrderToRefine)[i][1] = elem_glob;
          }
        }
      }      
    }
  }

  // flag all elements with ElemPos <= nelemcoarse for coarsening
  nelemcoarse = nelemtot_glob*fracCoarsen;

  if (pnelemcoarse != NULL) (*pnelemcoarse) = nelemcoarse;

  if (nelemcoarse > 0){
    for (egrp=0; egrp<negrp; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
        RI = RefIndicator->GenArray[egrp].iValue[elem];
        if (ElemPos[egrp][elem] <= nelemcoarse) {
          /*If we are doing robustness adaptation with penalty projection, 
            only the positive indicators should be considered*/
          if ((ErrIndicator->GenArray[egrp].rValue[elem][0] >= 0) || 
              ((ErrIndicator->GenArray[egrp].rValue[elem][0] < 0) && 
               (AdaptRobust != xfe_True) && (AdaptOn != xfe_AdaptOnPenalty))){
            RI[0] = -1;
          }
        }      
      }
    }
  }


  
  xf_Release2( (void **) ElemPos);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateRefIndicator
static int 
xf_CreateRefIndicator(xf_All *All, enum xfe_Bool AdaptIsotropic, 
                      enum xfe_Bool OrderRefFlag, enum xfe_Bool ConsiderCoarsen,
                      xf_Vector *AdaptIndicator, xf_Vector **pRefIndicator, 
                      int ***pOrderToRefine, int *nelemref, int *nelemcoarse)
{
  /* This function creates and fills the Refinement 
   indicator vector*/
  int ierr, negrp, egrp, elem, *RI, nelemtot, i;
  enum xfe_Bool AdaptSmoothRef;
  enum xfe_Bool ElemPicking;
  enum xfe_AdaptOnType AdaptOn;
  real frac, fracCoarsen, growth;
  xf_Vector *RefIndicator;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  // Find refinement indicator vector (integer for each element)
  ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem,
                                1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, 
                                xfe_True,  xfe_False, NULL, 
                                pRefIndicator, NULL));
  if (ierr != xf_OK) return ierr;
  
  RefIndicator = (*pRefIndicator);
  
  /* Determine what we are adapting on */
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptOn", 
                                     xfe_AdaptOnName, (int ) xfe_AdaptOnLast, 
                                     (int *) &AdaptOn));
  if (ierr != xf_OK) return ierr;
  
  // determine the element flagging method
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "AdaptFixedGrowth", 
                                     &ElemPicking));
  if (ierr != xf_OK) return ierr;
  
  //Uniform refinement//
  if (AdaptOn == xfe_AdaptOnUniform){ 
    // flag all elements for uniform refinement
    for (egrp=0; egrp<negrp; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
        RI = RefIndicator->GenArray[egrp].iValue[elem];
        RI[0] = 1;
      } 
    return xf_OK; // return immediately
  }
  
  //Fixed growth mode
  if (ElemPicking){ 
    // not yet supported for order refinement
    if (OrderRefFlag) return xf_Error(xf_NOT_SUPPORTED);
    
    // determine growth factor
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, 
                                       "AdaptFixedGrowthFactor", &growth));
    if (ierr != xf_OK) return ierr;
    
    if (AdaptIsotropic){
      //for now let's assume the same interpolation order in each element
      frac = (growth - 1.0)/(pow(2.0,Mesh->Dim) - 1.0);
    }
    else{
      /*for now let's assume single-cut refinement and 
       calculate a fraction based on that assumption.
       Note: This does not assume that the refined 
       elements are the worst amongst the ones marked 
       for refinement. Also this is the maximum number 
       of elements that have to be refined.*/
      frac = (growth - 1.0);
    }
  }
  else{
    // determine fixed fraction
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, 
                                       "AdaptFixedFraction", &frac));
    if (ierr != xf_OK) return ierr;
  }

  // Coarsening fraction if considering it
  if (ConsiderCoarsen){
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, 
                                       "AdaptCoarsenFraction", &fracCoarsen));
    if (ierr != xf_OK) return ierr;
  }
  else fracCoarsen = 0.;
  
  // Convert error indicator to refinement indicator
  ierr = xf_Error(xf_ErrInd2RefInd(All, AdaptIndicator, frac, fracCoarsen,
                                   RefIndicator, pOrderToRefine, nelemref, 
                                   nelemcoarse));
  if (ierr != xf_OK) return ierr;
  
  // smooth mesh refinement if desired
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "AdaptSmoothRef", 
                                     &AdaptSmoothRef));
  if (ierr != xf_OK) return ierr;
  
  if (AdaptSmoothRef){
    ierr = xf_Error(xf_AdaptSmoothRef(All, RefIndicator));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_PreAdaptHang
static int 
xf_PreAdaptHang(xf_All **pAll, xf_Vector **pAdaptIndicator, 
                xf_Vector **pRefIndicator, int **OrderToRefine, 
                int nelemref)
{
  int ierr, myRank, nProc;
  xf_Vector *RefIndicator_Glob, *AdaptIndicator_Glob;
  xf_DataSet *DataSet_Glob;
  xf_Mesh *Mesh_Glob;
  xf_Data *D;
  xf_All *All = (*pAll);
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  /* For now, serialize before adapting in parallel */
  if (nProc > 1){
    
    if (OrderToRefine != NULL){
      ierr = xf_Error(xf_MPI_Allreduce(OrderToRefine[0], nelemref*2, xfe_SizeInt, xfe_MPI_MAX));
      if (ierr != xf_OK) return ierr;
    }
    
    // RefIndicator is not in the data set
    if (myRank == 0){
      ierr = xf_Error(xf_CreateVector( &RefIndicator_Glob));
      if (ierr != xf_OK) return ierr;
    }
    
    ierr = xf_Error(xf_UnParallelizeVector(All->Mesh, (*pRefIndicator), RefIndicator_Glob));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0)
      (*pRefIndicator) = RefIndicator_Glob;
    else{
      ierr = xf_Error(xf_DestroyVector((*pRefIndicator), xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    
    // DataSet
    if (myRank == 0){
      ierr = xf_Error(xf_CreateDataSet( &DataSet_Glob));
      if (ierr != xf_OK) return ierr;
    }
    
    ierr = xf_Error(xf_UnParallelizeDataSet( All, All->DataSet, DataSet_Glob));
    if (ierr != xf_OK) return ierr;
    
    //this destroys AdaptIndicator in all processors
    ierr = xf_Error(xf_DestroyDataSet(All->DataSet));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0){
      All->DataSet = DataSet_Glob;
      ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "AdaptIndicator", xfe_Vector, &D));
      if (ierr != xf_OK) return ierr;
      D->ReadWrite = xfe_True;
      (*pAdaptIndicator) = (xf_Vector *) D->Data;
      
      ierr = xf_Error(xf_DataSetAdd(All->DataSet, "RefIndicator", xfe_Vector, xfe_True, 
                                    (void *)(*pRefIndicator), NULL));
      if (ierr != xf_OK) return ierr;
    }
    
    // Mesh
    if (myRank == 0){    
      ierr = xf_Error(xf_CreateMesh(&Mesh_Glob));
      if (ierr != xf_OK) return ierr;
    }
    
    ierr = xf_Error(xf_UnParallelizeMesh(All->Mesh, Mesh_Glob));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DestroyMesh(All->Mesh));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0) All->Mesh = Mesh_Glob;
  }  
  /* End of Serialization */
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_PostAdaptHang
static int 
xf_PostAdaptHang(xf_All *All, xf_Vector *RefIndicator, int **OrderToRefine)
{
  int ierr, myRank, nProc;
  xf_Mesh *Mesh_Glob;
  xf_DataSet *DataSet_Glob;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  xf_Release2((void **)OrderToRefine);
  
  if (nProc > 1){
    // Mesh
    if (myRank == 0) Mesh_Glob = All->Mesh;
    
    ierr = xf_Error(xf_CreateMesh(&All->Mesh));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ParallelizeMesh(Mesh_Glob, All->Mesh, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0){
      ierr = xf_Error(xf_DestroyMesh(Mesh_Glob));
      if (ierr != xf_OK) return ierr;
    }
    
    // DataSet
    if (myRank == 0){ 
      ierr = xf_Error(xf_DataSetRemove(All->DataSet, "RefIndicator", xfe_False));
      if (ierr != xf_OK) return ierr;
      
      DataSet_Glob = All->DataSet;
    }
    
    ierr = xf_Error(xf_CreateDataSet(&All->DataSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ParallelizeDataSet( All, DataSet_Glob, All->DataSet));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0){
      //this destroys AdaptIndicator and RefIndicator
      ierr = xf_Error(xf_DestroyDataSet(DataSet_Glob));
      if (ierr != xf_OK) return ierr;
    }
    
    // re-parallelize EqnSet
    ierr = xf_Error(xf_ReParallelizeEqnSet(All->EqnSet));
    if (ierr != xf_OK) return ierr;
    
  }
  else {
    ierr = xf_Error(xf_DestroyVector(RefIndicator, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
  
}

/******************************************************************/
//   FUNCTION Definition: xf_AdaptAll
int 
xf_AdaptAll(xf_All *All, int iAdapt, enum xfe_Bool *pDoneAdapt, 
            xf_Vector *Indicator)
{
  int ierr, egrp, elem, myRank, nProc, **OrderToRefine=NULL;
  int nelemtot, nelemref, nelemcoarse, s, sr, n, nn;
  enum xfe_Bool AdaptIsotropic, found, ProcViz, AdaptRobustIndBdown, ConsiderPref;
  enum xfe_AdaptOnType AdaptOn;
  enum xfe_AdaptMechType AdaptMechanics;
  char AdaptOutput[xf_MAXSTRLEN], AdaptVariableSet[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN], OutputFile[xf_MAXSTRLEN];
  char AdaptScalar[xf_MAXSTRLEN];
  real OutputError, AdaptTolerance, ElemVol, ResTolerance;
  real AdaptRobustCFLAmpFactor, MaxCFLAchieved;
  xf_Vector *AdaptIndicator, *RefIndicator, *U, *RefIndicatorP;
  xf_Data *D;
  xf_SolverData *SolverData;
  
  (*pDoneAdapt) = xfe_False;
  
  /* Delete non-essential data -- for memory cleanup */
  ierr = xf_Error(xf_DataSetDeleteNonEssential(All->DataSet));
  if (ierr != xf_OK) return ierr;
  
  /* Determine what we are adapting on */
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptOn", 
                                     xfe_AdaptOnName, (int ) xfe_AdaptOnLast, 
                                     (int *) &AdaptOn));
  if (ierr != xf_OK) return ierr;
  
  //should we consider P-refinement
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                     "AdaptIncludeP", 
                                     &ConsiderPref));
  if (ierr != xf_OK) return ierr;
  
  /* Is isotropic adaptation requested? */
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                     "AdaptIsotropic", &AdaptIsotropic));
  if (ierr != xf_OK) return ierr;
  
  /* ProcViz for debugging purposes */
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ProcViz", &ProcViz));
  if (ierr != xf_OK) return ierr;
  
  if (ProcViz){
    ierr = xf_Error(xf_ProcViz(All));
    if (ierr != xf_OK) return ierr;
  }
  
  /* Determine the adaptation mechanics */
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptMechanics", 
                                     xfe_AdaptMechName, (int ) xfe_AdaptMechLast, 
                                     (int *) &AdaptMechanics));
  if (ierr != xf_OK) return ierr;
  
  /* Determine adaptation tolerance */
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, 
                                     "AdaptTolerance", &AdaptTolerance));
  if (ierr != xf_OK) return ierr;
  
  
  /* Locate the adaptive indicator */
  if (Indicator == NULL){
    ierr = xf_Error(xf_FindAdaptIndicator(All, xfe_False, &AdaptIndicator));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_SetZeroVector(AdaptIndicator));
    if (ierr != xf_OK) return ierr;
    
    
    /*----------------------------------*/
    /* Calculate the adaptive indicator */
    /*----------------------------------*/
    
    if (AdaptOn == xfe_AdaptOnOutput){
      
      /* Determine which output we are adapting on */
      ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptOutput", AdaptOutput));
      if (ierr != xf_OK) return ierr;
      
      /* Determine entropy on which we are adapting*/
      ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptVariableSet", AdaptVariableSet));
      if (ierr != xf_OK) return ierr;
      
      /* Call error estimation if adapting on output; note:
       StoreFineSpace = !AdaptIsotropic
       ReuseFineSpace =  xfe_False
       */
      ierr = xf_Error(xf_ErrEstOutput(All, AdaptOutput, AdaptVariableSet, !AdaptIsotropic,
                                      xfe_False,  AdaptIndicator, xfe_False, &OutputError));
      if (ierr == xf_SOLVER_ERROR){
        xf_printf("xf_SOLVER_ERROR occurred in Error Estimation.\nSwitching to Residual-based adaptation.\n");
        AdaptOn = xfe_AdaptOnResidual;
        AdaptIsotropic = xfe_True;
      }
      else if (ierr != xf_OK) return ierr;
      else {
        xf_printf(" Output error estimate = %.10E\n", OutputError);
        
        if (fabs(OutputError) <= AdaptTolerance){
          xf_printf("\nAdaptation tolerance (%.10E) met. Exiting.\n", AdaptTolerance);
          (*pDoneAdapt) = xfe_True;
        }
      }
      
    }
    else if (AdaptOn == xfe_AdaptOnResidual){
      
      /* Call error estimation if adapting on residual; note:
       StoreFineSpace =  xfe_False
       ReuseFineSpace =  xfe_False
       */
      ierr = xf_Error(xf_ErrEstOutput(All, NULL, NULL, xfe_False, xfe_False,
                                      AdaptIndicator, xfe_False, &OutputError));
      if (ierr != xf_OK) return ierr;
      
      xf_printf(" Residual norm estimate = %.10E\n", OutputError);
      
      if (fabs(OutputError) <= AdaptTolerance){
        xf_printf("\nAdaptation tolerance (%.10E) met. Exiting.\n", AdaptTolerance);
        (*pDoneAdapt) = xfe_True;
      }
    }
    else if (AdaptOn == xfe_AdaptOnScalar){
      
      /* Determine which output we are adapting on */
      ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptScalar", AdaptScalar));
      if (ierr != xf_OK) return ierr;
      
      /* Call scalar-based adaptive indicator */
      ierr = xf_Error(xf_AdaptIndicatorScalar(All, AdaptScalar, AdaptIndicator));
      if (ierr != xf_OK) return ierr;
      
      xf_printf("\nNot checking adaptation tolerance since adapting on scalar.\n");
      
    }
    else if (AdaptOn == xfe_AdaptOnTest){
      
      /* A test adaptive indicator */
      ierr = xf_Error(xf_AdaptIndicatorTest(All, AdaptIndicator));
      if (ierr != xf_OK) return ierr;
      
    }
  }//Indicator == NULL
  else {
    AdaptIndicator = Indicator;
  }
  
  /* If SavePrefix is None or NULL, will not write anything */
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;
  
  if (xf_NotNull(SavePrefix)){
    /* Write out .xfa file before performing adaptation */
    sprintf(OutputFile, "%s_A%02d.xfa\0", SavePrefix, iAdapt);
    ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
    if (ierr!=xf_OK) return ierr;
  }
  
  if (*pDoneAdapt) return xf_OK;
  
  
  
  /*------------------------*/
  /* Perform the adaptation */
  /*------------------------*/
  
  
  ierr = xf_Error(xf_GetnElem(All->Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Allreduce(&nelemtot, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  xf_printf(" Current number of elements = %d.\n", nelemtot);
  
  xf_printf(" Adapting the mesh.\n");
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (AdaptMechanics == xfe_AdaptMechHangNode){
    /*Process:
     0-Create the refinement indicator
     1-Serialize the mesh and data.
     2-Calculate the direction of refinement 
     (includes p-ref as an option)
     3-Clean-up error estimation vectors
     4-Separate h and p refinement indicators
     5-Refine mesh/order using hanging-node/vorder.
     6-Parallelize mesh and data.
     */
    
    //0-Create refinement indicator
    ierr = xf_Error(xf_CreateRefIndicator(All, AdaptIsotropic, xfe_False, xfe_False,
                                          AdaptIndicator, &RefIndicator, 
                                          &OrderToRefine, &nelemref, NULL));
    if (ierr != xf_OK) return ierr;
    
    //1-Serialize the mesh and data.
    ierr = xf_Error(xf_PreAdaptHang(&All, &AdaptIndicator, &RefIndicator, 
                                    OrderToRefine, nelemref));
    if (ierr != xf_OK) return ierr;
    
    
    
    //2-Calculate the direction of refinement (if requested;includes p-ref)
    ierr = xf_Error(xf_OptimizeHangAnisotropic(All, AdaptIndicator, RefIndicator, 
                                               OrderToRefine, nelemref, AdaptIsotropic));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0){  
      //3-Clean-up error estimation vectors
      ierr = xf_Error(xf_ErrEstOutput(All, NULL, NULL, 
                                      xfe_False,xfe_True, NULL, xfe_False, NULL));
      if (ierr != xf_OK) return ierr;
      
      //4-Separate Refinement indicators
      if (ConsiderPref){
        ierr = xf_Error(xf_FindVector(All, "RefIndicatorP", xfe_LinkageGlobElem,
                                      1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, 
                                      xfe_SizeInt, xfe_False,  xfe_True, &D, 
                                      &RefIndicatorP, NULL));
        if (ierr != xf_OK) return ierr;
        D->ReadWrite = xfe_True;
        
        for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++){
          for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
            if (RefIndicator->GenArray[egrp].iValue[elem][0] < 0){
              //negative values correspond to p-refinement
              RefIndicatorP->GenArray[egrp].iValue[elem][0] = 
              -RefIndicator->GenArray[egrp].iValue[elem][0];
              RefIndicator->GenArray[egrp].iValue[elem][0] = 0;
            }
            else
              RefIndicatorP->GenArray[egrp].iValue[elem][0] = 0;
          }
        }
      }
      //5-Refine mesh using hanging-node/p-refinement.
      if (ConsiderPref){
        //reseting pointer after parallelization 
        ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "RefIndicatorP", xfe_Vector, &D));
        if (ierr != xf_OK) return ierr;
        RefIndicatorP = (xf_Vector *) D->Data;
        
        ierr = xf_Error(xf_AdaptOrderRef(All, SavePrefix, RefIndicatorP, NULL, NULL, NULL));
        if (ierr != xf_OK) return ierr;
        // destroy RefIndicator
        ierr = xf_Error(xf_DataSetRemove(All->DataSet, "RefIndicatorP", xfe_False));
        if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
      if (ierr != xf_OK) return ierr;
    }//myRank == 0
    xf_MPI_Barrier();
    //6-Parallelize the mesh and data.
    ierr = xf_Error(xf_PostAdaptHang(All, RefIndicator, OrderToRefine));
    if (ierr != xf_OK) return ierr;
  }
  else if (AdaptMechanics == xfe_AdaptMechOrderRef){
    /*Process:
     0-Create the refinement indicator
     1-Refine the mesh using order increase, all in parallel
     */
    
    // only adapt isotropically
    if (!AdaptIsotropic) return xf_Error(xf_NOT_SUPPORTED);
    
    // 0-Create the refinement indicator
    ierr = xf_Error(xf_CreateRefIndicator(All, xfe_True, xfe_True, xfe_True,
                                          AdaptIndicator, &RefIndicator, 
                                          NULL, &nelemref, &nelemcoarse));
    if (ierr != xf_OK) return ierr;
    
    xf_printf("%d elements chosen for p coarsening.\n", nelemcoarse);
    xf_printf("%d elements chosen for p refinement.\n", nelemref);
    
    // 1 refine the mesh
    ierr = xf_Error(xf_AdaptOrderRef(All, SavePrefix, RefIndicator, NULL, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    
    // destroy RefIndicator
    ierr = xf_Error(xf_DestroyVector(RefIndicator, xfe_True));
    if (ierr != xf_OK) return ierr;
    
  }
  else return xf_Error(xf_NOT_SUPPORTED);
  
  ierr = xf_Error(xf_GetnElem(All->Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Allreduce(&nelemtot, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  xf_printf(" New number of elements = %d.\n", nelemtot);
  
  if (xf_NotNull(SavePrefix)){
    /* Write out .xfa file after performing adaptation */
    sprintf(OutputFile, "%s_B%02d.xfa\0", SavePrefix, iAdapt+1);
    ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
    if (ierr!=xf_OK) return ierr;
  }
  
  return xf_OK;
}



/*-----------------------*/
/*  UNSTEADY ADAPTATION  */
/*-----------------------*/


/******************************************************************/
//   FUNCTION Definition: xf_ReadTemporalErrorSerial
static int
xf_ReadTemporalErrorSerial(const char *fname, const char *OutputName,
                           int nSlab, real *ErrIndTime, int *sdofTime)
{
  // This is called on the root processor to perform the actual file reading
  // ErrIndTime must be pre-allocated
  int ierr;
  int iSlab;
  int i, col, nOutput=0;
  char line[xf_MAXLINELEN];
  char **OutputNames=NULL;
  real *rbuf = NULL;
  FILE *fid;
  
  if ((fid = fopen(fname, "r")) == NULL)
    return xf_Error(xf_FILE_READ_ERROR);
  
  // find header with output names, col = column of our output
  col = -1;
  do{
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) continue;
    
    if (strncmp(line, "% iSlab DOF", 11)==0){ /* header line */
      ierr = xf_Error(xf_ScanXStringAlloc(line+11, xf_MAXSTRLEN, &nOutput, &OutputNames));
      if (ierr != xf_OK) return ierr;
      // check if matches the output we're looking for
      for (i=0; i<nOutput; i++)
        if (strcmp(OutputNames[i], OutputName) == 0){
          col = i+2; // adding 2 to account for iSlab and DOF
          break; 
        }
      xf_Release2( (void **) OutputNames);
    }
    if (col > 0) break;
  } while(feof(fid) == 0);	
  if (col < 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // allocate a real buffer
  ierr = xf_Error(xf_Alloc( (void **) &rbuf, nOutput+2, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // read in error estimates on each slab
  iSlab = 0;
  do{
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) continue;
    // each line contains nOutput+2 reals/integers
    ierr = xf_Error(xf_ScanReal(line, nOutput+2, rbuf));
    if (ierr != xf_OK) return ierr;
    if (iSlab >= nSlab) return xf_Error(xf_OUT_OF_BOUNDS);
    sdofTime[iSlab] = (int) rbuf[1]; // read in as a real, but dof is an int
    ErrIndTime[iSlab++] = rbuf[col];
  } while(feof(fid) == 0);
  if (iSlab != nSlab) return xf_Error(xf_FILE_READ_ERROR);
  
  // close file
  fclose(fid);
  
  // release buffer
  xf_Release( (void *) rbuf);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadTemporalError
static int
xf_ReadTemporalError(const char *fname, const char *OutputName,
                     int nSlab, real **pErrIndTime, int **psdofTime)
{
  /*
   
   PURPOSE:
   
   Reads temporally-localized error indicator from a text file, for one
   Output (the file may have more outputs)
   
   INPUTS:
   
   fname       : file name from which to read
   OutputName  : name of output error to read
   nSlab       : number of time slabs
   
   OUTPUTS: 
   
   (*pErrIndTime)  : nSlab vector of temporally-localized indicators
   (*psdofTime)    : nSlab vector of spatial degrees of freedom
   
   RETURN: Error code
   
   */
  int ierr, myRank;
  real *ErrIndTime;
  int *sdofTime;
  
  /* Determine myRank */
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  // all procs allocate
  ierr = xf_Error(xf_Alloc( (void **) pErrIndTime, nSlab, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) psdofTime, nSlab, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ErrIndTime = (*pErrIndTime);
  sdofTime = (*psdofTime);
  
  // only root reads
  if (myRank == 0)
    ierr = xf_Error(xf_ReadTemporalErrorSerial(fname, OutputName, nSlab, ErrIndTime, sdofTime));
  if (xf_PError(&ierr, 0) != xf_OK) return ierr;
  
  // broadcast ErrIndTime from root
  ierr = xf_Error(xf_MPI_Bcast((void *) ErrIndTime, nSlab*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  // broadcast sdofTime from root
  ierr = xf_Error(xf_MPI_Bcast((void *) sdofTime, nSlab*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadSpatialError
static int
xf_ReadSpatialError(xf_All *All, const char *SavePrefix, const char *OutputName,
                    xf_Vector *ErrIndSpace, xf_Vector **pBasisOrderElem)
{
  /*
   
   PURPOSE:
   
   Reads spatially-localized error indicators from files, one .data
   file for each time slab.  If the state is also present on disk,
   reads dof information from the state vectors.
   
   Also called in non-dynamics case, to allocate and fill in pBasisOrderElem
   
   INPUTS:
   
   All         : All structure
   SavePrefix  : prefix for saved files
   OutputName  : name of output of interest
   
   OUTPUTS: 
   
   ErrIndSpace  : vector of nSlab indicators for each element
   (*pBasisOrderElem) : basis and order information for each element (at each time)
   
   RETURN: Error code
   
   */
  int ierr, nSlab, iSlab, i, j;
  int Order;
  char fname[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  enum xfe_BasisType Basis;
  xf_DataSet *DataSet = NULL;
  xf_Data *D;
  xf_Vector *ErrInd, *U, *BasisOrderElem;
  
  nSlab = ErrIndSpace->StateRank;
  
  // allocate a vector of element basis, order
  ierr = xf_Error(xf_FindVector(All, "BasisOrderElem", xfe_LinkageGlobElem, 2*nSlab, NULL, 0, 0, 
                                NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_True,
                                xfe_False, NULL, pBasisOrderElem, NULL));
  if (ierr != xf_OK) return ierr;
  
  BasisOrderElem = (*pBasisOrderElem);
  
  
  for (iSlab=0; iSlab<nSlab; iSlab++){
    if (nSlab > 1){
      // create data set for reading
      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;
      // read .data from file 
      sprintf(fname, "%s_ErrS_%s%d.data\0", SavePrefix, OutputName, iSlab);
      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, fname, DataSet));
      if (ierr != xf_OK) return ierr;
      
      D = DataSet->Head;
      ErrInd = (xf_Vector *) D->Data;
      
      // take absolute values, since errors may be signed
      ierr = xf_Error(xf_VectorAbs(ErrInd));
      if (ierr != xf_OK) return ierr;
      
      // copy over error indicators into ErrIndSpace
      for (i=0; i<ErrIndSpace->nArraySelf; i++){
        for (j=0; j<ErrIndSpace->GenArray[i].n; j++){
          if (ErrIndSpace->GenArray[i].r != nSlab) return xf_Error(xf_CODE_LOGIC_ERROR);
          ErrIndSpace->GenArray[i].rValue[j][iSlab] = ErrInd->GenArray[i].rValue[j][0];
        } // j
      } // i
      
      // Destroy DataSet
      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;
    }
    
    // read state (U) .data from file (multiple vectors in one dataset)
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;
    sprintf(Title, "%s_U%d.data\0", SavePrefix, iSlab+1);
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, Title, DataSet));
    if (ierr != xf_OK) return xf_Error(ierr);
    
    D = DataSet->Head;
    U = (xf_Vector *) D->Data;
    
    for (i=0; i<BasisOrderElem->nArraySelf; i++){
      Basis = U->Basis[i];
      for (j=0; j<BasisOrderElem->GenArray[i].n; j++){
        Order = xf_InterpOrder(U, i, j);
        BasisOrderElem->GenArray[i].iValue[j][0*nSlab+iSlab] = Basis;
        BasisOrderElem->GenArray[i].iValue[j][1*nSlab+iSlab] = Order;
      } // j
    } // i
    
    // Destroy DataSet
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;
        
  } // iSlab  
  
  // BasisOrderElem needs to have data in halos when running in parallel
  if (All->Mesh->ParallelInfo != NULL){
    // begin communication of halo data
    ierr = xf_Error(xf_HaloExchangeVectorBegin(BasisOrderElem));
    if (ierr != xf_OK) return ierr;
    // end communication of halo data
    ierr = xf_Error(xf_HaloExchangeVectorEnd(BasisOrderElem));
    if (ierr != xf_OK) return ierr;
  }
  

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_OptimizeRefUnsteadyOneShot
static int 
xf_OptimizeRefUnsteadyOneShot(xf_All *All, xf_TimeHistData *TimeHistData, 
                              const char *OutputName, xf_Vector **pRefIndicatorSpace,
                              int *RefIndicatorTime, xf_Vector **pBasisOrderElem) 
{
  /*
   PURPOSE:
   
   Fills in RefIndicatorSpace and RefIndicatorTime (the adaptive
   refinement indicators), using output error estimates in space and
   time.  A "one-shot" approach is taken, in which a single mesh is
   used for all time, so that the spatial adaptation targets the
   worst-behaving elements after marginalizing (summing) over the time
   coordinate.  Similarly for the temporal adaptation.
   
   Modified after initial trials so that fixed fraction is not applied
   independently to the space direction and the time direction.
   Instead, a fixed fraction or fixed growth of *space-time* elements
   is enforced.
   
   INPUTS:
   
   All               : xf_All structure
   TimeHistData      : time history data
   OutputName        : output on which we are adapting
   
   OUTPUTS: 
   
   RefIndicatorSpace : spatial error indicator (preallocated)
   RefIndicatorTime  : temporal error indicator (preallocated)
   
   RETURN: Error code
   
   */
  int ierr, i, nref=0;
  int nTime, iTime, ielemtot;
  int egrp, elem, ielem;
  int posTime, posSpace;
  int Budget, Marked, Freed;
  int *StepPos = NULL;
  int *Pos2Step = NULL;
  int **ElemPos = NULL;
  int *LocPos = NULL;
  int **LocPos2Elem = NULL;
  int *sdofTime = NULL;
  int nelemtot, nelemtot_glob;
  int Order, EISrank, ddof, nn, iEIS;
  int AdaptOrderMin, AdaptOrderMax;
  int ndoftot, itot;
  int BudgetCoarsen;
  enum xfe_AdaptMechType AdaptMechanics;
  enum xfe_Bool Growth, Found, done;
  enum xfe_Bool DynamicSpatialRef;
  enum xfe_ShapeType Shape;
  enum xfe_BasisType Basis;
  char SavePrefix[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  real frac, errspace, fac, fracCoarsen;
  real facSlabFOM;
  real *ErrIndTime = NULL;
  xf_Vector *ErrIndSpaceTot = NULL;
  xf_Vector *BasisOrderElem = NULL;
  xf_Vector *RefIndicatorSpace = NULL;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  nTime = TimeHistData->nTime;
  
  // zero out RefIndicatorTime
  for (i=0; i<nTime; i++) RefIndicatorTime[i] = 0;
  
  // Get number of elements (total too)
  ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;
  
  // nelemtot_glob = global number of elements
  nelemtot_glob = nelemtot;
  ierr = xf_Error(xf_MPI_Allreduce(&nelemtot_glob, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  // determine the element flagging method
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "AdaptFixedGrowth", &Growth));
  if (ierr != xf_OK) return ierr;
  
  // pull off SavePrefix
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;
  
  // Read temporally-localized error estimate -> ErrIndTime
  sprintf(OutputFile, "%s_TemporalError.txt\0", SavePrefix);
  ierr = xf_Error(xf_ReadTemporalError(OutputFile, OutputName, nTime, &ErrIndTime, &sdofTime));
  if (ierr != xf_OK) return ierr;
  
  // modify ErrIndTime to be a benefit function (bang for buck)
  // bang = error addressed, buck = degrees of freedom associated with slab refinement
  for (i=0; i<nTime; i++) ErrIndTime[i] /= ((real) sdofTime[i]);
  
  // total number of degrees of freedom, pre-adaptation
  for (i=0,ndoftot=0; i<nTime; i++) ndoftot += sdofTime[i];
  xf_printf("Total number of degrees of freedom before adapt = %d\n", ndoftot);
  
  // determine the DOF budget 
  if (Growth){
    // fixed growth
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "AdaptFixedGrowthFactor", &frac));
    if (ierr != xf_OK) return ierr;
    // Calculate Budget = total # additional degrees of freedom
    Budget = (frac-1.0)*ndoftot; // fixed growth
  }
  else{
    // fixed fraction .. treat same as dof growth
    xf_printf("Applying fixed fraction to total degrees of freedom.\n");
    
    // determine fixed fraction
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "AdaptFixedFraction", &frac));
    if (ierr != xf_OK) return ierr;
    
    // Calculate Budget = total # additional degrees of freedom
    Budget = (frac-1.0)*ndoftot; // fixed growth
  }
  
  // minimum order for p-adaptation
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "AdaptOrderMin", &AdaptOrderMin));
  if (ierr != xf_OK) return ierr;

  // maximum order for p-adaptation
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "AdaptOrderMax", &AdaptOrderMax));
  if (ierr != xf_OK) return ierr;

  // Coarsening fraction
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "AdaptCoarsenFraction", &fracCoarsen));
  if (ierr != xf_OK) return ierr;

  // budget for coarsening = # dofs we want freed up during coarsening phase
  BudgetCoarsen = ndoftot*fracCoarsen;


  // allocate a vector for storing the time step position in sorted array
  ierr = xf_Error(xf_Alloc( (void **) &StepPos, 2*nTime, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  Pos2Step = StepPos+nTime;
  
  // Sort ErrIndTime (ascending)
  ierr = xf_Error(xf_SortRealParallel(ErrIndTime, nTime, xfe_True, StepPos));
  if (ierr != xf_OK) return ierr;
  
  // Form Pos2Step[i] = time slab index of ith lowest error
  for (i=0; i<nTime; i++) Pos2Step[StepPos[i]] = i;  
  
  /*   for (i=0; i<nTime; i++) */
  /*     xf_printf("%d TimeSlab = %d: %.10E\n", i, Pos2Step[i], ErrIndTime[i]); */
  
  // are we doing dynamic spatial refinement?  (only order ref supported for now)
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "DynamicSpatialRef", 
                                     &DynamicSpatialRef));
  if (ierr != xf_OK) return ierr;
  
  // rank of spatial error indicator vector: either 1 or nTime values per element
  EISrank = ((DynamicSpatialRef) ? nTime : 1);
  
  // Allocate refinement indicator in space
  ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem, EISrank, NULL, 0, 0, 
                                NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_True,
                                xfe_False, NULL, pRefIndicatorSpace, NULL));
  if (ierr != xf_OK) return ierr;
  RefIndicatorSpace = (*pRefIndicatorSpace);
  
  // zero out RefIndicatorSpace
  ierr = xf_Error(xf_SetZeroVector(RefIndicatorSpace));
  if (ierr != xf_OK) return ierr;
  
  
  
  // Locate spatial error indicator
  sprintf(Title, "ErrIndSpaceTot_%s", OutputName);
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, EISrank, NULL, 0, 0, 
                                NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,
                                xfe_True, NULL, &ErrIndSpaceTot, &Found));
  if (ierr != xf_OK) return ierr;
  if ((!DynamicSpatialRef) && (!Found)) return xf_Error(xf_NOT_FOUND);
  
  // in dynamic spatial refinement, read in individual error indicators at each time
  // calling not just fr dynamicspatialref, because will need BasisOrderElem in general
  ierr = xf_Error(xf_ReadSpatialError(All, SavePrefix, OutputName, ErrIndSpaceTot, pBasisOrderElem));
  if (ierr != xf_OK) return ierr;
  BasisOrderElem = (*pBasisOrderElem);
  
  
  /* Determine the adaptation mechanics */
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptMechanics", 
                                     xfe_AdaptMechName, (int ) xfe_AdaptMechLast, 
                                     (int *) &AdaptMechanics));
  if (ierr != xf_OK) return ierr;


  // no coarsening if using hanging-node (not yet supported)
  if (AdaptMechanics == xfe_AdaptMechHangNode){
    fracCoarsen = 0.0;
    BudgetCoarsen = 0;
  }
  
  
  // Modify ErrIndSpaceTot to a benefit function: divide by additional dof to be introduced
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      if (AdaptMechanics == xfe_AdaptMechHangNode){ // for hanging node
        ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_Ref2nElem(Shape, 1, &nref));
        if (ierr != xf_OK) return ierr;
        Basis = BasisOrderElem->GenArray[egrp].iValue[elem][0];
        Order = BasisOrderElem->GenArray[egrp].iValue[elem][1];
        ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
        if (ierr != xf_OK) return ierr;
        fac = (real) (nTime/EISrank*(nref-1)*nn);
        if (fac <= 0.) return xf_Error(xf_OUT_OF_BOUNDS);
        for (i=0; i<EISrank; i++)
          ErrIndSpaceTot->GenArray[egrp].rValue[elem][i] /= fac;
        
      }
      else if (AdaptMechanics == xfe_AdaptMechOrderRef){ // for order refinement
        for (i=0; i<EISrank; i++){
          Basis = BasisOrderElem->GenArray[egrp].iValue[elem][0*EISrank+i];
          Order = BasisOrderElem->GenArray[egrp].iValue[elem][1*EISrank+i];
          ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_Order2nNode(Basis, Order+1, &ddof));
          if (ierr != xf_OK) return ierr;
          ddof -= nn;
          fac = (real) ddof*nTime/EISrank;
          if (fac <= 0.) return xf_Error(xf_OUT_OF_BOUNDS);
          ErrIndSpaceTot->GenArray[egrp].rValue[elem][i] /= fac;
        } 
      }
      else return xf_Error(xf_NOT_SUPPORTED);
    } // elem
  } // egrp
  
  // Obtain ElemPos[egrp][elem] = global position number
  ierr = xf_Error(xf_ErrInd2ElemPos(All, ErrIndSpaceTot, &ElemPos));
  if (ierr != xf_OK) return ierr;
  
  // Form LocPos2Elem[nelemtot][4] = locally-sorted positions
  ierr = xf_Error(xf_Alloc2( (void ***) &LocPos2Elem, EISrank*nelemtot, 4, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &LocPos, EISrank*nelemtot, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SortIntPos(ElemPos[0], EISrank*nelemtot, LocPos, xfe_True));
  if (ierr != xf_OK) return ierr;
  for (egrp=0, itot=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      for (i=0; i<EISrank; i++, itot++){
        LocPos2Elem[itot][0] = ElemPos[0][itot]; // already sorted
        LocPos2Elem[LocPos[itot]][1] = egrp;
        LocPos2Elem[LocPos[itot]][2] = elem;
        LocPos2Elem[LocPos[itot]][3] = i;
      }
    }
  xf_Release( (void *) LocPos);
  xf_Release2( (void **) ElemPos);
  
  /*   for (ielemtot=0; ielemtot<nelemtot; ielemtot++) */
  /*     xf_printf("%d %d %d : %.10E\n", LocPos2Elem[ielemtot][0],  */
  /* 	      LocPos2Elem[ielemtot][1], LocPos2Elem[ielemtot][2], */
  /* 	      ErrIndSpaceTot->GenArray[LocPos2Elem[ielemtot][1]].rValue[LocPos2Elem[ielemtot][2]][0]); */
  
  /*------------*/
  /* COARSENING */
  /*------------*/

  Freed = 0;
  facSlabFOM = 2.0;  // upon coarsening, time slabs free up only half their dof
  
  if (BudgetCoarsen > 0){

    xf_printf("Beginning coarsening.\n");

    // set posTime=0, posSpace=0   (indices of lowest figure-of-merit)
    posTime  = 0;
    posSpace = 0;
    ielem    = 0;
    // begin while !done loop
    done = xfe_False;
    xf_printf("Coarsening Budget (# dofs we want freed) = %d\n", BudgetCoarsen);
    while (!done){
      // calculate next lowest space elem error
      if ((ielem >= 0) && (ielem < EISrank*nelemtot) && (LocPos2Elem[ielem][0] == posSpace)){
        egrp = LocPos2Elem[ielem][1];
        elem = LocPos2Elem[ielem][2];
        iEIS = LocPos2Elem[ielem][3];
        errspace = ErrIndSpaceTot->GenArray[egrp].rValue[elem][iEIS];
      }
      else{
        egrp = elem = iEIS = -1;
        errspace = -1.0;
      }
      // broadcast space elem error to all procs
      ierr = xf_Error(xf_MPI_Allreduce(&errspace, 1, xfe_SizeReal, xfe_MPI_MAX));
      if (ierr != xf_OK) return ierr;
      if ((posTime < 0) && (errspace < 0)) return xf_Error(xf_CODE_LOGIC_ERROR); // sanity check
    
      // compare, choose, increment Freed
      if ((posTime >= 0) && (ErrIndTime[posTime]*facSlabFOM < errspace)){
        // calculate index of next lowest time slab error
        iTime = Pos2Step[posTime];
        // coarsening of time slab is chosen
        RefIndicatorTime[iTime] = -1;
        Freed += ((Growth) ? 0.5*sdofTime[iTime] : 1); // freeing up 0.5 if coarsen by factor of 2
        posTime++;
        xf_printf("Coarsening time slab (iTime=%d). Freed = %d, BudgetCoarsen = %d\n", 
                  iTime, Freed, BudgetCoarsen);
      }
      else if (posSpace >= 0){
        // refinement of elem is chosen
        posSpace++;
        ddof = 0;
        if (egrp >= 0){ // only processor that contains egrp,elem coarsens
          RefIndicatorSpace->GenArray[egrp].iValue[elem][iEIS] = -1; // indicates coarsening
          ielem++;
          if (AdaptMechanics == xfe_AdaptMechHangNode){
            return xf_Error(xf_NOT_SUPPORTED);
          }
          else if (AdaptMechanics == xfe_AdaptMechOrderRef){
            Basis = BasisOrderElem->GenArray[egrp].iValue[elem][0*EISrank+iEIS];
            Order = BasisOrderElem->GenArray[egrp].iValue[elem][1*EISrank+iEIS];
            ierr = xf_Error(xf_Order2nNode(Basis, Order, &ddof));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_Order2nNode(Basis, max(AdaptOrderMin, Order-1), &nn));
            if (ierr != xf_OK) return ierr;
            ddof -= nn;
            ddof *= nTime/EISrank;
          }
          else return xf_Error(xf_NOT_SUPPORTED);
        }
        // broadcast number of new elements per spatial element to all procs
        ierr = xf_Error(xf_MPI_Allreduce(&ddof, 1, xfe_SizeInt, xfe_MPI_MAX));
        if (ierr != xf_OK) return ierr;
        Freed += ((Growth) ? ddof : 1);
        xf_printf("Choosing spatial coarsening. Freed = %d, BudgetCoarsen = %d\n", 
                  Freed, BudgetCoarsen);
      }
      else return xf_Error(xf_CODE_LOGIC_ERROR);
    
      // check if over BudgetCoarsen
      if ((Freed >= BudgetCoarsen) 
          || ((posTime>=nTime) && (posSpace>=(EISrank*nelemtot_glob)))) 
        done = xfe_True;
    } // end while not done
  
    xf_printf("Finished coarsening: nelemtot_glob=%d, nTime=%d, Freed = %d, BudgetCoarsen = %d\n", 
              nelemtot_glob, nTime, Freed, BudgetCoarsen);

  } // end if BudgetCoarsen > 0

  /*------------*/
  /* REFINEMENT */
  /*------------*/

  Budget += Freed;  // add freed-up dofs to refinement budget

  if (Budget > 0){

    xf_printf("Beginning refinement.\n");

    // set posTime=nTime-1, posSpace=nelemtot_glob-1
    posTime  = nTime-1;
    posSpace = EISrank*nelemtot_glob-1;
    ielem    = EISrank*nelemtot-1;
    // begin while !done loop
    done = xfe_False;
    Marked = 0;
    xf_printf("Refinement Budget (# dofs we can mark) = %d\n", Budget);
    while (!done){
      // calculate next highest space elem error
      if ((ielem >= 0) && (LocPos2Elem[ielem][0] == posSpace)){
        egrp = LocPos2Elem[ielem][1];
        elem = LocPos2Elem[ielem][2];
        iEIS = LocPos2Elem[ielem][3];
        errspace = ErrIndSpaceTot->GenArray[egrp].rValue[elem][iEIS];
      }
      else{
        egrp = elem = iEIS = -1;
        errspace = -1.0;
      }
      // broadcast space elem error to all procs
      ierr = xf_Error(xf_MPI_Allreduce(&errspace, 1, xfe_SizeReal, xfe_MPI_MAX));
      if (ierr != xf_OK) return ierr;
      if ((posTime < 0) && (errspace < 0)) return xf_Error(xf_CODE_LOGIC_ERROR); // sanity check
    
      // compare, choose, increment Marked
      if ((posTime >= 0) && (ErrIndTime[posTime] > errspace)){
        // calculate index of next highest time slab error
        iTime = Pos2Step[posTime];
        // refinement in time slab is chosen
        RefIndicatorTime[iTime] = 1;
        Marked += ((Growth) ? sdofTime[iTime] : 1);
        posTime--;
        xf_printf("Choosing temporal refinement (iTime=%d). Marked = %d, Budget = %d\n", 
                  iTime, Marked, Budget);
      }
      else if (posSpace >= 0){
        // refinement of elem is chosen
        posSpace--;
        ddof = 0;
        if (egrp >= 0){ // only processor that contains egrp,elem refines
          RefIndicatorSpace->GenArray[egrp].iValue[elem][iEIS] = 1;
          ielem--;
          if (AdaptMechanics == xfe_AdaptMechHangNode){
            ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_Ref2nElem(Shape, 1, &nref)); // assuming hanging node
            if (ierr != xf_OK) return ierr;
            Basis = BasisOrderElem->GenArray[egrp].iValue[elem][0];
            Order = BasisOrderElem->GenArray[egrp].iValue[elem][1];
            ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
            if (ierr != xf_OK) return ierr;
            ddof = (nref-1)*nTime*nn;
          }
          else if (AdaptMechanics == xfe_AdaptMechOrderRef){
            Basis = BasisOrderElem->GenArray[egrp].iValue[elem][0*EISrank+iEIS];
            Order = BasisOrderElem->GenArray[egrp].iValue[elem][1*EISrank+iEIS];
            ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_Order2nNode(Basis, min(AdaptOrderMax, Order+1), &ddof));
            if (ierr != xf_OK) return ierr;
            ddof -= nn;
            ddof *= nTime/EISrank;
          }
          else return xf_Error(xf_NOT_SUPPORTED);	
          //printf("  (egrp=%d, elem=%d, iTime=%d)\n", egrp, elem, iEIS);
        }
        // broadcast number of new elements per spatial element to all procs
        ierr = xf_Error(xf_MPI_Allreduce(&ddof, 1, xfe_SizeInt, xfe_MPI_MAX));
        if (ierr != xf_OK) return ierr;
        Marked += ((Growth) ? ddof : 1);
        xf_printf("Choosing spatial refinement. Marked = %d, Budget = %d\n", Marked, Budget);
      }
      else return xf_Error(xf_CODE_LOGIC_ERROR);
    
      // check if over Budget
      if ((Marked >= Budget) || ((posTime<0) && (posSpace<0))) done = xfe_True;
    } // end while not done
  
    xf_printf("Finished refining: nelemtot_glob=%d, nTime=%d, Marked = %d, Budget = %d\n", 
              nelemtot_glob, nTime, Marked, Budget);

  } // Budget

  
  // release memory
  xf_Release( (void *) sdofTime);
  xf_Release( (void *) ErrIndTime);
  xf_Release( (void *) StepPos);
  xf_Release2((void **) LocPos2Elem);
  
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SmoothTimeHistData
static int 
xf_SmoothTimeHistData(xf_TimeHistData *TimeHistData, real fac) 

{
  int ierr, i;
  int nTime;
  real *Time = NULL;
  
  // number of time steps (slabs)
  nTime = TimeHistData->nTime;

  // Allocate Time array
  ierr = xf_Error(xf_Alloc( (void **) &Time, nTime+1, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // copy over time
  for (i=0; i<nTime; i++) Time[i] = TimeHistData->Time[i];
  Time[nTime] = Time[nTime-1] + TimeHistData->TimeStep[nTime-1];

  // copy back smoothed version
  for (i=1; i<nTime; i++)
    TimeHistData->Time[i] = 0.25*fac*Time[i-1] + (1.-0.5*fac)*Time[i] + 0.25*fac*Time[i+1];

  // set modified time step
  for (i=0; i<nTime; i++){
    TimeHistData->TimeStep[i] = -TimeHistData->Time[i];
    TimeHistData->TimeStep[i] += (i==nTime-1) ? Time[nTime] : TimeHistData->Time[i+1];
  }

  // release time
  xf_Release( (void *) Time);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_OptimizeTimeHistData
static int 
xf_AdaptTimeHistData(xf_TimeHistData *TimeHistData, int *RefIndicatorTime,
                     xf_TimeHistData **pOldTimeHistData) 
{
  int ierr, i, j, jstart, jend;
  int nTime;
  int nTimeNew, iNew;
  int *TimeIsSet = NULL;
  real dN, Nreal, Nnext, fac, Time, dt, EndTime;
  real *dt_desired = NULL;
  real *n_desired = NULL;
  xf_TimeHistData *OldTimeHistData = NULL;
  
  // current number of time steps (slabs)
  nTime = TimeHistData->nTime;

  // end time
  EndTime = TimeHistData->Time[nTime-1] + TimeHistData->TimeStep[nTime-1];

  // Allocate a dt_desired array
  ierr = xf_Error(xf_Alloc( (void **) &dt_desired, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Fill in dt_desired using current dt (TimeHistData) and RefIndicatorTime
  for (i=0; i<nTime; i++){
    fac = 1.0;
    fac = (RefIndicatorTime[i] ==  1) ? 0.5 : fac; // time step grows in refinement 
    fac = (RefIndicatorTime[i] == -1) ? 2.0 : fac; // time step shrinks in coarsening     
    dt_desired[i] = TimeHistData->TimeStep[i] * fac;
  } // i

  // Add up dt/dt_desired to get total number of desired slabs -> Nreal
  for (i=0, Nreal=0.; i<nTime; i++) Nreal += TimeHistData->TimeStep[i]/dt_desired[i];

  // round (via truncation) to nearest integer Nreal->nTimeNew
  nTimeNew = (int) (Nreal + 0.5+MEPS); // MEPS is to make sure we round up, for unit tests

  // make sure we have at least one new time step
  nTimeNew = max(1, nTimeNew);

  // Scale dt_desired by Nreal/nTimeNew
  if (((real) nTimeNew) != Nreal)
    for (i=0; i<nTime; i++) dt_desired[i] *= Nreal/((real) nTimeNew);

  // Create OldTimeHistData that will store original time history
  ierr = xf_Error(xf_CreateTimeHistData(pOldTimeHistData));
  if (ierr != xf_OK) return ierr;
  OldTimeHistData = (*pOldTimeHistData);

  // Store away the original time history
  OldTimeHistData->nTime      = nTime;
  OldTimeHistData->Time       = TimeHistData->Time;
  OldTimeHistData->TimeStep   = TimeHistData->TimeStep;
  OldTimeHistData->TimeScheme = TimeHistData->TimeScheme;
    
  // Allocate new time history (nTimeNew slabs)
  TimeHistData->nTime = nTimeNew;
  ierr = xf_Error(xf_Alloc( (void **) &TimeHistData->Time, nTimeNew, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &TimeHistData->TimeStep, nTimeNew, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &TimeHistData->TimeScheme, nTimeNew,
                           sizeof(enum xfe_TimeSchemeType)));
  if (ierr != xf_OK) return ierr;

  // also reallocate space for outputs
  if (TimeHistData->nOutput > 0){
    ierr = xf_Error(xf_ReAllocCopy2( (void ***) &TimeHistData->OutputValues, 
                                    TimeHistData->nOutput, nTime, 
                                    TimeHistData->nOutput, nTimeNew, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  // convert dt_desired to n_desired = dt/dt_desired
  for (i=0; i<nTime; i++) dt_desired[i] = OldTimeHistData->TimeStep[i]/dt_desired[i];
  n_desired = dt_desired;


  // Allocate at TimeIsSet flag, and set it to zero
  ierr = xf_Error(xf_Alloc( (void **) &TimeIsSet, nTimeNew, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nTimeNew; i++) TimeIsSet[i] = 0;

  // set initial time
  Time = TimeHistData->Time[0] = OldTimeHistData->Time[0];
  TimeIsSet[0] = 1; // flag as set

  // Set times of new time history
  Nreal = 0.;
  for (i=0; i<nTime; i++){ // loop over old time slabs
    dt = OldTimeHistData->TimeStep[i]; // old time step
    Nnext = Nreal + n_desired[i]; // Nnext is a real number
    
    jstart = Nreal; // truncated
    if (jstart < Nreal) jstart++; // jstart is now >= Nreal
    jend   = Nnext; // truncated
    if (jend >= nTimeNew) jend = nTimeNew-1;
    
    //xf_printf("i=%d, Nreal=%.5f, Nnext=%.5f, jstart=%d, jend=%d\n", 
    //          i, Nreal, Nnext, jstart, jend);

    // sanity check
    if ((jstart>0) && (!TimeIsSet[jstart-1])) return xf_Error(xf_CODE_LOGIC_ERROR);

    for (j=jstart; j<=jend; j++){
      if (TimeIsSet[j]) continue;
      // factor into time slab
      fac = (((real) j) - Nreal)/(Nnext-Nreal);
      // set new time
      TimeHistData->Time[j] = OldTimeHistData->Time[i] + fac*dt;
      TimeIsSet[j] = 1;
      //xf_printf("  Just created a time node (j=%d) at Time = %.5f\n", j, TimeHistData->Time[j]); 
    }
    Nreal = Nnext;
  } // i

  // check that all times have been set
  for (i=0; i<nTimeNew; i++) 
    if (!TimeIsSet[i]) return xf_Error(xf_CODE_LOGIC_ERROR);

  // set time steps of new time history
  for (i=0; i<(nTimeNew-1); i++) 
    TimeHistData->TimeStep[i] = TimeHistData->Time[i+1] - TimeHistData->Time[i];
  TimeHistData->TimeStep[nTimeNew-1] = EndTime - TimeHistData->Time[nTimeNew-1];



/*   // Set times and time steps of new time history */
/*   Time = TimeHistData->Time[0] = OldTimeHistData->Time[0]; // set initial time */
/*   Nreal = 0.; */
/*   for (i=0, iNew=1; i<nTime; i++){ // loop over old time slabs */
/*     dt = OldTimeHistData->TimeStep[i]; // old time step */
/*     Nnext = Nreal + n_desired[i]; // Nnext is a real number */
/*     xf_printf("i = %d, iNew=%d, dt = %.5f, Nreal = %.5f, n_desired[i] = %.5f\n",  */
/*               i, iNew, dt, Nreal, n_desired[i]); */
/*     if ((Nnext < iNew) && (i==(nTime-1)) ){  */
/*       // machine precision may prevent us from exactly reaching the final time ... correct this */
/*       dN = ((real) iNew) - Nnext; */
/*       xf_printf("Correcting for machine precision in AdaptTimeHistData, dN = %.10E\n", dN); */
/*       if (dN > 1e-10) return xf_Error(xf_CODE_LOGIC_ERROR); */
/*       Nnext = (real) iNew; */
/*     } */
/*     if (Nnext >= iNew){ */
/*       while (iNew <= (Nnext+MEPS)){ */
/*         Time += (iNew - Nreal)/n_desired[i] * dt; */
/*         TimeHistData->TimeStep[iNew-1] = Time - TimeHistData->Time[iNew-1]; */
/*         if (iNew < nTimeNew) TimeHistData->Time[iNew] = Time; */
/*         xf_printf("  iNew=%d, Time=%.5f, Nreal=%.5f, Nnext=%.5f\n", iNew, Time, Nreal, Nnext); */
/*         xf_printf("    Nnext = %.15E\n", Nnext); */
/*         Nreal = (real) iNew; */
/*         iNew++; */
/*       } */
/*       if (Nnext > Nreal){ */
/*         Time += (Nnext-Nreal)/n_desired[i]*dt; */
/*         /\* if (iNew == nTimeNew){ // set last time step *\/ */
/* /\*           TimeHistData->TimeStep[iNew-1] = Time - TimeHistData->Time[iNew-1]; *\/ */
/* /\*           xf_printf("  +iNew=%d, Time=%.5f, Nreal=%.5f, Nnext=%.5f\n", iNew, Time, Nreal, Nnext); *\/ */
/* /\*           iNew++; *\/ */
/* /\*         } *\/ */
/*       } */
/*     } */
/*     else  */
/*       Time += dt; */
/*     Nreal = Nnext; */
/*   } // i */

/*   // sanity check */
/*   if (iNew != (nTimeNew+1)){ */
/*     xf_printf("iNew = %d, nTimeNew+1=%d\n", iNew, nTimeNew+1); */
/*     return xf_Error(xf_CODE_LOGIC_ERROR); */
/*   } */

  // set time scheme (assume constant, equal to first time slab of old time history)
  for (i=0; i<nTimeNew; i++){
    TimeHistData->TimeScheme[i] = OldTimeHistData->TimeScheme[0];
  } // i
  
  // free memory
  xf_Release( (void *) dt_desired);
  xf_Release( (void *) TimeIsSet);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SplitTimeHistData
int 
xf_SplitTimeHistData(xf_TimeHistData *TimeHistData, int *RefIndicatorTime, int iSplit) 
{
  int ierr;
  int i, j, k, nSplit;
  int nTime, nAdd, nTimeNew;
  enum xfe_Bool Uniform;
  real *Time = NULL;
  real *TimeStep = NULL;
  enum xfe_TimeSchemeType *TimeScheme = NULL;
  
  // number of time steps (slabs)
  nTime = TimeHistData->nTime;
  
  // set default flags
  Uniform = xfe_False;

  // number of additional time steps = sum(RefIndicatorTime)
  if (RefIndicatorTime == NULL){
    if (iSplit == -1){ // flag for uniform refinement
      nAdd = nTime;
      Uniform = xfe_True;
    }
    else nAdd = ((iSplit >= 0) && (iSplit < nTime));
  }
  else{
    for (i=0, nAdd=0; i<nTime; i++){
      // this function is not indended for coarsening -- see AdaptTimeHistData instead
      if (RefIndicatorTime[i] < 0) return xf_Error(xf_INPUT_ERROR);
      nAdd += RefIndicatorTime[i];
    }
  }

  if (nAdd == 0){
    xf_printf("Warning, no time steps selected for adaptation in SplitTimeHistData.\n");
    return xf_OK;
  }
  
  // allocate new vectors
  TimeHistData->nTime = nTimeNew = nTime + nAdd;
  ierr = xf_Error(xf_Alloc( (void **) &Time, nTimeNew, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &TimeStep, nTimeNew, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &TimeScheme, nTimeNew,
                           sizeof(enum xfe_TimeSchemeType)));
  if (ierr != xf_OK) return ierr;
  
  // copy over data, accounting for refinement
  for (i=0,j=0; i<nTime; i++){
    nSplit = ((RefIndicatorTime == NULL) ? ((i == iSplit) || Uniform) : RefIndicatorTime[i]);
    for (k=0; k<=nSplit; k++,j++){
      TimeStep[j]   = TimeHistData->TimeStep[i] / ( (real) nSplit+1);
      Time[j]       = TimeHistData->Time[i] + ( (real) k) * TimeStep[j];
      TimeScheme[j] = TimeHistData->TimeScheme[i]; // BDF2 handling corrected below
    } // j
  } // i
  if (j != nTimeNew) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  // special cases
  if ((nTimeNew > 2) && (TimeScheme[1] == xfe_TimeSchemeBDF1)
      && (TimeScheme[2] == xfe_TimeSchemeBDF2)) TimeScheme[1] = xfe_TimeSchemeBDF2;
  
  // reallocate space for outputs
  if (TimeHistData->nOutput > 0){
    ierr = xf_Error(xf_ReAllocCopy2( (void ***) &TimeHistData->OutputValues, 
                                    TimeHistData->nOutput, nTime, 
                                    TimeHistData->nOutput, nTimeNew, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  
  // destroy existing structures
  xf_Release(  (void  *) TimeHistData->Time);
  xf_Release(  (void  *) TimeHistData->TimeStep);
  xf_Release(  (void  *) TimeHistData->TimeScheme);
  
  // point to new structures
  TimeHistData->Time       = Time;
  TimeHistData->TimeStep   = TimeStep;
  TimeHistData->TimeScheme = TimeScheme;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AdaptAllUnsteady
int 
xf_AdaptAllUnsteady(xf_All *All, int iAdapt, const char *SavePrefixNext,
                    xf_TimeHistData *TimeHistData, enum xfe_Bool *pDoneAdapt)
{
  int ierr;
  int i, nTime;
  int *RefIndicatorTime = NULL;
  enum xfe_AdaptOnType AdaptOn;
  enum xfe_AdaptMechType AdaptMechanics;
  char AdaptOutput[xf_MAXSTRLEN];
  char AdaptVariableSet[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  char *s = NULL;
  real OutputError, AdaptTolerance;
  real facSmooth;
  xf_Vector *RefIndicatorSpace;
  xf_Vector *BasisOrderElem = NULL;
  xf_DataSet *DataSet = NULL;
  xf_Output *Output;
  xf_TimeHistData *OldTimeHistData = NULL;
  
  // pull off SavePrefix
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;
  
  /* Write out .xfa file before performing adaptation */
  sprintf(OutputFile, "%s_PostSolve.xfa\0", SavePrefix, iAdapt);
  ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
  if (ierr!=xf_OK) return ierr;
  
  // note, this flag must be defined on input
  if (*pDoneAdapt) return xf_OK;
  
  /* Determine what we are adapting on */
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptOn", 
                                     xfe_AdaptOnName, (int ) xfe_AdaptOnLast, 
                                     (int *) &AdaptOn));
  if (ierr != xf_OK) return ierr;
  
  /* Determine the adaptation mechanics */
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptMechanics", 
                                     xfe_AdaptMechName, (int ) xfe_AdaptMechLast, 
                                     (int *) &AdaptMechanics));
  if (ierr != xf_OK) return ierr;
  
  /* Determine adaptation tolerance */
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "AdaptTolerance", &AdaptTolerance));
  if (ierr != xf_OK) return ierr;
  
  
  // allocate temporal indicator vector
  nTime = TimeHistData->nTime;
  ierr = xf_Error(xf_Alloc( (void **) &RefIndicatorTime, nTime, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  
  if (AdaptOn == xfe_AdaptOnOutput){
    /* Determine which output we are adapting on */
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptOutput", AdaptOutput));
    if (ierr != xf_OK) return ierr;
    
    /* Determine entropy on which we are adapting*/
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptVariableSet", 
                                   AdaptVariableSet));
    if (ierr != xf_OK) return ierr;
    
    if (xf_NotNull(AdaptOutput)){
      // Output error estimate is already computed and stored in Output->ErrEst
      ierr = xf_Error(xf_FindOutput(All->EqnSet, AdaptOutput, &Output));
      if (ierr != xf_OK) return ierr;
      OutputError = Output->ErrEst;
      
      // print out info and exit if adaptation tolerance met
      xf_printf(" Unsteady output error estimate for %s = %.10E\n", AdaptOutput, OutputError);  
      if (fabs(OutputError) <= AdaptTolerance){
        xf_printf("\nAdaptation tolerance (%.10E) met. Exiting.\n", AdaptTolerance);
        
        (*pDoneAdapt) = xfe_True;
      }
    }
    // either name of output or "Output" if adapting on variable set
    s = (xf_NotNull(AdaptVariableSet) ? xfe_AdaptOnName[AdaptOn] : AdaptOutput);
    // Determine time indicator and space indicator for refinement
    ierr = xf_Error(xf_OptimizeRefUnsteadyOneShot(All, TimeHistData, s,
                                                  &RefIndicatorSpace, RefIndicatorTime,
                                                  &BasisOrderElem));
    if (ierr != xf_OK) return ierr;
  }
  else if ((AdaptOn == xfe_AdaptOnInterpol) ||
           (AdaptOn == xfe_AdaptOnResidual)){
    // Adaptation governed by interpolation error or residual
    // Error estimate already computed during forward run 
    //ierr = xf_Error(xf_ApplyTimeSchemeAdapt(All, SavePrefix, AdaptOn, TimeHistData));
    //if (ierr != xf_OK) return ierr;
    // Determine time indicator and space indicator for refinement
    ierr = xf_Error(xf_OptimizeRefUnsteadyOneShot(All, TimeHistData, xfe_AdaptOnName[AdaptOn],
                                                  &RefIndicatorSpace, RefIndicatorTime,
                                                  &BasisOrderElem));
    if (ierr != xf_OK) return ierr;
  }
  else if (AdaptOn == xfe_AdaptOnUniform){
    // Allocate refinement indicator in space
    ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                  NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_True,
                                  xfe_False, NULL, &RefIndicatorSpace, NULL));
    if (ierr != xf_OK) return ierr;
    
    // set indicators to refine all elements in space and all time slabs
    ierr = xf_Error(xf_SetConstVector(RefIndicatorSpace, 1, 0.));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nTime; i++) RefIndicatorTime[i] = 1;
  }
  else return xf_Error(xf_NOT_SUPPORTED);
  
  // return here if done with adapting (met tolerance)
  if (*pDoneAdapt) return xf_OK;
  
  /* Delete non-essential data -- for memory cleanup */
  ierr = xf_Error(xf_DataSetDeleteNonEssential(All->DataSet));
  if (ierr != xf_OK) return ierr;  
  
  /* Reload initial condition for projection purposes (may want to use
   to initialize steady presolve on next adapt iter).*/
  // read .data from file (multiple vectors in one dataset)
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;
  sprintf(Title, "%s_U%d.data\0", SavePrefix, 0);
  ierr = xf_ReadDataSetBinary(All->Mesh, NULL, Title, DataSet);
  if (ierr == xf_NOT_FOUND){
    xf_printf("Not reloading initial condition (not found) prior to adaptation.\n");
  }
  else if (ierr == xf_OK){
    /* Delete any existing states */
    ierr = xf_Error(xf_DataSetRemove(All->DataSet, "State", xfe_True));
    if (ierr != xf_OK) return ierr;
    
    // store U vector in a new data structure with Title="State"
    ierr = xf_Error(xf_DataSetAdd(All->DataSet, "State", xfe_Vector,
                                  xfe_True, (void *) DataSet->Head->Data, NULL));
    if (ierr != xf_OK) return ierr;
    
    DataSet->Head->Data = NULL;
  }
  else return xf_Error(ierr);
  ierr = xf_Error(xf_DestroyDataSet(DataSet)); 
  if (ierr != xf_OK) return ierr;
  
  /* Perform the temporal adaptation */
  ierr = xf_Error(xf_AdaptTimeHistData(TimeHistData, RefIndicatorTime, &OldTimeHistData));
  if (ierr != xf_OK) return ierr;

  // obtain factor for time history smoothing (0=none, 1=full)
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "TimeHistorySmoothFactor", &facSmooth));
  if (ierr != xf_OK) return ierr;

  // smooth new time history, if desired
  if (facSmooth != 0.){
    ierr = xf_Error(xf_SmoothTimeHistData(TimeHistData, facSmooth));
    if (ierr != xf_OK) return ierr;
  }

  
  /* Perform the spatial adaptation */
  if (AdaptMechanics == xfe_AdaptMechHangNode){
    /* Hanging-node adaptation */
    ierr = xf_Error(xf_AdaptHang(All, RefIndicatorSpace));
    if (ierr != xf_OK) return ierr;
    // Release memory
    ierr = xf_Error(xf_DestroyVector(RefIndicatorSpace, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  else if (AdaptMechanics == xfe_AdaptMechOrderRef){
    
    /* Order refinement */
    ierr = xf_Error(xf_AdaptOrderRef(All, SavePrefixNext, RefIndicatorSpace, 
                                     TimeHistData, OldTimeHistData, BasisOrderElem));
    if (ierr != xf_OK) return ierr;
    // Release memory
    ierr = xf_Error(xf_DestroyVector(RefIndicatorSpace, xfe_True));
    if (ierr != xf_OK) return ierr;
    
  }
  else return xf_Error(xf_NOT_SUPPORTED);
  
  
  // Release memory  
  ierr = xf_Error(xf_DestroyVector(BasisOrderElem, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyTimeHistData(OldTimeHistData));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) RefIndicatorTime);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_PreAdaptRobust
static int 
xf_PreAdaptRobust(xf_All *All, xf_KeyValue *pKeyValueOrig)
{
  int ierr, i;
  xf_KeyValue KeyValue;
  char *Key, KeyOrig[xf_MAXSTRLEN], ValueOrig[xf_MAXSTRLEN];
  real AdaptRobustFixedFraction;
  
  KeyValue = All->Param->KeyValue;
  
  ierr = xf_Error(xf_InitKeyValue(pKeyValueOrig));
  if (ierr != xf_OK) return ierr;
  
  /*************************/
  //backup of the parameters
  /*************************/
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptOn", ValueOrig));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_AddKeyValue(pKeyValueOrig, "AdaptOn", ValueOrig, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptIsotropic", ValueOrig));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_AddKeyValue(pKeyValueOrig, "AdaptIsotropic", ValueOrig, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptIncludeP", ValueOrig));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_AddKeyValue(pKeyValueOrig, "AdaptIncludeP", ValueOrig, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptFixedFraction", ValueOrig));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_AddKeyValue(pKeyValueOrig, "AdaptFixedFraction", ValueOrig, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix",
                                 ValueOrig));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_AddKeyValue(pKeyValueOrig, "SavePrefix", 
                                 ValueOrig, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  //not writing the in-between meshes
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, 
                                 "SavePrefix", "None"));
  if (ierr != xf_OK) return ierr;
  
  //setting up the parameters for adapting isotropically
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, 
                                 "AdaptIsotropic", "True"));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, 
                                 "AdaptIncludeP", "False"));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, 
                                     "AdaptRobustFixedFraction", 
                                     &AdaptRobustFixedFraction));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, 
                                     "AdaptFixedFraction",
                                     AdaptRobustFixedFraction));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_PostAdaptRobust
static int 
xf_PostAdaptRobust(xf_All *All, xf_KeyValue *pKeyValueOrig)
{
  int ierr, i;
  char ValueOrig[xf_MAXSTRLEN];
  
  for (i=0; i<(*pKeyValueOrig).nKey; i++){
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, (*pKeyValueOrig).Key[i], 
                                   (*pKeyValueOrig).Value[i]));
    if (ierr != xf_OK) return ierr;
  } //i
  
  ierr = xf_Error(xf_DestroyKeyValue(pKeyValueOrig));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AdaptRobust
int
xf_AdaptRobust(xf_All *All, xf_Vector *U, int *piAdaptRobust, 
               enum xfe_Bool *pAdapted, xf_SolverData *SolverData)
{
  int ierr, AdaptRobustMaxIter, AdaptRobustIter0, sr, egrp, elem, n, nn, s;
  enum xfe_Bool DoneAdapt, AdaptRobust, WriteIntermediate, AdaptIndBdown, found;
  enum xfe_AdaptOnType AdaptRobustInd;
  char filename[xf_MAXSTRLEN], SavePrefix[xf_MAXSTRLEN];
  real AdaptRobustCFLAmpFactor, ResTolerance, ElemVol, P;
  xf_KeyValue KeyValueOrig;
  xf_Vector *dt, *AdaptIndicator, *Gp, *dU, *R, *AdaptIndicatorBdown, *EG;
  xf_JacobianMatrix *R_U;
  xf_Data *D;
  
  (*pAdapted) = xfe_False;
  
  //get parameters
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                     "AdaptRobust", &AdaptRobust));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                     "AdaptRobustWriteInterm", 
                                     &WriteIntermediate));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, 
                                    "AdaptRobustMaxIter", 
                                    &AdaptRobustMaxIter));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, 
                                    "iAdaptRobust", &AdaptRobustIter0));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, 
                                     "AdaptRobustIndicator", 
                                     xfe_AdaptOnName, (int) xfe_AdaptOnLast, 
                                     (int *)&AdaptRobustInd));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue,
                                     "AdaptRobustMaxCFLAmpFactor",
                                     &AdaptRobustCFLAmpFactor));
  if (ierr != xf_OK) return ierr;
  
  /* Breakdown of the Adaptive indicator sum */
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                     "AdaptRobustIndBdown", 
                                     &AdaptIndBdown));
  if (ierr != xf_OK) return ierr;
  
  sr = All->EqnSet->StateRank;
  
  if (AdaptRobust == xfe_True && (*piAdaptRobust) < 
      AdaptRobustMaxIter+AdaptRobustIter0){
    
    //backup and change the apropriate parameters
    ierr = xf_Error(xf_PreAdaptRobust(All, &KeyValueOrig));
    if (ierr != xf_OK) return ierr;
    
    xf_printf(" Adapting the mesh to try to converge.\n iAdaptRobust = %d\n",
              (*piAdaptRobust));
    
    ierr = xf_Error(xf_FindAdaptIndicator(All, xfe_False, &AdaptIndicator));
    if (ierr != xf_OK) return ierr;
    AdaptIndicator->SolverRole = xfe_SolverRoleOther;
    
    if (AdaptIndBdown){
      ierr = xf_Error(xf_FindVector(All, "AdaptRobustIndBdown", 
                                    xfe_LinkageGlobElem, sr, U->StateName, 
                                    0, 0, NULL, NULL, NULL, NULL, NULL,
                                    xfe_SizeReal, xfe_False,  xfe_True, &D, 
                                    &AdaptIndicatorBdown, NULL));
      if (ierr != xf_OK) return ierr;
      D->ReadWrite = xfe_True;
      AdaptIndicatorBdown->SolverRole = xfe_SolverRoleOther;
    }
    
    switch (AdaptRobustInd) {
      case xfe_AdaptOnPenalty:
        ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                              xfe_True, NULL, &R_U, NULL));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, 
                                             xfe_True, NULL, &R, &found));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_CalculateResidual(All, U, R, R_U, SolverData);
        if (ierr != xf_OK) return ierr;
        
        //find the search direction dU
        ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_False, 
                                             xfe_True, NULL, &dU, NULL));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue,
                                           "ResidualTolerance",
                                           &ResTolerance));
        if (ierr != xf_OK) return ierr;
        
        //locate the timestep vector
        ierr = xf_Error(xf_FindVector(All, "TimeStep", xfe_LinkageGlobElem, 1, 
                                      NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, 
                                      xfe_SizeReal, xfe_False,  xfe_True, NULL, 
                                      &dt, &found));
        if (ierr != xf_OK) return ierr;
        
        // use maximum CFL achieved and amplify it by a factor
        ierr = xf_CalculateArtificialTimeStep(All, U, SolverData->MaxCFLAchieved*
                                              AdaptRobustCFLAmpFactor, dt);
        if (ierr != xf_OK) return ierr;
        
        //add mass matrix with a larger dt than the maximum previously achieved
        ierr = xf_Error(xf_AddMassMatrix(All, 0, dt, U, R, R_U,NULL));
        if (ierr != xf_OK) return ierr;
        
        SolverData->LinResTol = 1e-6;
        //calculate dU
        ierr = xf_SolveLinearSystem(All, R_U, R, xfe_False, -1, SolverData, dU);
        if (ierr != xf_OK) return ierr;
        
        //find the penalty gradient vector
        ierr = xf_Error(xf_FindSimilarVector(All, U, "GradPenalty", xfe_True, 
                                             xfe_False, NULL, &Gp, NULL));
        if (ierr != xf_OK) return ierr;
        
        //calculate the penalty function gradient
        ierr = xf_Error(xf_CalculatePenaltyGradient(All, U, Gp));
        if (ierr != xf_OK) return ierr;
        
        //compute element-wise dot products.
        for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++){
          for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
            AdaptIndicator->GenArray[egrp].rValue[elem][0] = 0.0;
            nn = U->GenArray[egrp].r/sr;
            for (s = 0; s < sr; s++){
              if (AdaptIndBdown)
                AdaptIndicatorBdown->GenArray[egrp].rValue[elem][s] = 0.0;
              for (n = 0; n < nn; n++){
                AdaptIndicator->GenArray[egrp].rValue[elem][0] += 
                dU->GenArray[egrp].rValue[elem][n*sr+s]*
                Gp->GenArray[egrp].rValue[elem][n*sr+s];
                if (AdaptIndBdown){
                  AdaptIndicatorBdown->GenArray[egrp].rValue[elem][s] += 
                  dU->GenArray[egrp].rValue[elem][n*sr+s]*
                  Gp->GenArray[egrp].rValue[elem][n*sr+s];
                }
              }
            }
            //normalization
            if (AdaptIndBdown){
              for (s = 0; s < sr; s++){
                AdaptIndicatorBdown->GenArray[egrp].rValue[elem][s] /= 
                AdaptIndicator->GenArray[egrp].rValue[elem][0];
              }
            }
          }
        }
        
        break;
      case xfe_AdaptOnResidual:
        //use the absolute values of the overall residual on each cell
        ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, 
                                             xfe_True, NULL, &R, &found));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_CalculateResidual(All, U, R, NULL, SolverData);
        if (ierr != xf_OK) return ierr;
        
        // get element geometry
        ierr = xf_Error(xf_FindElemGeom(All, &EG));
        if (ierr != xf_OK) return ierr;
        
        for (egrp = 0; egrp < All->Mesh->nElemGroup; egrp++){
          for (elem = 0; elem < All->Mesh->ElemGroup[egrp].nElem; elem++){
            ElemVol  = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
            AdaptIndicator->GenArray[egrp].rValue[elem][0] = 0.0;
            nn = U->GenArray[egrp].r/sr;
            for (s = 0; s < sr; s++){
              if (AdaptIndBdown)
                AdaptIndicatorBdown->GenArray[egrp].rValue[elem][s] = 0.0;
              for (n = 0; n < nn; n++){
                AdaptIndicator->GenArray[egrp].rValue[elem][0] += 
                fabs(R->GenArray[egrp].rValue[elem][n*sr+s]/ElemVol);
                if (AdaptIndBdown){
                  AdaptIndicatorBdown->GenArray[egrp].rValue[elem][s] += 
                  fabs(R->GenArray[egrp].rValue[elem][n*sr+s]/ElemVol);
                }
              }
            }
            //normalization
            if (AdaptIndBdown){
              for (s = 0; s < sr; s++){
                AdaptIndicatorBdown->GenArray[egrp].rValue[elem][s] /= 
                AdaptIndicator->GenArray[egrp].rValue[elem][0];
              }
            }
          }
        }
        EG = NULL;
        break;
      default:
        return xf_Error(xf_NOT_SUPPORTED);
        break;
    }
    
    //performing adaptation
    ierr = xf_Error(xf_AdaptAll(All, (*piAdaptRobust), &DoneAdapt, AdaptIndicator));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
    if (ierr != xf_OK) return ierr;
    
    (*piAdaptRobust)++;
    
    if (WriteIntermediate){
      /* If SavePrefix is None or NULL, will not write anything */
      ierr = xf_Error(xf_GetKeyValue(KeyValueOrig, "SavePrefix", SavePrefix));
      if (ierr != xf_OK) return ierr;
      
      if (xf_NotNull(SavePrefix)){
        /* Write out .xfa file*/
        sprintf(filename, "%s_iRA%02d.xfa\0", SavePrefix, (*piAdaptRobust));
        ierr = xf_Error(xf_WriteAllBinary(All, filename));
        if (ierr!=xf_OK) return ierr;
      }
    }
    
    //Copying back the values to All structure
    ierr = xf_Error(xf_PostAdaptRobust(All, &KeyValueOrig));
    if (ierr != xf_OK) return ierr;
    
    if (AdaptIndBdown){
    	AdaptIndicatorBdown->SolverRole = xfe_SolverRoleNone;
    }
    
    (*pAdapted) = xfe_True;
  }
  return xf_OK;
}



#if( UNIT_TEST==1 )
#include "xf_Adapt.test.in"
#endif
