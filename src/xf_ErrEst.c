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
 FILE:  xf_ErrEst.c
 
 This file contains functions for error estimation.
 
 */


#include "xf_AllStruct.h"
#include "xf_All.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_Math.h"
#include "xf_Memory.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Residual.h"
#include "xf_Solver.h"
#include "xf_AdaptStruct.h"
#include "xf_MPI.h"
#include "xf_MeshTools.h"
#include "xf_MeshToolsStruct.h"
#include "xf_Output.h"
#include "xf_EqnSet.h"
#include "xf_LinearSolver.h"
#include "xf_Basis.h"
#include "xf_AdaptHang.h"
#include "xf_AllPull.h"

/******************************************************************/
//   FUNCTION Definition: xf_PreFineSpaceSolve
int 
xf_PreFineSpaceSolve(xf_All *All, xf_KeyValue *pKeyValueOrig)
{
  int ierr, i;
  xf_KeyValue KeyValue;
  char *Key, KeyOrig[xf_MAXSTRLEN], ValueOrig[xf_MAXSTRLEN];
  real CFLStart;
  
  KeyValue = All->Param->KeyValue;
  
  ierr = xf_Error(xf_InitKeyValue(pKeyValueOrig));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<KeyValue.nKey; i++){
    Key = KeyValue.Key[i];
    if (strncmp(Key, "FineSpace_", 10) == 0){
      sprintf(KeyOrig, "%s\0", Key+10);
      ierr = xf_Error(xf_GetKeyValue(KeyValue, KeyOrig, ValueOrig));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_AddKeyValue(pKeyValueOrig, KeyOrig, ValueOrig, xfe_True));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_SetKeyValue(KeyValue, KeyOrig, KeyValue.Value[i]));
      if (ierr != xf_OK) return ierr;
    }
  } //i
  
  //turning off the robustness adaptation and penalization
  ierr = xf_Error(xf_GetKeyValue(KeyValue, "AdaptRobust", ValueOrig));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_AddKeyValue(pKeyValueOrig, "AdaptRobust", ValueOrig, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SetKeyValue(KeyValue, "AdaptRobust", "False"));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValue(KeyValue, "PenalizeResidual", ValueOrig));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_AddKeyValue(pKeyValueOrig, "PenalizeResidual", ValueOrig, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_SetKeyValue(KeyValue, "PenalizeResidual", "False"));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_PostFineSpaceSolve
int 
xf_PostFineSpaceSolve(xf_All *All, xf_KeyValue *pKeyValueOrig)
{
  int ierr, i;
  char ValueOrig[xf_MAXSTRLEN];
  
  for (i=0; i<(*pKeyValueOrig).nKey; i++){
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, (*pKeyValueOrig).Key[i], 
                                   (*pKeyValueOrig).Value[i]));
    if (ierr != xf_OK) return ierr;
  } //i
  
  ierr = xf_Error(xf_GetKeyValue((*pKeyValueOrig), "AdaptRobust", ValueOrig));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "AdaptRobust", ValueOrig));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValue((*pKeyValueOrig), "PenalizeResidual", ValueOrig));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "PenalizeResidual", ValueOrig));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyKeyValue(pKeyValueOrig));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ZeroVectorBoundingBox
static int 
xf_ZeroVectorBoundingBox(xf_All *All, xf_Vector *R, real *BB)
{
  /*
   
   PURPOSE: 
   
   Zeros out vector R on elements that are outside bounding box BB.
   
   INPUTS:
   
   All : All file
   R : vector to be zeroed out
   BB : bounding box = xmin xmax ymin ymax zmin zmax
   
   OUTPUTS: 
   
   R : modified   
   
   RETURNS: Error Code
   
  */
  int ierr;
  int egrp, elem;
  int i, dim, k, r, d;
  int *Node, nnode;
  int nn, nvec[xf_MAXQ1NODE];
  enum xfe_Bool Inside;
  real *x;
  xf_Mesh *Mesh;
  xf_GenArray *ga;

  Mesh = All->Mesh;
  dim = Mesh->Dim;

  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    ga = R->GenArray+egrp;

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      Node  = Mesh->ElemGroup[egrp].Node[elem];
      nnode = Mesh->ElemGroup[egrp].nNode;

      Inside = xfe_True;
      for (i=0; (i<nnode) && (Inside); i++){
	x = Mesh->Coord[Node[i]];
	for (d=0; d<dim; d++)
	  Inside = ((Inside) && (x[d] <= BB[2*d+1]) && (x[d] >= BB[2*d]));
      } // i

      if (Inside) continue;
      
      // zero out R on this elem
      r = ((ga->vr==NULL) ? ga->r : ga->vr[elem]);
      for (k=0; k<r; k++) ga->rValue[elem][k] = 0.;
    }
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_HRefnSubElem
static int 
xf_HRefnSubElem(enum xfe_ShapeType Shape, int *pnsub)
{ 
  /*
   
   PURPOSE: 
   
   Determines number of sub-elements corresponding to uniform h
   refinement of element of shape Shape.
   
   INPUTS:
   
   Shape : geometrical shape of element in question
   
   OUTPUTS: 
   
   (*pnsub) : number of sub elements
   
   RETURNS: Error Code
   
   */
  
  switch (Shape){
  case xfe_Point:
    (*pnsub) = 0;
    break;
  case xfe_Segment:
    (*pnsub) = 2;
    break;
  case xfe_Quadrilateral:
    (*pnsub) = 4;
    break;
  case xfe_Hexahedron:
    (*pnsub) = 8;
  case xfe_Triangle:
  case xfe_Tetrahedron:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_VectorRankMult
static int 
xf_VectorRankMult(xf_Vector *U, int egrp, int nsub)
{
 /*
    
   PURPOSE: 
  
   Expands storage of an interpolated real vector U to be nsub times
   the usual rank.  Used for storing, on the original mesh, states of
   uniformly-refined elements.
  
   INPUTS:
   
   U    : vector for which storage will be expanded the storage
   egrp : element group (i.e. general array index) to work with
   nsub : new rank will be nsub times the original rank
   
   OUTPUTS: 

   U    : rank is expanded (unless it already is expanded coming in)   
   
   RETURNS: Error Code
   
  */

  int ierr, nn, nelem, j, k;
  int rexpected, rcurrent;
  int *vrnew = NULL;
  enum xfe_Bool VariableOrder;
  xf_GenArray *ga = NULL;

  // number of elements
  nelem = U->GenArray[egrp].n;

  // expected rank
  ierr = xf_Error(xf_Order2nNode(U->Basis[egrp], U->Order[egrp], &nn));
  if (ierr != xf_OK) return ierr;
  rexpected = nn*U->StateRank*nsub;

  // current rank
  rcurrent = U->GenArray[egrp].r;

  if (rexpected == rcurrent) // nothing to do
    return xf_OK;
  else{

    ga = U->GenArray + egrp; // quick pointer to general array of interest

    VariableOrder = (ga->vr != NULL);  // are we dealing with a variable order?

    if (VariableOrder){ // yes, variable order
      
      // create new rank vector, vrnew
      ierr = xf_Error(xf_Alloc( (void **) &vrnew, nelem, sizeof(int)));
      if (ierr != xf_OK) return ierr;

      // set vrnew
      for (j=0; j<nelem; j++) vrnew[j] = ga->vr[j]*nsub;

      // variable rank reallocate
      ierr = xf_Error(xf_VReAllocCopy2( (void ***) &ga->rValue, nelem, ga->vr,
					nelem, vrnew, sizeof(real)));
      if (ierr != xf_OK) return ierr;

      // zero out newly-allocated space
      for (j=0; j<nelem; j++)
	for (k=ga->vr[j]; k<vrnew[j]; k++) ga->rValue[j][k] = 0.;

      // destroy vr
      xf_Release( (void *) ga->vr);

      // set vr = vrnew
      ga->vr = vrnew;

    }
    else{  // constant order on each elem

      // constant rank reallocate
      ierr = xf_Error(xf_ReAllocCopy2( (void ***) &ga->rValue, nelem, rcurrent,
				       nelem, rexpected, sizeof(real)));
      if (ierr != xf_OK) return ierr;

    }
    
    // set constant rank
    ga->r = rexpected;
  }
  
  return xf_OK;

}


/******************************************************************/
//   FUNCTION Definition: xf_TraceRefine
static int 
xf_TraceRefine(xf_All *All, int egrp, int elem, xf_Vector **pTracer)
{
 /*
    
   PURPOSE: 

   Creates a writable tracer (integer vector) that is 1 on (egrp,elem)
   of All, and zero on all other elements.  Then uniformly refines
   (egrp,elem).  The returned tracer corresponds to the refined mesh.
  
   INPUTS:
   
   All        : all structure
   egrp, elem : element to refine
   
   OUTPUTS: 

   (*pTracer) : tracer vector
   
   RETURNS: Error Code
   
  */
  int ierr;
  enum xfe_Bool found;
  xf_Vector *Tracer = NULL;
  xf_Vector *RefInd = NULL;
  xf_Data *D = NULL;

  // Create a writeable tracer vecotr on All
  ierr = xf_Error(xf_FindVector(All, "Tracer", xfe_LinkageGlobElem, 
				1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, 
				xfe_False, xfe_True, &D, &Tracer, NULL));
  if (ierr != xf_OK) return ierr;         
  D->ReadWrite = xfe_True;
            
  // Set tracer to zero
  ierr = xf_Error(xf_SetZeroVector(Tracer));
  if (ierr != xf_OK) return ierr;
            
  // Mark requested element
  Tracer->GenArray[egrp].iValue[egrp][elem] = 1;
            
  // Create a refinement indicator on All
  ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem, 
				1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, 
				xfe_False, xfe_False, NULL, &RefInd, NULL));
  if (ierr != xf_OK) return ierr;
  
  // set refinement indicator to zero
  ierr = xf_Error(xf_SetZeroVector(RefInd));
  if (ierr != xf_OK) return ierr;
  
  // turn off verbosity
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", "Low"));
  if (ierr != xf_OK) return ierr;
  
  // do not write VOrder vector file
  ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "AdaptWriteVOrder", xfe_False));
  if (ierr != xf_OK) return ierr;
            
  // mark egrp,elem for refinement
  RefInd->GenArray[egrp].iValue[elem][0] = 1; // indicates uniform refinement
              
  // perform the adaptation
  ierr = xf_Error(xf_AdaptHang(All, RefInd));
  if (ierr != xf_OK) return ierr;
  
  // destroy refinement indicator
  ierr = xf_Error(xf_DestroyVector(RefInd, xfe_True));
  if (ierr != xf_OK) return ierr;

  // find tracer on refined All -> point to it via pTracer
  found = xfe_False;
  ierr = xf_Error(xf_FindVector(All, "Tracer", xfe_LinkageGlobElem, 
				1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, 
				xfe_False, xfe_True, NULL, pTracer, &found));
  if (ierr != xf_OK) return ierr;
  if (found != xfe_True) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  return xf_OK;

}

/******************************************************************/
//   FUNCTION Definition: xf_HRefCalculateResidual
static int 
xf_HRefCalculateResidual(enum xfe_Bool HRef, int OrderIncrement, xf_All *All, 
                         xf_Vector *U, xf_Vector *R, xf_Vector *V)
{
  /*
   
   PURPOSE: 
   
   Wrapper for CalculateResidual; performs residual calculation on a
   uniformly-refined mesh if HRef is True.  Data for all sub-elements
   of an original element is stored associated with the original
   element in a rank-increased array as part of R.
   
   INPUTS:
   
   HRef : True means perform calculation on uniformly refined mesh.
   OrderIncrement : fine space order increment (to account for p-dep of residual)
   All  : All file
   U    : state
   R    : residual (pre-allocated to standard size regardless of HRef)
   V    : auxiliary vector to inject onto the refined space
   
   OUTPUTS: 
   
   R    : filled-in and possibly rank-increased residual vector R
   
   RETURNS: Error Code
   
  */
  int ierr;
  int isub, nsub;
  int egrp, elem, ie, k, r0;
  char StateTitle[] = "HRefCalculateResidual_State";
  char AuxTitle[] = "HRefCalculateResidual_Aux";
  enum xfe_ShapeType Shape;
  xf_SolverData *SolverData = NULL;
  xf_Data *D  = NULL;
  xf_Data *DV = NULL;
  xf_Data *DSmall   = NULL;
  xf_Data *DWallDist = NULL;
  xf_Vector *USmall = NULL;
  xf_Vector *VSmall = NULL;
  xf_Vector *Tracer = NULL;
  xf_Vector *RSmall = NULL;
  xf_GenArray *gaR = NULL;
  xf_All *AllSmall = NULL;

  // standard residual calculation (all we need if HRef == False)
  
  // create/allocate SolverData on All
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;
  
  // account for possible p-dependence of residual
  SolverData->ResidualOrderIncrement = -OrderIncrement;

  // Calculate residual on All
  ierr = xf_Error(xf_CalculateResidual(All, U, R, NULL, SolverData));
  if (ierr != xf_OK) return ierr;

  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;
  
  if (!HRef) return xf_OK; // done

  // Add input U as a writeable data piece in All->Data
  ierr = xf_Error(xf_DataSetAdd(All->DataSet, StateTitle, xfe_Vector, 
				xfe_True, (void *) U, &D));
  if (ierr != xf_OK) return ierr;

  if (V != NULL){  // same with V if it is not NULL
    ierr = xf_Error(xf_DataSetAdd(All->DataSet, AuxTitle, xfe_Vector, 
				  xfe_True, (void *) V, &DV));
    if (ierr != xf_OK) return ierr;
  }

  // make certain auxiliary vectors writeable so they get transferred to AllSmall
  ierr = xf_FindDataByTitle(All->DataSet, "WallDistance", xfe_Vector, &DWallDist);
  if (ierr == xf_OK) DWallDist->ReadWrite = xfe_True;
  else if (ierr != xf_NOT_FOUND) return xf_Error(ierr);
  

  // Loop over element groups
  for (egrp=0; egrp<All->Mesh->nElemGroup; egrp++){

    // How many sub-elements per original element do we get in this element group?
    ierr = xf_Error(xf_Basis2Shape(All->Mesh->ElemGroup[egrp].QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_HRefnSubElem(Shape, &nsub));
    if (ierr != xf_OK) return ierr;

    // Expand storage in R based on number of sub-elements, if necessary
    ierr = xf_Error(xf_VectorRankMult(R, egrp, nsub));
    if (ierr != xf_OK) return ierr;

    // same with V if not NULL
    if (V != NULL){
      ierr = xf_Error(xf_VectorRankMult(V, egrp, nsub));
      if (ierr != xf_OK) return ierr;
    }

    // Loop over elements
    for (elem=0; elem<All->Mesh->ElemGroup[egrp].nElem; elem++){

      // Pull small xfa
      ierr = xf_Error(xf_CreateAll(&AllSmall, xfe_False));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_AllPull(All, AllSmall, egrp, elem));
      if (ierr != xf_OK) return ierr;
      
      // Create an integer tracer vector and uniform refine central elem on small mesh
      ierr = xf_Error(xf_TraceRefine(AllSmall, egrp, 0, &Tracer));
      if (ierr != xf_OK) return ierr;
      
      // find USmall = state on small (via data name)
      ierr = xf_Error(xf_FindDataByTitle(AllSmall->DataSet, StateTitle, xfe_Vector, &DSmall));
      if (ierr != xf_OK) return ierr;
      USmall = (xf_Vector *) DSmall->Data;
      
      // same with V if not NULL
      if (V != NULL){
	ierr = xf_Error(xf_FindDataByTitle(AllSmall->DataSet, AuxTitle, xfe_Vector, &DSmall));
	if (ierr != xf_OK) return ierr;
	VSmall = (xf_Vector *) DSmall->Data;
      }
      
      // find residual vector, RSmall, that is similar to USmall
      ierr = xf_Error(xf_FindSimilarVector(AllSmall, USmall, "SmallResidual", xfe_False, 
					   xfe_True, NULL, &RSmall, NULL));
      if (ierr != xf_OK) return ierr;
      
      // create/allocate SolverData on AllSmall
      ierr = xf_Error(xf_CreateSolverData(AllSmall, &SolverData));
      if (ierr != xf_OK) return ierr;
      
      // account for possible p-dependence of residual
      SolverData->ResidualOrderIncrement = -OrderIncrement;

      // Calculate residual on small (use local SolverData)
      ierr = xf_Error(xf_CalculateResidual(AllSmall, USmall, RSmall, NULL, SolverData));
      if (ierr != xf_OK) return ierr;
      
      // destroy SolverData
      ierr = xf_Error(xf_DestroySolverData(SolverData));
      if (ierr != xf_OK) return ierr;

      // Loop over (sub) elements of AllSmall; if tracer is set, set R from RSmall
      for (ie=0, isub=0; ie<AllSmall->Mesh->ElemGroup[egrp].nElem; ie++){
	if (Tracer->GenArray[egrp].iValue[ie][0] == 1){
	  gaR = RSmall->GenArray + egrp;
	  r0 = ((gaR->vr == NULL) ? gaR->r : gaR->vr[ie]);
	  for (k=0; k<r0; k++)
	    R->GenArray[egrp].rValue[elem][isub*r0 + k] = gaR->rValue[ie][k];
	  if (V != NULL)
	    for (k=0; k<r0; k++)
	      V->GenArray[egrp].rValue[elem][isub*r0 + k] = VSmall->GenArray[egrp].rValue[ie][k];
	  isub++;
	}
      } // ie
      
      // Error if number of sub-elements is not consistent with the pre-determined value for this egrp
      if (isub != nsub) return xf_Error( xf_CODE_LOGIC_ERROR);

      // Set AllSmall->EqnSet to NULL and destroy AllSmall
      AllSmall->EqnSet = NULL;
      ierr = xf_Error(xf_DestroyAll(AllSmall));
      if (ierr != xf_OK) return ierr;

    } // elem

  } // egrp

  // Destroy newly-added data set in All corresponding to U ... do not destroy U though
  D->Data = NULL;
  ierr = xf_Error(xf_DataSetRemove(All->DataSet, StateTitle, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  if (V != NULL){  // same with V if it is not NULL
    DV->Data = NULL;
    ierr = xf_Error(xf_DataSetRemove(All->DataSet, AuxTitle, xfe_False));
    if (ierr != xf_OK) return ierr;
  }

  // make certain auxiliary vectors not-writeable
  ierr = xf_FindDataByTitle(All->DataSet, "WallDistance", xfe_Vector, &DWallDist);
  if (ierr == xf_OK) DWallDist->ReadWrite = xfe_False;
  else if (ierr != xf_NOT_FOUND) return xf_Error(ierr);

  return xf_OK;

}


/******************************************************************/
//   FUNCTION Definition: xf_LocalizeError
static int 
xf_LocalizeError(enum xfe_Bool HRef, xf_All *All, xf_Vector *Rh, 
		 xf_Vector *AdjFine, xf_Vector *ElemIndicator, 
		 enum xfe_Bool ElemIndSign, real *pOutputError)
{
  /*
   
   PURPOSE: 
   
   Computes output error estimate and localizes contributions to elements.
   
   INPUTS:
   
   HRef : True means perform calculation on uniformly refined mesh
   All : All file
   Rh : fine space residual vector, Rh(UH)
   AdjFine : fine space adjoint solution
   
   OUTPUTS: 
   
   ElemIndicator : elemental error indicator (optional)
   (*pOutputError) : output error (optional)
   
   
   RETURNS: Error Code
   
   */
  int ierr;
  int sr, r, k, i;
  int egrp, elem;
  int isub, nsub;
  int r0, i0;
  enum xfe_Verbosity Verbosity;
  enum xfe_Bool VolSpecificRes;
  enum xfe_Bool UseBoundingBox;
  char BoundingBoxString[xf_MAXSTRLEN];
  enum xfe_ShapeType Shape;
  real BoundingBox[6];
  real *evec, ElemVol;
  real SumIndicator;
  real *EI, *EAdj, *ER;
  xf_Mesh *Mesh;
  xf_Vector *EG;
  xf_GenArray *ga;
  
  // Determine verbosity for printing purposes
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", xfe_VerbosityName, 
                                     (int ) xfe_VerbosityLast, (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // Determine if bounding box is used to restrict err est to a certain region
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "ErrEstBoundingBox", BoundingBoxString));
  if (ierr != xf_OK) return ierr;
  UseBoundingBox = xfe_False;
  if (xf_NotNull(BoundingBoxString)){
    // Bounding box is on; zero out residuals on elements outside box
    UseBoundingBox = xfe_True;
    ierr = xf_Error(xf_ScanReal(BoundingBoxString, 2*All->Mesh->Dim, BoundingBox));
    if (ierr != xf_OK) return ierr;
  }

  
  if (pOutputError != NULL){
    if (AdjFine != NULL){
      // dJ = - Psih^T * Rh
      ierr = xf_Error(xf_VectorDot(Rh, AdjFine, pOutputError));
      if (ierr != xf_OK) return ierr;
      (*pOutputError) *=-1;
      
      if (Verbosity != xfe_VerbosityLow)
        xf_printf("%45s %.10E\n", "Output error estimate (JH - Jh) =", (*pOutputError));

      if (UseBoundingBox){
	// zero out R and repeat output calculation
	ierr = xf_Error(xf_ZeroVectorBoundingBox(All, Rh, BoundingBox));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VectorDot(Rh, AdjFine, pOutputError));
	if (ierr != xf_OK) return ierr;
	(*pOutputError) *=-1;
	if (Verbosity != xfe_VerbosityLow)
	  xf_printf("%45s %.10E\n", "Bounding-box output error estimate (JH - Jh) =", (*pOutputError));
      }
    }
    else{
      ierr = xf_Error(xf_VectorNorm(Rh, 1, pOutputError));
      if (ierr != xf_OK) return ierr;
      
      if (Verbosity != xfe_VerbosityLow)
        xf_printf("%45s %.10E\n", "Residual norm estimate =", (*pOutputError));
    }
  }
  
  if (ElemIndicator == NULL) return xf_OK; // no indicator is requested
  
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;
  
  ierr = xf_Error(xf_Alloc( (void **) &evec, sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  if (ierr != xf_OK) return ierr;
  
  // loop over elements and calculate elemental contributions to error
  SumIndicator = 0;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    if (HRef){
      // How many sub-elements per original element do we get in this element group?
      ierr = xf_Error(xf_Basis2Shape(All->Mesh->ElemGroup[egrp].QBasis, &Shape));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_HRefnSubElem(Shape, &nsub));
      if (ierr != xf_OK) return ierr;
    }
    else nsub = 1;
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      EI = ElemIndicator->GenArray[egrp].rValue[elem];
      SumIndicator -= EI[0];
      ga = Rh->GenArray+egrp; 
      r0 = ((ga->vr==NULL) ? ga->r : ga->vr[elem]); // possibly-augmented rank
      for (k=0; k<sr; k++) evec[k] = 0.;
      ER = Rh->GenArray[egrp].rValue[elem]; // residual
      
      if ((r0%nsub) != 0) return xf_Error(xf_INPUT_ERROR);
      r = r0/nsub; // original rank (e.g. on each subelement)

      for (isub=0; isub<nsub; isub++){ // loop over sub-elements

	i0 = isub*r; // starting index into adjoint and residual

	if (AdjFine != NULL){ // adjoint-weighted residual
	  EAdj = AdjFine->GenArray[egrp].rValue[elem]; // adjoint
	  for (i=0; i<r; i++) evec[i%sr] += EAdj[i0+i]*ER[i0+i]; // each equation separately
	}
	else // residual only (need fabs; otherwise sum is 0 due to Galerkin orthog.)
	  for (i=0; i<r; i++){
	    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
					       "VolumeSpecificResidual", 
					       &VolSpecificRes));
	    if (ierr != xf_OK) return ierr;
	    if (!VolSpecificRes)
	      evec[i%sr] += fabs(ER[i0+i]); // each equation separately
	    else {
	      ElemVol  = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
	      evec[i%sr] += fabs(ER[i0+i])/ElemVol;
	    } 
	  }
      
	if (ElemIndSign == xfe_True)
	  for (k=0; k<sr; k++) EI[0] += evec[k]; // sum signed values over equations
	else
	  for (k=0; k<sr; k++) EI[0] += fabs(evec[k]); // sum abs values over equations
	
	SumIndicator += EI[0];

      } // isub

    } // elem
  } // egrp
  
  if ((pOutputError != NULL) && (Verbosity != xfe_VerbosityLow)){
    // reduce-sum
    ierr = xf_Error(xf_MPI_Allreduce(&SumIndicator, 1, xfe_SizeReal, xfe_MPI_SUM));
    if (ierr != xf_OK) return ierr;
    xf_printf("%45s %.10E\n", "Sum(|element output error indicator|) =", SumIndicator);
  }
  
  EG = NULL;
  xf_Release(evec);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ErrEstOutput
int 
xf_ErrEstOutput(xf_All *All, char *OutputName, char *VariableSet, 
                enum xfe_Bool StoreFineSpace, enum xfe_Bool ReuseFineSpace,
                xf_Vector *ElemIndicator, enum xfe_Bool ElemIndSign, 
                real *pOutputError)
{
  int ierr, i;
  int nPsi;
  int OrderIncrement;
  //int *OrderVec = NULL;
  char FineTitlePrimal[] = "FinePrimal_h";
  char FineTitleDual[]   = "FineDual_h";
  char ValueOrig[xf_MAXSTRLEN];
  enum xfe_Bool Found, EntropyFlag = xfe_False, ResidualFlag = xfe_False;
  enum xfe_Bool UsePrimal = xfe_False, UseDual = xfe_False, Combined = xfe_False;
  enum xfe_Bool HRef = xfe_False;
  enum xfe_Verbosity Verbosity;
  real OutputError, fac;
  xf_KeyValue KeyValueOrig;
  xf_Output *Output;
  xf_Vector **Psi;
  xf_Vector *U, *Adj, *AdjRh = NULL, *UH = NULL;
  xf_Vector *UFine = NULL, *AdjFine = NULL, *AdjCoarse = NULL, *Rh;
  xf_Vector **pAdjRh = NULL;
  xf_SolverData *SolverData = NULL;
  xf_Data *D;

  int iRes, nRes, nResTerm;
  char ResTermName[xf_MAXSTRLEN];
  enum xfe_ResTermType ResType;
  xf_Vector *Rh2 = NULL;
  xf_ResTerm ResTerm, *pResTerm = NULL;
  
  //Options for just removing the old fine-space solutions
  if ((StoreFineSpace == xfe_False) && (ReuseFineSpace == xfe_True) 
      && (ElemIndicator == NULL)){
    ierr = xf_Error(xf_DataSetRemove(All->DataSet, FineTitlePrimal, xfe_False));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DataSetRemove(All->DataSet, FineTitleDual, xfe_False));
    if (ierr != xf_OK) return ierr;
    
    return xf_OK;
  }
  
  // Determine verbosity for printing purposes
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", xfe_VerbosityName, 
                                     (int ) xfe_VerbosityLast, (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  if ((OutputName == NULL) && (VariableSet == NULL)){
    ResidualFlag = xfe_True; // adapting on residual only
    if (Verbosity != xfe_VerbosityLow)
      xf_printf("Using fine-space residual to form error indicator.\n");
  }
  else{
    ResidualFlag = xfe_False;
    // both variables should be passed in
    if ((OutputName == NULL) || (VariableSet == NULL)) return xf_Error(xf_INPUT_ERROR);
    
    // Only one of OutputName or VariableSet is relevant
    if ((xf_NotNull(OutputName)) && (xf_NotNull(VariableSet))) return xf_Error(xf_INPUT_ERROR);

    Combined = xfe_False;
    
    if (xf_NotNull(OutputName)){
      ierr = xf_Error(xf_FindOutput(All->EqnSet, OutputName, &Output));
      if (ierr != xf_OK) return ierr;
      if (Output->nSumOutput > 0) Combined = xfe_True;
      
      if (Verbosity != xfe_VerbosityLow)
        xf_printf("Estimating error in output = %s\n", OutputName);
      EntropyFlag = xfe_False;
    }
    else if (xf_NotNull(VariableSet)){
      if (Verbosity != xfe_VerbosityLow)
        xf_printf("Estimating error in entropy flux integral output for %s\n", VariableSet);
      EntropyFlag = xfe_True;
    }
    else return xf_Error(xf_INPUT_ERROR);
  }
  
  // Determine if we're using the primal and/or dual error estimate form
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ErrEstUsePrimal", &UsePrimal));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ErrEstUseDual", &UseDual));
  if (ierr != xf_OK) return ierr;
  
  if (ResidualFlag) UseDual = xfe_False; // dual residual does not exist in this case
  
  if ((!UsePrimal) && (!UseDual)) return xf_Error(xf_INPUT_ERROR); // must use at least one
  
  /* Pull off order increment for error estimation */
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement", &OrderIncrement));
  if (ierr != xf_OK) return ierr;

  /* Is error estimate computed on a (uniformly) h-refined mesh? */
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ErrEstHRef", &HRef));
  if (ierr != xf_OK) return ierr;
  if (HRef) xf_printf("Using h-refinement in error estimate.\n");

  
  /* Pull off state and adjoint for the adapted output.  Error if
   required adjoint does not exist. */
  ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D, NULL));
  if (ierr != xf_OK) return ierr;
  U = (xf_Vector *) D->Data;
  
  if ((!EntropyFlag) && (!ResidualFlag)){
    ierr = xf_Error(xf_FindAdjointVectors(All, U, OutputName, xfe_False,
                                          xfe_False, &nPsi, &Psi, &Found));
    if (ierr != xf_OK) return ierr;
    if ((!Found) || (nPsi != 1)) return xf_Error(xf_NOT_FOUND);
    Adj = Psi[0];
    xf_Release( (void *) Psi);
  }
  
  // If it is a combined output, only the primal state should be used
  if (ReuseFineSpace && !Combined){
    
    /*-------------------------------------*/
    /* Locate UFine and AdjFine if reusing */
    /*-------------------------------------*/
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, FineTitlePrimal, xfe_Vector, &D));
    if (ierr != xf_OK) return ierr;
    UFine = (xf_Vector *) D->Data;
    
    if (UseDual) return xf_Error(xf_NOT_SUPPORTED); // need to add code logic for AdjRh calculation
    
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, FineTitleDual, xfe_Vector, &D));
    if (ierr != xf_OK) return ierr;
    AdjFine = (xf_Vector *) D->Data;
    
  }
  else{
    /*-------------------------------------------------*/
    /* Projection of state and adjoint into fine space */
    /*-------------------------------------------------*/ 
    
    if (Combined){
      ierr = xf_Error(xf_FindDataByTitle(All->DataSet, FineTitlePrimal, xfe_Vector, &D));
      if (ierr != xf_OK) return ierr;
      UFine = (xf_Vector *) D->Data;
    }
    else {// Not combined
     
      if (Verbosity != xfe_VerbosityLow)
        xf_printf("Using fine-space order increment of %d\n", OrderIncrement);
      
      // Copy and then project the state and adjoint into a finer space      
      ierr = xf_Error(xf_CreateVector(&UFine));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_CopyVector(All->Mesh, U, UFine);
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, UFine, UFine->Basis, 
							     xfe_BasisLast, OrderIncrement));
      if (ierr != xf_OK) return ierr;
    }
    
    
    if ((!EntropyFlag) && (!ResidualFlag)){
      ierr = xf_Error(xf_CreateVector(&AdjFine));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_CopyVector(All->Mesh, Adj, AdjFine);
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, AdjFine, AdjFine->Basis, 
							     xfe_BasisLast, OrderIncrement));
      if (ierr != xf_OK) return ierr;
    }
    
    
    /*---------------------*/
    /* Fine space solution */
    /*---------------------*/
    
    if (!ResidualFlag){ // not necessary if just using residual
      
      if (!Combined){
        // Set desired solver flags for the fine space solve
        ierr = xf_Error(xf_PreFineSpaceSolve(All, &KeyValueOrig));
        if (ierr != xf_OK) return ierr;
        
        // Solve primal on fine space
        if (Verbosity != xfe_VerbosityLow)
          xf_printf("\nIterating fine-space primal problem.\n");
        
        ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &UFine));
        if (ierr != xf_OK) return ierr;
      }
      // for debugging, dump out fine U solution
      /*     ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "Vector", UFine, "UFine.data")); */
      /*     if (ierr != xf_OK) return ierr; */
      
      if (!EntropyFlag){
        
        if (UseDual){
          // locate adjoint residual vector
          /* Do not need to create a separate vector; SolveAdjoints
           returns a pointer to this (same as adjoint residual) */
          pAdjRh = &AdjRh;
        }
        else pAdjRh = NULL;
        
        
        // Solve adjoint on fine space
        if (Verbosity != xfe_VerbosityLow)
          xf_printf("\nIterating fine-space adjoint problem.\n");
        /* note, if AdjRh != NULL is passed in, it will be returned with 
         the adjoint residual */
        ierr = xf_Error(xf_SolveAdjoints(All, 0.0, 1.0, xfe_False, UFine, 1, 
                                         NULL, &AdjFine, pAdjRh, xfe_False, 
                                         xfe_True));
        if (ierr != xf_OK) return ierr;
        
        // subtract off coarse adjoint
        ierr = xf_Error(xf_CreateVector(&AdjCoarse));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_CopyVector(All->Mesh, AdjFine, AdjCoarse);
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ProjectVector(All, Adj, xfe_False, AdjCoarse));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_SetVector(AdjCoarse, xfe_Sub, AdjFine));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_DestroyVector(AdjCoarse, xfe_True));
        if (ierr != xf_OK) return ierr;
        
      }
      else{
        // Set AdjFine to Entropy variables corresponding to UFine
        ierr = xf_Error(xf_CreateVector(&AdjFine));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_CopyVector(All->Mesh, UFine, AdjFine);
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ChangeVariableSet(All, VariableSet, AdjFine));
        if (ierr != xf_OK) return ierr;
        
        // subtract off coarse entropy
        ierr = xf_Error(xf_CreateVector(&AdjCoarse));
        if (ierr != xf_OK) return ierr;
  
        ierr = xf_CopyVector(All->Mesh, AdjFine, AdjCoarse);
        if (ierr != xf_OK) return ierr;
  
        ierr = xf_Error(xf_ProjectVector(All, U, xfe_False, AdjCoarse));
        if (ierr != xf_OK) return ierr;
  
        ierr = xf_Error(xf_ChangeVariableSet(All, VariableSet, AdjCoarse));
        if (ierr != xf_OK) return ierr;
  
        ierr = xf_Error(xf_SetVector(AdjCoarse, xfe_Sub, AdjFine));
        if (ierr != xf_OK) return ierr;
  
        ierr = xf_Error(xf_DestroyVector(AdjCoarse, xfe_True));
        if (ierr != xf_OK) return ierr;
      }      
      
      if (!Combined){
        // Reset solver flags to original
        ierr = xf_Error(xf_PostFineSpaceSolve(All, &KeyValueOrig));
        if (ierr != xf_OK) return ierr;
      }
      
    } //
    
  } // end else have to recalculate fine space quantities
  
  /*------------------------------------------*/
  /* Fine space weighted residual calculation */
  /*------------------------------------------*/
  
  // initialize output error to 0
  if (pOutputError != NULL) (*pOutputError) = 0.0;
  
  // both flags true indicate a 50/50 split
  fac = ((UsePrimal && UseDual) ? 0.5 : 1.0);
  
  // UH = coarse solution on fine space
  ierr = xf_Error(xf_CreateVector(&UH));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_CopyVector(All->Mesh, UFine, UH);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_ProjectVector(All, U, xfe_False, UH));
  if (ierr != xf_OK) return ierr;
  
  // set ElemIndicator to zero
  if (ElemIndicator != NULL){
    ierr = xf_Error(xf_SetZeroVector(ElemIndicator));
    if (ierr != xf_OK) return ierr;
  }
  
  // Primal residual error estimate
  if (UsePrimal){
    
    // Calculate fine-space residual
    ierr = xf_Error(xf_FindSimilarVector(All, UH, "Residual", xfe_False, xfe_True, NULL, &Rh, NULL));
    if (ierr != xf_OK) return ierr;
        
    // Calculate fine-space residual Rh(UH)
    if (Verbosity != xfe_VerbosityLow)
      xf_printf("\nCalculating fine-space primal residual.\n");
    ierr = xf_Error(xf_HRefCalculateResidual(HRef, OrderIncrement, All, UH, Rh, AdjFine));
    if (ierr != xf_OK) return ierr;
    
    // Error localization
    ierr = xf_Error(xf_LocalizeError(HRef, All, Rh, AdjFine, ElemIndicator, ElemIndSign,
                                     ((pOutputError == NULL) ? NULL : &OutputError)));
    if (ierr != xf_OK) return ierr;


    // is a residual-component breakdown of the error estimate requested?
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "ErrEstResBreakdown", ResTermName));
    if (ierr != xf_OK) return ierr;

    if (xf_NotNull(ResTermName)){
      
      xf_printf("Calculating error estimates from individual %s residual components.\n", ResTermName);
      
      // extra fine-space residual
      ierr = xf_Error(xf_FindSimilarVector(All, Rh, "Residual2", xfe_False, xfe_False, NULL, &Rh2, NULL));
      if (ierr != xf_OK) return ierr;
      
      // current number of residual terms
      nResTerm = All->EqnSet->ResTerms->nResTerm;

      // pull off residual type (enumerated)
      ierr = xf_Error(xf_Value2Enum(ResTermName, xfe_ResTermName, xfe_ResTermLast, (int *) &ResType));
      if( ierr != xf_OK ) return ierr;

      // count how many residual terms are of the requested type
      for (iRes=0, nRes=0; iRes<nResTerm; iRes++)
	nRes += (All->EqnSet->ResTerms->ResTerm[iRes].Type == ResType);
      All->EqnSet->ResTerms->nResTerm = nRes;
      
      // store current residual terms for safe-keeping
      pResTerm = All->EqnSet->ResTerms->ResTerm;

      // Temporarily allocate new residual terms
      ierr = xf_Error(xf_Alloc( (void **) &All->EqnSet->ResTerms->ResTerm, nRes, sizeof(xf_ResTerm)));
      if (ierr != xf_OK) return ierr;
      for (iRes=0, nRes=0; iRes<nResTerm; iRes++)
	if (pResTerm[iRes].Type == ResType){
	  ierr = xf_Error(xf_InitKeyValue(&All->EqnSet->ResTerms->ResTerm[nRes].KeyValue));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_CopyResTerm(pResTerm+iRes, All->EqnSet->ResTerms->ResTerm+nRes));
	  if (ierr != xf_OK) return ierr;
	  nRes += 1;
	}
      if (nRes != All->EqnSet->ResTerms->nResTerm) return xf_Error(xf_CODE_LOGIC_ERROR);
            
      // calculate the required residuals
      ierr = xf_Error(xf_HRefCalculateResidual(HRef, OrderIncrement, All, UH, Rh, NULL)); 
      if (ierr != xf_OK) return ierr;    

      ierr = xf_Error(xf_HRefCalculateResidual(HRef, OrderIncrement, All, UFine, Rh2, NULL)); 
      if (ierr != xf_OK) return ierr;

      // use residual difference for error estimate
      ierr = xf_Error(xf_SetVector(Rh2, xfe_Sub, Rh));
      if (ierr != xf_OK) return ierr;
      
      // error estimate and localization
      OutputError = 0.;       // make sure the output error is reinitialized to zero
      ierr = xf_Error(xf_SetZeroVector(ElemIndicator)); // same with indicator
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_LocalizeError(HRef, All, Rh, AdjFine, ElemIndicator, ElemIndSign,
				       &OutputError));
      if (ierr != xf_OK) return ierr;
      
      xf_printf(" %s-residual term error estimate (%d residual terms): %.10E\n",
		xfe_ResTermName[ResType], nRes, OutputError);
      
      // destroy EqnSet residual terms
      for (i=0; i<All->EqnSet->ResTerms->nResTerm; i++){
	ierr = xf_Error(xf_DestroyResTerm(All->EqnSet->ResTerms->ResTerm+i));
	if (ierr != xf_OK) return ierr;
      }
      xf_Release((void *) All->EqnSet->ResTerms->ResTerm);

      // reset residual terms in All->EqnSet
      All->EqnSet->ResTerms->ResTerm = pResTerm;
      All->EqnSet->ResTerms->nResTerm = nResTerm;
      
      ierr = xf_Error(xf_DestroyVector(Rh2, xfe_True));
      if (ierr != xf_OK) return ierr;  
    }

    if (pOutputError != NULL) (*pOutputError) += fac*OutputError;    

  }
  
  // Dual residual error estimate
  if (UseDual){
    
    if (ReuseFineSpace) return xf_Error(xf_NOT_SUPPORTED); // due to AdjRh calculation
    
    if (HRef) return xf_Error(xf_NOT_SUPPORTED); // need calcs on refined space

    if (AdjRh == NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    // Subtract UH from UFine
    ierr = xf_Error(xf_SetVector(UH, xfe_Sub, UFine));
    if (ierr != xf_OK) return ierr;
    
    if (Verbosity != xfe_VerbosityLow)
      xf_printf("\nLocalizing dual residual.\n");
    
    // Error localization using adjoint residual (computed above)
    ierr = xf_Error(xf_LocalizeError(HRef, All, AdjRh, UFine, ElemIndicator, ElemIndSign, 
                                     ((pOutputError == NULL) ? NULL : &OutputError)));
    if (ierr != xf_OK) return ierr;
    
    if (pOutputError != NULL) (*pOutputError) += fac*OutputError;

    if (xf_NotNull(ResTermName)){
      xf_printf("Dual error estimate with individual residuals is not supported.\n");
      return xf_Error(xf_NOT_SUPPORTED);
    }

  }
  
  
  // multiply indicator by fac if not 1.0
  if ((ElemIndicator != NULL) && (fac != 1.0)){
    ierr = xf_Error(xf_VectorMult(ElemIndicator, fac));
    if (ierr != xf_OK) return ierr;
  }
    
  // Release vectors or store them for future use
  ierr = xf_Error(xf_DestroyVector(UH, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Destroy residual because it might have a corrupted rank (when HRef == True)
  ierr = xf_Error(xf_DataSetRemove(All->DataSet, "Residual", xfe_False));
  if (ierr != xf_OK) return ierr;
  
  if (StoreFineSpace){
    if (!ReuseFineSpace){ // add only when first time (when not reusing)
      // point to UFine and AdjFine from All->DataSet (make writeable)
      if (!Combined){
        ierr = xf_Error(xf_DataSetAdd(All->DataSet, FineTitlePrimal, xfe_Vector, xfe_True, 
                                      (void *) UFine, NULL));
        if (ierr != xf_OK) return ierr;
      }
      if (AdjFine != NULL){
        ierr = xf_Error(xf_DataSetAdd(All->DataSet, FineTitleDual, xfe_Vector, xfe_True, 
                                      (void *) AdjFine, NULL));
        if (ierr != xf_OK) return ierr;
      }
      else return xf_Error(xf_INPUT_ERROR); // should not be storing finespace w/o AdjFine
    }
  }
  else{
    // destroy if not storing
    if (ReuseFineSpace || Combined){
      // destroy in All->DataSet
      ierr = xf_Error(xf_DataSetRemove(All->DataSet, FineTitlePrimal, xfe_False));
      if (ierr != xf_OK) return ierr;
      if (AdjFine != NULL){
        ierr = xf_Error(xf_DataSetRemove(All->DataSet, FineTitleDual, xfe_False));
        if (ierr != xf_OK) return ierr;
      }
      else return xf_Error(xf_INPUT_ERROR); // should not be reusing finespace w/o AdjFine
    }
    else{
      // destroy locally
      ierr = xf_Error(xf_DestroyVector(UFine, xfe_True));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_DestroyVector(AdjFine, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    if (Combined){
      // destroy locally
      ierr = xf_Error(xf_DestroyVector(AdjFine, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SeparateErrEst

int 
xf_SeparateErrEst(xf_All *All, const char *CombinedOutputName, 
                  xf_Vector *Eh, xf_Vector *UH_h)
{
  int ierr, iOutput, myRank, nProc;
  enum xfe_Bool VariableWeights;
  real SumRelativeErrors;
  xf_Output *CombinedOutput, *Output;
  xf_Vector *JprimeH_h;
  char FineTitlePrimal[] = "FinePrimal_h";
  xf_Data *D;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindOutput(All->EqnSet, CombinedOutputName, &CombinedOutput));
	if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "VariableWeights", &VariableWeights));
  if (ierr != xf_OK) return ierr;
  
  if (CombinedOutput->nSumOutput <=1) return xf_INPUT_ERROR;
  
  //create the sensitivity vector
  ierr = xf_Error(xf_CreateVector(&JprimeH_h));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_CopyVector(All->Mesh, UH_h, JprimeH_h);
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0)
    SumRelativeErrors = 0.0;
  
  for (iOutput = 0; iOutput < CombinedOutput->nSumOutput; iOutput++){
    ierr = xf_Error(xf_FindOutput(All->EqnSet, CombinedOutput->SumOutputNames[iOutput], &Output));
    if (ierr != xf_OK) return ierr;
    //calculate output and its derivative
    ierr = xf_Error(xf_CalculateOutput(All, CombinedOutput->SumOutputNames[iOutput], 
                                       UH_h, &(Output->Value), JprimeH_h, xfe_Set));
    if (ierr != xf_OK) return ierr;
    
    //Error estimation
    ierr = xf_Error(xf_VectorDot(Eh, JprimeH_h, &(Output->ErrEst)));
    if (ierr != xf_OK) return ierr;
    
    SumRelativeErrors += fabs(Output->ErrEst)/CombinedOutput->SumOutputErrTols[iOutput];
    
    if (!VariableWeights){
      CombinedOutput->SumOutputWeights[iOutput] = sign(Output->ErrEst)*
      fabs(CombinedOutput->SumOutputWeights[iOutput]);
    }
  } 
  
  ierr = xf_Error(xf_MPI_Bcast(&SumRelativeErrors, sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  
  for (iOutput = 0; iOutput < CombinedOutput->nSumOutput; iOutput++){
    ierr = xf_Error(xf_FindOutput(All->EqnSet, CombinedOutput->SumOutputNames[iOutput], &Output));
    if (ierr != xf_OK) return ierr;
    
    if (VariableWeights){
      //Computing weights
      CombinedOutput->SumOutputWeights[iOutput] = sign(Output->ErrEst)*(fabs(Output->ErrEst)/
                                                                        CombinedOutput->SumOutputErrTols[iOutput])/SumRelativeErrors;
    }
    xf_printf("Error Estimate for %s Output \t = %1.10e \t Weight \t = %1.3f\n", 
              CombinedOutput->SumOutputNames[iOutput], Output->ErrEst,
              CombinedOutput->SumOutputWeights[iOutput]);
  }
  
  ierr = xf_Error(xf_DestroyVector(JprimeH_h, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ErrEstSolution

int 
xf_ErrEstSolution(xf_All *All, xf_Vector *U, const char *CombinedOutputName, 
                  enum xfe_Bool StoreFineSpace, enum xfe_Bool ReuseFineSpace)
{
  int ierr, i, OrderIncrement;
  // int *OrderVec;
  char FineTitlePrimal[] = "FinePrimal_h";
  char FineTitleDual[] = "FineDual_h"; //not used here but created for output error estimation
  enum xfe_Verbosity Verbosity;
  enum xfe_Bool Found;
  real LinTol;
  xf_SolverData *SolverData;
  xf_JacobianMatrix *R_Uh;
  xf_Vector *Uh, *Eh, *Rh;
  xf_Data *D;
  xf_KeyValue KeyValueOrig;
  
  // Determine verbosity for printing purposes
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", xfe_VerbosityName, 
                                     (int ) xfe_VerbosityLast, (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  /* Order increment for error estimation */
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement", &OrderIncrement));
  if (ierr != xf_OK) return ierr;

  
  /*   ierr = xf_Error(xf_Alloc( (void **) &OrderVec, All->Mesh->nElemGroup, sizeof(int))); */
  /*   if (ierr != xf_OK) return ierr; */
  
  if (ReuseFineSpace){
    /*-------------------------------------*/
    /* Locate Uh if reusing */
    /*-------------------------------------*/
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, FineTitlePrimal, xfe_Vector, &D));
    if (ierr != xf_OK) return ierr;
    Uh = (xf_Vector *) D->Data;
    
    //for (i=0; i<All->Mesh->nElemGroup; i++) OrderVec[i] = Uh->Order[i];
  }
  else{
    /*-------------------------------------------------*/
    /* Projection of state into fine space */
    /*-------------------------------------------------*/
        
    if (Verbosity != xfe_VerbosityLow)
      xf_printf("Using fine-space order increment of %d\n", OrderIncrement);
    
    //for (i=0; i<All->Mesh->nElemGroup; i++) OrderVec[i] = U->Order[i] + OrderIncrement;
    
    // Copy and then project the state into a finer space
    
    ierr = xf_Error(xf_CreateVector(&Uh));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_CopyVector(All->Mesh, U, Uh));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Uh, Uh->Basis, 
							   xfe_BasisLast, OrderIncrement));
    if (ierr != xf_OK) return ierr;
    
  }
  // Solve sensitivity on fine space
  if (Verbosity != xfe_VerbosityLow)
    xf_printf("\nComputing solution sensitivity.\n");
  
  //backup of the linear solver tolerance
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "MinLinResDecreaseFactor", &LinTol));
  if (ierr != xf_OK) return ierr;
  //Solve the linear system to machine precision
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "MinLinResDecreaseFactor", 1e-12));
  if (ierr != xf_OK) return ierr;
  /****************************************************************************/
  // Solve for the approximate discretization error Eh
  /****************************************************************************/
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;
  
  // account for possible p-dependence of residual
  SolverData->ResidualOrderIncrement = -OrderIncrement;
  
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, Uh, NULL,
                                        xfe_True, NULL, &R_Uh, &Found));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, Uh, "Residual", xfe_False, xfe_True, NULL, &Rh, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_CalculateResidual(All, Uh, Rh, R_Uh, SolverData);
  if (ierr != xf_OK) return ierr;
  
  // Reset ResidualOrderIncrement
  SolverData->ResidualOrderIncrement = 0;
  

  ierr = xf_Error(xf_FindSimilarVector(All, Uh, "dU", xfe_True, xfe_True, NULL, &Eh, &Found));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SetZeroVector(Eh));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_SolveLinearSystem(All, R_Uh, Rh, xfe_False, -1, SolverData, Eh);
  if (ierr != xf_OK) return ierr; 
  //Setting the linear tolerance back to original value
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "MinLinResDecreaseFactor", LinTol));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;
  
  //Cleaning up and storing what is requested
  if (StoreFineSpace){
    if (!ReuseFineSpace){ // add only when first time (when not reusing)
      // point to UFine and AdjFine from All->DataSet (make writeable)
      ierr = xf_FindDataByTitle(All->DataSet, FineTitlePrimal, xfe_Vector, &D);
      if (ierr == xf_NOT_FOUND){
        ierr = xf_Error(xf_DataSetAdd(All->DataSet, FineTitlePrimal, xfe_Vector, xfe_True, 
                                      (void *) Uh, NULL));
        if (ierr != xf_OK) return ierr;
      }
      else {
        //delete old and add new
        ierr = xf_Error(xf_DataSetRemove(All->DataSet, FineTitlePrimal, xfe_False));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_DataSetAdd(All->DataSet, FineTitlePrimal, xfe_Vector, xfe_True, 
                                      (void *) Uh, NULL));
        if (ierr != xf_OK) return ierr;
      }
    }
  }
  else{
    // destroy if not storing
    if (ReuseFineSpace){
      // destroy in All->DataSet
      ierr = xf_Error(xf_DataSetRemove(All->DataSet, FineTitlePrimal, xfe_False));
      if (ierr != xf_OK) return ierr;
    }
    else{
      // destroy locally
      ierr = xf_Error(xf_DestroyVector(Uh, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  //Separating the error estimations for each output
  ierr = xf_Error(xf_SeparateErrEst(All, CombinedOutputName, Eh, Uh));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DataSetRemove(All->DataSet, "dU", xfe_False));
  if (ierr != xf_OK) return ierr;
  
  //xf_Release((void *)OrderVec);
  
  return xf_OK;
}


#if( UNIT_TEST==1 )
#include "xf_ErrEst.test.in"
#endif
