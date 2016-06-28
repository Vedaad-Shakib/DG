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
  FILE:  xf_MeshMotion.c

  This file contains functions for working with the mesh Motion
  structure.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_MeshMotionAnalytical.h"
#include "xf_MeshMotionIO.h"
#include "xf_Math.h"
#include "xf_MeshMotionGCL.h"



/******************************************************************/
//   FUNCTION Definition: xf_InitMeshMotion
static int 
xf_InitMeshMotion( xf_MeshMotion *Motion){

  int ierr;
  
  Motion->Type = xfe_MotionLast;
  Motion->Data = NULL;
  Motion->Active = xfe_True;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateMeshMotion
int 
xf_CreateMeshMotion( xf_MeshMotion **pMotion){

  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pMotion, 1, sizeof(xf_MeshMotion)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_InitMeshMotion((*pMotion)));
  if (ierr != xf_OK) return ierr;
   
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyMeshMotion
int 
xf_DestroyMeshMotion( xf_MeshMotion *Motion){
  
  int i, ierr;

  if (Motion == NULL) return xf_OK;

  switch (Motion->Type){
  case xfe_Motion_Analytical:
    ierr = xf_Error(xf_DestroyAnaMotionsSet((xf_AnaMotionsSet *) Motion->Data));
    if (ierr != xf_OK) return ierr;
    break;

  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  xf_Release( (void *) Motion);
 
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_InitMotionData
void
xf_InitMotionData( xf_MotionData *MData){
  
  MData->npoint      = 0;
  MData->dim         = 0;
  MData->x           = NULL;
  MData->vg          = NULL;
  MData->G           = NULL;
  MData->g           = NULL;
  MData->gb          = NULL;
  MData->gbigb_X     = NULL;
  MData->Ginv        = NULL;
  MData->GCLVector   = NULL;
}
/******************************************************************/
//   FUNCTION Definition: xf_CreateMotionData
int 
xf_CreateMotionData( xf_All *All, xf_MotionData **pMData){
  
  int ierr;
  enum xfe_Bool UseGCL;
  xf_MotionData *MData;
  
  ierr = xf_Error(xf_Alloc( (void **) pMData, 1, sizeof(xf_MotionData)));
  if (ierr != xf_OK) return ierr;

  MData = (*pMData);

  xf_InitMotionData(MData);

  // is GCL on?
  MData->GCLVector = NULL; // assume not, unless ...
  if (All != NULL){
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
    if (ierr != xf_OK) return ierr;
    
    // if GCL is on, find GCL vector (should exist)
    if (UseGCL){
      ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &MData->GCLVector));
      if (ierr != xf_OK) return ierr;
    }
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyMotionData
void 
xf_DestroyMotionData( xf_MotionData *MData){

  if (MData != NULL){
    
    xf_Release( (void *) MData->x);
    xf_Release( (void *) MData->vg);
    xf_Release( (void *) MData->G);
    xf_Release( (void *) MData->g);
    xf_Release( (void *) MData->gb);
    xf_Release( (void *) MData->gbigb_X);
    xf_Release( (void *) MData->Ginv);
    
    xf_Release( (void *) MData);

  }
}

/******************************************************************/
//   FUNCTION Definition: xf_AllocMotionData
int 
xf_AllocMotionData( unsigned int AllocFlag, int npoint, int dim, xf_MotionData *MData)
{
  int ierr;

  if (MData == NULL) return xf_Error(xf_INPUT_ERROR);

  MData->npoint = npoint; 
  MData->dim     = dim;

  // allocations are not independent
  AllocFlag |= xfb_MD_x; // always need this
  if (AllocFlag & xfb_MD_Ginv) AllocFlag |= xfb_MD_G;
  
  if (AllocFlag & xfb_MD_x){
    ierr = xf_Error(xf_ReAlloc((void **) &MData->x,  npoint*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  if (AllocFlag & xfb_MD_vg){
    ierr = xf_Error(xf_ReAlloc((void **) &MData->vg, npoint*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  if (AllocFlag & xfb_MD_G){
    ierr = xf_Error(xf_ReAlloc((void **) &MData->G,  npoint*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  if (AllocFlag & xfb_MD_g){
    ierr = xf_Error(xf_ReAlloc((void **) &MData->g,  npoint, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  if (AllocFlag & xfb_MD_gb){
    ierr = xf_Error(xf_ReAlloc((void **) &MData->gb,  npoint, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  if (AllocFlag & xfb_MD_gbigb_X){
    ierr = xf_Error(xf_ReAlloc((void **) &MData->gbigb_X, npoint*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  if (AllocFlag & xfb_MD_Ginv){
    ierr = xf_Error(xf_ReAlloc((void **) &MData->Ginv,  npoint*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReAllocMotionData
static int 
xf_ReAllocMotionData( enum xfe_Bool AllocAll, xf_MotionData *MData)
{
  int ierr;
  int npoint, dim;
  unsigned int AllocFlag;


  if (MData == NULL) return xf_Error(xf_INPUT_ERROR);

  npoint = MData->npoint;
  dim = MData->dim;

  AllocFlag = 0;
  AllocFlag |= ((AllocAll) || (MData->x       != NULL))*xfb_MD_x;
  AllocFlag |= ((AllocAll) || (MData->vg      != NULL))*xfb_MD_vg;
  AllocFlag |= ((AllocAll) || (MData->G       != NULL))*xfb_MD_G;
  AllocFlag |= ((AllocAll) || (MData->g       != NULL))*xfb_MD_g;
  AllocFlag |= ((AllocAll) || (MData->gb      != NULL))*xfb_MD_gb;
  AllocFlag |= ((AllocAll) || (MData->gbigb_X != NULL))*xfb_MD_gbigb_X;
  AllocFlag |= ((AllocAll) || (MData->Ginv    != NULL))*xfb_MD_Ginv;

  ierr = xf_Error(xf_AllocMotionData(AllocFlag, npoint, dim, MData));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionMap
int 
xf_MeshMotionMap( int egrp, int elem, xf_BasisData *PhiData, xf_MeshMotion *Motion, 
		  int npoint, int dim, real Time, const real *X, xf_MotionData *MData)
{
  int ierr, i;
  enum xfe_Bool FirstAlloc;
  enum xfe_Bool UseGCL;
  real *gbigb_X = NULL, *rtemp;

  // set mapping to identity if motion is not active
  if (!Motion->Active) return xf_Error(xf_NOT_SUPPORTED);

  // Are we using a GCL?
  UseGCL = ((MData->GCLVector != NULL) && (PhiData != NULL));

  // (re)allocate space for data, if need more points or higher dim
  if ((npoint > MData->npoint) || (dim > MData->dim)){
    FirstAlloc = (MData->npoint == 0); // useful behavior: means allocate all
    MData->npoint = max(MData->npoint, npoint);
    MData->dim    = max(MData->dim,    dim);
    ierr = xf_Error(xf_ReAllocMotionData(FirstAlloc, MData));
    if (ierr != xf_OK) return ierr;
  }

  // temporarily set gbigb_X to NULL if GCL is on
  // NEW: let's use the analytical gbigb_X
  //if (UseGCL) swap(MData->gbigb_X, gbigb_X, rtemp);

  // Apply different types of mesh motion, depending on type
  switch (Motion->Type){
  case xfe_Motion_Analytical:
    ierr = xf_Error(xf_MeshMotionMap_Analytical((xf_AnaMotionsSet *) Motion->Data, 
						npoint, dim, Time, X, MData));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  // reset gbigb_X from NULL if GCL is on
  // NEW: let's use the analytical gbigb_X
  //if (UseGCL) swap(MData->gbigb_X, gbigb_X, rtemp);

  // calculate inverse of mapping if desired
  if (MData->Ginv != NULL)
    for (i=0; i<npoint; i++){
      ierr = xf_Error(xf_MatDetInv(MData->G+i*dim*dim, dim, NULL, MData->Ginv+i*dim*dim));
      if (ierr != xf_OK) return ierr;
    } // i

  if (UseGCL){
    // calculate gb and gbigb_X if GCL is on
    ierr = xf_Error(xf_MeshMotionMap_GCL( MData->GCLVector, egrp, elem, PhiData, npoint,
					  dim, MData->gb, NULL));
    if (ierr != xf_OK) return ierr;
  }
  else{
    // if not using GCL, set gb = g
    if ((MData->gb != NULL) && (MData->g != NULL))
      for (i=0; i<npoint; i++) MData->gb[i] = MData->g[i];
  }

  
  return xf_OK;
}






