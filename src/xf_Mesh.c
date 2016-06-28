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
  FILE:  xf_Mesh.c

  This file contains functions for working with the Mesh data structure.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_MeshTools.h"
#include "xf_Quad.h"
#include "xf_MeshMotion.h"
#include "xf_MeshMotionIO.h"
#include "xf_Math.h"


/* Hash structure for reading ASCII grid files */
struct xf_FaceHash
{
  int nVisit;          // number of times this face has been visited
  int BFlag;           // boundary flag (True if this is a boundary face)
  int Group;           // BFGroup or ElemGroup # associated with this face
  int Elem;            // For IFaces, elem # at first visit
  int Face;            // Either BFace # or local face # at first visit
  int nfnode;          // number of nodes on this face
  int *svec;           // vector of nodes on this face, sorted
  struct xf_FaceHash *Next;   // pointer to next structure
};

typedef struct xf_FaceHash xf_FaceHash;


/******************************************************************/
//   FUNCTION Definition: xf_InitMesh
static void
xf_InitMesh( xf_Mesh *Mesh)
{
  // Initializes all data to 0 or NULL
  Mesh->Dim   = 0;
  Mesh->nNode = 0;
  Mesh->Coord = NULL;
  Mesh->nIFace = 0;
  Mesh->IFace  = NULL;
  Mesh->nBFaceGroup = 0;  
  Mesh->BFaceGroup  = NULL;
  Mesh->nElemGroup = 0;
  Mesh->ElemGroup  = NULL;
  Mesh->nPeriodicGroup = 0;
  Mesh->PeriodicGroup = NULL;
  Mesh->ParallelInfo = NULL;
  Mesh->BackgroundMesh = NULL;
  Mesh->Motion = NULL;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateMesh
int 
xf_CreateMesh( xf_Mesh **Mesh){

  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) Mesh, 1, sizeof(xf_Mesh)));
  if (ierr != xf_OK) return ierr;

  xf_InitMesh(*Mesh);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateCutFaceData
static int 
xf_CreateCutFaceData( xf_CutFaceData **pCutFaceData){

  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pCutFaceData, 1, sizeof(xf_CutFaceData)));
  if (ierr != xf_OK) return ierr;

  (*pCutFaceData)->QuadData = NULL;
  (*pCutFaceData)->OrigFace.Group = 0;
  (*pCutFaceData)->OrigFace.Number = 0;
  (*pCutFaceData)->GeomIndex = 0;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyCutFaceData
static int 
xf_DestroyCutFaceData( xf_CutFaceData *CutFaceData){
  if (CutFaceData == NULL) return xf_OK;

  return xf_Error(xf_DestroyQuadData(CutFaceData->QuadData));

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateCutElemData
static int 
xf_CreateCutElemData( int nelem, xf_CutElemData **pCutElemData){

  int ierr, elem;
  
  ierr = xf_Error(xf_Alloc((void **) pCutElemData, nelem, sizeof(xf_CutElemData)));
  if (ierr != xf_OK) return ierr;

  for (elem=0; elem<nelem; elem++){
    (*pCutElemData)[elem].QuadData = NULL;
    (*pCutElemData)[elem].QBasis = 0;
    (*pCutElemData)[elem].QCoord = NULL;
    (*pCutElemData)[elem].GeomIndex = 0;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyCutElemData
static int 
xf_DestroyCutElemData( xf_CutElemData *CutElemData){
  if (CutElemData == NULL) return xf_OK;

  return xf_Error(xf_DestroyQuadData(CutElemData->QuadData));

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyBFaceGroup
int 
xf_DestroyBFaceGroup( xf_BFaceGroup *BFaceGroup){
  int ierr, bface;

  for (bface=0; bface<BFaceGroup->nBFace; bface++){
    if (BFaceGroup->BFace[bface].CutFaceData != NULL){
      ierr = xf_Error(xf_DestroyCutFaceData(BFaceGroup->BFace[bface].CutFaceData));
      if (ierr != xf_OK) return ierr;
    }
  } // bface

  xf_Release((void *)BFaceGroup->Title);
  xf_Release((void *)BFaceGroup->BFace);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyElemGroup
static int 
xf_DestroyElemGroup( xf_ElemGroup *ElemGroup){
  int ierr, elem;

  if (ElemGroup->CutElemData != NULL){
    for (elem=0; elem<ElemGroup->nElem; elem++){
      ierr = xf_Error(xf_DestroyCutElemData(ElemGroup->CutElemData + elem));
      if (ierr != xf_OK) return ierr;
    } // elem
    xf_Release((void *)ElemGroup->CutElemData);
  }

  xf_Release( (void  *)ElemGroup->nFace);
  xf_Release2((void **)ElemGroup->Node);
  xf_Release2((void **)ElemGroup->Face);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateMeshParallelInfo
static int 
xf_CreateMeshParallelInfo( xf_MeshParallelInfo **pParallelInfo){
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pParallelInfo, 1, sizeof(xf_MeshParallelInfo)));
  if (ierr != xf_OK) return ierr;

  (*pParallelInfo)->nIFaceRegular = -1;
  (*pParallelInfo)->ElemLoc2Glob  = NULL;
  (*pParallelInfo)->IFaceLoc2Glob = NULL;
  (*pParallelInfo)->NodeLoc2Glob  = NULL;
  (*pParallelInfo)->nSendElem     = NULL;
  (*pParallelInfo)->SendElem      = NULL;
  (*pParallelInfo)->nRecvElem     = NULL;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyMeshParallelInfo
static int 
xf_DestroyMeshParallelInfo( xf_Mesh *Mesh, xf_MeshParallelInfo *ParallelInfo){
  int egrp;

  for (egrp=0; egrp<2*Mesh->nElemGroup; egrp++) // halos have this allocated too
    xf_Release( (void  *) ParallelInfo->ElemLoc2Glob[egrp]);

  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    xf_Release2((void **) ParallelInfo->SendElem[egrp]);

  xf_Release ( (void  *) ParallelInfo->ElemLoc2Glob);
  xf_Release ( (void  *) ParallelInfo->IFaceLoc2Glob);
  xf_Release ( (void  *) ParallelInfo->NodeLoc2Glob);
  xf_Release2( (void **) ParallelInfo->nSendElem);
  xf_Release ( (void  *) ParallelInfo->SendElem);
  xf_Release2( (void **) ParallelInfo->nRecvElem);  

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyPeriodicGroup
static int 
xf_DestroyPeriodicGroup( xf_PeriodicGroup *PeriodicGroup){
  
  xf_Release( (void  *) PeriodicGroup->PeriodicNode);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyMesh
int 
xf_DestroyMesh( xf_Mesh *Mesh){
  
  int ierr, i, negrp;

  if (Mesh == NULL) return xf_OK;

  xf_Release2((void **)Mesh->Coord);

  for (i=0; i<Mesh->nIFace; i++){
    if (Mesh->IFace[i].CutFaceData != NULL){
      ierr = xf_Error(xf_DestroyCutFaceData(Mesh->IFace[i].CutFaceData));
      if (ierr != xf_OK) return ierr;
    }
  }
  xf_Release((void *)Mesh->IFace);
  
  for (i=0; i<Mesh->nBFaceGroup; i++){
    ierr = xf_Error(xf_DestroyBFaceGroup(Mesh->BFaceGroup + i));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *)Mesh->BFaceGroup);

  // in parallel runs the number of element groups is twice nElemGroup
  negrp = Mesh->nElemGroup;
  if (Mesh->ParallelInfo != NULL) negrp *= 2;

  for (i=0; i<negrp; i++){
    ierr = xf_Error(xf_DestroyElemGroup(Mesh->ElemGroup + i));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *)Mesh->ElemGroup);

  for (i=0; i<Mesh->nPeriodicGroup; i++){
    ierr = xf_Error(xf_DestroyPeriodicGroup(Mesh->PeriodicGroup + i));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *)Mesh->PeriodicGroup);

  if (Mesh->ParallelInfo != NULL){
    ierr = xf_Error(xf_DestroyMeshParallelInfo(Mesh, Mesh->ParallelInfo));
    if (ierr != xf_OK) return ierr;
    xf_Release((void *)Mesh->ParallelInfo);
  }

  // destroy background mesh too
  if (Mesh->BackgroundMesh != NULL){
    ierr = xf_Error(xf_DestroyMesh((xf_Mesh *) Mesh->BackgroundMesh));
    if (ierr != xf_OK) return ierr;
  }

  // destroy mesh motion structure
  if (Mesh->Motion != NULL){
    ierr = xf_Error(xf_DestroyMeshMotion(Mesh->Motion));
    if (ierr != xf_OK) return ierr;
  }

  xf_Release((void *)Mesh);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_InitElemGroup
void
xf_InitElemGroup( xf_ElemGroup *ElemGroup){

  ElemGroup->CutFlag = xfe_False;
  ElemGroup->QBasis = xfe_BasisLast;
  ElemGroup->QOrder = 0;
  ElemGroup->nElem = 0;
  ElemGroup->nFace = NULL;
  ElemGroup->Face = NULL;
  ElemGroup->nNode = 0;
  ElemGroup->Node = NULL;
  ElemGroup->CutElemData = NULL;
}

/******************************************************************/
//   FUNCTION Definition: xf_InitIFace
void
xf_InitIFace( xf_IFace *IFace){

  IFace->ElemGroupL = 0;
  IFace->ElemGroupR = 0;
  IFace->ElemL = 0;
  IFace->ElemR = 0;  
  IFace->FaceL = 0;
  IFace->FaceR = 0;
  IFace->HangNumber = 0;
  IFace->OrientL = -1;
  IFace->OrientR = -1;
  IFace->CutFaceData = NULL;
}

/******************************************************************/
//   FUNCTION Definition: xf_InitBFace
void
xf_InitBFace( xf_BFace *BFace){

  BFace->ElemGroup = 0;
  BFace->Elem = 0;
  BFace->Face = 0;
  BFace->Orient = -1;
  BFace->CutFaceData = NULL;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateNodeHash
static int 
xf_CreateNodeHash(int nNode, xf_FaceHash ***pNode2Hash){
  int ierr, i;

  ierr = xf_Error(xf_Alloc( (void **) pNode2Hash, nNode, sizeof(xf_FaceHash *)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<nNode; i++) (*pNode2Hash)[i] = NULL;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AddFaceToHash
static int 
xf_AddFaceToHash(xf_FaceHash **Node2Hash, int nfnode, int *fvec, 
		 enum xfe_Bool BFlag, int Group, int Elem, int Face, 
		 xf_FaceHash **phinfo, enum xfe_Bool *Exists){
  int ierr, k, n0, *svec;
  xf_FaceHash *hinfo, **prev;

  if (nfnode <= 0) return xf_Error(xf_OUT_OF_BOUNDS);

  // Allocate svec, copy over data
  ierr = xf_Error(xf_Alloc((void **) &svec, nfnode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<nfnode; k++) svec[k] = fvec[k];

  // sort nodes in svec
  ierr = xf_Error(xf_SortInt(nfnode, svec));
  if (ierr != xf_OK) return ierr;

  // check if nodes already exist in hash
  n0 = svec[0];

  if (Node2Hash[n0] == NULL){
    (*Exists) = xfe_False;
    prev = Node2Hash + n0;
  }
  else{
    hinfo = Node2Hash[n0];
    do {
      (*Exists) = xfe_True;
      if (nfnode != hinfo->nfnode) (*Exists) = xfe_False;
      else{
	for (k=0; k<nfnode; k++)
	  if (hinfo->svec[k] != svec[k]){
	    (*Exists) = xfe_False;
	    break;
	  }
      }
      if (*Exists){
	hinfo->nVisit++;
	(*phinfo) = hinfo;
	xf_Release(svec);
	return xf_OK;
      }
      prev = &(hinfo->Next);
      hinfo = hinfo->Next;
    } while (hinfo != NULL);
  }
  

  ierr = xf_Error(xf_Alloc((void **) &hinfo, 1, sizeof(xf_FaceHash)));
  if (ierr != xf_OK) return ierr;
  hinfo->nVisit = 0;
  hinfo->BFlag = BFlag;
  hinfo->Group = Group;
  hinfo->Elem  = Elem;
  hinfo->Face  = Face;
  hinfo->nfnode = nfnode;
  ierr = xf_Error(xf_Alloc((void **) &hinfo->svec, nfnode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<nfnode; k++) hinfo->svec[k] = svec[k];
  hinfo->Next = NULL;

  (*prev)   = hinfo;
  (*phinfo) = hinfo;

  xf_Release(svec);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DeleteFaceFromHash
static int 
xf_DeleteFaceFromHash(xf_FaceHash **Node2Hash, int nfnode, int *fvec){

  int ierr, k, n0, *svec;
  xf_FaceHash *hinfo, *prev;
  enum xfe_Bool found;

  if (nfnode <= 0) return xf_Error(xf_OUT_OF_BOUNDS);

  // Allocate svec, copy over data
  ierr = xf_Error(xf_Alloc((void **) &svec, nfnode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<nfnode; k++) svec[k] = fvec[k];

  // sort nodes in fvec
  ierr = xf_Error(xf_SortInt(nfnode, svec));
  if (ierr != xf_OK) return ierr;

  // check if nodes already exist in hash
  n0 = svec[0];

  if (Node2Hash[n0] == NULL) return xf_Error(xf_NOT_FOUND);

  hinfo = Node2Hash[n0];
  prev = NULL;
  do{
    found = xfe_True;
    if (nfnode != hinfo->nfnode) found = xfe_False;
    else{
      for (k=0; k<nfnode; k++)
	if (hinfo->svec[k] != svec[k]){
	  found = xfe_False;
	  break;
	}
    }

    if (found){
      if (prev != NULL) 
	prev->Next = hinfo->Next;
      else
	Node2Hash[n0] = Node2Hash[n0]->Next;
      xf_Release((void *) hinfo->svec);
      xf_Release((void *) hinfo);
      xf_Release(svec);
      return xf_OK;
    }
    prev = hinfo;
    hinfo = hinfo->Next;
  } while (hinfo != NULL);

  xf_Release(svec);
  return xf_Error(xf_NOT_FOUND);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_HashIsEmpty
static enum xfe_Bool
xf_HashIsEmpty(xf_FaceHash **Node2Hash, int nNode)
{
  int i;

  // check if hash contains any data
  for (i=0; i<nNode; i++) if (Node2Hash[i] != NULL) return xfe_False;

  // hash is empty
  return xfe_True;
}

/******************************************************************/
//   FUNCTION Definition: xf_AddInteriorFace
void
xf_AddInteriorFace(xf_Mesh *Mesh, int egL, int elemL, int faceL, 
		   int egR, int elemR, int faceR)
{
  /* Adds an interior face to Mesh, with the left element info in
     hinfo, and the right element info in (egrp,elem,face).  Note:
     assumes Mesh->IFace has been adequately allocated! */
  xf_InitIFace(Mesh->IFace+Mesh->nIFace);
  Mesh->IFace[Mesh->nIFace].ElemGroupL = egL;
  Mesh->IFace[Mesh->nIFace].ElemL      = elemL;
  Mesh->IFace[Mesh->nIFace].FaceL      = faceL;
  Mesh->IFace[Mesh->nIFace].ElemGroupR = egR;
  Mesh->IFace[Mesh->nIFace].ElemR      = elemR;
  Mesh->IFace[Mesh->nIFace].FaceR      = faceR;
  Mesh->ElemGroup[egL].Face[elemL][faceL].Group  = xf_INTERIORFACE;
  Mesh->ElemGroup[egL].Face[elemL][faceL].Number = Mesh->nIFace;
  Mesh->ElemGroup[egR].Face[elemR][faceR].Group  = xf_INTERIORFACE;
  Mesh->ElemGroup[egR].Face[elemR][faceR].Number = Mesh->nIFace;
  Mesh->nIFace++;
}


/******************************************************************/
//   FUNCTION Definition: xf_AddHangingFace
static int
xf_AddHangingFace(xf_Mesh *Mesh, xf_FaceHash **Node2Hash, 
		  enum xfe_ShapeType Shape, int nnode, int *fvec)
{
  /* This is not an optimized function.  Primary mode of operation for
     non-conforming faces is to generate them during adaptation, not
     to read them in from a text file.  This is for testing purposes
     only.*/
  int ierr, i, k, egrp0, elem0, face0, nface0;
  int nfnode, pos1, pos2, nface, *nFaceNew;
  int iface1, iface2;
  int f0[2], f1[2], f2[2], fvec1[2], fvec2[2];
  int nvec[xf_MAXQ1FACENODE];
  enum xfe_Bool exists;
  xf_FaceHash *hinfo0, *hinfo1, *hinfo2;

  // Supporting only 2d now
  if (Shape != xfe_Segment) return xf_Error(xf_NOT_SUPPORTED);

  f0[0] = fvec[0]; f0[1] = fvec[2];  // coarse face
  f1[0] = fvec[0]; f1[1] = fvec[1];  // fine face 1
  f2[0] = fvec[1]; f2[1] = fvec[2];  // fine face 2


  // find f0
  ierr = xf_Error(xf_AddFaceToHash(Node2Hash, 2, f0, xfe_False, 
				   -1, -1, -1, &hinfo0, &exists));
  if (ierr != xf_OK) return ierr;
  
  if (!exists){ 
    xf_printf("Error, nonconforming face = %d, %d not found in hash.\n", f0[0], f0[1]);
    return xf_Error(xf_FILE_READ_ERROR);
  }

  
  // find f1
  ierr = xf_Error(xf_AddFaceToHash(Node2Hash, 2, f1, xfe_False, 
				   -1, -1, -1, &hinfo1, &exists));
  if (ierr != xf_OK) return ierr;
  
  if (!exists){ 
    xf_printf("Error, nonconforming face = %d, %d not found in hash.\n", f1[0], f1[1]);
    return xf_Error(xf_FILE_READ_ERROR);
  }

  // find f2
  ierr = xf_Error(xf_AddFaceToHash(Node2Hash, 2, f2, xfe_False, 
				   -1, -1, -1, &hinfo2, &exists));
  if (ierr != xf_OK) return ierr;
  
  if (!exists){ 
    xf_printf("Error, nonconforming face = %d, %d not found in hash.\n", f2[0], f2[1]);
    return xf_Error(xf_FILE_READ_ERROR);
  }

  // coarse element info
  egrp0 = hinfo0->Group;
  elem0 = hinfo0->Elem;
  face0 = hinfo0->Face;

  // number of faces in the original group
  ierr = xf_Error(xf_Basis2nFace(Mesh->ElemGroup[egrp0].QBasis, &nface0));
  if (ierr != xf_OK) return ierr;

  // Q1 nodes on coarse face
  ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp0].QBasis, 
				   Mesh->ElemGroup[egrp0].QOrder, face0, 
				   &nfnode, nvec));
  if (ierr != xf_OK) return ierr;

  if (Mesh->ElemGroup[egrp0].Node[elem0][nvec[0]] == fvec[0]){
    pos1 = 1; fvec1[0] = fvec[0]; fvec1[1] = fvec[1];
    pos2 = 2; fvec2[0] = fvec[1]; fvec2[1] = fvec[2];
  }
  else if (Mesh->ElemGroup[egrp0].Node[elem0][nvec[0]] == fvec[2]){
    pos1 = 2; fvec1[0] = fvec[2]; fvec1[1] = fvec[1];
    pos2 = 1; fvec2[0] = fvec[1]; fvec2[1] = fvec[0];
  }
  else {
    xf_printf("%d %d %d\n", Mesh->ElemGroup[egrp0].Node[elem0][nvec[0]], fvec[0], fvec[2] );
    return xf_Error(xf_OUT_OF_BOUNDS); // too non-conforming for our taste
  }


  // increment # faces in elem0
  nface = Mesh->ElemGroup[egrp0].nFace[elem0];
  ierr = xf_Error(xf_Alloc( (void **) &nFaceNew, Mesh->ElemGroup[egrp0].nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<Mesh->ElemGroup[egrp0].nElem; i++) nFaceNew[i] = Mesh->ElemGroup[egrp0].nFace[i];
  nFaceNew[elem0] += 1;
  ierr = xf_Error(xf_VReAllocCopy2( (void ***) &Mesh->ElemGroup[egrp0].Face,
				    Mesh->ElemGroup[egrp0].nElem,
				    Mesh->ElemGroup[egrp0].nFace,
				    Mesh->ElemGroup[egrp0].nElem,
				    nFaceNew, sizeof(xf_Face)));
  if (ierr != xf_OK) return ierr;
  xf_Release( (void *) nFaceNew);
  
  Mesh->ElemGroup[egrp0].nFace[elem0] += 1;

  // Add interior faces f1 and f2 and set hangflag
  if (pos1 == 1)
    xf_AddInteriorFace(Mesh, hinfo1->Group, hinfo1->Elem, hinfo1->Face, egrp0, elem0, face0);
  else
    xf_AddInteriorFace(Mesh, hinfo2->Group, hinfo2->Elem, hinfo2->Face, egrp0, elem0, face0);
  Mesh->IFace[iface1 = Mesh->nIFace-1].HangNumber = nface0*1 + face0;

  if (pos2 == 2)
    xf_AddInteriorFace(Mesh, hinfo2->Group, hinfo2->Elem, hinfo2->Face, egrp0, elem0, nface);
  else
    xf_AddInteriorFace(Mesh, hinfo1->Group, hinfo1->Elem, hinfo1->Face, egrp0, elem0, nface);
  Mesh->IFace[iface2 = Mesh->nIFace-1].HangNumber = nface0*2 + face0;

  
  //Mesh->ElemGroup[egrp0].Face[elem0][face0] = iface1; // face in pos 1
  //Mesh->ElemGroup[egrp0].Face[elem0][nface] = iface2; // face in pos 2

  // set orientations of fine faces wrt coarse elem
  Mesh->IFace[iface1].OrientR = fvec1[1]<fvec1[0];
  Mesh->IFace[iface2].OrientR = fvec2[1]<fvec2[0];

  if (elem0 == 5){
    xf_printf("1: orient = %d, eL = %d, 2: orient = %d, eL = %d\n", 
	      Mesh->IFace[iface1].OrientR, Mesh->IFace[iface1].ElemL,
	      Mesh->IFace[iface2].OrientR, Mesh->IFace[iface2].ElemL);
  }


  // delete faces from hash
  ierr = xf_Error(xf_DeleteFaceFromHash(Node2Hash, 2, f0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DeleteFaceFromHash(Node2Hash, 2, f1));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DeleteFaceFromHash(Node2Hash, 2, f2));
  if (ierr != xf_OK) return ierr;


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AddPeriodicFaces
static int
xf_AddPeriodicFaces(xf_Mesh *Mesh, xf_FaceHash **Node2Hash, xf_PeriodicGroup *PG)
{
  /* Adds internal faces to Mesh according to remaining faces in
     Node2Hash and periodic groups (e.g. matching node pairs) in
     PG. */
  int ierr, i, k;
  int fvec[xf_MAXQ1FACENODE];
  int *MapB2A = NULL;
  enum xfe_Bool exists, skip;
  xf_FaceHash *hinfoB, *hinfoA, *hinfoBNext;

  // terminology: group A = first nodes in pair, group B = second nodes
  // note, periodic node pairs must be consistent!

  // create node mapping
  ierr = xf_Error(xf_Alloc((void **) &MapB2A, Mesh->nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<Mesh->nNode; i++) MapB2A[i] = -1;
  for (i=0; i<PG->nPeriodicNode; i++)
    MapB2A[PG->PeriodicNode[2*i+1]] = PG->PeriodicNode[2*i];
  

  // loop over group B nodes
  for (i=0; i<PG->nPeriodicNode; i++){
    hinfoB = Node2Hash[PG->PeriodicNode[2*i+1]];
    // loop over hash faces for current node in Group B
    while (hinfoB != NULL){
      hinfoBNext = hinfoB->Next;  // pointer to next

      for (k=0; k<hinfoB->nfnode; k++) fvec[k] = MapB2A[hinfoB->svec[k]];
      
      // check if all nodes in fvec are in group B
      for (k=0,skip=xfe_False; (k<hinfoB->nfnode)&&(!skip); k++) skip = skip || (fvec[k] < 0);
      
      if (!skip){
	/* Attempt to add face (fvec nodes) to hash list */
	ierr = xf_Error(xf_AddFaceToHash(Node2Hash, hinfoB->nfnode, fvec, xfe_False, 
					 hinfoB->Group, hinfoB->Elem, hinfoB->Face, 
					 &hinfoA, &exists));
	if (ierr != xf_OK) return ierr;

	if (!exists){
	  xf_printf("Periodic matching failed:\n");
	  for (k=0; k<hinfoB->nfnode; k++) xf_printf("%d ", hinfoB->svec[k]);
	  xf_printf("\n should match \n");
	  for (k=0; k<hinfoB->nfnode; k++) xf_printf("%d ", fvec[k]);
	  xf_printf("\n which does not exist in the remaining face hash. \n");
	  return xf_Error(xf_INPUT_ERROR);
	}

	// add new interior face as a result of matching hinfoA and hinfoB
	xf_AddInteriorFace(Mesh, hinfoA->Group, hinfoA->Elem, hinfoA->Face, 
			   hinfoB->Group, hinfoB->Elem, hinfoB->Face);


	// delete both faces from hash
	ierr = xf_Error(xf_DeleteFaceFromHash(Node2Hash, hinfoA->nfnode, hinfoA->svec));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_DeleteFaceFromHash(Node2Hash, hinfoB->nfnode, hinfoB->svec));
	if (ierr != xf_OK) return ierr;
      }

      hinfoB = hinfoBNext;

    } // end while hinfoB != NULL

  } // end for i looping over node pairs

  // release memory
  xf_Release((void *) MapB2A);

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ReadGriFile
int 
xf_ReadGriFile( const char *InputFile, char **InputStrings, xf_Mesh *Mesh)
{
  int ierr, nNode, nElemTot, dim, nLocFaceMax;
  int i, k, nBFG, iBFG, nBFace, nfnode;
  int iPeriodicGroup, iHangGroup, nHangGroup;
  int *nvec, *fvec;
  int nElemCurrent, nElem, QOrder, nface, nnode;
  int egrp, elem, face, iBFace, nIFaceMax;
  int iString = 0, nHang, iHang;
  char *line, line0[xf_MAXLONGLINELEN], QBasisString[xf_MAXSTRLEN];
  char title[xf_MAXSTRLEN];

  FILE *fgri = NULL;

  enum xfe_Bool exists, WarnedFlag = xfe_False;
  enum xfe_Bool PeriodicFlag, HangFlag, MoreToFile;
  enum xfe_BasisType QBasis;
  enum xfe_ShapeType Shape;
  xf_FaceHash **Node2Hash, *hinfo;
  xf_PeriodicGroup *PG;
  real Volume;
  
  // Check input
  if (((InputFile == NULL) && (InputStrings == NULL)) || 
      ((InputFile != NULL) && (InputStrings != NULL))) 
    return xf_Error(xf_INPUT_ERROR);

  if (InputFile != NULL){
    // Open InputFile
    fgri = fopen(InputFile, "r");
    if(!fgri){
      xf_printf("File %s not found.\n\n", InputFile);
      return xf_Error(xf_FILE_READ_ERROR);
    }
    xf_printf("Reading .gri file.\n");
  }
  else{
    // Start at the first string
    iString = 0;
  }
  

  /* Read in first line: nNode nElemTot [dim] */
  ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
  if (ierr != xf_OK) return ierr;


  ierr = sscanf(line, "%d %d %d", &nNode, &nElemTot, &dim);
  if (ierr != 3){
    ierr = sscanf(line, "%d %d", &nNode, &nElemTot);
    if (ierr != 2) return xf_Error(xf_FILE_READ_ERROR);
    dim = 3; // no specified type implies dim =3
  }
  Mesh->Dim = dim;


  /* Read in nodes */
  Mesh->nNode = nNode;
    
  ierr = xf_Error(xf_Alloc2((void ***) &Mesh->Coord, nNode, dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for(i=0; i<nNode; i++){
    ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ScanReal(line, dim, Mesh->Coord[i]));
    if (ierr != xf_OK) return ierr;
  }

  /* Create a node-based hash for looking up faces given nodes */
  ierr = xf_Error(xf_CreateNodeHash(nNode, &Node2Hash));
  if (ierr != xf_OK) return ierr;


  /* Read in # bounday face groups (nBFG) */
  ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
  if (ierr != xf_OK) return ierr;
  ierr = sscanf(line, "%d", &nBFG);
  if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  Mesh->nBFaceGroup = nBFG;
  ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup, nBFG, sizeof(xf_BFaceGroup)));
  if (ierr != xf_OK) return ierr;

 
  /* Read in BFaces, add to hash list */
  fvec = NULL;
  for (iBFG=0; iBFG<nBFG; iBFG++){
   
    /* Each BFG must begin with a line:  nBFace [nfnode] [BFGTitle] */
    ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
  
    ierr = sscanf(line, "%d %d %s", &nBFace, &nfnode, title);
    if (ierr != 3){
      sprintf(title, "BoundaryGroup%d\0", iBFG);
      ierr = sscanf(line, "%d %d", &nBFace, &nfnode);
      if (ierr != 2){
	ierr = sscanf(line, "%d", &nBFace);
	if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);
	nfnode = dim; // Assume tri/tet if nfnode not specified
      }
    }
    //xf_printf("nBFace = %d, nfnode = %d, title = %s\n", nBFace, nfnode, title);
    
    /* Create name for bfgroup */
    ierr = xf_Error(xf_AllocString(&Mesh->BFaceGroup[iBFG].Title, xf_MAXSTRLEN, title));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReAlloc((void **) &fvec, nfnode, sizeof(int)));
    if (ierr != xf_OK)  return ierr;
 
    Mesh->BFaceGroup[iBFG].nBFace = nBFace;
    
    /* Initialize .BFace */
    ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup[iBFG].BFace, nBFace, sizeof(xf_BFace)));
    if (ierr != xf_OK)  return ierr;


    /* Read boundary faces */
    for (iBFace=0; iBFace<nBFace; iBFace++){

      ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_ScanInt(line, nfnode, fvec));
      if (ierr != xf_OK) {
	xf_printf("nfnode = %d, line = %s\n", nfnode, line);
	return ierr;
      }

      for (i=0; i<nfnode; i++){
	fvec[i]--;  // convert to 0-starting index
	if ((fvec[i] < 0) || (fvec[i] >= nNode)){
	  xf_printf("Node index = %d out of range when reading boundary faces.\n", fvec[i]+1);
	  xf_printf("Note, the indexing should start at 1.\n");
	  return xf_Error(xf_FILE_READ_ERROR);
	}
      }

      ierr = xf_Error(xf_AddFaceToHash(Node2Hash, nfnode, fvec, xfe_True, iBFG, 
				       -1, iBFace, &hinfo, &exists));
      if (ierr != xf_OK) return ierr;
      
      if (exists){ // face should not exist in hash list
	xf_printf("Error, boundary face repeated:\n");
	for (i=0; i<nfnode; i++) xf_printf("%d ", fvec[i]+1);
	xf_printf("\n (1-base numbering)\n");
	return xf_Error(xf_FILE_READ_ERROR);
      }
      
    } // iBFace
  } // iBFG


  /* Read in elements. Each element group (ElemGroup) must begin with
     a line: nElem [Q [BasisName]].  Element groups of different Q and
     shapes are allowed. */

  Mesh->nElemGroup = 0;
  Mesh->ElemGroup = NULL;
  
  /* Allocate IFace array (more than necessary, will resize at end) */
  ierr = xf_Error(xf_GetLocFaceMax(&nLocFaceMax));
  if (ierr != xf_OK) return ierr;

  // sanity check
  if (nLocFaceMax > 20) 
    xf_printf("Warning, nLocFaceMax=%d; remove this warning if this is ok.\n", nLocFaceMax);
  
  nIFaceMax = nElemTot*nLocFaceMax; // not dividing by 2 to allow non-conforming case
  ierr = xf_Error(xf_Alloc( (void **) &Mesh->IFace, nIFaceMax, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;

  Mesh->nIFace = 0; // running total of number of interior faces
  nElemCurrent = 0; // running total of number of elements
  
  nvec = NULL;
  while (nElemCurrent < nElemTot){

    /* read number of elements and geometry order of group */
    ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
    
    ierr = sscanf(line, "%d %d %s", &nElem, &QOrder, QBasisString);
    if (ierr != 3){
      if (dim == 1)
	QBasis = xfe_SegLagrange;
      else if (dim == 2)
	QBasis = xfe_TriLagrange;
      else
	QBasis = xfe_TetLagrange;
      ierr = sscanf(line, "%d %d", &nElem, &QOrder);
      if (ierr != 2){
	QOrder = 1; // no specified order implies Q = 1
	ierr = sscanf(line, "%d", &nElem);
	if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);
      }
    }
    else{
      ierr = xf_Error(xf_Value2Enum(QBasisString, xfe_BasisName, xfe_BasisLast, (int *) &QBasis));
      if (ierr != xf_OK) return ierr;
    }
    
    if (nElem <= 0) return xf_Error(xf_FILE_READ_ERROR);  
    if (nElemCurrent+nElem > nElemTot) return xf_Error(xf_FILE_READ_ERROR);
  

    // number of faces for this element group
    ierr = xf_Error(xf_Basis2nFace(QBasis, &nface));
    if (ierr != xf_OK) return ierr;
    
    // number of nodes for this element group
    ierr = xf_Error(xf_Order2nNode(QBasis, QOrder, &nnode));
    if (ierr != xf_OK) return ierr;
    
    // ReAlloc nvec and fvec
    ierr = xf_Error(xf_ReAlloc((void **) &nvec, nnode, sizeof(int)));
    if (ierr != xf_OK)  return ierr;

    ierr = xf_Error(xf_ReAlloc((void **) &fvec, nnode, sizeof(int)));
    if (ierr != xf_OK)  return ierr;


    /* Allocate new element groups */
    egrp = Mesh->nElemGroup;
    Mesh->nElemGroup += 1;
    ierr = xf_Error(xf_ReAlloc((void **) &Mesh->ElemGroup, Mesh->nElemGroup, 
			       sizeof(xf_ElemGroup)));
    if (ierr!=xf_OK) return ierr;

    xf_InitElemGroup(Mesh->ElemGroup + egrp);

    Mesh->ElemGroup[egrp].QBasis = QBasis;
    Mesh->ElemGroup[egrp].QOrder = QOrder;
    Mesh->ElemGroup[egrp].nElem = nElem;
    ierr = xf_Error(xf_Alloc((void **) &Mesh->ElemGroup[egrp].nFace, nElem, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    for (elem=0; elem<nElem; elem++) Mesh->ElemGroup[egrp].nFace[elem] = nface;
    Mesh->ElemGroup[egrp].nNode = nnode;
    ierr = xf_Error(xf_VAlloc2((void ***) &Mesh->ElemGroup[egrp].Face, nElem, 
			       Mesh->ElemGroup[egrp].nFace, sizeof(xf_Face)));
    if (ierr!=xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc2((void ***) &Mesh->ElemGroup[egrp].Node, nElem, nnode, sizeof(int)));
    if (ierr!=xf_OK) return ierr;


    // Loop over elements and read them in, node by node
    for (elem=0; elem<nElem; elem++){

      // this can be a really long line
      ierr = xf_Error(xf_LongLineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_ScanInt(line, nnode, nvec));
      if (ierr != xf_OK) return ierr;

      for (i=0; i<nnode; i++){
	nvec[i] -= 1;  // to C numbering
	if ((nvec[i] < 0) || (nvec[i] >= nNode)){
	  xf_printf("Node index = %d out of range when reading elements.\n", nvec[i]+1);
	  xf_printf("Note, the indexing should start at 1.\n");
	  return xf_Error(xf_FILE_READ_ERROR);
	}
      }

      // make sure element volume (jacobian at quad points) is positive
      if ((QOrder == 1) && ((QBasis == xfe_TriLagrange) || (QBasis == xfe_TetLagrange))){
	ierr = xf_LinearElemJacobian(QBasis, QOrder, nvec, Mesh->Coord, &Volume);
	if ((ierr == xf_NOT_SUPPORTED) || (ierr == xf_INPUT_ERROR)){
	  if ((!WarnedFlag) && (InputFile != NULL))
	    xf_printf("** In ReadGri: NOT CHECKING %s VOLUMES **\n", xfe_BasisName[QBasis]);
	  WarnedFlag = xfe_True;
	}
	else if (ierr != xf_OK) return xf_Error(ierr);
	if ((ierr == xf_OK) && (Volume <= 0.)){
	  if (QBasis == xfe_TriLagrange){
	    xf_printf("Swapping nodes to fix negative volume in element = %d\n", elem);
	    swap(nvec[0], nvec[1], k);
	  }
	  else{
	    xf_printf("Negative element Jacobian (volume) error.\n");
	    xf_printf("ElemGroup=%d, Elem=%d, detJ = %.10E\n", egrp, elem, Volume);
	    return xf_Error(xf_NEGATIVE_JACOBIAN);
	  }
	}
      }
      else{
	if ((!WarnedFlag) && (InputFile != NULL))
	  xf_printf("** In ReadGri: NOT CHECKING VALIDITY OF Q>1 VOLUMES **\n");
	WarnedFlag = xfe_True;
      }

      // Add nodes
      for (k = 0; k<nnode; k++) Mesh->ElemGroup[egrp].Node[elem][k] = nvec[k];

      
      /* loop over the faces and check the hash table */
      for (face=0; face< Mesh->ElemGroup[egrp].nFace[elem]; face++){

	// local nodes on face
	ierr = xf_Error(xf_Q1NodesOnFace(QBasis, QOrder, face, &nfnode, fvec));
	if (ierr != xf_OK) return ierr;
	
	// convert to global nodes
	for (k=0; k<nfnode; k++) fvec[k] = nvec[fvec[k]]; 

	/* Attempt to add face (fvec nodes) to hash list */
	ierr = xf_Error(xf_AddFaceToHash(Node2Hash, nfnode, fvec, xfe_False, 
					 egrp, elem, face, &hinfo, &exists));
	if (ierr != xf_OK) return ierr;

	if (exists){ // face already exists in hash list
	  if (hinfo->nVisit != 1){
	    xf_printf("Error, more than two elements share a face, or a\n");
	    xf_printf("boundary face is referenced by more than one element.\n");
	    return xf_Error(xf_FILE_READ_ERROR);
	  }
	  
	  // link elem to bface or to iface

	  if (hinfo->BFlag == xfe_True){ // boundary face
	    xf_InitBFace(Mesh->BFaceGroup[hinfo->Group].BFace + hinfo->Face);
	    Mesh->BFaceGroup[hinfo->Group].BFace[hinfo->Face].ElemGroup = egrp;
	    Mesh->BFaceGroup[hinfo->Group].BFace[hinfo->Face].Elem      = elem;
	    Mesh->BFaceGroup[hinfo->Group].BFace[hinfo->Face].Face      = face;
	    Mesh->ElemGroup[egrp].Face[elem][face].Group  = hinfo->Group;
	    Mesh->ElemGroup[egrp].Face[elem][face].Number = hinfo->Face;
	  }
	  else{ // interior face
	    xf_AddInteriorFace(Mesh, hinfo->Group, hinfo->Elem, hinfo->Face, egrp, elem, face);
	  }

	  ierr = xf_Error(xf_DeleteFaceFromHash(Node2Hash, nfnode, fvec));
	  if (ierr != xf_OK) return ierr;
	}

      } // face
     
    } // elem

    nElemCurrent += nElem;
 
  } // end while (nElemCurrent < nElemTot)

  
  // Check number of elements
  if (nElemCurrent != nElemTot){
    xf_printf("Error in number of elements.\n");
    xf_printf("nElemCurrent = %d, nElemTot = %d\n", nElemCurrent, nElemTot);
    return xf_Error(xf_FILE_READ_ERROR);
  }  

  // Do we have any more lines?
  ierr = xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line);
  MoreToFile = (ierr == xf_OK);

  PeriodicFlag = xfe_False;
  HangFlag     = xfe_False;
  if (MoreToFile){
    // check for periodic or hanging node flag
    ierr = sscanf(line, "%d %s", &k, title);
    if (ierr == 2){
      if (strncmp(title, "PeriodicGroup", 13) == 0){
	PeriodicFlag = xfe_True; Mesh->nPeriodicGroup = k;
      }
      else if (strncmp(title, "HangGroup", 9) == 0){
	HangFlag = xfe_True; nHangGroup = k;
      }
    }
  }

  // Read in periodic info if given
  if (PeriodicFlag){
    // allocate periodic groups
    ierr = xf_Error(xf_Alloc((void **) &Mesh->PeriodicGroup, Mesh->nPeriodicGroup, 
			     sizeof(xf_PeriodicGroup)));
    if (ierr != xf_OK) return ierr;

    for (iPeriodicGroup=0; iPeriodicGroup<Mesh->nPeriodicGroup; iPeriodicGroup++){
      PG = Mesh->PeriodicGroup+iPeriodicGroup; // pointer to current periodic group
      ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
      if (ierr != xf_OK) return ierr;
      ierr = sscanf(line, "%d %s", &PG->nPeriodicNode, title);
      if (ierr != 2) return xf_Error(xf_FILE_READ_ERROR);
      ierr = xf_Error(xf_Value2Enum(title, xfe_PeriodicityName, xfe_PeriodicityLast, 
				    (int *) &PG->Periodicity));
      if (ierr != xf_OK) return ierr;

      // only Translational currently supported
      if (PG->Periodicity != xfe_PeriodicityTranslational)
	return xf_Error(xf_NOT_SUPPORTED);

      // allocate memory for node pairs
      ierr = xf_Error(xf_Alloc((void **) &PG->PeriodicNode, 2*PG->nPeriodicNode, sizeof(int)));
      if (ierr != xf_OK) return ierr;

      // read in node pairs
      for(i=0; i<PG->nPeriodicNode; i++){
	ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
	if (ierr != xf_OK) return ierr;
	ierr = sscanf(line, "%d %d", PG->PeriodicNode+2*i, PG->PeriodicNode+2*i+1);
	if (ierr != 2) return xf_Error(xf_FILE_READ_ERROR);
	for (k=0; k<2; k++) PG->PeriodicNode[2*i+k]--; // convert to C numbering
      }
      
      // Generate new ifaces using remaining hash faces and periodic node pairs
      ierr = xf_Error(xf_AddPeriodicFaces(Mesh, Node2Hash, PG));
      if (ierr != xf_OK) return ierr;

    } // end for iPeriodicGroup

  } // end if PeriodicFlag

  // Read in hanging-node info if given
  if (HangFlag){
    // Read in non-conforming element information if the hash is not empty
    //while (!xf_HashIsEmpty(Node2Hash,Mesh->nNode)){
    for (iHangGroup=0; iHangGroup<nHangGroup; iHangGroup++){
      ierr = xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line);
      if (ierr == xf_OK){ // read line or break and error out when hash is checked
	ierr = sscanf(line, "%d %s", &nHang, title);
	if (ierr != 2) break;
      }
      else break;

      ierr = xf_Error(xf_Value2Enum(title, xfe_ShapeName, xfe_ShapeLast, (int *) &Shape));
      if (ierr != xf_OK) return ierr;
    
      switch(Shape){ // determine how many non-conforming nodes to read
      case xfe_Segment: nnode = 3; break;
      default:
	return xf_Error(xf_NOT_SUPPORTED);
	break;
      }
    
      ierr = xf_Error(xf_ReAlloc((void **) &fvec, nnode, sizeof(int)));
      if (ierr != xf_OK)  return ierr;

      for (iHang=0; iHang<nHang; iHang++){
	ierr = xf_Error(xf_LineFromFileOrStrings(fgri, InputStrings, &iString, line0, &line));
	if (ierr != xf_OK) return ierr;
      
	ierr = xf_Error(xf_ScanInt(line, nnode, fvec));
	if (ierr != xf_OK) return ierr;

	for (i=0; i<nnode; i++) fvec[i]--; // convert to C numbering

	ierr = xf_Error(xf_AddHangingFace(Mesh, Node2Hash, Shape, nnode, fvec));
	if (ierr != xf_OK) return ierr;
      
      } // iHang  
    } // end for-loop reading non-conforming faces 
  } // end HangFlag    
    
  // Make sure no faces are left in the hash; print out remaining faces if any
  for (i=0, face=0; i<Mesh->nNode; i++){
    if (Node2Hash[i] != NULL){
      hinfo = Node2Hash[i];
      while (hinfo != NULL){
	for (k=0; k<hinfo->nfnode; k++)
	  xf_printf("%d ", hinfo->svec[k]+1);
	xf_printf("\n");
	face++;
	hinfo = hinfo->Next;
      }
    }
  } // i
  if (face != 0){
    xf_printf("Mesh connectivity error: the above %d face(s) remain(s) in the hash.\n\n",
	      face);
    return xf_Error(xf_FILE_READ_ERROR);
  }

  // Resize IFace
  if (Mesh->nIFace > nIFaceMax) return xf_Error(xf_FILE_READ_ERROR);
  ierr = xf_Error(xf_ReAlloc( (void **) &Mesh->IFace, Mesh->nIFace, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;

  // update face orientation info  
  ierr = xf_Error(xf_UpdateFaceOrient(Mesh));
  if (ierr != xf_OK) return ierr;

    
  xf_Release((void *) fvec);
  xf_Release((void *) nvec);
  xf_Release((void *) Node2Hash);

  if (fgri != NULL) fclose(fgri);

  return xf_OK;
}



/******************************************************************/
int
PeriodMatch(real *xp1, real *xp2, int orient)
{
   real eps = 1.e-8;

   if(orient==1) //horizontal
   {
      if(fabs(xp1[1] - xp2[1]) < eps)
         return 1;
      else
         return 0;
   }
   
   if(orient==2) //vertical 
   {
      if(fabs(xp1[0] - xp2[0]) < eps)
         return 1;
      else
         return 0;
   }
}
/******************************************************************/
//   FUNCTION Definition: xf_WriteGriFile
int 
xf_WriteGriFile( xf_Mesh *Mesh, char *OutputFile )
{
  int ierr, i, j, dim, d;
  int nelemtot;
  int ibfgrp, ibface, egrp, elem, face;
  int nfnode, fvec[xf_MAXQ1FACENODE];
  xf_PeriodicGroup *PG;
  FILE *fgri;

  //~~~by Yu
  //hook-in for periodical boundary
  char input[xf_MAXSTRLEN];
  enum xfe_Bool WhePeriodic, inList;
  int offset, ipair, pair[4], tmp, numIn;
  int nPeriodicNode, *PeriodNodeList[2];

  if ((fgri =fopen(OutputFile ,"w"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);

  // total number of elements
  ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;

  // header + nodes
  fprintf(fgri, "%d %d %d\n", Mesh->nNode, nelemtot, (dim = Mesh->Dim));
  for (i=0; i<Mesh->nNode; i++){
    for (d=0; d<dim; d++)
      fprintf(fgri, "%.15E ", Mesh->Coord[i][d]);
    fprintf(fgri, "\n");
  }
 
  //~~~by Yu 
  //only support fully peridoical box case
  xf_printf("Are we doing periodical square:\n");
  scanf("%s", input); 
  if(strcmp(input, "yes") == 0)
  {//if we do periodical square
     fprintf(fgri, "%d\n", 0);
     if(Mesh->nBFaceGroup != 4) 
        return xf_NOT_SUPPORTED;

     for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
     
        if(strcmp(Mesh->BFaceGroup[ibfgrp].Title, "left") == 0)
           pair[0] = ibfgrp;
        if(strcmp(Mesh->BFaceGroup[ibfgrp].Title, "right") == 0)
           pair[1] = ibfgrp;
        if(strcmp(Mesh->BFaceGroup[ibfgrp].Title, "top") == 0)
           pair[2] = ibfgrp;
        if(strcmp(Mesh->BFaceGroup[ibfgrp].Title, "bottom") == 0)
           pair[3] = ibfgrp;
     } 
        
     //match node list
     nPeriodicNode = Mesh->BFaceGroup[pair[0]].nBFace +1 ; 
     ierr = xf_Error(xf_Alloc((void **) &PeriodNodeList[0], 2*nPeriodicNode, sizeof(int))); 
     if (ierr != xf_OK) return ierr;
     nPeriodicNode = Mesh->BFaceGroup[pair[2]].nBFace +1 ; 
     ierr = xf_Error(xf_Alloc((void **) &PeriodNodeList[1], 2*nPeriodicNode, sizeof(int))); 
     if (ierr != xf_OK) return ierr;

     numIn = 0;
     for(ipair=0; ipair < 4; ipair ++){ if(ipair==2) numIn=0;
     for(ibface=0; ibface<Mesh->BFaceGroup[pair[ipair]].nBFace; ibface++)
     {
        egrp = Mesh->BFaceGroup[pair[ipair]].BFace[ibface].ElemGroup;
        elem = Mesh->BFaceGroup[pair[ipair]].BFace[ibface].Elem;
        face = Mesh->BFaceGroup[pair[ipair]].BFace[ibface].Face;
      
        ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, Mesh->ElemGroup[egrp].QOrder, 
				         face, &nfnode, fvec));
        if (ierr != xf_OK) return ierr;

      
        for (i=0; i<nfnode; i++)
        {   
           tmp = Mesh->ElemGroup[egrp].Node[elem][fvec[i]]+1;
           inList = xfe_False;
           for (j=0; j<numIn; j++)
              if(PeriodNodeList[ipair/2][j] == tmp)
              {
                 inList = xfe_True;
                 break;
              }

           if(!inList){
              PeriodNodeList[ipair/2][numIn] = tmp;
              numIn++;
           }
        }
     }
     }

     //match node coordinates
     nPeriodicNode = Mesh->BFaceGroup[pair[0]].nBFace + 1;
     for(j=0; j<nPeriodicNode; j++)
     {
        for(i=nPeriodicNode; i<2*nPeriodicNode; i++)
        if(PeriodMatch(Mesh->Coord[PeriodNodeList[0][j]-1], Mesh->Coord[PeriodNodeList[0][i]-1], 1) == 1)
        {
           tmp = PeriodNodeList[0][i];
           PeriodNodeList[0][i] = PeriodNodeList[0][j+nPeriodicNode];
           PeriodNodeList[0][j+nPeriodicNode] = tmp;
        }
     }
     nPeriodicNode = Mesh->BFaceGroup[pair[2]].nBFace + 1;
     offset = 0;
     for(j=offset; j<offset+nPeriodicNode; j++)
     {
        for(i=offset+nPeriodicNode; i<offset+2*nPeriodicNode; i++)
        if(PeriodMatch(Mesh->Coord[PeriodNodeList[1][j]-1], Mesh->Coord[PeriodNodeList[1][i]-1], 2) == 1)
        {
           tmp = PeriodNodeList[1][i];
           PeriodNodeList[1][i] = PeriodNodeList[1][j+nPeriodicNode];
           PeriodNodeList[1][j+nPeriodicNode] = tmp;
        }
     }
     
  }
  else
  {

  // bfacegroups
  fprintf(fgri, "%d\n", Mesh->nBFaceGroup);
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    if (Mesh->BFaceGroup[ibfgrp].nBFace > 0){
      egrp = Mesh->BFaceGroup[ibfgrp].BFace[0].ElemGroup;
      elem = Mesh->BFaceGroup[ibfgrp].BFace[0].Elem;
      face = Mesh->BFaceGroup[ibfgrp].BFace[0].Face;
      // local nodes on face
      ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, Mesh->ElemGroup[egrp].QOrder, 
				       face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;
    }
    else nfnode = 0;
    fprintf(fgri, "%d %d %s\n", Mesh->BFaceGroup[ibfgrp].nBFace, nfnode, Mesh->BFaceGroup[ibfgrp].Title);
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      egrp = Mesh->BFaceGroup[ibfgrp].BFace[ibface].ElemGroup;
      elem = Mesh->BFaceGroup[ibfgrp].BFace[ibface].Elem;
      face = Mesh->BFaceGroup[ibfgrp].BFace[ibface].Face;
      // local nodes on face
      ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, Mesh->ElemGroup[egrp].QOrder, 
				       face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<nfnode; i++)
	fprintf(fgri, "%d ", Mesh->ElemGroup[egrp].Node[elem][fvec[i]]+1);
      fprintf(fgri, "\n");
    }
  }

  }//else

  // element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    fprintf(fgri, "%d %d %s\n", Mesh->ElemGroup[egrp].nElem, Mesh->ElemGroup[egrp].QOrder,
	    xfe_BasisName[Mesh->ElemGroup[egrp].QBasis]);
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      for (i=0; i<Mesh->ElemGroup[egrp].nNode; i++)
	fprintf(fgri, "%d ", Mesh->ElemGroup[egrp].Node[elem][i]+1);
      fprintf(fgri, "\n");
    }
  }
  
  if(strcmp(input, "yes") == 0)
  {
     fprintf(fgri, "2 PeriodicGroup\n");
     //left - right
     nPeriodicNode = Mesh->BFaceGroup[pair[0]].nBFace + 1;
     fprintf(fgri, "%d Translational\n", nPeriodicNode);
     for(i=0; i<nPeriodicNode; i++)
        fprintf(fgri, "%d %d\n", PeriodNodeList[0][i], PeriodNodeList[0][i+nPeriodicNode]);
     
     nPeriodicNode = Mesh->BFaceGroup[pair[2]].nBFace + 1;
     fprintf(fgri, "%d Translational\n", nPeriodicNode);
     offset = 0;
     for(i=offset; i<offset+nPeriodicNode; i++)
        fprintf(fgri, "%d %d\n", PeriodNodeList[1][i], PeriodNodeList[1][i+nPeriodicNode]);
  }

  // periodic groups
  if (Mesh->nPeriodicGroup > 0){
    fprintf(fgri, "%d PeriodicGroup\n", Mesh->nPeriodicGroup);
    for (i=0; i<Mesh->nPeriodicGroup; i++){
      PG = Mesh->PeriodicGroup+i; 
      fprintf(fgri, "%d %s\n", PG->nPeriodicNode, xfe_PeriodicityName[PG->Periodicity]);
      for (j=0; j<PG->nPeriodicNode; j++)
        fprintf(fgri, "%d %d\n", PG->PeriodicNode[2*j+0]+1, PG->PeriodicNode[2*j+1]+1);
    } // i
  }
  
  fclose(fgri);

  xf_Release((void *) PeriodNodeList[0]);
  xf_Release((void *) PeriodNodeList[1]);
  return xf_OK;
}


// Gmsh support functions
#include "xf_MeshGmsh.c"
// Fluent msh support functions
#include "xf_MeshFmsh.c"
// Bamg support functions
#include "xf_MeshBamg.c"


/******************************************************************/
//   FUNCTION Definition: xf_WriteFaceDataBinary
static int 
xf_WriteFaceDataBinary( xf_Face Face, FILE *fid){
  int ierr;

  if (fwrite(&Face.Group, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  if (fwrite(&Face.Number, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadFaceDataBinary
static int 
xf_ReadFaceDataBinary( FILE *fid, xf_Face *Face){
  int ierr;

  if (fread(&Face->Group, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  if (fread(&Face->Number, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteCutFaceDataBinary
static int 
xf_WriteCutFaceDataBinary( xf_CutFaceData *CutFaceData, FILE *fid){
  int ierr;

  ierr = xf_Error(xf_WriteQuadDataBinary(CutFaceData->QuadData, fid));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_WriteFaceDataBinary(CutFaceData->OrigFace, fid));
  if (ierr != xf_OK) return ierr;
  
  if (fwrite(&CutFaceData->GeomIndex, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadCutFaceDataBinary
static int 
xf_ReadCutFaceDataBinary( FILE *fid, xf_CutFaceData *CutFaceData){
  int ierr;

  ierr = xf_Error(xf_CreateQuadData(&CutFaceData->QuadData));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_ReadQuadDataBinary(fid, CutFaceData->QuadData));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_ReadFaceDataBinary(fid, &CutFaceData->OrigFace));
  if (ierr != xf_OK) return ierr;
  
  if (fread(&CutFaceData->GeomIndex, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteCutElemDataBinary
static int 
xf_WriteCutElemDataBinary( xf_CutElemData *CutElemData, FILE *fid){
  int ierr, nnode, i, dim;

  ierr = xf_Error(xf_WriteQuadDataBinary(CutElemData->QuadData, fid));
  if (ierr != xf_OK) return ierr;

  // basis of interpolation element
  ierr = xf_Error(xf_WriteStringBinary(xfe_BasisName[CutElemData->QBasis], fid));
  if (ierr != xf_OK) return ierr;
  
  // number of nodes for interpolation element (QOrder = 1)
  ierr = xf_Error(xf_Order2nNode(CutElemData->QBasis, 1, &nnode));
  if (ierr != xf_OK) return ierr;

  // dimension of coords
  ierr = xf_Error(xf_Basis2Dim(CutElemData->QBasis, &dim));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<nnode; i++)
    if (fwrite(CutElemData->QCoord[i], sizeof(real), dim, fid) != dim) 
      return xf_Error(xf_FILE_WRITE_ERROR);
    
  if (fwrite(&CutElemData->GeomIndex, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadCutElemDataBinary
static int 
xf_ReadCutElemDataBinary( FILE *fid, xf_CutElemData *CutElemData){
  int ierr, nnode, i, dim;

  ierr = xf_Error(xf_CreateQuadData(&CutElemData->QuadData));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_ReadQuadDataBinary(fid, CutElemData->QuadData));
  if (ierr != xf_OK) return ierr;

  // basis of interpolation element
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BasisName, xfe_BasisLast, 
				    (int *) &CutElemData->QBasis));
  if (ierr != xf_OK) return ierr;
  
  // number of nodes for interpolation element (QOrder = 1)
  ierr = xf_Error(xf_Order2nNode(CutElemData->QBasis, 1, &nnode));
  if (ierr != xf_OK) return ierr;

  // dimension of coords
  ierr = xf_Error(xf_Basis2Dim(CutElemData->QBasis, &dim));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc2((void ***) &CutElemData->QCoord, nnode, dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<nnode; i++)
    if (fread(CutElemData->QCoord[i], sizeof(real), dim, fid) != dim) 
      return xf_Error(xf_FILE_READ_ERROR);
    
  if (fread(&CutElemData->GeomIndex, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteMeshBinarySerial
static int 
xf_WriteMeshBinarySerial( xf_Mesh *Mesh, FILE *fid){

  int ierr, rev, len, i, j, k;
  int si, sr, iflag;
  enum xfe_Bool flag;
  xf_PeriodicGroup *PG;

  si = sizeof(int);
  sr = sizeof(real);

  rev = 2;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // write dim and nodes
  if (fwrite(&Mesh->Dim, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  if (fwrite(&Mesh->nNode, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  for (i=0; i<Mesh->nNode; i++)
    if (fwrite(Mesh->Coord[i], sr, Mesh->Dim, fid) != Mesh->Dim) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  

  // write IFaces
  if (fwrite(&Mesh->nIFace, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  for (i=0; i<Mesh->nIFace; i++){
    ierr = 0;
    ierr += fwrite(&Mesh->IFace[i].ElemGroupL, si, 1, fid);
    ierr += fwrite(&Mesh->IFace[i].ElemL     , si, 1, fid);
    ierr += fwrite(&Mesh->IFace[i].FaceL     , si, 1, fid);
    ierr += fwrite(&Mesh->IFace[i].ElemGroupR, si, 1, fid);
    ierr += fwrite(&Mesh->IFace[i].ElemR     , si, 1, fid);
    ierr += fwrite(&Mesh->IFace[i].FaceR     , si, 1, fid);
    ierr += fwrite(&Mesh->IFace[i].OrientL   , si, 1, fid);
    ierr += fwrite(&Mesh->IFace[i].OrientR   , si, 1, fid);
    ierr += fwrite(&Mesh->IFace[i].HangNumber, si, 1, fid);
    if (ierr != 9) return xf_Error(xf_FILE_WRITE_ERROR);
    
    iflag = (Mesh->IFace[i].CutFaceData != NULL);
    if (fwrite(&iflag, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
    if (iflag){
      ierr = xf_Error(xf_WriteCutFaceDataBinary(Mesh->IFace[i].CutFaceData, fid));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // write BFaces
  if (fwrite(&Mesh->nBFaceGroup, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  for (i=0; i<Mesh->nBFaceGroup; i++){
    ierr = xf_Error(xf_WriteStringBinary(Mesh->BFaceGroup[i].Title, fid));
    if (ierr != xf_OK) return ierr;
    if (fwrite(&Mesh->BFaceGroup[i].nBFace, si, 1, fid) != 1) 
      return xf_Error(xf_FILE_WRITE_ERROR);
    for (j=0; j<Mesh->BFaceGroup[i].nBFace; j++){
      ierr = 0;
      ierr += fwrite(&Mesh->BFaceGroup[i].BFace[j].ElemGroup, si, 1, fid);
      ierr += fwrite(&Mesh->BFaceGroup[i].BFace[j].Elem     , si, 1, fid);
      ierr += fwrite(&Mesh->BFaceGroup[i].BFace[j].Face     , si, 1, fid);
      ierr += fwrite(&Mesh->BFaceGroup[i].BFace[j].Orient   , si, 1, fid);
      if (ierr != 4) return xf_Error(xf_FILE_WRITE_ERROR);

      iflag = (Mesh->BFaceGroup[i].BFace[j].CutFaceData != NULL);
      if (fwrite(&iflag, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
      if (iflag){
	ierr = xf_Error(xf_WriteCutFaceDataBinary(Mesh->BFaceGroup[i].BFace[j].CutFaceData, fid));
	if (ierr != xf_OK) return ierr;
      }
    }
  }
  

  // write ElemGroups
  if (fwrite(&Mesh->nElemGroup, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  for (i=0; i<Mesh->nElemGroup; i++){
    ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[Mesh->ElemGroup[i].CutFlag], fid));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_WriteStringBinary(xfe_BasisName[Mesh->ElemGroup[i].QBasis], fid));
    if (ierr != xf_OK) return ierr;
    ierr = 0;
    ierr += fwrite(&Mesh->ElemGroup[i].QOrder, si, 1, fid);
    ierr += fwrite(&Mesh->ElemGroup[i].nElem , si, 1, fid);
    if (ierr != 2) return xf_Error(xf_FILE_WRITE_ERROR);

    for (j=0,ierr=0; j<Mesh->ElemGroup[i].nElem; j++)
      ierr += fwrite(Mesh->ElemGroup[i].nFace+j, si, 1, fid);
    if (ierr != Mesh->ElemGroup[i].nElem) return xf_Error(xf_FILE_WRITE_ERROR);

    ierr = fwrite(&Mesh->ElemGroup[i].nNode , si, 1, fid);
    if (ierr != 1) return xf_Error(xf_FILE_WRITE_ERROR);

    // is cut-cell data present?
    iflag = (Mesh->ElemGroup[i].CutElemData != NULL);
    if (fwrite(&iflag, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
    
    for (j=0; j<Mesh->ElemGroup[i].nElem; j++){
      
      for (k=0; k<Mesh->ElemGroup[i].nFace[j]; k++){
	ierr = xf_Error(xf_WriteFaceDataBinary(Mesh->ElemGroup[i].Face[j][k], fid));
	if (ierr != xf_OK) return ierr;
      }

      if (fwrite(Mesh->ElemGroup[i].Node[j], si, Mesh->ElemGroup[i].nNode, fid) != Mesh->ElemGroup[i].nNode) 
	return xf_Error(xf_FILE_WRITE_ERROR);

      if (iflag){
	ierr = xf_Error(xf_WriteCutElemDataBinary(Mesh->ElemGroup[i].CutElemData+j, fid));
	if (ierr != xf_OK) return ierr;
      }
    }
  }


  // write periodic groups
  if (fwrite(&Mesh->nPeriodicGroup, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  for (i=0; i<Mesh->nPeriodicGroup; i++){
    PG = Mesh->PeriodicGroup+i;
    ierr = xf_Error(xf_WriteStringBinary(xfe_PeriodicityName[PG->Periodicity], fid));
    if (ierr != xf_OK) return ierr;
    ierr = fwrite(&PG->nPeriodicNode, si, 1, fid);
    if (ierr != 1) return xf_Error(xf_FILE_WRITE_ERROR);
    ierr = fwrite(PG->PeriodicNode , si, 2*PG->nPeriodicNode, fid);
    if (ierr != 2*PG->nPeriodicNode) return xf_Error(xf_FILE_WRITE_ERROR);
  }

  // write background mesh
  flag = (Mesh->BackgroundMesh != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteMeshBinarySerial(Mesh->BackgroundMesh, fid));
    if (ierr != xf_OK) return ierr;
  }

  // write Motion
  flag = (Mesh->Motion != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteMeshMotionBinarySerial(Mesh->Motion, fid));
    if (ierr != xf_OK) return ierr;
  }


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadMeshBinarySerial
static int 
xf_ReadMeshBinarySerial( FILE *fid, xf_Mesh *Mesh){

  int ierr, rev, len, i, j, k;
  int si, sr, iflag, nface;
  enum xfe_Bool flag;
  xf_PeriodicGroup *PG;

  si = sizeof(int);
  sr = sizeof(real);

  // read + check revision number
  rev = 0;
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev > 2) return xf_Error(xf_FILE_READ_ERROR);

  // read dim and nodes
  if (fread(&Mesh->Dim, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  if (fread(&Mesh->nNode, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  ierr = xf_Error(xf_Alloc2((void ***) &Mesh->Coord, Mesh->nNode, Mesh->Dim, sr));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<Mesh->nNode; i++)
    if (fread(Mesh->Coord[i], sr, Mesh->Dim, fid) != Mesh->Dim) 
      return xf_Error(xf_FILE_READ_ERROR);
  

  // read IFaces
  if (fread(&Mesh->nIFace, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  ierr = xf_Error(xf_Alloc( (void **) &Mesh->IFace, Mesh->nIFace, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<Mesh->nIFace; i++){
    xf_InitIFace(Mesh->IFace+i);
    ierr = 0;
    ierr += fread(&Mesh->IFace[i].ElemGroupL, si, 1, fid);
    ierr += fread(&Mesh->IFace[i].ElemL     , si, 1, fid);
    ierr += fread(&Mesh->IFace[i].FaceL     , si, 1, fid);
    ierr += fread(&Mesh->IFace[i].ElemGroupR, si, 1, fid);
    ierr += fread(&Mesh->IFace[i].ElemR     , si, 1, fid);
    ierr += fread(&Mesh->IFace[i].FaceR     , si, 1, fid);
    ierr += fread(&Mesh->IFace[i].OrientL   , si, 1, fid);
    ierr += fread(&Mesh->IFace[i].OrientR   , si, 1, fid);
    ierr += fread(&Mesh->IFace[i].HangNumber, si, 1, fid);
    if (ierr != 9) return xf_Error(xf_FILE_READ_ERROR);
    
    if (fread(&iflag, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);

    if (iflag){
      ierr = xf_Error(xf_CreateCutFaceData(&Mesh->IFace[i].CutFaceData));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReadCutFaceDataBinary(fid, Mesh->IFace[i].CutFaceData));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // read BFaces
  if (fread(&Mesh->nBFaceGroup, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup, Mesh->nBFaceGroup, 
			   sizeof(xf_BFaceGroup)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<Mesh->nBFaceGroup; i++){
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &Mesh->BFaceGroup[i].Title));
    if (ierr != xf_OK) return ierr;
    if (fread(&Mesh->BFaceGroup[i].nBFace, si, 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    
    ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup[i].BFace, 
			     Mesh->BFaceGroup[i].nBFace, sizeof(xf_BFace)));
    if (ierr != xf_OK)  return ierr;

    for (j=0; j<Mesh->BFaceGroup[i].nBFace; j++){
      xf_InitBFace(Mesh->BFaceGroup[i].BFace + j);
      ierr = 0;
      ierr += fread(&Mesh->BFaceGroup[i].BFace[j].ElemGroup, si, 1, fid);
      ierr += fread(&Mesh->BFaceGroup[i].BFace[j].Elem     , si, 1, fid);
      ierr += fread(&Mesh->BFaceGroup[i].BFace[j].Face     , si, 1, fid);
      ierr += fread(&Mesh->BFaceGroup[i].BFace[j].Orient   , si, 1, fid);
      if (ierr != 4) return xf_Error(xf_FILE_READ_ERROR);

      if (fread(&iflag, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (iflag){
	ierr = xf_Error(xf_CreateCutFaceData(&Mesh->BFaceGroup[i].BFace[j].CutFaceData));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_WriteCutFaceDataBinary(Mesh->BFaceGroup[i].BFace[j].CutFaceData, fid));
	if (ierr != xf_OK) return ierr;
      }
    }
  }
  

  // read ElemGroups
  if (fread(&Mesh->nElemGroup, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  ierr = xf_Error(xf_Alloc((void **) &Mesh->ElemGroup, Mesh->nElemGroup, 
			   sizeof(xf_ElemGroup)));
  if (ierr!=xf_OK) return ierr;

  for (i=0; i<Mesh->nElemGroup; i++){
    xf_InitElemGroup(Mesh->ElemGroup + i);

    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, 
				      (int *) &Mesh->ElemGroup[i].CutFlag));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BasisName, xfe_BasisLast, 
				      (int *) &Mesh->ElemGroup[i].QBasis));
    if (ierr != xf_OK) return ierr;

    ierr = 0;
    ierr += fread(&Mesh->ElemGroup[i].QOrder, si, 1, fid);
    ierr += fread(&Mesh->ElemGroup[i].nElem , si, 1, fid);
    if (ierr != 2) return xf_Error(xf_FILE_READ_ERROR);

    ierr = xf_Error(xf_Alloc((void **) &Mesh->ElemGroup[i].nFace, 
			     Mesh->ElemGroup[i].nElem, sizeof(int)));
    if (ierr!=xf_OK) return ierr;

    if (rev == 0){
      ierr = fread(&nface , si, 1, fid);
      if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);
      for (j=0; j<Mesh->ElemGroup[i].nElem; j++) Mesh->ElemGroup[i].nFace[j] = nface;
    }
    else{
      for (j=0, ierr=0; j<Mesh->ElemGroup[i].nElem; j++)
	ierr += fread(Mesh->ElemGroup[i].nFace+j, si, 1, fid);
      if (ierr != Mesh->ElemGroup[i].nElem) return xf_Error(xf_FILE_READ_ERROR);
    }


    ierr = fread(&Mesh->ElemGroup[i].nNode , si, 1, fid);
    if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);

    // allocate Face and Node
    ierr = xf_Error(xf_VAlloc2((void ***) &Mesh->ElemGroup[i].Face, Mesh->ElemGroup[i].nElem, 
			       Mesh->ElemGroup[i].nFace, sizeof(xf_Face)));
    if (ierr!=xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc2((void ***) &Mesh->ElemGroup[i].Node, Mesh->ElemGroup[i].nElem, 
			      Mesh->ElemGroup[i].nNode, sizeof(int)));
    if (ierr!=xf_OK) return ierr;

    // is cut-cell data present?
    if (fread(&iflag, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
    
    if (iflag){
      ierr = xf_Error(xf_CreateCutElemData(Mesh->ElemGroup[i].nElem, 
					   &Mesh->ElemGroup[i].CutElemData));
      if (ierr != xf_OK) return ierr;
    }
    else
      Mesh->ElemGroup[i].CutElemData = NULL;

    for (j=0; j<Mesh->ElemGroup[i].nElem; j++){
      
      for (k=0; k<Mesh->ElemGroup[i].nFace[j]; k++){
	ierr = xf_Error(xf_ReadFaceDataBinary(fid, Mesh->ElemGroup[i].Face[j]+k));
	if (ierr != xf_OK) return ierr;
      }

      if (fread(Mesh->ElemGroup[i].Node[j], si, Mesh->ElemGroup[i].nNode, fid) != 
	  Mesh->ElemGroup[i].nNode) 
	return xf_Error(xf_FILE_READ_ERROR);

      if (iflag){
	ierr = xf_Error(xf_ReadCutElemDataBinary(fid, Mesh->ElemGroup[i].CutElemData+j));
	if (ierr != xf_OK) return ierr;
      }
    }
  }


  // read periodic groups
  if (fread(&Mesh->nPeriodicGroup, si, 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  if (Mesh->nPeriodicGroup > 0){
    ierr = xf_Error(xf_Alloc((void **) &Mesh->PeriodicGroup, Mesh->nPeriodicGroup,
			     sizeof(xf_PeriodicGroup)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<Mesh->nPeriodicGroup; i++){
      PG = Mesh->PeriodicGroup+i;
      ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_PeriodicityName, xfe_PeriodicityLast, 
					(int *) &PG->Periodicity));
      if (ierr != xf_OK) return ierr;
      if (fread(&PG->nPeriodicNode, si, 1, fid) != 1) 
	return xf_Error(xf_FILE_READ_ERROR);
      ierr = xf_Error(xf_Alloc((void **) &PG->PeriodicNode, 2*PG->nPeriodicNode,
			       sizeof(int)));
      if (ierr != xf_OK) return ierr;
      if (fread(PG->PeriodicNode , si, 2*PG->nPeriodicNode, fid) != 2*PG->nPeriodicNode) 
	return xf_Error(xf_FILE_READ_ERROR);
    }
  }
    
  // read background mesh
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_CreateMesh((xf_Mesh **) &Mesh->BackgroundMesh));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadMeshBinarySerial(fid, (xf_Mesh *) Mesh->BackgroundMesh));
    if (ierr != xf_OK) return ierr;
  }
  else
    Mesh->BackgroundMesh = NULL;

  // read mesh motion
  if (rev >= 2){
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
    if (ierr != xf_OK) return ierr;
    if (flag){
      ierr = xf_Error(xf_CreateMeshMotion(&Mesh->Motion));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ReadMeshMotionBinarySerial(fid, (xf_MeshMotion *) Mesh->Motion));
      if (ierr != xf_OK) return ierr;
    }
    else
      Mesh->Motion = NULL;
  }

  return xf_OK;
}



/*-----------------*/
/* Parallelization */
/*-----------------*/

#include "xf_MeshParallel.c"


/******************************************************************/
//   FUNCTION Definition: xf_WriteMeshBinary
int 
xf_WriteMeshBinary( xf_Mesh *Mesh, FILE *fid)
{
  int myRank, nProc;
  enum xfe_Bool ParallelFlag;
  int ierr, terr;
  xf_MeshMotion *Motion = NULL;
  xf_Mesh *Mesh_Glob = NULL;
  
  // Get myRank
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if (Mesh->ParallelInfo != NULL)
    ParallelFlag = xfe_True;
  else
    ParallelFlag = xfe_False;
  
  // unparallelize Mesh -> Mesh_Glob
  if (ParallelFlag){
    if (myRank == 0){    
      ierr = xf_Error(xf_CreateMesh(&Mesh_Glob));
      if (ierr != xf_OK) return ierr;
    }

    ierr = xf_Error(xf_UnParallelizeMesh(Mesh, Mesh_Glob));
    if (ierr != xf_OK) return ierr;
  }
  else
    Mesh_Glob = Mesh;
    
  // root writes Mesh_Glob
  if (myRank == 0)
    terr = xf_Error(xf_WriteMeshBinarySerial(Mesh_Glob, fid));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  /* during unparallelization, some structures were swapped
     pointer-wise.  Swap these back and destroy Mesh_Glob. */
  if ((ParallelFlag) && (myRank == 0)){
    swap(Mesh_Glob->Motion, Mesh->Motion, Motion);
    ierr = xf_Error(xf_DestroyMesh(Mesh_Glob));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadMeshBinary
int 
xf_ReadMeshBinary( FILE *fid, xf_Mesh *Mesh)
{
  int myRank;
  int ierr, terr;
  xf_Mesh *Mesh_Glob = NULL;
  
  // Get myRank
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  if (myRank == 0){    
    terr = xf_Error(xf_CreateMesh(&Mesh_Glob));
    if (terr == xf_OK)
      terr = xf_Error(xf_ReadMeshBinarySerial(fid, Mesh_Glob));
  }
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  ierr = xf_Error(xf_ParallelizeMesh(Mesh_Glob, Mesh, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyMesh(Mesh_Glob));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

