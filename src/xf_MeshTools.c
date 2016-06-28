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
 FILE:  xf_MeshTools.c
 
 This file contains functions for performing calculations with the Mesh
 geometry or connectivity.
 
 */

#include <stdlib.h>

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_Basis.h"
#include "xf_BasisFcn.h"
#include "xf_Memory.h"
#include "xf_Data.h"
#include "xf_Math.h"
#include "xf_Param.h"
#include "xf_MathLapack.h"
#include "xf_Quad.h"
#include "xf_MeshToolsStruct.h"
#include "xf_MPI.h"

/* Hash list for faces */
typedef struct
{
  int n0, n1, n2;
  int next;
}
xf_FaceHash;

// tetrahedron edge info
const int E2N[6][2] = {{0,1}, {0,2}, {1,2}, {0,3}, {1,3}, {2,3}}; // edge-to-node
const int OE[6] = {5, 4, 3, 2, 1, 0}; // opposite edge



/******************************************************************/
//   FUNCTION Definition: xf_GetnElem
int 
xf_GetnElem(xf_Mesh *Mesh, int **nElem, int *nelemtot)
{
  int ierr, nelem, egrp;
  
  if (nElem != NULL){
    ierr = xf_Error(xf_Alloc( (void **) nElem, Mesh->nElemGroup, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  if (nelemtot != NULL) (*nelemtot) = 0;
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    nelem = Mesh->ElemGroup[egrp].nElem;
    if (nelemtot != NULL) (*nelemtot) += nelem;
    if (nElem    != NULL) (*nElem)[egrp] = nelem;
  } // egrp
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_EgrpElem2Index
int 
xf_EgrpElem2Index(xf_Mesh *Mesh, int egrp0, int elem0, int *pie)
{
  int egrp, ie;

  if ((egrp0<0) || (egrp0>Mesh->nElemGroup)) return xf_Error(xf_OUT_OF_BOUNDS);

  ie=elem0;
  for (egrp=0; egrp<egrp0; egrp++) ie += Mesh->ElemGroup[egrp].nElem;
  
  (*pie) = ie;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Index2EgrpElem
int 
xf_Index2EgrpElem(xf_Mesh *Mesh, int ie0, int *pegrp, int *pelem)
{
  int egrp, elem, ie, nelem;

  if (ie<0) return xf_Error(xf_OUT_OF_BOUNDS);

  for (egrp=0, ie=0; egrp<Mesh->nElemGroup; egrp++){
    nelem = Mesh->ElemGroup[egrp].nElem;
    if (ie+nelem > ie0){
      elem = ie0-ie;
      break;
    }
    else if (egrp == Mesh->nElemGroup-1) return xf_Error(xf_OUT_OF_BOUNDS);
    ie += Mesh->ElemGroup[egrp].nElem;
  }
  
  (*pegrp) = egrp;
  (*pelem) = elem;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_NeighborAcrossFace
int 
xf_NeighborAcrossFace(xf_Mesh *Mesh, int egL, int eL, int faceL, 
                      int *egR, int *eR, int *faceR)
{
  
  int iiface;
  xf_Face Face;
  
  Face = Mesh->ElemGroup[egL].Face[eL][faceL];
  
  if ( Face.Group == xf_INTERIORFACE){
    iiface = Face.Number;
    *egR   = Mesh->IFace[iiface].ElemGroupR;
    *eR    = Mesh->IFace[iiface].ElemR;
    if (faceR != NULL) *faceR = Mesh->IFace[iiface].FaceR;
    
    if ((*egR == egL) && (*eR == eL)){
      *egR   = Mesh->IFace[iiface].ElemGroupL;
      *eR    = Mesh->IFace[iiface].ElemL;
      if (faceR != NULL) *faceR = Mesh->IFace[iiface].FaceL;
    }
  }
  else{ // note, xf_NULLFACE falls into this category
    *eR    = -1;
    *egR   = -1;
    if (faceR != NULL) *faceR = -1;
  }
  
  return xf_OK; 
}



/******************************************************************/
//   FUNCTION Definition: xf_CommonFace
int 
xf_CommonFace(xf_Mesh *Mesh, int egL, int eL, int egR, int eR,
              int *pfaceL, int *pfaceR)
{
  int faceL, faceR;
  int nfaceL, nfaceR;
  xf_Face *FaceL, *FaceR;
  
  nfaceL = Mesh->ElemGroup[egL].nFace[eL];
  nfaceR = Mesh->ElemGroup[egR].nFace[eR];
  FaceL  = Mesh->ElemGroup[egL].Face[eL];
  FaceR  = Mesh->ElemGroup[egR].Face[eR];
  
  for (faceL=0; faceL<nfaceL; faceL++){
    if (FaceL[faceL].Group != xf_INTERIORFACE) continue;
    for (faceR=0; faceR<nfaceR; faceR++){
      if (FaceR[faceR].Group != xf_INTERIORFACE) continue;
      if (FaceR[faceR].Number == FaceL[faceL].Number){
        (*pfaceL) = faceL;
        (*pfaceR) = faceR;
        return xf_OK;
      }
    } // faceR
  } // faceL
  return xf_NOT_FOUND;
  
}



/******************************************************************/
//   FUNCTION Definition: xf_IsElemOnLeft
int 
xf_IsElemOnLeft(xf_IFace IFace, int egrp, int elem, enum xfe_Bool *pIamL)
{
  if ((IFace.ElemGroupL == egrp) && (IFace.ElemL == elem))
    (*pIamL) = xfe_True;
  else if ((IFace.ElemGroupR == egrp) && (IFace.ElemR == elem))
    (*pIamL) = xfe_False;
  else return xf_Error(xf_MESH_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FaceElements
void 
xf_FaceElements(xf_Mesh *Mesh, int ibfgrp, int ibface, 
                int *egrpL, int *elemL, int *faceL, 
                int *egrpR, int *elemR, int *faceR)
{
  xf_IFace IFace;
  xf_BFace BFace;
  
  if (ibfgrp == -1){ // Interior face
    IFace = Mesh->IFace[ibface];
    if (egrpL != NULL) (*egrpL) = IFace.ElemGroupL;
    if (egrpR != NULL) (*egrpR) = IFace.ElemGroupR;
    if (elemL != NULL) (*elemL) = IFace.ElemL;
    if (elemR != NULL) (*elemR) = IFace.ElemR;
    if (faceL != NULL) (*faceL) = IFace.FaceL;
    if (faceR != NULL) (*faceR) = IFace.FaceR;
  }
  else{             // Boundary face
    BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
    if (egrpL != NULL) (*egrpL) =  BFace.ElemGroup;
    if (egrpR != NULL) (*egrpR) = -1;
    if (elemL != NULL) (*elemL) =  BFace.Elem;
    if (elemR != NULL) (*elemR) = -1;
    if (faceL != NULL) (*faceL) = BFace.Face;
    if (faceR != NULL) (*faceR) = -1;
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_GetFaceOrient
int 
xf_GetFaceOrient(xf_Mesh *Mesh, int egrp, int elem, int face, 
                 int *porient)
{
  xf_Face Face;
  xf_IFace *IFace;
  
  Face = Mesh->ElemGroup[egrp].Face[elem][face];
  
  if ( Face.Group == xf_INTERIORFACE){ // interior
    IFace = Mesh->IFace + Face.Number;
    if ((IFace->ElemGroupL == egrp) &&
        (IFace->ElemL      == elem))
      (*porient) = IFace->OrientL;
    else if ((IFace->ElemGroupR == egrp) &&
             (IFace->ElemR      == elem))
      (*porient) = IFace->OrientR;
    else
      return xf_Error(xf_MESH_ERROR);
  }
  else if ( Face.Group >= 0){
    (*porient) = Mesh->BFaceGroup[Face.Group].BFace[Face.Number].Orient;
  }
  else return xf_Error(xf_NOT_SUPPORTED);
  
  return xf_OK; 
}


/******************************************************************/
//   FUNCTION Definition: xf_DetermineFaceOrient_Segment
static int 
xf_DetermineFaceOrient_Segment(int *gnode, int *pfaceorient)
{
  (*pfaceorient) = 0;
  return xf_OK; 
}

/******************************************************************/
//   FUNCTION Definition: xf_DetermineFaceOrient_Triangle
static int 
xf_DetermineFaceOrient_Triangle(int *gnode, int *pfaceorient)
{
  (*pfaceorient) = (gnode[1] < gnode[0]);
  return xf_OK; 
}

/******************************************************************/
//   FUNCTION Definition: xf_DetermineFaceOrient_Quadrilateral
static int 
xf_DetermineFaceOrient_Quadrilateral(int *gnode, int *pfaceorient)
{
  (*pfaceorient) = (gnode[1] < gnode[0]);
  return xf_OK; 
}

/******************************************************************/
//   FUNCTION Definition: xf_DetermineFaceOrient_Tetrahedron
static int 
xf_DetermineFaceOrient_Tetrahedron(int *gnode, int *pfaceorient)
{
  int i0, i1, i2, itemp;
  
  // index of max node
  i2 = 0; 
  if (gnode[ 1] > gnode[i2]) i2 = 1;
  if (gnode[ 2] > gnode[i2]) i2 = 2;
  
  // other 2 indices
  i0 = (i2+1)%3; 
  i1 = (i2+2)%3;
  
  // i0 becomes index of min, i1 is next
  if (gnode[i0] > gnode[i1]) swap(i0,i1, itemp);
  
  (*pfaceorient) = i0;
  if (i1 != ((i0+1)%3)) (*pfaceorient) += 3;
  
  return xf_OK; 
}



/******************************************************************/
//   FUNCTION Definition: xf_DetermineFaceOrient_Hexahedron
static int 
xf_DetermineFaceOrient_Hexahedron(int *gnode, int *pfaceorient)
{
  int i0;
  
  // index of min node
  i0 = 0; 
  if (gnode[1] < gnode[i0]) i0 = 1;
  if (gnode[2] < gnode[i0]) i0 = 2;
  if (gnode[3] < gnode[i0]) i0 = 3;
  
  // first piece of info is loc of min node (4 possibilities)
  (*pfaceorient) = i0;
  
  // next piece of info is whether i0+1 node is > or < i0-1 node
  if (gnode[(i0+1)%4] > gnode[(i0+3)%4]) (*pfaceorient) += 4;
  
  return xf_OK; 
}


/******************************************************************/
//   FUNCTION Definition: xf_DetermineFaceOrient_Shape
static int 
xf_DetermineFaceOrient_Shape( enum xfe_ShapeType Shape, int *nvec, int *pfaceorient)
{
  switch (Shape){
    case xfe_Segment:   
      return xf_DetermineFaceOrient_Segment(nvec, pfaceorient);
      break;
    case xfe_Triangle:   
      return xf_DetermineFaceOrient_Triangle(nvec, pfaceorient);
      break;
    case xfe_Tetrahedron:  
      return xf_DetermineFaceOrient_Tetrahedron(nvec, pfaceorient);
      break;
    case xfe_Quadrilateral:   
      return xf_DetermineFaceOrient_Quadrilateral(nvec, pfaceorient);
      break;
    case xfe_Hexahedron:  
      return xf_DetermineFaceOrient_Hexahedron(nvec, pfaceorient);
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DetermineFaceOrient
int 
xf_DetermineFaceOrient(xf_Mesh *Mesh, int *NodeMap, int egrp, int elem, 
		       int face, int *pfaceorient)
{
  int ierr, i, QOrder, nfnode;
  int *Node;
  int nvec[xf_MAXQ1FACENODE];
  enum xfe_ShapeType Shape;
  enum xfe_BasisType QBasis;
  
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  
  // determine Shape of element
  ierr = xf_Error(xf_Basis2Shape(QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  Node = Mesh->ElemGroup[egrp].Node[elem];
  
  // Q1 nodes on face
  ierr = xf_Error(xf_Q1NodesOnFace(QBasis, QOrder, face, &nfnode, nvec));
  if (ierr != xf_OK) return ierr;
  
  if (nfnode > xf_MAXQ1FACENODE) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // convert to global node numbers
  for (i=0; i<nfnode; i++) nvec[i] = Node[nvec[i]];
  
  // use a node map if one is provided
  if (NodeMap != NULL)
    for (i=0; i<nfnode; i++) nvec[i] = NodeMap[nvec[i]];
  
  // call generic function
  ierr = xf_Error(xf_DetermineFaceOrient_Shape(Shape, nvec, pfaceorient));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReOrderHangNodes
static int 
xf_ReOrderHangNodes(enum xfe_ShapeType FShape, int nfnode0, 
                    int *fvec0, int nfine, int *posfine, int *nvfine,
                    int *pnedge, int *elist)
{
  /*
   PURPOSE:
   
   Reorders a set of node lists in nvfine to be consistent with the
   viewpoint from the coarse side of a hanging face, where each node
   list contains global node numbers on a subface of a hanging face.
   
   INPUTS:
   
   FShape : shape of face to consider
   nfnode0 : number of Q1 nodes on the face
   fvec0 : nfnode0 Q1 nodes as seen from coarse hang face, in desired order
   nfine : number of subfaces
   posfine : positions of subfaces
   nvfine : global node numbers on each subface ([nfnode0*nfine], unrolled),
   in some arbitrary order on each subface.
   
   OUTPUTS: 
   
   nvfine : sorted as described above
   (*pnedge) : number of bisected edges (optional)
   elist : list of three-tuples corresponding to global node numbers of
   bisected edges, e.g. [n0 n1 nm ... ], where nm is the global 
   node number of the node between n0 and n1.
   
   RETURN:
   
   Error Code
   */
  
  int i, j;
  int ipos1, ipos2;
  int t;
  int fvec[9];
  int *pvfine[xfe_QuadPosLast];
  int *pvf;
  int QU[4][4] = {{0,1,4,3}, {1,2,5,4}, {3,4,7,6}, {4,5,8,7}};
  
  switch (FShape){
    case xfe_Segment:
      if (nfine != 2) return xf_Error(xf_INPUT_ERROR);
      if (posfine[0] == xfe_SegPosLeft) ipos1 = 0;
      else if (posfine[1] == xfe_SegPosLeft) ipos1 = 1;
      else return xf_Error(xf_INPUT_ERROR);
      ipos2 = 1-ipos1;
      
      /* xf_printf("fvec0 = %d %d\n", fvec0[0], fvec0[1]); */
      /*     xf_printf("ipos1=%d: nvfine = %d %d\n", ipos1, nvfine[ipos1*nfnode0+0], */
      /* 	      nvfine[ipos1*nfnode0+1]); */
      /*     xf_printf("ipos2=%d: nvfine = %d %d\n", ipos2, nvfine[ipos2*nfnode0+0], */
      /* 	      nvfine[ipos2*nfnode0+1]); */
      
      if (nvfine[ipos1*nfnode0+1] == fvec0[0]){
        swap(nvfine[ipos1*nfnode0+0], nvfine[ipos1*nfnode0+1], t);
      }
      else if (nvfine[ipos1*nfnode0+0] != fvec0[0])
        return xf_Error(xf_INPUT_ERROR);
      
      if (nvfine[ipos2*nfnode0+0] == fvec0[1]){
        swap(nvfine[ipos2*nfnode0+0], nvfine[ipos2*nfnode0+1], t);
      }
      else if (nvfine[ipos2*nfnode0+1] != fvec0[1])
        return xf_Error(xf_INPUT_ERROR);
      
      if (pnedge != NULL) return xf_Error(xf_INPUT_ERROR);
      
      break;
      case xfe_Quadrilateral:
      if (nfnode0 != 4) return xf_Error(xf_INPUT_ERROR);
      // set fvec to -1s
      for (i=0; i<9; i++) fvec[i] = -1;
      
      /*     for (i=0; i<4; i++) xf_printf("fvec0[%d] = %d\n", i, fvec0[i]); */
      /*     for (i=0;i<nfine; i++) */
      /*       for (j=0; j<4; j++) xf_printf("nvfine[%d][%d] = %d\n", i, j, nvfine[i*4+j]); */
      
      if (posfine[0] <= xfe_QuadPosNE){ // uniform
        //if (nfine != 4) return xf_Error(xf_INPUT_ERROR); 
        // corners
        fvec[0] = fvec0[0]; fvec[2] = fvec0[1]; fvec[6] = fvec0[3]; fvec[8] = fvec0[2];
        // initialize pvfine to NULL
        for (i=0; i<xfe_QuadPosLast; i++) pvfine[i] = NULL;
        // order the positions
        for (i=0; i<nfine; i++) pvfine[posfine[i]] = nvfine + i*nfnode0;
        // check that all positions are present
        /*       if ((pvfine[xfe_QuadPosSW] == NULL) || (pvfine[xfe_QuadPosSE] == NULL) || */
        /* 	  (pvfine[xfe_QuadPosNW] == NULL) || (pvfine[xfe_QuadPosNE] == NULL)) */
        /* 	return xf_Error(xf_INPUT_ERROR); */
        // get center node, 4
        if ((pvfine[xfe_QuadPosSW] != NULL) && (pvfine[xfe_QuadPosNE] != NULL)){
          for (i=0; i<16; i++) 
            if (pvfine[xfe_QuadPosSW][i/4] == pvfine[xfe_QuadPosNE][i%4]){
              fvec[4] = pvfine[xfe_QuadPosSW][i/4];
              break;
            }
        }
        else if ((pvfine[xfe_QuadPosSE] != NULL) && (pvfine[xfe_QuadPosNW] != NULL)){
          for (i=0; i<16; i++) 
            if (pvfine[xfe_QuadPosSE][i/4] == pvfine[xfe_QuadPosNW][i%4]){
              fvec[4] = pvfine[xfe_QuadPosSE][i/4];
              break;
            }
        }
        
        if ((pvfine[xfe_QuadPosSW] != NULL) && (pvfine[xfe_QuadPosSE] != NULL))
          for (i=0; i<16; i++) // edge node, 1
            if (pvfine[xfe_QuadPosSW][i/4] == pvfine[xfe_QuadPosSE][i%4])
              if ((j=pvfine[xfe_QuadPosSW][i/4]) != fvec[4]){
                fvec[1] = j;
                break;
              }
        
        if ((pvfine[xfe_QuadPosSE] != NULL) && (pvfine[xfe_QuadPosNE] != NULL))
          for (i=0; i<16; i++) // edge node, 5
            if (pvfine[xfe_QuadPosSE][i/4] == pvfine[xfe_QuadPosNE][i%4])
              if ((j=pvfine[xfe_QuadPosSE][i/4]) != fvec[4]){
                fvec[5] = j;
                break;
              }
        
        if ((pvfine[xfe_QuadPosNE] != NULL) && (pvfine[xfe_QuadPosNW] != NULL))
          for (i=0; i<16; i++) // edge node, 7
            if (pvfine[xfe_QuadPosNE][i/4] == pvfine[xfe_QuadPosNW][i%4])
              if ((j=pvfine[xfe_QuadPosNE][i/4]) != fvec[4]){
                fvec[7] = j;
                break;
              }
        
        if ((pvfine[xfe_QuadPosNW] != NULL) && (pvfine[xfe_QuadPosSW] != NULL))
          for (i=0; i<16; i++) // edge node, 3
            if (pvfine[xfe_QuadPosNW][i/4] == pvfine[xfe_QuadPosSW][i%4])
              if ((j=pvfine[xfe_QuadPosNW][i/4]) != fvec[4]){
                fvec[3] = j;
                break;
              }
        
        // check if we missed any nodes
        //for (i=0; i<9; i++) if (fvec[i] < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
        
        // reassign nvfine via pvfine
        for (i=xfe_QuadPosSW; i<=xfe_QuadPosNE; i++)
          if (pvfine[i] != NULL)
            for (j=0; j<4; j++)
              pvfine[i][j] = fvec[QU[i-xfe_QuadPosSW][j]];
        
        // set edge list
        if (pnedge != NULL){
          j = 0;
          if (fvec[1] >= 0){
            elist[3*j+0] = fvec[0]; elist[3*j+1] = fvec[2]; elist[3*j+2] = fvec[1];
            j++;
          }
          if (fvec[5] >= 0){
            elist[3*j+0] = fvec[2]; elist[3*j+1] = fvec[8]; elist[3*j+2] = fvec[5];
            j++;
          }
          if (fvec[7] >= 0){
            elist[3*j+0] = fvec[8]; elist[3*j+1] = fvec[6]; elist[3*j+2] = fvec[7];
            j++;
          }
          if (fvec[3] >= 0){
            elist[3*j+0] = fvec[6]; elist[3*j+1] = fvec[0]; elist[3*j+2] = fvec[3];
            j++;
          }
          (*pnedge) = j;
        }
        
      }
      else if (posfine[0] <= xfe_QuadPosRight){ // left/right
        if (nfine != 2) return xf_Error(xf_INPUT_ERROR); 
        // corners
        fvec[0] = fvec0[0]; fvec[2] = fvec0[1]; fvec[3] = fvec0[3]; fvec[5] = fvec0[2];
        // initialize pvfine to NULL
        for (i=0; i<xfe_QuadPosLast; i++) pvfine[i] = NULL;
        // order the positions
        for (i=0; i<nfine; i++) pvfine[posfine[i]] = nvfine + i*nfnode0;
        // check that all positions are present
        if ((pvfine[xfe_QuadPosLeft] == NULL) || (pvfine[xfe_QuadPosRight] == NULL))
          return xf_Error(xf_INPUT_ERROR);
        // use Left to set fvec
        pvf = pvfine[xfe_QuadPosLeft];
        for (i=0; i<4; i++){
          if ((pvf[i]==fvec[0]) && (pvf[(i+3)%4]==fvec[3])){
            fvec[1]=pvf[(i+1)%4];
            fvec[4]=pvf[(i+2)%4];
            break;
          }
          if ((pvf[i]==fvec[0]) && (pvf[(i+1)%4]==fvec[3])){
            fvec[1]=pvf[(i+3)%4];
            fvec[4]=pvf[(i+2)%4];
            break;
          }
        }
        // check against Right
        pvf = pvfine[xfe_QuadPosRight];
        for (i=0; i<4; i++){
          if ((pvf[i]==fvec[2]) && (pvf[(i+3)%4]==fvec[5])){
            if ((fvec[1]!=pvf[(i+1)%4]) || (fvec[4]!=pvf[(i+2)%4]))
              return xf_Error(xf_CODE_LOGIC_ERROR);
          }
          if ((pvf[i]==fvec[2]) && (pvf[(i+1)%4]==fvec[5])){
            if ((fvec[1]!=pvf[(i+3)%4]) || (fvec[4]!=pvf[(i+2)%4]))
              return xf_Error(xf_CODE_LOGIC_ERROR);
          }
        }
        
        // check if we missed any nodes
        for (i=0; i<6; i++) if (fvec[i] < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
        // reassign nvfine via pvfine
        pvf = pvfine[xfe_QuadPosLeft];
        pvf[0] = fvec[0]; pvf[1] = fvec[1]; pvf[2] = fvec[4]; pvf[3] = fvec[3];
        pvf = pvfine[xfe_QuadPosRight];
        pvf[0] = fvec[1]; pvf[1] = fvec[2]; pvf[2] = fvec[5]; pvf[3] = fvec[4];
        
        // set edge list
        if (pnedge != NULL){
          (*pnedge) = 2;
          elist[3*0+0] = fvec[0]; elist[3*0+1] = fvec[2]; elist[3*0+2] = fvec[1];
          elist[3*1+0] = fvec[5]; elist[3*1+1] = fvec[3]; elist[3*1+2] = fvec[4];
        }
        
      }
      else{ // bottom/top
        if (nfine != 2) return xf_Error(xf_INPUT_ERROR); 
        // corners
        fvec[0] = fvec0[0]; fvec[1] = fvec0[1]; fvec[4] = fvec0[3]; fvec[5] = fvec0[2];
        // initialize pvfine to NULL
        for (i=0; i<xfe_QuadPosLast; i++) pvfine[i] = NULL;
        // order the positions
        for (i=0; i<nfine; i++) pvfine[posfine[i]] = nvfine + i*nfnode0;
        // check that all positions are present
        if ((pvfine[xfe_QuadPosBottom] == NULL) || (pvfine[xfe_QuadPosTop] == NULL))
          return xf_Error(xf_INPUT_ERROR);
        // use Bottom to set fvec
        pvf = pvfine[xfe_QuadPosBottom];
        for (i=0; i<4; i++){
          if ((pvf[i]==fvec[0]) && (pvf[(i+1)%4]==fvec[1])){
            fvec[2]=pvf[(i+3)%4];
            fvec[3]=pvf[(i+2)%4];
            break;
          }
          if ((pvf[i]==fvec[0]) && (pvf[(i+3)%4]==fvec[1])){
            fvec[2]=pvf[(i+1)%4];
            fvec[3]=pvf[(i+2)%4];
            break;
          }
        }
        // check against Top
        pvf = pvfine[xfe_QuadPosTop];
        for (i=0; i<4; i++){
          if ((pvf[i]==fvec[4]) && (pvf[(i+1)%4]==fvec[5])){
            if ((fvec[2]!=pvf[(i+3)%4]) || (fvec[3]!=pvf[(i+2)%4]))
              return xf_Error(xf_CODE_LOGIC_ERROR);
          }
          if ((pvf[i]==fvec[4]) && (pvf[(i+3)%4]==fvec[5])){
            if ((fvec[2]!=pvf[(i+1)%4]) || (fvec[3]!=pvf[(i+2)%4]))
              return xf_Error(xf_CODE_LOGIC_ERROR);
          }
        }
        
        // check if we missed any nodes
        for (i=0; i<6; i++) if (fvec[i] < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
        // reassign nvfine via pvfine
        pvf = pvfine[xfe_QuadPosBottom];
        pvf[0] = fvec[0]; pvf[1] = fvec[1]; pvf[2] = fvec[3]; pvf[3] = fvec[2];
        pvf = pvfine[xfe_QuadPosTop];
        pvf[0] = fvec[2]; pvf[1] = fvec[3]; pvf[2] = fvec[5]; pvf[3] = fvec[4];
        
        // set edge list
        if (pnedge != NULL){
          (*pnedge) = 2;
          elist[3*0+0] = fvec[1]; elist[3*0+1] = fvec[5]; elist[3*0+2] = fvec[3];
          elist[3*1+0] = fvec[4]; elist[3*1+1] = fvec[0]; elist[3*1+2] = fvec[2];
        }
        
      }
      
      break;
      case xfe_Triangle:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
      case xfe_Tetrahedron:
      case xfe_Hexahedron:
      return xf_Error(xf_INPUT_ERROR); // not a face!
      break;
      default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CheckHangFace
int 
xf_CheckHangFace(xf_Mesh *Mesh, int egrp, int elem, int face, 
                 int *phang, int *pface0, int *pfacepos, int *pelempos)
{
  int ierr, nface0, hang, face0, fpos, epos;
  enum xfe_ShapeType Shape;
  xf_Face *Face;
  xf_IFace *IFace;
  
  hang  = 0;
  face0 = face;
  fpos  = 0;
  epos  = 0;
  
  Face = Mesh->ElemGroup[egrp].Face[elem]+face;
  
  if ((Face->Group == xf_INTERIORFACE) && 
      ((hang = Mesh->IFace[Face->Number].HangNumber) != 0)){
    IFace = Mesh->IFace + Face->Number;
    //this mean the coarse side of a hang face
    if ( ((hang > 0) && (IFace->ElemGroupR == egrp) && (IFace->ElemR == elem)) ||
        ((hang < 0) && (IFace->ElemGroupL == egrp) && (IFace->ElemL == elem)) ){
      if (hang < 0) hang = -hang;
      
      // number of faces in the original group
      ierr = xf_Error(xf_Basis2nFace(Mesh->ElemGroup[egrp].QBasis, &nface0));
      if (ierr != xf_OK) return ierr;
      
      face0  = hang % nface0; // original face is stored in modulus of hang
      fpos   = hang / nface0; // position of fine face w.r.t coarse face
      
      if (pelempos != NULL){
        // shape of elem
        ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
        if (ierr != xf_OK) return ierr;
        
        // determine elem pos using face0 and fpos
        ierr = xf_Error(xf_FacePos2ElemPos(Shape, face0, fpos, &epos));
        if (ierr != xf_OK) return ierr;
      }
      
    }
    else
      hang = 0; // fine side of a hang face is not treated specially
  }
  
  // set requested values
  if (phang    != NULL) (*phang) = hang;
  if (pface0   != NULL) (*pface0) = face0;
  if (pfacepos != NULL) (*pfacepos)  = fpos;
  if (pelempos != NULL) (*pelempos)  = epos;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetCoarseOrients
int 
xf_SetCoarseOrients(xf_Mesh *Mesh, int egrp, int elem, int face0, 
                    enum xfe_Bool SetFlag, int *pnedge, int *elist)
{
  int ierr, i, nface, nfine, face;
  enum xfe_BasisType QBasis, QBasisN;
  int QOrder, QOrderN, nfnode0;
  int *posfine, *facevec, *nvfine;
  int orient, orientN, faceN, nfnode;
  int hang, f0, pos, k, egN, eN, ifine;
  enum xfe_Bool IamL;
  enum xfe_ShapeType Shape, FShape;
  xf_Face Face;
  xf_IFace *IFace;
  int fvec[xf_MAXQ1NODE], fvec0[xf_MAXQ1NODE];
  
  // number of faces
  nface = Mesh->ElemGroup[egrp].nFace[elem];
  
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  
  // get Q1 nodes on face0
  ierr = xf_Error(xf_Q1NodesOnFace(QBasis, QOrder, face0, &nfnode0, fvec0));
  if (ierr != xf_OK) return ierr;
  
  // convert to global node numbers
  for (i=0; i<nfnode0; i++) fvec0[i] = Mesh->ElemGroup[egrp].Node[elem][fvec0[i]];
  
  nfine = 0;
  
  ierr = xf_Error(xf_Alloc( (void **) &posfine, nface, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &facevec, nface, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &nvfine, nfnode0*nface, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  for (face=0; face<nface; face++){
    
    Face = Mesh->ElemGroup[egrp].Face[elem][face];
    
    /*     xf_printf("  In SCO: egrp=%d, elem=%d, face=%d, Face.Group = %d, Face.Number=%d\n", */
    /* 	      egrp, elem, face, Face.Group, Face.Number); */
    
    if (Face.Group != xf_INTERIORFACE) continue;
    
    IFace = Mesh->IFace + Face.Number;
    if ((IFace->ElemGroupL == egrp) &&
        (IFace->ElemL      == elem)){
      IamL = xfe_True;
      egN = IFace->ElemGroupR;
      eN  = IFace->ElemR;
      QBasisN = Mesh->ElemGroup[IFace->ElemGroupR].QBasis;
      QOrderN = Mesh->ElemGroup[IFace->ElemGroupR].QOrder;
      faceN   = IFace->FaceR;
      orientN = IFace->OrientR;
      orient = IFace->OrientL;
    }
    else if ((IFace->ElemGroupR == egrp) &&
             (IFace->ElemR      == elem)){
      IamL = xfe_False;
      egN = IFace->ElemGroupL;
      eN  = IFace->ElemL;
      QBasisN = Mesh->ElemGroup[IFace->ElemGroupL].QBasis;
      QOrderN = Mesh->ElemGroup[IFace->ElemGroupL].QOrder;
      faceN   = IFace->FaceL;
      orientN = IFace->OrientL;
      orient = IFace->OrientR;
    }
    else
      return xf_Error(xf_MESH_ERROR);
    
    hang = IFace->HangNumber;
    
    // get original face
    ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, &f0, &pos, NULL));
    if (ierr != xf_OK) return ierr;
    
    if (f0 != face0) continue; // not the face0 we're interested in
    
    /*  xf_printf(" In SetCoarseOrients: egrp = %d, elem=%d, face=%d, face0=%d, pos=%d\n", */
    /*     	      egrp, elem, face, face0, pos); */
    /*  xf_printf("   egN = %d, eN = %d, faceN = %d, orientN = %d, Face.Number = %d\n", */
    /*  	      egN, eN, faceN, orientN, Face.Number); */
    
    // get nodes on fine face, in local numbering of neighbor element
    ierr = xf_Error(xf_Q1NodesOnFaceNeighbor(QBasisN, QOrderN, faceN, orientN, 
                                             QBasis, QOrder, face0, orient, 
                                             &nfnode, fvec));
    if (ierr != xf_OK) return ierr;
    
    if (nfnode != nfnode0) return xf_Error(xf_MESH_ERROR);
    
    posfine[nfine] = pos;
    facevec[nfine] = Face.Number;
    for (k=0; k<nfnode0; k++) // convert to global numbers
      nvfine[nfine*nfnode0+k] = Mesh->ElemGroup[egN].Node[eN][fvec[k]];
    nfine++;
    
  } // face
  
  if (nfine <= 0) return xf_Error(xf_INPUT_ERROR);
  
  // get element Shape
  ierr = xf_Error(xf_Basis2Shape(QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  // get face shape on face0
  ierr = xf_Error(xf_FaceShape(Shape, face0, &FShape));
  if (ierr != xf_OK) return ierr;
  
  // reorder fine face nodes to a standard pattern (i.e. point of view
  // of coarse elem), using fvec0
  ierr = xf_Error(xf_ReOrderHangNodes(FShape, nfnode0, fvec0, nfine, 
                                      posfine, nvfine, pnedge, elist));
  if (ierr != xf_OK) return ierr;
  
  // loop over each fine face and determine the orientation
  if (SetFlag){
    for (ifine=0; ifine<nfine; ifine++){
      ierr = xf_Error(xf_DetermineFaceOrient_Shape(Shape, nvfine+ifine*nfnode0, &orient));
      if (ierr != xf_OK) return ierr;
      
      IFace = Mesh->IFace + facevec[ifine];
      
      IamL = ((IFace->ElemGroupL == egrp) && (IFace->ElemL == elem));
      
      if (IamL){ // set orient and also determine orient of neighbor
        IFace->OrientL = orient;
        //xf_printf(" neighbor = %d, %d\n", IFace->ElemGroupR, IFace->ElemR);
        ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NULL, IFace->ElemGroupR, IFace->ElemR,
                                               IFace->FaceR, &IFace->OrientR));
        if (ierr != xf_OK) return ierr;
      }
      else{
        IFace->OrientR = orient;
        //xf_printf(" neighbor = %d, %d\n", IFace->ElemGroupL, IFace->ElemL);
        ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NULL, IFace->ElemGroupL, IFace->ElemL,
                                               IFace->FaceL, &IFace->OrientL));
        if (ierr != xf_OK) return ierr;
      }
      
      /*       xf_printf("   IamL = %d, orientL = %d, orientR =%d\n", IamL, */
      /* 		IFace->OrientL, IFace->OrientR); */
      if (IFace->OrientL == IFace->OrientR) return xf_Error(xf_NOT_SUPPORTED);
      
    } // ifine
  }
  
  xf_Release( (void *) posfine);
  xf_Release( (void *) facevec);
  xf_Release( (void *) nvfine);
  
  
  return xf_OK; 
}

/******************************************************************/
//   FUNCTION Definition: xf_MakePeriodicNodeMap
static int
xf_MakePeriodicNodeMap(xf_Mesh *Mesh, int **pNodeMap)
{
  // NodeMap[j] = minimum global node number out of two possibilities (when periodic)
  int ierr;
  int i, iPG, *N;
  int *NodeMap;
  xf_PeriodicGroup *PG;

  if (Mesh->nPeriodicGroup <= 0) return xf_OK; // nothing to do

  ierr = xf_Error(xf_Alloc((void **) pNodeMap, Mesh->nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  NodeMap = (*pNodeMap);

  for (i=0; i<Mesh->nNode; i++) NodeMap[i] = i; // default node map (identity)

  for (iPG=0; iPG<Mesh->nPeriodicGroup; iPG++){
    PG = Mesh->PeriodicGroup+iPG;
    for (i=0; i<PG->nPeriodicNode; i++){
      N = PG->PeriodicNode + 2*i;
      NodeMap[N[0]] = NodeMap[N[1]] = min(NodeMap[N[0]],NodeMap[N[1]]);
    }
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UpdateFaceOrient
int 
xf_UpdateFaceOrient(xf_Mesh *Mesh)
{
  int ierr;
  int iiface;
  int ibfgrp, ibface;
  int *NodeMap = NULL;
  xf_IFace IFace;
  xf_BFace BFace;

  // make a node map if have periodic groups
  if (Mesh->nPeriodicGroup > 0){
    ierr = xf_Error(xf_MakePeriodicNodeMap(Mesh, &NodeMap));
    if (ierr != xf_OK) return ierr;
  }
  
  for (iiface=0; iiface<Mesh->nIFace; iiface++){
    IFace = Mesh->IFace[iiface];
    
    if (IFace.HangNumber >= 0){ // Do not determine orientation on coarse hanging faces
      ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NodeMap, IFace.ElemGroupL, IFace.ElemL,
                                             IFace.FaceL, &Mesh->IFace[iiface].OrientL));
      if (ierr != xf_OK) return ierr;
    }
    
    if (IFace.HangNumber <= 0){ // Do not determine orientation on coarse hanging faces
      ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NodeMap, IFace.ElemGroupR, IFace.ElemR,
                                             IFace.FaceR, &Mesh->IFace[iiface].OrientR));
      if (ierr != xf_OK) return ierr;
    }

    if ((IFace.HangNumber == 0) && (Mesh->Dim > 1)){
      if (Mesh->IFace[iiface].OrientL ==  Mesh->IFace[iiface].OrientR)
        return xf_Error(xf_OUT_OF_BOUNDS);
    }
  } // iiface
  
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      
      BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
      
      ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NodeMap, BFace.ElemGroup, BFace.Elem, BFace.Face, 
                                             &Mesh->BFaceGroup[ibfgrp].BFace[ibface].Orient));
      if (ierr != xf_OK) return ierr;
      
    } // ibface
  } // ibfgrp
  

  xf_Release( (void *) NodeMap);
  
  return xf_OK; 
}


/******************************************************************/
//   FUNCTION Definition: xf_GetRelativeOrient
int 
xf_GetRelativeOrient(enum xfe_ShapeType FShape, int orient1, int orient2, 
                     int *orientrel)
{
  int cycle1, flip1;
  int cycle2, flip2;
  int cycle, flip;
  int sg;
  
  switch (FShape){
    case xfe_Segment:
      (*orientrel) = !(orient1==orient2);
      break;
    case xfe_Quadrilateral:
      cycle1 = orient1%4;
      flip1  = orient1/4;
      cycle2 = orient2%4;
      flip2  = orient2/4;
      
      flip   =  !(flip1==flip2);
      sg = (flip == 1) ? 1 : -1;
      cycle  =  (4+cycle2+sg*cycle1)%4;
      
      (*orientrel) = 4*flip + cycle;
      
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK; 
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateJacobianData
static int 
xf_CreateJacobianData( xf_JacobianData **pJData)
{
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pJData, 1, sizeof(xf_JacobianData)));
  if (ierr != xf_OK) return ierr;
  
  (*pJData)->nq        = 0;
  (*pJData)->nqmax     = 0;
  (*pJData)->dim       = 0;
  (*pJData)->detJ      = NULL;
  (*pJData)->J         = NULL;
  (*pJData)->iJ        = NULL;
  (*pJData)->GPhi      = NULL;
  (*pJData)->GPhi_Basis= xfe_BasisLast;
  (*pJData)->GPhi_Order= -1;
  (*pJData)->T         = NULL;
  (*pJData)->Tsize     = 0;
  (*pJData)->AllocFlag = xfb_None;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyJacobianData
int 
xf_DestroyJacobianData( xf_JacobianData *JData)
{
  if (JData == NULL) return xf_OK;
  
  xf_Release( (void *) JData->detJ);
  xf_Release( (void *) JData->J);
  xf_Release( (void *) JData->iJ);
  xf_Release( (void *) JData->GPhi);
  xf_Release( (void *) JData->T);
  xf_Release( (void *) JData);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemJacobian
int 
xf_ElemJacobian(xf_Mesh *Mesh, int egrp, int elem, int nq, real *xq,
                unsigned int AllocFlag,  enum xfe_Bool PointsChanged,
		xf_JacobianData **pJData)
{
  int ierr, Order, nn, dim;
  int iq, i, j, k, l, *Node;
  int iqdim2, off0, off1;
  enum xfe_BasisType Basis;
  enum xfe_Bool ReSize;
  enum xfe_Bool Want_detJ,  Want_J,  Want_iJ;
  enum xfe_Bool Stale_detJ, Stale_J, Stale_iJ;
  xf_JacobianData *JD;
  real A[9], *GPhi, *G, *pdetJ, *piJ;
  real *T, *V, c, t;
  
  if ((*pJData) == NULL){
    ierr = xf_Error(xf_CreateJacobianData(pJData));
    if (ierr != xf_OK) return ierr;
  }
  
  JD = (*pJData);
  
  Basis = Mesh->ElemGroup[egrp].QBasis;
  Order = Mesh->ElemGroup[egrp].QOrder;
  
  // if Jacobian is constant on elem, only evaluate it once
  if ((Order == 1) &&
      (Basis != xfe_HexLagrange) &&
      (Basis != xfe_QuadLagrange))
    nq = 1;
  
  if (Mesh->ElemGroup[egrp].CutFlag) return xf_Error(xf_NOT_SUPPORTED);
  
  // determine nn and dim
  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Basis2Dim(Basis, &dim));
  if (ierr != xf_OK) return ierr;
  
  
  // alloc JD->GPhi if necessary
  if ((JD->GPhi == NULL) || (PointsChanged) || 
      (JD->GPhi_Basis != Basis) || (JD->GPhi_Order != Order)){
    ierr = xf_Error(xf_ReAlloc((void **) &JD->GPhi, 2*nn*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_GetGrads(Basis, Order, nq, xq, JD->GPhi+nn*nq*dim));
    if (ierr != xf_OK) return ierr;

    // Store GPhi in a form more convenient for memory access
    for (iq=0; iq<nq; iq++)
      for (k=0, off0=iq*nn*dim, l=0; k<dim; k++)
	for (i=0, off1=nn*nq*dim+nq*nn*k+iq*nn; i<nn; i++, l++)
	  JD->GPhi[off0+l] = JD->GPhi[off1+i];

    JD->GPhi_Basis = Basis;
    JD->GPhi_Order = Order;
  }
  GPhi = JD->GPhi;
  
  // Reallocate Tsize if necessary)
  if (nn*dim > JD->Tsize){
    JD->Tsize = nn*dim;
    ierr = xf_Error(xf_ReAlloc((void **) &JD->T, JD->Tsize, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  // set nq
  JD->nq = nq;
  
  ReSize = xfe_False;
  if ((nq > JD->nqmax) || (dim != JD->dim)){
    // resizing will be needed
    JD->nqmax = nn*nq;
    JD->dim   = dim;
    ReSize = xfe_True;
  }
  
  Want_detJ = (AllocFlag & xfb_detJ);  
  Want_J    = (AllocFlag & xfb_J);
  Want_iJ   = (AllocFlag & xfb_iJ);
  
  Stale_detJ = !(JD->AllocFlag & xfb_detJ);
  Stale_J    = !(JD->AllocFlag & xfb_J   );
  Stale_iJ   = !(JD->AllocFlag & xfb_iJ  );
  
  // obtain detJ, J, and iJ; reallocate if necessary
  if (Want_detJ && (ReSize || Stale_detJ)){
    ierr = xf_Error(xf_ReAlloc((void **) &JD->detJ, JD->nqmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  
  if (Want_J && (ReSize || Stale_J)){
    ierr = xf_Error(xf_ReAlloc((void **) &JD->J, dim*dim*JD->nqmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  
  if (Want_iJ && (ReSize || Stale_iJ)){
    ierr = xf_Error(xf_ReAlloc((void **) &JD->iJ, dim*dim*JD->nqmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  
  // set AllocFlag
  JD->AllocFlag = AllocFlag;
  
  // Copy over elem geometry nodes for quick access
  Node = Mesh->ElemGroup[egrp].Node[elem];
  for (k=0,l=0; k<dim; k++)
    for (i=0; i<nn; i++,l++)
      JD->T[l] = Mesh->Coord[Node[i]][k];

  // loop over points
  for (iq=0; iq<nq; iq++){

    // compute A = Jacobian
    for (j=0; j<dim; j++){
      T = JD->T + j*nn;
      G = GPhi + iq*nn*dim;
      for (k=0; k<dim; k++){
	for (i=0, t=0.; i<nn; i++) t += T[i]*G[i];
	A[j*dim+k] = t;
	G += nn;
      }
    }
    
    iqdim2 = iq*dim*dim;
    if (Want_J   ) for (i=0; i<dim*dim; i++) JD->J[iqdim2+i] = A[i];
    if (Want_detJ) pdetJ = JD->detJ+iq;       else pdetJ = NULL;
    if (Want_iJ )  piJ   = JD->iJ+iqdim2;     else piJ   = NULL;
    ierr = xf_Error(xf_MatDetInv(A, dim, pdetJ, piJ));
    if (ierr != xf_OK) return ierr;
    
    if ((pdetJ != NULL) && (*pdetJ <= 0)){
      xf_printf("Error, nonpositive Jacobian in ElemJacobian (egrp=%d, elem=%d).\n", 
                egrp, elem);
      return xf_Error(xf_NEGATIVE_JACOBIAN);
    }
    
  } // iq
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpolFace_Triangle
static int 
xf_ScaleHangInterpolFace_Triangle(int face0, int pos, int nq, real *xelem)
{
  // pos != 0; see doc under xf_ScaleHangInterpolFace
  int i;
  real cx, cy;
  
  switch(face0){
    case 0:
      cx = ((pos == xfe_SegPosLeft) ? .5 : 0.); 
      cy = ((pos == xfe_SegPosRight) ? .5 : 0.);
      for (i=0; i<nq; i++){
        xelem[2*i+0] = .5*xelem[2*i+0] + cx;
        xelem[2*i+1] = .5*xelem[2*i+1] + cy;
      }
      break;
      case 1:
      cx = ((pos == xfe_SegPosRight) ? .5 : 0.);
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0] + cx;
      break;
      case 2:
      cy = ((pos == xfe_SegPosLeft) ? .5 : 0.);
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1] + cy;
      break;
      default:
      return xf_Error(xf_OUT_OF_BOUNDS);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpolFace_Quadrilateral
static int 
xf_ScaleHangInterpolFace_Quadrilateral(int face0, int pos, int nq, real *xelem)
{
  // pos != 0; see doc under xf_ScaleHangInterpolFace
  int i;
  real cx, cy;
  
  switch(face0){
    case 0:
      cx = ((pos == xfe_SegPosRight) ? .5 : 0.);
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0] + cx;
      break;
      case 1:
      cy = ((pos == xfe_SegPosRight) ? .5 : 0.);
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1] + cy;
      break;
      case 2:
      cx = ((pos == xfe_SegPosLeft) ? .5 : 0.);
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0] + cx;
      break;
      case 3:
      cy = ((pos == xfe_SegPosLeft) ? .5 : 0.);
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1] + cy;
      break;
      default:
      return xf_Error(xf_OUT_OF_BOUNDS);
      break;
  }
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpolFace_Hexahedron
static int 
xf_ScaleHangInterpolFace_Hexahedron(int face0, int pos, int nq, real *xelem)
{
  // pos != 0; see doc under xf_ScaleHangInterpolFace
  int i, j;
  int iref[2]; // two axes involved in scaling
  real mref[2]; // multiplicative scalings
  real aref[2]; // additive constants
  real x0;
  
  // how are ref values on the quad scaled?
  switch(pos){
    case  xfe_QuadPosNone:
      return xf_Error(xf_INPUT_ERROR);
      break;
    case  xfe_QuadPosSW:
      mref[0] = mref[1] = 0.5;
      aref[0] = 0.; aref[1] = 0.;
      break;
    case  xfe_QuadPosSE:
      mref[0] = mref[1] = 0.5;
      aref[0] = 0.5; aref[1] = 0.;
      break;
    case  xfe_QuadPosNW:
      mref[0] = mref[1] = 0.5;
      aref[0] = 0.; aref[1] = 0.5;
      break;
    case  xfe_QuadPosNE:
      mref[0] = mref[1] = 0.5;
      aref[0] = 0.5; aref[1] = 0.5;
      break;
    case  xfe_QuadPosLeft:
      mref[0] = 0.5; mref[1] = 1.0;
      aref[0] = 0.; aref[1] = 0.;
      break;
    case  xfe_QuadPosRight:
      mref[0] = 0.5; mref[1] = 1.0;
      aref[0] = 0.5; aref[1] = 0.;
      break;
    case  xfe_QuadPosBottom:
      mref[0] = 1.0; mref[1] = 0.5;
      aref[0] = 0.; aref[1] = 0.;
      break;
    case  xfe_QuadPosTop:
      mref[0] = 1.0; mref[1] = 0.5;
      aref[0] = 0.; aref[1] = 0.5;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  // determine axes involved based on face0
  switch(face0){
    case 0: iref[0] = 1; iref[1] = 0; break;
    case 1: iref[0] = 0; iref[1] = 2; break;
    case 2: iref[0] = 1; iref[1] = 2; break;
    case 3: iref[0] = 0; iref[1] = 2; aref[0] = 1.-mref[0]-aref[0]; break;
    case 4: iref[0] = 1; iref[1] = 2; aref[0] = 1.-mref[0]-aref[0]; break;
    case 5: iref[0] = 0; iref[1] = 1; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
  }
  
  // scale the involved coordinates
  for (i=0; i<nq; i++) 
    for (j=0; j<2; j++){
      x0 = xelem[3*i+iref[j]];
      x0 = mref[j]*x0 + aref[j];
      xelem[3*i+iref[j]] = x0;
    } // j
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpolFace
static int 
xf_ScaleHangInterpolFace(enum xfe_ShapeType Shape, int face0, int pos, 
                         int nq, real *xelem)
{
  /*
   PURPOSE:
   
   Scale element ref coords obtained from RefFace2Interpol when called
   on the coarse side of a hanging face.  The input xelem is actually
   in the ref-space of a finer, mirror sub-element.  This function
   scales xelem to convert it to the ref coords of the original coarse
   element.  All input coords must lie on face0.
   
   INPUTS:
   
   Shape : element shape
   face0 : original face number
   pos : position of hanging face on the original face
   nq : number of points at which to calculate
   xelem : unrolled ref elem coords in mirror sub-element
   
   OUTPUTS: 
   
   xelem : unrolled ref elem coords on coarse elem
   
   RETURN:
   
   Error Code
   */
  
  // should only be called for non-conforming elements
  if (pos == 0) return xf_Error(xf_INPUT_ERROR);
  
  switch (Shape){
    case xfe_Triangle:   
      return xf_ScaleHangInterpolFace_Triangle(face0, pos, nq, xelem);
      break;
    case xfe_Tetrahedron:  
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    case xfe_Quadrilateral:
      return xf_ScaleHangInterpolFace_Quadrilateral(face0, pos, nq, xelem);
      break;
    case xfe_Hexahedron:
      return xf_ScaleHangInterpolFace_Hexahedron(face0, pos, nq, xelem);
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpol_Segment
static int 
xf_ScaleHangInterpol_Segment(int pos, int nq, real *xelem)
{
  int i;
  switch(pos){
  case xfe_SegPosNone:
    break;
  case xfe_SegPosLeft:
    for (i=0; i<nq; i++) xelem[i] = .5*xelem[i+0];
    break;
  case xfe_SegPosRight:
    for (i=0; i<nq; i++) xelem[i] = .5*xelem[i+0] + .5;
    break;
  default:
    return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpol_Triangle
static int 
xf_ScaleHangInterpol_Triangle(int pos, int nq, real *xelem)
{
  int i;
  switch(pos){
    case xfe_TriPosNone:
      break;
    case xfe_TriPosLeft:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0];
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1];
      break;
      case xfe_TriPosRight:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0]+.5;
      break;
      case xfe_TriPosCenter:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5-.5*xelem[2*i+0];
      for (i=0; i<nq; i++) xelem[2*i+1] = .5-.5*xelem[2*i+1];
      break;
      case xfe_TriPosTop:
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1]+.5;
      break;
      default:
      return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpol_Quadrilateral
static int 
xf_ScaleHangInterpol_Quadrilateral(int pos, int nq, real *xelem)
{
  int i;
  
  switch(pos){
    case xfe_QuadPosNone:
      break;
    case xfe_QuadPosSW:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0];
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1];
      break;
      case xfe_QuadPosSE:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0]+.5;
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1];
      break;
      case xfe_QuadPosNW:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0];
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1]+.5;
      break;
      case xfe_QuadPosNE:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0]+.5;
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1]+.5;
      break;
      case xfe_QuadPosLeft:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0];
      break;
      case xfe_QuadPosRight:
      for (i=0; i<nq; i++) xelem[2*i+0] = .5*xelem[2*i+0]+.5;
      break;
      case xfe_QuadPosBottom:
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1];
      break;
      case xfe_QuadPosTop:
      for (i=0; i<nq; i++) xelem[2*i+1] = .5*xelem[2*i+1]+.5;
      break;
      default:
      return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpol_Hexahedron
static int 
xf_ScaleHangInterpol_Hexahedron(int pos, int nq, real *xelem)
{
  int i, j;
  real mref[3]; // multiplicative scalings
  real aref[3]; // additive constants
  
  // this relies on the order in xf_MeshHangStruct.h
  if (pos == xfe_HexPosNone){
    mref[0] = 1.0; mref[1] = 1.0; mref[2] = 1.0;
    aref[0] = 0. ; aref[1] = 0. ; aref[2] = 0. ;
  }
  else if (pos <= xfe_HexPos111){  // uniform
    i = pos - xfe_HexPos000;
    mref[0] = 0.5; mref[1] = 0.5; mref[2] = 0.5;
    aref[0] = 0.5* (i % 2);
    aref[1] = 0.5* ((i/2) % 2);
    aref[2] = 0.5* (i / 4);
  }
  else if (pos >= xfe_HexPos200){ // SliceYZ
    i = pos - xfe_HexPos200;
    mref[0] = 1.0; mref[1] = 0.5; mref[2] = 0.5;
    aref[0] = 0.0; aref[1] = 0.5* (i%2); aref[2] = 0.5* (i/2);
  }
  else if (pos >= xfe_HexPos020){ // SliceXZ
    i = pos - xfe_HexPos020;
    mref[0] = 0.5; mref[1] = 1.0; mref[2] = 0.5;
    aref[0] = 0.5* (i%2); aref[1] = 0.0; aref[2] = 0.5* (i/2);
  }
  else if (pos >= xfe_HexPos002){ // SliceXY
    i = pos - xfe_HexPos002;
    mref[0] = 0.5; mref[1] = 0.5; mref[2] = 1.0;
    aref[0] = 0.5* (i%2); aref[1] = 0.5* (i/2); aref[2] = 0.0; 
  }
  else{
    switch(pos){
      case xfe_HexPos022: // SliceX, left
        mref[0] = 0.5; mref[1] = 1.0; mref[2] = 1.0;
        aref[0] = 0. ; aref[1] = 0. ; aref[2] = 0. ;
        break;
      case xfe_HexPos122: // SliceX, right
        mref[0] = 0.5; mref[1] = 1.0; mref[2] = 1.0;
        aref[0] = 0.5; aref[1] = 0. ; aref[2] = 0. ;
        break;
      case xfe_HexPos202: // SliceY, front
        mref[0] = 1.0; mref[1] = 0.5; mref[2] = 1.0;
        aref[0] = 0. ; aref[1] = 0. ; aref[2] = 0. ;
        break;
      case xfe_HexPos212: // SliceY, back
        mref[0] = 1.0; mref[1] = 0.5; mref[2] = 1.0;
        aref[0] = 0. ; aref[1] = 0.5; aref[2] = 0. ;
        break;
      case xfe_HexPos220: // SliceZ, bottom
        mref[0] = 1.0; mref[1] = 1.0; mref[2] = 0.5;
        aref[0] = 0. ; aref[1] = 0. ; aref[2] = 0. ;
        break;
      case xfe_HexPos221: // SliceZ, top
        mref[0] = 1.0; mref[1] = 1.0; mref[2] = 0.5;
        aref[0] = 0. ; aref[1] = 0. ; aref[2] = 0.5;
        break;
      default:
        return xf_Error(xf_INPUT_ERROR);
        break;
    } 
  }
  
  // scale the involved coordinates
  for (i=0; i<nq; i++) 
    for (j=0; j<3; j++)
      xelem[3*i+j] = mref[j]*xelem[3*i+j] + aref[j];
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpol
int 
xf_ScaleHangInterpol(enum xfe_ShapeType Shape, int pos, int nq, real *xelem)
{
  // return immediately if not refined
  if (pos == 0) return xf_OK;
  
  switch (Shape){
    case xfe_Segment:   
      return xf_ScaleHangInterpol_Segment(pos, nq, xelem);
      break;
    case xfe_Triangle:   
      return xf_ScaleHangInterpol_Triangle(pos, nq, xelem);
      break;
    case xfe_Tetrahedron:  
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    case xfe_Quadrilateral:
      return xf_ScaleHangInterpol_Quadrilateral(pos, nq, xelem);
      break;
    case xfe_Hexahedron:
      return xf_ScaleHangInterpol_Hexahedron(pos, nq, xelem);
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScaleHangInterpolInv
int 
xf_ScaleHangInterpolInv(enum xfe_ShapeType Shape, int pos, int nq, real *xelem)
{
  int ierr, i, j, dim;
  real x[3], b[3], A[9], iA[9];
  
  // return immediately if not refined
  if (pos == 0) return xf_OK;
  
  // get dim
  ierr = xf_Error(xf_Shape2Dim(Shape, &dim));
  if (ierr != xf_OK) return ierr;
  
  /*
   determine A an b in:  xH = A*xh + b
   where xH = coord in parent, xh = coord in child
   */
  
  // first see where origin of child maps to -- this is b
  for (j=0; j<dim; j++) x[j] = 0.;
  ierr = xf_Error(xf_ScaleHangInterpol(Shape, pos, 1, x));
  if (ierr != xf_OK) return ierr; 
  for (j=0; j<dim; j++) b[j] = x[j];
  
  // next get columns of A
  for (i=0; i<dim; i++){
    for (j=0; j<dim; j++) x[j] = 0.;
    x[i] = 1.0;
    ierr = xf_Error(xf_ScaleHangInterpol(Shape, pos, 1, x));
    if (ierr != xf_OK) return ierr; 
    for (j=0; j<dim; j++) A[j*dim+i] = x[j] - b[j];
  }
  
  // compute inverse of A -> iA
  ierr = xf_Error(xf_MatDetInv(A, dim, NULL, iA));
  if (ierr != xf_OK) return ierr;
  
  
  // map all requested points: xh = iA*(xH-b)
  for (i=0; i<nq; i++){
    for (j=0; j<dim; j++) x[j] = xelem[dim*i+j] - b[j];
    xf_MxV_Set(iA, x, dim, dim, xelem+dim*i);
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Pos2Hang
int 
xf_Pos2Hang(enum xfe_ShapeType FShape, int nface, int face0, 
            int pos, int *phang)
{ 
  (*phang) = nface * pos + face0;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefFace2Elem_Segment
static int 
xf_RefFace2Elem_Segment(int face, int faceorient, int nq, 
                        real *xface, real *xelem)
{
  
  if (faceorient < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  
  if (face == 0)
    xelem[0] = 0.;
  else if (face == 1)
    xelem[0] = 1.0;
  else return xf_Error(xf_INPUT_ERROR);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_RefFace2Elem_Triangle
static int 
xf_RefFace2Elem_Triangle(int face, int faceorient, int nq, 
                         real *xface, real *xelem)
{
  int i, n0, n1, itemp;
  const real xn[3] = {0., 1., 0.};
  const real yn[3] = {0., 0., 1.};
  real x0, y0, x1, y1;
  real b0, b1, *xref;
  
  if (faceorient < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // use face and faceorient to determine player nodes
  n0 = (face+1)%3;
  n1 = (face+2)%3;
  if (faceorient == 1) swap(n0,n1, itemp);
  
  // coords of player nodes
  x0 = xn[n0]; y0 = yn[n0];
  x1 = xn[n1]; y1 = yn[n1];
  
  // loop over points and set xelem
  for (i=0; i<nq; i++){ 
    // barycentric coords
    b1 = xface[i];  
    b0 = 1.0-b1;  
    
    xref = xelem + 2*i;
    xref[0] = b0*x0 + b1*x1;
    xref[1] = b0*y0 + b1*y1;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefFace2Elem_Quadrilateral
static int 
xf_RefFace2Elem_Quadrilateral(int face, int faceorient, int nq, 
                              real *xface, real *xelem)
{
  int ierr, i, nv[2], itemp;
  const real xn[4] = {0., 1., 0., 1.};
  const real yn[4] = {0., 0., 1., 1.};
  real x0, y0, x1, y1;
  real b0, b1, *xref;
  
  if (faceorient < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // use face and faceorient to determine player nodes
  ierr = xf_Error(xf_Q1NodesOnFace(xfe_QuadLagrange, 1, face, &itemp, nv));
  if (ierr != xf_OK) return ierr;
  if (faceorient == 1) swap(nv[0],nv[1], itemp);
  
  // coords of player nodes
  x0 = xn[nv[0]]; y0 = yn[nv[0]];
  x1 = xn[nv[1]]; y1 = yn[nv[1]];
  
  // loop over points and set xelem
  for (i=0; i<nq; i++){ 
    // barycentric coords
    b1 = xface[i];  
    b0 = 1.0-b1;  
    
    xref = xelem + 2*i;
    xref[0] = b0*x0 + b1*x1;
    xref[1] = b0*y0 + b1*y1;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefFace2Elem_Tetrahedron
static int 
xf_RefFace2Elem_Tetrahedron(int face, int faceorient, int nq, 
                            real *xface, real *xelem)
{
  int i, n0, n1, n2, itemp, cycle, flip;
  int F2N[4][3] = {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}};
  const real xn[4] = {0., 1., 0., 0.};
  const real yn[4] = {0., 0., 1., 0.};
  const real zn[4] = {0., 0., 0., 1.};
  real x0, y0, z0, x1, y1, z1, x2, y2, z2;
  real b0, b1, b2, *xref;
  
  if (faceorient < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // use face and faceorient to determine player nodes
  n0 = F2N[face][0];
  n1 = F2N[face][1];
  n2 = F2N[face][2];
  
  cycle = faceorient%3;
  flip  = faceorient/3;
  if (cycle == 1) cycle3(n0,n1,n2, itemp);
  if (cycle == 2) cycle3(n2,n1,n0, itemp);
  if (flip) swap(n1,n2, itemp);
  
  // coords of player nodes
  x0 = xn[n0]; y0 = yn[n0]; z0 = zn[n0];
  x1 = xn[n1]; y1 = yn[n1]; z1 = zn[n1];
  x2 = xn[n2]; y2 = yn[n2]; z2 = zn[n2];
  
  // loop over points and set xelem
  for (i=0; i<nq; i++){ 
    // barycentric coords
    b1 = xface[2*i+0];
    b2 = xface[2*i+1];
    b0 = 1.0-b1-b2;
    
    xref = xelem + 3*i;
    xref[0] = b0*x0 + b1*x1 + b2*x2;
    xref[1] = b0*y0 + b1*y1 + b2*y2;
    xref[2] = b0*z0 + b1*z1 + b2*z2;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_RefFace2Elem_Hexahedron
static int 
xf_RefFace2Elem_Hexahedron(int face, int faceorient, int nq, 
                           real *xface, real *xelem)
{
  int ierr, i, j, nv0[4], nv[4], itemp, cycle, flip;
  const real xn[8] = {0., 1., 0., 1., 0., 1., 0., 1.};
  const real yn[8] = {0., 0., 1., 1., 0., 0., 1., 1.};
  const real zn[8] = {0., 0., 0., 0., 1., 1., 1., 1.};
  real x[4], y[4], z[4], X, Y, b[4];
  real *xref;
  
  if (faceorient < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // use face and faceorient to determine player nodes
  ierr = xf_Error(xf_Q1NodesOnFace(xfe_HexLagrange, 1, face, &itemp, nv0));
  if (ierr != xf_OK) return ierr;
  cycle = faceorient%4;
  flip  = faceorient/4;
  for (i=0; i<4; i++) nv[i] = nv0[(i+cycle)%4];
  if (flip) swap(nv[1], nv[3], itemp);
  
  // coords of player nodes
  for (j=0; j<4; j++){
    x[j] = xn[nv[j]]; 
    y[j] = yn[nv[j]];
    z[j] = zn[nv[j]];
  }
  
  // loop over points and set xelem
  for (i=0; i<nq; i++){ 
    // shape functions
    X = xface[2*i+0];
    Y = xface[2*i+1];
    b[0] = (1.-X)*(1.-Y);
    b[1] =     X *(1.-Y);
    b[2] =     X *    Y ;
    b[3] = (1.-X)*    Y ;
    
    xref = xelem + 3*i;
    for (j=0, xref[0]=0.; j<4; j++) xref[0] += b[j]*x[j];
    for (j=0, xref[1]=0.; j<4; j++) xref[1] += b[j]*y[j];
    for (j=0, xref[2]=0.; j<4; j++) xref[2] += b[j]*z[j];
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefFace2Elem
static int 
xf_RefFace2Elem(enum xfe_ShapeType Shape, int face, int faceorient,
                int nq, real *xface, real *xelem){
  
  switch (Shape){
    case xfe_Segment:
      return xf_RefFace2Elem_Segment(face, faceorient, nq, xface, xelem);
      break;
    case xfe_Triangle:   
      return xf_RefFace2Elem_Triangle(face, faceorient, nq, xface, xelem);
      break;
    case xfe_Tetrahedron:  
      return xf_RefFace2Elem_Tetrahedron(face, faceorient, nq, xface, xelem);
      break;
    case xfe_Quadrilateral:
      return xf_RefFace2Elem_Quadrilateral(face, faceorient, nq, xface, xelem);
      break;
    case xfe_Hexahedron:
      return xf_RefFace2Elem_Hexahedron(face, faceorient, nq, xface, xelem);
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefFace2Interpol
int 
xf_RefFace2Interpol(xf_Mesh *Mesh, int egrp, int elem, int face, int faceorient,
                    int nq, real *xface, real *xelem)
{
  int ierr, hang, face0, fpos, epos;
  enum xfe_ShapeType Shape;
  
  // determine Shape of element
  ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  face0 = face; // original face number ... may be different for non-conforming elems
  
  // check if on the coarse side of a hanging face
  ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, &face0, 
                                   &fpos, &epos));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_RefFace2Elem(Shape, face0, faceorient, nq, xface, xelem));
  if (ierr != xf_OK) {
    xf_printf("egrp = %d, elem = %d, face0 = %d, hang = %d\n", egrp, elem, face0, hang);
    return ierr;
  }
  
  
  if (hang != 0){
    // scale to coarse elem interpol coords if a non-conforming elem
    ierr = xf_Error(xf_ScaleHangInterpol(Shape, epos, nq, xelem));
    if (ierr != xf_OK) return ierr; 
  }
  
  if (Mesh->ElemGroup[egrp].CutFlag){
    // convert ref space coords to interpolation coords
    return xf_Error(xf_NOT_SUPPORTED);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Ref2GlobElem
int 
xf_Ref2GlobElem(xf_Mesh *Mesh, int egrp, int elem, xf_BasisData **pPhiData, 
                enum xfe_Bool PointsChanged, int npoint, real *xref, real *xglob)
{
  int ierr, dim, nn, d, ipoint, n;
  int ipd, ipnn, ig, *Node;
  enum xfe_BasisType QBasis;
  int QOrder;
  real *Phi, **Coord, val;
  xf_BasisData *PhiDataLoc = NULL, **pPhi;
  
  if (pPhiData == NULL)
    pPhi = &PhiDataLoc;
  else
    pPhi = pPhiData;
  
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  
  if ((PointsChanged) || ((*pPhi) == NULL) || 
      ((*pPhi)->Basis != QBasis) || ((*pPhi)->Order != QOrder)){
    ierr = xf_Error(xf_EvalBasis(QBasis, QOrder, PointsChanged, 
                                 npoint, xref, xfb_Phi, pPhi));
    if (ierr != xf_OK) return ierr;
  }
  
  dim   = Mesh->Dim;
  Coord = Mesh->Coord;
  
  Phi = (*pPhi)->Phi;
  nn  = (*pPhi)->nn;
  
  if (nn != Mesh->ElemGroup[egrp].nNode) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  Node = Mesh->ElemGroup[egrp].Node[elem];
  
  for (ipoint=0; ipoint<npoint; ipoint++){
    ipd  = ipoint*dim;
    ipnn = ipoint*nn;
    for (d=0; d<dim; d++) xglob[ipd+d] = 0.0;
    for (n=0; n<nn; n++){
      ig  = Node[n];
      val = Phi[ipnn+n];
      for (d=0; d<dim; d++) xglob[ipd+d] += val*Coord[ig][d];
    }
  }
  
  /* Destroy Basis Data */
  if (pPhiData == NULL){
    ierr = xf_Error(xf_DestroyBasisData(PhiDataLoc, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Ref2GlobFaceUsingTable
int 
xf_Ref2GlobFaceUsingTable(xf_Mesh *Mesh, int egrp, int elem, int face,
			  int faceorient, enum xfe_Bool PointsChanged, 
			  xf_BasisData **pBasisData, xf_BasisTable *BasisTable,
			  int npoint, real *xref, real *xglob)
{
  int ierr, hang;
  enum xfe_ShapeType Shape;
  
  // if PointsChanged all stored bases need re-evaluating
  if (PointsChanged){
    ierr = xf_Error(xf_ForceReCalcBasisTable(BasisTable));
    if (ierr != xf_OK) return ierr;
  }

  // determine Shape of element
  ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
  if (ierr != xf_OK) return ierr;

  // Special check for non-conforming elements
  ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, NULL, NULL, NULL));
  if (ierr != xf_OK) return ierr;

  if (hang != 0){
    // if pBasisData points to a table, reset it to NULL
    if (((*pBasisData) != NULL) && ((*pBasisData)->InTable))
      (*pBasisData) = NULL;
    
    // we're on the coarse element of a hanging face; do not use table
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, pBasisData,
				    PointsChanged, npoint, xref, xglob));
    if (ierr != xf_OK) return ierr;
    (*pBasisData)->InTable = xfe_False;

    return xf_OK;
  }

  // Call Ref2GlobElem using table
  ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, 
				  &BasisTable->BasisData[Shape][face][faceorient],
				  PointsChanged, npoint, xref, xglob));
  if (ierr != xf_OK) return ierr;
  
  // if (*pBasisData) points to non-table data, destroy it
  if ((*pBasisData) != NULL){
    if ((*pBasisData)->InTable == xfe_False){
      ierr = xf_Error(xf_DestroyBasisData((*pBasisData), xfe_True));
      if (ierr != xf_OK) return ierr;
    }
  }

  // point (*pBasisData) to the table
  (*pBasisData) = BasisTable->BasisData[Shape][face][faceorient];
  (*pBasisData)->InTable = xfe_True;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Glob2RefElem
int 
xf_Glob2RefElem(xf_Mesh *Mesh, int egrp, int elem, real *xglob, 
		real intol, enum xfe_Bool Verbose,
                real *xref, enum xfe_Bool *pconverged)
{
  int ierr, dim, d;
  enum xfe_Bool converged, NonLinear;
  enum xfe_BasisType QBasis;
  enum xfe_ShapeType Shape;
  int QOrder, iNewton, iTry;
  int maxTry    = 5;     // maximum number of tries of Newton
  int maxNewton = 50;    // maximum number of Newton iterations
  real mintol   = 1e-12; // minimum absolute residual tolerance 
  real reltol   = 1e-12; // relative residual tolerance;
  real maxdxref = 0.5;   // maximum update in xref
  real omega, rnorm, tol;
  real R[3], xr[3], x0[3], x1[3], dxref[3];
  xf_JacobianData *JData = NULL;
  
  // set tolerance
  if (tol > 0.) mintol = intol;

  // basis and order of element
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  
  // is basis/order combination nonlinear?
  
  //NonLinear = ((QOrder > 1) || (xf_Q1BasisNotLinear(QBasis)));

  //Yu
  //temporarily use linear option for cube-type elements
  NonLinear = (QOrder > 1);
  
  // determine element Shape 
  ierr = xf_Error(xf_Basis2Shape(QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  // dimension
  dim = Mesh->Dim;
  
  // initialize xref to centroid
  ierr = xf_Error(xf_ShapeCentroid(Shape, xref));
  if (ierr != xf_OK) return ierr;
  
  
  if (!NonLinear){
    // Linear mapping, solve for ref coords in one shot
    
    // element Jacobian
    ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, 1, xref, xfb_iJ, xfe_True, &JData));
    if (ierr != xf_OK) return ierr;
    
    // where does xref=zero map to?  Store in R
    for (d=0; d<dim; d++) xref[d] = 0.;
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True,
                                    1, xref, R));
    if (ierr != xf_OK) return ierr;
    
    // R = xglob - x(0)
    for (d=0; d<dim; d++) R[d] = xglob[d]-R[d];
    
    // xref = iJ*R
    xf_MxV_Set(JData->iJ, R, dim, dim, xref);
    
    // status is converged
    converged = xfe_True;
    
  }
  else{  
    // Use Newton
    
    // -- first determine appropriate tolerance --
    // where does xref=0 map to?
    for (d=0; d<dim; d++) xr[d] = 0.;
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True,
                                    1, xr, x0));
    if (ierr != xf_OK) return ierr;
    // where does xref=1 map to?
    for (d=0; d<dim; d++) xr[d] = 1.;
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True,
                                    1, xr, x1));
    if (ierr != xf_OK) return ierr;
    // difference
    for (d=0, rnorm=0.; d<dim; d++) rnorm += (x1[d]-x0[d])*(x1[d]-x0[d]);
    rnorm = sqrt(rnorm);
    // tol based on rnorm
    tol = max(mintol, reltol*rnorm);
    
    // begin Newton loop
    converged = xfe_False;
    for (iTry=0; (iTry<maxTry) && (!converged); iTry++){
      if (iTry > 0)
        for (d=0; d<dim; d++) xref[d] = ((real) rand())/((real) RAND_MAX);
      for (iNewton=0; iNewton<maxNewton; iNewton++){
        
        // calculate residual: R = x(xref) - xglob
        ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True,
                                        1, xref, R));
        if (ierr != xf_OK) return ierr;
        for (d=0; d<dim; d++) R[d] -= xglob[d];
        
        // check tolerance on residual
        for (d=0, rnorm=0.; d<dim; d++) rnorm += R[d]*R[d];
        rnorm = sqrt(rnorm);
        
        //xf_printf("iNewton = %d, rnorm = %.10E\n", iNewton, rnorm);
        
        if (rnorm < tol){ // converged!
          converged = xfe_True;
          break;
        }
        
        // calculate residual linearization inverse, iR_xref
        ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, 1, xref, xfb_iJ, xfe_True, &JData));
        if (ierr != xf_OK) return ierr;
        
        // determine state update: dxref = -iR_xref*R
        xf_MxV_Neg(JData->iJ, R, dim, dim, dxref);
        
        // limit state update
        for (d=0,omega=1.0; d<dim; d++)
          if (fabs(dxref[d]) > maxdxref)
            omega = min(omega, maxdxref/fabs(dxref[d]));
        
        // apply state update
        for (d=0; d<dim; d++) xref[d] += omega*dxref[d];
        
      } // iNewton
    } // iTry
  }
 
  // print out warning if not converged
  if ((!converged) && (Verbose))
    xf_printf("Warning: Newton not converged in Glob2RefElem (rnorm = %.6E).\n", rnorm);
  
  // return convergence status if requested
  if (pconverged != NULL) (*pconverged) = converged;
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_CreateNormalData
static int 
xf_CreateNormalData( xf_NormalData **pNData)
{
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pNData, 1, sizeof(xf_NormalData)));
  if (ierr != xf_OK) return ierr;
  
  (*pNData)->nq     = 0;
  (*pNData)->ndqmax = 0;
  (*pNData)->dim    = 0;
  (*pNData)->n      = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyNormalData
int 
xf_DestroyNormalData( xf_NormalData *NData)
{
  if (NData == NULL) return xf_OK;
  
  xf_Release( (void *) NData->n);
  xf_Release( (void *) NData);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemNormal_Segment
static int 
xf_ElemNormal_Segment(xf_Mesh *Mesh, int egrp, int elem, int face, 
                      int faceorient, int nq, real *xface, 
                      real *nvec)
{
  if (nq != 1) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  if (face == 0)      // left face
    nvec[0] = -1.0;
  else if (face == 1) // right face
    nvec[0] =  1.0;
  else
    return xf_Error(xf_INPUT_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemNormal_TriangleQuad
static int 
xf_ElemNormal_TriangleQuad(xf_Mesh *Mesh, int egrp, int elem, int face, 
                           int faceorient, int nq, real *xface, 
                           real *nvec)
{
  int ierr, iq, j, d, itemp, nnodeq1;
  int QOrder, nn, fvec0[2], *fvec, *Node;
  enum xfe_BasisType QBasis;
  real *x0, *x1, *xn;
  real s, x_s[2];
  real *GPhi, **coord;
  
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  
  // make sure Lagrange
  if ((QBasis != xfe_TriLagrange) && (QBasis != xfe_QuadLagrange))
    return xf_Error(xf_NOT_SUPPORTED);
  
  // set Node pointer to elem's nodes
  Node = Mesh->ElemGroup[egrp].Node[elem];
  
  // if Q1, quick calc and return
  if (QOrder == 1){
    
    if (nq != 1) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    ierr = xf_Error(xf_Q1NodesOnFace(QBasis, 1, face, &itemp, fvec0));
    if (ierr != xf_OK) return ierr;
    
    x0 = Mesh->Coord[Node[fvec0[0]]];
    x1 = Mesh->Coord[Node[fvec0[1]]];
    
    nvec[0] =  (x1[1]-x0[1]);
    nvec[1] = -(x1[0]-x0[0]);
    
    return xf_OK;
  }
  
  // number of nodes per face
  nn = (QOrder+1);
  
  // alloc vars
  ierr = xf_Error(xf_Alloc((void **) &fvec, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &GPhi, nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &coord, nn, sizeof(real *)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) &xn, nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // set xn
  ierr = xf_Error(xf_LagrangeNodes(xfe_SegLagrange, QOrder, NULL, xn, NULL));
  if (ierr != xf_OK) return ierr;
  
  // determine loc nodes on the face as seen by ref element
  ierr = xf_Error(xf_NodesOnFace(QBasis, QOrder, face, &itemp, fvec));
  if (ierr != xf_OK) return ierr;
  
  // set coord
  for (j=0; j<nn; j++) coord[j] = Mesh->Coord[Node[fvec[j]]];
  
  // loop over points
  for (iq=0; iq<nq; iq++){
    if (faceorient == 0)
      s = xface[iq];
    else
      s = 1-xface[iq];

    // compute gradients
    ierr = xf_Error(xf_BasisLagrange1D(s, xn, nn, NULL, GPhi, NULL));
    if (ierr != xf_OK) return ierr;
    
    // calculate x_s
    for (d=0; d<2; d++)
      for (j=0, x_s[d] = 0.0; j<nn; j++) x_s[d] += GPhi[j]*coord[j][d];
    
    // take cross product, store in nvec
    nvec[2*iq+0] =  x_s[1];
    nvec[2*iq+1] = -x_s[0];
    
  } // iq
  
  // release memory
  xf_Release( (void *) fvec);
  xf_Release( (void *) GPhi);
  xf_Release( (void *) coord);
  xf_Release( (void *) xn);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemNormal_Tetrahedron
static int 
xf_ElemNormal_Tetrahedron(xf_Mesh *Mesh, int egrp, int elem, int face, 
                          int faceorient, int nq, real *xface, 
                          real *nvec)
{
  int ierr, iq, j, d, cycle, flip, itemp;
  int QOrder, nn, fvec0[3], *fvec, *Node;
  enum xfe_BasisType QBasis;
  real *x0, *x1, *x2, *x3;
  real b0, b1, b2, rtemp;
  real x_xi[3], x_eta[3], xr[2];
  real *GPhi, **coord;
  
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  
  // make sure Lagrange
  if (QBasis != xfe_TetLagrange) return xf_Error(xf_NOT_SUPPORTED);
  
  // set Node pointer to elem's nodes
  Node = Mesh->ElemGroup[egrp].Node[elem];
  
  // if Q1, quick calc and return
  if (QOrder == 1){
    
    if (nq != 1) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    // determine loc nodes on the face as seen by ref element
    ierr = xf_Error(xf_Q1NodesOnFace(QBasis, 1, face, &itemp, fvec0));
    if (ierr != xf_OK) return ierr;
    
    x0 = Mesh->Coord[Node[fvec0[0]]];
    x1 = Mesh->Coord[Node[fvec0[1]]];
    x2 = Mesh->Coord[Node[fvec0[2]]];
    
    x_xi[0] = x1[0]-x0[0];  x_eta[0] = x2[0]-x0[0];
    x_xi[1] = x1[1]-x0[1];  x_eta[1] = x2[1]-x0[1];
    x_xi[2] = x1[2]-x0[2];  x_eta[2] = x2[2]-x0[2];
    
    xf_CrossProduct(x_xi, x_eta, nvec);
    
    return xf_OK;
  }
  
  // number of nodes per face
  nn = (QOrder+1)*(QOrder+2)/2;
  
  // alloc vars
  ierr = xf_Error(xf_Alloc((void **) &fvec, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &GPhi, 2*nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &coord, nn, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  
  
  // determine loc nodes on the face as seen by ref element
  ierr = xf_Error(xf_NodesOnFace(QBasis, QOrder, face, &itemp, fvec));
  if (ierr != xf_OK) return ierr;
  
  // set coord
  for (j=0; j<nn; j++) coord[j] = Mesh->Coord[Node[fvec[j]]];
  
  cycle = faceorient%3;
  flip  = faceorient/3;
  
  // loop over points
  for (iq=0; iq<nq; iq++){
    // barycentric coords
    b1 = xface[2*iq+0];
    b2 = xface[2*iq+1];
    b0 = 1.0-b1-b2;
    
    // cycle and flip using faceorient
    if (flip) swap(b1,b2, rtemp);
    if (cycle == 1) cycle3(b2,b1,b0, rtemp);
    if (cycle == 2) cycle3(b0,b1,b2, rtemp);
    
    // convert barycentric to xr
    xr[0] = b1;
    xr[1] = b2;
    
    // compute gradients
    if (xf_Grad_TriLagrange(QOrder, xr, GPhi, nn) != 0) 
      return xf_Error(xf_BASIS_FCN_ERROR);
    
    // calculate x_xi, x_eta
    for (d=0; d<3; d++){
      for (j=0, x_xi[d] =0.; j<nn; j++) x_xi[d]  += GPhi[j+ 0]*coord[j][d];
      for (j=0, x_eta[d]=0.; j<nn; j++) x_eta[d] += GPhi[j+nn]*coord[j][d];
    }
    
    // take cross product, store in nvec
    xf_CrossProduct(x_xi, x_eta, nvec+3*iq);
    
  } // iq
  
  // release memory
  xf_Release( (void *) fvec);
  xf_Release( (void *) GPhi);
  xf_Release( (void *) coord);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemNormal_Hexahedron
static int 
xf_ElemNormal_Hexahedron(xf_Mesh *Mesh, int egrp, int elem, int face, 
                         int faceorient, int nq, real *xface, 
                         real *nvec)
{
  int ierr, iq, j, d, cycle, flip, itemp;
  int QOrder, nn, *fvec, *Node;
  enum xfe_BasisType QBasis;
  real b0, b1, b2, rtemp;
  real x_xi[3], x_eta[3], xr[2];
  real *GPhi, **coord;
  
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  
  // make sure Lagrange
  if (QBasis != xfe_HexLagrange) return xf_Error(xf_NOT_SUPPORTED);
  
  // set Node pointer to elem's nodes
  Node = Mesh->ElemGroup[egrp].Node[elem];
  
  // No quick calc for Q1 as normal may still depend on xface
  
  // number of nodes per face
  nn = (QOrder+1)*(QOrder+1);
  
  // alloc vars
  ierr = xf_Error(xf_Alloc((void **) &fvec, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &GPhi, 2*nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &coord, nn, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  
  
  // determine loc nodes on the face as seen by ref element
  ierr = xf_Error(xf_NodesOnFace(QBasis, QOrder, face, &itemp, fvec));
  if (ierr != xf_OK) return ierr;
  
  // set coord
  for (j=0; j<nn; j++) coord[j] = Mesh->Coord[Node[fvec[j]]];
  
  cycle = faceorient%4;
  flip  = faceorient/4;
  
  // loop over points
  for (iq=0; iq<nq; iq++){
    
    // face ref coords
    xr[0] = xface[2*iq+0];
    xr[1] = xface[2*iq+1];
    
    // cycle and flip using faceorient
    if (flip) swap(xr[0],xr[1], rtemp);
    if ((cycle%2) == 1) swap(xr[0],xr[1], rtemp);
    if ((cycle == 1) || (cycle == 2)) xr[0] = 1.-xr[0];
    if ((cycle == 2) || (cycle == 3)) xr[1] = 1.-xr[1];
    
    // compute gradients
    if (xf_Grad_TensorLagrange(2, QOrder, xr, GPhi, nn) != 0) 
      return xf_Error(xf_BASIS_FCN_ERROR);
    
    // calculate x_xi, x_eta
    for (d=0; d<3; d++){
      for (j=0, x_xi[d] =0.; j<nn; j++) x_xi[d]  += GPhi[j+ 0]*coord[j][d];
      for (j=0, x_eta[d]=0.; j<nn; j++) x_eta[d] += GPhi[j+nn]*coord[j][d];
    }
    
    // take cross product, store in nvec
    xf_CrossProduct(x_xi, x_eta, nvec+3*iq);
    
  } // iq
  
  // release memory
  xf_Release( (void *) fvec);
  xf_Release( (void *) GPhi);
  xf_Release( (void *) coord);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemNormal
int 
xf_ElemNormal(xf_Mesh *Mesh, int egrp, int elem, int face, 
              int faceorient, int nq, real *xface, 
              xf_NormalData **pNData)
{
  int ierr, dim;
  enum xfe_BasisType QBasis;
  enum xfe_ShapeType Shape;
  
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  
  // determine dim
  ierr = xf_Error(xf_Basis2Dim(QBasis, &dim));
  if (ierr != xf_OK) return ierr;
  
  // determine shape
  ierr = xf_Error(xf_Basis2Shape(QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  // check if Q1 (note a hexahedron's normal can vary even if Q1)
  if ((Mesh->ElemGroup[egrp].QOrder == 1) &&
      (QBasis != xfe_HexLagrange)) nq = 1;
  
  // check/(re)allocate (*pNData)
  if ((*pNData) == NULL){
    ierr = xf_Error(xf_CreateNormalData(pNData));
    if (ierr != xf_OK) return ierr;
  }
  (*pNData)->nq  = nq;
  (*pNData)->dim = dim;
  
  if (nq*dim > (*pNData)->ndqmax){
    (*pNData)->ndqmax = nq*dim;
    ierr = xf_Error(xf_ReAlloc((void **) &(*pNData)->n, nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  
  switch(Shape){
    case xfe_Segment:
      return xf_Error(xf_ElemNormal_Segment(Mesh, egrp, elem, face, faceorient, 
                                            nq, xface, (*pNData)->n));
      break;
    case xfe_Triangle:
      return xf_Error(xf_ElemNormal_TriangleQuad(Mesh, egrp, elem, face, faceorient, 
                                                 nq, xface, (*pNData)->n));
      break;
    case xfe_Tetrahedron:
      return xf_Error(xf_ElemNormal_Tetrahedron(Mesh, egrp, elem, face, faceorient, 
                                                nq, xface, (*pNData)->n));
      break;
    case xfe_Quadrilateral:
      return xf_Error(xf_ElemNormal_TriangleQuad(Mesh, egrp, elem, face, faceorient, 
                                                 nq, xface, (*pNData)->n));
      break;
    case xfe_Hexahedron:
      return xf_Error(xf_ElemNormal_Hexahedron(Mesh, egrp, elem, face, faceorient, 
                                               nq, xface, (*pNData)->n));
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK; 
}



/******************************************************************/
//   FUNCTION Definition: xf_IFaceNormal
int 
xf_IFaceNormal(xf_Mesh *Mesh, xf_IFace IFace, int nq, 
               real *xq, xf_NormalData **pNData)
{
  int ierr, i;
  int QOrderL, QOrderR;
  
  QOrderL = Mesh->ElemGroup[IFace.ElemGroupL].QOrder;
  QOrderR = Mesh->ElemGroup[IFace.ElemGroupR].QOrder;
  
  if (((QOrderL <= QOrderR) && (IFace.HangNumber == 0)) || (IFace.HangNumber > 0)){
    // note, hang>0 means R is a coarse hang elem -- easier to eval normal on L
    ierr = xf_Error(xf_ElemNormal(Mesh, IFace.ElemGroupL, IFace.ElemL, 
                                  IFace.FaceL, IFace.OrientL, nq, xq, pNData));
    if (ierr != xf_OK){
      xf_printf("L=(%d,%d,%d), R =(%d,%d,%d), hang = %d\n", IFace.ElemGroupL, IFace.ElemL, IFace.FaceL,
                IFace.ElemGroupR, IFace.ElemR, IFace.FaceR, IFace.HangNumber);
      return ierr;
    }
  }
  else{
    // conversely, hang<0 means L is a coarse hang elem -- eval normal on R
    ierr = xf_Error(xf_ElemNormal(Mesh, IFace.ElemGroupR, IFace.ElemR, 
                                  IFace.FaceR, IFace.OrientR, nq, xq, pNData));
    if (ierr != xf_OK) return ierr;
    for (i = 0; i < (*pNData)->nq*(*pNData)->dim; i++)
      (*pNData)->n[i] = -(*pNData)->n[i];
  }
  return xf_OK; 
}


/******************************************************************/
//   FUNCTION Definition: xf_BFaceNormal
int 
xf_BFaceNormal(xf_Mesh *Mesh, xf_BFace BFace, int nq, 
               real *xq, xf_NormalData **pNData, real *nvec)
{
  int ierr, i, dim, d, iq;
  xf_CutFaceData *CutFaceData;
  
  /* for embedded faces, call geometry kernel */
  if ((CutFaceData = BFace.CutFaceData) != NULL){
    //if (CutFaceData->OrigFace < 0) // signifies embedded
    return xf_Error(xf_NOT_SUPPORTED);
  }
  
  ierr = xf_Error(xf_ElemNormal(Mesh, BFace.ElemGroup, BFace.Elem, BFace.Face, 
                                BFace.Orient, nq, xq, pNData));
  if (ierr != xf_OK) return ierr;
  
  // provide normals at all quad points if requested
  if (nvec != NULL){
    dim = Mesh->Dim;
    for (d=0; d<dim; d++)
      for (iq=0;iq<nq; iq++) 
        nvec[iq*dim+d] = (*pNData)->n[iq*dim*((*pNData)->nq!=1)+d];
  }
  
  return xf_OK; 
}


/******************************************************************/
//   FUNCTION Definition: xf_LinearElemJacobian
int 
xf_LinearElemJacobian(enum xfe_BasisType QBasis, int QOrder, int *nvec, 
                      real **Coord, real *detJ){
  
  int ierr, dim, d0, d1;
  real A[9];
  int iQuadHex[3] = {1,2,4};
  
  if (QOrder != 1) return xf_Error(xf_NOT_SUPPORTED);
  
  // determine dim
  ierr = xf_Error(xf_Basis2Dim(QBasis, &dim));
  if (ierr != xf_OK) return ierr;
  
  switch (QBasis){
  case xfe_SegLagrange:
    (*detJ) = Coord[nvec[1]][0] - Coord[nvec[0]][0];
    break;
  case xfe_TriLagrange:
  case xfe_TetLagrange:
    for (d0=0; d0<dim; d0++)
      for (d1=0; d1<dim; d1++)
	A[d0*dim+d1] = Coord[nvec[d1+1]][d0] - Coord[nvec[0]][d0];
    ierr = xf_Error(xf_MatDetInv(A, dim, detJ, NULL));
    break;
  case xfe_QuadLagrange:
  case xfe_HexLagrange:
    for (d0=0; d0<dim; d0++)
      for (d1=0; d1<dim; d1++)
	A[d0*dim+d1] = Coord[nvec[iQuadHex[d1]]][d0] - Coord[nvec[0]][d0];
    ierr = xf_Error(xf_MatDetInv(A, dim, detJ, NULL));
    break;
  default:
    return xf_NOT_SUPPORTED;
    break;
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindElemGeom
int 
xf_FindElemGeom(xf_All *All, xf_Vector **pEG)
{
  int ierr, iq, d, dim;
  int egrp, elem, face, nqElem, nqFace;
  int QuadOrderElem, QuadOrderFace;
  enum xfe_Bool Found, HaloFlag, QuadChanged;
  xf_Face Face;
  xf_IFace IFace;
  xf_BFace BFace;
  real ElemVol, SurfArea, nmag;
  real *xqElem, *xqFace, *nv;
  xf_QuadData *QuadDataElem, *QuadDataFace;
  xf_JacobianData *JData;
  xf_NormalData *NData;
  xf_Vector *EG;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  
  ierr = xf_Error(xf_FindVector(All, "ElemGeom", xfe_LinkageGlobElem, xfe_ElemGeomLast, 
                                xfe_ElemGeomName, 0, 0, NULL, NULL, NULL, NULL, NULL, 
				xfe_SizeReal, xfe_False, xfe_True, NULL, pEG, &Found));
  if (ierr != xf_OK) return ierr;
  
  // return immediately if vector exists (no need to recreate)
  if (Found) return xf_OK;
  
  EG = (*pEG);
  
  // loop over elements and fill in EG
  QuadDataElem  = NULL;
  QuadDataFace  = NULL;
  JData         = NULL;
  NData         = NULL;
  for (egrp=0; egrp<EG->nArray; egrp++){
    
    HaloFlag = (egrp >= Mesh->nElemGroup);
    
    if (HaloFlag)
      xf_printf("Warning; not computing surface areas for elems in halo group = %d.\n", egrp);
    
    QuadOrderElem = Mesh->Dim*(Mesh->ElemGroup[egrp].QOrder-1);
    QuadOrderFace = (Mesh->Dim-1)*(Mesh->ElemGroup[egrp].QOrder-1);
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // quad points on elem
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrderElem, 
                                  &QuadDataElem, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nqElem = QuadDataElem->nquad;
      xqElem = QuadDataElem->xquad;
      
      /* Compute geometry Jacobian; if not constant, compute at quad
       points.  Note if jacobian is constant, only one Jacobian will
       be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nqElem, xqElem, xfb_detJ, 
				      QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;
      
      // compute element volume
      for (iq=0, ElemVol=0.; iq<nqElem; iq++) 
        ElemVol += QuadDataElem->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      EG->GenArray[egrp].rValue[elem][xfe_EGVolume] = ElemVol;
      
      // compute surface area if not halo element
      SurfArea = 0.;
      if (!HaloFlag){
        for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
          Face = Mesh->ElemGroup[egrp].Face[elem][face];
          if (Face.Group == xf_NULLFACE)
            continue;
          // quad points on face
          ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face, QuadOrderFace, 
                                      &QuadDataFace, &QuadChanged));
          if (ierr != xf_OK) return ierr;
          
          nqFace = QuadDataFace->nquad;
          xqFace = QuadDataFace->xquad;
          
          if (Face.Group == xf_INTERIORFACE){ // Interior Face
            IFace = Mesh->IFace[Face.Number];
            ierr = xf_Error(xf_IFaceNormal(Mesh, IFace, nqFace, xqFace, &NData));
            if (ierr != xf_OK) return ierr;
          }
          else if (Face.Group >= 0) { // Boundary Face
            BFace = Mesh->BFaceGroup[Face.Group].BFace[Face.Number];
            ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, nqFace, xqFace, &NData, NULL));
            if (ierr != xf_OK) return ierr;
          }
          else{
            return xf_Error(xf_NOT_SUPPORTED);
          }
          
          // add to surface area
          for (iq=0; iq<nqFace; iq++){
            nv = NData->n + iq*dim*(NData->nq!=1);
            for (d=0, nmag=0.; d<dim; d++) nmag += nv[d]*nv[d];
            SurfArea += QuadDataFace->wquad[iq]*sqrt(nmag);
          }
          
        } // face
      }
      EG->GenArray[egrp].rValue[elem][xfe_EGSurfArea] = SurfArea;
      
    } // elem
    
  } // egrp
  
  /* Destroy Normal Data */
  ierr = xf_Error(xf_DestroyNormalData(NData));
  if (ierr != xf_OK) return ierr;
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataElem));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataFace));
  if (ierr != xf_OK) return ierr;
  
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_FindElemHMetric
int 
xf_FindElemHMetric(xf_All *All, enum xfe_Bool MakeContinuous, xf_Vector **pEM)
{
  int ierr, i, k, l, dim, dim2;
  int myRank, nProc;
  int egrp, negrp, elem, negrphalo;
  int nn, iglob;
  int nvec[xf_MAXQ1NODE];
  int *OrderVec, *NodeC = NULL, *Node;
  enum xfe_Bool Found;
  enum xfe_Bool IsotropicFlag;
  enum xfe_Bool UseHydraulicDiam;
  enum xfe_BasisType *BasisVec;
  real maxE;
  real M[9], H[9], E[3];
  real xref0[xf_MAXQ1NODE*3];
  real Volume, SurfA;
  real **NodeH = NULL, *EH;
  xf_JacobianData *JData;
  xf_Vector *EM;
  xf_Vector *EG;
  xf_Data *D;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  dim2 = dim*dim;

  /* Number of element groups */
  negrp = Mesh->nElemGroup;

  /* Halo element groups are present in parallel runs */  
  negrphalo = ((Mesh->ParallelInfo == NULL) ? Mesh->nElemGroup : 2*Mesh->nElemGroup);
    
  // allocate Q1 Basis and Order vectors
  ierr = xf_Error(xf_Alloc( (void **) &BasisVec, negrphalo, sizeof(enum xfe_BasisType)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &OrderVec, negrphalo, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // request Q1 Lagrange basis if want continuous metric, else Q0
  for (egrp=0; egrp<negrphalo; egrp++){
    ierr = xf_Error(xf_Basis2UniformLagrange(Mesh->ElemGroup[egrp%negrp].QBasis, BasisVec+egrp));
    if (ierr != xf_OK) return ierr;
    if (BasisVec[egrp] != Mesh->ElemGroup[egrp].QBasis) return xf_Error(xf_NOT_SUPPORTED);
    OrderVec[egrp] = (int) MakeContinuous;
  }
  
  
  // Locate/create an element H metric vector (make parallel so can exchange halos)
  ierr = xf_Error(xf_FindVector(All, "ElemHMetric", xfe_LinkageGlobElem, dim2, 
                                NULL, 0, 0, BasisVec, OrderVec, NULL, NULL, NULL, xfe_SizeReal, 
                                xfe_True, xfe_True, &D, pEM, &Found));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_False;  // can set this to True for debugging .. but dangerous if left on
  
  
  xf_Release( (void  *) BasisVec);
  xf_Release( (void  *) OrderVec);
  
  // return immediately if vector exists (no need to recreate)
  if (Found) return xf_OK;
  
  EM = (*pEM);
  
  // Determine if we're calculating an isotropic metric
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "IsotropicHMetric", &IsotropicFlag));
  if (ierr != xf_OK) return ierr;
  
  // Determine if we're using the hydraulic diameter
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "IsotropicHMetricIsHD", &UseHydraulicDiam));
  if (ierr != xf_OK) return ierr;
  
  if (MakeContinuous){
    // create a nodal vector for averaging the elemental h metrics
    ierr = xf_Error(xf_Alloc2( (void ***) &NodeH, Mesh->nNode, dim2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (k=0; k<Mesh->nNode*dim2; k++) (*NodeH)[k] = 0.;
    
    // create a nodal vector for counting node hits
    ierr = xf_Error(xf_Alloc( (void **) &NodeC, Mesh->nNode, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (k=0; k<Mesh->nNode; k++) NodeC[k] = 0;
  }
  
  // need element geometry (volume, surfarea) if using the hydraulic diameter
  if (UseHydraulicDiam){
    ierr = xf_Error(xf_FindElemGeom(All, &EG));
    if (ierr != xf_OK) return ierr;
  }
  
  JData = NULL;
  
  // loop over self elements and add to NodeH
  
  for (egrp=0; egrp<negrp; egrp++){
    // pull off Q1 nodes
    ierr = xf_Error(xf_Q1Nodes(Mesh->ElemGroup[egrp].QBasis, Mesh->ElemGroup[egrp].QOrder, &nn, nvec));
    if (ierr != xf_OK) return ierr;

    if (EM->Order[egrp] == 0) nn = 1;
    
    // get ref space Q1 Lagrange nodes
    ierr = xf_Error(xf_LagrangeNodes(Mesh->ElemGroup[egrp].QBasis, EM->Order[egrp], 
				     NULL, xref0, NULL));
    if (ierr != xf_OK) return ierr;

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      Node = Mesh->ElemGroup[egrp].Node[elem];
      
      // calculate elem Jacobian at Q1 nodes
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nn, xref0, xfb_J, xfe_True, &JData));
      if (ierr != xf_OK) return ierr;

      for (i=0; i<nn; i++){ // loop over Q1 nodes

        // Create M = J * J^T
        k = ((JData->nq == 1) ? 0 : i*dim2);
        xf_MxMT_Set(JData->J+k, JData->J+k, dim, dim, dim, M);
        
        // Calculate H = sqrt(M); first calc eigenvalues/vectors of M
        ierr = xf_Error(xf_EigSym(dim, M, E));
        if (ierr != xf_OK) return ierr;
        
        // set E = sqrt(E)
        for (k=0; k<dim; k++){
          if (E[k] < 0.0) return xf_Error(xf_CODE_LOGIC_ERROR);
          E[k] = sqrt(E[k]);
        }
        
        // Assemble H = V^T*E*V, V = eigvectors, stored in column-major format in M
        for (k=0; k<dim; k++)
          for (l=0; l<dim; l++) M[k*dim+l] *= sqrt(E[k]);
        
        xf_MTxM_Set(M, M, dim, dim, dim, H);
        
        // use H = eye(dim)*max(E) when isotropic
        if (IsotropicFlag){
          if (UseHydraulicDiam){
            // use hydraulic diameter
            Volume = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
            SurfA  = EG->GenArray[egrp].rValue[elem][xfe_EGSurfArea];
            if (SurfA <= 0.) return xf_Error(xf_OUT_OF_BOUNDS);
            maxE = 2.0*Volume/SurfA;
          }
          else{
            // use maximum stretching length
            maxE = E[0];
            for (k=1; k<dim; k++) maxE = max(maxE, E[k]);
          }
          for (k=0; k<dim2; k++) H[k] = 0.;
          for (k=0; k<dim; k++) H[k*dim+k] = maxE;
        }
        
	if (MakeContinuous){
	  // add Jacobian to NodeH
	  iglob = Node[nvec[i]];
	  NodeC[iglob] += 1;
	  for (k=0; k<dim2; k++) NodeH[iglob][k] += H[k];
	}
	else{
	  EH = EM->GenArray[egrp].rValue[elem];
	  for (k=0; k<dim2; k++) EH[k] = H[k];
	}
      }
      
    } // elem
    
  } // egrp

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  if (MakeContinuous){

    // parallel sum NodeC and NodeH
    if (Mesh->ParallelInfo != NULL){
      ierr = xf_Error(xf_ParallelSumNodes(Mesh, (void *) NodeC   ,    1, xfe_SizeInt));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ParallelSumNodes(Mesh, (void *) NodeH[0], dim2, xfe_SizeReal));
      if (ierr != xf_OK) return ierr;
    }
    
    // average NodeH at nodes
    for (i=0; i<Mesh->nNode; i++)
      for (k=0; k<dim2; k++) NodeH[i][k] /= ((real) NodeC[i]);
  
  
    // loop over self elements and fill in EM
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      // pull off Q1 nodes
      ierr = xf_Error(xf_Q1Nodes(Mesh->ElemGroup[egrp].QBasis, Mesh->ElemGroup[egrp].QOrder, &nn, nvec));
      if (ierr != xf_OK) return ierr;
    
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
	Node = Mesh->ElemGroup[egrp].Node[elem];
      
	EH = EM->GenArray[egrp].rValue[elem];
      
	// Set EH using NodeH
	for (i=0; i<nn; i++){
	  iglob = Node[nvec[i]];
	  for (k=0; k<dim2; k++)
	    EH[i*dim2 + k] = NodeH[iglob][k];
	}
      
      } // elem
    
    } // egrp
  
  
    // exchange halo information if parallel
    ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
    if (ierr != xf_OK) return ierr;
  
    if (nProc > 1){
      // begin communication of halo data
      ierr = xf_Error(xf_HaloExchangeVectorBegin(EM));
      if (ierr != xf_OK) return ierr;
      // end communication of halo data
      ierr = xf_Error(xf_HaloExchangeVectorEnd(EM));
      if (ierr != xf_OK) return ierr;
    }
  }  
  
  xf_Release2((void **) NodeH);
  xf_Release ((void  *) NodeC);
  
  return xf_OK;
}






/******************************************************************/
//   FUNCTION Definition: xf_ElemSize
int 
xf_ElemSize(xf_All *All, int egrp, int elem, xf_Vector *EG, real *h)
{
  int ierr, dim;
  real ElemVol, SurfArea, fac;
  enum xfe_ShapeType Shape;
  
  ElemVol  = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
  SurfArea = EG->GenArray[egrp].rValue[elem][xfe_EGSurfArea];
  
  ierr = xf_Error(xf_Basis2Shape(All->Mesh->ElemGroup[egrp].QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  
  switch (Shape){
    case xfe_Segment:       fac = 2.0; break;
    case xfe_Triangle:      fac = 2.0; break;
    case xfe_Quadrilateral: fac = 4.0; break;
    case xfe_Tetrahedron:   fac = 3.0; break;
    case xfe_Hexahedron:    fac = 6.0; break;
    default: return xf_Error(xf_NOT_SUPPORTED); break;
  }
  
  if (SurfArea <= 0.){ // in case SurfArea was not set
    dim = All->Mesh->Dim;
    if (dim == 1) (*h) = ElemVol;
    else if (dim == 2) (*h) = sqrt(ElemVol);
    else (*h) = pow(ElemVol, 1.0/((real) dim));
  }
  else{
    (*h) = fac*ElemVol/SurfArea;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetRefineCoords
int 
xf_GetRefineCoords(enum xfe_ShapeType Shape, int p, int *nnode, 
                   real **pcoord, int *nsplit, int **pvsplit,
                   int *nbound, int **pvbound)
{
  int ierr, i, j, k, l, nn, ns, ix, iy, iz;
  int d0, d1, d, dx, dy, dz;
  int *v, n[8];
  int H2T[6][4] = {{0,1,2,4},
    {1,2,4,5},
    {6,5,4,2},
    {1,3,2,5},
    {2,5,3,6},
    {7,5,6,3}};
  
  if (p <= 0) return xf_Error(xf_INPUT_ERROR); 
  
  switch (Shape){
    case xfe_Triangle:
      nn = (*nnode) = (p+1)*(p+2)/2;
      ns = (*nsplit) = p*p;
      ierr = xf_Error(xf_ReAlloc((void **) pcoord, 2*nn, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc((void **) pvsplit, 3*ns, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_LagrangeNodesEqual(Shape, p, (*pcoord)));
      if (ierr != xf_OK) return ierr;
      k = 0; i = 0;
      for (iy=0; iy<p; iy++){
        dy = p+1-iy;
        for (ix=0; ix<(p-iy); ix++){
          v = (*pvsplit) + 3*k;
          v[0] = i; v[1] = i+1; v[2] = i+dy;
          k++;
          if (ix != p-iy-1){
            v = (*pvsplit) + 3*k;
            v[0] = i+1; v[1] = i+dy+1; v[2] = i+dy;
            k++;
          }
          i++;
        } // ix
        i++;
      } // iy
      
      if (nbound != NULL){
        (*nbound) = 3*p;
        ierr = xf_Error(xf_ReAlloc((void **) pvbound, (*nbound), sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        v = (*pvbound) + 0;
        v[0] = p;
        d0 = p; d1 = -1;
        for (i=1, d=d0; i<p; i++, d+=d1) v[i] = v[i-1]+d;
        v = (*pvbound) + p;
        v[0] = (p+1)*(p+2)/2-1;
        d0 = -2; d1 = -1;
        for (i=1, d=d0; i<p; i++, d+=d1) v[i] = v[i-1]+d;
        v = (*pvbound) + 2*p;
        v[0] = 0;
        d0 = 1; d1 = 0;
        for (i=1, d=d0; i<p; i++, d+=d1) v[i] = v[i-1]+d;
      }
      
      break;
      case xfe_Quadrilateral:
      nn = (*nnode) = (p+1)*(p+1);
      ns = (*nsplit) = 2*p*p;
      ierr = xf_Error(xf_ReAlloc((void **) pcoord, 2*nn, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc((void **) pvsplit, 3*ns, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_LagrangeNodesEqual(Shape, p, (*pcoord)));
      if (ierr != xf_OK) return ierr;
      for (iy=0, k=0, i=0; iy<p; iy++, i++){
        for (ix=0; ix<p; ix++, i++){
          v = (*pvsplit) + 3*k;
          v[0] = i; v[1] = i+1; v[2] = i+p+1;
          v = (*pvsplit) + 3*(k+1);
          v[0] = i+1; v[1] = i+p+2; v[2] = i+p+1;
          k+=2;
        } // ix
      } // iy
      
      if (nbound != NULL){
        (*nbound) = 4*p;
        ierr = xf_Error(xf_ReAlloc((void **) pvbound, (*nbound), sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        v = (*pvbound) + 0;
        v[0] = 0;
        for (i=1; i<p; i++) v[i] = v[i-1]+1;
        v = (*pvbound) + p;
        v[0] = p;
        for (i=1; i<p; i++) v[i] = v[i-1]+p+1;
        v = (*pvbound) + 2*p;
        v[0] = (p+2)*p;
        for (i=1; i<p; i++) v[i] = v[i-1]-1;
        v = (*pvbound) + 3*p;
        v[0] = (p+1)*p;
        for (i=1; i<p; i++) v[i] = v[i-1]-p-1;
        
      }
      
      break;
      
      case xfe_Tetrahedron: // refined into more tets
      nn = (*nnode) = (p+1)*(p+2)*(p+3)/6;
      ierr = xf_Error(xf_ReAlloc((void **) pcoord, 3*nn, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_LagrangeNodesEqual(Shape, p, (*pcoord)));
      if (ierr != xf_OK) return ierr;
      
      if ((nsplit != NULL) && (pvsplit != NULL)){
        ns = (*nsplit) = p*p*p;
        ierr = xf_Error(xf_ReAlloc((void **) pvsplit, 4*ns, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        k = 0; i = 0;
        for (iz=0; iz<p; iz++,i++){
          for (iy=0; iy<(p-iz); iy++,i++){ 
            dz = (p-iz+1)*(p-iz+2)/2 - iy;
            for (ix=0; ix<(p-iz-iy); ix++,i++){
              dy = p+1-iz-iy;
              dx = 1;
              v = (*pvsplit) + 4*k;
              v[0]=i; v[1]=i+dx; v[2]=i+dy; v[3]=i+dz; k++;
              if ((ix+iy+iz) < (p-1)){ // not at endpoint
                v = (*pvsplit) + 4*k;
                v[0]=i+dx; v[1]=i+dx+dz; v[2]=i+dy; v[3]=i+dz; k++;
                v = (*pvsplit) + 4*k;
                v[0]=i+dy; v[1]=i+dx+dz; v[2]=i+dy+dz-1; v[3]=i+dz; k++;
                v = (*pvsplit) + 4*k;
                v[0]=i+dx; v[1]=i+dx+dy; v[2]=i+dy; v[3]=i+dx+dz; k++;
                v = (*pvsplit) + 4*k;
                v[0]=i+dy; v[1]=i+dx+dz; v[2]=i+dx+dy; v[3]=i+dy+dz-1; k++;
              }
              if ((ix+iy+iz) < (p-2)){ // deep in center; one more tet
                v = (*pvsplit) + 4*k;
                v[0]=i+dx+dy; v[1]=i+dx+dy+dz-1; v[2]=i+dx+dz; v[3]=i+dy+dz-1; k++;
              }
            } // ix
          } // iy
        } // iz
      }
      
      // NEED CODE HERE FOR BOUNDARY FACES
      
      break;
      
      
      case xfe_Hexahedron: // refined into tets
      nn = (*nnode) = (p+1)*(p+1)*(p+1);
      ierr = xf_Error(xf_ReAlloc((void **) pcoord, 3*nn, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_LagrangeNodesEqual(Shape, p, (*pcoord)));
      if (ierr != xf_OK) return ierr;
      
      if ((nsplit != NULL) && (pvsplit != NULL)){
        ns = (*nsplit) = 6*p*p*p;
        ierr = xf_Error(xf_ReAlloc((void **) pvsplit, 4*ns, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        dx = 1;
        dy = p+1;
        dz = (p+1)*(p+1);
        for (iz=0, k=0; iz<p; iz++){
          for (iy=0; iy<p; iy++){
            for (ix=0; ix<p; ix++){
              i = dz*iz+dy*iy+dx*ix; // corner node number
              n[0] = i; n[1] = i+dx; n[2] = i+dy; n[3] = i+dx+dy; 
              for (j=0; j<4; j++) n[4+j] = n[j]+dz;
              v = (*pvsplit) + 4*k;
              for (j=0; j<6; j++)
                for (l=0; l<4; l++)
                  v[4*j+l] = n[H2T[j][l]];
              k += 6;
            }
          }
        } // iz
      }
      
      // NEED CODE HERE FOR BOUNDARY FACES
      
      break;
      case xfe_Segment:
      nn = (*nnode) = (p+1);
      ns = (*nsplit) = p;
      ierr = xf_Error(xf_ReAlloc((void **) pcoord, nn, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc((void **) pvsplit, 2*ns, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_LagrangeNodesEqual(Shape, p, (*pcoord)));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<p; i++){
        v = (*pvsplit) + 2*i;
        v[0] = i; v[1] = i+1;
      } // i
      break;
      
      default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_GetRefineCoordsOnFace
int 
xf_GetRefineCoordsOnFace(enum xfe_ShapeType Shape, int face, int p, 
                         int *nnode, real **pcoord, int *nsplit, 
                         int **pvsplit, int *nbound, int **pvbound)
{
  int ierr, dim;
  enum xfe_ShapeType FShape;
  real *fcoord = NULL;
  
  if (p <= 0) return xf_Error(xf_INPUT_ERROR); 
  
  // get shape on face
  ierr = xf_Error(xf_FaceShape(Shape, face, &FShape));
  if (ierr != xf_OK) return ierr;
  
  // refinement coords on face
  ierr = xf_Error(xf_GetRefineCoords(FShape, p, nnode, &fcoord, nsplit, 
                                     pvsplit, nbound, pvbound));
  if (ierr != xf_OK) return ierr;
  
  // get dim
  ierr = xf_Error(xf_Shape2Dim(Shape, &dim));
  if (ierr != xf_OK) return ierr;
  
  // allocate (*pcoord)
  ierr = xf_Error(xf_ReAlloc( (void **) pcoord, (*nnode)*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // convert to elem ref coords (orientation is not relevant)
  ierr = xf_Error(xf_RefFace2Elem(Shape, face, 0, (*nnode), fcoord, (*pcoord)));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) fcoord);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_EdgeHashAdd
int
xf_EdgeHashAdd(int n0, int n1, int *node2ehash, 
               xf_EdgeHash *hash, int *nhash)
{
  int ihash, prev;
  
  if (n0 > n1){
    printf("Error, n0 must not be > n1 when adding to hash.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  if (node2ehash[n0] == -1){
    node2ehash[n0] = ihash = (*nhash);
    (*nhash)++;
  }
  else{
    ihash = node2ehash[n0];
    
    do{
      if (hash[ihash].n1 == n1) return xf_OK; // already exists
      prev = ihash;
    } while ( (ihash = hash[ihash].next) != -1);
    
    ihash = (*nhash);
    hash[prev].next = ihash;
    (*nhash)++;
  }
  
  hash[ihash].n0 = n0;
  hash[ihash].n1 = n1;
  hash[ihash].next = -1;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_EdgeHashCheck
int
xf_EdgeHashCheck(int n0, int n1, int *node2ehash, 
                 xf_EdgeHash *hash, int *ihash)
{
  (*ihash) = node2ehash[n0];
  while ((*ihash) >= 0){
    if (hash[(*ihash)].n1 == n1) return xf_OK;
    (*ihash) = hash[(*ihash)].next;
  }
  return xf_NOT_FOUND;
}



/******************************************************************/
//   FUNCTION Definition: xf_FaceHashAdd
static int
xf_FaceHashAdd(int n0, int n1, int n2, int *node2fhash, 
               xf_FaceHash *hash, int *nhash)
{
  /*
   
   PURPOSE: Adds a face (n0-n1-n2) to the hash
   
   INPUTS:
   
   n0, n1: nodes to add
   node2fhash, hash: hash info
   
   OUTPUTS: 
   
   nhash: resulting number of faces in hash
   
   RETURNS: Error Code
   
   */
  int k;
  int ihash, prev;
  
  // sort n0, n1, n2
  if (n0 > n1) swap(n0, n1, k);
  if (n1 > n2) swap(n1, n2, k);
  if (n0 > n1) swap(n0, n1, k);
  
  
  if (node2fhash[n0] == -1){
    node2fhash[n0] = ihash = (*nhash);
    (*nhash)++;
  }
  else{
    ihash = node2fhash[n0];
    
    do{
      if ((hash[ihash].n1 == n1) &&
          (hash[ihash].n2 == n2)) return xf_OK; // already exists
      prev = ihash;
    } while ( (ihash = hash[ihash].next) != -1);
    
    ihash = (*nhash);
    hash[prev].next = ihash;
    (*nhash)++;
  }
  
  hash[ihash].n0 = n0;
  hash[ihash].n1 = n1;
  hash[ihash].n2 = n2;
  hash[ihash].next = -1;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_TriangleNormalOrient
static void 
xf_TriangleNormalOrient(int *vnode, int n3, real d3, real *coord, real *xref)
{
  /*
   This function orients the nodes of a triangle in 3-space such that
   the triangle normal is in a consistent direction in reference to
   another point (not on the triangle).  xref is exclusively for n3.
   */
  int d;
  real *x0, *x1, *x2, n[3];
  real v1[3], v2[3], vout[3], dp, sg;
  
  // compute current normal, n
  x0 = coord + 3*vnode[0];
  x1 = coord + 3*vnode[1];
  x2 = coord + 3*vnode[2];
  
  for (d=0; d<3; d++) v1[d] = x1[d] - x0[d];
  for (d=0; d<3; d++) v2[d] = x2[d] - x0[d];
  
  xf_CrossProduct(v1, v2, n); // normal = v1 x v2
  
  // compute vout = (x(n3) - x(vnode(0)))*sign(d3)
  sg = ((d3 >= 0.0) ? 1.0 : -1.0);
  for (d=0; d<3; d++) vout[d] = (xref[3*n3+d] - x0[d])*sg;
  
  // dot n with vout : reorient if negative
  xf_DotProduct(vout, n, 3, &dp);
  
  if (dp < 0.0) swap(vnode[1], vnode[2], d);
  
}

/******************************************************************/
//   FUNCTION Definition: xf_IntersectElemWithPlane
int 
xf_IntersectElemWithPlane(xf_Mesh *Mesh, int egrp, int elem, int refine, 
                          const real *CutPlane, const real *distin, 
                          enum xfe_Bool PointsChanged,  
                          xf_BasisData **pPhiData, int *pnnode, real **pcoord, 
                          int *pntri, int **pvtri, int *pnbound, int **pvbound)
{
  int ierr, dim, k, i, j, s;
  int itet, nitet, *i2tet;
  int nnode, ntet, *vtet, *n, *Node;
  int e, n0, n1, n3, iedge, nedge;
  int nin, nie, vin[4], vie[6], eglob[6];
  int m, m1, m2, nface, n01[2];
  int *node2i, *edge2i, *node2ehash, *node2fhash;
  int itri, *edge2b = NULL;
  enum xfe_Bool hit;
  enum xfe_ShapeType Shape;
  real *xref, *xglob, *dist, d[4];
  real d0, d1, fac, *x0, *x1;
  xf_EdgeHash *edgehash;
  xf_FaceHash *facehash;
  
  
  if ((dim = Mesh->Dim) == 2){
    
    if (CutPlane == NULL) return xf_Error(xf_INPUT_ERROR);
    
    // in 2D, we need a line-element intersection; relatively straightforward
    Node  = Mesh->ElemGroup[egrp].Node[elem];
    nnode = Mesh->ElemGroup[egrp].nNode;
    
    // allocate pcoord
    (*pnnode) = 0;
    ierr = xf_Error(xf_ReAlloc((void **) pcoord, 4, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // reference-space coords of nodes
    xref = NULL;
    ierr = xf_Error(xf_LagrangeNodes(Mesh->ElemGroup[egrp].QBasis, 1, NULL, NULL, &xref));
    if (ierr != xf_OK) return ierr;
    
    // first check for node intersections
    ierr = xf_Error(xf_Alloc((void **) &dist, nnode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    for (k=0,nin=0; k<nnode; k++){
      xglob = Mesh->Coord[Node[k]];
      for (i=0, dist[k]=CutPlane[2]; i<2; i++) dist[k] += CutPlane[i]*xglob[i];
      if (fabs(dist[k]) < MEPS){
        vin[nin++] = k;
        dist[k] = 0.0;
        for (i=0; i<2; i++) (*pcoord)[2*(*pnnode)+i] = xref[3*k+i];
        (*pnnode)++;
      }
    } // k
    
    
    // number of faces (= edges)
    ierr = xf_Error(xf_Basis2nFace(Mesh->ElemGroup[egrp].QBasis, &nedge));
    if (ierr != xf_OK) return ierr;
    
    // next check for edge intersections
    for (iedge=0; iedge<nedge; iedge++){
      ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, 1, iedge, &k, n01));
      if (ierr != xf_OK) return ierr;
      n0 = n01[0]; n1 = n01[1];
      if ((dist[n0] == 0.) || (dist[n1] == 0.)) continue;
      if (dist[n0]*dist[n1] < 0.){
        d0 = fabs(dist[n0]);
        d1 = fabs(dist[n1]);
        fac = d0/(d0+d1);
        x0 = xref + 2*n0;
        x1 = xref + 2*n1;
        for (i=0; i<2; i++)
          (*pcoord)[2*(*pnnode)+i] = fac*x1[i] + (1.-fac)*x0[i];
        (*pnnode)++;
      }
    }
    
    // must have hit 2 nodes, 1 node + 1 edge, or 2 edges
    if ((*pnnode) != 2){
      xf_printf("nnode = %d\n", (*pnnode)); fflush(stdout);
      return xf_Error(xf_INTERSECT_ERROR);
    }
    
    
    xf_Release ( (void *) dist);
    xf_Release ( (void *) xref);
    
    return xf_OK;
  } // end 2D case  
  
  
  
  // pull off element shape
  ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  // get refined coords
  xref = NULL;
  vtet = NULL;
  ierr = xf_Error(xf_GetRefineCoords(Shape, refine, &nnode, &xref, 
                                     &ntet, &vtet, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // obtain global coords: xref -> xglob
  ierr = xf_Error(xf_Alloc((void **) &xglob, 3*nnode, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, pPhiData, PointsChanged,
                                  nnode, xref, xglob));
  if (ierr != xf_OK) return ierr;
  
  
  // allocate more than enough space for (*pcoord) and (*pvtri)
  ierr = xf_Error(xf_ReAlloc((void **) pcoord, 3*4*ntet, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc((void **) pvtri, 3*2*ntet, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // allocate memory for node-to-intersection index
  ierr = xf_Error(xf_Alloc((void **) &node2i, nnode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<nnode; k++) node2i[k] = -1;
  
  (*pnnode) = 0; // total number of nodes we're going to return
  
  // mark each node as above/below/on plane (calculate distance)
  ierr = xf_Error(xf_Alloc((void **) &dist, nnode, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<nnode; k++){
    if (CutPlane == NULL){
      dist[k] = distin[k];
    }
    else{
      for (i=0, dist[k]=CutPlane[3]; i<3; i++) dist[k] += CutPlane[i]*xglob[3*k+i];
    }
    if (fabs(dist[k]) < MEPS) {
      dist[k] = 0.0;
      for (i=0; i<3; i++) (*pcoord)[3*(*pnnode)+i] = xref[3*k+i];
      node2i[k] = (*pnnode)++;
    }
  } // k
  
  // loop over ntet sub-tetrahedra and flag possible intersections
  ierr = xf_Error(xf_Alloc((void **) &i2tet, ntet, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (itet=0, nitet=0; itet<ntet; itet++){
    n = vtet + 4*itet;
    for (k=0; k<4; k++) d[k] = dist[n[k]];
    
    // check for hit
    if (d[0] == 0.) 
      hit = xfe_True;
    else
      for (k=1, s = sign(d[0]), hit = xfe_False; k<4; k++)
        if ((sign(d[k]) != s) || (d[k] == 0.)){
          hit = xfe_True;
          break;
        }
    
    if (hit) i2tet[nitet++] = itet;
  } // itet
  
  // allocate edge hash: node2ehash + edgehash
  ierr = xf_Error(xf_Alloc( (void **) &node2ehash, nnode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nnode; i++) node2ehash[i] = -1;
  ierr = xf_Error(xf_Alloc( (void **) &edgehash, 6*ntet, sizeof(xf_EdgeHash)));
  if (ierr != xf_OK) return ierr;
  
  // build an edge hash using intersected elements
  for (i=0, nedge=0; i<nitet; i++){
    itet = i2tet[i];
    n = vtet + 4*itet;
    for (e=0; e<6; e++){
      n0 = n[E2N[e][0]]; n1 = n[E2N[e][1]];
      if (n1 < n0) swap(n0,n1, k);
      ierr = xf_Error(xf_EdgeHashAdd(n0, n1, node2ehash, edgehash, &nedge));
      if (ierr != xf_OK) return ierr;
    } // e
  } // i
  
  
  // allocate memory for edge intersection index
  ierr = xf_Error(xf_Alloc((void **) &edge2i, nedge, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // calculate intersection for each edge in hash
  for (iedge=0; iedge<nedge; iedge++){
    edge2i[iedge] = -1;
    n0 = edgehash[iedge].n0; n1 = edgehash[iedge].n1;
    if ((dist[n0]*dist[n1] <= 0.0) && (node2i[n0] == -1) && (node2i[n1] == -1)){
      d0 = fabs(dist[n0]);
      d1 = fabs(dist[n1]);
      fac = d0/(d0+d1);
      x0 = xref + 3*n0;
      x1 = xref + 3*n1;
      for (i=0; i<3; i++)
        (*pcoord)[3*(*pnnode)+i] = fac*x1[i] + (1.-fac)*x0[i];
      edge2i[iedge] = (*pnnode)++;
    }
  } // iedge
  
  
  // allocate face hash: node2fhash + facehash
  ierr = xf_Error(xf_Alloc( (void **) &node2fhash, nnode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nnode; i++) node2fhash[i] = -1;
  ierr = xf_Error(xf_Alloc( (void **) &facehash, 2*nedge, sizeof(xf_FaceHash)));
  if (ierr != xf_OK) return ierr;
  
  // loop over intersected tets, add to vtri
  (*pntri) = 0;
  nface = 0;
  for (i=0; i<nitet; i++){
    
    itet = i2tet[i];
    n = vtet + 4*itet;
    nin = 0; // number of node intersections
    nie = 0; // number of edge intersections
    for (k=0, n3=n[0]; k<4; k++){
      if (node2i[n[k]] >= 0) vin[nin++] = k;
      else n3 = n[k]; // n3 will be a node not on the surface
    }
    for (e=0; e<6; e++){
      n0 = n[E2N[e][0]]; n1 = n[E2N[e][1]];
      if (n1 < n0) swap(n0, n1, k);
      ierr = xf_Error(xf_EdgeHashCheck(n0, n1, node2ehash, edgehash, &iedge));
      if (ierr != xf_OK) return ierr;
      eglob[e] = iedge;
      if (edge2i[iedge] >= 0) vie[nie++] = e;
    }
    if (nin == 3){
      // add vin[:] to face hash; if not already in, add vin[:] to vtri
      k = nface;
      ierr = xf_Error(xf_FaceHashAdd(n[vin[0]], n[vin[1]], n[vin[2]], 
                                     node2fhash, facehash, &nface));
      if (ierr != xf_OK) return ierr;
      if (k != nface){ // not added yet
        for (k=0; k<3; k++)
          (*pvtri)[3*(*pntri) + k] = node2i[n[vin[k]]];
        xf_TriangleNormalOrient((*pvtri) + 3*(*pntri), n3, dist[n3], (*pcoord), xref);
        (*pntri)++;
      }
    }
    else if (nie == 0){
      continue; // no edge intersections for this tet
    }
    else if (nie == 1){ // one edge, 2 node intersections
      if (nin != 2) return xf_Error(xf_INTERSECT_ERROR);
      for (k=0; k<2; k++)
        (*pvtri)[3*(*pntri) + k] = node2i[n[vin[k]]];
      (*pvtri)[3*(*pntri) + 2] = edge2i[eglob[vie[0]]];
      xf_TriangleNormalOrient((*pvtri) + 3*(*pntri), n3, dist[n3], (*pcoord), xref);
      //for (k=0; k<3; k++) (*pvtri)[3*(*pntri) + k] = 0; // TEMPORARY
      (*pntri)++;
    }
    else if (nie == 2){ // 2 edge, 1 node intersections
      if (nin != 1) return xf_Error(xf_INTERSECT_ERROR);
      for (k=0; k<2; k++)
        (*pvtri)[3*(*pntri) + k] = edge2i[eglob[vie[k]]];
      (*pvtri)[3*(*pntri) + 2] = node2i[n[vin[0]]];
      xf_TriangleNormalOrient((*pvtri) + 3*(*pntri), n3, dist[n3], (*pcoord), xref);
      //for (k=0; k<3; k++) (*pvtri)[3*(*pntri) + k] = 0; // TEMPORARY
      (*pntri)++;
    }
    else if (nie == 3){ // 3 edge intersections
      if (nin != 0) return xf_Error(xf_INTERSECT_ERROR);
      for (k=0; k<3; k++)
        (*pvtri)[3*(*pntri) + k] = edge2i[eglob[vie[k]]];
      xf_TriangleNormalOrient((*pvtri) + 3*(*pntri), n3, dist[n3], (*pcoord), xref);
      //for (k=0; k<3; k++) (*pvtri)[3*(*pntri) + k] = 0; // TEMPORARY
      (*pntri)++;
    }
    else if (nie == 4){ // two triangles
      e = vie[0];
      for (m = 1; m<4; m++) // m will be opposite edge
        if (OE[e] == vie[m]) break;
      if ((OE[e] != vie[m]) && (m == 3)) return xf_Error(xf_INTERSECT_ERROR);
      m1 = m%3+1;
      m2 = (m+1)%3+1;
      (*pvtri)[3*(*pntri) + 0] = edge2i[eglob[vie[0 ]]];
      (*pvtri)[3*(*pntri) + 1] = edge2i[eglob[vie[m ]]];
      (*pvtri)[3*(*pntri) + 2] = edge2i[eglob[vie[m1]]];
      xf_TriangleNormalOrient((*pvtri) + 3*(*pntri), n3, dist[n3], (*pcoord), xref);
      //for (k=0; k<3; k++) (*pvtri)[3*(*pntri) + k] = 0; // TEMPORARY
      (*pntri)++;
      (*pvtri)[3*(*pntri) + 0] = edge2i[eglob[vie[0 ]]];
      (*pvtri)[3*(*pntri) + 1] = edge2i[eglob[vie[m ]]];
      (*pvtri)[3*(*pntri) + 2] = edge2i[eglob[vie[m2]]];
      xf_TriangleNormalOrient((*pvtri) + 3*(*pntri), n3, dist[n3], (*pcoord), xref);
      //for (k=0; k<3; k++) (*pvtri)[3*(*pntri) + k] = 0; // TEMPORARY
      (*pntri)++;
    }
    else
      return xf_Error(xf_INTERSECT_ERROR);
    
  } // i
  
  // reallocate memory
  ierr = xf_Error(xf_ReAlloc((void **) pcoord, 3*(*pnnode), sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc((void **) pvtri, 3*(*pntri), sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  if (pnbound != NULL) (*pnbound) = 0;
  
  // loops on boundary -- for now, segments on boundary
  if ((pnbound != NULL) && (pvbound != NULL)){
    
    // re-allocate edge hash: node2ehash + edgehash
    ierr = xf_Error(xf_ReAlloc( (void **) &node2ehash, (*pnnode), sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<(*pnnode); i++) node2ehash[i] = -1;
    ierr = xf_Error(xf_ReAlloc( (void **) &edgehash, 3*(*pntri), sizeof(xf_EdgeHash)));
    if (ierr != xf_OK) return ierr;
    
    // flag to determine which edges are boundaries
    ierr = xf_Error(xf_Alloc( (void **) &edge2b, 3*(*pntri), sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<3*(*pntri); i++) edge2b[i] = 1; // start out as boundary
    
    // build edge hash
    nedge = 0;
    for (itri=0; itri<(*pntri); itri++){
      n = (*pvtri) + 3*itri;
      for (e=0; e<3; e++){
        n0 = n[(e+1)%3]; n1 = n[(e+2)%3];
        if (n1 < n0) swap(n0,n1, k);
        k = nedge;
        ierr = xf_Error(xf_EdgeHashAdd(n0, n1, node2ehash, edgehash, &nedge));
        if (ierr != xf_OK) return ierr;
        if (nedge == k){ // no edge was added; means shared by 2 triangles
          ierr = xf_Error(xf_EdgeHashCheck(n0, n1, node2ehash, edgehash, &k));
          if (ierr != xf_OK) return ierr;
          edge2b[k] = 0; // not a boundary
        }
      } // e
    } // itri
    
    // count how many boundary edges we have
    (*pnbound) = 0;
    for (iedge=0; iedge<nedge; iedge++)
      if (edge2b[iedge] == 1) (*pnbound)++;
    
    // allocate space for boundary edge list
    ierr = xf_Error(xf_ReAlloc( (void **) pvbound, 2*(*pnbound), sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    // add edge boundaries
    (*pnbound) = 0;
    for (iedge=0; iedge<nedge; iedge++)
      if (edge2b[iedge] == 1){
        (*pvbound)[2*(*pnbound)+0] = edgehash[iedge].n0;
        (*pvbound)[2*(*pnbound)+1] = edgehash[iedge].n1;	
        (*pnbound)++;
      }
    
    // this stores the total number of nodes
    (*pnbound) *= 2;
    
    xf_Release( (void *) edge2b);
    
  } // end boundary
  
  
  // release memory
  xf_Release( (void *) xref);
  xf_Release( (void *) xglob);
  xf_Release( (void *) dist);
  xf_Release( (void *) vtet);
  xf_Release( (void *) node2i);
  xf_Release( (void *) edge2i);
  xf_Release( (void *) i2tet);
  xf_Release( (void *) node2ehash);
  xf_Release( (void *) edgehash);
  xf_Release( (void *) node2fhash);
  xf_Release( (void *) facehash);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CheckVolumes
int 
xf_CheckVolumes(xf_Mesh *Mesh, int *OrderVec, enum xfe_Bool ErrorOutFlag,
                int nstoremax, int *pnstore, real *xstore)
{
  int ierr, d;
  int egrp, elem, face;
  int QuadOrder, Order;
  int nq, nqFace, iq;
  int orient, pnq, pnqFace;
  int istore;
  enum xfe_Bool QuadChanged;
  real minJ, maxJ, val;
  real *xq, *xelem, *xqFace, *xglob;
  xf_QuadData *QuadData, *QuadDataFace;
  xf_JacobianData *JData, *JDataFace;
  xf_Face Face;
  xf_IFace IFace;
  xf_BFace BFace;
  xf_BasisData *PhiData;
  
  QuadData      = NULL;
  QuadDataFace  = NULL;
  PhiData   = NULL;
  JData     = NULL;
  JDataFace = NULL;
  xelem     = NULL;
  xglob     = NULL;
  
  pnq     = -1;
  pnqFace = -1;
  minJ =  1e30;
  maxJ = -1e30;
  
  istore = 0;
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    Order = ((OrderVec != NULL) ? OrderVec[egrp] : 1);
    
    ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order, &QuadOrder));
    if (ierr != xf_OK) return ierr;
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      /* Pull off quad points for the element; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **) &xglob, Mesh->Dim*nq, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        pnq = nq;
      }
      
      // obtain glob coordinates
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &PhiData, xfe_True,
                                      nq, xq, xglob));
      if (ierr != xf_OK) return ierr;
      
      /* Compute geometry Jacobian; if not constant, compute at quad
       points.  Note if jacobian is constant, only one Jacobian will
       be computed/returned. */
      ierr = xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData);
      if (ierr == xf_NEGATIVE_JACOBIAN){
        for (iq=0; iq<JData->nq; iq++)
          if (JData->detJ[iq] <= 0.0){
            xf_printf("Elem interior:\n");
            xf_printf(" |J| = %.6E, egrp=%d, elem=%d\n xref = ", JData->detJ[iq], egrp, elem);
            for (d=0; d<Mesh->Dim; d++) xf_printf(" %.5E", xq[iq*Mesh->Dim+d]);
            xf_printf("\n");
            xf_printf(" xglob = ");
            for (d=0; d<Mesh->Dim; d++) xf_printf(" %.5E", xglob[iq*Mesh->Dim+d]);
            xf_printf("\n");
            if (istore < nstoremax){
              for (d=0; d<Mesh->Dim; d++)
                xstore[istore*Mesh->Dim + d] = xglob[iq*Mesh->Dim+d];
              istore++;
            }
          }
        if (ErrorOutFlag) return xf_Error(ierr);
      }
      else if (ierr != xf_OK) return ierr;
      
      for (iq=0; iq<JData->nq; iq++){
        val = JData->detJ[iq];
        minJ = min(minJ, val); maxJ = max(maxJ, val);
      }
      
      /* Compute geometry Jacobian on faces as well. */
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
        // quad points on face
        ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face, QuadOrder, 
                                    &QuadDataFace, &QuadChanged));
        if (ierr != xf_OK) return ierr;
        
        nqFace = QuadDataFace->nquad;
        xqFace = QuadDataFace->xquad;
        
        Face = Mesh->ElemGroup[egrp].Face[elem][face];
        
        if (Face.Group == xf_INTERIORFACE){ // Interior Face
          IFace = Mesh->IFace[Face.Number];
          orient = ((IFace.ElemGroupL==egrp)&&(IFace.ElemL==elem) ? IFace.OrientL : IFace.OrientR);
        }
        else if (Face.Group >= 0) { // Boundary Face
          BFace = Mesh->BFaceGroup[Face.Group].BFace[Face.Number];
          orient = BFace.Orient;
        }
        else return xf_Error(xf_NOT_SUPPORTED);
        
        if (nqFace > pnqFace){
          ierr = xf_Error(xf_ReAlloc( (void **) &xelem, Mesh->Dim*nqFace, sizeof(real)));
          if (ierr != xf_OK) return ierr;
          pnqFace = nqFace;
        }
        
        // convert xqFace to xelem
        ierr = xf_Error(xf_RefFace2Interpol(Mesh, egrp, elem, face, orient,
                                            nqFace, xqFace, xelem));
        if (ierr != xf_OK) return ierr;
        
        // obtain glob coordinates
        ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &PhiData, xfe_True,
                                        nqFace, xelem, xglob));
        if (ierr != xf_OK) return ierr;
        
        /* Compute geometry Jacobian */
        ierr = xf_ElemJacobian(Mesh, egrp, elem, nqFace, xelem, xfb_detJ, xfe_True, &JDataFace);
        if (ierr == xf_NEGATIVE_JACOBIAN){
          for (iq=0; iq<min(nqFace, JDataFace->nq); iq++)
            if (JDataFace->detJ[iq] <= 0.0){
              xf_printf("On face:\n");
              xf_printf(" |J| = %.6E, egrp=%d, elem=%d, face=%d\n xref = ", 
                        JDataFace->detJ[iq], egrp, elem, face);
              for (d=0; d<Mesh->Dim; d++) xf_printf(" %.5E", xelem[iq*Mesh->Dim+d]);
              xf_printf("\n");
              xf_printf(" xglob = ");
              for (d=0; d<Mesh->Dim; d++) xf_printf(" %.5E", xglob[iq*Mesh->Dim+d]);
              xf_printf("\n");
              if (istore < nstoremax){
                for (d=0; d<Mesh->Dim; d++)
                  xstore[istore*Mesh->Dim + d] = xglob[iq*Mesh->Dim+d];
                istore++;
              }
            }
          if (ErrorOutFlag) return xf_Error(ierr);
        }
        else if (ierr != xf_OK) return ierr;
        
        for (iq=0; iq<JDataFace->nq; iq++){
          val = JDataFace->detJ[iq];
          minJ = min(minJ, val); maxJ = max(maxJ, val);
        }
        
      } // face
      
    } // elem
    
  } // egrp
  
  if (pnstore != NULL) (*pnstore) = istore;
  
  xf_printf("Volume check finished: minJ = %.6E, maxJ = %.6E\n", minJ, maxJ);
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy geometry Jacobian Data on Face*/
  ierr = xf_Error(xf_DestroyJacobianData(JDataFace));
  if (ierr != xf_OK) return ierr;
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  // Only destroy QuadDataFace if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataFace));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) xelem);
  xf_Release( (void *) xglob);
  
  return xf_OK;
}

/******************************************************************/
//  FUNCTION Definition: xf_Transfinite2D
void 
xf_Transfinite2D(int dim, real *xref, real **x)
{
  int i, j, k, d, s;
  real X, Y;
  real R[9];
  real XV[3], YV[3];

  X = xref[0];
  Y = xref[1];
  
  // construct transfinite interpolation coefficients
  XV[0] = 1.-X; XV[1] = 1.; XV[2] = X;
  YV[0] = 1.-Y; YV[1] = 1.; YV[2] = Y;

  for (j=0, k=0; j<3; j++)
    for (i=0; i<3; i++)
      R[k++] = YV[j]*XV[i];
  R[4] = 0.;

  for (d=0; d<dim; d++)
    for (k=0, x[4][d]=0.0, s=-1; k<9; k++, s=-s)
      x[4][d] += R[k]*x[k][d]*s;

}


/******************************************************************/
//  FUNCTION Definition: xf_Transfinite3D
void 
xf_Transfinite3D(int dim, real *xref, real **x)
{
  int i, j, k, l, d, s;
  real X, Y, Z;
  real R[27];
  real XV[3], YV[3], ZV[3];

  X = xref[0];
  Y = xref[1];
  Z = xref[2];
  
  // construct transfinite interpolation coefficients
  XV[0] = 1.-X; XV[1] = 1.; XV[2] = X;
  YV[0] = 1.-Y; YV[1] = 1.; YV[2] = Y;
  ZV[0] = 1.-Z; ZV[1] = 1.; ZV[2] = Z;

  for (k=0, l=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++)
	R[l++] = ZV[k]*YV[j]*XV[i];
  R[13] = 0.;

  for (d=0; d<dim; d++)
    for (k=0, x[13][d]=0.0, s=1; k<27; k++, s=-s)
      x[13][d] += R[k]*x[k][d]*s;

  /* Comparison to  John Burkardt's matlab code:                     */
  /*   x =        ( 1.0 - r ) * ( 1.0 - s ) * ( 1.0 - t ) * x000 ... */
  /*            -               ( 1.0 - s ) * ( 1.0 - t ) * xr00 ... */
  /*            +         r   * ( 1.0 - s ) * ( 1.0 - t ) * x100 ... */
  /*            - ( 1.0 - r )               * ( 1.0 - t ) * x0s0 ... */
  /*            +                             ( 1.0 - t ) * xrs0 ... */
  /*            -         r                 * ( 1.0 - t ) * x1s0 ... */
  /*            + ( 1.0 - r ) *         s   * ( 1.0 - t ) * x010 ... */
  /*            -                       s   * ( 1.0 - t ) * xr10 ... */
  /*            +         r   *         s   * ( 1.0 - t ) * x110 ... */
  /*            - ( 1.0 - r ) * ( 1.0 - s )               * x00t ... */
  /*            +               ( 1.0 - s )               * xr0t ... */
  /*            -         r   * ( 1.0 - s )               * x10t ... */
  /*            + ( 1.0 - r )                             * x0st ... */
  /*            +         r                               * x1st ... */
  /*            - ( 1.0 - r ) *         s                 * x01t ... */
  /*            +                       s                 * xr1t ... */
  /*            -         r   *         s                 * x11t ... */
  /*            + ( 1.0 - r ) * ( 1.0 - s ) *         t   * x001 ... */
  /*            -               ( 1.0 - s ) *         t   * xr01 ... */
  /*            +         r   * ( 1.0 - s ) *         t   * x101 ... */
  /*            - ( 1.0 - r )                   *     t   * x0s1 ... */
  /*            +                                     t   * xrs1 ... */
  /*            -         r                 *         t   * x1s1 ... */
  /*            + ( 1.0 - r ) *         s   *         t   * x011 ... */
  /*            -                       s   *         t   * xr11 ... */
  /*            +         r   *         s   *         t   * x111;    */
}




// functions for element searching
#include "xf_MeshToolsElemSearch.c"

// functions for mesh curving
#include "xf_MeshToolsCurving.c"


#if( UNIT_TEST==1 )
#include "xf_MeshTools.test.in"
#endif

