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
 FILE:  xf_AdaptHang.c
 
 This file contains functions for hanging-node (non-conforming)
 adaptation.
 
 */


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Mesh.h"
#include "xf_MeshTools.h"
#include "xf_Math.h"
#include "xf_AdaptStruct.h"
#include "xf_EqnSet.h"


/* Maximum structure sizes for refined subelements from one coarse
 element.  Used to avoid dynamically allocating and reallocating
 memory.
 */
#define ha_MAX_VERT 27
#define ha_MAX_ELEM 8
#define ha_MAX_NODEQ1 8
#define ha_MAX_IFACE 12
#define ha_MAX_BFACE 24
#define ha_MAX_EDGE 12


// turn on to 1 for debugging
#define ah_DEBUG 0


// integer pair list
typedef struct
{
  int n;      // number of pairs in list
  int n0;     // max number for allocation purposes
  int *Pairs; // 2*n entries are valid
}
xf_IntPairList;



/******************************************************************/
//   FUNCTION Definition: xf_AddIntPair
static int 
xf_AddIntPair(int i0, int i1, xf_IntPairList *pList)
{
  /*
   PURPOSE:
   
   Adds an integer pair (i0, i1) to the list (*pList)
   
   INPUTS:
   
   i0, i1 : pair of integers
   (*pList) : integer pair list so far
   
   OUTPUTS: 
   
   (*pList) : modified integer pair list
   
   RETURN:
   
   Error Code
   */
  int ierr, n, n0;
  
  if ((*pList).n >= (*pList).n0){
    n0 = max(2*(*pList).n0, (*pList).n0+10);
    ierr = xf_Error(xf_ReAlloc( (void **) &(*pList).Pairs, 2*n0, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    (*pList).n0 = n0;
  }
  
  n = (*pList).n++;
  (*pList).Pairs[2*n+0] = i0;
  (*pList).Pairs[2*n+1] = i1;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemRefEffectOnFace
static int 
xf_ElemRefEffectOnFace(enum xfe_ShapeType Shape, int elemref, int face,
                       int *pfacerefloc)
{
  /*
   PURPOSE:
   
   Returns local face refinement from point of view of element, as
   implied by the elemental refinement, elemref.
   
   INPUTS:
   
   Shape : which element shape to consider
   elemref : element refinement
   face : face on which to consider the effect (original number)
   
   OUTPUTS: 
   
   (*pfacerefloc) : implied refinement of local face
   
   RETURN:
   
   Error Code
   */
  
  // return immediately if no refinement
  if (elemref == 0){
    (*pfacerefloc) = 0;
    return xf_OK;
  }
  
  switch (Shape){
    case xfe_Segment:
      // segment refinement has no effect on face
      (*pfacerefloc) = 0;
      break;
    case xfe_Quadrilateral:
      (*pfacerefloc) = xfe_SegRefNone;
      switch (elemref){
        case xfe_QuadRefUniform:
          (*pfacerefloc) = xfe_SegRefUniform;
          break;
        case xfe_QuadRefHoriz:
          if ((face == 1) || (face == 3)) (*pfacerefloc) = xfe_SegRefUniform;
          break;
        case xfe_QuadRefVert:
          if ((face == 0) || (face == 2)) (*pfacerefloc) = xfe_SegRefUniform;
          break;
        default:
          return xf_Error(xf_INPUT_ERROR);
          break;
      }
      break;
    case xfe_Triangle:
    case xfe_Tetrahedron:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    case xfe_Hexahedron:
      (*pfacerefloc) = xfe_QuadRefNone;
      switch (elemref){
        case xfe_HexRefUniform:
          (*pfacerefloc) = xfe_QuadRefUniform;
          break;
        case xfe_HexRefSliceX:
          if (face == 0) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 1) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 3) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 5) (*pfacerefloc) = xfe_QuadRefVert;
          break;
        case xfe_HexRefSliceY:
          if (face == 0) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 2) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 4) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 5) (*pfacerefloc) = xfe_QuadRefHoriz;
          break;
        case xfe_HexRefSliceZ:
          if (face == 1) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 2) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 3) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 4) (*pfacerefloc) = xfe_QuadRefHoriz;
          break;
        case xfe_HexRefSliceXY:
          if (face == 0) (*pfacerefloc) = xfe_QuadRefUniform;
          if (face == 1) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 2) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 3) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 4) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 5) (*pfacerefloc) = xfe_QuadRefUniform;
          break;
        case xfe_HexRefSliceXZ:
          if (face == 0) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 1) (*pfacerefloc) = xfe_QuadRefUniform;
          if (face == 2) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 3) (*pfacerefloc) = xfe_QuadRefUniform;
          if (face == 4) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 5) (*pfacerefloc) = xfe_QuadRefVert;
          break;
        case xfe_HexRefSliceYZ:
          if (face == 0) (*pfacerefloc) = xfe_QuadRefVert;
          if (face == 1) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 2) (*pfacerefloc) = xfe_QuadRefUniform;
          if (face == 3) (*pfacerefloc) = xfe_QuadRefHoriz;
          if (face == 4) (*pfacerefloc) = xfe_QuadRefUniform;
          if (face == 5) (*pfacerefloc) = xfe_QuadRefHoriz;
          break;
        default:
          return xf_Error(xf_INPUT_ERROR);
          break;
      }
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_FaceRefEffectOnElem
static int 
xf_FaceRefEffectOnElem(enum xfe_ShapeType Shape, int *frefvec, int *pelemref)
{
  /*
   PURPOSE:
   
   Returns simplest element refinement option that satisfies the given
   face refinements.
   
   INPUTS:
   
   Shape : which element shape to consider
   frefvec : vector (over all original faces) of desired face refinements
   
   OUTPUTS: 
   
   (*pelemref) : simplest consistent refinement of element
   
   RETURN:
   
   Error Code
   */
  
  int ierr, k, nface;
  enum xfe_Bool trivial;
  
  // number of faces in Shape
  ierr = xf_Error(xf_Shape2nFace(Shape, &nface));
  if (ierr != xf_OK) return ierr;
  
  // return immediately if no refinement on faces
  trivial = xfe_True;
  for (k=0; k<nface; k++) trivial = (trivial && (frefvec[k] == 0));
  if (trivial){
    (*pelemref) = 0;
    return xf_OK;
  }
  
  switch (Shape){
    case xfe_Segment:
      // should not be called for an element
      return xf_Error(xf_INPUT_ERROR);
      break;
    case xfe_Quadrilateral:
      if ((frefvec[0] == 0) && (frefvec[2] == 0)) (*pelemref) = xfe_QuadRefHoriz;
      else if ((frefvec[1] == 0) && (frefvec[3] == 0)) (*pelemref) = xfe_QuadRefVert;
      else (*pelemref) = xfe_QuadRefUniform;
      break;
    case xfe_Hexahedron:
      if ((frefvec[2] == xfe_QuadRefNone ) && (frefvec[4] == xfe_QuadRefNone   ) &&
          (frefvec[0] != xfe_QuadRefVert ) && (frefvec[0] != xfe_QuadRefUniform) &&
          (frefvec[1] != xfe_QuadRefHoriz) && (frefvec[1] != xfe_QuadRefUniform) &&
          (frefvec[3] != xfe_QuadRefHoriz) && (frefvec[3] != xfe_QuadRefUniform) &&
          (frefvec[5] != xfe_QuadRefHoriz) && (frefvec[5] != xfe_QuadRefUniform))
        (*pelemref) = xfe_HexRefSliceX;
      else if ((frefvec[1] == xfe_QuadRefNone ) && (frefvec[3] == xfe_QuadRefNone   ) &&
               (frefvec[0] != xfe_QuadRefHoriz) && (frefvec[0] != xfe_QuadRefUniform) &&
               (frefvec[2] != xfe_QuadRefHoriz) && (frefvec[2] != xfe_QuadRefUniform) &&
               (frefvec[4] != xfe_QuadRefHoriz) && (frefvec[4] != xfe_QuadRefUniform) &&
               (frefvec[5] != xfe_QuadRefVert ) && (frefvec[5] != xfe_QuadRefUniform))
        (*pelemref) = xfe_HexRefSliceY;
      else if ((frefvec[0] == xfe_QuadRefNone ) && (frefvec[5] == xfe_QuadRefNone   ) &&
               (frefvec[1] != xfe_QuadRefVert ) && (frefvec[1] != xfe_QuadRefUniform) &&
               (frefvec[2] != xfe_QuadRefVert ) && (frefvec[2] != xfe_QuadRefUniform) &&
               (frefvec[3] != xfe_QuadRefVert ) && (frefvec[3] != xfe_QuadRefUniform) &&
               (frefvec[4] != xfe_QuadRefVert ) && (frefvec[4] != xfe_QuadRefUniform))
        (*pelemref) = xfe_HexRefSliceZ;
      else if ((frefvec[1] != xfe_QuadRefHoriz) && (frefvec[1] != xfe_QuadRefUniform) &&
               (frefvec[2] != xfe_QuadRefHoriz) && (frefvec[2] != xfe_QuadRefUniform) &&
               (frefvec[3] != xfe_QuadRefHoriz) && (frefvec[3] != xfe_QuadRefUniform) &&
               (frefvec[4] != xfe_QuadRefHoriz) && (frefvec[4] != xfe_QuadRefUniform))
        (*pelemref) = xfe_HexRefSliceXY;
      else if ((frefvec[0] != xfe_QuadRefVert ) && (frefvec[0] != xfe_QuadRefUniform) &&
               (frefvec[2] != xfe_QuadRefVert ) && (frefvec[2] != xfe_QuadRefUniform) &&
               (frefvec[4] != xfe_QuadRefVert ) && (frefvec[4] != xfe_QuadRefUniform) &&
               (frefvec[5] != xfe_QuadRefHoriz) && (frefvec[5] != xfe_QuadRefUniform))
        (*pelemref) = xfe_HexRefSliceXZ;
      else if ((frefvec[0] != xfe_QuadRefHoriz) && (frefvec[0] != xfe_QuadRefUniform) &&
               (frefvec[1] != xfe_QuadRefVert ) && (frefvec[1] != xfe_QuadRefUniform) &&
               (frefvec[3] != xfe_QuadRefVert ) && (frefvec[3] != xfe_QuadRefUniform) &&
               (frefvec[5] != xfe_QuadRefVert ) && (frefvec[5] != xfe_QuadRefUniform))
        (*pelemref) = xfe_HexRefSliceYZ;
      else
        (*pelemref) = xfe_HexRefUniform;
      
      break;
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
//   FUNCTION Definition: xf_Face0RefEffectOnFace
static int 
xf_Face0RefEffectOnFace(enum xfe_ShapeType FShape, int fref0, int pos, 
                        int *pfref)
{
  /*
   PURPOSE:
   
   Returns refinement on subface of a hanging face, implied by the
   subface position (pos) and by the refinement of the original face
   (fref0).
   
   Note: this function relies on order of positions in xf_MeshHangStruct.h
   
   INPUTS:
   
   FShape : which face shape to consider
   fref0  : refinement of the original coarse face
   pos    : position of subface in original face
   
   OUTPUTS: 
   
   (*pfref) : required refinement on the subface
   
   RETURN:
   
   Error Code
   */
  (*pfref) = fref0;
  
  // return immediately if no refinement or not hanging
  if ((fref0 == 0) || (pos == 0)) return xf_OK;
  
  switch (FShape){
    case xfe_Segment:
      // pos != 0 means original segment is refined
      // we do not need to refine any subsegments
      (*pfref) = xfe_SegRefNone;
      break;
    case xfe_Quadrilateral:
      // first four positions mean face0 is already uniformly refined
      if (pos <= xfe_QuadPosNE) 
        (*pfref) = xfe_QuadRefNone; // nothing else to do
      else if (pos <= xfe_QuadPosRight){ // face0 is split vertically
        if ((fref0 == xfe_QuadRefHoriz) || (fref0 == xfe_QuadRefUniform))
          (*pfref) = xfe_QuadRefHoriz;
        else (*pfref) = xfe_QuadRefNone;
      }  
      else if (pos <= xfe_QuadPosTop){ // face0 is split horizontally
        if ((fref0 == xfe_QuadRefVert) || (fref0 == xfe_QuadRefUniform))
          (*pfref) = xfe_QuadRefVert;
        else (*pfref) = xfe_QuadRefNone;
      }
      else return xf_Error(xf_INPUT_ERROR);
      break;
    case xfe_Triangle:
      return xf_Error(xf_NOT_SUPPORTED);
    case xfe_Tetrahedron:
    case xfe_Hexahedron:
      return xf_Error(xf_INPUT_ERROR); // these are not face shapes
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_RotateFaceRef
static int 
xf_RotateFaceRef(enum xfe_ShapeType FShape, int facerefloc, int orient,
                 int *pfaceref)
{
  /*
   PURPOSE:
   
   Rotates face refinement option from element loc space (facerefloc)
   to the face space (faceref), using prescribed orientation number.
   
   INPUTS:
   
   FShape : face shape to consider
   facerefloc : face refinement option as seen from loc elem point of view
   orient : face orientation number (stored in grid structure)
   
   OUTPUTS: 
   
   (*pfaceref) : face refinement option in face space
   
   RETURN:
   
   Error Code
   */
  
  // return immediately if no refinement
  if (facerefloc == 0){
    (*pfaceref) = 0;
    return xf_OK;
  }
  
  switch (FShape){
    case xfe_Segment:
      (*pfaceref) = facerefloc;
      break;
    case xfe_Quadrilateral:
      (*pfaceref) = facerefloc;
      if ((orient == 1) || (orient == 3) || (orient == 4) || (orient == 6)){
        // horizontal becomes vertical and vice-versa
        if (facerefloc == xfe_QuadRefHoriz) (*pfaceref) = xfe_QuadRefVert ;
        if (facerefloc == xfe_QuadRefVert ) (*pfaceref) = xfe_QuadRefHoriz;
      }
      break;
    case xfe_Triangle:
    case xfe_Tetrahedron:
    case xfe_Hexahedron:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_InvRotateFaceRef
static int 
xf_InvRotateFaceRef(enum xfe_ShapeType FShape, int faceref, int orient,
                    int *pfacerefloc)
{
  /*
   PURPOSE:
   
   Rotates face refinement option from face space (faceref) to element
   loc space (*pfacerefloc).
   
   INPUTS:
   
   FShape : face shape to consider
   faceref : face refinement option in face space
   orient : face orientation number (stored in grid structure)
   
   OUTPUTS: 
   
   (*pfacerefloc) : face refinement option as seen from loc elem point of view
   
   
   RETURN:
   
   Error Code
   */
  
  // return immediately if no refinement
  if (faceref == 0){
    (*pfacerefloc) = 0;
    return xf_OK;
  }
  
  switch (FShape){
    case xfe_Segment:
      (*pfacerefloc) = faceref;
      break;
    case xfe_Quadrilateral:
      (*pfacerefloc) = faceref;
      if ((orient == 1) || (orient == 3) || (orient == 4) || (orient == 6)){
        // horizontal becomes vertical and vice-versa
        if (faceref == xfe_QuadRefHoriz) (*pfacerefloc) = xfe_QuadRefVert ;
        if (faceref == xfe_QuadRefVert ) (*pfacerefloc) = xfe_QuadRefHoriz;
      }
      break;
    case xfe_Triangle:
    case xfe_Tetrahedron:
    case xfe_Hexahedron:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefCompare
static int 
xf_RefCompare(enum xfe_ShapeType Shape, int ref1, int ref2, 
              enum xfe_RelType *prel)
{
  /*
   PURPOSE:
   
   Determine how two refinement options are related.  The relation of
   ref1 w.r.t ref2 is returned:
   
   ref1 == ref2:  -> rel = Equal
   ref1 in ref2:  -> rel = Subset
   ref2 in ref1:  -> rel = Superset
   
   INPUTS:
   
   Shape : which element shape to consider
   ref1, ref2 : the two refinement options to compare
   
   OUTPUTS: 
   
   (*prel) : relation between ref1 and ref2, as described above
   
   RETURN:
   
   Error Code
   */
  
  if (ref1 == ref2) {
    (*prel) = xfe_RelEqual;
    return xf_OK;
  }
  if (ref1 == 0){
    (*prel) = xfe_RelSubset;
    return xf_OK;
  }
  if (ref2 == 0){
    (*prel) = xfe_RelSuperset;
    return xf_OK;
  }
  
  switch (Shape){
    case xfe_Segment:
      return xf_Error(xf_INPUT_ERROR); // should not get here
      break;
    case xfe_Quadrilateral:
      if (ref2 == xfe_QuadRefUniform)
        (*prel) = xfe_RelSubset;
      else if (ref1 == xfe_QuadRefUniform)
        (*prel) = xfe_RelSuperset;
      else (*prel) = xfe_RelDisjoint;
      break;
    case xfe_Hexahedron:
      if (ref2 == xfe_QuadRefUniform)
        (*prel) = xfe_RelSubset;
      else if (ref1 == xfe_QuadRefUniform)
        (*prel) = xfe_RelSuperset;
      else{
        // both differently refined and neither is uniform
        if (((ref1==xfe_HexRefSliceX) && ((ref2==xfe_HexRefSliceXY) || (ref2==xfe_HexRefSliceXZ))) ||
            ((ref1==xfe_HexRefSliceY) && ((ref2==xfe_HexRefSliceXY) || (ref2==xfe_HexRefSliceYZ))) ||
            ((ref1==xfe_HexRefSliceZ) && ((ref2==xfe_HexRefSliceXZ) || (ref2==xfe_HexRefSliceYZ))))
          (*prel) = xfe_RelSubset;
        else if (((ref2==xfe_HexRefSliceX) && ((ref1==xfe_HexRefSliceXY) || (ref1==xfe_HexRefSliceXZ))) ||
                 ((ref2==xfe_HexRefSliceY) && ((ref1==xfe_HexRefSliceXY) || (ref1==xfe_HexRefSliceYZ))) ||
                 ((ref2==xfe_HexRefSliceZ) && ((ref1==xfe_HexRefSliceXZ) || (ref1==xfe_HexRefSliceYZ))))
          (*prel) = xfe_RelSuperset;
        else
          (*prel) = xfe_RelDisjoint;
      }
      break;
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
//   FUNCTION Definition: xf_LeastCommonRef
static int 
xf_LeastCommonRef(enum xfe_ShapeType Shape, int ref1, int ref2,
                  int *pref)
{
  /*
   PURPOSE:
   
   Returns the simplest shape refinement that includes both ref1
   and ref2.
   
   INPUTS:
   
   Shape : which element shape to consider
   ref1, ref2 : the two refinement options to merge
   
   OUTPUTS: 
   
   (*pref) : the least common refinement option
   
   RETURN:
   
   Error Code
   */
  int ierr;
  enum xfe_RelType rel;
  
  // if equal, return right away
  if (ref1 == ref2){
    (*pref) = ref1;
    return xf_OK;
  }
  
  // 0 is always assumed to mean no refinement
  if (ref1 == 0){
    (*pref) = ref2; 
    return xf_OK;
  }
  else if (ref2 == 0){
    (*pref) = ref1; 
    return xf_OK;
  }
  
  // consider each shape separately
  switch (Shape){
    case xfe_Segment:
      if ((ref1 != xfe_SegRefUniform) || (ref2 != xfe_SegRefUniform))
        return xf_Error(xf_INPUT_ERROR);
      (*pref) = xfe_SegRefUniform; // simple in this case
      break;
    case xfe_Quadrilateral:
      /* Both refs are nonzero and they are not equal.  Either one of
       them is uniform, or one is horiz and one is vert, in which case
       the inclusive common ref is uniform. */
      (*pref) = xfe_QuadRefUniform;
      break;
    case xfe_Hexahedron:
      if ((ref1 == xfe_HexRefUniform) || (ref2 == xfe_HexRefUniform))
        (*pref) = xfe_HexRefUniform;
      else{
        ierr = xf_Error(xf_RefCompare(Shape, ref1, ref2, &rel));
        if (ierr != xf_OK) return ierr;
        if (rel == xfe_RelSubset) (*pref) = ref2;
        else if (rel == xfe_RelSuperset) (*pref) = ref1;
        else if (rel == xfe_RelDisjoint){
          // disjoint and neither is uniform
          if ( ((ref1 == xfe_HexRefSliceX) && (ref2 == xfe_HexRefSliceY)) ||
              ((ref1 == xfe_HexRefSliceY) && (ref2 == xfe_HexRefSliceX)) )
            (*pref) = xfe_HexRefSliceXY;
          else if ( ((ref1 == xfe_HexRefSliceX) && (ref2 == xfe_HexRefSliceZ)) ||
                   ((ref1 == xfe_HexRefSliceZ) && (ref2 == xfe_HexRefSliceX)) )
            (*pref) = xfe_HexRefSliceXZ;
          else if ( ((ref1 == xfe_HexRefSliceY) && (ref2 == xfe_HexRefSliceZ)) ||
                   ((ref1 == xfe_HexRefSliceZ) && (ref2 == xfe_HexRefSliceY)) )
            (*pref) = xfe_HexRefSliceYZ;
          else
            (*pref) = xfe_HexRefUniform;
        }
        else return xf_Error(xf_CODE_LOGIC_ERROR);
      }
      break;
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
//   FUNCTION Definition: xf_Pos2FaceRef
static int 
xf_Pos2FaceRef(enum xfe_ShapeType FShape, int pos, int *pfacerefloc)
{
  /*
   PURPOSE:
   
   Returns the face refinement implied by hanging-face position, pos.
   This refinement is how the original face looks like from the coarse
   element.  Note, this function relies on some of the enumerated types
   in xfe_MeshHangStruct.h.
   
   INPUTS:
   
   FShape : shape of face
   pos : position of hanging face on coarse face
   
   OUTPUTS: 
   
   (*pfacerefloc) : implied refinement level on original face
   
   RETURN:
   
   Error Code
   */
  
  // pos = 0 means not a hanging face
  if (pos == 0){
    (*pfacerefloc) = 0;
    return xf_OK;
  }
  
  switch (FShape){
    case xfe_Segment:
      (*pfacerefloc) = xfe_SegRefUniform; // only option for a nonzero pos
      break;
    case xfe_Quadrilateral:
      // first four positions are for uniform refinement
      if (pos <= xfe_QuadPosNE) (*pfacerefloc) = xfe_QuadRefUniform;
      // next comes the vertical split (two positions)
      else if (pos <= xfe_QuadPosRight) (*pfacerefloc) = xfe_QuadRefVert;
      // next comes the horizontal split (two positions)
      else if (pos <= xfe_QuadPosTop) (*pfacerefloc) = xfe_QuadRefHoriz;
      else return xf_Error(xf_INPUT_ERROR);
      break;
    case xfe_Triangle:
      return xf_Error(xf_NOT_SUPPORTED);
    case xfe_Tetrahedron:
    case xfe_Hexahedron:
      return xf_Error(xf_INPUT_ERROR); // these are not face shapes
      break;
    default:
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_EnsureRefConsistency
static int 
xf_EnsureRefConsistency(xf_All *All, xf_Vector *RefIndicator, int *IFaceRef,
                        xf_Vector *ElemFaceRef)
{
  /*
   PURPOSE:
   
   Ensures that desired element refinements in RefIndicator are
   consistent across faces and that no "second degree" hanging faces
   are created.  Does this by iteratively looping over elements, and
   adjusting requested refinement based on the implied face refinements.
   
   INPUTS:
   
   All : All structure
   RefIndicator : element-based refinement indicator (stores initial request)
   
   OUTPUTS: 
   
   IFaceRef : will store refinement indicator on each face (in face ref space)
   ElemFaceRef : will store refinement indicator for each element original face
   from the point of view of the local reference elem
   
   RETURN:
   
   Error Code
   */
  
  int ierr, k, iter, negrp, egrp, elem, face;
  int hang, face0, pos, iiface, faceorig;
  int ref, ref0, frefloc, fref, freforig, orient;
  int *frefvec, nfacemax, nface0;
  enum xfe_ShapeType Shape, FShape;
  enum xfe_Bool changed;
  enum xfe_RelType rel;
  int maxrefiter = 1000; // maximum number of iters to ensure consistency
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  // nothing to do if in 1D
  if (Mesh->Dim == 1) return xf_OK;
  
  frefvec = NULL;
  nfacemax = -1;
  
  changed = xfe_True;
  iter = 0;
  
  while (changed){
    
    changed = xfe_False;
    
    /* Loop over elements and use RefIndicator to update IFaceRef.
     Note, IFaceRef stores the required refinement of each face
     (fine face if hanging). */
    for (egrp=0; egrp<negrp; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
        ref = RefIndicator->GenArray[egrp].iValue[elem][0];
        if (ref != 0){
          // get element Shape
          ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
          if (ierr != xf_OK) return ierr;
          
          // loop over faces
          for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
            
            // skip non-interior faces
            if (Mesh->ElemGroup[egrp].Face[elem][face].Group != xf_INTERIORFACE) continue;
            
            // get interior face number
            iiface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
            
            // get original face (while checking if hanging face)
            ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, &face0, &pos, NULL));
            if (ierr != xf_OK) return ierr;
            
            // get effect of ref on face0
            ierr = xf_Error(xf_ElemRefEffectOnFace(Shape, ref, face0, &frefloc));
            if (ierr != xf_OK) return ierr;
            
            // if ref does not affect face, continue
            if (frefloc <= 0) continue;
            
            // get Shape of face = Shape of face0
            ierr = xf_Error(xf_FaceShape(Shape, face0, &FShape));
            if (ierr != xf_OK) return ierr;
            
            // frefloc is now the refinement on face0
            // we want refinement on face (different if on coarse hang side)
            if (hang != 0){ // on coarse side of hang face
              // (frefloc, pos) -> new frefloc
              ierr = xf_Error(xf_Face0RefEffectOnFace(FShape, frefloc, pos, &frefloc));
              if (ierr != xf_OK) return ierr;
            }
            
            // pull off face orientation
            ierr = xf_Error(xf_GetFaceOrient(Mesh, egrp, elem, face, &orient));
            if (ierr != xf_OK) return ierr;
	    	    
            // transform effect to iface loc space via orient
            ierr = xf_Error(xf_RotateFaceRef(FShape, frefloc, orient, &fref));
            if (ierr != xf_OK) return ierr;
            
            freforig = IFaceRef[iiface];
            
            // update IFaceRef with least common superset refinement
            ierr = xf_Error(xf_LeastCommonRef(FShape, freforig, fref, IFaceRef + iiface));
            if (ierr != xf_OK) return ierr;
            
            changed = (changed || (freforig != IFaceRef[iiface]));
            
            //  - if on coarse side of a hang face, use ref, pos to
            //    determine how iface needs to be refined (if at all)
            //    don't forget to use orient to convert to iface loc space
            //    (note, in 2D, this step can be skipped, since a segment
            //     has only one refinement mode ... uniform)
            //
            //  3/24/09: taking care of this with Face0RefEffectOnFace call above
            //
            /* if ((Shape != xfe_Triangle) && (Shape != xfe_Quadrilateral)) */
            /* 	      return xf_Error(xf_NOT_SUPPORTED); */
            
            
          } // face
        }
      } // elem
    } // egrp
    
    
    // *** NEED A PARALLEL COMMUNICATION of IFaceRef HERE ... not if doing this in serial! ***
    
    
    /* Loop over elements and use IFaceRef to set ElemFaceRef.  Check
     for consistency with RefIndicator: if necessary, modify
     RefIndicator.  In particular, flag new elements here. */
    for (egrp=0; egrp<negrp; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
        
        // get element Shape
        ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
        if (ierr != xf_OK) return ierr;
        
        // number of original faces in Shape
        ierr = xf_Error(xf_Shape2nFace(Shape, &nface0));
        if (ierr != xf_OK) return ierr;
        
        if (nface0 > nfacemax){
          nfacemax = nface0;
          ierr = xf_Error(xf_ReAlloc( (void **) &frefvec, nfacemax, sizeof(int)));
          if (ierr != xf_OK) return ierr;
        }
        
        // original refinement request
        ref0 = RefIndicator->GenArray[egrp].iValue[elem][0];
        
        // initialize frefvec to zero = no refinement
        for (k=0; k<nface0; k++) frefvec[k] = 0; 
        
        // loop over faces and fill in frefvec
        for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
          
          // skip non-interior faces
          if (Mesh->ElemGroup[egrp].Face[elem][face].Group != xf_INTERIORFACE) continue;
          
          // get interior face number
          iiface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
          
          // continue if no refinement on iiface
          if (IFaceRef[iiface] <= 0) continue;
          
          // get original face (while checking if hanging face)
          ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, &face0, &pos, NULL));
          if (ierr != xf_OK) return ierr;
          
          // get Shape of face0
          ierr = xf_Error(xf_FaceShape(Shape, face0, &FShape));
          if (ierr != xf_OK) return ierr;
          
          if (hang){   
            /* Here, we are on coarse side of a hanging face, and the
             subface wants a refinement.  This means face0 needs to
             be refined again to maintain consistency.  
             
             Just refine according to hang pos.  Ignore comment below.
             
             [The safest option is to set the face0 refinement to
             uniform; however, if the current coarse-side refinement
             of face0 is the same as the requested subface
             refinement, we can use that instead (and this does not
             have to be uniform). So use a LeastCommonRef
             instead.] */
            
            // pos implies a coarse-side refinement level of face0
            ierr = xf_Error(xf_Pos2FaceRef(FShape, pos, &frefloc));
            if (ierr != xf_OK) return ierr;
            
            if (frefloc == 0) return xf_Error(xf_CODE_LOGIC_ERROR); // otherwise not a hang!
            
            frefvec[face0] = frefloc;
          }
          else{
            // not a hanging face; check if elem implies a refinement on face0
            ierr = xf_Error(xf_ElemRefEffectOnFace(Shape, ref0, face0, &frefloc));
            if (ierr != xf_OK) return ierr;
            
            if (frefloc != 0){
              /* if so, we are in a position where elem wants face0
               refined -- however, IFaceRef[iiface] might indicate
               more refinement.  For example, in 3D, two hexes
               splitting a face, each in opposite directions, means
               that uniform refinement of the common face is required.
               So we need to use IFaceRef[iiface] here. */
              
              // pull off face orientation
              ierr = xf_Error(xf_GetFaceOrient(Mesh, egrp, elem, face, &orient));
              if (ierr != xf_OK) return ierr;
              
              // on a normal face, inverse rotate the face ref into elem loc coords
              ierr = xf_Error(xf_InvRotateFaceRef(FShape, IFaceRef[iiface], orient, frefvec+face0));
              if (ierr != xf_OK) return ierr;
            }
            // no else here: we allow for one level of non-conformity
          }
        } // face
        
        // Determine new ref based on frefvec
        ierr = xf_Error(xf_FaceRefEffectOnElem(Shape, frefvec, &ref));
        if (ierr != xf_OK) return ierr;
        
        
        // Take the least common refinement option with the original request
        ierr = xf_Error(xf_LeastCommonRef(Shape, ref0, ref, &ref));
        if (ierr != xf_OK) return ierr;
        
        // make sure refinement request did not get coarser
        ierr = xf_Error(xf_RefCompare(Shape, ref0, ref, &rel));
        if (ierr != xf_OK) return ierr;
        
        if ((rel != xfe_RelEqual) && (rel != xfe_RelSubset)){
          xf_printf("Error: Refinement request got coarser when ensuring consistency.\n");
          return xf_Error(xf_CODE_LOGIC_ERROR);
        }
        
        // did the element refinement change?  If so, update it.
        changed = (changed || (ref != ref0));
        RefIndicator->GenArray[egrp].iValue[elem][0] = ref;
        
      } // elem
    } // egrp
    
    iter++;
    if (iter > maxrefiter){
      xf_printf("Could not satisfy hanging-node refinement consistency in a \n");
      xf_printf("reasonable number of iterations.  Input grid may not be valid.\n");
      return xf_Error(xf_OUT_OF_BOUNDS);
    }
  }
  
  xf_Release( (void *) frefvec);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_FillRefPriority
static int 
xf_FillRefPriority(xf_All *All, xf_Vector *RefIndicator, int *IFaceRef,
                   xf_Vector *RefPriority, int *pnPriority)
{
  /*
   PURPOSE:
   
   Assigns to each element a refinement priority that is used to order
   the element refinements.  This priority is necessary to maintain at
   most one degree of non-conformity at all times throughout the
   refinement.
   
   INPUTS:
   
   All : All structure
   RefIndicator : element-based refinement indicator (assumed consistent)
   IFaceRef : associated face refinement vector
   
   OUTPUTS: 
   
   RefPriority : will contain priority number for each element (higher number
   elements should get refined first)
   (*pnPriority) : number of priority levels
   
   RETURN:
   
   Error Code
   */
  int ierr, iter, nIFaceRegular, iiface, changed;
  int refL, refR;
  enum xfe_Bool MeshIsParallel;
  int maxpriorityiter = 1000; // maximum number of priority iterations
  xf_IFace *IFace;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  
  // Regular number of ifaces (for parallel)
  nIFaceRegular = -1;
  if (MeshIsParallel = (Mesh->ParallelInfo != NULL))
    nIFaceRegular = Mesh->ParallelInfo->nIFaceRegular;
  
  // initialize priorities to zero
  ierr = xf_Error(xf_SetZeroVector(RefPriority)); 
  if (ierr != xf_OK) return ierr;
  
  changed = 1;
  iter = 0;
  
  while (changed){
    
    changed = 0;
    
    // Loop over interior faces (backwards, to affect halo first)
    for (iiface=Mesh->nIFace-1; iiface>=0; iiface--){
      
      if (IFaceRef[iiface] != 0){ // iiface is marked for refinement
        IFace = Mesh->IFace + iiface;
        
        if (ah_DEBUG)
          xf_printf("Setting priority,iiface=%d (%d %d, %d %d): curpL=%d, curpR=%d, hang=%d\n",
                    iiface, IFace->ElemGroupL, IFace->ElemL,IFace->ElemGroupR,  IFace->ElemR,
                    RefPriority->GenArray[IFace->ElemGroupL].iValue[IFace->ElemL][0],
                    RefPriority->GenArray[IFace->ElemGroupR].iValue[IFace->ElemR][0],
                    IFace->HangNumber);
        
        // on iter, only consider iiface between two elements of priority iter
        if ((RefPriority->GenArray[IFace->ElemGroupL].iValue[IFace->ElemL][0] < iter) ||
            (RefPriority->GenArray[IFace->ElemGroupR].iValue[IFace->ElemR][0] < iter))
          continue;
        
        refL = RefIndicator->GenArray[IFace->ElemGroupL].iValue[IFace->ElemL][0];
        refR = RefIndicator->GenArray[IFace->ElemGroupR].iValue[IFace->ElemR][0];
        
        if (IFace->HangNumber > 0){
          // R element should be marked for refinement, and has priority
          if (refR==0)
            return xf_Error(xf_CODE_LOGIC_ERROR); // should have marked this elem before
          if (ah_DEBUG)
            xf_printf("setting priority on (%d,%d) to %d\n", IFace->ElemGroupR, IFace->ElemR, iter+1);
          RefPriority->GenArray[IFace->ElemGroupR].iValue[IFace->ElemR][0] = iter+1;
          changed = 1;
        }
        else if (IFace->HangNumber < 0){
          // L element should be marked for refinement, and has priority
          if (refL==0)
            return xf_Error(xf_CODE_LOGIC_ERROR); // should have marked this elem before
          if (ah_DEBUG)
            xf_printf("setting priority on (%d,%d) to %d\n", IFace->ElemGroupL, IFace->ElemL, iter+1);
          RefPriority->GenArray[IFace->ElemGroupL].iValue[IFace->ElemL][0] = iter+1;
          changed = 1;
        }
      }
      
      if (MeshIsParallel && (iiface == nIFaceRegular)){
        // begin halo exchange of RefPriority
        ierr = xf_Error(xf_HaloExchangeVectorBegin(RefPriority));
        if (ierr != xf_OK) return ierr;
      }
    } // iiface
    
    // end halo exchange of RefPriority
    if (MeshIsParallel){
      ierr = xf_Error(xf_HaloExchangeVectorEnd(RefPriority));
      if (ierr != xf_OK) return ierr;
    }
    
    // reduce-max changed over all processors
    if (MeshIsParallel){
      ierr = xf_Error(xf_MPI_Allreduce(&changed, 1, xfe_SizeInt, xfe_MPI_MAX));
      if (ierr != xf_OK) return ierr;
    }
    
    if (changed == 1){ // increment iteration number
      iter++;
      if (iter > maxpriorityiter){
        xf_printf("Exceeded max number of iterations when filling priorities.\n");
        return xf_Error(xf_OUT_OF_BOUNDS);
      }
    }
    
  } // end while changed
  
  // number of priority levels is the number of iterations taken
  (*pnPriority) = iter;
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_Ref2nElem
int 
xf_Ref2nElem(enum xfe_ShapeType Shape, int ref, int *pnelem)
{ 
  // return immediately if no refinement
  if (ref == 0){
    (*pnelem) = 1;
    return xf_OK;
  }
  
  switch (Shape){
    case xfe_Point:
      (*pnelem) = 0;
      break;
    case xfe_Segment:
      (*pnelem) = 2;
      break;
    case xfe_Quadrilateral:
      if (ref == xfe_QuadRefUniform) (*pnelem) = 4;
      else (*pnelem) = 2;
      break;
    case xfe_Hexahedron:
      if (ref == xfe_HexRefUniform) 
        (*pnelem) = 8;
      else if ((ref == xfe_HexRefSliceXY) ||
               (ref == xfe_HexRefSliceXZ) ||
               (ref == xfe_HexRefSliceYZ))
        (*pnelem) = 4;
      else
        (*pnelem) = 2;
      break;
    case xfe_Triangle:
    case xfe_Tetrahedron:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    default:
      xf_printf("Shape = %d\n", Shape);
      return xf_Error(xf_UNKNOWN_SHAPE);
      break;
  }
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_Ref2nIFace
static int 
xf_Ref2nIFace(enum xfe_ShapeType Shape, int ref, int *pniface)
{
  /*
   PURPOSE:
   
   Returns number of element interior faces that result from a
   refinement of Shape.
   
   INPUTS:
   
   Shape : shape of element in question
   ref  : refinement indicator for this element
   
   OUTPUTS: 
   
   (*pniface) : number of element interior faces
   
   RETURN:
   
   Error Code
   */
  
  
  // return immediately if no refinement
  if (ref == 0){
    (*pniface) = 0;
    return xf_OK;
  }
  
  switch (Shape){
    case xfe_Segment:
      (*pniface) = 1;
      break;
    case xfe_Quadrilateral:
      if (ref == xfe_QuadRefUniform) (*pniface) = 4;
      else (*pniface) = 1;
      break;
    case xfe_Hexahedron:
      if (ref == xfe_HexRefUniform) 
        (*pniface) = 12;
      else if ((ref == xfe_HexRefSliceXY) ||
               (ref == xfe_HexRefSliceXZ) ||
               (ref == xfe_HexRefSliceYZ))
        (*pniface) = 4;
      else
        (*pniface) = 1;
      break;
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
//   FUNCTION Definition: xf_FlipPosQuad
static int 
xf_FlipPosQuad(int pos, int *flippos) 
{
  /*
   PURPOSE:
   
   Flips a position number on a quadrilateral (imagine flipping a piece
   of paper and looking at the positions from behind).  Flip is done
   about diagonal going through lower-left corner.
   
   INPUTS:
   
   pos : position as seen from the front
   
   OUTPUTS: 
   
   (*flippos) : flipped position (as seen from the back)
   
   RETURN:
   
   Error Code
   */
  
  (*flippos) = pos;
  
  switch (pos){
    case xfe_QuadPosNone:
    case xfe_QuadPosSW:
    case xfe_QuadPosNE:
      (*flippos) = pos;
      break;
    case xfe_QuadPosSE:
      (*flippos) = xfe_QuadPosNW;
      break;
    case xfe_QuadPosNW:
      (*flippos) = xfe_QuadPosSE;
      break;
    case xfe_QuadPosLeft:
      (*flippos) = xfe_QuadPosBottom;
      break;
    case xfe_QuadPosRight:
      (*flippos) = xfe_QuadPosTop;
      break;
    case xfe_QuadPosBottom:
      (*flippos) = xfe_QuadPosLeft;
      break;
    case xfe_QuadPosTop:
      (*flippos) = xfe_QuadPosRight;
      break;
    default:
      return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CyclePosQuad
static int 
xf_CyclePosQuad(int pos, int *cyclepos) 
{
  /*
   PURPOSE:
   
   Cycles a position number counterclockwise (imagine rotating a piece
   of paper counterclockwise)
   
   INPUTS:
   
   pos : position as seen originally
   
   OUTPUTS: 
   
   (*cyclepos) : cycled position (as seen after paper turned, you stay
   stationary)
   
   RETURN:
   
   Error Code
   */
  
  (*cyclepos) = pos;
  
  switch (pos){
    case xfe_QuadPosNone:
      (*cyclepos) = pos;
      break;
    case xfe_QuadPosSW:
      (*cyclepos) = xfe_QuadPosSE;
      break;
    case xfe_QuadPosSE:
      (*cyclepos) = xfe_QuadPosNE;
      break;
    case xfe_QuadPosNE:
      (*cyclepos) = xfe_QuadPosNW;
      break;
    case xfe_QuadPosNW:
      (*cyclepos) = xfe_QuadPosSW;
      break;
    case xfe_QuadPosLeft:
      (*cyclepos) = xfe_QuadPosBottom;
      break;
    case xfe_QuadPosRight:
      (*cyclepos) = xfe_QuadPosTop;
      break;
    case xfe_QuadPosBottom:
      (*cyclepos) = xfe_QuadPosRight;
      break;
    case xfe_QuadPosTop:
      (*cyclepos) = xfe_QuadPosLeft;
      break;
    default:
      return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ConvertPos
static int 
xf_ConvertPos( enum xfe_ShapeType FShape, int pos1, int orient1, 
              int orient2, int *pos2)
{
  /*
   PURPOSE:
   
   Converts a position number as seen from the fine side of a hanging
   face to the position number as seen from the coarse side of a
   hanging face.  This function relies on the order of some of the
   enumerated types in xf_MeshHangStruct.h.
   
   INPUTS:
   
   FShape : face shape type
   pos1 : position as seen from the fine side
   orient1 : orientation of the ref face on the fine side
   orient2 : orientation of the ref face on the coarse side
   
   OUTPUTS: 
   
   (*pos2) : position as seen from coarse side of hanging face
   
   RETURN:
   
   Error Code
   */
  int ierr, i;
  int orientrel;
  int cycle, flip;
  
  // return immediately if trivial
  if ((pos1 == 0) || (orient1 == orient2)){
    (*pos2) = pos1;
    return xf_OK;
  } 
  
  switch (FShape){
    case xfe_Segment:
      (*pos2) = xfe_SegPosLast-pos1;
      break;
    case xfe_Quadrilateral:
      // obtain orientation of fine side relative to the coarse side
      ierr = xf_Error(xf_GetRelativeOrient(FShape, orient1, orient2, &orientrel));
      if (ierr != xf_OK) return ierr;
      
      if (ah_DEBUG)
        xf_printf("orient1=%d, orient2=%d, orientrel = %d\n", orient1, orient2, orientrel);
      
      cycle = orientrel%4;
      flip  = orientrel/4;
      
      // flip, then cycle, the position
      (*pos2) = pos1;
      if (flip){
        ierr = xf_Error(xf_FlipPosQuad((*pos2), pos2));
        if (ierr != xf_OK) return ierr;
      }
      for (i=0; i<cycle; i++){
        ierr = xf_Error(xf_CyclePosQuad((*pos2), pos2));
        if (ierr != xf_OK) return ierr;
      }
      
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_PosCompare
static int 
xf_PosCompare(enum xfe_ShapeType Shape, int pos1, int pos2, 
              enum xfe_RelType *prel, enum xfe_Bool *pcompatible)
{
  /*
   PURPOSE:
   
   Determine how two refinement positions are related.  The relation of
   pos1 w.r.t pos2 is returned:
   
   pos1 == pos2:       -> rel = Equal
   pos1 contains pos2: -> rel = Superset
   pos2 contains pos1: -> rel = Subset
   
   Relies on order of enumerated types in xf_MeshHangStruct.h
   
   INPUTS:
   
   Shape : which element shape to consider
   pos1, pos2 : the two positions to compare
   
   OUTPUTS: 
   
   (*prel) : relation between pos1 and pos2, as described above
   (*pcompatible) : Even if disjoint, pos1 and pos2 can be compatible
   (e.g. both arising from a uniform refinement).
   This flag indicates if pos1 and pos2 are indeed
   compatible.
   
   RETURN:
   
   Error Code
   */
  int posa, posb, k;
  
  (*pcompatible) = xfe_True;
  
  if (pos1 == pos2) {
    (*prel) = xfe_RelEqual;
    return xf_OK;
  }
  if (pos1 == 0){
    (*prel) = xfe_RelSuperset;
    return xf_OK;
  }
  if (pos2 == 0){
    (*prel) = xfe_RelSubset;
    return xf_OK;
  }
  
  switch (Shape){
    case xfe_Segment:
      (*prel) = xfe_RelDisjoint;
      break;
    case xfe_Quadrilateral:
      // incompatible if pos1 is left/right and pos2 is up/down or vice-versa
      if (  (pos1 >= xfe_QuadPosLeft ) && (pos2 >= xfe_QuadPosLeft )  && 
          ((pos1 <= xfe_QuadPosRight) || (pos2 <= xfe_QuadPosRight)) &&
          ((pos1 >  xfe_QuadPosRight) || (pos2 >  xfe_QuadPosRight)) ){
        (*pcompatible) = xfe_False;
      }
      else{
        // sort in increasing order: posa < posb
        posa = pos1; posb = pos2;
        if (posa > posb) swap(posa, posb, k);
        
        if ( ( (posb==xfe_QuadPosBottom) && ((posa==xfe_QuadPosSW) || (posa==xfe_QuadPosSE)) ) ||
            ( (posb==xfe_QuadPosTop   ) && ((posa==xfe_QuadPosNW) || (posa==xfe_QuadPosNE)) ) ||
            ( (posb==xfe_QuadPosLeft  ) && ((posa==xfe_QuadPosSW) || (posa==xfe_QuadPosNW)) ) ||
            ( (posb==xfe_QuadPosRight ) && ((posa==xfe_QuadPosSE) || (posa==xfe_QuadPosNE)) ) )
          (*prel) = xfe_RelSubset;
        else
          (*prel) = xfe_RelDisjoint;
        
        if (posa != pos1){ // we swapped the order above, so now swap the relation
          if ((*prel) == xfe_RelSubset  ) (*prel) = xfe_RelSuperset;
          else if ((*prel) == xfe_RelSuperset) (*prel) = xfe_RelSubset;
        }
      }
      break;
    case xfe_Hexahedron:
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
//   FUNCTION Definition: xf_RelativePos
static int
xf_RelativePos( enum xfe_ShapeType FShape, int pos1, int pos2, 
               int *relpos)
{
  /*
   PURPOSE:
   
   Computes the relative position of pos2 within pos1.  pos1 must be a
   superset of pos2, otherwise an error occurs.  Should not be called
   for trivial cases (e.g. no or equal refinement).
   
   INPUTS:
   
   FShape : face shape type
   pos1 : first position
   pos2 : second position
   
   OUTPUTS: 
   
   (*relpos) : relative position of pos2 within pos1
   
   RETURN:
   
   Error Code
   */
  int ierr;
  enum xfe_RelType rel;
  enum xfe_Bool compatible;
  
  // should not be called for trivial cases
  if ((pos1 == 0) || (pos2 == 0) || (pos1 == pos2))
    return xf_Error(xf_INPUT_ERROR);
  
  // check that pos2 is actually within pos1
  ierr = xf_Error(xf_PosCompare(FShape, pos1, pos2, &rel, &compatible));
  if (ierr != xf_OK) return ierr;
  
  if (rel != xfe_RelSuperset){
    xf_printf("pos1 = %d is not superset of pos2 = %d\n", pos1, pos2);
    return xf_Error(xf_INPUT_ERROR);
  }
  
  
  switch (FShape){
    case xfe_Segment:
      // should not be called for segments
      return xf_Error(xf_INPUT_ERROR);
      break;
    case xfe_Quadrilateral:
      
      if (pos1==xfe_QuadPosBottom)
        (*relpos) = ((pos2==xfe_QuadPosSW) ? xfe_QuadPosLeft : xfe_QuadPosRight);
      else if (pos1==xfe_QuadPosTop   )
        (*relpos) = ((pos2==xfe_QuadPosNW) ? xfe_QuadPosLeft : xfe_QuadPosRight);
      else if (pos1==xfe_QuadPosLeft  )
        (*relpos) = ((pos2==xfe_QuadPosSW) ? xfe_QuadPosBottom : xfe_QuadPosTop);
      else if (pos1==xfe_QuadPosRight )
        (*relpos) = ((pos2==xfe_QuadPosSE) ? xfe_QuadPosBottom : xfe_QuadPosTop);
      else return xf_Error(xf_INPUT_ERROR);
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetCommonNodeFlag
static int
xf_GetCommonNodeFlag( enum xfe_ShapeType FShape, int pos, int *CommonNode)
{
  /*
   PURPOSE:
   
   Returns a CommonNode vector that indicates, via a nonzero integer
   value, nodes on FShape that are common between the shape and a
   refined child in position pos.  The node ordering is assumed to be
   as returned by Q1NodesOnFace.
   
   INPUTS:
   
   FShape : face shape type
   pos : position of child
   
   OUTPUTS: 
   
   CommonNode : integer vector over face q1 nodes, containing nonzero values
   at nodes common to parent and child 
   
   RETURN:
   
   Error Code
   */
  int ierr, i;
  enum xfe_RelType rel;
  enum xfe_Bool compatible;
  
  // should not be called for a case
  if (pos == 0) return xf_Error(xf_INPUT_ERROR);
  
  switch (FShape){
    case xfe_Segment:
      CommonNode[0] = CommonNode[1] = 0;
      if (pos==xfe_SegPosLeft)
        CommonNode[0] = 1;
      else
        CommonNode[1] = 1;
      break;
    case xfe_Quadrilateral:
      for (i=0; i<4; i++) CommonNode[i] = 0;
      switch (pos){
        case xfe_QuadPosSW: CommonNode[0] = 1; break;
        case xfe_QuadPosSE: CommonNode[1] = 1; break;
        case xfe_QuadPosNE: CommonNode[2] = 1; break;
        case xfe_QuadPosNW: CommonNode[3] = 1; break;
        case xfe_QuadPosLeft:
          CommonNode[0] = 1;
          CommonNode[3] = 1;
          break;
        case xfe_QuadPosRight:
          CommonNode[1] = 1;
          CommonNode[2] = 1;
          break;
        case xfe_QuadPosBottom:
          CommonNode[0] = 1;
          CommonNode[1] = 1;
          break;
        case xfe_QuadPosTop:
          CommonNode[2] = 1;
          CommonNode[3] = 1;
          break;
        default:
          return xf_Error(xf_INPUT_ERROR);
          break;
      }
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_HangingFaceAdd
static int 
xf_HangingFaceAdd(xf_Mesh *Mesh, xf_IFace IFace, int FineIsL, int egF, 
                  int eF, int faceF, int orientF, int posF, int *piifacenew)
{
  /*
   PURPOSE:
   
   Adds a hanging face to Mesh.  The hanging face is based on iiface.
   (egF, eF) is the fine element for the new hanging face.  The coarse
   element is on the other side of iiface, according to FineIsL.  On
   the coarse element, the hang number is set appropriately using posF
   and face0C.
   
   INPUTS:
   
   Mesh : Mesh structure
   IFace : original interior face
   FineIsL : true is fine element is on the left of iiface
   egF, eF : fine element 
   faceF : local face number on fine element
   orientF : orientation of face relative to fine element
   posF : position of fine element relative to its parent
   
   OUTPUTS: 
   
   (*piifacenew) : interior face number of new face
   Also,  new face is added, or iiface is made hanging if this is the
   first hanging face being added.
   
   RETURN:
   
   Error Code
   */
  int ierr, i, nn, fvec0[8];
  int egC, eC, faceC, face0C, orientC, nface0C;
  int posC, hangC, hang, iiface, iifacenew;
  enum xfe_ShapeType ShapeC, FShapeC;
  
  
  if (FineIsL){
    egC     = IFace.ElemGroupR;
    eC      = IFace.ElemR;
    face0C  = IFace.FaceR; // local face number on the coarse element
    orientC = IFace.OrientR;
  }
  else{
    egC     = IFace.ElemGroupL;
    eC      = IFace.ElemL;
    face0C  = IFace.FaceL;
    orientC = IFace.OrientL;
  }
  
  
  if (ah_DEBUG){
    xf_printf("In HangingFaceAdd.\n");
    xf_printf("egC = %d, eC = %d, eF = %d, face0C = %d\n", egC, eC, eF, face0C);
    xf_printf("orientF = %d, orientC =%d\n", orientF, orientC);
    if (eC == 13)
      xf_printf("TEST: %d\n", Mesh->ElemGroup[0].Face[13][1].Group);
    // get Q1 nodes on face0C
    ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egC].QBasis, 
                                     Mesh->ElemGroup[egC].QOrder, 
                                     face0C, &nn, fvec0));
    if (ierr != xf_OK) return ierr;
    // print out global node numbers
    for (i=0; i<nn; i++) 
      xf_printf("Coarse node %d = %d\n", fvec0[i], Mesh->ElemGroup[egC].Node[eC][fvec0[i]]);
  }
  
  // obtain ShapeC
  ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egC].QBasis, &ShapeC));
  if (ierr != xf_OK) return ierr;
  
  
  // number of original faces in ShapeC
  ierr = xf_Error(xf_Shape2nFace(ShapeC, &nface0C));
  if (ierr != xf_OK) return ierr;
  
  // get Shape of face0C
  ierr = xf_Error(xf_FaceShape(ShapeC, face0C, &FShapeC));
  if (ierr != xf_OK) return ierr;
  
  // convert posF to posC
  ierr = xf_Error(xf_ConvertPos(FShapeC, posF, orientF, orientC, &posC));
  if (ierr != xf_OK) return ierr;
  if (ah_DEBUG) xf_printf("posF = %d, posC =%d\n", posF, posC);
  
  // convert posC to hangC
  ierr = xf_Error(xf_Pos2Hang(FShapeC, nface0C, face0C, posC, &hangC));
  if (ierr != xf_OK) return ierr;
  
  
  
  // is this the first hanging face we're putting on face0C?
  if (Mesh->ElemGroup[egC].Face[eC][face0C].Group != xf_INTERIORFACE) 
    return xf_Error(xf_INPUT_ERROR); // must be on an interior face
  iiface = Mesh->ElemGroup[egC].Face[eC][face0C].Number;
  hang = Mesh->IFace[iiface].HangNumber;
  
  if ((hang == 0) || 
      ((!FineIsL) && (hang > 0)) || 
      (( FineIsL) && (hang < 0))){ // first time: overwrite iiface
    iifacenew = iiface;
    faceC = face0C;
  }
  else{ // hanging face already exists, add a new one
    iifacenew = Mesh->nIFace++;
    faceC = Mesh->ElemGroup[egC].nFace[eC]++;
  }
  
  if (ah_DEBUG){
    xf_printf("egF = %d, egC = %d, eF = %d, eC = %d, faceF = %d, faceC = %d, face0C = %d\n",
              egF, egC, eF, eC, faceF, faceC, face0C);  fflush(stdout);
    xf_printf("Mesh->ElemGroup[egC].nFace[eC] = %d\n", Mesh->ElemGroup[egC].nFace[eC]);
  }
  
  
  Mesh->ElemGroup[egF].Face[eF][faceF].Group  = xf_INTERIORFACE;
  Mesh->ElemGroup[egF].Face[eF][faceF].Number = iifacenew;
  Mesh->ElemGroup[egC].Face[eC][faceC].Group  = xf_INTERIORFACE;
  Mesh->ElemGroup[egC].Face[eC][faceC].Number = iifacenew;
  
  // note, orientations are consistent placeholders only
  if (FineIsL){
    Mesh->IFace[iifacenew].ElemGroupL = egF;
    Mesh->IFace[iifacenew].ElemL      = eF;
    Mesh->IFace[iifacenew].FaceL      = faceF;
    Mesh->IFace[iifacenew].OrientL    = orientF;
    Mesh->IFace[iifacenew].ElemGroupR = egC;
    Mesh->IFace[iifacenew].ElemR      = eC;
    Mesh->IFace[iifacenew].FaceR      = faceC;
    Mesh->IFace[iifacenew].OrientR    = orientC;
    Mesh->IFace[iifacenew].HangNumber = hangC;
  }
  else{
    Mesh->IFace[iifacenew].ElemGroupR = egF;
    Mesh->IFace[iifacenew].ElemR      = eF;
    Mesh->IFace[iifacenew].FaceR      = faceF;
    Mesh->IFace[iifacenew].OrientR    = orientF;
    Mesh->IFace[iifacenew].ElemGroupL = egC;
    Mesh->IFace[iifacenew].ElemL      = eC;
    Mesh->IFace[iifacenew].FaceL      = faceC;
    Mesh->IFace[iifacenew].OrientL    = orientC;
    Mesh->IFace[iifacenew].HangNumber = -hangC;
  }
  
  if (ah_DEBUG) xf_printf("iifacenew = %d\n", iifacenew);  fflush(stdout);
  (*piifacenew) = iifacenew;
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetIFaceValues
static void 
xf_SetIFaceValues(xf_IFace *IFace, int *val)
{
  /*
   PURPOSE:
   
   Sets interior face values from an integer vector
   
   INPUTS:
   
   val : values (6 of them) [ElemL, FaceL, OrientL, ElemR, FaceR, OrientR] 
   
   OUTPUTS: 
   
   IFace : pointer to interior face that is set
   
   RETURN: None
   */
  IFace->ElemL   = val[0]; 
  IFace->FaceL   = val[1]; 
  IFace->OrientL = val[2];
  IFace->ElemR   = val[3]; 
  IFace->FaceR   = val[4];
  IFace->OrientR = val[5];
}


/******************************************************************/
//   FUNCTION Definition: xf_RefineGenericElem_Seg
static int 
xf_RefineGenericElem_Seg(int ref, int *pnNodeQ1, real *coord, 
                         int *pnvert, int *vert2corner,
                         int *pnelem, int *elem2vert, int *elem2pos,
                         int *pniface, xf_IFace *ifacelist, int *nbface, 
                         xf_BFace *bfacelist, int *bface2pos)
{
  /*
   PURPOSE:
   
   Performs a generic refinement of a Seg according to ref.
   
   INPUTS:
   
   ref  : refinement indicator
   
   OUTPUTS: 
   
   (*pnNodeQ1) : number of Q1 nodes for Shape (and hence for subelements)
   coord : unrolled array of vertex coordinates
   (*pnvert) : number of vertices
   vert2corner : mapping from vertices to corners of original element
   (*pnelem) : number of subelements
   elem2vert: unrolled array of Q1 nodes for each element
   elem2pos : position of each new elem w.r.t original elem
   (*pniface) : number of interior faces
   ifacelist : list of element-interior faces
   (*nbface) : number of element-boundary faces
   bfacelist : list of element-boundary faces
   bface2pos : index for each bface identifying position w.r.t. original face
   
   RETURN:
   
   Error Code
   */
  int k;
  
  int eUniform[4] = {0,1, 1,2};
  int epUniform[2] = {xfe_SegPosLeft, xfe_SegPosRight};
  int vcUniform[3] = {0,-1,1}; // vert2corner
  int beUniform[2] = {0,1}; // bface.elem
  int bfUniform[2] = {0,1}; // bface.face
  int bpUniform[2] = {xfe_SegPosNone, xfe_SegPosNone};
  int efoUniform[1][6] = {{0,1,0,1,0,0}};
  real xUniform[3] = {0., .5, 1.};
  
  (*pnNodeQ1) = 2;
  switch (ref){
    case xfe_QuadRefUniform:
      (*pnvert) = 3;
      for (k=0; k<3; k++) coord[k] = xUniform[k];
      (*pnelem) = 2;
      for (k=0; k<4; k++) elem2vert[k] = eUniform[k];
      for (k=0; k<2; k++) elem2pos[k] = epUniform[k];
      for (k=0; k<3; k++) vert2corner[k] = vcUniform[k];
      (*pniface) = 1;
      for (k=0; k<(*pniface); k++) xf_SetIFaceValues(ifacelist + k, efoUniform[k]);
      (*nbface) = 2;
      for (k=0; k<(*nbface); k++){
        bfacelist[k].Elem = beUniform[k];
        bfacelist[k].Face = bfUniform[k];
        bface2pos[k] = bpUniform[k];
      }
      break;
    default:
      return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefineGenericElem_Quad
static int 
xf_RefineGenericElem_Quad(int ref, int *pnNodeQ1, real *coord, 
                          int *pnvert, int *vert2corner,
                          int *pnelem, int *elem2vert, int *elem2pos,
                          int *pniface, xf_IFace *ifacelist, int *nbface, 
                          xf_BFace *bfacelist, int *bface2pos)
{
  /*
   PURPOSE:
   
   Performs a generic refinement of a Quad according to ref.
   
   INPUTS:
   
   ref  : refinement indicator
   
   OUTPUTS: 
   
   (*pnNodeQ1) : number of Q1 nodes for Shape (and hence for subelements)
   coord : unrolled array of vertex coordinates
   (*pnvert) : number of vertices
   vert2corner : mapping from vertices to corners of original element
   (*pnelem) : number of subelements
   elem2vert: unrolled array of Q1 nodes for each element
   elem2pos : position of each new elem w.r.t original elem
   (*pniface) : number of interior faces
   ifacelist : list of element-interior faces
   (*nbface) : number of element-boundary faces
   bfacelist : list of element-boundary faces
   bface2pos : index for each bface identifying position w.r.t. original face
   
   RETURN:
   
   Error Code
   */
  int k;
  
  int eUniform[16] = {0,1,3,4, 1,2,4,5, 3,4,6,7, 4,5,7,8};
  int epUniform[4] = {xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE};
  int vcUniform[9] = {0,-1,1,-1,-1,-1,2,-1,3}; // vert2corner
  int beUniform[8] = {0,1,1,3,3,2,2,0}; // bface.elem
  int bfUniform[8] = {0,0,1,1,2,2,3,3}; // bface.face
  int bpUniform[8] = {xfe_SegPosLeft, xfe_SegPosRight, xfe_SegPosLeft, xfe_SegPosRight,
    xfe_SegPosLeft, xfe_SegPosRight, xfe_SegPosLeft, xfe_SegPosRight};
  int efoUniform[4][6] = {{0,1,0,1,3,1}, {0,2,1,2,0,0}, {1,2,1,3,0,0}, {2,1,0,3,3,1}};
  real xUniform[18] = {0,0, .5,0, 1,0, 0,.5, .5,.5, 1,.5, 0,1, .5,1, 1,1};
  
  
  int eHoriz[8] = {0,1,2,3, 2,3,4,5};
  int epHoriz[2] = {xfe_QuadPosBottom, xfe_QuadPosTop};
  int vcHoriz[6] = {0,1,-1,-1,2,3};
  int beHoriz[6] = {0,0,1,1,1,0}; // bface.elem
  int bfHoriz[6] = {0,1,1,2,3,3}; // bface.face
  int bpHoriz[6] = {xfe_SegPosNone, xfe_SegPosLeft, xfe_SegPosRight,
    xfe_SegPosNone, xfe_SegPosLeft, xfe_SegPosRight}; // bface2pos 
  real xHoriz[12] = {0,0, 1,0, 0,.5, 1,.5, 0,1, 1,1};
  
  int eVert[16] = {0,1,3,4, 1,2,4,5};
  int epVert[2] = {xfe_QuadPosLeft, xfe_QuadPosRight};
  int vcVert[6] = {0,-1,1,2,-1,3};
  int beVert[6] = {0,1,1,1,0,0}; // bface.elem
  int bfVert[6] = {0,0,1,2,2,3}; // bface.face
  int bpVert[6] = {xfe_SegPosLeft, xfe_SegPosRight,xfe_SegPosNone,
    xfe_SegPosLeft, xfe_SegPosRight,xfe_SegPosNone}; // bface2pos  
  real xVert[12] = {0,0, .5,0, 1,0, 0,1, .5,1, 1,1};
  
  
  (*pnNodeQ1) = 4;
  switch (ref){
    case xfe_QuadRefUniform:
      (*pnvert) = 9;
      for (k=0; k<18; k++) coord[k] = xUniform[k];
      (*pnelem) = 4;
      for (k=0; k<16; k++) elem2vert[k] = eUniform[k];
      for (k=0; k<4; k++) elem2pos[k] = epUniform[k];
      for (k=0; k<9; k++) vert2corner[k] = vcUniform[k];
      (*pniface) = 4;
      for (k=0; k<4; k++) xf_SetIFaceValues(ifacelist + k, efoUniform[k]);
      (*nbface) = 8;
      for (k=0; k<8; k++){
        bfacelist[k].Elem = beUniform[k]; 	
        bfacelist[k].Face = bfUniform[k];
        bface2pos[k] = bpUniform[k];
      }
      break;
    case xfe_QuadRefHoriz:
      (*pnvert) = 6;
      for (k=0; k<12; k++) coord[k] = xHoriz[k];
      (*pnelem) = 2;
      for (k=0; k<2; k++) elem2pos[k] = epHoriz[k];
      for (k=0; k<8; k++) elem2vert[k] = eHoriz[k];
      for (k=0; k<6; k++) vert2corner[k] = vcHoriz[k];
      (*pniface) = 1;
      ifacelist[0].ElemL = 0; ifacelist[0].FaceL = 2; ifacelist[0].OrientL = 1;
      ifacelist[0].ElemR = 1; ifacelist[0].FaceR = 0; ifacelist[0].OrientR = 0;
      (*nbface) = 6;
      for (k=0; k<6; k++){
        bfacelist[k].Elem = beHoriz[k]; 	
        bfacelist[k].Face = bfHoriz[k];
        bface2pos[k] = bpHoriz[k];
      }
      break;
    case xfe_QuadRefVert:
      (*pnvert) = 6;
      for (k=0; k<12; k++) coord[k] = xVert[k];
      (*pnelem) = 2;
      for (k=0; k<2; k++) elem2pos[k] = epVert[k];
      for (k=0; k<8; k++) elem2vert[k] = eVert[k];
      for (k=0; k<6; k++) vert2corner[k] = vcVert[k];
      (*pniface) = 1;
      ifacelist[0].ElemL = 0; ifacelist[0].FaceL = 1; ifacelist[0].OrientL = 0;
      ifacelist[0].ElemR = 1; ifacelist[0].FaceR = 3; ifacelist[0].OrientR = 1;
      (*nbface) = 6;
      for (k=0; k<6; k++){
        bfacelist[k].Elem = beVert[k]; 	
        bfacelist[k].Face = bfVert[k];
        bface2pos[k] = bpVert[k];
      }
      break;
    default:
      return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefineGeneric_Hex
static int 
xf_RefineGenericElem_Hex(int ref, int *pnNodeQ1, real *coord, 
                         int *pnvert, int *vert2corner,
                         int *pnelem, int *elem2vert, int *elem2pos,
                         int *pniface, xf_IFace *ifacelist, int *pnbface, 
                         xf_BFace *bfacelist, int *bface2pos, int *pnedge,
                         int *elist)
{
  /*
   PURPOSE:
   
   Performs a generic refinement of a Hex according to ref.
   
   INPUTS:
   
   ref  : refinement indicator (enumerated)
   
   OUTPUTS: 
   
   (*pnNodeQ1) : number of Q1 nodes for Shape (and hence for subelements)
   coord : unrolled array of vertex coordinates
   (*pnvert) : number of vertices
   vert2corner : mapping from vertices to corners of original element
   (*pnelem) : number of subelements
   elem2vert: unrolled array of Q1 nodes for each element
   elem2pos : position of each new elem w.r.t original elem
   (*pniface) : number of interior faces
   ifacelist : list of element-interior faces
   (*pnbface) : number of element-boundary faces
   bfacelist : list of element-boundary faces
   bface2pos : index for each bface identifying position w.r.t. original face
   (*pnedge) : number of bisected edges 
   (*elist)  : list of three-tuples for each bisected edge, containing [n0, n1, nm],
   where [n0, n1] are the local node numbers of the edge vertices,
   and nm is the local node number of the midpoint.
   
   RETURN:
   
   Error Code
   */
  int k, ierr;
  
  // **** Uniform ****
  
  // for elem2vert
  int eUniform[64] = { 0, 1, 3, 4, 9,10,12,13,  1, 2, 4, 5,10,11,13,14,  
    3, 4, 6, 7,12,13,15,16,  4, 5, 7, 8,13,14,16,17,
    9,10,12,13,18,19,21,22, 10,11,13,14,19,20,22,23,
    12,13,15,16,21,22,24,25, 13,14,16,17,22,23,25,26};
  // for elem2pos
  int epUniform[8] = {xfe_HexPos000, xfe_HexPos100, xfe_HexPos010, xfe_HexPos110, 
    xfe_HexPos001, xfe_HexPos101, xfe_HexPos011, xfe_HexPos111};
  // for vert2corner
  int vcUniform[27] = {0,-1,1,-1,-1,-1,2,-1,3, -1,-1,-1,-1,-1,-1,-1,-1,-1, 
    4,-1,5,-1,-1,-1,6,-1,7};
  // for bfacelist.elem
  int beUniform[24] = {0,2,1,3, 0,1,4,5, 1,3,5,7, 3,2,7,6, 2,0,6,4, 4,5,6,7};
  // for bfacelist.face
  int bfUniform[24] = {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5};
  // for bface2pos (relative to elem being refined)
  int bpUniform[24] = {xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE,
    xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE,
    xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE,
    xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE,
    xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE,
    xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE};
  // for ifacelist: elem, face, orient (L and R)
  int efoUniform[12][6] = {{0,2,0, 1,4,5}, {1,3,5, 3,1,0}, {3,4,5, 2,2,0}, {2,1,0, 0,3,5},
    {0,5,0, 4,0,4}, {1,5,0, 5,0,4}, {2,5,0, 6,0,4}, {3,5,0, 7,0,4},
    {4,2,0, 5,4,5}, {5,3,5, 7,1,0}, {7,4,5, 6,2,0}, {6,1,0, 4,3,5}};
  real xUniform[81] = {0,0, 0, .5,0, 0, 1,0, 0, 0,.5, 0, .5,.5, 0, 1,.5, 0, 0,1, 0, .5,1, 0, 1,1, 0,
    0,0,.5, .5,0,.5, 1,0,.5, 0,.5,.5, .5,.5,.5, 1,.5,.5, 0,1,.5, .5,1,.5, 1,1,.5,
    0,0, 1, .5,0, 1, 1,0, 1, 0,.5, 1, .5,.5, 1, 1,.5, 1, 0,1, 1, .5,1, 1, 1,1, 1};
  int elUniform[36] = {0,2,1, 2,8,5, 8,6,7, 6,0,3, 0,18,9, 2,20,11, 8,26,17, 6,24,15,
    18,20,19, 20,26,23, 26,24,25, 24,18,21};
  
  // **** SliceX ****
  // for elem2vert
  int eSliceX[16] = {0,1,3,4,6,7,9,10,  1,2,4,5,7,8,10,11};
  // for elem2pos
  int epSliceX[2] = {xfe_HexPos022, xfe_HexPos122};
  // for vert2corner
  int vcSliceX[12] = {0,-1,1,2,-1,3, 4,-1,5,6,-1,7};
  // for bfacelist.elem
  int beSliceX[10] = {0,1,0,1,1,1,0,0,0,1};
  // for bfacelist.face
  int bfSliceX[10] = {0,0, 1,1, 2, 3,3, 4, 5,5};
  // for bface2pos (relative to elem being refined)
  int bpSliceX[10] = {xfe_QuadPosBottom, xfe_QuadPosTop, xfe_QuadPosLeft, xfe_QuadPosRight,
    xfe_QuadPosNone, xfe_QuadPosLeft, xfe_QuadPosRight, xfe_QuadPosNone,
    xfe_QuadPosLeft, xfe_QuadPosRight};
  // for ifacelist: elem, face, orient (L and R)
  int efoSliceX[1][6] = {{0,2,0, 1,4,5}};   
  // coordinates, unrolled
  real xSliceX[36] = {0,0,0, .5,0,0, 1,0,0, 0,1,0, .5,1,0, 1,1,0,
    0,0,1, .5,0,1, 1,0,1, 0,1,1, .5,1,1, 1,1,1};
  int elSliceX[12] = {0,2,1, 3,5,4, 6,8,7, 9,11,10};
  
  // **** SliceY ****
  // for elem2vert
  int eSliceY[16] = {0,1,2,3,6,7,8,9,  2,3,4,5,8,9,10,11};
  // for elem2pos
  int epSliceY[2] = {xfe_HexPos202, xfe_HexPos212};
  // for vert2corner
  int vcSliceY[12] = {0,1,-1,-1,2,3, 4,5,-1,-1,6,7};
  // for bfacelist.elem
  int beSliceY[10] = {0,1,0,0,1,1,1,0,0,1};
  // for bfacelist.face
  int bfSliceY[10] = {0,0, 1, 2,2, 3, 4,4, 5,5};
  // for bface2pos (relative to elem being refined)
  int bpSliceY[10] = {xfe_QuadPosLeft, xfe_QuadPosRight, xfe_QuadPosNone, xfe_QuadPosLeft,
    xfe_QuadPosRight, xfe_QuadPosNone, xfe_QuadPosLeft, xfe_QuadPosRight,
    xfe_QuadPosBottom, xfe_QuadPosTop};
  // for ifacelist: elem, face, orient (L and R)
  int efoSliceY[1][6] = {{0,3,5, 1,1,0}};   
  // coordinates, unrolled
  real xSliceY[36] = {0,0,0, 1,0,0, 0,.5,0, 1,.5,0, 0,1,0, 1,1,0,
    0,0,1, 1,0,1, 0,.5,1, 1,.5,1, 0,1,1, 1,1,1};
  int elSliceY[12] = {0,4,2, 1,5,3, 6,10,8, 7,11,9};
  
  // **** SliceZ ****
  // for elem2vert
  int eSliceZ[16] = {0,1,2,3,4,5,6,7, 4,5,6,7,8,9,10,11};
  // for elem2pos
  int epSliceZ[2] = {xfe_HexPos220, xfe_HexPos221};
  // for vert2corner
  int vcSliceZ[12] = {0,1,2,3, -1,-1,-1,-1, 4,5,6,7};
  // for bfacelist.elem
  int beSliceZ[10] = {0,0,1,0,1,0,1,0,1,1};
  // for bfacelist.face
  int bfSliceZ[10] = {0, 1,1, 2,2, 3,3, 4,4, 5};
  // for bface2pos (relative to elem being refined)
  int bpSliceZ[10] = {xfe_QuadPosNone, xfe_QuadPosBottom, xfe_QuadPosTop, xfe_QuadPosBottom,
    xfe_QuadPosTop, xfe_QuadPosBottom, xfe_QuadPosTop, xfe_QuadPosBottom,
    xfe_QuadPosTop, xfe_QuadPosNone};
  // for ifacelist: elem, face, orient (L and R)
  int efoSliceZ[1][6] = {{0,5,0, 1,0,4}};
  // coordinates, unrolled
  real xSliceZ[36] = {0,0, 0, 1,0, 0, 0,1, 0, 1,1, 0, 
    0,0,.5, 1,0,.5, 0,1,.5, 1,1,.5,
    0,0, 1, 1,0, 1, 0,1, 1, 1,1, 1};
  int elSliceZ[12] = {0,8,4, 1,9,5, 2,10,6, 3,11,7};
  
  // **** SliceXY ****
  // for elem2vert
  int eSliceXY[32] = {0,1,3,4,9,10,12,13,  1,2,4,5,10,11,13,14,
    3,4,6,7,12,13,15,16, 4,5,7,8,13,14,16,17};
  // for elem2pos
  int epSliceXY[4] = {xfe_HexPos002, xfe_HexPos102, xfe_HexPos012, xfe_HexPos112};
  // for vert2corner
  int vcSliceXY[18] = {0,-1,1,-1,-1,-1,2,-1,3, 4,-1,5,-1,-1,-1,6,-1,7};
  // for bfacelist.elem
  int beSliceXY[16] = {0,2,1,3,0,1,1,3,3,2,2,0,0,1,2,3};
  // for bfacelist.face
  int bfSliceXY[16] = {0,0,0,0, 1,1, 2,2, 3,3, 4,4, 5,5,5,5};
  // for bface2pos (relative to elem being refined)
  int bpSliceXY[16] = {xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE,
    xfe_QuadPosLeft, xfe_QuadPosRight, xfe_QuadPosLeft, xfe_QuadPosRight,
    xfe_QuadPosLeft, xfe_QuadPosRight, xfe_QuadPosLeft, xfe_QuadPosRight,
    xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE};
  // for ifacelist: elem, face, orient (L and R)
  int efoSliceXY[4][6] = {{0,2,0, 1,4,5}, {1,3,5, 3,1,0}, {3,4,5, 2,2,0}, {2,1,0, 0,3,5}};   
  // coordinates, unrolled
  real xSliceXY[54] = {0,0,0, .5,0,0, 1,0,0, 0,.5,0, .5,.5,0, 1,.5,0, 0,1,0, .5,1,0, 1,1,0,
    0,0,1, .5,0,1, 1,0,1, 0,.5,1, .5,.5,1, 1,.5,1, 0,1,1, .5,1,1, 1,1,1};
  int elSliceXY[24] = {0,2,1, 2,8,5, 8,6,7, 6,0,3, 9,11,10, 11,17,14, 17,15,16, 15,9,12};
  
  // **** SliceXZ ****
  // for elem2vert
  int eSliceXZ[32] = {0,1,3,4,6,7,9,10,  1,2,4,5,7,8,10,11,
    6,7,9,10,12,13,15,16, 7,8,10,11,13,14,16,17};
  // for elem2pos
  int epSliceXZ[4] = {xfe_HexPos020, xfe_HexPos120, xfe_HexPos021, xfe_HexPos121};
  // for vert2corner
  int vcSliceXZ[18] = {0,-1,1,2,-1,3,-1,-1,-1,-1,-1,-1,4,-1,5,6,-1,7};
  // for bfacelist.elem
  int beSliceXZ[16] = {0,1,0,1,2,3,1,3,1,0,3,2,0,2,2,3};
  // for bfacelist.face
  int bfSliceXZ[16] = {0,0, 1,1,1,1, 2,2, 3,3,3,3, 4,4, 5,5};
  // for bface2pos (relative to elem being refined)
  int bpSliceXZ[16] = {xfe_QuadPosBottom, xfe_QuadPosTop, xfe_QuadPosSW, xfe_QuadPosSE, 
    xfe_QuadPosNW, xfe_QuadPosNE, xfe_QuadPosBottom, xfe_QuadPosTop,
    xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE,
    xfe_QuadPosBottom, xfe_QuadPosTop, xfe_QuadPosLeft, xfe_QuadPosRight};
  // for ifacelist: elem, face, orient (L and R)
  int efoSliceXZ[4][6] = {{0,2,0, 1,4,5}, {1,5,0, 3,0,4}, {3,4,5, 2,2,0}, {2,0,4, 0,5,0}};
  // coordinates, unrolled
  real xSliceXZ[54] = {0,0, 0, .5,0, 0, 1,0, 0,  0,1, 0, .5,1, 0, 1,1, 0,
    0,0,.5, .5,0,.5, 1,0,.5,  0,1,.5, .5,1,.5, 1,1,.5,
    0,0, 1, .5,0, 1, 1,0, 1,  0,1, 1, .5,1, 1, 1,1, 1};
  int elSliceXZ[24] = {0,2,1, 3,5,4, 0,12,6, 2,14,8, 3,15,9, 5,17,11, 12,14,13, 15,17,16};
  
  // **** SliceYZ ****
  // for elem2vert
  int eSliceYZ[32] = {0,1,2,3,6,7,8,9,  2,3,4,5,8,9,10,11,
    6,7,8,9,12,13,14,15, 8,9,10,11,14,15,16,17};
  // for elem2pos
  int epSliceYZ[4] = {xfe_HexPos200, xfe_HexPos210, xfe_HexPos201, xfe_HexPos211};
  // for vert2corner
  int vcSliceYZ[18] = {0,1,-1,-1,2,3,-1,-1,-1,-1,-1,-1,4,5,-1,-1,6,7};
  // for bfacelist.elem
  int beSliceYZ[16] = {0,1,0,2,0,1,2,3,1,3,1,0,3,2,2,3};
  // for bfacelist.face
  int bfSliceYZ[16] = {0,0, 1,1, 2,2,2,2, 3,3, 4,4,4,4, 5,5};
  // for bface2pos (relative to elem being refined)
  int bpSliceYZ[16] = {xfe_QuadPosLeft, xfe_QuadPosRight, xfe_QuadPosBottom, xfe_QuadPosTop,
    xfe_QuadPosSW, xfe_QuadPosSE, xfe_QuadPosNW, xfe_QuadPosNE, 
    xfe_QuadPosBottom, xfe_QuadPosTop, xfe_QuadPosSW, xfe_QuadPosSE, 
    xfe_QuadPosNW, xfe_QuadPosNE, xfe_QuadPosBottom, xfe_QuadPosTop};
  // for ifacelist: elem, face, orient (L and R)
  int efoSliceYZ[4][6] = {{0,3,5, 1,1,0}, {1,5,0, 3,0,4}, {3,1,0, 2,3,5}, {2,0,4, 0,5,0}};
  // coordinates, unrolled
  real xSliceYZ[54] = {0,0, 0, 1,0, 0, 0,.5, 0, 1,.5, 0, 0,1, 0, 1,1, 0,
    0,0,.5, 1,0,.5, 0,.5,.5, 1,.5,.5, 0,1,.5, 1,1,.5,
    0,0, 1, 1,0, 1, 0,.5, 1, 1,.5, 1, 0,1, 1, 1,1, 1};
  int elSliceYZ[24] = {0,4,2, 1,5,3, 0,12,6, 1,13,7, 4,16,10, 5,17,11, 12,16,14, 13,17,15};
  
  
  (*pnNodeQ1) = 8;
  switch (ref){
    case xfe_HexRefUniform:
      (*pnvert)  = 27;
      (*pnelem)  = 8;
      (*pniface) = 12;
      (*pnbface) = 24;
      (*pnedge)  = 12;
      for (k=0; k<(*pnvert)*3; k++) coord[k]       = xUniform[k];
      for (k=0; k<(*pnelem)*8; k++) elem2vert[k]   = eUniform[k];
      for (k=0; k<(*pnelem)  ; k++) elem2pos[k]    = epUniform[k];
      for (k=0; k<(*pnvert)  ; k++) vert2corner[k] = vcUniform[k];
      for (k=0; k<(*pniface); k++) xf_SetIFaceValues(ifacelist + k, efoUniform[k]);
      for (k=0; k<(*pnbface); k++){
        bfacelist[k].Elem = beUniform[k]; 	
        bfacelist[k].Face = bfUniform[k];
        bface2pos[k]      = bpUniform[k];
      }
      for (k=0; k<(*pnedge)*3; k++) elist[k] = elUniform[k];
      break;
    case xfe_HexRefSliceX:
      (*pnvert)  = 12;
      (*pnelem)  = 2;
      (*pniface) = 1;
      (*pnbface)  = 10;
      (*pnedge)  = 4;
      for (k=0; k<(*pnvert)*3; k++) coord[k]       = xSliceX[k];
      for (k=0; k<(*pnelem)*8; k++) elem2vert[k]   = eSliceX[k];
      for (k=0; k<(*pnelem)  ; k++) elem2pos[k]    = epSliceX[k];
      for (k=0; k<(*pnvert)  ; k++) vert2corner[k] = vcSliceX[k];
      for (k=0; k<(*pniface); k++) xf_SetIFaceValues(ifacelist + k, efoSliceX[k]);
      for (k=0; k<(*pnbface); k++){
        bfacelist[k].Elem = beSliceX[k]; 	
        bfacelist[k].Face = bfSliceX[k];
        bface2pos[k]      = bpSliceX[k];
      }
      for (k=0; k<(*pnedge)*3; k++) elist[k] = elSliceX[k];
      break;
    case xfe_HexRefSliceY:
      (*pnvert)  = 12;
      (*pnelem)  = 2;
      (*pniface) = 1;
      (*pnbface) = 10;
      (*pnedge)  = 4;
      for (k=0; k<(*pnvert)*3; k++) coord[k]       = xSliceY[k];
      for (k=0; k<(*pnelem)*8; k++) elem2vert[k]   = eSliceY[k];
      for (k=0; k<(*pnelem)  ; k++) elem2pos[k]    = epSliceY[k];
      for (k=0; k<(*pnvert)  ; k++) vert2corner[k] = vcSliceY[k];
      for (k=0; k<(*pniface); k++) xf_SetIFaceValues(ifacelist + k, efoSliceY[k]);
      for (k=0; k<(*pnbface); k++){
        bfacelist[k].Elem = beSliceY[k]; 	
        bfacelist[k].Face = bfSliceY[k];
        bface2pos[k]      = bpSliceY[k];
      }
      for (k=0; k<(*pnedge)*3; k++) elist[k] = elSliceY[k];
      break;
    case xfe_HexRefSliceZ:
      (*pnvert)  = 12;
      (*pnelem)  = 2;
      (*pniface) = 1;
      (*pnbface) = 10;
      (*pnedge)  = 4;
      for (k=0; k<(*pnvert)*3; k++) coord[k]       = xSliceZ[k];
      for (k=0; k<(*pnelem)*8; k++) elem2vert[k]   = eSliceZ[k];
      for (k=0; k<(*pnelem)  ; k++) elem2pos[k]    = epSliceZ[k];
      for (k=0; k<(*pnvert)  ; k++) vert2corner[k] = vcSliceZ[k];
      for (k=0; k<(*pniface); k++) xf_SetIFaceValues(ifacelist + k, efoSliceZ[k]);
      for (k=0; k<(*pnbface); k++){
        bfacelist[k].Elem = beSliceZ[k]; 	
        bfacelist[k].Face = bfSliceZ[k];
        bface2pos[k]      = bpSliceZ[k];
      }
      for (k=0; k<(*pnedge)*3; k++) elist[k] = elSliceZ[k];
      break;
    case xfe_HexRefSliceXY:
      (*pnvert)  = 18;
      (*pnelem)  = 4;
      (*pniface) = 4;
      (*pnbface) = 16;
      (*pnedge)  = 8;
      for (k=0; k<(*pnvert)*3; k++) coord[k]       = xSliceXY[k];
      for (k=0; k<(*pnelem)*8; k++) elem2vert[k]   = eSliceXY[k];
      for (k=0; k<(*pnelem)  ; k++) elem2pos[k]    = epSliceXY[k];
      for (k=0; k<(*pnvert)  ; k++) vert2corner[k] = vcSliceXY[k];
      for (k=0; k<(*pniface); k++) xf_SetIFaceValues(ifacelist + k, efoSliceXY[k]);
      for (k=0; k<(*pnbface); k++){
        bfacelist[k].Elem = beSliceXY[k]; 	
        bfacelist[k].Face = bfSliceXY[k];
        bface2pos[k]      = bpSliceXY[k];
      }
      for (k=0; k<(*pnedge)*3; k++) elist[k] = elSliceXY[k];
      break;
    case xfe_HexRefSliceXZ:
      (*pnvert)  = 18;
      (*pnelem)  = 4;
      (*pniface) = 4;
      (*pnbface) = 16;
      (*pnedge)  = 8;
      for (k=0; k<(*pnvert)*3; k++) coord[k]       = xSliceXZ[k];
      for (k=0; k<(*pnelem)*8; k++) elem2vert[k]   = eSliceXZ[k];
      for (k=0; k<(*pnelem)  ; k++) elem2pos[k]    = epSliceXZ[k];
      for (k=0; k<(*pnvert)  ; k++) vert2corner[k] = vcSliceXZ[k];
      for (k=0; k<(*pniface); k++) xf_SetIFaceValues(ifacelist + k, efoSliceXZ[k]);
      for (k=0; k<(*pnbface); k++){
        bfacelist[k].Elem = beSliceXZ[k]; 	
        bfacelist[k].Face = bfSliceXZ[k];
        bface2pos[k]      = bpSliceXZ[k];
      }
      for (k=0; k<(*pnedge)*3; k++) elist[k] = elSliceXZ[k];
      break;
    case xfe_HexRefSliceYZ:
      (*pnvert)  = 18;
      (*pnelem)  = 4;
      (*pniface) = 4;
      (*pnbface) = 16;
      (*pnedge)  = 8;
      for (k=0; k<(*pnvert)*3; k++) coord[k]       = xSliceYZ[k];
      for (k=0; k<(*pnelem)*8; k++) elem2vert[k]   = eSliceYZ[k];
      for (k=0; k<(*pnelem)  ; k++) elem2pos[k]    = epSliceYZ[k];
      for (k=0; k<(*pnvert)  ; k++) vert2corner[k] = vcSliceYZ[k];
      for (k=0; k<(*pniface); k++) xf_SetIFaceValues(ifacelist + k, efoSliceYZ[k]);
      for (k=0; k<(*pnbface); k++){
        bfacelist[k].Elem = beSliceYZ[k]; 	
        bfacelist[k].Face = bfSliceYZ[k];
        bface2pos[k]      = bpSliceYZ[k];
      }
      for (k=0; k<(*pnedge)*3; k++) elist[k] = elSliceYZ[k];
      break;
      
      
    default:
      return xf_Error(xf_INPUT_ERROR);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_RefineGenericElem
static int 
xf_RefineGenericElem(enum xfe_ShapeType Shape, int ref, int *pnNodeQ1, 
                     real *coord, int *pnvert, int *vert2corner,
                     int *pnelem, int *elem2vert, int *elem2pos,
                     int *pniface, xf_IFace *ifacelist, int *nbface, 
                     xf_BFace *bfacelist, int *bface2pos, int *pnedge,
                     int *elist)
{
  /*
   PURPOSE:
   
   Performs a generic refinement of Shape according to ref.
   
   INPUTS:
   
   Shape : shape which to refine
   ref  : refinement indicator
   
   OUTPUTS: 
   
   (*pnNodeQ1) : number of Q1 nodes for Shape (and hence for subelements)
   coord : unrolled array of vertex coordinates
   (*pnvert) : number of vertices
   vert2corner : mapping from vertices to corners of original element
   (*pnelem) : number of subelements
   elem2vert: unrolled array of Q1 nodes for each element
   elem2pos : position of each new elem w.r.t original elem
   (*pniface) : number of interior faces
   ifacelist : list of element-interior faces
   (*nbface) : number of element-boundary faces
   bfacelist : list of element-boundary faces
   bface2pos : index for each bface identifying position w.r.t. original face
   (*pnedge) : number of bisected edges (relevant for 3D only)
   (*elist)  : list of three-tuples for each bisected edge, containing [n0, n1, nm],
   where [n0, n1] are the local node numbers of the edge vertices,
   and nm is the local node number of the midpoint.
   
   RETURN:
   
   Error Code
   */
  int ierr;
  
  switch (Shape){
    case xfe_Segment:
      (*pnedge) = 0;
      ierr = xf_Error(xf_RefineGenericElem_Seg(ref, pnNodeQ1, coord, pnvert, vert2corner, 
                                               pnelem, elem2vert, elem2pos, pniface,
                                               ifacelist, nbface, bfacelist, bface2pos));
      break;
    case xfe_Quadrilateral:
      (*pnedge) = 0;
      ierr = xf_Error(xf_RefineGenericElem_Quad(ref, pnNodeQ1, coord, pnvert, vert2corner, 
                                                pnelem, elem2vert, elem2pos, pniface,
                                                ifacelist, nbface, bfacelist, bface2pos));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_Hexahedron:
      ierr = xf_Error(xf_RefineGenericElem_Hex(ref, pnNodeQ1, coord, pnvert, vert2corner, 
                                               pnelem, elem2vert, elem2pos, pniface,
                                               ifacelist, nbface, bfacelist, bface2pos,
                                               pnedge, elist));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_RefineElemHang
static int 
xf_RefineElemHang(xf_Mesh *Mesh, int egrp, int elem, int ref, int nnew, 
                  int *vnew, int *vpos, int nNodeOrig, xf_IntPairList *pSameNode,
                  int *node2ehash, xf_EdgeHash *edgehash, int *pnedgehash)
{
  /*
   PURPOSE:
   
   Refines (egrp,elem) in Mesh and returns a valid, connected Mesh
   structure as a result.  If adapting several elements, this function
   should be called on the highest priority elements first, to ensure
   that we can always produce a valid mesh (e.g. do not try to adapt a
   fine element of a hanging face before adapting the coarse element).
   Mesh must have structures already over-allocated to accomodate the
   new elements, faces, nodes, etc.
   
   INPUTS:
   
   Mesh : Mesh structure
   egrp, elem : element to refine
   ref  : refinement indicator for this element
   nnew : number of new elements to be created (mainly for code logic check)
   vnew : list of new element numbers to use
   node2ehash : used for accessing the edge bisection hash (3D specific, modified))
   edgehash : contains node numbers of bisected edges (modified)
   pnedgehash : number of edges in hash (modified)
   
   
   OUTPUTS: 
   
   vpos : list of position numbers of new elems w.r.t old elem
   pSameNode : pairs of identical nodes with different indices (not used)
   
   RETURN:
   
   Error Code
   */
  
  int ierr, dim, i, j, k, d, nn, nn1;
  int nvert, nelem, niface, nbface;
  int ielem, face, ibfgrp, ibface;
  int nelem0, face0, hang, pos, face0C;
  int eL, faceL, egR, iiface, iifacenew;
  int eI, faceI, face0I, bpos, orientI, eR;
  int nNode, nNodeQ1, QOrder, nFaceOrig, nface0;
  int QBasisF, QOrderF, faceF, orientF, egF, eF;
  int QBasisC, QOrderC, faceC, orientC, egC, eC;
  int relpos, hangnew;
  int *elem2glob, nNode0, iglob, ivert;
  int nFaceToOrient = 0;
  int orientN, QBasisI, QBasisN, QOrderI, QOrderN, faceN;
  int nfacemax = ha_MAX_BFACE;
  int elem2vert[ha_MAX_ELEM*ha_MAX_NODEQ1];
  int elem2pos[ha_MAX_ELEM];
  int bface2pos[ha_MAX_BFACE];
  int vert2corner[ha_MAX_VERT];
  int vert2glob[ha_MAX_VERT];
  int corner2glob[xf_MAXQ1NODE];
  int nToOrient, nToOrientMax, **ToOrient;
  int nvecQ1[xf_MAXQ1NODE], vvec[xf_MAXQ1NODE], nvec[xf_MAXQ1NODE];
  int FaceToOrient[ha_MAX_BFACE][2];
  int CommonNode[xf_MAXQ1FACENODE];
  int  hangOrig[ha_MAX_BFACE];
  int   posOrig[ha_MAX_BFACE];
  int face0Orig[ha_MAX_BFACE];
  int nedge, elist[ha_MAX_EDGE*3];
  int n0, n1, ihash;
  enum xfe_Bool IamL, compatible;
  enum xfe_BasisType QBasis, QBasisLag;
  enum xfe_ShapeType Shape, FShape;
  enum xfe_RelType rel;
  real coord[ha_MAX_VERT*3];
  real *xref0, *xref, *xglob;
  xf_IFace ifacelist[ha_MAX_IFACE], IFace;
  xf_BFace bfacelist[ha_MAX_BFACE];
  xf_Face *FaceOrig, Face;
  xf_BasisData *PhiData, *PhiData1;
  
  
  dim = Mesh->Dim;
  
  // should not be called with no refinement request
  if (ref == 0) return xf_Error(xf_INPUT_ERROR);
  
  // get element Shape
  ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  // number of faces in Shape
  ierr = xf_Error(xf_Shape2nFace(Shape, &nface0));
  if (ierr != xf_OK) return ierr;
  
  // number of faces in egrp,elem and their .Face values
  nFaceOrig = Mesh->ElemGroup[egrp].nFace[elem];
  ierr = xf_Error(xf_Alloc( (void **) &FaceOrig, nFaceOrig, sizeof(xf_Face)));
  if (ierr != xf_OK) return ierr;
  for (face=0; face<nFaceOrig; face++) 
    FaceOrig[face] = Mesh->ElemGroup[egrp].Face[elem][face];
  
  
  // determine hang, pos, etc. on original faces
  if (nFaceOrig > ha_MAX_BFACE) return xf_Error(xf_CODE_LOGIC_ERROR);
  for (face=0; face<nFaceOrig; face++){
    // get original face (while checking if hanging face)
    ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, hangOrig+face, 
                                     face0Orig+face, posOrig+face, NULL));
    if (ierr != xf_OK) return ierr;
    
  }
  
  
  // vert2corner = mapping from vertices to corners of original element
  for (k=0; k<ha_MAX_VERT; k++) vert2corner[k] = -1;
  
  // get generic refinement information
  ierr = xf_Error(xf_RefineGenericElem(Shape, ref, &nNodeQ1, coord, &nvert, vert2corner, 
                                       &nelem, elem2vert, elem2pos, &niface, ifacelist, 
                                       &nbface, bfacelist, bface2pos, &nedge, elist));
  if (ierr != xf_OK) return ierr;
  
  
  // check number of elements
  if (nelem != nnew) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  // number of nodes per element
  nNode = Mesh->ElemGroup[egrp].nNode;
  
  // make sure Basis used to represent element shape is Lagrange
  QBasis = Mesh->ElemGroup[egrp].QBasis;
  ierr = xf_Error(xf_Basis2Lagrange(QBasis, &QBasisLag));
  if (ierr != xf_OK) return ierr;
  if (QBasis != QBasisLag) return xf_Error(xf_NOT_SUPPORTED);
  
  // Get ref-space coords of element nodes
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  xref0 = NULL;
  ierr = xf_Error(xf_LagrangeNodes(QBasis, QOrder, NULL, NULL, &xref0));
  if (ierr != xf_OK) return ierr;
  
  // Q1 basis at xref0
  PhiData1 = NULL;
  ierr = xf_Error(xf_EvalBasis(QBasis, 1, xfe_True, nNode, xref0, xfb_Phi, &PhiData1));
  if (ierr != xf_OK) return ierr;
  
  
  // allocate space for mapped nodes in ref space and glob space
  ierr = xf_Error(xf_Alloc( (void **) &xref, dim*nNode, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &xglob, nelem*dim*nNode, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  /*-------*/
  /* Nodes */
  /*-------*/
  
  // allocate space for storing new element numbers
  ierr = xf_Error(xf_Alloc( (void **) &elem2glob, nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  PhiData = NULL;
  
  //nelem0 = Mesh->ElemGroup[egrp].nElem; // current number of elements
  
  /* Map all nodes to physical space */
  for (i=0; i<nelem; i++){
    
    // evaluate original elem ref coords of new elem ref nodes
    for (k=0; k<nNode; k++){
      for (d=0; d<dim; d++) xref[k*dim+d] = 0;
      for (j=0; j<nNodeQ1; j++)
        for (d=0; d<dim; d++)
          xref[k*dim+d] += coord[dim*elem2vert[nNodeQ1*i+j]+d]*PhiData1->Phi[k*nNodeQ1+j];
    }
    
    // map xref to physical space -> xglob
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &PhiData, xfe_True,
                                    nNode, xref, xglob+i*dim*nNode));
    if (ierr != xf_OK) return ierr;
    
    ielem = vnew[i]; //((i == 0) ? elem : nelem0+i-1);
    elem2glob[i] = ielem;
    Mesh->ElemGroup[egrp].nFace[ielem] = nface0;
    vpos[ielem] = elem2pos[i]; // position of new element
    
    // initialize original faces to null
    for (j=0; j<nface0; j++){
      Mesh->ElemGroup[egrp].Face[ielem][j].Group  = xf_NULLFACE;
      Mesh->ElemGroup[egrp].Face[ielem][j].Number = -1; // indicates not set
    }
  } // i
  
  
  // determine original Q1 nodes
  ierr = xf_Error(xf_Q1Nodes(QBasis, QOrder, &nn, nvecQ1));
  if (ierr != xf_OK) return ierr;
  if (nn != nNodeQ1) return xf_Error(xf_CODE_LOGIC_ERROR);
  for (k=0; k<nNodeQ1; k++)
    corner2glob[k] = Mesh->ElemGroup[egrp].Node[elem][nvecQ1[k]];
  
  // vert2glob = mapping from each vert to an existing or a new global node
  for (k=0; k<nvert; k++) vert2glob[k] = -1;
  
  // first set Q1 nodes of original element
  for (k=0; k<nvert; k++)
    if (vert2corner[k] >= 0){
      if (ah_DEBUG)
        xf_printf("vert %d is on corner %d, where glob =%d\n",
                  k, vert2corner[k], corner2glob[vert2corner[k]]);
      vert2glob[k] = corner2glob[vert2corner[k]];
    }
  
  // set Q1 nodes from existing bisected edges using hash
  for (i=0; i<nedge; i++){
    // two edge nodes
    n0 = vert2glob[elist[3*i+0]];
    n1 = vert2glob[elist[3*i+1]];
    if (n0 > n1) swap(n0, n1, j);
    if (n0 < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
    ierr = xf_EdgeHashCheck(n0, n1, node2ehash, edgehash, &ihash);
    if (ierr == xf_OK) // edge bisection exists, set glob node number
      vert2glob[elist[3*i+2]] = edgehash[ihash].idata;
  }
  
  
  // set Q1 nodes on any new hanging faces while looping over faces
  
  /*-------*/
  /* Faces */
  /*-------*/
  
  
  if (ah_DEBUG){
    xf_printf("\n\n &&& nFaceOrig = %d &&&\n", nFaceOrig);
    xf_printf("    Starting out: Mesh->nIFace = %d, ref = %d\n", Mesh->nIFace, ref);
  }
  
  nToOrient    = 0;
  nToOrientMax = 4*nbface;
  ierr = xf_Error(xf_Alloc2( (void ***) &ToOrient, nToOrientMax, 3, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // add boundary faces (tie to neighbors)
  for (face=0; face<nFaceOrig; face++){ // loop over original elem faces
    
    Face  = FaceOrig[face];
    hang  = hangOrig[face];
    face0 = face0Orig[face];
    pos   = posOrig[face];
    
    if ((ibfgrp = Face.Group) == xf_NULLFACE){
      /* Do not do anything in the case of null original faces --
       faces of the subelements are already initialized to
       xf_NULLFACE before this loop.*/
      continue;
    } 
    else if (ibfgrp >= 0){ // on domain boundary
      face0 = face;
      // get corresponding face or faces from bfacelist
      ibface = Face.Number;
      for (i=0; i<nbface; i++){
        if ( bfacelist[i].Face == face0){ // refined elements share same face numbering
          eL    = elem2glob[bfacelist[i].Elem];
          faceL = bfacelist[i].Face;
          Mesh->BFaceGroup[ibfgrp].BFace[ibface].ElemGroup = egrp;
          Mesh->BFaceGroup[ibfgrp].BFace[ibface].Elem      = eL;
          Mesh->BFaceGroup[ibfgrp].BFace[ibface].Face      = faceL;
          Mesh->BFaceGroup[ibfgrp].BFace[ibface].Orient    = 0; // irrelevant for bfaces
          Mesh->ElemGroup[egrp].Face[eL][faceL].Group  = ibfgrp;
          Mesh->ElemGroup[egrp].Face[eL][faceL].Number = ibface;
          if (ah_DEBUG)
            xf_printf(" Added boundary to eL=%d, faceL=%d, ibfgrp=%d, ibface=%d\n",
                      eL, faceL, ibfgrp, ibface);
          if (ibface == Face.Number) ibface = Mesh->BFaceGroup[ibfgrp].nBFace;
          else { // else, need to augment number of bfaces
            Mesh->BFaceGroup[ibfgrp].nBFace++;
            ibface++;
          }
        }
      } // i
    }
    else if (ibfgrp == xf_INTERIORFACE){ // on an interior face
      iiface = Face.Number;
      IFace = Mesh->IFace[iiface];
      /*    // get original face (while checking if hanging face) */
      /*       ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, &face0, &pos, NULL)); */
      /*       if (ierr != xf_OK) return ierr; */
      
      // Is egrp,elem on left?
      ierr = xf_Error(xf_IsElemOnLeft(IFace, egrp, elem, &IamL));
      if (ierr != xf_OK) return ierr;
      
      // shape of element
      ierr = xf_Error(xf_Basis2Shape((IamL) ? Mesh->ElemGroup[IFace.ElemGroupL].QBasis : 
                                     Mesh->ElemGroup[IFace.ElemGroupR].QBasis, &Shape));
      if (ierr != xf_OK) return ierr;
      
      // get Shape of face = Shape of face0
      ierr = xf_Error(xf_FaceShape(Shape, face0, &FShape));
      if (ierr != xf_OK) return ierr;
      
      if (ah_DEBUG){
        xf_printf("\n[[ looking into original face = %d]]\n", face);
        xf_printf("   [ iiface = %d, faceL = %d, faceR = %d, HangNumber=%d]\n", 
                  iiface, IFace.FaceL, IFace.FaceR, IFace.HangNumber);
      }
      
      // egR, eR, faceR, orientR, orientL
      
      // get corresponding face or faces from bfacelist
      for (i=0; i<nbface; i++){
        if (ah_DEBUG) xf_printf("i = %d, face0=%d, bfacelist.Face = %d\n", i, face0, bfacelist[i].Face);
        if (bfacelist[i].Face == face0){
          bpos = bface2pos[i];
          ielem = bfacelist[i].Elem;
          eI    = elem2glob[ielem];
          face0I = bfacelist[i].Face;
          bfacelist[i].Face *= -1; // so we don't enter here again
          
          if (ah_DEBUG)
            xf_printf("\n\n**** eI = %d, egrp=%d, face0=%d, face0I = %d, pos=%d, bpos=%d ***\n\n",
                      eI, egrp, face0, face0I, pos, bpos);
          
          // make sure not refining a fine side of a hanging face
          if (bpos != 0)
            if ( ((IamL) && (IFace.HangNumber>0)) || ((!IamL) && (IFace.HangNumber<0)) ){
              xf_printf("Error, trying to refine the fine-elem side of a hanging face.\n");
              // should have caught this in consistency/priority checks
              return xf_Error(xf_CODE_LOGIC_ERROR);
            }
          
          // obtain relation between the positions
          ierr = xf_Error(xf_PosCompare(FShape, pos, bpos, &rel, &compatible));
          if (ierr != xf_OK) return ierr;
          
          // check if the proposed position (bpos) is compatible with the existing position (pos)
          if (!compatible){
            xf_printf("Error in adapting egrp=%d, elem=%d\n", egrp, elem);
            xf_printf("bpos = %d, is not compatible with existing pos = %d\n", bpos, pos);
            return xf_Error(xf_CODE_LOGIC_ERROR);
          }
          
          
          if ((pos == 0) && (bpos == 0)){ // face0 remains untouched
            if (IamL){
              if (face0I != IFace.FaceL) return xf_Error(xf_MESH_ERROR);
              Mesh->IFace[iiface].ElemGroupL = egrp;
              Mesh->IFace[iiface].ElemL = eI;
              Mesh->IFace[iiface].FaceL = face0I;
            }
            else{
              if (face0I != IFace.FaceR) return xf_Error(xf_MESH_ERROR);
              Mesh->IFace[iiface].ElemGroupR = egrp;
              Mesh->IFace[iiface].ElemR = eI;
              Mesh->IFace[iiface].FaceR = face0I;
            }
            iifacenew = iiface;
            Mesh->ElemGroup[egrp].Face[eI][face0I].Group  = xf_INTERIORFACE;
            Mesh->ElemGroup[egrp].Face[eI][face0I].Number = iifacenew;
          }
          else if (bpos == 0){ // face0 remains a hanging face
            // determine hanging face number in new element
            if (Mesh->ElemGroup[egrp].Face[eI][face0I].Number < 0)
              faceI = face0I;  // indicates face0I has not yet been filled
            else
              faceI = Mesh->ElemGroup[egrp].nFace[eI]++; // use next face
            
            if (IamL){
              Mesh->IFace[iiface].ElemGroupL = egrp;
              Mesh->IFace[iiface].ElemL = eI;
              Mesh->IFace[iiface].FaceL = faceI;
            }
            else{
              Mesh->IFace[iiface].ElemGroupR = egrp;
              Mesh->IFace[iiface].ElemR = eI;
              Mesh->IFace[iiface].FaceR = faceI;
            }
            
            iifacenew = iiface;
            Mesh->ElemGroup[egrp].Face[eI][faceI].Group  = xf_INTERIORFACE;
            Mesh->ElemGroup[egrp].Face[eI][faceI].Number = iifacenew;
            if (ah_DEBUG)
              xf_printf("Want to enter here again 1, hangnumber = %d\n", Mesh->IFace[iiface].HangNumber);
            bfacelist[i].Face *= -1; // so we do enter here again
          }
          else if (pos == 0){ // face0 becomes a hanging face
            if (IamL){
              orientI = IFace.OrientL;
              face0C  = IFace.FaceR;
              egC = IFace.ElemGroupR;
              eC  = IFace.ElemR;
            }
            else{
              orientI = IFace.OrientR;
              face0C  = IFace.FaceL;
              egC = IFace.ElemGroupL;
              eC  = IFace.ElemL;
            }
            
            if (ah_DEBUG)
              xf_printf(" pos = 0, but face0C = %d, IamL = %d, hang = %d, eC = %d, eI = %d, (orientL=%d, orientR=%d)\n", 
                        face0C, IamL, IFace.HangNumber, eC, eI, IFace.OrientL, IFace.OrientR);
            
            ierr = xf_Error(xf_HangingFaceAdd(Mesh, IFace, IamL, egrp, eI, face0I, 
                                              orientI, bpos, &iifacenew));
            if (ierr != xf_OK) return ierr;
            
            // so we re-orient this face (do not need to ... will reorient in DetermineFaceOrient
            FaceToOrient[nFaceToOrient][0] = IamL;
            FaceToOrient[nFaceToOrient][1] = iifacenew;
            if (nFaceToOrient++ >= nfacemax) return xf_Error(xf_CODE_LOGIC_ERROR);
            
            ToOrient[nToOrient][0] = egC;
            ToOrient[nToOrient][1] = eC;
            ToOrient[nToOrient][2] = face0C;
            nToOrient++;
            
            if (nToOrient > nToOrientMax) return xf_Error(xf_CODE_LOGIC_ERROR);
          }
          else if (pos == bpos){ // face0 is no longer a hanging face
            if (IamL){
              faceI = face0I;
              Mesh->IFace[iiface].ElemGroupL = egrp;
              Mesh->IFace[iiface].ElemL = eI;
              Mesh->IFace[iiface].FaceL = faceI;
              orientC = Mesh->IFace[iiface].OrientL;
              orientF = Mesh->IFace[iiface].OrientR;
              QBasisC = Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupL].QBasis;
              QBasisF = Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupR].QBasis;
              QOrderC = Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupL].QOrder;
              QOrderF = Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupR].QOrder;
              faceC = face0I;
              faceF = Mesh->IFace[iiface].FaceR;
              egF = Mesh->IFace[iiface].ElemGroupR;
              eF  = Mesh->IFace[iiface].ElemR;
            }
            else{
              faceI = face0I;
              Mesh->IFace[iiface].ElemGroupR = egrp;
              Mesh->IFace[iiface].ElemR = eI;
              Mesh->IFace[iiface].FaceR = faceI;
              orientF = Mesh->IFace[iiface].OrientL;
              orientC = Mesh->IFace[iiface].OrientR;
              QBasisF = Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupL].QBasis;
              QBasisC = Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupR].QBasis;
              QOrderF = Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupL].QOrder;
              QOrderC = Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupR].QOrder;
              faceF = Mesh->IFace[iiface].FaceL;
              faceC = face0I;
              egF = Mesh->IFace[iiface].ElemGroupL;
              eF  = Mesh->IFace[iiface].ElemL;
            }
            if (ah_DEBUG)
              xf_printf("no longer hanging: face0=%d, eF=%d, faceF=%d, eI = %d, faceC=%d\n",
                        face0, eF, faceF, eI, faceC);
            if (ah_DEBUG)
              xf_printf(" orientF =%d, orientC = %d\n",
                        orientF, orientC);
            
            Mesh->IFace[iiface].HangNumber = 0;
            iifacenew = iiface;
            Mesh->ElemGroup[egrp].Face[eI][faceI].Group  = xf_INTERIORFACE;
            Mesh->ElemGroup[egrp].Face[eI][faceI].Number = iifacenew;
            
            // so we re-orient this face
            FaceToOrient[nFaceToOrient][0] = IamL;
            FaceToOrient[nFaceToOrient][1] = iifacenew;
            if (nFaceToOrient++ >= nfacemax) return xf_Error(xf_CODE_LOGIC_ERROR);
            
            // get verts on face; do this by calling q1nodesonface with q=1
            ierr = xf_Error(xf_Q1NodesOnFace(QBasisC, 1, faceC, &nn, vvec));
            if (ierr != xf_OK) return ierr;
            // use orientF, orientC, etc. to get globnodes on face
            ierr = xf_Error(xf_Q1NodesOnFaceNeighbor(QBasisF, QOrderF, faceF, orientF, 
                                                     QBasisC, QOrderC, faceC, orientC, 
                                                     &nn1, nvec));
            if (ierr != xf_OK) return ierr;
            if (ah_DEBUG)
              for (k=0; k<nn; k++)
                xf_printf("Q1onFN[%d] = %d\n", k, nvec[k]);
            
            if (nn != nn1) return xf_Error(xf_CODE_LOGIC_ERROR);
            // set new vert2glob and check consistency of those already set
            for (k=0; k<nn; k++){
              ivert = elem2vert[nNodeQ1*ielem+vvec[k]];
              if (ah_DEBUG)
                xf_printf("egF = %d, eF = %d, k = %d, nvec[k] = %d\n", egF, eF, k, nvec[k]);
              iglob = Mesh->ElemGroup[egF].Node[eF][nvec[k]];
              if (ah_DEBUG)
                xf_printf("ivert = %d, iglob = %d, vert2glob = %d\n", ivert, iglob, vert2glob[ivert]);
              if (vert2glob[ivert] == -1)
                vert2glob[ivert] = iglob;
              else if (vert2glob[ivert] != iglob)
                return xf_Error(xf_CODE_LOGIC_ERROR); 
            }
          }
          else if (rel == xfe_RelSuperset){ // pos contains bpos: face becomes a hanging face
            if (IamL){
              orientI = IFace.OrientL;
              face0C  = IFace.FaceR;
              egC = IFace.ElemGroupR;
              eC  = IFace.ElemR;
              orientN = IFace.OrientR;
              QBasisI = Mesh->ElemGroup[IFace.ElemGroupL].QBasis;
              QBasisN = Mesh->ElemGroup[IFace.ElemGroupR].QBasis;
              QOrderI = Mesh->ElemGroup[IFace.ElemGroupL].QOrder;
              QOrderN = Mesh->ElemGroup[IFace.ElemGroupR].QOrder;
              faceN = Mesh->IFace[iiface].FaceR;
            }
            else{
              orientI = IFace.OrientR;
              face0C  = IFace.FaceL;
              egC = IFace.ElemGroupL;
              eC  = IFace.ElemL;
              orientN = IFace.OrientL;
              QBasisI = Mesh->ElemGroup[IFace.ElemGroupR].QBasis;
              QBasisN = Mesh->ElemGroup[IFace.ElemGroupL].QBasis;
              QOrderI = Mesh->ElemGroup[IFace.ElemGroupR].QOrder;
              QOrderN = Mesh->ElemGroup[IFace.ElemGroupL].QOrder;
              faceN = Mesh->IFace[iiface].FaceL;
            }
            
            // determine relative position of bpos within pos
            ierr = xf_Error(xf_RelativePos(FShape, pos, bpos, &relpos));
            if (ierr != xf_OK) return ierr;
            
            if (ah_DEBUG)
              xf_printf(" pos = 0, face0C = %d, IamL = %d, hangnumber = %d, eC = %d, eI = %d\n", 
                        face0C, IamL, IFace.HangNumber, eC, eI);
            
            ierr = xf_Error(xf_HangingFaceAdd(Mesh, IFace, IamL, egrp, eI, face0I, 
                                              orientI, relpos, &iifacenew));
            if (ierr != xf_OK) return ierr;
            
            // so we re-orient this face
            FaceToOrient[nFaceToOrient][0] = IamL;
            FaceToOrient[nFaceToOrient][1] = iifacenew;
            if (nFaceToOrient++ >= nfacemax) return xf_Error(xf_CODE_LOGIC_ERROR);
            
            // mark coarse side for reorientation
            ToOrient[nToOrient][0] = egC;
            ToOrient[nToOrient][1] = eC;
            ToOrient[nToOrient][2] = face0C;
            nToOrient++;
            if (nToOrient > nToOrientMax) return xf_Error(xf_CODE_LOGIC_ERROR);
            
            // set vert2glob; first obtain flag of common nodes between pos/bpos
            ierr = xf_Error(xf_GetCommonNodeFlag(FShape, relpos, CommonNode));
            if (ierr != xf_OK) return ierr;
            
            // get verts on face; do this by calling q1nodesonface with q=1
            ierr = xf_Error(xf_Q1NodesOnFace(QBasisI, 1, face0I, &nn, vvec));
            if (ierr != xf_OK) return ierr;
            // use orientF, orientC, etc. to get globnodes on face
            ierr = xf_Error(xf_Q1NodesOnFaceNeighbor(QBasisN, QOrderN, faceN, orientN, 
                                                     QBasisI, QOrderI, face0I, orientI, 
                                                     &nn1, nvec));
            if (ierr != xf_OK) return ierr;
            if (nn != nn1) return xf_Error(xf_CODE_LOGIC_ERROR);
            // set new vert2glob and check consistency of those already set
            for (k=0; k<nn; k++){
              ivert = elem2vert[nNodeQ1*ielem+vvec[k]];
              if (!CommonNode[k]) continue;
              iglob = Mesh->ElemGroup[egC].Node[eC][nvec[k]];
              if (vert2glob[ivert] == -1)
                vert2glob[ivert] = iglob;
              else if (vert2glob[ivert] != iglob)
                return xf_Error(xf_CODE_LOGIC_ERROR);
            }
            
          }
          else if (rel == xfe_RelSubset){ // bpos contains pos: face remains hanging, but downgraded
            // determine hanging face number in new element
            if (Mesh->ElemGroup[egrp].Face[eI][face0I].Number < 0)
              faceI = face0I;  // indicates face0I has not yet been filled
            else
              faceI = Mesh->ElemGroup[egrp].nFace[eI]++; // use next face
            
            
            if (IamL){
              Mesh->IFace[iiface].ElemGroupL = egrp;
              Mesh->IFace[iiface].ElemL = eI;
              Mesh->IFace[iiface].FaceL = faceI;
              orientI = IFace.OrientL;
              orientN = IFace.OrientR;
              QBasisI = Mesh->ElemGroup[IFace.ElemGroupL].QBasis;
              QBasisN = Mesh->ElemGroup[IFace.ElemGroupR].QBasis;
              QOrderI = Mesh->ElemGroup[IFace.ElemGroupL].QOrder;
              QOrderN = Mesh->ElemGroup[IFace.ElemGroupR].QOrder;
              faceN = Mesh->IFace[iiface].FaceR;
              egF = Mesh->IFace[iiface].ElemGroupR;
              eF  = Mesh->IFace[iiface].ElemR;
              egC = Mesh->IFace[iiface].ElemGroupL;
              eC  = Mesh->IFace[iiface].ElemL;
            }
            else{
              Mesh->IFace[iiface].ElemGroupR = egrp;
              Mesh->IFace[iiface].ElemR = eI;
              Mesh->IFace[iiface].FaceR = faceI;
              orientI = IFace.OrientR;
              orientN = IFace.OrientL;
              QBasisI = Mesh->ElemGroup[IFace.ElemGroupR].QBasis;
              QBasisN = Mesh->ElemGroup[IFace.ElemGroupL].QBasis;
              QOrderI = Mesh->ElemGroup[IFace.ElemGroupR].QOrder;
              QOrderN = Mesh->ElemGroup[IFace.ElemGroupL].QOrder;
              faceN = Mesh->IFace[iiface].FaceL;
              egF = Mesh->IFace[iiface].ElemGroupL;
              eF  = Mesh->IFace[iiface].ElemL;
              egC = Mesh->IFace[iiface].ElemGroupR;
              eC  = Mesh->IFace[iiface].ElemR;
            }
            
            // determine relative position of pos within bpos
            ierr = xf_Error(xf_RelativePos(FShape, bpos, pos, &relpos));
            if (ierr != xf_OK) return ierr;
            
            // make sure face is hanging already
            if (hang == 0) return xf_Error(xf_CODE_LOGIC_ERROR);
            
            // get new hanging number
            ierr = xf_Error(xf_Pos2Hang(FShape, nface0, face0, relpos, &hangnew));
            if (ierr != xf_OK) return ierr;
            
            // make appropriate sign
            if (IFace.HangNumber < 0) hangnew *= -1;
            
            iifacenew = iiface;
            Mesh->ElemGroup[egrp].Face[eI][faceI].Group  = xf_INTERIORFACE;
            Mesh->ElemGroup[egrp].Face[eI][faceI].Number = iifacenew;
            Mesh->IFace[iifacenew].HangNumber = hangnew;
            
            // mark coarse side for reorientation
            ToOrient[nToOrient][0] = egC;
            ToOrient[nToOrient][1] = eC;
            ToOrient[nToOrient][2] = face0I;
            nToOrient++;
            if (nToOrient > nToOrientMax) return xf_Error(xf_CODE_LOGIC_ERROR);
            
            // set vert2glob; first obtain flag of common nodes between pos/bpos
            ierr = xf_Error(xf_GetCommonNodeFlag(FShape, relpos, CommonNode));
            if (ierr != xf_OK) return ierr;
            
            // get verts on face; do this by calling q1nodesonface with q=1
            ierr = xf_Error(xf_Q1NodesOnFace(QBasisI, 1, face0I, &nn, vvec));
            if (ierr != xf_OK) return ierr;
            // use orientF, orientC, etc. to get globnodes on face
            ierr = xf_Error(xf_Q1NodesOnFaceNeighbor(QBasisN, QOrderN, faceN, orientN, 
                                                     QBasisI, QOrderI, face0I, orientI, 
                                                     &nn1, nvec));
            if (ierr != xf_OK) return ierr;
            if (nn != nn1) return xf_Error(xf_CODE_LOGIC_ERROR);
            if (ah_DEBUG)
              for (k=0; k<nn; k++) xf_printf("iglob[%d]=%d\n", k,  Mesh->ElemGroup[egF].Node[eF][nvec[k]]);
            // set new vert2glob and check consistency of those already set
            for (k=0; k<nn; k++){
              ivert = elem2vert[nNodeQ1*ielem+vvec[k]];
              if (!CommonNode[k]) continue;
              iglob = Mesh->ElemGroup[egF].Node[eF][nvec[k]];
              if (ah_DEBUG)
                xf_printf("k=%d, vert2glob[%d] = %d, iglob=%d\n", k, ivert, vert2glob[ivert], iglob);
              if (vert2glob[ivert] == -1)
                vert2glob[ivert] = iglob;
              else if (vert2glob[ivert] != iglob)
                return xf_Error(xf_CODE_LOGIC_ERROR); 
            }
            
            if (ah_DEBUG)
              xf_printf("Want to enter here again 1.5, hangnumber = %d\n", Mesh->IFace[iiface].HangNumber);
            bfacelist[i].Face *= -1; // so we do enter here again
            
          }
          else{ // at this point, bpos != pos (but they are compatible -- we already checked this)
            if (ah_DEBUG) xf_printf("Want to enter here again 2.\n");
            bfacelist[i].Face *= -1; // so we do enter here again
          }
          
        }
      } // i
      
    } // end else if interior face
    else return xf_Error(xf_NOT_SUPPORTED); // not null, interior, or boundary?
  } // end face loop over original faces
  
  
  // Add any non-assigned vert2globs as new nodes; increment Mesh->nNode
  for (ivert=0; ivert<nvert; ivert++){
    if (vert2glob[ivert] == -1){
      iglob = Mesh->nNode;
      // map coord[ivert] to physical space -> xglob
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &PhiData, xfe_True,
                                      1, coord+dim*ivert, Mesh->Coord[iglob]));
      if (ierr != xf_OK) return ierr;
      vert2glob[ivert] = iglob;
      Mesh->nNode++;
    }
  } // ivert
  
  // Point .Node of each subelement to the approriate nodes (existing or new)
  for (i=0; i<nelem; i++){
    
    ielem = vnew[i]; //((i == 0) ? elem : nelem0+i-1);
    
    // initialize .Node to -1 = not set
    for (k=0; k<nNode; k++) Mesh->ElemGroup[egrp].Node[ielem][k] = -1;
    
    // map Q1 nodes of this subelement to glob nodes indicated by verts
    for (k=0; k<nNodeQ1; k++){
      Mesh->ElemGroup[egrp].Node[ielem][nvecQ1[k]] = vert2glob[elem2vert[nNodeQ1*i+k]];
      if (ah_DEBUG) xf_printf("  Pointing elem=%d Q1 node %d to %d\n",
                              ielem, k, vert2glob[elem2vert[nNodeQ1*i+k]]);
    }
    
    // point rest of new element nodes to new global nodes; add these nodes too
    for (k=0; k<nNode; k++)
      if (Mesh->ElemGroup[egrp].Node[ielem][k] == -1){
        Mesh->ElemGroup[egrp].Node[ielem][k] = Mesh->nNode;
        for (d=0; d<dim; d++)
          Mesh->Coord[Mesh->nNode][d] = xglob[i*dim*nNode+k*dim+d];
        Mesh->nNode++;
      }
    
  } //i
  
  
  // add element-interior ifaces (and set orientations)
  if (ah_DEBUG)
    xf_printf("Adding interior faces, niface = %d, Mesh->nIFace = %d\n", niface, Mesh->nIFace);
  for (i=0; i<niface; i++){
    eL = elem2glob[ifacelist[i].ElemL];
    eR = elem2glob[ifacelist[i].ElemR];
    xf_AddInteriorFace(Mesh, egrp, eL, ifacelist[i].FaceL, 
                       egrp, eR, ifacelist[i].FaceR);
    ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NULL, egrp, elem2glob[ifacelist[i].ElemL], 
                                           ifacelist[i].FaceL,
                                           &Mesh->IFace[Mesh->nIFace-1].OrientL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NULL, egrp, elem2glob[ifacelist[i].ElemR], 
                                           ifacelist[i].FaceR,
                                           &Mesh->IFace[Mesh->nIFace-1].OrientR));
    if (ierr != xf_OK) return ierr;
    if (ah_DEBUG)
      xf_printf("Called (I) DetermineFaceOrient on egrp=%d, elem=(%d,%d), face=(%d,%d): orient = (%d,%d)\n",
                egrp, elem2glob[ifacelist[i].ElemL], elem2glob[ifacelist[i].ElemR],
                ifacelist[i].FaceL, ifacelist[i].FaceR, Mesh->IFace[Mesh->nIFace-1].OrientL,
                Mesh->IFace[Mesh->nIFace-1].OrientR);
    
  } // i
  
  // set orientations of any extra interior faces we marked above
  if (ah_DEBUG) xf_printf("nFaceToOrient = %d\n", nFaceToOrient);
  for (i=0; i<nFaceToOrient; i++){
    IamL   = FaceToOrient[i][0];
    iiface = FaceToOrient[i][1];
    if (IamL){
      ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NULL, Mesh->IFace[iiface].ElemGroupL, 
                                             Mesh->IFace[iiface].ElemL, Mesh->IFace[iiface].FaceL,
                                             &Mesh->IFace[iiface].OrientL));
      if (ierr != xf_OK) return ierr;
      if (ah_DEBUG)
        xf_printf("Called (L) DetermineFaceOrient on egrp=%d, elem=%d, face=%d: orient = %d\n",
                  Mesh->IFace[iiface].ElemGroupL, Mesh->IFace[iiface].ElemL, Mesh->IFace[iiface].FaceL,
                  Mesh->IFace[iiface].OrientL);
    }
    else{
      ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NULL, Mesh->IFace[iiface].ElemGroupR, 
                                             Mesh->IFace[iiface].ElemR, Mesh->IFace[iiface].FaceR,
                                             &Mesh->IFace[iiface].OrientR));
      if (ierr != xf_OK) return ierr;
      if (ah_DEBUG)
        xf_printf("Called (R) DetermineFaceOrient on egrp=%d, elem=%d, face=%d: orient = %d\n",
                  Mesh->IFace[iiface].ElemGroupR, Mesh->IFace[iiface].ElemR, Mesh->IFace[iiface].FaceR,
                  Mesh->IFace[iiface].OrientR);
    }
  } // i
  
  
  /* Set orientations of any neighboring coarse hanging faces we created */
  for (i=0; i<nToOrient; i++){
    egC    = ToOrient[i][0];
    eC     = ToOrient[i][1];
    face0C = ToOrient[i][2];
    if (ah_DEBUG) xf_printf(" Calling orient on egC = %d, eC = %d, face0C = %d\n", egC, eC, face0C);
    
    ierr = xf_Error(xf_SetCoarseOrients(Mesh, egC, eC, face0C, xfe_True, NULL, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  // Add new bisected edges while verifying consistency
  for (i=0; i<nedge; i++){
    // two edge nodes
    n0 = vert2glob[elist[3*i+0]];
    n1 = vert2glob[elist[3*i+1]];
    if (n0 > n1) swap(n0, n1, j);
    if (n0 < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
    ierr = xf_EdgeHashCheck(n0, n1, node2ehash, edgehash, &ihash);
    if (ierr == xf_OK){ // edge bisection exists, check glob node number
      if (vert2glob[elist[3*i+2]] != edgehash[ihash].idata)
        return xf_Error(xf_CODE_LOGIC_ERROR);
    }
    else{
      ihash = (*pnedgehash);
      ierr = xf_Error(xf_EdgeHashAdd(n0, n1, node2ehash, edgehash, pnedgehash));
      if (ierr != xf_OK) return ierr;
      if (ihash == (*pnedgehash)) return xf_Error(xf_CODE_LOGIC_ERROR);
      edgehash[ihash].idata = vert2glob[elist[3*i+2]];
    }
  } // i
  
  
  // release memory
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) xref0);
  xf_Release( (void *) xref);
  xf_Release( (void *) xglob);
  xf_Release( (void *) elem2glob);
  xf_Release( (void *) FaceOrig);
  xf_Release2( (void **) ToOrient);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReOrientElem
static int 
xf_ReOrientElem(xf_Mesh *Mesh, xf_IntPairList ElemReOrient)
{
  /*
   This function was written for reorienting elements after node
   numbering changes.  However, it should no longer be called, as node
   renumbering (based on same node pairs) has been eliminated.  The
   code was left in case such cabability is required again.
   */
  int ierr, i;
  int egrp, elem, face, face0;
  int nface, nface0;
  int *pOrient;
  enum xfe_ShapeType Shape;
  enum xfe_Bool IamL;
  xf_IFace *IFace;
  xf_Face Face;
  
  if (ElemReOrient.n <= 0) return xf_Error(xf_INPUT_ERROR);
  
  return xf_Error(xf_NOT_SUPPORTED); // should never get here
  
  // first reorient non-hang element faces
  for (i=0; i<ElemReOrient.n; i++){
    egrp = ElemReOrient.Pairs[2*i+0];
    elem = ElemReOrient.Pairs[2*i+1];
    
    nface = Mesh->ElemGroup[egrp].nFace[elem];
    
    // get element Shape
    ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    
    // number of original faces in Shape
    ierr = xf_Error(xf_Shape2nFace(Shape, &nface0));
    if (ierr != xf_OK) return ierr;
    
    for (face=0; face<nface; face++){
      Face = Mesh->ElemGroup[egrp].Face[elem][face];
      if (Face.Group >= 0) continue; // only care about interior faces
      IFace = Mesh->IFace + Face.Number;
      if (IFace->HangNumber != 0) continue; // only non-hanging faces first
      ierr = xf_Error(xf_IsElemOnLeft(*IFace, egrp, elem, &IamL));
      if (ierr != xf_OK) return ierr;
      pOrient = ((IamL) ?  &IFace->OrientL : &IFace->OrientR);
      ierr = xf_Error(xf_DetermineFaceOrient(Mesh, NULL, egrp, elem, face, pOrient));
      if (ierr != xf_OK) return ierr;
    } // face
  }
  
  // next reorient hanging element faces
  for (i=0; i<ElemReOrient.n; i++){
    egrp = ElemReOrient.Pairs[2*i+0];
    elem = ElemReOrient.Pairs[2*i+1];
    
    nface = Mesh->ElemGroup[egrp].nFace[elem];
    
    // get element Shape
    ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    
    // number of original faces in Shape
    ierr = xf_Error(xf_Shape2nFace(Shape, &nface0));
    if (ierr != xf_OK) return ierr;
    
    for (face=0; face<nface; face++){
      Face = Mesh->ElemGroup[egrp].Face[elem][face];
      if (Face.Group >= 0) continue; // only care about interior faces
      IFace = Mesh->IFace + Face.Number;
      ierr = xf_Error(xf_IsElemOnLeft(*IFace, egrp, elem, &IamL));
      if (ierr != xf_OK) return ierr;
      if ( ( IFace->HangNumber == 0) ||
          ((IFace->HangNumber < 0) && (!IamL)) ||
          ((IFace->HangNumber > 0) && (IamL)) )
        continue; // only hanging faces this time
      // get original face
      ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, NULL, &face0, NULL, NULL));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetCoarseOrients(Mesh, egrp, elem, face0, xfe_True, NULL, NULL));
      if (ierr != xf_OK) return ierr;
    } // face
  }
  
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PruneMesh
static int 
xf_PruneMesh(xf_Mesh *Mesh, int *nElem0, int **nFace0, int nNodeOrig,
             xf_IntPairList SameNode)
{
  /*
   PURPOSE:
   
   Prunes and reallocates structures in Mesh, which must be valid on
   input.  Unused Nodes and Interior/boundary faces are removed and
   connectivity is remapped.
   
   INPUTS:
   
   Mesh : mesh structure
   nElem0 : nElem0[egrp] = original (overallocated) number of elems
   nFace0 : nFace0[egrp][elem] = orig. overallocated # of faces per elem
   nNodeOrig : number of original nodes in the mesh
   SameNode : pairs of identical nodes with different indices (not used)
   
   OUTPUTS: 
   
   None: Mesh is pruned
   
   RETURN:
   
   Error Code
   */
  int ierr, i, j, d;
  int negrp, egrp, elem;
  int nNode0, nNode;
  int ip, isafety;
  enum xfe_Bool done;
  int *nodeflag = NULL;
  xf_IntPairList ElemReOrient;
  
  
  negrp = Mesh->nElemGroup;
  
  ElemReOrient.n  = 0;
  ElemReOrient.n0 = 0;
  ElemReOrient.Pairs = NULL;
  
  
  /*-------*/
  /* Nodes */
  /*-------*/
  
  // nodeflag is used to flag all used nodes
  nNode0 = Mesh->nNode;
  ierr = xf_Error(xf_Alloc((void **) &nodeflag, nNode0, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  if (SameNode.n > 0){
    if (ah_DEBUG) xf_printf("SameNode.n = %d\n", SameNode.n);
    
    return xf_Error(xf_NOT_SUPPORTED);
    
    // initialize nodeflag to identity map
    for (i=0; i<nNode0; i++) nodeflag[i] = i;
    
    
    // Use SameNodePair to modify nodeflag
    done = xfe_False;
    isafety = 0;
    while (!done){
      done = xfe_True;
      for (ip=0; ip<SameNode.n; ip++){
        j = min(nodeflag[SameNode.Pairs[2*ip+0]], nodeflag[SameNode.Pairs[2*ip+1]]);
        if ((j != nodeflag[SameNode.Pairs[2*ip+0]]) || 
            (j != nodeflag[SameNode.Pairs[2*ip+1]])){
          done = xfe_False;
          for (i=0; i<2; i++) nodeflag[SameNode.Pairs[2*ip+i]] = j;
        }
      }
      if (isafety++ > 2*SameNode.n) return xf_Error(xf_INFINITE_WHILE_ERROR);
    }
    
    // point to new nodes
    for (egrp=0; egrp<negrp; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
        for (i=0; i<Mesh->ElemGroup[egrp].nNode; i++){
          j = nodeflag[Mesh->ElemGroup[egrp].Node[elem][i]];
          if (j != Mesh->ElemGroup[egrp].Node[elem][i]){
            Mesh->ElemGroup[egrp].Node[elem][i] = j;
            ierr = xf_Error(xf_AddIntPair(egrp, elem, &ElemReOrient));
            if (ierr != xf_OK) return ierr;
          }
        }
  }
  
  // (re)set nodeflag
  for (i=0; i<nNode0; i++) nodeflag[i] = -1;
  
  // mark all used nodes
  for (egrp=0; egrp<negrp; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      for (i=0; i<Mesh->ElemGroup[egrp].nNode; i++)
        nodeflag[Mesh->ElemGroup[egrp].Node[elem][i]] = 1;
  
  // set nodeflag[i] = new node number
  for (i=0, nNode=0; i<nNode0; i++)
    if (nodeflag[i] != -1) nodeflag[i] = nNode++;
  
  if (nNode > Mesh->nNode) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  // copy over Mesh->Coord and reallocate
  for (i=0; i<nNode0; i++)
    if ((nodeflag[i] != -1) && ((j=nodeflag[i]) != i))
      for (d=0; d<Mesh->Dim; d++) Mesh->Coord[j][d] = Mesh->Coord[i][d];
  
  ierr = xf_Error(xf_ReAllocCopy2( (void ***) &Mesh->Coord, nNode0, Mesh->Dim, 
                                  nNode, Mesh->Dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  Mesh->nNode = nNode; // set new number of nodes
  
  // point to new nodes
  for (egrp=0; egrp<negrp; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      for (i=0; i<Mesh->ElemGroup[egrp].nNode; i++)
        Mesh->ElemGroup[egrp].Node[elem][i] = nodeflag[Mesh->ElemGroup[egrp].Node[elem][i]];
  
  
  // ReOrient elems that were affected by node change
  if (ElemReOrient.n > 0){
    ierr = xf_Error(xf_ReOrientElem(Mesh, ElemReOrient));
    if (ierr != xf_OK) return ierr;
  }
  
  
  /*----------*/
  /* Elements */
  /*----------*/
  
  for (egrp=0; egrp<negrp; egrp++){
    
    if (ah_DEBUG)
      xf_printf("egrp %d: nelem = %d, nelem0 = %d\n", egrp, 
                Mesh->ElemGroup[egrp].nElem, nElem0[egrp]);
    
    // not overallocating elements curently, so should match exactly
    if (Mesh->ElemGroup[egrp].nElem != nElem0[egrp])
      return xf_Error(xf_OUT_OF_BOUNDS);
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      if (Mesh->ElemGroup[egrp].nFace[elem] > nFace0[egrp][elem])
        return xf_Error(xf_OUT_OF_BOUNDS);
    
    
    if (Mesh->ElemGroup[egrp].nElem != nElem0[egrp]){
      
      // since not overallocating elements, should not get here
      return xf_Error(xf_CODE_LOGIC_ERROR);
      
      // Mesh->ElemGroup[egrp].Node
      ierr = xf_Error(xf_ReAllocCopy2( (void ***) &Mesh->ElemGroup[egrp].Node,
                                      nElem0[egrp], Mesh->ElemGroup[egrp].nNode,
                                      Mesh->ElemGroup[egrp].nElem, 
                                      Mesh->ElemGroup[egrp].nNode, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      
      // Mesh->ElemGroup[egrp].nFace
      ierr = xf_Error(xf_ReAlloc((void **) &Mesh->ElemGroup[egrp].nFace, 
                                 Mesh->ElemGroup[egrp].nElem, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
    }      
    
    // Mesh->ElemGroup[egrp].Face  
    ierr = xf_Error(xf_VReAllocCopy2( (void ***) &Mesh->ElemGroup[egrp].Face,
                                     nElem0[egrp], nFace0[egrp], Mesh->ElemGroup[egrp].nElem, 
                                     Mesh->ElemGroup[egrp].nFace, sizeof(xf_Face)));
    if (ierr != xf_OK) return ierr;  
    
  } // egrp
  
  
  xf_Release( (void *) nodeflag);
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_MapDataHang
static int 
xf_MapDataHang(xf_All *All, int *nElemOld, xf_Vector *RefIndicator,
               int **OldElem2nNew, int ***OldElem2New, int **NewElem2Pos)
{
  /*
   PURPOSE:
   
   Maps data from old mesh to new mesh and deletes data that it cannot
   map.
   
   Also maps point value outputs in All to new elems
   
   INPUTS:
   
   All : all structure
   nElemOld : nElemOld[egrp] = number of elems per group in old mesh
   RefIndicator : element-based refinement indicator (assumed consistent)
   OldElem2nNew : OldElem2nNew[egrp][elem] = # elems per old elem
   OldElem2New : OldElem2New[egrp][elem][i] = # ith new elem per old elem
   NewElem2Pos : NewElem2Pos[egrp][elem] = pos of each new elem
   
   OUTPUTS: 
   
   None, data is reallocated + mapped or deleted
   
   RETURN:
   
   Error Code
   */
  
  int ierr, sr, r, nn, i, k, ref;
  int pOrder, Order, nelem, nref, npos;
  int newelem, pos, egrp, elem;
  int iOutput, j, dim;
  int *nElem = NULL, rga, *vr = NULL;
  enum xfe_BasisType Basis;
  enum xfe_ShapeType Shape;
  enum xfe_Bool Interpolated;
  enum xfe_Bool VariableOrder;
  enum xfe_Bool deleteflag;
  enum xfe_Bool inside, done;
  real *v0 = NULL, **TT = NULL;
  real *xref, xref0[3];
  xf_DataSet *DataSet;
  xf_Data *D, *N;
  xf_Vector *V;
  xf_Output *Output;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim = Mesh->Dim;
  DataSet = All->DataSet;
  D = DataSet->Head;
  
  while (D != NULL){
    deleteflag = xfe_True;
    if ((D->Type == xfe_Vector) && (D->ReadWrite == xfe_True)){
      V = (xf_Vector *) D->Data;
      Interpolated = ((V->Basis != NULL) && (V->Order != NULL));
      
      /** CASE 1: Interpolated data **/
      if ((V->Linkage == xfe_LinkageGlobElem) && (Interpolated)){
        
        // only map vectors with data
        if ((V->nArray <= 0) || (V->GenArray[0].Size != xfe_SizeReal)) continue;
        
        sr = V->StateRank;
        
        // Are we dealing with a variable-order vector?
        VariableOrder = ((V->nComp != NULL) && (V->vOrder != NULL));
        
        // if variable order, reallocate vOrder and set nComp
        if (VariableOrder){
          
          // get number of elements
          ierr = xf_Error(xf_GetnElem(Mesh, &nElem, NULL));
          if (ierr != xf_OK) return ierr;
          
          // reallocate V->vOrder
          ierr = xf_Error(xf_VReAllocCopy2( (void ***) &V->vOrder, Mesh->nElemGroup,
                                           V->nComp, Mesh->nElemGroup, 
                                           nElem, sizeof(xf_Face)));
          
          // set V->nComp
          for (egrp=0; egrp<Mesh->nElemGroup; egrp++) V->nComp[egrp] = nElem[egrp];
          
          // release nElem
          xf_Release( (void *) nElem);
          
        }
        
        // map data
        for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
          Basis = V->Basis[egrp];
          
          // get element Shape
          ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
          if (ierr != xf_OK) return ierr;
          
          // reallocate V
          nelem = Mesh->ElemGroup[egrp].nElem;
          if (VariableOrder){ // variable order case
            
            // data should have variable rank information
            if (V->GenArray[egrp].vr == NULL) return xf_Error(xf_INPUT_ERROR);
            
            // allocate a new vr
            ierr = xf_Error(xf_Alloc((void **)&vr, nelem, sizeof(int)));
            if (ierr != xf_OK) return ierr;
            
            // fill in the new vr
            for (elem=0; elem<nElemOld[egrp]; elem++){
              vr[elem] = V->GenArray[egrp].vr[elem];
              ref = RefIndicator->GenArray[egrp].iValue[elem][0];
              if (ref == 0) continue; // this elem was not refined
              // loop over subelements
              ierr = xf_Error(xf_Ref2nElem(Shape, ref, &nref));
              if (ierr != xf_OK) return ierr;
              
              for (i=0; i<nref; i++){
                newelem = OldElem2New[egrp][elem][i];
                vr[newelem] = vr[elem]; // same rank on subelements
              } // i
            } // elem
            
            // set V->GenArray[egrp].r = max(vr)
            V->GenArray[egrp].r=0;
            for (j=0; j<V->GenArray[egrp].n; j++)
              V->GenArray[egrp].r = max(V->GenArray[egrp].r, vr[j]);
            
            // reallocate rValue using vr (copy over old data)
            ierr = xf_Error(xf_VReAllocCopy2((void ***) &V->GenArray[egrp].rValue, V->GenArray[egrp].n, 
                                             V->GenArray[egrp].vr, nelem, vr, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            
            // release existing V->GenArray[i].vr
            xf_Release( (void *) V->GenArray[egrp].vr);
            
            // set new vr
            V->GenArray[egrp].vr = vr;
          }
          else{
            ierr = xf_Error(xf_ReAllocCopy2((void ***) &V->GenArray[egrp].rValue, V->GenArray[egrp].n, 
                                            V->GenArray[egrp].r, nelem, V->GenArray[egrp].r, 
                                            sizeof(real)));
            if (ierr != xf_OK) return ierr;
          }
          V->GenArray[egrp].n = nelem; // set new size
          
          pOrder = -1;
          
          // begin loop over old elements
          for (elem=0; elem<nElemOld[egrp]; elem++){
            
            ref = RefIndicator->GenArray[egrp].iValue[elem][0];
            if (ref == 0) continue; // this elem was not refined
            
            // Order
            Order = xf_InterpOrder(V, egrp, elem);
            
            // obtain r = # unknowns * sr per element
            ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
            if (ierr != xf_OK) return ierr;
            r = nn*sr;
            
            // rank as seen by V
            rga = ((V->GenArray[egrp].vr==NULL) ? V->GenArray[egrp].r : V->GenArray[egrp].vr[elem]);
            
            // verify that current .r == r, and that Size == real
            if (rga                    != r           ) return xf_Error(xf_OUT_OF_BOUNDS);
            if (V->GenArray[egrp].Size != xfe_SizeReal) return xf_Error(xf_OUT_OF_BOUNDS);
            
            if (pOrder != Order){
              pOrder = Order;
              // calculate projection matrices for each refinement type
              ierr = xf_Error(xf_Shape2nPos(Shape, &npos));
              if (ierr != xf_OK) return ierr;
              ierr = xf_Error(xf_ReAlloc2((void ***) &TT, npos, nn*nn, sizeof(real)));
              if (ierr != xf_OK) return ierr;
              for (pos=0; pos<npos; pos++){
                ierr = xf_Error(xf_ComputeTransferMatrix(Basis, Order, Basis, Order, 
                                                         Shape, pos, TT[pos]));
                if (ierr != xf_OK) return ierr;
              }
              
              // reallocate v0 vector
              ierr = xf_Error(xf_ReAlloc((void **) &v0, r, sizeof(real)));
              if (ierr != xf_OK) return ierr;
            }
            
            // ** project data in each element **
            
            // save original element data in v0
            for (k=0; k<r; k++) v0[k] = V->GenArray[egrp].rValue[elem][k];
            
            // number of elements introduced in refinement
            ierr = xf_Error(xf_Ref2nElem(Shape, ref, &nref));
            if (ierr != xf_OK) return ierr;
            if (nref != OldElem2nNew[egrp][elem]) return xf_Error(xf_OUT_OF_BOUNDS);
            
            // transfer onto each new element
            for (i=0; i<nref; i++){
              newelem = OldElem2New[egrp][elem][i];
              pos = NewElem2Pos[egrp][newelem];
              xf_MxM_Set(TT[pos], v0, nn, nn, sr, V->GenArray[egrp].rValue[newelem]);
              // set new order if variable
              if (VariableOrder) V->vOrder[egrp][newelem] = V->vOrder[egrp][elem];
              
            } // i
          } // elem
          
        } // egrp
        
        deleteflag = xfe_False;
      }
      /** CASE 2: Non-interpolated data **/
      else if ((V->Linkage == xfe_LinkageGlobElem) && (!Interpolated)&& D->ReadWrite ==  xfe_True){
        
        // only map vectors with data
        if (V->nArray <= 0) continue;
        
        if (V != RefIndicator){ // do not map or destroy RefIndicator
          
          // map data (just transfer same to each elem)
          for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
            
            // get element Shape
            ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
            if (ierr != xf_OK) return ierr;
            
            // obtain r = # unknowns
            r = V->GenArray[egrp].r;
            
            // reallocate V (either real or int)
            nelem = Mesh->ElemGroup[egrp].nElem;
            if (V->GenArray[egrp].Size == xfe_SizeReal)
              ierr = xf_Error(xf_ReAllocCopy2((void ***) &V->GenArray[egrp].rValue, V->GenArray[egrp].n, 
                                              r, nelem, r, sizeof(real)));
            else if (V->GenArray[egrp].Size == xfe_SizeInt)
              ierr = xf_Error(xf_ReAllocCopy2((void ***) &V->GenArray[egrp].iValue, V->GenArray[egrp].n, 
                                              r, nelem, r, sizeof(int)));
            else return xf_Error(xf_OUT_OF_BOUNDS);
            
            if (ierr != xf_OK) return ierr;
            
            V->GenArray[egrp].n = nelem; // set new size
            
            // project data in each element
            for (elem=0; elem<nElemOld[egrp]; elem++){
              ref = RefIndicator->GenArray[egrp].iValue[elem][0];
              if (ref == 0) continue; // this elem was not refined
              
              // number of elements introduced in refinement
              ierr = xf_Error(xf_Ref2nElem(Shape, ref, &nref));
              if (ierr != xf_OK) return ierr;
              if (nref != OldElem2nNew[egrp][elem]) return xf_Error(xf_OUT_OF_BOUNDS);
              
              // transfer onto each new element
              for (i=0; i<nref; i++){
                newelem = OldElem2New[egrp][elem][i];
                if (V->GenArray[egrp].Size == xfe_SizeReal){
                  for (k=0; k<r; k++) 
                    V->GenArray[egrp].rValue[newelem][k] = V->GenArray[egrp].rValue[elem][k];
                }
                else{
                  for (k=0; k<r; k++) 
                    V->GenArray[egrp].iValue[newelem][k] = V->GenArray[egrp].iValue[elem][k];
                }
              } // i
            } // elem
            
          } // egrp
        }
        
        deleteflag = xfe_False;
      }
      
    }
    N = D->Next;
    if (deleteflag) xf_DestroyDataInSet(DataSet, D);
    D = N;
  }
  
  // map point outputs to new elems
  if ((All->EqnSet != NULL) && (All->EqnSet->Outputs != NULL)){
    for (iOutput=0; iOutput<All->EqnSet->Outputs->nOutput; iOutput++){
      Output = All->EqnSet->Outputs->Output + iOutput;
      if (Output->Type == xfe_PointValue){
        egrp = Output->egrp;
        elem = Output->elem;
        xref = Output->xref;
        // is this elem refined?
        ref = RefIndicator->GenArray[egrp].iValue[elem][0];
        if (ref == 0) continue; // this elem was not refined
        
        // number of elements introduced in refinement
        ierr = xf_Error(xf_Ref2nElem(Shape, ref, &nref));
        if (ierr != xf_OK) return ierr;
        if (nref != OldElem2nNew[egrp][elem]) return xf_Error(xf_OUT_OF_BOUNDS);
        
        // get element Shape
        ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
        if (ierr != xf_OK) return ierr;
        
        // check that parent contains point
        ierr = xf_Error(xf_InsideShape(Shape, xref, MEPS, &inside));
        if (ierr != xf_OK) return ierr;
        if (!inside) return xf_Error(xf_INPUT_ERROR);
        
        // point goes into first elem that contains it
        done = xfe_False;
        for (i=0; i<nref; i++){
          newelem = OldElem2New[egrp][elem][i];
          pos = NewElem2Pos[egrp][newelem];
          // copy xref into xref0
          for (j=0; j<dim; j++) xref0[j] = xref[j];
          // transform xref0 to child coord
          ierr = xf_Error(xf_ScaleHangInterpolInv(Shape, pos, 1, xref0));
          if (ierr != xf_OK) return ierr;
          // are we inside the child elem?
          ierr = xf_Error(xf_InsideShape(Shape, xref0, MEPS, &inside));
          if (ierr != xf_OK) return ierr;
          if (inside){
            done = xfe_True;
            Output->elem = newelem;
            for (j=0; j<dim; j++) xref[j] = xref0[j];
            break;
          }
        } // i
        if (!done) return xf_Error(xf_CODE_LOGIC_ERROR);
        xf_Release( (void * ) Output->elemLocal); // will get recalculated
        Output->elemLocal = NULL;
        //xf_printf("Adapted! Output->elem = %d\n", Output->elem);
      } // if PointValue
    } // iOutput
  }
  
  // Release memory
  xf_Release( (void *) v0);
  xf_Release2( (void **) TT);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FillHangEdgeHash
static int 
xf_FillHangEdgeHash(xf_Mesh *Mesh, int nedgemax, int *node2ehash, 
                    xf_EdgeHash *edgehash, int *pnedgehash)
{
  /*
   PURPOSE:
   
   Fills a pre-allocated edge hash based on hanging node faces.  Each
   edge in Mesh that is split as part of a hanging face gets put into
   the hash, with the original Q1 nodes of the parent edge and the
   midpoint Q1 node.  Specific to 3D -- in 2D, there are no edges, only
   1d faces.
   
   INPUTS:
   
   Mesh : mesh structure
   nedgemax : maximum number of edges for the hash (what was preallocated)
   node2hash : link into hash from the nodes (preallocated -- modified)
   edgehash : hash list of split edges (preallocated -- modified)
   (*pnedgehash) : original number of edges in hash
   
   OUTPUTS: 
   
   node2hash : link into hash from the nodes (preallocated -- modified)
   edgehash : hash list of split edges (preallocated -- modified)
   (*pnedgehash) : new number of edges in hash
   
   RETURN:
   
   Error Code
   */
  
  int ierr, i, j;
  int n0, n1;
  int egrp, elem, face;
  int nface0, face0;
  int ihash;
  int nedge;
  int elist[ha_MAX_EDGE*3];
  int vface0[xf_MAXLOCFACE];
  enum xfe_Bool IamL;
  enum xfe_ShapeType Shape;
  xf_Face Face;
  xf_IFace *IFace;
  
  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    // get element Shape
    ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    
    // number of original faces in Shape
    ierr = xf_Error(xf_Shape2nFace(Shape, &nface0));
    if (ierr != xf_OK) return ierr;
    
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // reset flag
      for (i=0; i<nface0; i++) vface0[i] = 0;
      
      // loop over faces
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
        Face = Mesh->ElemGroup[egrp].Face[elem][face];
        if (Face.Group != xf_INTERIORFACE) continue; // only care about interior faces
        IFace = Mesh->IFace + Face.Number;
        ierr = xf_Error(xf_IsElemOnLeft(*IFace, egrp, elem, &IamL));
        if (ierr != xf_OK) return ierr;
        if ( ( IFace->HangNumber == 0) ||
            ((IFace->HangNumber < 0) && (!IamL)) ||
            ((IFace->HangNumber > 0) && (IamL)) )
          continue; // only care about hanging faces
        // get original face
        ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, NULL, &face0, NULL, NULL));
        if (ierr != xf_OK) return ierr;
        vface0[face0] = 1; // want to pull of edges on this face0
      } // face
      
      // pull off split edges and add to hash
      for (face0=0; face0<nface0; face0++)
        if (vface0[face0]){
          ierr = xf_Error(xf_SetCoarseOrients(Mesh, egrp, elem, face0,
                                              xfe_False, &nedge, elist));
          if (ierr != xf_OK) return ierr;
          for (i=0; i<nedge; i++){
            n0 = elist[3*i+0]; // two edge nodes
            n1 = elist[3*i+1];
            if (n0 > n1) swap(n0, n1, j);
            ierr = xf_EdgeHashCheck(n0, n1, node2ehash, edgehash, &ihash);
            if (ierr == xf_OK){
              // make sure data (i.e. middle node) agrees
              if (edgehash[ihash].idata != elist[3*i+2]) return xf_Error(xf_MESH_ERROR);
            }
            else{
              // add to hash
              ihash = (*pnedgehash);
              ierr = xf_Error(xf_EdgeHashAdd(n0, n1, node2ehash, edgehash, pnedgehash));
              if (ierr != xf_OK) return ierr;
              if (ihash == (*pnedgehash)) return xf_Error(xf_CODE_LOGIC_ERROR);
              edgehash[ihash].idata = elist[3*i+2];
            }
          } // i
        }
    } // elem
  } // egrp
  
  if ((*pnedgehash) > nedgemax) return xf_Error(xf_MESH_ERROR); // can try increasing nedgemax
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RefineHangPriority
static int 
xf_RefineHangPriority(xf_All *All, xf_Vector *RefIndicator, int *IFaceRef,
                      xf_Vector *RefPriority, int nPriority)
{
  /*
   PURPOSE:
   
   Refines All->Mesh according to desired element requests in
   RefIndicator and according to the refinement priority in
   RefPriority, starting with priority nPriority and working down.
   
   INPUTS:
   
   All : All structure
   RefIndicator : element-based refinement indicator (assumed consistent)
   IFaceRef : associated face refinement vector
   
   RefPriority : priority number for each element (higher priority
   number elements will get refined first)
   
   nPriority : number of priority levels
   
   OUTPUTS: 
   
   None: Mesh is refined
   
   RETURN:
   
   Error Code
   */
  int ierr, i, iiface, hang, pos, k;
  int face0, elemnew, frefloc;
  int nIFaceNew, nface, nNodeNew;
  int negrp, egrp, elem, face, nelem, ref;
  int ibfgrp, ibface, nbfgrp;
  int iPriority, Priority;
  int nNodeOrig, nedgemax, nedgehash, nelemtot;
  int *nElem, *nElemNew, **OldElem2nNew, ***OldElem2New;
  int **nFaceNew, *nBFaceNew, *OldIFace2nNew;
  int **NewElem2Pos;
  int *node2ehash;
  enum xfe_ShapeType Shape, FShape;
  xf_IntPairList SameNode;
  xf_EdgeHash *edgehash;
  xf_IFace *IFace;
  xf_Face Face;
  xf_BFace BFace;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  negrp = Mesh->nElemGroup;
  nbfgrp = Mesh->nBFaceGroup;
  
  
  if (Mesh->Dim == 3){
    // total number of elements
    ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
    if (ierr != xf_OK) return ierr;
    
    // allocate edge hash: node2ehash + edgehash
    ierr = xf_Error(xf_Alloc( (void **) &node2ehash, Mesh->nNode, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<Mesh->nNode; i++) node2ehash[i] = -1;
    nedgemax = 8*(nelemtot+1);
    ierr = xf_Error(xf_Alloc( (void **) &edgehash, nedgemax, sizeof(xf_EdgeHash)));
    if (ierr != xf_OK) return ierr;
    
    // fill in edge hash using existing hanging faces
    nedgehash = 0;
    ierr = xf_Error(xf_FillHangEdgeHash(Mesh, nedgemax, node2ehash, edgehash, &nedgehash));
    if (ierr != xf_OK) return ierr;
  }
  else{
    node2ehash = NULL;
    edgehash = NULL;
    nedgemax = 0;
    nedgehash = 0;
  }
  
  
  /*------------------------*/
  /* Estimate new Mesh size */
  /*------------------------*/
  
  // OldIFace2nNew[iiface] = # of new interior faces per old interior face
  ierr = xf_Error(xf_Alloc( (void **) &OldIFace2nNew, Mesh->nIFace, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // total number of new interior faces (not including element-interior ones)
  nIFaceNew = 0; 
  
  for (iiface=0; iiface<Mesh->nIFace; iiface++){
    // get L element Shape
    ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[Mesh->IFace[iiface].ElemGroupL].QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    
    egrp = Mesh->IFace[iiface].ElemGroupL;
    elem = Mesh->IFace[iiface].ElemL;
    face = Mesh->IFace[iiface].FaceL;
    
    // get original face
    ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, &face0, &pos, NULL));
    if (ierr != xf_OK) return ierr;
    
    // get Shape of face = Shape of face0
    ierr = xf_Error(xf_FaceShape(Shape, face0, &FShape));
    if (ierr != xf_OK) return ierr;
    
    // number of faces introduced in refinement
    ierr = xf_Error(xf_Ref2nElem(FShape, IFaceRef[iiface], &nface));
    if (ierr != xf_OK) return ierr;
    
    OldIFace2nNew[iiface] = nface;
    nIFaceNew += nface;
  }
  
  
  // number of elements in each group
  ierr = xf_Error(xf_GetnElem(Mesh, &nElem, NULL));
  if (ierr != xf_OK) return ierr;
  
  // OldElem2nNew[egrp][elem] = # new elements per old element
  ierr = xf_Error(xf_VAlloc2( (void ***) &OldElem2nNew, negrp, nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // nElemNew[egrp] = number of new elements in each group
  ierr = xf_Error(xf_Alloc( (void **) &nElemNew, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (egrp=0; egrp<negrp; egrp++){
    // get element Shape
    ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    
    nElemNew[egrp] = 0;
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      ref = RefIndicator->GenArray[egrp].iValue[elem][0];
      
      // number of elements introduced in refinement
      ierr = xf_Error(xf_Ref2nElem(Shape, ref, &nelem));
      if (ierr != xf_OK) return ierr;
      
      OldElem2nNew[egrp][elem] = nelem;
      if (ah_DEBUG) xf_printf("elem=%d, o2nN = %d\n", nelem, nelem);
      nElemNew[egrp] += nelem;
      
      // number of element-interior faces introduced in refinement
      ierr = xf_Error(xf_Ref2nIFace(Shape, ref, &nface));
      if (ierr != xf_OK) return ierr;
      
      nIFaceNew += nface;
    } // elem
  } // egrp
  
  
  // OldElem2New[egrp][elem][k] = elem index of kth subelement for old elem
  ierr = xf_Error(xf_Alloc( (void **) &OldElem2New, negrp, sizeof(int **)));
  if (ierr != xf_OK) return ierr;
  
  for (egrp=0; egrp<negrp; egrp++){
    if (ah_DEBUG) xf_printf("nElemNew[%d] = %d\n", egrp, nElemNew[egrp]); fflush(stdout);
    ierr = xf_Error(xf_VAlloc2( (void ***) OldElem2New+egrp, nElem[egrp], 
                               OldElem2nNew[egrp], sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    elemnew = Mesh->ElemGroup[egrp].nElem; // new elements added at end
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      OldElem2New[egrp][elem][0] = elem; // first sub-element retains number
      for (k=1; k<OldElem2nNew[egrp][elem]; k++)
        OldElem2New[egrp][elem][k] = elemnew++;
    } // elem
    
    // sanity check
    if (elemnew != nElemNew[egrp]) return xf_Error(xf_CODE_LOGIC_ERROR);
    
  } //egrp
  
  
  // NewElem2Pos[egrp][elem] = pos of each new elem (to be filled in during adapt)
  ierr = xf_Error(xf_VAlloc2( (void ***) &NewElem2Pos, negrp, nElemNew, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (egrp=0; egrp<negrp; egrp++)
    for (elem=0; elem<nElemNew[egrp]; elem++) NewElem2Pos[egrp][elem]=0;
  
  
  // nFaceNew[egrp][elem] = # faces in each new element (overestimate)
  ierr = xf_Error(xf_VAlloc2( (void ***) &nFaceNew, negrp, nElemNew, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (egrp=0; egrp<negrp; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      nFaceNew[egrp][elem] = 0;
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
        Face = Mesh->ElemGroup[egrp].Face[elem][face];
        if (Face.Group >= 0) nface = 1;                // boundary
        else if (Face.Group == xf_NULLFACE) nface = 1; // null face
        else nface = OldIFace2nNew[Face.Number];       // interior
        nFaceNew[egrp][elem] += nface;
      }
      for (k=1; k<OldElem2nNew[egrp][elem]; k++)
        nFaceNew[egrp][OldElem2New[egrp][elem][k]] = nFaceNew[egrp][elem];
    }
  } // egrp
  
  
  // nBFaceNew[ibfgrp] = # bfaces in each group in new mesh
  ierr = xf_Error(xf_Alloc( (void **) &nBFaceNew, nbfgrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
    nBFaceNew[ibfgrp] = 0;
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      
      BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
      
      egrp = BFace.ElemGroup;
      elem = BFace.Elem;
      face = BFace.Face;
      
      // element refinement
      ref = RefIndicator->GenArray[egrp].iValue[elem][0];
      
      // get original face
      ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, face, &hang, &face0, &pos, NULL));
      if (ierr != xf_OK) return ierr;
      
      // get the element shape
      ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
      if (ierr != xf_OK) return ierr;
      
      // get effect of ref on face
      ierr = xf_Error(xf_ElemRefEffectOnFace(Shape, ref, face0, &frefloc));
      if (ierr != xf_OK) return ierr;
      
      if (frefloc == 0) // no refinement on boundary face
        nBFaceNew[ibfgrp] += 1;
      else{
        // get Shape of face = Shape of face0
        ierr = xf_Error(xf_FaceShape(Shape, face0, &FShape));
        if (ierr != xf_OK) return ierr;
        
        // number of faces introduced in refinement
        ierr = xf_Error(xf_Ref2nElem(FShape, frefloc, &nface));
        if (ierr != xf_OK) return ierr;
        nBFaceNew[ibfgrp] += nface;
      }
    } // ibface
    if (ah_DEBUG) xf_printf("nBFaceNew[%d] = %d\n", ibfgrp, nBFaceNew[ibfgrp]);
  } // ibfgrp
  
  
  // (over)estimate total number of nodes in new mesh
  nNodeNew = Mesh->nNode;
  for (egrp=0; egrp<negrp; egrp++) 
    nNodeNew += nElemNew[egrp]*Mesh->ElemGroup[egrp].nNode;
  if (nNodeNew < Mesh->nNode){
    xf_printf("Number of nodes decreased during refinement.  Make this a warning if ok.\n");
    return xf_Error(xf_MESH_ERROR);
  }
  
  
  // DEBUG
  if (ah_DEBUG){
    xf_printf("nIFaceNew = %d\n", nIFaceNew);
    xf_printf("nNodeNew = %d\n", nNodeNew);
    xf_printf("nelemNew[0] = %d\n", nElemNew[0]);
  }
  
  
  /*----------------------------*/
  /* Reallocate mesh structures */
  /*----------------------------*/
  
  // Mesh->Coord
  nNodeOrig = Mesh->nNode;
  ierr = xf_Error(xf_ReAllocCopy2( (void ***) &Mesh->Coord, Mesh->nNode, Mesh->Dim,
                                  nNodeNew, Mesh->Dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // Mesh->IFace
  ierr = xf_Error(xf_ReAlloc( (void **) &Mesh->IFace, nIFaceNew, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;
  
  for (iiface=Mesh->nIFace; iiface<nIFaceNew; iiface++) xf_InitIFace(Mesh->IFace+iiface);
  
  // Mesh->BFacegroup[ibfgrp].BFace
  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
    ierr = xf_Error(xf_ReAlloc( (void **) &Mesh->BFaceGroup[ibfgrp].BFace, 
                               nBFaceNew[ibfgrp], sizeof(xf_BFace)));
    if (ierr != xf_OK) return ierr;
    for (ibface=Mesh->BFaceGroup[ibfgrp].nBFace; ibface<nBFaceNew[ibfgrp]; ibface++) 
      xf_InitBFace(Mesh->BFaceGroup[ibfgrp].BFace+ibface);
  }
  
  for (egrp=0; egrp<negrp; egrp++){
    
    // Mesh->ElemGroup[egrp].Node
    ierr = xf_Error(xf_ReAllocCopy2( (void ***) &Mesh->ElemGroup[egrp].Node,
                                    nElem[egrp], Mesh->ElemGroup[egrp].nNode,
                                    nElemNew[egrp], Mesh->ElemGroup[egrp].nNode,
                                    sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    // Mesh->ElemGroup[egrp].nFace
    ierr = xf_Error(xf_ReAlloc((void **) &Mesh->ElemGroup[egrp].nFace, nElemNew[egrp], 
                               sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    
    for (elem=nElem[egrp]; elem<nElemNew[egrp]; elem++)
      Mesh->ElemGroup[egrp].nFace[elem] = 0;
    
    
    // Mesh->ElemGroup[egrp].Face  
    ierr = xf_Error(xf_VReAllocCopy2( (void ***) &Mesh->ElemGroup[egrp].Face,
                                     nElem[egrp], Mesh->ElemGroup[egrp].nFace, 
                                     nElemNew[egrp], nFaceNew[egrp], sizeof(xf_Face)));
    if (ierr != xf_OK) return ierr;  
    
    // set new number of elements
    Mesh->ElemGroup[egrp].nElem = nElemNew[egrp];
    
    if (Mesh->ElemGroup[egrp].CutElemData != NULL)
      return xf_Error(xf_NOT_SUPPORTED);
  } // egrp
  
  
  // prepare a node pair list for duplicated nodes (happens in 3D on edges)
  SameNode.n     = 0;
  SameNode.n0    = 0;
  SameNode.Pairs = NULL;
  
  /*-----------------------------*/
  /* Adapt using priority levels */
  /*-----------------------------*/
  
  for (iPriority=nPriority; iPriority >= 0; iPriority--){
    
    for (egrp=0; egrp<negrp; egrp++){
      
      for (elem=0; elem<nElem[egrp]; elem++){
        
        ref      = RefIndicator->GenArray[egrp].iValue[elem][0];
        Priority = RefPriority->GenArray[egrp].iValue[elem][0];
        
        if ((Priority == iPriority) && (ref != 0)){
          
          ierr = xf_Error(xf_RefineElemHang(Mesh, egrp, elem, ref, OldElem2nNew[egrp][elem], 
                                            OldElem2New[egrp][elem], NewElem2Pos[egrp], 
                                            nNodeOrig, &SameNode, node2ehash, edgehash,
                                            &nedgehash));
          if (ierr != xf_OK) return ierr;
        }
        
      } // elem
      
    } // egrp
    
  } // iPriority
  
  // did we add too many edges to hash?
  if (nedgehash > nedgemax){
    xf_printf("nedgehash = %d, nedgemax = %d\n", nedgehash, nedgemax);
    return xf_Error(xf_CODE_LOGIC_ERROR);
  }
  
  // no longer need edgehash
  xf_Release( (void *) node2ehash);
  xf_Release( (void *) edgehash);
  
  if (ah_DEBUG){
    xf_printf("\nBefore prunning\n");
    xf_printf("Mesh->nIFace = %d\n", Mesh->nIFace);
    xf_printf("Mesh->nNode0 = %d\n", Mesh->nNode);
  }
  
  if (Mesh->nNode > nNodeNew) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // prune and consolidate unused nodes and other over-allocated Mesh structures
  ierr = xf_Error(xf_PruneMesh(Mesh, nElemNew, nFaceNew, nNodeOrig, SameNode));
  if (ierr != xf_OK) return ierr;
  
  if (ah_DEBUG){
    xf_printf("\nAfter prunning\n");
    xf_printf("Mesh->nIFace = %d\n", Mesh->nIFace);
    xf_printf("Mesh->nNode = %d\n", Mesh->nNode);
  }
  
  
  // map data from old to new structures, delete data that cannot map
  ierr = xf_Error(xf_MapDataHang(All, nElem, RefIndicator, OldElem2nNew,
                                 OldElem2New, NewElem2Pos));
  if (ierr != xf_OK) return ierr;
  
  
  // PRINT MESH INFO
  if (ah_DEBUG){
    for (egrp=0; egrp<negrp; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
        xf_printf("Elem (%d,%d), nNode = %d, nFace = %d:\n", egrp, elem,
                  Mesh->ElemGroup[egrp].nNode, Mesh->ElemGroup[egrp].nFace[elem]);
        xf_printf("  Nodes = ");
        for (k=0; k<Mesh->ElemGroup[egrp].nNode; k++) 
          xf_printf("%d ", Mesh->ElemGroup[egrp].Node[elem][k]);
        xf_printf("\n");
        xf_printf("  Faces = ");
        for (k=0; k<Mesh->ElemGroup[egrp].nFace[elem]; k++)
          xf_printf("(%d,%d) ", Mesh->ElemGroup[egrp].Face[elem][k].Group,
                    Mesh->ElemGroup[egrp].Face[elem][k].Number);
        xf_printf("\n");
      }
    for (iiface=0; iiface<Mesh->nIFace; iiface++){
      IFace = Mesh->IFace+iiface;
      xf_printf("iiface (%d): L = (%d,%d,%d), R = (%d,%d,%d)\n", iiface,
                IFace->ElemGroupL, IFace->ElemL, IFace->FaceL,
                IFace->ElemGroupR, IFace->ElemR, IFace->FaceR);
    } 
    for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
      for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
        BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
        xf_printf("ibfgrp=%d, ibface=%d, (%d,%d,%d)\n",ibfgrp, ibface,
                  BFace.ElemGroup, BFace.Elem, BFace.Face);
      }
    }
  }
  
  // Parallelization ... make sure parallel structures are valid
  
  // Release memory
  
  for (egrp=0; egrp<negrp; egrp++)
    xf_Release2( (void **) OldElem2New[egrp]);
  xf_Release(  (void  *) OldElem2New);
  xf_Release2( (void **) OldElem2nNew);
  xf_Release(  (void  *) OldIFace2nNew);
  xf_Release2( (void **) NewElem2Pos);
  
  xf_Release(  (void  *) nElem);
  xf_Release(  (void  *) nElemNew);
  xf_Release2( (void **) nFaceNew);
  xf_Release(  (void  *) nBFaceNew);
  
  xf_Release( (void *) SameNode.Pairs);
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AdaptHang
int 
xf_AdaptHang(xf_All *All, xf_Vector *RefIndicator)
{
  int myRank, nProc;
  int ierr, negrp, egrp, elem, i;
  int nPriority;
  int *IFaceRef, *rvec;
  enum xfe_Bool ParallelFlag;
  xf_Vector *ElemFaceRef, *RefPriority;
  xf_Vector *RefIndicator_Glob;
  xf_DataSet *DataSet_Glob;
  xf_Mesh *Mesh_Glob;
  xf_Mesh *Mesh;
  
  
  // Get myRank
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ParallelFlag = ((nProc > 1) && (All->Mesh->ParallelInfo != NULL));
  
  // Return immediately if no elements flagged for refinement
  for (egrp=0, i=0; egrp<All->Mesh->nElemGroup; egrp++)
    for (elem=0; elem<All->Mesh->ElemGroup[egrp].nElem; elem++)
      i += RefIndicator->GenArray[egrp].iValue[elem][0];
  if (ParallelFlag){
    ierr = xf_Error(xf_MPI_Allreduce(&i, 1, xfe_SizeInt, xfe_MPI_SUM));
    if (ierr != xf_OK) return ierr;
  }
  if (i <= 0){
    xf_printf("Warning, no elements flagged for refinement in AdaptHang input.\n");
    return xf_OK;
  }
  
  
  /* For now, serialize before adapting in parallel */
  
  
  if (ParallelFlag){
    
    // RefIndicator
    if (myRank == 0){
      ierr = xf_Error(xf_CreateVector( &RefIndicator_Glob));
      if (ierr != xf_OK) return ierr;
    }
    
    ierr = xf_Error(xf_UnParallelizeVector(All->Mesh, RefIndicator, RefIndicator_Glob));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0) RefIndicator = RefIndicator_Glob;
    
    // DataSet
    if (myRank == 0){
      ierr = xf_Error(xf_CreateDataSet( &DataSet_Glob));
      if (ierr != xf_OK) return ierr;
    }
    
    ierr = xf_Error(xf_UnParallelizeDataSet( All, All->DataSet, DataSet_Glob));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DestroyDataSet(All->DataSet));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0) All->DataSet = DataSet_Glob;    
    
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
  
  /* end of serialization */
  
  
  if (myRank == 0 || !ParallelFlag){
    
    Mesh = All->Mesh;
    negrp = Mesh->nElemGroup;
    
    // DEBUG
    if (ah_DEBUG)
      for (egrp=0; egrp<negrp; egrp++)
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
          xf_printf("IN: egrp = %d, elem = %d, RefIndicator = %d\n",
                    egrp, elem, RefIndicator->GenArray[egrp].iValue[elem][0]);
    
    // IFaceRef will store refinement level of each face
    ierr = xf_Error(xf_Alloc((void **) &IFaceRef, Mesh->nIFace, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<Mesh->nIFace; i++) IFaceRef[i] = 0; // initialize to zero
    
    
    // ***** NOT SURE WE NEED ELEMFACEREF ... HAVEN'T USED IT YET! ****
    
    // rvec = # original faces per element group
    ierr = xf_Error(xf_Alloc( (void **) &rvec, negrp, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (egrp=0; egrp<negrp; egrp++){
      ierr = xf_Error(xf_Basis2nFace(Mesh->ElemGroup[egrp].QBasis, rvec+egrp));
      if (ierr != xf_OK) return ierr;
    }
    
    // ElemFaceRef will store refinement level on each *original* element faces
    ierr = xf_Error(xf_FindVector(All, "ElemFaceRef", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                  NULL, NULL, NULL, NULL, rvec, xfe_SizeInt, xfe_False,  
                                  xfe_False, NULL, &ElemFaceRef, NULL));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_SetZeroVector(ElemFaceRef));  // initialize ElemFaceRef to zero
    if (ierr != xf_OK) return ierr;
    
    // Augment refinement request to satisfy consistency
    ierr = xf_Error(xf_EnsureRefConsistency(All, RefIndicator, IFaceRef, ElemFaceRef));
    if (ierr != xf_OK) return ierr;
    
    // DEBUG
    if (ah_DEBUG)
      for (egrp=0; egrp<negrp; egrp++)
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
          xf_printf("After EnsureRefConsistency: egrp = %d, elem = %d, RefIndicator = %d\n",
                    egrp, elem, RefIndicator->GenArray[egrp].iValue[elem][0]);
    
    // Create a refinement priority vector (higher number will indicate higher priority for ref)
    ierr = xf_Error(xf_FindVector(All, "RefPriority", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                  NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, ParallelFlag, 
                                  xfe_False, NULL, &RefPriority, NULL));
    if (ierr != xf_OK) return ierr;
    
    // Fill in the refinement priority vector
    ierr = xf_Error(xf_FillRefPriority(All, RefIndicator, IFaceRef, RefPriority, &nPriority));
    if (ierr != xf_OK) return ierr;
    
    // DEBUG
    if (ah_DEBUG)
      for (egrp=0; egrp<negrp; egrp++)
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
          xf_printf("Post-priority: egrp = %d, elem = %d, RefPriority = %d\n",
                    egrp, elem, RefPriority->GenArray[egrp].iValue[elem][0]);
    
    // Adapt all elements based on a priority list
    ierr = xf_Error(xf_RefineHangPriority(All, RefIndicator, IFaceRef, RefPriority, nPriority));
    if (ierr != xf_OK) return ierr;
    
    // Snap points to geometry if have one
    ierr = xf_Error(xf_SnapToGeom(All));
    if (ierr != xf_OK) return ierr;

    // *** CALL PARALLEL REPARTITIONING AT END ***
        
    
    // Release memory
    ierr = xf_Error(xf_DestroyVector(ElemFaceRef, xfe_True));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyVector(RefPriority, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    xf_Release( (void *) IFaceRef);
    xf_Release( (void *) rvec);
    
  }
  
  /* Parallelize after adaptation */
  
  if (ParallelFlag){
    
    // RefIndicator
    if (myRank == 0){
      ierr = xf_Error(xf_DestroyVector(RefIndicator_Glob, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    
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
    if (myRank == 0) DataSet_Glob = All->DataSet;
    
    ierr = xf_Error(xf_CreateDataSet(&All->DataSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ParallelizeDataSet( All, DataSet_Glob, All->DataSet));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0){
      ierr = xf_Error(xf_DestroyDataSet(DataSet_Glob));
      if (ierr != xf_OK) return ierr;
    }
    
    // re-parallelize EqnSet
    ierr = xf_Error(xf_ReParallelizeEqnSet(All->EqnSet));
    if (ierr != xf_OK) return ierr;
    
  }
  
  return xf_OK;
}


#if( UNIT_TEST==1 )
#include "xf_AdaptHang.test.in"
#endif
