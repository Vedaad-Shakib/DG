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

#ifndef _xf_MeshTools_h
#define _xf_MeshTools_h 1

/*
  FILE:  xf_MeshTools.h

  This file contains the headers for functions dealing with the MeshTools structure.

*/

#include "xf_MeshToolsStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_GetnElem
extern int 
xf_GetnElem(xf_Mesh *Mesh, int **nElem, int *nelemtot);
/*
PURPOSE:

  Returns vector of number of elements (allocated in this function) in
  each group, as well as the total number of elements.

INPUTS:

  Mesh  : mesh structure

OUTPUTS: 

  (*nElem) : vector containing # of elems in each group (optional)
  (*nelemtot) : total number of elems

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_EgrpElem2Index
extern int 
xf_EgrpElem2Index(xf_Mesh *Mesh, int egrp0, int elem0, int *pie);
/*
PURPOSE:

  Calculates global (to the processor if in parallel) index of element
  (egrp0, elem0)
  
INPUTS:

  Mesh  : mesh structure
  egrp0, elem0 : element in question

OUTPUTS: 

  (*pie) : global index

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Index2EgrpElem
extern int 
xf_Index2EgrpElem(xf_Mesh *Mesh, int ie0, int *pegrp, int *pelem);
/*
PURPOSE:

  Calculates element group/number corresponding to global index ie0.
  
INPUTS:

  Mesh  : mesh structure
  ie0   : global index

OUTPUTS: 

  (*pegrp), (*pelem) : element group and number corresponding to ie0

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_NeighborAcrossFace
extern int 
xf_NeighborAcrossFace(xf_Mesh *Mesh, int egL, int eL, int faceL, 
		      int *egR, int *eR, int *faceR);
/*
PURPOSE:

  Determines element info of element lying across a local face of the
  input element.

INPUTS:

  Mesh  : mesh structure
  egL   : element group number of input element
  eL    : element number of input element
  faceL : face number of input element, across which to look

OUTPUTS: 

  egR   : element group number of neighbor (-1 if face is on boundary)
  eR    : element number of neighbor (-1 if face is on boundary)
  faceR : face number of neighbor (-1 if face is on boundary) [optional]

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_CommonFace
extern int 
xf_CommonFace(xf_Mesh *Mesh, int egL, int eL, int egR, int eR,
	      int *pfaceL, int *pfaceR);
/*
PURPOSE:

  Determines (first) shared interior face between element L and R
  
INPUTS:

  Mesh  : mesh structure
  egL, eL : left elem
  egR, eR : right elem
  
OUTPUTS: 

  faceL, faceR : local face numbers of shared interior face

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_IsElemOnLeft
extern int 
xf_IsElemOnLeft(xf_IFace IFace, int egrp, int elem, enum xfe_Bool *pIamL);
/*
PURPOSE:

  Checks if egrp,elem is on Left side on IFace.  Returns error if elem
  is not on either side of the face.

INPUTS:

  IFace : interior face structure
  egrp, elem : element in question

OUTPUTS: 

  (*IamL) : true if egrp,elem is on Left side of IFace, false otherwise

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FaceElements
extern  void 
xf_FaceElements(xf_Mesh *Mesh, int ibfgrp, int ibface, 
		int *egrpL, int *elemL, int *faceL, 
		int *egrpR, int *elemR, int *faceR);
/*
PURPOSE:

  Returns element(s) adjacent to a face indexed by ibfgrp and ibface.
  The face is regarded as interior if ibfgrp == -1, in which case
  ibface = iiface.  The face outputs are optional.

INPUTS:

  Mesh : mesh structure
  ibfgrp, ibface : interior or boundary face in question

OUTPUTS: 

  egrpL, elemL, faceL : left element
  egrpR, elemR, faceR : right elem (or -1 if on boundary)
  

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_GetFaceOrient
extern int 
xf_GetFaceOrient(xf_Mesh *Mesh, int egrp, int elem, int face, 
		 int *porient);
/*
PURPOSE:

  Returns orientation of face relative to elem.  Uses stored value in
  IFace or BFace.

INPUTS:

  Mesh  : mesh structure
  egrp, elem, face : element face for which to return the orientation

OUTPUTS: 

  (*porient) : returned orientation

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_DetermineFaceOrient
extern int 
xf_DetermineFaceOrient(xf_Mesh *Mesh, int *NodeMap, int egrp, int elem, 
		       int face, int *pfaceorient);
/*
PURPOSE:

  Determines orientation of face relative to elem, using node ordering.

INPUTS:

  Mesh  : mesh structure
  NodeMap : optional node map to use for a different global node numbering
            (e.g. useful when have periodic faces)
  egrp, elem, face : element face for which to calculate the orientation

OUTPUTS: 

  (*pfaceorient) : returned orientation

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CheckHangFace
extern int 
xf_CheckHangFace(xf_Mesh *Mesh, int egrp, int elem, int face, 
		 int *phang, int *pface0, int *pfacepos, int *pelempos);
/*
PURPOSE:

  Checks if the face in question is a hanging face, and if (egrp,elem)
  is the coarse element for that face.  If so, a nonzero hang value is
  returned, along with the orginal face number (if requested) and the
  position index (if requested).

INPUTS:

  Mesh : mesh structure
  egrp, elem : element for which to calculate
  face : face of elem on which xface is defined
  
OUTPUTS: 

  (*phang)  : nonzero if face is hanging and elem is on coarse side (optional)
  (*pface0) : original face number (optional)
  (*pfacepos) : position of hanging face on the original face (optional)
  (*pelempos) : a possible element position of a mirror subelement (optional)
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_SetCoarseOrients
extern int 
xf_SetCoarseOrients(xf_Mesh *Mesh, int egrp, int elem, int face0,
		    enum xfe_Bool SetFlag, int *pnedge, int *elist);
/*
PURPOSE:

  Calculates the coarse-side hanging face orientations.  (egrp,elem)
  must be the coarse element on the side of a hanging face, which was
  face0 originally.

INPUTS:

  Mesh : mesh structure
  egrp, elem : element for which to calculate orientations
  face0 : original face of elem which should now be a hanging face
  SetFlag : True to actually set orientations
  
OUTPUTS: 

  (*pnedge) : number of bisected edges (optional)
  elist : list of three-tuples corresponding to global node numbers of
          bisected edges, e.g. [n0 n1 nm ... ], where nm is the global 
	  node number of the node between n0 and n1.
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_UpdateFaceOrient
extern int 
xf_UpdateFaceOrient(xf_Mesh *Mesh);
/*
PURPOSE:

  Updates the IFace and BFace orientation info w.r.t elements. The
  orientation is stored for each face/elem combination, and indicates
  how the ref space of the face is positioned w.r.t. the element ref
  space.  Presently, the face ref space is based on a global node
  numbering of the face nodes, while the element likes to think of the
  nodes on each face as ordered so that the right-hand rule normal
  points outward.  The faceorient index enables switching between
  these two points of view.

INPUTS:

  Mesh  : mesh structure

OUTPUTS: 

  None : IFace.OrientL/R and BFace.Orient are modified for each face

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetRelativeOrient
extern int 
xf_GetRelativeOrient(enum xfe_ShapeType FShape, int orient1, int orient2, 
		     int *orientrel);
/*
PURPOSE:

  Returns relative orientation orient2 w.r.t orient1.

  Suppose elem1 and elem2 share a common face.  elem1 sees it via
  orient1, and elem2 sees it via orient2.  This function calculates
  orientrel, which is how elem2 sees elem1's face space.

  Put another way, let O1 be a mapping of a coordinate from the
  face-specific space to as-seen by elem1 -- this is described by
  orient1.  Similarly for elem2: O2 corresponds to orient2.  Then
  Orel, corresponding to orientrel is given by:

      Orel = O2 * inv(O1)

INPUTS:

  FShape  : shape of face
  orient1 : orientation of face from vantage point 1
  orient2 : orientation of face from vantage point 2


OUTPUTS: 

  orientrel : orientation of vantage point 1 from vantage point 2

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_DestroyJacobianData
extern int 
xf_DestroyJacobianData( xf_JacobianData *JData);
/*
PURPOSE:

  Destroys a JacobianData structure

INPUTS:

  JacobianData: structure to destroy
  
OUTPUTS: 

  None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ElemJacobian
extern int 
xf_ElemJacobian(xf_Mesh *Mesh, int egrp, int elem, int nq, real *xq,
		unsigned int AllocFlag, enum xfe_Bool PointsChanged,
		xf_JacobianData **pJData);
/*
PURPOSE:

  Calculates geometry jacobian for an element at a specified number of
  points.  If the element is geometrically linear, only one Jacobian
  calculation is performed.  The data is returned in JData, and
  requests for detJ (determinant of J), J (the Jacobian), iJ (inverse
  Jacobian) are controlled with AllocFlag -- see BasisStruct.h for
  details.

INPUTS:

  Mesh : mesh structure
  egrp, elem : element for which to calculate
  nq : number of points at which to calculate
  xq : unrolled point coords
  AllocFlag : 3 bit flag describing which of detJ, J, iJ to calculate
  PointsChanged : True if points have changed since last call
  
OUTPUTS: 

  JacobianData: structure containing jacobian info
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ScaleHangInterpol
extern int 
xf_ScaleHangInterpol(enum xfe_ShapeType Shape, int pos, 
		     int nq, real *xelem);
/*
PURPOSE:

  Scales element ref coords from a fine hanging-node refinement
  element, to the element ref coords on the parent element.


INPUTS:

  Shape : element shape
  pos : position of fine element within parent
  nq : number of points at which to calculate
  xelem : unrolled ref elem coords in fine element
  
OUTPUTS: 

  xelem : unrolled ref elem coords on parent coarse elem
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ScaleHangInterpolInv
extern int 
xf_ScaleHangInterpolInv(enum xfe_ShapeType Shape, int pos, 
			int nq, real *xelem);
/*
PURPOSE:

  Inveres of ScaleHangInterpol.  Scales element ref coords from a
  parent of a hanging-node refinement element,to the element ref
  coords on the child element.

INPUTS:

  Shape : element shape
  pos : position of fine element within parent
  nq : number of points at which to calculate
  xelem : unrolled ref elem coords in parent (coarse) element
  
OUTPUTS: 

  xelem : unrolled ref elem coords on child fine elem
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_Pos2Hang
extern int 
xf_Pos2Hang(enum xfe_ShapeType FShape, int nface, int face0, int pos, 
	    int *phang);
/*
PURPOSE:

  Converts nface,face0,pos -> hang, where hang is a hanging face
  number index.

INPUTS:

  FShape : face shape
  nface : number of faces for the element in question
  face0 : original face number in the element
  pos : position of hanging face on the original face
  
OUTPUTS: 

  (*phang) : index identifying both the original face and the hanging position

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_RefFace2Interpol
extern int 
xf_RefFace2Interpol(xf_Mesh *Mesh, int egrp, int elem, int face, int faceorient,
		    int nq, real *xface, real *xelem);
/*
PURPOSE:

  Converts a set of nq points from face reference space (xface) to element
  interpolation space (xelem).

INPUTS:

  Mesh : mesh structure
  egrp, elem : element for which to calculate
  face : face of elem on which xface is defined
  faceorient : index describing the face's orientation within the 
               reference element.  
  nq : number of points at which to calculate
  xface : unrolled point coords on face
  
OUTPUTS: 

  xelem : unrolled point coords on face; note that at each point,
          dim(xelem) = 1+dim(xface)
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_Ref2GlobElem
extern int 
xf_Ref2GlobElem(xf_Mesh *Mesh, int egrp, int elem, xf_BasisData **pPhiData, 
		enum xfe_Bool PointsChanged, int npoint, real *xref, real *xglob);
/*
PURPOSE:

  Converts a set of npoint points from elem reference space (xref) to
  global space (xglob).

INPUTS:

  Mesh : mesh structure
  egrp, elem : element for which to calculate
  pPhiData : pointer to a BasisData structure (optional).  Use for recurring
             calls for speedup.
  PointsChanged : only if passing pPhiData; force recalculation.
  npoint : number of points
  xref : ref elem coords
  
OUTPUTS: 

  xglob : global coords
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Ref2GlobFaceUsingTable
extern int 
xf_Ref2GlobFaceUsingTable(xf_Mesh *Mesh, int egrp, int elem, int face,
			  int faceorient, enum xfe_Bool PointsChanged, 
			  xf_BasisData **pBasisData, xf_BasisTable *BasisTable,
			  int npoint, real *xref, real *xglob);
/*
PURPOSE:

  Wrapper for Ref2GlobElem that does not recalculate the basis
  (usually the expensive step) for repeated calls with points on
  faces.  The basis evaluations for a given set of points are stored
  for all local face numbers and orientations via a Table.  Note that
  the points passed in, xref, are in the element reference space.

INPUTS:

  Mesh : mesh structure
  egrp, elem, face : element and face for which to calculate
  faceorient : face orientation
  PointsChanged : only if passing pPhiData; force recalculation.
  BasisTable : table where basis evaluations will be stored
  npoint : number of points
  xref : ref elem coords
  
OUTPUTS: 

  (*pBasisData) : structure pointing to the evaluated basis data 
  BasisTable: modified if a new Shape, face, or faceorient was calculated.
  xglob : global coords
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_Glob2RefElem
extern int 
xf_Glob2RefElem(xf_Mesh *Mesh, int egrp, int elem, real *xglob, 
		real tol, enum xfe_Bool Verbose,
		real *xref, enum xfe_Bool *pconverged);
/*
PURPOSE:

  Inverts a possibly-nonlinear reference-to-global map to compute
  element reference coordinates given global coordinates.

  Works for 1 set of coordinates at a time.  The global coordinates
  can be outside of the element.

INPUTS:

  Mesh : mesh structure
  egrp, elem : element for which to calculate
  xglob : global coords
  tol : tolerance for Newton (optional, negative means use built-in)
  Verbose : flag indicating whether or not to be verbose.
  
OUTPUTS: 

  xface : ref elem coords
  (*pconverged) : optional flag; False if inverse mapping not converged
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyNormalData
extern int 
xf_DestroyNormalData( xf_NormalData *NData);
/*
PURPOSE:

  Destroys a NormalData structure

INPUTS:

  NData: structure to destroy
  
OUTPUTS: 

  None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ElemNormal
extern int 
xf_ElemNormal(xf_Mesh *Mesh, int egrp, int elem, int face, 
              int faceorient, int nq, real *xface, 
              xf_NormalData **pNData);
/*
PURPOSE:

  Computes face normals at nq points, xq (face ref space) for element
  egrp,elem.  The normal points outward from elem.  If the face is
  straight, only one normal calculation is performed (regardless of
  nq) and one normal is returned.

INPUTS:

  Mesh : mesh structure
  egrp, elem : group and number of element in question
  face : local face number on element
  faceorient : orientation of face
  nq : number of points at which to calculate
  xface : unrolled point coords on face
  
OUTPUTS: 

  (*pNData) : NormalData structure containing the normals
                
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_IFaceNormal
extern int 
xf_IFaceNormal(xf_Mesh *Mesh, xf_IFace IFace, int nq, real *xq, 
	       xf_NormalData **pNData);
/*
PURPOSE:

  Computes interior face normals at nq points, xq (face ref space).
  The interior normals point from elemL to elemR.  If the face is
  straight, only one normal calculation is performed (regardless of
  nq) and one normal is returned.

INPUTS:

  Mesh : mesh structure
  IFace : interior face on which to compute the normal
  nq : number of points at which to calculate
  xface : unrolled point coords on face
  
OUTPUTS: 

  (*pNData) : NormalData structure containing the normals
                
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_BFaceNormal
extern int 
xf_BFaceNormal(xf_Mesh *Mesh, xf_BFace BFace, int nq, real *xq, 
	       xf_NormalData **pNData, real *nvec);
/*
PURPOSE:

  Computes boundary face normals at nq points, xq (face ref space).
  The interior normals point outward from elem.  If the face is
  straight, only one normal calculation is performed (regardless of
  nq) and one normal is returned.

INPUTS:

  Mesh : mesh structure
  BFace : boundary face on which to compute the normal
  nq : number of points at which to calculate
  xface : unrolled point coords on face
  
OUTPUTS: 

  (*pNData) : NormalData structure containing the normals
  nvec : optional, normal at all quad points.  Memory must
         be already allocated.

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_LinearElemJacobian
extern int 
xf_LinearElemJacobian(enum xfe_BasisType QBasis, int QOrder, int *nvec, 
		      real **Coord, real *detJ);
/*
PURPOSE:

  Computes Jacobian determinant of a linear (Q1) element.

INPUTS:

  QBasis: basis for element interpolation
  QOrder: q-order (should be 1)
  nvec  : nodes corresponding to the element in question (into Coord)
  Coord : coordinate list of all nodes
  
OUTPUTS: 

  (*detJ) : computed Jacobian determinant

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_FindElemGeom
extern int 
xf_FindElemGeom(xf_All *All, xf_Vector **pEG);
/*
PURPOSE:

  Finds/creates an element geometry vector.

INPUTS:

  All   : All structure
  
OUTPUTS: 

  EG  : vector containing element geometry data

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_FindElemHMetric
extern int 
xf_FindElemHMetric(xf_All *All, enum xfe_Bool MakeContinuous, xf_Vector **pEM);
/*
PURPOSE:

  Finds/creates a possibly-continuous element metric (h-size) vector.

INPUTS:

  All   : All structure
  MakeContinuous : if True, returned metric will be continuous
  
OUTPUTS: 

  EM : vector containing element metric data

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ElemSize
extern int 
xf_ElemSize(xf_All *All, int egrp, int elem, xf_Vector *EG, real *h);
/*
PURPOSE:

  Computes hydraulic-diameter element size factor*(Volume/SurfArea),
  where factor depends on shape of element.

INPUTS:

  All   : All structure
  egrp, elem : element in question
  EG  : vector containing element geometry data (from FindElemGeom)
  
OUTPUTS: 

  h : element size (Area/Perimeter)

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_GetRefineCoords
extern int 
xf_GetRefineCoords(enum xfe_ShapeType Shape, int p, int *nnode, 
		   real **pcoord, int *nsplit, int **pvsplit,
		   int *nbound, int **pvbound);
/*
PURPOSE:

  Returns ref coords and a sub-elem list for uniform refinement of
  Shape.

INPUTS:

  Shape : shape to refine
  p     : refinement order (p > 0)
  
OUTPUTS: 

  nnode : number of nodes in refinement
  pcoord : pointer to coord list (is reallocated)
  nsplit : number of subelements
  vsplit : list of subelement nodes (unrolled)
  nbound : number of subelement boundary edges/faces
  vbound : list of subelement boundary edges/faces

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_GetRefineCoordsOnFace
extern int 
xf_GetRefineCoordsOnFace(enum xfe_ShapeType Shape, int face, int p, 
			 int *nnode, real **pcoord, int *nsplit, 
			 int **pvsplit, int *nbound, int **pvbound);
/*
PURPOSE:

  Returns ref coords and a sub-elem list for uniform refinement of
  a face of Shape.

INPUTS:

  Shape : shape to refine
  face  : face of shape to refine
  p     : refinement order (p > 0)
  
OUTPUTS: 

  nnode : number of nodes in refinement
  pcoord : pointer to coord list in elem ref (is reallocated)
  nsplit : number of subfaces
  vsplit : list of subface nodes on face (unrolled)
  nbound : number of subface boundary edges
  vbound : list of subface boundary edges

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_EdgeHashAdd
extern int
xf_EdgeHashAdd(int n0, int n1, int *node2ehash, 
	       xf_EdgeHash *hash, int *nhash);
/*

PURPOSE: Adds an edge (n0-n1) to the hash
  
INPUTS:

  n0, n1: nodes to add
  node2ehash, hash: hash info

OUTPUTS: 

  nhash: resulting number of edges in hash

RETURNS: Error Code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_EdgeHashCheck
extern int
xf_EdgeHashCheck(int n0, int n1, int *node2ehash, 
		 xf_EdgeHash *hash, int *ihash);
/*

PURPOSE:  Checks if edge n0-n1 exists in hash table
  
INPUTS: 
  n0, n1: nodes to check
  node2ehash, hash: hash info

OUTPUTS: 

  ihash: edge number if it exists

RETURNS: xf_NOT_FOUND if edge does not exist

*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetRefineCoordsOnFace
extern int
xf_IntersectElemWithPlane(xf_Mesh *Mesh, int egrp, int elem, int refine,
			  const real *CutPlane, const real *distin, 
			  enum xfe_Bool PointsChanged,  
			  xf_BasisData **pPhiData, int *pnnode, real **pxref, 
			  int *pntri, int **pvtri, int *pnbound, int **pvbound); 
/*
PURPOSE:

  Intersects egrp,elem with a cut-plane (CutPlane).  Returns
  intersections in the form of a set of triangles.  Valid only for 3D
  meshes.

INPUTS:

  Mesh : mesh for egrp,elem
  egrp, elem : element information
  refine : refinement level for elem; if > 1, elem is uniformly refined
           refine times and then intersected with the plane.  The number
	   of intersecting triangles will be larger.  This is useful for
	   plotting.
  CutPlane : 4 real numbers a,b,c,d: ax + by + cz + d = 0
  distin : used if CutPlane == NULL; precomputed distances at nodes to
           be used instead of calculated distances to the plane
  PointsChanged: true if we expect the ref points (pxref) to have changed
                 from a previous call
  pPhiData : basis data for ref-to-glob transformations.
  
OUTPUTS: 

  (*pnnode) : number of nodes in refinement
  (*pxref)  : reference coords of the intersections (3*nnode)
  (*pntri)  : number of triangles resulting from the intersection
  (*pvtri)  : list of triangles (3*ntri)
  (*pnbound): number of boundary edges
  (*pvbound): list of boundary edges

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_CheckVolumes
extern int 
xf_CheckVolumes(xf_Mesh *Mesh, int *OrderVec, enum xfe_Bool ErrorOutFlag,
		int nstoremax, int *pnstore, real *xstore);
/*
PURPOSE:

  Checks Mesh for negative volumes at quad points used to interpolate
  a polynomial of order 1 (default) or one specified by the orders in
  OrderVec.

INPUTS:

  Mesh : mesh structure
  OrderVec : interpolation orders to use for quad points (optional)
  ErrorOutFlag : return with an error if nonpositive volumes are found
  nstoremax : max number of negative vol points to store (glob coords)

OUTPUTS: 

  (*pnstore) : number of points for which storing glob coords (neg vol)
  xstore : locations of negative volume points

RETURN:

  Error Code
*/


/******************************************************************/
//  FUNCTION Prototype: xf_Transfinite2D
extern void 
xf_Transfinite2D(int dim, real *xref, real **x);
/*

PURPOSE: 

  Determines position of interior node given corners and edges

  List of nodes is:

    6 7 8
    3 4 5
    0 1 2

  Node 4 is the one whose position is set.  Its position on input is
  not used.
  
INPUTS:

  dim  : dimension of coordinate vector
  xref : 2d reference coords, xref[0], xref[1]
  x    : x[i][d] = coordinate d of node i
         0 <= d < dim,  0 <= i < 9

OUTPUTS:  

  x    : x[4][:] is set using edge and corner values

RETURNS: None

*/

/******************************************************************/
//  FUNCTION Prototype: xf_Transfinite3D
extern void 
xf_Transfinite3D(int dim, real *xref, real **x);
/*

PURPOSE: 

  Determines position of interior node given corners and edges, and
  face centers.

  A list of 27 nodes is used.  This is a 3x3x3 cube, with i,j,k the
  indices from fastest to slowest running. The 1,1,1 node is the one
  that is set, based on all the other nodes; its position on input is
  not used.  xref is the reference coordinate of the i,j,k node.
  
INPUTS:

  dim  : dimension of coordinate vector
  xref : 3d reference coords, xref[0], xref[1], xref[2]
  x    : x[i][d] = coordinate d of node i
         0 <= d < dim,  0 <= i < 27

OUTPUTS:  

  x    : x[13][:] is set using edge and corner values

RETURNS: None

*/



/******************************************************************/
//   FUNCTION Prototype: xf_DestroyElemSearchStructure
extern int
xf_DestroyElemSearchStructure(xf_ElemSearchStruct *ESS);
/*
PURPOSE:

  Destroys allocated search vector structure

INPUTS:
 
  ESS: element search structure

OUTPUTS: 

RETURN:

  Error Code
*/




/******************************************************************/
//   FUNCTION Prototype: xf_BuildElemSearchStructure
extern int
xf_BuildElemSearchStructure(xf_All *All, xf_ElemSearchStruct *ESS);
/*
PURPOSE:

  Creates a bin-interval-based search structure (set of search
  vectors) for quick element lookup.

INPUTS:
 
  All   : All structure

OUTPUTS: 

  ESS: pointer to element search structure

RETURN:

  Error Code
*/
  
/******************************************************************/
//   FUNCTION Prototype: xf_FindElemUsingSearchStructure
extern int
xf_FindElemUsingSearchStructure(xf_All *All, const int np, real *xpoint,
				xf_ElemSearchStruct *ESS, int *pegrp,
				int *pelem, real *xref);
/*
PURPOSE:

  Finds an element containing a point (xpoint) using element search
  vector acceleration.

INPUTS:
 
  All    : All structure
  xpoint : Coords of point that we are trying to match with an elem 
  ESS    : element search structure

OUTPUTS: 
  
  (*pegrp) : element group
  (*pelem) : element number
  xref     : reference coords in elem corresponding to xpoint

RETURN:

  xf_NOT_FOUND : if element containing xpoint was not found
  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_CurveMeshBoundary
extern int 
xf_CurveMeshBoundary(xf_Mesh *Mesh, xf_Geom *Geom, 
		     const char *BFGTitles, int Q);
/*
PURPOSE:

  Curves, to high order Q, mesh elements that are adjacent to the
  boundaries in BFGTitles.  The curved elements are placed into
  separate element groups if Q is different from existing one.

INPUTS:
 
  Mesh : mesh structure
  Geom : geometry onto which points are projected (NULL means no curving)
  BFGTitles : titles of boundaries to consider, space separated list
  Q : desired order; <=0 means keep orders as they are

OUTPUTS: 

  None : mesh is changed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_SnapToGeom
extern int 
xf_SnapToGeom(xf_All *All);
/*
PURPOSE:

  Snaps elements adjacent to boundaries to the true geometry.
  Operates on those boundaries that have a geometry component in
  All->Geom.  Wrapper for CurveMeshBoundary.  Exits quietly when there
  is nothing to do.

INPUTS:
 
  All : All structure

OUTPUTS: 

  None : mesh is possibly changed

RETURN:

  Error Code
*/




#endif // end ifndef _xf_MeshTools_h