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

#ifndef _xf_Basis_h
#define _xf_Basis_h 1

/*
  FILE:  xf_Basis.h

  This file contains the headers for functions dealing with the Basis structure.

*/

#include "xf_MeshToolsStruct.h"
#include "xf_BasisStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_Basis2Shape
extern int 
xf_Basis2Shape(enum xfe_BasisType Basis, enum xfe_ShapeType *Shape);
/*
PURPOSE:

  Converts a Basis enumerated type into a Shape enumerated type

INPUTS:

  Basis: input enumerated type

OUTPUTS: 

  Shape: output enumerated type

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Q1BasisNotLinear
extern int 
xf_Q1BasisNotLinear(enum xfe_BasisType Basis);
/*
PURPOSE:

  Returns True if a Q1 version of the basis is not necessarily linear

INPUTS:

  Basis: input enumerated type

OUTPUTS:  

  None

RETURN:

  True if basis can be nonlinear for Q1
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Shape2Dim
extern int 
xf_Shape2Dim(enum xfe_ShapeType Shape, int *Dim);
/*
PURPOSE:

  Determines dimension given a shape enumerated type

INPUTS:

  Shape: input enumerated type

OUTPUTS: 

  Dim : dimension

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Basis2Dim
extern int 
xf_Basis2Dim(enum xfe_BasisType Basis, int *Dim);
/*
PURPOSE:

  Determines dimension given a basis enumerated type

INPUTS:

  Basis: input enumerated type

OUTPUTS: 

  Dim : dimension

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_FaceShape
extern int 
xf_FaceShape(enum xfe_ShapeType Shape, int face, 
	     enum xfe_ShapeType *FShape);
/*
PURPOSE:

  Returns the shape type of a face of an element of Shape.

INPUTS:

  Shape: element shape type
  face : face number of element

OUTPUTS: 

  FShape : face shape type

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Shape2nFace
extern int
xf_Shape2nFace(enum xfe_ShapeType Shape, int *nface);
/*
PURPOSE:

  Returns the number of faces on Shape.

INPUTS:

  Shape : which element shape to consider

OUTPUTS: 

  nface: number of faces (e.g. 3 for a triangle)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Basis2nFace
extern int
xf_Basis2nFace(enum xfe_BasisType Basis, int *nface);
/*
PURPOSE:

  Returns the number of faces corresponding to interpolation using
  the basis Basis

INPUTS:

  Basis: input enumerated type

OUTPUTS: 

  nface: number of faces (e.g. 3 for a triangle)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetLocFaceMax
extern int
xf_GetLocFaceMax(int *nLocFaceMax);
/*
PURPOSE:

  Returns the max number of local faces taken over all defined shapes

INPUTS:

  None

OUTPUTS: 

  nLocFaceMax: max number of local faces over all shapes

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Shape2nPos
extern int 
xf_Shape2nPos(enum xfe_ShapeType Shape, int *npos);
/*
PURPOSE:

  Returns the total number of hanging node refinement positions
  possible for a given shape.

INPUTS:

  Shape : which element shape to consider

OUTPUTS: 

  npos : number of hanging node refinement positions.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FacePos2ElemPos
extern int 
xf_FacePos2ElemPos(enum xfe_ShapeType Shape, int face0, 
		   int facepos, int *elempos);
/*
PURPOSE:

  Returns *an* element position number (in reference to a hanging face
  refinement) corresponding to a face position number on original
  face0, for an element with shape Shape. The element position number
  may not be unique, but there should always be at least one element
  refinement position for a given face0 and facepos.

INPUTS:

  Shape   : which element shape to consider
  face0   : original face number
  facepos : face position in a hanging face refinement

OUTPUTS: 

  elempos : element position number in a hanging face refinement

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_InsideShape
extern int 
xf_InsideShape(enum xfe_ShapeType Shape, real *xref, real tol,
	       enum xfe_Bool *inside);
/*
PURPOSE:

  Determines whether the reference coordinate xref is inside the
  specified shape.

INPUTS:

  Shape: input element shape
  xref : input reference coordinate to query
  tol : tolerance on identifying inside (ref elements are ~1.0 in size)

OUTPUTS: 
 
  inside : True if xref is inside Shape within a tolerance of tol

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Order2nNode
extern int 
xf_Order2nNode(enum xfe_BasisType Basis, int p, int *nnode);
/*
PURPOSE:

  Returns the number of nodes corresponding to a Basis and order(p).

INPUTS:

  Basis: input enumerated type
  p: input numerical order

OUTPUTS: 

  nnode: number of nodes associated with Basis/p

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Q1Nodes
extern int 
xf_Q1Nodes(enum xfe_BasisType Basis, int p,
	   int *nnode, int *nvec);
/*
PURPOSE:

  Returns the Q1 nodes on an element interpolated with
  functions of type Basis using order p

INPUTS:

  Basis: type of basis
  p: order of basis

OUTPUTS: 

  nnode: number of Q1 nodes on the element
  nvec: vector of nodes (must be preallocated; not reallocated here)

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Q1NodesOnFace
extern int 
xf_Q1NodesOnFace(enum xfe_BasisType Basis, int p, int face,
		 int *nfnode, int *fvec);
/*
PURPOSE:

  Returns the Q1 nodes on a given face of an element interpolated with
  functions of type Basis using order p

INPUTS:

  Basis: type of basis
  p: order of basis
  face: face number

OUTPUTS: 

  nfnode: number of Q1 nodes on the face
  fvec: vector of nodes (must be preallocated; not reallocated here)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_Q1NodesOnFaceNeighbor
extern int 
xf_Q1NodesOnFaceNeighbor(enum xfe_BasisType Basis1, int p1, int face1,
			 int orient1, enum xfe_BasisType Basis2, int p2, 
			 int face2, int orient2, int *nfnode, int *fvec);
/*
PURPOSE:

  Similar to Q1NodesOnFace, but nodes are returned in Basis1 element's
  local ordering, but as they would look from the neighbor, Basis2
  element.

INPUTS:

  Basis1, p1: basis type and order on element 1
  Basis2, p2: basis type and order on element 2
  face1: face on element 1
  face2: face on element 2
  orient1 : orientation on element 1
  orient2 : orientation on element 2
  
OUTPUTS: 

  nfnode: number of Q1 nodes on the face
  fvec: vector of nodes (must be preallocated; not reallocated here)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_NodesOnFace
extern int 
xf_NodesOnFace(enum xfe_BasisType Basis, int p, int face,
	       int *nfnode, int *fvec);
/*
PURPOSE:

  Returns the nodes on a given face of an element interpolated with
  functions of type Basis using order p.  The node ordering should be
  consistent with the standard face shape ordering, with a positive
  orientation when viewed from the outside of the element.

INPUTS:

  Basis: type of basis
  p: order of basis
  face: face number

OUTPUTS: 

  nfnode: number of nodes on the face
  fvec: vector of nodes (must be preallocated; not reallocated here)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_LagrangeNodes
extern int 
xf_LagrangeNodes(enum xfe_BasisType Basis, int p, int *pnn, real *xn, real **pxn);
/*
PURPOSE:

  Returns reference-element Lagrange nodes of specified order p, for a
  given Basis.  The basis must be of Lagrange type.  The ref
  coordinates are stored in xn, unrolled along dimension first
  (x0,y0,z0, x1,y1,z1 ...).

INPUTS:

  Basis : basis for which Lagrange nodes are requested
  p: order of basis

OUTPUTS: 

  pnn : number of nodes (optional)
  xn: vector of Lagrange nodes, unrolled as described above.
      must be allocated before call.  Alternatively, use pxn.
  pxn : if (xn == NULL), (*pxn) is allocated and filled with
        the Lagrange nodes.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_LagrangeNodesEqual
extern int 
xf_LagrangeNodesEqual(enum xfe_ShapeType Shape, int p, real *xn);
/*
PURPOSE:

  Returns equally-spaced reference-element Lagrange nodes of specified
  order p, for a given Shape.  The ref coordinates are stored in xn,
  unrolled along dimension first (x0,y0,z0, x1,y1,z1 ...).

INPUTS:

  Shape : shape on which equally-spaced Lagrange nodes are requested
  p: order of basis

OUTPUTS: 

  xn: vector of Lagrange nodes, unrolled as described above.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Basis2Lagrange
extern int 
xf_Basis2Lagrange(enum xfe_BasisType Basis, enum xfe_BasisType *BasisLag);
/*
PURPOSE:

  Returns a Lagrange Basis corresponding to a given Basis
  (i.e. defined on same Shape).

INPUTS:

  Basis : input basis

OUTPUTS: 

  BasisLag : output LagrangeBasis

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_Basis2UniformLagrange
extern int 
xf_Basis2UniformLagrange(enum xfe_BasisType Basis, enum xfe_BasisType *BasisLag);
/*
PURPOSE:

  Returns a uniform-spaced (including end-points) Lagrange Basis
  corresponding to a given Basis (i.e. defined on same Shape).

INPUTS:

  Basis : input basis

OUTPUTS: 

  BasisLag : output LagrangeBasis

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Shape2UniformLagrange
extern int 
xf_Shape2UniformLagrange(enum xfe_ShapeType Shape, enum xfe_BasisType *BasisLag);
/*
PURPOSE:

  Returns a default uniform-spaced (including end-points) Lagrange
  Basis corresponding to a given shape.

INPUTS:

  Shape : input shape

OUTPUTS: 

  BasisLag : output LagrangeBasis

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyBasisData
extern int 
xf_DestroyBasisData( xf_BasisData *BasisData, enum xfe_Bool ForceFlag);
/*
PURPOSE:

  Destroys a BasisData structure

INPUTS:

  BasisData: structure to destroy
  ForceFlag: if false and BasisData->InTable==True, BasisData is not
             destroyed

OUTPUTS: 

  None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_BasisLagrange1D
extern int
xf_BasisLagrange1D(real x, const real *xnode, int nnode, real *phi, 
		   real *gphi, real *hphi);
/*
PURPOSE:

  Computes 1D Lagrange basis functions at point x (scalar) using nnode
  Lagrange nodes, whose positions are stored in xnode.  Can compute
  gradients and Hessians as well.

INPUTS:

  x  : scalar location at which basis function should be computed
  xnode : Lagrange node positions (nnode scalar values)
  nnode : number of Lagrange nodes

OUTPUTS: 

  phi   : nnode Lagrange basis functions, evaluated at x (optional)
          must be pre-allocated if not NULL
  gphi  : nnode Lagrange basis gradients, evaluated at x (optional)
          must be pre-allocated if not NULL
  hphi  : nnode Lagrange basis Hessians, evaluated at x (optional)
          must be pre-allocated if not NULL

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Grad_TensorLagrange
extern int 
xf_Grad_TensorLagrange(int dim, int p, const real *x, real *gphi, 
		       int off);
/*
PURPOSE:

  Returns Lagrange basis gradients of order p at point xyz for a
  tensor product element (quad or hex, as determined by dim).

INPUTS:

  dim : dimension (2 or 3)
  p : desired order of basis
  xyz: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim].  The dim individual gradients
         are offset by n in memory

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetGrads
extern int 
xf_GetGrads(enum xfe_BasisType Basis, int Order, int nq, real *xq, 
	    real *GPhi);
/*
PURPOSE:

  Fills in GPhi with ref-space gradients of basis funcs of (Basis,
  Order) evaluated at nq points xq

INPUTS:

  Basis  : basis of fcn grad desired
  Order  : order of fcn grad desired
  nq     : number of points at which to evaluate
  xq     : locations of points at which to evaluate

OUTPUTS: 

  GPhi  : basis fcn gradients, unrolled first along nn(Order), then
          along nq, and finally along dim.  Specifically
	   
	  phi0_x(x0)..phin_x(x0), phi0_x(x1)..phin_x(x1), ...
	  phi0_y(x0)..phin_y(x0), phi0_y(x1)..phin_y(x1), ...
	  ...

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_EvalBasis
extern int 
xf_EvalBasis(enum xfe_BasisType Basis, int Order, enum xfe_Bool QuadChanged,
	     int nq, real *xq, unsigned int AllocFlag, xf_BasisData **pBasisData);
/*
PURPOSE:

  Evaluates basis functions (and gradients) specific to Basis, Order
  at specified interpolation coordinates (xq).  Info is stored in
  (*pBasisData).  Does not reallocate or recalculate if (*pBasisData)
  already contains the required info.

INPUTS:

  Basis:  basis type to evaluate
  Order:  order of basis
  QuadChanged:  signifies that xq is newer than what (*pBasisData) was
                evaluated at previously
  nq : number of points at which to evaluate
  xq : coords of points at which to evaluate
  AllocFlag: 3-bit flag describing which of Phi, GPhi, gPhi to create
             or allocate.  See BasisStruct.h for details.
  
OUTPUTS: 

  (*pBasisData) : structure containing evaluated basis data (can be
                  passed in as null, and will get reallocated).

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_EvalPhysicalGrad
extern int 
xf_EvalPhysicalGrad(xf_BasisData *BasisData, xf_JacobianData *JData);
/*
PURPOSE:

  Computes physical gradients given ref-space gradients and Jacobian
  data.  BasisData->GPhi and BasisData->gPhi must be allocated, and
  JData->iJ must exist and be valid.  JData->nq must equal 1
  (signifying a constant Jacobian for all quad points) or
  BasisData->nq.

INPUTS:

  BasisData : structure containing valid GPhi and space for gPhi.
              gPhi is calculated and filled in
  JData: Jacobian data, specifically containing the inverse Jacobian
         matrix, iJ
  
OUTPUTS: 

  None, BasisData->gPhi is modified

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CreateBasisTable
extern int 
xf_CreateBasisTable( xf_BasisTable **pBasisTable);
/*
PURPOSE:

  Allocates space for a BasisData lookup table.  Specifically,
  allocates a 3d array of pointers:
  
  BasisTable->BasisData[iShape][iface][ifaceorient]

  is a (*BasisData) pointer, initialized to NULL.  The ranges for
  the array sizes are the maximum values over all defined elements.

INPUTS:

  (*pBasisTable) : structure containing BasisData lookup table
  
OUTPUTS: 

  None, (p*BasisTable) is allocated as described above

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ForceReCalcBasisTable
extern int 
xf_ForceReCalcBasisTable( xf_BasisTable *BasisTable);
/*
PURPOSE:

  Sets "ReCalc" flag in BasisTable forcing a re-calculation on next
  call to evaluate.

INPUTS:

  BasisTable : structure containing table to force re-calculation.
  
OUTPUTS: 

  None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyBasisTable
extern int 
xf_DestroyBasisTable( xf_BasisTable *BasisTable);
/*
PURPOSE:

  Destroys a BasisData lookup table.  Loops over 3d array and
  Destroys individual BasisData structures.  Releases 3d array
  and then itself.

INPUTS:

  BasisTable : structure containing table to destroy
  
OUTPUTS: 

  None

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_EvalBasisOnFaceUsingTable
extern int 
xf_EvalBasisOnFaceUsingTable(xf_Mesh *Mesh, int egrp, int elem, int face, 
			     int faceorient, enum xfe_BasisType Basis, int Order, 
			     enum xfe_Bool QuadChanged, int nq, real *xq, 
			     unsigned int AllocFlag, xf_BasisData **pBasisData,
			     xf_BasisTable *PhiTable, real **pxelem);
/*
PURPOSE:

  Evaluates basis functions (and gradients) specific to Basis, Order
  at specified coordinates (xq) which are defined on a face of the
  element.  Info is stored in (*pBasisData) which, in the generic
  element case (i.e. not cut) will point to data that is stored in
  PhiTable[Shape][face][faceorient].  PhiTable is an array of pointers
  that must be pre-allocated.  It is used to prevent recalculation of
  common basis fcn values.  When releasing data at the end of the
  calling function, be sure to destroy PhiTable and to pass ForceFlag
  == False when destroying (*pBasisData), to prevent double-releasing.

INPUTS:

  Mesh : mesh structure
  egrp, elem : element for which to calculate the basis
  face : face of element on which xq is defined
  faceorient : index describing the face's orientation within the 
               reference element.  
  Basis:  basis type to evaluate
  Order:  order of basis
  QuadChanged:  signifies that xq is newer than what was passed in
                previous call
  nq : number of points at which to evaluate
  xq : coords of points (on face) at which to evaluate
  AllocFlag: 3-bit flag describing which of Phi, GPhi, gPhi to create
             or allocate.  See BasisStruct.h for details.
  PhiTable: pre-allocated table of BasisDataPointers that stores
            previously-calculated common basis functions 

OUTPUTS: 

  (*pBasisData) : structure pointing to the evaluated basis data 
  PhiTable: modified if a new Shape, face, or faceorient was calculated.
  (*pxelem) : (optional) element-ref space coords at points 

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindMassMatrixData
extern int 
xf_FindMassMatrixData( xf_All *All, int egrp, enum xfe_BasisType Basis,
		       int Order, xf_Matrix **pM);
/*
PURPOSE:

  Locates Matrix data structure in All->Dataset containing the mass
  matrix specific to egrp, Basis, Order.  If no mass matrix is found,
  a data node is allocated in All->DataSet, memory is allocated for
  the mass matrix, and the mass matrix is computed.

INPUTS:

  All  : all structure
  egrp : element group in question
  Basis, Order : basis and order of basis functions for the mass matrix
  
  
OUTPUTS: 

  M  : matrix data structure containing the generic or specific
       mass matrices

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindGenMassMatrixData
extern int 
xf_FindGenMassMatrixData( xf_All *All, int egrp, enum xfe_BasisType Basis1,
			  int Order1, enum xfe_BasisType Basis2, int Order2,
			  xf_Matrix **pM);
/*
PURPOSE:

  Locates Matrix data structure in All->Dataset containing the general
  mass matrix specific to egrp, Basis1, Order1, Basis2, Order2.  If no
  matrix is found, a data node is allocated in All->DataSet, memory is
  allocated for the general mass matrix, and the general mass matrix
  is computed.

INPUTS:

  All  : all structure
  egrp : element group in question
  Basis1, Order1, Basis2, Order2 : basis and order of basis functions 
                                   for the general mass matrix
  
  
OUTPUTS: 

  pM : matrix data structure containing the generic or specific
       general mass matrices

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_FindInvMassMatrixData
extern int 
xf_FindInvMassMatrixData( xf_All *All, int egrp, enum xfe_BasisType Basis,
			  int Order, xf_Matrix **piM);
/*
PURPOSE:

  Locates Matrix data structure in All->Dataset containing the inverse
  mass matrix specific to egrp, Basis, Order.  If no matrix is found,
  a data node is allocated in All->DataSet, memory is allocated for
  the inverse mass matrix, and the inverse mass matrix is computed.

INPUTS:

  All  : all structure
  egrp : element group in question
  Basis, Order : basis and order of basis functions for the mass matrix
  
  
OUTPUTS: 

  iM : matrix data structure containing the generic or specific
       inverse mass matrices

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ElemMassMatrix
extern int 
xf_ElemMassMatrix( xf_All *All, int egrp, int elem, enum xfe_BasisType Basis,
		   int Order, xf_Matrix *Min, real *detJin, real **pMM, real *fac);
/*
PURPOSE:

  Returns pointer to a possibly-generic (i.e unscaled) Mass matrix MM,
  along with a scalar factor, fac, which makes the mass matrix
  element-specific.  Specifically, the true scaled mass matrix is
  MM*fac.

INPUTS:

  All  : all structure
  egrp, elem : element in question
  Basis, Order : basis and order of basis functions for the mass matrix
  Min  : input matrix data structure containing the generic or specific
         mass matrices.  Optional.  If not given, a separate call to 
	 FindMassMatrixData will be made.
  detJin : input geometry-Jacobian determinant of elem, for constant-J elems.  
           If not given, the linear elem volume will be calculated to prepare
	   the appropriate scaling, fac.
  
OUTPUTS: 

  (*pMM) : pointer to mass matrix.  Do not deallocate this.
  (*fac) : factor that makes the mass matrix elem-specific, in case of a
           generic (*pMM)

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ElemGenMassMatrix
extern int 
xf_ElemGenMassMatrix( xf_All *All, int egrp, int elem, enum xfe_BasisType Basis1,
		      int Order1, enum xfe_BasisType Basis2, int Order2,
		      xf_Matrix *Min, real *detJin, real **pMM, real *fac);
/*
PURPOSE:

  Returns pointer to a possibly-generic (i.e unscaled) general Mass
  matrix MM, along with a scalar factor, fac, which makes the general
  mass matrix element-specific.  Specifically, the true scaled general
  mass matrix is MM*fac.

INPUTS:

  All  : all structure
  egrp, elem : element in question
  Basis1, Order1, Basis2, Order2 : basis and order of basis functions 
                                   for the general mass matrix
  Min  : input matrix data structure containing the generic or specific
         general mass matrices.  Optional.  If not given, a separate call to 
	 FindGenMassMatrixData will be made.
  detJin : input geometry-Jacobian determinant of elem, for constant-J elems.  
           If not given, the linear elem volume will be calculated to prepare
	   the appropriate scaling, fac.
  
OUTPUTS: 

  (*pMM) : pointer to general mass matrix.  Do not deallocate or alter.
  (*fac) : factor that makes the inverse mass matrix elem-specific, 
           in case of a generic (*pMM)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ElemInvMassMatrix
extern int 
xf_ElemInvMassMatrix( xf_All *All, int egrp, int elem, enum xfe_BasisType Basis,
		      int Order, xf_Matrix *iMin, real *detJin, real **piMM, 
		      real *fac);
/*
PURPOSE:

  Returns pointer to a possibly-generic (i.e unscaled) inverse Mass
  matrix iMM, along with a scalar factor, fac, which makes the inverse
  mass matrix element-specific.  Specifically, the true scaled inverse
  mass matrix is iMM*fac.

INPUTS:

  All  : all structure
  egrp, elem : element in question
  Basis, Order : basis and order of basis functions for the mass matrix
  iMin  : input matrix data structure containing the generic or specific
         inverse mass matrices.  Optional.  If not given, a separate call to 
	 FindInvMassMatrixData will be made.
  detJin : input geometry-Jacobian determinant of elem, for constant-J elems.  
           If not given, the linear elem volume will be calculated to prepare
	   the appropriate scaling, fac.
  
OUTPUTS: 

  (*piMM) : pointer to inverse mass matrix.  Do not deallocate or alter.
  (*fac) : factor that makes the inverse mass matrix elem-specific, 
           in case of a generic (*piMM)

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ComputeTransferMatrix
extern int 
xf_ComputeTransferMatrix(enum xfe_BasisType Basis1, int Order1,
			 enum xfe_BasisType Basis2, int Order2, 
			 enum xfe_ShapeType EShape, int pos, real *T);
/*
PURPOSE:

  Computes a transfer matrix suitable for converting data interpolated
  with Basis1,Order1 (elem1) to Basis2,Order2 (elem2).  elem2 may be a
  sub-element of elem1, in which case EShape and pos are used to
  discern how elem2 is positioned inside of elem1.

INPUTS:

  Basis1, Order1 : basis and order on elem1 (possibly background element)
  Basis2, Order2 : basis and order on elem2 (possibly a sub element)
  EShape: shape of element
  pos : position of fine element within the background element  
  
OUTPUTS: 

  T  : Transfer matrix

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindTransferMatrix
extern int 
xf_FindTransferMatrix( xf_DataSet *DataSet, enum xfe_BasisType Basis1,
		       int Order1, enum xfe_BasisType Basis2, 
		       int Order2, xf_Matrix **pT);
/*
PURPOSE:

  Returns a Transfer matrix suitable for converting data interpolated
  with Basis1,Order1 to data interpolated with Basis2,Order2.  This
  matrix will be a restriction, prolongation, or inter-basis transfer
  matrix.  If DataSet is passed in, DataSet is searched for the
  matrix; if not found, the matrix is saved in DataSet.  (*pT) should
  not be released following use in this case.  On the other hand, if
  DataSet==NULL is passed in, (*pT) is re-allocated, and should be
  released by the calling function when finished.

INPUTS:

  DataSet : DataSet structure (optional, see above discussion)
  Basis1, Order1 : current basis and order
  Basis2, Order2 : desired basis and order
  
  
OUTPUTS: 

  (*pT)  : Transfer matrix data structure

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectOnElemQR
extern int 
xf_ProjectOnElemQR(xf_Mesh *Mesh, int egrp, int elem, int sr, 
		   xf_QuadData *QuadData, enum xfe_BasisType Basis, 
		   int Order, enum xfe_Bool Regenerate, real **pQ, 
		   real **pR, int **pP, real *uin, real *U);
/*
PURPOSE:

  Projects point values u at quadrature points stored in QuadData to a
  polynomial of Basis and Order.  Projection coefficients of the
  polynomials are stored in U, which must be adequately allocated.
  Projection is performed using QR factorization.

INPUTS:

  Mesh : Mesh structure
  egrp, elem : element on which to do the projection
  sr : state rank in u
  QuadData : quadrature data structure containing points, weights
             NOTE: number of quad points, nq, must be equal to or more
	           than the number of basis functions for Basis,Order
  Basis, Order : desired projection basis and order
  Regenerate : True to regenerate Q,P,R matrices
  pQ, pR, pP : pointers to Q,P,R matrices so that the calling function can
               save time by passing in previous matrices if things have
	       not changed.
  u : unrolled state at all the quad points
  
OUTPUTS: 

  U : coefficients of polynomial functions in projection
  (*pQ), (*pR), (*pP) : altered Q,P,R matrices

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectOnElemQR_OnlyEnergy
extern int 
xf_ProjectOnElemQR_OnlyEnergy(xf_Mesh *Mesh, int egrp, int elem, int sr, 
		              xf_QuadData *QuadData, enum xfe_BasisType Basis, 
		              int Order, enum xfe_Bool Regenerate, real **pQ, 
		              real **pR, int **pP, real *uin, real *Uout);

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectOnElemLeastSquares
extern int 
xf_ProjectOnElemLeastSquares(xf_Mesh *Mesh, int egrp, int elem, int sr, 
			     xf_QuadData *QuadData, enum xfe_BasisType Basis, 
			     int Order, enum xfe_Bool Regenerate, real **pQ, 
			     real **pR, int **pP, real *uin, real *U);
/*
PURPOSE:

  Projects point values u at quadrature points stored in QuadData to a
  polynomial of Basis and Order.  Projection coefficients of the
  polynomials are stored in U, which must be adequately allocated.
  Projection is performed using least-squares.

INPUTS:

  Mesh : Mesh structure
  egrp, elem : element on which to do the projection
  sr : state rank in u
  QuadData : quadrature data structure containing points, weights
             NOTE: number of quad points, nq, must be equal to or more
	           than the number of basis functions for Basis,Order
  Basis, Order : desired projection basis and order
  Regenerate : True to regenerate Q,P,R matrices
  pQ, pR, pP : pointers to Q,P,R matrices so that the calling function can
               save time by passing in previous matrices if things have
	       not changed.
  u : unrolled state at all the quad points
  
OUTPUTS: 

  U : coefficients of polynomial functions in projection
  (*pQ), (*pR), (*pP) : altered Q,P,R matrices

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ProjectOnElemLeastSquares
extern int 
xf_ProjectOnElemLeastSquares_OnlyEnergy(xf_Mesh *Mesh, int egrp, int elem, int sr, 
			     xf_QuadData *QuadData, enum xfe_BasisType Basis, 
			     int Order, enum xfe_Bool Regenerate, real **pQ, 
			     real **pR, int **pP, real *uin, real *U);

/******************************************************************/
//   FUNCTION Prototype: xf_ShapeIsoFace
extern int 
xf_ShapeIsoFace(xf_All *All, int egrp, int elem, int face, int nq, real *xq, real *phi);
/*
PURPOSE:

  Evaluates face-isolated shape function on elem in question, at the
  desired points.  This function is nonzero only on the face in
  question and in the interior of the domain.  It is also positive
  everywhere.  Useful for interpolating a quantity from a face to the
  element interior.


INPUTS:

  All : all structure
  egrp, elem : element in question
  face : face on elem for which to calculate iso shape function
  nq : number of points
  xq : elem ref space locations of the points
  
OUTPUTS: 

  phi : shape function at all nq points

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_ShapeCentroid
extern int 
xf_ShapeCentroid(enum xfe_ShapeType Shape, real *xref);
/*
PURPOSE:

  Evaluates centroid coordinates of a given shape

INPUTS:

  Shape : element shape
  
OUTPUTS: 

  xref : reference coordinates of centroid

RETURN:

  Error Code
*/


#endif // end ifndef _xf_Basis_h
