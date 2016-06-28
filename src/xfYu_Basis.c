/*------------------------------------------------------------------*/
/* this routine is inherited from xf_Basis.c but targets to add     */
/* the implementation to support 2D and 3D Legendre basis by Yu     */
/* Authoer: Yu Lv;  Email: ylv@stanford.edu                         */
/*------------------------------------------------------------------*/

/*
  FILE:  xfYu_Basis.c

  This file contains functions for working with the Bases.

*/


#include "xf_AllStruct.h"
#include "xf_Math.h"
#include "xf_Memory.h"
#include "xf_BasisFcn.h"
#include "xf_MeshTools.h"
#include "xf_QuadRule.h"
#include "xf_IO.h"
#include "xf_Quad.h" 
#include "xf_Param.h" 
#include "xf_Data.h" 
#include "xf_DataMath.h" 

/******************************************************************/
//   FUNCTION Definition: xf_Basis2Shape
int 
xf_Basis2Shape(enum xfe_BasisType Basis, enum xfe_ShapeType *Shape){

  switch (Basis){
  case xfe_SegLagrange:
  case xfe_SegLagrangeGauss:
    (*Shape) = xfe_Segment;
    break;
  case xfe_TriLagrange:
  case xfe_TriHierarch:
    (*Shape) = xfe_Triangle;
    break;
  case xfe_TetLagrange:
  case xfe_TetHierarch:
    (*Shape) = xfe_Tetrahedron;
    break;
  case xfe_QuadLagrange: 
  case xfe_QuadLagrangeGauss:
  case xfe_QuadLegendre:
    (*Shape) = xfe_Quadrilateral;
    break;
  case xfe_HexLagrange:
  case xfe_HexLagrangeGauss:
  case xfe_HexLegendre:
    (*Shape) = xfe_Hexahedron;
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Q1BasisNotLinear
enum xfe_Bool
xf_Q1BasisNotLinear(enum xfe_BasisType Basis)
{
  enum xfe_Bool rval;

  switch (Basis){
  case xfe_QuadLagrange: 
  case xfe_QuadLagrangeGauss:
  case xfe_QuadLegendre:
  case xfe_HexLagrange:
  case xfe_HexLagrangeGauss:
  case xfe_HexLegendre:
    rval = xfe_True;
    break;
  default:
    rval = xfe_False;
    break;
  }

  return rval;
}


/******************************************************************/
//   FUNCTION Definition: xf_Shape2Dim
int 
xf_Shape2Dim(enum xfe_ShapeType Shape, int *Dim)
{
  switch (Shape){
  case xfe_Segment:
    (*Dim) = 1;
    break;
  case xfe_Quadrilateral:
  case xfe_Triangle:
    (*Dim) = 2;
    break;
  case xfe_Tetrahedron:
  case xfe_Hexahedron: 
    (*Dim) = 3;
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Basis2Dim
int 
xf_Basis2Dim(enum xfe_BasisType Basis, int *Dim){

  switch (Basis){
  case xfe_SegLagrange:
  case xfe_SegLagrangeGauss:
    (*Dim) = 1;
    break;
  case xfe_TriLagrange:
  case xfe_TriHierarch:
  case xfe_QuadLagrange: 
  case xfe_QuadLagrangeGauss:
  case xfe_QuadLegendre:
    (*Dim) = 2;
    break;
  case xfe_TetLagrange:
  case xfe_TetHierarch:
  case xfe_HexLagrange:
  case xfe_HexLagrangeGauss:
  case xfe_HexLegendre:
    (*Dim) = 3;
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FaceShape
int 
xf_FaceShape(enum xfe_ShapeType Shape, int face, 
	     enum xfe_ShapeType *FShape)
{
  
  switch (Shape){
  case xfe_Segment:
    (*FShape) = xfe_Point;
    break;
  case xfe_Triangle:
    (*FShape) = xfe_Segment;
    break;
  case xfe_Tetrahedron:
    (*FShape) = xfe_Triangle;
    break;
  case xfe_Quadrilateral:
    (*FShape) = xfe_Segment;
    break;
  case xfe_Hexahedron:
    (*FShape) = xfe_Quadrilateral;
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Shape2nFace
int 
xf_Shape2nFace(enum xfe_ShapeType Shape, int *nface)
{
  
  (*nface) = 0;

  switch (Shape){
  case xfe_Point:   /* Point      */
    (*nface) = 0;
    break;
  case xfe_Segment:   /* Segment      */
    (*nface) = 2;
    break;
  case xfe_Triangle:   /* Triangle      */
    (*nface) = 3;
    break;
  case xfe_Tetrahedron:   /* Tetrahedron   */
    (*nface) = 4;
    break;
  case xfe_Quadrilateral:  /* Quadrilateral */
    (*nface) = 4;
    break;
  case xfe_Hexahedron:   /* Hexahedron    */
    (*nface) = 6;
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Basis2nFace
int 
xf_Basis2nFace(enum xfe_BasisType Basis, int *nface){

  int ierr;
  enum xfe_ShapeType Shape;

  ierr = xf_Error(xf_Basis2Shape(Basis, &Shape));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Shape2nFace(Shape, nface));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GetLocFaceMax
int
xf_GetLocFaceMax(int *nLocFaceMax){

  int ierr;
  enum xfe_ShapeType Shape;
  int nface;

  (*nLocFaceMax) = 0;

  for (Shape=0; Shape < xfe_ShapeLast; Shape++){
    ierr = xf_Error(xf_Shape2nFace(Shape, &nface));
    if (ierr != xf_OK) return ierr;
    if (nface > (*nLocFaceMax)) (*nLocFaceMax) = nface;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Shape2nFaceOrient
static int 
xf_Shape2nFaceOrient(enum xfe_ShapeType Shape, int *nfaceorient){
  
  /* returns number of possible face orientations w.r.t the reference
     element.  For example, a triangle face on a tetrahedron can have
     its origin ("0") node in 3 possible locations and the face node
     ordering could be CW or CCW as viewed from outside the tet -- for
     a total of 2*3 = 6 possible face orientations. */

  (*nfaceorient) = 0;

  switch (Shape){
  case xfe_Point:
    (*nfaceorient) = 1;
    break;
  case xfe_Segment:
    (*nfaceorient) = 1;
    break;
  case xfe_Triangle:
    (*nfaceorient) = 2;
    break;
  case xfe_Tetrahedron:   
    (*nfaceorient) = 6;
    break;
  case xfe_Quadrilateral:
    (*nfaceorient) = 2;
    break;
  case xfe_Hexahedron:
    (*nfaceorient) = 8;
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetLocFaceOrientMax
static int
xf_GetLocFaceOrientMax(int *nLocFaceOrientMax){

  int ierr;
  enum xfe_ShapeType Shape;
  int nfaceorient;

  (*nLocFaceOrientMax) = 0;

  for (Shape=0; Shape < xfe_ShapeLast; Shape++){
    ierr = xf_Error(xf_Shape2nFaceOrient(Shape, &nfaceorient));
    if (ierr != xf_OK) return ierr;
    if (nfaceorient > (*nLocFaceOrientMax)) 
      (*nLocFaceOrientMax) = nfaceorient;
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Shape2nPos
int 
xf_Shape2nPos(enum xfe_ShapeType Shape, int *npos)
{
  switch (Shape){
  case xfe_Segment:
    (*npos) = xfe_SegPosLast;
    break;
  case xfe_Quadrilateral:
    (*npos) = xfe_QuadPosLast;
    break;
  case xfe_Triangle:
    (*npos) = xfe_TriPosLast;
    break;
  case xfe_Hexahedron: 
    (*npos) = xfe_HexPosLast;
    break;
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
//   FUNCTION Definition: xf_FacePosElemPos
int 
xf_FacePos2ElemPos(enum xfe_ShapeType Shape, int face0, 
		   int facepos, int *elempos)
{
  int F2EHex[6][xfe_QuadPosLast] = {
    {0,xfe_HexPos000,xfe_HexPos010,xfe_HexPos100,xfe_HexPos110,xfe_HexPos202,xfe_HexPos212,xfe_HexPos022,xfe_HexPos122},
    {0,xfe_HexPos000,xfe_HexPos100,xfe_HexPos001,xfe_HexPos101,xfe_HexPos022,xfe_HexPos122,xfe_HexPos220,xfe_HexPos221},
    {0,xfe_HexPos100,xfe_HexPos110,xfe_HexPos101,xfe_HexPos111,xfe_HexPos202,xfe_HexPos212,xfe_HexPos220,xfe_HexPos221},
    {0,xfe_HexPos110,xfe_HexPos010,xfe_HexPos111,xfe_HexPos011,xfe_HexPos122,xfe_HexPos022,xfe_HexPos220,xfe_HexPos221},
    {0,xfe_HexPos010,xfe_HexPos000,xfe_HexPos011,xfe_HexPos001,xfe_HexPos212,xfe_HexPos202,xfe_HexPos220,xfe_HexPos221},
    {0,xfe_HexPos001,xfe_HexPos101,xfe_HexPos011,xfe_HexPos111,xfe_HexPos022,xfe_HexPos122,xfe_HexPos202,xfe_HexPos212},
  };
  

  (*elempos) = 0;
  if (facepos == 0) return xf_OK;

  switch (Shape){
  case xfe_Segment:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  case xfe_Quadrilateral:
    switch (face0){
    case 0:
      if (facepos==xfe_SegPosLeft) (*elempos) = xfe_QuadPosSW;
      else (*elempos) = xfe_QuadPosSE;
      break;
    case 1:
      if (facepos==xfe_SegPosLeft) (*elempos) = xfe_QuadPosSE;
      else (*elempos) = xfe_QuadPosNE;
      break;
    case 2:
      if (facepos==xfe_SegPosLeft) (*elempos) = xfe_QuadPosNE;
      else (*elempos) = xfe_QuadPosNW;
      break;
    case 3:
      if (facepos==xfe_SegPosLeft) (*elempos) = xfe_QuadPosNW;
      else (*elempos) = xfe_QuadPosSW;
      break;
    default:
      return xf_Error(xf_INPUT_ERROR);
      break;
    }
    break;
  case xfe_Triangle:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  case xfe_Hexahedron: 
    if ((face0<0) || (face0>5) || (facepos<0) || (facepos>=xfe_QuadPosLast))
      return xf_Error(xf_INPUT_ERROR);
    (*elempos) = F2EHex[face0][facepos];
    break;
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
//   FUNCTION Definition: xf_InsideShape
int 
xf_InsideShape(enum xfe_ShapeType Shape, real *xref, real tol,
	       enum xfe_Bool *inside)
{
  switch (Shape){
  case xfe_Segment:
    (*inside) = ((xref[0] > -tol) && (xref[0] < 1.0+tol));
    break;
  case xfe_Quadrilateral:
    (*inside) = ( (xref[0] > -tol) && (xref[0] < 1.0+tol) &&
		  (xref[1] > -tol) && (xref[1] < 1.0+tol) );
    break;
  case xfe_Triangle:
    (*inside) = ( (xref[0] > -tol) && (xref[1] > -tol) &&
		  ((xref[0]+xref[1]) < 1.0+tol) );
    break;
  case xfe_Tetrahedron:
    (*inside) = ( (xref[0] > -tol) && (xref[1] > -tol) && (xref[2] > -tol) &&
		  ((xref[0]+xref[1]+xref[2]) < 1.0+tol) );
    break;
  case xfe_Hexahedron: 
    (*inside) = ( (xref[0] > -tol) && (xref[0] < 1.0+tol) &&
		  (xref[1] > -tol) && (xref[1] < 1.0+tol) &&
		  (xref[2] > -tol) && (xref[2] < 1.0+tol) );
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }
  return xf_OK;
}






/******************************************************************/
//   FUNCTION Definition: xf_Order2nNode
int 
xf_Order2nNode(enum xfe_BasisType Basis, int p, int *nnode){

  switch (Basis){
  case xfe_SegLagrange:
  case xfe_SegLagrangeGauss:
    (*nnode) = (p+1);
    break;
  case xfe_TriLagrange:
  case xfe_TriHierarch:
    (*nnode) = (p+1)*(p+2)/2;
    break;
  case xfe_TetLagrange:
  case xfe_TetHierarch:
    (*nnode) = (p+1)*(p+2)*(p+3)/6;
    break;
  case xfe_QuadLagrange: 
  case xfe_QuadLagrangeGauss:
  case xfe_QuadLegendre:
    (*nnode) = (p+1)*(p+1);
    break;
  case xfe_HexLagrange:
  case xfe_HexLagrangeGauss:
  case xfe_HexLegendre:
    (*nnode) = (p+1)*(p+1)*(p+1);
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Q1Nodes
int 
xf_Q1Nodes(enum xfe_BasisType Basis, int p,
	   int *nnode, int *nvec){
  switch (Basis){
  case xfe_SegLagrange:
    (*nnode) = 2;
    nvec[0] = 0; nvec[1] = p;
    break;
  case xfe_TriLagrange:
    (*nnode) = 3;
    nvec[0] = 0; nvec[1] = p; nvec[2] = (p+1)*(p+2)/2-1;
    break;
  case xfe_TetLagrange:
    (*nnode) = 4;
    nvec[0] = 0;
    nvec[1] = p;
    nvec[2] = (p+1)*(p+2)/2-1;
    nvec[3] = (p+1)*(p+2)*(p+3)/6-1;
    break;
  case xfe_QuadLagrange: 
    (*nnode) = 4;
    nvec[0] = 0; nvec[1] = p;
    nvec[2] = (p+1)*p; nvec[3] = (p+2)*p;
    break;
  case xfe_HexLagrange:
    (*nnode) = 8;
    nvec[0] = 0; nvec[1] = p; nvec[2] = (p+1)*(p+1)-(p+1); nvec[3] = (p+1)*(p+1)-1;
    nvec[7] = (p+1)*(p+1)*(p+1)-1;
    nvec[6] = nvec[7]-p; nvec[4] = nvec[7]-(p+1)*(p+1)+1; nvec[5] = nvec[4]+p;
    break;
  case xfe_TetHierarch:
  case xfe_TriHierarch:
  case xfe_QuadLagrangeGauss:
  case xfe_HexLagrangeGauss:
  case xfe_QuadLegendre:
  case xfe_HexLegendre:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Q1NodesOnFace
int 
xf_Q1NodesOnFace(enum xfe_BasisType Basis, int p, int face,
		 int *nfnode, int *fvec){
  int n[8];

  switch (Basis){
    case xfe_SegLagrange:
    (*nfnode) = 1;
    switch (face){
    case 0: fvec[0] = 0; break;
    case 1: fvec[0] = p; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    break;
  case xfe_TriLagrange:
    (*nfnode) = 2;
    switch (face){
    case 0: fvec[0] = p; fvec[1] = (p+1)*(p+2)/2-1; break;
    case 1: fvec[0] = (p+1)*(p+2)/2-1; fvec[1] = 0; break;
    case 2: fvec[0] = 0; fvec[1] = p; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    break;
  case xfe_TetLagrange:
    n[0] = 0;
    n[1] = p;
    n[2] = (p+1)*(p+2)/2-1;
    n[3] = (p+1)*(p+2)*(p+3)/6-1;
    (*nfnode) = 3;
    switch (face){
    case 0: fvec[0] = n[1]; fvec[1] = n[2]; fvec[2] = n[3]; break;
    case 1: fvec[0] = n[0]; fvec[1] = n[3]; fvec[2] = n[2]; break;
    case 2: fvec[0] = n[0]; fvec[1] = n[1]; fvec[2] = n[3]; break;
    case 3: fvec[0] = n[0]; fvec[1] = n[2]; fvec[2] = n[1]; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    break;
  case xfe_QuadLagrange: 
    (*nfnode) = 2;
    switch (face){
    case 0: fvec[0] = 0; fvec[1] = p; break;
    case 1: fvec[0] = p; fvec[1] = (p+2)*p; break;
    case 2: fvec[0] = (p+2)*p; fvec[1] = (p+1)*p; break;
    case 3: fvec[0] = (p+1)*p; fvec[1] = 0; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    break;
  case xfe_HexLagrange:
    n[0] = 0; n[1] = p; n[2] = (p+1)*(p+1)-(p+1); n[3] = (p+1)*(p+1)-1;
    n[7] = (p+1)*(p+1)*(p+1)-1;
    n[6] = n[7]-p; n[4] = n[7]-(p+1)*(p+1)+1; n[5] = n[4]+p;
    (*nfnode) = 4;
    switch (face){
    case 0: 
      fvec[0] = n[0]; fvec[1] = n[2]; fvec[2] = n[3]; fvec[3] = n[1]; 
      break;
    case 1: 
      fvec[0] = n[0]; fvec[1] = n[1]; fvec[2] = n[5]; fvec[3] = n[4]; 
      break;
    case 2: 
      fvec[0] = n[1]; fvec[1] = n[3]; fvec[2] = n[7]; fvec[3] = n[5]; 
      break;
    case 3: 
      fvec[0] = n[3]; fvec[1] = n[2]; fvec[2] = n[6]; fvec[3] = n[7]; 
      break;
    case 4: 
      fvec[0] = n[2]; fvec[1] = n[0]; fvec[2] = n[4]; fvec[3] = n[6]; 
      break;
    case 5: 
      fvec[0] = n[4]; fvec[1] = n[5]; fvec[2] = n[7]; fvec[3] = n[6]; 
      break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    break;

  case xfe_TetHierarch:
  case xfe_TriHierarch:
  case xfe_QuadLagrangeGauss:
  case xfe_HexLagrangeGauss:
  case xfe_QuadLegendre:
  case xfe_HexLegendre:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_RotateNodes
static int 
xf_RotateNodes(enum xfe_ShapeType FShape, int orient, int nnode, 
	       int *nvec1, int *nvec2)
{
  // Rotates face node order to as-seen by face
  int cycle, flip, i;
  switch (FShape){
  case xfe_Segment:
    if (orient == 0){
      nvec2[0] = nvec1[0]; nvec2[1] = nvec1[1];
    }
    else{
      nvec2[0] = nvec1[1]; nvec2[1] = nvec1[0];
    }
    break;
  case xfe_Quadrilateral:
    cycle = orient%4;
    flip  = orient/4;
    for (i=0; i<4; i++) nvec2[i] = nvec1[(i+cycle)%4];
    if (flip) swap(nvec2[1], nvec2[3], i);
    break;
  case xfe_Triangle:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  case xfe_Tetrahedron:
  case xfe_Hexahedron:
    return xf_Error(xf_INPUT_ERROR);
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_InvRotateNodes
static int 
xf_InvRotateNodes(enum xfe_ShapeType FShape, int orient, int nnode, 
		  int *nvec1, int *nvec2)
{
  // Rotates face node order to as-seen by element
  int cycle, flip, i;
  int ivec[4];

  switch (FShape){
  case xfe_Segment:
    if (orient == 0){
      nvec2[0] = nvec1[0]; nvec2[1] = nvec1[1];
    }
    else{
      nvec2[0] = nvec1[1]; nvec2[1] = nvec1[0];
    }
    break;
  case xfe_Quadrilateral:
    cycle = orient%4;
    flip  = orient/4;
    for (i=0; i<4; i++) ivec[i] = nvec1[i];
    if (flip) swap(ivec[1], ivec[3], i);
    for (i=0; i<4; i++) nvec2[i] = ivec[(i+4-cycle)%4];
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
//   FUNCTION Definition: xf_Q1NodesOnFaceNeighbor
int 
xf_Q1NodesOnFaceNeighbor(enum xfe_BasisType Basis1, int p1, int face1,
			 int orient1, enum xfe_BasisType Basis2, int p2, 
			 int face2, int orient2, int *nfnode, int *fvec)
{
  int ierr;
  int fvec1[8], fvec2[8];
  enum xfe_ShapeType Shape1, FShape;
  
  // first get nodes from point of view of Basis1
  ierr = xf_Error(xf_Q1NodesOnFace(Basis1, p1, face1, nfnode, fvec1));
  if (ierr != xf_OK) return ierr;

  // get Shape1
  ierr = xf_Error(xf_Basis2Shape(Basis1, &Shape1));
  if (ierr != xf_OK) return ierr;

  // Face shape
  ierr = xf_Error(xf_FaceShape(Shape1, face1, &FShape));
  if (ierr != xf_OK) return ierr;
  
  // Rotate nodes to face ref
  ierr = xf_Error(xf_RotateNodes(FShape, orient1, (*nfnode), fvec1, fvec2));
  if (ierr != xf_OK) return ierr;

  // inverse rotate nodes from face ref to Basis2 view
  ierr = xf_Error(xf_InvRotateNodes(FShape, orient2, (*nfnode), fvec2, fvec));
  if (ierr != xf_OK) return ierr;
  

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_NodesOnFace
int 
xf_NodesOnFace(enum xfe_BasisType Basis, int p, int face,
	       int *nfnode, int *fvec){
  int i0, d, d0, d1, i, j, k, jx, jy, jz;
  int i1, i2, ix, iy, iz, D0, D1, count, q, r;

  switch (Basis){
  case xfe_TriLagrange:
    (*nfnode) = (p+1);
    switch(face){
    case 0: i0 = p; d0 = p; d1 = -1; break;
    case 1: i0 = (p+1)*(p+2)/2-1; d0 = -2; d1 = -1; break;
    case 2: i0 = 0;  d0 = 1; d1 = 0; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    fvec[0] = i0;
    for (i=1, d=d0; i<p+1; i++, d+=d1) fvec[i] = fvec[i-1]+d;
    break;
  case xfe_TetLagrange:
    (*nfnode) = (p+1)*(p+2)/2;
    switch(face){
    case 0:  ix=p; iy=0; iz=0; D0 = -1; D1 = -2; break;
    case 1:  ix=0; iy=0; iz=0; D0 =  2; D1 =  1; break;
    case 2:  ix=0; iy=0; iz=0; D0 =  0; D1 =  2; break;
    case 3:  ix=0; iy=0; iz=0; D0 =  1; D1 =  0; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    count = 0;
    for (i2=0; i2<p+1; i2++){
      jx = ix; jy = iy; jz = iz;
      for (i1=0; i1<(p+1-i2); i1++){
	q = p-jz;
	r = q-jy;
	i = (p+1)*(p+2)*(p+3)/6-(q+1)*(q+2)*(q+3)/6 + (q+1)*(q+2)/2-(r+1)*(r+2)/2 + jx;
	fvec[count] = i;
	count++;
	d = D0;
	if (d < 0){ jx--; d = -d;}
	if (d == 0) jx++;
	else if (d == 1) jy++;
	else if (d == 2) jz++;
      }
      d = D1;
      if (d < 0){ ix--; d = -d;}
      if (d == 0) ix++;
      else if (d == 1) iy++;
      else if (d == 2) iz++;
    }

    break;
  case xfe_QuadLagrange: 
    (*nfnode) = (p+1);
    switch (face){
    case 0: i0 = 0;       d =    1; break;
    case 1: i0 = p;       d =  p+1; break;
    case 2: i0 = p*(p+2); d =   -1; break;
    case 3: i0 = p*(p+1); d = -p-1; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    for (i=0; i<p+1; i++)fvec[i] = i0+i*d;
    break;
  case xfe_HexLagrange:
    (*nfnode) = (p+1)*(p+1);
    switch (face){
    case 0: i0 = 0;             d0 =   p+1 ; d1 =          1 ; break;
    case 1: i0 = 0;             d0 =     1 ; d1 = (p+1)*(p+1); break;
    case 2: i0 = p;             d0 =   p+1 ; d1 = (p+1)*(p+1); break;
    case 3: i0 = (p+1)*(p+1)-1; d0 =    -1 ; d1 = (p+1)*(p+1); break;
    case 4: i0 = p*(p+1);       d0 = -(p+1); d1 = (p+1)*(p+1); break;
    case 5: i0 = p*(p+1)*(p+1); d0 =     1 ; d1 =        p+1 ; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    for (j=0,k=0; j<p+1; j++)
      for (i=0; i<p+1; i++,k++)
	fvec[k] = i0 + i*d0 + j*d1;
    break;
  case xfe_QuadLagrangeGauss:
  case xfe_QuadLegendre:
  case xfe_HexLagrangeGauss:
  case xfe_HexLegendre:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  case xfe_TetHierarch:
  case xfe_TriHierarch:
    return xf_Error(xf_NOT_SUPPORTED);
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_LagrangeNodes
int 
xf_LagrangeNodes(enum xfe_BasisType Basis, int p, int *pnn, real *xn0, 
		 real **pxn){
  int ierr, ix, iy, iz, k, nn, dim;
  real *x1 = NULL, dx, *xn;
  
  // number of nodes
  ierr = xf_Error(xf_Order2nNode(Basis, p, &nn));
  if (ierr != xf_OK) return ierr;

  if (pnn != NULL) (*pnn) = nn;

  // allocate xn if necessary
  if (xn0 == NULL){
    ierr = xf_Error(xf_Basis2Dim(Basis, &dim));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReAlloc((void **) pxn, dim*nn, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    xn = (*pxn);
  }
  else
    xn = xn0;


  // allocate x1 for tensor product nodes
  if ((p>0) && ((Basis == xfe_QuadLagrange) || (Basis == xfe_HexLagrange))){
    ierr = xf_Error(xf_Alloc((void **) &x1, p+1, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (k=0, dx = 1./((real) p); k<p+1; k++) x1[k] = k*dx;
  }
  if ((p>0) && ((Basis == xfe_QuadLagrangeGauss) || (Basis == xfe_HexLagrangeGauss))){
    ierr = xf_Error(xf_QuadLine(2*p+1, &k, &x1, NULL));
    if (ierr != xf_OK) return ierr;
    if (k != p+1) return xf_Error(xf_BASIS_FCN_ERROR);
  }

  // fill xn appropriate to Basis
  switch (Basis){
  case xfe_SegLagrange:
  case xfe_SegLagrangeGauss:
    if (p == 0){
      xn[0] = 0.0;
    }
    else{
      for (ix=0; ix<=p; ix++)
	xn[ix] = ( (real) ix) / ( (real) p);
    }   
    break;

  case xfe_TriLagrange:
    if (p == 0){
      xn[0] = xn[1] = 0.0;
    }
    else{
      k = 0;
      for (iy=0; iy<=p; iy++){ 
	for (ix=0; ix<=(p-iy); ix++){ 
	  xn[2*k+0] = ( (real) ix) / ( (real) p);
	  xn[2*k+1] = ( (real) iy) / ( (real) p);
	  k++;
	} // ix
      } // iy
    }   
    break;
  case xfe_QuadLagrange:
  case xfe_QuadLagrangeGauss:
    if (p == 0){
      xn[0] = xn[1] = 0.0;
    }
    else{
      for (iy=0,k=0; iy<p+1; iy++)
	for (ix=0; ix<p+1; ix++,k++){
	  xn[2*k+0] = x1[ix];
	  xn[2*k+1] = x1[iy];
	}
    }
    break;
  case xfe_TetLagrange:
    if (p == 0){
      xn[0] = xn[1] = xn[2] = 0.0;
    }
    else{
      k = 0;
      for (iz=0; iz<=p; iz++){ 
	for (iy=0; iy<=(p-iz); iy++){ 
	  for (ix=0; ix<=(p-iz-iy); ix++){ 
	    xn[3*k+0] = ( (real) ix) / ( (real) p);
	    xn[3*k+1] = ( (real) iy) / ( (real) p);
	    xn[3*k+2] = ( (real) iz) / ( (real) p);
	    k++;
	  } // ix
	} // iy
      } // iz
    }
    break;
  case xfe_HexLagrange:
  case xfe_HexLagrangeGauss:
    if (p == 0){
      xn[0] = xn[1] = 0.0;
    }
    else{
      for (iz=0,k=0; iz<p+1; iz++)
	for (iy=0; iy<p+1; iy++)
	  for (ix=0; ix<p+1; ix++,k++){
	    xn[3*k+0] = x1[ix];
	    xn[3*k+1] = x1[iy];
	    xn[3*k+2] = x1[iz];
	  }
    }
    break;
  case xfe_TriHierarch:
  case xfe_TetHierarch:
  case xfe_QuadLegendre:
  case xfe_HexLegendre:
    xf_printf("Error, calling LagrangeNodes function with a non-Lagrange basis.\n");
    return xf_Error(xf_INPUT_ERROR);
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }

  xf_Release( (void *) x1);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_LagrangeNodesEqual
int 
xf_LagrangeNodesEqual(enum xfe_ShapeType Shape, int p, real *xn){
  int ierr;
  enum xfe_BasisType Basis;
  
  switch (Shape){
  case xfe_Segment:
    Basis = xfe_SegLagrange;
    break;
  case xfe_Triangle:
    Basis = xfe_TriLagrange;
    break;
  case xfe_Quadrilateral:
    Basis = xfe_QuadLagrange;
    break;
  case xfe_Tetrahedron:
    Basis = xfe_TetLagrange;
    break;
  case xfe_Hexahedron:
    Basis = xfe_HexLagrange;
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }

  ierr = xf_Error(xf_LagrangeNodes(Basis, p, NULL, xn, NULL));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Basis2Lagrange
int 
xf_Basis2Lagrange(enum xfe_BasisType Basis, enum xfe_BasisType *BasisLag)
{  
  switch (Basis){
  case xfe_SegLagrange:
    (*BasisLag) = xfe_SegLagrange;
    break;
  case xfe_TriLagrange:
  case xfe_TriHierarch:
    (*BasisLag) = xfe_TriLagrange;
    break;
  case xfe_TetLagrange:
  case xfe_TetHierarch:
    (*BasisLag) = xfe_TetLagrange;
    break;
  case xfe_QuadLagrange: 
    (*BasisLag) = xfe_QuadLagrange;
    break;
  case xfe_QuadLagrangeGauss:
    (*BasisLag) = xfe_QuadLagrangeGauss;
    break;
  case xfe_HexLagrange:
    (*BasisLag) = xfe_HexLagrange;
    break;
  case xfe_HexLagrangeGauss:
    (*BasisLag) = xfe_HexLagrangeGauss;
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Basis2UniformLagrange
int 
xf_Basis2UniformLagrange(enum xfe_BasisType Basis, enum xfe_BasisType *BasisLag)
{  
  switch (Basis){
  case xfe_SegLagrange:
  case xfe_SegLagrangeGauss:
    (*BasisLag) = xfe_SegLagrange;
    break;
  case xfe_TriLagrange:
  case xfe_TriHierarch:
    (*BasisLag) = xfe_TriLagrange;
    break;
  case xfe_TetLagrange:
  case xfe_TetHierarch:
    (*BasisLag) = xfe_TetLagrange;
    break;
  case xfe_QuadLagrange: 
  case xfe_QuadLagrangeGauss:
  case xfe_QuadLegendre:
    (*BasisLag) = xfe_QuadLagrange;
    break;
  case xfe_HexLagrange:
  case xfe_HexLagrangeGauss:
  case xfe_HexLegendre:
    (*BasisLag) = xfe_HexLagrange;
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Shape2UniformLagrange
int 
xf_Shape2UniformLagrange(enum xfe_ShapeType Shape, enum xfe_BasisType *BasisLag)
{  
  switch (Shape){
  case xfe_Segment:
    (*BasisLag) = xfe_SegLagrange;
    break;
  case xfe_Triangle:
    (*BasisLag) = xfe_TriLagrange;
    break;
  case xfe_Quadrilateral:
    (*BasisLag) = xfe_QuadLagrange;
    break;
  case xfe_Tetrahedron:
    (*BasisLag) = xfe_TetLagrange;
    break;
  case xfe_Hexahedron:
    (*BasisLag) = xfe_HexLagrange;
    break;
  default:
    return xf_Error(xf_UNKNOWN_SHAPE);
    break;
  }  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateBasisData
static int 
xf_CreateBasisData( xf_BasisData **pBasisData)
{
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pBasisData, 1, sizeof(xf_BasisData)));
  if (ierr != xf_OK) return ierr;

  (*pBasisData)->InTable   = xfe_False;
  (*pBasisData)->Basis     = xfe_BasisLast;
  (*pBasisData)->Needs_ReCalc = xfe_False;
  (*pBasisData)->Order     = -1;
  (*pBasisData)->nn        = 0;
  (*pBasisData)->nq        = 0;
  (*pBasisData)->nnqmax     = 0;
  (*pBasisData)->dim       = 0;
  (*pBasisData)->Phi       = NULL;
  (*pBasisData)->GPhi      = NULL;
  (*pBasisData)->gPhi      = NULL;
  (*pBasisData)->HPhi      = NULL;
  (*pBasisData)->AllocFlag = xfb_None;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyBasisData
int 
xf_DestroyBasisData( xf_BasisData *BasisData, enum xfe_Bool ForceFlag)
{
  if (BasisData == NULL) return xf_OK;

  // do not destroy if this BasisData exists in a table
  if ((!ForceFlag) && BasisData->InTable) return xf_OK;

  xf_Release( (void *) BasisData->Phi);
  xf_Release( (void *) BasisData->GPhi);
  xf_Release( (void *) BasisData->gPhi);
  xf_Release( (void *) BasisData->HPhi);
  xf_Release( (void *) BasisData);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BasisLagrange1D
int
xf_BasisLagrange1D(real x, const real *xnode, int nnode, real *phi, real *gphi, 
		   real *hphi)
{
  // phi, gphi, hphi must be allocated prior to calling this function (if not NULL)
  int ierr, i, j, k, l;
  real pj, g, h;

  for (j=0; j<nnode; j++){
    if (phi != NULL){
      pj = 1.0;
      for (i=0  ; i<j    ; i++) pj *= (x-xnode[i])/(xnode[j]-xnode[i]);
      for (i=j+1; i<nnode; i++) pj *= (x-xnode[i])/(xnode[j]-xnode[i]);
      phi[j] = pj;
    }
    if (gphi != NULL){
      gphi[j] = 0.0;
      for (i=0  ; i<nnode; i++) 
	if (i != j){
	  g = 1.0/(xnode[j]-xnode[i]);
	  for (k=0; k<nnode; k++)
	    if ((k != i) && (k != j))
	      g *= (x-xnode[k])/(xnode[j]-xnode[k]);
	  gphi[j] += g;
	}
    }
    if (hphi != NULL){ // Hessian
      hphi[j] = 0.0;
      for (i=0  ; i<nnode; i++)
	if (i != j){
	  for (k=0; k<nnode; k++)
	    if ((k != i) && (k != j)){
	      h = 1.0/(xnode[j]-xnode[i]) * 1.0/(xnode[j]-xnode[k]);
	      for (l=0; l<nnode; l++)
		if ((l != i) && (l != j) && (l != k))
		  h *= (x-xnode[l])/(xnode[j]-xnode[l]);
	      hphi[j] += h;
	    }	  
	}
    }
  } // j
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BasisLegendre1D
int
xf_BasisLegendre1D(real x, const int p, real *phi, real *gphi,
                  real *hphi)
{
  // phi, gphi, hphi must be allocated prior to calling this function (if not NULL)
  int ierr;
  real xx;

  xx = x;

  if(phi != NULL){
     ierr = xf_Error(xf_Shape_LineLegendre(p, &xx, phi));
     if(ierr != xf_OK) return ierr;
  }

  if(gphi != NULL){
     ierr = xf_Error(xf_Grad_LineLegendre(p, &xx, gphi));
     if(ierr != xf_OK) return ierr;
  }

  if(hphi != NULL){
     ierr = xf_Error(xf_Hess_LineLegendre(p, &xx, hphi));
     if(ierr != xf_OK) return ierr;
  }

   return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BasisLagrange2D
static int
xf_BasisLagrange2D(const real *x, const real *xnode, int nnode, 
		   real *phi, real *gphi, real *hphi, int off)
{
  // phi and gphi must be allocated prior to calling this function (if not NULL)
  int ierr, ix, iy, ny;
  real py, gy, hy;
  real *phi1x, *phi1y;
  real *gphi1x, *gphi1y;
  real *hphi1x, *hphi1y;

  // allocate phi1x, phi1y
  ierr = xf_Error(xf_Alloc((void **) &phi1x, 2*nnode, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  phi1y = phi1x + nnode;

  // allocate gradients if necessary
  gphi1x = NULL;
  gphi1y = NULL;
  if ((gphi != NULL) || (hphi != NULL)){
    ierr = xf_Error(xf_Alloc((void **) &gphi1x, 2*nnode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    gphi1y = gphi1x + nnode;
  }

  // allocate Hessians if necessary
  hphi1x = NULL;
  hphi1y = NULL;
  if (hphi != NULL){
    ierr = xf_Error(xf_Alloc((void **) &hphi1x, 2*nnode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    hphi1y = hphi1x + nnode;
  }

  // call BasisLagrange1D with x -> phi1x
  ierr = xf_Error(xf_BasisLagrange1D(x[0], xnode, nnode, phi1x, gphi1x, hphi1x));
  if (ierr != xf_OK) return ierr;
  // call BasisLagrange1D with y -> phi1y
  ierr = xf_Error(xf_BasisLagrange1D(x[1], xnode, nnode, phi1y, gphi1y, hphi1y));
  if (ierr != xf_OK) return ierr;

  // set phi[iy*nn+ix] = phi1x[ix]*phi1y[iy]
  if (phi != NULL){
    for (iy=0; iy<nnode; iy++)
      for (ix=0, ny=iy*nnode, py = phi1y[iy]; ix<nnode; ix++)
 	phi[ny+ix] = phi1x[ix]*py;
  }

  if (gphi != NULL) // Gradients
    for (iy=0; iy<nnode; iy++)
      for (ix=0, ny=iy*nnode, py = phi1y[iy], gy = gphi1y[iy]; ix<nnode; ix++){
	gphi[ny+ix+  0] = gphi1x[ix]*py;
	gphi[ny+ix+off] = gy*phi1x[ix];
      }

  if (hphi != NULL) // Hessians
    for (iy=0; iy<nnode; iy++)
      for (ix=0, ny=iy*nnode, py=phi1y[iy], gy=gphi1y[iy], hy=hphi1y[iy]; ix<nnode; ix++){
	hphi[ny+ix+  0  ]                     = hphi1x[ix]*py;
	hphi[ny+ix+off  ] = hphi[ny+ix+off*2] = gy*gphi1x[ix];
	hphi[ny+ix+off*3]                     = hy*phi1x[ix];
      }

  xf_Release( (void *)  phi1x);
  xf_Release( (void *) gphi1x);
  xf_Release( (void *) hphi1x);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BasisLegendre2D
// added by Yu on Oct 2013
static int
xf_BasisLegendre2D(const real *x, const int p, real *phi, real *gphi,
                   real *hphi, int off)
{
    // phi and gphi must be allocated prior to calling this function (if not NULL)
    int ierr, ix, iy;
    int ny, nnode;
    real py, gy, hy;
    real *phi1x, *phi1y;
    real *gphi1x, *gphi1y;
    real *hphi1x, *hphi1y;
    
    nnode = p+1;
    
    // allocate phi1x, phi1y
    ierr = xf_Error(xf_Alloc((void **) &phi1x, 2*nnode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    phi1y = phi1x + nnode;
    
    // allocate gradients if necessary
    gphi1x = NULL;
    gphi1y = NULL;
    if ((gphi != NULL) || (hphi != NULL)){
        ierr = xf_Error(xf_Alloc((void **) &gphi1x, 2*nnode, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        gphi1y = gphi1x + nnode;
    }

    // allocate Hessians if necessary
    hphi1x = NULL;
    hphi1y = NULL;
    if (hphi != NULL){
        ierr = xf_Error(xf_Alloc((void **) &hphi1x, 2*nnode, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        hphi1y = hphi1x + nnode;
    }
    
    // call BasisLegendre1D with x -> phi1x
    ierr = xf_Error(xf_BasisLegendre1D(x[0], p, phi1x, gphi1x, hphi1x));
    if (ierr != xf_OK) return ierr;
    // call BasisLegendre1D with y -> phi1y
    ierr = xf_Error(xf_BasisLegendre1D(x[1], p, phi1y, gphi1y, hphi1y));
    if (ierr != xf_OK) return ierr;
    
    // set phi[iy*nn+ix] = phi1x[ix]*phi1y[iy]
    if (phi != NULL){
        for (iy=0; iy<nnode; iy++)
            for (ix=0, ny=iy*nnode, py = phi1y[iy]; ix<nnode; ix++)
                phi[ny+ix] = phi1x[ix]*py;
    }
    
    if (gphi != NULL) // Gradients
        for (iy=0; iy<nnode; iy++)
            for (ix=0, ny=iy*nnode, py = phi1y[iy], gy = gphi1y[iy]; ix<nnode; ix++){
                gphi[ny+ix+  0] = gphi1x[ix]*py;
                gphi[ny+ix+off] = gy*phi1x[ix];
            }
    
    if (hphi != NULL) // Hessians
        for (iy=0; iy<nnode; iy++)
            for (ix=0, ny=iy*nnode, py=phi1y[iy], gy=gphi1y[iy], hy=hphi1y[iy]; ix<nnode; ix++){
                hphi[ny+ix+  0  ]                     = hphi1x[ix]*py;
                hphi[ny+ix+off  ] = hphi[ny+ix+off*2] = gy*gphi1x[ix];
                hphi[ny+ix+off*3]                     = hy*phi1x[ix];
            }
    
    xf_Release( (void *)  phi1x);
    xf_Release( (void *) gphi1x);
    xf_Release( (void *) hphi1x);
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BasisLagrange3D
static int
xf_BasisLagrange3D(const real *x, real *xnode, int nnode, real *phi, 
		   real *gphi, real *hphi, int off)
{
  // phi and gphi must be allocated prior to calling this function (if not NULL)
  int ierr, ix, iy, ny, iz, nz, n;
  real py, gy, hy, pz, gz, hz; 
  real *phi1x, *phi1y, *phi1z;
  real *gphi1x, *gphi1y, *gphi1z;
  real *hphi1x, *hphi1y, *hphi1z;

  // allocate phi1x, ...
  ierr = xf_Error(xf_Alloc((void **) &phi1x, 3*nnode, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  phi1y = phi1x +   nnode;
  phi1z = phi1x + 2*nnode;

  // allocate gradients if necessary
  gphi1x = NULL;
  gphi1y = NULL;
  gphi1z = NULL;
  if ((gphi != NULL) || (hphi != NULL)) {
    ierr = xf_Error(xf_Alloc((void **) &gphi1x, 3*nnode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    gphi1y = gphi1x +   nnode;
    gphi1z = gphi1x + 2*nnode;
  }

  // allocate Hessians if necessary
  hphi1x = NULL;
  hphi1y = NULL;
  hphi1z = NULL;
  if (hphi != NULL){
    ierr = xf_Error(xf_Alloc((void **) &hphi1x, 3*nnode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    hphi1y = hphi1x +   nnode;
    hphi1z = hphi1x + 2*nnode;
  }

  // call BasisLagrange1D with x -> phi1x
  ierr = xf_Error(xf_BasisLagrange1D(x[0], xnode, nnode, phi1x, gphi1x, hphi1x));
  if (ierr != xf_OK) return ierr;
  // call BasisLagrange1D with y -> phi1y
  ierr = xf_Error(xf_BasisLagrange1D(x[1], xnode, nnode, phi1y, gphi1y, hphi1y));
  if (ierr != xf_OK) return ierr;
  // call BasisLagrange1D with z -> phi1z
  ierr = xf_Error(xf_BasisLagrange1D(x[2], xnode, nnode, phi1z, gphi1z, hphi1z));
  if (ierr != xf_OK) return ierr;


  // set phi[iy*nn+ix] = phi1x[ix]*phi1y[iy]*phi1z[iz]
  if (phi != NULL)
    for (iz=0; iz<nnode; iz++)
      for (iy=0, nz=iz*nnode*nnode, pz = phi1z[iz]; iy<nnode; iy++)
	for (ix=0, ny=iy*nnode, py = phi1y[iy]; ix<nnode; ix++)
	  phi[nz+ny+ix] = phi1x[ix]*py*pz;
  
  // Gradients
  if (gphi != NULL)
    for (iz=0; iz<nnode; iz++)
      for (iy=0, nz=iz*nnode*nnode, pz = phi1z[iz], gz = gphi1z[iz]; iy<nnode; iy++)
	for (ix=0, ny=iy*nnode, py = phi1y[iy], gy = gphi1y[iy]; ix<nnode; ix++){
	  n = nz+ny+ix;
	  gphi[n+    0] = gphi1x[ix]*py*pz;
	  gphi[n+  off] = phi1x[ix]*gy*pz;
	  gphi[n+2*off] = phi1x[ix]*py*gz;
	}
  
  // Hessians
  if (hphi != NULL)
    for (iz=0; iz<nnode; iz++)
      for (iy=0, nz=iz*nnode*nnode, pz=phi1z[iz], gz=gphi1z[iz], hz=hphi1z[iz]; iy<nnode; iy++)
	for (ix=0, ny=iy*nnode, py=phi1y[iy], gy=gphi1y[iy], hy=hphi1y[iy]; ix<nnode; ix++){
	  n = nz+ny+ix;
	  hphi[n+    0] =                 hphi1x[ix]*py*pz;
	  hphi[n+  off] = hphi[n+3*off] = gphi1x[ix]*gy*pz;
	  hphi[n+2*off] = hphi[n+6*off] = gphi1x[ix]*py*gz;
	  hphi[n+4*off] =                  phi1x[ix]*hy*pz;
	  hphi[n+5*off] = hphi[n+7*off] =  phi1x[ix]*gy*gz;
	  hphi[n+8*off] =                  phi1x[ix]*py*hz;
	}


  xf_Release( (void *)  phi1x);
  xf_Release( (void *) gphi1x);
  xf_Release( (void *) hphi1x);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BasisLegendre3D
// added by Yu on Oct 2013
static int
xf_BasisLegendre3D(const real *x, const int p, real *phi,
                   real *gphi, real *hphi, int off)
{
    // phi and gphi must be allocated prior to calling this function (if not NULL)
    int ierr, ix, iy, iz;
    int nnode, ny, nz, n;
    real py, gy, hy, pz, gz, hz;
    real *phi1x, *phi1y, *phi1z;
    real *gphi1x, *gphi1y, *gphi1z;
    real *hphi1x, *hphi1y, *hphi1z;
    
    nnode = p+1;
    
    // allocate phi1x, ...
    ierr = xf_Error(xf_Alloc((void **) &phi1x, 3*nnode, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    phi1y = phi1x +   nnode;
    phi1z = phi1x + 2*nnode;
    
    // allocate gradients if necessary
    gphi1x = NULL;
    gphi1y = NULL;
    gphi1z = NULL;
    if ((gphi != NULL) || (hphi != NULL)) {
        ierr = xf_Error(xf_Alloc((void **) &gphi1x, 3*nnode, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        gphi1y = gphi1x +   nnode;
        gphi1z = gphi1x + 2*nnode;
    }
    
    // allocate Hessians if necessary
    hphi1x = NULL;
    hphi1y = NULL;
    hphi1z = NULL;
    if (hphi != NULL){
        ierr = xf_Error(xf_Alloc((void **) &hphi1x, 3*nnode, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        hphi1y = hphi1x +   nnode;
        hphi1z = hphi1x + 2*nnode;
    }
    
    // call BasisLegendre1D with x -> phi1x
    ierr = xf_Error(xf_BasisLegendre1D(x[0], p, phi1x, gphi1x, hphi1x));
    if (ierr != xf_OK) return ierr;
    // call BasisLegendre1D with y -> phi1y
    ierr = xf_Error(xf_BasisLegendre1D(x[1], p, phi1y, gphi1y, hphi1y));
    if (ierr != xf_OK) return ierr;
    // call BasisLegendre1D with z -> phi1z
    ierr = xf_Error(xf_BasisLegendre1D(x[2], p, phi1z, gphi1z, hphi1z));
    if (ierr != xf_OK) return ierr;
    
    
    // set phi[iy*nn+ix] = phi1x[ix]*phi1y[iy]*phi1z[iz]
    if (phi != NULL)
        for (iz=0; iz<nnode; iz++)
            for (iy=0, nz=iz*nnode*nnode, pz = phi1z[iz]; iy<nnode; iy++)
                for (ix=0, ny=iy*nnode, py = phi1y[iy]; ix<nnode; ix++)
                    phi[nz+ny+ix] = phi1x[ix]*py*pz;
    
    // Gradients
    if (gphi != NULL)
        for (iz=0; iz<nnode; iz++)
            for (iy=0, nz=iz*nnode*nnode, pz = phi1z[iz], gz = gphi1z[iz]; iy<nnode; iy++)
                for (ix=0, ny=iy*nnode, py = phi1y[iy], gy = gphi1y[iy]; ix<nnode; ix++){
                    n = nz+ny+ix;
                    gphi[n+    0] = gphi1x[ix]*py*pz;
                    gphi[n+  off] = phi1x[ix]*gy*pz;
                    gphi[n+2*off] = phi1x[ix]*py*gz;
                }
    
    // Hessians
    if (hphi != NULL)
        for (iz=0; iz<nnode; iz++)
            for (iy=0, nz=iz*nnode*nnode, pz=phi1z[iz], gz=gphi1z[iz], hz=hphi1z[iz]; iy<nnode; iy++)
                for (ix=0, ny=iy*nnode, py=phi1y[iy], gy=gphi1y[iy], hy=hphi1y[iy]; ix<nnode; ix++){
                    n = nz+ny+ix;
                    hphi[n+    0] =                 hphi1x[ix]*py*pz;
                    hphi[n+  off] = hphi[n+3*off] = gphi1x[ix]*gy*pz;
                    hphi[n+2*off] = hphi[n+6*off] = gphi1x[ix]*py*gz;
                    hphi[n+4*off] =                  phi1x[ix]*hy*pz;
                    hphi[n+5*off] = hphi[n+7*off] =  phi1x[ix]*gy*gz;
                    hphi[n+8*off] =                  phi1x[ix]*py*hz;
                }
    
    
    xf_Release( (void *)  phi1x);
    xf_Release( (void *) gphi1x);
    xf_Release( (void *) hphi1x);
    
    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Shape_TensorLagrange
static int 
xf_Shape_TensorLagrange(int dim, int p, const real *x, real *phi)
{
  int ierr, nnode, i;
  real *xnode, dx;

  // return immediately if p == 0
  if (p == 0) { (*phi) = 1.; return xf_OK;}

  nnode = p+1; // for 1d

  // construct xnode
  ierr = xf_Error(xf_Alloc((void **) &xnode, nnode, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for (i=0, dx = 1./((real) p); i<nnode; i++) xnode[i] = i*dx;

  // call BasisLagrange
  if (dim == 1)
    ierr = xf_Error(xf_BasisLagrange1D(x[0], xnode, nnode, phi, NULL, NULL));
  else if (dim == 2)
    ierr = xf_Error(xf_BasisLagrange2D(x, xnode, nnode, phi, NULL, NULL, 0));
  else if (dim == 3)
    ierr = xf_Error(xf_BasisLagrange3D(x, xnode, nnode, phi, NULL, NULL, 0));
  else
    return xf_Error(xf_NOT_SUPPORTED);
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) xnode);
 
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Shape_TensorLegendre
// p: polynomial order
// added by Yu on Oct 2013
static int
xf_Shape_TensorLegendre(int dim, const int p, const real *x, real *phi)
{
    int ierr, i;
    
    // return immediately if p == 0
    if (p == 0) { (*phi) = 1.; return xf_OK;}
    
    // call BasisLegendre
    if (dim == 1)
        ierr = xf_Error(xf_BasisLegendre1D(x[0], p, phi, NULL, NULL));
    else if (dim == 2)
        ierr = xf_Error(xf_BasisLegendre2D(x, p, phi, NULL, NULL, 0));
    else if (dim == 3)
        ierr = xf_Error(xf_BasisLegendre3D(x, p, phi, NULL, NULL, 0));
    else
        return xf_Error(xf_NOT_SUPPORTED);
    if (ierr != xf_OK) return ierr;
    
    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Shape_TensorLagrangeGauss
static int 
xf_Shape_TensorLagrangeGauss(int dim, int p, const real *x, real *phi)
{
  int ierr, nnode, i;
  real *xnode;

  // return immediately if p == 0
  if (p == 0) {(*phi) = 1.; return xf_OK;}

  // construct xnode
  xnode = NULL;
  ierr = xf_Error(xf_QuadLine(2*p+1, &nnode, &xnode, NULL));
  if (ierr != xf_OK) return ierr;
  if (nnode != p+1) return xf_Error(xf_BASIS_FCN_ERROR);

  // call BasisLagrange
  if (dim == 1)
    ierr = xf_Error(xf_BasisLagrange1D(x[0], xnode, nnode, phi, NULL, NULL));
  else if (dim == 2)
    ierr = xf_Error(xf_BasisLagrange2D(x, xnode, nnode, phi, NULL, NULL, 0));
  else if (dim == 3)
    ierr = xf_Error(xf_BasisLagrange3D(x, xnode, nnode, phi, NULL, NULL, 0));
  else
    return xf_Error(xf_NOT_SUPPORTED);
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) xnode);
 
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Grad_TensorLagrange
int 
xf_Grad_TensorLagrange(int dim, int p, const real *x, real *gphi, int off)
{
  int ierr, nnode, i;
  real *xnode, dx;

  // return immediately if p == 0
  if (p == 0){
    gphi[0] = gphi[off] = 0.;
    if (dim == 3) gphi[2*off] = 0.;
    return xf_OK;
  }

  nnode = p+1; // for 1d

  // construct xnode
  ierr = xf_Error(xf_Alloc((void **) &xnode, nnode, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for (i=0, dx = 1./((real) p); i<nnode; i++) xnode[i] = i*dx;

  // call BasisLagrange
  if (dim == 1)
    ierr = xf_Error(xf_BasisLagrange1D(x[0], xnode, nnode, NULL, gphi, NULL));
  else if (dim == 2)
    ierr = xf_Error(xf_BasisLagrange2D(x, xnode, nnode, NULL, gphi, NULL, off));
  else if (dim == 3)
    ierr = xf_Error(xf_BasisLagrange3D(x, xnode, nnode, NULL, gphi, NULL, off));
  else
    return xf_Error(xf_NOT_SUPPORTED);
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) xnode);
 
  return xf_OK;

}

/******************************************************************/
//   FUNCTION Definition: xf_Grad_TensorLegendre
// added by Yu on Oct 2013
int
xf_Grad_TensorLegendre(int dim, int p, const real *x, real *gphi, int off)
{
    int ierr, i;
    
    // return immediately if p == 0
    if (p == 0){
        gphi[0] = gphi[off] = 0.;
        if (dim == 3) gphi[2*off] = 0.;
        return xf_OK;
    }
    
    // call BasisLegendre
    if (dim == 1)
        ierr = xf_Error(xf_BasisLegendre1D(x[0], p, NULL, gphi, NULL));
    else if (dim == 2)
        ierr = xf_Error(xf_BasisLegendre2D(x, p, NULL, gphi, NULL, off));
    else if (dim == 3)
        ierr = xf_Error(xf_BasisLegendre3D(x, p, NULL, gphi, NULL, off));
    else
        return xf_Error(xf_NOT_SUPPORTED);
    if (ierr != xf_OK) return ierr;
    
    return xf_OK;
    
}

/******************************************************************/
//   FUNCTION Definition: xf_Grad_TensorLagrangeGauss
static int 
xf_Grad_TensorLagrangeGauss(int dim, int p, const real *x, real *gphi, int off)
{
  int ierr, nnode, i;
  real *xnode;

  // return immediately if p == 0
  if (p == 0){
    gphi[0] = gphi[off] = 0.;
    if (dim == 3) gphi[2*off] = 0.;
    return xf_OK;
  }

  // construct xnode
  xnode = NULL;
  ierr = xf_Error(xf_QuadLine(2*p+1, &nnode, &xnode, NULL));
  if (ierr != xf_OK) return ierr;
  if (nnode != p+1) return xf_Error(xf_BASIS_FCN_ERROR);

  // call BasisLagrange
  if (dim == 1)
    ierr = xf_Error(xf_BasisLagrange1D(x[0], xnode, nnode, NULL, gphi, NULL));
  else if (dim == 2)
    ierr = xf_Error(xf_BasisLagrange2D(x, xnode, nnode, NULL, gphi, NULL, off));
  else if (dim == 3)
    ierr = xf_Error(xf_BasisLagrange3D(x, xnode, nnode, NULL, gphi, NULL, off));
  else
    return xf_Error(xf_NOT_SUPPORTED);
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) xnode);
 
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Hess_TensorLagrange
static int 
xf_Hess_TensorLagrange(int dim, int p, const real *x, real *hphi, int off)
{
  int ierr, nnode, i;
  real *xnode, dx;

  // return immediately if p == 0
  if (p == 0){
    for (i=0; i<dim*dim; i++) hphi[i*off] = 0.;
    return xf_OK;
  }

  nnode = p+1; // for 1d

  // construct xnode
  ierr = xf_Error(xf_Alloc((void **) &xnode, nnode, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for (i=0, dx = 1./((real) p); i<nnode; i++) xnode[i] = i*dx;

  // call BasisLagrange
  if (dim == 1)
    ierr = xf_Error(xf_BasisLagrange1D(x[0], xnode, nnode, NULL, NULL, hphi));
  else if (dim == 2)
    ierr = xf_Error(xf_BasisLagrange2D(x, xnode, nnode, NULL, NULL, hphi, off));
  else if (dim == 3)
    ierr = xf_Error(xf_BasisLagrange3D(x, xnode, nnode, NULL, NULL, hphi, off));
  else
    return xf_Error(xf_NOT_SUPPORTED);
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) xnode);
 
  return xf_OK;

}

/******************************************************************/
//   FUNCTION Definition: xf_Hess_TensorLegendre
// added by Yu on Oct 2013
static int
xf_Hess_TensorLegendre(int dim, int p, const real *x, real *hphi, int off)
{
    int ierr, i;
    
    // return immediately if p == 0
    if (p == 0){
        for (i=0; i<dim*dim; i++) hphi[i*off] = 0.;
        return xf_OK;
    }
    
    // call BasisLegendre
    if (dim == 1)
        ierr = xf_Error(xf_BasisLegendre1D(x[0], p, NULL, NULL, hphi));
    else if (dim == 2)
        ierr = xf_Error(xf_BasisLegendre2D(x, p, NULL, NULL, hphi, off));
    else if (dim == 3)
        ierr = xf_Error(xf_BasisLegendre3D(x, p, NULL, NULL, hphi, off));
    else
        return xf_Error(xf_NOT_SUPPORTED);
    if (ierr != xf_OK) return ierr;
    
    return xf_OK;
    
}

/******************************************************************/
//   FUNCTION Definition: xf_Hess_TensorLagrangeGauss
static int 
xf_Hess_TensorLagrangeGauss(int dim, int p, const real *x, real *hphi, int off)
{
  int ierr, nnode, i;
  real *xnode;

  // return immediately if p == 0
  if (p == 0){
    for (i=0; i<dim*dim; i++) hphi[i*off] = 0.;
    return xf_OK;
  }

  // construct xnode
  xnode = NULL;
  ierr = xf_Error(xf_QuadLine(2*p+1, &nnode, &xnode, NULL));
  if (ierr != xf_OK) return ierr;
  if (nnode != p+1) return xf_Error(xf_BASIS_FCN_ERROR);

  // call BasisLagrange
  if (dim == 1)
    ierr = xf_Error(xf_BasisLagrange1D(x[0], xnode, nnode, NULL, NULL, hphi));
  else if (dim == 2)
    ierr = xf_Error(xf_BasisLagrange2D(x, xnode, nnode, NULL, NULL, hphi, off));
  else if (dim == 3)
    ierr = xf_Error(xf_BasisLagrange3D(x, xnode, nnode, NULL, NULL, hphi, off));
  else
    return xf_Error(xf_NOT_SUPPORTED);
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) xnode);
 
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_GetShapes
// modified by Yu on Oct 2013
static int 
xf_GetShapes(enum xfe_BasisType Basis, int Order, int nq, real *xq, 
	     real *Phi){
  int ierr, iq, nn;

  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;

  ierr = 0;
  switch (Basis){
  case xfe_SegLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TensorLagrange(1, Order, xq+iq, Phi+nn*iq);
    break;
  case xfe_SegLagrangeGauss:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TensorLagrangeGauss(1, Order, xq+iq, Phi+nn*iq);
    break;
  case xfe_TriLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TriLagrange(Order, xq+2*iq, Phi+nn*iq);
    break;
  case xfe_TriHierarch:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TriHierarch(Order, xq+2*iq, Phi+nn*iq);
    break;
  case xfe_TetLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TetLagrange(Order, xq+3*iq, Phi+nn*iq);
    break;
  case xfe_TetHierarch:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TetHierarch(Order, xq+3*iq, Phi+nn*iq);
    break;
  case xfe_QuadLagrange: 
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TensorLagrange(2, Order, xq+2*iq, Phi+nn*iq);
    break;
  case xfe_QuadLagrangeGauss:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TensorLagrangeGauss(2, Order, xq+2*iq, Phi+nn*iq);
    break;
  case xfe_QuadLegendre:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TensorLegendre(2, Order, xq+2*iq, Phi+nn*iq);
    break;
  case xfe_HexLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TensorLagrange(3, Order, xq+3*iq, Phi+nn*iq);
    break;
  case xfe_HexLagrangeGauss:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TensorLagrangeGauss(3, Order, xq+3*iq, Phi+nn*iq);
    break;
  case xfe_HexLegendre:
    for (iq=0; iq<nq; iq++) ierr += xf_Shape_TensorLegendre(3, Order, xq+3*iq, Phi+nn*iq);
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }
  if (ierr != 0) return xf_Error(xf_BASIS_FCN_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetGrads
// modified by Yu on Oct 2013
int 
xf_GetGrads(enum xfe_BasisType Basis, int Order, int nq, real *xq, 
	    real *GPhi){
  int ierr, iq, nn, off;

  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;

  off = nn*nq; // offset in linear memory between grad wrt X, Y, Z

  ierr = 0;
  switch (Basis){
  case xfe_SegLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Grad_TensorLagrange(1, Order, xq+iq, GPhi+nn*iq, 0);
    break;
  case xfe_SegLagrangeGauss:
    for (iq=0; iq<nq; iq++) ierr += xf_Grad_TensorLagrangeGauss(1, Order, xq+iq, GPhi+nn*iq, 0);
    break;
  case xfe_TriLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Grad_TriLagrange(Order, xq+2*iq, GPhi+nn*iq, off);
    break;
  case xfe_TriHierarch:
    for (iq=0; iq<nq; iq++) ierr += xf_Grad_TriHierarch(Order, xq+2*iq, GPhi+nn*iq, off);
    break;
  case xfe_TetLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Grad_TetLagrange(Order, xq+3*iq, GPhi+nn*iq, off);
    break;
  case xfe_TetHierarch:
    for (iq=0; iq<nq; iq++) ierr += xf_Grad_TetHierarch(Order, xq+3*iq, GPhi+nn*iq, off);
    break;
  case xfe_QuadLagrange: 
    for (iq=0; iq<nq; iq++) 
      ierr += xf_Grad_TensorLagrange(2, Order, xq+2*iq, GPhi+nn*iq, off);
    break;
  case xfe_QuadLagrangeGauss:
    for (iq=0; iq<nq; iq++) 
      ierr += xf_Grad_TensorLagrangeGauss(2, Order, xq+2*iq, GPhi+nn*iq, off);
    break;
  case xfe_QuadLegendre:
    for (iq=0; iq<nq; iq++)
      ierr += xf_Grad_TensorLegendre(2, Order, xq+2*iq, GPhi+nn*iq, off);
    break;
  case xfe_HexLagrange:
    for (iq=0; iq<nq; iq++) 
      ierr += xf_Grad_TensorLagrange(3, Order, xq+3*iq, GPhi+nn*iq, off);
    break;
  case xfe_HexLagrangeGauss:
    for (iq=0; iq<nq; iq++) 
      ierr += xf_Grad_TensorLagrangeGauss(3, Order, xq+3*iq, GPhi+nn*iq, off);
    break;
  case xfe_HexLegendre:
    for (iq=0; iq<nq; iq++)
      ierr += xf_Grad_TensorLegendre(3, Order, xq+3*iq, GPhi+nn*iq, off);
    break;
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }
  if (ierr != 0) return xf_Error(xf_BASIS_FCN_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetHesss
// modifed by Yu on Oct 2013
static int 
xf_GetHesss(enum xfe_BasisType Basis, int Order, int nq, real *xq, 
	    real *HPhi){
  int ierr, iq, nn, off;

  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;

  off = nn*nq; // offset in linear memory between grad wrt X, Y, Z

  ierr = 0;
  switch (Basis){
  case xfe_SegLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Hess_TensorLagrange(1, Order, xq+iq, HPhi+nn*iq, 0);
    break;
  case xfe_TriLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Hess_TriLagrange(Order, xq+2*iq, HPhi+nn*iq, off);
    break;
  case xfe_TriHierarch:
    for (iq=0; iq<nq; iq++) ierr += xf_Hess_TriHierarch(Order, xq+2*iq, HPhi+nn*iq, off);
    break;
  case xfe_TetLagrange:
    for (iq=0; iq<nq; iq++) ierr += xf_Hess_TetLagrange(Order, xq+3*iq, HPhi+nn*iq, off);
    break;
  case xfe_TetHierarch:
    for (iq=0; iq<nq; iq++) ierr += xf_Hess_TetHierarch(Order, xq+3*iq, HPhi+nn*iq, off);
    break;
  case xfe_QuadLagrange: 
    for (iq=0; iq<nq; iq++) 
      ierr += xf_Hess_TensorLagrange(2, Order, xq+2*iq, HPhi+nn*iq, off);
    break;
  case xfe_QuadLagrangeGauss:
    for (iq=0; iq<nq; iq++) 
      ierr += xf_Hess_TensorLagrangeGauss(2, Order, xq+2*iq, HPhi+nn*iq, off);
    break;
  case xfe_QuadLegendre:
    for (iq=0; iq<nq; iq++)
      ierr += xf_Hess_TensorLegendre(2, Order, xq+2*iq, HPhi+nn*iq, off);
    break;
  case xfe_HexLagrange:
    for (iq=0; iq<nq; iq++) 
      ierr += xf_Hess_TensorLagrange(3, Order, xq+3*iq, HPhi+nn*iq, off);
    break;
  case xfe_HexLagrangeGauss:
    for (iq=0; iq<nq; iq++) 
      ierr += xf_Hess_TensorLagrangeGauss(3, Order, xq+3*iq, HPhi+nn*iq, off);
    break;
  case xfe_HexLegendre:
    for (iq=0; iq<nq; iq++)
      ierr += xf_Hess_TensorLegendre(3, Order, xq+3*iq, HPhi+nn*iq, off);
  default:
    return xf_Error(xf_UNKNOWN_BASIS);
    break;
  }
  if (ierr != 0) return xf_Error(xf_BASIS_FCN_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_NeedBasisReCalc
static enum xfe_Bool
xf_NeedBasisReCalc(enum xfe_BasisType Basis, int Order, enum xfe_Bool QuadChanged,
		   unsigned int AllocFlag, xf_BasisData *BD)
{
  enum xfe_Bool ReCalc;
  enum xfe_Bool  Want_Phi,  Want_GPhi,  Want_gPhi,  Want_HPhi;
  enum xfe_Bool Stale_Phi, Stale_GPhi, Stale_gPhi, Stale_HPhi;

  ReCalc = xfe_False;
  if (BD == NULL)
    ReCalc = xfe_True;
  else{
    if ((BD->Order != Order) || 
	(BD->Basis != Basis) ||
	(QuadChanged) ||
	(BD->Needs_ReCalc)) ReCalc = xfe_True;

    if (!ReCalc){
      Want_Phi  = (AllocFlag & xfb_Phi);  
      Want_GPhi = (AllocFlag & xfb_GPhi);
      Want_gPhi = (AllocFlag & xfb_gPhi);
      Want_HPhi = (AllocFlag & xfb_HPhi);
      
      Stale_Phi  = !(BD->AllocFlag & xfb_Phi);
      Stale_GPhi = !(BD->AllocFlag & xfb_GPhi);
      Stale_gPhi = !(BD->AllocFlag & xfb_gPhi);
      Stale_HPhi = !(BD->AllocFlag & xfb_HPhi);
      
      if (((Want_Phi  && Stale_Phi ) || 
	   (Want_GPhi && Stale_GPhi) ||
	   (Want_gPhi && Stale_gPhi) ||
	   (Want_HPhi && Stale_HPhi)))
	ReCalc = xfe_True;
    }
  }

  return ReCalc;
}


/******************************************************************/
//   FUNCTION Definition: xf_EvalBasis
int 
xf_EvalBasis(enum xfe_BasisType Basis, int Order, enum xfe_Bool QuadChanged,
	     int nq, real *xq, unsigned int AllocFlag, xf_BasisData **pBasisData)
{
  int ierr, nn, dim;
  enum xfe_Bool ReSize, ReCalc;
  enum xfe_Bool  Want_Phi,  Want_GPhi,  Want_gPhi,  Want_HPhi;
  enum xfe_Bool Stale_Phi, Stale_GPhi, Stale_gPhi, Stale_HPhi;
  xf_BasisData *BD;

  if ((*pBasisData) == NULL){
    ierr = xf_Error(xf_CreateBasisData(pBasisData));
    if (ierr != xf_OK) return ierr;
  }

  BD = (*pBasisData);

  ReCalc = xf_NeedBasisReCalc(Basis, Order, QuadChanged, AllocFlag, BD);

  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Basis2Dim(Basis, &dim));
  if (ierr != xf_OK) return ierr;
  
  // set Basis, Order, nn, and nq
  BD->Basis = Basis;
  BD->Order = Order;
  BD->nn = nn;
  BD->nq = nq;

  ReSize = xfe_False;
  if ((nn*nq > BD->nnqmax) || (dim != BD->dim)){
    // resizing will be needed
    BD->nnqmax = nn*nq;
    BD->dim   = dim;
    ReSize = xfe_True;
  }

  // sanity check
  if (ReSize && (!ReCalc)) return xf_Error(xf_CODE_LOGIC_ERROR);

  Want_Phi  = (AllocFlag & xfb_Phi);  
  Want_GPhi = (AllocFlag & xfb_GPhi);
  Want_gPhi = (AllocFlag & xfb_gPhi);
  Want_HPhi = (AllocFlag & xfb_HPhi);

  Stale_Phi  = !(BD->AllocFlag & xfb_Phi);
  Stale_GPhi = !(BD->AllocFlag & xfb_GPhi);
  Stale_gPhi = !(BD->AllocFlag & xfb_gPhi);
  Stale_HPhi = !(BD->AllocFlag & xfb_HPhi);

  // obtain basis + gradients; reallocate if necessary
  if (Want_Phi){
    if (ReSize || Stale_Phi){
      ierr = xf_Error(xf_ReAlloc((void **) &BD->Phi, BD->nnqmax, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    if (ReCalc || Stale_Phi){
      ierr = xf_Error(xf_GetShapes(Basis, Order, nq, xq, BD->Phi));
      if (ierr != xf_OK) return ierr;
    }
  }

  if (Want_GPhi){
    if (ReSize || Stale_GPhi){
      ierr = xf_Error(xf_ReAlloc((void **) &BD->GPhi, BD->dim*BD->nnqmax, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    if (ReCalc || Stale_GPhi){
      ierr = xf_Error(xf_GetGrads(Basis, Order, nq, xq, BD->GPhi));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  if (Want_gPhi && (ReSize || Stale_gPhi)){
    // just allocate here, since physical gradients will be element specific
    ierr = xf_Error(xf_ReAlloc((void **) &BD->gPhi, BD->dim*BD->nnqmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  // Hessians (second derivatives)
  if (Want_HPhi){
    if (ReSize || Stale_HPhi){
      ierr = xf_Error(xf_ReAlloc((void **) &BD->HPhi, BD->dim*BD->dim*BD->nnqmax, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    if (ReCalc || Stale_HPhi){
      ierr = xf_Error(xf_GetHesss(Basis, Order, nq, xq, BD->HPhi));
      if (ierr != xf_OK) return ierr;
    }
  }


  // set AllocFlag
  BD->AllocFlag = AllocFlag;

  // Just performed calculation, so set Needs_ReCalc to false
  if (BD->Needs_ReCalc) BD->Needs_ReCalc = xfe_False; 

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_EvalPhysicalGrad
int 
xf_EvalPhysicalGrad(xf_BasisData *BasisData, xf_JacobianData *JData)
{
  int nq, nn, dim, i0, d0, d1;

  /* error if GPhi, gPhi, iJ are Null or possibly out of date in terms
     of allocated size */
  if ((BasisData->GPhi == NULL) || (!(BasisData->AllocFlag & xfb_GPhi)) ||
      (BasisData->gPhi == NULL) || (!(BasisData->AllocFlag & xfb_gPhi)) ||
      (    JData->iJ   == NULL) || (!(    JData->AllocFlag & xfb_iJ  )))
    return xf_Error(xf_NULL_OR_STALE_POINTER);
  
  nq = JData->nq;
  if (BasisData->nq <= 0) return xf_Error(xf_OUT_OF_BOUNDS);
  if ((nq != 1) && (BasisData->nq != nq)) return xf_Error(xf_OUT_OF_BOUNDS);

  dim = JData->dim;

  if (dim != BasisData->dim) return xf_Error(xf_OUT_OF_BOUNDS);

  nn = BasisData->nn;

  if (JData->nq == 1)
    xf_MTxM_Set(JData->iJ, BasisData->GPhi, dim, dim, BasisData->nq*nn, BasisData->gPhi);
  else{
    for (d0=0; d0<dim; d0++){
      i0 = nq*nn*d0;
      xf_ColMult_Set(BasisData->GPhi, JData->iJ+d0, nq, nn, 
		     dim*dim, BasisData->gPhi+i0);
      for (d1=1; d1<dim; d1++)
	xf_ColMult_Add(BasisData->GPhi+nq*nn*d1, JData->iJ+dim*d1+d0, nq, nn, dim*dim, 
		       BasisData->gPhi+i0);
    } // d0
  }

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_CreateBasisTable
int 
xf_CreateBasisTable( xf_BasisTable **pBasisTable)
{
  int ierr;
  int Shape, face, faceorient;
  xf_BasisTable *BT;
  
  ierr = xf_Error(xf_Alloc((void **) pBasisTable, 1, sizeof(xf_BasisTable)));
  if (ierr != xf_OK) return ierr;

  BT = (*pBasisTable);
  
  BT->nShape = xfe_ShapeLast;

  ierr = xf_Error(xf_GetLocFaceMax(&BT->nface));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_GetLocFaceOrientMax(&BT->nfaceorient));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc3((void ****) &BT->BasisData, BT->nShape, BT->nface, 
			    BT->nfaceorient, sizeof(xf_BasisData *)));
  if (ierr != xf_OK) return ierr;

  for (Shape=0; Shape<BT->nShape; Shape++)
    for (face=0; face<BT->nface; face++)
      for (faceorient=0; faceorient<BT->nfaceorient; faceorient++)
	BT->BasisData[Shape][face][faceorient] = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ForceReCalcBasisTable
int 
xf_ForceReCalcBasisTable( xf_BasisTable *BasisTable)
{
  int ierr;
  int Shape, face, faceorient;
  xf_BasisData *BD;

  if (BasisTable == NULL) return xf_OK;
  
  for (Shape=0; Shape<BasisTable->nShape; Shape++)
    for (face=0; face<BasisTable->nface; face++)
      for (faceorient=0; faceorient<BasisTable->nfaceorient; faceorient++)
	if ( (BD = BasisTable->BasisData[Shape][face][faceorient]) != NULL)
	  BD->Needs_ReCalc = xfe_True;	 
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyBasisTable
int 
xf_DestroyBasisTable( xf_BasisTable *BasisTable)
{
  int ierr;
  int Shape, face, faceorient;
  xf_BasisData *BD;

  if (BasisTable == NULL) return xf_OK;
  
  for (Shape=0; Shape<BasisTable->nShape; Shape++)
    for (face=0; face<BasisTable->nface; face++)
      for (faceorient=0; faceorient<BasisTable->nfaceorient; faceorient++)
	if ( (BD = BasisTable->BasisData[Shape][face][faceorient]) != NULL){
	  ierr = xf_Error(xf_DestroyBasisData(BD, xfe_True));
	  if (ierr != xf_OK) return ierr;
	}
  
  xf_Release3( (void ***) BasisTable->BasisData);

  xf_Release( (void *) BasisTable);

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_EvalBasisOnFace
static int 
xf_EvalBasisOnFace(xf_Mesh *Mesh, int egrp, int elem, int face, 
		   int faceorient, int hang, enum xfe_BasisType Basis, 
		   int Order, enum xfe_Bool QuadChanged, int nq, real *xq, 
		   unsigned int AllocFlag, xf_BasisData **pBasisData,
		   real **pxelem)
{
  int ierr;
  real *xelem = NULL;
  enum xfe_Bool ReCalc;
  xf_BasisData *BD;

  // determine if need to recalculate BasisData
  ReCalc = xf_NeedBasisReCalc(Basis, Order, QuadChanged, AllocFlag, 
			      (*pBasisData));

  ReCalc = (ReCalc || (Mesh->ElemGroup[egrp].CutFlag || (hang != 0)));

  if (ReCalc){
    if (pxelem != NULL){
      ierr = xf_Error(xf_ReAlloc((void **) pxelem, nq*Mesh->Dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      xelem = (*pxelem);
    }
    else{
      ierr = xf_Error(xf_Alloc((void **) &xelem, nq*Mesh->Dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }

    // convert xq to xelem
    ierr = xf_Error(xf_RefFace2Interpol(Mesh, egrp, elem, face, faceorient,
					nq, xq, xelem));
    if (ierr != xf_OK) return ierr;

    // call EvalBasis (pass True for QuadChanged to ensure a recalc takes place)
    ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xelem, AllocFlag, 
				 pBasisData));
    if (ierr != xf_OK) return ierr;

    if (pxelem == NULL) xf_Release((void *) xelem);
  }
  else if (pxelem != NULL){
    // (*pxelem) must be calculated even if basisdata is in table
    ierr = xf_Error(xf_ReAlloc((void **) pxelem, nq*Mesh->Dim, sizeof(real)));
     if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_RefFace2Interpol(Mesh, egrp, elem, face, faceorient,
					nq, xq, (*pxelem)));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_EvalBasisOnFaceUsingTable
int 
xf_EvalBasisOnFaceUsingTable(xf_Mesh *Mesh, int egrp, int elem, int face, 
			     int faceorient, enum xfe_BasisType Basis, int Order, 
			     enum xfe_Bool QuadChanged, int nq, real *xq, 
			     unsigned int AllocFlag, xf_BasisData **pBasisData,
			     xf_BasisTable *BasisTable, real **pxelem)
{
  int ierr, hang;
  enum xfe_ShapeType Shape;

  if (Mesh->ElemGroup[egrp].CutFlag){ // do not use table for cut elements
    return xf_Error(xf_NOT_SUPPORTED);
/*     ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xelem, AllocFlag, pBasisData)); */
/*     if (ierr != xf_OK) return ierr; */
/*     (*pBasisData)->InTable = xfe_False; */
/*     return xf_OK; */
  }

  
  // if QuadChanged all stored bases need re-evaluating
  if (QuadChanged){
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
    ierr = xf_Error(xf_EvalBasisOnFace(Mesh, egrp, elem, face, faceorient, hang,
				       Basis, Order, QuadChanged, nq, xq, AllocFlag,
				       pBasisData, pxelem));
    if (ierr != xf_OK) return ierr;
    (*pBasisData)->InTable = xfe_False;
    return xf_OK;
  }

  // evaluate basis using table
  ierr = xf_Error(xf_EvalBasisOnFace(Mesh, egrp, elem, face, faceorient, 0,
				     Basis, Order, QuadChanged, nq, xq, AllocFlag, 
				     &BasisTable->BasisData[Shape][face][faceorient],
				     pxelem));
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
//   FUNCTION Definition: xf_ComputeGenMassMatrix
static int 
xf_ComputeGenMassMatrix(xf_Mesh *Mesh, int egrp, enum xfe_BasisType Basis1,
			int Order1, enum xfe_BasisType Basis2, int Order2,
			xf_Matrix *M){
  int ierr, nMatrix, i, j, iq, nn1, nn2;
  int pnq, nq, elem, QuadOrder;
  enum xfe_Bool GenericFlag, SameFlag, QuadChanged;
  real *MM, *xq, *wq, t, *phi1, *phi2;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData1;
  xf_BasisData *PhiData2;
  xf_JacobianData *JData;

  ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order1+Order2, &QuadOrder));
  if (ierr != xf_OK) return ierr;

  GenericFlag = (M->Linkage == xfe_LinkageElemGroup);
  SameFlag = ((Basis1 == Basis2) && (Order1 == Order2));

  QuadData = NULL;
  PhiData1 = NULL;
  PhiData2 = NULL;

  if (GenericFlag){
    MM = M->GenArray->rValue[0];  

    QuadData = NULL;
    ierr = xf_Error(xf_QuadElem(Mesh, egrp, 0, QuadOrder, &QuadData, &QuadChanged));
    if (ierr != xf_OK) return ierr;
    
    if (QuadData->Type != xfe_QuadDataGeneric) return xf_Error(xf_CODE_LOGIC_ERROR);
 
    nq = QuadData->nquad;
    xq = QuadData->xquad;
    wq = QuadData->wquad;

    // compute basis functions
    ierr = xf_Error(xf_EvalBasis(Basis1, Order1, QuadChanged, nq, xq, xfb_Phi, &PhiData1));
    if (ierr != xf_OK) return ierr;
    if (!SameFlag){
      ierr = xf_Error(xf_EvalBasis(Basis2, Order2, QuadChanged, nq, xq, xfb_Phi, &PhiData2));
      if (ierr != xf_OK) return ierr;
    }
    
    phi1 = phi2 = PhiData1->Phi;
    nn1  = nn2  = PhiData1->nn;
    if (!SameFlag){
      phi2 = PhiData2->Phi;
      nn2  = PhiData2->nn;
    }

    // form mass matrix
    for (i=0; i<nn1; i++)
      for (j=0; j<nn2; j++){
	t = 0.0;
	for (iq=0; iq<nq; iq++)
	  t += phi1[iq*nn1+i]*phi2[iq*nn2+j]*wq[iq];
	MM[i*nn2+j] = t;
      }
  }
  else{

    JData    = NULL;
    wq       = NULL;
    pnq = -1;
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      MM = M->GenArray->rValue[elem];
      /* Pull off quad points for the element; will not recalculate if
	 Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions
      ierr = xf_Error(xf_EvalBasis(Basis1, Order1, QuadChanged, nq, xq, xfb_Phi, &PhiData1));
      if (ierr != xf_OK) return ierr;
      if (!SameFlag){
	ierr = xf_Error(xf_EvalBasis(Basis2, Order2, QuadChanged, nq, xq, xfb_Phi, &PhiData2));
	if (ierr != xf_OK) return ierr;
      }

      /* Compute geometry Jacobian; if not constant, compute at quad
	 points.  Note if jacobian is constant, only one Jacobian will
	 be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;
      
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }
      
      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      phi1 = phi2 = PhiData1->Phi;
      nn1  = nn2  = PhiData1->nn;
      if (!SameFlag){
	phi2 = PhiData2->Phi;
	nn2  = PhiData2->nn;
      }

      // form mass matrix
      for (i=0; i<nn1; i++)
	for (j=0; j<nn2; j++){
	  t = 0.0;
	  for (iq=0; iq<nq; iq++)
	    t += phi1[iq*nn1+i]*phi2[iq*nn2+j]*wq[iq];
	  MM[i*nn2+j] = t;
	}
      
      pnq = nq;
    } // elem
      
    /* Destroy geometry Jacobian Data */
    ierr = xf_Error(xf_DestroyJacobianData(JData));
    if (ierr != xf_OK) return ierr;

    xf_Release((void *) wq);
  }

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(PhiData2, xfe_True));
  if (ierr != xf_OK) return ierr;
    
  
  return xf_OK;
}


/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_ComputeStiffMatrix */
/* static int  */
/* xf_ComputeStiffMatrix(xf_Mesh *Mesh, int egrp, enum xfe_BasisType Basis, */
/* 		      int Order, xf_Matrix *M){ */
/*   int ierr, nMatrix, i, j, iq, nn1, nn2; */
/*   int pnq, nq, elem, QuadOrder; */
/*   enum xfe_Bool GenericFlag, SameFlag, QuadChanged; */
/*   real *MM, *xq, *wq, t, *phi1, *phi2; */
/*   xf_QuadData *QuadData; */
/*   xf_BasisData *PhiData1; */
/*   xf_BasisData *PhiData2; */
/*   xf_JacobianData *JData; */

/*   ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order1+Order2, &QuadOrder)); */
/*   if (ierr != xf_OK) return ierr; */

/*   GenericFlag = (M->Linkage == xfe_LinkageElemGroup); */
/*   SameFlag = ((Basis1 == Basis2) && (Order1 == Order2)); */

/*   QuadData = NULL; */
/*   PhiData1 = NULL; */
/*   PhiData2 = NULL; */

/*   if (GenericFlag){ */
/*     MM = M->GenArray->rValue[0];   */

/*     QuadData = NULL; */
/*     ierr = xf_Error(xf_QuadElem(Mesh, egrp, 0, QuadOrder, &QuadData, &QuadChanged)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     if (QuadData->Type != xfe_QuadDataGeneric) return xf_Error(xf_CODE_LOGIC_ERROR); */
 
/*     nq = QuadData->nquad; */
/*     xq = QuadData->xquad; */
/*     wq = QuadData->wquad; */

/*     // compute basis functions */
/*     ierr = xf_Error(xf_EvalBasis(Basis1, Order1, QuadChanged, nq, xq, xfb_Phi, &PhiData1)); */
/*     if (ierr != xf_OK) return ierr; */
/*     if (!SameFlag){ */
/*       ierr = xf_Error(xf_EvalBasis(Basis2, Order2, QuadChanged, nq, xq, xfb_Phi, &PhiData2)); */
/*       if (ierr != xf_OK) return ierr; */
/*     } */
    
/*     phi1 = phi2 = PhiData1->Phi; */
/*     nn1  = nn2  = PhiData1->nn; */
/*     if (!SameFlag){ */
/*       phi2 = PhiData2->Phi; */
/*       nn2  = PhiData2->nn; */
/*     } */

/*     // form mass matrix */
/*     for (i=0; i<nn1; i++) */
/*       for (j=0; j<nn2; j++){ */
/* 	t = 0.0; */
/* 	for (iq=0; iq<nq; iq++) */
/* 	  t += phi1[iq*nn1+i]*phi2[iq*nn2+j]*wq[iq]; */
/* 	MM[i*nn2+j] = t; */
/*       } */
/*   } */
/*   else{ */
/*     JData    = NULL; */
/*     wq       = NULL; */
/*     pnq = -1; */
/*     for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){ */
/*       MM = M->GenArray->rValue[elem]; */
/*       /\* Pull off quad points for the element; will not recalculate if */
/* 	 Basis/Order have not changed. *\/ */
/*       ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged)); */
/*       if (ierr != xf_OK) return ierr; */
      
/*       nq = QuadData->nquad; */
/*       xq = QuadData->xquad; */
      
/*       // compute basis functions */
/*       ierr = xf_Error(xf_EvalBasis(Basis1, Order1, QuadChanged, nq, xq, xfb_Phi, &PhiData1)); */
/*       if (ierr != xf_OK) return ierr; */
/*       if (!SameFlag){ */
/* 	ierr = xf_Error(xf_EvalBasis(Basis2, Order2, QuadChanged, nq, xq, xfb_Phi, &PhiData2)); */
/* 	if (ierr != xf_OK) return ierr; */
/*       } */

/*       /\* Compute geometry Jacobian; if not constant, compute at quad */
/* 	 points.  Note if jacobian is constant, only one Jacobian will */
/* 	 be computed/returned. *\/ */
/*       ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, &JData)); */
/*       if (ierr != xf_OK) return ierr; */
      
/*       if (nq > pnq){ */
/* 	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real))); */
/* 	if (ierr != xf_OK) return ierr; */
/*       } */
      
/*       // form detJ-multiplied quad weight vector, wq */
/*       for (iq=0; iq<nq; iq++)  */
/* 	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)]; */
      
/*       phi1 = phi2 = PhiData1->Phi; */
/*       nn1  = nn2  = PhiData1->nn; */
/*       if (!SameFlag){ */
/* 	phi2 = PhiData2->Phi; */
/* 	nn2  = PhiData2->nn; */
/*       } */

/*       // form mass matrix */
/*       for (i=0; i<nn1; i++) */
/* 	for (j=0; j<nn2; j++){ */
/* 	  t = 0.0; */
/* 	  for (iq=0; iq<nq; iq++) */
/* 	    t += phi1[iq*nn1+i]*phi2[iq*nn2+j]*wq[iq]; */
/* 	  MM[i*nn2+j] = t; */
/* 	} */
      
/*       pnq = nq; */
/*     } // elem */
      
/*     /\* Destroy geometry Jacobian Data *\/ */
/*     ierr = xf_Error(xf_DestroyJacobianData(JData)); */
/*     if (ierr != xf_OK) return ierr; */

/*     xf_Release((void *) wq); */
/*   } */

/*   // Only destroy QuadData if points are generic */
/*   ierr = xf_Error(xf_DestroyGenericQuadData(QuadData)); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   /\* Destroy Basis Data *\/ */
/*   ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */
/*   ierr = xf_Error(xf_DestroyBasisData(PhiData2, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */
    
  
/*   return xf_OK; */
/* } */




/******************************************************************/
//   FUNCTION Definition: xf_FindMassMatrixData
int 
xf_FindMassMatrixData( xf_All *All, int egrp, enum xfe_BasisType Basis,
		       int Order, xf_Matrix **pM)
{
  int ierr, nMatrix, nn;
  enum xfe_Bool Found;
  enum xfe_Bool MeshIsLinear;
  enum xfe_LinkageType Linkage;
  const char Title[] = "MassMatrix";
  xf_Matrix *M;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // Determine whether mesh is forced linear
  MeshIsLinear = xfe_False;
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "MeshIsLinear", &MeshIsLinear);
  if ((ierr != xf_NOT_FOUND) && (ierr != xf_OK)) return ierr;

  // set Linkage and nMatrix depending on Q/cut-flag of egrp
  if ((Mesh->ElemGroup[egrp].QOrder > 1) || 
      (Mesh->ElemGroup[egrp].CutFlag) || 
      ((xf_Q1BasisNotLinear(Basis)) && (!MeshIsLinear)) ||
      ((xf_Q1BasisNotLinear(Mesh->ElemGroup[egrp].QBasis)) && (!MeshIsLinear)) ){
    Linkage = xfe_LinkageElem;
    nMatrix = Mesh->ElemGroup[egrp].nElem;
  }
  else{ // standard Q1
    Linkage = xfe_LinkageElemGroup;
    nMatrix = 1;
  }
  
  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_FindMatrix(All->DataSet, Title, Linkage, egrp,
				Order, Order, Basis, Basis, xfe_SizeReal,
				nMatrix, nn*nn, xfe_True, NULL, pM, &Found));
  if (ierr != xf_OK) return ierr;

  if (!Found){
    // matrix did not exist in DataSet; memory was allocated; fill in data
    M = (*pM);
    ierr = xf_Error(xf_ComputeGenMassMatrix(Mesh, egrp, Basis, Order, Basis, Order, M));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FindGenMassMatrixData
int 
xf_FindGenMassMatrixData( xf_All *All, int egrp, enum xfe_BasisType Basis1,
			  int Order1, enum xfe_BasisType Basis2, int Order2,
			  xf_Matrix **pM)
{
  int ierr, nMatrix, nn1, nn2;
  enum xfe_Bool Found;
  enum xfe_LinkageType Linkage;
  const char Title[] = "GenMassMatrix";
  xf_Matrix *M;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  // set Linkage and nMatrix depending on Q/cut-flag of egrp
  if ((Mesh->ElemGroup[egrp].QOrder > 1) || 
      (Mesh->ElemGroup[egrp].CutFlag) || 
      (xf_Q1BasisNotLinear(Basis1)) ||
      (xf_Q1BasisNotLinear(Basis2))){
    Linkage = xfe_LinkageElem;
    nMatrix = Mesh->ElemGroup[egrp].nElem;
  }
  else{ // standard Q1
    Linkage = xfe_LinkageElemGroup;
    nMatrix = 1;
  }
  
  ierr = xf_Error(xf_Order2nNode(Basis1, Order1, &nn1));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Order2nNode(Basis2, Order2, &nn2));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_FindMatrix(All->DataSet, Title, Linkage, egrp,
				Order1, Order2, Basis1, Basis2, xfe_SizeReal,
				nMatrix, nn1*nn2, xfe_True, NULL, pM, &Found));
  if (ierr != xf_OK) return ierr;

  if (!Found){
    // matrix did not exist in DataSet; memory was allocated; fill in data
    M = (*pM);
    ierr = xf_Error(xf_ComputeGenMassMatrix(Mesh, egrp, Basis1, Order1, 
					    Basis2, Order2, M));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindInvMassMatrixData
int 
xf_FindInvMassMatrixData( xf_All *All, int egrp, enum xfe_BasisType Basis,
			  int Order, xf_Matrix **piM)
{
  int ierr, nMatrix, nn;
  enum xfe_Bool Found;
  enum xfe_Bool MeshIsLinear;
  enum xfe_LinkageType Linkage;
  const char Title[] = "InvMassMatrix";
  xf_Matrix *M, *iM;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // Determine whether mesh is forced linear
  MeshIsLinear = xfe_False;
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "MeshIsLinear", &MeshIsLinear);
  if ((ierr != xf_NOT_FOUND) && (ierr != xf_OK)) return ierr;
  
  // set Linkage and nMatrix depending on Q/cut-flag of egrp
  if ((Mesh->ElemGroup[egrp].QOrder > 1) || 
      (Mesh->ElemGroup[egrp].CutFlag) ||
      ((xf_Q1BasisNotLinear(Basis)) && (!MeshIsLinear)) ||
      ((xf_Q1BasisNotLinear(Mesh->ElemGroup[egrp].QBasis)) && (!MeshIsLinear)) ){
    Linkage = xfe_LinkageElem;
    nMatrix = Mesh->ElemGroup[egrp].nElem;
  }
  else{ // standard Q1
    Linkage = xfe_LinkageElemGroup;
    nMatrix = 1;
  }
  
  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_FindMatrix(All->DataSet, Title, Linkage, egrp,
				Order, Order, Basis, Basis, xfe_SizeReal,
				nMatrix, nn*nn, xfe_True, NULL, piM, &Found));
  if (ierr != xf_OK) return ierr;

  if (!Found){
    // matrix did not exist in DataSet; memory was allocated; fill in data
    iM = (*piM);
    ierr = xf_Error(xf_FindMassMatrixData( All, egrp, Basis, Order, &M));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_MatrixInvert(M, nn, iM));
    if (ierr != xf_OK) return ierr;
  }


  return xf_OK;
}



/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_FindStiffMatrixData */
/* int  */
/* xf_FindStiffMatrixData( xf_All *All, int egrp, enum xfe_BasisType Basis, */
/* 		       int Order, xf_Matrix **pM) */
/* { */
/*   int ierr, nMatrix, nn; */
/*   enum xfe_Bool Found; */
/*   enum xfe_LinkageType Linkage; */
/*   const char Title[] = "StiffMatrix"; */
/*   xf_Matrix *M; */
/*   xf_Mesh *Mesh; */

/*   Mesh = All->Mesh; */
  
/*   // set Linkage and nMatrix depending on Q/cut-flag of egrp */
/*   if ((Mesh->ElemGroup[egrp].QOrder > 1) ||  */
/*       (Mesh->ElemGroup[egrp].CutFlag) ||  */
/*       (xf_Q1BasisNotLinear(Basis))){ */
/*     Linkage = xfe_LinkageElem; */
/*     nMatrix = Mesh->ElemGroup[egrp].nElem; */
/*   } */
/*   else{ // standard Q1 */
/*     Linkage = xfe_LinkageElemGroup; */
/*     nMatrix = 1; */
/*   } */
  
/*   ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn)); */
/*   if (ierr != xf_OK) return ierr; */

/*   ierr = xf_Error(xf_FindMatrix(All->DataSet, Title, Linkage, egrp, */
/* 				Order, Order, Basis, Basis, xfe_SizeReal, */
/* 				nMatrix, nn*nn, xfe_True, NULL, pM, &Found)); */
/*   if (ierr != xf_OK) return ierr; */

/*   if (!Found){ */
/*     // matrix did not exist in DataSet; memory was allocated; fill in data */
/*     M = (*pM); */
/*     ierr = xf_Error(xf_ComputeStiffMatrix(Mesh, egrp, Basis, Order, M)); */
/*     if (ierr != xf_OK) return ierr; */
/*   } */

/*   return xf_OK; */
/* } */


/******************************************************************/
//   FUNCTION Definition: xf_ElemMassMatrix
int 
xf_ElemMassMatrix( xf_All *All, int egrp, int elem, enum xfe_BasisType Basis,
		   int Order, xf_Matrix *Min, real *detJin, real **pMM, real *fac){
  int ierr;
  enum xfe_Bool GenericFlag;
  real detJ;
  xf_Matrix *M;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  if (Min == NULL){
    ierr = xf_Error(xf_FindMassMatrixData( All, egrp, Basis, Order, &M));
    if (ierr != xf_OK) return ierr;
  }
  else{
    M = Min;
  }

  GenericFlag = (M->Linkage == xfe_LinkageElemGroup);
  if (GenericFlag){
    if (detJin != NULL)
      detJ = (*detJin);
    else{
      // Get linear element volume
      ierr = xf_Error(xf_LinearElemJacobian(Mesh->ElemGroup[egrp].QBasis, 
					    Mesh->ElemGroup[egrp].QOrder,
					    Mesh->ElemGroup[egrp].Node[elem],
					    Mesh->Coord, &detJ));
      if (ierr != xf_OK) return ierr;
    }
    (*fac) = detJ;
    (*pMM) = M->GenArray->rValue[0];
  }
  else{
    (*fac) = 1.0;
    (*pMM) = M->GenArray->rValue[elem];
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemGenMassMatrix
int 
xf_ElemGenMassMatrix( xf_All *All, int egrp, int elem, enum xfe_BasisType Basis1,
		      int Order1, enum xfe_BasisType Basis2, int Order2,
		      xf_Matrix *Min, real *detJin, real **pMM, real *fac)
{
  int ierr;
  enum xfe_Bool GenericFlag;
  real detJ;
  xf_Matrix *M;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  if (Min == NULL){
    ierr = xf_Error(xf_FindGenMassMatrixData( All, egrp, Basis1, Order1, Basis2, 
					      Order2, &M));
    if (ierr != xf_OK) return ierr;
  }
  else{
    M = Min;
  }

  GenericFlag = (M->Linkage == xfe_LinkageElemGroup);
  if (GenericFlag){
    if (detJin != NULL)
      detJ = (*detJin);
    else{
      // Get linear element volume
      ierr = xf_Error(xf_LinearElemJacobian(Mesh->ElemGroup[egrp].QBasis, 
					    Mesh->ElemGroup[egrp].QOrder,
					    Mesh->ElemGroup[egrp].Node[elem],
					    Mesh->Coord, &detJ));
      if (ierr != xf_OK) return ierr;
    }
    (*fac) = detJ;
    (*pMM) = M->GenArray->rValue[0];
  }
  else{
    (*fac) = 1.0;
    (*pMM) = M->GenArray->rValue[elem];
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemInvMassMatrix
int 
xf_ElemInvMassMatrix( xf_All *All, int egrp, int elem, enum xfe_BasisType Basis,
		      int Order, xf_Matrix *iMin, real *detJin, real **piMM, 
		      real *fac){
  int ierr;
  enum xfe_Bool GenericFlag;
  real detJ;
  xf_Matrix *iM;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  if (iMin == NULL){
    ierr = xf_Error(xf_FindInvMassMatrixData( All, egrp, Basis, Order, &iM));
    if (ierr != xf_OK) return ierr;
  }
  else{
    iM = iMin;
  }

  GenericFlag = (iM->Linkage == xfe_LinkageElemGroup);
  if (GenericFlag){
    if (detJin != NULL)
      detJ = (*detJin);
    else{
      // Get linear element volume
      ierr = xf_Error(xf_LinearElemJacobian(Mesh->ElemGroup[egrp].QBasis, 
					    Mesh->ElemGroup[egrp].QOrder,
					    Mesh->ElemGroup[egrp].Node[elem],
					    Mesh->Coord, &detJ));
      if (ierr != xf_OK) return ierr;
    }
    (*fac) = 1.0/detJ;
    (*piMM) = iM->GenArray->rValue[0];
  }
  else{
    (*fac) = 1.0;
    (*piMM) = iM->GenArray->rValue[elem];
  }

  return xf_OK;
}



/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_ElemStiffMatrix */
/* int  */
/* xf_ElemStiffMatrix( xf_All *All, int egrp, int elem, enum xfe_BasisType Basis, */
/* 		    int Order, xf_Matrix *Min, real *detJin, real **pMM, real *fac) */
/* { */
/*   int ierr; */
/*   enum xfe_Bool GenericFlag; */
/*   real detJ; */
/*   xf_Matrix *M; */
/*   xf_Mesh *Mesh; */

/*   Mesh = All->Mesh; */
  
/*   if (Min == NULL){ */
/*     ierr = xf_Error(xf_FindStiffMatrixData( All, egrp, Basis, Order, &M)); */
/*     if (ierr != xf_OK) return ierr; */
/*   } */
/*   else{ */
/*     M = Min; */
/*   } */

/*   GenericFlag = (M->Linkage == xfe_LinkageElemGroup); */
/*   if (GenericFlag){ */
/*     if (detJin != NULL) */
/*       detJ = (*detJin); */
/*     else{ */
/*       // Get linear element volume */
/*       ierr = xf_Error(xf_LinearElemJacobian(Mesh->ElemGroup[egrp].QBasis,  */
/* 					    Mesh->ElemGroup[egrp].QOrder, */
/* 					    Mesh->ElemGroup[egrp].Node[elem], */
/* 					    Mesh->Coord, &detJ)); */
/*       if (ierr != xf_OK) return ierr; */
/*     } */
/*     (*fac) = detJ; */
/*     (*pMM) = M->GenArray->rValue[0]; */
/*   } */
/*   else{ */
/*     (*fac) = 1.0; */
/*     (*pMM) = M->GenArray->rValue[elem]; */
/*   } */

/*   return xf_OK; */
/* } */



/******************************************************************/
//   FUNCTION Definition: xf_ComputeTransferMatrix
int 
xf_ComputeTransferMatrix(enum xfe_BasisType Basis1, int Order1, 
			 enum xfe_BasisType Basis2, int Order2,
			 enum xfe_ShapeType EShape, int pos,  real *TT)
{
  int ierr, i, j, iq, n1, n2, nq, dim, *P;
  enum xfe_ShapeType Shape1, Shape2;
  real t, *M22, *xq, *wq, *phi1, *phi2;
  real *xq1, *xq2, *xq0 = NULL;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData1, *PhiData2;

  // determine interpolation Shape associated with Basis1
  ierr = xf_Error(xf_Basis2Shape(Basis1, &Shape1));
  if (ierr != xf_OK) return ierr;

  // determine interpolation Shape associated with Basis2
  ierr = xf_Error(xf_Basis2Shape(Basis2, &Shape2));
  if (ierr != xf_OK) return ierr;

  // verify that Shape1 == Shape2
  if (Shape1 != Shape2) return xf_Error(xf_INPUT_ERROR);

  // Create QuadData
  ierr = xf_Error(xf_CreateQuadData( &QuadData ));
  if (ierr != xf_OK) return ierr;

  // set Type, Shape, and Order
  QuadData->Type  = xfe_QuadDataGeneric;
  QuadData->Shape = Shape1;
  QuadData->Order = max(Order1 + Order2, 2*Order2);

  // obtain generic points
  ierr = xf_Error(xf_GetQuadPoints(QuadData));
  if (ierr != xf_OK) return ierr;
		       
  nq = QuadData->nquad;
  xq = QuadData->xquad;
  wq = QuadData->wquad;

  // if 2 is a hanging node refinement, convert xq to as-seen-by-parent=1
  if (pos != 0){
    dim = QuadData->qdim;
    ierr = xf_Error(xf_Alloc((void **) &xq0, nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (iq=0; iq<nq*dim; iq++) xq0[iq] = xq[iq];

    // xq0 is scaled to as seen by parent
    ierr = xf_Error(xf_ScaleHangInterpol(EShape, pos, nq, xq0));
    if (ierr != xf_OK) return ierr; 

    xq1 = xq0;
    xq2 = xq;
  }
  else{
    xq1 = xq;
    xq2 = xq;
  }
  
  // compute basis functions for Basis1
  PhiData1  = NULL;
  ierr = xf_Error(xf_EvalBasis(Basis1, Order1, xfe_True, nq, xq1, xfb_Phi, &PhiData1));
  if (ierr != xf_OK) return ierr;

  phi1 = PhiData1->Phi;
  n1   = PhiData1->nn;
  
  // compute basis functions for Basis2
  PhiData2  = NULL;
  ierr = xf_Error(xf_EvalBasis(Basis2, Order2, xfe_True, nq, xq2, xfb_Phi, &PhiData2));
  if (ierr != xf_OK) return ierr;

  phi2 = PhiData2->Phi;
  n2   = PhiData2->nn;  
  
  // form M22 = Basis2 mass matrix
  ierr = xf_Error(xf_Alloc((void **) &M22, n2*n2, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<n2; i++)
    for (j=0; j<n2; j++){
      t = 0.0;
      for (iq=0; iq<nq; iq++)
	t += phi2[iq*n2+i]*phi2[iq*n2+j]*wq[iq];
      M22[i*n2+j] = t;
    }

  // form TT = Basis1*Basis2 matrix 
  for (i=0; i<n2; i++)
    for (j=0; j<n1; j++){
      t = 0.0;
      for (iq=0; iq<nq; iq++)
	t += phi2[iq*n2+i]*phi1[iq*n1+j]*wq[iq];
      TT[i*n1+j] = t;
    }

  // set TT = inv(M22) * TT
  ierr = xf_Error(xf_Alloc( (void **) &P, n2, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ComputePLU(M22, n2, P));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SolvePLU_Matrix(M22, P, n2, n1, TT));
  if (ierr != xf_OK) return ierr;


  // Destroy QuadData
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  // Destroy Basis Data
  ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(PhiData2, xfe_True));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) M22);
  xf_Release( (void *) P);
  xf_Release( (void *) xq0);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_FindTransferMatrix
int 
xf_FindTransferMatrix( xf_DataSet *DataSet, enum xfe_BasisType Basis1,
		       int Order1, enum xfe_BasisType Basis2, 
		       int Order2, xf_Matrix **pT)
{
  int ierr, n1, n2, nMatrix;
  enum xfe_Bool Found;
  enum xfe_LinkageType Linkage;
  const char Title[] = "TransferMatrix";
  real *TT;
  xf_Matrix *T;

  // Transfer matrices have no linkage and only one matrix
  Linkage = xfe_LinkageNone;
  nMatrix = 1;
  
  // determine n1
  ierr = xf_Error(xf_Order2nNode(Basis1, Order1, &n1));
  if (ierr != xf_OK) return ierr;

  // determine n2
  ierr = xf_Error(xf_Order2nNode(Basis2, Order2, &n2));
  if (ierr != xf_OK) return ierr;

  // destroy (*pT) if DataSet is NULL in case creating a new one
  if (DataSet == NULL){
    ierr = xf_Error(xf_DestroyMatrix((*pT)));
    if (ierr != xf_OK) return ierr;
  }
  
  ierr = xf_Error(xf_FindMatrix(DataSet, Title, Linkage, 0, Order1, Order2, 
				Basis1, Basis2, xfe_SizeReal, nMatrix, n1*n2, 
				(DataSet != NULL), NULL, pT, &Found));
  if (ierr != xf_OK) return ierr;

  if (!Found){
    // matrix did not exist in DataSet; memory was allocated; fill in data
    T = (*pT);
    ierr = xf_Error(xf_ComputeTransferMatrix(Basis1, Order1, Basis2, Order2,
					     0, 0, T->GenArray->rValue[0]));
    if (ierr != xf_OK) return ierr;
  }


  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ProjectOnElemQR
int 
xf_ProjectOnElemQR(xf_Mesh *Mesh, int egrp, int elem, int sr, 
		   xf_QuadData *QuadData, enum xfe_BasisType Basis, 
		   int Order, enum xfe_Bool Regenerate, real **pQ, 
		   real **pR, int **pP, real *uin, real *U)
{
  int ierr, k;
  int iq, nq, nn;
  int *P, *P0 = NULL;
  enum xfe_Bool QuadChanged;
  real *u, *xq;
  real *Q, *R, *x, *Q0 = NULL, *R0 = NULL;
  xf_BasisData *PhiData;
  
  Q = ( (pQ == NULL) ? Q0 : (*pQ));
  R = ( (pR == NULL) ? R0 : (*pR));
  P = ( (pP == NULL) ? P0 : (*pP));

  Regenerate = ((Regenerate) || (Q == NULL) || (R == NULL) || (P == NULL));

  // determine nn = # unknowns for elements in this group
  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;
  
  nq = QuadData->nquad;
  xq = QuadData->xquad;
  
  // need at least as many quad points as nn = # basis functions
  if (nq < nn) return xf_Error(xf_INPUT_ERROR);

  PhiData     = NULL;
 
  // allocate data
  ierr = xf_Error(xf_Alloc( (void **) &u, nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &x, nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  if (Regenerate){
    // compute basis functions; always recompute (see weighting note below)
    ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xq, xfb_Phi, &PhiData));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReAlloc( (void **) &Q, nq*nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &R, nq*nn, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &P, nn, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    if (nn != PhiData->nn) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    /* weight Phi by wq (note, this changes PhiData, hence we
       force a recompute every time above) */
    xf_ColMult(PhiData->Phi, QuadData->wquad, nq, nn, 1);

    // QR factor PhiData->Phi
    ierr = xf_Error(xf_QRFactorHouseholder(PhiData->Phi, nq, nn, Q, R));
    if (ierr != xf_OK) return ierr;

    // Invert R
    ierr = xf_Error(xf_ComputePLU(R, nn, P));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
    
  // set u = uin
  for (k=0; k<nq*sr; k++) u[k] = uin[k];

  // weight u by wq
  xf_ColMult(u, QuadData->wquad, nq, sr, 1);
  
  // Set x = Q^T * u
  xf_MTxM_Set(Q, u, nq, nq, sr, x);
  
  // Set x = R^{-1} * x
  ierr = xf_Error(xf_SolvePLU_Matrix(R, P, nn, sr, x));
  if (ierr != xf_OK) return ierr;
  
  // have to do it this way until a compact QR function is written
  for (k=0; k<nn*sr; k++) U[k] = x[k];

  xf_Release((void *) u);
  xf_Release((void *) x);
  if (pQ == NULL) xf_Release((void *) Q); else (*pQ) = Q;
  if (pR == NULL) xf_Release((void *) R); else (*pR) = R;
  if (pP == NULL) xf_Release((void *) P); else (*pP) = P;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ProjectOnElemQR_OnlyEnergy
// expect to improve the affect of least-square fitting
int 
xf_ProjectOnElemQR_OnlyEnergy(xf_Mesh *Mesh, int egrp, int elem, int sr, 
                              xf_QuadData *QuadData, enum xfe_BasisType Basis, 
                              int Order, enum xfe_Bool Regenerate, real **pQ, 
                              real **pR, int **pP, real *uin, real *Uout)
{
    int ierr, k;
    int iq, nq, nn;
    int *P, *P0 = NULL;
    enum xfe_Bool QuadChanged;
    real *u, *xq, *U;
    real *Q, *R, *x, *Q0 = NULL, *R0 = NULL;
    xf_BasisData *PhiData;
    
    Q = ( (pQ == NULL) ? Q0 : (*pQ));
    R = ( (pR == NULL) ? R0 : (*pR));
    P = ( (pP == NULL) ? P0 : (*pP));
    
    Regenerate = ((Regenerate) || (Q == NULL) || (R == NULL) || (P == NULL));
    
    // determine nn = # unknowns for elements in this group
    ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
    if (ierr != xf_OK) return ierr;
    
    nq = QuadData->nquad;
    xq = QuadData->xquad;
    
    // need at least as many quad points as nn = # basis functions
    if (nq < nn) return xf_Error(xf_INPUT_ERROR);
    
    PhiData     = NULL;
    
    // allocate data
    ierr = xf_Error(xf_Alloc( (void **) &u, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &x, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &U, nn*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    if (Regenerate){
        // compute basis functions; always recompute (see weighting note below)
        ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xq, xfb_Phi, &PhiData));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ReAlloc( (void **) &Q, nq*nq, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &R, nq*nn, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &P, nn, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        if (nn != PhiData->nn) return xf_Error(xf_CODE_LOGIC_ERROR);
        
        /* weight Phi by wq (note, this changes PhiData, hence we
         force a recompute every time above) */
        xf_ColMult(PhiData->Phi, QuadData->wquad, nq, nn, 1);
        
        // QR factor PhiData->Phi
        ierr = xf_Error(xf_QRFactorHouseholder(PhiData->Phi, nq, nn, Q, R));
        if (ierr != xf_OK) return ierr;
        
        // Invert R
        ierr = xf_Error(xf_ComputePLU(R, nn, P));
        if (ierr != xf_OK) return ierr;
        
        /* Destroy Basis Data */
        ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
        if (ierr != xf_OK) return ierr;
    }
    
    // set u = uin
    for (k=0; k<nq*sr; k++) u[k] = uin[k];
    
    // weight u by wq
    xf_ColMult(u, QuadData->wquad, nq, sr, 1);
    
    // Set x = Q^T * u
    xf_MTxM_Set(Q, u, nq, nq, sr, x);
    
    // Set x = R^{-1} * x
    ierr = xf_Error(xf_SolvePLU_Matrix(R, P, nn, sr, x));
    if (ierr != xf_OK) return ierr;
    
    // have to do it this way until a compact QR function is written
    for (k=0; k<nn*sr; k++) U[k] = x[k];
 

    //only replace the energy polynomial
    for (k=0; k<nn; k++) Uout[k*sr+3] = U[k*sr+3];
    
    xf_Release((void *) u);
    xf_Release((void *) U);
    xf_Release((void *) x);
    if (pQ == NULL) xf_Release((void *) Q); else (*pQ) = Q;
    if (pR == NULL) xf_Release((void *) R); else (*pR) = R;
    if (pP == NULL) xf_Release((void *) P); else (*pP) = P;
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ProjectOnElemLeastSquares
int 
xf_ProjectOnElemLeastSquares(xf_Mesh *Mesh, int egrp, int elem, int sr, 
			     xf_QuadData *QuadData, enum xfe_BasisType Basis, 
			     int Order, enum xfe_Bool Regenerate, real **pQ, 
			     real **pR, int **pP, real *uin, real *U)
{
  int ierr, k;
  int iq, nq, nn;
  int *P, *P0 = NULL;
  enum xfe_Bool QuadChanged;
  real *u, *xq;
  real *Q, *R, *Q0 = NULL, *R0 = NULL;
  xf_BasisData *PhiData;
  
  Q = ( (pQ == NULL) ? Q0 : (*pQ));
  R = ( (pR == NULL) ? R0 : (*pR));
  P = ( (pP == NULL) ? P0 : (*pP));

  Regenerate = ((Regenerate) || (Q == NULL) || (R == NULL) || (P == NULL));

  // determine nn = # unknowns for elements in this group
  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
  if (ierr != xf_OK) return ierr;
  
  nq = QuadData->nquad;
  xq = QuadData->xquad;
  
  // need at least as many quad points as nn = # basis functions
  if (nq < nn) return xf_Error(xf_INPUT_ERROR);

  PhiData     = NULL;
 
  // allocate data
  ierr = xf_Error(xf_Alloc( (void **) &u, nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  if (Regenerate){
    // compute basis functions; always recompute (see weighting note below)
    ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xq, xfb_Phi, &PhiData));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReAlloc( (void **) &Q, nq*nn, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &R, nn*nn, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &P, nn, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    if (nn != PhiData->nn) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    // store Phi->Data in Q
    for (k=0; k<nq*nn; k++) Q[k] = PhiData->Phi[k];
    
    // multiply columns of Q by quad weights
    xf_ColMult(Q, QuadData->wquad, nq, nn, 1);

    // Compute PhiData->Phi^T*Q and store in R
    xf_MTxM(PhiData->Phi, Q, nn, nq, nn, xfe_Set, R);

    // Invert R
    ierr = xf_Error(xf_ComputePLU(R, nn, P));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
    
  // set u = uin
  for (k=0; k<nq*sr; k++) u[k] = uin[k];

  // Set U = Q^T * u
  xf_MTxM_Set(Q, u, nn, nq, sr, U);
  
  // Set U = R^{-1} * U
  ierr = xf_Error(xf_SolvePLU_Matrix(R, P, nn, sr, U));
  if (ierr != xf_OK) return ierr;
    
  xf_Release((void *) u);
  if (pQ == NULL) xf_Release((void *) Q); else (*pQ) = Q;
  if (pR == NULL) xf_Release((void *) R); else (*pR) = R;
  if (pP == NULL) xf_Release((void *) P); else (*pP) = P;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ProjectOnElemLeastSquares_OnlyEnergy
// Only Project Modified Energy State onto Element
int 
xf_ProjectOnElemLeastSquares_OnlyEnergy(xf_Mesh *Mesh, int egrp, int elem, int sr, 
                                        xf_QuadData *QuadData, enum xfe_BasisType Basis, 
                                        int Order, enum xfe_Bool Regenerate, real **pQ, 
                                        real **pR, int **pP, real *uin, real *Uout)
{
    int ierr, k;
    int iq, nq, nn;
    int *P, *P0 = NULL;
    enum xfe_Bool QuadChanged;
    real *u, *xq, *U;
    real *Q, *R, *Q0 = NULL, *R0 = NULL;
    xf_BasisData *PhiData;
    
    Q = ( (pQ == NULL) ? Q0 : (*pQ));
    R = ( (pR == NULL) ? R0 : (*pR));
    P = ( (pP == NULL) ? P0 : (*pP));
    
    Regenerate = ((Regenerate) || (Q == NULL) || (R == NULL) || (P == NULL));
    
    // determine nn = # unknowns for elements in this group
    ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
    if (ierr != xf_OK) return ierr;
    
    nq = QuadData->nquad;
    xq = QuadData->xquad;
    
    // need at least as many quad points as nn = # basis functions
    if (nq < nn) return xf_Error(xf_INPUT_ERROR);
    
    PhiData     = NULL;
    
    // allocate data
    ierr = xf_Error(xf_Alloc( (void **) &u, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &U, nn*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    if (Regenerate){
        // compute basis functions; always recompute (see weighting note below)
        ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xq, xfb_Phi, &PhiData));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ReAlloc( (void **) &Q, nq*nn, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &R, nn*nn, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &P, nn, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        if (nn != PhiData->nn) return xf_Error(xf_CODE_LOGIC_ERROR);
        
        // store Phi->Data in Q
        for (k=0; k<nq*nn; k++) Q[k] = PhiData->Phi[k];
        
        // multiply columns of Q by quad weights
        xf_ColMult(Q, QuadData->wquad, nq, nn, 1);
        
        // Compute PhiData->Phi^T*Q and store in R
        xf_MTxM(PhiData->Phi, Q, nn, nq, nn, xfe_Set, R);
        
        // Invert R
        ierr = xf_Error(xf_ComputePLU(R, nn, P));
        if (ierr != xf_OK) return ierr;
        
        /* Destroy Basis Data */
        ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
        if (ierr != xf_OK) return ierr;
    }
    
    // set u = uin
    for (k=0; k<nq*sr; k++) u[k] = uin[k];
    
    // Set U = Q^T * u
    xf_MTxM_Set(Q, u, nn, nq, sr, U);
    
    // Set U = R^{-1} * U
    ierr = xf_Error(xf_SolvePLU_Matrix(R, P, nn, sr, U));
    if (ierr != xf_OK) return ierr;
    
    //only replace the energy polynomial
    for (k=0; k<nn; k++) Uout[k*sr+3] = U[k*sr+3];
    
    xf_Release((void *) u);
    xf_Release((void *) U);
    if (pQ == NULL) xf_Release((void *) Q); else (*pQ) = Q;
    if (pR == NULL) xf_Release((void *) R); else (*pR) = R;
    if (pP == NULL) xf_Release((void *) P); else (*pP) = P;
    
    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ShapeIsoFace
int 
xf_ShapeIsoFace(xf_All *All, int egrp, int elem, int face, int nq, real *xq, real *phi)
{
  int ierr, iq;
  enum xfe_BasisType QBasis;
  real x, y;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  QBasis = Mesh->ElemGroup[egrp].QBasis;

  switch (QBasis) {
    case xfe_SegLagrange:
      for (iq=0; iq<nq; iq++){
        x = xq[iq];
        if (face == 0) phi[iq] = 1.-x;
        else if (face == 1) phi[iq] = x;
        else return xf_Error(xf_INPUT_ERROR);
      }
      break;
    case xfe_TriLagrange:
      for (iq=0; iq<nq; iq++){
        x = xq[2*iq+0];
        y = xq[2*iq+1];
       if (face == 0) phi[iq] = 4.*x*y;
       else if (face == 1) phi[iq] = 4.*(1.-x-y)*y;
       else if (face == 2) phi[iq] = 4.*(1.-x-y)*x;
       else return xf_Error(xf_INPUT_ERROR);
      }
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ShapeCentroid
int 
xf_ShapeCentroid(enum xfe_ShapeType Shape, real *xref)
{
  int ierr, d, dim;

  // get dim
  ierr = xf_Error(xf_Shape2Dim(Shape, &dim));
  if (ierr != xf_OK) return ierr;

  // initialize xref to centroid
  switch (Shape){
  case xfe_Segment:
  case xfe_Quadrilateral:
  case xfe_Hexahedron: 
    for (d=0; d<dim; d++) xref[d] = 0.5; break;
  case xfe_Triangle:
  case xfe_Tetrahedron:
    for (d=0; d<dim; d++) xref[d] = 1.0/(dim+1.0); break;
  default: return xf_Error(xf_UNKNOWN_SHAPE); break;
  }

  return xf_OK;
}


#if( UNIT_TEST==1 )
#include "xf_Basis.test.in"
#endif

