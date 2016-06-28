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


#include "xf_Unit.h"
#include "xf_Mesh.h"
#include "xf_All.h"


TEST_xf_ElemRefEffectOnFace()
{
  int ierr, ref;

  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Quadrilateral, xfe_QuadRefUniform, 0, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_SegRefUniform);

  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Quadrilateral, xfe_QuadRefHoriz, 1, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_SegRefUniform);

  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Quadrilateral, xfe_QuadRefVert, 3, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_SegRefNone);

  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Hexahedron, xfe_HexRefUniform, 4, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_QuadRefUniform);

  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Hexahedron, xfe_HexRefSliceXY, 1, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_QuadRefVert);
  
  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Hexahedron, xfe_HexRefSliceXZ, 0, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_QuadRefHoriz);
  
  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Hexahedron, xfe_HexRefSliceYZ, 2, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_QuadRefUniform);
  
  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Hexahedron, xfe_HexRefSliceZ, 5, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_QuadRefNone);
  
  ierr = xf_Error(xf_ElemRefEffectOnFace(xfe_Hexahedron, xfe_HexRefSliceY, 4, &ref));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ref, xfe_QuadRefVert);

  return xf_OK;  
}


TEST_xf_FaceRefEffectOnElem()
{
  int ierr, ref;
  int frv0[] = {xfe_SegRefUniform, 0, xfe_SegRefUniform, 0};
  int frv1[] = {xfe_SegRefUniform, xfe_SegRefUniform, 0, 0};
  int frv2[] = {xfe_QuadRefUniform, 0, 0, 0, 0, xfe_QuadRefUniform};
  int frv3[] = {0, xfe_QuadRefHoriz, 0, xfe_QuadRefHoriz, 0, 0};
  int frv4[] = {0, 0, xfe_QuadRefVert, 0, xfe_QuadRefVert, 0};
  int frv5[] = {xfe_QuadRefHoriz, 0, 0, 0, 0, xfe_QuadRefVert};
  int frv6[] = {0, xfe_QuadRefHoriz, 0, xfe_QuadRefVert, 0, 0};
  int frv7[] = {0, 0, 0, 0, 0, xfe_QuadRefHoriz};

  ierr = xf_Error(xf_FaceRefEffectOnElem(xfe_Quadrilateral, frv0, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefVert);

  ierr = xf_Error(xf_FaceRefEffectOnElem(xfe_Quadrilateral, frv1, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefUniform);

  ierr = xf_Error(xf_FaceRefEffectOnElem(xfe_Hexahedron, frv2, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceXY);

  ierr = xf_Error(xf_FaceRefEffectOnElem(xfe_Hexahedron, frv3, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceZ);

  ierr = xf_Error(xf_FaceRefEffectOnElem(xfe_Hexahedron, frv4, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceY);

  ierr = xf_Error(xf_FaceRefEffectOnElem(xfe_Hexahedron, frv5, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceX);

  ierr = xf_Error(xf_FaceRefEffectOnElem(xfe_Hexahedron, frv6, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceXZ);

  ierr = xf_Error(xf_FaceRefEffectOnElem(xfe_Hexahedron, frv7, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceY);
  
  return xf_OK;  
}


TEST_xf_Face0RefEffectOnFace()
{
  int ierr, ref;
  
  ierr = xf_Error(xf_Face0RefEffectOnFace(xfe_Segment, xfe_SegRefUniform, xfe_SegPosLeft, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_SegRefNone);
  
  ierr = xf_Error(xf_Face0RefEffectOnFace(xfe_Quadrilateral, xfe_QuadRefUniform, xfe_QuadPosSW, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefNone);
  
  ierr = xf_Error(xf_Face0RefEffectOnFace(xfe_Quadrilateral, xfe_QuadRefUniform, xfe_QuadPosLeft, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefHoriz);
  
  ierr = xf_Error(xf_Face0RefEffectOnFace(xfe_Quadrilateral, xfe_QuadRefUniform, xfe_QuadPosTop, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefVert);
  
  ierr = xf_Error(xf_Face0RefEffectOnFace(xfe_Quadrilateral, xfe_QuadRefHoriz, xfe_QuadPosTop, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefNone);
  
  ierr = xf_Error(xf_Face0RefEffectOnFace(xfe_Quadrilateral, xfe_QuadRefVert, xfe_QuadPosTop, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefVert);
  
  return xf_OK;  
}


TEST_xf_RotateFaceRef()
{
  int ierr, ref;
  
  ierr = xf_Error(xf_RotateFaceRef(xfe_Segment, xfe_SegRefUniform, 1, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_SegRefUniform);
  
  ierr = xf_Error(xf_RotateFaceRef(xfe_Quadrilateral, xfe_QuadRefUniform, 1, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefUniform);

  ierr = xf_Error(xf_RotateFaceRef(xfe_Quadrilateral, xfe_QuadRefHoriz, 1, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefVert);

  ierr = xf_Error(xf_RotateFaceRef(xfe_Quadrilateral, xfe_QuadRefVert, 6, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefHoriz);

  ierr = xf_Error(xf_RotateFaceRef(xfe_Quadrilateral, xfe_QuadRefHoriz, 7, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefHoriz);

  return xf_OK;  
}



TEST_xf_InvRotateFaceRef()
{
  int ierr, ref;
  
  ierr = xf_Error(xf_InvRotateFaceRef(xfe_Segment, xfe_SegRefUniform, 1, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_SegRefUniform);
  
  ierr = xf_Error(xf_InvRotateFaceRef(xfe_Quadrilateral, xfe_QuadRefUniform, 1, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefUniform);

  ierr = xf_Error(xf_InvRotateFaceRef(xfe_Quadrilateral, xfe_QuadRefHoriz, 1, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefVert);

  ierr = xf_Error(xf_InvRotateFaceRef(xfe_Quadrilateral, xfe_QuadRefVert, 6, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefHoriz);

  ierr = xf_Error(xf_InvRotateFaceRef(xfe_Quadrilateral, xfe_QuadRefHoriz, 7, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefHoriz);
  
  return xf_OK;  
}



TEST_xf_RefCompare()
{
  int ierr;
  enum xfe_RelType rel;
  
  ierr = xf_Error(xf_RefCompare(xfe_Quadrilateral, xfe_QuadRefHoriz, xfe_QuadRefVert, &rel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(rel, xfe_RelDisjoint);
  
  ierr = xf_Error(xf_RefCompare(xfe_Quadrilateral, xfe_QuadRefHoriz, xfe_QuadRefUniform, &rel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(rel, xfe_RelSubset);
  
  ierr = xf_Error(xf_RefCompare(xfe_Hexahedron, xfe_HexRefSliceX, xfe_HexRefUniform, &rel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(rel, xfe_RelSubset);

  ierr = xf_Error(xf_RefCompare(xfe_Hexahedron, xfe_HexRefSliceXY, xfe_HexRefSliceX, &rel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(rel, xfe_RelSuperset);

  ierr = xf_Error(xf_RefCompare(xfe_Hexahedron, xfe_HexRefSliceXZ, xfe_HexRefSliceZ, &rel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(rel, xfe_RelSuperset);

  ierr = xf_Error(xf_RefCompare(xfe_Hexahedron, xfe_HexRefSliceY, xfe_HexRefSliceXY, &rel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(rel, xfe_RelSubset);
  
  return xf_OK;  
}


TEST_xf_LeastCommonRef()
{
  int ierr, ref;
  
  ierr = xf_Error(xf_LeastCommonRef(xfe_Quadrilateral, xfe_QuadRefHoriz, xfe_QuadRefVert, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefUniform);
  
  ierr = xf_Error(xf_LeastCommonRef(xfe_Quadrilateral, xfe_QuadRefHoriz, xfe_QuadRefHoriz, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefHoriz);
  
  ierr = xf_Error(xf_LeastCommonRef(xfe_Hexahedron, xfe_HexRefSliceX, xfe_HexRefSliceXY, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceXY);
  
  ierr = xf_Error(xf_LeastCommonRef(xfe_Hexahedron, xfe_HexRefSliceXY, xfe_HexRefSliceYZ, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefUniform);
  
  ierr = xf_Error(xf_LeastCommonRef(xfe_Hexahedron, xfe_HexRefSliceZ, xfe_HexRefSliceX, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceXZ);
  
  ierr = xf_Error(xf_LeastCommonRef(xfe_Hexahedron, xfe_HexRefNone, xfe_HexRefSliceZ, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_HexRefSliceZ);

  return xf_OK;  
}


TEST_xf_Pos2FaceRef()
{
  int ierr, ref;
  
  ierr = xf_Error(xf_Pos2FaceRef(xfe_Quadrilateral, xfe_QuadPosSW, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefUniform);

  ierr = xf_Error(xf_Pos2FaceRef(xfe_Quadrilateral, xfe_QuadPosLeft, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefVert);

  ierr = xf_Error(xf_Pos2FaceRef(xfe_Quadrilateral, xfe_QuadPosTop, &ref));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(ref, xfe_QuadRefHoriz);
  
  return xf_OK;  
}

TEST_xf_Ref2nElem()
{
  int ierr, nelem;
  
  ierr = xf_Error(xf_Ref2nElem(xfe_Quadrilateral, xfe_QuadRefUniform, &nelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nelem, 4);
  
  ierr = xf_Error(xf_Ref2nElem(xfe_Quadrilateral, xfe_QuadRefHoriz, &nelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nelem, 2);
  
  ierr = xf_Error(xf_Ref2nElem(xfe_Hexahedron, xfe_HexRefUniform, &nelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nelem, 8);
  
  ierr = xf_Error(xf_Ref2nElem(xfe_Hexahedron, xfe_HexRefSliceX, &nelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nelem, 2);
  
  ierr = xf_Error(xf_Ref2nElem(xfe_Hexahedron, xfe_HexRefSliceZ, &nelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nelem, 2);
  
  ierr = xf_Error(xf_Ref2nElem(xfe_Hexahedron, xfe_HexRefSliceY, &nelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nelem, 2);
  
  ierr = xf_Error(xf_Ref2nElem(xfe_Hexahedron, xfe_HexRefSliceXZ, &nelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nelem, 4);
  
  ierr = xf_Error(xf_Ref2nElem(xfe_Hexahedron, xfe_HexRefSliceYZ, &nelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nelem, 4);

  return xf_OK;  
}


TEST_xf_Ref2nIFace()
{
  int ierr, niface;
  
  ierr = xf_Error(xf_Ref2nIFace(xfe_Quadrilateral, xfe_QuadRefUniform, &niface));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(niface, 4);
  
  ierr = xf_Error(xf_Ref2nIFace(xfe_Quadrilateral, xfe_QuadRefHoriz, &niface));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(niface, 1);
  
  ierr = xf_Error(xf_Ref2nIFace(xfe_Hexahedron, xfe_HexRefUniform, &niface));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(niface, 12);
  
  ierr = xf_Error(xf_Ref2nIFace(xfe_Hexahedron, xfe_HexRefSliceX, &niface));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(niface, 1);
  
  ierr = xf_Error(xf_Ref2nIFace(xfe_Hexahedron, xfe_HexRefSliceZ, &niface));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(niface, 1);
  
  ierr = xf_Error(xf_Ref2nIFace(xfe_Hexahedron, xfe_HexRefSliceY, &niface));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(niface, 1);
  
  ierr = xf_Error(xf_Ref2nIFace(xfe_Hexahedron, xfe_HexRefSliceXZ, &niface));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(niface, 4);
  
  ierr = xf_Error(xf_Ref2nIFace(xfe_Hexahedron, xfe_HexRefSliceYZ, &niface));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(niface, 4);

  return xf_OK;  
}


TEST_xf_ConvertPos()
{
  int ierr, pos2;
  
  ierr = xf_Error(xf_ConvertPos(xfe_Quadrilateral, xfe_QuadPosNone, 7, 3, &pos2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(pos2, xfe_QuadPosNone);

  ierr = xf_Error(xf_ConvertPos(xfe_Quadrilateral, xfe_QuadPosSW, 4, 1, &pos2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(pos2, xfe_QuadPosSE);

  ierr = xf_Error(xf_ConvertPos(xfe_Quadrilateral, xfe_QuadPosNW, 4, 1, &pos2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(pos2, xfe_QuadPosNE);

  ierr = xf_Error(xf_ConvertPos(xfe_Quadrilateral, xfe_QuadPosLeft, 4, 1, &pos2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(pos2, xfe_QuadPosRight);

  ierr = xf_Error(xf_ConvertPos(xfe_Quadrilateral, xfe_QuadPosNE, 4, 2, &pos2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(pos2, xfe_QuadPosSW);

  ierr = xf_Error(xf_ConvertPos(xfe_Quadrilateral, xfe_QuadPosNW, 1, 2, &pos2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(pos2, xfe_QuadPosSW);

  ierr = xf_Error(xf_ConvertPos(xfe_Quadrilateral, xfe_QuadPosLeft, 5, 0, &pos2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(pos2, xfe_QuadPosRight);

  return xf_OK;  
}


TEST_xf_PosCompare()
{
  int ierr;
  enum xfe_RelType rel;
  enum xfe_Bool compat;
  
  ierr = xf_Error(xf_PosCompare(xfe_Quadrilateral, xfe_QuadPosLeft, xfe_QuadPosRight, &rel, &compat));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(compat, xfe_True);
  xf_AssertEqual(rel, xfe_RelDisjoint);

  ierr = xf_Error(xf_PosCompare(xfe_Quadrilateral, xfe_QuadPosSW, xfe_QuadPosSW, &rel, &compat));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(compat, xfe_True);
  xf_AssertEqual(rel, xfe_RelEqual);

  ierr = xf_Error(xf_PosCompare(xfe_Quadrilateral, xfe_QuadPosSW, xfe_QuadPosNE, &rel, &compat));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(compat, xfe_True);
  xf_AssertEqual(rel, xfe_RelDisjoint);
  
  ierr = xf_Error(xf_PosCompare(xfe_Quadrilateral, xfe_QuadPosLeft, xfe_QuadPosTop, &rel, &compat));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(compat, xfe_False);
  xf_AssertEqual(rel, xfe_RelDisjoint);
  
  ierr = xf_Error(xf_PosCompare(xfe_Quadrilateral, xfe_QuadPosBottom, xfe_QuadPosRight, &rel, &compat));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(compat, xfe_False);
  xf_AssertEqual(rel, xfe_RelDisjoint);

  ierr = xf_Error(xf_PosCompare(xfe_Quadrilateral, xfe_QuadPosTop, xfe_QuadPosNW, &rel, &compat));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(compat, xfe_True);
  xf_AssertEqual(rel, xfe_RelSuperset);

  ierr = xf_Error(xf_PosCompare(xfe_Quadrilateral, xfe_QuadPosSW, xfe_QuadPosLeft, &rel, &compat));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(compat, xfe_True);
  xf_AssertEqual(rel, xfe_RelSubset);

  return xf_OK;  
}


TEST_xf_RelativePos()
{
  int ierr, relpos;
  
  ierr = xf_Error(xf_RelativePos(xfe_Quadrilateral, xfe_QuadPosLeft, xfe_QuadPosSW, &relpos));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(relpos, xfe_QuadPosBottom);

  ierr = xf_Error(xf_RelativePos(xfe_Quadrilateral, xfe_QuadPosRight, xfe_QuadPosNE, &relpos));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(relpos, xfe_QuadPosTop);

  ierr = xf_Error(xf_RelativePos(xfe_Quadrilateral, xfe_QuadPosBottom, xfe_QuadPosSE, &relpos));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(relpos, xfe_QuadPosRight);

  ierr = xf_Error(xf_RelativePos(xfe_Quadrilateral, xfe_QuadPosTop, xfe_QuadPosNW, &relpos));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(relpos, xfe_QuadPosLeft);

  return xf_OK;  
}



TEST_xf_RefineGenericElem_Quad()
{
  int ierr;
  int nNodeQ1, nvert, nelem, niface, nbface;
  int elem2vert[ha_MAX_ELEM*ha_MAX_NODEQ1];
  int elem2pos[ha_MAX_ELEM];
  int vert2corner[ha_MAX_VERT];
  int bface2pos[ha_MAX_BFACE];
  real coord[ha_MAX_VERT*3];
  xf_IFace ifacelist[ha_MAX_IFACE];
  xf_BFace bfacelist[ha_MAX_BFACE];

  ierr = xf_Error(xf_RefineGenericElem_Quad(xfe_QuadRefUniform, &nNodeQ1, coord, &nvert,
					    vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					    ifacelist, &nbface, bfacelist, bface2pos));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 9);   xf_AssertEqual(nelem, 4);
  xf_AssertEqual(niface, 4);   xf_AssertEqual(nbface, 8);
  xf_AssertEqual(vert2corner[2], 1);
  xf_AssertEqual(elem2vert[2*4+1], 4);
  xf_AssertEqual(elem2pos[2], xfe_QuadPosNW);
  xf_AssertEqual(ifacelist[3].ElemL, 2);
  xf_AssertEqual(bfacelist[3].Elem, 3);
  xf_AssertEqual(bface2pos[0], xfe_SegPosLeft);

  ierr = xf_Error(xf_RefineGenericElem_Quad(xfe_QuadRefHoriz, &nNodeQ1, coord, &nvert,
					    vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					    ifacelist, &nbface, bfacelist, bface2pos));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 6);   xf_AssertEqual(nelem, 2);
  xf_AssertEqual(niface, 1);   xf_AssertEqual(nbface, 6);
  xf_AssertEqual(vert2corner[2], -1);
  xf_AssertEqual(elem2vert[1*4+2], 4);
  xf_AssertEqual(elem2pos[1], xfe_QuadPosTop);
  xf_AssertEqual(ifacelist[0].ElemR, 1);
  xf_AssertEqual(bfacelist[2].Elem, 1);
  xf_AssertEqual(bface2pos[5], xfe_SegPosRight);

  return xf_OK;  
}


TEST_xf_RefineGenericElem_Hex()
{
  int ierr;
  int nNodeQ1, nvert, nelem, niface, nbface, nedge;
  int elem2vert[ha_MAX_ELEM*ha_MAX_NODEQ1];
  int elem2pos[ha_MAX_ELEM];
  int vert2corner[ha_MAX_VERT];
  int bface2pos[ha_MAX_BFACE];
  int elist[ha_MAX_EDGE*3];
  real coord[ha_MAX_VERT*3];
  xf_IFace ifacelist[ha_MAX_IFACE];
  xf_BFace bfacelist[ha_MAX_BFACE];

  ierr = xf_Error(xf_RefineGenericElem_Hex(xfe_HexRefUniform, &nNodeQ1, coord, &nvert,
					   vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					   ifacelist, &nbface, bfacelist, bface2pos,
					   &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 27);   xf_AssertEqual(nelem, 8);
  xf_AssertEqual(niface, 12);   xf_AssertEqual(nbface, 24);
  xf_AssertEqual(vert2corner[18], 4);
  xf_AssertEqual(elem2vert[5*8+2], 13);
  xf_AssertEqual(elem2pos[6], xfe_HexPos011);
  xf_AssertEqual(ifacelist[9].ElemR, 7);
  xf_AssertEqual(bfacelist[11].Elem, 7);
  xf_AssertEqual(bface2pos[10], xfe_QuadPosNW);
  xf_AssertEqual(nedge, 12);

  ierr = xf_Error(xf_RefineGenericElem_Hex(xfe_HexRefSliceX, &nNodeQ1, coord, &nvert,
					   vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					   ifacelist, &nbface, bfacelist, bface2pos,
					   &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 12);   xf_AssertEqual(nelem, 2);
  xf_AssertEqual(niface, 1);   xf_AssertEqual(nbface, 10);
  xf_AssertEqual(vert2corner[9], 6);
  xf_AssertEqual(elem2vert[1*8+2], 4);
  xf_AssertEqual(elem2pos[1], xfe_HexPos122);
  xf_AssertEqual(ifacelist[0].ElemR, 1);
  xf_AssertEqual(bfacelist[1].Elem, 1);
  xf_AssertEqual(bface2pos[9], xfe_QuadPosRight);
  xf_AssertEqual(nedge, 4);

  ierr = xf_Error(xf_RefineGenericElem_Hex(xfe_HexRefSliceY, &nNodeQ1, coord, &nvert,
					   vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					   ifacelist, &nbface, bfacelist, bface2pos,
					   &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 12);   xf_AssertEqual(nelem, 2);
  xf_AssertEqual(niface, 1);   xf_AssertEqual(nbface, 10);
  xf_AssertEqual(vert2corner[7], 5);
  xf_AssertEqual(elem2vert[0*8+2], 2);
  xf_AssertEqual(elem2pos[0], xfe_HexPos202);
  xf_AssertEqual(ifacelist[0].ElemL, 0);
  xf_AssertEqual(bfacelist[1].Elem, 1);
  xf_AssertEqual(bface2pos[0], xfe_QuadPosLeft);
  xf_AssertEqual(nedge, 4);

  ierr = xf_Error(xf_RefineGenericElem_Hex(xfe_HexRefSliceZ, &nNodeQ1, coord, &nvert,
					   vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					   ifacelist, &nbface, bfacelist, bface2pos,
					   &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 12);   xf_AssertEqual(nelem, 2);
  xf_AssertEqual(niface, 1);   xf_AssertEqual(nbface, 10);
  xf_AssertEqual(vert2corner[10], 6);
  xf_AssertEqual(elem2vert[1*8+2], 6);
  xf_AssertEqual(elem2pos[1], xfe_HexPos221);
  xf_AssertEqual(ifacelist[0].ElemL, 0);
  xf_AssertEqual(bfacelist[1].Elem, 0);
  xf_AssertEqual(bface2pos[3], xfe_QuadPosBottom);
  xf_AssertEqual(nedge, 4);

  ierr = xf_Error(xf_RefineGenericElem_Hex(xfe_HexRefSliceXY, &nNodeQ1, coord, &nvert,
					   vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					   ifacelist, &nbface, bfacelist, bface2pos,
					   &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 18);   xf_AssertEqual(nelem, 4);
  xf_AssertEqual(niface, 4);   xf_AssertEqual(nbface, 16);
  xf_AssertEqual(vert2corner[6], 2);
  xf_AssertEqual(elem2vert[2*8+5], 13);
  xf_AssertEqual(elem2pos[3], xfe_HexPos112);
  xf_AssertEqual(ifacelist[1].ElemR, 3);
  xf_AssertEqual(bfacelist[3].Elem, 3);
  xf_AssertEqual(bface2pos[6], xfe_QuadPosLeft);
  xf_AssertEqual(nedge, 8);

  ierr = xf_Error(xf_RefineGenericElem_Hex(xfe_HexRefSliceXZ, &nNodeQ1, coord, &nvert,
					   vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					   ifacelist, &nbface, bfacelist, bface2pos,
					   &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 18);   xf_AssertEqual(nelem, 4);
  xf_AssertEqual(niface, 4);   xf_AssertEqual(nbface, 16);
  xf_AssertEqual(vert2corner[14], 5);
  xf_AssertEqual(elem2vert[2*8+1], 7);
  xf_AssertEqual(elem2pos[2], xfe_HexPos021);
  xf_AssertEqual(ifacelist[3].ElemR, 0);
  xf_AssertEqual(bfacelist[1].Elem, 1);
  xf_AssertEqual(bface2pos[0], xfe_QuadPosBottom);
  xf_AssertEqual(nedge, 8);

  ierr = xf_Error(xf_RefineGenericElem_Hex(xfe_HexRefSliceYZ, &nNodeQ1, coord, &nvert,
					   vert2corner, &nelem, elem2vert, elem2pos, &niface, 
					   ifacelist, &nbface, bfacelist, bface2pos,
					   &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(nvert, 18);   xf_AssertEqual(nelem, 4);
  xf_AssertEqual(niface, 4);   xf_AssertEqual(nbface, 16);
  xf_AssertEqual(vert2corner[12], 4);
  xf_AssertEqual(elem2vert[2*8+6], 14);
  xf_AssertEqual(elem2pos[1], xfe_HexPos210);
  xf_AssertEqual(ifacelist[1].ElemL, 1);
  xf_AssertEqual(bfacelist[11].Elem, 0);
  xf_AssertEqual(bface2pos[12], xfe_QuadPosNW);
  xf_AssertEqual(nedge, 8);



  return xf_OK;  
}


// Helper function for below
static int 
xf_FindRefIndicatorTest(xf_All *All, xf_Vector **pRefIndicator){
  int ierr;
  
  // Create a RefIndicator
  ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_False,  xfe_True, NULL, 
				pRefIndicator, NULL));
  xf_AssertEqual(ierr, xf_OK);

  // zero out vector
  ierr = xf_Error(xf_SetZeroVector(*pRefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK; 
}


TEST_xf_AdaptHang_Quad()
{
  int ierr;
  xf_Vector *RefIndicator;
  xf_All *All;

  // pull off a quad mesh
  ierr = xf_Error(xf_UnitQ1QuadAll(&All));
  xf_AssertEqual(ierr, xf_OK);
  
  // Find RefIndicator
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // Set RefIndicator to horiz
  RefIndicator->GenArray[0].iValue[0][0] = xfe_QuadRefHoriz;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // Find RefIndicator
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  
  // Set RefIndicator to uniform for elem 0
  RefIndicator->GenArray[0].iValue[0][0] = xfe_QuadRefUniform;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // verify mesh is correct
  xf_AssertEqual(All->Mesh->ElemGroup[0].nElem, 5);

  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;
}


TEST_xf_AdaptHang_HexQ1()
{
  int ierr;
  xf_Vector *RefIndicator;
  xf_All *All;

  // pull off a single-element hex mesh
  ierr = xf_Error(xf_UnitQ1HexAll(&All));
  xf_AssertEqual(ierr, xf_OK);
  

  // Set RefIndicator to SliceX
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefSliceX;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);


  // Set RefIndicator to Uniform refinement for elem 0
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefUniform;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);


  // Set RefIndicator to Uniform refinement for elem 0
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefUniform;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  
  // Set RefIndicator to SliceZ refinement for elem 3
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[3][0] = xfe_HexRefSliceZ;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  
  // Set RefIndicator to SliceY refinement for elem 2
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[2][0] = xfe_HexRefSliceY;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // verify mesh is correct
  xf_AssertEqual(All->Mesh->ElemGroup[0].nElem, 21);

  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;
}



TEST_xf_AdaptHang_HexQ1_Harder()
{
  int ierr;
  xf_Vector *RefIndicator;
  xf_All *All;

  // pull off a single-element hex mesh
  ierr = xf_Error(xf_UnitQ1HexAll(&All));
  xf_AssertEqual(ierr, xf_OK);
  
  // Set RefIndicator to SliceX
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefSliceX;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);


  // Set RefIndicator to SliceY for elem 0
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefSliceY;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);


  // Set RefIndicator to SliceZ for elem 1
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[1][0] = xfe_HexRefSliceZ;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);


  // Set RefIndicator to SliceYZ for elem 0
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefSliceYZ;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);


  // Set RefIndicator to SliceXZ for elem 1
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[1][0] = xfe_HexRefSliceXZ;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);


  // Set RefIndicator to SliceXY for elem 14
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[14][0] = xfe_HexRefSliceXY;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // verify mesh is correct
  xf_AssertEqual(All->Mesh->ElemGroup[0].nElem, 27);


  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;
}



TEST_xf_AdaptHang_HexQ1_CommonNode()
{
  int ierr;
  xf_Vector *RefIndicator;
  xf_All *All;

  // pull off a single-element hex mesh
  ierr = xf_Error(xf_UnitQ1HexAll(&All));
  xf_AssertEqual(ierr, xf_OK);
  
  // Set RefIndicator to SliceXY
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefSliceXY;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  

  // Set RefIndicator to SliceZ for elems 0 and 2
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefSliceZ;
  RefIndicator->GenArray[0].iValue[3][0] = xfe_HexRefSliceZ;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // verify mesh is correct
  xf_AssertEqual(All->Mesh->ElemGroup[0].nElem, 6);
  xf_AssertEqual(All->Mesh->nNode, 25); // this is key

/*   // write out file */
/*   ierr = xf_Error(xf_WriteAllBinary(All, "test.xfa")); */
/*   xf_AssertEqual(ierr, xf_OK); */

  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;
}


TEST_xf_AdaptHang_HexQ2()
{
  int ierr;
  xf_Vector *RefIndicator;
  xf_All *All;

  // pull off a single-element hex mesh
  ierr = xf_Error(xf_UnitQ2HexAll(&All));
  xf_AssertEqual(ierr, xf_OK);

  // Set RefIndicator to SliceX
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefSliceX;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // Set RefIndicator to SliceY for elem 0
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[0][0] = xfe_HexRefSliceY;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // Set RefIndicator to SliceZ for elem 1
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[1][0] = xfe_HexRefSliceZ;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);


  // Set RefIndicator to SliceYZ for elem 4
  ierr = xf_Error(xf_FindRefIndicatorTest(All, &RefIndicator));
  xf_AssertEqual(ierr, xf_OK);
  RefIndicator->GenArray[0].iValue[4][0] = xfe_HexRefSliceYZ;

  /* Hanging-node adaptation */
  ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
  xf_AssertEqual(ierr, xf_OK);

  // verify mesh is correct
  xf_AssertEqual(All->Mesh->ElemGroup[0].nElem, 10);

  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;
}



/* TEST_xf_AdaptHang_Custom() */
/* { */
/*   int ierr; */
/*   enum xfe_Bool Found; */
/*   xf_All *All; */
/*   xf_Vector *RefIndicator; */
  
/*   /\* Create .xfa structure *\/ */
/*   ierr = xf_Error(xf_CreateAll(&All, xfe_False)); */
/*   if (ierr != xf_OK) return ierr; */

/*   /\* Read .xfa file*\/ */
/*   ierr = xf_Error(xf_ReadAllBinary("Test.xfa", All)); */
/*   if (ierr!=xf_OK) return ierr; */
  
/*   // find RefIndicator */
/*   ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem, 1, NULL, 0, 0,  */
/* 				NULL, NULL, NULL, xfe_SizeInt, xfe_False,  xfe_True, NULL,  */
/* 				&RefIndicator, &Found)); */
/*   xf_AssertEqual(ierr, xf_OK); */
/*   xf_AssertEqual(Found, xfe_True); */

/*   // call adaptation */
/*   ierr = xf_Error(xf_AdaptHang(All, RefIndicator)); */
/*   xf_AssertEqual(ierr, xf_OK); */
  
/*   // write out file */
/*   ierr = xf_Error(xf_WriteAllBinary(All, "test.xfa")); */
/*   xf_AssertEqual(ierr, xf_OK); */

/*   // clean up */
/*   ierr = xf_Error(xf_DestroyAll(All)); */
/*   xf_AssertEqual(ierr, xf_OK); */

/*   return xf_OK;  */
/* } */
