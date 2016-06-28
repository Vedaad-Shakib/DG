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



TEST_xf_EgrpElem2Index()
{
  int ierr;
  int ie, elem, egrp;
  xf_Mesh *Mesh;

  // Q1/Q2 mesh (hex)
  ierr = xf_Error(xf_UnitQ1Q2HexMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  // test (egrp,elem) -> index
  ierr = xf_Error(xf_EgrpElem2Index(Mesh, 0, 0, &ie));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ie, 0);
  ierr = xf_Error(xf_EgrpElem2Index(Mesh, 0, 1, &ie));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ie, 1);

  // test index -> (egrp, elem)
  ierr = xf_Error(xf_Index2EgrpElem(Mesh, 0, &egrp, &elem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(egrp, 0);
  xf_AssertEqual(elem, 0);
  ierr = xf_Error(xf_Index2EgrpElem(Mesh, 1, &egrp, &elem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(egrp, 1);
  xf_AssertEqual(elem, 0);

  // destroy mesh
  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  // Q1 mesh of 8 elements
  ierr = xf_Error(xf_UnitBoxQ1TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  // test (egrp,elem) -> index
  ierr = xf_Error(xf_EgrpElem2Index(Mesh, 0, 5, &ie));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(ie, 5);

  // test index -> (egrp, elem)
  ierr = xf_Error(xf_Index2EgrpElem(Mesh, 7, &egrp, &elem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(egrp, 0);
  xf_AssertEqual(elem, 7);

  // destroy mesh
  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_GetRelativeOrient()
{
  int ierr, orientrel;
  
  ierr = xf_Error(xf_GetRelativeOrient(xfe_Segment, 0, 1, &orientrel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(orientrel, 1);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 0, 0, &orientrel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(orientrel, 0);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 2, 0, &orientrel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(orientrel, 2);
  
  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 2, 4, &orientrel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(orientrel, 6);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 5, 6, &orientrel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(orientrel, 1);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 6, 5, &orientrel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(orientrel, 3);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 6, 1, &orientrel));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(orientrel, 7);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 6, 2, &orientrel));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(orientrel, 4);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 1, 2, &orientrel));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(orientrel, 1);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 2, 1, &orientrel));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(orientrel, 3);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 0, 5, &orientrel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(orientrel, 5);

  ierr = xf_Error(xf_GetRelativeOrient(xfe_Quadrilateral, 5, 0, &orientrel));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(orientrel, 5);

  return xf_OK;  
}




TEST_xf_ReOrderHangNodes()
{
  int ierr, k;
  int nedge;
  int elist[12];
  
  // segment on a quad/tri
  int fvec0[] = {15, 21};
  int pvf0[] = {xfe_SegPosRight, xfe_SegPosLeft};
  int nvf0[] = {21, 4, 15, 4};
  int nvf0_true[] = {4, 21, 15, 4};
  
  // quad on a hex : uniform
  int fvec1[] = {1,4,9,16};
  int pvf1[] = {xfe_QuadPosNW, xfe_QuadPosSE, xfe_QuadPosSW, xfe_QuadPosNE};
  int nvf1[] = {30,25,16,17, 13,30,5,4, 5,1,17,30, 25,30,13,9};
  int nvf1_true[] = {17,30,25,16, 5,4,13,30, 1,5,30,17, 30,13,9,25};
  int elist1[] = {1,4,5, 4,9,13, 9,16,25, 16,1,17};

  // quad on a hex : left/right
  int fvec2[] = {1,4,9,16};
  int pvf2[] = {xfe_QuadPosRight, xfe_QuadPosLeft};
  int nvf2[] = {4,9,25,5, 16,25,5,1};
  int nvf2_true[] = {5,4,9,25, 1,5,25,16};
  int elist2[] = {1,4,5, 9,16,25};

  // quad on a hex : top/bottom
  int fvec3[] = {3,17,9,1};
  int pvf3[] = {xfe_QuadPosTop, xfe_QuadPosBottom};
  int nvf3[] = {9,5,8,1, 8,3,17,5};
  int nvf3_true[] = {8,5,9,1, 3,17,5,8};
  int elist3[] = {17,9,5, 1,3,8};

  ierr = xf_Error(xf_ReOrderHangNodes(xfe_Segment, 2, fvec0, 2, pvf0, nvf0, NULL, NULL));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertIntVectorEqual(nvf0, nvf0_true, 4);

  ierr = xf_Error(xf_ReOrderHangNodes(xfe_Quadrilateral, 4, fvec1, 4, pvf1, nvf1, &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(nvf1, nvf1_true, 16);
  xf_AssertEqual(nedge, 4);
  xf_AssertIntVectorEqual(elist, elist1, 12);

  ierr = xf_Error(xf_ReOrderHangNodes(xfe_Quadrilateral, 4, fvec2, 2, pvf2, nvf2, &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(nvf2, nvf2_true, 8);
  xf_AssertEqual(nedge, 2);
  xf_AssertIntVectorEqual(elist, elist2, 6);

  ierr = xf_Error(xf_ReOrderHangNodes(xfe_Quadrilateral, 4, fvec3, 2, pvf3, nvf3, &nedge, elist));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(nvf3, nvf3_true, 8);
  xf_AssertEqual(nedge, 2);
  xf_AssertIntVectorEqual(elist, elist3, 6);


  return xf_OK;  
}


TEST_xf_ScaleHangInterpolFace()
{
  int ierr;
  // quad, face 1
  real x0[] = {1,.5, 1,.7};
  real x0_true[] = {1,.25, 1,.35};
  // hex, face 0, pos SW
  real x1[] = {.6,.8,0, .4,.2,0};
  real x1_true[] = {.3,.4,0, .2,.1,0};
  // hex, face 0, pos SE
  real x2[] = {.6,.8,0, .4,.7,0};
  real x2_true[] = {.3,.9,0, .2,.85,0};
  // hex, face 1, pos NW
  real x3[] = {.5,0,.9, .3,0,.2};
  real x3_true[] = {.25,0,.95, .15,0,.6};
  // hex, face 2, pos NE
  real x4[] = {1,.1,.8, 1,.2,.3};
  real x4_true[] = {1,.55,.9, 1,.6,.65};
  // hex, face 3, pos SW
  real x5[] = {.4,1,.2, .6,1,.4};
  real x5_true[] = {.7,1,.1, .8,1,.2};
  // hex, face 4, pos Left
  real x6[] = {0,.8,.2, 0,.7,.3};
  real x6_true[] = {0,.9,.2, 0,.85,.3};
  // hex, face 5, pos Top
  real x7[] = {.2,.9,1, .3,.4,1};
  real x7_true[] = {.2,.95,1, .3,.7,1};
  // hex, face 4, pos Right
  real x8[] = {0,.6,.2, 0,.4,.3};
  real x8_true[] = {0,.3,.2, 0,.2,.3};
  // hex, face 4, pos Top
  real x9[] = {0,.6,.2, 0,.4,.3};
  real x9_true[] = {0,.6,.6, 0,.4,.65};

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Quadrilateral, 1, xfe_SegPosLeft, 2, x0));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x0, x0_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 0, xfe_QuadPosSW, 2, x1));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x1, x1_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 0, xfe_QuadPosSE, 2, x2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x2, x2_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 1, xfe_QuadPosNW, 2, x3));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x3, x3_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 2, xfe_QuadPosNE, 2, x4));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x4, x4_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 3, xfe_QuadPosSW, 2, x5));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x5, x5_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 4, xfe_QuadPosLeft, 2, x6));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x6, x6_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 5, xfe_QuadPosTop, 2, x7));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x7, x7_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 4, xfe_QuadPosRight, 2, x8));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x8, x8_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolFace(xfe_Hexahedron, 4, xfe_QuadPosTop, 2, x9));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x9, x9_true, 4, UTOL0);

  return xf_OK;  
}


TEST_xf_ScaleHangInterpol()
{
  int ierr;
  // quad, SE
  real x0[] = {1,.5, .8,.7};
  real x0_true[] = {1,.25, .9,.35};
  // quad, Left
  real x1[] = {.6,.8, .4,.2};
  real x1_true[] = {.3,.8, .2,.2};
  // quad, Top
  real x2[] = {.3,.7, .4,.6};
  real x2_true[] = {.3,.85, .4,.8,0};
  // hex, Uniform 101
  real x3[] = {.5,.2,.9, .3,0,.2};
  real x3_true[] = {.75,.1,.95, .65,0,.6};
  // hex, SliceX right, 122
  real x4[] = {1,.1,.8, .8,.2,.3};
  real x4_true[] = {1,.1,.8, .9,.2,.3};
  // hex, SliceY front, 202
  real x5[] = {.4,1,.2, .6,.5,.4};
  real x5_true[] = {.4,.5,.2, .6,.25,.4};
  // hex, SliceZ bottom, 220
  real x6[] = {0,.8,.2, .1,.7,.3};
  real x6_true[] = {0,.8,.1, .1,.7,.15};
  // hex, SliceXY 102
  real x7[] = {.2,.9,1, .3,.4,.5};
  real x7_true[] = {.6,.45,1, .65,.2,.5};
  // hex, SliceXZ 121
  real x8[] = {0,.6,.2, .8,.4,.3};
  real x8_true[] = {.5,.6,.6, .9,.4,.65};
  // hex, SliceYZ 201
  real x9[] = {0,.6,.4, .5,.4,.7};
  real x9_true[] = {0,.3,.7, .5,.2,.85};

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Quadrilateral, xfe_QuadPosSE, 2, x0));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x0, x0_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Quadrilateral, xfe_QuadPosLeft, 2, x1));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x1, x1_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Quadrilateral, xfe_QuadPosTop, 2, x2));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x2, x2_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Hexahedron, xfe_HexPos101, 2, x3));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x3, x3_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Hexahedron, xfe_HexPos122, 2, x4));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x4, x4_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Hexahedron, xfe_HexPos202, 2, x5));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x5, x5_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Hexahedron, xfe_HexPos220, 2, x6));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x6, x6_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Hexahedron, xfe_HexPos102, 2, x7));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x7, x7_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Hexahedron, xfe_HexPos121, 2, x8));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x8, x8_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpol(xfe_Hexahedron, xfe_HexPos201, 2, x9));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x9, x9_true, 6, UTOL0);


  return xf_OK;  
}

TEST_xf_ScaleHangInterpolInv()
{
  int ierr;
  // quad, SE
  real x0[] = {1,.5, .8,.7};
  real x0_true[] = {1,.25, .9,.35};
  // quad, Left
  real x1[] = {.6,.8, .4,.2};
  real x1_true[] = {.3,.8, .2,.2};
  // quad, Top
  real x2[] = {.3,.7, .4,.6};
  real x2_true[] = {.3,.85, .4,.8,0};
  // hex, Uniform 101
  real x3[] = {.5,.2,.9, .3,0,.2};
  real x3_true[] = {.75,.1,.95, .65,0,.6};
  // hex, SliceX right, 122
  real x4[] = {1,.1,.8, .8,.2,.3};
  real x4_true[] = {1,.1,.8, .9,.2,.3};
  // hex, SliceY front, 202
  real x5[] = {.4,1,.2, .6,.5,.4};
  real x5_true[] = {.4,.5,.2, .6,.25,.4};
  // hex, SliceZ bottom, 220
  real x6[] = {0,.8,.2, .1,.7,.3};
  real x6_true[] = {0,.8,.1, .1,.7,.15};
  // hex, SliceXY 102
  real x7[] = {.2,.9,1, .3,.4,.5};
  real x7_true[] = {.6,.45,1, .65,.2,.5};
  // hex, SliceXZ 121
  real x8[] = {0,.6,.2, .8,.4,.3};
  real x8_true[] = {.5,.6,.6, .9,.4,.65};
  // hex, SliceYZ 201
  real x9[] = {0,.6,.4, .5,.4,.7};
  real x9_true[] = {0,.3,.7, .5,.2,.85};

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Quadrilateral, xfe_QuadPosSE, 2, x0_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x0, x0_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Quadrilateral, xfe_QuadPosLeft, 2, x1_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x1, x1_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Quadrilateral, xfe_QuadPosTop, 2, x2_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x2, x2_true, 4, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Hexahedron, xfe_HexPos101, 2, x3_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x3, x3_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Hexahedron, xfe_HexPos122, 2, x4_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x4, x4_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Hexahedron, xfe_HexPos202, 2, x5_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x5, x5_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Hexahedron, xfe_HexPos220, 2, x6_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x6, x6_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Hexahedron, xfe_HexPos102, 2, x7_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x7, x7_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Hexahedron, xfe_HexPos121, 2, x8_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x8, x8_true, 6, UTOL0);

  ierr = xf_Error(xf_ScaleHangInterpolInv(xfe_Hexahedron, xfe_HexPos201, 2, x9_true));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(x9, x9_true, 6, UTOL0);


  return xf_OK;  
}






TEST_xf_RefFace2Elem_Triangle()
{
  int ierr;
  real xface[2] = {0.1, 0.9};
  real xelem[4];
  real true0[4] = {0.9, 0.1, 0.1, 0.9};
  real true1[4] = {0.1, 0.9, 0.9, 0.1};
  real true2[4] = {0.1, 0.0, 0.9, 0.0};

  ierr = xf_Error(xf_RefFace2Elem_Triangle(0, 0, 2, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xelem, true0, 4, UTOL1);
  
  ierr = xf_Error(xf_RefFace2Elem_Triangle(0, 1, 2, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(xelem, true1, 4, UTOL1);

  ierr = xf_Error(xf_RefFace2Elem_Triangle(2, 0, 2, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(xelem, true2, 4, UTOL1);

  return xf_OK;  
}

TEST_xf_RefFace2Elem_Quadrilateral()
{
  int ierr;
  real xface[2] = {0.1, 0.7};
  real xelem[4];
  real true0[4] = {0.1, 0.0, 0.7, 0.0};
  real true1[4] = {1.0, 0.1, 1.0, 0.7};
  real true2[4] = {0.9, 1.0, 0.3, 1.0};
  real true3[4] = {0.0, 0.9, 0.0, 0.3};
  real true4[4] = {1.0, 0.9, 1.0, 0.3};

  ierr = xf_Error(xf_RefFace2Elem_Quadrilateral(0, 0, 2, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xelem, true0, 4, UTOL1);
  
  ierr = xf_Error(xf_RefFace2Elem_Quadrilateral(1, 0, 2, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(xelem, true1, 4, UTOL1);

  ierr = xf_Error(xf_RefFace2Elem_Quadrilateral(2, 0, 2, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(xelem, true2, 4, UTOL1);

  ierr = xf_Error(xf_RefFace2Elem_Quadrilateral(3, 0, 2, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(xelem, true3, 4, UTOL1);

  ierr = xf_Error(xf_RefFace2Elem_Quadrilateral(1, 1, 2, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(xelem, true4, 4, UTOL1);

  return xf_OK;  
}

TEST_xf_RefFace2Elem_Hexahedron()
{
  int ierr;
  real xface[2] = {0.1, 0.2};
  real xelem[3];
  real true0[3] = {0.2, 0.1, 0.0};
  real true1[3] = {0.1, 0.2, 0.0};
  real true2[3] = {0.1, 0.8, 0.0};
  real true3[3] = {1.0, 0.9, 0.8};

  ierr = xf_Error(xf_RefFace2Elem_Hexahedron(0, 0, 1, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xelem, true0, 3, UTOL1);
  
  ierr = xf_Error(xf_RefFace2Elem_Hexahedron(0, 4, 1, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xelem, true1, 3, UTOL1);

  ierr = xf_Error(xf_RefFace2Elem_Hexahedron(0, 1, 1, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xelem, true2, 3, UTOL1);

  ierr = xf_Error(xf_RefFace2Elem_Hexahedron(2, 2, 1, xface, xelem));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xelem, true3, 3, UTOL1);


  return xf_OK;  
}


TEST_xf_Ref2GlobElem()
{
  int ierr;
  real xref0[4] = {0.1,0.0, 0.0,0.5};
  real xref1[4] = {0.0,0.0, 0.5,0.5};
  real xglob[4];
  real true0[4] = {1.0,0.1, 0.5,0.5};
  real true1[4] = {1.0,0.0, 0.5,1.0};
  xf_Mesh *Mesh;
  xf_BasisData *PhiData;

  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_Ref2GlobElem(Mesh, 0, 1, NULL, xfe_True, 
				  2, xref0, xglob));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xglob, true0, 4, UTOL1);


  PhiData = NULL;
  ierr = xf_Error(xf_Ref2GlobElem(Mesh, 0, 1, &PhiData, xfe_True, 
				  2, xref0, xglob));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xglob, true0, 4, UTOL1);


  ierr = xf_Error(xf_Ref2GlobElem(Mesh, 0, 1, &PhiData, xfe_True, 
				  2, xref1, xglob));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xglob, true1, 4, UTOL1);

  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_Glob2RefElem()
{
  int ierr;
  enum xfe_Bool converged;
  real xglob0[2] = {0.5,0.5};
  real xglob1[2] = {0.5,1.0};
  real xglob2[2] = {1.0,0.05};
  real xglob3[3] = {.125,.01,.02};
  real xref[2];
  real true0[2] = {0.0,0.5};
  real true1[2] = {0.5,0.5};
  real true2[2] = {0.5,0.0};
  real true3[3] = {.5,0,0};
  xf_Mesh *Mesh;

  // Q1 mesh
  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_Glob2RefElem(Mesh, 0, 1, xglob0, -1.0, xfe_True, xref, &converged));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(converged, xfe_True);  
  xf_AssertRealVectorWithin(xref, true0, 2, UTOL1);

  ierr = xf_Error(xf_Glob2RefElem(Mesh, 0, 1, xglob1, -1.0, xfe_True, xref, NULL));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertRealVectorWithin(xref, true1, 2, UTOL1);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  // Q2 mesh (tri)
  ierr = xf_Error(xf_UnitBoxQ2TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_Glob2RefElem(Mesh, 0, 0, xglob2, -1.0, xfe_True, xref, &converged));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(converged, xfe_True);  
  xf_AssertRealVectorWithin(xref, true2, 2, UTOL3);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  // Q1/Q2 mesh (hex)
  ierr = xf_Error(xf_UnitQ1Q2HexMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_Glob2RefElem(Mesh, 0, 0, xglob3, -1.0, xfe_True, xref, &converged));
  xf_AssertEqual(ierr, xf_OK);  
  xf_AssertEqual(converged, xfe_True);  
  xf_AssertRealVectorWithin(xref, true3, 3, UTOL3);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}




TEST_xf_ElemNormal_Q1Triangle()
{
  int ierr;
  real s[1] = {0.5}, n_true[2];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 0, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 2);
  n_true[0] = 1.0; n_true[1] = 1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 1, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 2);
  n_true[0] = -1.0; n_true[1] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_ElemNormal_Q2Triangle()
{
  int ierr;
  real s[1] = {0.5}, n_true[2];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitTwoQ2TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 0, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 2);
  n_true[0] = 1.0; n_true[1] = 1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 1, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 2);
  n_true[0] = -1.0; n_true[1] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_ElemNormal_Q1Tet()
{
  int ierr;
  real s[2] = {0.0, 0.0}, n_true[3];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitQ1TetMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 0, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 3);
  n_true[0] = 1.0; n_true[1] = 1.0; n_true[2] = 1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 1, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = -1.0; n_true[1] = 0.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 2, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = -1.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 3, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = 0.0; n_true[2] = -1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_ElemNormal_Q2Tet()
{
  int ierr;
  real s[2] = {0.0, 0.0}, n_true[3];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitQ2TetMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 0, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 3);
  n_true[0] = 1.0; n_true[1] = 1.0; n_true[2] = 1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 1, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = -1.0; n_true[1] = 0.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 2, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = -1.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 3, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = 0.0; n_true[2] = -1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_ElemNormal_Q1Quad()
{
  int ierr;
  real s[1] = {0.5}, n_true[2];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitQ1QuadMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 0, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 2);
  n_true[0] = 0.0; n_true[1] = -1.0; 
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 1, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 1.0; n_true[1] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 2, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = 1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 3, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = -1.0; n_true[1] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}

TEST_xf_ElemNormal_Q2Quad()
{
  int ierr;
  real s[1] = {0.5}, n_true[2];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitQ2QuadMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 0, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 2);
  n_true[0] = 0.0; n_true[1] = -1.0; 
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 1, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 1.0; n_true[1] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 2, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = 1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 3, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = -1.0; n_true[1] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_ElemNormal_Q1Hex()
{
  int ierr;
  real s[2] = {0.0, 0.0}, n_true[3];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitQ1HexMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 0, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 3);
  n_true[0] = 0.0; n_true[1] = 0.0; n_true[2] = -1.0; 
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 1, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = -1.0; n_true[2] = 0.0; 
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 2, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 1.0; n_true[1] = 0.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 3, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = 1.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 4, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = -1.0; n_true[1] = 0.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 5, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = 0.0; n_true[2] = 1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_ElemNormal_Q2Hex()
{
  int ierr;
  real s[2] = {0.0, 0.0}, n_true[3];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitQ2HexMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 0, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 3);
  n_true[0] = 0.0; n_true[1] = 0.0; n_true[2] = -4.0; 
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 1, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = -4.0; n_true[2] = 0.0; 
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 2, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 4.0; n_true[1] = 0.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 3, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = 4.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 4, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = -4.0; n_true[1] = 0.0; n_true[2] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_ElemNormal(Mesh, 0, 0, 5, 0, 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  n_true[0] = 0.0; n_true[1] = 0.0; n_true[2] = 4.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 3, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_IFaceNormal()
{
  int ierr;
  real s[1] = {0.5}, n_true[2];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_IFaceNormal(Mesh, Mesh->IFace[0], 1, s, &NData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 2);
  n_true[0] = 1.0; n_true[1] = 1.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_BFaceNormal()
{
  int ierr;
  real s[1] = {0.5}, n_true[2];
  xf_NormalData *NData;
  xf_Mesh *Mesh;

  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  NData = NULL;

  ierr = xf_Error(xf_BFaceNormal(Mesh, Mesh->BFaceGroup[0].BFace[1], 1, s, &NData, NULL));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(NData->nq, 1);
  xf_AssertEqual(NData->dim, 2);
  n_true[0] = 1.0; n_true[1] = 0.0;
  xf_AssertRealVectorWithin(NData->n, n_true, 2, UTOL1);

  ierr = xf_Error(xf_DestroyNormalData(NData));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}

TEST_xf_GetRefineCoords()
{
  int ierr, nn, ns, nb;
  int *vs, *vb;
  int vs1[3] = {0,1,2};
  int vb1[3] = {1,2,0};
  int vs2[12] = {0,1,3, 1,4,3, 1,2,4, 3,4,5};
  int vb2[6] = {2,4,5,3,0,1};
  real *x;
  real x1[6] ={0,0, 1,0, 0,1};
  real x2[12] ={0,0, .5,0, 1,0, 0,.5, .5,.5, 0,1};

  x  = NULL;
  vs = NULL;
  vb = NULL;

  ierr = xf_Error(xf_GetRefineCoords(xfe_Triangle, 1, &nn, &x, &ns, &vs, &nb, &vb));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nn, 3);
  xf_AssertRealVectorWithin(x, x1, 6, UTOL1);
  xf_AssertEqual(ns, 1);
  xf_AssertIntVectorEqual(vs, vs1, 3);
  xf_AssertEqual(nb, 3);
  xf_AssertIntVectorEqual(vb, vb1, 3);

  ierr = xf_Error(xf_GetRefineCoords(xfe_Triangle, 2, &nn, &x, &ns, &vs, &nb, &vb));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nn, 6);
  xf_AssertRealVectorWithin(x, x2, 12, UTOL1);
  xf_AssertEqual(ns, 4);
  xf_AssertIntVectorEqual(vs, vs2, 12);
  xf_AssertEqual(nb, 6);
  xf_AssertIntVectorEqual(vb, vb2, 6);

  xf_Release( (void *) x);
  xf_Release( (void *) vs);
  xf_Release( (void *) vb);

  return xf_OK;  
}


TEST_xf_GetRefineCoordsTet()
{
  int ierr, nn, ns, nb;
  int *vs;
  int vs1[32] = {0,1,3,6, 1,7,3,6, 3,7,8,6, 1,4,3,7, 
		 3,7,4,8, 1,2,4,7, 3,4,5,8, 6,7,8,9};
  real *x;
  real x1[30] = {0,0,0, .5,0,0, 1,0,0, 0,.5,0, .5,.5,0, 0,1,0,
	 	 0,0,.5, .5,0,.5, 0,.5,.5, 0,0,1};
  x  = NULL;
  vs = NULL;

  ierr = xf_Error(xf_GetRefineCoords(xfe_Tetrahedron, 2, &nn, &x, &ns, &vs, NULL, NULL));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nn, 10);
  xf_AssertRealVectorWithin(x, x1, 30, UTOL1);
  xf_AssertEqual(ns, 8);
  xf_AssertIntVectorEqual(vs, vs1, 32);

  xf_Release( (void *) x);
  xf_Release( (void *) vs);

  return xf_OK;  
}

TEST_xf_GetRefineCoordsHex()
{
  int ierr, nn, ns, nb;
  int *vs;
  int vs1[24] = {0,1,2,4, 1,2,4,5, 6,5,4,2, 1,3,2,5, 2,5,3,6, 7,5,6,3};
  real *x;
  real x1[24] ={0,0,0, 1,0,0, 0,1,0, 1,1,0,
		0,0,1, 1,0,1, 0,1,1, 1,1,1};
  x  = NULL;
  vs = NULL;

  ierr = xf_Error(xf_GetRefineCoords(xfe_Hexahedron, 1, &nn, &x, &ns, &vs, NULL, NULL));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nn, 8);
  xf_AssertRealVectorWithin(x, x1, 24, UTOL1);
  xf_AssertEqual(ns, 6);
  xf_AssertIntVectorEqual(vs, vs1, 24);

  xf_Release( (void *) x);
  xf_Release( (void *) vs);

  return xf_OK;  
}



TEST_xf_FindElemGeom()
{
  int ierr;
  real Perim;
  xf_Vector *EG;
  xf_All *All;

  ierr = xf_Error(xf_UnitTwoQ1TriangleAll(&All));
  xf_AssertEqual(ierr, xf_OK);

  Perim = 2.0+sqrt(2.0);

  // obtain elem geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(EG->GenArray[0].rValue[0][xfe_EGVolume]  , 0.5,   UTOL1);
  xf_AssertWithin(EG->GenArray[0].rValue[0][xfe_EGSurfArea], Perim, UTOL1);
  xf_AssertWithin(EG->GenArray[0].rValue[1][xfe_EGVolume]  , 0.5,   UTOL1);
  xf_AssertWithin(EG->GenArray[0].rValue[1][xfe_EGSurfArea], Perim, UTOL1);
  
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  // Same test for Q2
  ierr = xf_Error(xf_UnitTwoQ2TriangleAll(&All));
  xf_AssertEqual(ierr, xf_OK);

  Perim = 2.0+sqrt(2.0);

  // obtain elem geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(EG->GenArray[0].rValue[0][xfe_EGVolume]  , 0.5,   UTOL1);
  xf_AssertWithin(EG->GenArray[0].rValue[0][xfe_EGSurfArea], Perim, UTOL1);
  xf_AssertWithin(EG->GenArray[0].rValue[1][xfe_EGVolume]  , 0.5,   UTOL1);
  xf_AssertWithin(EG->GenArray[0].rValue[1][xfe_EGSurfArea], Perim, UTOL1);
  
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_FindElemHMetric()
{
  int ierr, i;
  real H[4] = {1.0, 0.0, 0.0, 1.0};
  xf_Vector *EM;
  xf_All *All;

  ierr = xf_Error(xf_UnitQ2QuadAll(&All));
  xf_AssertEqual(ierr, xf_OK);


  // obtain elem geometry vector
  ierr = xf_Error(xf_FindElemHMetric(All, xfe_True, &EM));
  xf_AssertErrorOK(ierr);
  
  for (i=0; i<4; i++){
    xf_AssertRealVectorWithin(EM->GenArray[0].rValue[0]+i*4, H, 4, UTOL1);
  }
  
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_ElemSize()
{
  int ierr;
  real Perim, h;
  xf_Vector *EG;
  xf_All *All;

  ierr = xf_Error(xf_UnitTwoQ1TriangleAll(&All));
  xf_AssertEqual(ierr, xf_OK);

  Perim = 2.0+sqrt(2.0);

  // obtain elem geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_ElemSize(All, 0, 0, EG, &h));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(h, 2.*0.5/Perim, UTOL1);

  ierr = xf_Error(xf_ElemSize(All, 0, 1, EG, &h));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(h, 2.*0.5/Perim, UTOL1);
  
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_FindElemUsingSearchStructure()
{
  int ierr, egrp, elem;
  real xglob[2] = {1.3, 1.4};
  real xref[2];
  xf_ElemSearchStruct ESS;
  xf_Vector *EG;
  xf_All *All;

  ierr = xf_Error(xf_UnitBoxQ1Quad9All(&All));
  xf_AssertEqual(ierr, xf_OK);

  // build element search structure
  ierr = xf_Error(xf_BuildElemSearchStructure(All, &ESS));
  xf_AssertEqual(ierr, xf_OK);
  
  // find element that contains the point
  ierr = xf_FindElemUsingSearchStructure(All, xglob, &ESS, &egrp, &elem, xref);
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(egrp, 0);
  xf_AssertEqual(elem, 4);
  xf_AssertWithin(xref[0], 0.3, UTOL1);
  xf_AssertWithin(xref[1], 0.4, UTOL1);
  
  // destroy element search structure
  ierr = xf_Error(xf_DestroyElemSearchStructure(&ESS));
  if (ierr != xf_OK) return ierr;

  // destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}
