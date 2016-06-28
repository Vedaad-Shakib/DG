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

TEST_xf_Dist2Segment()
{
  int ierr;
  real xseg[4] = {1,1,5,3}, d, xp[2];
  real x0[2] = {1,6};
  real x1[2] = {1,1};
  real x2[2] = {5,3};
  real x3[2] = {2,1};
  real x4[2] = {0,0};

  ierr = xf_Error(xf_Dist2Segment(xseg, x0, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d,     2.*sqrt(5.), UTOL3);
  xf_AssertWithin(xp[0], 3., UTOL3);
  xf_AssertWithin(xp[1], 2., UTOL3);

  ierr = xf_Error(xf_Dist2Segment(xseg, x1, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 0., UTOL3);
  xf_AssertWithin(xp[0], 1., UTOL3);
  xf_AssertWithin(xp[1], 1., UTOL3);

  ierr = xf_Error(xf_Dist2Segment(xseg, x2, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 0., UTOL3);
  xf_AssertWithin(xp[0], 5., UTOL3);
  xf_AssertWithin(xp[1], 3., UTOL3);

  ierr = xf_Error(xf_Dist2Segment(xseg, x3, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 1./sqrt(5.), UTOL3);
  xf_AssertWithin(xp[0], 1.8, UTOL3);
  xf_AssertWithin(xp[1], 1.4, UTOL3);

  ierr = xf_Error(xf_Dist2Segment(xseg, x4, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, sqrt(2.0), UTOL3);
  xf_AssertWithin(xp[0], 1., UTOL3);
  xf_AssertWithin(xp[1], 1., UTOL3);

  return xf_OK;
}


TEST_xf_Dist2Segment3d()
{
  int ierr;
  real xseg[6] = {0,0,1, 0,0,4}, d, xp[3];
  real x0[3] = {3,4,1} , xp0[3] = {0,0,1};
  real x1[3] = {-8,6,4}, xp1[3] = {0,0,4};
  real x2[3] = {1,1,0},  xp2[3] = {0,0,1};
  real x3[3] = {0,0,1},  xp3[3] = {0,0,1};
  real x4[3] = {0,0,4},  xp4[3] = {0,0,4};
  real x5[3] = {5,12,2}, xp5[3] = {0,0,2};

  ierr = xf_Error(xf_Dist2Segment3d(xseg, x0, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 5., UTOL3);
  xf_AssertRealVectorWithin(xp, xp0, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Segment3d(xseg, x1, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 10., UTOL3);
  xf_AssertRealVectorWithin(xp, xp1, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Segment3d(xseg, x2, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, sqrt(3.), UTOL3);
  xf_AssertRealVectorWithin(xp, xp2, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Segment3d(xseg, x3, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 0., UTOL3);
  xf_AssertRealVectorWithin(xp, xp3, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Segment3d(xseg, x4, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 0., UTOL3);
  xf_AssertRealVectorWithin(xp, xp4, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Segment3d(xseg, x5, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 13., UTOL3);
  xf_AssertRealVectorWithin(xp, xp5, 3, UTOL3);

  return xf_OK;
}


TEST_xf_Dist2Triangle()
{
  int ierr;
  real xtri[9] = {0,0,0, 2,0,0, 0,2,0}, d, xp[3];
  real x0[3] = {0,0,0},      xp0[3] = {0,0,0};
  real x1[3] = {2,0,0},      xp1[3] = {2,0,0};
  real x2[3] = {0,2,0},      xp2[3] = {0,2,0};
  real x3[3] = {.3,.3,1.5},  xp3[3] = {.3,.3,0};
  real x4[3] = {.2,.3,-1.7}, xp4[3] = {.2,.3,0};
  real x5[3] = {2,2,1},      xp5[3] = {1,1,0};
  real x6[3] = {1,-2,0},     xp6[3] = {1,0,0};
  real x7[3] = {-1,1,1},     xp7[3] = {0,1,0};

  ierr = xf_Error(xf_Dist2Triangle(xtri, x0, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 0., UTOL3);
  xf_AssertRealVectorWithin(xp, xp0, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Triangle(xtri, x1, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 0., UTOL3);
  xf_AssertRealVectorWithin(xp, xp1, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Triangle(xtri, x2, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 0., UTOL3);
  xf_AssertRealVectorWithin(xp, xp2, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Triangle(xtri, x3, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 1.5, UTOL3);
  xf_AssertRealVectorWithin(xp, xp3, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Triangle(xtri, x4, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 1.7, UTOL3);
  xf_AssertRealVectorWithin(xp, xp4, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Triangle(xtri, x5, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, sqrt(3), UTOL3);
  xf_AssertRealVectorWithin(xp, xp5, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Triangle(xtri, x6, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, 2.0, UTOL3);
  xf_AssertRealVectorWithin(xp, xp6, 3, UTOL3);

  ierr = xf_Error(xf_Dist2Triangle(xtri, x7, &d, xp));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(d, sqrt(2.0), UTOL3);
  xf_AssertRealVectorWithin(xp, xp7, 3, UTOL3);

  return xf_OK;
}


TEST_xf_Dist2Face()
{
  int ierr;
  // Flat Tri
  real x1a[3] = {0.2, 0.2, 1.0},       xp1a[3] = {0.2, 0.2, 0};
  real x1b[3] = {0.3, 0.2, 1.0},       xp1b[3] = {0.3, 0.2, 0};
  real x1c[3] = {0.9, 0.01, 1.0},      xp1c[3] = {0.9, 0.01, 0.0};
  real x1d[3] = {-0.6, 1.3, 0.8},      xp1d[3] = {0.0, 1.3, 0.0}; 
  real x1e[3] = {0.4, -6.0, -8.0},     xp1e[3] = {0.4, 0.0, 0.0}; 
  real x1f[3] = {10, 10, 0.0},         xp1f[3] = {1,1,0};
  real X1[9] = {0,0,0, 2,0,0, 0,2,0};
  // Flat Quad
  real x2a[3] = {0.2, 0.2, 1.0},       xp2a[3] = {0.2, 0.2, 0};
  real x2b[3] = {0.3, 0.2, 1.0},       xp2b[3] = {0.3, 0.2, 0};
  real x2c[3] = {0.3, 2.6, -0.8},      xp2c[3] = {0.3, 2.0, 0};
  real X2[12] = {0,0,0, 2,0,0, 0,2,0, 2,2,0};
  // Curved Tri
  real x3a[3] = {0.1, 1.0, 1.0},       xp3a[3] = {0.1,1.0,0.2};
  real X3[18] = {0,0,0, 1,0,0, 2,0,0, 0,1,.2, 1,1,.2, 0,2,0};
  // Curved Quad
  real x4a[3] = {-0.3, 1.0, 0.5},      xp4a[3] = {0,1.0,.1};
  real X4[27] = {0,0,0, 1,0,.1, 2,0,0, 0,1,.1, 1,1,.3, 2,1,.1, 0,2,0, 1,2,.1, 2,2,0};
  // Flat Seg
  real x5a[2] = {0,5.0/3.0},           xp5a[3] = {0.8, 0.6};
  real X5[4] = {0,0, 4,3};
  // Curved Seg
  real x6a[2] = {0.5,-1.0},            xp6a[3] = {0.5,-1./8.};
  real X6[6] = {0,0, 1,0, 2,1};
  // Flat Seg, interp D
  real x7a[2] = {0.2, 0.7};
  real X7[4] = {0,0, 1,1};
  // Curved Tri, interp D
  real x8a[3] = {0.1, 1.3, .7};
  real X8[18] = {0,0,0, 1,0,0, 2,0,0, 0,1,.2, 1,1,.2, 0,2,.4};

  real dist, xp[3];
  
  // distance to flat tri face from various global points
  ierr = xf_Error(xf_Dist2Face(3, x1a, xfe_TriLagrange, 1, 3, X1, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 1.0, UTOL3);
  xf_AssertRealVectorWithin(xp, xp1a, 3, UTOL3);
  ierr = xf_Error(xf_Dist2Face(3, x1b, xfe_TriLagrange, 1, 3, X1, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 1.0, UTOL3);
  xf_AssertRealVectorWithin(xp, xp1b, 3, UTOL3);
  ierr = xf_Error(xf_Dist2Face(3, x1c, xfe_TriLagrange, 1, 3, X1, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 1.0, UTOL3);
  xf_AssertRealVectorWithin(xp, xp1c, 3, UTOL3);
  ierr = xf_Error(xf_Dist2Face(3, x1d, xfe_TriLagrange, 1, 3, X1, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 1.0, UTOL3);  // edge 1 is closest 
  xf_AssertRealVectorWithin(xp, xp1d, 3, UTOL3);
  ierr = xf_Error(xf_Dist2Face(3, x1e, xfe_TriLagrange, 1, 3, X1, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 10.0, UTOL3); // edge 2 is closest 
  xf_AssertRealVectorWithin(xp, xp1e, 3, UTOL3);
  ierr = xf_Error(xf_Dist2Face(3, x1f, xfe_TriLagrange, 1, 3, X1, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 9.0*sqrt(2.0), UTOL3); // edge 0 is closest 
  xf_AssertRealVectorWithin(xp, xp1f, 3, UTOL3);


  // distance to flat quad face
  ierr = xf_Error(xf_Dist2Face(3, x2a, xfe_QuadLagrange, 1, 4, X2, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 1.0, UTOL3);
  xf_AssertRealVectorWithin(xp, xp2a, 3, UTOL3);
  ierr = xf_Error(xf_Dist2Face(3, x2b, xfe_QuadLagrange, 1, 4, X2, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 1.0, UTOL3);
  xf_AssertRealVectorWithin(xp, xp2b, 3, UTOL3);
  ierr = xf_Error(xf_Dist2Face(3, x2c, xfe_QuadLagrange, 1, 4, X2, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 1.0, UTOL3);
  xf_AssertRealVectorWithin(xp, xp2c, 3, UTOL3);

  // distance to curved tri face
  ierr = xf_Error(xf_Dist2Face(3, x3a, xfe_TriLagrange, 2, 6, X3, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 0.8, UTOL3);
  xf_AssertRealVectorWithin(xp, xp3a, 3, UTOL3);

  // distance to curved quad face
  ierr = xf_Error(xf_Dist2Face(3, x4a, xfe_QuadLagrange, 2, 9, X4, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 0.5, UTOL3);
  xf_AssertRealVectorWithin(xp, xp4a, 3, UTOL3);

  // distance to flat seg face
  ierr = xf_Error(xf_Dist2Face(2, x5a, xfe_SegLagrange, 1, 2, X5, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 4.0/3.0, UTOL3);
  xf_AssertRealVectorWithin(xp, xp5a, 2, UTOL3);

  // distance to curved seg face
  ierr = xf_Error(xf_Dist2Face(2, x6a, xfe_SegLagrange, 2, 3, X6, &dist, xp));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(dist, 7.0/8.0, 0.1);
  xf_AssertRealVectorWithin(xp, xp6a, 2, 0.1);

  return xf_OK;  
}




/* A .gri file for a Q2 quad wedge */
static char *WedgeQ2QuadGri[] =
{
  "15 2 2",
  "-1 1", "0 1.414", "1 1",
  "-2 2", "0 2.828", "2 2",
  "-3 3", "0 4.242", "3 3",
  "-4 4", "0 5.656", "4 4",
  "-5 5", "0 7.070", "5 5",
  "4",
  "2 2 Left",
  "1 7","7 13",
  "2 2 Right",
  "3 9","9 15",
  "1 2 Bottom",
  "1 3",
  "1 2 Top",
  "13 15",
  "2 2 QuadLagrange",
  "1 2 3 4 5 6 7 8 9",
  "7 8 9 10 11 12 13 14 15",
  "\0",
};



TEST_xf_CalculateDistFcn_Wedge()
{
  int ierr;
  xf_All *All;

  // Create All
  ierr = xf_Error(xf_CreateAll(&All, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Read Wedge
  ierr = xf_Error( xf_ReadGriFile(NULL, WedgeQ2QuadGri, All->Mesh) );
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "WriteDistFcn", "True"));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "DistFcnOrder", "2"));
  if (ierr != xf_OK) return ierr;

  // Calculate distance function
  ierr = xf_Error(xf_CalculateDistFcn(All));
  if (ierr != xf_OK) return ierr;

  /* ierr = xf_Error(xf_WriteAllBinary(All, "dist.xfa")); */
  /*   if (ierr != xf_OK) return ierr;   */

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;
}



/*------------------*/
/*   LEGACY TESTS   */
/*------------------*/


/* TEST_xf_Dist2FaceResidual() */
/* { */
/*   int ierr, i, j; */
/*   int nf = 6; */
/*   real xi[2] = {0.5, 0.5}; */
/*   real x0[3] = {.5, .51, .2}; */
/*   real x1[3] = {0.3, 0.3, 1.0}; */
/*   real X[18] = {0,0,0, .5,0,.1, 1,0,.01, 0,.5,.11, .5,.51,.2, 0,1,-.02}; */
/*   real D[6] = {.5, .6, .57, .4, .55, .45}; */
/*   real dist; */
/*   real R[2], R0[2]; */
/*   real R_xi[4], R_xi0[4]; */
/*   real tol, eps, fd, an; */
/*   enum xfe_BasisType FBasis = xfe_TriLagrange; */
/*   xf_BasisData *PhiData = NULL; */

/*   // compute basis shape, grad, hess at xi */
/*   ierr = xf_Error(xf_EvalBasis(FBasis, 2, xfe_True, 1, xi, xfb_All, &PhiData)); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   // compute distance residual at x0 (which is on the face) */
/*   ierr = xf_Error(xf_Dist2FaceResidual(3, x0, nf, X, D, PhiData, R, R_xi, &dist)); */
/*   if (ierr != xf_OK) return ierr; */
/*   xf_AssertWithin(dist, .55, UTOL1); */
/*   xf_AssertWithin(R[0], 0.0, UTOL1); */
/*   xf_AssertWithin(R[1], 0.0, UTOL1); */

/*   // compute distance residual at x1 (off the face), get ready for pinging */
/*   ierr = xf_Error(xf_Dist2FaceResidual(3, x1, nf, X, D, PhiData, R0, R_xi0, &dist)); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   // compute distance residual at x1 (off the face), get ready for pinging */
/*   ierr = xf_Error(xf_Dist2FaceResidual(3, x1, nf, X, D, PhiData, R, R_xi, &dist)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // ping test */
/*   eps = 1e-5; */
/*   tol = 10.0*eps*eps; */
/*   for (j=0; j<2; j++){ */
/*     xi[j] += eps; */

/*     ierr = xf_Error(xf_EvalBasis(FBasis, 2, xfe_True, 1, xi, xfb_All, &PhiData)); */
/*     if (ierr != xf_OK) return ierr; */

/*     ierr = xf_Error(xf_Dist2FaceResidual(3, x1, nf, X, D, PhiData, R, R_xi, &dist)); */
/*     if (ierr != xf_OK) return ierr; */

/*     for (i=0; i<2; i++){ */
/*       fd = (R[i]-R0[i])/eps; */
/*       an = (R_xi[2*i+j]+R_xi0[2*i+j])*0.5; */
/*       if (fabs(fd-an) > tol){ */
/* 	xf_printf("ping failure[%d, %d]: fd = %.12E, an = %.12E, tol = %.6E\n", */
/* 		  i, j, fd, an, tol); */
/* 	return xf_Error(xf_PING_FAILED); */
/*       } */
/*     } */
/*     xi[j] -= eps; */
/*   } */

/*   // destroy basis */
/*   ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */

/*   return xf_OK;   */
/* } */



/* TEST_xf_CalculateDistFcnElems() */
/* { */
/*   int ierr; */
/*   int nface = 2; */
/*   xf_Face FaceList[2] = {{0,0}, {0,1}}; */
/*   enum xfe_BasisType Basis = xfe_TriLagrange; */
/*   int DistFcnOrder = 1; */
/*   real *Dist; */
/*   real Dist_true[3] = {1,2,1}; */
/*   xf_Vector *DistFcn; */
/*   xf_Vector *ElemChanged; */
/*   xf_Data *D; */
/*   xf_All *All; */

/*   ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False)); */
/*   if (ierr != xf_OK) return ierr; */

/*   if (All->Mesh->nElemGroup != 1) return xf_Error(xf_CODE_LOGIC_ERROR); */

/*   // look for DistFcn vector */
/*   ierr = xf_Error(xf_FindVector(All, "DistFcn", xfe_LinkageGlobElem, 1, */
/* 				NULL, 0, 0, &Basis, &DistFcnOrder, NULL, xfe_SizeReal, */
/* 				xfe_True, xfe_True, &D, &DistFcn, NULL)); */
/*   if (ierr != xf_OK) return ierr; */
/*   D->ReadWrite = xfe_True; */

/*   // set distances to a very large number (will take minimum) */
/*   ierr = xf_Error(xf_SetConstVector(DistFcn, 0, 1e30)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // look for ElemChanged vector (make parallel so can exchange halo) */
/*   ierr = xf_Error(xf_FindVector(All, "ElemChanged", xfe_LinkageGlobElem, 1, */
/* 				NULL, 0, 0, NULL, NULL, NULL, xfe_SizeInt, */
/* 				xfe_True, xfe_False, NULL, &ElemChanged, NULL)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // calculate distances for each elem to faces in FaceList */
/*   ierr = xf_Error(xf_CalculateDistFcnElems(All, nface, FaceList, DistFcn, ElemChanged)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // check distance function at elem 2 */
/*   Dist = DistFcn->GenArray[0].rValue[2]; */
/*   xf_AssertRealVectorWithin(Dist, Dist_true, 3, UTOL3); */
  
/*   // destroy temporary vectors */
/*   ierr = xf_Error(xf_DestroyVector(ElemChanged, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // Destroy All */
/*   ierr = xf_Error(xf_DestroyAll(All)); */
/*   xf_AssertEqual(ierr, xf_OK); */

/*   return xf_OK; */

/* } */


