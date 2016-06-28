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


/******************************************************************/
//   FUNCTION Definition: xf_UnitrMatrix
static int
xf_UnitrMatrix (int n, int r, xf_Matrix **pM)
{
  int ierr;

  ierr = xf_Error(xf_CreateMatrix(pM));
  if (ierr != xf_OK) return ierr;

  // allocate GenArray structure
  ierr = xf_Error(xf_Alloc((void **) &(*pM)->GenArray, 1, sizeof(xf_GenArray)));
  if (ierr != xf_OK) return ierr;

  (*pM)->GenArray[0].n = n;
  (*pM)->GenArray[0].r = r;
  (*pM)->GenArray[0].Size = xfe_SizeReal;
  (*pM)->GenArray[0].vr = NULL;

  ierr = xf_Error(xf_AllocGenArray((*pM)->GenArray));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


TEST_xf_VectorDot()
{
  int ierr, i, r=2;
  int n[1] = {1};
  int **vr0, **vr;
  xf_Vector *A, *B;
  real *rA, *rB;
  real dp, dp_true = 11.0;
  
  ierr = xf_Error(xf_VAlloc2((void ***) &vr0, 1, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  vr0[0][0] = 2;

  for (i=0; i<2; i++){

    vr = (i==1) ? vr0 : NULL;

    ierr = xf_Error(xf_UnitrVector(1, n, &r, vr, &A));
    xf_AssertEqual(ierr, xf_OK);
    
    ierr = xf_Error(xf_UnitrVector(1, n, &r, vr, &B));
    xf_AssertEqual(ierr, xf_OK);
  
    rA = A->GenArray[0].rValue[0];
    rA[0] = 1.0;  rA[1] = 2.0;
    
    rB = B->GenArray[0].rValue[0];
    rB[0] = 3.0;  rB[1] = 4.0;
  
    ierr = xf_Error(xf_VectorDot(A, B, &dp));
    xf_AssertEqual(ierr, xf_OK);
  
    xf_AssertWithin(dp, dp_true, UTOL0);

    ierr = xf_Error(xf_DestroyVector(A, xfe_True));
    xf_AssertEqual(ierr, xf_OK);

    ierr = xf_Error(xf_DestroyVector(B, xfe_True));
    xf_AssertEqual(ierr, xf_OK);

  } // i

  xf_Release2( (void **) vr0);

  return xf_OK;  
}

TEST_xf_VectorMult()
{
  int ierr, i, r=2;
  int n[1] = {1};
  int **vr0, **vr;
  xf_Vector *A;
  real *rA;
  real rtrue[2] = {2.0, 4.0};
  
  ierr = xf_Error(xf_VAlloc2((void ***) &vr0, 1, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  vr0[0][0] = 2;
  
  for (i=0; i<2; i++){

    vr = (i==1) ? vr0 : NULL;
  
    ierr = xf_Error(xf_UnitrVector(1, n, &r, vr, &A));
    xf_AssertEqual(ierr, xf_OK);

    rA = A->GenArray[0].rValue[0];
    rA[0] = 1.0;  rA[1] = 2.0;
    
    ierr = xf_Error(xf_VectorMult(A, 2.0));
    xf_AssertEqual(ierr, xf_OK);
    
    xf_AssertRealVectorWithin(rA, rtrue, 2, UTOL0);
    
    ierr = xf_Error(xf_DestroyVector(A, xfe_True));
    xf_AssertEqual(ierr, xf_OK);
  }    

  xf_Release2( (void **) vr0);

  return xf_OK;  
}

TEST_xf_VectorMultSet()
{
  int ierr, i, r=2;
  int n[1] = {1};
  int **vr0, **vr;
  xf_Vector *A, *B;
  real *rA, *rB;
  real rSet[2] = { 6, 8};
  real rNeg[2] = {-6,-8};
  real rAdd[2] = { 7,10};
  real rSub[2] = {-5,-6};
  
  ierr = xf_Error(xf_VAlloc2((void ***) &vr0, 1, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  vr0[0][0] = 2;
  
  for (i=0; i<2; i++){

    vr = (i==1) ? vr0 : NULL;
  
    ierr = xf_Error(xf_UnitrVector(1, n, &r, vr, &A));
    xf_AssertEqual(ierr, xf_OK);
    
    ierr = xf_Error(xf_UnitrVector(1, n, &r, vr, &B));
    xf_AssertEqual(ierr, xf_OK);
    
    rA = A->GenArray[0].rValue[0];
    rB = B->GenArray[0].rValue[0];
    rB[0] = 3.0;  rB[1] = 4.0;
  
    rA[0] = 1.0;  rA[1] = 2.0;
    ierr = xf_Error(xf_VectorMultSet(B, 2.0, xfe_Set, A));
    xf_AssertEqual(ierr, xf_OK);
    xf_AssertRealVectorWithin(rA, rSet, 2, UTOL0);
    
    rA[0] = 1.0;  rA[1] = 2.0;
    ierr = xf_Error(xf_VectorMultSet(B, 2.0, xfe_Neg, A));
    xf_AssertEqual(ierr, xf_OK);
    xf_AssertRealVectorWithin(rA, rNeg, 2, UTOL0);
    
    rA[0] = 1.0;  rA[1] = 2.0;
    ierr = xf_Error(xf_VectorMultSet(B, 2.0, xfe_Add, A));
    xf_AssertEqual(ierr, xf_OK);
    xf_AssertRealVectorWithin(rA, rAdd, 2, UTOL0);
    
    rA[0] = 1.0;  rA[1] = 2.0;
    ierr = xf_Error(xf_VectorMultSet(B, 2.0, xfe_Sub, A));
    xf_AssertEqual(ierr, xf_OK);
    xf_AssertRealVectorWithin(rA, rSub, 2, UTOL0);
    
    
    ierr = xf_Error(xf_DestroyVector(A, xfe_True));
    xf_AssertEqual(ierr, xf_OK);
    
    ierr = xf_Error(xf_DestroyVector(B, xfe_True));
    xf_AssertEqual(ierr, xf_OK);

  }

  xf_Release2( (void **) vr0);

  return xf_OK;  
}



TEST_xf_VectorNorm()
{
  int ierr, i, r=2;
  int n[1] = {1};
  int **vr0, **vr;
  xf_Vector *A;
  real *rA, rnorm;
  
  ierr = xf_Error(xf_VAlloc2((void ***) &vr0, 1, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  vr0[0][0] = 2;
  
  for (i=0; i<2; i++){

    vr = (i==1) ? vr0 : NULL;
  
    ierr = xf_Error(xf_UnitrVector(1, n, &r, vr, &A));
    xf_AssertEqual(ierr, xf_OK);
    
    rA = A->GenArray[0].rValue[0];
    rA[0] = 3.0;  rA[1] = 4.0;
    
    ierr = xf_Error(xf_VectorNorm(A, 2, &rnorm));
    xf_AssertEqual(ierr, xf_OK);
    
    xf_AssertWithin(rnorm, 5.0, UTOL0);
    
    ierr = xf_Error(xf_DestroyVector(A, xfe_True));
    xf_AssertEqual(ierr, xf_OK);

  }

  xf_Release2( (void **) vr0);

  return xf_OK;  
}

TEST_xf_MatrixInvert()
{
  int ierr, i;
  xf_Matrix *M, *iM;
  real rM[4] = {2,4,1,3};
  real riM[4] = {1.5,-2,-.5,1};
  
  ierr = xf_Error(xf_UnitrMatrix(1, 4, &M));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_UnitrMatrix(1, 4, &iM));
  xf_AssertEqual(ierr, xf_OK);

  for (i=0; i<4; i++) M->GenArray[0].rValue[0][i] = rM[i];

  ierr = xf_Error(xf_MatrixInvert(M, 2, iM));
  xf_AssertEqual(ierr, xf_OK);
  
  xf_AssertRealVectorWithin(iM->GenArray[0].rValue[0], riM, 4, UTOL0);

  ierr = xf_Error(xf_DestroyMatrix(M));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_DestroyMatrix(iM));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}

TEST_xf_VectorSetSVD()
{
  int ierr, i, j, k, l, i0, i1;
  int n = 2, r[2] = {2,2};
  int nelem[2] = {1, 2};
  xf_Vector *V;
  xf_VectorSet *VS;
  real *rV;
  real W[2], W_true[2], U[12];
  real U_true[12] = {  0.408248290463863,  0.000000000000000, 
		       0.612372435695795,  0.500000000000000, 
		       0.612372435695795, -0.500000000000000, 
		       0.204124145231932, -0.500000000000000, 
		      -0.204124145231932, -0.500000000000000, 
		       0,                   0 };
  
  // create a vectorset
  ierr = xf_Error(xf_UnitrVectorSet(n, 2, nelem, r, &VS));
  xf_AssertEqual(ierr, xf_OK);

  V = VS->Vector + 0;
  rV = V->GenArray[0].rValue[0];  rV[0] =  1.0;  rV[1] =  2.0;
  rV = V->GenArray[1].rValue[0];  rV[0] =  1.0;  rV[1] =  0.0;
  rV = V->GenArray[1].rValue[1];  rV[0] = -1.0;  rV[1] =  0.0;

  V = VS->Vector + 1;
  rV = V->GenArray[0].rValue[0];  rV[0] = 1.0;  rV[1] =  1.0;
  rV = V->GenArray[1].rValue[0];  rV[0] = 2.0;  rV[1] =  1.0;
  rV = V->GenArray[1].rValue[1];  rV[0] = 0.0;  rV[1] =  0.0;

  ierr = xf_Error(xf_VectorSetSVD(VS, n, W));
  xf_AssertEqual(ierr, xf_OK);

  i0 = 0; i1 = 1;
  if (W[0] < W[1]) swap(i0, i1, i);

  W_true[i0] = sqrt(12.); W_true[i1] = sqrt(2);
  xf_AssertRealVectorWithin(W, W_true, 2, UTOL1);

  // unroll left singular vectors
  V = VS->Vector + 0;
  k = 0;
  for (i=0; i<2; i++)
    for (j=0; j<V->GenArray[i].n; j++)
      for (l=0; l<V->GenArray[i].r; l++)
	U[2*(k++)+i0] = V->GenArray[i].rValue[j][l];
  V = VS->Vector + 1;
  k = 0;
  for (i=0; i<2; i++)
    for (j=0; j<V->GenArray[i].n; j++)
      for (l=0; l<V->GenArray[i].r; l++)
	U[2*(k++)+i1] = V->GenArray[i].rValue[j][l];

  xf_AssertRealVectorWithin(U, U_true, 12, UTOL1);

  ierr = xf_Error(xf_DestroyVectorSet(VS));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_ProjectVector()
{
  int ierr, i, j, sr = 1;
  int n[2] = {2, 3}, r[2] = {1,1};
  enum xfe_BasisType Basis[2] = {xfe_TriLagrange, xfe_QuadLagrange};
  int OrderA[2] = {2, 3};
  int OrderB[2] = {3, 4};
  int **vOrderA = NULL, **vOrderA0, vOrderA1[2][3] = {{1,2,-1},{2,3,1}};
  int **vOrderB = NULL, **vOrderB0, vOrderB1[2][3] = {{2,3,-1},{3,4,1}};
  int VOrder0[2][3] = {{2,4,-1},{4,4,3}};
  xf_Vector *A, *A0, *B, *VOrder;
  real rnorm;

  ierr = xf_Error(xf_VAlloc2((void ***) &vOrderA0, 2, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VAlloc2((void ***) &vOrderB0, 2, n, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<2; i++)
    for (j=0; j<n[i]; j++){
      vOrderA0[i][j] = vOrderA1[i][j];
      vOrderB0[i][j] = vOrderB1[i][j];
    }

  // allocate VOrder
  ierr = xf_Error(xf_UnitiVector(2, n, r, NULL, &VOrder));
  xf_AssertEqual(ierr, xf_OK);
  for (i=0; i<2; i++)
    for (j=0; j<n[i]; j++)
      VOrder->GenArray[i].iValue[j][0] = VOrder0[i][j];

  
  for (i=0; i<2; i++){

    vOrderA = (i==1) ? vOrderA0 : NULL;
    vOrderB = (i==1) ? vOrderB0 : NULL;

    // create A
    ierr = xf_Error(xf_UnitrVectorInterp(2, n, sr, Basis, OrderA, vOrderA, &A));
    xf_AssertEqual(ierr, xf_OK);
    
    // Fill A with pseudo-random values
    ierr = xf_Error(xf_VectorRand(A, 3));
    xf_AssertEqual(ierr, xf_OK);
    
    // copy A to A0
    ierr = xf_Error(xf_CreateVector(&A0));
    xf_AssertEqual(ierr, xf_OK);

    ierr = xf_Error(xf_CopyVector(NULL, A, A0));
    xf_AssertEqual(ierr, xf_OK);  
    
    // create B
    ierr = xf_Error(xf_UnitrVectorInterp(2, n, sr, Basis, OrderB, vOrderB, &B));
    xf_AssertEqual(ierr, xf_OK);
    
    // project A to B (B is higher-order, no losses expected)
    ierr = xf_Error(xf_ProjectVector( NULL, A, xfe_False, B));
    xf_AssertEqual(ierr, xf_OK);
    
    // project B back to A
    ierr = xf_Error(xf_ProjectVector( NULL, B, xfe_False, A));
    xf_AssertEqual(ierr, xf_OK);

    // check if A == A0
    ierr = xf_Error(xf_SetVector(A0, xfe_Sub, A));
    xf_AssertEqual(ierr, xf_OK);
    ierr = xf_Error(xf_VectorNorm(A, 2, &rnorm));
    xf_AssertEqual(ierr, xf_OK);
    xf_AssertWithin(rnorm, 0.0, UTOL2);

    // set A to A0
    ierr = xf_Error(xf_SetVector(A0, xfe_Set, A));
    xf_AssertEqual(ierr, xf_OK);

    // project A in place using B as a template
    ierr = xf_Error(xf_ProjectVectorInPlace_Vector(NULL, NULL, A, B));
    xf_AssertEqual(ierr, xf_OK);
    
    // project A back to original order (use A0 as template)
    ierr = xf_Error(xf_ProjectVectorInPlace_Vector(NULL, NULL, A, A0));
    xf_AssertEqual(ierr, xf_OK);
    
    // project A to (v)OrderB using OrderSet
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderSet(NULL, NULL, A, NULL, xfe_BasisLast, 
						     vOrderB, OrderB, -1));
    xf_AssertEqual(ierr, xf_OK);
    
    // project A back to original order (use A0 as template)
    ierr = xf_Error(xf_ProjectVectorInPlace_Vector(NULL, NULL, A, A0));
    xf_AssertEqual(ierr, xf_OK);
    
    // project A to +1 order using OrderIncrement
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(NULL, NULL, A, NULL, 
							   xfe_BasisLast, 1));
    xf_AssertEqual(ierr, xf_OK);
    
    // project A to -1 order using OrderIncrement (back to original)
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(NULL, NULL, A, NULL,
							   xfe_BasisLast, -1));
    xf_AssertEqual(ierr, xf_OK);
    
    // project A using VOrder
    ierr = xf_Error(xf_ProjectVectorInPlace_VOrder(NULL, NULL, A, NULL, xfe_BasisLast, 
						   VOrder));
    xf_AssertEqual(ierr, xf_OK);
    
    // project A back to original order (use A0 as template)
    ierr = xf_Error(xf_ProjectVectorInPlace_Vector(NULL, NULL, A, A0));
    xf_AssertEqual(ierr, xf_OK);
    
    // make A0 look like a variable order vector so we can compare it to A
    ierr = xf_Error(xf_PrepVectorVariableOrder(A0));
    if (ierr != xf_OK) return ierr;

    // check that A == A0
    ierr = xf_Error(xf_SetVector(A0, xfe_Sub, A));
    xf_AssertEqual(ierr, xf_OK);
    ierr = xf_Error(xf_VectorNorm(A, 2, &rnorm));
    xf_AssertEqual(ierr, xf_OK);
    xf_AssertWithin(rnorm, 0.0, UTOL2);
    
    // destroy vectors
    ierr = xf_Error(xf_DestroyVector(A,  xfe_True));
    xf_AssertEqual(ierr, xf_OK);
    ierr = xf_Error(xf_DestroyVector(A0, xfe_True));
    xf_AssertEqual(ierr, xf_OK);
    ierr = xf_Error(xf_DestroyVector(B,  xfe_True));
    xf_AssertEqual(ierr, xf_OK);

  } // i

  ierr = xf_Error(xf_DestroyVector(VOrder,  xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  xf_Release2( (void **) vOrderA0);
  xf_Release2( (void **) vOrderB0);

  return xf_OK;  
}




