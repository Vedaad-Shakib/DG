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

#include <limits.h>
extern "C" {
#include "xf_Unit.h"
#include "xf_All.h"
#include "xf_Basis.h"
#include "xf_Data.h"
#include "xf_Math.h"
}
#include "gtest/gtest.h"



TEST(MxV_Add_Test, All) {
    EXPECT_EQ(1, 1);
}

/*
TEST_xf_PowInt()
{
  real y;
  
  y = xf_PowInt(3.0, 4);
  xf_AssertWithin(y, 81.0, UTOL0);

  return xf_OK;
}

TEST_xf_MxV()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real u[2] = {2, 1};
  real v[3] = {-1, 0, 2};
  real v0[3] = {3, 10, 18};
  real v1[3] = {-1, 0, 2};
  real v2[3] = { 4,  10,  16};
  real v3[3] = {-4, -10, -16};

  xf_MxV(A, u, 3, 2, xfe_Add, v);
  xf_AssertRealVectorWithin(v, v0, 3, UTOL0);

  xf_MxV(A, u, 3, 2, xfe_Sub, v);
  xf_AssertRealVectorWithin(v, v1, 3, UTOL0);

  xf_MxV(A, u, 3, 2, xfe_Set, v);
  xf_AssertRealVectorWithin(v, v2, 3, UTOL0);

  xf_MxV(A, u, 3, 2, xfe_Neg, v);
  xf_AssertRealVectorWithin(v, v3, 3, UTOL0);

  return xf_OK;  
}

TEST_xf_cMxV_Add()
{
  real c = 2.0;
  real A[6] = {1, 2, 3, 4, 5, 6};
  real u[2] = {2, 1};
  real v[3] = {-1, 0, 2};
  real v_true[3] = {7, 20, 34};

  xf_cMxV_Add(c, A, u, 3, 2, v);
  xf_AssertRealVectorWithin(v, v_true, 3, UTOL0);

  return xf_OK;  
}


TEST_xf_MTxV()
{
  real A[6] = {1, 3, 5, 2, 4, 6};
  real u[2] = {2, 1};
  real v[3] = {-1, 0, 2};
  real v0[3] = {3, 10, 18};
  real v1[3] = {-1, 0, 2};

  xf_MTxV_Add(A, u, 3, 2, v);
  xf_AssertRealVectorWithin(v, v0, 3, UTOL0);

  xf_MTxV_Sub(A, u, 3, 2, v);
  xf_AssertRealVectorWithin(v, v1, 3, UTOL0);

  return xf_OK;  
}

TEST_xf_MxM()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 2, 1, 3, 1, 2};
  real C[4] = {2, 3, -1, -10};
  real C0[4] = {-4, -11, -16, -45};
  real C1[4] = {2, 3, -1, -10};
  real C2[4] = {6, 14, 15, 35};
  real C3[4] = {-6, -14, -15, -35};

  xf_MxM(A, B, 2, 3, 2, xfe_Sub, C);
  xf_AssertRealVectorWithin(C, C0, 4, UTOL0);

  xf_MxM(A, B, 2, 3, 2, xfe_Add, C);
  xf_AssertRealVectorWithin(C, C1, 4, UTOL0);

  xf_MxM(A, B, 2, 3, 2, xfe_Set, C);
  xf_AssertRealVectorWithin(C, C2, 4, UTOL0);

  xf_MxM(A, B, 2, 3, 2, xfe_Neg, C);
  xf_AssertRealVectorWithin(C, C3, 4, UTOL0);

  return xf_OK;  
}

TEST_xf_nMxM()
{
  real A[12] = { 1, 2, 3,  4, 5, 6,
	       -1, 2, 3, -4, 5, 6};
  real B[12] = {1, 2, 1, 3, 1, 2,
	       2, 0, 0, 1, 3, -1};
  real C[8];
  real C0[8] = {6, 14, 15, 35, 
		7, -1, 10, -1};
  real C1[8] = {12, 28, 30, 70, 
		14, -2, 20, -2};
  real Ap[8] = {1,2,3,4, -1,-2,-3,-4};
  real Bp[4] = {0,1,2,3};
  real Cp[4];
  real Cp0[4] = {2,4,-8,-18};

  xf_nMxM_Set(2, A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C0, 8, UTOL0);

  xf_nMxM_Add(2, A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C1, 8, UTOL0);

  xf_nMxM_Set(2, Ap, Bp, 2, 2, 1, Cp);
  xf_AssertRealVectorWithin(Cp, Cp0, 4, UTOL0);


  return xf_OK;  
}

TEST_xf_cMxM_Add()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 2, 1, 3, 1, 2};
  real C[4] = {2, 3, -1, -10};
  real C0[4] = {-1, -4, -8.5, -27.5};

  xf_cMxM_Add(-0.5, A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C0, 4, UTOL0);

  return xf_OK;  
}


TEST_xf_MTxM()
{
  real A[6] = {1, 4, 2, 5, 3, 6};
  real B[6] = {1, 2, 2, 1, 3, 1};
  real C[4] = {1, 1, 0, 0};
  real C0[4] = {15, 8, 32, 19};
  real C1[4] = {1, 1, 0, 0};
  real C2[4] = {14, 7, 32, 19};
  real C3[4] = {-14, -7, -32, -19};

  xf_MTxM(A, B, 2, 3, 2, xfe_Add, C);
  xf_AssertRealVectorWithin(C, C0, 4, UTOL0);

  xf_MTxM(A, B, 2, 3, 2, xfe_Sub, C);
  xf_AssertRealVectorWithin(C, C1, 4, UTOL0);

  xf_MTxM(A, B, 2, 3, 2, xfe_Set, C);
  xf_AssertRealVectorWithin(C, C2, 4, UTOL0);

  xf_MTxM(A, B, 2, 3, 2, xfe_Neg, C);
  xf_AssertRealVectorWithin(C, C3, 4, UTOL0);

  return xf_OK;  
}


TEST_xf_MTxwM()
{
  real A[6] = {1, 4, 2, 5, 3, 6};
  real w[3] = {3, 4, 5};
  real B[6] = {1, 2, 2, 1, 3, 1};
  real C[4];
  real C0[4] = {64, 29, 142, 74};
  real C1[4] = {128, 58, 284, 148};

  xf_MTxwM_Set(A, w, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C0, 4, UTOL0);

  xf_MTxwM_Add(A, w, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C1, 4, UTOL0);

  return xf_OK;  
}


TEST_xf_MxMT_Set()
{
  real A[6] = {1, 4, 2, 5, 3, 6};
  real B[6] = {1, 2, 2, 1, 3, 1};
  real C[4];
  real C0[4] = {13, 15, 23, 20};

  xf_MxMT_Set(A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C0, 4, UTOL0);

  return xf_OK;  
}



TEST_xf_ColMult()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real A_true[6] = {1, 2, 9, 12, 25, 30};

  xf_ColMult(A, v, 3, 2, 2);
  xf_AssertRealVectorWithin(A, A_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_ColMult_Set()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {1, 2, 9, 12, 25, 30}, B[6];

  xf_ColMult_Set(A, v, 3, 2, 2, B);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;  
}


TEST_xf_ColMult_Add()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 1, 0, 1, 0, 1};
  real B_true[6] = {2, 3, 9, 13, 25, 31};

  xf_ColMult_Add(A, v, 3, 2, 2, B);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_ColMult_Sub()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 1, 0, 1, 0, 1};
  real B_true[6] = {0, -1, -9, -11, -25, -29};

  xf_ColMult_Sub(A, v, 3, 2, 2, B);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;  
}


TEST_xf_ColcMult()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real A_true[6] = {1.5, 3.0, 13.5, 18, 37.5, 45};

  xf_ColcMult(A, v, 3, 2, 2, 1.5);
  xf_AssertRealVectorWithin(A, A_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_ColcMult_Set()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {1.5, 3.0, 13.5, 18, 37.5, 45}, B[6];

  xf_ColcMult_Set(A, v, 3, 2, 2, 1.5, B);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_ColcMult_Add()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {0.5, 2.0, 12.5, 17, 36.5, 44};
  real B[6] = {-1, -1, -1, -1, -1, -1};

  xf_ColcMult_Add(A, v, 3, 2, 2, 1.5, B);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_2ColMult_Set()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real w[6] = {2, 1, 4, 3, 0, 5};
  real B_true[6] = {2, 4, 36, 48, 0, 0}, B[6];

  xf_2ColMult_Set(A, v, w, 3, 2, 2, 2, B);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_2ColcMult_Set()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real w[6] = {2, 1, 4, 3, 0, 5};
  real B_true[6] = {1, 2, 18, 24, 0, 0}, B[6];

  xf_2ColcMult_Set(A, v, w, 3, 2, 2, 2, 0.5, B);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_2ColcMult_Add()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real w[6] = {2, 1, 4, 3, 0, 5};
  real B_true[6] = {0, 0, 15, 20, -5, -6};
  real B[6] = {-1, -2, -3, -4, -5, -6};

  xf_2ColcMult_Add(A, v, w, 3, 2, 2, 2, 0.5, B);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_MatDetInv()
{
  int ierr;
  real A[4] = {1,2,3,4};
  real iA[9], detA;
  real iA_true[4] = {-2,1,1.5,-0.5};
  real B[9] = {1,0,-1,0,1,1,1,2,0};
  real iB[9], detB;
  real iB_true[9] = {2,2,-1,-1,-1,1,1,2,-1};

  ierr = xf_Error(xf_MatDetInv(A, 2, &detA, iA));
  xf_AssertEqual(ierr, 0);
  xf_AssertWithin(detA, -2.0, 1e-15);
  xf_AssertRealVectorWithin(iA, iA_true, 4, UTOL0);

  ierr = xf_Error(xf_MatDetInv(B, 3, &detB, iB));
  xf_AssertEqual(ierr, 0);
  xf_AssertWithin(detB, -1.0, 1e-15);
  xf_AssertRealVectorWithin(iB, iB_true, 9, UTOL0);

  return xf_OK;  
}

TEST_xf_DotProduct()
{
  real a[3] = {1,2,3};
  real b[3] = {4,5,6};
  real c;
  
  xf_DotProduct(a, b, 3, &c);
  xf_AssertWithin(c, 32.0, 1e-15);

  return xf_OK;  
}


TEST_xf_CrossProduct()
{
  real a[3] = {1,2,3};
  real b[3] = {4,5,6};
  real c_true[3] = {-3,6,-3}, c[3];
  
  xf_CrossProduct(a, b, c);
  xf_AssertRealVectorWithin(c, c_true, 3, UTOL0);

  return xf_OK;  
}

TEST_xf_ComputeSolvePLUT()
{
  int ierr, P[3];
  real A[9] = {1, 2, -1, 1, -1, -2, 0, 2, 1};
  real b[3] = {1, 2, 3};
  real u0[3] = {20.0, -4.0, 11.0}, u[3];
  real u1[3] = {-7.0, 8.0, 12.0};

  ierr = xf_Error(xf_ComputePLU(A, 3, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_SolvePLU(A, P, b, 3, u, NULL));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(u, u0, 3, UTOL1);

  ierr = xf_Error(xf_SolvePLUT(A, P, b, 3, u, NULL));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(u, u1, 3, UTOL1);

  return xf_OK;
}


TEST_xf_SolvePLU_Matrix()
{
  int ierr, P[2];
  real A[4] = {1, 2, 3, 4};
  real B[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {2, 1, 0, -.5, .5, 1.5};

  ierr = xf_Error(xf_ComputePLU(A, 2, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_SolvePLU_Matrix(A, P, 2, 3, B));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;
}

TEST_xf_SolvePLU_MatrixR()
{
  int ierr, P[2];
  real A[4] = {1, 2, 3, 4};
  real B[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {1, 0, 0, 1, -1, 2};

  ierr = xf_Error(xf_ComputePLU(A, 2, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_SolvePLU_MatrixR(A, P, 2, 3, B));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(B, B_true, 6, UTOL0);

  return xf_OK;
}

TEST_xf_ComputeSolveBlockPLUT()
{
  int ierr, k, P[4];
  real A[16] =  {1,2,4,3, 3,4,2,1, 1,3,3,2, 2,4,1,4};
  real b[4] = {1, 6, 3, 4};
  real c[4] = {9, 3, 6, 1};
  real u0[4] = {0.85, 1.35, -0.65, -0.15}, u[4];
  real u1[4] = {2.45, 2.2, -3.0, 0.25};

  ierr = xf_Error(xf_ComputeBlockPLU(A, 2, 2, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_SolveBlockPLU(A, 2, 2, P, b, xfe_Set, u));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(u, u0, 4, UTOL0);

  ierr = xf_Error(xf_SolveBlockPLUT(A, 2, 2, P, c, xfe_Set, u));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(u, u1, 4, UTOL0);

  return xf_OK;
}

TEST_xf_ComputeSolveBlockPLU2()
{
  int ierr, k, P[4];
  real A[16] =  {1,2,3,4, 4,3,2,1, 1,3,2,4, 3,2,1,4};
  real b[4] = {1, 6, 3, 4};
  real u_true[4] = {-0.85, -1.35, 0.65, 0.15}, u[4];

  ierr = xf_Error(xf_ComputeBlockPLU(A, 1, 4, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_SolveBlockPLU(A, 1, 4, P, b, xfe_Neg, u));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(u, u_true, 4, UTOL0);

  return xf_OK;
}


TEST_xf_ComputeSolveBlockPLU_Matrix()
{
  int ierr, k, P[4];
  real A[16] =  {1,2,4,3, 3,4,2,1, 1,3,3,2, 2,4,1,4};
  real B[8] = {1, 2, 6, 12, 3, 6, 4, 8};
  real U_true[8] = {0.85, 1.7, 1.35, 2.7, -0.65, -1.3, -0.15, -.3}, U[8];
  real iA[16], I[16];

  ierr = xf_Error(xf_ComputeBlockPLU(A, 2, 2, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_SolveBlockPLU_Matrix(A, 2, 2, P, B, 1, xfe_Set, U));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(U, U_true, 8, UTOL0);

  ierr = xf_Error(xf_SolveBlockPLU_Matrix(A, 2, 2, P, B, 1, xfe_Sub, U));
  ierr = xf_Error(xf_SolveBlockPLU_Matrix(A, 2, 2, P, B, 1, xfe_Add, U));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(U, U_true, 8, UTOL0);

  // BlockIdentity is tested here
  xf_BlockIdentity(2, 2, I);
  ierr = xf_Error(xf_SolveBlockPLU_Matrix(A, 2, 2, P, I, 2, xfe_Set, iA));
  xf_AssertEqual(ierr, 0);
  xf_BlockMxBlockM(iA, 2, 2, 2, B, 1, xfe_Set, U);
  xf_AssertRealVectorWithin(U, U_true, 8, UTOL0);

  return xf_OK;
}

TEST_xf_PLUMxV_Set()
{
  int ierr, P[2];
  real A[4] = {1, 2, 3, 4};
  real u[2] = {1, 2};
  real b_true[2] = {5, 11}, b[2];

  ierr = xf_Error(xf_ComputePLU(A, 2, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_PLUMxV_Set(A, P, u, 2, b, NULL));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(b, b_true, 2, UTOL0);

  return xf_OK;
}

TEST_xf_PLUMTxV_Set()
{
  int ierr, P[3];
  real A[9] = {1, 2, -1, 1, -1, -2, 0, 2, 1};
  real u[3] = {1, 2, 3};
  real b0[3] = {3.0, 6.0, -2.0}, b[3];

  ierr = xf_Error(xf_ComputePLU(A, 3, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_PLUMTxV_Set(A, P, u, 3, b, NULL));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(b, b0, 3, UTOL0);

  return xf_OK;
}

TEST_xf_BlockPLUMTxV()
{
  int ierr, k, P[6];
  real A[36] = {1,2,6,5, 3,4,4,3, 5,6,2,1, 
		1,3,2,5, 2,4,3,1, 5,6,4,6, 
		4,3,2,5, 2,5,3,2, 6,1,6,1};
  real u[6] = {1, 2, 3, 4, 5, 6};
  real b0[6] = {91, 56, 90, 81, 72, 65}, b[6];
  real b1[6] = {56, 86, 57, 63, 106, 61};

  ierr = xf_Error(xf_ComputeBlockPLU(A, 3, 2, P));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_BlockPLUMxV(A, 3, 2, P, u, xfe_Set, b));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(b, b0, 6, UTOL2);

  ierr = xf_Error(xf_BlockPLUMTxV(A, 3, 2, P, u, xfe_Set, b));
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(b, b1, 6, UTOL2);

  return xf_OK;
}

TEST_xf_BlockMxV()
{
  real A[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  real u[4] = {2, 1, -1, -2};
  real v[2] = {0, 0};
  real z[2] = {0, 0};
  real v_true[2] = {-13, -13};

  xf_BlockMxV(A, 1, 2, 2, u, xfe_Set, v);
  xf_AssertRealVectorWithin(v, v_true, 2, UTOL0);

  xf_BlockMxV(A, 1, 2, 2, u, xfe_Sub, v);
  xf_AssertRealVectorWithin(v, z, 2, UTOL0);

  xf_BlockMxV(A, 1, 2, 2, u, xfe_Add, v);
  xf_AssertRealVectorWithin(v, v_true, 2, UTOL0);

  return xf_OK;
}

TEST_xf_BlockMTxV()
{
  real A[8] = {1, 3, 2, 4, 5, 7, 6, 8};
  real u[4] = {2, 1, -1, -2};
  real v[2] = {0, 0};
  real z[2] = {0, 0};
  real v_true[2] = {-13, -13};

  xf_BlockMTxV(A, 1, 2, 2, u, xfe_Set, v);
  xf_AssertRealVectorWithin(v, v_true, 2, UTOL0);

  xf_BlockMTxV(A, 1, 2, 2, u, xfe_Sub, v);
  xf_AssertRealVectorWithin(v, z, 2, UTOL0);

  xf_BlockMTxV(A, 1, 2, 2, u, xfe_Add, v);
  xf_AssertRealVectorWithin(v, v_true, 2, UTOL0);

  return xf_OK;
}


TEST_xf_BlockMxM()
{
  real A[32] = {1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7, 
		5,6,7,8, 6,7,8,9, 7,8,9,10, 8,9,10,11};
  real U[8] = {1,0, 0,1, 1,0, 0,1};
  real B_true[16] = {4,6,8,10, 6,8,10,12, 12,14,16,18, 14,16,18,20};
  real B[16] = {0};
  real Z[16] = {0};

  xf_BlockMxM(A, 2, 2, 4, U, 2, xfe_Set, B);
  xf_AssertRealVectorWithin(B, B_true, 16, UTOL0);

  xf_BlockMxM(A, 2, 2, 4, U, 2, xfe_Sub, B);
  xf_AssertRealVectorWithin(B, Z, 16, UTOL0);

  xf_BlockMxM(A, 2, 2, 4, U, 2, xfe_Neg, B);
  xf_BlockMxM(A, 2, 2, 4, U, 2, xfe_Add, B);
  xf_AssertRealVectorWithin(B, Z, 16, UTOL0);

  return xf_OK;
}


TEST_xf_MxBlockM()
{
  real A[32] = {1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7, 
		5,6,7,8, 6,7,8,9, 7,8,9,10, 8,9,10,11};
  real U[2] = {1, 1};
  real B_true[16] = {6,8,10,12, 8,10,12,14, 10,12,14,16, 12,14,16,18};
  real B[16] = {0};
  real Z[16] = {0};

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Set, B);
  xf_AssertRealVectorWithin(B, B_true, 16, UTOL0);

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Sub, B);
  xf_AssertRealVectorWithin(B, Z, 16, UTOL0);

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Neg, B);
  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Add, B);
  xf_AssertRealVectorWithin(B, Z, 16, UTOL0);

  return xf_OK;
}


TEST_xf_MTxBlockM()
{
  real A[32] = {1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7, 
		5,6,7,8, 6,7,8,9, 7,8,9,10, 8,9,10,11};
  real U[2] = {1, 1};
  real B_true[16] = {6,8,10,12, 8,10,12,14, 10,12,14,16, 12,14,16,18};
  real B[16] = {0};
  real Z[16] = {0};

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Set, B);
  xf_AssertRealVectorWithin(B, B_true, 16, UTOL0);

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Sub, B);
  xf_AssertRealVectorWithin(B, Z, 16, UTOL0);

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Neg, B);
  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Add, B);
  xf_AssertRealVectorWithin(B, Z, 16, UTOL0);

  return xf_OK;
}

TEST_xf_BlockMxBlockM()
{
  real A[32] = {1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7, 
		5,6,7,8, 6,7,8,9, 7,8,9,10, 8,9,10,11};
  real U[48] = {4,3,2,1, 5,4,3,2, 6,5,4,3, 
		7,6,5,4, 8,7,6,5, 9,8,7,6,
		3,2,1,4, 4,3,2,5, 5,4,3,6,
		6,5,4,7, 7,6,5,8, 8,7,6,9};
  real B_true[24] = {94 ,  106,  158, 170, 118, 130,
		     198,  210,  142, 154, 238, 250,
		     222,  234,  286, 298, 278, 290,
		     358,  370,  334, 346, 430, 442};
  real B[24] = {0};
  real Z[24] = {0};

  xf_BlockMxBlockM(A, 2, 2, 4, U, 3, xfe_Set, B);
  xf_AssertRealVectorWithin(B, B_true, 24, UTOL0);

  xf_BlockMxBlockM(A, 2, 2, 4, U, 3, xfe_Sub, B);
  xf_AssertRealVectorWithin(B, Z, 24, UTOL0);

  xf_BlockMxBlockM(A, 2, 2, 4, U, 3, xfe_Neg, B);
  xf_BlockMxBlockM(A, 2, 2, 4, U, 3, xfe_Add, B);
  xf_AssertRealVectorWithin(B, Z, 24, UTOL0);

  return xf_OK;
}

TEST_xf_OutProd_Add()
{
  real u[3] = {1, 2};
  real v[3] = {4, -5};
  real A[4] = {-1, -2, -3, -4};
  real A_true[4] = {3, -7, 5, -14};

  xf_OutProd_Add(u, v, 2, 2, A);
  xf_AssertRealVectorWithin(A, A_true, 4, UTOL0);

  return xf_OK;  
}

TEST_xf_OutProd_Sub()
{
  real u[2] = {1, 2};
  real v[2] = {4, -5};
  real A[4] = {-1, -2, -3, -4};
  real A_true[4] = {-5, 3, -11, 6};

  xf_OutProd_Sub(u, v, 2, 2, A);
  xf_AssertRealVectorWithin(A, A_true, 4, UTOL0);

  return xf_OK;  
}

TEST_xf_BlockOutProd_Add()
{
  real u[2] = {-2, 3};
  real v[4] = {1, 2, 3, 4};
  real A[8] = {0, 0, -1, 0, 0, 0, 0, 1};
  real A_true[8] = {-2, -4, 2, 6, -6, -8, 9, 13};

  xf_BlockOutProd_Add(u, v, 1, 2, 2, A);
  xf_AssertRealVectorWithin(A, A_true, 8, UTOL0);

  return xf_OK;  
}

TEST_xf_BlockOutProd_Sub()
{
  real u[2] = {-2, 3};
  real v[4] = {1, 2, 3, 4};
  real A[8] = {0, 0, 0, 0, 0, 1, 0, 1};
  real A_true[8] = {2, 4, -3, -6, 6, 9, -9, -11};

  xf_BlockOutProd_Sub(u, v, 1, 2, 2, A);
  xf_AssertRealVectorWithin(A, A_true, 8, UTOL0);

  return xf_OK;  
}


TEST_xf_QRFactorHouseholder()
{
  int ierr;
  real A[6] = {2.4, 1.8, 3.2, 2.4, 0.0, 2.0};
  real Q_true[9] = {-.6, 0, -.8, -.8, 0, 0.6, 0, 1.0, 0};
  real R_true[6] = {-4.0, -3.0, 0, 2.0, 0, 0};
  real Q[9], R[6];
  
  ierr = xf_QRFactorHouseholder(A, 3, 2, Q, R);
  xf_AssertEqual(ierr, 0);

  xf_AssertRealVectorWithin(Q, Q_true, 9, UTOL0);
  xf_AssertRealVectorWithin(R, R_true, 6, UTOL0);

 return xf_OK;
}


TEST_xf_EigSymmetricQR()
{
  int ierr;
  real tol = 1e-8;
  real A[9] = {1.64, -0.48, 0, -0.48, 1.36, 0, 0, 0, 3.0};
  real EG[3], EG_true[3] = {3,2,1};
  real EV[9], EV_true[9] = {0,0,-1, -0.8,0.6,0, -0.6,-0.8,0};
  
  ierr = xf_EigSymmetricQR(A, 3, tol, 100, EG, EV);
  xf_AssertEqual(ierr, 0);

  xf_AssertRealVectorWithin(EG, EG_true, 3, 10.*tol);
  xf_AssertRealVectorWithin(EV, EV_true, 9, 10.*tol);

 return xf_OK;
}

TEST_xf_SVDGolubReinsch()
{
  int ierr, i;
  real A[6] = {1,1,2,1,1,2};
  real U_true[6] = {-0.426401432711221,  0.0,
		    -0.639602149066831, -0.707106781186548,
		    -0.639602149066831,  0.707106781186547};
  real B[25] = {0};
  real W[5];
  real V[25], V_true[4];
  real W_true[5], B_true[25], C[25];

  W_true[0] = sqrt(11.); W_true[1] = 1.;
  V_true[0] = -0.5*sqrt(2.); V_true[1] = -0.5*sqrt(2.);
  V_true[2] = -0.5*sqrt(2.); V_true[3] =  0.5*sqrt(2.);
  ierr = xf_SVDGolubReinsch(A, 3, 2, xfe_True, W, V);
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(W, W_true, 2, UTOL1);
  xf_AssertRealVectorWithin(V, V_true, 4, UTOL1);
  xf_AssertRealVectorWithin(A, U_true, 6, UTOL1);

  for (i=0; i<5; i++){
    if (i > 0) B[6*i-1] = 1;
    if (i < 4) B[6*i+1] = 1;
    B[6*i] = 2;
  }
  for (i=0; i<25; i++) B_true[i] = B[i];
  W_true[0] = 3.732050807568878;
  W_true[1] = 3.0; W_true[2] = 2.0; W_true[3] = 1.0;
  W_true[4] = 0.267949192431122;
  ierr = xf_SVDGolubReinsch(B, 5, 5, xfe_True, W, V);
  xf_AssertEqual(ierr, 0);
  xf_AssertRealVectorWithin(W, W_true, 5, UTOL1);
  
  for (i=0; i<25; i++) B[i] *= W[i%5]; // B = B*diag(W)
  xf_MxMT_Set(B, V, 5, 5, 5, C); // C = B*V^T
  xf_AssertRealVectorWithin(C, B_true, 25, UTOL1);
       
  return xf_OK;
}


TEST_xf_QRBulgeChase()
{
  int ierr, k, np, nk;
  real shift[2] = {0.0, 1.0};
  real alpha[5] = {2,2,2,2,2};
  real beta[5] = {17, 1,1,1,1}; // first entry is irrelevant
  real Q[25];
  real alpha0[5] = { 3.263157894736842E+00, 2.826729745712597E+00,
		     2.453590620420127E+00, 4.565217391304347E-01,
		     1.000000000000000E+00}; // from Matlab
  real beta0[5] = {1.700000000000000E+01, 4.965253227398214E-01,
		   6.643481123188533E-01, 6.152596390471698E-01,
		   0.000000000000000E+00};
  real Qtrue[10] = { 6.882472016116853E-01, -3.647702873143482E-01,
		    8.595848852150008E-02, 3.686048903872426E-01,
		    -5.000000000000002E-01, 6.882472016116853E-01,
		    9.727207661715961E-02, 2.344322414222703E-02,
		     -5.160468465421397E-01, 5.000000000000002E-01};
  // only first 10 entries for Q

  np = 2;  nk = 3;
  ierr = xf_Error(xf_QRBulgeChase(shift, nk, np, alpha, beta, Q));
  xf_AssertEqual(ierr, xf_OK);

  xf_AssertRealVectorWithin(beta, beta0, 5, UTOL2);
  xf_AssertRealVectorWithin(alpha, alpha0, 5, UTOL2);
  for (k=0; k<10; k++) Q[k] = fabs(Q[k]);
  for (k=0; k<10; k++) Qtrue[k] = fabs(Qtrue[k]);
  xf_AssertRealVectorWithin(Q, Qtrue, 10, UTOL2);

  return xf_OK;
}

TEST_xf_QRBulgeChase2()
{
  int ierr, k, np, nk;
  real shift[2] = {0.0, 1.0};
  real alpha[8] = {2,2,2,2,2,2,2,2};
  real beta[8] = {17, 1,1,1,1,1,1,1}; // first entry is irrelevant
  real Q[64];

  np = 1;  nk = 7;
  ierr = xf_Error(xf_QRBulgeChase(shift, nk, np, alpha, beta, Q));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;
}

TEST_xf_SortIntPos()
{
  int ierr, i;
  int u[5] = {1, -3, 8, 2, 7};
  int u_true[5] = {-3, 1, 2, 7, 8};
  int pos[5];
  int pos_true[5] = {1, 0, 4, 2, 3};

  ierr = xf_Error(xf_SortIntPos(u, 5, pos, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
  
  xf_AssertIntVectorEqual(u, u_true, 5);
  xf_AssertIntVectorEqual(pos, pos_true, 5); 

  return xf_OK;
}


TEST_xf_SortReal()
{
  int ierr, i;
  real u[5] = {.2, -3.0, 7.2, 1.5, 7.1};
  real u_true[5] = {-3.0, .2, 1.5, 7.1, 7.2};
  int pos[5];
  int pos_true[5] = {1, 0, 4, 2, 3};

  ierr = xf_Error(xf_SortReal(u, 5, pos, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
  
  xf_AssertRealVectorWithin(u, u_true, 5, UTOL0);
  xf_AssertIntVectorEqual(pos, pos_true, 5); 

  return xf_OK;
}


TEST_xf_MergeSort()
{
  int ierr, i, j;
  int nlist = 3;
  real uList0[3][5] = {{.1,.3,.5,.8,0},
		       {.01,.2,.4,0,0},
		       {.5,.6,.65,.7,.9}};
  real **uList;
  int **iList;
  int iList_true[3][5] = {{1,3,5,10,0},
			  {0,2,4,0,0},
			  {6,7,8,9,11}};

  int N[3] = {4,3,5};

  // allocate data
  ierr = xf_Error(xf_Alloc2((void ***) &uList, nlist, 5, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc2((void ***) &iList, nlist, 5, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // fill in uList
  for (i=0; i<nlist; i++)
    for (j=0; j<N[i]; j++) uList[i][j] = uList0[i][j];


  ierr = xf_Error(xf_MergeSort(nlist, uList, N, iList));
  xf_AssertEqual(ierr, xf_OK);
  
  for (i=0; i<nlist; i++)
    xf_AssertIntVectorEqual(iList[i], iList_true[i], N[i]); 

  xf_Release2((void **) uList);
  xf_Release2((void **) iList);

  return xf_OK;
}


TEST_xf_QuadraticRoots()
{
  int  ierr;
  int  nroots;
  real roots[3];
  real coeff0[3] = {1.0, -2.0, 1.0};
  int  nroots0   = 1;
  real roots0[1] = {1.};
  real coeff1[3] = {10.0, -1.0, -2.0};
  int  nroots1   = 2;
  real roots1[2] = {-.4, .5};
  real coeff2[3] = {5.0, 1.0, 1.0};
  int  nroots2   = 0;
  

  ierr = xf_Error(xf_QuadraticRoots(coeff0, &nroots, roots));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nroots, nroots0);
  xf_AssertRealVectorWithin(roots, roots0, nroots0, UTOL0);


  ierr = xf_Error(xf_QuadraticRoots(coeff1, &nroots, roots));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nroots, nroots1);
  xf_AssertRealVectorWithin(roots, roots1, nroots1, UTOL0);


  ierr = xf_Error(xf_QuadraticRoots(coeff2, &nroots, roots));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nroots, nroots2);


  return xf_OK;
}


TEST_xf_RealNorm()
{
  real A[3] = {0.3, 0.4, 1.2};

  xf_AssertWithin(xf_RealNorm(A, 3), 1.3, UTOL0);

  return xf_OK;  
}

*/
