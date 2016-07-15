/*******************************************************************************                              
 **                                                                                                           
 ** "math_unittest.cc": Unittests for the source file xf_Math.c
 ** Author: Vedaad Shakib
 **                                                                                                          
 *******************************************************************************/

#include <limits.h>
#include <math.h>
#include "gtest/gtest.h"
#include "UnitTest.h"

extern "C" {
#include "xf.h"
#include "xf_Memory.h"
#include "xf_AllStruct.h"
#include "xf_All.h"
#include "xf_Basis.h"
#include "xf_Data.h"
#include "xf_Math.h"
}

using namespace std;

TEST(xf_PowInt, All)
{
  real y;
  
  y = xf_PowInt(3.0, 4);
  EXPECT_DOUBLE_EQ(y, 81.0);
}

TEST(xf_MxV, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real u[2] = {2, 1};
  real v[3] = {-1, 0, 2};
  real v0[3] = {3, 10, 18};
  real v1[3] = {-1, 0, 2};
  real v2[3] = { 4,  10,  16};
  real v3[3] = {-4, -10, -16};

  xf_MxV(A, u, 3, 2, xfe_Add, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v0, 3, UTOL0));

  xf_MxV(A, u, 3, 2, xfe_Sub, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v1, 3, UTOL0));

  xf_MxV(A, u, 3, 2, xfe_Set, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v2, 3, UTOL0));

  xf_MxV(A, u, 3, 2, xfe_Neg, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v3, 3, UTOL0));

    
}

TEST(xf_cMxV_Add, All)
{
  real c = 2.0;
  real A[6] = {1, 2, 3, 4, 5, 6};
  real u[2] = {2, 1};
  real v[3] = {-1, 0, 2};
  real v_true[3] = {7, 20, 34};

  xf_cMxV_Add(c, A, u, 3, 2, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v_true, 3, UTOL0));

    
}


TEST(xf_MTxV, All)
{
  real A[6] = {1, 3, 5, 2, 4, 6};
  real u[2] = {2, 1};
  real v[3] = {-1, 0, 2};
  real v0[3] = {3, 10, 18};
  real v1[3] = {-1, 0, 2};

  xf_MTxV_Add(A, u, 3, 2, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v0, 3, UTOL0));

  xf_MTxV_Sub(A, u, 3, 2, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v1, 3, UTOL0));

    
}

TEST(xf_MxM, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 2, 1, 3, 1, 2};
  real C[4] = {2, 3, -1, -10};
  real C0[4] = {-4, -11, -16, -45};
  real C1[4] = {2, 3, -1, -10};
  real C2[4] = {6, 14, 15, 35};
  real C3[4] = {-6, -14, -15, -35};

  xf_MxM(A, B, 2, 3, 2, xfe_Sub, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C0, 4, UTOL0));

  xf_MxM(A, B, 2, 3, 2, xfe_Add, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C1, 4, UTOL0));

  xf_MxM(A, B, 2, 3, 2, xfe_Set, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C2, 4, UTOL0));

  xf_MxM(A, B, 2, 3, 2, xfe_Neg, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C3, 4, UTOL0));

    
}

TEST(xf_nMxM, All)
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
  EXPECT_TRUE(AssertRealVectorWithin(C, C0, 8, UTOL0));

  xf_nMxM_Add(2, A, B, 2, 3, 2, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C1, 8, UTOL0));

  xf_nMxM_Set(2, Ap, Bp, 2, 2, 1, Cp);
  EXPECT_TRUE(AssertRealVectorWithin(Cp, Cp0, 4, UTOL0));


    
}

TEST(xf_cMxM_Add, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 2, 1, 3, 1, 2};
  real C[4] = {2, 3, -1, -10};
  real C0[4] = {-1, -4, -8.5, -27.5};

  xf_cMxM_Add(-0.5, A, B, 2, 3, 2, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C0, 4, UTOL0));

    
}


TEST(xf_MTxM, All)
{
  real A[6] = {1, 4, 2, 5, 3, 6};
  real B[6] = {1, 2, 2, 1, 3, 1};
  real C[4] = {1, 1, 0, 0};
  real C0[4] = {15, 8, 32, 19};
  real C1[4] = {1, 1, 0, 0};
  real C2[4] = {14, 7, 32, 19};
  real C3[4] = {-14, -7, -32, -19};

  xf_MTxM(A, B, 2, 3, 2, xfe_Add, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C0, 4, UTOL0));

  xf_MTxM(A, B, 2, 3, 2, xfe_Sub, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C1, 4, UTOL0));

  xf_MTxM(A, B, 2, 3, 2, xfe_Set, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C2, 4, UTOL0));

  xf_MTxM(A, B, 2, 3, 2, xfe_Neg, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C3, 4, UTOL0));

    
}


TEST(xf_MTxwM, All)
{
  real A[6] = {1, 4, 2, 5, 3, 6};
  real w[3] = {3, 4, 5};
  real B[6] = {1, 2, 2, 1, 3, 1};
  real C[4];
  real C0[4] = {64, 29, 142, 74};
  real C1[4] = {128, 58, 284, 148};

  xf_MTxwM_Set(A, w, B, 2, 3, 2, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C0, 4, UTOL0));

  xf_MTxwM_Add(A, w, B, 2, 3, 2, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C1, 4, UTOL0));

    
}


TEST(xf_MxMT_Set, All)
{
  real A[6] = {1, 4, 2, 5, 3, 6};
  real B[6] = {1, 2, 2, 1, 3, 1};
  real C[4];
  real C0[4] = {13, 15, 23, 20};

  xf_MxMT_Set(A, B, 2, 3, 2, C);
  EXPECT_TRUE(AssertRealVectorWithin(C, C0, 4, UTOL0));

    
}



TEST(xf_ColMult, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real A_true[6] = {1, 2, 9, 12, 25, 30};

  xf_ColMult(A, v, 3, 2, 2);
  EXPECT_TRUE(AssertRealVectorWithin(A, A_true, 6, UTOL0));

    
}

TEST(xf_ColMult_Set, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {1, 2, 9, 12, 25, 30}, B[6];

  xf_ColMult_Set(A, v, 3, 2, 2, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

    
}


TEST(xf_ColMult_Add, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 1, 0, 1, 0, 1};
  real B_true[6] = {2, 3, 9, 13, 25, 31};

  xf_ColMult_Add(A, v, 3, 2, 2, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

    
}

TEST(xf_ColMult_Sub, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 1, 0, 1, 0, 1};
  real B_true[6] = {0, -1, -9, -11, -25, -29};

  xf_ColMult_Sub(A, v, 3, 2, 2, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

    
}


TEST(xf_ColcMult, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real A_true[6] = {1.5, 3.0, 13.5, 18, 37.5, 45};

  xf_ColcMult(A, v, 3, 2, 2, 1.5);
  EXPECT_TRUE(AssertRealVectorWithin(A, A_true, 6, UTOL0));

    
}

TEST(xf_ColcMult_Set, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {1.5, 3.0, 13.5, 18, 37.5, 45}, B[6];

  xf_ColcMult_Set(A, v, 3, 2, 2, 1.5, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

    
}

TEST(xf_ColcMult_Add, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {0.5, 2.0, 12.5, 17, 36.5, 44};
  real B[6] = {-1, -1, -1, -1, -1, -1};

  xf_ColcMult_Add(A, v, 3, 2, 2, 1.5, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

    
}

TEST(xf_2ColMult_Set, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real w[6] = {2, 1, 4, 3, 0, 5};
  real B_true[6] = {2, 4, 36, 48, 0, 0}, B[6];

  xf_2ColMult_Set(A, v, w, 3, 2, 2, 2, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

    
}

TEST(xf_2ColcMult_Set, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real w[6] = {2, 1, 4, 3, 0, 5};
  real B_true[6] = {1, 2, 18, 24, 0, 0}, B[6];

  xf_2ColcMult_Set(A, v, w, 3, 2, 2, 2, 0.5, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

    
}

TEST(xf_2ColcMult_Add, All)
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real v[6] = {1, 2, 3, 4, 5, 6};
  real w[6] = {2, 1, 4, 3, 0, 5};
  real B_true[6] = {0, 0, 15, 20, -5, -6};
  real B[6] = {-1, -2, -3, -4, -5, -6};

  xf_2ColcMult_Add(A, v, w, 3, 2, 2, 2, 0.5, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

    
}

TEST(xf_MatDetInv, All)
{
  int ierr;
  real A[4] = {1,2,3,4};
  real iA[9], detA;
  real iA_true[4] = {-2,1,1.5,-0.5};
  real B[9] = {1,0,-1,0,1,1,1,2,0};
  real iB[9], detB;
  real iB_true[9] = {2,2,-1,-1,-1,1,1,2,-1};

  ierr = xf_Error(xf_MatDetInv(A, 2, &detA, iA));
  EXPECT_EQ(ierr, 0);
  EXPECT_DOUBLE_EQ(detA, -2.0);
  EXPECT_TRUE(AssertRealVectorWithin(iA, iA_true, 4, UTOL0));

  ierr = xf_Error(xf_MatDetInv(B, 3, &detB, iB));
  EXPECT_EQ(ierr, 0);
  EXPECT_DOUBLE_EQ(detB, -1.0);
  EXPECT_TRUE(AssertRealVectorWithin(iB, iB_true, 9, UTOL0));

    
}

TEST(xf_DotProduct, All)
{
  real a[3] = {1,2,3};
  real b[3] = {4,5,6};
  real c;
  
  xf_DotProduct(a, b, 3, &c);
  EXPECT_DOUBLE_EQ(c, 32.0);

    
}


TEST(xf_CrossProduct, All)
{
  real a[3] = {1,2,3};
  real b[3] = {4,5,6};
  real c_true[3] = {-3,6,-3}, c[3];
  
  xf_CrossProduct(a, b, c);
  EXPECT_TRUE(AssertRealVectorWithin(c, c_true, 3, UTOL0));

    
}

TEST(xf_ComputeSolvePLUT, All)
{
  int ierr, P[3];
  real A[9] = {1, 2, -1, 1, -1, -2, 0, 2, 1};
  real b[3] = {1, 2, 3};
  real u0[3] = {20.0, -4.0, 11.0}, u[3];
  real u1[3] = {-7.0, 8.0, 12.0};

  ierr = xf_Error(xf_ComputePLU(A, 3, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_SolvePLU(A, P, b, 3, u, NULL));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(u, u0, 3, UTOL1));

  ierr = xf_Error(xf_SolvePLUT(A, P, b, 3, u, NULL));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(u, u1, 3, UTOL1));

  
}


TEST(xf_SolvePLU_Matrix, All)
{
  int ierr, P[2];
  real A[4] = {1, 2, 3, 4};
  real B[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {2, 1, 0, -.5, .5, 1.5};

  ierr = xf_Error(xf_ComputePLU(A, 2, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_SolvePLU_Matrix(A, P, 2, 3, B));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

  
}

TEST(xf_SolvePLU_MatrixR, All)
{
  int ierr, P[2];
  real A[4] = {1, 2, 3, 4};
  real B[6] = {1, 2, 3, 4, 5, 6};
  real B_true[6] = {1, 0, 0, 1, -1, 2};

  ierr = xf_Error(xf_ComputePLU(A, 2, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_SolvePLU_MatrixR(A, P, 2, 3, B));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 6, UTOL0));

  
}

TEST(xf_ComputeSolveBlockPLUT, All)
{
  int ierr, k, P[4];
  real A[16] =  {1,2,4,3, 3,4,2,1, 1,3,3,2, 2,4,1,4};
  real b[4] = {1, 6, 3, 4};
  real c[4] = {9, 3, 6, 1};
  real u0[4] = {0.85, 1.35, -0.65, -0.15}, u[4];
  real u1[4] = {2.45, 2.2, -3.0, 0.25};

  ierr = xf_Error(xf_ComputeBlockPLU(A, 2, 2, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_SolveBlockPLU(A, 2, 2, P, b, xfe_Set, u));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(u, u0, 4, UTOL0));

  ierr = xf_Error(xf_SolveBlockPLUT(A, 2, 2, P, c, xfe_Set, u));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(u, u1, 4, UTOL0));

  
}

TEST(xf_ComputeSolveBlockPLU2, All)
{
  int ierr, k, P[4];
  real A[16] =  {1,2,3,4, 4,3,2,1, 1,3,2,4, 3,2,1,4};
  real b[4] = {1, 6, 3, 4};
  real u_true[4] = {-0.85, -1.35, 0.65, 0.15}, u[4];

  ierr = xf_Error(xf_ComputeBlockPLU(A, 1, 4, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_SolveBlockPLU(A, 1, 4, P, b, xfe_Neg, u));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(u, u_true, 4, UTOL0));

  
}


TEST(xf_ComputeSolveBlockPLU_Matrix, All)
{
  int ierr, k, P[4];
  real A[16] =  {1,2,4,3, 3,4,2,1, 1,3,3,2, 2,4,1,4};
  real B[8] = {1, 2, 6, 12, 3, 6, 4, 8};
  real U_true[8] = {0.85, 1.7, 1.35, 2.7, -0.65, -1.3, -0.15, -.3}, U[8];
  real iA[16], I[16];

  ierr = xf_Error(xf_ComputeBlockPLU(A, 2, 2, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_SolveBlockPLU_Matrix(A, 2, 2, P, B, 1, xfe_Set, U));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(U, U_true, 8, UTOL0));

  ierr = xf_Error(xf_SolveBlockPLU_Matrix(A, 2, 2, P, B, 1, xfe_Sub, U));
  ierr = xf_Error(xf_SolveBlockPLU_Matrix(A, 2, 2, P, B, 1, xfe_Add, U));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(U, U_true, 8, UTOL0));

  // BlockIdentity is tested here
  xf_BlockIdentity(2, 2, I);
  ierr = xf_Error(xf_SolveBlockPLU_Matrix(A, 2, 2, P, I, 2, xfe_Set, iA));
  EXPECT_EQ(ierr, 0);
  xf_BlockMxBlockM(iA, 2, 2, 2, B, 1, xfe_Set, U);
  EXPECT_TRUE(AssertRealVectorWithin(U, U_true, 8, UTOL0));

  
}

TEST(xf_PLUMxV_Set, All)
{
  int ierr, P[2];
  real A[4] = {1, 2, 3, 4};
  real u[2] = {1, 2};
  real b_true[2] = {5, 11}, b[2];

  ierr = xf_Error(xf_ComputePLU(A, 2, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_PLUMxV_Set(A, P, u, 2, b, NULL));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(b, b_true, 2, UTOL0));

  
}

TEST(xf_PLUMTxV_Set, All)
{
  int ierr, P[3];
  real A[9] = {1, 2, -1, 1, -1, -2, 0, 2, 1};
  real u[3] = {1, 2, 3};
  real b0[3] = {3.0, 6.0, -2.0}, b[3];

  ierr = xf_Error(xf_ComputePLU(A, 3, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_PLUMTxV_Set(A, P, u, 3, b, NULL));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(b, b0, 3, UTOL0));

  
}

TEST(xf_BlockPLUMTxV, All)
{
  int ierr, k, P[6];
  real A[36] = {1,2,6,5, 3,4,4,3, 5,6,2,1, 
		1,3,2,5, 2,4,3,1, 5,6,4,6, 
		4,3,2,5, 2,5,3,2, 6,1,6,1};
  real u[6] = {1, 2, 3, 4, 5, 6};
  real b0[6] = {91, 56, 90, 81, 72, 65}, b[6];
  real b1[6] = {56, 86, 57, 63, 106, 61};

  ierr = xf_Error(xf_ComputeBlockPLU(A, 3, 2, P));
  EXPECT_EQ(ierr, 0);

  ierr = xf_Error(xf_BlockPLUMxV(A, 3, 2, P, u, xfe_Set, b));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(b, b0, 6, UTOL2));

  ierr = xf_Error(xf_BlockPLUMTxV(A, 3, 2, P, u, xfe_Set, b));
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(b, b1, 6, UTOL2));

  
}

TEST(xf_BlockMxV, All)
{
  real A[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  real u[4] = {2, 1, -1, -2};
  real v[2] = {0, 0};
  real z[2] = {0, 0};
  real v_true[2] = {-13, -13};

  xf_BlockMxV(A, 1, 2, 2, u, xfe_Set, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v_true, 2, UTOL0));

  xf_BlockMxV(A, 1, 2, 2, u, xfe_Sub, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, z, 2, UTOL0));

  xf_BlockMxV(A, 1, 2, 2, u, xfe_Add, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v_true, 2, UTOL0));

  
}

TEST(xf_BlockMTxV, All)
{
  real A[8] = {1, 3, 2, 4, 5, 7, 6, 8};
  real u[4] = {2, 1, -1, -2};
  real v[2] = {0, 0};
  real z[2] = {0, 0};
  real v_true[2] = {-13, -13};

  xf_BlockMTxV(A, 1, 2, 2, u, xfe_Set, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v_true, 2, UTOL0));

  xf_BlockMTxV(A, 1, 2, 2, u, xfe_Sub, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, z, 2, UTOL0));

  xf_BlockMTxV(A, 1, 2, 2, u, xfe_Add, v);
  EXPECT_TRUE(AssertRealVectorWithin(v, v_true, 2, UTOL0));

  
}


TEST(xf_BlockMxM, All)
{
  real A[32] = {1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7, 
		5,6,7,8, 6,7,8,9, 7,8,9,10, 8,9,10,11};
  real U[8] = {1,0, 0,1, 1,0, 0,1};
  real B_true[16] = {4,6,8,10, 6,8,10,12, 12,14,16,18, 14,16,18,20};
  real B[16] = {0};
  real Z[16] = {0};

  xf_BlockMxM(A, 2, 2, 4, U, 2, xfe_Set, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 16, UTOL0));

  xf_BlockMxM(A, 2, 2, 4, U, 2, xfe_Sub, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, Z, 16, UTOL0));

  xf_BlockMxM(A, 2, 2, 4, U, 2, xfe_Neg, B);
  xf_BlockMxM(A, 2, 2, 4, U, 2, xfe_Add, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, Z, 16, UTOL0));

  
}


TEST(xf_MxBlockM, All)
{
  real A[32] = {1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7, 
		5,6,7,8, 6,7,8,9, 7,8,9,10, 8,9,10,11};
  real U[2] = {1, 1};
  real B_true[16] = {6,8,10,12, 8,10,12,14, 10,12,14,16, 12,14,16,18};
  real B[16] = {0};
  real Z[16] = {0};

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Set, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 16, UTOL0));

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Sub, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, Z, 16, UTOL0));

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Neg, B);
  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Add, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, Z, 16, UTOL0));

  
}


TEST(xf_MTxBlockM, All)
{
  real A[32] = {1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7, 
		5,6,7,8, 6,7,8,9, 7,8,9,10, 8,9,10,11};
  real U[2] = {1, 1};
  real B_true[16] = {6,8,10,12, 8,10,12,14, 10,12,14,16, 12,14,16,18};
  real B[16] = {0};
  real Z[16] = {0};

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Set, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 16, UTOL0));

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Sub, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, Z, 16, UTOL0));

  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Neg, B);
  xf_MxBlockM(U, 1, 2, A, 2, 4, xfe_Add, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, Z, 16, UTOL0));

  
}

TEST(xf_BlockMxBlockM, All)
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
  EXPECT_TRUE(AssertRealVectorWithin(B, B_true, 24, UTOL0));

  xf_BlockMxBlockM(A, 2, 2, 4, U, 3, xfe_Sub, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, Z, 24, UTOL0));

  xf_BlockMxBlockM(A, 2, 2, 4, U, 3, xfe_Neg, B);
  xf_BlockMxBlockM(A, 2, 2, 4, U, 3, xfe_Add, B);
  EXPECT_TRUE(AssertRealVectorWithin(B, Z, 24, UTOL0));

  
}

TEST(xf_OutProd_Add, All)
{
  real u[3] = {1, 2};
  real v[3] = {4, -5};
  real A[4] = {-1, -2, -3, -4};
  real A_true[4] = {3, -7, 5, -14};

  xf_OutProd_Add(u, v, 2, 2, A);
  EXPECT_TRUE(AssertRealVectorWithin(A, A_true, 4, UTOL0));

    
}

TEST(xf_OutProd_Sub, All)
{
  real u[2] = {1, 2};
  real v[2] = {4, -5};
  real A[4] = {-1, -2, -3, -4};
  real A_true[4] = {-5, 3, -11, 6};

  xf_OutProd_Sub(u, v, 2, 2, A);
  EXPECT_TRUE(AssertRealVectorWithin(A, A_true, 4, UTOL0));

    
}

TEST(xf_BlockOutProd_Add, All)
{
  real u[2] = {-2, 3};
  real v[4] = {1, 2, 3, 4};
  real A[8] = {0, 0, -1, 0, 0, 0, 0, 1};
  real A_true[8] = {-2, -4, 2, 6, -6, -8, 9, 13};

  xf_BlockOutProd_Add(u, v, 1, 2, 2, A);
  EXPECT_TRUE(AssertRealVectorWithin(A, A_true, 8, UTOL0));

    
}

TEST(xf_BlockOutProd_Sub, All)
{
  real u[2] = {-2, 3};
  real v[4] = {1, 2, 3, 4};
  real A[8] = {0, 0, 0, 0, 0, 1, 0, 1};
  real A_true[8] = {2, 4, -3, -6, 6, 9, -9, -11};

  xf_BlockOutProd_Sub(u, v, 1, 2, 2, A);
  EXPECT_TRUE(AssertRealVectorWithin(A, A_true, 8, UTOL0));

    
}


TEST(xf_QRFactorHouseholder, All)
{
  int ierr;
  real A[6] = {2.4, 1.8, 3.2, 2.4, 0.0, 2.0};
  real Q_true[9] = {-.6, 0, -.8, -.8, 0, 0.6, 0, 1.0, 0};
  real R_true[6] = {-4.0, -3.0, 0, 2.0, 0, 0};
  real Q[9], R[6];
  
  ierr = xf_QRFactorHouseholder(A, 3, 2, Q, R);
  EXPECT_EQ(ierr, 0);

  EXPECT_TRUE(AssertRealVectorWithin(Q, Q_true, 9, UTOL0));
  EXPECT_TRUE(AssertRealVectorWithin(R, R_true, 6, UTOL0));

 
}


TEST(xf_EigSymmetricQR, All)
{
  int ierr;
  real tol = 1e-8;
  real A[9] = {1.64, -0.48, 0, -0.48, 1.36, 0, 0, 0, 3.0};
  real EG[3], EG_true[3] = {3,2,1};
  real EV[9], EV_true[9] = {0,0,-1, -0.8,0.6,0, -0.6,-0.8,0};
  
  ierr = xf_EigSymmetricQR(A, 3, tol, 100, EG, EV);
  EXPECT_EQ(ierr, 0);

  EXPECT_TRUE(AssertRealVectorWithin(EG, EG_true, 3, 10.*tol));
  EXPECT_TRUE(AssertRealVectorWithin(EV, EV_true, 9, 10.*tol));

 
}

TEST(xf_SVDGolubReinsch, All)
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
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(W, W_true, 2, UTOL1));
  EXPECT_TRUE(AssertRealVectorWithin(V, V_true, 4, UTOL1));
  EXPECT_TRUE(AssertRealVectorWithin(A, U_true, 6, UTOL1));

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
  EXPECT_EQ(ierr, 0);
  EXPECT_TRUE(AssertRealVectorWithin(W, W_true, 5, UTOL1));
  
  for (i=0; i<25; i++) B[i] *= W[i%5]; // B = B*diag(W)
  xf_MxMT_Set(B, V, 5, 5, 5, C); // C = B*V^T
  EXPECT_TRUE(AssertRealVectorWithin(C, B_true, 25, UTOL1));
       
  
}


TEST(xf_QRBulgeChase, All)
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
  EXPECT_EQ(ierr, xf_OK);

  EXPECT_TRUE(AssertRealVectorWithin(beta, beta0, 5, UTOL2));
  EXPECT_TRUE(AssertRealVectorWithin(alpha, alpha0, 5, UTOL2));
  for (k=0; k<10; k++) Q[k] = fabs(Q[k]);
  for (k=0; k<10; k++) Qtrue[k] = fabs(Qtrue[k]);
  EXPECT_TRUE(AssertRealVectorWithin(Q, Qtrue, 10, UTOL2));

  
}

TEST(xf_QRBulgeChase2, All)
{
  int ierr, k, np, nk;
  real shift[2] = {0.0, 1.0};
  real alpha[8] = {2,2,2,2,2,2,2,2};
  real beta[8] = {17, 1,1,1,1,1,1,1}; // first entry is irrelevant
  real Q[64];

  np = 1;  nk = 7;
  ierr = xf_Error(xf_QRBulgeChase(shift, nk, np, alpha, beta, Q));
  EXPECT_EQ(ierr, xf_OK);

  
}

TEST(xf_SortIntPos, All)
{
  int ierr, i;
  int u[5] = {1, -3, 8, 2, 7};
  int u_true[5] = {-3, 1, 2, 7, 8};
  int pos[5];
  int pos_true[5] = {1, 0, 4, 2, 3};

  ierr = xf_Error(xf_SortIntPos(u, 5, pos, xfe_True));
  EXPECT_EQ(ierr, xf_OK);
  
  EXPECT_TRUE(AssertIntVectorEqual(u, u_true, 5));
  EXPECT_TRUE(AssertIntVectorEqual(pos, pos_true, 5)); 

  
}


TEST(xf_SortReal, All)
{
  int ierr, i;
  real u[5] = {.2, -3.0, 7.2, 1.5, 7.1};
  real u_true[5] = {-3.0, .2, 1.5, 7.1, 7.2};
  int pos[5];
  int pos_true[5] = {1, 0, 4, 2, 3};

  ierr = xf_Error(xf_SortReal(u, 5, pos, xfe_True));
  EXPECT_EQ(ierr, xf_OK);
  
  EXPECT_TRUE(AssertRealVectorWithin(u, u_true, 5, UTOL0));
  EXPECT_TRUE(AssertIntVectorEqual(pos, pos_true, 5)); 

  
}


TEST(xf_MergeSort, All)
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
  EXPECT_EQ(ierr, 0);
  ierr = xf_Error(xf_Alloc2((void ***) &iList, nlist, 5, sizeof(int)));
  EXPECT_EQ(ierr, 0);
  
  // fill in uList
  for (i=0; i<nlist; i++)
    for (j=0; j<N[i]; j++) uList[i][j] = uList0[i][j];


  ierr = xf_Error(xf_MergeSort(nlist, uList, N, iList));
  EXPECT_EQ(ierr, xf_OK);
  
  for (i=0; i<nlist; i++)
    EXPECT_TRUE(AssertIntVectorEqual(iList[i], iList_true[i], N[i])); 

  xf_Release2((void **) uList);
  xf_Release2((void **) iList);

  
}


TEST(xf_QuadraticRoots, All)
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
  EXPECT_EQ(ierr, xf_OK);
  EXPECT_EQ(nroots, nroots0);
  EXPECT_TRUE(AssertRealVectorWithin(roots, roots0, nroots0, UTOL0));


  ierr = xf_Error(xf_QuadraticRoots(coeff1, &nroots, roots));
  EXPECT_EQ(ierr, xf_OK);
  EXPECT_EQ(nroots, nroots1);
  EXPECT_TRUE(AssertRealVectorWithin(roots, roots1, nroots1, UTOL0));


  ierr = xf_Error(xf_QuadraticRoots(coeff2, &nroots, roots));
  EXPECT_EQ(ierr, xf_OK);
  EXPECT_EQ(nroots, nroots2);


  
}


TEST(xf_RealNorm, All)
{
  real A[3] = {0.3, 0.4, 1.2};

  EXPECT_TRUE(AssertWithin(xf_RealNorm(A, 3), 1.3, UTOL0));
}

