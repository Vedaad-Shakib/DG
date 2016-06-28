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
#include "xf_All.h"
#include "xf_Residual.h"
#include "xf_Solver.h"

// Functions for setting up and running a case within a unit test
#include "xf_UnitRun.c"


TEST_xf_GivensQRStep()
{
  int ierr;
  int nH[3] = {2, 3, 4};
  real **H;
  real g[3] = {3,13,6}; // rhs
  real F[2]; // Givens rotation: c,s
  real Ftrue[2] = { 0.8, -0.6};
  real Htrue[2] = { 5.0,  0.0};
  real gtrue[2] = {10.2,  8.6};
  
  // alloc H
  ierr = xf_Error(xf_VAlloc2( (void ***) &H, 3, nH, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  H[0][0] = 4.; H[0][1] = 3.;
  H[1][0] = 1.; H[1][1] = 1.; H[1][2] = 0.;
  H[2][0] = 2.; H[2][1] = 5.; H[2][2] = 3.; H[2][3] = 0.;

  ierr = xf_Error(xf_GivensQRStep(0, H, g, F));
  xf_AssertEqual(ierr, xf_OK);

  xf_AssertRealVectorWithin(F   , Ftrue, 2, UTOL1);
  xf_AssertRealVectorWithin(H[0], Htrue, 2, UTOL1);
  xf_AssertRealVectorWithin(g   , gtrue, 2, UTOL1);

  // release H
  xf_Release2((void **)  H);

  return xf_OK;  
}


TEST_xf_SolveUpper()
{
  int ierr;
  int nH[3] = {2, 3, 4};
  real **H;
  real g[3] = {3,13,6}; // rhs
  real y[3], y0[3] = {1,-3,-2}; // solution
  
  // alloc H
  ierr = xf_Error(xf_VAlloc2( (void ***) &H, 3, nH, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  H[0][0] = 4.; H[0][1] = 0.;
  H[1][0] = 1.; H[1][1] = 1.; H[1][2] = 0.;
  H[2][0] = 2.; H[2][1] = 5.; H[2][2] = 3.; H[2][3] = 0.;

  ierr = xf_Error(xf_SolveUpper(H, g, 3, y));
  xf_AssertEqual(ierr, xf_OK);

  xf_AssertRealVectorWithin(y, y0, 3, UTOL1);

  // release H
  xf_Release2((void **)  H);

  return xf_OK;  
}


TEST_xf_LinearIterCG()
{
  /* 
     Tests Conjugate-Gradient solver on 16 x 16 system:

     [2 1        ]     [- 1]
     [1 2 1      ]     [- 2]
     [  1 2 1    ] U + [- 3] = 0
     [     ...   ]     [-..]
     [      1 2 1]     [-16]
     
  */

  int ierr, i, k;
  int r = 16;
  int nelem[1] = {1};
  enum xfe_Bool done;
  enum xfe_LinearStatusType Status;
  xf_Vector *V, *W, *U, *R;
  xf_VectorSet *VS;
  xf_LinearSolverData *LinearSolverData;
  real *rV, *rW;
  real U0[16] = {0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8}; // true solution

  // create a vectorset, VS to be used for the CG work vectors
  ierr = xf_Error(xf_UnitrVectorSet(3, 1, nelem, &r, &VS));
  xf_AssertEqual(ierr, xf_OK);

  // Create R = linear residual vector
  ierr = xf_Error(xf_UnitrVector(1, nelem, &r, NULL, &R));
  xf_AssertEqual(ierr, xf_OK);
  rV = R->GenArray[0].rValue[0];
  for (i=0; i<r; i++) rV[i] = -((real) i) - 1.0;

  // Create U = solution vector
  ierr = xf_Error(xf_UnitrVector(1, nelem, &r, NULL,  &U));
  xf_AssertEqual(ierr, xf_OK);
  
  LinearSolverData = NULL;
  done = xfe_False;
  while (!done){
    
    ierr = xf_Error(xf_LinearIterCG(NULL, VS, 1e-12, xfe_VerbosityMedium, xfe_False, 
				    xfe_True, xfe_False, xfe_False,
				    R, U, &V, &W, &LinearSolverData, &Status));
    xf_AssertEqual(ierr, xf_OK);
    done = (Status == xfe_LinearConverged);
    if (!done){
      // W = A*V
      rV = V->GenArray[0].rValue[0];
      rW = W->GenArray[0].rValue[0];
      for (k=0; k<r; k++){
	rW[k] = 2.0*rV[k];
	if (k>0  ) rW[k] += rV[k-1];
	if (k<r-1) rW[k] += rV[k+1];
      }
      xf_AssertEqual((LinearSolverData->iIter < 100), xfe_True);
    }
  }

  rV = U->GenArray[0].rValue[0];
  //for (i=0; i<r; i++) xf_printf("  %.15E\n", rV[i]);
  xf_AssertRealVectorWithin(rV, U0, 16, UTOL2);

  ierr = xf_Error(xf_DestroyVectorSet(VS));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyVector(R, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyVector(U, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_ILUConsistency()
{
  /* Compares action of R_U to action of M and N for ILU preconditioner */
  int ierr, icase;
  int iTranspose;
  int dim, vdim[] = {1, 2, 2, 3};
  int Order, vOrder[] = {2, 2, 2, 1};
  enum xfe_Bool vTransposeFlag[2] = {xfe_False, xfe_True};
  enum xfe_Bool TransposeFlag;
  enum xfe_BasisType vBasis[] = {xfe_SegLagrange, xfe_TriLagrange, xfe_TriLagrange, xfe_HexLagrange};
  enum xfe_BasisType Basis;
  real rnorm;
  xf_Vector *U, *R, *X, *Y0, *Y1;
  xf_SolverData *SolverData;
  xf_JacobianMatrix *R_U = NULL;
  xf_All *All;

  // loop over cases
  for (icase=0; icase<4; icase++){

    Basis = vBasis[icase];
    Order = vOrder[icase];
    dim   = vdim[icase];
    
    // loop over False/True for TransposeFlag
    for (iTranspose=0; iTranspose<2; iTranspose++){
      TransposeFlag = vTransposeFlag[iTranspose];

      //xf_printf("icase=%d, TransposeFlag=%d\n", icase, iTranspose); fflush(stdout);

      // Get appropriate All structure Includes Param and EqnSet
      if (icase== 0)
	ierr = xf_Error(xf_UnitIntervalAllEuler(&All));
      else if (icase == 1)
	ierr = xf_Error(xf_UnitBoxQ1TriangleAllCNS(&All, xfe_UnitMotionNone));
      else if (icase == 2)
	ierr = xf_Error(xf_UnitSkewedQ1Tri5AllEuler(&All));
      else
	ierr = xf_Error(xf_UnitBoxQ1Hex2EG27All(&All));
      xf_AssertEqual(ierr, xf_OK);
    
      // Load dynamic library, register eqnset, initialize solution, test variable order
      ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &U));
      xf_AssertEqual(ierr, xf_OK);

      // locate Residual vector
      ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, xfe_True, 
					   NULL, &R, NULL));
      xf_AssertEqual(ierr, xf_OK);

      // locate Jacobian vector
      ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
					    xfe_True, NULL, &R_U, NULL));
      xf_AssertEqual(ierr, xf_OK);
  
      // create solver data for residual evaluation
      ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
      xf_AssertEqual(ierr, xf_OK);

      // calculate residual and Jacobian
      ierr = xf_Error(xf_CalculateResidual(All, U, R, R_U, SolverData));
      xf_AssertEqual(ierr, xf_OK);
  
      // Create a perturbation vector, X
      ierr = xf_Error(xf_FindSimilarVector(All, U, "XTestILU", xfe_False, xfe_True, 
					   NULL, &X, NULL));
      xf_AssertEqual(ierr, xf_OK);
      ierr = xf_Error(xf_VectorRand(X, 11));
      xf_AssertEqual(ierr, xf_OK);

      // Apply un-preconditioned R_U to X -> Y0
      ierr = xf_Error(xf_FindSimilarVector(All, U, "Y0TestILU", xfe_False, xfe_True, 
					   NULL, &Y0, NULL));
      xf_AssertEqual(ierr, xf_OK);
      ierr = xf_Error(xf_Jacobian_Mult(All, R_U, X, xfe_Set, TransposeFlag, SolverData, Y0));
      xf_AssertEqual(ierr, xf_OK);
  
      // precondition the Jacobian with ILU
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "ILUOrdering", "MDF"));
      xf_AssertEqual(ierr, xf_OK);
      ierr = xf_Error(xf_Jacobian_Precondition(All, R_U, xfe_PreconditionerILU0, xfe_False));
      xf_AssertEqual(ierr, xf_OK);

      // Apply preconditioned R_U (via M and N) to X -> Y1
      ierr = xf_Error(xf_FindSimilarVector(All, U, "Y1TestILU", xfe_False, xfe_True, 
					   NULL, &Y1, NULL));
      xf_AssertEqual(ierr, xf_OK);
      ierr = xf_Error(xf_Jacobian_Mult(All, R_U, X, xfe_Set, TransposeFlag, SolverData, Y1));
      xf_AssertEqual(ierr, xf_OK);

      // compare Y0 and Y1
      ierr = xf_Error(xf_SetVector(Y0, xfe_Sub, Y1));
      xf_AssertEqual(ierr, xf_OK);
      ierr = xf_Error(xf_VectorNorm(Y1, 2, &rnorm));
      xf_AssertEqual(ierr, xf_OK);
      xf_AssertWithin(rnorm, 0.0, UTOL3);

      // destroy solver data
      ierr = xf_Error(xf_DestroySolverData(SolverData));
      xf_AssertEqual(ierr, xf_OK);

      // Destroy All
      ierr = xf_Error(xf_DestroyAll(All));
      xf_AssertEqual(ierr, xf_OK);

    } // iTransposeFlag

  } // icase

  return xf_OK;  
}



TEST_xf_ILUvsLine()
{
  /* Compares ILU0 and line solver for a 1D case (should be equivalent
     in that both should be exact solves) */
  int ierr;
  int iTranspose;
  int Order = 2;
  enum xfe_Bool vTransposeFlag[2] = {xfe_False, xfe_True};
  enum xfe_Bool TransposeFlag;
  enum xfe_BasisType Basis = xfe_SegLagrange;
  real rnorm;
  xf_Vector *U, *R, *XLine, *XILU;
  xf_SolverData *SolverData;
  xf_JacobianMatrix *R_U = NULL;
  xf_All *All;

  // loop over False/True for TransposeFlag
  for (iTranspose=0; iTranspose<2; iTranspose++){
    TransposeFlag = vTransposeFlag[iTranspose];

    // Includes Param and EqnSet
    ierr = xf_Error(xf_UnitIntervalAllEuler(&All));
    xf_AssertEqual(ierr, xf_OK);

    // Load dynamic library, register eqnset, initialize solution, use variable order
    ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &U));
    xf_AssertEqual(ierr, xf_OK);

    // set default preconditioner as line
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Preconditioner", "LineJacobi"));
    xf_AssertEqual(ierr, xf_OK);

    // locate Residual vector
    ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, xfe_True, 
					 NULL, &R, NULL));
    xf_AssertEqual(ierr, xf_OK);

    // locate Jacobian vector
    ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
					  xfe_True, NULL, &R_U, NULL));
    xf_AssertEqual(ierr, xf_OK);
  
    // create solver data for residual evaluation
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    xf_AssertEqual(ierr, xf_OK);

    // calculate residual and Jacobian
    ierr = xf_Error(xf_CalculateResidual(All, U, R, R_U, SolverData));
    xf_AssertEqual(ierr, xf_OK);
  
    // Create a perturbation vector, XLine
    ierr = xf_Error(xf_FindSimilarVector(All, U, "XTestLine", xfe_False, xfe_True, 
					 NULL, &XLine, NULL));
    xf_AssertEqual(ierr, xf_OK);
    ierr = xf_Error(xf_VectorRand(XLine, 11));
    xf_AssertEqual(ierr, xf_OK);

    // Create a perturbation vector, XILU
    ierr = xf_Error(xf_FindSimilarVector(All, U, "XTestILU", xfe_False, xfe_True, 
					 NULL, &XILU, NULL));
    xf_AssertEqual(ierr, xf_OK);
    ierr = xf_Error(xf_VectorRand(XILU, 11));
    xf_AssertEqual(ierr, xf_OK);

  
    // precondition the Jacobian with Line
    ierr = xf_Error(xf_Jacobian_Precondition(All, R_U, xfe_PreconditionerLineJacobi, xfe_False));
    xf_AssertEqual(ierr, xf_OK);
    // Solve the system with Line
    ierr = xf_Error(xf_Jacobian_SolveM(All, R_U, XLine, xfe_Set, TransposeFlag, SolverData));
    xf_AssertEqual(ierr, xf_OK);

    // re-calculate the Jacobian
    ierr = xf_Error(xf_CalculateResidual(All, U, R, R_U, SolverData));
    xf_AssertEqual(ierr, xf_OK);
  
    // precondition the Jacobian with ILU
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "ILUOrdering", "Lex"));
    xf_AssertEqual(ierr, xf_OK);
    ierr = xf_Error(xf_Jacobian_Precondition(All, R_U, xfe_PreconditionerILU0, xfe_False));
    xf_AssertEqual(ierr, xf_OK);
    // Solve the system with ILU
    ierr = xf_Error(xf_Jacobian_SolveM(All, R_U, XILU, xfe_Set, TransposeFlag, SolverData));
    xf_AssertEqual(ierr, xf_OK);
 

    // compare XLine and XILU
    ierr = xf_Error(xf_SetVector(XLine, xfe_Sub, XILU));
    xf_AssertEqual(ierr, xf_OK);
    ierr = xf_Error(xf_VectorNorm(XILU, 2, &rnorm));
    xf_AssertEqual(ierr, xf_OK);
    xf_AssertWithin(rnorm, 0.0, UTOL4);

    // destroy solver data
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    xf_AssertEqual(ierr, xf_OK);

    // Destroy All
    ierr = xf_Error(xf_DestroyAll(All));
    xf_AssertEqual(ierr, xf_OK);


  } // iTransposeFlag

  return xf_OK;  
}

