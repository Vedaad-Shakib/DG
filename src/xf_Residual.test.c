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
#include "xf_Data.h"
#include "xf_Solver.h"

// Functions for setting up and running a case within a unit test
#include "xf_UnitRun.c"



/** Several comprehensive ping tests **/


TEST_xf_CalculateResidual_Scalar_Motion()
{
  /* Residual pinging for scalar eqn with motion. */

  int ierr, i;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  char *ConvKeyValue[] = {"VelocityFcn", "Constant", "VelocityData", "1.0 0.0", "\0"};
  char *DiffKeyValue[] = {"Unused", "Unused", "\0"};
  xf_Vector *U;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_UnitMotionPlunge2));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;
  
  // Inactivate existing residual terms, add new ones
  for (i=0; i<All->EqnSet->ResTerms->nResTerm; i++)
    All->EqnSet->ResTerms->ResTerm[i].Active = xfe_False;

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermConv, ConvKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, DiffKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Set Time
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Time", "0.37"));
  if (ierr != xf_OK) return ierr;
  
  // Ping Residual
  ierr = xf_Error(xf_PingResidual(All, U, 1e-2, 1e-1));
  xf_AssertErrorOK(ierr);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_CalculateResidual_CNS_Motion()
{

  /* Residual pinging for compressible Navier-Stokes with mesh
     motion. */

  int ierr, i;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  real perturb = 1e-3;
  xf_Vector *U, *dU;
  xf_Vector *GCL, *dGCL;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAllCNS(&All, xfe_UnitMotionPlunge1));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  
  // Set Time
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Time", "0.1"));
  if (ierr != xf_OK) return ierr;

  // Set UseGCL
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "UseGCL", "True"));
  if (ierr != xf_OK) return ierr;

  // Initialize GCL vector
  ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U, -1, xfe_True, &GCL));
  if(ierr != xf_OK) return ierr;


  // Find a dU vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_False, xfe_True, 
                                       NULL, &dU, NULL));
  if (ierr != xf_OK) return ierr;
  
  // Find a dGCL vector
  ierr = xf_Error(xf_FindSimilarVector(All, GCL, "dGCL", xfe_False, xfe_True, 
                                       NULL, &dGCL, NULL));
  if (ierr != xf_OK) return ierr;


  // perturb the state
  ierr = xf_Error(xf_VectorRand(dU, 17));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMult(dU, perturb));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetVector(dU, xfe_Add, U));
  if (ierr != xf_OK) return ierr;
  

  // also create a small perturbation to the initial GCL vector
  ierr = xf_Error(xf_VectorRand(dGCL, 11));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMult(dGCL, perturb));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetVector(dGCL, xfe_Add, GCL));
  if (ierr != xf_OK) return ierr;

  // Ping Residual
  ierr = xf_Error(xf_PingResidual(All, U, 1e-3, 1e-1));
  xf_AssertErrorOK(ierr);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;
}


TEST_xf_CalculateResidual_CNS_SA()
{
  /* Residual pinging for compressible Navier-Stokes with SA
     turbulence model. */

  int ierr, i;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  xf_Vector *U, *dU;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAllCNS_SA(&All, xfe_UnitMotionNone));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // perturb the state so that it is not constant (uninteresting)

  // Find a dU vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_False, xfe_True, 
				       NULL, &dU, NULL));
  if (ierr != xf_OK) return ierr;
  
  // create a small perturbation of the IC: dU(0)
  ierr = xf_Error(xf_VectorRand(dU, 17)); // dU in [0,1]
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorAdd(dU, -0.5)); // dU in [-0.5, 0.5]
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMult(dU, 1.0e-1)); // dU in [-0.05, 0.05]
  if (ierr != xf_OK) return ierr;

  // perturb with dU(0)
  ierr = xf_Error(xf_SetVector(dU, xfe_Add, U));
  if (ierr != xf_OK) return ierr;

  // Set Time
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Time", "0.37"));
  if (ierr != xf_OK) return ierr;
  
  // Ping Residual
  ierr = xf_Error(xf_PingResidual(All, U, 1e-5, 4e-1));
  xf_AssertErrorOK(ierr);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_CalculateResidual_CNS_Hex_RegStab()
{

  /* Residual pinging for a compressible NS flow in a hex Q1 Q2 box.
     Regularity stabilization included. */

  int ierr, i;
  enum xfe_BasisType Basis = xfe_HexLagrange;
  char *StabKeyValue[] = {"Stabilization", "Resolution", "\0"};
  int Order = 1;
  xf_Vector *U;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitQ1Q2HexAllCNS(&All));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // add stabilization term
  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // First test using a constant stabilization indicator
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", "Const"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitchValue", "1.0"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_PingResidual(All, U, 1e-3, 1e-1));
  xf_AssertErrorOK(ierr);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_CalculateResidual_ScalarDiff()
{

  /* Residual pinging for a scalar convection-diffusion-reaction with
     jump stabilization */

  int ierr, i;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  char *ConvKeyValue[] = {"VelocityFcn", "Constant", "VelocityData", "1.0 0.0", "\0"};
  char *DiffKeyValue[] = {"Unused", "Unused", "\0"};
  char *SourceKeyValue[] = {"SourceFcn", "Arrhenius", "SourceData", "5.0 0.2 .15 .24", "\0"};
  //char *StabKeyValue[] = {"Stabilization", "Jump", "\0"};
  xf_Vector *U;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_UnitMotionNone));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // add regularity scalar key
  ierr = xf_AddKeyValue(&All->EqnSet->KeyValue, "RegularityScalar", "Scalar", xfe_True);
  xf_AssertEqual( ((ierr == xf_OK) || (ierr == xf_OVERWROTE)), 1);

  // Inactivate existing residual terms, add new ones
  for (i=0; i<All->EqnSet->ResTerms->nResTerm; i++)
    All->EqnSet->ResTerms->ResTerm[i].Active = xfe_False;

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermConv, ConvKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, DiffKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermSource, SourceKeyValue));
  xf_AssertEqual(ierr, xf_OK);

 /*  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue)); */
/*   xf_AssertEqual(ierr, xf_OK); */


  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // set U to a perturbed value
  ierr = xf_Error(xf_VectorRand(U, 17));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorMult(U, 1e-3));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorAdd(U, 1.0));
  xf_AssertEqual(ierr, xf_OK);
  
  ierr = xf_Error(xf_PingResidual(All, U, 1e-2, 1e-1));
  xf_AssertErrorOK(ierr);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_CalculateResidual_Euler_JumpStab()
{

  /* Residual pinging for a compressible Euler flow in a box.  Jump
     stabilization included. */

  int ierr, i;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  char *StabKeyValue[] = {"Stabilization", "Jump", "\0"};
  int Order = 2;
  xf_Vector *U;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAllEuler(&All, xfe_UnitMotionNone));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // add stabilization term
  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // First test using a constant stabilization indicator
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", "Const"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitchValue", "1.0"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_PingResidual(All, U, 1e-3, 1e-1));
  xf_AssertErrorOK(ierr);

  // Second, test using a square stabilization indicator
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", "Square"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_PingResidual(All, U, 1e-4, 2e-1)); // more nonlinear
  xf_AssertErrorOK(ierr);


  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_CalculateResidual_Euler_JumpStab_IP()
{

  /* Residual pinging for a compressible Euler flow in a box.  Jump
     stabilization included.  Interior Penalty (non-default) viscous
     discretization is tested here.*/

  int ierr, i;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  char *StabKeyValue[] = {"Stabilization", "Jump", "\0"};
  int Order = 2;
  xf_Vector *U;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ2TriangleAllEuler(&All));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // add stabilization term
  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // Use IP viscous discretization
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "DiffusionDiscretization", "IP"));
  if (ierr != xf_OK) return ierr;

  // First test using a constant stabilization indicator
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", "Const"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitchValue", "1.0"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_PingResidual(All, U, 1e-3, 2e-1));
  xf_AssertErrorOK(ierr);

  // Second, test using a square stabilization indicator
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", "Square"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_PingResidual(All, U, 1e-4, 2e-1)); // more nonlinear
  xf_AssertErrorOK(ierr);


  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_CalculateResidual_Euler_RegStab()
{

  /* Residual pinging for a compressible Euler flow in a box.
     Resolution stabilization included. */

  int ierr, i;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  char *StabKeyValue[] = {"Stabilization", "Resolution", "\0"};
  xf_Vector *U, *V;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAllEuler(&All, xfe_UnitMotionNone));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // add stabilization term
  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
  xf_AssertEqual(ierr, xf_OK);
  
  // set U to a perturbed value
  ierr = xf_Error(xf_DuplicateVector(Mesh, U, &V));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorRand(V, 17));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorMultSet(V, 1e-2, xfe_Add, U));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_DestroyVector(V, xfe_True));
  xf_AssertEqual(ierr, xf_OK);



  // First test using a constant stabilization indicator
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", "Const"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitchValue", "1.0"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_PingResidual(All, U, 1e-3, 0.25));  
  xf_AssertErrorOK(ierr);

  // Second, test using a Linear stabilization indicator
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", "Linear"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_PingResidual(All, U, 1e-5, 1e-1)); // more nonlinear
  xf_AssertErrorOK(ierr);


  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_CalculateResidual_CNS_JumpStab()
{

  /* Residual pinging for a compressible NS flow in a box.  Jump
     stabilization included.  Interior Penalty (non-default) viscous
     discretization is tested here.*/

  int ierr, i;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  char *StabKeyValue[] = {"Stabilization", "Resolution", "\0"};
  int Order = 2;
  xf_Vector *U;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ2TriangleAllCNS(&All));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // add stabilization term
  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // First test using a constant stabilization indicator
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", "Const"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitchValue", "1.0"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_PingResidual(All, U, 1e-3, 2e-1));
  xf_AssertErrorOK(ierr);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_ResidualOrderIndep()
{
  /* Residual order independence test for various equations */

  int ierr, i;
  int iEqn, nEqn = 4;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  int OrderIncrement = 1;
  char *StabKeyValue[] = {"Stabilization", "Resolution", "\0"};
  real ddp, dpH, dph;
  real vtol[4] = {1e-10, 1e-8, 1e-8, 1e-8};
  xf_Vector *UH, *dUH, *VH, *RH, *Rh;
  xf_SolverData *SolverData;
  xf_All *All;
  xf_Mesh *Mesh;

  for (iEqn=0; iEqn<nEqn; iEqn++){

    // Obtain appropriate All
    if (iEqn==0){
      // Scalar
      ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
    }
    else if (iEqn==1){
      // Compressible Navier-Stokes
      ierr = xf_Error(xf_UnitBoxQ1TriangleAllCNS(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
    }
    else if (iEqn==2){
      // Compressible Navier-Stokes with SA RANS model
      ierr = xf_Error(xf_UnitBoxQ1TriangleAllCNS_SA(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
    }
    else if (iEqn==3){
      // Compressible Navier-Stokes with SA RANS model ...
      ierr = xf_Error(xf_UnitBoxQ1TriangleAllCNS_SA(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
      // ... with a shock-capturing stabilization term
      ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
      xf_AssertEqual(ierr, xf_OK);
    }
  
    Mesh = All->Mesh;

    // Load dynamic library, register eqnset, initialize solution (minus -> variable order)
    ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &UH));
    xf_AssertEqual(ierr, xf_OK);

    // perturb the state so that it is not constant

    // Find a dUH vector
    ierr = xf_Error(xf_FindSimilarVector(All, UH, "dUH", xfe_False, xfe_True, 
                                         NULL, &dUH, NULL));
    if (ierr != xf_OK) return ierr;
  
    // create a small (not actually random) perturbation of the IC: dUH(0)
    ierr = xf_Error(xf_VectorRand(dUH, 17)); // dUH in [0,1]
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_VectorAdd(dUH, -0.5)); // dUH in [-0.5, 0.5]
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_VectorMult(dUH, 0.1)); // dUH in [-0.05, 0.05]
    if (ierr != xf_OK) return ierr;

    // perturb UH with dUH(0)
    ierr = xf_Error(xf_SetVector(dUH, xfe_Add, UH));
    if (ierr != xf_OK) return ierr;

    // Create a VH vector
    ierr = xf_Error(xf_FindSimilarVector(All, UH, "VH", xfe_False, xfe_True, 
                                         NULL, &VH, NULL));
    if (ierr != xf_OK) return ierr;
  
    // make VH a (not actually) random vector
    ierr = xf_Error(xf_VectorRand(VH, 17)); // VH in [0,1]
    if (ierr != xf_OK) return ierr;

    // Create a residual vector, RH
    ierr = xf_Error(xf_FindSimilarVector(All, UH, "RH", xfe_False, xfe_True, 
                                         NULL, &RH, NULL));
    if (ierr != xf_OK) return ierr;

    // create solver data for residual evaluation
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    xf_AssertEqual(ierr, xf_OK);

    // Calculate residual
    ierr = xf_Error(xf_CalculateResidual(All, UH, RH, NULL, SolverData));
    if (ierr == xf_NEED_LAPACK){
      xf_printf("Need LAPACK to complete this test; continuing without an error.\n");
      return xf_OK;
    }
    if (ierr != xf_OK) return ierr;

    // destroy solver data
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    xf_AssertEqual(ierr, xf_OK);
  
    // Compute inner product with an arbitrary test function
    ierr = xf_Error(xf_VectorDot(VH, RH, &dpH));
    if (ierr != xf_OK) return ierr;


    // inject UH to higher order
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, UH, NULL, 
                                                           xfe_BasisLast, OrderIncrement));
    if (ierr != xf_OK) return ierr;

    // inject VH to higher order
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, VH, NULL, 
                                                           xfe_BasisLast, OrderIncrement));
    if (ierr != xf_OK) return ierr;

  
    // Create a residual vector, Rh (on order-incremented space)
    ierr = xf_Error(xf_FindSimilarVector(All, UH, "Rh", xfe_False, xfe_True, 
                                         NULL, &Rh, NULL));
    if (ierr != xf_OK) return ierr;

    // create solver data for residual evaluation
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    xf_AssertEqual(ierr, xf_OK);

    // want to evaluate the residual on the original order
    SolverData->ResidualOrderIncrement = -OrderIncrement;

    // Calculate residual
    ierr = xf_Error(xf_CalculateResidual(All, UH, Rh, NULL, SolverData));
    if (ierr != xf_OK) return ierr;

    // destroy solver data
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    xf_AssertEqual(ierr, xf_OK);
  
    // Galerkin orthogonality should hold if the residual p-independence is working
    ierr = xf_Error(xf_VectorDot(VH, Rh, &dph));
    if (ierr != xf_OK) return ierr;
  
    ddp = fabs(dpH - dph);

    xf_printf("iEqn=%d: dpH = %.10E, dph = %.10E\n", iEqn, dpH, dph);
    xf_AssertWithin(ddp, 0.0, vtol[iEqn]);

    // Destroy All
    ierr = xf_Error(xf_DestroyAll(All));
    xf_AssertEqual(ierr, xf_OK);

  } // iEqn

  return xf_OK;  

}


