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
#include "xf_Solver.h"

// Functions for setting up and running a case within a unit test
#include "xf_UnitRun.c"



TEST_xf_CalculateOutput_CNS_Motion()
{
  /*
    Tests calculation and linearization of outputs for CNS with motion
  */

  int ierr;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  int iOut, nOut = 5; // number of outputs to test
  const char *OutputNames[] = {"Point", "BoundaryInt", "BoundaryIntFS", 
                               "BoundaryIntNoFlux", "BoundaryIntTemp"}; 
  const char *OutputName;
  real perturb = 1e-2;
  xf_Vector *U, *dU;
  xf_Vector *GCL, *dGCL;
  xf_All *All;

  // get All with CNS eqnset
  ierr = xf_Error(xf_UnitBoxQ1TriangleAllCNS(&All, xfe_UnitMotionPlunge1));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &U)); // negative -> variable order
  xf_AssertEqual(ierr, xf_OK);

  // set verbosity
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", "High"));
  if (ierr != xf_OK) return ierr;

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

  // loop over outputs to test
  for (iOut = 0; iOut < nOut; iOut++){
    
    OutputName = OutputNames[iOut];

     // Create a J_GCL vector  for this output  
    ierr = xf_Error(xf_InitMeshMotionGCLLinearization(All, U, OutputName, -1, xfe_True, NULL));
    if (ierr != xf_OK) return ierr;

    // Ping output with respect to GCL
    xf_printf("Pinging Output = %s with respect to GCL\n", OutputName);
    ierr = xf_Error(xf_PingOutputGCL(All, OutputName, U, GCL));
    if (ierr != xf_OK) xf_printf(" Ping w.r.t GCL; Output = %s\n", OutputName);
    xf_AssertEqual(ierr, xf_OK);

    // Ping output with respect to state
    xf_printf("Pinging Output = %s with respect to U\n", OutputName);
    ierr = xf_Error(xf_PingOutput(All, OutputName, U));
    if (ierr != xf_OK) xf_printf(" Ping w.r.t U; Output = %s\n", OutputName);
    xf_AssertEqual(ierr, xf_OK);
    
  } // iOut

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;
}


TEST_xf_CalculateOutput_BoundaryScalar()
{
  /*
    Tests calculation and linearization of a boundary scalar integral output
  */

  int ierr;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  real J;
  xf_Vector *U;
  xf_All *All;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ2TriangleAllCNS(&All));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Solve system
  ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Calculate pressure integral output
  ierr = xf_Error(xf_CalculateOutput(All, "PressureIntegral", U, &J, NULL, xfe_Set));
  xf_AssertEqual(ierr, xf_OK);

  // Ping output
  ierr = xf_Error(xf_PingOutput(All, "PressureIntegral", U));
  xf_AssertEqual(ierr, xf_OK);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}




TEST_xf_CalculateOutput()
{
  /*
    Tests calculation and linearization of outputs
  */

  int ierr;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  real J;
  xf_Vector *U;
  xf_All *All;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Solve system
  ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Calculate heat flux out of right side
  ierr = xf_Error(xf_CalculateOutput(All, "HeatFlow", U, &J, NULL, xfe_Set));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(J, 0.1, UTOL2);

  // Ping output
  ierr = xf_Error(xf_PingOutput(All, "HeatFlow", U));
  xf_AssertEqual(ierr, xf_OK);

  
  // Calculate heat flux moment out of right side
  ierr = xf_Error(xf_CalculateOutput(All, "HeatFlowMoment", U, &J, NULL, xfe_Set));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertWithin(J, 0.1, UTOL2);
    
  // Ping output
  ierr = xf_Error(xf_PingOutput(All, "HeatFlowMoment", U));
  xf_AssertEqual(ierr, xf_OK);


  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}




TEST_xf_IncrementUnsteadyOutputs()
{
  /*
    Tests calculation and linearization of unsteady outputs
  */

  int ierr;
  int Order = 2;
  int iOrder, jOrder, OrderTime;
  int iScheme, nScheme = 1; // number of schemes to test
  int iOutput, nOutput = 3;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  const char *OutputNames[] = {"HeatFlowUnsteadyPoint", 
			       "HeatFlowUnsteadyInt",
			       "HeatFlowUnsteadySquare"}; 
  const char *OutputName;
  char StateName[xf_MAXSTRLEN];
  enum xfe_TimeSchemeType TimeScheme;
  enum xfe_TimeSchemeType TimeSchemes[] = {xfe_TimeSchemeDG1,
					   xfe_TimeSchemeDG2}; 
  real Times[][2] = {{0.6, 0.2}, // Time, TimeStep pairs
		     {0.6, 0.15},
		     {0.6, 0.12}};
  real tolvec[] = {UTOL5, UTOL5, 3e-7}; // one for each output
  real J, J0, dp;
  real Time, TimeStep;
  real perturb = 1.0e-3;
  xf_Vector *U0, *U, *dU, **Uj;
  xf_Vector **Value_Uj, **Value_Uj0;
  xf_Output *Output;
  xf_All *All;

  // Load scalar diffusion case
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);
  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);
  
  // copy IC over to U0
  ierr = xf_Error(xf_CreateVector(&U0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_CopyVector(All->Mesh, U, U0);
  if (ierr != xf_OK) return ierr;

  // Loop over time schemes
  for (iScheme=0; iScheme<nScheme; iScheme++){

    // current time scheme (must be DGinTime)
    TimeScheme = TimeSchemes[iScheme];
    
    // determine order in time
    ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
    if (ierr != xf_OK) return ierr;
    
    // Time and TimeStep
    Time     = Times[iScheme][0];
    TimeStep = Times[iScheme][1];
    
    // create a dU vector
    sprintf(StateName, "StateTemp");
    ierr = xf_Error(xf_FindSimilarVector(All, U, StateName, xfe_True, xfe_True,  
					 NULL, &dU, NULL));
    if (ierr != xf_OK) return ierr;

    // Create nodal U vectors
    ierr = xf_Error(xf_Alloc( (void **) &Uj, OrderTime+1, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      sprintf(StateName, "State_%d", iOrder);
      ierr = xf_Error(xf_FindSimilarVector(All, U, StateName, xfe_True, xfe_True,  
					   NULL, Uj + iOrder, NULL));
      if (ierr != xf_OK) return ierr;
      // set Uj vectors to U0
      ierr = xf_Error(xf_SetVector(U0, xfe_Set, Uj[iOrder]));
      if (ierr != xf_OK) return ierr;
    }
    // Add to Uj vectors quasi-random perturbations
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      ierr = xf_Error(xf_VectorRand(U, 2+3*iOrder));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorMultSet(U, perturb, xfe_Add, Uj[iOrder]));
      if (ierr != xf_OK) return ierr;
    }
    
    // Create Value_Uj (and Value_Uj0) vectors
    ierr = xf_Error(xf_Alloc( (void **) &Value_Uj, 2*(OrderTime+1), sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (iOrder = 0; iOrder < 2*(OrderTime+1); iOrder++){
      sprintf(StateName, "Value_Uj_%d", iOrder);
      ierr = xf_Error(xf_FindSimilarVector(All, U, StateName, xfe_True, xfe_True,  
					   NULL, Value_Uj + iOrder, NULL));
      if (ierr != xf_OK) return ierr;
    }
    Value_Uj0 = Value_Uj + (OrderTime+1);
  
    // loop over outputs
    for (iOutput=0; iOutput<nOutput; iOutput++){
    
      // set OutputName
      OutputName = OutputNames[iOutput];

      // find output structure
      ierr = xf_Error(xf_FindOutput(All->EqnSet, OutputName, &Output));
      if (ierr != xf_OK) return ierr;
      
      // calculate unsteady output (contributions) and linearization
      Output->Value = J0 = 0;
      for (iOrder = 0; iOrder < OrderTime+1; iOrder++){
	ierr = xf_Error(xf_SetZeroVector(Value_Uj0[iOrder]));
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_IncrementUnsteadyOutputs(All, OutputName, TimeScheme, 
						  OrderTime+1, Uj, U, Time, TimeStep, 
						  xfe_False, xfe_False, &J0, Value_Uj0));
      if (ierr != xf_OK) return ierr; 
      
      // perturb time-nodal states, each in turn
      for (iOrder = 0; iOrder < OrderTime+1; iOrder++){
	ierr = xf_Error(xf_VectorRand(dU, 3+5*iOrder));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VectorMult(dU, perturb));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VectorMultSet(dU, 1.0, xfe_Add, Uj[iOrder]));
	if (ierr != xf_OK) return ierr;
	
	// calculate unsteady output (contributions) and linearization again
	Output->Value = J = 0;
	for (jOrder = 0; jOrder < OrderTime+1; jOrder++){
	  ierr = xf_Error(xf_SetZeroVector(Value_Uj[jOrder]));
	  if (ierr != xf_OK) return ierr;
	}
	ierr = xf_Error(xf_IncrementUnsteadyOutputs(All, OutputName, TimeScheme, 
						    OrderTime+1, Uj, U, Time, TimeStep, 
						    xfe_False, xfe_False, &J, Value_Uj));
	if (ierr != xf_OK) return ierr; 
	
	// check linearization (problem is linear, except for square output)
	ierr = xf_Error(xf_VectorDot(dU, Value_Uj[iOrder], &dp));
	if (ierr != xf_OK) return ierr;
	xf_AssertWithin(J-J0, dp, tolvec[iOutput]); 
	
	// reset Uj back to original
	ierr = xf_Error(xf_VectorMultSet(dU, 1.0, xfe_Sub, Uj[iOrder]));
	if (ierr != xf_OK) return ierr;
      } // iOrder
      
    } // iOutput

    // release memory
    xf_Release( (void **) Uj);
    xf_Release( (void **) Value_Uj);
    
  } // iScheme
  
  // Destroy U0
  ierr = xf_Error(xf_DestroyVector(U0, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;
}
