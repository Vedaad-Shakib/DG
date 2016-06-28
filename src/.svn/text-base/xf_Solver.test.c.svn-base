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
#include "xf_Penalty.h"

// Functions for setting up and running a case within a unit test
#include "xf_UnitRun.c"


TEST_xf_Unsteady_Adjoint_MotionGCL()
{
  /*
    This test runs an unsteady simulation with mesh motion, using the GCL. 
    The compressible Euler equations are tested.  The computation is carried out
    for a specified number of time steps and spatial discretization order.  
    Next, the unsteady-adjoint is computed. The unsteady adjoint solution is 
    tested by using Psi(0) = sum(i>=1) [R(i)_U(0)]^T*Psi(i) to compute the 
    sensitivity of the final-time output to the initial condition U(0) 
    (see description in xf_Solver.h).  This sensitivity is verified by actually 
    perturbing U(0) and re-running the forward-time simulation.

    The time schemes used are DG1 and DG2 in time.
  */

  int ierr, j, nPsi;
  int Order = 2;
  int nTimeStep = 4;
  int iOut, nOut = 3; // number of outputs to test
  int iEqn, nEqn = 2; // number of eqns
  int iScheme, nScheme = 3; // number of schemes to test
  enum xfe_Bool LinearFlag;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  const char *OutputNames[] = {"Point", "BoundaryInt", "BoundaryIntNoFlux"}; 
  const char *OutputName;
  const char *TimeSchemes[] = {"BDF1", "DG1", "DG2"}; 
  const char SavePrefix0[] = "DELETEME_Unsteady_Adjoint_MotionGCL";
  const char *SavePrefix;
  char LogOutput[xf_MAXSTRLEN];
  real J, Jp, dJ, dJAdjoint, dJAdjointGCL;
  real perturb = 1.0e-3, tol, rnorm;
  xf_Vector *U, **Psi, *dU, *PsiGCL=NULL;
  xf_TimeHistData *TimeHistData = NULL;
  xf_All *All;
  xf_Vector *GCL=NULL, *dGCL=NULL;


  // system is NOT linear
  LinearFlag = xfe_False;

  // need a SavePrefix if not linear, so that states are written out
  SavePrefix = ((LinearFlag) ? NULL : SavePrefix0);

  // tolerance for checks
  tol = perturb*perturb; 

  // loop over equation sets to test
  for (iEqn = 0; iEqn < nEqn; iEqn++){ // run iEqn=nEqn for custom.job

    // loop over outputs to test
    for (iOut = 0; iOut < nOut; iOut++){

      /* Loop over the given TimeSchemes */    
      for (iScheme = 0; iScheme < nScheme; iScheme++){

        OutputName = OutputNames[iOut];
      
        xf_printf("Output = %s, TimeScheme = %s\n", OutputName, TimeSchemes[iScheme]);

        if (iEqn == 0){
          // Euler, plunge motion
          xf_printf("\nEuler equation set, ");
          ierr = xf_Error(xf_UnitBoxQ1TriangleAllEuler(&All, xfe_UnitMotionPlunge1));
          xf_AssertEqual(ierr, xf_OK);
        }
        else if (iEqn == 1){
          // Compressible Navier-Stokes, plunge motion
          xf_printf("\nCompressible Navier-Stokes equation set, ");
          ierr = xf_Error(xf_UnitBoxQ1TriangleAllCNS(&All, xfe_UnitMotionPlunge2));
          xf_AssertEqual(ierr, xf_OK);
        }
        else{
          // custom case
          ierr = xf_Error(xf_ReadAllFromJobFile("custom.job", xfe_True, &All));
          if (ierr != xf_OK) return ierr;
          // load dynamic library
          ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
          if (ierr != xf_OK) return ierr;
          All->EqnSet->Dim = All->Mesh->Dim;
          // register eqnset
          ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
          if (ierr != xf_OK) return ierr;
          // get state
          ierr = xf_Error(xf_FindOrCreatePrimalState(All, xfe_False, NULL, &U));
          if (ierr != xf_OK) return ierr;  
          // set outputname
          ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "LogOutput", LogOutput));
          if (ierr != xf_OK) return ierr;
          OutputName = LogOutput;
        }
            
        if (iEqn < nEqn){

          // Load dynamic library, register eqnset, initialize solution
          ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &U)); // use variable order
          xf_AssertEqual(ierr, xf_OK);

          // Set unsteady parameters
          ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "TimeScheme", 
                                         TimeSchemes[iScheme]));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Time", "0.0"));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "EndTime", "1.0"));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nTimeStep",
                                            nTimeStep));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "LinearFlag", 
                                             LinearFlag));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval",
                                            1));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "ResidualTolerance", 
                                         "1e-12"));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "AdjointResidualTolerance", 
                                         "1e-12"));
          if (ierr != xf_OK) return ierr; 
          ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "LogOutput", OutputName));
          if (ierr != xf_OK) return ierr;
      
          ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "AdjointOutputs", OutputName));
          if (ierr != xf_OK) return ierr;

          ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "UseGCL", "True")); 
          if (ierr != xf_OK) return ierr;
        }

        // Initialize GCL vector
        ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U, -1, xfe_True, &GCL));
        if(ierr != xf_OK) return ierr;

        // Need time-history data
        ierr = xf_Error(xf_CreateUniformTimeHistData(All, NULL, &TimeHistData));
        if (ierr != xf_OK) return ierr;
      
        // Apply unsteady time scheme (Time finishes = EndTime)
        ierr = xf_Error(xf_ApplyTimeScheme(All, SavePrefix, xfe_False, &U, TimeHistData));
        if (ierr != xf_OK) return ierr;

        // Unsteady output, J, is stored in the output structure
        ierr = xf_Error(xf_ReadStoredOutput(All->EqnSet, OutputName, &J, NULL));
        if (ierr != xf_OK) return ierr;

        // locate adjoint (set to zero)
        ierr = xf_Error(xf_FindAdjointVectors(All, U, OutputName, xfe_False, xfe_True, 
                                              &nPsi, &Psi, NULL));
        xf_AssertEqual(ierr, xf_OK);

        /* Run unsteady adjoint solve back to Time=0.  Note, TimeWeights in
           TimeHistData remains NULL, which means that Jtot is the output at
           the final time.
        */
        ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, SavePrefix, U, nPsi,
                                                  Psi, TimeHistData));
        if (ierr != xf_OK) return ierr;

        // locate single GCL adjoint (1 output)
        ierr = xf_Error(xf_InitMeshMotionGCLAdjoint(All, Psi[0], -1, xfe_False, &PsiGCL));
        if (ierr != xf_OK) return ierr;

        // Find a dU vector
        ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_False, xfe_True, 
                                             NULL, &dU, NULL));
        if (ierr != xf_OK) return ierr;

        // Find a dGCL vector
        ierr = xf_Error(xf_FindSimilarVector(All, GCL, "dGCL", xfe_False, xfe_True, 
                                             NULL, &dGCL, NULL));
        if (ierr != xf_OK) return ierr;

        // create a small perturbation of the IC: dU(0)
        ierr = xf_Error(xf_VectorRand(dU, 17));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_VectorMult(dU, perturb));
        if (ierr != xf_OK) return ierr;

        // also create a small perturbation to the initial GCL vector
        ierr = xf_Error(xf_VectorRand(dGCL, 11));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_VectorMult(dGCL, perturb));
        if (ierr != xf_OK) return ierr;
      
        // calculate dJAdjoint = Psi(0)^T*dU(0) 
        ierr = xf_Error(xf_VectorDot(Psi[0], dU, &dJAdjoint ));
        if (ierr != xf_OK) return ierr;

        // print out norm of GCL adjoint
        ierr = xf_Error(xf_VectorNorm(PsiGCL, 1, &rnorm));
        if (ierr != xf_OK) return ierr;
        xf_printf("Norm of PsiGCL = %.10E\n", rnorm);


        // add GCL contribution to dJAdjoint
        ierr = xf_Error(xf_VectorDot(PsiGCL, dGCL, &dJAdjointGCL ));
        if (ierr != xf_OK) return ierr;
        dJAdjoint += dJAdjointGCL;

        // Set CFL to 1
        ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "CFL", "1"));
        if (ierr != xf_OK) return ierr;

        // re-initialize U(0) and perturb with dU(0)
        ierr = xf_Error(xf_InitState(All, U));
        if (ierr != xf_OK) return ierr;

        ierr = xf_Error(xf_SetVector(dU, xfe_Add, U));
        if (ierr != xf_OK) return ierr;

        // Also reinitialize GCL and perturb with dGCL
        ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U, -1, xfe_True, &GCL));
        if(ierr != xf_OK) return ierr;

        ierr = xf_Error(xf_SetVector(dGCL, xfe_Add, GCL));
        if (ierr != xf_OK) return ierr;


        // Apply unsteady time scheme again (Time finishes at EndTime)
        ierr = xf_Error(xf_ApplyTimeScheme(All, NULL, xfe_False, &U, TimeHistData));
        if (ierr != xf_OK) return ierr;

        // Unsteady perturbed output, Jp, is stored in the output structure
        ierr = xf_Error(xf_ReadStoredOutput(All->EqnSet, OutputName, &Jp, NULL));
        if (ierr != xf_OK) return ierr;

        // Compute dJ = Jp - J
        dJ = Jp - J;
      
        xf_printf("Predicted dJ: %15.8E \n", dJAdjoint);
        xf_printf("Actual    dJ: %15.8E \n", dJ);
        xf_printf("J = %.10E, Jp = %.10E\n", J, Jp);

        // Verify dJ answer is the same (to within tolerance)
        xf_AssertWithin(dJ, dJAdjoint, tol); 

        xf_Release( (void *) Psi);

        // Destroy TimeHistData
        ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
        if (ierr != xf_OK) return ierr;
    
        // Destroy All
        ierr = xf_Error(xf_DestroyAll(All));
        xf_AssertEqual(ierr, xf_OK);

      } // iScheme

    } // iOut

  } // iEqn

  return xf_OK;

}



TEST_xf_ApplyTimeScheme_ConvRate()
{
  /*
    Tests the convergence rate of various unsteady solvers, for both
    the scalar and the Euler equations.
  */

  int ierr;
  int Order = 1;
  int nTimeStep = 32;
  int iEqn, nEqn = 1; // number of equation sets to test
  int iScheme, nScheme = 8; // number of schemes to test
  int iRefine, nRefine = 4; // number of uniform temporal refinements
  enum xfe_Bool LinearFlag;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  char *TempFcnName = NULL, ScalarFcnName[] = "OscillatingState";
  char *TempFcnData = NULL, ScalarFcnData[] = "1.0 0.2 6 0.01";
  const char *OutputNames[] = {"HeatFlowUnsteady"}; 
  const char *OutputName;
  const char *TimeSchemes[] = {"BDF1", "BDF2", "DG1", "DG2", "ESDIRK3", "ESDIRK4", "DIRK3", "DIRK4"}; 
  const char *SavePrefix = NULL;
  real J[4][xfe_TimeSchemeLast]; // nRefine x nScheme
  real Jexact ,err0, err1, rate;
  real rates[8] = {6.4, 2.4, 3.0, 5.0, 3.0, 4.0, 2.8, 4.0};
  real tolerance = 0.15;
  xf_Vector *U, *dU;
  xf_TimeHistData *TimeHistData = NULL;
  xf_All *All;

  // loop over equation sets to test
  for (iEqn = 0; iEqn < nEqn; iEqn++){

    OutputName = OutputNames[iEqn];

    // Load appropriate case
    if (iEqn == 0){
      // Scalar diffusion
      ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
      // Load dynamic library, register eqnset, initialize solution
      ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
      xf_AssertEqual(ierr, xf_OK);
      // system is linear
      LinearFlag = xfe_True;
      // make boundary condition unsteady
      TempFcnName = All->EqnSet->BCs->BC[0].Function;
      TempFcnData = All->EqnSet->BCs->BC[0].Data;
      All->EqnSet->BCs->BC[0].Function = ScalarFcnName;
      All->EqnSet->BCs->BC[0].Data = ScalarFcnData;
    }
    else{
      return xf_Error(xf_NOT_SUPPORTED);
    }

    // loop over refinements
    for (iRefine = 0; iRefine < nRefine; iRefine++){
      
      // refine temporal grid after the first iteration
      if (iRefine > 0) nTimeStep *= 2;
      
      /* Loop over the given TimeSchemes */
      for (iScheme = 0; iScheme < nScheme; iScheme++){
    
	// re-initialize U(0) for schemes beyond the first one
	if (iScheme > 0){
	  ierr = xf_Error(xf_InitState(All, U));
	  if (ierr != xf_OK) return ierr;
	}

	// Set unsteady parameters
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "TimeScheme", 
				       TimeSchemes[iScheme]));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Time", "0.0"));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "EndTime", "1.0"));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nTimeStep",
					  nTimeStep));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "LinearFlag", 
					   LinearFlag));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "ResidualTolerance", 
				       "1e-12"));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "LogOutput", OutputName));
	if (ierr != xf_OK) return ierr;
      

	// Need time-history data
	ierr = xf_Error(xf_CreateUniformTimeHistData(All, NULL, &TimeHistData));
	if (ierr != xf_OK) return ierr;
      
	// Apply unsteady time scheme (Time finishes = EndTime)
	ierr = xf_Error(xf_ApplyTimeScheme(All, SavePrefix, xfe_False, &U, TimeHistData));
	if (ierr != xf_OK) return ierr;

      
	// Unsteady output, J, is stored in the output structure
	ierr = xf_Error(xf_ReadStoredOutput(All->EqnSet, OutputName, &J[iRefine][iScheme], NULL));
	if (ierr != xf_OK) return ierr;

	// Destroy TimeHistData
	ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
	if (ierr != xf_OK) return ierr;

      } // iScheme

    } // iRefine

    /* Loop over the given TimeSchemes */
    for (iScheme = 0; iScheme < nScheme; iScheme++){
      
      // Verify convergence rate
      Jexact = J[nRefine-1][iScheme]; // last one in refinement sequence
      err0 = J[0][iScheme]-Jexact;
      err1 = J[1][iScheme]-Jexact;
      rate = log(fabs(err0/err1))/log(2.0);
      xf_printf("iEqn = %d, Scheme = %s, rate = %.3f (expected %.3f)\n", iEqn, TimeSchemes[iScheme], rate, rates[iScheme]);
      xf_AssertWithin(rate, rates[iScheme] , tolerance); 
    } // iScheme

    if (iEqn==0){
      All->EqnSet->BCs->BC[0].Function = TempFcnName;
      All->EqnSet->BCs->BC[0].Data = TempFcnData;
    }
    
    
    // Destroy All
    ierr = xf_Error(xf_DestroyAll(All));
    xf_AssertEqual(ierr, xf_OK);

  } // iEqn


  return xf_OK;  
}



TEST_xf_ScalarDiff_Steady_Adjoint()
{
  /*
    This test runs a steady scalar diffusion problem: heat conduction
    from a hot source on the left to a sink on the right of a square
    domain.  The computed integrated heat flux, "HeatFlow", is
    verified against the analytical value.  Next, the steady-state
    adjoint is computed: its accuracy is verified by adding an
    arbitrary source to the equations, re-solving for U, and checking
    that:

       dJ = J(U(S)) - J(U(0)) = Psi^T*S

    Where J is the output, U is the solution, Psi is the adjoint, and
    S is the arbitrary source.  Since the underlying primal problem is
    linear, the above equality should hold up to round-off.
  */

  int ierr, nPsi;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 1;
  real J, Jp, dJ, dJAdjoint, vnorm;
  xf_Vector *U, **Psi, *S;
  xf_All *All;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_UnitMotionNone));
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
  xf_AssertWithin(J, 0.1, UTOL2); // diffusion constant is 0.1

  // locate adjoint (set to zero)
  ierr = xf_Error(xf_FindAdjointVectors(All, U, "HeatFlow", xfe_False, xfe_True, 
					&nPsi, &Psi, NULL));
  xf_AssertEqual(ierr, xf_OK);

  // solve for adjoint
  ierr = xf_Error(xf_SolveAdjoints(All, 0.0, 1.0, xfe_False, U, nPsi, NULL, Psi, NULL, xfe_False, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
 
  
  // Find a source vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Source", xfe_False, xfe_True, NULL, &S, NULL));
  if (ierr != xf_OK) return ierr;

  // Set source to small varying perturbed values throughout the domain
  ierr = xf_Error(xf_VectorRand(S, 17));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMult(S, 1.0e-4));
  if (ierr != xf_OK) return ierr;

  // Solve system with source
  ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, S, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Compute output with perturbed solution
  ierr = xf_Error(xf_CalculateOutput(All, "HeatFlow", U, &Jp, NULL, xfe_Set));
  xf_AssertEqual(ierr, xf_OK);

  // Compute dJ = Jp - J
  dJ = Jp - J;

  // Compute Psi^T * S, which should be dJ (exactly since problem is linear)
  ierr = xf_Error(xf_VectorDot(Psi[0], S, &dJAdjoint ));
  if (ierr != xf_OK) return ierr;
  
  // Verify dJ answer is the same
  xf_AssertWithin(dJ, dJAdjoint, UTOL3); 

  xf_Release( (void *) Psi);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);


  return xf_OK;  
}


TEST_xf_Unsteady_Adjoint()
{
  /*
    This test runs an unsteady simulation. In the scalar case, the
    equation is pure diffusion: heat conduction from a hot source
    (temperature = 1) on the left to a sink (temperature = 0) on the
    right of a square domain, with the intial condition consisting of
    temperature = 0.5 everywhere.  The compressible Euler equations
    are also tested.  The computation is carried out for a specified
    number of time steps and spatial discretization order.  Next, the
    unsteady-adjoint is computed -- when the problem is linear, we do
    not write/read state vectors to/from disk.  The output of interest
    is the integrated heat flux at the final time.  The unsteady
    adjoint solution is tested by using Psi(0) = sum(i>=1)
    [R(i)_U(0)]^T*Psi(i) to compute the sensitivity of the final-time
    output to the initial condition U(0) (see description in
    xf_Solver.h).  This sensitivity is verified by actually perturbing
    U(0) and re-running the forward-time simulation.

    Several time shcemes are tested, including both multistep and
    DGinTime.

    When the underlying primal problem is linear, the sensitivity test
    should hold to round-off.  When it is not, it is tested to a
    prescribed relative tolerance.
  */

  int ierr, nPsi;
  int Order = 2;
  int nTimeStep = 4;
  int iEqn, nEqn = 2; // number of equation sets to test
  int iScheme, nScheme = 4; // number of schemes to test
  enum xfe_Bool LinearFlag;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  const char *OutputNames[] = {"HeatFlowUnsteady", "Point"}; 
  const char *OutputName;
  const char *TimeSchemes[] = {"BDF1", "BDF2", "DG1", "DG2"}; 
  const char SavePrefix0[] = "DELETEME_Unsteady_Adjoint";
  const char *SavePrefix;
  real J, Jp, dJ, dJAdjoint;
  real perturb = 1.0e-3, tol;
  xf_Vector *U, **Psi, *dU;
  xf_TimeHistData *TimeHistData = NULL;
  xf_All *All;

  // loop over equation sets to test
  for (iEqn = 0; iEqn < nEqn; iEqn++){

    OutputName = OutputNames[iEqn];

    // Load appropriate case
    if (iEqn == 0){
      // Scalar diffusion
      ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
      // Load dynamic library, register eqnset, initialize solution
      ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
      xf_AssertEqual(ierr, xf_OK);
      // system is linear
      LinearFlag = xfe_True;
      tol = UTOL5;
    }
    else{
      // Euler
      ierr = xf_Error(xf_UnitBoxQ1TriangleAllEuler(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
      // Load dynamic library, register eqnset, initialize solution
      ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &U)); // use variable order
      xf_AssertEqual(ierr, xf_OK);
      // system is NOT linear
      LinearFlag = xfe_False;
      tol = 0.1*perturb*perturb;
    }

    // need a SavePrefix if not linear, so that states are written out
    SavePrefix = ((LinearFlag) ? NULL : SavePrefix0);

    /* Loop over the given TimeSchemes */
    
    for (iScheme = 0; iScheme < nScheme; iScheme++){
    
      // re-initialize U(0) for schemes beyond the first one
      if (iScheme > 0){
	ierr = xf_Error(xf_InitState(All, U));
	if (ierr != xf_OK) return ierr;
      }

      // Set unsteady parameters
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "TimeScheme", 
				     TimeSchemes[iScheme]));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Time", "0.0"));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "EndTime", "1.0"));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nTimeStep",
					nTimeStep));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "LinearFlag", 
					 LinearFlag));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval",
					1));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "ResidualTolerance", 
				     "1e-12"));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "AdjointResidualTolerance", 
				     "1e-12"));
      if (ierr != xf_OK) return ierr; 
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "LogOutput", OutputName));
      if (ierr != xf_OK) return ierr;
     
      
      // Need time-history data
      ierr = xf_Error(xf_CreateUniformTimeHistData(All, NULL, &TimeHistData));
      if (ierr != xf_OK) return ierr;
      
      // Apply unsteady time scheme (Time finishes = EndTime)
      ierr = xf_Error(xf_ApplyTimeScheme(All, SavePrefix, xfe_False, &U, TimeHistData));
      if (ierr != xf_OK) return ierr;

      
      // Unsteady output, J, is stored in the output structure
      ierr = xf_Error(xf_ReadStoredOutput(All->EqnSet, OutputName, &J, NULL));
      if (ierr != xf_OK) return ierr;

      // locate adjoint (set to zero)
      ierr = xf_Error(xf_FindAdjointVectors(All, U, OutputName, xfe_False, xfe_True, 
					    &nPsi, &Psi, NULL));
      xf_AssertEqual(ierr, xf_OK);

      /* Run unsteady adjoint solve back to Time=0.  Note, TimeWeights in
	 TimeHistData remains NULL, which means that Jtot is the output at
	 the final time.
      */
      ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, SavePrefix, U, nPsi,
						Psi, TimeHistData));
      if (ierr != xf_OK) return ierr;

      // Find a dU vector
      ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_False, xfe_True, 
					   NULL, &dU, NULL));
      if (ierr != xf_OK) return ierr;

      // create a small perturbation of the IC: dU(0)
      ierr = xf_Error(xf_VectorRand(dU, 17));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorMult(dU, perturb));
      if (ierr != xf_OK) return ierr;

      // calculate dJAdjoint = Psi(0)^T*dU(0)
      ierr = xf_Error(xf_VectorDot(Psi[0], dU, &dJAdjoint ));
      if (ierr != xf_OK) return ierr;

      // re-initialize U(0) and perturb with dU(0)
      ierr = xf_Error(xf_InitState(All, U));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_SetVector(dU, xfe_Add, U));
      if (ierr != xf_OK) return ierr;

      // Apply unsteady time scheme again (Time finishes at EndTime)
      ierr = xf_Error(xf_ApplyTimeScheme(All, NULL, xfe_False, &U, TimeHistData));
      if (ierr != xf_OK) return ierr;

      // Unsteady perturbed output, Jp, is stored in the output structure
      ierr = xf_Error(xf_ReadStoredOutput(All->EqnSet, OutputName, &Jp, NULL));
      if (ierr != xf_OK) return ierr;

      // Compute dJ = Jp - J
      dJ = Jp - J;

      // Verify dJ answer is the same (to within tolerance)
      xf_AssertWithin(dJ, dJAdjoint, tol); 

      xf_Release( (void *) Psi);
      
      // Destroy TimeHistData
      ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
      if (ierr != xf_OK) return ierr;

    } // iScheme
    
    // Destroy All
    ierr = xf_Error(xf_DestroyAll(All));
    xf_AssertEqual(ierr, xf_OK);

  } // iEqn


  return xf_OK;  
}

// TEST WAS HERE

TEST_xf_Unsteady_Adjoint_ErrEst()
{
  /*
    Tests error estimation using unsteady adjoints.  Only temporal
    error estimation is checked.
  */

  int ierr, nPsi;
  int Order = 3;
  int nTimeStep = 4;
  int UErrEstnSubSlab;
  int UErrEstOrderIncrement;
  int iEqn, nEqn = 2; // number of equation sets to test
  int iScheme, nScheme = 2; // number of schemes to test
  int OrderInc; 
  enum xfe_Bool LinearFlag;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  const char *OutputNames[] = {"HeatFlowUnsteady", "Point"}; 
  const char *OutputName;
  const char *TimeSchemes[] = {"DG1", "DG2"}; 
  const char SavePrefix0[] = "DELETEME_Unsteady_Adjoint";
  const char *SavePrefix;
  real JH, Jh, dJ, dJEst;
  real perturb = 1.0e-3, tol;
  xf_Vector *U, **Psi, *dU;
  xf_TimeHistData *TimeHistData = NULL;
  xf_All *All;

  // loop over equation sets to test
  for (iEqn = 0; iEqn < nEqn; iEqn++){

    OutputName = OutputNames[iEqn];

    // Load appropriate case
    if (iEqn == 0){
      // Scalar diffusion
      ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
      // Load dynamic library, register eqnset, initialize solution
      ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
      xf_AssertEqual(ierr, xf_OK);
      // system is linear
      LinearFlag = xfe_True;
      tol = 1e-12;
    }
    else{
      // Euler
      ierr = xf_Error(xf_UnitBoxQ1TriangleAllEuler(&All, xfe_UnitMotionNone));
      xf_AssertEqual(ierr, xf_OK);
      // Load dynamic library, register eqnset, initialize solution
      ierr = xf_Error(xf_InitializeTestRun(All, Basis, -Order, &U)); // variable order
      xf_AssertEqual(ierr, xf_OK);
      // system is NOT linear
      LinearFlag = xfe_False;
      tol = 1e-8;
    }

    // need a SavePrefix for error estimation so that states are written out
    SavePrefix = SavePrefix0;

    // Find a dU vector
    ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_False, xfe_True, 
					 NULL, &dU, NULL));
    if (ierr != xf_OK) return ierr;

    /* Loop over the given TimeSchemes */
    for (iScheme = 0; iScheme < nScheme; iScheme++){

      /* Loop over order increase versus slab-refine choice */
      for (OrderInc=0; OrderInc<2; OrderInc++){

	// only try OrderInc for DG1
	if (OrderInc && (iScheme != 0)) continue;
	
	UErrEstnSubSlab = ((OrderInc) ? 1 : 2);
	UErrEstOrderIncrement = ((OrderInc) ? 1 : 0);

	// create a small perturbation of the IC: dU(0)
	ierr = xf_Error(xf_VectorRand(dU, 17));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VectorMult(dU, perturb));
	if (ierr != xf_OK) return ierr;
	
	// initialize U(0) and perturb with dU(0)
	ierr = xf_Error(xf_InitState(All, U));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_SetVector(dU, xfe_Add, U));
	if (ierr != xf_OK) return ierr;
	
	
	// Set unsteady parameters
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "TimeScheme", 
				       TimeSchemes[iScheme]));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Time", "0.0"));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "EndTime", "1.0"));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nTimeStep",
					  nTimeStep));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "LinearFlag", 
					   LinearFlag));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval",
					  1));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "ResidualTolerance", 
				       "1e-12"));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "AdjointResidualTolerance", 
				       "1e-12"));
	if (ierr != xf_OK) return ierr; 
	ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "LogOutput", OutputName));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "UErrEstOn", 
					   xfe_True));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "UErrEstOrderIncrement",
					  UErrEstOrderIncrement));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement",
					  0));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "UErrEstnSubSlab",
					  UErrEstnSubSlab));
	if (ierr != xf_OK) return ierr;
	  
	// Need time-history data
	ierr = xf_Error(xf_CreateUniformTimeHistData(All, NULL, &TimeHistData));
	if (ierr != xf_OK) return ierr;
	
	// Apply unsteady time scheme (Time finishes = EndTime)
	ierr = xf_Error(xf_ApplyTimeScheme(All, SavePrefix, xfe_False, &U, TimeHistData));
	if (ierr != xf_OK) return ierr;
	
	
	// Unsteady output, JH, is stored in the output structure
	ierr = xf_Error(xf_ReadStoredOutput(All->EqnSet, OutputName, &JH, NULL));
	if (ierr != xf_OK) return ierr;
	
	// locate adjoint (set to zero)
	ierr = xf_Error(xf_FindAdjointVectors(All, U, OutputName, xfe_False, xfe_True, 
					      &nPsi, &Psi, NULL));
	xf_AssertEqual(ierr, xf_OK);
	
	/* Run unsteady adjoint solve back to Time=0.  Error estimation
	   is performed here.
	*/
	ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, SavePrefix, U, nPsi,
						  Psi, TimeHistData));
	if (ierr != xf_OK) return ierr;
	
	// unsteady error estimate is stored in the output structure
	ierr = xf_Error(xf_ReadStoredOutput(All->EqnSet, OutputName, NULL, &dJEst));
	if (ierr != xf_OK) return ierr;
	
	// re-initialize U(0) (with perturbation)
	ierr = xf_Error(xf_InitState(All, U));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetVector(dU, xfe_Add, U));
	if (ierr != xf_OK) return ierr;
	
	// Destroy TimeHistData
	ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
	if (ierr != xf_OK) return ierr;
	
	// prepare a finer unsteady discretization
	if (OrderInc){
	  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "TimeScheme",
					 TimeSchemes[iScheme+1]));
	  if (ierr != xf_OK) return ierr;
	}
	else{
	  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nTimeStep",
					    nTimeStep*UErrEstnSubSlab));
	  if (ierr != xf_OK) return ierr;
 	}
	
	// create a time history with the finer discretization
	ierr = xf_Error(xf_CreateUniformTimeHistData(All, NULL, &TimeHistData));
	if (ierr != xf_OK) return ierr;
	
	// Apply a finer unsteady time scheme (Time still finishes at EndTime)
	ierr = xf_Error(xf_ApplyTimeScheme(All, NULL, xfe_False, &U, TimeHistData));
	if (ierr != xf_OK) return ierr;
	
	// Unsteady finer, Jh, is stored in the output structure
	ierr = xf_Error(xf_ReadStoredOutput(All->EqnSet, OutputName, &Jh, NULL));
	if (ierr != xf_OK) return ierr;
	
	// Compute error: dJ = JH - Jh
	dJ = JH - Jh;
	
	/*
 	xf_printf("TimeScheme = %s, OrderInc = %d, dJEst = %.10E, dJ = %.10E, tol = %.10E\n", 
 		  TimeSchemes[iScheme], OrderInc, dJEst, dJ, tol); 
  */
	
	// Verify dJ answer is the same to within the prescribed tolerance
	xf_AssertWithin(dJ, dJEst, tol);
	
	xf_Release( (void *) Psi);
	
	// Destroy TimeHistData
	ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
	if (ierr != xf_OK) return ierr;
	
      } // OrderInc
      
    } // iScheme
    
    // Destroy All
    ierr = xf_Error(xf_DestroyAll(All));
    xf_AssertEqual(ierr, xf_OK);

  } // iEqn


  return xf_OK;  
}


TEST_xf_DGTimeInterpolate()
{
  /*
    Tests DG state interpolation
  */

  int ierr;
  int Order = 2;
  int iOrder, OrderTime;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  char StateName[xf_MAXSTRLEN];
  enum xfe_TimeSchemeType TimeScheme;
  real norm, xt;
  xf_Vector *Utemp, *U, **Uj;
  xf_All *All;

  // Load scalar diffusion case
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_UnitMotionNone));
  xf_AssertEqual(ierr, xf_OK);
  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // DG1
  TimeScheme = xfe_TimeSchemeDG1;
  OrderTime  = 1;

  // Create nodal U vectors
  ierr = xf_Error(xf_Alloc( (void **) &Uj, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(StateName, "State_%d", iOrder);
    ierr = xf_Error(xf_FindSimilarVector(All, U, StateName, xfe_True, xfe_True,  
					 NULL, Uj + iOrder, NULL));
    if (ierr != xf_OK) return ierr;
    // set Uj vectors to U*(iOrder+1)
    ierr = xf_Error(xf_VectorMultSet(U, (real) iOrder+1, xfe_Set, Uj[iOrder]));
    if (ierr != xf_OK) return ierr;
  }
  
  // loate Utemp
  sprintf(StateName, "StateTemp");
  ierr = xf_Error(xf_FindSimilarVector(All, U, StateName, xfe_True, xfe_True,  
				       NULL, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;
  if (ierr != xf_OK) return ierr;
  
  // Interpolate to left time node
  xt = 0.;
  ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Uj, 1, &xt, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;
  // check if same as Uj[0]
  ierr = xf_Error(xf_SetVector(Uj[0], xfe_Sub, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorNorm(Utemp, 1, &norm));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(norm, 0., UTOL5);

  // Interpolate to right time node
  xt = 1.0;
  ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Uj, 1, &xt, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;
  // check if same as Uj[1]
  ierr = xf_Error(xf_SetVector(Uj[1], xfe_Sub, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorNorm(Utemp, 1, &norm));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(norm, 0., UTOL5);

  // Interpolate to 60% from left to right
  xt = 0.6;
  ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Uj, 1, &xt, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;
  // check if same as 0.4*Uj[0] + 0.6*Uj[1]
  ierr = xf_Error(xf_VectorMultSet(Uj[0], 0.4, xfe_Sub, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(Uj[1], 0.6, xfe_Sub, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorNorm(Utemp, 1, &norm));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(norm, 0., UTOL5);
  
  // release memory
  xf_Release( (void **) Uj);
    
  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}




TEST_xf_ErrEstSpaceTimeAniso()
{
  /*
    Tests DG space-time anisotropy measure (jumps)
  */

  int ierr;
  int Order = 2;
  int iOrder, OrderTime;
  enum xfe_BasisType Basis = xfe_QuadLagrange;
  char StateName[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  enum xfe_TimeSchemeType TimeScheme;
  real norm, xt;
  xf_Vector *U, **Uj, *SpaceTimePref;
  xf_All *All;

  // Load scalar convection diffusion case
  ierr = xf_Error(xf_UnitBoxQ1Quad9All(&All));
  xf_AssertEqual(ierr, xf_OK);
  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // DG1
  TimeScheme = xfe_TimeSchemeDG1;
  OrderTime  = 1;

  // Create nodal U vectors
  ierr = xf_Error(xf_Alloc( (void **) &Uj, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(StateName, "State_%d", iOrder);
    ierr = xf_Error(xf_FindSimilarVector(All, U, StateName, xfe_True, xfe_True,  
					 NULL, Uj + iOrder, NULL));
    if (ierr != xf_OK) return ierr;
    // set Uj vectors to U
    ierr = xf_Error(xf_SetVector(U, xfe_Set, Uj[iOrder]));
    if (ierr != xf_OK) return ierr;
  }

  // create a temporal variation
  ierr = xf_Error(xf_VectorMult(Uj[0], 0.8));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMult(Uj[1], 1.2));
  if (ierr != xf_OK) return ierr;

  // locate space-time pref vector
  sprintf(Title, "SpaceTimePref_State");
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, NULL, 
				NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, xfe_True, NULL, 
				&SpaceTimePref, NULL));
  if (ierr != xf_OK) return ierr;

  // estimate anisotropy
  ierr = xf_Error(xf_ErrEstSpaceTimeAniso(All, TimeScheme, Uj, U, 3.0, 0.5, SpaceTimePref, NULL));
  if (ierr != xf_OK) return ierr;

  // verify pref on a few elems
  xf_AssertWithin(SpaceTimePref->GenArray[0].rValue[0][0], 4./9., UTOL1);
  xf_AssertWithin(SpaceTimePref->GenArray[0].rValue[0][1],   1.0, UTOL1);
  xf_AssertWithin(SpaceTimePref->GenArray[0].rValue[0][2], 4./9., UTOL1);

  // release memory
  xf_Release( (void **) Uj);
    
  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_CalcPenaltyGradient_Euler_Quad()
{
  
  /* Penalty function Grdient pinging a compressible Euler flow in a quad Q1 box.*/
  
  int ierr, egrp, elem, r;
  enum xfe_BasisType Basis = xfe_QuadLagrange;
  int Order = 0;
  real Uorig, delta, eps, error1, error2, Pplus, Pminus, g1, g2;
  xf_Vector *U, *P, *GP;
  xf_Data *D;
  xf_All *All;
  xf_Mesh *Mesh;
  
  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1Quad9EulerPenaltyAll(&All));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;
  
  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);
  
  ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D,NULL));
  xf_AssertEqual(ierr, xf_OK);
  
  U = (xf_Vector *)D->Data;
  
  ierr = xf_Error(xf_FindVector(All, "PenaltyIntegral", xfe_LinkageGlobElem, 
                                1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, 
                                xfe_SizeReal, xfe_False,  xfe_True, NULL, 
                                &P, NULL));
  xf_AssertEqual(ierr, xf_OK);
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "GradPenalty", xfe_False, 
                                       xfe_True, NULL, &GP, NULL));
  xf_AssertEqual(ierr, xf_OK);
  
  
  //calculate the penalty function gradient
  ierr = xf_Error(xf_CalculatePenaltyGradient(All, U, GP));
  xf_AssertEqual(ierr, xf_OK);
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      for (r = 0; r < U->GenArray[egrp].r; r++){
        //Saving original U value
        Uorig = U->GenArray[egrp].rValue[elem][r];
        //Perturbing the state
        eps = 1e-4;
        U->GenArray[egrp].rValue[elem][r] = Uorig + eps;
        
        ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
        xf_AssertEqual(ierr, xf_OK);
        
        Pplus = P->GenArray[egrp].rValue[elem][0];
        
        U->GenArray[egrp].rValue[elem][r] = Uorig - eps;
        
        ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
        xf_AssertEqual(ierr, xf_OK);
        
        Pminus = P->GenArray[egrp].rValue[elem][0];
        
        g1 = (Pplus-Pminus)/(2.0*eps);
        error1 = fabs(GP->GenArray[egrp].rValue[elem][r]-g1);
        
        eps *= 0.5;//dividing epsilon by half and verifying
        U->GenArray[egrp].rValue[elem][r] = Uorig + eps;
        
        ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
        xf_AssertEqual(ierr, xf_OK);
        
        Pplus = P->GenArray[egrp].rValue[elem][0];
        
        U->GenArray[egrp].rValue[elem][r] = Uorig - eps;
        
        ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
        xf_AssertEqual(ierr, xf_OK);
        
        Pminus = P->GenArray[egrp].rValue[elem][0];
        
        g2 = (Pplus-Pminus)/(2.0*eps);
        error2 = fabs(GP->GenArray[egrp].rValue[elem][r]-g2);
        
        xf_AssertWithin(error1/error2,4.0,0.1);
        /* printf("GP: %1.10e g1: %1.10e g2: %1.10e error1: %1.10e error2: %1.10e ratio: %1.10e\n", 
                                       GP->GenArray[egrp].rValue[elem][r], g1, g2, error1, error2, error1/error2); */
        //Setting the U value back to original
        U->GenArray[egrp].rValue[elem][r] = Uorig;
      }
    }
  }
  
  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;  
}

 TEST_xf_CalcPenaltyHessian_Euler_Quad()
{  
  /* Penalty function Hessian pinging a compressible Euler flow in a quad Q1 box.*/

  int ierr, egrp, elem, i, j, sr, I, J, si, sj, nn;
  enum xfe_BasisType Basis = xfe_QuadLagrange;
  int Order = 1, index = 0;
  real UorigI, UorigJ, delta, eps, error1, error2;
  real Pip1jp1, Pip1jm1, Pim1jm1, Pim1jp1, Pij;
  real Pip1j, Pijp1, Pim1j, Pijm1;
  real HpApprox1, HpApprox2, ratio;
  xf_Vector *U, *P;
  xf_Data *D;
  xf_All *All;
  xf_Mesh *Mesh;
  xf_Matrix *Hp;
  
  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1Quad9EulerPenaltyAll(&All));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;
  
  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);
  
  ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D,NULL));
  xf_AssertEqual(ierr, xf_OK);
  
  U = (xf_Vector *)D->Data;
  
  ierr = xf_Error(xf_FindVector(All, "PenaltyIntegral", xfe_LinkageGlobElem, 
                                1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, 
                                xfe_SizeReal, xfe_False,  xfe_True, NULL, 
                                &P, NULL));
  xf_AssertEqual(ierr, xf_OK);
  
  sr = All->EqnSet->StateRank;
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    //allocate a Hessian for this element group
    ierr = xf_Error(xf_FindPenaltyHessian(All, egrp, U, &Hp));
    xf_AssertEqual(ierr, xf_OK);
    
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      xf_CalculatePenaltyHessian(All, U, Hp, egrp, elem);
      xf_AssertEqual(ierr, xf_OK);
      nn = U->GenArray[egrp].r/sr;
      //xf_printf("rank %d\n",U->GenArray[egrp].r);
      for (i = 0; i < U->GenArray[egrp].r; i++){
        for (j = 0; j < U->GenArray[egrp].r; j++){
          eps = 1e-2;
          //Saving original U value
          UorigI = U->GenArray[egrp].rValue[elem][i];
          UorigJ = U->GenArray[egrp].rValue[elem][j];
          
          I = (int)i/sr;
          J = (int)j/sr;
          si = i-I*sr;
          sj = j-J*sr;
          //xf_printf("I: %d J: %d si: %d sj: %d\n",I,J,si,sj);
          
          //****************************************************************************
          //getting Pim1jp1
          U->GenArray[egrp].rValue[elem][i] = UorigI - eps;
          U->GenArray[egrp].rValue[elem][j] = UorigJ + eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pim1jp1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pijp1
          U->GenArray[egrp].rValue[elem][j] = UorigJ + eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pijp1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pip1jp1
          U->GenArray[egrp].rValue[elem][i] = UorigI + eps;
          U->GenArray[egrp].rValue[elem][j] = UorigJ + eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pip1jp1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pim1j
          U->GenArray[egrp].rValue[elem][i] = UorigI - eps;
          
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pim1j = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          
          //****************************************************************************
          //getting Pij
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pij = P->GenArray[egrp].rValue[elem][0];
          
          //****************************************************************************
          //getting Pip1j
          U->GenArray[egrp].rValue[elem][i] = UorigI + eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pip1j = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          
          //****************************************************************************
          //getting Pim1jm1
          U->GenArray[egrp].rValue[elem][i] = UorigI - eps;
          U->GenArray[egrp].rValue[elem][j] = UorigJ - eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pim1jm1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pijm1
          U->GenArray[egrp].rValue[elem][j] = UorigJ - eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pijm1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pip1jm1
          U->GenArray[egrp].rValue[elem][i] = UorigI + eps;
          U->GenArray[egrp].rValue[elem][j] = UorigJ - eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pip1jm1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //diagonal
          if(i == j){
            HpApprox1 = (Pip1j-2.0*Pij+Pim1j)/(eps*eps);
          }
          else{
            HpApprox1 = (Pip1jp1-Pim1jp1-Pip1jm1+Pim1jm1)/(4.0*eps*eps);
          }
          
          eps *= 0.5;
          
          //****************************************************************************
          //getting Pim1jp1
          U->GenArray[egrp].rValue[elem][i] = UorigI - eps;
          U->GenArray[egrp].rValue[elem][j] = UorigJ + eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pim1jp1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pijp1
          U->GenArray[egrp].rValue[elem][j] = UorigJ + eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pijp1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pip1jp1
          U->GenArray[egrp].rValue[elem][i] = UorigI + eps;
          U->GenArray[egrp].rValue[elem][j] = UorigJ + eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pip1jp1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pim1j
          U->GenArray[egrp].rValue[elem][i] = UorigI - eps;
          
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pim1j = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          
          //****************************************************************************
          //getting Pij
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pij = P->GenArray[egrp].rValue[elem][0];
                    
          //****************************************************************************
          //getting Pip1j
          U->GenArray[egrp].rValue[elem][i] = UorigI + eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pip1j = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          
          //****************************************************************************
          //getting Pim1jm1
          U->GenArray[egrp].rValue[elem][i] = UorigI - eps;
          U->GenArray[egrp].rValue[elem][j] = UorigJ - eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pim1jm1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pijm1
          U->GenArray[egrp].rValue[elem][j] = UorigJ - eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pijm1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //****************************************************************************
          //getting Pip1jm1
          U->GenArray[egrp].rValue[elem][i] = UorigI + eps;
          U->GenArray[egrp].rValue[elem][j] = UorigJ - eps;
          
          ierr = xf_Error(xf_CalculatePenalty(All, U, P, NULL));
          xf_AssertEqual(ierr, xf_OK);
          
          Pip1jm1 = P->GenArray[egrp].rValue[elem][0];
          
          U->GenArray[egrp].rValue[elem][i] = UorigI;
          U->GenArray[egrp].rValue[elem][j] = UorigJ;
          
          //diagonal
          if(i == j){
            HpApprox2 = (Pip1j-2.0*Pij+Pim1j)/(eps*eps);
          }
          else{
            HpApprox2 = (Pip1jp1-Pim1jp1-Pip1jm1+Pim1jm1)/(4.0*eps*eps);
          }          
          
          
          error1 = fabs(Hp->GenArray[0].rValue[elem][I*nn*sr*sr + J*sr*sr + si*sr + sj]-HpApprox1);
          error2 = fabs(Hp->GenArray[0].rValue[elem][I*nn*sr*sr + J*sr*sr + si*sr + sj]-HpApprox2);
          ratio = error1/error2;
          xf_AssertWithin(ratio,4.0,0.1);
          index++;
          //xf_printf("index: %d Hp: %1.10e HpApprox1: %1.10e HpApprox2: %1.10e error1: %1.10e error2: %1.10e ratio: %1.10e\n", index,
          //       Hp->GenArray[0].rValue[elem][I*nn*sr*sr + J*sr*sr + si*sr + sj], HpApprox1, HpApprox2, error1, error2, ratio);
          
        }
      }
      //xf_printf("another element\n");
    }
    
    ierr = xf_Error(xf_DestroyMatrix(Hp));
    xf_AssertEqual(ierr, xf_OK);
  }
  
  return xf_OK; 
}

