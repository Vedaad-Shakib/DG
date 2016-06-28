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

/*
  FILE:  xf_SolverpMG.c

  This file contains p-Multigrid solver functions.

*/



/******************************************************************/
//   FUNCTION Definition: xf_MGSmooth
static int
xf_MGSmooth(xf_All *All, real c, enum xfe_Bool LinearFlag,
	    enum xfe_Bool ErrorSolve, xf_Vector *S, int nSmooth, xf_Vector *U)
{ 
/*
PURPOSE:

  Wrapper for a "smoothing" call to a nonlinear solver.

INPUTS:

  All: All structure
  c : constant in front of M*U product
  LinearFlag : If True, the Jacobian matrix will not be recalculated if
               it already exists.
  ErrorSolve : If True and S != NULL, source is modified on first iteration
  S : source vector (or NULL)
  nSmooth : number of smoothing iterations
  U : state vector

OUTPUTS:

  U : modified state vector

RETURN:

  Error code

*/
  int ierr;
  int nIterNonlinear;
  int iIterNonlinear;
  char Verbosity[xf_MAXSTRLEN];
  real CFL;
  real CFLIncreaseFactor;

  if (nSmooth == 0) return xf_OK; // nothing to do

  // pull off several parameters that do not want to change  
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", &nIterNonlinear));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", &iIterNonlinear));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFL));
  if (ierr != xf_OK) return ierr;
  /*  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFLIncreaseFactor", 
      &CFLIncreaseFactor)); */
  /*   if (ierr != xf_OK) return ierr; */
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "Verbosity", Verbosity));
  if (ierr != xf_OK) return ierr;


  // set nIterNonlinear to nSmooth
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", nSmooth));
  if (ierr != xf_OK) return ierr;
  /* // do not increase CFL during smoothing */
  /*   ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFLIncreaseFactor", 1.0)); */
  /*   if (ierr != xf_OK) return ierr; */
  // medium verbosity during smoothing
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", "Medium"));
  if (ierr != xf_OK) return ierr;

  // call the smoother
  ierr = xf_Error(xf_SolveNonlinearSystem_Newton(All, c, LinearFlag, ErrorSolve, S, U, NULL));
  if (ierr != xf_OK) return ierr;

  // reset the parameters that we did not want to change  
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", nIterNonlinear));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", iIterNonlinear));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFL));
  if (ierr != xf_OK) return ierr;
  /* ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFLIncreaseFactor", 
     CFLIncreaseFactor)); */
  /*   if (ierr != xf_OK) return ierr; */
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", Verbosity));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_MGCoarseSolve
static int
xf_MGCoarseSolve(xf_All *All, real c, enum xfe_Bool LinearFlag,
		 xf_Vector *S, xf_Vector *U)
{ 
/*
PURPOSE:

  Wrapper for a "smoothing" call to a nonlinear solver.

INPUTS:

  All: All structure
  c : constant in front of M*U product
  LinearFlag : If True, the Jacobian matrix will not be recalculated if
               it already exists.
  S : source vector (or NULL)
  U : state vector

OUTPUTS:

  U : modified state vector

RETURN:

  Error code

*/
  int ierr;
  int nIterNonlinear;
  int iIterNonlinear;
  char nCoarseIter[xf_MAXSTRLEN];
  char CoarseLinearSolver[xf_MAXSTRLEN], CoarsePreconditioner[xf_MAXSTRLEN];
  char LinearSolver[xf_MAXSTRLEN], Preconditioner[xf_MAXSTRLEN];
  char Verbosity[xf_MAXSTRLEN];
  real CFL;

  // pull off coarse-level settings
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "PMG_nCoarseIter", nCoarseIter));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "PMG_CoarseLinearSolver", 
				 CoarseLinearSolver));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "PMG_CoarsePreconditioner", 
				 CoarsePreconditioner));
  if (ierr != xf_OK) return ierr;

  // pull off several parameters that do not want to change  
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", &nIterNonlinear));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", &iIterNonlinear));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "Verbosity", Verbosity));
  if (ierr != xf_OK) return ierr; 
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "LinearSolver", LinearSolver));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "Preconditioner", Preconditioner));
  if (ierr != xf_OK) return ierr;


  // set nIterNonlinear to nCoarseIter
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "nIterNonlinear", nCoarseIter));
  if (ierr != xf_OK) return ierr;
  // set LinearSolver
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "LinearSolver", CoarseLinearSolver));
  if (ierr != xf_OK) return ierr;
  // set preconditioner
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Preconditioner", CoarsePreconditioner));
  if (ierr != xf_OK) return ierr;
  // medium verbosity during coarse solve
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", "Medium"));
  if (ierr != xf_OK) return ierr;

  // call the solver
  ierr = xf_Error(xf_SolveNonlinearSystem_Newton(All, c, LinearFlag, xfe_True, S, U, NULL));
  if (ierr != xf_OK) return ierr;

  // reset the parameters that we did not want to change  
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", nIterNonlinear));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", iIterNonlinear));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", Verbosity));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "LinearSolver", LinearSolver));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Preconditioner", Preconditioner));
  if (ierr != xf_OK) return ierr;


  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_MGResidual
static int
xf_MGResidual(xf_All *All, real c, xf_Vector *S, xf_Vector *U, 
	      xf_Vector **pR, enum xfe_Bool *pRewind)
{ 
/*
PURPOSE:

  Calculates the MG residual:

               (*pR) = c*M*U + R(U) + S

INPUTS:

  All: All structure
  c : constant in front of M*U product
  S : source vector (or NULL)
  U : state vector

OUTPUTS:

  (*pR) : MG residual
  (*pRewind) : should we rewind the calculation (i.e. error occured?) 
               (optional)

RETURN:

  Error code

*/
  int ierr;

  // locate a residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, xfe_True, 
				       NULL, pR, NULL));
  if (ierr != xf_OK) return ierr;

  // first, R(U)
  ierr = xf_CalculateResidual(All, U, (*pR), NULL, NULL);
  if (pRewind != NULL){
    if (xf_CheckSolverError(ierr, xfe_True, NULL, NULL, NULL, pRewind)) 
      return xf_Error(xf_SOLVER_ERROR);
    if (*pRewind) return xf_OK; // return if solution needs rewind
  }
  if (ierr != xf_OK) return ierr;

  // next, add source if present
  if (S != NULL){
    ierr = xf_Error(xf_SetVector(S, xfe_Add, (*pR)));
    if (ierr != xf_OK) return ierr;
  }

  // finally, add c*M*U to residual
  ierr = xf_Error(xf_AddMassMatrix(All, c, NULL, U, (*pR), NULL, NULL));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SolveNonlinearSystem_pMultigrid
static int
xf_SolveNonlinearSystem_pMultigrid(xf_All *All, real c, enum xfe_Bool LinearFlag,
				   xf_Vector *S, xf_Vector *U)
{ 
/*
PURPOSE:

  Solves:

                   c*M*U + R(U) + S = 0

  using a p-Multigrid method.  The smoother on each level consists of
  a call to the Newton solver, which encompasses single preconditioned
  iterative steps (e.g. block and line).

INPUTS:

  All: All structure
  c : constant in front of M*U product
  LinearFlag : If True, the Jacobian matrix will not be recalculated if
               it already exists.
  S : source vector (or NULL)
  U : state vector

OUTPUTS:

  U : modified state vector

RETURN:

  Error code

*/
  int ierr;
  int nPreSmooth, nPostSmooth;
  int iLevel, nLevel;
  int iMG, nMG, iIter;
  int egrp, negrp;
  int WriteInterval;
  int **Orders = NULL;
  int *CoarseOrders = NULL;
  char CoarseOrdersString[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  enum xfe_Bool Rewind, UpdateFlag, Converged = xfe_False;
  enum xfe_Bool LimitFlag = xfe_False, LimitFlagTot = xfe_False;
  enum xfe_Verbosity Verbosity;
  enum xfe_MultigridCycleType MGCycle;
  real ResidualTolerance;
  xf_Vector *R = NULL, *USafe = NULL;
  xf_Vector **Ui = NULL;
  xf_Vector **Si = NULL;
  xf_SolverData *SolverData = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  if (LinearFlag){
    xf_printf("LinearFlag == True not supported in p-Multigrid.\n");
    xf_printf("Continuing with LinearFlag = False.\n");
    LinearFlag = xfe_False;
  }

  /*** Pull off parameters ***/

  // multigrid cycle type
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "PMG_Cycle", 
				     xfe_MultigridCycleName, (int ) xfe_MultigridCycleLast, 
				     (int *) &MGCycle));
  if (ierr != xf_OK) return ierr;

  // number of levels
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "PMG_nLevel", &nLevel));
  if (ierr != xf_OK) return ierr;

  // list of coarse orders
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "PMG_CoarseOrders", 
				 CoarseOrdersString));
  if (ierr != xf_OK) return ierr;

  // current nonlinear iteration (i.e. MG iteration)
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", &iIter);
  if (ierr != xf_OK) return ierr;

  // number of MG iterations = number of nonlinear iterations
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", &nMG));
  if (ierr != xf_OK) return ierr;

  // Steady write interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "SteadyWriteInterval", 
				    &WriteInterval));
  if (ierr != xf_OK) return ierr;
  // set to -1 so that nonlinear solver does not write too often
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "SteadyWriteInterval", 
				    -1));
  if (ierr != xf_OK) return ierr;
  // SavePrefix
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;

  // number of pre-smoothing iterations
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "PMG_nPreSmooth", &nPreSmooth));
  if (ierr != xf_OK) return ierr;

  // number of post-smoothing iterations
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "PMG_nPostSmooth", &nPostSmooth));
  if (ierr != xf_OK) return ierr;

  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // pull off residual tolerance
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "ResidualTolerance", &ResidualTolerance);
  if (ierr != xf_OK) return ierr;

  
  /*** Construct orders for each element group on each level ***/

  // Variable Order is not yet supported
  if (U->vOrder != NULL) return xf_Error(xf_NOT_SUPPORTED);

  negrp = Mesh->nElemGroup;

  if (nLevel > 0){  // Use nLevel
    if (xf_NotNull(CoarseOrdersString))
      xf_printf("Warning, both nLevel and CoarseOrders specified.  Using nLevel.\n");
    
    ierr = xf_Error(xf_Alloc2((void ***) &Orders, nLevel, negrp, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    // fine orders from U
    for (egrp=0; egrp<negrp; egrp++) Orders[0][egrp] = U->Order[egrp];

    // coarse orders obtained by stepping down by 1
    for (iLevel=1; iLevel<nLevel; iLevel++)
      for (egrp=0; egrp<negrp; egrp++) 
	Orders[iLevel][egrp] = max(Orders[iLevel-1][egrp]-1, 0);
    
  }                 // Use CoarseOrders
  else if (xf_NotNull(CoarseOrdersString)){
    ierr = xf_Error(xf_ScanXIntAlloc( CoarseOrdersString, &nLevel, &CoarseOrders));
    if (ierr != xf_OK) return ierr;
    nLevel = nLevel+1; // to account for fine level

    ierr = xf_Error(xf_Alloc2((void ***) &Orders, nLevel, negrp, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    // fine orders from U
    for (egrp=0; egrp<negrp; egrp++) Orders[0][egrp] = U->Order[egrp];

    // coarse orders set directly
    for (iLevel=1; iLevel<nLevel; iLevel++){
      if (CoarseOrders[iLevel-1] < 0) return xf_Error(xf_INPUT_ERROR);
      for (egrp=0; egrp<negrp; egrp++) 
	Orders[iLevel][egrp] = CoarseOrders[iLevel-1];
    }
  }
  else{
    xf_printf("Need to specify nLevel or CoarseOrders for p-Multigrid.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  if (Verbosity != xfe_VerbosityLow)
    for (iLevel=0; iLevel<nLevel; iLevel++)
      for (egrp=0; egrp<negrp; egrp++) 
	xf_printf("iLevel = %d, egrp = %d, Order = %d\n", iLevel, egrp, Orders[iLevel][egrp]);
  


  /*** Allocate vectors on coarse levels ***/
 
  // pointer to state vectors on all levels
  ierr = xf_Error(xf_Alloc( (void **) &Ui, nLevel, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  // pointer to source vectors on all levels
  ierr = xf_Error(xf_Alloc( (void **) &Si, nLevel, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  
  // set fine level pointers
  Ui[0] = U;
  Si[0] = S;

  // allocate other vectors on coarse levels, and zero them out
  for (iLevel=1; iLevel<nLevel; iLevel++){
    sprintf(Title, "State_Level%d", iLevel);

    // create a shell of a vector: state
    ierr = xf_Error(xf_CreateVector(Ui+iLevel));
    if (ierr != xf_OK) return ierr;

    // Copy state vector from U to iLevel
    ierr = xf_CopyVector(All->Mesh, U, Ui[iLevel]);
    if (ierr != xf_OK) return ierr;

    // project state to required order
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderSet(All->Mesh, NULL, Ui[iLevel], U->Basis, 
						     xfe_BasisLast, NULL, Orders[iLevel], -1));
    if (ierr != xf_OK) return ierr;

    // create a shell of a vector: source
    ierr = xf_Error(xf_CreateVector(Si+iLevel));
    if (ierr != xf_OK) return ierr;
    
    // copy projected state to create a coarse source vector
    ierr = xf_CopyVector(All->Mesh, Ui[iLevel], Si[iLevel]);
    if (ierr != xf_OK) return ierr;

  } // iLevel


  /*** MG Iterations ***/

  if (MGCycle != xfe_MultigridCycleVCycle) return xf_Error(xf_NOT_SUPPORTED);

  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;

  // pull off CFL, CFLIncreaseFactor, CFLDecreaseFactor, etc.
  ierr = xf_Error(xf_FindCFLData(All->Param->KeyValue, SolverData));
  if (ierr != xf_OK) return ierr;
  
  // locate "safe" U vector
  ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], "USafe_PMG", xfe_False,
				       xfe_True, NULL, &USafe, NULL));
  if (ierr != xf_OK) return ierr;

  // Set Safe quantities
  Rewind = xfe_False;
  ierr = xf_Error(xf_SetVector(Ui[0], xfe_Set, USafe));
  if (ierr != xf_OK) return ierr;
  SolverData->CFLSafe = SolverData->CFL;

  // VCycles
  for (iMG=0; iMG<nMG; iMG++, iIter++){

    // Check if a Rewind was requested
    if (Rewind){
      // Restore Ui[0] to USafe
      xf_printf("Restoring U to USafe\n");
      ierr = xf_Error(xf_SetVector(USafe, xfe_Set, Ui[0]));
      if (ierr != xf_OK) return ierr;
      // Decrease CFL
      SolverData->CFL = min(SolverData->CFLSafe, 
			    SolverData->CFL/SolverData->CFLDecreaseFactor);
      xf_printf("Rewind: decreasing CFL, now  = %.10E, min = %.10E\n", 
		SolverData->CFL, SolverData->CFLMin);
      if (SolverData->CFL < SolverData->CFLMin){
	xf_printf("CFL is below CFLMin.  Exiting.\n");
	return xf_Error(xf_SOLVER_ERROR); // show's over
      }
    }

    LimitFlagTot = xfe_False;
    Rewind = xfe_False;

    // write log entry (also print to stdout)
    sprintf(Title, "%% -- Multigrid iteration %d --\n", iMG);
    ierr = xf_Error(xf_WriteLogLine(All, Title));
    if (ierr != xf_OK) return ierr;


    // Downward leg
    for (iLevel=0; iLevel<(nLevel-1); iLevel++){

      if (Verbosity != xfe_VerbosityLow) xf_printf("Down: iLevel = %d\n", iLevel);

      // pre-smooth on iLevel (ErrorSolve if not on finest level)
      ierr = xf_Error(xf_MGSmooth(All, c, LinearFlag, (iLevel != 0), 
				  Si[iLevel], nPreSmooth, Ui[iLevel]));
      if (ierr != xf_OK) return ierr;

      // calculate MG residual (c*M*U + R(U) + S) on iLevel
      ierr = xf_Error(xf_MGResidual(All, c, Si[iLevel], Ui[iLevel], &R, &Rewind));
      if (ierr != xf_OK) return ierr;
      if (Rewind) break;

      // compute residual norm for convergence check
      if (iLevel == 0){
	ierr = xf_Error(xf_VectorNorm(R, 1, &SolverData->ResNorm));
	if (ierr != xf_OK) return ierr;
	if (SolverData->ResNorm < ResidualTolerance){
	  Converged = xfe_True;
	  break;
	}
      }

      // restrict residual to source of iLevel+1
      ierr = xf_Error(xf_ProjectVector(All, R, xfe_True, Si[iLevel+1]));
      if (ierr != xf_OK) return ierr;

      // restrict state to iLevel+1
      ierr = xf_Error(xf_ProjectVector(All, Ui[iLevel], xfe_False, Ui[iLevel+1]));
      if (ierr != xf_OK) return ierr;

      /* Note: Incorporation of residual at restricted state is taken
	 care of by ErrorSolve */
      
    } // iLevel

    if (Rewind) continue;
    if (Converged) break;

    // solve/smooth coarse level problem    
    if (Verbosity != xfe_VerbosityLow) xf_printf("On coarsest level\n");
    ierr = xf_Error(xf_MGCoarseSolve(All, c, LinearFlag, Si[nLevel-1], Ui[nLevel-1]));
    if (ierr != xf_OK) return ierr;


    // Upward leg
    for (iLevel=nLevel-2; iLevel>=0; iLevel--){

      if (Verbosity != xfe_VerbosityLow) xf_printf("Up: iLevel = %d\n", iLevel);

      // -- prolongate state correction from iLevel+1 to iLevel --
      // first, find a vector on iLevel
      ierr = xf_Error(xf_FindSimilarVector(All, Ui[iLevel], "Residual",
					   xfe_False, xfe_True, NULL, &R, NULL));
      if (ierr != xf_OK) return ierr;
      // next, project Ui[iLevel] to Si[iLevel+1]
      ierr = xf_Error(xf_ProjectVector(All, Ui[iLevel], xfe_False, Si[iLevel+1]));
      if (ierr != xf_OK) return ierr;
      // next, subtract Si[iLevel+1] from Ui[iLevel+1]
      ierr = xf_Error(xf_SetVector(Si[iLevel+1], xfe_Sub, Ui[iLevel+1]));
      if (ierr != xf_OK) return ierr;
      // next, prolongate Ui[iLevel+1] to R
      ierr = xf_Error(xf_ProjectVector(All, Ui[iLevel+1], xfe_False, R));
      if (ierr != xf_OK) return ierr;
      // finally, add R(=dU) to Ui[iLevel].  but ...
      // (U + dU) may be non-physical; check what fraction of dU we can add
      ierr = xf_Error(xf_UpdateState(All, Ui[iLevel], R, &LimitFlag,
				     &UpdateFlag, NULL));
      if (ierr != xf_OK) return ierr;
      LimitFlagTot = (LimitFlagTot || LimitFlag);
      Rewind       = (UpdateFlag == xfe_False);
      if (Rewind) break;
      
      // post-smooth on iLevel
      ierr = xf_Error(xf_MGSmooth(All, c, LinearFlag, xfe_False, Si[iLevel], 
				  nPostSmooth, Ui[iLevel]));
      if (ierr != xf_OK) return ierr;

    } // iLevel

    // adjust CFL and store USafe if update was not limited
    if ((!Rewind) && (!LimitFlagTot)){
      SolverData->CFLSafe = SolverData->CFL;
      ierr = xf_Error(xf_SetVector(Ui[0], xfe_Set, USafe));
      if (ierr != xf_OK) return ierr;
      SolverData->CFL *= SolverData->CFLIncreaseFactor;
      SolverData->CFL = min(SolverData->CFL, SolverData->CFLMax);
    }

    // store parameters
    ierr = xf_SetKeyValueReal(All->Param->KeyValue, "CFL", SolverData->CFL);
    if (ierr != xf_OK) return ierr;

    // set iIterNonlinear
    ierr = xf_SetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", iIter);
    if (ierr != xf_OK) return ierr;

    /* Write out U to hard disk if at requested interval */
    if ((iMG > 0) && (WriteInterval > 0) && (iMG % WriteInterval) == 0){
      sprintf(OutputFile, "%s_MGIter%d.data\0", SavePrefix, iMG);
      ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "State", U, OutputFile));
      if (ierr != xf_OK) return ierr;
    }
    
  } // iMG

  if ((Converged) && (Verbosity != xfe_VerbosityLow))
    xf_printf("p-Multigrid nonlinear solver converged to tolerance.\n");
 
  // reset steady write interval
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "SteadyWriteInterval", 
				    WriteInterval));
  if (ierr != xf_OK) return ierr;

  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;

  /* Destroy coarse vectors */
  for (iLevel=1; iLevel<nLevel; iLevel++){
    ierr = xf_Error(xf_DestroyVector(Ui[iLevel], xfe_True));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyVector(Si[iLevel], xfe_True));
    if (ierr != xf_OK) return ierr;
  } // iLevel

  /*** Release memory ***/
  xf_Release2( (void **) Orders      );
  xf_Release(  (void  *) CoarseOrders);
  xf_Release(  (void  *) Ui);
  xf_Release(  (void  *) Si);

  return xf_OK;
}
