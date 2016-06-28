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
 FILE:  xf_SteadyContinue.c
 
 Program used for performing equation-set continuation in steady-state
 runs.  For example, can start with benign IC's parameters, and work
 way to a more demanding parameter choice.
 
 */

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_All.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Param.h"
#include "xf_Solver.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_ParamDefault.h"
#include "xf_EqnSetHook.h"
#include "xf_Arg.h"
#include "xf_Residual.h"
#include <time.h>  // for timing execution

/* Enumerated spacing types */
enum xfe_SpacingType{
  xfe_SpacingLinear,
  xfe_SpacingLog,
  xfe_SpacingLast
};

/* corresponding names */
static char *xfe_SpacingName[xfe_SpacingLast] = {
  "Linear",
  "Log",
};


/******************************************************************/
//   FUNCTION Definition: xf_SetEqnSet
static int
xf_SetEqnSet( xf_EqnSet *StartEqnSet, xf_EqnSet *EndEqnSet,
	      real s, enum xfe_SpacingType Spacing, xf_EqnSet *EqnSet)
{
  int ierr, i;
  char *Key, *StartValue;
  char EndValue[xf_MAXSTRLEN];
  xf_KeyValue StartKeyValue, EndKeyValue;
  real rStart, rEnd, rMid;

  // parameters are stored as key-value pairs
  StartKeyValue = StartEqnSet->KeyValue;
  EndKeyValue   = EndEqnSet->KeyValue;
  
  xf_printf("\nSetting new eqnset parameters\n", Key, rMid);

  // loop over key-value pairs
  for (i=0; i<StartKeyValue.nKey; i++){
    Key = StartKeyValue.Key[i];
    StartValue = StartKeyValue.Value[i];
    ierr = xf_Error(xf_GetKeyValue(EndKeyValue, Key, EndValue));
    if (ierr != xf_OK) return ierr;
    if (strcmp(StartValue, EndValue) != 0){ // values are different
      // well, they should be real, since that's the only difference we can handle
      ierr = xf_Error(xf_GetKeyValueReal(StartKeyValue, Key, &rStart));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_GetKeyValueReal(EndKeyValue, Key, &rEnd));
      if (ierr != xf_OK) return ierr;
      // interpolate between start and end
      if (Spacing == xfe_SpacingLinear)
	rMid = (1.-s)*rStart + s*rEnd;
      else
	rMid = exp((1.-s)*log(rStart) + s*log(rEnd));
      // set value
      xf_printf("Setting %s to %.8E\n", Key, rMid);
      ierr = xf_Error(xf_SetKeyValueReal(EqnSet->KeyValue, Key, rMid));
      if (ierr != xf_OK) return ierr;
    }
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: main
int
main(int argc, char *argv[])
{

  int ierr;
  int myRank, nProc;
  int nfail, nfailmax, niter;
  char *ArgIn[] = {"job", "NULL", ".job file describing the case",
		   "starteqn", "NULL", "starting .eqn file",
		   "endeqn", "NULL", "ending .eqn file",
		   "spacing", "Log", "continuation spacing; Linear or Log",
		   "nspace", "8", "nominal number of runs after first one",
		   "dsmin", "1e-3", "s in [0,1]; min ds below which give up",
		   "nfailmax", "10", "max # of failed solves before give up",
		   "strictrestol", "True", "res tol not achieved -> run failed",
		   "\0"};
  char jobFile[xf_MAXSTRLEN];
  char StartEqnSetFile[xf_MAXSTRLEN];
  char EndEqnSetFile[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];

  real s, ds, ds0, dsmin;
  real ResidualTolerance;
  real ResNorm, CFLStart;

  int nSpace;
  enum xfe_SpacingType Spacing;
  
  enum xfe_Bool RestartFlag;
  enum xfe_Bool StrictResTol;
  enum xfe_Bool done, success, Fail;

  clock_t clock_Start, clock_End;

  xf_Vector *U=NULL, *USafe=NULL, *R=NULL;
  xf_SolverData *SolverData;

  xf_EqnSet *StartEqnSet = NULL;
  xf_EqnSet   *EndEqnSet = NULL;
  xf_KeyValue KeyValue;
  xf_All *All; 
  
  
  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  xf_printf("\n");
  xf_printf("=== xf_SteadyContinue: Parameter continuation solver wrapper ===\n");
  xf_printf("\n");
          

  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  // Get jobFile
  ierr = xf_GetKeyValue(KeyValue, "job", jobFile);
  if (ierr != xf_OK) return ierr;

  // Get StartEqnSetFile 
  ierr = xf_GetKeyValue(KeyValue, "starteqn", StartEqnSetFile);
  if (ierr != xf_OK) return ierr;

  // Get EndEqnSetFile 
  ierr = xf_GetKeyValue(KeyValue, "endeqn", EndEqnSetFile);
  if (ierr != xf_OK) return ierr;

  // Get Spacing enumerated type
  ierr = xf_Error(xf_GetKeyValueEnum(KeyValue, "spacing", xfe_SpacingName, 
				     xfe_SpacingLast, (int *) &Spacing));
  if (ierr != xf_OK) return ierr;
  
  // Get nSpace
  ierr = xf_GetKeyValueInt(KeyValue, "nspace", &nSpace);
  if (ierr != xf_OK) return ierr;

  // Get nfailmax
  ierr = xf_GetKeyValueInt(KeyValue, "nfailmax", &nfailmax);
  if (ierr != xf_OK) return ierr;

  // Get dsmin
  ierr = xf_GetKeyValueReal(KeyValue, "dsmin", &dsmin);
  if (ierr != xf_OK) return ierr;

  // Get StrictResTol
  ierr = xf_GetKeyValueBool(KeyValue, "strictrestol", &StrictResTol);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  /*-------------------------------*/
  /* Read in .job file into All    */
  /*-------------------------------*/

  ierr = xf_Error(xf_ReadAllFromJobFile(jobFile, xfe_False, &All));
  if (ierr != xf_OK) return ierr;
  

  /*----------------------------------*/
  /* Read in start/end equation sets  */
  /*----------------------------------*/

  ierr = xf_Error(xf_CreateEqnSet(&StartEqnSet));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReadEqnSetFile(StartEqnSetFile, NULL, StartEqnSet));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_CreateEqnSet(&EndEqnSet));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReadEqnSetFile(EndEqnSetFile, NULL, EndEqnSet));
  if (ierr != xf_OK) return ierr;
  
  /* Also read StartEqnSet into All->EqnSet */
  ierr = xf_Error(xf_ReadEqnSetFile(StartEqnSetFile, NULL, All->EqnSet));
  if (ierr != xf_OK) return ierr;

  
  /*----------------------*/
  /* Load EqnSet Library  */
  /*----------------------*/

  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;
  All->EqnSet->Dim = All->Mesh->Dim;


  /*-----------------------*/
  /* Register equation set */
  /*-----------------------*/
  
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;


  /*---------------------*/
  /* Get params from All */
  /*---------------------*/
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "Restart", &RestartFlag));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;

  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "ResidualTolerance", &ResidualTolerance);
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFLStart));
  if (ierr != xf_OK) return ierr;

  
  /*-----------------------*/
  /* Create a state vector */
  /*-----------------------*/
  
  ierr = xf_Error(xf_FindOrCreatePrimalState(All, RestartFlag, NULL, &U));
  if (ierr != xf_OK) return ierr;  


  /*--------------------------------------*/
  /* Call nonlinear solver at StartEqnSet */
  /*--------------------------------------*/
  
  xf_printf("\n *** Initial solve, no continuation yet ***\n\n");
  clock_Start = clock(); // start timer          
  ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &U));
  if (ierr != xf_OK){
    xf_printf("Error occured during initial solve!  Cannot proceed with continuation.\n");
    return ierr;
  }
  clock_End = clock(); // end timer
  xf_printf("Initial steady solve CPU time = %.10E\n", 
	    ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC);

  // locate "safe" U vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "ContinuationUSafe", xfe_False, xfe_True, 
				       NULL, &USafe, NULL));
  if (ierr != xf_OK) return ierr;
  
  // locate Residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, xfe_True, 
				       NULL, &R, NULL));
  if (ierr != xf_OK) return ierr;

  /*-------------------*/
  /* Create SolverData */
  /*-------------------*/

  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;


  /*------------*/
  /* Begin Loop */
  /*------------*/

  s = 0.;                // s=0 means at start, s=1 means at end
  ds = ds0 = 1./nSpace;  // delta s
  
  done    = xfe_False;
  success = xfe_False;
  nfail   = 0;
  niter   = 0;

  while (!done){

    // save U in USafe
    ierr = xf_Error(xf_SetVector(U, xfe_Set, USafe));
    if (ierr != xf_OK) return ierr;

    // increment s towards end by ds
    s = s + ds;

    // if s is close to 1 or past it, just set it to 1 and prepare done flag
    if ((1.-s) < 0.1*ds){
      s = 1.0;
      done = xfe_True;
    }

    // set EqnSet using s
    ierr = xf_Error(xf_SetEqnSet(StartEqnSet, EndEqnSet, s, Spacing, All->EqnSet));
    if (ierr != xf_OK) return ierr;

    /* Set CFL */
    ierr = xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart);
    if (ierr != xf_OK) return ierr;


    // try solve
    xf_printf("\n *** Continuation solve at s=%.4f, ds=%.4f ***\n", s, ds);
    clock_Start = clock(); // start timer
    ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &U));

    Fail = xfe_False;

    // break out if user requests a halt
    if (xf_CheckUserHalt(NULL)) break;

    // Did solve fail?
    if (ierr != xf_OK){
      Fail = xfe_True; // yes if it threw an error
    }
    else if (StrictResTol){
      // also, residual not below tolerance can mean fail if user deems it so
      ierr = xf_CalculateResidual(All, U, R, NULL, SolverData);
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorNorm(R, 1, &ResNorm));
      if (ierr != xf_OK) return ierr;
      if (ResNorm > ResidualTolerance) Fail = xfe_True;
    }
    
    if (Fail){
      // solve not successful: lower ds, recall U from USafe, continue
      nfail++;
      s = s - ds;
      ds = ds/2.0;
      ierr = xf_Error(xf_SetVector(USafe, xfe_Set, U));
      if (ierr != xf_OK) return ierr;

      xf_printf("Solve failed. Lowering ds to %.4f.\n", ds);
      done = xfe_False;

      // do not allow ds to drop too small or nfail to get too big
      if ((ds < dsmin) || (nfail > nfailmax)){
	xf_printf("ds=%.3f, nfail=%d -> giving up\n", ds, nfail);
	xf_printf("Will exit with last successful state at s=%.3f\n", s);
	done = xfe_True;
	success = xfe_False;
      }

      continue;
    }
    clock_End = clock(); // end timer
    xf_printf("Continuation  solve CPU time = %.10E\n", 
	      ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC);
    
    // if successful, set ds = ds0 and keep going
    ds = ds0;
    if (done) success = xfe_True;
    else{
      // write out intermediate file
      sprintf(OutputFile, "%s_iter%d_s%.3f.xfa\0", SavePrefix, niter, s);
      ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
      if (ierr!=xf_OK) return ierr;
      niter++;
    }

    xf_printf("Solve succeeded! ds is %.4f.\n", ds);

  } // end while


  /*----------------------*/
  /* Write out final file */
  /*----------------------*/

  if (success)
    sprintf(OutputFile, "%s_Success.xfa\0", SavePrefix);
  else
    sprintf(OutputFile, "%s_Fail.xfa\0", SavePrefix);
  xf_printf("Writing %s\n", SavePrefix);
  ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
  if (ierr!=xf_OK) return ierr;


  /*----------*/
  /* Clean up */
  /*----------*/

  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;

  xf_printf("xf_SteadyContinue finished.\n");


  /* MPI finalize */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return xf_OK;

}
