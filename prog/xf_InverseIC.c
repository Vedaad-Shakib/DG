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
  FILE:  xf_InverseIC.c

  This program performs Hessian-based model reduction of initial
  conditions.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_MeshTools.h"
#include "xf_Param.h"
#include "xf_Basis.h"
#include "xf_EqnSetHook.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_Residual.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Output.h"
#include "xf_Solver.h"
#include "xf_LinearSolver.h"
#include "xf_Arg.h"
#include "xf_InverseIC_Common.h"
#include <time.h>


/******************************************************************/
//   FUNCTION Definition: xf_ApplyHessian
static int 
xf_ApplyHessian(xf_All *All, xf_JacobianMatrix *R_U, int nOut, xf_Vector **J_U0,
		xf_Vector *V, xf_Vector *Utemp, 
		xf_TimeHistData *TimeHistData, int nOutput, char *SumOutputName,
		int *TimeFlag, enum xfe_Bool MReg, real betaM, enum xfe_Bool LReg, real betaL,
		real beta, xf_SolverData *SolverData, enum xfe_Bool *pBreakFlag, xf_Vector *W,
		real *J)
{
  int ierr, i, iOut;
  int nTime, iTime;
  real val;

  (*pBreakFlag) = xfe_False;

  /*---------------------------------*/
  /* Set W = "operator applied to V" */
  /*---------------------------------*/

  if (J_U0 != NULL){

    // J(j) = J_U0(j) * V
    for (iOut=0; iOut<nOut; iOut++){
      ierr = xf_Error(xf_VectorDot(J_U0[iOut], V, J+iOut));
      if (ierr != xf_OK) return ierr;
    } // iOut

    // W = J_U0(j) * J(j)
    for (iOut=0; iOut<nOut; iOut++){
      ierr = xf_Error(xf_VectorMultSet(J_U0[iOut], J[iOut], 
				       ((iOut==0) ? xfe_Set : xfe_Add), W));
      if (ierr != xf_OK) return ierr;
    }

  }
  else{


    // First copy V into Utemp so that we do not overwrite V (or U)
    ierr = xf_Error(xf_SetVector(V, xfe_Set, Utemp));
    if (ierr != xf_OK) return ierr;


    // Unsteady Forward Solve
    ierr = xf_Error(xf_ApplyTimeScheme(All, NULL, xfe_False, &Utemp, TimeHistData));
    if (ierr != xf_OK) return ierr;

    // if Halted, break
    if (xf_CheckUserHalt(NULL)){ (*pBreakFlag) = xfe_True; return xf_OK;}

    nTime = TimeHistData->nTime;
    
    // Weight time-history outputs for adjoint solve
    TimeHistData->nSumWeight = nOutput;
    ierr = xf_Error(xf_ReAlloc2((void ***) &TimeHistData->SumOutputWeights,
				nTime, nOutput, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (iTime=1; iTime<nTime; iTime++){
      val = (TimeFlag == NULL) ? 1.0 : ((real) TimeFlag[iTime]);
      for (i=0; i<nOutput; i++)
	TimeHistData->SumOutputWeights[iTime][i] = val*TimeHistData->OutputValues[i][iTime];
    }
    
    // Set TimeWeights to all 1s
    ierr = xf_Error(xf_ReAlloc((void **) &TimeHistData->TimeWeights,
			       nTime, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (iTime=0; iTime<nTime; iTime++) TimeHistData->TimeWeights[iTime] = 1.0;


    // make W look like an adjoint vector associated with SumOutput
    if (W->OutputName != NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
    W->OutputName = SumOutputName; 

    // Initialize W to zero
    ierr = xf_Error(xf_SetZeroVector(W));
    if (ierr != xf_OK) return ierr;

    // Unsteady Adjoint Solve
    ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, NULL, Utemp, 1, &W, TimeHistData));
    if (ierr != xf_OK) return ierr;
    W->OutputName = NULL; // back to normal
  }


  /** Add regularization **/

  if (MReg){
    /* Mass matrix to approximate a continuous L2 norm */
    //Utemp = V
    ierr = xf_Error(xf_SetVector(V, xfe_Set, Utemp));
    if (ierr != xf_OK) return ierr;
    // Utemp = M*Utemp
    ierr = xf_Error(xf_MultMassMatrix(All, 1.0, Utemp));
    if (ierr != xf_OK) return ierr;
    // W += betaM*Utemp
    ierr = xf_Error(xf_VectorMultSet(Utemp, betaM, xfe_Add, W));
    if (ierr != xf_OK) return ierr;
  }
    
  if (LReg){

    /* Laplace regularization */

    /* Utemp = R_U * V  */
    ierr = xf_Error(xf_Jacobian_Mult(All, R_U, V, xfe_Set, xfe_False, SolverData, Utemp));
    if (ierr != xf_OK) return ierr;

    /* W += betaL*Utemp  */
    ierr = xf_Error(xf_VectorMultSet(Utemp, betaL, xfe_Add, W));
    if (ierr != xf_OK) return ierr;

    /* Utemp = (R_U)^T * V  */
    ierr = xf_Error(xf_Jacobian_Mult(All, R_U, V, xfe_Set, xfe_True, SolverData, Utemp));
    if (ierr != xf_OK) return ierr;

    /* W += betaL*Utemp  */
    ierr = xf_Error(xf_VectorMultSet(Utemp, betaL, xfe_Add, W));
    if (ierr != xf_OK) return ierr;
  }
    
  if ((!MReg) && (!LReg)){
    // A simple discrete L2 norm: W += beta*V
    ierr = xf_Error(xf_VectorMultSet(V, beta, xfe_Add, W));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int len, ierr, i, j, outerr;
  int myRank, nProc;
  int nOutput = 0;
  int iIter, nCG, nTime, iTime, maxiter_CG, nTimeMask;
  enum xfe_Bool RestartFlag, done, MReg, LReg, PreconditionFlag;
  enum xfe_Bool CGRestart, SavePoint, Precond, BreakFlag = xfe_False;
  enum xfe_Bool ParamFlag = xfe_False, SingleOutFlag = xfe_False;
  enum xfe_Bool Inside = xfe_False;
  enum xfe_TimeSchemeType TimeScheme;
  enum xfe_LinearStatusType Status;
  char *ArgIn[] = {"job", "NULL", ".job file name to read (run parameters)",
		   "yobs", "yobs.txt", "name file containing observations",
		   "beta", "1e-6", "regularization coefficient for L2 reg",
		   "betaM", "1e-6", "regularization coefficient for Mass reg",
		   "betaL", "1e-6", "regularization coefficient for Laplace reg",
		   "tol", "1e-8", "tolerance for linear solution",
		   "sol", "UIC.data", "file where inverse solution is written",
		   "MReg", "False", "regularization using Mass matrix (and precond)",
		   "LReg", "False", "regularization using Laplace matrix",
		   "maxiter", "200", "maximum number of CG iterations",
		   "Restart", "False", "if True, restart from CG.data and CG.txt",
		   "SavePoint", "True", "if True, writes CG.data/txt on every iter",
		   "Precond", "True", "if True, preconditioned CG is used",
		   "TimeMask", "NULL", "time mask file for fewer outputs",
		   "Prefix", "NULL", "prefix for any files written out",
		   "ParamFile", "NULL", "parameters for parametric inversion",
		   "epsFD", "1e-4", "finite difference epsilon (param inversion)",
		   "OutputFile", "NULL", "to use single outputs (param inversion)",
		   "nTry", "1", "number of random restarts (param inversion)",
		   "\0"};
  char jobFile[xf_MAXSTRLEN] = "NULL"; 
  char yobsFile[xf_MAXSTRLEN];
  char solFile[xf_MAXSTRLEN] = "NULL";
  char Prefix[xf_MAXSTRLEN] = "NULL";
  char TimeMaskFile[xf_MAXSTRLEN] = "NULL";
  char ParamFile[xf_MAXSTRLEN] = "NULL";
  char OutputFile[xf_MAXSTRLEN] = "NULL";
  char s[xf_MAXSTRLEN], line[xf_MAXLONGLINELEN], Title[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN] = "NULL";
  char SumOutputName[] = "SumOutputs";
  char *DiffKeyValue[] = {"Unused", "Unused", "\0"};

  int iNewton = 0;
  int iTry, nTry;
  int nNewtonMax = 20; // max number of Newton iterations before restart
  int iParam;
  int nParam;                // number of parameters for MCMC
  int *P = NULL;             // for PLU
  char **ParamName = NULL;   // parameter names
  real *ParamRange = NULL;   // min/max ranges of MCMC parameters
  real *ParamValue = NULL;   // current parameter values
  real *ParamValueMin = NULL;
  real *muR = NULL, *muJ = NULL;
  real *dmu = NULL;
  real *epsParam = NULL;
  real paramOmega = 0.2;     // Newton under-relaxation
  real dParam, epsFD, rnorm, randval;

  int iOut;
  int nOut = 0;             // total number of outputs (sensors and times)
  char **OutName = NULL;    // sensor names for each output reading
  int *OutTimeIndex = NULL; // time index at which sensor reading was taken
  real *OutValue = NULL;    // sensor output readings
  xf_Vector **J_U0 = NULL;  // sensitivity of outputs w.r.t ICs
  real *J = NULL;
  real cost, costmin=1e30;

  clock_t clock_Start, clock_End;

  FILE *fid;
  real beta, tol, **ObservedOutputs = NULL;
  real betaM, betaL, val;
  real ViscosityLReg, ViscosityOrig;
  int *TimeMask = NULL, *TimeFlag = NULL;
  xf_KeyValue KeyValueArg;
  xf_TimeHistData *TimeHistData = NULL;
  xf_LinearSolverData *LinearSolverData = NULL;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Vector *U = NULL, *V = NULL, *W = NULL;
  xf_Vector *LinR = NULL, *Utemp = NULL, *UIC = NULL;
  xf_Vector *R0 = NULL;
  xf_Vector *U_mu = NULL, *R = NULL;
  xf_VectorSet  *CGSet = NULL;
  xf_Output *SumOutput = NULL;
  xf_ResTerms *ResTermsLReg = NULL;
  xf_ResTerms *ResTermsOrig = NULL;
  xf_JacobianMatrix *R_U = NULL;
  xf_SolverData *SolverData = NULL;
  xf_All *All;

  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  xf_printf("\n");
  xf_printf("=== Linear Inverse Solver for Initial Conditions using Hessian ===\n");
  xf_printf("\n");
    
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValueArg));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValueArg);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  /* Get jobFile */
  ierr = xf_GetKeyValue(KeyValueArg, "job", jobFile);
  if (ierr != xf_OK) return ierr;

  /* Get yobsFile */
  ierr = xf_GetKeyValue(KeyValueArg, "yobs", yobsFile);
  if (ierr != xf_OK) return ierr;

  /* Get solFile */
  ierr = xf_GetKeyValue(KeyValueArg, "sol", solFile);
  if (ierr != xf_OK) return ierr;

  /* Get regularizations */
  ierr = xf_GetKeyValueReal(KeyValueArg, "beta", &beta);
  if (ierr != xf_OK) return ierr;
  if (beta <= 0) return xf_Error(xf_INPUT_ERROR);

  ierr = xf_GetKeyValueReal(KeyValueArg, "betaM", &betaM);
  if (ierr != xf_OK) return ierr;
  if (betaM <= 0) return xf_Error(xf_INPUT_ERROR);

  ierr = xf_GetKeyValueReal(KeyValueArg, "betaL", &betaL);
  if (ierr != xf_OK) return ierr;
  if (betaL <= 0) return xf_Error(xf_INPUT_ERROR);

  /* Get tolerance */
  ierr = xf_GetKeyValueReal(KeyValueArg, "tol", &tol);
  if (ierr != xf_OK) return ierr;
  if (tol <= 0) return xf_Error(xf_INPUT_ERROR);

  /* Get MReg */
  ierr = xf_GetKeyValueBool(KeyValueArg, "MReg", &MReg);
  if (ierr != xf_OK) return ierr;

  /* Get LReg */
  ierr = xf_GetKeyValueBool(KeyValueArg, "LReg", &LReg);
  if (ierr != xf_OK) return ierr;

  /* Get maxiter_CG */
  ierr = xf_GetKeyValueInt(KeyValueArg, "maxiter", &maxiter_CG);
  if (ierr != xf_OK) return ierr;

  /* CGRestart */
  ierr = xf_GetKeyValueBool(KeyValueArg, "Restart", &CGRestart);
  if (ierr != xf_OK) return ierr;

  /* SavePoint */
  ierr = xf_GetKeyValueBool(KeyValueArg, "SavePoint", &SavePoint);
  if (ierr != xf_OK) return ierr;

  /* Precond */
  ierr = xf_GetKeyValueBool(KeyValueArg, "Precond", &Precond);
  if (ierr != xf_OK) return ierr;

  /* Get Prefix */
  ierr = xf_GetKeyValue(KeyValueArg, "Prefix", Prefix);
  if (ierr != xf_OK) return ierr;

  /* Get TimeMaskFile */
  ierr = xf_GetKeyValue(KeyValueArg, "TimeMask", TimeMaskFile);
  if (ierr != xf_OK) return ierr;

  /* Get ParamFile */
  ierr = xf_GetKeyValue(KeyValueArg, "ParamFile", ParamFile);
  if (ierr != xf_OK) return ierr;

  /* Get OutputFile */
  ierr = xf_GetKeyValue(KeyValueArg, "OutputFile", OutputFile);
  if (ierr != xf_OK) return ierr;

  /* Get epsFD */
  ierr = xf_GetKeyValueReal(KeyValueArg, "epsFD", &epsFD);
  if (ierr != xf_OK) return ierr;

  /* Get nTry */
  ierr = xf_GetKeyValueInt(KeyValueArg, "nTry", &nTry);
  if (ierr != xf_OK) return ierr;

  // destroy key-value from arg list
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValueArg));
  if (ierr!=xf_OK) return ierr;


  /*-----------------------------------*/
  /*   Read in ParamFile if not NULL   */
  /*-----------------------------------*/

  ParamFlag = xfe_False;
  nParam = 0;
  if (xf_NotNull(ParamFile)){
    ParamFlag = xfe_True; // will be doing parameteric inversion
    ierr = xf_Error(xf_ReadICParamFile(ParamFile, &nParam, &ParamName, &ParamRange, 
				       &ParamValue));
    if (ierr != xf_OK) return ierr;
    
    // allocate memory
    ierr = xf_Error(xf_Alloc((void **) &P, nParam, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &muR, nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &dmu, nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &muJ, nParam*nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &ParamValueMin, nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nParam; i++) ParamValueMin[i] = 0.0;

    // create epsilon for finite differencing of gradients
    ierr = xf_Error(xf_Alloc( (void **) &epsParam, nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (iParam=0; iParam<nParam; iParam++){
      epsParam[iParam] = epsFD*(ParamRange[2*iParam+1] - ParamRange[2*iParam]);
      if (epsParam[iParam] <= 0) return xf_Error(xf_OUT_OF_BOUNDS);
    }
  }


  /*---------------------------------------*/
  /* Read All from .job, including eqnset. */
  /* Parallelization takes place here.     */
  /*---------------------------------------*/

  ierr = xf_Error(xf_ReadAllFromJobFile(jobFile, xfe_True, &All));
  if (ierr != xf_OK) return ierr;
  
  /* Load EqnSet Library  */  
  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;
  
  /* Register EqnSet and check/set default eqnset parameters. */
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;


  /* Set verbosity to low */
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", "Low"));
  if (ierr != xf_OK) return ierr;

  /* Save away SavePrefix */
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;
  /* Set SavePrefix to None so we do not write profusely to .log file */
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "SavePrefix", "None"));
  if (ierr != xf_OK) return ierr;

  /* Perform input error checks */
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
				     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
				     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;
  if ((TimeScheme == xfe_TimeSchemeSteady) || (TimeScheme >= xfe_TimeSchemeLast))
    return xf_Error(xf_INPUT_ERROR);

  if (All->EqnSet->Outputs == NULL) return xf_Error(xf_INPUT_ERROR);

  /*---------------------------------*/
  /* Read in OuptutFile if not NULL  */
  /*---------------------------------*/

  SingleOutFlag = xfe_False;
  nOut = 0;
  if (xf_NotNull(OutputFile) ){
    SingleOutFlag = xfe_True; // will be doing parameteric inversion
    // Read in OutputFile
    ierr = xf_Error(xf_ReadICOutputFile(OutputFile, &nOut, &OutName, &OutTimeIndex, 
					&OutValue, NULL));
    if (ierr != xf_OK) return ierr;

    // create/load sensitivity vectors
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_LoadSensitivityVectors(All, nOut, OutName, OutTimeIndex, &J_U0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
    if (ierr != xf_OK) return ierr;
    
    // allocate output vector
    ierr = xf_Error(xf_Alloc( (void **) &J, nOut, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }


  /*----------------------------------*/
  /* Read in Observation and TimeMask */
  /*----------------------------------*/

  /* Read observation file on root */
  if (!SingleOutFlag){
    outerr = xf_OK;
    if (myRank == 0){

      /* Read in yobs */
      if ((fid = fopen(yobsFile, "r")) == NULL) outerr = xf_Error(xf_FILE_READ_ERROR);
      else{
	/* First line of .out file: nTime (includes t=0), nOutput */
	if (fgets(line, xf_MAXLONGLINELEN, fid) == NULL) outerr = xf_Error(xf_FILE_READ_ERROR);
	if (sscanf(line, "%d %d", &nTime, &nOutput) != 2) outerr = xf_Error(xf_FILE_READ_ERROR);
	else{
	  ierr = xf_Error(xf_Alloc2((void ***) &ObservedOutputs,
				    nTime, nOutput, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	
	  for (iTime=0; (iTime<nTime) && (outerr == xf_OK); iTime++){
	    if (fgets(line, xf_MAXLONGLINELEN, fid) == NULL) outerr = xf_Error(xf_FILE_READ_ERROR);
	    // each line contains nOutput reals
	    ierr = xf_Error(xf_ScanReal(line, nOutput, ObservedOutputs[iTime]));
	    if (ierr != xf_OK) outerr = ierr;
	  }
	}
	fclose(fid);
      }
    } // end if myRank == 0
    if (xf_PError(&outerr, 0) != xf_OK) return outerr;
  

    // Broadcast ObservedOutputs
    ierr = xf_Error(xf_MPI_Bcast( (void *) &nTime,   sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_MPI_Bcast( (void *) &nOutput, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;

    if (myRank > 0){
      ierr = xf_Error(xf_Alloc2((void ***) &ObservedOutputs,
				nTime, nOutput, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }

    ierr = xf_Error(xf_MPI_Bcast((void *) ObservedOutputs[0], 
				 nTime*nOutput*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;


    if (xf_NotNull(TimeMaskFile)){
      /* Read time mask file on root */
      outerr = xf_OK;
      if (myRank == 0){
      
	/* Read in time mask */
	if ((fid = fopen(TimeMaskFile, "r")) == NULL) outerr = xf_Error(xf_FILE_READ_ERROR);
	else{
	  /* First line of .out file: nTimeMask */
	  if (fgets(line, xf_MAXLONGLINELEN, fid) == NULL) outerr = xf_Error(xf_FILE_READ_ERROR);
	  if (sscanf(line, "%d", &nTimeMask) != 1) outerr = xf_Error(xf_FILE_READ_ERROR);
	  else{
	    ierr = xf_Error(xf_Alloc((void **) &TimeMask, nTimeMask, sizeof(int)));
	    if (ierr != xf_OK) return ierr;
	  
	    for (iTime=0; (iTime<nTimeMask) && (outerr == xf_OK); iTime++){
	      if (fgets(line, xf_MAXLONGLINELEN, fid) == NULL) outerr = xf_Error(xf_FILE_READ_ERROR);
	      // each line contains nOutput reals
	      ierr = xf_Error(xf_ScanInt(line, 1, TimeMask+iTime));
	      if (ierr != xf_OK) outerr = ierr;
	    }
	  }
	  fclose(fid);
	}
      } // end if myRank == 0
      if (xf_PError(&outerr, 0) != xf_OK) return outerr;
  

      // Broadcast nTimeMask
      ierr = xf_Error(xf_MPI_Bcast( (void *) &nTimeMask, sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
    
      if (myRank > 0){
	ierr = xf_Error(xf_Alloc((void **) &TimeMask, nTimeMask, sizeof(int)));
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_MPI_Bcast((void *) TimeMask, nTimeMask*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;

      // Create TimeFlag
      ierr = xf_Error(xf_Alloc((void **) &TimeFlag, nTime, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      // Fill TimeFlag
      for (i=0; i<nTime; i++) TimeFlag[i] = 0;
      for (i=0; i<nTimeMask; i++){
	if ((TimeMask[i] < 0) || (TimeMask[i] >= nTime)) return xf_Error(xf_OUT_OF_BOUNDS);
	TimeFlag[TimeMask[i]] = 1;
      }
    }

    // Error checks
    if (nOutput != All->EqnSet->Outputs->nOutput){
      xf_printf("Error, number of outputs in yobs file ! = # in case eqnset.\n");
      return xf_Error(xf_INPUT_ERROR);
    }
    
    if (nOutput <= 0) return xf_Error(xf_INPUT_ERROR);
    
    for (i=0; i<nOutput; i++){
      if (All->EqnSet->Outputs->Output[i].Type == xfe_SumOutput){
	xf_printf("SumOutputs not allowed.\n"); // because we'll be creating one
	return xf_Error(xf_INPUT_ERROR);
      }
    }
    
    // Will need TimeHistData for unsteady runs, with all outputs logged
    line[0] = '\0';
    for (i=0; i<nOutput; i++){
      sprintf(s, "%s ", All->EqnSet->Outputs->Output[i].Name);
      strcat(line, s);
    }
    ierr = xf_Error(xf_CreateUniformTimeHistData(All, line, &TimeHistData));
    if (ierr != xf_OK) return ierr;

    /*---------------------------------*/
    /* Create a sum-output for adjoint */
    /*---------------------------------*/
    
    ierr = xf_Error(xf_ReAllocOutputs(All->EqnSet->Outputs, nOutput+1));
    if (ierr != xf_OK) return ierr;
    
    SumOutput = All->EqnSet->Outputs->Output + nOutput;
    
    // set Name
    ierr = xf_Error(xf_AllocString(&SumOutput->Name, xf_MAXSTRLEN, SumOutputName));
    if (ierr != xf_OK) return ierr;
    
    // set Type + nSumOutput
    SumOutput->Type = xfe_SumOutput;
    SumOutput->nSumOutput= nOutput;
    
    // Allocate and set SumOutputNames
    ierr = xf_Error(xf_Alloc2((void ***) &SumOutput->SumOutputNames, 
			      nOutput, xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    
    for (i=0; i<nOutput; i++)
      strcpy(SumOutput->SumOutputNames[i], All->EqnSet->Outputs->Output[i].Name);
    
    // Allocate SumOutputWeights
    SumOutput->SumOutputWeights     = NULL;
    ierr = xf_Error(xf_Alloc((void **) &SumOutput->SumOutputWeights, 
			     nOutput, sizeof(real)));
    if (ierr != xf_OK) return ierr;  
  }    


  /*-------------------------------------*/
  /* Find or create an initial vector, U */
  /*-------------------------------------*/

  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "Restart", &RestartFlag));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_FindOrCreatePrimalState(All, RestartFlag, NULL, &U));
  if (ierr != xf_OK) return ierr;

  /* If RestartFlag, copy U to initial condition, UIC */
  if (RestartFlag){
    ierr = xf_Error(xf_CreateVector(&UIC));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_CopyVector(All->Mesh, U, UIC);
    if (ierr != xf_OK) return ierr;
  }
  
  // start timer
  clock_Start = clock();

  /*----------------------------------------------------*/
  /* Create LinR =                                      */
  /*                                                    */
  /* [R(i)_U(0)]^T * [R(i)_U(j)]^-T * J_U(j) * yobs(j)  */
  /*                                                    */
  /* This is -(C*A^-1*F)^T * yobs in Omar Bashir's work */
  /*                                                    */
  /* The minus sign is because LinR goes on the LHS     */
  /*----------------------------------------------------*/

  xf_printf("Creating rhs vector for inverse problem: -(C*A^-1*F)^T * yobs\n");
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Inverse_Residual", xfe_True,
				       xfe_False, NULL, &LinR, NULL));
  if (ierr != xf_OK) return ierr;

  if (SingleOutFlag){

    xf_printf("Computing LinR using adjoint sensitivities.\n");
    
    // LinR = - J_U0(j) * yobs(j)
    for (iOut=0; iOut<nOut; iOut++){
      ierr = xf_Error(xf_VectorMultSet(J_U0[iOut], -OutValue[iOut], 
				       ((iOut==0) ? xfe_Set : xfe_Add), LinR));
      if (ierr != xf_OK) return ierr;
    }

  }
  else{

    /*-------------------------------------------------*/
    /* Run forward solver once to fill in TimeHistData */
    /*-------------------------------------------------*/
    
    xf_printf("Running forward solver once before any inverse attempts.\n");
    
    ierr = xf_Error(xf_ApplyTimeScheme(All, NULL, xfe_False, &U, TimeHistData));
    if (ierr != xf_OK) return ierr;
    
    if (nTime != TimeHistData->nTime){
      xf_printf("# time steps taken in forward solve (%d) != nTime from yobs file (%d).\n",
		TimeHistData->nTime, nTime);
      return xf_Error(xf_INPUT_ERROR);
    }


    // Weight time-history outputs for adjoint solve
    TimeHistData->nSumWeight = nOutput;
    ierr = xf_Error(xf_ReAlloc2((void ***) &TimeHistData->SumOutputWeights,
				nTime, nOutput, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (iTime=0; iTime<nTime; iTime++)
      for (i=0; i<nOutput; i++){
	val = ((TimeFlag == NULL) ? 1.0 : ((real) TimeFlag[iTime]));
	TimeHistData->SumOutputWeights[iTime][i] = val*ObservedOutputs[iTime][i];
	if ((TimeFlag != NULL) && (TimeFlag[iTime] == 1))
	  xf_printf("iTime = %d, iOutput = %d, Observed = %.10E\n",
		    iTime, i, ObservedOutputs[iTime][i]);
      }

    // No longer need ObservedOutputs
    xf_Release( (void *) ObservedOutputs);

  
    // Set TimeWeights to all 1s
    ierr = xf_Error(xf_ReAlloc((void **) &TimeHistData->TimeWeights,
			       nTime, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (iTime=0; iTime<nTime; iTime++) TimeHistData->TimeWeights[iTime] = 1.0;


    // make LinR look like an adjoint vector associated with SumOutput
    if (LinR->OutputName != NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
    LinR->OutputName = SumOutputName; 

    // Initialize LinR to zero
    ierr = xf_Error(xf_SetZeroVector(LinR));
    if (ierr != xf_OK) return ierr;

    // Unsteady Adjoint Solve
    ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, NULL, U, 1, &LinR, TimeHistData));
    if (ierr != xf_OK) return ierr;
  
    LinR->OutputName = NULL; // back to normal

    ierr = xf_Error(xf_VectorMult(LinR, -1.0)); // because LinR goes on LHS
    if (ierr != xf_OK) return ierr;

  }


  /*---------------------------*/
  /* Solve for U in:           */
  /* (H + beta*I)*U + LinR = 0 */
  /*---------------------------*/

  /* If RestartFlag, set U = UIC */
  if (RestartFlag){
    ierr = xf_SetVector(UIC, xfe_Set, U);
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_DestroyVector(UIC, xfe_True)); // no longer need it
    if (ierr != xf_OK) return ierr;
  }


  // Create a VectorSet for storing work Vectors (4 for CG, extra for Param)
  ierr = xf_Error(xf_FindSimilarVectorSet(All, U, 4+nParam, "CGSet", xfe_True, 
					  xfe_False, NULL, &CGSet));
  if (ierr != xf_OK) return ierr;

  // set vectors for parametric inversion
  if (ParamFlag){
    V = CGSet->Vector + 0;
    W = CGSet->Vector + 1;
    R = CGSet->Vector + 2;
    U_mu = CGSet->Vector + 4;
  }

  // Create a temporary vector for applying Hessian
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Inverse_Utemp", xfe_True,
				       xfe_False, NULL, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;

  /*--------------------------------------------*/
  /* If requesting Laplace regularization:      */
  /* Create a diffusion residual term + matrix  */
  /*--------------------------------------------*/

  if (LReg){
    
    // save pointer to original res terms
    ResTermsOrig = All->EqnSet->ResTerms;

    // create new ResTerms for just diffusion
    ierr = xf_Error(xf_Alloc((void **) &ResTermsLReg, 1, sizeof(xf_ResTerms)));
    if (ierr != xf_OK) return ierr;
    
    ResTermsLReg->nResTerm = 0;
    ResTermsLReg->ResTerm  = NULL;

    // add a diffusion term to new residual ResTermsLReg
    ierr = xf_Error(xf_AddResTerm(ResTermsLReg, xfe_ResTermDiff, DiffKeyValue));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_SetZeroVector(Utemp)); // Utemp = 0
    if (ierr != xf_OK) return ierr;

    
    All->EqnSet->ResTerms = ResTermsLReg; // set to new res terms

    // get orig viscosity
    ierr = xf_Error(xf_GetKeyValueReal(All->EqnSet->KeyValue, "Viscosity", &ViscosityOrig));
    if (ierr != xf_OK) return ierr;

    // set new viscosity
    ViscosityLReg = 1.0;
    ierr = xf_Error(xf_SetKeyValueReal(All->EqnSet->KeyValue, "Viscosity", ViscosityLReg));
    if (ierr != xf_OK) return ierr;


    
    // create/allocate SolverData
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    if (ierr != xf_OK) return ierr;


    // locate Jacobian matrix
    ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, "R_UDiff",
					  xfe_True,  NULL, &R_U, NULL));
    if (ierr != xf_OK) return ierr;
  

    // R0 = R(0)  and  R_U = R_U(0)
    ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual_0", xfe_True,
					 xfe_False, NULL, &R0, NULL));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_CalculateResidual(All, Utemp, R0, R_U, SolverData));
    if (ierr != xf_OK) return ierr;


    All->EqnSet->ResTerms = ResTermsOrig; // back to original
    ierr = xf_Error(xf_SetKeyValueReal(All->EqnSet->KeyValue, "Viscosity", ViscosityOrig));
    if (ierr != xf_OK) return ierr;
  }    
  


  // write log to see what's going on
  if (myRank == 0){
    if (xf_NotNull(Prefix)) sprintf(Title, "%s_ICLog.txt", Prefix);
    else sprintf(Title, "ICLog.txt");
    if ((fid = fopen(Title, "a")) != NULL){
      fprintf(fid, "\n*** [Re]starting CG iterations ***\n");
      fclose(fid);
    }
  }


  /* Iterations begin here */
  
  iIter = 0;
  iNewton = 0;
  iTry = 0;
  LinearSolverData = NULL;
  done = xfe_False; 
  // preconditioning allowed only if using Mass/Laplace matrix regularization
  PreconditionFlag = (MReg || LReg) && (Precond); 
  while (!done){

    xf_printf("\n== Starting solution iteration %d ==\n\n", iIter++);

    if (iIter > maxiter_CG){
      xf_printf("Maximum number of iterations exceeded.  Stopping.\n");
      break;
    }

    if (!ParamFlag){
      // Call iteration routine
      ierr = xf_Error(xf_LinearIterCG(All, CGSet, tol, xfe_VerbosityMedium, PreconditionFlag, 
				      !RestartFlag, SavePoint, CGRestart,
				      LinR,U, &V, &W, &LinearSolverData, &Status));
      if (ierr != xf_OK) return ierr;
      
      done = (Status == xfe_LinearConverged);
      if (done) break; // converged


      // Apply preconditioner if necessary
      if (Status == xfe_LinearPrecondition){
	
	xf_printf("\n Applying Preconditioner\n");
	
	// W = V
	ierr = xf_Error(xf_SetVector(V, xfe_Set, W));
	if (ierr != xf_OK) return ierr;
	// W = M{-1}*W
	ierr = xf_Error(xf_MultInvMassMatrix(All, 1.0, NULL, W));
	if (ierr != xf_OK) return ierr;
	continue;
      }

      xf_printf("\n CG Residual norm = %.15E\n", sqrt(LinearSolverData->rnorm2));
      
      // write log to see what's going on
      if (myRank == 0){
	if (xf_NotNull(Prefix)) sprintf(Title, "%s_ICLog.txt", Prefix);
	else sprintf(Title, "ICLog.txt");
	if ((fid = fopen(Title, "a")) != NULL){
	  fprintf(fid, "On CG Iteration %d, residual norm = %.15E\n", 
		  iIter, sqrt(LinearSolverData->rnorm2));
	  fclose(fid);
	}
      }
      
      /* Set W = "Hessian operator applied to V" (regularization included) */
      ierr = xf_Error(xf_ApplyHessian(All, R_U, nOut, J_U0, V, Utemp, TimeHistData, nOutput,
				      SumOutputName, TimeFlag, MReg, betaM, LReg, betaL,
				      beta, SolverData, &BreakFlag, W, J));
      if (ierr != xf_OK) return ierr;

      if (BreakFlag) break;

    }
    else{
      /* Parametric inversion via Newton iteration */

      // form IC -> store in V;  also gradient -> U_mu
      ierr = xf_Error(xf_ICParam2State(All, nParam, ParamName, ParamValue, 
				       epsParam, V, U_mu));
      if (ierr != xf_OK) return ierr;

      // calculate optimization problem residual: R = (H*V + LinR)
      ierr = xf_Error(xf_ApplyHessian(All, R_U, nOut, J_U0, V, Utemp, TimeHistData, nOutput,
				      SumOutputName, TimeFlag, MReg, betaM, LReg, betaL,
				      beta, SolverData, &BreakFlag, R, J));
      if (ierr != xf_OK) return ierr;
      if (BreakFlag) break;
      ierr = xf_Error(xf_SetVector(LinR, xfe_Add, R));
      if (ierr != xf_OK) return ierr;

      // project R to muR
      for (iParam=0; iParam<nParam; iParam++){
	ierr = xf_Error(xf_VectorDot(U_mu+iParam, R, muR+iParam));
	if (ierr != xf_OK) return ierr;
      }
      
      // convergence check, restart IC if necessary
      for (i=0, rnorm=0.; i<nParam; i++) rnorm += muR[i]*muR[i];
      rnorm = sqrt(rnorm);

      // write log to see what's going on
      xf_printf("\n Newton Residual norm = %.15E\n", rnorm);
      xf_printf(" ParamValue = ");
      for (i=0; i<nParam; i++) xf_printf(" %.5E", ParamValue[i]);
      xf_printf("\n");

      if (rnorm < tol){  // converged ... but maybe to a local optimum?

        xf_printf("Converged: rnorm = %.10E\n", rnorm);
	
	// calculate cost function (no regularization)
	for (iOut=0, cost=0.; iOut<nOut; iOut++)
	  cost += 0.5*xf_PowInt(J[iOut]-OutValue[iOut], 2);
	if (cost < costmin){
	  xf_printf("New minimum found, cost = %.10E\n", cost);
	  // test if inside parameter range
	  Inside = xfe_True;
          for (i=0; i<nParam; i++)
	    Inside = Inside && (ParamValue[i] > ParamRange[2*i]) && (ParamValue[i] < ParamRange[2*i+1]);
          if (!Inside) 
	    xf_printf("Oops, not inside the parameter domain.  Not considering this minimum.\n");
	  else{
	    for (i=0; i<nParam; i++) ParamValueMin[i] = ParamValue[i];		
	    costmin = cost;
          }
	}
	if (++iTry < nTry){
	  // generate random new param values
	  for (i=0; i<nParam; i++){
	    ierr = xf_Error(xf_RandUniform(1, &randval)); // note, this is parallel-safe
	    if (ierr != xf_OK) return ierr;
	    ParamValue[i] = ParamRange[2*i]*(1.0-randval) + ParamRange[2*i+1]*randval;
	  } // i
	}
	else{
	  // form IC from ParamValueMin -> store in U
	  ierr = xf_Error(xf_ICParam2State(All, nParam, ParamName, ParamValueMin, 
					   epsParam, U, NULL));
	  if (ierr != xf_OK) return ierr;
	  done = xfe_True; 
	  break;
	}
      }
	
      if (myRank == 0){
	if (xf_NotNull(Prefix)) sprintf(Title, "%s_ICLog.txt", Prefix);
	else sprintf(Title, "ICLog.txt");
	if ((fid = fopen(Title, "a")) != NULL){
	  fprintf(fid, "On Newton Iteration %d (iNewton =%d), residual norm = %.15E\n", 
		  iIter, iNewton, rnorm);
	  fclose(fid);
	}
      }

      // restart IC if taking too long
      iNewton++;
      if (iNewton > nNewtonMax){
	xf_printf("Newton not converging.  Restarting IC with random guess\n");
	iNewton = 0;
	for (i=0; i<nParam; i++){
	  ierr = xf_Error(xf_RandUniform(1, &randval)); // this is parallel-safe
	  if (ierr != xf_OK) return ierr;
	  ParamValue[i] = ParamRange[2*i]*(1.0-randval) + ParamRange[2*i+1]*randval;
	} // i
      }

      // apply Hessian to U_mu and form projected Jacobian: muJ
      for (j=0; j<nParam; j++){
	// W = H*U_mu[iParam]
	ierr = xf_Error(xf_ApplyHessian(All, R_U, nOut, J_U0, U_mu+j, Utemp, TimeHistData, nOutput,
					SumOutputName, TimeFlag, MReg, betaM, LReg, betaL,
					beta, SolverData, &BreakFlag, W, J));
	if (ierr != xf_OK) return ierr;
	// dot products to set columns of muJ
	for (i=0; i<nParam; i++){
	  ierr = xf_Error(xf_VectorDot(U_mu+i, W, muJ+i*nParam+j));
	  if (ierr != xf_OK) return ierr;
	}
	  
      } // j

      // Solve muJ*dmu + muR = 0
      ierr = xf_Error(xf_ComputePLU(muJ, nParam, P));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<nParam; i++) muR[i] *= -1.0;
      ierr = xf_Error(xf_SolvePLU(muJ, P, muR, nParam, dmu, NULL));
      if (ierr != xf_OK) return ierr;

      // sanity check on param range
      for (i=0; i<nParam; i++){
	dParam = paramOmega*(ParamRange[2*i+1] - ParamRange[2*i]); // positive
	dmu[i] = max(dmu[i], -dParam);
	dmu[i] = min(dmu[i], dParam);
	xf_printf("dmu[%d] = %.6E\n", i, dmu[i]);
      }

      // take an iteration of Newton
      for (i=0; i<nParam; i++) ParamValue[i] += dmu[i];

    }
    
    // if Halted, break
    if (xf_CheckUserHalt(NULL)) break;

  } // end CG linear solver loop


  // Solution is stored in U
  xf_printf("\nFinished solving inverse problem.\n");

  if (ParamFlag){
    xf_printf("\nParametric solution = ");
    for (i=0; i<nParam; i++) xf_printf(" %.8E", ParamValueMin[i]);
    xf_printf("\n");
    xf_printf("Cost function = %.8E\n", costmin);
  }


  /*---------------------*/
  /* Release some memory */
  /*---------------------*/

  // Destroy Time history 
  ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
  if (ierr != xf_OK) return ierr;

  // destroy LinR
  ierr = xf_Error(xf_DestroyVector(LinR, xfe_True));
  if (ierr != xf_OK) return ierr;

  // destroy Utemp
  ierr = xf_Error(xf_DestroyVector(Utemp, xfe_True));
  if (ierr != xf_OK) return ierr;

  // No longer need CGSet
  ierr = xf_Error(xf_DestroyVectorSet(CGSet));
  if (ierr != xf_OK) return ierr;
  
  // Destroy the temporarily-created SumOutput
  if (nOutput > 0){
    ierr = xf_Error(xf_ReAllocOutputs(All->EqnSet->Outputs, nOutput));
    if (ierr != xf_OK) return ierr;
  }

  // Destroy ResTermsLReg and R0 if was doing Laplace regularization
  if (LReg){
    ierr = xf_Error(xf_DestroyResTerms(ResTermsLReg));
    if (ierr != xf_OK) return ierr;

    // destroy R0
    ierr = xf_Error(xf_DestroyVector(R0, xfe_True));
    if (ierr != xf_OK) return ierr;

    // destroy SolverData
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    if (ierr != xf_OK) return ierr;

  }

  // end timer
  clock_End = clock(); // end timer
  xf_printf("CPU time = %.10E\n", ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC);



  /*--------------------------------------*/
  /* Write out Initial Condition solution */
  /*--------------------------------------*/

  xf_printf("\nWriting initial condition U to %s.\n", solFile);

  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DataSetAdd(DataSet, "State", xfe_Vector,
				xfe_True, (void *) U, &D));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, solFile));
  if (ierr != xf_OK) return ierr;
  D->Data = NULL;
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;


  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;

  /* Destroy sensitivity vectors */
  for (iOut=0; iOut<nOut; iOut++){
    ierr = xf_Error(xf_DestroyVector(J_U0[iOut], xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release( (void *) J_U0);
  
  // Release memory
  xf_Release( (void *) TimeMask);
  xf_Release( (void *) TimeFlag);
  xf_Release2( (void **) ParamName);
  xf_Release(  (void  *) ParamRange);
  xf_Release(  (void  *) ParamValue);
  xf_Release(  (void  *) ParamValueMin);
  xf_Release(  (void  *) P);
  xf_Release(  (void  *) muR);
  xf_Release(  (void  *) dmu);
  xf_Release(  (void  *) muJ);
  xf_Release(  (void  *) epsParam);
  xf_Release2( (void **) OutName);
  xf_Release(  (void  *) OutTimeIndex);
  xf_Release(  (void  *) OutValue);
  xf_Release(  (void  *) J);
 
  xf_printf("xf_InverseIC finished.\n");

  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
