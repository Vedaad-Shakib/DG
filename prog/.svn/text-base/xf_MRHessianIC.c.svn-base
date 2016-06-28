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
  FILE:  xf_MRHessianIC.c

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
#include "xf_EigSolver.h"
#include "xf_MathLapack.h"
#include "xf_Arg.h"

// minimum number of Lanczos vectors (in case only want 1 or 2 eigs)
#define MIN_LANCZOS 1
#define MAXITER_LANCZOS 5000

/******************************************************************/
//   FUNCTION Definition: xf_DumpVectorSet
static int 
xf_DumpVectorSet(xf_VectorSet *VS, int n, const char *fname)
{ 
  int ierr, egrp, elem, j, r;
  real **EV;
  xf_Vector *V;
  FILE *fid;

  
  // open file
  if ((fid = fopen(fname, "w")) == NULL) return xf_Error(xf_FILE_WRITE_ERROR);
  
  ierr = xf_Error(xf_Alloc( (void **) &EV, n, sizeof(real *)));
  if (ierr != xf_OK) return ierr;

  V = VS->Vector + 0;
  for (egrp=0; egrp<V->nArraySelf; egrp++){
    for (elem=0; elem<V->GenArray[egrp].n; elem++){
      for (j=0; j<n; j++) EV[j] = VS->Vector[j].GenArray[egrp].rValue[elem];
      for (r=0; r<V->GenArray[egrp].r; r++){
	for (j=0; j<n; j++)
	  fprintf(fid, "%.15E ", EV[j][r]);
	fprintf(fid, "\n");
      }
    }
  }
  fclose(fid);

  xf_Release( (void *) EV);
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_SetBasisSetToIdentity
static int 
xf_SetBasisSetToIdentity(xf_All *All, int nBasis, xf_VectorSet *BasisSet)
{
  /* Does what the name suggests */
  int ierr, j, k, r, sr, nelemtot, ndof;
  int egrp, elem, nn;
  xf_Vector *U;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  xf_printf(" Setting BasisSet to Identity.\n");

  
  U = BasisSet->Vector + 0;
  r = U->GenArray[0].r;
  sr = U->StateRank;
  ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;

  ndof = nelemtot*r;

  if (nBasis != ndof){
    xf_printf("Error, nBasis = %d is not ndof = %d in setting BasisSet = I\n",
	      nBasis, ndof);
    return xf_Error(xf_INPUT_ERROR);
  }
  
  egrp = 0; elem = 0; j = 0;
  for (k=0; k<nBasis; k++){
    U = BasisSet->Vector + k;
    if (egrp >= U->nArraySelf) return xf_Error(xf_OUT_OF_BOUNDS);
    ierr = xf_Error(xf_SetZeroVector(U));
    if (ierr != xf_OK) return ierr;  
    U->GenArray[egrp].rValue[elem][j] = 1.0;

    ierr = xf_Error(xf_Order2nNode(U->Basis[egrp], U->Order[egrp], &nn));
    if (ierr != xf_OK) return ierr;
    if ((++j) == nn*sr){
      j = 0;
      if ((++elem) == U->GenArray[egrp].n){
	elem = 0;
	egrp++;
      }
    }
  }
  
  xf_printf(" done.\n");

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int len, ierr, i, j, k, ndof, r;
  int outerr;
  int myRank, nProc;
  int nOutput;
  int iIter, nLanczos, nIC, nIC0, ICstart, ICend, nTime, iTime;
  int WriteHess = 0, SnapStart = 1, WriteInterval, SaveInterval;
  int iIC, nSnap, iSnap, nBasis, nTimeStep, nTimeMask;
  int nEigRestart = 0, nev_restart = 0;
  enum xfe_Bool ParallelFlag, RestartFlag, done, JustPOD, CheckEig;
  enum xfe_Bool EigRestart, LanczosRestart, SavePoint, PODFlag;
  enum xfe_Bool RetakeSnap, SkipSnap, SmallerPOD;
  enum xfe_TimeSchemeType TimeScheme;
  enum xfe_EigStatusType Status;
  char *ArgIn[] = {"job", "NULL", ".job file name to read (run parameters)",
		   "data", "IC.data", "name of IC data file",
		   "nIC", "10", "number of IC (eig) vectors to calculate",
		   "nIC0", "-1", "if > 0, # original ICs; useful if want subset of all ICs",
		   "tol", "1e-6", "relative tolerance on eigenvalues",
		   "WriteHess", "0", "set to n to dump n cols of H to Hessian.txt",
		   "nBasis", "1", "number of basis fcns for ROM",
		   "POD", "True", "if True, POD will be performed",
		   "JustPOD", "False", "snapshots already on disk; just do POD",
		   "CheckEig", "False", "if True, Hessian eigvectors will be checked",
		   "SnapStart", "1", "# WriteIntervals at which snapshot storage begins",
		   "Verbosity", "Medium", "Low/Medium/High -- amount of printouts to screen",
		   "Restart", "False", "if True, immediate restart from Lanczos.data/txt",
		   "EigRestart", "False", "if True, restart from IC.data and Eig.txt",
		   "SaveInterval", "1", "how often to write Lanczos.data/txt; < 0 to not write",
		   "RetakeSnap", "False", "if True, snapshots are retaken using IC data file",
		   "SkipSnap", "False", "if True, snapshot taking is skipped",
		   "SmallerPOD", "False", "pull off smaller ROM.m from existing Basis.data",
		   "TimeMask", "NULL", "time mask file for fewer outputs",
		   "\0"};
  char jobFile[xf_MAXSTRLEN] = "NULL"; 
  char ICFile[xf_MAXSTRLEN];
  char SavePrefix0[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN];
  char DataFile[xf_MAXSTRLEN];
  char TimeMaskFile[xf_MAXSTRLEN] = "NULL";
  char Verbosity[xf_MAXSTRLEN] = "Medium";
  char SumOutputName[] = "SumOutputs";
  char s[xf_MAXSTRLEN], line[xf_MAXLINELEN];
  char *DiffKeyValue[] = {"Unused", "Unused", "\0"};
  int *TimeMask = NULL, *TimeFlag = NULL;
  FILE *fid;
  real EigTol, *E;
  real *Rr, *Mr, *Ar, *Br, *Cr, val, val0;
  real ViscosityLReg, ViscosityOrig;
  xf_KeyValue KeyValueArg;
  xf_TimeHistData *TimeHistData = NULL;
  xf_EigSolverData *EigSolverData;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Vector *U, *V, *W, *R, *R0, *UT;
  xf_VectorSet  *LanczosSet, *ICSet = NULL;
  xf_VectorSet  *SnapSet, *BasisSet;
  xf_Output *SumOutput;
  xf_ResTerms *ResTermsLReg = NULL;
  xf_ResTerms *ResTermsOrig = NULL;
  xf_All *All;

  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if (nProc > 1) ParallelFlag = xfe_True;
  else  ParallelFlag = xfe_False;

    
  xf_printf("\n");
  xf_printf("=== Linear Model Reduction for Initial Conditions using Hessian ===\n");
  xf_printf("\n");
    
        
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValueArg));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValueArg);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  /* Get jobFile */
  ierr = xf_Error(xf_GetKeyValue(KeyValueArg, "job", jobFile));
  if (ierr != xf_OK) return ierr;

  /* Get ICFile */
  ierr = xf_Error(xf_GetKeyValue(KeyValueArg, "data", ICFile));
  if (ierr != xf_OK) return ierr;

  /* Get number of ICs eigenvectors */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "nIC", &nIC));
  if (ierr != xf_OK) return ierr;
  if (nIC <= 0) return xf_Error(xf_INPUT_ERROR); // nothing to work with

  /* Get number of original IC eigenvectors (useful if just doing POD with existing snapshots) */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "nIC0", &nIC0));
  if (ierr != xf_OK) return ierr;

  /* Get eigenvalue tolerance */
  ierr = xf_Error(xf_GetKeyValueReal(KeyValueArg, "tol", &EigTol));
  if (ierr != xf_OK) return ierr;
  if (EigTol <= 0) return xf_Error(xf_INPUT_ERROR);

  /* Get number of Basis vectors for POD eigenvectors */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "nBasis", &nBasis));
  if (ierr != xf_OK) return ierr;
  if (nBasis <= 0) return xf_Error(xf_INPUT_ERROR); // need nBasis > 0

  /* PODFlag */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValueArg, "POD", &PODFlag));
  if (ierr != xf_OK) return ierr;
    
  /* JustPOD */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValueArg, "JustPOD", &JustPOD));
  if (ierr != xf_OK) return ierr;

  /* CheckEig */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValueArg, "CheckEig", &CheckEig));
  if (ierr != xf_OK) return ierr;

  /* SnapStart */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "SnapStart", &SnapStart));
  if (ierr != xf_OK) return ierr;
  if (SnapStart < 0) return xf_Error(xf_INPUT_ERROR);

  /* WriteHess */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "WriteHess", &WriteHess));
  if (ierr != xf_OK) return ierr;

  /* Verbosity */
  ierr = xf_Error(xf_GetKeyValue(KeyValueArg, "Verbosity", Verbosity));
  if (ierr != xf_OK) return ierr;

  /* LanczosRestart */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValueArg, "Restart", &LanczosRestart));
  if (ierr != xf_OK) return ierr;

  /* EigRestart */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValueArg, "EigRestart", &EigRestart));
  if (ierr != xf_OK) return ierr;

  /* SaveInterval */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "SaveInterval", &SaveInterval));
  if (ierr != xf_OK) return ierr;

  /* RetakeSnap */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValueArg, "RetakeSnap", &RetakeSnap));
  if (ierr != xf_OK) return ierr;

  /* SkipSnap */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValueArg, "SkipSnap", &SkipSnap));
  if (ierr != xf_OK) return ierr;

  /* SmallerPOD */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValueArg, "SmallerPOD", &SmallerPOD));
  if (ierr != xf_OK) return ierr;

  /* Get TimeMaskFile */
  ierr = xf_GetKeyValue(KeyValueArg, "TimeMask", TimeMaskFile);
  if (ierr != xf_OK) return ierr;

  // destroy key-value from arg list
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValueArg));
  if (ierr!=xf_OK) return ierr;


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
  
  /* Find or create an initial vector, U */
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "Restart", &RestartFlag));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_FindOrCreatePrimalState(All, RestartFlag, NULL, &U));
  if (ierr != xf_OK) return ierr;

  /* Set verbosity */
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", Verbosity));
  if (ierr != xf_OK) return ierr;

  // do not write a .log file
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "WriteLog", "False"));
  if (ierr != xf_OK) return ierr;


  // pull off and save SavePrefix
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix0));
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


  // save away original UnsteadyWriteInterval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 
				    &WriteInterval));
  if (ierr != xf_OK) return ierr;

  if (All->EqnSet->Outputs == NULL) return xf_Error(xf_INPUT_ERROR);

  nOutput = All->EqnSet->Outputs->nOutput;

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
  
  /*-------------------*/
  /* Read TimeMaskFile */
  /*-------------------*/

  if (xf_NotNull(TimeMaskFile)){
    /* Read time mask file on root */
    outerr = xf_OK;
    if (myRank == 0){
      
      /* Read in yobs */
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

  }



 
  /*----------------------------------------------*/
  /* Write Hessian and exit if user requests this */
  /*----------------------------------------------*/

  if (WriteHess){

    if (ParallelFlag) return xf_Error(xf_NOT_SUPPORTED); // since we use identity

    // count ndof = # degrees of freedom
    for (i=0, ndof=0; i<U->nArraySelf; i++)
      ndof += U->GenArray[i].n*U->GenArray[i].r;

    if (U->nArraySelf > 1) return xf_Error(xf_NOT_SUPPORTED); // for now

    xf_printf("ndof = %d, WriteHess = %d\n", ndof, WriteHess);

    ierr = xf_Error(xf_FindSimilarVectorSet(All, U, WriteHess, "LanczosSet", 
					    xfe_True, xfe_True, &D, &LanczosSet));
    if (ierr != xf_OK) return ierr; 

    for (k=0; k<WriteHess; k++){

      xf_printf("\n\n\n k = %d \n\n\n", k);

      // Set U = e(k) -- assumes only one element group
      ierr = xf_Error(xf_SetZeroVector(U));
      if (ierr != xf_OK) return ierr;
      r = U->GenArray[0].r;
      U->GenArray[0].rValue[k/r][k%r] = 1.0;

      // Unsteady Forward Solve
      ierr = xf_Error(xf_ApplyTimeScheme(All, NULL, xfe_False, &U, TimeHistData));
      if (ierr != xf_OK) return ierr;
      
      nTime = TimeHistData->nTime;
      
      // Create TimeFlag
      if (TimeMask != NULL){
	ierr = xf_Error(xf_ReAlloc((void **) &TimeFlag, nTime, sizeof(int)));
	if (ierr != xf_OK) return ierr;
	// Fill TimeFlag
	for (i=0; i<nTime; i++) TimeFlag[i] = 0;
	for (i=0; i<nTimeMask; i++){
	  if ((TimeMask[i] < 0) || (TimeMask[i] >= nTime)) return xf_Error(xf_OUT_OF_BOUNDS);
	  TimeFlag[TimeMask[i]] = 1;
	}
      }

      // Weight time-history outputs for adjoint solve
      TimeHistData->nSumWeight = nOutput;
      ierr = xf_Error(xf_ReAlloc2((void ***) &TimeHistData->SumOutputWeights,
				  nTime, nOutput, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      for (iTime=1; iTime<nTime; iTime++)
	for (i=0; i<nOutput; i++){
	  val = ((TimeFlag == NULL) ? 1.0 : ((real) TimeFlag[iTime]));
	  TimeHistData->SumOutputWeights[iTime][i] = val*TimeHistData->OutputValues[i][iTime];
	}
      
      // Set TimeWeights to all 1s
      ierr = xf_Error(xf_ReAlloc((void **) &TimeHistData->TimeWeights,
				 nTime, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      for (iTime=0; iTime<nTime; iTime++) TimeHistData->TimeWeights[iTime] = 1.0;
      

      W = LanczosSet->Vector+k;
      
      // make W look like an adjoint vector associated with SumOutput
      if (W->OutputName != NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
      W->OutputName = SumOutputName; 
      
      // Initialize W to zero
      ierr = xf_Error(xf_SetZeroVector(W));
      if (ierr != xf_OK) return ierr;
      
      // Unsteady Adjoint Solve
      ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, NULL, U, 1, &W, TimeHistData));
      if (ierr != xf_OK) return ierr;
      
    } // k

    if (myRank == 0){
      ierr = xf_Error(xf_DumpVectorSet(LanczosSet, WriteHess, "Hessian.txt"));
      if (ierr != xf_OK) return ierr;
      xf_printf("Wrote Hessian.txt\n");
    }
    
    return xf_OK;
  }


  if (SmallerPOD) JustPOD = xfe_True;

  if (!JustPOD){

    // Create a VectorSet for storing Lanczos Vectors
    nLanczos = max(2*nIC+2, MIN_LANCZOS);
    ierr = xf_Error(xf_FindSimilarVectorSet(All, U, nLanczos+1, "LanczosSet", 
					    xfe_True, xfe_False, NULL, &LanczosSet));
    if (ierr != xf_OK) return ierr; 


    // Allocate vector for eigenvalues
    ierr = xf_Error(xf_Alloc( (void **) &E, nLanczos+1, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    
    // load existing IC.data and Eig.txt if asking for a restart
    if (EigRestart){

      xf_printf("Restarting from existing IC.data and Eig.txt\n");

      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, ICFile, DataSet));
      if (ierr != xf_OK) return ierr;

      D = DataSet->Head;
      if (strncmp(D->Title, "ICs", 3) != 0) return xf_Error(xf_INPUT_ERROR);
      ICSet = (xf_VectorSet *) D->Data;

      nEigRestart = ICSet->nVector;

      xf_printf("nEigRestart = %d\n", nEigRestart);

      if (nEigRestart >= nIC) return xf_Error(xf_INPUT_ERROR);

      // modify nLanczos -- will not need 2*nIC+2
      //nLanczos = min(nLanczos, nEigRestart + 2*(nIC-nEigRestart) + 2);

      // Copy ICs to LanczosSet
      for (i=0; i<nEigRestart; i++){
	ierr = xf_Error(xf_SetVector(ICSet->Vector+i, xfe_Set, LanczosSet->Vector+i));
	if (ierr != xf_OK) return ierr;
      }
	
      // destroy dataset
      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;

      xf_pprintf("Here.\n"); fflush(stdout);


      // Read in Eig.txt for eigenvalues
      outerr = xf_OK;
      if (myRank == 0){
	if ((fid = fopen("Eig.txt", "r")) == NULL) outerr = xf_Error(xf_FILE_READ_ERROR);
	else{
	  for (i=0; (i<nEigRestart) && (outerr==xf_OK); i++){
	    if (fgets(line, xf_MAXLINELEN, fid) == NULL) outerr = xf_Error(xf_FILE_READ_ERROR);
	    if (sscanf(line, "%lf", E+i) != 1) outerr = xf_Error(xf_FILE_READ_ERROR);
	  }
	}
	fclose(fid);
      }
      // broadcast outerr and quit if not OK
      if (xf_PError(&outerr, 0) != xf_OK) return outerr;

      // broadcast eigs to all procs
      ierr = xf_Error(xf_MPI_Bcast((void *) E, nEigRestart*sizeof(real), 0));
      if (ierr != xf_OK) return ierr;

    } // end if restart requested
    

    // Create a VectorSet for storing eigenvectors
    ierr = xf_Error(xf_FindSimilarVectorSet(All, U, nIC, "ICs", xfe_True, 
					    xfe_False, NULL, &ICSet));
    if (ierr != xf_OK) return ierr; 


    /*------------------------*/
    /* Start eigensolver loop */
    /*------------------------*/

    iIter = 0;
    EigSolverData = NULL;
    done = xfe_False;
    while (!done){

      // Write a status file with iteration number (for debugging)
      if (myRank == 0){
	if ((fid = fopen("Status.txt", "w")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
	fprintf(fid, "Hessian eigenvector calculation: on iteration # %d\n", iIter);
	fclose(fid);
      }


      // only ask for an IC restart on first iteration
      nev_restart = ((iIter == 0) ? nEigRestart : 0);

      xf_printf("\n== Starting Lanczos Iteration %d ==\n\n", iIter++);

      if (iIter > MAXITER_LANCZOS){
	xf_printf("Maximum number of calls to EigIterLanczos exceeded.  Stopping.\n");
	break;
      }

      // Should we save a checkpoint?
      SavePoint = ((SaveInterval > 0) && ((iIter % SaveInterval) == 0));

      // Call iteration routine
      ierr = xf_Error(xf_EigIterLanczos(All, nIC, nLanczos, LanczosSet, EigTol, 
					xfe_VerbosityHigh, nev_restart, SavePoint, 
					LanczosRestart, E, ICSet, &V, &W, 
					&EigSolverData, &Status));
      if (ierr != xf_OK) return ierr;

      done = ((Status == xfe_EigConverged) || (Status == xfe_EigForceExit));
      if (done) break; // converged

      /*---------------------------------*/
      /* Set W = "operator applied to V" */
      /*---------------------------------*/

      xf_printf(" Applying Hessian\n");

      // First copy V into U so that we do not overwrite V
      ierr = xf_Error(xf_SetVector(V, xfe_Set, U));
      if (ierr != xf_OK) return ierr;

      // Unsteady Forward Solve
      ierr = xf_Error(xf_ApplyTimeScheme(All, NULL, xfe_False, &U, TimeHistData));
      if (ierr != xf_OK) return ierr;

      nTime = TimeHistData->nTime;

      // Create TimeFlag
      if (TimeMask != NULL){
	ierr = xf_Error(xf_ReAlloc((void **) &TimeFlag, nTime, sizeof(int)));
	if (ierr != xf_OK) return ierr;
	// Fill TimeFlag
	for (i=0; i<nTime; i++) TimeFlag[i] = 0;
	for (i=0; i<nTimeMask; i++){
	  if ((TimeMask[i] < 0) || (TimeMask[i] >= nTime)) return xf_Error(xf_OUT_OF_BOUNDS);
	  TimeFlag[TimeMask[i]] = 1;
	}
      }
    
      // Weight time-history outputs for adjoint solve
      TimeHistData->nSumWeight = nOutput;
      ierr = xf_Error(xf_ReAlloc2((void ***) &TimeHistData->SumOutputWeights,
				  nTime, nOutput, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      for (iTime=1; iTime<nTime; iTime++)
	for (i=0; i<nOutput; i++){
	  val = ((TimeFlag == NULL) ? 1.0 : ((real) TimeFlag[iTime]));
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
      ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, NULL, U, 1, &W, TimeHistData));
      if (ierr != xf_OK) return ierr;
    
      W->OutputName = NULL; // back to normal

      // Check if halted by user
      if (xf_CheckUserHalt(NULL)){
	xf_printf("WARNING: computation halted by user. Do not trust the Eigenvalues/vectors.\n");
	return xf_Error(xf_NOT_SUPPORTED); // for now
	break;
      }
    
    } // end Lanczos eiengsolver loop


    /*---------------------*/
    /* Release some memory */
    /*---------------------*/

    // No longer need LanczosSet
    ierr = xf_Error(xf_DestroyVectorSet(LanczosSet));
    if (ierr != xf_OK) return ierr;
  


    /*-----------------------*/
    /* Print out eigenvalues */
    /*-----------------------*/

    xf_printf("\nFinished calculating eigenvalues.\n");
    xf_printf("Eigenvalues:\n");
    for (i=0; i<nIC; i++) xf_printf(" Eig[%d] = %.15E\n", i, E[i]);

    if (myRank == 0){ // also print to file Eig.txt
      if ((fid = fopen("Eig.txt", "w")) == NULL) return xf_Error(xf_FILE_WRITE_ERROR);
      for (i=0; i<nIC; i++)
	fprintf(fid, "%.15E\n", E[i]);
      if (fclose(fid)!= 0) return xf_Error(xf_FILE_READ_ERROR);
    }



    /*---------------------------------*/
    /* Write out eigenvector vectorset */
    /*---------------------------------*/

    xf_printf("\nWriting initial condition VectorSet to %s.\n", ICFile);

    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DataSetAdd(DataSet, "ICs", xfe_VectorSet,
				  xfe_True, (void *) ICSet, &D));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, ICFile));
    if (ierr != xf_OK) return ierr;
    D->Data = NULL;
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;




    // Check eigenvectors if desired

    if (CheckEig){
      for (iIC=0; iIC<nIC; iIC++){
	
	xf_printf("\nTesting whether V = IC # %d is an eigenvector of H\n", iIC);

	// copy initial condition into U
	ierr = xf_Error(xf_SetVector(ICSet->Vector+iIC, xfe_Set, U));
	if (ierr != xf_OK) return ierr;

	// Compute initial norm val0 = ||U||
	ierr = xf_Error(xf_VectorNorm(U, 2, &val0));
	if (ierr != xf_OK) return ierr;
	
	xf_printf(" Initial ||V|| = %.10E\n", val0);

	xf_printf(" Applying Hessian\n");

	// Unsteady Forward Solve
	ierr = xf_Error(xf_ApplyTimeScheme(All, NULL, xfe_False, &U, TimeHistData));
	if (ierr != xf_OK) return ierr;

	nTime = TimeHistData->nTime;

	// Create TimeFlag
	if (TimeMask != NULL){
	  ierr = xf_Error(xf_ReAlloc((void **) &TimeFlag, nTime, sizeof(int)));
	  if (ierr != xf_OK) return ierr;
	  // Fill TimeFlag
	  for (i=0; i<nTime; i++) TimeFlag[i] = 0;
	  for (i=0; i<nTimeMask; i++){
	    if ((TimeMask[i] < 0) || (TimeMask[i] >= nTime)) return xf_Error(xf_OUT_OF_BOUNDS);
	    TimeFlag[TimeMask[i]] = 1;
	  }
	}

	// Weight time-history outputs for adjoint solve
	TimeHistData->nSumWeight = nOutput;
	ierr = xf_Error(xf_ReAlloc2((void ***) &TimeHistData->SumOutputWeights,
				    nTime, nOutput, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	for (iTime=1; iTime<nTime; iTime++)
	  for (i=0; i<nOutput; i++){
	    val = ((TimeFlag == NULL) ? 1.0 : ((real) TimeFlag[iTime]));
	    TimeHistData->SumOutputWeights[iTime][i] = val*TimeHistData->OutputValues[i][iTime];
	  }
    
	// Set TimeWeights to all 1s
	ierr = xf_Error(xf_ReAlloc((void **) &TimeHistData->TimeWeights,
				   nTime, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	for (iTime=0; iTime<nTime; iTime++) TimeHistData->TimeWeights[iTime] = 1.0;


	// make U look like an adjoint vector associated with SumOutput
	if (U->OutputName != NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
	U->OutputName = SumOutputName; 

	// Initialize U to zero
	ierr = xf_Error(xf_SetZeroVector(U));
	if (ierr != xf_OK) return ierr;

	// Unsteady Adjoint Solve
	ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, NULL, ICSet->Vector+iIC, 1, 
						  &U, TimeHistData));
	if (ierr != xf_OK) return ierr;
    
	U->OutputName = NULL; // back to normal

	// Test whether U - E[iIC]*ICSet->Vector+iIC ~ 0
	ierr = xf_Error(xf_VectorMultSet(ICSet->Vector+iIC, E[iIC], xfe_Sub, U));
	if (ierr != xf_OK) return ierr;

	// val = ||U||
	ierr = xf_Error(xf_VectorNorm(U, 2, &val));
	if (ierr != xf_OK) return ierr;

	xf_printf(" Discrete L2 error ||H*V - lambda*V|| = %.10E\n", val);

      } // iIC
    }


    /*---------------------*/
    /* Release some memory */
    /*---------------------*/

    // Destroy Time history 
    ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
    if (ierr != xf_OK) return ierr;
  
    // Destroy the temporarily-created SumOutput
    ierr = xf_Error(xf_ReAllocOutputs(All->EqnSet->Outputs, nOutput));
    if (ierr != xf_OK) return ierr;

    xf_Release( (void *) E);

    
   


  } // JustPOD
  //------------------


  if (((!JustPOD) && (!SkipSnap)) || (RetakeSnap)){
    
    // Check if have ICSet -- read in if do not have it
    if (ICSet == NULL){
      // create dataset for ICSet reading
      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;
      
      xf_printf("Reading in IC.data\n"); fflush(stdout);
      
      // read in ICs.data
      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, ICFile, DataSet));
      if (ierr != xf_OK) return ierr;
      
      xf_printf("Done\n"); fflush(stdout);
      
      D = DataSet->Head;
      if (strncmp(D->Title, "ICs", 3) != 0) return xf_Error(xf_INPUT_ERROR);
      ICSet = (xf_VectorSet *) D->Data;

      // destroy DataSet
      D->Data = NULL; // so we do not destroy ICSet!
      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;
    }


    /*--------------------------------------------*/
    /* Take snapshots with each initial condition */
    /*--------------------------------------------*/

    // Create a time history data structure    
    ierr = xf_Error(xf_CreateUniformTimeHistData(All, NULL, &TimeHistData));
    if (ierr != xf_OK) return ierr;
    
    for (iIC=0; iIC<nIC; iIC++){
      
      xf_printf("\nTaking time snapshots with IC # %d\n\n", iIC);
      
      // prefix for saving away unsteady snapshot .data files
      sprintf(SavePrefix, "%s_IC%d", SavePrefix0, iIC);
      
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
      if (ierr != xf_OK) return ierr;
      
      // copy initial condition into U
      ierr = xf_Error(xf_SetVector(ICSet->Vector+iIC, xfe_Set, U));
      if (ierr != xf_OK) return ierr;
      
      
      // apply forward-time solve
      ierr = xf_Error(xf_ApplyTimeScheme(All, SavePrefix, xfe_False, &U, TimeHistData));
      if (ierr != xf_OK) return ierr;
      
      // Set Time to first value in history (for next solve)
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", TimeHistData->Time[0]));
      if (ierr != xf_OK) return ierr;
      
    } // iIC
    
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix0));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Time history */
    ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
    if (ierr != xf_OK) return ierr;
    
  
    /*---------------*/
    /* Destroy ICSet */
    /*---------------*/

    ierr = xf_Error(xf_DestroyVectorSet(ICSet));
    if (ierr != xf_OK) return ierr;


  }// end of taking snapshots



  if (PODFlag){
  

    if (!SmallerPOD){

      /*----------------------------------------------------*/
      /* Create a large vectorset for storing all snapshots */
      /*----------------------------------------------------*/

      ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nTimeStep", &nTimeStep));
      if (ierr != xf_OK) return ierr;

      xf_printf("SnapStart = %d, WriteInterval = %d\n", SnapStart, WriteInterval);

      if (WriteInterval <= 0) return xf_Error(xf_INPUT_ERROR);
      nSnap = nIC*(nTimeStep/WriteInterval + 1 - SnapStart);

      if (nSnap <= 0){
	xf_printf("Error, nSnap = %d.  Maybe decrease SnapStart?\n", nSnap);
	return xf_Error(xf_OUT_OF_BOUNDS);
      }

      ierr = xf_Error(xf_CreateVectorSet(nSnap, &SnapSet));
      if (ierr != xf_OK) return ierr;


      /*----------------*/
      /* Load snapshots */
      /*----------------*/
  
      /* Determine which initial conditions to load (ICs ordered in INCREASING eigenvalues) */
      ICstart = 0;
      ICend   = nIC;

      if (nIC0 > 0){
	ICend   = nIC0;
	ICstart = ICend-nIC;
	if (ICstart < 0){
	  xf_printf("Error, ICstart = %d < 0\n", ICstart);
	  return xf_Error(xf_OUT_OF_BOUNDS);
	}
	xf_printf("ICstart = %d, ICend = %d\n", ICstart, ICend);
      } 

      /* Loop over snapshots and read them into SnapSet */
      xf_printf("Loading snapshots ... ");

      iSnap = 0;
      for (iIC=ICstart; iIC<ICend; iIC++){ // loop over initial conditions
	for (iTime=SnapStart*WriteInterval; iTime<=nTimeStep; iTime+=WriteInterval, iSnap++){ // and time steps

	  sprintf(DataFile, "%s_IC%d_U%d.data", SavePrefix0, iIC, iTime);

	  ierr = xf_Error(xf_CreateDataSet(&DataSet));
	  if (ierr != xf_OK) return ierr;

	  ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, DataFile, DataSet));
	  if (ierr!=xf_OK) return ierr;

	  ierr = xf_Error(xf_FindPrimalState(DataSet, 0, &D, NULL));
	  if (ierr != xf_OK) return ierr;
	  UT = (xf_Vector *) D->Data;
	  SnapSet->Vector[iSnap] = (*UT);
	  xf_InitVector(UT);
    
	  ierr = xf_Error(xf_DestroyDataSet(DataSet));
	  if (ierr != xf_OK) return ierr;

	} // iTime
    
      } // iIC
      xf_printf("done.\n");



      /*-----------------------------*/
      /* Calculate POD Basis vectors */
      /*-----------------------------*/

      xf_printf("Calculating basis vectors from snapshots (POD) ...\n");
      xf_printf("  nSnap = %d\n", nSnap);
      xf_printf("  nBasis = %d\n", nBasis);
      if (nBasis <= 0){
	xf_printf(" Need nBasis > 0.\n");
	return xf_Error(xf_INPUT_ERROR);
      }

      ierr = xf_Error(xf_FindSimilarVectorSet(All, SnapSet->Vector+0, nBasis, "BasisSet",
					      xfe_True, xfe_False, NULL, &BasisSet));
      if (ierr != xf_OK) return ierr;

      /* Perform POD -> basis vectors. */
      ierr = xf_Error(xf_VectorSetPODMass(All, SnapSet, nSnap, nBasis, BasisSet));
      if (ierr != xf_OK) return ierr;
      xf_printf("done.\n");

      // destroy SnapSet
      ierr = xf_Error(xf_DestroyVectorSet(SnapSet));
      if (ierr != xf_OK) return ierr;
      SnapSet = NULL;


      /* Write out Basis.data */
      xf_printf("\nWriting Basis VectorSet to Basis.data.\n");

      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_DataSetAdd(DataSet, "BasisSet", xfe_VectorSet,
				    xfe_True, (void *) BasisSet, &D));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, "Basis.data"));
      if (ierr != xf_OK) return ierr;
      D->Data = NULL;
      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;

    } // end if !SmallerPOD
    else{
      /* At this point, we need to load an existing Basis.data */

      // create dataset for Basis reading
      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;

      xf_printf("Reading in Basis.data\n"); fflush(stdout);

      // read in Basis.data
      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, "Basis.data", DataSet));
      if (ierr != xf_OK) return ierr;
      
      xf_printf("Done\n"); fflush(stdout);
      
      D = DataSet->Head;
      if (strncmp(D->Title, "BasisSet", 8) != 0) return xf_Error(xf_INPUT_ERROR);
      BasisSet = (xf_VectorSet *) D->Data;
      D->Data = NULL; // so we do not destroy ICSet!

      // do we have enough basis vectors?
      xf_printf("BasisSet has %d vectors.  We need nBasis = %d for the ROM.\n", 
		BasisSet->nVector, nBasis);
      if (nBasis > BasisSet->nVector)
	return xf_Error(xf_OUT_OF_BOUNDS);
      else
	xf_printf("OK.\n");

      // destroy DataSet
      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;
    }



    /*----------------------*/
    /* Compute ROM matrices */
    /*----------------------*/

  
    /** Rr = V^T * R(0) **/

    ierr = xf_Error(xf_Alloc((void **) &Rr, nBasis, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    // create residual vector for computing R(0)
    ierr = xf_Error(xf_FindSimilarVector(All, BasisSet->Vector+0, "Residual_0",
					 xfe_True, xfe_False, NULL, &R0, NULL));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_SetZeroVector(U)); // U = 0
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_CalculateResidual(All, U, R0, NULL, NULL));
    if (ierr != xf_OK) return ierr;

    // Dot R(0) with BasisSet
    for (i=0; i<nBasis; i++){
      ierr = xf_Error(xf_VectorDot(BasisSet->Vector+i, R0, Rr+i));
      if (ierr != xf_OK) return ierr;
    }
  
  

    /**  Mr = V^T * M * V  **/

    ierr = xf_Error(xf_Alloc((void **) &Mr, nBasis*nBasis, sizeof(real)));
    if (ierr != xf_OK) return ierr;

  
    for (i=0; i<nBasis; i++){

      // U = BasisSet(i)
      ierr = xf_Error(xf_SetVector(BasisSet->Vector+i, xfe_Set, U));
      if (ierr != xf_OK) return ierr;

      // U = M * U
      ierr = xf_Error(xf_MultMassMatrix(All, 1.0, U));
      if (ierr != xf_OK) return ierr;

      // Mr(i,i:nBasis) = M(i:nBasis,i) = V^T * U
      for (j=i; j<nBasis; j++){
	ierr = xf_Error(xf_VectorDot(BasisSet->Vector+j, U, Mr + i*nBasis+j));
	if (ierr != xf_OK) return ierr;
	if (i != j) Mr[j*nBasis+i] = Mr[i*nBasis+j]; // Mr is symmetric
      }
    }



    /**  Ar = V^T * R_U * V  **/

    ierr = xf_Error(xf_Alloc((void **) &Ar, nBasis*nBasis, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    // create residual vector for computing R(U)
    ierr = xf_Error(xf_FindSimilarVector(All, BasisSet->Vector+0, "Residual",
					 xfe_True, xfe_False, NULL, &R, NULL));
    if (ierr != xf_OK) return ierr;

    for (i=0; i<nBasis; i++){

      // R = R(BasisSet(i))
      ierr = xf_Error(xf_CalculateResidual(All, BasisSet->Vector+i, R, NULL, NULL));
      if (ierr != xf_OK) return ierr;

      // R -= R0
      ierr = xf_Error(xf_VectorMultSet(R0, -1.0, xfe_Add, R));
      if (ierr != xf_OK) return ierr;

      // Ar(:,i) = V^T * R
      for (j=0; j<nBasis; j++){
	ierr = xf_Error(xf_VectorDot(BasisSet->Vector+j, R, Ar + j*nBasis+i));
	if (ierr != xf_OK) return ierr;
      }

    }


    /**  For Laplace regularization   **/
    /**  Br = V^T * R_U(Laplace) * V  **/

    ierr = xf_Error(xf_Alloc((void **) &Br, nBasis*nBasis, sizeof(real)));
    if (ierr != xf_OK) return ierr;

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
    
    All->EqnSet->ResTerms = ResTermsLReg; // set to new res terms

    // get orig viscosity
    ierr = xf_Error(xf_GetKeyValueReal(All->EqnSet->KeyValue, "Viscosity", &ViscosityOrig));
    if (ierr != xf_OK) return ierr;

    // set new viscosity
    ViscosityLReg = 1.0;
    ierr = xf_Error(xf_SetKeyValueReal(All->EqnSet->KeyValue, "Viscosity", ViscosityLReg));
    if (ierr != xf_OK) return ierr;

    // R0 = R(0)
    ierr = xf_Error(xf_SetZeroVector(U)); // U = 0
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_CalculateResidual(All, U, R0, NULL, NULL));
    if (ierr != xf_OK) return ierr;

    
    for (i=0; i<nBasis; i++){

      // R = R(BasisSet(i))
      ierr = xf_Error(xf_CalculateResidual(All, BasisSet->Vector+i, R, NULL, NULL));
      if (ierr != xf_OK) return ierr;

      // R -= R0
      ierr = xf_Error(xf_VectorMultSet(R0, -1.0, xfe_Add, R));
      if (ierr != xf_OK) return ierr;

      // Br(:,i) = V^T * R
      for (j=0; j<nBasis; j++){
	ierr = xf_Error(xf_VectorDot(BasisSet->Vector+j, R, Br + j*nBasis+i));
	if (ierr != xf_OK) return ierr;
      }

    }

    All->EqnSet->ResTerms = ResTermsOrig; // back to original
    ierr = xf_Error(xf_SetKeyValueReal(All->EqnSet->KeyValue, "Viscosity", ViscosityOrig));
    if (ierr != xf_OK) return ierr;


    // Destroy ResTermsLReg
    ierr = xf_Error(xf_DestroyResTerms(ResTermsLReg));
    if (ierr != xf_OK) return ierr;


    // destroy R0
    ierr = xf_Error(xf_DestroyVector(R0, xfe_True));
    if (ierr != xf_OK) return ierr;




    /** Cr = J_U^T*V **/

    ierr = xf_Error(xf_Alloc((void **) &Cr, nOutput*nBasis, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_SetZeroVector(U)); // U = 0
    if (ierr != xf_OK) return ierr;

    // loop over outputs
    for (i=0; i<nOutput; i++){

      // calculate J_U for output i (store it in R)
      ierr = xf_Error(xf_CalculateOutput(All, All->EqnSet->Outputs->Output[i].Name, 
					 U, &val, R, xfe_Set));
      if (ierr != xf_OK) return ierr;
    
      // Dot J_U with BasisSet
      for (j=0; j<nBasis; j++){
	ierr = xf_Error(xf_VectorDot(BasisSet->Vector+j, R, Cr+i*nBasis+j));
	if (ierr != xf_OK) return ierr;
      }
    } // i
    
    
  



    // destroy R
    ierr = xf_Error(xf_DestroyVector(R, xfe_True));
    if (ierr != xf_OK) return ierr;



    /* Destroy basis vectorset */
    ierr = xf_Error(xf_DestroyVectorSet(BasisSet));
    if (ierr != xf_OK) return ierr;


    /*----------------------------*/
    /* Write ROM matrices to file */
    /*----------------------------*/
  
    if (myRank == 0){ // root does the writing

      if ((fid = fopen("ROM.m", "w")) == NULL) return xf_Error(xf_FILE_WRITE_ERROR);

      fprintf(fid, "Rr = [\n");
      for (i=0; i<nBasis; i++){
	fprintf(fid, "%.15E\n", Rr[i]);
      }
      fprintf(fid, "];\n");
    
      fprintf(fid, "Mr = [\n");
      for (i=0; i<nBasis; i++){
	for (j=0; j<nBasis; j++){
	  fprintf(fid, "%.15E ", Mr[i*nBasis+j]);
	}
	fprintf(fid, "\n"); 
      }
      fprintf(fid, "];\n"); 
    
      fprintf(fid, "Ar = [\n");
      for (i=0; i<nBasis; i++){
	for (j=0; j<nBasis; j++){
	  fprintf(fid, "%.15E ", Ar[i*nBasis+j]);
	}
	fprintf(fid, "\n"); 
      }
      fprintf(fid, "];\n"); 

      fprintf(fid, "Br = [\n");
      for (i=0; i<nBasis; i++){
	for (j=0; j<nBasis; j++){
	  fprintf(fid, "%.15E ", Br[i*nBasis+j]);
	}
	fprintf(fid, "\n"); 
      }
      fprintf(fid, "];\n"); 

      fprintf(fid, "Cr = [\n");
      for (i=0; i<nOutput; i++){
	for (j=0; j<nBasis; j++)
	  fprintf(fid, "%.15E ", Cr[i*nBasis+j]);
	fprintf(fid, "\n"); 
      }
      fprintf(fid, "];\n");

      fclose(fid);
    }


    // Free memory
    xf_Release( (void *) Rr);
    xf_Release( (void *) Mr);
    xf_Release( (void *) Ar);
    xf_Release( (void *) Br);
    xf_Release( (void *) Cr);

  } // end if PODFlag
  else xf_printf("Not taking snapshots or performing POD\n");


  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;
  
  // Release memory
  xf_Release( (void *) TimeMask);
  xf_Release( (void *) TimeFlag);
 
 
  xf_printf("xf_MRHessianIC finished.\n");

  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
