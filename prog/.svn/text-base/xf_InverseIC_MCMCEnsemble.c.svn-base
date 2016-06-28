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
  FILE:  xf_InverseIC_MCMCEnsemble.c

  This program performs probabalistic initial condition inversion
  using the ensemble Markov-chain Monte Carlo approach (strecth move)
  with adjoint-based output calculation.

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
#include "xf_MathLapack.h"

#include "xf_Arg.h"
#include "xf_InverseIC_Common.h"
#include <time.h>
#include <stdlib.h>

// Enumerated type for MCMC proposal distribution type
enum xfe_MCMCPropType{
  xfe_MCMCPropIsoBox,       // isotropic box in parameter space
  xfe_MCMCPropLast
};

static char *xfe_MCMCPropName[xfe_MCMCPropLast] = {
  "IsoBox",
};

/******************************************************************/
// FUNCTION Definition: xf_WriteInitialSampleFile
static int
xf_WriteInitialSampleFile(int nInitialSample, int nParam, real **ParamInitialSample,
		char **ParamName, real *ParamRange, real a, const char *InitialSampleFile)
{
/*
PURPOSE:  Writes initial samples to a file

INPUTS:

  nInitialSample: number of initial samples (= nWalker)
  nParam: number of parameters
  ParamInitialSample: nWalker by nParam array of subsamples
  ParamName: names of parameters
  ParamRange: range of parameters
  a: scale of domain for scaling factor Z
  IntialSampleFile: name of file to which initial samples are written
  
OUTPUTS:  None  

RETURN: Error code

*/

  int ierr = xf_OK;
  int myRank;
  int iParam, iInitialSample;
  FILE *fid;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // Open file for writing (error handled in parallel)
  ierr = xf_Error(xf_fopen(InitialSampleFile, "w", &fid));
  if (ierr != xf_OK) return ierr;

  // Only root writes
  if (myRank == 0){
    // Write out header
    fprintf(fid, "%% Initial samples \n");
    fprintf(fid, "%% nInitialSample = %d\n", nInitialSample);
    fprintf(fid, "%% a = %.8E\n", a);
    fprintf(fid, "%% nParam = %d\n", nParam);
    
    fprintf(fid, "%% ParamMin = ");    
    for (iParam=0; iParam<nParam; iParam++)
      fprintf(fid, " %.6E", ParamRange[2*iParam+0]);
    fprintf(fid, "\n");
    
    fprintf(fid, "%% ParamMax = ");
    for (iParam=0; iParam<nParam; iParam++)
      fprintf(fid, " %.6E", ParamRange[2*iParam+1]);
    fprintf(fid, "\n");
    
    fprintf(fid, "%%");
    for (iParam=0; iParam<nParam; iParam++)
      fprintf(fid, " %s", ParamName[iParam]);
    fprintf(fid, "\n");
    
    for (iInitialSample=0; iInitialSample<nInitialSample; iInitialSample++){
      for (iParam=0; iParam<nParam; iParam++)
        fprintf(fid, " %20.15E ", ParamInitialSample[iInitialSample][iParam]);
      fprintf(fid, "\n");
    } // iInitialSample
  }
  
  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_WriteSubSampleFile
static int
xf_WriteSubSampleFile(int nTotalSubSample, int nSubSample, int nParam, 
    real **ParamSubSample, char **ParamName, real *ParamRange, real a, 
    const char *SubSampleFile, int countWrite)
{
/*
PURPOSE:  Writes MCMC samples to a file

INPUTS:

  nTotalSubSample: total number of sub-samples
  nSubSample: number of sub-samples to write
  nParam: number of parameters
  ParamSubSample: nSubSample by nParam array of subsamples
  ParamName: names of parameters
  ParamRange: range of parameters
  a: scale of domain for scaling factor Z
  SubSampleFile: name of file to which sub-samples are written
  countWrite: how many times writing into a file
  
OUTPUTS:  None  

RETURN: Error code

*/

  int ierr = xf_OK;
  int myRank;
  int iParam, iSubSample;
  FILE *fid;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // Open file for writing (error handled in parallel)
  ierr = xf_Error(xf_fopen(SubSampleFile, "a", &fid));
  if (ierr != xf_OK) return ierr;

  // Only root writes
  if (myRank == 0){
    if (countWrite == 1) {
    // Write out header
    fprintf(fid, "%% MCMC sub-samples \n");
    fprintf(fid, "%% nSubSample = %d\n", nTotalSubSample);
    fprintf(fid, "%% a = %.8E\n", a);
    fprintf(fid, "%% nParam = %d\n", nParam);
    
    fprintf(fid, "%% ParamMin = ");    
    for (iParam=0; iParam<nParam; iParam++)
      fprintf(fid, " %.6E", ParamRange[2*iParam+0]);
    fprintf(fid, "\n");
    
    fprintf(fid, "%% ParamMax = ");
    for (iParam=0; iParam<nParam; iParam++)
      fprintf(fid, " %.6E", ParamRange[2*iParam+1]);
    fprintf(fid, "\n");
    
    fprintf(fid, "%%");
    for (iParam=0; iParam<nParam; iParam++)
      fprintf(fid, " %s", ParamName[iParam]);
    fprintf(fid, "\n");
    }
    
    for (iSubSample=0; iSubSample<nSubSample; iSubSample++){
      for (iParam=0; iParam<nParam; iParam++)
        fprintf(fid, " %20.15E ", ParamSubSample[iSubSample][iParam]);
      fprintf(fid, "\n");
    } // iSubSample
  }
  
  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
// FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  // Arguments read in on the command line
  char *ArgIn[] = {"job", "NULL", ".job file (time discretization, .eqn file etc.)",
		   "MCMCParamFileEnsemble", "MCMCParamFileEnsemble.txt", "name of file containing MCMC parameters",
		   "OutputFile", "OutputFile.txt", "name of file containing output observations",
                   "nWalker", "100", "number of walkers",
		   "r", "10000", "how many steps each walker take",
		   "PropType", "IsoBox", "proposal type (IsoBox/AnisoEllipse)",
		   "a", "1.01", "scale domain for scaling factor Z, a>1",
                   "fac", "1.0", "scale the initial positions of walkers with respect to the domain (0<fac<=1)",
                   "SubSampleFile", "SubSamples.txt", "file to which subsamples are written",
		   "InitialSampleFile", "InitialSamples.txt", "file to which initial samples are written",
                   "SubSampleEvery", "1000", "subsample every how many samples? (Multiple of nWalker)",
		   "WriteEvery", "10000", "how often to write subsamples into a file (Multiple of SubSampleEvery)",
		   "PrintEvery", "100000", "how often to print",
		   "\0"};

  int ierr, myRank, nProc, index;
  int count, iter;
  int iOut;
  int nOut;                                  // Number of outputs
  real a;                                    // Set the domain for scaling factor, a > 1
  int r;                                     // How many steps each walker takes (r = nSample/nWalkers)
  int k;                                     // Walker k index
  int j;                                     // Walker j index
  real Z;                                    // Scaling factor
  real ksi;                                  // CDF(z)
  real G;                                    // Variable to normalized g(Z)
  int countWrite;                            // Count how many times writing samples to a file	
  real fac;                                  // Scale the initial position of walkers with respect to the domain

  int iParam;
  int nParam;                                // Number of parameters
  char **ParamName = NULL;                   // Parameter names
  long long int iSample;
  int nWalker;                               // Number of walkers
  int nAccept;                               // Count accepted samples
  int WriteEvery;                            // Frequency of writing subsamples into a file
  int PrintEvery;                            // Frequency of printing
  int SubSampleEvery;                        // Frequency of sub-sampling
  int nSample;                               // Total samples generated
  int nSubSample;                            // Total subsamples written each time
  int nTotalSubSample;                       // Total subsamples after all iterations
  int MaxCountWrite;                         // Total number of times writing samples to a file

  enum xfe_MCMCPropType PropType;            // Proposal type	
  enum xfe_Bool *Accept;                     // Flag the accepted proposed parameter values
  real *ParamRange = NULL;                   // Range of the domain (specified in MCMCParamFileEnsemble)
  real *ParamWindowRange = NULL;             // Range of the initial positions of all walkers (specified in MCMCParamFileEnsemble)
  real **ParamSubSample = NULL;              // Subsamples (size: nSubSample x nParam)
  real **ParamValueCurrent = NULL;           // Current parameter values (size: nWalker x nParam)
  real **ParamValueNext = NULL;              // Next parameter values (size: nWalker x nParam)
  real **ParamValueInitial = NULL;           // Initial parameter values (size: nWalker x nParam)
  real *ParamMin;                            // Minimum parameter values
  real *ParamMax;                            // Maximum parameter values

  char **OutName = NULL;                     // Sensor names for each output reading
  int *OutTimeIndex = NULL;                  // Time index at which sensor reading was taken
  real *OutValue = NULL;                     // Sensor output readings
  real *OutStdErr = NULL;                    // Sensor output measurement errors (stdev)
  real *J = NULL;                            // Outputs computed at each MCMC iteration
  real *oldpow = NULL;
  real newpow, logalpha, randval;
  int WalkerSetStart;                        // Index for set of walkers
  int iSubSample;
  real *epsParam = NULL;                     // Will store parameter epsilons
  real *J_mu = NULL;                         // Stores sensitivites of outputs to params (aniso)
	
  char jobFile[xf_MAXSTRLEN]               = "NULL"; 
  char MCMCParamFileEnsemble[xf_MAXSTRLEN] = "NULL";
  char OutputFile[xf_MAXSTRLEN]            = "NULL";
  char FileName[xf_MAXSTRLEN]              = "NULL";
  char SubSampleFile[xf_MAXSTRLEN]         = "NULL";
  char InitialSampleFile[xf_MAXSTRLEN]     = "NULL";
  
  real *ParamWindowMin;                      // Minimum parameter values of initial positions
  real *ParamWindowMax;                      // Maximum parameter values of initial positions
  int countOutput;                           // Count how many times outputs were calculated
  int multiple;

  FILE *fid;
  
  xf_KeyValue KeyValueArg;
  xf_TimeHistData *TimeHistData = NULL;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Vector *U0 = NULL;                      // An initial condition state vector
  xf_Vector **J_U0 = NULL;                   // Sensitivity of outputs w.r.t ICs
  xf_SolverData *SolverData = NULL;
  xf_All *All;

  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  xf_printf("\n");
  xf_printf("=== Linear Inverse Solver for Initial Conditions using MCMC ===\n");
  xf_printf("\n");

  // Initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValueArg));
  if (ierr != xf_OK) return ierr;

  // Parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValueArg);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  /* Get jobFile */
  ierr = xf_GetKeyValue(KeyValueArg, "job", jobFile);
  if (ierr != xf_OK) return ierr;

  /* Get MCMCParamFileEnsemble*/
  ierr = xf_GetKeyValue(KeyValueArg, "MCMCParamFileEnsemble", MCMCParamFileEnsemble);
  if (ierr != xf_OK) return ierr;

  /* Get OutputFile */
  ierr = xf_GetKeyValue(KeyValueArg, "OutputFile", OutputFile);
  if (ierr != xf_OK) return ierr;

  /* How many walkers */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "nWalker", &nWalker));
  if (ierr != xf_OK) return ierr;

  /* How many steps taken by each walker */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "r", &r));
  if (ierr != xf_OK) return ierr;

  /* Approximate the initial positions of the walkers with respect to the domain */
  ierr = xf_Error(xf_GetKeyValueReal(KeyValueArg, "fac", &fac));
  if (ierr != xf_OK) return ierr;

  /* Proposal type */
  ierr = xf_Error(xf_GetKeyValueEnum(KeyValueArg, "PropType", 
    xfe_MCMCPropName, (int ) xfe_MCMCPropLast, (int *) &PropType));
  if (ierr != xf_OK) return ierr;

  /* Proposal window for scale factor Z */
  ierr = xf_Error(xf_GetKeyValueReal(KeyValueArg, "a", &a));
  if (ierr != xf_OK) return ierr;

  /* Get SubSampleFile */
  ierr = xf_GetKeyValue(KeyValueArg, "SubSampleFile", SubSampleFile);
  if (ierr != xf_OK) return ierr;

  /* Get InitialSampleFile */
  ierr = xf_GetKeyValue(KeyValueArg, "InitialSampleFile", InitialSampleFile);
  if (ierr != xf_OK) return ierr;

  /* How often to sub-sample */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "SubSampleEvery", &SubSampleEvery));
  if (ierr != xf_OK) return ierr;

  /* How often to write */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "WriteEvery", &WriteEvery));
  if (ierr != xf_OK) return ierr;

  /* How often to print */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "PrintEvery", &PrintEvery));
  if (ierr != xf_OK) return ierr;

  // Destroy key-value from arg list
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValueArg));
  if (ierr!=xf_OK) return ierr;

  // Read in OutputFile
  ierr = xf_Error(xf_ReadICOutputFile(OutputFile, &nOut, &OutName, &OutTimeIndex, 
    &OutValue, &OutStdErr));
  if (ierr != xf_OK) return ierr;

  /*---------------------------------------*/
  /* Read All from .job, including EqnSet. */
  /*---------------------------------------*/

  ierr = xf_Error(xf_ReadAllFromJobFile(jobFile, xfe_True, &All));
  if (ierr != xf_OK) return ierr;
  
  /* Load dynamic EqnSet Library  */  
  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;
  
  /* Register EqnSet and check/set default eqnset parameters. */
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;

  /*------------------------------------------------*/
  /* Check for or create output sensitivity vectors */
  /* Also, Load these vectors into memory           */
  /*------------------------------------------------*/

  ierr = xf_Error(xf_LoadSensitivityVectors(All, nOut, OutName,
     OutTimeIndex, &J_U0));
  if (ierr != xf_OK) return ierr;

  // Read in MCMCParamFileEnsemble
  ierr = xf_Error(xf_ReadICParamFileEnsemble(MCMCParamFileEnsemble, 
    &nParam, &ParamName, &ParamRange, &ParamWindowRange));
  if (ierr != xf_OK) return ierr;

  /*-----------------*/
  /* MCMC Iterations */
  /*-----------------*/

  // Calculate number of samples and subsamples
  nSample = r*nWalker;
  nTotalSubSample = nSample/SubSampleEvery*nWalker;
  nSubSample = WriteEvery;
  
  MaxCountWrite = nTotalSubSample/WriteEvery+1;
  countOutput = 0;

  // Total subsamples is always multiple of nWalker
  if ((nTotalSubSample % nWalker) != 0) {
    multiple = nTotalSubSample/nWalker;
    nTotalSubSample = multiple*nWalker;
  }

  /* Check whether a>1 or not. If not, code will be terminated.*/
  if ((a<1) || (a==1)) {
    xf_printf("a has to be greater than 1. \n");
    return xf_Error(xf_INPUT_ERROR);
  }
  else if ((fac>1) || (fac<=0)) {
    xf_printf("fac has to be positive and not greater than one. \n");
    return xf_Error(xf_INPUT_ERROR);
  }
  else if ((SubSampleEvery % nWalker) != 0) {
    xf_printf("SubSampleEvery has to be multiple of nWalker. \n");
    return xf_Error(xf_INPUT_ERROR);
  }
  else if ((WriteEvery % SubSampleEvery) != 0) {
   xf_printf("WriteEvery has to be multiple of SubSampleEvery. \n");
   return xf_Error(xf_INPUT_ERROR);
  }


  // Allocate memory for subamples
  ierr = xf_Error(xf_Alloc2( (void ***) &ParamSubSample, nSubSample, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Allocate memory for ParamValueCurrent
  ierr = xf_Error(xf_Alloc2( (void ***) &ParamValueCurrent, nWalker, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Allocate memory for ParamValueInitial
  ierr = xf_Error(xf_Alloc2( (void ***) &ParamValueInitial, nWalker, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Allocate memory for ParamValueNext
  ierr = xf_Error(xf_Alloc2( (void ***) &ParamValueNext, nWalker, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Allocate memory for outputs 
  ierr = xf_Error(xf_Alloc( (void **) &J, nOut, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // Allocate memory for oldpow
  ierr = xf_Error(xf_Alloc( (void **) &oldpow, nWalker, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Allocate memory for Accept 
  ierr = xf_Error(xf_Alloc( (void **) &Accept, nWalker, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Allocate memory for ParamMin
  ierr = xf_Error(xf_Alloc( (void **) &ParamMin, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // Allocate memory for ParamMax
  ierr = xf_Error(xf_Alloc( (void **) &ParamMax, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Allocate memory for ParamWindowMin
  ierr = xf_Error(xf_Alloc( (void **) &ParamWindowMin, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Allocate memory for ParamWindowMax
  ierr = xf_Error(xf_Alloc( (void **) &ParamWindowMax, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;



  /* Determine initial parameter values for all walkers k */
  for (iParam=0; iParam<nParam; iParam++) {
    
    // Minimum and maximum parameter values
    ParamMin[iParam] = ParamRange[2*iParam+0];
    ParamMax[iParam] = ParamRange[2*iParam+1];
  
    // Minimum and maximum initial parameter values
    ParamWindowMin[iParam] = ParamWindowRange[2*iParam+0];
    ParamWindowMax[iParam] = ParamWindowRange[2*iParam+1];
	
    for (k=0; k<nWalker; k++) {

      // Randomly pick the initial parameter values of all walkers k.
      ierr = xf_Error(xf_RandUniform(1, &randval));
      if (ierr != xf_OK) return ierr;
        
      // Calculate the next parameter value. This will make sure that not all initial
      // positions will be the same. (The algorithm will not work if all walkers start
      // at the same position)
      if (ParamMax[iParam] == 0)
	ParamValueCurrent[k][iParam] = fac*ParamMin[iParam]*randval;
      else
	ParamValueCurrent[k][iParam] = fac*(ParamMax[iParam]*randval+ParamMin[iParam]);

      // Check whether initial parameter values are within the domain.
      // If not, then repick initial parameter values.			
      while ((ParamValueCurrent[k][iParam] < ParamWindowMin[iParam]) || 
             (ParamValueCurrent[k][iParam] > ParamWindowMax[iParam])) {  
        ierr = xf_Error(xf_RandUniform(1, &randval));
        if (ierr != xf_OK) return ierr;

	if (ParamMax[iParam] == 0)
	  ParamValueCurrent[k][iParam] = fac*ParamMin[iParam]*randval;
	else
	  ParamValueCurrent[k][iParam] = fac*(ParamMax[iParam]*randval+ParamMin[iParam]);
      }

      // Initial paramater values give the first nWalker samples
      ParamValueInitial[k][iParam] = ParamValueCurrent[k][iParam];
    } // k
  } // iParam

  // Total number of samples accepted at this point
  nAccept = nWalker;

  /* Write initial parameter values to a file */
  if (SubSampleEvery > nWalker) {
    xf_printf("# Initial samples = %d \n", nWalker); 
    ierr = xf_Error(xf_WriteInitialSampleFile(nWalker, nParam, ParamValueInitial,
                                              ParamName, ParamRange, a, InitialSampleFile));
    if (ierr != xf_OK) return ierr;
    xf_printf("# Initial parameter values were written to %s. \n", InitialSampleFile);
  }

  /* Store initial values to ParamSubSample (only if no sub-sampling) */
  if (SubSampleEvery == nWalker) {
    for (iParam=0; iParam<nParam; iParam++) {
      for (k=0; k<nWalker; k++)
        ParamSubSample[k][iParam] = ParamValueInitial[k][iParam];
    }

    xf_printf("Initial sample = %d \n", nWalker);
    xf_printf("No sub-sampling \n");
  }


  /* Destroy ParamValueInitial */
  xf_Release2((void **) ParamValueInitial);

  /* Calculate outputs and oldpow for first nWalker samples */
  for (k=0; k<nWalker; k++) {
    // Calculate outputs at first nWalker samples: store in J
    ierr = xf_Error(xf_CalculateLinearOutputs(All, nParam, ParamName,
					      ParamValueCurrent[k], nOut, J_U0, epsParam, J, J_mu, &U0));
    if (ierr != xf_OK) return ierr;

    // Precalculate oldpow (see MCMC code below)
    for (iOut=0, oldpow[k]=0.; iOut<nOut; iOut++)
      oldpow[k] -= 0.5*xf_PowInt((OutValue[iOut]-J[iOut])/OutStdErr[iOut], 2);
  }

  /* Begin MCMC Iterations */
	
  iSample = nWalker;
  iSubSample = 0;
  WalkerSetStart = 1;
  count = 0;
  countWrite = 0;
	
  // No sub-sampling
  if (SubSampleEvery == nWalker) {
    iSubSample = nWalker;
    WalkerSetStart = 2;
  }	

  for (iter=1; iter<r; iter++) {

    //Reset Accept to true
    for (k=0; k<nWalker; k++)
      Accept[k] = xfe_True;

    // Generate proposal and check whether it's in the domain
    // Note: all walkers have to be moved from t to t+1 at the same time
    for (k=0; k<nWalker; k++) {

      // Pick walker j randomly (walker j cannot be the same as walker k)
      ierr = xf_Error(xf_RandUniform(1, &randval));
      if (ierr != xf_OK) return ierr;
      
      j = nWalker*randval;
      while (j == k) {
	ierr = xf_Error(xf_RandUniform(1, &randval));
	if (ierr != xf_OK) return ierr;

	j= nWalker*randval;
      }
		
      // Choose scaling factor Z independently for each walker
      ierr = xf_Error(xf_RandUniform(1, &randval));
      if (ierr != xf_OK) return ierr;

      ksi = randval;  // CDF(z)
      G = 2*sqrt(a)-2/sqrt(a);  // Variable to normalized g(Z)
      Z = 0.25*xf_PowInt((ksi*G*sqrt(a)+2),2)/a;  // Scaling factor

      // Calculate the next parameter values of each walker using stretch move
      for (iParam=0; iParam<nParam; iParam++) {
	ParamValueNext[k][iParam] = ParamValueCurrent[j][iParam] + 
	  Z*(ParamValueCurrent[k][iParam]-ParamValueCurrent[j][iParam]);
      }

      // modify oldpow to incorporate Z^(npar-1)
      oldpow[k] -= (nParam-1.)*log(Z);

      // Check whether the next parameter value is within range. If not in
      // the domain, the proposed next paramater values will be flagged 
      // as false
      for (iParam=0; iParam<nParam; iParam++) {
	if(ParamValueNext[k][iParam] < ParamMin[iParam] || 
	   ParamValueNext[k][iParam] > ParamMax[iParam]) {
	  Accept[k] = xfe_False;
	}
      }
    } // k

    // Acceptance or rejection proposal
    for (k=0; k<nWalker; k++) {

      // Print out iteration info
      if ((iSample % PrintEvery) == 0){
	xf_printf("MCMC sample %d: acceptance ratio = %.6f\n", iSample, 
		  (((real) nAccept) / ((real) iSample)) );  fflush(stdout);
      }

      // If the proposal isn't in domain, then the proposal gets rejected
      if (Accept[k] == xfe_False) {
	for (iParam=0; iParam<nParam; iParam++)
	  ParamValueNext[k][iParam] = ParamValueCurrent[k][iParam];
        countOutput++;
      }
      // If the proposal is in domain, then the proposal will be evaluated
      else {

	// Evaluate outputs at munext; store in Jnew
	ierr = xf_Error(xf_CalculateLinearOutputs(All, nParam, ParamName,
						  ParamValueNext[k], nOut, J_U0, epsParam, J, J_mu, &U0));
	if (ierr != xf_OK) return ierr;

	// Calculate newpow
	for (iOut=0, newpow=0; iOut<nOut; iOut++)
	  newpow -= 0.5*xf_PowInt((OutValue[iOut]-J[iOut])/OutStdErr[iOut], 2);

	// Calculate logaplha
	logalpha = newpow - oldpow[k];

	// Acceptance step: flip MCMC coin
	ierr = xf_Error(xf_RandUniform(1, &randval));
	if (ierr != xf_OK) return ierr;

	if (log(randval) < logalpha){
	  nAccept++;
	  oldpow[k] = newpow;
	} 
	else {
	  for (iParam=0; iParam<nParam; iParam++)
	    ParamValueNext[k][iParam] = ParamValueCurrent[k][iParam];
	}
      } // End of if statement for checking proposal in domain or not
			
        // Update ParamSubSample
      if (iSample >= WalkerSetStart*SubSampleEvery-nWalker && 
	  iSample <= WalkerSetStart*SubSampleEvery-1) {
	for (iParam=0; iParam<nParam; iParam++)
	  ParamSubSample[iSubSample][iParam] = ParamValueNext[k][iParam];
    
	iSubSample++;					
	count++; //Keep track of stored samples
					
	// Reset count and WalkerSetStart
	if (count == nWalker)
	  WalkerSetStart++;

      } //iSample

      // Update ParamValueCurrent
      for (iParam=0; iParam<nParam; iParam++)  
	ParamValueCurrent[k][iParam] = ParamValueNext[k][iParam];
  
      iSample++;
			
      // Write sub-samples
      if ((nTotalSubSample < WriteEvery) && (iSample == nSample)) {
	xf_printf("#sub-samples %d: writing sub-samples to %s\n",
		  iSubSample, SubSampleFile);

	ierr = xf_Error(xf_WriteSubSampleFile(nTotalSubSample, iSubSample, nParam,
					      ParamSubSample, ParamName, ParamRange, a, SubSampleFile, countWrite));  
	if (ierr != xf_OK) return ierr;
      }

      if ((nTotalSubSample % WriteEvery) == 0 && (iSubSample % WriteEvery) == 0 
	  && (iSubSample != 0)) {
	countWrite++;
	xf_printf("# sub-samples %d: writing sub-samples to %s\n", 
		  countWrite*iSubSample, SubSampleFile); 
 
	ierr = xf_Error(xf_WriteSubSampleFile(nTotalSubSample, iSubSample, nParam, 
					      ParamSubSample, ParamName, ParamRange, a, SubSampleFile, countWrite));
	if (ierr != xf_OK) return ierr;
         
	iSubSample = 0;
      }

      if ((nTotalSubSample % WriteEvery) != 0) {
	countWrite++;
         
	if (countWrite == MaxCountWrite) {
	  xf_printf("# sub-samples %d: writing sub-samples to %s\n", 
		    (nTotalSubSample - (countWrite-1)*WriteEvery), SubSampleFile); 
 
	  ierr = xf_Error(xf_WriteSubSampleFile(nTotalSubSample, 
						(nTotalSubSample - (countWrite-1)*WriteEvery), nParam, 
						ParamSubSample, ParamName, ParamRange, a, SubSampleFile, countWrite));
	  if (ierr != xf_OK) return ierr;
	}

	if ((iSubSample % WriteEvery) == 0 && (iSubSample != 0)) {
	  xf_printf("# sub-samples %d: writing sub-samples to %s\n", 
		    countWrite*iSubSample, SubSampleFile); 
 
	  ierr = xf_Error(xf_WriteSubSampleFile(nTotalSubSample, iSubSample, nParam, 
						ParamSubSample, ParamName, ParamRange, a, SubSampleFile, countWrite));
	  if (ierr != xf_OK) return ierr;
         
	  iSubSample = 0;
	}

      }

    } // k
		
    count = 0;
  } // iter

    /* Acceptance rate */
  xf_printf("MCMC sample %d: acceptance ratio = %.6f\n", iSample,
	    (((real) nAccept) / ((real) iSample)) );

  /* Number of outputs calculated */
  xf_printf("Outputs were calculated %d times. \n", (nSample - countOutput));

  // Print out info on chain
  xf_printf("MCMC chain finished.\n");

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
  xf_Release2((void **) ParamName);
  xf_Release2((void **) ParamValueNext);
  xf_Release2((void **) ParamValueCurrent);
  xf_Release2((void **) ParamSubSample);
  xf_Release2((void **) OutName);
  xf_Release((void  *) ParamRange);
  xf_Release((void  *) ParamWindowRange);
  xf_Release((void  *) OutTimeIndex);
  xf_Release((void  *) OutValue);
  xf_Release((void  *) OutStdErr);
  xf_Release((void  *) J);
  xf_Release((void  *) Accept);
  xf_Release((void  *) ParamMin);
  xf_Release((void  *) ParamMax);
  xf_Release((void  *) oldpow);
  xf_Release((void  *) ParamWindowMax);
  xf_Release((void  *) ParamWindowMin);


  xf_printf("xf_InverseIC_MCMC finished.\n");

  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
  
}
