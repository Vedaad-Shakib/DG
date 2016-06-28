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
  FILE:  xf_InverseIC_MCMC.c

  This program performs probabalistic initial condition inversion
  using the Markov-chain Monte Carlo approach with adjoint-based
  output calculation.

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

// enumerated type for MCMC proposal distribution type
enum xfe_MCMCPropType{
  xfe_MCMCPropIsoBox,       // isotropic box in parameter space
  xfe_MCMCPropAnisoEllipse, // anisotropic ellipse (uses derivative info)
  xfe_MCMCPropLast
};

static char *xfe_MCMCPropName[xfe_MCMCPropLast] = {
  "IsoBox",
  "AnisoEllipse",
};

/******************************************************************/
//   FUNCTION Definition: xf_WriteSampleFile
static int
xf_WriteSampleFile(int nSample, int nParam, real **ParamSample,
		   char **ParamName, real *ParamRange, 
		   enum xfe_MCMCPropType PropType, real PropScale,
		   real PropDelta,  const char *SampleFile)
{
/*
PURPOSE:  Writes MCMC samples to a file

INPUTS:

  nSample: number of samples
  nParam: number of parameters
  ParamSample: nSample by nParam array of samples
  ParamName: names of parameters
  ParamRange: range of parameters
  PropType: proposal type
  PropScale: proposal scale
  PropDelta: proposal delta
  SampleFile: name of file to which samples are written
  
OUTPUTS:  None  

RETURN: Error code

*/

  int ierr = xf_OK;
  int myRank;
  int iParam, iSample;
  FILE *fid;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // open file for writing (error handled in parallel)
  ierr = xf_Error(xf_fopen(SampleFile, "w", &fid));
  if (ierr != xf_OK) return ierr;

  // only root writes
  if (myRank == 0){
    // write out header
    fprintf(fid, "%% MCMC samples \n");
    fprintf(fid, "%% nSample = %d\n", nSample);
    fprintf(fid, "%% PropType = %s\n", xfe_MCMCPropName[PropType]);
    fprintf(fid, "%% PropScale = %.8E\n", PropScale);
    fprintf(fid, "%% PropDelta = %.8E\n", PropScale);
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
    for (iSample=0; iSample<nSample; iSample++){
      for (iParam=0; iParam<nParam; iParam++)
	fprintf(fid, " %20.15E ", ParamSample[iSample][iParam]);
      fprintf(fid, "\n");
    } // iSample
  }
  
  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SensorCoverage
static int
xf_SensorCoverage(xf_All *All, int nOut, char **OutName, int *OutTimeIndex,
		  xf_Vector **Psij)
{
  /*
    PURPOSE: Finds the number of sensors that give the maximum coverage

    INPUTS: 
      All : all structure
      nOut: number of outputs
      OutName: output names
      OutTimeIndex : output time indices
      Psij: sensitivity vectors (basically, Psi(t=0))
  
    OUTPUTS:  Maximum coverage given by the sensors picked  

    RETURN: Error code

  */
	
  real ElemVol = 0;
  int *PickedFlag = NULL;
  int nPicked = 0;
  int j, k;
  int jmax = 0;
  enum xfe_ShapeType Shape;
  real MaxCoverage = -1;
  real *Coverage = NULL;
  real *AvgPsi = NULL;
  real average = 0.;
  real rmax;
  real xref[3], xglob[3];
  int ierr;
  int egrp;
  int elem;
  int *Flag = NULL;
  xf_Vector *ElemFlag;
  xf_Vector *Psi = NULL;
  xf_Vector *EG;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  //Initialize the PickedFlag[j]=0 for all outputs 
  ierr = xf_Error(xf_Alloc( (void **) &PickedFlag, nOut, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (j=0; j<nOut ; j++) PickedFlag[j]=-1; // -1 indicates not picked

  //Intialize ElemFlag=0 for all elements (1 integer per element)
  ierr = xf_Error(xf_FindVector(All, "ElemFlag", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_False, xfe_True,
				NULL, &ElemFlag, NULL));
  if (ierr != xf_OK) return ierr;

  // set ElemFlag to -1
  ierr = xf_Error(xf_SetConstVector( ElemFlag, -1, -1.0));
  if (ierr != xf_OK) return ierr;

  //Find element geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  if (ierr != xf_OK) return ierr;
 
  // Allocate Coverage
  ierr = xf_Error(xf_Alloc( (void **) &Coverage, nOut, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  // Allocate AvgPsi
  ierr = xf_Error(xf_Alloc( (void **) &AvgPsi, nOut, sizeof(real)));
  if (ierr != xf_OK) return ierr;


  // Order the outputs (via PickedFlag) in terms of coverage they provide
  while(nPicked < nOut){

    /* Zero out coverage */
    for (j=0; j<nOut ; j++) Coverage[j]=0.;

    // Loop over elements in mesh
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

	// ElemFlag != -1 means elem is already covered
	Flag = ElemFlag->GenArray[egrp].iValue[elem];			
	if (Flag[0] != -1) continue;
			
	// Loop over outputs
	for (j=0; j<nOut; j++) {

	  // initialize AvgPsi
	  AvgPsi[j] = 0.;
	  
	  // do not consider outputs that are already picked
	  if (PickedFlag[j] >= 0) continue;

	  //Pick off Psi_j from input vector of Psi's
	  Psi = Psij[j];
	  
	  // Calculate element average of abs(Psi_j)
	  average=0.;
	  for (k=0; k<Psi->GenArray[egrp].r; k++)
	    average += fabs(Psi->GenArray[egrp].rValue[elem][k]);
	  average /= ((real) Psi->GenArray[egrp].r);

	  // store average value
	  AvgPsi[j] = average;

	} // j

	// pick output of max average Psi value on elem
	rmax = -1.0;
	jmax = -1;
	for (j=0; j<nOut; j++){
	  if (PickedFlag[j] >= 0) continue;
	  if (AvgPsi[j] > rmax){
	    jmax = j;
	    rmax = AvgPsi[j];
	  }
	}
	if (jmax < 0) return xf_Error(xf_CODE_LOGIC_ERROR);

	// element volume
	ElemVol = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];

	// this element now belongs to output jmax
	ElemFlag->GenArray[egrp].iValue[elem][0] = jmax;

	// elem centroid
	// determine element Shape 
	ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ShapeCentroid(Shape, xref));
	if (ierr != xf_OK) return ierr;
  
	ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True,
					1, xref, xglob));
	if (ierr != xf_OK) return ierr;	


	// add coverage to coverage list; include AvgPsi through rmax
	// idea: coverage depends on strength*area
	Coverage[jmax] += ElemVol*rmax; //exp(-xglob[0]);

      } // elem
    } // egrp

        
    // reduce-sum coverage over all processors
    ierr = xf_Error(xf_MPI_Allreduce(Coverage, nOut, xfe_SizeReal, xfe_MPI_SUM));
    if (ierr != xf_OK) return ierr;


    // pick output of maximum coverage
    MaxCoverage = 0.;
    jmax = -1;
    for (j=0; j<nOut; j++){
      if (PickedFlag[j] >= 0) continue; // output j already picked
      if (Coverage[j]>MaxCoverage){					
	MaxCoverage=Coverage[j];
	jmax = j;
      }
    }
    if ((jmax < 0) || (MaxCoverage == 0.)){
      xf_printf("-- Outputs remain, but none provide any more coverage --\n");
      break;
    }

    //Set PickedFlag[jmax] = nPicked and add nPicked by 1
    PickedFlag[jmax] = nPicked; // low nPicked means high coverage
    nPicked++;
    xf_printf("Picked output %3d ( Name = %5s, TimeIndex = %4d) (nPicked = %3d)\n", jmax, 
	      OutName[jmax], ((OutTimeIndex == NULL) ? -1 : OutTimeIndex[jmax]), nPicked);

    // Loop over elements, flag ones associated with jmax
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	Flag = ElemFlag->GenArray[egrp].iValue[elem];

	// Flag < -1 means elem was chosen on a prior round
	if (Flag[0] < -1) continue;
			
	// if element does not belong to output jmax, make it available (uncovered)
	if (Flag[0] != jmax) Flag[0] = -1;
	// else make it chosen on this round (not touched again)
	else Flag[0] = -2-jmax; // always < -1

      } // end for elem
    } // end for egrp
      
  } // end while

  // reset ElemFlag for writing
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      Flag = ElemFlag->GenArray[egrp].iValue[elem];
      Flag[0] += 2;
      Flag[0] *= -1;
    }
  }
  
  // Write out ElemFlag
  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ElemFlag", ElemFlag, "ElemFlag.data"));
  if (ierr != xf_OK) return ierr;


  xf_Release((void *) PickedFlag);
  xf_Release((void *) Coverage);
  xf_Release((void *) AvgPsi);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SensorCoverageWrapper
static int
xf_SensorCoverageWrapper(xf_All *All, int nOut, char **OutName, int *OutTimeIndex,
			 xf_Vector **Psij, enum xfe_Bool Consolidate, real DecayRate)
{
  /*
    PURPOSE: Wrapper for xf_SensorCoverage

    INPUTS: 
      All : all structure
      nOut: number of outputs
      OutName: output names
      OutTimeIndex : output time indices
      Psij: sensitivity vectors (basically, Psi(t=0))
      Consolidate: True to consolidate outputs by name
      DecayRate : indicates time value of information (negative means earlier is better)
  
    OUTPUTS:  Maximum coverage given by the sensors picked  

    RETURN: Error code

  */

  int ierr;
  int iter, j, j0;
  int *ConsFlag = NULL;
  int nTimeStep;
  real decay;

  if (Consolidate == xfe_False){
    // no consolidation, standard call
    ierr = xf_Error(xf_SensorCoverage(All, nOut, OutName, OutTimeIndex, Psij));
    if (ierr != xf_OK) return ierr;
  }
  else{ // consolidate by name
    
    xf_printf("Consolidating outputs by name.\n");
    
    // get number of time steps for application of decay
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nTimeStep", &nTimeStep));
    if (ierr != xf_OK) return ierr;

    // Allocate ConsFlag
    ierr = xf_Error(xf_Alloc( (void **) &ConsFlag, nOut, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (j=0; j<nOut; j++) ConsFlag[j] = 0;

    for (iter=0; iter<nOut+1; iter++){
      // pick first output not yet consolidated
      for (j0=iter; j0<nOut; j0++) if (ConsFlag[j0] == 0) break;
      // no more outputs means we are done
      if (j0 >= nOut) break;
      
      // flag output j0 as consolidated
      ConsFlag[j0] = 1;
      
      // use absolute values of sensitivity vectors
      ierr = xf_Error(xf_VectorAbs(Psij[j0]));
      if (ierr != xf_OK) return ierr;

      // decay value depends on output time index
      decay = exp(DecayRate*((real)OutTimeIndex[j0])/((real)nTimeStep));

      // multiply sensitivity by decay
      ierr = xf_Error(xf_VectorMult(Psij[j0], decay));
      if (ierr != xf_OK) return ierr;
      //xf_printf("j0 = %d, decay = %.10E\n", j0, decay);

      // consolidate outputs of same name
      for (j=j0+1; j<nOut; j++){
	if (strcmp(OutName[j0],OutName[j]) == 0){

	  if (ConsFlag[j] != 0) return xf_Error(xf_CODE_LOGIC_ERROR);

	  ConsFlag[j] = 1;
	  
	  // use absolute values of sensitivity vectors
	  ierr = xf_Error(xf_VectorAbs(Psij[j]));
	  if (ierr != xf_OK) return ierr;
      
	  // decay value depends on output time index
	  decay = exp(DecayRate*((real)OutTimeIndex[j])/((real)nTimeStep));

	  //xf_printf("  j = %d, decay = %.10E\n", j, decay);

	  // add to j0 sensitivity
	  ierr = xf_Error(xf_VectorMultSet(Psij[j], decay, xfe_Add, Psij[j0]));
	  if (ierr != xf_OK) return ierr;
	}
      }

      // store total psi and name in location iter if not already there
      if (j0 != iter){
	ierr = xf_Error(xf_SetVector(Psij[j0], xfe_Set, Psij[iter]));
	if (ierr != xf_OK) return ierr;
	strcpy(OutName[iter], OutName[j0]);
      }

    }
    if (iter > nOut){
      xf_printf("Bug in consolidation!\n"); 
      return xf_Error(xf_CODE_LOGIC_ERROR);
    }


    xf_printf("Number of consolidated outputs (i.e. sensors) = %d\n", iter);
    for (j=0; j<iter; j++)
      xf_printf("  Sensor %3d: name = %s\n", j, OutName[j]); 
    xf_printf("\n");

    // call xf_SensorCoverage
    ierr = xf_Error(xf_SensorCoverage(All, iter, OutName, NULL, Psij));
    if (ierr != xf_OK) return ierr;

  }

  xf_Release( (void **) ConsFlag);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  // arguments read in on the command line
  char *ArgIn[] = {"job", "NULL", ".job file (time discretization, .eqn file etc.)",
		   "MCMCParamFile", "MCMCParamFile.txt", "name of file containing MCMC parameters",
		   "OutputFile", "OutputFile.txt", "name of file containing output observations",
		   "nSample", "1000", "number of MCMC samples to take",
		   "PropType", "IsoBox", "proposal type (IsoBox/AnisoEllipse)",
		   "PropScale", "0.1", "proposal window scale factor for IsoBox",
		   "PropDelta", "1.0", "proposal delta for AnisoEllipse",
		   "PropEps", "0.001", "PropEps*PropScale*range = epsilon for param FD in AnisoEllipse",
		   "SampleFile", "Samples.txt", "file to which samples are written",
		   "PrintEvery", "100", "print info every how many samples?",
		   "WriteEvery", "10000", "write Sample file every how many samples?",
		   "nBruteForce", "0", "if > 1, samples evaluated on n^param grid",
		   "OrderOut", "False", "if True, outputs will be ordered in terms of coverage",
		   "Consolidate", "True", "True to consolidate outputs by name into sensors",
		   "DecayRate", "-1.0", "Reflects time value of information for output consolidation",
		   "\0"};
  int ierr, myRank, nProc;
  int i, j, k, l;
  int nAccept;
  int iSample;
  int nSample;               // number of samples in MCMC
  int PrintEvery = 1;        // printing interval during MCMC chain
  int WriteEvery = 10000;    // file writing interval during MCMC chain
  int nBruteForce = 0;       // must be greater than 1 to initiate brute force sampling
  enum xfe_Bool Accept;
  enum xfe_Bool OrderOut = xfe_False;

  int iParam;
  int nParam;                // number of parameters for MCMC
  char **ParamName = NULL;   // parameter names
  real *ParamRange = NULL;   // min/max ranges of MCMC parameters
  real *ParamValue = NULL;   // current parameter values
  real PropScale;            // proposal window scale factor
  real PropDelta;            // proposal delta
  real PropEps;              // finite differencing epsilon factor for AnisoEllipse
  real **ParamSample = NULL; // parameter samples
  real *epsParam = NULL;     // will store parameter epsilons
  real *muprev, *munext;
  real oldpow, newpow, logalpha, randval, lnqxy, lnqyx;

  int iOut;
  int nOut;                 // number of outputs
  char **OutName = NULL;    // sensor names for each output reading
  int *OutTimeIndex = NULL; // time index at which sensor reading was taken
  real *OutValue = NULL;    // sensor output readings
  real *OutStdErr = NULL;   // sensor output measurement errors (stdev)
  real *J = NULL;           // outputs computed at each MCMC iteration
  real *J_mu = NULL;        // stores sensitivites of outputs to params (aniso)
  real *H = NULL;           // stores Hessian matrix (aniso)
  real *oldH = NULL, *newH = NULL, *oldR = NULL, *tempH = NULL;
  real *w = NULL, *v = NULL;// temporary storage for computing Hessian
  real dParam;              // used in param range for proposal calculation
  enum xfe_Bool Consolidate = xfe_False;
  real DecayRate;
  int countOutput;          // how many times outputs were calculated
  
  char jobFile[xf_MAXSTRLEN]       = "NULL"; 
  char MCMCParamFile[xf_MAXSTRLEN] = "NULL";
  char OutputFile[xf_MAXSTRLEN]    = "NULL";
  char FileName[xf_MAXSTRLEN]      = "NULL";
  char SampleFile[xf_MAXSTRLEN]    = "NULL";
  enum xfe_MCMCPropType PropType; // proposal type
  FILE *fid;
  xf_KeyValue KeyValueArg;
  xf_TimeHistData *TimeHistData = NULL;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Vector *U0 = NULL;  // an initial condition state vector
  xf_Vector **J_U0 = NULL;         // sensitivity of outputs w.r.t ICs
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

  /* Get MCMCParamFile*/
  ierr = xf_GetKeyValue(KeyValueArg, "MCMCParamFile", MCMCParamFile);
  if (ierr != xf_OK) return ierr;

  /* Get OutputFile */
  ierr = xf_GetKeyValue(KeyValueArg, "OutputFile", OutputFile);
  if (ierr != xf_OK) return ierr;

  /* Get SampleFile */
  ierr = xf_GetKeyValue(KeyValueArg, "SampleFile", SampleFile);
  if (ierr != xf_OK) return ierr;
  
  /* Number of samples */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "nSample", &nSample));
  if (ierr != xf_OK) return ierr;

  /* How often to print */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "PrintEvery", &PrintEvery));
  if (ierr != xf_OK) return ierr;

  /* How often to write */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "WriteEvery", &WriteEvery));
  if (ierr != xf_OK) return ierr;

  /* Brute Force sampling? */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValueArg, "nBruteForce", &nBruteForce));
  if (ierr != xf_OK) return ierr;

  /* Proposal type */
  ierr = xf_Error(xf_GetKeyValueEnum(KeyValueArg, "PropType", 
				     xfe_MCMCPropName, (int ) xfe_MCMCPropLast, 
				     (int *) &PropType));
  if (ierr != xf_OK) return ierr;

  /* Proposal window scale factor */
  ierr = xf_Error(xf_GetKeyValueReal(KeyValueArg, "PropScale", &PropScale));
  if (ierr != xf_OK) return ierr;

  /* Proposal delta */
  ierr = xf_Error(xf_GetKeyValueReal(KeyValueArg, "PropDelta", &PropDelta));
  if (ierr != xf_OK) return ierr;

  /* Proposal epsilon scale factor */
  ierr = xf_Error(xf_GetKeyValueReal(KeyValueArg, "PropEps", &PropEps));
  if (ierr != xf_OK) return ierr;

  // Get OrderOut
  ierr = xf_GetKeyValueBool(KeyValueArg, "OrderOut", &OrderOut);
  if (ierr != xf_OK) return ierr;

  /* Consolidate */
  ierr = xf_GetKeyValueBool(KeyValueArg, "Consolidate", &Consolidate);
  if (ierr != xf_OK) return ierr;

  /* DecayRate */
  ierr = xf_Error(xf_GetKeyValueReal(KeyValueArg, "DecayRate", &DecayRate));
  if (ierr != xf_OK) return ierr;

  // destroy key-value from arg list
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

  ierr = xf_Error(xf_LoadSensitivityVectors(All, nOut, OutName, OutTimeIndex, &J_U0));
  if (ierr != xf_OK) return ierr;

  /*-------------------------------------------------------*/
  /* Are we down-selecting set of outputs using coverage?  */
  /* If so, call xf_SensorCoverage                         */
  /*-------------------------------------------------------*/

  if (OrderOut){
    xf_printf("Finding the optimum number of sensors.\n");
    xf_printf("Will print out list of sensors in order of coverage.\n");
    ierr = xf_Error(xf_SensorCoverageWrapper(All, nOut, OutName, OutTimeIndex, J_U0, 
					     Consolidate, DecayRate));
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
    xf_Release2( (void **) OutName);
    xf_Release(  (void  *) OutTimeIndex);
    xf_Release(  (void  *) OutValue);
    xf_Release(  (void  *) OutStdErr);
    /* MPI finalize (no effect in serial) */
    ierr = xf_Error(xf_MPI_Finalize());
    if (ierr != xf_OK) return ierr;
    return xf_OK; // done here
  }


  // Read in MCMCParamFile
  ierr = xf_Error(xf_ReadICParamFile(MCMCParamFile, &nParam, &ParamName, &ParamRange, 
				     &ParamValue));
  if (ierr != xf_OK) return ierr;

  // are we doing brute force?  If so, alter number of samples
  if (nBruteForce > 1) nSample = xf_PowInt(nBruteForce,nParam);
  

  /*-------------------*/
  /*  MCMC iterations  */
  /*-------------------*/

  // allocate memory for samples
  ierr = xf_Error(xf_Alloc2( (void ***) &ParamSample, nSample, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // memory for outputs 
  ierr = xf_Error(xf_Alloc( (void **) &J, nOut, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // create finite-differnce epsilons, J_mu, H, w (aniso)
  if (PropType == xfe_MCMCPropAnisoEllipse){
    ierr = xf_Error(xf_Alloc( (void **) &epsParam, nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (iParam=0; iParam<nParam; iParam++){
      epsParam[iParam] = PropScale*PropEps*(ParamRange[2*iParam+1] - ParamRange[2*iParam]);
      if (epsParam[iParam] <= 0) return xf_Error(xf_OUT_OF_BOUNDS);
    }
    ierr = xf_Error(xf_Alloc( (void **) &J_mu, nOut*nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &H, 3*nParam*nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    oldH = H; newH = H+nParam*nParam; oldR = H+2*nParam*nParam;
    ierr = xf_Error(xf_Alloc( (void **) &w, 2*max(nOut,nParam), sizeof(real)));
    if (ierr != xf_OK) return ierr;
    v = w + max(nOut,nParam);
  }

  countOutput = 0;
  
  // are we doing brute force?  If so, pregenerate samples
  if (nBruteForce > 1){
    xf_printf("Running brute force.\n");
    xf_printf("Will print out samples and posterior probability every iteration.\n");
    for (iSample=0; iSample<nSample; iSample++){
      for (iParam=0, k=iSample; iParam<nParam; iParam++){
	for (i=0, l=1; i<nParam-iParam-1; i++) l *= nBruteForce; 
	j = k/l;
	ParamSample[iSample][iParam] = ParamRange[2*iParam] + 
	  ((real) j) / ((real) nBruteForce - 1.0) *(ParamRange[2*iParam+1] - ParamRange[2*iParam]);
	k = k%l;
      } // iParam
    } // iSample
    nAccept = 0;
  }
  else{
    // initial paramater values give first sample
    for (iParam=0; iParam<nParam; iParam++) ParamSample[0][iParam] = ParamValue[iParam];
    
    // Calculate outputs at first sample: store in J
    ierr = xf_Error(xf_CalculateLinearOutputs(All, nParam, ParamName, ParamSample[0], 
					      nOut, J_U0, epsParam, J, J_mu, &U0));
    if (ierr != xf_OK) return ierr;
    
    // precalculate oldpow (see MCMC code below)
    for (iOut=0, oldpow=0.; iOut<nOut; iOut++)
      oldpow -= 0.5*xf_PowInt((OutValue[iOut]-J[iOut])/OutStdErr[iOut], 2);
    nAccept = 1;

    // calculate and decompose derivative matrix if anisotropic
    if (PropType == xfe_MCMCPropAnisoEllipse){
      // oldH = 0.5/sigma^2 * J_mu' * J_mu
      for (iOut=0; iOut<nOut; iOut++) 
	w[iOut] = 0.5/(OutStdErr[iOut]*OutStdErr[iOut]);
      xf_MTxwM_Set(J_mu, w, J_mu, nParam, nOut, nParam, oldH);
      // CholDecomp(oldH) -> oldR upper triangular part is used (U)
      for (i=0; i<nParam*nParam; i++) oldR[i] = oldH[i];
      ierr = xf_Error(xf_CholDecomp(nParam, oldR, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
  }
    
  // begin MCMC iterations
  for (iSample=nAccept; iSample<nSample; iSample++){
    
    // print out iteration info
    if ((nBruteForce <= 1) && ((iSample % PrintEvery) == 0)){
      xf_printf("MCMC sample %d: acceptance ratio = %.6f\n", iSample, 
		(((real) nAccept) / ((real) iSample)) );  fflush(stdout);
    }

    // pointer to sample
    munext = ParamSample[iSample  ];

    if (nBruteForce > 1){
      Accept = xfe_True;
    }
    else{
      // pointer to previous sample
      muprev = ParamSample[iSample-1];
      
      Accept = xfe_True;

      // generate proposal (store in munext); also check if in domain
      if (PropType == xfe_MCMCPropIsoBox){
	for (iParam=0; iParam<nParam; iParam++){
	  dParam = ParamRange[2*iParam+1] - ParamRange[2*iParam];
	  ierr = xf_Error(xf_RandUniform(1, &randval));
	  if (ierr != xf_OK) return ierr;
	  munext[iParam] = muprev[iParam] + PropScale*dParam*0.5*(2*randval-1);
	} // iParam
      }
      else if (PropType == xfe_MCMCPropAnisoEllipse){
	// Solve oldR*w=randn(nParam,1)
	ierr = xf_Error(xf_RandNormal(nParam, w));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SolveU(oldR, nParam, w));
	if (ierr != xf_OK) return ierr;
	// munext = muprev + w
	for (iParam=0; iParam<nParam; iParam++)
	  munext[iParam] = muprev[iParam] + PropDelta*w[iParam];
	/* 	// iso clip */
	/* 	for (iParam=0; iParam<nParam; iParam++){ */
	/* 	  dParam = 0.5*PropScale*(ParamRange[2*iParam+1] - ParamRange[2*iParam]); */
	/* 	  munext[iParam] = min(munext[iParam], muprev[iParam] + dParam); */
	/* 	  munext[iParam] = max(munext[iParam], muprev[iParam] - dParam); */
	/* 	} */
      }
      else return xf_Error(xf_NOT_SUPPORTED);
      // check if munext is in parameter domain
      for (iParam=0; ((iParam<nParam) && Accept); iParam++){
	if ((munext[iParam] < ParamRange[2*iParam  ]) || 
	    (munext[iParam] > ParamRange[2*iParam+1]) ) {
          Accept = xfe_False;
        }
      } // iParam
    }

    if(!Accept) {
         countOutput++;
    }

    if (Accept){
      // the proposal is within the parameter domain

      // Evaluate outputs at munext; store in Jnew
      ierr = xf_Error(xf_CalculateLinearOutputs(All, nParam, ParamName, munext, 
						nOut, J_U0, epsParam, J, J_mu, &U0));
      if (ierr != xf_OK) return ierr;

      /*
	Probability of accepting proposal is:
      
	alpha = p(y)*q(y,x)/(p(x)*q(x,y))
	= p(y)/p(x)   for IsoBox proposal
	
	p(x) ~ exp(oldpow)
	p(y) ~ exp(newpow)

	x = muprev
	y = munext

	oldpow = - 0.5 * sum((OutValue[iOut] - Jold[iOut]).^2)/sigma[iOut]^2;
	newpow = - 0.5 * sum((OutValue[iOut] - Jnew[iOut]).^2)/sigma[iOut]^2;

	In the isotropic case, q(x,y) = q(y,x).
	In the anisotropic case ...
	q(x,y) ~ exp(-.5/delta^2 * (x-y)^T * H(x) * (x-y) )
	q(y,x) ~ exp(-.5/delta^2 * (y-x)^T * H(y) * (y-x) )
	(delta = PropDelta)
      */
      
      for (iOut=0, newpow=0.; iOut<nOut; iOut++)
	newpow -= 0.5*xf_PowInt((OutValue[iOut]-J[iOut])/OutStdErr[iOut], 2);
	
      if (nBruteForce > 1){
	// print out sample and power
	for (iParam=0; iParam<nParam; iParam++)
	  xf_printf("%.10E ", munext[iParam]);
	xf_printf("%.10E\n", newpow);
      }
      else{

	if (PropType == xfe_MCMCPropIsoBox){ // isotropic case
	  logalpha = newpow - oldpow;
	}
	else{ // anisotropic case
	  // newH = 0.5/sigma^2 * J_mu' * J_mu
	  for (iOut=0; iOut<nOut; iOut++) 
	    w[iOut] = 0.5/(OutStdErr[iOut]*OutStdErr[iOut]);
	  xf_MTxwM_Set(J_mu, w, J_mu, nParam, nOut, nParam, newH);
	  // w = x-y = muprev - munext
	  for (iParam=0; iParam<nParam; iParam++) w[iParam] = muprev[iParam] - munext[iParam];
	  // v = oldH*w
	  xf_MxV_Set(oldH, w, nParam, nParam, v);
	  // lnqxy = w^T * v
	  xf_DotProduct(w, v, nParam, &lnqxy);
	  // lnqxy *= (-0.5/delta^2)
	  lnqxy *= (-0.5/(PropDelta*PropDelta));
	  // w = y-x = munex - muprev
	  for (iParam=0; iParam<nParam; iParam++) w[iParam] = munext[iParam] - muprev[iParam];
	  // v = newH*w
	  xf_MxV_Set(newH, w, nParam, nParam, v);
	  // lnqyx = w^T * v
	  xf_DotProduct(w, v, nParam, &lnqyx);
	  // lnqyx *= (-0.5/delta^2)
	  lnqyx *= (-0.5/(PropDelta*PropDelta));
	  // compute logalpha
	  logalpha = newpow - oldpow + lnqyx - lnqxy;
	}
	

	// Acceptance step: flip MCMC coin
	ierr = xf_Error(xf_RandUniform(1, &randval));
	if (ierr != xf_OK) return ierr;
	if (log(randval) < logalpha){
	  nAccept++;
	  // proposal already stored in correct place
	  oldpow = newpow;
	  if (PropType == xfe_MCMCPropAnisoEllipse){
	    swap(oldH, newH, tempH); // oldH now contains newH data
	    // CholDecomp(oldH) -> oldR upper triangular part is used (U)
	    for (i=0; i<nParam*nParam; i++) oldR[i] = oldH[i];
	    ierr = xf_Error(xf_CholDecomp(nParam, oldR, xfe_True));
	    if (ierr != xf_OK) return ierr;
	  }
	}
	else Accept = xfe_False;

	
	/* for (iParam=0; iParam<nParam; iParam++) xf_printf(" %.6E", munext[iParam]); */
	/* 	xf_printf(" %% %d \n", Accept); */
      }
    }
    
    if (!Accept){
      // did not accept; use previous sample as next sample
      for (iParam=0; iParam<nParam; iParam++) munext[iParam] = muprev[iParam];
    }

    // Write out samples stored in ParamSample
    if ((nBruteForce <= 1) && ((((iSample+1) % WriteEvery) == 0) || (iSample == nSample-1))){
      xf_printf("# samples = %d: writing samples to %s\n", iSample+1, SampleFile); 
      ierr = xf_Error(xf_WriteSampleFile(iSample+1, nParam, ParamSample, ParamName, 
					 ParamRange, PropType, PropScale, PropDelta, SampleFile));
      if (ierr != xf_OK) return ierr;
    }
    
  } // iSample

  // print out info on chain
  xf_printf("MCMC chain finished.\n");

  /* How many times outputs were calculated */
  xf_printf("Outputs were calculated %d times. \n", (nSample - countOutput));

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
  xf_Release2( (void **) ParamName);
  xf_Release(  (void  *) ParamRange);
  xf_Release(  (void  *) ParamValue);
  xf_Release2( (void **) ParamSample);
  xf_Release2( (void **) OutName);
  xf_Release(  (void  *) OutTimeIndex);
  xf_Release(  (void  *) OutValue);
  xf_Release(  (void  *) OutStdErr);
  xf_Release(  (void  *) J);
  xf_Release(  (void  *) epsParam);
  xf_Release(  (void  *) J_mu);
  xf_Release(  (void  *) H);
  xf_Release(  (void  *) w);

  xf_printf("xf_InverseIC_MCMC finished.\n");

  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return xf_OK;

}
