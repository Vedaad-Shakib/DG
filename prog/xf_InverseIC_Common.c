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
  FILE:  xf_InverseIC_Common.c

  This file contains routines used by more than on IC inversion program.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_Param.h"
#include "xf_Math.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Solver.h"



/******************************************************************/
//   FUNCTION Definition: xf_ReadICParamFileSerial
static int
xf_ReadICParamFileSerial(const char *ICParamFile, int *pnParam, 
			   char ***pParamName, real **pParamRange, 
			   real **pParamValue)
{
/*
PURPOSE: Serial version of ReadICParamFile (see below)
*/
  int ierr, iParam, nParam;
  char line0[xf_MAXLINELEN], line1[xf_MAXLINELEN], *line;
  FILE *fid;
  
  if ((fid = fopen(ICParamFile, "r")) == NULL)
    return xf_Error(xf_FILE_READ_ERROR);

  // Look for first non-blank non-comment line and read the number of parameters 
  // (one integer on this line)
  do{
    // Read in line
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    line = line0;
    
    // If blank line or comment then continue
    if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;
    
    // Read number of parameters
    ierr = sscanf(line, "%d", pnParam);
    if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);
    break;
  } while(feof(fid) == 0);
  
  nParam = (*pnParam);
  
  // Allocate memory dynamically for names, ranges, starting values	
  ierr = xf_Error(xf_Alloc2((void ***) pParamName, nParam, xf_MAXSTRLEN, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) pParamRange, 2*nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) pParamValue, nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Read in parameter names, ranges, and starting values
  iParam = 0;
  do{
    // Read in line
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL) break;
    line = line0;
    
    // If blank line or comment then continue
    if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;
    
    // Read in param name, min range, max range, starting value
    ierr = sscanf(line, "%s %lf %lf %lf", (*pParamName)[iParam], (*pParamRange)+2*iParam, 
		  (*pParamRange)+2*iParam+1, (*pParamValue)+iParam);
    if (ierr != 4) return xf_Error(xf_FILE_READ_ERROR);

    // parameter range needs to be (min,max) not other way around
    if ((*pParamRange)[2*iParam+1] < (*pParamRange)[2*iParam]) return xf_Error(xf_OUT_OF_BOUNDS);

    iParam++;
  } while(feof(fid) == 0);

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_ReadICParamFile
int
xf_ReadICParamFile(const char *ICParamFile, int *pnParam, 
		     char ***pParamName, real **pParamRange, 
		     real **pParamValue)
{
  int ierr, myRank;
  int nParam;

  // check parallel
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root does the reading
  if (myRank == 0)
    ierr = xf_Error(xf_ReadICParamFileSerial(ICParamFile, pnParam, pParamName,
					       pParamRange, pParamValue));
  if (xf_PError(&ierr, 0) != xf_OK) return ierr;  
  
  // broadcast nParam
  ierr = xf_Error(xf_MPI_Bcast((void *) pnParam, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  nParam = (*pnParam);

  // other procs allocate
  if (myRank > 0){
    // Allocate memory dynamically for names, ranges, starting values	
    ierr = xf_Error(xf_Alloc2((void ***) pParamName, nParam, xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) pParamRange, 2*nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) pParamValue, nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  // broadcast from root
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pParamName)[0], nParam*xf_MAXSTRLEN*sizeof(char), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pParamRange), 2*nParam*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pParamValue), nParam*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_ReadICParamFileEnsembleSerial
static int
xf_ReadICParamFileEnsembleSerial(const char *ICParamFileEnsemble, int *pnParam, 
			   char ***pParamName, real **pParamRange, real **pParamWindowRange)
{
/*
PURPOSE: Serial version of ReadICParamFileEnsemble (see below)
*/
  int ierr, iParam, nParam;
  char line0[xf_MAXLINELEN], line1[xf_MAXLINELEN], *line;
  FILE *fid;
  
  if ((fid = fopen(ICParamFileEnsemble, "r")) == NULL)
    return xf_Error(xf_FILE_READ_ERROR);

  // Look for first non-blank non-comment line and read the number of parameters 
  // (one integer on this line)
  do{
    // Read in line
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    line = line0;
    
    // If blank line or comment then continue
    if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;
    
    // Read number of parameters
    ierr = sscanf(line, "%d", pnParam);
    if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);
    break;
  } while(feof(fid) == 0);
  
  nParam = (*pnParam);
  
  // Allocate memory dynamically for names, domain ranges, and initial param ranges	
  ierr = xf_Error(xf_Alloc2((void ***) pParamName, nParam, xf_MAXSTRLEN, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) pParamRange, 2*nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) pParamWindowRange, 2*nParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Read in parameter names, domain ranges, and initial param ranges
  iParam = 0;
  do{
    // Read in line
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL) break;
    line = line0;
    
    // If blank line or comment then continue
    if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;
    
    // Read in param name, min range and max range
    ierr = sscanf(line, "%s %lf %lf %lf %lf", (*pParamName)[iParam], (*pParamRange)+2*iParam, 
		  (*pParamRange)+2*iParam+1, (*pParamWindowRange)+2*iParam, (*pParamWindowRange)+2*iParam+1);
    if (ierr != 5) return xf_Error(xf_FILE_READ_ERROR);

    // domain and initial parameter ranges need to be (min,max) not other way around
    if ((*pParamRange)[2*iParam+1] < (*pParamRange)[2*iParam]) return xf_Error(xf_OUT_OF_BOUNDS);
    if ((*pParamWindowRange)[2*iParam+1] < (*pParamWindowRange)[2*iParam]) return xf_Error(xf_OUT_OF_BOUNDS);

    iParam++;
  } while(feof(fid) == 0);

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_ReadICParamFileEnsemble
int
xf_ReadICParamFileEnsemble(const char *ICParamFileEnsemble, int *pnParam, 
		     char ***pParamName, real **pParamRange, real **pParamWindowRange)
{
  int ierr, myRank;
  int nParam;

  // check parallel
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root does the reading
  if (myRank == 0)
    ierr = xf_Error(xf_ReadICParamFileEnsembleSerial(ICParamFileEnsemble,
			pnParam, pParamName, pParamRange, pParamWindowRange));
  if (xf_PError(&ierr, 0) != xf_OK) return ierr;  
  
  // broadcast nParam
  ierr = xf_Error(xf_MPI_Bcast((void *) pnParam, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  nParam = (*pnParam);

  // other procs allocate
  if (myRank > 0){
    // Allocate memory dynamically for names, domain ranges, and initial param ranges	
    ierr = xf_Error(xf_Alloc2((void ***) pParamName, nParam, xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) pParamRange, 2*nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) pParamWindowRange, 2*nParam, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  // broadcast from root
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pParamName)[0], nParam*xf_MAXSTRLEN*sizeof(char), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pParamRange), 2*nParam*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pParamWindowRange), 2*nParam*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_ICParam2State
int
xf_ICParam2State(xf_All *All, int nParam, char **ParamName, real *ParamValue,
		 real *epsParam, xf_Vector *U, xf_Vector *U_mu)
{
  int ierr;
  int iOut;
  int iParam, nICParam;
  int iHeader, nHeader, iHeader0;
  enum xfe_Bool found;
  xf_IC *IC = NULL;
  xf_EqnSet *EqnSet;
  char **ICHeader = NULL;
  const char *sParam = NULL;
  char s[xf_MAXSTRLEN];
  real *ICParam = NULL;

  EqnSet = All->EqnSet;

  // sanity checks
  if (nParam <= 0) return xf_Error(xf_INPUT_ERROR);
  if ((EqnSet->ICs == NULL) || (EqnSet->ICs[0].nIC <= 0)) 
    return xf_Error(xf_INPUT_ERROR);
  
  IC = EqnSet->ICs[0].IC;
  
  // pull off initial condition header
  ierr = xf_Error(xf_ScanXStringAlloc(IC->Header, xf_MAXSTRLEN, &nHeader, &ICHeader));
  if (ierr != xf_OK) return ierr;

  // Pull off initial condition parameters
  ierr = xf_Error(xf_ScanXRealAlloc(IC->Data, &nICParam, &ICParam));
  if (ierr != xf_OK) return ierr;

  if (nICParam != nHeader) return xf_Error(xf_INPUT_ERROR);

  // set IC parameters in All->EqnSet
  for (iParam=0; iParam<nParam; iParam++){
    sParam = ParamName[iParam];
    for (iHeader=0; iHeader<nHeader; iHeader++){
      if (strcmp(ICHeader[iHeader], sParam) == 0){
	ICParam[iHeader] = ParamValue[iParam];
      }
    } // iHeader
  } // iParam
  

  // store data into a string: IC->Data
  xf_Release((void *) IC->Data);
  ierr = xf_Error(xf_Alloc( (void **) &IC->Data, 20*nHeader, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  sprintf(IC->Data, "%.10E ", ICParam[0]);
  for (iHeader=1; iHeader<nHeader; iHeader++){
    sprintf(s, "%.10E ", ICParam[iHeader]);
    strcat(IC->Data, s);
  } // iHeader

  /* Create vector with appropriate initial condition, U */
  ierr = xf_Error(xf_InitState(All, U));
  if (ierr != xf_OK) return ierr;
  

  // calculate sensitivities if requested
  if (U_mu != NULL){
    for (iParam=0; iParam<nParam; iParam++){
      sParam = ParamName[iParam];
      iHeader0 = -1;
      for (iHeader=0; iHeader<nHeader; iHeader++){
	if (strcmp(ICHeader[iHeader], sParam) == 0){
	  ICParam[iHeader0 = iHeader] += epsParam[iParam];
	  break;
	}
      } // iHeader
      if (iHeader0 < 0) return xf_Error(xf_NOT_FOUND);

      // Store header in IC->Data
      sprintf(IC->Data, "%.10E ", ICParam[0]);
      for (iHeader=1; iHeader<nHeader; iHeader++){
	sprintf(s, "%.10E ", ICParam[iHeader]);
	strcat(IC->Data, s);
      } // iHeader
      
      /* Initialize state */
      ierr = xf_Error(xf_InitState(All, U_mu+iParam));
      if (ierr != xf_OK) return ierr;

      // take finite difference
      ierr = xf_Error(xf_SetVector(U, xfe_Sub, U_mu+iParam));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorMult(U_mu+iParam, 1.0/epsParam[iParam]));
      if (ierr != xf_OK) return ierr;

      ICParam[iHeader0] -= epsParam[iParam];

    } // iParam
  }

  xf_Release2( (void **) ICHeader);
  xf_Release ( (void  *) ICParam);

  return xf_OK;

}





/******************************************************************/
//   FUNCTION Definition: xf_ReadICOutputFileSerial
static int
xf_ReadICOutputFileSerial(const char *OutputFile, int *pnOut, char ***pOutName,
			  int **pOutTimeIndex, real **pOutValue, real **pOutStdErr)
{
/*
  PURPOSE: Serial version of ReadICOutputFile (see below)
*/

  int ierr, nOut, iOut;
  char line0[xf_MAXLINELEN], line1[xf_MAXLINELEN], *line;
  real val;
  FILE *fid;
  
  if ((fid = fopen(OutputFile, "r")) == NULL)
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Look for first non-blank non-comment line and read the number of outputs
  // (one integer on this line)
  do{
    // Read in line
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);

    line = line0;
    // If blank line or comment then continue
    if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;
    
    // Read number of outputs
    ierr = sscanf(line, "%d", pnOut);
    if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);
    break;
  } while(feof(fid) == 0);
  
  nOut = (*pnOut);
  
  // Allocate memory dynamically for output name, time index, value, and standard error
  ierr = xf_Error(xf_Alloc2((void ***) pOutName, nOut, xf_MAXSTRLEN, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) pOutTimeIndex, nOut, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) pOutValue, nOut, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  if (pOutStdErr != NULL){
    ierr = xf_Error(xf_Alloc((void **) pOutStdErr, nOut, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  // Read in parameter names, ranges, and starting values
  iOut = 0;  
  do{
    // Read in line
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL) break;
    line = line0;
    
    // If blank line or comment then continue
    if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;
    
    // Read in output name, time index, output value, and standard error
    ierr = sscanf(line, "%s %d %lf %lf", (*pOutName)[iOut], (*pOutTimeIndex)+iOut, 
		  (*pOutValue)+iOut, &val);
    if (pOutStdErr != NULL) (*pOutStdErr)[iOut] = val;
    if (ierr != 4) return xf_Error(xf_FILE_READ_ERROR);
    iOut++;
  } while(feof(fid) == 0);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadICOutputFile
int
xf_ReadICOutputFile(const char *OutputFile, int *pnOut, char ***pOutName,
		    int **pOutTimeIndex, real **pOutValue, real **pOutStdErr)
{
  int ierr, myRank;
  int nOut;

  // check parallel
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root does the reading
  if (myRank == 0)
    ierr = xf_Error(xf_ReadICOutputFileSerial(OutputFile, pnOut, pOutName, 
					      pOutTimeIndex, pOutValue, pOutStdErr));
  if (xf_PError(&ierr, 0) != xf_OK) return ierr;  
  
  // broadcast nOut
  ierr = xf_Error(xf_MPI_Bcast((void *) pnOut, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  nOut = (*pnOut);

  // other procs allocate
  if (myRank > 0){    
    ierr = xf_Error(xf_Alloc2((void ***) pOutName, nOut, xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) pOutTimeIndex, nOut, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) pOutValue, nOut, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    if (pOutStdErr != NULL){
      ierr = xf_Error(xf_Alloc((void **) pOutStdErr, nOut, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
  }
  // broadcast from root
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pOutName)[0], nOut*xf_MAXSTRLEN*sizeof(char), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pOutTimeIndex), nOut*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Bcast((void *) (*pOutValue), nOut*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;  
  if (pOutStdErr != NULL){
    ierr = xf_Error(xf_MPI_Bcast((void *) (*pOutStdErr), nOut*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_CreateSensitivityVector
static int
xf_CreateSensitivityVector(xf_All *All, const char *OutputName, 
			   int TimeIndex, const char *FileName)
{
/*
PURPOSE:

  Calculates an initial condition sensitivity vector for an output
  (OutputName at TimeIndex) and writes it to disk.

INPUTS:

  All : All structure
  OutputName : name of output for which to create sensitivity
  TimeIndex : time index at which output is taken
  FileName : name of file to which vector is written
  
OUTPUTS: 
  
  None, a vector data file is written to disk

RETURN: Error code

*/
  int ierr;
  int nPsi;
  int nTimeStep;
  real Time, EndTime, NewEndTime;
  char SavePrefix[xf_MAXSTRLEN];
  xf_Vector *U = NULL;
  xf_Vector **Psi = NULL;
  xf_TimeHistData *TimeHistData = NULL;

  /* set new temporal discretization */

  // number of time steps from .job file
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nTimeStep", &nTimeStep));
  if (ierr != xf_OK) return ierr;
  // current start time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
  if (ierr != xf_OK) return ierr;
  // current end time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "EndTime", &EndTime));
  if (ierr != xf_OK) return ierr;
  // new number of time steps
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nTimeStep", TimeIndex));
  if (ierr != xf_OK) return ierr;
  // set new end time
  NewEndTime = Time + (EndTime - Time) * ((real) TimeIndex) / ((real) nTimeStep);
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "EndTime", NewEndTime)); 
  if (ierr != xf_OK) return ierr;

  /* Get SavePrefix */
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;

  // create a uniform time history
  ierr = xf_Error(xf_CreateUniformTimeHistData(All, NULL, &TimeHistData));
  if (ierr != xf_OK) return ierr;
  
  // create a state for reference
  ierr = xf_Error(xf_FindOrCreatePrimalState(All, xfe_False, NULL, &U));
  if (ierr != xf_OK) return ierr;

  // locate adjoint vector
  ierr = xf_Error(xf_FindAdjointVectors(All, U, OutputName, xfe_False,
					xfe_True, &nPsi, &Psi, NULL));
  if (ierr != xf_OK) return ierr;
  if (nPsi != 1) return xf_Error(xf_CODE_LOGIC_ERROR);

  // Call unsteady adjoint solver
  ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, SavePrefix, U, nPsi,
					    Psi, TimeHistData));
  if (ierr != xf_OK) return ierr;

  // write out adjoint (sensitivity actually)
  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "Sensitivity", Psi[0], FileName));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) Psi);


  // reset old values for time discretization
  // reset number of time steps
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "nTimeStep", nTimeStep));
  if (ierr != xf_OK) return ierr;
  // reset end time
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "EndTime", EndTime)); 
  if (ierr != xf_OK) return ierr;

  
  /* Destroy Time history */
  ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
  if (ierr != xf_OK) return ierr;

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_LoadSensitivityVectors
int
xf_LoadSensitivityVectors(xf_All *All, int nOut, char **OutName, 
			  int *OutTimeIndex, xf_Vector ***pJ_U0)
{
/*
PURPOSE:

  Check for or create output sensitivity vectors
  Also, Load these vectors into memory           

INPUTS:

  All : All structure
  nOut : number of outputs
  OutName : names of outputs
  OutTimeIndex : time indices of outputs
  
OUTPUTS: 
  
  J_U0 : sensitivity vectors

RETURN: Error code

*/
  int ierr;
  int iOut;
  char FileName[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN];
  xf_Vector **J_U0;
  FILE *fid;
  xf_DataSet *DataSet;
  xf_Data *D;

  /* Get SavePrefix */
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
  if (ierr != xf_OK) return ierr;


  // allocate J_U0
  ierr = xf_Error(xf_Alloc( (void **) pJ_U0, nOut, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  J_U0 = (*pJ_U0);

  for (iOut=0; iOut<nOut; iOut++){
    sprintf(FileName, "Sens_%s_%s_%d.data", SavePrefix, OutName[iOut], OutTimeIndex[iOut]);  

    // does the file exist?  Try opening it ...
    ierr = xf_fopen(FileName, "rb", &fid);

    // if file does not exist, create it
    if (ierr == xf_NOT_FOUND){ 
      ierr = xf_Error(xf_CreateSensitivityVector(All, OutName[iOut], OutTimeIndex[iOut],
						 FileName));
      if (ierr != xf_OK) return ierr;
    }
    else if (ierr != xf_OK) return ierr;

    // load vector into memory
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, FileName, DataSet));
    if (ierr != xf_OK) return ierr;
    D = DataSet->Head; // use first vector in data set
    J_U0[iOut] = (xf_Vector *) D->Data;
    D->Data = NULL;
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;
    
  } // iOut

  return xf_OK;

}
  



/******************************************************************/
//   FUNCTION Definition: xf_CalculateLinearOutputs
int
xf_CalculateLinearOutputs(xf_All *All, int nParam, char **ParamName, real *ParamValue,
			  int nOut, xf_Vector **J_U0, real *epsParam, real *J, 
			  real *J_mu, xf_Vector **pU0)
{
  int ierr;
  int iOut;
  int iParam, nICParam;
  int iHeader, nHeader, iHeader0;
  enum xfe_Bool found;
  xf_IC *IC = NULL;
  xf_EqnSet *EqnSet;
  char **ICHeader = NULL;
  const char *sParam = NULL;
  char s[xf_MAXSTRLEN];
  real *ICParam = NULL;

  EqnSet = All->EqnSet;

  // sanity checks
  if (nParam <= 0) return xf_Error(xf_INPUT_ERROR);
  if ((EqnSet->ICs == NULL) || (EqnSet->ICs[0].nIC <= 0)) 
    return xf_Error(xf_INPUT_ERROR);
  
  IC = EqnSet->ICs[0].IC;
  
  // pull off initial condition header
  ierr = xf_Error(xf_ScanXStringAlloc(IC->Header, xf_MAXSTRLEN, &nHeader, &ICHeader));
  if (ierr != xf_OK) return ierr;

  // Pull off initial condition parameters
  ierr = xf_Error(xf_ScanXRealAlloc(IC->Data, &nICParam, &ICParam));
  if (ierr != xf_OK) return ierr;

  if (nICParam != nHeader) return xf_Error(xf_INPUT_ERROR);

  // set IC parameters in All->EqnSet
  for (iParam=0; iParam<nParam; iParam++){
    sParam = ParamName[iParam];
    for (iHeader=0; iHeader<nHeader; iHeader++){
      if (strcmp(ICHeader[iHeader], sParam) == 0){
	ICParam[iHeader] = ParamValue[iParam];
      }
    } // iHeader
  } // iParam
  

  // store data into a string: IC->Data
  xf_Release((void *) IC->Data);
  ierr = xf_Error(xf_Alloc( (void **) &IC->Data, 20*nHeader, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  sprintf(IC->Data, "%.10E ", ICParam[0]);
  for (iHeader=1; iHeader<nHeader; iHeader++){
    sprintf(s, "%.10E ", ICParam[iHeader]);
    strcat(IC->Data, s);
  } // iHeader

  /* Create vector with appropriate initial condition, U0 */
  if ((*pU0) == NULL){
    /* First call, create U0 vector */
    ierr = xf_Error(xf_FindOrCreatePrimalState(All, xfe_False, NULL, pU0));
    if (ierr != xf_OK) return ierr;
  }
  else{
    /* Initialize state */
    ierr = xf_Error(xf_InitState(All, (*pU0)));
    if (ierr != xf_OK) return ierr;
  }
  
  // take dot product of U0 and J_U0 to create J
  for (iOut=0; iOut<nOut; iOut++){
    ierr = xf_Error(xf_VectorDot(J_U0[iOut], (*pU0), J+iOut));
    if (ierr != xf_OK) return ierr;
    //xf_printf(" iOut=%d, J = %.10E\n", iOut, J[iOut]); fflush(stdout);
  } // iOut


  // calculate sensitivities if requested
  if (J_mu != NULL){
    for (iParam=0; iParam<nParam; iParam++){
      sParam = ParamName[iParam];
      iHeader0 = -1;
      for (iHeader=0; iHeader<nHeader; iHeader++){
	if (strcmp(ICHeader[iHeader], sParam) == 0){
	  ICParam[iHeader0 = iHeader] += epsParam[iParam];
	  break;
	}
      } // iHeader
      if (iHeader0 < 0) return xf_Error(xf_NOT_FOUND);

      // Store header in IC->Data
      sprintf(IC->Data, "%.10E ", ICParam[0]);
      for (iHeader=1; iHeader<nHeader; iHeader++){
	sprintf(s, "%.10E ", ICParam[iHeader]);
	strcat(IC->Data, s);
      } // iHeader
      
      /* Initialize state */
      ierr = xf_Error(xf_InitState(All, (*pU0)));
      if (ierr != xf_OK) return ierr;

      // calculate new outputs
      for (iOut=0; iOut<nOut; iOut++){
	ierr = xf_Error(xf_VectorDot(J_U0[iOut], (*pU0), J_mu + nParam*iOut + iParam));
	if (ierr != xf_OK) return ierr;
      }

      // take finite difference
      for (iOut=0; iOut<nOut; iOut++){
	J_mu[nParam*iOut + iParam] -= J[iOut];
	J_mu[nParam*iOut + iParam] /= epsParam[iParam];
      }
	
      ICParam[iHeader0] -= epsParam[iParam];
    } // iParam
  }

  xf_Release2( (void **) ICHeader);
  xf_Release ( (void  *) ICParam);

  return xf_OK;
}


