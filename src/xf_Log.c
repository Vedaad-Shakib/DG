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
  FILE:  xf_Log.c

  This file contains functions for writing to a .log file (and stdout)

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Output.h"
#include "xfYu_Statistics.h"

// maxline length for log file
#define xf_LOGLINELEN 1600


/******************************************************************/
//   FUNCTION Definition: xf_WriteLogHeader
int 
xf_WriteLogHeader( const xf_All *All, const char *PreHeader)
{
  int ierr, myRank, i, nOutput;
  enum xfe_Bool WriteLog, IsUnsteady, PenalizeResidual;
  enum xfe_Verbosity Verbosity;
  char line[xf_LOGLINELEN], LogOutput[xf_MAXLINELEN], s[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN], LogFile[xf_MAXSTRLEN];
  char **OutputList;
  FILE *fid;


  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  // Check if penalization is on
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "PenalizeResidual", 
                            &PenalizeResidual);
  if (ierr != xf_OK) return ierr;  

  // if verbosity is low, return immediately
  if (Verbosity == xfe_VerbosityLow) return xf_OK;


  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  if (myRank == 0){
    if (PreHeader != NULL)
      sprintf(line, "%s\n", PreHeader);
    else
      sprintf(line, "\0");
    if (PenalizeResidual)
      sprintf(s, "%7s %10s %7s %16s %16s %16s", "%> Iter", "CFL", "UFrac", "Rnorm", "(1+P)", "mu");
    else
      sprintf(s, "%7s %10s %7s %16s", "%> Iter", "CFL", "UFrac", "Rnorm");
    strcat(line, s);

    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "LogOutput", LogOutput));
    if (ierr != xf_OK) return ierr;

    if (xf_NotNull(LogOutput)){

      OutputList = NULL;

      ierr = xf_Error(xf_ScanXStringAlloc(LogOutput, xf_MAXSTRLEN, &nOutput, &OutputList));
      if (ierr != xf_OK) return ierr;
      
      for (i=0; i<nOutput; i++){
	ierr = xf_Error(xf_IsOutputUnsteady(All->EqnSet, OutputList[i], &IsUnsteady));
	if (ierr != xf_OK) return ierr;
	if (IsUnsteady) continue; // do not log unsteady outputs
	sprintf(s, " %16s", OutputList[i]);
	strcat(line, s);
	if (strlen(line) >= (xf_LOGLINELEN-1)){
	  xf_printf("String line for logging purposes exceeded.  Increase xf_LOGLINELEN in xf_Log.c.\n");
	  return xf_Error(xf_OUT_OF_BOUNDS);
	}
      }

      xf_Release2( (void **) OutputList);
    }

    strcat(line, "\n");
    
    xf_printf(line);

    /* write to log file only if WriteLog is True */
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "WriteLog", &WriteLog));
    if (ierr != xf_OK) return ierr;
    
    /* If SavePrefix is None or NULL, do not write anything */
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
    if (ierr != xf_OK) return ierr;
    
    if ((WriteLog) && (xf_NotNull(SavePrefix))){
      sprintf(LogFile, "%s.log", SavePrefix);
      // open .log file for append
      if ((fid=fopen(LogFile,"a"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);
      fprintf(fid, "%s", line);
      fclose(fid);
    }

  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteLogEntry
int 
xf_WriteLogEntry( xf_All *All, xf_SolverData *SolverData, xf_Vector *U)
{
  int ierr, myRank, i, j, nOutput, nProc;
  enum xfe_Bool WriteLog, ParallelFlag, IsUnsteady;
  enum xfe_Verbosity Verbosity;
  char line[xf_LOGLINELEN], LogOutput[xf_MAXLINELEN], s[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN], LogFile[xf_MAXSTRLEN];
  char **OutputList;
  FILE *fid;
  real Value, *ValueSet;


  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // if verbosity is low, return immediately
  if (Verbosity == xfe_VerbosityLow) return xf_OK;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  ParallelFlag = (nProc > 1);

  if (myRank == 0){ // print header and pull off LogOutput
    if (SolverData->PenalizeResidual)
      sprintf(line, "%7d %10.2E %7.3f %16.8E %16.8E %16.8E", SolverData->iIter, SolverData->CFL,
              SolverData->UpdateFrac, SolverData->ResNorm, SolverData->ResPenalty, SolverData->mu);
    else
      //sprintf(line, "%7d %10.2E %7.3f %16.8E", SolverData->iIter, SolverData->CFL,
      //        SolverData->UpdateFrac, SolverData->ResNorm);
      if(SolverData->StabRequired)
         i = 1;
      else
         i = 0;
       sprintf(line, "%7d %16.8E %16.8E %16.8E %d", SolverData->iIter, SolverData->currenttime,
              SolverData->currenttimestep, SolverData->ResNorm, i);

    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "LogOutput", LogOutput));
    if (ierr != xf_OK) return ierr;
  }

  if (ParallelFlag){ // broadcast LogOutput if in parallel
    ierr = xf_Error(xf_MPI_Bcast((void *) LogOutput, xf_MAXSTRLEN*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
  }

  // each proc calls CalculateOutput if LogOutput != NULL
  // old stuff not used any more
  /*  
  if (xf_NotNull(LogOutput)){
    
    OutputList = NULL;
    
    ierr = xf_Error(xf_ScanXStringAlloc(LogOutput, xf_MAXSTRLEN, &nOutput, &OutputList));
    if (ierr != xf_OK) return ierr;
    
    for (i=0; i<nOutput; i++){
      // calculate output i (if steady)
      ierr = xf_Error(xf_IsOutputUnsteady(All->EqnSet, OutputList[i], &IsUnsteady));
      if (ierr != xf_OK) return ierr;
      if (IsUnsteady) continue; // do not want to overwrite an accumulating output
      ierr = xf_Error(xf_CalculateOutput(All, OutputList[i], U, &Value, NULL, xfe_Set));
      if (ierr != xf_OK) return ierr;
      if (myRank == 0){
	sprintf(s, " %16.8E", Value);
	strcat(line, s);
      }
    } // i
    
    xf_Release2( (void **) OutputList);
  }
  */


  
  if (myRank == 0){
    strcat(line, "\n");
      
    xf_printf(line);

    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "WriteLog", &WriteLog));
    if (ierr != xf_OK) return ierr;

    /* If SavePrefix is None or NULL, do not write anything */
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
    if (ierr != xf_OK) return ierr;
    
    if ((WriteLog) && (xf_NotNull(SavePrefix))){
      sprintf(LogFile, "%s.log", SavePrefix);
      // open .log file for appen
      if ((fid=fopen(LogFile,"a"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);
      fprintf(fid, "%s", line);
      fclose(fid);
    }
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteLogLine
int 
xf_WriteLogLine( xf_All *All, const char *line)
{
  int ierr, myRank, nProc;
  enum xfe_Bool WriteLog;
  enum xfe_Verbosity Verbosity;
  char SavePrefix[xf_MAXSTRLEN], LogFile[xf_MAXSTRLEN];
  char **OutputList;
  FILE *fid;

  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // if verbosity is low, return immediately
  if (Verbosity == xfe_VerbosityLow) return xf_OK;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // root prints out line
  if (myRank == 0){

    xf_printf(line);

    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "WriteLog", &WriteLog));
    if (ierr != xf_OK) return ierr;

    /* If SavePrefix is None or NULL, do not write anything */
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
    if (ierr != xf_OK) return ierr;
    
    if ((WriteLog) && (xf_NotNull(SavePrefix))){
      sprintf(LogFile, "%s.log", SavePrefix);
      // open .log file for append
      if ((fid=fopen(LogFile,"a"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);
      fprintf(fid, "%s", line);
      fclose(fid);
    }
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_LogUnsteadyOutputs
int 
xf_LogUnsteadyOutputs( xf_All *All)
{
  int ierr, myRank, i, nOutput, nOutputUnsteady=0, nProc;
  enum xfe_Bool WriteLog;
  enum xfe_Verbosity Verbosity;
  char line[xf_LOGLINELEN], LogOutput[xf_MAXLINELEN], s[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN], LogFile[xf_MAXSTRLEN];
  char **OutputList;
  xf_Output *Output;
  FILE *fid;

  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // if verbosity is low, return immediately
  if (Verbosity == xfe_VerbosityLow) return xf_OK;
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // pull off LogOutput
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "LogOutput", LogOutput));
  if (ierr != xf_OK) return ierr;

  // print out unsteady outputs
  if (xf_NotNull(LogOutput)){
    
    OutputList = NULL;
    
    ierr = xf_Error(xf_ScanXStringAlloc(LogOutput, xf_MAXSTRLEN, 
					&nOutput, &OutputList));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0) // header on root
      sprintf(line, "%% Unsteady Outputs:\n");

    for (i=0; i<nOutput; i++){
      // find output structure
      ierr = xf_Error(xf_FindOutput(All->EqnSet, OutputList[i], &Output));
      if (ierr != xf_OK) return ierr;

      // print out output i (if unsteady)
      if ((Output->Type != xfe_SumOutput) && (Output->TimeNorm == xfe_TimeNormNone)) continue;
      nOutputUnsteady++;
      if (myRank == 0){
	sprintf(s, "%s = %20.12E\n", Output->Name, Output->Value);
	strcat(line, s);
      }
    } // i
    
    xf_Release2( (void **) OutputList);
  }
    
  if (nOutputUnsteady == 0) return xf_OK;

  // root does the printing
  if (myRank == 0){

    xf_printf(line);

    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "WriteLog", &WriteLog));
    if (ierr != xf_OK) return ierr;

    /* If SavePrefix is None or NULL, do not write anything */
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
    if (ierr != xf_OK) return ierr;
    
    if ((WriteLog) && (xf_NotNull(SavePrefix))){
      sprintf(LogFile, "%s.log", SavePrefix);
      // open .log file for append
      if ((fid=fopen(LogFile,"a"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);
      fprintf(fid, "%s", line);
      fclose(fid);
    }
  }

  return xf_OK;
}




