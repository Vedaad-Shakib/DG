/*------------------------------------------------------------------*/
/* XFLOW: A discontinuous Galerkin finite element software library. */
/*                                                                  */
/*                    Copyright  2007-2008                          */
/*           Krzysztof J. Fidkowski, kfid@alum.mit.edu              */
/*                                                                  */
/*                    Copyright  2008-2011                          */
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
 FILE:  xf_Exec.c
 
 This file contains the main function that executes xflow according
 to a job file.
 
 */

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_Memory.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_EqnSet.h"
#include "xf_ParamDefault.h"
#include "xf_EqnSetHook.h"
#include "xf_SolverStruct.h"
#include "xfYu_Solver.h"
#include "xf_Solver.h"
#include "xfYu_Model.h"
#include "xfYu_Limiter.h"
#include "xf_Log.h"
#include "xf_Output.h"
#include "xf_Math.h"
#include "xf_Adapt.h"
#include "xf_Partition.h"
#include "xfYu_EntropyBounding.h"
#include "xfYu_AdaptSolver.h"
#include <time.h>  // for timing execution

/******************************************************************/
//   FUNCTION Definition: xf_SanityCheck
static int
xf_SanityCheck(xf_All *All)
{
  int ierr, i, ibfgrp;
  int AdaptIter, nPsi;
  int UEIWriteInterval;
  xf_KeyValue KeyValue;
  enum xfe_Bool found, DynamicSpatialRef;
  enum xfe_AdaptOnType AdaptOn;
  char AdjointOutputs[xf_MAXSTRLEN], AdaptOutput[xf_MAXSTRLEN], s[xf_MAXSTRLEN];
  char AdaptVariableSet[xf_MAXSTRLEN];
  char **OutputNames = NULL;
  
  KeyValue = All->Param->KeyValue;
  
  /* Adaptation sanity check*/
  ierr = xf_Error(xf_GetKeyValueInt(KeyValue, "AdaptIter", &AdaptIter));
  if (ierr != xf_OK) return ierr;
  if (AdaptIter >= 1){
    ierr = xf_Error(xf_GetKeyValueEnum(KeyValue, "AdaptOn", xfe_AdaptOnName, 
                                       (int ) xfe_AdaptOnLast, (int *) &AdaptOn));
    if (ierr != xf_OK) return ierr;
    
    if (AdaptOn == xfe_AdaptOnOutput){
      /* Determine which output we are adapting on */
      ierr = xf_Error(xf_GetKeyValue(KeyValue, "AdaptOutput", AdaptOutput));
      if (ierr != xf_OK) return ierr;
      
      // pull off adjoint outputs
      ierr = xf_Error(xf_GetKeyValue(KeyValue, "AdaptVariableSet", AdaptVariableSet));
      if (ierr != xf_OK) return ierr;
      
      if (!xf_NotNull(AdaptVariableSet)){
        
        // pull off adjoint outputs
        ierr = xf_Error(xf_GetKeyValue(KeyValue, "AdjointOutputs", AdjointOutputs));
        if (ierr != xf_OK) return ierr;
        
        if (xf_NotNull(AdjointOutputs)){
          ierr = xf_Error(xf_ScanXStringAlloc(AdjointOutputs, xf_MAXSTRLEN, &nPsi, &OutputNames));
          if (ierr != xf_OK) return ierr;
          
          /* Make sure AdaptOutput exists in AdjointOutputs */
          for (i=0, found = xfe_False; i<nPsi; i++)
            if (strcmp(OutputNames[i], AdaptOutput) == 0) found = xfe_True;
        }
        else found = xfe_False;
        if (!found){
          xf_printf("Warning, output adaptation requested but AdaptOutput = %s\n", AdaptOutput);
          xf_printf("  is not in AdjointOutputs = %s.\n", AdjointOutputs);
          xf_printf("  Augmenting AdjointOutputs and continuing.\n");
          sprintf(s, " %s", AdaptOutput);
          strcat(AdjointOutputs, s);
          ierr = xf_Error(xf_SetKeyValue(KeyValue, "AdjointOutputs", AdjointOutputs));
          if (ierr != xf_OK) return ierr;
        }
      }
      
    }
  }
  
  xf_Release2( (void **) OutputNames);
  
  /* Unsteady dynamic refinement requires writing out of unsteady
   error indicators at every time step. */
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "DynamicSpatialRef", &DynamicSpatialRef));
  if (ierr != xf_OK) return ierr;
  
  if (DynamicSpatialRef){
    ierr = xf_Error(xf_GetKeyValueInt(KeyValue, "UErrEstIndicatorWriteInterval", &UEIWriteInterval));
    if (ierr != xf_OK) return ierr;
    if (UEIWriteInterval != 1){
      xf_printf("Setting UErrEstIndicatorWriteInterval=1 because DynamicSpatialRef == True\n");
      ierr = xf_Error(xf_SetKeyValueInt(KeyValue, "UErrEstIndicatorWriteInterval", 1));
      if (ierr != xf_OK) return ierr;
    }
  }

  /* Warn user if any boundary face groups will be skipped because of
     zero measure (indicated in title) */
  for (ibfgrp=0; ibfgrp<All->Mesh->nBFaceGroup; ibfgrp++)
    if (strncmp(All->Mesh->BFaceGroup[ibfgrp].Title, "ZeroMeasure", 11) == 0)
      xf_printf("NOTICE: will skip boundary face group %s during residual calc.\n",
		All->Mesh->BFaceGroup[ibfgrp].Title);
  
  return xf_OK;
  
}


/******************************************************************/
//   FUNCTION Definition: xf_main
int
xf_main(int argc, char *argv[])
{
  int len, i, ierr, iRestart;
  int myRank, nProc, iOrder, SequenceOrderAdd;
  int iEqnSet, nEqnSet, nPsi, iAdapt, AdaptIter;
  int WriteInterval;
  
  enum xfe_Bool RestartFlag, PingResidual, ParallelFlag;
  enum xfe_Bool SnapshotRestart, SearchFlag, DoneAdapt, DumpSystem; 
  enum xfe_Bool ProcViz, ScaleStateUsingIC, LinearFlag, Unsteady;
  enum xfe_TimeSchemeType TimeScheme;
  
  char jobFile[xf_MAXSTRLEN];
  
  char OutputFile[xf_MAXSTRLEN];
  char EqnSetFileRoot[xf_MAXSTRLEN];
  char EqnSetFile[xf_MAXSTRLEN];
  char SavePrefix[xf_MAXSTRLEN] = "None";
  char SavePrefixNext[xf_MAXSTRLEN] = "None";
  char SavePrefix0[xf_MAXSTRLEN] = "None";
  char SaveRoot[xf_MAXSTRLEN];
  char EqnSetLibName[xf_MAXSTRLEN] = "NotSet";
  char PingOutput[xf_MAXSTRLEN];
  char AdjointOutputs[xf_MAXSTRLEN];
  char LogOutput[xf_MAXLINELEN];
  char InputTimeHist[xf_MAXLINELEN] = "None";
  
  char value[xf_MAXSTRLEN];
  char PreHeader[xf_MAXSTRLEN];
  
  real CFLStart;
  
  xf_TimeHistData *TimeHistData = NULL;
  xf_DataSet *DataSet_Glob, *DataSet;
  xf_Data *StateData, *GammaDat;
  xf_Vector *State, **Psi, *GammaVec;
  xf_EqnSet *EqnSet;
  xf_ICs *ICsOrig = NULL;
  xf_All *All, *All_Glob; 
  Yu_Model Model;
  Yu_Limiter **Limiter;
  FILE *fid;
  
  clock_t clock_Start, clock_End;
  
  
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ParallelFlag = (nProc > 1);
  
  
  /* Print Info and check arguments. */
  if (myRank == 0){ 
    
    xf_printf("\n");
    xf_printf("+-----------------------------------+\n");
    xf_printf("| xflow: A DG Finite Element Solver |\n");  
    xf_printf("+-----------------------------------+\n");
    if (ParallelFlag) xf_printf("Parallel run on %d processors\n", nProc);
    xf_printf("\n");
    
    /* Check number of arguments */
    if( argc != 2 ){
      xf_printf("Usage:\n");
      xf_printf("xflow <jobfile>\n");
      xf_printf("\n");
      xf_printf("Where <jobfile> is the name of the job file.\n");
      xf_printf("\n");
      return xf_Error(xf_INPUT_ERROR);
    }
    
    // root sets jobFile
    strcpy(jobFile, argv[1]);
  }
  
  /*-----------------------------------------*/
  /* Read in model.yu for his problems. OK?  */
  /*-----------------------------------------*/
  ierr = xf_Error(PullinModel(&Model));
  if (ierr != xf_OK) return ierr;
  /*------------------------------------------*/
  /* Read in .job file into All.  Do not read */
  /* in EqnSet from file yet. Parallelization */
  /* takes place in this function.            */
  /*------------------------------------------*/

  ierr = xf_Error(xf_ReadAllFromJobFile(jobFile, xfe_False, &All));
  if (ierr != xf_OK) return ierr;
  All->Model = &Model; 

  /* Sanity check on parameters (to avoid errors further along) */
  ierr = xf_Error(xf_SanityCheck(All));
  if (ierr != xf_OK) return ierr;

  /*---------------------*/
  /* Get params from All */
  /*---------------------*/
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "Restart", &RestartFlag));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "SnapshotRestart", &SnapshotRestart));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFLStart));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix0));
  if (ierr != xf_OK) return ierr;
  strcpy(SavePrefix, SavePrefix0);
  
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "EqnSetFile", EqnSetFileRoot));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "iRestart", &iRestart);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "DumpSystem", &DumpSystem));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ScaleStateUsingIC", &ScaleStateUsingIC));
  if (ierr != xf_OK) return ierr;

  if ((ScaleStateUsingIC) && (!RestartFlag)){
    xf_printf("Cannot scale state using IC when Restart is False.\n");
    xf_printf("Continuing without rescaling.\n");
    ScaleStateUsingIC = xfe_False;
  }
    
    /*-------------------------------*/
    /* Begin Loop Over Equation Sets */
    /*-------------------------------*/
    
    ierr = xf_GetKeyValueInt(All->Param->KeyValue, "nEqnSet", &nEqnSet);
    if (ierr != xf_OK) return ierr;
    
    // Set EqnSetFile name
    if (nEqnSet == 1)
        strcpy(EqnSetFile, EqnSetFileRoot);
    else
        sprintf(EqnSetFile, "%s_%d.eqn", EqnSetFileRoot, iEqnSet);
    
    
    // Read EqnSetFile
    if (xf_NotNull(EqnSetFile)){
        
        ierr = xf_Error(xf_CreateEqnSet(&EqnSet));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ReadEqnSetFile(EqnSetFile, NULL, EqnSet));
        if (ierr != xf_OK) return ierr;
    }
    else 
        EqnSet = NULL;
    
    /*---------------------------------------------*/
    /* Update All->EqnSet if read EqnSet from file */
    /*---------------------------------------------*/
    
    // create an equation set if taking snapshots
    if (nEqnSet > 1){
        if (EqnSet == NULL){
            xf_printf("Error, need eqnset from file for snapshots.\n");
            return xf_Error(xf_INPUT_ERROR);
        }
        ierr = xf_Error(xf_DestroyEqnSet(All->EqnSet, xfe_True));
        if (ierr != xf_OK) return ierr;
        
        All->EqnSet = EqnSet;
        ICsOrig = NULL;
    }
    else if (EqnSet != NULL){
        
        // pull off original ICs if need to scale via ICs
        if (ScaleStateUsingIC){
            ICsOrig = All->EqnSet->ICs;
            All->EqnSet->ICs = NULL;
        }
        else ICsOrig = NULL;
        
        /* All->EqnSet is updated with values (params, ResTerms, ICs, BCs,
         Outputs) from EqnSet. */
        ierr = xf_Error(xf_MergeEqnSet(All->EqnSet, EqnSet));
        if (ierr != xf_OK) return ierr;
        
        /* No longer need the EqnSet read in from the .eqn file */
        ierr = xf_Error(xf_DestroyEqnSet(EqnSet, xfe_True));
        if (ierr != xf_OK) return ierr;
    }
   
    /*-----------------------*/
    /* Pull in Limiter Params*/
    /*-----------------------*/
    if(Model.LimiterFlag){
    ierr = xf_Error(xf_FullinLimiterStruct(All, &Model, &Limiter));
    if (ierr != xf_OK) return ierr;
    }

    /*-----------------------*/
    /* perform initialization*/
    /*-----------------------*/
    if(Model.GammaVaryFlag)
    {
       //first set gamma all to constant same as that specified in Model.yu
       ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, NULL, "HeatCapacityRatio", Model.GammaInit));
       if(ierr != xf_OK) return ierr;
       //initialize with constant gamma
       ierr = xf_Error(xfYu_FindOrCreatePrimalState(All, Model, RestartFlag, NULL, &State));
       if (ierr != xf_OK) return ierr;
       //update the gamma according to the given initial condtions
      // ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
      // if(ierr == xf_NOT_FOUND)
      // {
      //    xf_printf("Cannot find heat capacity ratio...\n");
      //    return ierr;
      // }
      // else
      //    GammaVec = (xf_Vector *) GammaDat->Data;
      // ierr = xf_Error(Yu_GammaVectorUpdate(All, &Model, GammaVec, State));
      // if (ierr != xf_OK) return ierr;
       //finally, re-initialize variable using the updated gamma variable
      // ierr = xf_Error(xfYu_FindOrCreatePrimalState(All, Model, SearchFlag, NULL, &State));
       ierr = xf_Error(Yu_EnergyCorrection(All, &Model, State));
       if (ierr != xf_OK) return ierr;
    }
    else
    {
     
       ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, NULL, "HeatCapacityRatio", Model.GammaInit));
       if(ierr != xf_OK) return ierr;
       ierr = xf_Error(xfYu_FindOrCreatePrimalState(All, Model, RestartFlag, NULL, &State));
       if (ierr != xf_OK) return ierr;
       //ierr = xf_Error(xf_FindOrCreatePrimalState(All, 0, NULL, &State));
       //if (ierr != xf_OK) return ierr;

    //specify another global-mesh-link vector output
    ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, State, "MaxPressure", 0.0));
    if(ierr != xf_OK) return ierr;
    /*-----------------------*/
    /*  Prepare for Stepping */
    /*-----------------------*/

    //for flexible time integration
    ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, State, "ElemMaxCharSpeed", 0.0));
    if(ierr != xf_OK) return ierr;
    ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, State, "ElemMinFaceLen", 0.0));
    if(ierr != xf_OK) return ierr;

    //one-time call to obtain length scale info for time step evaluation
    ierr = xf_Error(Yu_MinFaceLengthScale(All, &Model, State));
    if(ierr != xf_OK) return ierr;

    if(Model.AVmodel){
       ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, State, "AVmodel", 0.0));
       if(ierr != xf_OK) return ierr;
    }
    }//code logic!!

    /*-----------------------*/
    /* Solving Process       */
    /*-----------------------*/
    // create a default uniform time history
    ierr = xf_Error(xf_CreateUniformTimeHistData(All, LogOutput, &TimeHistData));
    if (ierr != xf_OK) return ierr;
   
    clock_Start = clock(); // start timer
    xf_printf("Starting timer\n");

    if(Model.Dyn_p_Adapt)
    {
       ierr = xf_Error(xfYu_ApplyUnsteadyAdapt(All, &Model, NULL, SavePrefix, xfe_False, State, TimeHistData));
       if (ierr != xf_OK) return ierr;
    }
    else
    {
       if(Model.LimiterFlag)
          ierr = xf_Error(xfYu_ApplyTimeScheme(All, &Model, Limiter, SavePrefix, xfe_False, State, TimeHistData));
       else
          ierr = xf_Error(xfYu_ApplyTimeScheme(All, &Model, NULL, SavePrefix, xfe_False, State, TimeHistData));
       //ierr = xf_Error(xf_ApplyTimeScheme(All, SavePrefix, xfe_False, &State, TimeHistData));
       if (ierr != xf_OK) return ierr;
    }
    
    clock_End = clock(); // end timer
    
    xf_printf("Forward CPU time = %.10E\n", ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC);
 
    /*----------------------*/
    /* Destroy Time history */
    /*----------------------*/
    ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
    if (ierr != xf_OK) return ierr;
    
    /*----------------------*/
    /* Writing Results      */
    /*----------------------*/
    sprintf(SaveRoot, "%s_0\0", SavePrefix);
    sprintf(OutputFile, "%s.xfa\0", SaveRoot);
    ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
    if (ierr!=xf_OK) return ierr;

    if(Model.LimiterFlag){
    ierr = xf_Error(DestroyLimiterStruct(All, Limiter));
    if (ierr!=xf_OK) return ierr;
    }

    ierr = xf_Error(DestroyModel(&Model));
    if (ierr!=xf_OK) return ierr;
   
    ierr = xf_Error(xf_DestroyAll(All));
    if (ierr!=xf_OK) return ierr;

    xf_printf("finished by YU!~\n");

    return xf_OK;
    //~~~~~~~~~~~~~~~~end of the code~~~~~~~~~~~~~~~~~~//
    
    for (iEqnSet=0; iEqnSet<nEqnSet; iEqnSet++)
    {
    /*----------------------*/
    /* Load EqnSet Library  */
    /*----------------------*/
    
    strncpy(EqnSetLibName, All->EqnSet->EqnSetLibrary, xf_MAXSTRLEN);
    ierr = xf_Error(xf_LoadEqnSetLibrary(EqnSetLibName));
    if (ierr != xf_OK) return ierr;
    
    /*----------------*/
    /* Set EqnSet Dim */
    /*----------------*/
    
    All->EqnSet->Dim = All->Mesh->Dim;
    
    
    /*-----------------------------*/
    /* Perform grid initialization */
    /*-----------------------------*/
    
    /* Increment iRestart if RestartFlag is True*/
    if ((RestartFlag) && (iEqnSet == 0)){
      ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "iRestart", &iRestart));
      if (ierr != xf_OK) return ierr;
      
      iRestart++;
      
      ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "iRestart", iRestart));
      if (ierr != xf_OK) return ierr;
    }
    else if (!RestartFlag){
      ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "iRestart", 0));
      if (ierr != xf_OK) return ierr;
    }
    
    
    /* Determine time scheme */
    ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
                                       xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
                                       (int *) &TimeScheme));
    if (ierr != xf_OK) return ierr;
    
    Unsteady = (TimeScheme != xfe_TimeSchemeSteady); // is time scheme unsteady?
    
    
    /* Determine if adaptation iterations are requested */
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "AdaptIter", &AdaptIter));
    if (ierr != xf_OK) return ierr;
    if (AdaptIter < 0) return xf_Error(xf_OUT_OF_BOUNDS);

    /* determine if mesh motion is active; set flag for quick access */
    if (All->Mesh->Motion != NULL){
      ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "MeshMotionActive", 
					 &All->Mesh->Motion->Active));
      if (ierr != xf_OK) return ierr;
    }
    
    
    /*-----------------------------*/
    /* Start Adaptation Iterations */
    /*-----------------------------*/
    
    for (iAdapt=0; iAdapt<=AdaptIter; iAdapt++){
      
      if (AdaptIter > 0) xf_printf("\nStarting adaptation iteration %d.\n", iAdapt);
      
      // Modify SavePrefix if performing unsteady adaptation
      if ((AdaptIter > 0) && (Unsteady)){
        sprintf(SavePrefix, "%s_UA%d", SavePrefix0, iAdapt);  
        sprintf(SavePrefixNext, "%s_UA%d", SavePrefix0, iAdapt+1);  
        ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "SavePrefix", SavePrefix));
        if (ierr != xf_OK) return ierr;
      }
      
      /* Set CFL */
      ierr = xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart);
      if (ierr != xf_OK) return ierr;
      
      
      /*-----------------------------------------*/
      /* (Re)Register equation based on ResTerms */
      /*-----------------------------------------*/
      
      /* StateRank, StateName, PosInState, FcnParamI, and FcnParamR set
       here.  Default eqnset parameters also set here. */
      ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
      if (ierr != xf_OK) return ierr;
      
      
      /*-----------------------------------------------*/
      /* Find or create+initialize primal state vector */
      /*-----------------------------------------------*/
      
      // if doing steady adaptation, state should already be there for iAdapt>0
      if ((iAdapt > 0) && (!Unsteady)){
        ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &StateData, NULL));
        if (ierr != xf_OK) return ierr;  
        State = (xf_Vector *) StateData->Data;
      }
      else{
        // search for existing State if any of the following is true
        SearchFlag = (RestartFlag ||                           // user wants restart
                      ((iEqnSet > 0) && (SnapshotRestart)));   // taking snapshots
        ierr = xf_Error(xf_FindOrCreatePrimalState(All, SearchFlag, ICsOrig, &State));
        if (ierr != xf_OK) return ierr;  
      }
      
      
      
      /* Ping Residual if requested */
      ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "PingResidual", &PingResidual));
      if (ierr != xf_OK) return ierr;
      if (PingResidual){
        ierr = xf_Error(xf_PingResidual(All, State, 1e-5, -1.0));
        if (ierr != xf_OK) return ierr;
      }
      
      /* Ping Output if requested */
      ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "PingOutput", PingOutput));
      if (ierr != xf_OK) return ierr;
      if (xf_NotNull(PingOutput)){
        ierr = xf_Error(xf_PingOutput(All, PingOutput, State));
        if (ierr != xf_OK) return ierr;
      }
      
      /* Create a processor ID vector, if requested*/
      ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ProcViz", &ProcViz));
      if (ierr != xf_OK) return ierr;
      if (ProcViz){
        ierr = xf_Error(xf_ProcViz(All));
        if (ierr != xf_OK) return ierr;
      }
      
      
      /* Write header */
      if (nEqnSet > 1)
        sprintf(PreHeader, "%s %d", "%% Snapshot", iEqnSet);
      else
        sprintf(PreHeader, "%%");
      ierr = xf_Error(xf_WriteLogHeader(All, PreHeader));
      if (ierr != xf_OK) return ierr;
      
      
      /*---------------------------*/
      /* Forward and Adjoint Solve */
      /*---------------------------*/
      
      ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdjointOutputs", 
                                     AdjointOutputs));
      if (ierr != xf_OK) return ierr;
      
      if (!Unsteady){
        
        /*---------------*/
        /* Steady Solver */
        /*---------------*/
        
        // Pull off SequenceOrderAdd
        ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "SequenceOrderAdd", 
                                          &SequenceOrderAdd));
        if (ierr != xf_OK) return ierr;
        
        // Begin order sequencing loop (only once through if no sequencing)
        for (iOrder = 0; iOrder <= SequenceOrderAdd; iOrder++){
          
          /* Set CFL */
          if (SequenceOrderAdd > 0){
            ierr = xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart);
            if (ierr != xf_OK) return ierr;
          }
          
          
          /* Call nonlinear solver directly for steady calculations */
          
          xf_printf("Calling steady solver.\n");

	  clock_Start = clock(); // start timer          
          ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &State));
          if (ierr != xf_OK){
            xf_printf("Error occured during Nonlinear solve. Attempting to exit nicely.\n");
            sprintf(SavePrefix, "%s_ERROR", SavePrefix0);
          }
          if (ierr != xf_OK) break;
	  clock_End = clock(); // end timer
	  xf_printf("Steady solve CPU time = %.10E\n", ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC);

          
          // Call Adjoint solver
          if (xf_NotNull(AdjointOutputs)){
            ierr = xf_Error(xf_FindAdjointVectors(All, State, AdjointOutputs, xfe_False,
                                                  xfe_True, &nPsi, &Psi, NULL));
            if (ierr != xf_OK) return ierr;
            xf_printf("\nCalling steady adjoint solver.\n");
	    clock_Start = clock(); // start timer          
           // ierr = xf_Error(xf_SolveAdjoints(All, 0.0, 1.0, xfe_False, State, nPsi, NULL, Psi, NULL,xfe_True));
            if (ierr != xf_OK) return ierr;
	    clock_End = clock(); // end timer
	    xf_printf("Steady adjoint solve CPU time = %.10E\n", 
		      ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC);
            xf_Release( (void *) Psi);
          }
          
          // save current order solution
          if (SequenceOrderAdd > 0){
            sprintf(SaveRoot, "%s_Order%d\0", SavePrefix, State->Order[0] );
            sprintf(OutputFile, "%s.xfa\0", SaveRoot);
            ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
            if (ierr!=xf_OK) return ierr;
          }
          
          // do not sequence order on adaptation iterations besides the first one
          if ((SequenceOrderAdd > 0) && (iOrder != SequenceOrderAdd) && (iAdapt==0)){
            ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, All->DataSet, State, NULL, 
                                                                   xfe_BasisLast, 1));
            if (ierr != xf_OK) return ierr;
            xf_printf("\nIncrementing order by 1: Order[0] is now %d.\n", State->Order[0] );
          }
        } // iOrder
        
      }
      else{
        
        /*-----------------*/
        /* Unsteady Solver */
        /*-----------------*/
        
        xf_printf("Calling unsteady solver.\n");
        
        // create time history data, only if iAdapt==0
        if (iAdapt == 0){
          ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "LogOutput", LogOutput));
          if (ierr != xf_OK) return ierr;
          ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "InputTimeHist", InputTimeHist));
          if (ierr != xf_OK) return ierr;
          if (xf_NotNull(InputTimeHist)){ // read time history from file
            ierr = xf_Error(xf_ReadTimeHistData(InputTimeHist, LogOutput, &TimeHistData));
            if (ierr != xf_OK) return ierr;
          }
          else{ // create a default uniform time history
            ierr = xf_Error(xf_CreateUniformTimeHistData(All, LogOutput, &TimeHistData));
            if (ierr != xf_OK) return ierr;
          }
        }
        
        // check flag settings for unsteady adjoint solve
        if (xf_NotNull(AdjointOutputs)){
          // Unsteady write interval
          ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, 
                                            "UnsteadyWriteInterval", 
                                            &WriteInterval));
          if (ierr != xf_OK) return ierr;
          
          // Is the system linear (as prescribed by the user)?
          ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, 
                                             "LinearFlag", &LinearFlag));
          if (ierr != xf_OK) return ierr;
          
          if ((!LinearFlag) && (WriteInterval != 1)){
            xf_printf("*** Setting WriteInterval=1 for unsteady adjoint ***\n");
            ierr = xf_SetKeyValueInt(All->Param->KeyValue, 
                                     "UnsteadyWriteInterval", 1);
            if (ierr != xf_OK) return ierr;
          }
        }
        
        /* Call time scheme if unsteady */
        
        clock_Start = clock(); // start timer
        ierr = xf_Error(xf_ApplyTimeScheme(All, SavePrefix, xfe_False, &State, TimeHistData));
        if (ierr != xf_OK) return ierr;
        clock_End = clock(); // end timer
        xf_printf("Forward CPU time = %.10E\n", ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC);
        
        // Call unsteady adjoint solver
        if (xf_NotNull(AdjointOutputs)){
          ierr = xf_Error(xf_FindAdjointVectors(All, State, AdjointOutputs, xfe_False,
                                                xfe_True, &nPsi, &Psi, NULL));
          if (ierr != xf_OK) return ierr;
          xf_printf("\nCalling unsteady adjoint solver.\n");
          clock_Start = clock(); // start timer
          ierr = xf_Error(xf_ApplyTimeSchemeAdjoint(All, SavePrefix, State, nPsi,
                                                    Psi, TimeHistData));
          if (ierr != xf_OK) return ierr;
          clock_End = clock(); // end timer
          xf_printf("Adjoint CPU time = %.10E\n", ((real) (clock_End-clock_Start)) / CLOCKS_PER_SEC);
          xf_Release( (void *) Psi);
        }
        
      }
      
      // dump system if requested
      if (DumpSystem){
        ierr = xf_Error(xf_DumpSystemSparse(All, State));
        if (ierr != xf_OK) return ierr;
      }
      
      /*-----------------------------------*/
      /* (Error Estimation and) Adaptation */
      /*-----------------------------------*/
      
      if (AdaptIter > 0){
        
        // break adaptation if user halt requested
        if (xf_CheckUserHalt(NULL)){
          xf_printf("Adaptation stopped by user halt.\n");
          break;
        }
        
        if (!Unsteady){
          if (iAdapt < AdaptIter){
            xf_printf("\nCalling steady adaptation.\n");
            ierr = xf_Error(xf_AdaptAll(All, iAdapt, &DoneAdapt, NULL));
            if (ierr != xf_OK) return ierr;
          }
        }
        else{
          DoneAdapt = (iAdapt == AdaptIter);
          xf_printf("\nCalling unsteady adaptation.\n");
          ierr = xf_Error(xf_AdaptAllUnsteady(All, iAdapt, SavePrefixNext, TimeHistData, &DoneAdapt));
          if (ierr != xf_OK) return ierr;
        }
        
        if (DoneAdapt) iAdapt = AdaptIter; // to force end of iterations
      }
      
      
    } // iAdapt (end of adaptation iteration loop)
    
    
    /*--------------------------------------------*/
    /* Write DataSet Snapshot if Taking Snapshots */
    /*--------------------------------------------*/
    if (nEqnSet > 1){      
      sprintf(OutputFile, "%s_%d.data\0", SavePrefix, iEqnSet);
      ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, All->DataSet, NULL, OutputFile));
      if (ierr != xf_OK) return ierr;
    }
    
    
    /*----------------------------*/
    /* Close equation set library */
    /*----------------------------*/
    xf_CloseEqnSetLibrary(&All); 
    
    /*----------------------*/
    /* Destroy Time history */
    /*----------------------*/
    ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
    if (ierr != xf_OK) return ierr;
    
    
  } // iEqnSet
  
  /*-----------*/
  /* Write All */
  /*-----------*/
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "iRestart", &iRestart));
  if (ierr != xf_OK) return ierr;
  sprintf(SaveRoot, "%s_%d\0", SavePrefix, iRestart);
  sprintf(OutputFile, "%s.xfa\0", SaveRoot);
  ierr = xf_Error(xf_WriteAllBinary(All, OutputFile));
  if (ierr!=xf_OK) return ierr;
 
  /*----------*/
  /* Clean up */
  /*----------*/
  
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;
  
  xf_printf("xflow finished.\n");
  
  
  return xf_OK;
  
}


#ifndef XF_REGRESSION

/******************************************************************/
//   FUNCTION Definition: main
int
main(int argc, char *argv[])
{
  int ierr, runerr;
  
  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;
  
  runerr = xf_Error(xf_main(argc, argv));
  
  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;
  
  return runerr;
}

#endif
