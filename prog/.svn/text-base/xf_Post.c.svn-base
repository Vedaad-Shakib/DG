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
  FILE:  xf_Post.c

  This program performs post-processing (output calculation).

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_Param.h"
#include "xf_Output.h"
#include "xf_Residual.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_EqnSetHook.h"
#include "xf_ErrEst.h"
#include "xf_Solver.h"
#include "xf_DataMath.h"
#include "xf_MeshTools.h"
#include "xf_MPI.h"
#include "xf_Arg.h"

/******************************************************************/
//   FUNCTION Definition: xf_PrintMeshStats
static int
xf_PrintMeshStats( xf_All *All)
{
  /* Prints out information on mesh */
  int ierr;
  int egrp, elem;
  real Volume;
  xf_Vector *EG = NULL;

  xf_printf("** Mesh statistics **\n");

  // obtain element geometry info
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  if (ierr != xf_OK) return ierr;

  Volume = 0.;
  for (egrp=0; egrp<All->Mesh->nElemGroup; egrp++)
    for (elem=0; elem<All->Mesh->ElemGroup[egrp].nElem; elem++)
      Volume += EG->GenArray[egrp].rValue[elem][xfe_EGVolume];

  xf_printf("Domain volume = %.15E\n", Volume);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int len, ierr, i, ibatch, nbatch;
  int nOutput, myRank, nProc, nPsi;
  enum xfe_Bool BatchMode, Found;
  char *ArgIn[] = {"xfa", "NULL", ".xfa file name to read (contains mesh)",
		   "data", "NULL", "alternate .data file to use (optional)",
		   "eqn", "NULL", "alternate eqnset file to use (optional)",
		   "output", "NULL", "name of output to calculate (all = everything)",
		   "f", "NULL", "file name to which output is written (optional)",
		   "batch", "0", "# data files to process in batch mode (name_#.data)",
		   "state", "NULL", "alternate name of state to use (optional)",
		   "l2err", "NULL", "exact .xfa or .data for L2 error calc (optional)",
		   "dump", "NULL", "dump file name for distributions (optional)",
		   "errest", "False", "if True, adjoint-based output error est. is computed",
 		   "dostats", "False", "if True, mesh stats will be printed",
		   "\0"};
  char xfaFile[xf_MAXSTRLEN];  
  char EqnSetFile[xf_MAXSTRLEN];  
  char dataRoot[xf_MAXSTRLEN];
  char dataFile[xf_MAXSTRLEN];
  char outFile[xf_MAXSTRLEN];
  char OutputName[xf_MAXSTRLEN];
  char StateName[xf_MAXSTRLEN];
  char exactFile[xf_MAXSTRLEN];
  char DumpFile[xf_MAXSTRLEN];
  char *AltName, *ext, *stemp;
  enum xfe_Bool ErrEst = xfe_False;
  enum xfe_Bool dostats = xfe_False;
  FILE *fout;
  real Value, OutputError, ErrorSum;
  xf_DataSet *ExactDataSet = NULL;
  xf_Data *D;
  xf_EqnSet *EqnSet;
  xf_Vector *U, **Psi = NULL, *ExactU = NULL;
  xf_Vector *AdaptIndicator;
  xf_KeyValue KeyValue;
  xf_Output *Output;
  xf_All *All, *ExactAll = NULL;
  
  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;


  xf_printf("\n");
  xf_printf("=== xf_Post: Post-Processing  ===\n");
  xf_printf("\n");
          

  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  // Get xfaFile
  ierr = xf_GetKeyValue(KeyValue, "xfa", xfaFile);
  if (ierr != xf_OK) return ierr;

  // Get dataRoot 
  ierr = xf_GetKeyValue(KeyValue, "data", dataRoot);
  if (ierr != xf_OK) return ierr;

  // Get EqnSetFile 
  ierr = xf_GetKeyValue(KeyValue, "eqn", EqnSetFile);
  if (ierr != xf_OK) return ierr;

  // Get OutputName
  ierr = xf_GetKeyValue(KeyValue, "output", OutputName);
  if (ierr != xf_OK) return ierr;

  // Get outFile
  ierr = xf_GetKeyValue(KeyValue, "f", outFile);
  if (ierr != xf_OK) return ierr;

  // Get nbatch
  ierr = xf_GetKeyValueInt(KeyValue, "batch", &nbatch);
  if (ierr != xf_OK) return ierr;

  // Get StateName
  ierr = xf_GetKeyValue(KeyValue, "state", StateName);
  if (ierr != xf_OK) return ierr;

  // Get exactFile
  ierr = xf_GetKeyValue(KeyValue, "l2err", exactFile);
  if (ierr != xf_OK) return ierr;

  // Get DumpFile
  ierr = xf_GetKeyValue(KeyValue, "dump", DumpFile);
  if (ierr != xf_OK) return ierr;

  // Get ErrEst
  ierr = xf_GetKeyValueBool(KeyValue, "errest", &ErrEst);
  if (ierr != xf_OK) return ierr;

  // Get DoStats
  ierr = xf_GetKeyValueBool(KeyValue, "dostats", &dostats);
  if (ierr != xf_OK) return ierr;


  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;

  // set BatchMode flag
  BatchMode = (nbatch > 0);
  if (!BatchMode) nbatch = 1;

  
  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Read .xfa or .gri file */
  len = strlen(xfaFile);
  if (len < 4) return xf_Error(xf_INPUT_ERROR);
  ext = xfaFile + len - 4; // pointer to extension
  if (strncmp(ext, ".xfa", 4) == 0){
    ierr = xf_Error(xf_ReadAllBinary(xfaFile, All));
    if (ierr!=xf_OK) return ierr;
  }
  else if (strncmp(ext, ".gri", 4) == 0){
    ierr = xf_Error(xf_ReadGriFile(xfaFile, NULL, All->Mesh));
    if (ierr!=xf_OK) return ierr;
    if (!xf_NotNull(dataRoot)) nbatch = 0;
  }

  /* Read eqnset file if specified*/
  if (xf_NotNull(EqnSetFile)){
    
    ierr = xf_Error(xf_CreateEqnSet(&EqnSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadEqnSetFile(EqnSetFile, NULL, EqnSet));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_DestroyEqnSet(All->EqnSet, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    All->EqnSet = EqnSet;
    All->EqnSet->Dim = All->Mesh->Dim;
  }


  // Dynamically load eqnset library
  if (All->EqnSet->EqnSetLibrary != NULL){
    ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
    if (ierr != xf_OK) return ierr;
    
    // Register the equation set
    ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
    if (ierr != xf_OK) return ierr;
  }  

  // print out stats if requested
  if (dostats){
    ierr = xf_Error(xf_PrintMeshStats(All));
    if (ierr != xf_OK) return ierr;
  }

  // open output file for writing
  if (myRank == 0){
    if (xf_NotNull(outFile)){
      if ((fout = fopen(outFile, "w")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    }
    else
      fout = stdout;
  }
  
  xf_printf("BatchMode = %d\n", BatchMode);


  // load exact .xfa or .data
  if (xf_NotNull(exactFile)){
    if ((len = strlen(exactFile)) < 5) return xf_Error(xf_FILE_READ_ERROR);
    if (strncmp(exactFile + len - 4, ".xfa", 4) == 0 ){ // read ExactAll
      ierr = xf_Error(xf_CreateAll(&ExactAll, xfe_False));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReadAllBinary(exactFile, ExactAll));
      if (ierr!=xf_OK) return ierr;
      ExactDataSet = ExactAll->DataSet;
    }
    else if (strncmp(exactFile + len - 5, ".data", 5) == 0 ){ // read ExactDataSet
      ierr = xf_Error(xf_CreateDataSet(&ExactDataSet));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, exactFile, ExactDataSet));
      if (ierr!=xf_OK) return ierr;
    }
    else if (!BatchMode) return xf_Error(xf_FILE_READ_ERROR);
     
    if (ExactDataSet != NULL){
      ierr = xf_Error(xf_FindPrimalState(ExactDataSet, 0, &D, NULL));
      if (ierr != xf_OK) return ierr;
      ExactU = (xf_Vector *) D->Data;
    }
  }


  // loop over batches
  for (ibatch=0; ibatch<nbatch; ibatch++){

    // load alternate .data file (required in batch-mode)
    if (xf_NotNull(dataRoot)){
      if (BatchMode)
	sprintf(dataFile, "%s_%d.data", dataRoot, ibatch);
      else
	strcpy(dataFile, dataRoot);

      xf_printf("Reading DataSet file: %s\n", dataFile);
      
      ierr = xf_Error(xf_DestroyDataSet(All->DataSet));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_CreateDataSet(&All->DataSet));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, dataFile, All->DataSet));
      if (ierr!=xf_OK) return ierr;
    }

    // find state
    AltName = ((xf_NotNull(StateName)) ? StateName : NULL);
    ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D, AltName));
    if (ierr != xf_OK) return ierr;
    U = (xf_Vector *) D->Data;

    // calculate output
    if (xf_NotNull(OutputName)){
      if (strncmp(OutputName, "all", 3) == 0){
	// using all outputs
	nOutput = All->EqnSet->Outputs->nOutput;
	Output = All->EqnSet->Outputs->Output;
      }
      else{
	// using only one output -- find it
	nOutput = 1;
	ierr = xf_Error(xf_FindOutput(All->EqnSet, OutputName, &Output));
	if (ierr != xf_OK) return ierr;
      }
      
      // loop over outputs and calculate each one
      for (i=0; i<nOutput; i++){
	
	// Set DumpFile if using one (or an alternate one)
	stemp = Output[i].DumpFile;
	if (xf_NotNull(DumpFile)) Output[i].DumpFile = DumpFile;
	
	ierr = xf_Error(xf_CalculateOutput(All, Output[i].Name, U, &Value, NULL, xfe_Set));
	if (ierr != xf_OK) return ierr;
	
	// reset DumpFile
	if (xf_NotNull(DumpFile)) Output[i].DumpFile = stemp;

	// root prints to file
	if (myRank == 0) fprintf(fout, "%.15E\n", Value);

	// call adjoint-based error estimation if desired
	if (ErrEst){
	  
	  // determine if we need to calculate an adjoint
	  ierr = xf_Error(xf_FindAdjointVectors(All, U, Output[i].Name, xfe_False,
						xfe_False, &nPsi, &Psi, &Found));
	  if (ierr != xf_OK) return ierr;
    
	  if (!Found){ // adjoint did not exist, need to solve for it
	    xf_printf("\nDid not find adjoint for %s; calling adjoint solver.\n", Output[i].Name);
	    ierr = xf_Error(xf_SolveAdjoints(All, 0.0, 1.0, xfe_False, U, nPsi, NULL, Psi, NULL, xfe_True,xfe_True));
	    if (ierr != xf_OK) return ierr;
	  }
	  xf_Release( (void *) Psi);
	  
	  // first get AdaptIndicator
	  ierr = xf_Error(xf_FindVector(All, "AdaptIndicator", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
					NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  
					xfe_False, NULL, &AdaptIndicator, NULL));
	  if (ierr != xf_OK) return ierr;

	  // next call error estimation in the output
	  xf_printf("Output.Name = %s\n", Output[i].Name);
	  ierr = xf_Error(xf_ErrEstOutput(All, Output[i].Name, "None", xfe_False,
					  xfe_False,  AdaptIndicator, xfe_False, &OutputError));
	  if (ierr != xf_OK) return ierr;
	  
	  ierr = xf_Error(xf_VectorNorm(AdaptIndicator, 1, &ErrorSum));
	  if (ierr != xf_OK) return ierr;

	  xf_printf(" Sum of error indicators = \n%.10E\n", ErrorSum);

	  xf_printf(" Output error estimate = \n%.10E\n", OutputError);
	
	  // destroy AdaptIndicator
	  xf_Error(xf_DestroyVector(AdaptIndicator, xfe_True));
          if (ierr != xf_OK) return ierr;


	}

      }
    }

    // pull off exact solution for batch mode
    if ((xf_NotNull(exactFile)) && (BatchMode)){
      ierr = xf_Error(xf_DestroyDataSet(ExactDataSet));
      if (ierr != xf_OK) return ierr;

      sprintf(dataFile, "%s_%d.data", exactFile, ibatch);

      ierr = xf_Error(xf_CreateDataSet(&ExactDataSet));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, dataFile, ExactDataSet));
      if (ierr!=xf_OK) return ierr;

      ierr = xf_Error(xf_FindPrimalState(ExactDataSet, 0, &D, NULL));
      if (ierr != xf_OK) return ierr;
      ExactU = (xf_Vector *) D->Data;      
    }

    // calculate error norm
    if (ExactU != NULL){

      ierr = xf_Error(xf_OutputExactErrorNorm(All, U, ExactU, &Value));
      if (ierr != xf_OK) return ierr;
      // print Value to output file
      if (myRank == 0) fprintf(fout, "%.15E\n", Value);
    }

  } // ibatch
  
  // close output file if wrote to one
  if (xf_NotNull(outFile)) fclose(fout);


  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;


  // destroy exact structure
  if (ExactAll != NULL){
    ierr = xf_Error(xf_DestroyAll(ExactAll));
    if (ierr!=xf_OK) return ierr;
    ExactDataSet = NULL;
  }
  if (ExactDataSet != NULL){
    ierr = xf_Error(xf_DestroyDataSet(ExactDataSet));
    if (ierr!=xf_OK) return ierr;
  }

  xf_printf("xf_Post finished.\n");

  /* MPI finalize */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
