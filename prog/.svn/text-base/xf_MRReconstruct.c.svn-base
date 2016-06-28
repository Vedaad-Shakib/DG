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
  FILE:  xf_MRReconstruct.c

  This program reconstructs FEM solution(s) from reduced coefficients

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Param.h"
#include "xf_Arg.h"
#include "xf_MRStruct.h"
#include "xf_Solver.h"




/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, N;
  int iRun, nRun, reverse;
  char *ArgIn[] = {"out", "online.out", "online output file",
		   "basis", "NULL", "basis set .data file",
		   "root", "NULL", "root for reconstructed solution .data files",
		   "reverse", "0", "if 1, V^T*M*root.data is computed.",
		   "xfa", "NULL", "Need an .xfa file if reverse==1",
		   "\0"};
  char outFile[xf_MAXSTRLEN];
  char basisFile[xf_MAXSTRLEN];
  char solRoot[xf_MAXSTRLEN];
  char solFile[xf_MAXSTRLEN];
  char xfaFile[xf_MAXSTRLEN];
  char line[xf_MAXLINELEN];
  FILE *fid;
  real *a, dp;
  xf_KeyValue KeyValue;
  xf_Vector *U;
  xf_VectorSet *BasisSet;
  xf_DataSet *DataSet, *SolDataSet;
  xf_Data *D;
  xf_All *All = NULL;
  
  xf_printf("\n");
  xf_printf("=== Model Reduction: Solution Reconstruction ===\n");
  xf_printf("\n");
    
      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  // Get online outFile 
  ierr = xf_GetKeyValue(KeyValue, "out", outFile);
  if (ierr != xf_OK) return ierr;

  // Get basisFile
  ierr = xf_GetKeyValue(KeyValue, "basis", basisFile);
  if (ierr != xf_OK) return ierr;

  // Get reverse flag
  ierr = xf_GetKeyValueInt(KeyValue, "reverse", &reverse);
  if (ierr != xf_OK) return ierr;

  // Get solRoot
  ierr = xf_GetKeyValue(KeyValue, "root", solRoot);
  if (ierr != xf_OK) return ierr;

  // Get xfaFile
  ierr = xf_GetKeyValue(KeyValue, "xfa", xfaFile);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  // read basis data set
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, basisFile, DataSet));
  if (ierr!=xf_OK) return ierr;

  // DataSet should contain one vectorset = BasisSet
  D = DataSet->Head;
  BasisSet = NULL;
  while (D != NULL){
    if ((D->Type == xfe_VectorSet) && (strcmp(D->Title, "BasisSet") == 0)){
      BasisSet = (xf_VectorSet *) D->Data;
      break;
    }
    D = D->Next;
  }
  if (BasisSet == NULL) return xf_Error(xf_INPUT_ERROR); // not found
  

  // If reverse flag, read root.data and compute V^T * M * U
  if (reverse){

    // Need All file to get Mass matrix
    if (!xf_NotNull(xfaFile)){
      xf_printf("Need .xfa file to perform inverse reconstruction.\n");
      return xf_Error(xf_INPUT_ERROR);
    }

    /* Create .xfa structure */
    ierr = xf_Error(xf_CreateAll(&All, xfe_False));
    if (ierr != xf_OK) return ierr;
    
    /* Read .xfa file */
    ierr = xf_Error(xf_ReadAllBinary(xfaFile, All));
    if (ierr!=xf_OK) return ierr;

    sprintf(solFile, "%s.data", solRoot);
    
    ierr = xf_Error(xf_CreateDataSet(&SolDataSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, solFile, SolDataSet));
    if (ierr!=xf_OK) return ierr;
    
    ierr = xf_Error(xf_FindPrimalState(SolDataSet, 0, &D, NULL));
    if (ierr != xf_OK) return ierr;
    U = (xf_Vector *) D->Data;
	
    // Multiply U by Mass matrix do be consistent with continuous norm
    ierr = xf_Error(xf_MultMassMatrix(All, 1.0, U));
    if (ierr != xf_OK) return ierr;

    for (i=0; i<BasisSet->nVector; i++){
      ierr = xf_Error(xf_VectorDot(BasisSet->Vector+i, U, &dp));
      if (ierr != xf_OK) return ierr;
      xf_printf("%.15E\n", dp);
    }
    
    // destroy DataSet
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;
    
    // destroy SolDataSet
    ierr = xf_Error(xf_DestroyDataSet(SolDataSet));
    if (ierr != xf_OK) return ierr;

    /* Destroy .xfa structure */
    ierr = xf_Error(xf_DestroyAll(All));
    if (ierr!=xf_OK) return ierr;


    return xf_OK; // execution terminates here
  }

  // create a solution vector, U
  ierr = xf_Error(xf_CreateVector(&U));
  if (ierr != xf_OK) return ierr;

  // copy fist BasisSet->Vector into U
  ierr = xf_CopyVector( ((All == NULL) ? NULL : All->Mesh), BasisSet->Vector+0, U);
  if (ierr != xf_OK) return ierr;
  
  // create a shell SolDataSet
  ierr = xf_Error(xf_CreateDataSet(&SolDataSet));
  if (ierr != xf_OK) return ierr;

  // add U to SolDataSet (so we can read/write it)
  U->SolverRole = xfe_SolverRolePrimalState;
  ierr = xf_Error(xf_DataSetAdd(SolDataSet, "State", xfe_Vector, xfe_True,
				(void *) U, NULL));
  if (ierr != xf_OK) return ierr;


  /* Open the .out file */
  if ((fid = fopen(outFile, "r")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  
  /* First line of .out file: nRun, N */
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d %d", &nRun, &N) != 2) return xf_Error(xf_FILE_READ_ERROR);
  
  // check N vs. BasisSet->nVector
  if (N > BasisSet->nVector){
    xf_printf("Error, number of basis vectors in basis .data file < N in .out file.\n");
    xf_printf("N = %d, nVector = %d\n", N, BasisSet->nVector);
    return xf_Error(xf_INPUT_ERROR);
  }
  else if (N < BasisSet->nVector){
    xf_printf("Warning, number of basis vectors in basis .data file != N in .out file.\n");
    xf_printf("N = %d, nVector = %d\n", N, BasisSet->nVector);
    xf_printf("Continuing.\n");
  }

  // allocate solution vector
  ierr = xf_Error(xf_Alloc((void **) &a, N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Loop over runs
  for (iRun=0; iRun<nRun; iRun++){
      
    // read line; should be "% Run .."
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (strncmp(line, "% Run", 1)!=0) return xf_Error(xf_FILE_READ_ERROR);

    // read N coefficients -> a
    for (i=0; i<N; i++){
      if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
      if (sscanf(line, "%lf", a+i) != 1) return xf_Error(xf_FILE_READ_ERROR);
    }
    
    // assemble U
    ierr = xf_Error(xf_SetZeroVector(U));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<N; i++){
      ierr = xf_Error(xf_VectorMultSet(BasisSet->Vector+i, a[i], xfe_Add, U));
      if (ierr != xf_OK) return ierr;
    }


    // write U
    sprintf(solFile, "%s_%d.data\0", solRoot, iRun);
    ierr = xf_Error(xf_WriteDataSetBinary(NULL, SolDataSet, NULL, solFile));
    if (ierr!=xf_OK) return ierr;
    
    // check for EOF
    if (feof(fid)){
      xf_printf("Warning, fewer runs in .out file than expected\n");
      nRun = iRun+1;
    }

  } // iRun

  // close outFile
  fclose(fid);
 

  
  // destroy DataSet
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;

  // destroy SolDataSet
  ierr = xf_Error(xf_DestroyDataSet(SolDataSet));
  if (ierr != xf_OK) return ierr;

  // Release memory
  xf_Release( (void *) a);

  xf_printf("xf_MRReconstruct finished.\n");

  return xf_OK;
}
