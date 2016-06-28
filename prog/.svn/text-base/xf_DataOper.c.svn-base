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
  FILE:  xf_DataOper.c

  This program performs a mathematical operation on two data vectors

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_MeshTools.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Math.h"
#include "xf_Quad.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Arg.h"
#include "xf_Solver.h"
#include "xf_EqnSetHook.h"


/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, i1, i2;
  int myRank, nProc;
  int OrderInc;
  int ibatch[3], iData, nData, ind;
  enum xfe_Bool dobatch = xfe_False;
  enum xfe_Bool Analytical, TakeDiff;
  enum xfe_Bool DotFlag = xfe_False;
  enum xfe_AddType oper;
  char *ArgIn[] = {"data1", "NULL", "first .data file",
		   "data2", "NULL", "second .data file",
		   "factor1", "1.0", "multiplicative factor for data1",
		   "factor2", "1.0", "multiplicative factor for data2",
		   "oper", "Add", "Add/Sub/Set/Neg/Dot (data2 used for Set/Neg)",
		   "xfa", "NULL", ".xfa file (optional)",
		   "UseMass", "False", "True to use Mass matrix in computation",
		   "Refine", "True", "True to order refine data1->out",
		   "Reconstruct", "False", "True to reconstruct data1->out",
		   "OrderInc", "1", "order increment for reconstructdion",
		   "VariableSet", "None", "Convert to new variable set, data1->out",
		   "batch", "NULL", "<start step end> string for batch mode",
		   "out", "out.data", "output .data file",
		   "\0"};
  char xfaFile[xf_MAXSTRLEN]; 
  char dataFile[xf_MAXSTRLEN];
  char data1File[xf_MAXSTRLEN];
  char data2File[xf_MAXSTRLEN];
  char outFile0[xf_MAXSTRLEN];
  char outFile[xf_MAXSTRLEN];
  char soper[xf_MAXSTRLEN];
  char line[xf_MAXLINELEN];
  char VariableSet[xf_MAXLINELEN];
  char batch[xf_MAXSTRLEN];
  enum xfe_Bool UseMass = xfe_False;
  enum xfe_Bool Reconstruct = xfe_False;
  enum xfe_Bool Refine = xfe_False;
  enum xfe_DataType DataType;
  real dp;
  real factor1, factor2;
  FILE *fid;
  xf_KeyValue KeyValue;
  xf_Vector *U1, *U2;
  xf_VectorSet *V1, *V2;
  xf_DataSet *DataSet1, *DataSet2;
  xf_Data *D1, *D2, *D;
  xf_All *All = NULL;
  xf_Mesh *Mesh = NULL;

  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  
  xf_printf("\n");
  xf_printf("=== Data Operation ===\n");
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

  // Get data1
  ierr = xf_GetKeyValue(KeyValue, "data1", data1File);
  if (ierr != xf_OK) return ierr;

  // Get data2
  ierr = xf_GetKeyValue(KeyValue, "data2", data2File);
  if (ierr != xf_OK) return ierr;

  // Get outFile
  ierr = xf_GetKeyValue(KeyValue, "out", outFile0);
  if (ierr != xf_OK) return ierr;

  // Get oper
  ierr = xf_GetKeyValue(KeyValue, "oper", soper);
  if (ierr != xf_OK) return ierr;

  // Get UseMass
  ierr = xf_GetKeyValueBool(KeyValue, "UseMass", &UseMass);
  if (ierr != xf_OK) return ierr;

  // Get Reconstruct
  ierr = xf_GetKeyValueBool(KeyValue, "Reconstruct", &Reconstruct);
  if (ierr != xf_OK) return ierr;

  // Get Refine
  ierr = xf_GetKeyValueBool(KeyValue, "Refine", &Refine);
  if (ierr != xf_OK) return ierr;

  // Get OrderInc
  ierr = xf_GetKeyValueInt(KeyValue, "OrderInc", &OrderInc);
  if (ierr != xf_OK) return ierr;

  // Get factor1
  ierr = xf_GetKeyValueReal(KeyValue, "factor1", &factor1);
  if (ierr != xf_OK) return ierr;

  // Get factor2
  ierr = xf_GetKeyValueReal(KeyValue, "factor2", &factor2);
  if (ierr != xf_OK) return ierr;

  // Get VariableSet
  ierr = xf_GetKeyValue(KeyValue, "VariableSet", VariableSet);
  if (ierr != xf_OK) return ierr;

  // Get batch
  ierr = xf_GetKeyValue(KeyValue, "batch", batch);
  if (ierr != xf_OK) return ierr;

  /* Get operation */
  if (strncmp(soper, "Dot", 3) == 0)
    DotFlag = xfe_True;
  else{
    ierr = xf_Error(xf_GetKeyValueEnum(KeyValue, "oper", xfe_AddName, 
				       xfe_AddLast, (int *) &oper));
    if (ierr != xf_OK) return ierr;
  }

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;

  // determine if doing batch mode
  dobatch = xfe_False;
  nData = 1;
  if (xf_NotNull(batch)){
    ierr = xf_Error(xf_ScanInt(batch, 3, ibatch));
    if (ierr != xf_OK) return ierr;
    dobatch = xfe_True;
    nData = 1 + (ibatch[2]-ibatch[0])/ibatch[1];
  }


  /* Read .xfa if requested */
  if (xf_NotNull(xfaFile)){
    /* Create .xfa structure */
    ierr = xf_Error(xf_CreateAll(&All, xfe_False));
    if (ierr != xf_OK) return ierr;
    
    /* Read .xfa file*/
    ierr = xf_Error(xf_ReadAllBinary(xfaFile, All));
    if (ierr!=xf_OK) return ierr;
    Mesh = All->Mesh;
  }


  
  // loop over data files (allows for batch mode)
  for (iData=0; iData<nData; iData++){
    
    // data file name
    if (dobatch){
      ind = ibatch[0] + ibatch[1]*iData;    // index of data file
      sprintf(dataFile, "%s%d.data", data1File, ind);
      sprintf( outFile, "%s%d.data",  outFile0, ind);
    }
    else{
      sprintf(dataFile, "%s", data1File);
      sprintf( outFile, "%s",  outFile0);
    }


    /* Read data1 */
    xf_printf("Reading first .data file\n");
    ierr = xf_Error(xf_CreateDataSet(&DataSet1));
    if (ierr != xf_OK) return ierr;
  
    ierr = xf_Error(xf_ReadDataSetBinary(Mesh, NULL, dataFile, DataSet1));
    if (ierr!=xf_OK) return ierr;
  
    // use first piece of data
    D1 = DataSet1->Head;
    DataType = D1->Type;
  

    /* Read data2 */
    if (xf_NotNull(data2File)){
      xf_printf("Reading second .data file\n");
      ierr = xf_Error(xf_CreateDataSet(&DataSet2));
      if (ierr != xf_OK) return ierr;
    
      ierr = xf_Error(xf_ReadDataSetBinary(Mesh, NULL, data2File, DataSet2));
      if (ierr!=xf_OK) return ierr;
    
      // use first piece of data
      D2 = DataSet2->Head;
      if (D2->Type != DataType) return xf_Error(xf_INPUT_ERROR);
    }
    else{
      DataSet2 = NULL;
      D2 = NULL;
    }
  
    if (xf_NotNull(VariableSet)){
      // Variable change takes precedence
      if (All == NULL) return xf_Error(xf_INPUT_ERROR); // need xfa
    
      // load dynamic library
      ierr = xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary);
      if (ierr != xf_OK) return ierr;
      // register equation set
      ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
      if (ierr != xf_OK) return ierr;

      // loop over data in DataSet1
      D = DataSet1->Head;
      while (D != NULL){
	if (D->Type == xfe_Vector){
	  U1 = (xf_Vector *) D->Data;

	  // transform vector, overwrite U1
	  ierr = xf_Error(xf_ChangeVariableSet(All, VariableSet, U1));
	  if (ierr != xf_OK) return ierr;
	
	  // apply factor1
	  ierr = xf_Error(xf_VectorMult(U1, factor1));
	  if (ierr != xf_OK) return ierr;

	}
	D = D->Next;
      }

      // write outFile
      ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet1, NULL, outFile));
      if (ierr != xf_OK) return ierr;

    }
    else if (Refine){
      if (DataType != xfe_Vector) return xf_Error(xf_NOT_SUPPORTED);
      U1 = (xf_Vector *) D1->Data;
    
      // increment order
      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(NULL, NULL, U1, NULL, 
                                                             xfe_BasisLast, OrderInc));
      if (ierr != xf_OK) return ierr;
      
      // save U1 to outFile
      ierr = xf_Error(xf_DumpVectorBinary(NULL, D1->Title, U1, outFile));
      if (ierr != xf_OK) return ierr;

    }
    else if (Reconstruct){
      if (DataType != xfe_Vector) return xf_Error(xf_NOT_SUPPORTED);
      U1 = (xf_Vector *) D1->Data;
    
      ierr = xf_Error(xf_ReconstructVector(All, OrderInc, U1, xfe_True, &U2));
      if (ierr != xf_OK) return ierr;
     
      // save U2 to outFile
      ierr = xf_Error(xf_DumpVectorBinary(NULL, "Vector", U2, outFile));
      if (ierr != xf_OK) return ierr;

      // destroy U2
      ierr = xf_Error(xf_DestroyVector(U2, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    else if (DataType == xfe_Vector){
      U1 = (xf_Vector *) D1->Data;
      U2 = (xf_Vector *) D2->Data;

      ierr = xf_Error(xf_VectorMult(U1, factor1));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_VectorMult(U2, factor2));
      if (ierr != xf_OK) return ierr;

      // perform the operation 
      if (DotFlag){
	ierr = xf_Error(xf_VectorDot(U1, U2, &dp));
	if (ierr != xf_OK) return ierr;
	xf_printf("Dot product = %.15E\n", dp);
      }
      else{
	// standard vector operation (result is stored in U1)
	ierr = xf_Error(xf_VectorMultSet(U2, 1.0, oper, U1));
	if (ierr != xf_OK) return ierr;  
      }

      // save U1 to outFile
      ierr = xf_Error(xf_DumpVectorBinary(NULL, "Vector", U1, outFile));
      if (ierr != xf_OK) return ierr;
    }
    else if (DataType == xfe_VectorSet){
      V1 = (xf_VectorSet *) D1->Data;
      V2 = (xf_VectorSet *) D2->Data;
      if (DotFlag){
	if (UseMass) xf_printf("V2^T * M * V1 =\n");
	else xf_printf("V2^T * V1 =\n");

	if (UseMass){
	  // V1 = M * V1
	  for (i1=0; i1<V1->nVector; i1++){
	    ierr = xf_Error(xf_MultMassMatrix(All, 1.0, V1->Vector+i1));
	    if (ierr != xf_OK) return ierr;
	  }
	}

	// answer(i2,i1) = V2(i2) dot V1(i1)      
	for (i2=0; i2<V2->nVector; i2++){		
	  for (i1=0; i1<V1->nVector; i1++){
	    ierr = xf_Error(xf_VectorDot(V2->Vector+i2, V1->Vector+i1, 
					 &dp));
	    if (ierr != xf_OK) return ierr;
	    xf_printf("%.15E ", dp);
	  }
	  xf_printf("\n");
	}
      }
      else return xf_Error(xf_NOT_SUPPORTED);

    }
    else return xf_Error(xf_NOT_SUPPORTED);


    // destroy DataSet1
    ierr = xf_Error(xf_DestroyDataSet(DataSet1));
    if (ierr != xf_OK) return ierr;
    
    // destroy DataSet2
    ierr = xf_Error(xf_DestroyDataSet(DataSet2));
    if (ierr != xf_OK) return ierr;

  } // iData in batch mode

  if (All != NULL){
    /* Destroy .xfa structure */
    ierr = xf_Error(xf_DestroyAll(All));
    if (ierr!=xf_OK) return ierr;
  }
    
  xf_printf("xf_DataOper finished.\n");

  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
