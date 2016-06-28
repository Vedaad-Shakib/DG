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
  FILE:  xf_DataPeel.c

  This program peels off a data vector from an xfa

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Arg.h"


/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i;
  int start, end, iv[2], nVector0;
  enum xfe_Bool found;
  char *ArgIn[] = {"title", "NULL", "title of data to peel",
		   "out", "NULL", "data will be saved in <out>",
		   "outtitle", "NULL", "if not null, new name for data (optional)",
		   "xfa", "NULL", ".xfa file containing the all structure (option 1)",
		   "dataset", "NULL", ".data set file containing vectors (option 2)",
		   "vecsetrange", "NULL", "pull off these vectors from set (0-based)",
		   "\0"};
  char outFile[xf_MAXSTRLEN];
  char outtitle[xf_MAXSTRLEN] = "NULL";
  char xfaFile[xf_MAXSTRLEN];
  char datasetFile[xf_MAXSTRLEN] = "NULL";
  char vecsetrange[xf_MAXSTRLEN] = "NULL";
  char title[xf_MAXSTRLEN];
  char line[xf_MAXLINELEN];
  char *newtitle = NULL;
  xf_KeyValue KeyValue;
  xf_DataSet *DataSet;
  xf_DataSet *DataSetNew;
  xf_Data *D, *DataU;
  xf_Vector *U, *Vector0;
  xf_VectorSet *VS;
  xf_All *All;
  
  xf_printf("\n");
  xf_printf("=== Data Peeling ===\n");
  xf_printf("\n");
    
      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);

  // Get title
  ierr = xf_GetKeyValue(KeyValue, "title", title);
  if (ierr != xf_OK) return ierr;

  // Get out
  ierr = xf_GetKeyValue(KeyValue, "out", outFile);
  if (ierr != xf_OK) return ierr;

  // Get outtitle
  ierr = xf_GetKeyValue(KeyValue, "outtitle", outtitle);
  if (ierr != xf_OK) return ierr;

  // Get xfaFile
  ierr = xf_GetKeyValue(KeyValue, "xfa", xfaFile);
  if (ierr != xf_OK) return ierr;

  // Get datasetFile
  ierr = xf_GetKeyValue(KeyValue, "dataset", datasetFile);
  if (ierr != xf_OK) return ierr;

  // Get vecsetrange
  ierr = xf_GetKeyValue(KeyValue, "vecsetrange", vecsetrange);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  if (xf_NotNull(datasetFile)){
    // create dataset for reading
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;
  
    // read in data
    xf_printf("Reading data set: %s\n", datasetFile);
    ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, datasetFile, DataSet));
    if (ierr != xf_OK) return ierr;
  }
  else{
    /* Create .xfa structure */
    ierr = xf_Error(xf_CreateAll(&All, xfe_False));
    if (ierr != xf_OK) return ierr;
    
    /* Read .xfa file */
    xf_printf("Reading xfa file: %s\n", xfaFile);
    ierr = xf_Error(xf_ReadAllBinary(xfaFile, All));
    if (ierr!=xf_OK) return ierr;

    DataSet = All->DataSet; // just a pointer
  }

  /* Locate title */
  found = xfe_False;
  D = DataSet->Head;
  while (D != NULL){
    if ((strcmp(D->Title, title) == 0)){
      found = xfe_True;
      break;
    }
    D = D->Next;
  }
  if (!found){
    xf_printf("Could not find title = %s in Dataset.\n", title);
    return xf_Error(xf_NOT_FOUND);
  }

  // new title
  newtitle = (xf_NotNull(outtitle) ? outtitle : D->Title);

  // set timeindex if necessary
  if (D->Type == xfe_Vector){
    U = (xf_Vector *) D->Data;
    U->TimeIndex = 0;
    U->SolverRole = xfe_SolverRolePrimalState;
  }

  // restrict output to just a few vectors if requested
  if ((D->Type == xfe_VectorSet) && (xf_NotNull(vecsetrange))){
    VS = (xf_VectorSet *) D->Data;
    ierr = xf_Error(xf_ScanInt(vecsetrange, 2, iv));
    if (ierr != xf_OK) return ierr;
    start = iv[0];
    end = iv[1];
    nVector0 = VS->nVector;
    Vector0 = VS->Vector;
    VS->nVector = end-start+1;
    VS->Vector = VS->Vector+start;
  }


  /* Create a wrapper dataset for writing */
  ierr = xf_Error(xf_CreateDataSet(&DataSetNew));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DataSetAdd(DataSetNew, newtitle, D->Type, D->ReadWrite, D->Data, &DataU));
  if (ierr != xf_OK) return ierr;

  
  /* Write data set in binary form */
  xf_printf("Writing data file: %s\n", outFile);
  ierr = xf_Error(xf_WriteDataSetBinary(NULL, DataSetNew, NULL, outFile));
  if (ierr != xf_OK) return ierr;

  // reset changes to vector set
  if ((D->Type == xfe_VectorSet) && (xf_NotNull(vecsetrange))){
    VS = (xf_VectorSet *) D->Data;
    VS->nVector = nVector0;
    VS->Vector = Vector0;
  }

  // destroy DataSetNew
  DataU->Data = NULL;
  ierr = xf_Error(xf_DestroyDataSet(DataSetNew));
  if (ierr != xf_OK) return ierr;
    
  if (xf_NotNull(datasetFile)){
    // destroy DataSet
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;
  }
  else{
    /* Destroy .xfa structure */
    ierr = xf_Error(xf_DestroyAll(All));
    if (ierr!=xf_OK) return ierr;
  }

  xf_printf("xf_DataPeel finished.\n");

  return xf_OK;
}
