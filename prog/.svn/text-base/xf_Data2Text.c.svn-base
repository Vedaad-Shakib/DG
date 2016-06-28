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
  FILE:  xf_Data2Text.c

  Writes a plain text file format for a given .data file

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



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr;
  int ind, iData, nData;
  int ibatch[3];
  int iVector, nVector;
  enum xfe_Bool dobatch;
  char *ArgIn[] = {"inroot", "NULL", "will read <inroot>.data file",
		   "outroot", "NULL", "will write <outroot>.txt file (NULL->inroot)",
		   "batch", "NULL", "<start step end> string for batch mode",
		   "\0"};
  char inroot[xf_MAXSTRLEN];
  char dataFile[xf_MAXSTRLEN];
  char outroot[xf_MAXSTRLEN];
  char outbase[xf_MAXSTRLEN];
  char outFile[xf_MAXSTRLEN];
  char batch[xf_MAXSTRLEN];
  xf_KeyValue KeyValue;
  xf_Vector *U;
  xf_DataSet *DataSet;
  xf_Data *D;

  
  xf_printf("\n");
  xf_printf("=== Data to text conversion ===\n");
  xf_printf("\n");
    
      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  // Get data
  ierr = xf_GetKeyValue(KeyValue, "inroot", inroot);
  if (ierr != xf_OK) return ierr;

  // Get outFile
  ierr = xf_GetKeyValue(KeyValue, "outroot", outroot);
  if (ierr != xf_OK) return ierr;

  // Get batch
  ierr = xf_GetKeyValue(KeyValue, "batch", batch);
  if (ierr != xf_OK) return ierr;

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

  // set outroot to inroot if NULL
  if (!xf_NotNull(outroot)) sprintf(outroot, "%s\0", inroot);

  // loop over data files (allows for batch mode)
  for (iData=0; iData<nData; iData++){

    // data file name
    if (dobatch){
      ind = ibatch[0] + ibatch[1]*iData;    // index of data file
      sprintf(dataFile, "%s%d.data", inroot, ind);
      sprintf(outbase, "%s%d", outroot, ind);
    }
    else{
      sprintf(dataFile, "%s.data", inroot);
      sprintf(outbase, "%s", outroot);
    }

    /* Read data */
    xf_printf("Reading %s\n", dataFile);
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, dataFile, DataSet));
    if (ierr!=xf_OK) return ierr;
  
    // count number of vectors
    D = DataSet->Head;
    nVector=0;
    while (D != NULL){
      if (D->Type == xfe_Vector) nVector++;
      D = D->Next;
    }

    if (nVector == 0) return xf_Error(xf_INPUT_ERROR);

    // write out each vector
    iVector=0;
    D = DataSet->Head;
    while (D != NULL){
      if (D->Type != xfe_Vector) continue;
      if (nVector > 1)
	sprintf(outFile, "%s_%d.txt", outbase, iVector);
      else
	sprintf(outFile, "%s.txt", outbase);

      U = (xf_Vector *) D->Data;
      ierr = xf_Error(xf_VectorTextOut(outFile, xfe_True, U));
      if (ierr != xf_OK) return ierr;

      iVector++;
      D = D->Next;
    }

  } // iData


  // destroy DataSet
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;
    
  xf_printf("xf_Data2Text finished.\n");

  return xf_OK;
}
