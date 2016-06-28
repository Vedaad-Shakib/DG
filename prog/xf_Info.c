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
  FILE:  xf_Info.c

  This program prints info about a .xfa file.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"


/******************************************************************/
//   FUNCTION Definition: xf_MeshInfo
static int 
xf_MeshInfo(xf_Mesh *Mesh, int Verbosity)
{
  int ibfgrp, egrp;

  xf_printf("---- Mesh Info ----\n");
  
  xf_printf(" Dim         = %d\n", Mesh->Dim);
  xf_printf(" nNode       = %d\n", Mesh->nNode);
  xf_printf(" nIFace      = %d\n", Mesh->nIFace);
  xf_printf(" nBFaceGroup = %d\n", Mesh->nBFaceGroup);
  xf_printf(" %8s %8s %s\n", "[Group]", "[nBFace]", "[Title]");
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    xf_printf(" %8d %8d %s\n", ibfgrp, Mesh->BFaceGroup[ibfgrp].nBFace,
	      Mesh->BFaceGroup[ibfgrp].Title);
  }
  xf_printf(" nElemGroup  = %d\n", Mesh->nElemGroup);
  xf_printf(" %8s %15s %8s %7s %7s\n", "[Group]", "[QBasis]", "[QOrder]", 
	    "[nElem]", "[nNode]");
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    xf_printf(" %8d %15s %8d %7d %7d\n", egrp, 
	      xfe_BasisName[Mesh->ElemGroup[egrp].QBasis],
	      Mesh->ElemGroup[egrp].QOrder, Mesh->ElemGroup[egrp].nElem,
	      Mesh->ElemGroup[egrp].nNode);
  }
  xf_printf(" Mesh Motion = ");
  if (Mesh->Motion != NULL){
    xf_printf("%s ", xfe_MotionName[Mesh->Motion->Type]);
    if (Mesh->Motion->Active) xf_printf("(Active)\n");
    else xf_printf("(Not active)\n");
  }
  else{
    xf_printf("NULL\n");
  }
  xf_printf("\n");
  return xf_OK;
}

 
/******************************************************************/
//   FUNCTION Definition: xf_GeomInfo
static int 
xf_GeomInfo(xf_Geom *Geom, int Verbosity)
{
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParamInfo
static int 
xf_ParamInfo(xf_Param *Param, int Verbosity)
{
  xf_printf("---- Param Info ----\n");
  xf_DumpKeyValue(Param->KeyValue, " %s = %s\n", stdout);
  xf_printf("\n");
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetInfo
static int 
xf_EqnSetInfo(xf_EqnSet *EqnSet, int Verbosity)
{
  int i, n;
  xf_Output *Output;

  xf_printf("---- EqnSet Info ----\n");
  xf_printf(" EqnSetLibrary = %s\n", EqnSet->EqnSetLibrary);
  xf_printf(" Parameters:\n");
  xf_DumpKeyValue(EqnSet->KeyValue, "  %s = %s\n", stdout);
  xf_printf(" nResTerms = ");
  if (EqnSet->ResTerms != NULL){
    xf_printf("%d", (n = EqnSet->ResTerms[0].nResTerm));
    for (i=0; i<n; i++)
      xf_printf(" [%s]", xfe_ResTermName[EqnSet->ResTerms[0].ResTerm[i].Type]);
    xf_printf("\n");
  }
  else xf_printf("NULL\n");
  xf_printf(" IC[0] = ");
  if ((EqnSet->ICs != NULL) && (EqnSet->ICs[0].nIC > 0))
    xf_printf("%s, %s\n", EqnSet->ICs[0].IC[0].Type, EqnSet->ICs[0].IC[0].Data);
  else xf_printf("NULL\n");
  xf_printf(" nBC = ");
  if (EqnSet->BCs != NULL){
    xf_printf("%d\n", (n = EqnSet->BCs[0].nBC));
    for (i=0; i<n; i++)
      xf_printf("  %15s : %s, %s\n", EqnSet->BCs[0].BC[i].BFGTitle,
		EqnSet->BCs[0].BC[i].Type, EqnSet->BCs[0].BC[i].Data);
    xf_printf("\n");
  }
  else xf_printf("NULL\n");
  xf_printf(" nOutput = ");
  if (EqnSet->Outputs != NULL){
    xf_printf("%d\n", (n = EqnSet->Outputs[0].nOutput));
    Output = EqnSet->Outputs[0].Output;
    for (i=0; i<n; i++)
      xf_printf("  %15s : %s, UsesFlux = %s, nBFG = %d\n", Output[i].Name, 
		xfe_OutputName[Output[i].Type],	xfe_BoolName[Output[i].UsesFlux], 
		Output[i].nBFG);
  }
  else xf_printf("NULL\n");

  xf_printf("\n");

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int len, ierr, Verbosity;

  char InputFile[xf_MAXSTRLEN];
  char *pext;
  xf_DataSet *DataSet;
  xf_All *All;

  xf_printf("\n");
  xf_printf("=== xf_Info ===\n");
  xf_printf("\n");

  /* Check number of arguments */
  if ( (argc != 2) && (argc != 3)){
    xf_printf("Usage:\n");
    xf_printf("xf_Info <xfafile> [<Verbosity>]\n");
    xf_printf("\n");
    xf_printf("<xfafile> is the name of the xfa file to read.\n");
    xf_printf("<Verbosity> = 0 or 1 is the level of output detail. Optional. 0 is default.\n");
    xf_printf("\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  Verbosity = 0;
  if ( (argc == 3) && (sscanf(argv[2], "%d", &Verbosity) != 1) ) return xf_Error(xf_INPUT_ERROR);
  if ( (Verbosity != 0) && (Verbosity != 1) ) return xf_Error(xf_INPUT_ERROR);

  /* Get InputFile name */
  strcpy(InputFile, argv[1]);
  len = strlen(InputFile);
  pext = InputFile + len - 4; // pointer to extension
  if (len < 4){
    xf_printf("Error, InputFile requires extension.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  if ( (len >= 5) && ((strncmp(pext-1, ".data", 5) == 0)) ){
    // asking for info on .data file
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;

    xf_printf("Reading in %s\n", InputFile);
      
    // read in .data file
    ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, InputFile, DataSet));
    if (ierr != xf_OK) return ierr;
    
    // print out info on dataset
    ierr = xf_Error(xf_DataSetInfo(DataSet));
    if (ierr != xf_OK) return ierr;

    // destroy DataSet
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;

    // done
    return xf_OK;
  
  }
  else if ((strncmp(pext, ".xfa", 4) != 0)){
    xf_printf("Error, unrecognized extension on InputFile.\n");
    return xf_Error(xf_INPUT_ERROR);
  }

  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Read .xfa file*/
  ierr = xf_Error(xf_ReadAllBinary(InputFile, All));
  if (ierr!=xf_OK) return ierr;


  /* Print info for each chunk */
  ierr = xf_Error(xf_MeshInfo(All->Mesh, Verbosity));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_GeomInfo(All->Geom, Verbosity));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DataSetInfo(All->DataSet));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_EqnSetInfo(All->EqnSet, Verbosity));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_ParamInfo(All->Param, Verbosity));
  if (ierr != xf_OK) return ierr;


  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;

  xf_printf("xf_Info finished.\n");

  return xf_OK;
}
