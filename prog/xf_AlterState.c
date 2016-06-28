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
  FILE:  xf_AlterState.c

  This program alters a state based on a given equation-specific function.

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
#include "xf_Arg.h"
#include <time.h>
#include <stdlib.h>



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  // arguments read in on the command line
  char *ArgIn[] = {"xfa", "NULL", ".xfa file containing mesh and primal state",
		   "out", "NULL", ".data file to which modified state is written",
		   "FcnName", "NULL", "Alteration function name",
		   "FcnParam", "NULL", "function parameters (as string)",
		   "\0"};
  int ierr;
  int nFcnParam;
  char xfaFile[xf_MAXSTRLEN] = "NULL"; 
  char outFile[xf_MAXSTRLEN] = "NULL";
  char FcnName[xf_MAXSTRLEN] = "NULL";
  char FcnParamString[xf_MAXSTRLEN] = "NULL";

  int *IParam = NULL;
  real *RParam = NULL;

  xf_KeyValue KeyValueArg;
  xf_Data *D;
  xf_Vector *U;
  xf_All *All;


  xf_printf("\n");
  xf_printf("=== State alteration ===\n");
  xf_printf("\n");
    
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValueArg));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValueArg);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  /* Get xfaFile */
  ierr = xf_GetKeyValue(KeyValueArg, "xfa", xfaFile);
  if (ierr != xf_OK) return ierr;

  /* Get outFile */
  ierr = xf_GetKeyValue(KeyValueArg, "out", outFile);
  if (ierr != xf_OK) return ierr;

  /* Get FcnName */
  ierr = xf_GetKeyValue(KeyValueArg, "FcnName", FcnName);
  if (ierr != xf_OK) return ierr;

  /* Get FcnParam */
  ierr = xf_GetKeyValue(KeyValueArg, "FcnParam", FcnParamString);
  if (ierr != xf_OK) return ierr;
  
  // destroy key-value from arg list
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValueArg));
  if (ierr!=xf_OK) return ierr;

  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  /* Read .xfa file */
  xf_printf("Reading xfa file: %s\n", xfaFile);
  ierr = xf_Error(xf_ReadAllBinary(xfaFile, All));
  if (ierr!=xf_OK) return ierr;

  // load dynamic library
  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;

  // register EqnSet
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;

  // Get State, U
  ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D, NULL));
  if (ierr != xf_OK) return ierr;
  U = (xf_Vector *) D->Data;

  // call function to modify state
  ierr = xf_Error(xf_AlterState(All, FcnName, FcnParamString, U));
  if (ierr != xf_OK) return ierr;

  // write out modified state
  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, D->Title, U, outFile));
  if (ierr != xf_OK) return ierr;
  

  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;

  // Release memory
  xf_Release( (void *) RParam);
  xf_Release( (void *) IParam);
  

  xf_printf("xf_AlterState finished.\n");


  return xf_OK;

}
