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
  FILE:  xf_ICSensitivity.c

  This program calculates the sensitivity of an initial condition
  field to a parameter governing the initial conditions.  The initial
  condition is assumed to be parameterized by inputs specified in the
  equation-set file.  This program computes the change in the initial
  condition that results from a perturbation of one parameter, using
  finite differences.

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
//   FUNCTION Definition: xf_CalculateSensitivity
static int
xf_CalculateSensitivity(xf_All *All, char *ParamName, real eps, 
			xf_Vector **pdU)
{
/*
PURPOSE:

  Calculates state perturbation that results from a parameter perturbation of epsilon.

INPUTS:

  All : All structure
  ParamName : name of parameter
  eps : how much are we perturbing the parameter?
  
OUTPUTS: 

  (*pdU) : state vector perturbation


RETURN: Error code

*/
  int ierr;
  int nICParam;
  int iHeader, nHeader;
  enum xfe_Bool found = xfe_False;
  xf_IC *IC = NULL;
  xf_EqnSet *EqnSet;
  char **ICHeader = NULL;
  char s[xf_MAXSTRLEN];
  real *ICParam = NULL;
  xf_Vector *U = NULL;

  EqnSet = All->EqnSet;

  // sanity checks
  if ((EqnSet->ICs == NULL) || (EqnSet->ICs[0].nIC <= 0)) 
    return xf_Error(xf_INPUT_ERROR);
  
  IC = EqnSet->ICs[0].IC;
  
  /* Create vector with appropriate initial condition, U */
  ierr = xf_Error(xf_FindOrCreatePrimalState(All, xfe_False, NULL, &U));
  if (ierr != xf_OK) return ierr;
  
  /* Initialize state */
  ierr = xf_Error(xf_InitState(All, U));
  if (ierr != xf_OK) return ierr;

  // pull off initial condition header
  ierr = xf_Error(xf_ScanXStringAlloc(IC->Header, xf_MAXSTRLEN, &nHeader, &ICHeader));
  if (ierr != xf_OK) return ierr;

  // Pull off initial condition parameters
  ierr = xf_Error(xf_ScanXRealAlloc(IC->Data, &nICParam, &ICParam));
  if (ierr != xf_OK) return ierr;

  if (nICParam != nHeader) return xf_Error(xf_INPUT_ERROR);

  // set IC parameters in All->EqnSet
  for (iHeader=0, found = xfe_False; iHeader<nHeader; iHeader++){
    if (strncmp(ICHeader[iHeader], ParamName, strlen(ParamName)) == 0){
      ICParam[iHeader] += eps;
      found = xfe_True;
      break;
    }
  } // iHeader

  if (!found) return xf_Error(xf_NOT_FOUND);

  // store data into a string: IC->Data
  xf_Release((void *) IC->Data);
  ierr = xf_Error(xf_Alloc( (void **) &IC->Data, 20*nHeader, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  sprintf(IC->Data, "%.10E ", ICParam[0]);
  for (iHeader=1; iHeader<nHeader; iHeader++){
    sprintf(s, "%.10E ", ICParam[iHeader]);
    strcat(IC->Data, s);
  } // iHeader

  xf_Release2( (void **) ICHeader);
  xf_Release ( (void  *) ICParam);


  // find delta U vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_True,
				       xfe_True, NULL, pdU, NULL));
  /* Initialize state */
  ierr = xf_Error(xf_InitState(All, (*pdU)));
  if (ierr != xf_OK) return ierr;

  /* Subtract off U from dU*/
  ierr = xf_Error(xf_SetVector(U, xfe_Sub, (*pdU)));
  if (ierr != xf_OK) return ierr;
   
  /* Scale to get sensitivity */
  ierr = xf_Error(xf_VectorMult((*pdU), 1.0/eps));
  if (ierr != xf_OK) return ierr;
      
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  // arguments read in on the command line
  char *ArgIn[] = {"xfa", "NULL", "All file to read",
		   "eqn", "NULL", "Alternate equation set file to read (optional)",
		   "ParamName", "NULL", "name of parameter for sensitivity calc",
		   "eps", "1e-4", "epsilon for finite difference calc",
		   "out", "dU.data", "name of file for sensitivity (1/eps included)",
		   "\0"};
  int ierr;
  char xfaFile[xf_MAXSTRLEN]   = "NULL"; 
  char eqnFile[xf_MAXSTRLEN]   = "NULL";
  char outFile[xf_MAXSTRLEN]   = "NULL";
  char ParamName[xf_MAXSTRLEN] = "NULL";
  real eps;
  xf_KeyValue KeyValueArg;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Vector *dU = NULL;  // sensitivity vector
  xf_All *All;
  xf_EqnSet *EqnSet;

  xf_printf("\n");
  xf_printf("===  Initial Condition Sensitivity Calculation ===\n");
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

  /* Get eqnFile */
  ierr = xf_GetKeyValue(KeyValueArg, "eqn", eqnFile);
  if (ierr != xf_OK) return ierr;

  /* Get outFile */
  ierr = xf_GetKeyValue(KeyValueArg, "out", outFile);
  if (ierr != xf_OK) return ierr;

  /* Get ParamName*/
  ierr = xf_GetKeyValue(KeyValueArg, "ParamName", ParamName);
  if (ierr != xf_OK) return ierr;

  /* epsilon value */
  ierr = xf_Error(xf_GetKeyValueReal(KeyValueArg, "eps", &eps));
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

  /* Read EqnSetFile if specified */
  if (xf_NotNull(eqnFile)){

    xf_printf("Reading alternate EqnSet file: %s\n", eqnFile);

    ierr = xf_Error(xf_CreateEqnSet(&EqnSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadEqnSetFile(eqnFile, NULL, EqnSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DestroyEqnSet(All->EqnSet, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    All->EqnSet = EqnSet;
    All->EqnSet->Dim = All->Mesh->Dim;
  }

  /* Load dynamic EqnSet Library  */  
  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;
  
  /* Register EqnSet and check/set default eqnset parameters. */
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;
    


  // Compute sensitivity
  ierr = xf_Error(xf_CalculateSensitivity(All, ParamName, eps, &dU));
  if (ierr != xf_OK) return ierr;

  // Write out vector dU
  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "dU", dU, outFile));
  if (ierr != xf_OK) return ierr;


  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;


  xf_printf("xf_ICSensitivity finished.\n");

  return xf_OK;
}
