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
  FILE:  xf_EigSolve.c

  This program performs Hessian-based model reduction of initial
  conditions.

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
#include "xf_EigSolver.h"
#include "xf_Arg.h"

// minimum number of Lanczos vectors (in case only want 1 or 2 eigs)
#define MIN_LANCZOS 1
#define MAXITER_LANCZOS 100


/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int len, ierr, i, j, k, ndof, r;
  int myRank, nProc;
  int nOutput;
  int iIter, nLanczos, neig;
  enum xfe_Bool ParallelFlag, RestartFlag, done;
  enum xfe_EigStatusType Status;
  char *ArgIn[] = {"job", "NULL", ".job file name to read (run parameters)",
		   "neig", "10", "number of eigenvalues to calculate",
		   "tol", "1e-6", "relative tolerance on eigenvalues",
		   "\0"};
  char jobFile[xf_MAXSTRLEN] = "NULL"; 
  real EigTol, *E;
  xf_KeyValue KeyValueArg;
  xf_EigSolverData *EigSolverData;
  xf_Vector *U, *V, *W;
  xf_VectorSet  *LanczosSet;
  xf_All *All;

  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if (nProc > 1) ParallelFlag = xfe_True;
  else  ParallelFlag = xfe_False;


  xf_printf("\n");
  xf_printf("=== Eigensolver for general matrices (hardcoded) ===\n");
  xf_printf("\n");
  
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValueArg));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValueArg);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
      
  /* Get jobFile */
  ierr = xf_GetKeyValue(KeyValueArg, "job", jobFile);
  if (ierr != xf_OK) return ierr;
  
  /* Get number of ICs eigenvectors */
  ierr = xf_GetKeyValueInt(KeyValueArg, "neig", &neig);
  if (ierr != xf_OK) return ierr;
  if (neig <= 0) return xf_Error(xf_INPUT_ERROR); // nothing to work with

  /* Get eigenvalue tolerance */
  ierr = xf_GetKeyValueReal(KeyValueArg, "tol", &EigTol);
  if (ierr != xf_OK) return ierr;
  if (EigTol <= 0) return xf_Error(xf_INPUT_ERROR);
  
  
  // destroy key-value from arg list
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValueArg));
  if (ierr!=xf_OK) return ierr;


  /*---------------------------------------*/
  /* Read All from .job, including eqnset. */
  /* Parallelization takes place here.     */
  /*---------------------------------------*/

  ierr = xf_Error(xf_ReadAllFromJobFile(jobFile, xfe_True, &All));
  if (ierr != xf_OK) return ierr;
  
  /* Load EqnSet Library  */  
  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;
  
  /* Register EqnSet and check/set default eqnset parameters. */
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;
  
  /* Find or create an initial vector, U */
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "Restart", &RestartFlag));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_FindOrCreatePrimalState(All, RestartFlag, NULL, &U));
  if (ierr != xf_OK) return ierr;

  /* Set verbosity to low */
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", "Low"));
  if (ierr != xf_OK) return ierr;



  // Create a VectorSet for storing Lanczos Vectors
  nLanczos = max(2*neig+2, MIN_LANCZOS);
  ierr = xf_Error(xf_FindSimilarVectorSet(All, U, nLanczos+1, "LanczosSet", 
					  xfe_True, xfe_False, NULL, &LanczosSet));
  if (ierr != xf_OK) return ierr; 

  // Allocate vector for eigenvalues
  ierr = xf_Error(xf_Alloc( (void **) &E, nLanczos+1, sizeof(real)));
  if (ierr != xf_OK) return ierr;


  /*------------------------*/
  /* Start eigensolver loop */
  /*------------------------*/

  iIter = 0;
  EigSolverData = NULL;
  done = xfe_False;
  while (!done){

    xf_printf("\n== Starting Lanczos Iteration %d ==\n\n", iIter++);

    if (iIter > MAXITER_LANCZOS){
      xf_printf("Maximum number of calls to EigIterLanczos exceeded.  Stopping.\n");
      break;
    }

    // Call iteration routine
    ierr = xf_Error(xf_EigIterLanczos(NULL, neig, nLanczos, LanczosSet, EigTol, 
				      xfe_VerbosityHigh, 0, xfe_False, xfe_False, 
				      E, NULL, &V, &W, &EigSolverData, &Status));
    if (ierr != xf_OK) return ierr;

    done = ((Status == xfe_EigConverged) || (Status == xfe_EigForceExit));
    if (done) break; // converged

    /*---------------------------------*/
    /* Set W = "operator applied to V" */
    /*---------------------------------*/

    xf_printf(" Applying Mass Matrix\n");

    // W = V
    ierr = xf_Error(xf_SetVector(V, xfe_Set, W));
    if (ierr != xf_OK) return ierr;
    // W = M*W
    ierr = xf_Error(xf_MultMassMatrix(All, 1.0, W));
    if (ierr != xf_OK) return ierr;

    // Check if halted by user
    if (xf_CheckUserHalt(NULL)){
      xf_printf("WARNING: computation halted by user. Do not trust the Eigenvalues/vectors.\n");
      break;
    }
    
  } // end Lanczos eiengsolver loop


  // Print out eigenvalues
  xf_printf("\nFinished calculating eigenvalues.\n");
  xf_printf("Eigenvalues:\n");
  for (i=0; i<neig; i++) xf_printf(" Eig[%d] = %.15E\n", i, E[i]);
  xf_Release( (void *) E);



  /*----------------*/
  /* Release memory */
  /*----------------*/

  // No longer need LanczosSet
  ierr = xf_Error(xf_DestroyVectorSet(LanczosSet));
  if (ierr != xf_OK) return ierr;


  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;
  
 
  xf_printf("xf_EigSolve finished.\n");

  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
