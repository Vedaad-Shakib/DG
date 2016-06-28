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
  FILE:  xf_MROnline.c

  This program performs the online portion of model reduction.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_MeshTools.h"
#include "xf_Param.h"
#include "xf_Basis.h"
#include "xf_EqnSetHook.h"
#include "xf_Residual.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_Math.h"
#include "xf_Arg.h"
#include "xf_MRStruct.h"
#include "xf_MRCommon.h"


// for storing the run list
typedef struct
{
  int iterm; // residual term to set
  char resKey[xf_MAXSTRLEN];   // residual key to set
  char resValue[xf_MAXSTRLEN]; // residual value to set
}
xf_RunListType;


/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int len, ierr, i, j, k;
  int N, M, Mmax, iRun, nRun;
  int nNonLinear, iNewton, nNewton;
  int iterm, nOutput;
  int *IParam, *P;
  char *ArgIn[] = {"rom", "NULL", "Reduced Model (.rom) file name",
		   "inp", "online.inp", "online input file with run info",
		   "out", "online.out", "online output file",
		   "iterm", "0", "nonlinear param Res term # for a single cmd-line run",
		   "key", "NULL", "ResTerm key for cmd-line run (!NULL -> inp not used)",
		   "value", "NULL", "Corresponding ResTerm value",
		   "\0"};
  char romFile[xf_MAXSTRLEN];  
  char inpFile[xf_MAXSTRLEN];
  char outFile[xf_MAXSTRLEN];
  char resKey[xf_MAXSTRLEN], resValue[xf_MAXSTRLEN];
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line[xf_MAXLINELEN];
  FILE *fid, *finp, *fout;
  real fac, tol, Rnorm, Rnorm0, *RParam;
  real *a, *b, *da, *R, *R_a, *u, *s, *s_u, *T;
  xf_KeyValue KeyValue;
  xf_EqnSet *EqnSet;
  xf_ReducedModel *RM;
  xf_RunListType *RunList;
  
  xf_printf("\n");
  xf_printf("=== Model Reduction: Online Computation ===\n");
  xf_printf("\n");
    
      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  xf_printf("nKey = %d\n", KeyValue.nKey);
  for (i=0; i<KeyValue.nKey; i++)
    xf_printf("%d : Key = %s, Value = %s\n", i, KeyValue.Key[i], KeyValue.Value[i]);
    
  // Get romFile
  ierr = xf_GetKeyValue(KeyValue, "rom", romFile);
  if (ierr != xf_OK) return ierr;

  // Get inpFile 
  ierr = xf_GetKeyValue(KeyValue, "inp", inpFile);
  if (ierr != xf_OK) return ierr;

  // Get outFile 
  ierr = xf_GetKeyValue(KeyValue, "out", outFile);
  if (ierr != xf_OK) return ierr;

  
  // create reduced model
  ierr = xf_Error(xf_CreateReducedModel(&RM));
  if (ierr != xf_OK) return ierr;

  // read reduced model
  if ((fid = fopen(romFile, "rb")) == NULL)  return xf_Error(xf_FILE_READ_ERROR);
  ierr = xf_Error(xf_ReadReducedModelBinary(fid, RM));
  if (ierr != xf_OK) return ierr;
  if (fclose(fid)!= 0) return xf_Error(xf_FILE_READ_ERROR);


  // set EqnSet
  EqnSet = RM->EqnSet;

  // Dynamically load eqnset library
  ierr = xf_Error(xf_LoadEqnSetLibrary(EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;

  // Register the equation set
  ierr = xf_Error(xf_EqnSetRegister(EqnSet));
  if (ierr != xf_OK) return ierr;

  // Get resKey
  ierr = xf_GetKeyValue(KeyValue, "key", resKey);
  if (ierr != xf_OK) return ierr;

  if (xf_NotNull(resKey)){
    // parameters specified on command-line
    nRun = 1;

    // allocate RunList
    ierr = xf_Error(xf_Alloc((void **) &RunList, nRun, sizeof(xf_RunListType)));
    if (ierr != xf_OK) return ierr;

    // Get iterm
    ierr = xf_GetKeyValueInt(KeyValue, "iterm", &iterm);
    if (ierr != xf_OK) return ierr;

    // Get resValue
    ierr = xf_GetKeyValue(KeyValue, "value", resValue);
    if (ierr != xf_OK) return ierr;

    RunList[0].iterm = iterm;
    strcpy(RunList[0].resKey, resKey);
    strcpy(RunList[0].resValue, resValue);
  }
  else{
    // parameters specified in input file

    // open input file
    if ((finp = fopen(inpFile, "r")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    
    // read first line : nRun
    if (fgets(line, xf_MAXLINELEN, finp) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (sscanf(line, "%d", &nRun) != 1) return xf_Error(xf_FILE_READ_ERROR);
    
    // allocate RunList
    ierr = xf_Error(xf_Alloc((void **) &RunList, nRun, sizeof(xf_RunListType)));
    if (ierr != xf_OK) return ierr;
    
    for (iRun=0; iRun<nRun; iRun++){

      /* iterm */
      if (fgets(line, xf_MAXLINELEN, finp) == NULL) return xf_Error(xf_FILE_READ_ERROR);
      if (strncmp(line, "iterm", 5)!=0) return xf_Error(xf_FILE_READ_ERROR);
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      sscanf(value, "%d", &RunList[iRun].iterm);

      /* key */
      if (fgets(line, xf_MAXLINELEN, finp) == NULL) return xf_Error(xf_FILE_READ_ERROR);
      if (strncmp(line, "key", 3)!=0) return xf_Error(xf_FILE_READ_ERROR);
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      strcpy(RunList[iRun].resKey, value);

      /* value */
      if (fgets(line, xf_MAXLINELEN, finp) == NULL) return xf_Error(xf_FILE_READ_ERROR);
      if (strncmp(line, "value", 5)!=0) return xf_Error(xf_FILE_READ_ERROR);
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      strcpy(RunList[iRun].resValue, value);

      if (feof(finp)) return xf_Error(xf_FILE_READ_ERROR);
    } // iRun

    fclose(finp);

  }  


  if (EqnSet->StateRank != 1) return xf_Error(xf_NOT_SUPPORTED); // for now

  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(NULL, EqnSet, &IParam, &RParam, NULL, NULL));
  //if (ierr != xf_OK) return ierr; // Do not return with an error


  // rank of reduced model
  N = RM->N;

  // # of outputs
  nOutput = RM->nOutput;

  // allocate solution vector a
  ierr = xf_Error(xf_Alloc((void **) &a, N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // allocate da for updates
  ierr = xf_Error(xf_Alloc((void **) &da, N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // allocate output vector, b
  ierr = xf_Error(xf_Alloc((void **) &b, nOutput, sizeof(real)));
  if (ierr != xf_OK) return ierr;


  // allocate reduced residual
  ierr = xf_Error(xf_Alloc((void **) &R, N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // allocate reduced Jacobian
  ierr = xf_Error(xf_Alloc((void **) &R_a, N*N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // nNonLinear
  nNonLinear = RM->nNonLinear;

  
  xf_printf("nNonLinear = %d, N = %d\n", nNonLinear, N);

  // Mmax
  for (i=0, Mmax=0; i<nNonLinear; i++) Mmax = max(Mmax, RM->M[i]);

  // allocate u for storing D*a
  ierr = xf_Error(xf_Alloc((void **) &u, Mmax, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // allocate s for storing nonlinear evaluations
  ierr = xf_Error(xf_Alloc((void **) &s, Mmax, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // allocate s_u (derivatives)
  ierr = xf_Error(xf_Alloc((void **) &s_u, Mmax, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // allocate temporary matrix T
  ierr = xf_Error(xf_Alloc((void **) &T, Mmax*N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // allocate permutation vector P
  ierr = xf_Error(xf_Alloc((void **) &P, N, sizeof(int)));
  if (ierr != xf_OK) return ierr;


  // open output file
  if ((fout = fopen(outFile, "w")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  
  // write nRun, N
  fprintf(fout, "%d %d\n", nRun, N);


  // loop over runs
  for (iRun=0; iRun<nRun; iRun++){

    // set params in RM
    ierr = xf_SetKeyValue(RM->ResTerm[RunList[iRun].iterm].KeyValue, 
			  RunList[iRun].resKey, RunList[iRun].resValue);
    if (ierr != xf_OK) return ierr;
    
    fprintf(fout, "%% Run %d\n", iRun);
    
    xf_printf("Run %d\n", iRun);

    /*--------------*/
    /* Solve system */
    /*--------------*/

    // initialize a
    for (i=0; i<N; i++) a[i] = 0.0;
    
    // outer newton iteration
    nNewton = 20;
    tol = 1e-10;
    for (iNewton=0; iNewton<nNewton; iNewton++){

      // linear contribution to R, R_a
      for (k=0; k<N; k++) R[k] = RM->L[k];  // R = L
      xf_MxV_Add(RM->A, a, N, N, R);        // R += A*a
      for (k=0; k<N*N; k++) R_a[k] = RM->A[k];

      // nonlinear contribution to R, R_a
      for (i=0; i<nNonLinear; i++){
	/*ierr = xf_Error(xf_GetKeyValue(RM->ResTerm[iterm].KeyValue, "SourceData", resValue)); */
	/*if (ierr != xf_OK) return ierr; */
	/*xf_printf("i = %d, Value = %s\n", i, resValue); */

	M = RM->M[i];
	xf_MxV_Set(RM->D[i], a, M, N, u);  // u = D*a
      
	/* 	xf_printf("u = [\n"); */
	/* 	for (k=0; k<M; k++) xf_printf(" %.10E\n", u[k]); */
	/* 	xf_printf("]\n"); */

	if (RM->ResTerm->Type != xfe_ResTermSource) return xf_Error(xf_NOT_SUPPORTED);
	ierr = xf_Error(xf_EqnSetSourceS(EqnSet, RM->ResTerm+i, 1, IParam, RParam, M, 
					 u, NULL, NULL, NULL, NULL, s, s_u, NULL, NULL));
	if (ierr != xf_OK) return ierr;
	xf_MxV_Add(RM->E[i], s, N, M, R);  // R += E*s
	for (k=0; k<M*N; k++) T[k] = RM->D[i][k]; // T = D
	xf_ColMult(T, s_u, M, N, 1); // multiply columns of T by s_u
	xf_MxM_Add(RM->E[i], T, N, M, N, R_a); // R_a += E*T
      }

      // check for convergence
      for (k=0, Rnorm = 0.; k<N; k++) Rnorm += R[k]*R[k];
      Rnorm = sqrt(Rnorm);
      xf_printf("iNewton = %d, Rnorm = %.10E\n", iNewton, Rnorm);
      if (i == 0) Rnorm0 = Rnorm;
      if ((Rnorm < tol*Rnorm0) && (iNewton > 4)) break;

      // solve for da
      ierr = xf_Error(xf_ComputePLU(R_a, N, P));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SolvePLU(R_a, P, R, N, da, NULL));
      if (ierr != xf_OK) return ierr;

      if (iNewton < 4) fac = 0.5;
      else if (iNewton < 8) fac = 0.8;
      else fac = 1.0;
      

      for (k=0; k<N; k++) a[k] -= fac*da[k]; // a = a - da
    } // iNewton

    // write a to online.out
    for (i=0; i<N; i++) fprintf(fout, "%.15E\n", a[i]);

    // compute outputs
    xf_MxV_Set(RM->F, a, nOutput, N, b);
    for (i=0; i<nOutput; i++)
      xf_printf("Output %d: %.15E\n", i, b[i] + RM->F0[i]);


    if (Rnorm > tol){
      xf_printf("Error: Newton failed to converge on iRun = %d.\n", iRun);
      return xf_Error(xf_NOT_CONVERGED);
    }

  } // iRun

  
  // close output file
  fclose(fout);

  // release memory
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) a);
  xf_Release( (void *) b);
  xf_Release( (void *) R);
  xf_Release( (void *) R_a);
  xf_Release( (void *) P);
  xf_Release( (void *) da);
  xf_Release( (void *) u);
  xf_Release( (void *) s);
  xf_Release( (void *) s_u);
  xf_Release( (void *) T);
  xf_Release( (void *) RunList);

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  /* Destroy Reduced Model */
  ierr = xf_Error(xf_DestroyReducedModel(RM, xfe_True));
  if (ierr != xf_OK) return ierr;

  xf_printf("xf_MROnline finished.\n");

  return xf_OK;
}
