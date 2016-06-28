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
  FILE:  xf_MRCommon.c

  This file contains common routines for model reduction.

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
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_Math.h"
#include "xf_MRStruct.h"


/******************************************************************/
//   FUNCTION Definition: xf_AllocIPoint
int 
xf_AllocIPoint(xf_IPointType *IP, int M, int dim)
{
  int ierr;

  IP->M = M;

  IP->Dim = dim;

  ierr = xf_Error(xf_Alloc((void **) &IP->proc, M, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) &IP->egrp, M, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &IP->elem, M, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &IP->node, M, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) &IP->xref, M*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyIPoint
static int 
xf_DestroyIPoint(xf_IPointType *IP)
{

  xf_Release( (void *) IP->proc);
  xf_Release( (void *) IP->egrp);
  xf_Release( (void *) IP->elem);
  xf_Release( (void *) IP->node);
  xf_Release( (void *) IP->xref);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateReducedModel
int 
xf_CreateReducedModel(xf_ReducedModel **pRM)
{
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pRM, 1, sizeof(xf_ReducedModel)));
  if (ierr != xf_OK) return ierr;

  (*pRM)->N  = 0;
  (*pRM)->A  = NULL;
  (*pRM)->L  = NULL;
  (*pRM)->F  = NULL;
  (*pRM)->F0 = NULL;
  (*pRM)->nNonLinear = 0;
  (*pRM)->M = NULL;
  (*pRM)->E = NULL;
  (*pRM)->D = NULL;
  (*pRM)->z = NULL;
  (*pRM)->ResTerm = NULL;
  (*pRM)->EqnSet  = NULL;
  (*pRM)->nSnap   = 0;
  (*pRM)->UN      = NULL;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AllocReducedModelLinear
int 
xf_AllocReducedModelLinear(xf_ReducedModel *RM, int N, int nOutput)
{
  int ierr;

  RM->N = N;

  ierr = xf_Error(xf_Alloc((void **) &RM->A, N*N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) &RM->L, N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  RM->nOutput = nOutput;

  ierr = xf_Error(xf_Alloc((void **) &RM->F, nOutput*N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) &RM->F0, nOutput, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AllocReducedModelNonLinear
int 
xf_AllocReducedModelNonLinear(xf_ReducedModel *RM, int nNonLinear)
{
  int ierr, i;

  RM->nNonLinear = nNonLinear;

  ierr = xf_Error(xf_Alloc( (void **) &RM->M, nNonLinear, sizeof(int))); 
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &RM->E, nNonLinear, sizeof(real *))); 
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &RM->D, nNonLinear, sizeof(real *))); 
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &RM->z, nNonLinear, sizeof(xf_IPointType))); 
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &RM->ResTerm, nNonLinear, sizeof(xf_ResTerm))); 
  if (ierr != xf_OK) return ierr;

  for (i=0; i<nNonLinear; i++){
    RM->ResTerm[i].Type = xfe_ResTermUnknown;
    /* Initialize key-value structure */
    ierr = xf_Error(xf_InitKeyValue(&RM->ResTerm[i].KeyValue));
    if (ierr != xf_OK) return ierr;
    RM->ResTerm[i].Active = xfe_False;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyReducedModel
int 
xf_DestroyReducedModel(xf_ReducedModel *RM, enum xfe_Bool DestroyEqnSet)
{
  int i, ierr;

  if (RM == NULL) return xf_OK;

  xf_Release( (void *) RM->A);
  xf_Release( (void *) RM->L);
  xf_Release( (void *) RM->F);
  xf_Release( (void *) RM->F0);

  for (i=0; i<RM->nNonLinear; i++){ 
    xf_Release( (void *) RM->E[i]);
    xf_Release( (void *) RM->D[i]);
    ierr = xf_Error(xf_DestroyIPoint(RM->z+i));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyResTerm(RM->ResTerm+i));
    if (ierr != xf_OK) return ierr;
  }

  xf_Release( (void *) RM->M);
  xf_Release( (void *) RM->E);
  xf_Release( (void *) RM->D);
  xf_Release( (void *) RM->z);
  xf_Release( (void *) RM->ResTerm);

  if (DestroyEqnSet){
    ierr = xf_Error(xf_DestroyEqnSet(RM->EqnSet, xfe_True));
    if (ierr != xf_OK) return ierr;
  }

  xf_Release( (void *) RM->UN);

  xf_Release( (void *) RM);


  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteIPointBinary
static int 
xf_WriteIPointBinary( xf_IPointType *z, FILE *fid)
{
  int M, Dim;

  // M
  if (fwrite(&z->M, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  M = z->M;

  // Dim
  if (fwrite(&z->Dim, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  Dim = z->Dim;

  // egrp
  if (fwrite(z->egrp, sizeof(int), M, fid) != M) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // elem
  if (fwrite(z->elem, sizeof(int), M, fid) != M) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // node
  if (fwrite(z->node, sizeof(int), M, fid) != M)
    return xf_Error(xf_FILE_WRITE_ERROR);

  // xref
  if (fwrite(z->xref, sizeof(real), M*Dim, fid) != M*Dim)
    return xf_Error(xf_FILE_WRITE_ERROR);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadIPointBinary
static int 
xf_ReadIPointBinary( FILE *fid, xf_IPointType *z)
{
  int ierr, M, Dim;

  // M
  if (fread(&z->M, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  M = z->M;

  // Dim
  if (fread(&z->Dim, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  Dim = z->Dim;

  // allocate z
  ierr = xf_Error(xf_AllocIPoint( z, M, Dim));
  if (ierr != xf_OK) return ierr;

  // egrp
  if (fread(z->egrp, sizeof(int), M, fid) != M) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // elem
  if (fread(z->elem, sizeof(int), M, fid) != M) 
    return xf_Error(xf_FILE_READ_ERROR);

  // node
  if (fread(z->node, sizeof(int), M, fid) != M)
    return xf_Error(xf_FILE_READ_ERROR);

  // xref
  if (fread(z->xref, sizeof(real), M*Dim, fid) != M*Dim)
    return xf_Error(xf_FILE_READ_ERROR);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteReducedModelBinary
int 
xf_WriteReducedModelBinary( xf_ReducedModel *RM, FILE *fid)
{
  int ierr, rev, nOutput;
  int N, nNonLinear, i, M;
  
  rev = 1;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // N
  if (fwrite(&RM->N, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  N = RM->N;

  if (N <= 0) return xf_Error(xf_OUT_OF_BOUNDS);

  // nOutput
  if (fwrite(&RM->nOutput, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  nOutput = RM->nOutput;
  
  // A
  if (fwrite(RM->A, sizeof(real), N*N, fid) != N*N) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // L
  if (fwrite(RM->L, sizeof(real), N, fid) != N) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // F and F0
  if (nOutput > 0){
    if (fwrite(RM->F, sizeof(real), N*nOutput, fid) != N*nOutput) 
      return xf_Error(xf_FILE_WRITE_ERROR);
    if (fwrite(RM->F0, sizeof(real), nOutput, fid) != nOutput) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }
  

  // nNonLinear
  if (fwrite(&RM->nNonLinear, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  nNonLinear = RM->nNonLinear;

  // M
  if (nNonLinear > 0){
    if (fwrite(RM->M, sizeof(int), nNonLinear, fid) != nNonLinear) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }

  for (i=0; i<nNonLinear; i++){ // loop over nonlinear terms

    M = RM->M[i];

    if (M <= 0) return xf_Error(xf_OUT_OF_BOUNDS);
    
    // E[i]
    if (fwrite(RM->E[i], sizeof(real), N*M, fid) != N*M)
      return xf_Error(xf_FILE_WRITE_ERROR);
    
    // D[i]
    if (fwrite(RM->D[i], sizeof(real), M*N, fid) != M*N) 
      return xf_Error(xf_FILE_WRITE_ERROR);

    // z[i]
    ierr = xf_Error(xf_WriteIPointBinary(RM->z + i, fid));
    if (ierr != xf_OK) return ierr;
    
    // ResTerm[i]
    ierr = xf_Error(xf_WriteResTermBinary(RM->ResTerm + i, fid));
    if (ierr != xf_OK) return ierr;

  } // i

  // EqnSet
  ierr = xf_Error(xf_WriteEqnSetBinarySerial(RM->EqnSet, fid));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadReducedModelBinary
int 
xf_ReadReducedModelBinary( FILE *fid, xf_ReducedModel *RM)
{
  int ierr, rev, nOutput;
  int N, nNonLinear, i, M;
  
  // revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev > 1) return xf_Error(xf_FILE_READ_ERROR);

  // N
  if (fread(&RM->N, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  N = RM->N;

  if (N <= 0) return xf_Error(xf_OUT_OF_BOUNDS);

  // nOutput
  if (rev >= 1){
    if (fread(&RM->nOutput, sizeof(int), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    nOutput = RM->nOutput;
  }
  else nOutput = 0;
  
  // allocate linear portion of reduced model
  ierr = xf_Error(xf_AllocReducedModelLinear(RM, N, nOutput));
  if (ierr != xf_OK) return ierr;

  // A
  if (fread(RM->A, sizeof(real), N*N, fid) != N*N) 
    return xf_Error(xf_FILE_READ_ERROR);

  // L
  if (fread(RM->L, sizeof(real), N, fid) != N) 
    return xf_Error(xf_FILE_READ_ERROR);

  // F and F0
  if (nOutput > 0){
    if (fread(RM->F, sizeof(real), N*nOutput, fid) != N*nOutput) 
    return xf_Error(xf_FILE_READ_ERROR);
    if (fread(RM->F0, sizeof(real), nOutput, fid) != nOutput) 
    return xf_Error(xf_FILE_READ_ERROR);
  }


  // nNonLinear
  if (fread(&RM->nNonLinear, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  nNonLinear = RM->nNonLinear;

  if (nNonLinear < 0) return xf_Error(xf_OUT_OF_BOUNDS);

  // allocate nonlinear portion of reduced model
  ierr = xf_Error(xf_AllocReducedModelNonLinear(RM, nNonLinear));
  if (ierr != xf_OK) return ierr;

  // M
  if (nNonLinear > 0){
    if (fread(RM->M, sizeof(int), nNonLinear, fid) != nNonLinear) 
      return xf_Error(xf_FILE_READ_ERROR);
  }

  for (i=0; i<nNonLinear; i++){ // loop over nonlinear terms

    M = RM->M[i];

    if (M <= 0) return xf_Error(xf_OUT_OF_BOUNDS);

    // allocate E
    ierr = xf_Error(xf_Alloc((void **) RM->E+i, N*M, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // E[i]
    if (fread(RM->E[i], sizeof(real), N*M, fid) != N*M)
      return xf_Error(xf_FILE_READ_ERROR);
      
    // allocate D
    ierr = xf_Error(xf_Alloc((void **) RM->D+i, M*N, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // D[i]
    if (fread(RM->D[i], sizeof(real), M*N, fid) != M*N) 
      return xf_Error(xf_FILE_READ_ERROR);

    // z[i]
    ierr = xf_Error(xf_ReadIPointBinary(fid, RM->z + i));
    if (ierr != xf_OK) return ierr;

    // ResTerm[i]
    ierr = xf_Error(xf_ReadResTermBinary(fid, RM->ResTerm + i));
    if (ierr != xf_OK) return ierr;

  } // i

  // alloc EqnSet
  ierr = xf_Error(xf_CreateEqnSet(&(RM->EqnSet)));
  if (ierr != xf_OK) return ierr;
  
  // read EqnSet
  ierr = xf_Error(xf_ReadEqnSetBinary(fid, RM->EqnSet));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteReducedModelAscii
int 
xf_WriteReducedModelAscii( xf_ReducedModel *RM, FILE *fid)
{
  int ierr, nOutput;
  int N, nNonLinear, i, j, k, M;
  

  // N
  fprintf(fid, "N = %d; %% # basis vectors\n", RM->N);
  
  N = RM->N;

  if (N <= 0) return xf_Error(xf_OUT_OF_BOUNDS);

  // nOutput
  fprintf(fid, "nOutput = %d; %% # outputs\n", RM->nOutput);

  nOutput = RM->nOutput;
  
  // A
  fprintf(fid, "%% Reduced matrix \nA = [\n");
  for (i=0; i<N; i++){
    for (j=0; j<N; j++) fprintf(fid, "%.15E ", RM->A[i*N+j]);
    fprintf(fid, "\n");
  }
  fprintf(fid,"];\n");

  // L
  fprintf(fid, "%% Load Vector \nL = [\n");
  for (i=0; i<N; i++)
    fprintf(fid, "%.15E\n", RM->L[i]);
  fprintf(fid,"];\n");

  // F and F0
  if (nOutput > 0){
    fprintf(fid, "%% Output Matrix \nF = [\n");
    for (i=0; i<nOutput; i++){
      for (j=0; j<N; j++) fprintf(fid, "%.15E ", RM->F[i*N+j]);
      fprintf(fid, "\n");
    }
    fprintf(fid,"];\n");
    
    fprintf(fid, "%% Output Vector \nF0 = [\n");
    for (i=0; i<nOutput; i++)
      fprintf(fid, "%.15E\n", RM->F0[i]);
    fprintf(fid,"];\n");
  }
  

  // nNonLinear
  fprintf(fid, "nNonLinear = %d; %% # nonlinear terms (usually 1)\n", RM->nNonLinear);

  nNonLinear = RM->nNonLinear;

  for (i=0; i<nNonLinear; i++){ // loop over nonlinear terms

    // M
    M = RM->M[i];
    fprintf(fid, "M = %d; %% # basis functions/points for nonlinear term\n", M);

    if (M <= 0) return xf_Error(xf_OUT_OF_BOUNDS);
    
    // E[i]
    fprintf(fid, "%% E matrix \nE = [\n");
    for (j=0; j<N; j++){
      for (k=0; k<M; k++) fprintf(fid, "%.15E ", RM->E[i][j*M+k]);
      fprintf(fid, "\n");
    }
    fprintf(fid,"];\n");
    
    // D[i]
    fprintf(fid, "%% D matrix \nD = [\n");
    for (j=0; j<M; j++){
      for (k=0; k<N; k++) fprintf(fid, "%.15E ", RM->D[i][j*N+k]);
      fprintf(fid, "\n");
    }
    fprintf(fid,"];\n");
  } // i

  // nSnap
  fprintf(fid, "nSnap = %d; %% # snapshots \n",RM->nSnap);
  fprintf(fid, "%% UN matrix \nUN = [\n");
  for (i=0; i<RM->nSnap; i++){
    for (j=0; j<N; j++) fprintf(fid, "%.15E ", RM->UN[i*N+j]);
    fprintf(fid, "\n");
  }
  fprintf(fid,"];\n");


  return xf_OK;
}

