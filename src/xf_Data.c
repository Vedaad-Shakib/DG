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
 FILE:  xf_Data.c
 
 This file contains functions for working with the Data structure.
 
 */


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_EqnSet.h"
#include "xf_Basis.h"
#include "xf_MeshTools.h"
#include "xf_Math.h"
#include "xf_Solver.h"
#include "xf_Line.h"
#include "xf_LinearSolver.h"
#include "xf_DataMath.h"
#include "xfYu_Model.h"


/******************************************************************/
//   FUNCTION Definition: xf_InterpOrder
int
xf_InterpOrder( const xf_Vector *U, int egrp, int elem){
  
  return ((U->vOrder==NULL) ? U->Order[egrp] : U->vOrder[egrp][elem]);
  
}

/******************************************************************/
//   FUNCTION Definition: xf_InitGenArray
static void 
xf_InitGenArray( xf_GenArray *ga ){
  
  ga->n = 0;
  ga->r = 0;
  ga->vr = NULL;
  ga->iValue  = NULL;
  ga->rValue  = NULL;
  ga->ParallelInfo = NULL;
  
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateGenArray
static int 
xf_CreateGenArray( xf_GenArray **pGenArray ){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pGenArray, 1, sizeof(xf_GenArray)));
  if (ierr != xf_OK) return ierr;
  
  xf_InitGenArray(*pGenArray);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateArrayParallelInfo
static int 
xf_CreateArrayParallelInfo( xf_ArrayParallelInfo **pParallelInfo){
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pParallelInfo, 1, sizeof(xf_ArrayParallelInfo)));
  if (ierr != xf_OK) return ierr;
  
  (*pParallelInfo)->HaloFlag  = 0;
  (*pParallelInfo)->Request   = NULL;
  (*pParallelInfo)->nSendElem = NULL;
  (*pParallelInfo)->SendElem  = NULL;
  (*pParallelInfo)->nRecvElem = NULL;
  (*pParallelInfo)->iSendBuf  = NULL;
  (*pParallelInfo)->rSendBuf  = NULL;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyArrayParallelInfo
static int 
xf_DestroyArrayParallelInfo( xf_ArrayParallelInfo *ParallelInfo){
  
  if (ParallelInfo == NULL) return xf_OK;
  
  xf_Release( (void   *) ParallelInfo->Request);
  xf_Release2((void  **) ParallelInfo->iSendBuf);
  xf_Release2((void  **) ParallelInfo->rSendBuf);
  xf_Release( (void   *) ParallelInfo);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyGenArray
static int 
xf_DestroyGenArray( xf_GenArray *GenArray )
{
  int ierr;
  
  xf_Release(  (void *)  GenArray->vr);
  xf_Release2( (void **) GenArray->iValue);
  xf_Release2( (void **) GenArray->rValue);
  ierr = xf_Error(xf_DestroyArrayParallelInfo(GenArray->ParallelInfo));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) GenArray);
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyGenArrays
static int 
xf_DestroyGenArrays( int nArray, xf_GenArray *GenArrays )
{
  int ierr, i;
  
  for (i=0; i<nArray; i++){
    xf_Release ( (void  *) GenArrays[i].vr);
    xf_Release2( (void **) GenArrays[i].iValue);
    xf_Release2( (void **) GenArrays[i].rValue);
    ierr = xf_Error(xf_DestroyArrayParallelInfo(GenArrays[i].ParallelInfo));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release( (void *) GenArrays);
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_InitVector
void 
xf_InitVector(xf_Vector *V)
{
  
  V->Linkage = xfe_LinkageNone;
  V->LinkageIndex = 0;
  V->SolverRole = xfe_SolverRoleNone;
  V->StateRank = 0;
  V->StateName = NULL;
  V->OutputName = NULL;
  V->TimeIndex = 0;
  V->MGIndex = 0;
  V->ParallelFlag = xfe_False;
  V->HaloInTransit = xfe_False;
  V->nArray = 0;
  V->nArraySelf = 0;
  V->Order = NULL;
  V->nComp = NULL;
  V->vOrder = NULL;
  V->Basis = NULL;
  V->Size = xfe_SizeReal;
  V->GenArray = NULL;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateVector
int 
xf_CreateVector(xf_Vector **pV){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pV, 1, sizeof(xf_Vector)));
  if (ierr != xf_OK) return ierr;
  
  xf_InitVector((*pV));
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyVector
int 
xf_DestroyVector(xf_Vector *V, enum xfe_Bool DestroySelf)
{  
  int ierr;
  
  if (V == NULL) return xf_OK;
  
  xf_Release2((void **) V->StateName);
  xf_Release( (void * ) V->Order);
  xf_Release( (void * ) V->nComp);
  xf_Release2((void **) V->vOrder);
  xf_Release( (void * ) V->Basis);
  xf_Release( (void * ) V->OutputName);
  
  ierr = xf_Error(xf_DestroyGenArrays(V->nArray, V->GenArray));
  if (ierr != xf_OK) return ierr;
  
  if (DestroySelf) xf_Release((void *) V);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateMatrix
int 
xf_CreateMatrix(xf_Matrix **pM){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pM, 1, sizeof(xf_Matrix)));
  if (ierr != xf_OK) return ierr;
  
  (*pM)->Linkage = xfe_LinkageNone;
  (*pM)->LinkageIndex = 0;
  (*pM)->Order1 = -1;
  (*pM)->Order2 = -1;
  (*pM)->Basis1 = xfe_TriLagrange;
  (*pM)->Basis2 = xfe_TriLagrange;
  (*pM)->GenArray = NULL;
  (*pM)->P = NULL;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyMatrix
int 
xf_DestroyMatrix(xf_Matrix *M){
  if (M == NULL) return xf_OK;
  xf_DestroyGenArray(M->GenArray);
  xf_Release2((void **)M->P);
  xf_Release((void *) M);
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_CreateJacobianMatrix
static int 
xf_CreateJacobianMatrix( xf_JacobianMatrix **pR_U){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pR_U, 1, sizeof(xf_JacobianMatrix)));
  if (ierr != xf_OK) return ierr;
  
  (*pR_U)->Preconditioner = xfe_PreconditionerNone;
  (*pR_U)->U         = NULL;
  (*pR_U)->LineSet   = NULL;
  (*pR_U)->ILUData   = NULL;
  (*pR_U)->negrp     = 0;
  (*pR_U)->negrphalo = 0;
  (*pR_U)->StateRank = 0;
  (*pR_U)->Basis  = NULL;
  (*pR_U)->Order  = NULL;
  (*pR_U)->nvec   = NULL;
  (*pR_U)->vnvec  = NULL;
  (*pR_U)->P      = NULL;
  (*pR_U)->egrpN  = NULL;
  (*pR_U)->elemN  = NULL;
  (*pR_U)->faceN  = NULL;
  (*pR_U)->Value  = NULL;
  (*pR_U)->R_Uc   = NULL;
  (*pR_U)->ProjectionNeeded  = xfe_False;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyJacobianMatrix
static int 
xf_DestroyJacobianMatrix(xf_JacobianMatrix *R_U){
  int ierr, i, negrp, negrphalo;
  xf_JacobianMatrix *R_Uc;
  
  if (R_U == NULL) return xf_OK;
  
  ierr = xf_Error(xf_DestroyLineSet(R_U->LineSet));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyILUData(R_U->ILUData));
  if (ierr != xf_OK) return ierr;

  negrp     = R_U->negrp;
  negrphalo = R_U->negrphalo;
  
  if (R_U->P != NULL){
    for (i=0; i<negrp; i++) xf_Release2( (void **) R_U->P[i]);
    xf_Release( (void *) R_U->P);
  }
  
  xf_Release( (void  *) R_U->Basis);
  xf_Release( (void  *) R_U->Order);
  xf_Release( (void  *) R_U->nvec);
  xf_Release2((void **) R_U->vnvec);
  
  if (R_U->egrpN != NULL){
    for (i=0; i<negrphalo; i++) xf_Release2( (void **) R_U->egrpN[i]);
    xf_Release( (void *) R_U->egrpN);
  }
  
  if (R_U->elemN != NULL){
    for (i=0; i<negrphalo; i++) xf_Release2( (void **) R_U->elemN[i]);
    xf_Release( (void *) R_U->elemN);
  }
  
  if (R_U->faceN != NULL){
    for (i=0; i<negrphalo; i++) xf_Release2( (void **) R_U->faceN[i]);
    xf_Release( (void *) R_U->faceN);
  }
  
  if (R_U->Value != NULL){
    for (i=0; i<negrphalo; i++){
      if (R_U->Value[i] != NULL){
        xf_Release( (void *) R_U->Value[i][0][0]);
        xf_Release2( (void **) R_U->Value[i]);
      }
    }
    xf_Release( (void *) R_U->Value);
  }
  
  
  if ((R_Uc = R_U->R_Uc) != NULL){
    R_Uc->LineSet = NULL;
    R_Uc->ILUData = NULL;
    xf_DestroyJacobianMatrix(R_Uc);
  }
  
  xf_Release((void *) R_U);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateVectorSet
int 
xf_CreateVectorSet(int nVector, xf_VectorSet **pVS){
  
  int ierr, i;
  
  // Allocate space for self
  ierr = xf_Error(xf_Alloc((void **) pVS, 1, sizeof(xf_VectorSet)));
  if (ierr != xf_OK) return ierr;
  
  if (nVector < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // set number of vectors
  (*pVS)->nVector = nVector;
  
  // allocate space for vectors
  ierr = xf_Error(xf_Alloc( (void **) &(*pVS)->Vector, nVector, sizeof(xf_Vector)));
  if (ierr != xf_OK) return ierr;
  
  // initialize all vectors
  for (i=0; i<nVector; i++) xf_InitVector((*pVS)->Vector + i);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_TrimVectorSet
int 
xf_TrimVectorSet(xf_VectorSet **pVS, int nVectorNew)
{  
  int ierr, i;
  
  if (nVectorNew > (*pVS)->nVector) return xf_Error(xf_OUT_OF_BOUNDS);
  
  for (i=nVectorNew; i<(*pVS)->nVector; i++){
    ierr = xf_Error(xf_DestroyVector((*pVS)->Vector + i, xfe_False));
    if (ierr != xf_OK) return ierr;
  }
  
  // set number of vectors
  (*pVS)->nVector = nVectorNew;
  
  // reallocate space for vectors
  ierr = xf_Error(xf_ReAlloc( (void **) &(*pVS)->Vector, nVectorNew, sizeof(xf_Vector)));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyVectorSet
int 
xf_DestroyVectorSet(xf_VectorSet *VS)
{  
  int ierr, i;
  
  if (VS == NULL) return xf_OK;
  
  for (i=0; i<VS->nVector; i++){
    ierr = xf_Error(xf_DestroyVector(VS->Vector + i, xfe_False));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *) VS->Vector);
  xf_Release((void *) VS);
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateData
static int 
xf_CreateData( xf_Data **pData){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pData, 1, sizeof(xf_Data)));
  if (ierr != xf_OK) return ierr;
  
  (*pData)->Title     = NULL;
  (*pData)->Type      = xfe_DataOther;
  (*pData)->ReadWrite = xfe_False;
  (*pData)->Data      = NULL;
  (*pData)->Prev      = NULL;
  (*pData)->Next      = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyData
static int 
xf_DestroyData( xf_Data *Data)
{
  int ierr;
  xf_Data *Prev, *Next;
  
  if (Data == NULL) return xf_OK;
  
  switch(Data->Type){
      
    case xfe_Vector:
      ierr = xf_Error(xf_DestroyVector((xf_Vector *) Data->Data, xfe_True));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_Matrix:
      ierr = xf_Error(xf_DestroyMatrix((xf_Matrix *) Data->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_JacobianMatrix:
      ierr = xf_Error(xf_DestroyJacobianMatrix((xf_JacobianMatrix *) Data->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_VectorSet:
      ierr = xf_Error(xf_DestroyVectorSet((xf_VectorSet *) Data->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_DataOther:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    default:
      xf_printf("Error destroying Data.  Unknown type.\n");
      return xf_Error(xf_UNKNOWN_TYPE);
      break;
  }
  
  xf_Release((void *)Data->Title);
  
  Prev = Data->Prev;
  Next = Data->Next;
  
  xf_Release((void *)Data);
  
  /* Keep double linked list tied together */
  if (Prev != NULL) Prev->Next = Next;
  if (Next != NULL) Next->Prev = Prev;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateDataSet
int 
xf_CreateDataSet( xf_DataSet **pDataSet){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pDataSet, 1, sizeof(xf_DataSet)));
  if (ierr != xf_OK) return ierr;
  
  (*pDataSet)->Head = NULL;
  (*pDataSet)->Tail = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDataSet
int 
xf_DestroyDataSet( xf_DataSet *DataSet){
  
  int ierr;
  xf_Data *Data, *D;
  
  if (DataSet == NULL) return xf_OK;
  
  // Destroy linked list
  Data = DataSet->Head;
  while (Data != NULL){
    D = Data->Next;
    ierr = xf_Error(xf_DestroyData(Data));
    if (ierr != xf_OK) return ierr;
    Data = D;
  }
  
  xf_Release( (void *) DataSet);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDataInSet
int 
xf_DestroyDataInSet( xf_DataSet *DataSet, xf_Data *Data)
{
  int ierr;
  xf_Data *P, *N;
  
  P = Data->Prev;
  N = Data->Next;
  
  // make sure head/tail are linked correctly to the DataSet
  if ((P == NULL) && (DataSet->Head != Data)) return xf_Error(xf_OUT_OF_BOUNDS);
  if ((N == NULL) && (DataSet->Tail != Data)) return xf_Error(xf_OUT_OF_BOUNDS);
  
  ierr = xf_Error(xf_DestroyData(Data));
  if (ierr != xf_OK) return ierr;
  if (P == NULL) DataSet->Head = N; // keep Head valid
  if (N == NULL) DataSet->Tail = P; // keep Tail valid
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DataSetAdd
int 
xf_DataSetAdd( xf_DataSet *DataSet, const char Title[],
              const enum xfe_DataType Type, const enum xfe_Bool ReadWrite, 
              void *Data, xf_Data **pD){
  
  int ierr;
  xf_Data *D, *Tail;
  
  // allocate memory
  ierr = xf_Error(xf_CreateData(&D));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_AllocString(&D->Title, xf_MAXSTRLEN, Title));
  if (ierr != xf_OK) return ierr;
  
  D->Type      = Type;
  D->ReadWrite = ReadWrite;
  D->Data      = Data;
  
  
  if (DataSet->Tail == NULL){
    // first time adding to list
    DataSet->Head = D;
    DataSet->Tail = D;
  }
  else{
    // link to end of list
    Tail = DataSet->Tail;
    D->Prev = Tail;
    Tail->Next = D;
    DataSet->Tail = D;
  }
  
  if (pD != NULL) (*pD) = D;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DataSetRemove
int 
xf_DataSetRemove( xf_DataSet *DataSet, const char Title[], enum xfe_Bool RootFlag){
  
  int ierr, len;
  enum xfe_Bool match;
  xf_Data *Data, *N;
  
  len = strlen(Title);
  
  // Destroy linked list
  Data = DataSet->Head;
  while (Data != NULL){
    N = Data->Next;
    
    if (RootFlag)
      match = (strncmp(Data->Title, Title, len) == 0);
    else
      match = (strcmp(Data->Title, Title) == 0);
    
    if (match){
      ierr = xf_Error(xf_DestroyDataInSet(DataSet, Data));
      if (ierr != xf_OK) return ierr;
    }
    Data = N;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindDataByTitle
int 
xf_FindDataByTitle( xf_DataSet *DataSet, const char Title[], 
                   enum xfe_DataType DataType, xf_Data **pD)
{
  
  int ierr;
  enum xfe_Bool found;
  xf_Data *D;
  
  //len = strlen(Title);
  
  (*pD) = NULL;
  found = xfe_False;
  
  D = DataSet->Head;
  while (D != NULL){
    if ((strcmp(D->Title, Title) == 0) && (D->Type == DataType)){
      if (found) return xf_Error(xf_MULTIPLE_MATCHES);
      (*pD) = D;
      found = xfe_True;
    }
    D = D->Next;
  }
  if (!found) return xf_NOT_FOUND;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DataSetMerge
int 
xf_DataSetMerge( xf_DataSet *DataSetFrom, xf_DataSet *DataSetTo)
{
  
  int ierr;
  xf_Data *D, *D2;
  
  
  D = DataSetFrom->Head;
  
  while (D != NULL){
    
    xf_printf("Merging %s.\n", D->Title);
    
    // remove existing data with same name from DataSetTo 
    ierr = xf_Error(xf_DataSetRemove(DataSetTo, D->Title, xfe_False));
    if (ierr != xf_OK) return ierr;
    
    // add the data to DataSetTo
    ierr = xf_Error(xf_DataSetAdd(DataSetTo, D->Title, D->Type,
                                  D->ReadWrite, D->Data, &D2));
    if (ierr != xf_OK) return ierr;
    
    D->Data = NULL; // this is the part that is destructive to DataSetFrom
    
    D = D->Next;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AllocGenArray
int 
xf_AllocGenArray(xf_GenArray *ga){
  
  int ierr;
  int s;
  void ***target = NULL;
  
  switch (ga->Size){
    case xfe_SizeInt:
      s = sizeof(int);
      target = (void ***) &(ga->iValue);
      ga->rValue = NULL;
      break;
    case xfe_SizeReal:
      s = sizeof(real);
      target = (void ***) &(ga->rValue);
      ga->iValue = NULL;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  if (ga->vr == NULL){
    ierr = xf_Error(xf_Alloc2(target, ga->n, ga->r, s));
    if (ierr != xf_OK) return ierr;
  }
  else{
    ierr = xf_Error(xf_VAlloc2(target, ga->n, ga->vr, s));
    if (ierr != xf_OK) return ierr;
  }
  
  ga->ParallelInfo = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelPrepVector
static int 
xf_ParallelPrepVector( xf_Mesh *Mesh, xf_Vector *V)
{
  /*
   PURPOSE:
   
   Prepares V for parallel computation.  Returns immediately if not
   running in parallel.  Otherwise sets V->ParallelFlag to True and
   allocates any required send buffers for communication.
   
   INPUTS:
   
   Mesh : associated mesh, for linkage purposes
   V : vector to prepare for parallel computation
   
   OUTPUTS: 
   
   None
   
   RETURN:
   
   Error Code
   */
  int ierr, i, j, egrp;
  int myRank, nProc, iProc;
  int *nSend;
  enum xfe_Bool HaloFlag;
  xf_ArrayParallelInfo *PInfo;
  
  V->ParallelFlag = xfe_False;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // if not parallel, return immediately
  if (nProc == 1) return xf_OK;
  
  /*   if ((V->MeshParallelInfo = Mesh->ParallelInfo) == NULL) */
  /*     return xf_Error(xf_PARALLEL_ERROR); */
  
  if (V->Linkage == xfe_LinkageGlobElem){
    
    V->ParallelFlag = xfe_True;
    nSend = NULL;
    
    for (i=0; i<V->nArray; i++){
      
      HaloFlag = (i >= Mesh->nElemGroup);
      egrp = i % Mesh->nElemGroup;
      
      ierr = xf_Error(xf_CreateArrayParallelInfo(&V->GenArray[i].ParallelInfo));
      if (ierr != xf_OK) return ierr;
      PInfo = V->GenArray[i].ParallelInfo;
      
      PInfo->HaloFlag = HaloFlag;
      ierr = xf_Error(xf_Alloc( &PInfo->Request, nProc, xf_MPI_RequestSize()));
      if (ierr != xf_OK) return ierr;
      
      if (!HaloFlag){
        fflush(stdout);
        PInfo->nSendElem = Mesh->ParallelInfo->nSendElem[egrp];
        PInfo->SendElem  = Mesh->ParallelInfo->SendElem[egrp];
        
        // Allocate send buffers
        ierr = xf_Error(xf_ReAlloc( (void **) &nSend, nProc, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        for (iProc=0; iProc<nProc; iProc++){
          if (V->GenArray[i].vr == NULL)
            nSend[iProc] = PInfo->nSendElem[iProc]*V->GenArray[i].r;
          else // account for variable orders
            for (j=0, nSend[iProc]=0; j<PInfo->nSendElem[iProc]; j++)
              nSend[iProc] += V->GenArray[i].vr[PInfo->SendElem[iProc][j]];
        }
        
        switch (V->GenArray[i].Size){
          case xfe_SizeInt:
            ierr = xf_Error(xf_VAlloc2((void ***) &PInfo->iSendBuf, nProc, nSend, sizeof(int)));
            if (ierr != xf_OK) return ierr;
            break;
          case xfe_SizeReal:
            ierr = xf_Error(xf_VAlloc2((void ***) &PInfo->rSendBuf, nProc, nSend, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            break;
          default:
            return xf_Error(xf_NOT_SUPPORTED);
            break;
        }
      }
      else
        PInfo->nRecvElem = Mesh->ParallelInfo->nRecvElem[egrp];
      
    } // i
    
    xf_Release( (void *) nSend);
    
  }
  else return xf_Error(xf_NOT_SUPPORTED);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindPrimalState
int 
xf_FindPrimalState( xf_DataSet *DataSet, int TimeIndex, xf_Data **pState,
                   char *AltName)
{
  int ierr, len;
  char root[xf_MAXSTRLEN] = "State";
  xf_Data *D;
  xf_Vector *V;
  
  D = DataSet->Head;
  
  // use alternate name if specified
  if (AltName != NULL) strcpy(root, AltName);
  
  len = strlen(root);
  
  while (D != NULL){
    if ((strncmp(D->Title, root, len) == 0) &&
        (D->Type == xfe_Vector)){
      V = (xf_Vector *) D->Data;
      if ((V->TimeIndex  == TimeIndex) && 
          (V->SolverRole == xfe_SolverRolePrimalState)){
        (*pState) = D; // found matching state, return
        return xf_OK;
      }
    }
    D = D->Next;
  }
  
  return xf_NOT_FOUND;
}

/******************************************************************/
//   FUNCTION Definition: xf_CompatibleLinkageVector
static enum xfe_Bool
xf_CompatibleLinkageVector( xf_Mesh *Mesh, xf_Vector *A){
  // returns True if A has enough (or more) size for linkage to Mesh
  int i, nArray, ierr;
  enum xfe_Bool compatible;
  enum xfe_Bool ParallelFlag;
  
  compatible = xfe_True;
  
  ParallelFlag = A->ParallelFlag;
  if ((ParallelFlag) && (Mesh->ParallelInfo == NULL)) compatible = xfe_False;  
  if (A->GenArray == NULL) compatible = xfe_False;
  
  switch (A->Linkage){
      
    case xfe_LinkageGlobElem:
      nArray = ((Mesh->ParallelInfo == NULL) ? Mesh->nElemGroup : 2*Mesh->nElemGroup);
      if (A->nArray != nArray) compatible = xfe_False; 
      for (i=0; (i<A->nArray)&&(compatible); i++){
        if (A->GenArray[i].n < Mesh->ElemGroup[i].nElem) compatible = xfe_False;
      }
      break;
      
    case xfe_LinkageFace:
      nArray = Mesh->nBFaceGroup+1;
      if (A->nArray != nArray) compatible = xfe_False; 
      if (A->GenArray[0].n < Mesh->nIFace) compatible = xfe_False;
      for (i=1; (i<A->nArray)&&(compatible); i++){
        if (A->GenArray[i].n < Mesh->BFaceGroup[i-1].nBFace) compatible = xfe_False;
      }
      break;
      
    default:
      compatible = xfe_False;
      break;
  }
  
  return compatible; 
}


/******************************************************************/
//   FUNCTION Definition: xf_CompatibleSizeVector
static enum xfe_Bool
xf_CompatibleSizeVector( xf_Vector *A, int negrp, int StateRank,
                        enum xfe_BasisType *Basis, int *Order, 
                        int *nComp, int **vOrder, int *rvec,
                        enum xfe_SizeType Size, 
                        enum xfe_Bool ParallelFlag){
  // returns True if A matches requested size
  int ierr, i, j, r;
  enum xfe_Bool compatible;
  enum xfe_Bool Interpolated, VariableOrder;
  xf_GenArray *gA;
  
  Interpolated  = ((Basis != NULL) && ( Order != NULL));
  VariableOrder = ((nComp != NULL) && (vOrder != NULL));
  
  compatible = xfe_True;
  
  // check size of .r (max r in case of variable order)
  if ((Interpolated) && ((A->Basis == NULL) || (A->Order == NULL)))
    compatible = xfe_False;
  else{
    
    // check size of vectors for variable order
    if (VariableOrder){
      if ((A->nComp == NULL) || (A->vOrder == NULL)) compatible = xfe_False;
      else{
        for (i=0; (i<A->nArray)&&(compatible); i++){
          if (A->nComp[i] != nComp[i]) compatible = xfe_False;
          else{
            for (j=0; (j<A->nComp[i])&&(compatible); j++)
              if (A->vOrder[i][j] != vOrder[i][j]) compatible = xfe_False;
          }
        } // i
      }
    }
    
    for (i=0; (i<A->nArray)&&(compatible); i++){ 
      gA = A->GenArray+i;
      if (Interpolated){
        if (VariableOrder){
          if (A->nComp[i] != nComp[i]) compatible = xfe_False;
          if (gA->n  != nComp[i]) compatible = xfe_False;
          if ((gA->n>0) && (gA->vr == NULL)) compatible = xfe_False;
          
        }
        else{
          if ((A->Basis[i] != Basis[i%negrp]) ||
              (A->Order[i] != Order[i%negrp])){
            compatible = xfe_False;
          }
        }
      }
      else{ // if not interpolated, r is rvec[i] or StateRank
        r = ((rvec == NULL) ? StateRank : rvec[i%negrp]);
        if (gA->r != r) compatible = xfe_False;
      }
    } // i
  }
  
  if (A->Size != Size) compatible = xfe_False;
  if (A->ParallelFlag != ParallelFlag) compatible = xfe_False;
  
  return compatible; 
}


/******************************************************************/
//   FUNCTION Definition: xf_CompatibleVectors
enum xfe_Bool
xf_CompatibleVectors( const xf_Vector *A, const xf_Vector *B){
  int ierr, i, j;
  enum xfe_Bool compatible;
  enum xfe_Bool VariableOrderA, VariableOrderB;
  xf_GenArray *gA, *gB;
  
  if (A->nArray != B->nArray) {
    //xf_printf("nArray: %d vs %d\n", A->nArray, B->nArray);
    return xfe_False;
  }
  if ((A->GenArray == NULL) || 
      (B->GenArray == NULL)){
    //xf_printf("GenArray mismatch\n");
    return xfe_False;
  }
  
  VariableOrderA = ((A->nComp != NULL) && (A->vOrder != NULL));
  VariableOrderB = ((B->nComp != NULL) && (B->vOrder != NULL));
  if (VariableOrderA != VariableOrderB){
    //xf_printf("Variable order mismatch\n");
    return xfe_False;
  }
  
  compatible = xfe_True;
  
  for (i = 0; i<A->nArray; i++){
    gA = A->GenArray+i;
    gB = B->GenArray+i;
    if ((gA->n==0) && (gB->n==0)) continue;
    if ((gA->Size != gB->Size) || 
        (gA->n != gB->n) || 
        (gA->r != gB->r)){
      xf_pprintf("%d, %d: [%d, %d,%d] vs [%d, %d,%d]\n", i, A->nArray, gA->Size,
                 gA->n, gA->r, gB->Size, gB->n, gB->r);
      compatible = xfe_False; break;
    }
    if (VariableOrderA){
      if ((A->nComp[i]==0) && (B->nComp[i]==0)) continue;
      if ((A->nComp[i] != B->nComp[i]) ||
          (gA->vr == NULL) || (gB->vr == NULL)){
	//xf_printf("A->nComp[%d]=%d, B->nComp[%d]=%d, ga->vr==NULL (%d), gB->vr==NULL (%d)\n",
	//	  i, A->nComp[i], i, B->nComp[i], gA->vr==NULL, gB->vr==NULL);
        compatible = xfe_False;
      }
      for (j=0; (j<A->nComp[i])&&(compatible); j++)
        if (gA->vr[j] != gB->vr[j]){
	  //xf_printf("gA->vr[%d]=%d, gB->vr[%d]=%d\n", j, gA->vr[j], j, gB->vr[j]);
          compatible = xfe_False;
	}
      if (!compatible) break;
    }
  }
  
  return compatible; 
}



/******************************************************************/
//   FUNCTION Definition: xf_CopyGenArray
static int
xf_CopyGenArray(xf_GenArray *gA, xf_GenArray *gB)
{
  // gB must already be allocated, but not the data in it
  int ierr, i, j;
  
  // set scalars and initialize pointers
  gB->Size = gA->Size;
  gB->n = gA->n;
  gB->r = gA->r;
  gB->vr = NULL;
  gB->iValue = NULL;
  gB->rValue = NULL;
  
  // account for variable order (vector of r values)
  if (gA->vr != NULL){
    ierr = xf_Error(xf_Alloc( (void **) &gB->vr, gB->n, sizeof(int))); 
    if (ierr != xf_OK) return ierr;
    for (i=0; i<gB->n; i++) gB->vr[i] = gA->vr[i];
  }
  
  // Allocate data in gB
  ierr = xf_Error(xf_AllocGenArray(gB));
  if (ierr != xf_OK) return ierr;
  
  if (gA->iValue != NULL){
    for (i=0; i<gB->n; i++) 
      for (j=0; j<((gB->vr==NULL) ? gB->r : gB->vr[i]); j++)   
        gB->iValue[i][j] = gA->iValue[i][j];
  }
  
  if (gA->rValue != NULL){
    for (i=0; i<gB->n; i++)
      for (j=0; j<((gB->vr==NULL) ? gB->r : gB->vr[i]); j++)
        gB->rValue[i][j] = gA->rValue[i][j];
  }
  
  gB->ParallelInfo = NULL;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CopyVector
int
xf_CopyVector(xf_Mesh *Mesh, xf_Vector *U, xf_Vector *V)
{
  int ierr, i, j;
  
  if (U == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // set scalars
  V->Linkage = U->Linkage;
  V->LinkageIndex = U->LinkageIndex;
  V->SolverRole = U->SolverRole;
  V->StateRank = U->StateRank;
  V->TimeIndex = U->TimeIndex;
  V->MGIndex = U->MGIndex;
  V->ParallelFlag = U->ParallelFlag;
  V->HaloInTransit = U->HaloInTransit;
  V->nArray = U->nArray;
  V->nArraySelf = U->nArraySelf;
  V->Size = U->Size;
  
  // StateName
  if (U->StateName != NULL){
    ierr = xf_Error(xf_Alloc2((void ***) &V->StateName, U->StateRank, 
                              xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<U->StateRank; i++)
      strncpy(V->StateName[i], U->StateName[i], xf_MAXSTRLEN);
  }
  
  // OutputName
  if (U->OutputName != NULL){
    ierr = xf_Error(xf_Alloc((void **) &V->OutputName, xf_MAXSTRLEN,
                             sizeof(char)));
    if (ierr != xf_OK) return ierr;
    strncpy(V->OutputName, U->OutputName, xf_MAXSTRLEN);
  }
  
  // Basis
  if (U->Basis != NULL){
    ierr = xf_Error(xf_Alloc((void **) &V->Basis, V->nArray, 
                             sizeof(enum xfe_BasisType)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<V->nArray; i++) V->Basis[i] = U->Basis[i];
  }
  
  // Order
  if (U->Order != NULL){
    ierr = xf_Error(xf_Alloc((void **) &V->Order, V->nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<V->nArray; i++) V->Order[i] = U->Order[i];
  }
  
  // nComp and vOrder
  if (U->nComp != NULL){
    ierr = xf_Error(xf_Alloc((void **) &V->nComp, V->nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<V->nArray; i++) V->nComp[i] = U->nComp[i];
    
    // vOrder
    if (U->vOrder != NULL){
      ierr = xf_Error(xf_VAlloc2((void ***) &V->vOrder, V->nArray, V->nComp, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<V->nArray; i++) 
        for (j=0; j<V->nComp[i]; j++)
          V->vOrder[i][j] = U->vOrder[i][j];
    }
  }
  
  // GenArray
  if (U->GenArray != NULL){
    ierr = xf_Error(xf_Alloc((void **) &V->GenArray, V->nArray, 
                             sizeof(xf_GenArray)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<V->nArray; i++){
      ierr = xf_Error(xf_CopyGenArray(U->GenArray+i, V->GenArray+i));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // prepare V for parallel communication if necessary
  if ((U->ParallelFlag) && (Mesh != NULL)){
    ierr = xf_Error(xf_ParallelPrepVector(Mesh, V));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DuplicateVector
int
xf_DuplicateVector(xf_Mesh *Mesh, xf_Vector *U, xf_Vector **pV)
{
  int ierr;
  
  ierr = xf_Error(xf_CreateVector(pV));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_CopyVector(Mesh, U, (*pV));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AllocVectorGuts
static int 
xf_AllocVectorGuts(xf_All *All, enum xfe_LinkageType Linkage, int StateRank,
                   char **StateName, int TimeIndex, int MGIndex, 
                   enum xfe_BasisType *Basis, int *Order, int *nComp, int **vOrder,
                   int *rvec, enum xfe_SizeType Size, enum xfe_Bool ParallelFlag, 
                   xf_Vector *V)
{
  /*
   PURPOSE:
   
   Allocates vector data, including general arrays
   
   INPUTS:
   
   All : All structure 
   Title: Title of B to search for
   Linkage, StateRank, TimeIndex, MGIndex, Size: matching characteristics
   StateName : array of strings containing names of state components (optional)
   Basis, Order : if both not NULL, these specify the desired inteprolation
   type for each array.
   nComp, vOrder : number of components (e.g. elements) and variable order vector
   if doing variable order interpolation
   rvec : vector of ranks for each array, one for each array.
   Optional, and only used if not interpolated.  If NULL passed
   in and not interpolated, StateRank is used as the rank for
   all arrays.
   Size : integer or real
   ParallelFlag: if true, means that the created vector (*pB) will be prepped for
   parallel comm, assuming that Mesh is parallelized
   V : vector for which "guts" are allocated 
   
   OUTPUTS:  None
   
   RETURN:
   
   Error Code
   */
  
  int ierr, negrp, i, j, r;
  int ngroup, n;
  enum xfe_Bool Interpolated, VariableOrder;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  V->Linkage   = Linkage;
  V->StateRank = StateRank;
  V->TimeIndex = TimeIndex;
  V->MGIndex   = MGIndex;
  V->Size      = Size;
  
  Interpolated  = ((Basis != NULL) && ( Order != NULL));
  VariableOrder = ((nComp != NULL) && (vOrder != NULL));
  
  if (StateName != NULL){
    ierr = xf_Error(xf_Alloc2((void ***)&V->StateName, StateRank, 
                              xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<StateRank; i++) strcpy(V->StateName[i], StateName[i]);
  }
  
  switch (Linkage){
      
    case xfe_LinkageGlobElem:
      /* Create a separate array for each element group in the Mesh.  If
       Mesh is parallel, create arrays on halo groups regardless of
       ParallelFlag.  ParallelFlag=True only refers to preparing the
       vector for parallel communication (allocating buffers,
       etc.). */
      negrp = Mesh->nElemGroup;
      V->nArraySelf = negrp;
      V->nArray     = ((Mesh->ParallelInfo == NULL) ? negrp : 2*negrp);
      
      if (Interpolated){
        ierr = xf_Error(xf_Alloc((void **) &V->Basis, V->nArray, sizeof(enum xfe_BasisType)));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_Alloc((void **) &V->Order, V->nArray, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        // set Basis and Order
        for (i=0; i<V->nArray; i++){
          V->Basis[i] = Basis[i]; // used to be i%negrp
          V->Order[i] = Order[i]; // used to be i%negrp
        }
        
        // set variable order if desired
        if (VariableOrder){
          ierr = xf_Error(xf_Alloc((void **) &V->nComp, V->nArray, sizeof(int)));
          if (ierr != xf_OK) return ierr;
          for (i=0; i<V->nArray; i++) V->nComp[i] = nComp[i];
          ierr = xf_Error(xf_VAlloc2((void ***) &V->vOrder, V->nArray, V->nComp, sizeof(int)));
          if (ierr != xf_OK) return ierr;
          for (i=0; i<V->nArray; i++) 
            for (j=0; j<V->nComp[i]; j++)
              V->vOrder[i][j] = vOrder[i][j];
        }
      }

      // allocate GenArray structures
      ierr = xf_Error(xf_Alloc((void **) &V->GenArray, V->nArray, sizeof(xf_GenArray)));
      if (ierr != xf_OK) return ierr;
      
      for (i=0; i<V->nArray; i++){
        
        V->GenArray[i].Size = Size;
        V->GenArray[i].n    = Mesh->ElemGroup[i].nElem;
        V->GenArray[i].r    = 0;
        V->GenArray[i].vr   = NULL;
        
        if (Interpolated){
          ierr = xf_Error(xf_Order2nNode(V->Basis[i], V->Order[i], &V->GenArray[i].r));
          if (ierr != xf_OK) return ierr;
          V->GenArray[i].r *= StateRank;
          if (VariableOrder){  // account for variable order (vector of r values)
            ierr = xf_Error(xf_Alloc( (void **) &V->GenArray[i].vr, V->nComp[i], sizeof(int))); 
            if (ierr != xf_OK) return ierr;
            for (j=0; j<V->nComp[i]; j++){
              ierr = xf_Error(xf_Order2nNode(V->Basis[i], V->vOrder[i][j], V->GenArray[i].vr+j));
              if (ierr != xf_OK) return ierr;
              V->GenArray[i].vr[j] *= StateRank;
            }
          }
        }
        else
          V->GenArray[i].r = ((rvec == NULL) ? StateRank : rvec[i%negrp]);

        ierr = xf_Error(xf_AllocGenArray(V->GenArray+i));
        if (ierr != xf_OK) return ierr;
      }
      
      break;
      
    case xfe_LinkageFace:
      /* Create a separate array for the iface and each bface group in
       the Mesh.  No special creation if Mesh is parallel. */
      ngroup = Mesh->nBFaceGroup+1;
      V->nArraySelf = V->nArray = ngroup;
      
      // interpolated face arrays are currently not supported
      if (Interpolated) return xf_Error(xf_NOT_SUPPORTED);
      
      // allocate GenArray structures
      ierr = xf_Error(xf_Alloc((void **) &V->GenArray, V->nArray, sizeof(xf_GenArray)));
      if (ierr != xf_OK) return ierr;
      
      for (i=0; i<V->nArray; i++){
        r = ((rvec == NULL) ? StateRank : rvec[i]);
        
        n = ((i == 0) ? Mesh->nIFace : Mesh->BFaceGroup[i-1].nBFace);
        
        V->GenArray[i].Size = Size;
        V->GenArray[i].n    = n;
        V->GenArray[i].r    = r;
        V->GenArray[i].vr   = NULL;
        
        ierr = xf_Error(xf_AllocGenArray(V->GenArray+i));
        if (ierr != xf_OK) return ierr;
      }
      
      break;
      
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  // Prepare vector for parallel computation
  if (ParallelFlag){
    ierr = xf_Error(xf_ParallelPrepVector(Mesh, V));
    if (ierr != xf_OK) return ierr;    
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_NewGREVector
int 
xf_NewGREVector( xf_All *All, enum xfe_BasisType InBasis, 
                int InOrder, int StateRank, enum xfe_Bool ParallelFlag,
                xf_Vector **pV){
  
  int i, ierr, negrp;
  enum xfe_BasisType *Basis = NULL;
  enum xfe_Bool MeshIsParallel;
  int *Order = NULL;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  
  MeshIsParallel = (Mesh->ParallelInfo != NULL);
  
  // allocate memory for the V structure (not the data yet)
  ierr = xf_Error(xf_CreateVector(pV));
  if (ierr != xf_OK) return ierr;
  
  if (InOrder >= 0){
    negrp = Mesh->nElemGroup;
    if (MeshIsParallel) negrp *=2; // to account for halos
    ierr = xf_Error(xf_Alloc( (void **) &Basis, negrp, sizeof(enum xfe_BasisType)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<negrp; i++) Basis[i] = InBasis;
    
    ierr = xf_Error(xf_Alloc( (void **) &Order, negrp, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<negrp; i++) Order[i] = InOrder; 
  }
  
  ierr = xf_Error(xf_AllocVectorGuts(All, xfe_LinkageGlobElem, StateRank, NULL,
                                     0, 0, Basis, Order, NULL, NULL, NULL, 
                                     xfe_SizeReal, ParallelFlag, (*pV)));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) Basis);
  xf_Release( (void *) Order);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindVector
int 
xf_FindVector( xf_All *All, const char Title[],
              enum xfe_LinkageType Linkage, int StateRank, 
              char **StateName, int TimeIndex, int MGIndex, 
              enum xfe_BasisType *Basis, int *Order, int *nComp, int **vOrder,
              int *rvec, enum xfe_SizeType Size, enum xfe_Bool ParallelFlag,
              enum xfe_Bool AddToDataSet, xf_Data **pBData, xf_Vector **pB,
              enum xfe_Bool *Found){
  
  int ierr, len, i, j, r, negrp;
  enum xfe_Bool match, Interpolated, VariableOrder;
  xf_Data *D, *N;
  xf_Vector *B;
  xf_GenArray *gB;
  xf_DataSet *DataSet;
  xf_Mesh *Mesh;
  
  // pull Mesh and DataSet from All
  Mesh    = All->Mesh;
  DataSet = All->DataSet;
  negrp = Mesh->nElemGroup;
  
  // parallelization request valid only if Mesh is parallel
  ParallelFlag = ((ParallelFlag) && (Mesh->ParallelInfo != NULL));
  
  Interpolated  = ((Basis != NULL) && ( Order != NULL));
  VariableOrder = ((nComp != NULL) && (vOrder != NULL));
  
  len = strlen(Title);
  
  if (AddToDataSet){
    D = DataSet->Head;
    while (D != NULL){
      N = D->Next;
      if ((strlen(D->Title)==len) &&
          (strncmp(D->Title, Title, len) == 0) &&
          (D->Type == xfe_Vector)){
        B = (xf_Vector *) D->Data;
        if ((B != NULL) &&
            (B->Linkage   == Linkage) && 
            (B->StateRank == StateRank)&& 
            (B->TimeIndex == TimeIndex) &&
            (B->MGIndex   == MGIndex)){
          
          // check if B is compatible (adequate size) with Mesh
          match = xf_CompatibleLinkageVector(Mesh, B);
          
          // check if B is compatible with requested size/interpolation
          match = ((match) && (xf_CompatibleSizeVector(B, negrp, StateRank, Basis, Order, 
                                                       nComp, vOrder, rvec, Size, ParallelFlag)));
          
          if (match){ // this vector is the one we're looking for
            if (Found  != NULL) (*Found) = xfe_True;
            if (pBData != NULL) (*pBData) = D;
            if (pB     != NULL) (*pB) = B;
            return xf_OK; // done
          }
          
        }
        // new: always destroy data that had the same name but did not match
        ierr = xf_Error(xf_DestroyDataInSet(DataSet, D));
        if (ierr != xf_OK) return ierr;
      }
      D = N;
    }
  }
  
  if (Found != NULL) (*Found) = xfe_False;
  
  // at this point, no suitable B was found, so create one
  ierr = xf_Error(xf_CreateVector(&B));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_AllocVectorGuts(All, Linkage, StateRank, StateName, TimeIndex, 
                                     MGIndex, Basis, Order, nComp, vOrder,
                                     rvec, Size, ParallelFlag, B));
  if (ierr != xf_OK) return ierr;
  
  if (pB != NULL) (*pB) = B;
  
  if (AddToDataSet){ // store vector in a new data structure with Title
    ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_Vector,
                                  xfe_False, (void *) B, pBData));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindSimilarVector
int 
xf_FindSimilarVector( xf_All *All, const xf_Vector *A, const char Title[],
                     enum xfe_Bool ParallelFlag, enum xfe_Bool AddToDataSet, 
                     xf_Data **pBData, xf_Vector **pB, enum xfe_Bool *Found){
  
  int ierr, i, *rvec;
  
  // allocate and fill rvec = vector of array ranks
  ierr = xf_Error(xf_Alloc((void **) &rvec, A->nArray, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<A->nArray; i++) rvec[i] = A->GenArray[i].r;
  
  // Call FindVector with parameters taken from A
  ierr = xf_Error(xf_FindVector(All, Title, A->Linkage, A->StateRank, A->StateName,
                                A->TimeIndex, A->MGIndex, A->Basis, A->Order, 
                                A->nComp, A->vOrder, rvec, A->Size, ParallelFlag, 
                                AddToDataSet, pBData, pB, Found));
  if (ierr != xf_OK) return ierr;
  
  // release rvec
  xf_Release( (void *) rvec);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindSimilarVectors
int 
xf_FindSimilarVectors( xf_All *All, const xf_Vector *A, const char Title[],
                      int nVector, enum xfe_Bool ParallelFlag, 
                      enum xfe_Bool AddToDataSet, xf_Vector ***pB)
{
  
  int ierr, i;
  char Name[xf_MAXSTRLEN];
  
  ierr = xf_Error(xf_Alloc( (void **) pB, nVector, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < nVector; i++){
    sprintf(Name, "%s_%d", Title, i);
    ierr = xf_Error(xf_FindSimilarVector(All, A, Name, ParallelFlag, AddToDataSet,
                                         NULL, (*pB) + i, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FindSimilarVectorHO
int 
xf_FindSimilarVectorHO( xf_All *All, xf_Vector *A, const char Title[],
                       enum xfe_Bool ParallelFlag, enum xfe_Bool AddToDataSet, 
                       int *OrderIncrement, int *DesiredOrder, xf_Data **pBData, 
                       xf_Vector **pB)
{  
  int ierr, i, j, *Order = NULL, **vOrder = NULL;
  
  // vector must be interpolated
  if ((A->Basis == NULL) && (A->Order == NULL)){
    xf_printf("FindSimilarVectorHO should only be called for interpolated vectors.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  // allocate and fill Order = vector of desired orders
  if (DesiredOrder != NULL)
    Order = DesiredOrder;
  else{
    if (OrderIncrement == NULL) return xf_Error(xf_INPUT_ERROR);
    ierr = xf_Error(xf_Alloc((void **) &Order, A->nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<A->nArray; i++) Order[i] = A->Order[i] + (*OrderIncrement);
    // account for variable order
    if ((A->nComp != NULL) && (A->vOrder != NULL)){
      ierr = xf_Error(xf_VAlloc2((void ***) &vOrder, A->nArray, A->nComp, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<A->nArray; i++) 
        for (j=0; j<A->nComp[i]; j++)
          vOrder[i][j] = A->vOrder[i][j] + (*OrderIncrement);
    }
  }
  
  // Call FindVector with parameters taken from A
  ierr = xf_Error(xf_FindVector(All, Title, A->Linkage, A->StateRank, A->StateName,
                                A->TimeIndex, A->MGIndex, A->Basis, Order, A->nComp, 
                                vOrder, NULL, A->Size, ParallelFlag, AddToDataSet, 
                                pBData, pB, NULL));
  if (ierr != xf_OK) return ierr;
  
  // release Order
  if (DesiredOrder == NULL){
    xf_Release( (void * ) Order);
    xf_Release2((void **) vOrder);
  }
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_FindVectorSet
int 
xf_FindVectorSet( xf_All *All, int nVector, const char Title[],
                 enum xfe_LinkageType Linkage, int StateRank, 
                 char **StateName, int TimeIndex, int MGIndex, 
                 enum xfe_BasisType *Basis, int *Order, 
                 int *nComp, int **vOrder, int *rvec,
                 enum xfe_SizeType Size, enum xfe_Bool ParallelFlag,
                 enum xfe_Bool AddToDataSet, xf_Data **pVSData, 
                 xf_VectorSet **pVS, enum xfe_Bool *Found)
{
  int ierr, i, len, negrp;
  enum xfe_Bool match;
  xf_Data *D, *N;
  xf_VectorSet *VS;
  xf_DataSet *DataSet;
  xf_Mesh *Mesh;
  xf_Vector *B;
  
  // pull off Mesh and DataSet
  Mesh = All->Mesh;
  DataSet = All->DataSet;
  negrp = Mesh->nElemGroup;
  
  // need at least one vector requested
  if (nVector < 1) return xf_Error(xf_INPUT_ERROR);
  
  // parallelization request valid only if Mesh is parallel
  ParallelFlag = ((ParallelFlag) && (Mesh->ParallelInfo != NULL));
  
  len = strlen(Title);
  
  if (AddToDataSet){
    D = DataSet->Head;
    while (D != NULL){
      N = D->Next;
      if ((strcmp(D->Title, Title) == 0) &
          (D->Type == xfe_VectorSet)){
        VS = (xf_VectorSet *) D->Data;
        
        if (VS->nVector < 1) return xf_Error(xf_DATA_ERROR);
        B = VS->Vector + 0;
        if ((B->Linkage      == Linkage  ) && 
            (B->StateRank    == StateRank) && 
            (B->TimeIndex    == TimeIndex) &&
            (B->MGIndex      == MGIndex  )){
          
          // check if B is compatible (adequate size) with Mesh
          match = xf_CompatibleLinkageVector(Mesh, B);
          
          // check if B is compatible with requested size/interpolation
          match = ((match) && (xf_CompatibleSizeVector(B, negrp, StateRank, Basis, Order, 
                                                       nComp, vOrder, rvec, Size, ParallelFlag)));
          
          // also need to match the number of vectors
          match = ((match) && (VS->nVector >= nVector));
          
          if (match){ // All good.  Found what we're looking for.
            if (Found   != NULL) (*Found) = xfe_True;
            if (pVSData != NULL) (*pVSData) = D;
            if (pVS     != NULL) (*pVS) = VS;
            return xf_OK; // done
          }
          else{ // VS is very similar but incompatible.  Destroy it.
            ierr = xf_Error(xf_DestroyDataInSet(DataSet, D));
            if (ierr != xf_OK) return ierr;
          }
        }
      }
      D = N;
    }
  }
  
  if (Found != NULL) (*Found) = xfe_False;
  
  // at this point, no suitable VectorSet was found, so create one
  ierr = xf_Error(xf_CreateVectorSet(nVector, &VS));
  if (ierr != xf_OK) return ierr;
  
  
  // allocate each vector in the set
  for (i=0; i<nVector; i++){
    ierr = xf_Error(xf_AllocVectorGuts(All, Linkage, StateRank, StateName, TimeIndex, 
                                       MGIndex, Basis, Order, nComp, vOrder, 
                                       rvec, Size, ParallelFlag, VS->Vector+i));
    if (ierr != xf_OK) return ierr;
  }
  
  if (pVS  != NULL) (*pVS) = VS;
  
  
  if (AddToDataSet){ // store VectorSet in a new data structure with Title
    ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_VectorSet,
                                  xfe_False, (void *) VS, pVSData));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_FindSimilarVectorSet
int 
xf_FindSimilarVectorSet( xf_All *All, xf_Vector *A, int nVector, 
                        const char Title[], enum xfe_Bool ParallelFlag, 
                        enum xfe_Bool AddToDataSet, xf_Data **pVSData, 
                        xf_VectorSet **pVS)
{
  int ierr, i, *rvec;
  
  // need at least one vector requested
  if (nVector < 1) return xf_Error(xf_INPUT_ERROR);
  
  // allocate and fill rvec = vector of array ranks
  ierr = xf_Error(xf_Alloc((void **) &rvec, A->nArray, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<A->nArray; i++) rvec[i] = A->GenArray[i].r;
  
  // Call FindVectorSet with parameters taken from A
  ierr = xf_Error(xf_FindVectorSet(All, nVector, Title, A->Linkage, A->StateRank,
                                   A->StateName, A->TimeIndex, A->MGIndex, A->Basis, 
                                   A->Order, A->nComp, A->vOrder, rvec, A->Size, 
                                   ParallelFlag, AddToDataSet, pVSData, pVS, NULL));
  if (ierr != xf_OK) return ierr;
  
  // release rvec
  xf_Release( (void *) rvec);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FindSimilarVectorSetHO
int 
xf_FindSimilarVectorSetHO( xf_All *All, xf_Vector *A, int nVector, 
                          const char Title[], enum xfe_Bool ParallelFlag, 
                          enum xfe_Bool AddToDataSet, int OrderIncrement,
                          xf_Data **pVSData, xf_VectorSet **pVS)
{
  int ierr, i, j, *Order = NULL, **vOrder = NULL;
  
  // need at least one vector requested
  if (nVector < 1) return xf_Error(xf_INPUT_ERROR);
  
  // vector must be interpolated
  if ((A->Basis == NULL) && (A->Order == NULL)){
    xf_printf("FindSimilarVectorHO should only be called for interpolated vectors.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  // allocate and fill Order = vector of desired orders
  ierr = xf_Error(xf_Alloc((void **) &Order, A->nArray, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<A->nArray; i++) Order[i] = A->Order[i] + OrderIncrement;
  if ((A->nComp != NULL) && (A->vOrder != NULL)){
    ierr = xf_Error(xf_VAlloc2((void ***) &vOrder, A->nArray, A->nComp, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<A->nArray; i++) 
      for (j=0; j<A->nComp[i]; j++)
        vOrder[i][j] = A->vOrder[i][j] + OrderIncrement;
  }
  
  
  // Call FindVectorSet with parameters taken from A
  ierr = xf_Error(xf_FindVectorSet(All, nVector, Title, A->Linkage, A->StateRank, 
                                   A->StateName, A->TimeIndex, A->MGIndex, A->Basis, 
                                   Order, A->nComp, A->vOrder, NULL, A->Size, 
                                   ParallelFlag, AddToDataSet, pVSData, pVS, NULL));
  if (ierr != xf_OK) return ierr;
  
  // release rvec
  xf_Release( (void * ) Order);
  xf_Release2((void **) vOrder);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindMatrix
int 
xf_FindMatrix( xf_DataSet *DataSet, const char Title[],
              enum xfe_LinkageType Linkage, int LinkageIndex, 
              int Order1, int Order2, enum xfe_BasisType Basis1,
              enum xfe_BasisType Basis2, enum xfe_SizeType Size, 
              int n, int r, enum xfe_Bool AddToDataSet, 
              xf_Data **pMData, xf_Matrix **pM, enum xfe_Bool *Found)
{  
  int ierr, len;
  enum xfe_Bool match;
  xf_Data *D, *N;
  xf_Matrix *M;
  xf_GenArray *gM;
  
  if (Found != NULL) (*Found) = xfe_True;
  
  len = strlen(Title);
  
  if (AddToDataSet){
    D = DataSet->Head;
    while (D != NULL){
      N = D->Next;
      if ((strcmp(D->Title, Title) == 0) &
          (D->Type == xfe_Matrix)){
        M = (xf_Matrix *) D->Data;
        if ((M->Linkage  == Linkage) && 
            (M->LinkageIndex == LinkageIndex) &&
            (M->Order1 == Order1) &&
            (M->Order2 == Order2) &&
            (M->Basis1 == Basis1) &&
            (M->Basis2 == Basis2)){
          match = xfe_True;
          if ((gM = M->GenArray) != NULL){
            if ((gM->n != n) || (gM->r != r))
              match = xfe_False;
          }
          if (match){ // data ranks are the same
            if (pMData != NULL) (*pMData) = D;
            if (pM     != NULL) (*pM) = M;
            return xf_OK; // done
          }
          else{ // destroy if incorrect rank
            ierr = xf_Error(xf_DestroyDataInSet(DataSet, D));
            if (ierr != xf_OK) return ierr;
          }
        }
      }
      D = N;
    }
  }
  
  if (Found != NULL) (*Found) = xfe_False;
  
  // allocate memory for the Matrix structure (not the data yet)
  ierr = xf_Error(xf_CreateMatrix(&M));
  if (ierr != xf_OK) return ierr;
  
  M->Linkage      = Linkage;
  M->LinkageIndex = LinkageIndex;
  M->Order1       = Order1;
  M->Order2       = Order2;
  M->Basis1       = Basis1;
  M->Basis2       = Basis2;
  
  // allocate GenArray structure
  ierr = xf_Error(xf_Alloc((void **) &M->GenArray, 1, sizeof(xf_GenArray)));
  if (ierr != xf_OK) return ierr;
  
  M->GenArray[0].Size = Size;
  M->GenArray[0].n    = n;
  M->GenArray[0].r    = r;
  M->GenArray[0].vr   = NULL;
  
  ierr = xf_Error(xf_AllocGenArray(M->GenArray));
  if (ierr != xf_OK) return ierr;
  
  if (pM != NULL) (*pM) = M;
  
  if (AddToDataSet){ // store Matrix in a new data structure with Title
    ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_Matrix,
                                  xfe_False, (void *) M, pMData));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_RemapState
static int 
xf_RemapState(xf_Vector *V, xf_Vector *U)
{
  /*
   PURPOSE:
   
   Maps state components from V to U, according to the StateNames in
   both vectors.  Not all states are required to match, and the
   StateRanks of U and V need not be the same.  Only the StateNames
   that appear in both U and V are copied from V to U.
   
   INPUTS:
   
   V: vector from which to map
   U: vector which receives new values
   
   OUTPUTS: 
   
   U: modified with new values for some/all of its states
   
   RETURN:
   
   Error Code
   */
  
  int ierr, k, len;
  int egrp, elem, i, j, r, nn;
  int *P;
  enum xfe_Bool VariableOrder;
  char *s;
  xf_GenArray *gU, *gV;
  
  // make sure U and V are compatible vectors
  if (!xf_CompatibleVectors(V,U)) return xf_Error(xf_INCOMPATIBLE);
  
  VariableOrder = ((V->nComp != NULL) && (V->vOrder != NULL));
  
  // need state names for remap
  if ((V->StateName == NULL) || (U->StateName == NULL)){
    xf_printf("Error, vector(s) passed into RemapState has no StateName.\n");
    return xf_Error(xf_EQNSET_ERROR);
  }
  
  ierr = xf_Error(xf_Alloc((void **) &P, U->StateRank, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (k=0; k<U->StateRank; k++){
    s = U->StateName[k];
    len = strlen(s);
    P[k] = -1;
    for (j=0; j<V->StateRank; j++){
      if (strcmp(s, V->StateName[j]) == 0){
        P[k] = j;
        break;
      }
    } // j
    if (P[k] < 0)
      xf_printf("FYI: state = %s not found during remap.  Continuing.\n", U->StateName[k]);
  } // k
  
  
  for (egrp=0; egrp<V->nArray; egrp++){
    gV = V->GenArray+egrp;
    gU = U->GenArray+egrp;
    for (elem=0; elem<gU->n; elem++){
      
      // pull off number of values, r
      r = ( (VariableOrder) ? gU->vr[elem] : gU->r);
      
      // make sure # basis functions is compatible between U and V
      nn = r/U->StateRank;
      if (nn != r/V->StateRank) return xf_Error(xf_INPUT_ERROR);
      
      for (i=0; i<nn; i++){
        for (k=0; k<U->StateRank; k++){
          if (P[k] < 0) continue;  // this state is not in V
          gU->rValue[elem][i*U->StateRank+k] = gV->rValue[elem][i*V->StateRank+k];
        }
      } // i
    } // elem
  } // egrp
  
  xf_Release((void *) P);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FindAdjointVectors
int 
xf_FindAdjointVectors( xf_All *All, xf_Vector *U, const char *OutputList, 
                      enum xfe_Bool UseAllOutputs, enum xfe_Bool ZeroFlag, 
                      int *nPsi, xf_Vector ***pPsi, enum xfe_Bool *Found)
{
  int ierr, len, i, iAdjoint;
  char tail[xf_MAXSTRLEN] = "Adjoint";
  char AdjointName[xf_MAXSTRLEN];
  char **OutputNames;
  enum xfe_Bool WriteAdjoint, found;
  xf_Data *D;
  xf_Vector *Adj;
  
  if (Found != NULL) (*Found) = xfe_False;
  
  if (UseAllOutputs){
    /* Create an adjoint for every output */
    if (All->EqnSet->Outputs == NULL) return xf_Error(xf_INPUT_ERROR);
    (*nPsi) = All->EqnSet->Outputs->nOutput;
    
    ierr = xf_Error(xf_Alloc2((void ***) &OutputNames, (*nPsi), xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    
    for (i=0; i<(*nPsi); i++)
      strcpy(OutputNames[i], All->EqnSet->Outputs->Output[i].Name);
  }
  else{
    /* Only use outputs in OutputList */
    ierr = xf_Error(xf_ScanXStringAlloc(OutputList, xf_MAXSTRLEN, nPsi,
                                        &OutputNames));
    if (ierr != xf_OK) return ierr;
  }
  
  if ((*nPsi) <= 0) return xf_Error(xf_OUT_OF_BOUNDS); // no outputs
  
  // Allocate memory for adjoint pointer
  ierr = xf_Error(xf_Alloc( (void **) pPsi, (*nPsi), sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  
  // should we make the adjoint vector(s) writable?
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "WriteAdjoint", &WriteAdjoint);
  if (ierr != xf_OK) return ierr;
  
  if (Found != NULL) (*Found) = xfe_True;
  for (iAdjoint=0; iAdjoint<(*nPsi); iAdjoint++){
    
    // set adjoint name
    sprintf(AdjointName, "%s_Adjoint", OutputNames[iAdjoint], tail);
    
    ierr = xf_Error(xf_FindSimilarVector(All, U, AdjointName, xfe_True, xfe_True, 
                                         &D, (*pPsi) + iAdjoint, &found));
    if (ierr != xf_OK) return ierr;
    D->ReadWrite = WriteAdjoint;
    if ((Found != NULL) && (!found)) (*Found) = xfe_False;
    
    Adj = (*pPsi)[iAdjoint];
    
    // Set OutputName
    len = strlen(OutputNames[iAdjoint])+1;
    ierr = xf_Error(xf_ReAlloc( (void **) &Adj->OutputName, len, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    strcpy(Adj->OutputName, OutputNames[iAdjoint]);
    
    // Initialize adjoint to zero if requested or if just created
    if (ZeroFlag || (!found)){
      ierr = xf_Error(xf_SetZeroVector(Adj));
      if (ierr != xf_OK) return ierr;
    }
    
    // Set SolverRole
    Adj->SolverRole = xfe_SolverRoleAdjointState;
    
  } // iAdjoint
  
  xf_Release2((void **) OutputNames);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AllocateJacobianMatrix
static int 
xf_AllocateJacobianMatrix( xf_Mesh *Mesh, xf_Vector *U, int StateRank, enum xfe_BasisType *Basis, 
                          int *Order, enum xfe_Bool AllocValue, xf_JacobianMatrix **pR_U){
  
  int ierr, negrp, negrphalo, nelem, nface;
  int egrp, elem, face, i, ord;
  int egN, eN, faceN, nn, nN, r, rN;
  int *nFacep1 = NULL;
  enum xfe_Bool VariableOrder;
  xf_long pos, **ipos;
  real *rtemp;
  xf_JacobianMatrix *R_U;
  
  ierr = xf_Error(xf_CreateJacobianMatrix(pR_U));
  if (ierr != xf_OK) return ierr;
  
  R_U = (*pR_U);
  
  R_U->U = U; // store pointer to state if we ever need it (e.g. lean solvers)
  
  R_U->negrp = negrp = Mesh->nElemGroup; // kept consistent with Mesh->nElemGroup
  
  /* Halo element groups are present in parallel runs */  
  negrphalo = ((Mesh->ParallelInfo == NULL) ? Mesh->nElemGroup : 2*Mesh->nElemGroup);

  R_U->negrphalo = negrphalo;

  R_U->StateRank = StateRank;
  
  // Basis
  ierr = xf_Error(xf_Alloc((void **) &R_U->Basis, negrphalo, sizeof(enum xfe_BasisType)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<negrphalo; i++) R_U->Basis[i] = Basis[i%negrp];
  
  // Order
  ierr = xf_Error(xf_Alloc((void **) &R_U->Order, negrphalo, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<negrphalo; i++) R_U->Order[i] = Order[i%negrp];
  
  // Are we dealing with a variable order?
  VariableOrder = ( (U != NULL) && (U->nComp != NULL) && (U->vOrder != NULL));
  
  if (VariableOrder){ // allocate vnvec
    if (negrphalo != U->nArray) return xf_Error(xf_OUT_OF_BOUNDS);
    ierr = xf_Error(xf_VAlloc2((void ***) &R_U->vnvec, U->nArray, U->nComp, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  else R_U->vnvec = NULL;
  
  ierr = xf_Error(xf_Alloc((void **) &R_U->nvec, negrphalo, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &R_U->egrpN, negrphalo, sizeof(int **)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &R_U->elemN, negrphalo, sizeof(int **)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &R_U->faceN, negrphalo, sizeof(int **)));
  if (ierr != xf_OK) return ierr;
  
  if (AllocValue){
    ierr = xf_Error(xf_Alloc((void **) &R_U->Value, negrphalo, sizeof(real ***)));
    if (ierr != xf_OK) return ierr;
  }
  
  for (egrp=0; egrp<negrphalo; egrp++){
    
    nelem = Mesh->ElemGroup[egrp].nElem;
    
    ierr = xf_Error(xf_Order2nNode(Basis[egrp%negrp], Order[egrp%negrp], &nn));
    if (ierr != xf_OK) return ierr;
    
    R_U->nvec[egrp] = nn;
    r = nn*StateRank;
    
    ierr = xf_Error(xf_ReAlloc((void **) &nFacep1, nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (elem=0; elem<nelem; elem++) nFacep1[elem] = Mesh->ElemGroup[egrp].nFace[elem]+1;
    
    ierr = xf_Error(xf_VAlloc2((void ***) R_U->egrpN+egrp, nelem, 
                               Mesh->ElemGroup[egrp].nFace, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_VAlloc2((void ***) R_U->elemN+egrp, nelem,
                               Mesh->ElemGroup[egrp].nFace, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_VAlloc2((void ***) R_U->faceN+egrp, nelem, 
                               Mesh->ElemGroup[egrp].nFace, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_VAlloc2((void ***) &ipos, nelem, nFacep1, sizeof(xf_long)));
    if (ierr != xf_OK) return ierr;
    
    if (AllocValue){
      ierr = xf_Error(xf_VAlloc2((void ***) R_U->Value+egrp, nelem, nFacep1, sizeof(real *)));
      if (ierr != xf_OK) return ierr;
    }
    
    pos = 0;
    
    // Loop to figure our how much space we need
    for (elem=0; elem<nelem; elem++){
      if (VariableOrder){
        ierr = xf_Error(xf_Order2nNode(Basis[egrp%negrp], U->vOrder[egrp][elem], &nn));
        if (ierr != xf_OK) return ierr;
        R_U->vnvec[egrp][elem] = nn;
        r = nn*StateRank;
      }
      nface = Mesh->ElemGroup[egrp].nFace[elem];
      for (face=-1; face<nface; face++){
        ipos[elem][1+face] = pos;
        if (face == -1)
          pos += r*r;
        else{
          ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, &egN, &eN, &faceN));
          if (ierr != xf_OK) return ierr;
          R_U->egrpN[egrp][elem][face] = egN;
          R_U->elemN[egrp][elem][face] = eN;
          R_U->faceN[egrp][elem][face] = faceN;
          if (egN >= 0){ // egN < 0 means face is on a boundary
            ord = ((VariableOrder) ? U->vOrder[egN][eN] : Order[egN%negrp]);
            ierr = xf_Error(xf_Order2nNode(Basis[egN%negrp], ord, &nN));
            if (ierr != xf_OK) return ierr;
            rN = nN*StateRank;
            pos += r*rN;
          }
        }
      } // face
    } // elem
    
    if (AllocValue){
      ierr = xf_Error(xf_LongAlloc( (void **) &rtemp, pos, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      
      for (elem=0; elem<nelem; elem++){
        nface = Mesh->ElemGroup[egrp].nFace[elem];
        for (face=-1; face<nface; face++){
          R_U->Value[egrp][elem][1+face] = rtemp+ipos[elem][1+face];  
        }
      }
    }
    
    xf_Release2( (void **) ipos);
    
  } // egrp
  
  xf_Release( (void *) nFacep1);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AllocateJacobianMatrixCoarse
int 
xf_AllocateJacobianMatrixCoarse( xf_Mesh *Mesh, int CoarseOrder, 
                                xf_JacobianMatrix *R_U)
{
  int ierr, i, negrp;
  int *Order;
  
  if (R_U->R_Uc != NULL) return xf_Error(xf_INPUT_ERROR);
  
  negrp = R_U->negrp;
  
  ierr = xf_Error(xf_Alloc( (void **) &Order, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<negrp; i++) Order[i] = CoarseOrder;
  
  ierr = xf_Error(xf_AllocateJacobianMatrix(Mesh,  NULL, R_U->StateRank, R_U->Basis,
                                            Order, (R_U->Value != NULL), 
                                            (xf_JacobianMatrix **) &R_U->R_Uc));
  if (ierr != xf_OK) return ierr;
  
  xf_Release((void *) Order);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindJacobianMatrix
int 
xf_FindJacobianMatrix( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *U, 
                      const char *TitleIn, enum xfe_Bool AllocValue, 
                      xf_Data **pR_UData, xf_JacobianMatrix **pR_U, 
                      enum xfe_Bool *Found){
  int ierr, len, i, j;
  const char TitleDefault[] = "R_U";
  const char *Title;
  enum xfe_Bool match, VariableOrder;
  xf_Data *D, *N;
  xf_JacobianMatrix *R_U;
  
  if (Found != NULL) (*Found) = xfe_True;
  Title = ((TitleIn == NULL) ? TitleDefault : TitleIn);
  
  len = strlen(Title);
  
  // Are we dealing with a variable order?
  VariableOrder = ((U->nComp != NULL) && (U->vOrder != NULL));
  
  D = DataSet->Head;
  while (D != NULL){
    N = D->Next;
    if ((strcmp(D->Title, Title) == 0) &
        (D->Type == xfe_JacobianMatrix)){
      R_U = (xf_JacobianMatrix *) D->Data;
      if ((R_U->negrp == Mesh->nElemGroup)){
        match = xfe_True;
        for (i=0; i<R_U->negrp; i++){
          if (VariableOrder){ // account for variable orders
            if (R_U->vnvec==NULL) match = xfe_False;
            else{
              for (j=0; j<U->nComp[i]; j++){
                if (R_U->vnvec[i][j]*R_U->StateRank != U->GenArray[i].vr[j]){
                  match = xfe_False;
                  break;
                }
              } // j
            }
          }
          else{
            if (R_U->nvec[i]*R_U->StateRank != U->GenArray[i].r)
              match = xfe_False;
          }
          if (!match) break;
        } // i
        
        // all-reduce match using  min to check if any procs find false
        ierr = xf_Error(xf_MPI_Allreduce(&match, 1, xfe_SizeInt, xfe_MPI_MIN));
        if (ierr != xf_OK) return ierr;
        
        
        if (match){ // data ranks are the same
          R_U->U = U; // store state that found this Jacobian
          if (pR_UData != NULL) (*pR_UData) = D;
          if (pR_U     != NULL) (*pR_U) = R_U;
          return xf_OK; // done
        }
        else{ // destroy if incorrect rank
          ierr = xf_Error(xf_DestroyDataInSet(DataSet, D));
          if (ierr != xf_OK) return ierr;
        }
      }
    }
    D = N;
  }
  
  if (Found != NULL) (*Found) = xfe_False;
  
  // Jacobian was not found in the DataSet; allocate it
  ierr = xf_Error(xf_AllocateJacobianMatrix( Mesh, U, U->StateRank, U->Basis, U->Order,
                                            AllocValue, pR_U));
  if (ierr != xf_OK) return ierr;
  
  // store vector in a new data structure with Title
  ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_JacobianMatrix,
                                xfe_False, (void *) (*pR_U), pR_UData));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetZeroJacobian
int 
xf_SetZeroJacobian( xf_Mesh *Mesh, xf_JacobianMatrix *R_U)
{ 
  // R_U = 0
  int egrp, elem, face, k;
  int egN, eN, r, rN, nN;
  int negrphalo;
  xf_JacobianMatrix *R_Uc;
  
  /* Return immediately if null or if no value */
  if ((R_U == NULL) || (R_U->Value == NULL)) return xf_OK;
  
  /* Halo element groups are present in parallel runs */  
  negrphalo = ((Mesh->ParallelInfo == NULL) ? Mesh->nElemGroup : 2*Mesh->nElemGroup);
  
  for (egrp = 0; egrp<negrphalo; egrp++){
    r = R_U->StateRank * R_U->nvec[egrp];
    for (elem = 0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      if (R_U->vnvec != NULL) r = R_U->StateRank*R_U->vnvec[egrp][elem];
      for (face = -1; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
        if (face == -1){
          egN = egrp;
          eN  = elem;
        }
        else{
          egN = R_U->egrpN[egrp][elem][face];
          eN  = R_U->elemN[egrp][elem][face];
        }
        if (egN < 0) continue;
        nN = ((R_U->vnvec != NULL) ? R_U->vnvec[egN][eN] : R_U->nvec[egN]);
        rN = R_U->StateRank * nN;
        for (k=0; k<r*rN; k++) R_U->Value[egrp][elem][1+face][k] = 0.0;
      } // face
    } // elem
  } // egrp
  
  R_U->Preconditioner = xfe_PreconditionerNone;
  
  // If Jacobian was recalculated, make sure R_Uc is projected next time it is needed
  if (R_U != NULL) R_U->ProjectionNeeded = xfe_True;
  if ((R_Uc = R_U->R_Uc) != NULL){
    R_Uc->Preconditioner = xfe_PreconditionerNone;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetVectorDOF
int 
xf_GetVectorDOF(xf_Vector *V, int *ndof)
{
  int i, j, ierr;
  
  (*ndof) = 0;
  
  if ((V->Linkage == xfe_LinkageGlobElem) && (V->GenArray != NULL)){
    for (i=0,*ndof=0; i<V->nArraySelf; i++){
      if (V->GenArray[i].vr != NULL){
        for (j=0; j<V->GenArray[i].n; j++)
          *ndof += V->GenArray[i].vr[j]/V->StateRank;
      }
      else *ndof += V->GenArray[i].n*V->GenArray[i].r/V->StateRank; 
    }
  }
  
  // reduce-sum
  ierr = xf_Error(xf_MPI_Allreduce(ndof, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DataSetInfo
int 
xf_DataSetInfo(xf_DataSet *DataSet)
{
  int ierr, count, r, i, j, ndof;
  int myRank, nProc;
  char s[200], sline[400];
  xf_Data *D;
  xf_Vector *V;
  xf_VectorSet *VS;
  xf_Matrix *M;
  xf_JacobianMatrix *R_U;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0){ // only root does the printing
    D = DataSet->Head;
    xf_printf("---- DataSet Info ----\n");
    xf_printf("%12s %8s  %s\n", "[Title]", "[Type]", "[Other Info]");
    count = 0;
    while (D != NULL){
      sprintf(sline, "\0");
      ndof = 0;
      if ((D->Type == xfe_Vector) || (D->Type == xfe_VectorSet)){
        if (D->Type == xfe_VectorSet){
          VS = (xf_VectorSet *) D->Data;
          sprintf(s, "nVector = %d, ", VS->nVector);
          strcat(sline, s);
          if (VS->nVector > 0)
            V = VS->Vector + 0;
          else
            V = NULL;
        }
        else
          V = (xf_Vector *) D->Data;
        if (V != NULL){
          r = ((V->GenArray != NULL) ? V->GenArray[0].r : -1);
          sprintf(s, "Linkage=%s, Size=%s, TimeIndex=%d, varorder=%d, sr=%d, nArray=%d, .r[0]=%d",
                  xfe_LinkageName[V->Linkage], xfe_SizeName[V->Size], V->TimeIndex,
                  V->vOrder != NULL, V->StateRank, V->nArray, r);
          strcat(sline, s);
          ierr = xf_Error(xf_GetVectorDOF(V, &ndof));
          if (ierr != xf_OK) return ierr;
        }
      }
      else if (D->Type == xfe_Matrix){
        M = (xf_Matrix *) D->Data;
        if (M != NULL){
          r = ((M->GenArray != NULL) ? M->GenArray[0].r : -1);
          sprintf(s, "Linkage=%s, LinkageIndex=%d, Order1=%d, Order2=%d, Basis1=%s, Basis2=%s",
                  xfe_LinkageName[M->Linkage], M->LinkageIndex, M->Order1, M->Order2,
                  xfe_BasisName[M->Basis1], xfe_BasisName[M->Basis2]);
          strcat(sline, s);
          if (M->GenArray != NULL) ndof = M->GenArray[0].n*r; 
        } 
      }
      else if (D->Type == xfe_JacobianMatrix){
        R_U = (xf_JacobianMatrix *) D->Data;
        if ((R_U != NULL) && (R_U->negrp > 0)){
          sprintf(s, "negrp=%d, StateRank=%d, Basis[0]=%s, Order[0]=%d",
                  R_U->negrp, R_U->StateRank, 
                  ((R_U->Basis == NULL) ? "NULL" : xfe_BasisName[R_U->Basis[0]]),
                  ((R_U->Order == NULL) ? -1 : R_U->Order[0]));
          strcat(sline, s);
          for (i=0,ndof=0; i<R_U->negrp; i++)
            ndof += (R_U->nvec[i]*R_U->StateRank)*(R_U->nvec[i]*R_U->StateRank);
        }
      }
      
      xf_printf("%12s %8s  %s\n", D->Title, xfe_DataName[D->Type], sline);
      if (ndof > 0) xf_printf("%12s ndof = %d\n", "\0", ndof);
      
      D = D->Next;
      count++;
    }
    xf_printf("%d total piece(s) of data.\n\n", count);
    
  } // if root
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_DataSetDeleteNonEssential
int 
xf_DataSetDeleteNonEssential(xf_DataSet *DataSet)
{
  int ierr;
  enum xfe_Bool Remove;
  xf_Data *D, *N;
  xf_Vector *V;
  xf_VectorSet *VS;
  
  D = DataSet->Head;
  while (D != NULL){
    
    N = D->Next;
    Remove = xfe_True;
    
    if (D->Type == xfe_Vector){
      V = (xf_Vector *) D->Data;
      if ((V != NULL) && (V->SolverRole != xfe_SolverRoleNone))
        Remove = xfe_False;
    }
    if (D->Type == xfe_VectorSet){
      VS = (xf_VectorSet *) D->Data;
      if ((VS != NULL) && (VS->nVector > 0) && (VS->Vector[0].SolverRole != xfe_SolverRoleNone))
        Remove = xfe_False;
    }
    
    // delete data if not essential
    if (Remove){
      ierr = xf_Error(xf_DestroyDataInSet(DataSet, D));
      if (ierr != xf_OK) return ierr;
    }
    
    D = N;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BuildVOrder
int 
xf_BuildVOrder(xf_All *All, xf_Vector *U, xf_Vector **pV)
{
  int ierr, i, j;
  enum xfe_Bool Interpolated;
  
  Interpolated  = ((U->Basis != NULL) && ( U->Order != NULL));
  
  if (!Interpolated) return xf_Error(xf_INPUT_ERROR);
  
  ierr = xf_Error(xf_FindVector(All, "VOrder", xfe_LinkageGlobElem,
                                1, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, 
                                xfe_False,  xfe_False, NULL, 
                                pV, NULL));
  if (ierr != xf_OK) return ierr;
  
  // loop over arrays
  for (i=0; i<U->nArray; i++)
    for (j=0; j<U->GenArray[i].n; j++)
      (*pV)->GenArray[i].iValue[j][0] = xf_InterpOrder(U,i,j);
  
  return xf_OK;
}


/*-----------------*/
/* Parallelization */
/*-----------------*/

#include "xf_DataParallel.c"


/*-------------------*/
/* Binary read/write */
/*-------------------*/

/******************************************************************/
//   FUNCTION Definition: xf_WriteGenArrayBinary
static int 
xf_WriteGenArrayBinary( xf_GenArray *G, FILE *fid)
{
  int ierr, i, si, sr, rev, r;
  enum xfe_Bool flag;
  
  si = sizeof(int);
  sr = sizeof(real);
  
  rev = 1;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Size
  ierr = xf_Error(xf_WriteStringBinary(xfe_SizeName[G->Size], fid));
  if (ierr != xf_OK) return ierr;
  
  // n, r
  if (fwrite(&G->n, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  if (fwrite(&G->r, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // vr
  flag = (G->vr != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    if (fwrite(G->vr, si, G->n, fid) != G->n) return xf_Error(xf_FILE_WRITE_ERROR);
  }
  
  // iValue
  flag = (G->iValue != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    for (i=0; i<G->n; i++){
      r = ((G->vr == NULL) ? G->r : G->vr[i]);
      if (fwrite(G->iValue[i], si, r, fid) != r) 
        return xf_Error(xf_FILE_WRITE_ERROR);
    }
  }
  
  // rValue
  flag = (G->rValue != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    for (i=0; i<G->n; i++){
      r = ((G->vr == NULL) ? G->r : G->vr[i]);
      if (fwrite(G->rValue[i], sr, r, fid) != r) 
        return xf_Error(xf_FILE_WRITE_ERROR);
    }
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadGenArrayBinary
static int 
xf_ReadGenArrayBinary( FILE *fid, xf_GenArray *G){
  int ierr, i, si, sr, rev, r;
  enum xfe_Bool flag;
  
  si = sizeof(int);
  sr = sizeof(real);
  
  rev = 0;  // writer revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev > 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // Size
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_SizeName, xfe_SizeLast, 
                                    (int *) &G->Size));
  if (ierr != xf_OK) return ierr;
  
  // n, r
  if (fread(&G->n, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (fread(&G->r, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // vr
  G->vr = NULL;
  if (rev >= 1){
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
    if (ierr != xf_OK) return ierr;
    if (flag){
      ierr = xf_Error(xf_Alloc((void **) &G->vr, G->n, si));
      if (ierr != xf_OK) return ierr;
      if (fread(G->vr, si, G->n, fid) != G->n) return xf_Error(xf_FILE_READ_ERROR);
    }
  }
  
  
  // iValue
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    if (G->vr == NULL){
      ierr = xf_Error(xf_Alloc2((void ***) &G->iValue, G->n, G->r, si));
      if (ierr != xf_OK) return ierr;
    }
    else{
      ierr = xf_Error(xf_VAlloc2((void ***) &G->iValue, G->n, G->vr, si));
      if (ierr != xf_OK) return ierr;
    }
    for (i=0; i<G->n; i++){
      r = ((G->vr == NULL) ? G->r : G->vr[i]);
      if (fread(G->iValue[i], si, r, fid) != r) 
        return xf_Error(xf_FILE_READ_ERROR);
    }
  }
  else
    G->iValue = NULL;
  
  // rValue
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    if (G->vr == NULL){
      ierr = xf_Error(xf_Alloc2((void ***) &G->rValue, G->n, G->r, sr));
      if (ierr != xf_OK) return ierr;
    }
    else{
      ierr = xf_Error(xf_VAlloc2((void ***) &G->rValue, G->n, G->vr, sr));
      if (ierr != xf_OK) return ierr;
    } 
    for (i=0; i<G->n; i++){
      r = ((G->vr == NULL) ? G->r : G->vr[i]);
      if (fread(G->rValue[i], sr, r, fid) != r) 
        return xf_Error(xf_FILE_READ_ERROR);
    }
  }
  else
    G->rValue = NULL;
  
  // ParallelInfo
  G->ParallelInfo = NULL;
  
  return xf_OK;
}





/******************************************************************/
//   FUNCTION Definition: xf_WriteVectorBinarySerial
static int 
xf_WriteVectorBinarySerial( xf_Vector *V, FILE *fid){
  int ierr, i, rev, si;
  enum xfe_Bool flag;
  
  si = sizeof(int);
  
  rev = 1;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Linkage
  ierr = xf_Error(xf_WriteStringBinary(xfe_LinkageName[V->Linkage], fid));
  if (ierr != xf_OK) return ierr;
  
  // Linkage Index
  if (fwrite(&V->LinkageIndex, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // SolverRole
  ierr = xf_Error(xf_WriteStringBinary(xfe_SolverRoleName[V->SolverRole], fid));
  if (ierr != xf_OK) return ierr;
  
  // State Rank
  if (fwrite(&V->StateRank, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // State Name
  flag = (V->StateName != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    for (i=0; i<V->StateRank; i++){
      ierr = xf_Error(xf_WriteStringBinary(V->StateName[i], fid));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // Output Name
  flag = (V->OutputName != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteStringBinary(V->OutputName, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // Time Index
  if (fwrite(&V->TimeIndex, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // MG Index
  if (fwrite(&V->MGIndex, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // nArray
  if (fwrite(&V->nArray, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Basis
  flag = (V->Basis != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    for (i=0; i<V->nArray; i++){
      ierr = xf_Error(xf_WriteStringBinary(xfe_BasisName[V->Basis[i]], fid));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // Order
  flag = (V->Order != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    if (fwrite(V->Order, si, V->nArray, fid) != V->nArray) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }
  
  // nComp
  flag = (V->nComp != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    if (fwrite(V->nComp, si, V->nArray, fid) != V->nArray) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }
  
  // vOrder
  flag = (V->vOrder != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    for (i=0; i<V->nArray; i++){
      if (fwrite(V->vOrder[i], si, V->nComp[i], fid) != V->nComp[i]) 
        return xf_Error(xf_FILE_WRITE_ERROR);
    }
  }
  
  // GenArray
  flag = (V->GenArray != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    for (i=0; i<V->nArray; i++){
      ierr = xf_Error(xf_WriteGenArrayBinary(V->GenArray+i, fid));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteVectorBinary
static int 
xf_WriteVectorBinary( xf_Mesh *Mesh, xf_Vector *V, FILE *fid){
  int ierr, terr;
  int myRank, nProc;
  enum xfe_Bool ParallelFlag;
  xf_Vector *V_Glob = NULL;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if ((Mesh != NULL) && (Mesh->ParallelInfo != NULL))
    ParallelFlag = xfe_True;
  else
    ParallelFlag = xfe_False;
  
  // unparallelize V -> V_Glob
  if (ParallelFlag){
    if (myRank == 0){
      ierr = xf_Error(xf_Error(xf_CreateVector(&V_Glob)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_UnParallelizeVector(Mesh, V, V_Glob));
    if (ierr != xf_OK) return ierr;
  }
  else
    V_Glob = V;
  
  // Root writes V_Glob
  if (myRank == 0)
    terr = xf_Error(xf_WriteVectorBinarySerial(V_Glob, fid));
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);
  
  
  if ((ParallelFlag) && (myRank == 0)){
    ierr = xf_Error(xf_DestroyVector(V_Glob, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadVectorBinarySerial
static int 
xf_ReadVectorBinarySerial(  FILE *fid, xf_Vector *V){
  int ierr, i, rev, si;
  enum xfe_Bool flag;
  enum xfe_SizeType Size = xfe_SizeReal;
  
  si = sizeof(int);
  
  rev = 0;  // reader revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev > 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // Linkage
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_LinkageName, xfe_LinkageLast, 
                                    (int *) &V->Linkage));
  if (ierr != xf_OK) return ierr;
  
  // Linkage Index
  if (fread(&V->LinkageIndex, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // SolverRole
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_SolverRoleName, xfe_SolverRoleLast, 
                                    (int *) &V->SolverRole));
  if (ierr != xf_OK) return ierr;
  
  // State Rank
  if (fread(&V->StateRank, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // State Name
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_Alloc2((void ***) &V->StateName, V->StateRank, xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<V->StateRank; i++){
      ierr = xf_Error(xf_ReadStringBinary(fid, xf_MAXSTRLEN, V->StateName[i], NULL));
      if (ierr != xf_OK) return ierr;
    }
  }
  else
    V->StateName = NULL;
  
  // Output Name
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &V->OutputName));
    if (ierr != xf_OK) return ierr;
  }
  else
    V->OutputName = NULL;
  
  
  // Time Index
  if (fread(&V->TimeIndex, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // MG Index
  if (fread(&V->MGIndex, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // nArray
  if (fread(&V->nArray, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  V->nArraySelf = V->nArray;
  
  // Basis
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_Alloc((void **) &V->Basis, V->nArray, si));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<V->nArray; i++){
      ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BasisName, xfe_BasisLast, (int *) V->Basis+i));
      if (ierr != xf_OK) return ierr;
    }
  }
  else
    V->Basis = NULL;
  
  // Order
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_Alloc((void **) &V->Order, V->nArray, si));
    if (ierr != xf_OK) return ierr;
    if (fread(V->Order, si, V->nArray, fid) != V->nArray) 
      return xf_Error(xf_FILE_READ_ERROR);
  }
  else
    V->Order = NULL;
  
  // nComp
  if (rev >= 1){
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
    if (ierr != xf_OK) return ierr;
    if (flag){
      ierr = xf_Error(xf_Alloc((void **) &V->nComp, V->nArray, si));
      if (ierr != xf_OK) return ierr;
      if (fread(V->nComp, si, V->nArray, fid) != V->nArray) 
        return xf_Error(xf_FILE_READ_ERROR);
    }
    else
      V->nComp = NULL;
  }
  
  // vOrder
  if (rev >= 1){
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
    if (ierr != xf_OK) return ierr;
    if (flag){
      ierr = xf_Error(xf_VAlloc2((void ***) &V->vOrder, V->nArray, V->nComp, si));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<V->nArray; i++)
        if (fread(V->vOrder[i], si, V->nComp[i], fid) != V->nComp[i]) 
          return xf_Error(xf_FILE_READ_ERROR);
    }
    else
      V->vOrder = NULL;
  }
  
  // GenArray
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_Alloc((void **) &V->GenArray, V->nArray, sizeof(xf_GenArray)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<V->nArray; i++){
      ierr = xf_Error(xf_ReadGenArrayBinary(fid, V->GenArray+i));
      if (ierr != xf_OK) return ierr;
      if (i == 0) Size = V->GenArray[i].Size;
      else if (Size != V->GenArray[i].Size) return xf_Error(xf_FILE_READ_ERROR);
    }
    V->Size = Size;
  }
  else
    V->GenArray = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadVectorBinary
static int 
xf_ReadVectorBinary(  xf_Mesh *Mesh, FILE *fid, xf_Vector *V){
  int ierr, terr;
  int myRank;
  xf_Vector *V_Glob = NULL;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  // root reads the vector from the file
  if (myRank == 0){
    terr = xf_Error(xf_CreateVector(&V_Glob));
    if (terr == xf_OK)
      terr = xf_Error(xf_ReadVectorBinarySerial(fid, V_Glob));
  }
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);
  
  // in serial, V just gets all of V_Glob's contents
  ierr = xf_Error(xf_ParallelizeVector(Mesh, V_Glob, V));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(V_Glob, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteMatrixBinary
static int 
xf_WriteMatrixBinary( xf_Mesh *Mesh, xf_Matrix *M, FILE *fid){
  int ierr, rev, si;
  int myRank, nProc;
  enum xfe_Bool flag;
  
  si = sizeof(int);
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // matrix writing in parallel is not supported yet (have not needed it)
  if (nProc > 1) return xf_Error(xf_NOT_SUPPORTED);
  
  
  rev = 0;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Linkage
  ierr = xf_Error(xf_WriteStringBinary(xfe_LinkageName[M->Linkage], fid));
  if (ierr != xf_OK) return ierr;
  
  // Linkage Index
  if (fwrite(&M->LinkageIndex, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Order1, Order2
  if (fwrite(&M->Order1, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  if (fwrite(&M->Order2, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Basis1, Basis2
  ierr = xf_Error(xf_WriteStringBinary(xfe_BasisName[M->Basis1], fid));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_WriteStringBinary(xfe_BasisName[M->Basis2], fid));
  if (ierr != xf_OK) return ierr;
  
  // GenArray
  flag = (M->GenArray != NULL);
  ierr = xf_Error(xf_WriteStringBinary(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_WriteGenArrayBinary(M->GenArray, fid));
    if (ierr != xf_OK) return ierr;
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadMatrixBinary
static int 
xf_ReadMatrixBinary( xf_Mesh *Mesh, FILE *fid, xf_Matrix *M){
  int ierr, rev, si;
  int myRank, nProc;
  enum xfe_Bool flag;
  
  si = sizeof(int);
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // matrix reading in parallel is not supported yet (have not needed it)
  if (nProc > 1) return xf_Error(xf_NOT_SUPPORTED);
  
  rev = 0;  // writer revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // Linkage
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_LinkageName, xfe_LinkageLast, 
                                    (int *) &M->Linkage));
  if (ierr != xf_OK) return ierr;
  
  // Linkage Index
  if (fread(&M->LinkageIndex, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // Order1, Order2
  if (fread(&M->Order1, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (fread(&M->Order2, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // Basis1, Basis2
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BasisName, xfe_BasisLast, (int *) &M->Basis1));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BasisName, xfe_BasisLast, (int *) &M->Basis2));
  if (ierr != xf_OK) return ierr;
  
  
  // GenArray
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  if (flag){
    ierr = xf_Error(xf_Alloc((void **) &M->GenArray, 1, sizeof(xf_GenArray)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadGenArrayBinary(fid, M->GenArray));
    if (ierr != xf_OK) return ierr;
  }
  else
    M->GenArray = NULL;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteVectorSetBinary
static int 
xf_WriteVectorSetBinary( xf_Mesh *Mesh, xf_VectorSet *VS, FILE *fid){
  int ierr, i, si, rev;
  
  si = sizeof(int);
  
  rev = 0;  // writer revision number
  ierr = xf_Error(xf_fwrite(&rev, sizeof(int), 1, fid));
  if (ierr != xf_OK) return ierr;
  
  // number of vectors in set
  ierr = xf_Error(xf_fwrite(&VS->nVector, sizeof(int), 1, fid));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<VS->nVector; i++){
    ierr = xf_Error(xf_WriteVectorBinary(Mesh, VS->Vector+i, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadVectorSetBinary
static int 
xf_ReadVectorSetBinary( xf_Mesh *Mesh, FILE *fid, xf_VectorSet **pVS){
  int ierr, i, rev, si, nVector;
  enum xfe_Bool flag;
  
  si = sizeof(int);
  
  // read revision number
  ierr = xf_Error(xf_fread(fid, sizeof(int), 1, &rev));
  if (ierr != xf_OK) return ierr;
  
  // revision # check
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // nVector
  ierr = xf_Error(xf_fread(fid, sizeof(int), 1, &nVector));
  if (ierr != xf_OK) return ierr;
  
  // create (*pVS)
  ierr = xf_Error(xf_CreateVectorSet(nVector, pVS));
  if (ierr != xf_OK) return ierr;
  
  // read vectors
  for (i=0; i<(*pVS)->nVector; i++){
    ierr = xf_Error(xf_ReadVectorBinary(Mesh, fid, (*pVS)->Vector+i));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteDataBinary
static int 
xf_WriteDataBinary( xf_Mesh *Mesh, xf_Data *D, FILE *fid){
  int ierr, rev;
  
  rev = 0;  // writer revision number
  ierr = xf_Error(xf_fwrite(&rev, sizeof(int), 1, fid));
  if (ierr != xf_OK) return ierr;
  
  // Title
  ierr = xf_Error(xf_WriteStringBinaryParallel(D->Title, fid));
  if (ierr != xf_OK) return ierr;
  
  // Type
  ierr = xf_Error(xf_WriteStringBinaryParallel(xfe_DataName[D->Type], fid));
  if (ierr != xf_OK) return ierr;

  switch (D->Type){
    case xfe_Vector:
      ierr = xf_Error(xf_WriteVectorBinary(Mesh, (xf_Vector *) D->Data, fid));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_Matrix:
      ierr = xf_Error(xf_WriteMatrixBinary(Mesh, (xf_Matrix *) D->Data, fid));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_VectorSet:
      ierr = xf_Error(xf_WriteVectorSetBinary(Mesh, (xf_VectorSet *) D->Data, fid));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadDataBinary
static int 
xf_ReadDataBinary( xf_Mesh *Mesh, FILE *fid, xf_Data *D)
{
  int ierr, rev;
  xf_VectorSet *VS;
  xf_Vector *V;
  xf_Matrix *M;
  
  // read revision number
  ierr = xf_Error(xf_fread(fid, sizeof(int), 1, &rev));
  if (ierr != xf_OK) return ierr;
  
  // revision # check
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // Title
  ierr = xf_Error(xf_ReadStringBinaryParallel(fid, -1, NULL, &D->Title));
  if (ierr != xf_OK) return ierr;
  
  // Type
  ierr = xf_Error(xf_ReadEnumBinaryParallel(fid, xfe_DataName, xfe_DataLast, 
                                            (int *) &D->Type));
  if (ierr != xf_OK) return ierr;
  
  switch (D->Type){
    case xfe_Vector:
      ierr = xf_Error(xf_CreateVector(&V));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReadVectorBinary(Mesh, fid, V));
      if (ierr != xf_OK) return ierr;
      D->Data = V;
      break;
    case xfe_Matrix:
      ierr = xf_Error(xf_CreateMatrix(&M));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReadMatrixBinary(Mesh, fid, M));
      if (ierr != xf_OK) return ierr;
      D->Data = M;
      break;
    case xfe_VectorSet:
      ierr = xf_Error(xf_ReadVectorSetBinary(Mesh, fid, &VS));
      if (ierr != xf_OK) return ierr;
      D->Data = VS;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}





/******************************************************************/
//   FUNCTION Definition: xf_WriteDataSetBinary
int 
xf_WriteDataSetBinary( xf_Mesh *Mesh, xf_DataSet *DataSet, FILE *fidin,
                      const char *fname)
{
  int ierr, rev, count, si, myRank, len;
  xf_long pos, pos2;
  enum xfe_Bool flag, match, done;
  xf_Data *D;
  char Title[xf_MAXSTRLEN];
  FILE *fid = NULL;
  
  si = sizeof(int);
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank,NULL));
  if (ierr != xf_OK) return ierr;
  
  /* if fname is specified, file is opened */
  if (fname != NULL){
    ierr = xf_Error(xf_fopen(fname, "wb", &fid));
    if (ierr != xf_OK) return ierr;
  }
  else
    fid = fidin;
  
  rev = 0;  // writer revision number
  ierr = xf_Error(xf_fwrite(&rev, sizeof(int), 1, fid));
  if (ierr != xf_OK) return ierr;
  
  
  flag = (DataSet->Head != NULL);
  ierr = xf_Error(xf_WriteStringBinaryParallel(xfe_BoolName[flag], fid));
  if (ierr != xf_OK) return ierr;
  if (!flag) return xf_OK;  // if Head is NULL, no need to continue
  
  D = DataSet->Head;
  
  count = 0;
  
  // placeholder for number of Data nodes written
  ierr = xf_Error(xf_ftell(fid, &pos));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_fwrite(&count, sizeof(int), 1, fid));
  if (ierr != xf_OK) return ierr;
  
  done = xfe_False;
  
  while (!done){
    //make sure data D is the same on all procs
    if (myRank == 0) len = strlen(D->Title)+1;
    
    ierr = xf_Error(xf_MPI_Bcast((void *) &len, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // bcast Title on root
    if (myRank == 0){
      strcpy(Title, D->Title);
      match = xfe_True;
    }
    else
      match = xfe_False;
      
    ierr = xf_Error(xf_MPI_Bcast((void **)&Title, len*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    
    // If not root proc, search for data with Title
    if (myRank != 0){
      D = DataSet->Head;
      while (D != NULL){
        if (strcmp(D->Title, Title) == 0){
          match = xfe_True;
          break;
        }
        D = D->Next;
      }
    }
    
    // Find minimum value of match on all procs; should be 1 if data was found on all
    ierr = xf_Error(xf_MPI_Allreduce((int*)&match, 1, xfe_SizeInt, xfe_MPI_MIN));
    if (ierr != xf_OK) return ierr;
    
    //Write data to file; throw error if data was not found on all procs
    if (D->ReadWrite){
      if (!match){
       return xf_Error(xf_PARALLEL_ERROR); 
      } 
      ierr = xf_Error(xf_WriteDataBinary(Mesh, D, fid)); 
      if (ierr != xf_OK) return ierr;
      
      count++;
    }
      
    D = D->Next; //move on to next piece of data
    
    if (myRank == 0 && D == NULL) done = xfe_True;
    
    ierr = xf_Error(xf_MPI_Allreduce((int*)&done, 1, xfe_SizeInt, xfe_MPI_MAX));
    if (ierr != xf_OK) return ierr;
    
  } //end while !done
  
  // get current position
  ierr = xf_Error(xf_ftell(fid, &pos2));
  if (ierr != xf_OK) return ierr;
  
  // seek to pos
  ierr = xf_Error(xf_fseek(fid, pos, SEEK_SET));
  if (ierr != xf_OK) return ierr;
  
  // write count
  ierr = xf_Error(xf_fwrite(&count, sizeof(int), 1, fid));
  if (ierr != xf_OK) return ierr;
  
  // seek back to pos2
  ierr = xf_Error(xf_fseek(fid, pos2, SEEK_SET));
  if (ierr != xf_OK) return ierr;
  
  // write end
  ierr = xf_Error(xf_WriteStringBinaryParallel("END DATASET", fid));
  if (ierr != xf_OK) return ierr;
  
  // close file if opened one
  if (fname != NULL){
    ierr = xf_Error(xf_fclose(fid));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DumpVectorBinary
int 
xf_DumpVectorBinary( xf_Mesh *Mesh, const char *Name,
                    xf_Vector *V, const char *fname)
{
  int ierr;
  xf_DataSet *DataSet;
  xf_Data *D;
  
  /* Create a dataset for writing the vector */
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DataSetAdd(DataSet, Name, xfe_Vector,
                                xfe_True, (void *) V, &D));
  if (ierr != xf_OK) return ierr;
  
  // write out data set
  ierr = xf_Error(xf_WriteDataSetBinary(Mesh, DataSet, NULL, fname));
  if (ierr != xf_OK) return ierr;
  
  // destroy DataSet
  D->Data = NULL;
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: Yu_DumpMultiVectorBinary
int 
Yu_DumpMultiVectorBinary( xf_Mesh *Mesh, const char *Name1,
                    xf_Vector *V1, const char *Name2, xf_Vector *V2, 
                    const char *fname)
{
  int ierr;
  xf_DataSet *DataSet;
  xf_Data *D1, *D2;
  
  /* Create a dataset for writing the vector */
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;

  //add vector data one-by-one
  ierr = xf_Error(xf_DataSetAdd(DataSet, Name1, xfe_Vector,
                                xfe_True, (void *) V1, &D1));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DataSetAdd(DataSet, Name2, xfe_Vector,
                                xfe_True, (void *) V2, &D2));
  if (ierr != xf_OK) return ierr;
  
  // write out data set
  ierr = xf_Error(xf_WriteDataSetBinary(Mesh, DataSet, NULL, fname));
  if (ierr != xf_OK) return ierr;
  
  // destroy DataSet
  D1->Data = NULL;
  D2->Data = NULL;
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}
/******************************************************************/
//   FUNCTION Definition: Yu_DumpMulti3VectorBinary
int 
Yu_DumpMulti3VectorBinary( xf_Mesh *Mesh, const char *Name1,
                    xf_Vector *V1, const char *Name2, xf_Vector *V2, 
                    const char *Name3, xf_Vector *V3, const char *fname)
{
  int ierr;
  xf_DataSet *DataSet;
  xf_Data *D1, *D2, *D3;
  
  /* Create a dataset for writing the vector */
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;

  //add vector data one-by-one
  ierr = xf_Error(xf_DataSetAdd(DataSet, Name1, xfe_Vector,
                                xfe_True, (void *) V1, &D1));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DataSetAdd(DataSet, Name2, xfe_Vector,
                                xfe_True, (void *) V2, &D2));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DataSetAdd(DataSet, Name3, xfe_Vector,
                                xfe_True, (void *) V3, &D3));
  if (ierr != xf_OK) return ierr;
  
  // write out data set
  ierr = xf_Error(xf_WriteDataSetBinary(Mesh, DataSet, NULL, fname));
  if (ierr != xf_OK) return ierr;
  
  // destroy DataSet
  D1->Data = NULL;
  D2->Data = NULL;
  D3->Data = NULL;
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadDataSetBinary
int 
xf_ReadDataSetBinary( xf_Mesh *Mesh, FILE *fidin, const char *fname,
                     xf_DataSet *DataSet){
  int ierr, rev, i, count;
  int terr;
  enum xfe_Bool flag;
  char s[xf_MAXSTRLEN];
  FILE *fid = NULL;
  xf_Data *D, *P;
  
  /* if fname is specified, file is opened */
  if (fname != NULL){
    ierr = xf_fopen(fname, "rb", &fid);
    if (ierr != xf_OK) return ierr;
  }
  else
    fid = fidin;
  
  // read revision number
  ierr = xf_Error(xf_fread(fid, sizeof(int), 1, &rev));
  if (ierr != xf_OK) return ierr;
  
  // revision # check
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // Initialize Head and Tail
  DataSet->Head = DataSet->Tail = NULL;
  
  // Check if no data
  ierr = xf_Error(xf_ReadEnumBinaryParallel(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
  if (ierr != xf_OK) return ierr;
  
  if (!flag) return xf_OK; // no need to continue  
  
  // read number of data nodes
  ierr = xf_Error(xf_fread(fid, sizeof(int), 1, &count));
  if (ierr != xf_OK) return ierr;
  
  if (count < 0) return xf_Error(xf_OUT_OF_BOUNDS);
  
  P = NULL;
  for (i=0; i<count; i++){
    ierr = xf_Error(xf_CreateData(&D));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadDataBinary(Mesh, fid, D));
    if (ierr != xf_OK) return ierr;
    D->ReadWrite = xfe_True;
    
    if (P == NULL) DataSet->Head = D;
    else P->Next = D;
    D->Prev = P;
    P = D;
  }
  DataSet->Tail = P;
  
  // read end
  ierr = xf_Error(xf_ReadStringBinaryParallel(fid, xf_MAXSTRLEN, s, NULL));
  if (ierr != xf_OK) return ierr;
  
  // close file if opened one
  if (fname != NULL){
    ierr = xf_Error(xf_fclose(fid));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_FindDataByPointer
int 
xf_FindDataByPointer( xf_DataSet *DataSet, void *TargetData, xf_Data **Data)
{  
  xf_Data *D = DataSet->Head;
  enum xfe_Bool found = xfe_False;
  
  while (D != NULL) {
    if (TargetData == D->Data) {//same piece of memory
      found = xfe_True;
      (*Data) = D;
      break;
    }
    D = D->Next;
  }
  if (!found)
    return xf_NOT_FOUND;
  
  return xf_OK;
}

/*****************************************************************/
int
xfYu_FindOrCreatePrimalState( xf_All *All, Yu_Model Model, enum xfe_Bool RestartFlag, xf_ICs *ICsOrig,
                             xf_Vector **pU)
{
    int i, ierr, Order, StateRank;
    enum xfe_BasisType Basis;
    enum xfe_Bool VariableBasis, VariableOrder, ProjectionRequired;
    enum xfe_Bool VisBad = xfe_True;
    enum xfe_Bool UsingVOrder = xfe_False;
    char value[xf_MAXSTRLEN];
    char VOrderFile[xf_MAXSTRLEN];
    xf_Data *StateData;
    xf_Vector *U, *V;

    /*---------------------------------------------*/
    //specify the Basis type, order and stateRank
    for(i=0; i<xfe_BasisLast; i++)
        if(strcmp(Model.basis, xfe_BasisName[i])==0)
        {  Basis = i; break;}
    
    if(i == xfe_BasisLast) xf_Error(xf_OUT_OF_BOUNDS);
    
    Order = Model.order;
    StateRank = Model.nVars;
    if(Model.Dyn_p_Adapt)
       VariableOrder = xfe_True;
    else
       VariableOrder = xfe_False;
    
    /*-----------------------------------------------------*/
    /* Always create and initialize a primal state vector */
    /*-----------------------------------------------------*/   
    /* Create a new glob real elem Vec linked to Mesh */
    ierr = xf_Error(xf_NewGREVector(All, Basis, Order, StateRank, 
                                    xfe_True, pU));
    if (ierr != xf_OK) return ierr;
    U = (*pU);
    U->SolverRole = xfe_SolverRolePrimalState;  
    
   /* Are we asking for variable orders from a file? */
   ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "VOrderFile", VOrderFile));
   if (ierr != xf_OK) return ierr;
   
   /* If so, project to variable order */
   UsingVOrder = xfe_False;
   if (xf_NotNull(VOrderFile)){
      ierr = xf_Error(xf_ProjectVectors_VOrderFile(All->Mesh, All->DataSet, 1, &U, VOrderFile));
      if (ierr != xf_OK) return ierr;
      UsingVOrder = xfe_True;
   }
    /* Initialize basic state (no fancy solves or data altering) */
    ierr = xf_Error(xfYu_InitAlterState(All, &Model, Model.nameVars, Model.initVars, U, StateRank)); 
    if (ierr != xf_OK) return ierr;
   
    /*--------------------------------*/
    /* Find existing state and remap  */
    /*--------------------------------*/
    
    ierr = xf_FindPrimalState(All->DataSet, 0, &StateData, NULL);
    if (ierr == xf_NOT_FOUND){

        if (RestartFlag){
            /* A TimeIndex==0 primal state is required for a restart */
            xf_printf("Error, Restart requested but no compatible primal state was found.\n");
            return xf_Error(ierr);
        }
    }
    else if (ierr != xf_OK) return xf_Error(ierr);
    else
    {
        // state was found
        V = (xf_Vector *) StateData->Data; // V is the old state
        
        ProjectionRequired = xfe_False;
        
        // check if need to project by comparing found vector (V) to generated one (U)
        if (!xf_CompatibleVectors(U, V))
           ProjectionRequired = xfe_True;
       
        //only order variation supported
        //xf_printf("Only order variation is not supported for restarting!~\n");
        //VariableOrder = xfe_False;

        if (ProjectionRequired){
            
              if ((RestartFlag) && (VariableOrder)){
                // if asking to use variable order on a restart, make the
                // newly-created vector, U, look like the found one V in terms
                // of order.  This has the effect of starting from V when
                // combined with the subsequent projection. 
                ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, U, V));
                if (ierr != xf_OK) return ierr;
                xf_printf("Restarting from an existing variable order state.\n");
              }
              else{
                 xf_printf("Projecting state to a new Basis and/or Order.\n");
                 xf_printf("To prevent projection, specify the same basis/order.\n");
              }
            
              ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, V, U));
              if (ierr != xf_OK) return ierr;
        }
         
        ierr = xf_Error(xf_RemapState(V, U));
        if (ierr != xf_OK){
            if (RestartFlag){
                xf_printf("Error, state remap failed.  Most likely, this means that the\n");
                xf_printf("existing restart state is incompatible with the equation set.\n");
                return xf_Error(ierr);
            }
            else VisBad = xfe_True;
        } 

    }
    
    /* Delete any existing states */
    ierr = xf_Error(xf_DataSetRemove(All->DataSet, "State", xfe_True));
    if (ierr != xf_OK) return ierr;
    
    // store U vector in a new data structure with Title="State"
    ierr = xf_Error(xf_DataSetAdd(All->DataSet, "State", xfe_Vector,
                                  xfe_True, (void *) U, NULL));
    if (ierr != xf_OK) return ierr;
    
    return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_FindOrCreatePrimalState
int 
xf_FindOrCreatePrimalState( xf_All *All, enum xfe_Bool RestartFlag, xf_ICs *ICsOrig, 
                           xf_Vector **pU)
{
    int i, ierr, Order, StateRank;
    enum xfe_BasisType Basis;
    enum xfe_Bool VariableBasis, VariableOrder, ProjectionRequired;
    enum xfe_Bool VisBad = xfe_True;
    enum xfe_Bool UsingVOrder = xfe_False;
    char value[xf_MAXSTRLEN];
    char VOrderFile[xf_MAXSTRLEN];
    xf_Data *StateData;
    xf_Vector *U, *V;
    
    
    /*-------------------------------------------------*/
    /* Get interpolation Basis + Order, and StateRank  */
    /*-------------------------------------------------*/
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "InterpBasis", value));
    if (ierr != xf_OK) return ierr;
    if (strncmp(value, "Variable", 8) == 0){
        VariableBasis = xfe_True;
    }
    else{
        VariableBasis = xfe_False;
        ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "InterpBasis", 
                                           xfe_BasisName, (int ) xfe_BasisLast, (int *) &Basis));
        if (ierr != xf_OK) return ierr;
    }
    
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "InterpOrder", value));
    if (ierr != xf_OK) return ierr;
    if (strncmp(value, "Variable", 8) == 0){
        VariableOrder = xfe_True;
        Order = 0;
    }
    else{
        VariableOrder = xfe_False;
        ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "InterpOrder", (int *) &Order));
        if (ierr != xf_OK) return ierr;
    }
    
    StateRank = All->EqnSet->StateRank;
    
    
    /*-----------------------------------------------------*/
    /* Always create and initialize a primal state vector */
    /*-----------------------------------------------------*/
    
    /* Create a new glob real elem Vec linked to Mesh */
    ierr = xf_Error(xf_NewGREVector(All, Basis, Order, StateRank, 
                                    xfe_True, pU));
    if (ierr != xf_OK) return ierr;
    U = (*pU);
    U->SolverRole = xfe_SolverRolePrimalState;
    
    /* Are we asking for variable orders from a file? */
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "VOrderFile", VOrderFile));
    if (ierr != xf_OK) return ierr;
    
    /* If so, project to variable order */
    UsingVOrder = xfe_False;
    if (xf_NotNull(VOrderFile)){
        ierr = xf_Error(xf_ProjectVectors_VOrderFile(All->Mesh, All->DataSet, 1, &U, VOrderFile));
        if (ierr != xf_OK) return ierr;
        UsingVOrder = xfe_True;
    }
    
    /* Initialize basic state (no fancy solves or data altering) */
    ierr = xf_Error(xf_InitAlterState(All, NULL, NULL, U));
    if (ierr != xf_OK) return ierr;
    
    
    /*--------------------------------*/
    /* Find existing state and remap  */
    /*--------------------------------*/
    
    /* Note: we look for an existing state if we are "restarting"
     (reusing existing state as is) or "reinitializing from best
     guess" (e.g. for unsteady solves when we need to do a prior
     steady solve for initialization) */
    
    ierr = xf_FindPrimalState(All->DataSet, 0, &StateData, NULL);
    if (ierr == xf_NOT_FOUND){
        if (RestartFlag){
            /* A TimeIndex==0 primal state is required for a restart */
            xf_printf("Error, Restart requested but no compatible primal state was found.\n");
            return xf_Error(ierr);
        }
    }
    else if (ierr != xf_OK) return xf_Error(ierr);
    else{
        // state was found
        V = (xf_Vector *) StateData->Data; // V is the old state
        
        ProjectionRequired = xfe_False;
        
        // check if need to project by comparing found vector (V) to generated one (U)
        if (!xf_CompatibleVectors(U, V))
            ProjectionRequired = xfe_True;
       
        if (ProjectionRequired){
            
            if ((RestartFlag) && (VariableOrder)){
                /* if asking to use variable order on a restart, make the
                 newly-created vector, U, look like the found one V in terms
                 of order.  This has the effect of starting from V when
                 combined with the subsequent projection. */
                ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, U, V));
                if (ierr != xf_OK) return ierr;
                xf_printf("Restarting from an existing variable order state.\n");
            }
            else{
                xf_printf("Projecting state to a new Basis and/or Order.\n");
                xf_printf("To prevent projection, specify the same basis/order.\n");
            }
            ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, V, U));
            if (ierr != xf_OK) return ierr;
        }
        
        // remap V -> U
        ierr = xf_Error(xf_RemapState(V, U));
        if (ierr != xf_OK){
            if (RestartFlag){
                xf_printf("Error, state remap failed.  Most likely, this means that the\n");
                xf_printf("existing restart state is incompatible with the equation set.\n");
                return xf_Error(ierr);
            }
            else VisBad = xfe_True;
        }
        
        // rescale state if desired
        if ((ICsOrig != NULL) && (!VisBad)){
            xf_printf("Scaling state using ICs.\n");
            ierr = xf_ScaleState(All, ICsOrig, U);
            if (ierr != xf_OK){
                if (RestartFlag){
                    xf_printf("Error, state rescale failed\n");
                    return xf_Error(ierr);
                }
                else VisBad = xfe_True;
            }
        }
        
        /* If unsteady restart, locate/initialize any additional state
         vectors and check their order. */
    }
    
    if ((!RestartFlag) || (All->EqnSet->ICs[0].IC->PriorSteadySolve)){
        if (VariableOrder || VariableBasis){
            xf_printf("Cannot use Variable Basis or Order with Restart = False.");
            return xf_Error(xf_INPUT_ERROR);
        }
        
        /* Not restarting, so initialize state */
        ierr = xf_Error(xf_InitState(All, U));
        if (ierr != xf_OK) return ierr;
        
    }
    
    /* Delete any existing states */
    ierr = xf_Error(xf_DataSetRemove(All->DataSet, "State", xfe_True));
    if (ierr != xf_OK) return ierr;
    
    // store U vector in a new data structure with Title="State"
    ierr = xf_Error(xf_DataSetAdd(All->DataSet, "State", xfe_Vector,
                                  xfe_True, (void *) U, NULL));
    if (ierr != xf_OK) return ierr;
    
    return xf_OK;
}



