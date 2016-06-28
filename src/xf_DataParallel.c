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
 FILE:  xf_DataParallel.c
 
 This file contains parallelization functions for working with the
 Data structure.  It is not compiled itself into an object.  Rather,
 it is included in xf_Data.c.
 
 */


/******************************************************************/
//   FUNCTION Definition: xf_HaloExchangeVectorBegin
int 
xf_HaloExchangeVectorBegin( xf_Vector *V){
  
  int ierr, i, j, k, r, pos, nrecv_tot;
  int myRank, nProc, iProc;
  int nrecv, nsend, recvsize, sendsize, reqsize;
  int *SendElem;
  enum xfe_Bool HaloFlag;
  char *request;
  void *sbuf, *rbuf;
  xf_GenArray *ga;
  xf_ArrayParallelInfo *PInfo;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc == 1 || V->ParallelFlag == xfe_False) return xf_OK; // not parallel, so return immediately
  
  if (V->HaloInTransit){
    xf_printf("Error.  Attempting to halo-exchange a vector while halo is in transit.\n");
    return xf_Error(xf_CODE_LOGIC_ERROR);
  }
  
  // only handle glob elem linkage for now
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  
  reqsize = xf_MPI_RequestSize();
  
  for (i=0; i<V->nArray; i++){
    ga = V->GenArray + i;
    r = ga->r;
    if ((PInfo = ga->ParallelInfo) == NULL) return xf_Error(xf_PARALLEL_ERROR);
    
    HaloFlag = PInfo->HaloFlag;
    
    if (!HaloFlag){ // Normal groups send data
      for (iProc=0; iProc<nProc; iProc++){
        nsend = PInfo->nSendElem[iProc];
        if (nsend > 0){
          SendElem = PInfo->SendElem[iProc];
          if (ga->Size == xfe_SizeInt){
            // pack send buffer
            for (j=0, pos=0; j<nsend; j++){
              r = (ga->vr != NULL) ? ga->vr[SendElem[j]] : ga->r;
              for (k=0; k<r; k++)
                PInfo->iSendBuf[iProc][pos+k] = ga->iValue[SendElem[j]][k];
              pos += r;
            }
            sbuf = (void *) PInfo->iSendBuf[iProc];
            sendsize = pos*sizeof(int);
          }
          else{
            // pack send buffer
            for (j=0, pos=0; j<nsend; j++){
              r = (ga->vr != NULL) ? ga->vr[SendElem[j]] : ga->r;
              for (k=0; k<r; k++){
                PInfo->rSendBuf[iProc][pos+k] = ga->rValue[SendElem[j]][k];
              }
              pos += r;
            }
            sbuf = (void *) PInfo->rSendBuf[iProc];
            sendsize = pos*sizeof(real);
          }
          request = ((char *) PInfo->Request) + iProc*reqsize;
          ierr = xf_Error(xf_MPI_Isend(sbuf, sendsize, iProc, 0, (void *) request));
          if (ierr != xf_OK) return ierr;
        }
      } // iProc
    }
    else{ // Halo groups receive data
      nrecv_tot = 0;
      for (iProc=0; iProc<nProc; iProc++){
        nrecv = PInfo->nRecvElem[iProc];
        if (nrecv > 0){
          pos = nrecv*r;
          if (ga->vr != NULL) for (j=nrecv_tot,pos=0;j<nrecv_tot+nrecv;j++) pos += ga->vr[j];
          if (ga->Size == xfe_SizeInt){
            rbuf = (void *) ga->iValue[0 + nrecv_tot];
            recvsize = pos*sizeof(int);
          }
          else{
            rbuf = (void *) ga->rValue[0 + nrecv_tot];
            recvsize = pos*sizeof(real);
          }
          //xf_pprintf("Beginning recv from iProc= %d, nrecv = %d.\n", iProc, nrecv);
          request = ((char *) PInfo->Request) + iProc*reqsize;
          ierr = xf_Error(xf_MPI_Irecv(rbuf, recvsize, iProc, 0, (void *) request));
          if (ierr != xf_OK) return ierr;
          nrecv_tot += nrecv;
        }
      } // iProc
    }
    
  } // i
  
  // set in-transit flag
  V->HaloInTransit = xfe_True;
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_HaloExchangeVectorEnd
int 
xf_HaloExchangeVectorEnd( xf_Vector *V){
  
  int ierr, i, j , k;
  int myRank, nProc, iProc;
  int nrecv, nsend, reqsize;
  enum xfe_Bool HaloFlag;
  char *request;
  xf_GenArray *ga;
  xf_ArrayParallelInfo *PInfo;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc == 1 || !V->ParallelFlag) return xf_OK; // not parallel, so return immediately
  
  if (!V->HaloInTransit) return xf_OK; // to allow for multiple calls 
  
  // only handle glob elem linkage for now
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  
  reqsize = xf_MPI_RequestSize();
  
  for (i=0; i<V->nArray; i++){
    ga = V->GenArray + i;
    if ((PInfo = ga->ParallelInfo) == NULL) return xf_Error(xf_PARALLEL_ERROR);
    
    HaloFlag = PInfo->HaloFlag;
    
    
    if (!HaloFlag){
      for (iProc=0; iProc<nProc; iProc++){
        nsend = PInfo->nSendElem[iProc];
        if (nsend > 0){
          request = ((char *) PInfo->Request) + iProc*reqsize;
          ierr = xf_Error(xf_MPI_Wait((void *) request)); 
          if (ierr != xf_OK) return ierr;
        }
      } // iProc
    }
    else{
      for (iProc=0; iProc<nProc; iProc++){
        nrecv = PInfo->nRecvElem[iProc];
        if (nrecv > 0){
          request = ((char *) PInfo->Request) + iProc*reqsize;
          ierr = xf_Error(xf_MPI_Wait((void *) request)); 
          if (ierr != xf_OK) return ierr;
        }
      } // iProc
    }
  }
  
  // set in-transit flag to false
  V->HaloInTransit = xfe_False;
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_HaloReverseExchangeVectorBegin
int 
xf_HaloReverseExchangeVectorBegin( xf_Vector *V){
  
  int ierr, i, j, k, r, pos, nrecv_tot;
  int myRank, nProc, iProc;
  int nrecv, nsend, recvsize, sendsize, reqsize;
  int *SendElem;
  enum xfe_Bool HaloFlag;
  char *request;
  void *sbuf, *rbuf;
  xf_GenArray *ga;
  xf_ArrayParallelInfo *PInfo;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc == 1 || V->ParallelFlag == xfe_False) return xf_OK; // not parallel, so return immediately
  
  if (V->HaloInTransit){
    xf_printf("Error.  Attempting to halo-exchange a vector while halo is in transit.\n");
    return xf_Error(xf_CODE_LOGIC_ERROR);
  }
  
  // only handle glob elem linkage for now
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  
  reqsize = xf_MPI_RequestSize();
  
  for (i=0; i<V->nArray; i++){
    ga = V->GenArray + i;
    r = ga->r;
    if ((PInfo = ga->ParallelInfo) == NULL) return xf_Error(xf_PARALLEL_ERROR);
    
    HaloFlag = PInfo->HaloFlag;
    
    if (!HaloFlag){ // Normal groups receive data
      for (iProc=0; iProc<nProc; iProc++){
        nsend = PInfo->nSendElem[iProc];
        if (nsend > 0){
          SendElem = PInfo->SendElem[iProc];
          if (ga->Size == xfe_SizeInt){
            // size of receive buffer
            for (j=0, pos=0; j<nsend; j++){
              r = (ga->vr != NULL) ? ga->vr[SendElem[j]] : ga->r;
              pos += r;
            }
            sbuf = (void *) PInfo->iSendBuf[iProc];
            sendsize = pos*sizeof(int);
          }
          else{
            // pack send buffer
            for (j=0, pos=0; j<nsend; j++){
              r = (ga->vr != NULL) ? ga->vr[SendElem[j]] : ga->r;
              pos += r;
            }
            sbuf = (void *) PInfo->rSendBuf[iProc];
            sendsize = pos*sizeof(real);
          }
          request = ((char *) PInfo->Request) + iProc*reqsize;
          ierr = xf_Error(xf_MPI_Irecv(sbuf, sendsize, iProc, 0, (void *) request));
          if (ierr != xf_OK) return ierr;
        }
      } // iProc
    }
    else{ // Halo groups send data
      nrecv_tot = 0;
      for (iProc=0; iProc<nProc; iProc++){
        nrecv = PInfo->nRecvElem[iProc];
        if (nrecv > 0){
          pos = nrecv*r;
          if (ga->vr != NULL) for (j=nrecv_tot,pos=0;j<nrecv_tot+nrecv;j++) pos += ga->vr[j];
          if (ga->Size == xfe_SizeInt){
            rbuf = (void *) ga->iValue[0 + nrecv_tot];
            recvsize = pos*sizeof(int);
          }
          else{
            rbuf = (void *) ga->rValue[0 + nrecv_tot];
            recvsize = pos*sizeof(real);
          }
          request = ((char *) PInfo->Request) + iProc*reqsize;
          ierr = xf_Error(xf_MPI_Isend(rbuf, recvsize, iProc, 0, (void *) request));
          if (ierr != xf_OK) return ierr;
          nrecv_tot += nrecv;
        }
      } // iProc
    }
    
  } // i
  
  // set in-transit flag
  V->HaloInTransit = xfe_True;
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_HaloReverseExchangeVectorEnd
int 
xf_HaloReverseExchangeVectorEnd( xf_Vector *V){
  
  int ierr, i, j , k;
  int myRank, nProc, iProc;
  int nrecv, nsend, reqsize;
  int pos, r;
  int *SendElem;
  enum xfe_Bool HaloFlag;
  char *request;
  xf_GenArray *ga;
  xf_ArrayParallelInfo *PInfo;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc == 1 || !V->ParallelFlag) return xf_OK; // not parallel, so return immediately
  
  if (!V->HaloInTransit) return xf_OK; // to allow for multiple calls 
  
  // only handle glob elem linkage for now
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  
  reqsize = xf_MPI_RequestSize();
  
  for (i=0; i<V->nArray; i++){
    ga = V->GenArray + i;
    r = ga->r;
    if ((PInfo = ga->ParallelInfo) == NULL) return xf_Error(xf_PARALLEL_ERROR);
    
    HaloFlag = PInfo->HaloFlag;
    
    if (!HaloFlag){ // non-halo groups finish receiving
      for (iProc=0; iProc<nProc; iProc++){
        nsend = PInfo->nSendElem[iProc];
        if (nsend > 0){
          request = ((char *) PInfo->Request) + iProc*reqsize;
          ierr = xf_Error(xf_MPI_Wait((void *) request)); 
          if (ierr != xf_OK) return ierr;

          // distribute data (via addition) to correct locations in regular group
          SendElem = PInfo->SendElem[iProc];
          if (ga->Size == xfe_SizeInt){
            // size of receive buffer
            for (j=0, pos=0; j<nsend; j++){
              r = (ga->vr != NULL) ? ga->vr[SendElem[j]] : ga->r;
              for (k=0; k<r; k++)
                ga->iValue[SendElem[j]][k] += PInfo->iSendBuf[iProc][pos+k];
              pos += r;
            }
          }
          else{
            // pack send buffer
            for (j=0, pos=0; j<nsend; j++){
              r = (ga->vr != NULL) ? ga->vr[SendElem[j]] : ga->r;
              for (k=0; k<r; k++)
                ga->rValue[SendElem[j]][k] += PInfo->rSendBuf[iProc][pos+k];
              pos += r;
            }
          }
        }
      } // iProc
    }
    else{  // Halo groups finish sending
      for (iProc=0; iProc<nProc; iProc++){
        nrecv = PInfo->nRecvElem[iProc];
        if (nrecv > 0){
          request = ((char *) PInfo->Request) + iProc*reqsize;
          ierr = xf_Error(xf_MPI_Wait((void *) request)); 
          if (ierr != xf_OK) return ierr;
        }
      } // iProc
    }
  }
  
  // set in-transit flag to false
  V->HaloInTransit = xfe_False;
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_BcastVectorBasicInfo
int 
xf_BcastVectorBasicInfo(xf_Vector *V_Glob, xf_Vector *V, 
                        enum xfe_Bool DestroyGlob)
{
  int ierr, myRank, nProc;
  int ibuf[13], len;
  enum xfe_Bool StateNameFlag, OutputNameFlag, BasisFlag, OrderFlag;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // bcast basic vector info
  if (myRank == 0){
    ibuf[ 0] =  V_Glob->Linkage;
    ibuf[ 1] =  V_Glob->LinkageIndex;
    ibuf[ 3] =  V_Glob->SolverRole;
    ibuf[ 4] =  V_Glob->StateRank;
    ibuf[ 5] =  V_Glob->TimeIndex;
    ibuf[ 6] =  V_Glob->MGIndex;
    ibuf[ 7] =  V_Glob->nArray;
    ibuf[ 8] = (V_Glob->StateName != NULL);
    ibuf[ 9] = (V_Glob->OutputName != NULL);
    ibuf[10] = (V_Glob->Basis  != NULL);
    ibuf[11] = (V_Glob->Order  != NULL);
    ibuf[12] = V_Glob->Size;
  }
  else {
    xf_InitVector(V);
  }
  
  
  ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 13*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  V->Linkage      = ibuf[0];
  V->LinkageIndex = ibuf[1];
  V->SolverRole   = ibuf[3];
  V->StateRank    = ibuf[4];
  V->TimeIndex    = ibuf[5];
  V->MGIndex      = ibuf[6];
  V->nArray       = ibuf[7];
  StateNameFlag   = ibuf[8];
  OutputNameFlag  = ibuf[9];
  BasisFlag       = ibuf[10];
  OrderFlag       = ibuf[11];
  V->Size         = ibuf[12];
  
  
  if (StateNameFlag){
    if (myRank == 0){ // this is destructive to V_Glob
      V->StateName  = V_Glob->StateName;  
      if (DestroyGlob)
        V_Glob->StateName  = NULL;
    }
    else{
      ierr = xf_Error(xf_Alloc2((void ***) &V->StateName, V->StateRank, xf_MAXSTRLEN, sizeof(char)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_MPI_Bcast(V->StateName[0], V->StateRank*xf_MAXSTRLEN*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  if (OutputNameFlag){
    if (myRank == 0){  // this is destructive to V_Glob
      V->OutputName   = V_Glob->OutputName;  
      if (DestroyGlob)
        V_Glob->OutputName = NULL;
      
      len = strlen(V->OutputName)+1;
    }
    ierr = xf_Error(xf_MPI_Bcast((void *) &len, 1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    if (myRank > 0){
      ierr = xf_Error(xf_Alloc((void **) &V->OutputName, len, sizeof(char)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_MPI_Bcast(V->OutputName, len*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  if (BasisFlag){
    if (myRank == 0){  // this is destructive to V_Glob
      V->Basis        = V_Glob->Basis;   
      if (DestroyGlob)
        V_Glob->Basis = NULL;
    }
    else{
      ierr = xf_Error(xf_Alloc((void **) &V->Basis, V->nArray, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_MPI_Bcast(V->Basis, V->nArray*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  if (OrderFlag){
    if (myRank == 0){  // this is destructive to V_Glob
      V->Order        = V_Glob->Order; 
      if (DestroyGlob)
        V_Glob->Order = NULL;
    }
    else{
      ierr = xf_Error(xf_Alloc((void **) &V->Order, V->nArray, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_MPI_Bcast(V->Order, V->nArray*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BcastVector
static int 
xf_BcastVector(xf_Vector *V_proc0, xf_Vector *V)
{
  int ierr, i, j, myRank, nProc, ntot, ibuf[4];
  enum xfe_Bool VariableOrder, nCompFlag, vOrderFlag;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  //broadcast basic info
  ierr = xf_Error(xf_BcastVectorBasicInfo(V_proc0, V, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0) {
    xf_InitVector(V);
    ierr = xf_Error(xf_DestroyVector(V, xfe_True));
    if (ierr != xf_OK) return ierr;
    //point to V_proc0
    V = V_proc0;
  }
  else {
    ierr = xf_Error(xf_Alloc((void **)&V->GenArray, V->nArray, 
                             sizeof(xf_GenArray)));
    if (ierr != xf_OK) return ierr;          
  }
  //processor specific info
  if (myRank == 0){
    ibuf[0] = (V->nComp  != NULL);
    ibuf[1] = (V->vOrder != NULL);
  }
  ierr = xf_Error(xf_MPI_Bcast((void *)ibuf, 2*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  nCompFlag       = ibuf[0];
  vOrderFlag      = ibuf[1];
  
  if (nCompFlag){
    if (myRank != 0) {
      //nComp 
      ierr = xf_Error(xf_Alloc((void **) &V->nComp, V->nArray, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_MPI_Bcast(V->nComp, V->nArray*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
  }
  if (vOrderFlag){
    if (myRank != 0) {
      ierr = xf_Error(xf_VAlloc2((void ***) &V->vOrder, V->nArray, V->nComp, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    for (i=0,ntot=0; i<V->nArray; i++) ntot += V->nComp[i];
    ierr = xf_Error(xf_MPI_Bcast((void *)V->vOrder[0], ntot*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
  }
  
  //broadcast GenArrays
  for (i = 0; i < V->nArray; i++){
    if (myRank == 0){
      ibuf[0] = V->GenArray[i].Size;
      ibuf[1] = V->GenArray[i].n;
      ibuf[2] = V->GenArray[i].r;
      ibuf[3] = (V->GenArray[i].vr != NULL); // variable order flag
    }
    
    ierr = xf_Error(xf_MPI_Bcast((void *)ibuf, 4*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    //unwrap, allocate and bcast
    V->GenArray[i].Size = ibuf[0];
    V->GenArray[i].n    = ibuf[1];
    V->GenArray[i].r    = ibuf[2];
    VariableOrder       = ibuf[3];
    
    if (VariableOrder){    // take care of vr != NULL
      if (myRank != 0){
        ierr = xf_Error(xf_Alloc((void **) &V->GenArray[i].vr, V->GenArray[i].n, sizeof(int)));
        if (ierr != xf_OK) return ierr;          
      }
      ierr = xf_Error(xf_MPI_Bcast((void *) V->GenArray[i].vr, V->GenArray[i].n*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
    }
    else{
    	V->GenArray[i].vr = NULL;
    }
    
    if (myRank != 0){ // allocate on non-root processors
      if (V->GenArray[i].Size == xfe_SizeInt){
        if (VariableOrder){
          ierr = xf_Error(xf_VAlloc2((void ***)&V->GenArray[i].iValue, 
                                     V->GenArray[i].n, V->GenArray[i].vr,
                                     sizeof(int)));
        }
        else{
          ierr = xf_Error(xf_Alloc2((void ***)&V->GenArray[i].iValue, 
                                    V->GenArray[i].n, V->GenArray[i].r, 
                                    sizeof(int)));
        }
        if (ierr != xf_OK) return ierr;
        V->GenArray[i].rValue = NULL;
      }
      else{
        if (VariableOrder){
          ierr = xf_Error(xf_VAlloc2((void ***)&V->GenArray[i].rValue, 
                                     V->GenArray[i].n, V->GenArray[i].vr,
                                     sizeof(real)));
        }
        else{
          ierr = xf_Error(xf_Alloc2((void ***)&V->GenArray[i].rValue, 
                                    V->GenArray[i].n, V->GenArray[i].r, 
                                    sizeof(real)));
        }
        if (ierr != xf_OK) return ierr;
        V->GenArray[i].iValue = NULL;
      }
      V->GenArray[i].ParallelInfo = NULL;
    }
    // let's broadcast
    ntot = V->GenArray[i].n*V->GenArray[i].r;
    if (VariableOrder) for (j=0,ntot=0; j<V->GenArray[i].n; j++) ntot += V->GenArray[i].vr[j];
    if (V->GenArray[i].Size == xfe_SizeInt){
      ierr = xf_Error(xf_MPI_Bcast((void *)V->GenArray[i].iValue[0], ntot*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
    }
    else {
      ierr = xf_Error(xf_MPI_Bcast((void *)V->GenArray[i].rValue[0], ntot*sizeof(real), 0));
      if (ierr != xf_OK) return ierr;
    }
    
  }//nArray  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeVector
int 
xf_ParallelizeVector( xf_Mesh *Mesh, xf_Vector *V_Glob, xf_Vector *V)
{
  int ierr, len, j, k, r;
  int myRank, nProc;
  int negrp, egrp, nelem, size;
  int nelemtot;
  int ibuf[12];
  int *sindex, *nElem, **ElemList;
  enum xfe_Bool VariableOrder = xfe_False;
  void *rbuf, *sbuf;
  int  *rlen, *slen;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc == 1){ 
    /* This serial behavior is useful.  V inherits all of V_Glob's
     data, and V_Glob has its data initialized to null. */
    (*V) = (*V_Glob);
    xf_InitVector(V_Glob);
    return xf_OK;
  }
  
  if (Mesh == NULL) return xf_Error(xf_PARALLEL_ERROR);
  
  // bcast basic vector info
  ierr = xf_Error(xf_BcastVectorBasicInfo(V_Glob, V, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  // handle each linkage type individually
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  
  // bcast whether we are dealing with a variable order vector
  if (myRank == 0)
    VariableOrder = ((V_Glob->nComp != NULL) && (V_Glob->vOrder != NULL));
  ierr = xf_Error(xf_MPI_Bcast((void *)&VariableOrder, sizeof(enum xfe_Bool), 0));
  if (ierr != xf_OK) return ierr;
  
  negrp = Mesh->nElemGroup;
  if (V->nArray != negrp) return xf_Error(xf_PARALLEL_ERROR);
  
  V->nArraySelf = negrp;
  
  // allocate GenArray for V
  V->nArray = 2*negrp;
  ierr = xf_Error(xf_Alloc((void **) &V->GenArray, V->nArray, 
                           sizeof(xf_GenArray)));
  if (ierr != xf_OK) return ierr;
  
  if (VariableOrder){
    // Allocate and fill in nComp 
    ierr = xf_Error(xf_Alloc((void **) &V->nComp, V->nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (egrp=0; egrp<V->nArray; egrp++) V->nComp[egrp] = Mesh->ElemGroup[egrp].nElem;
    // Allocate vOrder
    ierr = xf_Error(xf_VAlloc2((void ***) &V->vOrder, V->nArray, V->nComp, sizeof(int)));
    if (ierr != xf_OK) return ierr;    
  }  
  
  
  // loop over element groups (halos too) and prepare to allocate each GenArray[egrp]
  for (egrp=0; egrp<V->nArray; egrp++){
    
    // Initialize GenArray
    xf_InitGenArray(V->GenArray+egrp);
    
    if (myRank == 0){
      ibuf[0] = V_Glob->GenArray[egrp%negrp].Size;
      ibuf[1] = V_Glob->GenArray[egrp%negrp].r;
    }
    
    // bcast Size and r
    ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 2*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    V->GenArray[egrp].Size = ibuf[0];
    V->GenArray[egrp].n    = Mesh->ElemGroup[egrp].nElem;
    V->GenArray[egrp].r    = ibuf[1];
    
    if (VariableOrder){
      ierr = xf_Error(xf_Alloc((void **) &V->GenArray[egrp].vr, Mesh->ElemGroup[egrp].nElem, 
                               sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    else V->GenArray[egrp].vr = NULL;
    
  } // egrp
  
  
  // Set Basis and Order on halo groups
  if (V->Basis != NULL){
    ierr = xf_Error(xf_ReAlloc((void **) &V->Basis, V->nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (egrp=0; egrp<negrp; egrp++) V->Basis[negrp+egrp] = V->Basis[egrp];
  }
  
  if (V->Order != NULL){
    ierr = xf_Error(xf_ReAlloc((void **) &V->Order, V->nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (egrp=0; egrp<negrp; egrp++)
      V->Order[negrp+egrp] = V->Order[egrp];
  }
  
  // set pointers to NULL
  nElem    = NULL;
  ElemList = NULL;
  
  // allocate memory on proc 0
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &nElem, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }
  
  // loop over element groups, pack + scatter data (halos included in loop)
  for (egrp=0; egrp<2*negrp; egrp++){
    
    if ((nelem = Mesh->ElemGroup[egrp].nElem) != V->GenArray[egrp].n)
      return xf_Error(xf_PARALLEL_ERROR);
    
    // initialize variable order values
    if (VariableOrder){
      V->Order[egrp] = 0;
      V->GenArray[egrp].r = 0;
    }
    
    // gather number of elems to inform proc 0 of # on each proc
    ierr = xf_Error(xf_MPI_Gather((void *) &nelem,(void *) nElem, 1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // calculate total number of elements in this group
    nelemtot = nelem;
    ierr = xf_Error(xf_MPI_Allreduce(&nelemtot, 1, xfe_SizeInt, xfe_MPI_SUM));
    if (ierr != xf_OK) return ierr;
    if (nelemtot == 0) continue; // nothing to do if no elements!
    
    // proc 0 re-allocates ElemList
    if (myRank == 0){
      ierr = xf_Error(xf_VReAlloc2( (void ***) &ElemList, nProc, nElem, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    
    // variable gather global element numbers onto proc 0
    rbuf = ((myRank == 0) ? (void *) ElemList[0] : NULL);
    ierr = xf_Error(xf_MPI_Gatherv((void *) Mesh->ParallelInfo->ElemLoc2Glob[egrp], 
                                   nelem, rbuf, nElem, sizeof(int), 0));
    if (ierr!=xf_OK) return ierr;
    
    if (VariableOrder){
      sindex = ((myRank == 0) ? ElemList[0] : NULL);
      // send list of variable orders to each processor
      sbuf = ((myRank == 0) ? (void *) V_Glob->vOrder[egrp%negrp] : NULL);
      rbuf = ((nelem > 0) ? (void *) V->vOrder[egrp] : NULL);
      ierr = xf_Error(xf_MPI_PScatterv(sbuf, sindex, nElem, rbuf, nelem, sizeof(int), 0));
      if (ierr!=xf_OK) return ierr;
      // send list of variable array ranks to each processor
      sbuf = ((myRank == 0) ? (void *) V_Glob->GenArray[egrp%negrp].vr : NULL);
      rbuf = ((nelem > 0) ? (void *) V->GenArray[egrp].vr : NULL);
      ierr = xf_Error(xf_MPI_PScatterv(sbuf, sindex, nElem, rbuf, nelem, sizeof(int), 0));
      if (ierr!=xf_OK) return ierr;
      
      // set V->Order to max of V->vOrder
      V->Order[egrp] = 0;
      for (j=0; j<V->nComp[egrp]; j++)
        V->Order[egrp] = max(V->Order[egrp], V->vOrder[egrp][j]);
      // set V->GenArray[egrp].r to max of V->GenArray[egrp].vr
      V->GenArray[egrp].r = 0;
      for (j=0; j<V->GenArray[egrp].n; j++)
        V->GenArray[egrp].r = max(V->GenArray[egrp].r, V->GenArray[egrp].vr[j]);
      // sanity check
      if (V->GenArray[egrp].n != V->nComp[egrp]) return xf_Error(xf_CODE_LOGIC_ERROR);
    }
    
    // Allocate general array on each processor (halos included in this loop)
    ierr = xf_Error(xf_AllocGenArray(V->GenArray+egrp));
    if (ierr != xf_OK) return ierr;
    
    if (egrp >= negrp){ // if Halo
      // zero out Halo
      for (j=0; j<V->GenArray[egrp].n; j++){
        r = ((V->GenArray[egrp].vr == NULL) ? V->GenArray[egrp].r : V->GenArray[egrp].vr[j]);
        if (V->GenArray[egrp].iValue != NULL)
          for (k=0; k<r; k++) V->GenArray[egrp].iValue[j][k] = 0;
        if (V->GenArray[egrp].rValue != NULL)
          for (k=0; k<r; k++) V->GenArray[egrp].rValue[j][k] = 0.;
      }
      continue;   // do not send data to halos
    }
    
    // pack + scatter data
    sindex = ((myRank == 0) ? ElemList[0] : NULL);
    if (V->GenArray[egrp].Size == xfe_SizeInt){
      sbuf = ((myRank == 0) ? (void *) V_Glob->GenArray[egrp%negrp].iValue[0] : NULL);
      rbuf = ((nelem > 0) ? (void *) V->GenArray[egrp].iValue[0] : NULL);
      size = sizeof(int);
    }
    else if (V->GenArray[egrp].Size == xfe_SizeReal){
      sbuf = ((myRank == 0) ? (void *) V_Glob->GenArray[egrp%negrp].rValue[0] : NULL);
      rbuf = ((nelem > 0) ? (void *) V->GenArray[egrp].rValue[0]: NULL);
      size = sizeof(real);
    }
    else return xf_Error(xf_NOT_SUPPORTED);
    
    if (VariableOrder){ // use double pack+scatter for variable orders
      slen   = ((myRank == 0) ? V_Glob->GenArray[egrp%negrp].vr : NULL);
      rlen   = ((nelem > 0) ? V->GenArray[egrp].vr : NULL);
      
      ierr = xf_Error(xf_MPI_DPScatterv(sbuf, slen, sindex, nElem, rbuf, 
                                        rlen, nelem, size, 0));
      
    }
    else{ // standard vector pack+scatter for constant orders
      ierr = xf_Error(xf_MPI_PScatterv(sbuf, sindex, nElem, rbuf, nelem, 
                                       V->GenArray[egrp].r*size, 0));
      if (ierr!=xf_OK) return ierr;
    }
    
  } // egrp
  
  // release memory
  xf_Release( (void  *) nElem);
  xf_Release2((void **) ElemList);  
  
  // parallel-prep V (ParallelFlag is set here)
  ierr = xf_Error(xf_ParallelPrepVector(Mesh, V));
  if (ierr != xf_OK) return ierr;
  
  /* Fill in data on halos */
  // begin communication of halo data
  ierr = xf_Error(xf_HaloExchangeVectorBegin(V));
  if (ierr != xf_OK) return ierr;
  // end communication of halo data
  ierr = xf_Error(xf_HaloExchangeVectorEnd(V));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeVectorSet
static int 
xf_ParallelizeVectorSet( xf_Mesh *Mesh, xf_VectorSet *VS_Glob, 
                        xf_VectorSet *VS)
{
  int ierr, i;
  int myRank, nProc;
  xf_Vector *V_Glob;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<VS->nVector; i++){
    V_Glob = ((myRank == 0) ? VS_Glob->Vector + i : NULL);
    ierr = xf_Error(xf_ParallelizeVector(Mesh, V_Glob, VS->Vector + i));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeData
static int 
xf_ParallelizeData( xf_All *All, xf_Data *D_Glob, xf_Data *D)
{
  int ierr, len, nVector;
  int myRank, nProc;
  xf_Mesh *Mesh;
  xf_Vector *V_Glob, *V;
  xf_VectorSet *VS_Glob, *VS;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // pull off pointer to Mesh
  Mesh = All->Mesh;
  
  if (myRank == 0){
    // set Title
    len = strlen(D_Glob->Title)+1;
    ierr = xf_Error(xf_Alloc((void **) &D->Title, len, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    strncpy(D->Title, D_Glob->Title, len);
    
    // set Type
    D->Type = D_Glob->Type;
  }
  
  // bcast len of Title
  ierr = xf_Error(xf_MPI_Bcast((void *) &len, 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // alloc Title on all procs > 0
  if (myRank > 0){
    ierr = xf_Error(xf_Alloc((void **) &D->Title, len, sizeof(char)));
    if (ierr != xf_OK) return ierr;
  }
  
  // bcast Title
  ierr = xf_Error(xf_MPI_Bcast((void *) D->Title, len*sizeof(char), 0));
  if (ierr != xf_OK) return ierr;
  
  // bcast Type
  ierr = xf_Error(xf_MPI_Bcast((void *) &D->Type, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  
  // call appropriate function
  switch (D->Type){
    case xfe_Vector:
      ierr = xf_Error(xf_CreateVector(&V));
      if (ierr != xf_OK) return ierr;
      
      V_Glob = ((myRank == 0) ? (xf_Vector *) D_Glob->Data : NULL);
      ierr = xf_Error(xf_ParallelizeVector(Mesh, V_Glob, V));
      if (ierr != xf_OK) return ierr;
      
      D->Data = V;
      break;
      
    case xfe_VectorSet:
      
      VS_Glob = ((myRank == 0) ? (xf_VectorSet *) D_Glob->Data : NULL);
      if (myRank == 0) nVector = VS_Glob->nVector;
      ierr = xf_Error(xf_MPI_Bcast((void *) &nVector, sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_CreateVectorSet(nVector, &VS));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ParallelizeVectorSet(Mesh, VS_Glob, VS));
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
//   FUNCTION Definition: xf_ParallelizeDataSet
int 
xf_ParallelizeDataSet( xf_All *All, xf_DataSet *DataSet_Glob, 
                      xf_DataSet *DataSet){
  
  int ierr, count, i;
  int myRank, nProc;
  xf_Data *D, *P, *D_Glob;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // initialize head and tail on DataSet
  DataSet->Head = DataSet->Tail = NULL;
  
  // proc 0 counts # of write-able data nodes
  if (myRank == 0){
    D_Glob = DataSet_Glob->Head;
    count = 0;
    while (D_Glob != NULL){
      if (D_Glob->ReadWrite) count++;
      D_Glob = D_Glob->Next;
    }
  }
  
  // bcast count to all procs
  ierr = xf_Error(xf_MPI_Bcast((void *) &count, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0) 
    D_Glob = DataSet_Glob->Head;
  else
    D_Glob = NULL;
  P = NULL;
  for (i=0; i<count; i++){
    
    // find i'th read-writeable data node in DataSet_Glob
    if (myRank == 0){
      while (!D_Glob->ReadWrite){
        D_Glob = D_Glob->Next;
        if (D_Glob == NULL) return xf_Error(xf_PARALLEL_ERROR);
      }
      
      xf_printf("Parallelizing %s\n", D_Glob->Title);
    }
    
    // create data node
    ierr = xf_Error(xf_CreateData(&D));
    if (ierr != xf_OK) return ierr;
    
    // parallelize D_Glob into D
    ierr = xf_Error(xf_ParallelizeData(All, D_Glob, D));
    if (ierr != xf_OK) return ierr;
    D->ReadWrite = xfe_True;
    
    if (P == NULL) DataSet->Head = D;
    else P->Next = D;
    D->Prev = P;
    P = D;
    if (myRank == 0) D_Glob = D_Glob->Next;
    
  }
  DataSet->Tail = P;
  
  return xf_OK;
}




/*-------------------*/
/* UnParallelization */
/*-------------------*/  

/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeVector
int 
xf_UnParallelizeVector( xf_Mesh *Mesh, xf_Vector *V, xf_Vector *V_Glob)
{
  int ierr, i, j;
  int myRank, nProc;
  int negrp, egrp, nelem, size;
  int ibuf[10];
  int *rindex, *nElem, **ElemList;
  enum xfe_Bool VariableOrder = xfe_False;
  void *rbuf, *sbuf;
  int  *rlen;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // should not call this in serial
  if (nProc == 1) return xf_Error(xf_PARALLEL_ERROR);
  
  // need mesh
  if (Mesh == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // variable order flag
  
  // make sure basic info agrees on all procs
  VariableOrder = ((V->nComp != NULL) && (V->vOrder != NULL));
  if (myRank == 0){
    ibuf[0] = V->Linkage;
    ibuf[1] = V->LinkageIndex;
    ibuf[3] = V->SolverRole;
    ibuf[4] = V->StateRank;
    ibuf[5] = V->TimeIndex;
    ibuf[6] = V->MGIndex;
    ibuf[7] = V->nArray;
    ibuf[8] = V->Size;
    ibuf[9] = VariableOrder;
  }
  
  ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 10*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  if (ibuf[0] != V->Linkage)      return xf_Error(xf_PARALLEL_ERROR);
  if (ibuf[1] != V->LinkageIndex) return xf_Error(xf_PARALLEL_ERROR);
  if (ibuf[3] != V->SolverRole)   return xf_Error(xf_PARALLEL_ERROR);
  if (ibuf[4] != V->StateRank)    return xf_Error(xf_PARALLEL_ERROR);
  if (ibuf[5] != V->TimeIndex)    return xf_Error(xf_PARALLEL_ERROR);
  if (ibuf[6] != V->MGIndex)      return xf_Error(xf_PARALLEL_ERROR);
  if (ibuf[7] != V->nArray)       return xf_Error(xf_PARALLEL_ERROR);
  if (ibuf[8] != V->Size)         return xf_Error(xf_PARALLEL_ERROR);
  if (ibuf[9] != VariableOrder)   return xf_Error(xf_PARALLEL_ERROR);
  
  // set V_Glob basic info
  if (myRank == 0){
    V_Glob->Linkage      = V->Linkage;
    V_Glob->LinkageIndex = V->LinkageIndex;
    V_Glob->SolverRole   = V->SolverRole;
    V_Glob->StateRank    = V->StateRank;
    
    if (V->StateName != NULL){
      ierr = xf_Error(xf_Alloc2((void ***) &V_Glob->StateName, V_Glob->StateRank,
                                xf_MAXSTRLEN, sizeof(char)));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<V_Glob->StateRank; i++)
        strcpy(V_Glob->StateName[i], V->StateName[i]);
    }
    else V_Glob->StateName = NULL;
    
    ierr = xf_Error(xf_AllocString(&V_Glob->OutputName, xf_MAXSTRLEN, V->OutputName));
    if (ierr != xf_OK) return ierr;
    
    V_Glob->TimeIndex    = V->TimeIndex; 
    V_Glob->MGIndex      = V->MGIndex; 
    V_Glob->Size         = V->Size; 
  }
  
  // handle each linkage type individually
  if (V->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  
  negrp = Mesh->nElemGroup;
  if (V->nArray/2 != negrp) return xf_Error(xf_PARALLEL_ERROR);
  
  // Basis, Order, nComp
  if (myRank == 0){
    if (V->Basis != NULL){
      ierr = xf_Error(xf_Alloc( (void **) &V_Glob->Basis, negrp, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<negrp; i++) V_Glob->Basis[i] = V->Basis[i];
    }
    else V->Basis = NULL;
    if (V->Order != NULL){
      ierr = xf_Error(xf_Alloc( (void **) &V_Glob->Order, negrp, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      for (i=0; i<negrp; i++) V_Glob->Order[i] = V->Order[i];
    }
    else V->Order = NULL;
    if (VariableOrder){
      ierr = xf_Error(xf_Alloc( (void **) &V_Glob->nComp, negrp, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // Allocate space for GenArray on V_Glob
  if (myRank == 0){
    V_Glob->nArray = negrp;
    ierr = xf_Error(xf_Alloc((void **) &V_Glob->GenArray, V_Glob->nArray, 
                             sizeof(xf_GenArray)));
    if (ierr != xf_OK) return ierr;
  }
  
  // set basic info about GenArray, including number of total elements
  for (egrp=0; egrp<negrp; egrp++){
    nelem = Mesh->ElemGroup[egrp].nElem;
    ierr = xf_Error(xf_MPI_Allreduce( (void *) &nelem, 1, xfe_SizeInt, xfe_MPI_SUM));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0){
      V_Glob->GenArray[egrp].Size = V->GenArray[egrp].Size;
      V_Glob->GenArray[egrp].r    = V->GenArray[egrp].r;
      V_Glob->GenArray[egrp].n    = nelem;
      V_Glob->GenArray[egrp].vr   = NULL;
      
      if (VariableOrder){
        // .vr
        ierr = xf_Error(xf_ReAlloc((void **) &V_Glob->GenArray[egrp].vr, nelem, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        // fill in nComp
        V_Glob->nComp[egrp] = nelem;
      }
      
    }
  } // egrp
  
  // .vOrder on V_Glob
  if ((myRank == 0) && (VariableOrder)){
    ierr = xf_Error(xf_VAlloc2((void ***) &V_Glob->vOrder, negrp, V_Glob->nComp, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
	
  
  
  // set pointers to NULL
  nElem    = NULL;
  ElemList = NULL;
  
  // allocate memory on proc 0
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &nElem, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }
  
  // loop over element groups, gather + unpack data
  for (egrp=0; egrp<negrp; egrp++){
    
    if ((nelem = Mesh->ElemGroup[egrp].nElem) != V->GenArray[egrp].n)
      return xf_Error(xf_PARALLEL_ERROR);
    
    // gather number of elems to inform proc 0 of # on each proc
    ierr = xf_Error(xf_MPI_Gather((void *) &nelem,(void *) nElem, 1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // proc 0 re-allocates ElemList
    if (myRank == 0){
      ierr = xf_Error(xf_VReAlloc2( (void ***) &ElemList, nProc, nElem, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    
    // variable gather global element numbers onto proc 0
    rbuf = ((myRank == 0) ? (void *) ElemList[0] : NULL);
    ierr = xf_Error(xf_MPI_Gatherv((void *) Mesh->ParallelInfo->ElemLoc2Glob[egrp], 
                                   nelem, rbuf, nElem, sizeof(int), 0));
    if (ierr!=xf_OK) return ierr;
    
    if (VariableOrder){
      rindex = ((myRank == 0) ? ElemList[0] : NULL);
      // root receives list of variable orders from each processor
      sbuf = ((nelem > 0) ? (void *) V->vOrder[egrp] : NULL);
      rbuf = ((myRank == 0) ? (void *) V_Glob->vOrder[egrp] : NULL);
      ierr = xf_Error(xf_MPI_PGatherv(sbuf, nelem, rbuf, rindex, nElem,  sizeof(int), 0));
      if (ierr!=xf_OK) return ierr;
      // root receives list of variable array ranks from each processor
      sbuf = ((nelem > 0) ? (void *) V->GenArray[egrp].vr : NULL);
      rbuf = ((myRank == 0) ? (void *) V_Glob->GenArray[egrp].vr : NULL);
      ierr = xf_Error(xf_MPI_PGatherv(sbuf, nelem, rbuf, rindex, nElem,  sizeof(int), 0));
      if (ierr!=xf_OK) return ierr;
      
      if (myRank == 0){
        // set V_Glob->Order to max of V_Glob->vOrder
        V_Glob->Order[egrp]=0;
        for (j=0; j<V_Glob->nComp[egrp]; j++)
          V_Glob->Order[egrp] = max(V_Glob->Order[egrp], V_Glob->vOrder[egrp][j]);
        // set V_Glob->GenArray[egrp].r to max of V_Glob->GenArray[egrp].vr
        V_Glob->GenArray[egrp].r = 0;
        for (j=0; j<V_Glob->GenArray[egrp].n; j++)
          V_Glob->GenArray[egrp].r = max(V_Glob->GenArray[egrp].r, V_Glob->GenArray[egrp].vr[j]);
      }
    }
    
    // Allocate general array on root (now that have V_Glob->GenaArray[egrp].vr if variable order)
    if (myRank == 0){
      ierr = xf_Error(xf_AllocGenArray(V_Glob->GenArray+egrp));
      if (ierr != xf_OK) return ierr;
    }
    
    
    // gather and unpack data
    rindex = ((myRank == 0) ? ElemList[0] : NULL);
    if (V->GenArray[egrp].Size == xfe_SizeInt){
      sbuf = ((V->GenArray[egrp].n == 0) ? NULL : (void *) V->GenArray[egrp].iValue[0]);
      rbuf = ((myRank == 0) ? (void *) V_Glob->GenArray[egrp].iValue[0] : NULL);
      size = sizeof(int);
    }
    else if (V->GenArray[egrp].Size == xfe_SizeReal){
      sbuf = ((V->GenArray[egrp].n == 0) ? NULL : (void *) V->GenArray[egrp].rValue[0]);
      rbuf = ((myRank == 0) ? (void *) V_Glob->GenArray[egrp].rValue[0] : NULL);
      size = sizeof(real);
    }
    else return xf_Error(xf_NOT_SUPPORTED);
    
    if (VariableOrder){ // double pack+gather for variable orders
      rlen   = ((myRank == 0) ? V_Glob->GenArray[egrp].vr : NULL);
      ierr = xf_Error(xf_MPI_DPGatherv(sbuf, V->GenArray[egrp].vr, nelem, rbuf, rlen, 
                                       rindex, nElem, size, 0));
      if (ierr!=xf_OK) return ierr;
    }
    else{ // standard vector pack+gather for constant orders
      ierr = xf_Error(xf_MPI_PGatherv(sbuf, nelem, rbuf, rindex, nElem, 
                                      V->GenArray[egrp].r*size, 0));
      if (ierr!=xf_OK) return ierr;
    }
    
  } // egrp
  
  // release memory
  xf_Release( (void  *) nElem);
  xf_Release2((void **) ElemList);  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeVectorSet
static int 
xf_UnParallelizeVectorSet( xf_Mesh *Mesh, xf_VectorSet *VS, 
                          xf_VectorSet *VS_Glob)
{
  int ierr, i;
  int myRank, nProc;
  xf_Vector *V_Glob; 
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<VS->nVector; i++){
    V_Glob = ((myRank == 0) ? VS_Glob->Vector + i : NULL);
    ierr = xf_Error(xf_UnParallelizeVector(Mesh, VS->Vector + i, V_Glob));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeData
static int 
xf_UnParallelizeData( xf_All *All, xf_Data *D, xf_Data *D_Glob)
{
  int ierr, len, nVector;
  int myRank, nProc;
  xf_Mesh *Mesh;
  xf_Vector *V_Glob;
  xf_VectorSet *VS, *VS_Glob;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // pull off pointer to Mesh
  Mesh = All->Mesh;
  
  if (myRank == 0){
    // set Title
    len = strlen(D->Title)+1;
    ierr = xf_Error(xf_Alloc((void **) &D_Glob->Title, len, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    strncpy(D_Glob->Title, D->Title, len);
    
    // set Type
    D_Glob->Type = D->Type;
  }
  
  // call appropriate function
  switch (D->Type){
    case xfe_Vector:
      if (myRank == 0){
        ierr = xf_Error(xf_CreateVector(&V_Glob));
        if (ierr != xf_OK) return ierr;
      }
      else V_Glob = NULL;
      
      ierr = xf_Error(xf_UnParallelizeVector(Mesh, (xf_Vector *) D->Data, V_Glob));
      if (ierr != xf_OK) return ierr;
      
      if (myRank == 0) D_Glob->Data = V_Glob;
      break;
      
    case xfe_VectorSet:
      VS = (xf_VectorSet *) D->Data;
      nVector = VS->nVector;
      if (myRank == 0){
        ierr = xf_Error(xf_CreateVectorSet(nVector, &VS_Glob));
        if (ierr != xf_OK) return ierr;
      }
      else VS_Glob = NULL;
      
      ierr = xf_Error(xf_UnParallelizeVectorSet(Mesh, VS, VS_Glob));
      if (ierr != xf_OK) return ierr;
      
      if (myRank == 0) D_Glob->Data = VS_Glob;
      break;
      
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeDataSet
int 
xf_UnParallelizeDataSet( xf_All *All, xf_DataSet *DataSet, xf_DataSet *DataSet_Glob)
{
  int ierr, len;
  int myRank, nProc;
  enum xfe_Bool done, match, rw;
  char *Title;
  char TitleBuffer[xf_MAXSTRLEN];
  xf_Data *D, *P_Glob, *D_Glob;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // initialize head and tail on DataSet_Glob
  if (myRank == 0)
    DataSet_Glob->Head = DataSet_Glob->Tail = NULL;
  
  /* We require that all procs have the same order to their
   write-able dataset lists. */  
  Title  = NULL;
  P_Glob = NULL;
  D_Glob = NULL;
  D = DataSet->Head;
  
  done = xfe_False;
  
  while (!done){
    if (myRank == 0){
      if (D == NULL)
        done = xfe_True;
      else{
        while (D->ReadWrite != xfe_True){ //make sure it is writeable
          D = D->Next;
          if (D == NULL){
            done = xfe_True;
	    break;
	  }
        }
	if(!done)
          strcpy(TitleBuffer, D->Title);
        match = xfe_True;
      }
    }
    else
      match = xfe_False;
    
    ierr = xf_Error(xf_MPI_Bcast((void *) &done, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    if (done)
      break;
    
    //make sure data D is the same on all procs
    if (myRank == 0) len = strlen(D->Title)+1;
    
    ierr = xf_Error(xf_MPI_Bcast((void *) &len, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // bcast Title on root
    ierr = xf_Error(xf_MPI_Bcast((void **)&TitleBuffer, len*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    
    // If not root proc, search for data with Title
    if (myRank != 0){
      D = DataSet->Head;
      while (D != NULL){
        if (strcmp(D->Title, TitleBuffer) == 0){
          match = xfe_True;
          break;
        }
        D = D->Next;
      }
    }
    
    // Find minimum value of match on all procs; should be 1 if data was found on all
    ierr = xf_Error(xf_MPI_Allreduce((int*)&match, 1, xfe_SizeInt, xfe_MPI_MIN));
    if (ierr != xf_OK) return ierr;
    
    rw = D->ReadWrite;//should be should be writeable in all procs 
    ierr = xf_Error(xf_MPI_Allreduce((int*)&rw, 1, xfe_SizeInt, xfe_MPI_MIN));
    if (ierr != xf_OK) return ierr;
    
    if (!rw)//
      return xf_Error(xf_PARALLEL_ERROR); 
    
    if (!match){
      return xf_Error(xf_PARALLEL_ERROR); 
    }
  
    // create a data node for Glob
    if (myRank == 0){
      ierr = xf_Error(xf_CreateData(&D_Glob));
      if (ierr != xf_OK) return ierr;
      D_Glob->ReadWrite = xfe_True;
    }
    
    // unparallelize D into D_Glob
    ierr = xf_Error(xf_UnParallelizeData(All, D, D_Glob));
    if (ierr != xf_OK) return ierr;
    
    // link D_Glob into DataSet_Glob
    if (myRank == 0){
      if (P_Glob == NULL) DataSet_Glob->Head = D_Glob;
      else P_Glob->Next = D_Glob;
      D_Glob->Prev = P_Glob;
      P_Glob = D_Glob;
    }
    D = D->Next;
  }
  
  // point tail of linked-list to last P written
  if (myRank == 0)  DataSet_Glob->Tail = P_Glob;
  
  // release memory
  xf_Release( (void *) Title);
  
  return xf_OK;
}


/*--------*/
/* Extras */
/*--------*/  

/******************************************************************/
//   FUNCTION Definition: xf_ProcViz
int
xf_ProcViz(xf_All *All)
{
  
  int ierr, egrp, elem, myRank;
  xf_Vector *ProcID;
  xf_Data *D;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  
  // obtain nProc = # domains
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindVector(All, "ProcID", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,
                                xfe_True, &D, &ProcID, NULL));
  if (ierr != xf_OK) return ierr;
  
  D->ReadWrite = xfe_True; // make data writable
  
  // set ProcID = proc # for each elem
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      ProcID->GenArray[egrp].rValue[elem][0] = (real) (myRank);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelSumNodes
int 
xf_ParallelSumNodes( xf_Mesh *Mesh, void *NodeData, int rank, enum xfe_SizeType Size)
{
  int ierr, i, off, k;
  int myRank, nProc, iProc;
  int DataSize;
  int nNode_Glob, sum;
  int nmax, *ibuf, *nNode, **NodeList;
  void *rbuf;
  void *NodeData_Glob = NULL;
  void *NodeData_Buf = NULL;
  real *rData = NULL, *rData_Glob = NULL, *rData_Buf = NULL;
  int  *iData = NULL, *iData_Glob = NULL, *iData_Buf = NULL;
  
  // return immediately if not parallel
  if (Mesh->ParallelInfo == NULL) return xf_OK;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  switch (Size){
    case xfe_SizeInt:
      DataSize = sizeof(int);
      iData = (int *) NodeData;
      break;
    case xfe_SizeReal:
      DataSize = sizeof(real);
      rData = (real *) NodeData;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  // set pointers to null
  nNode    = NULL;
  NodeList = NULL;
  
  // alloc memory on proc 0
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &nNode, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }
  
  // determine max glob node # + 1 on each proc
  nmax = 0;
  for (i=0; i<Mesh->nNode; i++)
    nmax = max(nmax, Mesh->ParallelInfo->NodeLoc2Glob[i]+1);
  
  // reduce-max onto proc 0 to obtain number of global mesh nodes
  rbuf = ( (myRank == 0) ? (void *) &nNode_Glob : NULL);
  ierr = xf_Error(xf_MPI_Reduce(&nmax, rbuf, 1, xfe_SizeInt, xfe_MPI_MAX, 0));
  if (ierr != xf_OK) return ierr;
  
  // allocate NodeData_Glob on proc 0
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &NodeData_Glob, nNode_Glob*rank, DataSize));
    if (ierr != xf_OK) return ierr;
    
    // zero out NodeData_Glob
    switch (Size){
      case xfe_SizeInt:
        iData_Glob = (int *) NodeData_Glob;
        for (i=0; i<nNode_Glob*rank; i++) iData_Glob[i] = 0;
        break;
      case xfe_SizeReal:
        rData_Glob = (real *) NodeData_Glob;
        for (i=0; i<nNode_Glob*rank; i++) rData_Glob[i] = 0.;
        break;
      default:
        return xf_Error(xf_NOT_SUPPORTED);
        break;
    }
  }
  
  // gather number of nodes to inform proc 0 of # on each proc
  ierr = xf_Error(xf_MPI_Gather((void *) &Mesh->nNode, (void *) nNode, 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // proc 0 allocates NodeList
  if (myRank == 0){
    ierr = xf_Error(xf_VAlloc2( (void ***) &NodeList, nProc, nNode, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  // variable gather global node numbers onto proc 0
  rbuf = ((myRank == 0) ? (void *) NodeList[0] : NULL);
  ierr = xf_Error(xf_MPI_Gatherv((void *) Mesh->ParallelInfo->NodeLoc2Glob, 
                                 Mesh->nNode, rbuf, nNode, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // allocate NodeData_Buf, whose size is the sum of # nodes over all procs
  if (myRank == 0){
    for (iProc=0, sum=0; iProc<nProc; iProc++) sum += nNode[iProc];
    ierr = xf_Error(xf_Alloc((void **) &NodeData_Buf , sum*rank, DataSize));
    if (ierr != xf_OK) return ierr;
    if (Size == xfe_SizeInt)
      iData_Buf  = (int *) NodeData_Buf;
    else
      rData_Buf  = (real *) NodeData_Buf;
  }
  
  if (myRank == 0){
    // proc 0 sets its own node data
    if (Size == xfe_SizeInt)
      for (i=0; i<nNode[0]; i++)
        for (k=0; k<rank; k++) iData_Glob[NodeList[0][i]*rank+k] += iData[i*rank+k];
    else
      for (i=0; i<nNode[0]; i++) 
        for (k=0; k<rank; k++) rData_Glob[NodeList[0][i]*rank+k] += rData[i*rank+k];
    
    // proc 0 receives data from all other procs and adds to NodeData_Glob
    off = nNode[0]*rank;
    for (iProc=1; iProc<nProc; iProc++){
      ierr = xf_Error(xf_MPI_Recv(NodeData_Buf+off*DataSize, nNode[iProc]*rank*DataSize, iProc, 0));
      if (ierr != xf_OK) return ierr;
      off += nNode[iProc]*rank;
    }
    off = nNode[0]*rank;
    for (iProc=1; iProc<nProc; iProc++){
      if (Size == xfe_SizeInt)
        for (i=0; i<nNode[iProc]; i++) 
          for (k=0; k<rank; k++) iData_Glob[NodeList[iProc][i]*rank+k] += iData_Buf[off+i*rank+k];
      else
        for (i=0; i<nNode[iProc]; i++) 
          for (k=0; k<rank; k++) rData_Glob[NodeList[iProc][i]*rank+k] += rData_Buf[off+i*rank+k];
      off += nNode[iProc]*rank;
    } // iProc
    
  }
  else{
    ierr = xf_Error(xf_MPI_Send(NodeData, Mesh->nNode*rank*DataSize, 0, 0));
    if (ierr != xf_OK) return ierr;
  }
  
  
  // proc 0 sends data back
  if (myRank == 0){
    // copy ordered data in NodeData_Buf
    for (iProc=0, off=0; iProc<nProc; iProc++){
      if (Size == xfe_SizeInt)
        for (i=0; i<nNode[iProc]; i++) 
          for (k=0; k<rank; k++) iData_Buf[off+i*rank+k] = iData_Glob[NodeList[iProc][i]*rank+k];
      else
        for (i=0; i<nNode[iProc]; i++) 
          for (k=0; k<rank; k++) rData_Buf[off+i*rank+k] = rData_Glob[NodeList[iProc][i]*rank+k];
      off += nNode[iProc]*rank;
    }
    // set own data
    if (Size == xfe_SizeInt)
      for (i=0; i<nNode[0]; i++) 
        for (k=0; k<rank; k++) iData[i*rank+k] = iData_Buf[i*rank+k];
    else
      for (i=0; i<nNode[0]; i++) 
        for (k=0; k<rank; k++) rData[i*rank+k] = rData_Buf[i*rank+k];
    // send data to each proc
    off = nNode[0]*rank;
    for (iProc=1; iProc<nProc; iProc++){
      ierr = xf_Error(xf_MPI_Send(NodeData_Buf+off*DataSize, nNode[iProc]*rank*DataSize, iProc, 0));
      if (ierr != xf_OK) return ierr;
      off += nNode[iProc]*rank;
    }
  }
  else{
    ierr = xf_Error(xf_MPI_Recv(NodeData, Mesh->nNode*rank*DataSize, 0, 0));
    if (ierr != xf_OK) return ierr;
  }
  
  
  // release memory
  xf_Release( (void  *) nNode);
  xf_Release2((void **) NodeList);
  xf_Release( (void  *) NodeData_Glob);
  xf_Release( (void  *) NodeData_Buf);
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BcastMesh
int 
xf_BcastDataSet(xf_DataSet **pDataSet)
{
  // Only broadcasting Vectors
  int ierr, myRank, nProc, ibuf[4], i, len, your_turn, ReadWrite;
  enum xfe_Bool Done, Transfer;
  char *Title;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_Vector *V, *V_proc0;
  
  V = NULL;
  V_proc0 = NULL;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc > 1) {
    xf_printf("Broadcasting data ...");fflush(stdout);
    
    // Create dataset on all processors except root
    if (myRank != 0) {
      ierr = xf_Error(xf_CreateDataSet(pDataSet));
      if (ierr != xf_OK) return ierr;
    }
    DataSet = (*pDataSet);
    
    if (myRank == 0)
      D = DataSet->Head;
    
    Done = xfe_False;
    
    while (!Done){
      if (myRank == 0) {
        while (D != NULL){
          if (D->Type == xfe_Vector){
            V_proc0 = (xf_Vector *)D->Data;
            if (V_proc0->Linkage == xfe_LinkageGlobElem){
              Transfer = xfe_True;
              len = strlen(D->Title)+1;
              ReadWrite = D->ReadWrite;
              break;
            }
          }
          Transfer = xfe_False;
          D = D->Next;
        }
        if (D == NULL){
          Done = xfe_True;
          Transfer = xfe_False;
          len = -1;
        }
        
        ibuf[0] = Done;
        ibuf[1] = Transfer;
        ibuf[2] = len;
        ibuf[3] = ReadWrite;
      }
      
      // Broadcast Done and Transfer
      ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 4*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
      Done = ibuf[0];
      Transfer = ibuf[1];
      len = ibuf[2];
      ReadWrite = ibuf[3];
      
      // finish if all vectors have been transferred
      if (Done) 
        break;
      
      // if Transfer, enter here
      if (Transfer){
        if (myRank == 0)
          Title = D->Title;
        else {
          ierr = xf_Error(xf_Alloc((void **)&Title, xf_MAXSTRLEN, 
                                   sizeof(char)));
          if (ierr != xf_OK) return ierr;
        }
        
        // Bcast data title
        ierr = xf_Error(xf_MPI_Bcast((void *)Title, len, 0));
        if (ierr != xf_OK) return ierr;
        
        if (myRank == 0)
          Title = NULL;
        
        ierr = xf_Error(xf_CreateVector(&V));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_BcastVector(V_proc0, V));
        if (ierr != xf_OK) return ierr;
        
        /* The root processor will be searching for the next vector
         while the other processors add the vector to the dataset*/
        if (myRank != 0){
          // add to dataset
          ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_Vector, ReadWrite, 
                                        (void *)V, NULL));
          if (ierr != xf_OK) return ierr;
          
          xf_Release((void *)Title);
        }
        else {// proc0
          D = D->Next;
        }
      }// Transfer
    }// Done
    
    xf_MPI_Barrier();
    
    xf_printf("done.\n");
  }// nProc>1
  
  return xf_OK;
}  
