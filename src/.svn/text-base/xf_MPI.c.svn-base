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
  FILE:  xf_MPI.c

  This file contains functions for inter-processor communication using
  a message passing interface (MPI)

*/


#include "xf_AllStruct.h"
#include "xf_MPIStruct.h"
#include "xf_Memory.h"
#include "xf_IO.h"
#include "string.h"
#include "mpi.h"


// global variables and definitions
MPI_Comm xf_MPI_COMM_WORLD;
MPI_Comm xf_MPI_COMM_SELF;

#define xf_MPI_Request MPI_Request


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Init
int 
xf_MPI_Init(int *argc, char ***argv) 
{
  int ierr;
  
  // initialize MPI
  ierr = MPI_Init(argc, argv);
  if (ierr != MPI_SUCCESS) return xf_Error(xf_MPI_ERROR);

  // duplicate communicator handles to create xflow-specific ones
  ierr = MPI_Comm_dup(MPI_COMM_WORLD, &xf_MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) return xf_Error(xf_MPI_ERROR);
  
  ierr = MPI_Comm_dup(MPI_COMM_SELF, &xf_MPI_COMM_SELF);
  if (ierr != MPI_SUCCESS) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Finalize
int 
xf_MPI_Finalize()
{
  int ierr;

  /* Finish up MPI */  
  ierr = MPI_Finalize();
  if (ierr != MPI_SUCCESS) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_GetRank
int 
xf_MPI_GetRank(int *myRank, int *nProc)
{
  int ierr;

  if (myRank != NULL){
    ierr = MPI_Comm_rank(xf_MPI_COMM_WORLD, myRank);
    if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);
  }
  if (nProc  != NULL){
    ierr = MPI_Comm_size(xf_MPI_COMM_WORLD, nProc);
    if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_RequestSize
int 
xf_MPI_RequestSize()
{
  return sizeof(xf_MPI_Request);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Barrier
int 
xf_MPI_Barrier()
{
  int ierr;

  ierr = MPI_Barrier(xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Wait
int 
xf_MPI_Wait( void *request)
{
  int ierr;
  MPI_Status status;

  ierr = MPI_Wait( (xf_MPI_Request *) request, &status);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Bcast
int 
xf_MPI_Bcast(void *buf, int size, int root)
{
  int ierr;

  ierr = MPI_Bcast(buf, size, MPI_CHAR, root, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Type2MPI
static int 
xf_Type2MPI(enum xfe_SizeType size, MPI_Datatype *mpisize, int *nbyte)
{
  switch (size){
  case xfe_SizeInt:
    (*mpisize) = MPI_INT;
    (*nbyte) = sizeof(int);
    break;
  case xfe_SizeReal:
    (*nbyte) = sizeof(real);
    if (sizeof(real) == sizeof(double))
      (*mpisize) = MPI_DOUBLE;
    else if (sizeof(real) == sizeof(float))
      (*mpisize) = MPI_FLOAT;
    else
      return xf_Error(xf_NOT_SUPPORTED);
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Op2MPI
static int 
xf_Op2MPI(enum xfe_MPI_Op op, MPI_Op *mpiop)
{
  switch (op){
  case xfe_MPI_MAX:
    (*mpiop) = MPI_MAX;
    break;
  case xfe_MPI_MIN:
    (*mpiop) = MPI_MIN;
    break;
  case xfe_MPI_SUM:
    (*mpiop) = MPI_SUM;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Reduce
int 
xf_MPI_Reduce(void *sbuf, void *rbuf, int count, enum xfe_SizeType size, 
	      enum xfe_MPI_Op op, int root)
{
  int ierr, nbyte;
  MPI_Op mpiop;
  MPI_Datatype mpisize;

  // obtain MPI datatype
  ierr = xf_Error(xf_Type2MPI(size, &mpisize, &nbyte));
  if (ierr != xf_OK) return ierr;

  // obtain MPI operation type
  ierr = xf_Error(xf_Op2MPI(op, &mpiop));
  if (ierr != xf_OK) return ierr;

  // MPI_Reduce (&sendbuf,&recvbuf,count,datatype,op,root,comm)   
  ierr = MPI_Reduce(sbuf, rbuf, count, mpisize, mpiop, root, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Allreduce
int 
xf_MPI_Allreduce(void *buf, int count, enum xfe_SizeType size, 
		 enum xfe_MPI_Op op)
{
  int ierr, nbyte;
  MPI_Op mpiop;
  MPI_Datatype mpisize;
  void *rbuf;

  // obtain MPI datatype
  ierr = xf_Error(xf_Type2MPI(size, &mpisize, &nbyte));
  if (ierr != xf_OK) return ierr;

  // obtain MPI operation type
  ierr = xf_Error(xf_Op2MPI(op, &mpiop));
  if (ierr != xf_OK) return ierr;

  // allocate rbuf
  ierr = xf_Error(xf_Alloc(&rbuf, count, nbyte));
  if (ierr != xf_OK) return ierr;

  // MPI_Allreduce (&sendbuf,&recvbuf,count,datatype,op,comm) 
  ierr = MPI_Allreduce(buf, rbuf, count, mpisize, mpiop, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);
  
  // copy data back into buf
  memcpy(buf, rbuf, nbyte*count);

  // destroy rbuf
  xf_Release(rbuf);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Send
int 
xf_MPI_Send(void *sbuf, int size, int dest, int tag)
{
  int ierr;

  //  MPI_Send (&buf,count,datatype,dest,tag,comm) 
  ierr = MPI_Send(sbuf, size, MPI_CHAR, dest, tag, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Recv
int 
xf_MPI_Recv(void *rbuf, int size, int source, int tag)
{
  int ierr;
  MPI_Status status;

  // MPI_Recv (&buf,count,datatype,source,tag,comm,&status) 
  ierr = MPI_Recv(rbuf, size, MPI_CHAR, source, tag, xf_MPI_COMM_WORLD, &status);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Isend
int 
xf_MPI_Isend(void *sbuf, int size, int dest, int tag, void *request)
{
  int ierr;

  // MPI_Isend (&buf,count,datatype,dest,tag,comm,&request) 
  ierr = MPI_Isend(sbuf, size, MPI_CHAR, dest, tag, xf_MPI_COMM_WORLD, request);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Irecv
int 
xf_MPI_Irecv(void *rbuf, int size, int source, int tag, void *request)
{
  int ierr;

  // MPI_Irecv (&buf,count,datatype,source,tag,comm,&request) 
  ierr = MPI_Irecv(rbuf, size, MPI_CHAR, source, tag, xf_MPI_COMM_WORLD, request);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_MPI_Scatter
int 
xf_MPI_Scatter(void *sbuf, void *rbuf, int size, int root)
{
  int ierr;

  ierr = MPI_Scatter(sbuf, size, MPI_CHAR, rbuf, size, MPI_CHAR, 
		     root, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Gather
int 
xf_MPI_Gather(void *sbuf, void *rbuf, int size, int root)
{
  int ierr;

  ierr = MPI_Gather(sbuf, size, MPI_CHAR, rbuf, size, MPI_CHAR, 
		    root, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Scatterv
int 
xf_MPI_Scatterv(void *sbuf, int *scount, void *rbuf, int rcount, 
		int size, int root)
{
  int ierr, i, tot;
  int myRank, nProc;
  int *displs;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // multiply scount by size, since we work in bytes
  if (myRank == root)
    for (i=0; i<nProc; i++) scount[i] *= size;
  rcount *= size; // also multiply rcount = receive count

  /* on root, calculate displs = see below */
  if (myRank == root){
    ierr = xf_Error(xf_Alloc((void **) &displs, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0, tot=0; i<nProc; i++){
      displs[i] = tot;
      tot += scount[i];
    }
  }
  else displs = NULL;
  
  /*
    sendbuf 	address of send buffer (choice, significant only at root)
    sendcounts 	integer array (of length group size) specifying the number 
    of elements to send to each processor
    displs      integer array (of length group size). Entry i specifies 
                the displacement (relative to sendbuf from which to take the 
		outgoing data to process i
    sendtype 	data type of send buffer elements (handle)
    recvcount 	number of elements in receive buffer (integer)
    recvtype 	data type of receive buffer elements (handle)
    root        rank of sending process (integer)
    comm 	communicator (handle) 
  */
  ierr = MPI_Scatterv(sbuf, scount, displs, MPI_CHAR, rbuf, rcount, MPI_CHAR, 
		      root, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  // divide by size to restore scount 
  if (myRank == root)
    for (i=0; i<nProc; i++) scount[i] /= size;

  // release memory on root
  if (myRank == root) xf_Release( (void *) displs);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Gatherv
int 
xf_MPI_Gatherv(void *sbuf, int scount, void *rbuf, int *rcount, 
	       int size, int root)
{
  int ierr, i, tot;
  int myRank, nProc;
  int *displs;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // multiply rcount by size, since we work in bytes
  if (myRank == root)
    for (i=0; i<nProc; i++) rcount[i] *= size;
  scount *= size; // also multiply scount = send count

  /* on root, calculate displs = see below */
  if (myRank == root){
    ierr = xf_Error(xf_Alloc((void **) &displs, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0, tot=0; i<nProc; i++){
      displs[i] = tot;
      tot += rcount[i];
    }
  }
  else displs = NULL;

  /*
    sendbuf 	starting address of send buffer (choice)
    sendcount 	number of elements in send buffer (integer)
    sendtype 	data type of send buffer elements (handle)
    recvcounts 	integer array (of length group size) containing the 
                number of elements that are received from each process 
		(significant only at root)
    displs 	integer array (of length group size). Entry i specifies the 
                displacement relative to recvbuf at which to place the incoming 
		data from process i (significant only at root)
    recvtype 	data type of recv buffer elements (significant only at root) (handle)
    root 	rank of receiving process (integer)
    comm 	communicator (handle) 
  */
  ierr = MPI_Gatherv(sbuf, scount, MPI_CHAR, rbuf, rcount, displs, 
		     MPI_CHAR, root, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  // divide by size to restore rcount 
  if (myRank == root)
    for (i=0; i<nProc; i++) rcount[i] /= size;
  
  if (myRank == root) xf_Release( (void *) displs);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_PScatterv
int 
xf_MPI_PScatterv(void *sbuf, int *sindex, int *scount, void *rbuf, 
		 int rcount, int size, int root)
{
  int ierr, i, j, tot;
  int myRank, nProc;
  int isent;
  void *buf;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // root allocs and packs buf according to sindex, scount
  if (myRank == root){
    for (i=0, tot=0; i<nProc; i++) tot += scount[i];

    ierr = xf_Error(xf_Alloc(&buf, tot, size));
    if (ierr != xf_OK) return ierr;

    isent  = 0;
    for (i=0; i<nProc; i++){
      for (j=0; j<scount[i]; j++)
	memcpy(buf+size*(isent+j), sbuf + size*sindex[isent + j], size);
      isent += scount[i];
    }
  }
 
  ierr = xf_Error(xf_MPI_Scatterv(buf, scount, rbuf, rcount, size, root));
  if (ierr != xf_OK) return ierr;

  // release memory on root
  if (myRank == root) xf_Release( buf );

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_PGatherv
int 
xf_MPI_PGatherv(void *sbuf, int scount, void *rbuf, 
		int *rindex, int *rcount, int size, int root)
{
  int ierr, i, j, tot;
  int myRank, nProc;
  int irecvd;
  void *buf;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  buf = NULL;

  // root allocs a packed receive buffer 
  if (myRank == root){
    for (i=0, tot=0; i<nProc; i++) tot += rcount[i];

    ierr = xf_Error(xf_Alloc(&buf, tot, size));
    if (ierr != xf_OK) return ierr;
  }

  // gather data into the packed receive buffer
  ierr = xf_Error(xf_MPI_Gatherv(sbuf, scount, buf, rcount, size, root));
  if (ierr != xf_OK) return ierr;

  // root redistributes its packed receive buffer into rbuf, using rindex
  if (myRank == root){
    irecvd  = 0;
    for (i=0; i<nProc; i++){
      for (j=0; j<rcount[i]; j++){
	memcpy(rbuf + size*rindex[irecvd + j], buf+size*(irecvd+j), size);
      }
      irecvd += rcount[i];
    }
  } 

  // release memory on root
  if (myRank == root) xf_Release( buf );

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_MPI_DPScatterv
int 
xf_MPI_DPScatterv(void *sbuf, int *slength, int *sindex, int *scount, 
		  void *rbuf, int *rlength, int rcount, int size, int root)
{
  int ierr, i, j, k, ind, tot;
  int myRank, nProc;
  int isent, imax;
  int rtcount, *stcount = NULL, *spos = NULL;
  void *buf = NULL;
  real *rtemp = NULL;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // root allocs and packs buf according to slength, sindex, scount
  if (myRank == root){

    // stcount[proc] = total # of size units sending to proc
    ierr = xf_Error(xf_Alloc( (void **) &stcount, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    // tot = total number of units summed over all procs
    isent = 0;
    for (i=0, tot=0; i<nProc; i++){ 
      stcount[i] = 0;
      for (j=0; j<scount[i]; j++)
	stcount[i] += slength[sindex[isent+j]];
      tot += size*stcount[i];
      isent += scount[i];
    }

    ierr = xf_Error(xf_Alloc(&buf, tot, size));
    if (ierr != xf_OK) return ierr;

    // imax = largest index that appears in sindex
    isent = 0;
    for (i=0, imax=0; i<nProc; i++){
      for (j=0; j<scount[i]; j++)
	imax = max(imax, sindex[isent+j]+1);
      isent += scount[i];
    }
    
    // spos[i] = send buffer start pos of ith index, 0 <= i < imax
    ierr = xf_Error(xf_Alloc( (void **) &spos, imax, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0, tot=0; i<imax; i++){
      spos[i] = tot;
      tot += slength[i];
    }


    isent = 0;
    k = 0;
    for (i=0; i<nProc; i++){
      for (j=0; j<scount[i]; j++){
	ind = sindex[isent + j];
	memcpy(buf+size*k, sbuf + size*spos[ind], size*slength[ind]);
	k += slength[ind];
      }
      isent += scount[i];
    }
  }

  // determine rtcount = total # of size units expecting on this proc
  rtcount = 0;
  for (j=0; j<rcount; j++) rtcount += rlength[j];
 
  ierr = xf_Error(xf_MPI_Scatterv(buf, stcount, rbuf, rtcount, size, root));
  if (ierr != xf_OK) return ierr;

  // release memory on root
  if (myRank == root){
    xf_Release( stcount );
    xf_Release( spos );
    xf_Release( buf );
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_DPGatherv
int 
xf_MPI_DPGatherv(void *sbuf, int *slength, int scount, void *rbuf, 
		 int *rlength, int *rindex, int *rcount, int size, int root)
{
  int ierr, i, j, k, imax, ind, tot;
  int myRank, nProc;
  int irecvd;
  int stcount, *rtcount = NULL, *rpos = NULL;
  void *buf;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  buf = NULL;

  // root allocs a packed receive buffer 
  if (myRank == root){

    // rtcount[proc] = total # of size units expecting from proc
    ierr = xf_Error(xf_Alloc( (void **) &rtcount, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    // tot = total number of units summed over all procs
    irecvd = 0;
    for (i=0, tot=0; i<nProc; i++){ 
      rtcount[i] = 0;
      for (j=0; j<rcount[i]; j++)
	rtcount[i] += rlength[rindex[irecvd+j]];
      tot += size*rtcount[i];
      irecvd += rcount[i];
    }

    // this is the packed receive buffer
    ierr = xf_Error(xf_Alloc(&buf, tot, size));
    if (ierr != xf_OK) return ierr;
  }

  // determine stcount = total # of size units sending from this proc
  stcount = 0;
  for (j=0; j<scount; j++) stcount += slength[j];

  // gather data into the packed receive buffer
  ierr = xf_Error(xf_MPI_Gatherv(sbuf, stcount, buf, rtcount, size, root));
  if (ierr != xf_OK) return ierr;

  // root redistributes its packed receive buffer into rbuf, using rindex
  if (myRank == root){

    // imax = largest index that appears in rindex
    irecvd = 0;
    for (i=0, imax=0; i<nProc; i++){
      for (j=0; j<rcount[i]; j++)
	imax = max(imax, rindex[irecvd+j]+1);
      irecvd += rcount[i];
    }
    
    // rpos[i] = receive buffer start pos of ith index, 0 <= i < imax
    ierr = xf_Error(xf_Alloc( (void **) &rpos, imax, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0, tot=0; i<imax; i++){
      rpos[i] = tot;
      tot += rlength[i];
    }

    irecvd = 0;
    k = 0;
    for (i=0; i<nProc; i++){
      for (j=0; j<rcount[i]; j++){
	ind = rindex[irecvd + j];
	memcpy(rbuf + size*rpos[ind], buf + size*k, size*rlength[ind]);
	k += rlength[ind];
      }
      irecvd += rcount[i];
    }
  } 

  // release memory on root
  if (myRank == root){
    xf_Release( rtcount );
    xf_Release( rpos );
    xf_Release( buf );
  }

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_MPI_Alltoall
int 
xf_MPI_Alltoall(void *sbuf, void *rbuf, int size)
{
  int ierr;

  ierr = MPI_Alltoall(sbuf, size, MPI_CHAR, rbuf, size, MPI_CHAR, 
		      xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Alltoallv
int 
xf_MPI_Alltoallv(void *sbuf, int *scount, void *rbuf, int *rcount, int size)
{
  int ierr, i, tot;
  int myRank, nProc;
  int *sdispls, *rdispls;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // multiply scount, rcount by size, since we work in bytes
  for (i=0; i<nProc; i++) scount[i] *= size;
  for (i=0; i<nProc; i++) rcount[i] *= size;

  /* Create sdispls, rdispls */
  ierr = xf_Error(xf_Alloc((void **) &sdispls, nProc, sizeof(int)));
  ierr = xf_Error(xf_Alloc((void **) &rdispls, nProc, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0, tot=0; i<nProc; i++){
    sdispls[i] = tot;
    tot += scount[i];
  }
  for (i=0, tot=0; i<nProc; i++){
    rdispls[i] = tot;
    tot += rcount[i];
  }

  /*  
  int MPI_Alltoallv( void *sendbuf, int *sendcnts,  int *sdispls,   MPI_Datatype sendtype, 
                     void *recvbuf, int *recvcnts,  int *rdispls,   MPI_Datatype recvtype, 
	             MPI_Comm comm )
  sendbuf 	starting address of send buffer (choice)
  sendcounts 	integer array equal to the group size specifying the number 
                of elements to send to each processor
  sdispls 	integer array (of length group size). Entry j specifies the 
                displacement (relative to sendbuf from which to take the 
		outgoing data destined for process j
  sendtype 	data type of send buffer elements (handle)
  recvcounts 	integer array equal to the group size specifying the maximum 
                number of elements that can be received from each processor
  rdispls 	integer array (of length group size). Entry i specifies the 
                displacement (relative to recvbuf at which to place the 
		incoming data from process i
  recvtype 	data type of receive buffer elements (handle)
  comm 	        communicator (handle) 
  */

 /*  xf_pprintf("Before MPI_AlltoAllv.\n"); */
/*   for (i=0; i<nProc; i++) xf_pprintf("scount[%d] = %d\n", i, scount[i]); */
/*   for (i=0; i<nProc; i++) xf_pprintf("sdispls[%d] = %d\n", i, sdispls[i]); */
/*   for (i=0; i<nProc; i++) xf_pprintf("rcount[%d] = %d\n", i, rcount[i]); */
/*   for (i=0; i<nProc; i++) xf_pprintf("rdispls[%d] = %d\n", i, rdispls[i]); */
  
  ierr = MPI_Alltoallv(sbuf, scount, sdispls, MPI_CHAR, rbuf, rcount, 
		       rdispls, MPI_CHAR, xf_MPI_COMM_WORLD);
  if(ierr != MPI_SUCCESS ) return xf_Error(xf_MPI_ERROR);

  /* xf_pprintf("After MPI_AlltoAllv.\n"); */

  // divide by size to restore scount, rcount
  for (i=0; i<nProc; i++) scount[i] /= size;
  for (i=0; i<nProc; i++) rcount[i] /= size;

  /* Release memory */
  xf_Release( (void *) sdispls);
  xf_Release( (void *) rdispls);

  return xf_OK;
}


