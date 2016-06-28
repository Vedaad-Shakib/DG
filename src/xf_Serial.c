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
  FILE:  xf_Serial.c

  This file contains mainly dummy functions as placeholders for
  parallel functions that either should not be called in serial or
  that return immediately when called in serial.

*/


#include "xf_AllStruct.h"
#include "xf_MPIStruct.h"
#include "xf_IO.h"


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Init
int 
xf_MPI_Init(int *argc, char ***argv) 
{
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Finalize
int 
xf_MPI_Finalize()
{
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_GetRank
int 
xf_MPI_GetRank(int *myRank, int *nProc)
{
  if (myRank != NULL) (*myRank) = 0;
  if (nProc  != NULL) (*nProc)  = 1;
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_RequestSize
int 
xf_MPI_RequestSize()
{
  return 0;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Barrier
int 
xf_MPI_Barrier()
{
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Wait
int 
xf_MPI_Wait( void *request)
{
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Bcast
int 
xf_MPI_Bcast(void *buf, int size, int root)
{
  int ierr;
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Reduce
int 
xf_MPI_Reduce(void *sbuf, void *rbuf, int count, enum xfe_SizeType size, 
	      enum xfe_MPI_Op op, int root)
{
  return xf_Error(xf_NOT_SUPPORTED); // for now
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Allreduce
int 
xf_MPI_Allreduce(void *buf, int count, enum xfe_SizeType size, 
		 enum xfe_MPI_Op op)
{
  return xf_OK; // nothing to be done
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Send
int 
xf_MPI_Send(void *sbuf, int size, int dest, int tag)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Recv
int 
xf_MPI_Recv(void *rbuf, int size, int source, int tag)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Isend
int 
xf_MPI_Isend(void *sbuf, int size, int dest, int tag, void *request)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Irecv
int 
xf_MPI_Irecv(void *rbuf, int size, int source, int tag, void *request)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Scatter
int 
xf_MPI_Scatter(void *sbuf, void *rbuf, int size, int root)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Gather
int 
xf_MPI_Gather(void *sbuf, void *rbuf, int size, int root)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Scatterv
int 
xf_MPI_Scatterv(void *sbuf, int *scount, void *rbuf, int rcount, 
		int size, int root)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Gatherv
int 
xf_MPI_Gatherv(void *sbuf, int scount, void *rbuf, int *rcount, 
	       int size, int root)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_PScatterv
int 
xf_MPI_PScatterv(void *sbuf, int *sindex, int *scount, void *rbuf, 
		 int rcount, int size, int root)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_PGatherv
int 
xf_MPI_PGatherv(void *sbuf, int scount, void *rbuf, 
		int *rindex, int *rcount, int size, int root)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_DPScatterv
int 
xf_MPI_DPScatterv(void *sbuf, int *slength, int *sindex, int *scount, 
		  void *rbuf, int *rlength, int rcount, int size, int root)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}

/******************************************************************/
//   FUNCTION Definition: xf_MPI_DPGatherv
int 
xf_MPI_DPGatherv(void *sbuf, int *slength, int scount, void *rbuf, 
		 int *rlength, int *rindex, int *rcount, int size, int root)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Alltoall
int 
xf_MPI_Alltoall(void *sbuf, void *rbuf, int size)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}


/******************************************************************/
//   FUNCTION Definition: xf_MPI_Alltoallv
int 
xf_MPI_Alltoallv(void *sbuf, int *scount, void *rbuf, int *rcount, int size)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}

/******************************************************************/
//   FUNCTION Definition: xf_PartitionMesh
int 
xf_PartitionMesh( xf_Mesh *Mesh, int ***pElem2Proc)
{
  xf_printf("Error: parallel function called in a serial library.\n");
  return xf_Error(xf_CODE_LOGIC_ERROR);
}

/******************************************************************/
//   FUNCTION Definition: xf_CPUWorkReport
int
xf_CPUWorkReport(xf_All *All, xf_Vector *State, enum xfe_Bool *load_imbalance)
{
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CPULoadBalance
int
xf_CPULoadBalance(xf_All *All, xf_Vector **pState, int nState)
{
  return xf_OK;
}

int
xfYu_CPULoadBalance(xf_All *All, xf_Vector **pState, int nState)
{
  return xf_OK;
}