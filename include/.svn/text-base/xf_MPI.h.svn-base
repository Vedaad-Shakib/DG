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

#ifndef _xf_MPI_h
#define _xf_MPI_h 1

/*
  FILE:  xf_MPI.h

  This file contains headers for message passing interface (MPI)
  functions.  In the serial case, these functions are defined in
  xf_Serial.c

*/

#include "xf_MPIStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Init
extern int 
xf_MPI_Init(int *argc, char ***argv);
/*
PURPOSE:

  Initializes MPI and xflow-specific MPI structures

INPUTS:

  argc, argv : command line arguments from main

OUTPUTS:

  None

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Finalize
extern int 
xf_MPI_Finalize();
/*
PURPOSE:

  Finalizes MPI communication (last function called by main program)

INPUTS:

  None

OUTPUTS:

  None

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MPI_GetRank
extern int 
xf_MPI_GetRank(int *myRank, int *nProc);
/*
PURPOSE:

  Returns the current processor rank as well as the number of
  processors.  Both outputs are optional.  

INPUTS:

  None

OUTPUTS:

  (*myRank) : current processor rank (not set if NULL passed in)
  (*nProc)  : number of processors (not set if NULL passed in)

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Definition: xf_MPI_RequestSize
int 
xf_MPI_RequestSize();
/*
PURPOSE:

  Returns size of an MPI_Request structure

INPUTS:
  
  None.

OUTPUTS:

  None

RETURN:

  size of MPI_Request structure
*/

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Barrier
int 
xf_MPI_Barrier();
/*
PURPOSE:

  Waits until all processes call this function before returning.

INPUTS:
  
  None.

OUTPUTS:

  None.

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Definition: xf_MPI_Wait
int 
xf_MPI_Wait(void *request);
/*
PURPOSE:

  Waits until communication request is completed.

INPUTS:
  
  request: communication request on which to wait

OUTPUTS:

  None.

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Bcast
extern int 
xf_MPI_Bcast(void *buf, int size, int root);
/*
PURPOSE:

  Broadcasts contents of buf from root to all procs

INPUTS:

  buf  : buffer to broadcast
  size : number of bytes to broadcast
  root : proc from which to broadcast

OUTPUTS:

  On procs != root, buf is modified with contents from root

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Reduce
extern int 
xf_MPI_Reduce(void *sbuf, void *rbuf, int count, enum xfe_SizeType size, 
	      enum xfe_MPI_Op op, int root);
/*
PURPOSE:

  Reduces data from all procs onto root based on operation op.

INPUTS:

  sbuf  : data to reduce, valid on each proc
  count : number of values to reduce
  size  : xf-specific size idendifier (int, real, etc.)
  root  : proc which does reduction

OUTPUTS:

  rbuf  : data on proc root where result will be stored

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Allreduce
extern int 
xf_MPI_Allreduce(void *buf, int count, enum xfe_SizeType size, 
		 enum xfe_MPI_Op op);
/*
PURPOSE:

  Reduces data from all procs onto all procs based on operation op.
  Reduction is in place: buf is modified on all procs

INPUTS:

  buf   : data to reduce, valid on each proc
  count : number of values to reduce
  size  : xf-specific size idendifier (int, real, etc.)

OUTPUTS:

  buf   : data on each proc with stored result (original value 
          is overwritten)

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Send
extern int 
xf_MPI_Send(void *sbuf, int size, int dest, int tag);
/*
PURPOSE:

  Performs a blocking send of contents of sbuf, with size bytes, to
  proc dest

INPUTS:

  sbuf : buffer to send
  size : number of bytes to send
  dest : destination proc
  tag  : a message identifier tag

OUTPUTS:

  None

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Recv
extern int 
xf_MPI_Recv(void *rbuf, int size, int source, int tag);
/*
PURPOSE:

  Performs a blocking receive of contents of rbuf, with size bytes,
  from proc source

INPUTS:

  rbuf   : buffer into which receive is stored
  size   : number of bytes to receive
  source : source proc
  tag    : a message identifier tag

OUTPUTS:

  None

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Isend
extern int 
xf_MPI_Isend(void *sbuf, int size, int dest, int tag, void *request);
/*
PURPOSE:

  Begins a nonblocking send of contents of sbuf, with size bytes, to
  proc dest

INPUTS:

  sbuf : buffer to send
  size : number of bytes to send
  dest : destination proc
  tag  : a message identifier tag

OUTPUTS:

  request : a request handle for future queries (e.g. wait)

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Irecv
extern int 
xf_MPI_Irecv(void *rbuf, int size, int source, int tag, void *request);
/*
PURPOSE:

  Begins a nonblocking receive of contents of rbuf, with size bytes,
  from proc source

INPUTS:

  rbuf   : buffer into which receive is stored
  size   : number of bytes to receive
  source : source proc
  tag    : a message identifier tag

OUTPUTS:

  request : a request handle for future queries (e.g. wait)

RETURN:

  Error code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Scatter
extern int 
xf_MPI_Scatter(void *sbuf, void *rbuf, int size, int root);
/*
PURPOSE:

  Scatters contents of sbuf to procs.
 
INPUTS:

  sbuf : Send buffer.  Relevant only on root.  Contains nProc*size
         bytes, ordered sequentially by proc.
  size : number of bytes to scatter
  root : proc from which to scatter

OUTPUTS:

  rbuf : Receive buffer.  Relevant on all procs.  Must be alloced to
         hold size bytes.  Will contain proc's portion of sbuf.

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Gather
extern int 
xf_MPI_Gather(void *sbuf, void *rbuf, int size, int root);
/*
PURPOSE:

  Gathers size bytes from sbuf on each proc into rbuf on root.
 
INPUTS:

  sbuf : Send buffer.  Relevant on all procs.  Contains size bytes
  size : number of bytes to gather
  root : proc that does the gathering

OUTPUTS:

  rbuf : Receive buffer.  Relevant only on root.  Must be alloced 
         for nProc*size bytes.  Will be filled in sequentially 
	 (0 to nProc-1) with data from sbufs.
  
RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Scatterv
extern int 
xf_MPI_Scatterv(void *sbuf, int *scount, void *rbuf, int rcount, 
		int size, int root);
/*
PURPOSE:

  Scatters contents of sbuf to procs.  Each proc can get a
  variable-length portion of sbuf, as dictated by scount[iProc] =
  rcount.
 
INPUTS:

  sbuf : Send buffer.  Relevant only on root.  Contains 
         sbuf * sum(scount) bytes, ordered sequentially by proc.
  scount : Relevant only at root.  scount[iProc] is the number of
           data items to send to iProc
  rcount : Relevant on all procs.  Number of data items to receive.
           This is the same as scount[iProc] (but must be defined 
	   on each proc!)
  size : number of bytes to scatter
  root : proc from which to scatter

OUTPUTS:

  rbuf : Receive buffer.  Relevant on all procs.  Must be alloced to
         hold rcount bytes.  Will contain proc's portion of sbuf.

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Gatherv
extern int 
xf_MPI_Gatherv(void *sbuf, int scount, void *rbuf, int *rcount, 
	       int size, int root);
/*
PURPOSE:

  Gathers size*scount bytes from sbuf on each proc into rbuf on root.
  scount = rcount[iProc] may be different on each proc.
 
INPUTS:

  sbuf : Send buffer.  Relevant on all procs.  Contains 
         size*scount bytes
  size : number of bytes to gather
  root : proc that does the gathering

OUTPUTS:

  rcount : rcount[iProc] is the # of data items to receive from iProc.
           Relevant only at root
  rbuf : Receive buffer.  Relevant only on root.  Must be alloced 
         for size*sum(rcount) bytes.  Will be filled in sequentially 
	 (0 to nProc-1) with data from sbufs.
  

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MPI_PScatterv
extern int 
xf_MPI_PScatterv(void *sbuf, int *sindex, int *scount, void *rbuf, 
		 int rcount, int size, int root);
/*
PURPOSE:

  Packs and scatters select contents of sbuf to procs.  Each proc can
  receive a variable-length portion of sbuf, as dictated by
  scount[iProc] = rcount.  Furthermore, the sent data in sbuf need not
  be contiguous: rather, only the values of sbuf indexed by sindex are
  sent.  The data is packed before it is sent, and on the receiving
  end, will have contiguous data.
 
INPUTS:

  sbuf : Send buffer.  Relevant only on root.  Size is arbitrary, but
         must be consistent with sindex.
  sindex : The first scount[0] values are entries of sbuf that need
           to be sent to proc 0.  The next scount[1] values are entries
	   of sbuf that need to be sent to proc1, etc.  sindex is
	   relevant only at root.
  scount : Relevant only at root.  scount[iProc] is the number of
           data items to send to iProc
  rcount : Relevant on all procs.  Number of data items to receive.
           This is the same as scount[iProc] (but must be defined 
	   on each proc!)
  size : number of bytes to scatter
  root : proc from which to scatter

OUTPUTS:

  rbuf : Receive buffer.  Relevant on all procs.  Must be alloced to
         hold rcount bytes.  Will contain proc's portion of sbuf.

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_PGatherv
extern int 
xf_MPI_PGatherv(void *sbuf, int scount, void *rbuf, 
		int *rindex, int *rcount, int size, int root);
/*
PURPOSE:

  Gathers and unpacks contents of sbuf from each proc into rbuf on
  proc root.  Each proc can send a variable-length sbuf.  proc root
  must be aware of how many values it is receiving from each proc, as
  well as where they belong in rbuf: this is the job of rcount (number
  of values from each proc) and rindex (position of these values in
  rbuf).

 
INPUTS:

  sbuf : Send buffer.  Relevant at all procs.  Size is scount.
  scount : Relevant at all procs.  Number of data items to send to root
  rindex : The first rcount[0] values are the locations in rbuf where proc
           0's data will go.  The next rcount[1] values are the locations
	   in rbuf where proc 1's data will go, etc.  Only relevant at root.
  rcount : Relevant only on root.  rcount[iProc] is the number of data 
           items to receive from each iProc. 
  size : number of bytes to gather
  root : proc that does the gathering

OUTPUTS:

  rbuf : Receive buffer.  Relevant only on root.  Must be alloced to
         hold enough data from all procs 

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_DPScatterv
extern int 
xf_MPI_DPScatterv(void *sbuf, int *slength, int *sindex, int *scount, 
		  void *rbuf, int *rlength, int rcount, int size, int root);
/*
PURPOSE:

  "Double-level" pack and scatter.  i.e. variable-size version of
  PScatterv.  Root contains a variable-length 2D array (e.g. created
  with VAlloc2) over all elements.  It needs to scatter the contents
  of this array to all the procs, each of which has scount[proc]
  elements (specified in the unrolled array sindex), with each element
  containing slength[elem] values.  Each proc knows how many elements
  it's expecting (rcount), and the number of values per element
  (rlength vector).  The result is put into rbuf, which is a
  once-dereferenced variable-length 2D array.
 
INPUTS:

  sbuf : Send buffer.  Relevant only on root.  Size is arbitrary, but
         must be consistent with sindex.
  slength : a lookup-length table that specifies how many values are
            in each element indexed by sindex
  sindex : The first scount[0] values are elements of sbuf that need
           to be sent to proc 0.  The next scount[1] values are elements
	   of sbuf that need to be sent to proc1, etc.  sindex is
	   relevant only at root.  Note, each element has a different
	   number of values, as indicated by slength.
  scount : Relevant only at root.  scount[iProc] is the number of
           data items to send to iProc
  rlength : Relevant on all procs.  Number of values per element.
           Vector of size rcount.
  rcount : Relevant on all procs.  Number of data items (elements) to
           receive.  This is the same as scount[iProc] (but must be
           defined on each proc!)
  size : number of bytes to scatter
  root : proc from which to scatter

OUTPUTS:

  rbuf : Receive buffer.  Relevant on all procs.  Must be alloced to
         hold rcount*sum(rlength) bytes.  Will contain proc's portion
         of sbuf.

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_DPGatherv
extern int 
xf_MPI_DPGatherv(void *sbuf, int *slength, int scount, void *rbuf, 
		 int *rlength, int *rindex, int *rcount, int size, int root);
/*
PURPOSE:

  "Double-level" gather and unpack.  Gathers and unpacks contents of
  sbuf from each proc into rbuf on proc root.  Each proc can send a
  variable-number of elements, and moreover, each element can have a
  variable number of values, as indicated in slength.  sbuf stores the
  unrolled 2D variable array that contains the send values.  Proc root
  must be aware of how many elements it is receiving from each proc,
  as well as where they belong in rbuf: this is the job of rcount
  (number of elements from each proc), rindex (position of these
  elements in rbuf), rlength (number of values per element in rbuf).

 
INPUTS:

  sbuf : Send buffer.  Relevant at all procs.  Size is scount.
  slength : Relevant on all procs.  Number of values per element.
            Vector of size scount.
  scount : Relevant at all procs.  Number of elements to send to root
  rlength : a lookup-length table that specifies how many values are
            in each element indexed by rindex
  rindex : The first rcount[0] elements are the locations in rbuf where proc
           0's data will go.  The next rcount[1] elements are the locations
	   in rbuf where proc 1's data will go, etc.  Only relevant at root.
  rcount : Relevant only on root.  rcount[iProc] is the number of elements
           to receive from each iProc. 
  size : number of bytes to gather
  root : proc that does the gathering

OUTPUTS:

  rbuf : Receive buffer.  Relevant only on root.  Must be alloced to
         hold enough data from all procs 

RETURN:

  Error code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Alltoall
extern int 
xf_MPI_Alltoall(void *sbuf, void *rbuf, int size);
/*
PURPOSE:

  Each proc sends contents of sbuf to all procs.  Each proc receives
  data into rbof from all procs.  Size of data is size.
 
INPUTS:

  sbuf : Send buffer.  Relevant on all procs.  Contains size bytes.
  size : number of bytes to send

OUTPUTS:

  rbuf : Receive buffer.  Relevant on all procs.  Contains size bytes.

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MPI_Alltoallv
extern int 
xf_MPI_Alltoallv(void *sbuf, int *scount, void *rbuf, int *rcount, int size);
/*
PURPOSE:

  Each proc sends contents of sbuf to all procs.  Each proc receives
  data into rbuf from all procs.  Size of data sent by proc1 to proc2
  is scount[proc2] from proc1's point of view, which equals
  rcount[proc1] from proc2's point of view.
 
INPUTS:

  sbuf : Send buffer.  Relevant on all procs.  Contains sum(scount) entries.
  scount : Relevant on all procs.  scount[proc] is the number of entries
           to send to proc.
  size : number of bytes to send

OUTPUTS:

  rbuf : Receive buffer.  Relevant on all procs.  Contains sum(rcount) entries.
  rcount : Relevant on all procs.  rcount[proc] is the number of entries
           to receive from proc.

RETURN:

  Error code
*/




#endif // end ifndef _xf_MPI_h
