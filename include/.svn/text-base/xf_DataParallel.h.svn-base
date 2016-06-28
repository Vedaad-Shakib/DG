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

#ifndef _xf_DataParallel_h
#define _xf_DataParallel_h 1

/*
  FILE:  xf_DataParallel.h

  This file contains the headers for parallelization functions dealing
  with parameters.

*/



/******************************************************************/
//   FUNCTION Prototype: xf_HaloExchangeVectorBegin
extern int 
xf_HaloExchangeVectorBegin( xf_Vector *V);
/*
PURPOSE:

  Begins a non-blocking exchange of halo data in V.  Sets
  V->HaloInTransit flag to True.  Returns immediately of not running
  in parallel.  If running in parallel but V->ParallelFlag == False,
  errors out.

INPUTS:
 
  V   : vector for which to begin the halo exchange 

OUTPUTS: 

  None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_HaloExchangeVectorEnd
extern int 
xf_HaloExchangeVectorEnd( xf_Vector *V);
/*
PURPOSE:

  Ends a non-blocking exchange of halo data in V (i.e. waits on all
  commrequests).  Only waits if V->HaloInTransit flag is True; sets
  this flag to False afterwards.  Returns immediately of not running
  in parallel.  If running in parallel but V->ParallelFlag == False,
  errors out.

INPUTS:
 
  V   : vector for which to end the halo exchange 

OUTPUTS: 

  None

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_HaloReverseExchangeVectorBegin
extern int 
xf_HaloReverseExchangeVectorBegin( xf_Vector *V);
/*
PURPOSE:

  Begins a non-blocking reverse exchange of halo data in V.  Sets
  V->HaloInTransit flag to True.  Returns immediately of not running
  in parallel.  If running in parallel but V->ParallelFlag == False,
  errors out.

  The exchane is in "reverse" in that data on halo elements is sent to
  the corresponding regular elements on other processors, which see
  their values incremented by the received halo data.

INPUTS:
 
  V   : vector for which to begin the halo exchange 

OUTPUTS: 

  None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_HaloReverseExchangeVectorEnd
extern int 
xf_HaloReverseExchangeVectorEnd( xf_Vector *V);
/*
PURPOSE:

  Ends a non-blocking reverse exchange of halo data in V (i.e. waits
  on all commrequests).  Only waits if V->HaloInTransit flag is True;
  sets this flag to False afterwards.  Returns immediately of not
  running in parallel.  If running in parallel but V->ParallelFlag ==
  False, errors out.

  The "reverse" exchange means that halos send their data to regular
  elements (values are added).

INPUTS:
 
  V   : vector for which to end the halo exchange 

OUTPUTS: 

  None

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Definition: xf_BcastVectorBasicInfo
extern int 
xf_BcastVectorBasicInfo(xf_Vector *V_Glob, xf_Vector *V, 
                        enum xfe_Bool DestroyGlob);
/*
 PURPOSE:
 
 Broadcasts the basic information in V_Glob residing in proc0 
 to V residing in all processors
 
 INPUTS:
 
 V_Glob : vector to broacast (relevant only on proc 0)
 DestroyGlob: if True, the basic info pointers in V_Glob 
              are set to NULL
 
 OUTPUTS: 
 
 V : vector residing on all processors
 Note that some information in V_Glob gets destroyed
 (look at source code)
 
 RETURN:
 
 Error Code
 */


/******************************************************************/
//   FUNCTION Prototype: xf_ParallelizeVector
extern int 
xf_ParallelizeVector( xf_Mesh *Mesh, xf_Vector *V_Glob, xf_Vector *V);
/*
PURPOSE:

  Parallelizes the Vector structure V_Glob: each processor receives a
  relevant local copy V.  The top-level memory for V must be
  allocated, e.g. via a call to CreateVector.  The parallelization is
  destructive to V_Glob.  In the serial case, V inherits all data from
  V_Glob, and the data in V_Glob is initialized to null.

INPUTS:

  Mesh : pointer to Mesh on all procs
  V_Glob : vector to parallelize (relevant only on proc 0)

OUTPUTS:

  V : the parallelized vector structure on all procs

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ParallelizeDataSet
extern int 
xf_ParallelizeDataSet( xf_All *All, xf_DataSet *DataSet_Glob,
		       xf_DataSet *DataSet);
/*
PURPOSE:

  Parallelizes the DataSet structure DataSet_Glob.  Each Data node in
  the linked list is parallelized.  All on each proc is required for
  the parallelization.

INPUTS:

  All : pointer to All on each proc. 
  DataSet_Glob : the global DataSet structure, residing on proc 0

OUTPUTS:

  DataSet : the parallelized DataSet structure on all procs

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_UnParallelizeVector
extern int 
xf_UnParallelizeVector( xf_Mesh *Mesh, xf_Vector *V, xf_Vector *V_Glob);
/*
PURPOSE:

  UnParallelizes the Vector structure V from each processor to obtain
  a global V_Glob.  V_Glob must be created.

INPUTS:

  Mesh : pointer to Mesh on all procs
  V : the parallelized vector structure on all procs

OUTPUTS:

  V_Glob : unparallelized vector (relevant only on proc 0)


RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_UnParallelizeDataSet
extern int 
xf_UnParallelizeDataSet( xf_All *All, xf_DataSet *DataSet, 
			 xf_DataSet *DataSet_Glob);
/*
PURPOSE:

  UnParallelizes the DataSet structure that resides in All->DataSet on
  all proc.  Every writeflag=true data node in the linked list is un
  parallelized.  The assembled data set is stored in DataSet_Glob on
  proc 0.

INPUTS:
 
  All : pieces of All on each proc. Relevant on all procs
  DataSet : the parallelized DataSet structure on all procs

OUTPUTS: 

  DataSet_Glob : the global DataSet structure, residing on proc 0

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_ProcViz
extern int
xf_ProcViz(xf_All *All);
/*
PURPOSE:

  Stores processor number in a real element vector, ProcID, which is
  made writable.  Useful for visualizing partitions.

INPUTS:

  All : all structure

OUTPUTS: 

  None

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ParallelSumNodes
extern int 
xf_ParallelSumNodes( xf_Mesh *Mesh, void *NodeData, int rank, enum xfe_SizeType Size);
/*
PURPOSE:

  Nodal data, NodeData, on each proc is set to the sum of the values
  from the same global node over all procs.  Currently this involves a
  global send to the root proc, and hence this function should not be
  called to often.

INPUTS:

  Mesh : mesh structure
  NodeData : nodal data
  rank : # data items per node
  Size : xfe_SizeInt or xfe_SizeReal

OUTPUTS: 

  None, NodeData is modified

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_BcastDataSet
extern int 
xf_BcastDataSet(xf_DataSet **pDataSet);
/*
 PURPOSE:
 
 Creates a copy of the DataSet data structure of the root processor 
 on the other processors. 
 
 INPUTS:
 
 pDataSet : Pointer to the DataSet structure on each processor.
 
 OUTPUTS: 
 
 None. pDataSet gets modified on non-root processors.
 
 RETURN:
 
 Error Code
 */

#endif // end ifndef _xf_DataParallel_h
