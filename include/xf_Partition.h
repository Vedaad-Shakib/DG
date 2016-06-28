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

#ifndef _xf_Partition_h
#define _xf_Partition_h 1

/*
  FILE:  xf_Partition.h

  This file contains headers for partition functions.  In the serial
  case, these functions are defined in xf_Serial.c

*/


/******************************************************************/
//   FUNCTION Prototype: xf_PartitionMesh
extern int 
xf_PartitionMesh( xf_Mesh *Mesh, int ***pElem2Proc, int **ElemWeight, 
                 int *ConnectWeight);
/*
PURPOSE:

  Determines partitioning for the input Mesh.  Does not perform any
  communication, but stores the proc number of each elem in
  (*pElem2Proc).

INPUTS:

  Mesh : mesh to partition
  ElemWeight : If not NULL, the element weights will be used in 
               the mesh partitioning
  ConnectWeight : If not NULL, the interelement connection (iface) 
                  weight will be used in the mesh partitioning.

OUTPUTS:

  (*pElem2Proc): (*pElem2Proc)[egrp][elem] is the proc number of egrp,elem

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Definition: xf_CPUWorkReport
extern int
xf_CPUWorkReport(xf_All *All, xf_Vector *State, enum xfe_Bool *load_imbalance);
/*
 PURPOSE:
 
 Write a reportof representative numbers of CPU work 
 in MATLAB format
 
 INPUTS:
 
 All : All structure containing mesh and data
 pState : pointer to the state vector
 
 OUTPUTS:
 
 A MATLAB file named [SavePrefix]_[myRank]of[nProc].m is written
 
 RETURN:
 
 Error code
 */

/******************************************************************/
//   FUNCTION Definition: xf_CPULoadBalance
extern int
xf_CPULoadBalance(xf_All *All, xf_Vector **pState, int nState);
/*
 PURPOSE:
 
 Repartitions the mesh in order to improve the balance of CPU work 
 load and improve the effectivity of line-based preconditioners
 
 INPUTS:
 
 All : All structure containing mesh and data
 pState : pointer to the state vector
 nState : number of state vectors (ex: multi-step time integration)
 
 OUTPUTS:
 
 None: All->Mesh gets modified for all processors and the 
       pState pointers get reset
 
 RETURN:
 
 Error code
 */

extern int
xfYu_CPULoadBalance(xf_All *All, xf_Vector **pState, int nState);

#endif // end ifndef _xf_Partition_h

