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

#ifndef _xf_Heap_h
#define _xf_Heap_h 1

/*
 FILE:  xf_Heap.h
 
 This file contains prototypes for heap functions.
 
*/

#include "xf_HeapStruct.h"


/******************************************************************/
//   FUNCTION Prototype: xf_BuildMinHeap
extern int
xf_BuildMinHeap(xf_HeapNode *Node, int N, xf_Heap **pHeap);
/*
PURPOSE:

  Builds a min-heap out of N heap nodes.  The heap points directly
  into the Nodes, which are assumed contiguous in memory.  Also, the
  Heap structure assumes responsibility for destroying the heap nodes
  via xf_DestroyHeap.
  
INPUTS:

  Node  : vector of heap nodes
  N     : number of nodes

OUTPUTS: 

  Heap  : built heap structure

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MinHeapTakeAwayRoot
extern void
xf_MinHeapTakeAwayRoot(xf_Heap *Heap, xf_HeapNode **pNode);
/*
PURPOSE:

  Takes away root from Heap and adjusts its structure to keep it a
  min-heap.  On return, (*pNode) points to what was the root node.
  
INPUTS:

  Heap  : heap structure

OUTPUTS: 
  
  (*pNode) : pointer to former root node (no longer root of modified Heap).

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_MinHeapChangeNode
extern void
xf_MinHeapChangeNode(xf_Heap *Heap, int ipos, real rval);
/*
PURPOSE:

  Changes rval of heap node with integer value of ipos (original
  position).  Adjusts heap to preserve min-heap property.
  
INPUTS:

  Heap  : heap structure
  ipos  : integer value identifying node that needs to be modified
  rval  : new real value given to the node in question

OUTPUTS: 

  None, Heap is modified

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyHeap
extern void
xf_DestroyHeap(xf_Heap *Heap);
/*
PURPOSE:

  Destroys Heap structure
  
INPUTS:

  Heap  : heap structure to destroy

OUTPUTS: 

  None, Heap is destroyed

RETURN:

  Error Code
*/

#endif // end ifndef _xf_Heap_h
