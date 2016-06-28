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

#ifndef _xf_HeapStruct_h
#define _xf_HeapStruct_h 1

/*
 FILE:  xf_HeapStruct.h
 
 This file defines structures for working with heaps.
 
*/


/* Structure for individual heap nodes */
typedef struct
{
  real rval;        // real value used for sorting
  int  ipos;        // integer value storing original position in heap
}
xf_HeapNode;


/* Structure for entire heap */
typedef struct
{
  xf_HeapNode *Node;  // heap nodes
  int N;              // current size of heap
  int N0;             // max size of heap (number of allocated nodes)
  int *ipos2Node;     // ipos2Node[ival] = current heap node of original position ipos
}
xf_Heap;

#endif // end ifndef _xf_HeapStruct_h



