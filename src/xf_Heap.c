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
 FILE:  xf_Heap.c
 
 This file contains functions for working with heaps (e.g. a min-heap).
 
*/

#include "xf.h"
#include "xf_Memory.h"
#include "xf_HeapStruct.h"

/******************************************************************/
//   FUNCTION Definition: xf_MinHeapCascadeDown
void
xf_MinHeapCascadeDown(xf_Heap *Heap, int i)
{  
  int lchild, rchild, minchild;
  xf_HeapNode *H, Htemp;

  H = Heap->Node;
  Htemp = H[i]; // copy over node in question into Htemp
  while ( (minchild=lchild=2*i+1) < Heap->N ){
    if ( (rchild=lchild+1) < Heap->N)
      minchild = (H[rchild].rval < H[lchild].rval) ? rchild : lchild; // minimum
    if (H[minchild].rval < Htemp.rval){ // swap needed
      H[i] = H[minchild];
      Heap->ipos2Node[H[minchild].ipos] = i;
      i = minchild;
    }
    else break; // done
  } // end while i

  // place node in question into correct position
  H[i] = Htemp;
  Heap->ipos2Node[Htemp.ipos] = i;
}

/******************************************************************/
//   FUNCTION Definition: xf_MinHeapCascadeUp
void
xf_MinHeapCascadeUp(xf_Heap *Heap, int i)
{  
  int parent;
  xf_HeapNode *H, Htemp;

  H = Heap->Node;
  Htemp = H[i]; // copy over node in question into Htemp
  while ( i > 0){
    parent = (i-1)/2;
    if (H[parent].rval > Htemp.rval){ // parent needs to move down
      H[i] = H[parent];
      Heap->ipos2Node[H[parent].ipos] = i;
      i = parent;
    }
    else break; // done
  } // end while i

  // place node in question into correct position
  H[i] = Htemp;
  Heap->ipos2Node[Htemp.ipos] = i;
}


/******************************************************************/
//   FUNCTION Definition: xf_BuildMinHeap
int
xf_BuildMinHeap(xf_HeapNode *Node, int N, xf_Heap **pHeap)
{  
  int ierr, i;
  xf_Heap *Heap;

  // Allocate Heap structure
  ierr = xf_Error(xf_Alloc( (void **) pHeap, 1, sizeof(xf_Heap)));
  if (ierr != xf_OK) return ierr;
  Heap = (*pHeap);

  // point into nodes
  Heap->Node = Node;

  // set N and N0 = size of heap
  Heap->N = Heap->N0 = N;

  // allocate and fill in the initial index array
  ierr = xf_Error(xf_Alloc( (void **) &Heap->ipos2Node, Heap->N0, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<Heap->N0; i++) Heap->ipos2Node[i] = Heap->Node[i].ipos = i;

  // put all nodes in the correct place, starting from the end
  for (i=(Heap->N-1)/2; i>=0; i--) xf_MinHeapCascadeDown(Heap, i);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MinHeapTakeAwayRoot
void
xf_MinHeapTakeAwayRoot(xf_Heap *Heap, xf_HeapNode **pH)
{  
  int N;
  xf_HeapNode *H, Htemp;

  N = Heap->N;
  H = Heap->Node;

  // swap root with last position in heap
  Heap->ipos2Node[H[0  ].ipos] = N-1;
  Heap->ipos2Node[H[N-1].ipos] = 0;
  swap(H[0], H[N-1], Htemp);

  // heap becomes smaller
  Heap->N--;

  // cascade root down
  xf_MinHeapCascadeDown(Heap, 0);
  
  // point output to former last position in heap
  (*pH) = H+N-1;

}

/******************************************************************/
//   FUNCTION Definition: xf_MinHeapChangeNode
void
xf_MinHeapChangeNode(xf_Heap *Heap, int ipos, real rval)
{  
  int i;
  xf_HeapNode *H;

  H = Heap->Node;

  // get heap node
  i = Heap->ipos2Node[ipos];

  // set new real value at the node in question
  H[i].rval = rval;

  // nothing else to do if node is no longer in heap
  if (i >= Heap->N) return;

  if (rval > H[i].rval) // heap below i might be affected
    xf_MinHeapCascadeDown(Heap, i);
  else                  // heap above i might be affected
    xf_MinHeapCascadeUp(Heap, i);

}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyHeap
void
xf_DestroyHeap(xf_Heap *Heap)
{  

  // destroy Nodes and position index
  xf_Release( (void *) Heap->Node);
  xf_Release( (void *) Heap->ipos2Node);

  // also destroy self
  xf_Release( (void *) Heap);

}



#if( UNIT_TEST==1 )
#include "xf_Heap.test.in"
#endif





