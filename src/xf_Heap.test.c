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


#include "xf_Unit.h"


TEST_xf_MinHeap()
{
  int ierr, i;
  int N = 11;
  //               0    1    2    3    4    5    6    7    8    9   10
  real rval[] = {1.7, 0.3, 4.2, 0.7, 0.1, 0.8, 2.7, 1.3, 1.31, 1.4, 0.2};
  real rmin[] = {0.01, 0.1, 0.2, 0.3, 0.7, 0.8, 1.3, 1.31, 1.4, 1.7, 4.2};
  int  ipos[] = { 6, 4, 10, 1, 3, 5, 7, 8, 9, 0, 2};
  xf_HeapNode *H;
  xf_Heap *Heap;

  // allocate Heap node data structure (will get destroyed with heap)
  ierr = xf_Error(xf_Alloc( (void **) &H, N, sizeof(xf_HeapNode)));
  xf_AssertEqual(ierr, xf_OK); 
  
  // copy over data into heap nodes
  for (i=0; i<N; i++) H[i].rval = rval[i];

  // build a heap
  ierr = xf_Error(xf_BuildMinHeap(H, N, &Heap));
  xf_AssertEqual(ierr, xf_OK); 

  // heap size should be N
  xf_AssertEqual(Heap->N, N);

  // change a node (this will be the new minimum)
  xf_MinHeapChangeNode(Heap, 6, 0.01);

  // pull off and check minimum multiple times
  for (i=0; i<N; i++){
    xf_MinHeapTakeAwayRoot(Heap, &H);
    xf_AssertWithin(H->rval, rmin[i], UTOL0);
    xf_AssertEqual(H->ipos, ipos[i]);
  } // i
  
  // nothing should be left in the heap
  xf_AssertEqual(Heap->N, 0);

  // destroy heap
  xf_DestroyHeap(Heap);

  return xf_OK;
}
