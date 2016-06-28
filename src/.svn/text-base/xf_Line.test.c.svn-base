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



TEST_xf_TreePushPop()
{
  int ierr, i, k, root = -1;
  int **Tree;
  int TT0[18] = {-1,1,2, 0,-1,-1, 0,3,4, 2,-1,5, 2,-1,-1, 3,-1,-1};
  int TT1[18] = {-1,1,3, 0,-1,-1, -1,-1,-1, 0,-1,5, 5,-1,-1, 3,-1,4};
  int TT2[18] = {-1,1,3, 0,-1,-1, -1,-1,-1, 0,-1,4, 3,-1,-1, -1,-1,-1};
  int TT3[18] = {-1,1,3, 0,-1,-1, -1,-1,-1, 0,-1,-1, -1,-1,-1, -1,-1,-1};
  int TT4[18] = {-1,-1,-1, -1,-1,3, -1,-1,-1, 1,-1,-1, -1,-1,-1, -1,-1,-1};
  int TT5[18] = {-1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,-1};
  real val[6] = {5,3,10,6,17,9};
  real g[3] = {3,13,6}; // rhs
  real y[3] = {1,-3,-2}; // solution
  
  // alloc Tree
  ierr = xf_Error(xf_Alloc2( (void ***) &Tree, 6, 3, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<3*6; k++) Tree[0][k] = -1;
  
  // push values onto tree
  for (i=0; i<6; i++){
    ierr = xf_Error(xf_TreePush( i, val, &root, Tree));
    if (ierr != xf_OK) return ierr;
  }
  xf_AssertIntVectorEqual(Tree[0], TT0, 18);

  // pop node 2
  ierr = xf_Error(xf_TreePop( 2, val, &root, Tree));
  if (ierr != xf_OK) return ierr;
  xf_AssertIntVectorEqual(Tree[0], TT1, 18);
  
  // pop node 5
  ierr = xf_Error(xf_TreePop( 5, val, &root, Tree));
  if (ierr != xf_OK) return ierr;
  xf_AssertIntVectorEqual(Tree[0], TT2, 18);

  // pop node 4
  ierr = xf_Error(xf_TreePop( 4, val, &root, Tree));
  if (ierr != xf_OK) return ierr;
  xf_AssertIntVectorEqual(Tree[0], TT3, 18);

  // pop node 0
  ierr = xf_Error(xf_TreePop( 0, val, &root, Tree));
  if (ierr != xf_OK) return ierr;
  xf_AssertIntVectorEqual(Tree[0], TT4, 18);
  xf_AssertEqual(root, 1);

  // pop node 1
  ierr = xf_Error(xf_TreePop( 1, val, &root, Tree));
  if (ierr != xf_OK) return ierr;
  xf_AssertIntVectorEqual(Tree[0], TT5, 18);
  xf_AssertEqual(root, 3);

  // pop node 3
  ierr = xf_Error(xf_TreePop( 1, val, &root, Tree));
  if (ierr != xf_OK) return ierr;
  xf_AssertIntVectorEqual(Tree[0], TT5, 18);
  xf_AssertEqual(root, -1);

  // release Tree
  xf_Release2((void **)  Tree);

  return xf_OK;  
}
