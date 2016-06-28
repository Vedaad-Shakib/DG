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

#ifndef _xf_LineStruct_h
#define _xf_LineStruct_h 1

/*
  FILE:  xf_LineStruct.h

  This file contains the xflow Line structures

*/


/* Line structure.  The face vector stores, for each element in the
   line, the local face number across which the next element in the
   line resides. */
typedef struct
{
  int nelem;  // number of elements on this line
  int *egrp;  // vector of element group numbers
  int *elem;  // vector of element numbers
  int *face;  // vector of face numbers
}
xf_Line;


/* LineSet structure */
typedef struct
{
  int nLine;       // Number of lines
  xf_Line *Line;   // pointers to Line structures
  int **Elem2Line; // Elem2Line[egrp][elem] = line # of element
  int ***Face2M;   // 1 if elem across face is part of preconditioner
  int negrp;       // number of elem groups (for alloc/destroy purposes)
}
xf_LineSet;


/* Used for enumerating entries in a static binary search tree */ 
enum xfe_TreeNode { 
  xfe_Parent, 
  xfe_ChildL, 
  xfe_ChildR, 
  xfe_TreeLast
};


#endif // end ifndef _xf_LineStruct_h
