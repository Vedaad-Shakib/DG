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

#ifndef _xf_Line_h
#define _xf_Line_h 1

/*
  FILE:  xf_Line.h

  This file contains the headers for functions in xf_Line.c

*/

#include "xf_LinearSolverStruct.h"



/******************************************************************/
//   FUNCTION Prototype: xf_DestroyLineSet
extern int 
xf_DestroyLineSet( xf_LineSet *LineSet);
/*
PURPOSE:

  Destroys a LineSet structure

INPUTS:

  LineSet : pointer to line set

OUTPUTS: 

  None.  LineSet is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_PreconditionerLineCheck
extern void
xf_PreconditionerLineCheck(enum xfe_PreconditionerType Preconditioner,
			   enum xfe_Bool *CRequired, enum xfe_Bool *SortLines);
/*
PURPOSE:

  Determines whether Preconditioner requires lines (and hence a
  connectivity vector), and whether the lines need to be sorted
  (e.g. for Line Gauss Seidel).

INPUTS:
  
  Preconditioner: type of preconditioner in question

OUTPUTS:

  (*CRequired): True if connectivity and lines are required
  (*SortLines): True if lines need to be sorted

RETURNS:

  None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_FindLineConnectivity
extern int
xf_FindLineConnectivity(xf_All *All, xf_Vector **pC);
/*
PURPOSE:

  Looks for a line connectivity vector in All->DataSet.  Creates
  one if not found.  Stores vector in C

INPUTS:

  All : All structure

OUTPUTS:

  (*pC) : elem-to-elem connectivity vector

RETURNS:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CreateLines
extern int
xf_CreateLines(xf_All *All, xf_Vector *C, xf_JacobianMatrix *R_U, 
               xf_LineSet **pLineSet);
/*
PURPOSE:
 
  Computes set of element lines according to elem-to-elem connectivity
  vector, C.  Stores this LineSet in R_U->LineSet.  If R is provided,
  lines are ordered from maximum residual to minimum residual norm (of
  elements contained within the lines).


INPUTS:

  All : All structure
  C   : elem-to-elem connectivity vector
  R_U : Jacobian matrix
  pLineSet : Pointer to the created lineset

OUTPUTS:

  None, R_U->LineSet is created/filled.

RETURNS:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_SortLines
extern int
xf_SortLines(xf_All *All, xf_Vector *C, xf_JacobianMatrix *R_U, xf_Vector *R);
/*
PURPOSE:

  Sorts lines from maximum R to minimum R norm (of elements contained
  within the lines).  Uses L1 norm of R on each element.  The sorting
  favors neighboring lines because as each line is taken into the
  ordering, its L1(R) is distributed among its neighbor lines, in
  percentages proportional to the line-to-line connectivity.


INPUTS:

  All : All structure
  C   : elem-to-elem connectivity vector
  R_U : Jacobian matrix
  R   : real vector (e.g. residual)

OUTPUTS:

  None

RETURNS:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeLineSet
extern int
xf_UnParallelizeLineSet(xf_Mesh *Mesh_Glob, xf_Vector *C_Glob, 
                        xf_LineSet **LineSet_Glob, xf_Mesh *Mesh,
                        xf_LineSet *LineSet);
/*
 PURPOSE:
 
 Gets LineSet from all processors and put it in LineSet_Glob 
 and connects the lines.
 
 
 INPUTS:
 
 Mesh_Glob: Global mesh (only used by proc 0)
 C_Glob: Global element coupling vector (only used by proc 0)
 LineSet_Glob: Global lineset. No need to pre-allocate.
 Mesh: Local mesh in each processor.
 LineSet: Local lineset in each processor.
  
 OUTPUTS:
 
 None, LineSet_Glob gets created.
 
 RETURNS:
 
 Error Code
 */


#endif // end ifndef _xf_Line_h
