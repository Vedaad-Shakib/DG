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

#ifndef _xf_LinQueue_h
#define _xf_LinQueue_h 1


/*
  FILE:  xf_LinQueue.h

  This file contains prototypes for xf_LinQueue.c

*/

#include "xf_LinQueueStruct.h"


/******************************************************************/
//   FUNCTION Prototype: xf_InitLinQueue
extern int
xf_InitLinQueue(xf_LinQueueData *LinQ);
/*
PURPOSE:
 
  Initializes linearization queue, LinQ (sets quantities to NULL)

INPUTS:

  LinQ : linearization queue data

OUTPUTS:  None

RETURNS:  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyLinQueue
extern void
xf_DestroyLinQueue(xf_LinQueueData *LinQ);
/*
PURPOSE:
 
  Destroys linearization queue.  Note, self is not destroyed.

INPUTS:

  LinQ : linearization queue data for which data is released.

OUTPUTS:  None

RETURNS:  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ClearLinQueue
extern void
xf_ClearLinQueue(xf_LinQueueData *LinQ);
/*
PURPOSE:
 
  Sets all terms in linearization queue to inactive.  Does not
  de-allocate.

INPUTS:

  LinQ : linearization queue data

OUTPUTS:  None

RETURNS:  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_AddToLinQueue
extern int
xf_AddToLinQueue(real *F_u, enum xfe_LinQTermType Term,
		 int nq, int dim, int sr2, int idim, real *w,
		 int dw, real fac, xf_LinQueueData *LinQ);
/*
PURPOSE:
 
  Adds a flux (gradient) linearization to a linearization queue.

INPUTS:

  F_u  : flux (gradient) linearization, nq*sr2[*dim]
  Term : type indicating what (Phi) F_u is sandwiched in-between
  nq   : number of quadrature points
  dim  : spatial dimension
  sr2  : state rank squared
  idim : if >= 0, indicates what component of flux gradient is affected
  w    : weight at each quadrature point applied during summation
  dw   : spacing in memory of weights
  fac  : constant factor included in incrementation of LinQ

OUTPUTS:

  LinQ : linearization queue data that is incremented

RETURNS:  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ApplyLinQueue
extern int
xf_ApplyLinQueue(xf_LinQueueData *LinQ, xf_BasisData *RowPhiData,
		 xf_BasisData *ColPhiData, int sr2, real *R_U);
/*
PURPOSE:
 
  Applies terms in LinQ to R_U.

INPUTS:

  LinQ : linearization queue data
  RowPhiData : basis functions/gradients for the row dimension
  ColPhiData : basis functions/gradients for the column dimension
  sr2  : state rank squared

OUTPUTS:

  R_U : residual linearization (added to)

RETURNS:  Error code

*/

#endif // end ifndef _xf_LinQueue_h
