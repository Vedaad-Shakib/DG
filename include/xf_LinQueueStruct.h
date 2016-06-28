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

#ifndef _xf_LinQueueStruct_h
#define _xf_LinQueueStruct_h 1


/*
  FILE:  xf_LinQueueStruct.h

  This file contains structures for residual linearization

*/

/* Supported linearization terms */
enum xfe_LinQTermType{
  xfe_LinQTerm_PhiPhi,    // Phi^T  *  F_u *  Phi
  xfe_LinQTerm_GPhiPhi,   // GPhi^T *  F_u *  Phi
  xfe_LinQTerm_PhiGPhi,   //  Phi^T *  F_u * GPhi
  xfe_LinQTerm_GPhiGPhi,  // GPhi^T *  F_u * GPhi
  xfe_LinQTerm_Last
};


/* Linearization queue data */
typedef struct
{
  real *F_u[xfe_LinQTerm_Last];              // stores the matrices
  enum xfe_Bool Active[xfe_LinQTerm_Last];   // indicates whether term is active
  int Size[xfe_LinQTerm_Last];               // indicates matrix sizes (for persistency)
  real *T;                                   // temporary storage
  int Tsize;                                 // temporary storage size
}
xf_LinQueueData;



#endif // end ifndef _xf_LinQueue_h
