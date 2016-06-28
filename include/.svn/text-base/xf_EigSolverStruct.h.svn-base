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

#ifndef _xf_EigSolverStruct_h
#define _xf_EigSolverStruct_h 1

/*
  FILE:  xf_EigSolverStruct.h

  This file contains the structures for the eigenvalue solver

*/

#include "xf.h"


/* Steady or unsteady calculations dictated by this type */
enum xfe_EigStatusType{
  xfe_EigMultiply,
  xfe_EigConverged,
  xfe_EigForceExit,
  xfe_EigStatusLast,
};

/* corresponding names */
static char *xfe_EigStatusName[xfe_EigStatusLast] = {
  "Multiply",
  "Converged",
  "ForceExit",
};



/* EigSolverData structure (not read or written) */
typedef struct
{
  int  iIter;             // current iteration number
  int  nRestart;          // number of implicit restarts taken
  
  int  nInvariant;        // number of times an invariant subspace was found
  enum xfe_Bool InvariantSubspaceFound;

  int ncv;
  real *H; /* Symmetric tridiagonal matrix. Subdiagonal stored from 1
	      to ncv-1. Main diagonal stored from ncv to 2*ncv-1.*/

}
xf_EigSolverData;
 

#endif // end ifndef _xf_EigSolverStruct_h
