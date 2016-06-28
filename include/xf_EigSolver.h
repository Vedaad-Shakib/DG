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

#ifndef _xf_EigSolver_h
#define _xf_EigSolver_h 1

/*
  FILE:  xf_EigSolver.h

  This file contains header for eigenvalue solver routines

*/

#include "xf_EigSolverStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_EigIterLanczos
extern int
xf_EigIterLanczos(xf_All *All, int nev, int ncv, xf_VectorSet *VS, real tol, 
		  enum xfe_Verbosity Verbosity, int nev_restart, 
		  enum xfe_Bool SavePoint, enum xfe_Bool LanczosRestart, 
		  real *E, xf_VectorSet *EV, xf_Vector **pV, xf_Vector **pW, 
		  xf_EigSolverData **pEigSolverData, 
		  enum xfe_EigStatusType *Status);
/*
PURPOSE:

  Performs an iteration of the Lanczos method for finding eigenvalues
  of a real, symmetric matrix.  This function does not know about the
  matrix, and returns control to the calling function (with an
  appropriately-set Status) whenever it requires a matrix-vector
  multiplication.  The calling function must then compute W = A*V, and
  call this function again.  When convergence is detected, Status is
  set appropriately.

  The algorithm used here is the implicitly-restarted Lanczos method,
  as described in the ARPACK documentation.  The top-level ARPACK
  function on which this function is modeled is dsaupd.f


INPUTS:

  All   : All structure. Optional. Only necessary if reading a restart
          Lanczos file or writing savepoints.

  nev   : number of desired eigenvalues and eigenvectors

  ncv   : number of available vectors in VS (must be > nev), not
          including a residual vector, which is the last one in VS.  A
          good practice is to use about 2*nev for ncv.

  VS    : VectorSet containing ncv+1 allocated vectors
  tol   : desired relative tolerance for eigenvalues (e.g. .01 is 1%)
  Verbosity : as defined in xf.h; this is the level of output printed
              to the screen.  High prints out a lot of debug information.

  nev_restart: if > 0, a restart is assumed in which the first
               nev_restart vectors in VS contain existing
               eigenvectors, and the corresponding eigenvalues are in
               E.  Only allowed on first call to this function.
  
  SavePoint : if True, Lanczos.data and Lanczos.txt save points will be
              written every Lanczos iteration

  LanczosRestart: if True, Lanczos.data and Lanczos.txt will be read in for
                  a restart.  Only checked on first call.

OUTPUTS:

  E : upon convergence, contains eigenvalues.  This vector must be
      allocated before the call.

  EV : If not NULL, upon convergence, contains eigenvectors.  This
       vectorset must be allocated before the call.  NULL indicates
       that eigenvector computation is not requested.

  V, W : if Status == xfe_Eig_Multiply, the operation W = A*V is
         requested. The function must then be called again for further
         iterations.  V and W are only pointers to vectors and need
	 not be allocated before the call.

  EigSolverData : internal storage required by this function through
                  subsequent calls.  Calling functions do not have to
                  allocate or release this -- just pass in the same
                  pointer every time.

  Status : Current status of the Lanczos factorization:
           xfe_EigMultiply : means the operation W = A*V is requested
	   xfe_EigConverged : eigs have converged to tolerance (tol)

RETURN:

  Error Code

*/

#endif // end ifndef _xf_EigSolver_h
