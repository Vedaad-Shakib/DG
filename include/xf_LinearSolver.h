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

#ifndef _xf_LinearSolver_h
#define _xf_LinearSolver_h 1

/*
  FILE:  xf_LinearSolver.h

  This file contains the headers for the linear solver functions

*/

#include "xf_LinearSolverStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_Jacobian_n
extern int
xf_Jacobian_n(xf_JacobianMatrix *R_U, int egrp, int elem);
/*
PURPOSE:

  Pulls of n (# unknowns) for egrp,elem, allowing for variable orders

INPUTS:

  R_U : Jacobian matrix
  egrp : element group (less than zero allowed, returns 0)
  elem : element number (optional, less than zero means use nvec)
  
OUTPUTS: None

RETURNS:

  number of unknowns stored in R_U->nvec, or R_U->vnvec[egrp][elem]
  or zero if egrp < 0

*/  

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyILUData
extern int
xf_DestroyILUData(xf_LinearSolverILUData *ILUData);
/*
PURPOSE:

  Destroys ILUData structure (including self)

INPUTS:

  ILUData : structure to destroy

OUTPUTS:

RETURNS:

  Error Code

*/ 



/******************************************************************/
//   FUNCTION Prototype: xf_Jacobian_Mult
extern int
xf_Jacobian_Mult(xf_All *All, xf_JacobianMatrix *R_U,
		 xf_Vector *X, enum xfe_AddType AddFlag, 
		 enum xfe_Bool TransposeFlag, 
		 xf_SolverData *SolverData, xf_Vector *Y);
/*
PURPOSE:

  Sets 
         Y @= R_U*X

  where @= is one of {=, +=, -=, =-} according to AddFlag.

INPUTS:

  All : All structure
  R_U : Jacobian matrix (contains preconditioner info, if any)
  X   : Input that gets multiplied by R_U 
  AddFlag: one of xfe_Set, xfe_Add, etc.
  TransposeFlag : if True, R_U^T is used instead of R_U
  SolverData : solver data structure

OUTPUTS:

  Y   : Resulting vector.

RETURNS:

  Error Code

*/  

/******************************************************************/
//   FUNCTION Prototype: xf_SolveLinearSystem
extern int
xf_SolveLinearSystem(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *R, 
		     enum xfe_Bool TransposeFlag, int nIter,
		     xf_SolverData *SolverData, xf_Vector *dU );
/*
PURPOSE:

  Solves:
                 R_U   (dU) + R = 0       (TransposeFlag == False)
                 R_U^T (dU) + R = 0       (TransposeFlag == True )

  using the specified LinearSolver (in All)

INPUTS:

  All : All structure
  R_U : Jacobian matrix
  R   : residual vector (see above)
  TransposeFlag: if true, R_U^T is used instead
  nIter : if > 0, requested number of linear iterations.  Will overwrite
          any linear-solver-specific number in All's parameters.
  SolverData : solver data structure

OUTPUTS:

  dU : linear-system solution vector

RETURNS:

   Error Code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_CalculateLinearResidual
extern int
xf_CalculateLinearResidual(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X, 
			   xf_Vector *R, enum xfe_Bool TransposeFlag, 
			   xf_SolverData *SolverData, xf_Vector *LinR);
/*
PURPOSE:

  Computes the linear residual:

  LinR =     R_U   (X) + R        (TransposeFlag == False)
             R_U^T (X) + R        (TransposeFlag == True )

  R is optional -- if not passed in, then the calculation is:

  LinR += R_U * X

INPUTS:

  All : All structure
  R_U : Jacobian matrix
  X   : vector that is multiplied by R_U
  R   : constant lhs vector (optional)
  TransposeFlag: if true, R_U^T is used instead
  SolverData : solver data structure

OUTPUTS:

  LinR : linear-system residual

RETURNS:

   Error Code

*/





/******************************************************************/
//   FUNCTION Prototype: xf_LinearIterCG
extern int
xf_LinearIterCG(xf_All *All, xf_VectorSet *VS, real tol, enum xfe_Verbosity Verbosity,
		enum xfe_Bool PreconditionFlag, enum xfe_Bool ZeroFlag,
		enum xfe_Bool SavePoint, enum xfe_Bool CGRestart,
		xf_Vector *R, xf_Vector *U, xf_Vector **pV, xf_Vector **pW, 
		xf_LinearSolverData **pLinearSolverData, 
		enum xfe_LinearStatusType *Status);
/*
PURPOSE:

  Solves:   Operator * U + R = 0

  using the Conjugate Gradient method.  "Operator", which must be
  symmetric, is a user-implemented matrix-vector multiplication,
  performed via a call-back.  That is, this function should be called
  multiple times and Status should be checked after each call.  If a
  multiplication is requested, W = Operator * V should be performed.
  The final solution is stored in U.

  Conjugate Gradient algorithm (A = Operator):

  U = 0
  p = r = R
  rnorm2 = <r,r>
  while (sqrt(rnorm2) > tol),
    w = A*p
    alpha = rnorm2/<w,p>
    U = U - alpha*p
    r = r - alpha*w
    rnorm2_new = <r,r>
    p = r + rnorm2_new/rnorm2*p
    rnorm2 = rnorm2_new
  end while


  Preconditioned Conjugate Gradient algorithm
  (A = Operator, M^{-1} = preconditioner)

  U = 0
  r = R
  z = M^{-1}*r
  p = z
  rnorm2 = <z,r>
  while (sqrt(rnorm2) > tol),
    w = A*p
    alpha = rnorm2/<w,p>
    U = U - alpha*p
    r = r - alpha*w
    z = M^{-1}*r
    rnorm2_new = <z,r>
    p = z + rnorm2_new/rnorm2*p
    rnorm2 = rnorm2_new
  end while
  

INPUTS:
  All   : All structure. Optional. Only necessary if reading a restart
          CG file or writing savepoints.

  VS    : VectorSet containing at least 3 allocated vectors.  Used for
          temporary storage.
  tol   : desired relative residual tolerance
  Verbosity : as defined in xf.h; this is the level of output printed
              to the screen.  High prints out a lot of debug information.

  PreconditionFlag : if True, the preconditioned version is used.  The calling
                     function must be prepared to apply M^{-1} when requested.

  ZeroFlag : if False, U is not set to zero at the beginning of the algorithm.
             Rather, the provided U serves as the initial guess for the iterations.

  SavePoint : if True, CG.data and CG.txt save points will be
              written every CG iteration
	      
  CGRestart: if True, CG.data and CG.txt will be read in for
	     a restart.  Only checked on first call.

  R : Input residual vector (see description above)
	      

OUTPUTS:

  U : upon convergence, contains the solution.  This vector must be
      allocated before the call.

  EV : If not NULL, upon convergence, contains eigenvectors.  This
       vectorset must be allocated before the call.  NULL indicates
       that eigenvector computation is not requested.

  V, W : if Status == xfe_Eig_Multiply, the operation W = Operator*V
         is requested. The function must then be called again for
         further iterations.  V and W are only pointers to vectors and
         need not be allocated before the call.

  LinearSolverData : internal storage required by this function
                     through subsequent calls.  Calling functions do
                     not have to allocate or release this -- just pass
                     in the same pointer every time.

  Status : Current status of the CG solution:
           xfe_LinearMultiply : means the operation W = A*V is requested
           xfe_LinearPrecondition : means the operation W = M^{-1}*V is requested
	   xfe_LinearConverged : system has converged to tolerance (tol)

RETURNS:

   Error Code

*/






#endif // end ifndef _xf_LinearSolver_h
