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

#ifndef _xf_LeanSolver_h
#define _xf_LeanSolver_h 1

/*
  FILE:  xf_LeanSolver.h

  This file contains the headers for the memory-lean solver functions

*/

#include "xf_SolverStruct.h"
#include "xf_SolverTools.h"



/******************************************************************/
//   FUNCTION Prototype: xf_Jacobian_SolveM_BlockJacobiLean
extern int
xf_Jacobian_SolveM_BlockJacobiLean(xf_All *All, xf_JacobianMatrix *R_U,
				   xf_Vector *X, enum xfe_AddType AddFlag,
				   enum xfe_Bool TransposeFlag,
				   xf_SolverData *SolverData);
/*
PURPOSE:

  X @= M^{-1} * X      
  X @= M^{-T} * X      (TransposeFlag == True)

  "@=" corresponds to AddFlag

  M is the block preconditioner part of R_U.  It is calculated on the
  fly in a memory-lean fashion.

INPUTS:

  All : All structure
  R_U : Jacobian matrix -- no allocated values, but containing rank and
        neighbor connectivity information
  X   : vector to be multiplied
  AddFlag : operation indicating set, neg, add, sub
  TransposeFlag : True if the transpose opeartion is requested
  SolverData : solver data structure

OUTPUTS: 

  X   : multiplied vector

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Jacobian_MultM_BlockJacobiLean
extern int
xf_Jacobian_MultM_BlockJacobiLean(xf_All *All, xf_JacobianMatrix *R_U,
				  xf_Vector *X, enum xfe_AddType AddFlag,
				  enum xfe_Bool TransposeFlag, 
				  xf_SolverData *SolverData, xf_Vector *Y);
/*
PURPOSE:

  Y @= M   * X      
  Y @= M^T * X      (TransposeFlag == True)

  "@=" corresponds to AddFlag

  M is the block preconditioner part of R_U.  It is calculated on the
  fly in a memory-lean fashion.

INPUTS:

  All : All structure
  R_U : Jacobian matrix -- no allocated values, but containing rank and
        neighbor connectivity information
  X   : vector to be multiplied
  AddFlag : operation indicating set, neg, add, sub
  TransposeFlag : True if the transpose opeartion is requested
  SolverData : solver data structure

OUTPUTS: 

  Y   : multiplied vector

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_Jacobian_MultN_BlockJacobiLean
extern int
xf_Jacobian_MultN_BlockJacobiLean(xf_All *All, xf_JacobianMatrix *R_U,
				  xf_Vector *X, enum xfe_AddType AddFlag,
				  enum xfe_Bool TransposeFlag, 
				  xf_SolverData *SolverData,
				  enum xfe_Bool HaloFlag,  xf_Vector *Y);
/*
PURPOSE:

  Y @= N   * X      
  Y @= N^T * X      (TransposeFlag == True)

  "@=" corresponds to AddFlag

  N is the off-diagonal part of R_U.  It is calculated on the
  fly in a memory-lean fashion.

INPUTS:

  All : All structure
  R_U : Jacobian matrix -- no allocated values, but containing rank and
        neighbor connectivity information
  X   : vector to be multiplied
  AddFlag : operation indicating set, neg, add, sub
  TransposeFlag : True if the transpose opeartion is requested
  SolverData : solver data structure
  HaloFlag : if True, only Halo effects are multiplied
             if False, effect of halo elements is skipped.

OUTPUTS: 

  Y   : multiplied vector

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_LeanSolverStepBlockElem
extern int
xf_LeanSolverStepBlockElem( xf_All *All, int egrp, int elem,
			    xf_SolverData *SolverData, xf_Vector *U,
			    xf_JacobianMatrix *R_U, xf_Vector *dU);
/*
PURPOSE:

  Applies one step of a memory-lean block solver to obtain a state
  update, dU, for one element only

INPUTS:

  All : All structure
  egrp, elem : element on which lean block step is applied
  SolverData : structure of useful solver data + params
  U : state vector
  R_U : Jacobian matrix without values, but with connectivity

OUTPUTS: 

  dU : state update

RETURN:

  Error Code
*/


/* /\******************************************************************\/ */
/* //   FUNCTION Prototype: xf_LeanSolverStep */
/* extern int  */
/* xf_LeanSolverStep( xf_All *All, enum xfe_NonlinearSolverType NonlinearSolver,  */
/* 		   xf_SolverData *SolverData, xf_Vector *U,  */
/* 		   xf_JacobianMatrix *R_U, xf_Vector *dU); */
/* /\* */
/* PURPOSE: */

/*   Applies a step of the nonlinear solver stored in NonlinearSolver. */
/*   The state update is returned in dU.  Does not perform a halo */
/*   transfer before calculations. */

/* INPUTS: */

/*   All : All structure */
/*   NonlinearSolver : type of nonlinear solver to use (has to be lean) */
/*   SolverData : structure of useful solver data + params */
/*   U : state vector */
/*   R_U : Jacobian matrix -- no allocated values, but containing rank and */
/*         neighbor connectivity information */


/* OUTPUTS:  */

/*   dU : state update */

/* RETURN: */

/*   Error Code */
/* *\/ */




#endif // end ifndef _xf_LeanSolver_h

