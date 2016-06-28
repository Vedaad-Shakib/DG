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

#ifndef _xf_SolverROBST_h
#define _xf_SolverROBST_h 1

/*
 FILE:  xf_SolverROBST.h
 
 This file contains the headers for the functions related 
 to the ROBST solver
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_Calculate_F
extern int
xf_Calculate_F(xf_All *All, xf_Vector *U, xf_SolverData *SolverData, 
               enum xfe_Bool UseResidual);

/*
 PURPOSE:
 
 Calculates the value of f(U) = R(U)^2 + P(U),
 where R and P are the residual and penalty functions recpectively
 
 INPUTS:
 
 All: All structure
 U: primal state 
 SolverData: structure that will receive the value of 
             SolverData->AugResidual = f(U) and
             SolverData->PenaltyFcn = P(U)
 UseResidual: If true, the residual R(U) is not recalculated
 
 OUTPUTS: 
 
 None: SolverData get modified
 
 RETURN:
 
 Error Code
 */


/******************************************************************/
//   FUNCTION Definition: xf_Calculate_dF
extern int
xf_Calculate_dF(xf_All *All, xf_Vector *U, xf_Vector *P, real *dF, 
                xf_SolverData *SolverData, enum xfe_Bool UseGradient,
                enum xfe_Bool Normalize);

/*
 PURPOSE:
 
 Calculates the derivative of f(U) in the direction of P in U space
 
 INPUTS:
 
 All: All structure
 U: primal state 
 P: Trial direction
 dF: this will receive the value of the derivative
 SolverData: structure that will be used only if UseGradient is false
 UseGradient: If true, the gradient is not recalculated
 Normalize: If true, the derivative is divided by the norm of G
 
 OUTPUTS: 
 
 (*dF): the dot product G.*P is stored here.
 
 RETURN:
 
 Error Code
 */

/******************************************************************/
//   FUNCTION Definition: xf_SolveNonlinearSystem_ROBST
extern int
xf_SolveNonlinearSystem_ROBST(xf_All *All, enum xfe_Bool LinearFlag, xf_Vector *S, 
                              xf_Vector *U);

/*
 PURPOSE:
 
 Solves R(U) + S = 0 using a Residual Optimization Back-Solving
 Technique
 
 INPUTS:
 
 All: All structure
 LinearFlag: If true, the Hessian is calculated only once
 S: Source vector
 U: State vector
 
 OUTPUTS: 
 
 None: The state vector gets modified
 
 RETURN:
 
 Error Code
 */

#endif // end ifndef _xf_SolverROBST_h
