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

#ifndef _xf_Penalty_h
#define _xf_Penalty_h 1

/*
 FILE:  xf_Penalty.h
 
 This file contains the headers for interior penalty functions.
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_FindPenaltyHessian
extern int
xf_FindPenaltyHessian(xf_All *All, int egrp, xf_Vector *U, xf_Matrix **pHp);
/*
 PURPOSE: 
 
 Finds (allocates) the memory the Hessian of the penalty 
 function for a given element group.
 
 INPUTS:
 
 All : All structure
 egrp: element group index
 U   : Primal state vector
 pHP : Pointer for the Hessian Matrix
 
 OUTPUTS:
 
 None : pHp gets modified.
 
 RETURN:
 
 Error code
 
 */



/******************************************************************/
//   FUNCTION Definition: xf_CalculatePenalty
extern int
xf_CalculatePenalty(xf_All *All, xf_Vector *U, xf_Vector *P, real *pP);
/*
 PURPOSE: 
 
 Calculates the penalty value for a given state vector U.
 The penalty function shape depends on the EqnSet used.
 
 INPUTS:
 
 All : All structure
 U   : Primal state vector
 P   : Vector with one penalty integral per element per unit volume
 *pP : Pointer for the sum of penalty function values
 
 OUTPUTS:
 
 P : Receives a penalty integral per unit volume for each element
 *pP: Penalty integral per unit volume for the Domain
 
 
 RETURN:
 
 Error code
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_CalculatePenaltyGradient
extern int
xf_CalculatePenaltyGradient(xf_All *All, xf_Vector *U, xf_Vector *Gp);
/*
 PURPOSE: 
 
 Calculates the gradient of the penalty function for a given state vector U.
 The penalty function shape depends on the EqnSet used.
 
 INPUTS:
 
 All : All structure
 U   : Primal state vector
 Gp  : Gradient of the penalty function
 
 OUTPUTS:
 
 None : GP gets changed
 
 RETURN:
 
 Error code
 
 */



/****************************************************************************/
//   FUNCTION Definition: xf_CalculatePenaltyHessian
extern int
xf_CalculatePenaltyHessian(xf_All *All, xf_Vector *U, xf_Matrix *Hp, int egrp, int elem);
/*
 PURPOSE: 
 
 Calculates the Hessian of the penalty function for a given state vector U.
 The penalty function shape depends on the EqnSet used.
 
 INPUTS:
 
 All : All structure
 U   : Primal state vector
 Hp  : Hessian matrix of the penalty function
 egrp, elem : coordinates to a element. 
 
 OUTPUTS:
 
 None : Hp gets changed
 
 RETURN:
 
 Error code
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_ApplyPenaltyHessian
extern int
xf_ApplyPenaltyHessian(xf_All *All, xf_Vector *U, xf_Vector *Y, 
                       xf_Vector *W, enum xfe_Bool CalcHessian, 
                       enum xfe_Bool InverseFlag, enum xfe_AddType AddFlag);

/*
 PURPOSE: 
 
 Computes W @= Hp*Y if InverseFlag = false. 
          W @= Hp^-1*Y otherwise
 
 INPUTS:
 
 All : All structure
 U   : Primal state vector
 Y   : vector to apply on
 W   : vector to put the result
 CalcHessian: If true, Hessian gets recalculated
 InverseFlag: Described in purpose
 AddFlag: Add, subtract or set.
 
 OUTPUTS:
 
 None : W gets changed
 
 RETURN:
 
 Error code
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_ApplyResidualHessian
extern int
xf_ApplyResidualHessian(xf_All *All, xf_Vector *U, xf_Vector *Y, 
                        xf_Vector *W, xf_SolverData *SolverData, 
                        enum xfe_AddType AddFlag);
/*
 PURPOSE: 
 
 Computes W @= 2*(dR_dU)^T*dR_dU*dU + d2Pe_dUidUj*dU
 This function assumes dR_dU is up to date with U
 
 INPUTS:
 
 All : All structure
 U   : Primal state vector
 Y   : vector to apply on
 W   : vector to put the result
 SolverData: Struture that stores solver parameters
 AddFlag: Add, subtract or set.
 
 OUTPUTS:
 
 None : W gets changed
 
 RETURN:
 
 Error code
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_ApplyHessian
extern int
xf_ApplyHessian(xf_All *All, xf_Vector *U, xf_Vector *dU, xf_Vector *W,
                xf_SolverData *SolverData, enum xfe_Bool CalcHessian);
/*
 PURPOSE: 
 
 Wrapper for ApplyResidualHessian and ApplyPenaltyHessian 
 It applies the residual Hessian and the Penalty Hessian
 respectively
 
 INPUTS:
 
 All : All structure
 U   : Primal state vector
 Y   : vector to apply on
 W   : vector to put the result
 SolverData: Struture that stores solver parameters
 AddFlag: Add, subtract or set.
 
 OUTPUTS:
 
 None : W gets changed
 
 RETURN:
 
 Error code
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_ApplyHessianPrecond
extern int
xf_ApplyHessianPrecond(xf_All *All, xf_Vector *U, xf_Vector *Y, xf_Vector *W,
                       xf_SolverData *SolverData, enum xfe_Bool CalcHessian);
/*
 PURPOSE: 
 
 This function applies the approximate inverse when the penalty 
 Hessian dominates. This approximate is computed using 
 Taylor expansion of Hp^-1 * (I + H*Hp)^-1
 W = (I - Hp^-1*Hr + (Hp^-1*Hr)*(Hp^-1*Hr))*Hp^-1*Y
 
 INPUTS:
 
 All : All structure
 U   : Primal state vector
 Y   : vector to apply on
 W   : vector to put the result
 SolverData: Struture that stores solver parameters
 CalcHessian: If true, Hessian gets recalculated
 
 OUTPUTS:
 
 None : W gets changed
 
 RETURN:
 
 Error code
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_PenalizeResidual
extern int
xf_PenalizeResidual(xf_All *All, xf_Vector *U, xf_Vector *R, 
                    xf_SolverData *SolverData);

/*
 PURPOSE: 
 
 This function applies the penalty matrix to R:
 R *= (I+Phi).
 Check AIAA 2011-3696 Paper.
 
 INPUTS:
 
 All : All structure
 U   : Primal state vector
 R   : Residual vector
 SolverData: Struture that stores solver parameters
 
 OUTPUTS:
 
 None : R gets updated
 
 RETURN:
 
 Error code
 
 */

#endif // end ifndef _xf_Penalty_h

