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

#ifndef _xf_ResidualStab_h
#define _xf_ResidualStab_h 1

/*
  FILE:  xf_ResidualStab.h

  This file contains the headers for stabilization residual-calculation functions

*/

#include "xf_SolverStruct.h"



/******************************************************************/
//   FUNCTION Prototype: xf_InitStabData
extern void
xf_InitStabData( xf_StabData *StabData);
/*
PURPOSE:

  Initializes stabilization data structure
 
INPUTS:

  StabData : pointer to structure

OUTPUTS:

  StabData is initialized
 
RETURNS:

  None

*/


/******************************************************************/
//   FUNCTION Prototype: xf_CalculateStabilization
extern int
xf_CalculateStabilization(xf_All *All, xf_Vector *U, 
			  enum xfe_Bool Need_Grad,
			  xf_SolverData *SolverData);
/*
PURPOSE:

  Determines if any stabilization terms are required, and, if so,
  calculates necessary data for them.
 
INPUTS:

 
  All     : All structure
  U       : State Vector
  Need_Grad : flag indicating that linearization is needed
  SolverData : contains stab structure

OUTPUTS:

  None, data in SolverData is modified
 
RETURNS:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_CalculateStabViscElem
extern int
xf_CalculateStabViscElem(xf_All *All, int egrp, int elem, int nq, real *xq, 
			 real *EM, const xf_StabData *StabData,
			 real *StabVisc, real *ResMetric, enum xfe_Bool *SkipDiffStab);
/*
PURPOSE:

  Calculates stabilization viscosity on an element, at nq points.
 
INPUTS:

  All : All structure
  egrp, elem : element in question
  nq : number of points
  xq : location of points (elem ref space)
  EM : element metric coefficients (Q1 interpolation)
  StabData : data structure containing stabilization information

OUTPUTS:

  StabVisc  : stabilization viscosities
  ResMetric    : resolution matrix (dim x dim) at each point
  SkipDiffStab : true if no stabilization is necessary for this elem
 
RETURNS:

  Error code

*/
  


/******************************************************************/
//   FUNCTION Prototype: xf_CalculateStabViscElem
extern int
xf_AddLinearizationStabViscElem(xf_All *All, int egrp, int elem, int sr, int nn, int nq,
				real *xq, real *StabVisc, const xf_StabData *StabData,
				real *gPhi, real *wq, real *AwStab, real *T,
				real *ER_U, real **ER_NU);
/*
PURPOSE:

  Adds linearization of stabilization viscosity to element residual:

  ER_U{n,k;m,a} += gPhi{i,n,q}*AwStab{i,q,k}*wq{q}*StabVisc_U{q,m,a}
 
INPUTS:

  All : All structure
  egrp, elem : element in question
  sr : state rank
  nn : order rank of state on elem
  nq : number of points
  xq : elem ref locations of points
  StabVisc : values of stabilization viscosities
  StabData : data structure containing stabilization information
  gPhi  : element basis gradients
  wq    : quadrature weights
  AwStab : A matrix multiplying vector of interest, not including stab viscosities
  T  : temporary work matrix

OUTPUTS:

  ER_U : residual linearization on elem w.r.t U on elem (added to)
  ER_NU : residual linearization on elem w.r.t U on neighbors (added to)

RETURNS:

  Error code

*/
  


/******************************************************************/
//   FUNCTION Prototype: xf_CalculateStabViscIFace
extern int
xf_CalculateStabViscIFace(xf_All *All, int iiface, int nq, real *xelemL,
			  real *xelemR, real *EML, real *EMR,
			  const xf_StabData *StabData,
			  real *StabViscL, real *StabViscR, 
			  real **pStabViscL_UL, real **pStabViscL_UR,
			  real **pStabViscR_UL, real **pStabViscR_UR,
			  real *StabPhiL, real *StabPhiR,
			  real *ResMetricL, real *ResMetricR,
			  enum xfe_Bool *SkipDiffStab);
/*
PURPOSE:

  Calculates stabilization viscosity on an interior face
 
INPUTS:

  All : All structure
  iiface : interior face number in question
  nq : number of points
  xelemL : L elem ref space locations of points
  xelemR : R elem ref space locations of points
  EML : element metric coefficients on L
  EMR : element metric coefficients on R
  StabData : data structure containing stabilization information
  
OUTPUTS:

  StabViscL : stabilization viscosities on left elem
  StabViscR : stabilization viscosities on right elem
  (*pStabViscL_UL) : linearization of L stab viscosity w.r.t UL
  (*pStabViscL_UR) : linearization of L stab viscosity w.r.t UR
  (*pStabViscR_UL) : linearization of R stab viscosity w.r.t UL
  (*pStabViscR_UR) : linearization of R stab viscosity w.r.t UR
  StabPhiL      : face viscosity basis interpolation values on L
  StabPhiR      : face viscosity basis interpolation values on R
  ResMetricL    : resolution matrix (dim x dim) on L elem at each point
  ResMetricR    : resolution matrix (dim x dim) on R elem at each point
  SkipDiffStab   : True if no stabilization is necessary on this face

RETURNS:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_CalculateStabViscBFace
extern int
xf_CalculateStabViscBFace(xf_All *All, int ibfgrp, int ibface, int nq, real *xelem,
			  real *EM, const xf_StabData *StabData,
			  real *StabVisc, real **pStabVisc_U, 
			  real *StabPhi, real *ResMetric, enum xfe_Bool *SkipDiffStab);
/*
PURPOSE:

  Calculates stabilization viscosity for points on a boundary face
 
INPUTS:

  All : All structure
  ibfgrp, ibface : boundary face number in question
  nq : number of points
  xelem : elem ref space locations of points
  EM : element metric coefficients
  StabData : data structure containing stabilization information
  
OUTPUTS:

  StabVisc       : stabilization viscosities at nq points
  (*pStabVisc_U) : linearization of stab viscosity w.r.t U
  StabPhi        : face viscosity basis interpolation values 
  ResMetric      : resolution matrix (dim x dim) at each point
  SkipDiffStab   : True if no stabilization is necessary on this face

RETURNS:

  Error code

*/




/******************************************************************/
//   FUNCTION Prototype: xf_CalculateStabViscBFace
extern int
xf_StabDiffA(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
	     const int *IParam, const real *RParam, const real *ResMetric, 
	     int nq, const real *u, const real *w, real *Aw, real *A, 
	     real *A_uw, enum xfe_Bool *ConstA, real *mu);
/*
PURPOSE:

  Evaluates the diffusion "A matrix" corresponding for stabilization.
  Note, the stabilization viscosity is not added in this A matrix, as
  it is added separately so that the effect of the stabilization
  viscosity linearization can be separated.  A non-Null w means that
  the product Aw = A*w, is desired.  A_uw is not changed, since the
  stabilization A matrix is constant w.r.t u.  All outputs are
  incremented, not set.

INPUTS:

  EqnSet  : equation set structure
  ResTerm : the diffusion residual term to evaluate
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  ResMetric : resolution metric [length]; used for artificial stabilization
  nq      : number of times to calculate A
  u       : nq state vectors, unrolled along StateRank first
  w       : supplementary vector for matrix-vector products. 
            Also at nq points, unrolled along StateRank first.

OUTPUTS: 

  Aw   : dim*nq vectors, unrolled along sr(= StateRank) first, then 
         along nq, then along dim,

 	 Aw{i,q,k} = Aw[i*nq*sr + q*sr + k]

  A    : dim*dim*nq matrices, each sr x sr:

         A{i,j,q,k,l} = A[(i*dim+j)*nq*sr*sr + q*sr*sr + k*sr + l]

	 This is the no-frills attached matrix A
  A_uw : not set

  ConstA : if True, means A is not a const w.r.t u.  In this case, A_uw
           will have been computed, assuming it was not passed in as NULL.

  mu : viscosity (or some measure thereof) at nq points.  Used for
       line creation.  This variable is added-to.

RETURN:

  Error Code
*/



#endif // end ifndef _xf_ResidualStab_h
