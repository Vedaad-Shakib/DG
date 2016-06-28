#ifndef _xf_ResidualDiff_h
#define _xf_ResidualDiff_h 1

/*
  FILE:  xf_ResidualDiff.h

  This file contains the headers for the diffusion-specific
  residual-calculation functions

*/

#include "xf_LinQueueStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_CreateDiffJumpData
extern int
xf_CreateDiffJumpData(xf_DiffJumpData **pDData);
/*
PURPOSE:
 
  Creates a data structure for storing diffusion-specific matrices for
  terms associated with inter-element jumps.  Initializes constants
  and pointers.

INPUTS:

  None

OUTPUTS:

  pDData : pointer to created diffusion data structure

RETURNS:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReAllocDiffJumpData
extern int
xf_ReAllocDiffJumpData(int sr, int nq, int dim, int nn, 
		       enum xfe_Bool Need_Grad, xf_DiffJumpData *DData);
/*
PURPOSE:
 
  Reallocates matrices in DData to accomodate sr, nq, dim, nn.
  Specific for inter-element jump diffusion terms.  DData must already
  have been creaated.

INPUTS:

  sr : state rank
  nq : number of quadrature points
  dim : dimension
  nn : number of basis functions
  Need_Grad : if true, matrices required for derivative computations 
              will also be allocated.

OUTPUTS:

  DData : pointer containing allocated diffusion data structure

RETURNS:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyDiffJumpData
extern void
xf_DestroyDiffJumpData(xf_DiffJumpData *DData);
/*
PURPOSE:
 
  Destroys DData and all its associated matrices.

INPUTS:

  DData: a diffusion-term data structure

OUTPUTS:

  None

RETURNS:

  Error code

*/



/******************************************************************/
//   FUNCTION Definition: xf_ComputeDiffA
extern int 
xf_ComputeDiffA(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm,
		int nDiff, const int *IParam, const real *RParam,
		const real *ResMetric, const real *StabVisc, int nq, 
		real *u, real *gu, real *w, enum xfe_Bool uIsRef, 
		const xf_MotionData *MD, real *Aw, real *A, 
		real *A_uw, real *AwStab, enum xfe_Bool *ConstA, real *CF);
/*
PURPOSE:
 
  Wrapper for equation-specific DiffA function.  Evaluates the total
  diffusion "A matrix" corresponding to the nDiff ResTerms, at state
  u, at nq points.  The A matrix multiplies grad U in the
  discretization to yield the visous flux.  A non-Null W means that
  the products AW = A*W, and A_UW = A_U*W (if A_uw != NULL) are
  desired.  The A-matrix outputs are set to zero at the beginning of
  the function.  The connectivity, CF, is added-to.

  All stabilization diffusion matrices (and outputs) are multiplied by
  the stabilization viscosity, StabVisc, before being added to the
  total sum.  Thus, if StabVisc == 0.0, the stabilization As are not
  even computed. A distinctive feature of this function is the output
  AwStab, which is the Aw product formed from only the stabilization A
  matrices multiplying w, not including the StabVisc multiplier.
  AwStab is required for computing the linearization (it multiplies
  StabVisc_U via the product rule).

  In the presence of ALE-based mesh motion (MD != NULL) the matrices
  are modified according to:
  
          Aref     = G^{-1} * Aphys * G^{-T}
          Awref    = G^{-1} * Aphys * G^{-T}*w
	  A_uwref  = G^{-1} * 1/g * A_uwphys * G^{-T}*w

  where "ref" refers to the reference space, and phys refers to the
  physical space.  G is the ref-to-glob coord mapping Jacobian, and g
  = det(G).  The factor of 1/g on A_uw accounts for the chain rule so
  that A_uwref is the linearization w.r.t the ref state uref =
  g*uphys.

INPUTS:

  EqnSet    : equation set structure
  ResTerm   : the diffusion residual term(s) to evaluate
  nDiff     : number of diffusion terms to evaluate
  IParam    : vector of integer params as requested by EqnSet->IParamKey
  RParam    : vector of real params as requested by EqnSet->RParamKey
  ResMetric : dimensionless resolution metric [dim^2]; used for artificial stabilization
  StabVisc  : stabilization viscosity optional; only if stabilization is being used
  nq        : number of times to calculate A
  u         : nq state vectors, unrolled along StateRank first
  gu        : nq state gradient vectors, unrolled along StateRank first
  w         : supplementary vector for matrix-vector products. 
              Also at nq points, unrolled along StateRank first.
  uIsRef    : flag indicating whether u is in physical (F) or reference (T) space.
  MD        : mesh motion data structure (optional, indicates motion is on)

OUTPUTS:

  Aw   : dim*nq vectors, unrolled along sr(= StateRank) first, then 
         along nq, then along dim,

 	 Aw{i,q,k} = Aw[i*nq*sr + q*sr + k]

	 If w is grad u, Aw is the viscous flux vector.
  A    : dim*dim*nq matrices, each sr x sr:

         A{i,j,q,k,l} = A[(i*dim+j)*nq*sr*sr + q*sr*sr + k*sr + l]

	 This is the no-frills attached matrix A
  A_uw : dim*nq matrices, each sr x sr:
  
         A_uw{i,q,k,a} = A_uw[i*nq*sr*sr + q*sr*sr + k*sr + l]
	               = A_u{i,j,q,k,l;a}* * w{j,q,l}  (sum over j, l)

	 Only computed if A is not const w.r.t U.

  AwStab : dim*nq vectors, unrolled along sr(= StateRank) first, then 
          along nq, then along dim,

 	  AwStab{i,q,k} = AwStab[i*nq*sr + q*sr + k]

	  This is the product of the stabilization A matrices with w,
	  not including the stabilization viscosity, StabVisc.  Only
	  computed if StabVisc != NULL.

  ConstA : if True, means A is not a const w.r.t U.  In this case, A_UW
           will have been computed, assuming it was not passed in as NULL.

  CF : viscosity (or some measure of coupling) at nq points.  Used for
       line creation.  This variable is added-to.

RETURNS:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_DiffFluxQJump
extern int
xf_DiffFluxQJump(xf_All *All, const int iiface, const xf_ResTerm *ResTerm, 
		 int nDiff, enum xfe_DiffDiscType DiffDisc, const int *IParam, 
		 const real *RParam, const xf_BasisData *PhiDataL, 
		 const xf_BasisData *PhiDataR, const xf_BasisData *ResPhiDataL,
                 const xf_BasisData *ResPhiDataR, const xf_StabData *StabData,
		 real ElemVolL, real ElemVolR, enum xfe_BasisType BasisL, 
		 enum xfe_BasisType BasisR, int OrderL, int OrderR, int nq, 
		 const real *wn, real *uL, real *uR, real *guL, real *guR, 
		 const xf_MotionData *MDL, const xf_MotionData *MDR, real *RL, 
		 real *RR, real *RL_UL, real *RL_UR, real *RR_UL, real *RR_UR, 
		 xf_LinQueueData *LinQ, xf_DiffJumpData *DData, real *CF);
/*
PURPOSE:
 
  Computes diffusion flux, Q, and any associated derivatives, for an
  inter-element face.  Derivatives are computed if the derivative
  structures in DData are not Null.  Residual and derivatives are
  updated with the integrated flux, int(Qn*Phi), using the input basis
  functions in PhiDataL/R.  Stabilization is included in Qn,
  corresponding to the diffusion discretization type in DiffDisc.

INPUTS:

  All : all structure
  IFace : interior face
  ResTerm : residual terms to consider
  nDiff : number of res terms (in sequence in ResTerm)
  DiffDisc : type of diffusion discretization to use
  IParam, RParam : integer/real parameters passed into eqnset functions
  PhiDataL, PhiDataR : basis functions + gradients from L and R elements
  ResPhiDataL, ResPhiDataR : possibly-different basis functions for residual evaluation
  StabData : structure containing data for stabilization (or NULL)
  ElemVolL, ElemVolR : volumes of L and R elements
  BasisL, BasisR : bases of elems L and R
  OrderL, OrderR : orders of elems L and R
  nq : number of quad points on face
  wn : normal vectors at nq points, with quad weights included
  uL/uR : left/right states at nq points
  guL, guR : left/right gradients at nq points
  MDL/MDR : left/right motion data structures (non-null indicates motion is on)

OUTPUTS:

  RL, RR : residuals that are incremented with the integrated Qn
  RL_UL, RR_UL : derivatives of residuals w.r.t. UL
  RL_UR, RR_UR : derivatives of residuals w.r.t. UR
  LinQ : updated linearization queue with some of the contributions to R_U
  DData: diffusion-term data structure containing Q and any derivatives
  CF : inter-element coupling at the nq points


RETURNS:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_CreateDiffBCData
extern int
xf_CreateDiffBCData(xf_DiffBCData **pDData);
/*
PURPOSE:
 
  Creates a data structure for storing diffusion-specific matrices for
  terms associated with boundary jumps.  Initializes constants
  and pointers.

INPUTS:

  None

OUTPUTS:

  pDData : pointer to created diffusion data structure

RETURNS:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReAllocDiffBCData
extern int
xf_ReAllocDiffBCData(int sr, int nq, int dim, int nn,
		     enum xfe_Bool Need_Grad, xf_DiffBCData *DData);
/*
PURPOSE:
 
  Reallocates matrices in DData to accomodate sr, nq, dim, nn.
  Specific for boundary jump diffusion terms.  DData must already have
  been creaated.

INPUTS:

  sr : state rank
  nq : number of quadrature points
  dim : dimension
  nn : number of basis functions
  Need_Grad : if true, matrices required for derivative computations 
              will also be allocated.

OUTPUTS:

  DData : pointer containing allocated diffusion data structure

RETURNS:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyDiffBCData
extern void
xf_DestroyDiffBCData(xf_DiffBCData *DData);
/*
PURPOSE:
 
  Destroys DData and all its associated matrices.

INPUTS:

  DData: a diffusion-term data structure

OUTPUTS:

  None

RETURNS:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_DiffFluxQBC
extern int
xf_DiffFluxQBC(xf_All *All, int ibfgrp, int ibface, const xf_ResTerm *ResTerm, 
	       int nDiff, enum xfe_DiffDiscType DiffDisc, const int *IParam, 
	       const real *RParam, const xf_BasisData *PhiData, 
               const xf_BasisData *ResPhiData, const xf_StabData *StabData, 
	       real ElemVol, enum xfe_BasisType Basis, int Order, int nq, const real *wn,  
	       const real *xglob, real *uI, real *uB, const real *uB_uI, real *guI, 
	       xf_MotionData *MD, xf_OutputEvalData *OutputEval, real *ER, real *ER_U, 
	       xf_LinQueueData *LinQ, xf_DiffBCData *DData, real *CF);
/*
PURPOSE:
 
  Computes diffusion flux, Q, and any associated derivatives, for a
  boundary face.  Derivatives are computed if the derivative
  structures in DData are not Null.  The integrated flux dotted with
  the normal, Qn, is added to the residual, ER.  Stabilization is
  included in Qn.  If Value (EV_U for linearization) is not NULL, an
  output derived from a boundary integral of a linear combination of
  the flux components is computed.

INPUTS:

  All : all structure
  ibfgrp, ibface : boundary face group number and number
  ResTerm : residual terms to consider
  nDiff : number of res terms (in sequence in ResTerm)
  DiffDisc : type of diffusion discretization to use
  IParam, RParam : integer/real parameters passed into eqnset functions
  PhiData : basis functions + gradients on elem
  ResPhiData : possibly-different basis functions for residual evaluation
  StabData : structure containing data for stabilization (or NULL)
  ElemVol : volume of element
  Basis : solution basis
  Order : solution order
  Res : resolution vector (for artificial viscosity) [length]
  nq : number of quad points on face
  wn : normal vectors at nq points, with quad weights included
  xglob : global coords of nq points
  uI/uB : interior/boundary states at nq points
  uB_uI : derivative of boundary states w.r.t interior states
  guI : interior gradients at nq points
  MD : motion data structure (non-null indicates motion is on) 
  OutputEval : output evaluation structure

OUTPUTS:

  ER : residual incremented with the integrated Qn
  ER_U : derivative of ER w.r.t. UI
  LinQ : updated linearization queue with some of the contributions to R_U
  DData: diffusion-term data structure containing Q and any derivatives
  CF : element-boundary coupling at the nq points

RETURNS:

  Error code

*/


#endif // end ifndef _xf_ResidualDiff_h
