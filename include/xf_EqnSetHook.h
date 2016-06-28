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

#ifndef _xf_EqnSetHook_h
#define _xf_EqnSetHook_h 1

/*
  FILE:  xf_EqnSetHook.h

  This file contains the headers for Equation-set specific functions.

*/


/* The first two prototypes are for opening/closing the library */

/******************************************************************/
//   FUNCTION Prototype: xf_LoadEqnSetLibrary
extern int
xf_LoadEqnSetLibrary(const char *LibName);


/******************************************************************/
//   FUNCTION Prototype: xf_CloseEqnSetLibrary
extern void
xf_CloseEqnSetLibrary();
/*
PURPOSE:

  Closes the equation set library

INPUTS:

  None;

OUTPUTS: 

  None;

RETURN:

  Error Code
*/


/* The remaining prototypes are for the equation-set specific hooks. */


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetRegister
extern int 
xf_EqnSetRegister(xf_EqnSet *EqnSet);
/*
PURPOSE:

  Interprets the residual terms that make up the equation.  These
  residual terms must be present in EqnSet->ResTerm.  Based on the
  terms present, this function decides which state vector components
  are going to participate in the solve, and appropriately fills in
  StateRank, StateName, and PositionInState (see EqnSetStruct.h for
  definitions).


INPUTS:

  EqnSet: the equation-set structure

OUTPUTS: 

  None

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetQuadOrderElem
extern int 
xf_EqnSetQuadOrderElem(const xf_EqnSet *EqnSet, int Order, int *QuadOrder);
/*
PURPOSE:

  Returns the required quadrature order for performing element
  integrals of the ResTerms stored in EqnSet, when the state is
  represented with order Order.

INPUTS:

  EqnSet : equation set structure
  Order : state interpolation order

OUTPUTS: 

  (*QuadOrder): contains the quadrature order

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetQuadOrderFace
extern int 
xf_EqnSetQuadOrderFace(const xf_EqnSet *EqnSet, int Order, int *QuadOrder);
/*
PURPOSE:

  Returns the required quadrature order for performing face
  integrals of the ResTerms stored in EqnSet, when the state is
  represented with order Order.

INPUTS:

  EqnSet : equation set structure
  Order : state interpolation order

OUTPUTS: 

  (*QuadOrder): contains the quadrature order

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetICState
extern int 
xf_EqnSetICState(const xf_EqnSet *EqnSet, const xf_IC *IC, 
		 const int *IParam, const real *RParam, int nq,
		 const real *xglob, const real *pTime, real *U);

/*
PURPOSE:

  Fills in state vector, u, at nq points, using initial condition IC
  and possibly the global positions in xglob

INPUTS:

  EqnSet : the equation-set structure
  IC     : initial-condition structure
  IParam : vector of integer params as requested by EqnSet->IParamKey
  RParam : vector of real params as requested by EqnSet->RParamKey
  nq     : number of points at which to evaluate ICs
  xglob  : locations of points in physical space
  pTime  : pointer to simulation time (optional)

OUTPUTS: 

  U : vector of data that gets initialized.  Must be pre-allocated [nq*sr]

RETURN:

  Error Code (e.g. if the initial condition type or data are not understood)
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetBCState
extern int 
xf_EqnSetBCState(const xf_EqnSet *EqnSet, const xf_BC *BC,
		 const int *IParam, const real *RParam, int nq, 
		 const real *n, const real *xglob, const real *pTime, 
		 const real *vg, const real *UI, real *UB, real *UB_UI);
/*
PURPOSE:

  Evaluates nq boundary states, UB, corresponding to interior states,
  UI, using boundary condition information in BC.  When boundary
  velocities are given, they are used in setting the BC state.

INPUTS:

  EqnSet  : equation set structure
  BC      : the boundary condition for which to evaluate UB
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of times to calculate UB
  n       : nq normal vectors
  xglob   : global positions of the nq points
  pTime   : pointer to simulation time (optional)
  vg      : boundary velocity at nq points [nq*dim]
  UI      : nq state vectors, unrolled along StateRank first (Interior)

OUTPUTS: 

  UB      : nq state vectors, unrolled along StateRank first (Boundary)
  UB_UI   : nq srxsr derivative matrices (boundary state w.r.t interior)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetFcnState
extern int 
xf_EqnSetFcnState(const xf_EqnSet *EqnSet, const char *Fcn, 
		  const char *Data, const int *IParam, 
		  const real *RParam, int nq, const real *xglob, 
		  const real *pTime, real *u);
/*
PURPOSE:

  Fills in state vector, u, at nq points, using Fcn with Data,
  and possibly the global positions in xglob

INPUTS:

  EqnSet : the equation-set structure
  Fcn    : name of function to use
  Data   : data string passed into function
  IParam : vector of integer params as requested by EqnSet->IParamKey
  RParam : vector of real params as requested by EqnSet->RParamKey
  nq     : number of points at which to evaluate ICs
  xglob  : locations of points in physical space
  pTime   : pointer to simulation time (optional)

OUTPUTS: 

  U : vector of data that gets set  Must be pre-allocated [nq*sr]

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetConvF
extern int 
xf_EqnSetConvF(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm,
	       int nConv, const int *IParam, const real *RParam, int nq, 
	       const real *U, real **AuxU, const real *xglob, 
	       real *F, real *F_U);
/*
PURPOSE:

  Evaluates the element-interior flux for the convection residual
  term(s) stored in ResTerm, at nq points.

INPUTS:

  EqnSet  : equation set structure
  ResTerm : the convection residual term(s) to evaluate
  nConv   : number of residual terms
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of times to calculate the interior flux
  U       : nq state vectors, unrolled along StateRank first
  AuxU    : pointers to any requested auxiliary vectors.  Each vector
            is interpolated at the nq points.
  xglob   : global positions of the nq points

OUTPUTS: 

  F   : EqnSet->Dim*nq flux vectors, unrolled along StateRank first,
        then along nq, then along dim. Specifically,

	F[d*StateRank*nq + StateRank*iq + k]

	stores the kth component of the dim=d flux at point iq.
	Storage for F is pre-allocated before the call.  

  F_U : only computed if F_U != NULL.  Stores flux Jacobians for
        EqnSet->Dim dimensions and nq points, unrolled along
        StateRank^2 first, then nq, then dim.  Specifically,
	  
	F_U [d*StateRank^2*nq + StateRank^2*iq + StateRank*k + a] 

	stores the derivative of dim=d flux term k w.r.t state vector
	term a, at point iq.  Storage for F_U is pre-allocated before
	the call.

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetConvFJump
extern int 
xf_EqnSetConvFJump(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
		   int nConv, const int *IParam, const real *RParam, 
		   int nq, const real *UL, const real *UR, real **AuxUL,
		   real **AuxUR, const real *n,  const real *xglob, 
		   const real *vg, real *F, real *F_UL, real *F_UR, real *C);
/*
PURPOSE:

  Evaluates an interior face jump flux for the convection residual
  term(s) stored in ResTerm, at nq points.  When an interface velocity
  is given, the jump flux takes vg into account by using F-U*vg
  instead of F in the (approximate Riemann) flux calculation.

INPUTS:

  EqnSet  : equation set structure
  ResTerm : the convection residual term(s) to evaluate
  nConv   : number of residual terms
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of times to calculate the jump flux
  UL      : nq Left state vectors, unrolled along StateRank first
  UR      : nq Right state vectors, unrolled along StateRank first
  AuxUL   : pointers to any requested auxiliary vectors on Left.  
            Each vector is interpolated at the nq points.
  AuxUR   : pointers to any requested auxiliary vectors on Right.
            Each vector is interpolated at the nq points.
  nq      : number of normal vectors (if 1 means normal is const)
  n       : nq normal vectors, unrolled along dim first
  xglob   : global positions of the nq points
  vg      : interface velocity at nq points [nq*dim]

OUTPUTS: 

  F   : nq flux vectors, unrolled along StateRank first, then along nq.
        Specifically,

	F[StateRank*iq + k]

	stores the kth component of the flux at point iq. Storage for F
	is pre-allocated before the call.

  F_UL : Flux Jacobian w.r.t UL.  Only computed if F_UL != NULL.  
  F_UR : Flux Jacobian w.r.t UR.  Only computed if F_UR != NULL.  

        Each stores flux Jacobians at nq points, unrolled along
        StateRank^2 first, then nq.  Specifically (e.g. for F_UL),
	  
	F_UL [StateRank^2*iq + StateRank*k + a] 

	stores the derivative of the kth flux term w.r.t state vector
	term a, at point iq.  Storage for F_UL and F_UR is
	pre-allocated before the call.
  C    : Connectivity value at each quad point (if not passed in as NULL)

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetConvFBC
extern int 
xf_EqnSetConvFBC(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, int nConv, 
		 const xf_BC *BC, const int *IParam, const real *RParam, 
		 int nq, const real *UI, const real *UB, const real *UB_UI, 
		 real **AuxUB, const real *n, const real *xglob, const real *vg,
		 real *F, real *F_U, real *C);
/*
PURPOSE:

  Evaluates a boundary face jump flux for the convection residual
  term(s) stored in ResTerm, at nq points.  When a boundary velocity
  is given, the returned flux is modified to take into account this
  motion via Fnew = F - UB*vg.

INPUTS:

  EqnSet  : equation set structure
  ResTerm : the convection residual term(s) to evaluate
  nConv   : number of residual terms
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of times to calculate the boundary flux
  UI      : nq state vectors, unrolled along StateRank first (Interior)
  UB      : nq state vectors, unrolled along StateRank first (Boundary)
  UB_UI   : nq srxsr derivative matrices (boundary state w.r.t interior)
  AuxUB   : pointers to any requested auxiliary vectors.  Each vector
            is interpolated at the nq points.
  n       : nq normal vectors, unrolled along dim first
  xglob   : global positions of the nq points
  vg      : boundary velocity at nq points [nq*dim]


OUTPUTS: 

  F   : nq flux vectors, unrolled along StateRank first, then along nq.
        Specifically,

	F[StateRank*iq + k]

	stores the kth component of the flux at point iq. Storage for F
	is pre-allocated before the call.

  F_U : only computed if F_U != NULL.  Stores flux Jacobians at nq
        points, unrolled along StateRank^2 first, then nq.
        Specifically,
	  
	F_U [StateRank^2*iq + StateRank*k + a] 

	stores the derivative of the kth flux term w.r.t state vector
	term a, at point iq.  Storage for F_U is pre-allocated before
	the call.
  C    : Connectivity value at each quad point (if not passed in as NULL)


RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetDiffA
extern int 
xf_EqnSetDiffA(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
	       const int *IParam, const real *RParam, const real *ResMetric, 
	       int nq, const real *U, const real *gU, const real *W, 
	       real *AW, real *A, real *A_UW, enum xfe_Bool *ConstA, real *mu);
/*
PURPOSE:

  Evaluates the diffusion "A matrix" corresponding to ResTerm, at
  state U, at nq points.  The A matrix multiplies grad U in the
  discretization to yield the visous flux.  A non-Null W means that
  the products AW = A*W, and A_UW = A_U*W are desired.  All outputs
  are incremented, not set.

INPUTS:

  EqnSet  : equation set structure
  ResTerm : the diffusion residual term to evaluate
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  ResMetric : resolution metric [length]; used for artificial stabilization
  nq      : number of times to calculate A
  U       : nq state vectors, unrolled along StateRank first
  gU      : nq state vector gradients, unrolled along StateRank first
  W       : supplementary vector for matrix-vector products. 
            Also at nq points, unrolled along StateRank first.

OUTPUTS: 

  AW   : dim*nq vectors, unrolled along sr(= StateRank) first, then 
         along nq, then along dim,

 	 AW{i,q,k} = AW[i*nq*sr + q*sr + k]

	 If W is grad U, AW is the viscous flux vector.
  A    : dim*dim*nq matrices, each sr x sr:

         A{i,j,q,k,l} = A[(i*dim+j)*nq*sr*sr + q*sr*sr + k*sr + l]

	 This is the no-frills attached matrix A
  A_UW : dim*nq matrices, each sr x sr:
  
         A_UW{i,q,k,a} = A_UW[i*nq*sr*sr + q*sr*sr + k*sr + l]
	               = A_U{i,j,q,k,l;a}* * W{j,q,l}  (sum over j, l)

	 Only computed if A is not const w.r.t U.

  ConstA : if True, means A is not a const w.r.t U.  In this case, A_UW
           will have been computed, assuming it was not passed in as NULL.

  mu : viscosity (or some measure thereof) at nq points.  Used for
       line creation.  This variable is added-to.

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetDiffFBC
extern int 
xf_EqnSetDiffFBC(const xf_EqnSet *EqnSet, const xf_BC *BC, const int *IParam, 
		 const real *RParam, int nq, const real *wn, const real *xglob, 
		 const real *gUI, real *gUB, real *gUB_gUI, 
		 enum xfe_Bool *AFlag, int *SetIndex);
/*
PURPOSE:

  Evaluates a diffusion boundary face jump flux at nq points.

INPUTS:

  EqnSet  : equation set structure
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of times to calculate the boundary flux
  gUI     : interior state gradient at the nq points, unrolled along
            StateRank first, then nq, then dim.
  wn      : normal vectors at nq points, unrolled along dim first
  xglob   : global positions of the nq points

OUTPUTS: 

  gUB     : If AFlag == False, boundary state gradient at nq points,
            unrolled just like gUI; if AFlag == True, viscous flux
            components (see below)
  gUB_gUI : derivatives of gUB w.r.t gUI; size = nq * dim*dim*sr*sr
  AFlag   : True means BC imposed condition on visc flux itself, not on
            grad u.  In this case, the entries of guB for which
            SetIndex is 1 are the set viscous flux components (already
            dotted with the normal); the remaining entries of the
            viscous flux vector are taken from AB*guI.  AFlag == False
            means that the BC imposed a condition on grad u; SetIndex
            is not used in this case.
  SetIndex: Only valid if AFlag == True.  SetIndex[k] == 1 means kth
            component of the diffusion flux is set by the bc.
  C    : Connectivity value at each quad point (if not passed in as NULL)

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetSourceS
extern int 
xf_EqnSetSourceS(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
		 int nSource, const int *IParam, const real *RParam, int nq, 
		 const real *U, const real *gU, real **qAuxU, const real *xglob, 
		 const real *pTime, real *S, real *S_U, real *S_gU, 
		 enum xfe_Bool *Nonzero_gu);
/*
PURPOSE:

  Evaluates a source term(s) stored in ResTerm, at nq points.

INPUTS:

  EqnSet  : equation set structure
  ResTerm : the source residual term(s) to evaluate
  nSource : number of residual terms
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of times to calculate the source term
  U       : nq state vectors, U[iq*sr + k]
  gU      : nq state vector gradients, gU[i*nq*sr + iq*sr + k] 
  AuxU    : pointers to any requested auxiliary vectors.  Each vector
            is interpolated at the nq points.
  xglob   : global positions of the nq points
  pTime   : pointer to simulation time


OUTPUTS: 

  S   : nq source term vectors, unrolled along StateRank first, then along nq.
        Specifically,

	S[StateRank*iq + k]

	stores the kth component of the flux at point iq. Storage for S
	must be pre-allocated before the call.

  S_U : only computed if S_U != NULL.  Stores derivatives of S at nq
        points, unrolled along StateRank^2 first, then nq.
        Specifically,
	  
	S_U [StateRank^2*iq + StateRank*k + a] 

	stores the derivative of the kth source term w.r.t state
	vector term a, at point iq.  Storage for S_U must be
	pre-allocated before the call.

 S_gU : only computed if S_gU != NULL.  Stores derivatives of S w.r.t
        gU at nq points.  i.e. the derivative of S{iq,k} w.r.t
        gU{i,iq,a} is
	  
	S_gU [i*nq*sr2 + iq*sr2 + k*sr + a]

	Storage for S_U must be pre-allocated.  Only computed if
	source term depends on the gradient of the state.  The calling
	function is notified of such a case by the flag Nonzero_gu.

  Nonzero_gu : set to True if S_gU was calculated; i.e. if S depends on gU.

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetScalar
extern int 
xf_EqnSetScalar(const xf_EqnSet *EqnSet, const char *Name, const int *IParam, 
		const real *RParam, int nq, const real *U, const real *gU, 
		real *s, real *s_U, int *nScalar, char ***pScalarNames, const real Gamma);
/*
PURPOSE:

  Computes a scalar quantity, s, and possibly a linearization, s_U,
  given a state vector, U, at nq points.

INPUTS:

  EqnSet  : equation set structure
  Name    : string containing the name of the scalar to calculate
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of points
  U       : state vector at the nq points, unrolled along
            StateRank first, then nq
  gU      : state gradients at all points, gU[i*nq*sr+iq*sr+k], optional

OUTPUTS: 

  s       : calculated scalar, at nq points (size nq)
  s_U     : linearization of scalar w.r.t state (size nq*StateRank),
            unrolled along StateRank first, then nq
  nScalar : number of scalars
  pScalarNames : names of scalars (do not deallocate this).  Note, if
                 these last two outputs are requested, no other outputs 
		 are returned.

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetVector
extern int 
xf_EqnSetVector(const xf_EqnSet *EqnSet, const char *Name, const int *IParam, 
		const real *RParam, int nq, const real *U, real *v, real *v_U);
/*
PURPOSE:

  Computes a vector quantity, v, and possibly a linearization, v_U,
  given a state vector, U, at nq points.

INPUTS:

  EqnSet  : equation set structure
  Name    : string containing the name of the vector to calculate
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of points
  U       : state vector at the nq points, unrolled along
            StateRank first, then nq

OUTPUTS: 

  v       : calculated vector, at nq points (size nq*dim) 
            v[0], v[1], (v[2]) is the vector at the first point
  v_U     : linearization of vector w.r.t state (size nq*dim*StateRank),
            unrolled along StateRank first, then dim, then nq

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetVariableChange
extern int 
xf_EqnSetVariableChange(const xf_EqnSet *EqnSet, const char *Name, const int *IParam, 
			const real *RParam, int nq, const real *U, real *V,
			int *nVSet, char ***pVSetNames);
/*
PURPOSE:

  Transforms variables in U from conservative to a desired other set,
  indicated by Name.

INPUTS:

  EqnSet  : equation set structure
  Name    : string containing the name of the variable set to change to
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of points
  U       : conservative-variable state vector at the nq points, unrolled
            along StateRank first, then nq

OUTPUTS: 

  V       : transformed-variable state vector at the nq points, unrolled 
            along StateRank first, then nq
  nVSet   : number of variable sets
  pScalarNames : names of variable sets (do not deallocate this).  Note, if
                 these last two outputs are requested, no other outputs 
		 are returned.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetScaleState
extern int 
xf_EqnSetScaleState(const xf_EqnSet *EqnSet, const xf_IC *ICOrig, const int *IParam, 
		    const real *RParam, int nq, real *U);
/*
PURPOSE:

  Scales variables in U based on new initial conditions.  Useful for
  restarting during parameter sequencing.

INPUTS:

  EqnSet  : equation set structure
  ICOrig  : original IC
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of points
  U       : unrolled state rank, with nq points

OUTPUTS: 

  U       : scaled state vector at the nq points (overwritten)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetAlterState
extern int 
xf_EqnSetAlterState(const xf_EqnSet *EqnSet, const char *AlterFcn, const int *IParam, 
		    const real *RParam, const real *xglob, const real *FcnParam, 
		    int nq, real *U);
/*
PURPOSE:

  Alters state U based on a prescribed alteration function, AlterFcn.

INPUTS:

  EqnSet  : equation set structure
  AlterFcn: which function to use for state change
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  xglob   : global coordinates at all points
  FcnParam: parameters specific to the function
  nq      : number of points
  U       : unrolled state rank, with nq points

OUTPUTS: 

  U       : altered state vector at the nq points (overwritten)

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetUpdateFraction
extern int 
xf_EqnSetUpdateFraction(const xf_EqnSet *EqnSet, const int *IParam, const real *RParam, 
			int nq, const real *u, const real *du, real *frac);
/*
PURPOSE:

  Evaluates maximum allowable update fraction that will keep u+frac*du
  physical (and possibly within a certain percentage of u) over nq
  points.  u and du are stored unrolled, with EqnSet->StateRank values
  for each of the nq points.

INPUTS:

  EqnSet  : equation set structure
  IParam  : vector of integer params as requested by EqnSet->IParamKey
  RParam  : vector of real params as requested by EqnSet->RParamKey
  nq      : number of points at which to check the update fraction
  u       : nq state vectors, unrolled along StateRank first
  du      : nq state update vectos


OUTPUTS: 

  frac    : maximum allowable update fraction

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetMaxCharSpeed
extern int 
xf_EqnSetMaxCharSpeed(const xf_EqnSet *EqnSet, int nq, const real *U, 
		   real **AuxU, const real *xglob, const int *IParam,
		   const real *RParam, real *vmax, real *dtmax, real *v, real *v_U);
/*
PURPOSE:

  Evaluates maximum characteristic speed over an element, by
  evaluating at the nq quad points.

INPUTS:

  EqnSet  : equation set structure
  nq      : number of points at which to check the max char speed
  U       : nq state vectors, unrolled along StateRank first
  AuxU    : pointers to any requested auxiliary vectors.  Each vector
            is interpolated at the nq points.
  xglob   : global positions of the nq points
  IParam  : integer parameters requested by the equation set
  RParam  : real parameters requested by the equation set


OUTPUTS: 

  vmax    : maximum characteristic speed
  dtmax   : maximum time step (set to > 0 to use this)
  v       : max characteristic speed at each point (optional)
  v_U     : linearization of v w.r.t U (optional)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetBCIsWall
extern int 
xf_EqnSetBCIsWall(const xf_BC *BC, enum xfe_Bool *IsWall);
/*
PURPOSE:

  Determines if BC is a Wall (e.g. for wall distance calculations)

INPUTS:

  BC      : boundary condition structure

OUTPUTS: 

  (*IsWall) : true if BC corresponds to a wall

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetPenaltyTerm
extern int 
xf_EqnSetPenaltyTerm(const xf_EqnSet *EqnSet, const real *RParam, 
                     const int *IParam, int nq, real *u, real *pPe, 
                     real *Pq, real *Pe_u, real *Pe_uu);

/*
 PURPOSE:
 
 Computes a penalty based on realizability constraints
 
 INPUTS:
 
 EqnSet  : equation set structure
 IParam  : vector of integer params as requested by EqnSet->IParamKey
 RParam  : vector of real params as requested by EqnSet->RParamKey
 nq      : number of points at which to check the update fraction
 u       : nq state vectors, unrolled along StateRank first
 
 OUTPUTS: 
 
 pPe     : total penalty value
 Pq      : penalty for each quadrature state
 Pe_u    : penalty gradient w.r.t. each of the quadrature states
 Pe_uu   : penalty Hessian w.r.t. each of the quadrature states
 
 RETURN:
 
 Error Code
 */


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetPerturbParam
extern int
xf_EqnSetPerturbParam(xf_EqnSet *EqnSet, xf_Sensitivity *Sensitivity, 
                      real *epsilon, enum xfe_Bool CorrectFlag);

/*
 PURPOSE:
 
 Perturbs a parameter defined in the xf_Sensitivity structure
 
 INPUTS:
 
 EqnSet      : equation set structure
 Sensitivity : structure containing sensitivity information
 epsilon     : if != NULL, it receives the value of the 
               eqnset-param-specific perturbation
 CorrectFlag : if True, the BCs and outputs are corrected by 
               epsilon
 
 OUTPUTS: 
 
 None: EqnSet gets modified
 
 RETURN:
 
 Error Code
 */



/******************************************************************/
//   FUNCTION Prototype: xf_EqnSetOutputDependentBCs
extern int 
xf_EqnSetOutputDependentBCs(xf_EqnSet *EqnSet, xf_KeyValue *Outputs,
                            real *Sensitivity, real *epsilon, 
                            enum xfe_Bool *Converged);
/*
PURPOSE:

  Sets output-dependent BCs using outputs given in "Outputs"

INPUTS:

  EqnSet : contains the BCs that are modified
  Outputs : key-value list containing the necessary outputs
  Sensitivity: if != NULL the bc will be corrected using 
               the sensitivity
  epsilon : Value of the BC parameter perturbation. 
            If Sensitivity == NULL, then epsilon gets calculated 
            according to the appropriate EqnSet
OUTPUTS: 
 
  Converged: used by certain output-dependent BCs.
  BCs in EqnSet are set

RETURN:

  Error Code
*/



#endif // end ifndef _xf_EqnSetHook_h
