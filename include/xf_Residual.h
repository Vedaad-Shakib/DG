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

#ifndef _xf_Residual_h
#define _xf_Residual_h 1

/*
  FILE:  xf_Residual.h

  This file contains the headers for the residual-calculation functions

*/

#include "xf_SolverStruct.h"
#include "xfYu_Model.h"

/******************************************************************/
//   FUNCTION Prototype: xf_FindSupportedVector
extern int
xf_FindSupportedVector(xf_All *All, const char Title[], xf_Vector **pV);
/*
PURPOSE:

  Searches for a supported vector with name Title
  (e.g. "WallDistance") in All.  If not found, attempts to create it.

INPUTS:
 
  All   : All structure
  Title : Name of vector to create

OUTPUTS: 

  (*pV)  : pointer to located/created vector

RETURN:

  Error Code
*/
  

/******************************************************************/
//   FUNCTION Prototype: xf_RetrieveFcnParams
extern int
xf_RetrieveFcnParams(xf_All *All, const xf_EqnSet *EqnSet, int **pIParam, 
		     real **pRParam, int *pnAuxU, xf_Vector ***pAuxU);
/*
PURPOSE:
 
  Retrieves integer and real parameters specified in
  EqnSet->IParamKey, and EqnSet->RParamKey.  These parameters are then
  passed to eqnset-specific functions.

INPUTS:

  All     : All->Param->KeyValue contains additional key value list 
  EqnSet  : the equation-set structure, with its own EqnSet->Param

OUTPUTS:

  pIParam : pointer to IParam list, gets reallocated
  pRParam : pointer to RParam list, gets reallocated
  (*pnAuxU) : number of auxiliary vectors requested by EqnSet
  (*AuxU) : pointer to auxiliary vectors from All->DataSet

RETURNS:

  Error code

*/
  

/******************************************************************/
//   FUNCTION Prototype: xf_SortEqnSetBCs
extern int
xf_SortEqnSetBCs(const xf_Mesh *Mesh, xf_BCs *BCs);
/*
PURPOSE:
 
  Arranges BCs->BC so that they match the boundary face groups in
  Mesh.  The ordering in BCs is changed if it was out of order.

INPUTS:

  Mesh: mesh structure
  BCs: boundary conditions structure

OUTPUTS:

  None: BCs->BC are rearranged if out of order

RETURNS:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ModMotionPreEqnCall
extern void
xf_ModMotionPreEqnCall(int nq, int dim, int sr, xf_MotionData *MD,
		       real *u, real *wn);
/*
PURPOSE: 

  Converts state and normal to physical space prior to a call to an
  EqnSet-specific function.

INPUTS: 
 
  nq   : number of points
  dim  : dimension
  sr   : state rank
  MD   : mesh motion data structure
  u    : reference-space state [nq*sr]
  wn   : reference-space normals [nq*dim]

OUTPUTS:

  u, wn : modified versions

RETURNS: Error code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ModMotionPhysGrad
extern void
xf_ModMotionPhysGrad(int nq, int dim, int sr, xf_MotionData *MD,
		     real *u, real *gu);
/*
PURPOSE: 

  Converts gradient from ALE reference space to physical space.

INPUTS: 
 
  nq   : number of points
  dim  : dimension
  sr   : state rank
  MD   : mesh motion data structure
  u    : physical-space state (necessary due to product rule)
  gu   : reference-space gradient [dim*nq*sr];

OUTPUTS:

  gu   : modified gradient (now in physical space)

RETURNS: Error code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ModMotionPostEqnCall
extern void
xf_ModMotionPostEqnCall(int nq, int dim, int sr, xf_MotionData *MD,
			real *u, real *wn);
/*
PURPOSE: 

  Converts state and normal from physical space to ref space after a
  call to an EqnSet-specific function.

INPUTS: 
 
  nq   : number of points
  dim  : dimension
  sr   : state rank
  MD   : mesh motion data structure
  u    : physical-space states [nq*sr]
  wn   : physical-space normals [nq*dim]

OUTPUTS:

  u, wn : modified versions

RETURNS: Error code
*/


/******************************************************************/
//   FUNCTION Prototpye: xf_CalculateResidualBFaces
extern int
xf_CalculateResidualBFaces(xf_All *All, const xf_Vector *U, xf_Vector *R, 
			   xf_JacobianMatrix *R_U, xf_OutputEvalData *OutputEval,
			   int nBFG, int *BFGs, xf_SolverData *SolverData);
/*

PURPOSE:

  Calculates the boundary contribution to the spatial residual, R,
  using All->EqnSet and the state vector U.  The Jacobian matrix, R_U,
  is computed if not passed as NULL.  In addition a flux output, can
  be calculated if OutputEval is not passed in as NULL.  This flux
  output consists of the integrated convective and diffusive fluxes
  weighted by a given FluxWeights and possibly a FluxMoments vector
  (size StateRank) to form a scalar.


INPUTS:

  All: All structure
  U : state vector
  OutputEval : output evaluation structure (see xf_SolverStruct.h).
               If not NULL, output evaluation is assumed.
  nBFG : number of boundary face groups over which to integrate
         (only relevant if BFGs is not NULL)
  BFGs : vector of boundary face groups for the integration (if NULL,
         all boundary face groups will be considered)
  SolverData: contains pointers to solver variables such as the
              connectivity vector or the regularity vector.  Some
	      of these (e.g. connectivity) may be modified.

OUTPUTS:

  R_U : Jacobian matrix
  R : Residual vector(s)

RETURN:

  Error code

*/


/******************************************************************/
//function: Implicit Time Integration for Detailed Chemistry
extern int
xf_CalculateResidualElems_DetailChem(xf_All *All, Yu_Model *pModel, xf_Vector *U, 
                                     xf_Vector *R, xf_JacobianMatrix *R_U, 
                                     xf_SolverData *SolverData);
/******************************************************************/
//   FUNCTION Prototype: xfYu_CalculateResidual
extern int
xfYu_CalculateResidual(xf_All *All, Yu_Model *Model, xf_Vector *U, xf_Vector *R, 
		     xf_JacobianMatrix *R_U, xf_SolverData *SolverData);

extern int
xf_CalculateResidual(xf_All *All, xf_Vector *U, xf_Vector *R, 
		     xf_JacobianMatrix *R_U, xf_SolverData *SolverData);

/******************************************************************/
//   FUNCTION Prototype: xf_CalculateResidualLeanElem
extern int
xf_CalculateResidualLeanElem(xf_All *All, int egrp, int elem, xf_Vector *U, 
			     real *ER, real *ER_EU, real **ER_NU, real **NR_EU, 
			     xf_JacobianMatrix *R_U, xf_SolverData *SolverData);
/*

PURPOSE:

  Calculates the spatial residual for one element, using a memory-lean
  approach.  This is computationally more expensive when aggregated
  over all the elements in the sense that face residual calculations
  are duplicated for neighboring elements.  All output vectors and
  matrices are set (not added to) in this function.  Note that not all
  outputs are required, but if given must be preallocated.

  Also note that parallel communication is assumed to have taken place
  before the call to this function.  That is, no halo exchanges take
  place here.  In the current setup where this function is called from
  the linear solver, the linear solver functions take care of the halo
  exchange.

INPUTS:

  All: All structure
  egrp, elem : element in question
  R_U : Jacobian matrix -- no values allocated, just for ranks and connectivity
  U : entire state vector (will need Basis, Order, etc.)

OUTPUTS:

  ER : residual on elem (optional)
  ER_EU : self jacobian on elem (optional)
  ER_NU : ER_NU[j] = jacobian of ER w.r.t neighboring elem j (optional)
          must be preallocated with adequate size!
  NR_EU : NR_EU[j] = jacobian of neighbor elem j residual w.r.t EU (optional)
          must be preallocated with adequate size!

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_MeshMotionGCL_PsiTxR_G
extern int
xf_MeshMotionGCL_PsiTxR_G(xf_All *All, xf_Vector *U, int nPsi, 
                          xf_Vector **Psi, xf_Vector **SGCL);
/*
PURPOSE:

  Calculates Psi^T * (Residual_gbar) using finite difference
  approximations.  The residual linearization with respect to gbar is
  calculated by perturbing gbar on each element and its neighbors in
  turn.  Multiple adjoints can be passed in via Psi, resulting in
  multiple incremented sources in SGCL.

INPUTS:

  All: All structure
  U  : state vector
  nPsi : number of adjoints
  Psi : adjoint vectors

OUTPUTS:

  SGCL : vector of source vectors that get incremented

RETURN:

  Error code

*/



#endif // end ifndef _xf_Residual_h
