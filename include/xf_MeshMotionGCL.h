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

#ifndef _xf_MeshMotionGCL_h
#define _xf_MeshMotionGCL_h 1

/*
  FILE:  xf_MeshMotionGCL.h

  This file contains the headers for functions dealing with the
  Geometric Conservation Law for mesh motion.

*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindMeshMotionGCLVector
extern int 
xf_FindMeshMotionGCLVector( xf_All *All, xf_Vector **pGCL);
/*
PURPOSE:

  Finds a GCL vector based purely on a title search.

INPUTS:

  All : all structure

OUTPUTS: 

  (*pGCL) : pointer to GCL vector

RETURN:

  Error Code: xf_NOT_FOUND if not found

*/


/******************************************************************/
//   FUNCTION Prototype: xf_InitMeshMotionGCLVector
extern int 
xf_InitMeshMotionGCLVector( xf_All *All, xf_Vector *U, int ind, 
                            enum xfe_Bool InitFlag, xf_Vector **pGCL);
/*
PURPOSE:

  Finds/creates and initializes a GCL vector consistent (in terms of
  basis and order) with the passed-in state, U.

INPUTS:

  All : all structure
  U : state vector (need this for order/basis information)
  ind : if >= 0, suffix on data name (to distinguish between multiple similar ones)
  InitFlag : if True, GCL vector will be initialized to ones even
             if the GCL vector is found

OUTPUTS: 

  (*pGCL) : pointer to GCL vector

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindMeshMotionGCLLinearization
extern int 
xf_FindMeshMotionGCLLinearization( xf_All *All, const char *Root, 
				   int ind, xf_Vector **pRoot_GCL);
/*
PURPOSE:

  Finds a vector used for linearizing some quantity (e.g. an output)
  whose name is in Root, with respect to the GCL state.

INPUTS:

  All : all structure
  Root : string containing name of quantity being linearized
  ind : index appended to name of vector in data set

OUTPUTS: 

  (*pRoot_GCL) : pointer to Root_GCL vector

RETURN:

  Error Code: xf_NOT_FOUND if not found

*/



/******************************************************************/
//   FUNCTION Prototype: xf_InitMeshMotionGCLLinearization
extern int 
xf_InitMeshMotionGCLLinearization( xf_All *All, xf_Vector *U, const char *OutputName,
                                   int ind, enum xfe_Bool InitFlag, xf_Vector **pJ_GCL);
/*
PURPOSE:

  Finds/Initializes a vector used for linearizing some quantity
  (e.g. an output) whose name is in Root, with respect to the GCL
  state.

INPUTS:

  All : all structure
  U : state or state-like vector (containing basis/order/etc. info)
  OutputName : output for which we want the linearization
  ind : index appended to name of vector in data set
  InitFlag : if True, linearization will be initialized to zero even
             if it is found

OUTPUTS: 

  (*J_GCL) : pointer to the linearization vector

RETURN:

  Error Code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_InitMeshMotionGCLAdjoint
extern int 
xf_InitMeshMotionGCLAdjoint( xf_All *All, xf_Vector *Psi, int ind, 
                             enum xfe_Bool InitFlag, xf_Vector **pPsiGCL);
/*
PURPOSE:

  Finds/creates and initializes a GCL adjoint vector consistent (in
  terms of basis and order) with the passed-in adjoint, Psi.

INPUTS:

  All : all structure
  Psi : adjoint vector (need this for order/basis/OutputName information)
  ind : if >= 0, suffix on data name (to distinguish between multiple similar ones)
  InitFlag : if True, GCL adjoint will be initialized to zero even
             if the GCL adjoint is found

OUTPUTS: 

  (*pPsiGCL) : pointer to GCL adjoint vector

RETURN:

  Error Code
*/




/******************************************************************/
//   FUNCTION Prototype: xf_MeshMotionGCLResidual
extern int 
xf_MeshMotionGCLResidual( xf_All *All, xf_Vector *GCL, enum xfe_Bool ZeroFlag, xf_Vector *RGCL);
/*
PURPOSE:

  Calculates residual for marching the GCL equation in time.

INPUTS:

  All : all structure
  GCL : GCL Vector
  ZeroFlag : if True, RGCL is zeroed out before computation
             if False, residual is added to RGCL

OUTPUTS: 

  RGCL : residual 

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MeshMotionGCLSolveSystem
extern int 
xf_MeshMotionGCLSolveSystem( xf_All *All, real c, xf_Vector *S, xf_Vector *GCL);
/*
PURPOSE:

  Solves the following system for the GCL Vector:

    c*M*GCL + RGCL(GCL) + S = 0

  Note that RGCL does not depend on GCL, so this is an easy solve
  (just a mass matrix inversion).

INPUTS:

  All : all structure
  c   : a constant multiplying the mass matrix
  S   : source vector

OUTPUTS: 

  GCL : GCL Vector that solves the above equation

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_MeshMotionMap_GCL
extern int 
xf_MeshMotionMap_GCL( xf_Vector *GCL, int egrp, int elem, xf_BasisData *PhiData, 
		      int npoint, int dim, real *gb, real *gbigb_X);
/*
PURPOSE:

  Interpolates the GCL vector (gb = gbar) and spatial derivatives
  (gbigb_X = 1/gbar * gbar_X).

INPUTS:

  GCL : discrete vector of unknowns in GCL vector approximation
  egrp, elem: element in question
  PhiData : basis functions + possibly gradients at desired points
  npoint : number of points
  dim : dimension

OUTPUTS: 

  gb : interpolated GCL vector [npoint]
  gbigb_X : interpolated GCL gradient * 1/gbar [npoint*dim]
            (dim = fastest-running index)

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DGTimeFindGCLVectors
extern int
xf_DGTimeFindGCLVectors(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                        xf_Vector ***pGCLj, xf_Data ***pDj);
/*

PURPOSE:

  Looks for GCL vectors at time nodes for a DG-in-time discretization.
  Note, (*pGCLj) should eventually be released by the calling
  function, and same with (*pDj).  Both are optional outputs.
  
INPUTS:

  All: All data struct
  TimeScheme: Timescheme being used; e.g. DG1 or DG2
  
OUTPUTS:

  (*pGCLj) : (allocated) vector of pointers to found or allocated (OrderTime+1) GCL vectors
             (optional)
  (*pDj)   : (allocated) corresponding pointers to the data within All->DataSet (optional)
  
RETURN:

  Error Code
  
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DGTimeInterpolateGCL
extern int
xf_DGTimeInterpolateGCL(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                        real tq, real *phiout);

/*

PURPOSE:

  Interpolates the already-computed GCLTime vector to the reference-space time tq.
  
INPUTS:

  All: All data struct
  TimeScheme: Timescheme being used; e.g. DG1 or DG2
  tq: Time to interpolate GCL at
  
OUTPUTS:

  GCL: Sets GCL = GCL(tq) in All data struct
  phiout : temporal basis functions at tq (optional)
  
RETURN:

  Error Code
  
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DGTimeUnsteadyResidualGCL
extern int
xf_DGTimeUnsteadyResidualGCL(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                             xf_Vector **GCLj, xf_Vector *GCLprev, xf_Vector *GCLtemp,
                             real Time, real TimeStep, xf_Vector **RGCLi);
/*
 PURPOSE: 

   Computes residual vectors associated with DG-in-time discretization
   of the GCL equation.  Similar to corresponding "U" function in
   xf_SolverDGTime.  Here, we make use of the above source GCL
   function and just add an extra temporal stiffness matrix term.

 INPUTS:
 
   All: All struct
   TimeScheme: Timescheme, e.g. DG1 or DG2
   GCLj : GCL vectors on current time slab
   GCLprev : GCL at end of previous time slab
   GCLtemp : temporary vector used in calculations
   Time : time at start of slab
   TimeStep: current timestep
   
 OUTPUTS:
 
   RGCLi: Filled residual at each temporal node
   
 RETURNS:
   
   Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DGInTimeGCLSolve
extern int 
xf_DGTimeGCLSolve(xf_All *All, enum xfe_TimeSchemeType TimeScheme, real Time, 
                  real TimeStep, xf_SolverData *SolverData, xf_Vector *GCLprev,
                  xf_Vector **GCLj);

/*
PURPOSE:

  Solves for the GCL coefficients at each temporal basis node associated with a 
  DGInTime scheme, for all elements. For an order-r DGInTime scheme, there are 
  (nn)x(r+1) of these coefficients for each element. They are computed by taking
  
    |GCL_0| =            |[M]^{-1}*S_0(x,t) + GCLprev|
    |GCL_1| =  [K]^{-1}  |[M]^{-1}*S_1(x,t)          |
                               ...
    |GCL_r| =            |[M]^{-1}*S_r(x,t)          |
    
  where each GCL_i is a length-[nn] vector and [K] is an [(r+1)nn] x [(r+1)nn]
  block matrix, with (r+1) unique values. The [K] values are calculated 
  analytically for DG1 and DG2 schemes from the temporal stiffness matrix AT, 
  and the [K] matrix itself is not actually formed. [M] is the element mass 
  matrix and S_0(x,t) is a source term related to the spatial part of the GCL PDE.


INPUTS:

  All: Global All struct
  TimeScheme: e.g. DG1 or DG2
  Time: current time
  TimeStep: current time step size
  SolverData: (not currently used)
  GCLprev: previous GCL

OUTPUTS:

  GCLj: solved GCL vector on time slab
  
RETURN:

  Error Code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_DGinTimeStepGCLAdjoint
extern int
xf_DGinTimeStepGCLAdjoint(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
                          real Time, real TimeStep, int iSlab, int nSlab,
                          xf_SolverData *SolverData, xf_Vector **Uj,
                          int nPsi,  xf_Vector ***Psii, xf_Vector **PsiGCLnext, 
                          xf_Vector ***PsiGCLi);
/*
PURPOSE:

  Takes a step of a DG in time unsteady adjoint solver for the GCL
  variable.  On input, the latest GCL adjoint (left node of future
  adjacent time slab) is stored in PsiGCLnext.  On output, the entire
  GCL adjoint in the current time slab is computed and stored in
  PsiGCLi.
  
INPUTS:

  All         : all structure
  TimeScheme  : what time scheme to use
  Time        : current time
  TimeStep    : delta t
  iSlab       : current time slab
  nSlab       : number of slabs
  SolverData  : solver data structure (contains CFL)
  Uj          : state vectors on current time slab
  nPsi        : number of adjoints
  Psii        : ALL the state adjoints on the current time slab
  PsiGCLnext  : GCL adjoint on left node of next time slab
  
OUTPUTS: 

  PsiGCLi     : ALL the GCL adjoints on the current time slab

RETURN: Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_MeshMotionGCLSolveAdjoints
extern int
xf_MeshMotionGCLSolveAdjoints(xf_All *All, real c, real d, int nPsi, 
			      xf_Vector **SGCL, xf_Vector **PsiGCL);
/*
PURPOSE:

  GCL version of xf_SolveAdjoints.  Solves, for each of the nPsi
  adjoint equations,

    (c*M + 0)^T PsiGCL + SGCL + d*J_GCL = 0

  where M is the mass matrix, the spatial Jacobian is zero, SGCL is an
  optional source term, and J_GCL is the linearization of the output
  with respect to the GCL state.  This function is only used in an
  unsteady calculation.

  Note, each GCL adjoint vector, PsiGCL, should have an associated
  output stored in its OutputName field.

INPUTS:

  All: All structure
  c : real constant that multiplies the Mass matrix added to the
      Jacobian (zero in the GCL case).
  d : real constant that multiplies the output linearization before it
      is added to the residual. 
  nPsi : number of adjoints
  SGCL : vector of source vectors for the adjoint equations (optional)

OUTPUTS:

  PsiGCL : vector of adjoint solution vectors. These vectors are overwritten,
           so that they do not need to be initialized on input

RETURN:

  Error code

*/


#endif // end ifndef _xf_MeshMotionGCL_h
