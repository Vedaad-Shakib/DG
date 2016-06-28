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

#ifndef _xf_Solver_h
#define _xf_Solver_h 1

/*
  FILE:  xf_Solver.h

  This file contains the headers for the top-level solver functions

*/

#include "xf_SolverStruct.h"
#include "xf_LinearSolverStruct.h"
#include "xf_SolverTools.h"
#include "xf_AdaptStruct.h"
#include "xfYu_Model.h"


/******************************************************************/
//   FUNCTION Prototype: xf_CreateSolverData
extern int 
xf_CreateSolverData( xf_All *All, xf_SolverData **pSolverData);
/*
PURPOSE:

  Creates (allocates) a solver data structure.

INPUTS:

  All : all structure

OUTPUTS: 

  (*pSolverData) : allocated structure

RETURN: Error code

*/

/******************************************************************/
//   FUNCTION Definition: xf_LimitPoints
extern int 
xf_LimitPoints( xf_Mesh *Mesh, int egrp, int elem, xf_Vector *U,
               xf_EqnSet *EqnSet, int *pnq, real **pxq);
/*
 PURPOSE:
 
 Creates (allocates) a set of quadature points and limit points on the faces.
 
 INPUTS:
 
 Mesh : mesh structure
 egrp, elem: element coordinates
 U : state vector
 EqnSet : EqnSet structure
 *pnq : pointer to the number of points
 **pxq : pointer to the points coordinates in reference space
 
 OUTPUTS: 
 
 None : pnq and pxq are modified
 
 RETURN: 
 
 Error code
 
 */


/******************************************************************/
//   FUNCTION Definition: xf_DestroySolverData
extern int 
xf_DestroySolverData( xf_SolverData *SolverData);
/*
PURPOSE:

  Destroys a solver data structure.

INPUTS:

  (*SolverData) : solverdata structure to destroy


OUTPUTS: 


RETURN: Error code

*/

/******************************************************************/
//   FUNCTION Definition: xf_CalculateArtificialTimeStep
extern int
xf_CalculateArtificialTimeStep(xf_All *All, xf_Vector *U, real CFL, 
                               xf_Vector *dt);
/*
 PURPOSE:
 
 Calculates a vector of dt's based on the CFL number and if 
 local time-step is requested
 
 INPUTS:
 
 All : All structure
 U   : Primal state
 CFL : Value of the CFL
 
 OUTPUTS: 
 
 dt : vector gets filled with values. Create it before 
 calling this function
 
 RETURN: Error code
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_AddPenaltyMatrix
extern int
xf_AddPenaltyMatrix(xf_All *All, xf_Vector *U, xf_Vector *R, 
                    xf_JacobianMatrix *R_U, xf_SolverData *SolverData);
/*
 PURPOSE:
 
 Calculates and adds a penalization matrix to the Jacobian 
 to account for feasibility constraints
 
 INPUTS:
 
 All : All structure
 U   : Primal state
 R_U : Jacobian
 SolverData: Structure storing parameters of the solver
 
 OUTPUTS: 
 
 None, the Jacobian gets modified.
 
 RETURN: Error code
 
 */

/******************************************************************/
//   FUNCTION Prototype: xf_PreconditionerLeanCheck
extern int
xf_PreconditionerLeanCheck(enum xfe_PreconditionerType Preconditioner,
			   enum xfe_Bool *pLeanFlag);
/*
PURPOSE:

  Checks if Preconditioner is memory-lean

INPUTS:

  Preconditioner : the preconditioner type

OUTPUTS: 

  (*pLeanFlag) : true if preconditioner is memory lean

RETURN: Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SolveNonlinearSystem
extern int
xf_SolveNonlinearSystem(xf_All *All, real c, enum xfe_Bool LinearFlag,
			xf_Vector *S, xf_Vector **pU);
/*
PURPOSE:

  Solves:

                   c*M*U + R(U) + S = 0

  using an appropriate nonlinear solver, specified in All->Param.  For
  standard steady-state runs, pass c=0.  For unsteady runs, c should
  be the appropriate coefficient from the time discretization.

  For robustness, artificial time-stepping is used in the course of
  the nonlinear solve, so that at each nonlinear iteration, the
  system looks like

               (c + 1/dta)*M*U + R(U) + S = 0,

  where dta is the (possibly element-local) artificial time step.
  During the course of the nonlinear solve, 1/dta -> 0.

INPUTS:

  All: All structure
  c : constant in front of M*U product
  LinearFlag : If True, the Jacobian matrix will not be recalculated if
               it already exists.
  S : source vector (or NULL)
  (*U) : state vector

OUTPUTS:

  (*U) : modified state vector

RETURN:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_SolveAdjoints
extern int
xf_SolveAdjoints(xf_All *All, real c, real d, enum xfe_Bool ReuseJacobian,
                 xf_Vector *U, int nPsi, xf_Vector **S, xf_Vector **Psi,
                 xf_Vector **LinR, enum xfe_Bool CalcSolError, 
                 enum xfe_Bool SetZero);
/*
PURPOSE:

  Solves nAdjoint adjoint equations, each of the form:

    (c*M + R_U)^T Psi + S + d*J_U = 0

  where M is the mass matrix, R_U is the Jacobian, S is an optional
  source term, and J_U is the linearization of the output.  The
  constant d will usually be 1 for steady adjoint calculations, but
  can be <> 1 for weighting during an unsteady-adjoint solve.
  Including c*M and S enables use of this function during an unsteady
  adjoint solve.  Note that Psi and S (if specified) must each be
  pointers to nPsi vectors.

  U is provided as a state vector for which the Jacobian, R_U, will be
  computed.  If ReuseJacobian is set to True, the Jacobian will not be
  recomputed (assuming it exists).  However, a state vector U should
  still be provided for vector sizing purposes, and for calculating
  outputs (e.g. if any are nonlinear).

  Note, each adjoint vector, Psi, should have an associated output
  stored in its OutputName field.

  An optional linear residual output, LinR, is provided in cases where
  the linear adjoint residual is requested.  Only supported for one
  adjoint vector (will be for the last one if multiple adjoints exist).


INPUTS:

  All: All structure
  c : real constant that multiplies the Mass matrix added to the
      Jacobian. Note, if ReuseJacobian == True, R_U is not recomputed
      and hence M is not re-added -- thus c is irrelevant in this
      case.
  d : real constant that multiplies the output linearization before it
      is added to the residual. 
  ReuseJacobian : if True R_U is not recomputed, assuming it is found
                  in the DataSet.
  U : state vector about which R_U is computed; also used for
      calculating the output and sizing vectors

  nPsi : number of adjoints
  S : vector of source vectors for the adjoint equations (optional)
  CalcSolError : If True and there is a combined output, a solution 
                 error estimate will be computed following Hartmann's 
                 approach
  SetZero : If True, the initial guess of Psi is set to zero.

OUTPUTS:

  Psi : vector of adjoint solution vectors. These vectors must be
        initialized properly (e.g. set to zero or a good guess) before
	calling this function.
  LinR : pointer to linear residual vector (optional -- only one
         supported, see above)

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteTimeHist
extern int 
xf_WriteTimeHist(xf_TimeHistData *TimeHistData, const char *fname);
/*
PURPOSE:

  Writes Time history data to a text file, fname

INPUTS:

  TimeHistData : structure to write out
  fname        : name of file to write

OUTPUTS: 

  None, TimeHistData is written out
		   
RETURN: Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_ApplyTimeScheme
//extern int
//xfYu_ApplyTimeScheme(xf_All *All, Yu_Model *Model, const char *SavePrefix, 
//		     enum xfe_Bool RestartFlag, xf_Vector **pU0,
//		     xf_TimeHistData *TimeHistData);

extern int
xf_ApplyTimeScheme(xf_All *All, const char *SavePrefix, 
		   enum xfe_Bool RestartFlag, xf_Vector **pU0,
		   xf_TimeHistData *TimeHistData);
/*
PURPOSE:

  Runs the unsteady time-stepping solver.  Pulls info from All.

  Each implicit time scheme is written as:
   
  M/dt * (c0*u^{n+1} + c1*u^n + c2*u^{n-1} + ... ) + R(u^{n+1}) = 0

  The coefficients are particular to the unsteady scheme.


INPUTS:

  All: All structure
  SavePrefix : prefix for writing unsteady write-interval files.
  RestartFlag: True if the run is restarted, in which case the appropriate 
               time-index states should be present in All->DataSet.  For non-
	       restarted runs, additional vectors are created and initialized
	       to the time-index=0 state.
  pU0: initial/time-index=0 state (pointer passed to be reset in case there 
                                   is re[artitioning)

  TimeHistData : structure used to accumulate history of the time
                 values run and any outputs at all times.  Optional:
                 can pass in as NULL.  If not NULL, the structure must
                 have been created before the call, and the desired
                 number of Outputs must be stored in nOutput, with
                 corresponding names in OutputNames.  The Time vector
                 and the OutputValues vector are reallocated in this
                 function.

OUTPUTS:

  U0: final time state data is stored here
  TimeHistData : If given as not NULL, modified according to the above
                 description.

RETURN:

  Error code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ApplyTimeSchemeAdjoint
extern int
xf_ApplyTimeSchemeAdjoint(xf_All *All, const char *SavePrefix, 
			  xf_Vector *U0, int nPsi, xf_Vector **Psi, 
			  xf_TimeHistData *TimeHistData);
/*
PURPOSE:

  Runs the unsteady adjoint time-stepping solver.  The primal unsteady
  solver (ApplyTimeScheme) must have been run prior to this to
  generate the time history data (TimeHistData) and primal state .data
  files at each time step (for the nonlinear case).

  A discrete unsteady adjoint solve is performed, which means that the
  entire discrete unsteady Jacobian matrix is transposed.  This is
  best illustrated by a BDF1 method, in which the adjoint at each time
  step, Psi(i) is given by the solution to:

  
   [ D(1) -M/dt(2)   0.0    ] [Psi(1)]   [w(1)*J_U(1)]
   [ 0.0    D(2)   -M/dt(3) ] [Psi(2)] + [w(2)*J_U(2)] = 0
   [ 0.0    0.0      D(3)   ] [Psi(3)]   [w(3)*J_U(3)]
   \___________ ____________/            \_____ _____/
               V                               V
          [R(i)_U(j)]^T                Jtot_U(j) (see below)

  where,
  D(i) = [M/dt(i) + Rs_U(i)]^T
  dt(i) = t(i)-t(i-1) = the time step right before time i.
  Rs(U(i)) = steady-state residual at time i
  Note, M^T = M (mass matrix is symmetric)

  The calculated adjoints Psi(i) would be useful as-is if the total
  output from the unsteady simulation, Jtot, were a simple sum of the
  outputs at each time step.  Alternatively, if the output only
  depends on the final, t=T, state, then only J_U(T) would be nonzero
  in the above system.  If the output is a weighted sum of outputs at
  each of the different times (as described below), then an
  appropriate output weighting has to be specified in the TimeHistData
  structure.

  When this function returns, the values stored in Psi are not the
  "t=1" values of the adjoints.  Rather, Psi is set to
  [R(i)_U(0)]^T*Psi(i) (summation implied), where R(i)_U(0) is the
  linearization of the unsteady residual at time i w.r.t. the initial
  condition, U(0).  For the BDF1 example above, only the time 1
  residual, R(1), depends on U(0):

         R(1) = M/dt(1)*[U(1) - U(0)] + Rs(U(1))
  
  so that [R(i)_U(0)]^T = [-M/dt(1), 0, 0, ...].  The reason for this
  multiplication is that [R(i)_U(0)]^T*Psi(i) is precisely the vector
  (size of U(0)) that represents the sensitivity of the output to the
  initial condition U(0).  The following derivation makes this clear:

    Jtot = w(j)*J(U(j)),   w(j) = weight at time j, j >= 1
    
    Jtot_U(j) = w(j)*J_U(j)  (can be computed at each time step)

    Jtot_U(0)^T = [Jtot_U(j)]^T * U(j)_U(0)
                                  ^^^^^^^^^
    R(i) = 0  implies   R(i)_U(j) * U(j)_U(0) + R(i)_U(0) = 0
                                    ^^^^^^^^^
    Jtot_U(0)^T = -[Jtot_U(j)]^T * [R(i)_U(j)]^-1 * R(i)_U(0)
                  ~~~~~~~~~~~~~~V~~~~~~~~~~~~~~~~
                             Psi(i)^T

    where [R(i)_U(j)]^T * Psi(i) + Jtot_U(j) = 0

    Thus, 
    Jtot_U(0)^T = Psi(i)^T * R(i)_U(0)
    Jtot_U(0)   = [R(i)_U(0)]^T*Psi(i)}^T
    

  Note, for schemes like the trapezoidal method, in which the unsteady
  residual R(i) depends not only on Rs(U(i)) but also on Rs(U(i-1)),
  etc., the off-diagonal terms in the unsteady Jacobian, R(i)_U(j)
  will now have Rs_U terms.  These terms complicate the unsteady
  adjoint solution because products Rs_U(i)*Psi(i) are required at
  times when the full steady Jacobian, Rs_U(i), is not stored (we do
  not want to store two full Jacobians at once for memory reasons).
  Probably the best solution to this is to have the CalculateResidual
  function compute the product Rs_U*V for an arbitrary V on the fly --
  this should also be faster than forming a full Jacobian.  BDF2, with
  BDF1 on the first time step, does not have this complication; the
  bandwidth on R(i)_U(j) just increases.

INPUTS:

  All: All structure

  c : real constant that multiplies the Mass matrix added to the
      Jacobian. Note, if ReuseJacobian == True, R_U is not recomputed
      and hence M is not re-added -- thus c is irrelevant in this
      case.

  ReuseJacobian : if True R_U is not recomputed, assuming it is found
                  in the DataSet.

  U0 : a model state vector for sizing purposes.  In a linear-primal
       case, this state vector is used throughout the computation,
       getting passed into the residual calculation and the output
       calculation -- thus, a reasonable state vector should still be
       passed in.  In the nonlinear case, the state vector will be
       read in from disk at each time step.

  nPsi : number of adjoints

  Psi : vector of adjoint solutions; these serve as the initial
        conditions, and the final "t=0 adjoint solutions" are stored
        here.  See description above. Each adjoint vector should
        have an associated output stored in its OutputName field.

  TimeHistData : structure containing time history, generated by the
                 primal unsteady solve, ApplyTimeScheme. Required.
                 Time weights can be supplied through this structure.

OUTPUTS:

  Psi : vector of adjoint solution vectors

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ApplyTimeSchemeAdapt
extern int
xf_ApplyTimeSchemeAdapt(xf_All *All, const char *SavePrefix, 
			enum xfe_AdaptOnType AdaptOn,
			xf_TimeHistData *TimeHistData);
/*
PURPOSE:

  Calculates non-output-based adaptive indicator for unsteady
  simulations.  The type of indicator is governed by AdaptOn.

INPUTS:

  All: All structure
  SavePrefix: prefix for saving files
  AdaptOn : what type of adaptation are we performing?  For example,
            this could be interpolation error or residual.
  TimeHistData : structure containing time history, generated by the
                 primal unsteady solve, ApplyTimeScheme. Required.

OUTPUTS:

  None.  Data is stored in the All structure and/or written to files.

RETURN:

  Error code

*/



/******************************************************************/
// FUNCTION Prototype: xf_DGTimeScheme2Order
extern int
xf_DGTimeScheme2Order(enum xfe_TimeSchemeType TimeScheme, int *OrderTime);
/*

PURPOSE:

  Determines Order in time for a given DG TimeScheme

INPUTS:

  TimeScheme  : time scheme

OUTPUTS: 
  
  OrderTime : order in time

RETURN: Error code

*/


/******************************************************************/
// FUNCTION Prototype: xf_DGTimeSchemeVars
extern int
xf_DGTimeSchemeVars(enum xfe_TimeSchemeType TimeScheme, real *AT, 
		    real *VP, int *nq, real *tq, real *wq, real *iAT);

/*
PURPOSE:

  Returns variables associated with a TimeScheme.  Made a function for
  this because these are used multiple times (e.g. forward and adjoint
  solve), and duplicate code is error prone.

  All outputs are optional: they can be passed in as NULL if not wanted.

INPUTS:

  TimeScheme  : temporal time scheme 
  
OUTPUTS: 

  AT  : temporal stiffness matrix: Mij = int phii*dphij/dt
  iAT : inverse temporal stiffness matrix
  VP  : influence of previous U on current time slab
  nq  : number of quad points for unsteady residual eval
  tq  : quad time points on [0,1]
  wq  : quad time weights (sum to 1)
  

RETURN: Error code

*/

/******************************************************************/
// FUNCTION Prototype: xf_DGTimeInterpolate
extern int
xf_DGTimeInterpolate(enum xfe_TimeSchemeType TimeScheme, xf_Vector **Uj,
		     int nt, real *xt, xf_Vector **Vj, real *phiout );
/*

PURPOSE:

  Interpolates state on a time slab, at specified nn points.

INPUTS:

  TimeScheme  : time scheme
  Uj          : states at time nodes (i.e. time Lagrange basis coeffs)
  nt          : number of time points at which to interpolate
  xt          : locations of time points

OUTPUTS: 
  
  Vj : interpolated vectors -- must be pre-allocated
  phiout : interpolating basis at last time point (optional)

RETURN: Error code

*/

/******************************************************************/
// FUNCTION Prototype: xf_DGTimeInterpolateState
extern int
xf_DGTimeInterpolateState(xf_All *All, enum xfe_TimeSchemeType TimeScheme, xf_Vector **Uj, real *xt, xf_Vector **Vj, real *phiout );
		     
/*

PURPOSE:

  Serves as a wrapper for xf_DGTimeInterpolate and xf_DGTimeInterpolateGCL, and is used
  if the state is only being interpolated to a single time. If the GCL is active,
  the GCL vector is automatically interpolated alongside the state.

INPUTS:

  All         : All struct
  TimeScheme  : time scheme
  Uj          : states at time nodes (i.e. time Lagrange basis coeffs)
  xt          : locations of time point

OUTPUTS: 
  
  Vj : interpolated vector -- must be pre-allocated
  phiout : interpolating basis at time point (optional)

RETURN: Error code

*/
  

/******************************************************************/
//   FUNCTION Prototype: xf_PingResidual
extern int
xf_PingResidual(xf_All *All, xf_Vector *U, real ep, real tol);
/*
PURPOSE:

  Pings R_U in CalculateResidual

INPUTS:

  All: All structure
  U: state vector
  ep  : epsilon used to perturb state
  tol : if > 0, this function will exit with error if
        the ping fails by tol*ep^2

OUTPUTS:

  ping output to stdout

RETURN:

  Error code
*/





#endif // end ifndef _xf_Solver_h
