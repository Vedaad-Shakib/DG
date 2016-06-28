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


#include "xf_Unit.h"
#include "xf_All.h"
#include "xf_Data.h"


// Functions for setting up and running a case within a unit test
#include "xf_UnitRun.c"

TEST_xf_ScalarDiff_ErrEstOutput()
{
  /*
    This tests output error estimation for a steady scalar diffusion
    problem: heat conduction from a hot temperature on the left to a
    sink on the right of a square domain.  A linear source term is
    also included in the domain so that the solution is not p=1.  The
    computed integrated heat flux, "HeatFlow", is the output of
    interest.  The error estimate at p=1 for this output is compared
    to the actual error between a p=2 and a p=2 solution.
  */

  int ierr, nPsi;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  char *ConvKeyValue[] = {"VelocityFcn", "Constant", "VelocityData", ".1 0.0", "\0"};
  char *SourceKeyValue[] = {"SourceFcn", "Linear", "SourceData", "0.1", "\0"};
  char *ErrEstKeyValue[] = { "ErrEstOrderIncrement", "1",
			     "FineSpace_nIterNonlinear", "10",
			     "FineSpace_LinearSolver", "GMRES",
			     "FineSpace_nIterLinear", "20",
			     "FineSpace_Preconditioner", "BlockJacobi",
			     "FineSpace_nIterAdjoint", "20",
			     "\0"};
  real JH, Jh, Jtrue, OutputError, OutputErrorTest;
  xf_Vector *U, **Psi;
  xf_All *All;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);

  // Set custom error estimation key values
  ierr = xf_Error(xf_AddKeyValueList(&All->Param->KeyValue, ErrEstKeyValue, xfe_True, xfe_False));
  if (ierr != xf_OK) return ierr;

  // Inactivate diffusion, add convection and source terms
  All->EqnSet->ResTerms->ResTerm[0].Active = xfe_False;

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermConv, ConvKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermSource, SourceKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);


  // Solve system
  ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Calculate heat flux out of right side
  ierr = xf_Error(xf_CalculateOutput(All, "HeatFlow", U, &JH, NULL, xfe_Set));
  xf_AssertEqual(ierr, xf_OK);

  // locate adjoint (set to zero)
  ierr = xf_Error(xf_FindAdjointVectors(All, U, "HeatFlow", xfe_False, xfe_True, 
					&nPsi, &Psi, NULL));
  xf_AssertEqual(ierr, xf_OK);

  // solve for adjoint
  ierr = xf_Error(xf_SolveAdjoints(All, 0.0, 1.0, xfe_False, U, nPsi, NULL, Psi, NULL, xfe_False, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
 
  xf_Release( (void *) Psi);

  // estimate the error
  ierr = xf_Error(xf_ErrEstOutput(All, "HeatFlow", "None", xfe_False, xfe_False, NULL, xfe_False, &OutputError));
  if (ierr != xf_OK) return ierr;

  // Galerkin orthogonality check by not doing any iterations on the adjoint problem
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "FineSpace_nIterAdjoint", "0"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ErrEstOutput(All, "HeatFlow", "None", xfe_False, xfe_False, NULL, xfe_False, &OutputErrorTest));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(OutputErrorTest, 0.0, UTOL3);

 
  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);
  

  /** Re-do solve with Order+1 **/

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order+1, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Inactivate diffusion, add convection and source terms
  All->EqnSet->ResTerms->ResTerm[0].Active = xfe_False;

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermConv, ConvKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermSource, SourceKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // Solve system
  ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &U));
  xf_AssertEqual(ierr, xf_OK);

  // Calculate heat flux out of right side
  ierr = xf_Error(xf_CalculateOutput(All, "HeatFlow", U, &Jh, NULL, xfe_Set));
  xf_AssertEqual(ierr, xf_OK);

  // Check if calculated output error correctly
  xf_AssertWithin(JH-Jh, OutputError, UTOL3);

  
  // Also, for kicks, test if Jh is close to the analytical solution
  //Jtrue = -2*0.1*(-sin(2.) + (-cos(2.)/sin(2.))*cos(2.)); // diffusion-reactuib
  Jtrue = 2*.1*exp(-2); // convection-reaction
  xf_AssertWithin(Jh, Jtrue, 1e-7); // for p=3 this should be enough


  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}
