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
#include "xf_Solver.h"
#include "xf_Quad.h"
#include "xf_QuadRule.h"


// Functions for setting up and running a case within a unit test
#include "xf_UnitRun.c"


TEST_xf_ResolutionIndicator()
{
  int ierr, k, nq, nn, nn1, sr, j;
  int Order, QuadOrder;
  enum xfe_BasisType Basis;
  real *xq, *wq, *s, *s1, *s_u, *s1_u1;
  real U[6] = {1.0, 1.1, 1.2, 1.1, 1.2, 1.23};
  real U1[3], *TT;
  real Reg0[1], Reg_U0[6];
  real Reg[1] , Reg_U[6];
  real eps, tol, fd, an;
  xf_Matrix *T;
  xf_BasisData *PhiData;
  xf_BasisData *PhiData1;

  sr = 1;
  Basis = xfe_TriLagrange;
  Order = 2;

  // quad points
  QuadOrder = Order+Order-1;
  xq = wq = NULL;
  ierr = xf_Error(xf_QuadTriangle(QuadOrder, &nq, &xq, &wq));

  // restriction matrix
  T = NULL;
  ierr = xf_Error(xf_FindTransferMatrix(NULL, Basis, Order, Basis, 
					Order-1, &T));
  xf_AssertEqual(ierr, xf_OK);

  // pull off real array of the matrix
  TT = T->GenArray->rValue[0];

  // compute basis functions (and grads) if quad or basis or order changed
  PhiData  = NULL;
  ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nq, xq, 
			       xfb_Phi, &PhiData));
  xf_AssertEqual(ierr, xf_OK);
  nn = PhiData->nn;

  PhiData1 = NULL;
  ierr = xf_Error(xf_EvalBasis(Basis, Order-1, xfe_True, nq, xq, 
			       xfb_Phi, &PhiData1));
  xf_AssertEqual(ierr, xf_OK);
  nn1 = PhiData1->nn;

  // allocate memory
  ierr = xf_Error(xf_Alloc((void **) &s , nq, sizeof(real)));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_Alloc((void **) &s1, nq, sizeof(real)));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_Alloc((void **) &s_u  , nq*sr, sizeof(real)));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_Alloc((void **) &s1_u1, nq*sr, sizeof(real)));
  xf_AssertEqual(ierr, xf_OK);
  
  // restrict U to get p-1 coefficients, U1
  xf_MxM_Set(TT, U, nn1, nn, sr, U1);
  
  // interpolate state at quad points
  xf_MxM_Set(PhiData->Phi,  U,  nq, nn, sr, s);   // p
  xf_MxM_Set(PhiData1->Phi, U1, nq, nn1, sr, s1); // p-1

  // set derivatives to identity
  for (k=0; k<nq*sr; k++) s_u[k] = 1.0;
  for (k=0; k<nq*sr; k++) s1_u1[k] = 1.0;

  // Obtain resolution indicator
  ierr = xf_Error(xf_ResolutionIndicator(s, s1, wq, PhiData->Phi, PhiData1->Phi,
					 Order, nq, sr, nn, nn1, s_u, s1_u1, 
					 TT, xfe_StabSwitchLog, 0.0, 1.1, Reg0, Reg_U0));
  xf_AssertEqual(ierr, xf_OK);

  // ping test
  eps = 1e-6;
  tol = eps; // ResolutionIndicator is a higly non-linear function
  for (j=0; j<6; j++){
    U[j] += eps;

    // restrict U to get p-1 coefficients, U1
    xf_MxM_Set(TT, U, nn1, nn, sr, U1);
    
    // interpolate state at quad points
    xf_MxM_Set(PhiData->Phi,  U,  nq, nn, sr, s);   // p
    xf_MxM_Set(PhiData1->Phi, U1, nq, nn1, sr, s1); // p-1

    ierr = xf_Error(xf_ResolutionIndicator(s, s1, wq, PhiData->Phi, PhiData1->Phi,
					   Order, nq, sr, nn, nn1, s_u, s1_u1, 
					   TT, xfe_StabSwitchLog, 0.0, 1.1, Reg, Reg_U));
    xf_AssertEqual(ierr, xf_OK);
    
    fd = (Reg[0]-Reg0[0])/eps;
    an = (Reg_U[j]+Reg_U0[j])*0.5;
    if (fabs(fd-an) > tol){
      xf_printf("ping failure[%d]: fd = %.12E, an = %.12E, tol = %.6E\n",
		j, fd, an, tol);
      xf_printf("R_U0 = %.15E\n", Reg_U0[j]);
      xf_printf("R_U  = %.15E\n", Reg_U[j]);
      
      return xf_Error(xf_PING_FAILED);
    }
    U[j] -= eps;
  }


  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
  
  /* Destroy transfer matrix */
  ierr = xf_Error(xf_DestroyMatrix(T));
  xf_AssertEqual(ierr, xf_OK);

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  xf_Release( (void *) s);
  xf_Release( (void *) s1);
  xf_Release( (void *) s_u);
  xf_Release( (void *) s1_u1);
  

  return xf_OK;
}



TEST_xf_CalculateRegStab_Euler()
{
  /*
    Tests regularity stabilization calculation and linearization
  */

  int ierr, i, irun, nrun;
  int r, n, egrp, elem;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  char *StabKeyValue[] = {"Stabilization", "Resolution", "\0"};
  char *StabSwitch[] = {"Const", "Linear", "Square"};
  real epv[]  = {1e-5 , 1e-5, 1e-5};
  real tolv[] = {1e-10, 1e-9, 1e-9};
  real fd, an, tol, ep;
  real *E, *E_U;
  xf_Vector *U, *V;
  xf_Vector *Reg0, *Reg_U0;
  xf_SolverData *SolverData;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAllEuler(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // add stabilization term
  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // set U to a perturbed value
  ierr = xf_Error(xf_DuplicateVector(Mesh, U, &V));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorRand(V, 17));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorMultSet(V, 1e-2, xfe_Add, U));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_DestroyVector(V, xfe_True));
  xf_AssertEqual(ierr, xf_OK);


  // create SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  xf_AssertEqual(ierr, xf_OK);

  // set stab-switch constant value for use below
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitchValue", "1.0"));
  if (ierr != xf_OK) return ierr;

  nrun = 3;
  for (irun=0; irun<nrun; irun++){

    // each run tests a different stabilization switch
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", StabSwitch[irun]));
    if (ierr != xf_OK) return ierr;
    
    // likewise, a different ping tolerance is associated with each run;
    ep  =  epv[irun];
    tol = tolv[irun];
    

    // calculate jump stabilization
    ierr = xf_Error(xf_CalculateStabilization(All, U, xfe_True, SolverData));
    xf_AssertEqual(ierr, xf_OK);
    
    // store initial vectors 
    ierr = xf_Error(xf_DuplicateVector(Mesh, SolverData->StabData.Reg, &Reg0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DuplicateVector(Mesh, SolverData->StabData.Reg_U, &Reg_U0));
    if (ierr != xf_OK) return ierr;    
    
    r = U->GenArray[0].r;
    
    // Ping test
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	for (n=0; n<r; n++){
	  U->GenArray[egrp].rValue[elem][n] += ep;
	  ierr = xf_Error(xf_CalculateStabilization(All, U, xfe_True, SolverData));
	  xf_AssertEqual(ierr, xf_OK);
	  E   = SolverData->StabData.Reg->GenArray[egrp].rValue[elem];
	  E_U = SolverData->StabData.Reg_U->GenArray[egrp].rValue[elem];
	  
	  fd = (E[0]   -   Reg0->GenArray[egrp].rValue[elem][0])/ep;
	  an = (E_U[n] + Reg_U0->GenArray[egrp].rValue[elem][n])*.5;
	  
	  if (fabs(fd-an) > tol){
	    xf_printf("ping failure[egrp=%d, elem=%d, n=%d]: fd = %.12E, an = %.12E, tol = %.6E\n",
		      egrp, elem, n, fd, an, tol);
	  }
	  xf_AssertWithin(fabs(fd-an), 0.0, tol);
	  U->GenArray[egrp].rValue[elem][n] -= ep;
	}
      } // ibface
    } // ibfgrp
    
    // Destroy vectors
    ierr = xf_Error(xf_DestroyVector(Reg0, xfe_True));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyVector(Reg_U0, xfe_True));
    if (ierr != xf_OK) return ierr;

  }

  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  xf_AssertEqual(ierr, xf_OK);


  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_CalculateJumpStab_ScalarDiff()
{
  /*
    Tests jump stabilization calculation and linearization
  */

  int ierr, i, irun, nrun;
  int ibfgrp, ibface, nibface, n;
  int r, egrpL, elemL, egrpR, elemR;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  char *ConvKeyValue[] = {"VelocityFcn", "Constant", "VelocityData", "1.0 0.0", "\0"};
  char *StabKeyValue[] = {"Stabilization", "Jump", "\0"};
  char *StabSwitch[] = {"Const", "Linear", "Square"};
  real epv[]  = {1e-5 , 1e-5, 1e-5};
  real tolv[] = {1e-10, 1e-2, 1e-6};
  real fd, an, tol, ep;
  real *J, *J_UL;
  xf_Vector *U;
  xf_Vector *Jump0, *Jump_UL0, *Jump_UR0;
  xf_SolverData *SolverData;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // add regularity scalar key
  ierr = xf_AddKeyValue(&All->EqnSet->KeyValue, "RegularityScalar", "Scalar", xfe_True);
  xf_AssertEqual( ((ierr == xf_OK) || (ierr == xf_OVERWROTE)), 1);

  // Inactivate existing residual terms, add convection and stabilization terms
  for (i=0; i<All->EqnSet->ResTerms->nResTerm; i++)
    All->EqnSet->ResTerms->ResTerm[i].Active = xfe_False;

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermConv, ConvKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // set U to a perturbed value
  ierr = xf_Error(xf_VectorRand(U, 17));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorMult(U, 1e-3));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorAdd(U, 1.0));
  xf_AssertEqual(ierr, xf_OK);

  
  // create SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  xf_AssertEqual(ierr, xf_OK);

  // set stab-switch constant value for use below
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitchValue", "1.0"));
  if (ierr != xf_OK) return ierr;

  nrun = 3;
  for (irun=0; irun<nrun; irun++){

    // each run tests a different stabilization switch
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", StabSwitch[irun]));
    if (ierr != xf_OK) return ierr;

    // likewise, a different ping tolerance is associated with each run;
    ep  =  epv[irun];
    tol = tolv[irun];

    // calculate jump stabilization
    ierr = xf_Error(xf_CalculateStabilization(All, U, xfe_True, SolverData));
    xf_AssertEqual(ierr, xf_OK);
    
    // store initial vectors
    ierr = xf_Error(xf_DuplicateVector(Mesh, SolverData->StabData.Jump, &Jump0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DuplicateVector(Mesh, SolverData->StabData.Jump_UL, &Jump_UL0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DuplicateVector(Mesh, SolverData->StabData.Jump_UR, &Jump_UR0));
    if (ierr != xf_OK) return ierr;

    r = U->GenArray[0].r;

    // Ping test

    // Loop over interior group and boundary face groups
    for (ibfgrp=-1; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
      nibface = ((ibfgrp == -1) ? Mesh->nIFace : Mesh->BFaceGroup[ibfgrp].nBFace);
      // loop over faces
      for (ibface=0; ibface<nibface; ibface++){
	xf_FaceElements(Mesh, ibfgrp, ibface, &egrpL, &elemL, NULL, &egrpR, &elemR, NULL);

	// ping L
	for (n=0; n<r; n++){
	  U->GenArray[egrpL].rValue[elemL][n] += ep;
	  ierr = xf_Error(xf_CalculateStabilization(All, U, xfe_True, SolverData));
	  xf_AssertEqual(ierr, xf_OK);
	  J    = SolverData->StabData.Jump->GenArray[1+ibfgrp].rValue[ibface];
	  J_UL = SolverData->StabData.Jump_UL->GenArray[1+ibfgrp].rValue[ibface];
	
	  fd = (J[0]    - Jump0->GenArray[1+ibfgrp].rValue[ibface][0])/ep;
	  an = (J_UL[n] + Jump_UL0->GenArray[1+ibfgrp].rValue[ibface][n])*.5;

	  if (fabs(fd-an) > tol){
	    xf_printf("ping failure[ibfgrp=%d, ibface=%d, n=%d]: fd = %.12E, an = %.12E, tol = %.6E\n",
		      ibfgrp, ibface, n, fd, an, tol);
	  }
	  xf_AssertWithin(fabs(fd-an), 0.0, tol);
	  U->GenArray[egrpL].rValue[elemL][n] -= ep;
	}
      } // ibface
    } // ibfgrp

  
    // Destroy vectors
    ierr = xf_Error(xf_DestroyVector(Jump0, xfe_True));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyVector(Jump_UL0, xfe_True));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyVector(Jump_UR0, xfe_True));
    if (ierr != xf_OK) return ierr;

  } // irun

  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  xf_AssertEqual(ierr, xf_OK);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;
}



TEST_xf_CalculateJumpStab_Euler()
{
  /*
    Tests jump stabilization calculation and linearization
  */

  int ierr, i, irun, nrun;
  int ibfgrp, ibface, nibface, n;
  int r, egrpL, elemL, egrpR, elemR;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 2;
  char *StabKeyValue[] = {"Stabilization", "Jump", "\0"};
  char *StabSwitch[] = {"Const", "Linear", "Square"};
  real epv[]  = {1e-5 , 1e-5, 1e-5};
  real tolv[] = {1e-10, 1e-4, 5e-6};
  real fd, an, tol, ep;
  real *J, *J_UL;
  xf_Vector *U, *V;
  xf_Vector *Jump0, *Jump_UL0, *Jump_UR0;
  xf_SolverData *SolverData;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAllEuler(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);
  
  Mesh = All->Mesh;

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  // add stabilization term
  ierr = xf_Error(xf_AddResTerm(All->EqnSet->ResTerms, xfe_ResTermDiff, StabKeyValue));
  xf_AssertEqual(ierr, xf_OK);

  // set U to a perturbed value
  ierr = xf_Error(xf_DuplicateVector(Mesh, U, &V));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorRand(V, 17));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_VectorMultSet(V, 1e-2, xfe_Add, U));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_DestroyVector(V, xfe_True));
  xf_AssertEqual(ierr, xf_OK);


  // create SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  xf_AssertEqual(ierr, xf_OK);

  // set stab-switch constant value for use below
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitchValue", "1.0"));
  if (ierr != xf_OK) return ierr;

  nrun = 3;
  for (irun=0; irun<nrun; irun++){

    // each run tests a different stabilization switch
    ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "StabSwitch", StabSwitch[irun]));
    if (ierr != xf_OK) return ierr;
    
    // likewise, a different ping tolerance is associated with each run;
    ep  =  epv[irun];
    tol = tolv[irun];
    

    // calculate jump stabilization
    ierr = xf_Error(xf_CalculateStabilization(All, U, xfe_True, SolverData));
    xf_AssertEqual(ierr, xf_OK);
    
    // store initial vectors 
    ierr = xf_Error(xf_DuplicateVector(Mesh, SolverData->StabData.Jump, &Jump0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DuplicateVector(Mesh, SolverData->StabData.Jump_UL, &Jump_UL0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DuplicateVector(Mesh, SolverData->StabData.Jump_UR, &Jump_UR0));
    if (ierr != xf_OK) return ierr;
    
    
    r = U->GenArray[0].r;
    
    // Ping test
    // Loop over interior group and boundary face groups
    for (ibfgrp=-1; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
      nibface = ((ibfgrp == -1) ? Mesh->nIFace : Mesh->BFaceGroup[ibfgrp].nBFace);
      // loop over faces
      for (ibface=0; ibface<nibface; ibface++){
	xf_FaceElements(Mesh, ibfgrp, ibface, &egrpL, &elemL, NULL, &egrpR, &elemR, NULL);
	
	// ping L
	for (n=0; n<r; n++){
	  U->GenArray[egrpL].rValue[elemL][n] += ep;
	  ierr = xf_Error(xf_CalculateStabilization(All, U, xfe_True, SolverData));
	  xf_AssertEqual(ierr, xf_OK);
	  J    = SolverData->StabData.Jump->GenArray[1+ibfgrp].rValue[ibface];
	  J_UL = SolverData->StabData.Jump_UL->GenArray[1+ibfgrp].rValue[ibface];
	  
	  fd = (J[0]    -    Jump0->GenArray[1+ibfgrp].rValue[ibface][0])/ep;
	  an = (J_UL[n] + Jump_UL0->GenArray[1+ibfgrp].rValue[ibface][n])*.5;
	  
	  if (fabs(fd-an) > tol){
	    xf_printf("ping failure[ibfgrp=%d, ibface=%d, n=%d]: fd = %.12E, an = %.12E, tol = %.6E\n",
		      ibfgrp, ibface, n, fd, an, tol);
	  }
	  xf_AssertWithin(fabs(fd-an), 0.0, tol);
	  U->GenArray[egrpL].rValue[elemL][n] -= ep;
	}
      } // ibface
    } // ibfgrp
    
    // Destroy vectors
    ierr = xf_Error(xf_DestroyVector(Jump0, xfe_True));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyVector(Jump_UL0, xfe_True));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyVector(Jump_UR0, xfe_True));
    if (ierr != xf_OK) return ierr;

  }

  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  xf_AssertEqual(ierr, xf_OK);



  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}
