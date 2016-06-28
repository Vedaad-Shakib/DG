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


TEST_xf_InitState()
{
  /*
    This tests state initialization via Interpolation and QR
    projection
  */

  int ierr, sr, nq, iq, k, itest, iBasis;
  enum xfe_BasisType BasisVec[2] = {xfe_TriLagrange, xfe_TriHierarch};
  enum xfe_BasisType Basis = xfe_TriLagrange;
  enum xfe_Bool QuadChanged;
  int Order, OrderStart=2, OrderEnd=5;
  int egrp, elem, dim;
  int *IParam;
  char *Type0, *Function0, *Data0;
  char Type1[] = "FullState", *Function1 = NULL, Data1[] = "0.5";
  char Type2[] = "FullState", Function2[] = "Test", Data2[] = "0.0";
  real *u, *ue, *xglob, *xq, *EU;
  real *RParam;
  real Time=0.;
  xf_Vector *U;
  xf_ICs *ICs;
  xf_IC  *IC;
  xf_BasisData *PhiData, *GeomPhiData;
  xf_QuadData *QuadData;
  xf_All *All;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;

  Order = OrderStart;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  Mesh = All->Mesh;
  EqnSet = All->EqnSet;
  dim = Mesh->Dim;

  ICs = EqnSet->ICs;
  xf_AssertEqual((ICs->nIC>0), 1);
  
  IC  = ICs->IC+0;

  // store original values
  Type0     = IC->Type;  
  Function0 = IC->Function;
  Data0     = IC->Data;

  // build eqnset-desired parameter lists for passing into EqnSetICState
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  xf_AssertEqual(ierr, xf_OK);

  PhiData     = NULL;
  GeomPhiData = NULL;
  QuadData    = NULL;
  u           = NULL;
  ue          = NULL;
  xglob       = NULL;
  sr          = U->StateRank;


  for (itest=0; itest<4; itest++){

    if ((itest%2) == 0){
      // first test
      IC->Type     = Type1     ;
      IC->Function = Function1 ; 
      IC->Data     = Data1     ;
    }
    else{
      // second test
      IC->Type     = Type2     ;
      IC->Function = Function2 ; 
      IC->Data     = Data2     ;
    }

    if (itest > 1){ // tests are repeated for interpolating ICs instead of doing QR
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "InterpolateIC", "True"));
      if (ierr != xf_OK) return ierr;
    }
    else{
      ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "InterpolateIC", "False"));
      if (ierr != xf_OK) return ierr;
    }

    for (iBasis=0; iBasis<2; iBasis++){  // testing several bases

      Basis = BasisVec[iBasis];

      for (Order = OrderStart; Order<=OrderEnd; Order++){ // testing several orders
	
	// project U to desired Basis and Order
	ierr = xf_Error(xf_ProjectVectorInPlace_OrderSet(All->Mesh, NULL, U, NULL, Basis, NULL, 
							 NULL, Order));
	xf_AssertEqual(ierr, xf_OK);
      
	// initialize U
	ierr = xf_Error(xf_InitState(All, U));
	xf_AssertEqual(ierr, xf_OK);

	QuadChanged = xfe_True;

	// Check if u is exact value on entire domain
	for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

	  for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	    // max quad points on elem
	    ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, -1, &QuadData, &QuadChanged));
	    xf_AssertEqual(ierr, xf_OK);
	    nq = QuadData->nquad;
	    xq = QuadData->xquad;

	    // compute basis functions
	    ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged,  nq, xq, xfb_Phi, &PhiData));
	    xf_AssertEqual(ierr, xf_OK);

	    // reallocate memory
	    ierr = xf_Error(xf_ReAlloc( (void **) &u, nq*sr, sizeof(real)));
	    xf_AssertEqual(ierr, xf_OK);

	    ierr = xf_Error(xf_ReAlloc( (void **) &ue, nq*sr, sizeof(real)));
	    xf_AssertEqual(ierr, xf_OK);

	    ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
	    xf_AssertEqual(ierr, xf_OK);

	    // Calculate U at xq -> u
	    EU = U->GenArray[egrp].rValue[elem];
	    xf_MxM_Set(PhiData->Phi, EU, nq, PhiData->nn, sr, u);

	    // obtain global coords of quad points
	    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
					    nq, xq, xglob));
	    xf_AssertEqual(ierr, xf_OK);

	    // determine exact values of u at xq
	    ierr = xf_Error(xf_EqnSetICState(EqnSet, IC, IParam, RParam, nq, xglob, &Time, ue));
	    xf_AssertEqual(ierr, xf_OK);

	    // check if got IC exactly
	    for (iq=0; iq<nq; iq++)
	      xf_AssertWithin(u[iq], ue[iq], UTOL5);
      
	  } // elem
    
	} // egrp
      } // Order
    }// iBasis
  }// itest

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
  
  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  xf_Release((void *) u);
  xf_Release((void *) ue);
  xf_Release((void *) xglob);

  
  xf_Release((void *) IParam);
  xf_Release((void *) RParam);

  // reset original
  IC->Type     = Type0     ;
  IC->Function = Function0 ; 
  IC->Data     = Data0     ;


  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}



TEST_xf_ChangeVariableSet()
{
  /*
    This tests the changing of variable set for a state vector.
  */

  int ierr, i, j, k;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order = 4;
  xf_Vector *U, *V;
  xf_GenArray *gU, *gV;
  xf_All *All;
  xf_Mesh *Mesh;

  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  Mesh = All->Mesh;

  // Copy Vector U to V
  ierr = xf_Error(xf_CreateVector(&V));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_CopyVector(Mesh, U, V));
  xf_AssertEqual(ierr, xf_OK);  

 
  // Change variable set in V using the "Test" variable set
  ierr = xf_Error(xf_ChangeVariableSet(All, "Test", V));
  xf_AssertEqual(ierr, xf_OK);

  // Test that V = - U
  for (i=0; i<U->nArray; i++){
    gU = U->GenArray+i;
    gV = V->GenArray+i;
    for (j=0; j<gU->n; j++){
      for (k=0; k<gU->r; k++)
	xf_AssertWithin(gU->rValue[j][k], -gV->rValue[j][k], UTOL4);
    } // j
  } // i
  
  // destroy V
  ierr = xf_Error(xf_DestroyVector(V, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}


TEST_xf_ReconstructVector()
{
  /*
    This tests state reconstruction.
  */

  int ierr;
  enum xfe_BasisType Basis = xfe_TriLagrange;
  int Order=2;
  char *Type0, *Function0, *Data0;
  char Type1[] = "FullState", Function1[] = "Test", Data1[] = "0.0";
  real norm;
  xf_Vector *U, *V;
  xf_ICs *ICs;
  xf_IC  *IC;
  xf_All *All;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  // Includes Param and EqnSet
  ierr = xf_Error(xf_UnitBoxQ1TriangleAll(&All, xfe_False));
  xf_AssertEqual(ierr, xf_OK);

  // Load dynamic library, register eqnset, initialize solution
  ierr = xf_Error(xf_InitializeTestRun(All, Basis, Order, &U));
  xf_AssertEqual(ierr, xf_OK);

  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;
  ICs    = EqnSet->ICs;
  xf_AssertEqual((ICs->nIC>0), 1);  
  IC     = ICs->IC+0;

  // store original IC values
  Type0     = IC->Type;  
  Function0 = IC->Function;
  Data0     = IC->Data;

  // set function to Test: should be a quadratic
  IC->Type     = Type1     ;
  IC->Function = Function1 ; 
  IC->Data     = Data1     ;

  // reconstruct into a cubic -> V
  ierr = xf_Error(xf_ReconstructVector(All, 1, U, xfe_True, &V));
  if (ierr != xf_OK) return ierr;

  // project U to cubic (injection)
  ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, U, NULL, 
							 xfe_BasisLast, 1));
  if (ierr != xf_OK) return ierr;

  // subtract U from V
  ierr = xf_Error(xf_SetVector(U, xfe_Sub, V));
  if (ierr != xf_OK) return ierr;

  // check norm of V (should be zero)
  ierr = xf_Error(xf_VectorNorm(V, 1, &norm));
  if (ierr != xf_OK) return ierr;
  xf_AssertWithin(norm, 0., UTOL5);

  // reset original
  IC->Type     = Type0     ;
  IC->Function = Function0 ; 
  IC->Data     = Data0     ;

  // Destroy All
  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, xf_OK);

  return xf_OK;  
}




