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

/*
 FILE:  xf_SolverTools.c
 
 This file contains helper functions for the solver.
 
 */

#include "xf_AllStruct.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Residual.h"
#include "xf_Line.h"
#include "xf_Basis.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Memory.h"
#include "xf_Quad.h"
#include "xf_MeshTools.h"
#include "xf_MathLapack.h"
#include "xf_EqnSetHook.h"
#include "xf_Solver.h"
#include "xf_LinearSolver.h"
#include "xf_MeshMotionGCL.h"
#include "xf_Output.h"
#include "xfYu_Model.h"

/******************************************************************/
//   FUNCTION Definition:  xf_CheckUserHalt
enum xfe_Bool
xf_CheckUserHalt(char *pfname)
{
  int ierr, myRank;
  enum xfe_Bool ReturnFlag = xfe_False;
  char fname0[]="STOP", *fname;
  char line[xf_MAXLINELEN];
  FILE *fid;
  
  fname = ((pfname == NULL) ? fname0 : pfname);
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0){
    if ( (fid = fopen(fname,"r")) != NULL ){
      if (fgets(line, xf_MAXLINELEN, fid) != NULL){
        if ( strncmp(line, "STOP",4) == 0 ) {
          xf_printf("User halt.\n");
          ReturnFlag = xfe_True;
        }
      }
      fclose(fid);
    }
  }
  
  // reduce-max ReturnFlag
  ierr = xf_Error(xf_MPI_Allreduce(&ReturnFlag, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  
  return ReturnFlag;
}


/******************************************************************/
//   FUNCTION Definition: xf_MultMassMatrix
int
xf_MultMassMatrix(xf_All *All, real c, xf_Vector *U)
{
  //  U = c*M*U
  int ierr, egrp, elem;
  int pOrder, Order, i, nn, sr;
  enum xfe_BasisType Basis;
  real *MM, fac, val, *EU, *y = NULL;
  xf_Matrix *M;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  sr = U->StateRank;
  
  if (c == 0.0) return xf_OK; // nothing to do
  
  

  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    pOrder = -1;
    Basis = U->Basis[egrp];

    //  loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);

      if (Order != pOrder){
	pOrder = Order;
	//  find mass matrix data (if generic, n==1)
	ierr = xf_Error(xf_FindMassMatrixData(All, egrp, Basis, Order, &M));
	if (ierr != xf_OK) return ierr;
      
	ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
	if (ierr != xf_OK) return ierr;
      
	ierr = xf_Error(xf_ReAlloc((void **) &y, nn*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      ierr = xf_Error(xf_ElemMassMatrix(All, egrp, elem, Basis, Order, M, NULL, &MM, &fac));
      if (ierr != xf_OK) return ierr;
      
      // y = c*M*U
      EU = U->GenArray[egrp].rValue[elem];
      for (i=0; i<nn*sr; i++) y[i] = 0.;
      xf_cMxM_Add(c*fac, MM, EU, nn, nn, sr, y);
      
      for (i=0; i<nn*sr; i++) EU[i] = y[i];
    } // elem
  } // egrp
  
  xf_Release((void *) y);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MultInvMassMatrix
int
xf_MultInvMassMatrix(xf_All *All, real c0, xf_Vector *dt, xf_Vector *U)
{
  //  U = c*M^{-1}*U
  int ierr, egrp, elem;
  int Order, pOrder, i, nn, sr;
  enum xfe_BasisType Basis;
  real *iMM, fac, val, *EU, *y = NULL;
  real c;
  xf_Matrix *iM;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  sr = U->StateRank;
  
  if (c0 == 0.0){
    ierr = xf_Error(xf_SetZeroVector(U));
    if (ierr != xf_OK) return ierr;
    return xf_OK; // nothing else to do
  }
  
  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    pOrder = -1;
    Basis = U->Basis[egrp];

    //  loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);
      
      if (pOrder != Order){
	pOrder = Order;

	//  find inverse mass matrix data (if generic, n==1)
	ierr = xf_Error(xf_FindInvMassMatrixData(All, egrp, Basis, Order, &iM));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_ReAlloc((void **) &y, nn*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      ierr = xf_Error(xf_ElemInvMassMatrix(All, egrp, elem, Basis, Order, iM, NULL, &iMM, &fac));
      if (ierr != xf_OK) return ierr;
      
      // include dt if not NULL
      c = (dt == NULL) ? c0 : c0*dt->GenArray[egrp].rValue[elem][0];
      
      // y = c*iM*U
      EU = U->GenArray[egrp].rValue[elem];
      for (i=0; i<nn*sr; i++) y[i] = 0.;
      xf_cMxM_Add(c*fac, iMM, EU, nn, nn, sr, y);
      for (i=0; i<nn*sr; i++) EU[i] = y[i];
    } // elem
  } // egrp
  
  xf_Release((void *) y);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MultStiffMatrix
int
xf_MultStiffMatrix(xf_All *All, real c, xf_Vector *U)
{
  //  U = c*K*U, K_ij = int(gPhi{i}*gphi{j})
  int ierr, k, dim, egrp, elem, pOrder, Order, QuadOrder;
  int nq, pnq, nn, sr, r, *IParam;
  int i, iq, d;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged;
  real *xq, *EU, *u, *gu, *y, *RParam;
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_QuadData *QuadData;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  sr   = All->EqnSet->StateRank;
  
  if (c == 0.0) return xf_OK; // nothing to do
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  QuadData    = NULL;
  PhiData     = NULL;
  JData       = NULL;
  u           = NULL;
  gu          = NULL;
  y           = NULL;
  pnq         = -1;
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    pOrder = -1;
    Basis = U->Basis[egrp];

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);

      if (Order != pOrder){
	pOrder = Order;
	// determine quadrature order
	ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*Order, &QuadOrder));
	if (ierr != xf_OK) return ierr;
      
	ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_ReAlloc((void **) &y, nn*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      
      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);
      
      if (Order != pOrder){
        pOrder = Order;
        // determine quadrature order
        ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*Order, &QuadOrder));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ReAlloc((void **) &y, nn*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
      }
      
      /* Pull off quad points for the element; will not recalculate if
       Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions (and grads) if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, 
                                   xfb_Phi | xfb_GPhi | xfb_gPhi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      /* Compute geometry Jacobian; if not constant, compute at quad
       points.  Note if jacobian is constant, only one Jacobian will
       be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ | xfb_iJ, 
                                      QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;
      
      // convert reference basis grads (GPhi) to physical grads, gPhi
      ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
      if (ierr != xf_OK) return ierr;
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **) &gu, nq*sr*dim, sizeof(real)));
        if (ierr != xf_OK) return ierr;
      }
      
      nn = PhiData->nn;
      EU = U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
      
      // interpolate gradient
      for (d=0; d<dim; d++)
        xf_MxM_Set(PhiData->gPhi+nn*nq*d, EU, nq, nn, sr, gu+nq*sr*d);
      
      
      // multiply gu by quad weigths and detJ
      for (iq=0; iq<nq; iq++)
        for (d=0; d<dim; d++)
          gu[nq*sr*d+iq] *= QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      
      // set y = sum_i sum_q gPhi{i,q,n}^T*gu{i,q,k}
      xf_MTxM_Set(PhiData->gPhi+nn*nq*0, gu+nq*sr*0, nn, nq, sr, y);
      for (d=1; d<dim; d++)
        xf_MTxM_Add(PhiData->gPhi+nn*nq*d, gu+nq*sr*d, nn, nq, sr, y);
      
      // set EU = c*y
      for (i=0; i<nn*sr; i++) EU[i] = c*y[i];
      
      pnq = nq;
    } // elem
  } // egrp
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  // Release memory
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) u);
  xf_Release( (void *) y);
  xf_Release( (void *) gu);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AddMassMatrix
int
xf_AddMassMatrix(xf_All *All, real c, xf_Vector *dt, xf_Vector *U, 
                 xf_Vector *R, xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
  //  R += c*M*U, R_U += (c+1/dt)*M
  int ierr, egrp, elem;
  int pOrderU, pOrderR, OrderU, OrderR, i, nnU, nnR, k, sr, sr2;
  enum xfe_BasisType BasisU, BasisR;
  enum xfe_Bool PenalizeResidual;
  real *MM, fac, val, *A, *EU, *ER, P;
  xf_Matrix *M;
  xf_Mesh *Mesh;
  xf_Vector *Pvec;
  
  // return immediately if nothing to do
  if ((c == 0.) && ((R_U == NULL) || (R_U->Value == NULL) || (dt == NULL))) return xf_OK;
  
  if (SolverData != NULL)
    PenalizeResidual = SolverData->PenalizeResidual;
  else 
    PenalizeResidual = xfe_False;
  
  if (PenalizeResidual) {
    Pvec = SolverData->Pvec;
  }
  
  Mesh = All->Mesh;
  sr = U->StateRank;
  sr2 = sr*sr;

  // consistency check in state rank
  if ((R != NULL) && (sr != R->StateRank)) return xf_Error(xf_INPUT_ERROR);
  
  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    BasisU = U->Basis[egrp];
    BasisR = ((R==NULL) ? BasisU : R->Basis[egrp]);
    pOrderU = -1;
    pOrderR = -1;
    
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      //get penalization factor
      if (PenalizeResidual)
        P = Pvec->GenArray[egrp].rValue[elem][0];
      else 
        P = 0.0;
      
      // get interpolation order
      OrderU = xf_InterpOrder(U, egrp, elem);
      OrderR = ((R==NULL) ? OrderU : xf_InterpOrder(R, egrp, elem));
      
      if ((pOrderU != OrderU) || (pOrderR != OrderR)){
        pOrderU = OrderU;
        pOrderR = OrderR;
        
        // find general mass matrix data (if generic, n==1)
        ierr = xf_Error(xf_FindGenMassMatrixData(All, egrp, BasisR, OrderR, BasisU, OrderU, &M));
        if (ierr != xf_OK) return ierr;
        
        // nnR
        ierr = xf_Error(xf_Order2nNode(BasisR, OrderR, &nnR));
        if (ierr != xf_OK) return ierr;
        
        // nnU
        ierr = xf_Error(xf_Order2nNode(BasisU, OrderU, &nnU));
        if (ierr != xf_OK) return ierr;
        
      }
      
      ierr = xf_Error(xf_ElemGenMassMatrix(All, egrp, elem, BasisR, OrderR, BasisU, OrderU, 
                                           M, NULL, &MM, &fac));
      if (ierr != xf_OK) return ierr;
      
      // R += c*M*U
      if ((c != 0.) && (R != NULL)){
        EU = U->GenArray[egrp].rValue[elem];
        ER = R->GenArray[egrp].rValue[elem];
        xf_cMxM_Add(c*fac, MM, EU, nnR, nnU, sr, ER);
      }
      
      if (dt == NULL)
        fac *= c;
      else
        fac *= (c + 1.0/dt->GenArray[egrp].rValue[elem][0]);
      
      // R_U += M*(c + 1/dt)
      if (R_U != NULL){
        A = R_U->Value[egrp][elem][0];
        for (i=0; i<nnR*nnU; i++){
          val = MM[i]*fac/(1.0+P);
          for (k=0; k<sr2; k+=(sr+1))
            A[i*sr2+k] += val;
        } // i
      }
    } // elem
  } // egrp
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_InitAlterState
int 
xf_InitAlterState(xf_All *All,  const char *FcnName, 
                  const char *FcnParam, xf_Vector *U)
{
  int ierr, i, j, k, dim, sr, nn, egrp, elem;
  int iq, nq, pnq, pOrder, Order, QuadOrder;
  int nFcnParam;
  int *IParam, *P;
  enum xfe_Bool QuadChanged, InterpolateIC;
  enum xfe_Bool Altering = xfe_False;
  enum xfe_BasisType Basis;
  enum xfe_BasisType *BasisNew;
  real *u, *xq, *xglob, x0[3] = {0.};
  real *xn, *EU, *RParam;
  real *Q, *R;
  real *FParam = NULL;
  real Time;
  xf_GenArray *ga;
  xf_BasisData *PhiData, *GeomPhiData;
  xf_QuadData *QuadData;
  xf_ICs *ICs;
  xf_IC  *IC;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  
  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;
  dim    = Mesh->Dim;
  sr     = EqnSet->StateRank;
  
  if (U->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  
  if (U->StateRank != sr) return xf_Error(xf_OUT_OF_BOUNDS);
  
  Altering = ((FcnName != NULL) && (xf_NotNull(FcnName)));
  
  // pull off parameters
  if (Altering && (FcnParam != NULL)){
    ierr = xf_Error(xf_ScanXRealAlloc(FcnParam, &nFcnParam, &FParam));
    if (ierr != xf_OK) return ierr;
  }
  
  // set U->StateName
  ierr = xf_Error(xf_ReAlloc2((void ***)&U->StateName, sr, xf_MAXSTRLEN, sizeof(char)));
  if (ierr != xf_OK) return ierr;
  
  for (i=0; i<sr; i++) strcpy(U->StateName[i], EqnSet->StateName[i]);
  
  // Pull off first initial condition
  if (!Altering){
    ICs = EqnSet->ICs;
    
    if (ICs->nIC <= 0){
      xf_printf("Error, no initial conditions present!\n");
      return xf_Error(xf_OUT_OF_BOUNDS);
    }
    IC  = ICs->IC+0;
  }
  
  // build eqnset-desired parameter lists for passing into EqnSetICState
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // determine simulation time
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time);
  if (ierr != xf_OK) return ierr;

  // Determine whether should interpolate IC if a function
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "InterpolateIC", &InterpolateIC);
  if (ierr != xf_OK) return ierr;
  
  
  if ((!Altering) && (IC->Function == NULL)){
    
    /*--------------------------*/
    /* No IC Function specified */
    /*--------------------------*/
    
    /* Only need to call EqnSetICState once, since the initial state
     is constant throughout all of the domain. */
    
    // allocate u[sr]
    ierr = xf_Error(xf_Alloc((void **) &u, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // fill in u using IC by calling eqnset-specific function
    ierr = xf_Error(xf_EqnSetICState(EqnSet, IC, IParam, RParam, 1, x0, &Time, u));
    if (ierr != xf_OK) return ierr;
    
    // fill in array data as if basis were Lagrange
    for (i=0; i<U->nArray; i++){
      ga = U->GenArray+i;
      for (j=0; j<ga->n; j++){  
	for (k=0; k<((ga->vr==NULL) ? ga->r : ga->vr[j]); k++)
	  ga->rValue[j][k] = u[k%sr];
      } // j
    } // i
    
    xf_Release((void *) u);

    // convert from Lagrange basis to actual Basis    
    ierr = xf_Error(xf_ConvertVectorFromLagrange(U));
    if (ierr != xf_OK) return ierr;
    
  }
  else if (InterpolateIC){
    
    /*---------------------------*/
    /* Interpolating IC function */
    /*---------------------------*/
    
    /* An IC function is specified.  Need to project/interpolate the
     initial condition onto the space spanned by the basis fcns of
     U.  Here we perform interpolation. */
    
    // allocate vector for storing bases
    ierr = xf_Error(xf_Alloc((void **) &BasisNew, U->nArray, sizeof(enum xfe_BasisType)));
    if (ierr != xf_OK) return ierr;
    
    for (i=0; i<U->nArray; i++){
      // get an appropriate Lagrange basis
      BasisNew[i] = U->Basis[i];
      ierr = xf_Error(xf_Basis2Lagrange(BasisNew[i], U->Basis + i));
      if (ierr != xf_OK) return ierr;
    }
    
    PhiData     = NULL;
    GeomPhiData = NULL;
    xn          = NULL;
    xglob       = NULL;
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

      Basis = U->Basis[egrp];
      pOrder = -1;

      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	
	// get interpolation order
	Order = xf_InterpOrder(U, egrp, elem);
	
	if (Order != pOrder){
	  pOrder = Order;

	  // determine nn = # unknowns for elements in this group
	  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
	  if (ierr != xf_OK) return ierr;
	  
	  // obtain Lagrange node locations
	  ierr = xf_Error(xf_LagrangeNodes(Basis, Order, NULL, NULL, &xn));
	  if (ierr != xf_OK) return ierr;
	  
	  // compute basis functions at Lagrange Nodes
	  ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, nn, xn, xfb_Phi, &PhiData));
	  if (ierr != xf_OK) return ierr;

	  // reallocate xglob
	  ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nn*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	
	}

	// global coordinates
	ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, xfe_True, 
					nn, xn, xglob));
	if (ierr != xf_OK) return ierr;

	if (Altering) return xf_Error(xf_NOT_SUPPORTED);

	// Call eqnset function to determine u
	EU = U->GenArray[egrp].rValue[elem]; //  U on elem [nn*sr]
	ierr = xf_Error(xf_EqnSetICState(EqnSet, IC, IParam, RParam, nn, xglob, &Time, EU));
	if (ierr != xf_OK) return ierr;

      } // elem
      
    } // egrp
    

    // convert from Lagrange basis to actual Basis
    ierr = xf_Error(xf_ProjectVectorInPlace_Basis(Mesh, NULL, U, BasisNew, xfe_BasisLast));
    if (ierr != xf_OK) return ierr;
    
    xf_Release((void *) BasisNew);
    
    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Geometry Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    xf_Release((void *) xglob);
    xf_Release((void *) xn);
  }
  else{
    
    /*-----------------------------------------*/
    /* Least-squares projection of IC function */
    /*-----------------------------------------*/
    
    /* Here we perform least-squares projection of a specified IC
     function. */
    QuadData    = NULL;
    PhiData     = NULL;
    GeomPhiData = NULL;
    u           = NULL;
    xglob       = NULL;
    Q           = NULL;
    R           = NULL;
    P           = NULL;
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      Basis = U->Basis[egrp];

      pOrder = -1;
      pOrder = -1;
      pnq = -1; // here so that we reallocate at least once for every element group
      
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

	// get interpolation order
	Order = xf_InterpOrder(U, egrp, elem);

	if (Order != pOrder){
	  pOrder = Order;

	  // determine nn = # unknowns for elements in this group
	  ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
	  if (ierr != xf_OK) return ierr;

	  // determine quadrature order (not using eqnset)
	  ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*Order, &QuadOrder));
	  if (ierr != xf_OK) return ierr;
	}
	
	/* Pull off quad points for the element; make sure have enough
	   quad points (more than nn). */
	ierr = xf_Error(xf_QuadElemAtLeast(Mesh, egrp, elem, QuadOrder, nn, &QuadData, &QuadChanged));
	if (ierr != xf_OK) return ierr;
	nq = QuadData->nquad;
	xq = QuadData->xquad;

	// re-allocate data if quad points increased
	if (nq > pnq){
	  ierr = xf_Error(xf_ReAlloc( (void **) &u, nq*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}

	// obtain global coords of quad points
	ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
					nq, xq, xglob));
	if (ierr != xf_OK) return ierr;

	if (Altering){
	  // need basis
	  ierr = xf_Error(xf_EvalBasis(U->Basis[egrp], Order, QuadChanged, 
				       nq, xq, xfb_Phi, &PhiData));
	  if (ierr != xf_OK) return ierr;
	  // interpolate current u if Altering
	  xf_MxM_Set(PhiData->Phi, U->GenArray[egrp].rValue[elem], nq, PhiData->nn, sr, u);  
	  // alter state in u
	  ierr = xf_Error(xf_EqnSetAlterState(EqnSet, FcnName, IParam, RParam, xglob, 
					      FParam, nq, u));
	  if (ierr != xf_OK) return ierr;
	}
	else{
	  // Call eqnset function to determine u
	  ierr = xf_Error(xf_EqnSetICState(EqnSet, IC, IParam, RParam, nq, xglob, &Time, u));
	  if (ierr != xf_OK) return ierr;
	}	

	// Project u to get coefficients EU
       	ierr = xf_Error(xf_ProjectOnElemQR(Mesh, egrp, elem, sr, QuadData, 
                                           Basis, Order, QuadChanged, &Q, &R, &P,
                                           u, U->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr;
        
        pnq = nq;
      } // elem
    } // egrp
    
    // Only destroy QuadData if points are generic
    ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Geometry Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    xf_Release((void *) u);
    xf_Release((void *) xglob);
    xf_Release((void *) Q);
    xf_Release((void *) R);
    xf_Release((void *) P);
  }
  
  xf_Release((void *) IParam);
  xf_Release((void *) RParam);
  xf_Release((void *) FParam);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_InitState
int 
xf_InitState(xf_All *All, xf_Vector *U)
{
  int ierr;
  real CFL;
  enum xfe_Bool UseGCL = xfe_False;
  xf_ICs *ICs;
  xf_IC *IC;
  xf_MeshMotion *Motion;  

  ICs = All->EqnSet->ICs;
  
  if (ICs->nIC <= 0){
    xf_printf("Error, no initial conditions present!\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  IC  = ICs->IC+0;
  
  
  // initialize with steady solve?  If so, run solver. (initial U assumed valid)
  if (IC->PriorSteadySolve){
    // do not want to change CFL, so store it 
    ierr = xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFL);
    if (ierr != xf_OK) return ierr;
    // do not want mesh motion on for the steady solve
    Motion = All->Mesh->Motion;
    All->Mesh->Motion = NULL;
    xf_printf("Running steady solve for initialization.\n");
    
    ierr = xf_Error(xf_SolveNonlinearSystem(All, 0, xfe_False, NULL, &U));
    if (ierr != xf_OK){
      xf_printf("Error occured during steady solve!\n");
      return ierr;
    }
    // reset CFL to original
    ierr = xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFL);
    if (ierr != xf_OK) return ierr;
    // reset Mesh motion
    All->Mesh->Motion = Motion;

  }
  else{
    // straight-up initialize if no prior solve
    ierr = xf_Error(xf_InitAlterState(All, NULL, NULL, U));
    if (ierr != xf_OK) return ierr;
  }
  
  
  // Alter steady state?
  if ((IC->AlterFunction != NULL) && (xf_NotNull(IC->AlterFunction)) ){
    ierr = xf_Error(xf_InitAlterState(All, IC->AlterFunction, IC->AlterData, U));
    if (ierr != xf_OK) return ierr;
  }

  // If using GCL, find GCL vector and initialize it
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U, -1, xfe_True, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_AlterState
int 
xf_AlterState(xf_All *All, const char *FcnName, const char *FcnParam,
              xf_Vector *U)
{
  return xf_Error(xf_InitAlterState(All, FcnName, FcnParam, U));
}



/******************************************************************/
//   FUNCTION Definition: xf_ChangeVariableSet
int 
xf_ChangeVariableSet(xf_All *All, const char *VariableSet, xf_Vector *U)
{
  int ierr, sr, nn, egrp, elem;
  int nq, pnq, Order, pOrder, QuadOrder;
  int *IParam = NULL, *P;
  enum xfe_Bool QuadChanged;
  enum xfe_BasisType Basis;
  real *u, *v, *xq;
  real *Q, *R;
  real *EU, *RParam = NULL;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  
  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;
  sr     = EqnSet->StateRank;
  
  if (U->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  
  if (U->StateRank != sr) return xf_Error(xf_OUT_OF_BOUNDS);
  
  // build eqnset-desired parameter lists 
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  /* Here we perform least-squares projection of the transformed
   variables to obtain the desired coefficients. */
  QuadData    = NULL;
  PhiData     = NULL;
  u           = NULL;
  v           = NULL;
  Q           = NULL;
  R           = NULL;
  P           = NULL;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    pOrder = -1;
    Basis = U->Basis[egrp];
    pnq = -1; // here so that we reallocate at least once for every element group

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);
      
      if (Order != pOrder){
	pOrder = Order;
	
	// determine nn = # unknowns for elements in this group
	ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
	if (ierr != xf_OK) return ierr;
	
	// determine quadrature order (not using eqnset)
	ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*Order, &QuadOrder));
	if (ierr != xf_OK) return ierr;
      }
          
      /* Pull off quad points for the element; make sure have enough
       quad points (more than nn). */
      ierr = xf_Error(xf_QuadElemAtLeast(Mesh, egrp, elem, QuadOrder, nn, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **) &u, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &v, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        pnq = nq;
      }
      
      // Calculate U at xq -> u
      EU = U->GenArray[egrp].rValue[elem];
      xf_MxM_Set(PhiData->Phi, EU, nq, PhiData->nn, sr, u);     
      
      // Call eqnset function to determine new variables
      ierr = xf_Error(xf_EqnSetVariableChange(EqnSet, VariableSet, IParam, RParam, nq, u, v,
                                              NULL, NULL));
      if (ierr != xf_OK) return ierr;
      
      // Project v to get coefficients -> EU
      ierr = xf_Error(xf_ProjectOnElemQR(Mesh, egrp, elem, sr, QuadData, 
                                         Basis, Order, QuadChanged, &Q, &R, &P, 
                                         v, EU));
      if (ierr != xf_OK) return ierr;
      
    } // elem
  } // egrp
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  xf_Release((void *) u);
  xf_Release((void *) v);
  xf_Release((void *) IParam);
  xf_Release((void *) RParam);
  xf_Release((void *) Q);
  xf_Release((void *) R);
  xf_Release((void *) P);
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScaleState
int 
xf_ScaleState(xf_All *All, xf_ICs *ICsOrig, xf_Vector *U)
{
  int ierr, sr, nn, egrp, elem;
  int nq, pnq, Order, QuadOrder;
  int *IParam = NULL, *P;
  enum xfe_Bool QuadChanged;
  enum xfe_BasisType Basis;
  real *u, *xq;
  real *Q, *R;
  real *EU, *RParam = NULL;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  xf_IC *ICOrig;
  
  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;
  sr     = EqnSet->StateRank;
  
  if (U->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
  if (U->StateRank != sr) return xf_Error(xf_OUT_OF_BOUNDS);
  if (ICsOrig == NULL)  return xf_Error(xf_INPUT_ERROR);
  if (ICsOrig->nIC <= 0){
    xf_printf("No original IC from which to scale.  Continuing without scaling.\n");
    return xf_OK;
  }
  ICOrig  = ICsOrig->IC+0;
  
  // build eqnset-desired parameter lists
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  /* Here we perform least-squares projection of the scaled
   variables to obtain the desired coefficients. */
  QuadData    = NULL;
  PhiData     = NULL;
  u           = NULL;
  Q           = NULL;
  R           = NULL;
  P           = NULL;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    Basis = U->Basis[egrp];
    
    pnq = -1; // here so that we reallocate at least once for every element group
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);

      // determine nn = # unknowns for elements in this group
      ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
      if (ierr != xf_OK) return ierr;
      
      // determine quadrature order (not using eqnset)
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;

      
      /* Pull off quad points for the element; make sure have enough
       quad points (more than nn). */
      ierr = xf_Error(xf_QuadElemAtLeast(Mesh, egrp, elem, QuadOrder, nn, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      
      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      // re-allocate data if quad points increased
      if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **) &u, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        pnq = nq;
      }
      
      // Calculate U at xq -> u
      EU = U->GenArray[egrp].rValue[elem];
      xf_MxM_Set(PhiData->Phi, EU, nq, PhiData->nn, sr, u);     
      
      // Call eqnset function to determine new variables (overwrite u)
      ierr = xf_Error(xf_EqnSetScaleState(EqnSet, ICOrig, IParam, RParam, nq, u));
      if (ierr != xf_OK) return ierr;
      
      // Project v to get coefficients -> EU
      ierr = xf_Error(xf_ProjectOnElemQR(Mesh, egrp, elem, sr, QuadData, Basis, Order, 
                                         QuadChanged, &Q, &R, &P, u, EU));
      if (ierr != xf_OK) return ierr;
      
    } // elem
  } // egrp
  
  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  xf_Release((void *) u);
  xf_Release((void *) IParam);
  xf_Release((void *) RParam);
  xf_Release((void *) Q);
  xf_Release((void *) R);
  xf_Release((void *) P);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReconstructVector
int 
xf_ReconstructVector(xf_All *All, int OrderIncrement, xf_Vector *U,
                     enum xfe_Bool AllocFlag, xf_Vector **pV)
{
  int ierr, sr, nn, dim, d;
  int iq, nq, iqtot, nqtot;
  int egrp, elem, face;
  int egrpN, elemN;
  int k, nskip;
  int nface, nqperelem;
  int Order, Orderh, QuadOrder;
  int *P;
  enum xfe_BasisType Basis, Basish;
  enum xfe_Bool SkipFlag, QuadChanged, converged;
  real *xq, xglob[3], xref[3], *pxref;
  real *u, *Q, *R;
  xf_QuadData *QuadData, *QuadDatah;
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_Vector *V;
  xf_EqnSet *EqnSet;
  xf_Mesh *Mesh;
  
  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;
  dim    = Mesh->Dim;
  sr     = EqnSet->StateRank;
  
  if (AllocFlag){
    // allocate HO vector
    ierr = xf_Error(xf_FindSimilarVectorHO(All, U, "VReconstruct", xfe_True, xfe_False,
                                           &OrderIncrement, NULL, NULL,  pV));
    if (ierr != xf_OK) return ierr;
  }
  V = (*pV);
  
  if (Mesh->ParallelInfo != NULL){
    // parallel communicate U so halos are updated
    ierr = xf_Error(xf_HaloExchangeVectorBegin(U));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
    if (ierr != xf_OK) return ierr;
  }    
  
  // create QuadDatah, used for least squares
  ierr = xf_Error(xf_CreateQuadData(&QuadDatah));
  if (ierr != xf_OK) return ierr;
  
  QuadData    = NULL;
  PhiData     = NULL;
  JData       = NULL;
  u           = NULL;
  Q           = NULL;
  R           = NULL;
  P           = NULL;
  nqtot = 0; // total nubmer of quad points allocated for QuadDatah
  
  
  // Loop over groups, perform least-squares projection into higher order
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    // fine space basis and order
    Basish = V->Basis[egrp];

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // get interpolation order
      Orderh = xf_InterpOrder(V, egrp, elem);
      
      // determine nn = # unknowns/elem for fine-space in this group
      ierr = xf_Error(xf_Order2nNode(Basish, Orderh, &nn));
      if (ierr != xf_OK) return ierr;
      
      // determine quadrature order (not using eqnset)
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*Orderh, &QuadOrder));
      if (ierr != xf_OK) return ierr;    

      iqtot = 0; // total number of quad points for least squares on this elem
      
      // number of valid faces (for deciding # quad points) -- NOT USED
      for (face=0, nface=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++)
        if (Mesh->ElemGroup[egrp].Face[elem][face].Group == xf_INTERIORFACE) nface++;
      // minimum number of quad points per elem
      nqperelem = 4*(nn/(nface+1))+1;
      
      // loop over self + neighbor elements
      for (face=-1; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
	if (face == -1){
	  egrpN = egrp;
	  elemN = elem;
	}
	else{
	  ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, &egrpN, &elemN, NULL));
	  if (ierr != xf_OK) return ierr;
	}

	// skip if on boundary
	if (egrpN < 0) continue;

	// original basis and order on egrpN
	Basis = U->Basis[egrpN];
	Order = xf_InterpOrder(U, egrpN, elemN);

	/* Pull off quad points for the element; make sure have enough
	   quad points (more than nn). */
	ierr = xf_Error(xf_QuadElemAtLeast(Mesh, egrpN, elemN, QuadOrder, nn,
					   &QuadData, &QuadChanged));
	if (ierr != xf_OK) return ierr;

	nq = QuadData->nquad;
	xq = QuadData->xquad;

	/* Compute geometry Jacobian; if not constant, compute at quad
	   points.  Note if jacobian is constant, only one Jacobian will
	   be computed/returned. */
	ierr = xf_Error(xf_ElemJacobian(Mesh, egrpN, elemN, nq, xq, xfb_detJ, QuadChanged, &JData));
	if (ierr != xf_OK) return ierr;

	// compute original basis functions if quad or basis or order changed
	ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, 
				     xfb_Phi, &PhiData));
	if (ierr != xf_OK) return ierr;

	// reallocate QuadDatah and u if necessary
	if (iqtot+nq > nqtot){
	  nqtot = iqtot+nq;
	  ierr = xf_Error(xf_ReAlloc( (void **) &QuadDatah->xquad, nqtot*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;  
	  ierr = xf_Error(xf_ReAlloc( (void **) &QuadDatah->wquad, nqtot*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;  
	  ierr = xf_Error(xf_ReAlloc( (void **) &u, nqtot*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}

	// interpolate U at xq -> u
	xf_MxM_Set(PhiData->Phi, U->GenArray[egrpN].rValue[elemN], 
		   nq, PhiData->nn, sr, u+iqtot*sr);

	// obtain ref coords of xq, from point of view of elem -> QuadDatah->xquad
	for (iq=0, nskip=0; iq<nq; iq++){
	  SkipFlag = xfe_False;
	  pxref = xq + iq*dim;
	  if (face >= 0){
	    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrpN, elemN, NULL, xfe_True, 
					    1, pxref, xglob));
	    if (ierr != xf_OK) return ierr;
	    ierr = xf_Error(xf_Glob2RefElem(Mesh, egrp, elem, xglob, -1.0, xfe_True,
					    xref, &converged));
	    if (ierr != xf_OK) return ierr;
	    if (!converged) SkipFlag = xfe_True; // will skip this point
	    pxref = xref;
	  }
	  if (!SkipFlag){ // store in QuadDatah
	    for (d=0; d<dim; d++) QuadDatah->xquad[iqtot*dim+d] = pxref[d];
	    // use detJ-multiplied quad weight
	    QuadDatah->wquad[iqtot] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
	    // copy over u to earlier entry 
	    if (nskip > 0) for (k=0; k<sr; k++) u[iqtot*sr+k] = u[iqtot*sr+nskip*sr+k];
	    iqtot++;
	  }
	  else nskip++;
	} // iq

      } // face
      
      QuadDatah->nquad = iqtot;
      
      // Project u to get coefficients EU
      ierr = xf_Error(xf_ProjectOnElemLeastSquares(Mesh, egrp, elem, sr, QuadDatah,
                                                   Basish, Orderh, xfe_True, &Q, &R, &P,
                                                   u, V->GenArray[egrp].rValue[elem]));
      if (ierr != xf_OK) return ierr;
      
    } // elem
  } // egrp
  
  
  // Destroy QuadData and QuadDatah
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadDatah));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  xf_Release((void *) u);
  xf_Release((void *) Q);
  xf_Release((void *) R);
  xf_Release((void *) P);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorTextOut
int
xf_VectorTextOut(const char *fname, enum xfe_Bool Project2Lagrange, 
                 xf_Vector *U)
{
  int ierr, i;
  int egrp, r, nelem, elem, k, Order;
  enum xfe_BasisType *Basis;
  FILE *fid;
  xf_GenArray *ga;
  
  if ( (fid = fopen(fname,"w")) == NULL )
    return xf_Error(xf_FILE_WRITE_ERROR);

  // if data is not interpolated set a default p=0 basis + order
  if ((U->Basis == NULL) || (U->Order == NULL)){
    xf_Release((void *) U->Basis);
    xf_Release((void *) U->Order);
    ierr = xf_Error(xf_Alloc((void **) &U->Basis, U->nArray, sizeof(enum xfe_BasisType)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &U->Order, U->nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<U->nArray; i++){
      // irrelevant for p = 0
      U->Basis[i] = xfe_SegLagrange;
      U->Order[i] = 0;
    }

  }
  
  if (Project2Lagrange){
    // convert U to a default Lagrange basis (to make .txt readers simpler)
    ierr = xf_Error(xf_Alloc( (void **) &Basis, U->nArray, sizeof(enum xfe_BasisType)));
    if (ierr != xf_OK) return ierr;
    for (egrp=0; egrp<U->nArray; egrp++){
      ierr = xf_Error(xf_Basis2UniformLagrange(U->Basis[egrp], Basis+egrp));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_ProjectVectorInPlace_Basis(NULL, NULL, U, Basis, xfe_BasisLast));
    if (ierr != xf_OK) return ierr;
    xf_Release( (void *) Basis);
  }

  // write out to text file
  fprintf(fid, "%% nElemGroup = %d\n", U->nArray);
  for (egrp=0; egrp<U->nArray; egrp++){
    ga = U->GenArray + egrp;
    r = ga->r;
    nelem = ga->n;
    fprintf(fid, "%% ElemGroup %d\n%%  nElem = %d\n%%  Basis = %s\n%%  max_Order = %d\n", 
	    egrp, nelem, xfe_BasisName[U->Basis[egrp]], U->Order[egrp]);
    fprintf(fid, "%%  max_r = %d\n", r);
    fprintf(fid, "%%  max_dof_per_elem = %d\n", r/U->StateRank);
    fprintf(fid, "%%  StateRank = %d\n", U->StateRank);
    fprintf(fid, "%%  StateNames = < ");
    if (U->StateName != NULL){
      for (k=0; k<U->StateRank; k++)
	fprintf(fid, "%s ", U->StateName[k]);
    }
    else{
      fprintf(fid, "Value ");
    }
    fprintf(fid, " >\n");
    fprintf(fid, "%%  The following numbers are: Order r; <values>\n\n");
    for (elem=0; elem<nelem; elem++){
      Order = xf_InterpOrder(U, egrp, elem);
      r = ((ga->vr == NULL) ? ga->r : ga->vr[elem]);
      fprintf(fid, "%d %d\n", Order, r);
      if (U->GenArray[egrp].rValue != NULL)
	for (k=0; k<r; k++)
	  fprintf(fid, "%.15E ", U->GenArray[egrp].rValue[elem][k]);
      else if (U->GenArray[egrp].iValue != NULL)
	for (k=0; k<r; k++)
	  fprintf(fid, "%d ", U->GenArray[egrp].iValue[elem][k]);
      else
	return xf_Error(xf_INPUT_ERROR);
      fprintf(fid, "\n");
    }
  }
  
  fclose(fid);
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DumpSystemSparse
int
xf_DumpSystemSparse(xf_All *All, xf_Vector *U)
{
  //  Dumps system to ascii text file in sparse format
  int ierr, egrp, elem, face, nface;
  int egN, eN, faceN, negrp, nn, nnN;
  int Order, pOrder;
  int myRank, nProc;
  int globindex;
  int *egrpOffset;
  int r, rN, row, col;
  int i, j, ij, k, sr, sr2;
  int iOutput, nOutput;
  char fname[xf_MAXSTRLEN], AdjointOutputs[xf_MAXSTRLEN];
  char **OutputNames = NULL;
  real *MM, fac, val, J, *A;
  FILE *fid;
  xf_JacobianMatrix *R_U;
  xf_Vector *R;
  xf_Vector *J_U;
  xf_Matrix *M;
  xf_Mesh *Mesh;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc != 1){
    xf_printf("Warning, dumping system to .txt not supported in parallel.");
    return xf_OK;
  }
  
  // locate Jacobian vector
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_True, NULL, &R_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  // locate Residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_False, xfe_True, 
                                       NULL, &R, NULL));
  if (ierr != xf_OK) return ierr;
  
  // calculate the residual
  ierr = xf_Error(xf_CalculateResidual(All, U, R, R_U, NULL));
  if (ierr != xf_OK);
  
  Mesh = All->Mesh;
  sr = All->EqnSet->StateRank;
  sr2 = sr*sr;
  negrp = Mesh->nElemGroup;
  
  ierr = xf_Error(xf_Alloc( (void **) &egrpOffset, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // loop over element groups, calculate global index offsets
  globindex = 0;
  for (egrp=0; egrp<negrp; egrp++){
    egrpOffset[egrp] = globindex;
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      nn = xf_Jacobian_n(R_U,egrp,elem);
      globindex += nn*sr;
    }
  }
  
  // loop over groups, dump out Jacobian
  if ((fid = fopen("A.txt", "w")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  for (egrp=0; egrp<negrp; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      nface = Mesh->ElemGroup[egrp].nFace[elem];
      for (face=-1; face<nface; face++){
	if (face == -1){
	  egN   = egrp;
	  eN    = elem;
	  faceN = face;
	}
	else{
	  egN   = R_U->egrpN[egrp][elem][face];
	  eN    = R_U->elemN[egrp][elem][face];
	  faceN = R_U->faceN[egrp][elem][face];
	}

	if (egN < 0) continue; // indicates a boundary face

	A = R_U->Value[egrp][elem][1+face];
	nn = xf_Jacobian_n(R_U,egrp,elem);
	r  = sr * nn;
	nnN = xf_Jacobian_n(R_U,egN,eN);
	rN  = sr * nnN;
	for (row=egrpOffset[egrp],i=0; i<elem; i++) row += sr*xf_Jacobian_n(R_U,egrp,i);
	for (col=egrpOffset[egN] ,i=0; i<eN  ; i++) col += sr*xf_Jacobian_n(R_U,egN, i);	

	for (i=0; i<nn; i++)
	  for (j=0; j<nnN; j++)
	    for (k=0, ij = i*nnN+j; k<sr2; k++)
	      if ((val = A[ij*sr2+k]) != 0.0)
		fprintf(fid, "%d %d %.15E\n", row+i*sr+(k/sr)+1, col+j*sr+(k%sr)+1, val);
      } // face
    } // elem
    if (egrp == (negrp-1)){
      fprintf(fid, "%% %d\n", globindex);
    }
  }
  fclose(fid);
  
  // loop over groups, dump out Mass Matrix
  if ((fid = fopen("M.txt", "w")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  for (egrp=0; egrp<negrp; egrp++){

    pOrder = -1;

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      Order = xf_InterpOrder(U, egrp, elem); // get interpolation order
      
      if (pOrder != Order){
        pOrder = Order;
        ierr = xf_Error(xf_FindMassMatrixData(All, egrp, U->Basis[egrp], 
                                              Order, &M));
        if (ierr != xf_OK) return ierr;
      }
      
      Order = xf_InterpOrder(U, egrp, elem); // get interpolation order

      if (pOrder != Order){
	pOrder = Order;
	ierr = xf_Error(xf_FindMassMatrixData(All, egrp, U->Basis[egrp], 
					      Order, &M));
	if (ierr != xf_OK) return ierr;
      }

      ierr = xf_Error(xf_ElemMassMatrix(All, egrp, elem, U->Basis[egrp], 
					Order, M, NULL, &MM, &fac));
      if (ierr != xf_OK) return ierr;

      nn = xf_Jacobian_n(R_U,egrp,elem);
      r  = sr * nn;
      for (row=egrpOffset[egrp],i=0; i<elem; i++) row += sr*xf_Jacobian_n(R_U,egrp,i);
      col = row;
      
      nn = xf_Jacobian_n(R_U,egrp,elem);
      r  = sr * nn;
      for (row=egrpOffset[egrp],i=0; i<elem; i++) row += sr*xf_Jacobian_n(R_U,egrp,i);
      col = row;
      
      for (i=0; i<nn; i++)
        for (j=0; j<nn; j++)
          if ((val = fac*MM[i*nn+j]) != 0.0)
            for (k=0; k<sr; k++)
              fprintf(fid, "%d %d %.15E\n", row+i*sr+k+1, col+j*sr+k+1, val);
      
    } // elem
    if (egrp == (negrp-1)){
      fprintf(fid, "%% %d\n", globindex);
    }
  } // egrp
  fclose(fid);
  
  xf_Release( (void *) egrpOffset);


  // loop over Adjoint outputs, dump out output sensitivities
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdjointOutputs", AdjointOutputs));
  if (ierr != xf_OK) return ierr;

  if (xf_NotNull(AdjointOutputs)){

    /* Only use outputs in AdjointOutputs */
    ierr = xf_Error(xf_ScanXStringAlloc(AdjointOutputs, xf_MAXSTRLEN, &nOutput, &OutputNames));
    if (ierr != xf_OK) return ierr;
    
    // create J_U
    ierr = xf_Error(xf_FindSimilarVector(All, U, "J_U", xfe_False, xfe_True, NULL, &J_U, NULL));
    if (ierr != xf_OK) return ierr;  
    

    for (iOutput=0; iOutput<nOutput; iOutput++){
      
      // compute J_U
      ierr = xf_Error(xf_CalculateOutput(All, OutputNames[iOutput], U, &J, J_U, xfe_Set));
      if (ierr != xf_OK) return ierr;

      sprintf(fname, "%s_U.txt", OutputNames[iOutput]);
      ierr = xf_Error(xf_VectorTextOut(fname, xfe_False, J_U));
      if (ierr != xf_OK) return ierr;

    } // iOutput
    
    xf_Release2((void **) OutputNames);  

  } // if AdjointOutputs not null

  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GetTimeSchemeInfo
int
xf_GetTimeSchemeInfo(enum xfe_TimeSchemeType TimeScheme, 
                     enum xfe_Bool *Explicit, enum xfe_Bool *MultiStep,
		     enum xfe_Bool *MultiStage)
{
  switch (TimeScheme){
    case xfe_TimeSchemeSteady:
    case xfe_TimeSchemeBDF1:
    case xfe_TimeSchemeBDF2:
    case xfe_TimeSchemeBDF3:
    case xfe_TimeSchemeTrapezoidal:
      if (Explicit   != NULL) (*Explicit)   = xfe_False;
      if (MultiStep  != NULL) (*MultiStep)  = xfe_True;
      if (MultiStage != NULL) (*MultiStage) = xfe_False;
      break;
    case xfe_TimeSchemeFE:
    case xfe_TimeSchemeRK4:
    case xfe_TimeSchemeRK2:
      if (Explicit   != NULL) (*Explicit)   = xfe_True;
      if (MultiStep  != NULL) (*MultiStep)  = xfe_True;
      if (MultiStage != NULL) (*MultiStage) = xfe_True;
      break;
    case xfe_TimeSchemeDG1:
    case xfe_TimeSchemeDG2:
      if (Explicit   != NULL) (*Explicit)   = xfe_False;
      if (MultiStep  != NULL) (*MultiStep)  = xfe_False;
      if (MultiStage != NULL) (*MultiStage) = xfe_False;
      break;
    case xfe_TimeSchemeESDIRK3:
    case xfe_TimeSchemeESDIRK4:
    case xfe_TimeSchemeESDIRK5:
    case xfe_TimeSchemeDIRK3:
    case xfe_TimeSchemeDIRK4:
      if (Explicit   != NULL) (*Explicit)   = xfe_False;
      if (MultiStep  != NULL) (*MultiStep)  = xfe_True;
      if (MultiStage != NULL) (*MultiStage) = xfe_True;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateTimeHistData
int 
xf_CreateTimeHistData( xf_TimeHistData **pTimeHistData)
{
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pTimeHistData, 1, sizeof(xf_TimeHistData)));
  if (ierr != xf_OK) return ierr;
  
  (*pTimeHistData)->nTime         = 0;
  (*pTimeHistData)->Time          = NULL;
  (*pTimeHistData)->TimeStep      = NULL;
  (*pTimeHistData)->TimeScheme    = NULL;
  (*pTimeHistData)->ConstTimeStep = xfe_False;
  (*pTimeHistData)->nOutput       = 0;
  (*pTimeHistData)->OutputNames   = NULL;
  (*pTimeHistData)->OutputValues  = NULL;
  (*pTimeHistData)->TimeWeights   = NULL;
  (*pTimeHistData)->nSumWeight       = 0;
  (*pTimeHistData)->SumOutputWeights = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyTimeHistData
int 
xf_DestroyTimeHistData( xf_TimeHistData *TimeHistData)
{
  if (TimeHistData == NULL) return xf_OK;
  
  xf_Release(  (void  *) TimeHistData->Time);
  xf_Release(  (void  *) TimeHistData->TimeStep);
  xf_Release(  (void  *) TimeHistData->TimeScheme);
  xf_Release2( (void **) TimeHistData->OutputNames);
  xf_Release2( (void **) TimeHistData->OutputValues);
  xf_Release(  (void  *) TimeHistData->TimeWeights);
  xf_Release2( (void **) TimeHistData->SumOutputWeights);
  
  xf_Release(  (void  *) TimeHistData);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadTimeHistDataSerial
static int 
xf_ReadTimeHistDataSerial( const char *fname, xf_TimeHistData *TimeHistData)
{
  /*
   
   PURPOSE: 
   
   Reads an ASCII TimeHistData file on one processor.
   
   INPUTS:
   
   fname : name of TimeHistData file
   
   OUTPUTS: 
   
   TimeHistData : filled-in time history.  Outermost data structure must be
   pre-allocated.
   
   RETURNS: Error Code
   
   */
  
  int ierr, itemp, i, nTime;
  char line0[xf_MAXLINELEN], line1[xf_MAXLINELEN], *line;
  char TimeSchemeMap[xfe_TimeSchemeLast][xf_MAXSTRLEN];
  char key[xf_MAXSTRLEN];
  enum xfe_Bool done;
  FILE *fid;
  
  if ((fid = fopen(fname, "r")) == NULL)
    return xf_Error(xf_FILE_READ_ERROR);
  
  // read in first line
  if (fgets(line0, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  done = xfe_False;
  do{
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL) continue;
    if (strncmp(line0, "% Time schemes", 14)==0){ // start of time scheme map
      for (i=0; i<xfe_TimeSchemeLast; i++){
        if (fgets(line0, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
        line = line0+1;
        /* if blank line or comment then error */
        if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) return xf_Error(xf_FILE_READ_ERROR);
        /* read key and value from line */
        ierr = xf_Error(xf_ReadKey(line, "=", key, TimeSchemeMap[i], xf_MAXSTRLEN));
        if (ierr != xf_OK) return ierr;
      }
      done = xfe_True;
      break;
    }
  } while(feof(fid) == 0);
  if (!done) return xf_Error(xf_FILE_READ_ERROR);
  // header line
  if (fgets(line0, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  
  // preallocate TimeHistData (initial size of nTime)
  nTime = 10;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->Time, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeStep, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeScheme, nTime,
                             sizeof(enum xfe_TimeSchemeType)));
  if (ierr != xf_OK) return ierr;
  
  // read through time lines in file
  TimeHistData->nTime = 0;
  do{
    if (fgets(line0, xf_MAXLINELEN, fid) == NULL) continue;
    i = TimeHistData->nTime++;
    // reallocate TimeHistData if needed
    if (TimeHistData->nTime > nTime){
      nTime *= 2;
      ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->Time, nTime, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeStep, nTime, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeScheme, nTime,
                                 sizeof(enum xfe_TimeSchemeType)));
      if (ierr != xf_OK) return ierr;
    }
    // fill Time, TimeStep, TimeScheme (not super robust because using int for TimeScheme)
    ierr = sscanf(line0, "%d %lf %lf %d %s", &itemp, TimeHistData->Time+i,
                  TimeHistData->TimeStep+i, (int *) TimeHistData->TimeScheme+i, line1);
    if ((ierr != 5) & (ierr != 4)) return xf_Error(xf_FILE_READ_ERROR);
  } while(feof(fid) == 0);
  
  // final reallocate of TimeHistData
  nTime = TimeHistData->nTime;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->Time, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeStep, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeScheme, nTime,
                             sizeof(enum xfe_TimeSchemeType)));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadTimeHistData
int 
xf_ReadTimeHistData( const char *fname, const char *LogOutput,
                    xf_TimeHistData **pTimeHistData)
{
  int ierr, myRank, nProc;
  int nTime;
  xf_TimeHistData *TimeHistData;
  
  xf_printf("Reading TimeHistData file: %s\n", fname);
  
  // allocate structure for TimeHistData
  ierr = xf_Error(xf_CreateTimeHistData(pTimeHistData));
  if (ierr != xf_OK) return ierr;
  TimeHistData = (*pTimeHistData);
  
  // check parallel
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // root does the reading
  if (myRank == 0)
    ierr = xf_Error(xf_ReadTimeHistDataSerial(fname, TimeHistData));
  if (xf_PError(&ierr, 0) != xf_OK) return ierr;
  
  // broadcast nTime
  ierr = xf_Error(xf_MPI_Bcast((void *) &TimeHistData->nTime, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  nTime = TimeHistData->nTime;
  
  // other procs allocate space
  if (myRank > 0){
    ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->Time, nTime, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeStep, nTime, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeScheme, nTime,
                               sizeof(enum xfe_TimeSchemeType)));
    if (ierr != xf_OK) return ierr;
  }
  // broadcast from root
  ierr = xf_Error(xf_MPI_Bcast((void *) TimeHistData->Time, nTime*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Bcast((void *) TimeHistData->TimeStep, nTime*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Bcast((void *) TimeHistData->TimeScheme, 
                               nTime*sizeof(enum xfe_TimeSchemeType), 0));
  if (ierr != xf_OK) return ierr;
  
  // allocate room for outputs
  if ( (LogOutput != NULL) && (xf_NotNull(LogOutput)) ){
    ierr = xf_Error(xf_ScanXStringAlloc(LogOutput, xf_MAXSTRLEN, 
                                        &TimeHistData->nOutput, 
                                        &TimeHistData->OutputNames));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc2( (void ***) &TimeHistData->OutputValues, 
                                TimeHistData->nOutput, nTime, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
  
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateUniformTimeHistData
int 
xf_CreateUniformTimeHistData( xf_All *All, const char *LogOutput,
			      xf_TimeHistData **pTimeHistData)
{
  int ierr, nTimeStep, nTime, i;
  enum xfe_Bool MultiStep, MultiStage;
  enum xfe_TimeSchemeType TimeScheme;
  real StartTime, EndTime, TimeStep;
  xf_TimeHistData *TimeHistData;
  
  // allocate structure for TimeHistData
  ierr = xf_Error(xf_CreateTimeHistData(pTimeHistData));
  if (ierr != xf_OK) return ierr;
  TimeHistData = (*pTimeHistData);
  
  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
                                     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
                                     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;
  
  // Start Time = Current Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &StartTime));
  if (ierr != xf_OK) return ierr;
 
  // End Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "EndTime", &EndTime));
  if (ierr != xf_OK) return ierr;
  
  // Number of time steps
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nTimeStep", &nTimeStep));
  if (ierr != xf_OK) return ierr;
  
  // Make sure input is valid
  if ((EndTime <= StartTime) || (nTimeStep < 0)) return xf_Error(xf_INPUT_ERROR);
  if (nTimeStep == 0){
    xf_printf("Warning, nTimeStep=0, not allocating time history structure.\n");
    return xf_OK; // nothing to do
  }
  
  // uniform spacing means constant time step
  TimeHistData->ConstTimeStep = xfe_True;
  
  // Calculate time step
  TimeStep = (EndTime - StartTime)/nTimeStep;
  
  // Is Time Scheme MultiStep or MultiStage?
  ierr = xf_Error(xf_GetTimeSchemeInfo(TimeScheme, NULL, &MultiStep, &MultiStage));
  if (ierr != xf_OK) return ierr;
  
  if (MultiStep || MultiStage) 
    nTime = nTimeStep+1; // number of time nodes for multistep
  else
    nTime = nTimeStep;   // number of time slabs for DG in time
  
  
  // set nTime, allocate Time, TimeStep, TimeScheme
  TimeHistData->nTime = nTime;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->Time, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeStep, nTime, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &TimeHistData->TimeScheme, nTime,
                             sizeof(enum xfe_TimeSchemeType)));
  if (ierr != xf_OK) return ierr;
  
  // fill Time, TimeStep, TimeScheme
  for (i=0; i<nTime; i++) TimeHistData->Time[i]       = StartTime + i*TimeStep;
  for (i=0; i<nTime; i++) TimeHistData->TimeStep[i]   = TimeStep;
  for (i=0; i<nTime; i++) TimeHistData->TimeScheme[i] = TimeScheme;
  
  // special first time step rule if MultiStep
  if (MultiStep){
    if (TimeScheme == xfe_TimeSchemeBDF2)
      TimeHistData->TimeScheme[1] = xfe_TimeSchemeBDF1; // i=0 not used
    if (TimeScheme == xfe_TimeSchemeBDF3){
      if (nTime < 3) return xf_Error(xf_INPUT_ERROR);
      TimeHistData->TimeScheme[1] = xfe_TimeSchemeBDF1;
      TimeHistData->TimeScheme[2] = xfe_TimeSchemeBDF2;
    }
  }
  
  // if LogOutput given: fill in nOutput, OutputNames; allocate OutputValues
  if ( (LogOutput != NULL) && (xf_NotNull(LogOutput)) ){
    ierr = xf_Error(xf_ScanXStringAlloc(LogOutput, xf_MAXSTRLEN, 
                                        &TimeHistData->nOutput, 
                                        &TimeHistData->OutputNames));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc2( (void ***) &TimeHistData->OutputValues, 
                                TimeHistData->nOutput, nTime, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<TimeHistData->nOutput*nTime; i++) TimeHistData->OutputValues[0][i] = 0.;
  }
  
  return xf_OK;
}

/******************************************************************/
// counterpart of xf_InitAlterState
// sr: number of ranks
// only support uniform initialization
int
xfYu_InitAlterState(xf_All *All, Yu_Model *Model, char **FcnName, real *FcnParam, xf_Vector *U, 
                    const int sr)
{
    int ierr, i, j, k, dim, nn, egrp, elem;
    int iq, nq, pnq, pOrder, Order, QuadOrder;
    int nFcnParam, InterpolateInitData;
    int *IParam, *P;
    enum xfe_Bool QuadChanged, InterpolateIC;
    enum xfe_Bool Altering = xfe_False;
    enum xfe_BasisType Basis;
    enum xfe_BasisType *BasisNew; 
    real *pu, Gamma;
    real *u, *xq, *xglob, x0[3] = {0.};
    real *xn, *EU, *RParam;
    real *Q, *R;
    real *FParam = NULL;
    real Time;
    xf_GenArray *ga;
    xf_BasisData *PhiData, *GeomPhiData;
    xf_QuadData *QuadData;
    xf_Data *GammaDat;
    xf_Vector *GammaVec;
    xf_ICs *ICs;
    xf_IC  *IC;
    xf_Mesh *Mesh;
  

    Mesh   = All->Mesh;
    dim    = Mesh->Dim;
    InterpolateInitData = Model->InterpolateInitData;

    //find Gamma Vector
    ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
    if(ierr == xf_NOT_FOUND)
    {
       xf_printf("Cannot find heat capacity ratio...\n");
       return ierr;
    }
    else
       GammaVec = (xf_Vector *) GammaDat->Data;

    if (U->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);
    
    // set U->StateName
    ierr = xf_Error(xf_ReAlloc2((void ***)&U->StateName, sr, xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    
    // pull in variable name
    for (i=0; i<sr; i++) strcpy(U->StateName[i], FcnName[i]);
    
    // determine simulation time
    ierr = xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time);
    if (ierr != xf_OK) return ierr;
    
    // Determine whether should interpolate IC if a function
    ierr = xf_GetKeyValueBool(All->Param->KeyValue, "InterpolateIC", &InterpolateIC);
    if (ierr != xf_OK) return ierr;
   

    //start from user defined file 
    if(InterpolateInitData == 2)  
    {
       //init Yu_datafile struct in Model
       //ierr = xf_Error(xf_Alloc( (void **) &Model->init_datafile, 1, sizeof(Yu_datafile)));
       //if (ierr != xf_OK) return ierr;
       
       ierr = xf_Error(init_data_read_from_file(&Model->init_datafile));
       if (ierr != xf_OK) return ierr;
 
    }
    else
       Model->init_datafile = NULL;

    /*--------------------------*/
    /* No IC Function specified */
    /*--------------------------*/
    
    if (!InterpolateInitData){
    /* Only need to call EqnSetICState once, since the initial state
     is constant throughout all of the domain. */
    
    // allocate u[sr]
    ierr = xf_Error(xf_Alloc((void **) &u, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    //pull in IC state
    for(i=0; i<sr; i++)
        u[i] = FcnParam[i];
   
    // fill in array data as if basis were Lagrange
    for (i=0; i<U->nArray; i++){
        ga = U->GenArray+i;
        for (j=0; j<ga->n; j++){  
            for (k=0; k<((ga->vr==NULL) ? ga->r : ga->vr[j]); k++)
            {
                ga->rValue[j][k] = u[k%sr];
            }
        } // j
    } // i
    
    xf_Release((void *) u);
    
    // allocate vector for storing  bases
    ierr = xf_Error(xf_Alloc((void **) &BasisNew, U->nArray, sizeof(enum xfe_BasisType)));
    if (ierr != xf_OK) return ierr;
    
    for (i=0; i<U->nArray; i++){
        // get an appropriate Lagrange basis
        BasisNew[i] = U->Basis[i];
        //ierr = xf_Error(xf_Basis2Lagrange(BasisNew[i], U->Basis + i));
        ierr = xf_Error(xf_Basis2UniformLagrange(BasisNew[i], U->Basis + i));
        if (ierr != xf_OK) return ierr;
    }
    
    // convert from Lagrange basis to actual Basis
    ierr = xf_Error(xf_ProjectVectorInPlace_Basis(Mesh, NULL, U, BasisNew, xfe_BasisLast));
    if (ierr != xf_OK) return ierr;
    
    xf_Release((void *) BasisNew);
    
    }
    else 
    {
        
    /*-----------------------------------------*/
    /* Least-squares projection of IC function */
    /*-----------------------------------------*/
        
    /* Here we perform least-squares projection of a specified IC
    function. */
    QuadData    = NULL;
    PhiData     = NULL;
    GeomPhiData = NULL;
    u           = NULL;
    xglob       = NULL;
    Q           = NULL;
    R           = NULL;
    P           = NULL;
        
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
        Basis = U->Basis[egrp];
            
        pOrder = -1;
        pOrder = -1;
        pnq = -1; // here so that we reallocate at least once for every element group
            
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
               
            // get latest gamma value
            Gamma = GammaVec->GenArray[egrp].rValue[elem][0]; 

            // get interpolation order
            Order = xf_InterpOrder(U, egrp, elem);
                
            if (Order != pOrder){
                pOrder = Order;
                    
                // determine nn = # unknowns for elements in this group
                ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
                if (ierr != xf_OK) return ierr;
                    
                // determine quadrature order (not using eqnset)
                ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*Order, &QuadOrder));
                if (ierr != xf_OK) return ierr;
            }
                
            /* Pull off quad points for the element; make sure have enough
            quad points (more than nn). */
            ierr = xf_Error(xf_QuadElemAtLeast(Mesh, egrp, elem, QuadOrder, nn, &QuadData, &QuadChanged));
            if (ierr != xf_OK) return ierr;
            nq = QuadData->nquad;
            xq = QuadData->xquad;
                
            // re-allocate data if quad points increased
            if (nq > pnq){
                ierr = xf_Error(xf_ReAlloc( (void **) &u, nq*sr, sizeof(real)));
                if (ierr != xf_OK) return ierr;
                ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
                if (ierr != xf_OK) return ierr;
            }
                
            // obtain global coords of quad points
            ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
                                            nq, xq, xglob));
            if (ierr != xf_OK) return ierr;

            // Call eqnset function to determine u
            //ierr = xf_Error(xf_EqnSetICState(EqnSet, IC, IParam, RParam, nq, xglob, &Time, u));
            //if (ierr != xf_OK) return ierr;
            ierr = xf_Error(InitializeDataInterpolation(Model, nq, sr, dim, xglob, u, Gamma));
            if (ierr != xf_OK) return ierr;
                
            // Project u to get coefficients EU
            ierr = xf_Error(xf_ProjectOnElemQR(Mesh, egrp, elem, sr, QuadData, 
                                               Basis, Order, QuadChanged, &Q, &R, &P,
                                               u, U->GenArray[egrp].rValue[elem]));
            if (ierr != xf_OK) return ierr;
                
            pnq = nq;
            
        } // elem
    } // egrp
        
    // Only destroy QuadData if points are generic
    ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
    if (ierr != xf_OK) return ierr;
        
    /* Destroy Geometry Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
        
    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
        
    xf_Release((void *) u);
    xf_Release((void *) xglob);
    xf_Release((void *) Q);
    xf_Release((void *) R);
    xf_Release((void *) P);

    }  

    //free data from file
    if(InterpolateInitData == 2 && Model->init_datafile != NULL)
    {
       if(Model->init_datafile->data != NULL)
       xf_Release((void *) Model->init_datafile->data);
       
       xf_Release((void *) Model->init_datafile);
    }

    return xf_OK;
}





#if( UNIT_TEST==1 )
#include "xf_SolverTools.test.in"
#endif


