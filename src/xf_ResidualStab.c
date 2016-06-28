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
  FILE:  xf_ResidualStab.c

  This file contains stabilization residual-calculation functions.

*/
#include <stdlib.h>

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Quad.h"
#include "xf_Basis.h"
#include "xf_MeshTools.h"
#include "xf_EqnSetHook.h"
#include "xf_Line.h"
#include "xf_Residual.h"
#include "xf_ResidualDiff.h"


/******************************************************************/
//   FUNCTION Definition: xf_InitStabData
void
xf_InitStabData( xf_StabData *StabData)
{
  StabData->StabType    = xfe_StabilizationNone;
  StabData->Reg         = NULL;
  StabData->Reg_U       = NULL;

  StabData->Jump         = NULL;
  StabData->Jump_UL      = NULL;
  StabData->Jump_UR      = NULL;
}



/******************************************************************/
//   FUNCTION Definition: xf_ResolutionIndicator
static int 
xf_ResolutionIndicator(const real *s, const real *s1, const real *wq, 
		       const real *Phi, const real *Phi1, int Order, int nq, 
		       int sr, int nn, int nn1, const real *s_u, 
		       const real *s1_u1, const real *TT, 
		       enum xfe_StabSwitchType StabSwitch, real StabSwitchValue, 
		       real StabSwitchFactor, real *EReg, real *EReg_U)
{
/* 
PURPOSE

  Computes regularity estimate for an interpolated scalar quantity s,
  provided at Order and Order-1 (s1 is at Order-1).  The values are
  provided at quadrature points, and wq are the weights.  The
  regularity of the scalar is estimated by determining the ratio of
  the "energy" in s-s1 vs the energy in s1.  The energy is measured as
  an integral of the square of the quantity of interest.  This ratio
  is converted to a log scale and clipped with a smooth sine variation
  to a hardcoded range.  A linearization is computed if EReg_U is not
  NULL.  The Order to Order-1 restriction matrix, TT, must be provided
  in this case.

INPUTS:

  s    : scalar at Order; given at nq points
  s1   : scalar at Order-1; given at nq points
  wq   : quadrature weights
  Phi  : nn basis functions (Order) at nq points
  Phi1 : nn1 basis functions (Order-1) at nq points
  Order : scalar interpolation order
  nq   : number of quadrature points
  sr   : state rank
  nn   : number of basis functions at Order
  nn1  : number of basis functions at Order-1
  Need_Grad : true if gradients
  s_u  : linearization of s w.r.t state; only provided if linearization
         is requested.  Size is nq x sr
  s1_u1: linearization of s1 w.r.t state; only provided if linearization
         is requested.  Size is nq x sr
  TT   : restriction matrix from Order to Order-1
  StabSwitch      : type of switch to use to convert to indicator
  StabSwitchValue : const stab switch requires a value
  StabSwitchFactor : multiplicative factor for stabswitch

OUTPUTS:

  EReg : regularity estimate (between 0 and 1)
  EReg_U : linearization of regularity estimate (optional)

RETURNS:  Error code

*/
  int ierr, iq, k;
  enum xfe_Bool Need_Grad;
  real Num, Den, S, S0, DS, val;
  real *Num_U = NULL, *Den_U = NULL, *S_U = NULL, *rtemp = NULL;
  real sNum, sDen;
  real f, f0;

  Need_Grad = (EReg_U != NULL);

  // Allocate memory
  if (Need_Grad){
    ierr = xf_Error(xf_Alloc( (void **) &Num_U, nn*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &Den_U, nn*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &S_U  , nn*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &rtemp, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  

  // Numerator: Num = wq{q}*2*(s{q}-s1{q})*(s{q}-s1{q})
  for (iq=0, Num=0.; iq<nq; iq++) Num += wq[iq]*(s[iq]-s1[iq])*(s[iq]-s1[iq]);
  sNum = ((Num < 0.0) ? -1.0 : 1.0); // negative quad weights can make Num < 0
  Num = Num * sNum; // need Num positive

  if (Need_Grad){
    /* 
       Num_U{n,k} =  wq{q}*2*(s{q}-s1{q})*s_u{q,k}*Phi{q,n}
       -wq{q}*2*(s{q}-s1{q})*s1_u1{q,k}*Phi1{q,n1}*TT{n1,n}
    */

    // rtemp{q,k} = wq{q}*2*s{q}*s_u{q,k}
    xf_2ColcMult_Set(s_u, s,  wq, nq, sr, 1, 1,  2.0*sNum, rtemp);
    // rtemp{q,k} += -wq{q}*2*s1{q}*s_u{q,k}
    xf_2ColcMult_Add(s_u, s1, wq, nq, sr, 1, 1, -2.0*sNum, rtemp);
    // Num_U += Phi{q,n}*rtemp{q,k}
    xf_MTxM_Set(Phi, rtemp, nn, nq, sr, Num_U);

    // rtemp{q,k} = wq{q}*2*s{q}*s1_u1{q,k}
    xf_2ColcMult_Set(s1_u1, s,  wq, nq, sr, 1, 1,  2.0*sNum, rtemp);
    // rtemp{q,k} += -wq{q}*2*s1{q}*s1_u1{q,k}
    xf_2ColcMult_Add(s1_u1, s1, wq, nq, sr, 1, 1, -2.0*sNum, rtemp);
    // S_U{n1,k} = Phi1{q,n1}*rtemp{q,k}
    xf_MTxM_Set(Phi1, rtemp, nn1, nq, sr, S_U);
    // Num_U{n,k} -= TT{n1,n}*S_U{n1,k}
    xf_MTxM_Sub(TT, S_U, nn, nn1, sr, Num_U);
	
  }

  // Denominator: Den = wq{q}*s{q}*s{q}
  for (iq=0, Den=0.; iq<nq; iq++) Den += wq[iq]*s[iq]*s[iq];
  sDen = ((Den < 0.0) ? -1.0 : 1.0);
  Den = Den * sDen; // need Den positive

  if (Need_Grad){
    /* 
       Den_U{n,k} = wq{q}*2*s{q}*s_u{q,k}*Phi{q,n}
    */
    // rtemp{q,k} = wq{q}*2*s{q}*s_u{q,k}
    xf_2ColcMult_Set(s_u, s, wq, nq, sr, 1, 1, 2.0*sDen, rtemp);
    // Den_U{n,k} = Phi{q,n}*rtemp{q,k}
    xf_MTxM_Set(Phi, rtemp, nn, nq, sr, Den_U);
  }


  /* Transform Num/Den into an indicator value using a nonlinear switch */
      
  switch (StabSwitch){

  case xfe_StabSwitchConst: // A constant indicator value
    EReg[0] = StabSwitchValue;
    if (Need_Grad) for (k=0; k<nn*sr; k++) EReg_U[k] = 0.0;
    break;

  case xfe_StabSwitchLinear: // A linear power indicator
    if (Order == 1)
      f0 = .1; // .05
    else if (Order == 2)
      f0 = .04; //.02;
    else if (Order == 3)
      f0 = .03; //.015;
    else
      f0 = .01;

    EReg[0] = (Num/Den)/f0;
    if (Need_Grad)
      for (k=0; k<nn*sr; k++) EReg_U[k] = Num_U[k]/Den/f0 - EReg[0]/Den*Den_U[k];
    break;

  case xfe_StabSwitchSquare: // A square power indicator
/*     if (Order == 1) */
/*       f0 = 5e-1; // 5e-2; //.01; // .05 */
/*     else if (Order == 2) */
/*       f0 = 5e-2; // 2e-3; //2e-3; //.02; */
/*     else if (Order == 3) */
/*       f0 = 5e-3; // 3e-5; // 3e-5; //5e-5; //.00005; //.015; */
/*     else if (Order == 4) */
/*       f0 = 5e-4; // 5e-6; */
/*     else if (Order == 5) */
/*       f0 = 5e-5; // 1e-6; */
/*     else if (Order == 6) */
/*       f0 = 5e-6; // 5e-7; */
/*     else if (Order == 7) */
/*       f0 = 5e-7; // 1e-7; */
/*     else if (Order == 8) */
/*       f0 = 5e-8; */
/*     else */
/*       f0 = 5e-9*pow(10., 8.-(real) Order); */
/*     //f0 = 1e-6*pow(40., 4.-(real) Order); */
    
    //f0 = 0.37 * exp( -2. * ((real) Order));
    f0 = 0.5*pow(10., -(real) Order);
    

	
    f = Num/Den;
    if (Need_Grad)
      for (k=0; k<nn*sr; k++) EReg_U[k] = Num_U[k]/Den - f/Den*Den_U[k];

    EReg[0] = (f/(f+f0))*(f/f0);
    if (Need_Grad)
      for (k=0; k<nn*sr; k++) EReg_U[k] = (2.*f/(f+f0)/f0 - EReg[0]/(f+f0))*EReg_U[k];
    break;


  case xfe_StabSwitchLog: // A highly-nonlinear log indicator

    
    // Hardcoded values from G. Barter's thesis work
    S0 = -(log(10.0)*4.0 + 4.25*log(Order));
    DS = log(10.0)*0.5;

    /*
      S = log(Num/Den) = ln(Num) - ln(Den);
      S_U = Num_U/Num - Den_U/Den
    */
    if (Need_Grad) for (k=0; k<nn*sr; k++) S_U[k] = 0.0;
    if ((Num == 0.0) && (Den == 0.0))
      S = 0.0;
    else if (Num == 0.0)
      S = S0-2.0*DS;
    else if (Den == 0.0)
      S = S0+2.0*DS;
    else{
      S = log(Num) - log(Den);
      if (Need_Grad) 
	for (k=0; k<nn*sr; k++) S_U[k] = Num_U[k]/Num - Den_U[k]/Den;
    }

    /*  xf_printf("S = %.10E, S0 =%.10E, Num = %.10E, Den = %.10E\n",  */
    /* 	    S, S0, Num, Den); */
      
    /* 
       Set EReg:
       { 0.0                        S <= S0 - DS
       EReg =  { .5*(1+sin(pi/2*(S-S0)/DS)) |S-S0| < DS
       { 1.0                        S >= S0 + DS

       When EReg != 0,
       EReg_U = 0.5*cos(pi/2*(S-S0)/DS)*pi/(2*DS)*S_U
    */
    if (Need_Grad) for (k=0; k<nn*sr; k++) EReg_U[k] = 0.0;
    if (S <= (S0-DS))
      EReg[0] = 0.0;
    else if (S >= (S0+DS))
      EReg[0] = 1.0;
    else{
      EReg[0] = 0.5*(1.0+sin(0.5*PI*(S-S0)/DS));
      val = 0.5*cos(0.5*PI*(S-S0)/DS)*0.5*PI/DS;
      if (Need_Grad) for (k=0; k<nn*sr; k++) EReg_U[k] = val*S_U[k];
    }
    break;

  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  // multiply by factor
  if (StabSwitchFactor != 1.0){
    EReg[0] *= StabSwitchFactor;
    if (Need_Grad) for (k=0; k<nn*sr; k++) EReg_U[k] *= StabSwitchFactor;
  }
      

  // Release memory
  xf_Release( (void *) Num_U);
  xf_Release( (void *) Den_U);  
  xf_Release( (void *) S_U);  
  xf_Release( (void *) rtemp);

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateRegStab
static int
xf_CalculateRegStab(xf_All *All, xf_Vector *U, enum xfe_Bool Need_Grad, 
		    xf_SolverData *SolverData)
{
  int ierr, i, k, n, sr, sr2, iq, nq, pnq, dim;
  int pOrder, Order, Order1, QuadOrder, nn, nn1;
  int ResidualOrderIncrement;
  int egrp, elem, nResTerm;
  int negrp, negrphalo;
  int *IParam;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged;
  enum xfe_Verbosity Verbosity;
  enum xfe_StabSwitchType StabSwitch;
  char RegularityScalar[xf_MAXSTRLEN] = "None";
  real *RParam, *xq, *wq, *u, *u1, *s, *s1;
  real *s_u, *s1_u1;
  real *EReg, *EU, *EU1, *TT;
  real *EReg_U;
  real *lam, *lam_u, *lav_U, lav, sumw, fac;
  real StabSwitchValue, ERegMax = -1.0;
  real StabSwitchFactor;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData, *PhiData1;
  xf_Matrix *T = NULL;
  xf_JacobianData *JData;
  xf_Data *D;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  xf_ResTerms *ResTerms;
  xf_StabData *StabData;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;

  EqnSet = All->EqnSet;
  sr     = EqnSet->StateRank;
  sr2    = sr*sr;

  if (SolverData == NULL) return xf_OK; // need StabData to proceed

  StabData = &(SolverData->StabData);

  // error if no regularity estimate necessary
  if (StabData->StabType != xfe_StabilizationResolution) return xf_Error(xf_INPUT_ERROR);

  // Residual order increase
  ResidualOrderIncrement = ((SolverData == NULL) ? 0 : SolverData->ResidualOrderIncrement);
  
  // handle nonzero residual order increment by temporarily projecting the state
  if (ResidualOrderIncrement != 0){
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, U, NULL, 
                                                           xfe_BasisLast, ResidualOrderIncrement));
    if (ierr != xf_OK) return ierr;
  }

  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // determine stabilization switch type and value (only applicable for certain switches)
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "StabSwitch", 
				     xfe_StabSwitchName, (int ) xfe_StabSwitchLast, 
				     (int *) &StabSwitch));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "StabSwitchValue", &StabSwitchValue));
  if (ierr != xf_OK) return ierr;

  // determine stabilization factor (multiplicative)
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "StabSwitchFactor", &StabSwitchFactor));
  if (ierr != xf_OK) return ierr;



  // pull off name of scalar used for regularity estimate
  ierr = xf_Error(xf_GetKeyValue(EqnSet->KeyValue, "RegularityScalar", RegularityScalar));
  if (ierr != xf_OK) return ierr;
  
  // locate Regularity vector
  ierr = xf_Error(xf_FindVector(All, "StabRegVisc", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, xfe_True, &D, 
				&StabData->Reg, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True; // TEMPORARY


  // locate Regularity linearization vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "StabRegVisc_U", xfe_False, xfe_True, 
				       NULL, &StabData->Reg_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
    

  // initialize vars to NULL
  QuadData    = NULL;
  PhiData     = NULL;
  PhiData1    = NULL;
  JData       = NULL;
  EU1         = NULL; 
  EReg_U      = NULL;
  u           = NULL;
  u1          = NULL;
  wq          = NULL;
  s           = NULL;
  s1          = NULL;
  s_u         = NULL;
  s1_u1       = NULL;
  lam         = NULL;
  lam_u       = NULL;
  lav_U       = NULL;

  pnq         = -1;  // previous number of quad points

  // loop over element groups (including halo)
  negrp = Mesh->nElemGroup;
  negrphalo = ((Mesh->ParallelInfo == NULL) ? Mesh->nElemGroup : 2*Mesh->nElemGroup);

  for (egrp=0; egrp<negrphalo; egrp++){

    if ((Mesh->ParallelInfo != NULL) && (egrp == Mesh->nElemGroup)){
      // wait for end communication of halo data: state, if starting on halo
      ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
      if (ierr != xf_OK) return ierr;
    }

    // Determine Basis and Order from the state, U
    Basis = U->Basis[egrp%negrp];
    pOrder = -1;

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // get interpolation order
      Order = xf_InterpOrder(U, egrp, elem);
      
      if (Order != pOrder){
	pOrder = Order;
	
	// determine required integration order
	// 8/14/2011: not passing in EqnSet due to order being too high
	ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, max(3*Order, 0), 
					    &QuadOrder));
	if (ierr != xf_OK) return ierr;
	
	// number of nodes at Order
	ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_ReAlloc( (void **) &lav_U, nn*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;

	// locate transfer matrix from p to p-1
	if (Order > 0){
	  ierr = xf_Error(xf_FindTransferMatrix(All->DataSet, Basis, Order, Basis, Order-1, &T));
	  if (ierr != xf_OK) return ierr;
	  // pull off real array of the matrix
	  TT = T->GenArray->rValue[0];
	  
	  // number of nodes at Order-1
	  ierr = xf_Error(xf_Order2nNode(Basis, Order-1, &nn1));
	  if (ierr != xf_OK) return ierr;
	  
	  // reallocate memory for EU1
	  ierr = xf_Error(xf_ReAlloc( (void **) &EU1, nn1*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}
      }

      EReg = StabData->Reg->GenArray[egrp].rValue[elem];
      if (Need_Grad)
	EReg_U = StabData->Reg_U->GenArray[egrp].rValue[elem];

      if (Order == 0){
	EReg[0] = 0.0;
	if (EReg_U != NULL) for (k=0; k<nn*sr; k++) EReg_U[k] = 0.0;
	continue; // p-1 not valid when p == 0
      }

      Order1 = Order-1;

      /* Pull off quad points for the element; will not recalculate if
	 Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;
      
      // compute basis functions if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_EvalBasis(Basis, Order1, QuadChanged, nq, xq, xfb_Phi, &PhiData1));
      if (ierr != xf_OK) return ierr;
     
      /* Compute geometry Jacobian; if not constant, compute at quad
         points.  Note if jacobian is constant, only one Jacobian will
         be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;

      nn = PhiData->nn;  // number of interpolation nodes for p
      if (nn1 != PhiData1->nn) return xf_Error(xf_CODE_LOGIC_ERROR); // sanity check

      // re-allocate data if quad points increased
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &u1, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **)  &s, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &s1, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &lam, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	if (Need_Grad){
	  ierr = xf_Error(xf_ReAlloc( (void **)   &s_u, nq*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &s1_u1, nq*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &lam_u, nq*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}
      }


      EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]

      // restrict EU to get p-1 coefficients, EU1
      xf_MxM_Set(TT, EU, nn1, nn, sr, EU1);
	
      // interpolate state at quad points
      xf_MxM_Set(PhiData->Phi,  EU,  nq, nn, sr, u);   // p
      xf_MxM_Set(PhiData1->Phi, EU1, nq, nn1, sr, u1); // p-1 
   
      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];

      // call eqnset specific function for scalar using u and u1
      ierr = xf_Error(xf_EqnSetScalar(EqnSet, RegularityScalar, IParam, 
				      RParam, nq, u, NULL, s, s_u, NULL, NULL, 0.0));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_EqnSetScalar(EqnSet, RegularityScalar, IParam, 
				      RParam, nq, u1, NULL, s1, s1_u1, NULL, NULL, 0.0));
      if (ierr != xf_OK) return ierr;


      // Obtain resolution indicator
      ierr = xf_Error(xf_ResolutionIndicator(s, s1, wq, PhiData->Phi, PhiData1->Phi,
					     Order, nq, sr, nn, nn1, s_u, s1_u1, 
					     TT, StabSwitch, StabSwitchValue, 
					     StabSwitchFactor, EReg, EReg_U));
      if (ierr != xf_OK) return ierr;


      // Keep track of maximum for printing out
      ERegMax = max(EReg[0], ERegMax);

      /*
	Make artificial viscosity by adding lambda*1/p
      */
      ierr = xf_EqnSetMaxCharSpeed(EqnSet, nq, u, NULL, NULL, IParam, RParam, NULL, NULL, lam, lam_u);
      if (ierr != xf_OK) return ierr;
      
      // Determine average lav, lav_U
      for (iq=0, sumw=0.; iq<nq; iq++) sumw += wq[iq];
      lav = 0.0;
      fac = 1./((real) Order);
      if (Need_Grad) for (k=0; k<nn*sr; k++) lav_U[k] = 0.0;
      for (iq=0; iq<nq; iq++){
	lav += fac*wq[iq]*lam[iq]/sumw;
	if (Need_Grad){
	  for (n=0; n<nn; n++)
	    for (k=0; k<sr; k++)
	      lav_U[n*sr+k] += fac*wq[iq]*lam_u[iq*sr+k]*PhiData->Phi[iq*nn+n]/sumw;
	}
      } // iq
            
      if (Need_Grad)
	for (k=0; k<nn*sr; k++) EReg_U[k] = EReg_U[k]*lav + EReg[0]*lav_U[k];
      EReg[0] *= lav;
      
      pnq = nq;
    } // elem

  } // egrp

  // restore the state to the original order
  if (ResidualOrderIncrement != 0){
    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, U, NULL, 
                                                           xfe_BasisLast, -ResidualOrderIncrement));
    if (ierr != xf_OK) return ierr;
  }

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True));
  if (ierr != xf_OK) return ierr;

   /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) EU1);
  xf_Release( (void *) u);
  xf_Release( (void *) u1);
  xf_Release( (void *) wq);
  xf_Release( (void *) s);
  xf_Release( (void *) s1);
  xf_Release( (void *) s_u);
  xf_Release( (void *) s1_u1);
  xf_Release( (void *) lam);
  xf_Release( (void *) lam_u);
  xf_Release( (void *) lav_U);


  return xf_OK;

}




/******************************************************************/
//   FUNCTION Definition: xf_CreateElemJumpVisc
static int
xf_CreateElemJumpVisc(xf_All *All, xf_StabData *StabData)
{
  int ierr, negrp, nn, i;
  int egrp, elem, face, nface;
  int *OrderVec;
  enum xfe_BasisType *BasisVec;
  real fval;
  real *xref, *phi;
  xf_Data *D;
  xf_Face Face;
  xf_Vector *EJV;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  
  
  // allocate Q1 Basis and Order vectors
  negrp = Mesh->nElemGroup;
  ierr = xf_Error(xf_Alloc( (void **) &BasisVec, negrp, sizeof(enum xfe_BasisType)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &OrderVec, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // request Q1 Lagrange basis
  for (egrp=0; egrp<negrp; egrp++){
    ierr = xf_Error(xf_Basis2UniformLagrange(Mesh->ElemGroup[egrp].QBasis, BasisVec+egrp));
    if (ierr != xf_OK) return ierr;
    if (BasisVec[egrp] != Mesh->ElemGroup[egrp].QBasis) return xf_Error(xf_NOT_SUPPORTED);
    OrderVec[egrp] = 2;
  }
  
  ierr = xf_Error(xf_FindVector(All, "ElemJumpVisc", xfe_LinkageGlobElem, 1, 
				NULL, 0, 0, BasisVec, OrderVec, NULL, NULL, NULL, xfe_SizeReal, 
				xfe_False, xfe_True, &D, &EJV, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True;
  ierr = xf_Error(xf_SetZeroVector(EJV));
  if (ierr != xf_OK) return ierr;
 
  xf_Release( (void  *) BasisVec);
  xf_Release( (void  *) OrderVec);

  xref = NULL;
  phi  = NULL;

  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    // determine nn = # unknowns for order 2 lagrange interpolation for this group
    ierr = xf_Error(xf_Order2nNode(Mesh->ElemGroup[egrp].QBasis, 2, &nn));
    if (ierr != xf_OK) return ierr;

    // get ref space order 2 Lagrange nodes
    ierr = xf_Error(xf_LagrangeNodes(Mesh->ElemGroup[egrp].QBasis, 2, NULL, NULL, &xref));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReAlloc( (void **) &phi, nn, sizeof(real)));
    if (ierr != xf_OK) return ierr;


    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      nface = Mesh->ElemGroup[egrp].nFace[elem];
      for (face=0; face<nface; face++){
	Face = Mesh->ElemGroup[egrp].Face[elem][face];
	if (Face.Group < -1) continue;
	// value of stabilization viscosity on face
	fval = StabData->Jump->GenArray[1+Face.Group].rValue[Face.Number][0];
		
	// isoface shape functions on face
	ierr = xf_Error(xf_ShapeIsoFace(All, egrp, elem, face, nn, xref, phi));
	if (ierr != xf_OK) return ierr;

	for (i=0; i<nn; i++) EJV->GenArray[egrp].rValue[elem][i] += phi[i]*fval;

      } // face
    }
  } // egrp

  xf_Release( (void *) xref);
  xf_Release( (void *) phi);

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateJumpStab
static int
xf_CalculateJumpStab(xf_All *All, xf_Vector *U, enum xfe_Bool Need_Grad, 
		     xf_SolverData *SolverData)
{
  int ierr, i, k, sr, sr2, iq, nq, pnq, dim;
  int Order, OrderL, OrderR, QuadOrder, nL, nR, n;
  int egrp, rmax, nn;
  int egrpL, elemL, faceL;
  int egrpR, elemR, faceR;
  int ibfgrp, ibface, nibface;
  int OrientL, OrientR;
  int nIFaceRegular;
  int *IParam;
  enum xfe_BasisType BasisL, BasisR;
  enum xfe_Bool QuadChanged;
  enum xfe_Bool MeshIsParallel;
  enum xfe_Verbosity Verbosity;
  enum xfe_StabSwitchType StabSwitch;
  char RegularityScalar[xf_MAXSTRLEN] = "None";
  real ln10;
  real StabSwitchValue;
  real StabSwitchFactor;
  real *RParam, *xq, *wq, *uL, *uR, *sL, *sR;
  real *sL_uL, *sR_uR, *lamL, *lamR, *lamL_uL, *lamR_uR;
  real *wn, *xelemL, *xglob, *uR_uL, *sR_uL = NULL;
  real *UL, *UR;
  real *Jump, *Jump_UL = NULL, *Jump_UR = NULL;
  real Num, Den, sjmp, savg, sumw, lam;
  real *lam_UL, *lam_UR;
  real *Num_UL = NULL, *Den_UL = NULL;
  real *Num_UR = NULL, *Den_UR = NULL;
  real S0, DS, S, val, fac, f0, f, pw;
  real JumpMax = -1.0;
  real Time;
  xf_IFace IFace;
  xf_BFace BFace;
  xf_BasisTable *PhiTable;
  xf_QuadData *QuadData;
  xf_BasisData *PhiDataL, *PhiDataR;
  xf_BasisData *GeomPhiData;
  xf_NormalData *NData;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  xf_ResTerms *ResTerms;
  xf_StabData *StabData;
  xf_BC *BC;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;

  EqnSet = All->EqnSet;
  sr     = EqnSet->StateRank;
  sr2    = sr*sr;

  if (SolverData == NULL) return xf_OK; // need StabData to proceed

  StabData = &(SolverData->StabData);

  // error if no jump stab necessary
  if (StabData->StabType != xfe_StabilizationJump) return xf_Error(xf_INPUT_ERROR);

  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // determine stabilization switch type and value (only applicable for certain switches)
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "StabSwitch", 
				     xfe_StabSwitchName, (int ) xfe_StabSwitchLast, 
				     (int *) &StabSwitch));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "StabSwitchValue", &StabSwitchValue));
  if (ierr != xf_OK) return ierr;

  // determine stabilization factor (multiplicative)
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "StabSwitchFactor", &StabSwitchFactor));
  if (ierr != xf_OK) return ierr;

  // determine Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
  if (ierr != xf_OK) return ierr;


  // pull off name of scalar used for regularity estimate
  ierr = xf_Error(xf_GetKeyValue(EqnSet->KeyValue, "RegularityScalar", RegularityScalar));
  if (ierr != xf_OK) return ierr;

  // locate Jump viscosity vector
  ierr = xf_Error(xf_FindVector(All, "StabJumpVisc", xfe_LinkageFace, 1, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, xfe_True, NULL, 
				&StabData->Jump, NULL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetZeroVector(StabData->Jump));
  if (ierr != xf_OK) return ierr;

  // determine max size for Jump_UL, Jump_UR vectors
  for (egrp=0, rmax=0; egrp<U->nArray; egrp++){
    // number of nodes at Order
    ierr = xf_Error(xf_Order2nNode(U->Basis[egrp], U->Order[egrp], &nn));
    if (ierr != xf_OK) return ierr;
    rmax = max(rmax, nn*U->StateRank);
  }


  // locate Jump linearization vectors
  ierr = xf_Error(xf_FindVector(All, "StabJumpVisc_UL", xfe_LinkageFace, rmax, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, xfe_True, NULL, 
				&StabData->Jump_UL, NULL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetZeroVector(StabData->Jump_UL));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_FindVector(All, "StabJumpVisc_UR", xfe_LinkageFace, rmax, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, xfe_True, NULL, 
				&StabData->Jump_UR, NULL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetZeroVector(StabData->Jump_UR));
  if (ierr != xf_OK) return ierr;

  


  // Sort EqnSet->BCs to match boundary face groups
  ierr = xf_Error(xf_SortEqnSetBCs(Mesh, EqnSet->BCs+0));
  if (ierr != xf_OK) return ierr;
  BC = EqnSet->BCs[0].BC;


  // Allocate memory for jump calculations
  if (Need_Grad){
    ierr = xf_Error(xf_Alloc( (void **) &sR_uL, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &Num_UL, rmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &Den_UL, rmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &Num_UR, rmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &Den_UR, rmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &lam_UL, rmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &lam_UR, rmax, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  /* Create a basis table, PhiTable, that will store computed basis
     functions specific to each [element shape, face in element,
     orientation of face] combination, for quick lookup. */
  ierr = xf_Error(xf_CreateBasisTable(&PhiTable));
  if (ierr != xf_OK) return ierr;

  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
    
  // natural log of 10 (used below)
  ln10 = log(10.0);


  // initialize vars to NULL
  QuadData    = NULL;
  PhiDataL    = NULL;
  PhiDataR    = NULL;
  uL          = NULL;
  uR          = NULL;
  sL          = NULL;
  sR          = NULL;
  sL_uL       = NULL;
  sR_uR       = NULL;
  lamL        = NULL;
  lamR        = NULL;
  lamL_uL     = NULL;
  lamR_uR     = NULL;
  xelemL      = NULL;
  wn          = NULL;
  xglob       = NULL;
  uR_uL       = NULL;
  GeomPhiData = NULL;
  NData       = NULL;


  pnq         = -1;  // previous number of quad points

  // number of "regular" interior faces (not adjacent to other procs)
  nIFaceRegular = -1;
  if (MeshIsParallel = (Mesh->ParallelInfo != NULL))
    nIFaceRegular = Mesh->ParallelInfo->nIFaceRegular;

  // Loop over interior group and boundary face groups
  for (ibfgrp=-1; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    
    nibface = ((ibfgrp == -1) ? Mesh->nIFace : Mesh->BFaceGroup[ibfgrp].nBFace);

    // loop over faces
    for (ibface=0; ibface<nibface; ibface++){
      
      if (ibfgrp == -1){ // Interior face
	
	if (MeshIsParallel && (ibface == nIFaceRegular)){
	  // wait for end communication of halo data: state
	  ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
	  if (ierr != xf_OK) return ierr;
	}

	IFace = Mesh->IFace[ibface];
	// elements on L and R
	egrpL   = IFace.ElemGroupL;
	egrpR   = IFace.ElemGroupR;
	elemL   = IFace.ElemL;
	elemR   = IFace.ElemR;
	faceL   = IFace.FaceL;
	faceR   = IFace.FaceR;
	OrientL = IFace.OrientL;
	OrientR = IFace.OrientR;
      }
      else{             // Boundary face
	BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
	// element on R = elem on L (for minimal changes below)
	egrpL   = egrpR   = BFace.ElemGroup;
	elemL   = elemR   = BFace.Elem;
	faceL   = faceR   = BFace.Face;
	OrientL = OrientR = BFace.Orient;
      }

      BasisL = U->Basis[egrpL];
      BasisR = U->Basis[egrpR];
      OrderL = xf_InterpOrder(U, egrpL, elemL);
      OrderR = xf_InterpOrder(U, egrpR, elemR);      

      /* Pointers to stabilization vectors.  Note, Jump_UR is not used
	 for boundary faces, but the memory is still accessed (and
	 altered) below for minimal conditional checks. */
      Jump    = StabData->Jump->GenArray[ibfgrp+1].rValue[ibface];
      if (Need_Grad){
	Jump_UL = StabData->Jump_UL->GenArray[ibfgrp+1].rValue[ibface];
	Jump_UR = StabData->Jump_UR->GenArray[ibfgrp+1].rValue[ibface];
      }
    
    
      // determine required integration order
      if (ibfgrp == -1){
	ierr = xf_Error(xf_GetQuadOrderIFace(Mesh, EqnSet, IFace, max(OrderL,OrderR), 
					     &QuadOrder));
	if (ierr != xf_OK) return ierr;
      }
      else{
	ierr = xf_Error(xf_GetQuadOrderBFace(Mesh, EqnSet, BFace, OrderL, &QuadOrder));
	if (ierr != xf_OK) return ierr;
      }
    
      /* Pull off quad points for the iface; will not recalculate if
	 Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadFace(Mesh, egrpL, elemL, faceL, 
				  QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;
      wq = QuadData->wquad;


      // compute basis functions if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrpL, elemL, faceL, OrientL,
						   BasisL, OrderL, QuadChanged, nq, xq, 
						   xfb_Phi, &PhiDataL, PhiTable, &xelemL));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrpR, elemR, faceR, OrientR,
						   BasisR, OrderR, QuadChanged, nq, xq, 
						   xfb_Phi, &PhiDataR, PhiTable, NULL));
      if (ierr != xf_OK) return ierr;
      
      nL = PhiDataL->nn;
      nR = PhiDataR->nn;
      
      UL = U->GenArray[egrpL].rValue[elemL]; // U on elemL [nL*sr]
      UR = U->GenArray[egrpR].rValue[elemR]; // U on elemR [nR*sr]
      

      // re-allocate data if quad points increased
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc( (void **) &wn, nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &uL, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &uR, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &sL, nq  , sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &sR, nq  , sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &lamL, nq  , sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &lamR, nq  , sizeof(real)));
	if (ierr != xf_OK) return ierr;
	
	if (Need_Grad){
	  ierr = xf_Error(xf_ReAlloc( (void **) &sL_uL, nq*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &sR_uR, nq*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &lamL_uL, nq*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &lamR_uR, nq*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &uR_uL, nq*sr2, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}
      }


      // interpolate state 
      xf_MxM_Set(PhiDataL->Phi, UL, nq, nL, sr, uL);
      xf_MxM_Set(PhiDataR->Phi, UR, nq, nR, sr, uR);
      
      // on boundaries determine boundary state -> uR
      if (ibfgrp >= 0){
	/* Compute normal(s) at quad points. */
	ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, nq, xq, &NData, wn));
	if (ierr != xf_OK) return ierr;

	// obtain global coords of quad points
	ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrpL, elemL, &GeomPhiData, 
					xfe_True, nq, xelemL, xglob));
	if (ierr != xf_OK) return ierr;

	// we do not support mesh motion for the jump indicator
	if ((Mesh->Motion != NULL) && (Mesh->Motion->Active)) 
	  return xf_Error(xf_NOT_SUPPORTED);

	// determine boundary state -> uR and derivative w.r.t uL -> uR_uL
	ierr = xf_Error(xf_EqnSetBCState(EqnSet, BC+ibfgrp, IParam, RParam,
					 nq, wn, xglob, &Time, NULL, uL, uR, uR_uL));
	if (ierr != xf_OK) return ierr;
      }
      
    
      // call eqnset specific function for scalar using uL and uR
      ierr = xf_Error(xf_EqnSetScalar(EqnSet, RegularityScalar, IParam, 
				      RParam, nq, uL, NULL, sL, sL_uL, NULL, NULL, 0.0));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_EqnSetScalar(EqnSet, RegularityScalar, IParam, 
				      RParam, nq, uR, NULL, sR, sR_uR, NULL, NULL, 0.0));
      if (ierr != xf_OK) return ierr;
      
      /* 
	 Compute jump regularity estimate:
	 
	 Num   = sum_q  sq{q}*(sL{q}-sR{q})^2 / [.25*sq{q}*(sL{q}+sR{q})^2]

	 Num_UL{n,k} =  sum_q blah/blah
	 
      */

      // zero out quantities before adding
      Num = 0.;
      Den = 0.;
      if (Need_Grad)
	for (i=0; i<rmax; i++) Num_UL[i] = Num_UR[i] = Den_UL[i] = Den_UR[i] = 0.;
      
      // loop over quad points and add
      for (iq=0; iq<nq; iq++){
	sjmp =     (sL[iq]-sR[iq]);
	savg = 0.5*(sL[iq]+sR[iq]);

	f = wq[iq]*sjmp*sjmp/(savg*savg);
	Num += f;

	if (Need_Grad){
	  for (n=0; n<nL; n++)
	    for (k=0; k<sr; k++){
	      Num_UL[n*sr+k] +=   2.*wq[iq]*sjmp*sL_uL[iq*sr+k]*PhiDataL->Phi[iq*nL+n]/(savg*savg)
		                - 2.*f/savg*0.5*sL_uL[iq*sr+k]*PhiDataL->Phi[iq*nL+n];
	    }
	  if (ibfgrp == -1){ // for interior faces
	    for (n=0; n<nR; n++)
	      for (k=0; k<sr; k++){
		Num_UR[n*sr+k] += - 2.*wq[iq]*sjmp*sR_uR[iq*sr+k]*PhiDataR->Phi[iq*nR+n]/(savg*savg)
                                  - 2.*f/savg*0.5*sR_uR[iq*sr+k]*PhiDataR->Phi[iq*nR+n];
	      }
	  }
	  else{  // for boundary faces
	    xf_MTxV(uR_uL+iq*sr2, sR_uR+iq*sr, sr, sr, xfe_Set, sR_uL); // chain rule
	    for (n=0; n<nL; n++)
	      for (k=0; k<sr; k++){
		Num_UL[n*sr+k] += - 2.*wq[iq]*sjmp*sR_uL[k]*PhiDataL->Phi[iq*nL+n]/(savg*savg)
                                  - 2.*f/savg*0.5*sR_uL[k]*PhiDataL->Phi[iq*nL+n];
	      }
	  }
	}
      } // iq

      // Num should not be negative, unless negative quad weights really messed things up
      if (Num < 0.0) return xf_Error(xf_OUT_OF_BOUNDS); 

      Jump[0] = sqrt(Num);
      if (Need_Grad){
	if (Jump[0] > 0.0){
	  for (k=0; k<nL*sr; k++) Jump_UL[k] = .5/sqrt(Num)*Num_UL[k];
	  for (k=0; k<nR*sr; k++) Jump_UR[k] = .5/sqrt(Num)*Num_UR[k];
	}
	else{
	  for (k=0; k<nL*sr; k++) Jump_UL[k] = 0.0;
	  for (k=0; k<nR*sr; k++) Jump_UR[k] = 0.0;
	}
      }



      /** Below is a different formulation that is not currently used **/

      /* 	 Compute jump regularity estimate: */
	 
      /* 	 Num   = sum_q      sq{q}*(sL{q}-sR{q})^2 */
      /* 	 Den   = sum_q  .25*sq{q}*(sL{q}+sR{q})^2 */
      
      /* 	 Num_UL{n,k} = sum_q   2*(sL{q}-sR{q})*sL_uL{q,k}*PhiL{q,n}; */
      /* 	 Den_UL{n,k} = sum_q  .5*(sL{q}+sR{q})*sL_uL{q,k}*PhiL{q,n}; */
	 
      /*    *\/ */

      /*    // zero out quantities before adding */
      /*    Num = 0.; */
      /*    Den = 0.; */
      /*    if (Need_Grad) */
      /* 	for (i=0; i<rmax; i++) Num_UL[i] = Num_UR[i] = Den_UL[i] = Den_UR[i] = 0.; */
      
      /*    // loop over quad points and add */
      /*    for (iq=0; iq<nq; iq++){ */
      /* 	sjmp =     (sL[iq]-sR[iq]); */
      /* 	savg = 0.5*(sL[iq]+sR[iq]); */
      /* 	Num += wq[iq]*sjmp*sjmp; */
      /* 	Den += wq[iq]*savg*savg; */
      /* 	if (Need_Grad){ */
      /* 	  for (n=0; n<nL; n++) */
      /* 	    for (k=0; k<sr; k++){ */
      /* 	      Num_UL[n*sr+k] += 2*wq[iq]*sjmp*sL_uL[iq*sr+k]*PhiDataL->Phi[iq*nL+n]; */
      /* 	      Den_UL[n*sr+k] +=   wq[iq]*savg*sL_uL[iq*sr+k]*PhiDataL->Phi[iq*nL+n]; */
      /* 	    } */
      /* 	  if (ibfgrp == -1){ // for interior faces */
      /* 	    for (n=0; n<nR; n++) */
      /* 	      for (k=0; k<sr; k++){ */
      /* 		Num_UR[n*sr+k] -= 2*wq[iq]*sjmp*sR_uR[iq*sr+k]*PhiDataR->Phi[iq*nR+n]; */
      /* 		Den_UR[n*sr+k] +=   wq[iq]*savg*sR_uR[iq*sr+k]*PhiDataR->Phi[iq*nR+n]; */
      /* 	      } */
      /* 	  } */
      /* 	  else{  // for boundary faces */
      /* 	    xf_MTxV(uR_uL+iq*sr2, sR_uR+iq*sr, sr, sr, xfe_Set, sR_uL); // chain rule */
      /* 	    for (n=0; n<nL; n++) */
      /* 	      for (k=0; k<sr; k++){ */
      /* 		Num_UL[n*sr+k] -= 2*wq[iq]*sjmp*sR_uL[k]*PhiDataL->Phi[iq*nL+n]; */
      /* 		Den_UL[n*sr+k] +=   wq[iq]*savg*sR_uL[k]*PhiDataL->Phi[iq*nL+n]; */
      /* 	      } */
      /* 	  } */
      /* 	} */
      /*    } // iq */

      /*    // TEMPORARY */
      /*    if (Need_Grad) for (k=0; k<rmax; k++) Jump_UL[k] = Jump_UR[k] = 0.0; */
      /*    Jump[0] = sqrt(Num/Den); */
      /*    if (Need_Grad){ */
      /* 	//for (k=0; k<nL*sr; k++) Jump_UL[k] = (Num_UL[k]/Den - Jump[0]*Den_UL[k]/Den); */
      /* 	//for (k=0; k<nR*sr; k++) Jump_UR[k] = (Num_UR[k]/Den - Jump[0]*Den_UR[k]/Den); */
      /* 	if (Jump[0] > MEPS){ */
      /* 	  for (k=0; k<nL*sr; k++) Jump_UL[k] = .5*(Num_UL[k]/(Den*Jump[0]) - Jump[0]*Den_UL[k]/Den); */
      /* 	  for (k=0; k<nR*sr; k++) Jump_UR[k] = .5*(Num_UR[k]/(Den*Jump[0]) - Jump[0]*Den_UR[k]/Den); */
      /* 	} */
      /*    } */


      // Interpolation order
      Order = max(min(OrderL, OrderR), 1);


      /* Transform Jump into an indicator value using a nonlinear switch */
      
      switch (StabSwitch){

      case xfe_StabSwitchConst: // A constant indicator value
	Jump[0] = StabSwitchValue;
	if (Need_Grad){
	  for (k=0; k<nL*sr; k++) Jump_UL[k] = 0.0;
	  for (k=0; k<nR*sr; k++) Jump_UR[k] = 0.0;
	}
	break;

      case xfe_StabSwitchLinear: // A linear power indicator
	if (Order == 1)
	  f0 = .1; // .05
	else if (Order == 2)
	  f0 = .04; //.02;
	else if (Order == 3)
	  f0 = .03; //.015;
	else
	  f0 = .01;

	Jump[0] = Jump[0]/f0;
	if (Need_Grad){
	  for (k=0; k<nL*sr; k++) Jump_UL[k] = Jump_UL[k]/f0;
	  for (k=0; k<nR*sr; k++) Jump_UR[k] = Jump_UR[k]/f0;
	}
	break;


      case xfe_StabSwitchSquare: // A square power indicator
	if (Order == 1)
	  f0 = .1; // .05
	else if (Order == 2)
	  f0 = .04; //.02;
	else if (Order == 3)
	  f0 = .03; //.015;
	else
	  f0 = .01;
	
	f = Jump[0];
	Jump[0] = f*f/(f+f0)/f0;
	if ((Need_Grad) && (f > 0.0)){
	  for (k=0; k<nL*sr; k++) Jump_UL[k] = Jump[0]*(2./f - 1./(f+f0))*Jump_UL[k];
	  for (k=0; k<nR*sr; k++) Jump_UR[k] = Jump[0]*(2./f - 1./(f+f0))*Jump_UR[k];
	}
	break;


      case xfe_StabSwitchLog: // A highly-nonlinear log indicator

	// Hardcoded values from G. Barter's thesis work:
	S0 = -(2.25 + 3.0*log10( (real) Order));
	DS = 0.5;
      
	/*
	  S = log(Jump[0])
	  S_U = Jump_U/Jump[0]
	*/
	if (Need_Grad) for (k=0; k<rmax; k++) Jump_UL[k] = Jump_UR[k] = 0.0;
	if (Jump[0] <= 0.0)
	  S = 0.0;
	else{
	  S = log10(Jump[0]);
	  if (Need_Grad){
	    for (k=0; k<nL*sr; k++) Jump_UL[k] = (Jump_UL[k]/Jump[0])/ln10;
	    for (k=0; k<nR*sr; k++) Jump_UR[k] = (Jump_UR[k]/Jump[0])/ln10;
	  }
	}

	/*
	  Set Jump:
	          { 0.0                        S <= S0 - DS
	  Jump =  { .5*(1+sin(pi/2*(S-S0)/DS)) |S-S0| < DS
	          { 1.0                        S >= S0 + DS
	  
	  When Jump != 0,1,
	  Jump_U = 0.5*cos(pi/2*(S-S0)/DS)*pi/(2*DS)*S_U
	*/
	
	if (S <= (S0-DS)){
	  Jump[0] = 0.0;
	  if (Need_Grad) for (k=0; k<rmax; k++) Jump_UL[k] = Jump_UR[k] = 0.0;
	}
	else if (S >= (S0+DS)){
	  Jump[0] = 1.0;
	  if (Need_Grad) for (k=0; k<rmax; k++) Jump_UL[k] = Jump_UR[k] = 0.0;
	}
	else{
	  Jump[0] = 0.5*(1.0+sin(0.5*PI*(S-S0)/DS));
	  val = 0.5*cos(0.5*PI*(S-S0)/DS)*0.5*PI/DS;
	  if (Need_Grad){
	    for (k=0; k<nL*sr; k++) Jump_UL[k] *= val;
	    for (k=0; k<nR*sr; k++) Jump_UR[k] *= val;
	  }
	}

	break;

      default:
	return xf_Error(xf_NOT_SUPPORTED);
	break;
      }
      
      // multiply by factor
      if (StabSwitchFactor != 1.0){
	Jump[0] *= StabSwitchFactor;
	if (Need_Grad) for (k=0; k<nL*sr; k++) Jump_UL[k] *= StabSwitchFactor;
	if (Need_Grad) for (k=0; k<nR*sr; k++) Jump_UR[k] *= StabSwitchFactor;
      }


      // Keep track of maximum for printing out
      JumpMax = max(Jump[0], JumpMax);


    
      /*
	Make artificial viscosity per length by adding lambda/p
      */
      ierr = xf_EqnSetMaxCharSpeed(EqnSet, nq, uL, NULL, NULL, IParam, RParam, NULL, NULL, lamL, lamL_uL);
      if (ierr != xf_OK) return ierr;
      ierr = xf_EqnSetMaxCharSpeed(EqnSet, nq, uR, NULL, NULL, IParam, RParam, NULL, NULL, lamR, lamR_uR);
      if (ierr != xf_OK) return ierr;
      
      // Determine lam, lam_UL, and lam_UR
      for (iq=0, sumw=0.; iq<nq; iq++) sumw += wq[iq];
      lam = 0.0;
      fac = 1.0/((real) Order);

      if (Need_Grad) for (k=0; k<rmax; k++) lam_UL[k] = lam_UR[k] = 0.0;
      for (iq=0; iq<nq; iq++){
	lam += 0.5*fac*wq[iq]*(lamL[iq] + lamR[iq])/sumw;
	if (Need_Grad){
	  for (n=0; n<nL; n++)
	    for (k=0; k<sr; k++)
	      lam_UL[n*sr+k] += 0.5*fac*wq[iq]*lamL_uL[iq*sr+k]*PhiDataL->Phi[iq*nL+n]/sumw;

	  if (ibfgrp == -1){ // for interior faces
	    for (n=0; n<nR; n++)
	      for (k=0; k<sr; k++)
		lam_UR[n*sr+k] += 0.5*fac*wq[iq]*lamR_uR[iq*sr+k]*PhiDataR->Phi[iq*nR+n]/sumw;
	  }
	  else{  // for boundary faces
	    xf_MTxV(uR_uL+iq*sr2, lamR_uR+iq*sr, sr, sr, xfe_Set, sR_uL); // chain rule
	    for (n=0; n<nL; n++)
	      for (k=0; k<sr; k++){
		lam_UL[n*sr+k] += 0.5*fac*wq[iq]*sR_uL[k]*PhiDataL->Phi[iq*nL+n]/sumw;
	      }
	  }
	}
      } // iq

    
      if (Need_Grad){
	for (k=0; k<nL*sr; k++) Jump_UL[k] = Jump_UL[k]*lam + Jump[0]*lam_UL[k];
	for (k=0; k<nR*sr; k++) Jump_UR[k] = Jump_UR[k]*lam + Jump[0]*lam_UR[k];
      }
      Jump[0] *= lam;

      pnq = nq;
    } // ibface
  } // ibfgrp


  // Print out maximum jump value
  if (Verbosity != xfe_VerbosityLow)
    xf_printf("JumpMax = %.10E\n", JumpMax);

  
  /* Interpolate jump viscosities to elements for visualization */
  ierr = xf_Error(xf_CreateElemJumpVisc(All, StabData));
  if (ierr != xf_OK) return ierr;



  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiDataL, xfe_False));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(PhiDataR, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Normal Data */
  ierr = xf_Error(xf_DestroyNormalData(NData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Table */
  ierr = xf_Error(xf_DestroyBasisTable(PhiTable));
  if (ierr != xf_OK) return ierr;


  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) uL);
  xf_Release( (void *) uR);
  xf_Release( (void *) sL);
  xf_Release( (void *) sR);
  xf_Release( (void *) sL_uL);
  xf_Release( (void *) sR_uR);
  xf_Release( (void *) lamL_uL);
  xf_Release( (void *) lamR_uR);
  xf_Release( (void *) lamL);
  xf_Release( (void *) lamR);
  xf_Release( (void *) sR_uL);
  xf_Release( (void *) Num_UL);
  xf_Release( (void *) Den_UL);
  xf_Release( (void *) Num_UR);
  xf_Release( (void *) Den_UR);
  xf_Release( (void *) lam_UL);
  xf_Release( (void *) lam_UR);
  xf_Release( (void *) wn);
  xf_Release( (void *) uR_uL);
  xf_Release( (void *) xelemL);
  xf_Release( (void *) xglob);

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateStabilization
int
xf_CalculateStabilization(xf_All *All, xf_Vector *U, enum xfe_Bool Need_Grad, 
			  xf_SolverData *SolverData)
{
  int ierr, i;
  int nResTerm;
  char Value[xf_MAXSTRLEN];
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  xf_ResTerms *ResTerms;
  xf_StabData *StabData;

  
  Mesh = All->Mesh;

  EqnSet = All->EqnSet;

  if (SolverData == NULL) return xf_OK; // need StabData to proceed

  StabData = &(SolverData->StabData);

  // initialize flags
  SolverData->StabRequired = xfe_False;
  StabData->StabType  = xfe_StabilizationNone;
  StabData->Reg  = NULL;      
  StabData->Jump = NULL;

  // determine if regularity estimate is required
  ResTerms = EqnSet->ResTerms;
  nResTerm = ResTerms->nResTerm;
  for (i=0; i<nResTerm; i++){
    if ((ResTerms->ResTerm[i].Type == xfe_ResTermDiff) &&
	(ResTerms->ResTerm[i].Active)){
      ierr = xf_GetKeyValue(ResTerms->ResTerm[i].KeyValue, "Stabilization", Value);
      if (ierr != xf_OK) continue;
      
      ierr = xf_Error(xf_Value2Enum(Value, xfe_StabilizationName, xfe_StabilizationLast, 
				    (int *) &StabData->StabType));
      if (ierr != xf_OK) return ierr;
      
      
      SolverData->StabRequired = xfe_True;
    }
  }
  
  // calculate appropriate stabilization terms
  if (SolverData->StabRequired){
    switch (StabData->StabType){
    case xfe_StabilizationResolution:
      ierr = xf_Error(xf_CalculateRegStab(All, U, Need_Grad, SolverData));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_StabilizationJump:
      ierr = xf_Error(xf_CalculateJumpStab(All, U, Need_Grad, SolverData));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break; 
    }
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateStabViscElem
int
xf_CalculateStabViscElem(xf_All *All, int egrp, int elem, int nq, real *xq, 
			 real *EM, const xf_StabData *StabData,
			 real *StabVisc, real *ResMetric, enum xfe_Bool *SkipDiffStab)
{
  int ierr, k, dim, face, nface;
  int iq;
  real *phi, fval;
  xf_Face Face;
  xf_BasisData *PhiData = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;
  
  switch (StabData->StabType){
  case xfe_StabilizationResolution:
    for (iq=0; iq<nq; iq++)
      StabVisc[iq]  = StabData->Reg->GenArray[egrp].rValue[elem][0];
    break;
  case xfe_StabilizationJump:
    for (iq=0; iq<nq; iq++) StabVisc[iq] = 0.; // zero out StabVisc

    ierr = xf_Error(xf_Alloc( (void **) &phi, nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    // loop over faces, add to StabVisc
    nface = Mesh->ElemGroup[egrp].nFace[elem];
    for (face=0; face<nface; face++){
      Face = Mesh->ElemGroup[egrp].Face[elem][face];
      if (Face.Group < -1) continue;
      // value of stabilization viscosity on face
      fval = StabData->Jump->GenArray[1+Face.Group].rValue[Face.Number][0];

      // isoface shape functions on face
      ierr = xf_Error(xf_ShapeIsoFace(All, egrp, elem, face, nq, xq, phi));
      if (ierr != xf_OK) return ierr;

      for (iq=0; iq<nq; iq++) StabVisc[iq] += phi[iq]*fval;
    } // face

    xf_Release( (void *) phi);
   
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  // can we skip diffusion addition?
  (*SkipDiffStab) = xfe_True;
  for (iq=0; iq<nq; iq++) if (StabVisc[iq] != 0.) (*SkipDiffStab) = xfe_False;

  // calculate Q1 basis functions at xq
  ierr = xf_Error(xf_EvalBasis(Mesh->ElemGroup[egrp].QBasis, 1, xfe_True, nq, xq, xfb_Phi, &PhiData));
  if (ierr != xf_OK) return ierr;

  // Interpolate ElemHMetric to fill in ResMetric
  xf_MxM_Set(PhiData->Phi, EM, nq, PhiData->nn, dim*dim, ResMetric); 

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AddLinearizationStabViscElem
int
xf_AddLinearizationStabViscElem(xf_All *All, int egrp, int elem, int sr, int nn, int nq,
				real *xq, real *StabVisc, const xf_StabData *StabData,
				real *gPhi, real *wq, real *AwStab, real *T,
				real *ER_U, real **ER_NU)
{
  int ierr, i, k, d, iq, face, nface, dim;
  enum xfe_Bool IamL;
  real *Reg_U;
  real *Jump_UL, *Jump_UR;
  real *Jump_EU, *Jump_NU;
  real fval;
  real *phi;
  xf_Face Face;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim  = Mesh->Dim;
    
  switch (StabData->StabType){
  case xfe_StabilizationResolution:
    
    if (ER_U == NULL) return xf_OK; // nothing to do if don't want self derivative

    // T{n,k} = gPhi{i,n,q}*AwStab{i,q,k}*wq{q}
    xf_MTxwM_Set(gPhi+nn*nq*0, wq, AwStab+0*nq*sr, nn, nq, sr, T);
    for (i=1; i<dim; i++)
      xf_MTxwM_Add(gPhi+nn*nq*i, wq, AwStab+i*nq*sr, nn, nq, sr, T);

    // ER_U{n,k;m,a} += T{n,k}*Reg_U{m,a}
    Reg_U = StabData->Reg_U->GenArray[egrp].rValue[elem];
      
    xf_BlockOutProd_Add(T, Reg_U, nn, sr, nn, ER_U);

    break;

  case xfe_StabilizationJump:
    
    // ER_U{n,k;m,a} += gPhi{i,n,q}*AwStab{i,q,k}*wq{q}*StabVisc_U{q,m,a}

    // Contribution to StabVisc due to a face is: StabVisc(face)*phi(face,x)

    ierr = xf_Error(xf_Alloc( (void **) &phi, nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    // loop over faces, add to ER_U
    nface = Mesh->ElemGroup[egrp].nFace[elem];
    for (face=0; face<nface; face++){
      Face = Mesh->ElemGroup[egrp].Face[elem][face];
      if (Face.Group < -1) continue;
      Jump_UL = StabData->Jump_UL->GenArray[1+Face.Group].rValue[Face.Number];
      Jump_UR = StabData->Jump_UR->GenArray[1+Face.Group].rValue[Face.Number];
      if (Face.Group == -1){ // interior face
	ierr = xf_Error(xf_IsElemOnLeft(Mesh->IFace[Face.Number], egrp, elem, &IamL));
	if (ierr != xf_OK) return ierr;
	if (IamL){
	  Jump_EU = Jump_UL;
	  Jump_NU = Jump_UR;
	}
	else{
	  Jump_EU = Jump_UR;
	  Jump_NU = Jump_UL;
	}
      }
      else{                  // boundary face
	Jump_EU = Jump_UL;
	Jump_NU = NULL;
      }

      // value of face viscosity
      fval = StabData->Jump->GenArray[1+Face.Group].rValue[Face.Number][0];

      // isoface shape functions due to face
      ierr = xf_Error(xf_ShapeIsoFace(All, egrp, elem, face, nq, xq, phi));
      if (ierr != xf_OK) return ierr;

      // multiply phi by quad weights
      for (iq=0; iq<nq; iq++) phi[iq] *= wq[iq];
      
      // T{n,k} = gPhi{i,n,q}*AwStab{i,q,k}*phi{q}
      xf_MTxwM_Set(gPhi+nn*nq*0, phi, AwStab+0*nq*sr, nn, nq, sr, T);
      for (i=1; i<dim; i++)
	xf_MTxwM_Add(gPhi+nn*nq*i, phi, AwStab+i*nq*sr, nn, nq, sr, T);
      
      // ER_U{n,k;m,a} += T{n,k}*Jump_U{m,a}

      // add to ER_U
      if (ER_U != NULL) 
	xf_BlockOutProd_Add(T, Jump_EU, nn, sr, nn, ER_U);

      // add to ER_NU
      if ((ER_NU != NULL) && (Jump_NU != NULL))
	xf_BlockOutProd_Add(T, Jump_NU, nn, sr, nn, ER_NU[face]);

    } // face

    xf_Release( (void *) phi);

    break;

  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
 
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateStabViscIFace
int
xf_CalculateStabViscIFace(xf_All *All, int iiface, int nq, real *xelemL,
			  real *xelemR, real *EML, real *EMR,
			  const xf_StabData *StabData,
			  real *StabViscL, real *StabViscR, 
			  real **pStabViscL_UL, real **pStabViscL_UR,
			  real **pStabViscR_UL, real **pStabViscR_UR,
			  real *StabPhiL, real *StabPhiR,
			  real *ResMetricL, real *ResMetricR,
			  enum xfe_Bool *SkipDiffStab)
{
  int ierr, k, dim, iq;
  int egrpL, egrpR, elemL, elemR;
  real fval;
  xf_IFace IFace;
  xf_BasisData *PhiData = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;
  IFace = Mesh->IFace[iiface];

  // elements on L and R
  egrpL = IFace.ElemGroupL;
  egrpR = IFace.ElemGroupR;
  elemL = IFace.ElemL;
  elemR = IFace.ElemR;

  
  switch (StabData->StabType){
  case xfe_StabilizationResolution:
    
    for (iq=0; iq<nq; iq++)
      StabViscL[iq]  = StabData->Reg->GenArray[egrpL].rValue[elemL][0];
    for (iq=0; iq<nq; iq++)
      StabViscR[iq]  = StabData->Reg->GenArray[egrpR].rValue[elemR][0];

    for (iq=0; iq<nq; iq++) StabPhiL[iq] = 1.0;
    for (iq=0; iq<nq; iq++) StabPhiR[iq] = 1.0;

    (*pStabViscL_UL) = StabData->Reg_U->GenArray[egrpL].rValue[elemL];
    (*pStabViscL_UR) = NULL;
    (*pStabViscR_UL) = NULL;
    (*pStabViscR_UR) = StabData->Reg_U->GenArray[egrpR].rValue[elemR];
    break;

  case xfe_StabilizationJump:
    // stab visc is the same for both L and R

    // isoface shape functions on left face
    ierr = xf_Error(xf_ShapeIsoFace(All, egrpL, elemL, IFace.FaceL, nq, xelemL, StabPhiL));
    if (ierr != xf_OK) return ierr;
    for (iq=0; iq<nq; iq++) StabPhiR[iq] = StabPhiL[iq]; // same on right face

    // face viscosity value
    fval = StabData->Jump->GenArray[0].rValue[iiface][0];

    // interpolated viscosity values
    for (iq=0; iq<nq; iq++) StabViscL[iq] = StabPhiL[iq]*fval;
    for (iq=0; iq<nq; iq++) StabViscR[iq] = StabPhiR[iq]*fval;
    
    // derivatives are of the scalar face viscosities
    (*pStabViscL_UL) = StabData->Jump_UL->GenArray[0].rValue[iiface];
    (*pStabViscL_UR) = StabData->Jump_UR->GenArray[0].rValue[iiface];
    (*pStabViscR_UL) = (*pStabViscL_UL);
    (*pStabViscR_UR) = (*pStabViscL_UR);
    
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  // can we skip diffusion addition?
  (*SkipDiffStab) = xfe_True;
  for (iq=0; iq<nq; iq++) 
    if ((StabViscL[iq] != 0.) || (StabViscR[iq] != 0.)) (*SkipDiffStab) = xfe_False;

  // calculate Q1 basis functions at xq
  ierr = xf_Error(xf_EvalBasis(Mesh->ElemGroup[egrpL].QBasis, 1, xfe_True, nq, xelemL, xfb_Phi, &PhiData));
  if (ierr != xf_OK) return ierr;

  // Interpolate ElemHMetric to fill in ResMetric (metric is continuous, so using just L is ok)
  xf_MxM_Set(PhiData->Phi, EML, nq, PhiData->nn, dim*dim, ResMetricL); 
  for (iq=0; iq<nq*dim*dim; iq++) ResMetricR[iq] = ResMetricL[iq];

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateStabViscBFace
int
xf_CalculateStabViscBFace(xf_All *All, int ibfgrp, int ibface, int nq, real *xelem,
			  real *EM, const xf_StabData *StabData,
			  real *StabVisc, real **pStabVisc_U, 
			  real *StabPhi, real *ResMetric, enum xfe_Bool *SkipDiffStab)
{
  int ierr, k, dim, iq;
  int egrp, elem;
  real fval;
  xf_BFace BFace;
  xf_BasisData *PhiData = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;
  BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];

  // elements on L and R
  egrp = BFace.ElemGroup;
  elem = BFace.Elem;
  
  switch (StabData->StabType){
  case xfe_StabilizationResolution:
    for (iq=0; iq<nq; iq++)
      StabVisc[iq]  = StabData->Reg->GenArray[egrp].rValue[elem][0];
    for (iq=0; iq<nq; iq++) StabPhi[iq] = 1.0;
    (*pStabVisc_U) = StabData->Reg_U->GenArray[egrp].rValue[elem];
    break;

  case xfe_StabilizationJump:
    // isoface shape functions on face
    ierr = xf_Error(xf_ShapeIsoFace(All, egrp, elem, BFace.Face, nq, xelem, StabPhi));
    if (ierr != xf_OK) return ierr;

    // face viscosity value
    fval = StabData->Jump->GenArray[1+ibfgrp].rValue[ibface][0];

    // interpolated viscosity values
    for (iq=0; iq<nq; iq++) StabVisc[iq] = StabPhi[iq]*fval;

    // derivatives of scalar face viscosity
    (*pStabVisc_U) = StabData->Jump_UL->GenArray[1+ibfgrp].rValue[ibface];
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  
  // can we skip diffusion addition?
  (*SkipDiffStab) = xfe_True;
  for (iq=0; iq<nq; iq++) 
    if (StabVisc[iq] != 0.) (*SkipDiffStab) = xfe_False;

  // calculate Q1 basis functions at xq
  ierr = xf_Error(xf_EvalBasis(Mesh->ElemGroup[egrp].QBasis, 1, xfe_True, nq, xelem, xfb_Phi, &PhiData));
  if (ierr != xf_OK) return ierr;

  // Interpolate ElemHMetric to fill in ResMetric 
  xf_MxM_Set(PhiData->Phi, EM, nq, PhiData->nn, dim*dim, ResMetric); 

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_StabDiffA
int 
xf_StabDiffA(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
	     const int *IParam, const real *RParam, const real *ResMetric, 
	     int nq, const real *u, const real *w, real *Aw, real *A, 
	     real *A_uw, enum xfe_Bool *ConstA, real *mu)
{
  int i, j, dim, dim2;
  int ii, jj, ij2, iqsr, iqsr2;
  int k, sr, sr2;
  int iq;
  real res;

  dim  = EqnSet->Dim;
  dim2 = dim*dim;
  sr   = EqnSet->StateRank;
  sr2  = sr*sr;

  /*** Laplacian viscosity matrix ***/
    
  for (i=0; i<dim; i++){
    for (j=0; j<dim; j++){
      ii  = i*nq*sr;
      jj  = j*nq*sr;
      ij2 = (i*dim+j)*nq*sr2;
      
      for (iq=0; iq<nq; iq++){
	iqsr  = iq*sr;
	iqsr2 = iq*sr2;

	res = ResMetric[iq*dim2 + i*dim+j];

	if (A != NULL)
	  for (k=0; k<sr; k++)
	    A[ij2 + iqsr2 + k*(sr+1)] += res;
	if (Aw != NULL)
	  for (k=0; k<sr; k++)
	    Aw[ii + iqsr + k] += res * w[jj + iqsr + k];
      }
    }
  }
  
  if (mu != NULL)
    for (iq=0; iq<nq; iq++) mu[iq] += 1.0;

  if (ConstA != NULL) (*ConstA) = xfe_True;

  return xf_OK;
}




#if( UNIT_TEST==1 )
#include "xf_ResidualStab.test.in"
#endif
