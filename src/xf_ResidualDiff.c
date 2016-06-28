/*
  FILE:  xf_ResidualDiff.c

  This file contains the residual-calculation functions specific to
  diffusion terms.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_Math.h"
#include "xf_MeshTools.h"
#include "xf_EqnSetHook.h"
#include "xf_ResidualStab.h"
#include "xf_MeshMotionStruct.h"
#include "xf_LinQueue.h"


/******************************************************************/
//   FUNCTION Definition: xf_JumpPenaltyEta
int
xf_JumpPenaltyEta(enum xfe_DiffDiscType DiffDisc, 
		  enum xfe_BasisType Basis, int Order, int nFace, 
		  real *eta)
{	  
  /* Coefficient for stabilization */

  int ierr;
  enum xfe_ShapeType Shape;
  int MaxOrderTriIP = 5;
  int MaxOrderTetIP = 5;
  int MaxOrderQuadIP = 8;
  int MaxOrderHexIP = 8;
  // real etaTriIP[6] = {1.0, 4.0, 11.0, 14.0, 17.0, 20.0};
  real etaTriIP[6] = {1.0, 4.0, 11.0, 20.0, 30.0, 40.0};
  //real etaQuadIP[9] = {1.0, 1.5, 3.0, 4.0, 8.0, 10.0, 20.0, 24.0, 50.0};
  //real etaQuadIP[9] = {1.0, 4.0, 8.0, 16.0, 20.0, 30.0, 35.0, 45.0, 50.0};
  real etaQuadIP[9] = {1.0, 4.0, 12.0, 12.0, 20.0, 30.0, 35.0, 45.0, 50.0};
  real etaTetIP[6] = {1.0, 5.0, 20.0, 30.0, 40.0, 50.0};
  real etaHexIP[9] = {1.0, 2.5, 5.0, 8.0, 10.0, 14.0, 20.0, 24.0, 50.0};

  ierr = xf_Error(xf_Basis2Shape(Basis, &Shape));
  if (ierr != xf_OK) return ierr;

  switch (DiffDisc){
  case xfe_DiffDiscIP:
    switch (Shape){
    case xfe_Triangle: 
      (*eta) = etaTriIP[min(Order,MaxOrderTriIP)]; break;
    case xfe_Quadrilateral: 
      (*eta) = etaQuadIP[min(Order,MaxOrderQuadIP)]; break;
    case xfe_Tetrahedron: 
      (*eta) = etaTetIP[min(Order,MaxOrderTetIP)]; break;
    case xfe_Hexahedron: 
      (*eta) = etaHexIP[min(Order,MaxOrderHexIP)]; break;
    default: 
      return xf_Error(xf_NOT_SUPPORTED); break;
    }
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  (*eta) *= nFace;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateDiffJumpData
int
xf_CreateDiffJumpData(xf_DiffJumpData **pDData)
{
  int ierr;
  xf_DiffJumpData *DData;
  
  ierr = xf_Error(xf_Alloc( (void **) pDData, sizeof(xf_DiffJumpData), 1));
  if (ierr != xf_OK) return ierr;

  DData = (*pDData);

  DData->Need_Grad = xfe_True;
  DData->ConstAL   = xfe_False;
  DData->ConstAR   = xfe_False;
  DData->alpha     = 0.0;
  DData->duL_uL = DData->duL_uR = 0.0;
  DData->duR_uL = DData->duR_uR = 0.0;

  DData->dunL     = NULL; DData->dunR     = NULL;
  DData->ALguL    = NULL; DData->ARguR    = NULL;
  DData->ALdunL   = NULL; DData->ARdunR   = NULL;
  DData->AL       = NULL; DData->AR       = NULL;
  DData->A_uLguL  = NULL; DData->A_uRguR  = NULL;
  DData->A_uLdunL = NULL; DData->A_uRdunR = NULL;
  DData->N        = NULL; 
  DData->Qn       = NULL; 
  DData->Qn_u     = NULL;
  DData->Qn_gu    = NULL;
  DData->AwStabL  = NULL; DData->AwStabR  = NULL;
  DData->T        = NULL;
  DData->Aw       = NULL;
  DData->DnL      = NULL; DData->DnR      = NULL;
  DData->SLL      = NULL; DData->SLR      = NULL;
  DData->SRL      = NULL; DData->SRR      = NULL;
  DData->PL       = NULL; DData->PR       = NULL;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReAllocDiffJumpData
int
xf_ReAllocDiffJumpData(int sr, int nq, int dim, int nn, 
		       enum xfe_Bool Need_Grad, xf_DiffJumpData *DData)
{
  int ierr;
  int sr2, Tsize;

  sr2 = sr*sr;
  DData->Need_Grad = Need_Grad;

  ierr = xf_Error(xf_ReAlloc( (void **) &DData->dunL, dim*nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->dunR, dim*nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ALguL, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ARguR, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ALdunL, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ARdunR, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn, sr*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->N , nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->DnL, 2*dim*nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  DData->DnR = DData->DnL + dim*nn*sr;

  // BR2 needs these even without gradient
  Tsize = nn*nq;
  Tsize = max(Tsize, nn*sr);
  Tsize = max(Tsize, nn*nn);
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->T, Tsize, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->SLL, 4*dim*nn*nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  DData->SLR = DData->SLL +   dim*nn*nn;
  DData->SRL = DData->SLL + 2*dim*nn*nn;
  DData->SRR = DData->SLL + 3*dim*nn*nn;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->PL, 2*dim*nn*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  DData->PR = DData->PL + dim*nn*nq;
	
  if (Need_Grad){
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AL, sr2*nq*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AR, sr2*nq*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uLguL , sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uRguR , sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uLdunL, sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uRdunR, sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_u, nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_gu, dim*nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AwStabL, dim*nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AwStabR, dim*nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Aw, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDiffJumpData
void
xf_DestroyDiffJumpData(xf_DiffJumpData *DData)
{
  xf_Release( (void *) DData->dunL);     xf_Release( (void *) DData->dunR);
  xf_Release( (void *) DData->AL);       xf_Release( (void *) DData->AR);
  xf_Release( (void *) DData->ALguL);    xf_Release( (void *) DData->ARguR);
  xf_Release( (void *) DData->ALdunL);   xf_Release( (void *) DData->ARdunR);
  xf_Release( (void *) DData->A_uLguL);  xf_Release( (void *) DData->A_uRguR);
  xf_Release( (void *) DData->A_uLdunL); xf_Release( (void *) DData->A_uRdunR);
  xf_Release( (void *) DData->N);
  xf_Release( (void *) DData->Qn);
  xf_Release( (void *) DData->Qn_u);
  xf_Release( (void *) DData->Qn_gu);
  xf_Release( (void *) DData->AwStabL);  xf_Release( (void *) DData->AwStabR);
  xf_Release( (void *) DData->T);
  xf_Release( (void *) DData->Aw);
  xf_Release( (void *) DData->DnL);
  xf_Release( (void *) DData->SLL);
  xf_Release( (void *) DData->PL);
    

  xf_Release( (void *) DData);
}

/******************************************************************/
//   FUNCTION Definition: xf_DiffFluxUJump
static int
xf_DiffFluxUJump(enum xfe_DiffDiscType DiffDisc, real *alpha){
  
  switch (DiffDisc){
  case xfe_DiffDiscIP:
  case xfe_DiffDiscBR2:
    (*alpha) = 0.0;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ModMotionDiffA
static void
xf_ModMotionDiffA(int nq, int dim, int sr, const xf_MotionData *MD, real *A)
{
/*
PURPOSE: 

  Modifies diffusive flux matrix to take into account mesh motion

INPUTS: 
 
  nq     : number of points
  dim    : dimension
  sr     : state rank
  MD     : mesh motion data structure
  A      : diffusive flux linearization matrix [dim*dim*nq*sr*sr] (optional)

OUTPUTS:

  A      : set to G^{-1} * A * G^{-T} * g/gb

RETURNS: Error code
*/
  int d, iq, dim2, sr2;

  dim2 = dim*dim;
  sr2  = sr*sr;

  /* A{i,j,q,k,a} = Ginv{q,i,ii}*A{ii,jj,q,k,a}*Ginv{q,jj,j} * g{iq}/gb{iq} */
  for (d=0; d<dim; d++)
    for (iq=0; iq<nq; iq++)
      xf_dMxM(MD->Ginv+iq*dim2, 1.0, dim, sr2, dim*nq*sr2, A+iq*sr2+d*nq*sr2);
  for (d=0; d<dim; d++)
    for (iq=0; iq<nq; iq++)
      xf_dMxM(MD->Ginv+iq*dim2, MD->g[iq]/MD->gb[iq], dim, sr2, nq*sr2, A+iq*sr2+d*dim*nq*sr2);

}

/******************************************************************/
//   FUNCTION Definition: xf_ModMotionDiffAw
static void
xf_ModMotionDiffAw(int nq, int dim, int sr, const xf_MotionData *MD, real *Aw)
{
/*
PURPOSE: 

  Modifies diffusive flux matrix-vector product to take into account mesh motion

INPUTS: 
 
  nq     : number of points
  dim    : dimension
  sr     : state rank
  MD     : mesh motion data structure
  Aw     : A multiplied by some vector w [dim*nq*sr] (optional)

OUTPUTS:

  Aw     : set to g/gb * G^{-1} * Aw

RETURNS: Error code
*/
  int d, iq, dim2;

  dim2 = dim*dim;

  /* Aw{i,q,k} = Ginv{q,i,ii}*Aw{ii,q,k} * g{iq}/gb{iq} */
  for (iq=0; iq<nq; iq++) 
    xf_dMxM(MD->Ginv+iq*dim2, MD->g[iq]/MD->gb[iq], dim, sr, nq*sr, Aw+iq*sr);
  
}

/******************************************************************/
//   FUNCTION Definition: xf_ModMotionDiffA_uw
static void
xf_ModMotionDiffA_uw(int nq, int dim, int sr, enum xfe_Bool uIsRef,
                     const xf_MotionData *MD, real *A_uw)
{
/*
PURPOSE: 

  Modifies diffusive flux matrix linearization to take into account mesh motion

INPUTS: 
 
  nq     : number of points
  dim    : dimension
  sr     : state rank
  uIsRef : true if A was computed from a reference u (i.e. uX)
  MD     : mesh motion data structure
  A_uw   : linearization of A multiplied by w [dim*nq*sr*sr] (optional)

OUTPUTS:

  A_uw   : set to g/gb * gb^{-1}*G^{-1} * A_uw
                  ----   -------
		    ^       ^ from chain rule derivative transform (only if uIsRef)
		    |_ from GCL

RETURNS: Error code
*/
  int d, iq, dim2, sr2;
  real fac;
  
  dim2 = dim*dim;
  sr2  = sr*sr;

  /* A_uw{i,q,k,a} = (g{q}/gb{q})*(1/gb{q})*Ginv{q,i,ii}*A_uw{ii,q,k,a} */
  for (iq=0; iq<nq; iq++){
    fac = MD->g[iq]/MD->gb[iq];
    if (uIsRef) fac /= MD->gb[iq];
    xf_dMxM(MD->Ginv+iq*dim2, fac, dim, sr2, nq*sr2, A_uw+iq*sr2);
  }
    
}



/******************************************************************/
//   FUNCTION Definition: xf_ComputeDiffA
int 
xf_ComputeDiffA(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm,
		int nDiff, const int *IParam, const real *RParam,
		const real *ResMetric, const real *StabVisc, int nq, 
		real *u, real *gu, real *w, enum xfe_Bool uIsRef, 
		const xf_MotionData *MD, real *Aw, real *A, 
		real *A_uw, real *AwStab, enum xfe_Bool *ConstA, real *CF)
{
  int ierr, d, dim, sr, sr2, iDiff, k, iq;
  enum xfe_Bool StabFlag = xfe_False;
  enum xfe_Bool SkipStab = xfe_True;
  enum xfe_Bool ConstATemp = xfe_False;
  enum xfe_Bool ConstAReturn = xfe_False;
  enum xfe_Bool HaveViscModel = xfe_False;
  enum xfe_Bool MotionOn;
  char Value[xf_MAXSTRLEN];
  real *CFtemp = NULL;

  dim = EqnSet->Dim;
  sr  = EqnSet->StateRank;
  sr2 = sr*sr;

  if (nDiff <= 0) return xf_Error(xf_CODE_LOGIC_ERROR);

  // MD passed-in means motion is on
  MotionOn = (MD != NULL);

  // Can we skip stabilization
  SkipStab = xfe_True;
  if (StabVisc != NULL)
    for (iq=0, SkipStab = xfe_True; iq<nq; iq++) if (StabVisc[iq] != 0.) SkipStab = xfe_False;
  
  // if AwStab != NULL and Aw == NULL error (for now)
  if ((AwStab != NULL) && (Aw == NULL)) return xf_Error(xf_INPUT_ERROR);

  // Initialize ConstA
  if (ConstA != NULL) (*ConstA) = xfe_False;
  ConstAReturn = xfe_False;

  // Initialize Aw, A, etc. to zero
  if (Aw     != NULL) for (k=0; k<dim*nq*sr     ; k++) Aw[k]     = 0.0;
  if (A      != NULL) for (k=0; k<dim*dim*nq*sr2; k++) A[k]      = 0.0;
  if (A_uw   != NULL) for (k=0; k<dim*nq*sr2    ; k++) A_uw[k]   = 0.0;
  if (AwStab != NULL) for (k=0; k<dim*nq*sr     ; k++) AwStab[k] = 0.0;

  // transform state to physical if input is reference
  if (uIsRef && MotionOn) xf_ColDiv(u, MD->gb, nq, sr, 1);

  // if mesh motion is on, transform w -> G^{-T}*w,
  if (MotionOn)
    for (iq=0; iq<nq; iq++)
      xf_dMxMT(MD->Ginv+iq*dim*dim, 1.0, dim, sr, nq*sr, w+iq*sr);

  // Handle stabilization terms first
  if (!SkipStab){
    // For adding to CF, allocate a temp vector and add at end
    if (CF != NULL){
      ierr = xf_Error(xf_Alloc( (void **) &CFtemp, nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      for (k=0; k<nq; k++) CFtemp[k] = 0.0;
    }

    // Loop over diffusion terms
    for (iDiff=0; iDiff<nDiff; iDiff++){
      
      // Determine if this diffusion term is for stabilization
      ierr = xf_GetKeyValue(ResTerm[iDiff].KeyValue, "Stabilization", Value);
      StabFlag = (ierr == xf_OK);
      if (!StabFlag) continue;

      // Determine if have an equation-set specific viscosity model
      ierr = xf_GetKeyValue(ResTerm[iDiff].KeyValue, "ViscModel", Value); 
      HaveViscModel = (ierr == xf_OK);


      // Add stabilization A; if have a ViscModel, use eqn-set function
      if (HaveViscModel){
	ierr = xf_Error(xf_EqnSetDiffA(EqnSet, ResTerm+iDiff, IParam, RParam,
				       ResMetric, nq, u, gu, w, Aw, A, A_uw, 
				       &ConstATemp, CF));
	if (ierr != xf_OK) return ierr;
      }
      else{ // use default EqnSet-general Laplace viscosity
	ierr = xf_Error(xf_StabDiffA(EqnSet, ResTerm+iDiff, IParam, RParam,
				     ResMetric, nq, u, w, Aw, A, A_uw,
				     &ConstATemp, CFtemp));
	if (ierr != xf_OK) return ierr;
      }

      if (ConstA != NULL) ConstAReturn = (ConstAReturn) && ConstATemp;
    }

    // Store Aw in AwStab
    if (AwStab != NULL) for (k=0; k<dim*nq*sr; k++) AwStab[k] = Aw[k];

    // Multiply matrices by StabVisc
    if (Aw    != NULL) 
      for (d=0; d<dim    ; d++) xf_ColMult(Aw  +nq*sr *d, StabVisc, nq, sr , 1);    
    if (A     != NULL) 
      for (d=0; d<dim*dim; d++) xf_ColMult(A   +nq*sr2*d, StabVisc, nq, sr2, 1);    
    if (A_uw  != NULL) 
      for (d=0; d<dim    ; d++) xf_ColMult(A_uw+nq*sr2*d, StabVisc, nq, sr2, 1);

    // add to CF and release CFtemp
    if (CF != NULL){
      for (iq=0; iq<nq; iq++) CF[iq] += StabVisc[iq]*CFtemp[iq];
      xf_Release((void *) CFtemp);
    }
  }

  // Add non-stabilization terms next
  for (iDiff=0; iDiff<nDiff; iDiff++){
    // Determine if this diffusion term is for stabilization
    StabFlag = xfe_False;
    ierr = xf_GetKeyValue(ResTerm[iDiff].KeyValue, "Stabilization", Value);
    if (ierr == xf_OK) StabFlag = xfe_True;
    
    if (StabFlag) continue;
    
    ierr = xf_Error(xf_EqnSetDiffA(EqnSet, ResTerm+iDiff, IParam, RParam,
				   ResMetric, nq, u, gu, w, Aw, A, A_uw, 
				   &ConstATemp, CF));
    if (ierr != xf_OK) return ierr;
    
    if (ConstA != NULL) ConstAReturn = (ConstAReturn) && ConstATemp;
  }

  if (ConstA != NULL) (*ConstA) = ConstAReturn;


  // if motion on, apply mesh motion modifications to A, Aw, AwStab, A_uw
  if (MotionOn){
    if (A  != NULL) xf_ModMotionDiffA(nq, dim, sr, MD, A);
    if (Aw != NULL) xf_ModMotionDiffAw(nq, dim, sr, MD, Aw);
    if ((AwStab != NULL) && (!SkipStab)) 
      xf_ModMotionDiffAw(nq, dim, sr, MD, AwStab);
    if ((A_uw != NULL) && (!ConstAReturn))
      xf_ModMotionDiffA_uw(nq, dim, sr, uIsRef, MD, A_uw);
  }

  // if mesh motion is on, transform w back to original: w -> G^T*w,
  if (MotionOn)
    for (iq=0; iq<nq; iq++)
      xf_dMxMT(MD->G+iq*dim*dim, 1.0, dim, sr, nq*sr, w+iq*sr);

  // transform state back to reference if original input was reference
  if (uIsRef && MotionOn) xf_ColMult(u, MD->gb, nq, sr, 1);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DiffFluxQJump
int
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
		 xf_LinQueueData *LinQ, xf_DiffJumpData *DData, real *CF)
{
  int ierr, i, j, k, iq, dim, ii, sr, sr2;
  int egrpL, egrpR, elemL, elemR, nFace, Order;
  int nL, nR, n, n2;
  int ResOrderL, ResOrderR, RnL, RnR;
  enum xfe_Bool Need_Grad;
  enum xfe_Bool MotionOn;
  enum xfe_BasisType Basis;
  real ViscStabFactor;
  real *T, nval, eta, FaceArea, ihL, ihR;
  real *ResMetricL, *ResMetricR;
  real *StabViscL = NULL, *StabViscL_UL = NULL, *StabViscL_UR = NULL;
  real *StabViscR = NULL, *StabViscR_UL = NULL, *StabViscR_UR = NULL;
  real *iMML, *iMMR, facL, facR;
  real *wL = NULL, *wR = NULL;
  xf_LinQueueData *LinQLL, *LinQLR, *LinQRL, *LinQRR;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  xf_IFace IFace;

  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;  
  IFace  = Mesh->IFace[iiface];

  dim  = Mesh->Dim;
  sr   = EqnSet->StateRank;
  sr2  = sr*sr;

  nL = PhiDataL->nn;
  nR = PhiDataR->nn;

  // residual orders (possibly different from state orders)
  ResOrderL = ResPhiDataL->Order;
  ResOrderR = ResPhiDataR->Order;
  RnL = ResPhiDataL->nn;
  RnR = ResPhiDataR->nn;

  // determine stabilization factor
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "ViscStabFactor", &ViscStabFactor);
  if (ierr == xf_NOT_FOUND)  ViscStabFactor = 1.0; // backwards-compatibility
  else if (ierr != xf_OK) return xf_Error(ierr); 
  
  // element group info
  egrpL = IFace.ElemGroupL;  egrpR = IFace.ElemGroupR;
  elemL = IFace.ElemL;       elemR = IFace.ElemR;

  // pull off temporary matrices
  T = DData->T;

  // MD passed-in means motion is on
  MotionOn = ((MDL != NULL) || (MDR != NULL));

  // pull off data required for stabilization
  if (StabData != NULL){
    StabViscL    = StabData->StabViscL;
    StabViscR    = StabData->StabViscR;
    StabViscL_UL = StabData->StabViscL_UL;
    StabViscL_UR = StabData->StabViscL_UR;
    StabViscR_UL = StabData->StabViscR_UL;
    StabViscR_UR = StabData->StabViscR_UR;
    ResMetricL   = StabData->ResMetricL;
    ResMetricR   = StabData->ResMetricR;

  }
  else StabViscL = StabViscR = NULL;

  /* Point to residual linearization queues */
  LinQLL=LinQ+0;  LinQLR=LinQ+1;  LinQRL=LinQ+2;  LinQRR=LinQ+3;

  /* For a general diffusion discretization:

     Qn{q,k} = 0.5*(ALguL{i,q,k}+ARguR{i,q,k})*wn{q,i} - (visc. stab) 

     where (visc. stab) depends on the paticular form used.
     Quad weights are included through the use of wn.

     Normally (taking L as an example),
         ALguL{i,q,k} = AL{i,j,q,k,a}*guL{j,q,a};
     however, if mesh motion is on, use
         ALguL{i,q,k} = AL{i,j,q,k,a}*(guL{j,q,a} - uL{q,a}*gbigb_X{q,j})
	                              \____ use this as w in Aw  _____/
				           [store in DData->dunL]
  */

  // zero out coupling
  if (CF != NULL) for (iq=0;iq<nq; iq++) CF[iq] = 0.;

  // wL is the vector that AL multiplies
  wL = DData->dunL;  // using this space as temporary storage
  for (k=0;k<dim*nq*sr;k++) wL[k] = guL[k]; // set wL = gradient
  if (MotionOn)      // modify state gradient in presence of mesh motion
    for (j=0; j<dim; j++)
      xf_ColMult_Sub(uL, MDL->gbigb_X+j, nq, sr, dim, wL+nq*sr*j);

  // calculate ALguL [and AL, A_uLguL] at quad points
  ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm, nDiff, IParam, RParam,
				  ResMetricL, StabViscL, nq, uL, guL, wL, 
				  xfe_True, MDL, DData->ALguL, 
				  DData->AL, DData->A_uLguL, DData->AwStabL,
				  &DData->ConstAL, CF));
  if (ierr != xf_OK) return ierr;

  // wR is the vector that AR multiplies
  wR = DData->dunR;  // using this space as temporary storage
  for (k=0;k<dim*nq*sr;k++) wR[k] = guR[k]; // set wR = gradient
  if (MotionOn)      // modify state gradient in presence of mesh motion
    for (j=0; j<dim; j++)
      xf_ColMult_Sub(uR, MDR->gbigb_X+j, nq, sr, dim, wR+nq*sr*j);

  // same for R
  ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm, nDiff, IParam, RParam,
				  ResMetricR, StabViscR, nq, uR, guR, wR, 
				  xfe_True, MDR, DData->ARguR,
				  DData->AR, DData->A_uRguR, DData->AwStabR,
				  &DData->ConstAR, CF));
  if (ierr != xf_OK) return ierr;

  // average coupling
  if (CF != NULL) for (iq=0;iq<nq; iq++) CF[iq] *= 0.5;

 
  // flux term first: Qn{q,k} = 0.5*(ALguL{i,q,k}+ARguR{i,q,k})*wn{q,i}
  xf_ColcMult_Set(DData->ALguL+0*nq*sr, wn+0, nq, sr, dim, 0.5, DData->Qn);
  xf_ColcMult_Add(DData->ARguR+0*nq*sr, wn+0, nq, sr, dim, 0.5, DData->Qn);
  for (i=1; i<dim; i++){
    xf_ColcMult_Add(DData->ALguL+i*nq*sr, wn+i, nq, sr, dim, 0.5, DData->Qn);
    xf_ColcMult_Add(DData->ARguR+i*nq*sr, wn+i, nq, sr, dim, 0.5, DData->Qn);
  }

  /* Derivatives of flux term */

  /*
    Qn_guL{j,q,k;a} = 0.5*AL{i,j,q,k,a}*wn{i,q}
    Qn_uL{q,k;a}    = 0.5*A_uLguL{i,q,k,a}*wn{i,q}
    RL_UL{n,k;m,a} -= PhiL{n,q}*[Qn_uL{q,k;a}*PhiL{q,m}
	                         + Qn_guL{i,q,k;a}*gPhiL{i,q,m}
	       (if MotionOn)     - Qn_guL{i,q,k;a}*gbigb_X{q,i}*PhiL{q,m}]
    RR_UL{n,k;m,a} += PhiR{n,q}*[Qn_uL{q,k;a}*PhiL{q,m}
	                         + Qn_guL{i,q,k;a}*gPhiL{i,q,m}
	       (if MotionOn)     - Qn_guL{i,q,k;a}*gbigb_X{q,i}*PhiL{q,m}]
  */
  if ((RL_UL != NULL) || (RR_UL != NULL)){
    for (j=0; j<dim; j++){ // Qn_guL
      ii = j*nq*sr2;
      xf_ColcMult_Set(DData->AL+j*nq*sr2, wn+0, nq, sr2, dim, 0.5, DData->Qn_gu+ii);
      for (i=1; i<dim; i++)
	xf_ColcMult_Add(DData->AL+(i*dim+j)*nq*sr2, wn+i, nq, sr2, dim, 0.5, DData->Qn_gu+ii);
    }
    if (!DData->ConstAL){ // Qn_uL
      xf_ColcMult_Set(DData->A_uLguL, wn+0, nq, sr2, dim, 0.5, DData->Qn_u);
      for (i=1; i<dim; i++)
	xf_ColcMult_Add(DData->A_uLguL+i*nq*sr2, wn+i, nq, sr2, dim, 0.5, DData->Qn_u);
    }

    if (RL_UL != NULL){ // RL_UL
      // Qn_uL contribution, only if AL is not constant
      if (!DData->ConstAL){
	ierr = xf_Error(xf_AddToLinQueue(DData->Qn_u, xfe_LinQTerm_PhiPhi, nq, dim, 
					 sr2, -1, NULL, 0, -1.0, LinQLL));
	if (ierr != xf_OK) return ierr;
      }
      // Qn_guL contribution
      ierr = xf_Error(xf_AddToLinQueue(DData->Qn_gu, xfe_LinQTerm_PhiGPhi, nq, dim, 
				       sr2, -1, NULL, 0, -1.0, LinQLL));
      if (ierr != xf_OK) return ierr;
      // mesh motion contribution
      if (MotionOn){
	for (i=0; i<dim; i++){
	  ierr = xf_Error(xf_AddToLinQueue(DData->Qn_gu+i*nq*sr2, xfe_LinQTerm_PhiPhi, nq, 
					   dim, sr2, -1, MDL->gbigb_X+i, dim, 1.0, LinQLL));
	  if (ierr != xf_OK) return ierr;
	}
      }
    }
    if (RR_UL != NULL){ // RR_UL
      // Qn_uL contribution, only if AL is not constant
      if (!DData->ConstAL){
	ierr = xf_Error(xf_AddToLinQueue(DData->Qn_u, xfe_LinQTerm_PhiPhi, nq, dim, 
					 sr2, -1, NULL, 0, 1.0, LinQRL));
	if (ierr != xf_OK) return ierr;
      }
      // Qn_guL contribution
      ierr = xf_Error(xf_AddToLinQueue(DData->Qn_gu, xfe_LinQTerm_PhiGPhi, nq, dim, 
				       sr2, -1, NULL, 0, 1.0, LinQRL));
      if (ierr != xf_OK) return ierr;
      // mesh motion contribution
      if (MotionOn){
	for (i=0; i<dim; i++){
	  ierr = xf_Error(xf_AddToLinQueue(DData->Qn_gu+i*nq*sr2, xfe_LinQTerm_PhiPhi, nq, 
					   dim, sr2, -1, MDL->gbigb_X+i, dim, -1.0, LinQRL));
	  if (ierr != xf_OK) return ierr;
	}
      }
    }
  }

  /*
    Qn_guR{j,q,k;a} = 0.5*AR{i,j,q,k,a}*wn{i,q}
    Qn_uR{q,k;a}    = 0.5*A_uRguR{i,q,k,a}*wn{i,q}
    RL_UR{n,k;m,a} -= PhiL{n,q}*[Qn_uR{q,k;a}*PhiR{q,m}
  	                         + Qn_guR{i,q,k;a}*gPhiR{i,q,m}
	       (if MotionOn)     - Qn_guR{i,q,k;a}*gbigb_X{q,i}*PhiR{q,m}]
    RR_UR{n,k;m,a} += PhiR{n,q}*[Qn_uR{q,k;a}*PhiR{q,m}
 	                         + Qn_guR{i,q,k;a}*gPhiR{i,q,m}
	       (if MotionOn)     - Qn_guR{i,q,k;a}*gbigb_X{q,i}*PhiR{q,m}]
  */
  if ((RL_UR != NULL) || (RR_UR != NULL)){
    for (j=0; j<dim; j++){ // Qn_guR
      ii = j*nq*sr2;
      xf_ColcMult_Set(DData->AR+j*nq*sr2, wn+0, nq, sr2, dim, 0.5, DData->Qn_gu+ii);
      for (i=1; i<dim; i++)
	xf_ColcMult_Add(DData->AR+(i*dim+j)*nq*sr2, wn+i, nq, sr2, dim, 0.5, DData->Qn_gu+ii);
    }
    if (!DData->ConstAR){ // Qn_uR
      xf_ColcMult_Set(DData->A_uRguR, wn+0, nq, sr2, dim, 0.5, DData->Qn_u);
      for (i=1; i<dim; i++)
	xf_ColcMult_Add(DData->A_uRguR+i*nq*sr2, wn+i, nq, sr2, dim, 0.5, DData->Qn_u);
    }

    if (RL_UR != NULL){ // RL_UR
      // Qn_uR contribution, only if AR is not constant
      if (!DData->ConstAR){
	ierr = xf_Error(xf_AddToLinQueue(DData->Qn_u, xfe_LinQTerm_PhiPhi, nq, dim, 
					 sr2, -1, NULL, 0, -1.0, LinQLR));
	if (ierr != xf_OK) return ierr;
      }
      // Qn_guR contribution
      ierr = xf_Error(xf_AddToLinQueue(DData->Qn_gu, xfe_LinQTerm_PhiGPhi, nq, dim, 
				       sr2, -1, NULL, 0, -1.0, LinQLR));
      if (ierr != xf_OK) return ierr;
      // mesh motion contribution
      if (MotionOn){
	for (i=0; i<dim; i++){
	  ierr = xf_Error(xf_AddToLinQueue(DData->Qn_gu+i*nq*sr2, xfe_LinQTerm_PhiPhi, nq, 
					   dim, sr2, -1, MDR->gbigb_X+i, dim, 1.0, LinQLR));
	  if (ierr != xf_OK) return ierr;
	}
      }
    }
    if (RR_UR != NULL){ // RR_UR
      // Qn_uR contribution, only if AR is not constant
      if (!DData->ConstAR){
	ierr = xf_Error(xf_AddToLinQueue(DData->Qn_u, xfe_LinQTerm_PhiPhi, nq, dim, 
					 sr2, -1, NULL, 0, 1.0, LinQRR));
	if (ierr != xf_OK) return ierr;
      }
      // Qn_guR contribution
      ierr = xf_Error(xf_AddToLinQueue(DData->Qn_gu, xfe_LinQTerm_PhiGPhi, nq, dim, 
				       sr2, -1, NULL, 0, 1.0, LinQRR));
      if (ierr != xf_OK) return ierr;
      // mesh motion contribution
      if (MotionOn){
	for (i=0; i<dim; i++){
	  ierr = xf_Error(xf_AddToLinQueue(DData->Qn_gu+i*nq*sr2, xfe_LinQTerm_PhiPhi, nq, 
					   dim, sr2, -1, MDR->gbigb_X+i, dim, -1.0, LinQRR));
	  if (ierr != xf_OK) return ierr;
	}
      }
    }
  }

  Need_Grad = ((RL_UL != NULL) || (RL_UR != NULL) || (RR_UL != NULL) || (RR_UR != NULL));

  /* Also include linearization of stabilization term */
  if ((StabViscL != NULL) && (!StabData->SkipDiffStab) && (Need_Grad) ){
    // RL_UL{n,k;m,a} -= PhiL{n,q}*0.5*AwStabL{i,q,k}*StabPhiL{q}*wn{i,q}*StabViscL_UL{m,a}
    // RR_UL{n,k;m,a} += PhiR{n,q}*0.5*AwStabL{i,q,k}*StabPhiL{q}*wn{i,q}*StabViscL_UL{m,a}
    // RL_UR{n,k;m,a} -= PhiL{n,q}*0.5*AwStabL{i,q,k}*StabPhiL{q}*wn{i,q}*StabViscL_UR{m,a}
    // RR_UR{n,k;m,a} += PhiR{n,q}*0.5*AwStabL{i,q,k}*StabPhiL{q}*wn{i,q}*StabViscL_UR{m,a}
    xf_2ColcMult_Set(DData->AwStabL+0*nq*sr, wn+0, StabData->StabPhiL, nq, sr, dim, 1, 0.5, DData->Aw);
    for (i=1; i<dim; i++)
      xf_2ColcMult_Add(DData->AwStabL+i*nq*sr, wn+i, StabData->StabPhiL, nq, sr, dim, 1, 0.5, DData->Aw);

    xf_MTxM_Set(PhiDataL->Phi, DData->Aw, nL, nq, sr, T);
    if ((RL_UL != NULL) && (StabViscL_UL != NULL))
      xf_BlockOutProd_Sub(T, StabViscL_UL, nL, sr, nL, RL_UL);
    if ((RL_UR != NULL) && (StabViscL_UR != NULL))
      xf_BlockOutProd_Sub(T, StabViscL_UR, nL, sr, nR, RL_UR);

    xf_MTxM_Set(PhiDataR->Phi, DData->Aw, nR, nq, sr, T);
    if ((RR_UL != NULL) && (StabViscL_UL != NULL))
      xf_BlockOutProd_Add(T, StabViscL_UL, nR, sr, nL, RR_UL);
    if ((RR_UR != NULL) && (StabViscL_UR != NULL))
      xf_BlockOutProd_Add(T, StabViscL_UR, nR, sr, nR, RR_UR);
  }
  if ((StabViscR != NULL) && (!StabData->SkipDiffStab) && (Need_Grad) ){
    // RL_UL{n,k;m,a} -= PhiL{n,q}*0.5*AwStabR{i,q,k}*StabPhiR{q}*wn{i,q}*StabViscR_UL{m,a}
    // RR_UL{n,k;m,a} += PhiR{n,q}*0.5*AwStabR{i,q,k}*StabPhiR{q}*wn{i,q}*StabViscR_UL{m,a}
    // RL_UR{n,k;m,a} -= PhiL{n,q}*0.5*AwStabR{i,q,k}*StabPhiR{q}*wn{i,q}*StabViscR_UR{m,a}
    // RR_UR{n,k;m,a} += PhiR{n,q}*0.5*AwStabR{i,q,k}*StabPhiR{q}*wn{i,q}*StabViscR_UR{m,a}
    xf_2ColcMult_Set(DData->AwStabR+0*nq*sr, wn+0, StabData->StabPhiR, nq, sr, dim, 1, 0.5, DData->Aw);
    for (i=1; i<dim; i++)
      xf_2ColcMult_Add(DData->AwStabR+i*nq*sr, wn+i, StabData->StabPhiR, nq, sr, dim, 1, 0.5, DData->Aw);

    xf_MTxM_Set(PhiDataL->Phi, DData->Aw, nL, nq, sr, T);
    if ( (RL_UL != NULL) && (StabViscR_UL != NULL) )
      xf_BlockOutProd_Sub(T, StabViscR_UL, nL, sr, nL, RL_UL);
    if ( (RL_UR != NULL) && (StabViscR_UR != NULL) )
      xf_BlockOutProd_Sub(T, StabViscR_UR, nL, sr, nR, RL_UR);

    xf_MTxM_Set(PhiDataR->Phi, DData->Aw, nR, nq, sr, T);
    if ( (RR_UL != NULL) && (StabViscR_UL != NULL) )
      xf_BlockOutProd_Add(T, StabViscR_UL, nR, sr, nL, RR_UL);
    if ( (RR_UR != NULL) && (StabViscR_UR != NULL) )
      xf_BlockOutProd_Add(T, StabViscR_UR, nR, sr, nR, RR_UR);
  }


  // Need u-flux: uhat = 0.5*(uL+uR) + alpha*(uL-uR)
  ierr = xf_Error(xf_DiffFluxUJump(DiffDisc, &DData->alpha));
  if (ierr != xf_OK) return ierr;
  

  /* Need jumps, so calculate:
     dunL = (uL - uhat)*nL = (0.5-alpha) * (uL - uR) * nL 
     dunR = (uR - uhat)*nR = (0.5+alpha) * (uR - uL) * nR = (0.5+alpha) * (uL - uR) * nL
     Note: quad weights are included
     
     Derivatives defined according to:
     dunL_uL = duL_uL*nL
     dunL_uR = duL_uR*nL
     dunR_uL = duR_uL*nL (yes, this is nL)
     dunR_uR = duR_uR*nL (same here)
  */
  DData->duL_uL =  (0.5-DData->alpha); DData->duL_uR = -(0.5-DData->alpha);
  DData->duR_uR = -(0.5+DData->alpha); DData->duR_uL =  (0.5+DData->alpha);
  for (i=0; i<dim; i++)
    for (iq=0;iq<nq; iq++)
      for (k=0, nval=wn[dim*iq+i]; k<sr; k++){
	DData->dunL[(nq*i+iq)*sr+k] = DData->duL_uL*(uL[iq*sr+k]-uR[iq*sr+k])*nval;
	DData->dunR[(nq*i+iq)*sr+k] = DData->duR_uL*(uL[iq*sr+k]-uR[iq*sr+k])*nval;
      }

  // calculate ALdu [and A_uLdu]
  ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm, nDiff, IParam, RParam, ResMetricL, 
				  StabViscL, nq, uL, guL, DData->dunL, xfe_True, MDL,
				  DData->ALdunL, NULL, DData->A_uLdunL, DData->AwStabL, 
				  &DData->ConstAL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // same for R
  ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm, nDiff, IParam, RParam, ResMetricR, 
				  StabViscR, nq, uR, guR, DData->dunR, xfe_True, MDR, 
				  DData->ARdunR, NULL, DData->A_uRdunR, DData->AwStabR, 
				  &DData->ConstAR, NULL));
  if (ierr != xf_OK) return ierr;

  // Calculate DData->N = normalized wn; and FaceArea
  for (iq=0, FaceArea=0.; iq<nq; iq++){
    for (i=0, nval=0.; i<dim; i++) nval += wn[dim*iq+i]*wn[dim*iq+i];
    FaceArea += (nval = sqrt(nval));
    for (i=0; i<dim; i++) DData->N[dim*iq+i] = wn[dim*iq+i]/nval;
    if (CF != NULL) CF[iq] *= nval; // face area term included in Coupling here
  }


  /*-----------------------*/
  /* Viscous Stabilization */
  /*    (BR2, IP, etc.)    */
  /*-----------------------*/

  if (DiffDisc == xfe_DiffDiscIP){

    /* IP = Interior Penalty viscous discretization

       Stabilization terms:

       Qn{q,k} -= eta*N{i,q}*0.5*[ ihL*ALdunL{i,q,k} + ihR*ARdunR{i,q,k}]
       Qn_uL{q,k;a} -= eta*N{i,q}*0.5 * [ ihL*A_uLdunL{i,q,k,a} + ihL*AL{i,j,q,k,a}*wn{q,j}*duL_uL
                                                                + ihR*AR{i,j,q,k,a}*wn{q,j}*duR_uL]
       Qn_uR{q,k;a} -= eta*N{i,q}*0.5 * [                         ihL*AL{i,j,q,k,a}*wn{q,j}*duL_uR
                                        + ihR*A_uRdunR{i,q,k,a} + ihR*AR{i,j,q,k,a}*wn{q,j}*duR_uR ]
       N{i,q} = normalized wn{i,q}
    */

    // need eta for IP stabilization
    Basis = min(Mesh->ElemGroup[egrpL].QBasis, Mesh->ElemGroup[egrpR].QBasis);
    Order = max(ResOrderL, ResOrderR); // accounts for possible residual p-dependence
    nFace = max(Mesh->ElemGroup[egrpL].nFace[elemL], Mesh->ElemGroup[egrpR].nFace[elemR]);
    ierr = xf_Error(xf_JumpPenaltyEta(DiffDisc, Basis, Order, nFace, &eta));
    if (ierr != xf_OK) return ierr;

    // ih = average (1/h) normal to face
    ihL = FaceArea/ElemVolL;
    ihR = FaceArea/ElemVolR;


    //  Qn{q,k} (see above)
    for (i=0; i<dim; i++){
      xf_ColcMult_Add(DData->ALdunL+sr*nq*i, DData->N+i, nq, sr, dim,
		      -eta*0.5*ihL, DData->Qn);
      xf_ColcMult_Add(DData->ARdunR+sr*nq*i, DData->N+i, nq, sr, dim,
		      -eta*0.5*ihR, DData->Qn);
    }
  
    /* Add Q term to R:
       RL{n,k} -= PhiL{n,q}*Qn{q,k},  sum over q
       RR{n,k} += PhiR{n,q}*Qn{q,k},  sum over q (+ because nR = -nL)
    */
    if (RL != NULL)
      xf_MTxM_Sub(PhiDataL->Phi, DData->Qn, nL, nq, sr, RL);
    if (RR != NULL)
      xf_MTxM_Add(PhiDataR->Phi, DData->Qn, nR, nq, sr, RR);
  
  
    /*
      Qn_uL{q,k;a} (see above)
      RL_UL{n,k;m,a} -= PhiL{n,q}*[Qn_uL{q,k;a}*PhiL{q,m}]
      RR_UL{n,k;m,a} += PhiR{n,q}*[Qn_uL{q,k;a}*PhiL{q,m}]
    */
    if ((RL_UL != NULL) || (RR_UL != NULL)){
      for (k=0; k<nq*sr2; k++) DData->Qn_u[k] = 0.;
      for (i=0; i<dim; i++){
	for (j=0; j<dim; j++){
	  xf_2ColcMult_Add(DData->AL+(i*dim+j)*nq*sr2, wn+j, DData->N+i, nq, sr2,
			   dim, dim, -eta*0.5*ihL*DData->duL_uL, DData->Qn_u);
	  xf_2ColcMult_Add(DData->AR+(i*dim+j)*nq*sr2, wn+j, DData->N+i, nq, sr2,
			   dim, dim, -eta*0.5*ihR*DData->duR_uL, DData->Qn_u);
	} // j
	if (!DData->ConstAL)
	  xf_ColcMult_Add(DData->A_uLdunL+sr2*nq*i, DData->N+i, nq,
			  sr2, dim, -eta*0.5*ihL, DData->Qn_u);
      } // i
      if (RL_UL != NULL){ // RL_UL
	for (n=0; n<nL; n++){
	  xf_ColMult_Set(PhiDataL->Phi, PhiDataL->Phi+n, nq, nL, nL, T);
	  xf_MTxM_Sub(T, DData->Qn_u, nL, nq, sr2, RL_UL+n*nL*sr2);
	} // n
      }
      if (RR_UL != NULL){ // RR_UL
	for (n=0; n<nR; n++){
	  xf_ColMult_Set(PhiDataL->Phi, PhiDataR->Phi+n, nq, nL, nR, T);
	  xf_MTxM_Add(T, DData->Qn_u, nL, nq, sr2, RR_UL+n*nL*sr2);
	} // n
      }
    }

    /*
      Qn_uR{q,k;a} (see above)
      RL_UR{n,k;m,a} -= PhiL{n,q}*[Qn_uR{q,k;a}*PhiR{q,m}]
      RR_UR{n,k;m,a} += PhiR{n,q}*[Qn_uR{q,k;a}*PhiR{q,m}]
    */
    if ((RL_UR != NULL) || (RR_UR != NULL)){
      for (k=0; k<nq*sr2; k++) DData->Qn_u[k] = 0.;
      for (i=0; i<dim; i++){
	for (j=0; j<dim; j++){
	  xf_2ColcMult_Add(DData->AL+(i*dim+j)*nq*sr2, wn+j, DData->N+i, nq, sr2,
			   dim, dim, -eta*0.5*ihL*DData->duL_uR, DData->Qn_u);
	  xf_2ColcMult_Add(DData->AR+(i*dim+j)*nq*sr2, wn+j, DData->N+i, nq, sr2,
			   dim, dim, -eta*0.5*ihR*DData->duR_uR, DData->Qn_u);
	} // j
	if (!DData->ConstAR)
	  xf_ColcMult_Add(DData->A_uRdunR+sr2*nq*i, DData->N+i, nq,
			  sr2, dim, -eta*0.5*ihR, DData->Qn_u);
      } // i
      if (RL_UR != NULL){ // RL_UR
	for (n=0; n<nL; n++){
	  xf_ColMult_Set(PhiDataR->Phi, PhiDataL->Phi+n, nq, nR, nL, T);
	  xf_MTxM_Sub(T, DData->Qn_u, nR, nq, sr2, RL_UR+n*nR*sr2);
	} // n
      }
      if (RR_UR != NULL){ // RR_UR
	for (n=0; n<nR; n++){
	  xf_ColMult_Set(PhiDataR->Phi, PhiDataR->Phi+n, nq, nR, nR, T);
	  xf_MTxM_Add(T, DData->Qn_u, nR, nq, sr2, RR_UR+n*nR*sr2);
	} // n
      }
    }

    /* Include linearization of stabilization term */
    if ((StabViscL != NULL) && (!StabData->SkipDiffStab) && (Need_Grad)){
      // RL_UL{n,k;m,a} += PhiL{n,q}*eta*N{i,q}*0.5*ihL*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UL{m,a}
      // RR_UL{n,k;m,a} -= PhiR{n,q}*eta*N{i,q}*0.5*ihL*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UL{m,a}
      // RL_UR{n,k;m,a} += PhiL{n,q}*eta*N{i,q}*0.5*ihL*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UR{m,a}
      // RR_UR{n,k;m,a} -= PhiR{n,q}*eta*N{i,q}*0.5*ihL*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UR{m,a}
      // Note, quad weights are included in AwStabL (due to dunL)
      xf_2ColcMult_Set(DData->AwStabL+0*nq*sr, DData->N+0, StabData->StabPhiL,
		       nq, sr, dim, 1, 0.5*eta*ihL, DData->Aw);
      for (i=1; i<dim; i++)
	xf_2ColcMult_Add(DData->AwStabL+i*nq*sr, DData->N+i, StabData->StabPhiL,
			 nq, sr, dim, 1, 0.5*eta*ihL, DData->Aw);

      xf_MTxM_Set(PhiDataL->Phi, DData->Aw, nL, nq, sr, T);
      if ( (RL_UL != NULL) && (StabViscL_UL != NULL))
	xf_BlockOutProd_Add(T, StabViscL_UL, nL, sr, nL, RL_UL);
      if ( (RL_UR != NULL) && (StabViscL_UR != NULL))
	xf_BlockOutProd_Add(T, StabViscL_UR, nL, sr, nR, RL_UR);

      xf_MTxM_Set(PhiDataR->Phi, DData->Aw, nR, nq, sr, T);
      if ( (RR_UL != NULL) && (StabViscL_UL != NULL))
	xf_BlockOutProd_Sub(T, StabViscL_UL, nR, sr, nL, RR_UL);
      if ( (RR_UR != NULL) && (StabViscL_UR != NULL))
	xf_BlockOutProd_Sub(T, StabViscL_UR, nR, sr, nR, RR_UR);

    }
    if ((StabViscR != NULL) && (!StabData->SkipDiffStab) && (Need_Grad)){
      // RL_UL{n,k;m,a} += PhiL{n,q}*eta*N{i,q}*0.5*ihR*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UL{m,a}
      // RR_UL{n,k;m,a} -= PhiR{n,q}*eta*N{i,q}*0.5*ihR*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UL{m,a}
      // RL_UR{n,k;m,a} += PhiL{n,q}*eta*N{i,q}*0.5*ihR*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UR{m,a}
      // RR_UR{n,k;m,a} -= PhiR{n,q}*eta*N{i,q}*0.5*ihR*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UR{m,a}
      // Note, quad weights are included in AwStabR (due to dunR)
      xf_2ColcMult_Set(DData->AwStabR+0*nq*sr, DData->N+0, StabData->StabPhiR,
		       nq, sr, dim, 1, 0.5*eta*ihR, DData->Aw);
      for (i=1; i<dim; i++)
	xf_2ColcMult_Add(DData->AwStabR+i*nq*sr, DData->N+i, StabData->StabPhiR,
			 nq, sr, dim, 1, 0.5*eta*ihR, DData->Aw);

      xf_MTxM_Set(PhiDataL->Phi, DData->Aw, nL, nq, sr, T);
      if ( (RL_UL != NULL) && (StabViscR_UL != NULL) )
	xf_BlockOutProd_Add(T, StabViscR_UL, nL, sr, nL, RL_UL);
      if ( (RL_UR != NULL) && (StabViscR_UR != NULL) )
	xf_BlockOutProd_Add(T, StabViscR_UR, nL, sr, nR, RL_UR);

      xf_MTxM_Set(PhiDataR->Phi, DData->Aw, nR, nq, sr, T);
      if ( (RR_UL != NULL) && (StabViscR_UL != NULL) )
	xf_BlockOutProd_Sub(T, StabViscR_UL, nR, sr, nL, RR_UL);
      if ( (RR_UR != NULL) && (StabViscR_UR != NULL) )
	xf_BlockOutProd_Sub(T, StabViscR_UR, nR, sr, nR, RR_UR);
    }

    // Coupling is mu*eta*ih*|wn|, where mu*|wn| is already stored in CF
    if (CF != NULL) for (iq=0;iq<nq; iq++) CF[iq] *= eta*0.5*(ihL+ihR);

  }
  else if (DiffDisc == xfe_DiffDiscBR2){
 
    /* BR2 = Second form of Bassi & Rebay discretization

       Stabilization delta terms:

       Qn{q,k}   -= eta*dn{q,k} 
       dn{q,k}    = 0.5*(dnL{q,k} + dnR{q,k})
       dnL{q,k}   = ResPhiL{q,n}*DnL{i,n,k}*wn{i,q}
       dnR{q,k}   = ResPhiR{q,n}*DnR{i,n,k}*wn{i,q}
       DnL{i,n,k} = ResiML{n,m} * 0.5*ResPhiL{g,m}*ALdunL{i,g,k}
       DnR{i,n,k} = ResiMR{n,m} * 0.5*ResPhiR{g,m}*ARdunR{i,g,k}
       
       DnL_UL{i,n,k;o,a} = ResiML{n,m} * 0.5*ResPhiL{g,m}*(A_uLdunL{i,g,k,a} + AL{i,j,g,k,a}*wn{g,j}*duL_uL)*PhiL{g,o}
       DnL_UR{i,n,k;o,a} = ResiML{n,m} * 0.5*ResPhiL{g,m}*(                    AL{i,j,g,k,a}*wn{g,j}*duL_uR)*PhiR{g,o}
       DnR_UL{i,n,k;o,a} = ResiMR{n,m} * 0.5*ResPhiR{g,m}*(                    AR{i,j,g,k,a}*wn{g,j}*duR_uL)*PhiL{g,o}
       DnR_UR{i,n,k;o,a} = ResiMR{n,m} * 0.5*ResPhiR{g,m}*(A_uRdunR{i,g,k,a} + AR{i,j,g,k,a}*wn{g,j}*duR_uR)*PhiR{g,o}
       
       Since,
       RL{n,k} += PhiL{n,q}*eta*0.5*(ResPhiL{q,m}*DnL{i,m,k} + ResPhiR{q,m}*DnR{i,m,k})*wn{i,q}
                = SLL{i,n,m}*DnL{i,m,k} + SLR{i,n,m}*DnR{i,m,k}
       RR{n,k} -= PhiR{n,q}*eta*0.5*(ResPhiL{q,m}*DnL{i,m,k} + ResPhiR{q,m}*DnR{i,m,k})*wn{i,q}
                = SRL{i,n,m}*DnL{i,m,k} + SRR{i,n,m}*DnR{i,m,k}
       then
       RL_UL{n,k;o,a} += SLL{i,n,m}*DnL_UL{i,m,k;o,a} + SLR{i,n,m}*DnR_UL{i,m,k;o,a}
       RL_UR{n,k;o,a} += SLL{i,n,m}*DnL_UR{i,m,k;o,a} + SLR{i,n,m}*DnR_UR{i,m,k;o,a}
       RR_UL{n,k;o,a} -= SRL{i,n,m}*DnL_UL{i,m,k;o,a} + SRR{i,n,m}*DnR_UL{i,m,k;o,a}
       RR_UR{n,k;o,a} -= SRL{i,n,m}*DnL_UR{i,m,k;o,a} + SRR{i,n,m}*DnR_UR{i,m,k;o,a}

       N{i,q} = normalized wn{i,q}
       Note, in using the mass matrices, iMML must be multiplied by facL and iMMR by facR

       Note, ResPhiL, ResPhiR, ResiML, and ResiMR refer to quantities
       computed at a possibly different order than the state
       approximation.  This is useful for error estimation when the
       order of the residual needs to remain constant.

    */


    /* First add existing Qn to RL and RR, since we will be dealing
       directly with RL and RR.  Note that derivatives like Qn_uL,
       Qn_guL, etc. have already been added above. */
    if (RL != NULL)
      xf_MTxM_Sub(PhiDataL->Phi, DData->Qn, nL, nq, sr, RL);
    if (RR != NULL)
      xf_MTxM_Add(PhiDataR->Phi, DData->Qn, nR, nq, sr, RR);
   

    /* Now can deal with delta terms.  First pull off inverse mass
       matrices, iML, iMR.  These are at the residual order. */
    ierr = xf_Error(xf_ElemInvMassMatrix(All, egrpL, elemL, BasisL,
					 ResOrderL, NULL, NULL, &iMML, &facL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ElemInvMassMatrix(All, egrpR, elemR, BasisR,
					 ResOrderR, NULL, NULL, &iMMR, &facR));
    if (ierr != xf_OK) return ierr;
    
    /* Calculate delta coefficients DnL, DnR (note use of facL and facR for mass matrix scaling) */
    for (i=0; i<dim; i++){
      xf_MTxM_Set(ResPhiDataL->Phi, DData->ALdunL+i*nq*sr, RnL, nq, sr, T);
      xf_MxM_Set(iMML, T, RnL, RnL, sr, DData->DnL+i*RnL*sr);
    }
    for (k=0; k<dim*RnL*sr; k++) DData->DnL[k] *= 0.5*facL;
    for (i=0; i<dim; i++){
      xf_MTxM_Set(ResPhiDataR->Phi, DData->ARdunR+i*nq*sr, RnR, nq, sr, T);
      xf_MxM_Set(iMMR, T, RnR, RnR, sr, DData->DnR+i*RnR*sr);
    }
    for (k=0; k<dim*RnR*sr; k++) DData->DnR[k] *= 0.5*facR;
    
    // Determine eta
    nFace = max(Mesh->ElemGroup[egrpL].nFace[elemL], Mesh->ElemGroup[egrpR].nFace[elemR]);
    eta = ViscStabFactor * ((real) nFace);
    // ihL and ihR used to set coupling
    ihL = FaceArea/ElemVolL;
    ihR = FaceArea/ElemVolR;


    /* Calculate:
         SLL{i,n,m} = 0.5*eta*PhiL{n,q}*ResPhiL{q,m}*wn{i,q}
         SLR{i,n,m} = 0.5*eta*PhiL{n,q}*ResPhiR{q,m}*wn{i,q}
         SRL{i,n,m} = 0.5*eta*PhiR{n,q}*ResPhiL{q,m}*wn{i,q}
         SRR{i,n,m} = 0.5*eta*PhiR{n,q}*ResPhiR{q,m}*wn{i,q}
    */
    n2 = max(nL,nR)*max(RnL,RnR); // for ease of indexing into S**
    for (i=0; i<dim; i++){
      xf_ColcMult_Set(PhiDataL->Phi, wn+i, nq,  nL, dim, 0.5*eta, T); // row
      xf_MTxM_Set(T, ResPhiDataL->Phi, nL, nq, RnL, DData->SLL+i*n2); // col

      xf_ColcMult_Set(PhiDataL->Phi, wn+i, nq,  nL, dim, 0.5*eta, T); // row
      xf_MTxM_Set(T, ResPhiDataR->Phi, nL, nq, RnR, DData->SLR+i*n2); // col

      xf_ColcMult_Set(PhiDataR->Phi, wn+i, nq,  nR, dim, 0.5*eta, T); // row
      xf_MTxM_Set(T, ResPhiDataL->Phi, nR, nq, RnL, DData->SRL+i*n2); // col

      xf_ColcMult_Set(PhiDataR->Phi, wn+i, nq,  nR, dim, 0.5*eta, T); // row
      xf_MTxM_Set(T, ResPhiDataR->Phi, nR, nq, RnR, DData->SRR+i*n2); // col
    }

    /* Add to RL and RR:
         RL{n,k} += SLL{i,n,m}*DnL{i,m,k} + SLR{i,n,m}*DnR{i,m,k}
	 RR{n,k} -= SRL{i,n,m}*DnL{i,m,k} + SRR{i,n,m}*DnR{i,m,k}
    */
    if (RL != NULL){
      for (i=0; i<dim; i++){
	xf_MxM_Add(DData->SLL+i*n2, DData->DnL+i*RnL*sr, nL, RnL, sr, RL);
	xf_MxM_Add(DData->SLR+i*n2, DData->DnR+i*RnR*sr, nL, RnR, sr, RL);
      }
    }
    if (RR != NULL){
      for (i=0; i<dim; i++){
	xf_MxM_Sub(DData->SRL+i*n2, DData->DnL+i*RnL*sr, nR, RnL, sr, RR);
	xf_MxM_Sub(DData->SRR+i*n2, DData->DnR+i*RnR*sr, nR, RnR, sr, RR);
      }
    }

    
    /* SLL contribution to derivatives:
         PL{i,q,n} = SLL{i,n,m}*iMML{m,o}*ResPhiL{q,o}
         RL_UL{n,k;o,a} += PL{i,n,q}*0.5*(A_uLdunL{i,q,k,a} + AL{i,j,q,k,a}*wn{q,j}*duL_uL)*PhiL{q,o}
	 RL_UR{n,k;o,a} += PL{i,n,q}*0.5*(                    AL{i,j,q,k,a}*wn{q,j}*duL_uR)*PhiR{q,o}
       SRL contribution to derivatives:
         PR{i,q,n} = SRL{i,n,m}*iMML{m,o}*ResPhiL{q,o}
         RR_UL{n,k;o,a} -= PR{i,n,q}*0.5*(A_uLdunL{i,q,k,a} + AL{i,j,q,k,a}*wn{q,j}*duL_uL)*PhiL{q,o}
	 RR_UR{n,k;o,a} -= PR{i,n,q}*0.5*(                    AL{i,j,q,k,a}*wn{q,j}*duL_uR)*PhiR{q,o}
    */
    if (Need_Grad){
      for (i=0; i<dim; i++){
	// T{n,o} = SLL{i,n,m}*iMML{m,o}
	xf_MxM_Set(DData->SLL+i*n2, iMML, nL, RnL, RnL, T);
	for (k=0; k<nL*RnL; k++) T[k] *= facL; // facL is included here
	// PL{i,q,n} = ResPhiL{q,o}*T{o,n}
	xf_MxMT_Set(ResPhiDataL->Phi, T, nq, RnL, nL, DData->PL+i*nL*nq);
      }

      for (i=0; i<dim; i++){
	// T{n,o} = SRL{i,n,m}*iMML{m,o}
	xf_MxM_Set(DData->SRL+i*n2, iMML, nR, RnL, RnL, T);
	for (k=0; k<nR*RnL; k++) T[k] *= facL; // facL is included here
	// PR{q,n} = ResPhiL{q,o}*T{o,n}
	xf_MxMT_Set(ResPhiDataL->Phi, T, nq, RnL, nR, DData->PR+i*nR*nq);
      }	


      if ((RL_UL != NULL) || (RR_UL != NULL)){
	for (i=0; i<dim; i++){
	  // Qn_u{q,k,a} = 0.5*(A_uLdunL{i,q,k,a} + AL{i,j,q,k,a}*wn{q,j}*duL_uL)
	  for (k=0; k<nq*sr2; k++) DData->Qn_u[k] = 0.;
	  for (j=0; j<dim; j++){
	    xf_ColcMult_Add(DData->AL+(i*dim+j)*nq*sr2, wn+j, nq, sr2,
			     dim, 0.5*DData->duL_uL, DData->Qn_u);
	  } // j
	  if (!DData->ConstAL)
	    for (k=0; k<nq*sr2; k++) DData->Qn_u[k] += 0.5*DData->A_uLdunL[sr2*nq*i+k];
	  if (RL_UL != NULL){ // RL_UL
	    for (n=0; n<nL; n++){
	      xf_ColMult_Set(PhiDataL->Phi, DData->PL+i*nL*nq+n, nq, nL, nL, T); // T is nq x nL
	      xf_MTxM_Add(T, DData->Qn_u, nL, nq, sr2, RL_UL+n*nL*sr2);
	    } // n
	  }
	  if (RR_UL != NULL){ // RR_UL
	    for (n=0; n<nR; n++){
	      xf_ColMult_Set(PhiDataL->Phi, DData->PR+i*nR*nq+n, nq, nL, nR, T); // T is nq x nL
	      xf_MTxM_Sub(T, DData->Qn_u, nL, nq, sr2, RR_UL+n*nL*sr2);
	    } // n
	  }
	} // i
      }
      

      if ((RL_UR != NULL) || (RR_UR != NULL)){
	for (i=0; i<dim; i++){
	  // Qn_u{q,k,a} = 0.5*(                    AL{i,j,q,k,a}*wn{q,j}*duL_uR)
	  for (k=0; k<nq*sr2; k++) DData->Qn_u[k] = 0.;
	  for (j=0; j<dim; j++){
	    xf_ColcMult_Add(DData->AL+(i*dim+j)*nq*sr2, wn+j, nq, sr2,
			    dim, 0.5*DData->duL_uR, DData->Qn_u);
	  } // j
	  
	  if (RL_UR != NULL){ // RL_UR
	    for (n=0; n<nL; n++){
	      xf_ColMult_Set(PhiDataR->Phi, DData->PL+i*nL*nq+n, nq, nR, nL, T); // T is nq x nR
	      xf_MTxM_Add(T, DData->Qn_u, nR, nq, sr2, RL_UR+n*nR*sr2);
	    } // n
	  }
	  if (RR_UR != NULL){ // RR_UR
	    for (n=0; n<nR; n++){
	      xf_ColMult_Set(PhiDataR->Phi, DData->PR+i*nR*nq+n, nq, nR, nR, T); // T is nq x nR
	      xf_MTxM_Sub(T, DData->Qn_u, nR, nq, sr2, RR_UR+n*nR*sr2);
	    } // n
	  }
	} // i
      }

      /* Linearization of stabilization viscosity:
  	   RL_UL{n,k;o,a} += PL{i,n,q}*0.5*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UL{m,a}
	   RL_UR{n,k;o,a} += PL{i,n,q}*0.5*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UR{m,a}
  	   RR_UL{n,k;o,a} -= PR{i,n,q}*0.5*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UL{m,a}
	   RR_UR{n,k;o,a} -= PR{i,n,q}*0.5*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UR{m,a}
      */
      if ((StabViscL != NULL) && (StabViscL[0] != 0.0)){
	// assumption: StabVisc = 0 implies StabVisc_U = 0 (valid for continuous slope switch)
	for (i=0; i<dim; i++){
	  for (k=0; k<nq*sr; k++) DData->Aw[k] = DData->AwStabL[i*nq*sr+k]*0.5;
	  xf_ColMult(DData->Aw, StabData->StabPhiL, nq, sr, 1);
	  xf_MTxM_Set(DData->PL+i*nL*nq, DData->Aw, nL, nq, sr, T); // T is nL x sr
	  if ( (RL_UL != NULL) && (StabViscL_UL != NULL))
	    xf_BlockOutProd_Add(T, StabViscL_UL, nL, sr, nL, RL_UL);
	  if ( (RL_UR != NULL) && (StabViscL_UR != NULL))
	    xf_BlockOutProd_Add(T, StabViscL_UR, nL, sr, nR, RL_UR);
	  xf_MTxM_Set(DData->PR+i*nR*nq, DData->Aw, nR, nq, sr, T); // T is nR x sr
	  if ( (RR_UL != NULL) && (StabViscL_UL != NULL))
	    xf_BlockOutProd_Sub(T, StabViscL_UL, nR, sr, nL, RR_UL);
	  if ( (RR_UR != NULL) && (StabViscL_UR != NULL))
	    xf_BlockOutProd_Sub(T, StabViscL_UR, nR, sr, nR, RR_UR);
	}
      }
    }

     
    /* SLR contribution to derivatives:
         PL{i,q,n} = SLR{i,n,m}*iMMR{m,o}*ResPhiR{q,o}
         RL_UL{n,k;o,a} += PL{i,n,q}*0.5*(                    AR{i,j,q,k,a}*wn{q,j}*duR_uL)*PhiL{q,o}
	 RL_UR{n,k;o,a} += PL{i,n,q}*0.5*(A_uRdunR{i,q,k,a} + AR{i,j,q,k,a}*wn{q,j}*duR_uR)*PhiR{q,o}
       SRR contribution to derivatives:
         PR{i,q,n} = SRR{i,n,m}*iMMR{m,o}*ResPhiR{q,o}
         RR_UL{n,k;o,a} -= PR{i,n,q}*0.5*(                    AR{i,j,q,k,a}*wn{q,j}*duR_uL)*PhiL{q,o}
	 RR_UR{n,k;o,a} -= PR{i,n,q}*0.5*(A_uRdunR{i,q,k,a} + AR{i,j,q,k,a}*wn{q,j}*duR_uR)*PhiR{q,o}
    */
    if (Need_Grad){
      for (i=0; i<dim; i++){
	// T{n,o} = SLR{i,n,m}*iMMR{m,o}
	xf_MxM_Set(DData->SLR+i*n2, iMMR, nL, RnR, RnR, T);
	for (k=0; k<nL*RnR; k++) T[k] *= facR; // facR is included here
	// PL{i,q,n} = ResPhiR{q,o}*T{o,n}
	xf_MxMT_Set(ResPhiDataR->Phi, T, nq, RnR, nL, DData->PL+i*nL*nq);
      }

      for (i=0; i<dim; i++){
	// T{n,o} = SRR{i,n,m}*iMMR{m,o}
	xf_MxM_Set(DData->SRR+i*n2, iMMR, nR, RnR, RnR, T);
	for (k=0; k<nR*RnR; k++) T[k] *= facR; // facR is included here
	// PR{i,q,n} = ResPhiR{q,o}*T{o,n}
	xf_MxMT_Set(ResPhiDataR->Phi, T, nq, RnR, nR, DData->PR+i*nR*nq);
      }

      if ((RL_UL != NULL) || (RR_UL != NULL)){
	for (i=0; i<dim; i++){
	  // Qn_u{q,k,a} = 0.5*(                     AR{i,j,q,k,a}*wn{q,j}*duR_uL)
	  for (k=0; k<nq*sr2; k++) DData->Qn_u[k] = 0.;
	  for (j=0; j<dim; j++){
	    xf_ColcMult_Add(DData->AR+(i*dim+j)*nq*sr2, wn+j, nq, sr2,
			    dim, 0.5*DData->duR_uL, DData->Qn_u);
	  } // j
	  
	  if (RL_UL != NULL){ // RL_UL
	    for (n=0; n<nL; n++){
	      xf_ColMult_Set(PhiDataL->Phi, DData->PL+i*nL*nq+n, nq, nL, nL, T); // T is nq x nL
	      xf_MTxM_Add(T, DData->Qn_u, nL, nq, sr2, RL_UL+n*nL*sr2);
	    } // n
	  }
	  if (RR_UL != NULL){ // RR_UL
	    for (n=0; n<nR; n++){
	      xf_ColMult_Set(PhiDataL->Phi, DData->PR+i*nR*nq+n, nq, nL, nR, T); // T is nq x nL
	      xf_MTxM_Sub(T, DData->Qn_u, nL, nq, sr2, RR_UL+n*nL*sr2);
	    } // n
	  }
	} // i
      }

      if ((RL_UR != NULL) || (RR_UR != NULL)){
	for (i=0; i<dim; i++){
	  // Qn_u{q,k,a} = 0.5*(A_uRdunR{i,q,k,a} + AR{i,j,q,k,a}*wn{q,j}*duR_uR)
	  for (k=0; k<nq*sr2; k++) DData->Qn_u[k] = 0.;
	  for (j=0; j<dim; j++){
	    xf_ColcMult_Add(DData->AR+(i*dim+j)*nq*sr2, wn+j, nq, sr2,
			    dim, 0.5*DData->duR_uR, DData->Qn_u);
	  } // j
	  if (!DData->ConstAL)
	    for (k=0; k<nq*sr2; k++) DData->Qn_u[k] += 0.5*DData->A_uRdunR[sr2*nq*i+k];
	  if (RL_UR != NULL){ // RL_UR
	    for (n=0; n<nL; n++){
	      xf_ColMult_Set(PhiDataR->Phi, DData->PL+i*nL*nq+n, nq, nR, nL, T); // T is nq x nR
	      xf_MTxM_Add(T, DData->Qn_u, nR, nq, sr2, RL_UR+n*nR*sr2);
	    } // n
	  }
	  if (RR_UR != NULL){ // RR_UR
	    for (n=0; n<nR; n++){
	      xf_ColMult_Set(PhiDataR->Phi, DData->PR+i*nR*nq+n, nq, nR, nR, T); // T is nq x nR
	      xf_MTxM_Sub(T, DData->Qn_u, nR, nq, sr2, RR_UR+n*nR*sr2);
	    } // n
	  }
	} // i
      }

      /* Linearization of stabilization viscosity:
  	   RL_UL{n,k;o,a} += PL{i,n,q}*0.5*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UL{m,a}
	   RL_UR{n,k;o,a} += PL{i,n,q}*0.5*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UR{m,a}
  	   RR_UL{n,k;o,a} -= PR{i,n,q}*0.5*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UL{m,a}
	   RR_UR{n,k;o,a} -= PR{i,n,q}*0.5*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UR{m,a}
      */
      if ((StabViscL != NULL) && (StabViscL[0] != 0.0)){
	// assumption: StabVisc = 0 implies StabVisc_U = 0 (valid for continuous slope switch)
	for (i=0; i<dim; i++){
	  for (k=0; k<nq*sr; k++) DData->Aw[k] = DData->AwStabR[i*nq*sr+k]*0.5;
	  xf_ColMult(DData->Aw, StabData->StabPhiR, nq, sr, 1);
	  xf_MTxM_Set(DData->PL+i*nL*nq, DData->Aw, nL, nq, sr, T); // T is nL x sr
	  if ( (RL_UL != NULL) && (StabViscR_UL != NULL))
	    xf_BlockOutProd_Add(T, StabViscR_UL, nL, sr, nL, RL_UL);
	  if ( (RL_UR != NULL) && (StabViscR_UR != NULL))
	    xf_BlockOutProd_Add(T, StabViscR_UR, nL, sr, nR, RL_UR);
	  xf_MTxM_Set(DData->PR+i*nR*nq, DData->Aw, nR, nq, sr, T); // T is nR x sr
	  if ( (RR_UL != NULL) && (StabViscR_UL != NULL))
	    xf_BlockOutProd_Sub(T, StabViscR_UL, nR, sr, nL, RR_UL);
	  if ( (RR_UR != NULL) && (StabViscR_UR != NULL))
	    xf_BlockOutProd_Sub(T, StabViscR_UR, nR, sr, nR, RR_UR);
	}
      }
    }
     
    // Coupling is mu*eta*ih*|wn|, where mu*|wn| is already stored in CF
    if (CF != NULL) for (iq=0;iq<nq; iq++) CF[iq] *= eta*0.5*(ihL+ihR);    

  }
  else return xf_Error(xf_NOT_SUPPORTED); // only handle IP or BR2



  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateDiffBCData
int
xf_CreateDiffBCData(xf_DiffBCData **pDData)
{
  int ierr;
  xf_DiffBCData *DData;
  
  ierr = xf_Error(xf_Alloc( (void **) pDData, sizeof(xf_DiffBCData), 1));
  if (ierr != xf_OK) return ierr;

  DData = (*pDData);

  DData->Need_Grad = xfe_True;
  DData->ConstAB   = xfe_False;
  DData->alpha     = 0.0;

  DData->dun       = NULL; 
  DData->du_uI     = NULL;
  DData->du_gb     = NULL;
  DData->ABgu      = NULL; 
  DData->ABdun     = NULL; 
  DData->AB        = NULL; 
  DData->A_uBgu    = NULL; 
  DData->A_uBdun   = NULL; 
  DData->AB_uIdun  = NULL;
  DData->N         = NULL; 
  DData->Qn        = NULL; 
  DData->Qn_uI     = NULL; 
  DData->Qn_guI    = NULL; 
  DData->Qn_gb     = NULL; 
  DData->Qn_ggb    = NULL;
  DData->AwStab    = NULL; 
  DData->guB       = NULL; 
  DData->guB_guI   = NULL; 
  DData->Aw        = NULL; 
  DData->T         = NULL; 
  DData->Dn        = NULL; 
  DData->S         = NULL; 
  DData->P         = NULL; 
  DData->Vw        = NULL; 

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReAllocDiffBCData
int
xf_ReAllocDiffBCData(int sr, int nq, int dim, int nn, 
		     enum xfe_Bool Need_Grad, xf_DiffBCData *DData)
{
  int ierr;
  int sr2, Tsize;

  sr2 = sr*sr;
  DData->Need_Grad = Need_Grad;

  ierr = xf_Error(xf_ReAlloc( (void **) &DData->dun, dim*nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ABgu, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ABdun, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn, sr*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->guB, dim*nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->N , nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Dn, dim*nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // BR2 needs these even without R_U request
  Tsize = dim*nq*sr2;
  Tsize = max(Tsize, nn*nq);
  Tsize = max(Tsize, nn*sr2);
  Tsize = max(Tsize, nn*nn);
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->T, Tsize, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->S, dim*nn*nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->P, dim*nn*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // for output calculation
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Vw, dim*nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

	
  if (Need_Grad){ 
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->du_uI, nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->du_gb, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AB, sr2*nq*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uBgu , sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uBdun, sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AB_uIdun, sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_uI, nq*nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_guI, dim*nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_gb, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_ggb, dim*nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AwStab, dim*nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->guB_guI, sr2*nq*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Aw, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDiffBCData
void
xf_DestroyDiffBCData(xf_DiffBCData *DData){

  xf_Release( (void *) DData->dun);
  xf_Release( (void *) DData->du_uI);     
  xf_Release( (void *) DData->du_gb);      
  xf_Release( (void *) DData->ABgu);   
  xf_Release( (void *) DData->ABdun);    
  xf_Release( (void *) DData->AB);  
  xf_Release( (void *) DData->A_uBgu); 
  xf_Release( (void *) DData->A_uBdun);
  xf_Release( (void *) DData->AB_uIdun);
  xf_Release( (void *) DData->N);
  xf_Release( (void *) DData->Qn);
  xf_Release( (void *) DData->Qn_uI);   
  xf_Release( (void *) DData->Qn_guI);
  xf_Release( (void *) DData->Qn_gb);   
  xf_Release( (void *) DData->Qn_ggb);
  xf_Release( (void *) DData->AwStab);  
  xf_Release( (void *) DData->guB);
  xf_Release( (void *) DData->guB_guI);
  xf_Release( (void *) DData->Aw);
  xf_Release( (void *) DData->T);
  xf_Release( (void *) DData->Dn);
  xf_Release( (void *) DData->S);
  xf_Release( (void *) DData->P);
  xf_Release( (void *) DData->Vw);

  xf_Release( (void *) DData);
}

/******************************************************************/
//   FUNCTION Definition: xf_DiffFluxUBC
static int
xf_DiffFluxUBC(enum xfe_DiffDiscType DiffDisc, real *alpha){
  
  switch(DiffDisc){
  case xfe_DiffDiscIP:
  case xfe_DiffDiscBR2:
    (*alpha) = 0.0;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SetPrescribedQBC
static int
xf_SetPrescribedQBC(const int *SetIndex, int nq, int sr, int dim,
		    const real *guB, const real *guB_guI, real *Qn, 
                    real *Aw, real *Qn_g, real *Qn_uI, real *Qn_guI,
		    enum xfe_Bool ZeroFlag)
{
  // Account for BC-imposed values of viscous flux
  // Zero flag means we just want to zero out components that are Set
  int k, iq, i, l, ii, sr2;

  sr2 = sr*sr;

  for (k=0; k<sr; k++)
    if (SetIndex[k]){
      for (iq=0;iq<nq; iq++){
	if (Qn != NULL){
	  if (ZeroFlag) Qn[iq*sr+k] = 0.0;
	  else Qn[iq*sr+k] = guB[iq*sr+k];
	}
	if (Aw != NULL)
	  Aw[iq*sr+k] = 0.0;
        if (Qn_g != NULL)
          for (i=0; i<dim; i++) Qn_g[i*nq*sr+iq*sr+k] = 0.;
	if (Qn_guI != NULL)
	  for (i=0; i<dim; i++)
	    for (l=0, ii = (i*nq+iq)*sr2+k*sr; l<sr; l++){
	      if (ZeroFlag) Qn_guI[ii+l] = 0.0;
	      else Qn_guI[ii+l] = guB_guI[ii+l];
	    }
	if (Qn_uI != NULL)
	  for (l=0; l<sr; l++) Qn_uI[iq*sr2+k*sr+l] = 0.0;
      }
    }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DiffFluxQBC
int
xf_DiffFluxQBC(xf_All *All, int ibfgrp, int ibface, const xf_ResTerm *ResTerm, 
	       int nDiff, enum xfe_DiffDiscType DiffDisc, const int *IParam, 
	       const real *RParam, const xf_BasisData *PhiData, 
               const xf_BasisData *ResPhiData, const xf_StabData *StabData, 
	       real ElemVol, enum xfe_BasisType Basis, int Order, int nq, const real *wn,  
	       const real *xglob, real *uI, real *uB, const real *uB_uI, real *guI, 
	       const xf_MotionData *MD, xf_OutputEvalData *OutputEval, real *ER, real *ER_U, 
	       xf_LinQueueData *LinQ, xf_DiffBCData *DData, real *CF)
{
  int ierr, i, j, k, l, dim, ii, jj, sr, sr2;
  int egrp, elem, iq, *SetIndex, nn, n, n2;
  int ResOrder, Rnn;
  int *FluxMoments = NULL;
  enum xfe_Bool AFlag, Identity_gu_guI = xfe_False;
  enum xfe_Bool MotionOn;
  real eta, fac, nval, FaceArea, ih, val, t;
  const real *gu, *gu_guI;
  real ViscStabFactor;
  real *T, xfac, rtemp;
  real *FluxWeights = NULL, *EV_U = NULL, *EV_G = NULL, *Value = NULL;
  real *ResMetric;
  real *StabVisc = NULL, *StabVisc_U = NULL;
  real *iMM;
  real *w = NULL;
  xf_BFace BFace;
  xf_BC   *BC;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;

  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;
  BC     = EqnSet->BCs[0].BC+ibfgrp;

  dim  = Mesh->Dim;
  sr   = EqnSet->StateRank;
  sr2  = sr*sr;

  // determine stabilization factor
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "ViscStabFactor", &ViscStabFactor);
  if (ierr == xf_NOT_FOUND)  ViscStabFactor = 1.0; // backwards-compatibility
  else if (ierr != xf_OK) return xf_Error(ierr); 
  
  // output information
  if (OutputEval != NULL){
    Value       = OutputEval->Value;
    EV_U        = OutputEval->EV_U;
    EV_G        = OutputEval->EV_G;
    FluxWeights = OutputEval->FluxWeights;
    FluxMoments = OutputEval->FluxMoments;
  }

  // element group info
  BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
  egrp  = BFace.ElemGroup;
  elem  = BFace.Elem;

  // pull off temporary matrices
  T = DData->T;

  nn = PhiData->nn;

  // residual order (possibly different from state order)
  ResOrder = ResPhiData->Order;
  Rnn = ResPhiData->nn;

  // MD passed-in means motion is on
  MotionOn = (MD != NULL);

  // allocate temp memory
  ierr = xf_Error(xf_Alloc((void **) &SetIndex, sr, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<sr; k++) SetIndex[k] = 0;


  // pull off data required for stabilization
  if (StabData != NULL){
    StabVisc    = StabData->StabVisc;
    StabVisc_U  = StabData->StabVisc_U;
    ResMetric   = StabData->ResMetric;
  }
  else StabVisc = NULL;


  /* For a general diffusion discretization:

     Qn{q,k} = ABgu{i,q,k}*wn{q,i} - (visc. stab)

     where (visc. stab) depends on the paticular form used

     AB = A matrix calculated using uB
     gu = grad u as prescribed by BCs

     In the presence of mesh motion, instead of gu, use:
       gu{j,q,a} - uB{q,a}*gbigb_X{q,j}
     but only when the viscous BC flux is not set directly

     Quad weights are included through the use of wn
  */

  // obtain viscous BC flux or grad u, using eqnset
  // TODO: Normal passed in should be transformed, before and after, when Motion is on
  //       (OK for now because currently no visc flux BCs explicitly use the normal)
  ierr = xf_Error(xf_EqnSetDiffFBC(EqnSet, BC, IParam, RParam, nq, wn, 
				   xglob, guI, DData->guB, DData->guB_guI, 
				   &AFlag, SetIndex));
  if (ierr != xf_OK) return ierr;

  /* AFlag == True means BC imposed condition on visc flux itself, not
     on grad u.  In this case, the entries of guB for which SetIndex
     is 1 are the set viscous flux components (already dotted with the
     normal); the remaining entries of the viscous flux vector are
     taken from AB*guI.  AFlag == False means that the BC imposed a
     condition on grad u; SetIndex is not used in this case. */
  if (AFlag){
    gu = guI;
    Identity_gu_guI = xfe_True;
  }
  else{
    Identity_gu_guI = xfe_False;
    gu = DData->guB;
    gu_guI = DData->guB_guI;
  }
  
  // zero out coupling
  if (CF != NULL) for (iq=0;iq<nq; iq++) CF[iq] = 0.;

  // w is the vector that AB multiplies
  w = DData->dun;  // using this space as temporary storage
  for (k=0;k<dim*nq*sr;k++) w[k] = gu[k]; // set w = gradient
  if (MotionOn)      // modify state gradient in presence of mesh motion
    for (j=0; j<dim; j++)
      xf_ColMult_Sub(uB, MD->gbigb_X+j, nq, sr, dim, w+nq*sr*j);

  // Calculate AB*gu, A_uB*gu ; uB is now physical
  ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm, nDiff, IParam, RParam, ResMetric, 
				  StabVisc, nq, uB, guI, w, xfe_False, MD, DData->ABgu, 
				  DData->AB, DData->A_uBgu, DData->AwStab, 
				  &DData->ConstAB, CF));
  if (ierr != xf_OK) return ierr;


  // flux term first: Qn{q,k} = ABgu{i,q,k}*wn{q,i}
  xf_ColMult_Set(DData->ABgu+0*nq*sr, wn+0, nq, sr, dim, DData->Qn);
  for (i=1; i<dim; i++)
    xf_ColMult_Add(DData->ABgu+i*nq*sr, wn+i, nq, sr, dim, DData->Qn);

  /* Derivatives of flux term */
  if ((ER_U != NULL) || (EV_U != NULL)){
    /* 
       Qn_guI{j,q,k;a} = AB{i,j,q,k,l}*wn{q,i} * gu_guI{j,d,q,l,a}
       Qn_uI{q,k;a} = (A_uBgu{i,q,k,l}-AB{i,j,q,k,l}*gbigb_X{q,j})*wn{q,i} * uB_uI{q,l,a}
                                      \____ only if MotionOn ___/
       Below is for GCL:
       Qn_ggb{j,q,k;} = AB{i,j,q,k,l}*wn{q,i} * (-uB{q,l}/gb{q})
       Qn_gb{q,k;} = -Qn{q,k}/gb{q} + Qn_uI{q,k;a}*(-uI{q,a}/gb{q})
    */

    // Qn_guI -- first use gu_guI = identity
    for (j=0; j<dim; j++){ 
      jj = j*nq*sr2;
      xf_ColMult_Set(DData->AB+j*nq*sr2, wn+0, nq, sr2, dim, DData->Qn_guI+jj);
      for (i=1; i<dim; i++)
	xf_ColMult_Add(DData->AB+(i*dim+j)*nq*sr2, wn+i, nq, sr2, dim, DData->Qn_guI+jj);
    }
    
    // GCL specific, Qn_ggb{j,q,k;} = Qn_guI{j,q,k;l}*(-uB{q,l}/gb{q})
    if (EV_G != NULL){
      for (j=0; j<dim; j++)
        for (iq=0; iq<nq; iq++){
          jj = j*nq*sr+iq*sr;
          xf_MxV_Set(DData->Qn_guI+j*nq*sr2+iq*sr2, uB+iq*sr, sr, sr, DData->Qn_ggb+jj);
          for (k=0; k<sr; k++) DData->Qn_ggb[jj+k] *= -1./MD->gb[iq];
        } 
    }

    if ((!DData->ConstAB) || (MotionOn)){
      // Qn_uI
      for (k=0; k<nq*sr2; k++) T[k] = 0.;
      if (!DData->ConstAB){
	for (i=0; i<dim; i++)
	  xf_ColMult_Add(DData->A_uBgu+i*nq*sr2, wn+i, nq, sr2, dim, T);
      }
      if (MotionOn){
	for (i=0; i<dim; i++)
	  for (j=0; j<dim; j++)
	    xf_2ColcMult_Add(DData->AB+(i*dim+j)*nq*sr2, wn+i, MD->gbigb_X+j, nq, 
			     sr2, dim, dim, -1., T);
      }
      xf_nMxM_Set(nq, T, uB_uI, sr, sr, sr, DData->Qn_uI);

      // GCL specific, Qn_gb{q,k;} = -Qn{q,k}/gb{q} + Qn_uI{q,k;a}*(-uI{q,a}/gb{q})
      if (EV_G != NULL){
        for (iq=0; iq<nq; iq++){
          xf_MxV_Set(DData->Qn_uI+iq*sr2, uI+iq*sr, sr, sr, DData->Qn_gb+iq*sr);
          for (k=0; k<sr; k++)
            DData->Qn_gb[iq*sr+k] = -(DData->Qn[iq*sr+k] + DData->Qn_gb[iq*sr+k])/MD->gb[iq];
        }
      }
    }
    
    if (!Identity_gu_guI){
      // carry out Qn_guI{j,q,k;a} *= gu_guI{j,d,q,l,a}
      for (j=0; j<dim; j++){
	jj = j*nq*sr2;
	xf_nMxM_Set(nq, DData->Qn_guI, gu_guI+jj, sr, sr, sr, T+jj);
	for (i=1; i<dim; i++){
	  ii = i*nq*sr2;
	  xf_nMxM_Add(nq, DData->Qn_guI+ii, gu_guI+i*dim*nq*sr2+jj, sr, sr, sr, T+jj);
	}
      } // j
      for (i=0; i<dim*nq*sr2; i++) DData->Qn_guI[i] = T[i];
    }

    
    /* Multiply AwStab by quad weights/normals and add to DData->Aw.
       Do this in preparation for linearizing the stabilization terms,
       before setting prescribed QBCs.  */
    if (StabData != NULL){
      xf_2ColMult_Set(DData->AwStab+0*nq*sr, wn+0, StabData->StabPhi, nq, sr, dim, 1, DData->Aw);
      for (i=1; i<dim; i++)
	xf_2ColMult_Add(DData->AwStab+i*nq*sr, wn+i, StabData->StabPhi, nq, sr, dim, 1, DData->Aw);
    }

    // before adding to ER_U, zero out derivatives for prescribed flux entries
    if (AFlag){
      ierr = xf_Error(xf_SetPrescribedQBC(SetIndex, nq, sr, dim, DData->guB,
					  DData->guB_guI, NULL, DData->Aw, NULL,
					  DData->Qn_uI, DData->Qn_guI, xfe_False));
      if (ierr != xf_OK) return ierr;
      if (EV_G != NULL){ // zero out GCL-specific entries too
        ierr = xf_Error(xf_SetPrescribedQBC(SetIndex, nq, sr, dim, NULL,
                                            NULL, NULL, DData->Qn_gb,
                                            DData->Qn_ggb, NULL, NULL, xfe_True));
        if (ierr != xf_OK) return ierr;
      }
    }
     
    /*
      Add to ER_U
      ER_U{n,k;m,a} -= Phi{n,q}*[Qn_uI{q,k;a}*Phi{q,m}
                                 + Qn_guI{j,q,k;a}*gPhi{j,q,m}]
    */
    if (ER_U != NULL){
      if ((!DData->ConstAB) || (MotionOn)){
	ierr = xf_Error(xf_AddToLinQueue(DData->Qn_uI, xfe_LinQTerm_PhiPhi, 
					 nq, dim, sr2, -1, NULL, 0, -1.0, LinQ));
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_AddToLinQueue(DData->Qn_guI, xfe_LinQTerm_PhiGPhi, 
				       nq, dim, sr2, -1, NULL, 0, -1.0, LinQ));
      if (ierr != xf_OK) return ierr;

      /* Also include linearization of stabilization term */
      if ((StabVisc != NULL) && (!StabData->SkipDiffStab)){
	// ER_U{n,k;m,a} -= Phi{n,q}*AwStab{i,q,k}*wn{i,q}*StabVisc_U{m,a}
	xf_MTxM_Set(PhiData->Phi, DData->Aw, nn, nq, sr, T);
	xf_BlockOutProd_Sub(T, StabVisc_U, nn, sr, nn, ER_U);
      }
    }

    
    /* 
       Add to EV_U:
       EV_U{n,l} -= FluxWeights{k}* [Qn_uI{q,k;l}*Phi{q,n}
                  + Qn_guI{i,q,k;l}*gPhi{i,q,n}] [*xglob{q,FluxMoments{k}}]
       EV_G{n} -= FluxWeights{k}* [Qn_gb{q,k;}*Phi{q,n}
                  + Qn_ggb{i,q,k;}*gPhi{i,q,n}] [*xglob{q,FluxMoments{k}}]
                  ** adj consistency fix: Qn_ggb is now zero **
    */
    if (EV_U != NULL){
      for (n=0; n<nn; n++)
	for (k=0; k<sr; k++)
          for (iq=0; iq<nq; iq++){
            xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
            rtemp =  FluxWeights[k]*PhiData->Phi[iq*nn+n]*xfac;
            for (l=0; l<sr; l++) EV_U[n*sr+l] -= rtemp * DData->Qn_uI[iq*sr2+k*sr+l];
            if (EV_G != NULL) EV_G[n] -= rtemp * DData->Qn_gb[iq*sr+k];
            for (i=0; i<dim; i++){
              rtemp = FluxWeights[k]*PhiData->gPhi[nn*nq*i+iq*nn+n]*xfac;
              for (l=0; l<sr; l++) EV_U[n*sr+l] -= rtemp * DData->Qn_guI[nq*sr2*i+iq*sr2+k*sr+l];
              // no dependence on grad of gb anymore
              // if (EV_G != NULL) EV_G[n] -= rtemp * DData->Qn_ggb[nq*sr*i+iq*sr+k];
            } // i
          } // iq
      
      /* Also include linearization of stabilization term (GCL not needed here, U is ref state)*/
      if ((StabVisc != NULL) && (!StabData->SkipDiffStab)){
	// EV_U{m,l} -= FluxWeights{k}*AwStab{i,q,k}*wn{i,q}*StabVisc_U{m,l}
	//              [*xglob{q,FluxMoments{k}}]
	// Note, Aw{q,k} = AwStab{i,q,k}*wn{i,q}
	for (iq=0, val=0.; iq<nq; iq++)
	  for (k=0, ii=iq*sr; k<sr; k++){
	    xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
	    val += FluxWeights[k]*DData->Aw[ii+k]*xfac;
	  }
	for (i=0; i<nn*sr; i++) EV_U[i] -= val*StabVisc_U[i];
      }
    }
  } // end of adding flux derivatives; diffusion stabilization is next
			     		     
  // NOTE: if boundary is deforming, might want to multiply uB by g here
			     
  // need bc flux in u: uhat = uB + alpha*(uI-uB);
  ierr = xf_Error(xf_DiffFluxUBC(DiffDisc, &DData->alpha));
  if (ierr != xf_OK) return ierr;

  // need jump in u: dun = (uI - uhat)*wn = (1-alpha)*(uI-uB)*wn
  fac = 1.-DData->alpha;
  for (i=0; i<dim; i++)
    for (iq=0;iq<nq; iq++)
      for (k=0, nval=wn[dim*iq+i]; k<sr; k++)
	DData->dun[(nq*i+iq)*sr+k] = fac*(uI[iq*sr+k]-uB[iq*sr+k])*nval;

  // du_uI{q,k,a} = (1-alpha) * (eye{q,k,a} - uB_uI{q,k,a})
  if (DData->du_uI != NULL)
    for (iq=0;iq<nq; iq++){
      for (k=0; k<sr2; k++)  DData->du_uI[iq*sr2+k] = -fac*uB_uI[iq*sr2+k];
      for (k=0; k<sr2; k+=(sr+1))  DData->du_uI[iq*sr2+k] += fac;
    }

  // du_gb{q,k;} = du_uB{q,k;a}*uB_gb{q,a;}
  // where  du_uB{q,k;a} = -(1-alpha)*eye{k,a} for each q
  // and    uB_gb{q,a;} = (uB{q,a} - uB_uI{q,a,l}*uI{q,l} )/gb{q}
  if (EV_G != NULL){
    xf_cV_Add(uB, -fac, nq*sr, xfe_Set, DData->du_gb); // du_gb{q,a} = -(1-alpha)*uB{q,a}
    xf_nMxM_Set(nq, uB_uI, uI, sr, sr, 1, T);   // T{q,a} = uB_uI{q,a,l}*uI{q,l}
    xf_cV_Add(T, fac, nq*sr, xfe_Add, DData->du_gb);   // du_gb{q,a} += (1-alpha)*T{q,a}
    xf_ColDiv(DData->du_gb, MD->gb, nq, sr, 1);        // du_gb{q,a} /= gb{q}
  }


  // Calculate AB, AB*dun, A_uB*dun ; uB is now physical
  ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm, nDiff, IParam, RParam, ResMetric,
				  StabVisc, nq, uB, guI, DData->dun, xfe_False, MD, 
				  DData->ABdun, NULL, DData->A_uBdun, DData->AwStab, 
				  &DData->ConstAB, NULL));
  if (ierr != xf_OK) return ierr;

  // Calculate AB_uI*dun = A_uB*dun * uB_uI
  if ((!DData->ConstAB) && (DData->Need_Grad)){
    for (i=0; i<dim; i++)
      xf_nMxM_Set(nq, DData->A_uBdun+sr2*nq*i, uB_uI, sr, sr, sr, 
		  DData->AB_uIdun+sr2*nq*i);
  }

  // Calculate DData->N = normalized wn; and FaceArea
  for (iq=0, FaceArea=0.; iq<nq; iq++){
    for (i=0, nval=0.; i<dim; i++) nval += wn[dim*iq+i]*wn[dim*iq+i];
    FaceArea += (nval = sqrt(nval));
    for (i=0; i<dim; i++) DData->N[dim*iq+i] = wn[dim*iq+i]/nval;
    if (CF != NULL) CF[iq] *= nval; // face area term included in Coupling here
  }


  /*-----------------------*/
  /* Viscous Stabilization */
  /*    (BR2, IP, etc.)    */
  /*-----------------------*/
  
  if (DiffDisc == xfe_DiffDiscIP){

    /* IP = Interior Penalty viscous discretization

       Stabilization terms:

       Qn{q,k} -= eta*ih*N{i,q}*ABdun{i,q,k}
       Qn_uI{q,k;a} -= eta*ih*N{i,q} * [AB_uI{i,q,k,a}
                     + AB{i,j,q,k,l}*du_uI{q,l,a}*wn{q,j}]
       N{i,q} = normalized wn{i,q}
    */

    if (EV_G != NULL) return xf_Error(xf_NOT_SUPPORTED); // for now

    // For stabilization, need eta (use ResOrder to account for possible p-dependence)
    ierr = xf_Error(xf_JumpPenaltyEta(DiffDisc, Mesh->ElemGroup[egrp].QBasis, 
				      ResOrder, Mesh->ElemGroup[egrp].nFace[elem], &eta));
    if (ierr != xf_OK) return ierr;
        
    // ih = average (1/h) normal to face
    ih = FaceArea/ElemVol;
    
    // Qn{q,k} -= eta*ih*N{i,q}*ABdun{i,q,k}
    for (i=0; i<dim; i++)
      xf_ColcMult_Add(DData->ABdun+sr*nq*i, DData->N+i, nq, sr, dim,
		      -eta*ih, DData->Qn);
    
    
    // Store flux components set by BC in Qn
    if (AFlag){
      ierr = xf_Error(xf_SetPrescribedQBC(SetIndex, nq, sr, dim, DData->guB,
					  DData->guB_guI, DData->Qn, NULL, NULL,
					  NULL, NULL, xfe_False));
      if (ierr != xf_OK) return ierr;
    }

    /* Add Q term to R:
       ER{n,k} -= Phi{n,q}*Qn{q,k},  sum over q
    */
    if (ER != NULL)
      xf_MTxM_Sub(PhiData->Phi, DData->Qn, nn, nq, sr, ER);
    
    /* Add Q term to Value:
       Value -= FluxWeights{k}*Qn{q,k}[*xglob{q,FluxMoments{k}}], sum over q,k
    */
    if (Value != NULL)
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++){
	  xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
	  (*Value) -= FluxWeights[k]*DData->Qn[iq*sr+k]*xfac;
	}
    
    if ((ER_U != NULL) || (EV_U != NULL)){
      /*
	Qn_uI{q,k;a} -= eta*ih*N{i,q} * [AB_uIdun{i,q,k,a}
 	              + AB{i,j,q,k,l}*du_uI{q,l,a}*wn{q,j}]
      */
      for (k=0; k<nq*sr2; k++) DData->Qn_uI[k] = 0.;
      for (i=0; i<dim; i++){
	for (j=0; j<dim; j++){
	  xf_nMxM_Set(nq, DData->AB+(i*dim+j)*nq*sr2, DData->du_uI, sr, sr, sr, T);
	  xf_2ColcMult_Add(T, wn+j, DData->N+i, nq, sr2, dim, dim,
			   -eta*ih, DData->Qn_uI); // T is nq x sr2
	} // j
	if (!DData->ConstAB)
	  xf_ColcMult_Add(DData->AB_uIdun+sr2*nq*i, DData->N+i, nq,
			  sr2, dim, -eta*ih, DData->Qn_uI);
      } // i
      
      // Set Aw{q,k} = N{i,q}*AwStab{i,q,k}*wn{i,q}
      // Note, quad weights (wn) are included in AwStab (due to dun)
      if (StabData != NULL){
	xf_2ColcMult_Set(DData->AwStab+0*nq*sr, DData->N+0, StabData->StabPhi, 
			 nq, sr, dim, 1, eta*ih, DData->Aw);
	for (i=1; i<dim; i++)
	  xf_2ColcMult_Add(DData->AwStab+i*nq*sr, DData->N+i, StabData->StabPhi, 
			   nq, sr, dim, 1, eta*ih, DData->Aw);
      }
      
      // before adding to ER_U, zero out derivatives for prescribed flux entries
      if (AFlag){
	ierr = xf_Error(xf_SetPrescribedQBC(SetIndex, nq, sr, dim, DData->guB,
					    DData->guB_guI, NULL, DData->Aw, NULL,
					    DData->Qn_uI, NULL, xfe_False));
	if (ierr != xf_OK) return ierr;
      }


      /*
	Add to ER_U:
	ER_U{n,k;m,a} -= Phi{n,q}*[Qn_uI{q,k;a}*Phi{q,m}]
      */
      if (ER_U != NULL){
	for (n=0; n<nn; n++){
	  xf_ColMult_Set(PhiData->Phi, PhiData->Phi+n, nq, nn, nn, T);
	  xf_MTxM_Sub(T, DData->Qn_uI, nn, nq, sr2, ER_U+n*nn*sr2);
	} // n
	
	/* Also include linearization of stabilization term */
	if ((StabVisc != NULL) && (!StabData->SkipDiffStab)){	  
	  // ER_U{n,k;m,a} += Phi{n,q}*eta*ih*N{i,q}*AwStab{i,q,k}*wn{i,q}*StabVisc_U{m,a}
	  // Note, Aw{q,k} = N{i,q}*AwStab{i,q,k}*wn{i,q}
	  xf_MTxM_Set(PhiData->Phi, DData->Aw, nn, nq, sr, T);
	  xf_BlockOutProd_Add(T, StabVisc_U, nn, sr, nn, ER_U);
	}
      }

      
      /*
	Add to EV_U:
	EV_U{m,l} -= FluxWeights{k}* [Qn_uI{q,k;l}*Phi{q,m}][*xglob{q,FluxMoments{k}}]
      */
      if (EV_U != NULL){
	//xf_MTxM_Set(PhiData->Phi, DData->Qn_uI, nn, nq, sr2, T);
	for (n=0; n<nn; n++)
	  for (k=0; k<sr; k++)
	    for (l=0; l<sr; l++)
	      for (iq=0; iq<nq; iq++){
		xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
		EV_U[n*sr+l] -= FluxWeights[k]*DData->Qn_uI[iq*sr2+k*sr+l]*PhiData->Phi[iq*nn+n]*xfac;
	      }
	//EV_U[n*sr+l] -= FluxWeights[k]*T[n*sr2+k*sr+l];
	
	/* Also include linearization of stabilization term */
	if ((StabVisc != NULL) && (!StabData->SkipDiffStab)){
	  // EV_U{m,l} += FluxWeights{k}*eta*ih*N{i,q}*AwStab{i,q,k}*wn{i,q}*StabVisc_U{m,l}*xfac
	  // Note, Aw{q,k} = N{i,q}*AwStab{i,q,k}*wn{i,q}
	  for (iq=0, val=0.; iq<nq; iq++){
	    for (k=0, ii=iq*sr, t=0.; k<sr; k++){
	      xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
	      t += FluxWeights[k]*DData->Aw[ii+k] * xfac;
	    }
	    val += t*eta*ih;
	  }
	  for (i=0; i<nn*sr; i++) EV_U[i] += val*StabVisc_U[i];
	}
      }
    }
    
    // Coupling is mu*eta*ih*|wn|, where mu*|wn| is already stored in CF
    if (CF != NULL) for (iq=0;iq<nq; iq++) CF[iq] *= eta*ih;

  }
  else if (DiffDisc == xfe_DiffDiscBR2){
    
    /* BR2 = Second form of Bassi & Rebay discretization

       Stabilization delta terms:

       Qn{q,k}   -= eta*dn{q,k}
       dn{q,k}    = ResPhi{q,n}*Dn{i,n,k}*wn{i,q}
       Dn{i,n,k}  = ResiM{n,m} * ResPhi{q,m}*ABdun{i,q,k}
       Dn_UI{i,n,k;o,a} = ResiM{n,m} * ResPhi{q,m}*(AB_uIdun{i,q,k;a} + AB{i,j,q,k,l}*du_uI{q,l;a}*wn{q,j})*Phi{q,o}
                                                    \_____________ ABdun_uI{i,q,k;a} ____________________/
       Since,
       ER{n,k} += Phi{n,q}*eta*ResPhi{q,m}*Dn{i,m,k}*wn{i,q} = S{i,n,m}*Dn{i,m,k}
       then
       ER_UI{n,k;o,a} += S{i,n,m}*Dn_UI{i,m,k;o,a}

       N{i,q} = normalized wn{i,q}
       Note, in using the mass matrix, iMM must be multiplied by fac

       Note, ResPhi and ResiM refer to quantities computed at a
       possibly different order than the state approximation.  This is
       useful for error estimation when the order of the residual
       needs to remain constant.

       GCL linearization:
       Dn_G{i,n,k;o} = ResiM{n,m}*ResPhi{q,m}*( 
                       -ABdun{i,q,k}/gb{q} + AB_uIdun{i,q,k;a}*(-uI{q,a}/gb{q})
                       + AB{i,j,q,k,l}*du_gb{q,l}*wn{q,j} ) * Phi{q,o}
       where  du_gb{q,l} = du_uB{q,l;a}*uB_gb{q,a;}
              uB_gb{q,a;} = (uB{q,a} - uB_uI{q,a,l}*uI{q,l} )/gb{q} 

       old calculations ... (not quite right)
                       + A_uBdun{i,q,k,a} * uB_gb{q,a;}
                       + AB{i,j,q,k,l}*du_uB{q,l,a}*uB_gb{q,a;}*wn{q,j})*Phi{q,o}
                     = ResiM{n,m}*ResPhi{q,m}*(-ABdun{i,q,k}/gb{q} 
                       + (A_uBdun{i,q,k,a} + AB{i,j,q,k,l}*du_uB{q,l,a}*wn{q,j})
                         * (uB{q,a} - uB_uI{q,a,l}*uI{q,l} )/gb{q} * Phi{q,o}
                     = ResiM{n,m}*ResPhi{q,m}*(-ABdun{i,q,k}/gb{q} + ABdun_uI{i,q,k;a}*(-uI{q,a}/gb{q}))*Phi{q,o}
                                                                     \_defined above_/

    */


    
    /* First add existing Qn to ER, since we will be dealing
       directly with ER.  Note that derivatives like Qn_uI
       Qn_guI, etc. have already been added above. */

    // Store flux components set by BC in Qn
    if (AFlag){
      ierr = xf_Error(xf_SetPrescribedQBC(SetIndex, nq, sr, dim, DData->guB,
					  DData->guB_guI, DData->Qn, NULL, NULL,
					  NULL, NULL, xfe_False));
      if (ierr != xf_OK) return ierr;
    }

    /* Add Q term to R:
       ER{n,k} -= Phi{n,q}*Qn{q,k},  sum over q
    */
    if (ER != NULL)
      xf_MTxM_Sub(PhiData->Phi, DData->Qn, nn, nq, sr, ER);

    /* Add Q term to Value:
       Value -= FluxWeights{k}*Qn{q,k}[*xglob{q,FluxMoments{k}}], sum over q,k
    */
    if (Value != NULL)
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++){
	  xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
	  (*Value) -= FluxWeights[k]*DData->Qn[iq*sr+k]*xfac;
	}

    /* Now can deal with delta terms.  First pull off inverse mass
       matrix, iMM. This is at the residual order. */
    ierr = xf_Error(xf_ElemInvMassMatrix(All, egrp, elem, Basis, ResOrder, 
                                         NULL, NULL, &iMM, &fac));
    if (ierr != xf_OK) return ierr;
    
    /* Calculate delta coefficients Dn (note use of fac for mass matrix scaling) */
    for (i=0; i<dim; i++){
      xf_MTxM_Set(ResPhiData->Phi, DData->ABdun+i*nq*sr, Rnn, nq, sr, T);
      xf_MxM_Set(iMM, T, Rnn, Rnn, sr, DData->Dn+i*Rnn*sr);
    }
    for (k=0; k<dim*Rnn*sr; k++) DData->Dn[k] *= fac;

    // delta is zero for flux components prescribed by BC -> zero out appropriate cols of Dn
    if (AFlag){
      ierr = xf_Error(xf_SetPrescribedQBC(SetIndex, dim*Rnn, sr, dim, NULL, NULL,
					  DData->Dn, NULL, NULL, NULL, NULL, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    
    // Determine eta
    eta = ViscStabFactor * 0.5*((real) Mesh->ElemGroup[egrp].nFace[elem]);
    // ih used to set coupling
    ih = FaceArea/ElemVol;

    /* Calculate:
         S{i,n,m} = eta*Phi{n,q}*ResPhi{q,m}*wn{i,q}
    */
    n2 = nn*Rnn;
    for (i=0; i<dim; i++){
      xf_ColcMult_Set(PhiData->Phi, wn+i, nq,  nn, dim, eta, T);   // row
      xf_MTxM_Set(T, ResPhiData->Phi, nn, nq, Rnn, DData->S+i*n2); // col
    }    

    /* Add to Qn (so that the returned Qn is actually a valid flux):
         Qn{q,k}   -= eta*dn{q,k}
         dn{q,k}    = ResPhi{q,n}*Dn{i,n,k}*wn{i,q}
    */
    for (i=0; i<dim; i++){
      // T{q,k}  = ResPhi{q,n}*Dn{i,n,k}
      xf_MxM_Set(ResPhiData->Phi, DData->Dn+i*Rnn*sr, nq, Rnn, sr, T);
      // Qn{q,k} += -etaT{q,k}*wn{i,q}
      xf_ColcMult_Add(T, wn+i, nq, sr, dim, -eta, DData->Qn);
    } 

    /* Add to ER:
         ER{n,k} += S{i,n,m}*Dn{i,m,k}
    */
    if (ER != NULL) 
      for (i=0; i<dim; i++)
	xf_MxM_Add(DData->S+i*n2, DData->Dn+i*Rnn*sr, nn, Rnn, sr, ER);
    

    /* Add to Value:
       Value += FluxWeights{k}*eta*dn{q,k}[*xglob{q,FluxMoments{k}}], sum over q,k
       Value += Vw{i,n,k}*Dn{i,n,k}  (restatement using more convenient form)
    */
    if ((Value != NULL) || (EV_U != NULL)){
      // Vw{i,n,k} = FluxWeights{k}*eta*Phi{q,n}*wn{i,q}[*xglob{q,FluxMoments{k}}]
      for (i=0; i<dim; i++){
	// construct T{q,k} = eta*x{q,k}*wn{i,q}*FluxWeights{k}
	for (iq=0; iq<nq; iq++)
	  for (k=0; k<sr; k++){
	    xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
	    T[iq*sr + k] = wn[iq*dim+i]*xfac*eta*FluxWeights[k];
	  }
	// set Vw{i,n,k} = ResPhi{n,q}*T{q,k}
	xf_MTxM_Set(PhiData->Phi, T, Rnn, nq, sr, DData->Vw+i*Rnn*sr);
      }
    }

    if (Value != NULL){
      // Value += Vw{i,n,k}*Dn{i,n,k} 
      xf_DotProduct(DData->Vw, DData->Dn, dim*Rnn*sr, &val);
      (*Value) += val;
    }

    
    /* Contribution to derivatives:
         P{i,n,q} = S{i,n,m}*iMM{m,o}*Phi{q,o}
         ER_U{n,k;o,a} += P{i,n,q}*(AB_uIdun{i,q,k,a} + AB{i,j,q,k,l}*du_uI{q,l,a}*wn{q,j})*Phi{q,o}
	 EV_U{o,a} += Vw{i,n,k}*ResiMM{n,m} * ResPhi{q,m}* ABdun_uI{i,q,k;a} * Phi{q,o}
         EV_G{o}   += Vw{i,n,k}*ResiMM{n,m} * ResPhi{q,m}*
                      (-ABdun{i,q,k} - ABdun_uI{i,q,k;a}*uI{q,a})/gb{q} * Phi{q,o}
         EV_G{o}   += Vw{i,n,k}*ResiMM{n,m} * ResPhi{q,m}*
                      ( -ABdun{i,q,k}/gb{q} + AB_uIdun{i,q,k;a}*(-uI{q,a}/gb{q})
                       +  AB{i,j,q,k,l}*du_gb{q,l}*wn{q,j}  ) * Phi{q,o}
         ABdun_uI{i,q,k;a} = AB_uIdun{i,q,k;a} + AB{i,j,q,k,l}*du_uI{q,l,a}*wn{q,j}
    */
    if ((ER_U != NULL) || (EV_U != NULL)){
      for (i=0; i<dim; i++){
	// T{n,o} = S{i,n,m}*ResiMM{m,o}
	xf_MxM_Set(DData->S+i*n2, iMM, nn, Rnn, Rnn, T);
	for (k=0; k<nn*Rnn; k++) T[k] *= fac; // fac is included here
	// P{i,q,n} = ResPhi{q,o}*T{o,n}
	xf_MxMT_Set(ResPhiData->Phi, T, nq, Rnn, nn, DData->P+i*nn*nq);
      }

      for (i=0; i<dim; i++){
	// Qn_uI{q,k,a} = (AB_uIdun{i,q,k,a} + AB{i,j,q,k,l}*du_uI{q,l,a}*wn{q,j})
	for (k=0; k<nq*sr2; k++) DData->Qn_uI[k] = 0.;
	for (j=0; j<dim; j++){
	  xf_nMxM_Set(nq, DData->AB+(i*dim+j)*nq*sr2, DData->du_uI, sr, sr, sr, T);
	  xf_ColcMult_Add(T, wn+j, nq, sr2, dim, 1.0, DData->Qn_uI);
	} // j
	if (!DData->ConstAB)
	  for (k=0; k<nq*sr2; k++) DData->Qn_uI[k] += DData->AB_uIdun[sr2*nq*i+k];      

        if (EV_G != NULL){
          /* Qn_gb{q,k;} = ( -ABdun{i,q,k}/gb{q} + AB_uIdun{i,q,k;a}*(-uI{q,a}/gb{q}) 
                             +  AB{i,j,q,k,l}*du_gb{q,l}*wn{q,j}  ) */
          for (k=0; k<nq*sr; k++) DData->Qn_gb[k] = 0.;
          // First, -ABdun{i,q,k}/gb{q}
          for (iq=0; iq<nq; iq++)
            for (k=0; k<sr; k++)
              DData->Qn_gb[iq*sr+k] = -DData->ABdun[i*nq*sr+iq*sr+k]/MD->gb[iq];
          // Second, AB_uIdun{i,q,k;a}*(-uI{q,a}/gb{q}) 
          if (!DData->ConstAB){
            for (iq=0; iq<nq; iq++){
              xf_MxV_Set(DData->AB_uIdun + sr2*nq*i + iq*sr2, uI+iq*sr, sr, sr, T);
              for (k=0; k<sr; k++) DData->Qn_gb[iq*sr+k] -= T[k]/MD->gb[iq];
            }
          }
          // Third, AB{i,j,q,k,l}*du_gb{q,l}*wn{q,j}
          for (j=0; j<dim; j++){
            //for (iq=0; iq<nq; iq++)
            //  xf_MxV_Set(DData->AB+(i*dim+j)*nq*sr2+iq*sr2, DData->du_gb+iq*sr, sr, sr, T+iq*sr);
            xf_nMxM_Set(nq, DData->AB+(i*dim+j)*nq*sr2, DData->du_gb, sr, sr, 1, T);
            //for (iq=0; iq<nq*sr2; iq++) xf_printf("AB[%d]=%.10E\n", iq, DData->AB[(i*dim+j)*nq*sr2+iq]);
            xf_ColcMult_Add(T, wn+j, nq, sr, dim, 1.0, DData->Qn_gb);
          } // j
          
          /* // Qn_gb{q,k;} = (-ABdun{i,q,k} - Qn_uI{q,k;a}*uI{q,a})/gb{q} */
          /* for (iq=0; iq<nq; iq++){ */
          /*   xf_MxV_Set(DData->Qn_uI+iq*sr2, uI+iq*sr, sr, sr, DData->Qn_gb+iq*sr); */
          /*   for (k=0; k<sr; k++) */
          /*     DData->Qn_gb[iq*sr+k] = -(DData->ABdun[i*nq*sr+iq*sr+k] + DData->Qn_gb[iq*sr+k])/MD->gb[iq]; */
          /* } */
          
        }

	// Set Aw{q,k} = AwStab{i,q,k}*Phi{q,k} (for stabilization, before zeroing out Aw)
	if (StabData != NULL)
	  xf_ColMult_Set(DData->AwStab+i*nq*sr, StabData->StabPhi, nq, sr, 1, DData->Aw);

	// before adding to ER_U, zero out derivatives for prescribed flux entries
	// entries in Aw{q,k} are zeroed out as well
	if (AFlag){
	  ierr = xf_Error(xf_SetPrescribedQBC(SetIndex, nq, sr, dim, DData->guB,
					      DData->guB_guI, NULL, DData->Aw, NULL,
					      DData->Qn_uI, NULL, xfe_True));
	  if (ierr != xf_OK) return ierr;
          if (EV_G != NULL){ // zero out GCL-specific entries too
            ierr = xf_Error(xf_SetPrescribedQBC(SetIndex, nq, sr, dim, NULL,
                                                NULL, NULL, DData->Qn_gb,
                                                NULL, NULL, NULL, xfe_True));
            if (ierr != xf_OK) return ierr;
          }
	}

	/*
	  Add to ER_U:
	  ER_U{n,k;o,a} += P{i,n,q}*[Qn_uI{q,k;a}*Phi{q,o}]
	*/
	if (ER_U != NULL){
	  for (n=0; n<nn; n++){
	    xf_ColMult_Set(PhiData->Phi, DData->P+i*nn*nq+n, nq, nn, nn, T);
	    xf_MTxM_Add(T, DData->Qn_uI, nn, nq, sr2, ER_U+n*nn*sr2);
	  } // n
	
	  /* Also include linearization of stabilization term */
	  if ((StabVisc != NULL) && (!StabData->SkipDiffStab) && (StabVisc_U != NULL)){
	    // ER_U{n,k;m,a} += P{i,n,q}*AwStab{i,q,k}*StabPhi{q}*StabVisc_U{m,a}
	    //                = P{i,n,q}*Aw{q,k}                 *StabVisc_U{m,a}
	    xf_MTxM_Set(DData->P+i*nn*nq, DData->Aw, nn, nq, sr, T);  // T is nn x sr
	    xf_BlockOutProd_Add(T, StabVisc_U, nn, sr, nn, ER_U);
	  }
	}
	
	/*
	  Add to EV_U:
	  EV_U{o,a} += Vw{i,n,k}*ResiMM{n,m} * ResPhi{q,m}* Qn_uI{q,k;a} * Phi{q,o}
          EV_G{o}   += Vw{i,n,k}*ResiMM{n,m} * ResPhi{q,m}* Qn_gb{q,k;} * Phi{q,o}
	*/
	if (EV_U != NULL){
	  // T{q,n} = ResPhi{q,m}*iMM{m,n}
	  xf_MxMT_Set(ResPhiData->Phi, iMM, nq, Rnn, Rnn, T); // T is nq x Rnn
	  // Aw{q,k} = T{q,n}*Vw{i,n,k}
	  xf_MxM_Set(T, DData->Vw+i*Rnn*sr, nq, Rnn, sr, DData->Aw);
	  for (k=0; k<nq*sr; k++) DData->Aw[k] *= fac; // fac is included here
	  // T{q,a} = Aw{q,k}*Qn_uI{q,k;a}
	  for (iq=0; iq<nq; iq++)
	    xf_MTxV(DData->Qn_uI+iq*sr2, DData->Aw+iq*sr, sr, sr, xfe_Set, T+iq*sr);
	  // EV_U{o,a} += Phi{o,q}*T{q,a}
	  xf_MTxM_Add(PhiData->Phi, T, nn, nq, sr, EV_U);

          if (EV_G != NULL){ // GCL-specific
            // T{q} = Aw{q,k}*Qn_gb{q,k}
            for (iq=0; iq<nq; iq++)
              xf_DotProduct(DData->Aw+iq*sr, DData->Qn_gb+iq*sr, sr, T+iq);
            // EV_G{o} += Phi{o,q}*T{q}
            xf_MTxM_Add(PhiData->Phi, T, nn, nq, 1, EV_G); 
          }

	
	  /* Also include linearization of stabilization term */
	  if ((StabVisc != NULL) && (!StabData->SkipDiffStab)){
	    // EV_U{o,a} += Vw{i,n,k}*ResiMM{n,m}*ResPhi{q,m}*AwStab{i,q,k}*StabVisc_U{o,a}
	    //            = Aw{q,k}                          *AwStab{i,q,k}*StabVisc_U{o,a} 
	    xf_DotProduct(DData->Aw, DData->AwStab+i*nq*sr, nq*sr, &val);
	    for (k=0; k<nn*sr; k++) EV_U[k] += val*StabVisc_U[k];
	  }
	}
  
      } // i
    }
       
    // Coupling is mu*eta*ih*|wn|, where mu*|wn| is already stored in CF
    if (CF != NULL) for (iq=0;iq<nq; iq++) CF[iq] *= eta*ih;
  }
  else return xf_Error(xf_NOT_SUPPORTED); // only support IP or BR2

  xf_Release( (void *) SetIndex);

  return xf_OK;
}

