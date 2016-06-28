/*------------------------------------------------------------------*
 * This code is heritaged from XFLOW, but have totally different    *
 * functionality. Author: Yu Lv Email:lvyu@umich.edu                *
 * Knowledgements: Kfid for orignal contribution
 *------------------------------------------------------------------*/

/*
  FILE:  xf_Residual.c

  This file contains the residual-calculation functions.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_Math.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Quad.h"
#include "xf_Basis.h"
#include "xf_MeshTools.h"
#include "xf_EqnSetHook.h"
#include "xf_Line.h"
#include "xf_ResidualDiff.h"
#include "xf_ResidualStab.h"
#include "xf_MeshDistance.h"
#include "xf_All.h"
#include "xf_Solver.h"
#include "xf_MeshMotion.h"
#include "xf_MeshMotionGCL.h"
#include "xf_LinQueue.h"
#include "xf_LinearSolver.h"
#include "xfYu_Model.h"
#include "Yu_DiffusDiscret.c"
//#include "xfYu_AdaptSolver.h"

static Yu_Model *Model;


/* Structure for storing  "static" data in elem residual calculation.
   Static data refers to any data that can be passed repeatedly into
   the function to avoid excessive reallocations and calculations. 
typedef struct
{
  int             iConv[2], iDiff[2], iSource[2];
  xf_QuadData     *QuadData;
  xf_BasisData    *PhiData;
  xf_JacobianData *JData;
  xf_BasisData    *GeomPhiData;
  real            *xglob;
  real            *T, *wq, *u, *gu;
  real            *F, *F_u, *Aw, *A, *A_uw, *AwStab;
  real            *StabVisc;
  real            *ResMetric;
  real            *S, *S_u, *S_gu;
  int             pnq, nqConv, nqDiff, Tsize;
  int             *IParam;
  real            *RParam;
  int             nAuxU;
  xf_BasisData    **AuxPhiData;
  real            **Auxu;  
  xf_Vector       **AuxU;
  xf_Vector       *EG;
  xf_Vector       *EM;
  real            Time;
  xf_MotionData   *MD;
  real            *gum; 
  xf_LinQueueData  LinQ[1];
}
xf_StaticDataElem;
*/
typedef struct
{
    xf_QuadData     *QuadData;
    xf_BasisData    *PhiData;
    xf_JacobianData *JData;
    xf_BasisData    *GeomPhiData;
    real            *xglob;
    real            *wq, *u, *F, *S;
    real            *gu, *Aw;
    int             pnq;
    int             nqConv;
    int             nqDiff;
    xf_Vector       *EG;
    real            Time;
}
xf_StaticDataElem;

/* Structure for storing "static" data in iface residual calculation. 
typedef struct
{
  int             iConv[2], iDiff[2], iSource[2];
  xf_QuadData     *QuadData;
  xf_BasisData    *PhiDataL, *PhiDataR;
  xf_BasisTable   *PhiTable;
  xf_NormalData   *NData;
  xf_JacobianData *JDataL, *JDataR;
  xf_BasisData    *GeomPhiData;
  xf_BasisTable   *GeomPhiTable;
  real            *xglob;
  real            *wn;
  real            *xelemL, *xelemR;
  real            *uL, *uR;
  real            *guL, *guR;
  real            *F, *F_uL, *F_uR;
  real            *CF;
  xf_DiffJumpData *DData;
  int             pnq, nqConv, nqDiff, nnDiff;
  int             *IParam;
  real            *RParam;
  int             nAuxU;
  xf_BasisData    **AuxPhiDataL, **AuxPhiDataR;
  real            **AuxuL, **AuxuR;  
  xf_Vector       **AuxU;
  xf_Vector       *EG;
  xf_Vector       *EM;
  real            *StabViscL,  *StabViscR;
  real            *StabPhiL,   *StabPhiR;
  real            *ResMetricL, *ResMetricR;
  real            Time;
  xf_MotionData   *MDL, *MDR;
  xf_LinQueueData  LinQ[4];
  xf_BasisData    *ResPhiDataL, *ResPhiDataR;
  xf_BasisTable   *ResPhiTable;
}
xf_StaticDataIFace;
*/

typedef struct
{
    xf_QuadData     *QuadData;
    xf_BasisData    *PhiDataL, *PhiDataR;
    xf_BasisTable   *PhiTable;
    xf_NormalData   *NData;
    xf_JacobianData *JDataL, *JDataR;
    xf_BasisData    *GeomPhiData;
    xf_BasisTable   *GeomPhiTable;
    xf_DiffJumpData *DData;
    real            *xglob;
    real            *wn;
    real            *xelemL, *xelemR;
    real            *uL, *uR;
    real            *guL, *guR;
    real            *F;
    int             pnq, nqConv, nqDiff, nnDiff;
    xf_Vector       *EG;
    real            Time;
}
xf_StaticDataIFace;

/* Structure for storing "static" data in bface residual calculation. 
typedef struct
{
  int             iConv[2], iDiff[2], iSource[2];
  xf_QuadData     *QuadData;
  xf_BasisData    *PhiData;
  xf_BasisTable   *PhiTable;
  xf_NormalData   *NData;
  xf_JacobianData *JData;
  xf_BasisData    *GeomPhiData;
  xf_BasisTable   *GeomPhiTable;
  real            *xglob;
  real            *T, *T2, *wn;
  real            *xelem;
  real            *uI;
  real            *guI;
  real            *uB;
  real            *uB_uI;
  real            *F, *F_uI;
  real            *CF;
  xf_DiffBCData   *DData;
  int             pnq, nqConv, nqDiff, nnDiff, Tsize, T2size;
  int             *IParam;
  real            *RParam;
  int             nAuxU;
  xf_BasisData    **AuxPhiData;
  real            **Auxu;
  xf_Vector       **AuxU;
  xf_Vector       *EG;
  xf_Vector       *EM;
  real            *StabVisc;
  real            *StabPhi;
  real            *ResMetric;
  real            Time;
  xf_MotionData   *MD;
  xf_LinQueueData  LinQ[1];
  xf_BasisData    *ResPhiData;
  xf_BasisTable   *ResPhiTable;
}
xf_StaticDataBFace;
*/
typedef struct
{
    xf_QuadData     *QuadData;
    xf_BasisData    *PhiData;
    xf_BasisTable   *PhiTable;
    xf_NormalData   *NData;
    xf_JacobianData *JData;
    xf_BasisData    *GeomPhiData;
    xf_BasisTable   *GeomPhiTable;
    real            *xglob;
    real            *wn;
    real            *xelem;
    real            *uI, *guI;
    real            *uB;
    xf_DiffBCData   *DData;
    real            *F;
    int             pnq, nqConv, nqDiff, nnDiff;
    xf_Vector       *EG;
    real            Time;
}
xf_StaticDataBFace;

/******************************************************************/
//   FUNCTION Definition: xf_CheckRecoverable
static enum xfe_Bool
xf_CheckRecoverable(int InputError, int *ReturnError)
{
  if (InputError == xf_OK) return xfe_True;
  
  (*ReturnError) = InputError;

  if (InputError == xf_NON_PHYSICAL){
    return xfe_True;
  }
  else // not recoverable
    return xfe_False;
}

/******************************************************************/
//special data structure for detailed chemistry ODE solver
typedef struct
{
    xf_QuadData     *QuadData;
    xf_BasisData    *PhiData;
    xf_JacobianData *JData;
    xf_BasisData    *GeomPhiData;
    real            *xglob;
    real            *wq, *u, *F, *S;
    int             pnq;
    xf_Vector       *EG;
    real            Time;
}
xf_StaticDataElem_DetailChem;

/**********************************************************************/
//   FUNCTION Definition: xf_CreateStaticDataElem (for detail chemistry)
static int
xf_CreateStaticDataElem_DetailChem(xf_All *All, xf_SolverData *SolverData, 
                                   xf_StaticDataElem_DetailChem **pSD)
{
    /*
     PURPOSE: Initializes element static data
     
     INPUTS: 
     
     All : all structure to determine residual terms
     SolverData : solver data structure
     
     OUTPUTS: 
     
     (*pSD) : static data that is created and initialized
     
     RETURNS: Error code
     
     */
    int ierr, iAux;
    xf_StaticDataElem_DetailChem *SD;
    
    // allocate memory
    ierr = xf_Error(xf_Alloc( (void **) pSD, 1, sizeof(xf_StaticDataElem_DetailChem)));
    if (ierr != xf_OK) return ierr;
    SD = (*pSD);
    
    // initialize variables to NULL
    SD->QuadData    =  NULL;
    SD->PhiData     =  NULL;
    SD->JData       =  NULL;
    SD->wq	  =  NULL;
    SD->GeomPhiData =  NULL;
    SD->xglob       =  NULL;
    //SD->T	          =  NULL;
    SD->u	          =  NULL; 
    //SD->gu	  =  NULL; 
    SD->F	          =  NULL; 
    //SD->F_u	  =  NULL; 
    //SD->Aw	  =  NULL; 
    //SD->A	          =  NULL; 
    //SD->A_uw        =  NULL; 
    //SD->StabVisc    =  NULL; 
    //SD->ResMetric   =  NULL; 
    //SD->AwStab      =  NULL; 
    SD->S	          =  NULL; 
    //SD->S_u	  =  NULL; 
    //SD->S_gu	  =  NULL; 
    SD->pnq         =  -1; 
    //SD->nqConv      =  -1; 
    //SD->nqDiff      =  -1; 
    //SD->Tsize       =  -1;
    SD->EG          =  NULL;
    //SD->EM          =  NULL;
    //SD->gum         =  NULL;
    
    // determine Time
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &SD->Time));
    if (ierr != xf_OK) return ierr;
    
    // build eqnset-desired parameter lists for passing into functions
    //ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &SD->IParam, &SD->RParam, 
    //			       &SD->nAuxU, &SD->AuxU));
    //if (ierr != xf_OK) return ierr;
    
    // Allocate vector of basis data for interpolating auxiliary vectors
    /*
     ierr = xf_Error(xf_Alloc( (void **) &SD->AuxPhiData, SD->nAuxU, sizeof(xf_BasisData *)));
     if (ierr != xf_OK) return ierr;
     ierr = xf_Error(xf_Alloc( (void **) &SD->Auxu, SD->nAuxU, sizeof(real *)));
     if (ierr != xf_OK) return ierr;
     for (iAux=0; iAux<SD->nAuxU; iAux++){
     SD->AuxPhiData[iAux] = NULL;
     SD->Auxu[iAux]       = NULL;
     }
     
     if ( (SolverData != NULL) && (SolverData->StabRequired) ){
     // obtain elem geometry vector
     ierr = xf_Error(xf_FindElemGeom(All, &SD->EG));
     if (ierr != xf_OK) return ierr;
     // obtain elem metric vector
     ierr = xf_Error(xf_FindElemHMetric(All, xfe_True, &SD->EM));
     if (ierr != xf_OK) return ierr;
     }
     
     // initialize motion data
     SD->MD = NULL;
     if ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active)){
     ierr = xf_Error(xf_CreateMotionData(All, &SD->MD));
     if (ierr != xf_OK) return ierr;
     }
     
     // initialize data structure for residual linearization
     ierr = xf_Error(xf_InitLinQueue(SD->LinQ));
     if (ierr != xf_OK) return ierr;
     */
    
    // everything must go back with (*pSD)
    (*pSD) = SD;
    
    return xf_OK;
}



/***********************************************************************/
//   FUNCTION Definition: xf_DestroyStaticDataElem (for detail chemistry)
static int
xf_DestroyStaticDataElem_DetailChem(xf_StaticDataElem_DetailChem *SD)
{
    /*
     PURPOSE: Destroys static element structure
     
     INPUTS: 
     
     SD : static data
     
     OUTPUTS: None, SD and all of its contents are destroyed
     
     RETURNS: Error code
     
     */
    
    int ierr, iAux;
    
    if (SD == NULL) return xf_OK;
    
    // Only destroy QuadData if points are generic
    ierr = xf_Error(xf_DestroyGenericQuadData(SD->QuadData));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(SD->PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Geometry Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(SD->GeomPhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy geometry Jacobian Data */
    ierr = xf_Error(xf_DestroyJacobianData(SD->JData));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Auxiliary Vector Data */
    /*
     for (iAux=0; iAux<SD->nAuxU; iAux++){
     ierr = xf_Error(xf_DestroyBasisData(SD->AuxPhiData[iAux], xfe_True));
     if (ierr != xf_OK) return ierr;
     xf_Release( (void *) SD->Auxu[iAux]);
     }
     xf_Release( (void *) SD->AuxPhiData);
     xf_Release( (void *) SD->Auxu);
     xf_Release( (void *) SD->AuxU);
     */
    
    /* Destroy mesh motion data */
    //xf_DestroyMotionData(SD->MD);
    
    // Destroy residual linearization structure
    //xf_DestroyLinQueue(SD->LinQ);
    
    // Release all other memory
    //xf_Release( (void *) SD->IParam);
    //xf_Release( (void *) SD->RParam);
    xf_Release( (void *) SD->wq);
    xf_Release( (void *) SD->xglob);
    //xf_Release( (void *) SD->T);
    xf_Release( (void *) SD->u);
    //xf_Release( (void *) SD->gu);
    xf_Release( (void *) SD->F);
    //xf_Release( (void *) SD->F_u);
    //xf_Release( (void *) SD->Aw);
    //xf_Release( (void *) SD->A);
    //xf_Release( (void *) SD->A_uw);
    //xf_Release( (void *) SD->StabVisc);
    //xf_Release( (void *) SD->ResMetric);
    //xf_Release( (void *) SD->AwStab);
    xf_Release( (void *) SD->S);
    //xf_Release( (void *) SD->S_u);
    //xf_Release( (void *) SD->S_gu);
    //xf_Release( (void *) SD->gum);
    
    // Destroy self
    xf_Release( (void *) SD);
    
    return xf_OK;
}

/***********************************************************************/
//   FUNCTION Definition: xf_CalculateResidualElem (for detail chemistry)
static int
xf_CalculateResidualElem_DetailChem(xf_All *All, int egrp, int elem, xf_Vector *U, 
                                    real *ER, real *ER_U, real **ER_NU, 
                                    xf_StaticDataElem_DetailChem **pSD, 
                                    xf_SolverData *SolverData)
{
    /*
     PURPOSE: 
     
     Calculates the residual and residual Jacobian (if requested)
     associated with an element interior integration on (egrp, elem).
     
     INPUTS: 
     
     All: All structure
     egrp, elem : element in question
     U : state vector on all elements
     pSD : pointer to static data (optional, can pass in as NULL)
     useful for avoiding reallocations when calling multiple times
     SolverData : solver data structure
     
     OUTPUTS:
     
     ER : element residual (must be preallocated)
     ER_U : element residual Jacobian (optional, but if given must be preallocated)
     ER_NU : element jacobian w.r.t neighbor states (optional, but preallocated if given)
     
     RETURNS: Error code
     
     */
    int ierr, sr, sr2, iq, nq, nn, d, n, i, j;
    int dim, Order, QuadOrder, iAux, nAuxU;
    int ResOrder, ResidualOrderIncrement;
    int Tsize, pnq, nqConv, nqDiff;
    int  nConv,  nDiff,  nSource;
    int *iConv, *iDiff, *iSource;
    int *IParam;
    enum xfe_BasisType Basis;
    enum xfe_Bool QuadChanged, ConstA, SkipDiffusion, SkipDiffStab, StabRequired;
    enum xfe_Bool Need_gu = xfe_False, Nonzero_gu = xfe_False;
    enum xfe_Bool Need_Grad;
    enum xfe_Bool MotionOn;
    real *EU, *F, *F_u, *Aw, *A, *A_uw, *AwStab, Gamma;
    real *S, *S_u, *S_gu, **Auxu = NULL;
    real *RParam, *xq, *wq, *xglob, *u, *gu, *T, *gum = NULL, *w = NULL;
    real Time;
    real *StabVisc = NULL;
    xf_ResTerm *ResTerm;
    xf_QuadData *QuadData;
    xf_BasisData *PhiData, *GeomPhiData, **AuxPhiData;
    xf_JacobianData *JData;
    xf_StabData *StabData = NULL;
    xf_Vector *EG, **AuxU, *V, *GammaVec;
    xf_MotionData *MD = NULL;
    xf_LinQueueData *LinQ = NULL;
    xf_StaticDataElem_DetailChem *SD = NULL;
    xf_Mesh   *Mesh;
    xf_Data   *GammaDat;
    //xf_EqnSet *EqnSet;
    
    // General information
    Mesh    = All->Mesh;
    dim     = Mesh->Dim;  
    //EqnSet  = All->EqnSet;
    //sr      = EqnSet->StateRank;
    sr      = Model->nVars;
    sr2     = sr*sr;
    
    
    // Determine Basis and Order from the state, U
    Basis = U->Basis[egrp];
    Order = xf_InterpOrder(U, egrp, elem);
    
    // do we have source terms? quick check to see if we can skip p=0
    /*
     nSource = 0;
     for (i=0; i<EqnSet->ResTerms->nResTerm; i++)
     if ((EqnSet->ResTerms->ResTerm[i].Type == xfe_ResTermSource) &&
     (EqnSet->ResTerms->ResTerm[i].Active)) 
     nSource++;
     // skip Order = 0 elements if no source
     if ((Order == 0) && (nSource == 0)) return xf_OK;
     */
    //for classic finite volume (no need to compute this term
    //if (Order == 0) return xf_OK;
    
    // Create and initialize static data if not passed in
    if ((pSD == NULL) || ((*pSD) == NULL)){
        ierr = xf_Error(xf_CreateStaticDataElem_DetailChem(All, SolverData, (pSD != NULL) ? pSD : &SD));
        if (ierr != xf_OK) return ierr;
    }
    if (pSD != NULL) SD = (*pSD);
    
    /*
     ResTerm = EqnSet->ResTerms->ResTerm;
     
     // Are we doing mesh motion?
     MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));
     
     // Pull off variables from StaticData  
     iConv   = SD->iConv;
     iDiff   = SD->iDiff;
     iSource = SD->iSource;
     nConv   =   iConv[1]-  iConv[0];
     nDiff   =   iDiff[1]-  iDiff[0]; 
     nSource = iSource[1]-iSource[0];
     */
    
    //find Gamma Vector
    ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
    if(ierr == xf_NOT_FOUND)
    {
        xf_printf("Cannot find heat capacity ratio...\n");
        return ierr;
    }
    else
        GammaVec = (xf_Vector *) GammaDat->Data;
    //if found, directly get the gamma value at this element
    Gamma = GammaVec->GenArray[egrp].rValue[elem][0];
    
    QuadData    = SD->QuadData;
    PhiData     = SD->PhiData;
    JData       = SD->JData;
    wq          = SD->wq;
    GeomPhiData = SD->GeomPhiData;
    xglob       = SD->xglob;
    //T           = SD->T;
    u           = SD->u;
    //gu          = SD->gu;
    F           = SD->F;
    //F_u         = SD->F_u;
    //Aw          = SD->Aw;
    //A           = SD->A;
    //A_uw        = SD->A_uw;
    //AwStab      = SD->AwStab;
    S           = SD->S;
    //S_u         = SD->S_u;
    //S_gu        = SD->S_gu;
    pnq         = SD->pnq;        // previous # quad points
    //nqConv      = SD->nqConv;     // for allocation purposes
    //nqDiff      = SD->nqDiff;
    //Tsize       = SD->Tsize;      // for temporary matrix allocation
    //IParam      = SD->IParam;
    //RParam      = SD->RParam;
    //nAuxU       = SD->nAuxU;      // number of auxiliary vectors
    //AuxU        = SD->AuxU;       // auxiliary vectors
    //AuxPhiData  = SD->AuxPhiData; // auxiliary vector basis data
    //Auxu        = SD->Auxu;       // auxiliary vector data
    EG          = SD->EG;         // element geometry vector
    Time        = SD->Time;       // simulation time
    //LinQ        = SD->LinQ;       // residual linearization structure
    //if (MotionOn){
    //  MD          = SD->MD;       // mesh motion data
    //  gum         = SD->gum;      // physical gradient for mesh motion
    //}
    
    // will we need the gradient of u?
    //if ((nDiff > 0) || (nSource > 0)) Need_gu = xfe_True;
    
    // do we need a linearization
    //Need_Grad = ((ER_U != NULL) || (ER_NU != NULL));
    
    // do we need stabilization
    //StabRequired = ( (SolverData != NULL) && (SolverData->StabRequired) );
    
    // Residual order increase
    //ResidualOrderIncrement = ((SolverData == NULL) ? 0 : SolverData->ResidualOrderIncrement);
    
    // Residual order
    //ResOrder = (ResidualOrderIncrement == 0) ? Order : max(Order + ResidualOrderIncrement, 0);
    
    // determine required integration order 
    
    ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order, &QuadOrder));
    if (ierr != xf_OK) return ierr;
    
    if(!Model->TwistFlag)
        QuadOrder -= Mesh->Dim;
    
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
    
    nn = PhiData->nn; // number of interpolation nodes
    
    // re-allocate data if quad points increased
    if (nq > pnq){
        ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
        if (ierr != xf_OK) return ierr;
        /*
         if (Need_gu){
         ierr = xf_Error(xf_ReAlloc( (void **) &gu, nq*sr*dim, sizeof(real)));
         if (ierr != xf_OK) return ierr;
         }
         if (StabRequired){
         ierr = xf_Error(xf_ReAlloc( (void **) &SD->StabVisc, nq, sizeof(real)));
         if (ierr != xf_OK) return ierr;
         ierr = xf_Error(xf_ReAlloc( (void **) &SD->ResMetric, nq*dim*dim, sizeof(real)));
         if (ierr != xf_OK) return ierr;
         }
         */
    }
    
    // obtain global coords of quad points
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
                                    nq, xq, xglob));
    if (ierr != xf_OK) return ierr;
    
    // obtain transformation map if doing mesh motion
    /*if (MotionOn){
     ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, All->Mesh->Motion, 
     nq, dim, Time, xglob, MD));
     if (ierr != xf_OK) return ierr;
     }*/
    
    EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
    
    // interpolate state and gradient at quad points
    xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u);      
    /*
     if (Need_gu)
     for (d=0; d<dim; d++)
     xf_MxM_Set(PhiData->gPhi+nn*nq*d, EU, nq, nn, sr, gu+nq*sr*d);
     
     // interpolate any auxiliary vectors, using AuxPhiData
     for (iAux=0; iAux<nAuxU; iAux++){
     V  = AuxU[iAux];
     if (nq > pnq){ // reallocate Auxu[iAux] if necessary
     ierr = xf_Error(xf_ReAlloc( (void **) Auxu+iAux, nq*V->StateRank, sizeof(real)));
     if (ierr != xf_OK) return ierr;
     }
     ierr = xf_Error(xf_EvalBasis(V->Basis[egrp], xf_InterpOrder(V,egrp,elem), QuadChanged, 
     nq, xq, xfb_Phi, AuxPhiData+iAux));
     if (ierr != xf_OK) return ierr;
     xf_MxM_Set(AuxPhiData[iAux]->Phi, V->GenArray[egrp].rValue[elem], nq, 
     AuxPhiData[iAux]->nn, V->StateRank, Auxu[iAux]);
     }*/
    
    
    // form detJ-multiplied quad weight vector, wq
    for (iq=0; iq<nq; iq++) 
        wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
    
    // pull off data required for stabilization
    /*
     if ( StabRequired ){
     StabData = &(SolverData->StabData);
     ierr = xf_Error(xf_CalculateStabViscElem(All, egrp, elem, nq, xq, SD->EM->GenArray[egrp].rValue[elem], 
     StabData, SD->StabVisc, SD->ResMetric, &SkipDiffStab));
     if (ierr != xf_OK) return ierr;
     StabVisc = SD->StabVisc;
     }
     else {
     StabVisc = NULL;
     SkipDiffStab = xfe_False;
     }*/
    
    /* Clear residual linearization queue */
    //if (ER_U != NULL) xf_ClearLinQueue(LinQ);
     
    /*--------------*/
    /* SOURCE TERMS */
    /*--------------*/
    
    // if (nSource > 0){
    // source term for sponge layer

    if(Model->Sponge != NULL){
    
        if (nq > pnq){	// realloc S, S_u if necessary
            ierr = xf_Error(xf_ReAlloc( (void **) &S, sr*nq, sizeof(real)));
            if (ierr != xf_OK) return ierr;
        }

        ierr = xf_Error(SpongeSource(nq, sr, dim, xglob, u, S, Gamma,  Model->dt_size, Time,  
                                     Model->Sponge->GenArray[egrp].rValue[elem][0], Model->SpongeParam));
        if (ierr != xf_OK) return ierr;

        // multiply S by quad weights*J
        xf_ColMult(S, wq, nq, sr, 1); // S is modified here

        xf_MTxM_Add(PhiData->Phi, S, nn, nq, sr, ER);


    }

    if(Model->ChemSource && Model->DetailChem) {
       
        if (nq > pnq){	// realloc S, S_u if necessary
            ierr = xf_Error(xf_ReAlloc( (void **) &S, sr*nq, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            /*    if (ER_U != NULL){
             ierr = xf_Error(xf_ReAlloc( (void **) &S_u, sr*sr*nq, sizeof(real)));
             if (ierr != xf_OK) return ierr;
             if (Need_gu){
             ierr = xf_Error(xf_ReAlloc( (void **) &S_gu, sr*sr*nq*dim, sizeof(real)));
             if (ierr != xf_OK) return ierr;
             }
             }
             */
        }
        
        // realloc temporary matrix for ER_U construction, if necessary
        /* if ((ER_U != NULL)  && (nn*nq > Tsize)){
         Tsize = nn*nq;
         ierr = xf_Error(xf_ReAlloc( (void **) &T, Tsize, sizeof(real)));
         if (ierr != xf_OK) return ierr;
         }
         
         // transform state to physical
         if (MotionOn) xf_ModMotionPreEqnCall(nq, dim, sr, MD, u, NULL);
         */
        // calculate S [and S_u] at quad points
        //ierr = xf_Error(xf_EqnSetSourceS(EqnSet, ResTerm+iSource[0], nSource, IParam, RParam, 
        //				     nq, u, gu, Auxu, xglob, &Time, S, S_u, S_gu, &Nonzero_gu));
        ierr = xf_Error(ChemicalSource(nq, sr, dim, u, xglob, S, Gamma, Model->DetailChem, Model->dt_size));
        if (ierr != xf_OK) return ierr;
        
        // transform state to ref
        //if (MotionOn) xf_ModMotionPostEqnCall(nq, dim, sr, MD, u, NULL);
        
        
        // multiply quad weights by g if MeshMotion is on (last time quad weights used)
        //if (MotionOn) for (iq=0; iq<nq; iq++) wq[iq] *= MD->g[iq];
        
        // multiply S by quad weights*J
        xf_ColMult(S, wq, nq, sr, 1); // S is modified here
        
        // Add to R:
        //   ER{n,k} += sum_q Phi{q,n}^T * S{q,k}*wq{q}
        
        xf_MTxM_Add(PhiData->Phi, S, nn, nq, sr, ER);
        
        // Add to ER_U:
        //   ER_U{n,k;m,a} += sum_q Phi{q,n}^T * S_u{q,k;a}*Phi{q,m}
        /* 
         if (ER_U != NULL){
         
         ierr = xf_Error(xf_AddToLinQueue(S_u, xfe_LinQTerm_PhiPhi, nq, dim, 
         sr2, -1, wq, 1, 1.0, LinQ));
         if (ierr != xf_OK) return ierr;
         
         // Nonzero_gu indicates that S depends on the gradient of u, gu
         //ER_U{n,k;m,a} += sum_i sum_q Phi{q,n}^T * S_gu{i,q,k;a}*gPhi{i,q,m}
         //(note, this is a dual-inconsistent discretization at this point)
         
         if (Nonzero_gu){
         ierr = xf_Error(xf_AddToLinQueue(S_gu, xfe_LinQTerm_PhiGPhi, nq, dim,
         sr2, -1, wq, 1, 1.0, LinQ));
         if (ierr != xf_OK) return ierr;
         }
         }
         */ 
    } // end if nSource > 0
   // else
   //     return xf_Error(xf_NOT_SUPPORTED);  //not consistent
    
    // apply linearizations queued up
    //if (ER_U != NULL){
    //  ierr = xf_Error(xf_ApplyLinQueue(LinQ, PhiData, PhiData, sr2, ER_U));
    //  if (ierr != xf_OK) return ierr;
    //}
    
    pnq = nq; // set previous quad point # for next element
    
    // Store possibly-altered or resized data back in StaticData
    SD->QuadData    =  QuadData;
    SD->PhiData     =  PhiData;  
    SD->JData       =  JData;
    SD->wq          =  wq;
    SD->GeomPhiData =  GeomPhiData;
    SD->xglob       =  xglob;
    //SD->T	          =  T;
    SD->u	          =  u;
    //SD->gu	  =  gu;
    SD->F	          =  F;
    //SD->F_u	  =  F_u;
    //SD->Aw	  =  Aw;
    //SD->A	          =  A;
    //SD->A_uw        =  A_uw;
    //SD->AwStab      =  AwStab;
    SD->S	          =  S;
    //SD->S_u	  =  S_u;
    //SD->S_gu        =  S_gu;
    SD->pnq         =  pnq;    
    //SD->nqConv      =  nqConv; 
    //SD->nqDiff      =  nqDiff;
    //SD->Tsize       =  Tsize;  
    //SD->AuxPhiData  =  AuxPhiData; 
    //SD->Auxu        =  Auxu;
    //SD->LinQ[0]     =  *LinQ; // not really necessary
    //if (MotionOn){
    //  SD->MD        =  MD;
    //  SD->gum       =  gum;
    //}
    
    if (pSD == NULL){
        // Delete StaticData that we just created
        ierr = xf_Error(xf_DestroyStaticDataElem_DetailChem(SD));
        if (ierr != xf_OK) return ierr;
    }
    
    return xf_OK;
    
}




/************************************************************************/
//   FUNCTION Definition: xf_CalculateResidualElems (for detail chemistry)
//   Time integration for detailed chemistry using implicit ODE solver
int
xf_CalculateResidualElems_DetailChem(xf_All *All, Yu_Model *pModel, xf_Vector *U, 
                                     xf_Vector *R, xf_JacobianMatrix *R_U, 
                                     xf_SolverData *SolverData)
{
    /*
     PURPOSE: 
     
     Calculates the residual and residual Jacobian (if requested)
     associated with an element interior integration on all elements.
     
     INPUTS: 
     
     All: All structure
     U : state vector on all elements
     SolverData : solver data structure
     
     OUTPUTS:
     
     R : residual
     R_U : element Jacobian (optional, or may be given but Value may not exist)
     
     RETURNS: Error code
     
     */
    
    
    int ierr, egrp, elem;
    real *ER, *ER_U;
    real **ER_NU = NULL;
    xf_StaticDataElem_DetailChem *StaticData = NULL;
    xf_Mesh   *Mesh;
    
    //plug in model parameter from Yu
    Model = pModel;
    Mesh = All->Mesh;
    
    // initialize R to zero
    ierr = xf_Error(xf_SetZeroVector(R));
    if (ierr != xf_OK) return ierr;
    
    
    // loop over element groups
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
        
        // loop over elements
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
            
            /* Prepare pointers to residual and Jacobian */
            ER = R->GenArray[egrp].rValue[elem];
            if ((R_U == NULL) || (R_U->Value == NULL)){
                ER_U  = NULL;
                ER_NU = NULL;
            }
            else{
                ER_U  = R_U->Value[egrp][elem][0];
                ER_NU = R_U->Value[egrp][elem]+1;
            }
            
            /* Calculate residual on elem, passing in StaticData  */
            ierr = xf_Error(xf_CalculateResidualElem_DetailChem(All, egrp, elem, U, ER, ER_U, 
                                                                ER_NU, &StaticData, SolverData));
            if (ierr != xf_OK) return ierr;
            
        } // elem
        
    } // egrp
    
    
    // Delete StaticData
    ierr = xf_Error(xf_DestroyStaticDataElem_DetailChem(StaticData));
    if (ierr != xf_OK) return ierr;
    
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ProcessResTerms
static int
xf_ProcessResTerms(xf_ResTerms *ResTerms, int iConv[], int iDiff[], int iSource[])
{
/*
PURPOSE:
 
  Arranges ResTerms so that they are in sequence: first Conv, then
  Diff, then Source.  Returns index ranges of all terms; i.e.

    ResTerms->ResTerm[iConv[0] .. iConv[1]-1] are convection terms
    ResTerms->ResTerm[iDiff[0] .. iDiff[1]-1] are diffusion terms
    ...
  
  The number of (e.g.) convection terms is iConv[1]-iConv[0].  If, for
  example, no convection terms are present, iConv[1] == iConv[0].
  Only Active residual terms are considered.

INPUTS:

  ResTerms: set of residual terms

OUTPUTS:

  iConv  : [start,end+1] range for convection terms (preallocated 2-int vector)
  iDiff  : [start,end+1] range for diffusion terms (preallocated 2-int vector)
  iSource: [start,end+1] range for source terms (preallocated 2-int vector)

RETURNS:

  Error code

*/
  int ierr, nResTerm, pos, i;
  xf_ResTerm *ResTermNew;

  nResTerm = ResTerms->nResTerm;

  if (nResTerm == 0){
    iConv[0] = iConv[1] = 0;
    iDiff[0] = iDiff[1] = 0;
    iSource[0] = iSource[1] = 0;
    return xf_OK;
  }

  ierr = xf_Error(xf_Alloc((void **) &ResTermNew, nResTerm, sizeof(xf_ResTerm)));
  if (ierr != xf_OK) return ierr;

  pos = 0;

  iConv[0] = pos;
  for (i=0; i<nResTerm; i++)
    if ((ResTerms->ResTerm[i].Type == xfe_ResTermConv) &&
	(ResTerms->ResTerm[i].Active)){
      ResTermNew[pos] = ResTerms->ResTerm[i];
      pos++;
    }
  iConv[1] = pos;

  iDiff[0] = pos;
  for (i=0; i<nResTerm; i++)
    if ((ResTerms->ResTerm[i].Type == xfe_ResTermDiff) &&
	(ResTerms->ResTerm[i].Active)){
      ResTermNew[pos] = ResTerms->ResTerm[i];
      pos++;
    }
  iDiff[1] = pos;

  iSource[0] = pos;
  for (i=0; i<nResTerm; i++)
    if ((ResTerms->ResTerm[i].Type == xfe_ResTermSource) &&
	(ResTerms->ResTerm[i].Active)){
      ResTermNew[pos] = ResTerms->ResTerm[i];
      pos++;
    }
  iSource[1] = pos;

  for (i=0; i<nResTerm; i++)
    if (!ResTerms->ResTerm[i].Active){
      ResTermNew[pos] = ResTerms->ResTerm[i];
      pos++;
    }
  
  /* Replace original pointer vec with new pointer vec */
  xf_Release( (void *) ResTerms->ResTerm);
  ResTerms->ResTerm = ResTermNew;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateSupportedVector
static int
xf_CreateSupportedVector(xf_All *All, const char Title[])
{
/*
PURPOSE:

  Creates a supported vector with name Title.  An example is the
  "WallDistance" vector.

INPUTS:
 
  All   : All structure
  Title : Name of vector to create

OUTPUTS: 

  None  : appropriate vector is stored in All->DataSet if created.

RETURN:

  Error Code
*/
  int ierr;

  if (strncmp(Title, "WallDistance", 12) == 0){
    // create a wall distance vector
    ierr = xf_Error(xf_CalculateDistFcn(All));
    if (ierr != xf_OK) return ierr;
  }
  else{
    xf_printf("Auxiliary vector = %s is not supported.\n", Title);
    return xf_Error(xf_NOT_SUPPORTED);
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_FindSupportedVector
int
xf_FindSupportedVector(xf_All *All, const char Title[], xf_Vector **pV)
{
  int ierr;
  xf_Data *D;

  ierr = xf_FindDataByTitle(All->DataSet, Title, xfe_Vector, &D);
  if (ierr == xf_NOT_FOUND){ 
    // try creating vector if supported (e.g. wall distance)
    ierr = xf_Error(xf_CreateSupportedVector(All, Title));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, Title, xfe_Vector, &D));
  }
  if (ierr != xf_OK) return xf_Error(ierr);
  
  (*pV) = (xf_Vector *) D->Data;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_RetrieveFcnParams
int
xf_RetrieveFcnParams(xf_All *All, const xf_EqnSet *EqnSet, 
		     int **pIParam, real **pRParam, int *pnAuxU, 
		     xf_Vector ***pAuxU)
{  
  int ierr, i, iAux;
  char *Name;
  xf_KeyValue *pKeyValue;
  xf_Data *D;

  pKeyValue = ((All == NULL) ? NULL : &All->Param->KeyValue);

  // Re-allocate memory for IParam and RParam
  ierr = xf_Error(xf_Alloc((void **) pIParam, EqnSet->nIParam, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc((void **) pRParam, EqnSet->nRParam, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<EqnSet->nIParam; i++){
    ierr = xf_GetKeyValueInt(EqnSet->KeyValue, EqnSet->IParamKey[i], (*pIParam)+i);
    if (ierr == xf_STRING_ERROR) // try a boolean read
      ierr = xf_GetKeyValueBool(EqnSet->KeyValue, EqnSet->IParamKey[i], 
				(enum xfe_Bool *) (*pIParam)+i);
    if ((ierr == xf_NOT_FOUND) && (pKeyValue != NULL)){
      ierr = xf_GetKeyValueInt((*pKeyValue), EqnSet->IParamKey[i], (*pIParam)+i);
      if (ierr == xf_STRING_ERROR) // try a boolean read
	ierr = xf_Error(xf_GetKeyValueBool((*pKeyValue), EqnSet->IParamKey[i], 
					   (enum xfe_Bool *) (*pIParam)+i));
      if (ierr != xf_OK){
	xf_printf("Error, could not find desired key = %s.\n", EqnSet->IParamKey[i]);
	return xf_Error(ierr);
      }
    }
    else
      if (ierr != xf_OK) return ierr;
  }

  for (i=0; i<EqnSet->nRParam; i++){
    ierr = xf_GetKeyValueReal(EqnSet->KeyValue, EqnSet->RParamKey[i], (*pRParam)+i);
    if ((ierr == xf_NOT_FOUND) && (pKeyValue != NULL)){
      ierr = xf_Error(xf_GetKeyValueReal((*pKeyValue), EqnSet->RParamKey[i], (*pRParam)+i));
      if (ierr != xf_OK){
	xf_printf("Error, could not find desired key = %s.\n", EqnSet->RParamKey[i]);
	return xf_Error(ierr);
      }
    }
    else
      if (ierr != xf_OK) return ierr;
  }

  // Auxiliary vectors to be passed into eqnset functions
  if (pAuxU != NULL) (*pAuxU) = NULL;
  if ((pnAuxU != NULL) &&  (((*pnAuxU) = EqnSet->nAuxU) != 0)){

    if ((pAuxU == NULL) || (All == NULL)) return xf_Error(xf_INPUT_ERROR);

    ierr = xf_Error(xf_Alloc((void **) pAuxU, (*pnAuxU), sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    
    for (iAux=0; iAux<(*pnAuxU); iAux++){
      (*pAuxU)[iAux] = NULL;
      Name = EqnSet->AuxUNames[iAux];

      ierr = xf_FindDataByTitle(All->DataSet, Name, xfe_Vector, &D);
      if (ierr == xf_NOT_FOUND){ 
	// try creating vector if supported (e.g. wall distance)
	ierr = xf_Error(xf_CreateSupportedVector(All, Name));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_FindDataByTitle(All->DataSet, Name, xfe_Vector, &D));
      }
      if (ierr != xf_OK) return xf_Error(ierr);
      
      (*pAuxU)[iAux] = (xf_Vector *) D->Data;

    } // iAux
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateStaticDataElem
static int
xf_CreateStaticDataElem(xf_All *All, xf_SolverData *SolverData, xf_StaticDataElem **pSD)
{
/*
PURPOSE: Initializes element static data

INPUTS: 
 
  All : all structure to determine residual terms
  SolverData : solver data structure

OUTPUTS: 

  (*pSD) : static data that is created and initialized

RETURNS: Error code

*/
  int ierr, iAux;
  xf_StaticDataElem *SD;
    
  // allocate memory
  ierr = xf_Error(xf_Alloc( (void **) pSD, 1, sizeof(xf_StaticDataElem)));
  if (ierr != xf_OK) return ierr;
  SD = (*pSD);
  
  // initialize variables to NULL
  SD->QuadData    =  NULL;
  SD->PhiData     =  NULL;
  SD->JData       =  NULL;
  SD->wq	  =  NULL;
  SD->GeomPhiData =  NULL;
  SD->xglob       =  NULL;
  //SD->T	          =  NULL;
  SD->u	          =  NULL; 
  SD->gu	  =  NULL; 
  SD->F	          =  NULL; 
  //SD->F_u	  =  NULL; 
  SD->Aw	  =  NULL; 
  //SD->A	          =  NULL; 
  //SD->A_uw        =  NULL; 
  //SD->StabVisc    =  NULL; 
  //SD->ResMetric   =  NULL; 
  //SD->AwStab      =  NULL; 
  SD->S	          =  NULL; 
  //SD->S_u	  =  NULL; 
  //SD->S_gu	  =  NULL; 
  SD->pnq         =  -1; 
  SD->nqConv      =  -1; 
  SD->nqDiff      =  -1; 
  //SD->Tsize       =  -1;
  SD->EG          =  NULL;
  //SD->EM          =  NULL;
  //SD->gum         =  NULL;

  // determine Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &SD->Time));
  if (ierr != xf_OK) return ierr;

  // build eqnset-desired parameter lists for passing into functions
  //ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &SD->IParam, &SD->RParam, 
  //			       &SD->nAuxU, &SD->AuxU));
  //if (ierr != xf_OK) return ierr;
  
  // Allocate vector of basis data for interpolating auxiliary vectors
  /*
  ierr = xf_Error(xf_Alloc( (void **) &SD->AuxPhiData, SD->nAuxU, sizeof(xf_BasisData *)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &SD->Auxu, SD->nAuxU, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  for (iAux=0; iAux<SD->nAuxU; iAux++){
    SD->AuxPhiData[iAux] = NULL;
    SD->Auxu[iAux]       = NULL;
  }
   
  if ( (SolverData != NULL) && (SolverData->StabRequired) ){
    // obtain elem geometry vector
    ierr = xf_Error(xf_FindElemGeom(All, &SD->EG));
    if (ierr != xf_OK) return ierr;
    // obtain elem metric vector
    ierr = xf_Error(xf_FindElemHMetric(All, xfe_True, &SD->EM));
    if (ierr != xf_OK) return ierr;
  }
    
  // initialize motion data
  SD->MD = NULL;
  if ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active)){
    ierr = xf_Error(xf_CreateMotionData(All, &SD->MD));
    if (ierr != xf_OK) return ierr;
  }

  // initialize data structure for residual linearization
  ierr = xf_Error(xf_InitLinQueue(SD->LinQ));
  if (ierr != xf_OK) return ierr;
  */
   
  // everything must go back with (*pSD)
  (*pSD) = SD;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyStaticDataElem
static int
xf_DestroyStaticDataElem(xf_StaticDataElem *SD)
{
/*
PURPOSE: Destroys static element structure

INPUTS: 
 
  SD : static data

OUTPUTS: None, SD and all of its contents are destroyed

RETURNS: Error code

*/

  int ierr, iAux;
  
  if (SD == NULL) return xf_OK;

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(SD->QuadData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(SD->PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(SD->GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(SD->JData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Auxiliary Vector Data */
  /*
  for (iAux=0; iAux<SD->nAuxU; iAux++){
    ierr = xf_Error(xf_DestroyBasisData(SD->AuxPhiData[iAux], xfe_True));
    if (ierr != xf_OK) return ierr;
    xf_Release( (void *) SD->Auxu[iAux]);
  }
  xf_Release( (void *) SD->AuxPhiData);
  xf_Release( (void *) SD->Auxu);
  xf_Release( (void *) SD->AuxU);
  */
   
  /* Destroy mesh motion data */
  //xf_DestroyMotionData(SD->MD);

  // Destroy residual linearization structure
  //xf_DestroyLinQueue(SD->LinQ);

  // Release all other memory
  //xf_Release( (void *) SD->IParam);
  //xf_Release( (void *) SD->RParam);
  xf_Release( (void *) SD->wq);
  xf_Release( (void *) SD->xglob);
  //xf_Release( (void *) SD->T);
  xf_Release( (void *) SD->u);
  xf_Release( (void *) SD->gu);
  xf_Release( (void *) SD->F);
  //xf_Release( (void *) SD->F_u);
  xf_Release( (void *) SD->Aw);
  //xf_Release( (void *) SD->A);
  //xf_Release( (void *) SD->A_uw);
  //xf_Release( (void *) SD->StabVisc);
  //xf_Release( (void *) SD->ResMetric);
  //xf_Release( (void *) SD->AwStab);
  xf_Release( (void *) SD->S);
  //xf_Release( (void *) SD->S_u);
  //xf_Release( (void *) SD->S_gu);
  //xf_Release( (void *) SD->gum);

  // Destroy self
  xf_Release( (void *) SD);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ModMotionConvF
static void
xf_ModMotionConvF(real *u, int nq, int dim, int sr, xf_MotionData *MD,
		  real *F, real *F_u)
{
/*
PURPOSE: 

  Modifies convective flux F and linearization F_u (if not NULL) to
  take into account motion of the mesh prescribed in the structure MD.

INPUTS: 
 
  u   : state vector at nq points [nq*sr]
  nq  : number of points
  dim : dimension
  sr  : state rank
  MD  : mesh motion data structure
  F   : convective flux [dim*nq*sr]
  F_u : convective flux linearization [dim*nq*sr*sr]

OUTPUTS:

  F, F_u : modified versions

RETURNS: Error code
*/
  int d, iq, k, i, dim2, sr2;

  dim2 = dim*dim;
  sr2  = sr*sr;

  // F{d,q,k} -= u{q,k}*vg{q,d}
  for (d=0; d<dim; d++) 
    xf_ColMult_Sub(u, MD->vg+d, nq, sr, dim, F+nq*sr*d);

  // multiply F{:,q,k} by the matrix g*Ginv{q,:,:}
  for (iq=0; iq<nq; iq++) 
    xf_dMxM(MD->Ginv+iq*dim2, MD->g[iq], dim, sr, nq*sr, F+iq*sr);

  /* modify F_u to incorporate the above changes to F */
  if (F_u != NULL){
    // F_u -= I*vg
    for (d=0; d<dim; d++)
      for (iq=0; iq<nq; iq++)
	for (k=0,i=(d*nq+iq)*sr2; k<sr*sr; k+=(sr+1))
	  F_u[i+k] -= MD->vg[iq*dim+d];
    // multiply F_u by g*Ginv/gb -- note, the 1/gb converts derivatives w.r.t u to w.r.t uX
    for (iq=0; iq<nq; iq++) 
      xf_dMxM(MD->Ginv+iq*dim2, MD->g[iq]/MD->gb[iq], dim, sr2, nq*sr2, F_u+iq*sr2);
  }

}



/******************************************************************/
//   FUNCTION Definition: xf_ModMotionPreEqnCall
void
xf_ModMotionPreEqnCall(int nq, int dim, int sr, xf_MotionData *MD,
		       real *u, real *wn)
{

  // divide state by gb: i.e. uX -> u = uX/gb (physical state)
  if (u  != NULL) xf_ColDiv(u, MD->gb, nq, sr, 1);
  // transform normal: multiply wn by g*G^{-T}
  if (wn != NULL) xf_ndMTxVc(nq, MD->Ginv, MD->g, dim, wn);

}


/******************************************************************/
//   FUNCTION Definition: xf_ModMotionPhysGrad
void
xf_ModMotionPhysGrad(int nq, int dim, int sr, xf_MotionData *MD,
		     real *u, real *gu)
{
  // convert reference gu to physical one
  int d, iq;

  // gu -= u*gb*gbigb_X
  for (d=0; d<dim; d++)
    xf_2ColcMult_Add(u, MD->gbigb_X+d, MD->gb, nq, sr, dim, 1, -1.0, gu+nq*sr*d);
  //xf_ColcMult_Add(u, MD->gig_X+d, nq, sr, dim, -1.0, gu+nq*sr*d);

  // gu{:,q,k} = (1/gb) * G^{-T} * gu
  for (iq=0; iq<nq; iq++) 
    xf_dMxMT(MD->Ginv+iq*dim*dim, 1./MD->gb[iq], dim, sr, nq*sr, gu+iq*sr);
  //xf_dMxM(MD->Ginv+iq*dim*dim, 1./MD->g[iq], dim, sr, nq*sr, gu+iq*sr);

}




/******************************************************************/
//   FUNCTION Definition: xf_ModMotionPostEqnCall
void
xf_ModMotionPostEqnCall(int nq, int dim, int sr, xf_MotionData *MD,
			real *u, real *wn)
{
  // multiply state by gb: i.e. u -> uX = gb*u (back to reference)
  if (u != NULL) xf_ColMult(u, MD->gb, nq, sr, 1);
  // transform normal back to reference: multiply wn by (1/g)*G^T
  if (wn != NULL) xf_ndMTxVic(nq, MD->G, MD->g, dim, wn);

}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualElem
static int
xf_CalculateResidualElem(xf_All *All, int egrp, int elem, xf_Vector *U, real *ER, 
			 real *ER_U, real **ER_NU, xf_StaticDataElem **pSD, 
			 xf_SolverData *SolverData)
{
/*
PURPOSE: 

  Calculates the residual and residual Jacobian (if requested)
  associated with an element interior integration on (egrp, elem).

INPUTS: 
 
  All: All structure
  egrp, elem : element in question
  U : state vector on all elements
  pSD : pointer to static data (optional, can pass in as NULL)
        useful for avoiding reallocations when calling multiple times
  SolverData : solver data structure

OUTPUTS:

  ER : element residual (must be preallocated)
  ER_U : element residual Jacobian (optional, but if given must be preallocated)
  ER_NU : element jacobian w.r.t neighbor states (optional, but preallocated if given)

RETURNS: Error code

*/
  int ierr, sr, sr2, iq, nq, nn, d, n, i, j;
  int dim, Order, QuadOrder, iAux, nAuxU;
  int ResOrder, ResidualOrderIncrement;
  int Tsize, pnq, nqConv, nqDiff;
  int  nConv,  nDiff,  nSource;
  int *iConv, *iDiff, *iSource;
  int *IParam;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged, ConstA, SkipDiffusion, SkipDiffStab, StabRequired;
  enum xfe_Bool Need_gu = xfe_False, Nonzero_gu = xfe_False;
  enum xfe_Bool Need_Grad;
  enum xfe_Bool MotionOn;
  real *EU, *F, *F_u, *Aw, *A, *A_uw, *AwStab, Gamma;
  real *S, *S_u, *S_gu, **Auxu = NULL;
  real *RParam, *xq, *wq, *xglob, *xq_glob, *u, *gu, *T, *gum = NULL, *w = NULL;
  real Time;
  real *StabVisc = NULL;
  real *AVelem;
  xf_ResTerm *ResTerm;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData, *GeomPhiData, **AuxPhiData;
  xf_JacobianData *JData;
  xf_StabData *StabData = NULL;
  xf_Vector *EG, **AuxU, *V, *GammaVec;
  xf_MotionData *MD = NULL;
  xf_LinQueueData *LinQ = NULL;
  xf_StaticDataElem *SD = NULL;
  xf_Mesh   *Mesh;
  xf_Data   *GammaDat;
  //xf_EqnSet *EqnSet;

  // General information
  Mesh    = All->Mesh;
  dim     = Mesh->Dim;  
  //EqnSet  = All->EqnSet;
  //sr      = EqnSet->StateRank;
  sr      = Model->nVars;
  sr2     = sr*sr;


  // Determine Basis and Order from the state, U
  Basis = U->Basis[egrp];
  Order = xf_InterpOrder(U, egrp, elem);

  // do we have source terms? quick check to see if we can skip p=0
  /*
  nSource = 0;
  for (i=0; i<EqnSet->ResTerms->nResTerm; i++)
    if ((EqnSet->ResTerms->ResTerm[i].Type == xfe_ResTermSource) &&
	(EqnSet->ResTerms->ResTerm[i].Active)) 
      nSource++;
  // skip Order = 0 elements if no source
  if ((Order == 0) && (nSource == 0)) return xf_OK;
  */
  //for classic finite volume (no need to compute this term
    if (Order == 0) return xf_OK;
  
  // Create and initialize static data if not passed in
  if ((pSD == NULL) || ((*pSD) == NULL)){
    ierr = xf_Error(xf_CreateStaticDataElem(All, SolverData, (pSD != NULL) ? pSD : &SD));
    if (ierr != xf_OK) return ierr;
  }
  if (pSD != NULL) SD = (*pSD);

  /*
  ResTerm = EqnSet->ResTerms->ResTerm;

  // Are we doing mesh motion?
  MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));

  // Pull off variables from StaticData  
  iConv   = SD->iConv;
  iDiff   = SD->iDiff;
  iSource = SD->iSource;
  nConv   =   iConv[1]-  iConv[0];
  nDiff   =   iDiff[1]-  iDiff[0]; 
  nSource = iSource[1]-iSource[0];
  */

  //find Gamma Vector
  ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
  if(ierr == xf_NOT_FOUND)
  {
     xf_printf("Cannot find heat capacity ratio...\n");
     return ierr;
  }
  else
     GammaVec = (xf_Vector *) GammaDat->Data;
  //if found, directly get the gamma value at this element
  Gamma = GammaVec->GenArray[egrp].rValue[elem][0];

  QuadData    = SD->QuadData;
  PhiData     = SD->PhiData;
  JData       = SD->JData;
  wq          = SD->wq;
  GeomPhiData = SD->GeomPhiData;
  xglob       = SD->xglob;
  //T           = SD->T;
  u           = SD->u;
  gu          = SD->gu;
  F           = SD->F;
  //F_u         = SD->F_u;
  Aw          = SD->Aw;
  //A           = SD->A;
  //A_uw        = SD->A_uw;
  //AwStab      = SD->AwStab;
  S           = SD->S;
  //S_u         = SD->S_u;
  //S_gu        = SD->S_gu;
  pnq         = SD->pnq;        // previous # quad points
  nqConv      = SD->nqConv;     // for allocation purposes
  nqDiff      = SD->nqDiff;
  //Tsize       = SD->Tsize;      // for temporary matrix allocation
  //IParam      = SD->IParam;
  //RParam      = SD->RParam;
  //nAuxU       = SD->nAuxU;      // number of auxiliary vectors
  //AuxU        = SD->AuxU;       // auxiliary vectors
  //AuxPhiData  = SD->AuxPhiData; // auxiliary vector basis data
  //Auxu        = SD->Auxu;       // auxiliary vector data
  EG          = SD->EG;         // element geometry vector
  Time        = SD->Time;       // simulation time
  //LinQ        = SD->LinQ;       // residual linearization structure
  //if (MotionOn){
  //  MD          = SD->MD;       // mesh motion data
  //  gum         = SD->gum;      // physical gradient for mesh motion
  //}
  
  // will we need the gradient of u?
  //if ((nDiff > 0) || (nSource > 0)) Need_gu = xfe_True;

  // do we need a linearization
  //Need_Grad = ((ER_U != NULL) || (ER_NU != NULL));

  // do we need stabilization
  //StabRequired = ( (SolverData != NULL) && (SolverData->StabRequired) );

  // Residual order increase
  //ResidualOrderIncrement = ((SolverData == NULL) ? 0 : SolverData->ResidualOrderIncrement);
  
  // Residual order
  //ResOrder = (ResidualOrderIncrement == 0) ? Order : max(Order + ResidualOrderIncrement, 0);

  // determine required integration order 
  
  ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order, &QuadOrder));
  if (ierr != xf_OK) return ierr;

  if(!Model->TwistFlag)
     QuadOrder -= Mesh->Dim;

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
  
  nn = PhiData->nn; // number of interpolation nodes
  
  // re-allocate data if quad points increased
  if (nq > pnq){
    ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
   
    //for diffusion implementation (May 2013)
    if (Model->DiffFlag){
      ierr = xf_Error(xf_ReAlloc( (void **) &gu, nq*sr*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    /*
    if (StabRequired){
      ierr = xf_Error(xf_ReAlloc( (void **) &SD->StabVisc, nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc( (void **) &SD->ResMetric, nq*dim*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    */
  }
  
  // obtain global coords of quad points
  ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
				  nq, xq, xglob));
  if (ierr != xf_OK) return ierr;

  // obtain transformation map if doing mesh motion
  /*if (MotionOn){
    ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, All->Mesh->Motion, 
				      nq, dim, Time, xglob, MD));
    if (ierr != xf_OK) return ierr;
  }*/
  
  EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
  
  // interpolate state and gradient at quad points
  xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u);      
  
  if (Model->DiffFlag)
    for (d=0; d<dim; d++)
      xf_MxM_Set(PhiData->gPhi+nn*nq*d, EU, nq, nn, sr, gu+nq*sr*d);
  /*
  // interpolate any auxiliary vectors, using AuxPhiData
  for (iAux=0; iAux<nAuxU; iAux++){
    V  = AuxU[iAux];
    if (nq > pnq){ // reallocate Auxu[iAux] if necessary
      ierr = xf_Error(xf_ReAlloc( (void **) Auxu+iAux, nq*V->StateRank, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_EvalBasis(V->Basis[egrp], xf_InterpOrder(V,egrp,elem), QuadChanged, 
				 nq, xq, xfb_Phi, AuxPhiData+iAux));
    if (ierr != xf_OK) return ierr;
    xf_MxM_Set(AuxPhiData[iAux]->Phi, V->GenArray[egrp].rValue[elem], nq, 
	       AuxPhiData[iAux]->nn, V->StateRank, Auxu[iAux]);
  }*/
  

  // form detJ-multiplied quad weight vector, wq
  for (iq=0; iq<nq; iq++) 
    wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
  
  // pull off data required for stabilization
  /*
  if ( StabRequired ){
    StabData = &(SolverData->StabData);
    ierr = xf_Error(xf_CalculateStabViscElem(All, egrp, elem, nq, xq, SD->EM->GenArray[egrp].rValue[elem], 
					     StabData, SD->StabVisc, SD->ResMetric, &SkipDiffStab));
    if (ierr != xf_OK) return ierr;
    StabVisc = SD->StabVisc;
  }
  else {
    StabVisc = NULL;
    SkipDiffStab = xfe_False;
  }*/

  /* Clear residual linearization queue */
  //if (ER_U != NULL) xf_ClearLinQueue(LinQ);

  /*------------------*/
  /* CONVECTION TERMS */
  /*------------------*/

  //if ((nConv > 0) & (Order != 0)){

    if (nq > nqConv){ // realloc F, F_u if necessary	
      nqConv = nq;
      ierr = xf_Error(xf_ReAlloc( (void **) &F, sr*nq*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
  //    if (ER_U != NULL){
  //	ierr = xf_Error(xf_ReAlloc( (void **) &F_u, sr*sr*nq*dim, sizeof(real)));
  //	if (ierr != xf_OK) return ierr;
  //    }
    }	

    // transform state to physical
    //if (MotionOn) xf_ModMotionPreEqnCall(nq, dim, sr, MD, u, NULL);

    // calculate F [and F_u] at quad points
    //ierr = xf_Error(xf_EqnSetConvF(EqnSet, ResTerm+iConv[0], nConv, IParam, 
    //				   RParam, nq, u, Auxu, xglob, F, F_u));
    ierr = xf_Error(ConvFluxInterior(nq, sr, dim, u, F, Gamma));
    if (ierr != xf_OK) return ierr;

    // apply mesh motion modifications
    //if (MotionOn) xf_ModMotionConvF(u, nq, dim, sr, MD, F, F_u);

    // transform state to ref
    //if (MotionOn) xf_ModMotionPostEqnCall(nq, dim, sr, MD, u, NULL);
	
    // multiply F by quad weights*J
    for (d=0; d<dim; d++)
      xf_ColMult(F+nq*sr*d, wq, nq, sr, 1); // F is modified here
	
    /* Add to R:
       ER{n,k} -= sum_i sum_q gPhi{i,q,n}^T * F{i,q,k}*wq{q}
    */
    for (d=0; d<dim; d++)
      xf_MTxM_Sub(PhiData->gPhi+nn*nq*d, F+nq*sr*d, nn, nq, sr, ER);
	
    /* Add to ER_U:
       ER_U{n,k;m,a} -= sum_i sum_q gPhi{i,q,n}^T * F_u{i,q,k;a}*Phi{q,m}
    */
    //if (ER_U != NULL){
    //  ierr = xf_Error(xf_AddToLinQueue(F_u, xfe_LinQTerm_GPhiPhi, nq, dim, 
	//			       sr2, -1, wq, 1, -1.0, LinQ));
    //  if (ierr != xf_OK) return ierr;
    //
    //}

  //} // end if nConv > 0

      
  /*-----------------*/
  /* DIFFUSION TERMS */
  /*-----------------*/
      
  /* Skip diffusion terms if only have a stabilization term and
     the stabilization viscosity is 0. */
  if(Model->DiffFlag){
//  SkipDiffusion = ((nDiff == 1) && SkipDiffStab);

//  if ((nDiff > 0) && (!SkipDiffusion) && (Order != 0)){  

    if (nq > nqDiff){ // realloc Aw, A, A_uw if necessary	
      nqDiff = nq;
      ierr = xf_Error(xf_ReAlloc( (void **) &Aw, sr*nq*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
/*      if (Need_Grad){
	ierr = xf_Error(xf_ReAlloc( (void **) &A, sr2*nq*dim*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &A_uw, sr2*nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }
      if (StabVisc != NULL){
	ierr = xf_Error(xf_ReAlloc( (void **) &AwStab, sr*nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }
      if (MotionOn){
	ierr = xf_Error(xf_ReAlloc( (void **) &gum, sr*nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }
*/    }	
    // realloc temporary matrix for R_U construction, if necessary
/*    if ((Need_Grad) && (max(nn,sr)*nq > Tsize)){
      Tsize = max(nn,sr)*nq;
      ierr = xf_Error(xf_ReAlloc( (void **) &T, Tsize, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }

    if (MotionOn){ // Modify state gradient in presence of mesh motion
      // gum = gu
      for (i=0; i<nq*sr*dim; i++) gum[i] = gu[i];
      // gum -= u*gbigb_X  -- u here is in ref space
      for (d=0; d<dim; d++)
        xf_ColcMult_Add(u, MD->gbigb_X+d, nq, sr, dim, -1.0, gum+nq*sr*d);
      w = gum;
    }
    else 
*/      w = gu;

    // calculate Aw [and A, A_uw] at quad points
//    ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm+iDiff[0], nDiff, IParam, RParam, SD->ResMetric, 
//				    StabVisc, nq, u, gu, w, xfe_True, MD, Aw, A, A_uw, AwStab, 
//				    &ConstA, NULL));
//    if (ierr != xf_OK) return ierr;
    if(Model->AVmodel){  
       AVelem = Model->AVmodel_data->GenArray[egrp].rValue[elem];
       //hook-in to avoid the activation around boundary
       //otherwise, we might have terrible instabilty there
       /* 
       for (iq=0; iq<nq; iq++){
           xq_glob = xglob + dim * iq;
           if(fabs(xq_glob[0])>0.0)
           {
              for(d=1; d<sr+1; d++)
                 AVelem[d] = 0.0;
          
              break;
           }
        }*/
    }
    else 
       AVelem = NULL;

    ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, u, gu, w, Aw, NULL, Gamma, AVelem));
    if (ierr != xf_OK) return ierr;

    // multiply Aw by quad weights
    for (d=0; d<dim; d++)
      xf_ColMult(Aw+nq*sr*d, wq, nq, sr, 1);
    
    // Add to R:
    //   ER{n,k} += sum_i sum_q gPhi{i,q,n}^T * Aw{i,q,k}*wq{q}
    
    for (d=0; d<dim; d++)
      xf_MTxM_Add(PhiData->gPhi+nn*nq*d, Aw+nq*sr*d, nn, nq, sr, ER);

    // Add to ER_U:
    //   ER_U{n,k;m,a} += sum_i sum_q wq{q} * gPhi{i,q,n}^T *
    //                    [ A{i,j,q,k,a}*gPhi{j,q,m} + A_uw{i,q,k,a}*Phi{q,m} 
	//		 -A{i,j,q,k,a}*gbigb_X{q,j}*Phi{q,m}  (if MotionOn) ]
    //   A_uw{i,q,k,a} = A_u{i,j,q,k,l; a} * w{j,q,l}  (sum over j,l)
/*  we do explicit solving; no need for linearization 
    
    if (ER_U != NULL){

      // A contribution
      ierr = xf_Error(xf_AddToLinQueue(A, xfe_LinQTerm_GPhiGPhi, nq, dim, 
				       sr2, -1, wq, 1, 1.0, LinQ));
      if (ierr != xf_OK) return ierr;

      // A_uw contribution
      if (!ConstA){
	ierr = xf_Error(xf_AddToLinQueue(A_uw, xfe_LinQTerm_GPhiPhi, nq, dim, 
					 sr2, -1, wq, 1, 1.0, LinQ));
	if (ierr != xf_OK) return ierr;
      }

      // mesh motion contribution (A is overwritten here)
      if (MotionOn){
        // A{i,0,q,k,a} = sum_j A{i,j,q,k,a}*gbigb_X{q,j}
        for (i=0; i<dim; i++){
          xf_ColMult(A+nq*sr2*i*dim, MD->gbigb_X+0, nq, sr2, dim);
          for (j=1; j<dim; j++)
            xf_ColMult_Add(A+nq*sr2*(i*dim+j), MD->gbigb_X+j, nq, sr2, dim, A+nq*sr2*i*dim);
	  
          ierr = xf_Error(xf_AddToLinQueue(A+nq*sr2*i*dim, xfe_LinQTerm_GPhiPhi, nq, dim,
        				   sr2, i, wq, 1, -1.0, LinQ));
          if (ierr != xf_OK) return ierr;
        } // i
      }
    }
*/    
    // Linearization of stabilization
//    if ((StabRequired) && ((ER_U != NULL) || (ER_NU != NULL)) ){
//
      // ER_U{n,k;m,a} += gPhi{i,n,q}*AwStab{i,q,k}*wq{q}*StabVisc_U{q,m,a}
//      ierr = xf_Error(xf_AddLinearizationStabViscElem(All, egrp, elem, sr, nn, nq, xq, StabVisc, 
//						      StabData, PhiData->gPhi, wq, AwStab,
//						      T, ER_U, ER_NU));
//      if (ierr != xf_OK) return ierr;
//    }
	
  } // end if nDiff > 0
  
   
  /*--------------*/
  /* SOURCE TERMS */
  /*--------------*/
  
 // if (nSource > 0){
  if((Model->ChemSource && !Model->DetailChem) || (Model->MMSSource)) {

    if (nq > pnq){	// realloc S, S_u if necessary
      ierr = xf_Error(xf_ReAlloc( (void **) &S, sr*nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
  /*    if (ER_U != NULL){
	ierr = xf_Error(xf_ReAlloc( (void **) &S_u, sr*sr*nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	if (Need_gu){
	  ierr = xf_Error(xf_ReAlloc( (void **) &S_gu, sr*sr*nq*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}
      }
  */
    }
	
    // realloc temporary matrix for ER_U construction, if necessary
   /* if ((ER_U != NULL) && (nn*nq > Tsize)){
      Tsize = nn*nq;
      ierr = xf_Error(xf_ReAlloc( (void **) &T, Tsize, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }

    // transform state to physical
    if (MotionOn) xf_ModMotionPreEqnCall(nq, dim, sr, MD, u, NULL);
    */
    // calculate S [and S_u] at quad points
    //ierr = xf_Error(xf_EqnSetSourceS(EqnSet, ResTerm+iSource[0], nSource, IParam, RParam, 
    //				     nq, u, gu, Auxu, xglob, &Time, S, S_u, S_gu, &Nonzero_gu));
    ierr = xf_Error(ChemicalSource(nq, sr, dim, u, xglob, S, Gamma, Model->DetailChem, Model->dt_size));
    if (ierr != xf_OK) return ierr;

    // transform state to ref
    //if (MotionOn) xf_ModMotionPostEqnCall(nq, dim, sr, MD, u, NULL);


    // multiply quad weights by g if MeshMotion is on (last time quad weights used)
    //if (MotionOn) for (iq=0; iq<nq; iq++) wq[iq] *= MD->g[iq];
	
    // multiply S by quad weights*J
    xf_ColMult(S, wq, nq, sr, 1); // S is modified here
	
    // Add to R:
    //   ER{n,k} += sum_q Phi{q,n}^T * S{q,k}*wq{q}
    
    xf_MTxM_Add(PhiData->Phi, S, nn, nq, sr, ER);
	
    // Add to ER_U:
    //   ER_U{n,k;m,a} += sum_q Phi{q,n}^T * S_u{q,k;a}*Phi{q,m}
   /* 
    if (ER_U != NULL){
      
      ierr = xf_Error(xf_AddToLinQueue(S_u, xfe_LinQTerm_PhiPhi, nq, dim, 
				       sr2, -1, wq, 1, 1.0, LinQ));
      if (ierr != xf_OK) return ierr;

      // Nonzero_gu indicates that S depends on the gradient of u, gu
	 //ER_U{n,k;m,a} += sum_i sum_q Phi{q,n}^T * S_gu{i,q,k;a}*gPhi{i,q,m}
	 //(note, this is a dual-inconsistent discretization at this point)
      
      if (Nonzero_gu){
	ierr = xf_Error(xf_AddToLinQueue(S_gu, xfe_LinQTerm_PhiGPhi, nq, dim,
					 sr2, -1, wq, 1, 1.0, LinQ));
	if (ierr != xf_OK) return ierr;
      }
    }
    */
  } // end if nSource > 0

  // apply linearizations queued up
  //if (ER_U != NULL){
  //  ierr = xf_Error(xf_ApplyLinQueue(LinQ, PhiData, PhiData, sr2, ER_U));
  //  if (ierr != xf_OK) return ierr;
  //}

  pnq = nq; // set previous quad point # for next element

  // Store possibly-altered or resized data back in StaticData
  SD->QuadData    =  QuadData;
  SD->PhiData     =  PhiData;  
  SD->JData       =  JData;
  SD->wq          =  wq;
  SD->GeomPhiData =  GeomPhiData;
  SD->xglob       =  xglob;
  //SD->T	          =  T;
  SD->u	          =  u;
  SD->gu	  =  gu;
  SD->F	          =  F;
  //SD->F_u	  =  F_u;
  SD->Aw	  =  Aw;
  //SD->A	          =  A;
  //SD->A_uw        =  A_uw;
  //SD->AwStab      =  AwStab;
  SD->S	          =  S;
  //SD->S_u	  =  S_u;
  //SD->S_gu        =  S_gu;
  SD->pnq         =  pnq;    
  SD->nqConv      =  nqConv; 
  SD->nqDiff      =  nqDiff;
  //SD->Tsize       =  Tsize;  
  //SD->AuxPhiData  =  AuxPhiData; 
  //SD->Auxu        =  Auxu;
  //SD->LinQ[0]     =  *LinQ; // not really necessary
  //if (MotionOn){
  //  SD->MD        =  MD;
  //  SD->gum       =  gum;
  //}
  
  if (pSD == NULL){
    // Delete StaticData that we just created
    ierr = xf_Error(xf_DestroyStaticDataElem(SD));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;

}




/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualElems
static int
xf_CalculateResidualElems(xf_All *All, xf_Vector *U, xf_Vector *R, 
			  xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
/*
PURPOSE: 

  Calculates the residual and residual Jacobian (if requested)
  associated with an element interior integration on all elements.

INPUTS: 
 
  All: All structure
  U : state vector on all elements
  SolverData : solver data structure

OUTPUTS:

  R : residual
  R_U : element Jacobian (optional, or may be given but Value may not exist)

RETURNS: Error code

*/


  int ierr, egrp, elem;
  real *ER, *ER_U;
  real **ER_NU = NULL;
  xf_StaticDataElem *StaticData = NULL;
  xf_Mesh   *Mesh;

  Mesh = All->Mesh;

  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      /* Prepare pointers to residual and Jacobian */
      ER = R->GenArray[egrp].rValue[elem];
      if ((R_U == NULL) || (R_U->Value == NULL)){
	ER_U  = NULL;
	ER_NU = NULL;
      }
      else{
	ER_U  = R_U->Value[egrp][elem][0];
	ER_NU = R_U->Value[egrp][elem]+1;
      }

      /* Calculate residual on elem, passing in StaticData  */
      ierr = xf_Error(xf_CalculateResidualElem(All, egrp, elem, U, ER, ER_U, ER_NU, 
					       &StaticData, SolverData));
      if (ierr != xf_OK) return ierr;
      
    } // elem

  } // egrp


  // Delete StaticData
  ierr = xf_Error(xf_DestroyStaticDataElem(StaticData));
  if (ierr != xf_OK) return ierr;


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateStaticDataIFace
static int
xf_CreateStaticDataIFace(xf_All *All, xf_SolverData *SolverData, xf_StaticDataIFace **pSD)
{
/*
PURPOSE: 

  Creates and initializes iface static data

INPUTS: 
 
  All : all structure to determine residual terms
  SolverData : solver data structure

OUTPUTS:

  (*pSD) : static data that is created and initialized

RETURNS: Error code

*/
  int ierr, k, iAux;
  xf_StaticDataIFace *SD;

  // allocate memory
  ierr = xf_Error(xf_Alloc( (void **) pSD, 1, sizeof(xf_StaticDataIFace)));
  if (ierr != xf_OK) return ierr;
  SD = (*pSD);
  
  //ierr = xf_Error(xf_ProcessResTerms(All->EqnSet->ResTerms, SD->iConv, 
  //				     SD->iDiff, SD->iSource));
  //if (ierr != xf_OK) return ierr;
  
  // initialize variables to NULL
  SD->QuadData    =  NULL;
  SD->PhiDataL    =  NULL;   SD->PhiDataR    =  NULL;
  SD->NData       =  NULL;
  SD->JDataL      =  NULL;   SD->JDataR      =  NULL;
  SD->GeomPhiData =  NULL;
  SD->wn	  =  NULL;
  SD->xglob       =  NULL;
  SD->xelemL      =  NULL;   SD->xelemR      =  NULL;
  SD->uL          =  NULL;   SD->uR          =  NULL; 
  SD->guL	  =  NULL;   SD->guR         =  NULL; 
  SD->F	          =  NULL; 
  //SD->F_uL        =  NULL;   SD->F_uR        =  NULL; 
  //SD->CF          =  NULL;
  SD->DData	  =  NULL;
  SD->pnq         =  -1; 
  SD->nqConv      =  -1; 
  SD->nqDiff      =  -1; 
  SD->nnDiff      =  -1; 
  SD->EG          =  NULL;
  //SD->EM          =  NULL;
  //SD->StabViscL   =  NULL;   SD->StabViscR   =  NULL;
  //SD->StabPhiL    =  NULL;   SD->StabPhiR    =  NULL;
  //SD->ResMetricL  =  NULL;   SD->ResMetricR  =  NULL;
  //SD->ResPhiDataL =  NULL;   SD->ResPhiDataR =  NULL;
  //SD->ResPhiTable =  NULL;

  // determine Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &SD->Time));
  if (ierr != xf_OK) return ierr;

  // build eqnset-desired parameter lists for passing into functions
  //ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &SD->IParam, &SD->RParam, 
  //				       &SD->nAuxU, &SD->AuxU));
  //if (ierr != xf_OK) return ierr;

  // create structure for storing jump data 
  ierr = xf_Error(LYDG_CreateDiffJumpData(&SD->DData));
  if (ierr != xf_OK) return ierr;

  // Allocate vector of basis data for interpolating auxiliary vectors
    /*
  ierr = xf_Error(xf_Alloc( (void **) &SD->AuxPhiDataL, SD->nAuxU, sizeof(xf_BasisData *)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &SD->AuxPhiDataR, SD->nAuxU, sizeof(xf_BasisData *)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &SD->AuxuL, SD->nAuxU, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &SD->AuxuR, SD->nAuxU, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  for (iAux=0; iAux<SD->nAuxU; iAux++){
    SD->AuxPhiDataL[iAux] = NULL;
    SD->AuxuL[iAux]       = NULL;
    SD->AuxPhiDataR[iAux] = NULL;
    SD->AuxuR[iAux]       = NULL;
  }*/

  /* Create a basis table, PhiTable, that will store computed basis
     functions specific to each [element shape, face in element,
     orientation of face] combination, for quick lookup. */
  ierr = xf_Error(xf_CreateBasisTable(&SD->PhiTable));
  if (ierr != xf_OK) return ierr;

  /* Create another basis table, GeomPhiTable, for storing the
     geometry approximation basis functions evaluated at
     faces/orientations/shapes. */
  ierr = xf_Error(xf_CreateBasisTable(&SD->GeomPhiTable));
  if (ierr != xf_OK) return ierr;

  // find element geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &SD->EG));
  if (ierr != xf_OK) return ierr;

  /*
  if ( (SolverData != NULL) && (SolverData->StabRequired) ){
    // obtain elem metric vector
    ierr = xf_Error(xf_FindElemHMetric(All, xfe_True, &SD->EM));
    if (ierr != xf_OK) return ierr;
  }

  // initialize motion data
  SD->MDL = SD->MDR = NULL;
  if ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active)){
    ierr = xf_Error(xf_CreateMotionData(All, &SD->MDL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_CreateMotionData(All, &SD->MDR));
    if (ierr != xf_OK) return ierr;
  }

  // initialize data structure for residual linearization
  for (k=0; k<4; k++){
    ierr = xf_Error(xf_InitLinQueue(SD->LinQ+k));
    if (ierr != xf_OK) return ierr;
  }*/

  /* Create another basis table, ResPhiTable, that will store computed
     basis functions specific to each [element shape, face in element,
     orientation of face] at an incremented order. */
  //if ( (SolverData != NULL) && (SolverData->ResidualOrderIncrement != 0)) {
  //  ierr = xf_Error(xf_CreateBasisTable(&SD->ResPhiTable));
  //  if (ierr != xf_OK) return ierr;
  //}

  // everything must go back with (*pSD)
  (*pSD) = SD;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyStaticDataIFace
static int
xf_DestroyStaticDataIFace(xf_StaticDataIFace *SD)
{
/*
PURPOSE: Destroys static iface structure

INPUTS: 
 
  SD : static data

OUTPUTS: None, SD and all of its contents are destroyed

RETURNS: Error code

*/

  int ierr, k, iAux;

  if (SD == NULL) return xf_OK;

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(SD->QuadData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(SD->PhiDataL, xfe_False));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyBasisData(SD->PhiDataR, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Table */
  ierr = xf_Error(xf_DestroyBasisTable(SD->PhiTable));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(SD->GeomPhiData, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Destroy Geometry Basis Table */
  ierr = xf_Error(xf_DestroyBasisTable(SD->GeomPhiTable));
  if (ierr != xf_OK) return ierr;

  /* Destroy Normal Data */
  ierr = xf_Error(xf_DestroyNormalData(SD->NData));
  if (ierr != xf_OK) return ierr;

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(SD->JDataL));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyJacobianData(SD->JDataR));
  if (ierr != xf_OK) return ierr;

  /* Destroy Diff Data */
  xf_DestroyDiffJumpData(SD->DData);

  /* Destroy Auxiliary Vector Data */
    /*
  for (iAux=0; iAux<SD->nAuxU; iAux++){
    ierr = xf_Error(xf_DestroyBasisData(SD->AuxPhiDataL[iAux], xfe_True));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyBasisData(SD->AuxPhiDataR[iAux], xfe_True));
    if (ierr != xf_OK) return ierr;
    xf_Release( (void *) SD->AuxuL[iAux]);
    xf_Release( (void *) SD->AuxuR[iAux]);
  }
  xf_Release( (void *) SD->AuxPhiDataL);
  xf_Release( (void *) SD->AuxPhiDataR);
  xf_Release( (void *) SD->AuxuL);
  xf_Release( (void *) SD->AuxuR);
  xf_Release( (void *) SD->AuxU);
  */
  /* Destroy mesh motion data */
  //xf_DestroyMotionData(SD->MDL);
  //xf_DestroyMotionData(SD->MDR);

  // Destroy residual linearization structures
  //for (k=0; k<4; k++) xf_DestroyLinQueue(SD->LinQ+k);

  /* Destroy Residual Basis Data */
  //ierr = xf_Error(xf_DestroyBasisData(SD->ResPhiDataL, xfe_False));
  //if (ierr != xf_OK) return ierr;

  //ierr = xf_Error(xf_DestroyBasisData(SD->ResPhiDataR, xfe_False));
  //if (ierr != xf_OK) return ierr;

  /* Destroy Residual Basis Table */
  //ierr = xf_Error(xf_DestroyBasisTable(SD->ResPhiTable));
  //if (ierr != xf_OK) return ierr;


  // Release memory
  //xf_Release( (void *) SD->IParam);
  //xf_Release( (void *) SD->RParam);
  xf_Release( (void *) SD->wn);
  xf_Release( (void *) SD->xglob);
  xf_Release( (void *) SD->xelemL);
  xf_Release( (void *) SD->xelemR);
  xf_Release( (void *) SD->uL);
  xf_Release( (void *) SD->uR);
  xf_Release( (void *) SD->guL);
  xf_Release( (void *) SD->guR);
  xf_Release( (void *) SD->F);
  //xf_Release( (void *) SD->F_uL);
  //xf_Release( (void *) SD->F_uR);
  //xf_Release( (void *) SD->CF);
  //xf_Release( (void *) SD->StabViscL);
  //xf_Release( (void *) SD->StabPhiL);
  //xf_Release( (void *) SD->ResMetricL);

  // Destroy self
  xf_Release( (void *) SD);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ModMotionConvFJump
static void
xf_ModMotionConvFJump(int nq, int dim, int sr, xf_MotionData *MD, real *F_u)
{
/*
PURPOSE: 

  Modifies convective jump flux F (result of Riemann solver) and
  linearizations F_uL and F_uR (if not NULL) to take into account
  motion of the mesh prescribed in the structure MD.  Currently the
  transformation is incorporated into the normal vector so that a
  change to F is not required.

INPUTS: 
 
  nq   : number of points
  dim  : dimension
  sr   : state rank
  MD   : mesh motion data structure
  (not currently passed in) F    : convective jump flux [nq*sr]
  F_u  : convective jump flux linearization with respect to state [nq*sr*sr]

OUTPUTS:

  F_u : modified version

RETURNS: Error code
*/

  int iq, k, i, sr2;
  real gb;

  sr2  = sr*sr;

  /* multiply F_u 1/gb */
  if (F_u != NULL){
    for (iq=0; iq<nq; iq++)
      for (k=0, gb=MD->gb[iq], i=iq*sr2; k<sr2; k++) F_u[i+k] /= gb;
  }

}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualIFace
static int
xf_CalculateResidualIFace(xf_All *All, int iiface, xf_Vector *U, real *RL,
			  real *RR, real *RL_UL, real *RR_UR, real *RL_UR, 
			  real *RR_UL, xf_StaticDataIFace **pSD, xf_SolverData *SolverData)
{
/*
PURPOSE: 

  Calculates the residuals and residual Jacobians (if requested)
  associated with the left and right elements due to integration on
  interior face IFace

INPUTS: 
 
  All: All structure
  iiface: number of interior face in question
  U : state vector on all elements
  pSD : pointer to static data (optional, can pass in as NULL)
        useful for avoiding reallocations when calling multiple times
  SolverData : solver data structure

OUTPUTS:

  RL, RR : residuals on left and right elements (must be preallocated)
  RL_UL, RR_UR :  self Jacobians (optional, but if given must be preallocated)
  RL_UR, RR_UL : cross Jacobians (optional, but if given must be preallocated)


RETURNS: Error code

*/

  int ierr, sr, sr2, iq, nq, nn, nL, nR, d, n, m, i, j, k;
  int dim, OrderL, OrderR, QuadOrder;
  int ResOrderL, ResOrderR, MaxOrder;
  int egrpL, elemL, faceL;
  int egrpR, elemR, faceR;
  int pnq, nqConv, nqDiff, nnDiff;
  int ResidualOrderIncrement, NewOrder, Rnn;
  int  nConv,  nDiff,  nSource;
  int *iConv, *iDiff, *iSource;
  int iAux, nAuxU;
  int *IParam;
  enum xfe_BasisType BasisL, BasisR;
  enum xfe_DiffDiscType DiffDisc;
  enum xfe_Bool QuadChanged, AnyR_U;
  enum xfe_Bool Need_gu = xfe_False;
  enum xfe_Bool SkipDiffusion, SkipDiffStab;
  enum xfe_Bool StabRequired;
  enum xfe_Bool MotionOn;
  real *RParam, *xq, *wq, *wn, *xglob, *xelemL, *xelemR;
  real *UL, *UR, *uL, *uR, *guL, *guR, GammaL, GammaR;
  real *uLghost, *uRghost;
  real **AuxuL = NULL, **AuxuR = NULL;
  real *F, *F_uL, *F_uR, *CF, cval, ElemVolL, ElemVolR;
  real *phiL, *gphiL, t, LRU[25];
  real Time;
  xf_ResTerm *ResTerm;
  xf_QuadData *QuadData;
  xf_BasisTable *PhiTable, *GeomPhiTable;
  xf_BasisData *PhiDataL, *PhiDataR, *GeomPhiData;
  xf_BasisTable *ResPhiTable;
  xf_BasisData *ResPhiDataL, *ResPhiDataR;
  xf_BasisData **AuxPhiDataL, **AuxPhiDataR;
  xf_NormalData *NData;
  xf_JacobianData *JDataL, *JDataR;
  xf_DiffJumpData *DData;
  xf_Vector *EG, *C, *V, **AuxU, *GammaVec;
  xf_MotionData *MDL = NULL, *MDR = NULL;
  xf_LinQueueData *LinQ = NULL, *LinQLL, *LinQLR, *LinQRL, *LinQRR;
  xf_StaticDataIFace *SD = NULL;
  xf_StabData *StabData = NULL;
  xf_IFace IFace;
  xf_Mesh   *Mesh;
  xf_Data   *GammaDat;
  //xf_EqnSet *EqnSet;

  //added in Feb 2014
  real MaxCharSpeed;
  xf_Data *MaxC_Dat;
  xf_Vector *MaxC_Vec, *AdaptIndicator;

  real tmp1, tmp;

  // General information
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  sr  = Model->nVars;
  sr2 = sr*sr;

  // Interior face structure
  IFace = Mesh->IFace[iiface];
  

  // Create and initialize static data if not passed in
  if ((pSD == NULL) || ((*pSD) == NULL)){
    ierr = xf_Error(xf_CreateStaticDataIFace(All, SolverData, (pSD != NULL) ? pSD : &SD));
    if (ierr != xf_OK) return ierr;
  }
  if (pSD != NULL) SD = (*pSD);
 
  //ResTerm = EqnSet->ResTerms->ResTerm;

  // Are we doing mesh motion?
  //MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));
  
  // Pull off variables from StaticData  
  //iConv   = SD->iConv;
  //iDiff   = SD->iDiff;
  //iSource = SD->iSource;
  //nConv   =   iConv[1]-  iConv[0];
  //nDiff   =   iDiff[1]-  iDiff[0]; 
  //nSource = iSource[1]-iSource[0];
  
  //find Gamma Vector
  ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
  if(ierr == xf_NOT_FOUND)
  {
     xf_printf("Cannot find heat capacity ratio...\n");
     return ierr;
  }
  else
     GammaVec = (xf_Vector *) GammaDat->Data;

  //find data/vector for maximum char speed
  ierr = xf_FindDataByTitle(All->DataSet, "ElemMaxCharSpeed", xfe_Vector, &MaxC_Dat);
  if(ierr == xf_NOT_FOUND)
  {
     xf_printf("Cannot find maximum characteristic speed...\n");
     return ierr;
  }
  else
     MaxC_Vec = (xf_Vector *) MaxC_Dat->Data;

  QuadData    = SD->QuadData;
  PhiDataL    = SD->PhiDataL;       PhiDataR    = SD->PhiDataR;
  PhiTable    = SD->PhiTable;
  NData       = SD->NData;
  JDataL      = SD->JDataL;         JDataR      = SD->JDataR;
  GeomPhiData = SD->GeomPhiData;
  GeomPhiTable= SD->GeomPhiTable;
  wn          = SD->wn;
  xglob       = SD->xglob;
  xelemL      = SD->xelemL;         xelemR      = SD->xelemR;
  uL          = SD->uL;             uR          = SD->uR;
  guL         = SD->guL;            guR         = SD->guR;
  F           = SD->F;
  //F_uL        = SD->F_uL;           F_uR        = SD->F_uR;
  //CF          = SD->CF;
  DData       = SD->DData;
  pnq         = SD->pnq;        
  nqConv      = SD->nqConv;     
  nqDiff      = SD->nqDiff;
  nnDiff      = SD->nnDiff;
  //IParam      = SD->IParam;
  //RParam      = SD->RParam;
  //nAuxU       = SD->nAuxU;      
  //AuxPhiDataL = SD->AuxPhiDataL;    AuxPhiDataR = SD->AuxPhiDataR;
  //AuxuL       = SD->AuxuL;          AuxuR       = SD->AuxuR;
  //AuxU        = SD->AuxU;      
  EG          = SD->EG;      
  Time        = SD->Time;      // simulation time
  //LinQ        = SD->LinQ;      // residual linearization structure
  //if (MotionOn){
  //  MDL       = SD->MDL;       // mesh motion data  
  //  MDR       = SD->MDR;       // mesh motion data  
  //}
  //ResPhiDataL = SD->ResPhiDataL;    ResPhiDataR  = SD->ResPhiDataR;
  //ResPhiTable = SD->ResPhiTable;
  
  // if we have diffusion terms
  /*
    if (nDiff > 0){
    Need_gu = xfe_True;  // we will need the gradient, and we need the type of discretization
    ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "DiffusionDiscretization", 
				       xfe_DiffDiscName, (int ) xfe_DiffDiscLast, 
				       (int *) &DiffDisc));
    if (ierr != xf_OK) return ierr;
  }

  // Set Connectivity vector if required
  if ((SolverData != NULL) && (SolverData->CRequired))
    C = SolverData->C;
  else
    C = NULL;

  // Stabilization required ?
  StabRequired = ((SolverData != NULL) && (SolverData->StabRequired));

  // Residual order increase
  ResidualOrderIncrement = ((SolverData == NULL) ? 0 : SolverData->ResidualOrderIncrement);
*/
    
  // elements on L and R
  egrpL = IFace.ElemGroupL;
  egrpR = IFace.ElemGroupR;
  elemL = IFace.ElemL;
  elemR = IFace.ElemR;
  faceL = IFace.FaceL;
  faceR = IFace.FaceR;
  
  BasisL = U->Basis[egrpL];
  BasisR = U->Basis[egrpR];
  OrderL = xf_InterpOrder(U, egrpL, elemL);
  OrderR = xf_InterpOrder(U, egrpR, elemR);
  MaxOrder = max(OrderL, OrderR);

  //Gamma vector is found, then go to evaluate gamma value for left and right elements
  GammaL = GammaVec->GenArray[egrpL].rValue[elemL][0];
  GammaR = GammaVec->GenArray[egrpR].rValue[elemR][0];

  // residual orders (may be different from state)
  /*if (ResidualOrderIncrement != 0){
    // do not allow order to go below 0
    ResOrderL = max(OrderL + ResidualOrderIncrement, 0); 
    ResOrderR = max(OrderR + ResidualOrderIncrement, 0);
    MaxOrder = max(ResOrderL, ResOrderR);
  }*/

  // determine required integration order
  ierr = xf_Error(xf_GetQuadOrderIFace(Mesh, NULL, IFace, MaxOrder, &QuadOrder));
  if (ierr != xf_OK) return ierr;

  if(!Model->TwistFlag)
     QuadOrder -= (Mesh->Dim - 1);

  /* Pull off quad points for the iface; will not recalculate if
     Basis/Order have not changed. */
  ierr = xf_Error(xf_QuadFace(Mesh, egrpL, elemL, faceL, 
			      QuadOrder, &QuadData, &QuadChanged));
  if (ierr != xf_OK) return ierr;
    
  nq = QuadData->nquad;
  xq = QuadData->xquad;
  wq = QuadData->wquad;
  
  // compute basis functions if quad or basis or order changed
  ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrpL, elemL, faceL, IFace.OrientL,
					       BasisL, OrderL, QuadChanged, nq, xq, 
					       xfb_Phi | xfb_GPhi | xfb_gPhi, 
					       &PhiDataL, PhiTable, &xelemL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrpR, elemR, faceR, IFace.OrientR,
					       BasisR, OrderR, QuadChanged, nq, xq, 
					       xfb_Phi | xfb_GPhi | xfb_gPhi, 
					       &PhiDataR, PhiTable, &xelemR));
  if (ierr != xf_OK) return ierr;

  // Also compute lower/higher order basis functions for residual eval, if needed
  /*if (ResidualOrderIncrement != 0){
    if (ResPhiTable == NULL) return xf_Error(xf_CODE_LOGIC_ERROR); // sanity check
    ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrpL, elemL, faceL, IFace.OrientL,
                                                 BasisL, ResOrderL, QuadChanged, nq, xq, 
                                                 xfb_Phi, &ResPhiDataL, ResPhiTable, NULL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrpR, elemR, faceR, IFace.OrientR,
                                                 BasisR, ResOrderR, QuadChanged, nq, xq, 
                                                 xfb_Phi, &ResPhiDataR, ResPhiTable, NULL));
    if (ierr != xf_OK) return ierr;
  }*/ 


  /* C ompute normal(s) at quad points.  If face is straight, only
     one normal will be computed/returned. */
  ierr = xf_Error(xf_IFaceNormal(Mesh, IFace, nq, xq, &NData));
  if (ierr != xf_OK) return ierr;

  if (Model->DiffFlag){
    /* Compute geometry Jacobian; if not constant, compute at quad
       points.  Note if jacobian is constant, only one Jacobian will
       be computed/returned. */
    ierr = xf_Error(xf_ElemJacobian(Mesh, egrpL, elemL, nq, xelemL, 
				    xfb_detJ | xfb_iJ, xfe_True, &JDataL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ElemJacobian(Mesh, egrpR, elemR, nq, xelemR, 
				    xfb_detJ | xfb_iJ, xfe_True, &JDataR));
    if (ierr != xf_OK) return ierr;
      
    // convert reference basis grads (GPhi) to physical grads, gPhi
    ierr = xf_Error(xf_EvalPhysicalGrad(PhiDataL, JDataL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_EvalPhysicalGrad(PhiDataR, JDataR));
    if (ierr != xf_OK) return ierr;
  }
    
  nL = PhiDataL->nn;
  nR = PhiDataR->nn;
  nn = max(nL, nR);
  //if (ResidualOrderIncrement != 0){  // account for possibly-different residual order
  //  Rnn = max(ResPhiDataL->nn, ResPhiDataR->nn);
  //  nn = max(nn, Rnn);
  //}

  UL = U->GenArray[egrpL].rValue[elemL]; // U on elemL [nL*sr]
  UR = U->GenArray[egrpR].rValue[elemR]; // U on elemR [nR*sr]

  // re-allocate data if quad points increased
  if (nq > pnq){
    ierr = xf_Error(xf_ReAlloc( (void **) &wn, nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &uL, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &uR, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    if (Model->DiffFlag){
      ierr = xf_Error(xf_ReAlloc( (void **) &guL, nq*sr*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc( (void **) &guR, nq*sr*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
/*
    if (C != NULL){
      ierr = xf_Error(xf_ReAlloc( (void **) &CF, nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }

    if (StabRequired){
      ierr = xf_Error(xf_ReAlloc( (void **) &SD->StabViscL, 2*nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      SD->StabViscR = SD->StabViscL + nq;
      ierr = xf_Error(xf_ReAlloc( (void **) &SD->StabPhiL, 2*nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      SD->StabPhiR = SD->StabPhiL + nq;
      ierr = xf_Error(xf_ReAlloc( (void **) &SD->ResMetricL, 2*nq*dim*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      SD->ResMetricR = SD->ResMetricL + nq*dim*dim;
    }*/

  }

  // construct wn = weighted normals = normals multiplied by quad weights
  for (d=0; d<dim; d++)
    for (iq=0;iq<nq; iq++) 
      wn[iq*dim+d] = NData->n[iq*dim*(NData->nq!=1)+d]*wq[iq];

  // pull off element volumes
  ElemVolL = EG->GenArray[egrpL].rValue[elemL][xfe_EGVolume];
  ElemVolR = EG->GenArray[egrpR].rValue[elemR][xfe_EGVolume];

  // obtain global coords of quad points (use L, should be same as R if mesh is valid)
  ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrpL, elemL, &GeomPhiData,
				  xfe_True, nq, xelemL, xglob));
  if (ierr != xf_OK) return ierr;
/*   ierr = xf_Error(xf_Ref2GlobFaceUsingTable(Mesh, egrpL, elemL, faceL, IFace.OrientL, */
/* 					    QuadChanged, &GeomPhiData, GeomPhiTable, */
/* 					    nq, xelemL, xglob)); */
/*   if (ierr != xf_OK) return ierr; */



  // obtain transformation map if doing mesh motion
  /*if (MotionOn){
    ierr = xf_Error(xf_MeshMotionMap( egrpL, elemL, PhiDataL, All->Mesh->Motion, 
				      nq, dim, Time, xglob, MDL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_MeshMotionMap( egrpR, elemR, PhiDataR, All->Mesh->Motion, 
				      nq, dim, Time, xglob, MDR));
    if (ierr != xf_OK) return ierr;
  }*/

  // interpolate state and gradient at quad points
  xf_MxM_Set(PhiDataL->Phi, UL, nq, nL, sr, uL);
  xf_MxM_Set(PhiDataR->Phi, UR, nq, nR, sr, uR);
  if (Model->DiffFlag)
    for (d=0; d<dim; d++){
      xf_MxM_Set(PhiDataL->gPhi+nL*nq*d, UL, nq, nL, sr, guL+nq*sr*d);
      xf_MxM_Set(PhiDataR->gPhi+nR*nq*d, UR, nq, nR, sr, guR+nq*sr*d);
    }


  //use the velocity C1-continuity indicator
   if(Model->Dyn_p_Adapt)
   {
      ierr = xf_Error(Yu_CreateOrFindAdaptIndicator(All, xfe_False, NUM_VAR_DYM_P, &AdaptIndicator));
      if (ierr != xf_OK) return ierr;
     
     
      for(iq=0; iq<nq; iq++)
      {
         tmp = 0.;
         for (d=0; d<dim; d++)
            tmp += wn[iq*dim + d] * fabs(guL[d*nq*sr + iq*sr + 1] - guR[d*nq*sr + iq*sr + 1]); 
    
         tmp1 = 0.;
         for (d=0; d<dim; d++)
            tmp1 += wn[iq*dim + d] * fabs(guL[d*nq*sr + iq*sr + 2] - guR[d*nq*sr + iq*sr + 2]);

         AdaptIndicator->GenArray[egrpL].rValue[elemL][0] += sqrt(tmp*tmp + tmp1*tmp1);
         AdaptIndicator->GenArray[egrpR].rValue[elemR][0] += sqrt(tmp*tmp + tmp1*tmp1);
      }
   }
/*
  // zero out Connectivity
  if (C != NULL){
    C->GenArray[egrpL].rValue[elemL][faceL] = 0.;
    C->GenArray[egrpR].rValue[elemR][faceR] = 0.;
  }


  // pull off data required for stabilization
  if (StabRequired){
    StabData = &(SolverData->StabData);
    ierr = xf_Error(xf_CalculateStabViscIFace(All, iiface, nq, xelemL, xelemR,
					      SD->EM->GenArray[egrpL].rValue[elemL],
					      SD->EM->GenArray[egrpR].rValue[elemR],
					      StabData, SD->StabViscL, SD->StabViscR,
					      &StabData->StabViscL_UL, &StabData->StabViscL_UR,
					      &StabData->StabViscR_UL, &StabData->StabViscR_UR,
					      SD->StabPhiL, SD->StabPhiR,
					      SD->ResMetricL, SD->ResMetricR, &SkipDiffStab));
    if (ierr != xf_OK) return ierr;
    StabData->StabViscL    = SD->StabViscL;
    StabData->StabViscR    = SD->StabViscR;
    StabData->StabPhiL     = SD->StabPhiL;
    StabData->StabPhiR     = SD->StabPhiR;
    StabData->ResMetricL   = SD->ResMetricL;
    StabData->ResMetricR   = SD->ResMetricR;
    StabData->SkipDiffStab = SkipDiffStab;
  }
  else{
    StabData = NULL;
    SkipDiffStab = xfe_False;
  }

  // interpolate any auxiliary vectors, using AuxPhiData
  for (iAux=0; iAux<nAuxU; iAux++){
    V  = AuxU[iAux];
    if (nq > pnq){ // reallocate Auxu[iAux] if necessary
      ierr = xf_Error(xf_ReAlloc( (void **) AuxuL+iAux, nq*V->StateRank, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc( (void **) AuxuR+iAux, nq*V->StateRank, sizeof(real)));
      if (ierr != xf_OK) return ierr;

    }
    ierr = xf_Error(xf_EvalBasis(V->Basis[egrpL], xf_InterpOrder(V,egrpL,elemL), xfe_True,
				 nq, xelemL, xfb_Phi, AuxPhiDataL+iAux));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_EvalBasis(V->Basis[egrpR], xf_InterpOrder(V,egrpR,elemR), xfe_True,
				 nq, xelemR, xfb_Phi, AuxPhiDataR+iAux));
    if (ierr != xf_OK) return ierr;

    xf_MxM_Set(AuxPhiDataR[iAux]->Phi, V->GenArray[egrpR].rValue[elemR], nq, 
	       AuxPhiDataR[iAux]->nn, V->StateRank, AuxuR[iAux]);
    xf_MxM_Set(AuxPhiDataL[iAux]->Phi, V->GenArray[egrpL].rValue[elemL], nq, 
	       AuxPhiDataL[iAux]->nn, V->StateRank, AuxuL[iAux]);
	    
  }
  
  // will we need any Jacobians?
  AnyR_U = ((RL_UL != NULL) || (RR_UR != NULL) || (RL_UR != NULL) || (RR_UL != NULL));
*/
  /* Clear residual linearization queues */
  //LinQLL=LinQ+0;  LinQLR=LinQ+1;  LinQRL=LinQ+2;  LinQRR=LinQ+3;
  //if (RL_UL != NULL) xf_ClearLinQueue(LinQLL);
  //if (RL_UR != NULL) xf_ClearLinQueue(LinQLR);
  //if (RR_UL != NULL) xf_ClearLinQueue(LinQRL);
  //if (RR_UR != NULL) xf_ClearLinQueue(LinQRR);
      
  /*------------------*/
  /* CONVECTION TERMS */
  /*------------------*/
  //if (nConv > 0){  
    
    // realloc F if necessary
    if (nq > nqConv){ 
      nqConv = nq;
      ierr = xf_Error(xf_ReAlloc( (void **) &F, sr*nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;

    //  if (AnyR_U){
	//ierr = xf_Error(xf_ReAlloc( (void **) &F_uL, sr*sr*nq, sizeof(real)));
	//if (ierr != xf_OK) return ierr;
	//ierr = xf_Error(xf_ReAlloc( (void **) &F_uR, sr*sr*nq, sizeof(real)));
	//if (ierr != xf_OK) return ierr;
    //  }
    }
        
    // transform state and normal to physical space
    //if (MotionOn){
    // xf_ModMotionPreEqnCall(nq, dim, sr, MDL, uL, wn  );
    //  xf_ModMotionPreEqnCall(nq, dim, sr, MDR, uR, NULL);
    //}

    // calculate F [and F_u] at quad points
    //ierr = xf_Error(xf_EqnSetConvFJump(EqnSet, ResTerm+iConv[0], nConv, IParam, RParam, 
    //				       nq, uL, uR, AuxuL, AuxuR, wn, xglob, 
    //				       (MotionOn) ? MDL->vg : NULL, F, F_uL, F_uR, CF));
    if(fabs(GammaL - GammaR) < Model->Gammathreshold)
    {
    ierr = xf_Error(ConvFluxInteriorFace(nq, sr, dim, uL, uR, wn, F, GammaR, &MaxCharSpeed));
    if (ierr != xf_OK) return ierr;

    //update max speed for elemL & elemR
    if(MaxC_Vec->GenArray[egrpL].rValue[elemL][0] < MaxCharSpeed)
       MaxC_Vec->GenArray[egrpL].rValue[elemL][0] = MaxCharSpeed;

    if(MaxC_Vec->GenArray[egrpR].rValue[elemR][0] < MaxCharSpeed)
       MaxC_Vec->GenArray[egrpR].rValue[elemR][0] = MaxCharSpeed;
    // Modify fluxes if motion is on 
    //if (MotionOn){
    //  xf_ModMotionConvFJump(nq, dim, sr, MDL, F_uL);
    //  xf_ModMotionConvFJump(nq, dim, sr, MDR, F_uR);
    //}

    // transform state and normal back to reference space
    //if (MotionOn){
    //  xf_ModMotionPostEqnCall(nq, dim, sr, MDL, uL,   wn);
    //  xf_ModMotionPostEqnCall(nq, dim, sr, MDR, uR, NULL);
    //}

    /* Add to R:
       RL_{n,k} += sum_q PhiL_{n,q}*F_{k,q}
       RR_{n,k} -= sum_q PhiR_{n,q}*F_{k,q}
    */
    if (RL != NULL) xf_MTxM_Add(PhiDataL->Phi, F, nL, nq, sr, RL);
    if (RR != NULL) xf_MTxM_Sub(PhiDataR->Phi, F, nR, nq, sr, RR);
    }
    else
    {
    
       ierr = xf_Error(xf_Alloc( (void **) &uLghost, nq*sr, sizeof(real)));
       if (ierr != xf_OK) return ierr;
       ierr = xf_Error(xf_Alloc( (void **) &uRghost, nq*sr, sizeof(real)));
       if (ierr != xf_OK) return ierr;

       //construct ghost solution according to double flux model
       ierr = xf_Error(Yu_GhostStateConstruct(nq, sr, dim, uL, uR, uLghost, uRghost, GammaL, GammaR));
       if( ierr != xf_OK) return ierr;
    
       ierr = xf_Error(ConvFluxInteriorFace(nq, sr, dim, uL, uRghost, wn, F, GammaL, &MaxCharSpeed));
       if (ierr != xf_OK) return ierr;
       if (RL != NULL) xf_MTxM_Add(PhiDataL->Phi, F, nL, nq, sr, RL);
       
       ierr = xf_Error(ConvFluxInteriorFace(nq, sr, dim, uLghost, uR, wn, F, GammaR, &MaxCharSpeed));
       if (ierr != xf_OK) return ierr;
       if (RR != NULL) xf_MTxM_Sub(PhiDataR->Phi, F, nR, nq, sr, RR);

       xf_Release( (void *) uLghost);
       xf_Release( (void *) uRghost);
    }
    /* Add to R_U:
       RL_UL{n,k;m,a} += sum_q PhiL_{n,q}*F_uL{k,q;a}*PhiL_{m,q}
       RL_UR{n,k;m,a} += sum_q PhiL_{n,q}*F_uR{k,q;a}*PhiR_{m,q}
       RR_UL{n,k;m,a} -= sum_q PhiR_{n,q}*F_uL{k,q;a}*PhiL_{m,q}
       RR_UR{n,k;m,a} -= sum_q PhiR_{n,q}*F_uR{k,q;a}*PhiR_{m,q}
    */
/*
    if (RL_UL != NULL){
      ierr = xf_Error(xf_AddToLinQueue(F_uL, xfe_LinQTerm_PhiPhi, nq, dim, 
				       sr2, -1, NULL, 0, 1.0, LinQLL));
      if (ierr != xf_OK) return ierr;
    }
    if (RL_UR != NULL){
      ierr = xf_Error(xf_AddToLinQueue(F_uR, xfe_LinQTerm_PhiPhi, nq, dim, 
				       sr2, -1, NULL, 0, 1.0, LinQLR));
      if (ierr != xf_OK) return ierr;
    }
    if (RR_UL != NULL){
      ierr = xf_Error(xf_AddToLinQueue(F_uL, xfe_LinQTerm_PhiPhi, nq, dim, 
				       sr2, -1, NULL, 0, -1.0, LinQRL));
      if (ierr != xf_OK) return ierr;
    }
    if (RR_UR != NULL){
      ierr = xf_Error(xf_AddToLinQueue(F_uR, xfe_LinQTerm_PhiPhi, nq, dim, 
				       sr2, -1, NULL, 0, -1.0, LinQRR));
      if (ierr != xf_OK) return ierr;
    }
    
    if (C != NULL){ // add to connectivity
      for (iq=0, cval=0.; iq<nq; iq++) cval += CF[iq];
      C->GenArray[egrpL].rValue[elemL][faceL] += fabs(cval);
      C->GenArray[egrpR].rValue[elemR][faceR] += fabs(cval);
    }
*/
  //} // end if nConv > 0

  /*-----------------*/
  /* DIFFUSION TERMS */
  /*-----------------*/
  
  // Skip diffusion terms if only have a stabilization term and
  //   the regularity flags are 0. 
  //SkipDiffusion = ((nDiff == 1) && SkipDiffStab);

  if(Model->DiffFlag && !Model->AVmodel){
  //if(Model->DiffFlag){// && !Model->AVmodel){

//  if ((nDiff > 0) && (!SkipDiffusion)){

    if ((nq > nqDiff) || (nn > nnDiff)){  // realloc DiffData if necesasry
      nqDiff = max(nq, nqDiff); 
      nnDiff = max(nn, nnDiff);
      ierr = xf_Error(LYDG_ReAllocDiffJumpData(sr, nqDiff, dim, nnDiff,
					       xfe_False, DData));
      if (ierr != xf_OK) return ierr;
    }

    //return xf_NOT_SUPPORTED;

     
    //   Compute viscous flux integrated on IFace: 
	// 
    //   RL{n,k} -= int( PhiL{n} * Qn{k} )
    //   RR{n,k} += int( PhiR{n} * Qn{k} )
	// 
    //   And derivatives w.r.t UL and UR:
	// 
    //   RL_UL{n,k} -= int( PhiL{n} * Qn_UL{k} )
    //   RL_UR{n,k} -= int( PhiL{n} * Qn_UR{k} )
    //   RR_UL{n,k} += int( PhiR{n} * Qn_UL{k} )
    //   RR_UR{n,k} += int( PhiR{n} * Qn_UR{k} )
	// 
    //   The viscous flux dotted with the normal (nL) is:
	// 
    //   Qn = 0.5*(ALguL+ARguR)*nL + [stabilization terms]
	//  
    //   where the stabilization terms depend on the type of
    //   viscous disrectization (BR2, etc.)
	// 
    //   Note: dunL, dunR, AL, AR, etc. are filled in in DData
/*    
    ierr = xf_Error(xf_DiffFluxQJump(All, iiface, ResTerm+iDiff[0], nDiff, DiffDisc,
				     IParam, RParam, PhiDataL, PhiDataR, 
                                     (ResidualOrderIncrement == 0) ? PhiDataL : ResPhiDataL,
                                     (ResidualOrderIncrement == 0) ? PhiDataR : ResPhiDataR,
                                     StabData, ElemVolL, ElemVolR, BasisL, BasisR, OrderL, OrderR, 
				     nq, wn, uL, uR, guL, guR, MDL, MDR, RL, RR, RL_UL, 
				     RL_UR, RR_UL, RR_UR, LinQ, DData, CF));
    if (ierr != xf_OK) return ierr;
*/   
    ierr = xf_Error(LYDG_FaceViscTerm(All, Model, iiface, PhiDataL, PhiDataR, PhiDataL,
                                      PhiDataR, ElemVolL, ElemVolR, BasisL, BasisR,
                                      OrderL, OrderR, nq, wn, uL, uR, guL, guR, RL, RR,
                                      GammaL, GammaR, DData));
    if (ierr != xf_OK) return ierr;

    // Dual consistency term
    // Add A*(u-uhat)*n term to RL and RR:
    //
    //   RL{n,k} -= gPhiL{i,q,n} * ALdunL{i,q,k},  q,i summed
    //   RR{n,k} -= gPhiR{i,q,n} * ARdunR{i,q,k},  q,i summed
    //
    //   Note: normal and quad weights are included in dunL and dunR
   
    if (RL != NULL)
      for (d=0; d<dim; d++)
	xf_MTxM_Sub(PhiDataL->gPhi+nL*nq*d, DData->ALdunL+sr*nq*d, nL, nq, sr, RL);
    if (RR != NULL)
      for (d=0; d<dim; d++)
	xf_MTxM_Sub(PhiDataR->gPhi+nR*nq*d, DData->ARdunR+sr*nq*d, nR, nq, sr, RR);

    // Add derivatives of A*(u-uhat)*n term to R_U:
    //
    //   RL_UL{n,k;m,a} -= gPhiL{i,q,n} * [ A_uLdunL{i,q,k,a} * PhiL{q,m}
    //   + AL{i,j,q,k,a}*PhiL{q,m}*wn{q,j}*duL_uL ]
    //   RL_UR{n,k;m,a} -= gPhiL{i,q,n} * [
    //   + AL{i,j,q,k,a}*PhiR{q,m}*wn{q,j}*duL_uR ]
    //   RR_UL{n,k;m,a} -= gPhiR{i,q,n} * [
    //   + AR{i,j,q,k,a}*PhiL{q,m}*wn{q,j}*duR_uL ]
    //   RR_UR{n,k;m,a} -= gPhiR{i,q,n} * [ A_uRdunR{i,q,k,a} * PhiR{q,m}
    //   + AR{i,j,q,k,a}*PhiR{q,m}*wn{q,j}*duR_uR ]
    //
    //   Note: nL is stored with quad weights in wn;
    //   duL_uL, etc. are scalars; q and i are summed
    
/*
    if ((RL_UL != NULL) || (RL_UR != NULL)){
      for (i=0; i<dim; i++){

	// Set AL{i,0,q,k,a} = sum_j AL{i,j,q,k,a}*wn{q,j}  (AL is overwritten here)
	xf_ColMult(DData->AL+nq*sr2*i*dim, wn+0, nq, sr2, dim); 
	for (j=1; j<dim; j++)
	  xf_ColMult_Add(DData->AL+nq*sr2*(i*dim+j), wn+j, nq, sr2, dim, DData->AL+nq*sr2*i*dim); 
	
	// Contribution of AL to RL_UL linearization
	if (RL_UL != NULL){
	  ierr = xf_Error(xf_AddToLinQueue(DData->AL+nq*sr2*i*dim, xfe_LinQTerm_GPhiPhi, 
					   nq, dim, sr2, i, NULL, 0, -DData->duL_uL, LinQLL));
	  if (ierr != xf_OK) return ierr;
	}

	// Contribution of AL to RL_UR linearization
	if (RL_UR != NULL){
	  ierr = xf_Error(xf_AddToLinQueue(DData->AL+nq*sr2*i*dim, xfe_LinQTerm_GPhiPhi, 
					   nq, dim, sr2, i, NULL, 0, -DData->duL_uR, LinQLR));
	  if (ierr != xf_OK) return ierr;
	}
      } // i
    }
*/
    // Contribution of A_uLdunL to RL_UL
 /*   
    if ((RL_UL != NULL) && (!DData->ConstAL)){
      ierr = xf_Error(xf_AddToLinQueue(DData->A_uLdunL, xfe_LinQTerm_GPhiPhi, 
				       nq, dim, sr2, -1, NULL, 0, -1.0, LinQLL));
      if (ierr != xf_OK) return ierr;
    }

    if ((RR_UL != NULL) || (RR_UR != NULL)){
      for (i=0; i<dim; i++){

	// Set AR{i,0,q,k,a} = sum_j AR{i,j,q,k,a}*wn{q,j}  (AR is overwritten here)
	xf_ColMult(DData->AR+nq*sr2*i*dim, wn+0, nq, sr2, dim); 
	for (j=1; j<dim; j++)
	  xf_ColMult_Add(DData->AR+nq*sr2*(i*dim+j), wn+j, nq, sr2, dim, DData->AR+nq*sr2*i*dim); 
	
	// Contribution of AR to RR_UL linearization
	if (RR_UL != NULL){
	  ierr = xf_Error(xf_AddToLinQueue(DData->AR+nq*sr2*i*dim, xfe_LinQTerm_GPhiPhi, 
					   nq, dim, sr2, i, NULL, 0, -DData->duR_uL, LinQRL));
	  if (ierr != xf_OK) return ierr;
	}

	// Contribution of AR to RR_UR linearization
	if (RR_UR != NULL){
	  ierr = xf_Error(xf_AddToLinQueue(DData->AR+nq*sr2*i*dim, xfe_LinQTerm_GPhiPhi, 
					   nq, dim, sr2, i, NULL, 0, -DData->duR_uR, LinQRR));
	  if (ierr != xf_OK) return ierr;
	}
      } // i
    }

    // Contribution of A_uRdunR to RR_UR
    if ((RR_UR != NULL) && (!DData->ConstAR)){
      ierr = xf_Error(xf_AddToLinQueue(DData->A_uRdunR, xfe_LinQTerm_GPhiPhi, 
				       nq, dim, sr2, -1, NULL, 0, -1.0, LinQRR));
      if (ierr != xf_OK) return ierr;
    }

    
    // Include linearization of stabilization term
    if ( ((RL_UL != NULL) || (RL_UR != NULL)) && (StabData != NULL) ){
      // assumption: StabVisc = 0 implies StabVisc_U = 0 (valid for continuous slope switch)
      // RL_UL{n,k;m,a} -= gPhiL{i,n,q}*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UL{m,a}
      // RL_UR{n,k;m,a} -= gPhiL{i,n,q}*AwStabL{i,q,k}*StabPhiL{q}*StabViscL_UR{m,a}
      xf_MTxwM_Set(PhiDataL->gPhi+nL*nq*0, StabData->StabPhiL, DData->AwStabL+0*nq*sr, nL, nq, sr, DData->T);
      for (i=1; i<dim; i++)
	xf_MTxwM_Add(PhiDataL->gPhi+nL*nq*i, StabData->StabPhiL, DData->AwStabL+i*nq*sr, nL, nq, sr, DData->T);
      if ((RL_UL != NULL) && (StabData->StabViscL_UL != NULL))
	xf_BlockOutProd_Sub(DData->T, StabData->StabViscL_UL, nL, sr, nL, RL_UL);
      if ((RL_UR != NULL) && (StabData->StabViscL_UR != NULL))
	xf_BlockOutProd_Sub(DData->T, StabData->StabViscL_UR, nL, sr, nL, RL_UR);
    }
    
    if ( ((RR_UL != NULL) || (RR_UR != NULL)) && (StabData != NULL) ){
      // assumption: StabVisc = 0 implies StabVisc_U = 0 (valid for continuous slope switch)
      // RR_UL{n,k;m,a} -= gPhiR{i,n,q}*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UL{m,a}
      // RR_UR{n,k;m,a} -= gPhiR{i,n,q}*AwStabR{i,q,k}*StabPhiR{q}*StabViscR_UR{m,a}
      xf_MTxwM_Set(PhiDataR->gPhi+nR*nq*0, StabData->StabPhiR, DData->AwStabR+0*nq*sr, nR, nq, sr, DData->T);
      for (i=1; i<dim; i++)
	xf_MTxwM_Add(PhiDataR->gPhi+nR*nq*i, StabData->StabPhiR, DData->AwStabR+i*nq*sr, nR, nq, sr, DData->T);
      if ((RR_UL != NULL) && (StabData->StabViscR_UL != NULL))
	xf_BlockOutProd_Sub(DData->T, StabData->StabViscR_UL, nR, sr, nR, RR_UL);
      if ((RR_UR != NULL) && (StabData->StabViscR_UR != NULL))
	xf_BlockOutProd_Sub(DData->T, StabData->StabViscR_UR, nR, sr, nR, RR_UR);
    }

    if (C != NULL){ // add to connectivity
      for (iq=0, cval=0.; iq<nq; iq++) cval += CF[iq];
      C->GenArray[egrpL].rValue[elemL][faceL] += fabs(cval);
      C->GenArray[egrpR].rValue[elemR][faceR] += fabs(cval);
    }
*/
  } // end if nDiff > 0
 
  // apply linearizations queued up
/*  
  if (RL_UL != NULL){
    ierr = xf_Error(xf_ApplyLinQueue(LinQLL, PhiDataL, PhiDataL, sr2, RL_UL));
    if (ierr != xf_OK) return ierr;
  }
  if (RL_UR != NULL){
    ierr = xf_Error(xf_ApplyLinQueue(LinQLR, PhiDataL, PhiDataR, sr2, RL_UR));
    if (ierr != xf_OK) return ierr;
  }
  if (RR_UL != NULL){
    ierr = xf_Error(xf_ApplyLinQueue(LinQRL, PhiDataR, PhiDataL, sr2, RR_UL));
    if (ierr != xf_OK) return ierr;
  }
  if (RR_UR != NULL){
    ierr = xf_Error(xf_ApplyLinQueue(LinQRR, PhiDataR, PhiDataR, sr2, RR_UR));
    if (ierr != xf_OK) return ierr;
  }
*/

  pnq = nq; // set previous quad point # for next iface
  
  // Store possibly-altered or resized data back in StaticData
  SD->QuadData    = QuadData;
  SD->PhiDataL    = PhiDataL;       SD->PhiDataR    = PhiDataR;
  SD->PhiTable    = PhiTable;	          
  SD->NData       = NData;	          
  SD->JDataL      = JDataL;         SD->JDataR      = JDataR;
  SD->GeomPhiData = GeomPhiData;        
  SD->GeomPhiTable= GeomPhiTable; 
  SD->wn          = wn;	          
  SD->xglob       = xglob;	          
  SD->xelemL      = xelemL;         SD->xelemR      = xelemR;
  SD->uL          = uL;             SD->uR          = uR;
  SD->guL         = guL;            SD->guR         = guR;
  SD->F           = F;	          
  //SD->F_uL        = F_uL;           SD->F_uR        = F_uR;
  //SD->CF          = CF;	          
  SD->DData       = DData;	          
  SD->pnq         = pnq;                
  SD->nqConv      = nqConv;             
  SD->nqDiff      = nqDiff;	          
  SD->nnDiff      = nnDiff;	          
  //SD->IParam      = IParam;	          
  //SD->RParam      = RParam;	          
  //SD->nAuxU       = nAuxU;              
  //SD->AuxPhiDataL = AuxPhiDataL;    SD->AuxPhiDataR = AuxPhiDataR;
  //SD->AuxuL       = AuxuL;          SD->AuxuR       = AuxuR;
  //SD->AuxU        = AuxU;      
  SD->EG          = EG;
  //if (MotionOn){
  //  SD->MDL = MDL;
  //  SD->MDR = MDR;
  //}
  //SD->ResPhiDataL  = ResPhiDataL;   SD->ResPhiDataR = ResPhiDataR;
  //SD->ResPhiTable  = ResPhiTable;

  if (pSD == NULL){
    // Delete StaticData that we just created
    ierr = xf_Error(xf_DestroyStaticDataIFace(SD));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualIFaces
static int
xf_CalculateResidualIFaces(xf_All *All, xf_Vector *U, xf_Vector *R, 
			   xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
/*
PURPOSE: 

  Calculates the residuals and residual Jacobians (if requested)
  associated with the left and right elements due to integration on
  all interior faces

INPUTS: 
 
  All: All structure
  U : state vector on all elements
  SolverData : solver data structure

OUTPUTS:

  R : residual vector
  R_U : Jacobian matrix (optional: can be null, or R_U->Value can be null)

RETURNS: Error code

*/

  int ierr;
  int nIFaceRegular, iiface;
  int egrpL, egrpR, elemL, elemR, faceL, faceR;
  int ReturnError = xf_OK;
  real *RL = NULL, *RR = NULL;
  real *RL_UL = NULL, *RL_UR = NULL, *RR_UR = NULL, *RR_UL = NULL;
  enum xfe_Bool MeshIsParallel;
  enum xfe_Bool UseGCL;
  enum xfe_Bool MotionOn;
  xf_StaticDataIFace *StaticData = NULL;
  xf_Vector *GCL = NULL, *GammaVec;
  xf_IFace IFace;
  xf_Mesh   *Mesh;
  xf_Data   *GammaDat;
  ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
  if(ierr == xf_NOT_FOUND) {xf_printf("Cannot find heat capacity ratio...\n"); return ierr;}
  else
     GammaVec = (xf_Vector *) GammaDat->Data;

  Mesh = All->Mesh;
  
  nIFaceRegular = -1;
  if (MeshIsParallel = (Mesh->ParallelInfo != NULL))
    nIFaceRegular = Mesh->ParallelInfo->nIFaceRegular;

  // determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  
  // Is mesh motion on?
  MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));
  
  if (UseGCL && MotionOn){
    // find primary GCL vector
    ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
    if (ierr != xf_OK) return ierr;
  }

  // Loop over interior faces
  for (iiface=0; iiface<Mesh->nIFace; iiface++){

    if (MeshIsParallel && (iiface == nIFaceRegular)){
      // wait for end communication of halo data: state
      ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
      if (ierr != xf_OK) return ierr;

      if (UseGCL && MotionOn){
	// same with GCL
	ierr = xf_Error(xf_HaloExchangeVectorEnd(GCL));
	if (ierr != xf_OK) return ierr;
      }
    }

    IFace = Mesh->IFace[iiface];
    egrpL = IFace.ElemGroupL;
    egrpR = IFace.ElemGroupR;
    elemL = IFace.ElemL;
    elemR = IFace.ElemR;
    faceL = IFace.FaceL;
    faceR = IFace.FaceR;

    /* Prepare pointers to residual and Jacobian */
    RL = R->GenArray[egrpL].rValue[elemL];
    RR = R->GenArray[egrpR].rValue[elemR];

    if ((R_U == NULL) || (R_U->Value == NULL)){
      RL_UL = RR_UR = RL_UR = RR_UL = NULL;
    }
    else{
      RL_UL = R_U->Value[egrpL][elemL][0];
      RL_UR = R_U->Value[egrpL][elemL][1+faceL];
      RR_UR = R_U->Value[egrpR][elemR][0];
      RR_UL = R_U->Value[egrpR][elemR][1+faceR];
    }
    
    /* Calculate residual on IFace, passing in StaticData  */
    ierr = xf_Error(xf_CalculateResidualIFace(All, iiface, U, RL, RR, RL_UL, RR_UR,
					      RL_UR, RR_UL, &StaticData, SolverData));
    if (ierr != xf_OK){  // for parallel's sake, do not return immediately if recoverable
      if (!xf_CheckRecoverable(ierr, &ReturnError)) return ierr;
      continue; // recoverable error occured, move on
    }

  } // iiface

  // Delete StaticData
  ierr = xf_Error(xf_DestroyStaticDataIFace(StaticData));
  if (ierr != xf_OK) return ierr;

  return ReturnError;
}


/******************************************************************/
//   FUNCTION Definition: xf_SortEqnSetBCs
int
xf_SortEqnSetBCs(const xf_Mesh *Mesh, xf_BCs *BCs)
{
  int ierr, nBC, nbfgrp, iBC, ibfgrp;
  enum xfe_Bool found;
  xf_BC *BCnew;

  nBC = BCs->nBC;
  nbfgrp = Mesh->nBFaceGroup;

  if (nBC != nbfgrp) return xf_Error(xf_BOUNDARY_CONDITION_ERROR);
  
  ierr = xf_Error(xf_Alloc((void **) &BCnew, nBC, sizeof(xf_BC)));
  if (ierr != xf_OK) return ierr;

  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
    found = xfe_False;
    for (iBC=0; iBC<nBC; iBC++)
      if (strcmp(Mesh->BFaceGroup[ibfgrp].Title, BCs->BC[iBC].BFGTitle) == 0){
	BCnew[ibfgrp] = BCs->BC[iBC];
	found = xfe_True;
      }
    if (!found){
      xf_printf("Error. Was not able to find Mesh bfg %s in the EqnSet BC list.\n", 
		Mesh->BFaceGroup[ibfgrp].Title);
      return xf_Error(xf_BOUNDARY_CONDITION_ERROR);
    }
  }
  
  xf_Release( (void *) BCs->BC);
  BCs->BC = BCnew;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateStaticDataBFace
static int
xf_CreateStaticDataBFace(xf_All *All, xf_SolverData *SolverData, xf_StaticDataBFace **pSD)
{
/*
PURPOSE: 

  Creates and initializes bface static data

INPUTS: 
 
  All : all structure to determine residual terms
  SolverData : solver data structure

OUTPUTS:

  (*pSD) : static data that is created and initialized

RETURNS: Error code

*/
  int ierr, iAux;
  xf_StaticDataBFace *SD;

  // allocate memory
  ierr = xf_Error(xf_Alloc( (void **) pSD, 1, sizeof(xf_StaticDataBFace)));
  if (ierr != xf_OK) return ierr;
  SD = (*pSD);
  
  //ierr = xf_Error(xf_ProcessResTerms(All->EqnSet->ResTerms, SD->iConv, 
  //	 			     SD->iDiff, SD->iSource));
  //if (ierr != xf_OK) return ierr;
  
  // initialize variables to NULL
  SD->QuadData    =  NULL;
  SD->PhiData     =  NULL; 
  SD->NData       =  NULL;
  SD->JData       =  NULL; 
  SD->GeomPhiData =  NULL;
  SD->wn	  =  NULL;
  SD->xglob       =  NULL;
  SD->xelem       =  NULL;
  SD->uI          =  NULL; 
  SD->guI	  =  NULL; 
  SD->uB          =  NULL; 
  //SD->uB_uI       =  NULL; 
  SD->F	          =  NULL; 
  //SD->F_uI        =  NULL; 
  //SD->CF          =  NULL;
  SD->DData	  =  NULL;
  SD->pnq         =  -1; 
  SD->nqConv      =  -1; 
  SD->nqDiff      =  -1; 
  SD->nnDiff      =  -1; 
  SD->EG          =  NULL;
  //SD->EM          =  NULL;
  //SD->StabVisc    =  NULL;  
  //SD->StabPhi     =  NULL;  
  //SD->ResMetric   =  NULL;
  //SD->ResPhiData  =  NULL;
  //SD->ResPhiTable =  NULL;

  // determine Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &SD->Time));
  if (ierr != xf_OK) return ierr; 

  // build eqnset-desired parameter lists for passing into functions
  /*
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &SD->IParam, &SD->RParam, 
				       &SD->nAuxU, &SD->AuxU));
  if (ierr != xf_OK) return ierr;
*/
  // create structure for storing jump data 
  ierr = xf_Error(LYDG_CreateDiffBCData(&SD->DData));
  if (ierr != xf_OK) return ierr;
/*
  // Allocate vector of basis data for interpolating auxiliary vectors
  ierr = xf_Error(xf_Alloc( (void **) &SD->AuxPhiData, SD->nAuxU, sizeof(xf_BasisData *)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &SD->Auxu, SD->nAuxU, sizeof(real *)));
  if (ierr != xf_OK) return ierr;
  for (iAux=0; iAux<SD->nAuxU; iAux++){
    SD->AuxPhiData[iAux] = NULL;
    SD->Auxu[iAux]       = NULL;
  }*/

  /* Create a basis table, PhiTable, that will store computed basis
     functions specific to each [element shape, face in element,
     orientation of face] combination, for quick lookup. */
  ierr = xf_Error(xf_CreateBasisTable(&SD->PhiTable));
  if (ierr != xf_OK) return ierr;

  /* Create another basis table, GeomPhiTable, for storing the
     geometry approximation basis functions evaluated at
     faces/orientations/shapes. */
  ierr = xf_Error(xf_CreateBasisTable(&SD->GeomPhiTable));
  if (ierr != xf_OK) return ierr;

  // find element geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &SD->EG));
  if (ierr != xf_OK) return ierr;

/*
  if ( (SolverData != NULL) && (SolverData->StabRequired) ){
    // obtain elem metric vector
    ierr = xf_Error(xf_FindElemHMetric(All, xfe_True, &SD->EM));
    if (ierr != xf_OK) return ierr;
  }
  
  // initialize motion data
  SD->MD = NULL;
  if ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active)){
    ierr = xf_Error(xf_CreateMotionData(All, &SD->MD));
    if (ierr != xf_OK) return ierr;
  }

  // initialize data structure for residual linearization
  ierr = xf_Error(xf_InitLinQueue(SD->LinQ));
  if (ierr != xf_OK) return ierr;
*/
  
  /* Create another basis table, ResPhiTable, that will store computed
     basis functions specific to each [element shape, face in element,
     orientation of face] at an incremented order. */
  //if ( (SolverData != NULL) && (SolverData->ResidualOrderIncrement != 0)) {
  //  ierr = xf_Error(xf_CreateBasisTable(&SD->ResPhiTable));
  //  if (ierr != xf_OK) return ierr;
  //}


  // everything must go back with (*pSD)
  (*pSD) = SD;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyStaticDataBFace
static int
xf_DestroyStaticDataBFace(xf_StaticDataBFace *SD)
{
/*
PURPOSE: Destroys static iface structure

INPUTS: 
 
  SD : static data

OUTPUTS: None, SD and all of its contents are destroyed

RETURNS: Error code

*/

  int ierr, iAux;
  
  if (SD == NULL) return xf_OK;

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(SD->QuadData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(SD->PhiData, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Table */
  ierr = xf_Error(xf_DestroyBasisTable(SD->PhiTable));
  if (ierr != xf_OK) return ierr;

  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(SD->GeomPhiData, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Destroy Geometry Basis Table */
  ierr = xf_Error(xf_DestroyBasisTable(SD->GeomPhiTable));
  if (ierr != xf_OK) return ierr;

  /* Destroy Normal Data */
  ierr = xf_Error(xf_DestroyNormalData(SD->NData));
  if (ierr != xf_OK) return ierr;

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(SD->JData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Diff Data */
  xf_DestroyDiffBCData(SD->DData);

  /* Destroy Auxiliary Vector Data */
  /*for (iAux=0; iAux<SD->nAuxU; iAux++){
    ierr = xf_Error(xf_DestroyBasisData(SD->AuxPhiData[iAux], xfe_True));
    if (ierr != xf_OK) return ierr;
    xf_Release( (void *) SD->Auxu[iAux]);
  }
  xf_Release( (void *) SD->AuxPhiData);
  xf_Release( (void *) SD->Auxu);
  xf_Release( (void *) SD->AuxU);
*/
  /* Destroy mesh motion data */
  //xf_DestroyMotionData(SD->MD);

  // Destroy residual linearization structure
  //xf_DestroyLinQueue(SD->LinQ);
  
  /* Destroy Residual Basis Data */
  //ierr = xf_Error(xf_DestroyBasisData(SD->ResPhiData, xfe_False));
  //if (ierr != xf_OK) return ierr;

  /* Destroy Residual Basis Table */
  //ierr = xf_Error(xf_DestroyBasisTable(SD->ResPhiTable));
  //if (ierr != xf_OK) return ierr;


  // Release memory
  //xf_Release( (void *) SD->IParam);
  //xf_Release( (void *) SD->RParam);
  xf_Release( (void *) SD->wn);
  xf_Release( (void *) SD->xglob);
  xf_Release( (void *) SD->xelem);
  xf_Release( (void *) SD->uI);
  xf_Release( (void *) SD->guI);
  xf_Release( (void *) SD->uB);
  //xf_Release( (void *) SD->uB_uI);
  xf_Release( (void *) SD->F);
  //xf_Release( (void *) SD->F_uI);
  //xf_Release( (void *) SD->CF);
  //xf_Release( (void *) SD->StabVisc);
  //xf_Release( (void *) SD->StabPhi);
  //xf_Release( (void *) SD->ResMetric);

  // Destroy self
  xf_Release( (void *) SD);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ModMotionConvFBC
static void
xf_ModMotionConvFBC(real *uB, real *uB_uI, real *wn, int nq, int dim, 
		    int sr, xf_MotionData *MD, real *F, real *F_uI)
{
/*
PURPOSE: 

  Modifies convective BC linearization F_uI (if not NULL) to take into
  account motion of the mesh prescribed in the structure MD.  The flux
  is already modified in the eqn-set specific function.

INPUTS: 
 
  uB   : boundary state in physical space [nq*sr]
  uB_uI: linearization of bc state wrt interior state [nq*sr2]
  wn   : normal vector in physical space [nq*dim]
  nq   : number of points
  dim  : dimension
  sr   : state rank
  MD   : mesh motion data structure
  F    : convective BC flux [nq*sr]
  F_uI : convective BC flux linearization w.r.t the interior state [nq*sr*sr]

OUTPUTS:

  F, F_uI : modified versions

RETURNS: Error code
*/

  int iq, d, k, a, i, i2, sr2;
  real gb;

  sr2  = sr*sr;

  /* multiply F_uI by 1/gb */
  if (F_uI != NULL){
    for (iq=0; iq<nq; iq++)
      for (k=0, gb=MD->gb[iq], i=iq*sr2; k<sr2; k++) F_uI[i+k] /= gb;
  }

}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualBFace
static int
xf_CalculateResidualBFace(xf_All *All, int ibfgrp, int ibface, const xf_Vector *U, 
			  real *ER, real *ER_U, xf_OutputEvalData *OutputEval, 
			  xf_StaticDataBFace **pSD, xf_SolverData *SolverData)
{
/*

PURPOSE:

  Calculates the boundary contribution to the spatial residual, ER,
  due to integration on BFace = ibfgrp, ibface.  The residual
  Jacobian, R_U, is computed if not passed as NULL.  In addition a
  flux output, Value, can be calculated if not passed in as NULL.
  This flux output consists of the integrated convective and diffusive
  fluxes weighted by a given FluxWeights vector (size StateRank) to
  form a scalar.


INPUTS:

  All: All structure
  ibfgrp, ibface : boundary face in question
  U : state vector
  OutputEval : output evaluation structure.  If not NULL, output
               evaluation is assumed.
  pSD : pointer to static data (optional, can pass in as NULL)
        useful for avoiding reallocations when calling multiple times
  SolverData: contains pointers to solver variables such as the
              connectivity vector or the regularity vector.  Some
	      of these (e.g. connectivity) may be modified.

OUTPUTS:

  ER : Residual on element adjacent to BFace
  ER_U : associated residual Jacobian

RETURN:

  Error code

*/

  int ierr, sr, sr2, i, j, k, l, iq, nq, nn, d, n;
  int dim, Order, QuadOrder, ResOrder;
  int egrp, elem, face;
  int pnq, nqConv, nqDiff, nnDiff;
  int ResidualOrderIncrement, NewOrder, nnmax;
  int  nConv,  nDiff,  nSource;
  int *iConv, *iDiff, *iSource;
  int iAux, nAuxU;
  int *IParam;
  int *FluxMoments = NULL;
  int localBCflag;
  enum xfe_BasisType Basis;
  enum xfe_DiffDiscType DiffDisc;
  enum xfe_Bool QuadChanged;
  enum xfe_Bool Need_guI = xfe_False;
  enum xfe_Bool Need_Grad;
  enum xfe_Bool SkipDiffusion, SkipDiffStab;
  enum xfe_Bool StabRequired;
  enum xfe_Bool MotionOn;
  char Title[xf_MAXSTRLEN];
  real cval, ElemVol, **Auxu = NULL;
  real *RParam, *xq, *wq, *wn, *xglob, *xelem, *CF;
  real *EU, *uI, *guI, *uB, *uB_uI, *F, *F_uI, Gamma;
  real *Value = NULL, *EV_U = NULL, *EV_G = NULL, *FluxWeights = NULL;
  real xfac, rtemp, Nmag, Time;
  xf_ResTerm *ResTerm;
  xf_BFace BFace;
  xf_BC *BC;
  xf_QuadData *QuadData;
  xf_BasisTable *PhiTable, *GeomPhiTable;
  xf_BasisData *PhiData, *GeomPhiData, **AuxPhiData;
  xf_BasisTable *ResPhiTable;
  xf_BasisData *ResPhiData;
  xf_NormalData *NData;
  xf_JacobianData *JData;
  xf_DiffBCData *DData;
  xf_Vector *EG, *C, *V, **AuxU, *GammaVec;
  xf_Data *GammaDat;
  xf_MotionData *MD = NULL;
  xf_LinQueueData *LinQ = NULL;
  xf_StaticDataBFace *SD = NULL;
  xf_StabData *StabData = NULL;
  xf_Mesh   *Mesh;
  //xf_EqnSet *EqnSet;
  FILE *fidDump = NULL;

  //added in Feb 2014
  real MaxCharSpeed;
  xf_Data *MaxC_Dat;
  xf_Vector *MaxC_Vec;

  // General information
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  sr   = Model->nVars;
  sr2 = sr*sr;
  BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];

  // output information
  if (OutputEval != NULL){
    Value       = OutputEval->Value;
    EV_U        = OutputEval->EV_U;
    EV_G        = OutputEval->EV_G;
    FluxWeights = OutputEval->FluxWeights;
    FluxMoments = OutputEval->FluxMoments;
    fidDump     = OutputEval->fidDump;
  }

  // Create and initialize static data if not passed in
  if ((pSD == NULL) || ((*pSD) == NULL)){
    ierr = xf_Error(xf_CreateStaticDataBFace(All, SolverData, (pSD != NULL) ? pSD : &SD));
    if (ierr != xf_OK) return ierr;
  }
  if (pSD != NULL) SD = (*pSD);

   
  //ResTerm = EqnSet->ResTerms->ResTerm;

  // Are we doing mesh motion?
  //MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));

  //find Gamma Vector
  ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
  if(ierr == xf_NOT_FOUND)
  {
     xf_printf("Cannot find heat capacity ratio...\n");
     return ierr;
  }
  else
     GammaVec = (xf_Vector *) GammaDat->Data;

  //find the data/vector for maximum characteristic speed
  ierr = xf_FindDataByTitle(All->DataSet, "ElemMaxCharSpeed", xfe_Vector, &MaxC_Dat);
  if(ierr == xf_NOT_FOUND)
  {
     xf_printf("Cannot find maximum characteristic speed...\n");
     return ierr;
  }
  else
     MaxC_Vec = (xf_Vector *) MaxC_Dat->Data;

  // Pull off variables from StaticData  
  /*iConv   = SD->iConv;
  iDiff   = SD->iDiff;
  iSource = SD->iSource;
  nConv   =   iConv[1]-  iConv[0];
  nDiff   =   iDiff[1]-  iDiff[0]; 
  nSource = iSource[1]-iSource[0];
  */
  QuadData    = SD->QuadData;
  PhiData     = SD->PhiData;
  PhiTable    = SD->PhiTable;
  NData       = SD->NData;
  JData       = SD->JData;
  GeomPhiData = SD->GeomPhiData;
  GeomPhiTable= SD->GeomPhiTable;
  wn          = SD->wn;
  xglob       = SD->xglob;
  xelem       = SD->xelem;
  uI          = SD->uI;
  guI         = SD->guI;
  uB          = SD->uB;
  //uB_uI       = SD->uB_uI;
  F           = SD->F;
  //F_uI        = SD->F_uI;
  //CF          = SD->CF;
  DData       = SD->DData;
  pnq         = SD->pnq;        
  nqConv      = SD->nqConv;     
  nqDiff      = SD->nqDiff;
  nnDiff      = SD->nnDiff;
  //IParam      = SD->IParam;
  //RParam      = SD->RParam;
  //nAuxU       = SD->nAuxU;      
  //AuxPhiData  = SD->AuxPhiData;
  //Auxu        = SD->Auxu;     
  //AuxU        = SD->AuxU;      
  EG          = SD->EG;
  Time        = SD->Time;       // simulation time
  //LinQ        = SD->LinQ;       // residual linearization structure
  //if (MotionOn)
  //  MD        = SD->MD;         // mesh motion data
  //ResPhiData  = SD->ResPhiData;
  //ResPhiTable = SD->ResPhiTable;

  // if we have diffusion terms
  /*if (nDiff > 0){
    Need_guI = xfe_True; // we will need the gradient, and we need the type of discretization
    ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "DiffusionDiscretization", 
				       xfe_DiffDiscName, (int ) xfe_DiffDiscLast, 
				       (int *) &DiffDisc));
    if (ierr != xf_OK) return ierr;
  }

  // Set Connectivity vector if required
  if ((SolverData != NULL) && (SolverData->CRequired))
    C = SolverData->C;
  else
    C = NULL;

  // Stabilization required?
  StabRequired = ((SolverData != NULL) && (SolverData->StabRequired));

  // Residual order increase
  ResidualOrderIncrement = ((SolverData == NULL) ? 0 : SolverData->ResidualOrderIncrement);
*/

  // Sort EqnSet->BCs to match boundary face groups
  //ierr = xf_Error(xf_SortEqnSetBCs(Mesh, EqnSet->BCs+0));
  //if (ierr != xf_OK) return ierr;
  //BC = EqnSet->BCs[0].BC;
  //match up Yu's new BC structure in Yu_Model
  for(i=0; i<Model->nBCs; i++)
     if(strcmp(Mesh->BFaceGroup[ibfgrp].Title, Model->nameBCs[i]) == 0)
     {  localBCflag = i; break;  }

  if(i==Model->nBCs) {
     printf("Boundary name does not match model defining.\n");
     xf_Error(xf_BOUNDARY_CONDITION_ERROR);
  }
  
  /* Determine whether we need to calculate gradients, _U */
  Need_Grad = ((ER_U != NULL) || (EV_U != NULL) || (EV_G != NULL));

 
  egrp = BFace.ElemGroup;
  elem = BFace.Elem;
  face = BFace.Face;
      
  Basis = U->Basis[egrp];
  Order = xf_InterpOrder(U, egrp, elem);

  //if found, directly get the gamma value at this element
  Gamma = GammaVec->GenArray[egrp].rValue[elem][0];
  
  // Residual order
  //ResOrder = (ResidualOrderIncrement == 0) ? Order : max(Order + ResidualOrderIncrement, 0);
  
  // determine required integration order
  ierr = xf_Error(xf_GetQuadOrderBFace(Mesh, NULL, BFace, Order, &QuadOrder));
  if (ierr != xf_OK) return ierr;
  
  if(!Model->TwistFlag)
     QuadOrder -= (Mesh->Dim - 1);

  /* Pull off quad points for the bface; will not recalculate if
     Basis/Order have not changed. */
  ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face, QuadOrder, 
			      &QuadData, &QuadChanged));
  if (ierr != xf_OK) return ierr;
  
  nq = QuadData->nquad;
  xq = QuadData->xquad;
  wq = QuadData->wquad;
      
  // compute basis functions if quad or basis or order changed
  ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrp, elem, face, BFace.Orient,
					       Basis, Order, QuadChanged, nq, xq, 
					       xfb_Phi | xfb_GPhi | xfb_gPhi, 
					       &PhiData, PhiTable, &xelem));
  if (ierr != xf_OK) return ierr;

  // also compute lower/higher order basis functions for residual eval
  /*if (ResidualOrderIncrement != 0){
    if (ResPhiTable == NULL) return xf_Error(xf_CODE_LOGIC_ERROR); // sanity check
    NewOrder = max(Order + ResidualOrderIncrement, 0);
    ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrp, elem, face, BFace.Orient,
                                                 Basis, NewOrder, QuadChanged, nq, xq, 
                                                 xfb_Phi, &ResPhiData, ResPhiTable, NULL));
    if (ierr != xf_OK) return ierr;
  }*/


  /* Compute normal(s) at quad points.  If face is straight, only
     one normal will be computed/returned. */
  ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, nq, xq, &NData, NULL));
  if (ierr != xf_OK) return ierr;

  if (Model->DiffFlag){
    /* Compute geometry Jacobian; if not constant, compute at quad
       points.  Note if jacobian is constant, only one Jacobian will
       be computed/returned. */
    ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xelem, 
				    xfb_detJ | xfb_iJ, xfe_True, &JData));
    if (ierr != xf_OK) return ierr;
       
    // convert reference basis grads (GPhi) to physical grads, gPhi
    ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
    if (ierr != xf_OK) return ierr;

  }

  nn = PhiData->nn;
  nnmax = nn = PhiData->nn;
  //if (ResidualOrderIncrement != 0)  // account for possibly-different residual order
  //  nnmax = max(nn, ResPhiData->nn);
  
  EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]

  // re-allocate data if quad points increased
  if (nq > pnq){
    ierr = xf_Error(xf_ReAlloc( (void **) &wn, nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &uI, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &uB, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    if (Model->DiffFlag){
      ierr = xf_Error(xf_ReAlloc( (void **) &guI, nq*sr*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
/*    if (Need_Grad){
      ierr = xf_Error(xf_ReAlloc( (void **) &uB_uI, nq*sr2, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    if (C != NULL){
      ierr = xf_Error(xf_ReAlloc( (void **) &CF, nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    if (StabRequired){
      ierr = xf_Error(xf_ReAlloc( (void **) &SD->StabVisc, nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc( (void **) &SD->StabPhi, nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReAlloc( (void **) &SD->ResMetric, nq*dim*dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }*/

  }
      
  // pull off element volume
  ElemVol = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];

  // obtain global coords of quad points
  ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, xfe_True,
				  nq, xelem, xglob));
  if (ierr != xf_OK) return ierr;
/*   ierr = xf_Error(xf_Ref2GlobFaceUsingTable(Mesh, egrp, elem, face, BFace.Orient, */
/* 					    QuadChanged, &GeomPhiData, GeomPhiTable, */
/* 					    nq, xelem, xglob)); */
/*   if (ierr != xf_OK) return ierr; */


  // obtain transformation map if doing mesh motion
  /*if (MotionOn){
    ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, All->Mesh->Motion, 
				      nq, dim, Time, xglob, MD));
    if (ierr != xf_OK) return ierr;
  }*/
      
  // interpolate state and gradient at quad points
  xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, uI);
  if (Model->DiffFlag){
    for (d=0; d<dim; d++)
      xf_MxM_Set(PhiData->gPhi+nn*nq*d, EU, nq, nn, sr, guI+nq*sr*d);
  }
/*
  // interpolate any auxiliary vectors, using AuxPhiData
  for (iAux=0; iAux<nAuxU; iAux++){
    V  = AuxU[iAux];
    if (nq > pnq){ // reallocate Auxu[iAux] if necessary
      ierr = xf_Error(xf_ReAlloc( (void **) Auxu+iAux, nq*V->StateRank, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_EvalBasis(V->Basis[egrp], xf_InterpOrder(V,egrp,elem), xfe_True, 
				 nq, xelem, xfb_Phi, AuxPhiData+iAux));
    if (ierr != xf_OK) return ierr;
    xf_MxM_Set(AuxPhiData[iAux]->Phi, V->GenArray[egrp].rValue[elem], nq, 
	       AuxPhiData[iAux]->nn, V->StateRank, Auxu[iAux]);
  }*/


  /* Initially, set wn = just the normals.  Will include quad
     weights after call to EqnSetBCState, since some weights may
     be negative while EqnSetBCState needs to have the true
     outward-pointing normal at each point.
  */
  for (d=0; d<dim; d++)
    for (iq=0;iq<nq; iq++) 
      wn[iq*dim+d] = NData->n[iq*dim*(NData->nq!=1)+d];

  // transform state and normal to physical space
  //if (MotionOn) xf_ModMotionPreEqnCall(nq, dim, sr, MD, uI, wn);
  
  // obtain boundary state and derivative at quad points
  //ierr = xf_Error(xf_EqnSetBCState(EqnSet, BC+ibfgrp, IParam, RParam, nq, 
  //				   wn, xglob, &Time, (MotionOn) ? MD->vg : NULL, 
  //				   uI, uB, uB_uI));
  //ierr = xf_Error(ConvFluxBoundaryState(nq, sr, dim, uI, uB, wn, xglob, Time, Model->typeBCs[localBCflag],
  //                                      Model->paraBCs[localBCflag], Model->BCFunSpf, Gamma));
  ierr = xf_Error(ConvFluxBoundaryState(Model, nq, uI, uB, wn, xglob, Time, localBCflag, Gamma));
  if (ierr != xf_OK) return ierr;

  // transform state and normal back to reference space
  //if (MotionOn){
  //  xf_ModMotionPostEqnCall(nq, dim, sr, MD, uI, wn);
  //  if (uB_uI != NULL) xf_ColDiv(uB_uI, MD->gb, nq, sr*sr, 1);  // uB_uI now w.r.t ref uI
  //}
  
  // multiply wn by quad weights for further use in integration
  for (iq=0;iq<nq; iq++) 
    for (d=0; d<dim; d++)
      wn[iq*dim+d] *= wq[iq];
           
  // zero out Connectivity
  //if (C != NULL) C->GenArray[egrp].rValue[elem][face] = 0.;


  // pull off data required for stabilization
  /*if ( (SolverData != NULL) && (SolverData->StabRequired) ){
    StabData = &(SolverData->StabData);
    ierr = xf_Error(xf_CalculateStabViscBFace(All, ibfgrp, ibface, nq, xelem, 
					      SD->EM->GenArray[egrp].rValue[elem],
					      StabData, SD->StabVisc, &StabData->StabVisc_U, 
					      SD->StabPhi, SD->ResMetric, &SkipDiffStab));
    if (ierr != xf_OK) return ierr;
    StabData->StabVisc     = SD->StabVisc;
    StabData->StabPhi      = SD->StabPhi;
    StabData->ResMetric    = SD->ResMetric;
    StabData->SkipDiffStab = SkipDiffStab;
  }
  else{
    StabData = NULL;
    SkipDiffStab = xfe_False;
  }*/
   
  /* Clear residual linearization queue */
  //if (ER_U != NULL) xf_ClearLinQueue(LinQ);

  /*------------------*/
  /* CONVECTION TERMS */
  /*------------------*/
  //if (nConv > 0){
	
    if (nq > nqConv){ // realloc F, F_uI if necessary
      nqConv = nq;
      ierr = xf_Error(xf_ReAlloc( (void **) &F, sr*nq, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    //  if (Need_Grad){
	//ierr = xf_Error(xf_ReAlloc( (void **) &F_uI, sr*sr*nq, sizeof(real)));
	//if (ierr != xf_OK) return ierr;
    //  }
    }
	
    // transform state and normal to physical space
    //if (MotionOn){
    //  xf_ModMotionPreEqnCall(nq, dim, sr, MD, uI, wn);
    //  if (uB_uI != NULL) xf_ColMult(uB_uI, MD->gb, nq, sr*sr, 1); // uB_uI is now w.r.t phys uI
    //}

    // calculate F [and F_uI] at quad points
    //ierr = xf_Error(xf_EqnSetConvFBC(EqnSet, ResTerm+iConv[0], nConv, BC+ibfgrp, IParam, 
    //				     RParam, nq, uI, uB, uB_uI, Auxu, wn, xglob, 
    //				     (MotionOn) ? MD->vg : NULL, F, F_uI, CF));
    ierr = xf_Error(ConvFluxBoundaryFace(nq, sr, dim, uI, uB, wn, F, Gamma, &MaxCharSpeed));
    if (ierr != xf_OK) return ierr;

    //update MaxCharSpeed for boundary element
    if(MaxC_Vec->GenArray[egrp].rValue[elem][0] < MaxCharSpeed)
       MaxC_Vec->GenArray[egrp].rValue[elem][0] = MaxCharSpeed;

    // Modify fluxes if motion is on 
    //if (MotionOn) xf_ModMotionConvFBC(uB, uB_uI, wn, nq, dim, sr, MD, F, F_uI);

    // transform state and normal back to reference space
    //if (MotionOn){
    //  xf_ModMotionPostEqnCall(nq, dim, sr, MD, uI, wn);
    //  if (uB_uI != NULL) xf_ColDiv(uB_uI, MD->gb, nq, sr*sr, 1);  // uB_uI now w.r.t ref uI
    //}

    /* Add to R:
       ER_{n,k} += sum_q Phi_{n,q}*F_{q,k}
    */
    if (ER != NULL)
      xf_MTxM_Add(PhiData->Phi, F, nn, nq, sr, ER);
	
    /* Add to R_U:
       ER_U{n,k;m,a} += sum_q Phi_{n,q}*F_uI{q,k;a}*Phi_{m,q}
    */
    //if (ER_U != NULL){
    //  ierr = xf_Error(xf_AddToLinQueue(F_uI, xfe_LinQTerm_PhiPhi, nq, dim, 
	//			       sr2, -1, NULL, 0, 1.0, LinQ));
    //  if (ierr != xf_OK) return ierr;
    //}
	
    //if (C != NULL){ // add to connectivity; note, wquad included via wn above
    //  for (iq=0, cval=0.; iq<nq; iq++) cval += CF[iq];
    //  C->GenArray[egrp].rValue[elem][face] += fabs(cval);
    //}


    // OUTPUT
    /* Add F term to Value:
       Value += FluxWeights{k}*F{q,k}[*xglob{q,FluxMoments{k}}], sum over q,k
    */
    if (Value != NULL)
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++){
	  xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
	  (*Value) += FluxWeights[k]*F[iq*sr+k]*xfac;
	}
	
    /* Add F_uI to EV_U:
       EV_U{n,l} += FluxWeights{k}*F_uI{q,k;l}*Phi{q,n}[*xglob{q,FluxMoments{k}}]
       EV_G{n}   += FluxWeights{k}*F_uI{q,k;l}*(-uI{q,l})*Phi{q,n}[*xglob{q,FluxMoments{k}}]
       Note: above summed over q,k
    */
    /*
    if (EV_U != NULL){
      // xf_MTxM_Set(PhiData->Phi, F_uI, nn, nq, sr2, T2);
      for (n=0; n<nn; n++)
	for (k=0; k<sr; k++)
	  for (l=0; l<sr; l++)
	    for (iq=0; iq<nq; iq++){
	      xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
              rtemp = FluxWeights[k]*F_uI[iq*sr2+k*sr+l]*PhiData->Phi[iq*nn+n]*xfac;
	      EV_U[n*sr+l] += rtemp;
              if (EV_G != NULL) EV_G[n] += rtemp*(-uI[iq*sr+l]/MD->gb[iq]); // GCL linearization
              // note, in above with MM on, uI is the reference state
	    }
    }*/
	  
  //} // end if nConv > 0

  /*-----------------*/
  /* DIFFUSION TERMS */
  /*-----------------*/

  // Skip diffusion terms if only have a stabilization term and
  //   the regularity flag is 0. 
  //SkipDiffusion = ((nDiff == 1) && SkipDiffStab);

  if(Model->DiffFlag && !Model->AVmodel){
  //if(Model->DiffFlag){ // && !Model->AVmodel){

//  if ((nDiff > 0) && (!SkipDiffusion)){

    if ((nq > nqDiff) || (nnmax > nnDiff)){  // realloc DiffData if necesasry
      nqDiff = max(nq,    nqDiff); 
      nnDiff = max(nnmax, nnDiff);
      ierr = xf_Error(xf_ReAllocDiffBCData(sr, nqDiff, dim, nnDiff, 
					   xfe_False, DData));
      if (ierr != xf_OK) return ierr;
    }
    
    // compute Q flux; in the process DData is filled in
/*
    //modify BC here
    ierr = xf_Error(xf_DiffFluxQBC(All, ibfgrp, ibface, ResTerm+iDiff[0], nDiff, 
				   DiffDisc, IParam, RParam, PhiData, 
                                   (ResidualOrderIncrement == 0) ? PhiData : ResPhiData,
                                   StabData, ElemVol, Basis, Order, nq, wn, xglob, uI, uB, 
				   uB_uI, guI, MD, OutputEval, ER, ER_U, LinQ, DData, CF));
    if (ierr != xf_OK) return ierr;
*/
    ierr = xf_Error(LYDG_BoundaryViscTerm(All, Model, ibfgrp, ibface, PhiData, PhiData,
                                          ElemVol, Basis, Order, nq, wn, xglob, uI, uB, 
                                          guI, ER, DData, Gamma, OutputEval));
    if (ierr != xf_OK) return ierr;
	
    // Add AB*(u-uhat)*n term to ER:
    //
    //   ER_{n,k} -= gPhi{i,q,n} * ABdun{i,q,k},  q,i summed
	//   
    //   Note: normal and quad weights are included in dun
   
    if (ER != NULL)
      for (d=0; d<dim; d++)
	xf_MTxM_Sub(PhiData->gPhi+nn*nq*d, DData->ABdun+sr*nq*d, nn, nq, sr, ER);

    // Add derivatives of AB*(u-uhat)*n term to R_U:
    //
    //   ER_U{n,k;m,a} -= gPhi{i,q,n} * [ AB_uIdun{i,q,k,a} * Phi{q,m}
    //   + AB{i,j,q,k,l}*du_uI{q,l,a}*Phi{q,m}*wn{q,j} ]
    //   Note: wn stores the normal with quad weights;
    //   du_uI is a scalar; q and i are summed
    
/*	
    if (ER_U != NULL){
      for (i=0; i<dim; i++){
	// Set AB{i,0,q,k,a} = sum_j AB{i,j,q,k,a}*wn{q,j}  (AB is overwritten here)
	xf_ColMult(DData->AB+nq*sr2*i*dim, wn+0, nq, sr2, dim); 
	for (j=1; j<dim; j++)
	  xf_ColMult_Add(DData->AB+nq*sr2*(i*dim+j), wn+j, nq, sr2, dim, DData->AB+nq*sr2*i*dim); 

	// multiply AB{i,0,q,k,a} by DData->du_uI (a matrix)
	xf_nMxM_Set(nq, DData->AB+nq*sr2*i*dim, DData->du_uI, sr, sr, sr, DData->T);
	
	// Contribution of AB (now in DData->T) to R_U linearization
	ierr = xf_Error(xf_AddToLinQueue(DData->T, xfe_LinQTerm_GPhiPhi, 
					 nq, dim, sr2, i, NULL, 0, -1.0, LinQ));
	if (ierr != xf_OK) return ierr;
      } // i
      if (!DData->ConstAB){
	ierr = xf_Error(xf_AddToLinQueue(DData->AB_uIdun, xfe_LinQTerm_GPhiPhi, 
					 nq, dim, sr2, -1, NULL, 0, -1.0, LinQ));
	if (ierr != xf_OK) return ierr;
      }
    } // end if R_U != NULL

    // Include linearization of stabilization term 
    if ( (ER_U != NULL) && (StabData != NULL) ){
      // assumption: StabVisc = 0 implies StabVisc_U = 0 (valid for continuous slope switch)
      // ER_U{n,k;m,a} -= gPhiL{i,n,q}*AwStab{i,q,k}*StabVisc_U{m,a}
      xf_MTxwM_Set(PhiData->gPhi+nn*nq*0, StabData->StabPhi, DData->AwStab+0*nq*sr, nn, nq, sr, DData->T);
      for (i=1; i<dim; i++)
	xf_MTxwM_Add(PhiData->gPhi+nn*nq*i, StabData->StabPhi, DData->AwStab+i*nq*sr, nn, nq, sr, DData->T);
      xf_BlockOutProd_Sub(DData->T, StabData->StabVisc_U, nn, sr, nn, ER_U);
    }
	
    if (C != NULL){ // add to connectivity
      for (iq=0, cval=0.; iq<nq; iq++) cval += CF[iq];
      C->GenArray[egrp].rValue[elem][face] += fabs(cval);
    }
*/    
  } // end if nDiff > 0


  // TODO: make this a separate function
  // Dump output data to a file if requested
/*  if (fidDump != NULL){
    // header (for each bface)
    fprintf(fidDump, "%% ibfgrp = %d (%s), ibface = %d, nq = %d\n", ibfgrp,
	    Mesh->BFaceGroup[ibfgrp].Title, ibface, nq);
    fprintf(fidDump, "%% %20s %20s", "x", "y");
    if (dim == 3) fprintf(fidDump, " %20s", "z");
    fprintf(fidDump, " %20s", "wq*Nmag");
    if (nConv > 0){
      for (k=0; k<sr; k++){
	sprintf(Title, "%s%s", "FConv:",  EqnSet->StateName[k]);
	fprintf(fidDump, " %20s", Title);
      }
    }
    if (nDiff> 0){
      for (k=0; k<sr; k++){
	sprintf(Title, "%s%s", "FDiff:",  EqnSet->StateName[k]);
	fprintf(fidDump, " %20s", Title);
      }
    }
    fprintf(fidDump, "\n");
    for (iq=0;iq<nq; iq++) { // loop over quad points
      fprintf(fidDump, " %20.12E %20.12E", xglob[dim*iq+0], xglob[dim*iq+1]);
      if (dim == 3) fprintf(fidDump, " %20.12E", xglob[dim*iq+2]);

      // Nmag = magnitude of normal
      for (d=0, Nmag=0.; d<dim; d++) 
	Nmag += NData->n[iq*dim*(NData->nq!=1)+d] * NData->n[iq*dim*(NData->nq!=1)+d];
      Nmag = sqrt(Nmag);

      fprintf(fidDump, " %20.12E", wq[iq]*Nmag);
      if (nConv > 0){
	for (k=0; k<sr; k++)
	  fprintf(fidDump, " %20.12E", F[iq*sr+k]/(wq[iq]*Nmag));
      }
      if (nDiff> 0){
	for (k=0; k<sr; k++)
	  fprintf(fidDump, " %20.12E", DData->Qn[iq*sr+k]/(wq[iq]*Nmag));
      }
      fprintf(fidDump, "\n");
    } // iq
  }
  
  // apply queued linearization
  if (ER_U != NULL){
    ierr = xf_Error(xf_ApplyLinQueue(LinQ, PhiData, PhiData, sr2, ER_U));
    if (ierr != xf_OK) return ierr;
  }
*/    
  pnq = nq; // set previous quad point # for next bface

  // Store possibly-altered or resized data back in StaticData
  SD->QuadData    = QuadData;
  SD->PhiData     = PhiData;
  SD->PhiTable    = PhiTable;
  SD->NData       = NData;
  SD->JData       = JData;
  SD->GeomPhiData = GeomPhiData;
  SD->GeomPhiTable= GeomPhiTable;
  SD->wn          = wn;
  SD->xglob       = xglob;
  SD->xelem       = xelem;
  SD->uI          = uI;
  SD->guI         = guI;
  SD->uB          = uB;
  //SD->uB_uI       = uB_uI;
  SD->F           = F;
  //SD->F_uI        = F_uI;
  //SD->CF          = CF;
  SD->DData       = DData;
  SD->pnq         = pnq;        
  SD->nqConv      = nqConv;     
  SD->nqDiff      = nqDiff;
  SD->nnDiff      = nnDiff;
  //SD->IParam      = IParam;
  //SD->RParam      = RParam;
  //SD->nAuxU       = nAuxU;      
  //SD->AuxPhiData  = AuxPhiData;
  //SD->Auxu        = Auxu;     
  //SD->AuxU        = AuxU;      
  SD->EG          = EG;       
  //if (MotionOn) SD->MD = MD;
  //SD->ResPhiData  = ResPhiData;
  //SD->ResPhiTable = ResPhiTable;

  if (pSD == NULL){
    // Delete StaticData that we just created
    ierr = xf_Error(xf_DestroyStaticDataBFace(SD));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualBFaces
int
xf_CalculateResidualBFaces(xf_All *All, const xf_Vector *U, xf_Vector *R, 
			   xf_JacobianMatrix *R_U, xf_OutputEvalData *OutputEval,
			   int nBFG, int *BFGs, xf_SolverData *SolverData)
{
  int ierr;
  int ibfg, nbfgrp, ibfgrp, ibface;
  int egrp, elem;
  int ReturnError = xf_OK;
  real *ER = NULL, *ER_U = NULL;
  xf_StaticDataBFace *StaticData = NULL;
  xf_BFace BFace;
  xf_Mesh *Mesh;


  Mesh = All->Mesh;

  /* OutputEval != NULL means output calculation is requested */
  if ((OutputEval != NULL) && (OutputEval->Value != NULL)) (*(OutputEval->Value)) = 0.0;

  // Loop over boundary groups and faces
  nbfgrp = ((BFGs == NULL) ? Mesh->nBFaceGroup : nBFG);
  
  for (ibfg=0; ibfg<nbfgrp; ibfg++){

    ibfgrp = ((BFGs == NULL) ? ibfg : BFGs[ibfg]);

    // skip zero-measure boundary groups
    if (strncmp(Mesh->BFaceGroup[ibfgrp].Title, "ZeroMeasure", 11) == 0) continue;

    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      
      BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];

      egrp = BFace.ElemGroup;
      elem = BFace.Elem;

      // Prepare pointers to residual, Jacobian, and output
      if (R != NULL)
	ER = R->GenArray[egrp].rValue[elem];
      if ((R_U != NULL) && (R_U->Value != NULL))
	ER_U = R_U->Value[egrp][elem][0];
      if (OutputEval != NULL){
	if (OutputEval->Value_U != NULL)
	  OutputEval->EV_U = (OutputEval->Value_U)->GenArray[egrp].rValue[elem];
	else
	  OutputEval->EV_U = NULL;
	if (OutputEval->Value_G != NULL)
	  OutputEval->EV_G = (OutputEval->Value_G)->GenArray[egrp].rValue[elem];
	else
	  OutputEval->EV_G = NULL;
      }
    
      /* Calculate residual on BFace, passing in StaticData  */
      ierr = xf_Error(xf_CalculateResidualBFace(All, ibfgrp, ibface, U, ER, ER_U, OutputEval, 
						&StaticData, SolverData));
      if (ierr != xf_OK){  // for parallel's sake, do not return immediately if recoverable
	if (!xf_CheckRecoverable(ierr, &ReturnError)) return ierr;
	continue; // recoverable error occured, move on
      }
    } // ibface
  } // ibfg


  // Delete StaticData
  ierr = xf_Error(xf_DestroyStaticDataBFace(StaticData));
  if (ierr != xf_OK) return ierr;

  return ReturnError;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidual
int
xf_CalculateResidual(xf_All *All, xf_Vector *U, xf_Vector *R, 
                       xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
    int ierr;
    enum xfe_Bool UseGCL;
    enum xfe_Bool MotionOn;
    int ReturnError = xf_OK;
    xf_Vector *C = NULL;
    xf_Vector *GCL = NULL;
    xf_Vector *AVmodel_data;
    xf_Vector *GammaVec;
    xf_Data   *GammaDat;
   
    ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
    if(ierr == xf_NOT_FOUND) {xf_printf("Cannot find heat capacity ratio...\n"); return ierr;}
    else
       GammaVec = (xf_Vector *) GammaDat->Data;
  
    if(Model->AVmodel){
       AVmodel_data = Model->AVmodel_data;
       ierr = xf_Error(xf_HaloExchangeVectorBegin(AVmodel_data));
       if (ierr != xf_OK) return ierr;
    }

    // determine if using a Geometric Conservation Law
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
    if (ierr != xf_OK) return ierr;
    
    // Is mesh motion on?
    MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));
    
    // begin communication of halo data: state
    ierr = xf_Error(xf_HaloExchangeVectorBegin(U));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_HaloExchangeVectorBegin(GammaVec));
    if (ierr != xf_OK) return ierr;

    if (UseGCL && MotionOn){
        // find primary GCL vector
        ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
        if (ierr != xf_OK) return ierr;
        
        // begin communication of halo data: GCL
        ierr = xf_Error(xf_HaloExchangeVectorBegin(GCL));
        if (ierr != xf_OK) return ierr;
    }
    
    // initialize R to zero
    ierr = xf_Error(xf_SetZeroVector(R));
    if (ierr != xf_OK) return ierr;
    
    // initialize R_U to zero
    if (R_U != NULL){
        ierr = xf_Error(xf_SetZeroJacobian(All->Mesh, R_U));
        if (ierr != xf_OK) return ierr;
    }
    
    // zero out connectivity if it is required
    if ((SolverData != NULL) && (SolverData->CRequired)){
        ierr = xf_Error(xf_SetZeroVector(SolverData->C));
        if (ierr != xf_OK) return ierr;
        C = SolverData->C;
    }
    
    /* Calculate terms for stabilization (returns immediately if not required) */
    ierr = xf_Error(xf_CalculateStabilization(All, U, R_U != NULL, SolverData));
    if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
   
    /* Calculate residual contribution from boundary faces */
    ierr = xf_Error(xf_CalculateResidualBFaces(All, U, R, R_U, NULL, 0, NULL, SolverData));
    if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
    
    /* Calculate residual contribution from element interiors */
    ierr = xf_Error(xf_CalculateResidualElems(All, U, R, R_U, SolverData));
    if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
   
    ierr = xf_Error(xf_HaloExchangeVectorEnd(GammaVec));
    if (ierr != xf_OK) return ierr;
  
    if(Model->AVmodel)
    {
       ierr = xf_Error(xf_HaloExchangeVectorEnd(AVmodel_data));
       if (ierr != xf_OK) return ierr;
    }

    /* Calculate residual contribution from internal faces (wait for end
     of halo exchange occurs here, unless the wait already occured) */
    ierr = xf_Error(xf_CalculateResidualIFaces(All, U, R, R_U, SolverData));
    if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
    
    // make sure halo is exchanged 
    ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
    if (ierr != xf_OK) return ierr;
    
    if (UseGCL && MotionOn){
        // also make sure GCL is exchanged
        ierr = xf_Error(xf_HaloExchangeVectorEnd(GCL));
        if (ierr != xf_OK) return ierr;
    }
    
    // create lines if required
    if ((SolverData != NULL) && (SolverData->CRequired) && (R_U != NULL)){
        ierr = xf_Error(xf_CreateLines(All, C, R_U, NULL));
        if (ierr != xf_OK) return ierr;
        
        if (SolverData->SortLines){
            ierr = xf_Error(xf_SortLines(All, C, R_U, R));
            if (ierr != xf_OK) return ierr;
        }
    }
    
    return ReturnError;
    
}


/******************************************************************/
//   FUNCTION Definition: xfYu_CalculateResidual
int
xfYu_CalculateResidual(xf_All *All, Yu_Model *pModel, xf_Vector *U, xf_Vector *R, 
		       xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
  int ierr;
  enum xfe_Bool UseGCL;
  enum xfe_Bool MotionOn;
  int ReturnError = xf_OK;
  xf_Vector *C = NULL;
  xf_Vector *GCL = NULL;
  xf_Vector *GammaVec;
  xf_Vector *AVmodel_data;
  xf_Data   *GammaDat;

  //plug in model parameter from Yu
  Model = pModel;
  
  ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
  if(ierr == xf_NOT_FOUND) {xf_printf("Cannot find heat capacity ratio...\n"); return ierr;}
  else
     GammaVec = (xf_Vector *) GammaDat->Data;
   
  // determine if using a Geometric Conservation Law
  //ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  //if (ierr != xf_OK) return ierr;
  
  // Is mesh motion on?
  //MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));

  // begin communication of halo data: state
  ierr = xf_Error(xf_HaloExchangeVectorBegin(U));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_HaloExchangeVectorBegin(GammaVec));
  if (ierr != xf_OK) return ierr;
  if(Model->AVmodel){
     AVmodel_data = Model->AVmodel_data;
     ierr = xf_Error(xf_HaloExchangeVectorBegin(AVmodel_data));
     if (ierr != xf_OK) return ierr;
  }
  /*
  if (UseGCL && MotionOn){
    // find primary GCL vector
    ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
    if (ierr != xf_OK) return ierr;
    
    // begin communication of halo data: GCL
    ierr = xf_Error(xf_HaloExchangeVectorBegin(GCL));
    if (ierr != xf_OK) return ierr;
  }*/

  // initialize R to zero
  ierr = xf_Error(xf_SetZeroVector(R));
  if (ierr != xf_OK) return ierr;
 
  // initialize R_U to zero
  if (R_U != NULL){
    ierr = xf_Error(xf_SetZeroJacobian(All->Mesh, R_U));
    if (ierr != xf_OK) return ierr;
  }
 
  // zero out connectivity if it is required
  /*
  if ((SolverData != NULL) && (SolverData->CRequired)){
    ierr = xf_Error(xf_SetZeroVector(SolverData->C));
    if (ierr != xf_OK) return ierr;
    C = SolverData->C;
  }
  */

  /* Calculate terms for stabilization (returns immediately if not required) */
  //ierr = xf_Error(xf_CalculateStabilization(All, U, R_U != NULL, SolverData));
  //if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;

  /* Calculate residual contribution from boundary faces */
  ierr = xf_Error(xf_CalculateResidualBFaces(All, U, R, R_U, NULL, 0, NULL, SolverData));
  if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;

  /* Calculate residual contribution from element interiors */
  ierr = xf_Error(xf_CalculateResidualElems(All, U, R, R_U, SolverData));
  if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;

  //make sure gamma exchange is done
  ierr = xf_Error(xf_HaloExchangeVectorEnd(GammaVec));
  if (ierr != xf_OK) return ierr;

  if(Model->AVmodel)
  {
     ierr = xf_Error(xf_HaloExchangeVectorEnd(AVmodel_data));
     if (ierr != xf_OK) return ierr;
  }

  /* Calculate residual contribution from internal faces (wait for end
      of halo exchange occurs here, unless the wait already occured) */
  ierr = xf_Error(xf_CalculateResidualIFaces(All, U, R, R_U, SolverData));
  if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
  
  // make sure halo is exchanged 
  ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
  if (ierr != xf_OK) return ierr;

  //if (UseGCL && MotionOn){
  //  // also make sure GCL is exchanged
  //  ierr = xf_Error(xf_HaloExchangeVectorEnd(GCL));
  //  if (ierr != xf_OK) return ierr;
  //} 
  
  // create lines if required
  /*
  if ((SolverData != NULL) && (SolverData->CRequired) && (R_U != NULL)){
    ierr = xf_Error(xf_CreateLines(All, C, R_U, NULL));
    if (ierr != xf_OK) return ierr;

    if (SolverData->SortLines){
      ierr = xf_Error(xf_SortLines(All, C, R_U, R));
      if (ierr != xf_OK) return ierr;
     }
  } */

  return ReturnError;

}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualLeanElem
int
xf_CalculateResidualLeanElem(xf_All *All, int egrp, int elem, xf_Vector *U, 
			     real *ER, real *ER_EU, real **ER_NU, real **NR_EU, 
			     xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
  int ierr, k, sr, nn, r, r2;
  int nface, face, egN, eN;
  int ibfgrp, iiface, ibface;
  int *nvecN;
  enum xfe_Bool IamL;
  real *RL, *RR, *RL_UL, *RR_UR, *RL_UR, *RR_UL;
  xf_IFace IFace;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // pull off info about elem and data
  sr = R_U->StateRank;
  nn = ((R_U->vnvec == NULL) ? R_U->nvec[egrp] : R_U->vnvec[egrp][elem]);
  r  = nn*sr;  r2 = r*r;
  nface = Mesh->ElemGroup[egrp].nFace[elem];

  // nvecN[face] = # basis functions in adjacent element; negative if face is on boundary
  ierr = xf_Error(xf_Alloc( (void **) &nvecN, nface, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (face=0; face<nface; face++){
    egN = R_U->egrpN[egrp][elem][face];
    eN  = R_U->elemN[egrp][elem][face];
    if (egN < 0) nvecN[face] = 0;
    else nvecN[face] = ((R_U->vnvec == NULL) ? R_U->nvec[egN] : R_U->vnvec[egN][eN]);
  }
  
  // zero out residual and Jacobians
  if (ER    != NULL) for (k=0; k<r ; k++) ER[k]    = 0.;
  if (ER_EU != NULL) for (k=0; k<r2; k++) ER_EU[k] = 0.;
  if (ER_NU != NULL)
    for (face=0; face<nface; face++)
      for (k=0; k< r*nvecN[face]*sr; k++) ER_NU[face][k] = 0;
  
  if (NR_EU != NULL)
    for (face=0; face<nface; face++)
      for (k=0; k<nvecN[face]*sr*r; k++) NR_EU[face][k] = 0;
  
  // element interior contribution
  ierr = xf_Error(xf_CalculateResidualElem(All, egrp, elem, U, ER, ER_EU, ER_NU, NULL, SolverData));
  if (ierr != xf_OK) return ierr;


  // interior and boundary face contributions
  for (face=0; face<nface; face++){

    if ((ibfgrp = Mesh->ElemGroup[egrp].Face[elem][face].Group) == xf_INTERIORFACE){

      /** Interior Face **/

      iiface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
      IFace  = Mesh->IFace[iiface];

      // check which side elem is on
      ierr = xf_Error(xf_IsElemOnLeft(IFace, egrp, elem, &IamL));
      if (ierr != xf_OK) return ierr;

      // set RL, RR, etc. appropriately
      if (IamL){
	RL = ER;
	RR = NULL;
	RL_UL = ER_EU;
	RL_UR = ((ER_NU == NULL) ? NULL : ER_NU[face]);
	RR_UR = NULL;
	RR_UL = ((NR_EU == NULL) ? NULL : NR_EU[face]);
      }
      else{
	RR = ER;
	RL = NULL;
	RR_UR = ER_EU;
	RR_UL = ((ER_NU == NULL) ? NULL : ER_NU[face]);
	RL_UL = NULL;
	RL_UR = ((NR_EU == NULL) ? NULL : NR_EU[face]);
      }
       
      // call face residual calculation
      ierr = xf_Error(xf_CalculateResidualIFace(All, iiface, U, RL, RR, RL_UL, RR_UR,
						RL_UR, RR_UL, NULL, SolverData));
      if (ierr != xf_OK) return ierr;
    }
    else if (ibfgrp >= 0){

      /** Boundary Face **/

      ibface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
      
      // skip zero-measure boundary groups
      if (strncmp(Mesh->BFaceGroup[ibfgrp].Title, "ZeroMeasure", 11) != 0){
       
	// call face residual calculation
	ierr = xf_Error(xf_CalculateResidualBFace(All, ibfgrp, ibface, U, ER, ER_EU, NULL, 
						  NULL, SolverData));
	if (ierr != xf_OK) return ierr;
      }
    }
    // do nothing for NULL Faces

  } // face
  

  xf_Release( (void *) nvecN);     

  return xf_OK;

}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualLeanElemFast
static int
xf_CalculateResidualLeanElemFast(xf_All *All, int egrp, int elem, xf_Vector *U, 
                                 real *ER, real *ER_EU, real **ER_NU, real **NR_EU, 
                                 xf_JacobianMatrix *R_U, xf_SolverData *SolverData,
                                 xf_StaticDataElem **pSDE, xf_StaticDataIFace **pSDI,
                                 xf_StaticDataBFace **pSDB)
{
  /* This is a version of CalculateResidualLeanElem that re-uses
     static data if such data is passed in. */
  int ierr, k, sr, nn, r, r2;
  int nface, face, egN, eN;
  int ibfgrp, iiface, ibface;
  int *nvecN;
  enum xfe_Bool IamL;
  real *RL, *RR, *RL_UL, *RR_UR, *RL_UR, *RR_UL;
  xf_IFace IFace;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // pull off info about elem and data
  sr = R_U->StateRank;
  nn = ((R_U->vnvec == NULL) ? R_U->nvec[egrp] : R_U->vnvec[egrp][elem]);
  r  = nn*sr;  r2 = r*r;
  nface = Mesh->ElemGroup[egrp].nFace[elem];

  // nvecN[face] = # basis functions in adjacent element; negative if face is on boundary
  ierr = xf_Error(xf_Alloc( (void **) &nvecN, nface, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (face=0; face<nface; face++){
    egN = R_U->egrpN[egrp][elem][face];
    eN  = R_U->elemN[egrp][elem][face];
    if (egN < 0) nvecN[face] = 0;
    else nvecN[face] = ((R_U->vnvec == NULL) ? R_U->nvec[egN] : R_U->vnvec[egN][eN]);
  }
  
  // zero out residual and Jacobians
  if (ER    != NULL) for (k=0; k<r ; k++) ER[k]    = 0.;
  if (ER_EU != NULL) for (k=0; k<r2; k++) ER_EU[k] = 0.;
  if (ER_NU != NULL)
    for (face=0; face<nface; face++)
      for (k=0; k< r*nvecN[face]*sr; k++) ER_NU[face][k] = 0;
  
  if (NR_EU != NULL)
    for (face=0; face<nface; face++)
      for (k=0; k<nvecN[face]*sr*r; k++) NR_EU[face][k] = 0;
  
  // element interior contribution
  ierr = xf_Error(xf_CalculateResidualElem(All, egrp, elem, U, ER, ER_EU, ER_NU, pSDE, SolverData));
  if (ierr != xf_OK) return ierr;


  // interior and boundary face contributions
  for (face=0; face<nface; face++){

    if ((ibfgrp = Mesh->ElemGroup[egrp].Face[elem][face].Group) == xf_INTERIORFACE){

      /** Interior Face **/

      iiface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
      IFace  = Mesh->IFace[iiface];

      // check which side elem is on
      ierr = xf_Error(xf_IsElemOnLeft(IFace, egrp, elem, &IamL));
      if (ierr != xf_OK) return ierr;

      // set RL, RR, etc. appropriately
      if (IamL){
	RL = ER;
	RR = NULL;
	RL_UL = ER_EU;
	RL_UR = ((ER_NU == NULL) ? NULL : ER_NU[face]);
	RR_UR = NULL;
	RR_UL = ((NR_EU == NULL) ? NULL : NR_EU[face]);
      }
      else{
	RR = ER;
	RL = NULL;
	RR_UR = ER_EU;
	RR_UL = ((ER_NU == NULL) ? NULL : ER_NU[face]);
	RL_UL = NULL;
	RL_UR = ((NR_EU == NULL) ? NULL : NR_EU[face]);
      }
       
      // call face residual calculation
      ierr = xf_Error(xf_CalculateResidualIFace(All, iiface, U, RL, RR, RL_UL, RR_UR,
						RL_UR, RR_UL, pSDI, SolverData));
      if (ierr != xf_OK) return ierr;
    }
    else if (ibfgrp >= 0){

      /** Boundary Face **/

      ibface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
      
      // skip zero-measure boundary groups
      if (strncmp(Mesh->BFaceGroup[ibfgrp].Title, "ZeroMeasure", 11) != 0){
       
	// call face residual calculation
	ierr = xf_Error(xf_CalculateResidualBFace(All, ibfgrp, ibface, U, ER, ER_EU, NULL, 
						  pSDB, SolverData));
	if (ierr != xf_OK) return ierr;
      }
    }
    // do nothing for NULL Faces

  } // face
  

  xf_Release( (void *) nvecN);     

  return xf_OK;

}


/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionGCL_PsiTxR_G
int
xf_MeshMotionGCL_PsiTxR_G(xf_All *All, xf_Vector *U, int nPsi, 
                          xf_Vector **Psi, xf_Vector **SGCL)
{
  // Calculates Psi^T * (Residual_gbar) using finite difference approximations
  int ierr, k, sr;
  int egrp, elem, face, nface;
  int rmax, nn, r;
  int nN, nnN;
  int egN, eN;
  int iAdjoint;
  enum xfe_Bool Found;
  real ep = 1e-5; // epsilon for finite-differencing gbar
  real *ER0, *ER;
  real *EGCLN, *EPsi, *ESGCLN;
  real dp;
  xf_JacobianMatrix *R_U;
  xf_SolverData *SolverData;
  xf_StaticDataElem *SDE = NULL;
  xf_StaticDataIFace *SDI = NULL;
  xf_StaticDataBFace *SDB = NULL;
  xf_Vector *GCL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // locate current GCL state
  ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
  if (ierr != xf_OK) return ierr;

  // begin communication: make sure Halos for U and GCL are up-to-date
  ierr = xf_Error(xf_HaloExchangeVectorBegin(U));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_HaloExchangeVectorBegin(GCL));
  if (ierr != xf_OK) return ierr;

  // locate Jacobian matrix
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_False, NULL, &R_U, &Found));
  if (ierr != xf_OK) return ierr;
  if (!Found) return xf_Error(xf_NOT_FOUND);

  // state rank
  sr = R_U->StateRank;
  
  // get rmax = max(nn*sr)
  rmax=0;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      rmax = max(  sr*xf_Jacobian_n(R_U,egrp,elem), rmax);

  // allocate ER and ER0 = residual vectors on a single element
  ierr = xf_Error(xf_Alloc( (void **) &ER0, 2*rmax, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ER = ER0 + rmax;

  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;
  
  // end halo communication
  ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_HaloExchangeVectorEnd(GCL));
  if (ierr != xf_OK) return ierr;

  // zero out Halo entries in SGCL
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    ierr = xf_Error(xf_SetZeroHaloVector(SGCL[iAdjoint]));
    if (ierr != xf_OK) return ierr;
  }

  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      nn = ((R_U->vnvec == NULL) ? R_U->nvec[egrp] : R_U->vnvec[egrp][elem]);
      r  = nn*sr;
      nface = Mesh->ElemGroup[egrp].nFace[elem];

      // calculate ER0 = spatial residual on elem
      ierr = xf_Error(xf_CalculateResidualLeanElemFast(All, egrp, elem, U, ER0, NULL, 
                                                       NULL, NULL, R_U, SolverData,
                                                       &SDE, &SDI, &SDB));
      if (ierr != xf_OK) return ierr;

      // loop over self and neighbors
      for (face=-1; face<nface; face++){
	if (face==-1){ // self
	  egN = egrp;
	  eN  = elem;
	}
	else{          // neighbor
	  egN = R_U->egrpN[egrp][elem][face];
	  eN  = R_U->elemN[egrp][elem][face];
	}
	
	if (egN < 0) continue; // boundary face
	
	// # of basis functions in adjacent element
	nnN = ((R_U->vnvec == NULL) ? R_U->nvec[egN] : R_U->vnvec[egN][eN]);

	// GCL on neighbor
	EGCLN = GCL->GenArray[egN].rValue[eN];
			
	// loop over basis functions on adjacent element
	for (nN = 0; nN < nnN; nN++){
	  
	  // increment GCL on adjacent element
	  EGCLN[nN] += ep;

	  // calculate ER = spatial residual on elem
	  ierr = xf_Error(xf_CalculateResidualLeanElemFast(All, egrp, elem, U, ER, NULL, 
                                                           NULL, NULL, R_U, SolverData,
                                                           &SDE, &SDI, &SDB));
	  if (ierr != xf_OK) return ierr;
	  
	  // Set ER{k} = (ER{k}-ER0{k})/ep
	  for (k=0; k<r; k++) ER[k] = (ER[k]-ER0[k])/ep;

	  // Add ER{k} dot Psi[iAdjoint]{k} to ESGCLN{nN}
	  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){

	    // Adjoint, Psi, on current element
	    EPsi = Psi[iAdjoint]->GenArray[egrp].rValue[elem];

	    // SGCL on neighbor element
	    ESGCLN = SGCL[iAdjoint]->GenArray[egN].rValue[eN];

	    // take dot product
	    for (k=0, dp=0.; k<r; k++) dp += EPsi[k]*ER[k];

	    // increment SGCL on neighbor
	    ESGCLN[nN] += dp;

	  } // iAdjoint

	  // decrement GCL on adjacent element (return to orig value)
	  EGCLN[nN] -= ep;
	  
	} // nN

      } // face

    } // elem

  } // egrp

  // Reverse-halo-exchange SGCL to include SGCL from halos
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    ierr = xf_Error(xf_HaloReverseExchangeVectorBegin(SGCL[iAdjoint]));
    if (ierr != xf_OK) return ierr;
  }
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    ierr = xf_Error(xf_HaloReverseExchangeVectorEnd(SGCL[iAdjoint]));
    if (ierr != xf_OK) return ierr;
  }


  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;

  // destroy static data
  ierr = xf_Error(xf_DestroyStaticDataElem(SDE));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyStaticDataIFace(SDI));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyStaticDataBFace(SDB));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) ER0);

  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_CalculateAdjointResidual
int
xf_CalculateAdjointResidual(xf_All *All, xf_Vector *U, xf_Vector *Adj, xf_Vector *R, 
			    xf_JacobianMatrix *R_U, xf_SolverData *SolverDataIn)
{
  int ierr;
  int ReturnError = xf_OK;
  xf_SolverData *SolverData;
  xf_Vector *C = NULL;

  // begin communication of halo data: state
  ierr = xf_Error(xf_HaloExchangeVectorBegin(U));
  if (ierr != xf_OK) return ierr;

  // initialize R to zero
  ierr = xf_Error(xf_SetZeroVector(R));
  if (ierr != xf_OK) return ierr;
 
  // initialize R_U to zero
  if (R_U != NULL){
    ierr = xf_Error(xf_SetZeroJacobian(All->Mesh, R_U));
    if (ierr != xf_OK) return ierr;
  }
  
  // zero out connectivity if it is required
  if ((SolverDataIn != NULL) && (SolverDataIn->CRequired)){
    ierr = xf_Error(xf_SetZeroVector(SolverDataIn->C));
    if (ierr != xf_OK) return ierr;
    C = SolverDataIn->C;
  }

  // create a dummy solver data if passed in NULL (e.g. for stabilization)
  if (SolverDataIn == NULL){
    ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
    if (ierr != xf_OK) return ierr;
  }
  else SolverData = SolverDataIn;

  /* Calculate terms for stabilization (returns immediately if not required) */
  ierr = xf_Error(xf_CalculateStabilization(All, U, R_U != NULL, SolverData));
  if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;

  /* Calculate residual contribution from boundary faces */
  ierr = xf_Error(xf_CalculateResidualBFaces(All, U, R, R_U, NULL, 0, NULL, SolverData));
  if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;

  /* Calculate residual contribution from element interiors */
  ierr = xf_Error(xf_CalculateResidualElems(All, U, R, R_U, SolverData));
  if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;

  /* Calculate residual contribution from internal faces (wait for end
     of halo exchange occurs here, unless the wait already occured) */
  ierr = xf_Error(xf_CalculateResidualIFaces(All, U, R, R_U, SolverData));
  if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;

  // create lines if required
  if ((SolverData != NULL) && (SolverData->CRequired) && (R_U != NULL)){
    ierr = xf_Error(xf_CreateLines(All, C, R_U, NULL));
    if (ierr != xf_OK) return ierr;

    if (SolverData->SortLines){
      ierr = xf_Error(xf_SortLines(All, C, R_U, R));
      if (ierr != xf_OK) return ierr;
    }
  }

  // destroy SolverData if created it
  if (SolverDataIn == NULL){
    ierr = xf_Error(xf_DestroySolverData(SolverData));
    if (ierr != xf_OK) return ierr;
  }

  return ReturnError;

}




#if( UNIT_TEST==1 )
#include "xf_Residual.test.in"
#endif
