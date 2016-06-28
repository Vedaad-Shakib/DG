/* this is an addon to the xflow for combustion simulation
 * Author: YU LV
 * Date: Sep 2012
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
#include "xfYu_Solver.h"
#include "xfYu_Model.h"
#include "xf_MeshMotion.h"

/* static pointer to global model structure*/
static Yu_Model *Model;

/* for element interior residual */
typedef struct
{
    int             iConv[2], iDiff[2], iSource[2];
    xf_QuadData     *QuadData;
    xf_BasisData    *PhiData;
    xf_JacobianData *JData;
    xf_BasisData    *GeomPhiData;
    real            *xglob;
    real            *T, *wq, *u, *gu;
    real            *F;     //convective term
    real            *S;     //source term
    real            *ResMetric;
    int             pnq, nqConv, nqDiff, Tsize;
    xf_Vector       *EG;
    xf_Vector       *EM;
    real            Time;
    real            *gum; 
}
xf_StaticDataElem;

typedef struct
{
    int             iConv[2], iDiff[2], iSource[2];
    xf_QuadData     *QuadData;
    xf_BasisData    *PhiDataL, *PhiDataR;
    xf_BasisTable   *PhiTable;
    xf_NormalData   *NData;
    xf_JacobianData *JDataL, *JDataR;
    xf_BasisData    *GeomPhiData;
    real            *xglob;
    real            *T, *wn;
    real            *xelemL, *xelemR;
    real            *uL, *uR;
    real            *guL, *guR;
    real            *F;    //convective flux
    int             pnq, nqConv, nqDiff, nnDiff, Tsize;
    int             *IParam;
    real            *RParam;
    xf_Vector       *EG;
    xf_Vector       *EM;
    real            *ResMetricL, *ResMetricR;
    real            Time;
    xf_MotionData   *MD;
}
xf_StaticDataIFace;

/* Structure for storing "static" data in bface residual calculation. */
typedef struct
{
    xf_QuadData     *QuadData;
    xf_BasisData    *PhiData;
    xf_BasisTable   *PhiTable;
    xf_NormalData   *NData;
    xf_JacobianData *JData;
    xf_BasisData    *GeomPhiData;
    real            *xglob;
    real            *T, *T2, *wn;
    real            *xelem;
    real            *uI;
    real            *guI;
    real            *uB;
    real            *F;
    int             pnq, Tsize, T2size;
    xf_Vector       *EG;
    xf_Vector       *EM;
    real            Time;
}
xf_StaticDataBFace;

/******************************************************************/
//   FUNCTION Definition: xf_CreateStaticDataElem
static int
xf_CreateStaticDataElem(xf_All *All, xf_SolverData *SolverData, xf_StaticDataElem **pSD)
{
    
    int ierr;
    xf_StaticDataElem *SD;

    // allocate memory
    ierr = xf_Error(xf_Alloc( (void **) pSD, 1, sizeof(xf_StaticDataElem)));
    if (ierr != xf_OK) return ierr;
    SD = (*pSD);
    
    // initialize variables to NULL
    SD->QuadData    =  NULL;
    SD->PhiData     =  NULL;
    SD->JData       =  NULL;
    SD->wq	        =  NULL;
    SD->GeomPhiData =  NULL;
    SD->xglob       =  NULL;
    SD->T	        =  NULL;
    SD->u	        =  NULL; 
    SD->gu	        =  NULL; 
    SD->F	        =  NULL; 
    SD->S	        =  NULL;
    SD->ResMetric   =  NULL;
    SD->pnq         =  -1; 
    SD->nqConv      =  -1; 
    SD->nqDiff      =  -1; 
    SD->Tsize       =  -1;
    SD->EG          =  NULL;
    SD->EM          =  NULL;
    SD->gum         =  NULL;
    
    // determine Time
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &SD->Time));
    if (ierr != xf_OK) return ierr;
   
    // obtain elem geometry vector
    //ierr = xf_Error(xf_FindElemGeom(All, &SD->EG));
    //if (ierr != xf_OK) return ierr;

    // everything must go back with (*pSD)
    (*pSD) = SD;
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyStaticDataElem
static int
xf_DestroyStaticDataElem(xf_StaticDataElem *SD)
{
    int ierr;

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

    // Release all other memory
    xf_Release( (void *) SD->wq);
    xf_Release( (void *) SD->xglob);
    xf_Release( (void *) SD->T);
    xf_Release( (void *) SD->u);
    xf_Release( (void *) SD->gu);
    xf_Release( (void *) SD->F);
    xf_Release( (void *) SD->ResMetric);
    xf_Release( (void *) SD->S);
    xf_Release( (void *) SD->gum);
    
    // Destroy self
    xf_Release( (void *) SD);
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualElem
static int
xf_CalculateResidualElem(xf_All *All, int egrp, int elem, xf_Vector *U, real *ER, 
                         xf_StaticDataElem **pSD, xf_SolverData *SolverData)
{
    int ierr, sr, sr2, iq, nq, nn, d, n, i, j;
    int dim, Order, QuadOrder;
    int Tsize, pnq, nConv,  nDiff,  nSource;
    int *iConv, *iDiff, *iSource;
    enum xfe_BasisType Basis;
    enum xfe_Bool QuadChanged, ConstA;
    real *EU, *F, *S, Gamma;
    real *xq, *wq, *xglob, *u, *gu, *T, *gum = NULL, *w = NULL;
    real Time, h;
    xf_QuadData *QuadData;
    xf_BasisData *PhiData, *GeomPhiData;
    xf_JacobianData *JData;
    xf_Vector *EG, *V, *GammaVec;
    xf_StaticDataElem *SD = NULL;
    xf_Mesh   *Mesh;
    xf_Data   *GammaDat;

    //pull in general information
    Mesh    = All->Mesh;
    dim     = Mesh->Dim; 
    sr      = Model->nVars;
    sr2     = sr * sr;
    
    // Determine Basis and Order from the state, U
    Basis = U->Basis[egrp];
    Order = xf_InterpOrder(U, egrp, elem);
    
    // Create and initialize static data if not passed in
    if ((pSD == NULL) || ((*pSD) == NULL)){
        ierr = xf_Error(xf_CreateStaticDataElem(All, SolverData, (pSD != NULL) ? pSD : &SD));
        if (ierr != xf_OK) return ierr;
    }
    if (pSD != NULL) SD = (*pSD);
   
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
    T           = SD->T;
    u           = SD->u;
    gu          = SD->gu;
    F           = SD->F;
    S           = SD->S;
    pnq         = SD->pnq;        // previous # quad points
    Tsize       = SD->Tsize;      // for temporary matrix allocation
    EG          = SD->EG;         // element geometry vector
    Time        = SD->Time;       // simulation time
    
    // determine required integration order
    ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order, &QuadOrder));
    if (ierr != xf_OK) return ierr;
    
    /* Pull off quad points for the element; will not recalculate if
     Basis/Order have not changed. */
    ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
    if (ierr != xf_OK) return ierr;
   
    // check element size
    //ierr = xf_Error(xf_ElemSize(All, egrp, elem, EG, &h)); // A/perim
    //if (ierr != xf_OK) return ierr;

    //printf("%lf\n", h);

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
        ierr = xf_Error(xf_ReAlloc( (void **) &gu, nq*sr*dim, sizeof(real)));
        if (ierr != xf_OK) return ierr;
    }

    // obtain global coords of quad points
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
                                    nq, xq, xglob));
    if (ierr != xf_OK) return ierr;
    
    EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
    
    // interpolate state and gradient at quad points
    xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u); 
    
    // form detJ-multiplied quad weight vector, wq
    for (iq=0; iq<nq; iq++) 
        wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
    
    /*------------------*/
    /* CONVECTION TERMS */
    /*------------------*/
    
    // allocate memory
    ierr = xf_Error(xf_ReAlloc( (void **) &F, sr*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // calculate F [and F_u] at quad points  
    // modify later........................
    //ierr = xf_Error(xf_EqnSetConvF(EqnSet, ResTerm+iConv[0], nConv, IParam, 
    //                               RParam, nq, u, Auxu, xglob, F, F_u));
    //if (ierr != xf_OK) return ierr;
    //Yu's add-on
    ierr = xf_Error(ConvFluxInterior(nq, sr, dim, u, F, Gamma));
    if(ierr != xf_OK) return ierr;


    // multiply F by quad weights*J
    for (d=0; d<dim; d++)
        xf_ColMult(F+nq*sr*d, wq, nq, sr, 1); // F is modified here
	
    /* Add to R:
     ER{n,k} -= sum_i sum_q gPhi{i,q,n}^T * F{i,q,k}*wq{q}
     */
    for (d=0; d<dim; d++)
    { 
       xf_MTxM_Sub(PhiData->gPhi+nn*nq*d, F+nq*sr*d, nn, nq, sr, ER);
    }
    
    /*--------------*/
    /* SOURCE TERMS */
    /*--------------*/
   
    /*
    // allocate memory
    ierr = xf_Error(xf_ReAlloc( (void **) &S, sr*nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // calculate S [and S_u] at quad points
    ierr = xf_Error(xf_EqnSetSourceS(EqnSet, ResTerm+iSource[0], nSource, IParam, RParam, 
                                     nq, u, gu, Auxu, xglob, &Time, S, S_u, S_gu, &Nonzero_gu));
    if (ierr != xf_OK) return ierr;
    
    // multiply S by quad weights*J
    xf_ColMult(S, wq, nq, sr, 1); // S is modified here
    
    // Add to R:
    // ER{n,k} += sum_q Phi{q,n}^T * S{q,k}*wq{q}
     
    xf_MTxM_Add(PhiData->Phi, S, nn, nq, sr, ER);
    */

    // end of residual evalation
    pnq = nq; // set previous quad point # for next element
    
    // Store possibly-altered or resized data back in StaticData
    SD->QuadData    =  QuadData;
    SD->PhiData     =  PhiData;  
    SD->JData       =  JData;
    SD->wq          =  wq;
    SD->GeomPhiData =  GeomPhiData;
    SD->xglob       =  xglob;
    SD->T	        =  T;
    SD->u	        =  u;
    SD->gu	        =  gu;
    SD->F	        =  F;
    SD->S	        =  S;
    SD->pnq         =  pnq;
    
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
    int ierr, egrp, elem;
    real *ER;
    xf_StaticDataElem *StaticData = NULL;
    xf_Mesh   *Mesh;
    
    Mesh = All->Mesh;
    
    // loop over element groups
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
        
        // loop over elements
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
            
            /* Prepare pointers to residual and Jacobian */
            ER = R->GenArray[egrp].rValue[elem];
            
            /* Calculate residual on elem, passing in StaticData  */
            ierr = xf_Error(xf_CalculateResidualElem(All, egrp, elem, U, ER, 
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
    int ierr;
    xf_StaticDataIFace *SD;
    
    // allocate memory
    ierr = xf_Error(xf_Alloc( (void **) pSD, 1, sizeof(xf_StaticDataIFace)));
    if (ierr != xf_OK) return ierr;
    SD = (*pSD);
    
    // initialize variables to NULL
    SD->QuadData    =  NULL;
    SD->PhiDataL    =  NULL;   SD->PhiDataR    =  NULL;
    SD->NData       =  NULL;
    SD->JDataL      =  NULL;   SD->JDataR      =  NULL;
    SD->GeomPhiData =  NULL;
    SD->wn	        =  NULL;
    SD->xglob       =  NULL;
    SD->xelemL      =  NULL;   SD->xelemR      =  NULL;
    SD->T	        =  NULL;
    SD->uL          =  NULL;   SD->uR          =  NULL; 
    SD->guL	        =  NULL;   SD->guR         =  NULL; 
    SD->F	        =  NULL; 
    SD->pnq         =  -1;
    SD->EG          =  NULL;
    SD->EM          =  NULL;
    
    // determine Time
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &SD->Time));
    if (ierr != xf_OK) return ierr;
    
    /* Create a basis table, PhiTable, that will store computed basis
     functions specific to each [element shape, face in element,
     orientation of face] combination, for quick lookup. */
    ierr = xf_Error(xf_CreateBasisTable(&SD->PhiTable));
    if (ierr != xf_OK) return ierr;
    
    // find element geometry vector
    ierr = xf_Error(xf_FindElemGeom(All, &SD->EG));
    if (ierr != xf_OK) return ierr;
    
    // everything must go back with (*pSD)
    (*pSD) = SD;
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyStaticDataIFace
static int
xf_DestroyStaticDataIFace(xf_StaticDataIFace *SD)
{
    int ierr;
    
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
    ierr = xf_Error(xf_DestroyBasisData(SD->GeomPhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Normal Data */
    ierr = xf_Error(xf_DestroyNormalData(SD->NData));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy geometry Jacobian Data */
    ierr = xf_Error(xf_DestroyJacobianData(SD->JDataL));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DestroyJacobianData(SD->JDataR));
    if (ierr != xf_OK) return ierr;
    
    //release memory
    xf_Release( (void *) SD->wn);
    xf_Release( (void *) SD->xglob);
    xf_Release( (void *) SD->xelemL);
    xf_Release( (void *) SD->xelemR);
    xf_Release( (void *) SD->uL);
    xf_Release( (void *) SD->uR);
    xf_Release( (void *) SD->guL);
    xf_Release( (void *) SD->guR);
    xf_Release( (void *) SD->F);
    xf_Release( (void *) SD->T);
    
    // Destroy self
    xf_Release( (void *) SD);
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualIFace
static int
xf_CalculateResidualIFace(xf_All *All, int iiface, xf_Vector *U, real *RL,
                          real *RR, xf_StaticDataIFace **pSD, xf_SolverData *SolverData)
{
    int ierr, sr, sr2, iq, nq, nn, nL, nR, d, n, i, j;
    int dim, OrderL, OrderR, QuadOrder;
    int egrpL, elemL, faceL;
    int egrpR, elemR, faceR;
    int pnq;
    enum xfe_BasisType BasisL, BasisR;
    enum xfe_Bool QuadChanged;
    real *xq, *wq, *wn, *xglob, *xelemL, *xelemR;
    real *UL, *UR, *uL, *uR, *guL, *guR, *F;
    real Time, ElemVolL, ElemVolR;
    xf_QuadData *QuadData;
    xf_BasisTable *PhiTable;
    xf_BasisData *PhiDataL, *PhiDataR, *GeomPhiData;
    xf_NormalData *NData;
    xf_JacobianData *JDataL, *JDataR;
    xf_Vector *EG, *C, *V;
    xf_StaticDataIFace *SD = NULL;
    xf_IFace IFace;
    xf_Mesh  *Mesh;
    
    // General information
    Mesh = All->Mesh;
    dim  = Mesh->Dim;
    sr   = Model->nVars;
    sr2  = sr*sr;
    
    // Interior face structure
    IFace = Mesh->IFace[iiface];
    
    // Create and initialize static data if not passed in
    if ((pSD == NULL) || ((*pSD) == NULL)){
        ierr = xf_Error(xf_CreateStaticDataIFace(All, SolverData, (pSD != NULL) ? pSD : &SD));
        if (ierr != xf_OK) return ierr;
    }
    if (pSD != NULL) SD = (*pSD);
    
    // pull in variables
    QuadData    = SD->QuadData;
    PhiDataL    = SD->PhiDataL;       PhiDataR    = SD->PhiDataR;
    PhiTable    = SD->PhiTable;
    NData       = SD->NData;
    JDataL      = SD->JDataL;         JDataR      = SD->JDataR;
    GeomPhiData = SD->GeomPhiData;
    wn          = SD->wn;
    xglob       = SD->xglob;
    xelemL      = SD->xelemL;         xelemR      = SD->xelemR;
    uL          = SD->uL;             uR          = SD->uR;
    guL         = SD->guL;            guR         = SD->guR;
    F           = SD->F;
    pnq         = SD->pnq; 
    EG          = SD->EG;      
    Time        = SD->Time;       // simulation time
    
    // Set Connectivity vector if required
    //if ((SolverData != NULL) && (SolverData->CRequired))
    //    C = SolverData->C;
    //else
    //    C = NULL;
    
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
    
    // determine required integration order
    ierr = xf_Error(xf_GetQuadOrderIFace(Mesh, NULL, IFace, max(OrderL,OrderR), 
                                         &QuadOrder));
    if (ierr != xf_OK) return ierr;
    
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
    
    /* Compute normal(s) at quad points.  If face is straight, only
     one normal will be computed/returned. */
    ierr = xf_Error(xf_IFaceNormal(Mesh, IFace, nq, xq, &NData));
    if (ierr != xf_OK) return ierr;
    
    nL = PhiDataL->nn;
    nR = PhiDataR->nn;
    nn = max(nL, nR);
    
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
        
        //if (C != NULL){
        //    ierr = xf_Error(xf_ReAlloc( (void **) &CF, nq, sizeof(real)));
        //    if (ierr != xf_OK) return ierr;
        //}
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
    
    // interpolate state and gradient at quad points
    xf_MxM_Set(PhiDataL->Phi, UL, nq, nL, sr, uL);
    xf_MxM_Set(PhiDataR->Phi, UR, nq, nR, sr, uR);
    
    // zero out Connectivity
    if (C != NULL){
        C->GenArray[egrpL].rValue[elemL][faceL] = 0.;
        C->GenArray[egrpR].rValue[elemR][faceR] = 0.;
    }
    
    /*------------------*/
    /* CONVECTION TERMS */
    /*------------------*/

    ierr = xf_Error(xf_ReAlloc( (void **) &F, sr*nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // calculate F [and F_u] at quad points for the left state(flag ==0)
    // need to modify ......
    //ierr = xf_Error(xf_EqnSetConvFJump(EqnSet, ResTerm+iConv[0], nConv, IParam, RParam, 
    //                                   nq, uL, uR, AuxuL, AuxuR, wn, xglob, 
    //                                   (MotionOn) ? MD->vg : NULL, F, F_uL, F_uR, CF, &flag)); //0 for left state
    //if (ierr != xf_OK) return ierr;
    ierr = xf_Error(ConvFluxInteriorFace(nq, sr, dim, uL, uR, wn, F));
    if (ierr != xf_OK) return ierr;
    
    /* Add to R:
     RL_{n,k} += sum_q PhiL_{n,q}*F_{k,q}
     RR_{n,k} -= sum_q PhiR_{n,q}*F_{k,q}
     */
    if (RL != NULL) xf_MTxM_Add(PhiDataL->Phi, F, nL, nq, sr, RL);
    if (RR != NULL) xf_MTxM_Sub(PhiDataR->Phi, F, nR, nq, sr, RR);
    
    //add to connectivity
    //if (C != NULL){ 
    //    for (iq=0, cval=0.; iq<nq; iq++) cval += CF[iq];
    //    C->GenArray[egrpL].rValue[elemL][faceL] += fabs(cval);
    //    C->GenArray[egrpR].rValue[elemR][faceR] += fabs(cval);
    //}
    
    pnq = nq; // set previous quad point # for next iface
    
    // Store possibly-altered or resized data back in StaticData
    SD->QuadData    = QuadData;
    SD->PhiDataL    = PhiDataL;       SD->PhiDataR    = PhiDataR;
    SD->PhiTable    = PhiTable;	          
    SD->NData       = NData;	          
    SD->JDataL      = JDataL;         SD->JDataR      = JDataR;
    SD->GeomPhiData = GeomPhiData; 
    SD->wn          = wn;	          
    SD->xglob       = xglob;	          
    SD->xelemL      = xelemL;         SD->xelemR      = xelemR;
    SD->uL          = uL;             SD->uR          = uR;
    SD->guL         = guL;            SD->guR         = guR;
    SD->F           = F;	
    SD->pnq         = pnq;  
    SD->EG          = EG;
    
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
    int ierr;
    int nIFaceRegular, iiface;
    int egrpL, egrpR, elemL, elemR, faceL, faceR;
    int ReturnError = xf_OK;
    real *RL = NULL, *RR = NULL;
    enum xfe_Bool MeshIsParallel;
    xf_StaticDataIFace *StaticData = NULL;
    xf_IFace IFace;
    xf_Mesh   *Mesh;
    
    Mesh = All->Mesh;
    
    nIFaceRegular = -1;
    if (MeshIsParallel = (Mesh->ParallelInfo != NULL))
        nIFaceRegular = Mesh->ParallelInfo->nIFaceRegular;
    
    // Loop over interior faces
    for (iiface=0; iiface<Mesh->nIFace; iiface++){
        
        if (MeshIsParallel && (iiface == nIFaceRegular)){
            // wait for end communication of halo data: state
            ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
            if (ierr != xf_OK) return ierr;
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
        
        /* Calculate residual on IFace, passing in StaticData  */
        ierr = xf_Error(xf_CalculateResidualIFace(All, iiface, U, RL, RR, &StaticData, SolverData));
        if (ierr != xf_OK){  // for parallel's sake, do not return immediately if recoverable
            //if (!xf_CheckRecoverable(ierr, &ReturnError)) return ierr;
            //continue; // recoverable error occured, move on
            return ierr;
        }
        
    } // iiface
    
    // Delete StaticData
    ierr = xf_Error(xf_DestroyStaticDataIFace(StaticData));
    if (ierr != xf_OK) return ierr;
    
    return ReturnError;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateStaticDataBFace
static int
xf_CreateStaticDataBFace(xf_All *All, xf_SolverData *SolverData, xf_StaticDataBFace **pSD)
{
    int ierr;
    xf_StaticDataBFace *SD;
    
    // allocate memory
    ierr = xf_Error(xf_Alloc( (void **) pSD, 1, sizeof(xf_StaticDataBFace)));
    if (ierr != xf_OK) return ierr;
    SD = (*pSD);
    
    // initialize variables to NULL
    SD->QuadData    =  NULL;
    SD->PhiData     =  NULL; 
    SD->NData       =  NULL;
    SD->JData       =  NULL; 
    SD->GeomPhiData =  NULL;
    SD->wn	      =  NULL;
    SD->xglob       =  NULL;
    SD->xelem       =  NULL; 
    SD->T	          =  NULL;
    SD->T2          =  NULL;
    SD->uI          =  NULL; 
    SD->guI	      =  NULL; 
    SD->uB          =  NULL; 
    SD->F	          =  NULL; 
    SD->pnq         =  -1; 
    SD->Tsize       =  -1;
    SD->T2size      =  -1;
    SD->EG          =  NULL;
    SD->EM          =  NULL;
    
    // determine Time
    ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &SD->Time));
    if (ierr != xf_OK) return ierr; 
    
    /* Create a basis table, PhiTable, that will store computed basis
     functions specific to each [element shape, face in element,
     orientation of face] combination, for quick lookup. */
    ierr = xf_Error(xf_CreateBasisTable(&SD->PhiTable));
    if (ierr != xf_OK) return ierr;
    
    // find element geometry vector
    ierr = xf_Error(xf_FindElemGeom(All, &SD->EG));
    if (ierr != xf_OK) return ierr;
    
    // everything must go back with (*pSD)
    (*pSD) = SD;
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyStaticDataBFace
static int
xf_DestroyStaticDataBFace(xf_StaticDataBFace *SD)
{
    int ierr;
    
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
    ierr = xf_Error(xf_DestroyBasisData(SD->GeomPhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy Normal Data */
    ierr = xf_Error(xf_DestroyNormalData(SD->NData));
    if (ierr != xf_OK) return ierr;
    
    /* Destroy geometry Jacobian Data */
    ierr = xf_Error(xf_DestroyJacobianData(SD->JData));
    if (ierr != xf_OK) return ierr;
    
    // Release memory
    xf_Release( (void *) SD->T);
    xf_Release( (void *) SD->T2);
    xf_Release( (void *) SD->wn);
    xf_Release( (void *) SD->xglob);
    xf_Release( (void *) SD->xelem);
    xf_Release( (void *) SD->uI);
    xf_Release( (void *) SD->guI);
    xf_Release( (void *) SD->uB);
    xf_Release( (void *) SD->F);
    
    // Destroy self
    xf_Release( (void *) SD);
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculateResidualBFace
static int
xf_CalculateResidualBFace(xf_All *All, int ibfgrp, int ibface, const xf_Vector *U, 
                          real *ER, xf_OutputEvalData *OutputEval, 
                          xf_StaticDataBFace **pSD, xf_SolverData *SolverData)
{
    int ierr, sr, sr2, i, j, k, l, iq, nq, nn, d, n;
    int dim, Order, QuadOrder;
    int egrp, elem, face;
    int pnq, Tsize, T2size;
    int localBCflag;
    enum xfe_BasisType Basis;
    enum xfe_DiffDiscType DiffDisc;
    enum xfe_Bool QuadChanged;
    char Title[xf_MAXSTRLEN];
    real *T, *T2, cval, ElemVol, Time;
    real *xq, *wq, *wn, *xglob, *xelem;
    real *EU, *uI, *guI, *uB, *F;
    xf_BFace BFace;
    xf_QuadData *QuadData;
    xf_BasisTable *PhiTable;
    xf_BasisData *PhiData, *GeomPhiData;
    xf_NormalData *NData;
    xf_JacobianData *JData;
    xf_StaticDataBFace *SD = NULL;
    xf_Vector *EG, *C, *V;
    xf_Mesh   *Mesh;
    
    // General information
    Mesh = All->Mesh;
    dim  = Mesh->Dim;
    sr  = Model->nVars;
    sr2 = sr*sr;
    
    BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
    
    // Create and initialize static data if not passed in
    if ((pSD == NULL) || ((*pSD) == NULL)){
        ierr = xf_Error(xf_CreateStaticDataBFace(All, SolverData, (pSD != NULL) ? pSD : &SD));
        if (ierr != xf_OK) return ierr;
    }
    if (pSD != NULL) SD = (*pSD);
    
    //pull in variables from StaticData
    QuadData    = SD->QuadData;
    PhiData     = SD->PhiData;
    PhiTable    = SD->PhiTable;
    NData       = SD->NData;
    JData       = SD->JData;
    GeomPhiData = SD->GeomPhiData;
    T           = SD->T;
    T2          = SD->T2;
    wn          = SD->wn;
    xglob       = SD->xglob;
    xelem       = SD->xelem;
    uI          = SD->uI;
    guI         = SD->guI;
    uB          = SD->uB;
    F           = SD->F;
    pnq         = SD->pnq;
    Tsize       = SD->Tsize;      
    T2size      = SD->T2size;
    EG          = SD->EG;
    Time        = SD->Time;       // simulation time
    
    // Set Connectivity vector if required
    //if ((SolverData != NULL) && (SolverData->CRequired))
    //    C = SolverData->C;
    //else
    //    C = NULL;
    
    // Sort EqnSet->BCs to match boundary face groups
    //ierr = xf_Error(xf_SortBCs(Mesh, EqnSet->BCs+0));
    //if (ierr != xf_OK) return ierr;
    //BC = EqnSet->BCs[0].BC;
    //match up Yu's new BC structure in Yu_Model
    for(i=0; i<Model->nBCs; i++)
       if(strcmp(Mesh->BFaceGroup[ibfgrp].Title, Model->nameBCs[i]) == 0)
       {  localBCflag = i;  break; }

    if(i==Model->nBCs) {
       printf("Boundary name does not match model defining.\n");
       xf_Error(xf_BOUNDARY_CONDITION_ERROR);
    }
    
    egrp = BFace.ElemGroup;
    elem = BFace.Elem;
    face = BFace.Face;
    
    Basis = U->Basis[egrp];
    Order = xf_InterpOrder(U, egrp, elem);
    
    // determine required integration order
    ierr = xf_Error(xf_GetQuadOrderBFace(Mesh, NULL, BFace, Order, &QuadOrder));
    if (ierr != xf_OK) return ierr;
   
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
    
    /* Compute normal(s) at quad points.  If face is straight, only
     one normal will be computed/returned. */
    ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, nq, xq, &NData, NULL));
    if (ierr != xf_OK) return ierr;
    
    nn = PhiData->nn;
   
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
        //if (C != NULL){
        //    ierr = xf_Error(xf_ReAlloc( (void **) &CF, nq, sizeof(real)));
        //    if (ierr != xf_OK) return ierr;
        //}
    }
    
    // pull off element volume
    ElemVol = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
    
    // obtain global coords of quad points
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, xfe_True, 
                                    nq, xelem, xglob));
    if (ierr != xf_OK) return ierr;
    
    // interpolate state and gradient at quad points
    xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, uI);
    
    /* Initially, set wn = just the normals.  Will include quad
     weights after call to EqnSetBCState, since some weights may
     be negative while EqnSetBCState needs to have the true
     outward-pointing normal at each point.
     */
    for (d=0; d<dim; d++)
        for (iq=0;iq<nq; iq++) 
            wn[iq*dim+d] = NData->n[iq*dim*(NData->nq!=1)+d];
    
    // obtain boundary state and derivative at quad points
    // need to modify later......
    // deine boundary state ......
    //ierr = xf_Error(xf_EqnSetBCState(EqnSet, BC+ibfgrp, IParam, RParam, nq, 
    //                                 wn, xglob, &Time, (MotionOn) ? MD->vg : NULL, 
    //                                 uI, uB, uB_uI));
    ierr = xf_Error(ConvFluxBoundaryState(nq, sr, dim, uI, uB, wn, xglob, Model->typeBCs[localBCflag], Model->paraBCs[localBCflag]));
    if (ierr != xf_OK) return ierr;
    
    // multiply wn by quad weights for further use in integration
    for (iq=0;iq<nq; iq++) 
        for (d=0; d<dim; d++)
            wn[iq*dim+d] *= wq[iq];
    
    // zero out Connectivity
    //if (C != NULL) C->GenArray[egrp].rValue[elem][face] = 0.;
    
    /*------------------*/
    /* CONVECTION TERMS */
    /*------------------*/
    // realloc F, F_uI if necessary
    ierr = xf_Error(xf_ReAlloc( (void **) &F, sr*nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(ConvFluxBoundaryFace(nq, sr, dim, uI, uB, wn, F));
    if (ierr != xf_OK) return ierr;
    
    /* Add to R:
     ER_{n,k} += sum_q Phi_{n,q}*F_{q,k}
     */
    if (ER != NULL)
        xf_MTxM_Add(PhiData->Phi, F, nn, nq, sr, ER);
    
    
    pnq = nq; // set previous quad point # for next bface
    
    // Store possibly-altered or resized data back in StaticData
    SD->QuadData    = QuadData;
    SD->PhiData     = PhiData;
    SD->PhiTable    = PhiTable;
    SD->NData       = NData;
    SD->JData       = JData;
    SD->GeomPhiData = GeomPhiData;
    SD->T           = T;
    SD->T2          = T2;
    SD->wn          = wn;
    SD->xglob       = xglob;
    SD->xelem       = xelem;
    SD->uI          = uI;
    SD->guI         = guI;
    SD->uB          = uB;
    SD->F           = F;
    SD->pnq         = pnq;
    SD->EG          = EG;
    if (pSD == NULL){
        // Delete StaticData that we just created
        ierr = xf_Error(xf_DestroyStaticDataBFace(SD));
        if (ierr != xf_OK) return ierr;
    }
    
    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xfYu_CalculateResidualBFaces
int
xfYu_CalculateResidualBFaces(xf_All *All, const xf_Vector *U, xf_Vector *R, 
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
            }
            
            /* Calculate residual on BFace, passing in StaticData  */
            ierr = xf_Error(xf_CalculateResidualBFace(All, ibfgrp, ibface, U, ER, OutputEval, 
                                                      &StaticData, SolverData));
            if (ierr != xf_OK){  // for parallel's sake, do not return immediately if recoverable
            //      if (!xf_CheckRecoverable(ierr, &ReturnError)) return ierr;
            //      continue; // recoverable error occured, move on
                return ierr;
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
xfYu_CalculateResidual(xf_All *All, Yu_Model *pModel, xf_Vector *U, xf_Vector *R, 
                     xf_JacobianMatrix *R_U, xf_SolverData *SolverData)
{
    int ierr;
    int ReturnError = xf_OK;
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
 
    //Yu's model should be pulled in
    Model = pModel;

    // zero out connectivity if it is required
    // not connectivity is really need in Yu cases
    //if ((SolverData != NULL) && (SolverData->CRequired)){
    //    ierr = xf_Error(xf_SetZeroVector(SolverData->C));
    //    if (ierr != xf_OK) return ierr;
    //    C = SolverData->C;
    //}
    
    /* Calculate terms for stabilization (returns immediately if not required) */
    //ierr = xf_Error(xf_CalculateStabilization(All, U, R_U != NULL, SolverData));
    //if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
    
    /* Calculate residual contribution from boundary faces */
    ierr = xf_Error(xfYu_CalculateResidualBFaces(All, U, R, R_U, NULL, 0, NULL, SolverData));
    //if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
    if (ierr != xf_OK) return ierr;
    
    /* Calculate residual contribution from element interiors */
    ierr = xf_Error(xf_CalculateResidualElems(All, U, R, R_U, SolverData));
    //if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
    if (ierr != xf_OK) return ierr;
    
    /* Calculate residual contribution from internal faces (wait for end
     of halo exchange occurs here, unless the wait already occured) */
    ierr = xf_Error(xf_CalculateResidualIFaces(All, U, R, R_U, SolverData));
    //if ((ierr != xf_OK) && (!xf_CheckRecoverable(ierr, &ReturnError)) ) return ierr;
    if (ierr != xf_OK) return ierr;
    
    // make sure halo is exchanged 
    ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
    if (ierr != xf_OK) return ierr;
    
    // create lines if required
    //if ((SolverData != NULL) && (SolverData->CRequired) && (R_U != NULL)){
    //    ierr = xf_Error(xf_CreateLines(All, C, R_U));
    //    if (ierr != xf_OK) return ierr;
        
    //    if (SolverData->SortLines){
    //        ierr = xf_Error(xf_SortLines(All, C, R_U, R));
    //        if (ierr != xf_OK) return ierr;
    //    }
    //}
    
    return ReturnError;
    
}
