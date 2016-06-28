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
  FILE:  xf_SolverDGTime.c

  This file contains functions for DG in time solvers.

*/

#include "xf_MeshMotionGCL.h"
#include "xf_Partition.h"

/* Global variable (hack) for increasing the number of quadrature
   points in time-slab integrations.  NOTE, if setting this above 0,
   increase xf_MAXDGTIMENODE in xf_SolverStruct.h -- specifically,
   xf_MAXDGTIMENODE should be 3 + # extra time quad points. */
int nExtraTimeQuadPoints = 0;


/* Data structure for error estimation/adaptation */
typedef struct
{  
  int SpaceOrderIncrement;        // spatial order increment
  int TimeOrderIncrement;         // temporal order increment
  int OrderTime;                  // temporal order (coarse)
  int OrderTimeh;                 // fine temporal order
  int nSlab;                      // number of time slabs

  const char *AdaptOnName;        // what are we adapting on?
  const char *SavePrefix;         // prefix for saving files
  char AdaptVariableSet[xf_MAXSTRLEN]; // used for entropy adjoint

  int *sdofTime;                  // spatial dof per time slab

  real **ErrIndTime;              // temporally-localized err ind
  real OutputError;               // output error
  real OutputErrorSpace;          // spatial contribution to output error
  real OutputErrorTime;           // temporal contribution to output error

  enum xfe_Verbosity Verbosity;   // verbosity level
  enum xfe_TimeSchemeType TimeScheme;  // coarse time scheme
  enum xfe_TimeSchemeType TimeSchemeh; // fine time scheme
  enum xfe_SpaceTimeAnisoType UErrEstAnisoMeasure;  // anisotropy measure
  enum xfe_AdaptOnType AdaptOn;   // what are we adapting on?
  enum xfe_Bool UseGCL; // are we using a geometric conservation law?


  xf_Vector *ErrIndElemSTot;  // total spatial-loc err ind (space)
  xf_Vector *ErrIndElem    ;  // slab spatial-loc err ind (both space and time)
  xf_Vector *ErrIndElemS   ;  // slab spatial-loc err ind (space)
  xf_Vector *ErrIndElemT   ;  // spatial-loc err ind (time)
  xf_Vector *SpaceTimePref ;  // spatial vs temporal preference
  xf_Vector **SUj          ;  // state storage vectors
  xf_Vector *SUprev        ;  // previous state storage

  xf_Vector **SVj          ;  // extra state storage vectors (for entropy adjoint)
  xf_Vector *SVprev        ;  // extra previous state storage(for entropy adjoint)

}
xf_DGTimeErrEstData;


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeDGTimeErrEstData
static int 
xf_UnParallelizeDGTimeErrEstData(xf_All *All, xf_DGTimeErrEstData *ErrEstData)
{
  int ierr, terr, i, myRank, nProc;
  xf_Data *D;
  xf_Vector *Temp = NULL;
  xf_Mesh *Mesh = All->Mesh;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  //we only need to unparalellize the vectors
  if (nProc > 1){
    /* ============ */
    //ErrIndElemSTot
    /* ============ */
    if (myRank == 0){//Create a temporary vector for the global information
      ierr = xf_Error(xf_CreateVector(&Temp));
      if (ierr != xf_OK) return ierr;
    }
    //Unparallelize into Temp
    ierr = xf_Error(xf_UnParallelizeVector(Mesh, ErrEstData->ErrIndElemSTot, Temp));
    if (ierr != xf_OK) return ierr;
    //destroy old
    terr = xf_FindDataByPointer(All->DataSet, 
                                (void*)ErrEstData->ErrIndElemSTot, &D);
    if (terr == xf_NOT_FOUND) {//safe to destroy
      ierr = xf_Error(xf_DestroyVector(ErrEstData->ErrIndElemSTot, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    else if (terr == xf_OK){
      ierr = xf_Error(xf_DestroyDataInSet(All->DataSet, D));
      if (ierr != xf_OK) return ierr;
    }
    else return terr;

    //reset the pointer to the temporary (note: Temp is NULL in non-root processors)
    ErrEstData->ErrIndElemSTot = Temp;
    /* ========== */
    //ErrIndElem
    /* ========== */
    if (myRank == 0){
      ierr = xf_Error(xf_CreateVector(&Temp));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_UnParallelizeVector(Mesh, ErrEstData->ErrIndElem, Temp));
    if (ierr != xf_OK) return ierr;
    //destroy old
    terr = xf_FindDataByPointer(All->DataSet, 
                                (void*)ErrEstData->ErrIndElem, &D);
    if (terr == xf_NOT_FOUND) {//safe to destroy
      ierr = xf_Error(xf_DestroyVector(ErrEstData->ErrIndElem, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    else if (terr == xf_OK){
      ierr = xf_Error(xf_DestroyDataInSet(All->DataSet, D));
      if (ierr != xf_OK) return ierr;
    }
    else return terr;
    ErrEstData->ErrIndElem = Temp;
    /* ========== */
    //ErrIndElemS
    /* ========== */
    if (myRank == 0){
      ierr = xf_Error(xf_CreateVector(&Temp));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_UnParallelizeVector(Mesh, ErrEstData->ErrIndElemS, Temp));
    if (ierr != xf_OK) return ierr;
    //destroy old
    terr = xf_FindDataByPointer(All->DataSet, 
                                (void*)ErrEstData->ErrIndElemS, &D);
    if (terr == xf_NOT_FOUND) {//safe to destroy
      ierr = xf_Error(xf_DestroyVector(ErrEstData->ErrIndElemS, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    else if (terr == xf_OK){
      ierr = xf_Error(xf_DestroyDataInSet(All->DataSet, D));
      if (ierr != xf_OK) return ierr;
    }
    else return terr;
    ErrEstData->ErrIndElemS = Temp;
    /* ========== */
    //ErrIndElemT
    /* ========== */
    if (myRank == 0){
      ierr = xf_Error(xf_CreateVector(&Temp));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_UnParallelizeVector(Mesh, ErrEstData->ErrIndElemT, Temp));
    if (ierr != xf_OK) return ierr;
    //destroy old
    terr = xf_FindDataByPointer(All->DataSet, 
                                (void*)ErrEstData->ErrIndElemT, &D);
    if (terr == xf_NOT_FOUND) {//safe to destroy
      ierr = xf_Error(xf_DestroyVector(ErrEstData->ErrIndElemT, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    else if (terr == xf_OK){
      ierr = xf_Error(xf_DestroyDataInSet(All->DataSet, D));
      if (ierr != xf_OK) return ierr;
    }
    else return terr;
    ErrEstData->ErrIndElemT = Temp;
    /* ========== */
    //SpaceTimePref
    /* ========== */
    if (myRank == 0){
      ierr = xf_Error(xf_CreateVector(&Temp));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_UnParallelizeVector(Mesh, ErrEstData->SpaceTimePref, Temp));
    if (ierr != xf_OK) return ierr;
    //destroy old
    terr = xf_FindDataByPointer(All->DataSet, 
                                (void*)ErrEstData->SpaceTimePref, &D);
    if (terr == xf_NOT_FOUND) {//safe to destroy
      ierr = xf_Error(xf_DestroyVector(ErrEstData->SpaceTimePref, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    else if (terr == xf_OK){
      ierr = xf_Error(xf_DestroyDataInSet(All->DataSet, D));
      if (ierr != xf_OK) return ierr;
    }
    else return terr;
    ErrEstData->SpaceTimePref = Temp;
    
    
    //SUj
    if (ErrEstData->SUj != NULL){
      for (i = 0; i < ErrEstData->OrderTimeh+1; i++){
        ierr = xf_Error(xf_DestroyVector(ErrEstData->SUj[i], xfe_True));
        if (ierr != xf_OK) return ierr;
        ErrEstData->SUj[i] = NULL;
      }
    }
    //SVj
    if (ErrEstData->SVj != NULL){
      for (i = 0; i < ErrEstData->OrderTimeh+1; i++){
        ierr = xf_Error(xf_DestroyVector(ErrEstData->SVj[i], xfe_True));
        if (ierr != xf_OK) return ierr;
        ErrEstData->SVj[i] = NULL;
      }
    }
    ErrEstData->SUprev = NULL;
    ErrEstData->SVprev = NULL;
    //NOTE: we are still keeping ErrEstData->SVj[i] and ErrEstData->SUj[i]
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeDGTimeErrEstData
static int 
xf_ParallelizeDGTimeErrEstData(xf_All *All, xf_DGTimeErrEstData *ErrEstData)
{
  int ierr, myRank, nProc;
  char Title[xf_MAXSTRLEN];
  xf_Vector *Temp = NULL;
  xf_Mesh *Mesh = All->Mesh;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  //everything stays the same, only xf_Vectors need to be parallelize
  if (nProc > 1){
    /* ============ */
    //ErrIndElemSTot
    /* ============ */
    ierr = xf_Error(xf_CreateVector(&Temp));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ParallelizeVector(Mesh, ErrEstData->ErrIndElemSTot, Temp));
    if (ierr != xf_OK) return ierr;
    //ErrIndElemSTot should be null in non-root processors
    ierr = xf_Error(xf_DestroyVector(ErrEstData->ErrIndElemSTot, xfe_True));
    if (ierr != xf_OK) return ierr;
    ErrEstData->ErrIndElemSTot = Temp;
    sprintf(Title, "ErrIndSpaceTot_%s", ErrEstData->AdaptOnName);
    ierr = xf_Error(xf_DataSetAdd(All->DataSet, Title, xfe_Vector, xfe_True, 
                                  (void*)ErrEstData->ErrIndElemSTot, NULL));
    if (ierr != xf_OK) return ierr;
    
    
    /* ============ */
    //ErrIndElem
    /* ============ */
    ierr = xf_Error(xf_CreateVector(&Temp));
    if (ierr != xf_OK) return ierr;    
    ierr = xf_Error(xf_ParallelizeVector(Mesh, ErrEstData->ErrIndElem, Temp));
    if (ierr != xf_OK) return ierr;
    //ErrIndElem should be null in non-root processors
    ierr = xf_Error(xf_DestroyVector(ErrEstData->ErrIndElem, xfe_True));
    if (ierr != xf_OK) return ierr;
    ErrEstData->ErrIndElem = Temp;
    
    /* ============ */
    //ErrIndElemS
    /* ============ */
    ierr = xf_Error(xf_CreateVector(&Temp));
    if (ierr != xf_OK) return ierr;    
    ierr = xf_Error(xf_ParallelizeVector(Mesh, ErrEstData->ErrIndElemS, Temp));
    if (ierr != xf_OK) return ierr;
    //ErrIndElemS should be null in non-root processors
    ierr = xf_Error(xf_DestroyVector(ErrEstData->ErrIndElemS, xfe_True));
    if (ierr != xf_OK) return ierr;
    ErrEstData->ErrIndElemS = Temp;
    
    /* ============ */
    //ErrIndElemT
    /* ============ */
    ierr = xf_Error(xf_CreateVector(&Temp));
    if (ierr != xf_OK) return ierr;    
    ierr = xf_Error(xf_ParallelizeVector(Mesh, ErrEstData->ErrIndElemT, Temp));
    if (ierr != xf_OK) return ierr;
    //ErrIndElemT should be null in non-root processors
    ierr = xf_Error(xf_DestroyVector(ErrEstData->ErrIndElemT, xfe_True));
    if (ierr != xf_OK) return ierr;
    ErrEstData->ErrIndElemT = Temp;
    
    /* ============ */
    //SpaceTimePref
    /* ============ */
    ierr = xf_Error(xf_CreateVector(&Temp));
    if (ierr != xf_OK) return ierr;    
    ierr = xf_Error(xf_ParallelizeVector(Mesh, ErrEstData->SpaceTimePref, Temp));
    if (ierr != xf_OK) return ierr;
    //SpaceTimePref should be null in non-root processors
    ierr = xf_Error(xf_DestroyVector(ErrEstData->SpaceTimePref, xfe_True));
    if (ierr != xf_OK) return ierr;
    ErrEstData->SpaceTimePref = Temp;
    
  }
  
  return xf_OK;
}

/******************************************************************/
// FUNCTION Definition: xf_DGTimeScheme2Order
int
xf_DGTimeScheme2Order(enum xfe_TimeSchemeType TimeScheme, int *OrderTime)
{
  switch (TimeScheme){
  case xfe_TimeSchemeDG1:  
    (*OrderTime) = 1; 
    break;
  case xfe_TimeSchemeDG2:
    (*OrderTime) = 2; 
    break;
  default: 
    return xf_Error(xf_NOT_SUPPORTED); 
    break;
  }
  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_DGTimeSchemeVars
int
xf_DGTimeSchemeVars(enum xfe_TimeSchemeType TimeScheme, real *AT, 
		    real *VP, int *nq, real *tq, real *wq, real *iAT)
{
  int iq, nq0;

  // 2 points
  real tq2[] = {0.211324865405187, 0.788675134594813};
  real wq2[] = {0.500000000000000, 0.500000000000000};
  
  // 3 points
  real tq3[] = {0.112701665379258, 0.500000000000000, 0.887298334620742};
  real wq3[] = {0.277777777777778, 0.444444444444444, 0.277777777777778};
  
  // 4 points
  real tq4[] = {0.069431844202974, 0.330009478207572, 0.669990521792428, 0.930568155797026};
  real wq4[] = {0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727};

  // 5 points
  real tq5[] = {0.046910077030668, 0.230765344947158, 0.500000000000000, 0.769234655052841,
                0.953089922969332};
  real wq5[] = {0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683,
                0.118463442528095};

  real *tqv[5], *wqv[5];

  // set pointers to points and weights
  tqv[0] = NULL;   wqv[0] = NULL;
  tqv[1] = NULL;   wqv[1] = NULL;
  tqv[2] =  tq2;   wqv[2] =  wq2;
  tqv[3] =  tq3;   wqv[3] =  wq3;
  tqv[4] =  tq4;   wqv[4] =  wq4;
  tqv[5] =  tq5;   wqv[5] =  wq5;

  switch (TimeScheme){
  case xfe_TimeSchemeDG1:  
    if (AT != NULL){
      AT[0] =  0.5; AT[1] = 0.5;
      AT[2] = -0.5; AT[3] = 0.5;
    }
    if(iAT != NULL){
      iAT[0] = 1.; iAT[1] = -1.;
      iAT[2] = 1.; iAT[3] = 1.;
    }
    if (VP != NULL){
      VP[0] = -1.0; VP[1] = 0.;
    }
    nq0 = 2;
    //if (tq != NULL){ tq[0] = (1-1./sqrt(3))/2.; tq[1] = (1+1./sqrt(3))/2.;}
    //if (wq != NULL){ wq[0] =               0.5; wq[1] =               0.5;}
    break;
  case xfe_TimeSchemeDG2:
    if (AT != NULL){
      AT[0] =  1./2.; AT[1] =  2./3.; AT[2] = -1./6.;
      AT[3] = -2./3.; AT[4] =     0.; AT[5] =  2./3.;
      AT[6] =  1./6.; AT[7] = -2./3.; AT[8] =  1./2.;
    }
    if(iAT != NULL){
      iAT[0] = 1.; iAT[1] = -1./2.; iAT[2] = 1.;
      iAT[3] = 1.; iAT[4] = 5./8.;  iAT[5] = -1./2.;
      iAT[6] = 1.; iAT[7] = 1.;     iAT[8] = 1.;
    }
    if (VP != NULL){
      VP[0] = -1.0; VP[1] = 0.; VP[2] = 0.;
    }
    nq0 = 3;
    //if (tq != NULL){tq[0] = (1.-sqrt(0.6))/2.; tq[1] =      0.5; tq[2] = (1.+sqrt(0.6))/2.;}
    //if (wq != NULL){wq[0] =          5.0/18.0; wq[1] = 8.0/18.0; wq[2] = 5.0/18.0;}
    break;
  default: 
    return xf_Error(xf_NOT_SUPPORTED); 
    break;
  }

  // add quad point delta
  nq0 += nExtraTimeQuadPoints;

  // error check number of points
  if ((nq0<2) || (nq0>5)) return xf_Error(xf_NOT_SUPPORTED);

  // set number of points
  if (nq != NULL) (*nq) = nq0;

  // set points and weights
  if (tq != NULL)  for (iq=0; iq<nq0; iq++) tq[iq] = tqv[nq0][iq];
  if (wq != NULL)  for (iq=0; iq<nq0; iq++) wq[iq] = wqv[nq0][iq];


  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_DGTimeInterpolate
int
xf_DGTimeInterpolate(enum xfe_TimeSchemeType TimeScheme, xf_Vector **Uj,
		     int nt, real *xt, xf_Vector **Vj, real *phiout )
{
  int ierr, OrderTime, i, it;
  real dx, xnode[xf_MAXDGTIMENODE];
  real phi[xf_MAXDGTIMENODE];

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;

  // assume a Lagrange basis
  for (i=0, dx=1./((real) OrderTime); i<OrderTime+1; i++) xnode[i] = i*dx;
  
  // loop over time points
  for (it=0; it<nt; it++){
    // pull off time basis functions at time point
    ierr = xf_Error(xf_BasisLagrange1D(xt[it], xnode, OrderTime+1, phi, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    
    // interpolate Uj
    for (i=0; i<OrderTime+1; i++){
      ierr = xf_Error(xf_VectorMultSet(Uj[i], phi[i], ((i==0) ? xfe_Set : xfe_Add), Vj[it]));
      if (ierr != xf_OK) return ierr;
    }

    // on last time node, store phiout, if requested
    if ((it==(nt-1)) && (phiout != NULL)){
      for (i=0; i<OrderTime+1; i++) phiout[i] = phi[i];
    }

  } // it

  return xf_OK;
}

/******************************************************************/
// FUNCTION Definition: xf_DGTimeInterpolateState
int
xf_DGTimeInterpolateState(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                          xf_Vector **Uj, real *xt, xf_Vector **Vj, real *phiout)
{
  int ierr;
  enum xfe_Bool UseGCL;
  
  
  //Determine if using a Geometric Conservation Law
  if(All != NULL){
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
    if (ierr != xf_OK) return ierr;
  }
  
  //If state is requested at only a single time, interpolate GCL to that time as well
  if((All != NULL) && UseGCL){
    ierr = xf_Error(xf_DGTimeInterpolateGCL(All, TimeScheme, (*xt), phiout));
    if (ierr != xf_OK) return ierr;
  }
    
  //Interpolate state vectors to temporal point, as requested
  ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Uj, 1, xt, Vj, phiout));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_DGTimeInterpolateDown
static int
xf_DGTimeInterpolateDown(enum xfe_TimeSchemeType TimeScheme, xf_Vector **Uj,
			 int nt, real *xt, int OrderDecrement, xf_Vector **Vj)
{
  int ierr, OrderTime, OrderTimeH, i, j, it;
  real dx, xnode[xf_MAXDGTIMENODE];
  real phi[xf_MAXDGTIMENODE];
  real T[6] = {2./3., 2./3., -1./3., -1./3., 2./3., 2./3.};
  real S[3];

  if (OrderDecrement == 0)
    return xf_Error(xf_DGTimeInterpolate(TimeScheme, Uj, nt, xt, Vj, NULL));

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;

  // for now support only DG2 to DG1
  if (OrderTime != 2) return xf_Error(xf_NOT_SUPPORTED);
  if (OrderDecrement != 1) return xf_Error(xf_NOT_SUPPORTED);
  OrderTimeH = 1;

  // assume a Lagrange basis
  for (i=0, dx=1./((real) OrderTimeH); i<OrderTimeH+1; i++) xnode[i] = i*dx;
  
  // loop over time points
  for (it=0; it<nt; it++){
    // pull off time basis functions at time point
    ierr = xf_Error(xf_BasisLagrange1D(xt[it], xnode, OrderTimeH+1, phi, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    for (j=0; j<OrderTime+1; j++) S[j] = 0.;
    for (i=0; i<OrderTimeH+1; i++)
      for (j=0; j<OrderTime+1; j++)
	S[j] += T[i*(OrderTime+1)+j]*phi[i];
    
    // interpolate Uj
    for (j=0; j<OrderTime+1; j++){
      ierr = xf_Error(xf_VectorMultSet(Uj[j], S[j], ((j==0) ? xfe_Set : xfe_Add), Vj[it]));
      if (ierr != xf_OK) return ierr;
    }

  } // it

  return xf_OK;
}



/******************************************************************/
// FUNCTION Definition: xf_DGTimeInject
static int
xf_DGTimeInject(int OrderTimeH, int OrderTimeh, xf_Vector **Vj)
{
/*
PURPOSE:

  Injects Vj from temporal order OrderTimeH to OrderTimeh

INPUTS:
 
  OrderTimeH   : coarse (current) order
  OrderTimeh   : fine (desired) order
  Vj           : vector to project
  
OUTPUTS: 

  Vj           : projected vector

RETURN: Error code

*/
  int ierr;
  
  if ((OrderTimeH == 1) && (OrderTimeh == 2)){
  
    // Vj[2] = Vj[1]
    ierr = xf_Error(xf_VectorMultSet(Vj[1], 1.0, xfe_Set, Vj[2]));
    if (ierr != xf_OK) return ierr;
    
    // Vj[1] = 0.5(Vj[0] + Vj[2]);
    ierr = xf_Error(xf_VectorMultSet(Vj[0], 0.5, xfe_Set, Vj[1]));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_VectorMultSet(Vj[2], 0.5, xfe_Add, Vj[1]));
    if (ierr != xf_OK) return ierr;

  }
  else if ((OrderTimeH == 2) && (OrderTimeh == 1)){
    // "injecting" to a lower order is only valid if Uj is a low-order polynomial (in t)

    // Vj[1] = Vj[2]
    ierr = xf_Error(xf_VectorMultSet(Vj[2], 1.0, xfe_Set, Vj[1]));
    if (ierr != xf_OK) return ierr;
    
  }
  else return xf_Error(xf_NOT_SUPPORTED);

  
  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_DGTimeReconstruct
static int
xf_DGTimeReconstruct(xf_All *All, int OrderTimeH, int OrderTimeh, xf_Vector **Vj, 
		     xf_Vector *Vnpin, enum xfe_Bool LeftRadau)
{
/*
PURPOSE:

  Reconstructs Vj from temporal order OrderTimeH to OrderTimeh taking
  advantage of super convergence of left or right Radau points.  Uses
  only the right-node solution from the previous time slab (for Right
  Radau points, forward solve) or the left-node solution from the next
  time slab (for Left Radau points, adjoint solve).

INPUTS:
 
  All          : all structure
  OrderTimeH   : coarse (current) order
  OrderTimeh   : fine (desired) order
  Vj           : vector to reconstruct
  Vnp          : if LeftRadau == False, this should be the right-node
                 solution from the previous time slab (e.g. Uprev)
		 if LeftRadau == True, this should be the left-node
		 solution from the next time slab (e.g. Adjnext)
  LeftRadau    : boolean indicating forward versus adjoint solve
  
OUTPUTS: 

  Vj           : reconstructed vector

RETURN: Error code

*/
  int ierr;
  real x;
  xf_Vector *Vnp = NULL;
  
  if ((OrderTimeH != 1) && (OrderTimeh != 2))
    return xf_Error(xf_NOT_SUPPORTED);
  
  // forward case is easy to put in, but do not need it yet
  if (LeftRadau != xfe_True) return xf_Error(xf_NOT_SUPPORTED);

  // project Vnpin to space of Vj[0] if incompatible
  if (!xf_CompatibleVectors(Vnpin, Vj[0])){
    ierr = xf_Error(xf_FindSimilarVector(All, Vj[0], "DGReconTemp", 
					 xfe_False, xfe_True, NULL, 
					 &Vnp, NULL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ProjectVector(All, Vnpin, xfe_False, Vnp));
    if (ierr != xf_OK) return ierr;
  }
  else Vnp = Vnpin;

  // Vj[2] = DG1 interpolation of (Vj[0], Vj[1]) at x = 2/3
  x = 2./3.;
  ierr = xf_Error(xf_DGTimeInterpolate(xfe_TimeSchemeDG1, Vj, 1, &x, Vj+2, NULL));
  if (ierr != xf_OK) return ierr;

  // Vj[1] = Vj[0]*1/8 + Vj[2]*9/8 + Vnp*(-1/4)
  ierr = xf_Error(xf_VectorMultSet(Vj[0],  1.0/8.0, xfe_Set, Vj[1]));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(Vj[2],  9.0/8.0, xfe_Add, Vj[1]));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(Vnp  , -1.0/4.0, xfe_Add, Vj[1]));
  if (ierr != xf_OK) return ierr;

  // Vj[2] = Vnp
  ierr = xf_Error(xf_VectorMultSet(Vnp  , 1.0, xfe_Set, Vj[2]));
  if (ierr != xf_OK) return ierr;
    
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DGinTimeUnsteadyResidual
static int
xf_DGinTimeUnsteadyResidual(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
			    xf_Vector **Uj, xf_Vector *Uprev, real Time, 
			    real TimeStep, xf_SolverData *SolverData, 
			    xf_Vector *Utemp, xf_Vector *Rtemp, 
			    xf_Vector **Ri, enum xfe_Bool *pRewind)
{
/*
PURPOSE:

  Calculates the unsteady DG residual for the GCL.


INPUTS:

  All         : all structure
  Uj          : states on the current time slab
  Uprev       : state vector at the right node of the previous time slab
  TimeScheme  : temporal time scheme
  Time        : current time
  TimeStep    : delta t
  SolverData  : solver data structure (contains CFL)
  Utemp, Rtemp: temporary state and residual vectors
  
  
OUTPUTS: 

  Ri          : unsteady residual
  (*pRewind)  : True if a solution rewind is requested (e.g. non-physical)

RETURN: Error code

*/
  int ierr;
  int iOrder, jOrder, OrderTime;
  int iq, nq;
  real tq[xf_MAXDGTIMENODE], wq[xf_MAXDGTIMENODE], phi[xf_MAXDGTIMENODE]; 
  real Atime[xf_MAXDGTIMENODE*xf_MAXDGTIMENODE], Vprev[xf_MAXDGTIMENODE];

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;

  // pull off variables for this time scheme
  ierr = xf_Error(xf_DGTimeSchemeVars(TimeScheme, Atime, Vprev, &nq, tq, wq, NULL));
  if (ierr != xf_OK) return ierr;
 
  /* 
     DG1:
     R0 = M*( 1/2*U0+1/2*U1)-M'*Uprev + int_0^dt phi0*Rspace(U) dt
     R1 = M*(-1/2*U0+1/2*U1)          + int_0^dt phi1*Rspace(U) dt

     DG2:     
     R0 = M*( 1/2*U0+2/3*U1-1/6*U2)-M'*Uprev + int_0^dt phi0*Rspace(U) dt
     R1 = M*(-2/3*U0        2/3*U2)          + int_0^dt phi1*Rspace(U) dt
     R2 = M*( 1/6*U0-2/3*U1+1/2*U2)          + int_0^dt phi2*Rspace(U) dt

     M' = possibly multi-order mass matrix if U0/U0 and Uprev have different orders
     
     Construction of the unsteady residual with Utemp, Rtemp (e.g. DG1):
   
     Utemp = U at first quadrature point
     Rtemp = R(Utemp)
     Ri[0] += factor00*Rtemp;
     Ri[1] += factor01*Rtemp;
     
     Utemp = U at second quadrature point
     Rtemp = R(Utemp)
     Ri[0] += factor10*Rtemp;
     Ri[1] += factor11*Rtemp;
     
     Then add in M*U etc. contributions to Ri[0] and Ri[1] (one at
     a time), using only Utemp.
     
  */
  
  // initialize Ri to zero
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    ierr = xf_Error(xf_SetZeroVector(Ri[iOrder]));
    if (ierr != xf_OK) return ierr;
  }
  
  // loop over the quadrature points
  for (iq=0; iq<nq; iq++){
    // Set Time to that at quadrature point, Time + tq[iq]*TimeStep
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time",
				       Time+tq[iq]*TimeStep));
    if (ierr != xf_OK) return ierr;
    
    // Calculate U at the quad point Utemp = U(tq), and phii(tq)
    ierr = xf_Error(xf_DGTimeInterpolateState(All,TimeScheme, Uj, tq+iq, &Utemp, phi));
    if (ierr != xf_OK) return ierr;
	
    // Calculate the spatial residual at the iq quadrature state -> Rtemp
    ierr = xf_CalculateResidual(All, Utemp, Rtemp, NULL, SolverData);
    if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, &Uprev,
				    OrderTime+1, Uj, &TimeStep, pRewind))
      return xf_Error(xf_SOLVER_ERROR);
    if (*pRewind){ if (SolverData->iIter==0) return xf_Error(xf_NON_PHYSICAL);
      else break;}
	
    // Ri[i] += wq(iq)*phii(tq)*TimeStep*Rtemp
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      ierr = xf_Error(xf_VectorMultSet(Rtemp, wq[iq]*phi[iOrder]*TimeStep,
				       xfe_Add, Ri[iOrder]));
      if (ierr != xf_OK) return ierr;
    } 
  }
  
  // Take care of temporal stiffness matrix and influence of Uprev
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    // temporal stiffness matrix
    for (jOrder = 0; jOrder <= OrderTime; jOrder++){
      ierr = xf_Error(xf_VectorMultSet(Uj[jOrder],  Atime[iOrder*(OrderTime+1)+jOrder], 
				       ((jOrder==0) ? xfe_Set : xfe_Add), Utemp));
      if (ierr != xf_OK) return ierr;
    }
    // add M*Utemp to temporal residual Ri[iOrder]
    ierr = xf_Error(xf_AddMassMatrix(All, 1.0, NULL, Utemp, Ri[iOrder], NULL, NULL));
    if (ierr != xf_OK) return ierr;
    // influence of Uprev
    if (Vprev[iOrder] != 0.){
      // add M'*Uprev to temporal residual Ri[iOrder]
      ierr = xf_Error(xf_AddMassMatrix(All, Vprev[iOrder], NULL, Uprev, Ri[iOrder], NULL, NULL));
      if (ierr != xf_OK) return ierr;
    }
  }

  return xf_OK;
}



/******************************************************************/
// FUNCTION Definition: xf_DG1ApproxSolver
static int
xf_DG1ApproxSolver(xf_All *All, real Time, real TimeStep, 
		   xf_SolverData *SolverData, xf_Vector **U, 
		   xf_Vector **R, xf_Vector **X, 
		   enum xfe_Bool TransposeFlag)
{
/*

PURPOSE

  This function solves the linear system of equations arising from a
  DG(r=1) discretization approximately, using Thomas Ritcher's method.
  The system is:

  [ M/2 + dt*A/3]*X[0] + [M/2 + dt*A/6]*X[1] = -R[0]
  [-M/2 + dt*A/6]*X[0] + [M/2 + dt*A/3]*X[1] = -R[1]

  The solution steps are:

    [M+dt/sqrt(6)*A]*Y    +  R[0] + R[1] + A*M^{-1}*dt*(2*R[1]-R[0])/3 = 0
    [M+dt/sqrt(6)*A]*X[1] - M*Y = 0
    [M+dt*2/3    *A]*X[0] + 2*R[0] + [M+dt/3*A]X[1] = 0

  The transpose problem is (note M^T = M):

  [M/2 + dt*A^T/3]*X[0] + [-M/2 + dt*A^T/6]*X[1] = -R[0]
  [M/2 + dt*A^T/6]*X[0] + [ M/2 + dt*A^T/3]*X[1] = -R[1]

  And the transpose solution steps are:

    [M+dt/sqrt(6)*A^T]*Y    = -R[1] - R[0] - A^T*M^{-1}*dt*(2*R[0]-R[1])/3
    [M+dt/sqrt(6)*A^T]*X[0] = M*Y
    [M+dt*2/3    *A^T]*X[1] = -2*R[1] - [M+dt/3*A^T]X[0]


INPUTS:

  All         : All structure
  Time        : Time at start of time slab
  TimeStep    : Delta t = time at end of slab - time at start
  SolverData  : solver data structure (contains CFL)
  U           : The states at the two ends of the time slab
  R           : (negative) right hand side of the equation (two vectors)
  Transpose   : If true, solve the transpose of the problem (see above)

OUTPUTS: 

  X           : solution vector

RETURN: Error code

*/
  int ierr;
  enum xfe_Bool Found, Rewind;
  enum xfe_Verbosity Verbosity;
  enum xfe_Bool UseGCL;
  xf_Vector  *Utemp, *Rtemp, *Vtemp;
  xf_JacobianMatrix *R_U;
  
  //Determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;

  /* Allocate space for vectors and matrices */

  // locate temporary state vector, Utemp
  ierr = xf_Error(xf_FindSimilarVector(All, U[0], "Utemp_DG", xfe_True, xfe_True,
				       NULL, &Utemp, &Found));
  if (ierr != xf_OK) return ierr;
  if (!Found) return xf_Error(xf_CODE_LOGIC_ERROR); // vector should already be in All
  
  // locate array for temporary state vector, Rtemp
  ierr = xf_Error(xf_FindSimilarVector(All, Utemp, "Rtemp_DG", xfe_True, xfe_True,
				       NULL, &Rtemp, &Found));
  if (ierr != xf_OK) return ierr;
  if (!Found) return xf_Error(xf_CODE_LOGIC_ERROR); // vector should already be in All



  // locate Jacobian matrix
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, Utemp, NULL,
					xfe_True, NULL, &R_U, &Found));
  if (ierr != xf_OK) return ierr;
  

  // If TransposeFlag == true, switch the input R[0] and R[1]
  if (TransposeFlag) swap(R[0], R[1], Vtemp);
 
  /* Begin linear solve */

  // Set Time to slab midpoint, Time + TimeStep/2
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time+TimeStep/2.0));
  if (ierr != xf_OK) return ierr;
     
  // Calculate the time slab midpoint state for calculation of Jacobian matrix
  ierr = xf_Error(xf_VectorMultSet(U[0], 0.5, xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(U[1], 0.5, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;

  // Calculate Jacobian with Utemp, 
  
  if(UseGCL){
    // First, interpolate GCL to time of residual evaluation, tref=0.5
    ierr = xf_Error(xf_DGTimeInterpolateGCL(All, xfe_TimeSchemeDG1, 0.5, NULL));
    if(ierr != xf_OK) return ierr;
  }
  
  ierr = xf_CalculateResidual(All, Utemp, Rtemp, R_U, SolverData);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);

  // build rhs (on left) of first equation in three-step system
  ierr = xf_Error(xf_VectorMultSet(R[1], 2.*TimeStep/3., xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[0], -1.*TimeStep/3., xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MultInvMassMatrix(All, 1., NULL, Utemp));
  if (ierr != xf_OK) return ierr;

  // Calculate A*Utemp, store in Rtemp
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Utemp, xfe_Set, TransposeFlag, 
				   SolverData, Rtemp));
  if (ierr != xf_OK) return ierr;
      
  // add Ri[0] and Ri[1]
  ierr = xf_Error(xf_VectorMultSet(R[0], 1., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[1], 1., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;

  // the entire equation is multiplied by this factor
  ierr = xf_Error(xf_VectorMult(Rtemp, sqrt(6.0)/TimeStep));
  if (ierr != xf_OK) return ierr;

  // add c*M to Jacobian
  ierr = xf_Error(xf_AddMassMatrix(All, sqrt(6.0)/TimeStep, NULL, U[0], NULL, R_U, NULL));
  if (ierr != xf_OK) return ierr;

  // solve the 1st linear system, and save the result in Utemp
  ierr = xf_SolveLinearSystem(All, R_U, Rtemp, TransposeFlag, -1, SolverData, Utemp);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);

  // calculate the rhs (on left) of the 2nd linear system
  ierr = xf_Error(xf_MultMassMatrix(All, -sqrt(6.0)/TimeStep, Utemp));
  if (ierr != xf_OK) return ierr;

  // solve the 2nd linear system, and save the result in X[1]
  ierr = xf_SolveLinearSystem(All, R_U, Utemp, TransposeFlag, -1, SolverData, X[1]);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);

  // calculate the 3rd linear system

  // Recalculate the Jacobian using the time slab midpoint
  ierr = xf_Error(xf_VectorMultSet(U[0], 0.5, xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(U[1], 0.5, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
      
  // Calculate Jacobian with Utemp
  if(UseGCL){
    // First, interpolate GCL to time of residual evaluation, tref=0.5
    ierr = xf_Error(xf_DGTimeInterpolateGCL(All, xfe_TimeSchemeDG1, 0.5, NULL));
    if(ierr != xf_OK) return ierr;
  }
  
  ierr = xf_CalculateResidual(All, Utemp, Rtemp, R_U, SolverData);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);

  // set Rtemp = X[1]
  ierr = xf_Error(xf_VectorMultSet(X[1], 1.0, xfe_Set, Rtemp));
  if (ierr != xf_OK) return ierr;

  // Utemp = A*Rtemp  (note Rtemp and X[1] are the same at this point)
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Rtemp, xfe_Set, TransposeFlag, 
				   SolverData, Utemp));
  if (ierr != xf_OK) return ierr;
	  
  // continue building rhs (on left) of third system
  ierr = xf_Error(xf_VectorMult(Utemp, TimeStep/3.0));
  if (ierr != xf_OK) return ierr;

  // Rtemp = M*Rtemp
  ierr = xf_Error(xf_MultMassMatrix(All, 1.0, Rtemp));
  if (ierr != xf_OK) return ierr;
  
  // Utemp += Rtemp
  ierr = xf_Error(xf_VectorMultSet(Rtemp, 1.0, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  
  // Utemp += 2*R[0]
  ierr = xf_Error(xf_VectorMultSet(R[0], 2.0, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
      
  // whole equation is multiplied by 3/(2*dt)
  ierr = xf_Error(xf_VectorMult(Utemp, 3/(2*TimeStep)));
  if (ierr != xf_OK) return ierr;

  // add M*3/(2*dt) to the Jacobian
  ierr = xf_Error(xf_AddMassMatrix(All, 3/(2*TimeStep), NULL, Utemp, NULL, R_U, NULL));
  if (ierr != xf_OK) return ierr;
		      
  // solve the 3rd linear system, and save the result in X[0]
  ierr = xf_SolveLinearSystem(All, R_U, Utemp, TransposeFlag, -1, SolverData, X[0]);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);


  // if TransposeFlag == True, switch the vector in X
  if (TransposeFlag) swap(X[0], X[1], Vtemp);

  return xf_OK;
}



/******************************************************************/
// FUNCTION Definition: xf_DG2ApproxSolver
static int
xf_DG2ApproxSolver(xf_All *All, real Time, real TimeStep, 
		   xf_SolverData *SolverData, xf_Vector **U, 
		   xf_Vector **R, xf_Vector **X, 
		   enum xfe_Bool TransposeFlag)
{
/*

PURPOSE

  This function solves the linear system of equations arising from a
  DG(r=2) discretization approximately, using Thomas Ritcher's method.
  The system is:

  [ 1/2*M+2/15*dtA]*X0 + [ 2/3*M+1/15*dtA]*X1 + [-1/6*M-1/30*dtA]*X2 = -R0
  [-2/3*M+1/15*dtA]*X0 + [       8/15*dtA]*X1 + [ 2/3*M+1/15*dtA]*X2 = -R1
  [ 1/6*M-1/30*dtA]*X0 + [-2/3*M+1/15*dtA]*X1 + [ 1/2*M+2/15*dtA]*X2 = -R2

  where dtA = dt*A and A = Jacobian at slab midpoint and U=.25*U0+.5*U1+.25*U2.
  The solution steps are:
  
    [M+   60^{-1/3}*dtA]*B  + M*G0 = 0
    [M+   60^{-1/3}*dtA]*Y  - M*B  = 0
    [M+   60^{-1/3}*dtA]*X2 - M*Y  = 0
    [M+sqrt(3/20)*dtA]*Z  + M*G1 = 0
    [M+sqrt(3/20)*dtA]*X1 - M*Z  = 0
    [M+      4/15*dtA]*X0 + M*G2 = 0

  where

    M*G0 = R0+R1+R2 + dtA*M^{-1}*[-2/5*R0+1/10*R1+3/5*R2
                                  +dtA*M^{-1}*(1/20*R0-1/40*R1+3/20*R2)]
    M*G1 = 3/2*R0+9/8*R1+1/2*M*X2 + dtA*[M^{-1}*(-3/20*R0+3/10*R1+1/40*dtA*X2)
                                         + 1/4*X2]
    M*G2 = 2*R0+M*(4/3*X1-1/3*X2) + dtA*[2/15*X1-1/15*X2]

  The transpose problem has dtA^T instead of dtA, and M^T=M.  The
  transpose solution steps are the same as the forward ones with the
  roles of R0/R2 and U0/U2 interchanged.


INPUTS:

  All         : All structure
  Time        : Time at start of time slab
  TimeStep    : Delta t = time at end of slab - time at start
  SolverData  : solver data structure (contains CFL)
  U           : The states at the Lagrange nodes of the time slab
  R           : (negative) right hand side of the equation (three vectors)
  Transpose   : If true, solve the transpose of the problem (see above)

OUTPUTS: 

  X           : solution vector

RETURN: Error code

*/
  int ierr;
  enum xfe_Bool Found, Rewind;
  enum xfe_Verbosity Verbosity;
  enum xfe_Bool UseGCL;
  real c;
  xf_Vector  *Utemp, *Rtemp, *Vtemp;
  xf_JacobianMatrix *R_U;
  
  //Determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;

  /* Allocate space for vectors and matrices */

  // locate temporary state vector, Utemp
  ierr = xf_Error(xf_FindSimilarVector(All, U[0], "Utemp_DG", xfe_True, xfe_True,
				       NULL, &Utemp, &Found));
  if (ierr != xf_OK) return ierr;
  if (!Found) return xf_Error(xf_CODE_LOGIC_ERROR); // vector should already be in All

  
  // locate array for temporary state vector, Rtemp
  ierr = xf_Error(xf_FindSimilarVector(All, Utemp, "Rtemp_DG", xfe_True, xfe_True,
				       NULL, &Rtemp, &Found));
  if (ierr != xf_OK) return ierr;
  if (!Found) return xf_Error(xf_CODE_LOGIC_ERROR); // vector should already be in All



  // locate Jacobian matrix
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, Utemp, NULL,
					xfe_True, NULL, &R_U, &Found));
  if (ierr != xf_OK) return ierr;
  
  // If TransposeFlag == true, switch the input R[0] and R[2]
  if (TransposeFlag) swap(R[0], R[2], Vtemp);

  // Set Time to slab midpoint, Time + TimeStep/2
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time+TimeStep/2.0));
  if (ierr != xf_OK) return ierr;

  /*** Begin linear solve ***/

  /*-- Equation 1:  [M+   60^{-1/3}*dtA]*B  + M*G0 = 0 --*/
  
  // Calculate the approximate Jacobian at time slab midpoint
  ierr = xf_Error(xf_VectorMultSet(U[0], 0.25, xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(U[1], 0.5 , xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(U[2], 0.25, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  
  if(UseGCL){
    // First, interpolate GCL to time of residual evaluation, tref=0.5
    ierr = xf_Error(xf_DGTimeInterpolateGCL(All, xfe_TimeSchemeDG2, 0.5, NULL));
    if(ierr != xf_OK) return ierr;
  }
  
  ierr = xf_CalculateResidual(All, Utemp, Rtemp, R_U, SolverData);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);

  // build rhs (on left) = M*G0 of first equation
  // M*G0 = R0+R1+R2 + dtA*M^{-1}*[-2/5*R0+1/10*R1+3/5*R2
  //                               +dtA*M^{-1}*(1/20*R0-1/40*R1+3/20*R2)]
  ierr = xf_Error(xf_VectorMultSet(R[0],  1./20., xfe_Set, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[1], -1./40., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[2],  3./20., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MultInvMassMatrix(All, TimeStep, NULL, Rtemp));
  if (ierr != xf_OK) return ierr;
  // Utemp = A*Rtemp
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Rtemp, xfe_Set, TransposeFlag, 
				   SolverData, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[0], -2./5. , xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[1],  1./10., xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[2],  3./5. , xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MultInvMassMatrix(All, TimeStep, NULL, Utemp));
  if (ierr != xf_OK) return ierr;
  // Rtemp = A*Utemp
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Utemp, xfe_Set, TransposeFlag, 
				   SolverData, Rtemp));
  ierr = xf_Error(xf_VectorMultSet(R[0], 1., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[1], 1., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[2], 1., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;

  // the entire equation is multiplied by this factor
  c = pow(60.,1./3.)/TimeStep;
  ierr = xf_Error(xf_VectorMult(Rtemp, c));
  if (ierr != xf_OK) return ierr;
  // add c*M to Jacobian
  ierr = xf_Error(xf_AddMassMatrix(All, c, NULL, U[0], NULL, R_U, NULL));
  if (ierr != xf_OK) return ierr;

  // solve the 1st linear system, and save the result in Utemp
  ierr = xf_SolveLinearSystem(All, R_U, Rtemp, TransposeFlag, -1, SolverData, Utemp);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);


  /*-- Equation 2:  [M+   60^{-1/3}*dtA]*Y  - M*B  = 0 --*/

  // calculate the rhs (on left) of the 2nd linear system, multiplied by c
  ierr = xf_Error(xf_MultMassMatrix(All, -c, Utemp));
  if (ierr != xf_OK) return ierr;

  // solve the 2nd linear system, and save the result in Rtemp
  ierr = xf_SolveLinearSystem(All, R_U, Utemp, TransposeFlag, -1, SolverData, Rtemp);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);


  /*-- Equation 3:  [M+   60^{-1/3}*dtA]*X2 - M*Y  = 0 --*/
  
  // calculate the rhs (on left) of the 2nd linear system, multiplied by c
  ierr = xf_Error(xf_MultMassMatrix(All, -c, Rtemp));
  if (ierr != xf_OK) return ierr;

  // solve the 3rd linear system, and save the result in X[2]
  ierr = xf_SolveLinearSystem(All, R_U, Rtemp, TransposeFlag, -1, SolverData, X[2]);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);


  /*-- Equation 4:   [M+sqrt(3/20)*dtA]*Z  + M*G1 = 0 --*/
 
  // Calculate the approximate Jacobian at time slab midpoint
  ierr = xf_Error(xf_VectorMultSet(U[0], 0.25, xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(U[1], 0.5 , xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(U[2], 0.25, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  
  if(UseGCL){
    // First, interpolate GCL to time of residual evaluation, tref=0.5
    ierr = xf_Error(xf_DGTimeInterpolateGCL(All, xfe_TimeSchemeDG2, 0.5, NULL));
    if(ierr != xf_OK) return ierr;
  }
  
  ierr = xf_CalculateResidual(All, Utemp, Rtemp, R_U, SolverData);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);

  // build rhs (on left) = M*G1
  // M*G1 = 3/2*R0+9/8*R1+1/2*M*X2 + dtA*[M^{-1}*(-3/20*R0+3/10*R1+1/40*dtA*X2)
  //                                      + 1/4*X2]
  ierr = xf_Error(xf_VectorMultSet(X[2], 1./40.*TimeStep, xfe_Set, Rtemp));
  if (ierr != xf_OK) return ierr;
  // Utemp = A*Rtemp
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Rtemp, xfe_Set, TransposeFlag, 
				   SolverData, Utemp));
  ierr = xf_Error(xf_VectorMultSet(R[0], -3./20., xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[1],  3./10., xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MultInvMassMatrix(All, TimeStep, NULL, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(X[2], 1./4.*TimeStep, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  // Rtemp = A*Utemp
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Utemp, xfe_Set, TransposeFlag, 
				   SolverData, Rtemp));
  ierr = xf_Error(xf_VectorMultSet(R[0],  3./2., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[1],  9./8., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetVector(X[2], xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MultMassMatrix(All, 1.0, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(Utemp, 1./2., xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;

  // the entire equation is multiplied by this factor
  c = sqrt(20./3.)/TimeStep;
  ierr = xf_Error(xf_VectorMult(Rtemp, c));
  if (ierr != xf_OK) return ierr;
  // add c*M to Jacobian
  ierr = xf_Error(xf_AddMassMatrix(All, c, NULL, U[0], NULL, R_U, NULL));
  if (ierr != xf_OK) return ierr;

  // solve the 4th linear system, and save the result in Utemp
  ierr = xf_SolveLinearSystem(All, R_U, Rtemp, TransposeFlag, -1, SolverData, Utemp);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);


  /*-- Equation 5:  [M+sqrt(3/20)*dtA]*X1 - M*Z  = 0 --*/

  // calculate the rhs (on left) of the 2nd linear system, multiplied by c
  ierr = xf_Error(xf_MultMassMatrix(All, -c, Utemp));
  if (ierr != xf_OK) return ierr;

  // solve the 5th linear system, and save the result in X[1]
  ierr = xf_SolveLinearSystem(All, R_U, Utemp, TransposeFlag, -1, SolverData, X[1]);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);


  /*-- Equation 6: [M+      4/15*dtA]*X0 + M*G2 = 0 --*/
 
  // Calculate the approximate Jacobian at time slab midpoint
  ierr = xf_Error(xf_VectorMultSet(U[0], 0.25, xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(U[1], 0.5 , xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(U[2], 0.25, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  
  if(UseGCL){
    // First, interpolate GCL to time of residual evaluation, tref=0.5
    ierr = xf_Error(xf_DGTimeInterpolateGCL(All, xfe_TimeSchemeDG2, 0.5, NULL));
    if(ierr != xf_OK) return ierr;
  }
  
  ierr = xf_CalculateResidual(All, Utemp, Rtemp, R_U, SolverData);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);

  // build rhs (on left) = M*G2
  // M*G2 = 2*R0+M*(4/3*X1-1/3*X2) + dtA*[2/15*X1-1/15*X2]

  ierr = xf_Error(xf_VectorMultSet(X[2], -1./15.*TimeStep, xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(X[1],  2./15.*TimeStep, xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  // Rtemp = A*Utemp
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Utemp, xfe_Set, TransposeFlag, 
				   SolverData, Rtemp));
  ierr = xf_Error(xf_VectorMultSet(X[1],  4./3., xfe_Set, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(X[2], -1./3., xfe_Add, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MultMassMatrix(All, 1.0, Utemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(Utemp, 1.0, xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorMultSet(R[0], 2.0, xfe_Add, Rtemp));
  if (ierr != xf_OK) return ierr;

  // the entire equation is multiplied by this factor
  c = 15./(4.*TimeStep);
  ierr = xf_Error(xf_VectorMult(Rtemp, c));
  if (ierr != xf_OK) return ierr;
  // add c*M to Jacobian
  ierr = xf_Error(xf_AddMassMatrix(All, c, NULL, U[0], NULL, R_U, NULL));
  if (ierr != xf_OK) return ierr;

  // solve the 6th linear system, and save the result in X[0]
  ierr = xf_SolveLinearSystem(All, R_U, Rtemp, TransposeFlag, -1, SolverData, X[0]);
  if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
				  0, NULL, &TimeStep, &Rewind)) 
    return xf_Error(xf_SOLVER_ERROR);
  if (Rewind) return xf_Error(xf_REWIND);

  // if TransposeFlag == True, switch the vectors in X
  if (TransposeFlag) swap(X[0], X[2], Vtemp);

  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_DGApproxSolver
static int
xf_DGApproxSolver(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
		  real Time, real TimeStep, xf_SolverData *SolverData, 
		  xf_Vector **U, xf_Vector **R, xf_Vector **X, 
		  enum xfe_Bool TransposeFlag)
{
/*
PURPOSE:

  Wrapper for order-dependent DGinTime approximate solvers.

INPUTS:

  All         : all structure
  TimeScheme  : what time scheme to use
  Time        : current time
  TimeStep    : delta t
  SolverData  : solver data structure (contains CFL)
  U           : The states at the Lagrange nodes of the time slab
  R           : (negative) right hand side of the equation
  Transpose   : If true, solve the transpose of the problem
  
OUTPUTS: 

  X           : solution vector

RETURN: Error code

*/
  
  switch (TimeScheme){
  case xfe_TimeSchemeDG1:  
    return xf_Error(xf_DG1ApproxSolver(All, Time, TimeStep, SolverData,
				       U, R, X, TransposeFlag));
    break;
  case xfe_TimeSchemeDG2:
    return xf_Error(xf_DG2ApproxSolver(All, Time, TimeStep, SolverData,
				       U, R, X, TransposeFlag));
    break;
  default: 
    return xf_Error(xf_NOT_SUPPORTED); 
    break;
  }

}

/******************************************************************/
//   FUNCTION Definition: xf_DGinTimeStep
static int
xf_DGinTimeStep(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
		real Time, real TimeStep, enum xfe_Bool FineFlag,
		xf_SolverData *SolverData, xf_Vector *Uprev, 
		xf_Vector **Uj, enum xfe_Bool *Redo)
{
/*
PURPOSE:

  Takes a step of a DG in time unsteady solver.  On input, the latest
  state (right node of previous time slab) is stored in Uprev.  On
  output, the entire state in the time slab is returned in the vector
  of state vectors, Uj.
  
  An automatic time-step reduction is performed for robustness in the
  event of a recoverable solver error during either the residual
  calculation, the linear system solution, or the state update.  A
  warning is printed when this reduction, controlled by
  TimeStepDecreaseFactor occurs.

INPUTS:

  All         : all structure
  TimeScheme  : what time scheme to use
  Time        : current time
  TimeStep    : delta t
  FineFlag    : If True, this is a fine space solve (not necessarily to convergence)
  SolverData  : solver data structure (contains CFL)
  Uprev       : state vector at the right node of the previous time slab
  
OUTPUTS: 

  Uj          : pointer to ALL the states on the current time slab (0=left)
  Redo        : on error, this flag is set to True to redo the time step

RETURN: Error code

*/
  int ierr, nIter, iq, nq;
  int iOrder, OrderTime;
  enum xfe_Bool UpdateFlag, Halted, Rewind, Converged;
  enum xfe_Bool Found, LeanFlag;
  enum xfe_Bool LimitFlag;
  enum xfe_PreconditionerType Preconditioner;
  enum xfe_Verbosity Verbosity;
  char StateName[xf_MAXSTRLEN];
  real ResidualTolerance, rtemp, omega;
  xf_Vector *U, *Utemp, *Rtemp;
  xf_Vector **Wj, **Ri;
  xf_JacobianMatrix *R_U;
  
  (*Redo)   = xfe_False;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;

  // set U to Uj[0] (in pointer fashion) for allocation convenience
  U = Uj[0];
  
  // locate preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
				     xfe_PreconditionerName, (int ) xfe_PreconditionerLast, 
				     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;

  // Is the preconditioner memory-lean?
  ierr = xf_Error(xf_PreconditionerLeanCheck(Preconditioner, &LeanFlag));
  if (ierr != xf_OK) return ierr;

  // cannot use lean preconditioner for the approximate factorization
  if (LeanFlag){
    xf_printf("Lean preconditioners not supported for DG in Time.\n");
    return xf_Error(xf_NOT_SUPPORTED);
  }

  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // number of nonlinear iterations
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", &nIter);
  if (ierr != xf_OK) return ierr;

  // locate Jacobian matrix
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
					xfe_True, NULL, &R_U, &Found));
  if (ierr != xf_OK) return ierr;
  

  // pull off time stepping + CFL quantities
  ierr = xf_Error(xf_FindCFLData(All->Param->KeyValue, SolverData));
  if (ierr != xf_OK) return ierr;

  // pull off residual tolerance
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "ResidualTolerance", &ResidualTolerance);
  if (ierr != xf_OK) return ierr;

  // Under-relaxation factor
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "DGTimeUnderRelax", &omega);
  if (ierr != xf_OK) return ierr;

  // locate array of required update vectors, Wj
  ierr = xf_Error(xf_Alloc( (void **) &Wj, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(StateName, "Update_DG_%d", iOrder);
    ierr = xf_Error(xf_FindSimilarVector(All, U, StateName, xfe_True, xfe_True,  
					 NULL, Wj + iOrder, NULL));
    if (ierr != xf_OK) return ierr;
  }

  // locate array of required residual vectors, Ri
  ierr = xf_Error(xf_Alloc( (void **) &Ri, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(StateName, "Residual_DG_%d", iOrder);
    ierr = xf_Error(xf_FindSimilarVector(All,U,StateName, xfe_True, xfe_True,
					 NULL, Ri + iOrder,NULL));
    if (ierr != xf_OK) return ierr;
  }


  // locate temporary state vector, Utemp
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Utemp_DG", xfe_True, xfe_True,
				       NULL, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;
  
  // locate temporary residual vector, Rtemp
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Rtemp_DG", xfe_True, xfe_True,
				       NULL, &Rtemp, NULL));
  if (ierr != xf_OK) return ierr;

  // initialize Uj to Uprev
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    // using project to allow for varying orders
    ierr = xf_Error(xf_ProjectVector(All, Uprev, xfe_False, Uj[iOrder]));
    if (ierr != xf_OK) return ierr;
  }
  
  // initialize Wj to zero
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    ierr = xf_Error(xf_SetZeroVector(Wj[iOrder]));
    if (ierr != xf_OK) return ierr;
  }
  
  /*** Start Newton outer steps ***/

  Converged = xfe_False;
  for (SolverData->iIter=0; SolverData->iIter < nIter; SolverData->iIter++){
      
    if (Halted = xf_CheckUserHalt(NULL)) break;
    
    /* Construct unsteady residual, Ri */
    ierr = xf_Error(xf_DGinTimeUnsteadyResidual(All, TimeScheme, Uj, Uprev, Time, 
						TimeStep, SolverData, Utemp, Rtemp, 
						Ri, &Rewind));
    if (ierr != xf_OK) return ierr;
    if (Rewind) break;
    
    // compute residual norm
    for (iOrder=0, SolverData->ResNorm=0.; iOrder<=OrderTime; iOrder++){
      ierr = xf_Error(xf_VectorNorm(Ri[iOrder], 1, &rtemp));
      if (ierr != xf_OK) return ierr;
      SolverData->ResNorm += rtemp;
    }
    
    // print out residual norm to stdout
    if (Verbosity != xfe_VerbosityLow)
      xf_printf("  DG in Time iteration %d: |R| = %.10E\n", 
		SolverData->iIter, SolverData->ResNorm);
    
    // convergence check
    if (SolverData->ResNorm < ResidualTolerance){
      Converged = xfe_True;
      break;
    }

    // Solve DG linear system using approximate factorization
    ierr = xf_Error(xf_DGApproxSolver(All, TimeScheme, Time, TimeStep, 
				      SolverData, Uj, Ri, Wj, xfe_False));
    if (ierr == xf_REWIND){
      (*Redo) = xfe_True;
      break;
    }
    if (ierr != xf_OK) return ierr;

    // Up to now we calculated the update vectors Wj[j]
    // Update the state vectors, Uj[j], and check them (are they physical?)
    
    // first under-relax the updates (1.0 is default)
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      ierr = xf_Error(xf_VectorMult(Wj[iOrder], omega));
      if (ierr != xf_OK) return ierr;
    } // iOrder
    
    // now do the update
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      ierr = xf_Error(xf_UpdateState(All, Uj[iOrder], Wj[iOrder], &LimitFlag,
				     &UpdateFlag, SolverData));
      if (ierr != xf_OK) return ierr;
	
      if (xf_CheckSolverErrorUnsteady(((UpdateFlag) ? xf_OK : xf_NO_UPDATE), xfe_True,
				      SolverData, 1, &Uprev, OrderTime+1, Uj,
				      &TimeStep, &Rewind))
	return xf_Error(xf_SOLVER_ERROR);
      if (Rewind) break;

    }
    if (Rewind) break;
          
  } // iiter

  if (Rewind) (*Redo) = xfe_True; // need to do time step over  

  // notify if converged
  if ((Converged) && (Verbosity != xfe_VerbosityLow)){
    xf_printf("DG in Time step converged to tolerance.\n");
  }
  
  // notify if not converged
  if ((!Converged) && (!FineFlag)){
    if (Verbosity != xfe_VerbosityLow)
      xf_printf("Warning, DG in Time step not converged at Time = %.10E, TimeStep = %.10E\n",
		Time, TimeStep);
    (*Redo) = xfe_True;
  }

  // notify if redoing
  if (*Redo) xf_printf("Asking to redo this time slab.\n");
  else{
    // write log entry (also print to stdout) only if not redoing time step
    ierr = xf_Error(xf_WriteLogEntry(All, SolverData, Uj[OrderTime]));
    if (ierr != xf_OK) return ierr;
  }

  // memory cleanup
  xf_Release( (void *) Ri);
  xf_Release( (void *) Wj);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ErrEstSpaceTimeAniso
static int
xf_ErrEstSpaceTimeAniso(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
			xf_Vector **Uj, xf_Vector *Uprev, real Time, 
			real TimeStep, xf_Vector *SpaceTimePref,
			xf_Vector *AdaptIndicator)
{
/*

PURPOSE:

  Computes preferred direction of space versus time refinement.
  Currently based on solution jumps.  Result is stored in
  SpaceTimePref, under the following convention:
  
    OLD:
        0 indicates spatial refinement is prefered
        1 indicates temporal refinement is prefered

    NEW:
        value indicates percentage of error attributable to 
	temporal resolution (1 = 100% temporal error)


INPUTS:

  All         : All structure
  TimeScheme  : DG Time scheme
  Uj          : state on current time slab
  Uprev       : final-time (right-hand node) state on previous slab
  Time        : time at start of slab
  TimeStep    : slab length

OUTPUTS: 

  SpaceTimePref : real vector indicating refinement precentage preference
  AdaptIndicator : solution jump adaptive indicator for each element in slab.
		   
RETURN: Error code

*/
  int ierr;
  int egrp, elem, face;
  int dim, sr, d, k, iq, nq;
  int Order, Orderp, OrderL, OrderR, QuadOrder;
  int nqmax = -1;
  int nsize = -1;
  int iOrder, OrderTime;
  int *IParam = NULL;
  enum xfe_Bool QuadChanged;
  enum xfe_Bool MotionOn;
  real *xq, *wq;
  real *xelemL = NULL;
  real *xelemR = NULL;
  real *xelem  = NULL;
  real *xglob  = NULL;
  real *uL = NULL;
  real *uR = NULL;
  real *wn = NULL;
  real *UM = NULL;
  real *RParam = NULL;
  real TimeMid, Num, Den, dubarTime, dubarSpace;
  real *CI;
  real CI1[2] = {0.5, 0.5};
  real CI2[3] = {0.25, 0.5, 0.25};
  xf_IFace IFace;
  xf_BFace BFace;
  xf_Face Face;
  xf_QuadData *QuadDataElem = NULL;
  xf_QuadData *QuadDataFace = NULL;
  xf_BasisData *PhiDataElemL = NULL;
  xf_BasisData *PhiDataElemR = NULL;
  xf_BasisData *PhiDataL = NULL;
  xf_BasisData *PhiDataR = NULL;
  xf_BasisData *PhiData  = NULL;
  xf_BasisData *GeomPhiData = NULL;
  xf_NormalData *NData = NULL;
  xf_BasisTable *PhiTable = NULL;
  xf_MotionData *MD = NULL;
  xf_BC *BC;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  
  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;
  dim    = Mesh->Dim;
  sr     = EqnSet->StateRank;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;
  
  CI = CI1;
  if (TimeScheme == xfe_TimeSchemeDG2) CI = CI2;

  // Sort EqnSet->BCs to match boundary face groups
  ierr = xf_Error(xf_SortEqnSetBCs(Mesh, EqnSet->BCs+0));
  if (ierr != xf_OK) return ierr;
  BC = EqnSet->BCs[0].BC;

  // Set Time to that at 0.5 slab, Time + 0.5*TimeStep
  TimeMid = Time+0.5*TimeStep;
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", TimeMid));
  if (ierr != xf_OK) return ierr;

  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam,  NULL, NULL));
  if (ierr != xf_OK) return ierr;

  // Are we doing mesh motion?
  MotionOn =  ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));
  if (MotionOn){
    ierr = xf_Error(xf_CreateMotionData(All, &MD));
    if (ierr != xf_OK) return ierr;
  }

  // zero out SpaceTimePref
  ierr = xf_Error(xf_SetZeroVector(SpaceTimePref));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_CreateBasisTable(&PhiTable));
  if (ierr != xf_OK) return ierr;

  // need halos in correct place for anisotropy measure
  for (iOrder=0; iOrder<=OrderTime; iOrder++){
    // begin communication of halo data: state
    ierr = xf_Error(xf_HaloExchangeVectorBegin(Uj[iOrder]));
    if (ierr != xf_OK) return ierr;
    // immediately end halo communication
    ierr = xf_Error(xf_HaloExchangeVectorEnd(Uj[iOrder]));
    if (ierr != xf_OK) return ierr;
  }

  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
 
      // order of interpolation of U
      Order  = xf_InterpOrder(Uj[0], egrp, elem);
      Orderp = xf_InterpOrder(Uprev, egrp, elem);
      Order = max(Order, Orderp);

      /*-- compute temporal direction average jump, dubarTime --*/

      // quad points on elem
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, Order, 
				  &QuadDataElem, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      nq = QuadDataElem->nquad;
      xq = QuadDataElem->xquad;
      wq = QuadDataElem->wquad;
      // re-allocate uL, uR if necessary
      if (nq > nqmax){
	nqmax = nq;
	ierr = xf_Error(xf_ReAlloc( (void **) &uL, nqmax*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &uR, nqmax*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &wn, nqmax*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nqmax*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }
      // interpolate Uj[0] at xq
      Order = xf_InterpOrder(Uj[0], egrp, elem);
      ierr = xf_Error(xf_EvalBasis(Uj[0]->Basis[egrp], Order, QuadChanged, 
				   nq, xq, xfb_Phi, &PhiDataElemL));
      if (ierr != xf_OK) return ierr;
      xf_MxM_Set(PhiDataElemL->Phi, Uj[0]->GenArray[egrp].rValue[elem], nq, PhiDataElemL->nn, sr, uL);
      // interpolate Uprev at xq
      Orderp = xf_InterpOrder(Uprev, egrp, elem);
      ierr = xf_Error(xf_EvalBasis(Uprev->Basis[egrp], Orderp, QuadChanged, 
				   nq, xq, xfb_Phi, &PhiDataElemR));
      if (ierr != xf_OK) return ierr;
      xf_MxM_Set(PhiDataElemR->Phi, Uprev->GenArray[egrp].rValue[elem], nq, PhiDataElemR->nn, sr, uR);

      // loop over quad points on elem, form dubarTime = Num/Den
      for (iq=0, Num=Den=0.; iq<nq; iq++){
	for (k=0; k<sr; k++) Num += wq[iq]*fabs(uL[iq*sr+k]-uR[iq*sr+k]);
	Den += wq[iq];
      } // iq
      dubarTime = Num/Den;


      /*-- compute spatial direction average jump, dubarSpace --*/
      Num = Den = 0.;
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
	
	// Face
	Face = Mesh->ElemGroup[egrp].Face[elem][face];
	
	// quadrature order
	if (Face.Group == xf_INTERIORFACE){
	  IFace = Mesh->IFace[Face.Number];
	  OrderL = xf_InterpOrder(Uj[0], IFace.ElemGroupL, IFace.ElemL);
	  OrderR = xf_InterpOrder(Uj[0], IFace.ElemGroupR, IFace.ElemR);
	  QuadOrder = max(OrderL, OrderR);
	}
	else if (Face.Group >= 0) { // Boundary Face
	  QuadOrder = xf_InterpOrder(Uj[0], egrp, elem);
	}
	else{
	  return xf_Error(xf_NOT_SUPPORTED);
	}

	// quad points on face
	ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face, QuadOrder, 
				    &QuadDataFace, &QuadChanged));
	if (ierr != xf_OK) return ierr;
	nq = QuadDataFace->nquad;
	xq = QuadDataFace->xquad;
	wq = QuadDataFace->wquad;

	// reallocate uL, uR if necessary
	if (nq > nqmax){
	  nqmax = nq;
	  ierr = xf_Error(xf_ReAlloc( (void **) &uL, nqmax*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &uR, nqmax*sr, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &wn, nqmax*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nqmax*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}
	// interpolate state inside and outside elem
	if (Face.Group == xf_INTERIORFACE){  /** Interior Face **/
	  // compute basis functions on L and R
	  OrderL = xf_InterpOrder(Uj[0], IFace.ElemGroupL, IFace.ElemL);
	  ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, IFace.ElemGroupL, IFace.ElemL,
						       IFace.FaceL, IFace.OrientL,
						       Uj[0]->Basis[IFace.ElemGroupL],
						       OrderL, QuadChanged,
						       nq, xq, xfb_Phi, &PhiDataL, PhiTable, &xelemL));
	  if (ierr != xf_OK) return ierr;
	  OrderR = xf_InterpOrder(Uj[0], IFace.ElemGroupR, IFace.ElemR);
	  ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, IFace.ElemGroupR, IFace.ElemR,
						       IFace.FaceR, IFace.OrientR,
						       Uj[0]->Basis[IFace.ElemGroupR],
						       OrderR, QuadChanged,
						       nq, xq, xfb_Phi, &PhiDataR, PhiTable, &xelemR));
	  if (ierr != xf_OK) return ierr; 
	  // reallocate UM
	  if (max(PhiDataL->nn,PhiDataR->nn) > nsize){
	    nsize = max(PhiDataL->nn,PhiDataR->nn);
	    ierr = xf_Error(xf_ReAlloc( (void **) &UM, nsize*sr, sizeof(real)));
	    if (ierr != xf_OK) return ierr;
	  }
	  // interpolate state at quad points (use state at half time slab)
	  for (k=0; k<PhiDataL->nn*sr; k++) // L
	    for (iOrder=0, UM[k]=0.; iOrder<=OrderTime; iOrder++)
	      UM[k] += CI[iOrder]*Uj[iOrder]->GenArray[IFace.ElemGroupL].rValue[IFace.ElemL][k];
	  xf_MxM_Set(PhiDataL->Phi, UM, nq, PhiDataL->nn, sr, uL);
	  for (k=0; k<PhiDataR->nn*sr; k++) // R
	    for (iOrder=0, UM[k]=0.; iOrder<=OrderTime; iOrder++)
	      UM[k] += CI[iOrder]*Uj[iOrder]->GenArray[IFace.ElemGroupR].rValue[IFace.ElemR][k];
	  xf_MxM_Set(PhiDataR->Phi, UM, nq, PhiDataR->nn, sr, uR);
	}
	else if (Face.Group >= 0) { /** Boundary Face **/
	  
	  BFace = Mesh->BFaceGroup[Face.Group].BFace[Face.Number];
	  // compute basis functions if quad or basis or order changed
	  Order = xf_InterpOrder(Uj[0], egrp, elem);
	  ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrp, elem, face, BFace.Orient,
						       Uj[0]->Basis[egrp], Order,
						       QuadChanged, nq, xq, xfb_Phi, 
						       &PhiData, PhiTable, &xelem));
	  if (ierr != xf_OK) return ierr;
	  // reallocate UM
	  if (PhiData->nn > nsize){
	    nsize = PhiData->nn;
	    ierr = xf_Error(xf_ReAlloc( (void **) &UM, nsize*sr, sizeof(real)));
	    if (ierr != xf_OK) return ierr;
	  }
	  // interpolate state at quad points (use state at half time slab)
	  for (k=0; k<PhiData->nn*sr; k++)
	    for (iOrder=0, UM[k]=0.; iOrder<=OrderTime; iOrder++)
	      UM[k] += CI[iOrder]*Uj[iOrder]->GenArray[egrp].rValue[elem][k];
	  xf_MxM_Set(PhiData->Phi, UM, nq, PhiData->nn, sr, uL);
	  
	  /* Compute normal(s) at quad points.  If face is straight, only
	     one normal will be computed/returned. */
	  ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, nq, xq, &NData, NULL));
	  if (ierr != xf_OK) return ierr;
	  for (d=0; d<dim; d++)
	    for (iq=0;iq<nq; iq++) 
	      wn[iq*dim+d] = NData->n[iq*dim*(NData->nq!=1)+d];

	  // obtain global coords of quad points
	  ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, xfe_True, 
					  nq, xelem, xglob));
	  if (ierr != xf_OK) return ierr;

	  // obtain transformation map if doing mesh motion
	  if (MotionOn){
	    ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, All->Mesh->Motion, 
					      nq, dim, TimeMid, xglob, MD));
	    if (ierr != xf_OK) return ierr;
	    // transform state + normal to physical
	    if (MotionOn) xf_ModMotionPreEqnCall(nq, dim, sr, MD, uL, wn);
	  }

	  // obtain boundary state and derivative at quad points
	  ierr = xf_Error(xf_EqnSetBCState(EqnSet, BC+Face.Group, IParam, RParam,
					   nq, wn, xglob, &TimeMid, 
					   (MotionOn) ? MD->vg : NULL, uL, uR, NULL));
	  if (ierr != xf_OK) return ierr;

	  // transform state and normal back to reference space
	  if (MotionOn) xf_ModMotionPostEqnCall(nq, dim, sr, MD, uL, wn);

	}
	else{
	  return xf_Error(xf_NOT_SUPPORTED);
	}

	// loop over quad points on elem, form dubarSpace = Num/Den
	for (iq=0; iq<nq; iq++){
	  for (k=0; k<sr; k++) Num += wq[iq]*fabs(uL[iq*sr+k]-uR[iq*sr+k]);
	  Den += wq[iq];
	} // iq
	
      } // face
      dubarSpace = Num/Den;
      
      // set SpaceTimePref based on dubarTime and dubarSpace
      /* xf_printf("egrp=%d, elem=%d, dubarTime = %.6E, dubarSpace = %.6E\n", */
      /*       		egrp, elem, dubarTime, dubarSpace);  */
      
      // Option 1: Binary switch
      //SpaceTimePref->GenArray[egrp].iValue[elem][0] = (int) (dubarSpace < dubarTime);
      // Option 2: Continuous switch, MEPS = safety to avoid dividing by zero
      SpaceTimePref->GenArray[egrp].rValue[elem][0] = dubarTime / (dubarSpace + dubarTime + MEPS);
      if (AdaptIndicator != NULL)
	AdaptIndicator->GenArray[egrp].rValue[elem][0] += (dubarTime + dubarSpace);

    } // elem
  } // egrp
   

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiDataElemL, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(PhiDataElemR, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  /* Destroy Normal Data */
  ierr = xf_Error(xf_DestroyNormalData(NData));
  if (ierr != xf_OK) return ierr;

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataElem));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataFace));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Table */
  ierr = xf_Error(xf_DestroyBasisTable(PhiTable));
  if (ierr != xf_OK) return ierr;

  /* Destroy mesh motion data */
  xf_DestroyMotionData(MD);

  // Release memory
  xf_Release( (void *) uL);
  xf_Release( (void *) uR);
  xf_Release( (void *) wn);
  xf_Release( (void *) UM);
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) xelemL);
  xf_Release( (void *) xelemR);
  xf_Release( (void *) xelem);
  xf_Release( (void *) xglob);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SpaceLocalizeAdd
static int
xf_SpaceLocalizeAdd(xf_Vector *A, xf_Vector *B, real c, 
		    real *OutputError, xf_Vector *ErrIndElem)
{
/*

PURPOSE:

  Spatially localizes an inner product between two vectors (and a
  multiplicative scalar):

       (*OutputError) += c * A^T * B

  The result of localizing this inner product is an element-based
  indicator.  The indicator is signed (no absolute values), and is
  added to ErrIndElem.  If A is NULL, then an absolute-value L1 norm
  indicator is computed, consisting of |c*B|.

INPUTS:

  A           : first vector of inner product
  B           : second vector of inner product
  c           : multiplicative constant

OUTPUTS:

  OutputError : gets incremented with inner product value
  ErrIndElem : gets incremented with *localized* inner product
		   
RETURN: Error code

*/
  int ierr, k;
  int egrp, elem;
  real dp;
  real *EA, *EB;
  xf_GenArray *ga;

  if (A != NULL){
    // take inner product -- checks if vectors are compatible
    ierr = xf_Error(xf_VectorDot(A, B, &dp));
    if (ierr != xf_OK) return ierr;
    
    // add inner product (times c) to OutputError
    if (OutputError != NULL) (*OutputError) += dp*c;
  }
  
  // Localization (note, no absolute value signs at this point)
  for (egrp=0; egrp<B->nArraySelf; egrp++){
    for (elem=0; elem<B->GenArray[egrp].n; elem++){
      EB = B->GenArray[egrp].rValue[elem];
      if (A == NULL){ // unweighted residual, use abs for L1 norm
	ga = B->GenArray+egrp;
	for (k=0,dp=0.; k<((ga->vr==NULL) ? ga->r : ga->vr[elem]); k++)
	  dp += fabs(c*EB[k]);
      }
      else{ // weighted residual, no abs
	EA = A->GenArray[egrp].rValue[elem];
	ga = A->GenArray+egrp;
	for (k=0,dp=0.; k<((ga->vr==NULL) ? ga->r : ga->vr[elem]); k++)
	  dp += c*EA[k]*EB[k];
      }
      ErrIndElem->GenArray[egrp].rValue[elem][0] += dp;
    } // elem
  } // egrp
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DGTimeBuilddAdj
static int
xf_DGTimeBuilddAdj(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
		    xf_Vector **Adji, const char *Root, xf_Vector ***pdAdji)
{
/*

PURPOSE:

  Builds delta adjoint vectors for DG time error estimation. These
  vectors are stored in All with prefix Root.  Spatial and temporal
  order increment is taken into account in this function.

INPUTS:

  All         : all structure
  TimeScheme  : time scheme
  Adji        : adjoints on time slab (OrderTime+1 of them) 

OUTPUTS: 

  (*pdAdji)   : allocated adjoint delta vectors (should be released by caller)
		   
RETURN: Error code

*/
  int ierr, i;
  int OrderTime;
  int SpaceOrderIncrement = 0;
  int TimeOrderIncrement = 0;
  real xt[xf_MAXDGTIMENODE];
  xf_Vector **dAdji = NULL;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;
   
  // time nodes
  for (i=0; i<=OrderTime; i++)
    xt[i] = ((real) (i)) / ((real) (OrderTime));

  // locate array of adjoint delta vectors
  ierr = xf_Error(xf_FindSimilarVectors(All, Adji[0], Root, OrderTime+1,
                                        xfe_True, xfe_True, pdAdji));
  if (ierr != xf_OK) return ierr;
  dAdji = (*pdAdji);

  // dAdji = Adji
  for (i=0; i<OrderTime+1; i++){
    ierr = xf_Error(xf_SetVector(Adji[i], xfe_Set, dAdji[i]));
    if (ierr != xf_OK) return ierr;
  } // i

  // No longer taking out coarse adjoint

/*   // temporal order increment */
/*   ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UErrEstOrderIncrement",  */
/* 				    &TimeOrderIncrement)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // spatial order increment */
/*   ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement",  */
/* 				    &SpaceOrderIncrement)); */
/*   if (ierr != xf_OK) return ierr; */


/*   // only take out projection if incrementing order in both space and time */
/*   if ((SpaceOrderIncrement > 0) && (TimeOrderIncrement > 0)){ */
    
/*     if (TimeOrderIncrement > 0){ */
/*       // project dAdj down temporally */
/*       ierr = xf_Error(xf_DGTimeInterpolateDown(TimeScheme, Adji, OrderTime+1, xt, */
/*                                                TimeOrderIncrement, dAdji)); */
/*       if (ierr != xf_OK) return ierr; */
/*     } */

/*     if (SpaceOrderIncrement > 0){ */
/*       for (i=0; i<OrderTime+1; i++){ */
/*         // project dAdj down spatially */
/*         ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, dAdji[i], NULL, */
/*                                                                xfe_BasisLast, -SpaceOrderIncrement)); */
/*         if (ierr != xf_OK) return ierr; */
/*         // Project dAdj up spatially (info already lost as desired) */
/*         ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, dAdji[i], NULL, */
/*                                                                xfe_BasisLast, SpaceOrderIncrement)); */
/*         if (ierr != xf_OK) return ierr; */
/*       } // i */
/*     } */
    
/*     // now set dAdji = Adji - dAdji */
/*     for (i=0; i<OrderTime+1; i++){ */
/*       // dAdji *= -1 */
/*       ierr = xf_Error(xf_VectorMult(dAdji[i], -1.0)); */
/*       if (ierr != xf_OK) return ierr; */
/*       // dAdji += Adji */
/*       ierr = xf_Error(xf_SetVector(Adji[i], xfe_Add, dAdji[i])); */
/*       if (ierr != xf_OK) return ierr; */
/*     } // i */
/*   } */

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DGTimeLocalizeError
static int
xf_DGTimeLocalizeError(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
                       xf_Vector **Ri, xf_Vector **dAdji, xf_Vector *Utemp,
                       xf_Vector *SpaceTimePref, real *pOutputError, 
                       real *pOutputErrorSpace, real *pOutputErrorTime, 
                       xf_Vector *ErrIndElem, xf_Vector *ErrIndElemS, 
                       xf_Vector *ErrIndElemT)
{
/*

PURPOSE:

  Localizes adjoint-weighted-residual error to
  elemental/spatial/temporal contributions.


INPUTS:

  All         : all structure
  TimeScheme  : time scheme
  Ri          : residuals on time slab (OrderTime+1 of them)
  dAdji       : adjoint deltas on time slab (OrderTime+1 of them) 
  SpaceTimePref : real vector indicating refinement precentage preference
                  (temporal versus spatial)

OUTPUTS: 
  
  (*OutputError)     : scalar output error estimate on time slab.  Pure
                       adjoint-residual inner product -- can be a positive
                       or negative number.  Gets incremented.
  (*OutputErrorSpace): spatial contribution to output error (gets incremented)
  (*OutputErrorTime) : temporal contribution to output error (gets incremented)
  (*ErrIndElem)      : total output error indicator localized to each 
                       element.  Gets incremented.
  (*ErrIndElemS)     : spatial portion of output error indicator localized to each 
                       element.  Gets incremented.
  (*ErrIndElemT)     : temporal portion of output error indicator localized to each 
                       element.  Gets incremented.
		   
RETURN: Error code

*/

  int ierr, i, k, r;
  int egrp, elem;
  int OrderTime;
  int SpaceOrderIncrement = 0;
  int TimeOrderIncrement = 0;
  enum xfe_Bool UseGCL = xfe_False;
  real xt[xf_MAXDGTIMENODE];
  real dp, fac, *ER, *EU;
  real OutputError, OutputErrorSpace, OutputErrorTime;

  OutputError      = 0.;
  OutputErrorSpace = 0.;
  OutputErrorTime  = 0.;
  

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;
  
  // temporal order increment
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UErrEstOrderIncrement", 
				    &TimeOrderIncrement));
  if (ierr != xf_OK) return ierr;

  // spatial order increment
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement", 
				    &SpaceOrderIncrement));
  if (ierr != xf_OK) return ierr;
 
  // time nodes
  for (i=0; i<=OrderTime; i++)
    xt[i] = ((real) (i)) / ((real) (OrderTime));
  
  // dot product yields output error contribution
  if (dAdji != NULL){
    for (i=0; i<OrderTime+1; i++){
      ierr = xf_Error(xf_VectorDot(dAdji[i], Ri[i], &dp));
      if (ierr != xf_OK) return ierr;  
      OutputError -= dp;
    } // i
  }

  // localize output error -> ErrIndElem
  if (ErrIndElem != NULL){
    for (i=0; i<OrderTime+1; i++){
      ierr = xf_Error(xf_SpaceLocalizeAdd(((dAdji == NULL) ? NULL : dAdji[i]), 
					  Ri[i], -1.0, NULL, ErrIndElem));
      if (ierr != xf_OK) return ierr; 
    } // i
  }
  
  if (SpaceTimePref == NULL){ // no anisotropy measure provided, use projection
    
    // Spatial error
    for (i=0; i<OrderTime+1; i++){
      if (dAdji != NULL){
	// interpolate dAdji down in time -> Utemp
	ierr = xf_Error(xf_DGTimeInterpolateDown(TimeScheme, dAdji, 1, xt+i,
						 TimeOrderIncrement, &Utemp));
	if (ierr != xf_OK) return ierr;
      }
      // localize dot product
      ierr = xf_Error(xf_SpaceLocalizeAdd(((dAdji == NULL) ? NULL : Utemp), 
					  Ri[i], -1.0, &OutputErrorSpace, 
					  ErrIndElemS));
      if (ierr != xf_OK) return ierr;  
    } //  i
  
    // Temporal error
    for (i=0; i<OrderTime+1; i++){
      if (dAdji != NULL){
	// Project dAdji down spatially
	ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, dAdji[i], NULL, 
							       xfe_BasisLast, -SpaceOrderIncrement));
	if (ierr != xf_OK) return ierr;
	// Project dAdji up spatially (info already lost as desired)
	ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, dAdji[i], NULL, 
							       xfe_BasisLast, SpaceOrderIncrement));
	if (ierr != xf_OK) return ierr;
      }
      // localize dot product
      ierr = xf_Error(xf_SpaceLocalizeAdd(((dAdji == NULL) ? NULL : dAdji[i]), 
					  Ri[i], -1.0, &OutputErrorTime,
					  ErrIndElemT));
      if (ierr != xf_OK) return ierr;
    }// i

  }
  
  if (SpaceTimePref != NULL){
    //desired anisotropy is given
    for (i=0; i<OrderTime+1; i++){

      // first, spatial fraction of error
      if (ErrIndElemS != NULL){
	for (egrp=0; egrp<All->Mesh->nElemGroup; egrp++){
	  for (elem=0; elem<All->Mesh->ElemGroup[egrp].nElem; elem++){
	    // fac = temporal error fraction, assume 1-fac = spatial
	    fac = SpaceTimePref->GenArray[egrp].rValue[elem][0];
	    ER  = Ri[i]->GenArray[egrp].rValue[elem];
	    EU  = Utemp->GenArray[egrp].rValue[elem];
	    r = ((Utemp->GenArray[egrp].vr == NULL) ? 
		 Utemp->GenArray[egrp].r : Utemp->GenArray[egrp].vr[elem]);
	    for (k=0; k<r; k++) EU[k] = ER[k]*(1.0-fac);
	  } //elem
	} //egrp
	ierr = xf_Error(xf_SpaceLocalizeAdd(((dAdji == NULL) ? NULL : dAdji[i]), 
					    Utemp, -1.0, &OutputErrorSpace, 
					    ErrIndElemS));
	if (ierr != xf_OK) return ierr;
      }

      // second, temporal fraction of error
      if (ErrIndElemT != NULL){
	for (egrp=0; egrp<All->Mesh->nElemGroup; egrp++){
	  for (elem=0; elem<All->Mesh->ElemGroup[egrp].nElem; elem++){
	    // fac = temporal error fraction
	    fac = SpaceTimePref->GenArray[egrp].rValue[elem][0];
	    ER  = Ri[i]->GenArray[egrp].rValue[elem];
	    EU  = Utemp->GenArray[egrp].rValue[elem];
	    r = ((Utemp->GenArray[egrp].vr == NULL) ? 
		 Utemp->GenArray[egrp].r : Utemp->GenArray[egrp].vr[elem]);
	    for (k=0; k<r; k++) EU[k] = ER[k]*fac;
	  } //elem
	} //egrp
	ierr = xf_Error(xf_SpaceLocalizeAdd(((dAdji == NULL) ? NULL : dAdji[i]), 
					    Utemp, -1.0, &OutputErrorTime, 
					    ErrIndElemT));
	if (ierr != xf_OK) return ierr;
      }
    } // i
  }

  if (pOutputError      != NULL) (*pOutputError     ) += OutputError;
  if (pOutputErrorSpace != NULL) (*pOutputErrorSpace) += OutputErrorSpace;
  if (pOutputErrorTime  != NULL) (*pOutputErrorTime ) += OutputErrorTime;


  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ErrEstDGTimeSlab
static int
xf_ErrEstDGTimeSlab(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
		    real Time, real TimeStep, xf_SolverData *SolverData,
		    xf_Vector **Uj, xf_Vector *Uprev, xf_Vector **GCLj,
                    xf_Vector *GCLprev, xf_Vector **Adji, xf_Vector **AdjGCLi,
		    xf_Vector *SpaceTimePref, real *pOutputError, 
		    real *pOutputErrorSpace, real *pOutputErrorTime, 
		    xf_Vector *ErrIndElem, xf_Vector *ErrIndElemS, 
		    xf_Vector *ErrIndElemT)
{
/*

PURPOSE:

  Computes output-based error estimate for one output in a DG-in-time
  unsteady discretization, using the adjoint-weighted residual method.

  GCL contribution to output error estimate is added if GCL vectors
  are provided.

  If Adjoint is not given, unweighted residual error estimate will be
  used (as an adaptive indicator).

INPUTS:

  All         : all structure
  TimeScheme  : time scheme
  Time        : time at start of slab
  TimeStep    : size of time slab
  SolverData  : solver data structure
  Uj          : states on time slab (OrderTime+1 of them)
  Uprev       : state on the previous time slab
  GCLj        : GCL vectors on time slab (OrderTime+1 of them)
  GCLprev     : GCL vector on the previous time slab
  Adji        : adjoints on time slab (OrderTime+1 of them) 
  AdjGCLi     : GCL adjoints on time slab (OrderTime+1 of them) 
  SpaceTimePref : real vector indicating refinement precentage preference
                  (temporal versus spatial)

OUTPUTS: 
  
  (*OutputError)     : scalar output error estimate on time slab.  Pure
                       adjoint-residual inner product -- can be a positive
                       or negative number.  Gets incremented.
  (*OutputErrorSpace): spatial contribution to output error (gets incremented)
  (*OutputErrorTime) : temporal contribution to output error (gets incremented)
  (*ErrIndElem)      : total output error indicator localized to each 
                       element.  Gets incremented.
  (*ErrIndElemS)     : spatial portion of output error indicator localized to each 
                       element.  Gets incremented.
  (*ErrIndElemT)     : temporal portion of output error indicator localized to each 
                       element.  Gets incremented
		   
RETURN: Error code

*/
  int ierr, i;
  int egrp, elem;
  int OrderTime;
  int SpaceOrderIncrement;
  enum xfe_Bool Rewind = xfe_False;
  enum xfe_Bool UseGCL = xfe_False;
  xf_Vector *Rtemp, *Utemp, *GCLtemp = NULL;
  xf_Vector **Ri = NULL;
  xf_Vector **RGCLi = NULL;
  xf_Vector **dAdji = NULL;
  xf_Vector **dAdjGCLi = NULL;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;

  // spatial order increment
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement", 
				    &SpaceOrderIncrement));
  if (ierr != xf_OK) return ierr;

  // are we also doing a GCL-adjoint-weighted GCL-residual?
  // TODO: also use GCL for residual-only error estimate (no adjoint)
  UseGCL = ((GCLj != NULL) && (GCLprev != NULL) && (AdjGCLi != NULL));

  // look for a temporary state vector
  ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], "Utemp_DG_HO", xfe_True, xfe_True, 
				       NULL, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;

  // same with GCL
  if (UseGCL){
    ierr = xf_Error(xf_FindSimilarVector(All, GCLj[0], "GCLtemp_DG_HO", xfe_True, xfe_True, 
                                         NULL, &GCLtemp, NULL));
    if (ierr != xf_OK) return ierr;
  }

  // look for a temporary residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], "Rtemp_DG_HO", xfe_True, xfe_True, 
				       NULL, &Rtemp, NULL));
  if (ierr != xf_OK) return ierr;
  
  // Build array of adjoint delta vectors
  if (Adji != NULL){
    ierr = xf_Error(xf_DGTimeBuilddAdj(All, TimeScheme, Adji, "Adjoint_DG_HO", &dAdji));
    if (ierr != xf_OK) return ierr;
  }
  // same with GCL
  if (AdjGCLi != NULL){
    ierr = xf_Error(xf_DGTimeBuilddAdj(All, TimeScheme, AdjGCLi, "GCLAdjoint_DG_HO", &dAdjGCLi));
    if (ierr != xf_OK) return ierr;
  }

  // account for possible p-dependence of residual
  SolverData->ResidualOrderIncrement = -SpaceOrderIncrement;
  
  // locate array of residual vectors
  ierr = xf_Error(xf_FindSimilarVectors(All, Uj[0], "Residual_DG_HO", OrderTime+1,
					xfe_True, xfe_True, &Ri));
  if (ierr != xf_OK) return ierr;

  // calculate unsteady residual
  ierr = xf_Error(xf_DGinTimeUnsteadyResidual(All, TimeScheme, Uj, Uprev, Time, 
					      TimeStep, SolverData, Utemp, Rtemp, 
					      Ri, &Rewind));
  if (ierr != xf_OK) return ierr;
  if (Rewind) return xf_Error(xf_NOT_SUPPORTED);

  // Also calculate GCL residual
  if (UseGCL){
    ierr = xf_Error(xf_FindSimilarVectors(All, GCLj[0], "GCLResidual_DG_HO", OrderTime+1,
                                          xfe_True, xfe_True, &RGCLi));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DGTimeUnsteadyResidualGCL(All, TimeScheme, GCLj, GCLprev, GCLtemp,
                                                 Time, TimeStep, RGCLi));
    if (ierr != xf_OK) return ierr;
  }    

  // done calculating residuals for error estimation; reset this variable
  SolverData->ResidualOrderIncrement = 0;


  // Localize error (and add to OutputError, etc.)
  ierr = xf_Error(xf_DGTimeLocalizeError(All, TimeScheme, Ri, (Adji==NULL) ? NULL: dAdji,
                                         Utemp, SpaceTimePref,
                                         pOutputError, pOutputErrorSpace, pOutputErrorTime,
                                         ErrIndElem, ErrIndElemS, ErrIndElemT));
  if (ierr != xf_OK) return ierr;
  //same with GCL
  if (UseGCL){
    ierr = xf_Error(xf_DGTimeLocalizeError(All, TimeScheme, RGCLi, (AdjGCLi==NULL) ? NULL : dAdjGCLi,
                                           GCLtemp, SpaceTimePref,
                                           pOutputError, pOutputErrorSpace, pOutputErrorTime,
                                           ErrIndElem, ErrIndElemS, ErrIndElemT));
    if (ierr != xf_OK) return ierr;
  }

  //Reset Time to start of time slab
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) Ri);
  xf_Release( (void *) dAdji);
  xf_Release( (void *) RGCLi);
  xf_Release( (void *) dAdjGCLi);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteTemporalError
static int
xf_WriteTemporalError(int nSlab, int nPsi, xf_Vector **Psi, const char *OutName,
		      real **ErrIndTime, int *sdofTime, const char *fname)
{
/*

PURPOSE:

  Writes temporally-localized error indicator to a text file.

INPUTS:

  nSlab       : number of time slabs
  nPsi        : number of adjoint vectors
  Psi         : all of the adjoint vectors
  OutName     : if Psi is not specified, name to use for writing error
  ErrIndTime  : nSlab x nPsi array of temporally-localized indicators
  sdofTime    : nSlab vector of spatial degrees of freedom per slab
  fname       : file name to write

OUTPUTS: 

  None, fname is written (over-written if exists).  
		   
RETURN: Error code

*/
  int ierr, myRank;
  int iSlab, iAdjoint;
  const char *s = NULL;
  FILE *fid;

  /* Determine myRank */
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // open file for writing (error handled in parallel)
  ierr = xf_Error(xf_fopen(fname, "w", &fid));
  if (ierr != xf_OK) return ierr;

  // only root writes
  if (myRank == 0){
    // write out header
    fprintf(fid, "%% Temporally-localized error estimates\n");
    fprintf(fid, "%% iSlab");
    if (sdofTime != NULL) fprintf(fid, " DOF");
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      s = ((Psi != NULL) ? Psi[iAdjoint]->OutputName : OutName);
      fprintf(fid, " %s", s);
    }
    fprintf(fid, "\n");
    
    // write out output error data
    for (iSlab=0; iSlab<nSlab; iSlab++){
      fprintf(fid, " %d", iSlab);
      if (sdofTime != NULL) fprintf(fid, " %d", sdofTime[iSlab]);
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) 
	fprintf(fid, " %.12E", ErrIndTime[iSlab][iAdjoint]);
      fprintf(fid, "\n");
    } // iSlab
  }
  
  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_ApplyAniso2ErrInd
static int
xf_ApplyAniso2ErrInd(xf_All *All, xf_Vector *ErrIndElem, 
		     xf_Vector *ErrIndElemS, xf_Vector *ErrIndElemT)
{
/*

PURPOSE:

  Sets the spatial (ErrIndElemS) and temporal (ErrIndElemT) output
  error contributions using existing values in ErrIndElemS and
  ErrIndElemT to set the anisotropy, which scales the existing output
  error in ErrIndElem.

INPUTS:

  All             : all structure
  (*ErrIndElem)  : total output error indicator localized to each 
                    element.

OUTPUTS: 
  
  (*ErrIndElemS) : spatial portion of output error indicator localized to each 
                    element. (Also an input for the purpose of anisotropy).
  (*ErrIndElemT) : temporal portion of output error indicator localized to each 
                    element. (Also an input for the purpose of anisotropy).
		   
RETURN: Error code

*/
  int egrp, elem;
  real EI, ES, ET;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  for (egrp=0; egrp<All->Mesh->nElemGroup; egrp++){
    for (elem=0; elem<All->Mesh->ElemGroup[egrp].nElem; elem++){
      EI = ErrIndElem->GenArray[egrp].rValue[elem][0]; // error indicator
      ES = fabs(ErrIndElemS->GenArray[egrp].rValue[elem][0]); // spatial error indicator
      ET = fabs(ErrIndElemT->GenArray[egrp].rValue[elem][0]); // temporal error indicator
      // use MEPS to avoid division by zero when errors are all zero
      ErrIndElemS->GenArray[egrp].rValue[elem][0] = EI*ES/(ES+ET+MEPS); // set new spatial
      ErrIndElemT->GenArray[egrp].rValue[elem][0] = EI*ET/(ES+ET+MEPS); // set new temporal
    } //elem
  } //egrp
  
  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_InitUErrEst
static int
xf_InitUErrEst(xf_All *All, xf_TimeHistData *TimeHistData, 
	       const char *SavePrefix,
	       xf_DGTimeErrEstData *EData,
	       enum xfe_Bool *pUErrEstOn)
{
/*
PURPOSE:

  Initializes unsteady error estimation variables

INPUTS:

  All         : all structure
  TimeHistData: time history data
  SavePrefix  : prefix for saving files
  
OUTPUTS: 

  EData       : unsteady error estimation structure (initialized)
  UErrEstOn   : flag indicating whether err est should proceed

RETURN: Error code

*/
  int ierr;
  int iOrder;
  char Title[xf_MAXSTRLEN];
  xf_Data *D;

  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &EData->Verbosity));
  if (ierr != xf_OK) return ierr;

  /* Determine what we are adapting on */
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptOn", 
                                     xfe_AdaptOnName, (int ) xfe_AdaptOnLast, 
                                     (int *) &EData->AdaptOn));
  if (ierr != xf_OK) return ierr;

  /* Determine entropy on which we are adapting*/
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "AdaptVariableSet", 
				 EData->AdaptVariableSet));
  if (ierr != xf_OK) return ierr;

  // Are we using a Geometric Conservation Law?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &EData->UseGCL));
  if (ierr != xf_OK) return ierr;

  
  /* Should we proceed with error estimation during forward solve? */
  (*pUErrEstOn) = xfe_False;
  // Yes if adapting on heuristic
  if ((EData->AdaptOn == xfe_AdaptOnInterpol) || 
      (EData->AdaptOn == xfe_AdaptOnResidual)){
    (*pUErrEstOn) = xfe_True;
  }
  // Yes if adapting on entropy variable
  if (xf_NotNull(EData->AdaptVariableSet)){
    if (EData->Verbosity != xfe_VerbosityLow)
      xf_printf("Estimating error in unsteady entropy flux integral output for %s\n", 
		EData->AdaptVariableSet);
    (*pUErrEstOn) = xfe_True;
  }
  if (!(*pUErrEstOn)) return xf_OK;

  // # of time slabs
  EData->nSlab = TimeHistData->nTime;
  
  // Set pointer to SavePrefix
  EData->SavePrefix = SavePrefix;
  
  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
				     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
				     (int *) &EData->TimeScheme));
  if (ierr != xf_OK) return ierr;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(EData->TimeScheme, &EData->OrderTime));
  if (ierr != xf_OK) return ierr;

  // What is the temporal order increment for error estimation/adaptation?
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UErrEstOrderIncrement", 
				    &EData->TimeOrderIncrement));
  if (ierr != xf_OK) return ierr;

  // fine time scheme
  EData->TimeSchemeh = EData->TimeScheme;
  if (EData->TimeOrderIncrement == 1){
    if (EData->TimeScheme != xfe_TimeSchemeDG1) return xf_Error(xf_NOT_SUPPORTED);
    EData->TimeSchemeh = xfe_TimeSchemeDG2;
  }
  else if (EData->TimeOrderIncrement != 0) return xf_Error(xf_NOT_SUPPORTED);

  // determine fine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(EData->TimeSchemeh, &EData->OrderTimeh));
  if (ierr != xf_OK) return ierr;

  // Spatial order increment
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement", 
				    &EData->SpaceOrderIncrement));
  if (ierr != xf_OK) return ierr;

  // Name of adaptive indicator
  EData->AdaptOnName = xfe_AdaptOnName[EData->AdaptOn];

  // create spatial error indicator vector
  sprintf(Title, "ErrIndSpaceTot_%s", EData->AdaptOnName);
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  xfe_True,
				&D, &EData->ErrIndElemSTot, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True; // make indicator writeable
  // zero out spatial error indicator (running total, with absolute values)
  ierr = xf_Error(xf_SetZeroVector(EData->ErrIndElemSTot));
  if (ierr != xf_OK) return ierr;
  // will need another error indicator vector, for local storage (not added to dataset)
  sprintf(Title, "ErrIndElem_%s", EData->AdaptOnName);
  ierr = xf_Error(xf_FindSimilarVector(All, EData->ErrIndElemSTot, Title, xfe_False, 
				       xfe_False,  NULL, &EData->ErrIndElem, NULL));
  if (ierr != xf_OK) return ierr;
  // and another
  sprintf(Title, "ErrIndElemS_%s", EData->AdaptOnName);
  ierr = xf_Error(xf_FindSimilarVector(All, EData->ErrIndElemSTot, Title, xfe_False, 
				       xfe_False,  NULL, &EData->ErrIndElemS, NULL));
  if (ierr != xf_OK) return ierr;


  // Allocate reals for storing temporally-localized adaptive indicator
  ierr = xf_Error(xf_Alloc2( (void ***) &EData->ErrIndTime, EData->nSlab, 1, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  // pull off UErrEstAnisoMeasure
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "UErrEstAnisoMeasure", 
				     xfe_SpaceTimeAnisoName, (int) xfe_SpaceTimeAnisoLast, 
				     (int *) &EData->UErrEstAnisoMeasure));
  if (ierr != xf_OK) return ierr;

  if (((EData->AdaptOn == xfe_AdaptOnInterpol) || 
       (EData->AdaptOn == xfe_AdaptOnResidual)) 
      && (EData->UErrEstAnisoMeasure != xfe_SpaceTimeAnisoJump)){
    xf_printf("Warning, AdaptOn = %s, but AnisoMeasure = %s\n", 
	      xfe_AdaptOnName[EData->AdaptOn],
	      xfe_SpaceTimeAnisoName[EData->UErrEstAnisoMeasure]);
    xf_printf("Setting AnisoMeasure = Jump.\n");
    EData->UErrEstAnisoMeasure = xfe_SpaceTimeAnisoJump;
  }

  // Allocate reals for storing spatial dof
  ierr = xf_Error(xf_Alloc( (void **) &EData->sdofTime, EData->nSlab, sizeof(int)));
  if (ierr != xf_OK) return ierr;


  // create spatial vs. temporal preference vector (not added to dataset)
  sprintf(Title, "SpaceTimePref_State");
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 
				0, NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, 
				xfe_False, NULL, &EData->SpaceTimePref, NULL));
  if (ierr != xf_OK) return ierr;
  // will need another temp error indicator vector
  sprintf(Title, "ErrIndElemT_%s", EData->AdaptOnName);
  ierr = xf_Error(xf_FindSimilarVector(All, EData->ErrIndElemSTot, Title, xfe_False, 
				       xfe_False,  NULL, &EData->ErrIndElemT, NULL));
  if (ierr != xf_OK) return ierr;

  // Allocate vector of state pointers, set to NULL
  ierr = xf_Error(xf_Alloc( (void **) &EData->SUj, EData->OrderTimeh+1, 
			    sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= EData->OrderTimeh; iOrder++)
    EData->SUj[iOrder] = NULL;
  EData->SUprev = NULL;

  // Allocate vector of extra state pointers, set to NULL
  ierr = xf_Error(xf_Alloc( (void **) &EData->SVj, EData->OrderTimeh+1, 
			    sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= EData->OrderTimeh; iOrder++)
    EData->SVj[iOrder] = NULL;
  EData->SVprev = NULL;

  // set output error to zero
  EData->OutputError      = 0.0;
  EData->OutputErrorSpace = 0.0;
  EData->OutputErrorTime  = 0.0;
 

  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_UpdateUErrEst
static int
xf_UpdateUErrEst(xf_All *All, int iSlab, real Time, real TimeStep,
		 xf_Vector **Uj, xf_Vector *Uprev, xf_Vector **GCLj,
                 xf_Vector *GCLprev, enum xfe_Bool SpatialDyno, 
                 xf_SolverData *SolverData, xf_DGTimeErrEstData *EData)
{
/*
PURPOSE:

  Updates unsteady error estimate.

INPUTS:

  All         : all structure
  iSlab       : slab number
  Time        : current time
  TimeStep    : current time step
  Uj          : forward state
  Uprev       : right-node state from previous time slab
  GCLj        : GCL state
  GCLprev     : right-node GCl state from previous time slab
  SpatialDyno : if True, indicates orders on current slab changed
                from orders on previous slab
  SolverData  : solver data structure
  EData       : unsteady error estimation structure
  
OUTPUTS: 

  None

RETURN: Error code
*/
  int ierr;
  int iOrder;
  int UEIWriteInterval;
  enum xfe_Bool Redo;
  char Title[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  xf_Vector *rvec = NULL;
  xf_Vector *GCL = NULL;
  xf_KeyValue KeyValueOrig;


  // zero quantities out
  EData->ErrIndTime[iSlab][0] = 0.0;
  ierr = xf_Error(xf_SetZeroVector(EData->ErrIndElem ));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetZeroVector(EData->ErrIndElemS));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetZeroVector(EData->ErrIndElemT));
  if (ierr != xf_OK) return ierr;
  
   
  // Estimate anisotropy (space vs. time) of error
  if (EData->UErrEstAnisoMeasure == xfe_SpaceTimeAnisoJump){
    
    rvec = ((EData->AdaptOn == xfe_AdaptOnInterpol) ? EData->ErrIndElemS : NULL);
    ierr = xf_Error(xf_ErrEstSpaceTimeAniso(All, EData->TimeScheme, Uj, Uprev, Time, 
					    TimeStep, EData->SpaceTimePref, rvec));
    if (ierr != xf_OK) return ierr;
    if (EData->AdaptOn == xfe_AdaptOnInterpol){  
      ierr = xf_Error(xf_SetVector(EData->ErrIndElemS, xfe_Set, EData->ErrIndElemT));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorVectorMult(EData->SpaceTimePref, EData->ErrIndElemT));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorMult(EData->SpaceTimePref, -1.0));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorAdd(EData->SpaceTimePref, 1.0));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorVectorMult(EData->SpaceTimePref, EData->ErrIndElemS));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  // output and residual-based error estimate requires evaluation of residuals
  // note, output here refers to entropy-adjoint adaptation
  if ((EData->AdaptOn == xfe_AdaptOnResidual) ||
      (EData->AdaptOn == xfe_AdaptOnOutput)) {

    // destroy SUj if spatial dyno
    if ((EData->SUj != NULL) && (SpatialDyno)){
      for (iOrder = 0; iOrder <= EData->OrderTimeh; iOrder++){
        ierr = xf_Error(xf_DestroyVector(EData->SUj[iOrder], xfe_True));
        if (ierr != xf_OK) return ierr;
        EData->SUj[iOrder] = NULL;
      }
    }

    // allocate SUj if necessary
    for (iOrder = 0; iOrder <= EData->OrderTimeh; iOrder++){
      if (EData->SUj[iOrder] == NULL){
	sprintf(Title, "HOState_%d", iOrder);
	ierr = xf_Error(xf_FindSimilarVectorHO(All, Uj[0], Title, xfe_True, 
					       xfe_False, &EData->SpaceOrderIncrement,
					       NULL, NULL,&EData->SUj[iOrder]));
	if (ierr != xf_OK) return ierr;
      }
    } // iOrder

    // allocate SUprev if necessary
    if (EData->SUprev == NULL){
      ierr = xf_Error(xf_FindSimilarVectorHO(All, Uprev, "HOState_Prev", xfe_True, 
					     xfe_False, &EData->SpaceOrderIncrement, 
					     NULL, NULL, &EData->SUprev));
      if (ierr != xf_OK) return ierr;
    }

    // project Uj to higher order -> SUj
    if (EData->SpaceOrderIncrement > 0){
      for (iOrder = 0; iOrder <= EData->OrderTime; iOrder++){
	ierr = xf_Error(xf_ProjectVector(All, Uj[iOrder], xfe_False, 
					 EData->SUj[iOrder]));
	if (ierr != xf_OK) return ierr;
      }
      // also project Uprev -> SUprev
      ierr = xf_Error(xf_ProjectVector(All, Uprev, xfe_False, 
				       EData->SUprev));
      if (ierr != xf_OK) return ierr;

      // with GCL, project spatially in place
      if (EData->UseGCL){
        // GCLj
        for (iOrder = 0; iOrder <= EData->OrderTime; iOrder++){
          ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLj[iOrder], NULL, 
                                                                 xfe_BasisLast, EData->SpaceOrderIncrement));
          if (ierr != xf_OK) return ierr;
        }
        // GCL (note, this vector persists, so here we just make it look like GCLj[0])
        ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL)); // vector must exist
        if(ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, GCL, GCLj[0]));
        if (ierr != xf_OK) return ierr;

        // GCLprev
        ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLprev, NULL, 
                                                               xfe_BasisLast, EData->SpaceOrderIncrement));
        if (ierr != xf_OK) return ierr;
      }

    }
      
    // inject into OrderTimeh
    if (EData->OrderTimeh != EData->OrderTime){
      ierr = xf_Error(xf_DGTimeInject(EData->OrderTime, EData->OrderTimeh, 
				      EData->SUj));
      if (ierr != xf_OK) return ierr;

      if (EData->UseGCL){   // same with GCL

        // make sure we have an extra vector (or more)
        for (iOrder = EData->OrderTime+1; iOrder <= EData->OrderTimeh; iOrder++){
          // note, these will be at appropriate (high) spatial order because using SUj[0] here
          ierr = xf_Error(xf_InitMeshMotionGCLVector(All, EData->SUj[0], iOrder, xfe_True, GCLj + iOrder));
          if (ierr != xf_OK) return ierr;
        }

        // inject GCL vectors to high temporal order
        ierr = xf_Error(xf_DGTimeInject(EData->OrderTime, EData->OrderTimeh, GCLj));
        if (ierr != xf_OK) return ierr;
      }
    }
      
    if (EData->AdaptOn == xfe_AdaptOnResidual){
      // residual error estimate
      ierr = xf_Error(xf_ErrEstDGTimeSlab(All, EData->TimeSchemeh, Time, TimeStep, 
                                          SolverData, EData->SUj, EData->SUprev, 
                                          GCLj, GCLprev, NULL, NULL, 
                                          EData->SpaceTimePref, NULL, NULL, NULL, 
                                          NULL, EData->ErrIndElemS, EData->ErrIndElemT));
      if (ierr != xf_OK) return ierr;
    }
    else{
      // entropy adjoint output error calculation

      // Allocate SVj and SVprev if necessary (not added to All)
      if (EData->SVprev == NULL){
	// first SVprev
	ierr = xf_Error(xf_FindSimilarVectorHO(All, Uprev, "EvarHOState_Prev", xfe_True, 
					       xfe_False, &EData->SpaceOrderIncrement, 
					       NULL, NULL, &EData->SVprev));
	if (ierr != xf_OK) return ierr;
	// next SVj
	for (iOrder = 0; iOrder <= EData->OrderTimeh; iOrder++){
	  sprintf(Title, "EvarHOState_%d", iOrder);
	  ierr = xf_Error(xf_FindSimilarVectorHO(All, Uj[0], Title, xfe_True, 
						 xfe_False, &EData->SpaceOrderIncrement,
						 NULL, NULL,&EData->SVj[iOrder]));
	  if (ierr != xf_OK) return ierr;
	} // iOrder
      }

      // copy SUj to SVj
      for (iOrder = 0; iOrder <= EData->OrderTimeh; iOrder++){
	ierr = xf_Error(xf_SetVector(EData->SUj[iOrder], xfe_Set, EData->SVj[iOrder]));
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_SetVector(EData->SUprev, xfe_Set, EData->SVprev));
      if (ierr != xf_OK) return ierr;

      // Set desired solver flags for the fine space solve
      ierr = xf_Error(xf_PreFineSpaceSolve(All, &KeyValueOrig));
      if (ierr != xf_OK) return ierr;
      // Perform DGinTimeStep (fine space iterations)
      ierr = xf_Error(xf_DGinTimeStep(All, EData->TimeSchemeh, Time, TimeStep, xfe_True,
				      SolverData, EData->SVprev, EData->SVj, &Redo));
      if (ierr != xf_OK) return ierr;
      // Handle Redo (eventually include sub-slabs, and just increase here)
      if (Redo) return xf_Error(xf_NOT_SUPPORTED); // for now
      // Reset solver flags to original
      ierr = xf_Error(xf_PostFineSpaceSolve(All, &KeyValueOrig));
      if (ierr != xf_OK) return ierr;
	
      // Convert SVj to entropy variables
      for (iOrder = 0; iOrder <= EData->OrderTimeh; iOrder++){
	ierr = xf_Error(xf_ChangeVariableSet(All, EData->AdaptVariableSet, 
					     EData->SVj[iOrder]));
        if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_ChangeVariableSet(All, EData->AdaptVariableSet, EData->SVprev));
      if (ierr != xf_OK) return ierr;

      // call xf_ErrEstDGTimeSlab with SVj as adjoint
      rvec = ((EData->UErrEstAnisoMeasure == xfe_SpaceTimeAnisoProj) ? NULL : 
	      EData->SpaceTimePref);
      ierr = xf_Error(xf_ErrEstDGTimeSlab(All, EData->TimeSchemeh, Time, TimeStep, 
					  SolverData, EData->SUj, EData->SUprev, NULL, NULL,
					  EData->SVj, NULL, rvec, &EData->OutputError, 
					  &EData->OutputErrorSpace, 
					  &EData->OutputErrorTime,
					  EData->ErrIndElem, EData->ErrIndElemS,
					  EData->ErrIndElemT));
      if (ierr != xf_OK) return ierr;

      // modify ErrIndElemS and ErrIndElemT to incorporate output error in ErrIndElem
      ierr = xf_Error(xf_ApplyAniso2ErrInd(All, EData->ErrIndElem, EData->ErrIndElemS, 
					   EData->ErrIndElemT));
      if (ierr != xf_OK) return ierr;
      
    }

    
    /* destroy SUprev if spatial dyno (we're doing this after error
       estimation, in preparation for the next time slab) */
    if ((EData->SUprev != NULL) && (SpatialDyno)){
      ierr = xf_Error(xf_DestroyVector(EData->SUprev, xfe_True));
      if (ierr != xf_OK) return ierr;
      EData->SUprev = NULL;
    }


    // undo changes to GCL vector (no loss since we injected in space and time)
    if (EData->UseGCL){ 

      // reverse temporal injection
      ierr = xf_Error(xf_DGTimeInject(EData->OrderTimeh, EData->OrderTime, GCLj));
      if (ierr != xf_OK) return ierr;

      // spatial projection

      // GCLj
      for (iOrder = 0; iOrder <= EData->OrderTime; iOrder++){
        ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLj[iOrder], NULL, 
                                                               xfe_BasisLast, -EData->SpaceOrderIncrement));
        if (ierr != xf_OK) return ierr;
      }
      // GCL (note, this vector persists, so here we just make it look like GCLj[0])
      ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, GCL, GCLj[0]));
      if (ierr != xf_OK) return ierr;
      
      // GCLprev
      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLprev, NULL, 
                                                             xfe_BasisLast, -EData->SpaceOrderIncrement));
      if (ierr != xf_OK) return ierr;

    }

  }

  // write out error indicators if at requested interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UErrEstIndicatorWriteInterval", 
				    &UEIWriteInterval));
  if (ierr != xf_OK) return ierr;
  if ((UEIWriteInterval > 0) && (iSlab % UEIWriteInterval) == 0){
    sprintf(OutputFile, "%s_ErrS_%s%d.data\0", EData->SavePrefix, EData->AdaptOnName, iSlab);
    ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ErrIndicatorSpace", 
					EData->ErrIndElemS, OutputFile));
    if (ierr != xf_OK) return ierr;
    sprintf(OutputFile, "%s_ErrT_%s%d.data\0", EData->SavePrefix, EData->AdaptOnName, iSlab);
    ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ErrIndicatorTime", 
					EData->ErrIndElemT, OutputFile));
    if (ierr != xf_OK) return ierr;
  }


  // add ErrIndElemS to running total
  ierr = xf_Error(xf_VectorAbs(EData->ErrIndElemS)); // always use conservative measure
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetVector(EData->ErrIndElemS, xfe_Add, 
			       EData->ErrIndElemSTot));
  if (ierr != xf_OK) return ierr;

  // add ErrIndElemT to running total
  ierr = xf_Error(xf_VectorAbs(EData->ErrIndElemT));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VectorNorm(EData->ErrIndElemT, 0, EData->ErrIndTime[iSlab]));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
// FUNCTION Definition: xf_FinishUErrEst
static int
xf_FinishUErrEst(xf_All *All, xf_DGTimeErrEstData *EData)
{
/*
PURPOSE:

  Finalizes unsteady error estimation.  Releases memory.

INPUTS:

  All         : all structure
  EData       : unsteady error estimation structure
  
OUTPUTS: 

  None

RETURN: Error code

*/
  int ierr;
  int iOrder, iSlab;
  char OutputFile[xf_MAXSTRLEN];
  real sumSpace, sumTime;

  // Write out total (summed) spatially-localized error indicators
  sprintf(OutputFile, "%s_ErrTot_%s.data\0", EData->SavePrefix, EData->AdaptOnName);
  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ErrIndicatorSpace", 
				      EData->ErrIndElemSTot, OutputFile));
  if (ierr != xf_OK) return ierr;

  // Write out table of temporally-localized error indicators
  sprintf(OutputFile, "%s_TemporalError.txt\0", EData->SavePrefix);
  ierr = xf_Error(xf_WriteTemporalError(EData->nSlab, 1, NULL, EData->AdaptOnName, 
					EData->ErrIndTime, EData->sdofTime, OutputFile));
  if (ierr != xf_OK) return ierr;

  // store calculated error estimates and print out values
  if (EData->Verbosity != xfe_VerbosityLow){
    // Calculate for printing
    ierr = xf_Error(xf_VectorNorm(EData->ErrIndElemSTot, 1, &sumSpace));
    if (ierr != xf_OK) return ierr;
    for (iSlab=0,sumTime=0.; iSlab<EData->nSlab; iSlab++) 
      sumTime += EData->ErrIndTime[iSlab][0];
    //xf_printf("AdaptOn = %s  sumSpace = %.10E  sumTime = %.10E\n", 
    //EData->AdaptOnName, sumSpace, sumTime);
    xf_printf("AdaptOn = %s\n", EData->AdaptOnName);
    xf_printf("ErrEst = %.10E\n", EData->OutputError);
    xf_printf("sumSpace = %.10E\n", sumSpace);
    xf_printf("sumTime = %.10E\n", sumTime);

    xf_printf("AdaptVariableSet = %s\n", EData->AdaptVariableSet);
  }
   
  // destroy ErrIndElemS (not ErrIndElemSTot, which is part of All->DataSet)
  ierr = xf_Error(xf_DestroyVector(EData->ErrIndElemS, xfe_True));
  if (ierr != xf_OK) return ierr;
  // destroy ErrIndElem
  ierr = xf_Error(xf_DestroyVector(EData->ErrIndElem, xfe_True));
  if (ierr != xf_OK) return ierr;
  // destroy ErrIndElemT
  ierr = xf_Error(xf_DestroyVector(EData->ErrIndElemT, xfe_True));
  if (ierr != xf_OK) return ierr;

  // destroy SpaceTimePref, which is also not part of All->DataSet
  ierr = xf_Error(xf_DestroyVector(EData->SpaceTimePref, xfe_True));
  if (ierr != xf_OK) return ierr;

  // destroy SUprev and SVprev
  ierr = xf_Error(xf_DestroyVector(EData->SUprev, xfe_True));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyVector(EData->SVprev, xfe_True));
  if (ierr != xf_OK) return ierr;

  // destroy SUj and SVj
  for (iOrder = 0; iOrder <= EData->OrderTimeh; iOrder++){
    if (EData->SUj != NULL){
      ierr = xf_Error(xf_DestroyVector(EData->SUj[iOrder], xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    if (EData->SVj != NULL){
      ierr = xf_Error(xf_DestroyVector(EData->SVj[iOrder], xfe_True));
      if (ierr != xf_OK) return ierr;
    }
  }

  // release memory
  xf_Release( (void *) EData->SUj);
  xf_Release( (void *) EData->SVj);
  xf_Release( (void  *) EData->sdofTime);
  xf_Release2((void **) EData->ErrIndTime);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CheckSpatialDyno
static int
xf_CheckSpatialDyno(xf_All *All, const char *SavePrefix, int iSlab, 
		    int nVector, xf_Vector **Uj, enum xfe_Bool *pSpatialDyno)
{
  int ierr, j;
  int OrderInc;
  enum xfe_Bool DynamicSpatialRef;
  enum xfe_AdaptMechType AdaptMechanics;
  char Title[xf_MAXSTRLEN];
  xf_DataSet *DataSet = NULL;
  xf_Vector *VOrder = NULL;
  xf_Data *D;

  
  // initially assume no dynamic spatial refinement
  (*pSpatialDyno) = xfe_False;

  // check to see if dynamic spatial refinement is requested
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "DynamicSpatialRef", &DynamicSpatialRef));
  if (ierr != xf_OK) return ierr;

  if (!DynamicSpatialRef) return xf_OK;

  // Determine dynamic spatial refinement order increment
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "DynamicSpatialRefOrderInc", &OrderInc));
  if (ierr != xf_OK) return ierr;
  
  // what kind of refinement/adaptation are we talking about?
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "AdaptMechanics", 
                                     xfe_AdaptMechName, (int ) xfe_AdaptMechLast, 
                                     (int *) &AdaptMechanics));
  if (ierr != xf_OK) return ierr;

  // currently only support order changes
  if (AdaptMechanics != xfe_AdaptMechOrderRef) return xf_Error(xf_NOT_SUPPORTED);

  // try projecting using VOrder in file ... if it exists
  sprintf(Title, "%s_VOrder%d.data\0", SavePrefix, iSlab);
  ierr = xf_ProjectVectors_VOrderFile(All->Mesh, All->DataSet, nVector, Uj, Title);
  if (ierr == xf_NOT_FOUND){ // VOrder does not exist, that's fine
    if (OrderInc == 0) return xf_OK;
  }
  else if (ierr != xf_OK) return xf_Error(ierr);

  // if got to here, that means we are doing dynamic spatial refinement
  (*pSpatialDyno) = xfe_True;

  // increase orders of vectors by one: useful for debugging error estimation
  if (OrderInc != 0){
    for (j=0; j<nVector; j++){
      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Uj[j], 
                                                             NULL, xfe_BasisLast, OrderInc));
      if (ierr != xf_OK) return ierr;
    }
  }

  return xf_OK;
}
 


/******************************************************************/
//   FUNCTION Definition: xf_ApplyTimeScheme_DG
static int
xf_ApplyTimeScheme_DG(xf_All *All, const char *SavePrefix,
		      enum xfe_Bool RestartFlag, xf_Vector **pU0,
		      xf_TimeHistData *TimeHistData,
		      xf_SolverData **pSolverData)
{
  int ierr, i, iStep, nTimeStep, iSlab;
  int WriteInterval, OrderTime, iOrder;
  int Redid = 0, RepartEvery, myRank, nProc;
  enum xfe_TimeSchemeType TimeScheme;
  enum xfe_Bool Redo = xfe_False, ProcViz;
  enum xfe_Bool UErrEstOn = xfe_False;
  enum xfe_Bool SpatialDyno = xfe_False;
  enum xfe_Bool UseGCL;
  char StateName[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  char PreHeader[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  real Time, TimeStep, CFLStart, J;
  real DecreaseFactor;
  xf_Vector **Uj, *U0;   // vector of all state vectors for one time slab
  xf_Vector **GCLj = NULL;
  xf_Vector *Uprev; // previous time slab end node state (just a pointer)
  xf_Vector *GCLprev = NULL;
  xf_Vector *UExtra = NULL; // extra state vector for temp/prev storage
  xf_Vector *GCLExtra = NULL; 
  xf_Vector *GCL = NULL;
  xf_Data *D = NULL;
  xf_DataSet *DataSet = NULL;
  xf_SolverData *SolverData = NULL;
  xf_DGTimeErrEstData ErrEstData;
  U0 = (*pU0);
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  SolverData = (*pSolverData);
  
  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
				     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
				     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;
  
  // need TimeHistData
  if (TimeHistData == NULL) return xf_Error(xf_INPUT_ERROR);

  // TimeHistData must be filled in
  nTimeStep = TimeHistData->nTime;
  if (nTimeStep <  0) return xf_Error(xf_INPUT_ERROR);
  if (nTimeStep == 0) return xf_OK; // nothing to do

  // Unsteady write interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 
				    &WriteInterval));
  if (ierr != xf_OK) return ierr;
  
  // Repartition interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "DD_RepartEveryNdt", 
                                    &RepartEvery));
  if (ierr != xf_OK) return ierr;
  if (RepartEvery <= 0) //no repartitioning
    RepartEvery = nTimeStep+2;
  
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ProcViz", 
                                     &ProcViz));
  if (ierr != xf_OK) return ierr;

   // determine linear residual decrease factor
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "MinLinResDecreaseFactor", 
				     &DecreaseFactor));
  if (ierr != xf_OK) return ierr;
  SolverData->LinResTol = DecreaseFactor;

  // Initial CFL for artificial time stepping
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFLStart));
  if (ierr != xf_OK) return ierr;

  // Determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;

  // allocate space for state vector pointers
  ierr = xf_Error(xf_Alloc( (void **) &Uj, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;

  // Uprev vector = end-node of previous time slab = U0 initially
  Uprev = U0; // no new space allocated

  // locate required state vectors; create if necessary
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(StateName, "State_%d", iOrder);
    ierr = xf_Error(xf_FindSimilarVector(All, U0, StateName, xfe_True, xfe_True,  
					 NULL, Uj + iOrder, NULL));
    if (ierr != xf_OK) return ierr;
    Uj[iOrder]->TimeIndex = iOrder;
    Uj[iOrder]->SolverRole = xfe_SolverRolePrimalState;
  }

  // Allocate/find similar variables for GCL
  if (UseGCL){

    // default GCL vector (contains the IC)
    ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U0, -1, xfe_False, &GCL));
    if (ierr != xf_OK) return ierr;
    
    // using max dg time node in case have to inject to high temporal order for errest
    ierr = xf_Error(xf_Alloc( (void **) &GCLj, xf_MAXDGTIMENODE, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;

    // Set GCLprev = GCLExtra = initial GCL vector (note that GCLprev is just a pointer)
    ierr = xf_Error(xf_FindSimilarVector(All, GCL, "GCLExtra", xfe_True, xfe_False,  
                                         NULL, &GCLExtra, NULL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_SetVector(GCL, xfe_Set, GCLExtra));
    if (ierr != xf_OK) return ierr;
    GCLprev = GCLExtra;

    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U0, iOrder, xfe_True, GCLj + iOrder));
      if (ierr != xf_OK) return ierr;
      GCLj[iOrder]->TimeIndex = iOrder;
    }

  }


  if (SavePrefix != NULL){ // only if requesting saves
    
    /* Create a dataset for writing the initial condition */
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;
    // add state to dataset
    Uprev->TimeIndex = 0; // so we can use as a restart .data file
    ierr = xf_Error(xf_DataSetAdd(DataSet, "State", xfe_Vector,
				  xfe_True, (void *) Uprev, NULL));
    if (ierr != xf_OK) return ierr;
    // add GCL vector to dataset
    if (UseGCL){
      ierr = xf_Error(xf_DataSetAdd(DataSet, GCLVectorTitle, xfe_Vector,
				    xfe_True, (void *) GCL, NULL));
      if (ierr != xf_OK) return ierr;
    }

    /* Write out iTime == 0 vector (initial condition) */
    sprintf(OutputFile, "%s_U%d.data\0", SavePrefix, 0);
    ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, OutputFile));
    if (ierr != xf_OK) return ierr;
  
    // Destroy the DataSet that we just used, but not the vectors
    for (D=DataSet->Head; D != NULL; D=D->Next) D->Data = NULL;
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;

  }

  // set time index of Uprev for later use in repartitioning
  Uprev->TimeIndex = -1;

  
  /* Create a dataset for writing the state at each WriteInterval */
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;
  
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(Title, "State_%d", iOrder);
    ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_Vector, xfe_True, 
                                  (void *) Uj[iOrder],NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  // also add GCL vectors to DataSet for writing if using the GCL
  if (UseGCL){
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      sprintf(Title, "GCL_%d", iOrder);
      ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_Vector, xfe_True, 
                                    (void *) GCLj[iOrder], NULL));
      if (ierr != xf_OK) return ierr;
    } // iOrder 
  }
  
  // Are we doing unsteady error estimation?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UErrEstOn", &UErrEstOn));
  if (ierr != xf_OK) return ierr;


  /* Initialize error estimation variables */
  if (UErrEstOn){
    ierr = xf_Error(xf_InitUErrEst(All, TimeHistData, SavePrefix, &ErrEstData,
                                   &UErrEstOn));
    if (ierr != xf_OK) return ierr;
  }
  

  // Begin loop over time slabs
  for (iSlab=0; iSlab<nTimeStep; iSlab++){
    //check if we are repartitioning the mesh
    if ((iSlab+1)%RepartEvery == 0 && nProc > 1 && iSlab > 0) {
      xf_printf("Rebalancing CPU-work load...\n");fflush(stdout);
      if (UErrEstOn){
        ierr = xf_Error(xf_UnParallelizeDGTimeErrEstData(All, &ErrEstData));
        if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_DestroySolverData(SolverData));
      if (ierr != xf_OK) return ierr;
      //rebalance the cpu work load
      ierr = xf_Error(xf_CPULoadBalance(All, Uj, OrderTime+1));
      if (ierr != xf_OK) return ierr;
      
      //reset pointers in DataSet since the data was reparallelized in xf_CPULoadBalance
      DataSet->Head = NULL;
      DataSet->Tail = NULL;
      
      for (iOrder = 0; iOrder <= OrderTime; iOrder++){
        sprintf(Title, "State_%d", iOrder);
        ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_Vector, xfe_True, (void *) Uj[iOrder], NULL));
        if (ierr != xf_OK) return ierr;
      }
      //find state from previous timestep (TimeIndex = -1)
      ierr = xf_Error(xf_FindPrimalState(All->DataSet, -1, &D, NULL));
      if (ierr != xf_OK) return ierr;
      U0 = (xf_Vector *)D->Data;
      (*pU0) = (xf_Vector *)D->Data;
      ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
      if (ierr != xf_OK) return ierr;
      
      if (UseGCL){
        //reset the pointer to the GCL vector
        ierr = xf_Error(xf_InitMeshMotionGCLVector(All, Uj[0], -1, xfe_False, &GCL));
        if (ierr != xf_OK) return ierr;
        
        //Destroy old GCLExtra
        ierr = xf_Error(xf_DestroyVector(GCLExtra, xfe_True));
        if (ierr != xf_OK) return ierr;
        
        // Set GCLprev = GCLExtra = initial GCL vector (note that GCLprev is just a pointer)
        ierr = xf_Error(xf_FindSimilarVector(All, GCL, "GCLExtra", xfe_True, xfe_False,  
                                             NULL, &GCLExtra, NULL));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_SetVector(GCL, xfe_Set, GCLExtra));
        if (ierr != xf_OK) return ierr;
        GCLprev = GCLExtra;
        
        //reset pointers and store in DataSet for writing out .data
        for (iOrder = 0; iOrder <= OrderTime; iOrder++){
          ierr = xf_Error(xf_InitMeshMotionGCLVector(All, Uj[0], iOrder, xfe_False, GCLj + iOrder));
          if (ierr != xf_OK) return ierr;
          
          sprintf(Title, "GCL_%d", iOrder);
          ierr = xf_Error(xf_DataSetAdd(DataSet, Title, xfe_Vector, xfe_True, 
                                        (void *) GCLj[iOrder], NULL));
          if (ierr != xf_OK) return ierr;
        }
      }
      
      if (ProcViz) {
        ierr = xf_Error(xf_FindDataByTitle(All->DataSet, "ProcID", xfe_Vector, &D));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_DataSetAdd(DataSet, D->Title, xfe_Vector, 
                                      xfe_True, D->Data, NULL));
        if (ierr != xf_OK) return ierr;
      }
      
      /* Destroy UExtra if it exists and spatial refinement changed */
      ierr = xf_Error(xf_DestroyVector(UExtra, xfe_True));
      if (ierr != xf_OK) return ierr;
      UExtra = NULL;
      
      ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], "UExtra", xfe_True, xfe_False,  
                                           NULL, &UExtra, NULL));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_SetVector(Uj[OrderTime], xfe_Set, UExtra));
      if (ierr != xf_OK) return ierr;
      Uprev = UExtra;
      
      if (UErrEstOn){
        ierr = xf_Error(xf_ParallelizeDGTimeErrEstData(All, &ErrEstData));
        if (ierr != xf_OK) return ierr;
      }
      
      xf_printf("done.\n");fflush(stdout);
    }
    // current slab start time and time step
    Time     = TimeHistData->Time[iSlab];
    TimeStep = TimeHistData->TimeStep[iSlab];
    
    // set Time
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
    if (ierr != xf_OK) return ierr;
    
    // Set initial CFL
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart));
    if (ierr != xf_OK) return ierr;
    
    // Write out header for this time slab
    sprintf(PreHeader, "%s %.8E %s %8E", "%% DG Forward: Time slab start =", Time,
	    " Step =", TimeStep);
    ierr = xf_Error(xf_WriteLogHeader(All, PreHeader));
    if (ierr != xf_OK) return ierr;

    // check for dynamic refinement, apply to Uj
    ierr = xf_Error(xf_CheckSpatialDyno(All, SavePrefix, iSlab, OrderTime+1, Uj, &SpatialDyno));
    if (ierr != xf_OK) return ierr;

    // apply dynamic refinement to GCL and GCLj as well
    if ((UseGCL) && (SpatialDyno)){
      for (iOrder=0; iOrder<=OrderTime; iOrder++){
        ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, GCLj[iOrder], Uj[iOrder]));
        if (ierr != xf_OK) return ierr;
      }
      
      //Find GCL vector in All and project it
      ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
      if(ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, GCL, Uj[0]));
      if(ierr != xf_OK) return ierr;
    }

    // calculate spatial degrees of freedom on this slab
    if ((UErrEstOn) && (ErrEstData.sdofTime != NULL)) {
      ierr = xf_Error(xf_GetVectorDOF(Uj[0], ErrEstData.sdofTime+iSlab));
      if (ierr != xf_OK) return ierr;
    }

    // GCL step comes first
    if(UseGCL){

      ierr = xf_Error(xf_DGTimeGCLSolve(All, TimeScheme, Time, TimeStep, SolverData,
                                        GCLprev, GCLj));
      if (ierr != xf_OK) return ierr;
    }

    // Solve for states in time slab, Uj
    ierr = xf_Error(xf_DGinTimeStep(All, TimeScheme, Time, TimeStep, xfe_True,
				    SolverData, Uprev, Uj, &Redo));
    if (ierr != xf_OK) return ierr;

    // if Redo==True, an error occured; we will lower the time step and try again
    if (Redo){
      ierr = xf_Error(xf_SplitTimeHistData(TimeHistData, NULL, iSlab));
      if (ierr != xf_OK) return ierr;
      nTimeStep = TimeHistData->nTime;
      Redid++; // increment counter
      iSlab--; // so that time slab remains the same on next loop
      continue;
    }

    /* Store time history data (write out outputs for right node) */
    ierr = xf_Error(xf_StoreTimeHistData(All, TimeHistData, iSlab, Uj[OrderTime]));
    if (ierr != xf_OK) return ierr;

    /* Update error estimate */
    if (UErrEstOn){
      ierr = xf_Error(xf_UpdateUErrEst(All, iSlab, Time, TimeStep, Uj, Uprev, 
				       GCLj, GCLprev, SpatialDyno, SolverData, &ErrEstData));
      if (ierr != xf_OK) return ierr;
    }
    
    /* Destroy UExtra if it exists and spatial refinement changed */
    if ((UExtra != NULL) && (SpatialDyno)){
      ierr = xf_Error(xf_DestroyVector(UExtra, xfe_True));
      if (ierr != xf_OK) return ierr;
      UExtra = NULL;
    }
    // same with GCL
    if (UseGCL){
      if ((GCLExtra != NULL) && (SpatialDyno)){
        ierr = xf_Error(xf_DestroyVector(GCLExtra, xfe_True));
        if (ierr != xf_OK) return ierr;
        GCLExtra = NULL;
      } 
    }
    
    /* Set UExtra = vector similar to Uj */
    if (UExtra == NULL){
      ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], "UExtra", xfe_True, xfe_False,  
                                           NULL, &UExtra, NULL));
      if (ierr != xf_OK) return ierr;
    }
    // same with GCL
    if (UseGCL){
      if (GCLExtra == NULL){
        ierr = xf_Error(xf_FindSimilarVector(All, GCLj[0], "GCLExtra", xfe_True, xfe_False,  
                                             NULL, &GCLExtra, NULL));
        if (ierr != xf_OK) return ierr;
      }
    }
   
    /* Calculate unsteady outputs (UExtra is temp storage here) */
    ierr = xf_Error(xf_IncrementUnsteadyOutputs(All, NULL, TimeScheme, OrderTime+1, 
                                                Uj, UExtra, Time, TimeStep, (iSlab==0), 
                                                (iSlab==nTimeStep-1), &J, NULL));
    if (ierr != xf_OK) return ierr;

    // set Uprev = UExtra = Uj[OrderTime] = the right node from the previous time slab
    ierr = xf_Error(xf_SetVector(Uj[OrderTime], xfe_Set, UExtra));
    if (ierr != xf_OK) return ierr;
    Uprev = UExtra;
    if(UseGCL){ // same with GCL
      ierr = xf_Error(xf_SetVector(GCLj[OrderTime], xfe_Set, GCLExtra));
      if (ierr != xf_OK) return ierr;
      GCLprev = GCLExtra;
    }

    if(UseGCL){
      // Interpolate GCL to right time node prior to writing out (for plotting purposes)
      ierr = xf_Error(xf_DGTimeInterpolateGCL(All, TimeScheme, 1., NULL));
      if (ierr != xf_OK) return ierr;
    }

    /* Write out Uj[OrderTime] to hard disk if at requested interval */
    if ((SavePrefix != NULL) && (((iSlab+1) % WriteInterval) == 0)){
      sprintf(OutputFile, "%s_U%d.data\0", SavePrefix, iSlab+1); // 1 = first slab
      ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, OutputFile));
      if (ierr != xf_OK) return ierr;
    }

    // break out if user requests a halt
    if (xf_CheckUserHalt(NULL)) break;

  } // iSlab

  // print message if redid any time slabs
  if (Redid > 0) xf_printf("Redid (i.e. split) %d time slabs for robustness.\n", Redid);

  /* Finish error estimate */
  if (UErrEstOn){
    ierr = xf_Error(xf_FinishUErrEst(All, &ErrEstData));
    if (ierr != xf_OK) return ierr;
  }

  // No need to destroy the data in DataSet because it stores pointers to Data in All->DataSet
  if (DataSet != NULL){
    for (D=DataSet->Head; D != NULL; D=D->Next) D->Data = NULL;
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;
  }

  /* Destroy UExtra if it exists */
  if (UExtra != NULL){
    ierr = xf_Error(xf_DestroyVector(UExtra, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  /* Destroy GCLExtra if it exists */
  if (GCLExtra != NULL){
    ierr = xf_Error(xf_DestroyVector(GCLExtra, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  (*pSolverData) = SolverData;
  // release memory
  xf_Release( (void *) Uj);
  xf_Release( (void *) GCLj);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DGinTimeStepAdjoint
static int
xf_DGinTimeStepAdjoint(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
		       real Time, real TimeStep, int iSlab, int nSlab,
		       xf_SolverData *SolverData, xf_Vector **Uj, 
		       enum xfe_Bool InitFlag, enum xfe_Bool FineFlag, 
		       int nPsi, xf_Vector **Psi, xf_Vector **Psinext, 
		       xf_Vector ***Psii, enum xfe_Bool *Redo)
{
/*
PURPOSE:

  Takes a step of a DG in time unsteady adjoint solver.  On input, the
  latest adjoint (left node of future adjacent time slab) is stored in
  Psinext.  On output, the entire adjoint in the current time slab is
  computed and stored in Psii.
  
INPUTS:

  All         : all structure
  TimeScheme  : what time scheme to use
  Time        : current time
  TimeStep    : delta t
  iSlab       : current time slab
  nSlab       : number of slabs
  SolverData  : solver data structure (contains CFL)
  Uj          : state vectors on current time slab
  InitFlag    : if True, Psii will be initialized before solving
  FineFlag    : if True, fine space solver variables will be used
  nPsi        : number of adjoints
  Psi         : generic adjoint array storing OutputName
  Psinext     : adjoint on left node of next time slab
  
OUTPUTS: 

  Psii        : pointer to ALL the adjoints on the current time slab (0=left)
  Redo        : on error, this flag is set to True to redo the time step

RETURN: Error code

*/

  int ierr, nIter, iAdjoint, iq, nq;
  int iOrder, jOrder, OrderTime;
  enum xfe_Bool Found, LeanFlag, TransposeFlag;
  enum xfe_Bool Halted, Rewind, Converged;
  enum xfe_Verbosity Verbosity;
  enum xfe_PreconditionerType Preconditioner;
  char *OutputName;
  char StateName[xf_MAXSTRLEN];
  xf_JacobianMatrix *R_U;
  xf_KeyValue KeyValueOrig;
  real AdjointResidualTolerance, rtemp, omega;
  real tq[xf_MAXDGTIMENODE], wq[xf_MAXDGTIMENODE]; 
  real Atime[xf_MAXDGTIMENODE*xf_MAXDGTIMENODE], Vprev[xf_MAXDGTIMENODE];
  real phi[xf_MAXDGTIMENODE];
  xf_Vector **AdjRj, **Adji, **dAdji, *Adjnext;
  xf_Vector *Rtemp, *Utemp;

  (*Redo)   = xfe_False;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;

  TransposeFlag = xfe_True;

  // pull off variables for this time scheme
  ierr = xf_Error(xf_DGTimeSchemeVars(TimeScheme, Atime, Vprev, &nq, tq, wq, NULL));
  if (ierr != xf_OK) return ierr;

  // check input args
  if (nPsi <= 0) return xf_Error(xf_INPUT_ERROR);

  // Set desired solver flags for the fine space solve
  if (FineFlag){
    ierr = xf_Error(xf_PreFineSpaceSolve(All, &KeyValueOrig));
    if (ierr != xf_OK) return ierr;
  }

  // locate preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
				     xfe_PreconditionerName, (int ) xfe_PreconditionerLast, 
				     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;

  // Is the preconditioner memory-lean?
  ierr = xf_Error(xf_PreconditionerLeanCheck(Preconditioner, &LeanFlag));
  if (ierr != xf_OK) return ierr;

  // cannot use lean preconditioner for the approximate factorization
  if (LeanFlag){
    xf_printf("Lean preconditioners not supported for Adjoint DG in Time.\n");
    return xf_Error(xf_NOT_SUPPORTED);
  }

  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // pull off adjoint residual tolerance (not relative in this case)
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "AdjointResidualTolerance", 
			    &AdjointResidualTolerance);
  if (ierr != xf_OK) return ierr;

  // Under-relaxation factor
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "DGTimeUnderRelax", &omega);
  if (ierr != xf_OK) return ierr;

  // number of Adjoint iterations
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "nIterAdjoint", &nIter);
  if (ierr != xf_OK) return ierr;


  // locate array of required update vectors, dAdji
  ierr = xf_Error(xf_Alloc( (void **) &dAdji, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(StateName, "Update_DG_%d", iOrder);
    ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], StateName, xfe_True, xfe_True,  
					 NULL, dAdji + iOrder, NULL));
    if (ierr != xf_OK) return ierr;
  }

  // locate array of required residual vectors, AdjRj
  ierr = xf_Error(xf_Alloc( (void **) &AdjRj, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(StateName, "Residual_DG_%d", iOrder);
    ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], StateName, xfe_True, xfe_True,
					 NULL, AdjRj + iOrder,NULL));
    if (ierr != xf_OK) return ierr;
  }

  // locate temporary state vector
  ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], "Utemp_DG", xfe_True, xfe_True, 
				       NULL, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;

  // locate temporary residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], "Rtemp_DG", xfe_True, xfe_True, 
				       NULL, &Rtemp, NULL));
  if (ierr != xf_OK) return ierr;
  
 
  // locate Jacobian and aux vectors; check size; create if necessary
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, Uj[0], NULL,
					!LeanFlag, NULL, &R_U, &Found));
  if (ierr != xf_OK) return ierr;

  /***  Loop over adjoints  ***/
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){

    // current adjoint vectors
    Adji    = Psii[iAdjoint]; 
    Adjnext = Psinext[iAdjoint];
    
    // associated output name (make sure not NULL)
    if ((OutputName = Psi[iAdjoint]->OutputName) == NULL) 
      return xf_Error(xf_INPUT_ERROR);
    
/*     if ((iSlab != nSlab-1) && (OutputName = Adjnext->OutputName) == NULL)  */
/*       return xf_Error(xf_INPUT_ERROR); */


    /*** Initialize Adji  ***/
    if (InitFlag){
      if (iSlab == nSlab-1){ 
	// on last time slab initizlize Adji to zero
	for (iOrder = 0; iOrder <= OrderTime; iOrder++){
	  ierr = xf_Error(xf_SetZeroVector(Adji[iOrder]));
	  if (ierr != xf_OK) return ierr;
	} // iOrder
      }
      else{
	// on all other time slabs initialize using Adjnext
	for (iOrder = 0; iOrder <= OrderTime; iOrder++){
	  // using project to allow for varying orders
	  ierr = xf_Error(xf_ProjectVector(All, Adjnext, xfe_False, Adji[iOrder]));
	  if (ierr != xf_OK) return ierr;
	} // iOrder
      }
    }

    
    /*** Begin solver loop ***/

    Halted    = xfe_False;
    Rewind    = xfe_False;
    Converged = xfe_False;

    for (SolverData->iIter=0; SolverData->iIter < nIter; SolverData->iIter++){

      if (Halted = xf_CheckUserHalt(NULL)) break;
	
      if (Rewind){
	xf_printf("DGTime Adjoint requested a Rewind.  Not currently supported.\n");
	return xf_Error(xf_NOT_SUPPORTED);
      }

      /*** Calculate RHS (on left) for adjoint equations  ***/
      
      // initialize AdjRj = 0
      for (iOrder = 0; iOrder <= OrderTime; iOrder++){
	ierr = xf_Error(xf_SetZeroVector(AdjRj[iOrder]));
	if (ierr != xf_OK) return ierr;
      } // iOrder
      
      // calculate output contributions first (Time used correctly in this function)
      ierr = xf_Error(xf_IncrementUnsteadyOutputs(All, OutputName, TimeScheme, 
						  OrderTime+1, Uj, Utemp, 
						  Time, TimeStep, (nSlab==0), 
						  (iSlab==nSlab-1), NULL, AdjRj));
      if (ierr != xf_OK) return ierr; 
      
      /* add contributions from stiffness matrix terms and influence of Adjnext*/
      for (jOrder = 0; jOrder <= OrderTime; jOrder++){
	// temporal stiffness matrix
	for (iOrder = 0; iOrder <= OrderTime; iOrder++){
	  ierr = xf_Error(xf_VectorMultSet(Adji[iOrder],  Atime[iOrder*(OrderTime+1)+jOrder], 
					   ((iOrder==0) ? xfe_Set : xfe_Add), Utemp));
	  if (ierr != xf_OK) return ierr;
	}
	// add M*Utemp to adjoint residual AdjRj[jOrder]
	ierr = xf_Error(xf_AddMassMatrix(All, 1.0, NULL, Utemp, AdjRj[jOrder], NULL, NULL));
	if (ierr != xf_OK) return ierr;
	// influence of Adjnext
	if ((iSlab != nSlab-1) && (Vprev[OrderTime-jOrder] != 0.)){
	  // add M'*Adjnext to AdjRj[jOrder]
	  ierr = xf_Error(xf_AddMassMatrix(All, Vprev[OrderTime-jOrder], NULL, Adjnext, 
					   AdjRj[jOrder], NULL, NULL));
	  if (ierr != xf_OK) return ierr;
	}
      }


      /* add contributions from quadrature points */

      for (iq=0; iq<nq; iq++){ // loop over (the two) quadrature points

	// Set Time to that at the quadrature point, Time + tq[iq]*TimeStep
	ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", 
					   Time+tq[iq]*TimeStep));
	if (ierr != xf_OK) return ierr;
	
	// Calculate U at the quad point Utemp = U(tq), and phii(tq); GCL is automatically interpolated as well
	ierr = xf_Error(xf_DGTimeInterpolateState(All,TimeScheme, Uj, tq+iq, &Utemp, phi));
	if (ierr != xf_OK) return ierr;
	
	// Calculate the Jacobian at the iq quadrature state
	ierr = xf_CalculateResidual(All, Utemp, Rtemp, R_U, SolverData);
	if (xf_CheckSolverErrorUnsteady(ierr, xfe_True, SolverData, 1, NULL, 
					OrderTime+1, NULL, &TimeStep, &Rewind)) 
	return xf_Error(xf_SOLVER_ERROR);
	if (Rewind) break;

	// Calculate Adjoint at the quad point, Utemp = Adj(tq)
	ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Adji, 1, tq+iq, &Utemp, NULL));
	if (ierr != xf_OK) return ierr;

	// Rtemp = R_U*Utemp
	ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Utemp, xfe_Set, 
					 TransposeFlag, SolverData, Rtemp));
	if (ierr != xf_OK) return ierr;

	// AdjRj[j] += wq(iq)*phij(tq)*TimeStep*Rtemp
	for (jOrder = 0; jOrder <= OrderTime; jOrder++){
	  ierr = xf_Error(xf_VectorMultSet(Rtemp, wq[iq]*phi[jOrder]*TimeStep,
					   xfe_Add, AdjRj[jOrder]));
	  if (ierr != xf_OK) return ierr;
	}

      }
      if (Rewind) break;


      // compute residual norm
      for (iOrder=0, SolverData->ResNorm=0.; iOrder<=OrderTime; iOrder++){
	ierr = xf_Error(xf_VectorNorm(AdjRj[iOrder], 1, &rtemp));
	if (ierr != xf_OK) return ierr;
	SolverData->ResNorm += rtemp;
      }
	  
      // print out residual norm to stdout
      if (Verbosity != xfe_VerbosityLow)
	xf_printf("  DG in Time Adjoint iteration %d: |R| = %.10E\n", 
		  SolverData->iIter, SolverData->ResNorm);

      // convergence check
      if (SolverData->ResNorm < AdjointResidualTolerance){
	Converged = xfe_True;
	break;
      }


      // Solve DG linear system using approximate factorization
      ierr = xf_Error(xf_DGApproxSolver(All, TimeScheme, Time, TimeStep, 
					SolverData, Uj, AdjRj, dAdji, TransposeFlag));
      if (ierr == xf_REWIND){
	(*Redo) = xfe_True;
	break;
      }
      if (ierr != xf_OK) return ierr;
      
      // under-relax the updates (1.0 is default)
      for (iOrder = 0; iOrder <= OrderTime; iOrder++){
	ierr = xf_Error(xf_VectorMult(dAdji[iOrder], omega));
 	if (ierr != xf_OK) return ierr;
      } // iOrder

      
      // update Adjoint, no need to check if physical
      for (iOrder = 0; iOrder <= OrderTime; iOrder++){
	ierr = xf_Error(xf_SetVector(dAdji[iOrder], xfe_Add, Adji[iOrder]));
	if (ierr != xf_OK) return ierr;
      } // iOrder
    
    } // iter
    
    // Reset Time to start of time slab
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
    if (ierr != xf_OK) return ierr;
    
    if (Rewind) (*Redo) = xfe_True; // need to do time step over
    
    // notify if converged
    if ((Converged) && (Verbosity != xfe_VerbosityLow))
      xf_printf("DGTime Adjoint step (%s) converged to tolerance.\n", OutputName);

    // notify if not converged (not a problem if on fine space)
    if ((!Converged) && (FineFlag == xfe_False)){
      if (Verbosity != xfe_VerbosityLow) 
	xf_printf("Warning, DGTime Adjoint (%s) not converged at Time = %.10E, TimeStep = %.10E\n",
		  OutputName, Time, TimeStep);
      (*Redo) = xfe_True;
    }

    // notify if redoing
    if (*Redo){
      xf_printf("Asking to redo this adjoint time slab.\n");
      break; // out of iAdjoint loop
    }

  } // iAdjoint

  // Reset solver flags to original
  if (FineFlag){
    ierr = xf_Error(xf_PostFineSpaceSolve(All, &KeyValueOrig));
    if (ierr != xf_OK) return ierr;
  }

  // memory cleanup
  xf_Release( (void *) AdjRj);
  xf_Release( (void *) dAdji);

  return xf_OK;

}



/******************************************************************/
//  FUNCTION Definition: xf_ApplyTimeSchemeAdjoint_DG
static int
xf_ApplyTimeSchemeAdjoint_DG(xf_All *All, const char *SavePrefix,
			     xf_Vector *U0, int nPsi, xf_Vector **Psi, 
			     xf_TimeHistData *TimeHistData)
{
  int ierr, nSlab, iSlab, iOrder, OrderTime, OrderTimeh, OrderTimeIh;
  int WriteInterval, iAdjoint, nSumWeight;
  int nSubSlab = 1, iSubSlab, nSubSlab0, Redid = 0;
  int SpaceOrderIncrement = 0;
  int TimeOrderIncrement = 0;
  int UEIWriteInterval = 0;
  int *sdofTime = NULL;
  enum xfe_Bool LinearFlag;
  enum xfe_Bool UErrEstOn = xfe_False;
  enum xfe_Bool ErrEstUseReconstruct = xfe_False;
  enum xfe_Bool UErrEstIterative = xfe_False;
  enum xfe_Bool UErrEstIterativeGS = xfe_False;
  enum xfe_Bool UErrEstUseReconstruct = xfe_False;
  enum xfe_Bool UErrEstConservative = xfe_False;
  enum xfe_Bool SpatialDyno = xfe_False;
  enum xfe_Bool btemp, Redo = xfe_False;
  enum xfe_Bool UseGCL = xfe_False;
  enum xfe_Bool MotionOn = xfe_False;
  enum xfe_TimeSchemeType TimeScheme;
  enum xfe_TimeSchemeType TimeSchemeh, TimeSchemeIh;
  enum xfe_SpaceTimeAnisoType UErrEstAnisoMeasure;
  enum xfe_Verbosity Verbosity;
  char Title[xf_MAXSTRLEN];
  char PreHeader[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  real Time, TimeStep, CFLStart;
  real SlabTimeStart, SlabTimeStep;
  real *TimeWeights, **SumOutputWeights;
  real *OutputError = NULL;
  real *OutputErrorSpace = NULL, *OutputErrorTime = NULL;
  real **ErrIndTime = NULL;
  real xt[xf_MAXDGTIMENODE];
  real LinResTol;
  real sumSpace, sumTime;
  void *voidpointer = NULL; // used for swapping D->Data
  xf_Vector ***Psii; // vector of all adjoint vectors (at each time index)
  xf_Vector ***PsiGCLi = NULL, **PsiGCL = NULL;
  xf_Vector **Uj = NULL, **SUj = NULL;
  xf_Vector **AdjiFine = NULL;
  xf_Vector **AdjGCLiFine = NULL;
  xf_Vector **Psinext;   // next time slab first node adjoints
  xf_Vector **PsiGCLnext = NULL;
  xf_Vector **PsiExtra;  // extra vectors for adjoint next storage
  xf_Vector **PsiGCLExtra = NULL;  
  xf_Vector **PsinextFine;   // next time slab first node adjoints (fine space)
  xf_Vector **PsinextIterative;
  xf_Vector *Uprev, *SUprev;
  xf_Vector **ErrIndElemSTot = NULL;
  xf_Vector **ErrIndElem = NULL;
  xf_Vector **ErrIndElemS = NULL;
  xf_Vector **ErrIndElemT = NULL;
  xf_Vector *SpaceTimePref = NULL;
  xf_Vector *ETemp = NULL, *rvec = NULL;
  xf_Vector *GCL = NULL; // current GCL vector
  xf_Vector **GCLj = NULL; // GCL vectors in All struct
  xf_Vector *GCLprev0 = NULL; // GCL vector on previous time slab
  xf_Vector *GCLprev = NULL; // another pointer to GCL on previous time slab
  xf_Vector **GCLjTemp = NULL;
  xf_Vector *J_GCL = NULL;
  xf_Output *Output;
  xf_Data *D, **Dj;
  xf_DataSet *DataSetPsi;
  xf_DataSet *DataSet = NULL;
  xf_DataSet *DataSetPrev = NULL;
  xf_SolverData *SolverData;


  // set convenient variables
  nSlab            = TimeHistData->nTime; // # of time slabs
  TimeWeights      = TimeHistData->TimeWeights;
  nSumWeight       = TimeHistData->nSumWeight;
  SumOutputWeights = TimeHistData->SumOutputWeights;

  // do not currently support SumOutput with DG Adjoint
  if ((nSumWeight > 0) && (SumOutputWeights != NULL))
    return xf_Error(xf_NOT_SUPPORTED);

  // Is the system linear (as prescribed by the user)?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "LinearFlag", &LinearFlag));
  if (ierr != xf_OK) return ierr;
  
  // Are we using a Geometric Conservation Law?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;

  // Are we doing mesh motion?
  MotionOn = ((All->Mesh->Motion != NULL) && (All->Mesh->Motion->Active));
  
  // Unsteady write interval
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UnsteadyWriteInterval", 
				    &WriteInterval));
  if (ierr != xf_OK) return ierr;

  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;


  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
				     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
				     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;


  // Are we doing unsteady error estimation?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UErrEstOn", &UErrEstOn));
  if (ierr != xf_OK) return ierr;

  // if so, pull off several other parameters
  if (UErrEstOn){
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UErrEstConservative", 
				       &UErrEstConservative));
    if (ierr != xf_OK) return ierr;
    // Are we iteratively correcting a coarse adjoint?
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UErrEstIterative", 
				       &UErrEstIterative));
    if (ierr != xf_OK) return ierr;
    // Are we propagating Psinext in iterative adjoint solution (Gauss-Seidel in time)?
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UErrEstIterativeGS", 
				       &UErrEstIterativeGS));
    if (ierr != xf_OK) return ierr;

    // default fine-space time schemes
    TimeSchemeh  = TimeScheme;
    TimeSchemeIh = TimeScheme;

    // number of time sub-slabs
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UErrEstnSubSlab", 
				      &nSubSlab));
    if (ierr != xf_OK) return ierr;
    if (nSubSlab < 1) return xf_Error(xf_INPUT_ERROR);

    // de we (also) want HO in time?
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UErrEstOrderIncrement", 
				      &TimeOrderIncrement));
    if (ierr != xf_OK) return ierr;

    if (TimeOrderIncrement > 0){
      // only support DG1 -> DG2
      if (TimeOrderIncrement > 1) return xf_Error(xf_NOT_SUPPORTED);
      if (TimeScheme != xfe_TimeSchemeDG1) return xf_Error(xf_NOT_SUPPORTED);
      // Are we reconstructing to get temporal high-order?
      ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UErrEstUseReconstruct", 
					 &UErrEstUseReconstruct));
      if ((UErrEstIterative) || (UErrEstUseReconstruct)){
	// will be solving DG1 and iterating/reconstructing DG2
	TimeSchemeh  = xfe_TimeSchemeDG1; // same as TimeScheme
	TimeSchemeIh = xfe_TimeSchemeDG2;
      }
      else{ // will be solving on DG2
	TimeSchemeh  = xfe_TimeSchemeDG2;
	TimeSchemeIh = xfe_TimeSchemeDG2;
      }
    }

    // Spatial order increment
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement", 
				      &SpaceOrderIncrement));
    if (ierr != xf_OK) return ierr;
    // Are we reconstructing to get spatial high-order?
    ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "ErrEstUseReconstruct", 
				       &ErrEstUseReconstruct));
    if (ierr != xf_OK) return ierr;
    // write interval for unsteady error indicator
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UErrEstIndicatorWriteInterval", 
				      &UEIWriteInterval));
    if (ierr != xf_OK) return ierr;
    
    // Allocate vector of pointers for error estimation vectors
    ierr = xf_Error(xf_Alloc( (void **) &ErrIndElemSTot, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &ErrIndElem, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &ErrIndElemS, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &ErrIndElemT, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      // create spatial error indicator vector
      sprintf(Title, "ErrIndSpaceTot_%s", Psi[iAdjoint]->OutputName);
      ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				    NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  xfe_True, &D, 
				    ErrIndElemSTot+iAdjoint, NULL));
      if (ierr != xf_OK) return ierr;
      D->ReadWrite = xfe_True; // make indicator writeable
      // zero out spatial error indicator (running total, with option of absolute values)
      ierr = xf_Error(xf_SetZeroVector(ErrIndElemSTot[iAdjoint]));
      if (ierr != xf_OK) return ierr;
      // will need another error indicator vector, for local storage (not added to dataset)
      sprintf(Title, "ErrIndElem_%s", Psi[iAdjoint]->OutputName);
      ierr = xf_Error(xf_FindSimilarVector(All, ErrIndElemSTot[iAdjoint], Title, xfe_False, 
					   xfe_False,  NULL, ErrIndElem+iAdjoint, NULL));
      if (ierr != xf_OK) return ierr;
      // will need another error indicator vector, for local storage (not added to dataset)
      sprintf(Title, "ErrIndElemS_%s", Psi[iAdjoint]->OutputName);
      ierr = xf_Error(xf_FindSimilarVector(All, ErrIndElemSTot[iAdjoint], Title, xfe_False, 
					   xfe_False,  NULL, ErrIndElemS+iAdjoint, NULL));
      if (ierr != xf_OK) return ierr;
      // will need another error indicator vector, for local storage (not added to dataset)
      sprintf(Title, "ErrIndElemT_%s", Psi[iAdjoint]->OutputName);
      ierr = xf_Error(xf_FindSimilarVector(All, ErrIndElemSTot[iAdjoint], Title, xfe_False, 
					   xfe_False,  NULL, ErrIndElemT+iAdjoint, NULL));
      if (ierr != xf_OK) return ierr;

    }
    // Allocate reals for storing temporally-localized error indicator
    ierr = xf_Error(xf_Alloc2( (void ***) &ErrIndTime, nSlab, nPsi, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    // Allocate ints for storing spatial dof for each time slab
    ierr = xf_Error(xf_Alloc( (void **) &sdofTime, nSlab, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    // Allocate reals for storing non-localized, unadulterated, output error estimates
    ierr = xf_Error(xf_Alloc( (void **) &OutputError, 3*nPsi, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    OutputErrorSpace = OutputError + nPsi;
    OutputErrorTime  = OutputError + nPsi*2;
    for (iAdjoint=0; iAdjoint<3*nPsi; iAdjoint++) OutputError[iAdjoint] = 0.;
    // pull off UErrEstAnisoMeasure
    ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "UErrEstAnisoMeasure", 
				       xfe_SpaceTimeAnisoName, (int) xfe_SpaceTimeAnisoLast, 
				       (int *) &UErrEstAnisoMeasure));
    if (ierr != xf_OK) return ierr;
    // create spatial vs. temporal preference vector (not added to dataset)
    sprintf(Title, "SpaceTimePref_State");
    ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				  NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, xfe_False, NULL, 
				  &SpaceTimePref, NULL));
    if (ierr != xf_OK) return ierr;
    // will need another temp error indicator vector
    sprintf(Title, "ETemp");
    ierr = xf_Error(xf_FindSimilarVector(All, ErrIndElemSTot[0], Title, xfe_False, 
					 xfe_False,  NULL, &ETemp, NULL));
    if (ierr != xf_OK) return ierr;
    // need PsinextFine pointers
    ierr = xf_Error(xf_Alloc( (void **) &PsinextFine, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) PsinextFine[iAdjoint] = NULL;
  }
  else{
    nSubSlab = 1;
    ErrEstUseReconstruct = xfe_False;
    SpaceOrderIncrement = 0;
    TimeSchemeh = TimeScheme;
    TimeSchemeIh = TimeScheme;
  }

  // determine fine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeSchemeh, &OrderTimeh));
  if (ierr != xf_OK) return ierr;

  // determine iterative/reconstructed fine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeSchemeIh, &OrderTimeIh));
  if (ierr != xf_OK) return ierr;
  
  // need PsiExtra pointers, initialize to NULL
  ierr = xf_Error(xf_Alloc( (void **) &PsiExtra, nPsi, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) PsiExtra[iAdjoint] = NULL;

  // same with GCL
  if (UseGCL){
    ierr = xf_Error(xf_Alloc( (void **) &PsiGCLExtra, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) PsiGCLExtra[iAdjoint] = NULL;
  }


  // Allocate vector of adjoint pointers
  ierr = xf_Error(xf_Alloc( (void **) &Psii, nPsi, sizeof(xf_Vector **)));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){ // same with GCL if using it
    ierr = xf_Error(xf_Alloc( (void **) &PsiGCLi, nPsi, sizeof(xf_Vector **)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &PsiGCL, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
  }


  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    // allocate space for adjoint vector pointers
    ierr = xf_Error(xf_Alloc( (void **) (Psii+iAdjoint), OrderTimeIh+1, 
			      sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    // locate required adjoint vectors; create if necessary
    for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
      sprintf(Title, "%s_Adjoint_%d", Psi[iAdjoint]->OutputName, iOrder);
      ierr = xf_Error(xf_FindSimilarVector(All, Psi[iAdjoint], Title, xfe_True, 
					   xfe_True,  NULL, Psii[iAdjoint] + iOrder, NULL));
      if (ierr != xf_OK) return ierr;
    } // iOrder

    // same with GCL
    if (UseGCL){
      ierr = xf_Error(xf_Alloc( (void **) (PsiGCLi+iAdjoint), OrderTimeIh+1, 
				sizeof(xf_Vector *)));
      if (ierr != xf_OK) return ierr;
      // find a GCL adjoint vector based on regular adjoint, Psi
      ierr = xf_Error(xf_InitMeshMotionGCLAdjoint(All, Psi[iAdjoint], -1, xfe_True, PsiGCLi[iAdjoint]));
      if (ierr != xf_OK) return ierr;
      PsiGCL[iAdjoint] = PsiGCLi[iAdjoint][0];
      // locate and initialize GCL adjoint vectors 
      for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
        ierr = xf_Error(xf_InitMeshMotionGCLAdjoint(All, Psi[iAdjoint], iOrder, xfe_True, PsiGCLi[iAdjoint]+iOrder));
        if (ierr != xf_OK) return ierr;
      }
    }

  } // iAdjoint


  // need Psinext pointers, initialize to Psi (input)
  ierr = xf_Error(xf_Alloc( (void **) &Psinext, nPsi, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) Psinext[iAdjoint] = Psi[iAdjoint];

  // same with GCL
  if (UseGCL){
    ierr = xf_Error(xf_Alloc( (void **) &PsiGCLnext, nPsi, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) PsiGCLnext[iAdjoint] = PsiGCL[iAdjoint];
  }


  // Project Psinext to HO if doing errest and not reconstructing spatially
  if ((UErrEstOn) && (!ErrEstUseReconstruct) && (!UErrEstIterative) && (SpaceOrderIncrement > 0)){
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Psinext[iAdjoint], 
					      NULL, xfe_BasisLast, SpaceOrderIncrement));
      if (ierr != xf_OK) return ierr;
    }
    if (UseGCL){
      // same with GCL
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
        ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, PsiGCLnext[iAdjoint], 
                                                               NULL, xfe_BasisLast, SpaceOrderIncrement));
        if (ierr != xf_OK) return ierr;
      }
    }
  }


  /* Create a dataset for writing the adjoint at each WriteInterval */
  ierr = xf_Error(xf_CreateDataSet(&DataSetPsi));
  if (ierr != xf_OK) return ierr;
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
      sprintf(Title, "%s_Adjoint_%d", Psi[iAdjoint]->OutputName, iOrder);
      ierr = xf_Error(xf_DataSetAdd(DataSetPsi, Title, xfe_Vector, xfe_True,
				    (void *) Psii[iAdjoint][iOrder], NULL));
      if (ierr != xf_OK) return ierr;
    } // iOrder
    if (UseGCL){
      for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
        sprintf(Title, "%s_GCLAdjoint_%d", Psi[iAdjoint]->OutputName, iOrder);
        ierr = xf_Error(xf_DataSetAdd(DataSetPsi, Title, xfe_Vector, xfe_True,
                                      (void *) PsiGCLi[iAdjoint][iOrder], NULL));
        if (ierr != xf_OK) return ierr;
      } // iOrder
    }
  } // iAdjoint

  
  // Allocate vector of state pointers
  ierr = xf_Error(xf_Alloc( (void **) &Uj, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  
  if(UseGCL){
    // Allocate vector of pointers to GCLj vectors in All struct
    ierr = xf_Error(xf_Alloc( (void **) &GCLj, OrderTimeIh+1, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    
    // Allocate vector of pointers to GCLj temporary vectors used in error estimation
    ierr = xf_Error(xf_Alloc( (void **) &GCLjTemp, OrderTimeIh+1, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;

    // ensure that GCL vector exists
    ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U0, -1, xfe_True, &GCL));
    if (ierr != xf_OK) return ierr;

    // ensure that GCLj (up to OrderTime) vectors exists
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      ierr = xf_Error(xf_InitMeshMotionGCLVector(All, U0, iOrder, xfe_True, GCLj + iOrder));
      if (ierr != xf_OK) return ierr;
    }
  }

  // need sub-slab vector points, even if not doing error estimation
  ierr = xf_Error(xf_Alloc( (void **) &SUj, OrderTimeIh+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;

  // need fine-space adjoint pointers, even if not doing error estimation
  ierr = xf_Error(xf_Alloc( (void **) &AdjiFine, OrderTimeIh+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){ // same with GCL
    ierr = xf_Error(xf_Alloc( (void **) &AdjGCLiFine, OrderTimeIh+1, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
  }
  
  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;

  // pull off adjoint residual tolerance (relative)
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "AdjointResidualTolerance", 
			    &LinResTol);
  if (ierr != xf_OK) return ierr;
  SolverData->LinResTol = LinResTol;

  // pull off CFL, CFLIncreaseFactor, CFLDecreaseFactor, etc.
  ierr = xf_Error(xf_FindCFLData(All->Param->KeyValue, SolverData));
  if (ierr != xf_OK) return ierr;

  // Initial CFL for artificial time stepping
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFLStart));
  if (ierr != xf_OK) return ierr;

  // original # of sub slabs requested (can change during solve for robustness)
  nSubSlab0 = nSubSlab;
  Redo = xfe_False;

  /*----------------------------------------*/
  /* Begin loop over time slabs (backwards) */
  /*----------------------------------------*/
  
  for (iSlab=nSlab-1; iSlab>=0; iSlab--){
    
    // break out if user requests a halt
    if (xf_CheckUserHalt(NULL)) break;

    // initialize Redo
    Redo = xfe_False;

    // Set Time to start of slab
    SlabTimeStart = TimeHistData->Time[iSlab];
    SlabTimeStep  = TimeHistData->TimeStep[iSlab];
    
    // Read in U from storage if problem is not linear, if doing error estimation, or if using GCL
    if (((LinearFlag) && (!UErrEstOn)) && (!UseGCL) && (!MotionOn))
      for (iOrder = 0; iOrder <= OrderTime; iOrder++)
	Uj[iOrder] = U0; // use the default input state
    else{
      // need SavePrefix
      if (SavePrefix == NULL) return xf_Error(xf_INPUT_ERROR);

      // read .data from file (multiple vectors in one dataset)
      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;
      sprintf(Title, "%s_U%d.data\0", SavePrefix, iSlab+1);
      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, Title, DataSet));
      if (ierr != xf_OK) return ierr;
      
      // use data vectors in order of data
      for (iOrder=0, D=DataSet->Head; iOrder<=OrderTime; iOrder++, D=D->Next)
	Uj[iOrder] = (xf_Vector *) D->Data;
      if ((!UseGCL) && (D != NULL)) return xf_Error(xf_CODE_LOGIC_ERROR); 
    }
	   
    if(UseGCL){
      // Find GCLj vectors in All struct
      ierr = xf_Error(xf_DGTimeFindGCLVectors(All, TimeScheme, NULL, &Dj)); 
      if(ierr != xf_OK) return ierr;

      // Read in GCL vectors from DataSet and swap with those in All->DataSet
      for (iOrder=0; iOrder<=OrderTime; iOrder++, D=D->Next){
        swap(Dj[iOrder]->Data, D->Data, voidpointer);
        GCLj[iOrder] = (xf_Vector *) Dj[iOrder]->Data;
      }
	      
      // Destroy Dj
      xf_Release( (void *) Dj);

      // sanity check ... should be at end of DataSet
      if (D != NULL) return xf_Error(xf_CODE_LOGIC_ERROR);       
    }

    // read in Uprev for error estimation
    if (UErrEstOn){
      // read .data from file 
      ierr = xf_Error(xf_CreateDataSet(&DataSetPrev));
      if (ierr != xf_OK) return ierr;
      sprintf(Title, "%s_U%d.data\0", SavePrefix, iSlab);
      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, Title, DataSetPrev));
      if (ierr != xf_OK) return ierr;
      
      // use last, or only in case of (iSlab==0), data vector
      D = DataSetPrev->Head;
      if (iSlab != 0) for (iOrder=0; (iOrder<OrderTime)&&(D!=NULL); iOrder++, D=D->Next);
      if (D == NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
      Uprev = (xf_Vector *) D->Data;
      
      // also need GCLprev0 for error estimation
      if (UseGCL){
        D=D->Next;
        if (iSlab != 0) for (iOrder=0; (iOrder<OrderTime)&&(D!=NULL); iOrder++, D=D->Next);
        if (D == NULL) return xf_Error(xf_CODE_LOGIC_ERROR); // GCL needs to be in there
        if (D->Next != NULL) return xf_Error(xf_CODE_LOGIC_ERROR); // nothing should come after the GCL (for now)
        GCLprev0 = (xf_Vector *) D->Data; // leave in Data, will get deleted with DataSetPrev
      }
    }

    // Estimate anisotropy (space vs. time) of error (e.g. via jumps)
    if (UErrEstOn && (UErrEstAnisoMeasure == xfe_SpaceTimeAnisoJump)){
      ierr = xf_Error(xf_ErrEstSpaceTimeAniso(All, TimeScheme, Uj, Uprev, 
					      SlabTimeStart, SlabTimeStep,
					      SpaceTimePref, NULL));
      if (ierr != xf_OK) return ierr;
    }
    
    // calculate spatial degrees of freedom on this slab
    if (sdofTime != NULL){
      ierr = xf_Error(xf_GetVectorDOF(Uj[0], sdofTime+iSlab));
      if (ierr != xf_OK) return ierr;
    }

    
    // Project Uj to HO if not planning to reconstruct spatially
    if ((UErrEstOn) && (!ErrEstUseReconstruct) && (!UErrEstIterative) && (SpaceOrderIncrement > 0)){
      
      for (iOrder = 0; iOrder <= OrderTime; iOrder++){
	ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Uj[iOrder], NULL, 
							       xfe_BasisLast, SpaceOrderIncrement));
        if (ierr != xf_OK) return ierr;
      }
      
      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Uprev, NULL, 
							     xfe_BasisLast, SpaceOrderIncrement));
      if (ierr != xf_OK) return ierr;

      // same with GCL
      if(UseGCL){
        // spatially project GCLj vectors [still OrderTime]
        for (iOrder = 0; iOrder <= OrderTime; iOrder++){
          ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLj[iOrder], NULL, 
                                                                 xfe_BasisLast, SpaceOrderIncrement));
          if (ierr != xf_OK) return ierr;
        } // iOrder

        // spatially project GCL (note, this vector persists, so here we just make it look like GCLj[0])
        ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, GCL, GCLj[0]));
        if (ierr != xf_OK) return ierr;

        // spatially project GCLprev0 vector
        ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLprev0, NULL, 
                                                               xfe_BasisLast, SpaceOrderIncrement));
        if (ierr != xf_OK) return ierr;
      }
    } 
    
    // project Psii to same space (orders) as Uj if different
    SpatialDyno = xfe_False;
    if ((nPsi > 0) && (!xf_CompatibleVectors(Psii[0][0], Uj[0]))) {
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
        xf_printf("projecting Adjoint vectors.\n");
	for (iOrder=0; iOrder<=OrderTimeIh; iOrder++){
	  // project Psii[iAdjoint][iOrder] to make it look like Uj[0]
	  ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, 
							 Psii[iAdjoint][iOrder], Uj[0]));
	  if (ierr != xf_OK) return ierr;
          // same with GCL; note that only order info is used from Uj[0]
          if (UseGCL){
            ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, 
                                                           PsiGCLi[iAdjoint][iOrder], Uj[0]));
            if (ierr != xf_OK) return ierr;
          }
	} // iOrder
        // also project mesh motion J_GCL vectors for output linearization calculation
        if (UseGCL){
          for (iOrder=-1; iOrder<=OrderTimeIh; iOrder++){
            ierr = xf_Error(xf_FindMeshMotionGCLLinearization(All, Psi[iAdjoint]->OutputName, iOrder, &J_GCL));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, 
                                                           J_GCL, Uj[0]));
            if (ierr != xf_OK) return ierr;
          }
        }
      } // iAdjoint
      if (!xf_CompatibleVectors(Psii[0][0], Uj[0])) return xf_Error(xf_CODE_LOGIC_ERROR);
      SpatialDyno = xfe_True;
    }


    if (UErrEstOn){
      // SUj = additional vectors for sub-slab states (not added to All->DataSet)
      for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
	sprintf(Title, "SubState_%d", iOrder);
	// create required state vectors
	ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], Title, xfe_True, xfe_False,
					     NULL, SUj + iOrder, NULL));
	if (ierr != xf_OK) return ierr;
      }

      // make sure GCLj has extra vectors if OrderTimeIh > OrderTime
      if ((UseGCL) && (OrderTimeIh > OrderTime)){
        for (iOrder = OrderTime+1; iOrder <= OrderTimeIh; iOrder++){
          // note, these will be at appropriate (high) spatial order because using Uj[0] here
          ierr = xf_Error(xf_InitMeshMotionGCLVector(All, Uj[0], iOrder, xfe_True, GCLj + iOrder));
          if (ierr != xf_OK) return ierr;
        }
      }

      // zero out element-localized error indicators
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
	ierr = xf_Error(xf_SetZeroVector(ErrIndElem[iAdjoint]));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetZeroVector(ErrIndElemS[iAdjoint]));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_SetZeroVector(ErrIndElemT[iAdjoint]));
	if (ierr != xf_OK) return ierr;
      }

      // zero out slab-localized error indicator
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++)
	ErrIndTime[iSlab][iAdjoint] = 0.0;
	
      // Fill temporary GCLjTemp vectors [OrderTimeIh, on original (macro) slub]
      if (UseGCL){	    
        for(iOrder=0; iOrder<=OrderTimeIh; iOrder++){
          // Allocate space for GCLjTemp vectors and copy GCLj data into them
          sprintf(Title, "GCLjTemp%d_DG%d", iOrder, OrderTimeIh);
          ierr = xf_Error(xf_FindSimilarVector(All, GCLj[iOrder], Title, xfe_True, 
                                               xfe_True, NULL, GCLjTemp+iOrder, NULL));
          if(ierr != xf_OK) return ierr;
          
          ierr = xf_Error(xf_SetVector(GCLj[iOrder], xfe_Set, GCLjTemp[iOrder]));
          if(ierr != xf_OK) return ierr;
        }
      }
	
    } // UErrEstOn
    

    // Loop over sub-slabs (backwards)
    for (iSubSlab=nSubSlab-1; iSubSlab>=0; iSubSlab--){
    
      // sub-slab time step
      TimeStep = SlabTimeStep/nSubSlab;

      // time at start of sub-slab
      Time = SlabTimeStart + TimeStep*iSubSlab;
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
      if (ierr != xf_OK) return ierr;

      // Write out header for this time slab
      sprintf(PreHeader, "%s %.8E %s %8E", "%% DG Adjoint: Time slab start =", Time,
	      " Step =", TimeStep);
      ierr = xf_Error(xf_WriteLogHeader(All, PreHeader));
      if (ierr != xf_OK) return ierr;
      
      if (UErrEstOn){
	// node positions of sub-slab within slab
	for (iOrder=0; iOrder<=OrderTimeh; iOrder++)
	  xt[iOrder] = ((real) (iSubSlab*OrderTimeh+iOrder)) / ((real) (nSubSlab*OrderTimeh));
	// Interpolate U to sub-slab nodes -> SU
	ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Uj, OrderTimeh+1, xt, SUj, NULL));
	if (ierr != xf_OK) return ierr;
	
	// Set SUprev
	SUprev = ((iSubSlab==0) ? Uprev : SUj[0]);

        if (UseGCL){
          // Temporally interpolate GCLjTemp (OrderTime) -> GCLj (OrderTimeh, possibly on sub-slabs)
          ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, GCLjTemp, OrderTimeh+1, xt, GCLj, NULL));
          if (ierr != xf_OK) return ierr;

          // Set GCLprev
          GCLprev = ((iSubSlab==0) ? GCLprev0: GCLj[0]);
        }        
      }
      else{
	for (iOrder=0; iOrder<=OrderTime; iOrder++) SUj[iOrder] = Uj[iOrder];
	SUprev = Uprev;
        // don't need to set GCL quantities here ...
      }

      // Set initial CFL
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart));
      if (ierr != xf_OK) return ierr;
      

      // Solve for adjoint states on current time (sub) slab
      ierr = xf_Error(xf_DGinTimeStepAdjoint(All, TimeSchemeh, Time, TimeStep,
					     iSlab*nSubSlab+iSubSlab, 
					     nSlab*nSubSlab, SolverData, SUj,
					     xfe_True, xfe_False, nPsi, Psi, 
					     Psinext, Psii, &Redo));
      if (ierr != xf_OK) return ierr;

      if (Redo) break; // error occured, try again with a lower time step

      if (UseGCL){
        // Solve for GCL adjoint on current time (sub) slab
        xf_printf("Solving GCL adjoint system.\n");
        ierr = xf_Error(xf_DGinTimeStepGCLAdjoint(All, TimeSchemeh, Time, TimeStep,
                                                  iSlab*nSubSlab+iSubSlab, 
                                                  nSlab*nSubSlab, SolverData, SUj,
                                                  nPsi, Psii, PsiGCLnext, PsiGCLi));
        if (ierr != xf_OK) return ierr; 
      }

      /*** Error Estimation ***/

      if (UErrEstOn){

	if (((UErrEstIterative) || (UErrEstUseReconstruct)) && (OrderTimeIh != OrderTimeh)){
	  // project SUj (i.e. Uj) temporally if iterating or time-reconstructing adjoint
	  for (iOrder=0; iOrder<=OrderTimeIh; iOrder++) // node positions
	    xt[iOrder] = ((real) (iSubSlab*OrderTimeIh+iOrder)) / ((real) (nSubSlab*OrderTimeIh));
	  // Interpolate U to nodes -> SU
	  ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Uj, OrderTimeIh+1, xt, SUj, NULL));
	  if (ierr != xf_OK) return ierr;
          // Same with GCL
          if (UseGCL){
            ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, GCLjTemp, OrderTimeIh+1, xt, GCLj, NULL));
            if (ierr != xf_OK) return ierr;
          }   
	}
	if ((ErrEstUseReconstruct || UErrEstIterative) && (SpaceOrderIncrement > 0)){
	  // project SUj spatially if reconstructing or iterating adjoint
	  for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
	    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, SUj[iOrder], NULL, 
								   xfe_BasisLast, SpaceOrderIncrement));
	    if (ierr != xf_OK) return ierr;
	  }
	  // also project SUprev = Uprev if on SubSlab 0
	  if (iSubSlab == 0){
	    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, SUprev, NULL, 
								   xfe_BasisLast, SpaceOrderIncrement));
	    if (ierr != xf_OK) return ierr;
	  }
          // same with GCL
          if (UseGCL){
            // GCLj
            for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
              ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLj[iOrder], NULL, 
                                                                     xfe_BasisLast, SpaceOrderIncrement));
              if (ierr != xf_OK) return ierr;
            }
            // GCL (note, this vector persists, so here we just make it look like GCLj[0])
            ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, GCL, GCLj[0]));
            if (ierr != xf_OK) return ierr;
            // GCLprev
            if (iSubSlab == 0){
              ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLprev, NULL, 
                                                                     xfe_BasisLast, SpaceOrderIncrement));
              if (ierr != xf_OK) return ierr;
            }
          }
	}

	// Loop over adjoint vectors
	for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
	  if (UErrEstIterative){
	    if (OrderTimeIh != OrderTimeh){
	      // temporally project adjoint
	      ierr = xf_Error(xf_DGTimeInject(OrderTimeh, OrderTimeIh, Psii[iAdjoint]));
	      if (ierr != xf_OK) return ierr;
              // same with GCL
              if (UseGCL){
                ierr = xf_Error(xf_DGTimeInject(OrderTimeh, OrderTimeIh, PsiGCLi[iAdjoint]));
                if (ierr != xf_OK) return ierr;
              }
	    }
	    if (SpaceOrderIncrement > 0){
	      // spatially project Psii(p)->AdjiFine(p+) and Psinext(p)->Psinext(p+)
	      for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
		ierr = xf_Error(xf_FindSimilarVector(All, Psii[iAdjoint][iOrder], "Irrelevant", 
						     xfe_True, xfe_False, NULL, 
						     AdjiFine+iOrder, NULL));
		if (ierr != xf_OK) return ierr;
		ierr = xf_Error(xf_SetVector(Psii[iAdjoint][iOrder], xfe_Set, AdjiFine[iOrder]));
		if (ierr != xf_OK) return ierr;
		ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, AdjiFine[iOrder], NULL, 
								       xfe_BasisLast, SpaceOrderIncrement));
		if (ierr != xf_OK) return ierr;
                // TODO: same with GCL
                if (UseGCL) return xf_Error(xf_NOT_SUPPORTED);
	      }
	      if (!UErrEstIterativeGS){ // project Psinext if not Gauss-Seidel
		ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Psinext[iAdjoint], NULL,
								       xfe_BasisLast, SpaceOrderIncrement));
		if (ierr != xf_OK) return ierr;
                // TODO: same with GCL
                if (UseGCL) return xf_Error(xf_NOT_SUPPORTED);
	      }
	    }
	    else{ // no spatial projection, AdjiFine = Psii
	      for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++)
		AdjiFine[iOrder] = Psii[iAdjoint][iOrder];
              // TODO: same with GCL
              if (UseGCL) return xf_Error(xf_NOT_SUPPORTED);
	    }
	    // Which Psinext vector to use for the iterative solver?
	    PsinextIterative = ( (UErrEstIterativeGS) ? PsinextFine+iAdjoint : Psinext+iAdjoint);
            // TODO: same with GCL
            if (UseGCL) return xf_Error(xf_NOT_SUPPORTED);

	    // take iterations of fine solver for adjoint
	    ierr = xf_Error(xf_DGinTimeStepAdjoint(All, TimeSchemeIh, Time, TimeStep,
						   iSlab*nSubSlab+iSubSlab, nSlab*nSubSlab, 
						   SolverData, SUj, xfe_False, xfe_True, 
						   1, Psi, PsinextIterative, &AdjiFine, &btemp));
	    if (ierr != xf_OK) return ierr;

            // TODO: Take iterations of GCL adjoint solver too
            if (UseGCL) return xf_Error(xf_NOT_SUPPORTED);
	    
	    // store first adjoint vector in PsinextFine
	    if (UErrEstIterativeGS){
	      if (PsinextFine[iAdjoint] != NULL){
		ierr = xf_Error(xf_DestroyVector(PsinextFine[iAdjoint], xfe_True));
		if (ierr != xf_OK) return ierr;
	      }
	      ierr = xf_Error(xf_CreateVector(PsinextFine+iAdjoint));
	      if (ierr != xf_OK) return ierr;
	      ierr = xf_CopyVector(All->Mesh, AdjiFine[0], PsinextFine[iAdjoint]);
	      if (ierr != xf_OK) return ierr;
              // TODO: same with GCL
              if (UseGCL) return xf_Error(xf_NOT_SUPPORTED);
	    }
            
	    if ((!UErrEstIterativeGS) && (SpaceOrderIncrement > 0)){
	      // bring Psinext back to original order if not GS
	      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Psinext[iAdjoint],
								     NULL, xfe_BasisLast, -SpaceOrderIncrement));
	      if (ierr != xf_OK) return ierr;
              // TODO: same with GCL
              if (UseGCL) return xf_Error(xf_NOT_SUPPORTED);
	    }
	  } // end if iterative
	  else if ((UErrEstUseReconstruct) || (ErrEstUseReconstruct)){
	    // reconstruct in time and/or space
	    if (UErrEstUseReconstruct){
	      // temporally reconstruct adjoint
	      if (OrderTimeIh == OrderTimeh) return xf_Error(xf_CODE_LOGIC_ERROR);
	      ierr = xf_Error(xf_DGTimeReconstruct(All, OrderTimeh, OrderTimeIh, Psii[iAdjoint], 
						   Psinext[iAdjoint], xfe_True));
	      if (ierr != xf_OK) return ierr;
              // same with GCL
              if (UseGCL){
                ierr = xf_Error(xf_DGTimeReconstruct(All, OrderTimeh, OrderTimeIh, PsiGCLi[iAdjoint], 
                                                     PsiGCLnext[iAdjoint], xfe_True));
                if (ierr != xf_OK) return ierr;
              }
	    }
	    if ((ErrEstUseReconstruct) && (SpaceOrderIncrement > 0)){
	      // spatially reconstruct adjoint
	      for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
		ierr = xf_Error(xf_ReconstructVector(All, SpaceOrderIncrement,
						     Psii[iAdjoint][iOrder], 
						     xfe_True, AdjiFine+iOrder));
		if (ierr != xf_OK) return ierr;
	      }
              // same with GCL
              if (UseGCL){
                for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
                  ierr = xf_Error(xf_ReconstructVector(All, SpaceOrderIncrement,
                                                       PsiGCLi[iAdjoint][iOrder], 
                                                       xfe_True, AdjGCLiFine+iOrder));
                  if (ierr != xf_OK) return ierr;
                }
              }
	    }
	    else{
	      // otherwise, set AdjiFine = Psii[iAdjoint] (pointer-wise)
	      for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++)
		AdjiFine[iOrder] = Psii[iAdjoint][iOrder];
              // same with GCL
              if (UseGCL){
                for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++)
                  AdjGCLiFine[iOrder] = PsiGCLi[iAdjoint][iOrder];
              }
	    }
	  }
	  else{
	    // otherwise, set AdjiFine = Psii[iAdjoint] (pointer-wise)
	    for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++)
	      AdjiFine[iOrder] = Psii[iAdjoint][iOrder];
            // same with GCL
            if (UseGCL){
              for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++)
                AdjGCLiFine[iOrder] = PsiGCLi[iAdjoint][iOrder];
            }
	  }
        
	  // compute error estimate contribution from sub-slab
	  rvec = ((UErrEstAnisoMeasure == xfe_SpaceTimeAnisoProj) ? NULL : SpaceTimePref);
	  ierr = xf_Error(xf_ErrEstDGTimeSlab(All, TimeSchemeIh, Time, TimeStep, SolverData,
					      SUj, SUprev, GCLj, GCLprev, AdjiFine, AdjGCLiFine, rvec, 
                                              OutputError+iAdjoint, OutputErrorSpace+iAdjoint, 
                                              OutputErrorTime+iAdjoint, ErrIndElem[iAdjoint], 
                                              ErrIndElemS[iAdjoint], ErrIndElemT[iAdjoint]));
	  if (ierr != xf_OK) return ierr;
        
	  // destroy AdjiFine if reconstructed
	  if ((UErrEstIterative || ErrEstUseReconstruct) && (SpaceOrderIncrement > 0)){
	    for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
	      ierr = xf_Error(xf_DestroyVector(AdjiFine[iOrder], xfe_True));
	      if (ierr != xf_OK) return ierr;
	    } // iOrder
            // same with GCL
            if (UseGCL)
              for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
                ierr = xf_Error(xf_DestroyVector(AdjGCLiFine[iOrder], xfe_True));
                if (ierr != xf_OK) return ierr;
              } // iOrder
	  }
	  
	} // iAdjoint

	// project SUj back to original order so it can be used on next sub-slab
	if ((ErrEstUseReconstruct || UErrEstIterative) && (SpaceOrderIncrement > 0) 
	    && (iSubSlab > 0)){
	  for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
	    ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, SUj[iOrder], NULL, 
						    xfe_BasisLast, -SpaceOrderIncrement));
	    if (ierr != xf_OK) return ierr;
	  }
          // same with GCL vector
          if (UseGCL){
            // GCL
            ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCL, NULL, 
                                                                   xfe_BasisLast, -SpaceOrderIncrement));
            if (ierr != xf_OK) return ierr;
            // GCLj
            for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
              ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, GCLj[iOrder], NULL, 
                                                                     xfe_BasisLast, -SpaceOrderIncrement));
              if (ierr != xf_OK) return ierr;
            }
          } // UseGCL
	}
      	
      } // end if UErrEstOn
      
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
	/* Destroy PsiExtra if it exists and spatial refinement changed */
	if ((PsiExtra[iAdjoint] != NULL) && (SpatialDyno)){
	  ierr = xf_Error(xf_DestroyVector(PsiExtra[iAdjoint], xfe_True));
	  if (ierr != xf_OK) return ierr;
	  PsiExtra[iAdjoint] = NULL;
	}
        /* Set PsiExtra = vector similar to Psii */
	if (PsiExtra[iAdjoint] == NULL){
	  ierr = xf_Error(xf_FindSimilarVector(All, Psii[iAdjoint][0], "PsiExtra", xfe_True, xfe_False,  
					       NULL, PsiExtra+iAdjoint, NULL));
	  if (ierr != xf_OK) return ierr;
	}

        if (UseGCL){
          // same with GCL
          if ((PsiGCLExtra[iAdjoint] != NULL) && (SpatialDyno)){
            ierr = xf_Error(xf_DestroyVector(PsiGCLExtra[iAdjoint], xfe_True));
            if (ierr != xf_OK) return ierr;
            PsiGCLExtra[iAdjoint] = NULL;
          }
          if (PsiGCLExtra[iAdjoint] == NULL){
            ierr = xf_Error(xf_FindSimilarVector(All, PsiGCLi[iAdjoint][0], "PsiGCLExtra", xfe_True, xfe_False,  
                                                 NULL, PsiGCLExtra+iAdjoint, NULL));
            if (ierr != xf_OK) return ierr;
          }          
        }
      }
      
      // set Psinext = Psii[0] = the left node from the previous time slab
      // Note, Psii could have been modified above, but Psii[0] remains unchanged
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
	// set by way of PsiExtra as Psinext consists only of pointers
	ierr = xf_Error(xf_SetVector(Psii[iAdjoint][0], xfe_Set, 
				     PsiExtra[iAdjoint]));
	if (ierr != xf_OK) return ierr;
	Psinext[iAdjoint] = PsiExtra[iAdjoint];
      }
      if (UseGCL){
        // same with GCL
        for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
          ierr = xf_Error(xf_SetVector(PsiGCLi[iAdjoint][0], xfe_Set,
                                       PsiGCLExtra[iAdjoint]));
          if (ierr != xf_OK) return ierr;
          PsiGCLnext[iAdjoint] = PsiGCLExtra[iAdjoint];
        }
      }
     
      } // iSubSlab

    // Redo time slab if necessary
    if (Redo){
      nSubSlab++; // increase nSubSlab (smaller time step per sub slab)
      Redid++; // increment counter
      iSlab++; // so that time slab remains the same on next loop
      continue;
    }
    else
      nSubSlab = nSubSlab0; // original number of sub slabs
    

    /* Write out Psi to hard disk if at requested interval (all
       adjoints are written to one dataset) */
    if ((SavePrefix != NULL) && ((iSlab % WriteInterval) == 0)){
      sprintf(OutputFile, "%s_Psi%d.data\0", SavePrefix, iSlab);
      ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSetPsi, NULL, OutputFile));
      if (ierr != xf_OK) return ierr;
    }

    
    if (UErrEstOn){

      // keep track of error indicators
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
	
	// modify ErrIndElemS and ErrIndElemT to incorporate output error in ErrIndElem
	ierr = xf_Error(xf_ApplyAniso2ErrInd(All, ErrIndElem[iAdjoint], ErrIndElemS[iAdjoint], 
					     ErrIndElemT[iAdjoint]));
	if (ierr != xf_OK) return ierr;

	// write out error indicators if at requested interval
	if ((UEIWriteInterval > 0) && (iSlab % UEIWriteInterval) == 0){
	  sprintf(OutputFile, "%s_Err_%s%d.data\0", SavePrefix, Psi[iAdjoint]->OutputName, iSlab);
	  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ErrIndicator", 
					      ErrIndElem[iAdjoint], OutputFile));
	  if (ierr != xf_OK) return ierr;
	  sprintf(OutputFile, "%s_ErrS_%s%d.data\0", SavePrefix, Psi[iAdjoint]->OutputName, iSlab);
	  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ErrIndicatorSpace", 
					      ErrIndElemS[iAdjoint], OutputFile));
	  if (ierr != xf_OK) return ierr;
	  sprintf(OutputFile, "%s_ErrT_%s%d.data\0", SavePrefix, Psi[iAdjoint]->OutputName, iSlab);
	  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ErrIndicatorTime", 
					      ErrIndElemT[iAdjoint], OutputFile));
	  if (ierr != xf_OK) return ierr;
	}

	// use absolute values if conservative
	if (UErrEstConservative){
	  ierr = xf_Error(xf_VectorAbs(ErrIndElemS[iAdjoint]));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_VectorAbs(ErrIndElemT[iAdjoint]));
	  if (ierr != xf_OK) return ierr;
	}
	// add ErrIndElemS to running total
	ierr = xf_Error(xf_SetVector(ErrIndElemS[iAdjoint], xfe_Add, ErrIndElemSTot[iAdjoint]));
	if (ierr != xf_OK) return ierr;
	// add ErrIndElemT to running total
	ierr = xf_Error(xf_VectorNorm(ErrIndElemT[iAdjoint], 0, ErrIndTime[iSlab]+iAdjoint));
	if (ierr != xf_OK) return ierr;
	if (!UErrEstConservative) ErrIndTime[iSlab][iAdjoint] = fabs(ErrIndTime[iSlab][iAdjoint]);
      } // iAdjoint

      // destroy SUj
      for (iOrder = 0; iOrder <= OrderTimeIh; iOrder++){
	ierr = xf_Error(xf_DestroyVector(SUj[iOrder], xfe_True));
	if (ierr != xf_OK) return ierr;
      }

    }

    // destroy DataSet and DataSetPrev (and the vectors included therein)
    ierr = xf_Error(xf_DestroyDataSet(DataSet)); // contains Uj
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyDataSet(DataSetPrev)); // contains Uprev
    if (ierr != xf_OK) return ierr;

  } // iSlab

  // print message if redid any time slabs
  if (Redid > 0) xf_printf("Redid (i.e. split) %d adjoint time slabs for robustness.\n", Redid);
    
  /* Following the loop, set the Psi vector to: [R(i)_U(0)]^T*Psii(i)
     = -M*Psii(0), which is the sensitivity of the output to the
     initial condition, U0 -- see description in xf_Solver.h.  Note,
     Psii(0) is already stored in Psinext, but there may be an order
     difference, so project. */
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    // Make Psi look like Psinext
    ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, 
						   Psi[iAdjoint], Psinext[iAdjoint]));
    if (ierr != xf_OK) return ierr;
    // copy over data
    ierr = xf_Error(xf_SetVector(Psinext[iAdjoint], xfe_Set, Psi[iAdjoint]));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_MultMassMatrix(All, -1.0, Psi[iAdjoint]));
    if (ierr != xf_OK) return ierr;
  }
  if (UseGCL){
    // same with GCL
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      // Make Psi look like Psinext
      ierr = xf_Error(xf_ProjectVectorInPlace_Vector(All->Mesh, All->DataSet, 
                                                     PsiGCL[iAdjoint], PsiGCLnext[iAdjoint]));
      if (ierr != xf_OK) return ierr;
      // copy over data
      ierr = xf_Error(xf_SetVector(PsiGCLnext[iAdjoint], xfe_Set, PsiGCL[iAdjoint]));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_MultMassMatrix(All, -1.0, PsiGCL[iAdjoint]));
      if (ierr != xf_OK) return ierr;
    }
  }


  if (UErrEstOn){
    if (!UErrEstConservative){ // take absolute values at very end
      for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
	ierr = xf_Error(xf_VectorAbs(ErrIndElemSTot[iAdjoint]));
	if (ierr != xf_OK) return ierr;
      }
    }

    // Write out total (summed) spatially-localized error indicators
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      sprintf(OutputFile, "%s_ErrTot_%s.data\0", SavePrefix, Psi[iAdjoint]->OutputName);
      ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ErrIndicatorSpace", 
					  ErrIndElemSTot[iAdjoint], OutputFile));
      if (ierr != xf_OK) return ierr;
    }
    // Write out table of temporally-localized error indicators
    sprintf(OutputFile, "%s_TemporalError.txt\0", SavePrefix);
    ierr = xf_Error(xf_WriteTemporalError(nSlab, nPsi, Psi, NULL, ErrIndTime, sdofTime, OutputFile));
    if (ierr != xf_OK) return ierr;
    // store calculated error estimates and print out values
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      ierr = xf_Error(xf_FindOutput(All->EqnSet, Psi[iAdjoint]->OutputName, &Output));
      if (ierr != xf_OK) return ierr;
      Output->ErrEst = OutputError[iAdjoint];
      if (Verbosity != xfe_VerbosityLow){
        // Calculate for printing
        ierr = xf_Error(xf_VectorNorm(ErrIndElemSTot[iAdjoint], 1, &sumSpace));
        if (ierr != xf_OK) return ierr;
        for (iSlab=0,sumTime=0.; iSlab<nSlab; iSlab++) sumTime += ErrIndTime[iSlab][iAdjoint];
	xf_printf("Output = %s  SpatialError = %.10E  TemporalError = %.10E\n", 
		  Psi[iAdjoint]->OutputName, OutputErrorSpace[iAdjoint], 
		  OutputErrorTime[iAdjoint]);
	xf_printf("Output = %s  Value = %.10E  ErrEst = %.10E\n  sumSpace = %.10E  sumTime = %.10E\n", 
		  Psi[iAdjoint]->OutputName, Output->Value, Output->ErrEst, sumSpace, sumTime);
      }
    }
  }
   
  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;
  
  // Destroy DataSetPsi, but not the contained Psii vectors (which are in All)
  for (D=DataSetPsi->Head; D != NULL; D=D->Next) D->Data = NULL;
  ierr = xf_Error(xf_DestroyDataSet(DataSetPsi));
  if (ierr != xf_OK) return ierr;

  // destroy vectors in error estimation
  if (UErrEstOn){
    // first destroy PsinextFine
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      if (PsinextFine[iAdjoint] != NULL){
	ierr = xf_Error(xf_DestroyVector(PsinextFine[iAdjoint], xfe_True));
	if (ierr != xf_OK) return ierr;
      }
    }
    xf_Release( (void *) PsinextFine);
    // destroy ErrIndElem (not ErrIndElemSTot, which is part of All->DataSet)
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      ierr = xf_Error(xf_DestroyVector(ErrIndElem[iAdjoint], xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    // destroy ErrIndElemS
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      ierr = xf_Error(xf_DestroyVector(ErrIndElemS[iAdjoint], xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    // destroy ErrIndElemT
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
      ierr = xf_Error(xf_DestroyVector(ErrIndElemT[iAdjoint], xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    // destroy SpaceTimePref, which is also not part of All->DataSet
    ierr = xf_Error(xf_DestroyVector(SpaceTimePref, xfe_True));
    if (ierr != xf_OK) return ierr;
    // same with ETemp
    ierr = xf_Error(xf_DestroyVector(ETemp, xfe_True));
    if (ierr != xf_OK) return ierr;

    xf_Release( (void *) ErrIndElemSTot);
    xf_Release( (void *) ErrIndElem);
    xf_Release( (void *) ErrIndElemS);
    xf_Release( (void *) ErrIndElemT);
  }

  // Destroy PsiExtra
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++)
    if (PsiExtra[iAdjoint] != NULL){
      ierr = xf_Error(xf_DestroyVector(PsiExtra[iAdjoint], xfe_True));
      if (ierr != xf_OK) return ierr;
    }
  xf_Release( (void *) PsiExtra);

  // Destroy PsiGCLExtra
  if (UseGCL){
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++)
      if (PsiGCLExtra[iAdjoint] != NULL){
        ierr = xf_Error(xf_DestroyVector(PsiGCLExtra[iAdjoint], xfe_True));
        if (ierr != xf_OK) return ierr;
      }
    xf_Release( (void *) PsiGCLExtra);
  }


    
  // release memory
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) xf_Release( (void *) Psii[iAdjoint]);
  xf_Release( (void *) Psii);
  xf_Release( (void *) Psinext);
  if (UseGCL){
    for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++) xf_Release( (void *) PsiGCLi[iAdjoint]);
    xf_Release( (void *) PsiGCLi);
    xf_Release( (void *) PsiGCL);
    xf_Release( (void *) PsiGCLnext);
  }
  xf_Release( (void *) Uj);
  xf_Release( (void *) SUj);
  xf_Release( (void *) AdjiFine);
  xf_Release( (void *) AdjGCLiFine);
  xf_Release( (void *) OutputError);
  xf_Release( (void *) sdofTime);
  xf_Release( (void *) GCLj);
  xf_Release( (void *) GCLjTemp);
  xf_Release2( (void **) ErrIndTime);


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ApplyTimeSchemeAdapt_DG
static int
xf_ApplyTimeSchemeAdapt_DG(xf_All *All, const char *SavePrefix,
			   enum xfe_AdaptOnType AdaptOn,
			   xf_TimeHistData *TimeHistData)
{
  int ierr, nSlab, iSlab, iOrder, OrderTime, OrderTimeh;
  int SpaceOrderIncrement = 0;
  int TimeOrderIncrement = 0;
  int *sdofTime = NULL;
  enum xfe_TimeSchemeType TimeScheme, TimeSchemeh;
  enum xfe_SpaceTimeAnisoType UErrEstAnisoMeasure;
  enum xfe_Verbosity Verbosity;
  const char *AdaptOnName = NULL;
  char Title[xf_MAXSTRLEN];
  char PreHeader[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  real Time, TimeStep;
  real sumSpace, sumTime;
  real **ErrIndTime = NULL;
  xf_Vector **Uj = NULL;
  xf_Vector *Uprev = NULL;
  xf_Vector *ErrIndElemSTot = NULL;
  xf_Vector *ErrIndElemS    = NULL;
  xf_Vector *ErrIndElemT   = NULL;
  xf_Vector *SpaceTimePref = NULL;
  xf_Vector *rvec = NULL;
  xf_Data *D;
  xf_DataSet *DataSet = NULL;
  xf_DataSet *DataSetPrev = NULL;
  xf_SolverData *SolverData;

  // set convenient variables
  nSlab            = TimeHistData->nTime; // # of time slabs

  // need SavePrefix
  if (SavePrefix == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // Time scheme
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "TimeScheme", 
				     xfe_TimeSchemeName, (int ) xfe_TimeSchemeLast, 
				     (int *) &TimeScheme));
  if (ierr != xf_OK) return ierr;

  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;

  // What is the temporal order increment for error estimation/adaptation?
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "UErrEstOrderIncrement", 
				    &TimeOrderIncrement));
  if (ierr != xf_OK) return ierr;

  // fine time scheme
  TimeSchemeh = TimeScheme;
  if (TimeOrderIncrement == 1){
    if (TimeScheme != xfe_TimeSchemeDG1) return xf_Error(xf_NOT_SUPPORTED);
    TimeSchemeh = xfe_TimeSchemeDG2;
  }
  else if (TimeOrderIncrement != 0) return xf_Error(xf_NOT_SUPPORTED);

  // determine fine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeSchemeh, &OrderTimeh));
  if (ierr != xf_OK) return ierr;


  // Spatial order increment
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "ErrEstOrderIncrement", 
				    &SpaceOrderIncrement));
  if (ierr != xf_OK) return ierr;

  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
				     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
				     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;

  // Name of adaptive indicator
  AdaptOnName = xfe_AdaptOnName[AdaptOn];

  // create spatial error indicator vector
  sprintf(Title, "ErrIndSpaceTot_%s", AdaptOnName);
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  xfe_True, &D, 
				&ErrIndElemSTot, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True; // make indicator writeable
  // zero out spatial error indicator (running total, with absolute values)
  ierr = xf_Error(xf_SetZeroVector(ErrIndElemSTot));
  if (ierr != xf_OK) return ierr;
  // will need another error indicator vector, for local storage (not added to dataset)
  sprintf(Title, "ErrIndElemS_%s", AdaptOnName);
  ierr = xf_Error(xf_FindSimilarVector(All, ErrIndElemSTot, Title, xfe_False, 
				       xfe_False,  NULL, &ErrIndElemS, NULL));
  if (ierr != xf_OK) return ierr;

  // Allocate reals for storing temporally-localized adaptive indicator
  ierr = xf_Error(xf_Alloc2( (void ***) &ErrIndTime, nSlab, 1, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  // Allocate ints for storing spatial dofs per time slab
  ierr = xf_Error(xf_Alloc( (void **) &sdofTime, nSlab, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // pull off UErrEstAnisoMeasure
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "UErrEstAnisoMeasure", 
				     xfe_SpaceTimeAnisoName, (int) xfe_SpaceTimeAnisoLast, 
				     (int *) &UErrEstAnisoMeasure));
  if (ierr != xf_OK) return ierr;

  if (((AdaptOn == xfe_AdaptOnInterpol) || (AdaptOn == xfe_AdaptOnResidual)) 
      && (UErrEstAnisoMeasure != xfe_SpaceTimeAnisoJump)){
    xf_printf("Warning, AdaptOn = %s, but AnisoMeasure = %s\n", xfe_AdaptOnName[AdaptOn],
	      xfe_SpaceTimeAnisoName[UErrEstAnisoMeasure]);
    xf_printf("Setting AnisoMeasure = Jump.\n");
    UErrEstAnisoMeasure = xfe_SpaceTimeAnisoJump;
  }

  // create spatial vs. temporal preference vector (not added to dataset)
  sprintf(Title, "SpaceTimePref_State");
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False, xfe_False, NULL, 
				&SpaceTimePref, NULL));
  if (ierr != xf_OK) return ierr;
  // will need another temp error indicator vector
  sprintf(Title, "ErrIndElemT_%s", AdaptOnName);
  ierr = xf_Error(xf_FindSimilarVector(All, ErrIndElemSTot, Title, xfe_False, 
				       xfe_False,  NULL, &ErrIndElemT, NULL));
  if (ierr != xf_OK) return ierr;
  
  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;

  // Allocate vector of state pointers
  ierr = xf_Error(xf_Alloc( (void **) &Uj, OrderTimeh+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;

 
  /*----------------------------------------*/
  /* Begin loop over time slabs (forwards)  */
  /*----------------------------------------*/
  
  for (iSlab=0; iSlab<nSlab; iSlab++){
    
    // break out if user requests a halt
    if (xf_CheckUserHalt(NULL)) break;
    
    // Set Time to start of slab
    Time =      TimeHistData->Time[iSlab];
    TimeStep  = TimeHistData->TimeStep[iSlab];
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
    if (ierr != xf_OK) return ierr;

    // Write out header for this time slab
    sprintf(PreHeader, "%s %.8E %s %8E", "%% DG Adapt: Time slab start =", Time,
	    " Step =", TimeStep);
    ierr = xf_Error(xf_WriteLogHeader(All, PreHeader));
    if (ierr != xf_OK) return ierr;

    // zero quantities out
    ErrIndTime[iSlab][0] = 0.0;
    ierr = xf_Error(xf_SetZeroVector(ErrIndElemS));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_SetZeroVector(ErrIndElemT));
    if (ierr != xf_OK) return ierr;
	
    // read .data from file (multiple vectors in one dataset)
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;
    sprintf(Title, "%s_U%d.data\0", SavePrefix, iSlab+1);
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, Title, DataSet));
    if (ierr != xf_OK) return ierr;
    
    // use data vectors in order of data
    for (iOrder=0, D=DataSet->Head; iOrder<=OrderTime; iOrder++, D=D->Next)	
      Uj[iOrder] = (xf_Vector *) D->Data;
    if (D != NULL) return xf_Error(xf_CODE_LOGIC_ERROR);

    // read in Uprev for error estimation
    ierr = xf_Error(xf_CreateDataSet(&DataSetPrev));
    if (ierr != xf_OK) return ierr;
    sprintf(Title, "%s_U%d.data\0", SavePrefix, iSlab);
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, Title, DataSetPrev));
    if (ierr != xf_OK) return ierr;
    
    // use last, or only in case of (iSlab==0), data vector
    D = DataSetPrev->Head;
    if (iSlab != 0) for (iOrder=0; (iOrder<OrderTime)&&(D!=NULL); iOrder++, D=D->Next);
    if (D == NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
    Uprev = D->Data;

    // calculate spatial degrees of freedom on this slab
    if (sdofTime != NULL){
      ierr = xf_Error(xf_GetVectorDOF(Uj[0], sdofTime+iSlab));
      if (ierr != xf_OK) return ierr;
    }

    
    // Estimate anisotropy (space vs. time) of error
    if (UErrEstAnisoMeasure == xfe_SpaceTimeAnisoJump){
      rvec = ((AdaptOn == xfe_AdaptOnInterpol) ? ErrIndElemS : NULL);
      ierr = xf_Error(xf_ErrEstSpaceTimeAniso(All, TimeScheme, Uj, Uprev, Time, 
					      TimeStep, SpaceTimePref, rvec));
      if (ierr != xf_OK) return ierr;
      if (AdaptOn == xfe_AdaptOnInterpol){
	ierr = xf_Error(xf_SetVector(ErrIndElemS, xfe_Set, ErrIndElemT));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VectorVectorMult(SpaceTimePref, ErrIndElemT));
  	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VectorMult(SpaceTimePref, -1.0));
 	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VectorAdd(SpaceTimePref, 1.0));
 	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_VectorVectorMult(SpaceTimePref, ErrIndElemS));
  	if (ierr != xf_OK) return ierr;
      }
    }

    // residual-based error estimate requires evaluation of residuals
    if (AdaptOn == xfe_AdaptOnResidual){
      // project Uj to higher order
      for (iOrder = 0; iOrder <= OrderTime; iOrder++){
	ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Uj[iOrder], NULL, 
							       xfe_BasisLast, SpaceOrderIncrement));
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_ProjectVectorInPlace_OrderIncrement(All->Mesh, NULL, Uprev, NULL, 
							     xfe_BasisLast, SpaceOrderIncrement));
      if (ierr != xf_OK) return ierr; 
      

      if (OrderTimeh != OrderTime){
	// additional vectors for temporal high-order (not added to All)
	for (iOrder = OrderTime+1; iOrder <= OrderTimeh; iOrder++){
	  sprintf(Title, "HOState_%d", iOrder);
	  // create required state vectors
	  ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], Title, xfe_True, xfe_False,
					       NULL, Uj + iOrder, NULL));
	  if (ierr != xf_OK) return ierr;
	}
	
	// inject into OrderTimeh
	ierr = xf_Error(xf_DGTimeInject(OrderTime, OrderTimeh, Uj));
	if (ierr != xf_OK) return ierr;
      }
      
      // residual error estimate
      ierr = xf_Error(xf_ErrEstDGTimeSlab(All, TimeScheme, Time, TimeStep, SolverData, Uj, 
					  Uprev, NULL, NULL, NULL, NULL, SpaceTimePref, NULL, NULL, NULL, 
					  NULL, ErrIndElemS, ErrIndElemT));
      if (ierr != xf_OK) return ierr;

      // destroy extra U vectors
      for (iOrder = OrderTime+1; iOrder <= OrderTimeh; iOrder++){
	ierr = xf_Error(xf_DestroyVector(Uj[iOrder], xfe_True));
	if (ierr != xf_OK) return ierr;
      }
    }

    // add ErrIndElemS to running total
    ierr = xf_Error(xf_VectorAbs(ErrIndElemS)); // always use conservative measure
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_SetVector(ErrIndElemS, xfe_Add, ErrIndElemSTot));
    if (ierr != xf_OK) return ierr;

    // add ErrIndElemT to running total
    ierr = xf_Error(xf_VectorAbs(ErrIndElemT));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_VectorNorm(ErrIndElemT, 0, ErrIndTime[iSlab]));
    if (ierr != xf_OK) return ierr;

    // destroy DataSet and DataSetPrev (and the vectors included therein)
    ierr = xf_Error(xf_DestroyDataSet(DataSet)); // contains Uj
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyDataSet(DataSetPrev)); // contains Uprev
    if (ierr != xf_OK) return ierr;

  } // iSlab

  // Write out total (summed) spatially-localized error indicators
  sprintf(OutputFile, "%s_ErrTot_%s.data\0", SavePrefix, AdaptOnName);
  ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "ErrIndicatorSpace", 
				      ErrIndElemSTot, OutputFile));
  if (ierr != xf_OK) return ierr;

  // Write out table of temporally-localized error indicators
  sprintf(OutputFile, "%s_TemporalError.txt\0", SavePrefix);
  ierr = xf_Error(xf_WriteTemporalError(nSlab, 1, NULL, AdaptOnName, ErrIndTime, sdofTime, OutputFile));
  if (ierr != xf_OK) return ierr;

  // store calculated error estimates and print out values
  if (Verbosity != xfe_VerbosityLow){
    // Calculate for printing
    ierr = xf_Error(xf_VectorNorm(ErrIndElemSTot, 1, &sumSpace));
    if (ierr != xf_OK) return ierr;
    for (iSlab=0,sumTime=0.; iSlab<nSlab; iSlab++) sumTime += ErrIndTime[iSlab][0];
    xf_printf("AdaptOn = %s  sumSpace = %.10E  sumTime = %.10E\n", 
	      AdaptOnName, sumSpace, sumTime);
  }
   
  // destroy ErrIndElemS (not ErrIndElemSTot, which is part of All->DataSet)
  ierr = xf_Error(xf_DestroyVector(ErrIndElemS, xfe_True));
  if (ierr != xf_OK) return ierr;
  // destroy ErrIndElemT
  ierr = xf_Error(xf_DestroyVector(ErrIndElemT, xfe_True));
  if (ierr != xf_OK) return ierr;

  // destroy SpaceTimePref, which is also not part of All->DataSet
  ierr = xf_Error(xf_DestroyVector(SpaceTimePref, xfe_True));
  if (ierr != xf_OK) return ierr;

  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) Uj);
  xf_Release( (void *) sdofTime);
  xf_Release2((void **) ErrIndTime);

  return xf_OK;
}

