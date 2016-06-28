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
 FILE:  xf_LineSearch.c
 
 This file contains functions related to line search.
 
 */


#include "xf.h"
#include "xf_AllStruct.h"
#include "xf_Solver.h"
#include "xf_SolverTools.h"
#include "xf_SolverStruct.h"
#include "xf_SolverROBST.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Penalty.h"
#include "xf_LineSearch.h"
#include "xf_Math.h"
#include "xf_Param.h"
#include "../dyn/CompressibleNS/xf_CompressibleNSStruct.h"

/******************************************************************/
//   FUNCTION Definition: xf_UpdateStateWithP
static int 
xf_UpdateStateWithP(xf_All *All, xf_Vector *U, xf_Vector *P, 
                    real atrial, xf_Vector *Utrial, real *Ftrial, 
                    real *dFtrial)
{
  int ierr;
  xf_SolverData *SolverData;
  
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;
  
  //Update Utrial
  ierr = xf_Error(xf_SetVector(U, xfe_Set, Utrial));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VectorMultSet(P, atrial, xfe_Add, Utrial));
  if (ierr != xf_OK) return ierr;
  
  if (Ftrial != NULL){
    //calculate f for Utrial
    ierr = xf_Error(xf_Calculate_F(All, Utrial, SolverData, xfe_False));
    if (ierr != xf_OK) return ierr;  
    Ftrial[0] = SolverData->AugResidual;
  }
  
  if (dFtrial != NULL){
    //calculate df for Utrial
    ierr = xf_Error(xf_Calculate_dF(All, Utrial, P, dFtrial, SolverData, xfe_False, xfe_True));
    if (ierr != xf_OK) return ierr;  
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_LineSearchTrialPoint
static int 
xf_LineSearchTrialPoint(real Lo[3], real Hi[3], real *aTrial)
{
  real aL, FL, dFL, aR, FR, dFR, curv;
  
  //ordering the values
  if (Lo[0] < Hi[0]){
    aL   = Lo[0];
    FL   = Lo[1];
    dFL  = Lo[2];
    aR   = Hi[0];
    FR   = Hi[1];
    dFR  = Hi[2];
  }
  else{
    aL   = Hi[0];
    FL   = Hi[1];
    dFL  = Hi[2];
    aR   = Lo[0];
    FR   = Lo[1];
    dFR  = Lo[2];
  }
  
  curv = 2.0*(FR-FL - dFL*(aR-aL))/xf_PowInt(aR-aL,2);
  
  if (fabs(curv) <= MEPS){ 
    xf_printf("Curvature = %1.10e\n",curv/2.0);
    return xf_SINGULAR;
  }
  
  aTrial[0] = -dFL/curv + aL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_LineSearchZoom
static int 
xf_LineSearchZoom(xf_All *All, xf_Vector *U, xf_Vector *P, 
                  real *aStar, xf_SolverData *SolverData, 
                  real dF, real Lo[3], real Hi[3])
{
  int ierr, it, nIter;
  real atrial, Ftrial, dFtrial, F = SolverData->AugResidual;
  xf_Vector *Utrial;
  enum xfe_Bool Update = xfe_False;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Utrial", xfe_True, 
                                       xfe_False, NULL, &Utrial, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "nIterLineSearchZoom", &nIter);
  if (ierr != xf_OK) return ierr;
  
  it = 0;
  while (!Update){
    it++;
    //checks for maximum iterations or an interval that is too-small
    if ((it > nIter) || (fabs(Hi[0]-Lo[0]) <= MEPS)){
      aStar[0] = ((Hi[0] > Lo[0]) ? Hi[0] : Lo[0]);
      Update = xfe_True;
      continue;
    }
    //step 1,2: interpolating and picking a trial point
    //get trial point
    ierr = xf_Error(xf_LineSearchTrialPoint(Lo, Hi, &atrial));
    if (ierr != xf_OK) return ierr;
    //calculate f and df
    ierr = xf_Error(xf_UpdateStateWithP(All, U, P, atrial, Utrial, &Ftrial, &dFtrial));
    if (ierr != xf_OK) return ierr;
    
    //step 3
    if ((Ftrial <= F + xf_mu1*atrial*dF) && (Ftrial <= Lo[1])){
      //step 4.1
      if (fabs(dFtrial) <= -xf_mu2*dF){
        //step 4.2
        aStar[0] = atrial;
        Update = xfe_True;
        continue;
      }
      else {
        //step 4.3
        if (dFtrial*(Hi[0]-Lo[0]) >= 0.0){
          Hi[0] = Lo[0];
          Hi[1] = Lo[1];
          Hi[2] = Lo[2];
        }
        //step 4.4
        Lo[0] = atrial;
        Lo[1] = Ftrial;
        Lo[2] = dFtrial;
        continue;
      }
    }
    else {
      Hi[0] = atrial;
      Hi[1] = Ftrial;
      Hi[2] = dFtrial;
      continue;
    }
  }
  
  //cleaning up
  ierr = xf_Error(xf_DestroyVector(Utrial, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_LineSearch
int 
xf_LineSearch(xf_All *All, xf_Vector *U, xf_Vector *P, 
              xf_SolverData *SolverData)
{
  int ierr;
  real aPrev, aCurr, FPrev, FCurr, dFPrev, dFCurr, IncFac, DecFac, proj, F, dF;
  real aStar, Lo[3], Hi[3], Pnorm;
  enum xfe_Bool Update, found;
  xf_Vector *UPrev, *UCurr, *Gp;
  
  //get norm of P
  ierr = xf_Error(xf_VectorNorm(P, 1, &Pnorm));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VectorMult(P, 1.0/Pnorm));
  if (ierr != xf_OK) return ierr;
  
  //extracting parameters
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, 
                            "UpdateFracAmplification", &IncFac);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, 
                            "UpdateFracReduction", &DecFac);
  if (ierr != xf_OK) return ierr;
  
  /****************************************************************************/
  //now let's check if we need to start small
  ierr = xf_Error(xf_FindSimilarVector(All, U, "GradPenalty", xfe_True, 
                                       xfe_True, NULL, &Gp, &found));
  if (ierr != xf_OK) return ierr;
  if (!found)
    return xf_CODE_LOGIC_ERROR;
  
  //calculate the projection
  ierr = xf_Error(xf_VectorDot(P, Gp, &proj));
  if (ierr != xf_OK) return ierr;
  
  if (proj > 0.0){
    //it's likely that we are going towards a non-physical region
    //start small
    aCurr = xf_mu1;
  }
  else 
    aCurr = Pnorm;
  /****************************************************************************/
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "UPrev", xfe_True, 
                                       xfe_False, NULL, &UPrev, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "UCurr", xfe_True, 
                                       xfe_False, NULL, &UCurr, NULL));
  if (ierr != xf_OK) return ierr;
  
  //Copy state and calculate F and dF
  //Note: UPrev is physical because U is assumed to be physical
  aPrev = 0.0;
  ierr = xf_Error(xf_UpdateStateWithP(All, U, P, aPrev, UPrev, &FPrev, &dFPrev));
  if (ierr != xf_OK) return ierr;
  F = FPrev;
  dF = dFPrev;
  
  //Start loop
  Update = xfe_False;
  while (!Update) {
    //Update UCurr
    ierr = xf_Error(xf_UpdateStateWithP(All, U, P, aCurr, UCurr, &FCurr, &dFCurr));
    if (ierr == xf_NON_PHYSICAL){ //first guess was too large
      aCurr *= DecFac;
      aPrev = 0.0;
      FPrev = F;
      dFPrev = dF;
      continue;
    }
    else if (ierr != xf_OK) return ierr;
    //Check sufficient descent
    if (FCurr > F + xf_mu1*aCurr*dF || FCurr > FPrev){//Does not satisfy
      //Zoom and return
      Lo[0] = aPrev;
      Lo[1] = FPrev;
      Lo[2] = dFPrev;
      Hi[0] = aCurr;
      Hi[1] = FCurr;
      Hi[2] = dFCurr;
      
      ierr = xf_Error(xf_LineSearchZoom(All, U, P, &aStar, SolverData, dF, Lo, Hi));
      if (ierr != xf_OK) return ierr;
      Update = xfe_True;
      continue;
    }
    else {
      //Check for curvature condition
      if (fabs(dFCurr) <= -xf_mu2*dF){
        aStar = aCurr;
        Update = xfe_True;
        continue;
      }
      if (dFCurr >= 0.0){
        //Zoom and return
        Lo[0] = aCurr;
        Lo[1] = FCurr;
        Lo[2] = dFCurr;
        Hi[0] = aPrev;
        Hi[1] = FPrev;
        Hi[2] = dFPrev;
        
        ierr = xf_Error(xf_LineSearchZoom(All, U, P, &aStar, SolverData, dF, Lo, Hi));
        if (ierr != xf_OK) return ierr;
        Update = xfe_True;
        continue;
      }
      else {
        aPrev = aCurr;
        FPrev = FCurr;
        dFPrev = dFCurr;
        //increase alpha and continue
        aCurr *= IncFac;
      }
    }
  }
  
  SolverData->UpdateFrac = aStar;
  
  //updating the state
  ierr = xf_Error(xf_VectorMultSet(P, SolverData->UpdateFrac, xfe_Add, U));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Calculate_F(All, U, SolverData, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  //rescaling P
  ierr = xf_Error(xf_VectorMult(P, Pnorm));
  if (ierr != xf_OK) return ierr;
  
  //cleaning up
  ierr = xf_Error(xf_DestroyVector(UCurr, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyVector(UPrev, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}
