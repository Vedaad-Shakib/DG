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
 FILE:  xf_SolverROBST.c
 
 This file contains functions for the ROBST type of solver
 
 */

#include "xf.h"
#include "xf_AllStruct.h"
#include "xf_SolverStruct.h"
#include "xf_Penalty.h"
#include "xf_LineSearch.h"
#include "xf_LinearSolver.h"
#include "xf_Memory.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Math.h"
#include "xf_Residual.h"
#include "xf_Solver.h"
#include "xf_SolverTools.h"
#include "xf_Log.h"
#include "../dyn/CompressibleNS/xf_CompressibleNSStruct.h"


/******************************************************************/
//   FUNCTION Definition: xf_SolveLinearSystemCG
static int
xf_SolveLinearSystemCG(xf_All *All, xf_Vector *U, xf_Vector *dU, 
                       xf_Vector *R, xf_SolverData *SolverData, 
                       real LinearTolerance)
{
  int ierr, nIter, it, i, n;
  enum xfe_LinearStatusType Status;
  enum xfe_Bool PreconditionFlag;
  xf_Vector *V, *W, *TempVector;
  xf_VectorSet *VS;
  xf_LinearSolverData *LinearSolverData;
  
  LinearSolverData = NULL;
  V = W = NULL;
  PreconditionFlag = xfe_False; //for now this is hard coded
  
  Status = xfe_LinearMultiply;
  
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "CG_MaxIter", 
                           &nIter);
  if (ierr != xf_OK) return ierr;
  
  it = 0;
  
  n = 3 + (int)PreconditionFlag;
  
  ierr = xf_Error(xf_CreateVectorSet(n, &VS));
  if (ierr != xf_OK) return ierr;
  
  for (i = 0; i < n; i++){
    ierr = xf_Error(xf_FindSimilarVector(All, U, "Vector", xfe_True, 
                                         xfe_False, NULL, &TempVector, NULL));
    if (ierr != xf_OK) return ierr;
    VS->Vector[i] = (*TempVector);
  }
  
  while (Status != xfe_LinearConverged && it < nIter) {
    ierr = xf_Error(xf_LinearIterCG(All, VS, LinearTolerance, xfe_VerbosityLow, 
                                    PreconditionFlag, xfe_False, 
                                    xfe_False, xfe_False, R, dU, &V, &W, 
                                    &LinearSolverData, &Status));
    if (ierr != xf_OK) return ierr;
    
    if (Status == xfe_LinearMultiply){
      ierr = xf_Error(xf_ApplyHessian(All, U, V, W, SolverData, 
                                      ((it == 0) ? xfe_True : xfe_False)));
      if (ierr != xf_OK) return ierr;
    }
    
    if (Status == xfe_LinearPrecondition){
      ierr = xf_Error(xf_ApplyHessianPrecond(All, U, V, W, SolverData, 
                                             ((it == 0) ? xfe_True : xfe_False)));
      if (ierr != xf_OK) return ierr;
    }
    it++;
  }
  
  xf_printf("LinearCG => Iter: %d |R|_2: %1.10e\n", it, LinearSolverData->rnorm2);
  
  ierr = xf_Error(xf_DestroyVectorSet(VS));
  if (ierr != xf_OK) return ierr;
  
  VS = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculateGradient

static int
xf_CalculateGradient(xf_All *All, xf_Vector *U, xf_Vector *G, 
                     xf_SolverData *SolverData)
{
  int ierr;
  xf_Vector *R, *GP;
  xf_JacobianMatrix *R_U;
  xf_Data *D;
  
  //Set the pointer to the residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_True, 
                                       xfe_True, &D, &R, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True;
  ierr = xf_Error(xf_FindSimilarVector(All, U, "GradPenalty", xfe_True, 
                                       xfe_True, &D, &GP, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True;
  //Find the Jacobian Matrix
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_False, NULL, &R_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  //Calculte residual and Jacobian correspondent to state U
  ierr = xf_CalculateResidual(All, U, R, R_U, SolverData);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VectorNorm(R, 1, &SolverData->ResNorm));
  if (ierr != xf_OK) return ierr;
  
  //computing the gradient without penalty terms
  //G = 2*R*dR/dU
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, R, xfe_Set, xfe_True, 
                                   SolverData, G));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VectorMult(G,2.0));
  if (ierr != xf_OK) return ierr;
  
  //add penalty gradient
  ierr = xf_Error(xf_CalculatePenaltyGradient(All, U, GP));
  if (ierr != xf_OK) return ierr;  
  
  ierr = xf_Error(xf_SetVector(GP, xfe_Add, G));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VectorNorm(G, 2, &SolverData->Gnorm2));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Calculate_F
int
xf_Calculate_F(xf_All *All, xf_Vector *U, xf_SolverData *SolverData, 
               enum xfe_Bool UseResidual)
{
  int ierr;
  real Rnorm_L2;
  enum xfe_Bool found;
  xf_Vector *R;
  
  //Set the pointer to the residual vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_True, 
                                       xfe_True, NULL, &R, &found));
  if (ierr != xf_OK) return ierr;
  
  if(!UseResidual){
    //Calculate Residual
    ierr = xf_CalculateResidual(All, U, R, NULL, SolverData);
    if (ierr != xf_OK) return ierr;
  }  
  
  //Calcuating F
  ierr = xf_Error(xf_VectorNorm(R, 2, &Rnorm_L2));
  if (ierr != xf_OK) return ierr;
  
  SolverData->AugResidual = Rnorm_L2*Rnorm_L2;
  
  //Add penalty term F += P(U)
  ierr = xf_Error(xf_CalculatePenalty(All, U, &SolverData->PenaltyFcn));
  if (ierr != xf_OK) return ierr;
  
  SolverData->AugResidual += SolverData->PenaltyFcn;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Calculate_dF
int
xf_Calculate_dF(xf_All *All, xf_Vector *U, xf_Vector *P, real *dF, 
                xf_SolverData *SolverData, enum xfe_Bool UseGradient,
                enum xfe_Bool Normalize)
{
  int ierr;
  real Pnorm, Gnorm;
  enum xfe_Bool found;
  xf_Vector *G;
  
  //Check inputs
  if (P == NULL)
    return xf_INPUT_ERROR;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Gradient", xfe_True, 
                                       xfe_True, NULL, &G, &found));
  if (ierr != xf_OK) return ierr;
  
  if (!found)
    return xf_CODE_LOGIC_ERROR;
  
  //Norm of P
  ierr = xf_Error(xf_VectorNorm(P, 1, &Pnorm));
  if (ierr != xf_OK) return ierr;
  
  
  if(!UseGradient){
    //Calculate gradient
    ierr = xf_Error(xf_CalculateGradient(All, U, G, SolverData));
    if (ierr != xf_OK) return ierr;
  }
  
  //Calculating dF
  ierr = xf_Error(xf_VectorDot(P, G, &(*dF)));
  if (ierr != xf_OK) return ierr;
  
  (*dF) = (*dF)/Pnorm;
  if (Normalize){
    ierr = xf_Error(xf_VectorNorm(G, 1, &Gnorm));
    if (ierr != xf_OK) return ierr;
    
    (*dF) = (*dF)/Gnorm;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ComputeSearchDirection
static int 
xf_ComputeSearchDirection(xf_All *All, xf_Vector *U, xf_Vector *P, 
                          xf_SolverData *SolverData)
{
  int ierr;
  real LinearTolerance;
  xf_Vector *R, *G, *dU;
  xf_JacobianMatrix *R_U;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Gradient", xfe_True, 
                                       xfe_True, NULL, &G, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_True, 
                                       xfe_True, NULL, &R, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_True, 
                                       xfe_True, NULL, &dU, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_False, NULL, &R_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  //CalculateGradient updates the Jacobian and the residual
  ierr = xf_Error(xf_CalculateGradient(All, U, G, SolverData));
  if (ierr != xf_OK) return ierr;
  
  LinearTolerance = MEPS;//SolverData->Gnorm2/1.0e6; //hard coded
  
  //initial estimate: dU = -G
  ierr = xf_Error(xf_SetVector(G, xfe_Neg, dU));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SolveLinearSystemCG(All, U, dU, G, SolverData, 
                                         LinearTolerance));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SetVector(dU, xfe_Set, P));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DumpR2PlusP
static int 
xf_DumpR2PlusP(xf_All *All, xf_Vector *U, xf_Vector *P, 
               xf_SolverData *SolverData)
{
  int ierr, i, SliceNPoints, Iter0;
  real MaxAlpha, DAlpha, Alpha, Pnorm;
  char fname[xf_MAXSTRLEN];
  xf_Vector *Utrial;
  FILE *fid;
  
  //getting parameters
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, 
                            "DumpR2PlusPMaxAlpha", &MaxAlpha);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, 
                           "DumpR2PlusPN", &SliceNPoints);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", 
                           &Iter0);
  if (ierr != xf_OK) return ierr;
  
  DAlpha = MaxAlpha/(SliceNPoints-1); //alpha increment
  
  sprintf(fname,"R2PlusP_it%d.txt",SolverData->iIter+Iter0);
  fid = fopen(fname,"w");
  
  //get norm of P
  ierr = xf_Error(xf_VectorNorm(P, 1, &Pnorm));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VectorMult(P, 1.0/Pnorm));
  if (ierr != xf_OK) return ierr;
  
  //Getting a trial state vector
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Utrial",
                                       xfe_True, xfe_False, 
                                       NULL, &Utrial, NULL));
  if (ierr != xf_OK) return ierr;
  
  
  for (i = 0; i < SliceNPoints; i++){
    //update state
    Alpha = i*DAlpha;
    
    ierr = xf_Error(xf_SetVector(U, xfe_Set, Utrial));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_VectorMultSet(P, Alpha, xfe_Add, Utrial));
    if (ierr != xf_OK) return ierr;
    
    //Calculate R2+P and print: Alpha R2 P R2+P
    ierr = xf_Error(xf_Calculate_F(All, Utrial, SolverData, xfe_False));
    if (ierr == xf_NON_PHYSICAL){
      fprintf(fid,"%1.16e nan nan nan\n",Alpha);
    }
    else if (ierr == xf_OK){
      fprintf(fid,"%1.16e %1.16e %1.16e %1.16e\n", Alpha, 
              SolverData->AugResidual-SolverData->PenaltyFcn,
              SolverData->PenaltyFcn, SolverData->AugResidual);
    }
    else
      return ierr;
  }
  
  //rescaling P
  ierr = xf_Error(xf_VectorMult(P, Pnorm));
  if (ierr != xf_OK) return ierr;
  
  //cleaning up
  ierr = xf_Error(xf_DestroyVector(Utrial, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  fclose(fid);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SolveNonlinearSystem_ROBST
int
xf_SolveNonlinearSystem_ROBST(xf_All *All, enum xfe_Bool LinearFlag, 
                              xf_Vector *S, xf_Vector *U)
{
  int ierr, major_it, *IParam, nIterNonlinear, Iter0;
  real ResidualTolerance, dF, *RParam, Gnorm_prev;
  real PenaltyFcnFactor, PenaltyReduceFactor;
  enum xfe_Bool Converged, WriteR2PlusP, ReducePenalty;
  xf_Vector *R, *G, *dU, *P;
  xf_JacobianMatrix *R_U;
  xf_SolverData *SolverData;
  xf_Data *D;
  
  ierr = xf_GetKeyValueBool(All->Param->KeyValue, "DumpR2PlusP", 
                            &WriteR2PlusP);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "nIterNonlinear", 
                           &nIterNonlinear);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_GetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", 
                           &Iter0);
  if (ierr != xf_OK) return ierr;
  
  //pull off penalty reduction factor
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "PenaltyReduceFactor",
                            &PenaltyReduceFactor);
  if (ierr != xf_OK) return ierr;
  
  // pull off residual tolerance
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "ResidualTolerance", 
                            &ResidualTolerance);
  if (ierr != xf_OK) return ierr;
  
  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_True, NULL, &R_U, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Residual", xfe_True, 
                                       xfe_True, &D, &R, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Gradient", xfe_True, 
                                       xfe_True, &D, &G, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "dU", xfe_True, xfe_True,
                                       NULL, &dU, NULL));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_FindSimilarVector(All, U, "Pk", xfe_True, xfe_True, 
                                       NULL, &P, NULL));
  if (ierr != xf_OK) return ierr;
  
  //setting to zero the search direction
  ierr = xf_Error(xf_SetZeroVector(dU));
  if (ierr != xf_OK) return ierr;
  
  //setting up the parameters
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, 
                                       &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  PenaltyFcnFactor = RParam[xfe_PenaltyFcnFactor];

  major_it = 0;
  SolverData->iIter = Iter0;
  //start
  Converged = xfe_False;
  while (!Converged && SolverData->iIter < nIterNonlinear + Iter0){
    ReducePenalty = xfe_False;
    
    ierr = xf_Error(xf_CalculateGradient(All, U, G, SolverData));
    if (ierr != xf_OK) return ierr;
    Gnorm_prev = SolverData->Gnorm2;
    
    ierr = xf_Error(xf_Calculate_F(All, U, SolverData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_WriteLogEntry(All, SolverData, U));
    if (ierr != xf_OK) return ierr;
    
    while (SolverData->iIter < nIterNonlinear + Iter0) {
      //check for user halt
      if (xf_CheckUserHalt(NULL)) break;
      
      //check for convergence
      if ((SolverData->Gnorm2 <= Gnorm_prev/1.0e1) ||//Hard coded
          (SolverData->ResNorm <= ResidualTolerance)){
        xf_printf("Nonlinear solver converged for PenaltyFcnFactor = %1.2e\n",
                  PenaltyFcnFactor);
        ReducePenalty = xfe_True;
        break;
      }
      //we should go this way...hehe
      ierr = xf_Error(xf_ComputeSearchDirection(All, U, P, SolverData));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Calculate_dF(All, U, P, &dF, SolverData, 
                                      xfe_False, xfe_True));
      if (ierr != xf_OK) return ierr;
      
      //check sign of df
      if (dF >= 0){
        xf_printf("Warning, dU is not a descent direction! Setting P to -G\n");
        ierr = xf_Error(xf_SetVector(G, xfe_Neg, P));
        if (ierr != xf_OK) return ierr;
      }
      
      if (WriteR2PlusP){
        ierr = xf_Error(xf_DumpR2PlusP(All, U, P, SolverData));
        if (ierr != xf_OK) return ierr;
      }
      
      ierr = xf_Error(xf_LineSearch(All, U, P, SolverData));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_CalculateGradient(All, U, G, SolverData));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_Calculate_F(All, U, SolverData, xfe_True));
      if (ierr != xf_OK) return ierr;
      
      SolverData->iIter++;
      ierr = xf_Error(xf_WriteLogEntry(All, SolverData, U));
      if (ierr != xf_OK) return ierr;
    }
    major_it++;
    //check for user halt
    if (xf_CheckUserHalt(NULL)) break;
    
    if ((SolverData->Gnorm2 <= ResidualTolerance || 
         SolverData->ResNorm <= ResidualTolerance)&& 
        PenaltyFcnFactor < ResidualTolerance){
      xf_printf("Nonlinear solver converged to tolerance.\n");
      Converged = xfe_True;
    }
    
    if (ReducePenalty){
      //Reduce the Penalty Factor
      PenaltyFcnFactor *= xf_PowInt(PenaltyReduceFactor,major_it);
      ierr = xf_Error(xf_SetKeyValueReal(All->EqnSet->KeyValue, 
                                         "PenaltyFcnFactor", 
                                         PenaltyFcnFactor));
      if (ierr != xf_OK) return ierr;
      xf_printf("Reducing PenaltyFcnFactor to %1.2e\n",PenaltyFcnFactor);
      continue;
    }
    
  }
  
  ierr = xf_SetKeyValueInt(All->Param->KeyValue, "iIterNonlinear", SolverData->iIter);
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;
  
  xf_Release( (void *) RParam);
  xf_Release( (void *) IParam);
  
  return xf_OK;
}
