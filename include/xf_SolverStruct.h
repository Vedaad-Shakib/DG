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

#ifndef _xf_SolverStruct_h
#define _xf_SolverStruct_h 1

/*
  FILE:  xf_SolverStruct.h

  This file contains the xflow Solver structures

*/


#include "xf.h"
#include "xf_DataStruct.h"

/* Maximum number of vectors used at once for the unsteady
   multistep methods */
#define xf_MAXMULTISTEP 4

// maximum number of time nodes for DG in time = order+1
#define xf_MAXDGTIMENODE 3

/* Steady or unsteady calculations dictated by this type */
enum xfe_TimeSchemeType{
  xfe_TimeSchemeSteady,
  xfe_TimeSchemeBDF1,
  xfe_TimeSchemeBDF2,
  xfe_TimeSchemeBDF3,
  xfe_TimeSchemeTrapezoidal,
  xfe_TimeSchemeFE,
  xfe_TimeSchemeRK4,
  xfe_TimeSchemeRK2,
  xfe_TimeSchemeDG1,
  xfe_TimeSchemeDG2,
  xfe_TimeSchemeESDIRK3,
  xfe_TimeSchemeESDIRK4,
  xfe_TimeSchemeESDIRK5,
  xfe_TimeSchemeDIRK3,
  xfe_TimeSchemeDIRK4,
  xfe_TimeSchemeLast
};

/* corresponding names */
static char *xfe_TimeSchemeName[xfe_TimeSchemeLast] = {
  "Steady",
  "BDF1",          // Backwards Euler
  "BDF2",          // Second-order backwards difference
  "BDF3",          // Third-order backwards difference
  "Trapezoidal",   // Crank-Nicholson
  "FE",            // Forward Euler
  "RK4",           // 4th order Runge-Kutta
  "RK2",           // 2nd order Runge-Kutta
  "DG1",           // DG in time, r=1
  "DG2",           // DG in time, r=2
  "ESDIRK3",       // ESDIRK3
  "ESDIRK4",       // ESDIRK4
  "ESDIRK5",       // ESDIRK4
  "DIRK3",         // DIRK3
  "DIRK4",         // DIRK4
};



/* Multistep data type */
typedef struct
{
  int nStep; // number of steps in this time scheme
             // more generally, number of state vectors used
  real alpha[xf_MAXMULTISTEP]; // coefficients on U
  real beta; // coefficient on residual R
}
xf_MultiStepData;


/* Multistep coefficients, relevant for multi-step implicit methods
   only:

   (delta t) * dU/dt = c(0)*U^{n+1} + c(1)*U^n + ...

   where c(i) = MultiStepData(i)

*/
static xf_MultiStepData MultiStepData[xfe_TimeSchemeLast] = {
  {0, {0.0}, 0.0},            		 // Steady
  {1, {1.0, -1.0}, 0.0},      		 // BDF1
  {2, {1.5, -2.0, 0.5}, 0.0}, 	         // BDF2
  {3, {11./6., -3.0, 1.5, -1./3.}, 0.0}, // BDF3
  {1, {2.0, -2.0}, 1.0},      		 // Trapezoidal
  {1, {1.0, -1.0}, 0.0},      		 // FE
  {1, {1.0, -1.0}, 0.0},      		 // RK4
  {1, {1.0, -1.0}, 0.0},      		 // RK2
  {1, {1.0, -1.0}, 0.0},      		 // DG1
  {2, {1.0, -1.0}, 0.0},       		 // DG2
  {1, {1.0, -1.0}, 0.0},      		 // ESDIRK3
  {1, {1.0, -1.0}, 0.0},                 // ESDIRK4
  {1, {1.0, -1.0}, 0.0},                 // ESDIRK5
  {1, {1.0, -1.0}, 0.0},                 // DIRK3
  {1, {1.0, -1.0}, 0.0}                  // DIRK4
};
  
/* Nonlinear Solver type */
enum xfe_NonlinearSolverType{
  xfe_NonlinearSolverNone,
  xfe_NonlinearSolverNewton,
  xfe_NonlinearSolverpMultigrid,
  xfe_NonlinearSolverLast
};

/* corresponding names */
static char *xfe_NonlinearSolverName[xfe_NonlinearSolverLast] = {
  "None",
  "Newton",
  "pMultigrid"
};

/* CFL Evolution type */
enum xfe_CFLEvolutionType{
  xfe_CFL_Exp,
  xfe_CFL_SERA,
  xfe_CFL_SERB,
  xfe_CFL_SERAB,
  xfe_CFL_RDM,
  xfe_CFL_mRDM,
  xfe_CFL_Last
};

/* corresponding names */
static char *xfe_CFLEvolutionName[xfe_CFL_Last] = {
  "Exp",
  "SER-A",
  "SER-B",
  "SER-AB",
  "RDM",
  "mRDM"
};

/* Diffusion discretization type */
enum xfe_DiffDiscType{
  xfe_DiffDiscIP,   // Interior penalty
  xfe_DiffDiscBR2,  // second form of Bassi & Rebay
  xfe_DiffDiscLast
};

/* corresponding names */
static char *xfe_DiffDiscName[xfe_DiffDiscLast] = {
  "IP",
  "BR2"
};

/* Multigrid cycle type */
enum xfe_MultigridCycleType{
  xfe_MultigridCycleVCycle,
  xfe_MultigridCycleFMG,
  xfe_MultigridCycleLast
};

/* corresponding names */
static char *xfe_MultigridCycleName[xfe_MultigridCycleLast] = {
  "VCycle",
  "FMG"
};


typedef struct
{
  int  egrp;
  int  elem;
  real value;
}
xf_ElemReal;




/* Numerically-based stabilization types */
enum xfe_StabilizationType{
  xfe_StabilizationNone,
  xfe_StabilizationResolution, // element resolution-based
  xfe_StabilizationJump,       // face jump-based
  xfe_StabilizationLast
};

/* corresponding names */
static char *xfe_StabilizationName[xfe_StabilizationLast] = {
  "None",
  "Resolution",
  "Jump"
};



/* Stabilization switch type */
enum xfe_StabSwitchType{
  xfe_StabSwitchNone,
  xfe_StabSwitchConst,        // Constant switch value over entire domain
  xfe_StabSwitchLinear,       // Benign linear function of jumps
  xfe_StabSwitchSquare,       // More aggressive square function
  xfe_StabSwitchLog,          // Highly-nonlinear log value
  xfe_StabSwitchLast,
};

/* corresponding names */
static char *xfe_StabSwitchName[xfe_StabSwitchLast] = {
  "None",
  "Const",
  "Linear",
  "Square",
  "Log"
};


/* Space-time anisotropy measure type */
enum xfe_SpaceTimeAnisoType{
  xfe_SpaceTimeAnisoJump,   // solution jumps (heuristic)
  xfe_SpaceTimeAnisoProj,  // order projection
  xfe_SpaceTimeAnisoLast
};

/* corresponding names */
static char *xfe_SpaceTimeAnisoName[xfe_SpaceTimeAnisoLast] = {
  "Jump",
  "Proj"
};

// Structure for storing stabilization terms
typedef struct
{
  enum xfe_StabilizationType StabType; // type of stabilization (see above)

  // for element-based regularity stabilization
  xf_Vector *Reg;            // element regularity estimate; also in All->DataSet
  xf_Vector *Reg_U;          // linearization of elem regularity estimate (also in All)
  
  // for face-based jump stabilization
  xf_Vector *Jump;            // face stabilization terms, also in All
  xf_Vector *Jump_UL;         // linearization of JumpStab w.r.t UL (also in All)
  xf_Vector *Jump_UR;         // linearization of JumpStab w.r.t UR (also in All)


  // For passing into functions; interior faces
  real *StabViscL;
  real *StabViscR;
  real *StabPhiL;
  real *StabPhiR;
  real *StabViscL_UL;
  real *StabViscL_UR;
  real *StabViscR_UL;
  real *StabViscR_UR;
  real *ResMetricL;
  real *ResMetricR;

  // boundary faces
  real *StabVisc;
  real *StabPhi;
  real *StabVisc_U;
  real *ResMetric;

  enum xfe_Bool SkipDiffStab;

}
xf_StabData;


/* SolverData structure (not read or written) */
typedef struct
{
  int  iIter;            // current iteration number
  real CFL;              // current CFL number
  real CFLPrev;             // initial CFL number
  real CFLSafe;          // last safe CFL number
  real CFLDecreaseFactor;// CFL decrease factor
  real CFLIncreaseFactor;// CFL decrease factor
  real CFLMax;           // maximum CFL
  real CFLMin;           // minimum CFL
  real MaxCFLAchieved;   // maximum CFL achieved
  real UpdateFrac;       // UpdateFraction
  real ResNorm;          // Primal residual norm
  real *AdjResNorm;      // Adjoint residual norms for adjoint runs
  real *Output;          // evaluated scalar outputs corresponding to EqnSet->Outputs
  real normdU_prev;      // last calculated |du|
  real LinResTol;        // relative tolerance on linear residual (for linear solve)
  real RnormPrev;        // previous L-2 norm of the primal residual 
  real Rnorm;            // current L-2 norm of the primal residual
  real ResPenaltyPrev;   // previous value of the primal residual taxation (1+P)
  real ResPenalty;       // current value of the primal residual taxation (1+P)
  real gfnormPrev;       // previous L-2 norm of the gradient of f (look at note for ROT) 
  real gfnorm;           // current L-2 norm of the gradient of f (look at note for ROT)
  real muPrev;           // previous penalty factor
  real mu;               // current penalty factor
  
  enum xfe_Bool CRequired; // True if elem-to-elem connectivity needs to be computed
  enum xfe_Bool SortLines; // True if lines need to be sorted
  xf_Vector *C;            // elem-to-elem connectivity vector; also in All->DataSet

  enum xfe_CFLEvolutionType CFLEvolution; //type o CFL evolution
  enum xfe_Bool PenalizeResidual; // if true, Residual gets taxed by ResPenalty
  
  enum xfe_Bool ReusedU; /*if True, the previous dU is used as initial guess 
                          in the linear solve.*/
  
  enum xfe_Bool StabRequired; // True if stabilization is required
  xf_StabData StabData;       // data for various forms of stabilization

  xf_Vector *dt;         // artificial time-step. Also in All->DataSet
  real c;                // constant in front of M*U product for true time stepping

  real TimeStepDecreaseFactor; // time step for robustness control in unsteady simulations
  
  xf_Vector *Pvec;       // pointer to vector of penalty values
  
  enum xfe_Bool SkipParallelExchange; // True for domain decomposition solver

  int ResidualOrderIncrement;  // residual is evaluated at state+increment order

  //particular for Yu's explicit solver
  real currenttime;
  real currenttimestep;
}
xf_SolverData;



/* Time history data (not read or written) */
typedef struct
{
  int nTime;            // number of time points or slabs
  real *Time;           // time value at each point [nTime]
  real *TimeStep;       // time step size associated with slab
  enum xfe_Bool ConstTimeStep; // True if time step is constant
  enum xfe_TimeSchemeType *TimeScheme; // time scheme at each point [nTime], 0 N/A
 
  /* In order to collect time histories of outputs, set nOutput > 0
     and OutputNames to the desired outputs */
  int nOutput;          // number of outputs with a time history
  char **OutputNames;   // names of outputs with a time history [nOutput]
  real **OutputValues;  // calculated output values [nOutput x nTime], 0 N/A

  real *TimeWeights;    // w(i) in Jtot = w(i)*J(U(i)).  NULL means
			// all w(i) are zero except w(T). [nTime], 0 N/A

  /* The following data provide extra weight flexibility for
     SumOutputs, in that the weights for a weighted-sum output can
     vary from time-step to time step.  Main reason for this is to
     support Hessian-based IC model reduction.
   */
  int nSumWeight;          // number of weights per time-point
  real **SumOutputWeights; // weights used for sum outputs in unsteady
			   // adjoint solve [nTime x nweight], 0 N/A
}
xf_TimeHistData;


/* Data for computing diffusion jump residuals */
typedef struct
{
  enum xfe_Bool Need_Grad;
  enum xfe_Bool ConstAL;
  enum xfe_Bool ConstAR;

  real alpha;
  real duL_uL, duL_uR, duR_uL, duR_uR;
  
  real *dunL    , *dunR;
  real *ALguL   , *ARguR;
  real *ALdunL  , *ARdunR;
  real *AL      , *AR;
  real *A_uLguL , *A_uRguR;
  real *A_uLdunL, *A_uRdunR;
  real *N;
  real *Qn;
  real *Qn_dp;
  real *Qn_u;
  real *Qn_gu;
  real *AwStabL , *AwStabR;
  real *T;
  real *Aw;
  real *DnL, *DnR;
  real *SLL, *SLR, *SRL, *SRR;
  real *PL, *PR;

  //July 2015; added by Yu
  real *RL_Fholder, *RR_Fholder;


}
xf_DiffJumpData;


/* Data for computing diffusion BC residuals */
typedef struct
{
  enum xfe_Bool Need_Grad;
  enum xfe_Bool ConstAB;

  real alpha;
  
  real *dun;
  real *du_uI;
  real *du_gb;
  real *ABgu;
  real *ABdun;
  real *AB;
  real *A_uBgu;
  real *A_uBdun;
  real *AB_uIdun;
  real *N;
  real *Qn;
  real *Qn_uI;
  real *Qn_guI;
  real *Qn_gb;
  real *Qn_ggb;
  real *AwStab;
  real *guB;
  real *guB_guI;
  real *T;
  real *Aw;
  real *Dn;
  real *S;
  real *P;
  real *Vw;
}
xf_DiffBCData;

/* Solution update methods */
enum xfe_StateUpdateMethodType{
  xfe_StateUpdateLineSearch,
  xfe_StateUpdateMaxPrimChange,
  xfe_StateUpdateMethodLast
};

/* corresponding names */
static char *xfe_StateUpdateMethodName[xfe_StateUpdateMethodLast] = {
  "LineSearch",
  "MaxPrimitiveChange"
};


#endif // end ifndef _xf_SolverStruct_h
