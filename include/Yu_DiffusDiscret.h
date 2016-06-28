/* header for Yu_DiffusDiscret.c */
// some notes on this subroutine
// Implemented based on BR2; fully explicit;
// Support variable thermodynamics
// Author: ylv@stanford.edu
// Update till May, 2013
//
/******************************************************************/
//   FUNCTION Definition: LYDG_CreateDiffJumpData
int
LYDG_CreateDiffJumpData(xf_DiffJumpData **pDData);

/******************************************************************/
//   FUNCTION Definition: LYDG_ReAllocDiffJumpData
int
LYDG_ReAllocDiffJumpData(int sr, int nq, int dim, int nn,
                         enum xfe_Bool Need_Grad, xf_DiffJumpData *DData);

/******************************************************************/
//   FUNCTION Definition: LYDG_DestroyDiffJumpData
void
LYDG_DestroyDiffJumpData(xf_DiffJumpData *DData);

/******************************************************************/
//   FUNCTION Definition: LYDG_InteriorViscTerm
int
LYDG_InteriorViscTerm(const Yu_Model *Model, int nDiff, int nq, real *u,
                      real *gu, real *w, real *Aw, real *A, const real gamma);

/******************************************************************/
//   FUNCTION Definition: LYDG_FaceViscTerm
int
LYDG_FaceViscTerm(xf_All *All, Yu_Model *Model, const int iiface, const xf_BasisData *PhiDataL,
                  const xf_BasisData *PhiDataR, real ElemVolL, real ElemVolR,
                  enum xfe_BasisType BasisL, enum xfe_BasisType BasisR, int OrderL,
                  int OrderR, int nq, const real *wn, real *uL, real *uR, real *guL,
                  real *guR, real *RL, real *RR, const real gammaL, const real gammaR);

/******************************************************************/
//   FUNCTION Definition: LYDG_CreateDiffBCData
int
xf_CreateDiffBCData(xf_DiffBCData **pDData);

/******************************************************************/
//   FUNCTION Definition: LYDG_ReAllocDiffBCData
int
LYDG_ReAllocDiffBCData(int sr, int nq, int dim, int nn,
                       enum xfe_Bool Need_Grad, xf_DiffBCData *DData);

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDiffBCData
void
LYDG_DestroyDiffBCData(xf_DiffBCData *DData);

/******************************************************************/
//   FUNCTION Definition:  LYDG_SetVisFluxBC
static int
LYDG_SetVisFluxBC(const int *SetIndex, int nq, int sr, int dim, const real *guB,
                  const real *Qn, real *Aw, enum xfe_Bool ZeroFlag);

/******************************************************************/
//   FUNCTION Definition:  LYDG_BoundaryViscTerm
int
LYDG_BoundaryViscTerm(xf_All *All, Yu_Model *Model, int ibfgrp, int ibface,
                      const xf_BasisData *PhiData, const xf_BasisData *ResPhiData,
                      real ElemVol, enum xfe_BasisType Basis, int Order, int nq,
                      const real *wn, const real *xglob, real *uI, real *uB,
                      real *guI, real *ER, xf_DiffBCData *DData const real gamma,
                      xf_OutputEvalData *OutputEval);

