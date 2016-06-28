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
  FILE:  xf_EqnSetHook.c

  This file contains "hooks" for equation-set specific functions.  It
  is used to dynamically load functions that reside in equation-set
  specific libraries.

*/

#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_IO.h"
#include "xf_Dynamic.h"


static xf_DLHandle LibHandle = NULL;



/* The following objects are pointers to the respective equation-set
   specific funcitons in the dynamically-loaded library.  The names
   are stripped of the xf_ prefix to prevent a namespace conflict with
   the hook functions defined further down in this file.  See
   xf_EqnSetHook.h for documentation on the individual functions. */

static int (*EqnSetRegister) (xf_EqnSet *EqnSet)
= NULL;

static int (*EqnSetQuadOrderElem) (const xf_EqnSet *EqnSet, int Order, int *QuadOrder)
= NULL;

static int (*EqnSetQuadOrderFace) (const xf_EqnSet *EqnSet, int Order, int *QuadOrder)
= NULL;

static int (*EqnSetICState) (const xf_EqnSet *EqnSet, const xf_IC *IC,
			     const int *IParam, const real *RParam, int nq,
			     const real *xglob, const real *pTime, real *u) 
= NULL;

static int (*EqnSetBCState) (const xf_EqnSet *EqnSet, const xf_BC *BC,
			     const int *IParam, const real *RParam, int nq, 
			     const real *n, const real *xglob, const real *pTime,
			     const real *vg, const real *UI, real *UB, real *UB_UI)
= NULL;

static int (*EqnSetFcnState) (const xf_EqnSet *EqnSet, const char *Fcn, const char *Data,
			      const int *IParam, const real *RParam, int nq,
			      const real *xglob, const real *pTime, real *u) 
= NULL;

static int (*EqnSetConvF) (const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
			   int nConv, const int *IParam, const real *RParam, int nq, 
			   const real *U, real **AuxU, const real *xglob, 
			   real *F, real *F_U)
= NULL;

static int (*EqnSetConvFJump) (const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
			       int nConv, const int *IParam, const real *RParam, int nq, 
			       const real *UL, const real *UR, real **AuxUL,
			       real **AuxUR, const real *n, const real *xglob, 
			       const real *vg, real *F, real *F_UL, real *F_UR, 
			       real *C)
= NULL;

static int (*EqnSetConvFBC) (const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, int nConv, 
			     const xf_BC *BC, const int *IParam, const real *RParam, 
			     int nq, const real *UI, const real *UB, const real *UB_UI,
			     real **AuxUB, const real *n, const real *xglob,
			     const real *vg, real *F, real *F_U, real *C)
= NULL;

static int (*EqnSetDiffA) (const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
			   const int *IParam, const real *RParam, const real *ResMetric, 
			   int nq, const real *U, const real *gU, const real *W, 
			   real *AW, real *A, real *A_UW, enum xfe_Bool *ConstA, real *mu)
= NULL;

static int (*EqnSetDiffFBC) (const xf_EqnSet *EqnSet, const xf_BC *BC, const int *IParam, 
			     const real *RParam, int nq, const real *n, const real *xglob, 
			     const real *gUI, real *gUB, real *gUB_gUI, 
			     enum xfe_Bool *AFlag, int *SetIndex)
= NULL;

static int (*EqnSetSourceS) (const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
			     int nSource, const int *IParam, const real *RParam, int nq, 
			     const real *U, const real *gU, real **qAuxU, const real *xglob,
			     const real *pTime, real *S, real *S_U, real *S_gU, 
			     enum xfe_Bool *Nonzero_gu)
= NULL;

static int (*EqnSetUpdateFraction) (const xf_EqnSet *EqnSet, const int *IParam, const real *RParam,
				 int nq, const real *u, const real *du, real *frac)
= NULL;

static int (*EqnSetMaxCharSpeed) (const xf_EqnSet *EqnSet, int nq, const real *U, 
			       real **AuxU, const real *xglob, const int *IParam,
			       const real *RParam, real *vmax, real *dtmax, real *v, real *vmax_U)
= NULL;

static int (*EqnSetScalar) (const xf_EqnSet *EqnSet, const char *Name, const int *IParam, 
			    const real *RParam, int nq, const real *U, const real *gU,
			    real *s, real *s_U, int *nScalar, char ***pScalarNames)
= NULL;

static int (*EqnSetVector) (const xf_EqnSet *EqnSet, const char *Name, const int *IParam, 
			    const real *RParam, int nq, const real *U, real *v, real *v_U)
= NULL;

static int (*EqnSetVariableChange) (const xf_EqnSet *EqnSet, const char *Name, const int *IParm, 
				    const real *RParam, int nq, const real *U, real *V,
				    int *nVSet, char ***pVSetNames)
= NULL;

static int (*EqnSetScaleState) (const xf_EqnSet *EqnSet, const xf_IC *ICOrig, const int *IParm, 
				const real *RParam, int nq, real *U)
= NULL;

static int (*EqnSetBCIsWall) (const xf_BC *BC, enum xfe_Bool *IsWall)
= NULL;

static int (*EqnSetAlterState) (const xf_EqnSet *EqnSet, const char *AlterFcn, const int *IParm, 
				const real *RParam, const real *qxglob, const real *FcnParam, 
				int nq, real *qU)
= NULL;

static int (*EqnSetPenaltyTerm) (const xf_EqnSet *EqnSet, const real *RParam, const int *IParam, 
                                 int nq, real *u, real *pPe, real *Pq, real *Pe_u,
                                 real *Pe_uu)
= NULL;

static int (*EqnSetPerturbParam) (xf_EqnSet *EqnSet, xf_Sensitivity *Sensitivity, 
                                  real *epsilon, enum xfe_Bool CorrectFlag)
= NULL;

static int (*EqnSetOutputDependentBCs) (xf_EqnSet *EqnSet, xf_KeyValue *Outputs,
                                        real *Sensitivity, real *epsilon, 
                                        enum xfe_Bool *Converged)
= NULL;



/* The following functions are used for dynamically opening and
   closing the equation-set library. */

/******************************************************************/
//   FUNCTION Definition: xf_LoadEqnSetLibrary
int
xf_LoadEqnSetLibrary(const char *LibName)
{

  xf_printf("Loading EqnSet Library = %s\n", LibName);

  if (!(LibHandle = xf_DLOpen(LibName))){
    xf_printf("Error during xf_DLOpen: %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetRegister = xf_DLSym(LibHandle, "xf_EqnSetRegister");
  if (!EqnSetRegister){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetQuadOrderElem = xf_DLSym(LibHandle, "xf_EqnSetQuadOrderElem");
  if (!EqnSetQuadOrderElem){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetQuadOrderFace= xf_DLSym(LibHandle, "xf_EqnSetQuadOrderFace");
  if (!EqnSetQuadOrderFace){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetICState = xf_DLSym(LibHandle, "xf_EqnSetICState");
  if (!EqnSetICState){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetBCState = xf_DLSym(LibHandle, "xf_EqnSetBCState");
  if (!EqnSetBCState){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetFcnState = xf_DLSym(LibHandle, "xf_EqnSetFcnState");
  /* OK if not present */

  EqnSetConvF = xf_DLSym(LibHandle, "xf_EqnSetConvF");
  /* OK if not present */

  EqnSetConvFJump = xf_DLSym(LibHandle, "xf_EqnSetConvFJump");
  /* OK if not present */

  EqnSetConvFBC = xf_DLSym(LibHandle, "xf_EqnSetConvFBC");
  /* OK if not present */

  EqnSetDiffA = xf_DLSym(LibHandle, "xf_EqnSetDiffA");
  /* OK if not present */

  EqnSetDiffFBC = xf_DLSym(LibHandle, "xf_EqnSetDiffFBC");
  /* OK if not present */

  EqnSetSourceS = xf_DLSym(LibHandle, "xf_EqnSetSourceS");
  /* OK if not present */

  EqnSetScalar = xf_DLSym(LibHandle, "xf_EqnSetScalar");
  /* OK if not present */

  EqnSetVector = xf_DLSym(LibHandle, "xf_EqnSetVector");
  /* OK if not present */

  EqnSetVariableChange = xf_DLSym(LibHandle, "xf_EqnSetVariableChange");
  /* OK if not present */

  EqnSetScaleState = xf_DLSym(LibHandle, "xf_EqnSetScaleState");
  /* OK if not present */

  EqnSetBCIsWall = xf_DLSym(LibHandle, "xf_EqnSetBCIsWall");
  /* OK if not present */

  EqnSetAlterState = xf_DLSym(LibHandle, "xf_EqnSetAlterState");
  /* OK if not present */

  EqnSetPenaltyTerm = xf_DLSym(LibHandle, "xf_EqnSetPenaltyTerm");
  if (!EqnSetPenaltyTerm){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetUpdateFraction = xf_DLSym(LibHandle, "xf_EqnSetUpdateFraction");
  if (!EqnSetUpdateFraction){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetMaxCharSpeed= xf_DLSym(LibHandle, "xf_EqnSetMaxCharSpeed");
  if (!EqnSetMaxCharSpeed){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  EqnSetPerturbParam = xf_DLSym(LibHandle, "xf_EqnSetPerturbParam");
  if (!EqnSetPerturbParam){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }
  
  EqnSetOutputDependentBCs = xf_DLSym(LibHandle, "xf_EqnSetOutputDependentBCs");
  if (!EqnSetOutputDependentBCs){
    xf_printf("object not found.  xf_DLSym error = %s\n", xf_DLError());
    return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CloseEqnSetLibrary
void
xf_CloseEqnSetLibrary()
{
  xf_DLClose(LibHandle);
}


/* The following functions are the "hook"/entry-point functions.  When
   called, they just pass their arguments to the dynamically-loaded
   functions. */

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetRegister
int 
xf_EqnSetRegister(xf_EqnSet *EqnSet){
  if (!EqnSetRegister) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetRegister(EqnSet);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetQuadOrderElem
int 
xf_EqnSetQuadOrderElem(const xf_EqnSet *EqnSet, int Order, int *QuadOrder){
  if (!EqnSetQuadOrderElem) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetQuadOrderElem(EqnSet, Order, QuadOrder);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetQuadOrderFace
int 
xf_EqnSetQuadOrderFace(const xf_EqnSet *EqnSet, int Order, int *QuadOrder){
  if (!EqnSetQuadOrderFace) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetQuadOrderFace(EqnSet, Order, QuadOrder);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetICState
int 
xf_EqnSetICState(const xf_EqnSet *EqnSet, const xf_IC *IC, 
		 const int *IParam, const real *RParam, int nq, 
		 const real *xglob, const real *pTime, real *u){
  if (!EqnSetICState) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetICState(EqnSet, IC, IParam, RParam, nq, xglob, pTime, u);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetBCState
int 
xf_EqnSetBCState(const xf_EqnSet *EqnSet, const xf_BC *BC,
		 const int *IParam, const real *RParam, int nq, 
		 const real *n, const real *xglob, const real *pTime, 
		 const real *vg, const real *UI, real *UB, real *UB_UI){
  if (!EqnSetBCState) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetBCState(EqnSet, BC, IParam, RParam, nq, n, xglob, 
		       pTime, vg, UI, UB, UB_UI);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetFcnState
int 
xf_EqnSetFcnState(const xf_EqnSet *EqnSet, const char *Fcn, 
		  const char *Data, const int *IParam, 
		  const real *RParam, int nq, const real *xglob, 
		  const real *pTime, real *u){
  if (!EqnSetFcnState) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetFcnState(EqnSet, Fcn, Data, IParam, RParam, nq, xglob, pTime, u);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetConvF
int 
xf_EqnSetConvF(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
	       int nConv, const int *IParam, const real *RParam, int nq, 
	       const real *U, real **AuxU, const real *xglob, 
	       real *F, real *F_U)
{
  if (!EqnSetConvF) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetConvF(EqnSet, ResTerm, nConv, IParam, RParam, nq, U, 
		     AuxU, xglob, F, F_U);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetConvFJump
int 
xf_EqnSetConvFJump(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
		   int nConv, const int *IParam, const real *RParam, int nq, 
		   const real *UL, const real *UR, real **AuxUL,
		   real **AuxUR, const real *n, const real *xglob, 
		   const real *vg, real *F, real *F_UL, real *F_UR, real *C)
{
  if (!EqnSetConvFJump) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetConvFJump(EqnSet, ResTerm, nConv, IParam, RParam, nq, 
			 UL, UR, AuxUL, AuxUR, n, xglob, vg, F, F_UL, F_UR, C);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetConvFBC
int 
xf_EqnSetConvFBC(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, int nConv,
		 const xf_BC *BC, const int *IParam, const real *RParam, 
		 int nq, const real *UI, const real *UB, const real *UB_UI, 
		 real **AuxUB, const real *n, const real *xglob, 
		 const real *vg, real *F, real *F_U, real *C)
{
  if (!EqnSetConvFBC) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetConvFBC(EqnSet, ResTerm, nConv, BC, IParam, RParam, nq, 
		       UI, UB, UB_UI, AuxUB, n, xglob, vg, F, F_U, C);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetDiffA
int 
xf_EqnSetDiffA(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
	       const int *IParam, const real *RParam, const real *ResMetric, 
	       int nq, const real *U, const real *gU, const real *W, 
	       real *AW, real *A, real *A_UW, enum xfe_Bool *ConstA, real *mu)
{
  if (!EqnSetDiffA) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetDiffA(EqnSet, ResTerm, IParam, RParam, ResMetric, nq, 
		     U, gU, W, AW, A, A_UW, ConstA, mu);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetDiffFBC
int 
xf_EqnSetDiffFBC(const xf_EqnSet *EqnSet, const xf_BC *BC, const int *IParam, 
		 const real *RParam, int nq, const real *wn, const real *xglob, 
		 const real *gUI, real *gUB, real *gUB_gUI, 
		 enum xfe_Bool *AFlag, int *SetIndex)
{
  if (!EqnSetDiffFBC) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetDiffFBC(EqnSet, BC, IParam, RParam, nq, wn, xglob, 
		       gUI, gUB, gUB_gUI, AFlag, SetIndex);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetSourceS
int 
xf_EqnSetSourceS(const xf_EqnSet *EqnSet, const xf_ResTerm *ResTerm, 
		 int nSource, const int *IParam, const real *RParam, int nq, 
		 const real *U, const real *gU, real **qAuxU, const real *xglob, 
		 const real *pTime, real *S, real *S_U, real *S_gU, 
		 enum xfe_Bool *Nonzero_gu)
{
  if (!EqnSetSourceS) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetSourceS(EqnSet, ResTerm, nSource, IParam, RParam, nq, 
		       U, gU, qAuxU, xglob, pTime, S, S_U, S_gU, Nonzero_gu);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetScalar
int 
xf_EqnSetScalar(const xf_EqnSet *EqnSet, const char *Name, const int *IParam, 
		const real *RParam, int nq, const real *U, const real *gU, 
		real *s, real *s_U, int *nScalar, char ***pScalarNames)
{
  if (!EqnSetScalar) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetScalar(EqnSet, Name, IParam, RParam, nq, U, gU, s, s_U, 
		      nScalar, pScalarNames);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetVector
int 
xf_EqnSetVector(const xf_EqnSet *EqnSet, const char *Name, const int *IParam, 
		const real *RParam, int nq, const real *U, real *v, real *v_U)
{
  if (!EqnSetVector) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetVector(EqnSet, Name, IParam, RParam, nq, U, v, v_U);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetVariableChange
int 
xf_EqnSetVariableChange(const xf_EqnSet *EqnSet, const char *Name, const int *IParam, 
			const real *RParam, int nq, const real *U, real *V,
			int *nVSet, char ***pVSetNames)
{
  if (!EqnSetVariableChange) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetVariableChange(EqnSet, Name, IParam, RParam, nq, U, V,
			      nVSet, pVSetNames);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetScaleState
int 
xf_EqnSetScaleState(const xf_EqnSet *EqnSet, const xf_IC *ICOrig, const int *IParam, 
		    const real *RParam, int nq, real *U)
{
  if (!EqnSetScaleState) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetScaleState(EqnSet, ICOrig, IParam, RParam, nq, U);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetAlterState
int 
xf_EqnSetAlterState(const xf_EqnSet *EqnSet, const char *AlterFcn, const int *IParam, 
		    const real *RParam, const real *xglob, const real *FcnParam, 
		    int nq, real *U)
{
  if (!EqnSetAlterState) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetAlterState(EqnSet, AlterFcn, IParam, RParam, xglob, FcnParam, nq, U);
}


/******************************************************************/
//   FUNCTION Definition: xf_EqnSetUpdateFraction
int 
xf_EqnSetUpdateFraction(const xf_EqnSet *EqnSet, const int *IParam, const real *RParam,
			const int nq, const real *u, const real *du, real *frac){
  if (!EqnSetUpdateFraction) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetUpdateFraction(EqnSet, IParam, RParam, nq, u, du, frac);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetMaxCharSpeed
int 
xf_EqnSetMaxCharSpeed(const xf_EqnSet *EqnSet, int nq, const real *U, 
		   real **AuxU, const real *xglob, const int *IParam,
		   const real *RParam, real *vmax, real *dtmax,  real *v,  real *v_U){
  if (!EqnSetMaxCharSpeed) return xf_Error(xf_DYNAMIC_LIBRARY_ERROR);
  return EqnSetMaxCharSpeed(EqnSet, nq, U, AuxU, xglob, IParam, RParam, vmax, dtmax, v, v_U);
}



/******************************************************************/
//   FUNCTION Definition: xf_EqnSetBCIsWall
int 
xf_EqnSetBCIsWall(const xf_BC *BC, enum xfe_Bool *IsWall)
{
  if (!EqnSetBCIsWall) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetBCIsWall(BC, IsWall);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetPenaltyTerm
int 
xf_EqnSetPenaltyTerm(const xf_EqnSet *EqnSet, const real *RParam, 
                     const int *IParam, int nq, real *u, real *pPe, 
                     real *Pq, real *Pe_u, real *Pe_uu)
{
  if (!EqnSetPenaltyTerm) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetPenaltyTerm(EqnSet, RParam, IParam, nq, u, pPe, Pq, Pe_u, Pe_uu);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetPerturbParam
int 
xf_EqnSetPerturbParam(xf_EqnSet *EqnSet, xf_Sensitivity *Sensitivity, 
                      real *epsilon, enum xfe_Bool CorrectFlag)
{
  if (!EqnSetPerturbParam) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetPerturbParam(EqnSet, Sensitivity, epsilon, CorrectFlag);
}

/******************************************************************/
//   FUNCTION Definition: xf_EqnSetOutputDependentBCs
int 
xf_EqnSetOutputDependentBCs(xf_EqnSet *EqnSet, xf_KeyValue *Outputs,
                            real *Sensitivity, real *epsilon, 
                            enum xfe_Bool *Converged)
{
  if (!EqnSetOutputDependentBCs) return xf_DYNAMIC_LIBRARY_ERROR;
  return EqnSetOutputDependentBCs(EqnSet, Outputs, Sensitivity, 
                                  epsilon, Converged);
}

