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
 FILE:  xf_LinearSolver.c
 
 This file contains linear solver functions, where the linear system
 of interest involves the Jacobian matrix, R_U.
 
 */

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_LinearSolverStruct.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Solver.h"
#include "xf_LeanSolver.h"


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_n
int
xf_Jacobian_n(xf_JacobianMatrix *R_U, int egrp, int elem)
{
  // pulls off n (# unknowns) for egrp,elem, allowing for variable orders
  int n;
  
  if (egrp < 0)
    n = 0;
  else if ((elem < 0) || (R_U->vnvec == NULL))
    n = R_U->nvec[egrp];
  else
    n = R_U->vnvec[egrp][elem];

  return n;
}


// ILU linear solver (temporarily included as a .c file)
#include "xf_LinearSolverILU.c"


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_Precondition_Block
static int
xf_Jacobian_Precondition_Block(xf_All *All, xf_JacobianMatrix *R_U)
{  
  int ierr, egrp, elem;
  xf_Mesh *Mesh;
  Mesh = All->Mesh;
  
  // loop over elements, PLU factor diagonal blocks
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      ierr = xf_Error(xf_ComputeBlockPLU(R_U->Value[egrp][elem][0], 
					 xf_Jacobian_n(R_U,egrp,elem),
                                         R_U->StateRank, R_U->P[egrp][elem]));
      if (ierr != xf_OK) return ierr;
    } // elem
  } // egrp
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_Precondition_Line
static int
xf_Jacobian_Precondition_Line(xf_All *All, xf_JacobianMatrix *R_U)
{  
  /*
   PURPOSE:
   
   Computes a quasi-LU factorization of the block-tridiagonal systems
   in R_U associated with the lines. The factorization is "quasi-LU" in
   the sense that the lower diagonal blocks remain unchanged (whereas
   in real LU these would be the blocks of L).  Thus, only the diagonal
   blocks are altered.  For example, for a block 3x3 system:
   
   [ A0  B0  0  ]      [ LU(A0)    B0      0     ]
   [ C0  A1  B1 ]  ->  [   C0    LU(A1')   B1    ]
   [ 0   C1  A2 ]      [   0       C1    LU(A2') ]
   
   where:   A1' = A1 - C0 * inv(A0 ) * B0                   
   A2' = A2 - C1 * inv(A1') * B1                   
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix
   
   OUTPUTS:
   
   R_U : Preconditioned Jacobian matrix
   
   RETURN:  Error code
   
   */
  int ierr, iLine, ie;
  int egrp, elem, face, egrpN, elemN, faceN;
  int *P, n0, n1, sr, Tsize = 0;
  real *R00, *R01, *R10, *R11, *T = NULL;
  xf_LineSet *LineSet;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  sr   = R_U->StateRank;
  
  if ((LineSet = R_U->LineSet) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // loop over lines
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    
    // loop over elements in iLine 
    for (ie=0; ie<LineSet->Line[iLine].nelem; ie++){
      
      // pull off info of elem ie on line
      egrp = LineSet->Line[iLine].egrp[ie];
      elem = LineSet->Line[iLine].elem[ie];
      face = LineSet->Line[iLine].face[ie];
      
      R00 = R_U->Value[egrp][elem][0];    // R_U(elem,elem)
      n0  = xf_Jacobian_n(R_U,egrp,elem); // r/StateRank for egrp,elem
      P   = R_U->P[egrp][elem];           // permutation vector for elem
      
      // PLU factor R00
      ierr = xf_Error(xf_ComputeBlockPLU(R00, n0, sr, P));
      if (ierr != xf_OK) return ierr;
      
      if (ie != (LineSet->Line[iLine].nelem-1)){ // if not last elem
        
        // pull off info of next elem
        egrpN = R_U->egrpN[egrp][elem][face];
        elemN = R_U->elemN[egrp][elem][face];
        faceN = R_U->faceN[egrp][elem][face]; // elem is on other side
        
        R01 = R_U->Value[egrp ][elem ][1+face ]; // R_U(elem ,elemN)
        R10 = R_U->Value[egrpN][elemN][1+faceN]; // R_U(elemN,elem )
        R11 = R_U->Value[egrpN][elemN][0      ]; // R_U(elemN,elemN)
	n1  = xf_Jacobian_n(R_U,egrpN,elemN);    // r/StateRank for egrpN,elemN
        
        // reallocate T if necessary; T is n0 x n1 blocks
        if (n0*n1*sr*sr > Tsize){
          Tsize = n0*n1*sr*sr;
          ierr = xf_Error(xf_ReAlloc( (void **) &T, Tsize, sizeof(real)));
          if (ierr != xf_OK) return ierr;
        }
        
        // T = R00^{-1} * R01
        ierr = xf_Error(xf_SolveBlockPLU_Matrix(R00, n0, sr, P, R01, n1, xfe_Set, T));
        if (ierr != xf_OK) return ierr;
        
        // R11 -= R10 * T
        ierr = xf_Error(xf_BlockMxBlockM(R10, n1, sr, n0, T, n1, xfe_Sub, R11));
        if (ierr != xf_OK) return ierr;
      }
      
    } // ie
    
  } // iLine
  
  // release memory
  xf_Release( (void *) T);
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_Precondition
static int
xf_Jacobian_Precondition(xf_All *All, xf_JacobianMatrix *R_U,
                         enum xfe_PreconditionerType Preconditioner,
                         enum xfe_Bool CoarseCorrectionFlag)
{
  /*
   PURPOSE:
   
   Preconditions R_U (in-place) according to Preconditioner.  R_U must
   either be not preconditioned or already preconditioned with
   Preconditioner.  Otherwise an error results.  Info on the
   preconditioning status is stored in R_U, so that this function can
   be called multiple times without worry of recomputing the same
   preconditioner.
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix
   Preconditioner : desired preconditioner
   
   OUTPUTS:
   
   R_U : Preconditioned Jacobian matrix
   
   RETURNS:
   
   Error Code
   
   */  
  
  int ierr;
  int egrp, negrp;
  int CoarseOrder;
  enum xfe_Bool LeanFlag;
  xf_JacobianMatrix *R_Uc;
  
  if (R_U->Preconditioner == Preconditioner) // already preconditioned
    return xf_OK;
  
  if (R_U->Preconditioner != xfe_PreconditionerNone)
    return xf_Error(xf_INPUT_ERROR); // some other preconditioner already applied
  
  // Is the preconditioner memory-lean?
  ierr = xf_Error(xf_PreconditionerLeanCheck(Preconditioner, &LeanFlag));
  if (ierr != xf_OK) return ierr;
  
  if (CoarseCorrectionFlag == xfe_True){
    
    if (LeanFlag == xfe_True) return xf_Error(xf_NOT_SUPPORTED);
    
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "CoarseCorrectionOrder", 
                                      &CoarseOrder)); 
    if (ierr != xf_OK) return ierr;
    
    if (R_U->R_Uc == NULL){
      ierr = xf_Error(xf_AllocateJacobianMatrixCoarse(All->Mesh, CoarseOrder, R_U));
      if (ierr != xf_OK) return ierr;
    }
    
    R_Uc = (xf_JacobianMatrix *) R_U->R_Uc;
    
    ierr = xf_Error(xf_ProjectJacobian(All, R_U, R_Uc));
    if (ierr != xf_OK) return ierr;
    
    R_Uc->LineSet = R_U->LineSet;
    R_Uc->Preconditioner = xfe_PreconditionerNone; // just to make sure
    
    ierr = xf_Error(xf_Jacobian_Precondition(All, R_Uc, Preconditioner, xfe_False));
    if (ierr != xf_OK) return ierr;
  }
  
  
  if (R_U->P == NULL){
    negrp = All->Mesh->nElemGroup;
    ierr = xf_Error(xf_Alloc((void **) &R_U->P, negrp, sizeof(int **)));
    if (ierr != xf_OK) return ierr;
    for (egrp=0; egrp<negrp; egrp++){
      if (R_U->vnvec == NULL)
	ierr = xf_Error(xf_Alloc2((void ***) &R_U->P[egrp], All->Mesh->ElemGroup[egrp].nElem,
				  R_U->nvec[egrp]*R_U->StateRank, sizeof(int)));
      else
	ierr = xf_Error(xf_VAlloc2((void ***) &R_U->P[egrp], All->Mesh->ElemGroup[egrp].nElem,
				   R_U->vnvec[egrp], R_U->StateRank*sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  switch (Preconditioner){
  case xfe_PreconditionerBlockJacobi:
    ierr = xf_Error(xf_Jacobian_Precondition_Block(All, R_U));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerBlockJacobiLean:
    // nothing to do for lean preconditioning
    break;
  case xfe_PreconditionerLineJacobi:
  case xfe_PreconditionerLineGS:
    ierr = xf_Error(xf_Jacobian_Precondition_Line(All, R_U));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerILU0:
    ierr = xf_Error(xf_Jacobian_Precondition_ILU0(All, R_U));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  
  R_U->Preconditioner = Preconditioner;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_BlockJacobi
static int
xf_Jacobian_SolveM_BlockJacobi(xf_All *All, xf_JacobianMatrix *R_U,
                               xf_Vector *X, enum xfe_AddType AddFlag,
                               enum xfe_Bool TransposeFlag)
{  
  int ierr, egrp, elem, nn;
  xf_Mesh *Mesh;
  Mesh = All->Mesh;
  
  // loop over elements, PLU solve diagonal blocks, apply to X
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      nn = xf_Jacobian_n(R_U,egrp,elem);
      if (!TransposeFlag){
        ierr = xf_Error(xf_SolveBlockPLU(R_U->Value[egrp][elem][0], nn,
                                         R_U->StateRank, R_U->P[egrp][elem],
                                         X->GenArray[egrp].rValue[elem], AddFlag,
                                         X->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr;
      }
      else{
        ierr = xf_Error(xf_SolveBlockPLUT(R_U->Value[egrp][elem][0], nn,
                                          R_U->StateRank, R_U->P[egrp][elem],
                                          X->GenArray[egrp].rValue[elem], AddFlag,
                                          X->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr;
      }
      
    } // elem
  } // egrp
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_OneLine
static int
xf_Jacobian_SolveM_OneLine(xf_Line *Line, xf_JacobianMatrix *R_U,
                           xf_Vector *X, enum xfe_AddType AddFlag,
                           enum xfe_Bool TransposeFlag, real *XT)
{  
  // X = M^{-1} * X
  int ierr, k, egrp, elem, face, r, sr;
  int egrpN, elemN, faceN, n0, n1;
  int iLine, ie, *P;
  real *R00, *R01, *R10, *X0, *X1;
  
  // pull of StateRank
  sr = R_U->StateRank;
  
  // only support Set for now
  if (AddFlag != xfe_Set) return xf_Error(xf_NOT_SUPPORTED);
  
  /*-------------------------------------------------*/
  /* (I) Set:   X' = inv(L) X                        */
  /*                                                 */
  /* inv(L) = [    I         0       0]   (3x3 e.g.) */
  /*          [-C0*A0'^-1    I       0]              */
  /*          [   0      -C1*A1'^-1  I]              */
  /*                                                 */
  /* Note that in the quasi-LU factorization of R_U, */
  /* the C blocks are stored on the lower diagonal,  */
  /* and the A' blocks are stored, PLU factored, on  */
  /* the diagonal blocks.                            */
  /*-------------------------------------------------*/
  
  
  // loop over elements in a line
  for (ie=0; ie<(Line->nelem-1); ie++){
    
    // pull off info of elem ie on line
    egrp = Line->egrp[ie];
    elem = Line->elem[ie];
    face = Line->face[ie];
    
    R00 = R_U->Value[egrp][elem][0];    // R_U(elem,elem); PLU factored
    n0  = xf_Jacobian_n(R_U,egrp,elem); // rank/sr for egrp,elem
    P   = R_U->P[egrp][elem];           // permutation vector for elem
    X0  = X->GenArray[egrp].rValue[elem]; // X(elem)
    
    // pull off info of next elem
    egrpN = R_U->egrpN[egrp][elem][face];
    elemN = R_U->elemN[egrp][elem][face];
    faceN = R_U->faceN[egrp][elem][face]; // elem is on other side
    
    R10 = R_U->Value[egrpN][elemN][1+faceN]; // R_U(elemN,elem )
    n1  = xf_Jacobian_n(R_U,egrpN,elemN);    // rank/sr for egrpN,elemN
    X1  = X->GenArray[egrpN].rValue[elemN];  // X(elemN)
    
    if (!TransposeFlag){
      // XT = R00^{-1} * X0
      ierr = xf_Error(xf_SolveBlockPLU(R00, n0, sr, P, X0, xfe_Set, XT));
      if (ierr != xf_OK) return ierr;
      
      // X1 -= R10 * XT
      ierr = xf_Error(xf_BlockMxV(R10, n1, sr, n0, XT, xfe_Sub, X1));
      if (ierr != xf_OK) return ierr;
    }
    else{
      // XT = R00^{-T} * X0
      ierr = xf_Error(xf_SolveBlockPLUT(R00, n0, sr, P, X0, xfe_Set, XT));
      if (ierr != xf_OK) return ierr;
      
      // X1 -= R01^T * XT
      R01 = R_U->Value[egrp][elem][1+face]; // R_U(elem ,elemN)
      ierr = xf_Error(xf_BlockMTxV(R01, n1, sr, n0, XT, xfe_Sub, X1));
      if (ierr != xf_OK) return ierr;
    }
    
  } // ie
  
  
  /*-----------------------------------------------*/
  /* (II) Solve for X via back substitution:       */
  /*                                               */
  /*             U X = X'                          */
  /*                                               */
  /*    U = [ A0'   B0    0 ]        (3x3 e.g.)    */
  /*        [ 0     A1'   B1]                      */
  /*        [ 0     0    A2']                      */
  /*                                               */
  /*-----------------------------------------------*/
  
  // loop over elements in a line
  for (ie=(Line->nelem-1); ie>=0; ie--){
    
    // pull off info of elem ie on line
    egrp = Line->egrp[ie];
    elem = Line->elem[ie];
    face = Line->face[ie];
    
    R00 = R_U->Value[egrp][elem][0];    // R_U(elem,elem); PLU factored
    n0  = xf_Jacobian_n(R_U,egrp,elem); // rank/sr for egrp
    P   = R_U->P[egrp][elem];           // permutation vector for elem
    X0  = X->GenArray[egrp].rValue[elem]; // X(elem)
    
    if (ie < Line->nelem-1){ // if not last element
      // pull off info of next elem
      egrpN = R_U->egrpN[egrp][elem][face];
      elemN = R_U->elemN[egrp][elem][face];
      faceN = R_U->faceN[egrp][elem][face];
      
      R01 = R_U->Value[egrp][elem][1+face];   // R_U(elem,elemN)
      n1  = xf_Jacobian_n(R_U,egrpN,elemN);   // rank/sr for egrpN,elemN
      X1  = X->GenArray[egrpN].rValue[elemN]; // X(elemN)
      
      if (!TransposeFlag){
        // X0 -= R01 * X1
        ierr = xf_Error(xf_BlockMxV(R01, n0, sr, n1, X1, xfe_Sub, X0));
        if (ierr != xf_OK) return ierr;
      }
      else{
        // X0 -= R10^T * X1
        R10 = R_U->Value[egrpN][elemN][1+faceN]; // R_U(elemN,elem )
        ierr = xf_Error(xf_BlockMTxV(R10, n0, sr, n1, X1, xfe_Sub, X0));
        if (ierr != xf_OK) return ierr;
      }
    }
    
    if (!TransposeFlag){
      // X0 @= R00^{-1} * X0,  @ is AddFlag (Set)
      ierr = xf_Error(xf_SolveBlockPLU(R00, n0, sr, P, X0, AddFlag, X0));
      if (ierr != xf_OK) return ierr;
    }
    else{
      // X0 @= R00^{-T} * X0,  @ is AddFlag (Set)
      ierr = xf_Error(xf_SolveBlockPLUT(R00, n0, sr, P, X0, AddFlag, X0));
      if (ierr != xf_OK) return ierr; 
    }
    
  } // ie
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_LineJacobi
static int
xf_Jacobian_SolveM_LineJacobi(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X, 
                              enum xfe_AddType AddFlag, enum xfe_Bool TransposeFlag)
{  
  // X = M^{-1} * X  or X = M^{-T} * X  if TransposeFlag == True
  int ierr, egrp, negrp, r, iLine;
  real *XT;
  xf_LineSet *LineSet;
  xf_Mesh *Mesh;
  
  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  if ((LineSet = R_U->LineSet) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // allocate XT
  for (egrp=0, r=0; egrp<negrp; egrp++) r = max(r,X->GenArray[egrp].r);
  ierr = xf_Error(xf_Alloc((void **) &XT, r, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // loop over lines
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    
    // X_line @= M_line^{-1} * X_line,  @ = AddFlag, _line = elems on iLine
    ierr = xf_Error(xf_Jacobian_SolveM_OneLine(LineSet->Line+iLine, R_U, X, 
                                               AddFlag, TransposeFlag, XT));
    if (ierr != xf_OK) return ierr;
    
  } // iLine
  
  // release memory
  xf_Release( (void *) XT);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_LineGS
static int
xf_Jacobian_SolveM_LineGS(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X, 
                          enum xfe_AddType AddFlag, enum xfe_Bool TransposeFlag)
{  
  // X = M^{-1} * X
  int ierr, k, egrp, elem, face, negrp, r, sr;
  int egrpN, elemN, faceN, pface, n0, n1;
  int iLine, iLineN, fi, ie;
  real *R01, *R10, *X0, *X1, *XT;
  xf_LineSet *LineSet;
  xf_Mesh *Mesh;
  
  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup;
  sr    = R_U->StateRank;
  
  // only support Set for now
  if (AddFlag != xfe_Set) return xf_Error(xf_NOT_SUPPORTED);
  
  if ((LineSet = R_U->LineSet) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // allocate XT
  for (egrp=0, r=0; egrp<negrp; egrp++) r = max(r,X->GenArray[egrp].r);
  ierr = xf_Error(xf_Alloc((void **) &XT, r, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // loop over lines
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    
    /*------------------------------------------------*/
    /* X_line -= M_off * X_off                        */
    /* _off = elems not on line and with lower line # */
    /*------------------------------------------------*/
    
    pface = -1; // initialize previous face number
    
    // loop over elements in a line
    for (ie=0; ie<LineSet->Line[iLine].nelem; ie++){
      
      // pull off info of elem ie on line
      egrp = LineSet->Line[iLine].egrp[ie];
      elem = LineSet->Line[iLine].elem[ie];
      face = LineSet->Line[iLine].face[ie];
      
      n0 = xf_Jacobian_n(R_U,egrp,elem);   // rank/sr for egrp,elem
      X0 = X->GenArray[egrp].rValue[elem]; // X(elem)
      
      // loop over faces of elem
      for (fi=0; fi<Mesh->ElemGroup[egrp].nFace[elem]; fi++){
        
        // do not consider previous or next elements
        if ((fi==face) || (fi==pface)) continue;
        
        egrpN = R_U->egrpN[egrp][elem][fi];
        elemN = R_U->elemN[egrp][elem][fi];
        
        // do not consider halo or boundary elements
        if ((egrpN >= negrp) || (egrpN < 0)) continue;
        
        // do not consider neighbors with iLineN >= iLine
        iLineN = LineSet->Elem2Line[egrpN][elemN];
        if (iLineN >= iLine) continue;
        
        R01 = R_U->Value[egrp ][elem ][1+fi];   // R_U(elem, elemN)
        n1  = xf_Jacobian_n(R_U,egrpN,elemN);   // rank/sr for egrpN,elemN
        X1  = X->GenArray[egrpN].rValue[elemN]; // X(elemN)
        
        // X0 -= R01*X1
        if (!TransposeFlag){
          ierr = xf_Error(xf_BlockMxV(R01, n0, sr, n1, X1, xfe_Sub, X0));
          if (ierr != xf_OK) return ierr;
        }
        else{
          faceN = R_U->faceN[egrp][elem][fi];
          R10 = R_U->Value[egrpN][elemN][1+faceN];   // R_U(elemN, elem)
          ierr = xf_Error(xf_BlockMTxV(R10, n0, sr, n1, X1, xfe_Sub, X0));
          if (ierr != xf_OK) return ierr;
        }
        
      } // fi
      
      // set previous face number
      if (face >= 0) pface = R_U->faceN[egrp][elem][face];
    } // ie
    
    /*--------------------------------*/
    /* X_line = M_line^{-1} * X_line  */
    /*  _line = elems on iLine        */
    /*--------------------------------*/
    
    ierr = xf_Error(xf_Jacobian_SolveM_OneLine(LineSet->Line+iLine, R_U, X, 
                                               xfe_Set, TransposeFlag, XT));
    if (ierr != xf_OK) return ierr;
    
  } // iLine
  
  // release memory
  xf_Release( (void *) XT);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM
static int
xf_Jacobian_SolveM(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X,
                   enum xfe_AddType AddFlag, enum xfe_Bool TransposeFlag,
                   xf_SolverData *SolverData)
{
  /*
   PURPOSE:
   
   Sets 
   X @= M^{-1}*X
   
   where @ is the AddFlag.  R_U must be preconditioned: R_U = M + N.
   Thus, M is the preconditioner part of R_U.
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix (contains preconditioner info)
   X   : Input that gets multiplied by M^{-1}
   AddFlag: one of xfe_Set, xfe_Add, etc.
   TransposeFlag: if True, M^{-T} is used instead of M^{-1}
   SolverData : solver data structure
   
   OUTPUTS:
   
   Y   : Resulting vector.
   
   RETURNS:
   
   Error Code
   
   */  
  
  int ierr;
  
  switch (R_U->Preconditioner){
  case xfe_PreconditionerNone:
    return xf_Error(xf_INPUT_ERROR); // R_U is not preconditioned
    break;
  case xfe_PreconditionerBlockJacobi:
    ierr = xf_Error(xf_Jacobian_SolveM_BlockJacobi(All, R_U, X, AddFlag,
						   TransposeFlag));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerBlockJacobiLean:
    ierr = xf_Error(xf_Jacobian_SolveM_BlockJacobiLean(All, R_U, X, AddFlag,
						       TransposeFlag, SolverData));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerLineJacobi:
    ierr = xf_Error(xf_Jacobian_SolveM_LineJacobi(All, R_U, X, AddFlag,
						  TransposeFlag));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerLineGS:
    ierr = xf_Error(xf_Jacobian_SolveM_LineGS(All, R_U, X, AddFlag,
					      TransposeFlag));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerILU0:
    ierr = xf_Error(xf_Jacobian_SolveM_ILU0(All, R_U, X, AddFlag,
					    TransposeFlag));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultM_BlockJacobi
static int
xf_Jacobian_MultM_BlockJacobi(xf_All *All, xf_JacobianMatrix *R_U,
                              xf_Vector *X, enum xfe_AddType AddFlag,
                              enum xfe_Bool TransposeFlag, xf_Vector *Y)
{  
  int ierr, egrp, elem, nn;
  xf_Mesh *Mesh;
  Mesh = All->Mesh;
  
  // loop over elements, multiply diagonal blocks
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      nn = xf_Jacobian_n(R_U,egrp,elem);
      if (!TransposeFlag){
        ierr = xf_Error(xf_BlockPLUMxV(R_U->Value[egrp][elem][0], nn,
                                       R_U->StateRank, R_U->P[egrp][elem],
                                       X->GenArray[egrp].rValue[elem], AddFlag,
                                       Y->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr;
      }
      else{
        ierr = xf_Error(xf_BlockPLUMTxV(R_U->Value[egrp][elem][0], nn,
                                        R_U->StateRank, R_U->P[egrp][elem],
                                        X->GenArray[egrp].rValue[elem], AddFlag,
                                        Y->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr;
      }
      
    } // elem
  } // egrp
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultM_LineJacobi
static int
xf_Jacobian_MultM_LineJacobi(xf_All *All, xf_JacobianMatrix *R_U,
                             xf_Vector *X, enum xfe_AddType AddFlag, 
                             enum xfe_Bool TransposeFlag, xf_Vector *Y)
{  
  // Y @= M * X  (or M^T * X if TransposeFlag)
  int ierr, k, egrp, elem, face, negrp, r, sr;
  int egrpN, elemN, faceN, fi, pface, n0, n1;
  int iLine, ie, *P;
  enum xfe_AddType AddFlag2;
  real *R00, *R01, *R10, *X0, *X1, *XT, *XS, *Y0, *Y1;
  xf_LineSet *LineSet;
  xf_Mesh *Mesh;
  
  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup;
  sr    = R_U->StateRank;
  
  if ((LineSet = R_U->LineSet) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  // allocate XT and XS
  for (egrp=0, r=0; egrp<negrp; egrp++) r = max(r,X->GenArray[egrp].r);
  ierr = xf_Error(xf_Alloc((void **) &XT, r, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &XS, r, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // loop over lines
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    
    
    /*--------------------------------------------------*/
    /*     Y @= M   * X                                 */ 
    /*  or Y @= M^T * X  (if TransposeFlag == True)     */ 
    /*                                                  */
    /* i.e. Y @= L * U * X, where (3x3 e.g.)            */
    /*                                                  */
    /*     L  = [    I         0       0]               */
    /*          [ C0*A0'^-1    I       0]               */
    /*          [   0       C1*A1'^-1  I]               */
    /*                                                  */
    /* and                                              */
    /*                                                  */
    /*     U  = [ LU(A0)      B0       0  ]             */
    /*          [   0      LU(A1')    B1  ]             */
    /*          [   0         0    LU(A2')]             */
    /*                                                  */
    /*  Note, M^T = U^T * L^T                           */
    /*--------------------------------------------------*/
    
    
    // loop over elements in a line (starting from end)
    for (ie=LineSet->Line[iLine].nelem-1; ie>=0; ie--){
      
      // pull off info of elem ie on line
      egrp = LineSet->Line[iLine].egrp[ie];
      elem = LineSet->Line[iLine].elem[ie];
      face = LineSet->Line[iLine].face[ie];
      
      R00 = R_U->Value[egrp][elem][0];      // R_U(elem,elem); PLU factored
      n0  = xf_Jacobian_n(R_U,egrp,elem);   // rank/sr for egrp,elem
      P   = R_U->P[egrp][elem];             // permutation vector for elem
      X0  = X->GenArray[egrp].rValue[elem]; // X(elem)
      Y0  = Y->GenArray[egrp].rValue[elem]; // Y(elem)
      
      if (!TransposeFlag){
        // XT = R00 * X0
        ierr = xf_Error(xf_BlockPLUMxV(R00, n0, sr, P, X0, xfe_Set, XT));
        if (ierr != xf_OK) return ierr;
      }
      else{
        // XT = R00^T * X0
        ierr = xf_Error(xf_BlockPLUMTxV(R00, n0, sr, P, X0, xfe_Set, XT));
        if (ierr != xf_OK) return ierr;
      }
      
      if (ie < (LineSet->Line[iLine].nelem-1)){ // if not last elem in line
        
        // pull off info of next elem
        egrpN = R_U->egrpN[egrp][elem][face];
        elemN = R_U->elemN[egrp][elem][face];
        faceN = R_U->faceN[egrp][elem][face]; // elem is on other side
        
        R01 = R_U->Value[egrp ][elem ][1+face ]; // R_U(elem ,elemN)
        R10 = R_U->Value[egrpN][elemN][1+faceN]; // R_U(elemN,elem )
        n1  = xf_Jacobian_n(R_U,egrpN,elemN);    // rank/sr for egrpN,elemN
        X1  = X->GenArray[egrpN].rValue[elemN];  // X(elemN)
        Y1  = Y->GenArray[egrpN].rValue[elemN];  // Y(elemN)
        
        if (!TransposeFlag){
          // XT += R01 * X1
          ierr = xf_Error(xf_BlockMxV(R01, n0, sr, n1, X1, xfe_Add, XT));
          if (ierr != xf_OK) return ierr;
          
          // XS = R00^{-1} * XT
          ierr = xf_Error(xf_SolveBlockPLU(R00, n0, sr, P, XT, xfe_Set, XS));
          if (ierr != xf_OK) return ierr;
          
          // Y1 &= R10 * XS,  & is AddFlag2
          ierr = xf_Error(xf_BlockMxV(R10, n1, sr, n0, XS, AddFlag2, Y1));
          if (ierr != xf_OK) return ierr;
        }
        else{
          
          // XT += R10^T * X1
          ierr = xf_Error(xf_BlockMTxV(R10, n0, sr, n1, X1, xfe_Add, XT));
          if (ierr != xf_OK) return ierr;
          
          // XS = R00^{-T} * XT
          ierr = xf_Error(xf_SolveBlockPLUT(R00, n0, sr, P, XT, xfe_Set, XS));
          if (ierr != xf_OK) return ierr;
          
          // Y1 &= R01^T * XS,  & is AddFlag2
          ierr = xf_Error(xf_BlockMTxV(R01, n1, sr, n0, XS, AddFlag2, Y1));
          if (ierr != xf_OK) return ierr;
        }
        
      }
      
      // Y0 @= XT,  @= AddFlag
      xf_V_Add(XT, n0*sr, AddFlag, Y0);
      
    } // ie
    
  } // iLine
  
  // release memory
  xf_Release( (void *) XT);
  xf_Release( (void *) XS);
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultM
static int
xf_Jacobian_MultM(xf_All *All, xf_JacobianMatrix *R_U,
                  xf_Vector *X, enum xfe_AddType AddFlag,
                  enum xfe_Bool TransposeFlag, 
                  xf_SolverData *SolverData, xf_Vector *Y)
{
  /*
   PURPOSE:
   
   Sets 
   Y @= M*X
   
   where @= is one of {=, +=, -=, =-} according to AddFlag.  R_U must
   be preconditioned: R_U = M + N.  Thus, M is the preconditioner part
   of R_U.
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix (contains preconditioner info)
   X   : Input that gets multiplied by M
   AddFlag: one of xfe_Set, xfe_Add, etc.
   TransposeFlag: if True, M^T is used instead of M
   SolverData : solver data structure
   
   OUTPUTS:
   
   Y   : Resulting vector.
   
   RETURNS:
   
   Error Code
   
   */  
  
  int ierr;
  
  switch (R_U->Preconditioner){
  case xfe_PreconditionerBlockJacobi:
    ierr = xf_Error(xf_Jacobian_MultM_BlockJacobi(All, R_U, X, AddFlag, 
						  TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerBlockJacobiLean:
    ierr = xf_Error(xf_Jacobian_MultM_BlockJacobiLean(All, R_U, X, AddFlag, TransposeFlag, 
						      SolverData, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerLineJacobi:
    ierr = xf_Error(xf_Jacobian_MultM_LineJacobi(All, R_U, X, AddFlag, 
						 TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerILU0:
    ierr = xf_Error(xf_Jacobian_MultM_ILU0(All, R_U, X, AddFlag, 
					   TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultN_BlockJacobi
static int
xf_Jacobian_MultN_BlockJacobi(xf_All *All, xf_JacobianMatrix *R_U,
                              xf_Vector *X, enum xfe_AddType AddFlag,
                              enum xfe_Bool TransposeFlag, xf_Vector *Y)
{  
  int ierr, egrp, negrp, elem, face;
  int egrpN, elemN, faceN, k, n, nN;
  enum xfe_AddType AddFlag2; 
  real *rY;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  // loop over elements, add to Y
  for (egrp=0; egrp<negrp; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      n = xf_Jacobian_n(R_U,egrp,elem);
      if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){ // Set to 0 if necessary
        rY = Y->GenArray[egrp].rValue[elem];
        for (k=0; k<n*R_U->StateRank; k++) rY[k] = 0.0;
      }
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
        egrpN = R_U->egrpN[egrp][elem][face];
        elemN = R_U->elemN[egrp][elem][face];
        if ((egrpN < 0) || (egrpN >= negrp)) continue; // X on halo is not here yet
	nN = xf_Jacobian_n(R_U,egrpN,elemN);
        if (!TransposeFlag){
          ierr = xf_Error(xf_BlockMxV(R_U->Value[egrp][elem][1+face], n,
                                      R_U->StateRank, nN, 
                                      X->GenArray[egrpN].rValue[elemN], AddFlag2,
                                      Y->GenArray[egrp].rValue[elem]));
          if (ierr != xf_OK) return ierr;
        }
        else{
          faceN = R_U->faceN[egrp][elem][face];
          ierr = xf_Error(xf_BlockMTxV(R_U->Value[egrpN][elemN][1+faceN], 
                                       n, R_U->StateRank, nN,
                                       X->GenArray[egrpN].rValue[elemN], AddFlag2,
                                       Y->GenArray[egrp].rValue[elem]));
          if (ierr != xf_OK) return ierr;
        }
      } // face
    } // elem
  } // egrp
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultN_LineJacobi
static int
xf_Jacobian_MultN_LineJacobi(xf_All *All, xf_JacobianMatrix *R_U,
                             xf_Vector *X, enum xfe_AddType AddFlag,
                             enum xfe_Bool TransposeFlag, xf_Vector *Y)
{  
  int ierr, k, egrp, elem, negrp;
  int egrpN, elemN, faceN, fi, n0, n1;
  int sr, *R_UegrpN, *R_UelemN;
  enum xfe_AddType AddFlag2; 
  real *X1, *Y0, *R01, *R10;
  xf_LineSet *LineSet;
  xf_Mesh *Mesh;
  
  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup;
  sr    = R_U->StateRank;
  
  if ((LineSet = R_U->LineSet) == NULL) return xf_Error(xf_INPUT_ERROR);
  if (LineSet->Face2M == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  for (egrp=0; egrp<negrp; egrp++){
        
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      R_UegrpN = R_U->egrpN[egrp][elem];
      R_UelemN = R_U->elemN[egrp][elem];
    
      n0 = xf_Jacobian_n(R_U,egrp,elem);   // rank/sr for egrp,elem

      Y0 = Y->GenArray[egrp].rValue[elem]; // Y(elem)
      if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg))
        for (k=0; k<n0*sr; k++) Y0[k] = 0.0;
      
      // loop over faces of elem
      for (fi=0; fi<Mesh->ElemGroup[egrp].nFace[elem]; fi++){
        
        if (LineSet->Face2M[egrp][elem][fi] != 0) continue;
        
        egrpN = R_UegrpN[fi];
        elemN = R_UelemN[fi];
        
        R01 = R_U->Value[egrp ][elem ][1+fi];   // R_U(elem, elemN)
        n1  = xf_Jacobian_n(R_U,egrpN,elemN);   // rank/sr for egrpN,elemN
        X1  = X->GenArray[egrpN].rValue[elemN]; // X(elemN)
        
        if (!TransposeFlag){
          // Y0 &= R01*X1, where & is AddFlag2
          ierr = xf_Error(xf_BlockMxV(R01, n0, sr, n1, X1, AddFlag2, Y0));
          if (ierr != xf_OK) return ierr;
        }
        else{
          // or, if TransposeFlag, Y0 &= R10^T*X1
          faceN = R_U->faceN[egrp][elem][fi];
          R10 = R_U->Value[egrpN][elemN][1+faceN];  // R_U(elemN, elem)
          ierr = xf_Error(xf_BlockMTxV(R10, n0, sr, n1, X1, AddFlag2, Y0));
          if (ierr != xf_OK) return ierr;
        }
        
      } // fi
    } // elem
  } // egrp
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultN_LineGS
static int
xf_Jacobian_MultN_LineGS(xf_All *All, xf_JacobianMatrix *R_U,
                         xf_Vector *X, enum xfe_AddType AddFlag,
                         enum xfe_Bool TransposeFlag, xf_Vector *Y)
{  
  int ierr, k, egrp, elem, face, negrp;
  int egrpN, elemN, faceN, fi, pface, n0, n1;
  int iLine, iLineN, ie, sr;
  enum xfe_AddType AddFlag2; 
  real *X1, *Y0, *R01, *R10;
  xf_LineSet *LineSet;
  xf_Mesh *Mesh;
  
  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup;
  sr    = R_U->StateRank;
  
  if ((LineSet = R_U->LineSet) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  // loop over lines
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    
    pface = -1; // initialize previous face number
    
    // loop over elements in a line
    for (ie=0; ie<LineSet->Line[iLine].nelem; ie++){
      
      // pull off info of elem ie on line
      egrp = LineSet->Line[iLine].egrp[ie];
      elem = LineSet->Line[iLine].elem[ie];
      face = LineSet->Line[iLine].face[ie];
      
      n0 = xf_Jacobian_n(R_U,egrp,elem);   // rank/sr for egrp,elem
      Y0 = Y->GenArray[egrp].rValue[elem]; // Y(elem)
      
      if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg))
        for (k=0; k<n0*sr; k++) Y0[k] = 0.0;
      
      // loop over faces of elem
      for (fi=0; fi<Mesh->ElemGroup[egrp].nFace[elem]; fi++){
        
        // do not consider previous or next elements
        if ((fi==face) || (fi==pface)) continue;
        
        egrpN = R_U->egrpN[egrp][elem][fi];
        elemN = R_U->elemN[egrp][elem][fi];
        
        // do not consider halo or boundary elements
        if ((egrpN >= negrp) || (egrpN < 0)) continue;
        
        // do not consider neighbors with iLineN < iLine
        iLineN = LineSet->Elem2Line[egrpN][elemN];
        if (iLineN < iLine) continue;
        
        R01 = R_U->Value[egrp ][elem ][1+fi];   // R_U(elem, elemN)
        n1  = xf_Jacobian_n(R_U,egrpN,elemN);   // rank/sr for egrpN,elemN
        X1  = X->GenArray[egrpN].rValue[elemN]; // X(elemN)
        
        if (!TransposeFlag){
          // Y0 &= R01*X1, where & is AddFlag2
          ierr = xf_Error(xf_BlockMxV(R01, n0, sr, n1, X1, AddFlag2, Y0));
          if (ierr != xf_OK) return ierr;
        }
        else{
          // or, if TransposeFlag, Y0 &= R10^T*X1
          faceN = R_U->faceN[egrp][elem][fi];
          R10 = R_U->Value[egrpN][elemN][1+faceN];  // R_U(elemN, elem)
          ierr = xf_Error(xf_BlockMTxV(R10, n0, sr, n1, X1, AddFlag2, Y0));
          if (ierr != xf_OK) return ierr;
        }
        
      } // fi
      
      // set previous face number
      if (face >= 0) pface = R_U->faceN[egrp][elem][face]; 
      
    } // ie
  } // iLine
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultN_Halo
static int
xf_Jacobian_MultN_Halo(xf_All *All, xf_JacobianMatrix *R_U,
                       xf_Vector *X, enum xfe_AddType AddFlag,
                       enum xfe_Bool TransposeFlag, 
                       xf_SolverData *SolverData, xf_Vector *Y)
{  
  int ierr, egrp, negrp, elem, face;
  int egrpN, elemN, faceN, n, nN;
  enum xfe_Bool LeanFlag;
  enum xfe_AddType AddFlag2; 
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  if (Mesh->ParallelInfo == NULL) return xf_OK; // Mesh is not parallel
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  // check if preconditioner is lean
  ierr = xf_Error(xf_PreconditionerLeanCheck(R_U->Preconditioner, &LeanFlag));
  if (ierr != xf_OK) return ierr;
  
  // if LeanFlag, call lean multiplication and return immediately
  if (LeanFlag)
    return xf_Error(xf_Jacobian_MultN_BlockJacobiLean(All, R_U, X, AddFlag, TransposeFlag, 
                                                      SolverData, xfe_True, Y));
  
  
  // loop over elements, finish multiplication by handling halo
  for (egrp=0; egrp<negrp; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      n = xf_Jacobian_n(R_U,egrp,elem);
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
        egrpN = R_U->egrpN[egrp][elem][face];
        elemN = R_U->elemN[egrp][elem][face];
        if (egrpN < negrp) continue; // only consider halo
        
	nN = xf_Jacobian_n(R_U,egrpN,elemN);
        if (!TransposeFlag){
          ierr = xf_Error(xf_BlockMxV(R_U->Value[egrp][elem][1+face], n,
                                      R_U->StateRank, nN, 
                                      X->GenArray[egrpN].rValue[elemN], AddFlag2,
                                      Y->GenArray[egrp].rValue[elem]));
          if (ierr != xf_OK) return ierr;
        }
        else{
          faceN = R_U->faceN[egrp][elem][face];
          ierr = xf_Error(xf_BlockMTxV(R_U->Value[egrpN][elemN][1+faceN], 
                                       n, R_U->StateRank, nN, 
                                       X->GenArray[egrpN].rValue[elemN], AddFlag2,
                                       Y->GenArray[egrp].rValue[elem]));
          if (ierr != xf_OK) return ierr;
        }
      } // face
    } // elem
  } // egrp
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultN
static int
xf_Jacobian_MultN(xf_All *All, xf_JacobianMatrix *R_U,
                  xf_Vector *X, enum xfe_AddType AddFlag,
                  enum xfe_Bool TransposeFlag, 
                  xf_SolverData *SolverData, xf_Vector *Y)
{
  /*
   PURPOSE:
   
   Sets 
   Y @= N*X
   
   where @= is one of {=, +=, -=, =-} according to AddFlag.  R_U must
   be preconditioned: R_U = M + N.  Thus, N is the non-preconditioner
   part of R_U.
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix (contains preconditioner info)
   X   : Input that gets multiplied by N
   AddFlag: one of xfe_Set, xfe_Add, etc.
   TransposeFlag : if True, R_U^T is used instead of R_U
   SolverData : solver data structure
   
   OUTPUTS:
   
   Y   : Resulting vector.
   
   RETURNS:
   
   Error Code
   
   */  
  
  int ierr;
  enum xfe_AddType AddFlag2; 
  
  if (!SolverData->SkipParallelExchange){
    // begin communication of halo data: X
    ierr = xf_Error(xf_HaloExchangeVectorBegin(X));
    if (ierr != xf_OK) return ierr;
  }
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  switch (R_U->Preconditioner){
  case xfe_PreconditionerBlockJacobi:
    ierr = xf_Error(xf_Jacobian_MultN_BlockJacobi(All, R_U, X, AddFlag, 
						  TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerBlockJacobiLean:
    ierr = xf_Error(xf_Jacobian_MultN_BlockJacobiLean(All, R_U, X, AddFlag, TransposeFlag, 
						      SolverData, xfe_False, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerLineJacobi:
    ierr = xf_Error(xf_Jacobian_MultN_LineJacobi(All, R_U, X, AddFlag, 
						 TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerLineGS:
    ierr = xf_Error(xf_Jacobian_MultN_LineGS(All, R_U, X, AddFlag, 
					     TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerILU0:
    ierr = xf_Error(xf_Jacobian_MultN_ILU0(All, R_U, X, AddFlag, 
					   TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  
  if (!SolverData->SkipParallelExchange){
    // end communication of halo data: X
    ierr = xf_Error(xf_HaloExchangeVectorEnd(X));
    if (ierr != xf_OK) return ierr;
  
    // finish multiplying halo
    ierr = xf_Error(xf_Jacobian_MultN_Halo(All, R_U, X, AddFlag2, 
					   TransposeFlag, SolverData, Y));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_Mult_None
static int
xf_Jacobian_Mult_None(xf_All *All, xf_JacobianMatrix *R_U,
                      xf_Vector *X, enum xfe_AddType AddFlag,
                      enum xfe_Bool TransposeFlag, xf_Vector *Y)
{  
  int ierr, egrp, negrp, elem, face;
  int egrpN, elemN, faceN, n, nN;
  enum xfe_AddType AddFlag2; 
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  // loop over elements, apply R_U to X
  for (egrp=0; egrp<negrp; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      n = xf_Jacobian_n(R_U,egrp,elem);
      if (!TransposeFlag){
        ierr = xf_Error(xf_BlockMxV(R_U->Value[egrp][elem][0], n,
                                    R_U->StateRank, n,
                                    X->GenArray[egrp].rValue[elem], AddFlag,
                                    Y->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr;
      }
      else{
        ierr = xf_Error(xf_BlockMTxV(R_U->Value[egrp][elem][0], n,
                                     R_U->StateRank, n,
                                     X->GenArray[egrp].rValue[elem], AddFlag,
                                     Y->GenArray[egrp].rValue[elem]));
        if (ierr != xf_OK) return ierr;
      }
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){
        egrpN = R_U->egrpN[egrp][elem][face];
        elemN = R_U->elemN[egrp][elem][face];
        if ((egrpN < 0) || (egrpN >= negrp)) continue; // X on halo is not here yet
	nN = xf_Jacobian_n(R_U,egrpN,elemN);
        if (!TransposeFlag){
          ierr = xf_Error(xf_BlockMxV(R_U->Value[egrp][elem][1+face], n,
                                      R_U->StateRank, nN, 
                                      X->GenArray[egrpN].rValue[elemN], AddFlag2,
                                      Y->GenArray[egrp].rValue[elem]));
          if (ierr != xf_OK) return ierr;
        }
        else{
          faceN = R_U->faceN[egrp][elem][face];
          ierr = xf_Error(xf_BlockMTxV(R_U->Value[egrpN][elemN][1+faceN], 
                                       n, R_U->StateRank, nN, 
                                       X->GenArray[egrpN].rValue[elemN], AddFlag2,
                                       Y->GenArray[egrp].rValue[elem]));
          if (ierr != xf_OK) return ierr;
        }
      } // face
    } // elem
  } // egrp
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_Mult
int
xf_Jacobian_Mult(xf_All *All, xf_JacobianMatrix *R_U,
                 xf_Vector *X, enum xfe_AddType AddFlag, 
                 enum xfe_Bool TransposeFlag, 
                 xf_SolverData *SolverData, xf_Vector *Y)
{
  
  int ierr;
  enum xfe_AddType AddFlag2; 
  
  // begin communication of halo data: X
  ierr = xf_Error(xf_HaloExchangeVectorBegin(X));
  if (ierr != xf_OK) return ierr;
  
  // if addflag is set or neg, initialize Y=0
  if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
    ierr = xf_Error(xf_SetZeroVector(Y));
    if (ierr != xf_OK) return ierr;
  }
  
  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  
  switch (R_U->Preconditioner){
  case xfe_PreconditionerNone:
    ierr = xf_Error(xf_Jacobian_Mult_None(All, R_U, X, AddFlag, TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerBlockJacobi:
    ierr = xf_Error(xf_Jacobian_MultM_BlockJacobi(All, R_U, X, AddFlag, TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Jacobian_MultN_BlockJacobi(All, R_U, X, AddFlag2, TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerBlockJacobiLean:
    ierr = xf_Error(xf_Jacobian_MultM_BlockJacobiLean(All, R_U, X, AddFlag, TransposeFlag, 
						      SolverData, Y));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Jacobian_MultN_BlockJacobiLean(All, R_U, X, AddFlag2, TransposeFlag, 
						      SolverData, xfe_False, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerLineJacobi:
    ierr = xf_Error(xf_Jacobian_MultM_LineJacobi(All, R_U, X, AddFlag, TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Jacobian_MultN_LineJacobi(All, R_U, X, AddFlag2, TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_PreconditionerILU0:
    ierr = xf_Error(xf_Jacobian_MultM_ILU0(All, R_U, X, AddFlag, TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Jacobian_MultN_ILU0(All, R_U, X, AddFlag2, TransposeFlag, Y));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  
  // end communication of halo data: X
  ierr = xf_Error(xf_HaloExchangeVectorEnd(X));
  if (ierr != xf_OK) return ierr;
  
  // finish multiplying halo
  ierr = xf_Error(xf_Jacobian_MultN_Halo(All, R_U, X, AddFlag2, TransposeFlag, SolverData, Y));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}





/******************************************************************/
//   FUNCTION Definition: xf_GeneralLinearStep
static int
xf_GeneralLinearStep(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X, 
                     xf_Vector *R, xf_Vector *InVtemp, enum xfe_AddType AddFlag, 
                     enum xfe_Bool TransposeFlag, real omega, 
                     xf_SolverData *SolverData, xf_Vector *Y )
{
  /*
   PURPOSE:
   
   Sets 
   Y @= omega*M^{-1} (R_U * X + R)
   
   where @= is one of {=, +=, -=, =-} according to AddFlag.  M is the
   preconditioner as indicated by R_U->Preconditioner (none means M
   is effectively the identity).  R can be NULL.  X can also be NULL.
   Both R and X NULL would be silly.  Y and [X or R] can be the same
   vectors.
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix (contains preconditioner info, if any)
   X   : Input that gets multiplied by R_U (can be NULL)
   R   : Input that is added (see above)
   InVtemp : a temporary-storage vector so don't have to look for one 
   here (can be NULL)
   AddFlag: one of xfe_Set, xfe_Add, etc.
   TransposeFlag: if True, R_U^T is used instead of R_U, and
   M^{-T} is used instead of M^{-1}
   omega : under-relaxation factor
   SolverData : solver data structure
   
   OUTPUTS:
   
   Y   : Resulting vector.  Can be same as X or R
   
   RETURNS:
   
   Error Code
   
   */  
  
  int ierr;
  enum xfe_AddType XFlag;
  enum xfe_Bool Preconditioned;
  xf_Vector *Vtemp;
  
  if ((R == NULL) && (X== NULL)) return xf_Error(xf_INPUT_ERROR);
  
  // check if R_U is preconditioned
  Preconditioned = (R_U->Preconditioner != xfe_PreconditionerNone);
  
  // find Vtemp if necessary
  if (InVtemp == NULL){
    ierr = xf_Error(xf_FindSimilarVector(All, Y, "GeneralLinearStep_Vtemp",
                                         xfe_False, xfe_True, NULL, &Vtemp, NULL));
    if (ierr != xf_OK) return ierr;
  }
  else
    Vtemp = InVtemp;
  
  // Set Vtemp to R and set XFlag = add or set depending on whether R is NULL
  if (R != NULL){
    // Vtemp = R
    ierr = xf_Error(xf_SetVector(R, xfe_Set, Vtemp));
    if (ierr != xf_OK) return ierr;
    XFlag = xfe_Add;
  }
  else XFlag = xfe_Set;
  
  
  // handle not preconditioned case first:  Y @= R_U*X + R
  if (!Preconditioned){
    
    if (X != NULL){
      // Vtemp &= R_U*X,  & is XFlag
      ierr = xf_Error(xf_Jacobian_Mult(All, R_U, X, XFlag, TransposeFlag, SolverData, Vtemp));
      if (ierr != xf_OK) return ierr;
    }
    
    // Y @= Vtemp,  @ is AddFlag
    ierr = xf_Error(xf_SetVector(Vtemp, AddFlag, Y));
    if (ierr != xf_OK) return ierr;
    
    return xf_OK; // Done, return
  }
  
  
  // At this point we're doing:  Y @= X + M^{-1}*(N*X + R)
  
  if (X != NULL){
    // Vtemp &= N*X,  & is XFlag
    ierr = xf_Error(xf_Jacobian_MultN(All, R_U, X, XFlag, TransposeFlag, SolverData, Vtemp));
    if (ierr != xf_OK) return ierr;
  }
  
  // Vtemp = M^{-1} * Vtemp
  ierr = xf_Error(xf_Jacobian_SolveM(All, R_U, Vtemp, xfe_Set, TransposeFlag, SolverData));
  if (ierr != xf_OK) return ierr;
  
  if (X != NULL){
    // Vtemp += X
    ierr = xf_Error(xf_SetVector(X, xfe_Add, Vtemp));
    if (ierr != xf_OK) return ierr;
  }
  
  
  // Y @= omega*Vtemp,  @ is AddFlag
  ierr = xf_Error(xf_VectorMultSet(Vtemp, omega, AddFlag, Y));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SolveLinearSystem_Iterative
static int
xf_SolveLinearSystem_Iterative(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *R, 
                               enum xfe_Bool TransposeFlag, int nIterIn, 
                               xf_SolverData *SolverData, xf_Vector *dU )
{
  /*
   
   PURPOSE:
   
   Solves linear system using a preconditioned iterative method
   (e.g. Jacobi or Gauss-Seidel).  The system is:
   
   R_U*dU + R = 0
   
   For this method, dU is initially set to 0, and then updated
   according to:
   
   dU^{m+1} = dU^{m} - M^{-1} [R_U*dU^{m} + R]
   
   where M is the preconditioner.  Note, a preconditioner is required
   for this function.  If we decompose R_U into R_U = M + N, the above
   expression can be simplified (fewer operations) into:
   
	 dU^{m+1} = -M^{-1} [N*dU^{m} + R]
   
   Note, this can also be derived from: M*dU^{m+1} + N*dU^{m} + R = 0
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix
   R   : Residual vector
   TransposeFlag : if True, R_U^T is used instead of R_U
   nIterIn : if >= 0, this # iterations is used
   SolverData : solver data structure
   
   
   OUTPUTS:
   
   dU : iteratively-obtained solution vector
   
   RETURNS:
   
   Error Code
   
   */
  int ierr, nIter, iiter;
  xf_Vector *Vtemp;
  enum xfe_Bool CheckHalt;
  enum xfe_Bool DoDomainBlockJacobi = xfe_False;
  enum xfe_PreconditionerType Preconditioner;
  real omega;

  // determine max number of linear iterations
  if (nIterIn < 0){
    ierr = xf_GetKeyValueInt(All->Param->KeyValue, "nIterLinear", &nIter);
    if (ierr != xf_OK) return ierr;
  }
  else nIter = nIterIn;  

  // special case of zero requested iterations (set dU = 0)
  if (nIter == 0){
    ierr = xf_Error(xf_SetZeroVector(dU));
    if (ierr != xf_OK) return ierr;
    return xf_OK;
  }

  // determine preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
                                     xfe_PreconditionerName, (int ) xfe_PreconditionerLast, 
                                     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;
  
  // need a preconditioner for iterative method
  if (Preconditioner == xfe_PreconditionerNone){
    xf_printf("Must specify a Preconditioner for the Jacobi linear solver.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  // precondition R_U
  ierr = xf_Error(xf_Jacobian_Precondition(All, R_U, Preconditioner, xfe_False));
  if (ierr != xf_OK) return ierr;

  
  // do we need to check for halt during iterations?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "CheckHaltLinear", &CheckHalt));
  if (ierr != xf_OK) return ierr;
  
  // Under-relaxation factor
  ierr = xf_GetKeyValueReal(All->Param->KeyValue, "IterativeUnderRelax", &omega);
  if (ierr != xf_OK) return ierr;

  
  /* set dU = - omega*M^{-1} R */
  
  // dU = -omega*R
  ierr = xf_Error(xf_VectorMultSet(R, omega, xfe_Neg, dU));
  if (ierr != xf_OK) return ierr;
  
  // dU = M^{-1}*dU
  ierr = xf_Error(xf_Jacobian_SolveM(All, R_U, dU, xfe_Set, TransposeFlag, SolverData));
  if (ierr != xf_OK) return ierr;
  
  
  // if only asking for one iteration, return
  if (nIter == 1) return xf_OK;
  
  // find Vtemp vector
  ierr = xf_Error(xf_FindSimilarVector(All, R, "Jacobi_Vtemp", xfe_False,
                                       xfe_True, NULL, &Vtemp, NULL));
  if (ierr != xf_OK) return ierr;

  // are we doing domain-decomposition block Jacobi?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "DD_DoDomainBlockJacobi", 
				     &DoDomainBlockJacobi));
  if (ierr != xf_OK) return ierr;

  // block Jacobi domain decomposition means not exchanging messages every iteration
  if (DoDomainBlockJacobi) SolverData->SkipParallelExchange = xfe_True;
  
  // perform remaining iterations (already did first one)
  for (iiter=1; iiter < nIter; iiter++){
    
    // dU -= omega*M^{-1} * (R_U*dU + R)
    ierr = xf_Error(xf_GeneralLinearStep(All, R_U, dU, R, Vtemp, xfe_Sub, 
                                         TransposeFlag, omega, SolverData, dU));
    if (ierr != xf_OK) return ierr;
    
    // check for halt
    if ((CheckHalt) && (xf_CheckUserHalt(NULL))) break;
    
  } // iiter

  SolverData->SkipParallelExchange = xfe_False;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_LinearStepCoarseCorrection
static int
xf_LinearStepCoarseCorrection(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X, 
                              xf_Vector *R, xf_Vector *InVtemp, 
                              enum xfe_Bool TransposeFlag, xf_SolverData *SolverData,
                              xf_Vector *Y )
{
  /*
   PURPOSE:
   
   Solves:
   
   R_U*Y - (R_U * X + R) = 0
   
   using a coarse correction procedure:
   
   LinR  = -(R_U * X + R)
	 LinRc = T^T * LinR       (T = prolongation operator, T^T = restriction) 
   solve R_Uc * Yc + LinRc = 0 using Iterative solver
   Y = T*Yc
   Y = Y - w*M^{-1} (R_U*Y + LinR) (nSmooth times)
   
	 
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix (contains preconditioner info, if any)
   X   : Input that gets multiplied by R_U (can be NULL)
   R   : Input that is added (see above)
   InVtemp : a temporary-storage vector so don't have to look for one 
   here (can be NULL)
   TransposeFlag: if True, R_U^T is used instead of R_U, and
   M^{-T} is used instead of M^{-1}
   SolverData : solver data structure
   
   OUTPUTS:
   
   Y   : Resulting vector.  Can be same as X or R
   
   RETURNS:
   
   Error Code
   
   */  
  
  int ierr, nIterCoarse;
  real omega, rnorm;
  int iSmooth, nSmooth;
  enum xfe_Bool AdditiveFlag = xfe_False;
  enum xfe_Bool DoDomainBlockJacobi = xfe_False;
  enum xfe_AddType XFlag;
  xf_Vector *LinR, *LinRc, *Yc, *Vtemp;
  xf_JacobianMatrix *R_Uc;
  
  if ((R == NULL) && (X== NULL)) return xf_Error(xf_INPUT_ERROR);
  
  if ((R_Uc = (xf_JacobianMatrix *) R_U->R_Uc) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // R_U needs to be preconditioned for smoothing
  if (R_U->Preconditioner == xfe_PreconditionerNone) return xf_Error(xf_INPUT_ERROR);
  
  // Is CoarseCorrection additive (Schwarz)?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "CoarseCorrectionIsAdditive", 
				     &AdditiveFlag));
  if (ierr != xf_OK) return ierr;
  

  // find LinR
  ierr = xf_Error(xf_FindSimilarVector(All, X, "LinearStepCoarseCorrection_LinR",
                                       xfe_False, xfe_True, NULL, &LinR, NULL));
  if (ierr != xf_OK) return ierr;
  
  
  /* LinR = -(R_U*X + R) */
  
  // LinR = -R
  if (R != NULL){
    ierr = xf_Error(xf_SetVector(R, xfe_Neg, LinR));
    if (ierr != xf_OK) return ierr;
  }
  else{
    ierr = xf_Error(xf_SetZeroVector(LinR));
    if (ierr != xf_OK) return ierr;
  }
  
  // LinR -= R_U*X
  if (X != NULL){
    ierr = xf_Error(xf_Jacobian_Mult(All, R_U, X, xfe_Sub, TransposeFlag, SolverData, LinR));
    if (ierr != xf_OK) return ierr;
  }
  
  
  if (LinR->Basis == NULL) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  // find LinRc
  ierr = xf_Error(xf_FindSimilarVectorHO(All, LinR, "LinearStepCoarseCorrection_LinRc", xfe_False, 
                                         xfe_True, NULL, R_Uc->Order, NULL, &LinRc));
  if (ierr != xf_OK) return ierr;
  
  /* LinRc = T^T * LinR */
  ierr = xf_Error(xf_ProjectVector(All, LinR, xfe_True, LinRc));
  if (ierr != xf_OK) return ierr;
  
  
  // find Yc
  ierr = xf_Error(xf_FindSimilarVector(All, LinRc, "LinearStepCoarseCorrection_Yc",
                                       xfe_True, xfe_True, NULL, &Yc, NULL));
  if (ierr != xf_OK) return ierr;
  
  // pull off number of coarse iterations
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "CoarseCorrectionIter", &nIterCoarse));
  if (ierr != xf_OK) return ierr;
  
  
  // do not do domain-decomposition block Jacobi on coarse level
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "DD_DoDomainBlockJacobi", 
				     &DoDomainBlockJacobi));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "DD_DoDomainBlockJacobi", 
				 "False"));
  if (ierr != xf_OK) return ierr;
  

  /* Solve R_Uc * Yc + LinRc = 0 using Iterative solver */
  ierr = xf_Error(xf_SolveLinearSystem_Iterative(All, R_Uc, LinRc, TransposeFlag, 
                                                 nIterCoarse, SolverData, Yc));
  if (ierr != xf_OK) return ierr;  

  // reset domain block Jacobi flag
  ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "DD_DoDomainBlockJacobi", 
				     DoDomainBlockJacobi));
  if (ierr != xf_OK) return ierr;


  
  /* Set Y = T*Yc */
  if (!AdditiveFlag){
    ierr = xf_Error(xf_ProjectVector(All, Yc, xfe_False, Y));
    if (ierr != xf_OK) return ierr;
  }
  
  
  /*   // TESTING */
  /*   ierr = xf_Error(xf_SetZeroVector(LinRc)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   ierr = xf_Error(xf_Jacobian_Mult(All, R_Uc, Yc, xfe_Add, TransposeFlag, LinRc)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   ierr = xf_Error(xf_VectorNorm(LinRc, 2, &rnorm)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   xf_printf("rnorm1 = %.10E\n", rnorm); */
  /*   ierr = xf_Error(xf_ProjectVector(All, Yc, xfe_False, Y)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   ierr = xf_Error(xf_SetZeroVector(LinR)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   ierr = xf_Error(xf_Jacobian_Mult(All, R_U, Y, xfe_Add, TransposeFlag, LinR)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   ierr = xf_Error(xf_ProjectVector(All, LinR, xfe_True, LinRc)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   ierr = xf_Error(xf_VectorNorm(LinRc, 2, &rnorm)); */
  /*   if (ierr != xf_OK) return ierr; */
  /*   xf_printf("rnorm2 = %.10E\n", rnorm); */
  
  /*   // TEMPORARY */
  /*   ierr = xf_Error(xf_SetZeroVector(Y)); */
  /*   if (ierr != xf_OK) return ierr; */
  
  
  // pull off number of smoothing iterations
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "CoarseCorrectionSmooth", &nSmooth));
  if (ierr != xf_OK) return ierr;
  
  // pull off under-relaxation factor
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CoarseCorrectionRelax", &omega));
  if (ierr != xf_OK) return ierr;


  // block Jacobi domain decomposition means not exchanging messages every iteration
  if (DoDomainBlockJacobi) SolverData->SkipParallelExchange = xfe_True;

  
  // find Vtemp if necessary
  if (InVtemp == NULL){
    ierr = xf_Error(xf_FindSimilarVector(All, X, "LinearStepCoarseCorrection_Vtemp",
                                         xfe_False, xfe_True, NULL, &Vtemp, NULL));
    if (ierr != xf_OK) return ierr;
  }
  else
    Vtemp = InVtemp;

  // If Additive, zero out Y
  if (AdditiveFlag){
    ierr = xf_Error(xf_SetZeroVector(Y)); 
    if (ierr != xf_OK) return ierr; 
  }
  
  // Vtemp = LinR
  ierr = xf_Error(xf_SetVector(LinR, xfe_Set, Vtemp));
  if (ierr != xf_OK) return ierr;
  
  /* Perform post-correction smoothing */
  
  for (iSmooth = 0; iSmooth < nSmooth; iSmooth++){
    
    if (iSmooth > 0){ // LinR = Vtemp
      ierr = xf_Error(xf_SetVector(Vtemp, xfe_Set, LinR));
      if (ierr != xf_OK) return ierr;
    }
    
    // Y = (1-w)*Y - w*M^{-1} * (R_U*Y + LinR)
    ierr = xf_Error(xf_Jacobian_MultN(All, R_U, Y, xfe_Add, TransposeFlag, SolverData, LinR));
    if (ierr != xf_OK) return ierr;
    
    // LinR = M^{-1} * LinR
    ierr = xf_Error(xf_Jacobian_SolveM(All, R_U, LinR, xfe_Set, TransposeFlag, SolverData));
    if (ierr != xf_OK) return ierr;
    
    // LinR *= omega
    ierr = xf_Error(xf_VectorMult(LinR, omega));
    if (ierr != xf_OK) return ierr;
    
    // Y *= (1-omega)
    ierr = xf_Error(xf_VectorMult(Y, 1.0-omega));
    if (ierr != xf_OK) return ierr;
    
    // Y -= LinR
    ierr = xf_Error(xf_SetVector(LinR, xfe_Sub, Y));
    if (ierr != xf_OK) return ierr;
  }

  // add coarse-grid correction if additive preconditioner
  if (AdditiveFlag){

    // Vtemp = T*Yc
    ierr = xf_Error(xf_ProjectVector(All, Yc, xfe_False, Vtemp));
    if (ierr != xf_OK) return ierr;

    // Y += Vtemp
    ierr = xf_Error(xf_SetVector(Vtemp, xfe_Add, Y));
    if (ierr != xf_OK) return ierr;

  }

  // reset SkipParallelExchange
  SolverData->SkipParallelExchange = xfe_True;
  
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_LinearStep
static int
xf_LinearStep(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X, 
              xf_Vector *R, xf_Vector *InVtemp, enum xfe_AddType AddFlag, 
              enum xfe_Bool TransposeFlag, xf_SolverData *SolverData, xf_Vector *Y )
{
  /*
   PURPOSE:
   
   Wrapper that calls either GeneralLinearStep, or
   LinearStepCoarseCorrection.
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix (contains preconditioner info, if any)
   X   : Input that gets multiplied by R_U (can be NULL)
   R   : Input that is added (see above)
   InVtemp : a temporary-storage vector so don't have to look for one 
   (can be NULL)
   AddFlag: one of xfe_Set, xfe_Add, etc.
   TransposeFlag: if True, R_U^T is used instead of R_U, and
   M^{-T} is used instead of M^{-1} (M is the preconditioner)
   SolverData : solver data structure
   
   OUTPUTS:
   
   Y   : Resulting vector.  Can be same as X or R
   
   RETURNS:
   
   Error Code
   
   */  
  int ierr;
  int nSmooth;
  xf_Vector *LinR;
  enum xfe_Bool CoarseCorrectionFlag;
  
  // determine if coarse correction is desired
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "CoarseCorrectionFlag", 
                                     &CoarseCorrectionFlag)); 
  if (ierr != xf_OK) return ierr;
  
  
  if (!CoarseCorrectionFlag){
    
    // pull off number of smoothing iterations
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nSmoothLinear", &nSmooth));
    if (ierr != xf_OK) return ierr;
    
    // Try multiple smoothing iterations 
    if ((nSmooth >= 2) && (AddFlag == xfe_Set)){
      
      // find LinR if necessary
      if (InVtemp == NULL){
        ierr = xf_Error(xf_FindSimilarVector(All, X, "LinearStep_Vtemp",
                                             xfe_True, xfe_True, NULL, &LinR, NULL));
        if (ierr != xf_OK) return ierr;
      }
      else
        LinR = InVtemp;
      
      if (nSmooth == 2){ // Below is optimized for two steps
        // LinR = N*X + R
        if (X != NULL){
          ierr = xf_Error(xf_Jacobian_MultN(All, R_U, X, xfe_Set, TransposeFlag, SolverData, LinR));
          if (ierr != xf_OK) return ierr;
        }
        if (R != NULL){
          ierr = xf_Error(xf_SetVector(R, xfe_Add, LinR));
          if (ierr != xf_OK) return ierr;
        }
        
        // LinR = M^{-1} * LinR
        ierr = xf_Error(xf_Jacobian_SolveM(All, R_U, LinR, xfe_Set, TransposeFlag, SolverData));
        if (ierr != xf_OK) return ierr;
        
        // Y = -N*LinR
        ierr = xf_Error(xf_Jacobian_MultN(All, R_U, LinR, xfe_Neg, TransposeFlag, SolverData, Y));
        if (ierr != xf_OK) return ierr;
        
        // Y += R
        if (R != NULL){
          ierr = xf_Error(xf_SetVector(R, xfe_Add, Y));
          if (ierr != xf_OK) return ierr;
        }
        
        // Y = M^{-1}*Y
        ierr = xf_Error(xf_Jacobian_SolveM(All, R_U, Y, xfe_Set, TransposeFlag, SolverData));
        if (ierr != xf_OK) return ierr;
        
        // Y += X
        if (X != NULL){
          ierr = xf_Error(xf_SetVector(X, xfe_Add, Y));
          if (ierr != xf_OK) return ierr;
        }
      }
      else{ // Below is for nSmooth > 2
        
        /* LinR = -(R_U*X + R) */
        
        // LinR = -R
        if (R != NULL){
          ierr = xf_Error(xf_SetVector(R, xfe_Neg, LinR));
          if (ierr != xf_OK) return ierr;
        }
        else{
          ierr = xf_Error(xf_SetZeroVector(LinR));
          if (ierr != xf_OK) return ierr;
        }
        
        // LinR -= R_U*X
        if (X != NULL){
          ierr = xf_Error(xf_Jacobian_Mult(All, R_U, X, xfe_Sub, TransposeFlag, SolverData, LinR));
          if (ierr != xf_OK) return ierr;
        }
        
        /* Solve R_U * Y + LinR = 0 using Iterative solver */
        ierr = xf_Error(xf_SolveLinearSystem_Iterative(All, R_U, LinR, TransposeFlag,
                                                       nSmooth, SolverData, Y));
        if (ierr != xf_OK) return ierr;
      }
      
      return xf_OK;
    }
    
    // Standard single linear step
    ierr = xf_Error(xf_GeneralLinearStep(All, R_U, X, R, InVtemp, xfe_Set,
                                         TransposeFlag, 1.0, SolverData, Y));
    if (ierr != xf_OK) return ierr;
  }
  else{
    // Coarse correction + smoothing
    if (AddFlag != xfe_Set) return xf_Error(xf_NOT_SUPPORTED);
    ierr = xf_Error(xf_LinearStepCoarseCorrection(All, R_U, X, R, InVtemp,
                                                  TransposeFlag, SolverData, Y));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SolveLinearSystem_Debug
static int
xf_SolveLinearSystem_Debug(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *R, 
                           enum xfe_Bool TransposeFlag, xf_SolverData *SolverData,
                           xf_Vector *dU )
{
  /*
   
   PURPOSE:
   
   Not an actual linear solver.  Contains several checks that are
   useful for debugging the preconditioned M and N functions.
   
   INPUTS:
   
   All : All structure
   R_U : Jacobian matrix
   R   : Residual vector
   TransposeFlag : if True, R_U^T is used instead of R_U
   SolverData : solver data structure
   
   dU  : a vector for testing
   
   OUTPUTS:
   
   None : see above
   
   RETURNS:
   
   Error Code
   
   */
  int ierr, i, j, k, r;
  real nv, ns, rnorm, diff;
  enum xfe_PreconditionerType Preconditioner;
  
  if (TransposeFlag) return xf_Error(xf_NOT_SUPPORTED);
  
  // determine preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
                                     xfe_PreconditionerName, (int ) xfe_PreconditionerLast, 
                                     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;
  
  // precondition R_U
  ierr = xf_Error(xf_Jacobian_Precondition(All, R_U, Preconditioner, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  // set dU = some vector (values nonzero, different, but irrelevant)
  ierr = xf_Error(xf_VectorRand(dU, 17));
  if (ierr != xf_OK) return ierr;
  
  // set R += R_U * dU
  ierr = xf_Error(xf_Jacobian_Mult(All, R_U, dU, xfe_Add, TransposeFlag, SolverData, R));
  if (ierr != xf_OK) return ierr;
  
  // print out R and Rnorm
  for (i=0; i<R->nArray; i++)
    for (j=0; j<R->GenArray[i].n; j++){
      r = ((R->GenArray[i].vr == NULL) ? R->GenArray[i].r : R->GenArray[i].vr[j]);
      for (k=0; k<r; k++)
        xf_printf("R[%d][%d][%d] = %.10E\n", i, j, k, R->GenArray[i].rValue[j][k]);
    }
  ierr = xf_Error(xf_VectorNorm(R, 2, &rnorm));
  if (ierr != xf_OK) return ierr;
  
  xf_printf("rnorm = %.15E\n", rnorm);
  
  // test M and M^{-1}
  
  // R = dU
  ierr = xf_Error(xf_SetVector(dU, xfe_Set, R));
  if (ierr != xf_OK) return ierr;
  
  // R = M^{-1}*R
  ierr = xf_Error(xf_Jacobian_SolveM(All, R_U, R, xfe_Set, TransposeFlag, SolverData));
  if (ierr != xf_OK) return ierr;
  
  // dU = M*R
  ierr = xf_Error(xf_Jacobian_MultM(All, R_U, R, xfe_Set, TransposeFlag, SolverData, dU));
  if (ierr != xf_OK) return ierr;
  
  // check dU
  for (i=0; i<dU->nArray; i++){
    nv = (real) dU->GenArray[i].n * dU->GenArray[i].r;
    ns = 0;
    for (j=0; j<dU->GenArray[i].n; j++){
      r = ((dU->GenArray[i].vr == NULL) ? dU->GenArray[i].r : dU->GenArray[i].vr[j]);
      for (k=0; k<r; k++){
        diff = dU->GenArray[i].rValue[j][k] - ns/nv;
        if (fabs(diff) > 1e-13)
          xf_printf("Warning: M*M^{-1} diff[%d][%d][%d] = %.10E\n", i, j, k, diff);
        ns += 1.0;
      }
    } // j
  } // i
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_GivensQRStep
static int
xf_GivensQRStep(int m, real **H, real *g, real *F)
{
  /*
   PURPOSE:
   
   Used in GMRES.  First applies m-1 Givens rotation form F to
   H[m]. Next, Calculates Givens rotation that sets the last two
   entries of H[m]: (H[m], H[m+1]) to (*, 0).  Applies this rotation
   to (g[m], g[m+1]).  Stores this rotation in F.  See Ref. [1] of
   SolveLinearSystem_GMRES for details.
   
   INPUTS:
   
   m : the column we're on
   H : Hessenberg matrix whose m'th column has just been added and
   contains m+2 entries (1 too many, hence the Givens rotations)
   g : associated right-hand side vector that also receives the rotations.
   
   OUTPUTS:
   
   H, g : modified as described above
   F    : stores newest Givens rotation in F[2*m] and F[2*m+1]
   
   RETURNS:
   
   Error Code
   */
  int i;
  real a, b, c, s, den;
  
  /* Apply m-1 rotations to column H[m] */
  for (i=0; i<m; i++){
    c = F[2*i+0]; 
    s = F[2*i+1];
    a = H[m][i]; 
    b = H[m][i+1];
    H[m][i  ] = c*a - s*b;
    H[m][i+1] = s*a + c*b;
  }
  
  /* Calculate mth rotation, store it in F */
  a = H[m][m];
  b = H[m][m+1];
  den = sqrt(a*a + b*b);
  F[2*m+0] = c =  a/den;
  F[2*m+1] = s = -b/den;
  
  /* Apply mth rotation to H */
  a = H[m][m]; 
  b = H[m][m+1];
  H[m][m  ] = c*a - s*b;
  H[m][m+1] = s*a + c*b;
  
  /* Apply mth rotation to g */
  a = g[m]; 
  b = g[m+1];
  g[m  ] = c*a - s*b;
  g[m+1] = s*a + c*b;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SolveUpper
static int
xf_SolveUpper(real **H, real *g, int m, real *y)
{
  /*
   
   PURPOSE:
   
   Used in GMRES.  Solves upper-triangular system H*y + g = 0;
   H is stored by columns.  Applied to first m  columns.
   
   INPUTS:
   
   H : Upper triangular matrix, stored by columns
   g : associated right-hand side vector
   m : size of system
   
   OUTPUTS:
   
   y : solution vector
   
   RETURNS:
   
   Error Code
   */
  int col, row;
  
  for (col=0; col<m; col++) y[col] = g[col];
  
  for (col=m-1; col>=0; col--){
    y[col] = -y[col]/H[col][col];
    for (row=col-1; row>=0; row--){
      y[row] += y[col]*H[col][row];
    }
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_SolveLinearSystem_GMRES
static int
xf_SolveLinearSystem_GMRES(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *R, 
                           enum xfe_Bool TransposeFlag, int nIterIn,
                           xf_SolverData *SolverData, xf_Vector *dU )
{
  /*
   PURPOSE:
   
   Solves linear system, R_U * dU + R = 0, using restarted GMRES:
   
   x0 = 0
   
   for iOuter = 0:nOuter-1,
   
   v0 = M^{-1}*(R_U * x0 + R)
   v0 = v0/||v0||  (beta = ||v0||)
   
   for j = 0:iInner-1,
   
   # Preconditioner (applied on left)
   w = M^{-1} R_U v_j
   
   # Modified Gram-Schmidt
   for i = 0:j,
   H[i,j] = (w, v_i)
	 w -= H(i,j) v_i
   end
   
   H[j+1,j] = ||w||
   v_{j+1} = w/||w||
   
   # QR factor H via Givens rotations on the fly (see [1])
   end 
   
   x0 += V*y, where V = {v_0, v_1, ...}, and y minimizes
   ||beta*e1 + H*y||  (see references)
   end
   
   Note, in the references, the system solved is usually A*x = b, and
   the residual is b-A*x, while in our case everything is on the left
   hand side.  The only difference that this introduces is that the
   minimized quantity is ||beta*e1 + H*y|| instead of ||beta*e1 -
   H*y||.  Hence, the SolveUpper function solves H*y + g = 0, instead
   of H*y - g = 0.
   
   References:
   
   [1] Saad and Schultz, SIAM J. Sci. Stat. Comput.  Vol 7, No. 3, 1986, 856-869
   [2] Sarkis and Szyld, Comp. Meth. in Appl. Mech. Eng.  Vol 196, 2007, 1612-1621
   
   INPUTS:
   
   All   : all structure
   R_U   : Jacobian matrix
   R     : "right-hand side", although actually appears on left (see above)
   TransposeFlag : if True, system solved will be R_U^T * dU + R = 0
   nIterIn : if >= 0, this number of outer iterations is used
   SolverData : solver data structure
   
   OUTPUTS:
   
   dU  : solution
   
   RETURN:
   
   Error Code
   
   */
  
  int ierr, i, nOuter, nInner;
  int iInner, iOuter;
  int *nH;
  enum xfe_Bool converged, FirstIteration, CheckHalt, CoarseCorrectionFlag;
  enum xfe_Verbosity Verbosity;
  real **H, *g, *y, *F, rnorm, rnorm0, normdU, normdU_prev, rcriterion=0.;
  real DecreaseFactor, fac, normdU0, LinResTol;
  xf_VectorSet *VS;
  xf_Vector *w, *v0, *x0, *Vtemp;
  enum xfe_PreconditionerType Preconditioner;
  
  // determine max number of outer iterations
  if (nIterIn < 0){
    ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nIterGMRESOuter", &nOuter));
    if (ierr != xf_OK) return ierr;
  }
  else nOuter = nIterIn;
  
  // determine max number of inner iterations
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "nIterGMRESInner", &nInner));
  if (ierr != xf_OK) return ierr;
  
  // do we need to check for halt during iterations?
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "CheckHaltLinear", &CheckHalt));
  if (ierr != xf_OK) return ierr;
  
  // determine verbosity
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  // determine residual decrease factor
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "MinLinResDecreaseFactor", 
                                     &DecreaseFactor));
  if (ierr != xf_OK) return ierr;
  
  // determine preconditioner
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Preconditioner", 
                                     xfe_PreconditionerName, (int ) xfe_PreconditionerLast, 
                                     (int *) &Preconditioner));
  if (ierr != xf_OK) return ierr;
  
  // determine if coarse correction is desired
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "CoarseCorrectionFlag", 
                                     &CoarseCorrectionFlag)); 
  if (ierr != xf_OK) return ierr;
  
  // precondition R_U
  ierr = xf_Error(xf_Jacobian_Precondition(All, R_U, Preconditioner, CoarseCorrectionFlag));
  if (ierr != xf_OK) return ierr;
  
  
  // find set of vectors for inner iterations (the v_j in the above description)
  ierr = xf_Error(xf_FindSimilarVectorSet(All, dU, nInner+1, "VSearch", 
                                          xfe_True, xfe_True, NULL, &VS));
  if (ierr != xf_OK) return ierr;
  
  // find a Vtemp vector
  ierr = xf_Error(xf_FindSimilarVector(All, dU, "GMRES_Vtemp", xfe_True,
                                       xfe_True, NULL, &Vtemp, NULL));
  if (ierr != xf_OK) return ierr;
  
  // find x0 vector
  ierr = xf_Error(xf_FindSimilarVector(All, dU, "x0", xfe_True, xfe_True, NULL, &x0, NULL));
  if (ierr != xf_OK) return ierr;
  
  // x0 = 0
  if (SolverData->ReusedU) {
    ierr = xf_Error(xf_SetVector(dU, xfe_Set, x0));
    if (ierr != xf_OK) return ierr;
  }
  else {
    ierr = xf_Error(xf_SetZeroVector(x0));
    if (ierr != xf_OK) return ierr;
  }
  
  // nH[i] = length of i'th column of H
  ierr = xf_Error(xf_Alloc( (void **) &nH, nInner, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nInner; i++) nH[i] = i+2;
  
  // allocate H via a variable-length alloc  # note, H is stored by COLUMNS
  ierr = xf_Error(xf_VAlloc2( (void ***) &H, nInner, nH, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // allocate g vector for storing rhs in reduced least-squares system
  ierr = xf_Error(xf_Alloc( (void **) &g, nInner+1, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // allocate y vector for storing solution in reduced least-squares system
  ierr = xf_Error(xf_Alloc( (void **) &y, nInner+1, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // allocate F vector for storing Givens rotations (c,s) at each iInner
  ierr = xf_Error(xf_Alloc( (void **) &F, 2*nInner, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // determine normdU_prev
  FirstIteration = xfe_False;
  if ((normdU_prev = SolverData->normdU_prev) < 0.) FirstIteration = xfe_True;
  
  // if LinResTolerance > 0, use this requested linear residual tolerance (relative)
  LinResTol = SolverData->LinResTol;
  
  /*----------------------------*/
  /* Begin outer iteration loop */
  /*----------------------------*/
  normdU0 = 0.; normdU = -1.0;
  converged = xfe_False;
  for (iOuter=0; iOuter<nOuter; iOuter++){
    
    // v0 points to VS[0]
    v0 = VS->Vector + 0;
    
    // v0 = M^{-1}*(R_U * x0 + R)
    ierr = xf_Error(xf_LinearStep(All, R_U, x0, R, Vtemp, xfe_Set, TransposeFlag, 
                                  SolverData, v0));
    if (ierr != xf_OK) return ierr;
    
    // rnorm = ||v0||
    ierr = xf_Error(xf_VectorNorm(v0, 2, &rnorm));
    if (ierr != xf_OK) return ierr;
    if (iOuter == 0){
      if ((rnorm0 = rnorm) == 0.){
        converged = xfe_True; // RHS was identically zero!
        break;
      }
    }
    
    //if (rnorm < MEPS) return xf_Error(xf_SINGULAR);
    
    // v0 *= 1/rnorm
    ierr = xf_Error(xf_VectorMult(v0, 1.0/rnorm));
    if (ierr != xf_OK) return ierr;
    
    // set g = [rnorm, 0, 0, ...]
    g[0] = rnorm;
    for (i=1; i<nInner+1; i++) g[i] = 0.0;
    
    /*----------------------------*/
    /* Begin inner iteration loop */
    /*----------------------------*/
    
    iInner = 0;
    while (iInner < nInner){
      
      // w points to VS[iInner + 1]
      w = VS->Vector + iInner+1;
      
      // w = M^{-1} R_U v_j
      ierr = xf_Error(xf_LinearStep(All, R_U, VS->Vector + iInner, NULL, Vtemp,
                                    xfe_Set, TransposeFlag, SolverData, w));
      if (ierr != xf_OK) return ierr;
      
      // Modified Gram-Schmidt
      for (i=0; i<=iInner; i++){
        
        // H[i,iInner] = (w, v_i)  # recall H is stored by COLUMNS
        ierr = xf_Error(xf_VectorDot(w, VS->Vector + i, H[iInner] + i ));
        if (ierr != xf_OK) return ierr;
        
        // w -= H(i,j) v_i
        ierr = xf_Error(xf_VectorMultSet(VS->Vector + i, H[iInner][i], xfe_Sub, w));
        if (ierr != xf_OK) return ierr;
        
      } // i
      
      // rnorm = ||w||
      ierr = xf_Error(xf_VectorNorm(w, 2, &rnorm));
      if (ierr != xf_OK) return ierr;
      
      //if (rnorm < MEPS) return xf_Error(xf_SINGULAR);
      
      // H[iInner+1,iInner] = rnorm
      H[iInner][iInner+1] = rnorm;
      
      // v_{j+1} = w/||w||
      ierr = xf_Error(xf_VectorMult(w, 1.0/rnorm));
      if (ierr != xf_OK) return ierr;
      
      // QR factor column that we just added to H
      ierr = xf_GivensQRStep(iInner, H, g, F);
      if(ierr != xf_OK) return ierr;
      
      rnorm = fabs(g[iInner+1]);
      
      iInner++;
      
      /* Convergence check */
      if (LinResTol > 0.0){
	rcriterion = max(LinResTol*rnorm0, MEPS);
      }
      else{
        // compute y = H^{-1}*g, where H has already been QR factored (in fact, H = R)
        ierr = xf_Error(xf_SolveUpper(H, g, iInner, y)); 
        if(ierr != xf_OK) return ierr;
        
        // compute approx |dU| so far
        for (i=0, normdU = normdU0*normdU0; i<iInner; i++) normdU += y[i]*y[i];
        normdU = sqrt(normdU);
        if (FirstIteration) normdU_prev = normdU;
        
        // stop if rnorm < rnorm0*(normdU/normdU_prev)^2 * DecreaseFactor
        fac = min(normdU/normdU_prev, 1.0);
        fac = min(fac*fac, DecreaseFactor);      
	rcriterion = max(rnorm0*fac, MEPS);
      }
      
      // converged if below criterion, so break
      if (rnorm < rcriterion){
	converged = xfe_True;
	break;
      }

    } // iInner
    
    if (Verbosity != xfe_VerbosityLow){
      xf_printf("   GMRES: iOuter = %d, iInner = %d, Linear norm = %.10E\n", 
                iOuter, iInner, rnorm);
    }
    
    // compute y = H^{-1}*g, where H has already been QR factored (in fact, H = R)
    ierr = xf_Error(xf_SolveUpper(H, g, iInner, y)); 
    if(ierr != xf_OK) return ierr;
    
    // set x0 += sum_i (VS[i]*y[i])
    for (i=0; i<iInner; i++){
      ierr = xf_Error(xf_VectorMultSet(VS->Vector + i, y[i], xfe_Add, x0));
      if (ierr != xf_OK) return ierr;
    }
    
    // store |dU| calculated so far
    normdU0 = normdU;
    
    // did we really converge?  Sometimes inner residual estimate is not accurate.
    if (converged){
      // v0 points to VS[0]
      v0 = VS->Vector + 0;
      // v0 = M^{-1}*(R_U * x0 + R)
      ierr = xf_Error(xf_LinearStep(All, R_U, x0, R, Vtemp, xfe_Set, TransposeFlag, 
				    SolverData, v0));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VectorNorm(v0, 2, &rnorm));  // this is the actual linear residual
      if (ierr != xf_OK) return ierr;
      converged = (rnorm < rcriterion);
    }

    // stop if converged
    if (converged) iOuter = nOuter;
    
    // check for halt
    if ((CheckHalt) && (xf_CheckUserHalt(NULL))) break;
    
  } // iOuter
  
  // set dU = x0
  ierr = xf_Error(xf_SetVector(x0, xfe_Set, dU));
  if (ierr != xf_OK) return ierr;
  
  SolverData->normdU_prev = normdU;
  
  // release memoroy
  xf_Release( (void  *) nH);
  xf_Release2((void **)  H);
  xf_Release( (void  *)  g);
  xf_Release( (void  *)  y);
  xf_Release( (void  *)  F);
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_SolveLinearSystem
int
xf_SolveLinearSystem(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *R, 
                     enum xfe_Bool TransposeFlag, int nIter, 
                     xf_SolverData *SolverData, xf_Vector *dU)
{
  int ierr; 
  enum xfe_LinearSolverType LinearSolver;
  
  // Linear solver type
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "LinearSolver", 
                                     xfe_LinearSolverName, (int ) xfe_LinearSolverLast, 
                                     (int *) &LinearSolver));
  if (ierr != xf_OK) return ierr;
  
  if (LinearSolver == xfe_GMRES){
    // call GMRES solver
    ierr = xf_SolveLinearSystem_GMRES(All, R_U, R, TransposeFlag, nIter, SolverData, dU);
    if (ierr != xf_OK) return ierr;
  }
  else if (LinearSolver == xfe_Iterative){
    // call Jacobi solver
    ierr = xf_SolveLinearSystem_Iterative(All, R_U, R, TransposeFlag, nIter, SolverData, dU);
    if (ierr != xf_OK) return ierr;
  }
  else if (LinearSolver == xfe_Debug){
    // call Debug solver
    ierr = xf_SolveLinearSystem_Debug(All, R_U, R, TransposeFlag, SolverData, dU);
    if (ierr != xf_OK) return ierr;
  }
  else
    return xf_Error(xf_NOT_SUPPORTED);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateLinearResidual
int
xf_CalculateLinearResidual(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X, 
                           xf_Vector *R, enum xfe_Bool TransposeFlag, 
                           xf_SolverData *SolverData, xf_Vector *LinR)
{
  int ierr; 
  
  if (R != NULL){
    // LinR = R
    ierr = xf_Error(xf_SetVector(R, xfe_Set, LinR));
    if (ierr != xf_OK) return ierr;
  }
  
  // LinR += R_U*X
  if (X != NULL){
    ierr = xf_Error(xf_Jacobian_Mult(All, R_U, X, xfe_Add, TransposeFlag, SolverData, LinR));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateLinearSolverData
static int
xf_CreateLinearSolverData(xf_LinearSolverData **pLinearSolverData)
{
  /*
   PURPOSE:
   
   Creates a structure for storing data required by a linear solver
   using call-back (currently just CG). 
   
   INPUTS:
   
   OUTPUTS:
   
   (*pLinearSolverData) : allocated LinearSolverData structure
   
   RETURN:
   
   Error Code
   
   */
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pLinearSolverData, 1,
                           sizeof(xf_LinearSolverData)));
  if (ierr != xf_OK) return ierr;
  
  (*pLinearSolverData)->iIter   = 0;
  (*pLinearSolverData)->PreconditionFlag = xfe_False;
  (*pLinearSolverData)->AfterPrecondition = xfe_False;
  (*pLinearSolverData)->rnorm2  = 0.0;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyLinearSolverData
static int
xf_DestroyLinearSolverData(xf_LinearSolverData *LinearSolverData)
{
  /*
   PURPOSE:  
   
   Destroys an LinearSolverData structure.
   
   INPUTS:  
   
   LinearSolverData : structure to be destroyed
   
   OUTPUTS: None
   
   RETURN:
   
   Error Code
   
   */
  xf_Release( (void *) LinearSolverData);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CGReadPoint
static int
xf_CGReadPoint(xf_All *All, xf_VectorSet *VS, xf_LinearSolverData *LinearSolverData)
{
  int ierr, i, nCG;
  int myRank, nProc;
  char line[xf_MAXLINELEN];
  FILE *fid;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_VectorSet *CGSet;
  
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  xf_printf("Restarting from existing CG.data and CG.txt\n");
  
  // Read in CG.txt
  if (myRank == 0){
    if ((fid = fopen("CG.txt", "r")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (sscanf(line, "%d", &LinearSolverData->iIter) != 1) return xf_Error(xf_FILE_READ_ERROR);
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (sscanf(line, "%d", (int *) &LinearSolverData->PreconditionFlag) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (sscanf(line, "%d", (int *) &LinearSolverData->AfterPrecondition) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (sscanf(line, "%lf", &LinearSolverData->rnorm2) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
    fclose(fid);
  }
  
  // broadcast info to all procs
  ierr = xf_Error(xf_MPI_Bcast((void *) &LinearSolverData->iIter, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &LinearSolverData->PreconditionFlag, 
                               sizeof(enum xfe_Bool), 0));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &LinearSolverData->AfterPrecondition, 
                               sizeof(enum xfe_Bool), 0));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Bcast((void *) &LinearSolverData->rnorm2, sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  
  
  // create dataset for CG
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;
  
  // read in CG.data
  ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, "CG.data", DataSet));
  if (ierr != xf_OK) return ierr;
  
  
  D = DataSet->Head;
  if (strncmp(D->Title, "CGSet", 5) != 0) return xf_Error(xf_INPUT_ERROR);
  CGSet = (xf_VectorSet *) D->Data;
  
  // Copy CG vectors to VS
  for (i=0; i<CGSet->nVector; i++){
    ierr = xf_Error(xf_SetVector(CGSet->Vector+i, xfe_Set, VS->Vector+i));
    if (ierr != xf_OK) return ierr;
  }
  
  // destroy DataSet
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CGSavePoint
static int
xf_CGSavePoint(xf_All *All, xf_VectorSet *VS, xf_LinearSolverData *LinearSolverData)
{
  int ierr, i, nCG;
  int myRank, nProc;
  char line[xf_MAXLINELEN];
  FILE *fid;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_VectorSet *CGSet;
  
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  nCG  = LinearSolverData->iIter;  // number of CG vectors for current run
  
  xf_printf("Saving CG.data and CG.txt: nIter = %d\n", nCG);
  
  
  // Write CG.txt
  if (myRank == 0){
    if ((fid = fopen("CG.txt", "w")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    fprintf(fid, "%d\n", nCG);
    fprintf(fid, "%d\n", LinearSolverData->PreconditionFlag);
    fprintf(fid, "%d\n", LinearSolverData->AfterPrecondition);
    fprintf(fid, "%.15E\n", LinearSolverData->rnorm2);
    fclose(fid);
  }
  
  // Write out CG.data, containing VS
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DataSetAdd(DataSet, "CGSet", xfe_VectorSet,
                                xfe_True, (void *) VS, &D));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, "CG.data"));
  if (ierr != xf_OK) return ierr;
  D->Data = NULL;
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_LinearIterCG
int
xf_LinearIterCG(xf_All *All, xf_VectorSet *VS, real tol, enum xfe_Verbosity Verbosity,
                enum xfe_Bool PreconditionFlag, enum xfe_Bool ZeroFlag,
                enum xfe_Bool SavePoint, enum xfe_Bool CGRestart,
                xf_Vector *R, xf_Vector *U, xf_Vector **pV, xf_Vector **pW, 
                xf_LinearSolverData **pLinearSolverData, 
                enum xfe_LinearStatusType *Status)
{
  int ierr, iIter;
  int ir = 0, ip = 1, iw = 2, iz = 3;
  enum xfe_Bool debug = xfe_False;
  static enum xfe_Bool Forming_r = xfe_False;
  real dp, val, rnorm, alpha, zdotr;
  xf_Vector *r, *p, *w, *z;
  xf_LinearSolverData *LinearSolverData;
  
  // print debug results when high verbosity is requested
  debug = (Verbosity == xfe_VerbosityHigh);
  
  if (VS->nVector < (3 + (int) PreconditionFlag)){
    xf_printf("Conjugate Gradient function requires at least 3 workspace functions.\n");
    xf_printf(" 4 if preconditioner is requested.\n");
    return xf_Error(xf_INPUT_ERROR);
  }
  
  // pull off work vectors
  r = VS->Vector + ir;
  p = VS->Vector + ip;
  w = VS->Vector + iw;
  z = ((PreconditionFlag) ? VS->Vector + iz : VS->Vector + ir);
  
  // on first call ...
  if (((*pLinearSolverData) == NULL) || (Forming_r)){
    
    if ((*pLinearSolverData) == NULL){
      // allocate LinearSolverData
      ierr = xf_Error(xf_CreateLinearSolverData(pLinearSolverData));
      if (ierr != xf_OK) return ierr;
      LinearSolverData = (*pLinearSolverData);
      
      /* If restarting from CG data, we pick up where we left off (which
       was right before we asked for a multiplication by A. */
      if ((CGRestart) && (All != NULL)){
        ierr = xf_Error(xf_CGReadPoint(All, VS, LinearSolverData));
        if (ierr != xf_OK) return ierr;
        (*Status) = xfe_LinearMultiply;
        (*pV) = VS->Vector + ip;
        (*pW) = VS->Vector + iw;
        return xf_OK;
      }
    }
    LinearSolverData = (*pLinearSolverData);
    
    if (!Forming_r){
      if (!ZeroFlag){
        // r = A*U
        Forming_r = xfe_True;
        (*Status) = xfe_LinearMultiply;
        (*pV) = U;
        (*pW) = VS->Vector + ir;
        return xf_OK;
      }
      else{
        // r = 0
        ierr = xf_Error(xf_SetZeroVector(r));
        if (ierr != xf_OK) return ierr;
      }
    }
    
    Forming_r = xfe_False;
    
    if (ZeroFlag){
      // U = 0
      ierr = xf_Error(xf_SetZeroVector(U));
      if (ierr != xf_OK) return ierr;
    }
    
    // r += R
    ierr = xf_Error(xf_SetVector(R, xfe_Add, r));
    if (ierr != xf_OK) return ierr;
    
    
    // rorm = sqrt(<r,r>)
    ierr = xf_Error(xf_VectorNorm(r, 2, &rnorm));
    if (ierr != xf_OK) return ierr;
    
    // U = 0 is solution if rnorm < tol
    if (rnorm < tol){
      if (debug) xf_printf("U0 satisfies tolerance in CG.  Exiting.\n");
      (*Status) = xfe_LinearConverged;
      ierr = xf_Error(xf_DestroyLinearSolverData(LinearSolverData));
      if (ierr != xf_OK) return ierr;
      return xf_OK;
    }
    
    
    // Request z = M^{-1}*r if need preconditioning
    LinearSolverData->PreconditionFlag = PreconditionFlag;
    if (PreconditionFlag){
      (*Status) = xfe_LinearPrecondition;
      (*pV) = r;
      (*pW) = z;
      return xf_OK;
    }
    
  }
  
  
  // increment iteration counter
  LinearSolverData = (*pLinearSolverData);
  iIter = LinearSolverData->iIter++;
  PreconditionFlag = LinearSolverData->PreconditionFlag;
  
  if (iIter == 0){
    // set p = z
    ierr = xf_Error(xf_SetVector(z, xfe_Set, p));
    if (ierr != xf_OK) return ierr;
    
    // rnorm2 = <z,r>, 
    ierr = xf_Error(xf_VectorDot(z, r, &zdotr));
    if (ierr != xf_OK) return ierr;
    LinearSolverData->rnorm2 = zdotr;
    
    // Request w = A*p
    (*Status) = xfe_LinearMultiply;
    (*pV) = VS->Vector + ip;
    (*pW) = VS->Vector + iw;
    return xf_OK;
  }
  
  
  if (!LinearSolverData->AfterPrecondition){
    
    // Back from computing  w = A*p
    
    // alpha = rnorm2 / <w,p>
    ierr = xf_Error(xf_VectorDot(w, p, &dp ));
    if (ierr != xf_OK) return ierr;
    alpha = LinearSolverData->rnorm2/dp;
    
    // U -= alpha*p
    ierr = xf_Error(xf_VectorMultSet(p, alpha, xfe_Sub, U));
    if (ierr != xf_OK) return ierr;
    
    // r -= alpha*w
    ierr = xf_Error(xf_VectorMultSet(w, alpha, xfe_Sub, r));
    if (ierr != xf_OK) return ierr;
    
    // Request z = M^{-1}*r if need preconditioning
    if (PreconditionFlag){
      LinearSolverData->AfterPrecondition = xfe_True;
      (*Status) = xfe_LinearPrecondition;
      (*pV) = r;
      (*pW) = z;
      return xf_OK;
    }
  }
  LinearSolverData->AfterPrecondition = xfe_False;
  
  
  // zdotr = <z,r>, 
  ierr = xf_Error(xf_VectorDot(z, r, &zdotr));
  if (ierr != xf_OK) return ierr;
  rnorm = sqrt(zdotr);
  
  /*   // rnorm = ||r|| */
  ierr = xf_Error(xf_VectorNorm(r, 2, &rnorm));
  if (ierr != xf_OK) return ierr;
  
  // rnorm < tol means convergence
  if (rnorm < tol){
    if (debug) xf_printf("CG converged in %d iterations.  Exiting.\n", iIter);
    (*Status) = xfe_LinearConverged;
    ierr = xf_Error(xf_DestroyLinearSolverData(LinearSolverData));
    if (ierr != xf_OK) return ierr;
    return xf_OK;
  }
  
  // p = z + rnorm*rnorm/rnorm2*p
  val = rnorm/sqrt(LinearSolverData->rnorm2); // to prevent underflow
  ierr = xf_Error(xf_VectorMult(p, val*val));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_SetVector(z, xfe_Add, p));
  if (ierr != xf_OK) return ierr;
  
  // rnorm2 = rnorm*rnrom
  LinearSolverData->rnorm2 = zdotr;
  
  
  /* Save away CG data for immediate restart */
  if ((SavePoint) && (All != NULL)){
    ierr = xf_Error(xf_CGSavePoint(All, VS, LinearSolverData));
    if (ierr != xf_OK) return ierr;
  }
  
  
  // Request w = A*p
  (*Status) = xfe_LinearMultiply;
  (*pV) = VS->Vector + ip;
  (*pW) = VS->Vector + iw;
  
  return xf_OK;
}


#if( UNIT_TEST==1 )
#include "xf_LinearSolver.test.in"
#endif
