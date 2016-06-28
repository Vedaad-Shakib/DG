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

#ifndef _xf_LinearSolverStruct_h
#define _xf_LinearSolverStruct_h 1

/*
  FILE:  xf_LinearSolverStruct.h

  This file contains the structures for the linear solver

*/

#include "xf.h"
#include "xf_LineStruct.h"


/* Linear Solver type; others could be Direct, ... */
enum xfe_LinearSolverType{
  xfe_LinearSolverNone,
  xfe_GMRES,
  xfe_Iterative,
  xfe_Debug,
  xfe_LinearSolverLast
};

/* corresponding names */
static char *xfe_LinearSolverName[xfe_LinearSolverLast] = {
  "None",
  "GMRES",
  "Iterative",
  "Debug",
};

/* Preconditioner type */
enum xfe_PreconditionerType{
  xfe_PreconditionerNone,
  xfe_PreconditionerPMG,
  xfe_PreconditionerBlockJacobi,
  xfe_PreconditionerBlockJacobiLean,
  xfe_PreconditionerLineJacobi,
  xfe_PreconditionerLineGS,
  xfe_PreconditionerILU0,
  xfe_PreconditionerLast
};

/* corresponding names */
static char *xfe_PreconditionerName[xfe_PreconditionerLast] = {
  "None",
  "PMG",
  "BlockJacobi",
  "BlockJacobiLean",
  "LineJacobi",
  "LineGS",
  "ILU0",
};


/* Status of linear solver  */
enum xfe_LinearStatusType{
  xfe_LinearMultiply,
  xfe_LinearPrecondition,
  xfe_LinearConverged,
  xfe_LinearStatusLast,
};

/* corresponding names */
static char *xfe_LinearStatusName[xfe_LinearStatusLast] = {
  "Multiply",
  "Precondition",
  "Converged",
};


/* LinearSolverData structure (not read or written) */
typedef struct
{
  int  iIter;   // current iteration number
  enum xfe_Bool PreconditionFlag; // true if preconditioner is used
  enum xfe_Bool AfterPrecondition; // true after M^{-1} was just applied

  real rnorm2;  // <r,r> (see documentation in CG function header)

}
xf_LinearSolverData;


/* ILU ordering type */ 
enum xfe_ILUOrderingType{
  xfe_ILUOrderingNone,
  xfe_ILUOrderingMDF,       // minimum discarded fill
  xfe_ILUOrderingLex,       // Lexicographical ordering
  xfe_ILUOrderingLast
};

/* corresponding names */
static char *xfe_ILUOrderingName[xfe_ILUOrderingLast] = {
  "None",
  "MDF",
  "Lex",
};


/* ILU structure */
typedef struct
{
  int nelemtot;     // total number of elements
  int **Index2Elem; // Index2Elem[i][0,1] = (egrp,elem) at ordered index i
  int **Elem2Index; // Elem2Index[egrp][elem] = ordered index of element
  real ***C;        // matrix for MDF weight calculation
  int rmax;         // max r/StateRank over all elements
  int nfacemax;     // max number of faces over all elements
}
xf_LinearSolverILUData;




#endif // end ifndef _xf_LinearSolverStruct_h
