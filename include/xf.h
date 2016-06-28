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

#ifndef _xf_h
#define _xf_h 1

/*
  FILE:  xf.h

  This file contains commonly used "defines" and enumerated types

*/

/* include files */
#include <stddef.h>
#include "xf_Error.h"

typedef double real;
typedef long long xf_long; // for 32/64-bit compatibility
#define PI 3.141592653589793238
#define MEPS 1e-15 // machine precision
#define xf_MAXSTRLEN 150  // maximum string length, in characters
#define xf_MAXLINELEN 1024  // maximum line length, in characters
#define xf_MAXLONGLINELEN 1000  // maximum long line length, in characters
#define xf_ISAFETY 10000  // to check for infinite while loops
#define xfb_All 127 // all-on bitwise flag (for up to 7 flags)
#define xfb_None  0 // all-off bitwise flag
#define xf_INTERIORFACE -1 // designates interior face groups (do not change)
#define xf_NULLFACE -2 // designates null face groups
//#define RSTRCT restrict // c-99 standard for pointer passing optimization
#define RSTRCT  // c-99 standard for pointer passing optimization

#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#define sign(a) (((a) < 0.0) ? -1 : 1)
#define swap(a,b, t)  {t = a; a = b; b = t;}
#define cycle3(a,b,c, t)  {t = a; a = b; b = c; c = t;}
#define cycle4(a,b,c,d, t)  {t = a; a = b; b = c; c = d; d = t;}

/* True/False.  False always needs to correspond to 0. */
enum xfe_Bool { 
  xfe_False, 
  xfe_True, 
  xfe_BoolLast
};
static char *xfe_BoolName[xfe_BoolLast] = {
  "False", 
  "True"
};

/* Set/Negative/Add/Subtract */
enum xfe_AddType { 
  xfe_Set, 
  xfe_Neg,
  xfe_Add, 
  xfe_Sub,
  xfe_AddLast
};
static char *xfe_AddName[xfe_AddLast] = {
  "Set",
  "Neg",
  "Add",
  "Sub"
};


/* Relation type for comparing quantities/sets */
enum xfe_RelType { 
  xfe_RelEqual, 
  xfe_RelSubset, 
  xfe_RelSuperset,
  xfe_RelDisjoint,
  xfe_RelOther,
  xfe_RelLast
};
static char *xfe_RelName[xfe_RelLast] = {
  "Equal",
  "Subset",
  "Superset",
  "Disjoint",
  "Other"
};


/* Int/Real/Bool/Long in an enumerated type */
enum xfe_SizeType { 
  xfe_SizeInt, 
  xfe_SizeReal,
  xfe_SizeBool,
  xfe_SizeLong,
  xfe_SizeLast
};
static char *xfe_SizeName[xfe_SizeLast] = {
  "Int", 
  "Real",
  "Bool",
  "Long",
};


/* Verbosity level */
enum xfe_Verbosity { 
  xfe_VerbosityLow, 
  xfe_VerbosityMedium, 
  xfe_VerbosityHigh, 
  xfe_VerbosityLast
};
static char *xfe_VerbosityName[xfe_VerbosityLast] = {
  "Low", 
  "Medium", 
  "High"
};


/* Error Macro: used to report error occurrences */
#define xf_Error(X) (xf_ErrorReport( __FILE__, __LINE__, #X, (X)))
#define xf_PError(X, Y) (xf_PErrorReport( __FILE__, __LINE__, #X, (X), Y))

/* Error Codes. xf_OK must correspond to 0. */
enum xf_Errors {
  xf_OK,                                      // 0
  xf_MEMORY_ERROR,                            // 1
  xf_OUT_OF_BOUNDS,                           // 2
  xf_NOT_FOUND,                               // 3
  xf_FILE_READ_ERROR,			      // 4
  xf_FILE_WRITE_ERROR,			      // 5
  xf_SYSTEM_ERROR,                            // 6
  xf_PARALLEL_ERROR,			      // 7
  xf_INFINITE_WHILE_ERROR,		      // 8
  xf_MESH_ERROR,                              // 9
  xf_LINE_ERROR,                              // 10
  xf_SINGULAR,				      // 11
  xf_NO_UPDATE,				      // 12
  xf_CFL_BELOW_MINIMUM,			      // 13
  xf_NOT_CONVERGED,			      // 14
  xf_NON_PHYSICAL,			      // 15
  xf_CUT_CELL_ERROR,			      // 16
  xf_ADAPT_ERROR,			      // 17
  xf_STRING_ERROR,			      // 18
  xf_OVERWROTE,				      // 19
  xf_UNKNOWN_SHAPE,                           // 20
  xf_UNKNOWN_BASIS,			      // 21
  xf_NOT_SUPPORTED,			      // 22
  xf_NEGATIVE_JACOBIAN,			      // 23
  xf_SORT_ERROR,			      // 24
  xf_UNKNOWN_TYPE,			      // 25
  xf_EQNSET_ERROR,			      // 26
  xf_INPUT_ERROR,			      // 27
  xf_INCOMPATIBLE,			      // 28
  xf_NULL_OR_STALE_POINTER,		      // 29
  xf_UNSUPPORTED_QUAD_ORDER,                  // 30
  xf_BASIS_FCN_ERROR,			      // 31
  xf_CODE_LOGIC_ERROR,			      // 32
  xf_ASSERT_FAILED,			      // 33
  xf_PING_FAILED,			      // 34
  xf_MPI_ERROR,				      // 35
  xf_SOLVER_ERROR,			      // 36
  xf_DATA_ERROR,			      // 37
  xf_TREE_ERROR,			      // 38
  xf_DYNAMIC_LIBRARY_ERROR,		      // 39
  xf_BOUNDARY_CONDITION_ERROR,                // 40
  xf_FORCE_QUIT,			      // 41
  xf_INTERSECT_ERROR,			      // 42
  xf_LAPACK_ERROR,			      // 43
  xf_END_OF_FILE,			      // 44
  xf_TEST_FAILED,			      // 45
  xf_MULTIPLE_MATCHES,			      // 46
  xf_NEED_LAPACK,			      // 47
  xf_NEED_ADAPT,			      // 48
  xf_REWIND				      // 49
};


#endif // end ifndef _xf_h
