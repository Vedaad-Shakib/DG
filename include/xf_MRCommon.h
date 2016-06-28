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

#ifndef _xf_MRCommon_h
#define _xf_MRCommon_h 1

/*
  FILE:  xf_MRCommon.h

  This file contains function prototypes for functions in xf_MRCommon.c

*/

/******************************************************************/
//   FUNCTION Prototype: xf_AllocIPoint
extern int 
xf_AllocIPoint(xf_IPointType *IP, int M, int dim);


/******************************************************************/
//   FUNCTION Prototype: xf_CreateReducedModel
extern int 
xf_CreateReducedModel(xf_ReducedModel **pRM);

/******************************************************************/
//   FUNCTION Prototype: xf_AllocReducedModelLinear
extern int 
xf_AllocReducedModelLinear(xf_ReducedModel *RM, int N, int nOutput);


/******************************************************************/
//   FUNCTION Prototype: xf_AllocReducedModelNonLinear
extern int 
xf_AllocReducedModelNonLinear(xf_ReducedModel *RM, int M);


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyReducedModel
extern int 
xf_DestroyReducedModel(xf_ReducedModel *RM, enum xfe_Bool DestroyEqnSet);


/******************************************************************/
//   FUNCTION Prototype: xf_WriteReducedModelBinary
extern int 
xf_WriteReducedModelBinary( xf_ReducedModel *RM, FILE *fid);

/******************************************************************/
//   FUNCTION Prototype: xf_ReadReducedModelBinary
extern int 
xf_ReadReducedModelBinary( FILE *fid, xf_ReducedModel *RM);

/******************************************************************/
//   FUNCTION Prototype: xf_WriteReducedModelAscii
extern int 
xf_WriteReducedModelAscii( xf_ReducedModel *RM, FILE *fid);

#endif // end ifndef _xf_MRCommon_h
