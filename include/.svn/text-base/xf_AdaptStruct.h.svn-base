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

#ifndef _xf_AdaptStruct_h
#define _xf_AdaptStruct_h 1

/*
  FILE:  xf_AdaptStruct.h

  This file contains the xflow Adaptation structures

*/


#include "xf.h"



/* The adaptation performed is driven by this type */
enum xfe_AdaptOnType{
  xfe_AdaptOnOutput,
  xfe_AdaptOnResidual,
  xfe_AdaptOnPenalty,
  xfe_AdaptOnScalar,
  xfe_AdaptOnUniform,
  xfe_AdaptOnInterpol,
  xfe_AdaptOnTest,
  xfe_AdaptOnLast
};

/* corresponding names */
static char *xfe_AdaptOnName[xfe_AdaptOnLast] = {
  "Output",
  "Residual",
  "Penalty",
  "Scalar",
  "Uniform",
  "Interpol",
  "Test"
};



/* These are the available adaptation mechanics (means of changing the mesh) */
enum xfe_AdaptMechType{
  xfe_AdaptMechHangNode,
  xfe_AdaptMechOrderRef,
  xfe_AdaptMechLast
};

/* corresponding names */
static char *xfe_AdaptMechName[xfe_AdaptMechLast] = {
  "HangNode",
  "OrderRef"
};

/* These are the available cost metrics for optimization-based refinement */
enum xfe_AdaptCostType{
  xfe_AdaptCostDeltaDOF,
  xfe_AdaptCostDOF,
  xfe_AdaptCostDeltaNonZeros,
  xfe_AdaptCostNonZeros,
  xfe_AdaptCostLast
};

/* corresponding names */
static char *xfe_AdaptCostName[xfe_AdaptCostLast] = {
  "DeltaDOF",
  "DOF",
  "DeltaNonZeros",
  "NonZeros"
};

/*---------------------------------*/
/* Hanging-node refinement options */
/*---------------------------------*/

/* Types of refinements defined for segments */
enum xfe_SegRefType{
  xfe_SegRefNone,
  xfe_SegRefUniform,
  xfe_SegRefLast
};

/* corresponding names */
static char *xfe_SegRefName[xfe_SegRefLast] = {
  "None",
  "Uniform"
};

/* Types of refinements defined for quads */
enum xfe_QuadRefType{
  xfe_QuadRefNone,
  xfe_QuadRefUniform,
  xfe_QuadRefHoriz,
  xfe_QuadRefVert,
  xfe_QuadRefLast
};

/* corresponding names */
static char *xfe_QuadRefName[xfe_QuadRefLast] = {
  "None",
  "Uniform",
  "Horiz",
  "Vert",
};


/* Types of refinements defined for hexes */
enum xfe_HexRefType{
  xfe_HexRefNone,
  xfe_HexRefUniform,
  xfe_HexRefSliceX,    // cut plane normal is in x direction
  xfe_HexRefSliceY,
  xfe_HexRefSliceZ,
  xfe_HexRefSliceXY,   // two cut planes, with normals in x and y
  xfe_HexRefSliceXZ,
  xfe_HexRefSliceYZ,
  xfe_HexRefLast
};

/* corresponding names */
static char *xfe_HexRefName[xfe_HexRefLast] = {
  "None",
  "Uniform",
  "SliceX",
  "SliceY",
  "SliceZ",
  "SliceXY",
  "SliceXZ",
  "SliceYZ",
};





#endif // end ifndef _xf_AdaptStruct_h
