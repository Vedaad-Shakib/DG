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

#ifndef _xf_GeomSpline_h
#define _xf_GeomSpline_h 1

/*
  FILE:  xf_GeomSpline.h

  This file contains the headers for functions dealing with spline
  geometry components.

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectToGeomComp_Spline
extern int 
xf_ProjectToGeomComp_Spline( xf_GeomCompSpline *GCS, int dim, real *x);
/*
PURPOSE:

  Projects x to spline geometry component GCS.

INPUTS:

  GCS      : Spline geometry component
  dim      : spatial dimension of x
  x        : on input, original location of point

OUTPUTS: 

  x        : on output, new (projected) location of point

RETURN:

  Error Code
*/




/******************************************************************/
//   FUNCTION Prototype: xf_PointsOnGeom_Spline
extern int 
xf_PointsOnGeom_Spline( xf_GeomCompSpline *GCS, int dim, 
			int np, real dl, enum xfe_GeomSpacingType Spacing, 
			int *pnout, real **px);
/*
PURPOSE:

  Returns points on spline Geometry component iComp, according to
  a specified spacing.

INPUTS:

  GCS      : pointer to spline geometry structure
  dim      : dimension
  np       : number of points requested
  dl       : average spacing requested
  Spacing  : type of spacing requested
  
OUTPUTS: 

  (*pnout) : actual number of points provided
  (*px)    : actual points, unrolled, allocated in this function

RETURN:

  Error Code
*/


#endif // end ifndef _xf_GeomSpline_h
