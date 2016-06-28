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

#ifndef _xf_GeomAnalytical_h
#define _xf_GeomAnalytical_h 1

/*
  FILE:  xf_GeomAnalytical.h

  This file contains the headers for functions dealing with analytical
  geometry components.

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectToGeomComp_Analytical
extern int 
xf_ProjectToGeomComp_Analytical( xf_GeomCompAnalytical *GCA, int dim, real *x);
/*
PURPOSE:

  Projects x to analytical geometry component GCA.

INPUTS:

  GCA      : Analytical geometry component
  dim      : spatial dimension of x
  x        : on input, original location of point

OUTPUTS: 

  x        : on output, new (projected) location of point

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_PointsOnGeom_Analytical
extern int 
xf_PointsOnGeom_Analytical( xf_GeomCompAnalytical *GCA, int dim, 
			    int np, real dl, enum xfe_GeomSpacingType Spacing, 
			    int *pnout, real **px);
/*
PURPOSE:

  Returns points on analytical Geometry component iComp, according to
  a specified spacing.

INPUTS:

  GCA      : pointer to analytical geometry structure
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



#endif // end ifndef _xf_GeomAnalytical_h
