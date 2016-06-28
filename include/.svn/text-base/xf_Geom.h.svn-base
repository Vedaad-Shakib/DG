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

#ifndef _xf_Geom_h
#define _xf_Geom_h 1

/*
  FILE:  xf_Geom.h

  This file contains the headers for functions dealing with the Geom structure.

*/

#include "xf_GeomStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_CreateGeom
extern int 
xf_CreateGeom( xf_Geom **Geom);
/*
PURPOSE:

  Creates an Geom structure and all of its children.  Memory is
  allocated and initial (zero) values are set.

INPUTS:

  Geom : pointer to Geom

OUTPUTS: 

  None: Geom is allocated

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyGeom
extern int 
xf_DestroyGeom( xf_Geom *Geom);
/*
PURPOSE:

  Destroys a Geom structure and all of its children.  Memory is
  released

INPUTS:

  Geom : pointer to Geom

OUTPUTS: 

  None: Geom is destroyed

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ProjectToGeom
extern int 
xf_ProjectToGeom( xf_Geom *Geom, int iComp, char *BFGTitle, 
		  int dim, int np, real *x);
/*
PURPOSE:

  Projects np points provided in x to geometry.  Overwrites data in x.

INPUTS:

  Geom     : Geometry structure
  iComp    : component index (optional)
  BFGTitle : boundary face group title (optional, need if iComp < 0)
  dim      : spatial dimension of provided coordinates
  np       : number of points
  x        : on input, original locations of points

OUTPUTS: 

  x        : on output, new (projected) locations of points

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_PointsOnGeom
extern int 
xf_PointsOnGeom( xf_Geom *Geom, int iComp, int dim, int np, real dl, 
		 enum xfe_GeomSpacingType Spacing, int *pnout, real **px);
/*
PURPOSE:

  Returns points on Geometry component iComp, according to a specified
  spacing.

INPUTS:

  Geom     : pointer to Geom structure
  iComp    : component index
  dim      : dimension (only dim=2 supported for now)
  np       : number of points requested
  dl       : average spacing requested
  Spacing  : type of spacing requested
  
OUTPUTS: 

  (*pnout) : actual number of points provided
  (*px)    : actual points, unrolled, allocated in this function

RETURN:

  Error Code
*/




#endif // end ifndef _xf_Geom_h
