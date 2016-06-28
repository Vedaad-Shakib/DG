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

#ifndef _xf_QuadRule_h
#define _xf_QuadRule_h 1

/*
  FILE:  xf_QuadRule.h

  This file contains the headers for functions that return quadrature points

*/

/******************************************************************/
//   FUNCTION Prototype: xf_QuadLine
extern int 
xf_QuadLine( int Order, int *pnq, real **pxq, real **pwq);
/*
PURPOSE:

  Returns quadrature points and weights for integrating up to Order on
  a reference unitsegment [0,1].  (*pxq) and (*pwq) are reallocated.

INPUTS:

  Order: desired order of integration

OUTPUTS: 

  (*pnq) : number of quad points
  (*pxq) : ref-space coords of quad points: x0,x1,..
  (*pwq) : corresponding weights (optional; can pass in as NULL)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_QuadTriangle
extern int 
xf_QuadTriangle( int Order, int *pnq, real **pxq, real **pwq);
/*
PURPOSE:

  Returns quadrature points and weights for integrating up to Order on
  a reference unit right triangle [0,1].  (*pxq) and (*pwq) are
  reallocated.

INPUTS:

  Order: desired order of integration

OUTPUTS: 

  (*pnq) : number of quad points
  (*pxq) : ref-space coords of quad points: x0,y0,x1,y1,..
  (*pwq) : corresponding weights

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_QuadTetrahedron
extern int 
xf_QuadTetrahedron( int Order, int *pnq, real **pxq, real **pwq);
/*
PURPOSE:

  Returns quadrature points and weights for integrating up to Order on
  a reference unit right tetrahedron [0,1].  (*pxq) and (*pwq) are
  reallocated.

INPUTS:

  Order: desired order of integration

OUTPUTS: 

  (*pnq) : number of quad points
  (*pxq) : ref-space coords of quad points: x0,y0,z0,x1,y1,z1,..
  (*pwq) : corresponding weights

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_QuadQuadrilateral
extern int 
xf_QuadQuadrilateral( int Order, int *pnq, real **pxq, real **pwq);
/*
PURPOSE:

  Returns quadrature points and weights for integrating up to Order on
  a reference unit square [0,1].  (*pxq) and (*pwq) are reallocated.

INPUTS:

  Order: desired order of integration

OUTPUTS: 

  (*pnq) : number of quad points
  (*pxq) : ref-space coords of quad points: x0,y0,x1,y1,..
  (*pwq) : corresponding weights

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_QuadHexahedron
extern int 
xf_QuadHexahedron( int Order, int *pnq, real **pxq, real **pwq);
/*
PURPOSE:

  Returns quadrature points and weights for integrating up to Order on
  a reference unit cube [0,1].  (*pxq) and (*pwq) are reallocated.

INPUTS:

  Order: desired order of integration

OUTPUTS: 

  (*pnq) : number of quad points
  (*pxq) : ref-space coords of quad points: x0,y0,z0,x1,y1,z1,..
  (*pwq) : corresponding weights

RETURN:

  Error Code
*/






#endif // end ifndef _xf_QuadRule_h
