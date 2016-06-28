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

#ifndef _xf_MeshDistance_h
#define _xf_MeshDistance_h 1

/*
  FILE:  xf_MeshDistance.h

  This file contains the headers for Mesh distance functions

*/

/******************************************************************/
//   FUNCTION Prototype: xf_Dist2Face
extern int
xf_Dist2Face(int dim, const real *x0in, enum xfe_BasisType FBasis,
	     int Order, int nf, real *X, real *pdist, real *xproj);
/*
PURPOSE:

  Calculates minimum distance between x0 and a Lagrange-interpolated
  face.

  Note, the face geometry is interpolated using some Lagrange basis
  and order.  This function uses subdivision of the face into
  subfaces, and hence is approximate for curved faces.  On the plus
  side, this function is more robust than the Newton version
  (commented out in the "legacy" section below).
  

INPUTS:
 
  dim  : dimension of x0
  x0in : point from which to compute distance
  FBasis : Lagrange basis of face
  Order : Order with which face is interpolated
  nf : number of Lagrange nodes interpolating face
  X  : coords of face Lagrange points
  
OUTPUTS:

  (*pdist)  : calculated minimum distance, as described above
  xproj     : projection of x0in to face -- i.e. point on face
              that is closest to x0in.
 
RETURNS: Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_CalculateDistFcn
extern int 
xf_CalculateDistFcn(xf_All *All);
/*
PURPOSE:

  Calculates a distance function from wall boundaries to element
  interiors (a wall is identified by the equation set), and stores
  this distance in an interpolated element vector, "WallDistance".

INPUTS:

  All : all structure

OUTPUTS: 

RETURN:

  Error Code
*/




#endif // end ifndef _xf_MeshDistance_h
