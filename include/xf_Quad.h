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

#ifndef _xf_Quad_h
#define _xf_Quad_h 1

/*
  FILE:  xf_Quad.h

  This file contains the headers for quadrature-specific functions.

*/

#include "xf_QuadStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_CreateQuadData
extern int 
xf_CreateQuadData( xf_QuadData **pQuadData);
/*
PURPOSE:

  Allocates memory for a QuadData structure

INPUTS:

  pQuadData: address of pointer to a QuadData structure

OUTPUTS: 

  None, (*pQuadData) is allocated and initialized

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyQuadData
extern int 
xf_DestroyQuadData( xf_QuadData *QuadData);
/*
PURPOSE:

  Destroys QuadData structure

INPUTS:

  QuadData: pointer to a QuadData structure

OUTPUTS: 

  None, QuadData is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyGenericQuadData
extern int 
xf_DestroyGenericQuadData( xf_QuadData *QuadData);
/*
PURPOSE:

  Destroys QuadData structure only if Type == xfe_QuadDataTypeGeneric

INPUTS:

  QuadData: pointer to a QuadData structure

OUTPUTS: 

  None, QuadData is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetQuadOrderElem
extern int 
xf_GetQuadOrderElem( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet, 
		     int egrp, int Order, int *QuadOrder);
/*
PURPOSE:

  Returns the appropriate quadrature order, after calling an
  eqnset-specific function and accounting for the geometry order.

INPUTS:

  Mesh : mesh structure
  EqnSet : equation set structure; Null to not perform eqnset call.
  egrp : element group
  Order : state interpolation order passed into the eqnset function.
          Order == -1 indicates maximum quadrature rule is desired

OUTPUTS: 

  (*QuadOrder): contains the quadrature order

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetQuadOrderIFace
extern int 
xf_GetQuadOrderIFace( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet, 
		      xf_IFace IFace, int Order, int *QuadOrder);
/*
PURPOSE:

  Returns the appropriate quadrature order, after calling an
  eqnset-specific function and accounting for the geometry order.

INPUTS:

  Mesh : mesh structure
  IFace : interior face to consider
  EqnSet : equation set structure; Null to not perform eqnset call.
  Order : state interpolation order passed into the eqnset function.
          Order == -1 indicates maximum quadrature rule is desired

OUTPUTS: 

  (*QuadOrder): contains the quadrature order

RETURN:

  Error Code
*/

/******************************************************************/
int
xf_GetQuadOrderAcrossFace( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet,
                            int egrpL, int egrpR, int Order, int *QuadOrder);
/*
 * One simple modified version of the above routine.
 */

/******************************************************************/
//   FUNCTION Prototype: xf_GetQuadOrderGeneralFace
extern int 
xf_GetQuadOrderGeneralFace( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet, 
			    int egrp, int Order, int *QuadOrder);
/*
PURPOSE:

  Returns the appropriate quadrature order for a general face of
  element group egrp, after calling an eqnset-specific function and
  accounting for the geometry order.

INPUTS:

  Mesh : mesh structure
  egrp : element group number
  EqnSet : equation set structure; Null to not perform eqnset call.
  Order : state interpolation order passed into the eqnset function.
          Order == -1 indicates maximum quadrature rule is desired

OUTPUTS: 

  (*QuadOrder): contains the quadrature order

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetQuadOrderBFace
extern int 
xf_GetQuadOrderBFace( const xf_Mesh *Mesh, const xf_EqnSet *EqnSet, 
		      xf_BFace BFace, int Order, int *QuadOrder);
/*
PURPOSE:

  Returns the appropriate quadrature order, after calling an
  eqnset-specific function and accounting for the geometry order.

INPUTS:

  Mesh : mesh structure
  BFace : boundary face to consider
  EqnSet : equation set structure; Null to not perform eqnset call.
  Order : state interpolation order passed into the eqnset function.
          Order == -1 indicates maximum quadrature rule is desired

OUTPUTS: 

  (*QuadOrder): contains the quadrature order

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_GetQuadPoints
extern int 
xf_GetQuadPoints( xf_QuadData *QuadData);
/*
PURPOSE:

  Computes quadrature points for a reference shape, adequate for
  integrating up to QuadData->Order.  The points and weights in
  QuadData are reallocated.

INPUTS:

  QuadData: pointer to an allocated QuadData structure, which has
            valid QuadData->Order.  The points and weights need
	    not be allocated.
  

OUTPUTS: 

  QuadData: contains quadrature points and weights in ref space

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_QuadElem
extern int 
xf_QuadElem( xf_Mesh *Mesh, int egrp, int elem, int Order, 
	     xf_QuadData **pQuadData, enum xfe_Bool *changed);
/*
PURPOSE:

  Returns a QuadData structure appropriate to the input element and
  the input Order.  Does not recalculate points if (*pQuadData)
  already contains suitable points.

INPUTS:

  Mesh : Mesh structure
  egrp, elem : input element info
  Order : requested integration order (-1 for maximum available rule)
  pQuadData: address of pointer to a QuadData structure
  

OUTPUTS: 

  (*pQuadData): contains pointer to suitable QuadData structure
  (*changed) : True if the points/weights in (*pQuadData) changed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_QuadElemAtLeast
extern int 
xf_QuadElemAtLeast( xf_Mesh *Mesh, int egrp, int elem, int Order, int nmin,
		    xf_QuadData **pQuadData, enum xfe_Bool *changed);
/*
PURPOSE:

  Wrapper for QuadElem.  Returns at least nn points.  Does this by
  upping the Order until have enough points.  Returns an error if
  maxes out quad rule without meeting nq > nn criterion.

INPUTS:

  Mesh : Mesh structure
  egrp, elem : input element info
  Order : requested integration order (-1 for maximum available rule)
  nmin : minimum # of quad points desired
  pQuadData: address of pointer to a QuadData structure
  

OUTPUTS: 

  (*pQuadData): contains pointer to suitable QuadData structure
  (*changed) : True if the points/weights in (*pQuadData) changed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_QuadElemJustAbove
extern int 
xf_QuadElemJustAbove( xf_Mesh *Mesh, int egrp, int elem, int nmin, 
		      xf_QuadData **pQuadData, enum xfe_Bool *changed);
/*
PURPOSE:

  Wrapper for QuadElem.  Returns at least nn points, but just barely.
  Does this by upping the Order (starting from 1) until have enough
  points.  Returns an error if maxes out quad rule without meeting nq
  > nn criterion.

INPUTS:

  Mesh : Mesh structure
  egrp, elem : input element info
  nmin : minimum # of quad points desired
  pQuadData: address of pointer to a QuadData structure
  

OUTPUTS: 

  (*pQuadData): contains pointer to suitable QuadData structure
  (*changed) : True if the points/weights in (*pQuadData) changed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_QuadFace
extern int 
xf_QuadFace( xf_Mesh *Mesh, int egrp, int elem, int face, int Order, 
	     xf_QuadData **pQuadData, enum xfe_Bool *changed);
/*
PURPOSE:

  Returns a QuadData structure appropriate to the input egrp,elem,face
  and the input Order.  Does not recalculate points if (*pQuadData)
  already contains suitable points.

INPUTS:

  Mesh : Mesh structure
  egrp,elem,face : input face info
  Order : requested integration order (-1 for maximum available rule)
  pQuadData: address of pointer to a QuadData structure
  

OUTPUTS: 

  (*pQuadData): contains pointer to suitable QuadData structure
  (*changed) : True if the points/weights in (*pQuadData) changed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_QuadIFace
extern int 
xf_QuadIFace( xf_Mesh *Mesh, xf_IFace IFace, int Order, 
	     xf_QuadData **pQuadData, enum xfe_Bool *changed);
/*
PURPOSE:

  Returns a QuadData structure appropriate to the input face and the
  input Order.  Does not recalculate points if (*pQuadData) already
  contains suitable points.  The points are defined in the face
  reference space.

INPUTS:

  Mesh : Mesh structure
  egrp, elem : input element info
  Order : requested integration order (-1 for maximum available rule)
  pQuadData: address of pointer to a QuadData structure
  

OUTPUTS: 

  (*pQuadData): contains pointer to suitable QuadData structure
  (*changed) : True if the points/weights in (*pQuadData) changed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteQuadDataBinary
extern int 
xf_WriteQuadDataBinary( xf_QuadData *QuadData, FILE *fid);
/*
PURPOSE:

  Writes QuadData to a binary file

INPUTS:

  QuadData : QuadData structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadQuadDataBinary
extern int 
xf_ReadQuadDataBinary( FILE *fid, xf_QuadData *QuadData);
/*
PURPOSE:

  Reads QuadData from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  QuadData : QuadData structure to read

RETURN:

  Error Code
*/


#endif // end ifndef _xf_Quad_h
