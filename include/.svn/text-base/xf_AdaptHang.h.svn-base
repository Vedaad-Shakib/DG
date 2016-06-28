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

#ifndef _xf_AdaptHang_h
#define _xf_AdaptHang_h 1

/*
  FILE:  xf_AdaptHang.h

  This file contains the headers for hanging-node adaptation mechanics.

*/



/******************************************************************/
//   FUNCTION Prototype: xf_Ref2nElem
extern int 
xf_Ref2nElem(enum xfe_ShapeType Shape, int ref, int *pnelem);
/*
PURPOSE:

  Returns number of sub-elements resulting from a refinement of Shape

INPUTS:

  Shape : shape of element in question
  ref  : refinement indicator for this element

OUTPUTS: 

  (*pnelem) : number of subelements

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_AdaptHang
extern int 
xf_AdaptHang(xf_All *All, xf_Vector *RefIndicator);

/*
PURPOSE:

  Adapts All->Mesh according to the desired refinement specified in
  the elemental vector RefIndicator.  Additional elements may be
  flagged and adapted to maintain the required one-level depth of
  non-conformity.  Regarding All->Data: vectors are reallocated and
  projected if interpolated, or destroyed otherwise (e.g. if a vector
  stored the areas, safest to destroy it and recreate for the new mesh
  when required); matrices, including the Jacobian, are destroyed
  (will be recreated when required on the new mesh).

INPUTS:

  All : All structure
  RefIndicator : refinement indicator (integer vector over all elements)

OUTPUTS: 

  None: All is modified.  In particular, structures are rallocated
        and some data may be deleted

RETURN:

  Error Code
*/




#endif // end ifndef _xf_AdaptHang_h
