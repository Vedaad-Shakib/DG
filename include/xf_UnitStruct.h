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

#ifndef _xf_UnitStruct_h
#define _xf_UnitStruct_h 1

/*
  FILE:  xf_UnitStruct.h

  This file contains structures used in unit testing.

*/

#include "xf.h"

/* Different mesh motions for unit tests */
enum xfe_UnitMotionType{
  xfe_UnitMotionNone,
  xfe_UnitMotionPlunge1, // plunge centered at 1,1
  xfe_UnitMotionPlunge2, // plunge centered at 1,0 (e.g. walls affected)
  xfe_UnitMotionLast
};


#endif // end ifndef _xf_UnitStruct_h
