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

#ifndef _xf_Dynamic_h
#define _xf_Dynamic_h 1

/*
  FILE:  xf_Dynamic.h

  This file contains headers for Dynamic library functions.

*/

#include "xf_DynamicStruct.h"


/******************************************************************/
//   FUNCTION Prototype: xf_DLError
extern char *
xf_DLError(void);
/*
PURPOSE:

  Returns string containing any error message regarding dynamic
  libraries.

INPUTS:

  None

OUTPUTS:

  None

RETURN:

  String containing error message
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DLOpen
extern xf_DLHandle 
xf_DLOpen(const char *LibName);
/*
PURPOSE:

  Loads a dynamic library

INPUTS:

  LibName : name of library to load

OUTPUTS:

  None

RETURN:

  Handle to library just opened
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DLSym
void *
xf_DLSym(xf_DLHandle LibHandle, const char *SymName);
/*
PURPOSE:

  Locates symbol SynName in library with handle LibHandle

INPUTS:

  LibHandle : handle of library to search
  SynName   : name of object for which to search

OUTPUTS:

  None

RETURN:

  pointer to object found

*/



/******************************************************************/
//   FUNCTION Prototype: xf_DLClose
extern void
xf_DLClose(xf_DLHandle LibHandle);
/*
PURPOSE:

  Closes a dynamic library

INPUTS:

  LibHandle: Handle of library to close

OUTPUTS:

  None

RETURN:

  None
*/


#endif // end ifndef _xf_Dynamic_h
