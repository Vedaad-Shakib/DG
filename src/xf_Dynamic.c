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

/*
  FILE:  xf_Dynamic.c

  This file contains functions for working with dynamic libraries.

*/


#include "xf_AllStruct.h"
#include "xf_DynamicStruct.h"
#include "xf_IO.h"
#include <dlfcn.h>


/******************************************************************/
//   FUNCTION Definition: xf_DLError
char *
xf_DLError(void)
{
  return dlerror();
}

/******************************************************************/
//   FUNCTION Definition: xf_DLOpen
xf_DLHandle 
xf_DLOpen(const char *LibName)
{
  return dlopen(LibName, RTLD_LAZY);
}

/******************************************************************/
//   FUNCTION Definition: xf_DLSym
void *
xf_DLSym(xf_DLHandle LibHandle, const char *SymName)
{
  return dlsym(LibHandle, SymName);
}

/******************************************************************/
//   FUNCTION Definition: xf_DLClose
void
xf_DLClose(xf_DLHandle LibHandle)
{
  int ierr;
  if ((LibHandle != NULL) && (dlclose(LibHandle) != 0))
    xf_printf("Warning, could not close dynamic library. Error = %s\n", dlerror());
}
