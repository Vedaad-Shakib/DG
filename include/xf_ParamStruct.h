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

#ifndef _xf_ParamStruct_h
#define _xf_ParamStruct_h 1

/*
  FILE:  xf_ParamStruct.h

  This file contains the xflow parameter structure

*/

#include "xf.h"


/* Key-Value structure */
typedef struct
{
  int nKey;
  /* Number of key-value pairs */
  
  char **Key;     
  /* Key[ikey] = string of key # ikey */
  
  char **Value;   
  /* Value[ikey] = string of value # ikey */
  
  int MaxStrLen;
  /* Maximum string length of keys + values */
  
  int DKey;
  /* Number of pairs to allocate at a time, to prevent frequent
     re-allocations*/

  int nKey0;
  /* Number of allocated key-value pairs.  Should always be >= nKey */

}
xf_KeyValue;


/* Parameter structure definition */
typedef struct
{ 

  xf_KeyValue KeyValue;
  /* Key-Value variables.  All non-equation-set-specific parameters
     are stored here; they are read and written. */
}
xf_Param;


#endif // end ifndef _xf_ParamStruct_h
