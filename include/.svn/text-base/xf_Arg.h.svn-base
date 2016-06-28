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

#ifndef _xf_Arg_h
#define _xf_Arg_h 1

/*
  FILE:  xf_Arg.h

  This file contains headers for functions in xf_Arg.c

*/

#include <stdio.h>


/******************************************************************/
//   FUNCTION Prototype: xf_ParseArg
extern int 
xf_ParseArg(char **Keys, int argc, char *argv[], xf_KeyValue *KeyValue);
/*
PURPOSE:

  Parses argument list from command line into KeyValue pairs acording
  to recognized Keys list passed in.  The argument list is expected to
  be of the form:

  ProgramName -key0 value0 -key1 value1 ...

  Note the leading '-' in front of every key.

INPUTS:

  Keys : string array of recognized Keys, Default values, and Documentation.
         
         e.g.

	 Keys = {"key0", "value0", "This is key0",
  	         "key1", "value1", "This is key1",
		 ...
		 "\0"};

         Note final null terminating string, "\0".

  argc : # arguments (from command line)
  argv : string array of arguments (from command line)

OUTPUTS: 

  KeyValue : key-value list containing any keys set via the arguments.
  Keys not set have their values = "NULL".

RETURN:

  Error Code

*/


#endif // end ifndef _xf_Arg_h
