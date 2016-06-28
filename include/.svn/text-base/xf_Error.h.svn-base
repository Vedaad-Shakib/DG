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

#ifndef _xf_Error_h
#define _xf_Error_h 1

/*
  FILE:  xf_Error.h

  This file contains the header for the error reporting functions.
  The error codes are defined in xf.h.

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ErrorReport
extern int xf_ErrorReport( char *file, int line, char *call, int ierr);
/*
PURPOSE:
  Function called by the macro xf_Error to report any error instances

INPUTS:
  file:  file name in which the error occurred
  line:  line number at which the error (call) occurred
  call:  string containing the name of the call that produced the error
  ierr:  error code

OUTPUTS: 
  None; prints an error message if neccessary

RETURN:
  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_PErrorReport
extern int xf_PErrorReport( char *file, int line, char *call, int *pierr, int root);
/*
PURPOSE:

  Function called by the macro xf_PError to report any error instances
  in parallel.  The parallel aspect is that (*pierr) is first
  broadcast to all processors from root before being checked.  A
  pointer is required because the error is modified during the
  broadcast so that all procs have the same error as on root.

INPUTS:
  file :  file name in which the error occurred
  line :  line number at which the error (call) occurred
  call :  string containing the name of the call that produced the error
  pierr:  error code
  root :  root processor number

OUTPUTS: 
  (*pierr): modified on other procs to have root's ierr

RETURN:
  Error Code
*/




#endif // end ifndef _xf_Error_h
