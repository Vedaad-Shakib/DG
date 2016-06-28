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

#ifndef _xf_Log_h
#define _xf_Log_h 1

/*
  FILE:  xf_Log.h

  This file contains the headers for functions in xf_Log.c

*/

/******************************************************************/
//   FUNCTION Prototype: xf_WriteLogHeader
extern int 
xf_WriteLogHeader( const xf_All *All, const char *PreHeader);
/*
PURPOSE:

  Writes Log header to file and stdout

INPUTS:

  All  : All structure
  PreHeader : optional string to write before header (can be NULL)

OUTPUTS: 

  None

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_WriteLogEntry
extern int 
xf_WriteLogEntry( xf_All *All, xf_SolverData *SolverData, xf_Vector *U);
/*
PURPOSE:

  Writes Log entry to file and stdout

INPUTS:

  All  : All structure
  SolverData :  contains iIter, CFL, etc.
  U : state vector

OUTPUTS: 

  None

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_WriteLogLine
extern int 
xf_WriteLogLine( xf_All *All, const char *line);
/*
PURPOSE:

  Writes a specified line to the log file

INPUTS:

  All  : All structure
  line : string to write, newline included

OUTPUTS: 

  None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_LogUnsteadyOutputs
extern int 
xf_LogUnsteadyOutputs( xf_All *All);
/*
PURPOSE:

  Writes unsteady outputs to the log file.  Does not calculate them,
  but instead only writes their Value fields.  Should be called at the
  end of an unsteady run, when these values are filled in.

INPUTS:

  All  : All structure

OUTPUTS: 

  None

RETURN:

  Error Code
*/


#endif // end ifndef _xf_Log_h
