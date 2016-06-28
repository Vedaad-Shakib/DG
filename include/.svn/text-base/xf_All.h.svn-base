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

#ifndef _xf_All_h
#define _xf_All_h 1

/*
  FILE:  xf_All.h

  This file contains the headers for functions dealing with the All structure.

*/

/******************************************************************/
//   FUNCTION Prototype: xf_CreateAll
extern int 
xf_CreateAll( xf_All **All, enum xfe_Bool DefaultFlag);
/*
PURPOSE:

  Creates an All structure and all of its children.  Memory is
  allocated and initial (zero) values are set.

INPUTS:

  All : pointer to All
  DefaultFlag : true to set default parameters

OUTPUTS: 

  None: All is allocated

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyAll
extern int 
xf_DestroyAll( xf_All *All);
/*
PURPOSE:

  Destroys an All structure and all of its children.  Memory is
  released

INPUTS:

  All : All structure

OUTPUTS: 

  None: All is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_CheckAll
extern int 
xf_CheckAll( xf_All *All);
/*
PURPOSE:

  Checks components of All structure for consistency

INPUTS:

  All : All structure

OUTPUTS: 

  None: All is checked

RETURN:

  Error Code (if any) from check
*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteAllBinary
extern int 
xf_WriteAllBinary( xf_All *All, const char *fname);
/*
PURPOSE:

  Writes All to a binary file

INPUTS:

  All : All structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ReadAllBinary
extern int 
xf_ReadAllBinary( const char *fname, xf_All *All);
/*
PURPOSE:

  Reads All from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  All : All structure to read.  Should have Mesh, Geom, etc. allocated.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadAllInputFile
extern int 
xf_ReadAllInputFile( const char *InputFile, xf_KeyValue *KeyValue, 
		     enum xfe_Bool DefaultFlag, xf_All **pAll);
/*
PURPOSE:

  Reads All from InputFile.  Extension on InputFile identifies the
  type of file.

INPUTS:

  InputFile : name of file to read (with extension)
  KeyValue : pointer to key value list of optional parameters (optional)
  DefaultFlag : if True, any existing key-value pairs in the read All
                file will be reset to default values.

OUTPUTS: 

  (*pAll) : All structure read in.  Should not be pre-allocated.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadAllFromJobFile
extern int 
xf_ReadAllFromJobFile( const char *jobFileIn, enum xfe_Bool ReadEqnSet,
		       xf_All **pAll);
/*
PURPOSE:

  Reads All from key-value pairs given in a job file.  All and
  children should NOT be pre-allocated (via call to CreateAll).  An
  eqnset file is read in only if ReadEqnSet is True.

INPUTS:

  jobFileIn : name of .job file to read (with or without extension)
  ReadEqnSet : if True, .eqn file will also be read

OUTPUTS: 

  (*pAll) : All structure read in.  Should not be pre-allocated.

RETURN:

  Error Code
*/



#endif // end ifndef _xf_All_h
