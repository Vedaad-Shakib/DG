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

#ifndef _xf_MeshMotionIO_h
#define _xf_MeshMotionIO_h 1

/*
  FILE:  xf_MeshMotionIO.h

  This file contains headers for mesh motion IO

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ParallelizeMeshMotion
extern int 
xf_ParallelizeMeshMotion( xf_MeshMotion *Motion);
/*
PURPOSE:

  Parallelizes the MeshMotion data structure.  Each proc gets a copy of the
  MeshMotion data structure that initially only resides on proc 0.

INPUTS:

  Motion : pointer to MeshMotion structure;
           structure must exist (via CreateMeshMotion) on all procs

OUTPUTS: 

  None: Motion is parallelized

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadMeshMotionFile
extern int 
xf_ReadMeshMotionFile(char *MotionFile, char **InputStrings, 
		      xf_MeshMotion *Motion);
/*
PURPOSE:

  Reads a MeshMotion data structure from a text file or from
  InputStrings that are assumed to represent lines of the file.

INPUTS:

  MotionFile : name of file to read (optional)
  InputStrings : list of strings representing lines of file
                 (must be specified if MotionFile is NULL)

OUTPUTS: 

  Motion : MeshMotion structure.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Definition: xf_WriteMeshMotionBinarySerial
extern int 
xf_WriteMeshMotionBinarySerial( xf_MeshMotion *Motion, FILE *fid);
/******************************************************************/
//   FUNCTION Prototype: xf_WriteMeshMotionBinary
extern int 
xf_WriteMeshMotionBinary(xf_MeshMotion *Motion, FILE *fid);
/*
PURPOSE:

  Writes Motion to a binary file

INPUTS:

  Motion : Motion structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ReadMeshMotionBinarySerial
extern int 
xf_ReadMeshMotionBinarySerial(  FILE *fid, xf_MeshMotion *Motion);
/******************************************************************/
//   FUNCTION Prototype: xf_ReadMeshMotionBinary
extern int 
xf_ReadMeshMotionBinary( FILE *fid, xf_MeshMotion *Motion);
/*
PURPOSE:

  Reads Motion from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  Motion : MeshMotion structure to read

RETURN:

  Error Code
*/


#endif // end ifndef _xf_MeshMotionIO_h


