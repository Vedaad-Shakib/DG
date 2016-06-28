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

#ifndef _xf_GeomIO_h
#define _xf_GeomIO_h 1

/*
  FILE:  xf_GeomIO.h

  This file contains headers for mesh motion IO

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ParallelizeGeom
extern int 
xf_ParallelizeGeom( xf_Geom *Geom);
/*
PURPOSE:

  Parallelizes the Geom data structure.  Each proc gets a copy of the
  Geom data structure that initially only resides on proc 0.

INPUTS:

  Geom : pointer to Geom structure;
           structure must exist (via CreateGeom) on all procs

OUTPUTS: 

  None: Geom is parallelized

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadGeomFile
extern int 
xf_ReadGeomFile(char *GeomFile, char **InputStrings, xf_Geom *Geom);
/*
PURPOSE:

  Reads a Geom data structure from a text file or from
  InputStrings that are assumed to represent lines of the file.

INPUTS:

  GeomFile : name of file to read (optional)
  InputStrings : list of strings representing lines of file
                 (must be specified if GeomFile is NULL)

OUTPUTS: 

  Geom : Geom structure.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteGeomBinary
extern int 
xf_WriteGeomBinary(xf_Geom *Geom, FILE *fid);
/*
PURPOSE:

  Writes Geom to a binary file

INPUTS:

  Geom : Geom structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadGeomBinary
extern int 
xf_ReadGeomBinary( FILE *fid, xf_Geom *Geom);
/*
PURPOSE:

  Reads Geom from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  Geom : Geom structure to read

RETURN:

  Error Code
*/


#endif // end ifndef _xf_GeomIO_h


