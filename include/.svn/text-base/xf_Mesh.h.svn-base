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

#ifndef _xf_Mesh_h
#define _xf_Mesh_h 1

/*
  FILE:  xf_Mesh.h

  This file contains the headers for functions dealing with the Mesh structure.

*/

#include "xf_MeshParallel.h"

/******************************************************************/
//   FUNCTION Prototype: xf_CreateMesh
extern int 
xf_CreateMesh( xf_Mesh **Mesh);
/*
PURPOSE:

  Creates an Mesh structure and all of its children.  Memory is
  allocated and initial (zero) values are set.

INPUTS:

  Mesh : pointer to Mesh

OUTPUTS: 

  None: Mesh is allocated

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_DestroyBFaceGroup
extern int 
xf_DestroyBFaceGroup( xf_BFaceGroup *BFaceGroup);
/*
PURPOSE:

  Destroys a boundary face group (self pointer not released)

INPUTS:

  BFaeGroup : pointer to boundary face group

OUTPUTS: 

  None: BFaceGroup is destroyed

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyMesh
extern int 
xf_DestroyMesh( xf_Mesh *Mesh);
/*
PURPOSE:

  Destroys a Mesh structure and all of its children.  Memory is
  released

INPUTS:

  Mesh : pointer to Mesh

OUTPUTS: 

  None: Mesh is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_InitElemGroup
extern void
xf_InitElemGroup( xf_ElemGroup *ElemGroup);
/*
PURPOSE:

  Initializes a pointer to an element group to default values.

INPUTS:

  ElemGroup : pointer to element group structure

OUTPUTS: 

  None

RETURN: None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_InitIFace
extern void
xf_InitIFace( xf_IFace *IFace);
/*
PURPOSE:

  Initializes a pointer to an IFace to default values.

INPUTS:

  IFace : pointer to element group structure

OUTPUTS: 

  None

RETURN: None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_InitBFace
extern void
xf_InitBFace( xf_BFace *BFace);
/*
PURPOSE:

  Initializes a pointer to a BFace structure to default values.

INPUTS:

  BFace : pointer to element group structure

OUTPUTS: 

  None

RETURN: None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_AddInteriorFace
extern void
xf_AddInteriorFace(xf_Mesh *Mesh, int egL, int elemL, int faceL, 
		   int egR, int elemR, int faceR);
/*
PURPOSE:

  Adds an interior face to Mesh, with specified left and right
  elements.  Note: assumes Mesh->IFace has been adequately allocated!


INPUTS:


  Mesh : input Mesh
  egL, elemL, faceL : left element info
  egR, elemR, faceR : right element info

OUTPUTS: 

  None, Mesh has one face added to it

RETURN: None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadGriFile
extern int 
xf_ReadGriFile( const char *InputFile, char **InputStrings, xf_Mesh *Mesh);
/*
PURPOSE:

  Reads a .gri InputFile into structure Mesh.  The input format for
  the .gri file is described in the UserGuide documentation.
  Alternately, this function can read the file from a sequence of
  InputStrings which contain the lines of a .gri file.  Note, one and
  only one of InputFile or InputStrings must be NULL.

INPUTS:

  InputFile: name of file to read (or NULL)
  InputStrings: strings containing lines of input file (or NULL)

OUTPUTS: 

  Mesh : pointer to Mesh structure; must have been allocated before
         the call to this function.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteGriFile
extern int 
xf_WriteGriFile( xf_Mesh *Mesh, char *OutputFile );
/*
PURPOSE:

  Writes a .gri OutputFile from structure Mesh.

INPUTS:

  Mesh : mesh structure to write

OUTPUTS: 

  OutpuFile : name of ascii text file to write to

RETURN:

  Error Code
*/


// Gmsh support functions
#include "xf_MeshGmsh.h"
// Fluent msh support functions
#include "xf_MeshFmsh.h"
// Bamg support functions
#include "xf_MeshBamg.h"


/******************************************************************/
//   FUNCTION Prototype: xf_WriteMeshBinary
extern int 
xf_WriteMeshBinary( xf_Mesh *Mesh, FILE *fid);
/*
PURPOSE:

  Writes Mesh to a binary file

INPUTS:

  Mesh : mesh structure to write
  fid : file to write to

OUTPUTS: 

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadMeshBinary
extern int 
xf_ReadMeshBinary( FILE *fid, xf_Mesh *Mesh);
/*
PURPOSE:

  Reads Mesh from a binary file

INPUTS:

  fid : file from which to read

OUTPUTS: 

  Mesh : Mesh structure to read

RETURN:

  Error Code
*/


#endif // end ifndef _xf_Mesh_h
