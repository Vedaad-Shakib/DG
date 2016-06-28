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

#ifndef _xf_MeshParallel_h
#define _xf_MeshParallel_h 1

/*
  FILE:  xf_MeshParallel.h

  This file contains the headers for parallelization functions dealing
  with the Mesh structure.

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ParallelizeMesh
extern int 
xf_ParallelizeMesh(xf_Mesh *Mesh_Glob, xf_Mesh *Mesh, int **ElemWeight,
                   int *ConnectWeight);
/*
PURPOSE:

  Parallelizes the Mesh data structure.  Each proc gets a submesh
  consisting of a chunk of relatively contiguous elements, as
  determined by the graph partitioner, and a halo of nearest-neighbor
  elements.  The entire mesh initially resides in Mesh_Glob on proc 0.

INPUTS:
 
  Mesh_Glob : the global Mesh structure, residing on proc 0
  ElemWeight : If not NULL, the element weights will be used in 
               the mesh partitioning
  ConnectWeight : If not NULL, the interelement connection (iface) 
                  weight will be used in the mesh partitioning.

OUTPUTS: 

  Mesh : pointer to resultant submesh on each proc.  Each proc must
         have this allocated before the call. 

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_UnParallelizeMesh
extern int 
xf_UnParallelizeMesh( xf_Mesh *Mesh, xf_Mesh *Mesh_Glob);
/*
PURPOSE:

  UnParallelizes the Mesh data structure.  Each proc inputs its submesh
  The entire assembled mesh is put into Mesh_Glob on proc 0.

INPUTS:
 
  Mesh : the submesh on each proc.  Relevant on all procs

OUTPUTS: 

  Mesh_Glob : the global Mesh structure, residing on proc 0

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_BcastMesh
extern int 
xf_BcastMesh(xf_Mesh **pMesh);
/*
 PURPOSE:
 
 Creates a copy of the Mesh data structure of the root processor 
 on the other processors. 
 
 INPUTS:
 
 pMesh : Pointer to the mesh structure on each processor.
 
 OUTPUTS: 
 
 None. pMesh gets modified on non-root processors.
 
 RETURN:
 
 Error Code
 */


#endif // end ifndef _xf_MeshParallel_h
