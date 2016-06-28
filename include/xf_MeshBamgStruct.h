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

#ifndef _xf_MeshBamgStruct_h
#define _xf_MeshBamgStruct_h 1

/*
  FILE:  xf_MeshBamgStruct.h

  This file contains structures for working with BAMG through xflow.

*/


/* This structure stores points for defining a BAMG geometry */
typedef struct
{
  int nVert;      // number of vertices
  real *Vert;     // coordinates of vertices
  real *Tangent;  // tangents at nodes

  int nEdge;      // number of edges
  int *Edge;      // 3*nEdge list: [n1, n2, boundary_flag]

  int nCorner;    // number of corners
  int *Corner;    // vertices that are corners

} 
xf_BamgGeomPoints;


#endif // end ifndef _xf_MeshBamgStruct_h
