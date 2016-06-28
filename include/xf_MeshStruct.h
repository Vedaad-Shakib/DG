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

#ifndef _xf_MeshStruct_h
#define _xf_MeshStruct_h 1

/*
  FILE:  xf_MeshStruct.h

  This file contains the xflow mesh data structure

*/

#include "xf.h"
#include "xf_BasisStruct.h"
#include "xf_QuadStruct.h"
#include "xf_MeshHangStruct.h"
#include "xf_MeshMotionStruct.h"
#include "xf_MeshBamgStruct.h"


/* Face structure to simplify indexing into groups and face #s */
typedef struct
{
  int Group;
  /* xf_INTERIORFACE : if face is interior  (negative, defined in xf.h)
     xf_NULLFACE     : if face is null      (negative, defined in xf.h)
     BFGroup #       : otherwise; */

  int Number;
  /* number of interior or boundary face within Group */

} 
xf_Face;


/* CutFaceData stores quadrature points and other info for cut faces */
typedef struct
{
  xf_QuadData *QuadData; /* quad points, weights */
  xf_Face OrigFace;  /* Background mesh face */
  int GeomIndex; /* Index (e.g.) for identifying the cut face in the
		    geometry kernel */
}
xf_CutFaceData;


/* Interior face structure: stores element connectivity info.  The
   designation of left (L) and right (R) is arbitrary (no
   right-hand-rule assumptions).
*/
typedef struct
{
  int ElemGroupL;  /* ElemGroup number of Left element */
  int ElemL;       /* Elem number within ElemGroupL of Left element */
  int FaceL;       /* Local face number within Left element */
  int ElemGroupR;  /* ElemGroup number of Right element */
  int ElemR;       /* Elem number within ElemGroupR of Right element */
  int FaceR;       /* Local face number within Right element */

  int OrientL;  /* Index specifying how the face ref space is oriented
		   w.r.t ref elem L */
  int OrientR;  /* Index specifying how the face ref space is oriented
		   w.r.t ref elem R */

  int HangNumber; 
  /* If nonzero, indicates that this interior face is part of a
     hanging-node refinement.  If negative (positive), means that the
     Left (Right) element is the coarse element that must deal with
     the "hanging" face. abs(HangNumber) indicates the manner in which
     the face is positioned relative to the element. */

  xf_CutFaceData *CutFaceData;
  /* Stores quad points and interpolation information if IFace is cut.
     Null for non-cut faces. */
} 
xf_IFace;


/* Boundary face structure: identifies a boundary face within the mesh */
typedef struct
{
  int ElemGroup;  /* ElemGroup number of adjacent element */
  int Elem;       /* Elem number within the ElemGroup */
  int Face;       /* Local face number within the element */

  int Orient;  /* Index specifying how the face ref space is oriented
		  w.r.t the ref elem */

  xf_CutFaceData *CutFaceData;
  /* Stores quad points and interpolation information if BFace is cut.
     Null for non-cut faces. */
} 
xf_BFace;


/* Boundary face group structure */
typedef struct
{
  char *Title;       /* Name of the boundary face group */
  int nBFace;        /* Number of boundary faces in this boundary face group */
  xf_BFace *BFace;   /* Array of boundary faces in this group */  
} 
xf_BFaceGroup;



/* CutElemData stores quadrature points and other info for cut elems */
typedef struct
{
  xf_QuadData *QuadData; /* quad points, weights */
  enum xfe_BasisType QBasis;  /* Basis for interpolation element (QOrder = 1) */
  real **QCoord; /* Coordinates of interpolation element nodes in
		   original element ref space [dim reals per node] */

  int GeomIndex; /* Index (e.g.) for identifying the cut elem in the
		    geometry kernel */
}
xf_CutElemData;


/* Element Group structure: all elements within a group should have
   the same CutFlag, geometry type (QBasis and QOrder), nFace,
   and nNode. */
typedef struct 
{

  enum xfe_Bool CutFlag;
  /* If true, element group consists of cut elements */

  enum xfe_BasisType QBasis;
  /* Type of basis used for interpolating the element geometry */
  
  int QOrder;
  /* Order of geometry basis (Q) */

  int nElem;
  /* Number of elements in this element group */

  int *nFace;
  /* nFace[elem] = number of faces in element elem = 0 to nElem-1 */

  xf_Face **Face;
  /*
    Face[elem][i] = face identifier of i'th face of element elem
    elem = 0 to nElem-1
    i    = 0 to nFace[elem]-1
    
    Face[elem][i].Group  = -1/ibfgrp for interior/boundary faces
    Face[elem][i].Number = iface or bface number
  */

  int nNode;
  /* Number of nodes for each element in this group */

  int **Node;
  /* Node[elem][i] = global node number of i'th node of element elem
     i = 0 to nElement-1
     j = 0 to nNode-1    
  */

  xf_CutElemData *CutElemData;
  /* Stores quad points and interpolation information if CutFlag==True
     (elems in group are cut).  nElem of these structures are
     allocated.  Null for non-cut groups. */

} 
xf_ElemGroup;


/* Mesh Parallel Info structure */
typedef struct
{
  int nIFaceRegular;         
  /* Number of regular interior faces in the mesh, where a regular
     iface is one for which both adjacent elements lie on the current
     processor.  The alternative would be a "halo" iface, for which
     one of the adjacent elements is in a halo group.  Regular IFaces
     are listed first in the local mesh IFace list, so nIFaceRegular
     is used to know when to pause to allow inter-processor data
     communication to complete, when computing the iface residuals. */

  int **ElemLoc2Glob;
  /* ElemLoc2Glob[egrp][elem] stores the element number in the global
     (unparallelized) mesh corresponding to the local element [egrp,
     elem].  Note, the element group number in the global mesh is just
     egrp. */
  
  int *IFaceLoc2Glob;
  /* IFaceLoc2Glob[iiface] stores the global-mesh interior face number
     of local interior face iiface. */

  int *NodeLoc2Glob;
  /* NodeLoc2Glob[iNode] stores the global-mesh node number of local
     node iNode. */

  int **nSendElem;
  /* nSendElem[egrp][iProc] stores the number of elements in egrp that
     serve as halo cells for processor iProc.  It is used jointly with
     SendElem. */

  int ***SendElem;
  /* SendElem[egrp][iProc][i], where 0 <= i < nSendElem[egrp][iProc],
     stores the element numbers of those elements in egrp that serve
     as halo cells for processor iProc.  When it comes time to
     exchange data, the current processor will send data for these
     nSendElem[egrp][iProc] elements to the appropriate halo element
     group (Mesh->nElemGroup + egrp) on processor iProc. */

  int **nRecvElem;
  /* nRecvElem[egrp][iProc] stores the number of elements that halo
     element group Mesh->nElemGroup + egrp expects to receive from
     element group egrp on processor iProc during a data exchange.  It
     equals iProc's value of nSendElem[egrp][currentProc].  NOTE: the
     halo group elements are assumed to be ordered such that all the
     elems residing on iProc 0 come before the elements on iProc 1,
     etc. */

} 
xf_MeshParallelInfo;

/* Periodicity types*/
enum xfe_PeriodicityType { 
  xfe_PeriodicityTranslational, 
  xfe_PeriodicityRotational,
  xfe_PeriodicityLast
};
static char *xfe_PeriodicityName[xfe_PeriodicityLast] = {
  "Translational",
  "Rotational",
};


/* Periodic Group structure */
typedef struct
{
  enum xfe_PeriodicityType Periodicity; /* e.g. translational, rotational */
  int nPeriodicNode;  /* number of periodic nodes */
  int *PeriodicNode;  /* unrolled vector of matching node pairs */  
} 
xf_PeriodicGroup;


/*--------------- Mesh data structure definition  ----------------*/
typedef struct
{
  
  int Dim;
  /* Dimension of mesh: 1, 2, or 3 */

  int nNode;
  /* Number of nodes in the mesh */

  real **Coord;
  /* Coordinates of the nodes:
     Coord[iNode][d], 0 <= iNode < nNode, 0 <= d < Dim
  */

  int nIFace;
  /* Number of interior faces */

  xf_IFace *IFace;
  /* Array of interior face connectivity structures */

  int nBFaceGroup;
  /* Number of boundary face groups, including embedded groups if
     present */

  xf_BFaceGroup *BFaceGroup;
  /* Array of boundary face group objects */

  int nElemGroup;
  /* Number of element groups */

  xf_ElemGroup *ElemGroup;
  /* Array of element group objects */

  int nPeriodicGroup;
  /* Number of periodic groups for periodic BCs */
  
  xf_PeriodicGroup *PeriodicGroup;
  /* Array of periodic group objects */
  
  xf_MeshParallelInfo *ParallelInfo;
  /* Stores how this mesh portion fits together with neighboring
     portions on other processors.  NULL for serial runs. */

  void *BackgroundMesh;
  /* Background mesh for cut-cell meshes */

  xf_MeshMotion *Motion;
  /* Structure defining mesh motion in time */

}
xf_Mesh;


#endif // end ifndef _xf_MeshStruct_h
