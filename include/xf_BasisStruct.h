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

#ifndef _xf_BasisStruct_h
#define _xf_BasisStruct_h 1

/*
  FILE:  xf_BasisStruct.h

  This file contains enumerated types of basis functions

*/


/* Simply for preventing excessive reallocation in certain routines.
   Can increase when necessary */
#define xf_MAXQ1FACENODE 4
#define xf_MAXQ1NODE 8
#define xf_MAXLOCFACE 6
#define xf_MAXEDGES 12


/* BasisType describes the class of functions used for general
   interpolation.  The interpolation could be of a solution variable
   or of the physical shape of the element (i.e. geometry).  Thanks to
   the DG setting, it is possible to use different basis types for
   geometry vs solution interpolation.  For example, the geometry of a
   triangle can be interpolated using TriLagrange basis functions,
   while the solution can be interpolated using QuadLagrangeGauss
   basis functions.  Of course, the geometry interpolation must be
   consistent with the shape of the element; for example, using a Quad
   class of functions for interpolating the geometry of a triangle
   does not make sense, as the extra nodes are not available.
*/
enum xfe_BasisType{

  xfe_SegLagrange,
   
  /* For a line segment element: Lagrange interpolation with equi-spaced
     nodes.  The reference coords of the nodes for order p
     interpolation are given by:
 
     for i=0:p
       node(i) = i/(p+1)
     end

  */

  xfe_SegLagrangeGauss,
  
  /* Similar to SegLagrange, except using Gauss points instead of
     equi-spaced points in each direction.
  */

  xfe_TriLagrange,

  /* For a triangle element: Lagrange interpolation with equi-spaced
     nodes.  The reference coords of the nodes for order p
     interpolation are given by:

     i=0; 
     for iy=0:p
       for ix=0:p-iy
         node(i) = (ix,iy)
	 i = i + 1
       end
     end

  */

  xfe_TriHierarch,

  /* For a triangle element: Hierarchical interpolation with basis
     functions as described by Solin etal [Solin, Segeth, DolevZel,
     "High-Order Finite Element Methods", 2003]. The basis is only
     hierarchical for p>=1, as p=0 is simply the constant function.
  */

  xfe_TetLagrange,

  /* For a tetrahedron element: Lagrange interpolation with equi-
     spaced nodes.  The reference coords of the nodes are given by:

     i=0; 
     for iz=0:p
       for iy=0:p-iz
         for ix=0:p-iz-iy
           node(i) = (ix,iy,iz)
	   i = i + 1
	 end
       end
     end

  */
  xfe_TetHierarch,

  /* For a tetrahedron element: Hierarchical interpolation with basis
     functions as described by Solin etal [ref above]. The basis is only
     hierarchical for p>=1, as p=0 is simply the constant function.
  */

  xfe_QuadLagrange,

  /* For a quadrilateral element: Lagrange interpolation with
     equispaced nodes.  The coords of the nodes are:

     i=0; 
     for iy=0:p
       for ix=0:p
         node(i) = (ix,iy)
	 i = i + 1
       end
     end
  
  */

  xfe_QuadLagrangeGauss,
  
  /* Similar to QuadLagrange, except using Gauss points instead of
     equi-spaced points in each direction.
  */

  xfe_QuadLegendre,

  /* Basis defined on 2D quads, tensor products basis, similar to
   * above two bases, but Lagrangian polynomial is replaced by 
   * Legendre polynomial, which simply Yu's entropy bounding idea
   */

  xfe_HexLagrange,

  /* For a hexahedron element: Lagrange interpolation with
     equispaced nodes.  The coords of the nodes are:

     i=0;
     for iz = 0:p
       for iy=0:p
         for ix=0:p
           node(i) = (ix,iy,iz)
	   i = i + 1
	 end
       end
     end
  
  */

  xfe_HexLagrangeGauss,

  /* Similar to HexLagrange, except using Gauss points instead of
     equi-spaced points in each direction.
  */ 

  xfe_HexLegendre,

  /* Extend xfe_QuadLegendre to 3D */

  xfe_BasisLast
};

/* corresponding names */
static char *xfe_BasisName[xfe_BasisLast] = {
  "SegLagrange",
  "SegLagrangeGauss",
  "TriLagrange",
  "TriHierarch",
  "TetLagrange",
  "TetHierarch",
  "QuadLagrange",
  "QuadLagrangeGauss",
  "QuadLegendre",
  "HexLagrange",
  "HexLagrangeGauss",
  "HexLegendre"
};
  


/* Physical shape types for reference elements */
enum xfe_ShapeType{
  xfe_Point,         /* Point */
  xfe_Segment,       /* Line segment  */
  xfe_Triangle,      /* Triangle      */
  xfe_Tetrahedron,   /* Tetrahedron   */
  xfe_Quadrilateral, /* Quadrilateral */
  xfe_Hexahedron,    /* Hexahedron    */
  xfe_ShapeLast
};

/* corresponding names */
static char *xfe_ShapeName[xfe_ShapeLast] = {
  "Point",
  "Segment",
  "Triangle",
  "Tetrahedron",
  "Quadrilateral",
  "Hexahedron",
};


/* Structure for storing basis functions at nq points.  Useful for
   simplifying reallocation bookkeeping. */
typedef struct
{
  enum xfe_Bool InTable; /* Flag to identify whether this basis exists
			    in a lookup table (for memory releasing) */ 
  enum xfe_BasisType Basis;  /* interpolation basis for which the
				values of Phi are computed */
  enum xfe_Bool Needs_ReCalc; /* If true, means basis will be recalculated
				 next time it it required.  Used to explicitly
				 force a basis recalculation. */
  int Order;  /* interpolation order for which the values of Phi are
		 computed. */

  int nn;     /* number of basis functions for which the vectors have
		 useful data. */
  int nq;     /* number of points (usually quad points) for which the
		 vectors have meaningful data */
  int nnqmax;  /* number of (points times basis functions) for which
		  the vectors are allocated (>= nq*nn) */
  int dim;    /* dimension for the gradient */
  real *Phi;  /* basis functions, unrolled along nn first, then nq:

		   phi0,phi1,..phin,  phi0,phi1,..phin, ...
		       q = 0               q = 1

	         total size = [nq*nn] */
  real *GPhi; /* ref-space gradients of basis functions, unrolled
		 along nn first, then nq, then dim:

		   phi0,phi1,..phin,  phi0,phi1,..phin, ...
		    q = 0, dim = 0     q = 1, dim = 0

		   phi0,phi1,..phin,  phi0,phi1,..phin, ...
		    q = 0, dim = 1     q = 1, dim = 1

		   ...

	         total size = [dim*nq*nn] */
  real *gPhi; /* global space gradients of basis functions (same
		 storage as GPhi)*/

  real *HPhi; /* ref-space Hessians of basis functions, unrolled
		 along nn first, then nq, then dim2, then dim1:


		 HPhi{i,j,q,n} = second derivative of Phi{q,n}
		                 w.r.t. xref{i} and xref{j}

	         total size = [dim*dim*nq*nn] */

  unsigned int AllocFlag;  
  /* bit flag for which of Phi, GPhi, gPhi, HPhi are allocated
     (consistently with nqmax).
     
     e.g.
     
     0  = 000:  none of [gPhi, GPhi, Phi] are alloced
     1  = 001:  Phi is alloced, others not
     6  = 110:  GPhi and gPhi are alloced, Phi is not

     use the bit masks (below) in practice
  */
}
xf_BasisData;

// bit masks for BasisData AllocFlag
#define xfb_Phi   1
#define xfb_GPhi  2
#define xfb_gPhi  4
#define xfb_HPhi  8


/* Structure for storing pointers to BasisData structures for lookup. */
typedef struct
{

  int nShape;      /* max number of shape types for which the pointer
		      array is allocated */
  int nface;       /* max number of faces for which the pointer array
		      is allocated */
  int nfaceorient; /* max number of face orientations for which the
		      pointer array is allocated */
  
  xf_BasisData ****BasisData;
  /* BasisTable->BasisData[Shape][face][faceorient] is a pointer to a
     BasisData structure defined for the specific Shape, face, and
     faceorient.  The BasisData structure stores info on what Basis,
     and Order it is defined for; if these change, the BasisData
     structure is updated (in the function EvalBasisOnFaceUsingTable)
     and the most recent structure is stored in the table. */

  real *gPhi; /* global space gradients of basis functions (same
		 storage as GPhi)*/

}
xf_BasisTable;



#endif // end ifndef _xf_BasisStruct_h