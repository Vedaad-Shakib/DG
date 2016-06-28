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

#ifndef _xf_MeshToolsStruct_h
#define _xf_MeshToolsStruct_h 1

/*
  FILE:  xf_MeshToolsStruct.h

  This file contains structures for working with MeshTools functions

*/


/* Structure for storing geometry Jacobian at possibly nq points.
   Useful for simplifying reallocation bookkeeping. */
typedef struct
{
  int nq;     /* number of points (usually quad points) for which the
		 vectors have meaningful data */
  int nqmax;  /* number of points (usually quad points) for which the
		 vectors are allocated (>= nq) */
  int dim;    /* dimension of geom Jacobian */
  real *detJ; /* determinant of Jacobian at the points [nq] */
  real *J;    /* Unrolled Jacobian matrices at the points [dim*dim*nq], 
		 stored as (2d e.g.):

		   J00,J01,J10,J11, J00,J01,J10,J11, ...
		       q = 0             q = 1

	      */
  real *iJ;   /* Unrolled inverse Jacobians (same storage as J) */

  real *GPhi;  /* gradients of geometry approximation basis functions */
  enum xfe_BasisType GPhi_Basis; /* basis used for gradients*/
  int GPhi_Order; /* order used for gradients*/

  real *T;  /* temporary storage (for speedup, to avoid reallocs) */
  int Tsize;  /* temporary storage size */


  unsigned int AllocFlag;  
  /* bit flag for which of detJ, J, iJ are
     allocated (consistently with nqmax).
     
     e.g.
     
     0  = 000:  none are alloced
     1  = 001:  detJ is alloced, others not
     6  = 110:  J and iJ are alloced, detJ is not

     use the bit masks (below) in practice
  */
}
xf_JacobianData;

// bit masks for JacobianData AllocFlag
#define xfb_detJ 1
#define xfb_J    2
#define xfb_iJ   4


/* Structure for storing geometry normal at possibly nq points.
   Useful for simplifying reallocation bookkeeping. */
typedef struct
{
  int nq;     /* number of points (usually quad points) for which the
		 vector has meaningful data */
  int dim;    /* dimension of geom normal */
  int ndqmax; /* size (dim * nqmax) for which the vector is allocated
		 (>= nq*dim) */
  real *n;    /* Unrolled normal at the points [dim*nq], 
		 stored as (2d e.g.):

		   nx0,ny0,  nx1,ny1, ...
		    q = 0     q = 1
	      */
}
xf_NormalData;


/* Element geometry data enumerated type */
enum xfe_ElemGeomType{
  xfe_EGVolume,
  xfe_EGSurfArea,
  xfe_ElemGeomLast
};

/* corresponding names */
static char *xfe_ElemGeomName[xfe_ElemGeomLast] = {
  "Volume",
  "SurfArea"
};

/* Element quality/statistical info enumerated type */
enum xfe_ElemStatType {
  xfe_ESVolume,
  xfe_ESSurfArea,
  xfe_ESAspectRatio,
  xfe_ESLast
};

/* corresponding names */
static char *xfe_ElemStatName[xfe_ESLast] = {
  "Volume",
  "SurfArea",
  "AspectRatio"
};

/* Hash list for edges */
typedef struct
{
  int n0, n1;
  int idata;
  int next;
}
xf_EdgeHash;


/* Search vectors for quick element lookup by position */
typedef struct
{
  int  dim;         // number of dimensions (0,1,2)
  int  nBin[3];     // number of bins in dim d, power of 2
  real *xBin[3];    // xBin[d][bin] = bin interval locations
  int  *nElem[3];   // nElem[d][bin] = number of elems in bin
  int  **Elem[3];   // Elem[d][bin][2*i+(0,1)] = ith (egrp,elem)
  int  *buf;        // buffer for using this structure
}
xf_ElemSearchStruct;


#endif // end ifndef _xf_MeshToolsStruct_h
