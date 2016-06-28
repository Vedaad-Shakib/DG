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

#ifndef _xf_QuadStruct_h
#define _xf_QuadStruct_h 1

/*
  FILE:  xf_QuadStruct.h

  This file contains the xflow Solver structures

*/


#include "xf.h"


/* This enumerated type is used to describe what kind of quadrature
   data is stored/pointed-to by an xf_QuadData stucture */
enum xfe_QuadDataType{
  xfe_QuadDataGeneric, /* The quad points are generic (e.g. for a
			  reference triangle) and can be released or
			  (re)alloc'ed at will. */
  xfe_QuadDataSpecific,/* The quad points are elem or face specific,
			  and should not be freed or (re)alloc'ed
			  during routine computations.*/
  xfe_QuadDataLast
};

static char *xfe_QuadDataName[xfe_QuadDataLast] = {
  "Generic", 
  "Specific"
};


/* Quadrature point data structure */
typedef struct
{
  enum xfe_QuadDataType Type;  /* Whether quad data is generic or
				  specific (see above) */
  
  enum xfe_ShapeType Shape; /* Element shape on which the quad rule is
			       defined.  Used in generic case. */

  int Order; /* Order up to which this quad rule can integrate. */

  int nquad;    /* number of quad points */
  int qdim;     /* dim of xquad and nvec */
  real *xquad;  /* quad point positions:
		     for a regular face: in face ref space
		     for a regular elem: in elem ref space
		     for a cut face: in face ref space
		     for an embedded face: in *elem* ref space
		     for a cut elem: in *interpolation* elem ref space
		   size = [nquad * qdim] */
  real *wquad;  /* quad point weights */
  real *nvec;   /* normals at quad points, for embedded faces
		   size = [nquad * qdim] */
  
}
xf_QuadData;


 

#endif // end ifndef _xf_QuadStruct_h
