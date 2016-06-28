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

/*
  FILE:  xf_GeomAnalytical.c

  This file contains functions for working with analytical geometry
  components.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Mesh.h"
#include "xf_Param.h"
#include "xf_Math.h"

/* List of understood object types */
enum xfe_GCAObjectType{
  xfe_GCAObjectNone,
  xfe_GCAObjectCylinder,
  xfe_GCAObjectLast
};

/* Corresponding names */
static char *xfe_GCAObjectName[xfe_GCAObjectLast] = {
  "None",
  "Cylinder"
};


/******************************************************************/
//   FUNCTION Definition: xf_GCAParams_Cylinder
static int 
xf_GCAParams_Cylinder( xf_GeomCompAnalytical *GCA, real **pxc, real *pR)
{
  int ierr;
  int d, nset;
  // Real parameters: enumerated type, names, defaults
  enum {XCenter, YCenter, Radius, RPLast};
  char *RPName[RPLast] = {"XCenter", "YCenter", "Radius"};
  real RPDef[RPLast] = {0., 0., 1.};
  real *RP;

  // fill in if do not have parameters in shortcut storage
  if (GCA->RParam == NULL){
    ierr = xf_Error(xf_AllocFillRParam(GCA->KeyValue, RPLast, RPName, RPDef, 
				       &GCA->RParam, &nset));
    if (ierr != xf_OK) return ierr;
  }
  RP = GCA->RParam;

  (*pxc) = RP+0;       // center
  (*pR)  = RP[Radius]; // radius

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GCAProject_Cylinder
static int 
xf_GCAProject_Cylinder( xf_GeomCompAnalytical *GCA, int dim, real *x)
{
  // Object-specific projection function
  int ierr, i;
  real *xc, R, r;

  // dimension check
  if (dim != 2) return xf_Error(xf_INPUT_ERROR);

  ierr = xf_Error(xf_GCAParams_Cylinder(GCA, &xc, &R));
  if (ierr != xf_OK) return ierr;
    
  // projection: maintains theta, radius is normalized to R
  r = xf_Distance(xc, x, 2);
  for (i=0; i<2; i++) x[i] = xc[i] + (x[i]-xc[i])*R/r;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ProjectToGeomComp_Analytical
int 
xf_ProjectToGeomComp_Analytical( xf_GeomCompAnalytical *GCA, int dim, real *x)
{
  int ierr;
  
  // Have we worked with this object before?
  if (GCA->Object <= 0){
    // determine object type
    ierr = xf_Error(xf_GetKeyValueEnum(GCA->KeyValue, "Object", xfe_GCAObjectName,
				       (int ) xfe_GCAObjectLast, (int *) &GCA->Object));
    if (ierr != xf_OK) return ierr;
  }

  if (GCA->Object <= 0) return xf_Error(xf_INPUT_ERROR);

  // call object-specific function
  switch (GCA->Object){
  case xfe_GCAObjectCylinder:
    ierr = xf_Error(xf_GCAProject_Cylinder(GCA, dim, x));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GCAPointsOn_Cylinder
static int 
xf_GCAPointsOn_Cylinder( xf_GeomCompAnalytical *GCA, int dim, 
		      int np, real dl, enum xfe_GeomSpacingType Spacing, 
		      int *pnout, real **px)
{
  // Object-specific "points on" function
  int ierr, i, N;
  real *xc, R;
  real theta, dtheta;

  // dimension check
  if (dim != 2) return xf_Error(xf_INPUT_ERROR);

  ierr = xf_Error(xf_GCAParams_Cylinder(GCA, &xc, &R));
  if (ierr != xf_OK) return ierr;
   
  // how many points do we need?
  if (np > 0) N = np;
  else N = 2.*PI*R/dl;
  (*pnout) = N;

  // determine dtheta
  dtheta = 2.*PI/N;

  // allocate memory for points
  ierr = xf_Error(xf_Alloc( (void **) px, 2*N, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // calculate points (clockwise ordering)
  for (i=0, theta=0.; i<N; i++, theta-=dtheta){
    (*px)[2*i+0] = xc[0]+R*cos(theta);
    (*px)[2*i+1] = xc[1]+R*sin(theta);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_PointsOnGeom_Analytical
int 
xf_PointsOnGeom_Analytical( xf_GeomCompAnalytical *GCA, int dim, 
			    int np, real dl, enum xfe_GeomSpacingType Spacing, 
			    int *pnout, real **px)
{
  int ierr;

  // Have we worked with this object before?
  if (GCA->Object <= 0){
    // determine object type
    ierr = xf_Error(xf_GetKeyValueEnum(GCA->KeyValue, "Object", xfe_GCAObjectName,
				       (int ) xfe_GCAObjectLast, (int *) &GCA->Object));
    if (ierr != xf_OK) return ierr;
  }

  if (GCA->Object <= 0) return xf_Error(xf_INPUT_ERROR);

  // call object-specific function
  switch (GCA->Object){
  case xfe_GCAObjectCylinder:
    ierr = xf_Error(xf_GCAPointsOn_Cylinder(GCA, dim, np, dl, Spacing, pnout, px));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}
