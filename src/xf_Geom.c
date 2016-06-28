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
  FILE:  xf_Geom.c

  This file contains functions for working with the Geometry data
  sructure.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Mesh.h"
#include "xf_Param.h"
#include "xf_Math.h"
#include "xf_GeomAnalytical.h"
#include "xf_GeomSpline.h"
#include "xf_GeomPanel.h"



/******************************************************************/
//   FUNCTION Definition: xf_CreateGeom
int 
xf_CreateGeom( xf_Geom **pGeom){
  
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pGeom, 1, sizeof(xf_Geom)));
  if (ierr != xf_OK) return ierr;

  (*pGeom)->Dim   = 0;
  (*pGeom)->nComp = 0;
  (*pGeom)->Comp  = NULL;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyGeomCompAnalytical
static int 
xf_DestroyGeomCompAnalytical( xf_GeomCompAnalytical *GeomCompAnalytical){
  
  int ierr;

  /* Initialize key-value structure */
  ierr = xf_Error(xf_DestroyKeyValue(&GeomCompAnalytical->KeyValue));
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) GeomCompAnalytical->RParam);
  xf_Release( (void *) GeomCompAnalytical->IParam);

  xf_Release( (void *) GeomCompAnalytical);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyGeomCompSpline
static int
xf_DestroyGeomCompSpline( xf_GeomCompSpline *GeomCompSpline){
  
  int ierr;

  xf_Release( (void *) GeomCompSpline->X );
  xf_Release( (void *) GeomCompSpline->Y );
  xf_Release( (void *) GeomCompSpline->S );
  xf_Release( (void *) GeomCompSpline->XS);
  xf_Release( (void *) GeomCompSpline->YS);

  xf_Release( (void *) GeomCompSpline);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyGeomCompPanel
static int
xf_DestroyGeomCompPanel( xf_GeomCompPanel *GeomCompPanel){
  
  int ierr;

  xf_Release( (void *) GeomCompPanel->Coord );
  xf_Release( (void *) GeomCompPanel->Panels);

  xf_Release( (void *) GeomCompPanel);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DestroyGeom
int 
xf_DestroyGeom( xf_Geom *Geom){
  
  int i, ierr;
  xf_GeomComp *Comp;

  if (Geom == NULL) return xf_OK;

  for (i=0; i<Geom->nComp; i++){
    
    Comp = Geom->Comp + i;
    
    xf_Release ((void *) Comp->Name);
    xf_Release ((void *) Comp->BFGTitle);
    
    switch (Comp->Type){
    case xfe_GeomCompNone:
      break;
    case xfe_GeomCompAnalytical:
      ierr = xf_Error(xf_DestroyGeomCompAnalytical( (xf_GeomCompAnalytical *) Comp->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompSpline:
      ierr = xf_Error(xf_DestroyGeomCompSpline( (xf_GeomCompSpline *) Comp->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompPanel:
      ierr = xf_Error(xf_DestroyGeomCompPanel( (xf_GeomCompPanel *) Comp->Data));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    }
  }

  xf_Release((void *) Geom->Comp);

  xf_Release((void *) Geom);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ProjectToGeomComp
static int 
xf_ProjectToGeomComp( xf_GeomComp *Comp, int dim, real *x)
{
/*
PURPOSE:

  Projects 1 point, x, to 1 geometry component, Comp.  Overwrites data
  in x.

INPUTS:

  Comp     : Geometry component
  dim      : spatial dimension of x
  x        : on input, original location of point

OUTPUTS: 

  x        : on output, new (projected) location of point

RETURN:

  Error Code
*/
  int ierr;

  switch (Comp->Type){
  case xfe_GeomCompNone:
    break;
  case xfe_GeomCompAnalytical:
    ierr = xf_Error(xf_ProjectToGeomComp_Analytical( (xf_GeomCompAnalytical *) Comp->Data, 
						     dim, x));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_GeomCompSpline:
    ierr = xf_Error(xf_ProjectToGeomComp_Spline( (xf_GeomCompSpline *) Comp->Data, 
						 dim, x));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_GeomCompPanel:
    ierr = xf_Error(xf_ProjectToGeomComp_Panel( (xf_GeomCompPanel *) Comp->Data, 
						dim, x));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ProjectToGeom
int 
xf_ProjectToGeom( xf_Geom *Geom, int iComp, char *BFGTitle, 
		  int dim, int np, real *x)
{  
  int ierr;
  int ip, i, iC;
  real *x0, xproj[3], x1[3]; 
  real dist, mindist;

  // check dimension of geometry for consistency
  if (Geom->Dim != dim) return xf_Error(xf_INPUT_ERROR);

  // loop over np points
  for (ip=0; ip<np; ip++){

    x0 = x + ip*dim;

    for (i=0; i<dim; i++) xproj[i] = x0[i];
    mindist = 1e30;

    // loop over components
    for (iC=0; iC<Geom->nComp; iC++){

      // if component not in consideration, continue
      if ((iC != iComp) && (strcmp(BFGTitle, Geom->Comp[iC].BFGTitle) != 0)) continue;
  
      // store x in x1
      for (i=0; i<dim; i++) x1[i] = x0[i];

      // project x1 to geometry component
      ierr = xf_Error(xf_ProjectToGeomComp(Geom->Comp+iC, dim, x1));
      if (ierr != xf_OK) return ierr;
  
      // determine if minimum distance
      dist = xf_Distance(x0, x1, dim);
      if (dist < mindist){
	for (i=0; i<dim; i++) xproj[i] = x1[i];
	mindist = dist;
      }
      
    } // iC

    // store final projected point in x
    for (i=0; i<dim; i++) x0[i] = xproj[i];

  } // ip      

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_PointsOnGeom
int 
xf_PointsOnGeom( xf_Geom *Geom, int iComp, int dim, int np, real dl, 
		 enum xfe_GeomSpacingType Spacing, int *pnout, real **px)
{
  int ierr;
  xf_GeomComp *Comp;

  // check dimension of geometry for consistency
  if (Geom->Dim != dim) return xf_Error(xf_INPUT_ERROR);
						
  // pointer to geometry component
  Comp = Geom->Comp+iComp;

  switch (Comp->Type){
  case xfe_GeomCompNone:
    break;
  case xfe_GeomCompAnalytical:
    ierr = xf_Error(xf_PointsOnGeom_Analytical( (xf_GeomCompAnalytical *) Comp->Data, 
						dim, np, dl, Spacing, pnout, px));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_GeomCompSpline:
    ierr = xf_Error(xf_PointsOnGeom_Spline( (xf_GeomCompSpline *) Comp->Data, 
					    dim, np, dl, Spacing, pnout, px));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}

