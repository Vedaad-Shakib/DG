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
  FILE:  xf_GeomPanel.c

  This file contains functions for working with panel geometry
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
#include "xf_QuadRule.h"
#include "xf_Basis.h"
#include "xf_MeshDistance.h"



/******************************************************************/
//   FUNCTION Definition: xf_ProjectToGeomComp_Panel
int 
xf_ProjectToGeomComp_Panel( xf_GeomCompPanel *GCP, int dim, real *x)
{
  int ierr;
  int nn, ip, j, k;
  real *X = NULL;
  real dist, distmin;
  real x0[3], xproj[3];

  // determine number of nodes per panel
  ierr = xf_Error(xf_Order2nNode(GCP->Basis, GCP->Order, &nn));
  if (ierr != xf_OK) return ierr;

  // allocate memory for passing node coordinates
  ierr = xf_Error(xf_Alloc((void **) &X, nn*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // copy x to x0
  for (k=0; k<dim; k++) x0[k] = x[k];

  // loop over panels
  distmin = 1e30;
  for (ip=0; ip<GCP->nPanel; ip++){

    // store node coords in a local structure
    for (j=0; j<nn; j++)
      for (k=0; k<dim; k++) X[j*dim+k] = GCP->Coord[GCP->Panels[ip][j]][k];

    // project x to panel
    ierr = xf_Error(xf_Dist2Face(dim, x0, GCP->Basis, GCP->Order, nn, X, &dist, xproj));
    if (ierr != xf_OK) return ierr;

    // check for minimum
    if (dist < distmin){
      for (k=0; k<dim; k++) x[k] = xproj[k]; // store minimum projection in x
      distmin = dist;
    }

  } // ip
  
  // release memory
  xf_Release( (void *) X);


  return xf_OK;
}


