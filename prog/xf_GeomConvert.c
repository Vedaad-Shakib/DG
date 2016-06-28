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
 FILE:  xf_GeomConvert.c
 
 This program converts xflow geometry to various formats
 
 */

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Geom.h"
#include "xf_Mesh.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_Arg.h"
#include "xf_Geom.h"
#include "xf_GeomIO.h"

/* Farfield structure */
typedef struct
{
  int dim;                 // spatial dimension
  char type[xf_MAXSTRLEN]; // string indicating type (e.g. Box)
  real center[3];          // farfield is centered about this point
  real dist;               // distance from center to farfield
  real space;              // mesh/point spacing desired at farfield
}
xf_Farfield;

/******************************************************************/
//   FUNCTION Definition: xf_GeomToBamgPoints
static int 
xf_GeomToBamgPoints( xf_Geom *Geom, real nearspace, 
		     xf_Farfield Farfield, xf_BamgGeomPoints **pBamgPoints)
{
  int ierr;
  int ibfgrp;
  int iComp;
  int i, n, k, kk, nf;
  real x, y, dx, dy, dd;
  real *xv = NULL;
  real IFB[4] = {1.0, 0.0, -1.0, 0.0};
  xf_BamgGeomPoints *BP;
  
  
  // allocate points structure
  ierr = xf_Error(xf_CreateBamgGeomPoints(pBamgPoints));
  if (ierr != xf_OK) return ierr;
  BP = (*pBamgPoints);
  
  // initialize boundary face group counter
  ibfgrp = 0;

  // only box farfield is suported
  if (strcmp(Farfield.type, "Box") != 0) return xf_Error(xf_NOT_SUPPORTED);
  
  // dimension needs to be 2
  if (Farfield.dim != 2) return xf_Error(xf_NOT_SUPPORTED);

  // add farfield points
  nf = 2.*Farfield.dist/Farfield.space+1; // number of points per edge
  
  BP->nVert = 4*(nf-1); // number of vertices
  ierr = xf_Error(xf_Alloc( (void **) &BP->Vert, 2*BP->nVert, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // vertex coordinates (bottom, right, top, left)
  x = Farfield.center[0]-Farfield.dist;
  y = Farfield.center[1]-Farfield.dist;
  dd = 2.*Farfield.dist/(nf-1.);
  for (i=0; i<BP->nVert; i++){
    BP->Vert[2*i+0] = x;
    BP->Vert[2*i+1] = y;
    if (i%(nf-1) == 0){
      k = i/(nf-1);
      dx = dd*IFB[ k   %4];
      dy = dd*IFB[(k+3)%4];
    }
    x = x + dx;
    y = y + dy;
  } // i

  // edges are just one big loop on farfield
  BP->nEdge = 4*(nf-1);
  ierr = xf_Error(xf_Alloc( (void **) &BP->Edge, 3*BP->nEdge, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<BP->nEdge; i++){
    if (i%(nf-1) == 0) ibfgrp++;
    BP->Edge[3*i+0] = i+1;
    BP->Edge[3*i+1] = (i+1)%BP->nEdge + 1;
    BP->Edge[3*i+2] = ibfgrp;
  }


  // add geometry component points
  for (iComp=0; iComp<Geom->nComp; iComp++){

    // get points on component iComp
    ierr = xf_Error(xf_PointsOnGeom(Geom, iComp, 2, -1, nearspace,
				    xfe_GeomSpacingDefaultD, &n, &xv));
    if (ierr != xf_OK) return ierr;

    // add to total number of points
    BP->nVert += n;

    // reallocate vertex coordinates
    ierr = xf_Error(xf_ReAlloc( (void **) &BP->Vert, 2*BP->nVert, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    // copy over xv
    for (k=0; k<2*n; k++) BP->Vert[(BP->nVert-n)*2+k] = xv[k];

    // add to edge counter (assume point list is a loop)
    BP->nEdge += n;
    
    // reallocate edge list
    ierr = xf_Error(xf_ReAlloc( (void **) &BP->Edge, 3*BP->nEdge, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    // increment boundary face group counter
    ibfgrp++;

    // add edges (assume point list is ordered correctly)
    for (k=0; k<n; k++){
      kk = 3*(BP->nEdge-n+k);
      BP->Edge[kk  ] = BP->nVert-n+ k     +1;  // first node of edge
      BP->Edge[kk+1] = BP->nVert-n+(k+1)%n+1;  // second node of edge
      BP->Edge[kk+2] = ibfgrp;           // boundary flag
    }

    // release x
    xf_Release( (void *) xv);
    
  } // iComp
  

  return xf_OK;

}




/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, len;
  int outExtLen;
  char *ArgIn[] = {"in", "NULL", "input geometry file (xflow)",
		   "out", "NULL", "output file (type by extension)",
		   "fartype", "Box", "type of farfield",
		   "farcenter", "0 0", "center of farfield [dim]",
		   "fardist", "50", "farfield characteristic distance",
		   "farspace", "10", "farfield spacing",
		   "nearspace", "0.1", "nearfield spacing",
		   "\0"};
  char inFile[xf_MAXSTRLEN];  
  char outFile[xf_MAXSTRLEN];
  char value[xf_MAXSTRLEN];
  char *outExt;
  real nearspace;
  xf_KeyValue KeyValue;
  xf_Geom *Geom;
  xf_Farfield Farfield;
  xf_BamgGeomPoints *BamgPoints;
  
  xf_printf("\n");
  xf_printf("=== xf_GeomConvert: Geometry Conversion  ===\n");
  xf_printf("\n");
  
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;
  
  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
  
  // Get inFile
  ierr = xf_GetKeyValue(KeyValue, "in", inFile);
  if (ierr != xf_OK) return ierr;
  
  // Get outFile
  ierr = xf_GetKeyValue(KeyValue, "out", outFile);
  if (ierr != xf_OK) return ierr;
  
  // Get fartype
  ierr = xf_GetKeyValue(KeyValue, "fartype", Farfield.type);
  if (ierr != xf_OK) return ierr;

  // Get farcenter
  ierr = xf_GetKeyValue(KeyValue, "farcenter", value);
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ScanXReal(value, &Farfield.dim, Farfield.center));
  if (ierr != xf_OK) return ierr;

  // Get fardist
  ierr = xf_GetKeyValueReal(KeyValue, "fardist", &Farfield.dist);
  if (ierr != xf_OK) return ierr;

  // Get farspace
  ierr = xf_GetKeyValueReal(KeyValue, "farspace", &Farfield.space);
  if (ierr != xf_OK) return ierr;

  // Get nearspace
  ierr = xf_GetKeyValueReal(KeyValue, "nearspace", &nearspace);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;

  // extension is required on outFile
  if ((len = strlen(outFile)) <= 0) return xf_Error(xf_FILE_READ_ERROR);
  for (outExtLen=0; ((outExtLen<len)&&(outFile[len-outExtLen-1]!='.')); outExtLen++);
  outExt = outFile + len - outExtLen; 
  
  // Create a geometry
  ierr = xf_Error(xf_CreateGeom(&Geom));
  if (ierr != xf_OK) return ierr;

  // Read the input geometry
  ierr = xf_Error(xf_ReadGeomFile(inFile, NULL, Geom));
  if (ierr != xf_OK) return ierr;

  // Write output file based on extension
  if ((outExtLen == 4) && (strncmp( outExt, "bamg", 4) == 0)){  // BAMG format

    // convert geometry + farfield to points
    ierr = xf_Error(xf_GeomToBamgPoints(Geom, nearspace, 
					Farfield, &BamgPoints));
    if (ierr != xf_OK) return ierr;
    
    // write BAMG points
    ierr = xf_Error( xf_WriteBamgGeometry_Points(BamgPoints, outFile) );
    if (ierr != xf_OK) return ierr;

    // destroy BAMG points
    ierr = xf_Error(xf_DestroyBamgGeomPoints(BamgPoints));
    if (ierr != xf_OK) return ierr;

  }
  else
    return xf_Error(xf_NOT_SUPPORTED);


  // destroy geometry
  ierr = xf_Error(xf_DestroyGeom(Geom));
  if (ierr != xf_OK) return ierr;

  
  xf_printf("xf_GeomConvert finished.\n");
  
  return xf_OK;
}
