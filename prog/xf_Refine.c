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
  FILE:  xf_Refine.c

  This program refines a mesh.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Param.h"
#include "xf_Mesh.h"
#include "xf_Math.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Arg.h"
#include "xf_AdaptHang.h"
#include "xf_MeshTools.h"
#include "xf_Solver.h"
#include "xf_Residual.h"
#include "xf_SolverTools.h"
#include "xf_Adapt.h"
#include "xf_ParamDefault.h"
#include "xf_Geom.h"
#include "xf_GeomIO.h"


/******************************************************************/
//   FUNCTION Definition: RefineTriTet
static int 
xf_RefineTriTet(xf_All *All, const char *outFile)
{
  int ierr, i, d, dim, nelemtot;
  int ibfgrp, ibface, egrp, elem, face;
  int nfnode, fvec[xf_MAXQ1FACENODE];
  int q, qq, fac, ix, iy;
  int *Node;
  enum xfe_BasisType QBasis;
  FILE *fgri = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim  = Mesh->Dim;

  ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;

  // open .gri file and write to it
  if ((fgri = fopen(outFile ,"w"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);

  fac = ((dim == 2) ? 2 : 5);
  fprintf(fgri, "%d %d %d\n", Mesh->nNode, fac*nelemtot, dim);
		  
  /* Write Node coordinates (these stay the same)*/
  for (i=0; i<Mesh->nNode; i++){
    for (d=0; d<dim; d++) fprintf(fgri, "%.15E ", Mesh->Coord[i][d]);
    fprintf(fgri, "\n");
  }
  
  // bfacegroups
  fprintf(fgri, "%d\n", Mesh->nBFaceGroup);
  fac = ((dim == 3) ? 2 : 1);
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    if (Mesh->BFaceGroup[ibfgrp].nBFace > 0){
      egrp = Mesh->BFaceGroup[ibfgrp].BFace[0].ElemGroup;
      elem = Mesh->BFaceGroup[ibfgrp].BFace[0].Elem;
      face = Mesh->BFaceGroup[ibfgrp].BFace[0].Face;
      // local nodes on face
      ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, Mesh->ElemGroup[egrp].QOrder, 
				       face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;
    }
    else nfnode = 0;
    fprintf(fgri, "%d %d %s\n", fac*Mesh->BFaceGroup[ibfgrp].nBFace, nfnode, Mesh->BFaceGroup[ibfgrp].Title);
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      egrp = Mesh->BFaceGroup[ibfgrp].BFace[ibface].ElemGroup;
      elem = Mesh->BFaceGroup[ibfgrp].BFace[ibface].Elem;
      face = Mesh->BFaceGroup[ibfgrp].BFace[ibface].Face;
      // local nodes on face
      ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, Mesh->ElemGroup[egrp].QOrder, 
				       face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;
      Node = Mesh->ElemGroup[egrp].Node[elem];
      if (dim == 2) // tri edges same as tet edges
	for (i=0; i<nfnode; i++)
	  fprintf(fgri, "%d ", Node[fvec[i]]+1);
      else{
	if (nfnode != 4) return xf_Error(xf_INPUT_ERROR);
	fprintf(fgri, "%d %d %d\n", Node[fvec[0]]+1, Node[fvec[2]]+1, Node[fvec[3]]+1);
	fprintf(fgri, "%d %d %d\n", Node[fvec[3]]+1, Node[fvec[1]]+1, Node[fvec[0]]+1);
      }
      fprintf(fgri, "\n");
    }
  } // ibfgrp

  // element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    QBasis = ((dim == 2) ? xfe_TriLagrange : xfe_TetLagrange);
    fac = ((dim == 2) ? 2 : 5);
    q = Mesh->ElemGroup[egrp].QOrder;
    qq = q+1;
    fprintf(fgri, "%d %d %s\n", fac*Mesh->ElemGroup[egrp].nElem, q, xfe_BasisName[QBasis]);
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      Node = Mesh->ElemGroup[egrp].Node[elem];
      if (dim == 2){
	for (iy=0; iy<qq; iy++) // elem 0
	  for (ix=0; ix<(qq-iy); ix++)
	    fprintf(fgri, "%d ", Node[iy*qq+ix]+1);
	fprintf(fgri, "\n");
	for (ix=qq-1; ix>=0; ix--) // elem 1
	  for (iy=qq-1-ix; iy<qq; iy++) 
	    fprintf(fgri, "%d ", Node[iy*qq+ix]+1);
	fprintf(fgri, "\n");
      }
      else{
	return xf_Error(xf_NOT_SUPPORTED);
      }
    }
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_SetBlayerIndicator
static int 
xf_SetBlayerIndicator(xf_All *All, const char *bgroups, int nblayer, 
		      xf_Vector *RefIndicator)
{
  int ierr;
  int nbgroup;
  int ibfgrp, ibg, ibface, egrp, elem;
  int ilayer, iiface, egrpL, elemL, egrpR, elemR;
  int ivalL, ivalR, ival;
  xf_IFace IFace;
  enum xfe_Bool found = xfe_False;
  char **bgroup = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // set indicator to all 0s
  ierr = xf_Error(xf_SetConstVector(RefIndicator, 0, 0));
  if (ierr != xf_OK) return ierr;

  // determine which boundaries we're looking at
  ierr = xf_Error(xf_ScanXStringAlloc(bgroups, xf_MAXSTRLEN, &nbgroup, &bgroup));
  if (ierr != xf_OK) return ierr;

  // loop over boundaries, flag elements next to boundaries
  for (ibg=0; ibg<nbgroup; ibg++){
    found = xfe_False;
    for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
      if (strcmp(Mesh->BFaceGroup[ibfgrp].Title, bgroup[ibg]) == 0){
	found = xfe_True;
	break;
      }
    }
    if (!found){
      xf_printf("Mesh does not have a boundary group entitiled %s\n", bgroup[ibg]);
      return xf_Error(xf_INPUT_ERROR);
    }
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      egrp = Mesh->BFaceGroup[ibfgrp].BFace[ibface].ElemGroup;
      elem = Mesh->BFaceGroup[ibfgrp].BFace[ibface].Elem;
      RefIndicator->GenArray[egrp].iValue[elem][0] = 1;
    } // ibface 
  }

  // flag elements on additional layers
  for (ilayer=1; ilayer<nblayer; ilayer++){
    for (iiface=0; iiface<Mesh->nIFace; iiface++){
      IFace = Mesh->IFace[iiface];
      egrpL = IFace.ElemGroupL;
      egrpR = IFace.ElemGroupR;
      elemL = IFace.ElemL;
      elemR = IFace.ElemR;
      ivalL = RefIndicator->GenArray[egrpL].iValue[elemL][0];
      ivalR = RefIndicator->GenArray[egrpR].iValue[elemR][0];
      if ((ivalL == 1) && (ivalR == 0))
	RefIndicator->GenArray[egrpR].iValue[elemR][0] = -1;
      else if ((ivalL == 0) && (ivalR == 1))
	RefIndicator->GenArray[egrpR].iValue[elemL][0] = -1;
    }
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
	if (RefIndicator->GenArray[egrp].iValue[elem][0] == -1)
	  RefIndicator->GenArray[egrp].iValue[elem][0] = 1;
  }
  

  xf_Release2( (void **) bgroup);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteRefined
static int 
xf_WriteRefined(xf_All *All, char *outFile)
{
  // Writes All in .gri or .xfa format
  int ierr;
  int len;
  char *outExt;

  if ((len = strlen(outFile)) < 4) return xf_Error(xf_FILE_WRITE_ERROR);
  outExt = outFile + len - 4; // pointer to extension
  
  if (strncmp(outExt, ".gri", 4) == 0){
    ierr = xf_Error(xf_WriteGriFile(All->Mesh, outFile));
    if (ierr != xf_OK) return ierr;
  }
  else if (strncmp(outExt, ".xfa", 4) == 0){
    ierr = xf_Error(xf_WriteAllBinary(All, outFile));
    if (ierr!=xf_OK) return ierr;
  }
  else return xf_Error(xf_NOT_SUPPORTED);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ChooseBFGTitles
static int 
xf_ChooseBFGTitles(xf_Mesh *Mesh, const char *Question, int *pnbf, 
		   char *BFGTitles)
{
/*
 PURPOSE:
 
   Queries user for list of boundary titles.  Returns space-separated
   list.
 
 INPUTS:
 
   Mesh : mesh structure
   Question : prompt question to ask
 
 OUTPUTS:
   
   (*pnbf) : number of boundary face groups selected
   BFGTitles : space-separated list of boundary-face group titles

 RETURN:  Error Code

*/
  int ierr, i;
  int nbf, *vbf = NULL;
  int ibfgrp;
  char instr[xf_MAXSTRLEN];


  // get list of wall boundaries
  xf_printf("The following boundary groups are available:\n\n");
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++)
    xf_printf(" %d : %s\n", ibfgrp, Mesh->BFaceGroup[ibfgrp].Title);
  xf_printf("\n%s\n", Question);
  xf_printf("Provide a space-separated list of numbers:\n> ");
  if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) return xf_INPUT_ERROR;

  // convert list of numbers to list of names
  ierr = xf_Error(xf_ScanXIntAlloc(instr, &nbf, &vbf));
  if (ierr != xf_OK) return ierr; 
  if (nbf <= 0) return xf_INPUT_ERROR;
  sprintf(BFGTitles, "\0");
  for (i=0; i<nbf; i++){
    sprintf(instr, "%s %s", BFGTitles, Mesh->BFaceGroup[vbf[i]].Title);
    strcpy(BFGTitles, instr);
  }
  if (strlen(BFGTitles) == 0) return xf_INPUT_ERROR;
 
 xf_Release( (void *) vbf);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WallDistanceMetric
static int 
xf_WallDistanceMetric(xf_All *All, real **pNM)
{
/*
 PURPOSE:
 
   Computes a node-based metric using the wall distance function
 
 INPUTS:
 
   All : all structure
 
 OUTPUTS: None, All->Mesh is modified

 RETURN:  Error Code
*/
  int ierr, i, dim, dim2;
  int DistFcnOrder;
  int ibfgrp, egrp, elem, j, nnode;
  int *counter = NULL;
  int choice, node, k;
  int nvec[xf_MAXQ1NODE];
  enum xfe_Bool done;
  enum {LinearMetric, LogMetric, PowerMetric};
  int MetricType;
  char instr[xf_MAXSTRLEN];
  char DistFcnWallBoundariesOrig[xf_MAXSTRLEN];
  char DistFcnWallBoundaries[xf_MAXSTRLEN];
  real *NM;
  real d, d2, h, h0, m, p;
  xf_Vector *WD;
  xf_Mesh *Mesh;
  xf_ElemGroup *EG;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  dim2 = dim*dim;

  xf_printf("\nSetting metric based on a wall distance function.\n");

  // set parameters for calculating wall distance
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "DistFcnOrder", &DistFcnOrder));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "DistFcnOrder", 1)); // harcoded
  if (ierr != xf_OK) return ierr;

  // get list of wall boundaries
  ierr = xf_Error(xf_ChooseBFGTitles(Mesh, "Which of these should be treated as walls?", 
				     NULL, DistFcnWallBoundaries));
  if (ierr != xf_OK) return ierr;

  // set this list (one string) as a parameter, read by CalculateDistFcn
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "DistFcnWallBoundaries", DistFcnWallBoundariesOrig));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "DistFcnWallBoundaries", DistFcnWallBoundaries));
  if (ierr != xf_OK) return ierr;

  // calculate wall distance
  ierr = xf_Error(xf_FindSupportedVector(All, "WallDistance", &WD));
  if (ierr != xf_OK) return ierr;

  // request functional form for metric
  done = xfe_False;
  while (!done){

    xf_printf("Choose a form for mesh size as a function of d = distance to wall:\n\n");
    xf_printf("   0  Linear: h(d) = h0 + m*d\n");
    xf_printf("   1  Log:    h(d) = h0 * [1 + log(1+m*d)]\n");
    xf_printf("   2  Power:  h(d) = h0 * [1 + m*d^p]\n");
    // TODO: add a mask capability
    xf_printf("\n> ");
    if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) continue;
    sscanf(instr, "%d", &choice);

    MetricType = choice;

    if (MetricType == LinearMetric){
      xf_printf("\nEnter h0 and m (space separated)\n");
      xf_printf("\n> ");
      if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) continue;
      if (sscanf(instr, "%lf %lf", &h0, &m) == 2) done = xfe_True;
      MetricType = LinearMetric;
    }
    else if (MetricType == LogMetric){
      xf_printf("\nEnter h0=spacing at wall, and d2=dist for doubling of h (space separated)\n");
      xf_printf("\n> ");
      if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) continue;
      if (sscanf(instr, "%lf %lf", &h0, &d2) == 2) done = xfe_True;
      MetricType = LogMetric;
      m = (exp(1)-1)/d2;
    }
    else if (MetricType == PowerMetric){
      xf_printf("\nEnter h0=spacing at wall (eg .01), p=power (eg 1.3), and d2=dist for doubling of h (eg 0.1) (space separated)\n");
      xf_printf("\n> ");
      if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) continue;
      if (sscanf(instr, "%lf %lf %lf", &h0, &p, &d2) == 3) done = xfe_True;
      MetricType = PowerMetric;
      m = 1./pow(d2,p);
    }
    else return xf_INPUT_ERROR;
  }

  // allocate (*pNM)
  ierr = xf_Error(xf_Alloc( (void **) pNM, Mesh->nNode*dim2, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  NM = (*pNM);
  for (j=0; j<Mesh->nNode*dim2; j++) NM[j] = 0.;

  // allocate a counter vector
  ierr = xf_Error(xf_Alloc( (void **) &counter, Mesh->nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (j=0; j<Mesh->nNode; j++) counter[j] = 0;


  // loop over elements (effectively nodes) and calculate metric
  for (egrp=0; egrp < Mesh->nElemGroup; egrp++){
    for (elem=0, EG=Mesh->ElemGroup+egrp; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      ierr = xf_Error(xf_Q1Nodes(EG->QBasis, EG->QOrder, &nnode, nvec));
      if (ierr != xf_OK) return ierr;
      for (j=0; j<nnode; j++){
	node = Mesh->ElemGroup[egrp].Node[elem][nvec[j]];
	if (counter[node] == 1) continue;
	// We assume that the distance function is approximated by a uniform-Lagrange basis
	d = WD->GenArray[egrp].rValue[elem][j];
	switch(MetricType){
	case LinearMetric:
	  h = h0 + m*d;
	  break;
	case LogMetric:
	  h = h0*(1+log(1 + m*d));
	  break;
	case PowerMetric:
	  h = h0*(1+m*pow(d,p));
	  break;
	default:
	  return xf_Error(xf_NOT_SUPPORTED);
	  break;
	}
	for (k=0; k<dim2; k+=(dim+1)) NM[node*dim2+k] = h;	
	counter[node] = 1;
      } // j

    } // elem
  } // egrp

  // reset parameters
  ierr = xf_Error(xf_SetKeyValueInt(All->Param->KeyValue, "DistFcnOrder", DistFcnOrder)); 
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "DistFcnWallBoundaries", DistFcnWallBoundariesOrig));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) counter);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ModifyMetricUsingPrimitive
static int 
xf_ModifyMetricUsingPrimitive(xf_All *All, xf_Vector *EM)
{
/*
 PURPOSE:
 
   Modifies metric using a scaling factor in the vicinity of a primitive
 
 INPUTS:
 
   All : all structure
 
 OUTPUTS: None, All->Mesh is modified

 RETURN:  Error Code
*/
  int ierr, dim, dim2;
  int egrp, elem, j, k;
  enum {PrimitiveCircle, PrimitiveDisk, PrimitivePoint};
  int Primitive;
  enum xfe_ShapeType Shape;
  char instr[xf_MAXSTRLEN];
  real xc, yc, rc, r, d;
  real factor, sigma, fac;
  real xref[3], xglob[3];
  xf_ElemGroup *EG;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  dim2 = dim*dim;

  xf_printf("Choose a primitive object to work with:\n");
  xf_printf("   0  Circle centered at xc,yc, and radius r\n");
  xf_printf("   1  Disk centered at xc,yc, and radius r\n");
  xf_printf("   2  Point at xc,yc\n");
  xf_printf("> ");
  if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) return xf_OK;
  if (sscanf(instr, "%d", (int *) &Primitive) != 1) return xf_OK;

  switch (Primitive){
  case PrimitiveCircle:
    xf_printf("Enter center (xc,yc) and radius (rc), three numbers space-separated:\n");
    xf_printf("> ");
    if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) return xf_OK;
    if (sscanf(instr, "%lf %lf %lf", &xc, &yc, &rc) != 3) return xf_OK;  
    break;
  case PrimitiveDisk:
    xf_printf("Enter center (xc,yc) and radius (rc), three numbers space-separated:\n");
    xf_printf("> ");
    if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) return xf_OK;
    if (sscanf(instr, "%lf %lf %lf", &xc, &yc, &rc) != 3) return xf_OK;  
    break;
  case PrimitivePoint:
    xf_printf("Enter point location (xc,yc), two numbers space-separated:\n");
    xf_printf("> ");
    if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) return xf_OK;
    if (sscanf(instr, "%lf %lf", &xc, &yc) != 2) return xf_OK;  
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  xf_printf("Enter local scaling factor (<1 = refinement) and spread (a length):\n");
  xf_printf("> ");
  if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) return xf_OK;
  if (sscanf(instr, "%lf %lf", &factor, &sigma) != 2) return xf_OK;  

  // loop over elements (effectively nodes) and calculate metric
  for (egrp=0; egrp < Mesh->nElemGroup; egrp++){
    for (elem=0, EG=Mesh->ElemGroup+egrp; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      // determine element shape
      ierr = xf_Error(xf_Basis2Shape(EG->QBasis, &Shape));
      if (ierr != xf_OK) return ierr;
      // shape centroid in ref coordinates
      ierr = xf_Error(xf_ShapeCentroid(Shape, xref));
      if (ierr != xf_OK) return ierr;
      // corresponding global coordinates
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True, 1, xref, xglob));
      if (ierr != xf_OK) return ierr;

      // determine distance to primitive shape
      switch (Primitive){
      case PrimitiveCircle:
	r = sqrt((xglob[0]-xc)*(xglob[0]-xc) + (xglob[1]-yc)*(xglob[1]-yc));
	d = fabs(r-rc);
	break;
      case PrimitiveDisk:
	r = sqrt((xglob[0]-xc)*(xglob[0]-xc) + (xglob[1]-yc)*(xglob[1]-yc));
	d = max(r-rc, 0.);
	break;
      case PrimitivePoint:
	r = sqrt((xglob[0]-xc)*(xglob[0]-xc) + (xglob[1]-yc)*(xglob[1]-yc));
	d = r;
	break;
      default:
	return xf_Error(xf_NOT_SUPPORTED);
	break;
      }

      // determine local factor
      fac = 1. + (factor-1.)*exp(-(d/sigma)*(d/sigma));
      
      // multiply elemental metric by fac
      for (k=0; k<dim2; k++) EM->GenArray[egrp].rValue[elem][k] *= fac;
    
    } // elem
  } // egrp

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_InteractiveRefine
static int 
xf_InteractiveRefine(xf_All *All, char *outFile)
{
/*
 PURPOSE:
 
   Refines a mesh interactively
 
 INPUTS:
 
   All : all structure
 
 OUTPUTS: None, All->Mesh is modified

 RETURN:  Error Code
*/

  int ierr;
  int nelemtot;
  int choice, Q; 
  char instr[xf_MAXSTRLEN];
  char BFGTitles[xf_MAXSTRLEN];
  enum xfe_Bool done;
  enum xfe_Bool IsotropicHMetric;
  enum xfe_Bool UseHydraulicDiam;
  real factor;
  real *NM = NULL;
  xf_Vector *EM;
  xf_Mesh *Mesh, *OldMesh = NULL;


  xf_printf("\n <<< Interactive refinement >>> \n");
  
  done = xfe_False;
  while (!done){

    Mesh = All->Mesh;

    // total number of elements
    ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
    if (ierr != xf_OK) return ierr;

    xf_printf("\nCurrent number of elements = %d.\n", nelemtot);
    xf_printf("\nCurrent number of boundary face groups = %d.\n", Mesh->nBFaceGroup);

    xf_printf("Choose option (Ctrl-c to break):\n\n");
    xf_printf("   0  Scale mesh isotropically (one scaling factor)\n");
    xf_printf("   1  Scale mesh anisotropically (one scaling factor)\n");
    xf_printf("   2  Set mesh size using wall distance\n");
    xf_printf("   3  Scale mesh in vicinity of a geometrical object (e.g. circle)\n");
    xf_printf("   4  Curve a high-order boundary to geometry\n");
    // TODO: isotropic ref with separate boundary factors
    xf_printf("  -1  Write output file\n");
    xf_printf("  -2  Write output file and exit\n");
    xf_printf("  -3  Undo last operation (restore previous mesh)\n");
    xf_printf("\n> ");
    if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) continue;
    sscanf(instr, "%d", &choice);

    if (choice == -3){ // Undo
      if (OldMesh == NULL){
	xf_printf("No previous mesh (can only undo one step).\n");
	continue;
      }
      xf_printf("Undoing last step -- restoring previous mesh.\n");
      // destroy existing mesh
      ierr = xf_Error(xf_DestroyMesh(Mesh));
      if (ierr != xf_OK) return ierr;
      // set pointer to previous mesh
      Mesh = All->Mesh = OldMesh;
      OldMesh = NULL;
      continue;
    }


    if ((choice == 0) || (choice == 1) || (choice == 3)){
      xf_printf("\nEnter global h scaling factor (less than 1 indicates refinement).\n");
      xf_printf("> ");
      if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) continue;
      sscanf(instr, "%lf", &factor);
      
      // choice 0 -> isotropic;  choice 1 -> anisotropic
      ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "IsotropicHMetric", 
				       &IsotropicHMetric));
      ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "IsotropicHMetric", 
					 (choice != 1)));
      if (ierr != xf_OK) return ierr;
      // We want to use the hhydraulic diameter
      ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "IsotropicHMetricIsHD", 
					 &UseHydraulicDiam));
      ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "IsotropicHMetricIsHD", 
					 xfe_True));
      if (ierr != xf_OK) return ierr;
      
      // compute current metric
      ierr = xf_Error(xf_FindElemHMetric(All, xfe_False, &EM));
      if (ierr != xf_OK) return ierr;

      if (choice == 3){
	ierr = xf_Error(xf_ModifyMetricUsingPrimitive(All, EM));
	if (ierr != xf_OK)
	  xf_printf("Error modifying metric; continuing.\n");
	  
      }
      /* multiply metric by factor.  Note, 6/sqrt(3) is exact for
	 equilateral triangles when using the hydraulic diameter. */
      ierr = xf_Error(xf_VectorMult(EM, factor*6.0/sqrt(3)));
      if (ierr != xf_OK) return ierr;
      
      // destroy old mesh
      ierr = xf_Error(xf_DestroyMesh(OldMesh));
      if (ierr != xf_OK) return ierr;

      // only have one tool for refining a mesh ... Bamg
      ierr = xf_Error(xf_RefineMeshBamg(All, (choice == 1), EM, NULL, &OldMesh));
      if (ierr != xf_OK) return ierr;
      Mesh = All->Mesh;

      // reset IsotropicHMetric flag and hydraulic diameter flag
      ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "IsotropicHMetric", 
					 IsotropicHMetric));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_SetKeyValueBool(All->Param->KeyValue, "IsotropicHMetricIsHD", 
					 UseHydraulicDiam));
      if (ierr != xf_OK) return ierr;

    }
    else if (choice == 2){

      // Compute a wall-distance-based metric
      ierr = xf_Error(xf_WallDistanceMetric(All, &NM));
      if (ierr != xf_OK) continue; // ignore errors

      // destroy old mesh
      ierr = xf_Error(xf_DestroyMesh(OldMesh));
      if (ierr != xf_OK) return ierr;

      // only have one tool for refining a mesh ... Bamg
      ierr = xf_Error(xf_RefineMeshBamg(All, xfe_False, NULL, NM, &OldMesh));
      if (ierr != xf_OK) return ierr;
      Mesh = All->Mesh;

    }   
    else if (choice == 4){

      // make sure we have geometry
      if ((All->Geom == NULL) || (All->Geom->nComp <= 0)){
	xf_printf("No geometry.\n");
	continue;
      }

      // ask for BFGTitles for curving
      ierr = xf_Error(xf_ChooseBFGTitles(Mesh, "Which of these should be used for curving?", 
					 NULL, BFGTitles));
      if (ierr != xf_OK) return ierr;
     
      // ask for desired Q
      xf_printf("What geometry order (Q) should the adjacent elements be?\n");
      xf_printf("\n> ");
      if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) continue;
      if (sscanf(instr, "%d", &Q) != 1) continue;

      // curve mesh
      ierr = xf_Error(xf_CurveMeshBoundary(All->Mesh, All->Geom, BFGTitles, Q));
      if (ierr != xf_OK) return ierr;

    }  
    else if ((choice == -1) || (choice == -2)){
      // always request to input an output file
      xf_printf("Provide output file name (include .gri or .xfa extension)\n> ");
      if (fgets(instr, xf_MAXSTRLEN, stdin) == NULL) continue;
      sscanf(instr, "%s", outFile);      
      if (strlen(outFile) < 4) continue;
      xf_printf("\nWriting %s.\n", outFile);
      ierr = xf_Error(xf_WriteRefined(All, outFile));
      if (ierr != xf_OK) continue;
      if (choice == -2){
	xf_printf("\nExiting.\n");
	done = xfe_True;
      }
    }
    else{
      xf_printf("Unrecognized choice.  Please try again.\n");
    }

  } // end while !done  

  // clean up
  xf_Release( (void *) NM);

  return xf_OK;

}

/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, len, nblayer;
  char *ArgIn[] = {"in", "NULL", "input file (.xfa or .gri)",
		   "out", "NULL", "output file (.xfa or .gri)",
		   "geom", "NULL", "optional geometry file",
		   "uniform", "False", "uniform refinement",
		   "bgroups", "None", "names of bgroups for proximity ref",
		   "nblayer", "1", "# of elem layers to refine close to broups",
		   "tritet", "False", "True to convert quad to tri, hex to tet",
		   "interact", "False", "True to enter interactive mode",
                   "orderref", "False", "True to do order refinement",
		   "\0"};
  char inFile[xf_MAXSTRLEN];  
  char outFile[xf_MAXSTRLEN];
  char bgroups[xf_MAXSTRLEN];
  char geomFile[xf_MAXSTRLEN];  
  char *inExt, *outExt;
  enum xfe_Bool uniform;
  enum xfe_Bool tritet;
  enum xfe_Bool interact;
  enum xfe_Bool orderref;
  enum xfe_Bool written = xfe_False;
  xf_Vector *RefIndicator;
  xf_TimeHistData *TimeHistData = NULL;
  xf_KeyValue KeyValue;
  xf_All *All;

  xf_printf("\n");
  xf_printf("=== xf_Refine: Mesh Refinement  ===\n");
  xf_printf("\n");
    
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  xf_printf("nKey = %d\n", KeyValue.nKey);
  for (i=0; i<KeyValue.nKey; i++)
    xf_printf("%d : Key = %s, Value = %s\n", i, KeyValue.Key[i], KeyValue.Value[i]);
    
  // Get inFile
  ierr = xf_GetKeyValue(KeyValue, "in", inFile);
  if (ierr != xf_OK) return ierr;

  // Get outFile
  ierr = xf_GetKeyValue(KeyValue, "out", outFile);
  if (ierr != xf_OK) return ierr;

  // Get uniform
  ierr = xf_GetKeyValueBool(KeyValue, "uniform", &uniform);
  if (ierr != xf_OK) return ierr;

  // Get tritet
  ierr = xf_GetKeyValueBool(KeyValue, "tritet", &tritet);
  if (ierr != xf_OK) return ierr;

  // Get interact
  ierr = xf_GetKeyValueBool(KeyValue, "interact", &interact);
  if (ierr != xf_OK) return ierr;

  // Get orderref
  ierr = xf_GetKeyValueBool(KeyValue, "orderref", &orderref);
  if (ierr != xf_OK) return ierr;

  // Get bgroups
  ierr = xf_GetKeyValue(KeyValue, "bgroups", bgroups);
  if (ierr != xf_OK) return ierr;

  // Get nblayer
  ierr = xf_GetKeyValueInt(KeyValue, "nblayer", &nblayer);
  if (ierr != xf_OK) return ierr;

  // Get geomFile
  ierr = xf_GetKeyValue(KeyValue, "geom", geomFile);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  // extensions are required
  if ((len = strlen(inFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
  inExt = inFile + len - 4; // pointer to extension

  if ((len = strlen(outFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
  outExt = outFile + len - 4; // pointer to extension


  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  /* Read input file */
  if (strncmp( inExt, ".xfa", 4) == 0){
    // read .xfa
    ierr = xf_Error(xf_ReadAllBinary(inFile, All));
    if (ierr!=xf_OK) return ierr;
  }
  else if (strncmp( inExt, ".gri", 4) == 0){
    // read .gri
    ierr = xf_Error(xf_ReadGriFile(inFile, NULL, All->Mesh));
    if (ierr!=xf_OK) return ierr;
    // set default params
    ierr = xf_AddKeyValueList(&All->Param->KeyValue, xf_DefaultParamList,
			      xfe_True, xfe_False);
    if ((ierr != xf_OK) && (ierr != xf_NOT_FOUND)) return xf_Error(ierr);
  }
  else if (strncmp( inExt, ".txt", 4) == 0){
    // read input time history
    ierr = xf_Error(xf_ReadTimeHistData(inFile, NULL, &TimeHistData));
    if (ierr != xf_OK) return ierr;
  }
  else return xf_Error(xf_NOT_SUPPORTED);

  // read geom file if specified
  if (xf_NotNull(geomFile)){
    xf_printf("Reading geometry file.\n");
    ierr = xf_Error(xf_ReadGeomFile(geomFile, NULL, All->Geom));
    if (ierr != xf_OK) return ierr;
  }


  if (TimeHistData != NULL){

    if (uniform){
      /* Perform the temporal adaptation */
      ierr = xf_Error(xf_SplitTimeHistData(TimeHistData, NULL, -1));
      if (ierr != xf_OK) return ierr;
    }

    // order refinement means DG1->DG2
    if (orderref){
      for (i=0; i<TimeHistData->nTime; i++)
        if (TimeHistData->TimeScheme[i] == xfe_TimeSchemeDG1)
          TimeHistData->TimeScheme[i] = xfe_TimeSchemeDG2;
    }
  }
  else if (uniform){ // uniform refinement

    // create a refinement indicator
    ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				  NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_False,  xfe_True, NULL, 
				  &RefIndicator, NULL));
    if (ierr != xf_OK) return ierr;
    
    // set indicator to all 1s
    ierr = xf_Error(xf_SetConstVector(RefIndicator, 1, 0));
    if (ierr != xf_OK) return ierr;

    /* Hanging-node adaptation (uniform) */
    ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
    if (ierr != xf_OK) return ierr;

  }
  else if (xf_NotNull(bgroups)){ // bgroup refinement
    // create a refinement indicator
    ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				  NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_False,  xfe_True, NULL, 
				  &RefIndicator, NULL));
    if (ierr != xf_OK) return ierr;
    
    // set indicators next to desired boundaries to 1
    ierr = xf_Error(xf_SetBlayerIndicator(All, bgroups, nblayer, RefIndicator));
    if (ierr != xf_OK) return ierr;

    /* Hanging-node adaptation (uniform) */
    ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
    if (ierr != xf_OK) return ierr;

  }
  else if (tritet){  // quad to tri or hex to tet

    if (strncmp(outExt, ".gri", 4) != 0)
      return xf_Error(xf_NOT_SUPPORTED);
    
    ierr = xf_Error(xf_RefineTriTet(All, outFile));
    if (ierr != xf_OK) return ierr;

    written = xfe_True;

  }
  else if (interact){
    
    ierr = xf_Error(xf_InteractiveRefine(All, outFile));
    if (ierr != xf_OK) return ierr;
    if ((len = strlen(outFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
    outExt = outFile + len - 4; // pointer to extension
    
  }
  else return xf_Error(xf_INPUT_ERROR);

  
  /* Write output file */
  if (TimeHistData != NULL){
    ierr = xf_Error(xf_WriteTimeHist(TimeHistData, outFile));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
    if (ierr != xf_OK) return ierr;
    xf_printf("Wrote %s\n", outFile);
    written = xfe_True;
  }
  if (!written){
    ierr = xf_Error(xf_WriteRefined(All, outFile));
    if (ierr != xf_OK) return ierr;
  }


  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;
  

  xf_printf("xf_Refine finished.\n");

  return xf_OK;
}
