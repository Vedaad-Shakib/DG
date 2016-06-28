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
  FILE:  xf_AdaptToVol.c

  This program adapts a mesh to a volume requirement (from another mesh)

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Param.h"
#include "xf_Mesh.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Arg.h"
#include "xf_AdaptHang.h"
#include "xf_MeshTools.h"
#include "xf_Solver.h"
#include "xf_SolverTools.h"
#include "xf_Adapt.h"
#include "xf_MPI.h"
#include "xf_Math.h"



/******************************************************************/
//   FUNCTION Definition: CreateVolList
static int 
xf_CreateVolList(xf_All *All, const char *vollist)
{
  // writes a volume list file
  int ierr;
  int egrp, elem, d, dim;
  enum xfe_Bool PointsChanged = xfe_False;
  enum xfe_BasisType QBasis;
  enum xfe_ShapeType Shape;
  FILE *fid;
  real xref[3], xglob[3], ElemVol;
  xf_BasisData *GeomPhiData = NULL;
  xf_Vector *EG;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim = Mesh->Dim;

  // find element geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  if (ierr != xf_OK) return ierr;

  // open file for writing
  ierr = xf_Error(xf_fopen(vollist, "w", &fid));
  if (ierr != xf_OK) return ierr;

  // loop over elements, write volumes
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    // basis and order of element
    QBasis = Mesh->ElemGroup[egrp].QBasis;
    // determine element Shape 
    ierr = xf_Error(xf_Basis2Shape(QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    // centroid in reference coords
    ierr = xf_Error(xf_ShapeCentroid(Shape, xref));
    if (ierr != xf_OK) return ierr;

    // perturb centroid so as to avoid trouble with uniformly-refined elems
    for (d=0; d<dim; d++) xref[d] += 0.000176234012*(d*d-2);  // pseudo-random

    PointsChanged = xfe_True;

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // element volume
      ElemVol = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
      
      // centroid in global coords
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, 
				      PointsChanged, 1, xref, xglob));
      if (ierr != xf_OK) return ierr;

      PointsChanged = xfe_False;
      
      // print out info
      for (d=0; d<dim; d++)
	fprintf(fid, "%.12E ", xglob[d]);
      fprintf(fid,"%.12E\n", ElemVol);

    } // elem
  } // egrp

  // close file
  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;

  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: VolListToRefIndicator
static int 
xf_VolListToRefIndicator(xf_All *All, const char *vollist,
			 xf_Vector *RefIndicator, int *pnMarked)
{
  int ierr, dim;
  int myRank, nProc, outerr;
  int egrp, elem;
  int count;
  enum xfe_Bool done, found, converged, inside;
  enum xfe_ShapeType Shape;
  char line[xf_MAXLONGLINELEN];
  real x[3], xref[3];
  real ElemVol, vol;
  xf_ElemSearchStruct ESS;
  FILE *fid;
  xf_Vector *EG, *ElemBoundBox = NULL;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  dim = Mesh->Dim;

  (*pnMarked) = 0;

  /* Determine myRank */
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // find element geometry vector
  ierr = xf_Error(xf_FindElemGeom(All, &EG));
  if (ierr != xf_OK) return ierr;

  // zero out indicator
  ierr = xf_Error(xf_SetZeroVector(RefIndicator));
  if (ierr != xf_OK) return ierr;

  // build element search structure
  ierr = xf_Error(xf_BuildElemSearchStructure(All, &ESS));
  if (ierr != xf_OK) return ierr;
    
  // open file for reading
  ierr = xf_Error(xf_fopen(vollist, "r", &fid));
  if (ierr != xf_OK) return ierr;

  // read in file line by line
  done = xfe_False;
  count = 0;
  while (!done){
    if (myRank == 0){
      if (fgets(line, xf_MAXLONGLINELEN, fid) == NULL) done = xfe_True;
      else{
	outerr = xf_OK;
	if (dim == 2)
	  ierr = sscanf(line, "%lf %lf %lf", x+0,x+1,&vol);
	else
	  ierr = sscanf(line, "%lf %lf %lf %lf", x+0,x+1,x+2,&vol);
	if (ierr != dim+1) outerr = xf_Error(xf_FILE_READ_ERROR);
      }
    }
    if (xf_PError(&outerr, 0) != xf_OK) return outerr;

    // are we done?
    ierr = xf_Error(xf_MPI_Bcast( (void *) &done, sizeof(enum xfe_Bool), 0));
    if (ierr != xf_OK) return ierr;
    if (done) break;
    
    count++;
    if ((count%100)==0) xf_printf("Read %7d lines \n", count);

    // not done yet
    ierr = xf_Error(xf_MPI_Bcast( (void *) x, 3*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_MPI_Bcast( (void *) &vol, sizeof(real), 0));
    if (ierr != xf_OK) return ierr;

    // find element that contains the point
    ierr = xf_FindElemUsingSearchStructure(All, x, &ESS, &egrp, &elem, xref);
    if (ierr == xf_NOT_FOUND){
      xf_printf("Warning, element containing point not found.  Continuing.\n");
      continue;
    }
    if (ierr != xf_OK) return ierr;
    found = ((egrp >= 0) && (elem >= 0));

/*     found = xfe_False; */
/*     for (egrp=0; egrp<Mesh->nElemGroup; egrp++){ */
/*       // determine element Shape  */
/*       ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape)); */
/*       if (ierr != xf_OK) return ierr; */
/*       for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){ */
/* 	// get reference coordinates */
/* 	ierr = xf_Error(xf_Glob2RefElem(Mesh, egrp, elem, x, 1e-4, xfe_False, */
/* 					xref, &converged)); */
/* 	if (ierr != xf_OK) return ierr; */
/* 	if (!converged) continue; // skip elements where glob2refelem fails */
/* 	// if inside, found=True, break */
/* 	ierr = xf_Error(xf_InsideShape(Shape, xref, MEPS, &inside)); */
/* 	if (ierr != xf_OK) return ierr; */
/* 	if (inside){ */
/* 	  found = xfe_True; */
/* 	  break; */
/* 	} */
/*       } // elem */
/*       if (found) break; */
/*     } // egrp */

    if (found){
      // mark element if requested volume is smaller
      ElemVol = EG->GenArray[egrp].rValue[elem][xfe_EGVolume];
      /*       xf_printf("-- egrp = %d, elem = %d, xglob = %.5E %.5E %.5E, vol=%.6E, ElemVol=%.6E\n",  */
      /* 		egrp, elem, x[0], x[1], x[2], vol, ElemVol); fflush(stdout); */
      /*       xf_printf("   xref  = %.5E %.5E %.5E\n", xref[0], xref[1], xref[2]); */
      
      /*       ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL,  */
      /* 				      xfe_True, 1, xref, x)); */
      /*       if (ierr != xf_OK) return ierr; */
      /*       xf_printf("   xflob = %.5E %.5E %.5E\n", x[0], x[1], x[2]); */
      
      if (ElemVol > 1.001*vol){
	if (RefIndicator->GenArray[egrp].iValue[elem][0] == 0){
	  (*pnMarked) = (*pnMarked)+1;
	  RefIndicator->GenArray[egrp].iValue[elem][0] = 1;
	}
      }
    }

    // reduce to max
    ierr = xf_Error(xf_MPI_Allreduce(pnMarked, 1, xfe_SizeInt, xfe_MPI_MAX));
    if (ierr != xf_OK) return ierr;
    
  } // until done

  // close file
  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;

  // destroy element search structure
  ierr = xf_Error(xf_DestroyElemSearchStructure(&ESS));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, len, nMarked, iter, maxIter;
  int myRank, nProc;
  char *ArgIn[] = {"in", "NULL", "input file (.xfa or .gri)",
		   "out", "NULL", "output file (.xfa or .gri)",
		   "vollist", "NULL", "text volume list",
		   "create", "False", "Create volume list",
		   "maxiter", "10", "Max number of iterations",
		   "\0"};
  char inFile[xf_MAXSTRLEN];  
  char outFile[xf_MAXSTRLEN];
  char vollist[xf_MAXSTRLEN];
  char *inExt, *outExt;
  enum xfe_Bool createvollist, done;
  xf_Vector *RefIndicator;
  xf_TimeHistData *TimeHistData = NULL;
  xf_KeyValue KeyValue;
  xf_All *All;


  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  xf_printf("\n");
  xf_printf("=== xf_AdaptToVol: Mesh adaptation based on volume ===\n");
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

  // Get vollist
  ierr = xf_GetKeyValue(KeyValue, "vollist", vollist);
  if (ierr != xf_OK) return ierr;

  // Get createvollist
  ierr = xf_GetKeyValueBool(KeyValue, "create", &createvollist);
  if (ierr != xf_OK) return ierr;

  // Get maxIter
  ierr = xf_GetKeyValueInt(KeyValue, "maxiter", &maxIter);
  if (ierr != xf_OK) return ierr;


  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  // extensions are required
  if ((len = strlen(inFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
  inExt = inFile + len - 4; // pointer to extension

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
  }
  else return xf_Error(xf_NOT_SUPPORTED);


  if (createvollist){
    /* Create volume list from existing All */
    ierr = xf_Error(xf_CreateVolList(All, vollist));
    if (ierr != xf_OK) return ierr;
  }
  else{
    /* Use volume list to adapt All, iteratively */
    done = xfe_False;
    for (iter=0; iter<maxIter; iter++){
      
      /* Delete non-essential data -- for memory cleanup */
      ierr = xf_Error(xf_DataSetDeleteNonEssential(All->DataSet));
      if (ierr != xf_OK) return ierr;

      // create/find a refinement indicator
      ierr = xf_Error(xf_FindVector(All, "RefIndicator", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				    NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_False,  xfe_True, NULL, 
				    &RefIndicator, NULL));
      if (ierr != xf_OK) return ierr;
      
      // fill the indicator using volume information
      ierr = xf_Error(xf_VolListToRefIndicator(All, vollist, RefIndicator, &nMarked));
      if (ierr != xf_OK) return ierr;
      
      xf_printf("Iteration %d, nMarked = %d\n", iter, nMarked);

      if (nMarked == 0){
	done = xfe_True;
	break; // nothing marked, so done!
      }
      
      /* Hanging-node adaptation (uniform) */
      ierr = xf_Error(xf_AdaptHang(All, RefIndicator));
      if (ierr != xf_OK) return ierr;
      
    } // iiter

    if (!done)
      xf_printf(" Volume adaptation not converged after %d iterations; nMarked = %d\n",
		maxIter, nMarked);
  }

  
  /* Write output file */
  if (xf_NotNull(outFile)){
    if ((len = strlen(outFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
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
  }

  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;
  

  xf_printf("xf_AdaptToVol finished.\n");  

  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
