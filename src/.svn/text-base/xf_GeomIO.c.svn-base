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
  FILE:  xf_GeomIO.c

  This file contains input/output functions pertinent to the Geometry
  data sructure.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_Mesh.h"
#include "xf_Basis.h"


/*---------------*/
/* Constructors  */
/*---------------*/  


/******************************************************************/
//   FUNCTION Definition: xf_CreateGeomCompAnalytical
static int 
xf_CreateGeomCompAnalytical( xf_GeomCompAnalytical **pGeomCompAnalytical){
  
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pGeomCompAnalytical, 1, sizeof(xf_GeomCompAnalytical)));
  if (ierr != xf_OK) return ierr;

  /* Initialize key-value structure */
  ierr = xf_Error(xf_InitKeyValue(&(*pGeomCompAnalytical)->KeyValue));
  if (ierr != xf_OK) return ierr;

  (*pGeomCompAnalytical)->Object = 0;
  (*pGeomCompAnalytical)->RParam = NULL;
  (*pGeomCompAnalytical)->IParam = NULL;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateGeomCompSpline
static int 
xf_CreateGeomCompSpline( xf_GeomCompSpline **pGeomCompSpline){
  
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pGeomCompSpline, 1, sizeof(xf_GeomCompSpline)));
  if (ierr != xf_OK) return ierr;

  (*pGeomCompSpline)->Order     = 0;
  (*pGeomCompSpline)->N         = 0;
  (*pGeomCompSpline)->X         = NULL;
  (*pGeomCompSpline)->Y         = NULL;
  (*pGeomCompSpline)->S         = NULL;
  (*pGeomCompSpline)->XS        = NULL;
  (*pGeomCompSpline)->YS        = NULL;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateGeomCompPanel
static int 
xf_CreateGeomCompPanel( xf_GeomCompPanel **pGeomCompPanel){
  
  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pGeomCompPanel, 1, sizeof(xf_GeomCompPanel)));
  if (ierr != xf_OK) return ierr;

  (*pGeomCompPanel)->nNode     = 0;
  (*pGeomCompPanel)->dim       = 0;
  (*pGeomCompPanel)->nPanel    = 0;
  (*pGeomCompPanel)->Basis     = 0;
  (*pGeomCompPanel)->Order     = 1;
  (*pGeomCompPanel)->Coord     = NULL;
  (*pGeomCompPanel)->Panels    = NULL;

  return xf_OK;
}



/*-----------------*/
/* Parallelization */
/*-----------------*/  


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeGeomCompAnalytical
static int 
xf_ParallelizeGeomCompAnalytical( xf_GeomCompAnalytical *GeomCompAnalytical){
  
  int ierr;
  int myRank, nProc;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // KeyValue
  ierr = xf_Error(xf_ParallelizeKeyValue(&GeomCompAnalytical->KeyValue));
  if (ierr != xf_OK) return ierr;

  // Object
  ierr = xf_Error(xf_MPI_Bcast((void *) &GeomCompAnalytical->Object, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeGeomCompSpline
static int 
xf_ParallelizeGeomCompSpline( xf_GeomCompSpline *GeomCompSpline){
  
  int ierr, N;
  int myRank, nProc;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // Order
  ierr = xf_Error(xf_MPI_Bcast((void *) &GeomCompSpline->Order, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  // N
  ierr = xf_Error(xf_MPI_Bcast((void *) &GeomCompSpline->N, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  if ((N = GeomCompSpline->N) <= 0) return xf_OK; // no points; nothing else to do

  // allocate X and Y on all procs other than root
  if (myRank > 0){
    ierr = xf_Error(xf_Alloc( (void **) &GeomCompSpline->X, GeomCompSpline->N, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &GeomCompSpline->Y, GeomCompSpline->N, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  // X
  ierr = xf_Error(xf_MPI_Bcast((void *) GeomCompSpline->X, N*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;

  // Y
  ierr = xf_Error(xf_MPI_Bcast((void *) GeomCompSpline->Y, N*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeGeom
int 
xf_ParallelizeGeom( xf_Geom *Geom)
{
  int ierr, i;
  int myRank, nProc;
  xf_GeomComp *pGeomComp;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if (nProc == 1) return xf_OK; // nothing to do

  // broadcast spatial dimension
  ierr = xf_Error(xf_MPI_Bcast((void *) &Geom->Dim, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  // broadcast number of components
  ierr = xf_Error(xf_MPI_Bcast((void *) &Geom->nComp, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  // all procs except for root allocate components
  if (myRank != 0){
    ierr = xf_Error(xf_Alloc( (void **) &(Geom->Comp), Geom->nComp,
			      sizeof(xf_GeomComp)));
    if (ierr != xf_OK) return ierr;
  }

  // Parallelize components
  for (i=0; i<Geom->nComp; i++){

    pGeomComp = Geom->Comp+i;

    // Name
    ierr = xf_Error(xf_ParallelizeString(&pGeomComp->Name));
    if (ierr != xf_OK) return ierr;

    // BFGTitle
    ierr = xf_Error(xf_ParallelizeString(&pGeomComp->BFGTitle));
    if (ierr != xf_OK) return ierr;

    // Type
    ierr = xf_Error(xf_MPI_Bcast((void *) &pGeomComp->Type, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;

    // Data -- depends what kind of component this is
    switch (pGeomComp->Type){
    case xfe_GeomCompNone:
      pGeomComp->Data = NULL;
      break;
    case xfe_GeomCompAnalytical:
      if (myRank > 0){ // create analytical geometry component on non-root procs
	ierr = xf_Error(xf_CreateGeomCompAnalytical((xf_GeomCompAnalytical **) &pGeomComp->Data));
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_ParallelizeGeomCompAnalytical((xf_GeomCompAnalytical *) pGeomComp->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompSpline:
      if (myRank > 0){ // create spline geometry component on non-root procs
	ierr = xf_Error(xf_CreateGeomCompSpline((xf_GeomCompSpline **) &pGeomComp->Data));
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_Error(xf_ParallelizeGeomCompSpline((xf_GeomCompSpline *) pGeomComp->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompPanel:
      if (myRank > 0){ // for now, paneling exists only on root, as it is only used by root
	pGeomComp->Data = NULL;
      }
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    } 
  }

  return xf_OK;
}


/*-----------------*/
/* Text File Input */
/*-----------------*/ 

/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomCompAnalytical
static int 
xf_ReadGeomCompAnalytical(FILE *fgeom, char **InputStrings, int *piString,
			  xf_GeomCompAnalytical **pGeomCompAnalytical)
{
  int ierr;
  enum xfe_Bool done;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line0[200], *line;
  xf_GeomCompAnalytical *GCA;


  ierr = xf_Error(xf_CreateGeomCompAnalytical(pGeomCompAnalytical));
  if (ierr != xf_OK) return ierr;

  GCA = (*pGeomCompAnalytical);

  /* Description of analytical object consists of just key-value pairs */
  done = xfe_False;
  do{
    ierr = xf_LineFromFileOrStrings(fgeom, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line

    if (strncmp(line, "ENDCOMPONENT", 12) == 0){
      done = xfe_True;
      break;
    }
    
    /* reading a key=value pair */
    ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
    if (ierr != xf_OK) return ierr;
    
    /* Store key and value in local list */
    ierr = xf_Error(xf_AddKeyValue(&GCA->KeyValue, key, value, xfe_True));
    if (ierr == xf_OVERWROTE){
      xf_printf("Error. Key %s assigned more than once in geom analytical term.\n", key);
      return xf_FILE_READ_ERROR;
    }
    if (ierr != xf_OK) return ierr;
    
  } while (!done);
    
  // ENDCOMPONENT line not found
  if (!done) return xf_Error(xf_FILE_READ_ERROR);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomCompSpline
static int 
xf_ReadGeomCompSpline(FILE *fgeom, char **InputStrings, int *piString,
		      xf_GeomCompSpline **pGeomCompSpline)
{
  int ierr;
  int k;
  enum xfe_Bool done;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line0[200], *line;
  xf_GeomCompSpline *GCS;

  ierr = xf_Error(xf_CreateGeomCompSpline(pGeomCompSpline));
  if (ierr != xf_OK) return ierr;

  GCS = (*pGeomCompSpline);

  // default values
  GCS->Order = 3;

  /* Read description of spline object */
  done = xfe_False;
  do{
    ierr = xf_LineFromFileOrStrings(fgeom, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line

    if (strncmp(line, "ENDCOMPONENT", 12) == 0){
      done = xfe_True;
      break;
    }

    // Order
    if (strncmp(line, "Order", 5) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      if (sscanf(value, "%d", &GCS->Order) != 1) return xf_Error(xf_FILE_READ_ERROR);
      continue;	      
    }
    
    // Points
    if (strncmp(line, "Points", 6) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      if (strncmp(value, "Inline", 6) != 0) return xf_Error(xf_NOT_SUPPORTED);

      // first line is number of points
      ierr = xf_Error(xf_LineFromFileOrStrings(fgeom, InputStrings, piString, line0, &line));
      if (ierr != xf_OK) return ierr;
      ierr = sscanf(line, "%d", &GCS->N);
      if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);

      // next are all x,y coordinates, one point per line
      ierr = xf_Error(xf_Alloc( (void **) &GCS->X, GCS->N, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_Alloc( (void **) &GCS->Y, GCS->N, sizeof(real)));
      if (ierr != xf_OK) return ierr;
      for (k=0; k<GCS->N; k++){
	ierr = xf_Error(xf_LineFromFileOrStrings(fgeom, InputStrings, piString, line0, &line));
	if (ierr != xf_OK) return ierr;

	ierr = sscanf(line, "%lf %lf", GCS->X+k, GCS->Y+k);
	if (ierr != 2) return xf_Error(xf_FILE_READ_ERROR);
      } // k

      continue;
    }

    xf_printf("Unrecognized line in spline geometry description.\n");
    return xf_Error(xf_FILE_READ_ERROR);
    
  } while (!done);
    
  // ENDCOMPONENT line not found
  if (!done) return xf_Error(xf_FILE_READ_ERROR);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadPanelNodeFile
static int 
xf_ReadPanelNodeFile(char *fname, xf_GeomCompPanel *GCP)
{
  int ierr, i;
  char line0[200], *line;
  FILE *fid = NULL;

  if ((fid = fopen(fname, "r")) == NULL){
    xf_printf("Could not find panel node file: %s\n", fname);
    return xf_Error(xf_FILE_READ_ERROR);
  }

  // first line is: nNode dim
  ierr = xf_Error(xf_LineFromFileOrStrings(fid, NULL, NULL, line0, &line));
  if (ierr != xf_OK) return ierr;
  ierr = sscanf(line, "%d %d", &GCP->nNode, &GCP->dim);
  if (ierr != 2) return xf_Error(xf_FILE_READ_ERROR);
  
  // allocate memory
  ierr = xf_Error(xf_Alloc2((void ***) &GCP->Coord, GCP->nNode, GCP->dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // next is all the nodes
  for (i=0; i<GCP->nNode; i++){
    ierr = xf_Error(xf_LineFromFileOrStrings(fid, NULL, NULL, line0, &line));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ScanReal(line, GCP->dim, GCP->Coord[i]));
    if (ierr != xf_OK) return ierr;
  } // i

  fclose(fid);
 
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadPanelPanelsFile
static int 
xf_ReadPanelPanelsFile(char *fname, xf_GeomCompPanel *GCP)
{
  int ierr, i, j, nn;
  char line0[200], *line;
  char BasisString[xf_MAXSTRLEN];
  FILE *fid = NULL;

  if ((fid = fopen(fname, "r")) == NULL){
    xf_printf("Could not find panel panels file: %s\n", fname);
    return xf_Error(xf_FILE_READ_ERROR);
  }

  // first line is: nPanel BasisString Order
  ierr = xf_Error(xf_LineFromFileOrStrings(fid, NULL, NULL, line0, &line));
  if (ierr != xf_OK) return ierr;
  ierr = sscanf(line, "%d %s %d", &GCP->nPanel, BasisString, &GCP->Order);
  if (ierr != 3) return xf_Error(xf_FILE_READ_ERROR);
  ierr = xf_Error(xf_Value2Enum(BasisString, xfe_BasisName, xfe_BasisLast, 
				(int *) &GCP->Basis));
  
  // how many nodes per panel do we expect?
  ierr = xf_Error(xf_Order2nNode(GCP->Basis, GCP->Order, &nn));
  if (ierr != xf_OK) return ierr;
  
  // allocate memory
  ierr = xf_Error(xf_Alloc2((void ***) &GCP->Panels, GCP->nPanel, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // next is all the panels
  for (i=0; i<GCP->nPanel; i++){
    ierr = xf_Error(xf_LineFromFileOrStrings(fid, NULL, NULL, line0, &line));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ScanInt(line, nn, GCP->Panels[i]));
    if (ierr != xf_OK) return ierr;
    for (j=0; j<nn; j++) GCP->Panels[i][j]--; // file is assumed 1-based; we want 0-based
  } // i

  fclose(fid);
 
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomCompPanel
static int 
xf_ReadGeomCompPanel(FILE *fgeom, char **InputStrings, int *piString,
		      xf_GeomCompPanel **pGeomCompPanel)
{
  int ierr;
  int k;
  enum xfe_Bool done;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line0[200], *line;
  xf_GeomCompPanel *GCP;

  ierr = xf_Error(xf_CreateGeomCompPanel(pGeomCompPanel));
  if (ierr != xf_OK) return ierr;

  GCP = (*pGeomCompPanel);

  // default values
  GCP->Order = 1;
  GCP->Basis = xfe_TriLagrange;

  /* Read description of panel object */
  done = xfe_False;
  do{
    ierr = xf_LineFromFileOrStrings(fgeom, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line

    if (strncmp(line, "ENDCOMPONENT", 12) == 0){
      done = xfe_True;
      break;
    }
    
    // Nodes
    if (strncmp(line, "Nodes", 5) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      // for now only support file names
      if (strncmp(value, "Inline", 6) == 0) return xf_Error(xf_NOT_SUPPORTED);
      ierr = xf_Error(xf_ReadPanelNodeFile(value, GCP));
      if (ierr != xf_OK) return ierr;
      continue;
    }

    // Panels
    if (strncmp(line, "Panels", 5) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      // for now only support file names
      if (strncmp(value, "Inline", 6) == 0) return xf_Error(xf_NOT_SUPPORTED);
      ierr = xf_Error(xf_ReadPanelPanelsFile(value, GCP));
      if (ierr != xf_OK) return ierr;
      continue;
    }

    xf_printf("Unrecognized line in panel geometry description.\n");
    return xf_Error(xf_FILE_READ_ERROR);
    
  } while (!done);
    
  // ENDCOMPONENT line not found
  if (!done) return xf_Error(xf_FILE_READ_ERROR);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomFileSerial
static int 
xf_ReadGeomFileSerial(char *GeomFile, char **InputStrings, xf_Geom *Geom)
{
  int ierr, iString = 0;
  int nComp, iComp;
  enum xfe_Bool done;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line0[200], *line;
  FILE *fgeom = NULL;

  /* Check input */
  if (((GeomFile == NULL) && (InputStrings == NULL)) || 
      ((GeomFile != NULL) && (InputStrings != NULL))) 
    return xf_Error(xf_INPUT_ERROR);

  if (GeomFile != NULL){
    /* Open file */
    if ((fgeom = fopen(GeomFile, "r")) == NULL){
      xf_printf("Could not find Geometry file: %s\n", GeomFile);
      xf_printf("Make sure the appropriate extension is included.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
  }
  else{
    // Start at the first string
    iString = 0;
  }

  /* Read dimension and number of components */
  nComp = -1;
  do{
    ierr = xf_LineFromFileOrStrings(fgeom, InputStrings, &iString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;
    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line

    // dimension must precede nComponent
    if (strncmp(line, "Dim", 3) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      if (sscanf(value, "%d", &Geom->Dim) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (Geom->Dim < 1) return xf_Error(xf_FILE_READ_ERROR);
    }

    if (strncmp(line, "nComponent", 10) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      if (sscanf(value, "%d", &nComp) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (nComp <= 0) return xf_Error(xf_FILE_READ_ERROR);
    }

  } while (nComp < 0);

  if (nComp < 0){
    xf_printf("Error: no geometry components specified in geometry file.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }

  // allocate components
  Geom->nComp = nComp;
  ierr = xf_Error(xf_Alloc( (void **) &(Geom->Comp), Geom->nComp,
			    sizeof(xf_GeomComp)));
  if (ierr != xf_OK) return ierr;


  /* Loop over components, read each one in turn */
  for (iComp=0; iComp<nComp; iComp++){
    
    // read component name and type
    done = xfe_False;
    do{
      ierr = xf_LineFromFileOrStrings(fgeom, InputStrings, &iString, line0, &line);
      if (ierr == xf_END_OF_FILE) break;
      if (ierr != xf_OK) return ierr;
      if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line

      // geometry component name is optional, must precede type
      if (strncmp(line, "ComponentName", 13) == 0){
	ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_AllocString(&Geom->Comp[iComp].Name, xf_MAXSTRLEN, value));
	if (ierr != xf_OK) return ierr;
      }

      // boundary face group title is optional, must precede type
      if (strncmp(line, "BFGTitle", 8) == 0){
	ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_AllocString(&Geom->Comp[iComp].BFGTitle, xf_MAXSTRLEN, value));
	if (ierr != xf_OK) return ierr;
      }

      // geometry component type is required
      if (strncmp(line, "Type", 4) == 0){
	ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
	if (ierr != xf_OK) return ierr;
	
	ierr = xf_Error(xf_Value2Enum(value, xfe_GeomCompName, xfe_GeomCompLast, 
				      (int *) &Geom->Comp[iComp].Type));
	if( ierr != xf_OK ) return ierr;
	done = xfe_True;
      }

    } while (!done);
    
    // inconsistent number of components
    if (!done) return xf_Error(xf_FILE_READ_ERROR);

    // Call appropriate reader
    switch(Geom->Comp[iComp].Type){
    case xfe_GeomCompAnalytical:
      ierr = xf_Error(xf_ReadGeomCompAnalytical(fgeom, InputStrings, &iString, 
						(xf_GeomCompAnalytical **) 
						&Geom->Comp[iComp].Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompSpline:
      ierr = xf_Error(xf_ReadGeomCompSpline(fgeom, InputStrings, &iString, 
					    (xf_GeomCompSpline **) 
					    &Geom->Comp[iComp].Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompPanel:
      ierr = xf_Error(xf_ReadGeomCompPanel(fgeom, InputStrings, &iString,
					   (xf_GeomCompPanel **) 
					   &Geom->Comp[iComp].Data));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    }
    
  }  // iComp


  if (fgeom != NULL) fclose(fgeom);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomFile
int 
xf_ReadGeomFile(char *GeomFile, char **InputStrings, xf_Geom *Geom)
{
  int ierr, terr, myRank;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root reads in geometry file
  if (myRank == 0)
    terr = xf_Error(xf_ReadGeomFileSerial(GeomFile, InputStrings, Geom));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  // parallelize Geom structure
  ierr = xf_Error(xf_ParallelizeGeom(Geom));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}




/*---------------------*/
/* Binary Input/Output */
/*---------------------*/ 


/******************************************************************/
//   FUNCTION Definition: xf_WriteGeomCompAnalyticalBinary
static int 
xf_WriteGeomCompAnalyticalBinary( xf_GeomCompAnalytical *GCA, FILE *fid)
{
  int ierr, rev;
    
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // KeyValue
  ierr = xf_Error(xf_WriteKeyValueBinary(GCA->KeyValue, fid));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomCompAnalyticalBinary
static int 
xf_ReadGeomCompAnalyticalBinary( FILE *fid, xf_GeomCompAnalytical *GCA )
{
  int ierr, rev, si, i;
  enum xfe_Bool flag;

  si = sizeof(int);

  rev = 0;  // writer revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // KeyValue
  ierr = xf_Error(xf_ReadKeyValueBinary(fid, &GCA->KeyValue));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteGeomCompSplineBinary
static int 
xf_WriteGeomCompSplineBinary( xf_GeomCompSpline *GCS, FILE *fid)
{
  int ierr, rev;
    
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // Order
  if (fwrite(&GCS->Order, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // N
  if (fwrite(&GCS->N, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // X and Y
  if (GCS->N > 0){
    if (fwrite(GCS->X, sizeof(real), GCS->N, fid) != GCS->N) 
      return xf_Error(xf_FILE_WRITE_ERROR);
    if (fwrite(GCS->Y, sizeof(real), GCS->N, fid) != GCS->N) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomCompSplineBinary
static int 
xf_ReadGeomCompSplineBinary( FILE *fid, xf_GeomCompSpline *GCS )
{
  int ierr, rev, si, i;
  enum xfe_Bool flag;

  si = sizeof(int);

  rev = 0;  // writer revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);

  // Order
  if (fread(&GCS->Order, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  // N  
  if (fread(&GCS->N, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  // X 
  if (GCS->N > 0){

    ierr = xf_Error(xf_Alloc( (void **) &GCS->X, GCS->N, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &GCS->Y, GCS->N, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    if (fread(GCS->X, sizeof(real), GCS->N, fid) != GCS->N) 
      return xf_Error(xf_FILE_READ_ERROR);
    if (fread(GCS->Y, sizeof(real), GCS->N, fid) != GCS->N) 
      return xf_Error(xf_FILE_READ_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteGeomCompPanelBinary
static int 
xf_WriteGeomCompPanelBinary( xf_GeomCompPanel *GCP, FILE *fid)
{
  int ierr, rev, nn;
    
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // nNode
  if (fwrite(&GCP->nNode, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // dim
  if (fwrite(&GCP->dim, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // Coord
  if (GCP->nNode*GCP->dim > 0){
    if (fwrite(GCP->Coord[0], sizeof(real), GCP->nNode*GCP->dim, fid) != GCP->nNode*GCP->dim) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }

  // nPanel
  if (fwrite(&GCP->nPanel, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Basis
  ierr = xf_Error(xf_WriteStringBinary(xfe_BasisName[GCP->Basis], fid));
  if (ierr != xf_OK) return ierr;

  // Order
  if (fwrite(&GCP->Order, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  ierr = xf_Error(xf_Order2nNode(GCP->Basis, GCP->Order, &nn));
  if (ierr != xf_OK) return ierr;

  // Panels
  if (GCP->nPanel*nn > 0){
    if (fwrite(GCP->Panels[0], sizeof(int), GCP->nPanel*nn, fid) != GCP->nPanel*nn) 
      return xf_Error(xf_FILE_WRITE_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomCompPanelBinary
static int 
xf_ReadGeomCompPanelBinary( FILE *fid, xf_GeomCompPanel *GCP )
{
  int ierr, rev, si, i, nn;
  enum xfe_Bool flag;

  si = sizeof(int);

  rev = 0;  // writer revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);

  // nNode
  if (fread(&GCP->nNode, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  // dim
  if (fread(&GCP->dim, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  // Coord
  if (GCP->nNode*GCP->dim > 0){
    ierr = xf_Error(xf_Alloc2( (void ***) &GCP->Coord, GCP->nNode, GCP->dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    if (fread(GCP->Coord[0], sizeof(real), GCP->nNode*GCP->dim, fid) != GCP->nNode*GCP->dim) 
      return xf_Error(xf_FILE_READ_ERROR);
  }

  // nPanel
  if (fread(&GCP->nPanel, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Basis
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BasisName, xfe_BasisLast, 
				    (int *) &GCP->Basis));
  if (ierr != xf_OK) return ierr;

  // Order
  if (fread(&GCP->Order, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);

  ierr = xf_Error(xf_Order2nNode(GCP->Basis, GCP->Order, &nn));
  if (ierr != xf_OK) return ierr;


  // Panels
  if (GCP->nPanel*nn > 0){
    ierr = xf_Error(xf_Alloc2( (void ***) &GCP->Panels, GCP->nPanel, nn, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    if (fread(GCP->Panels[0], sizeof(int), GCP->nPanel*nn, fid) != GCP->nPanel*nn) 
      return xf_Error(xf_FILE_READ_ERROR);
  }


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteGeomBinarySerial
static int 
xf_WriteGeomBinarySerial( xf_Geom *Geom, FILE *fid )
{
  int ierr, rev, si, i;
  int iComp;
  enum xfe_Bool flag;
  xf_GeomComp *GC;

  si = sizeof(int);

  rev = 1;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // dimension
  if (fwrite(&Geom->Dim, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // number of components
  if (fwrite(&Geom->nComp, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  /* Loop over components, write each one in turn */
  for (iComp=0; iComp<Geom->nComp; iComp++){
    
    GC = Geom->Comp+iComp;

    // Name
    ierr = xf_Error(xf_WriteStringBinary(GC->Name, fid));
    if (ierr != xf_OK) return ierr;

    // BFGTitle
    ierr = xf_Error(xf_WriteStringBinary(GC->BFGTitle, fid));
    if (ierr != xf_OK) return ierr;

    // Type
    ierr = xf_Error(xf_WriteStringBinary(xfe_GeomCompName[GC->Type], fid));
    if (ierr != xf_OK) return ierr;

    // Call appropriate writer for data
    switch(Geom->Comp[iComp].Type){
    case xfe_GeomCompAnalytical:
      ierr = xf_Error(xf_WriteGeomCompAnalyticalBinary((xf_GeomCompAnalytical *) GC->Data, fid));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompSpline:
      ierr = xf_Error(xf_WriteGeomCompSplineBinary((xf_GeomCompSpline *) GC->Data, fid));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompPanel:
      ierr = xf_Error(xf_WriteGeomCompPanelBinary((xf_GeomCompPanel *) GC->Data, fid));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    }

  }  // iComp
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomBinarySerial
static int 
xf_ReadGeomBinarySerial(  FILE *fid, xf_Geom *Geom)
{
  int ierr, rev, si, i;
  int iComp;
  char title[xf_MAXSTRLEN];
  enum xfe_Bool flag;
  xf_GeomComp *GC;

  si = sizeof(int);

  rev = 1;  // writer revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);

  // backwards-compatibility support for version 0 geometry
  if (rev == 0){
    // Type
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, title, NULL));
    if (ierr != xf_OK) return ierr;
    // SurfMesh
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
    if (ierr != xf_OK) return ierr;
    if (flag) return xf_Error(xf_NOT_SUPPORTED);
    // AuxData
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_BoolName, xfe_BoolLast, (int *) &flag));
    if (ierr != xf_OK) return ierr;
    if (flag) return xf_Error(xf_NOT_SUPPORTED);

    return xf_OK;
  }

  if (rev != 1) return xf_Error(xf_FILE_READ_ERROR);

  // dimension
  if (fread(&Geom->Dim, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // number of components
  if (fread(&Geom->nComp, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (Geom->nComp < 0) return xf_Error(xf_FILE_READ_ERROR);

  // allocate component vector
  ierr = xf_Error(xf_Alloc( (void **) &(Geom->Comp), Geom->nComp,
			    sizeof(xf_GeomComp)));
  if (ierr != xf_OK) return ierr;

  /* Loop over components, read each one in turn */
  for (iComp=0; iComp<Geom->nComp; iComp++){
    
    GC = Geom->Comp+iComp;

    // Name
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &GC->Name));
    if (ierr != xf_OK) return ierr;

    // BFGTitle
    ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &GC->BFGTitle));
    if (ierr != xf_OK) return ierr;

    // Type
    ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_GeomCompName, xfe_GeomCompLast, 
				      (int *) &GC->Type));
    if (ierr != xf_OK) return ierr;

    // Call appropriate reader for data
    switch(Geom->Comp[iComp].Type){
    case xfe_GeomCompAnalytical:
      ierr = xf_Error(xf_CreateGeomCompAnalytical((xf_GeomCompAnalytical **) &GC->Data));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReadGeomCompAnalyticalBinary(fid, (xf_GeomCompAnalytical *) GC->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompSpline:
      ierr = xf_Error(xf_CreateGeomCompSpline((xf_GeomCompSpline **) &GC->Data));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReadGeomCompSplineBinary(fid, (xf_GeomCompSpline *) GC->Data));
      if (ierr != xf_OK) return ierr;
      break;
    case xfe_GeomCompPanel:
      ierr = xf_Error(xf_CreateGeomCompPanel((xf_GeomCompPanel **) &GC->Data));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ReadGeomCompPanelBinary(fid, (xf_GeomCompPanel *) GC->Data));
      if (ierr != xf_OK) return ierr;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    }

  }  // iComp
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteGeomBinary
int 
xf_WriteGeomBinary(xf_Geom *Geom, FILE *fid)
{
  int ierr, terr, myRank;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root writes
  if (myRank == 0)
    terr = xf_Error(xf_WriteGeomBinarySerial(Geom, fid));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadGeomBinary
int 
xf_ReadGeomBinary(  FILE *fid, xf_Geom *Geom)
{
  int ierr, terr, myRank;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root reads 
  if (myRank == 0)
    terr = xf_Error(xf_ReadGeomBinarySerial(fid, Geom));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  // parallelize Geom
  ierr = xf_Error(xf_ParallelizeGeom(Geom));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

