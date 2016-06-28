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
 FILE:  xf_All.c
 
 This file contains functions for working with the All data structure.
 
 */

#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_Mesh.h"
#include "xf_MeshMotion.h"
#include "xf_MeshMotionIO.h"
#include "xf_Geom.h"
#include "xf_GeomIO.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_EqnSet.h"
#include "xf_String.h"
#include "xf_ParamDefault.h"
#include "xf_EqnSetHook.h"

/******************************************************************/
//   FUNCTION Definition: xf_CreateAll
int 
xf_CreateAll( xf_All **pAll, enum xfe_Bool DefaultFlag){
  
  int ierr;
  xf_All *All;
  
  ierr = xf_Error(xf_Alloc((void **) pAll, 1, sizeof(xf_All)));
  if (ierr != xf_OK) return ierr;
  
  All = (*pAll);
  
  ierr = xf_Error(xf_CreateMesh(&(All->Mesh)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_CreateGeom(&(All->Geom)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_CreateDataSet(&(All->DataSet)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_CreateParam(&(All->Param)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_CreateEqnSet(&(All->EqnSet)));
  if (ierr != xf_OK) return ierr;
  
  if (DefaultFlag){
    // Params default values
    ierr = xf_Error(xf_AddKeyValueList(&(*pAll)->Param->KeyValue, xf_DefaultParamList,
                                       xfe_True, xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyAll
int 
xf_DestroyAll( xf_All *All){
  
  int ierr;
  
  ierr = xf_Error(xf_DestroyMesh(All->Mesh));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyGeom(All->Geom));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyDataSet(All->DataSet));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyParam(All->Param));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_DestroyEqnSet(All->EqnSet, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  xf_Release(All);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CheckAll
int 
xf_CheckAll( xf_All *All){
  
  int ierr;
  
  // check Mesh for consistency
  
  // check agreement between EqnSet and Mesh
  
  // check data linkage
  
  xf_printf("\n** NOT CHECKING ALL **\n\n");
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteAllBinary
int 
xf_WriteAllBinary( xf_All *All, const char *fname){
  
  int ierr, rev, len;
  int si, sc, sl;
  xf_long pos, pos2;
  FILE *fid;
  
  si = sizeof(int);
  sc = sizeof(char);
  sl = sizeof(xf_long);
  
  
  // open file for writing
  ierr = xf_Error(xf_fopen(fname, "wb", &fid));
  if (ierr != xf_OK) return ierr;
  
  // write header
  ierr = xf_Error(xf_WriteStringBinaryParallel("XFlow Binary Restart File", fid));
  if (ierr != xf_OK) return ierr;
  
  rev = 0;  // writer revision number
  ierr = xf_Error(xf_fwrite(&rev, si, 1, fid));
  if (ierr != xf_OK) return ierr;
  
  // write mesh
  ierr = xf_Error(xf_WriteStringBinaryParallel("Mesh", fid));
  if (ierr != xf_OK) return ierr;
  
  // placeholder for end of chunk position
  ierr = xf_Error(xf_ftell(fid, &pos)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_WriteMeshBinary(All->Mesh, fid));
  if (ierr != xf_OK) return ierr;
  
  // get/write end-of-chunk position
  ierr = xf_Error(xf_ftell(fid, &pos2));  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos, SEEK_SET)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos2, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos2, SEEK_SET)); if (ierr != xf_OK) return ierr;
  
  // write geometry
  ierr = xf_Error(xf_WriteStringBinaryParallel("Geom", fid));
  if (ierr != xf_OK) return ierr;
  
  // placeholder for end of chunk position
  ierr = xf_Error(xf_ftell(fid, &pos)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_WriteGeomBinary(All->Geom, fid));
  if (ierr != xf_OK) return ierr;
  
  // get/write end-of-chunk position
  ierr = xf_Error(xf_ftell(fid, &pos2));  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos, SEEK_SET)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos2, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos2, SEEK_SET)); if (ierr != xf_OK) return ierr;
  
  // write DataSet
  ierr = xf_Error(xf_WriteStringBinaryParallel("DataSet", fid));
  if (ierr != xf_OK) return ierr;
  
  // placeholder for end of chunk position
  ierr = xf_Error(xf_ftell(fid, &pos)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, All->DataSet, fid, NULL));
  if (ierr != xf_OK) return ierr;
  
  // get/write end-of-chunk position
  ierr = xf_Error(xf_ftell(fid, &pos2));  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos, SEEK_SET)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos2, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos2, SEEK_SET)); if (ierr != xf_OK) return ierr;
  
  // write Param
  ierr = xf_Error(xf_WriteStringBinaryParallel("Param", fid));
  if (ierr != xf_OK) return ierr;
  
  // placeholder for end of chunk position
  ierr = xf_Error(xf_ftell(fid, &pos)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_WriteParamBinary(All->Param, fid));
  if (ierr != xf_OK) return ierr;
  
  // get/write end-of-chunk position
  ierr = xf_Error(xf_ftell(fid, &pos2));  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos, SEEK_SET)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos2, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos2, SEEK_SET)); if (ierr != xf_OK) return ierr;
  
  // write EqnSet  
  ierr = xf_Error(xf_WriteStringBinaryParallel("EqnSet", fid));
  if (ierr != xf_OK) return ierr;
  
  // placeholder for end of chunk position
  ierr = xf_Error(xf_ftell(fid, &pos)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_WriteEqnSetBinary(All->EqnSet, fid));
  if (ierr != xf_OK) return ierr;
  
  // get/write end-of-chunk position
  ierr = xf_Error(xf_ftell(fid, &pos2));  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos, SEEK_SET)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fwrite(&pos2, sl, 1, fid)); if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_fseek(fid, pos2, SEEK_SET)); if (ierr != xf_OK) return ierr;
  
  
  // write end
  ierr = xf_Error(xf_WriteStringBinaryParallel("END ALL", fid));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadAllBinary
int 
xf_ReadAllBinary( const char *fname, xf_All *All){
  
  int ierr, rev, si, sl;
  enum xfe_Bool done;
  xf_long pos;
  char title[xf_MAXSTRLEN];
  FILE *fid;
  
  si = sizeof(int);
  sl = sizeof(xf_long);
  
  // open file for reading
  ierr = xf_Error(xf_fopen(fname, "rb", &fid));
  if (ierr != xf_OK) return ierr;
  
  // read header (discard)
  ierr = xf_Error(xf_ReadStringBinaryParallel(fid, xf_MAXSTRLEN, title, NULL));
  if (ierr != xf_OK) return ierr;
  
  // read revision number
  ierr = xf_Error(xf_fread(fid, si, 1, &rev));
  if (ierr != xf_OK) return ierr;
  
  // revision # check
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  
  // read chunks until hit END
  done = xfe_False;
  while (!done){
    if (xf_feof(fid)) return xf_Error(xf_FILE_READ_ERROR);
    
    // read chunk title
    ierr = xf_Error(xf_ReadStringBinaryParallel(fid, xf_MAXSTRLEN, title, NULL));
    if (ierr != xf_OK) return ierr;
    
    if (strncmp(title, "END ALL", 7) == 0){
      done = xfe_True;
      break;
    }
    
    // read pos of end of chunk
    ierr = xf_Error(xf_fread(fid, sl, 1, &pos));
    if (ierr != xf_OK) return ierr;
    
    // read chunk
    if (strncmp(title, "Mesh", 4) == 0){
      ierr = xf_Error(xf_ReadMeshBinary(fid, All->Mesh));
      if (ierr != xf_OK) return ierr;
    }
    else if (strncmp(title, "Geom", 4) == 0){
      ierr = xf_Error(xf_ReadGeomBinary(fid, All->Geom));
      if (ierr != xf_OK) return ierr;
    }
    else if (strncmp(title, "DataSet", 7) == 0){
      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, fid, NULL, All->DataSet));
      if (ierr != xf_OK) return ierr;
    }
    else if (strncmp(title, "Param", 5) == 0){
      ierr = xf_Error(xf_ReadParamBinary(fid, All->Param));
      if (ierr != xf_OK) return ierr;
    }
    else if (strncmp(title, "EqnSet", 6) == 0){
      ierr = xf_Error(xf_ReadEqnSetBinary(fid, All->EqnSet));
      if (ierr != xf_OK) return ierr;
    }
    else{
      xf_printf("Skipping chunk entitled %s (not recognized) in ReadAll.\n", title);
      ierr = xf_Error(xf_fseek(fid, pos, SEEK_SET));
      if (ierr != xf_OK) return ierr;
    }
    
  }
  
  // close file
  ierr = xf_Error(xf_fclose(fid));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadAsciiFile
static int
xf_ReadAsciiFile(const char *InputFile, xf_KeyValue *KeyValue, xf_Mesh *Mesh)
{
  /*
   
   PURPOSE: 
   
   Reads an ASCII InputFile into a Mesh structure.  This function is
   called by all processors when in parallel.
   
   INPUTS:
   
   InputFile : name of input file
   KeyValue  : KeyValue list containing any useful parameters (from the .job file)
   
   OUTPUTS: 
   
   Mesh : Mesh structure (preallocated in that *Mesh exists)
   
   RETURNS: Error Code
   
   */
  int len, ierr, terr;
  int myRank, nProc;
  const char *pext;
  xf_Mesh *Mesh_Glob = NULL;
  
  if ((len = strlen(InputFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;  
  
  pext = InputFile + len - 4; // pointer to extension
  
  if (strncmp(pext, ".gri", 4) == 0){
    
    // root creates Mesh_Glob and reads the ascii file
    if (myRank == 0){
      terr = xf_Error(xf_CreateMesh(&Mesh_Glob));
      if (terr == xf_OK)
        terr = xf_Error( xf_ReadGriFile(InputFile, NULL, Mesh_Glob) );
    }
    ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    if (terr != xf_OK) return xf_Error(terr);
    
    // Mesh_Glob is parallelized
    ierr = xf_Error(xf_ParallelizeMesh(Mesh_Glob, Mesh, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    
    // no longer need Mesh_Glob
    ierr = xf_Error(xf_DestroyMesh(Mesh_Glob));
    if (ierr != xf_OK) return ierr;
  }
  else if (strncmp(pext, "gmsh", 4) == 0){
    
    // root creates Mesh_Glob and reads the ascii file
    if (myRank == 0){
      terr = xf_Error(xf_CreateMesh(&Mesh_Glob));
      if (terr == xf_OK){
        ierr = xf_Error( xf_ReadGmshFile(InputFile, Mesh_Glob) );
        if (ierr != xf_OK) return ierr;
      }
    }
    ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    if (terr != xf_OK) return xf_Error(terr);
    
    // Mesh_Glob is parallelized
    ierr = xf_Error(xf_ParallelizeMesh(Mesh_Glob, Mesh, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    
    // no longer need Mesh_Glob
    ierr = xf_Error(xf_DestroyMesh(Mesh_Glob));
    if (ierr != xf_OK) return ierr;
  }
  else{
    xf_printf("InputFile = %s not supported in ReadAsciiFile.\n", InputFile);
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadAllInputFile
int 
xf_ReadAllInputFile( const char *InputFile, xf_KeyValue *KeyValue, 
                    enum xfe_Bool DefaultFlag, xf_All **pAll)
{
  int ierr;
  int len;
  const char *pext;
  char BackgroundMesh[xf_MAXSTRLEN];
  xf_All *All;
  
  /*----------------------------*/
  /* Create all structure, All, */
  /* and children structures    */
  /*----------------------------*/ 
  
  ierr = xf_Error(xf_CreateAll(pAll, xfe_False));
  if (ierr != xf_OK) return ierr;
  All = (*pAll);
  
  
  /*----------------*/
  /* Read InputFile */
  /*----------------*/
  
  len = strlen(InputFile);
  if (len < 4){
    xf_printf("Error, InputFile requires an extension.\n");
    xf_printf(" strlen(InputFile) = %d < 4.\n", len);
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  pext = InputFile + len - 4; // pointer to extension
  
  if (strncmp(pext, ".xfa", 4) == 0 ){
    // Binary .xfa file
    ierr = xf_Error(xf_ReadAllBinary(InputFile, All));
    if (ierr!=xf_OK) return ierr;
  }
  else if( (strncmp(pext, ".gri", 4) == 0) ||
	   (strncmp(pext, ".fro", 4) == 0) ||
	   (strncmp(pext, "gmsh", 4) == 0) ){
    // Ascii (text) file
    ierr = xf_Error(xf_ReadAsciiFile(InputFile, KeyValue, All->Mesh));
    if (ierr != xf_OK) return ierr;
  }
  else if (strncmp(pext, ".ebg", 4) == 0){
    // Ascii embedded boundary group file
    if (KeyValue == NULL) return xf_Error(xf_INPUT_ERROR); // need KeyValue here
    
    ierr = xf_Error(xf_GetKeyValue(*KeyValue, "BackgroundMesh", BackgroundMesh));
    if (ierr != xf_OK){
      xf_printf("Error reading job file.  Please specify BackgroundMesh.\n");
      return ierr;
    }
    
    // Read background mesh
    ierr = xf_Error(xf_ReadAsciiFile(BackgroundMesh, KeyValue, All->Mesh));
    if (ierr != xf_OK) return ierr;
    
    // Read embedded boundary group info
    return xf_Error(xf_NOT_SUPPORTED);
    //ierr = xf_Error(xf_ReadEbgFile(InputFile, All->Mesh));
    if (ierr != xf_OK) return ierr;
  }
  else{
    xf_printf("Unrecognized extension on InputFile: %s\n", pext);
    return xf_Error(xf_INPUT_ERROR);
  }
  
  /*----------------*/
  /* Set Parameters */
  /*----------------*/
  
  /* Set default parameters: if DefaultFlag is true, all are
   overwritten to default. */
  ierr = xf_AddKeyValueList(&All->Param->KeyValue, xf_DefaultParamList,
                            DefaultFlag, xfe_True);
  if ((ierr != xf_OK) && (ierr != xf_NOT_FOUND)) return xf_Error(ierr);
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadJobFileKeyValue
static int 
xf_ReadJobFileKeyValue( const char *jobFileIn, xf_KeyValue *pKeyValueLoc)
{
  // Only one processor should enter into this function
  int ierr, len;
  int myRank, nProc;
  
  char jobFile[xf_MAXSTRLEN];
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line0[xf_MAXLINELEN], *line;
  
  FILE *fjob;
  
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // only allow root to enter here
  if (myRank != 0) return xf_Error(xf_PARALLEL_ERROR);
  
  strcpy(jobFile, jobFileIn);
  
  /* Check for .job extension.  If not present, append. */
  len = strlen(jobFile);
  if (len < 4) 
    strcat(jobFile, ".job");
  else if( (strcmp(&jobFile[len-4],".job")) != 0 )
    strcat(jobFile, ".job");
  
  /* Open the .job file */
  fjob = fopen(jobFile, "r");
  if( fjob == NULL ){ 
    xf_printf("\nJob file '%s' not found!\n", jobFile);
    xf_printf(".job extension is optional.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  
  /*-----------------------------------------*/
  /*  Read Key-Value pairs in the .job file  */
  /*-----------------------------------------*/
  
  do{
    
    /* read line of job file */
    if (fgets(line0, xf_MAXLINELEN, fjob) == NULL)
      continue;
    line = line0; // set pointer
    
    /* if blank line or comment, continue to next line */
    if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;
    
    /* read key and value from line */
    ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
    if (ierr != xf_OK) return ierr;
    
    /* Each line that is not blank or a commment must be of the form
     key = value, where key and value are strings.  Neither key nor
     value may contain the delimeter, "=".  They may, however,
     contain spaces, " ".*/
    
    /* Store key and value in local list */
    ierr = xf_Error(xf_AddKeyValue(pKeyValueLoc, key, value, xfe_True));
    if (ierr == xf_OVERWROTE){
      xf_printf("Error. The key %s is assigned more than once in the .job file.\n", key);
      return xf_FILE_READ_ERROR;
    }
    if (ierr != xf_OK) return ierr;
    
    
  } while(feof(fjob) == 0);
  
  /* Done with job file */
  fclose(fjob);
  
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadAllFromJobFile
int 
xf_ReadAllFromJobFile( const char *jobFileIn, enum xfe_Bool ReadEqnSet,
                      xf_All **pAll)
{
  int ierr, terr;
  int myRank, nProc;
  int len, nData, iData;
  enum xfe_Bool ParallelFlag;
  enum xfe_Bool RestartFlag, DefaultFlag;
  xf_KeyValue KeyValueLoc;
  
  char value[xf_MAXSTRLEN];
  char InputFile[xf_MAXSTRLEN];
  char EqnSetFile[xf_MAXSTRLEN];
  char MeshMotionFile[xf_MAXSTRLEN];
  char GeometryFile[xf_MAXSTRLEN];
  char DataMerge[xf_MAXSTRLEN] = "None";
  char **DataFiles;
  
  xf_DataSet *DataSet;
  xf_EqnSet *EqnSet;
  xf_All *All = NULL;
  
  FILE *fid;
  
  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ParallelFlag = (nProc > 1);
  
  ierr = xf_Error(xf_InitKeyValue(&KeyValueLoc));
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0) // root reads the jobFile
    terr = xf_Error(xf_ReadJobFileKeyValue(jobFileIn, &KeyValueLoc));
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr); // all procs exit due to file read error
  
  /* Parallelize KeyValueLoc so every proc gets a copy */
  ierr = xf_Error(xf_ParallelizeKeyValue(&KeyValueLoc));
  if (ierr != xf_OK) return ierr;
  
  
  /*-------------------------------------*/
  /* Check for necessary run information */
  /*-------------------------------------*/
  
  ierr = xf_GetKeyValue(KeyValueLoc, "InputFile", InputFile);
  if( ierr != xf_OK ){
    xf_printf("Error reading job file!\n");
    xf_printf("InputFile not set.  Add following line to job file:\n\n");
    xf_printf("InputFile = filename.abc \n\n");
    xf_printf("where filename.abc is your input file.\n");
    xf_printf("An extension is required.\n\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  ierr = xf_GetKeyValue(KeyValueLoc, "EqnSetFile", value);
  if( ierr != xf_OK ){
    xf_printf("Error reading job file!\n");
    xf_printf("EqnSetFile not set.  Add following line to job file:\n\n");
    xf_printf("EqnSetFile = filename\n\n");
    xf_printf("where filename is your Equation Set file (extension required).\n");
    xf_printf("To use the EqnSet in the input .xfa file, specify\n");
    xf_printf("EqnSetFile = None or NULL.\n\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  ierr = xf_GetKeyValue(KeyValueLoc, "SavePrefix", value);
  if (ierr != xf_OK){
    xf_printf("Error reading job file.  Please specify SavePrefix.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  ierr = xf_GetKeyValueBool(KeyValueLoc, "Restart", &RestartFlag);
  if( ierr != xf_OK ){
    xf_printf("Error reading job file!\n");
    xf_printf("Restart parameter not set.  Add following line to job file:\n");
    xf_printf("\n");
    xf_printf("Restart = True/False\n");
    xf_printf("\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  DefaultFlag = xfe_False;
  ierr = xf_GetKeyValue(KeyValueLoc, "Default", value);
  if( ierr == xf_OK ){
    ierr = xf_Error(xf_Value2Enum(value, xfe_BoolName, xfe_BoolLast, 
                                  (int *) &DefaultFlag));
    if( ierr != xf_OK ) return ierr;
  }
  
  /* If restarting, print message about how parameters are set */
  if (RestartFlag == xfe_True){ 
    if( DefaultFlag == xfe_True ){
      xf_printf("\nDefault has been set to True in the job file.\n");
      xf_printf("Parameters from the xfa file will be ignored.\n");
      xf_printf("Parameters not specified in the job file will take on\n");
      xf_printf("default values.\n\n");
    }
    else{
      xf_printf("\nOrder of preference for parameters (lowest to highest):\n");
      xf_printf("- Initially, parameters are set to default values.\n");
      xf_printf("- Next, all parameters in the .xfa file take precedence\n");
      xf_printf("- Finally, all parameters in the job file take precedence\n\n");
    }
  }
  
  
  /*----------------*/
  /* Read InputFile */
  /*----------------*/
  
  ierr = xf_Error(xf_ReadAllInputFile(InputFile, &KeyValueLoc, DefaultFlag, pAll));
  if (ierr != xf_OK) return ierr;
  All = (*pAll);
  
  
  
  /* Add local parameters from .job file */
  ierr = xf_MergeKeyValue(&All->Param->KeyValue, KeyValueLoc, 2);
  if ((ierr != xf_OK) && (ierr != xf_OVERWROTE)) return xf_Error(ierr);
  
  
  /*---------------------*/
  /* Check All structure */
  /*---------------------*/
  
  All->EqnSet->Dim = All->Mesh->Dim; // Set EqnSet dimension
  
  ierr = xf_Error(xf_CheckAll(All));
  if (ierr != xf_OK) return ierr;
  
  /* Release local memory */
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValueLoc));
  if (ierr!=xf_OK) return ierr;
  
  
  /*-------------------------------*/
  /* Read any additional data sets */
  /*-------------------------------*/
  
  ierr = xf_GetKeyValue(All->Param->KeyValue, "DataMerge", DataMerge);
  if (ierr != xf_OK) return ierr;
  
  if (xf_NotNull(DataMerge)){
    
    ierr = xf_Error(xf_ScanXStringAlloc(DataMerge, xf_MAXSTRLEN, &nData,
                                        &DataFiles));
    if (ierr != xf_OK) return ierr;
    
    for (iData=0; iData<nData; iData++){
      
      ierr = xf_Error(xf_CreateDataSet(&DataSet));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, DataFiles[iData], DataSet));
      if (ierr!=xf_OK) return ierr;
      
      ierr = xf_Error(xf_DataSetMerge(DataSet, All->DataSet));
      if (ierr != xf_OK) return ierr;
      
      // destroy dataset
      ierr = xf_Error(xf_DestroyDataSet(DataSet));
      if (ierr != xf_OK) return ierr;
    }
    xf_Release2( (void **) DataFiles);
  }

  /*----------------------------*/
  /*  Read EqnSet if requested  */
  /*----------------------------*/
 
  if (ReadEqnSet){
    
    ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "EqnSetFile", EqnSetFile));
    if (ierr != xf_OK) return ierr;
    
    // Read EqnSetFile
    if (xf_NotNull(EqnSetFile)){
      
      ierr = xf_Error(xf_CreateEqnSet(&EqnSet));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ReadEqnSetFile(EqnSetFile, NULL, EqnSet));
      if (ierr != xf_OK) return ierr;
      
      /* All->EqnSet is updated with values (params, ResTerms, ICs, BCs,
       Outputs) from EqnSet. */
      ierr = xf_Error(xf_MergeEqnSet(All->EqnSet, EqnSet));
      if (ierr != xf_OK) return ierr;
      
      /* No longer need the EqnSet read in from the .eqn file */
      ierr = xf_Error(xf_DestroyEqnSet(EqnSet, xfe_True));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  
  /*--------------------------------*/
  /*  Read MeshMotion if requested  */
  /*--------------------------------*/

  // Check if MeshMotion read if requested
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "MeshMotionFile", MeshMotionFile));
  if (ierr != xf_OK) return ierr;
  
  // Read MeshMotionFile
  if (xf_NotNull(MeshMotionFile)){
    
    ierr = xf_Error(xf_CreateMeshMotion(&All->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadMeshMotionFile(MeshMotionFile, NULL, All->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
  }

  /*------------------------------*/
  /*  Read Geometry if requested  */
  /*------------------------------*/

  // Check if geometry read if requested
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "GeometryFile", GeometryFile));
  if (ierr != xf_OK) return ierr;
  
  // Read GeometryFile
  if (xf_NotNull(GeometryFile)){
    ierr = xf_Error(xf_ReadGeomFile(GeometryFile, NULL, All->Geom));
    if (ierr != xf_OK) return ierr;
  }


  /*----------------*/
  /* Set EqnSet Dim */
  /*----------------*/
  
  All->EqnSet->Dim = All->Mesh->Dim;
  
  return xf_OK;
}
