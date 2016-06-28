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
  FILE:  xf_MeshMotion.c

  This file contains functions for input/output and parallelization of
  the mesh Motion structure.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_MeshMotionAnalytical.h"



/*-----------------*/
/* Parallelization */
/*-----------------*/  

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeAnaMotions
static int 
xf_ParallelizeAnaMotions( xf_AnaMotions *AnaMotions){
  
  int ierr, i;
  int myRank, nProc;
  xf_AnaMotion *pAnaMotion;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_MPI_Bcast((void *) &AnaMotions->nTerm, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  if (myRank > 0){ // allocate on non-root procs
    ierr = xf_Error(xf_AllocAnaMotions(AnaMotions, AnaMotions->nTerm));
    if (ierr != xf_OK) return ierr;
  }

  for (i=0; i<AnaMotions->nTerm; i++){
    
    pAnaMotion = AnaMotions->AnaMotion+i;

    // Name
    ierr = xf_Error(xf_ParallelizeString(&pAnaMotion->Name));
    if (ierr != xf_OK) return ierr;

    // MotionType
    ierr = xf_Error(xf_MPI_Bcast((void *) &pAnaMotion->MotionType, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;

    // key-values
    ierr = xf_Error(xf_ParallelizeKeyValue(&pAnaMotion->MotionKeyValue));
    if (ierr != xf_OK) return ierr;

  } // i

  // BlendType
  ierr = xf_Error(xf_MPI_Bcast((void *) &AnaMotions->BlendType, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;


  // key values
  ierr = xf_Error(xf_ParallelizeKeyValue(&AnaMotions->BlendKeyValue));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeAnaMotionsSet
static int 
xf_ParallelizeAnaMotionsSet( xf_AnaMotionsSet *AnaMotionsSet){
  
  int ierr, i;
  int myRank, nProc;
  xf_AnaMotion *pAnaMotion;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_MPI_Bcast((void *) &AnaMotionsSet->nMotion, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  if (myRank > 0){ // allocate on non-root procs
    ierr = xf_Error(xf_AllocAnaMotionsSet(AnaMotionsSet, AnaMotionsSet->nMotion));
    if (ierr != xf_OK) return ierr;
  }

  for (i=0; i<AnaMotionsSet->nMotion; i++){

    ierr = xf_Error(xf_ParallelizeAnaMotions(AnaMotionsSet->AnaMotions + i));
    if (ierr != xf_OK) return ierr;

  } // i
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeMeshMotion
int 
xf_ParallelizeMeshMotion( xf_MeshMotion *Motion)
{
  int ierr;
  int myRank, nProc;

  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if (nProc == 1) return xf_OK; // nothing to do

  // broadcast type
  ierr = xf_Error(xf_MPI_Bcast((void *) &Motion->Type, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  // Parallelize different types of mesh motion data
  switch (Motion->Type){
  case xfe_Motion_Analytical:
    if (myRank > 0){ // create analytical motions set on non-root procs
      ierr = xf_Error(xf_CreateAnaMotionsSet((xf_AnaMotionsSet **) &Motion->Data));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_ParallelizeAnaMotionsSet((xf_AnaMotionsSet *) Motion->Data));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}


/*-----------------*/
/* Text File Input */
/*-----------------*/ 


/******************************************************************/
//   FUNCTION Definition: xf_ReadMeshMotionAnalytical
static int 
xf_ReadMeshMotionAnalytical(FILE *fmm, char **InputStrings, int *piString,
			    xf_AnaMotionsSet **pAnaMotionsSet)
{
  int ierr, nTerm, iTerm, i;
  int nMotion, iMotion;
  char line0[xf_MAXLINELEN], *line, longvalue[xf_MAXLINELEN];
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  enum xfe_Bool ReachedEnd = xfe_False;
  xf_AnaMotionsSet *AnaMotionsSet = NULL;
  xf_AnaMotions *AnaMotions = NULL;
  xf_KeyValue *pKeyValue = NULL;

  ierr = xf_Error(xf_CreateAnaMotionsSet(pAnaMotionsSet));
  if (ierr != xf_OK) return ierr;

  AnaMotionsSet = (*pAnaMotionsSet);

  nMotion =  0;
  iMotion = -1;
  nTerm =  0;
  iTerm = -1;
  do{      
    /* read line of file */
    ierr = xf_LineFromFileOrStrings(fmm, InputStrings, piString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;

    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line
    
    /* Check if reached end of block */
    if (strncmp(line, "ENDBLOCK", 8) == 0){
      ReachedEnd = xfe_True; break;
    }
    if (strncmp(line, "STARTBLOCK", 10) == 0) break;
    
    if (strncmp(line, "nMotion", 7) == 0){

      if (iTerm != (nTerm-1)){
	xf_printf("Number of Analytical terms is less than nTerm.\n");
	return xf_Error(xf_FILE_READ_ERROR);
      }

      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%d", &nMotion) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (iMotion != -1){
	xf_printf("Lines in ANALYTICAL block out of order.\n");
	return xf_Error(xf_FILE_READ_ERROR);
      }
      
      ierr = xf_Error(xf_AllocAnaMotionsSet(AnaMotionsSet, nMotion));
      if (ierr != xf_OK) return ierr;

      continue;
    }

    if (strncmp(line, "nTerm", 5) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;
      
      if (sscanf(value, "%d", &nTerm) != 1) return xf_Error(xf_FILE_READ_ERROR);
      if (nMotion == 0){
	nMotion = 1; // nMotion flag is optional, hence this
	ierr = xf_Error(xf_AllocAnaMotionsSet(AnaMotionsSet, nMotion));
	if (ierr != xf_OK) return ierr;
      }
      iMotion++;
      iTerm = -1;
      AnaMotions = AnaMotionsSet->AnaMotions+iMotion;

      ierr = xf_Error(xf_AllocAnaMotions(AnaMotions, nTerm));
      if (ierr != xf_OK) return ierr;

      continue;
    }

    // sanity check: we need nMotion to be positive by now
    if (nMotion <= 0){
      xf_printf("Error reading ANALYTICAL block.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }

    if (nTerm < 0){
      xf_printf("Error reading ANALYTICAL block.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }

    if (strncmp(line, "Name", 4) == 0){
      iTerm++;  
      if (iTerm >= nTerm){
	xf_printf("Number of MotionNames exceeds nTerm.\n");
	return xf_Error(xf_FILE_READ_ERROR);
      }

      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_AllocString(&AnaMotions->AnaMotion[iTerm].Name, xf_MAXSTRLEN, value));
      if (ierr != xf_OK) return ierr;

      continue;
    }

    if (iTerm < 0){
      xf_printf("Unrecognized line in ANALYTICAL block before first MotionName.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }

    
    if (strncmp(line, "MotionType", 10) == 0){
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_Value2Enum(value, xfe_AnaMotionName, xfe_AnaMotionLast, 
				    (int *) &AnaMotions->AnaMotion[iTerm].MotionType));
      if( ierr != xf_OK ) return ierr;
      
      // any following key values will be added to MotionKeyValue
      pKeyValue = &AnaMotions->AnaMotion[iTerm].MotionKeyValue;

      continue;
    }

    if (strncmp(line, "BlendType", 9) == 0){

      if (iTerm != (nTerm-1)){
	xf_printf("Only one BlendType allowed; place after all terms in a motion.\n");
	return xf_Error(xf_FILE_READ_ERROR);
      }
	
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_Value2Enum(value, xfe_AnaBlendName, xfe_AnaBlendLast, 
				    (int *) &AnaMotions->BlendType));
      if( ierr != xf_OK ) return ierr;
      
      // any following key values will be added to BlendKeyValue
      pKeyValue = &AnaMotions->BlendKeyValue;

      continue;
    }

    if (pKeyValue == NULL) return xf_Error(xf_FILE_READ_ERROR);


    /* At this point, reading a key=value pair */
    ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
    if (ierr != xf_OK) return ierr;
      
    /* Store key and value in local list */
    ierr = xf_Error(xf_AddKeyValue(pKeyValue, key, value, xfe_True));
    if (ierr == xf_OVERWROTE){
      xf_printf("Error. Key %s is assigned more than once in ANALYTICAL term %d.\n", 
		key, iTerm);
      return xf_FILE_READ_ERROR;
    }
    if (ierr != xf_OK) return ierr;
    
  } while (1);

  if (iTerm != (nTerm-1)){
    xf_printf("Number of Analytical terms is less than nTerm.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }

  if (!ReachedEnd){
    xf_printf("Error. ANALYTICAL block not terminated with ENDBLOCK\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadMeshMotionFileSerial
static int 
xf_ReadMeshMotionFileSerial(char *MotionFile, char **InputStrings, 
			    xf_MeshMotion *Motion)
{
  int ierr, iString = 0;
  enum xfe_Bool done;
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line0[200], *line;
  FILE *fmm = NULL;

  /* Check input */
  if (((MotionFile == NULL) && (InputStrings == NULL)) || 
      ((MotionFile != NULL) && (InputStrings != NULL))) 
    return xf_Error(xf_INPUT_ERROR);

  if (MotionFile != NULL){
    /* Open file */
    if ((fmm = fopen(MotionFile, "r")) == NULL){
      xf_printf("Could not find Mesh Motion file: %s\n", MotionFile);
      xf_printf("Make sure the appropriate extension is included.\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
  }
  else{
    // Start at the first string
    iString = 0;
  }
  
  /* Read blocks */
  done = xfe_False;
  do{
    ierr = xf_LineFromFileOrStrings(fmm, InputStrings, &iString, line0, &line);
    if (ierr == xf_END_OF_FILE) break;
    if (ierr != xf_OK) return ierr;

    if (xf_TrimAndCheckBlank(&line, 200)) continue; // blank or comment line


    if (strncmp(line, "STARTBLOCK", 10) == 0){
      /* read in second token */
      ierr = xf_Error(xf_DesiredToken(line, " ", 1, value));
      if (ierr != xf_OK) return ierr;

      if (strncmp(value, "ANALYTICAL", 5) == 0){
	/* Read Analytical block */
	ierr = xf_Error(xf_ReadMeshMotionAnalytical(fmm, InputStrings, &iString, 
						    (xf_AnaMotionsSet **) &Motion->Data));
	if (ierr != xf_OK) return ierr;
	Motion->Type = xfe_Motion_Analytical;
      }
      else return xf_Error(xf_NOT_SUPPORTED);
    }

  } while (!done);

  if (fmm != NULL) fclose(fmm);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadMeshMotionFile
int 
xf_ReadMeshMotionFile(char *MotionFile, char **InputStrings, 
		      xf_MeshMotion *Motion)
{
  int ierr, terr, myRank;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root reads in mesh motion file
  if (myRank == 0)
    terr = xf_Error(xf_ReadMeshMotionFileSerial(MotionFile, InputStrings, Motion));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  // parallelize Mesh Motion structure
  ierr = xf_Error(xf_ParallelizeMeshMotion(Motion));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/*---------------------*/
/* Binary Input/Output */
/*---------------------*/ 

/******************************************************************/
//   FUNCTION Definition: xf_WriteAnaMotionBinary
static int 
xf_WriteAnaMotionBinary( xf_AnaMotion *AnaMotion, FILE *fid)
{
  int ierr, rev;
  
  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // Name
  ierr = xf_Error(xf_WriteStringBinary(AnaMotion->Name, fid));
  if (ierr != xf_OK) return ierr;

  // Motion Type
  ierr = xf_Error(xf_WriteStringBinary(xfe_AnaMotionName[AnaMotion->MotionType], fid));
  if (ierr != xf_OK) return ierr;

  // MotionKeyValue
  ierr = xf_Error(xf_WriteKeyValueBinary(AnaMotion->MotionKeyValue, fid));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadAnaMotionBinary
static int 
xf_ReadAnaMotionBinary( FILE *fid, xf_AnaMotion *AnaMotion)
{
  int ierr, rev;
  
  rev = 0;  // writer revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);

  // Name
  ierr = xf_Error(xf_ReadStringBinary(fid, -1, NULL, &AnaMotion->Name));
  if (ierr != xf_OK) return ierr;

  // MotionType
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_AnaMotionName, xfe_AnaMotionLast, 
				    (int *) &AnaMotion->MotionType));
  if (ierr != xf_OK) return ierr;

  // MotionKeyValue
  ierr = xf_Error(xf_ReadKeyValueBinary(fid, &AnaMotion->MotionKeyValue));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteAnaMotionsBinary
static int 
xf_WriteAnaMotionsBinary( xf_AnaMotions *AnaMotions, FILE *fid){
  int ierr, rev, i;

  rev = 0;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // nTerm
  if (fwrite(&AnaMotions->nTerm, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // AnaMotion
  for (i=0; i<AnaMotions->nTerm; i++){
    ierr = xf_Error(xf_WriteAnaMotionBinary(AnaMotions->AnaMotion + i, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  // Blend Type
  ierr = xf_Error(xf_WriteStringBinary(xfe_AnaBlendName[AnaMotions->BlendType], fid));
  if (ierr != xf_OK) return ierr;

  // BlendKeyValue
  ierr = xf_Error(xf_WriteKeyValueBinary(AnaMotions->BlendKeyValue, fid));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadAnaMotionsBinary
static int 
xf_ReadAnaMotionsBinary( FILE *fid, xf_AnaMotions *AnaMotions){
  int ierr, rev, i;
  
  rev = 0;  // writer revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);

  // nTerm
  if (fread(&AnaMotions->nTerm, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_READ_ERROR);
  
  // Allocate AnaMotion
  ierr = xf_Error(xf_AllocAnaMotions(AnaMotions, AnaMotions->nTerm));
  if (ierr != xf_OK) return ierr;

  // Read AnaMotion
  for (i=0; i<AnaMotions->nTerm; i++){
    ierr = xf_Error(xf_ReadAnaMotionBinary(fid, AnaMotions->AnaMotion + i));
    if (ierr != xf_OK) return ierr;
  }

  // BlendType
  ierr = xf_Error(xf_ReadEnumBinary(fid, xfe_AnaBlendName, xfe_AnaBlendLast, 
				    (int *) &AnaMotions->BlendType));
  if (ierr != xf_OK) return ierr;

  // BlendKeyValue
  ierr = xf_Error(xf_ReadKeyValueBinary(fid, &AnaMotions->BlendKeyValue));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteAnaMotionsSetBinary
static int 
xf_WriteAnaMotionsSetBinary( xf_AnaMotionsSet *AnaMotionsSet, FILE *fid){
  int ierr, rev, i;

  rev = 1;  // writer revision number
  if (fwrite(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // nMotion
  if (fwrite(&AnaMotionsSet->nMotion, sizeof(int), 1, fid) != 1) 
    return xf_Error(xf_FILE_WRITE_ERROR);

  // AnaMotions
  for (i=0; i<AnaMotionsSet->nMotion; i++){
    ierr = xf_Error(xf_WriteAnaMotionsBinary(AnaMotionsSet->AnaMotions + i, fid));
    if (ierr != xf_OK) return ierr;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadAnaMotionsSetBinary
static int 
xf_ReadAnaMotionsSetBinary( FILE *fid, xf_AnaMotionsSet *AnaMotionsSet){
  int ierr, rev, i;
  
  rev = 1;  // writer revision number
  if (fread(&rev, sizeof(int), 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev > 1) return xf_Error(xf_FILE_READ_ERROR);

  // nMotion
  if (rev >= 1){
    if (fread(&AnaMotionsSet->nMotion, sizeof(int), 1, fid) != 1) 
      return xf_Error(xf_FILE_READ_ERROR);
  }
  else AnaMotionsSet->nMotion = 1;
  
  // Allocate AnaMotionsSet
  ierr = xf_Error(xf_AllocAnaMotionsSet(AnaMotionsSet, AnaMotionsSet->nMotion));
  if (ierr != xf_OK) return ierr;

  // Read AnaMotions
  for (i=0; i<AnaMotionsSet->nMotion; i++){
    ierr = xf_Error(xf_ReadAnaMotionsBinary(fid, AnaMotionsSet->AnaMotions + i));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteMeshMotionBinarySerial
int 
xf_WriteMeshMotionBinarySerial( xf_MeshMotion *Motion, FILE *fid)
{
  int ierr, rev, si, i;
  enum xfe_Bool flag;

  si = sizeof(int);

  rev = 0;  // writer revision number
  if (fwrite(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);

  // Type
  if (fwrite(&Motion->Type, si, 1, fid) != 1) return xf_Error(xf_FILE_WRITE_ERROR);
  
  // Write different types of mesh motion data
  switch (Motion->Type){
  case xfe_Motion_Analytical:
    ierr = xf_Error(xf_WriteAnaMotionsSetBinary((xf_AnaMotionsSet *) Motion->Data, fid));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadMeshMotionBinarySerial
int 
xf_ReadMeshMotionBinarySerial(  FILE *fid, xf_MeshMotion *Motion)
{
  int ierr, rev, si, i;
  enum xfe_Bool flag;

  si = sizeof(int);

  rev = 0;  // writer revision number
  if (fread(&rev, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  if (rev != 0) return xf_Error(xf_FILE_READ_ERROR);
  
  // Type
  if (fread(&Motion->Type, si, 1, fid) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // Read different types of mesh motion data
  switch (Motion->Type){
  case xfe_Motion_Analytical:
    ierr = xf_Error(xf_CreateAnaMotionsSet((xf_AnaMotionsSet **)&Motion->Data));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadAnaMotionsSetBinary(fid, (xf_AnaMotionsSet *) Motion->Data));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_WriteMeshMotionBinary
int 
xf_WriteMeshMotionBinary(xf_MeshMotion *Motion, FILE *fid)
{
  int ierr, terr, myRank;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root writes
  if (myRank == 0)
    terr = xf_Error(xf_WriteMeshMotionBinarySerial(Motion, fid));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ReadMeshMotionBinary
int 
xf_ReadMeshMotionBinary(  FILE *fid, xf_MeshMotion *Motion)
{
  int ierr, terr, myRank;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root reads 
  if (myRank == 0)
    terr = xf_Error(xf_ReadMeshMotionBinarySerial(fid, Motion));

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  // parallelize Motion
  ierr = xf_Error(xf_ParallelizeMeshMotion(Motion));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}






