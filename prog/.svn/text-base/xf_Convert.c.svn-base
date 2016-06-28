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
 FILE:  xf_Convert.c
 
 This program converts to/from various mesh formats.
 
 */

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_Arg.h"


/******************************************************************/
//  FUNCTION Definition: xf_TetGen2Gri
static int 
xf_TetGen2Gri(const char *inFile, const char *outFile, const char *bnames)
{
  int ierr, k, len, inode, ibfgrp, nbfgrp;
  int nnode, nelem, elem, i0, i1, i2, i3, i4;
  int **B, bmin, bmax, bcur, ibface, nbface, nbfacetot;
  int nBFGTitle = 0;
  char **BFGTitle = NULL;
  char inRoot[xf_MAXSTRLEN], Title[xf_MAXSTRLEN];
  char snode[xf_MAXSTRLEN], sele[xf_MAXSTRLEN], sface[xf_MAXSTRLEN];
  real r0, r1, r2;
  FILE *fnode, *fele, *fface, *fgri;
  
  /* Determine inRoot */
  strcpy(inRoot, inFile);
  len = strlen(inFile);
  for (k=len-4; k<len; k++) inRoot[k] = '\0';
  
  /* Determine if any boundary names are given */
  if ((bnames != NULL) && (xf_NotNull(bnames))){
    ierr = xf_Error(xf_ScanXStringAlloc(bnames, xf_MAXSTRLEN, &nBFGTitle,
                                        &BFGTitle));
    if (ierr != xf_OK) return ierr;
  }
  
  /* Convert (root.node, root.ele, root.face) to root.gri */
  sprintf(snode, "%s.node", inRoot);
  sprintf(sele,  "%s.ele",  inRoot);
  sprintf(sface, "%s.face", inRoot);
  
  if ((fnode=fopen(snode,"r"))==NULL) return xf_Error(xf_FILE_READ_ERROR);
  if ((fele =fopen(sele ,"r"))==NULL) return xf_Error(xf_FILE_READ_ERROR);
  if ((fface=fopen(sface,"r"))==NULL) return xf_Error(xf_FILE_READ_ERROR);
  if ((fgri =fopen(outFile ,"w"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);
  
  ierr = fscanf(fnode, "%d %d %d %d", &nnode, &i0, &i1, &i2);
  if (ierr != 4) return xf_Error(xf_FILE_READ_ERROR);
  
  ierr = fscanf(fele, "%d %d %d", &nelem, &i0, &i1);
  if (ierr != 3) return xf_Error(xf_FILE_READ_ERROR);
  
  fprintf(fgri, "%d %d 3\n", nnode, nelem);
  for (inode=0; inode<nnode; inode++){
    ierr = fscanf(fnode, "%d %lf %lf %lf", &i0, &r0, &r1, &r2);
    if (ierr != 4) return xf_Error(xf_FILE_READ_ERROR);
    fprintf(fgri, "%.15E %.15E %.15E\n", r0, r1, r2);
  }
  
  ierr = fscanf(fface, "%d %d", &nbfacetot, &i0);
  if (ierr != 2) return xf_Error(xf_FILE_READ_ERROR);
  
  ierr = xf_Error(xf_Alloc2((void ***) &B, nbfacetot, 4, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (ibface=0; ibface<nbfacetot; ibface++){
    ierr = fscanf(fface, "%d %d %d %d %d", &i0, 
                  B[ibface]+0, B[ibface]+1, B[ibface]+2, B[ibface]+3);
    if (ierr != 5) return xf_Error(xf_FILE_READ_ERROR);
    bcur = B[ibface][3];
    if (bcur < 0) bcur = -bcur;
    B[ibface][3] = bcur;
    if (ibface == 0)
      bmin = bmax = bcur;
    else{
      bmin = min(bmin, bcur);
      bmax = max(bmax, bcur);
    }
  } // ibface
  
  xf_printf("bmin = %d, bmax = %d\n",bmin,bmax);
  nbfgrp = bmax-bmin+1;
  
  fprintf(fgri, "%d\n", nbfgrp);
  
  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
    nbface = 0;
    for (ibface=0; ibface<nbfacetot; ibface++)
      if (B[ibface][3] == bmin+ibfgrp) nbface++;
    
    if ((BFGTitle != NULL) && (ibfgrp < nBFGTitle))
      strcpy(Title, BFGTitle[ibfgrp]);
    else
      sprintf(Title, "Boundary_%d", ibfgrp);
    
    fprintf(fgri, "%d 3 %s\n", nbface, Title);
    for (ibface=0; ibface<nbfacetot; ibface++)
      if (B[ibface][3] == bmin+ibfgrp) 
        fprintf(fgri, "%d %d %d\n", B[ibface][0], B[ibface][1], B[ibface][2]);
  } // ibfgrp
  
  
  // Release memory
  xf_Release2( (void **) B);
  xf_Release2( (void **) BFGTitle);
  
  fprintf(fgri, "%d 1 TetLagrange\n", nelem);
  for (elem=0; elem<nelem; elem++){
    ierr = fscanf(fele, "%d %d %d %d %d", &i0, &i1, &i2, &i3, &i4); 
    if (ierr != 5) return xf_Error(xf_FILE_READ_ERROR);
    fprintf(fgri, "%d %d %d %d\n", i1, i2, i3, i4);
  } // elem
  
  
  fclose(fnode);
  fclose(fele);
  fclose(fface);
  fclose(fgri);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, len;
  int inExtLen, outExtLen;
  char *ArgIn[] = {"in", "NULL", "input mesh file (with extension)",
    "out", "NULL", "output mesh file (with extension)",
    "bnames", "NULL", "space-separated boundary names",
  "\0"};
  char inFile[xf_MAXSTRLEN];  
  char outFile[xf_MAXSTRLEN];
  char bnames[xf_MAXSTRLEN];
  char *inExt, *outExt;
  xf_KeyValue KeyValue;
  xf_All *All;
  xf_Mesh *Mesh;
  
  xf_printf("\n");
  xf_printf("=== xf_Convert: Mesh Conversion  ===\n");
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
  
  // Get bnames
  ierr = xf_GetKeyValue(KeyValue, "bnames", bnames);
  if (ierr != xf_OK) return ierr;
  
  // extensions are required
  if ((len = strlen(inFile)) <= 0) return xf_Error(xf_FILE_READ_ERROR);
  for (inExtLen=0; ((inExtLen<len)&&(inFile[len-inExtLen-1]!='.')); inExtLen++);
  inExt = inFile + len - inExtLen; 

  if ((len = strlen(outFile)) <= 0) return xf_Error(xf_FILE_READ_ERROR);
  for (outExtLen=0; ((outExtLen<len)&&(outFile[len-outExtLen-1]!='.')); outExtLen++);
  outExt = outFile + len - outExtLen; 
  
  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;
  
  /* Perform the conversion */
  
  if ( ( inExtLen == 3) && (strncmp( inExt, "tgn", 3) == 0) &&
       (outExtLen == 3) && (strncmp(outExt, "gri", 3) == 0)){
    // TetGen to gri can be performed in one shot
    ierr = xf_Error( xf_TetGen2Gri(inFile, outFile, bnames) );
    if (ierr != xf_OK) return ierr;
  }
  else{
    // All other conversions are through .xfa

    /* Create .xfa structure */
    ierr = xf_Error(xf_CreateAll(&All, xfe_False));
    if (ierr != xf_OK) return ierr;

    Mesh = All->Mesh;

    /* Read in using extension to identify the format */

    if ((inExtLen == 3) && (strncmp( inExt, "xfa", 3) == 0)){
      // binary .xfa format
      ierr = xf_Error(xf_ReadAllBinary(inFile, All));
      if (ierr!=xf_OK) return ierr;
    }
    else if ((inExtLen == 3) && (strncmp( inExt, "msh", 3) == 0)){
      // Fluent/Gambit format
      ierr = xf_Error( xf_ReadFmshFile(inFile, Mesh) );
      if (ierr != xf_OK) return ierr;
    }
    else if ((inExtLen == 3) && (strncmp( inExt, "gri", 3) == 0)){
      // xflow ASCII format
      ierr = xf_Error( xf_ReadGriFile(inFile, NULL, Mesh) );
      if (ierr != xf_OK) return ierr;
    }
    else if ((inExtLen == 4) && (strncmp( inExt, "gmsh", 4) == 0)){
      // Gmsh format
      ierr = xf_Error( xf_ReadGmshFile(inFile, Mesh) );
      if (ierr != xf_OK) return ierr;
    }
    else if ((inExtLen == 4) && (strncmp( inExt, "bamg", 4) == 0)){
      // BAMG format
      ierr = xf_Error( xf_ReadBamgFile(inFile, NULL, Mesh) );
      if (ierr != xf_OK) return ierr;
    } 
    else
      return xf_Error(xf_NOT_SUPPORTED);


    /* Write out using extension to identify the format */
    
    xf_printf("Writing %s.\n", outFile);

    if ((outExtLen == 3) && (strncmp( outExt, "xfa", 3) == 0)){
      // binary .xfa format
      ierr = xf_Error(xf_WriteAllBinary(All, outFile));
      if (ierr!=xf_OK) return ierr;
    }
    else if ((outExtLen == 3) && (strncmp( outExt, "msh", 3) == 0)){
      // Fluent/Gambit format
      return xf_Error(xf_NOT_SUPPORTED);
    }
    else if ((outExtLen == 3) && (strncmp( outExt, "gri", 3) == 0)){
      // xflow ASCII format
      ierr = xf_Error( xf_WriteGriFile(Mesh, outFile) );
      if (ierr != xf_OK) return ierr;
    }
    else if ((outExtLen == 4) && (strncmp( outExt, "gmsh", 4) == 0)){
      // Gmsh format
      ierr = xf_Error( xf_WriteGmshFile(Mesh, outFile) );
      if (ierr != xf_OK) return ierr;
    }
    else if ((outExtLen == 4) && (strncmp( outExt, "bamg", 4) == 0)){
      // BAMG format
      ierr = xf_Error( xf_WriteBamgFile(Mesh, xfe_False, NULL, outFile) );
      if (ierr != xf_OK) return ierr;
    }
    else
      return xf_Error(xf_NOT_SUPPORTED);
        
    /* Destroy .xfa structure */
    ierr = xf_Error(xf_DestroyAll(All));
    if (ierr!=xf_OK) return ierr;
  }
  
  xf_printf("xf_Convert finished.\n");
  
  return xf_OK;
}
