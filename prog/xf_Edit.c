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
  FILE:  xf_Edit.c

  This program allows interactive editing of a .xfa or .gri file

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
#include "xf_MeshTools.h"
#include "xf_Math.h"



/******************************************************************/
//  FUNCTION Definition: xf_RenameBoundaries
static int 
xf_RenameBoundaries(xf_All *All)
{
/*

PURPOSE: 

  Interactively renames boundaries in All
  
INPUTS:

  All : all structure

OUTPUTS:  

RETURNS: Error Code

*/
  int ierr, i, nbfgrp, ibfgrp;
  enum xfe_Bool done;
  char buf[xf_MAXSTRLEN];
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  done = xfe_False;

  while (!done){
    nbfgrp = Mesh->nBFaceGroup;
    
    xf_printf(" %8s %8s %s\n", "[Group]", "[nBFace]", "[Title]");
    for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
      xf_printf(" %8d %8d %s\n", ibfgrp, Mesh->BFaceGroup[ibfgrp].nBFace,
		Mesh->BFaceGroup[ibfgrp].Title);
    }

    xf_printf("Choose group number to rename (or -1 to finish)\n> ");

    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%d", &i) == 1) && (i >= -1) && (i < nbfgrp)){
      if (i == -1)
	done = xfe_True;
      else{
	xf_printf("Title of ibfgrp=%d is %s\n", i,  Mesh->BFaceGroup[i].Title);
	xf_printf("Choose new title:\n> ");
	if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	    (strlen(buf)>1)){
	  buf[strlen(buf)-1] = '\0'; // remove ending newline
	  xf_Release( (void *) Mesh->BFaceGroup[i].Title);
	  ierr = xf_Error(xf_AllocString(&Mesh->BFaceGroup[i].Title, xf_MAXSTRLEN, buf));
	  if (ierr != xf_OK) return ierr;
	}
	else{
	  xf_printf("Not understood.\n");
	}
      }
    }
    else
      xf_printf("Not understood.\n");

  } // end while ! done


  return xf_OK;
}


/******************************************************************/
//  FUNCTION Definition: xf_MergeBoundaries
static int 
xf_MergeBoundaries(xf_All *All)
{
/*

PURPOSE: 

  Interactively merges boundaries in All
  
INPUTS:

  All : all structure

OUTPUTS:  

RETURNS: Error Code

*/

  int ierr, i, nbfgrp, ibfgrp;
  int n, j, nBFaceTot;
  int ibface;
  int *vi = NULL;
  enum xfe_Bool done, valid;
  char buf[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  xf_BFace *BFaceNew = NULL;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  done = xfe_False;

  while (!done){
    nbfgrp = Mesh->nBFaceGroup;
    
    xf_printf(" %8s %8s %s\n", "[Group]", "[nBFace]", "[Title]");
    for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
      xf_printf(" %8d %8d %s\n", ibfgrp, Mesh->BFaceGroup[ibfgrp].nBFace,
		Mesh->BFaceGroup[ibfgrp].Title);
    }

    xf_printf("Choose group numbers to merge (or -1 to finish)\n> ");

    if (fgets(buf, xf_MAXSTRLEN, stdin) != NULL){
      vi = NULL;
      ierr = xf_Error(xf_ScanXIntAlloc(buf, &n, &vi));
      if (ierr != xf_OK) return ierr;
      
      if ((n == 1) && (vi[0] == -1)) done = xfe_True;

      if (n > 1){
	valid = xfe_True;
	for (i=0; i<n; i++) if ((vi[i]<0) || (vi[i]>=nbfgrp)) valid = xfe_False;

	if ((valid) && (n > 1)){
	  
	  xf_printf("Merging %d groups.\n", n);

	  // merge the n groups in vi
	  nBFaceTot = 0;
	  for (i=0; i<n; i++) nBFaceTot += Mesh->BFaceGroup[vi[i]].nBFace;

	  xf_printf("New total bfaces = %d\n", nBFaceTot);
	  
	  // Allocate new BFace
	  ierr = xf_Error(xf_Alloc( (void **) &BFaceNew, nBFaceTot, sizeof(xf_BFace)));
	  if (ierr != xf_OK) return ierr;

	  // Copy over bfaces from groups being merged
	  for (i=0,nBFaceTot=0; i<n; i++)
	    for (ibface=0; ibface<Mesh->BFaceGroup[vi[i]].nBFace; ibface++){
	      BFaceNew[nBFaceTot++] = Mesh->BFaceGroup[vi[i]].BFace[ibface];
	      xf_InitBFace(Mesh->BFaceGroup[vi[i]].BFace + ibface);
	    }

	  // choose title
	  strcpy(Title, Mesh->BFaceGroup[vi[0]].Title);
	  
	  xf_printf("Choose new title of merged group:\n> ");
	  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	      (strlen(buf)>1)){
	    buf[strlen(buf)-1] = '\0'; // remove ending newline
	    strcpy(Title, buf);
	  }
	  else{
	    xf_printf("Not understood, using Title = %s.\n", buf);
	  }
	  
	  // Destroy old bface groups
	  for (i=0; i<n; i++){
	    ierr = xf_Error(xf_DestroyBFaceGroup(Mesh->BFaceGroup + vi[i]));
	    if (ierr != xf_OK) return ierr;
	    Mesh->BFaceGroup[vi[i]].BFace = NULL;
	  }
	  
	  // New group is stored in position of min(vi[i]);
	  ibfgrp = vi[0];
	  for (i=1; i<n; i++) ibfgrp = min(ibfgrp, vi[i]);
	  ierr = xf_Error(xf_AllocString(&Mesh->BFaceGroup[ibfgrp].Title, xf_MAXSTRLEN, Title));
	  if (ierr != xf_OK) return ierr;
	  Mesh->BFaceGroup[ibfgrp].BFace  = BFaceNew;
	  Mesh->BFaceGroup[ibfgrp].nBFace = nBFaceTot;

	  // shift all other groups and reallocate .BFace
	  for (i=j=ibfgrp+1; i<nbfgrp; i++){
	    if (Mesh->BFaceGroup[i].BFace != NULL){ 
	      if (i != j) Mesh->BFaceGroup[j] = Mesh->BFaceGroup[i];
	      j++;
	    }
	  }
	  Mesh->nBFaceGroup = nbfgrp = nbfgrp-n+1; // new number of groups
	  xf_printf("New number of boundary groups = %d\n", nbfgrp);
	  ierr = xf_Error(xf_ReAlloc((void **) &Mesh->BFaceGroup, nbfgrp, sizeof(xf_BFaceGroup)));
	  if (ierr != xf_OK) return ierr;	  
	  
	}
	else xf_printf("Input error.\n");
	
      }
      else xf_printf("Not understood.\n");
      xf_Release( (void *) vi);

    }
    else xf_printf("Not understood.\n");
    
  } // end while ! done

  return xf_OK;
}

/******************************************************************/
//  FUNCTION Definition: xf_SplitBoundaries
static int 
xf_SplitBoundaries(xf_All *All)
{
/*

PURPOSE: 

  Interactively splits boundaries in All according to coordinates
  
INPUTS:

  All : all structure

OUTPUTS:  

RETURNS: Error Code

*/

  int ierr, i, nbfgrp, ibfgrp;
  int d, dim, inode, ibfgrpNew;
  int ibface, nbface1, nbface2;
  int n, nBFaceTot, nBFaceNew;
  int nfnode, egrp, elem, face;
  enum xfe_Bool done, valid;
  char buf[xf_MAXSTRLEN];
  char Title1[xf_MAXSTRLEN];
  char Title2[xf_MAXSTRLEN];
  int fvec[xf_MAXQ1NODE];
  int *bflag;
  real val, rval;
  xf_BFace BFace;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim = Mesh->Dim;

  done = xfe_False;

  while (!done){
    nbfgrp = Mesh->nBFaceGroup;
    
    xf_printf(" %8s %8s %s\n", "[Group]", "[nBFace]", "[Title]");
    for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
      xf_printf(" %8d %8d %s\n", ibfgrp, Mesh->BFaceGroup[ibfgrp].nBFace,
		Mesh->BFaceGroup[ibfgrp].Title);
    }

    xf_printf("Choose group number to split (or negative to finish)\n> ");

    valid = xfe_False;
    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%d", &ibfgrp) == 1) && (ibfgrp < nbfgrp))
      valid = xfe_True;
    if (!valid){ xf_printf("Invalid group number.\n"); continue;}
    
    if (ibfgrp < 0){
      done = xfe_True;
      continue;
    }
    
    xf_printf("Splitting group %d.\n", ibfgrp);
    
    xf_printf("Enter dimension for demarcation (0=x, 1=y, 2=z)\n> ");
    valid = xfe_False;
    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%d", &d) == 1) && (d >= 0) && (d < dim))
      valid = xfe_True;
    if (!valid){ xf_printf("Invalid dimension.\n"); continue;}

    xf_printf("Enter real value of demarcation\n> ");
    valid = xfe_False;
    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%lf", &val) == 1))
      valid = xfe_True;
    if (!valid){ xf_printf("Invalid demarcation value.\n"); continue;}

    // flag over bfaces
    nBFaceTot = Mesh->BFaceGroup[ibfgrp].nBFace;
    ierr = xf_Error(xf_Alloc( (void **) &bflag, nBFaceTot, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    nBFaceNew = 0;
    for (ibface=0; ibface< nBFaceTot; ibface++){

      // boundary face
      BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
      egrp = BFace.ElemGroup;
      elem = BFace.Elem;
      face = BFace.Face;
      
      // get Q1 nodes on face
      ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, 
				       Mesh->ElemGroup[egrp].QOrder, 
				       face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;
      
      // calculate rval = average x[d] for bface
      rval = 0.;
      for (i=0; i<nfnode; i++){
	inode = Mesh->ElemGroup[egrp].Node[elem][fvec[i]];
	rval += Mesh->Coord[inode][d];
      }
      rval = rval/((real) nfnode);

      if (rval < val)
	bflag[ibface] = 0;
      else{
	bflag[ibface] = 1;
	nBFaceNew++;
      }

    } // ibface

    if ((nBFaceNew == 0) || (nBFaceNew == nbfgrp)){
      xf_printf("All/None split.  Not doing anything.\n");
      continue;
    }
      

    xf_printf("Enter name of group with x[%d] < %.3E ... (%d faces)\n> ", 
	      d, val, nBFaceTot - nBFaceNew);
    valid = xfe_False;
    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%s", Title1) == 1))
      valid = xfe_True;
    if (!valid){ xf_printf("Invalid name.\n"); continue;}
      
    xf_printf("Enter name of group with x[%d] >= %.3E ... (%d faces)\n> ",
	      d, val, nBFaceNew);
    valid = xfe_False;
    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%s", Title2) == 1))
      valid = xfe_True;
    if (!valid){ xf_printf("Invalid name.\n"); continue;}


    Mesh->nBFaceGroup = nbfgrp = nbfgrp+1; // new number of groups
    xf_printf("New number of boundary groups = %d\n", nbfgrp);
    ierr = xf_Error(xf_ReAlloc((void **) &Mesh->BFaceGroup, nbfgrp, sizeof(xf_BFaceGroup)));
    if (ierr != xf_OK) return ierr;
    ibfgrpNew = nbfgrp-1;
    
    
    // new title for old group
    strcpy(Mesh->BFaceGroup[ibfgrp].Title, Title1);

    /* Create name for new group */
    ierr = xf_Error(xf_AllocString(&Mesh->BFaceGroup[ibfgrpNew].Title, xf_MAXSTRLEN, Title2));
    if (ierr != xf_OK) return ierr;

    // Allocate new BFace
    ierr = xf_Error(xf_Alloc( (void **) &Mesh->BFaceGroup[ibfgrpNew].BFace, 
			      nBFaceNew, sizeof(xf_BFace)));
    if (ierr != xf_OK) return ierr;
    
    // Copy over bfaces that we are splitting off, consolidate ones that stay
    nbface1 = nbface2 = 0;
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      if (bflag[ibface] == 1) 
	Mesh->BFaceGroup[ibfgrpNew].BFace[nbface2++] = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
      else{
	if (ibface != nbface1)
	  Mesh->BFaceGroup[ibfgrp   ].BFace[nbface1] = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
	nbface1++;
      }
    }

    /*     xf_printf("nBFaceTot = %d, nBFaceNew = %d, nbface1 = %d, nbface2 = %d\n", */
    /* 	      nBFaceTot, nBFaceNew, nbface1, nbface2); */

    // sanity check
    if (nbface1 != (nBFaceTot-nBFaceNew)) return xf_Error(xf_CODE_LOGIC_ERROR);
    if (nbface2 !=            nBFaceNew ) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    Mesh->BFaceGroup[ibfgrp   ].nBFace = nbface1;
    Mesh->BFaceGroup[ibfgrpNew].nBFace = nbface2;

    // ReAllocate old BFace
    ierr = xf_Error(xf_ReAlloc( (void **) &Mesh->BFaceGroup[ibfgrp].BFace, 
			      nBFaceTot-nBFaceNew, sizeof(xf_BFace)));
    if (ierr != xf_OK) return ierr;


    xf_Release((void *) bflag);
    
  } // end while ! done

  return xf_OK;
}




/******************************************************************/
//  FUNCTION Definition: xf_ChooseDeleteData
static int 
xf_ChooseDeleteData(xf_All *All)
{
/*

PURPOSE: 

  Interactively deletes Data from All->DataSet
  
INPUTS:

  All : all structure

OUTPUTS:  

RETURNS: Error Code

*/

  int ierr, i, k, count;
  enum xfe_Bool done;
  char buf[xf_MAXSTRLEN];
  char Title[xf_MAXSTRLEN];
  xf_DataSet *DataSet;
  xf_Data *Data;

  DataSet = All->DataSet;

  done = xfe_False;

  while (!done){

    Data = DataSet->Head;
    count = 0;
    while (Data != NULL){
      xf_printf("%4d: %s\n", count++, Data->Title);
      Data = Data->Next;
    }
    xf_printf("Choose data # to delete (or -1 to finish)\n> ");

    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%d", &k) == 1)){
      if ((k >= 0) && (k < count)){
	for (i=0, Data=DataSet->Head; i<k; i++) Data = Data->Next;
	xf_printf("Deleting %s\n", Data->Title);
	ierr = xf_Error(xf_DataSetRemove(DataSet, Data->Title, xfe_False));
	if (ierr != xf_OK) return ierr;
      }
      else if (k == -1)
	done = xfe_True;
      else
	xf_printf("Out of range.\n");
    }
    else
      xf_printf("Not understood.\n");    
  } // end while ! done

  return xf_OK;
}




/******************************************************************/
//  FUNCTION Definition: xf_WriteFile
static int 
xf_WriteFile(xf_All *All, const char *fname)
{
/*

PURPOSE: 

  Writes file name fname, or other
  
INPUTS:

  All : all structure
  fname : default file name (with extension)

OUTPUTS:  

RETURNS: Error Code

*/
  int ierr, len;
  char buf[xf_MAXSTRLEN];
  char outFile[xf_MAXSTRLEN];
  char *outExt;

  xf_printf("Choose file name to write, with extension [%s]:\n> ", fname);
  
  if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
      (strlen(buf)>4)){
    strcpy(outFile,buf);
    if ((len = strlen(outFile)) > 0) outFile[len-1] = '\0';
  }
  else{
    xf_printf("Using default = %s\n", fname);
    strcpy(outFile, fname);
  }


  xf_printf("outFile = %s\n", outFile);
  
  if ((len = strlen(outFile)) < 4) return xf_Error(xf_INPUT_ERROR);
  outExt = outFile + len - 4; // pointer to extension

  if (strncmp(outExt, ".gri", 4) == 0){
    ierr = xf_Error(xf_WriteGriFile(All->Mesh, outFile));
    if (ierr != xf_OK) return ierr;
  }
  else if (strncmp(outExt, ".xfa", 4) == 0){
    ierr = xf_Error(xf_WriteAllBinary(All, outFile));
    if (ierr!=xf_OK) return ierr;
  }
  else
    return xf_Error(xf_NOT_SUPPORTED);

  xf_printf("Wrote %s.\n", outFile);


  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, len;
  char inFile[xf_MAXSTRLEN];
  char buf[xf_MAXSTRLEN];
  char *inExt;
  enum xfe_Bool done;
  xf_All *All;

  xf_printf("\n");
  xf_printf("=== xf_Edit: interactive mesh/all editor  ===\n");
  xf_printf("\n");

  /* Check number of arguments */
  if( argc != 2 ){
    xf_printf("Usage:\n");
    xf_printf("xf_Edit <file>\n");
    xf_printf("\n");
    xf_printf("Where <file> is an .xfa or .gri file.\n");
    xf_printf("\n");
    return xf_Error(xf_INPUT_ERROR);
  }

  // root sets jobFile
  strcpy(inFile, argv[1]);

  // extension is required
  if ((len = strlen(inFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
  inExt = inFile + len - 4; // pointer to extension

  /*-------------*/
  /* Read inFile */
  /*-------------*/

  ierr = xf_Error(xf_ReadAllInputFile(inFile, NULL, xfe_False, &All));
  if (ierr != xf_OK) return ierr;


  /*-----------*/
  /* Main Menu */
  /*-----------*/

  done = xfe_False;
  while (!done){
    xf_printf("\n");
    xf_printf("+---------------+\n");
    xf_printf("|   MAIN MENU   |\n");
    xf_printf("+---------------+\n");
    xf_printf("\n");
    xf_printf("0. Exit\n");
    xf_printf("1. Rename boundary groups\n");
    xf_printf("2. Merge boundary groups\n");
    xf_printf("3. Split boundary groups\n");
    xf_printf("4. Check volumes\n");
    xf_printf("5. Choose and delete data\n");
    xf_printf("6. Write file\n");
    xf_printf("\n");
    xf_printf("Choose option:\n> ");
    
    if ((fgets(buf, xf_MAXSTRLEN, stdin) != NULL) && 
	(sscanf(buf, "%d", &i) == 1) && (i >= 0) && (i < 7)){
      switch(i){
	
      case 0: // Exit
	done = xfe_True;
	break;
	
      case 1:
	ierr = xf_Error(xf_RenameBoundaries(All));
	if (ierr != xf_OK) 
	  xf_printf("Error = %d occured during RenameBoundaries. Continuing.\n", ierr);
	break;

      case 2:
	ierr = xf_Error(xf_MergeBoundaries(All));
	if (ierr != xf_OK) return ierr;
	break;

      case 3:
	ierr = xf_Error(xf_SplitBoundaries(All));
	if (ierr != xf_OK) return ierr;
	break;
	
      case 4:
	ierr = xf_Error(xf_CheckVolumes(All->Mesh, NULL, xfe_False, 0, NULL, NULL));
	if (ierr != xf_OK) return ierr;
	break;

      case 5:
	ierr = xf_Error(xf_ChooseDeleteData(All));
	if (ierr != xf_OK) return ierr;
	break;

      case 6:
	ierr = xf_Error(xf_WriteFile(All, inFile));
	if (ierr != xf_OK) 
	  xf_printf("Error = %d occured during WriteFile. Continuing.\n", ierr);
	break;

      default:
	xf_printf("Out of range.\n");
	break;

      } // end switch i

    }
    else
      xf_printf("Not understood.\n");
  } 
  
  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;
  
  xf_printf("xf_Edit finished.\n");

  return xf_OK;
}
