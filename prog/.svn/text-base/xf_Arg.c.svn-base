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
  FILE:  xf_Arg.c

  This file contains argument parsing routines.

*/


#include <stdio.h>
#include "xf.h"
#include "xf_ParamStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_Param.h"
#include "xf_String.h"


/******************************************************************/
//   FUNCTION Definition: xf_ParseArg
int 
xf_ParseArg(char **Keys, int argc, char *argv[], xf_KeyValue *KeyValue)
{
  int ierr, outerr = xf_OK, i, nKey;
  int myRank, nProc;
  char s0[xf_MAXSTRLEN];
  char s1[xf_MAXSTRLEN];
  char *key;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // count number of keys
  i = 0;
  while (strlen(Keys[i]) > 0){
    i++;
    if (i > xf_ISAFETY) return xf_Error(xf_INPUT_ERROR);
  }
  if ((i % 3) != 0) return xf_Error(xf_INPUT_ERROR);
  nKey = i/3;

  if (nKey == 0){
    xf_printf("No keys given for parsing.\n");
    return xf_OK;
  }
  
  
  // proc 0 performs the parsing
  if (myRank == 0){
    if (argc < 1) outerr = xf_Error(xf_INPUT_ERROR);
    else if (argc == 1){ // no key-value pairs in arguments
      xf_printf("\nUsage:\n");
      for (i=0; i<nKey; i++){
	sprintf(s0, "-%s", Keys[3*i+0]);
	sprintf(s1, "[%s]", Keys[3*i+1]);
	xf_printf("%12s : %14s : %s\n", s0, s1, Keys[3*i+2]);
      }
      xf_printf("\n");
      outerr = xf_FORCE_QUIT;
    }
    else{
      // pull off Keys, store them in KeyValue
      for (i=0; i<nKey; i++){
	ierr = xf_AddKeyValue(KeyValue, Keys[3*i+0], Keys[3*i+1], xfe_True);
	if (ierr == xf_OVERWROTE){
	  xf_printf("Error, repeated key found in Keys input.\n");
	  outerr = xf_Error(xf_INPUT_ERROR);
	}
	else if (ierr != xf_OK) outerr = xf_Error(ierr);
      }
    }
    
    // proceed with parsing if no errors up to now
    if (outerr == xf_OK){
      for (i=0; i<argc; i++){
	if ((strlen(argv[i]) > 0) && (strncmp(argv[i], "-", 1) == 0) 
	    && (i < (argc-1))){
	  key = argv[i] + 1;
	  ierr = xf_SetKeyValue((*KeyValue), key, argv[i+1]);
	  // do not check for error; allow other keys to be present
	}
      }
      // print out summary
      xf_printf("nKey = %d\n", KeyValue->nKey);
      for (i=0; i<KeyValue->nKey; i++)
	xf_printf("%d : Key = %s, Value = %s\n", i, KeyValue->Key[i], KeyValue->Value[i]);
      xf_printf("\n");

    }
  } // end if myRank == 0

  // broadcast outerr and quit if not OK
  if (xf_PError(&outerr, 0) != xf_OK) return outerr;
 

  // parallelize KeyValue, so every proc gets a copy
  ierr = xf_Error(xf_ParallelizeKeyValue(KeyValue));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
