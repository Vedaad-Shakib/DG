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
  file   xf_Error.c

  This file contains the error handling functions.

*/

#include "xf.h"
#include "xf_IO.h"
#include "xf_MPI.h"


/******************************************************************/
//   FUNCTION Definition: xf_ErrorReport
int xf_ErrorReport( char *file, int line, char *call, int ierr){

  int myRank;

  if (ierr == xf_OK) return xf_OK;

  if (xf_MPI_GetRank(&myRank, NULL) != xf_OK) myRank = -1;

  xf_pprintf("Error %d has occured [myRank = %d].\n", ierr, myRank);
  xf_pprintf("   File : %s   Line : %d\n   Call : %s\n", file, line, call);
  fflush(stdout);
  return ierr;
}

/******************************************************************/
//   FUNCTION Definition: xf_PErrorReport
int xf_PErrorReport( char *file, int line, char *call, int *pierr, int root){
  int terr;

  // broadcast ierr from root
  terr = xf_Error(xf_MPI_Bcast( (void *) pierr, sizeof(int), root));
  if (terr != xf_OK) return terr;

  if ((*pierr) == xf_OK) return xf_OK;

  //xf_printf("Error %d has occured.\n", ierr);
  //xf_printf("   File : %s   Line : %d\n   Call : %s\n", file, line, call);
  //fflush(stdout);
  return (*pierr);
}


