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
  FILE:  xf_IO.c

  This file contains input-output routines.

*/



#include <stdio.h>
#include <stdarg.h>

#include "xf_AllStruct.h"
#include "xf_Error.h"
#include "xf_MPI.h"


/******************************************************************/
//   FUNCTION Definition: xf_printf
int 
xf_printf(const char *fmt, ...){
  int myRank;
  va_list argp;
  if (xf_MPI_GetRank(&myRank, NULL) != xf_OK) myRank = 0;

  if (myRank == 0){
    va_start(argp, fmt);
    vprintf(fmt, argp);
    va_end(argp);
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_printf
int 
xf_pprintf(char *fmt, ...){

  int myRank;
  char newfmt[xf_MAXLINELEN];
  va_list argp;

  if (xf_MPI_GetRank(&myRank, NULL) != xf_OK) myRank = -1;
  sprintf(newfmt, "[myRank = %d] %s", myRank, fmt);
  
  va_start(argp, fmt);
  vprintf(newfmt, argp);
  va_end(argp);
  fflush(stdout);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_fread
int 
xf_fread(FILE *fid, int size, int n, void *data )
{
  int ierr, myRank, terr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root performs the read
  if (myRank == 0)
    terr = ((fread(data, size, n, fid) != n) ? xf_FILE_READ_ERROR : xf_OK);

  // read status is broadcasted in case of error
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  // data is broadcasted
  ierr = xf_Error(xf_MPI_Bcast(data, n*size, 0));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_fwrite
int 
xf_fwrite(const void *data, int size, int n, FILE *fid)
{
  int ierr, myRank, terr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  // root performs the write
  if (myRank == 0)
    terr = ((fwrite(data, size, n, fid) != n) ? xf_FILE_WRITE_ERROR : xf_OK);

  // status is broadcasted in case of error
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_fopen
int 
xf_fopen(const char *fname, const char *mode, FILE **pfid)
{
  int ierr, myRank, terr = 0;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  (*pfid) = NULL;

  /* Root opens file */
  if (myRank == 0)
    terr = ((((*pfid) = fopen(fname, mode)) == NULL) ? xf_NOT_FOUND : xf_OK);
  
  // error from file open is broadcasted
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return terr; 

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_fclose
int 
xf_fclose(FILE *fid)
{
  int ierr, myRank, terr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;
  
  /* Root closes file */
  if (myRank == 0)
    terr = ((fclose(fid) != 0) ? xf_FILE_READ_ERROR : xf_OK);
  
  // error from file close is broadcasted
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr); 

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_feof
enum xfe_Bool 
xf_feof(FILE *fid)
{
  int ierr, myRank;
  enum xfe_Bool flag;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  /* Root checks feof */
  if (myRank == 0) flag = feof(fid);

  // flag is broadcasted
  ierr = xf_Error(xf_MPI_Bcast((void *) &flag, sizeof(enum xfe_Bool), 0));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_fseek
int 
xf_fseek(FILE *fid, xf_long pos, int whence)
{
  int ierr, myRank, terr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  /* Root seeks file */
  if (myRank == 0)
    terr = ((fseek(fid, (long) pos, whence) != 0) ? xf_FILE_READ_ERROR : xf_OK);
  
  // error is broadcasted
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr); 

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ftell
int 
xf_ftell(FILE *fid, xf_long *pos)
{
  int ierr, myRank, terr;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
  if (ierr != xf_OK) return ierr;

  /* Root ftells file */
  if (myRank == 0)
    terr = (((*pos = (xf_long) ftell(fid)) < 0) ? xf_FILE_READ_ERROR : xf_OK);

  // error is broadcasted
  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr != xf_OK) return xf_Error(terr); 

  return xf_OK;
}
