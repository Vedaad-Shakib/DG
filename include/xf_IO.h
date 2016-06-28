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

#ifndef _xf_IO_h
#define _xf_IO_h 1

/*
  FILE:  xf_IO.h

  This file contains headers for functions in xf_IO.c

*/

#include <stdio.h>


/******************************************************************/
//   FUNCTION Prototype: xf_printf
extern int xf_printf(const char *fmt, ...);
/*
PURPOSE:
  xflow specific print function

INPUTS:
  fmt: format (just as for printf)

OUTPUTS: 
  None; prints the buffer

RETURN:
  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_pprintf
extern int xf_pprintf(char *fmt, ...);
/*
PURPOSE:

  Parallel xflow specific print function.  Prints out [myRank = %d]
  before the string.  Also flushes the stdout buffer. Useful for
  debugging.

INPUTS:
  fmt: format (just as for printf)

OUTPUTS: 
  None; prints the buffer

RETURN:
  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_fread
extern int 
xf_fread(FILE *fid, int size, int n, void *data );
/*
PURPOSE:

  Wrapper for fread.  Safe to call in parallel.

INPUTS:

  fid : file pointer from which to read
  size : size of data to read
  n : number of units of data to read
  
OUTPUTS: 
  
  data : output from read stored in here

RETURN:
  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_fwrite
extern int 
xf_fwrite(const void *data, int size, int n, FILE *fid);
/*
PURPOSE:

  Wrapper for fwrite.  Safe to call in parallel.

INPUTS:

  data : stores data that will be written
  size : size of data
  n : number of units of data to write
  fid : file pointer to write to

  
OUTPUTS: 

RETURN:
  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_fopen
extern int 
xf_fopen(const char *fname, const char *mode, FILE **pfid);
/*
PURPOSE:

  Wrapper for fopen.  Safe to call in parallel.

INPUTS:

  fname : name of file to open
  mode : mode to open in, e.g. "r", "w", "rb", etc.
  
OUTPUTS: 
 
  fid : resulting file pointer

RETURN:
  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_fclose
extern int 
xf_fclose(FILE *fid);
/*
PURPOSE:

  Wrapper for fclose.  Safe to call in parallel.

INPUTS:

  fid : file pointer to close
  
OUTPUTS: 

RETURN:
  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_feof
extern enum xfe_Bool 
xf_feof(FILE *fid);
/*
PURPOSE:

  Wrapper for feof.  Safe to call in parallel.

INPUTS:

  fid : file pointer to check for feof
  
OUTPUTS: 

RETURN:
  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_fseek
extern int 
xf_fseek(FILE *fid, xf_long pos, int whence);
/*
PURPOSE:

  Wrapper for fseek.  Safe to call in parallel.

INPUTS:

  fid : file pointer to seek
  pos : position to seek to
  whence: e.g. SEEK_SET
  
OUTPUTS: 

RETURN:
  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ftell
extern int 
xf_ftell(FILE *fid, xf_long *pos);
/*
PURPOSE:

  Wrapper for ftell.  Safe to call in parallel.

INPUTS:

  fid : file pointer to tell
  
OUTPUTS: 

  pos : position returned by ftell

RETURN:
  Error Code
*/

#endif // end ifndef _xf_IO_h
