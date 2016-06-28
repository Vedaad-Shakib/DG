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

#ifndef _xf_Memory_h
#define _xf_Memory_h 1

/*
  FILE:  xf_Memory.h

  This file contains the headers for functions that handle memory
  management.

*/

/******************************************************************/
//   FUNCTION Prototype:  xf_Alloc
extern int
xf_Alloc( void **pchunk, int n, int size);
/*
PURPOSE:

  Allocates a block of memory. 
  n*size = total bytes allocated

INPUTS:

  n   : Number of elements to be allocated.
  size: Size of each element in bytes.

OUTPUTS:

  pchunk: pointer to memory chunk allocated

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype:  xf_LongAlloc
extern int
xf_LongAlloc( void **pchunk, xf_long n, int size);
/*
PURPOSE:

  Allocates a block of memory, with n of long type.
  n*size = total bytes allocated

INPUTS:

  n   : Number of elements to be allocated.
  size: Size of each element in bytes.

OUTPUTS:

  pchunk: pointer to memory chunk allocated

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype:  xf_Alloc2
extern int
xf_Alloc2( void ***pchunk, int n1, int n2, int size);
/*
PURPOSE:

  Allocs a block of memory for 2d array (a matrix).
  Block of memory of n1*n2*size bytes is allocated.

INPUTS:

  n1  : Number of elements in first index (i.e number of rows).
  n2  : Number of elements in second index (i.e. number of columns).
  size: Size of each element in bytes.

OUTPUTS:

  pchunk: Pointer to memory chunk allocated.

RETURN:

  Error code
*/

/******************************************************************/
//   FUNCTION Prototype:  xf_VAlloc2
extern int
xf_VAlloc2( void ***pchunk, int n1, int *n2, int size);
/*
PURPOSE:

  Allocs a block of memory for variable-length 2d array.  Imagine a
  matrix with n1 rows but with each row i having n2[i] columns,
  0<=i<=n1.  The array is allocated similarly to Alloc2, in that first
  a single chunk of size sum_i(n2[i]) is allocated and then a pointer
  of pointers is allocated into this chunk.  Thus, (*pchunk)[0] points
  to the unwound array.
 
INPUTS:

  n1  : Number of elements in first index (i.e number of rows).
  n2  : Vector number of columns for each row
  size: Size of each element in bytes.

OUTPUTS:

  pchunk: Pointer to memory chunk allocated.

RETURN:

  Error code
*/


/******************************************************************/
//   FUNCTION Prototype:  xf_Alloc3
extern int
xf_Alloc3( void ****pchunk, int n1, int n2, int n3, int size);
/*
PURPOSE:

  Allocs a block of memory for a 3-D array.
  Block of memory of n1*n2*n3*size bytes is allocated.

INPUTS:

  n1  : Number of elements in first index.
  n2  : Number of elements in second index.
  n3  : Number of elements in third index.
  size: Size of each element in bytes.

OUTPUTS:
  
  pchunk: Pointer to memory chunk allocated

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_Release
extern void 
xf_Release( void *chunk);
/*
PURPOSE:

  Frees a block of memory that was allocated using xf_Alloc or xf_ReAlloc.

INPUTS:

  chunk  : Pointer to block of memory to be freed

OUTPUTS:

  None

RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype:  xf_Release2
extern void 
xf_Release2( void **chunk);
/*
PURPOSE:

  Frees a block of memory that was allocated using xf_Alloc2 or xf_ReAlloc2.

INPUTS:
  chunk  : Pointer to pointer to block of memory to be freed.

OUTPUTS:

  None

RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype:  xf_Release3
extern void 
xf_Release3( void ***chunk);
/*
PURPOSE:

  Frees a block of memory that was allocated using xf_Alloc3.

INPUTS:
  chunk  : Triple pointer to block of memory as returned by xf_Alloc3.

OUTPUTS:

  None

RETURN:

  None

*/


/******************************************************************/
//   FUNCTION Prototype:  xf_ReAlloc
extern int
xf_ReAlloc( void **pchunk, int n, int size);
/*
PURPOSE: 

  Reallocates (resizes) a block of memory pointed to by pchunk.
  If pchunk is NULL, xf_Alloc is called.

INPUTS:

  pchunk : Pointer to block of memory.
  n   : Total number of elements requested
  size: Size of each element in bytes.

OUTPUTS:

  pchunk: points to reallocated memory


RETURN:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_ReAlloc2
extern int
xf_ReAlloc2( void ***pchunk, int n1, int n2, int size);
/*
PURPOSE:

  Reallocates a block of memory for a 2d array.  Block of memory of
  n1*n2*size bytes is allocated.  Note that the fcn frees memory that
  pchunk points to before allocating new memory.  Thus, data in pchunk
  is destroyed.


INPUTS:
  pchunk : Pointer to pointer to block of memory.
  n1  : Number of elements in first index.
  n2  : Number of elements in second index.
  size: Size of each element in bytes.

OUTPUTS: 
 
  pchunk : pointer to reallocated memory

RETURN:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_VReAlloc2
extern int
xf_VReAlloc2( void ***pchunk, int n1, int *n2, int size);
/*
PURPOSE:

  Reallocates a block of memory for a variable-length 2d array.  Block
  of memory of n1*sum(n2)*size bytes is allocated.  Note that the fcn
  frees memory that pchunk points to before allocating new memory.
  Thus, data in pchunk is destroyed.


INPUTS:
  pchunk : Pointer to pointer to block of memory.
  n1  : Number of elements in first index (i.e number of rows).
  n2  : Vector number of columns for each row
  size: Size of each element in bytes.

OUTPUTS: 
 
  pchunk : pointer to reallocated memory

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype:  xf_VReAllocCopy2
extern int
xf_VReAllocCopy2( void ***pchunk, int p1, int *p2, int n1, int *n2, int size);
/*
PURPOSE:

  Reallocates a block of memory for a 2-D variable array and copies
  data.  Block of memory of n1*sum(n2)*size bytes is allocated and old
  data is copied to new memory block where possible.

INPUTS:

  pchunk : (Pointer to) pointer to original block of memory.
  p1  : Original first index size (scalar)
  p2  : Original second index size vector
  n1  : New first index size (scalar)
  n2  : New second index size vector
  size: Size of each element in bytes.

OUTPUTS: 

  pchunk : pointer to reallocated memory

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype:  xf_ReAllocCopy2
extern int
xf_ReAllocCopy2( void ***pchunk, int p1, int p2, int n1, int n2, int size);
/*
PURPOSE:

  Reallocates a block of memory for a 2-D array and copies data.
  Block of memory of n1*n2*size bytes is allocated and old data is
  copied to new memory block where possible.

INPUTS:

  pchunk : (Pointer to) double pointer to original block of memory.
  p1  : Original first index size.
  p2  : Original second index size.
  n1  : New first index size.
  n2  : New second index size.
  size: Size of each element in bytes.

OUTPUTS: 

  pchunk : pointer to reallocated memory

RETURN:

  Error code

*/



#endif // end ifndef _xf_Memory_h
