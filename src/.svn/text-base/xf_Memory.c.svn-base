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
  FILE:  xf_Memory.c

  This file contains functions for memory handling

*/

#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "xf.h"
#include "xf_IO.h"

/******************************************************************/
//   FUNCTION Definition: xf_Alloc
int
xf_Alloc( void **pchunk, int n, int size)
{
  xf_long totalsize;
  
  (*pchunk) = NULL;
  
  totalsize = (xf_long) n*size;

  if (totalsize == 0) return xf_OK;
  
  if (totalsize < 0){
    xf_printf("Warning, requesting allocation of negative memory size.\n");
    return xf_OK;
  }
  
  if (((*pchunk) = (void *)malloc(totalsize)) == NULL)
    return xf_Error(xf_MEMORY_ERROR);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_LongAlloc
int
xf_LongAlloc( void **pchunk, xf_long n, int size)
{
  xf_long totalsize;
  
  (*pchunk) = NULL;
  
  totalsize = (xf_long) n*size;

  if (totalsize == 0) return xf_OK;
  
  if (totalsize < 0){
    xf_printf("Warning, requesting allocation of negative memory size.\n");
    return xf_OK;
  }
  
  if (((*pchunk) = (void *)malloc(totalsize)) == NULL)
    return xf_Error(xf_MEMORY_ERROR);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition:  xf_Alloc2
int
xf_Alloc2( void ***pchunk, int n1, int n2, int size)
{
  char *temp;
  int i;
  xf_long totalsize;
  
  (*pchunk) = NULL;
  
  totalsize = (xf_long) n1*n2*size;
  
  if (totalsize == 0) return xf_OK;

  if (totalsize < 0){
    xf_printf("Warning, requesting allocation of negative memory size.\n");
    return xf_OK;
  }

  if ((temp = (char *)malloc( totalsize)) == NULL)
    return xf_Error(xf_MEMORY_ERROR);

  if (((*pchunk) = (void **)malloc( n1*sizeof(char *) )) == NULL)
    return xf_Error(xf_MEMORY_ERROR);
  
  for(i = 0; i<n1; i++)
    (*pchunk)[i] = temp + (xf_long) i*n2*size;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition:  xf_VAlloc2
int
xf_VAlloc2( void ***pchunk, int n1, int *n2, int size)
{
  char *temp;
  int i;
  xf_long totalsize;
  
  (*pchunk) = NULL;
  
  for (i=0, totalsize=0; i<n1; i++) totalsize += (xf_long) n2[i]*size;
  
  if (totalsize == 0) return xf_OK;

  if (totalsize < 0){
    xf_printf("Warning, requesting allocation of negative memory size.\n");
    return xf_OK;
  }

  if ((temp = (char *)malloc( totalsize)) == NULL)
    return xf_Error(xf_MEMORY_ERROR);

  if (((*pchunk) = (void **)malloc( n1*sizeof(char *) )) == NULL)
    return xf_Error(xf_MEMORY_ERROR);
  
  for(i=0, totalsize=0; i<n1; i++){
    (*pchunk)[i] = temp + totalsize; 
    totalsize += (xf_long) n2[i]*size;
  }
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition:  xf_Alloc3
int
xf_Alloc3( void ****pchunk, int n1, int n2, int n3, int size)
{
  char *temp,  **temp2;
  int i, j;
  xf_long totalsize;

  (*pchunk) = NULL;
  
  totalsize = (xf_long) n1*n2*n3*size;
  
  if (totalsize <= 0) return xf_OK;
  
  if ((temp = (char *)malloc( totalsize)) == NULL)
    return xf_Error(xf_MEMORY_ERROR);
  
  if ((temp2 = (char **)malloc( n1*n2*sizeof(char *) )) == NULL)
    return xf_Error(xf_MEMORY_ERROR);
  
  if (((*pchunk) = (void ***)malloc( n1* sizeof(char **) )) == NULL)
    return xf_Error(xf_MEMORY_ERROR);
  
  
  for(i = 0; i<n1; i++){            
    for(j = 0; j<n2; j++)
      temp2[i*n2+j] = temp + (xf_long) (i*n2+j)*n3*size;
    
    (*pchunk)[i] = (void **) (temp2 + i*n2);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition:  xf_Release
void
xf_Release( void *chunk)
{
  if (chunk == NULL) return;
  free( (void *)chunk);
}

/******************************************************************/
//   FUNCTION Definition:  xf_Release2
void
xf_Release2(void **chunk)
{
  if (chunk == NULL) return;
  free( (void * ) chunk[0]);
  free( (void **) chunk   );
}

/******************************************************************/
//   FUNCTION Definition:  xf_Release3
void
xf_Release3(void ***chunk)
{
  if (chunk == NULL) return;
  free( (void * )  chunk[0][0]);
  free( (void **)  chunk[0]   );
  free( (void ***) chunk   );
}

/******************************************************************/
//   FUNCTION Definition:  xf_ReAlloc
int
xf_ReAlloc( void **pchunk, const int n, const int size)
{
  int merr;
  xf_long totalsize;
  
  if ((*pchunk) == NULL)
    return xf_Error(xf_Alloc(pchunk, n, size));
  
  totalsize = (xf_long) n*size;
  
  if (totalsize <= 0) {
    xf_Release((*pchunk));
    (*pchunk) = NULL;
    return xf_OK;
  }
  
  if (((*pchunk) = (char *)realloc( (void *)(*pchunk), totalsize)) == NULL) {
    merr = errno;
    xf_printf("In xf_ReAlloc: totalsize = %d, errno = %d\n", totalsize, merr);
    if (merr == ENOMEM) xf_printf("ENOMEM\n");
    if (merr == EAGAIN) xf_printf("EAGAIN\n");
    return xf_Error(xf_MEMORY_ERROR);
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition:  xf_ReAlloc2
int
xf_ReAlloc2( void ***pchunk, int n1, int n2, int size)
{
  xf_long totalsize;

  if ((*pchunk) == NULL)
    return xf_Error(xf_Alloc2(pchunk, n1, n2, size));
  
  totalsize = (xf_long) n1*n2*size;

  if (totalsize < 0)
    return xf_Error(xf_MEMORY_ERROR);
  
  xf_Release2((*pchunk));
  return xf_Error(xf_Alloc2( pchunk, n1, n2, size));

}

/******************************************************************/
//   FUNCTION Definition:  xf_VReAlloc2
int
xf_VReAlloc2( void ***pchunk, int n1, int *n2, int size)
{
  xf_long totalsize;

  xf_Release2((*pchunk));
  return xf_Error(xf_VAlloc2( pchunk, n1, n2, size));

}



/******************************************************************/
//   FUNCTION Definition:  xf_VReAllocCopy2
int
xf_VReAllocCopy2( void ***pchunk, int p1, int *p2, int n1, int *n2, int size)
{
   char **po;
   int ierr, i, m1, m2;
   xf_long totalsize;
   
   if ((*pchunk) == NULL)
     return xf_Error(xf_VAlloc2( pchunk, n1, n2, size));

   totalsize = 0;
   for (i=0; i<n1; i++) totalsize += (xf_long) n2[i]*size;

   if (totalsize < 0) return xf_Error(xf_MEMORY_ERROR);

   ierr = xf_Error(xf_VAlloc2( (void ***) &po, n1, n2, size));
   if (ierr != xf_OK) return ierr;

   if(p1<n1)
     m1 = p1;
   else
     m1 = n1;
   
   for (i=0; i<m1; i++){
     if(p2[i]<n2[i])
       m2 = p2[i];
     else
       m2 = n2[i];

     if (m2 <= 0) return xf_Error(xf_MEMORY_ERROR);
  
     po[i] = (char *)memcpy((void *)po[i], (*pchunk)[i], m2*size);
   }

   xf_Release2((*pchunk));

   (*pchunk) = (void **) po;

   return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition:  xf_ReAllocCopy2
int
xf_ReAllocCopy2( void ***pchunk, int p1, int p2, int n1, int n2, int size)
{
   char **po;
   int ierr, i, m1, m2;
   xf_long totalsize;
   
   if ((*pchunk) == NULL)
     return xf_Error(xf_Alloc2(pchunk, n1, n2, size));

   totalsize = (xf_long) n1*n2*size;

   if (totalsize < 0)
     return xf_Error(xf_MEMORY_ERROR);

   ierr = xf_Error(xf_Alloc2( (void ***) &po, n1, n2, size));
   if (ierr != xf_OK) return ierr;

   if(p1<n1)
     m1 = p1;
   else
     m1 = n1;

   if(p2<n2)
     m2 = p2;
   else
     m2 = n2;
  
   if ((m1<=0)||(m2<=0))
     return xf_Error(xf_MEMORY_ERROR);
   
   for (i=0; i<m1; i++)
     po[i] = (char *)memcpy((void *)po[i], (*pchunk)[i], m2*size);

   xf_Release2((*pchunk));

   (*pchunk) = (void **) po;

   return xf_OK;
}
