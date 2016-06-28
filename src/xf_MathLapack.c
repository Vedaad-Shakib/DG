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
  FILE:  xf_MathLapack.c

  This file contains top-level wrappers for math functions that are
  implemented in LAPACK.

*/

#include "xf.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"


/* Lapack prototypes */
void *dsteqr_(char *, int *, real *, real *, real *, int *, real *, int *);
void *dsyev_(char *, char *, int *, real *, int *, real *, real *, int *, int *);
void *dpotrf_(char *, int *, real *, int *, int *);
void *dgeev_(char *, char*, int *, real *, int *, real *, real *, real *, int *,
             real *, int *, real *, int *, int *);

/******************************************************************/
//   FUNCTION Definition: xf_EigSymTriDiag
int
xf_EigSymTriDiag(int n, real *D, real *E, real *Z)
{
  int ierr, ldz = 1;
  int k, info;
  char compz[] = "N";
  real *work;

  // allocate work
  ierr = xf_Error(xf_Alloc((void **) &work, 2*n, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // prepare to calculate eigenvectors
  if (Z != NULL){
    strcpy(compz, "I");
    ldz = n;
    // set Z = identity
    for (k=0; k<n*n; k++) Z[k] = 0.0;
    for (k=0; k<n*n; k+=(n+1)) Z[k] = 1.0;
  }

  /* DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )*/
  dsteqr_(compz, &n, D, E, Z, &ldz, work, &info);

  if (info != 0) return xf_Error(xf_LAPACK_ERROR);

  xf_Release( (void *) work);

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_EigSym
int
xf_EigSym(int n, real *A, real *E)
{
  int ierr, lwork;
  int info;
  char jobz[] = "V"; // always calculate eigenvectors
  char uplo[] = "U";
  real *work;

  // allocate work
  lwork = max(n*n,3*n); lwork = max(lwork, 1);
  ierr = xf_Error(xf_Alloc((void **) &work, lwork, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  /* DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO ) */
  dsyev_(jobz, uplo, &n, A, &n, E, work, &lwork, &info);

  if (info != 0) return xf_Error(xf_LAPACK_ERROR);

  xf_Release( (void *) work);

  return xf_OK;
}

/******************************************************************/
//      FUNCTION Definition: Yu_EigenVector
int 
Yu_EigenVector(int n, real *A, real *VL, real *VR)
{
   int i, j, ierr, lwork;
   int info;
   char jobvl[] = "V";
   char jobvr[] = "V";
   real *work, *wRe, *wIm;

   //allocate work
   lwork = max(1, 4*n);
   ierr = xf_Error(xf_Alloc((void **) &work, lwork, sizeof(real)));
   if (ierr != xf_OK) return ierr;

   ierr = xf_Error(xf_Alloc((void **) &wRe, n, sizeof(real)));
   if (ierr != xf_OK) return ierr;

   ierr = xf_Error(xf_Alloc((void **) &wIm, n, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   
   //SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
   dgeev_(jobvl, jobvr, &n, A, &n, wRe, wIm, VL, &n, VR, &n, work, &lwork, &info);

   for(i=0; i<4; i++)
   {
      for(j=0; j<4; j++)
         printf("%.12lf ", A[4*i + j]);
   printf("\n");
   }
   printf("\n");

   printf("%.15lf %.15lf %.15lf %.15lf", wRe[0], wRe[1], wRe[2], wRe[3]);
   printf("\n");
   if (info != 0) return xf_Error(xf_LAPACK_ERROR);

   xf_Release( (void *) work);
   xf_Release( (void *) wRe);
   xf_Release( (void *) wIm);

   return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CholDecomp
int
xf_CholDecomp(int n, real *A, enum xfe_Bool LUFlag)
{
  int info;
  char uploL[] = "L";
  char uploU[] = "U";
  char *uplo;

  uplo = ((LUFlag) ? uploL : uploU);

  // SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
  dpotrf_(uplo, &n, A, &n, &info);

  if (info != 0) return xf_Error(xf_LAPACK_ERROR);

  return xf_OK;
}


#if( UNIT_TEST==1 )
#include "xf_MathLapack.test.in"
#endif
