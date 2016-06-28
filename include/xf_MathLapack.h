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

#ifndef _xf_MathLapack_h
#define _xf_MathLapack_h 1

/*
  FILE:  xf_MathLapack.h

  This file contains the headers for math functions that use LAPACK.

*/

 
/* Note, LAPACK functions work with matrices in column-major format
   (Fortran-style).  Check the below documentation on each function
   regarding calling convention. */


/******************************************************************/
//   FUNCTION Prototype: xf_EigSymTriDiag
extern int
xf_EigSymTriDiag(int n, real *D, real *E, real *Z);
/*
PURPOSE: 

  Computes the eigenvalues and optionally the eigenvectors of a
  symmetric tri-diagonal matrix.  Uses the LAPACK function dsteqr.

INPUTS:

  n : size of matrix
  D : on-diagonal entries of matrix [n]
  E : off-diagonal entries of matrix [n-1] (note, matrix is symmetric)

OUTPUTS:

  D : stores n eigenvalues upon successful calculation
  E : used for storage and hence overwritten
  Z : eigenvectors stores in column-major format; optional output;
      pass NULL to not calculate eigenvectors

RETURN: Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_EigSymTriDiag
extern int
xf_EigSym(int n, real *A, real *E);
/*
PURPOSE: 

  Computes the eigenvalues and eigenvectors of a symmetric matrix, A.
  A is overwritten with the eigenvectors (column-major format).  Uses
  the LAPACK function dsyev.  The eigenvalues are in ascending order.

INPUTS:

  n : size of matrix
  A : symmetrix matrix of interest [nxn]

OUTPUTS:

  E : stores the eigenvalues (must be preallocated) in ascending order
  A : overwritten with eigenvectors

RETURN: Error code

*/

/******************************************************************/
//   FUNCTION  Prototype: Yu_EigenVector
extern int
Yu_EigenVector(int n, real *A, real *VL, real *VR);

/******************************************************************/
//   FUNCTION Prototype: xf_CholDecomp
extern int
xf_CholDecomp(int n, real *A, enum xfe_Bool LUFlag);
/*
PURPOSE: 

  Computes the Cholesky decomposition of a real symmetric positive
  definite matrix A.

  A = U^T * U     (LUFlag = True )
  A = L * L^T     (LUFlag = False)

INPUTS:

  n : size of matrix
  A : symmetrix matrix of interest [nxn]
  LUFLag : True to work with and return the upper portion of A. 
           False to work with and return the lower portion of A

OUTPUTS:

  A : overwritten with Cholesky decomposition

RETURN: Error code

*/



#endif // end ifndef _xf_MathLapack_h
