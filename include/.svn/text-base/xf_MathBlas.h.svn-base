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

#ifndef _xf_MathBlas_h
#define _xf_MathBlas_h 1

/*
  FILE:  xf_MathBlas.h

  This file contains headers for math functions that come in BLAS
  versions.

*/


/******************************************************************/
//   FUNCTION Prototypes: xf_MxM_Set, _Neg, _Add, _Sub
extern void
xf_MxM_Set(const real *A, const real *B, int rA, int n, int cB, real *C);
extern void
xf_MxM_Neg(const real *A, const real *B, int rA, int n, int cB, real *C);
extern void
xf_MxM_Add(const real *A, const real *B, int rA, int n, int cB, real *C);
extern void
xf_MxM_Sub(const real *A, const real *B, int rA, int n, int cB, real *C);
/*
PURPOSE: 

  Computes matrix-matrix product C @= A*B,
  where @= is {=, =-, +=, -=} for {Set, Neg, Add, Sub}

INPUTS:

  A : rA x n  matrix
  B :  n x cB matrix
  rA, n, cB : dimensions introduced above

OUTPUTS:

  C :  rA x cB resultant matrix

RETURN:

  None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MTxM_Set, _Neg, _Add, _Sub
extern void
xf_MTxM_Set(const real *A, const real *B, int cA, int n, int cB, real *C);
extern void
xf_MTxM_Neg(const real *A, const real *B, int cA, int n, int cB, real *C);
extern void
xf_MTxM_Add(const real * RSTRCT A, const real * RSTRCT B, int cA, int n, int cB, 
	    real * RSTRCT C);
extern void
xf_MTxM_Sub(const real * RSTRCT A, const real * RSTRCT B, int cA, int n, int cB, 
	    real * RSTRCT C);
/*
PURPOSE: 

  Computes matrix^T-matrix product C @= A^T*B
  where @= is {=, =-, +=, -=} for {Set, Neg, Add, Sub}

INPUTS:

  A :  n x cA  matrix  (note, A^T is cA x n)
  B :  n x cB matrix
  rA, n, cB : dimensions introduced above

OUTPUTS:

  C :  rA x rB resultant matrix

RETURN:

  None
*/


#endif // end ifndef _xf_MathBlas_h
