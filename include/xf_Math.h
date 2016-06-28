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

#ifndef _xf_Math_h
#define _xf_Math_h 1

/*
  FILE:  xf_Math.h

  This file contains the headers for math functions.

*/

#include <math.h>
#include <time.h>
#include "xf.h"
#include "xf_MathBlas.h"

/* Note, all matrices are assumed stored unrolled in row-major form,
   unless specified otherwise. */


/******************************************************************/
//   FUNCTION Prototype: xf_GetAddFlag2
extern enum xfe_AddType
xf_GetAddFlag2(enum xfe_AddType AddFlag);  
/*
  PURPOSE:

  Returns add type used for stringing multiple AddFlag operations
  together.  For example, if the operation V @= (A + B), @ = AddFlag,
  was broken up into two to avoid a temp vector, the two operations
  would be:  V @= A followed by  V &= B, where & = AddFlag2.

  INPUTS:  AddFlag : Input operation

  OUTPUTS:  None

  RETURN: AddFlag2 for stringing subsequent AddFlag operations
*/


/******************************************************************/
//   FUNCTION Prototype: xf_GetAddFlagNeg
extern enum xfe_AddType
xf_GetAddFlagNeg(enum xfe_AddType AddFlag);  
/*
  PURPOSE:

  Returns add type used for taking the negative operation of Add Flag:
  
  AddFlag -> returned NegAddFlag
  Add     ->          Sub
  Sub     ->          Add
  Set     ->          Neg
  Neg     ->          Set

  INPUTS:  AddFlag : Input operation

  OUTPUTS:  None

  RETURN: AddFlagNeg for opposite operation
*/


/******************************************************************/
//   FUNCTION Prototype: xf_PowInt
extern real
xf_PowInt(real x, int n);
/*
  PURPOSE:

  Calculates x^n via a simple loop.  If n < 0, returns 0 (useful for
  differentiation).

  INPUTS:  x : base
  n : power

  OUTPUTS:  None

  RETURN: x^n
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Distance
extern real
xf_Distance(real *x0, real *x1, int dim);
/*
PURPOSE: 

  Calculates distance between x0 and x1

INPUTS: 
 
  x0,x1 : points between which distance is desired
 
OUTPUTS:
 
RETURNS: distance

*/


/******************************************************************/
//   FUNCTION Prototype: xf_V_Add
extern void 
xf_V_Add(const real *u, int n, enum xfe_AddType AddFlag, real *v);
/*
  PURPOSE:

  Sets v @= u, for vectors u,v, where @ is the AddFlag (+/-/etc.)

  INPUTS:  

  u : a vector of size n
  n : size of vectors
  AddFlag : Operation to perform (add, sub, set, neg)

  OUTPUTS:
 
  v : modified output vector, size n

  RETURN: None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_cV_Add
extern void 
xf_cV_Add(const real *u, const real c, int n, enum xfe_AddType AddFlag, real *v);
/*
  PURPOSE:

  Sets v @= c*u, for vectors u,v, where @ is the AddFlag (+/-/etc.)

  INPUTS:  

  u : a vector of size n
  c : scaling constant
  n : size of vectors
  AddFlag : Operation to perform (add, sub, set, neg)

  OUTPUTS:
 
  v : modified output vector, size n

  RETURN: None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_MxV_Add, _Sub, _Set, _Neg
extern void
xf_MxV_Add(const real *A, const real *u, int rA, int cA, real *v);
extern void
xf_MxV_Sub(const real *A, const real *u, int rA, int cA, real *v);
extern void
xf_MxV_Set(const real *A, const real *u, int rA, int cA, real *v);
extern void
xf_MxV_Neg(const real *A, const real *u, int rA, int cA, real *v);
extern void
xf_MxV(const real *A, const real *u, int rA, int cA, 
       enum xfe_AddType AddFlag, real *v);
/*
  PURPOSE: 

  Computes matrix-vector product v @= A*u
  where @= is +=, -=, or = (AddFlag)

  INPUTS:

  A : rA x cA  matrix
  u : cA x 1 vector
  rA, cA : dimensions of A

  OUTPUTS:

  v : rA x 1 output vector

  RETURN:

  None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_cMxV_Add
extern void
xf_cMxV_Add(real c, real *A, real *u, int rA, int cA, real *v);
/*
  PURPOSE: 

  Computes matrix-vector product v += c*A*u, where c is a const

  INPUTS:

  c : scalar constant
  A : rA x cA  matrix
  u : cA x 1 vector
  rA, cA : dimensions of A

  OUTPUTS:

  v : rA x 1 output vector

  RETURN:

  None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MTxV
extern void
xf_MTxV(const real *A, const real *u, int cA, int rA, 
	enum xfe_AddType AddFlag, real *v);
/*
  PURPOSE: 

  Computes matrix-vector product v @= A^T*u
  where @= is +=, -=, or = (AddFlag)

  INPUTS:

  A : rA x cA  matrix
  u : rA x 1 vector
  rA, cA : dimensions of A

  OUTPUTS:

  v : cA x 1 output vector

  RETURN:

  None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ndMTxVc
extern void
xf_ndMTxVc(int n, const real *A, const real *f, int d, real *u);
/*
  PURPOSE: 

  Computes  u{i,:} = A{i,:,:}*f{i}*u{i,:},  where (:) = [0,..,d-1]
                                            and i = 0 .. n-1
  u and A are unrolled about d, d*d, respectively.

  INPUTS:
 
  n : number of products
  A : n unrolled dxd matrices
  f : n scalars
  d : size, <=3
  u : n unrolled dx1 vectors
  
  OUTPUTS:

  u : modified according to the above description

  RETURN: None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ndMTxVic
extern void
xf_ndMTxVic(int n, const real *A, const real *f, int d, real *u);
/*
  PURPOSE: 

  Computes  u{i,:} = A{i,:,:}*(1/f{i})*u{i,:},  where (:) = [0,..,d-1]
                                            and i = 0 .. n-1
  u and A are unrolled about d, d*d, respectively.

  INPUTS:
 
  n : number of products
  A : n unrolled dxd matrices
  f : n scalars
  d : size, <=3
  u : n unrolled dx1 vectors
  
  OUTPUTS:

  u : modified according to the above description

  RETURN: None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_nMxM_Set, _Add
extern void
xf_nMxM_Set(int n, const real *A, const real *B, int rA, int cA, 
	    int cB, real *C);
extern void
xf_nMxM_Add(int n, const real *A, const real *B, int rA, int cA, 
	    int cB, real *C);
/*
  PURPOSE: 

  Computes matrix-matrix products C[n] @= A[n]*B[n], where A[n],
  etc. are sequentially stored matrices of sizes detailed below.

  INPUTS:

  n : number of products to compute
  A : n * [rA x cA] matrices
  B : n * [cA x cB] matrices
  rA, cA, cB : dimensions introduced above

  OUTPUTS:

  C :  n * [rA x cB] resultant matrices

  RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_cMxM_Add
extern void
xf_cMxM_Add(real c, real *A, real *B, int rA, int n, int cB, real *C);
/*
  PURPOSE: 

  Computes matrix-matrix product C += c*A*B

  INPUTS:

  c : scalar constant
  A : rA x n  matrix
  B :  n x cB matrix
  rA, n, cB : dimensions introduced above

  OUTPUTS:

  C :  rA x cB resultant matrix

  RETURN:

  None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_MTxM_Set, _Add, _Sub
extern void
xf_MTxM(const real *A, const real *B, int cA, int n, int cB, 
	enum xfe_AddType AddFlag, real *C);
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


/******************************************************************/
//   FUNCTION Prototype: xf_MTxwM_Set, _Add, _Sub
extern void
xf_MTxwM_Set(const real *A, const real *w, const real *B, int cA, int n, int cB, real *C);
extern void
xf_MTxwM_Add(const real *A, const real *w, const real *B, int cA, int n, int cB, real *C);
/*
  PURPOSE: 

  Computes weighted matrix^T-matrix product C @= A^T*diag(w)*B
  where diag(w) is a nxn matrix with w on the main diagonal

  INPUTS:

  A :  n x cA matrix  (note, A^T is cA x n)
  w :  n x 1  weighting vector
  B :  n x cB matrix
  rA, n, cB : dimensions introduced above

  OUTPUTS:

  C :  rA x rB resultant matrix

  RETURN:

  None
*/




/******************************************************************/
//   FUNCTION Prototype: xf_MxMT_Set
extern void
xf_MxMT_Set(const real *A, const real *B, int rA, int n, int rB, real *C);
/*
  PURPOSE: 

  Computes matrix-matrix^T product C @= A*B^T
  where @= is {=, =-, +=, -=} for {Set, Neg, Add, Sub}

  INPUTS:

  A :  rA x n  matrix
  B :  rB x n  matrix  (note, B^T is n x rB)
  rA, n, rB : dimensions introduced above

  OUTPUTS:

  C :  rA x rB resultant matrix

  RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ColMult, _Set, _Add, _Sub
extern void
xf_ColMult(real *A, const real *v, int rA, int cA, int dv);
extern void
xf_ColDiv(real *A, const real *v, int rA, int cA, int dv);
extern void
xf_ColMult_Set(const real *A, const real *v, int rA, int cA, 
	       int dv, real *B);
extern void
xf_ColMult_Add(const real *A, const real *v, int rA, int cA, 
	       int dv, real *B);
extern void
xf_ColMult_Sub(const real *A, const real *v, int rA, int cA, 
	       int dv, real *B);
extern void
xf_ColcMult(real *A, const real *v, int rA, int cA, int dv, real f);
extern void
xf_ColcMult_Add(const real *A, const real *v, int rA, int cA, 
		int dv, real f, real *B);
extern void
xf_ColcMult_Set(const real *A, const real *v, int rA, int cA, 
		int dv, real f, real *B);
/*
  PURPOSE: 

  Multiplies columns of A entry-wise by column vector v.  In Matlab
  notation,

  B @= A.*repmat(v,1,cA)
  
  where @= is {=, =-, +=, -=} for {Set, Neg, Add, Sub}
  ColMult without a B input just sets A = A.*repmat(v,1,cA)

  To allow for the case when v is part of a matrix, an offset dv is
  specified, so that the column vector equals:

  column vector = [v, v+dv, v+2*dv, ..., v+rA*dv]^T

  ColcMult also multiplies A by a constant, f

  ColDiv divides instead of multiplying

  INPUTS:

  A :  rA x cA  matrix 
  v :  column vector with rA entries offset by dv in memory
  rA, cA : dimensions of A
  dv : offset for v storage (see above)
  f  : additional constant that multiplies A

  OUTPUTS:

  B :  stores resulting matrix

  RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_2ColMult_Set, xf_2ColcMult_Set
extern void
xf_2ColMult_Set(const real *A, const real *v, const real *w, int rA, 
		int cA, int dv, int dw, real *B);
extern void
xf_2ColMult_Add(const real *A, const real *v, const real *w, int rA, 
		int cA, int dv, int dw, real *B);
extern void
xf_2ColcMult_Set(const real *A, const real *v, const real *w, int rA, 
		 int cA, int dv, int dw, real f, real *B);
extern void
xf_2ColcMult_Add(const real *A, const real *v, const real *w, int rA, 
		 int cA, int dv, int dw, real f, real *B);
/*
  PURPOSE: 

  Multiplies columns of A entry-wise by column vectors v and w, and
  stores the result in B.  In Matlab notation,

  B = A.*repmat(v,1,cA).*repmat(w,1,cA)

  To allow for the case when v(w) is part of a matrix, an offset dv
  (dw) is specified, so that the column vector equals:

  column vector = [v, v+dv, v+2*dv, ..., v+rA*dv]^T  (similarly for w)

  2ColcMult_Set, etc. also multiply A by a constant, f.

  INPUTS:

  A :  rA x cA  matrix 
  v :  column vector with rA entries offset by dv in memory
  w :  column vector with rA entries offset by dw in memory
  rA, cA : dimensions of A
  dv : offset for v storage (see above)
  dw : offset for w storage (see above)
  f  : additional constant that multiplies A

  OUTPUTS:

  B :  stores resulting product

  RETURN:

  None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_3ColcMult_Add
extern void
xf_3ColcMult_Add(const real *A, const real *v, const real *w, 
		 const real *u, int rA, int cA, int dv, int dw, 
		 int du, real f, real *B);
/*
  PURPOSE: 

  Multiplies columns of A entry-wise by column vectors v, w, and u, and
  stores the result in B.  In Matlab notation,

  B += f*A.*repmat(v,1,cA).*repmat(w,1,cA).*repmat(u,1,cA)

  To allow for the case when v(w,u) is part of a matrix, an offset dv
  (dw,du) is specified, so that the column vector equals:

  column vector = [v, v+dv, v+2*dv, ..., v+rA*dv]^T  (similarly for w,u)


  INPUTS:

  A :  rA x cA  matrix 
  v :  column vector with rA entries offset by dv in memory
  w :  column vector with rA entries offset by dw in memory
  u :  column vector with rA entries offset by du in memory
  rA, cA : dimensions of A
  dv : offset for v storage (see above)
  dw : offset for w storage (see above)
  du : offset for u storage (see above)
  f  : additional constant that multiplies A

  OUTPUTS:

  B :  stores resulting product (incremented)

  RETURN:

  None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_MatDetInv
extern int
xf_MatDetInv(real *A, int d, real *detA, real *iA);
/*
  PURPOSE: 

  Calculates determinant and/or inverse of A for 1x1, 2x2, or 3x3 matrices.

  INPUTS:

  A :  d x d  matrix (d = 1, 2, or 3)
  d :  dim of A
  
  OUTPUTS:

  detA :  if not NULL, (*detA) will store the determinant of A
  iA :  if not NULL, (*iA) will store the inverse of A

  RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_dMxM
extern void
xf_dMxM(const real *D, real f, int m, int n, int offset, real *A);
/*
  PURPOSE: 

  Computes matrix-matrix product between D and a subset of A
  (e.g. when A has a non-standard structure).

  INPUTS:

  D : m x m matrix
  f : constant factor to use in the multiplication
  m : size of square matrix D (<=3)
  n : number of columns in A
  offset : index offset for rows of A
  A : mxn matrix where rows are separated by offset
  
  OUTPUTS:

  A : gets set to D*A in-place

  RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_dMxMT
extern void
xf_dMxMT(const real *D, real f, int m, int n, int offset, real *A);
/*
  PURPOSE: 

  Computes matrix-matrix product between A^T and D, when A has a
  non-standard structure.

  INPUTS:

  D : m x m matrix
  f : constant factor to use in the multiplication
  m : size of square matrix D (<= 3)
  n : number of columns in A
  offset : index offset for rows of A
  A : mxn matrix where rows are separated by offset
  
  OUTPUTS:

  A : gets set to (A^T*D)^T in-place

  RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DotProduct
extern void
xf_DotProduct(const real *a, const real *b, int n, real *dp);
/*
  PURPOSE: 

  Calculates the dot product dp = a.b

  INPUTS:

  a : [n] vector
  b : [n] vector
  
  OUTPUTS:

  dp : scalar dot product

  RETURN:

  None
*/

/******************************************************************/
//   FUNCTION Prototype: xf_CrossProduct
extern void
xf_CrossProduct(real *a, real *b, real *c);
/*
  PURPOSE: 

  Calculates the cross product c = a x b

  INPUTS:

  a : [3] vector
  b : [3] vector
  
  OUTPUTS:

  c : [3] vector = a cross b

  RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ComputePLU
extern int
xf_ComputePLU(real * RSTRCT A, int r, int * RSTRCT P);
/*
  PURPOSE: 

  Computes PLU of A and stores it in A.  A is a r x r matrix,
  The permutation vector is stored in P.

  INPUTS:

  A : r x r matrix
  r : size of A
  
  OUTPUTS:

  A : PLU factorization is stored here
  P : stores permutation 

  RETURN:

  xf_SINGULAR : if matrix is singular to machine precision

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ComputeBlockPLU
extern int
xf_ComputeBlockPLU(real  *A, int n, int s, int *P);
/*
  PURPOSE: 

  Computes PLU of A and stores it in A.  A is a n x n block matrix,
  where each block is s x s.  The permutation vector is stored in P.

  INPUTS:

  A : n x n block matrix, each block is s x s
  The n x n blocks are stored unrolled:
  A = [
  B00_0, ..., B00_s2; ... B0n_0, ..., B0n_s2;
  ..
  Bn0_0, ..., Bn0_s2; ... Bnn_0, ..., Bnn_s2;
  ]
  n : A is an n x n block matrix
  s : each block is s x s
  
  OUTPUTS:

  A : PLU factorization is stored here
  P : stores permutation 

  RETURN:

  xf_SINGULAR : if matrix is singular to machine precision

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SolveU
extern int
xf_SolveU(real *U, int r, real *u);
/*
  PURPOSE: 

  Sets u = U^{-1}*u, where U is an upper triangular matrix

  INPUTS:

  U  : upper triangular r x r matrix
  r  : size of square matrix U
  u  : on input, right hand side of equation (r x 1)

  OUTPUTS:

  u  : on output, solution to equation

  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SolvePLU
extern int
xf_SolvePLU(real * RSTRCT A, int * RSTRCT P, real * RSTRCT b, int r, 
	    real * RSTRCT u, real * RSTRCT yin);
/*
  PURPOSE: 

  Sets u = A^{-1}*b when A has been PLU factored. 

  INPUTS:

  A  : PLU factored r x r matrix,
  P  : r x 1 permutation vector (if NULL, P=identity is assumed)
  r  : A is an r x r matrix
  yin : optional, memory for temporary array to prevent many alloc/
  release calls if called repeatedly

  OUTPUTS:

  u : an r x 1 solution vector


  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SolvePLU_Matrix
extern int
xf_SolvePLU_Matrix(real *A, int *P, int r, int cB, real *B);
/*
  PURPOSE: 

  Sets  B = A^{-1}*B when A has been PLU factored. 

  INPUTS:

  A  : PLU factored r x r matrix,
  r  : A is an r x r block matrix
  cB : # columns in B

  OUTPUTS:

  B : an r x cB matrix that is set to A^{-1}*B


  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SolveBlockPLU
extern int
xf_SolveBlockPLU(real *A, int n, int s, int *P, real *b, 
		 enum xfe_AddType AddFlag, real *u);
/*
  PURPOSE: 

  Solves:   A*u = b, when A has been Block-PLU factored. 

  INPUTS:

  A : PLU n x n block matrix, each block is s x s
  n : A is an n x n block matrix
  s : each block is s x s
  P : stores permutation vector
  b : rhs
  AddFlag : whether to set, add-to, or subtract-from u


  OUTPUTS:

  u : is modified according to AddFlag with A\b


  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SolveBlockPLUT
extern int
xf_SolveBlockPLUT(real *A, int n, int s, int *P, real *b, 
		  enum xfe_AddType AddFlag, real *u);
/*
  PURPOSE: 

  Solves:  A^T*u = b, when A has been Block-PLU factored. 

  INPUTS:

  A : PLU n x n block matrix, each block is s x s
  n : A is an n x n block matrix
  s : each block is s x s
  P : stores permutation vector
  b : rhs
  AddFlag : whether to set, add-to, or subtract-from u


  OUTPUTS:

  u : is modified according to AddFlag with A^T\b


  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SolveBlockPLU_Matrix
extern int
xf_SolveBlockPLU_Matrix(real *A, int n, int s, int *P, real *B, 
			int m, enum xfe_AddType AddFlag, real *U);
/*
  PURPOSE: 

  Solves:   A*U = B, when A has been Block-PLU factored. 

  INPUTS:

  A : PLU n x n block matrix, each block is s x s
  n : A is an n x n block matrix
  s : each block is s x s
  P : stores permutation vector
  B : rhs matrix, of size n x m
  m : number of columns in B
  AddFlag : whether to set, add-to, or subtract-from u

  OUTPUTS:

  U : size n x m blocks; modified according to AddFlag with A\B

  RETURN:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_BlockPLUMxV
extern int
xf_BlockPLUMxV(real *A, int n, int s, int *P, real *u, 
	       enum xfe_AddType AddFlag, real *b);
/*
  PURPOSE: 

  Sets b = A*u (or +=/-= according to AddFlag) when A has been
  block-PLU factored.

  INPUTS:

  A : PLU n x n block matrix, each block is s x s
  n : A is an n x n block matrix
  s : each block is s x s
  P : stores permutation vector
  AddFlag : whether to set, add-to, or subtract-from u
  b : rhs

  OUTPUTS:

  u : is modified according to AddFlag with A*b

  RETURN:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_BlockPLUMTxV
extern int
xf_BlockPLUMTxV(real *A, int n, int s, int *P, real *u, 
		enum xfe_AddType AddFlag, real *b);
/*
  PURPOSE: 

  Sets b = A^T*u (or +=/-= according to AddFlag) when A has been
  block-PLU factored.

  INPUTS:

  A : PLU n x n block matrix, each block is s x s
  n : A is an n x n block matrix
  s : each block is s x s
  P : stores permutation vector
  AddFlag : whether to set, add-to, or subtract-from u
  b : rhs

  OUTPUTS:

  u : is modified according to AddFlag with A^T*b

  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_BlockMxV
extern int
xf_BlockMxV(real *A, int n, int s, int m, real *u, 
	    enum xfe_AddType AddFlag, real *b);
/*
  PURPOSE: 

  Sets b = A*u (or +=/-= according to AddFlag), where A can be
  rectangular.

  INPUTS:

  A : block matrix
  n,m : A is an n x m block matrix
  s : each block is s x s
  u : input vector, gets multiplied by A
  AddFlag : whether to set, add-to, or subtract-from u


  OUTPUTS:

  b : is modified according to AddFlag with A*u

  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_BlockMTxV
extern int
xf_BlockMTxV(real *A, int n, int s, int m, real *u, 
	     enum xfe_AddType AddFlag, real *b);
/*
  PURPOSE: 

  Sets b = A^T*u (or +=/-= according to AddFlag), where A can be
  rectangular.

  INPUTS:

  A : block matrix
  n,m : A^T is an n x m block matrix
  s : each block is s x s
  u : input vector, gets multiplied by A^T
  AddFlag : whether to set, add-to, or subtract-from u


  OUTPUTS:

  b : is modified according to AddFlag with A^T*u

  RETURN:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_BlockMxM
extern void
xf_BlockMxM(const real *A, int n, int s, int m, const real *U, 
	    int l, enum xfe_AddType AddFlag, real *B);
/*
  PURPOSE: 

  Sets B = A*U (or +=/-= according to AddFlag), where A is a block
  matrix, while U is not.  To make the multiplication work, think of
  replacing U(i,j) with eye(s)*U(i,j).

  INPUTS:

  A : block matrix
  n,m : A is an n x m block matrix
  s : each block is s x s
  U : input m x l matrix, gets multiplied by A
  l : U is m x l (not block)
  AddFlag : whether to set, add-to, or subtract-from B

  OUTPUTS:

  B : n x l block matrix; modified according to AddFlag with A*U

  RETURN:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_MxBlockM
extern void
xf_MxBlockM(const real *U, int l, int n, const real *A, int s,
	    int m, enum xfe_AddType AddFlag, real *B);
/*
  PURPOSE: 

  Sets B = U*A (or +=/-= according to AddFlag), where A is a block
  matrix, while U is not.  To make the multiplication work, think of
  replacing U(i,j) with eye(s)*U(i,j).

  INPUTS:

  A : block matrix
  n,m : A is an n x m block matrix
  s : each block is s x s
  U : input l x n matrix, multiplies A
  l : U is l x n (not block)
  AddFlag : whether to set, add-to, or subtract-from B

  OUTPUTS:

  B : n x l block matrix; modified according to AddFlag with U*A

  RETURN:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_MTxBlockM
extern void
xf_MTxBlockM(const real *U, int l, int n, const real *A, int s,
	     int m, enum xfe_AddType AddFlag, real *B);
/*
  PURPOSE: 

  Sets B = U^T*A (or +=/-= according to AddFlag), where A is a block
  matrix, while U is not.  To make the multiplication work, think of
  replacing U(i,j) with eye(s)*U(i,j).

  INPUTS:

  A : block matrix
  n,m : A is an n x m block matrix
  s : each block is s x s
  U : input n x l matrix, multiplies A
  l : U is n x l (not block)
  AddFlag : whether to set, add-to, or subtract-from B

  OUTPUTS:

  B : n x l block matrix; modified according to AddFlag with U*A

  RETURN:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_BlockMxBlockM
extern int
xf_BlockMxBlockM(const real *A, int n, int s, int m, const real *U, 
		 int l, enum xfe_AddType AddFlag, real *B);
/*
  PURPOSE: 

  Sets B = A*U (or +=/-= according to AddFlag), where A can be
  rectangular.

  INPUTS:

  A : block matrix
  n,m : A is an n x m block matrix
  s : each block is s x s
  U : input m x l block matrix, gets multiplied by A
  l : U is m x l blocks
  AddFlag : whether to set, add-to, or subtract-from B

  OUTPUTS:

  B : n x l block matrix; modified according to AddFlag with A*U

  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_ComputeBlockPLU
extern void
xf_BlockIdentity(int n, int s, real * RSTRCT A);
/*
  PURPOSE: 

  Returns a block identity matrix.

  INPUTS:

  n : A is an n x n block matrix
  s : each block is s x s
  
  OUTPUTS:

  A : n x n block identity matrix (on-diagonal blocks = I, all else is 0)

  RETURN: None

*/


/******************************************************************/
//   FUNCTION Purpose: xf_BlockOutProd_Add, _Sub
extern void
xf_BlockOutProd_Add(const real *u, const real *v, int n, int s, 
		    int m, real *A);
extern void
xf_BlockOutProd_Sub(const real *u, const real *v, int n, int s, 
		    int m, real *A);
/*
  PURPOSE: 

  Performs an outer product of u and v and adds the result to the
  block-format matrix A.

  INPUTS:

  u : an n x s block vector
  v : an m x s block vector
  n,m : A is an n x m block matrix
  s : each block is s x s

  OUTPUTS:

  A : block matrix modified with the outer product

  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_QRFactorHouseholder
extern int
xf_QRFactorHouseholder(real *A, int rA, int cA, real *Q, real *R);
/*
  PURPOSE: 

  Calculates a QR factorization of A, rA >= cA. Both Q (rA*rA) and
  R(rA*cA) must be pre-allocated. Q is returned as the full rA x rA
  matrix.  Only the first cA*cA entries of R are valid on return.

  INPUTS:

  A : rA*cA matrix
  rA, cA : # rows/columns in A matrix

  OUTPUTS:

  Q, R : factorization matrices as described above (must be preallocated)

  RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_EigSymmetricQR
extern int
xf_EigSymmetricQR(real *A, int r, real tol, int maxiter, 
		  real *EG, real *EV);
/*
  PURPOSE: 

  Calculates eigenvalues and eigenvectors of a real, symmetric matrix
  A.  Employs a QR factorization trick, whereby 

  A0 = Q0*R0, 
  A1 = R0*Q0 = Q1*R1,
  A2 = R1*Q1 = Q2*R2,
  ...

  AN approaches a diagonal matrix with eigs on the diagonal.
  Eigenvectors are then the columns of S = Q0*Q1*...*QN

  INPUTS:

  A : r*r matrix
  r : size of A matrix
  tol : tolerance to which off-diagonal entries of AN are checked
  maxiter : maximum number of iterations

  OUTPUTS:

  EG : vector of r eigenvalues, sorted from largest to smallest in mag
  EV : [r*r], rows are the corresponding eigenvectors

  RETURN:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_SVDGolubReinsch
extern int 
xf_SVDGolubReinsch(real *U, int m, int n, enum xfe_Bool NeedU, 
		   real *W, real *V);
/* 
   PURPOSE

   Computes Singular Value Decomposition of an m by n matrix U.

   Taken from:

   Singular Value Decomposition and Least Squares Solutions
   G. H. GOLUB and C. REINSCH
   Numer. Math. 14, 403--420 (1970)

   as implemented in:

   Algorithm 581, COLLECTED ALGORITHMS FROM ACM.
   ACM-TRANS. MATH. SOFTWARE, VOL.8, NO. 1, MAR., 1982, P. 84.
   TONY CHAN 
  
   Translated from FORTRAN into C by hand.

   INPUTS:

   U : an m by n matrix (this function changes U)
   m, n : size of A
   NeedU : if True, left singular vectors will be computed and
   returned in U

   OUTPUTS:

   U : left singular vectors if NeedU is True. Otherwise, U is used
   as temporary storage (data is altered).
   W : vector of n singular values (memory must be preallocated)
   V : columns are n right singular vectors; optional; only computed
   if V != NULL; memory must be preallocated.

   RETURNS:  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_QRBulgeChase
extern int 
xf_QRBulgeChase(real *shift, int nk, int np, real *alpha,
		real *beta, real *Q);
/* 
   PURPOSE

   Implicitly applies np shifts to a symmetric tridiagonal matrix,
   using QR factorization.  Implements:
  
   Q = I
   for j = 0:np-1,
   [Qj,R] = qr(H - shift(j)*I)
   H = Qj^T*H*Qj
   Q = Q*Qj
   end

   Translated from the ARPACK Fortran routine dsapps.f

   INPUTS:

   nk : size of square H matrix is np+nk
   np : number of shifts
   alpha : main diagonal of H
   beta  : subdiagonal of H
  
   OUTPUTS:

   Q : accumulated rotations.  Must be pre-allocated to hold
   (np+nk)*(np+nk) values
   alpha : modified diagonal of H
   beta  : modified sub-diagonal of H

   RETURNS:  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SortIntPos
extern int 
xf_SortIntPos(int *u, int n, int *pos, enum xfe_Bool InverseFlag);
/* 
PURPOSE

  Sorts n entries in u in ascending order and returns ranks in pos.

  Function modified from public-domain C implementation written by
  Darel Rex Finley

  In parallel, u is sorted in serial on each processor.

INPUTS:

  u : vector of integer entries to be sorted
  n : number of entries to be sorted
  InverseFlag : whether forward or inverse mapping is desired (see below)

OUTPUTS:

  pos : must be pre-allocated.  pos[i] stores the post-sorting rank of
  the ith (original) entry in u if InverseFlag == True.
  Otherwise, pos[i] stores the index of the rank i entry.

RETURNS:  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_SortRealParallel
extern int 
xf_SortRealParallel(real *u, int n, enum xfe_Bool SerialFlag, int *pos);
/* 
   PURPOSE

   Sorts n entries in u in ascending order and returns ranks in pos.
   Suitable for calling in parallel.

   INPUTS:

   u : vector of real entries to be sorted
   n : number of entries to be sorted
   SerialFlag : if true, force a serial sort (every proc sorts its own
                list indpendently)
  
   OUTPUTS:

   pos : must be pre-allocated.  pos[i] stores the post-sorting rank of
   the ith (original) entry in u.  The rank is over data on all
   procs, if parallel.

   RETURNS:  Error code

*/

// extern int 
// xf_SortRealParallel2(real *u, int n, enum xfe_Bool SerialFlag, int *pos, int nelemtot_glob, int *locorders, int *ElemOrder);

/******************************************************************/
//   FUNCTION Prototype: xf_QuadraticRoots
extern int 
xf_QuadraticRoots(real *coeff, int *nroots, real *roots);
/* 
   PURPOSE

   Applies the quadratic formula to the equation:

   coeff[0]*x^2 + coeff[1]*x + coeff[2] = 0.

   Returns the number of *real* roots, and their values, sorted in
   ascending order if two real roots exist.

   INPUTS:

   coeff : coefficients of quadratic (vector of 3 reals)
  
   OUTPUTS:

   nroots : number of real roots
   roots  : vector of real roots; valid entries are roots[0..nroots-1]

   RETURNS:  Error code

*/

/******************************************************************/
//   FUNCTION Definition: xf_SortInt
extern int
xf_SortInt(int n, int *a);
/* 
   PURPOSE
	
   Sort an integer array. Uses the Shell's method modified 
   from Numerical Recipes in C.
	
   INPUTS:
	
   n: number of entries in a
   a: integer array to be sorted
	
   OUTPUTS:
	
   None. "a" gets modified.
	
   RETURNS:  Error code
	
*/


/******************************************************************/
//   FUNCTION Prototype: xf_RandUniform
extern int 
xf_RandUniform(int n, real *pval);
/*
  PURPOSE:  

  Calculates n random real numbers betwen 0 and 1.  Parallel-safe (all
  procs will have same values).

  INPUTS: 
  
  n : number of values to calculate

  OUTPUTS: 

  (*pval) : vector of random real numbers

  RETURN: error code 
*/


/******************************************************************/
//   FUNCTION Prototype: xf_RandNormal
extern int 
xf_RandNormal(int n, real *pval);
/*
  PURPOSE:  

  Calculates n random real numbers drawn from a uniform distribution
  with mean 0 and variance 1.  Parallel-safe (all procs will have same
  values).

  INPUTS: 
  
  n : number of values to calculate

  OUTPUTS: 

  (*pval) : vector of random real numbers

  RETURN: error code 
*/

/******************************************************************/
//   FUNCTION Definition: xf_BinSearch
extern int xf_BinSearch(const int target, const int *set, 
                 int begin, int end, int *rank);
/*
 PURPOSE:  
 
 Performs a binary search for "target" in "set" ranging from "begin"
 to "end". 
 
 INPUTS: 
 
 target : integer to search for.
 set : ORDERED set of integers.
 begin, end: indices that define the range in "set"
 
 OUTPUTS: 
 
 (*rank): if NOT NULL will assume the position of "target" in "set".
 
 RETURN: error code 
*/


/******************************************************************/
//   FUNCTION Prototype: xf_RealNorm
extern real 
xf_RealNorm(const real * RSTRCT A, int r);

/*
 PURPOSE:  Calculates L2 norm of real-valued vector A 
 
 INPUTS: 
 
 A : pointer to vector
 r : rank of vector
 
 OUTPUTS: None
 
 RETURN: L2 norm of A
*/

/******************************************************************/
//   FUNCTION Definition: xf_Add2OrderedSet
extern int xf_Add2OrderedSet(const int entry, int *set_size, int **set,
                             int **orig_rank);
/*
 PURPOSE:  
 
 Adds "entry" to "set" such that "set" is kept in the crescent order. 
 
 INPUTS: 
 
 entry : integer to include to "set".
 set : ORDERED set of integers.
 orig_rank : if profided, it will store the original 
             order at which entry was received
 set_size: number of entries in "set".
 
 OUTPUTS: 
 
 set_size is updated.
 orig_rank is updated  
 set is updated
 
 RETURN: error code 
 */

#endif // end ifndef _xf_Math_h
