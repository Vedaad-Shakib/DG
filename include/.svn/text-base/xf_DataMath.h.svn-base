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

#ifndef _xf_DataMath_h
#define _xf_DataMath_h 1

/*
  FILE:  xf_DataMath.h

  This file contains the headers for functions dealing with parameters.

*/



/******************************************************************/
//   FUNCTION Prototype: xf_SetVector
extern int 
xf_SetVector( xf_Vector *B, enum xfe_AddType AddFlag, xf_Vector *A);
/*
PURPOSE:

  Sets  A @= B, where @= is one of {=, +=, -=, =-}

INPUTS:

  B: Vector
  AddFlag : see above

OUTPUTS: 

  A: Resultant Vector
  
RETURN:

  xf_INCOMPATIBLE if A and B are not compatible
  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_SetZeroVector
extern int 
xf_SetZeroVector( xf_Vector *A);
/*
PURPOSE:

  Sets A = 0;

INPUTS:

  A: Vector
  
OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_SetZeroHaloVector
extern int 
xf_SetZeroHaloVector( xf_Vector *A);
/*
PURPOSE:

  Sets A = 0 on Halo elements only

INPUTS:

  A: Vector
  
OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_SetConstVector
extern int 
xf_SetConstVector( xf_Vector *A, int ival, real rval);
/*
PURPOSE:

  Sets A = ival for integer vectors, rval for real vectors

INPUTS:

  A: Vector
  
OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_VectorAbs
extern int 
xf_VectorAbs( xf_Vector *A);
/*
PURPOSE:

  Sets A to the absolute-value of A, abs(A), component-wise.

INPUTS:

  A: Vector
  
OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_VectorMin
extern int 
xf_VectorMin( xf_Vector *A, int *imin, real *rmin);
/*
PURPOSE:

  Computes minimum value of A

INPUTS:

  A: Vector
  
OUTPUTS: 

  (*imin) : minimum integer value if integer vector
   - OR -
  (*rmin) : minimum real value if real vector
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_VectorMask
extern int 
xf_VectorMask( xf_Vector *A, xf_Vector *Mask, int MaskVal);
/*
PURPOSE:

  Components of A for which Mask is not equal MaskVal are set to zero

INPUTS:

  A    : real or integer Vector
  Mask : integer Mask vector
  MaskVal : mask value -- components with this as Mask are left alone
  
OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_VectorDot
extern int 
xf_VectorDot( const xf_Vector *A, const xf_Vector *B, real *pdp);
/*
PURPOSE:

  Sets dp = dot(A,B);

INPUTS:

  A, B: Vectors

OUTPUTS: 

  (*pdp): dot product
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_VectorMult
extern int 
xf_VectorMult( xf_Vector *A, real c);
/*
PURPOSE:

  Sets A *= c;

INPUTS:

  A: Vector
  c: real number

OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_VectorVectorMult
extern int 
xf_VectorVectorMult( const xf_Vector *B, xf_Vector *A);
/*
PURPOSE:

  Sets A *= B;

INPUTS:

  B: Vector
  A: Vector

OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_VectorInv
extern int 
xf_VectorInv( xf_Vector *A);
/*
PURPOSE:

  Sets A = 1/A;

INPUTS:

  A: Vector

OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_VectorAdd
extern int 
xf_VectorAdd( xf_Vector *A, real c);
/*
PURPOSE:

  Sets A += c;

INPUTS:

  A: Vector
  c: real number

OUTPUTS: 

  A: Modified Vector
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_VectorMultSet
extern int 
xf_VectorMultSet( xf_Vector *B, real c, enum xfe_AddType AddFlag, xf_Vector *A);
/*
PURPOSE:

  Sets  A &= c*B, where &= is one of {=, +=, -=, =-}

INPUTS:

  B: Vector
  c: real number
  AddFlag : see above

OUTPUTS: 

  A: Resultant Vector
  
RETURN:

  xf_INCOMPATIBLE if A and B are not compatible
  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_VectorNorm
extern int 
xf_VectorNorm( const xf_Vector *A, const int p, real *pnorm);
/*
PURPOSE:

  Computes p-norm of real vector in A.  In parallel, uses only non-halo
  groups for ElemGlob vectors and reduces rnorm over all procs

INPUTS:

  A : Vector
  p : norm type to compute (0, 1, or 2) ... (0 = sum of entries, no abs)
  
OUTPUTS: 

  (*pnorm) : resulting norm
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Definition: xf_VectorStateNorms
extern int 
xf_VectorStateNorms(xf_All *All, xf_Vector *A, real *StateNorms,
                    int p);
/*
 PURPOSE:
 
 Computes p-norm of a interpolated real vector A for 
 each state quantity that A entries represents.
 
 INPUTS:
 
 A : Vector
 p : norm type to compute (0, 1, or 2) ... (0 = sum of entries, no abs)
 
 OUTPUTS: 
 
 *StateNorms : array of size StateRank with the p-norms of each 
               state entry. Must be pre-allocated
 
 RETURN:
 
 Error Code
 */

/******************************************************************/
//   FUNCTION Prototype: xf_VectorRand
extern int 
xf_VectorRand(xf_Vector *A, int PseudoStep);
/*
PURPOSE:

  Sets entries of vector A to random real numbers between 0 and 1.  If
  PseudoStep > 0, random numbers are not used; rather the entries of A
  are set incrementally between 0 and 1 in intervals of PseudoStep.

INPUTS:

  V : Vector to set to random
  PseudoStep : if > 0, random values are not used (as described above)
  
OUTPUTS: 

  V : entries set to random

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MatrixInvert
extern int 
xf_MatrixInvert( xf_Matrix *M, int nn, xf_Matrix *iM);
/*
PURPOSE:

  Inverts nn x nn array(s) in M and stores it (them) in iM

INPUTS:

  M : Matrix to invert
  nn : size of M (see above)
  
OUTPUTS: 

  iM : Matrix storing inverse of each array in M.  Data for iM
       must be allocated before the call.
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ProjectVector
extern int 
xf_ProjectVector( xf_All *All, const xf_Vector *A, 
		  enum xfe_Bool TransposeFlag, xf_Vector *B);
/*
PURPOSE:

  Prolongates/restricts all data in A to B according to the
  (Basis,Order) values in A and B. The prolongation / restriction is
  performed in the element interpolation ref. space.  If TransposeFlag
  is True, the transpose of the reverse operator is used.

INPUTS:

  All : All structure (for quickly pulling off transfer matrices)
  A  : source vector
  B  : destination vector
  TransposeFlag : if True, transpose of reverse operator is used;
                  for example, useful when residual is being restricted
		  by the transpose of the prolongation operator.
  
OUTPUTS: 

  None.  Data in B is altered
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectVectorInPlace
extern int 
xf_ProjectVectorInPlace( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
			 enum xfe_BasisType *BasisVec, enum xfe_BasisType BasisScal, 
			 xf_Vector *VOrder, int *OrderVec, int OrderScal,
			 int OrderIncrement);
/*
PURPOSE:

  Prolongates/restricts all data in V to (Basis,Order).  The Basis and
  Order in V can be general.  

  BasisVec specifies the basis in each element group; if passed in as
  NULL, BasisScal is used for all groups.  If BasisScal ==
  xfe_BasisLast, the existing basis is used.

  The order can be specified in a variety of ways.  If the element
  integer vector VOrder is provided, with order information for each
  individual element, this is used.  If VOrder is NULL, OrderVec is
  used to specify the desired order in each element group.  If
  OrderVec is NULL, OrderScal is used for all elements.  If OrderScal
  < 0, then the Order Increment is applied to the existing order.  If
  OrderIncrement < 0, then the existing order is used.

  The prolongation / restriction is performed in the element
  interpolation ref. space.

INPUTS:

  Mesh : optional; if provided, input BasisVec and OrderVec are
         accessed modulo negrp, in case halo groups are
         traveresed. Safest to do this if passing in BasisVec and
         OrderVec as not NULL
  DataSet : if provided, DataSet is searched for existing transfer
            matrices, and any computed matrices are stored here
  V : Vector
  BasisVec  : desired basis vector
  BasisScal : desired basis (used if BasisVec == NULL)
  VOrder    : integer vector of element-wise desired orders
  OrderVec  : desired order vector (used if VOrder == NULL)
  OrderScal : desired order (used if OrderVec == NULL)
  OrderIncrement: desired order increment, used if OrderVec == NULL,
                  and if OrderScal < 0
  
OUTPUTS: 

  None.  Data in V is altered if for any i up to V->nArray, 
         V->Basis[i] != Basis or V->Order[i] != Order 
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectVectorInPlace_Basis
extern int 
xf_ProjectVectorInPlace_Basis( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
			       enum xfe_BasisType *BasisVec, enum xfe_BasisType BasisScal);
/*
PURPOSE:

  Projects vector V to a new basis

INPUTS:

  Mesh : optional (useful in parallel if basis/order not stored on halo in U)
  DataSet : if provided, DataSet is searched for existing transfer
            matrices, and any computed matrices are stored here
  BasisVec  : desired basis vector
  BasisScal : desired basis (used if BasisVec == NULL)
  
OUTPUTS: 

  None.  Data in V is altered, as is order information for V.
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectVectorInPlace_OrderSet
extern int 
xf_ProjectVectorInPlace_OrderSet( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
				  enum xfe_BasisType *BasisVec, enum xfe_BasisType BasisScal, 
				  int **vOrder, int *OrderVec, int OrderScal);
/*
PURPOSE:

  Projects vector V using the provided Basis information and requested
  order information.  Use BasisVec=Null, BasisScal = xfe_BasisLast to
  keep the same basis as already in V.

INPUTS:

  Mesh : optional (useful in parallel if basis/order not stored on halo in U)
  DataSet : if provided, DataSet is searched for existing transfer
            matrices, and any computed matrices are stored here
  BasisVec  : desired basis vector
  BasisScal : desired basis (used if BasisVec == NULL)
  vOrder    : desired elemental order vector
  OrderVec  : desired order vector (used if vOrder == NULL)
  OrderScal : desired order (used if OrderVec == NULL)
  
OUTPUTS: 

  None.  Data in V is altered, as is order information for V.
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ProjectVectorInPlace_OrderIncrement
extern int 
xf_ProjectVectorInPlace_OrderIncrement( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
					enum xfe_BasisType *BasisVec, 
					enum xfe_BasisType BasisScal, int OrderIncrement);
/*
PURPOSE:

  Projects vector V using the provided Basis information and
  OrderIncrement.  Use BasisVec=Null, BasisScal = xfe_BasisLast to
  keep the same basis as already in V.

INPUTS:

  Mesh : optional (useful in parallel if basis/order not stored on halo in U)
  DataSet : if provided, DataSet is searched for existing transfer
            matrices, and any computed matrices are stored here
  BasisVec  : desired basis vector
  BasisScal : desired basis (used if BasisVec == NULL)
  OrderIncrement: desired order increment (can be negative)
  
OUTPUTS: 

  None.  Data in V is altered, as is order information for V.
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ProjectVectorInPlace_Vector
extern int 
xf_ProjectVectorInPlace_Vector( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
				xf_Vector *U);
/*
PURPOSE:

  Projects vector V using the basis/order information in U as a guide.

INPUTS:

  Mesh : optional (useful in parallel if basis/order not stored on halo in U)
  DataSet : if provided, DataSet is searched for existing transfer
            matrices, and any computed matrices are stored here
  U : provides basis/order information 
  
OUTPUTS: 

  None.  Data in V is altered if for any i up to V->nArray, 
         V->Basis[i] != desired Basis or V->Order[i] != desired Order 
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ProjectVectorInPlace_VOrder
extern int 
xf_ProjectVectorInPlace_VOrder( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *V,
				enum xfe_BasisType *BasisVec, 
				enum xfe_BasisType BasisScal, xf_Vector *VOrder);
/*
PURPOSE:

  Projects vector V using the given basis info and order information
  stored in the integer vector VOrder

INPUTS:

  Mesh : optional (useful in parallel if basis/order not stored on halo in U)
  DataSet : if provided, DataSet is searched for existing transfer
            matrices, and any computed matrices are stored here
  BasisVec  : desired basis vector
  BasisScal : desired basis (used if BasisVec == NULL)
  VOrder    : integer vector of element-wise desired orders

  
OUTPUTS: 

  None.  Data in V is altered if for any i up to V->nArray, 
         V->Basis[i] != desired Basis or V->Order[i] != desired Order 
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototpye: xf_ProjectVectors_VOrderFile
extern int 
xf_ProjectVectors_VOrderFile( xf_Mesh *Mesh, xf_DataSet *DataSet, int nVector,
			      xf_Vector **Vi, const char *fname);
/*
PURPOSE:

  Projects vectors Vi using the basis information provided and order
  information from a file containing a dataset containing an integer
  vector specifying the desired order.

INPUTS:

  Mesh : optional (useful in parallel if basis/order not stored on halo in U)
  DataSet : if provided, DataSet is searched for existing transfer
            matrices, and any computed matrices are stored here
  nVector : number of vectors to project
  fname : file name containing integer vector of orders
  
OUTPUTS: 

  None.  Data in Vi vectors is altered if for any i up to Vi->nArray, 
         Vi->Basis[i] != desired Basis or Vi->Order[i] != desired Order 
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_ProjectJacobian
extern int 
xf_ProjectJacobian( xf_All *All, const xf_JacobianMatrix *R_UA,
		    xf_JacobianMatrix *R_UB);
/*
PURPOSE:

  Prolongates/restricts Jacobian from A to B according to the
  (Basis,Order) values in A and B. The prolongation / restriction is
  performed in the element interpolation ref. space.

INPUTS:

  All : All structure (for quickly pulling off transfer matrices)
  R_UA  : source Jacobian
  R_UB  : destination Jacobian
  
OUTPUTS: 

  None.  Data in B is altered
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_VectorSetSVD
extern int 
xf_VectorSetSVD( xf_VectorSet *VS, int n, real *W);
/*
PURPOSE:

  Computes the SVD of the matrix formed by the columns of the first n
  vectors in VS.  The singular values are stored in W, while the
  vectors in VS are modified to contain the left singular vectors.

INPUTS:

  VS : vector set
  n  : number of vectors in VS to consider
  
OUTPUTS: 

  VS : vectors become left singular vectors
  W  : singular values [n]
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_VectorSetPOD
extern int 
xf_VectorSetPOD( xf_VectorSet **pSnapSet, int nSnap, int nBasis,
		 enum xfe_Verbosity Verbosity);
/*
PURPOSE:

  Wrapper for VectorSetSVD.  After computing the SVD, orders the
  vectors in SnapSet so that the first nBasis vectors with the largest
  singular values appear first.  Debug output is printed according to
  Verbosity.

INPUTS:

  pSnapSet : pointer to snapshot set -- SVD will be calculated of this
  nSnap : number of snapshots
  nBasis : number of basis functions ( < nSnap)
  Verbosity : set to high to print debug output info
  
OUTPUTS: 

  (*pSnapSet) : first nBasis vectors contain the POD vectors, ordered
                along singular values (from largest to smallest)
                
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_VectorSetPODMass
extern int 
xf_VectorSetPODMass(xf_All *All, const xf_VectorSet *SnapSet, 
		    int nSnap, int nBasis, xf_VectorSet *BasisSet);
/*
PURPOSE:

  Computes Mass-matrix norm POD of a snapshot set by constructing a
  matrix = SnapSet^T*M*SnapSet and then calculating its
  eigenvectors. Warning, this is not a well-conditioned way to
  calculate the POD basis.  Use the function based on the SVD in
  DataMath instead.  This is only for testing, although probably ok if
  nBasis << nSnap (only interested in the highest-energy eigenvalues)

  If All is NULL, Mass matrix will not be used


INPUTS:

  All : all structure; required for mass matrix
  SnapSet : snapshot set -- POD will be calculated of this
  nSnap : number of snapshots
  nBasis : number of basis functions ( < nSnap)
  
OUTPUTS: 

  BasisSet : basis vector set
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_MatrixFrobNorm
extern int 
xf_MatrixFrobNorm(real *M, int n, int m, real *Fnorm);
/*
 PURPOSE:
 
 Computes the Frobenius norm of the matrix M
 
 
 INPUTS:
 
 M : matrix
 n : number of rows
 m : number of collumns
 
 OUTPUTS: 
 
 Fnorm : value of the Frobenius norm
 
 RETURN:
 
 Error Code
 */



/******************************************************************/
//   FUNCTION Prototype: xf_ConvertVectorFromLagrange
extern int 
xf_ConvertVectorFromLagrange(xf_Vector *U);
/*
 PURPOSE:
 
   Treats existing numbers in U as coefficients on Lagrange basis
   functions, regardless of what U->Basis is.  Converts these numbers
   to the actual basis in U->Basis such that the reprsented function
   is unchanged.

   Useful when we want to, for example, initialize a field vector to a
   constant, for an arbitrary basis.
 
 INPUTS:
 
   U : vector to convert
 
 OUTPUTS: 
 
   None; data in U is changed
 
 RETURN:
 
 Error Code
*/



#endif // end ifndef _xf_DataMath_h
