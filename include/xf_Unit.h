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

#if !defined(_xf_Unit_h) && !defined(_UnitTest_h)
#define _xf_Unit_h 1

/*
  FILE:  xf_Unit.h

  This file contains macros and function prototypes for unit testing.

*/

#include "xf.h"
#include "xf_IO.h"
#include "xf_AllStruct.h"
#include "xf_UnitStruct.h"
#include <math.h>

#define UTOL0 1e-15
#define UTOL1 1e-14
#define UTOL2 1e-13
#define UTOL3 1e-12
#define UTOL4 1e-11
#define UTOL5 1e-10

// global variables
int unit_itemp;

#define xf_AssertEqual(a, b) \
  if ((a) != (b)) return xf_Error(xf_ASSERT_FAILED)

#define xf_AssertErrorOK(ierr)			       \
  if (ierr != xf_OK) {				               \
    if (ierr == xf_NEED_LAPACK){                               \
      xf_printf("This test requires LAPACK. Continuing.\n");   \
      return xf_OK;                                            \
    }                                                          \
    else                                                       \
     return xf_Error(xf_ASSERT_FAILED);			       \
  }

#define xf_AssertWithin(a, b, c)			       \
  if (fabs((a)-(b))>(c)) {				       \
    xf_printf("a = %.15f, b = %.15f\n", a, b);		       \
    return xf_Error(xf_ASSERT_FAILED);			       \
  }

#define xf_AssertIntVectorEqual(A, B, n)				\
  if (unit_itemp = AssertIntVectorEqual(A,B,n)){			\
    xf_printf("A[%d] = %d, B[%d] = %d\n", unit_itemp-1,			\
	      A[unit_itemp-1], unit_itemp-1, B[unit_itemp-1]);		\
    return xf_Error(xf_ASSERT_FAILED);					\
  }  

#define xf_AssertRealVectorWithin(A, B, n, c)				\
  if (unit_itemp = AssertRealVectorWithin(A,B,n,c)){			\
    xf_printf("A[%d] = %.15f, B[%d] = %.15f\n", unit_itemp-1,		\
	      (A)[unit_itemp-1], unit_itemp-1, (B)[unit_itemp-1]);	\
    return xf_Error(xf_ASSERT_FAILED);					\
  }  


/* Prototypes */

/******************************************************************/
//   FUNCTION Prototype: AssertIntVectorEqual
extern int
AssertIntVectorEqual(int *A, int *B, int n);
/*
PURPOSE: 

  Returns nonzero error if any A[i] != B[i],  0<i<n-1

INPUTS:

  A : int vector
  B : int vector
  n : vector length
  
OUTPUTS:

  None

RETURN:

  Zero or i+1 of nonmatching index

*/


/******************************************************************/
//   FUNCTION Prototype: AssertRealVectorWithin
extern int
AssertRealVectorWithin(real *A, real *B, int n, real c);
/*
PURPOSE: 

  Returns nonzero error if any fabs(A[i]-B[i]) > c,  0<i<n-1


  A : real vector
  B : real vector
  n : vector length
  
OUTPUTS:

  None

RETURN:

  Zero or i+1 of nonmatching index

  None
*/

/******************************************************************/
//   FUNCTION Prototype: xfu_PingPerturb
extern real
xfu_PingPerturb(real eps, int ieps);
/*
PURPOSE: 

  Returns perturbation for pinging
  
INPUTS:

  eps  : baseline perturbation
  ieps : index for order of convergence check

OUTPUTS: none

RETURN: perturbation = eps*2^(-ieps)

*/

/******************************************************************/
//   FUNCTION Prototype: xfu_PingCheckRate
extern int
xfu_PingCheckRate(real *veps, real tol, real meps);
/*
PURPOSE: 

  Checks if the rate of convergence in veps is at least 2.0-tol.
  
INPUTS:

  veps : two numbers = differences between FD and analytical for two epsilons
  tol  : convergence rate okay if in [2-tol, infty)
  meps : machine zero cutoff (rate is okay if veps[0] is below meps)

OUTPUTS: none

RETURN: error code, e.g. PING_FAILED

*/


/******************************************************************/
//   FUNCTION Prototype: xf_UnitrVector
extern int
xf_UnitrVector(int nArray, int *n, int *r, int **vr, xf_Vector **pV);
/*
PURPOSE: 

  Creates and allocates a real Vector.

INPUTS:

  nArray : desired number of arrays in the vector
  n      : number of elements in each array [nArray]
  r      : rank of data within each element [nArray]
  vr     : variable rank for each element (optional)
  
OUTPUTS:

  (*pV)  : created and allocated Vector

RETURN:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_UnitiVector
extern int
xf_UnitiVector(int nArray, int *n, int *r, int **vr, xf_Vector **pV);
/*
PURPOSE: 

  Creates and allocates an integer Vector.

INPUTS:

  nArray : desired number of arrays in the vector
  n      : number of elements in each array [nArray]
  r      : rank of data within each element [nArray]
  vr     : variable rank for each element (optional)
  
OUTPUTS:

  (*pV)  : created and allocated Vector

RETURN:

  Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_UnitrVectorInterp
extern int
xf_UnitrVectorInterp(int nArray, int *n, int sr, enum xfe_BasisType *Basis,
		     int *Order, int **vOrder, xf_Vector **pV);
/*
PURPOSE: 

  Creates and allocates a real interpolated Vector.

INPUTS:

  nArray : desired number of arrays in the vector
  n      : number of elements in each array [nArray]
  sr     : state rank
  Basis  : vector of basis values, one for each array
  Order  : vector of order values, one for each array
  vOrder : variable order 2D array (optional)

OUTPUTS:

  (*pV)  : created and allocated Vector

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Definition: xf_UnitrVectorSet
extern int
xf_UnitrVectorSet (int nVector, int nArray, int *n, int *r, 
		   xf_VectorSet **pVS);
/*
PURPOSE: 

  Creates and allocates a real Vector Set.  Effectively a wrapper for
  xf_UnitrVector.

INPUTS:

  nVector : desired number of vectors
  nArray  : desired number of arrays in each vector
  n       : number of elements in each array [nArray]
  r       : rank of data within each element [nArray]
  
OUTPUTS:

  (*pVS)  : created and allocated Vector Set

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_UnitTwoQ1TriangleMesh
extern int
xf_UnitTwoQ1TriangleMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitTwoQ1TriangleAll(xf_All **pAll);
/*
PURPOSE: 

  Returns a two element mesh on the unit square.  The elements are Q1
  triangles.  All version returns an All structure with the mesh.

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_UnitTwoQ2TriangleMesh
extern int
xf_UnitTwoQ2TriangleMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitTwoQ2TriangleAll(xf_All **pAll);
/*
PURPOSE: 

  Returns a two element mesh on the unit square.  The elements are Q2
  triangles. All version returns an All structure with the mesh.

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_Unit*Tet*
extern int
xf_UnitQ1TetMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitQ1TetAll(xf_All **pAll);
extern int
xf_UnitQ2TetMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitQ2TetAll(xf_All **pAll);
/*
PURPOSE: 

  Returns a one element tetrahedron mesh.  The elements are Q1 or Q2
  depending on the call.  All version returns an All structure with
  the mesh.

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Unit*Quad*
extern int
xf_UnitQ1QuadMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitQ1QuadAll(xf_All **pAll);
extern int
xf_UnitQ2QuadMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitQ2QuadAll(xf_All **pAll);
/*
PURPOSE: 

  Returns a one element quadrilateral mesh.  The elements are Q1 or Q2
  depending on the call.  All version returns an All structure with
  the mesh.

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_Unit*Hex*
extern int
xf_UnitQ1HexMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitQ1HexAll(xf_All **pAll);
extern int
xf_UnitQ2HexMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitQ2HexAll(xf_All **pAll);
/*
PURPOSE: 

  Returns a one element hexahedron mesh.  The elements are Q1 or Q2
  depending on the call.  All version returns an All structure with
  the mesh.

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_UnitBoxQ1TriangleMesh
extern int
xf_UnitBoxQ1TriangleMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitBoxQ1TriangleAll(xf_All **pAll, enum xfe_MotionType MotionType);
extern int
xf_UnitBoxQ1TriangleAllEuler(xf_All **pAll, enum xfe_MotionType MotionType);
extern int
xf_UnitBoxQ1TriangleAllCNS(xf_All **pAll, enum xfe_MotionType MotionType);
extern int
xf_UnitBoxQ1TriangleAllCNS_SA(xf_All **pAll, enum xfe_MotionType MotionType);
/*
PURPOSE: 

  Returns an eight element mesh of a square.  The elements are Q1
  triangles.  All version return an All structure with the mesh, with
  filled in Param and an EqnSet Library (Scalar, Euler, CNS, etc.).

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  MotionType : enumerated type specifying motion, see xf_UnitStruct.h
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_UnitBoxQ2Triangle
extern int
xf_UnitBoxQ2TriangleMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitBoxQ2TriangleAllEuler(xf_All **pAll);
extern int
xf_UnitBoxQ2TriangleAllCNS(xf_All **pAll);
/*
PURPOSE: 

  Returns a two-element mesh of a square.  The elements are Q2 deformed
  triangles.  All version returns an All structure with the mesh, with
  filled in Param and an Euler or CNS EqnSet library.

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/




/******************************************************************/
//   FUNCTION Prototype: xf_UnitQ1Q2HexAllCNS
extern int
xf_UnitQ1Q2HexMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitQ1Q2HexAllCNS(xf_All **pAll);
/*
PURPOSE: 

  Returns a two-element mesh of a 3D box.  The elements are a Q2
  deformed hex adjacent to a Q1 hex.  An All structure or just the
  mesh, is returned with the mesh, with filled in Param and a CNS
  EqnSet library.

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_UnitBoxQ1Quad9Mesh
extern int
xf_UnitBoxQ1Quad9Mesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitBoxQ1Quad9All(xf_All **pAll);
/*
PURPOSE: 

  Returns a nine-element, Q1, quad mesh of a square. 

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_UnitBoxQ1Quad2EG9Mesh
extern int
xf_UnitBoxQ1Quad2EG9Mesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitBoxQ1Quad2EG9All(xf_All **pAll);
/*
PURPOSE: 

  Returns a nine-element, Q1, quad mesh of a square. The elements are
  divided into two groups (4 elements in the first, 5 in the second).

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/



/******************************************************************/
//   FUNCTION Prototype: xf_UnitBoxQ1Hex27Mesh
extern int
xf_UnitBoxQ1Hex27Mesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitBoxQ1Hex27All(xf_All **pAll);
/*
PURPOSE: 

  Returns a 27-element, Q1, hex mesh of a box. 

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_UnitBoxQ1Hex2EG27Mesh
extern int
xf_UnitBoxQ1Hex2EG27Mesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag);
extern int
xf_UnitBoxQ1Hex2EG27All(xf_All **pAll);
/*
PURPOSE: 

  Returns a 27-element, Q1, hex mesh of a box. The elements are
  divided into two groups (13 elements in the first, 14 in the
  second).

INPUTS:

  CreateFlag : if True, (*pMesh) will be allocated with a call to 
               CreateMesh.
  
OUTPUTS:

  (*pMesh) : allocated/created mesh structure
  (*pAll)  : allocated/created all structure

RETURN:

  None
*/

extern int
xf_UnitBoxQ1Quad9EulerPenaltyAll(xf_All **pAll);
/*
 PURPOSE: 
 
 Returns a nine-element, Q1, quad mesh of a square
 with a EqnSet using the penalty function. 
 
 INPUTS:
 
 CreateFlag : if True, (*pMesh) will be allocated with a call to 
 CreateMesh.
 
 OUTPUTS:
 
 (*pMesh) : allocated/created mesh structure
 (*pAll)  : allocated/created all structure
 
 RETURN:
 
 None
 */



/******************************************************************/
//   FUNCTION Prototype: xf_UnitIntervalAllEuler
extern int
xf_UnitIntervalAllEuler(xf_All **pAll);
/*
 PURPOSE: 
 
 Returns an 8-element, Q1 mesh of an interval with an EqnSet
 using the Euler equations.
 
 INPUTS: None
 
 OUTPUTS:
 
   (*pAll)  : allocated/created all structure
 
 RETURN: None
*/


/******************************************************************/
//   FUNCTION Prototype: xf_UnitSkewedQ1Tri5AllEuler
extern int
xf_UnitSkewedQ1Tri5AllEuler(xf_All **pAll);
/*
 PURPOSE: 
 
 Returns a 5-element, Q1, triangular mesh of a box with an EqnSet
 using the Euler equations.
 
 INPUTS: None
 
 OUTPUTS:
 
   (*pAll)  : allocated/created all structure
 
 RETURN: None
*/


#endif // end ifndef _xf_Unit_h

