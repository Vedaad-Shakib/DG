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

#ifndef _xf_DataStruct_h
#define _xf_DataStruct_h 1

/*
  FILE:  xf_DataStruct.h

  This file contains the xflow Data structures

*/


#include "xf.h"
#include "xf_LineStruct.h"
#include "xf_LinearSolverStruct.h"

/* Data representation type */
enum xfe_DataType{

  xfe_Vector,
  /* Vector data with linkage, possibly interpolated.  Solution,
     residual, and adjoint vectors are stored as this type.  Other
     mesh-linked data (e.g. node flags, face values) can also be
     stored as this type. */

  xfe_Matrix,    
  /* Matrix data with linkage.  Element-specific matrices (e.g. Mass
     Matrix) can be stored as this type. */

  xfe_JacobianMatrix,
  /* The element-based Jacobian Matrix is stored as its own data type
     that spans all the element groups.  A compact (nearest-neighbor)
     stencil is required.
   */
  xfe_VectorSet,
  /* 
     Set of Vectors under one common data node.
   */
  xfe_DataOther,          /* Other data type */
  xfe_DataLast
};

/* corresponding names */
static char *xfe_DataName[xfe_DataLast] = {
  "Vector",
  "Matrix",
  "JacobianMatrix",
  "VectorSet",
  "DataOther"
};


 
/* Linkage Type for DataWithLinkage */
enum xfe_LinkageType{
  xfe_LinkageNode,
  xfe_LinkageElem,
  xfe_LinkageFace,
  xfe_LinkageIFace,
  xfe_LinkageBFace,
  xfe_LinkageGlobElem,
  xfe_LinkageElemGroup,
  xfe_LinkageBFaceGroup,
  xfe_LinkageNone,
  xfe_LinkageLast
};

/* corresponding names */
static char *xfe_LinkageName[xfe_LinkageLast] = {
  "Node",
  "Elem",
  "Face",
  "IFace",
  "BFace",
  "GlobElem",
  "ElemGroup",
  "BFaceGroup",
  "None"
};


/* Solver Role Type */
enum xfe_SolverRoleType{
  xfe_SolverRolePrimalState,
  xfe_SolverRolePrimalRes,
  xfe_SolverRoleAdjointState,
  xfe_SolverRoleAdjointRes,
  xfe_SolverRoleOther,
  xfe_SolverRoleNone,
  xfe_SolverRoleLast
};

/* corresponding names */
static char *xfe_SolverRoleName[xfe_SolverRoleLast] = {
  "PrimalState",
  "PrimalRes",
  "AdjointState",
  "AdjointRes",
  "Other",
  "None",
};


/* Array Parallel Info structure */
typedef struct
{

  enum xfe_Bool HaloFlag;
  /* If True, the associated array is part of a halo group.  This
     means that none of the array elements need to be sent to any
     other processors and that the entire contents of the data serve
     as receive buffers for data incoming from adjacent processors. */

  void *Request;
  /* Request[iProc] stores a communicator-specific handle that is
     created upon initiation of a Send (normal group) or Recv (halo
     group) and that can be queried with Wait to determine when
     complete. */

  int *nSendElem;   // Only valid if HaloFlag == False
  int **SendElem;   // Only valid if HaloFlag == False
  int *nRecvElem;   // Only valid if HaloFlag == True
  /* These are soft-link pointers to the corresponding ParallelInfo
     fields in MeshStruct (dereferenced by one pointer in that they do
     not include the [egrp]), declared here for convenience when
     dealing with element arrays.  See the descriptions in
     MeshStruct.h.  NOTE: soft-link pointers means that these arrays
     are not allocated or released. */

  int  **iSendBuf;
  real **rSendBuf;
  /* SendBuf[iProc] points to a buffer of size nsend*GenArray->r,
     where nsend is the number of r-size pieces of data to send to
     iProc (e.g nSendElem[iProc] for element arrays).  At most one of
     these buffers is not NULL, depending on GenArray->Size.  Halo
     group arrays do not need send buffers, and hence keep both of
     these equal to NULL. Variable orders are supported, in which
     case the data pieces need not be all of constant size r. */
} 
xf_ArrayParallelInfo;


/* GenArray stores integer or real data as a 2d array */
typedef struct
{
  enum xfe_SizeType Size;
  /* int or real */

  int n;
  /* first dimension */
  
  int r;
  /* constant second dimension, or max if varying */

  int *vr;
  /* varying second dimension, or NULL if second dim is constant */

  int  **iValue;
  real **rValue;
  /* array values: Value[j][k]
         j = [0 .. n-1]
         k = [0 .. r[j]-1]

     Only one of iValue or rValue is allocated, depending on the Size
     type.  The other pointer is set to NULL.  This approach is used
     to allow simpler dereferencing than with a void pointer.
  */

  xf_ArrayParallelInfo *ParallelInfo;
  /* In parallel runs, contains structure storing send buffers for
     each processor as well as communication request handles.  NULL
     for serial runs. */
  
}
xf_GenArray;



typedef struct
{
  enum xfe_LinkageType Linkage;
  /* type of Mesh structure to which the data is associated */
  
  int LinkageIndex;
  /* e.g. ElemGroup or BFaceGroup number */

  enum xfe_SolverRoleType SolverRole;
  /* primal state, adjoint state, residual, other, SolverRoleNone */

  int StateRank;
  /* Number of state variables */

  char **StateName;
  /* Names of state variables (if any) */

  char *OutputName;
  /* Name of output if this vector is associated with one (NULL otherwise). */
  
  int TimeIndex;
  /* For unsteady runs: e.g. 0 is current time level  */

  int MGIndex;
  /* For multigrid runs: e.g. 0 is fine MG level */

  enum xfe_Bool ParallelFlag;
  /* If True, this vector can exchange data with other processors.
     This means the GenArray structures have a ParallelInfo field
     defined. */

  enum xfe_Bool HaloInTransit;
  /* Not read or written.  Used during halo communication do signify
     that the vector is in transit. */

  int nArray;
  /* number of arrays of storage (e.g. 1 or nElemGroup) */

  int nArraySelf;
  /* number of arrays on own (self) processor.  Equal to nArray for
     serial runs.  For parallel element arrays, this will usually be
     nArray/2, to not count halos. */

  enum xfe_BasisType *Basis;
  /* Interpolation basis type vector [nArray]; NULL if data is not
     interpolated */

  int *Order;
  /* Interpolation order vector [nArray]: NULL, if data is not
     interpolated */

  int *nComp;
  /* number of components in each array; e.g. for variable-order
     interpolation.  NULL if not doing variable order. */

  int **vOrder;
  /* Array of variable interpolation orders [nArray][nComp[egrp]]:
     NULL, if not interpolated or constant order per group. */

  enum xfe_SizeType Size;
  /* Size type (int, real, etc.) used in the GenArrays */

  xf_GenArray *GenArray;
  /* Arrays of values [nArray] */
}
xf_Vector;



typedef struct
{
  enum xfe_LinkageType Linkage;
  /* type of other structure to which the data is associated */

  int LinkageIndex;
  /* e.g. ElemGroup or BFaceGroup number */

  int Order1, Order2;
  /* Interpolation orders associated with matrix (for a Mass matrix
     these would be the same, but for an inter-order transfer matrix,
     these would be different).*/

  enum xfe_BasisType Basis1, Basis2;
  /* Interpolation basis type; valid only if data is interpolated
     (i.e. if InterpOrder >= 0) */

  xf_GenArray *GenArray;
  /* Array of values [1] */
  
  int **P;
  /*
   When the matrix is PLU-factored. P stores the permutation indices.
   */
}
xf_Matrix;



/* JacobianMatrix stores a real-precision, element-based Jacobian
   matrix that includes, for each element, blocks corresponding to
   self and to neighbors across all faces.  Note, a Jacobian matrix is
   inherently linked with a mesh.*/
typedef struct
{
  int Preconditioner;
  /* If the Jacobian is preconditioned, this specifies the
     Preconditioner that was computed (e.g. the block diagonals are
     PLU factored for Block). */

  xf_Vector *U;
  /* Pointer to state vector; useful for memory-lean
     preconditioners. */

  xf_LineSet *LineSet;
  /* If a line preconditioner is used, this will point to a LineSet
     structure for the current set of lines.  See SolverStruct.h for a
     definition. */

  xf_LinearSolverILUData *ILUData;
  /* Used when ILU is the preconditioner.  Stores ordering and other
     temporary variables. */


  int ***P;
  /* When on-diagonal blocks are PLU factored, P stores the
     permutation vector for each block.
         P[egrp][elem][k] = permutation vector , 0 <= k < r[egrp]
   */

  int negrp;
  /* Number of element groups */

  int negrphalo;
  /* Number of element groups including halo.  Different from negrp if
     run is parallel. Need this info for de-allocation purposes. */

  int StateRank;
  /* State rank of the data */

  enum xfe_BasisType *Basis;
  /* vector of bases for each group [negrp] */

  int *Order;
  /* vector of Orders for each group [negrp] */

  int *nvec;
  /* Rank/StateRank for each element group.  The self block for each
     element in group egrp is of size nvec[egrp]*StateRank x
     nvec[egrp]*StateRank. The off-diagonal blocks could be
     rectangular, nvec[egrp1]*StateRank x nvec[egrp2]*StateRank, if
     different interpolation orders/bases are used for different
     element groups
  */

  int **vnvec;
  /* Variable-order version of nvec.  One rank/StateRank per element. */

  int ***egrpN;
  int ***elemN;
  int ***faceN;
  /* For Jacobian ease of use, these arrays store the adjacent elements
     for each element:
         egrpN[egrp][elem][face] = adjacent egrp
	 elemN[egrp][elem][face] = adjacent elem
	 faceN[egrp][elem][face] = local face number on adjacent elem
  */

  real ****Value;  
  /* array of pointers to values
     value[egrp][elem][1+face][k]:

     egrp = [ 0..nElemGroup-1]
     elem = [ 0..nElem(egrp)-1]
     face = [-1..nFace(elem)-1]
            -1 corresponds to self-block
     k    = [ 0..rvec[egrp]*rvec[egR]]
            egR = group of element neighbor across face
  */

  void *R_Uc;
  /* Coarse Jacobian for coarse-grid correction (or NULL) */

  enum xfe_Bool ProjectionNeeded;
  /* If True, R_Uc will be recalculated next time it is required. */

}
xf_JacobianMatrix;



/* A VectorSet is just a group of identically-sized vectors*/
typedef struct
{
  int nVector;
  /* Number vectors*/

  xf_Vector *Vector;
  /* Array of vectors */
}
xf_VectorSet;

/*------------- Data structure definition  --------------*/
struct xf_Data
{
  
  char *Title;
  /* Null-terminated string, not necessarily unique */

  enum xfe_DataType Type;
  /* Type of data */

  enum xfe_Bool ReadWrite;
  /* If true, data should be read/written */

  void *Data;
  /* Pointer to structure that stores the data */

  struct xf_Data *Prev, *Next;
  /* Pointer to previous and next data structures */
  
};
typedef struct xf_Data xf_Data;


/*------------- DataSet structure definition  --------------*/
struct xf_DataSet
{

  xf_Data *Head;
  xf_Data *Tail;
  /* Pointer to head and tail of linked list */
  
};
typedef struct xf_DataSet xf_DataSet;



#endif // end ifndef _xf_DataStruct_h
