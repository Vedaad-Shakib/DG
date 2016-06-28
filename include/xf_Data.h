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

#ifndef _xf_Data_h
#define _xf_Data_h 1

/*
  FILE:  xf_Data.h

  This file contains the headers for functions dealing with data.

*/

#include "xf_DataParallel.h"
#include "xfYu_Model.h"
#include "xf_IO.h"

extern int 
xfYu_FindOrCreatePrimalState( xf_All *All, Yu_Model Model, enum xfe_Bool RestartFlag, xf_ICs *ICsOrig,
                             xf_Vector **pU);

/******************************************************************/
//   FUNCTION Prototype: xf_InterpOrder
extern int
xf_InterpOrder( const xf_Vector *U, int egrp, int elem);
/*
PURPOSE:

  Calculate order of an element in a vector, allowing for the
  possibility of variable orders.

INPUTS:

  U : vector
  egrp, elem : element in question (elem only required if variable)

OUTPUTS: None

RETURN:

  U->vOrder[egrp][elem] if vOrder exists; else U->Order[egrp]

*/

/******************************************************************/
//   FUNCTION Prototype: xf_InitVector
extern void 
xf_InitVector(xf_Vector *V);
/*
PURPOSE:

  Initializes a vector to default 0/NULL values.

INPUTS:

  V : pointer to Vector

OUTPUTS: 

  None. V is initialized

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_CreateVector
extern int 
xf_CreateVector(xf_Vector **pV);
/*
PURPOSE:

  Creates a vector and initializes it

INPUTS:

  pV : pointer to Vector

OUTPUTS: 

  (*pV) is created and initialized

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyVector
extern int 
xf_DestroyVector(xf_Vector *V, enum xfe_Bool DestroySelf);
/*
PURPOSE:

  Destroys a Vector structure. 

INPUTS:

  V : pointer to Vector
  DestroySelf : if True, the pointer to V is released as well.

OUTPUTS: 

  None. V is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CreateMatrix
extern int 
xf_CreateMatrix( xf_Matrix **pM);
/*
PURPOSE:

  Creates a Matrix structure

INPUTS:

  pM : pointer to Matrix

OUTPUTS: 

  None: Matrix is allocated

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyMatrix
extern int 
xf_DestroyMatrix( xf_Matrix *M);
/*
PURPOSE:

  Destroys a Matrix structure. 

INPUTS:

  M : pointer to Matrix

OUTPUTS: 

  None. M is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CreateVectorSet
extern int 
xf_CreateVectorSet(int nVector, xf_VectorSet **pVS);
/*
PURPOSE:

  Creates a VectorSet structure. 

INPUTS:

  nVector : number of vectors that VS will hold

OUTPUTS: 

  pVS : pointer to created VectorSet

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_TrimVectorSet
extern int 
xf_TrimVectorSet(xf_VectorSet **pVS, int nVectorNew);
/*
PURPOSE:

  Trims number of vectors in (*pVS) to be nVectorNew.

INPUTS:

  pVS : pointer to VectorSet
  nVectorNew : desired number of vectors

OUTPUTS: 

  None. (*pVS) is trimmed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyVectorSet
extern int 
xf_DestroyVectorSet(xf_VectorSet *VS);
/*
PURPOSE:

  Destroys a VectorSet structure. 

INPUTS:

  VS : pointer to VectorSet

OUTPUTS: 

  None. VS is destroyed

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_CreateDataSet
extern int 
xf_CreateDataSet( xf_DataSet **pDataSet);
/*
PURPOSE:

  Creates a DataSet structure for storing general types of data.  A
  doubly-linked list is initialized with Head = Tail = NULL.

INPUTS:

  pDataSet : pointer to DataSet

OUTPUTS: 

  None: DataSet is allocated

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyDataSet
extern int 
xf_DestroyDataSet( xf_DataSet *DataSet);
/*
PURPOSE:

  Destroys all Data in a DataSet structure.  Memory is released for
  DataSet itself.

INPUTS:

  DataSet : pointer to DataSet

OUTPUTS: 

  None: DataSet and all of its Data are destroyed

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyDataInSet
extern int 
xf_DestroyDataInSet( xf_DataSet *DataSet, xf_Data *Data);
/*
PURPOSE:

  Destroys Data in DataSet while keeping connectivity valid (prev and
  next pointers)

INPUTS:

  Data : pointer to a piece of data in DataSet

OUTPUTS: 

  None: Data in DataSet is destroyed

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DataSetAdd
extern int 
xf_DataSetAdd( xf_DataSet *DataSet, const char Title[], const enum xfe_DataType Type,
	       const enum xfe_Bool ReadWrite, const void *Data, xf_Data **pD);
/*
PURPOSE:

  Adds a Data node to the DataSet linked list.

INPUTS:

  DataSet : DataSet to which the Data should be added
  Title   : Title ascribed to this Data node
  Type    : Type of Data
  ReadWrite : True if this Data node should be read/written
  Data    : void pointer to actual data
  

OUTPUTS: 

  pD: (optional) pointer to newly-created xf_Data node

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DataSetRemove
extern int 
xf_DataSetRemove(xf_DataSet *DataSet, const char Title[], enum xfe_Bool RootFlag);
/*
PURPOSE:

  Removes Data node(s) from the DataSet linked list for which the
  Data->Title maches the given Title.

INPUTS:

  DataSet  : DataSet from which the node(s) should be removed
  Title    : All Data nodes with matching Title will be removed
  RootFlag : If this is true, all data nodes whose Data->Title strings begin 
             with Title are removed.  For example, if RootFlag = "State", 
	     a node with Data->Title = "StateFoo" would be removed.
	     If False, an exact match is required for removal.  

OUTPUTS: 

  None (Data from DataSet is removed)

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_FindDataByTitle
extern int 
xf_FindDataByTitle( xf_DataSet *DataSet, const char Title[], 
		    enum xfe_DataType DataType, xf_Data **pD);
/*
PURPOSE:

  Locates data of title = Title, and type == DataType, in DataSet.

INPUTS:

  DataSet  : DataSet to search
  Title    : Name of Data node to look for
  DataType : Whether to look for a Vector, Matrix, etc.

OUTPUTS: 

  (*pD)    : Data node, if found

RETURN:

  Error Code: 
    NOT_FOUND: data named Title does not exist, 
    MULTIPLE_MATCHES: multiple data nodes are named Title
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DataSetMerge
extern int 
xf_DataSetMerge( xf_DataSet *DataSetFrom, xf_DataSet *DataSetTo);
/*
PURPOSE:

  Adds a Data node to the DataSet linked list.

INPUTS:

  DataSetFrom : DataSet from which to copy data.  The data here is destroyed.
  
OUTPUTS: 

  DataSetTo   : DataSet to which to copy data.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_AllocGenArray
extern int 
xf_AllocGenArray( xf_GenArray *ga);
/*
PURPOSE:

  Allocates a general array of size ga->n*ga->r, or with variable r
  values if ga->vr is not NULL.  In the latter case, ga->n*sum(ga->vr)
  values are allocated.

INPUTS:

  ga->Size : int or real
  ga->n    : first array dimension
  ga->r    : second array dimension, or max if variable
  ga->vr   : variable second array dimension (a vector), or NULL

OUTPUTS: 

  ga : with allocated general array pointer

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CopyVector
extern int 
xf_CopyVector(xf_Mesh *Mesh, xf_Vector *U, xf_Vector *V);
/*
PURPOSE:

  Copies vector U to V.  Parallel communication send and receive
  buffers are not copied.

INPUTS:

  Mesh : mesh structure; can be NULL but should be passed in if 
         parallel communication structure is to be copied
  U : input vector
  
OUTPUTS: 

  V : copied vector
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DuplicateVector
extern int 
xf_DuplicateVector(xf_Mesh *Mesh, xf_Vector *U, xf_Vector **pV);
/*
PURPOSE:

  Wrapper for CopyVector.  Creates (*pV) first.

INPUTS:

  Mesh : mesh structure; can be NULL but should be passed in if 
         parallel communication structure is to be copied
  U : input vector
  
OUTPUTS: 

  (*pV) : copied vector (must not be allocated on input)
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_NewGREVector
extern int 
xf_NewGREVector( const xf_All *All, enum xfe_BasisType Basis, 
		 int Order, int StateRank, enum xfe_Bool ParallelFlag,
		 xf_Vector **pV);
/*
PURPOSE:

  Creates a new glob real elem vector sized according to All->Mesh,
  with an appropriately-dimensioned array of data allocated.  The
  vector is not placed in any DataSet structure.

INPUTS:

  All          : All structure with Mesh to link to
  Basis        : interpolation basis
  Order        : interpolation order (< 0 for not interpolated)
  StateRank    : # interpolated unknowns
  ParallelFlag : if True, vector will be prepared for parallel communication
  
OUTPUTS: 

  pV : pointer to vector just allocated
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_FindPrimalState
extern int 
xf_FindPrimalState( xf_DataSet *DataSet, int TimeIndex, xf_Data **pState,
		    char *AltName);
/*
PURPOSE:

  Looks for a primal state vector with matching TimeIndex.  Retruns the
  data pointer to the data node in the linked list, if a suitable
  state is found.

INPUTS:

  DataSet   : DataSet to look through
  TimeIndex : TimeIndex to match
  AltName   : alternate name for state (optional; can be NULL)
  
OUTPUTS: 

  pState: pointer to primal state if found; note, if multiple suitable
          states are present, only the first one is returned.
  
RETURN:

  xf_NOT_FOUND if could not find the primal state
  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_CompatibleVectors
extern enum xfe_Bool
xf_CompatibleVectors( const xf_Vector *A, const xf_Vector *B);
/*
PURPOSE:

  Checks if vectors A and B are compatible based on size of arrays.

INPUTS:

  A, B : the two vectors
  
OUTPUTS:  None
  
RETURN:

  True if vectors compatible, False otherwise
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindVector
extern int 
xf_FindVector( xf_All *All, const char Title[], enum xfe_LinkageType Linkage, 
	       int StateRank, char **StateName, int TimeIndex,
	       int MGIndex, enum xfe_BasisType *Basis, int *Order, 
	       int *nComp, int **vOrder, int *rvec,
	       enum xfe_SizeType Size, enum xfe_Bool ParallelFlag,  
	       enum xfe_Bool AddToDataSet, xf_Data **pBData, xf_Vector **pB,
	       enum xfe_Bool *Found);
/*
PURPOSE:

  Looks for a vector B in DataSet that matches the input
  characteristics: Linkage, StateRank, TimeIndex, MGIndex, etc.  If
  such a B is not found, it is created, and added to the DataSet if
  AddToDataSet is True.

INPUTS:

  All : All structure 
  Title: Title of B to search for
  Linkage, StateRank, TimeIndex, MGIndex, Size: matching characteristics
  StateName : array of strings containing names of state components (optional)
  Basis, Order : if both not NULL, these specify the desired inteprolation
                 type for each array.
  nComp, vOrder : number of components (e.g. elements) and variable order vector
                  if doing variable order interpolation
  rvec : vector of ranks for each array, one for each array.
         Optional, and only used if not interpolated.  If NULL passed
         in and not interpolated, StateRank is used as the rank for
         all arrays.
  ParallelFlag: if true, means that the created vector (*pB) will be prepped for
                parallel comm, assuming that Mesh is parallelized
  AddToDataSet: if true, B is added to DataSet if not found
  
OUTPUTS: 

  pBData: pointer to the xf_Data containing B (optional); not returned
         if pB does not exist in DataSet and AddToDataSet == False
  pB: pointer to vector B
  Found : True if matrix was found; false if not
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindSimilarVector
extern int 
xf_FindSimilarVector( xf_All *All, const xf_Vector *A, const char Title[],
		      enum xfe_Bool ParallelFlag, enum xfe_Bool AddToDataSet, 
		      xf_Data **pBData, xf_Vector **pB, enum xfe_Bool *Found);
/*
PURPOSE:

  Looks for a vector B in All->DataSet that is similar to vector A.
  Similar means same Linkage, StateRank, TimeIndex, MGIndex, and
  nArray; and also compatible data ranks in GenArray.  If such a B is
  not found, it is created, and added to the DataSet if AddToDataSet
  is True.

INPUTS:

  All : All structure 
  A: Model vector to try to match (as described above)
  Title: Title of B to search for
  ParallelFlag: if true, means that the created vector (*pB) will be prepped for
                parallel comm
  AddToDataSet: if true, B is added to DataSet if not found
  Found : True if vector was found; false if not (in which case it was created)
  
OUTPUTS: 

  pBData: pointer to the xf_Data containing B (optional); not returned
         if pB does not exist in DataSet and AddToDataSet == False
  pB: pointer to vector B
  
RETURN:

  xf_INPUT_ERROR if A exists in DataSet under Title
  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindSimilarVectors
extern int 
xf_FindSimilarVectors( xf_All *All, const xf_Vector *A, const char Title[],
		       int nVector, enum xfe_Bool ParallelFlag, enum xfe_Bool AddToDataSet, 
		       xf_Vector ***pB);
/*
PURPOSE:

  Wrapper for FindSimilarVector; seeks nVector vectors.

INPUTS:

  All : All structure 
  A: Model vector to try to match (as described above)
  Title: Title of B to search for
  nVector : number of vectors to create
  ParallelFlag: if true, means that the created vector (*pB) will be prepped for
                parallel comm
  AddToDataSet: if true, B is added to DataSet if not found
  OrderIncrement: the order sought for B will be A->Order + (*OrderIncrement)
  DesiredOrder: if not NULL, OrderIncrement is not used, and the order sought for 
                B will instead be DesiredOrder
  
OUTPUTS: 

  pBData: pointer to the xf_Data containing B (optional); not returned
         if pB does not exist in DataSet and AddToDataSet == False
  pB: pointer to vector B
  
RETURN:

  xf_INPUT_ERROR if A exists in DataSet under Title
  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_FindSimilarVectorHO
extern int 
xf_FindSimilarVectorHO( xf_All *All, xf_Vector *A, const char Title[],
			enum xfe_Bool ParallelFlag, enum xfe_Bool AddToDataSet, 
			int *OrderIncrement, int *DesiredOrder, xf_Data **pBData, 
			xf_Vector **pB);
/*
PURPOSE:

  Just like FindSimilarVector, but seeks a vector with Order increased
  by OrderIncrement (or set to DesiredOrder).  The vector must be
  interpolated.

INPUTS:

  All : All structure 
  A: Model vector to try to match (as described above)
  Title: Title of B to search for
  ParallelFlag: if true, means that the created vector (*pB) will be prepped for
                parallel comm
  AddToDataSet: if true, B is added to DataSet if not found
  OrderIncrement: the order sought for B will be A->Order + (*OrderIncrement)
  DesiredOrder: if not NULL, OrderIncrement is not used, and the order sought for 
                B will instead be DesiredOrder
  
OUTPUTS: 

  pBData: pointer to the xf_Data containing B (optional); not returned
         if pB does not exist in DataSet and AddToDataSet == False
  pB: pointer to vector B
  
RETURN:

  xf_INPUT_ERROR if A exists in DataSet under Title
  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_FindVectorSet
extern int 
xf_FindVectorSet( xf_All *All, int nVector, const char Title[],
		  enum xfe_LinkageType Linkage, int StateRank, 
		  char **StateName, int TimeIndex, int MGIndex, 
		  enum xfe_BasisType *Basis, int *Order, int *rvec,
		  enum xfe_SizeType Size, enum xfe_Bool ParallelFlag,
		  enum xfe_Bool AddToDataSet, xf_Data **pVSData, 
		  xf_VectorSet **pVS, enum xfe_Bool *Found);
/*
PURPOSE:

  Looks for a vectorset B in DataSet that matches the input
  characteristics: Linkage, StateRank, TimeIndex, MGIndex, etc.  If
  such a B is not found, it is created, and added to the DataSet if
  AddToDataSet is True.

INPUTS:

  All : All structure 
  nVector : desired number of vectors in set
  Title: Title of B to search for
  Linkage, StateRank, TimeIndex, MGIndex, Size: matching characteristics
  StateName : array of strings containing names of state components (optional)
  Basis, Order : if both not NULL, these specify the desired inteprolation
                 type for each array.
  rvec : vector of ranks for each array.  Optional, and only used if not 
         interpolated.  If NULL passed in and not interpolated, StateRank
	 is used as the rank for all arrays.
  ParallelFlag: if true, means that the created vector (*pB) will be prepped for
                parallel comm, assuming that Mesh is parallelized
  AddToDataSet: if true, B is added to DataSet if not found
  
OUTPUTS: 

  pBData: pointer to the xf_Data containing B (optional); not returned
         if pB does not exist in DataSet and AddToDataSet == False
  pB: pointer to vectorset B
  Found : True if matrix was found; false if not

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_FindSimilarVectorSet
extern int 
xf_FindSimilarVectorSet( xf_All *All, xf_Vector *A, int nVector, 
			 const char Title[], enum xfe_Bool ParallelFlag, 
			 enum xfe_Bool AddToDataSet, xf_Data **pVSData, 
			 xf_VectorSet **pVS);
/*
PURPOSE:

  Looks for a VectorSet VS in All->DataSet, whose vectors are similar
  to A.  Similar means same Linkage, StateRank, TimeIndex, MGIndex,
  and nArray; and also compatible data ranks in GenArray.  If such a VS
  is not found, it is created, and added to the DataSet if
  AddToDataSet is True.

INPUTS:

  All : All structure 
  A: Model vector to try to match (as described above)
  Title: Title of VS to search for
  ParallelFlag: if true, means that each vector in the the created vector 
                set (*pVS) will be prepped for parallel comm
  AddToDataSet: if true, VS is added to DataSet if not found
  OrderIncrement: the order sought for VS will be A->Order + OrderIncrement
  
OUTPUTS: 

  pVSData: pointer to the xf_Data containing VS (optional); not returned
         if pB does not exist in DataSet and AddToDataSet == False
  pVS: pointer to vector VS
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_FindSimilarVectorSetHO
extern int 
xf_FindSimilarVectorSetHO( xf_All *All, xf_Vector *A, int nVector, 
			   const char Title[], enum xfe_Bool ParallelFlag, 
			   enum xfe_Bool AddToDataSet, int OrderIncrement,
			   xf_Data **pVSData, xf_VectorSet **pVS);
/*
PURPOSE:

  Just like FindSimilarVectorSet but looks for a VectorSet with order
  increased by OrderIncrement.  The vectors must be interpolated.

INPUTS:

  All : All structure 
  A: Model vector to try to match (as described above)
  Title: Title of VS to search for
  ParallelFlag: if true, means that each vector in the the created vector 
                set (*pVS) will be prepped for parallel comm
  AddToDataSet: if true, VS is added to DataSet if not found
  
OUTPUTS: 

  pVSData: pointer to the xf_Data containing VS (optional); not returned
         if pB does not exist in DataSet and AddToDataSet == False
  pVS: pointer to vector VS
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindMatrix
extern int 
xf_FindMatrix( xf_DataSet *DataSet, const char Title[],
	       enum xfe_LinkageType Linkage, int LinkageIndex, 
	       int Order1, int Order2, enum xfe_BasisType Basis1,
	       enum xfe_BasisType Basis2, enum xfe_SizeType Size, 
	       int n, int r, enum xfe_Bool AddToDataSet, 
	       xf_Data **pMData, xf_Matrix **pM, enum xfe_Bool *Found);
/*
PURPOSE:

  Looks for a matrix M in DataSet that matches the input
  characteristics: Linkage, Order, Basis, size, etc.  If such an M is
  not found, it is created, and added to the DataSet if AddToDataSet
  is True.

INPUTS:

  DataSet : Data set to search
  Title: Title of M to search for
  Linkage, LinkageIndex, Order1, Order2, Basis1, Basis2,
    Size: matching matrix characteristics
  n, r : desired matrix dimensions
  AddToDataSet: if true, M is added to DataSet if not found
  
OUTPUTS: 

  pMData: pointer to the xf_Data containing M (optional); not returned
         if pM does not exist in DataSet and AddToDataSet == False
  pM: pointer to matrix M
  Found : True if matrix was found; false if not
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindOrCreatePrimalState
extern int 
xf_FindOrCreatePrimalState( xf_All *All, enum xfe_Bool RestartFlag, xf_ICs *ICsOrig,
			    xf_Vector **pU);

/*
PURPOSE:

  Looks for a Primal State vector, U, in All when RestartFlag is True.
  If found, projects U to appropriate Order, as dictated by the
  parameters in All.  If U is not found, it is created and initialized
  to the IC in All->EqnSet.

INPUTS:

  All       : All structure
  RestartFlag : If True, All will be searched for a suitable primal State.
  ICsOrig : original initial conditions, if want to scale the State
  
OUTPUTS: 

  U : pointer to primal state vector
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindAdjointVectors
extern int 
xf_FindAdjointVectors( xf_All *All, xf_Vector *U, const char *OutputList, 
		       enum xfe_Bool UseAllOutputs, enum xfe_Bool ZeroFlag, 
		       int *nPsi, xf_Vector ***pPsi, enum xfe_Bool *Found);
/*
PURPOSE:

  Looks for a set of adjoint vectors corresponding to adjoints listed
  sequentially in the string OutputList.  The adjoints are allocated
  to match the size of vector U (e.g. a state vector).  The number of
  adjoints is returned in (*nPsi), and (*pPsi) is allocated to be a
  vector of pointers to the adjoint vectors.

INPUTS:

  All : all file
  U : vector on which to base the size of the adjoint (e.g. a state vector)
  OutputList : space-separated list of output names for desired adjoints
  UseAllOutputs : if True, OutputList is ignored and all existing outputs
                  are asigned an adjoint
  ZeroFlag : if True, adjoints are zeroed out
  
OUTPUTS: 

  nPsi : number of adjoint vectors
  pPsi : vector of allocated pointers to the adjoint vectors
  Found : True if *all* adjoint vectors were found; false if one or
          more was not found (in which case it was created)
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Definition: xf_AllocateJacobianMatrix
extern int 
xf_AllocateJacobianMatrixCoarse( xf_Mesh *Mesh, int CoarseOrder, 
				 xf_JacobianMatrix *R_U);
/*
PURPOSE:

  Allocates a Jacobian Matrix, R_U->R_Uc, with the same Basis
  as used in R_U, but with Order set to CoarseOrder.

INPUTS:

  Mesh: Mesh structure 
  CoarseOrder : order to use for R_U->R_Uc
  R_U: Parent Jacobian matrix

OUTPUTS: 

  None : R_U->R_Uc is allocated
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_FindJacobianMatrix
extern int 
xf_FindJacobianMatrix( xf_Mesh *Mesh, xf_DataSet *DataSet, xf_Vector *U, 
		       const char *TitleIn, enum xfe_Bool AllocValue, 
		       xf_Data **pR_UData, xf_JacobianMatrix **pR_U, 
		       enum xfe_Bool *Found);
/*
PURPOSE:

  Looks for the Jacobian Matrix, R_U in DataSet that is suitable for
  the given Mesh and vector U.  The TimeIndex and MGIndex of R_U and U
  must match.

  If the data size of a suitable R_U is not adequate, R_U is
  reallocated.

  The Jacobian is always added to the DataSet, with Title "R_U"

INPUTS:

  Mesh: Mesh structure 
  DataSet: Data Set structure
  U: vector upon which R_U is sized
  TitleIn : optional title of Jacobian matrix.  If not specified, a
            default title will be used.
  AllocValue : if False, Value data will not be allocated.  This is
               useful if want to use R_U just for neighbor info.
  
OUTPUTS: 

  pR_UData: pointer to the xf_Data containing R_U (optional);
  pR_U: pointer to R_U
  Found : True if Jacobian matrix was found in the data
  
RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_SetZeroJacobian
extern int 
xf_SetZeroJacobian( xf_Mesh *Mesh, xf_JacobianMatrix *R_U);
/*
PURPOSE:

  Sets R_U = 0

INPUTS:

  Mesh : mesh structure
  R_U : jacobian matrix

OUTPUTS: 

  R_U is set to zero
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_GetVectorDOF
extern int 
xf_GetVectorDOF(xf_Vector *V, int *ndof);
/*
PURPOSE:

  Computes the number of degrees of freedom for a vector; does not
  include StateRank in the count.

INPUTS:

  V : vector of interest

OUTPUTS: 

  ndof : number of degrees of freedom
  
RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_DataSetInfo
extern int 
xf_DataSetInfo(xf_DataSet *DataSet);
/*
PURPOSE:

  Prints out information on a DataSet

INPUTS:

  DataSet : dataset for which to print out information

OUTPUTS: 

  None, only stdout information
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DataSetDeleteNonEssential
extern int 
xf_DataSetDeleteNonEssential(xf_DataSet *DataSet);
/*
PURPOSE:

  Deletes non-essential data from DataSet.  Only vectors or vector
  sets with a solver role are considered essential.

INPUTS:

  DataSet : dataset from which to delete data

OUTPUTS: 

  None
  
RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_WriteDataSetBinary
extern int 
xf_WriteDataSetBinary( xf_Mesh *Mesh, xf_DataSet *DataSet, FILE *fidin,
		       const char *fname);
/*
PURPOSE:

  Writes All->DataSet (or AltDataSet if passed in as not NULL) to file
  fname.  All may be parallelized, in which case All->DataSet (or
  AltDataSet) is temporarily unparallelized and written by the root
  processor.

INPUTS:

  Mesh  : only required if writing in parallel
  DataSet : dataset to write out
  fidin : file pointer to use
  fname : if not null, fidin will be ignored and file fname will 
          be opened and used for writing instaed

OUTPUTS: 

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DumpVectorBinary
extern int 
xf_DumpVectorBinary( xf_Mesh *Mesh, const char *Name,
		     xf_Vector *V, const char *fname);
/*
PURPOSE:

  Writes vector V in the .data file fname.

INPUTS:

  Mesh  : only required if writing in parallel
  Name  : Name of Vector within the data set that is written
  V     : vector to write out
  fname : file name to write

OUTPUTS: None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: Yu_DumpMultiVectorBinary
extern int 
Yu_DumpMultiVectorBinary( xf_Mesh *Mesh, const char *Name1,
		     xf_Vector *V1, const char *Name2,
                     xf_Vector *V2, const char *fname);
/*
 * Modified version 
 */

/******************************************************************/
//   FUNCTION Prototype: Yu_DumpMulti3VectorBinary
extern int 
Yu_DumpMulti3VectorBinary( xf_Mesh *Mesh, const char *Name1,
		     xf_Vector *V1, const char *Name2,
                     xf_Vector *V2, const char *Name3, 
                     xf_Vector *V3, const char *fname);
/*
 * Modified version 
 */

/******************************************************************/
//   FUNCTION Prototype: xf_ReadDataSetBinary
extern int 
xf_ReadDataSetBinary( xf_Mesh *Mesh, FILE *fidin, const char *fname,
		      xf_DataSet *DataSet);
/*
PURPOSE:

  Reads DataSet in from a file, and parallelizes it if running in
  parallel.  This function can read very large data files such
  as large vector sets as it parallelizes while reading.

INPUTS:

  Mesh  : only required if reading in parallel
  fidin : file pointer to use
  fname : if not null, fidin will be ignored and file fname will 
          be opened and used for reading instaed

OUTPUTS: 

  DataSet : dataset that was read in, valid on all procs (parallelized)

RETURN:

  Error Code; specifically xf_NOT_FOUND if file could not be opened
  (most likely because it does not exist).  This is a silent error (no
  messages to the screen, just the error code) so that this function
  can be used to check if a file exists.

*/


/******************************************************************/
//   FUNCTION Prototype: xf_BuildVOrder
extern int 
xf_BuildVOrder(xf_All *All, xf_Vector *U, xf_Vector **pV);
/*
PURPOSE:

  Builds an integer vector, V, that contains the elemental orders of
  the interpolated vector U.

INPUTS:

  All   : all structure
  U     : input vector (must be interpolated)

OUTPUTS: 

  V     : output integer vector of interpolation errors

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Definition: xf_FindDataByPointer
extern int 
xf_FindDataByPointer( xf_DataSet *DataSet, void *TargetData, xf_Data **Data);
/*
 PURPOSE:
 
 Loops through the DataSet and finds the data node that contains TargetData.
 
 INPUTS:
 
 DataSet   : dataset structure
 TargetData: pointer to data (xf_Vector, xf_Matrix, ect.)
 
 OUTPUTS: 
 
 Data      : if TargetData is part of the DataSet, then Data points to 
             corresponding data node.
 
 RETURN:
 
 Error Code
 */


#endif // end ifndef _xf_Data_h
