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
  FILE:  xf_Unit.c

  This file contains unit-test-specific functions.

*/
#include "xf.h"
#include <math.h>
#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_Mesh.h"
#include "xf_Param.h"
#include "xf_EqnSet.h"
#include "xf_Data.h"
#include "xf_All.h"
#include "xf_Memory.h"
#include "xf_EqnSetHook.h"
#include "xf_Solver.h"
#include "xf_Basis.h"
#include "xf_MeshMotion.h"
#include "xf_MeshMotionIO.h"
#include "xf_UnitStruct.h"

/******************************************************************/
//   FUNCTION Definition: AssertIntVectorEqual
int
AssertIntVectorEqual(int *A, int *B, int n)
{
  int i;

  for (i=0; i<n; i++)
    if (A[i] != B[i]) return (i+1);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: AssertRealVectorWithin
int
AssertRealVectorWithin(real *A, real *B, int n, real c)
{
  int i;

  for (i=0; i<n; i++)
    if (fabs(A[i]-B[i]) > c) return (i+1);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xfu_PingPerturb
extern real
xfu_PingPerturb(real eps, int ieps){
  int k;
  real ep;

  ep = eps;
  for (k=0; k<ieps; k++) ep *= 0.5;
  
  return ep;
}


/******************************************************************/
//   FUNCTION Definition: xfu_PingCheckRate
extern int
xfu_PingCheckRate(real *veps, real tol, real meps){
  real rate;

  //xf_printf("veps = [%.6E, %.6E]\n", veps[0], veps[1]);

  if (veps[0] <= meps) return xf_OK; // machine precision check
  rate = log(veps[0]/veps[1])/log(2.0);
  //xf_printf("rate = %.4f\n", rate);
  if (rate < (2.0-tol)){
    xf_printf("Ping failure: observed rate of %.5f is less than the expected rate of %.5f\n",
	      rate, 2.0-tol);
    return xf_Error(xf_PING_FAILED);
  }
  else
    return xf_OK;
}



// Utility functions

/******************************************************************/
//   FUNCTION Definition: xf_UnitirVector
static int
xf_UnitirVector (int nArray, int *n, int *r, int **vr, enum xfe_SizeType Size, 
		 xf_Vector **pV)
{
  int ierr, i, j;

  ierr = xf_Error(xf_CreateVector(pV));
  if (ierr != xf_OK) return ierr;

  (*pV)->nArraySelf = (*pV)->nArray = nArray;

  // allocate GenArray structures
  ierr = xf_Error(xf_Alloc((void **) &(*pV)->GenArray, (*pV)->nArray, sizeof(xf_GenArray)));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<nArray; i++){
    (*pV)->GenArray[i].Size = Size;
    (*pV)->GenArray[i].n    = n[i];
    (*pV)->GenArray[i].r    = r[i];

    (*pV)->GenArray[i].vr   = NULL;
    if (vr != NULL){ // variable order support
      ierr = xf_Error(xf_Alloc((void **) &(*pV)->GenArray[i].vr, n[i], sizeof(int)));
      if (ierr != xf_OK) return ierr;
      for (j=0; j<n[i]; j++) (*pV)->GenArray[i].vr[j]   = vr[i][j];
    }

    ierr = xf_Error(xf_AllocGenArray((*pV)->GenArray+i));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitiVector
int
xf_UnitiVector (int nArray, int *n, int *r, int **vr, xf_Vector **pV)
{  
  return xf_Error(xf_UnitirVector(nArray, n, r, vr, xfe_SizeInt, pV));
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitrVector
int
xf_UnitrVector (int nArray, int *n, int *r, int **vr, xf_Vector **pV)
{  
  return xf_Error(xf_UnitirVector(nArray, n, r, vr, xfe_SizeReal, pV));
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitrVectorInterp
int
xf_UnitrVectorInterp(int nArray, int *n, int sr, enum xfe_BasisType *Basis,
		     int *Order, int **vOrder, xf_Vector **pV)
{
  int ierr, i, j, nn;
  int **vr = NULL, *r = NULL;
  enum xfe_Bool VariableOrder;

  // are we doing variable orders?
  VariableOrder = (vOrder != NULL);

  // allocate r
  ierr = xf_Error(xf_Alloc( (void **) &r, nArray, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // if so, allocate and fill vr
  vr = NULL;
  if (VariableOrder){
    ierr = xf_Error(xf_VAlloc2( (void ***) &vr, nArray, n, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nArray; i++){
      r[i] = -1;
      for (j=0; j<n[i]; j++){
	ierr = xf_Error(xf_Order2nNode(Basis[i], vOrder[i][j], &nn));
	if (ierr != xf_OK) return ierr;
	vr[i][j] = nn*sr;	
	r[i] = max(r[i], vr[i][j]);
      } // j
    } // i
  }
  else{ // just use Order to determine r
    for (i=0; i<nArray; i++){
      ierr = xf_Error(xf_Order2nNode(Basis[i], Order[i], &nn));
      if (ierr != xf_OK) return ierr;
      r[i] = nn*sr;
    } // i
  }

  // create the vector
  ierr = xf_Error(xf_UnitrVector(nArray, n, r, vr, pV));
  if (ierr != xf_OK) return ierr;

  // add a linkage (some functions require this)
  (*pV)->Linkage = xfe_LinkageGlobElem;

  // Set StateRank
  (*pV)->StateRank = sr;

  // allocate and fill Basis
  ierr = xf_Error(xf_Alloc((void **) &(*pV)->Basis, nArray, sizeof(enum xfe_BasisType)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nArray; i++) (*pV)->Basis[i] = Basis[i];

  // allocate and fill Order
  ierr = xf_Error(xf_Alloc((void **) &(*pV)->Order, nArray, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nArray; i++) (*pV)->Order[i] = Order[i];

  // allocate and fill nComp, vOrder if VariableOrder
  if (VariableOrder){
    ierr = xf_Error(xf_Alloc((void **) &(*pV)->nComp, nArray, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nArray; i++) (*pV)->nComp[i] = n[i];
    ierr = xf_Error(xf_VAlloc2( (void ***) &(*pV)->vOrder, nArray, n, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<nArray; i++)
      for (j=0; j<n[i]; j++)
	(*pV)->vOrder[i][j] = vOrder[i][j];
  }

  xf_Release2( (void **) vr);
  xf_Release ( (void  *) r);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitrVectorSet
int
xf_UnitrVectorSet (int nVector, int nArray, int *n, int *r, 
		   xf_VectorSet **pVS)
{
  int ierr, i;
  xf_Vector *V;

  // create a vectorset
  ierr = xf_Error(xf_CreateVectorSet(nVector, pVS));
  if (ierr != xf_OK) return ierr;
  
  // allocate vectors in VS
  for (i=0; i<nVector; i++){
    ierr = xf_Error(xf_UnitrVector(nArray, n, r, NULL, &V));
    if (ierr != xf_OK) return ierr;
    (*pVS)->Vector[i] = *V;
    xf_InitVector(V);
    ierr = xf_Error(xf_DestroyVector(V, xfe_True));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/* A .gri file for a Two-Q1 triangle mesh */
static char *TwoQ1TriangleGri[] =
{
  "4 2 2",
  "0 0",
  "1 0",
  "0 1",
  "1 1",
  "1",
  "4 2 All",
  "1 2",
  "2 4",
  "4 3",
  "3 1",
  "2 1 TriLagrange",
  "1 2 3",
  "2 4 3",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitTwoQ1TriangleMesh
int
xf_UnitTwoQ1TriangleMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag)
{
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, TwoQ1TriangleGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitTwoQ1TriangleAll
int
xf_UnitTwoQ1TriangleAll(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&(*pAll)->Mesh, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/* A .gri file for a Two-Q2 triangle mesh */
static char *TwoQ2TriangleGri[] =
{
  "9 2 2",
  " 0  0",
  ".5  0",
  " 1  0",
  " 0 .5",
  ".5 .5",
  " 1 .5",
  " 0  1",
  ".5  1",
  " 1  1",
  "1",
  "4 2 All",
  "1 3",
  "3 9",
  "9 7",
  "7 1",
  "2 2 TriLagrange",
  "1 2 3 4 5 7",
  "3 6 9 5 8 7",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitTwoQ2TriangleMesh
int
xf_UnitTwoQ2TriangleMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag)
{
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, TwoQ2TriangleGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitTwoQ2TriangleAll
int
xf_UnitTwoQ2TriangleAll(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_UnitTwoQ2TriangleMesh(&(*pAll)->Mesh, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/* A .gri file for a Q1 tet mesh */
static char *OneQ1TetGri[] =
{
  "5 1 3",
  " 0 0 0", // extra node so that local/global node #s do not coincide
  " 0 0 0",
  " 1 0 0",
  " 0 1 0",
  " 0 0 1",
  "1",
  "4 3 All",
  "3 4 5",
  "4 2 5",
  "2 3 5",
  "2 4 3",
  "1 1 TetLagrange",
  "2 3 4 5",
  "\0",
};

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ1TetMesh
int
xf_UnitQ1TetMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag)
{
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, OneQ1TetGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ1TetAll
int
xf_UnitQ1TetAll(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_UnitQ1TetMesh(&(*pAll)->Mesh, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/* A .gri file for a Q2 tet mesh */
static char *OneQ2TetGri[] =
{
  "11 1 3",
  "0.0 0.0 0.0", // extra node so that local/global node #s do not coincide
  "0.0 0.0 0.0",
  "0.5 0.0 0.0",
  "1.0 0.0 0.0",
  "0.0 0.5 0.0",
  "0.5 0.5 0.0",
  "0.0 1.0 0.0",
  "0.0 0.0 0.5",
  "0.5 0.0 0.5",
  "0.0 0.5 0.5",
  "0.0 0.0 1.0",
  "1",
  "4 3 All",
  "4 7 11",
  "7 2 11",
  "2 4 11",
  "2 7 4",
  "1 2 TetLagrange",
  "2 3 4 5 6 7 8 9 10 11",
  "\0",
};

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ2TetMesh
int
xf_UnitQ2TetMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag)
{
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, OneQ2TetGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ2TetAll
int
xf_UnitQ2TetAll(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_UnitQ2TetMesh(&(*pAll)->Mesh, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/* A .gri file for a Q1 Quad mesh */
static char *OneQ1QuadGri[] =
{
  "5 1 2",
  "0.0 0.0", // extra node so that local/global node #s do not coincide
  "0.0 0.0",
  "1.0 0.0",
  "0.0 1.0",
  "1.0 1.0",
  "1",
  "4 2 All",
  "2 3",
  "3 5",
  "5 4",
  "4 2",
  "1 1 QuadLagrange",
  "2 3 4 5",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitQ1QuadMesh
int
xf_UnitQ1QuadMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag)
{
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, OneQ1QuadGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ1QuadAll
int
xf_UnitQ1QuadAll(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_UnitQ1QuadMesh(&(*pAll)->Mesh, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/* A .gri file for a Q2 Quad mesh */
static char *OneQ2QuadGri[] =
{
  "10 1 2",
  "0.0 0.0", // extra node so that local/global node #s do not coincide
  "0.0 0.0",
  "0.5 0.0",
  "1.0 0.0",
  "0.0 0.5",
  "0.5 0.5",
  "1.0 0.5",
  "0.0 1.0",
  "0.5 1.0",
  "1.0 1.0",
  "1",
  "4 2 All",
  "2 4",
  "4 10",
  "10 8",
  "8 2",
  "1 2 QuadLagrange",
  "2 3 4 5 6 7 8 9 10",
  "\0",
};

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ2QuadMesh
int
xf_UnitQ2QuadMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag)
{
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, OneQ2QuadGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ2QuadAll
int
xf_UnitQ2QuadAll(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_UnitQ2QuadMesh(&(*pAll)->Mesh, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/* A .gri file for a Q1 Hex mesh */
static char *OneQ1HexGri[] =
{
  "9 1 3",
  "0 0 0", // extra node so that local/global node #s do not coincide
  "0 0 0",
  "1 0 0",
  "0 1 0",
  "1 1 0",
  "0 0 1",
  "1 0 1",
  "0 1 1",
  "1 1 1",
  "1",
  "6 4 All",
  "2 3 7 6",
  "3 5 9 7",
  "5 4 8 9",
  "4 2 6 8",
  "2 4 5 3",
  "6 7 9 8",
  "1 1 HexLagrange",
  "2 3 4 5 6 7 8 9",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitQ1HexMesh
int
xf_UnitQ1HexMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag)
{
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, OneQ1HexGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ1HexAll
int
xf_UnitQ1HexAll(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_UnitQ1HexMesh(&(*pAll)->Mesh, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/* A .gri file for a Q2 Hex mesh */
static char *OneQ2HexGri[] =
{
  "28 1 3",
  "0 0 0", // extra node so that local/global node #s do not coincide
  "0 0 0",  "1 0 0",  "2 0 0",
  "0 1 0",  "1 1 0",  "2 1 0",
  "0 2 0",  "1 2 0",  "2 2 0",
  "0 0 1",  "1 0 1",  "2 0 1",
  "0 1 1",  "1 1 1",  "2 1 1",
  "0 2 1",  "1 2 1",  "2 2 1",
  "0 0 2",  "1 0 2",  "2 0 2",
  "0 1 2",  "1 1 2",  "2 1 2",
  "0 2 2",  "1 2 2",  "2 2 2",
  "1",
  "6 4 All",
  "2 4 22 20",
  "4 10 28 22",
  "10 8 26 28",
  "8 2 20 26",
  "2 8 10 4",
  "20 22 28 26",
  "1 2 HexLagrange",
  "2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitQ2HexMesh
int
xf_UnitQ2HexMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag)
{
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, OneQ2HexGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitQ2HexAll
int
xf_UnitQ2HexAll(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_UnitQ2HexMesh(&(*pAll)->Mesh, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/* A .gri file for a Q1 triangle mesh of a 2D box (8 elems) */
static char *BoxQ1TriangleGri[] =
{
  "9 8 2",
  "0 0", "1 0", "2 0",
  "0 1", "1 1", "2 1",
  "0 2", "1 2", "2 2",
  "4",
  "2 2 Left",
  "1 4","4 7",
  "2 2 Right",
  "3 6","6 9",
  "2 2 Bottom",
  "1 2","2 3",
  "2 2 Top",
  "9 8","8 7",
  "8 1 TriLagrange",
  "1 2 4", "2 5 4", "2 3 5", "3 6 5",
  "4 5 7", "5 8 7", "5 6 8", "6 9 8",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1TriangleMesh
int
xf_UnitBoxQ1TriangleMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag){
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1TriangleGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/* A Scalar .eqn file for a 2D Box mesh */
static char *Box2DScalarEqn[] =
{
  "EqnSetLibrary = libScalar.so",
  // Parameters
  "STARTBLOCK PARAM",
  "Viscosity = 0.1",
  "ENDBLOCK",
  // Residual Terms
  "STARTBLOCK RESTERM",
  "nResTerm = 1",
  "TermType = Diffusion",
  //"TermType = Convection",
  //"VelocityFcn = Constant",
  //"VelocityData = .1 0.0", 
  "ENDBLOCK",
  // Initial Conditions
  "STARTBLOCK IC",
  "ICType = FullState",
  "Data = 0.5",
  "ENDBLOCK",
  // Boundary Conditions
  "STARTBLOCK BC",
  "nBC = 4",
  "BFGTitle = Left",  "BCType = FullState", "Data = 1.0",
  "BFGTitle = Right", "BCType = FullState", "Data = 0.0",
  "BFGTitle = Bottom","BCType = HeatFlux",  "Data = 0.0",
  "BFGTitle = Top",   "BCType = HeatFlux",  "Data = 0.0",
  "ENDBLOCK",
  // Outputs
  "STARTBLOCK OUTPUT",
  "nOutput = 6",
  "OutputName = HeatFlow",
  "OutputType = BoundaryIntegral",
  "UsesFlux = True",
  "FluxComponentNames = Scalar",
  "FluxComponentWeights = 1.0",
  "BFGTitles = Right",
  "OutputName = HeatFlowMoment",
  "OutputType = BoundaryIntegral",
  "UsesFlux = True",
  "FluxComponentNames = Scalar",
  "FluxComponentWeights = 1.0",
  "FluxComponentMoments = 1",
  "BFGTitles = Right",
  "OutputName = HeatFlowUnsteady",
  "OutputType = BoundaryIntegral",
  "TimeNorm = Final",
  "UsesFlux = True",
  "FluxComponentNames = Scalar",
  "FluxComponentWeights = 1.0",
  "BFGTitles = Right",
  "OutputName = HeatFlowUnsteadyInt",
  "OutputType = BoundaryIntegral",
  "TimeNorm = Integral",
  "StartTime = 0.4",
  "EndTime = 0.8",
  "UsesFlux = True",
  "FluxComponentNames = Scalar",
  "FluxComponentWeights = 1.0",
  "BFGTitles = Right",
  "OutputName = HeatFlowUnsteadyPoint",
  "OutputType = BoundaryIntegral",
  "TimeNorm = Point",
  "StartTime = 0.7",
  "UsesFlux = True",
  "FluxComponentNames = Scalar",
  "FluxComponentWeights = 1.0",
  "BFGTitles = Right",
  "OutputName = HeatFlowUnsteadySquare",
  "OutputType = BoundaryIntegral",
  "TimeNorm = SquareIntegral",
  "StartTime = 0.3",
  "EndTime = 0.7",
  "UsesFlux = True",
  "FluxComponentNames = Scalar",
  "FluxComponentWeights = 1.0",
  "BFGTitles = Right",
  "ENDBLOCK",
  "\0",
};


/* A Scalar .eqn file for a 2D Box mesh */
static char *Box2DScalarConvEqn[] =
{
  "EqnSetLibrary = libScalar.so",
  // Parameters
  "STARTBLOCK PARAM",
  "Viscosity = 0.1",
  "ENDBLOCK",
  // Residual Terms
  "STARTBLOCK RESTERM",
  "nResTerm = 1",
  "TermType = Convection",
  "VelocityFcn = Constant",
  "VelocityData = 1.0 0.5", 
  "ENDBLOCK",
  // Initial Conditions
  "STARTBLOCK IC",
  "ICType = FullState",
  "Function = Gaussian",
  "Data = 50 0.5 50 0.5",
  "ENDBLOCK",
  // Boundary Conditions
  "STARTBLOCK BC",
  "nBC = 4",
  "BFGTitle = Left",  "BCType = FullState", "Data = 0.0",
  "BFGTitle = Bottom","BCType = FullState", "Data = 0.0",
  "BFGTitle = Right", "BCType = None",
  "BFGTitle = Top",   "BCType = None", 
  "ENDBLOCK",
  // Outputs
  "STARTBLOCK OUTPUT",
  "nOutput = 0",
  "ENDBLOCK",
  "\0",
};

static char *Box2DEulerPenaltyEqn[] =
{
  "EqnSetLibrary = libCompressibleNS.so",
  //Parameters
  "STARTBLOCK PARAM",
  "SpecificHeatRatio = 1.4",
  "GasConstant = 0.4",
  "PenaltyFcnFactor = 1e0",
  "PenaltyFcnPower = 1.0",
  "ENDBLOCK",
  //Residual Terms
  "STARTBLOCK RESTERM",
  "nResTerm = 1",
  "TermType = Convection",
  "InteriorFlux = Standard",
  "FluxFunction = Roe",
  "ENDBLOCK",
  //Initial Conditions
  "STARTBLOCK IC",
  "ICType = FullState",
  "Data = 1.0 -0.8 0.6 0.6",
  "ENDBLOCK",
  //Boundary Conditions
  "STARTBLOCK BC",
  "nBC = 4",
  //Top
  "BFGTitle = Top",
  "BCType = InviscidWall",
  //Bottom
  "BFGTitle = Bottom",
  "BCType = InviscidWall",
  //Left
  "BFGTitle = Left",
  "BCType = FullState",
  "Data = 1.0 1.0 0.0 7.64285714286",
  //Right
  "BFGTitle = Right",
  "BCType = FullState",
  "Data = 1.0 1.0 0.0 7.64285714286",
  "ENDBLOCK",
  "\0",
};

/* An Euler .eqn file for a 2D Box mesh */
static char *Box2DEulerEqn[] =
{
  "EqnSetLibrary = libCompressibleNS.so",
  // Parameters
  "STARTBLOCK PARAM",
  "SpecificHeatRatio = 1.4",
  "GasConstant = 1.0",
  "Viscosity = .1",
  "RegularityScalar = Density",
  "ENDBLOCK",
  // Residual Terms
  "STARTBLOCK RESTERM",
  "nResTerm = 1",
  "TermType = Convection",
  "InteriorFlux = Standard",
  "FluxFunction = Roe",
  "ENDBLOCK",
  // Initial Conditions
  "STARTBLOCK IC",
  "ICType = FullState",
  "Data = 1.0 0.59160797830996159 0.0 2.675",
  "ENDBLOCK",
  // Boundary Conditions
  "STARTBLOCK BC",
  "nBC = 4",
  "BFGTitle = Left",  "BCType = TtPta", "Data = 1.05 1.186212638044398 0.0",
  "BFGTitle = Right", "BCType = StaticP", "Data = 1.0",
  "BFGTitle = Bottom","BCType = InviscidWall",
  "BFGTitle = Top",   "BCType = InviscidWall",
  "ENDBLOCK",
  // Outputs
  "STARTBLOCK OUTPUT",
  "nOutput = 4",
  //Output 1
  "OutputName = EntropyNorm",
  "OutputType = DomainIntegral",
  "DomainNorm = None",
  "UsesFlux = False",
  "ScalarName = Entropy",
  //Output 2
  "OutputName = BoundaryInt",
  "OutputType = BoundaryIntegral",
  "BFGTitles = Right",
  "UsesFlux = True",
  "FluxComponentNames = XMomentum YMomentum",
  "FluxComponentWeights = 1.0 0.5",
  "TimeNorm = Final",
   //Output 3
  "OutputName = BoundaryIntNoFlux",
  "OutputType = BoundaryIntegral",
  "BFGTitles = Right",
  "UsesFlux = False",
  "ScalarName = Density",
  "TimeNorm = Final",
  //Output 4
  "OutputName = Point",
  "OutputType = PointValue",
  "TimeNorm = Final",
  "egrp = 0",
  "elem = 3",
  "xref = 0.1 0.8",
  "ScalarName = Density",
  "ENDBLOCK",
  "\0",
};

/* A CNS .eqn file for a 2D Box mesh */
static char *Box2DCNSEqn[] =
{
  "EqnSetLibrary = libCompressibleNS.so",
  // Parameters
  "STARTBLOCK PARAM",
  "SpecificHeatRatio = 1.4",
  "GasConstant = 1.0",
  "Viscosity = .1",
  "RegularityScalar = Density",
  "ENDBLOCK",
  // Residual Terms
  "STARTBLOCK RESTERM",
  "nResTerm = 2",
  "TermType = Convection",
  "InteriorFlux = Standard",
  "FluxFunction = Roe",
  "TermType = Diffusion",
  "ENDBLOCK",
  // Initial Conditions
  "STARTBLOCK IC",
  "ICType = FullState",
  "Data = 1.0 0.59160797830996159 0.01 2.675",
  "ENDBLOCK",
  // Boundary Conditions
  "STARTBLOCK BC",
  "nBC = 4",
  "BFGTitle = Left",  "BCType = FullState", "Data = 1.0 0.59160797830996159 0.01 2.675",
  //"BFGTitle = Left",  "BCType = TtPta", "Data = 1.05 1.186212638044398 -0.01",
  "BFGTitle = Right", "BCType = StaticP", "Data = 1.0",
  "BFGTitle = Bottom","BCType = NoSlipHeat", "Data = 0.0",
  "BFGTitle = Top",   "BCType = NoSlipTemp", "Data = 1.0",
  /* "BFGTitle = Left",  "BCType = FullState", "Data = 1.0 0.59160797830996159 0.01 2.675", */
  /* "BFGTitle = Right",  "BCType = FullState", "Data = 1.0 0.59160797830996159 0.01 2.675", */
  /* "BFGTitle = Bottom",  "BCType = FullState", "Data = 1.0 0.59160797830996159 0.01 2.675", */
  /* "BFGTitle = Top",  "BCType = FullState", "Data = 1.0 0.59160797830996159 0.01 2.675", */

  "ENDBLOCK",
  // Outputs
  "STARTBLOCK OUTPUT",
  "nOutput = 7",
  // Output 1
  "OutputName = EntropyNorm",
  "OutputType = DomainIntegral",
  "DomainNorm = None",
  "UsesFlux = False",
  "ScalarName = Entropy",
  // Output 2
  "OutputName = PressureIntegral",
  "OutputType = BoundaryIntegral",
  "UsesFlux = False",
  "ScalarName = Pressure",
  "BFGTitles = Bottom",
  // Output 3
  "OutputName = BoundaryInt",
  "OutputType = BoundaryIntegral",
  "BFGTitles = Bottom",
  "UsesFlux = True",
  "FluxComponentNames = XMomentum YMomentum",
  "FluxComponentWeights = 1.0 0.5",
  "TimeNorm = Final",
  // Output 4
  "OutputName = BoundaryIntNoFlux",
  "OutputType = BoundaryIntegral",
  "BFGTitles = Top",
  "UsesFlux = False",
  "ScalarName = Density",
  "TimeNorm = Final",
  // Output 5
  "OutputName = Point",
  "OutputType = PointValue",
  "TimeNorm = Final",
  "egrp = 0",
  "elem = 3",
  "xref = 0.1 0.8",
  "ScalarName = Density",
  // Output 6
  "OutputName = BoundaryIntFS",
  "OutputType = BoundaryIntegral",
  "BFGTitles = Left",
  "UsesFlux = True",
  "FluxComponentNames = XMomentum YMomentum",
  "FluxComponentWeights = 0.7 -0.2",
  "TimeNorm = Final",
  // Output 7
  "OutputName = BoundaryIntTemp",
  "OutputType = BoundaryIntegral",
  "BFGTitles = Top",
  "UsesFlux = True",
  "FluxComponentNames = XMomentum YMomentum",
  "FluxComponentWeights = 1.0 0.5",
  "TimeNorm = Final",
  "ENDBLOCK",
  "\0",
};

/* A CNS_SA (RANS) .eqn file for a 2D Box mesh */
static char *Box2DCNS_SAEqn[] =
{
  "EqnSetLibrary = libCompressibleNS.so",
  // Parameters
  "STARTBLOCK PARAM",
  "SpecificHeatRatio = 1.4",
  "GasConstant = 1.0",
  "Viscosity = .1",
  "SA_NonDim = 10.0",
  "RegularityScalar = Density",
  "ENDBLOCK",
  // Residual Terms
  "STARTBLOCK RESTERM",
  "nResTerm = 3",
  "TermType = Convection",
  "InteriorFlux = Standard",
  "FluxFunction = Roe",
  "TermType = Diffusion",
  "TermType = Source",
  "SourceFcn = Turbulence",
  "TurbModel = SA",
  "ENDBLOCK",
  // Initial Conditions
  "STARTBLOCK IC",
  "ICType = FullState",
  "Header = Density XMomentum YMomentum Energy TurbNutil",
  "Data = 1.1 0.6796022962302372 0.0118625022060646 3.21 1e-2",
  "ENDBLOCK",
  // Boundary Conditions
  "STARTBLOCK BC",
  "nBC = 4",
  "BFGTitle = Left",  "BCType = TtPta", "Header = Tt pt alpha TurbNutil", 
                                        "Data = 1.14545454545455 1.42345516565328 .02 2e-2",
  "BFGTitle = Right", "BCType = StaticP", "Data = 1.2",
  "BFGTitle = Bottom","BCType = NoSlipHeat", "Header = Q TurbNutil", "Data = 0.0 0.0",
  "BFGTitle = Top",   "BCType = NoSlipTemp", "Header = T TurbNutil", "Data = 1.0 0.0",
  "ENDBLOCK",
  // Outputs
  "STARTBLOCK OUTPUT",
  "nOutput = 2",
  // Output 1
  "OutputName = EntropyNorm",
  "OutputType = DomainIntegral",
  "DomainNorm = None",
  "UsesFlux = False",
  "ScalarName = Entropy",
  // Output 2
  "OutputName = PressureIntegral",
  "OutputType = BoundaryIntegral",
  "UsesFlux = False",
  "ScalarName = Pressure",
  "BFGTitles = Bottom",
  "ENDBLOCK",
  "\0",
};



/* A 2D test motion file, plunge about 1,1 */
static char *Box2DMMPlunge1[] =
{
  "STARTBLOCK ANALYTICAL",
  "nTerm = 1",
  "Name = PlungeY",
  "MotionType = Plunge",
  "YAmplitude = 0.2",
  "Frequency = 5.4",
  "Phase = 0.0",
  "BlendType = Cubic",
  "XCenter = 1.0",
  "YCenter = 1.0",
  "Radius = 0.2",
  "Dist = 0.8",
  "ENDBLOCK",
  "\0",
};

/* A 2D test motion file, plunge about 1,0.1 */
static char *Box2DMMPlunge2[] =
{
  "STARTBLOCK ANALYTICAL",
  "nTerm = 1",
  "Name = PlungeY",
  "MotionType = Plunge",
  "YAmplitude = 0.2",
  "Frequency = 5.75",
  "Phase = 0.0",
  "BlendType = Cubic",
  "XCenter = 1.0",
  "YCenter = 0.1",
  "Radius = 0.2",
  "Dist = 0.8",
  "ENDBLOCK",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_GetMotionStrings
static char**
xf_GetMotionStrings(enum xfe_MotionType MotionType)
{
  switch (MotionType){
  case xfe_UnitMotionNone:    return NULL;           break;
  case xfe_UnitMotionPlunge1: return Box2DMMPlunge1; break;
  case xfe_UnitMotionPlunge2: return Box2DMMPlunge2; break;
  default: return NULL; break;
  }
}



/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1TriangleAll
int
xf_UnitBoxQ1TriangleAll(xf_All **pAll, enum xfe_MotionType MotionType)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1TriangleGri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DScalarEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  (*pAll)->EqnSet->Dim = (*pAll)->Mesh->Dim; // Set EqnSet dimension

  // mesh motion
  if (MotionType != xfe_UnitMotionNone){
    ierr = xf_Error(xf_CreateMeshMotion(&(*pAll)->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadMeshMotionFile(NULL, xf_GetMotionStrings(MotionType), 
					  (*pAll)->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1TriangleAllEuler
int
xf_UnitBoxQ1TriangleAllEuler(xf_All **pAll, enum xfe_UnitMotionType MotionType)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1TriangleGri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DEulerEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  // mesh motion
  if (MotionType != xfe_UnitMotionNone){
    ierr = xf_Error(xf_CreateMeshMotion(&(*pAll)->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadMeshMotionFile(NULL, xf_GetMotionStrings(MotionType), 
					  (*pAll)->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1TriangleAllCNS
int
xf_UnitBoxQ1TriangleAllCNS(xf_All **pAll, enum xfe_UnitMotionType MotionType)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1TriangleGri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DCNSEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  // mesh motion
  if (MotionType != xfe_UnitMotionNone){
    ierr = xf_Error(xf_CreateMeshMotion(&(*pAll)->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadMeshMotionFile(NULL, xf_GetMotionStrings(MotionType), 
					  (*pAll)->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1TriangleAllCNS_SA
int
xf_UnitBoxQ1TriangleAllCNS_SA(xf_All **pAll, enum xfe_UnitMotionType MotionType)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1TriangleGri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DCNS_SAEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  // mesh motion
  if (MotionType != xfe_UnitMotionNone){
    ierr = xf_Error(xf_CreateMeshMotion(&(*pAll)->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReadMeshMotionFile(NULL, xf_GetMotionStrings(MotionType), 
					  (*pAll)->Mesh->Motion));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/* A .gri file for a small deformed Q2 triangle mesh of a 2D box (2 elems) */
static char *BoxQ2TriangleGri[] =
{
  "9 2 2",
  "0 0", "1 0.05", "2 0",
  "-.04 1", "1.03 1.06", "2 1.02",
  "0 2", "1 1.95", "2 2",
  "4",
  "1 2 Left",
  "1 7",
  "1 2 Right",
  "3 9",
  "1 2 Bottom",
  "1 3",
  "1 2 Top",
  "7 9",
  "2 2 TriLagrange",
  "1 2 3 4 5 7",
  "3 6 9 5 8 7",
  "\0",
};

/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ2TriangleMesh
int
xf_UnitBoxQ2TriangleMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag){
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ2TriangleGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ2TriangleAllEuler
int
xf_UnitBoxQ2TriangleAllEuler(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ2TriangleGri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DEulerEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ2TriangleAllCNS
int
xf_UnitBoxQ2TriangleAllCNS(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ2TriangleGri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DCNSEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}




/* A .gri file for a Q1/Q2 Hex mesh */
static char *Q1Q2HexGri[] =
{
"31 2 3",
"   0.00    0.0    0.",
"   0.25    0.0    0.",
"   0.50    0.0    0.",
"   0.00    0.5    0.",
"   0.25    0.5    0.",
"   0.50    0.5    0.",
"   0.00    0.0   .75",
"   0.25    0.0   .75",
"   0.50    0.0   .75",
"   0.00    0.5   .75",
"   0.25    0.5   .75",
"   0.50    0.5   .75",
"   0.125   0.01   0.02",
"   0.0     0.26  -0.01",
"   0.125   0.25   0.01",
"   0.25    0.25   0.0",
"   0.125   0.51   0.0",
"   0.0    -0.01   0.375",
"   0.125   0.0    0.374",
"   0.25    0.0    0.375",
"   0.0     0.25   0.375",
"   0.125   0.25   0.375",
"   0.25    0.25   0.375",
"   0.0     0.5    0.375",
"   0.125   0.51   0.375",
"   0.25    0.5    0.375",
"   0.125   0.0    0.75",
"   0.0    0.25    0.75",
"   0.125  0.24   0.76",
"   0.25   0.25   0.75",
"   0.125  0.51   0.75",
"6",
"1 4 Left",
" 1 4 10 7",
"1 4 Right",
" 3 6 12 9",
"2 4 Front",
" 1 2 8 7",
" 2 3 9 8",
"2 4 Bottom",
" 1 2 5 4",
" 2 3 6 5",
"2 4 Top",
" 7 8 11 10",
" 8 9 12 11",
"2 4 Back",
" 4 5 11 10",
" 5 6 12 11",
"1 2 HexLagrange",
"1 13 2 14 15 16 4 17 5 18 19 20 21 22 23 24 25 26 7 27 8 28 29 30 10 31 11",
"1 1 HexLagrange",
"2 3 5 6 8 9 11 12",
  "\0",
};




/* A CNS .eqn file for a 3D Box mesh */
static char *Box3DCNSEqn[] =
{
  "EqnSetLibrary = libCompressibleNS.so",
  // Parameters
  "STARTBLOCK PARAM",
  "SpecificHeatRatio = 1.4",
  "GasConstant = 1.0",
  "Viscosity = .1",
  "RegularityScalar = Density",
  "ENDBLOCK",
  // Residual Terms
  "STARTBLOCK RESTERM",
  "nResTerm = 2",
  "TermType = Convection",
  "InteriorFlux = Standard",
  "FluxFunction = Roe",
  "TermType = Diffusion",
  "ENDBLOCK",
  // Initial Conditions
  "STARTBLOCK IC",
  "ICType = FullState",
  "Data = 1.0 0.59160797830996159 0.01 .02 2.675",
  "ENDBLOCK",
  // Boundary Conditions
  "STARTBLOCK BC",
  "nBC = 6",
  "BFGTitle = Left",  "BCType = TtPta", "Data = 1.05 1.186212638044398 0.01 -.02",
  "BFGTitle = Right", "BCType = StaticP", "Data = 1.0",
  "BFGTitle = Bottom","BCType = NoSlipHeat", "Data = 0.0",
  "BFGTitle = Top",   "BCType = NoSlipTemp", "Data = 1.0",
  "BFGTitle = Front", "BCType = Symmetry",
  "BFGTitle = Back",  "BCType = NoSlipHeat", "Data = 1.0",
  "ENDBLOCK",
  // Outputs
  "STARTBLOCK OUTPUT",
  "nOutput = 1",
  "OutputName = EntropyNorm",
  "OutputType = DomainIntegral",
  "DomainNorm = None",
  "UsesFlux = False",
  "ScalarName = Entropy",
  "ENDBLOCK",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitQ1Q2HexMesh
int
xf_UnitQ1Q2HexMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag){
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, Q1Q2HexGri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitQ1Q2HexAllCNS
int
xf_UnitQ1Q2HexAllCNS(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, Q1Q2HexGri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box3DCNSEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/* A .gri file for a Q1 quad mesh of a 2D box (9 elems) */
static char *BoxQ1Quad9Gri[] =
{
  "16 9 2",
  "0 0", "1 0", "2 0", "3 0",
  "0 1", "1 1", "2 1", "3 1",
  "0 2", "1 2", "2 2", "3 2",
  "0 3", "1 3", "2 3", "3 3",
  "4",
  "3 2 Left",
  "1 5","5 9", "9 13",
  "3 2 Right",
  "4 8","8 12", "12 16",
  "3 2 Bottom",
  "1 2","2 3", "3 4",
  "3 2 Top",
  "16 15","15 14", "14 13",
  "9 1 QuadLagrange",
  "1 2 5 6", "2 3 6 7", "3 4 7 8",
  "5 6 9 10", "6 7 10 11", "7 8 11 12",
  "9 10 13 14", "10 11 14 15", "11 12 15 16",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Quad9Mesh
int
xf_UnitBoxQ1Quad9Mesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag){
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Quad9Gri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Quad9All
int
xf_UnitBoxQ1Quad9All(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Quad9Gri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DScalarEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  (*pAll)->EqnSet->Dim = (*pAll)->Mesh->Dim; // Set EqnSet dimension

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Quad9EulerPenaltyAll
int
xf_UnitBoxQ1Quad9EulerPenaltyAll(xf_All **pAll)
{
  int ierr;
  
  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Quad9Gri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;
  
  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DEulerPenaltyEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;
  
  (*pAll)->EqnSet->Dim = (*pAll)->Mesh->Dim; // Set EqnSet dimension
  
  return xf_OK;
}



/* A .gri file for a Q1 quad mesh of a 2D box (9 elems); 2 element groups */
static char *BoxQ1Quad2EG9Gri[] =
{
  "16 9 2",
  "0 0", "1 0", "2 0", "3 0",
  "0 1", "1 1", "2 1", "3 1",
  "0 2", "1 2", "2 2", "3 2",
  "0 3", "1 3", "2 3", "3 3",
  "4",
  "3 2 Left",
  "1 5","5 9", "9 13",
  "3 2 Right",
  "4 8","8 12", "12 16",
  "3 2 Bottom",
  "1 2","2 3", "3 4",
  "3 2 Top",
  "16 15","15 14", "14 13",
  "4 1 QuadLagrange",
  "1 2 5 6", "2 3 6 7", "3 4 7 8", "5 6 9 10", 
  "5 1 QuadLagrange",
  "6 7 10 11", "7 8 11 12",
  "9 10 13 14", "10 11 14 15", "11 12 15 16",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Quad2EG9Mesh
int
xf_UnitBoxQ1QuadMesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag){
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Quad2EG9Gri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Quad2EG9All
int
xf_UnitBoxQ1Quad2EG9All(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Quad2EG9Gri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DScalarEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  (*pAll)->EqnSet->Dim = (*pAll)->Mesh->Dim; // Set EqnSet dimension

  return xf_OK;
}


/* A .gri file for a Q1 hex mesh of a 3D box (27 elems) */
static char *BoxQ1Hex27Gri[] =
{
  "64 27 3",
  " 0 0 0",  " 1 0 0",  " 2 0 0",  " 3 0 0",  " 0 1 0",
  " 1 1 0",  " 2 1 0",  " 3 1 0",  " 0 2 0",  " 1 2 0",
  " 2 2 0",  " 3 2 0",  " 0 3 0",  " 1 3 0",  " 2 3 0",
  " 3 3 0",  " 0 0 1",  " 1 0 1",  " 2 0 1",  " 3 0 1",
  " 0 1 1",  " 1 1 1",  " 2 1 1",  " 3 1 1",  " 0 2 1",
  " 1 2 1",  " 2 2 1",  " 3 2 1",  " 0 3 1",  " 1 3 1",
  " 2 3 1",  " 3 3 1",  " 0 0 2",  " 1 0 2",  " 2 0 2",
  " 3 0 2",  " 0 1 2",  " 1 1 2",  " 2 1 2",  " 3 1 2",
  " 0 2 2",  " 1 2 2",  " 2 2 2",  " 3 2 2",  " 0 3 2",
  " 1 3 2",  " 2 3 2",  " 3 3 2",  " 0 0 3",  " 1 0 3",
  " 2 0 3",  " 3 0 3",  " 0 1 3",  " 1 1 3",  " 2 1 3",
  " 3 1 3",  " 0 2 3",  " 1 2 3",  " 2 2 3",  " 3 2 3",
  " 0 3 3",  " 1 3 3",  " 2 3 3",  " 3 3 3",
  "6",
  "9 4 Left",
  "   1    5   21   17",
  "  17   21   37   33",
  "  33   37   53   49",
  "   5    9   25   21",
  "  37   41   57   53",
  "   9   13   29   25",
  "  25   29   45   41",
  "  41   45   61   57",
  "  21   25   41   37",
  "9 4 Right",
  "   4    8   24   20",
  "  20   24   40   36",
  "  36   40   56   52",
  "   8   12   28   24",
  "  24   28   44   40",
  "  40   44   60   56",
  "  12   16   32   28",
  "  28   32   48   44",
  "  44   48   64   60",
  "9 4 Bottom",
  "   1    2    6    5",
  "   5    6   10    9",
  "   9   10   14   13",
  "   2    3    7    6",
  "   6    7   11   10",
  "  10   11   15   14",
  "   3    4    8    7",
  "   7    8   12   11",
  "  11   12   16   15",
  "9 4 Top",
  "  49   50   54   53",
  "  53   54   58   57",
  "  57   58   62   61",
  "  50   51   55   54",
  "  54   55   59   58",
  "  58   59   63   62",
  "  51   52   56   55",
  "  55   56   60   59",
  "  59   60   64   63",
  "9 4 Front",
  "   1    2   18   17",
  "  17   18   34   33",
  "  33   34   50   49",
  "   2    3   19   18",
  "  18   19   35   34",
  "  34   35   51   50",
  "   3    4   20   19",
  "  19   20   36   35",
  "  35   36   52   51",
  "9 4 Back",
  "  13   14   30   29",
  "  29   30   46   45",
  "  45   46   62   61",
  "  14   15   31   30",
  "  30   31   47   46",
  "  46   47   63   62",
  "  15   16   32   31",
  "  31   32   48   47",
  "  47   48   64   63",
  "27 1 HexLagrange",
  " 1  2  5  6 17 18 21 22",
  " 2  3  6  7 18 19 22 23",
  " 3  4  7  8 19 20 23 24",
  " 5  6  9 10 21 22 25 26",
  " 6  7 10 11 22 23 26 27",
  " 7  8 11 12 23 24 27 28",
  " 9 10 13 14 25 26 29 30",
  "10 11 14 15 26 27 30 31",
  "11 12 15 16 27 28 31 32",
  "17 18 21 22 33 34 37 38",
  "18 19 22 23 34 35 38 39",
  "19 20 23 24 35 36 39 40",
  "21 22 25 26 37 38 41 42",
  "22 23 26 27 38 39 42 43",
  "23 24 27 28 39 40 43 44",
  "25 26 29 30 41 42 45 46",
  "26 27 30 31 42 43 46 47",
  "27 28 31 32 43 44 47 48",
  "33 34 37 38 49 50 53 54",
  "34 35 38 39 50 51 54 55",
  "35 36 39 40 51 52 55 56",
  "37 38 41 42 53 54 57 58",
  "38 39 42 43 54 55 58 59",
  "39 40 43 44 55 56 59 60",
  "41 42 45 46 57 58 61 62",
  "42 43 46 47 58 59 62 63",
  "43 44 47 48 59 60 63 64",
  "\0",
};



/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Hex27Mesh
int
xf_UnitBoxQ1Hex27Mesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag){
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Hex27Gri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Hex27All
int
xf_UnitBoxQ1Hex27All(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Hex27Gri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DScalarEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  (*pAll)->EqnSet->Dim = (*pAll)->Mesh->Dim; // Set EqnSet dimension

  return xf_OK;
}


/* A .gri file for a Q1 hex mesh of a 2D box (27 elems): 2 element groups */
static char *BoxQ1Hex2EG27Gri[] =
{
  "64 27 3",
  " 0 0 0",  " 1 0 0",  " 2 0 0",  " 3 0 0",  " 0 1 0",
  " 1 1 0",  " 2 1 0",  " 3 1 0",  " 0 2 0",  " 1 2 0",
  " 2 2 0",  " 3 2 0",  " 0 3 0",  " 1 3 0",  " 2 3 0",
  " 3 3 0",  " 0 0 1",  " 1 0 1",  " 2 0 1",  " 3 0 1",
  " 0 1 1",  " 1 1 1",  " 2 1 1",  " 3 1 1",  " 0 2 1",
  " 1 2 1",  " 2 2 1",  " 3 2 1",  " 0 3 1",  " 1 3 1",
  " 2 3 1",  " 3 3 1",  " 0 0 2",  " 1 0 2",  " 2 0 2",
  " 3 0 2",  " 0 1 2",  " 1 1 2",  " 2 1 2",  " 3 1 2",
  " 0 2 2",  " 1 2 2",  " 2 2 2",  " 3 2 2",  " 0 3 2",
  " 1 3 2",  " 2 3 2",  " 3 3 2",  " 0 0 3",  " 1 0 3",
  " 2 0 3",  " 3 0 3",  " 0 1 3",  " 1 1 3",  " 2 1 3",
  " 3 1 3",  " 0 2 3",  " 1 2 3",  " 2 2 3",  " 3 2 3",
  " 0 3 3",  " 1 3 3",  " 2 3 3",  " 3 3 3",
  "6",
  "9 4 Left",
  "   1    5   21   17",
  "  17   21   37   33",
  "  33   37   53   49",
  "   5    9   25   21",
  "  37   41   57   53",
  "   9   13   29   25",
  "  25   29   45   41",
  "  41   45   61   57",
  "  21   25   41   37",
  "9 4 Right",
  "   4    8   24   20",
  "  20   24   40   36",
  "  36   40   56   52",
  "   8   12   28   24",
  "  24   28   44   40",
  "  40   44   60   56",
  "  12   16   32   28",
  "  28   32   48   44",
  "  44   48   64   60",
  "9 4 Bottom",
  "   1    2    6    5",
  "   5    6   10    9",
  "   9   10   14   13",
  "   2    3    7    6",
  "   6    7   11   10",
  "  10   11   15   14",
  "   3    4    8    7",
  "   7    8   12   11",
  "  11   12   16   15",
  "9 4 Top",
  "  49   50   54   53",
  "  53   54   58   57",
  "  57   58   62   61",
  "  50   51   55   54",
  "  54   55   59   58",
  "  58   59   63   62",
  "  51   52   56   55",
  "  55   56   60   59",
  "  59   60   64   63",
  "9 4 Front",
  "   1    2   18   17",
  "  17   18   34   33",
  "  33   34   50   49",
  "   2    3   19   18",
  "  18   19   35   34",
  "  34   35   51   50",
  "   3    4   20   19",
  "  19   20   36   35",
  "  35   36   52   51",
  "9 4 Back",
  "  13   14   30   29",
  "  29   30   46   45",
  "  45   46   62   61",
  "  14   15   31   30",
  "  30   31   47   46",
  "  46   47   63   62",
  "  15   16   32   31",
  "  31   32   48   47",
  "  47   48   64   63",
  "13 1 HexLagrange",
  " 1  2  5  6 17 18 21 22",
  " 2  3  6  7 18 19 22 23",
  " 3  4  7  8 19 20 23 24",
  " 5  6  9 10 21 22 25 26",
  " 6  7 10 11 22 23 26 27",
  " 7  8 11 12 23 24 27 28",
  " 9 10 13 14 25 26 29 30",
  "10 11 14 15 26 27 30 31",
  "11 12 15 16 27 28 31 32",
  "17 18 21 22 33 34 37 38",
  "18 19 22 23 34 35 38 39",
  "19 20 23 24 35 36 39 40",
  "21 22 25 26 37 38 41 42",
  "14 1 HexLagrange",
  "22 23 26 27 38 39 42 43",
  "23 24 27 28 39 40 43 44",
  "25 26 29 30 41 42 45 46",
  "26 27 30 31 42 43 46 47",
  "27 28 31 32 43 44 47 48",
  "33 34 37 38 49 50 53 54",
  "34 35 38 39 50 51 54 55",
  "35 36 39 40 51 52 55 56",
  "37 38 41 42 53 54 57 58",
  "38 39 42 43 54 55 58 59",
  "39 40 43 44 55 56 59 60",
  "41 42 45 46 57 58 61 62",
  "42 43 46 47 58 59 62 63",
  "43 44 47 48 59 60 63 64",
  "\0",
};



/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Hex2EG27Mesh
int
xf_UnitBoxQ1Hex2EG27Mesh(xf_Mesh **pMesh, enum xfe_Bool CreateFlag){
  int ierr;
  
  if (CreateFlag){
    ierr = xf_Error(xf_CreateMesh(pMesh));
    if (ierr != xf_OK) return ierr;
  }

  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Hex2EG27Gri, (*pMesh)) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnitBoxQ1Hex2EG27All
int
xf_UnitBoxQ1Hex2EG27All(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, BoxQ1Hex2EG27Gri, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box3DCNSEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  (*pAll)->EqnSet->Dim = (*pAll)->Mesh->Dim; // Set EqnSet dimension

  return xf_OK;
}


/* A .gri file for an 8-element interval */
static char *Interval8[] =
{
  "9 8 1",
  "0", "1", "2", "3", "4", "5", "6", "7", "8",
  "2",
  "1 1 Left",
  "1",
  "1 1 Right",
  "9",
  "8 1 SegLagrange",
  "1 2", "2 3", "3 4", "4 5", "5 6", "6 7", "7 8", "8 9",
  "\0",
};

/* An Euler .eqn file for 1d interval mesh */
static char *Interval8EulerEqn[] =
{
  "EqnSetLibrary = libCompressibleNS.so",
  // Parameters
  "STARTBLOCK PARAM",
  "SpecificHeatRatio = 1.4",
  "GasConstant = 1.0",
  "Viscosity = .1",
  "RegularityScalar = Density",
  "ENDBLOCK",
  // Residual Terms
  "STARTBLOCK RESTERM",
  "nResTerm = 1",
  "TermType = Convection",
  "InteriorFlux = Standard",
  "FluxFunction = Roe",
  "ENDBLOCK",
  // Initial Conditions
  "STARTBLOCK IC",
  "ICType = FullState",
  "Data = 1.0 0.59160797830996159 2.675",
  "ENDBLOCK",
  // Boundary Conditions
  "STARTBLOCK BC",
  "nBC = 2",
  "BFGTitle = Left",  "BCType = FullState", "Data = 1.0 0.5 2.675",
  "BFGTitle = Right", "BCType = FullState", "Data = 1.0 0.6 2.675",
  "ENDBLOCK",
  // Outputs
  "STARTBLOCK OUTPUT",
  "nOutput = 2",
  "OutputName = EntropyNorm",
  "OutputType = DomainIntegral",
  "DomainNorm = None",
  "UsesFlux = False",
  "ScalarName = Entropy",
  "OutputName = Point",
  "OutputType = PointValue",
  "TimeNorm = Final",
  "egrp = 0",
  "elem = 3",
  "xref = 0.1",
  "ScalarName = Density",
  "ENDBLOCK",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitIntervalAllEuler
int
xf_UnitIntervalAllEuler(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, Interval8, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Interval8EulerEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/* A .gri file for an 5-element skewed tri mesh */
static char *SkewedQ1Tri5[] =
{
  "6 5 2",
  "0 0", "2 0",  "1 1",  "0 2",  "1 2",  "2 2",
  "4",
  "1 2 Bottom",  "1 2",
  "1 2 Right",   "2 6",
  "2 2 Top",     "4 5",  "5 6",
  "1 2 Left",    "1 4",
  "5 1 TriLagrange",
  "1 5 4",  "1 3 5",  "1 2 3",  "2 5 3",  "2 6 5",
  "\0",
};


/******************************************************************/
//   FUNCTION Definition: xf_UnitSkewedQ1Tri5AllEuler
int
xf_UnitSkewedQ1Tri5AllEuler(xf_All **pAll)
{
  int ierr;

  ierr = xf_Error(xf_CreateAll(pAll, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Mesh
  ierr = xf_Error( xf_ReadGriFile(NULL, SkewedQ1Tri5, (*pAll)->Mesh) );
  if (ierr != xf_OK) return ierr;

  // EqnSet
  ierr = xf_Error( xf_ReadEqnSetFile(NULL, Box2DEulerEqn, (*pAll)->EqnSet) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
