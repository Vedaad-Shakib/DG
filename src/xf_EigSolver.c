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
  FILE:  xf_EigSolver.c

  This file contains eigenvalue solver functions for sparse systems.
  Note, some of these functions use Lapack-dependent functions from
  MathLapack.c.  Hence, Lapack is required for their use.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_EigSolverStruct.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_MathLapack.h"
#include "xf_Solver.h"
#include "xf_MPI.h"
#include "xf_String.h"


/******************************************************************/
//   FUNCTION Definition: xf_CreateEigSolverData
static int
xf_CreateEigSolverData(int ncv, xf_EigSolverData **pEigSolverData)
{
/*
PURPOSE:

  Creates a structure for storing data required by the eigensolver.
  Currently, only a Lanczos iteration is supported, and the required
  storage consists of a main (alpha) and sub (beta) diagonal of the
  symmetric, tridiagonal Lanczos matrix, H.

INPUTS:

  ncv : maximum expected number of Lanczos vectors; i.e. the Lanczos
        matrix will be at most ncv x ncv
  
OUTPUTS:

  (*pEigSolverData) : allocated EigSolverData structure

RETURN:

  Error Code

*/
  int ierr, k;

  ierr = xf_Error(xf_Alloc((void **) pEigSolverData, 1,
			   sizeof(xf_EigSolverData)));
  if (ierr != xf_OK) return ierr;
  
  (*pEigSolverData)->iIter      = 0;
  (*pEigSolverData)->nRestart   = 0;
  (*pEigSolverData)->nInvariant = 0;
  (*pEigSolverData)->ncv        = ncv;

  (*pEigSolverData)->InvariantSubspaceFound = xfe_False;

  ierr = xf_Error(xf_Alloc((void **) &(*pEigSolverData)->H, 2*ncv,
			   sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for (k=0; k<2*ncv; k++) (*pEigSolverData)->H[k] = 0.0;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyEigSolverData
static int
xf_DestroyEigSolverData(xf_EigSolverData *EigSolverData)
{
/*
PURPOSE:  

  Destroys an EigSolverData structure.

INPUTS:  
   
  EigSolverData : structure to be destroyed

OUTPUTS: None

RETURN:

  Error Code

*/
  xf_Release( (void *) EigSolverData->H);
  xf_Release( (void *) EigSolverData);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ApplyShifts
static int
xf_ApplyShifts(xf_VectorSet *VS, int nk, int np, real *shift,
	       real *alpha, real *beta)
{
/*
PURPOSE:

  Translated version of the function dsapps.f in ARPACK.

  Given the Lanczos factorization

     A*V_m = V_m*H_m + r_m*e_m^T,

  m = k+p, applies shifts implicitly resulting in

     A*(V_m*Q) = (V_m*Q)*(Q^T* H_m*Q) + r_m*e_m^T * Q

  where Q is an orthogonal matrix of size m x m. Q is the product of
  rotations resulting from np bulge chasing sweeps.  The updated
  Lanczos factorization becomes:

     A*Vnew_k = Vnew_k*Hnew_k + rnew_k*e_k^T.

  This new Lanczos factorization is only valid for the first k-1
  vectors in Vnew_k.  That is, to continue the Lanczos factorization,
  begin with multiplying Vnew_{k-1} by the operator.

  Note, the Lanczos vectors are stored in VS, and must have been all
  computed before this function is called.  The Lanczos matrix, H, is
  passed in as a main and sub diagonal (alpha and beta).

INPUTS:

  VS    : VectorSet containing np+nk Lanczos vectors (V_m)
  nk    : number of desired eigs (k in the above description)
  np    : number of shifts (p in the above description)
  shift : vector of shifts [np]
  alpha : main diagonal of H_m
  beta  : subdiagonal of H_m
  
OUTPUTS:

  VS    : modified Lanczos vectors (Vnew_k)
  alpha : main diagonal of Hnew_k
  beta  : sub-diagonal of Hnew_k
  

RETURN:

  Error Code

*/
  int ierr, i, j, iv, nm;
  real *Q;
  xf_Vector V0;
  
  nm = nk + np;

  // quick return if possible
  if (np == 0) return xf_OK;

  if ((nk <= 0) || (nm <= 0) || (np <= 0)) return xf_Error(xf_INPUT_ERROR);

  // Allocate Q [nm x nm]
  ierr = xf_Error(xf_Alloc( (void **) &Q, nm*nm, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // Calculate Q and apply it to [alpha,beta] in the process
  ierr = xf_Error(xf_QRBulgeChase(shift, nk, np, alpha, beta, Q));
  if (ierr != xf_OK) return ierr;

  /* Set V <-- V*Q (columns 1 to nk), in backward order; result is
     stored in vectors nm-nk to nm-1*/
  for (i=0; i<nk; i++){
    iv = nm-1-i;
    ierr = xf_Error(xf_VectorMult(VS->Vector+iv, Q[iv*nm+nk-i-1]));
    if (ierr != xf_OK) return ierr;
    for (j=1; j<nm-i; j++){
      ierr = xf_Error(xf_VectorMultSet(VS->Vector+iv-j, Q[(iv-j)*nm+nk-i-1],
				       xfe_Add, VS->Vector+iv));
      if (ierr != xf_OK) return ierr;
    }
  }

  // Copy vectors [nm-nk:nm-1] to [1:nk]
  for (i=0; i<nk; i++)
    swap(VS->Vector[i], VS->Vector[nm-nk+i], V0);
  
  xf_Release((void *) Q);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_OrthogonalizeVector
static int
xf_OrthogonalizeVector(xf_VectorSet *VS, int n, enum xfe_Bool debug,
		       xf_Vector *W, real *dalpha)
{
  /* This function orthogonalizes W against the first n vectors in VS.
     Optionally, dalpha = <VS+n-1,W> is returned. */

  int ierr, j;
  real *s;
    
  ierr = xf_Error(xf_Alloc( (void **) &s, n, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  for (j=0; j<n; j++){
    ierr = xf_Error(xf_VectorDot(VS->Vector+j, W, s+j));
    if (ierr != xf_OK) return ierr;
    if (debug)
      xf_printf("  In Orthogonalize, n=%d; s[%d] = %.10E\n", n, j, s[j]);
  }
  for (j=0; j<n; j++){
    ierr = xf_Error(xf_VectorMultSet(VS->Vector+j, s[j], xfe_Sub, W));
    if (ierr != xf_OK) return ierr;
  }

  if (dalpha != NULL) (*dalpha) = s[n-1];
  
  xf_Release( (void *) s);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SortEigenvalues
static void
xf_SortEigenvalues(int nev, int ncv, real *ritz, real *bounds, int *P)
{
  /* This function uses a bubble sort to find the nev largest
     magnitude eigenvalues from the vector ritz.  These nev largest
     values are placed at the end of ritz, sorted so that the largest
     value is in the last entry of ritz.  bounds and P are optional
     real and integer vectors associated with the eigenvalues; if
     provided, they are ordered to maintain correspondence with the
     ritz vector. */
  int np, i, j, k;
  real val;

  np = ncv - nev;
  for (i=ncv-1; i>=ncv-nev; i--){
    for (j=i-1, k=i; j>=0; j--){
      if (fabs(ritz[j]) > fabs(ritz[k])) k = j;
    }
    if (i != k){
      swap(ritz[i], ritz[k], val);
      if (bounds != NULL) swap(bounds[i], bounds[k], val);
      if (P != NULL) swap(P[i], P[k], j);
    }
  }
}

/******************************************************************/
//   FUNCTION Definition: xf_LanczosReadPoint
static int
xf_LanczosReadPoint(xf_All *All, xf_VectorSet *VS, xf_EigSolverData *EigSolverData)
{
  int ierr, i, ncv, nLanczos;
  int myRank, nProc, terr;
  char line[xf_MAXLINELEN];
  FILE *fid;
  real *alpha, *beta, *H;
  xf_DataSet *DataSet;
  xf_Data *D;
  xf_VectorSet *LanczosSet;


  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  xf_printf("Restarting from existing Lanczos.data and Lanczos.txt\n");
  
  ncv  = EigSolverData->ncv;  // number of Lanczos vectors for current run

  // set alpha and beta
  H = EigSolverData->H;
  beta  = H+0;
  alpha = H+ncv;

  // Read in Lanczos.txt for H matrix
  ierr = xf_Error(xf_fopen("Lanczos.txt", "r", &fid));
  if (ierr != xf_OK) return ierr;

  terr = xf_OK;
  if (myRank == 0){
    if ((fid = fopen("Lanczos.txt", "r")) == NULL) terr = xf_Error(xf_FILE_READ_ERROR);
    if (terr == xf_OK)
      if (fgets(line, xf_MAXLINELEN, fid) == NULL) terr = xf_Error(xf_FILE_READ_ERROR);
    if (terr == xf_OK)
      if (sscanf(line, "%d", &nLanczos) != 1) terr = xf_Error(xf_FILE_READ_ERROR);
    
    // error out if savepoint file has more Lanczos vectors than we can handle
    if (terr == xf_OK)
      if (nLanczos > ncv) terr = xf_Error(xf_OUT_OF_BOUNDS);

    if (terr == xf_OK)
      for (i=0; i<nLanczos; i++){
	if (fgets(line, xf_MAXLINELEN, fid) == NULL) terr = xf_Error(xf_FILE_READ_ERROR);
	if (terr != xf_OK) break;
	if (sscanf(line, "%lf %lf", alpha+i, beta+i) != 2) 
	  return xf_Error(xf_FILE_READ_ERROR);
      }
    fclose(fid);
  }

  // broadcast terr and quit if not OK
  if (xf_PError(&terr, 0) != xf_OK) return terr;
  
  // broadcast info to all procs
  ierr = xf_Error(xf_MPI_Bcast((void *) &nLanczos, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_MPI_Bcast((void *) alpha, nLanczos*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_MPI_Bcast((void *) beta, nLanczos*sizeof(real), 0));
  if (ierr != xf_OK) return ierr;

  EigSolverData->iIter = nLanczos;

  // create dataset for Lanczos
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;

  // read in Lanczos.data
  ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, "Lanczos.data", DataSet));
  if (ierr != xf_OK) return ierr;


  D = DataSet->Head;
  if (strncmp(D->Title, "LanczosSet", 10) != 0) return xf_Error(xf_INPUT_ERROR);
  LanczosSet = (xf_VectorSet *) D->Data;
  
  // # vectors in Lanczos.data must be equal to or more than nLanczos
  if (nLanczos > LanczosSet->nVector) return xf_Error(xf_OUT_OF_BOUNDS);

  // Copy first nLanczos Lanczos vectors to VS
  for (i=0; i<nLanczos; i++){
    ierr = xf_Error(xf_SetVector(LanczosSet->Vector+i, xfe_Set, VS->Vector+i));
    if (ierr != xf_OK) return ierr;
  }

  // destroy DataSet
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;

  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_LanczosSavePoint
static int
xf_LanczosSavePoint(xf_All *All, xf_VectorSet *VS, xf_EigSolverData *EigSolverData)
{
  int ierr, i, nLanczos;
  int myRank, nProc;
  char line[xf_MAXLINELEN];
  FILE *fid;
  real *alpha, *beta, *H;
  xf_DataSet *DataSet, *DataSet_Glob = NULL;
  xf_Data *D;
  xf_VectorSet *LanczosSet;


  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  nLanczos  = EigSolverData->iIter-1;  // number of Lanczos vectors for current run

  xf_printf("Saving Lanczos.data and Lanczos.txt: nLanczos = %d\n", nLanczos);
  

  // set alpha and beta
  H = EigSolverData->H;
  beta  = H+0;
  alpha = H+EigSolverData->ncv;

  // Write Lanczos.txt for H matrix
  if (myRank == 0){
    if ((fid = fopen("Lanczos.txt", "w")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    fprintf(fid, "%d\n", nLanczos);
    for (i=0; i<nLanczos; i++)
      fprintf(fid, "%.15E %.15E\n", alpha[i], beta[i]);
    fclose(fid);
  }
  
  // Write out Lanczos.data, containing VS
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_DataSetAdd(DataSet, "LanczosSet", xfe_VectorSet,
				xfe_True, (void *) VS, &D));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet, NULL, "Lanczos.data"));
  if (ierr != xf_OK) return ierr;
  D->Data = NULL;
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_EigIterLanczos
int
xf_EigIterLanczos(xf_All *All, int nev, int ncv, xf_VectorSet *VS, real tol, 
		  enum xfe_Verbosity Verbosity, int nev_restart, 
		  enum xfe_Bool SavePoint, enum xfe_Bool LanczosRestart,
		  real *E, xf_VectorSet *EV, xf_Vector **pV, 
		  xf_Vector **pW, xf_EigSolverData **pEigSolverData, 
		  enum xfe_EigStatusType *Status)
{
  /*
    Implicitly-Restarted Arnoldi-Method (symmetric version) modeled on
    ARPACK.  For documentation, see xf_EigSolver.h.
  */
  int ierr, iIter, i, j, k, np, iOrtho, nev0;
  int nconv, nptemp, nevbef, *P;
  int nRestartMax = 30;
  enum xfe_Bool debug;
  enum xfe_Bool ForceExit = xfe_False;
  real wnorm, rnorm, rnorm1, val, dalpha;
  real *alpha, *beta, *H, *s, *Z;
  real *ritz, *bounds, *workl;
  xf_EigSolverData *EigSolverData;
  xf_Vector *V, *W;

  // print debug results when high verbosity is requested
  debug = (Verbosity == xfe_VerbosityHigh);

  /* On first iteration, initialize EigSolverData, generate a random
     vector, ask for a mult that will place this vector into the range
     of A, and return */
  if ((*pEigSolverData) == NULL){
    
    // check input
    if ((nev <= 0) || (ncv <= nev) || (VS->nVector <= ncv))
      return xf_Error(xf_INPUT_ERROR);

    // allocate EigSolverData
    ierr = xf_Error(xf_CreateEigSolverData(ncv, pEigSolverData));
    if (ierr != xf_OK) return ierr;
    EigSolverData = (*pEigSolverData);

    // generate a random vector in VS->Vector[ncv]
    ierr = xf_Error(xf_VectorRand(VS->Vector + ncv, 0));
    if (ierr != xf_OK) return ierr;

    // Set Status and pV
    (*Status) = xfe_EigMultiply;
    (*pV) = VS->Vector + ncv;

    /* If restarting from Lanczos data, we pick up where we left off */
    if ((LanczosRestart) && (All != NULL)){
      ierr = xf_Error(xf_LanczosReadPoint(All, VS, EigSolverData));
      if (ierr != xf_OK) return ierr;
      (*pV) = VS->Vector + EigSolverData->iIter-1;
      (*pW) = VS->Vector + EigSolverData->iIter;
      return xf_OK;
    }

    /* if restarting from existing eigenvectors, we effectively have an
       invariant subspace. */
    if (nev_restart > 0){
      xf_printf("Initializing EigSolver to restart from %d existing vectors.\n",
		nev_restart);
      H = EigSolverData->H;
      beta  = H + 0;
      alpha = H + ncv;
      for (iIter=0; iIter<nev_restart; iIter++){
	alpha[iIter] = E[iIter];
	beta[iIter] = 0.0;
      }
      beta[nev_restart] = 0;
      EigSolverData->iIter = nev_restart;
      EigSolverData->InvariantSubspaceFound = xfe_True;
      (*pW) = VS->Vector + nev_restart;
    }
    else
      (*pW) = VS->Vector + 0;
    
    return xf_OK;
  }

  EigSolverData = (*pEigSolverData);
  W = (*pW);
  V = (*pV);
  iIter = EigSolverData->iIter++;
  H = EigSolverData->H;
  beta  = H + 0;
  alpha = H + ncv;
  ritz = E;
  ForceExit = xfe_False;
  nev0 = nev; // original # of requested eigs (extras will be padded with 0s)

  /* If an invariant subspace was found on the last iteration, W
     contains A*x, where x is a random vector.  Need to orthogonalize
     W against the first iIter vectors in VS.  */
  if (EigSolverData->InvariantSubspaceFound){
    
    ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm));
    if (ierr != xf_OK) return ierr;

    if (debug)
      xf_printf("Back after Invariant subspace found.  Reorthogonalizing. iIter = %d\n", iIter);
    ierr = xf_Error(xf_OrthogonalizeVector(VS, iIter, debug, W, NULL));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm1));
    if (ierr != xf_OK) return ierr;

    if (rnorm1 < 0.717*rnorm){ // possibly no other directions, try reortho again
      rnorm = rnorm1;
      ierr = xf_Error(xf_OrthogonalizeVector(VS, iIter, debug, W, NULL));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm1));
      if (ierr != xf_OK) return ierr;
      if ((rnorm1 < 0.717*rnorm) || (rnorm1 < 10.0*MEPS)){
	if (debug)
	  xf_printf("Forcing exit because rnorm = %.10E rnorm1 = %.10E\n", rnorm, rnorm1);
	ForceExit = xfe_True; // no other directions to search; exit.
      }
    }
  }

  // wnorm = ||W||
  ierr = xf_Error(xf_VectorNorm(W, 2, &wnorm));
  if (ierr != xf_OK) return ierr;

  if (debug) xf_printf("wnorm = %.10E\n", wnorm);
  if (wnorm < MEPS) ForceExit = xfe_True;

  if (ForceExit){
    /* A forced exit occurs when an invariant subspace of eigenvectors
       has been found and when that subspace contains all the nonzero
       eigenvalues.  i.e. all other search directions (operator
       applied to a random vector) result in a vector in the span of
       the invariant subspace. */
    if (debug) 
      xf_printf("Forced exit at iteration iIter = %d\n", iIter);
    ncv = min(nev,iIter); // so we can calculate the eigs of H with the same code
    nev = ncv; // # of nonzero eigenvalues returned
    rnorm = 0.0; // residual norm is zero
  }
  else{
    /* On first iteration, have only one vector, so ask for another
       mult.  If an invariant subspace was found, we're also basically
       starting over, so we need another mult. */
    if ((iIter == 0) || (EigSolverData->InvariantSubspaceFound)){
      beta[iIter] = 0.0;
      EigSolverData->InvariantSubspaceFound = xfe_False;

      // W = W / wnorm
      ierr = xf_Error(xf_VectorMult(W, 1.0/wnorm));
      if (ierr != xf_OK) return ierr;

      (*Status) = xfe_EigMultiply;
      (*pV) = VS->Vector + iIter;
      (*pW) = VS->Vector + iIter+1;
      return xf_OK;
    }

    // at this point, iIter >= 1, V = VS(iIter-1), W = VS(iIter)
  
    // alpha(iIter-1) = V^T * W
    ierr = xf_Error(xf_VectorDot(V, W, alpha+iIter-1));
    if (ierr != xf_OK) return ierr;

    // W = W - beta(iIter-1)*VS(iIter-2) - alpha(iIter-1)*VS(iIter-1)
    if (iIter > 1){
      ierr = xf_Error(xf_VectorMultSet(VS->Vector+iIter-2, beta[iIter-1], xfe_Sub, W));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_VectorMultSet(V, alpha[iIter-1], xfe_Sub, W));
    if (ierr != xf_OK) return ierr;
    
    // rnorm = ||W||
    ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm));
    if (ierr != xf_OK) return ierr;

    if (debug)
      xf_printf("iIter = %d; wnorm = %.10E, rnorm = %.15E\n", iIter, wnorm, rnorm);
  
    // if rnorm < 0.717*wnorm, reorthogonalize (following ARPACK)
    rnorm1 = rnorm;
    rnorm  = wnorm;
    for (iOrtho=0; iOrtho<2; iOrtho++){ // reorthogonalize up to two times
      if (rnorm1 < 0.717*rnorm){ 
	if (debug){
	  xf_printf("Reorthogonalizing: iOrtho = %d\n", iOrtho);
	  xf_printf("rnorm = %.10E, rnorm1 = %.10E\n", rnorm, rnorm1);
	}

	// Call orthogonalization function
	ierr = xf_Error(xf_OrthogonalizeVector(VS, iIter, debug, W, &dalpha));
	if (ierr != xf_OK) return ierr;

	alpha[iIter-1] += dalpha;
      
	rnorm = rnorm1;
      
	// rnorm1 = ||W||
	ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm1));
	if (ierr != xf_OK) return ierr;
      }
    } // iOrtho

    /* If two attempts at reorthogonalization did not help, an invariant
       subspace has been found.  Following ARPACK, generate a new random
       vector, A*rand; make sure to orthogonalize it 
    */
    if (debug) 
      xf_printf("rnorm = %.10E, rnorm1 = %.10E\n", rnorm, rnorm1);
    if (((rnorm1 < 0.717*rnorm) || (rnorm1 < 10.0*MEPS)) && (iIter < ncv)){
    
      if (debug)
	xf_printf("Invariant subspace found. iIter = %d\n", iIter);
    
      beta[iIter] = 0.0;

      EigSolverData->InvariantSubspaceFound = xfe_True;
      EigSolverData->iIter--; // so that iteration number stays the same
      if ((EigSolverData->nInvariant++) > nRestartMax){
	xf_printf("Lanczos stalled after trying to get away from an invariant subspace.\n");
	xf_printf(" %d attempts were made to try to find a new direction.\n", nRestartMax);
	return xf_Error(xf_NOT_CONVERGED);
      }
    
      // generate a random vector in VS->Vector[ncv]
      ierr = xf_Error(xf_VectorRand(VS->Vector + ncv, 0));
      if (ierr != xf_OK) return ierr;
    
      // Set pV and pW
      (*Status) = xfe_EigMultiply;
      (*pV) = VS->Vector + ncv;
      (*pW) = VS->Vector + iIter;
    
      return xf_OK;
    }

    if (iIter < ncv){
      beta[iIter] = rnorm; // set beta(iIter)
    
      // W = W / rnorm
      ierr = xf_Error(xf_VectorMult(W, 1.0/rnorm));
      if (ierr != xf_OK) return ierr;
    }


    /* Save away Lanczos data for immediate restart */
    if ((SavePoint) && (All != NULL)){
      ierr = xf_Error(xf_LanczosSavePoint(All, VS, EigSolverData));
      if (ierr != xf_OK) return ierr;
    }


    /* If don't have enough vectors, ask for another mult and return */
    if (iIter < ncv){
      (*Status) = xfe_EigMultiply;
      (*pV) = VS->Vector + iIter;
      (*pW) = VS->Vector + iIter+1;
      return xf_OK;
    }
  } // end else if no ForceExit


  /* Allocate memory for ritz value calculation */
  ierr = xf_Error(xf_Alloc( (void **) &Z, ncv*ncv, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc( (void **) &workl, ncv, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc( (void **) &bounds, ncv, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc( (void **) &P, ncv, sizeof(int)));
  if (ierr != xf_OK) return ierr; // permutation vector for sorting
  for (k=0; k<ncv; k++) P[k] = k;

  /* Copy H matrix into ritz (main diagonal) and workl (sub-diagonal) */
  for (k=0; k<ncv; k++)   ritz[k]  = alpha[k];
  for (k=0; k<ncv-1; k++) workl[k] = beta[k+1];

  if (debug){
    xf_printf("Before call to EigSymTriDiag:\n");
    for (k=0; k<ncv; k++){
      xf_printf("alpha[%d] = %.15E, beta[%d] = %.15E\n", 
		k, alpha[k], k, beta[k]); 
    }
  }


  /* Compute the eigenvalues and eigenvectors for the curent symmetric
     tridiagonal matrix, H.  Note, this is where MathLapack is
     required. */
  ierr = xf_Error(xf_EigSymTriDiag(ncv, ritz, workl, Z));
  if (ierr != xf_OK) return ierr;

  /* Compute the error bounds for the Ritz values (eigs of H).  These
     are just the last entries in the eigenvectors, scaled by the
     residual norm, rnorm.  Note, the eigenvectors in Z are stored in
     column-order (Fortran style). */
  for (k=0; k<ncv; k++)
    bounds[k] = fabs(Z[ncv*k+ncv-1]*rnorm);

  if (debug)
    for (k=0; k<ncv; k++){
      xf_printf("ritz[%d] = %.10E, bounds[%d] = %.10E\n", 
		k, ritz[k], k, bounds[k]); 
    }
  
 
  /* Select wanted Ritz values (e.g. largest in magnitude). The
     unwanted Ritz values (i.e. shifts) are placed into the first np =
     (ncv-nev) locations in ritz, while the wanted values are placed
     into locations [np+1:ncv]. The bounds get sorted as well. */
  np = ncv - nev;
  xf_SortEigenvalues(nev, ncv, ritz, bounds, P);
  
  if (debug){
    xf_printf("After Sorting:\n");
    for (k=0; k<ncv; k++){
      xf_printf("ritz[%d] = %.15E, bounds[%d] = %.15E\n", 
		k, ritz[k], k, bounds[k]); 
    }
  }

  /* Check convergence of wanted Ritz values */
  nconv = 0; // # converged values
  for (i=ncv-nev; i<ncv; i++)
    nconv += (bounds[i] <= tol*max(MEPS, fabs(ritz[i])));

  if (Verbosity >= xfe_VerbosityMedium){
    xf_printf("[iRestart = %d] Lanczos eigenvalue estimates and bounds:\n",
	      EigSolverData->nRestart);
    for (i=ncv-nev; i<ncv; i++){
      xf_printf(" [%d] Eigenvalue = %.10E, Bound = %.10E\n", 
		i-ncv+nev, ritz[i], bounds[i]);
    }
  }
  
  /* Check if user requested a halt: stop computation then */
  if (xf_CheckUserHalt("EIGSTOP")){
    nconv = ncv;
    ForceExit = xfe_True;
  }


  /* Exit successfully if have enough converged values.  Converged
     eigs/bounds are placed in ritz(0:nconv-1)/bounds(0:nconv-1).
     Eigenvectors are computed if desired. */
  if ((nconv >= nev) || (ForceExit)){
    if (ForceExit) nconv = ncv; // should be true anyway since rnorm == 0.0
    for (i=0; i<nconv; i++){
      swap(ritz[i], ritz[ncv-nev+i], val);
      swap(bounds[i], bounds[ncv-nev+i], val);
      swap(P[i], P[ncv-nev+i], j);
    }
    for (i=nconv; i<nev0; i++) E[i] = 0.0;

    if (EV != NULL){ // Computation of eigenvectors
      for (i=0; i<nconv; i++){
	ierr = xf_Error(xf_SetZeroVector(EV->Vector + i));
	if (ierr != xf_OK) return ierr;
	for (j=0; j<ncv; j++){
	  ierr = xf_Error(xf_VectorMultSet(VS->Vector + j, Z[P[i]*ncv + j], 
					   xfe_Add, EV->Vector+i));
	  if (ierr != xf_OK) return ierr;
	}
      }
      for (i=nconv; i<nev0; i++){ // zero out extra eigenvectors if ForceExit
	ierr = xf_Error(xf_SetZeroVector(EV->Vector + i));
	if (ierr != xf_OK) return ierr;
      }
    }

    if ((ForceExit) && (nev < nev0))
      (*Status) = xfe_EigForceExit;
    else
      (*Status) = xfe_EigConverged;

    ierr = xf_Error(xf_DestroyEigSolverData(EigSolverData));
    if (ierr != xf_OK) return ierr;
  }
  else{

    /* Error out if applied too many restarts */ 
    if (EigSolverData->nRestart++ > nRestartMax){
      xf_printf("Lanczos not converged after %d restarts. Giving up.\n", nRestartMax);
      return xf_Error(xf_NOT_CONVERGED);
    }

    
    /* Count # of unwated Ritz values with zero bounds
       estimates. Decrease np, the number of shifts to apply. */
    for (i=0, nptemp = np; i<nptemp; i++)
      if (bounds[i] <= MEPS) np--;
      //if (bounds[i] == 0.) np--;
    
    if (nptemp != np){
      nev += nptemp-np;
      xf_printf("Decreased # shifts due to zero bounds estimates. np = %d, nev = %d\n", 
		np, nev);
    }

    /* Adjust nev to prevent stagnation */
    nevbef = nev;
    nev = nev + min (nconv, np/2);
    if ((nev == 1) && (ncv >= 6))
      nev = ncv / 2;
    else if ((nev == 1) && (ncv > 2))
      nev = 2;
    np = ncv - nev;

    // if nev has increased, resort eigenvalues
    if (nevbef < nev){
      if (debug) xf_printf("Resorting eigenvalues: nevbef = %d, nev = %d\n", nevbef, nev);
      xf_SortEigenvalues(nev, ncv, ritz, bounds, P);
      if (debug){
	xf_printf("After ReSorting:\n");
	for (k=0; k<ncv; k++){
	  xf_printf("ritz[%d] = %.15E, bounds[%d] = %.15E\n",
		    k, ritz[k], k, bounds[k]);
	}
      }
    }
    
    if (debug){
      xf_printf("Before shifts:\n");
      for (k=0; k<ncv; k++){
	xf_printf("alpha[%d] = %.15E, beta[%d] = %.15E\n", 
		  k, alpha[k], k, beta[k]); 
      }
    }
    
    /* Apply shifts -- separate function to do this */
    ierr = xf_Error(xf_ApplyShifts(VS, nev, np, ritz, alpha, beta));
    if (ierr != xf_OK) return ierr;
    
    if (debug){
      xf_printf("After shifts:\n");
      for (k=0; k<ncv; k++){
	xf_printf("alpha[%d] = %.15E, beta[%d] = %.15E\n", 
		  k, alpha[k], k, beta[k]); 
      }
    }
      
    if (debug)
      for (j=0; j<iIter-1; j++){
	ierr = xf_Error(xf_VectorDot(VS->Vector+j, VS->Vector+j+1, &val));
	if (ierr != xf_OK) return ierr;
	xf_printf(" dp[%d] = %.10E\n", j, val);
      }

    /* Set pointers for next call and return with a multiply request */
    (*Status) = xfe_EigMultiply;
    EigSolverData->iIter = nev;
    EigSolverData->nInvariant = 0; // reset counter upon restart
    (*pV) = VS->Vector + nev-1;
    (*pW) = VS->Vector + nev;
  }

  /* Release memory */
  xf_Release( (void *) bounds);
  xf_Release( (void *) workl);
  xf_Release( (void *) Z); 
  xf_Release( (void *) P);

  return xf_OK;

}




/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_EigIterLanczosPrepRestart */
/* int */
/* xf_EigIterLanczosPrepRestart(int nev_current, int ncv, xf_VectorSet *VS, */
/* 			     xf_Vector **pV, xf_Vector **pW,  */
/* 			     xf_EigSolverData **pEigSolverData,  */
/* 			     enum xfe_EigStatusType *Status) */
/* { */
/*   /\* */
/*     Prepares for a restart of the Lanczos eigenvalue solver. */
/*   *\/ */
/*   int ierr, iIter, i, j, k, np, iOrtho, nev0; */
/*   int nconv, nptemp, nevbef, *P; */
/*   int nRestartMax = 30; */
/*   enum xfe_Bool debug; */
/*   enum xfe_Bool ForceExit = xfe_False; */
/*   real wnorm, rnorm, rnorm1, val, dalpha; */
/*   real *alpha, *beta, *H, *s, *Z; */
/*   real *ritz, *bounds, *workl; */
/*   xf_EigSolverData *EigSolverData; */
/*   xf_Vector *V, *W; */

/*   // print debug results when high verbosity is requested */
/*   debug = (Verbosity == xfe_VerbosityHigh); */

/*   /\* On first iteration, initialize EigSolverData, generate a random */
/*      vector, ask for a mult that will place this vector into the range */
/*      of A, and return *\/ */
/*   if ((*pEigSolverData) != NULL) return xf_Error(xf_INPUT_ERROR); */
  
/*   // check input */
/*   if ((ncv <= nev_current) || (VS->nVector <= ncv)) */
/*     return xf_Error(xf_INPUT_ERROR); */
  
/*   // allocate EigSolverData */
/*   ierr = xf_Error(xf_CreateEigSolverData(ncv, pEigSolverData)); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   H = EigSolverData->H; */
/*   beta  = H + 0; */
/*   alpha = H + ncv; */

/*   beta[iIter] = 0.0; */
  
/*   EigSolverData->InvariantSubspaceFound = xfe_True; */
/*   EigSolverData->iIter--; // so that iteration number stays the same */
/*   if ((EigSolverData->nInvariant++) > nRestartMax){ */
/*     xf_printf("Lanczos stalled after trying to get away from an invariant subspace.\n"); */
/*     xf_printf(" %d attempts were made to try to find a new direction.\n", nRestartMax); */
/*     return xf_Error(xf_NOT_CONVERGED); */
/*   } */
  
/*   // generate a random vector in VS->Vector[ncv] */
/*   ierr = xf_Error(xf_VectorRand(VS->Vector + ncv, 0)); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   // Set pV and pW */
/*   (*Status) = xfe_EigMultiply; */
/*   (*pV) = VS->Vector + ncv; */
/*   (*pW) = VS->Vector + iIter; */
  
  

  
/*   return xf_OK; */








/*   EigSolverData = (*pEigSolverData); */
/*   W = (*pW); */
/*   V = (*pV); */
/*   iIter = EigSolverData->iIter++; */
/*   H = EigSolverData->H; */
/*   beta  = H + 0; */
/*   alpha = H + ncv; */
/*   ritz = E; */
/*   ForceExit = xfe_False; */
/*   nev0 = nev; // original # of requested eigs (extras will be padded with 0s) */

/*   /\* If an invariant subspace was found on the last iteration, W */
/*      contains A*x, where x is a random vector.  Need to orthogonalize */
/*      W against the first iIter vectors in VS.  *\/ */
/*   if (EigSolverData->InvariantSubspaceFound){ */
    
/*     ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm)); */
/*     if (ierr != xf_OK) return ierr; */

/*     if (debug) */
/*       xf_printf("Back after Invariant subspace found.  Reorthogonalizing. iIter = %d\n", iIter); */
/*     ierr = xf_Error(xf_OrthogonalizeVector(VS, iIter, debug, W, NULL)); */
/*     if (ierr != xf_OK) return ierr; */

/*     ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm1)); */
/*     if (ierr != xf_OK) return ierr; */

/*     if (rnorm1 < 0.717*rnorm){ // possibly no other directions, try reortho again */
/*       rnorm = rnorm1; */
/*       ierr = xf_Error(xf_OrthogonalizeVector(VS, iIter, debug, W, NULL)); */
/*       if (ierr != xf_OK) return ierr; */
      
/*       ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm1)); */
/*       if (ierr != xf_OK) return ierr; */
/*       if ((rnorm1 < 0.717*rnorm) || (rnorm1 < 10.0*MEPS)){ */
/* 	if (debug) */
/* 	  xf_printf("Forcing exit because rnorm = %.10E rnorm1 = %.10E\n", rnorm, rnorm1); */
/* 	ForceExit = xfe_True; // no other directions to search; exit. */
/*       } */
/*     } */
/*   } */

/*   // wnorm = ||W|| */
/*   ierr = xf_Error(xf_VectorNorm(W, 2, &wnorm)); */
/*   if (ierr != xf_OK) return ierr; */

/*   if (debug) xf_printf("wnorm = %.10E\n", wnorm); */
/*   if (wnorm < MEPS) ForceExit = xfe_True; */

/*   if (ForceExit){ */
/*     /\* A forced exit occurs when an invariant subspace of eigenvectors */
/*        has been found and when that subspace contains all the nonzero */
/*        eigenvalues.  i.e. all other search directions (operator */
/*        applied to a random vector) result in a vector in the span of */
/*        the invariant subspace. *\/ */
/*     if (debug)  */
/*       xf_printf("Forced exit at iteration iIter = %d\n", iIter); */
/*     ncv = min(nev,iIter); // so we can calculate the eigs of H with the same code */
/*     nev = ncv; // # of nonzero eigenvalues returned */
/*     rnorm = 0.0; // residual norm is zero */
/*   } */
/*   else{ */
/*     /\* On first iteration, have only one vector, so ask for another */
/*        mult.  If an invariant subspace was found, we're also basically */
/*        starting over, so we need another mult. *\/ */
/*     if ((iIter == 0) || (EigSolverData->InvariantSubspaceFound)){ */
/*       beta[iIter] = 0.0; */
/*       EigSolverData->InvariantSubspaceFound = xfe_False; */

/*       // W = W / wnorm */
/*       ierr = xf_Error(xf_VectorMult(W, 1.0/wnorm)); */
/*       if (ierr != xf_OK) return ierr; */

/*       (*Status) = xfe_EigMultiply; */
/*       (*pV) = VS->Vector + iIter; */
/*       (*pW) = VS->Vector + iIter+1; */
/*       return xf_OK; */
/*     } */

/*     // at this point, iIter >= 1, V = VS(iIter-1), W = VS(iIter) */
  
/*     // alpha(iIter-1) = V^T * W */
/*     ierr = xf_Error(xf_VectorDot(V, W, alpha+iIter-1)); */
/*     if (ierr != xf_OK) return ierr; */

/*     // W = W - beta(iIter-1)*VS(iIter-2) - alpha(iIter-1)*VS(iIter-1) */
/*     if (iIter > 1){ */
/*       ierr = xf_Error(xf_VectorMultSet(VS->Vector+iIter-2, beta[iIter-1], xfe_Sub, W)); */
/*       if (ierr != xf_OK) return ierr; */
/*     } */
/*     ierr = xf_Error(xf_VectorMultSet(V, alpha[iIter-1], xfe_Sub, W)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     // rnorm = ||W|| */
/*     ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm)); */
/*     if (ierr != xf_OK) return ierr; */

/*     if (debug) */
/*       xf_printf("iIter = %d; wnorm = %.10E, rnorm = %.15E\n", iIter, wnorm, rnorm); */
  
/*     // if rnorm < 0.717*wnorm, reorthogonalize (following ARPACK) */
/*     rnorm1 = rnorm; */
/*     rnorm  = wnorm; */
/*     for (iOrtho=0; iOrtho<2; iOrtho++){ // reorthogonalize up to two times */
/*       if (rnorm1 < 0.717*rnorm){  */
/* 	if (debug){ */
/* 	  xf_printf("Reorthogonalizing: iOrtho = %d\n", iOrtho); */
/* 	  xf_printf("rnorm = %.10E, rnorm1 = %.10E\n", rnorm, rnorm1); */
/* 	} */

/* 	// Call orthogonalization function */
/* 	ierr = xf_Error(xf_OrthogonalizeVector(VS, iIter, debug, W, &dalpha)); */
/* 	if (ierr != xf_OK) return ierr; */

/* 	alpha[iIter-1] += dalpha; */
      
/* 	rnorm = rnorm1; */
      
/* 	// rnorm1 = ||W|| */
/* 	ierr = xf_Error(xf_VectorNorm(W, 2, &rnorm1)); */
/* 	if (ierr != xf_OK) return ierr; */
/*       } */
/*     } // iOrtho */

/*     /\* If two attempts at reorthogonalization did not help, an invariant */
/*        subspace has been found.  Following ARPACK, generate a new random */
/*        vector, A*rand; make sure to orthogonalize it  */
/*     *\/ */
/*     if (debug)  */
/*       xf_printf("rnorm = %.10E, rnorm1 = %.10E\n", rnorm, rnorm1); */
/*     if (((rnorm1 < 0.717*rnorm) || (rnorm1 < 10.0*MEPS)) && (iIter < ncv)){ */
    
/*       if (debug) */
/* 	xf_printf("Invariant subspace found. iIter = %d\n", iIter); */
    
/*       beta[iIter] = 0.0; */

/*       EigSolverData->InvariantSubspaceFound = xfe_True; */
/*       EigSolverData->iIter--; // so that iteration number stays the same */
/*       if ((EigSolverData->nInvariant++) > nRestartMax){ */
/* 	xf_printf("Lanczos stalled after trying to get away from an invariant subspace.\n"); */
/* 	xf_printf(" %d attempts were made to try to find a new direction.\n", nRestartMax); */
/* 	return xf_Error(xf_NOT_CONVERGED); */
/*       } */
    
/*       // generate a random vector in VS->Vector[ncv] */
/*       ierr = xf_Error(xf_VectorRand(VS->Vector + ncv, 0)); */
/*       if (ierr != xf_OK) return ierr; */
    
/*       // Set pV and pW */
/*       (*Status) = xfe_EigMultiply; */
/*       (*pV) = VS->Vector + ncv; */
/*       (*pW) = VS->Vector + iIter; */
    
/*       return xf_OK; */
/*     } */

/*     if (iIter < ncv){ */
/*       beta[iIter] = rnorm; // set beta(iIter) */
    
/*       // W = W / rnorm */
/*       ierr = xf_Error(xf_VectorMult(W, 1.0/rnorm)); */
/*       if (ierr != xf_OK) return ierr; */
/*     } */


/*     /\* If don't have enough vectors, ask for another mult and return *\/ */
/*     if (iIter < ncv){ */
/*       (*Status) = xfe_EigMultiply; */
/*       (*pV) = VS->Vector + iIter; */
/*       (*pW) = VS->Vector + iIter+1; */
/*       return xf_OK; */
/*     } */
/*   } // end else if no ForceExit */


/*   /\* Allocate memory for ritz value calculation *\/ */
/*   ierr = xf_Error(xf_Alloc( (void **) &Z, ncv*ncv, sizeof(real))); */
/*   if (ierr != xf_OK) return ierr; */

/*   ierr = xf_Error(xf_Alloc( (void **) &workl, ncv, sizeof(real))); */
/*   if (ierr != xf_OK) return ierr; */

/*   ierr = xf_Error(xf_Alloc( (void **) &bounds, ncv, sizeof(real))); */
/*   if (ierr != xf_OK) return ierr; */

/*   ierr = xf_Error(xf_Alloc( (void **) &P, ncv, sizeof(int))); */
/*   if (ierr != xf_OK) return ierr; // permutation vector for sorting */
/*   for (k=0; k<ncv; k++) P[k] = k; */

/*   /\* Copy H matrix into ritz (main diagonal) and workl (sub-diagonal) *\/ */
/*   for (k=0; k<ncv; k++)   ritz[k]  = alpha[k]; */
/*   for (k=0; k<ncv-1; k++) workl[k] = beta[k+1]; */

/*   if (debug){ */
/*     xf_printf("Before call to EigSymTriDiag:\n"); */
/*     for (k=0; k<ncv; k++){ */
/*       xf_printf("alpha[%d] = %.15E, beta[%d] = %.15E\n",  */
/* 		k, alpha[k], k, beta[k]);  */
/*     } */
/*   } */


/*   /\* Compute the eigenvalues and eigenvectors for the curent symmetric */
/*      tridiagonal matrix, H.  Note, this is where MathLapack is */
/*      required. *\/ */
/*   ierr = xf_Error(xf_EigSymTriDiag(ncv, ritz, workl, Z)); */
/*   if (ierr != xf_OK) return ierr; */

/*   /\* Compute the error bounds for the Ritz values (eigs of H).  These */
/*      are just the last entries in the eigenvectors, scaled by the */
/*      residual norm, rnorm.  Note, the eigenvectors in Z are stored in */
/*      column-order (Fortran style). *\/ */
/*   for (k=0; k<ncv; k++) */
/*     bounds[k] = fabs(Z[ncv*k+ncv-1]*rnorm); */

/*   if (debug) */
/*     for (k=0; k<ncv; k++){ */
/*       xf_printf("ritz[%d] = %.10E, bounds[%d] = %.10E\n",  */
/* 		k, ritz[k], k, bounds[k]);  */
/*     } */
  
 
/*   /\* Select wanted Ritz values (e.g. largest in magnitude). The */
/*      unwanted Ritz values (i.e. shifts) are placed into the first np = */
/*      (ncv-nev) locations in ritz, while the wanted values are placed */
/*      into locations [np+1:ncv]. The bounds get sorted as well. *\/ */
/*   np = ncv - nev; */
/*   xf_SortEigenvalues(nev, ncv, ritz, bounds, P); */
  
/*   if (debug){ */
/*     xf_printf("After Sorting:\n"); */
/*     for (k=0; k<ncv; k++){ */
/*       xf_printf("ritz[%d] = %.15E, bounds[%d] = %.15E\n",  */
/* 		k, ritz[k], k, bounds[k]);  */
/*     } */
/*   } */

/*   /\* Check convergence of wanted Ritz values *\/ */
/*   nconv = 0; // # converged values */
/*   for (i=ncv-nev; i<ncv; i++) */
/*     nconv += (bounds[i] <= tol*max(MEPS, fabs(ritz[i]))); */

/*   if (Verbosity >= xfe_VerbosityMedium){ */
/*     xf_printf("[iRestart = %d] Lanczos eigenvalue estimates and bounds:\n", */
/* 	      EigSolverData->nRestart); */
/*     for (i=ncv-nev; i<ncv; i++){ */
/*       xf_printf(" [%d] Eigenvalue = %.10E, Bound = %.10E\n",  */
/* 		i-ncv+nev, ritz[i], bounds[i]); */
/*     } */
/*   } */
  
/*   /\* Check if user requested a halt: stop computation then *\/ */
/*   if (xf_CheckUserHalt("EIGSTOP")){ */
/*     nconv = ncv; */
/*     ForceExit = xfe_True; */
/*   } */


/*   /\* Exit successfully if have enough converged values.  Converged */
/*      eigs/bounds are placed in ritz(0:nconv-1)/bounds(0:nconv-1). */
/*      Eigenvectors are computed if desired. *\/ */
/*   if ((nconv >= nev) || (ForceExit)){ */
/*     if (ForceExit) nconv = ncv; // should be true anyway since rnorm == 0.0 */
/*     for (i=0; i<nconv; i++){ */
/*       swap(ritz[i], ritz[ncv-nev+i], val); */
/*       swap(bounds[i], bounds[ncv-nev+i], val); */
/*       swap(P[i], P[ncv-nev+i], j); */
/*     } */
/*     for (i=nconv; i<nev0; i++) E[i] = 0.0; */

/*     if (EV != NULL){ // Computation of eigenvectors */
/*       for (i=0; i<nconv; i++){ */
/* 	ierr = xf_Error(xf_SetZeroVector(EV->Vector + i)); */
/* 	if (ierr != xf_OK) return ierr; */
/* 	for (j=0; j<ncv; j++){ */
/* 	  ierr = xf_Error(xf_VectorMultSet(VS->Vector + j, Z[P[i]*ncv + j],  */
/* 					   xfe_Add, EV->Vector+i)); */
/* 	  if (ierr != xf_OK) return ierr; */
/* 	} */
/*       } */
/*       for (i=nconv; i<nev0; i++){ // zero out extra eigenvectors if ForceExit */
/* 	ierr = xf_Error(xf_SetZeroVector(EV->Vector + i)); */
/* 	if (ierr != xf_OK) return ierr; */
/*       } */
/*     } */

/*     if ((ForceExit) && (nev < nev0)) */
/*       (*Status) = xfe_EigForceExit; */
/*     else */
/*       (*Status) = xfe_EigConverged; */

/*     ierr = xf_Error(xf_DestroyEigSolverData(EigSolverData)); */
/*     if (ierr != xf_OK) return ierr; */
/*   } */
/*   else{ */

/*     /\* Error out if applied too many restarts *\/  */
/*     if (EigSolverData->nRestart++ > nRestartMax){ */
/*       xf_printf("Lanczos not converged after %d restarts. Giving up.\n", nRestartMax); */
/*       return xf_Error(xf_NOT_CONVERGED); */
/*     } */

    
/*     /\* Count # of unwated Ritz values with zero bounds */
/*        estimates. Decrease np, the number of shifts to apply. *\/ */
/*     for (i=0, nptemp = np; i<nptemp; i++) */
/*       if (bounds[i] <= MEPS) np--; */
/*       //if (bounds[i] == 0.) np--; */
    
/*     if (nptemp != np){ */
/*       nev += nptemp-np; */
/*       xf_printf("Decreased # shifts due to zero bounds estimates. np = %d, nev = %d\n",  */
/* 		np, nev); */
/*     } */

/*     /\* Adjust nev to prevent stagnation *\/ */
/*     nevbef = nev; */
/*     nev = nev + min (nconv, np/2); */
/*     if ((nev == 1) && (ncv >= 6)) */
/*       nev = ncv / 2; */
/*     else if ((nev == 1) && (ncv > 2)) */
/*       nev = 2; */
/*     np = ncv - nev; */

/*     // if nev has increased, resort eigenvalues */
/*     if (nevbef < nev){ */
/*       if (debug) xf_printf("Resorting eigenvalues: nevbef = %d, nev = %d\n", nevbef, nev); */
/*       xf_SortEigenvalues(nev, ncv, ritz, bounds, P); */
/*       if (debug){ */
/* 	xf_printf("After ReSorting:\n"); */
/* 	for (k=0; k<ncv; k++){ */
/* 	  xf_printf("ritz[%d] = %.15E, bounds[%d] = %.15E\n", */
/* 		    k, ritz[k], k, bounds[k]); */
/* 	} */
/*       } */
/*     } */
    
/*     if (debug){ */
/*       xf_printf("Before shifts:\n"); */
/*       for (k=0; k<ncv; k++){ */
/* 	xf_printf("alpha[%d] = %.15E, beta[%d] = %.15E\n",  */
/* 		  k, alpha[k], k, beta[k]);  */
/*       } */
/*     } */
    
/*     /\* Apply shifts -- separate function to do this *\/ */
/*     ierr = xf_Error(xf_ApplyShifts(VS, nev, np, ritz, alpha, beta)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     if (debug){ */
/*       xf_printf("After shifts:\n"); */
/*       for (k=0; k<ncv; k++){ */
/* 	xf_printf("alpha[%d] = %.15E, beta[%d] = %.15E\n",  */
/* 		  k, alpha[k], k, beta[k]);  */
/*       } */
/*     } */
      
/*     if (debug) */
/*       for (j=0; j<iIter-1; j++){ */
/* 	ierr = xf_Error(xf_VectorDot(VS->Vector+j, VS->Vector+j+1, &val)); */
/* 	if (ierr != xf_OK) return ierr; */
/* 	xf_printf(" dp[%d] = %.10E\n", j, val); */
/*       } */

/*     /\* Set pointers for next call and return with a multiply request *\/ */
/*     (*Status) = xfe_EigMultiply; */
/*     EigSolverData->iIter = nev; */
/*     EigSolverData->nInvariant = 0; // reset counter upon restart */
/*     (*pV) = VS->Vector + nev-1; */
/*     (*pW) = VS->Vector + nev; */
/*   } */

/*   /\* Release memory *\/ */
/*   xf_Release( (void *) bounds); */
/*   xf_Release( (void *) workl); */
/*   xf_Release( (void *) Z);  */
/*   xf_Release( (void *) P); */

/*   return xf_OK; */

/* } */





#if( UNIT_TEST==1 )
#include "xf_EigSolver.test.in"
#endif
