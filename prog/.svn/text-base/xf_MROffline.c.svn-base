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
  FILE:  xf_MROffline.c

  This program performs the offline portion of model reduction.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_MeshTools.h"
#include "xf_SolverTools.h"
#include "xf_Param.h"
#include "xf_Basis.h"
#include "xf_EqnSetHook.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_Residual.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Quad.h"
#include "xf_Output.h"
#include "xf_Arg.h"
#include "xf_MRStruct.h"
#include "xf_MRCommon.h"
#include "xf_ParamDefault.h"



/******************************************************************/
//   FUNCTION Definition: xf_CalculateNonlinearSnap
static int 
xf_CalculateNonlinearSnap(xf_All *All, const xf_VectorSet *SnapSet, 
			  xf_EqnSet **EqnSets, int iResTerm, 
			  int OrderIncrement, int **pQuadOrder, 
			  xf_VectorSet **pSSnapSet)
{
  int ierr, nSnap, iSnap, sr, i, dim, nq;
  int egrp, elem, negrp, *QuadOrder, pnq;
  int *IParam, *rvec;
  enum xfe_Bool QuadChanged;
  char resValue[80];
  real *RParam, *u, *s, *xq, *EU, *ES;
  xf_Data *D;
  xf_Vector *U;
  xf_VectorSet *SSnapSet;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;
  xf_ResTerms *ResTerms;
  xf_ResTerm *ResTerm;

  Mesh = All->Mesh;
  sr   = All->EqnSet->StateRank;
  dim  = All->EqnSet->Dim;


  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(All, All->EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;

  // pull of prototypical (first) snapshot vector
  U = SnapSet->Vector+0; // all vectors in set are same basis/order

  // determine number of element groups
  if ((negrp = Mesh->nElemGroup) != SnapSet->Vector[0].nArraySelf)
    return xf_Error(xf_INCOMPATIBLE);


  // allocate vector of quad orders and ranks
  ierr = xf_Error(xf_Alloc( (void **) pQuadOrder, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  QuadOrder = (*pQuadOrder);
  ierr = xf_Error(xf_Alloc( (void **) &rvec, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // determine quad order and # quad points for each group
  QuadData = NULL;
  for (egrp=0; egrp<negrp; egrp++){
    ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*U->Order[egrp]+OrderIncrement, 
					QuadOrder+egrp));
    if (ierr != xf_OK) return ierr;
    rvec[egrp] = 0;
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      /* Pull off quad points for the element; will not recalculate in generic case */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder[egrp], &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;
      rvec[egrp] = max(rvec[egrp], QuadData->nquad*sr);
    }
    xf_printf("egrp = %d: rvec = %d\n", egrp, rvec[egrp]);
    xf_printf("  U->Order = %d\n", U->Order[egrp]);
    xf_printf("  OrderIncrement = %d\n", OrderIncrement);
  } // egrp
  
  // number of snapshots
  nSnap = SnapSet->nVector;
  if (nSnap <= 0) return xf_Error(xf_INPUT_ERROR); // no snapshots


  /* Allocate a non-interpolated snapshot vectorset.  This vector will
     store the nonlinear snapshots at the desired quad points. */
  ierr = xf_Error(xf_FindVectorSet(All, nSnap, "SSnapSet", xfe_LinkageGlobElem, sr, NULL, 
				   0, 0, NULL, NULL, rvec, xfe_SizeReal, xfe_False, 
				   xfe_False, NULL, pSSnapSet, NULL));
  if (ierr != xf_OK) return ierr;

  // release rvec
  xf_Release( (void *) rvec);
  
  SSnapSet = (*pSSnapSet);

  // loop over elements; calculate source at Lagrange nodes
  PhiData = NULL;
  u = NULL;
  s = NULL;
  pnq = -1;
  for (egrp=0; egrp<negrp; egrp++){
   
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      /* Pull off quad points for the element; will not recalculate in generic case */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder[egrp], &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;

      // compute basis functions for U vector
      ierr = xf_Error(xf_EvalBasis(U->Basis[egrp], U->Order[egrp], QuadChanged, 
				   nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;

      
      // re=allocate space for u and s
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc((void **) &u, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc((void **) &s, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }


      for (iSnap=0; iSnap<nSnap; iSnap++){
	EU = SnapSet->Vector[iSnap].GenArray[egrp].rValue[elem];
      
	// interpolate snapshot at quad points
	xf_MxM_Set(PhiData->Phi, EU, nq, PhiData->nn, sr, u);

	// pull off ResTerm for this snapshot (contains params)
	EqnSet = EqnSets[iSnap];
	ResTerms = EqnSet->ResTerms;
	ResTerm = ResTerms->ResTerm+iResTerm;
	if (ResTerm->Type != xfe_ResTermSource) return xf_Error(xf_NOT_SUPPORTED);

	/* Evaluate nonlinear source term at all points for this
	   snapshot.  Note, use of All->EqnSet means that all variable
	   parameters must be specified in ResTerm itself.  Could
	   change this in the future if want say viscosity to be a
	   parameter. */
	ierr = xf_Error(xf_EqnSetSourceS(All->EqnSet, ResTerm, 1, IParam, RParam, nq, 
					 u, NULL, NULL, NULL, NULL, s, NULL, NULL, NULL));
	if (ierr != xf_OK) return ierr;

	ES = SSnapSet->Vector[iSnap].GenArray[egrp].rValue[elem];
	for (i=0; i<nq*sr; i++)  ES[i] = s[i];
      }

      pnq = nq;
    } // elem

  } // egrp

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
      
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) u);
  xf_Release( (void *) s);


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Snapshot2BasisPODOld
static int 
xf_Snapshot2BasisPODOld(const xf_VectorSet *SnapSet, int nSnap, int nBasis,
			xf_VectorSet *BasisSet)
{
  int ierr, i, j;
  real dp;
  real *K, *EG, *EV;
  xf_Data *D;

  // allocate memory
  ierr = xf_Error(xf_Alloc( (void **) &K, nSnap*nSnap, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  xf_printf("  Calculating K ... ");
  for (i=0; i<nSnap; i++)
    for (j=i; j<nSnap; j++){
      ierr = xf_Error(xf_VectorDot(SnapSet->Vector+i, SnapSet->Vector+j, &dp));
      if (ierr != xf_OK) return ierr;
      dp = dp/( (real) nSnap);
      K[i*nSnap+j] = dp;
      K[j*nSnap+i] = dp; // K is symmetric
    } // j
  xf_printf("done.\n");
  
  xf_printf("K = [");
  for (i=0; i<nSnap; i++)
    for (j=0, xf_printf("\n"); j<nSnap; j++)
      xf_printf("%.6E ", K[i*nSnap+j]);
  xf_printf("\n]\n");


  // allocate memory for eigenvalues and eigenvectors
  ierr = xf_Error(xf_Alloc( (void **) &EG, nSnap, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc( (void **) &EV, nSnap*nSnap, sizeof(real)));
  if (ierr != xf_OK) return ierr;


  /* Calculate nBasis max eigs and eigvs of K. */
  xf_printf("  Calculating eigenvalues ... "); fflush(stdout);
  ierr = xf_Error(xf_EigSymmetricQR(K, nSnap, 1e-8, 25, EG, EV));
  if (ierr == xf_NOT_CONVERGED) xf_printf("Warning, eig not converged.\n");
  else if (ierr != xf_OK) return ierr;
  
  //for (i=0; i<nSnap; i++) xf_printf("\n %.10E", EG[i]);
  xf_printf("done.\n");

  
  /* Construct basis vectors */
  xf_printf("  Constructing Basis vectors ... ");
  for (i=0; i<nBasis; i++){
    ierr = xf_Error(xf_SetZeroVector(BasisSet->Vector+i));
    if (ierr != xf_OK) return ierr;
    for (j=0; j<nSnap; j++){
      ierr = xf_Error(xf_VectorMultSet(SnapSet->Vector+j, EV[i*nSnap+j], 
				       xfe_Add, BasisSet->Vector+i));
      if (ierr != xf_OK) return ierr;
    }
  }
  xf_printf("done.\n");

  /* Release memory */
  xf_Release( (void *) K);
  xf_Release( (void *) EG);
  xf_Release( (void *) EV);


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DumpVectorSet
static int 
xf_DumpVectorSet(xf_VectorSet *VS, int n, const char *fname)
{ 
  int ierr, egrp, elem, j, r;
  real **EV;
  xf_Vector *V;
  FILE *fid;

  
  // open file
  if ((fid = fopen(fname, "w")) == NULL) return xf_Error(xf_FILE_WRITE_ERROR);
  
  ierr = xf_Error(xf_Alloc( (void **) &EV, n, sizeof(real *)));
  if (ierr != xf_OK) return ierr;

  V = VS->Vector + 0;
  for (egrp=0; egrp<V->nArraySelf; egrp++){
    for (elem=0; elem<V->GenArray[egrp].n; elem++){
      for (j=0; j<n; j++) EV[j] = VS->Vector[j].GenArray[egrp].rValue[elem];
      for (r=0; r<V->GenArray[egrp].r; r++){
	for (j=0; j<n; j++)
	  fprintf(fid, "%.15E ", EV[j][r]);
	fprintf(fid, "\n");
      }
    }
  }
  fclose(fid);

  xf_Release( (void *) EV);
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_VectorMaxLoc
static int 
xf_VectorMaxLoc(xf_Mesh *Mesh, const int *QuadOrder, const xf_Vector *V, 
		xf_IPointType *z, int l)
{
  int ierr, i, j, k, r, dim;
  int myRank, nProc, iProc;
  int proc, egrp, elem, node;
  int procmax;
  int ibuf[3], *ibuf0;
  enum xfe_Bool QuadChanged;
  real val, valmax, rbuf[4], xref[3], *rbuf0, *xn;
  xf_QuadData *QuadData;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  // calculate maximum value on proc
  valmax = -1.0;
  for (i=0; i<V->nArraySelf; i++)
    for (j=0; j<V->GenArray[i].n; j++)
      for (r=0; r<V->GenArray[i].r; r++)
	if ((val = fabs(V->GenArray[i].rValue[j][r])) > valmax){
	  egrp = i;
	  elem = j;
	  node = r;
	  valmax = val;
	}

  if (valmax < 0) return xf_Error(xf_CODE_LOGIC_ERROR);

  // get dim
  dim = Mesh->Dim;
  
  // get coordinates of quad point node on egrp,elem
  QuadData = NULL;
  ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder[egrp], 
			      &QuadData, &QuadChanged));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<dim; k++) xref[k] = QuadData->xquad[node*dim+k];
  // Destroy QuadData (if points are generic)
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;



  // collect data into send buffers
  for (k=0; k<dim; k++) rbuf[k] = xref[k];
  rbuf[dim] = valmax;
  ibuf[0] = egrp; ibuf[1] = elem; ibuf[2] = node;

  // allocate memory for receive buffers on root
  rbuf0 = NULL;
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc( (void **) &rbuf0, nProc*(dim+1), sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  ibuf0 = NULL;
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc( (void **) &ibuf0, nProc*3, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }

  // gather data on root
  if (nProc > 1){
    ierr = xf_Error(xf_MPI_Gather( (void *) rbuf, (void *) rbuf0, (dim+1)*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_MPI_Gather( (void *) ibuf, (void *) ibuf0, 3*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
  }
  else{
    for (k=0; k<(dim+1); k++) rbuf0[k] = rbuf[k];
    for (k=0; k<3; k++) ibuf0[k] = ibuf[k];
  }

  // root picks out max value, fills in z
  if (myRank == 0){
    valmax  = -1.0;
    procmax = -1;
    for (iProc=0; iProc<nProc; iProc++){
      if ((val=rbuf0[(dim+1)*iProc+dim]) > valmax){
	valmax = val;
	procmax = iProc;
      }
    }
    if (procmax < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    z->proc[l] = procmax;
    z->egrp[l] = ibuf0[3*procmax+0];
    z->elem[l] = ibuf0[3*procmax+1];
    z->node[l] = ibuf0[3*procmax+2];
    for (k=0; k<dim; k++) 
      z->xref[dim*l+k] = rbuf0[(dim+1)*procmax+k];
  }  

  xf_Release( (void *) rbuf0);
  xf_Release( (void *) ibuf0);

  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateVectorSetAtPoint
static int 
xf_CalculateVectorSetAtPoint(const xf_VectorSet *BasisSet, int M,
			     xf_IPointType *z, int l, real *PA)
{
  int ierr, j;
  int myRank, nProc;
  int proc, egrp, elem, node;
  int ibuf[4];

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;


  // Broadcast proc, egrp, elem, node of z{l}
  if (myRank == 0){
    ibuf[0] = z->proc[l]; ibuf[1] = z->egrp[l];
    ibuf[2] = z->elem[l]; ibuf[3] = z->node[l];
  }    
  if (nProc > 1){
    ierr = xf_Error(xf_MPI_Bcast( (void *) ibuf, 4*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
  }
  proc = ibuf[0]; egrp = ibuf[1]; elem = ibuf[2]; node = ibuf[3];
    
  // calculate basis functions at z{l}; store in lth row of PA
  if (myRank == proc){      
    for (j=0; j<M; j++)
      PA[j] = BasisSet->Vector[j].GenArray[egrp].rValue[elem][node];
  }

  // send row of PA to root if in parallel
  if (proc != 0){	
    if (myRank == proc){
      // send PA row
      ierr = xf_Error(xf_MPI_Send( (void *) PA, M*sizeof(real), 0, 0));
      if (ierr != xf_OK) return ierr;
    }
    if (myRank == 0){
      // receive PA row
      ierr = xf_Error(xf_MPI_Recv( (void *) PA, M*sizeof(real), proc, 0));
      if (ierr != xf_OK) return ierr;
    }
  }

  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CoeffExpansionEIM
static int 
xf_CoeffExpansionEIM(xf_All *All, const xf_VectorSet *BasisSet, int M, 
		     const int *QuadOrder, xf_IPointType *z, 
		     xf_VectorSet *CardSet)
{
  int ierr, i, j, k, l, dim;
  int myRank, nProc;
  int *p;
  enum xfe_Bool ParallelFlag;
  real *xn, *PA, *P, *b, *sigma;
  xf_Vector *Phi, *PhiTemp;
  xf_Mesh *Mesh;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  ParallelFlag = (nProc > 1);

  Mesh = All->Mesh;
  dim = Mesh->Dim;

  // allocate interpolation point structure on root
  if (myRank == 0){
    ierr = xf_Error(xf_AllocIPoint( z, M, dim));
    if (ierr != xf_OK) return ierr;
  }

  /*
    Each interpolation point is stored as an element index and the
    reference coordinates in the element.

    EIM algorithm:  

    (BasisSet = Phi{i}, CardSet = Psi{i})
    (The Phi are assumed to be stored using a nodal basis)

    z{0} = arg max |Phi{0}| 
    for l = 1:M-1,
      L = [0:l-1]
      P = Phi{L}(z{L})    % rows correspond to the z{L} 
      b = Phi{l}(z{L})
      sigma = inv(P)*b
      z{l} = arg max |Phi{l}(x) - Phi{L}(x)*sigma{L}^T|
    end

    Psi{j}(x) = Phi{k}(x) * P{k,j}
    Psi{j}(z{i}) = delta{i,j} -> P{k,j} = inv( Phi{:}(z{:}) )

  */

  /* Find maximum value and location of |Phi{0}| */
  Phi = BasisSet->Vector + 0;
  if (Phi->StateRank != 1) return xf_Error(xf_NOT_SUPPORTED); // for now
  ierr = xf_Error(xf_VectorMaxLoc(Mesh, QuadOrder, Phi, z, 0));
  if (ierr != xf_OK) return ierr;
  
  // allocate memory for PA = all Phi evaluated at all interpolation points
  ierr = xf_Error(xf_Alloc( (void **) &PA, M*M, sizeof(real))); 
  if (ierr != xf_OK) return ierr;

  // allocate memory for P matrix
  ierr = xf_Error(xf_Alloc( (void **) &P, M*M, sizeof(real))); 
  if (ierr != xf_OK) return ierr;

  // allocate memory for permutation vector, p
  ierr = xf_Error(xf_Alloc( (void **) &p, M, sizeof(int))); 
  if (ierr != xf_OK) return ierr;

  // allocate memory for rhs vector
  ierr = xf_Error(xf_Alloc( (void **) &b, M, sizeof(real))); 
  if (ierr != xf_OK) return ierr;

  // allocate memory for solution vector, sigma
  ierr = xf_Error(xf_Alloc( (void **) &sigma, M, sizeof(real))); 
  if (ierr != xf_OK) return ierr;

  // allocate memory for a temporary vector, PhiTemp
  ierr = xf_Error(xf_FindSimilarVector(All, Phi, "PhiTemp", xfe_False, 
				       xfe_False, NULL, &PhiTemp, NULL));
  if (ierr != xf_OK) return ierr;

  // calculate basis functions at z0; store in first row of PA  
  ierr = xf_Error(xf_CalculateVectorSetAtPoint(BasisSet, M, z, 0, PA));
  if (ierr != xf_OK) return ierr;
  
  /* Loop from 1 to M-1 */
  for (l=1; l<M; l++){
    
    if (myRank == 0){
      // construct P matrix P = Phi{0:l-1}(z{0:l-1})
      for (i=0; i<l; i++)
	for (j=0; j<l; j++)
	  P[i*l+j] = PA[i*M+j];
	
      // construct rhs b = Phi{l}(z{0:l-1})
      for (i=0; i<l; i++)
	b[i] = PA[i*M+l];
      
      // PLU factor P
      ierr = xf_Error(xf_ComputePLU(P, l, p));
      if (ierr != xf_OK) return ierr;
      
      // solve for sigma
      ierr = xf_Error(xf_SolvePLU(P, p, b, l, sigma, NULL));
      if (ierr != xf_OK) return ierr;
    }

    // Broadcast solution vector, sigma, if in parallel
    if (ParallelFlag){
      ierr = xf_Error(xf_MPI_Bcast( (void *) sigma, l*sizeof(real), 0));
      if (ierr != xf_OK) return ierr;
    }

    // PhiTemp = Phi{l}(x) - Phi{L}(x)*sigma{L}^T
    ierr = xf_Error(xf_VectorMultSet(BasisSet->Vector+l, 1.0, xfe_Set, PhiTemp));
    if (ierr != xf_OK) return ierr;
    for (j=0; j<l; j++){
      ierr = xf_Error(xf_VectorMultSet(BasisSet->Vector+j, -sigma[j], xfe_Add, PhiTemp));
      if (ierr != xf_OK) return ierr;
    }
    
    // z{l} = arg max |PhiTemp|
    ierr = xf_Error(xf_VectorMaxLoc(Mesh, QuadOrder, PhiTemp, z, l));
    if (ierr != xf_OK) return ierr;

    // calculate basis functions at z{l}; store in lth row of PA  
    ierr = xf_Error(xf_CalculateVectorSetAtPoint(BasisSet, M, z, l, PA+l*M));
    if (ierr != xf_OK) return ierr;

  } // l


  if (myRank == 0){
    // P = inv(PA)
    ierr = xf_Error(xf_ComputePLU(PA, M, p));
    if (ierr != xf_OK) return ierr;
    for (k=0; k<M*M; k++     ) P[k] = 0.0;
    for (k=0; k<M*M; k+=(M+1)) P[k] = 1.0;
    ierr = xf_Error(xf_SolvePLU_Matrix(PA, p, M, M, P));
    if (ierr != xf_OK) return ierr;
  }

  // Broadcast solution matrix, P, if in parallel
  if (ParallelFlag){
    ierr = xf_Error(xf_MPI_Bcast( (void *) P, M*M*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
  }

  // Psi{j}(x) = Phi{k}(x) * P{k,j}
  for (j=0; j<M; j++){
    ierr = xf_Error(xf_SetZeroVector(CardSet->Vector+j));
    if (ierr != xf_OK) return ierr;
    for (k=0; k<M; k++){
      ierr = xf_Error(xf_VectorMultSet(BasisSet->Vector+k, P[k*M+j], xfe_Add, CardSet->Vector+j));
      if (ierr != xf_OK) return ierr;
    }
  }

  /*  // test whether CardSet is cardinal */
  /*   xf_printf("Cardinality test: \n"); */
  /*   for (l=0; l<M; l++){ */
  /*     xf_printf("\n"); */
  /*     for (j=0; j<M; j++) */
  /*       xf_printf("%.6E ", CardSet->Vector[j].GenArray[z->egrp[l]].rValue[z->elem[l]][z->node[l]]); */
  /*   } */
  /*   xf_printf("\n"); */

  xf_Release( (void *) PA);
  xf_Release( (void *) P);
  xf_Release( (void *) b);
  xf_Release( (void *) p);
  xf_Release( (void *) sigma);

  ierr = xf_Error(xf_DestroyVector(PhiTemp, xfe_True));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ScaleVectorSetByWeights
static int 
xf_ScaleVectorSetByWeights(xf_All *All, xf_VectorSet *VS, int nVector,
			   const int *QuadOrder, enum xfe_Bool InverseFlag)
{
  int ierr, sr, k, iq, nq, pnq;
  int egrp, elem, i;
  enum xfe_Bool QuadChanged;
  real *EU, *xq, *wq, fac;
  xf_QuadData *QuadData;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;

  if (nVector <= 0) return xf_Error(xf_INPUT_ERROR);

  Mesh = All->Mesh;
  sr   = VS->Vector[0].StateRank;

  if (sr != 1) return xf_Error(xf_NOT_SUPPORTED); // for now

  QuadData = NULL;
  JData    = NULL;
  wq       = NULL;
  pnq = -1;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      /* Pull off quad points for the element; will not recalculate in generic case */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder[egrp], &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;

      // element Jacobian
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;

      // re-allocate memory if quad points increased
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      // scale U
      for (i=0; i<nVector; i++){
	EU = VS->Vector[i].GenArray[egrp].rValue[elem];
	if (VS->Vector[i].GenArray[egrp].r != nq) return xf_Error(xf_OUT_OF_BOUNDS);
	for (iq=0; iq<nq; iq++){
	  fac = ((InverseFlag) ? 1.0/sqrt(fabs(wq[iq])) : sqrt(fabs(wq[iq])));
	  EU[iq] *= fac;
	}
      }

      pnq = nq;
    } // elem
  } // egrp


  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) wq);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_NonlinearReducedMatrixE
static int 
xf_NonlinearReducedMatrixE(xf_All *All, const xf_VectorSet *UBasisSet,
			   const xf_VectorSet *SBasisSet, int NU, int NS,
			   const int *QuadOrder, real *E)
{
  int ierr, dim, sr, k, iq, nq, pnq;
  int egrp, elem, nu, ns, nn;
  int Order;
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged;
  real *EU, *ES, *xq, *u, *wq, t;
  xf_Vector *U, *S;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;

  if ((NU <=0) || (NS <= 0)) return xf_Error(xf_INPUT_ERROR);

  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  sr   = UBasisSet->Vector[0].StateRank;

  if (sr != 1) return xf_Error(xf_NOT_SUPPORTED); // for now

  // representative vectors
  U = UBasisSet->Vector+0;
  S = SBasisSet->Vector+0;

  // zero out E
  for (k=0; k<NU*NS; k++) E[k] = 0.;

  QuadData = NULL;
  PhiData  = NULL;
  JData    = NULL;
  u        = NULL;
  wq       = NULL;
  pnq = -1;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    Basis = U->Basis[egrp]; Order = U->Order[egrp];

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      /* Pull off quad points for the element; will not recalculate in generic case */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder[egrp], &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;

      // compute basis functions for U vector
      ierr = xf_Error(xf_EvalBasis(U->Basis[egrp], U->Order[egrp], QuadChanged, 
				   nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;

      // element Jacobian
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;

      nn = PhiData->nn; // number of interpolation nodes

      // re-allocate memory if quad points increased
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      // add to entries of E matrix
      for (nu=0; nu<NU; nu++){
	EU = UBasisSet->Vector[nu].GenArray[egrp].rValue[elem];
	xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u); // interpolate u
	for (ns=0; ns<NS; ns++){
	  ES = SBasisSet->Vector[ns].GenArray[egrp].rValue[elem];
	  if (SBasisSet->Vector[ns].GenArray[egrp].r != nq) 
	    return xf_Error(xf_OUT_OF_BOUNDS);
	  for (iq=0, t=0.; iq<nq; iq++) t += wq[iq]*u[iq]*ES[iq];
	  E[nu*NS+ns] += 1.0*t;
	}
      }

      pnq = nq;
    } // elem
  } // egrp

  // reduce-sum E
  ierr = xf_Error(xf_MPI_Allreduce(E, NU*NS, xfe_SizeReal, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

  // Destroy Basis Data
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) u);
  xf_Release( (void *) wq);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_NonlinearReducedMatrixD
static int 
xf_NonlinearReducedMatrixD(const xf_VectorSet *BasisSet, xf_IPointType *z,
			   int N, int M, real *D)
{
  int ierr, m, n, k, dim, sr;
  int egrp, elem, proc;
  int myRank, nProc;
  int ibuf[4];
  real *EV, xref[3] = {0};
  xf_Vector *V;
  xf_BasisData *PhiData;

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if ((N <=0) || (M <= 0)) return xf_Error(xf_INPUT_ERROR);
  sr = BasisSet->Vector[0].StateRank;

  if (sr != 1) return xf_Error(xf_NOT_SUPPORTED); // for now
  
  // Evaluate the N BasisSet functions at the M points z
  PhiData = NULL;
  for (m=0; m<M; m++){ // loop over the M points

    // pull off the element info + reference coordinate
    if (myRank == 0){
      ibuf[0] = z->proc[m];
      ibuf[1] = z->egrp[m];
      ibuf[2] = z->elem[m];
      ibuf[3] = z->Dim;
      for (k=0; k<z->Dim; k++) xref[k] = z->xref[m*z->Dim+k];
    }

    // broadcast data if in parallel
    if (nProc > 1){
      ierr = xf_Error(xf_MPI_Bcast( (void *) ibuf, 4*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_MPI_Bcast( (void *) xref, 3*sizeof(real), 0));
      if (ierr != xf_OK) return ierr;
    }
    proc = ibuf[0]; egrp = ibuf[1]; elem = ibuf[2]; dim = ibuf[3];


    if (myRank == proc){

      // evaluate the FE basis
      V = BasisSet->Vector+0;
      ierr = xf_Error(xf_EvalBasis(V->Basis[egrp], V->Order[egrp], xfe_True, 
				   1, xref, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      // interpolate the BasisSet functions at xref
      for (n=0; n<N; n++){
	EV = BasisSet->Vector[n].GenArray[egrp].rValue[elem];
	xf_MxM_Set(PhiData->Phi, EV, 1, PhiData->nn, sr, D + m*N + n);
      } //n
    }

    // send row of D to root if in parallel
    if (proc != 0){	
      if (myRank == proc){
	// send row
	ierr = xf_Error(xf_MPI_Send( (void *) (D + m*N), N*sizeof(real), 0, 0));
	if (ierr != xf_OK) return ierr;
      }
      if (myRank == 0){
	// receive row
	ierr = xf_Error(xf_MPI_Recv( (void *) (D + m*N), N*sizeof(real), proc, 0));
	if (ierr != xf_OK) return ierr;
      }
    }

  } //m

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ResTermIsLinear
static enum xfe_Bool
xf_ResTermIsLinear(xf_ResTerm *ResTerm)
{
  // hardcode for now; later use eqnset
  return (ResTerm->Type != xfe_ResTermSource);
}

/******************************************************************/
//   FUNCTION Definition: mymain
int 
mymain(int argc, char *argv[])
{
  int len, ierr, terr, i, j;
  int myRank, nProc;
  int nSnap, iSnap, nBasis, nSBasis, iResTerm;
  int iNonLinear, nNonLinear, nOutput, N, M;
  int OrderIncrement, ibuf[6], *QuadOrder;
  enum xfe_Bool ParallelFlag, DebugFlag = xfe_False, UNFlag = xfe_False;
  enum xfe_Bool AsciiFlag = xfe_False;
  char *ArgIn[] = {"xfa", "NULL", ".xfa file name to read (contains mesh)",
		   "eqn", "NULL", "root name in [eqn]_#.eqn for each snapshot",
		   "data", "NULL", "root name for snapshot data files",
		   "K", "10", "number of snapshots to load",
		   "N", "5", "number of Basis (e.g. POD) vectors to use",
		   "M", "5", "number of nonlinear basis vectors to use",
		   "pplus", "0", "order increment for nonlinear snapshots (not used)",
		   "UN", "1", "set to >0 to calculate UN info (see xf_MRStruct.h)",
		   "ascii", "1", "set to >0 to write ROM file in ascii (.m) format",
		   "debug", "0", "set to >0 to print out debug info",
		   "\0"};
  char xfaFile[xf_MAXSTRLEN]; 
  char EqnSetRoot[xf_MAXSTRLEN]; 
  char EqnSetFile[xf_MAXSTRLEN]; 
  char DataRoot[xf_MAXSTRLEN];
  char DataFile[xf_MAXSTRLEN];
  char OutputFile[xf_MAXSTRLEN];
  FILE *fid;
  real dp, val, *Emat, *Dmat, xglob[3];
  xf_KeyValue KeyValue;
  xf_DataSet *DataSet;
  xf_Data *StateData, *D;
  xf_Vector *U, *R0;
  xf_VectorSet  *SnapSet,  *BasisSet, *Residuals;
  xf_VectorSet *SSnapSet, *SBasisSet, *SCardSet;
  xf_ReducedModel *RM;
  xf_IPointType *z;
  xf_ResTerm *ResTerm;
  xf_ResTerms *ResTerms;
  xf_EqnSet *EqnSet;
  xf_EqnSet **EqnSets; // vector of eqnsets for all snapshots
  xf_All *All;

 

  /* Determine myRank*/
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  if (nProc > 1) ParallelFlag = xfe_True;
  else  ParallelFlag = xfe_False;


  xf_printf("\n");
  xf_printf("=== Model Reduction: Offline Computation ===\n");
  xf_printf("\n");
  

  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  /* Parse arguments on root */
  if (myRank == 0)
    terr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);

  ierr = xf_Error(xf_MPI_Bcast((void *) &terr, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (terr == xf_FORCE_QUIT) return xf_OK;
  if (terr != xf_OK) return xf_Error(terr);

  /* Parallelize KeyValue so every proc gets a copy */
  ierr = xf_Error(xf_ParallelizeKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;


  xf_printf("nKey = %d\n", KeyValue.nKey);
  for (i=0; i<KeyValue.nKey; i++)
    xf_printf("%d : Key = %s, Value = %s\n", i, KeyValue.Key[i], KeyValue.Value[i]);
    
  /* Get DataRoot */
  ierr = xf_GetKeyValue(KeyValue, "data", DataRoot);
  if (ierr != xf_OK) return ierr;

  /* Get number of snapshots */
  ierr = xf_GetKeyValueInt(KeyValue, "K", &nSnap);
  if (ierr != xf_OK) return ierr;
  if (nSnap <= 0) return xf_Error(xf_INPUT_ERROR); // nothing to work with

  /* Get number of pod basis vectors */
  ierr = xf_GetKeyValueInt(KeyValue, "N", &nBasis);
  if (ierr != xf_OK) return ierr;
  if ((nBasis <= 0) || (nBasis > nSnap)) return xf_Error(xf_INPUT_ERROR);

  /* Get number of nonlinear basis vectors */
  ierr = xf_GetKeyValueInt(KeyValue, "M", &nSBasis);
  if (ierr != xf_OK) return ierr;
  if ((nBasis <= 0) || (nBasis > nSnap)) return xf_Error(xf_INPUT_ERROR);

  /* Get order increment for nonlinear term */
  ierr = xf_GetKeyValueInt(KeyValue, "pplus", &OrderIncrement);
  if (ierr != xf_OK) return ierr;
  if (OrderIncrement < 0) return xf_Error(xf_INPUT_ERROR);
    
  
  /* Get DebugFlag */
  ierr = xf_GetKeyValueInt(KeyValue, "debug", (int *) &DebugFlag);
  if (ierr != xf_OK) return ierr;
  
  /* Get UNFlag */
  ierr = xf_GetKeyValueInt(KeyValue, "UN", (int *) &UNFlag);
  if (ierr != xf_OK) return ierr;
  
  /* Get AsciiFlag */
  ierr = xf_GetKeyValueInt(KeyValue, "ascii", (int *) &AsciiFlag);
  if (ierr != xf_OK) return ierr;

  /* Get .xfa name */
  ierr = xf_GetKeyValue(KeyValue, "xfa", xfaFile);
  if (ierr != xf_OK) return ierr;
  len = strlen(xfaFile);
  if ((len < 4) || (strncmp(xfaFile+len-4, ".xfa", 4) != 0)){
    xf_printf("Error, xfaFile requires .xfa extension.\n");
    return xf_Error(xf_INPUT_ERROR);
  }

  
  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  /* Read .xfa file*/
  ierr = xf_Error(xf_ReadAllBinary(xfaFile, All));
  if (ierr!=xf_OK) return ierr;
  
  /* Set default parameters: do not overwrite. */
  ierr = xf_Error(xf_AddKeyValueList(&All->Param->KeyValue, xf_DefaultParamList,
				     xfe_False, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* EqnSetRoot */
  ierr = xf_GetKeyValue(KeyValue, "eqn", EqnSetRoot);
  if (ierr != xf_OK) return ierr;
  
  /* need eqnsets */
  if (!xf_NotNull(EqnSetRoot)) return xf_Error(xf_INPUT_ERROR);
    
  ierr = xf_Error(xf_Alloc((void **) &EqnSets, nSnap, sizeof(xf_EqnSet *)));
  if (ierr != xf_OK) return ierr;
  
  for (iSnap=0; iSnap<nSnap; iSnap++){
    sprintf(EqnSetFile, "%s_%d.eqn", EqnSetRoot, iSnap);
    
    ierr = xf_Error(xf_CreateEqnSet(&EqnSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadEqnSetFile(EqnSetFile, NULL, EqnSet));
    if (ierr != xf_OK) return ierr;
    
    if (iSnap == 0){ // 0th eqnset is put into All for output purposes
      ierr = xf_Error(xf_DestroyEqnSet(All->EqnSet, xfe_True));
      if (ierr != xf_OK) return ierr;
      
      All->EqnSet = EqnSet;
      All->EqnSet->Dim = All->Mesh->Dim;
    }
    
    EqnSets[iSnap] = EqnSet;
    
  } // iSnap
  
  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  /*----------------*/
  /* Load snapshots */
  /*----------------*/

  /* Create a VectorSet for the snapshots: SnapSet */
  ierr = xf_Error(xf_CreateVectorSet(nSnap, &SnapSet));
  if (ierr != xf_OK) return ierr;


  /* Loop over snapshots and read them into SnapSet */
  xf_printf("Loading snapshots ... ");

  for (iSnap=0; iSnap<nSnap; iSnap++){

    sprintf(DataFile, "%s_%d.data", DataRoot, iSnap);

    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, DataFile, DataSet));
    if (ierr != xf_OK) return ierr;
    

    // pull off state vector -> store in SnapSet[iSnap]
    ierr = xf_Error(xf_FindPrimalState(DataSet, 0, &D, NULL));
    if (ierr != xf_OK) return ierr;
    U = (xf_Vector *) D->Data;
    SnapSet->Vector[iSnap] = *U;
    xf_InitVector(U); // so that snapshot is not lost when DataSet is destroyed

    // destroy DataSet
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;

  } // iSnap
  xf_printf("done.\n");


  /*-----------------------------*/
  /* Calculate POD Basis vectors */
  /*-----------------------------*/

  xf_printf("Calculating basis vectors from snapshots ...\n");

  /* Create a VectorSet for holding the Basis vectors in All-DataSet */
  ierr = xf_Error(xf_FindSimilarVectorSet(All, SnapSet->Vector+0, nBasis, "BasisSet",
					  xfe_True, xfe_True, &D, &BasisSet));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True;


  /* Perform POD -> basis vectors. */
  xf_printf("Calculating basis vectors from snapshots (POD) ...\n");
  ierr = xf_Error(xf_VectorSetPODMass(All, SnapSet, nSnap, nBasis, BasisSet));
  if (ierr != xf_OK) return ierr;
  xf_printf("done.\n");


  /* Dynamically load eqnset library */
  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;

  /* Register the equation set: assume aspects of all eqnsets are
     covered with this one registration. */
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;

  // number of outputs
  nOutput = All->EqnSet->Outputs->nOutput;

  if (myRank == 0){
    // create reduced model
    ierr = xf_Error(xf_CreateReducedModel(&RM));
    if (ierr != xf_OK) return ierr;
    
    // Set RM->EqnSet
    RM->EqnSet = All->EqnSet;
    
    // allocate linear portion of reduced model
    ierr = xf_Error(xf_AllocReducedModelLinear(RM, nBasis, nOutput));
    if (ierr != xf_OK) return ierr;
  }
  else 
    RM = NULL;

  /*-----------------------------------*/
  /* Linear terms of the reduced model */
  /*-----------------------------------*/

  xf_printf("\nConstructing linear terms of the reduced model ... \n");

  /* Turn off non-linear terms */
  ResTerms = All->EqnSet->ResTerms;
  for (i=0; i<ResTerms->nResTerm; i++)
    if (!xf_ResTermIsLinear(ResTerms->ResTerm + i))
      ResTerms->ResTerm[i].Active = xfe_False;


  // create local residual vectors
  ierr = xf_Error(xf_FindSimilarVector(All, BasisSet->Vector+0, "Residual_0",
				       xfe_True, xfe_False, NULL, &R0, NULL));
  if (ierr != xf_OK) return ierr;

    
  // vector set of residuals for speed performance in A_ij construction
  ierr = xf_Error(xf_FindSimilarVectorSet(All, BasisSet->Vector+0, nBasis, "Residuals",
					  xfe_True, xfe_False, NULL, &Residuals));
  if (ierr != xf_OK) return ierr;

  // calculate R0 = Residual evaluated with 0 state (linear terms only)
  ierr = xf_Error(xf_SetZeroVector(Residuals->Vector+0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_CalculateResidual(All, Residuals->Vector+0, R0, NULL, NULL));
  if (ierr != xf_OK) return ierr;

  // construct L_i = (phi_i, R0)
  for (i=0; i<nBasis; i++){
    ierr = xf_Error(xf_VectorDot(BasisSet->Vector+i, R0, &dp));
    if (ierr != xf_OK) return ierr;
    if (myRank == 0) RM->L[i] = dp;
  }

  // calculate residuals corresponding to all BasisSet vectors: Rj = R(phi_j) - R0
  for (j=0; j<nBasis; j++){    
    ierr = xf_Error(xf_CalculateResidual(All, BasisSet->Vector+j, 
					 Residuals->Vector+j, NULL, NULL));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_VectorMultSet(R0, -1.0, xfe_Add, Residuals->Vector+j));
    if (ierr != xf_OK) return ierr;
  } // j

  // construct A_ij = (phi_i, Rj)
  for (i=0; i<nBasis; i++){
    for (j=0; j<nBasis; j++){
      ierr = xf_Error(xf_VectorDot(BasisSet->Vector+i, Residuals->Vector+j, &dp));
      if (ierr != xf_OK) return ierr;
      if (myRank == 0) RM->A[i*RM->N+j] = dp;
    }
  }

  // destroy Residuals
  ierr = xf_Error(xf_DestroyVectorSet(Residuals));
  if (ierr != xf_OK) return ierr;

  xf_printf("done.\n");

  if ((DebugFlag) && (myRank == 0)){
    xf_printf("A = [\n");
    for (i=0; i<nBasis; i++){
      for (j=0; j<nBasis; j++){
	xf_printf("%.10E ", RM->A[i*nBasis+j]);
      }
    xf_printf("\n");
    }
    xf_printf("]\n");
    
    xf_printf("L = [\n");
    for (i=0; i<nBasis; i++)
      xf_printf("%.10E\n", RM->L[i]);
    xf_printf("]\n");
  }


  /*---------*/
  /* Outputs */
  /*---------*/

  xf_printf("\nConstructing reduced outputs ... \n");

  ierr = xf_Error(xf_SetZeroVector(R0));
  if (ierr != xf_OK) return ierr;

  for (i=0; i<nOutput; i++){
    
    // zero contribution to output (if any)
    ierr = xf_Error(xf_CalculateOutput(All, All->EqnSet->Outputs->Output[i].Name, 
				       R0, &val, NULL, xfe_Set));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0) RM->F0[i] = val;

    for (j=0; j<nBasis; j++){
      ierr = xf_Error(xf_CalculateOutput(All, All->EqnSet->Outputs->Output[i].Name, 
					 BasisSet->Vector+j, &val, NULL, xfe_Set));
      if (ierr != xf_OK) return ierr;
      
      if (myRank == 0) RM->F[i*RM->N+j] = val;
    } // j

  } // i
  

  xf_printf("done.\n");


  /*----*/
  /* UN */
  /*----*/

  if (UNFlag){
    xf_printf("\nTaking inner product of snapshots with basis vectors ... \n");

    if (myRank == 0){ // allocate space for UN
      RM->nSnap = nSnap;
      ierr = xf_Error(xf_Alloc((void **) &RM->UN, nSnap*nBasis, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }

    // construct UN_ij = (snapshot_i, phi_j)
    for (i=0; i<nSnap; i++){

      // R0 = SnapSet->Vector+i
      ierr = xf_Error(xf_SetVector(SnapSet->Vector+i, xfe_Set, R0)); 
      if (ierr != xf_OK) return ierr;

      // Multiply R0 by Mass matrix do be consistent with continuous norm
      // i.e. the basis functions are orthonormal w.r.t (a,b) = a^T*M*b
      ierr = xf_Error(xf_MultMassMatrix(All, 1.0, R0)); 
      if (ierr != xf_OK) return ierr;

      for (j=0; j<nBasis; j++){	
	ierr = xf_Error(xf_VectorDot(R0, BasisSet->Vector+j, &dp));
	if (ierr != xf_OK) return ierr;
	if (myRank == 0) RM->UN[i*nBasis+j] = dp;
      }
    } // i

    xf_printf("done.\n");
  }


  // destroy R0
  ierr = xf_Error(xf_DestroyVector(R0, xfe_True));
  if (ierr != xf_OK) return ierr;


  
  /*--------------------------------------*/
  /* Nonlinear terms of the reduced model */
  /*--------------------------------------*/

  xf_printf("\nReducing the nonlinear terms ...\n");

  EqnSet = All->EqnSet;
  ResTerms = EqnSet->ResTerms;

  // determine number of nonlinear terms
  nNonLinear = 0;
  for (i=0; i<ResTerms->nResTerm; i++)
    if (!xf_ResTermIsLinear(ResTerms->ResTerm + i))
      nNonLinear++;

  if (nNonLinear > 1) return xf_Error(xf_NOT_SUPPORTED);

  if (myRank == 0){
    // allocate memory for nonlinear portion of RM
    ierr = xf_Error(xf_AllocReducedModelNonLinear(RM, nNonLinear));
    if (ierr != xf_OK) return ierr;
  }
 
  // loop over nonlinear terms
  QuadOrder = NULL;
  for (iNonLinear=0; iNonLinear<nNonLinear; iNonLinear++){

    // pull off ResTerm
    for (i=0, j=0; i<ResTerms->nResTerm; i++)
      if ( (!xf_ResTermIsLinear(ResTerms->ResTerm + i)) && (j == iNonLinear)){
	ResTerm = ResTerms->ResTerm + i;
	iResTerm = i;
	break;
      }
    if (i >= ResTerms->nResTerm) return xf_Error(xf_CODE_LOGIC_ERROR);

    // store ResTerm in RM
    if (myRank == 0){
      ierr = xf_Error(xf_CopyResTerm(ResTerm, RM->ResTerm + iNonLinear));
      if (ierr != xf_OK) return ierr;
    }

    N = nBasis;  M = nSBasis;
    if (myRank == 0) RM->M[iNonLinear] = M;

    /* Calculate nonlinear snapshots */
    xf_printf("Computing the nonlinear snapshots. \n");
    ierr = xf_Error(xf_CalculateNonlinearSnap(All, SnapSet, EqnSets, iResTerm, 
					      OrderIncrement, &QuadOrder, &SSnapSet));
    if (ierr != xf_OK) return ierr;



    // Allocate memory for basis vectors
    ierr = xf_Error(xf_FindSimilarVectorSet(All, SSnapSet->Vector+0, M, "SBasisSet",
					    xfe_True, xfe_False, &D, &SBasisSet));
    if (ierr != xf_OK) return ierr;


    /* Scale SSnapSet by quadrature weights */
    /*  ierr = xf_Error(xf_ScaleVectorSetByWeights(All, SSnapSet, nSnap, */
    /* 					       QuadOrder, xfe_False)); */
    /*     if (ierr != xf_OK) return ierr; */

    /* Perform POD on SSnapSet */
    xf_printf("Computing the nonlinear POD vectors. \n");
    ierr = xf_Error(xf_VectorSetPODMass(NULL, SSnapSet, nSnap, M, SBasisSet));
    if (ierr != xf_OK) return ierr;
    xf_printf("done.\n");

    /* Scale SSnapSet by 1/(quadrature weights) */
    /*     ierr = xf_Error(xf_ScaleVectorSetByWeights(All, SSnapSet, M, */
    /* 					       QuadOrder, xfe_True)); */
    /*     if (ierr != xf_OK) return ierr; */


    /* Trim SSnapSet to remove unused (nSnap - M) vectors */
    ierr = xf_Error(xf_TrimVectorSet(&SSnapSet, M));
    if (ierr != xf_OK) return ierr;

    // rename SSnapSet -> SCardSet (to save release/alloc)
    SCardSet = SSnapSet;
    SSnapSet = NULL;

    // calculate EIM points (RM->z) and expansion functions SCardSet
    xf_printf("Performing the coefficient expansion. \n");
    z = ((myRank == 0) ? RM->z+iNonLinear : NULL);
    ierr = xf_Error(xf_CoeffExpansionEIM(All, SBasisSet, M, QuadOrder, z, SCardSet));
    if (ierr != xf_OK) return ierr;
    
    // debug output
    if ((DebugFlag) && (nProc == 1)){
      xf_printf("z = [\n");
      for (i=0; i<M; i++){
	xf_printf("%d %d %d %d ", i, RM->z[iNonLinear].egrp[i], RM->z[iNonLinear].elem[i],
		  RM->z[iNonLinear].node[i]);
	ierr = xf_Error(xf_Ref2GlobElem(All->Mesh, RM->z[iNonLinear].egrp[i],
					RM->z[iNonLinear].elem[i], NULL, xfe_True,
					1, RM->z[iNonLinear].xref+i*RM->z[iNonLinear].Dim,
					xglob));
	if (ierr != xf_OK) return ierr;
	for (j=0; j<RM->z[iNonLinear].Dim; j++)
	  xf_printf("%.10E ", xglob[j]);
	xf_printf("\n");
      }
      xf_printf("]\n");
    }

    // allocate memory for E
    ierr = xf_Error(xf_Alloc( (void **) &Emat, N*M, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    // calculate E
    xf_printf("Calculating E. \n");
    ierr = xf_Error(xf_NonlinearReducedMatrixE(All, BasisSet, SCardSet, N, M, 
					       QuadOrder, Emat));
    if (ierr != xf_OK) return ierr;

    // only retain E on root
    if (myRank == 0)
      RM->E[iNonLinear] = Emat;
    else
      xf_Release( (void *) Emat);

      
    if ((DebugFlag) && (myRank == 0)){ // DEBUG
      xf_printf("E = [\n");
      for (i=0; i<N; i++){
	for (j=0; j<M; j++){
	  xf_printf("%.10E ", RM->E[iNonLinear][i*M+j]);
	}
	xf_printf("\n");
      }
      xf_printf("]\n");
    }

    // allocate memory for D
    ierr = xf_Error(xf_Alloc( (void **) &Dmat, M*N, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    // calculate D
    xf_printf("Calculating D. \n");
    z = ((myRank == 0) ? RM->z+iNonLinear : NULL);
    ierr = xf_Error(xf_NonlinearReducedMatrixD(BasisSet, z, N, M, Dmat));
    if (ierr != xf_OK) return ierr;

    // only retain D on root
    if (myRank == 0)
      RM->D[iNonLinear] = Dmat;
    else
      xf_Release( (void *) Dmat);

    if ((DebugFlag) && (myRank == 0)){ // DEBUG
      xf_printf("D = [\n");
      for (i=0; i<M; i++){
	for (j=0; j<N; j++){
	  xf_printf("%.10E ", RM->D[iNonLinear][i*N+j]);
	}
	xf_printf("\n");
      }
      xf_printf("]\n");
    }

    // destroy card set
    ierr = xf_Error(xf_DestroyVectorSet(SCardSet));
    if (ierr != xf_OK) return ierr;

    /* Destroy nonlinear basis vectors */
    ierr = xf_Error(xf_DestroyVectorSet(SBasisSet));
    if (ierr != xf_OK) return ierr;

  } // iNonLinear


  xf_Release( (void *) QuadOrder);

  xf_printf("Finished reducing the nonlinear terms.\n");

  // write reduced RM
  if (myRank == 0){ // root performs the writing
    
    // binary
    sprintf(OutputFile, "%s_N%d_M%d.rom\0", DataRoot, nBasis, nSBasis);
    if ((fid = fopen(OutputFile, "wb")) == NULL)  return xf_Error(xf_FILE_WRITE_ERROR);
    ierr = xf_Error(xf_WriteReducedModelBinary(RM, fid));
    if (ierr!=xf_OK) return ierr;
    if (fclose(fid)!= 0) return xf_Error(xf_FILE_WRITE_ERROR);
    // ascii
    sprintf(OutputFile, "%s_N%d_M%d.m\0", DataRoot, nBasis, nSBasis);
    if ((fid = fopen(OutputFile, "w")) == NULL)  return xf_Error(xf_FILE_WRITE_ERROR);
    ierr = xf_Error(xf_WriteReducedModelAscii(RM, fid));
    if (ierr!=xf_OK) return ierr;
    if (fclose(fid)!= 0) return xf_Error(xf_FILE_WRITE_ERROR);
  }

  /* Destroy snapshots (only from memory) */
  ierr = xf_Error(xf_DestroyVectorSet(SnapSet));
  if (ierr != xf_OK) return ierr;

  // Write a basis .data file
  sprintf(OutputFile, "%s_Basis_N%d.data\0", DataRoot, nBasis);
  ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, All->DataSet, NULL, OutputFile));
  if (ierr != xf_OK) return ierr;

  /* Destroy Reduced Model */
  ierr = xf_Error(xf_DestroyReducedModel(RM, xfe_False));
  if (ierr != xf_OK) return ierr;


  /* Destroy EqnSets */
  for (iSnap=0; iSnap<nSnap; iSnap++){
    ierr = xf_Error(xf_DestroyEqnSet(EqnSets[iSnap], xfe_True));
    if (ierr != xf_OK) return ierr;
  }
  All->EqnSet = NULL; // don't want to destroy the first eqnset twice

  xf_Release( (void *) EqnSets);

  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;
  
 
  xf_printf("xf_MROffline finished.\n");


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, runerr;

  /* Initialize parallel-run (no effect in serial) */
  ierr = xf_Error(xf_MPI_Init(&argc, &argv));
  if (ierr != xf_OK) return ierr;
  
  runerr = xf_Error(mymain(argc, argv));

  /* MPI finalize (no effect in serial) */
  ierr = xf_Error(xf_MPI_Finalize());
  if (ierr != xf_OK) return ierr;

  return runerr;
}

