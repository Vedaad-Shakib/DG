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
 FILE:  xf_LinearSolverILU.c
 
 This file contains linear solver functions relating to the incomplete
 LU (ILU) factorization.
 
*/

#include "xf_MeshTools.h"
#include "xf_Line.h"
#include "xf_Heap.h"

/******************************************************************/
//   FUNCTION Definition: xf_CreateILUData
static int
xf_CreateILUData(xf_All *All, xf_JacobianMatrix *R_U, xf_LinearSolverILUData **pILUData)
{  
  int ierr, sr;
  int negrp, egrp, elem, nelemtot;
  int *nElem = NULL;
  xf_LinearSolverILUData *ILUData;
  xf_Vector *CVector;
  xf_Mesh *Mesh;

  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup; // number of element groups
  sr    = R_U->StateRank;

  // allocate memory for the data structure
  ierr = xf_Error(xf_Alloc( (void **) pILUData, 1, sizeof(xf_LinearSolverILUData)));
  if (ierr != xf_OK) return ierr;
  ILUData = *pILUData;

  // determine total number of elements
  ierr = xf_Error(xf_GetnElem(Mesh, &nElem, &ILUData->nelemtot));
  if (ierr != xf_OK) return ierr;

  nelemtot = ILUData->nelemtot;

  // allocate Index2Elem [nelemtot x 2]
  ierr = xf_Error(xf_Alloc2( (void ***) &ILUData->Index2Elem, nelemtot, 2, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // allocate Elem2Index [negrp, nelem[egrp]]
  ierr = xf_Error(xf_VAlloc2( (void ***) &ILUData->Elem2Index, negrp, nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // allocate C matrix
  ierr = xf_Error(xf_FindLineConnectivity(All, &CVector));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &ILUData->C, negrp, sizeof(real **)));
  if (ierr != xf_OK) return ierr;
  for (egrp=0; egrp<negrp; egrp++)
    ILUData->C[egrp] = CVector->GenArray[egrp].rValue;

  // determine rmax = max(r), and nfacemax = max(nface) over all elements
  ILUData->rmax=0;
  ILUData->nfacemax=0;
  for (egrp=0; egrp<negrp; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      ILUData->rmax     = max(  sr*xf_Jacobian_n(R_U,egrp,elem), ILUData->rmax    );
      ILUData->nfacemax = max(Mesh->ElemGroup[egrp].nFace[elem], ILUData->nfacemax);
    }
  
  // release memory
  xf_Release( (void *) nElem);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyILUData
int
xf_DestroyILUData(xf_LinearSolverILUData *ILUData)
{  
  if (ILUData == NULL) return xf_OK;

  xf_Release2( (void **) ILUData->Index2Elem);
  xf_Release2( (void **) ILUData->Elem2Index);
  xf_Release ( (void  *) ILUData->C         );
  
  xf_Release( (void *) ILUData); // destroy self

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_JacobianNeighbor
static void
xf_JacobianNeighbor(xf_JacobianMatrix *R_U, int egrp, int elem, int face,
		    int *pegrpN, int *pelemN, int *pfaceN)
{  
  // returns info on neighbor using information in R_U
  if (pegrpN != NULL) (*pegrpN) = R_U->egrpN[egrp][elem][face];
  if (pelemN != NULL) (*pelemN) = R_U->elemN[egrp][elem][face];
  if (pfaceN != NULL) (*pfaceN) = R_U->faceN[egrp][elem][face];
}

/******************************************************************/
//   FUNCTION Definition: xf_ILUOrderElements_MDF
static int
xf_ILUOrderElements_MDF(xf_All *All, xf_JacobianMatrix *R_U)
{  
  /* Reference = Persson & Peraire 2008 */
  int ierr, sr, k, index, kie;
  int negrp, nelemtot;
  int ni, nj, rmax;
  int iegrp, ielem, iface, ifacei, ifacej, ifacek;
  int jegrp, jelem, jfacei, jfacek;
  int kegrp, kelem, kfacei, kfacej;
  int legrp, lelem, lfacek;
  int *P = NULL;
  real w = 0.0, dw;
  real * RSTRCT Aii = NULL;
  real * RSTRCT Aij = NULL;
  real ***C = NULL;
  real *R_Uij;
  xf_LinearSolverILUData *ILUData;
  xf_HeapNode *HeapNode, *H;
  xf_Heap *Heap;
  xf_Mesh *Mesh;

  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup;
  sr    = R_U->StateRank;

  // ILUData needs to be allocated prior to call
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);

  C        = ILUData->C; // connectivity matrix
  nelemtot = ILUData->nelemtot;  // total number of elements

  // allocate matrices/vectors
  rmax = ILUData->rmax;
  ierr = xf_Error(xf_Alloc( (void **) &Aii, 2*rmax*rmax, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  Aij = Aii+rmax*rmax;
  ierr = xf_Error(xf_Alloc( (void **) &P, rmax, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  

  /* Fill in C matrix */
  for (iegrp=0; iegrp<negrp; iegrp++){
    for (ielem=0; ielem<Mesh->ElemGroup[iegrp].nElem; ielem++){

      ni = xf_Jacobian_n(R_U,iegrp,ielem); // r/StateRank for iegrp,ielem

      // set Aii = PLU factorization of R_U{ii}
      xf_V_Add(R_U->Value[iegrp][ielem][0], ni*ni*sr*sr, xfe_Set, Aii);
      ierr = xf_Error(xf_ComputeBlockPLU(Aii, ni, sr, P));
      if (ierr != xf_OK) return ierr;

      // loop over neighbors
      for (ifacej=0; ifacej<Mesh->ElemGroup[iegrp].nFace[ielem]; ifacej++){

	// neighboring element
	xf_JacobianNeighbor(R_U, iegrp, ielem, ifacej, &jegrp, &jelem, NULL);

	// skip halos and boundaries
	if ((jegrp<0) || (jegrp>=negrp)) continue;
     
	nj = xf_Jacobian_n(R_U,jegrp,jelem); // r/StateRank

	// Aij = Aii^{-1}*R_U{ij}
	R_Uij = R_U->Value[iegrp][ielem][ifacej + 1];
        ierr = xf_Error(xf_SolveBlockPLU_Matrix(Aii, ni, sr, P, R_Uij, nj, xfe_Set, Aij));
        if (ierr != xf_OK) return ierr;

	// Cij = norm(Aij)
	C[iegrp][ielem][ifacej] = xf_RealNorm(Aij, ni*nj*sr*sr);

      } // ifacej
      
    } // ielem
  } // iegrp


  // allocate Heap node data structure (will get destroyed with heap)
  ierr = xf_Error(xf_Alloc( (void **) &HeapNode, nelemtot, sizeof(xf_HeapNode)));
  if (ierr != xf_OK) return ierr;

  /* Calculate initial weights */
  index = 0;
  for (kegrp=0; kegrp<negrp; kegrp++){
    for (kelem=0; kelem<Mesh->ElemGroup[kegrp].nElem; kelem++){
      w = 0.;
      // loop over neighbor pairs of k: i and j
      for (kfacei=0; kfacei<Mesh->ElemGroup[kegrp].nFace[kelem]; kfacei++){

	// (iegrp,ielem) = neighbor of k across kfacei
	xf_JacobianNeighbor(R_U, kegrp, kelem, kfacei, &iegrp, &ielem, &ifacek);

	// skip halos and boundaries
	if ((iegrp<0) || (iegrp>=negrp)) continue;

	for (kfacej=0; kfacej<Mesh->ElemGroup[kegrp].nFace[kelem]; kfacej++){

	  if (kfacej == kfacei) continue;

	  // (jegrp,jelem) = neighbor of k across kfacej
	  xf_JacobianNeighbor(R_U, kegrp, kelem, kfacej, &jegrp, &jelem, &jfacek);
	  
	  // skip halos and boundaries
	  if ((jegrp<0) || (jegrp>=negrp)) continue;

	  // are i and j connected?
	  ierr = xf_CommonFace(Mesh, iegrp, ielem, jegrp, jelem, &ifacej, &jfacei);
	  if (ierr == xf_NOT_FOUND){
	    // add to weight only if i and j are not connected
	    dw = C[jegrp][jelem][jfacek]*C[kegrp][kelem][kfacei];
	    w += dw*dw; // using Frobenius norm
	  }
	    
	} // kfacej
      } // kfacei
      H = HeapNode+(index++);                    // pointer to heap entry
      H->rval = sqrt(w);                         // store weight in heap
      // No need to set H->ipos since we are visiting (kegrp,kelem) in order
    } // kelem
  } // kegrp

  // sanity check
  if (index != nelemtot) return xf_Error(xf_OUT_OF_BOUNDS);


  /* Build a min-heap out of nodes in HeapNode */
  ierr = xf_Error(xf_BuildMinHeap(HeapNode, nelemtot, &Heap));
  if (ierr != xf_OK) return ierr;
  
  // initialize Index2Elem and Elem2Index
  for (k=0; k<nelemtot*2; k++) ILUData->Index2Elem[0][k] = -1;
  for (k=0; k<nelemtot  ; k++) ILUData->Elem2Index[0][k] = -1;

  
  /* Fill in Index2Elem and Elem2Index (i.e. order the elements) */
  for (index=0; index<nelemtot; index++){

    // pull off minimum weight (fill-in) element
    xf_MinHeapTakeAwayRoot(Heap, &H);
    ierr = xf_Error(xf_Index2EgrpElem(Mesh, H->ipos, &legrp, &lelem));
    if (ierr != xf_OK) return ierr;

    // store data in Index2Elem and Elem2Index
    ILUData->Index2Elem[index][0] = legrp;
    ILUData->Index2Elem[index][1] = lelem;
    ILUData->Elem2Index[legrp][lelem] = index;
    //xf_printf("index=%d: legrp=%d, lelem=%d, weight=%.8E\n", index, legrp, lelem, H->rval);
    
    // loop over neighbors of (legrp, lelem) -> k
    for (lfacek=0; lfacek<Mesh->ElemGroup[legrp].nFace[lelem]; lfacek++){
      
      // (kegrp,kelem) = neighbor of l across lfacek
      xf_JacobianNeighbor(R_U, legrp, lelem, lfacek, &kegrp, &kelem, NULL);

      // skip halos and boundaries
      if ((kegrp<0) || (kegrp>=negrp)) continue;

      if (ILUData->Elem2Index[kegrp][kelem] != -1) continue; // already ordered

      // recompute weight for elem k using only unordered elements
      dw = 0.;
      // loop over neighbor pairs of k: i and j
      for (kfacei=0; kfacei<Mesh->ElemGroup[kegrp].nFace[kelem]; kfacei++){

	// (iegrp,ielem) = neighbor of k across kfacei
	xf_JacobianNeighbor(R_U, kegrp, kelem, kfacei, &iegrp, &ielem, &ifacek);

	// skip halos and boundaries
	if ((iegrp<0) || (iegrp>=negrp)) continue;

	if (ILUData->Elem2Index[iegrp][ielem] != -1) continue; // already ordered

	for (kfacej=0; kfacej<Mesh->ElemGroup[kegrp].nFace[kelem]; kfacej++){

	  if (kfacej == kfacei) continue;  // i and j must be different

	  // (jegrp,jelem) = neighbor of k across kfacej
	  xf_JacobianNeighbor(R_U, kegrp, kelem, kfacej, &jegrp, &jelem, &jfacek);

	  // skip halos and boundaries
	  if ((jegrp<0) || (jegrp>=negrp)) continue;

	  if (ILUData->Elem2Index[jegrp][jelem] != -1) continue; // already ordered

	  // are i and j connected?
	  ierr = xf_CommonFace(Mesh, iegrp, ielem, jegrp, jelem, &ifacej, &jfacei);
	  if (ierr == xf_NOT_FOUND){
	    // add to weight only if i and j are not connected
	    dw = C[jegrp][jelem][jfacek]*C[kegrp][kelem][kfacei];
	    w += dw*dw; // using Frobenius norm
	  }

	} // kfacej
      } // kfacei

      w = sqrt(w);

      // pull off global index of (kegrp, kelem)
      ierr = xf_Error(xf_EgrpElem2Index(Mesh, kegrp, kelem, &kie));
      if (ierr != xf_OK) return ierr;

      // modify weight
      xf_MinHeapChangeNode(Heap, kie, w);
      
    } // lfacek    

  } // index  

  
  // destroy heap
  xf_DestroyHeap(Heap);

  // Release memory
  xf_Release( (void *) Aii);
  xf_Release( (void *) P);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ILUOrderElements_Lex
static int
xf_ILUOrderElements_Lex(xf_All *All, xf_JacobianMatrix *R_U)
{  
  /* Lexicographical ordering (as provided in the mesh) */
  int ierr, index;
  int k, legrp, lelem;
  xf_LinearSolverILUData *ILUData;
  xf_Mesh *Mesh;

  Mesh  = All->Mesh;

  // ILUData needs to be allocated prior to call
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // initialize Index2Elem and Elem2Index
  for (k=0; k<ILUData->nelemtot*2; k++) ILUData->Index2Elem[0][k] = -1;
  for (k=0; k<ILUData->nelemtot  ; k++) ILUData->Elem2Index[0][k] = -1;

  
  /* Fill in Index2Elem and Elem2Index (i.e. order the elements) */
  for (index=0; index<ILUData->nelemtot; index++){

    ierr = xf_Error(xf_Index2EgrpElem(Mesh, index, &legrp, &lelem));
    if (ierr != xf_OK) return ierr;

    // store data in Index2Elem and Elem2Index
    ILUData->Index2Elem[index][0] = legrp;
    ILUData->Index2Elem[index][1] = lelem;
    ILUData->Elem2Index[legrp][lelem] = index;
  } // iindex    

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ILUOrderElements
static int
xf_ILUOrderElements(xf_All *All, xf_JacobianMatrix *R_U)
{  
  int ierr;
  enum xfe_ILUOrderingType ILUOrdering;

  // determine type of ordering requested
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "ILUOrdering", 
                                     xfe_ILUOrderingName, (int ) xfe_ILUOrderingLast, 
                                     (int *) &ILUOrdering));
  if (ierr != xf_OK) return ierr;

  // call appropriate ordering function
  switch (ILUOrdering){
  case xfe_ILUOrderingMDF:
    ierr = xf_Error(xf_ILUOrderElements_MDF(All, R_U));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_ILUOrderingLex:
    ierr = xf_Error(xf_ILUOrderElements_Lex(All, R_U));
    if (ierr != xf_OK) return ierr;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Factorize_ILU0
static int
xf_Factorize_ILU0(xf_All *All, xf_JacobianMatrix *R_U)
{  
  int ierr, negrp;
  int sr, sr2, rmax, k;
  int iindex, jindex, kindex;
  int ni, nj, nk;
  int iegrp, ielem, ifacej, ifacek;
  int jegrp, jelem, jfacei, jfacek;
  int kegrp, kelem, kfacei, kfacej;
  int *P;
  real *Rii, *Rij, *Rji, *Rik, *Rjk;
  real *T    = NULL;
  real *iRii = NULL;
  xf_LinearSolverILUData *ILUData;
  xf_Mesh *Mesh;

  // ILUData needs to exist
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // store desired quantities
  Mesh     = All->Mesh;
  negrp    = Mesh->nElemGroup;
  sr       = R_U->StateRank;
  sr2      = sr*sr;
  rmax     = ILUData->rmax;

  // allocate a temporary matrix
  ierr = xf_Error(xf_ReAlloc( (void **) &T   , rmax*rmax, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  // alocate another matrix for storing inverses of block diagonals
  ierr = xf_Error(xf_ReAlloc( (void **) &iRii, rmax*rmax, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  

  // loop over all elements, i
  for (iindex=0; iindex<ILUData->nelemtot; iindex++){

    // obtain ordered element (iegrp, ielem)
    iegrp = ILUData->Index2Elem[iindex][0];
    ielem = ILUData->Index2Elem[iindex][1];
    
    Rii = R_U->Value[iegrp][ielem][0];    // on-diagonal block
    ni  = xf_Jacobian_n(R_U,iegrp,ielem); // r/StateRank for iegrp,ielem
    P   = R_U->P[iegrp][ielem];           // permutation vector for ielem

    // PLU factor Rii
    ierr = xf_Error(xf_ComputeBlockPLU(Rii, ni, sr, P));
    if (ierr != xf_OK) return ierr;

    // last element; nothing else to do
    if (iindex ==ILUData->nelemtot-1) continue;

    // Set iRii = (Rii)^{-1}
    xf_BlockIdentity(ni, sr, T);
    ierr = xf_Error(xf_SolveBlockPLU_Matrix(Rii, ni, sr, P, T, ni, xfe_Set, iRii));
    if (ierr != xf_OK) return ierr;

    // loop over neighbors, j, but consider only higher-ordered elements
    for (ifacej=0; ifacej<Mesh->ElemGroup[iegrp].nFace[ielem]; ifacej++){
      // obtain jegrp, jelem, jfacei
      xf_JacobianNeighbor(R_U, iegrp, ielem, ifacej, &jegrp, &jelem, &jfacei);

      // skip halos and boundaries
      if ((jegrp<0) || (jegrp>=negrp)) continue;

      // ordered index of j
      jindex = ILUData->Elem2Index[jegrp][jelem];
      if (jindex <= iindex) continue;  // this neighbor was already visited

      nj  = xf_Jacobian_n(R_U,jegrp,jelem); // r/StateRank for jegrp,jelem
      
      // pointers into Jacobian off-diagonal blocks
      Rij = R_U->Value[iegrp][ielem][ifacej + 1];
      Rji = R_U->Value[jegrp][jelem][jfacei + 1];

      // Set L of ILU factorization: Lji = Rji*Rii^{-1}  (stored back in Rji)
      ierr = xf_Error(xf_BlockMxBlockM(Rji, nj, sr, ni, iRii, ni, xfe_Set, T));
      if (ierr != xf_OK) return ierr;
      xf_V_Add(T, ni*nj*sr2, xfe_Set, Rji);

      // loop over elements k that are neighbors of both j and i (also include j=k)
      for (jfacek=-1; jfacek<Mesh->ElemGroup[jegrp].nFace[jelem]; jfacek++){

	if (jfacek == -1){ // self block (j=k)
	  kegrp = jegrp;
	  kelem = jelem;
	  ifacek = ifacej;
	}
	else{
	  // obtain kegrp, kelem, kfacej
	  xf_JacobianNeighbor(R_U, jegrp, jelem, jfacek, &kegrp, &kelem, &kfacej);
	  // skip halos and boundaries
	  if ((kegrp<0) || (kegrp>=negrp)) continue;
	  // ordered index of k
	  kindex = ILUData->Elem2Index[kegrp][kelem];
	  if (kindex <= iindex) continue;  // this neighbor was already visited
	  // are i and k connected?
	  ierr = xf_CommonFace(Mesh, iegrp, ielem, kegrp, kelem, &ifacek, &kfacei);
	  if (ierr == xf_NOT_FOUND) continue; // if not, nothing to do
	}
	  
	nk  = xf_Jacobian_n(R_U,kegrp,kelem); // r/StateRank for kegrp,kelem

	// pointers into Jacobian
	Rik = R_U->Value[iegrp][ielem][ifacek + 1];
	Rjk = R_U->Value[jegrp][jelem][jfacek + 1];
	
	// Rjk -= Lji*Rik
	ierr = xf_Error(xf_BlockMxBlockM(Rji, nj, sr, ni, Rik, nk, xfe_Sub, Rjk));
	if (ierr != xf_OK) return ierr;

      } // jfacek
    } // ifacej
  } // iindex


  // release memory
  xf_Release((void *)    T);
  xf_Release((void *) iRii);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_Precondition_ILU0
static int
xf_Jacobian_Precondition_ILU0(xf_All *All, xf_JacobianMatrix *R_U)
{  
  int ierr;

  // delete any existing R_U->ILUData
  ierr = xf_Error(xf_DestroyILUData(R_U->ILUData));
  if (ierr != xf_OK) return ierr;

  // create new R_U->ILUData
  ierr = xf_Error(xf_CreateILUData(All, R_U, &R_U->ILUData));
  if (ierr != xf_OK) return ierr;
    
  // order the elements (e.g. to minimize discarded fill)
  ierr = xf_Error(xf_ILUOrderElements(All, R_U));
  if (ierr != xf_OK) return ierr;

  // perform incomplete factorization
  ierr = xf_Error(xf_Factorize_ILU0(All, R_U));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_ILU0_Lower
static int
xf_Jacobian_SolveM_ILU0_Lower(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X)
{  
  // X = L^{-1}*X
  int ierr, sr, negrp;
  int ni, nj;
  int iindex, jindex;
  int iegrp, ielem, ifacej;
  int jegrp, jelem, jfacei;
  real *Rij, *Xi, *Xj;
  xf_LinearSolverILUData *ILUData;
  xf_Mesh *Mesh;

  // ILUData needs to exist
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);

  // store desired quantities
  Mesh     = All->Mesh;
  negrp    = Mesh->nElemGroup;
  sr       = R_U->StateRank;

  // loop forward over elements, i
  for (iindex=0; iindex<ILUData->nelemtot; iindex++){
    
    // obtain ordered element (iegrp, ielem)
    iegrp = ILUData->Index2Elem[iindex][0];
    ielem = ILUData->Index2Elem[iindex][1];

    ni  = xf_Jacobian_n(R_U,iegrp,ielem);    // r/StateRank for iegrp, ielem
    Xi  = X->GenArray[iegrp].rValue[ielem];  // X for ielem

    // loop over neighbors, j, but consider only lower-ordered elements
    for (ifacej=0; ifacej<Mesh->ElemGroup[iegrp].nFace[ielem]; ifacej++){

      // obtain jegrp, jelem, jfacei
      xf_JacobianNeighbor(R_U, iegrp, ielem, ifacej, &jegrp, &jelem, &jfacei);
      
      // skip halos and boundaries
      if ((jegrp<0) || (jegrp>=negrp)) continue;

      // get ordered index of j and consider only jindex lower than iindex
      jindex = ILUData->Elem2Index[jegrp][jelem];
      if (jindex >= iindex) continue;  

      // set pointers
      Rij = R_U->Value[iegrp][ielem][ifacej + 1]; // off-diag Jacobian
      nj  = xf_Jacobian_n(R_U,jegrp,jelem);       // r/StateRank for jegrp,jelem
      Xj  = X->GenArray[jegrp].rValue[jelem];     // X for jelem
    
      // Xi -= Rij*Xj
      ierr = xf_Error(xf_BlockMxV(Rij, ni, sr, nj, Xj, xfe_Sub, Xi));
      if (ierr != xf_OK) return ierr;

    } // ifacej

  } // iindex

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_ILU0_LowerT
static int
xf_Jacobian_SolveM_ILU0_LowerT(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X)
{  
  // X = L^{-T}*X
  int ierr, sr, negrp;
  int ni, nj;
  int iindex, jindex;
  int iegrp, ielem, ifacej;
  int jegrp, jelem, jfacei;
  real *Rji, *Xi, *Xj;
  xf_LinearSolverILUData *ILUData;
  xf_Mesh *Mesh;

  // ILUData needs to exist
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);

  // store desired quantities
  Mesh     = All->Mesh;
  negrp    = Mesh->nElemGroup;
  sr       = R_U->StateRank;

  // loop backward over elements, i
  for (iindex=ILUData->nelemtot-1; iindex>=0; iindex--){
    
    // obtain ordered element (iegrp, ielem)
    iegrp = ILUData->Index2Elem[iindex][0];
    ielem = ILUData->Index2Elem[iindex][1];

    ni  = xf_Jacobian_n(R_U,iegrp,ielem);    // r/StateRank for iegrp, ielem
    Xi  = X->GenArray[iegrp].rValue[ielem];  // X for ielem

    // loop over neighbors, j, but consider only higher-ordered elements
    for (ifacej=0; ifacej<Mesh->ElemGroup[iegrp].nFace[ielem]; ifacej++){

      // obtain jegrp, jelem, jfacei
      xf_JacobianNeighbor(R_U, iegrp, ielem, ifacej, &jegrp, &jelem, &jfacei);

      // skip halos and boundaries
      if ((jegrp<0) || (jegrp>=negrp)) continue;

      // get ordered index of j and consider only jindex higher than iindex
      jindex = ILUData->Elem2Index[jegrp][jelem];
      if (jindex <= iindex) continue;  

      // set pointers
      Rji = R_U->Value[jegrp][jelem][jfacei + 1]; // off-diag Jacobian
      nj  = xf_Jacobian_n(R_U,jegrp,jelem);       // r/StateRank for jegrp,jelem
      Xj  = X->GenArray[jegrp].rValue[jelem];     // X for jelem
    
      // Xi -= (Rji)^T*Xj
      ierr = xf_Error(xf_BlockMTxV(Rji, ni, sr, nj, Xj, xfe_Sub, Xi));
      if (ierr != xf_OK) return ierr;

    } // ifacej

  } // iindex

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_ILU0_Upper
static int
xf_Jacobian_SolveM_ILU0_Upper(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X)
{  
  // X = U^{-1}*X
  int ierr, sr, negrp;
  int ni, nj;
  int iindex, jindex;
  int iegrp, ielem, ifacej;
  int jegrp, jelem, jfacei;
  int *P;
  real *Rij, *Rii, *Xi, *Xj;
  xf_LinearSolverILUData *ILUData;
  xf_Mesh *Mesh;

  // ILUData needs to exist
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);

  // store desired quantities
  Mesh     = All->Mesh;
  negrp    = Mesh->nElemGroup;
  sr       = R_U->StateRank;

  // loop backwards over elements, i
  for (iindex=ILUData->nelemtot-1; iindex>=0; iindex--){
    
    // obtain ordered element (iegrp, ielem)
    iegrp = ILUData->Index2Elem[iindex][0];
    ielem = ILUData->Index2Elem[iindex][1];

    ni  = xf_Jacobian_n(R_U,iegrp,ielem);    // r/StateRank for iegrp, ielem
    Xi  = X->GenArray[iegrp].rValue[ielem];  // X for ielem

    // loop over neighbors, j, but consider only higher-ordered elements
    for (ifacej=0; ifacej<Mesh->ElemGroup[iegrp].nFace[ielem]; ifacej++){

      // obtain jegrp, jelem, jfacei
      xf_JacobianNeighbor(R_U, iegrp, ielem, ifacej, &jegrp, &jelem, &jfacei);
      
      // skip halos and boundaries
      if ((jegrp<0) || (jegrp>=negrp)) continue;

      // get ordered index of j and consider only jindex higher than iindex
      jindex = ILUData->Elem2Index[jegrp][jelem];
      if (jindex <= iindex) continue;  

      // set pointers
      Rij = R_U->Value[iegrp][ielem][ifacej + 1]; // off-diag Jacobian
      nj  = xf_Jacobian_n(R_U,jegrp,jelem);       // r/StateRank for jegrp,jelem
      Xj  = X->GenArray[jegrp].rValue[jelem];     // X for jelem
    
      // Xi -= Rij*Xj
      ierr = xf_Error(xf_BlockMxV(Rij, ni, sr, nj, Xj, xfe_Sub, Xi));
      if (ierr != xf_OK) return ierr;

    } // ifacej
    
    Rii = R_U->Value[iegrp][ielem][0]; // on-diag Jacobian
    P   = R_U->P[iegrp][ielem];        // permutation vector for i

    // Xi = Rii^{-1} * Xi
    ierr = xf_Error(xf_SolveBlockPLU(Rii, ni, sr, P, Xi, xfe_Set, Xi));
    if (ierr != xf_OK) return ierr;


  } // iindex

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_ILU0_UpperT
static int
xf_Jacobian_SolveM_ILU0_UpperT(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X)
{  
  // X = U^{-T}*X
  int ierr, sr, negrp;
  int ni, nj;
  int iindex, jindex;
  int iegrp, ielem, ifacej;
  int jegrp, jelem, jfacei;
  int *P;
  real *Rji, *Rii, *Xi, *Xj;
  xf_LinearSolverILUData *ILUData;
  xf_Mesh *Mesh;

  // ILUData needs to exist
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);

  // store desired quantities
  Mesh     = All->Mesh;
  negrp    = Mesh->nElemGroup;
  sr       = R_U->StateRank;

  // loop forward over elements, i
  for (iindex=0; iindex<ILUData->nelemtot; iindex++){
    
    // obtain ordered element (iegrp, ielem)
    iegrp = ILUData->Index2Elem[iindex][0];
    ielem = ILUData->Index2Elem[iindex][1];

    ni  = xf_Jacobian_n(R_U,iegrp,ielem);    // r/StateRank for iegrp, ielem
    Xi  = X->GenArray[iegrp].rValue[ielem];  // X for ielem

    // loop over neighbors, j, but consider only lower-ordered elements
    for (ifacej=0; ifacej<Mesh->ElemGroup[iegrp].nFace[ielem]; ifacej++){

      // obtain jegrp, jelem, jfacei
      xf_JacobianNeighbor(R_U, iegrp, ielem, ifacej, &jegrp, &jelem, &jfacei);

      // skip halos and boundaries
      if ((jegrp<0) || (jegrp>=negrp)) continue;

      // get ordered index of j and consider only jindex lower than iindex
      jindex = ILUData->Elem2Index[jegrp][jelem];
      if (jindex >= iindex) continue;  

      // set pointers
      Rji = R_U->Value[jegrp][jelem][jfacei + 1]; // off-diag Jacobian
      nj  = xf_Jacobian_n(R_U,jegrp,jelem);       // r/StateRank for jegrp,jelem
      Xj  = X->GenArray[jegrp].rValue[jelem];     // X for jelem
    
      // Xi -= (Rji^T)*Xj
      ierr = xf_Error(xf_BlockMTxV(Rji, ni, sr, nj, Xj, xfe_Sub, Xi));
      if (ierr != xf_OK) return ierr;

    } // ifacej
    
    Rii = R_U->Value[iegrp][ielem][0]; // on-diag Jacobian
    P   = R_U->P[iegrp][ielem];        // permutation vector for i

    // Xi = Rii^{-1} * Xi
    ierr = xf_Error(xf_SolveBlockPLUT(Rii, ni, sr, P, Xi, xfe_Set, Xi));
    if (ierr != xf_OK) return ierr;

  } // iindex

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_SolveM_ILU0
static int
xf_Jacobian_SolveM_ILU0(xf_All *All, xf_JacobianMatrix *R_U, xf_Vector *X, 
			enum xfe_AddType AddFlag, enum xfe_Bool TransposeFlag)
{ 
  // X = M^{-1}*X 
  int ierr;

  // only support Set for now
  if (AddFlag != xfe_Set) return xf_Error(xf_NOT_SUPPORTED);

  if (!TransposeFlag){
    // lower solve
    ierr = xf_Error(xf_Jacobian_SolveM_ILU0_Lower(All, R_U, X));
    if (ierr != xf_OK) return ierr;
    // upper solve
    ierr = xf_Error(xf_Jacobian_SolveM_ILU0_Upper(All, R_U, X));
    if (ierr != xf_OK) return ierr;
  }
  else{
    // transposed upper solve
    ierr = xf_Error(xf_Jacobian_SolveM_ILU0_UpperT(All, R_U, X));
    if (ierr != xf_OK) return ierr;
    // transposed lower solve
    ierr = xf_Error(xf_Jacobian_SolveM_ILU0_LowerT(All, R_U, X));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultM_ILU0
static int
xf_Jacobian_MultM_ILU0(xf_All *All, xf_JacobianMatrix *R_U,
		       xf_Vector *X, enum xfe_AddType AddFlag, 
		       enum xfe_Bool TransposeFlag, xf_Vector *Y)
{  
  // Y @= M * X  (or M^T * X if TransposeFlag)
  int ierr;
  int sr, sr2, rmax;
  int iindex, jindex, negrp;
  int ni, nj;
  int iegrp, ielem, ifacej;
  int jegrp, jelem, jfacei;
  int *P;
  enum xfe_AddType AddFlag2;
  real *Rii, *Rij, *Rji;
  real *Xi, *Zi, *Xj, *Yi, *Zj;
  xf_LinearSolverILUData *ILUData;
  xf_Vector *Z;
  xf_Mesh *Mesh;

  // ILUData needs to exist
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // store desired quantities
  Mesh     = All->Mesh;
  negrp    = Mesh->nElemGroup;
  sr       = R_U->StateRank;
  sr2      = sr*sr;
  rmax     = ILUData->rmax;

  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);

  // find Z = a temporary vector
  ierr = xf_Error(xf_FindSimilarVector(All, Y, "ILU_Temp", xfe_True,
                                       xfe_True, NULL, &Z, NULL));
  if (ierr != xf_OK) return ierr;
  
  /*
    Normal operation:
         Zi  = Uij*Xj      (j >= i)
	 Zi += Lij*Zj      (j <  i)
	 Yi @= Zi
    Transpose operation:
         Zi  = (Lji)^T*Xj  (j >= i)
	 Zi  = (Uii)^T*Zi + (Uji)^T*Zj  (j <  i)
	 Yi @= Zi
  */

  /* First, Z = U*X  */

  // loop forward over all elements, i
  for (iindex=0; iindex<ILUData->nelemtot; iindex++){

    // obtain ordered element (iegrp, ielem)
    iegrp = ILUData->Index2Elem[iindex][0];
    ielem = ILUData->Index2Elem[iindex][1];
    
    Rii = R_U->Value[iegrp][ielem][0];       // on-diagonal block
    ni  = xf_Jacobian_n(R_U,iegrp,ielem);    // r/StateRank for iegrp,ielem
    P   = R_U->P[iegrp][ielem];              // permutation vector for ielem
    Xi  = X->GenArray[iegrp].rValue[ielem];  // X for ielem
    Zi  = Z->GenArray[iegrp].rValue[ielem];  // Z for ielem

    if (!TransposeFlag){
      // Zi = Rii*Xi */
      ierr = xf_Error(xf_BlockPLUMxV(Rii, ni, sr, P, Xi, xfe_Set, Zi));
      if (ierr != xf_OK) return ierr;
    }
    else{
      // Zi = Xi, to account for I on the diagonal in L
      xf_V_Add(Xi, ni*sr, xfe_Set, Zi);
    }
    
    // include upper-diagonal blocks via a loop over higher-ordered neighbors
    for (ifacej=0; ifacej<Mesh->ElemGroup[iegrp].nFace[ielem]; ifacej++){

      // obtain jegrp, jelem
      xf_JacobianNeighbor(R_U, iegrp, ielem, ifacej, &jegrp, &jelem, &jfacei);

      // skip halos and boundaries
      if ((jegrp<0) || (jegrp>=negrp)) continue;

      // ordered index of j
      jindex = ILUData->Elem2Index[jegrp][jelem];
      if (jindex <= iindex) continue;  // this neighbor is of lower ordering

      nj  = xf_Jacobian_n(R_U,jegrp,jelem); // r/StateRank for jegrp,jelem

      // set pointers
      Xj  = X->GenArray[jegrp].rValue[jelem];      // X for jelem
    
      if (!TransposeFlag){
	// Zi += Rij*Xj */
	Rij = R_U->Value[iegrp][ielem][ifacej + 1];  // off-diag Jacobian
	ierr = xf_Error(xf_BlockMxV(Rij, ni, sr, nj, Xj, xfe_Add, Zi));
	if (ierr != xf_OK) return ierr;
      }
      else{
	// Zi += (Rji)^T*Xj */
	Rji = R_U->Value[jegrp][jelem][jfacei + 1];  // off-diag Jacobian
	ierr = xf_Error(xf_BlockMTxV(Rji, ni, sr, nj, Xj, xfe_Add, Zi));
	if (ierr != xf_OK) return ierr;
      }
      
    } // ifacej
  } // iindex


  /* Second, Z = L*Z  */

  // loop backward over all elements, i
  for (iindex=ILUData->nelemtot-1; iindex>=0; iindex--){

    // obtain ordered element (iegrp, ielem)
    iegrp = ILUData->Index2Elem[iindex][0];
    ielem = ILUData->Index2Elem[iindex][1];
    
    ni  = xf_Jacobian_n(R_U,iegrp,ielem);    // r/StateRank for iegrp,ielem
    Yi  = Y->GenArray[iegrp].rValue[ielem];  // Y for ielem
    Zi  = Z->GenArray[iegrp].rValue[ielem];  // Z for ielem
    
    if (!TransposeFlag){
      // Yi @= Zi
      xf_V_Add(Zi, ni*sr, AddFlag, Yi);
    }
    else{
      // Yi @= Rii^T*Zi */
      Rii = R_U->Value[iegrp][ielem][0];       // on-diagonal block
      P   = R_U->P[iegrp][ielem];              // permutation vector for ielem
      ierr = xf_Error(xf_BlockPLUMTxV(Rii, ni, sr, P, Zi, AddFlag, Yi));
      if (ierr != xf_OK) return ierr;
    }
    
    // include lower-diagonal blocks via a loop over lower-ordered neighbors
    for (ifacej=0; ifacej<Mesh->ElemGroup[iegrp].nFace[ielem]; ifacej++){

      // obtain jegrp, jelem
      xf_JacobianNeighbor(R_U, iegrp, ielem, ifacej, &jegrp, &jelem, &jfacei);

      // skip halos and boundaries
      if ((jegrp<0) || (jegrp>=negrp)) continue;

      // ordered index of j
      jindex = ILUData->Elem2Index[jegrp][jelem];
      if (jindex >= iindex) continue;  // this neighbor is of higher ordering

      nj  = xf_Jacobian_n(R_U,jegrp,jelem); // r/StateRank for jegrp,jelem

      // set pointers
      Zj  = Z->GenArray[jegrp].rValue[jelem];      // Z for jelem
      
      if (!TransposeFlag){
	// Yi @= Rij*Zj */
	Rij = R_U->Value[iegrp][ielem][ifacej + 1];  // off-diag Jacobian
	ierr = xf_Error(xf_BlockMxV(Rij, ni, sr, nj, Zj, AddFlag2, Yi));
	if (ierr != xf_OK) return ierr;
      }
      else{
	// Yi @= (Rji)^T*Zj */
	Rji = R_U->Value[jegrp][jelem][jfacei + 1];  // off-diag Jacobian
	ierr = xf_Error(xf_BlockMTxV(Rji, ni, sr, nj, Zj, AddFlag2, Yi));
	if (ierr != xf_OK) return ierr;
      }
      
    } // ifacej
  } // iindex

  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Jacobian_MultN_ILU0
static int
xf_Jacobian_MultN_ILU0(xf_All *All, xf_JacobianMatrix *R_U,
		       xf_Vector *X, enum xfe_AddType AddFlag,
		       enum xfe_Bool TransposeFlag, xf_Vector *Y)
{  
  // Y @= N * X  (or N^T * X if TransposeFlag)
  int ierr, k, off, negrp;
  int sr, sr2, rmax, nfacemax;
  int iindex, jindex, kindex;
  int ni, nj, nk;
  int iegrp, ielem, ifacej, ifacek;
  int jegrp, jelem, jfacei, jfacek;
  int kegrp, kelem, kfacei, kfacej;
  int *iFlag = NULL;
  enum xfe_AddType AddFlag2, AddFlag2Neg;
  enum xfe_Bool HaveUXk;
  real *iUXk=NULL, *UXk=NULL;
  real *Rkj, *Rik, *Rjk, *Rki, *Xj, *Yi;
  xf_LinearSolverILUData *ILUData;
  xf_Mesh *Mesh;

  // ILUData needs to exist
  if ((ILUData = R_U->ILUData) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // store desired quantities
  Mesh     = All->Mesh;
  negrp    = Mesh->nElemGroup;
  sr       = R_U->StateRank;
  sr2      = sr*sr;
  rmax     = ILUData->rmax;
  nfacemax = ILUData->nfacemax;

  // zero out Y if add flag is Set or Neg
  if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
    ierr = xf_Error(xf_SetZeroVector(Y));
    if (ierr != xf_OK) return ierr;
  }


  // AddFlag2 is used for additional operations: always either add or sub
  AddFlag2 = xf_GetAddFlag2(AddFlag);
  // determine opposite operation of AddFlag2
  AddFlag2Neg = xf_GetAddFlagNeg(AddFlag2);
  

  // allocate temporary memory
  ierr = xf_Error(xf_Alloc( (void **)  &iUXk, rmax*nfacemax, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **)   &UXk, rmax         , sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &iFlag,      nfacemax, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  /* 
     Note, N*X = -L*U*X for all L*U entries that fall outside the
     sparsity pattern of A.

     Specifically,  
       Yi  @= Nij*Xj       becomes
       Yi !@= Lik*Ukj*Xj,  where i,j are not connected (Aij = 0)

     The operation !@ is the negation of @ (e.g. Add -> Sub)

     For the transpose operation, we have N^T*X = -U^T*L^T*X:
       Yi  @= (N^T)ij*Xj       becomes
       Yi !@= (Uki)^T*(Ljk)^T*Xj,  where i,j are not connected (Aij = 0) 
  */

  
  // Primary loop is over k (the intermediate index) 
  for (kindex=0; kindex<ILUData->nelemtot; kindex++){

    // obtain ordered element (kegrp, kelem)
    kegrp = ILUData->Index2Elem[kindex][0];
    kelem = ILUData->Index2Elem[kindex][1];

    nk = xf_Jacobian_n(R_U,kegrp,kelem); // r/StateRank for k

    // zero out iFlag
    for (k=0; k<nfacemax; k++) iFlag[k] = 0;
    
    // loop over higher-ordered neighbors, j
    for (kfacej=0; kfacej<Mesh->ElemGroup[kegrp].nFace[kelem]; kfacej++){
      
      // obtain jegrp, jelem
      xf_JacobianNeighbor(R_U, kegrp, kelem, kfacej, &jegrp, &jelem, &jfacek);

      // skip halos and boundaries
      if ((jegrp<0) || (jegrp>=negrp)) continue;

      // ordered index of j
      jindex = ILUData->Elem2Index[jegrp][jelem];
      if (jindex <= kindex) continue;  // skip as j is of lower ordering than k

      nj = xf_Jacobian_n(R_U,jegrp,jelem); // r/StateRank for j

      // set UXk flag
      HaveUXk = xfe_False;

      // another loop over higher-ordered neighbors of k, i
      for (kfacei=0; kfacei<Mesh->ElemGroup[kegrp].nFace[kelem]; kfacei++){
	
	// obtain iegrp, ielem
	xf_JacobianNeighbor(R_U, kegrp, kelem, kfacei, &iegrp, &ielem, NULL);

	// skip halos and boundaries
	if ((iegrp<0) || (iegrp>=negrp)) continue;

	// ordered index of i
	iindex = ILUData->Elem2Index[iegrp][ielem];
	if (iindex <= kindex) continue;  // skip as i is of lower ordering than k

	// make sure i and j are not connected
	ierr = xf_CommonFace(Mesh, iegrp, ielem, jegrp, jelem, &ifacej, &jfacei);
	if (ierr != xf_NOT_FOUND) continue; // if connected, skip (we only consider Aij = 0)

	ni = xf_Jacobian_n(R_U,iegrp,ielem); // r/StateRank for i
	
	if (!HaveUXk){
	  Xj  = X->GenArray[jegrp].rValue[jelem];      // X for jelem
	  if (!TransposeFlag){
	    // compute UXk = Ukj*Xj
	    Rkj = R_U->Value[kegrp][kelem][kfacej + 1];  // Jacobian off-diag block
	    ierr = xf_Error(xf_BlockMxV(Rkj, nk, sr, nj, Xj, xfe_Set, UXk));
	    if (ierr != xf_OK) return ierr;
	  }
	  else{
	    // compute UXk = (Ljk)^T*Xj
	    Rjk = R_U->Value[jegrp][jelem][jfacek + 1];  // Jacobian off-diag block
	    ierr = xf_Error(xf_BlockMTxV(Rjk, nk, sr, nj, Xj, xfe_Set, UXk));
	    if (ierr != xf_OK) return ierr;
	  }
	  HaveUXk = xfe_True;
	}

	// iUXk += UXk
	off = kfacei*rmax;                             // offset into iUXk
	if (!iFlag[kfacei])                            // zero out if first add
	  for (k=0; k<nk*sr; k++) iUXk[off+k] = 0.;     
	for (k=0; k<nk*sr; k++) iUXk[off+k] += UXk[k]; // actual addition

	// mark i as visited via kfacei
	iFlag[kfacei] = 1;

      } // kfacei
      
    } // kfacej

    // now account for Lik via loop over higher-ordered neighbors of k, i
    for (kfacei=0; kfacei<Mesh->ElemGroup[kegrp].nFace[kelem]; kfacei++){
      
      // skip if no contributions from this neighbor
      if (!iFlag[kfacei]) continue;
      
      // obtain iegrp, ielem, ifacek
      xf_JacobianNeighbor(R_U, kegrp, kelem, kfacei, &iegrp, &ielem, &ifacek);
      
      // error if i is a halo or boundary (iFlag should not have been set)
      if ((iegrp<0) || (iegrp>=negrp)) return xf_Error(xf_CODE_LOGIC_ERROR);

      // error if ordered index of i is lower than k (why was iFlag set, then?)
      iindex = ILUData->Elem2Index[iegrp][ielem];
      if (iindex <= kindex) return xf_Error(xf_CODE_LOGIC_ERROR);
	            
      ni  = xf_Jacobian_n(R_U,iegrp,ielem); // r/StateRank for i

      // set offset and pointer into Y
      off = kfacei*rmax;                       // offset into iUXk
      Yi  = Y->GenArray[iegrp].rValue[ielem];  // Y for ielem

      if (!TransposeFlag){
	Rik = R_U->Value[iegrp][ielem][ifacek + 1];  // off-diagonal block
	// Yi !@= Lik*iUXk
	ierr = xf_Error(xf_BlockMxV(Rik, ni, sr, nk, iUXk+off, AddFlag2Neg, Yi));
	if (ierr != xf_OK) return ierr;
      }
      else{
	Rki = R_U->Value[kegrp][kelem][kfacei + 1];  // off-diagonal block
	// Yi !@= (Uki^T)*iUXk
	ierr = xf_Error(xf_BlockMTxV(Rki, ni, sr, nk, iUXk+off, AddFlag2Neg, Yi));
	if (ierr != xf_OK) return ierr;
      }
          
    } // kfacei
  } // kindex

  // release memory
  xf_Release( (void *)   UXk);
  xf_Release( (void *)  iUXk);
  xf_Release( (void *) iFlag);

  return xf_OK;

}

