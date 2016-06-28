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
 FILE:  xf_MeshParallel.c
 
 This file contains parallelization and unparallelization functions
 for the Mesh data structure.  It is not compiled itself into an
 object.  Rather, it is included in xf_Mesh.c.
 
 */


#include "xf_MPI.h"
#include "xf_Partition.h" 

/*-----------------*/
/* Parallelization */
/*-----------------*/  



/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeElems
static int 
xf_ParallelizeElems( int **Elem2Proc, xf_Mesh *Mesh_Glob, xf_Mesh *Mesh)
{
  int ierr, i, nelemtot;
  int myRank, nProc, iProc, jProc, ibuf[4];
  int negrp, nelem, nfacemax, egrp, elem, face;
  int *itemp, **LocHaloList, *nElem, *nHalo;
  int *egrpN, *elemN, *ProcFlag, *HaloProcNum;
  int **ElemList, **HaloList, **ProcList;
  int **EL2G, **HL2G, *sindex, *slen;
  int *nSendElem, ***pSendElem, *nRecvElem;
  void *sbuf, *rbuf;
  xf_ElemGroup *EG, *HG;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  negrp = Mesh->nElemGroup;
  
  // allocate element groups, including halo
  ierr = xf_Error(xf_Alloc((void **) &Mesh->ElemGroup, 2*Mesh->nElemGroup,
                           sizeof(xf_ElemGroup)));
  if (ierr!=xf_OK) return ierr;
  
  // allocate arrays in ParallelInfo
  ierr = xf_Error(xf_Alloc((void **) &Mesh->ParallelInfo->ElemLoc2Glob,
                           2*Mesh->nElemGroup, sizeof(int *))); // halo gets one too
  if (ierr!=xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc2((void ***) &Mesh->ParallelInfo->nSendElem,
                            Mesh->nElemGroup, nProc, sizeof(int)));
  if (ierr!=xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void **) &Mesh->ParallelInfo->SendElem,
                           Mesh->nElemGroup, sizeof(int **)));
  if (ierr!=xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc2((void ***) &Mesh->ParallelInfo->nRecvElem,
                            Mesh->nElemGroup, nProc, sizeof(int)));
  if (ierr!=xf_OK) return ierr;
  
  // set pointers to NULL
  LocHaloList  = NULL;
  nElem        = NULL;
  nHalo        = NULL;
  egrpN        = NULL;
  elemN        = NULL;
  ProcFlag     = NULL;
  ElemList     = NULL;
  HaloList     = NULL;
  ProcList     = NULL;
  HaloProcNum  = NULL;
  
  // allocate proc-0-specific memory
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &nElem, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &nHalo, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &ProcFlag, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }
  else{
    // alloc dummy proc > 0 variables
    ierr = xf_Error(xf_Alloc2( (void ***) &ElemList, 1, 1, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc2( (void ***) &HaloList, 1, 1, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc2( (void ***) &ProcList, 1, 1, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
  }
  
  nelemtot = 0; // counter to ensure all procs get at least one elem
  
  for (egrp=0; egrp<negrp; egrp++){
    
    EG = Mesh->ElemGroup + egrp;         // local elem group
    HG = Mesh->ElemGroup + negrp + egrp; // local halo group
    
    xf_InitElemGroup(EG);
    xf_InitElemGroup(HG);
    
    // local-to-global numbering on egrp and halo egrp
    EL2G = Mesh->ParallelInfo->ElemLoc2Glob + egrp;
    HL2G = Mesh->ParallelInfo->ElemLoc2Glob + negrp + egrp;
    
    // pointers to send and receive lists
    nSendElem = Mesh->ParallelInfo->nSendElem[egrp];
    pSendElem = Mesh->ParallelInfo->SendElem +egrp;
    nRecvElem = Mesh->ParallelInfo->nRecvElem[egrp];
    
    // bcast basic elemgroup info
    if (myRank == 0){
      ibuf[0] = Mesh_Glob->ElemGroup[egrp].CutFlag;
      ibuf[1] = Mesh_Glob->ElemGroup[egrp].QBasis;
      ibuf[2] = Mesh_Glob->ElemGroup[egrp].QOrder;
      ibuf[3] = Mesh_Glob->ElemGroup[egrp].nNode;
    }
    
    ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 4*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    HG->CutFlag = EG->CutFlag = ibuf[0];
    HG->QBasis  = EG->QBasis  = ibuf[1];
    HG->QOrder  = EG->QOrder  = ibuf[2];
    HG->nNode   = EG->nNode   = ibuf[3];
    
    
    // proc 0 determines ElemList, HaloList, ProcList for each proc
    if (myRank == 0){
      nelem = Mesh_Glob->ElemGroup[egrp].nElem;
      
      for (elem=0,nfacemax=0; elem<nelem; elem++) 
        nfacemax = max(nfacemax, Mesh_Glob->ElemGroup[egrp].nFace[elem]);
      
      ierr = xf_Error(xf_ReAlloc((void **) &egrpN, nfacemax, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      
      ierr = xf_Error(xf_ReAlloc((void **) &elemN, nfacemax, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      
      // in first loop over elems[egrp], count nElem[iProc], nHalo[iProc]
      for (iProc=0; iProc<nProc; iProc++) nElem[iProc] = 0;
      for (iProc=0; iProc<nProc; iProc++) nHalo[iProc] = 0;
      for (elem=0; elem<nelem; elem++){
        iProc = Elem2Proc[egrp][elem];
        nElem[iProc]++;
        for (face=0; face<Mesh_Glob->ElemGroup[egrp].nFace[elem]; face++){
          ierr = xf_Error(xf_NeighborAcrossFace(Mesh_Glob, egrp, elem, face, 
                                                egrpN+face, elemN+face, NULL));
          if (ierr != xf_OK) return ierr;
          if (egrpN[face] < 0) continue;
          ProcFlag[Elem2Proc[egrpN[face]][elemN[face]]] = 0;
        }
        // for Halo, count only distinct procs != iProc, adjacent to elem
        for (face=0; face<Mesh_Glob->ElemGroup[egrp].nFace[elem]; face++){
          if (egrpN[face] < 0) continue;
          jProc = Elem2Proc[egrpN[face]][elemN[face]];
          if ((jProc != iProc) && (ProcFlag[jProc] == 0)) 
            nHalo[jProc]++;
          ProcFlag[jProc]++;
        }
      } // elem
      
      /*       for (iProc=0; iProc<nProc; iProc++) */
      /*       	xf_printf("On iProc = %d, nElem = %d, nHalo = %d.\n", */
      /*       		  iProc, nElem[iProc], nHalo[iProc]); */
      
      // variable re-alloc ElemList, HaloList, ProcList
      ierr = xf_Error(xf_VReAlloc2( (void ***) &ElemList, nProc, nElem, sizeof(int) ));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VReAlloc2( (void ***) &HaloList, nProc, nHalo, sizeof(int) ));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_VReAlloc2( (void ***) &ProcList, nProc, nHalo, sizeof(int) ));
      if (ierr != xf_OK) return ierr;
      
      // in second loop over elems[egrp], construct ElemList, HaloList
      for (iProc=0; iProc<nProc; iProc++) nElem[iProc] = 0;
      for (iProc=0; iProc<nProc; iProc++) nHalo[iProc] = 0;
      for (elem=0; elem<nelem; elem++){
        iProc = Elem2Proc[egrp][elem];
        ElemList[iProc][nElem[iProc]++] = elem;
        for (face=0; face<Mesh_Glob->ElemGroup[egrp].nFace[elem]; face++){
          ierr = xf_Error(xf_NeighborAcrossFace(Mesh_Glob, egrp, elem, face, 
                                                egrpN+face, elemN+face, NULL));
          if (ierr != xf_OK) return ierr;
          if (egrpN[face] < 0) continue;
          ProcFlag[Elem2Proc[egrpN[face]][elemN[face]]] = 0;
        }
        // add only distinct procs != iProc, adjacent to elem
        for (face=0; face<Mesh_Glob->ElemGroup[egrp].nFace[elem]; face++){
          if (egrpN[face] < 0) continue;
          jProc = Elem2Proc[egrpN[face]][elemN[face]];
          if ((jProc != iProc) && (ProcFlag[jProc] == 0)){
            HaloList[jProc][nHalo[jProc]] = elem;
            ProcList[jProc][nHalo[jProc]] = iProc;
            nHalo[jProc]++;
          }
          ProcFlag[jProc]++;
        }
      } // elem
      
      /*       for (iProc=0; iProc<nProc; iProc++){ */
      /*       	xf_printf("ElemList[iProc = %d] = ", iProc); */
      /*       	for (elem=0; elem<nElem[iProc]; elem++) */
      /*       	  xf_printf("%d ", ElemList[iProc][elem]); */
      /*       	xf_printf("\n"); */
      /*       	xf_printf("HaloList[iProc = %d] = ", iProc); */
      /*       	for (elem=0; elem<nHalo[iProc]; elem++) */
      /*       	  xf_printf("%d ", HaloList[iProc][elem]); */
      /*       	xf_printf("\n"); */
      /*       	xf_printf("ProcList[iProc = %d] = ", iProc); */
      /*       	for (elem=0; elem<nHalo[iProc]; elem++) */
      /*       	  xf_printf("%d ", ProcList[iProc][elem]); */
      /*       	xf_printf("\n"); */
      /*       } */
      /*       fflush(stdout); */
    }
    
    
    // scatter nElem
    ierr = xf_Error(xf_MPI_Scatter((void *) nElem, (void *) &EG->nElem, 
                                   1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // scatter nHalo
    ierr = xf_Error(xf_MPI_Scatter((void *) nHalo, (void *) &HG->nElem, 
                                   1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // alloc ElemLoc2Glob for egrp and halo
    ierr = xf_Error(xf_Alloc((void **) EL2G, EG->nElem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) HL2G, HG->nElem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    // re-alloc HaloProcNum
    ierr = xf_Error(xf_ReAlloc((void **) &HaloProcNum, HG->nElem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    
    // variable-length scatter ElemList, HaloList, ProcList
    
    ierr = xf_Error(xf_MPI_Scatterv((void *) ElemList[0], nElem, (void *) (*EL2G), 
                                    EG->nElem, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    sbuf = ((HaloList == NULL) ? NULL: (void *) HaloList[0]);
    ierr = xf_Error(xf_MPI_Scatterv(sbuf, nHalo, (void *) (*HL2G), 
                                    HG->nElem, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    sbuf = ((ProcList == NULL) ? NULL: (void *) ProcList[0]);
    ierr = xf_Error(xf_MPI_Scatterv(sbuf, nHalo, (void *) HaloProcNum, 
                                    HG->nElem, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // use HaloProcNum to sort HL2G, form nRecvElem[iProc]
    for (iProc=0; iProc<nProc; iProc++) nRecvElem[iProc] = 0;
    for (i=0; i<HG->nElem; i++) nRecvElem[HaloProcNum[i]]++;
    
    ierr = xf_Error(xf_VReAlloc2( (void ***) &LocHaloList, nProc, nRecvElem, 
                                 sizeof(int) ));
    if (ierr != xf_OK) return ierr;
    
    for (iProc=0; iProc<nProc; iProc++) nRecvElem[iProc] = 0;
    for (i=0; i<HG->nElem; i++){
      iProc = HaloProcNum[i];
      LocHaloList[iProc][nRecvElem[iProc]++] = (*HL2G)[i];
    }
    
    if (LocHaloList != NULL) swap(LocHaloList[0], (*HL2G), itemp);    
    
    /* Each proc sends a list of its halo elements to neighboring
     procs.  The neighboring procs receiving this list use it to
     form their SendElem list.  This is done in two steps: */
    
    /* first, use an all-to-all scalar communication: each proc1 tells
     every other proc2 how many elems proc2 has in proc1's halo group */
    ierr = xf_Error(xf_MPI_Alltoall((void *) nRecvElem, (void *) nSendElem, 
                                    1*sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    // variable-alloc2 SendElem
    ierr = xf_Error(xf_VAlloc2( (void ***) pSendElem, nProc, nSendElem, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
    
    /* second, use a variable-length all-to-all vector communication:
     each proc1 sends a list of its halo elements to every
     appropriate neighboring proc2 (i.e. every proc2 for which
     proc1's nRecvElem[proc2] > 0) */
    sbuf = (((*pSendElem) == NULL) ? NULL: (void *) (*pSendElem)[0]);
    ierr = xf_Error(xf_MPI_Alltoallv((void *) (*HL2G), nRecvElem, 
                                     sbuf, nSendElem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    /* variable-length gather HL2G into HaloList to inform proc 0 of
     the new halo ordering */
    sbuf = ((HaloList == NULL) ? NULL: (void *) HaloList[0]);
    ierr = xf_Error(xf_MPI_Gatherv((void *) (*HL2G), HG->nElem, sbuf, 
                                   nHalo, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    /*-------*/
    /* .Node */
    /*-------*/ 
    
    // allocate group and halo Node
    ierr = xf_Error(xf_Alloc2((void ***) &EG->Node, EG->nElem, 
                              EG->nNode, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc2((void ***) &HG->Node, HG->nElem, 
                              HG->nNode, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    
    // proc 0 packs and sends off Node lists to each Proc
    sbuf = ((myRank == 0) ? (void *) Mesh_Glob->ElemGroup[egrp].Node[0] : NULL);
    rbuf = ((EG->Node == NULL) ? NULL : (void *) EG->Node[0]);
    sindex = ((ElemList == NULL) ? NULL : ElemList[0]);
    ierr = xf_Error(xf_MPI_PScatterv(sbuf, sindex, nElem, rbuf, EG->nElem, 
                                     EG->nNode*sizeof(int), 0));
    if (ierr!=xf_OK) return ierr;
    rbuf = ((HG->Node == NULL) ? NULL : (void *) HG->Node[0]);
    sindex = ((HaloList == NULL) ? NULL : HaloList[0]);
    ierr = xf_Error(xf_MPI_PScatterv(sbuf, sindex, nHalo, rbuf, HG->nElem, 
                                     HG->nNode*sizeof(int), 0));
    if (ierr!=xf_OK) return ierr;
    
    /*-------*/
    /* .Face */
    /*-------*/ 
    
    // allocate group and Halo nFace
    ierr = xf_Error(xf_Alloc((void **) &EG->nFace, EG->nElem, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &HG->nFace, HG->nElem, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    
    // proc 0 packs and sends off nFace list to each proc (group + halo)
    sbuf = ((myRank == 0) ? (void *) Mesh_Glob->ElemGroup[egrp].nFace : NULL);
    sindex = ((ElemList == NULL) ? NULL : ElemList[0]);
    ierr = xf_Error(xf_MPI_PScatterv(sbuf, sindex, nElem, EG->nFace, EG->nElem, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    sindex = ((HaloList == NULL) ? NULL : HaloList[0]);
    ierr = xf_Error(xf_MPI_PScatterv(sbuf, sindex, nHalo, HG->nFace, HG->nElem, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // allocate group and Halo Face
    ierr = xf_Error(xf_VAlloc2((void ***) &EG->Face, EG->nElem, 
                               EG->nFace, sizeof(xf_Face)));
    if (ierr!=xf_OK) return ierr;
    ierr = xf_Error(xf_VAlloc2((void ***) &HG->Face, HG->nElem, 
                               HG->nFace, sizeof(xf_Face)));
    if (ierr!=xf_OK) return ierr;
    
    // proc 0 packs and sends off Face lists to each Proc
    sbuf = ((myRank == 0) ? (void *) Mesh_Glob->ElemGroup[egrp].Face[0] : NULL);
    slen = ((myRank == 0) ? (void *) Mesh_Glob->ElemGroup[egrp].nFace : NULL);
    rbuf = ((EG->Face == NULL) ? NULL : (void *) EG->Face[0]);
    sindex = ((ElemList == NULL) ? NULL : ElemList[0]);
    ierr = xf_Error(xf_MPI_DPScatterv(sbuf, slen, sindex, nElem, rbuf, EG->nFace,
                                      EG->nElem, sizeof(xf_Face), 0));
    if (ierr!=xf_OK) return ierr;
    rbuf = ((HG->Face == NULL) ? NULL : (void *) HG->Face[0]);
    sindex = ((HaloList == NULL) ? NULL : HaloList[0]);
    ierr = xf_Error(xf_MPI_DPScatterv(sbuf, slen, sindex, nHalo, rbuf, HG->nFace,
                                      HG->nElem, sizeof(xf_Face), 0));
    if (ierr!=xf_OK) return ierr;
    
    
    /* each proc nulls out its Halo group .Face list; the necessary
     iface connectivity will be set during renumbering. */
    for (elem=0; elem<HG->nElem; elem++)
      for (face=0; face<HG->nFace[elem]; face++){
        HG->Face[elem][face].Group  = xf_NULLFACE;
	HG->Face[elem][face].Number = -1;
      }
    
    nelemtot += EG->nElem;
    
  } // egrp
  
  
  // reduce-min nelemtot
  ierr = xf_Error(xf_MPI_Allreduce(&nelemtot, 1, xfe_SizeInt, xfe_MPI_MIN));
  if (ierr != xf_OK) return ierr;
  
  if (nelemtot <= 0){
    xf_printf("Error, one or more procs has zero elements following partitioning.\n");
    xf_printf("Reduce number of processors.\n");
    return xf_Error(xf_PARALLEL_ERROR);
  }
  
  
  
  xf_Release2((void **) LocHaloList);
  xf_Release2((void **) ElemList);
  xf_Release2((void **) HaloList);
  xf_Release2((void **) ProcList);
  xf_Release( (void  *) nElem);
  xf_Release( (void  *) nHalo);
  xf_Release( (void  *) ProcFlag);
  xf_Release( (void  *) HaloProcNum);
  xf_Release( (void  *) egrpN);
  xf_Release( (void  *) elemN);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeIFaces
static int 
xf_ParallelizeIFaces( int **Elem2Proc, xf_Mesh *Mesh_Glob, xf_Mesh *Mesh){
  
  int ierr, iiface;
  int myRank, nProc, iProc, iProcL, iProcR;
  int *nIFace, **IFaceList, *IFL2G, *ibuf;
  int *nIFNorm, *nIFHalo, **IFaceHalo;
  void *sbuf;
  xf_IFace *iface;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  nIFace    = NULL;
  nIFNorm   = NULL;
  IFaceList = NULL;
  nIFHalo   = NULL;
  IFaceHalo = NULL;
  nIFNorm   = NULL;
  
  // proc 0 determines IFaceList to send to each proc
  if (myRank == 0){
    
    // allocate IFace counter and normal (non-halo) IFace counter
    ierr = xf_Error(xf_Alloc( (void **) &nIFace, nProc, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &nIFNorm, nProc, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
    
    // zero out counters
    for (iProc=0; iProc<nProc; iProc++) nIFace[iProc]  = 0;
    for (iProc=0; iProc<nProc; iProc++) nIFNorm[iProc] = 0;
    
    for (iiface=0; iiface<Mesh_Glob->nIFace; iiface++){
      iface = Mesh_Glob->IFace + iiface;
      iProcL = Elem2Proc[iface->ElemGroupL][iface->ElemL];
      iProcR = Elem2Proc[iface->ElemGroupR][iface->ElemR];
      if (iProcL == iProcR){  // normal iface on iProcL = iProcR
        nIFNorm[iProcL]++;
        nIFace[iProcL]++;
      }
      else{  // halo iface between iProcL and iProcR
        nIFace[iProcL]++;
        nIFace[iProcR]++;
      }
    } // iiface
    
    // alloc IFaceList as a variable-length 2d array
    ierr = xf_Error(xf_VAlloc2( (void ***) &IFaceList, nProc, nIFace, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
    
    // temporary pointers into halo ifaces, since these will come last in IFaceList
    ierr = xf_Error(xf_Alloc( (void **) &IFaceHalo, nProc, sizeof(int *) ));
    if (ierr != xf_OK) return ierr;
    for (iProc=0; iProc<nProc; iProc++) IFaceHalo[iProc] = IFaceList[iProc]+nIFNorm[iProc];
    
    // allocate and zero out counter for halo ifaces
    ierr = xf_Error(xf_Alloc( (void **) &nIFHalo, nProc, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
    for (iProc=0; iProc<nProc; iProc++) nIFHalo[iProc] = 0;
    
    // reset normal (non-halo) IFace counter
    for (iProc=0; iProc<nProc; iProc++) nIFNorm[iProc] = 0;
    
    for (iiface=0; iiface<Mesh_Glob->nIFace; iiface++){
      iface = Mesh_Glob->IFace + iiface;
      iProcL = Elem2Proc[iface->ElemGroupL][iface->ElemL];
      iProcR = Elem2Proc[iface->ElemGroupR][iface->ElemR];
      if (iProcL == iProcR){  // normal iface on iProcL = iProcR
        IFaceList[iProcL][nIFNorm[iProcL]] = iiface;
        nIFNorm[iProcL]++;
      }
      else{  // halo iface between iProcL and iProcR
        IFaceHalo[iProcL][nIFHalo[iProcL]] = iiface;
        IFaceHalo[iProcR][nIFHalo[iProcR]] = iiface;
        nIFHalo[iProcL]++;
        nIFHalo[iProcR]++;
      }
    } // iiface
    
    // error check
    for (iProc=0; iProc<nProc; iProc++)
      if ( (nIFNorm[iProc] + nIFHalo[iProc]) != nIFace[iProc])
        return xf_Error(xf_PARALLEL_ERROR);
    
    // release no-longer-necessary memory
    xf_Release( (void *) nIFHalo);
    xf_Release( (void *) IFaceHalo);
    
  }
  
  // scatter nIFNorm
  ierr = xf_Error(xf_MPI_Scatter((void *) nIFNorm, 
                                 (void *) &Mesh->ParallelInfo->nIFaceRegular, 
                                 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // scatter nIFace
  ierr = xf_Error(xf_MPI_Scatter((void *) nIFace, (void *) &Mesh->nIFace, 
                                 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // alloc IFaceLoc2Glob
  ierr = xf_Error(xf_Alloc((void **) &Mesh->ParallelInfo->IFaceLoc2Glob, 
                           Mesh->nIFace, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;
  IFL2G = Mesh->ParallelInfo->IFaceLoc2Glob;
  
  // variable-length scatter IFaceList
  ibuf = ((myRank == 0) ? IFaceList[0] : NULL);
  ierr = xf_Error(xf_MPI_Scatterv((void *) ibuf, nIFace, (void *) IFL2G, 
                                  Mesh->nIFace, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // allocate Mesh->IFace
  ierr = xf_Error(xf_Alloc( (void **) &Mesh->IFace, Mesh->nIFace, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;
  
  // pack and scatter Mesh->IFace
  sbuf = ((myRank == 0) ? (void *) Mesh_Glob->IFace : NULL);
  ibuf = ((myRank == 0) ?              IFaceList[0] : NULL);
  ierr = xf_Error(xf_MPI_PScatterv(sbuf, ibuf, nIFace, (void *) Mesh->IFace,
                                   Mesh->nIFace,  sizeof(xf_IFace), 0));
  if (ierr!=xf_OK) return ierr;
  
  
  /* Release memory on proc 0 */
  if (myRank == 0){
    xf_Release( (void *) nIFace);
    xf_Release( (void *) nIFNorm);
    xf_Release2( (void **) IFaceList);
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeBFaces
static int 
xf_ParallelizeBFaces( int **Elem2Proc, xf_Mesh *Mesh_Glob, xf_Mesh *Mesh){
  
  int ierr, ibfgrp, nbfgrp, ibface;
  int myRank, iProc, nProc, len;
  int *nBFace, **BFaceList, *ibuf;
  void *sbuf;
  xf_BFace *bface;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  nbfgrp = Mesh->nBFaceGroup;
  
  ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup, Mesh->nBFaceGroup, 
                           sizeof(xf_BFaceGroup)));
  if (ierr != xf_OK) return ierr;
  
  nBFace    = NULL;
  BFaceList = NULL;
  
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc( (void **) &nBFace, nProc, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
  }
  
  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
    
    // Bcast bfgroup title
    if (myRank == 0) len = strlen(Mesh_Glob->BFaceGroup[ibfgrp].Title)+1;
    ierr = xf_Error(xf_MPI_Bcast((void *) &len, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc( (void **) &Mesh->BFaceGroup[ibfgrp].Title, 
                             len, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    
    if (myRank == 0)
      strncpy(Mesh->BFaceGroup[ibfgrp].Title, Mesh_Glob->BFaceGroup[ibfgrp].Title, len);
    ierr = xf_Error(xf_MPI_Bcast((void *) Mesh->BFaceGroup[ibfgrp].Title, 
                                 len*sizeof(char), 0));
    if (ierr != xf_OK) return ierr;
    
    // proc 0 determines BFaceList to send to each proc
    if (myRank == 0){
      for (iProc=0; iProc<nProc; iProc++) nBFace[iProc] = 0;
      for (ibface=0; ibface<Mesh_Glob->BFaceGroup[ibfgrp].nBFace; ibface++){
        bface = Mesh_Glob->BFaceGroup[ibfgrp].BFace+ibface;
        nBFace[Elem2Proc[bface->ElemGroup][bface->Elem]]++;
      } // ibface
      
      // variable-length realloc
      ierr = xf_Error(xf_VReAlloc2( (void ***) &BFaceList, nProc, nBFace, sizeof(int) ));
      if (ierr != xf_OK) return ierr;
      
      for (iProc=0; iProc<nProc; iProc++) nBFace[iProc] = 0;
      for (ibface=0; ibface<Mesh_Glob->BFaceGroup[ibfgrp].nBFace; ibface++){
        bface = Mesh_Glob->BFaceGroup[ibfgrp].BFace+ibface;
        iProc = Elem2Proc[bface->ElemGroup][bface->Elem];
        BFaceList[iProc][nBFace[iProc]] = ibface;
        nBFace[iProc]++;
      } // ibface
    }
    
    // scatter nBFace
    ierr = xf_Error(xf_MPI_Scatter((void *) nBFace, (void *) &Mesh->BFaceGroup[ibfgrp].nBFace, 
                                   1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // allocate BFaces
    ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup[ibfgrp].BFace, 
                             Mesh->BFaceGroup[ibfgrp].nBFace, sizeof(xf_BFace)));
    if (ierr != xf_OK)  return ierr;
    
    // pack and scatter BFaces
    sbuf = ((myRank == 0) ? (void *) Mesh_Glob->BFaceGroup[ibfgrp].BFace : NULL);
    ibuf = ((myRank == 0) ?                                 BFaceList[0] : NULL);
    ierr = xf_Error(xf_MPI_PScatterv(sbuf, ibuf, nBFace, 
                                     (void *) Mesh->BFaceGroup[ibfgrp].BFace,
                                     Mesh->BFaceGroup[ibfgrp].nBFace, sizeof(xf_BFace), 0));
    if (ierr!=xf_OK) return ierr;
    
  } // ibfgrp
  
  
  /* Release memory on proc 0 */
  if (myRank == 0){
    xf_Release( (void *) nBFace);
    xf_Release2( (void **) BFaceList);
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeNodes
static int 
xf_ParallelizeNodes( xf_Mesh *Mesh_Glob, xf_Mesh *Mesh){
  
  int ierr, i, egrp, elem;
  int myRank, nProc;
  int nglob, nnode, negrp;
  int *NL2G, *NodeFlag, *nNode, **NodeList, *ibuf;
  char s[500], sout[500];
  void *sbuf;
  xf_ElemGroup *EG, *HG;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  negrp = Mesh->nElemGroup;
  
  // glob node number has already been broadcast
  nglob = Mesh->nNode;
  
  // allocate and initialize NodeFlag
  ierr = xf_Error(xf_Alloc( (void **) &NodeFlag, nglob, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nglob; i++) NodeFlag[i] = 0;
  
  
  // each proc counts unique nodes
  for (egrp=0; egrp<negrp; egrp++){
    EG = Mesh->ElemGroup + egrp;         // local elem group
    HG = Mesh->ElemGroup + negrp + egrp; // local halo group
    for (elem=0; elem<EG->nElem; elem++)
      for (i=0; i<EG->nNode; i++) NodeFlag[EG->Node[elem][i]] = 1;
    for (elem=0; elem<HG->nElem; elem++)
      for (i=0; i<HG->nNode; i++) NodeFlag[HG->Node[elem][i]] = 1;
  } // egrp
  
  for (i=0, nnode=0; i<nglob; i++) nnode += NodeFlag[i];
  Mesh->nNode = nnode;  // set local number of nodes
  
  // allocate NodeLoc2Glob, Coord on each proc
  ierr = xf_Error(xf_Alloc((void **) &Mesh->ParallelInfo->NodeLoc2Glob, nnode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  NL2G = Mesh->ParallelInfo->NodeLoc2Glob;
  
  ierr = xf_Error(xf_Alloc2((void ***) &Mesh->Coord, nnode, Mesh->Dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // fill in NL2G
  for (i=0, nnode=0; i<nglob; i++)
    if (NodeFlag[i]){
      NodeFlag[i] = nnode;
      NL2G[nnode++] = i;
    }
  
  // allocate nNode on 0
  nNode = NULL;
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc( (void **) &nNode, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  // gather number of nodes to proc 0
  ierr = xf_Error(xf_MPI_Gather((void *) &nnode, (void *) nNode, 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // variable-allocate node list on proc 0
  NodeList = NULL;
  if (myRank == 0){
    ierr = xf_Error(xf_VAlloc2( (void ***) &NodeList, nProc, nNode, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  // gather node lists onto proc 0
  ibuf = ((myRank == 0) ?                  NodeList[0] : NULL);
  ierr = xf_Error(xf_MPI_Gatherv((void *) NL2G, nnode, (void *) ibuf, 
                                 nNode, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // variable pack-and-scatter the node coordinates
  sbuf = ((myRank == 0) ? (void *) Mesh_Glob->Coord[0] : NULL);
  ibuf = ((myRank == 0) ?                  NodeList[0] : NULL);
  ierr = xf_Error(xf_MPI_PScatterv(sbuf, ibuf, nNode, (void *) Mesh->Coord[0], 
                                   Mesh->nNode, Mesh->Dim*sizeof(real), 0));
  if (ierr!=xf_OK) return ierr;
  
  // renumber nodes
  for (egrp=0; egrp<negrp; egrp++){
    EG = Mesh->ElemGroup + egrp;         // local elem group
    HG = Mesh->ElemGroup + negrp + egrp; // local halo group
    for (elem=0; elem<EG->nElem; elem++)
      for (i=0; i<EG->nNode; i++) EG->Node[elem][i] = NodeFlag[EG->Node[elem][i]];
    for (elem=0; elem<HG->nElem; elem++)
      for (i=0; i<HG->nNode; i++) HG->Node[elem][i] = NodeFlag[HG->Node[elem][i]];
  } // egrp
  
  
  /*   xf_pprintf("negrp = %d\n", negrp); */
  /*   for (egrp=0; egrp<2*negrp; egrp++){ */
  /*     EG = Mesh->ElemGroup + egrp; */
  /*     xf_pprintf("egrp = %d: nelem = %d, nNode=%d, nFace=%d\n", */
  /* 	    egrp, EG->nElem, EG->nNode, EG->nFace); */
  /*     for (elem=0; elem<EG->nElem; elem++){ */
  /*       xf_pprintf("elem = %d\n", elem); */
  /*       sprintf(sout, "Nodes = "); */
  /*       for (i=0; i<EG->nNode; i++) { */
  /* 	sprintf( s, " %d", EG->Node[elem][i]); */
  /* 	strcat(sout, s); */
  /*       } */
  /*       xf_pprintf("%s\n", sout); */
  /*       sprintf(sout, "Faces = "); */
  /*       for (i=0; i<EG->nFace; i++){ */
  /* 	sprintf( s, " [%d,%d]", EG->Face[elem][i].Group, EG->Face[elem][i].Number); */
  /* 	strcat(sout,s); */
  /*       } */
  /*       xf_pprintf("%s\n", sout); */
  /*       xf_pprintf("\n"); */
  /*     } */
  /*   } */
  
  // release NodeFlag
  xf_Release ((void  *) NodeFlag);
  xf_Release ((void  *) nNode);
  xf_Release2((void **) NodeList);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeRenumber
static int 
xf_ParallelizeRenumber( xf_Mesh *Mesh_Glob, xf_Mesh *Mesh){
  
  int ierr, i;
  int myRank, nProc, iProc;
  int egrp, negrp, elem, nelem;
  int iiface, ibfgrp, ibface;
  int *nSendElem, **SendElem;
  xf_IFace *iface;
  xf_BFace *bface;
  int **ElemFlag;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // construct ElemFlag[egrp][oldelem] = newelem #
  negrp = Mesh->nElemGroup;
  ierr = xf_Error(xf_Alloc( (void **) &ElemFlag, negrp, sizeof(int *)));
  if (ierr != xf_OK) return ierr;
  
  for (egrp=0; egrp<negrp; egrp++){
    if (myRank == 0) nelem = Mesh_Glob->ElemGroup[egrp].nElem;
    ierr = xf_Error(xf_MPI_Bcast((void *) &nelem, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc( (void **) ElemFlag + egrp, nelem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    for (i=0; i<nelem; i++) ElemFlag[egrp][i] = 0;
    
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      ElemFlag[egrp][Mesh->ParallelInfo->ElemLoc2Glob[egrp][elem]] = 1+elem;
    
    for (elem=0; elem<Mesh->ElemGroup[negrp+egrp].nElem; elem++)
      ElemFlag[egrp][Mesh->ParallelInfo->ElemLoc2Glob[negrp+egrp][elem]] = -1-elem;
  }
  
  // renumber IFace-Element connectivity
  for (iiface=0; iiface<Mesh->nIFace; iiface++){
    iface = Mesh->IFace + iiface;
    elem = ElemFlag[iface->ElemGroupL][iface->ElemL];
    if (elem == 0) 
      return xf_Error(xf_PARALLEL_ERROR);
    else if (elem > 0) // regular group
      iface->ElemL =  elem-1;
    else{ // halo group
      iface->ElemGroupL += negrp;
      iface->ElemL = -elem-1;
    }
    Mesh->ElemGroup[iface->ElemGroupL].Face[iface->ElemL][iface->FaceL].Group  = xf_INTERIORFACE;
    Mesh->ElemGroup[iface->ElemGroupL].Face[iface->ElemL][iface->FaceL].Number = iiface;
    
    elem = ElemFlag[iface->ElemGroupR][iface->ElemR];
    if (elem == 0) 
      return xf_Error(xf_PARALLEL_ERROR);
    else if (elem > 0) // regular group
      iface->ElemR =  elem-1;
    else{ // halo group
      iface->ElemGroupR += negrp;
      iface->ElemR = -elem-1;
    }
    Mesh->ElemGroup[iface->ElemGroupR].Face[iface->ElemR][iface->FaceR].Group  = xf_INTERIORFACE;
    Mesh->ElemGroup[iface->ElemGroupR].Face[iface->ElemR][iface->FaceR].Number = iiface;    
  } // iiface
  
  // renumber BFace-Element connectivity
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      bface = Mesh->BFaceGroup[ibfgrp].BFace + ibface;
      elem = ElemFlag[bface->ElemGroup][bface->Elem];
      if (elem == 0) 
        return xf_Error(xf_PARALLEL_ERROR);
      else if (elem > 0) // regular group
        bface->Elem =  elem-1;
      else{ // halo group
        return xf_Error(xf_CODE_LOGIC_ERROR); // should not have sent this bface
        bface->ElemGroup += negrp;
        bface->Elem = -elem-1;
      }
      Mesh->ElemGroup[bface->ElemGroup].Face[bface->Elem][bface->Face].Group  = ibfgrp;
      Mesh->ElemGroup[bface->ElemGroup].Face[bface->Elem][bface->Face].Number = ibface;
    } // ibface
  } // ibfgrp
  
  
  // renumber SendElem
  for (egrp=0; egrp<negrp; egrp++){
    nSendElem = Mesh->ParallelInfo->nSendElem[egrp];
    SendElem  = Mesh->ParallelInfo->SendElem[egrp];
    for (iProc=0; iProc<nProc; iProc++){
      for (i=0; i<nSendElem[iProc]; i++){
        elem = ElemFlag[egrp][SendElem[iProc][i]];
        if (elem <= 0) return xf_Error(xf_CODE_LOGIC_ERROR); // halos are not for sending
        SendElem[iProc][i] = elem-1;
      }
    }
  }
  
  // release ElemFlag
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    xf_Release((void *) ElemFlag[egrp]);
  xf_Release((void *) ElemFlag);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ParallelizeMesh
int 
xf_ParallelizeMesh(xf_Mesh *Mesh_Glob, xf_Mesh *Mesh, int **ElemWeight, 
                   int *ConnectWeight){
  
  int ierr;
  int myRank, nProc;
  int ibuf[10];
  int **Elem2Proc;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc == 1){ 
    /* This serial behavior is useful.  Mesh inherits all of
     Mesh_Glob's data, and Mesh_Glob has its data initalized to
     null. */
    (*Mesh) = (*Mesh_Glob);
    xf_InitMesh(Mesh_Glob);
    return xf_OK;
  }
  
  /* Determine on which proc each elem should reside.
   Elem2Proc[egrp][elem] = proc. */
  Elem2Proc = NULL;
  if (myRank == 0){
    xf_printf("Calling PartitionMesh ... ");
    ierr = xf_Error(xf_PartitionMesh(Mesh_Glob, &Elem2Proc, ElemWeight, 
                                     ConnectWeight));
    if (ierr != xf_OK) return ierr;
    xf_printf(" done.\n");
  }
  
  /*-----------------*/
  /* Mesh basic info */
  /*-----------------*/
  
  if (myRank == 0){
    ibuf[0] = Mesh_Glob->Dim;
    ibuf[1] = Mesh_Glob->nNode;
    ibuf[2] = Mesh_Glob->nIFace;
    ibuf[3] = Mesh_Glob->nBFaceGroup;
    ibuf[4] = Mesh_Glob->nElemGroup;
  }
  
  ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 5*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  Mesh->Dim         = ibuf[0];
  Mesh->nNode       = ibuf[1];
  Mesh->nIFace      = ibuf[2];
  Mesh->nBFaceGroup = ibuf[3];
  Mesh->nElemGroup  = ibuf[4];
  
  // create parallelinfo
  ierr = xf_Error(xf_CreateMeshParallelInfo( &Mesh->ParallelInfo));
  if (ierr != xf_OK) return ierr;
  
  /*-----------------------------*/
  /* Element Groups and Elements */
  /*-----------------------------*/
  ierr = xf_Error(xf_ParallelizeElems(Elem2Proc, Mesh_Glob, Mesh));
  if (ierr != xf_OK) return ierr;
  
  /*--------*/
  /* IFaces */
  /*--------*/
  ierr = xf_Error(xf_ParallelizeIFaces(Elem2Proc, Mesh_Glob, Mesh));
  if (ierr != xf_OK) return ierr;
  
  /*--------*/
  /* BFaces */
  /*--------*/
  ierr = xf_Error(xf_ParallelizeBFaces(Elem2Proc, Mesh_Glob, Mesh));
  if (ierr != xf_OK) return ierr;
  
  /*------------------*/
  /* Node Coordinates */
  /*------------------*/
  ierr = xf_Error(xf_ParallelizeNodes(Mesh_Glob, Mesh));
  if (ierr != xf_OK) return ierr;
  
  /*--------------------------------*/
  /* Renumber elems, ifaces, bfaces */
  /*--------------------------------*/
  ierr = xf_Error(xf_ParallelizeRenumber(Mesh_Glob, Mesh));
  if (ierr != xf_OK) return ierr;

  /*-------------*/
  /* Mesh Motion */
  /*-------------*/
  if (myRank == 0) ibuf[0] = (Mesh_Glob->Motion != NULL);
  ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  if (ibuf[0]){
    if (myRank == 0){
      Mesh->Motion = Mesh_Glob->Motion;
      Mesh_Glob->Motion = NULL;
    }
    else{
      ierr = xf_Error(xf_CreateMeshMotion(&Mesh->Motion));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_ParallelizeMeshMotion(Mesh->Motion));
    if (ierr != xf_OK) return ierr;
  }

  /*-----------------*/
  /* Periodic Groups */
  /*-----------------*/
  if (myRank == 0){ // root gets all info
    Mesh->nPeriodicGroup = Mesh_Glob->nPeriodicGroup;
    Mesh->PeriodicGroup  = Mesh_Glob->PeriodicGroup;
    Mesh_Glob->nPeriodicGroup = 0;
    Mesh_Glob->PeriodicGroup  = NULL;
  }
  
  /* Release Elem2Proc*/
  xf_Release2( (void **) Elem2Proc);
  
  return xf_OK;
}



/*-------------------*/
/* UnParallelization */
/*-------------------*/  


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeElems
static int 
xf_UnParallelizeElems( xf_Mesh *Mesh, xf_Mesh *Mesh_Glob){
  
  int ierr, i;
  int myRank, nProc;
  int egrp, elem, *nElem, **ElemList;
  int *ibuf, *NL2G, *rlen, **NodeLocal;
  xf_ElemGroup *EG, *EG_Glob;
  void *rbuf, *sbuf;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // set pointers to NULL
  nElem     = NULL;
  ElemList  = NULL;
  NodeLocal = NULL;
  
  // allocate element groups for Mesh_Glob
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &Mesh_Glob->ElemGroup, Mesh->nElemGroup,
                             sizeof(xf_ElemGroup)));
    if (ierr!=xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc((void **) &nElem, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }
  
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    EG = Mesh->ElemGroup + egrp;         // local elem group
    
    // initialize element group on Mesh_Glob, and set basic info
    if (myRank == 0){
      EG_Glob = Mesh_Glob->ElemGroup + egrp;
      xf_InitElemGroup(EG_Glob);
      
      EG_Glob->CutFlag = EG->CutFlag;
      EG_Glob->QBasis  = EG->QBasis;
      EG_Glob->QOrder  = EG->QOrder;
      EG_Glob->nNode   = EG->nNode;
    }
    
    // Reduce-sum number of elements onto proc 0
    rbuf = ( (myRank == 0) ? (void *) &EG_Glob->nElem : NULL);
    ierr = xf_Error(xf_MPI_Reduce(&EG->nElem, rbuf, 1, xfe_SizeInt, xfe_MPI_SUM, 0));
    if (ierr != xf_OK) return ierr;
    
    // allocate nFace and Node on Mesh_Glob
    if (myRank == 0){
      ierr = xf_Error(xf_Alloc((void **) &EG_Glob->nFace, EG_Glob->nElem, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      
      ierr = xf_Error(xf_Alloc2((void ***) &EG_Glob->Node, EG_Glob->nElem, 
                                EG_Glob->nNode, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
    }
    
    // each proc saves its local .Node numbering
    ierr = xf_Error(xf_ReAlloc2( (void ***) &NodeLocal, EG->nElem, EG->nNode, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (elem=0; elem<EG->nElem; elem++)
      for (i=0; i<EG->nNode; i++) NodeLocal[elem][i] = EG->Node[elem][i];
    
    // each proc renumbers its .Node back to global numbering (.Face done later)
    NL2G  = Mesh->ParallelInfo->NodeLoc2Glob;
    for (elem=0; elem<EG->nElem; elem++)
      for (i=0; i<EG->nNode; i++) EG->Node[elem][i] = NL2G[EG->Node[elem][i]];
    
    // gather number of elems to inform proc 0 of # on each proc
    ierr = xf_Error(xf_MPI_Gather((void *) &EG->nElem, (void *) nElem, 1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // proc 0 re-allocates ElemList
    if (myRank == 0){
      ierr = xf_Error(xf_VReAlloc2( (void ***) &ElemList, nProc, nElem, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    
    // variable gather global element numbers onto proc 0
    rbuf = ((myRank == 0) ? (void *) ElemList[0] : NULL);
    ierr = xf_Error(xf_MPI_Gatherv((void *) Mesh->ParallelInfo->ElemLoc2Glob[egrp], 
                                   EG->nElem, rbuf, nElem, sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // gather and unpack .Node
    ibuf = ((myRank == 0) ? ElemList[0] : NULL);
    rbuf = ((myRank == 0) ? (void *) Mesh_Glob->ElemGroup[egrp].Node[0] : NULL);
    sbuf = ((EG->Node == NULL) ? NULL : (void *) EG->Node[0]);
    ierr = xf_Error(xf_MPI_PGatherv(sbuf, EG->nElem, rbuf, ibuf, nElem, 
                                    EG->nNode*sizeof(int), 0));
    if (ierr!=xf_OK) return ierr;
    
    // gather and unpack .nFace
    rbuf = ((myRank == 0) ? (void *) Mesh_Glob->ElemGroup[egrp].nFace : NULL);
    ierr = xf_Error(xf_MPI_PGatherv(EG->nFace, EG->nElem, rbuf, ibuf, nElem, 
                                    sizeof(int), 0));
    if (ierr!=xf_OK) return ierr;
    
    // proc 0 allocates .Face on EG_Glob = Mesh_Glob->ElemGroup[egrp]
    if (myRank == 0){
      ierr = xf_Error(xf_VAlloc2((void ***) &EG_Glob->Face, EG_Glob->nElem,
                                 EG_Glob->nFace, sizeof(xf_Face)));
      if (ierr!=xf_OK) return ierr;
    }
    
    // gather and unpack .Face
    rbuf = ((myRank == 0) ? (void *) Mesh_Glob->ElemGroup[egrp].Face[0] : NULL);
    rlen = ((myRank == 0) ? (void *) Mesh_Glob->ElemGroup[egrp].nFace : NULL);
    sbuf = ((EG->Face == NULL) ? NULL : (void *) EG->Face[0]);
    ierr = xf_Error(xf_MPI_DPGatherv(sbuf, EG->nFace, EG->nElem, rbuf, 
                                     rlen, ibuf, nElem, sizeof(xf_Face), 0));
    if (ierr!=xf_OK) return ierr;
    
    // each proc restores its local .Node numbering
    for (elem=0; elem<EG->nElem; elem++)
      for (i=0; i<EG->nNode; i++) EG->Node[elem][i] = NodeLocal[elem][i];
  } // egrp
  
  // release memory
  xf_Release( (void  *) nElem);
  xf_Release2((void **) ElemList);
  xf_Release2((void **) NodeLocal);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeIFaces
static int 
xf_UnParallelizeIFaces( xf_Mesh *Mesh, xf_Mesh *Mesh_Glob){
  int ierr, i, iiface, fmax;
  int myRank, nProc;
  int *nIFace, **IFaceList, *ibuf;
  int **EL2G, negrp;
  void *rbuf;
  xf_IFace *iface;
  xf_IFace *IFaceLocal;
  xf_Face *Face;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // set pointers to null
  nIFace    = NULL;
  IFaceList = NULL;
  
  // alloc memory on proc 0
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &nIFace, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }
  
  // determine max glob iface+1 # on each proc
  fmax = 0;
  for (i=0; i<Mesh->nIFace; i++)
    fmax = max(fmax, Mesh->ParallelInfo->IFaceLoc2Glob[i]+1);
  
  // reduce-max onto proc 0 to obtain number of global ifaces
  rbuf = ( (myRank == 0) ? (void *) &Mesh_Glob->nIFace : NULL);
  ierr = xf_Error(xf_MPI_Reduce(&fmax, rbuf, 1, xfe_SizeInt, xfe_MPI_MAX, 0));
  if (ierr != xf_OK) return ierr;
  
  // allocate Faces on Mesh_Glob
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc( (void **) &Mesh_Glob->IFace, Mesh_Glob->nIFace, 
                             sizeof(xf_IFace)));
    if (ierr != xf_OK) return ierr;
  }
  
  // each proc saves its local IFaces
  ierr = xf_Error(xf_Alloc( (void **) &IFaceLocal, Mesh->nIFace, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;
  for (iiface=0; iiface<Mesh->nIFace; iiface++)
    IFaceLocal[iiface] = Mesh->IFace[iiface];
  
  // each proc sets the elem #s of its ifaces to their global values
  EL2G  = Mesh->ParallelInfo->ElemLoc2Glob;
  negrp = Mesh->nElemGroup;
  
  for (iiface=0; iiface<Mesh->nIFace; iiface++){
    iface = Mesh->IFace + iiface;
    iface->ElemL = EL2G[iface->ElemGroupL][iface->ElemL];
    iface->ElemR = EL2G[iface->ElemGroupR][iface->ElemR];
    if (iface->ElemGroupL >= negrp) iface->ElemGroupL -= negrp; // halo back to reg
    if (iface->ElemGroupR >= negrp) iface->ElemGroupR -= negrp; // halo back to reg
  } // iiface
  
  
  // gather number of ifaces to inform proc 0 of # on each proc
  ierr = xf_Error(xf_MPI_Gather((void *) &Mesh->nIFace, (void *) nIFace, 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // proc 0 allocates FaceList
  if (myRank == 0){
    ierr = xf_Error(xf_VAlloc2( (void ***) &IFaceList, nProc, nIFace, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  // variable gather global iface numbers onto proc 0
  rbuf = ((myRank == 0) ? (void *) IFaceList[0] : NULL);
  ierr = xf_Error(xf_MPI_Gatherv((void *) Mesh->ParallelInfo->IFaceLoc2Glob, 
                                 Mesh->nIFace, rbuf, nIFace, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // gather and unpack faces
  ibuf = ((myRank == 0) ? IFaceList[0] : NULL);
  rbuf = ((myRank == 0) ? (void *) Mesh_Glob->IFace : NULL);
  ierr = xf_Error(xf_MPI_PGatherv((void *) Mesh->IFace, Mesh->nIFace, rbuf, 
                                  ibuf, nIFace, sizeof(xf_IFace), 0));
  if (ierr!=xf_OK) return ierr;
  
  // loop over Mesh_Glob ifaces and make sure elems point to ifaces correctly
  if (myRank == 0){
    for (iiface=0; iiface<Mesh_Glob->nIFace; iiface++){
      iface = Mesh_Glob->IFace + iiface;
      Face = Mesh_Glob->ElemGroup[iface->ElemGroupL].Face[iface->ElemL] + iface->FaceL;
      if (Face->Group != xf_INTERIORFACE) return xf_Error(xf_PARALLEL_ERROR);
      Face->Number = iiface;
      Face = Mesh_Glob->ElemGroup[iface->ElemGroupR].Face[iface->ElemR] + iface->FaceR;
      if (Face->Group != xf_INTERIORFACE) return xf_Error(xf_PARALLEL_ERROR);
      Face->Number = iiface;
    } // iiface
  }
  
  // each proc restores its local IFaces
  for (iiface=0; iiface<Mesh->nIFace; iiface++)
    Mesh->IFace[iiface] = IFaceLocal[iiface];
  
  // release memory
  xf_Release( (void  *) nIFace);
  xf_Release2((void **) IFaceList);
  xf_Release( (void  *) IFaceLocal);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeBFaces
static int 
xf_UnParallelizeBFaces( xf_Mesh *Mesh, xf_Mesh *Mesh_Glob){
  
  int ierr, ibfgrp, ibface;
  int myRank, nProc, len;
  int *nBFace, *ibuf, **EL2G;
  void *rbuf;
  xf_BFace *bface;
  xf_BFace *BFaceLocal;
  xf_Face *Face;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  nBFace = NULL;
  BFaceLocal = NULL;
  
  // alloc memory on proc 0
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &Mesh_Glob->BFaceGroup, Mesh_Glob->nBFaceGroup, 
                             sizeof(xf_BFaceGroup)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc( (void **) &nBFace, nProc, sizeof(int) ));
    if (ierr != xf_OK) return ierr;
  }
  
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    
    // copy bfgroup title from proc 0
    if (myRank == 0) {
      len = strlen(Mesh->BFaceGroup[ibfgrp].Title)+1;
      ierr = xf_Error(xf_Alloc( (void **) &Mesh_Glob->BFaceGroup[ibfgrp].Title, 
                               len, sizeof(char)));
      if (ierr != xf_OK) return ierr;
      
      strncpy(Mesh_Glob->BFaceGroup[ibfgrp].Title, Mesh->BFaceGroup[ibfgrp].Title, len);
    }
    
    // each proc saves its local BFaces
    ierr = xf_Error(xf_ReAlloc( (void **) &BFaceLocal, Mesh->BFaceGroup[ibfgrp].nBFace, sizeof(xf_BFace)));
    if (ierr != xf_OK) return ierr;
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++)
      BFaceLocal[ibface] = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
    
    // each proc sets the elem #s of its bfaces to their global values
    EL2G = Mesh->ParallelInfo->ElemLoc2Glob;
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      bface = Mesh->BFaceGroup[ibfgrp].BFace+ibface;
      bface->Elem = EL2G[bface->ElemGroup][bface->Elem];
      if (bface->ElemGroup >= Mesh->nElemGroup)
        return xf_Error(xf_PARALLEL_ERROR);
    }
    
    // reduce-sum number of bfaces onto proc 0
    rbuf = ( (myRank == 0) ? (void *) &Mesh_Glob->BFaceGroup[ibfgrp].nBFace : NULL);
    ierr = xf_Error(xf_MPI_Reduce(&Mesh->BFaceGroup[ibfgrp].nBFace, 
                                  rbuf, 1, xfe_SizeInt, xfe_MPI_SUM, 0));
    if (ierr != xf_OK) return ierr;
    
    // allocate bfaces on proc 0
    if (myRank == 0){
      ierr = xf_Error(xf_Alloc((void **) &Mesh_Glob->BFaceGroup[ibfgrp].BFace, 
                               Mesh_Glob->BFaceGroup[ibfgrp].nBFace, sizeof(xf_BFace)));
      if (ierr != xf_OK)  return ierr;
    }
    
    // gather number of bfaces to inform proc 0 of # on each proc
    ierr = xf_Error(xf_MPI_Gather((void *) &Mesh->BFaceGroup[ibfgrp].nBFace, 
                                  (void *) nBFace, 1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // variable-gather bfaces (order will be based on proc)
    rbuf = ((myRank == 0) ? (void *) Mesh_Glob->BFaceGroup[ibfgrp].BFace : NULL);
    ierr = xf_Error(xf_MPI_Gatherv((void *) Mesh->BFaceGroup[ibfgrp].BFace, 
                                   Mesh->BFaceGroup[ibfgrp].nBFace, rbuf,
                                   nBFace, sizeof(xf_BFace), 0));
    if (ierr!=xf_OK) return ierr;
    
    // loop over Mesh_Glob bfaces and make sure elems point to bfaces correctly
    if (myRank == 0){
      for (ibface=0; ibface<Mesh_Glob->BFaceGroup[ibfgrp].nBFace; ibface++){
        bface = Mesh_Glob->BFaceGroup[ibfgrp].BFace+ibface;
        Face = Mesh_Glob->ElemGroup[bface->ElemGroup].Face[bface->Elem] + bface->Face;
        if (Face->Group != ibfgrp) return xf_Error(xf_PARALLEL_ERROR);
        Face->Number = ibface;
      }
    }
    
    // each proc restores its bfaces
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++)
      Mesh->BFaceGroup[ibfgrp].BFace[ibface] = BFaceLocal[ibface];
    
  } // ibfgrp
  
  // release memory
  xf_Release( (void  *) nBFace);
  xf_Release( (void  *) BFaceLocal);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeNodes
static int 
xf_UnParallelizeNodes( xf_Mesh *Mesh, xf_Mesh *Mesh_Glob){
  int ierr, i;
  int myRank, nProc;
  int nmax, *ibuf, *nNode, **NodeList;
  void *rbuf;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  // set pointers to null
  nNode    = NULL;
  NodeList = NULL;
  
  // alloc memory on proc 0
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc((void **) &nNode, nProc, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }
  
  // determine max glob node # + 1 on each proc
  nmax = 0;
  for (i=0; i<Mesh->nNode; i++)
    nmax = max(nmax, Mesh->ParallelInfo->NodeLoc2Glob[i]+1);
  
  // reduce-max onto proc 0 to obtain number of global mesh nodes
  rbuf = ( (myRank == 0) ? (void *) &Mesh_Glob->nNode : NULL);
  ierr = xf_Error(xf_MPI_Reduce(&nmax, rbuf, 1, xfe_SizeInt, xfe_MPI_MAX, 0));
  if (ierr != xf_OK) return ierr;
  
  // allocate Node coordinates on Mesh_Glob
  if (myRank == 0){
    ierr = xf_Error(xf_Alloc2((void ***) &Mesh_Glob->Coord, Mesh_Glob->nNode, 
                              Mesh_Glob->Dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }
  
  // gather number of nodes to inform proc 0 of # on each proc
  ierr = xf_Error(xf_MPI_Gather((void *) &Mesh->nNode, (void *) nNode, 1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // proc 0 allocates NodeList
  if (myRank == 0){
    ierr = xf_Error(xf_VAlloc2( (void ***) &NodeList, nProc, nNode, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  // variable gather global node numbers onto proc 0
  rbuf = ((myRank == 0) ? (void *) NodeList[0] : NULL);
  ierr = xf_Error(xf_MPI_Gatherv((void *) Mesh->ParallelInfo->NodeLoc2Glob, 
                                 Mesh->nNode, rbuf, nNode, sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  // gather and unpack node coordinates
  ibuf = ((myRank == 0) ? NodeList[0] : NULL);
  rbuf = ((myRank == 0) ? (void *) Mesh_Glob->Coord[0] : NULL);
  ierr = xf_Error(xf_MPI_PGatherv((void *) Mesh->Coord[0], Mesh->nNode, rbuf, 
                                  ibuf, nNode, Mesh->Dim*sizeof(real), 0));
  if (ierr!=xf_OK) return ierr;
  
  // release memory
  xf_Release( (void  *) nNode);
  xf_Release2((void **) NodeList);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeMesh
int 
xf_UnParallelizeMesh( xf_Mesh *Mesh, xf_Mesh *Mesh_Glob){
  
  int ierr;
  int myRank, nProc;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  /*-----------------*/
  /* Mesh basic info */
  /*-----------------*/
  
  if (myRank == 0){
    Mesh_Glob->Dim         = Mesh->Dim;
    Mesh_Glob->nNode       = Mesh->nNode;
    Mesh_Glob->nIFace      = Mesh->nIFace;
    Mesh_Glob->nBFaceGroup = Mesh->nBFaceGroup;
    Mesh_Glob->nElemGroup  = Mesh->nElemGroup;
  }
  
  /*-----------------------------*/
  /* Element Groups and Elements */
  /*-----------------------------*/
  ierr = xf_Error(xf_UnParallelizeElems(Mesh, Mesh_Glob));
  if (ierr != xf_OK) return ierr;
  
  /*--------*/
  /* IFaces */
  /*--------*/
  ierr = xf_Error(xf_UnParallelizeIFaces(Mesh, Mesh_Glob));
  if (ierr != xf_OK) return ierr;
  
  /*--------*/
  /* BFaces */
  /*--------*/
  ierr = xf_Error(xf_UnParallelizeBFaces(Mesh, Mesh_Glob));
  if (ierr != xf_OK) return ierr;
  
  /*------------------*/
  /* Node Coordinates */
  /*------------------*/
  ierr = xf_Error(xf_UnParallelizeNodes(Mesh, Mesh_Glob));
  if (ierr != xf_OK) return ierr;

  /*-------------*/
  /* Mesh Motion */
  /*-------------*/
  if (myRank == 0){
    Mesh_Glob->Motion = Mesh->Motion;
    Mesh->Motion = NULL;
  }
  
  /*-----------------*/
  /* Periodic groups */
  /*-----------------*/
  if (myRank == 0){
    Mesh_Glob->nPeriodicGroup = Mesh->nPeriodicGroup;
    Mesh_Glob->PeriodicGroup  = Mesh->PeriodicGroup;
    Mesh->nPeriodicGroup = 0;
    Mesh->PeriodicGroup  = NULL;
  }
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_BcastMesh
int 
xf_BcastMesh(xf_Mesh **pMesh)
{
  //Not broadcasting cut-elements and periodicgroup information
  int ierr, myRank, nProc, ibuf[10], size, i, j, k, your_turn;
  int *sbuf;
  char string[xf_MAXSTRLEN], filename[xf_MAXSTRLEN];
  xf_Mesh *Mesh;
  xf_Face FaceBuf;
  FILE *fid;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (nProc > 1){
    xf_printf("Broadcasting mesh ...");
    
    //Create mesh on all processors except root
    if (myRank != 0) {
      ierr = xf_Error(xf_CreateMesh(pMesh));
      if (ierr != xf_OK) return ierr;
    }
    Mesh = (*pMesh);
    
    //Get mesh basic info from root and broadcast
    if (myRank == 0){
      ibuf[0] = Mesh->Dim;
      ibuf[1] = Mesh->nNode;
      ibuf[2] = Mesh->nIFace;
      ibuf[3] = Mesh->nBFaceGroup;
      ibuf[4] = Mesh->nElemGroup;
    }
    
    ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 5*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    Mesh->Dim         = ibuf[0];
    Mesh->nNode       = ibuf[1];
    Mesh->nIFace      = ibuf[2];
    Mesh->nBFaceGroup = ibuf[3];
    Mesh->nElemGroup  = ibuf[4];
    
    //node coordinates
    if (myRank != 0) {
      ierr = xf_Error(xf_Alloc2((void ***)&Mesh->Coord, Mesh->nNode, 
                                Mesh->Dim, sizeof(real)));
      if (ierr != xf_OK) return ierr;
    }
    
    ierr = xf_Error(xf_MPI_Bcast((void *)Mesh->Coord[0], 
                                 Mesh->nNode*Mesh->Dim*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    
    //interior faces
    if (myRank != 0) {
      ierr = xf_Error(xf_Alloc((void **)&Mesh->IFace, Mesh->nIFace, 
                               sizeof(xf_IFace)));
      if (ierr != xf_OK) return ierr;
    }
    
    for (i = 0; i < Mesh->nIFace; i++){
      if (myRank == 0){
        ibuf[0] = Mesh->IFace[i].ElemGroupL;
        ibuf[1] = Mesh->IFace[i].ElemGroupR;
        ibuf[2] = Mesh->IFace[i].ElemL;
        ibuf[3] = Mesh->IFace[i].ElemR;
        ibuf[4] = Mesh->IFace[i].FaceL;
        ibuf[5] = Mesh->IFace[i].FaceR;
        ibuf[6] = Mesh->IFace[i].HangNumber;
        ibuf[7] = Mesh->IFace[i].OrientL;
        ibuf[8] = Mesh->IFace[i].OrientR;
      }
      else {
        xf_InitIFace(&(Mesh->IFace[i]));
      }
      
      ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 9*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
      
      Mesh->IFace[i].ElemGroupL = ibuf[0];
      Mesh->IFace[i].ElemGroupR = ibuf[1];
      Mesh->IFace[i].ElemL      = ibuf[2];
      Mesh->IFace[i].ElemR      = ibuf[3];
      Mesh->IFace[i].FaceL      = ibuf[4];
      Mesh->IFace[i].FaceR      = ibuf[5];
      Mesh->IFace[i].HangNumber = ibuf[6];
      Mesh->IFace[i].OrientL    = ibuf[7];
      Mesh->IFace[i].OrientR    = ibuf[8];
    }//nIFace
    
    //boundary faces
    if (myRank != 0) {
      ierr = xf_Error(xf_Alloc((void **)&Mesh->BFaceGroup, 
                               Mesh->nBFaceGroup, sizeof(xf_BFaceGroup)));
      if (ierr != xf_OK) return ierr;
    }
    
    for (i = 0; i < Mesh->nBFaceGroup; i++){
      //get original title
      if (myRank == 0){
        size = strlen(Mesh->BFaceGroup[i].Title)+1;
        strncpy(string, Mesh->BFaceGroup[i].Title, size);
      }
      
      ierr = xf_Error(xf_MPI_Bcast((void *)&size, 
                                   sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
      
      //broadcast title
      ierr = xf_Error(xf_MPI_Bcast((void *)string, 
                                   size*sizeof(char), 0));
      if (ierr != xf_OK) return ierr;
      
      //put title in mesh structure
      if (myRank != 0){
        ierr = xf_Error(xf_AllocString(&Mesh->BFaceGroup[i].Title, xf_MAXSTRLEN, 
                                       string));
        if (ierr != xf_OK) return ierr;
      }
      
      //nBFace
      ierr = xf_Error(xf_MPI_Bcast((void *)&Mesh->BFaceGroup[i].nBFace, 
                                   sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
      
      //BFaces
      if (myRank != 0){
        ierr = xf_Error(xf_Alloc((void **)&Mesh->BFaceGroup[i].BFace, 
                                 Mesh->BFaceGroup[i].nBFace, sizeof(xf_BFace)));
        if (ierr != xf_OK) return ierr;
      }
      for (j = 0; j < Mesh->BFaceGroup[i].nBFace; j++){
        if (myRank != 0)
          xf_InitBFace(&(Mesh->BFaceGroup[i].BFace[j]));
        
        if (myRank == 0){
          ibuf[0] = Mesh->BFaceGroup[i].BFace[j].ElemGroup;
          ibuf[1] = Mesh->BFaceGroup[i].BFace[j].Elem;
          ibuf[2] = Mesh->BFaceGroup[i].BFace[j].Face;
          ibuf[3] = Mesh->BFaceGroup[i].BFace[j].Orient;
        }
        ierr = xf_Error(xf_MPI_Bcast((void *) ibuf, 4*sizeof(int), 0));
        if (ierr != xf_OK) return ierr;
        
        Mesh->BFaceGroup[i].BFace[j].ElemGroup = ibuf[0];
        Mesh->BFaceGroup[i].BFace[j].Elem      = ibuf[1];
        Mesh->BFaceGroup[i].BFace[j].Face      = ibuf[2];
        Mesh->BFaceGroup[i].BFace[j].Orient    = ibuf[3];
      }//nBFace
    }//nBFaceGroup
    
    //element groups
    if (myRank != 0) {
      ierr = xf_Error(xf_Alloc((void **)&Mesh->ElemGroup, Mesh->nElemGroup, 
                               sizeof(xf_ElemGroup)));
      if (ierr != xf_OK) return ierr;
    }
    
    for (i = 0; i < Mesh->nElemGroup; i++){
      if (myRank != 0)
        xf_InitElemGroup(&(Mesh->ElemGroup[i]));
      
      if (myRank == 0){
        ibuf[0] = Mesh->ElemGroup[i].CutFlag;
        ibuf[1] = Mesh->ElemGroup[i].QBasis;
        ibuf[2] = Mesh->ElemGroup[i].QOrder;
        ibuf[3] = Mesh->ElemGroup[i].nElem;
        ibuf[4] = Mesh->ElemGroup[i].nNode;
      }
      
      ierr = xf_Error(xf_MPI_Bcast((void *)ibuf, 5*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
      
      Mesh->ElemGroup[i].CutFlag = ibuf[0];
      Mesh->ElemGroup[i].QBasis  = ibuf[1];
      Mesh->ElemGroup[i].QOrder  = ibuf[2];
      Mesh->ElemGroup[i].nElem   = ibuf[3];
      Mesh->ElemGroup[i].nNode   = ibuf[4];
      
      if (myRank != 0){
        ierr = xf_Error(xf_Alloc((void **)&Mesh->ElemGroup[i].nFace, 
                                 Mesh->ElemGroup[i].nElem, sizeof(int)));
        if (ierr != xf_OK) return ierr;
      }
      
      ierr = xf_Error(xf_MPI_Bcast((void *)Mesh->ElemGroup[i].nFace, 
                                   Mesh->ElemGroup[i].nElem*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
      
      if (myRank != 0){
        ierr = xf_Error(xf_VAlloc2((void ***)&Mesh->ElemGroup[i].Face, 
                                   Mesh->ElemGroup[i].nElem, 
                                   Mesh->ElemGroup[i].nFace, sizeof(xf_Face)));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_Alloc2((void ***)&Mesh->ElemGroup[i].Node, 
                                  Mesh->ElemGroup[i].nElem, 
                                  Mesh->ElemGroup[i].nNode, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
      }
      
      for (j = 0; j < Mesh->ElemGroup[i].nElem; j++)
        for (k = 0; k < Mesh->ElemGroup[i].nFace[j]; k++){
          if (myRank == 0)
            FaceBuf = Mesh->ElemGroup[i].Face[j][k];
          
          ierr = xf_Error(xf_MPI_Bcast((void *)&FaceBuf, sizeof(xf_Face), 0));
          if (ierr != xf_OK) return ierr;
          
          Mesh->ElemGroup[i].Face[j][k] = FaceBuf;
        }
      
      ierr = xf_Error(xf_MPI_Bcast((void *)Mesh->ElemGroup[i].Node[0], 
                                   Mesh->ElemGroup[i].nElem*
                                   Mesh->ElemGroup[i].nNode*sizeof(int), 0));
      if (ierr != xf_OK) return ierr;
    }//nElemGroup
    xf_MPI_Barrier();
    xf_printf("done.\n");
    
    
    
  }//if nProc > 1
  
  //debugging
  /*Mesh = (*pMesh);
    your_turn = 0;
    while (your_turn < nProc){
      if (myRank == your_turn){
        sprintf(filename,"BFaceGroups_mr%d-%d.txt",myRank,nProc);
        fid = fopen(filename, "w");
        for (i = 0; i < Mesh->nBFaceGroup; i++) {
          fprintf(fid,"Mesh->BFaceGroup[%d].Title = %s\n",i,Mesh->BFaceGroup[i].Title);
          for (j = 0; j < Mesh->BFaceGroup[i].nBFace; j++) {
            fprintf(fid,"Mesh->BFaceGroup[%d].BFace[%d].ElemGroup = %d\n",i,j,Mesh->BFaceGroup[i].BFace[j].ElemGroup);
            fprintf(fid,"Mesh->BFaceGroup[%d].BFace[%d].Elem = %d\n",i,j,Mesh->BFaceGroup[i].BFace[j].Elem);
            fprintf(fid,"Mesh->BFaceGroup[%d].BFace[%d].Face = %d\n",i,j,Mesh->BFaceGroup[i].BFace[j].Face);
            fprintf(fid,"Mesh->BFaceGroup[%d].BFace[%d].Orient = %d\n",i,j,Mesh->BFaceGroup[i].BFace[j].Orient);
          }
        }
        fclose(fid);
        your_turn++;
      }
      ierr = xf_Error(xf_MPI_Allreduce(&your_turn, 1, xfe_SizeInt, xfe_MPI_MAX));
      if (ierr != xf_OK) return ierr;
    }
    
    your_turn = 0;
    while (your_turn < nProc){
      if (myRank == your_turn){
        sprintf(filename,"IFaces_mr%d-%d.txt",myRank,nProc);
        fid = fopen(filename, "w");
        for (i = 0; i < Mesh->nIFace; i++) {
          fprintf(fid,"Mesh->IFace[%d].ElemGroupL = %d\n",i,Mesh->IFace[i].ElemGroupL);
          fprintf(fid,"Mesh->IFace[%d].ElemGroupR = %d\n",i,Mesh->IFace[i].ElemGroupR);
          fprintf(fid,"Mesh->IFace[%d].ElemL = %d\n",i,Mesh->IFace[i].ElemL);
          fprintf(fid,"Mesh->IFace[%d].ElemR = %d\n",i,Mesh->IFace[i].ElemR);
          fprintf(fid,"Mesh->IFace[%d].FaceL = %d\n",i,Mesh->IFace[i].FaceL);
          fprintf(fid,"Mesh->IFace[%d].FaceR = %d\n",i,Mesh->IFace[i].FaceR);
          fprintf(fid,"Mesh->IFace[%d].HangNumber = %d\n",i,Mesh->IFace[i].HangNumber);
          fprintf(fid,"Mesh->IFace[%d].OrientL = %d\n",i,Mesh->IFace[i].OrientL);
          fprintf(fid,"Mesh->IFace[%d].OrientR = %d\n",i,Mesh->IFace[i].OrientR);
        }
        fclose(fid);
        your_turn++;
      }
      ierr = xf_Error(xf_MPI_Allreduce(&your_turn, 1, xfe_SizeInt, xfe_MPI_MAX));
      if (ierr != xf_OK) return ierr;
    }
    
    your_turn = 0;
    while (your_turn < nProc){
      if (myRank == your_turn){
        sprintf(filename,"ElemGroups_mr%d-%d.txt",myRank,nProc);
        fid = fopen(filename, "w");
        for (i = 0; i < Mesh->nElemGroup; i++) {
          fprintf(fid,"Mesh->ElemGroup[%d].CutFlag = %d\n",i,Mesh->ElemGroup[i].CutFlag);
          fprintf(fid,"Mesh->ElemGroup[%d].QBasis = %d\n",i,Mesh->ElemGroup[i].QBasis);
          fprintf(fid,"Mesh->ElemGroup[%d].QOrder = %d\n",i,Mesh->ElemGroup[i].QOrder);
          fprintf(fid,"Mesh->ElemGroup[%d].nElem = %d\n",i,Mesh->ElemGroup[i].nElem);
          fprintf(fid,"Mesh->ElemGroup[%d].nNode = %d\n",i,Mesh->ElemGroup[i].nNode);
          for (j = 0; j < Mesh->ElemGroup[i].nElem; j++) {
            for (k = 0; k < Mesh->ElemGroup[i].nFace[j]; k++) {
              fprintf(fid,"Mesh->ElemGroup[%d].Face[%d][%d].Group = %d\n",i,j,k,Mesh->ElemGroup[i].Face[j][k].Group);
              fprintf(fid,"Mesh->ElemGroup[%d].Face[%d][%d].Number = %d\n",i,j,k,Mesh->ElemGroup[i].Face[j][k].Number);
            }
            for (k = 0; k < Mesh->ElemGroup[i].nNode; k++) {
              fprintf(fid,"Mesh->ElemGroup[%d].Node[%d][%d] = %d\n",i,j,k,Mesh->ElemGroup[i].Node[j][k]);
            }
          }
        }
        fclose(fid);
        your_turn++;
      }
      ierr = xf_Error(xf_MPI_Allreduce(&your_turn, 1, xfe_SizeInt, xfe_MPI_MAX));
      if (ierr != xf_OK) return ierr;
    }*/
     
  
  return xf_OK;
}

