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
 FILE:  xf_Line.c
 
 This file contains functions dealing with creation and ordering of
 lines of elements for line-implicit preconditioning.
 
 */

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_LinearSolverStruct.h"
#include "xf_Data.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_MeshTools.h"
#include "xf_Math.h"
#include "xf_MPI.h"
#include "xf_Basis.h"


/******************************************************************/
//   FUNCTION Definition: xf_AllocLine
static int 
xf_AllocLine( xf_Line *Line, int nelem ){
  
  int ierr;
  
  Line->nelem = nelem;
  
  ierr = xf_Error(xf_Alloc((void **) &Line->egrp, nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &Line->elem, nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &Line->face, nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReAllocLine
static int 
xf_ReAllocLine( xf_Line *Line, int nelem ){
  
  int ierr;
  
  Line->nelem = nelem;
  
  ierr = xf_Error(xf_ReAlloc((void **) &Line->egrp, nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc((void **) &Line->elem, nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc((void **) &Line->face, nelem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyLine
static int 
xf_DestroyLine( xf_Line *Line, enum xfe_Bool DestroySelf ){
  
  int ierr;
  
  if (Line == NULL) return xf_OK;
  
  xf_Release( (void *) Line->egrp);
  xf_Release( (void *) Line->elem);
  xf_Release( (void *) Line->face);
  
  if (DestroySelf) xf_Release( (void *) Line);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReverseLine
static int 
xf_ReverseLine(xf_Mesh *Mesh, xf_Line *Line){
  // this function reverses the direction of elements in a line
  int ierr, n, i, k, face;
  
  n = Line->nelem;
  for (i=0; i<n/2; i++){
    swap(Line->egrp[i], Line->egrp[n-i-1], k);
    swap(Line->elem[i], Line->elem[n-i-1], k);
  }
  
  // now reconnect the faces
  for (i=0; i<n-1; i++){
    ierr = xf_Error(xf_CommonFace(Mesh, Line->egrp[i], Line->elem[i],
                                  Line->egrp[i+1], Line->elem[i+1],
                                  Line->face+i, &face));
    if (ierr == xf_NOT_FOUND) return xf_Error(xf_LINE_ERROR);
    if (ierr != xf_OK) return ierr;
  } // i
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_CreateLineSet
static int 
xf_CreateLineSet( xf_LineSet **pLineSet ){
  
  int ierr;
  
  ierr = xf_Error(xf_Alloc((void **) pLineSet, 1, sizeof(xf_LineSet)));
  if (ierr != xf_OK) return ierr;
  
  (*pLineSet)->nLine     = 0;
  (*pLineSet)->Line      = NULL;
  (*pLineSet)->Elem2Line = NULL;
  (*pLineSet)->Face2M    = NULL;
  (*pLineSet)->negrp     = 0;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyLineSet
int 
xf_DestroyLineSet( xf_LineSet *LineSet ){
  
  int ierr, i, iLine;
  
  if (LineSet == NULL) return xf_OK;
  
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    ierr = xf_Error(xf_DestroyLine(LineSet->Line + iLine, xfe_False));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *) LineSet->Line);
  
  xf_Release2((void **) LineSet->Elem2Line);
  
  if (LineSet->Face2M != NULL){
    for (i=0; i<LineSet->negrp; i++) xf_Release2( (void **) LineSet->Face2M[i]);
    xf_Release( (void *) LineSet->Face2M);
  }
  
  xf_Release((void *) LineSet);
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AddElem2List
static int
xf_AddElem2List(int egrp, int elem, int face, int **pelist, 
                int *nelem, int *nelem0)
{
  int ierr;
  
  (*nelem)++;
  if ((*nelem) > (*nelem0)){
    (*nelem0) = max((*nelem), 3*(*nelem0));
    ierr = xf_Error(xf_ReAlloc((void **) pelist, 3*(*nelem0), sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  (*pelist)[3*(*nelem)-3] = egrp;
  (*pelist)[3*(*nelem)-2] = elem;
  (*pelist)[3*(*nelem)-1] = face;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_PreconditionerLineCheck
void
xf_PreconditionerLineCheck(enum xfe_PreconditionerType Preconditioner,
                           enum xfe_Bool *CRequired, enum xfe_Bool *SortLines)
{
  (*CRequired) = ((Preconditioner == xfe_PreconditionerLineJacobi) || 
                  (Preconditioner == xfe_PreconditionerLineGS    ));
  (*SortLines) =  (Preconditioner == xfe_PreconditionerLineGS);
}

/******************************************************************/
//   FUNCTION Definition: xf_FindLineConnectivity
int
xf_FindLineConnectivity(xf_All *All, xf_Vector **pC)
{
  int ierr, negrp, i, j, *rvec;
  int negrphalo, nface;
  
  // rvec = vector of ranks = nFace for each group
  negrp = All->Mesh->nElemGroup;
  negrphalo = ((All->Mesh->ParallelInfo == NULL) ? negrp : 2*negrp);
  ierr = xf_Error(xf_Alloc((void **) &rvec, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<negrp; i++) rvec[i] = 0;
  for (i=0; i<negrphalo; i++)
    for (j=0; j<All->Mesh->ElemGroup[i].nElem; j++)
      rvec[i%negrp] = max(rvec[i%negrp],All->Mesh->ElemGroup[i].nFace[j]);
  
  if (All->Mesh->ParallelInfo != NULL) {
    ierr = xf_Error(xf_MPI_Allreduce(rvec, negrp, xfe_SizeInt, xfe_MPI_MAX));
    if (ierr != xf_OK) return ierr;
  }
  
  ierr = xf_Error(xf_FindVector(All, "LineConnectivity", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                NULL, NULL, NULL, NULL, rvec, xfe_SizeReal, xfe_False, xfe_True,
                                NULL, pC, NULL));
  if (ierr != xf_OK) return ierr;
  
  // release rvec
  xf_Release( (void *) rvec);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_StoreLines4Viz
static int
xf_StoreLines4Viz(xf_All *All, xf_LineSet *LineSet)
{
  // Stores LineSet in a real vector "LineID"
  
  int ierr, iLine, egrp, elem, ie;
  xf_Vector *LineID;
  xf_Data *D;
  
  ierr = xf_Error(xf_FindVector(All, "LineID", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                NULL, NULL, NULL, NULL, NULL, xfe_SizeReal, xfe_False,  
                                xfe_True, &D, &LineID, NULL));
  if (ierr != xf_OK) return ierr;
  
  D->ReadWrite = xfe_True; // make data writeable
  
  // set LineID = line # for each elem
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    for (ie=0; ie<LineSet->Line[iLine].nelem; ie++){
      egrp = LineSet->Line[iLine].egrp[ie];
      elem = LineSet->Line[iLine].elem[ie];
      LineID->GenArray[egrp].rValue[elem][0] = (real) (iLine);
    }
  } // iLine
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateFace2M
static int
xf_CreateFace2M(xf_All *All, xf_LineSet *LineSet, xf_JacobianMatrix *R_U)
{
  /* Creates LineSet->Face2M, which stores information about the
   neighbors of each element: whether they are adjacent in the same
   line, or whether they are not but otherwise valid elements, or
   whether they are halo or boundary elements.  This structure is
   used in the solver for speedup in traversing through the
   lines. */
  
  int ierr, k, egrp, elem, face, negrp;
  int egrpN, elemN, faceN, fi, pface;
  int iLine, ie, *R_UegrpN, *R_UelemN;
  xf_Mesh *Mesh;
  
  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  if (LineSet->Face2M != NULL) return xf_Error(xf_INPUT_ERROR);
  
  // For speedup in traversing lines
  
  LineSet->negrp = negrp;
  ierr = xf_Error(xf_Alloc( (void **) &LineSet->Face2M, negrp, sizeof(int **)));
  if (ierr != xf_OK) return ierr;
  for (egrp=0; egrp<negrp; egrp++){
    ierr = xf_Error(xf_VAlloc2( (void ***) &LineSet->Face2M[egrp], 
                               Mesh->ElemGroup[egrp].nElem,
                               Mesh->ElemGroup[egrp].nFace, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  // loop over lines
  for (iLine=0; iLine<LineSet->nLine; iLine++){      
    pface = -1; // initialize previous face number
    // loop over elements in a line
    for (ie=0; ie<LineSet->Line[iLine].nelem; ie++){
      
      // pull off info of elem ie on line
      egrp = LineSet->Line[iLine].egrp[ie];
      elem = LineSet->Line[iLine].elem[ie];
      face = LineSet->Line[iLine].face[ie];
      
      // loop over faces of elem
      for (fi=0; fi<Mesh->ElemGroup[egrp].nFace[elem]; fi++){
        
        // default is neither next/previous nor boundary/halo
        LineSet->Face2M[egrp][elem][fi] = 0;
        
        // previous or next elements in line are marked with a 1
        if ((fi==face) || (fi==pface)) LineSet->Face2M[egrp][elem][fi] = 1;
        
        egrpN = R_U->egrpN[egrp][elem][fi];
        elemN = R_U->elemN[egrp][elem][fi];
        
        // halo or boundary elements are marked with a -1
        if ((egrpN >= negrp) || (egrpN < 0)) LineSet->Face2M[egrp][elem][fi] = -1;
      } // fi
      
      // set previous face number
      if (face >= 0) pface = R_U->faceN[egrp][elem][face]; 
      
    } // ie
  } // iLine 
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_LineCreation
static int
xf_LineCreation(xf_Mesh *Mesh, xf_Vector *C, xf_LineSet *LineSet)
{
  int ierr, negrp, *nElem, nelemtot, elem, k, j, egrp, face;
  int isafety, nface, SeedEgrp, SeedElem, SeedFace[2], step, tot;
  int nelem[2], egrpN, elemN, faceN, nfaceN, count, ielem, nLine0;
  enum xfe_Bool done, TerminatePath;
  real *CS, *CN, cval;
  int **Elem2Line, nelem0[2], *elist[2];
  xf_Line Line;
  
  negrp = Mesh->nElemGroup;
  nLine0 = 0;
  
  // Create Elem2Line
  ierr = xf_Error(xf_GetnElem(Mesh, &nElem, &nelemtot));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VAlloc2((void ***)&Elem2Line, negrp, nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // set Elem2Line values to -1
  for (elem=0; elem<nelemtot; elem++) Elem2Line[0][elem] = -1;
  
  /*** START OF STAGE I: line creation  ***/
  
  // find first seed element (just first element in a group with more than zero elements)
  SeedEgrp = 0;
  SeedElem = 0;
  done = xfe_True;
  while ((SeedEgrp < Mesh->nElemGroup) && (done)){
    while ((SeedElem < Mesh->ElemGroup[SeedEgrp].nElem) && (done)){
      if (Elem2Line[SeedEgrp][SeedElem] < 0) // not yet part of a line
        done = xfe_False;
      else
        SeedElem++;
    }
    if (done) { // moving on to next elem group
      SeedElem = 0;
      SeedEgrp++;
    }
  }
  
  // initialize path element lists
  nelem0[0] = nelem0[1] = 0;
  elist[0] = elist[1] = NULL;
  
  // begin loop over seed elements
  LineSet->nLine = 0;
  done = xfe_False;
  isafety = 0;
  while ((!done) && (isafety < nelemtot+1)){
    
    if (Elem2Line[SeedEgrp][SeedElem] >= 0) return xf_Error(xf_LINE_ERROR);
    
    // obtain connectivity for seed
    CS = C->GenArray[SeedEgrp].rValue[SeedElem];
    
    // obtain number of faces for the seed elem
    if ((nface = Mesh->ElemGroup[SeedEgrp].nFace[SeedElem]) < 2) 
      return xf_Error(xf_LINE_ERROR);
    
    // find two faces with largest C values, CS[SeedFace[0]] > CS[SeedFace[1]]
    for (k=1, SeedFace[0]=0; k<nface; k++)
      if (CS[k] > CS[SeedFace[0]]) SeedFace[0] = k;
    
    for (k=1, SeedFace[1]=(SeedFace[0]+1)%nface; k<(nface-1); k++)
      if (CS[face = (SeedFace[0]+1+k)%nface] > CS[SeedFace[1]]) SeedFace[1] = face;
    
    
    if (SeedFace[0] == SeedFace[1]) return xf_Error(xf_LINE_ERROR);
    
    // Perform forward and backward path searches from seed
    for (step=0; step<2; step++){ // 0 = forward, 1 = backward
      
      nelem[step] = 0; // initialize number of elements in this path
      
      if (step == 0){ // add seed as first elem in forward path
        ierr = xf_Error(xf_AddElem2List(SeedEgrp, SeedElem, SeedFace[step], 
                                        elist+step, nelem+step, nelem0+step));
        if (ierr != xf_OK) return ierr;
        Elem2Line[SeedEgrp][SeedElem] = LineSet->nLine;
      }
      
      // initialize start of path to seed
      egrp = SeedEgrp;
      elem = SeedElem;
      face = SeedFace[step];
      cval = CS[face];
      
      // start path
      TerminatePath = xfe_False;
      while ((!TerminatePath) && (isafety < nelemtot+1)){
        
        ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, 
                                              &egrpN, &elemN, &faceN));
        if (ierr != xf_OK) return ierr;
        
        // stop if hit domain boundary (or halo)
        if (TerminatePath = ((egrpN < 0) || (egrpN >= negrp))) break;
        
        // stop if neighbor is part of a line
        if (TerminatePath = (Elem2Line[egrpN][elemN] >= 0)) break;
        
        // stop if neighbor has 2 other faces of higher C
        CN = C->GenArray[egrpN].rValue[elemN];
        if (fabs(cval - CN[faceN]) > 1e-14) return xf_Error(xf_INPUT_ERROR);
        nfaceN = Mesh->ElemGroup[egrpN].nFace[elemN];
        for (k=1, count=0; k<nfaceN; k++)
          if (CN[(faceN+k)%nfaceN] > cval) count++;
        if (TerminatePath = (count > 1)) break;
        
        // set elem <- elemN
        egrp = egrpN;
        elem = elemN;
        
        // set face to highest connectivity, not counting faceN
        for (k=1, face=(faceN+1)%nfaceN; k<(nfaceN-1); k++)
          if (CN[j = (faceN+1+k)%nfaceN] > CN[face]) face = j;
        cval = CN[face];
        
        
        // add elemN to list
        if (step == 0) // for forward path, store next face
          ierr = xf_Error(xf_AddElem2List(egrpN, elemN, face, elist+step, 
                                          nelem+step, nelem0+step));
        else // for backward path, store lookback face of elemN
          ierr = xf_Error(xf_AddElem2List(egrpN, elemN, faceN, elist+step, 
                                          nelem+step, nelem0+step));
        if (ierr != xf_OK) return ierr;
        Elem2Line[egrpN][elemN] = LineSet->nLine;
        
        isafety++;
      }
      if (!TerminatePath) return xf_Error(xf_LINE_ERROR);
      
    } // end for step
    
    
    // allocate a Line
    ierr = xf_Error(xf_AllocLine(&Line, nelem[0] + nelem[1]));
    if (ierr != xf_OK) return ierr;
    
    // add path elements to line
    ielem=0;
    for (k=nelem[1]-1; k>=0; k--){ // backward path first
      Line.egrp[ielem] = egrp = elist[1][3*k+0];
      Line.elem[ielem] = elem = elist[1][3*k+1];
      Line.face[ielem] = face = elist[1][3*k+2];
      Elem2Line[egrp][elem] = LineSet->nLine;
      ielem++;
    }
    for (k=0; k<nelem[0]; k++){ // forward path second
      Line.egrp[ielem] = egrp = elist[0][3*k+0];
      Line.elem[ielem] = elem = elist[0][3*k+1];
      Line.face[ielem] = face = elist[0][3*k+2];
      Elem2Line[egrp][elem] = LineSet->nLine;
      ielem++;
    }
    if (ielem != (nelem[0] + nelem[1])) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    // set .face of last elem to -1, to signify end of line
    Line.face[ielem-1] = -1;
    
    tot += ielem;
    
    // increment line counter
    LineSet->nLine++;
    
    
    // reallocate LineSet->Line if necessary
    if (LineSet->nLine > nLine0){
      nLine0 = max(LineSet->nLine, 2*nLine0);
      ierr = xf_Error(xf_ReAlloc( (void **) &LineSet->Line, nLine0, sizeof(xf_Line)));
      if (ierr != xf_OK) return ierr;
    }
    
    // add Line to LineSet
    LineSet->Line[LineSet->nLine-1] = Line;
    
    // find next seed element
    done = xfe_True;
    while ((SeedEgrp < Mesh->nElemGroup) && (done)){
      while ((SeedElem < Mesh->ElemGroup[SeedEgrp].nElem) && (done)){
        if (Elem2Line[SeedEgrp][SeedElem] < 0) // not yet part of a line
          done = xfe_False;
        else
          SeedElem++;
      }
      if (done) { // moving on to next elem group
        SeedElem = 0;
        SeedEgrp++;
      }
    }
    
    isafety++;
  } // end while loop over 
  
  if (!done) return xf_Error(xf_LINE_ERROR); // seed elems remaining after loop
  
  /*** END OF STAGE I ***/
  
  // store Elem2Line in LineSet
  LineSet->Elem2Line = Elem2Line;
  
  // Release memory
  xf_Release( (void  *) nElem);
  
  xf_Release( (void  *) elist[0]);
  xf_Release( (void  *) elist[1]);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateLines
static int
xf_LineConnection(xf_Mesh *Mesh, xf_Vector *C, xf_LineSet *LineSet)
{
  int ierr, egrp, elem, nelemtot, negrp, iLine, ne, isafety, step;
  int **Elem2End, *nElem, ie, nface, k, face, iLineN, nelem[2], stepN;
  int egrpN, elemN, egrpP, elemP, nfaceN, faceN, iLineP;
  real *CS, *CN;
  xf_Line *pLine;
  
  negrp = Mesh->nElemGroup;
  
  /*** START OF STAGE II: line connection ***/
  
  ierr = xf_Error(xf_GetnElem(Mesh, &nElem, &nelemtot));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_VAlloc2((void ***) &Elem2End, 
                             negrp, nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // set Elem2End values to -1
  for (elem=0; elem<nelemtot; elem++) Elem2End[0][elem] = -1;
  
  /* Loop over Lines, set Elem2End to line# for start and end points */
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    pLine = LineSet->Line + iLine;
    ne    = pLine->nelem;
    Elem2End[pLine->egrp[0   ]][pLine->elem[0   ]] = iLine;
    Elem2End[pLine->egrp[ne-1]][pLine->elem[ne-1]] = iLine;
  } // iLine
  
  // loop over lines
  isafety = 0;
  for (iLine=0; iLine<LineSet->nLine; iLine++){
    // Continue if line no longer exists
    pLine = LineSet->Line + iLine;
    if ((ne=pLine->nelem) == 0) continue;
    // loop over endpoints of current line
    for (step=0; step<2; step++){
      // elem at endpoint
      ie = ((step==0) ? 0 : ne-1);
      egrp = pLine->egrp[ie];
      elem = pLine->elem[ie];
      
      // obtain connectivity for this elem
      CS = C->GenArray[egrp].rValue[elem];
      // obtain number of faces for this elem
      if ((nface=Mesh->ElemGroup[egrp].nFace[elem]) < 2) 
        return xf_Error(xf_LINE_ERROR);
      
      /* find face with highest connectivity such that:
       1) adjacent elem is an endpoint of another line
       2) face is a boundary face
       */
      for (k=0, face=-1; k<nface; k++){
        ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, k, 
                                              &egrpN, &elemN, NULL));
        if (ierr != xf_OK) return ierr;
        
        if ((egrpN < 0) || (egrpN >= negrp) ||
            (((iLineN = Elem2End[egrpN][elemN]) >= 0)&& (iLineN != iLine)))
          face = ( ((face==-1) || (CS[k] > CS[face])) ? k : face);
      } // k
      // if no such face exists or if a bface, cannot join, so continue
      if (face == -1) continue;
      ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, face, 
                                            &egrpN, &elemN, NULL));
      if (ierr != xf_OK) return ierr;
      if (egrpN < 0) continue;
      
      // now we have an endpoint of another line: call it (egrpN, elemN)
      if (egrpN >= negrp) continue; // cannot connect across halos
      iLineN = Elem2End[egrpN][elemN];
      
      // does egrpN, elemN want to be joined to elem?
      CN = C->GenArray[egrpN].rValue[elemN]; // connectivity for neighbor
      nfaceN = Mesh->ElemGroup[egrpN].nFace[elemN];
      for (k=0, faceN=-1; k<nfaceN; k++){ 
        ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrpN, elemN, k, 
                                              &egrpP, &elemP, NULL));
        if (ierr != xf_OK) return ierr;
        
        if ((egrpP < 0) || (egrpP >= negrp) ||
            (((iLineP = Elem2End[egrpP][elemP]) >= 0) && (iLineP != iLineN)))
          faceN = ( ((faceN==-1) || (CN[k] > CN[faceN])) ? k : faceN);
      } // k
      
      // if not a single candidate face exists on N, error
      if (faceN == -1) return xf_Error(xf_LINE_ERROR);
      
      ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrpN, elemN, faceN, 
                                            &egrpP, &elemP, NULL));
      if (ierr != xf_OK) return ierr;
      // if N wants another elem (or boundary), continue
      if ((egrpP != egrp) || (elemP != elem)) continue;
      
      /* At this point, lines iLine and iLineN need to be connected */
      
      // number of elements in the lines to be joined
      nelem[0] = LineSet->Line[iLine ].nelem;
      nelem[1] = LineSet->Line[iLineN].nelem;
      
      // Do we need to traverse iLine or iLineN in reverse?
      stepN = ((LineSet->Line[iLineN].egrp[0] != egrpN) ||
               (LineSet->Line[iLineN].elem[0] != elemN));
      
      if ((stepN) && ((LineSet->Line[iLineN].egrp[nelem[1]-1] != egrpN) ||
                      (LineSet->Line[iLineN].elem[nelem[1]-1] != elemN))){
        return xf_Error(xf_CODE_LOGIC_ERROR);
      }
      /*
       iLine     iLineN
       <-------||------->   step = 0, stepN = 0
       <-------||<-------   step = 0, stepN = 1
       ------->||------->   step = 1, stepN = 0
       ------->||<-------   step = 1, stepN = 1
       */
      if (step  == 0){ // reverse elems in iLine
        ierr = xf_Error(xf_ReverseLine(Mesh, LineSet->Line+iLine ));
        if (ierr != xf_OK) return ierr;
      }
      if (stepN == 1){ // reverse elems in iLineN
        ierr = xf_Error(xf_ReverseLine(Mesh, LineSet->Line+iLineN));
        if (ierr != xf_OK) return ierr;
      }
      
      // reallocate .egrp, .elem, .face in LineSet->Line+iLine
      ierr = xf_Error(xf_ReAllocLine(LineSet->Line+iLine, nelem[0] + nelem[1]));
      if (ierr != xf_OK) return ierr;
      
      // copy over egrp,elem,face from iLineN to iLine
      for (ie=0; ie<nelem[1]; ie++){
        LineSet->Line[iLine].egrp[nelem[0]+ie] = LineSet->Line[iLineN].egrp[ie];
        LineSet->Line[iLine].elem[nelem[0]+ie] = LineSet->Line[iLineN].elem[ie];
        LineSet->Line[iLine].face[nelem[0]+ie] = LineSet->Line[iLineN].face[ie];
      }
      // set connecting face between the two lines
      LineSet->Line[iLine].face[nelem[0]-1] = face;
      
      // set correct number of elements
      LineSet->Line[iLine ].nelem = nelem[0] + nelem[1];
      LineSet->Line[iLineN].nelem = 0;
      
      // release memory from iLineN
      ierr = xf_Error(xf_DestroyLine(LineSet->Line+iLineN, xfe_False));
      if (ierr != xf_OK) return ierr;
      
      // (egrp,elem) and (egrpN,elemN) are no longer endpoints
      Elem2End[egrp ][elem ] = -1;
      Elem2End[egrpN][elemN] = -1;
      
      // set new endpoints for current long line
      pLine = LineSet->Line + iLine;
      ne    = nelem[0] + nelem[1];
      Elem2End[pLine->egrp[0   ]][pLine->elem[0   ]] = iLine;
      Elem2End[pLine->egrp[ne-1]][pLine->elem[ne-1]] = iLine;
      
      // mark current line to be rechecked, and check for infinite loop
      step = 2; 
      iLine--;
      if (isafety++ > 2*LineSet->nLine) return xf_Error(xf_LINE_ERROR);
    } // step
  } // iLine
  
  // Remove stale lines (those with nelem=0)
  for (iLine=0, k=0; k<LineSet->nLine; k++){
    if (LineSet->Line[k].nelem == 0) continue;
    if (iLine != k) LineSet->Line[iLine] = LineSet->Line[k];
    iLine++;
  }
  
  if (iLine <= 0) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  // Re-allocate number of lines if merging occurred, 
  if (iLine != LineSet->nLine){
    ierr = xf_Error(xf_ReAlloc( (void **) &LineSet->Line, iLine, sizeof(xf_Line)));
    if (ierr != xf_OK) return ierr;
    
    // Stage II reduced # lines from LineSet->nLine to iLine
    LineSet->nLine = iLine;
    
    // change Elem2Line using new line numbering
    for (iLine=0, k=0; k<LineSet->nLine; k++){
      pLine = LineSet->Line+k;
      for (ie=0; ie<pLine->nelem; ie++){
        LineSet->Elem2Line[pLine->egrp[ie]][pLine->elem[ie]] = iLine;
      } // ie
    } // iLine
  }
  
  // release Elem2End
  xf_Release2((void **) Elem2End);
  xf_Release( (void  *) nElem);
  /***  END OF STAGE II  ***/
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_LineReNumber
/* this function is used only to facilitate parallel debugging 
 since it makes the line ordering unique */
static int
xf_LineReNumber(xf_Mesh *Mesh, xf_LineSet *LineSet)
{
  int ierr, i, *FirstElem, *LineRank, egrp, elem;
  xf_LineSet *LineSet_buffer;
  
  ierr = xf_Error(xf_Alloc((void**)&FirstElem, LineSet->nLine,
                           sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void**)&LineRank, LineSet->nLine, 
                           sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  for (i=0;i<LineSet->nLine;i++)
    FirstElem[i] = LineSet->Line[i].elem[0];
  
  ierr = xf_Error(xf_SortIntPos(FirstElem, LineSet->nLine, 
                                LineRank, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void**)&LineSet_buffer, 
                           1, sizeof(xf_LineSet)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void**)&LineSet_buffer->Line, 
                           LineSet->nLine, sizeof(xf_Line)));
  if (ierr != xf_OK) return ierr;
  
  for (i = 0; i < LineSet->nLine; i++) {
    LineSet_buffer->Line[LineRank[i]].nelem = LineSet->Line[i].nelem;
    //store pointers
    LineSet_buffer->Line[LineRank[i]].elem = LineSet->Line[i].elem;
    LineSet_buffer->Line[LineRank[i]].egrp = LineSet->Line[i].egrp;
    LineSet_buffer->Line[LineRank[i]].face = LineSet->Line[i].face;
  }
  for (i = 0; i < LineSet->nLine; i++) {
    LineSet->Line[i].nelem = LineSet_buffer->Line[i].nelem;
    LineSet->Line[i].elem = LineSet_buffer->Line[i].elem;
    LineSet->Line[i].egrp = LineSet_buffer->Line[i].egrp;
    LineSet->Line[i].face = LineSet_buffer->Line[i].face;
  }
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++) {
      LineSet->Elem2Line[egrp][elem] = LineRank[LineSet->Elem2Line[egrp][elem]];
    }
  }
  
  xf_Release((void*)LineSet_buffer->Line);
  xf_Release((void*)LineSet_buffer);
  xf_Release((void*)FirstElem);
  xf_Release((void*)LineRank);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateLines
int
xf_CreateLines(xf_All *All, xf_Vector *C, xf_JacobianMatrix *R_U, 
               xf_LineSet **pLineSet)
{
  int ierr;
  enum xfe_Bool LineViz;
  xf_LineSet *LineSet;
  xf_Mesh *Mesh;
  
  Mesh  = All->Mesh;
  
  
  ierr = xf_Error(xf_CreateLineSet(&LineSet));
  if (ierr != xf_OK) return ierr;
  //Stage I
  ierr = xf_Error(xf_LineCreation(Mesh, C, LineSet));
  if (ierr != xf_OK) return ierr;
  //Stage II
  ierr = xf_Error(xf_LineConnection(Mesh, C, LineSet));
  if (ierr != xf_OK) return ierr;
  
  if (R_U != NULL){
    // Create Face2M for speedup during solve
    ierr = xf_Error(xf_CreateFace2M(All, LineSet, R_U));
    if (ierr != xf_OK) return ierr;
    
    
    // Destroy any existing R_U->LineSet and store the one just created
    ierr = xf_Error(xf_DestroyLineSet(R_U->LineSet));
    if (ierr != xf_OK) return ierr;
    R_U->LineSet = LineSet;
  }
  
  ierr = xf_Error(xf_LineReNumber(All->Mesh, LineSet));
  if (ierr != xf_OK) return ierr;
  
  // create a LineID vizualization vector if asking for visualization
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "LineViz", &LineViz));
  if (ierr != xf_OK) return ierr;
  if (LineViz){
    ierr = xf_Error(xf_StoreLines4Viz(All, LineSet));
    if (ierr != xf_OK) return ierr;
  }
  
  if (pLineSet != NULL) {
    (*pLineSet) = LineSet;
  }
  
  return xf_OK;
  
}


/******************************************************************/
//   FUNCTION Definition: xf_TreePush
static int
xf_TreePush(int ipush, real *val, int *proot, int **Tree)
{
  // pushes index ipush with valpush = val[ipush] onto Tree
  int inode, root, i;
  enum xfe_TreeNode child;
  enum xfe_Bool done;
  real valpush, valnode;
  
  if ((*proot) < 0){
    (*proot) = ipush;
    return xf_OK;
  }
  
  valpush = val[ipush]; // value to push onto tree
  
  inode = root = (*proot);
  
  i = 0;
  
  done = xfe_False;
  while (!done){
    valnode = val[inode];
    child = ((valpush > valnode) ? xfe_ChildR : xfe_ChildL);
    if (Tree[inode][child] < 0){
      Tree[inode][child] = ipush;
      Tree[ipush][xfe_Parent] = inode;
      done = xfe_True;
    }
    else{
      inode = Tree[inode][child];
    }
  } // while !done
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_TreePop
static int
xf_TreePop(int ipop, real *val, int *root, int **Tree)
{
  // pops index ipop with valpop = val[ipop] off of Tree
  int ierr, inode;
  int parent, childL, childR;
  enum xfe_TreeNode child;
  real valpop, valnode;
  
  valpop = val[ipop]; // value to pop off of tree
  
  parent = Tree[ipop][xfe_Parent]; // parent of ipop
  childL = Tree[ipop][xfe_ChildL]; // left child
  childR = Tree[ipop][xfe_ChildR]; // right child
  
  // determine whether ipop is L or R child of parent
  if (parent >= 0){ // parent < 0 means ipop is root node
    if (Tree[parent][xfe_ChildL] == ipop) child = xfe_ChildL;
    else if (Tree[parent][xfe_ChildR] == ipop) child = xfe_ChildR;
    else return xf_Error(xf_TREE_ERROR);
  }
  
  // set -1 values for popped node
  Tree[ipop][xfe_Parent] = -1;
  Tree[ipop][xfe_ChildL] = -1;
  Tree[ipop][xfe_ChildR] = -1;
  
  // now deal with children
  if (childL >= 0){
    // left branch inherits ipop and right branch is attached to it
    Tree[childL][xfe_Parent] = parent;
    if (parent >= 0) 
      Tree[parent][child] = childL;
    else
      (*root) = childL;
    
    if (childR >= 0){ // attach right branch
      ierr = xf_Error(xf_TreePush( childR, val, &childL, Tree));
      if (ierr != xf_OK) return ierr;
    }
    
  }
  else if (childR >= 0){
    // no left child means right branch inherits ipop position
    Tree[childR][xfe_Parent] = parent;
    if (parent >= 0) 
      Tree[parent][child] = childR;
    else
      (*root) = childR;
  }
  else{ 
    // ipop had no children, set parent to -1
    if (parent >= 0) 
      Tree[parent][child] = -1;
    else if ((*root) = ipop) // ipop was sole entry in Tree
      (*root) = -1;
  }
  
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SortLines
int
xf_SortLines(xf_All *All, xf_Vector *C, xf_JacobianMatrix *R_U, 
             xf_Vector *R)
{
  int ierr, k, i, ie, fi, nN, nN0, negrp;
  int nLine, iLine, iLineN, egrp, elem, face, pface;
  int egrpN, elemN, maxLine, TreeRoot;
  int **Elem2Line, **Tree, *LineFlag, *NList, *Line2New;
  enum xfe_Bool LineViz, SortLines;
  real t, Ctot, Cval, *Rline, *RlineAdd;
  xf_Line *Line, *NewLines;
  xf_LineSet *LineSet;
  xf_Mesh *Mesh;
  
  // This code is somewhat buggy now, so not supporting it until resolve issues (parallel)
  return xf_Error(xf_NOT_SUPPORTED);
  
  // check if sorting of lines is requested
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "SortLines", &SortLines));
  if (ierr != xf_OK) return ierr;
  
  if (!SortLines) return(xf_OK);
  
  Mesh  = All->Mesh;
  negrp = Mesh->nElemGroup;
  
  if ((LineSet = R_U->LineSet) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  if ((Elem2Line = LineSet->Elem2Line) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // set nLine = # of lines
  nLine = LineSet->nLine;
  
  // allocate NewLines
  ierr = xf_Error(xf_Alloc( (void **) &NewLines, nLine, sizeof(xf_Line)));
  if (ierr != xf_OK) return ierr;
  
  
  // allocate Rline = L1 norm of R on each line
  ierr = xf_Error(xf_Alloc( (void **) &Rline, nLine, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // set Rline = 0
  for (iLine=0; iLine<nLine; iLine++) Rline[iLine] = 0.0;
  
  // loop over elems and fill in Rline
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      for (k=0, t=0.; k<R->GenArray[egrp].r; k++)
        t += fabs(R->GenArray[egrp].rValue[elem][k]);
      Rline[Elem2Line[egrp][elem]] += t;
    }
  
  // allocate space for modified Rline
  ierr = xf_Error(xf_Alloc( (void **) &RlineAdd, nLine, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // allocate memory for a static binary search tree over the lines
  ierr = xf_Error(xf_Alloc2( (void ***) &Tree, nLine, xfe_TreeLast, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (iLine=0; iLine<nLine; iLine++)
    for (k=0; k<xfe_TreeLast; k++) Tree[iLine][k] = -1; // fill tree with -1s
  
  // set root to -1
  TreeRoot = -1;
  
  // Loop over lines and push them onto Tree one at a time
  for (iLine=0; iLine<nLine; iLine++){
    ierr = xf_Error(xf_TreePush( iLine, Rline, &TreeRoot, Tree));
    if (ierr != xf_OK) return ierr;
  }
  
  // allocate Line2New: Line2New[oldline] = newline
  ierr = xf_Error(xf_Alloc( (void **) &Line2New, nLine, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // allocate LineFlag indicator over lines and set it to 0
  ierr = xf_Error(xf_Alloc( (void **) &LineFlag, nLine, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (iLine=0; iLine<nLine; iLine++) LineFlag[iLine] = 0;
  
  // initialize variables before main loop
  NList = NULL;
  nN0   = 0;
  
  // begin main reorder loop over lines
  for (iLine=0; iLine<nLine; iLine++){
    
    // find line with maximum Rline using Tree
    if ((maxLine = TreeRoot) < 0) return xf_Error(xf_TREE_ERROR);
    while ((k = Tree[maxLine][xfe_ChildR]) >= 0) maxLine = k;
    
    Line = LineSet->Line + maxLine;
    
    // pop maxLine off tree
    ierr = xf_Error(xf_TreePop( maxLine, Rline, &TreeRoot, Tree));
    if (ierr != xf_OK) return ierr;
    
    /* determine Ctot = total connectivity of line to neighbors */
    Ctot = 0.;
    nN   = 0; // number of neighbor lines
    pface = -1; // initialize previous face number
    for (ie=0; ie<Line->nelem; ie++){
      // pull off info of elem ie on line
      egrp = Line->egrp[ie]; 
      elem = Line->elem[ie]; 
      face = Line->face[ie];
      
      for (fi=0; fi<Mesh->ElemGroup[egrp].nFace[elem]; fi++){
        // do not consider previous or next elements
        if ((fi==face) || (fi==pface)) continue;
        
        egrpN = R_U->egrpN[egrp][elem][fi];
        elemN = R_U->elemN[egrp][elem][fi];
        
        Cval = C->GenArray[egrp].rValue[elem][fi];
        
        // do not consider halo or boundary elements for nN
        if ((egrpN >= negrp) || (egrpN < 0)){
          Ctot += Cval;
          continue;
        }
        
        iLineN = Elem2Line[egrpN][elemN];
        
        // only consider iLineN != maxLine and iLineN still on Tree
        if (iLineN == maxLine) continue;
        if ((Tree[iLineN][xfe_Parent] < 0) && (TreeRoot != iLineN)) continue;
        
        Ctot += Cval;
        
        if (LineFlag[iLineN] == 0) {
          LineFlag[iLineN] = 1;
          nN++;
        }
        
      } // fi
      
      // set previous face number
      if (face >= 0) pface = R_U->faceN[egrp][elem][face];
      
    } // ie
    
    // allocate NList = neighbor line list
    if (nN > nN0){
      ierr = xf_Error(xf_ReAlloc( (void **) &NList, (nN0 = nN), sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    
    
    /* add C/Ctot*MaxRline to Rline of each line in NList and fill NList */
    nN    = 0;
    pface = -1; // initialize previous face number
    for (ie=0; ie<Line->nelem; ie++){
      // pull off info of elem ie on line
      egrp = Line->egrp[ie]; 
      elem = Line->elem[ie]; 
      face = Line->face[ie];
      
      for (fi=0; fi<Mesh->ElemGroup[egrp].nFace[elem]; fi++){
        // do not consider previous or next elements
        if ((fi==face) || (fi==pface)) continue;
        
        egrpN = R_U->egrpN[egrp][elem][fi];
        elemN = R_U->elemN[egrp][elem][fi];
        
        // do not consider halo or boundary elements for NList
        if ((egrpN >= negrp) || (egrpN < 0)) continue;
        
        // line number of neighbor element
        iLineN = Elem2Line[egrpN][elemN];
        
        // only consider iLineN != maxLine and iLineN still on Tree
        if (iLineN == maxLine) continue;
        if ((Tree[iLineN][xfe_Parent] < 0) && (TreeRoot != iLineN)) continue;
        
        // store iLineN in NList if not already in there
        if (LineFlag[iLineN] == 1){
          LineFlag[iLineN] = 0 ;
          RlineAdd[iLineN] = 0.;
          NList[nN] = iLineN;
          nN++;
        }
        
        // add portion of Rline[maxLine] to neighbor
        RlineAdd[iLineN] += Rline[maxLine]*C->GenArray[egrp].rValue[elem][fi]/Ctot;
        
      } // fi
      
      // set previous face number
      if (face >= 0) pface = R_U->faceN[egrp][elem][face];
      
    } // ie
    
    // Pop each line in NList off Tree
    for (i=0; i<nN; i++){
      ierr = xf_Error(xf_TreePop( NList[i], Rline, &TreeRoot, Tree));
      if (ierr != xf_OK) return ierr;
    }
    
    // modify Rline using RlineAdd
    for (i=0; i<nN; i++) Rline[NList[i]] += RlineAdd[NList[i]];
    
    // Push each line in NList onto Tree (new Rline will be used)
    for (i=0; i<nN; i++){
      ierr = xf_Error(xf_TreePush( NList[i], Rline, &TreeRoot, Tree));
      if (ierr != xf_OK) return ierr;
    }
    
    // point iLine of NewLines to the maxLine
    NewLines[iLine] = LineSet->Line[maxLine];
    Line2New[maxLine] = iLine;
    
    // TEMPORARY
    //if ((nLine % 41) == 0) return xf_Error(-17);
    //NewLines[iLine] = LineSet->Line[(iLine*41)%nLine];
    //Line2New[(iLine*41)%nLine] = iLine;
  } // iLine (end main reorder loop)
  
  // Set LineSet->Line to NewLines
  xf_Release((void *) LineSet->Line); // release previous vector of pointers
  LineSet->Line = NewLines;
  
  // renumber Elem2Line
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      Elem2Line[egrp][elem] = Line2New[Elem2Line[egrp][elem]];
  
  // store LineID vizualization vector if asking for visualization
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "LineViz", &LineViz));
  if (ierr != xf_OK) return ierr;
  if (LineViz){
    ierr = xf_Error(xf_StoreLines4Viz(All, LineSet));
    if (ierr != xf_OK) return ierr;
  }
  
  // release memory
  xf_Release( (void  *) Rline);
  xf_Release( (void  *) RlineAdd);
  xf_Release2((void **) Tree );
  xf_Release( (void  *) NList );
  xf_Release( (void  *) LineFlag );
  xf_Release( (void  *) Line2New );
  
  return xf_OK;

}

/******************************************************************/
//   FUNCTION Definition: xf_UnParallelizeLineSet
int
xf_UnParallelizeLineSet(xf_Mesh *Mesh_Glob, xf_Vector *C_Glob, 
                        xf_LineSet **LineSet_Glob, xf_Mesh *Mesh,
                        xf_LineSet *LineSet)
{
  int ierr, myRank, nProc, *Proc2nLine, *Line2nElem;
  int i, nLineTot, nelemtot, nelemtot_Glob, *LineElem, *LineEgrp, *sbuf;
  int *LineFace, j, e, elem, egrp, *Proc2nElem, *nElem, negrp, *sbuf2, *sbuf3;
  xf_Line Line;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0){
    //allocate global lineset
    ierr = xf_Error(xf_CreateLineSet(LineSet_Glob));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc((void**)&Proc2nLine, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc((void**)&Proc2nElem, nProc, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  
  /*We need the number of lines in each processor 
  to figure out the displacements*/
  ierr = xf_Error(xf_MPI_Gather((void*)&LineSet->nLine, (void *)Proc2nLine, 
                                1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void**)&sbuf, LineSet->nLine, 
                           sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0){
    nLineTot = 0;
    for (i=0; i < nProc; i++){ 
      nLineTot += Proc2nLine[i];
    }
    ierr = xf_Error(xf_Alloc((void**)&Line2nElem, nLineTot, 
                             sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  else {
    //number of elements in each line (1 per processor)
    ierr = xf_Error(xf_Alloc((void**)&Line2nElem, LineSet->nLine, 
                             sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  //proc 0 puts its data in the global array
  for (i = 0; i < LineSet->nLine; i++)
    sbuf[i] = LineSet->Line[i].nelem;
  //xf_MPI_Gatherv calculates the displacements for us
  ierr = xf_Error(xf_MPI_Gatherv((void*)sbuf, LineSet->nLine, 
                                 (void*)Line2nElem, Proc2nLine, 
                                 sizeof(int), 0));
  if (ierr != xf_OK) return ierr;

  xf_Release((void*)sbuf);
  
  nelemtot_Glob = nelemtot;
  //let all procs know what is the total number of elements
  ierr = xf_Error(xf_MPI_Allreduce(&nelemtot_Glob, 1, xfe_SizeInt, 
                                   xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;
  
  if (myRank == 0) {
    ierr = xf_Error(xf_Alloc((void**)&LineElem, nelemtot_Glob, 
                             sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc((void**)&LineEgrp, nelemtot_Glob, 
                             sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc((void**)&LineFace, nelemtot_Glob, 
                             sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  else {
    ierr = xf_Error(xf_Alloc((void**)&LineElem, nelemtot, 
                             sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc((void**)&LineEgrp, nelemtot, 
                             sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc((void**)&LineFace, nelemtot, 
                             sizeof(int)));
    if (ierr != xf_OK) return ierr;
  }
  ierr = xf_Error(xf_Alloc((void**)&sbuf, nelemtot, 
                           sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void**)&sbuf2, nelemtot, 
                           sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  ierr = xf_Error(xf_Alloc((void**)&sbuf3, nelemtot, 
                           sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  //unroll lines
  e=0;
  for (i=0; i < LineSet->nLine; i++) {
    for (j=0; j < LineSet->Line[i].nelem; j++) {
      egrp = LineSet->Line[i].egrp[j];
      elem = LineSet->Line[i].elem[j];
      sbuf[e] = Mesh->ParallelInfo->ElemLoc2Glob[egrp][elem];
      sbuf2[e] = egrp;
      sbuf3[e] = LineSet->Line[i].face[j];
      e++;
    }
  }
  
  //number of element in each proc
  ierr = xf_Error(xf_MPI_Gather((void*)&nelemtot, (void *)Proc2nElem, 
                                1*sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  //gather unrolled lines from procs
  ierr = xf_Error(xf_MPI_Gatherv((void*)sbuf, nelemtot, 
                                 (void*)LineElem, Proc2nElem, 
                                 sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Gatherv((void*)sbuf2, nelemtot, 
                                 (void*)LineEgrp, Proc2nElem, 
                                 sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_MPI_Gatherv((void*)sbuf3, nelemtot, 
                                 (void*)LineFace, Proc2nElem, 
                                 sizeof(int), 0));
  if (ierr != xf_OK) return ierr;
  
  xf_Release((void*)sbuf);
  xf_Release((void*)sbuf2);
  xf_Release((void*)sbuf3);
  
  if (myRank == 0){
    //put lines in LineSet_Glob
    (*LineSet_Glob)->nLine = nLineTot;
    ierr = xf_Error(xf_Alloc((void**)&(*LineSet_Glob)->Line, nLineTot, 
                             sizeof(xf_Line)));
    if (ierr != xf_OK) return ierr;
    j = 0;
    for (i=0; i<(*LineSet_Glob)->nLine; i++) {
      ierr = xf_Error(xf_AllocLine(&Line, Line2nElem[i]));
      if (ierr != xf_OK) return ierr;
      (*LineSet_Glob)->Line[i] = Line;
      for (e=0; e<Line2nElem[i]; e++) {
        (*LineSet_Glob)->Line[i].elem[e] = LineElem[j];
        (*LineSet_Glob)->Line[i].egrp[e] = LineEgrp[j];
        (*LineSet_Glob)->Line[i].face[e] = LineFace[j];
        j++;
      }
    }
    
    //allocate Elem2Line
    ierr = xf_Error(xf_GetnElem(Mesh_Glob, &nElem, NULL));
    if (ierr != xf_OK) return ierr;
    negrp = Mesh_Glob->nElemGroup;
    ierr = xf_Error(xf_VAlloc2((void ***)&(*LineSet_Glob)->Elem2Line, 
                               negrp, nElem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    xf_Release((void*)nElem);    
    //connect lines
    ierr = xf_Error(xf_LineConnection(Mesh_Glob, C_Glob, 
                                      (*LineSet_Glob)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_LineReNumber(Mesh_Glob, (*LineSet_Glob)));
    if (ierr != xf_OK) return ierr;
    
    xf_Release((void*)Proc2nLine);
    xf_Release((void*)Proc2nElem);
  }

  xf_Release((void*)Line2nElem);
  xf_Release((void*)LineElem);
  xf_Release((void*)LineEgrp);
  xf_Release((void*)LineFace);

  return xf_OK;
}
#if( UNIT_TEST==1 )
#include "xf_Line.test.in"
#endif
