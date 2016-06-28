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
 FILE:  xf_AllPull.c
 
 This file contains functions for pulling the All Structure
 
 */

#include "xf.h"
#include "xf_AllStruct.h"
#include "xf_All.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_Math.h"
#include "xf_Memory.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Residual.h"
#include "xf_Solver.h"
#include "xf_AdaptStruct.h"
#include "xf_Mesh.h"
#include "xf_MeshStruct.h"
#include "xf_MeshTools.h"
#include "xf_Basis.h"
#include "xf_EqnSet.h"
#include "xf_AllPull.h"
#include "xf_AllPullStruct.h"
#include "xf_MPI.h"

// turn on to 1 for debugging
#define ah_DEBUG 0


/******************************************************************/
//   FUNCTION Definition: xf_CheckExist
int
xf_CheckExist(int num, int size, int *vector, int *index)
{
	
  /*
   
   PURPOSE: 
   
   Check existance of num in vector.
   
   INPUTS:
   
   num : value to be checked
   size : size of vector
   vector : vector of integers among which num will be checked
   out: -1 -> vector does not have num; >=0 -> vector has num at index
   
   OUTPUTS: 
   None. value of index is changed
   
   RETURNS: Error code		
   */
  int i;
	
  (*index) = -1;
	
  for (i = 0; i < size; i++){ 
    if (num == vector[i]){
      (*index) = i; 
      break;
    }
  }
	
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateMeshPullInfo
static int
xf_CreateMeshPullInfo(xf_All *All, xf_MeshPullInfo **pMeshPullInfo)
{
  int ierr, i, nGrp;
  xf_Mesh *Mesh;
  
  Mesh = All->Mesh;
  
  //allocate structure space
  ierr = xf_Error(xf_Alloc((void **)&(*pMeshPullInfo), 1, sizeof(xf_MeshPullInfo)));
  if (ierr != xf_OK) return ierr;
  
  //setting default values
  (*pMeshPullInfo)->Dim = Mesh->Dim;
  
  /* Nodes */
  (*pMeshPullInfo)->Nodes.nItem = 0;//initial value
  (*pMeshPullInfo)->Nodes.Item = NULL;
  
  /* IFaces */
  (*pMeshPullInfo)->IFaces.nItem = 0;//initial value
  (*pMeshPullInfo)->IFaces.Item = NULL;
  
  /* BFaces */
  (*pMeshPullInfo)->nBFaceGroup = Mesh->nBFaceGroup;
  
  ierr = xf_Error(xf_Alloc((void **)&((*pMeshPullInfo)->BFaceGroups), 
                           (*pMeshPullInfo)->nBFaceGroup, sizeof(xf_ConvTable)));
  if (ierr != xf_OK) return ierr;
  
  for (i = 0; i < (*pMeshPullInfo)->nBFaceGroup; i++){
    //number of faces in each group: initialize it to 0
    (*pMeshPullInfo)->BFaceGroups[i].nItem = 0;
    (*pMeshPullInfo)->BFaceGroups[i].Item = NULL;
  }
  
  /* Elems */
  (*pMeshPullInfo)->nElemGroup = Mesh->nElemGroup;

  // account for possible halo element groups
  if (Mesh->ParallelInfo != NULL) (*pMeshPullInfo)->nElemGroup *= 2;
  
  ierr = xf_Error(xf_Alloc((void **)&((*pMeshPullInfo)->ElemGroups), 
                           (*pMeshPullInfo)->nElemGroup, sizeof(xf_ConvTable)));
  if (ierr != xf_OK) return ierr;
  
  for (i = 0; i < (*pMeshPullInfo)->nElemGroup; i++){
    //number of faces in each group: initialize it to 0
    (*pMeshPullInfo)->ElemGroups[i].nItem = 0;
    (*pMeshPullInfo)->ElemGroups[i].Item = NULL;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UpdateConvTable
static int
xf_UpdateConvTable(xf_ConvTable **pConvTable, int Item)
{
  int ierr, index, nItem;
  
  //is the Item already listed?
  ierr = xf_Error(xf_CheckExist(Item, (*pConvTable)->nItem, 
                                (*pConvTable)->Item, &index));
  if (ierr != xf_OK) return ierr;
  
  if (index < 0){//not listed yet
    (*pConvTable)->nItem++;
    
    nItem = (*pConvTable)->nItem;
    
    //adjust the size of the conversion table
    ierr = xf_Error(xf_ReAlloc((void **)&((*pConvTable)->Item), 
                               nItem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    
    //update the table
    (*pConvTable)->Item[nItem-1] = Item;
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UpdateFaceConvTable
static int
xf_UpdateFaceConvTable(xf_Face Face, xf_MeshPullInfo **pMeshPullInfo)
{
  int ierr, fgrp, face, index, nItem;
  xf_MeshPullInfo *MeshPullInfo;
  xf_ConvTable *FaceConvTable;
  
  fgrp = Face.Group;
  face = Face.Number;
  MeshPullInfo = (*pMeshPullInfo);
  
  //Error Check
  if (fgrp == xf_NULLFACE){
    //    xf_printf("A NULLFACE should not be pulled\n");
    //    return xf_MESH_ERROR;
    return xf_OK; // no harm in pulling a NULLFACE
  }
  
  //to which conversion table should we point?
  if (fgrp == xf_INTERIORFACE)
    FaceConvTable = &(MeshPullInfo->IFaces);
  else{
    if (fgrp >= MeshPullInfo->nBFaceGroup)
      return xf_INPUT_ERROR;
    
    FaceConvTable = &(MeshPullInfo->BFaceGroups[fgrp]);
  }
  
  //update table
  ierr = xf_Error(xf_UpdateConvTable(&FaceConvTable, face));
  if (ierr != xf_OK) return ierr;
  
  FaceConvTable = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_UpdateElemConvTable
static int
xf_UpdateElemConvTable(int egrp, int elem, xf_MeshPullInfo **pMeshPullInfo)
{
  int ierr, index, nItem;
  xf_MeshPullInfo *MeshPullInfo;
  xf_ConvTable *ElemConvTable;
  
  MeshPullInfo = (*pMeshPullInfo);
  
  if (egrp >= MeshPullInfo->nElemGroup)
    return xf_INPUT_ERROR;
  
  //point to conversion table of egrp
  ElemConvTable = MeshPullInfo->ElemGroups+egrp;
  
  ierr = xf_Error(xf_UpdateConvTable(&ElemConvTable, elem));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_PullHangFacesFromFace
static int
xf_PullHangFacesFromFace(xf_Mesh *Mesh, int egrp, int elem, int face, 
                         xf_MeshPullInfo **pMeshPullInfo)
{
  int ierr, nface0, face0, iface, nface, hang;
  int egrpL, elemL, egrpR, elemR;
  xf_Face Face;
  
  // number of faces in the original group
  ierr = xf_Error(xf_Basis2nFace(Mesh->ElemGroup[egrp].QBasis, &nface0));
  if (ierr != xf_OK) return ierr;
  
  nface = Mesh->ElemGroup[egrp].nFace[elem];
  
  if (face >= nface)
    return xf_INPUT_ERROR;
  
  for (iface = 0; iface < nface; iface++){
    Face = Mesh->ElemGroup[egrp].Face[elem][iface];
    
    // skip null faces (those that connect halo elements to each other)
    if (Face.Group == xf_NULLFACE) continue;

    ierr = xf_Error(xf_CheckHangFace(Mesh, egrp, elem, iface, NULL, &face0, NULL, NULL));
    if (ierr != xf_OK) return ierr;
    
    if (face0 == face){//check if face is part of the original face
      //add this face to conversion table
      ierr = xf_Error(xf_UpdateFaceConvTable(Face, &(*pMeshPullInfo)));
      if (ierr != xf_OK) return ierr;
      
      //add elements to conversion table
      egrpL = Mesh->IFace[Face.Number].ElemGroupL;
      elemL = Mesh->IFace[Face.Number].ElemL;

      ierr = xf_Error(xf_UpdateElemConvTable(egrpL, elemL, &(*pMeshPullInfo)));
      if (ierr != xf_OK) return ierr;
      
      egrpR = Mesh->IFace[Face.Number].ElemGroupR;
      elemR = Mesh->IFace[Face.Number].ElemR;
      
      ierr = xf_Error(xf_UpdateElemConvTable(egrpR, elemR, &(*pMeshPullInfo)));
      if (ierr != xf_OK) return ierr;
    }
  }//iface
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteNodes2AllSmall
static int
xf_WriteNodes2AllSmall(xf_All *All, xf_All *All_Small, 
                       xf_MeshPullInfo *MeshPullInfo)
{
  int ierr, inode, i, j;
  
  All_Small->Mesh->nNode = MeshPullInfo->Nodes.nItem;
	
  ierr = xf_Error(xf_Alloc2((void ***) &All_Small->Mesh->Coord, 
                            All_Small->Mesh->nNode, 
                            All_Small->Mesh->Dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  for (i = 0; i < All_Small->Mesh->nNode; i++){
    for(j = 0; j < All_Small->Mesh->Dim; j++){ 
      inode = MeshPullInfo->Nodes.Item[i];
      All_Small->Mesh->Coord[i][j] = All->Mesh->Coord[inode][j];
    }
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteElems2AllSmall
static int
xf_WriteElems2AllSmall(xf_All *All, xf_All *All_Small, 
                       xf_MeshPullInfo *MeshPullInfo)
{
  int ierr, egrp, elem, elem_glob, node, node_glob, node_index;
  int face;
  xf_Face Face_glob;
  xf_Mesh *Mesh_Small, *Mesh;
  
  Mesh_Small = All_Small->Mesh;
  Mesh = All->Mesh;
  
  Mesh_Small->nElemGroup = MeshPullInfo->nElemGroup;
  
  //allocate groups
  ierr = xf_Error(xf_Alloc((void **)&Mesh_Small->ElemGroup, 
                           Mesh_Small->nElemGroup, sizeof(xf_ElemGroup)));
  if (ierr != xf_OK) return ierr;
  
  for (egrp = 0; egrp < Mesh_Small->nElemGroup; egrp++){
    if (MeshPullInfo->ElemGroups[egrp].nItem == 0){
      //CutFlag, CutElemData are set to NULL
      xf_InitElemGroup(&(Mesh_Small->ElemGroup[egrp]));
      //QBasis, QOrder, nNode
      Mesh_Small->ElemGroup[egrp].QBasis = Mesh->ElemGroup[egrp].QBasis;
      Mesh_Small->ElemGroup[egrp].QOrder = Mesh->ElemGroup[egrp].QOrder;
      Mesh_Small->ElemGroup[egrp].nNode = Mesh->ElemGroup[egrp].nNode;
    }
    else{
      Mesh_Small->ElemGroup[egrp] = Mesh->ElemGroup[egrp]; 
      
      //nElem
      Mesh_Small->ElemGroup[egrp].nElem = MeshPullInfo->ElemGroups[egrp].nItem;
      //allocating nFace
      ierr = xf_Error(xf_Alloc((void **)&Mesh_Small->ElemGroup[egrp].nFace, 
                               Mesh_Small->ElemGroup[egrp].nElem, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      //allocating **Node
      ierr = xf_Error(xf_Alloc2((void ***)&Mesh_Small->ElemGroup[egrp].Node, 
                                Mesh_Small->ElemGroup[egrp].nElem, 
                                Mesh_Small->ElemGroup[egrp].nNode, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      
      for (elem = 0; elem < Mesh_Small->ElemGroup[egrp].nElem; elem++){
        //global element number
        elem_glob = MeshPullInfo->ElemGroups[egrp].Item[elem];
        //nFace
        Mesh_Small->ElemGroup[egrp].nFace[elem] = Mesh->ElemGroup[egrp].nFace[elem_glob];
        //Node
        for (node = 0; node < Mesh_Small->ElemGroup[egrp].nNode; node++){
          node_glob = Mesh->ElemGroup[egrp].Node[elem_glob][node];
          //get index in local mesh
          ierr = xf_Error(xf_CheckExist(node_glob, MeshPullInfo->Nodes.nItem, 
                                        MeshPullInfo->Nodes.Item, &node_index));
          if (ierr != xf_OK) return ierr;
          
          Mesh_Small->ElemGroup[egrp].Node[elem][node] = node_index;
        }
      }
      
      //allocating **Face structures
      ierr = xf_Error(xf_VAlloc2((void ***)&Mesh_Small->ElemGroup[egrp].Face,
                                 Mesh_Small->ElemGroup[egrp].nElem, 
                                 Mesh_Small->ElemGroup[egrp].nFace, sizeof(xf_Face)));
      if (ierr != xf_OK) return ierr;
      //this loop is for the faces
      for (elem = 0; elem < Mesh_Small->ElemGroup[egrp].nElem; elem++){
        //loop through faces and translate information
        for (face = 0; face < Mesh_Small->ElemGroup[egrp].nFace[elem]; face++){
          //For now the faces will be set to be NULL. This will be fixed in the Faces Loop.
          Mesh_Small->ElemGroup[egrp].Face[elem][face].Group = xf_NULLFACE;
          Mesh_Small->ElemGroup[egrp].Face[elem][face].Number = 0;
        }
      }
    }
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteFaces2AllSmall
static int
xf_WriteFaces2AllSmall(xf_All *All, xf_All *All_Small, 
                       xf_MeshPullInfo *MeshPullInfo)
{
  int ierr, fgrp, face, face_glob;
  int ElemGroupL, ElemGroupR, FaceL, FaceR, ElemL, ElemR;
  xf_Mesh *Mesh_Small, *Mesh;
  
  Mesh_Small = All_Small->Mesh;
  Mesh = All->Mesh;
  
  //write internal faces
  Mesh_Small->nIFace = MeshPullInfo->IFaces.nItem;
  
  ierr = xf_Error(xf_Alloc((void **)&Mesh_Small->IFace, 
                           Mesh_Small->nIFace, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;
  
  for (face = 0; face < Mesh_Small->nIFace; face++){
    face_glob = MeshPullInfo->IFaces.Item[face];
    /*copy all values from global structure and change 
     only the pertinent to all_small*/
    Mesh_Small->IFace[face] = Mesh->IFace[face_glob];
    
    ElemGroupL = Mesh->IFace[face_glob].ElemGroupL;//no need to convert this
    ElemL = Mesh->IFace[face_glob].ElemL;
    FaceL = All->Mesh->IFace[face_glob].FaceL;//no need to convert this
		
    ierr = xf_Error(xf_CheckExist(ElemL, MeshPullInfo->ElemGroups[ElemGroupL].nItem, 
                                  MeshPullInfo->ElemGroups[ElemGroupL].Item, &ElemL));
    if (ierr != xf_OK) return ierr;
    //error check
    if (ElemL < 0){
      xf_printf("\nElemL should exist in MeshPullInfo\n");
      return xf_NOT_FOUND;
		}
    
    ElemGroupR = Mesh->IFace[face_glob].ElemGroupR;//no need to convert this
    ElemR = Mesh->IFace[face_glob].ElemR;
    FaceR = All->Mesh->IFace[face_glob].FaceR;//no need to convert this
		
    ierr = xf_Error(xf_CheckExist(ElemR, MeshPullInfo->ElemGroups[ElemGroupR].nItem, 
                                  MeshPullInfo->ElemGroups[ElemGroupR].Item, &ElemR));
    if (ierr != xf_OK) return ierr;
		//error check
    if (ElemR < 0){
      xf_printf("\nElemR should exist in MeshPullInfo\n");
      return xf_NOT_FOUND;
		}
    
    //Putting values in All_Small
    All_Small->Mesh->IFace[face] = All->Mesh->IFace[face_glob];
    
    All_Small->Mesh->IFace[face].ElemL = ElemL;
    All_Small->Mesh->IFace[face].ElemGroupL = ElemGroupL;
    All_Small->Mesh->IFace[face].FaceL = FaceL;
    //correcting the face number and group from the elem writting function
    All_Small->Mesh->ElemGroup[ElemGroupL].Face[ElemL][FaceL].Number = face;
    All_Small->Mesh->ElemGroup[ElemGroupL].Face[ElemL][FaceL].Group = xf_INTERIORFACE;
    
    All_Small->Mesh->IFace[face].ElemR = ElemR;
    All_Small->Mesh->IFace[face].ElemGroupR = ElemGroupR;		
    All_Small->Mesh->IFace[face].FaceR = FaceR;
    //correcting the face number and group from the elem writting function
    All_Small->Mesh->ElemGroup[ElemGroupR].Face[ElemR][FaceR].Number = face;
    All_Small->Mesh->ElemGroup[ElemGroupR].Face[ElemR][FaceR].Group = xf_INTERIORFACE;
  }
  
  //write boundary faces
  Mesh_Small->nBFaceGroup = MeshPullInfo->nBFaceGroup;
  
  ierr = xf_Error(xf_Alloc((void **)&Mesh_Small->BFaceGroup, 
                           Mesh_Small->nBFaceGroup, sizeof(xf_BFaceGroup)));
  if (ierr != xf_OK) return ierr;
  
  for (fgrp = 0; fgrp < Mesh_Small->nBFaceGroup; fgrp++){
    //Title: allocate and copy
    ierr = xf_Error(xf_AllocString(&Mesh_Small->BFaceGroup[fgrp].Title, 
                                   xf_MAXSTRLEN, Mesh->BFaceGroup[fgrp].Title));
    if (ierr != xf_OK) return ierr;
    
    //nBFace
    Mesh_Small->BFaceGroup[fgrp].nBFace = MeshPullInfo->BFaceGroups[fgrp].nItem;
    
    //BFace
    ierr = xf_Error(xf_Alloc((void **)&Mesh_Small->BFaceGroup[fgrp].BFace, 
                             Mesh_Small->BFaceGroup[fgrp].nBFace, 
                             sizeof(xf_BFace)));
    if (ierr != xf_OK) return ierr;
    
    for (face = 0; face < Mesh_Small->BFaceGroup[fgrp].nBFace; face++){
      face_glob = MeshPullInfo->BFaceGroups[fgrp].Item[face];
      
      //Orient and CutFaceData is set
      xf_InitBFace(Mesh_Small->BFaceGroup[fgrp].BFace+face);
      
      ElemGroupL = Mesh->BFaceGroup[fgrp].BFace[face_glob].ElemGroup;
      ElemL = Mesh->BFaceGroup[fgrp].BFace[face_glob].Elem;
      FaceL = Mesh->BFaceGroup[fgrp].BFace[face_glob].Face;
      //local mesh element number
      ierr = xf_Error(xf_CheckExist(ElemL, MeshPullInfo->ElemGroups[ElemGroupL].nItem, 
                                    MeshPullInfo->ElemGroups[ElemGroupL].Item, &ElemL));
      if (ierr != xf_OK) return ierr;
      //error check
      if (ElemL < 0){
        xf_printf("\nElem should exist in MeshPullInfo\n");
        return xf_NOT_FOUND;
      }
      
      
      /*copy all values from global structure and change 
       only the pertinent to all_small*/
      Mesh_Small->BFaceGroup[fgrp].BFace[face] = Mesh->BFaceGroup[fgrp].BFace[face_glob];
      
      Mesh_Small->BFaceGroup[fgrp].BFace[face].ElemGroup = ElemGroupL;
      Mesh_Small->BFaceGroup[fgrp].BFace[face].Elem = ElemL;
      Mesh_Small->BFaceGroup[fgrp].BFace[face].Face = FaceL;
      Mesh_Small->ElemGroup[ElemGroupL].Face[ElemL][FaceL].Group = fgrp;
      Mesh_Small->ElemGroup[ElemGroupL].Face[ElemL][FaceL].Number = face;
    }
    
  }
  
  Mesh = NULL;
  Mesh_Small = NULL;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteData2AllSmall
static int
xf_WriteData2AllSmall(xf_All *All, xf_All *All_Small, 
                      xf_MeshPullInfo *MeshPullInfo)
{
  int ierr, egrp, elem, elem_glob, r, nr,myRank,nProc;
  int *nComp, **vOrder, *rvec;
  xf_Data *D, *DSmall;
  xf_Vector *V, *VSmall;
  
  nComp = rvec = NULL;
  vOrder = NULL;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  if (All->DataSet->Head != NULL){
    D = All->DataSet->Head;
    while (D != NULL){
      if (D->Type == xfe_Vector){
        V = (xf_Vector *) D->Data;
        if ((V->Linkage == xfe_LinkageGlobElem &&
            V->SolverRole != xfe_SolverRoleNone)
            || D->ReadWrite == xfe_True){

          if (V->nComp != NULL && V->vOrder != NULL){//variable order
            if (All_Small->Mesh->nElemGroup != V->nArray) 
              return xf_Error(xf_CODE_LOGIC_ERROR);
            
            ierr = xf_Error(xf_Alloc((void **)&nComp, 
                                     All_Small->Mesh->nElemGroup, 
                                     sizeof(int)));
            if (ierr != xf_OK) return ierr;
            
            ierr = xf_Error(xf_Alloc((void **)&rvec, 
                                     All_Small->Mesh->nElemGroup, 
                                     sizeof(int)));
            if (ierr != xf_OK) return ierr;
            
            for (egrp = 0; egrp < All_Small->Mesh->nElemGroup; egrp++){
              nComp[egrp] = All_Small->Mesh->ElemGroup[egrp].nElem;
              rvec[egrp] = V->GenArray[egrp].r;
            }
            
            ierr = xf_Error(xf_VAlloc2((void ***)&vOrder, 
                                       All_Small->Mesh->nElemGroup, 
                                       nComp, sizeof(int)));
            if (ierr != xf_OK) return ierr;
            for (egrp = 0; egrp < All_Small->Mesh->nElemGroup; egrp++){
              for (elem = 0; elem < nComp[egrp]; elem++){
                elem_glob = MeshPullInfo->ElemGroups[egrp].Item[elem];
                vOrder[egrp][elem] = V->vOrder[egrp][elem_glob];
              }
            }
          }
                    
          //Creating vector to copy data to
          ierr = xf_Error(xf_FindVector(All_Small, D->Title, V->Linkage, 
                                        V->StateRank, V->StateName, V->TimeIndex, 
                                        V->MGIndex, V->Basis, V->Order, 
                                        nComp, vOrder, rvec, V->Size, xfe_False, 
                                        xfe_True, &DSmall, &VSmall, NULL));
          if (ierr != xf_OK) return ierr;

          xf_Release((void *)nComp);
          xf_Release((void *)rvec);
          xf_Release2((void **)vOrder);
          nComp = rvec = NULL;
          vOrder = NULL;
					
          DSmall->ReadWrite = D->ReadWrite;
          VSmall->SolverRole = V->SolverRole;
          VSmall->ParallelFlag = xfe_False;
					
          //looping through the arrays in V
          for (egrp = 0; egrp < All_Small->Mesh->nElemGroup; egrp++){
            for (elem = 0; elem < All_Small->Mesh->ElemGroup[egrp].nElem; elem++){
              //Global number of the element
              elem_glob = MeshPullInfo->ElemGroups[egrp].Item[elem];
              if (VSmall->GenArray[egrp].vr != NULL){ 
                //VSmall->GenArray[egrp].vr[elem] = V->GenArray[egrp].vr[elem_glob];
                //VSmall->GenArray[egrp].r = V->GenArray[egrp].r;
                nr = VSmall->GenArray[egrp].vr[elem];
              }
              else{
                //VSmall->GenArray[egrp].r = V->GenArray[egrp].r;
                nr = VSmall->GenArray[egrp].r;
              }
              for (r = 0; r < nr; r++) {
                if (VSmall->GenArray[egrp].Size == xfe_SizeReal)
                  VSmall->GenArray[egrp].rValue[elem][r] = V->GenArray[egrp].rValue[elem_glob][r];
                else if (VSmall->GenArray[egrp].Size == xfe_SizeInt)
                  VSmall->GenArray[egrp].iValue[elem][r] = V->GenArray[egrp].iValue[elem_glob][r];
              }
            }
          }
        }
      }
      D = D->Next;
    }
  }
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyMeshPullInfo
static int
xf_DestroyMeshPullInfo(xf_MeshPullInfo *MeshPullInfo)
{
  int ierr, i;
  
  if (MeshPullInfo != NULL){
    xf_Release((void *)MeshPullInfo->Nodes.Item);
    xf_Release((void *)MeshPullInfo->IFaces.Item);
    
    for (i = 0; i < MeshPullInfo->nBFaceGroup; i++)
      xf_Release((void *)MeshPullInfo->BFaceGroups[i].Item);
    
    xf_Release((void *)MeshPullInfo->BFaceGroups);
    
    for (i = 0; i < MeshPullInfo->nElemGroup; i++)
      xf_Release((void *)MeshPullInfo->ElemGroups[i].Item);
    xf_Release((void *)MeshPullInfo->ElemGroups);
  }
  
  xf_Release((void *)MeshPullInfo);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteMeshPullInfo
static int
xf_WriteMeshPullInfo(xf_MeshPullInfo *MeshPullInfo, char *filename)
{
  int i, j;
  FILE *fid;
  
  fid = fopen(filename, "w");
  
  fprintf(fid, "====MeshPullInfo====\n");
  
  //dimension
  fprintf(fid, "Mesh Dimension: %d\n\n",MeshPullInfo->Dim);
  
  //nodes
  fprintf(fid, "Number of Nodes: %d\n\n",MeshPullInfo->Nodes.nItem);
  for (i = 0; i < MeshPullInfo->Nodes.nItem; i++) {
    fprintf(fid, "Node[%d] => %d\n", i, MeshPullInfo->Nodes.Item[i]);
  }
  //ifaces
  fprintf(fid, "Number of IFaces: %d\n\n",MeshPullInfo->IFaces.nItem);
  for (i = 0; i < MeshPullInfo->IFaces.nItem; i++) {
    fprintf(fid, "IFace[%d] => %d\n", i, MeshPullInfo->IFaces.Item[i]);
  }
  //bfaces
  fprintf(fid, "Number of BFaceGroups: %d\n\n",MeshPullInfo->nBFaceGroup);
  for (i = 0; i < MeshPullInfo->nBFaceGroup; i++){
    fprintf(fid, "Number of BFaces in group %d: %d\n", i, MeshPullInfo->BFaceGroups[i].nItem);
    for (j = 0; j < MeshPullInfo->BFaceGroups[i].nItem; j++){
      fprintf(fid, "\tBFace[%d] => %d\n", j, MeshPullInfo->BFaceGroups[i].Item[j]);
    }
  }
  //elemgroups
  fprintf(fid, "Number of ElemGroups: %d\n\n",MeshPullInfo->nElemGroup);
  for (i = 0; i < MeshPullInfo->nElemGroup; i++){
    fprintf(fid, "Number of Elems in group %d: %d\n", i, MeshPullInfo->ElemGroups[i].nItem);
    for (j = 0; j < MeshPullInfo->ElemGroups[i].nItem; j++){
      fprintf(fid, "\tElem[%d] => %d\n", j, MeshPullInfo->ElemGroups[i].Item[j]);
    }
  }
  
  fclose(fid);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_AllPull
int
xf_AllPull(xf_All *All, xf_All *All_Small, int egrp, int elem)
{
  /*
   
   PURPOSE: 
   
   Copy the local mesh (elements and 1st-level neighbors) to All_Small.
   Point the EqnSet to the global EqnSet.
   Copy global data to local data.
   
   INPUTS:
   
   All : All Structure
   All_Small : All_Small structure
   elem : element to pull local mesh structure from
   egrp : element group to which elem belongs
   
   OUTPUTS: 
   None: All_Small structure is modified.
   
   
   RETURNS: Error Code
   
   */
  int ierr, i, j, k, nFace, fgrp, face, hang, n, myRank, nProc;
  int nbregrp, nbrelem, nbrface, face0, igrp, ielem, e, g, inode, iface; 
  int ElemGroupL, ElemGroupR, FaceL, FaceR, ElemL, ElemR, OrientL, OrientR;
  char fname[xf_MAXSTRLEN];
  xf_MeshPullInfo *MeshPullInfo;
  xf_ConvTable *NodesConvTable;
  xf_Face Face;
  
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;
  
  //Pulling mesh information
  //Centralizing all the pulling information
  ierr = xf_Error(xf_CreateMeshPullInfo(All, &MeshPullInfo));
  if (ierr != xf_OK) return ierr;
  
  //adding central element to conversion table
  //this guarantees that the central elem is number 0 in egrp in All_Small
  ierr = xf_Error(xf_UpdateElemConvTable(egrp, elem, &MeshPullInfo));
  if (ierr != xf_OK) return ierr;
  
  nFace = All->Mesh->ElemGroup[egrp].nFace[elem];
  
  //getting the faces of the central element, elements and element groups
  for (iface = 0;iface < nFace; iface++){
    //global face number
    Face = All->Mesh->ElemGroup[egrp].Face[elem][iface];
    fgrp = Face.Group;
    face = Face.Number;
		
    //interior faces
    if(fgrp == xf_INTERIORFACE){
      //Check if it is a "hanging" face
      hang = All->Mesh->IFace[face].HangNumber;
      //Note: pull all the hanging faces that have the same original face
      if (hang != 0){ 
        if (hang > 0){
          //guarantee to look from the coarse side
          nbregrp = All->Mesh->IFace[face].ElemGroupR;
          nbrelem = All->Mesh->IFace[face].ElemR;
          nbrface = All->Mesh->IFace[face].FaceR;
        }
        else {
          //guarantee to look from the coarse side
          nbregrp = All->Mesh->IFace[face].ElemGroupL;
          nbrelem = All->Mesh->IFace[face].ElemL;
          nbrface = All->Mesh->IFace[face].FaceL;
        }
        //get original face
        ierr = xf_Error(xf_CheckHangFace(All->Mesh, nbregrp, 
                                         nbrelem, nbrface, NULL, 
                                         &face0, NULL, NULL));
        if (ierr != xf_OK) return ierr;
        /*add all hanging faces and correspondent elements that are 
         part of the original face*/
        ierr = xf_Error(xf_PullHangFacesFromFace(All->Mesh, nbregrp, 
                                                 nbrelem, face0, 
                                                 &MeshPullInfo));
        if (ierr != xf_OK) return ierr;
      }
      else {//add face and elements to conversion table
        //add face to conversion table
        ierr = xf_Error(xf_UpdateFaceConvTable(Face, &MeshPullInfo));
        if (ierr != xf_OK) return ierr;
        
        //add neighbor to conversion table
        ierr = xf_Error(xf_NeighborAcrossFace(All->Mesh, egrp, elem, iface, 
                                              &nbregrp, &nbrelem,NULL));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_UpdateElemConvTable(nbregrp, nbrelem, &MeshPullInfo));
        if (ierr != xf_OK) return ierr;
        
      }
    }
    else if (fgrp != xf_NULLFACE){
      //only add face to conversion table
      ierr = xf_Error(xf_UpdateFaceConvTable(Face, &MeshPullInfo));
      if (ierr != xf_OK) return ierr;
    }
  }
  
  //getting the nodes and neighbors' boundary faces
  NodesConvTable = &(MeshPullInfo->Nodes);
  for (igrp = 0; igrp < MeshPullInfo->nElemGroup; igrp++){
    for (ielem = 0; ielem < MeshPullInfo->ElemGroups[igrp].nItem; ielem++){
      //global element number
      e = MeshPullInfo->ElemGroups[igrp].Item[ielem];
      //nodes loop
      for (inode = 0; inode < All->Mesh->ElemGroup[igrp].nNode; inode++){
        //global node number
        n = All->Mesh->ElemGroup[igrp].Node[e][inode];
        //add node to conversion table
        
        ierr = xf_Error(xf_UpdateConvTable(&NodesConvTable, n));
        if (ierr != xf_OK) return ierr;
      }
      //faces loop
      for (iface = 0; iface < All->Mesh->ElemGroup[igrp].nFace[e]; iface++){
        Face = All->Mesh->ElemGroup[igrp].Face[e][iface];
        if (Face.Group >= 0){
          ierr = xf_Error(xf_UpdateFaceConvTable(Face, &MeshPullInfo));
          if (ierr != xf_OK) return ierr;
        }
      }
    }
  }
  NodesConvTable = NULL;
  //order node conversion table for consistency
  ierr = xf_Error(xf_SortInt(MeshPullInfo->Nodes.nItem, 
                             MeshPullInfo->Nodes.Item));
  if (ierr != xf_OK) return ierr;
  
  //print MeshPullInfo
  //sprintf(fname,"MeshPullInfo_%d-%d_%d_%d.txt",myRank,nProc,egrp,elem);
  //xf_WriteMeshPullInfo(MeshPullInfo, fname);
  /****************************************************************************/
  //Writing info into All_Small
	All_Small->Mesh->ParallelInfo = NULL;
  All_Small->Mesh->Dim = MeshPullInfo->Dim;	
	
  //Writing Nodes
  ierr = xf_Error(xf_WriteNodes2AllSmall(All, All_Small, MeshPullInfo));
  if (ierr != xf_OK) return ierr;
  
  //Writing Elements
  ierr = xf_Error(xf_WriteElems2AllSmall(All, All_Small, MeshPullInfo));
  if (ierr != xf_OK) return ierr;
  
  //Writing Faces
  ierr = xf_Error(xf_WriteFaces2AllSmall(All, All_Small, MeshPullInfo));
  if (ierr != xf_OK) return ierr;
	
  //Pulling EqnSet								
  ierr = xf_Error(xf_DestroyEqnSet(All_Small->EqnSet, xfe_True));
  if (ierr != xf_OK) return ierr;
	
  All_Small->EqnSet = All->EqnSet;  //Only pointing the equation set.
	
  //Writing data
  ierr = xf_Error(xf_WriteData2AllSmall(All, All_Small, MeshPullInfo));
  if (ierr != xf_OK) return ierr;
  
  //Pulling Parameters								
  ierr = xf_Error(xf_CopyParam(All->Param,All_Small->Param));
  if (ierr != xf_OK) return ierr;
  
  //Releasing the conversion tables memory
	ierr = xf_Error(xf_DestroyMeshPullInfo(MeshPullInfo));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


#if( UNIT_TEST==1 )
#include "xf_AllPull.test.in"
#endif

