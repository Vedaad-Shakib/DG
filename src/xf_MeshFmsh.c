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
 FILE:  xf_MeshFmsh.c
 
 This file contains functions for working with the Fluent msh format.
 
 */

#include "xf_MeshStruct.h"
#include "xf_Math.h"

enum xfe_ZoneType {
  xfe_FluentInterior,
  xfe_FluentPeriodic,
  xfe_FluentPeriodicShadow,
  xfe_FluentOther
}; /* special treatment for interior and periodic faces, 
    "other" is treated as a boundary*/

typedef struct
{
  int nNode;
  /* Number of nodes */
  
  int nFace;
  /* Number of faces in this face group */
  
  int **Node;
  /* Nodes composing each face Node[0->nFace-1][0->nNode-1] */
  
  enum xfe_BasisType QBasis; 
  
  int First, Last;
  /* amplitude of indices in this face group */
  
  int ZoneId;
  /* zone id in msh file */
  
  int *cr, *cl;
  /* right and left cells for each face */
  
  enum xfe_ZoneType ZoneType;
  /* Zone Type for this face group */
  
  char *Title;
  /* zone name */
}
xf_FluentFaceGroup; //Please refer to Appendix B of Fluent manual

typedef struct
{
  int First, Last;
  /* amplitude of indices in this cell group */
  
  int nFace;
  
  enum xfe_BasisType QBasis; 
}
xf_FluentCellGroup;

typedef struct
{
  int nP;                  /* number of pairs */
  int nPprealloc;          /* preallocation size */
  real average_distance;   /* average cartesian distance between periodic nodes */
  int *im_node;            /* image node number */
  int *sh_node;            /* shadow node number */
  enum xfe_PeriodicityType Periodicity;
  enum xfe_Bool MergeNodes;/* if true, nodes from this pair will be merged */
}
xf_FluentPeriodicNodePair;

typedef struct
{
  int nP;                  /* number of pairs */
  int *im_fgrp;            /* image face group */
  int *im_face;            /* image face number */
  int *sh_fgrp;            /* shadow face group */
  int *sh_face;            /* shadow face number */
}
xf_FluentPeriodicFacePair;

/******************************************************************/
//  FUNCTION Definition: xf_FluentElemType2QBasis
static int 
xf_FluentElemType2QBasis(int elem_type, enum xfe_BasisType *QBasis)
{
  switch (elem_type) {
    case 1:
      (*QBasis) = xfe_TriLagrange;
      break;
    case 2:
      (*QBasis) = xfe_TetLagrange;
      break;
    case 3:
      (*QBasis) = xfe_QuadLagrange;
      break;
    case 4:
      (*QBasis) = xfe_HexLagrange;
      break;
    default:
      return xf_NOT_SUPPORTED;
      break;
  }
  return xf_OK;
}

/******************************************************************/
//  FUNCTION Definition: xf_FluentElemType2QBasis
static int 
xf_FluentFaceType2QBasis(int face_type, enum xfe_BasisType *QBasis)
{
  switch (face_type) {
    case 2:
      (*QBasis) = xfe_SegLagrange;
      break;
    case 3:
      (*QBasis) = xfe_TriLagrange;
      break;
    case 4:
      (*QBasis) = xfe_QuadLagrange;
      break;
    default:
      return xf_NOT_SUPPORTED;
      break;
  }
  return xf_OK;
}

/******************************************************************/
//  FUNCTION Definition: xf_FluentBC2ZoneType
static int 
xf_FluentBC2ZoneType(int bc_type, enum xfe_ZoneType *ZoneType)
{
  switch (bc_type) {
    case 2:
      (*ZoneType) = xfe_FluentInterior;
      break;
    case 12:
      (*ZoneType) = xfe_FluentPeriodic;
      break;
    case 8:
      (*ZoneType) = xfe_FluentPeriodicShadow;
      break;
    default:
      (*ZoneType) = xfe_FluentOther;
      break;
  }
  return xf_OK;
}


/******************************************************************/
//  FUNCTION Definition: xf_Hex2DecInt
static int 
xf_Hex2DecInt(char hexstr[])
{
  int i, d, dec=0, base=1;
  char c;
  
  
  for (i = strlen(hexstr)-1; i >= 0; i--){
    c = hexstr[i];
    if ((c == ' ') || (c == '\n')) continue; // skip trailing blanks
    if ((c >= '0') && (c <= '9')) d = c - '0';
    else if ((c >= 'a') && (c <= 'f')) d = c - 'a' + 10;
    else if ((c >= 'A') && (c <= 'F')) d = c - 'A' + 10;
    else continue; // do not consider unrecognized characters
    dec += d*base;
    base *= 16;		
  }
  
	
  return dec;
}



/******************************************************************/
//  FUNCTION Definition: xf_StitchElement
static int 
xf_StitchElement(xf_Mesh *Mesh, int egrp, int elem, int *fnodes, 
                 enum xfe_Bool OrientOut, int *face)
{
  int t, nodes[xf_MAXQ1FACENODE];
  
  switch (Mesh->ElemGroup[egrp].QBasis){
    case xfe_TriLagrange:
      //always use face in positive orientation
      for (t=0; t < 2; t++){
        if (OrientOut)
          nodes[t]=fnodes[t];
        else 
          nodes[t]=fnodes[1-t];
      }
      if (Mesh->ElemGroup[egrp].Node[elem][0] == -1){//first visit
        //Bottom face
        Mesh->ElemGroup[egrp].Node[elem][0] = nodes[1];
        Mesh->ElemGroup[egrp].Node[elem][1] = nodes[0];
        (*face) = 2;//following numbering in xf_Q1NodesOnFace
      }
      else {
        //Right or Left face
        if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[1]){
          Mesh->ElemGroup[egrp].Node[elem][2] = nodes[0];
          (*face) = 0;//following numbering in xf_Q1NodesOnFace
        }
        else{ 
          Mesh->ElemGroup[egrp].Node[elem][2] = nodes[1];
          (*face) = 1;//following numbering in xf_Q1NodesOnFace
        }
      }
      break;
    case xfe_QuadLagrange:
      //always use face in positive orientation
      for (t=0; t < 2; t++){
        if (OrientOut)
          nodes[t]=fnodes[t];
        else 
          nodes[t]=fnodes[1-t];
      }
      if (Mesh->ElemGroup[egrp].Node[elem][0] == -1){//first visit
        //South face
        Mesh->ElemGroup[egrp].Node[elem][0] = nodes[1];
        Mesh->ElemGroup[egrp].Node[elem][1] = nodes[0];
        (*face) = 0;//following numbering in xf_Q1NodesOnFace
      }
      else {
        //East face
        if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[1]){
          Mesh->ElemGroup[egrp].Node[elem][3] = nodes[0];
          (*face) = 1;//following numbering in xf_Q1NodesOnFace
        }
        //West face
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[0]){
          Mesh->ElemGroup[egrp].Node[elem][2] = nodes[1];
          (*face) = 3;//following numbering in xf_Q1NodesOnFace
        }
        //North face
        else{
          Mesh->ElemGroup[egrp].Node[elem][2] = nodes[0];
          Mesh->ElemGroup[egrp].Node[elem][3] = nodes[1];
          (*face) = 2;//following numbering in xf_Q1NodesOnFace
        }
      }
      break;
    case xfe_HexLagrange:
      //always use face in positive orientation
      for (t=0; t < 4; t++){
        if (OrientOut)
          nodes[t]=fnodes[t];
        else 
          nodes[t]=fnodes[3-t];
      }
      if (Mesh->ElemGroup[egrp].Node[elem][0] == -1){//first visit
        //Bottom face
        Mesh->ElemGroup[egrp].Node[elem][0] = nodes[0];
        Mesh->ElemGroup[egrp].Node[elem][1] = nodes[3];
        Mesh->ElemGroup[egrp].Node[elem][2] = nodes[1];
        Mesh->ElemGroup[egrp].Node[elem][3] = nodes[2];
        (*face) = 0;//following numbering in xf_Q1NodesOnFace
      }
      else {
        //check for different rotations
        //Front face (rotation 1)
        if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[0] &&
            Mesh->ElemGroup[egrp].Node[elem][1] == nodes[1]){
          Mesh->ElemGroup[egrp].Node[elem][4] = nodes[3];
          Mesh->ElemGroup[egrp].Node[elem][5] = nodes[2];
          (*face) = 1;//following numbering in xf_Q1NodesOnFace
        }
        //Front face (rotation 2)
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[3] &&
                 Mesh->ElemGroup[egrp].Node[elem][1] == nodes[0]){
          Mesh->ElemGroup[egrp].Node[elem][4] = nodes[2];
          Mesh->ElemGroup[egrp].Node[elem][5] = nodes[1];
          (*face) = 1;//following numbering in xf_Q1NodesOnFace
        }
        //Front face (rotation 3)
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[2] &&
                 Mesh->ElemGroup[egrp].Node[elem][1] == nodes[3]){
          Mesh->ElemGroup[egrp].Node[elem][4] = nodes[1];
          Mesh->ElemGroup[egrp].Node[elem][5] = nodes[0];
          (*face) = 1;//following numbering in xf_Q1NodesOnFace
        }
        //Front face (rotation 4)
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[1] &&
                 Mesh->ElemGroup[egrp].Node[elem][1] == nodes[2]){
          Mesh->ElemGroup[egrp].Node[elem][4] = nodes[0];
          Mesh->ElemGroup[egrp].Node[elem][5] = nodes[3];
          (*face) = 1;//following numbering in xf_Q1NodesOnFace
        }
        //check for different rotations
        //Back face (rotation 1)
        else if (Mesh->ElemGroup[egrp].Node[elem][2] == nodes[1] &&
                 Mesh->ElemGroup[egrp].Node[elem][3] == nodes[0]){
          Mesh->ElemGroup[egrp].Node[elem][6] = nodes[2];
          Mesh->ElemGroup[egrp].Node[elem][7] = nodes[3];
          (*face) = 3;//following numbering in xf_Q1NodesOnFace
        }
        //Back face (rotation 2)
        else if (Mesh->ElemGroup[egrp].Node[elem][2] == nodes[0] &&
                 Mesh->ElemGroup[egrp].Node[elem][3] == nodes[3]){
          Mesh->ElemGroup[egrp].Node[elem][6] = nodes[1];
          Mesh->ElemGroup[egrp].Node[elem][7] = nodes[2];
          (*face) = 3;//following numbering in xf_Q1NodesOnFace
        }
        //Back face (rotation 3)
        else if (Mesh->ElemGroup[egrp].Node[elem][2] == nodes[3] &&
                 Mesh->ElemGroup[egrp].Node[elem][3] == nodes[2]){
          Mesh->ElemGroup[egrp].Node[elem][6] = nodes[0];
          Mesh->ElemGroup[egrp].Node[elem][7] = nodes[1];
          (*face) = 3;//following numbering in xf_Q1NodesOnFace
        }
        //Back face (rotation 4)
        else if (Mesh->ElemGroup[egrp].Node[elem][2] == nodes[2] &&
                 Mesh->ElemGroup[egrp].Node[elem][3] == nodes[1]){
          Mesh->ElemGroup[egrp].Node[elem][6] = nodes[3];
          Mesh->ElemGroup[egrp].Node[elem][7] = nodes[0];
          (*face) = 3;//following numbering in xf_Q1NodesOnFace
        }
        //only check for sides and top
        //Right face (rotation 1)
        else if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[0] &&
                 Mesh->ElemGroup[egrp].Node[elem][3] == nodes[1]){
          (*face) = 2;//following numbering in xf_Q1NodesOnFace
        }
        //Right face (rotation 2)
        else if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[3] &&
                 Mesh->ElemGroup[egrp].Node[elem][3] == nodes[0]){
          (*face) = 2;//following numbering in xf_Q1NodesOnFace
        }
        //Right face (rotation 3)
        else if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[2] &&
                 Mesh->ElemGroup[egrp].Node[elem][3] == nodes[3]){
          (*face) = 2;//following numbering in xf_Q1NodesOnFace
        }
        //Right face (rotation 4)
        else if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[1] &&
                 Mesh->ElemGroup[egrp].Node[elem][3] == nodes[2]){
          (*face) = 2;//following numbering in xf_Q1NodesOnFace
        }
        //Left face (rotation 1)
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[0] &&
                 Mesh->ElemGroup[egrp].Node[elem][2] == nodes[3]){
          (*face) = 4;//following numbering in xf_Q1NodesOnFace
        }
        //Left face (rotation 2)
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[3] &&
                 Mesh->ElemGroup[egrp].Node[elem][2] == nodes[2]){
          (*face) = 4;//following numbering in xf_Q1NodesOnFace
        }
        //Left face (rotation 3)
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[2] &&
                 Mesh->ElemGroup[egrp].Node[elem][2] == nodes[1]){
          (*face) = 4;//following numbering in xf_Q1NodesOnFace
        }
        //Left face (rotation 4)
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[1] &&
                 Mesh->ElemGroup[egrp].Node[elem][2] == nodes[0]){
          (*face) = 4;//following numbering in xf_Q1NodesOnFace
        }
        else {//the only option left is the top face
          (*face) = 5;//following numbering in xf_Q1NodesOnFace
        }
      }
      break;
    case xfe_TetLagrange:
      //always use face in positive orientation
      for (t=0; t < 3; t++){
        if (OrientOut)
          nodes[t]=fnodes[t];
        else 
          nodes[t]=fnodes[2-t];
      }
      if (Mesh->ElemGroup[egrp].Node[elem][0] == -1){//first visit
        //Bottom face
        Mesh->ElemGroup[egrp].Node[elem][0] = nodes[2];
        Mesh->ElemGroup[egrp].Node[elem][1] = nodes[1];
        Mesh->ElemGroup[egrp].Node[elem][2] = nodes[0];
        (*face) = 3;//following numbering in xf_Q1NodesOnFace
      }
      else {//check 3 orientations
        if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[0] &&
            Mesh->ElemGroup[egrp].Node[elem][1] == nodes[1]){
          Mesh->ElemGroup[egrp].Node[elem][3] = nodes[2];
          (*face) = 2;//following numbering in xf_Q1NodesOnFace
        }
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[2] &&
                 Mesh->ElemGroup[egrp].Node[elem][1] == nodes[0]){
          Mesh->ElemGroup[egrp].Node[elem][3] = nodes[1];
          (*face) = 2;//following numbering in xf_Q1NodesOnFace
        }
        else if (Mesh->ElemGroup[egrp].Node[elem][0] == nodes[1] &&
                 Mesh->ElemGroup[egrp].Node[elem][1] == nodes[2]){
          Mesh->ElemGroup[egrp].Node[elem][3] = nodes[0];
          (*face) = 2;//following numbering in xf_Q1NodesOnFace
        }
        else if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[0] &&
                 Mesh->ElemGroup[egrp].Node[elem][2] == nodes[1]){
          (*face) = 0;//following numbering in xf_Q1NodesOnFace
        }
        else if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[2] &&
                 Mesh->ElemGroup[egrp].Node[elem][2] == nodes[0]){
          (*face) = 0;//following numbering in xf_Q1NodesOnFace
        }
        else if (Mesh->ElemGroup[egrp].Node[elem][1] == nodes[1] &&
                 Mesh->ElemGroup[egrp].Node[elem][2] == nodes[2]){
          (*face) = 0;//following numbering in xf_Q1NodesOnFace
        }
        else {
          (*face) = 1;//following numbering in xf_Q1NodesOnFace
        }
      }
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
  }
  
  return xf_OK;
}

/******************************************************************/
//  FUNCTION Definition: xf_GetZoneInfo
static int
xf_GetZoneInfo(char *line, char *zone_info)
{
  int i, d;
  
  d = 1;
  while (line[d] != '(') d++;
  
  for (i = 0; i < strlen(line) - d; i++){
    if (line[d+1+i] == ')') break;
    else
      zone_info[i] = line[d+1+i];
  }
  zone_info[i] = '\0';
  
  return xf_OK;
}

/******************************************************************/
//  FUNCTION Definition: xf_AddPeriodicNodePair
static int xf_AddPeriodicNodePair(xf_FluentPeriodicNodePair *NodePair, 
                                  int im_node, int sh_node)
{
  int ierr, dest, src, movesize, rank;
  
  if (NodePair == NULL)
    return xf_Error(xf_INPUT_ERROR);
  
  /* adding entry for the first time? */
  if (NodePair->nP == 0){
    NodePair->nP++;
    if (NodePair->nP > NodePair->nPprealloc){
      NodePair->nPprealloc = NodePair->nP;
      ierr = xf_Error(xf_ReAlloc((void **)&NodePair->im_node, 
                                 NodePair->nPprealloc, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ReAlloc((void **)&NodePair->sh_node, 
                                 NodePair->nPprealloc, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    NodePair->im_node[0] = im_node;
    NodePair->sh_node[0] = sh_node;
  }
  else{
    /* is it a new entry? */
    if (xf_BinSearch(im_node, NodePair->im_node, 0, 
                     NodePair->nP-1, &rank) == xf_NOT_FOUND){
      /* make sure that there is enough space in the 
       arrays for the new entry */
      if (NodePair->nP+1 > NodePair->nPprealloc){
        NodePair->nPprealloc = NodePair->nP+1;
        ierr = xf_Error(xf_ReAlloc((void **)&NodePair->im_node, 
                                   NodePair->nPprealloc, sizeof(int)));
        if (ierr != xf_OK) return ierr;
        
        ierr = xf_Error(xf_ReAlloc((void **)&NodePair->sh_node, 
                                   NodePair->nPprealloc, sizeof(int)));
        if (ierr != xf_OK) return ierr;
      }
      //make sure to keep the crescent order
      if (rank == NodePair->nP){
        movesize = 0;
      }
      else if (rank == -1){
        //move all one position forward
        rank = 0;
        src = 0;
        dest = 1;
        movesize = NodePair->nP;
      }
      else {
        src = rank;
        dest = rank+1;
        movesize = NodePair->nP-rank;
      }
      
      if (movesize > 0){
        if (memmove(NodePair->im_node+dest,NodePair->im_node+src,
                    movesize*sizeof(int)) == NULL)
          return xf_Error(xf_MEMORY_ERROR);
        if (memmove(NodePair->sh_node+dest,NodePair->sh_node+src,
                    movesize*sizeof(double)) == NULL)
          return xf_Error(xf_MEMORY_ERROR);
      }
      NodePair->im_node[rank] = im_node;
      NodePair->sh_node[rank] = sh_node;
      NodePair->nP++;
    }
  }
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadGmshFile
int
xf_ReadFmshFile( const char *InputFile, xf_Mesh *Mesh)
{
  int ierr, i, j, nElemTot, egrp, elem_type, face_type, nFaceTot, nFaceSection;
  int fgrp, n, nface, nInterior, nPeriodic, elem, face, fL, fR, choice, image, shadow;
  int image_face, shadow_face, nPeriodic2Rm, rank, src, dest;
  int nnodeL, nnodeR, fvecR[xf_MAXQ1FACENODE], *Node2RmIdx;
  char line[xf_MAXLINELEN], zone_info[xf_MAXLINELEN], zone_id[xf_MAXSTRLEN];
  char bc_string[xf_MAXSTRLEN];
  char hex0[xf_MAXSTRLEN], hex1[xf_MAXSTRLEN], hex2[xf_MAXSTRLEN], string[xf_MAXSTRLEN];
  enum xfe_Bool dim_set;
  char **info = NULL;
  int egL, egR, eL, eR, im, sh, im_node, sh_node, d, nNode2Rm, *Node2Rm;
  real x[3], average_distance, temp;
  xf_FluentFaceGroup *FaceGroup;
  xf_FluentCellGroup *CellGroup;
  xf_FluentPeriodicNodePair *NodePair = NULL;
  xf_FluentPeriodicFacePair *FacePair = NULL;
  FILE *fid;
  
  // Check input
  if (InputFile == NULL) 
    return xf_Error(xf_INPUT_ERROR);
  
  if (InputFile != NULL){
    // Open InputFile
    fid = fopen(InputFile, "r");
    if(!fid){
      xf_printf("File %s not found.\n\n", InputFile);
      return xf_Error(xf_FILE_READ_ERROR);
    }
    xf_printf("Reading Fluent .msh file.\n");
  }
  
  //Get Dimension
  dim_set = xfe_False;
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_True, "(2",line));
  if (ierr == xf_NOT_FOUND) 
    xf_printf("Warning: no dimension section found in file %s\n",InputFile);
  else if(ierr != xf_OK)
    xf_Error(xf_FILE_READ_ERROR);
  
  if (ierr == xf_OK){
    ierr = sscanf(line,"(2 %d)\n",&Mesh->Dim);
    if (ierr != 1) return xf_Error(xf_FILE_READ_ERROR);
    dim_set = xfe_True;
  }
  
  //Get number of nodes
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_True, "(10",line));
  if (ierr != xf_OK) xf_Error(xf_FILE_READ_ERROR);
  
  ierr = xf_Error(xf_GetZoneInfo(line,zone_info));
  if (ierr != xf_OK) return ierr;
  
  ierr = sscanf(zone_info,"%s %s %s %d %d\n",zone_id, hex0, hex1, &i, &j);
  if (ierr != 5) return xf_Error(xf_FILE_READ_ERROR);
  if (j != Mesh->Dim && dim_set == xfe_True){
    xf_printf("Inconsistent dimension in node section of file %s\n",InputFile);
    return xf_Error(xf_FILE_READ_ERROR);
  }
  else {
    Mesh->Dim = j;
    dim_set = xfe_True;
  }
  //nNode
  Mesh->nNode = xf_Hex2DecInt(hex1);
  ierr = xf_Error(xf_Alloc2((void ***)&Mesh->Coord, Mesh->nNode, Mesh->Dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  if (xf_Hex2DecInt(zone_id) == 0){
    //coordinates
    ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_False, "(10",line));
    if (ierr != xf_OK) xf_Error(xf_FILE_READ_ERROR);
    ierr = xf_Error(xf_GetZoneInfo(line,zone_info));
    if (ierr != xf_OK) return ierr;
    ierr = sscanf(zone_info,"%s %s %s %d %d\n",zone_id, hex0, hex1, &i, &j);
    if (ierr != 5) return xf_Error(xf_FILE_READ_ERROR);
  }
  if (xf_Hex2DecInt(hex0) != 1) {
    xf_printf("Node coordinates section does not start with node 1\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  //finding the node coordinates beginning
  for (i = strlen(line)-4;i < strlen(line);i++){
    if (line[i] == '(')
      break;
  }
  if (line[i] != '(')
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  
  for (i = 0; i < Mesh->nNode; i++){
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    ierr = xf_Error(xf_ScanReal(line, Mesh->Dim, x));
    if (ierr != xf_OK) return ierr;
    for (j = 0; j < Mesh->Dim; j++)
      Mesh->Coord[i][j] = x[j];
  }
  ierr = xf_SetFilePositionAfterString_Serial(fid, xfe_False, "(10",line);
  if (ierr != xf_NOT_FOUND) {
    xf_printf("Multiple node coordinates sections\n");
    return xf_Error(xf_NOT_SUPPORTED);
  }
  
  //Get number of cell groups and cell types
  CellGroup = NULL;
  rewind(fid);
  while (!feof(fid)) {
    ierr = xf_SetFilePositionAfterString_Serial(fid, xfe_False, "(12",line);
    if (ierr == xf_NOT_FOUND)
      break;
    if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR);
    
    ierr = xf_Error(xf_GetZoneInfo(line,zone_info));
    if (ierr != xf_OK) return ierr;
    
    ierr = sscanf(zone_info,"%s %s %s %d %d\n",zone_id, hex0, hex1, &i, &elem_type);
    if (ierr != 5) return xf_Error(xf_FILE_READ_ERROR);
    //declaration
    if (xf_Hex2DecInt(zone_id) == 0) {
      nElemTot = xf_Hex2DecInt(hex1);
      if (nElemTot <= 0){
        xf_printf("No cells in the mesh\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
    }
    else {
      egrp = Mesh->nElemGroup;
      Mesh->nElemGroup++;
      ierr = xf_Error(xf_ReAlloc((void **)&CellGroup, Mesh->nElemGroup, 
                                 sizeof(xf_FluentCellGroup)));
      if (ierr != xf_OK) return ierr;
      CellGroup[egrp].First = xf_Hex2DecInt(hex0)-1;//zero base
      CellGroup[egrp].Last = xf_Hex2DecInt(hex1)-1;//zero base
      ierr = xf_Error(xf_FluentElemType2QBasis(elem_type, 
                                               &CellGroup[egrp].QBasis));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_Basis2nFace(CellGroup[egrp].QBasis, 
                                     &CellGroup[egrp].nFace));
      if (ierr != xf_OK) return ierr;
    }
  }
  if (Mesh->nElemGroup == 0) {
    xf_printf("No cells in the mesh\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  //Get faces
  FaceGroup = NULL;
  rewind(fid);
  nFaceSection = 0;
  nInterior = nPeriodic = 0;
  while (!feof(fid)) {
    ierr = xf_SetFilePositionAfterString_Serial(fid, xfe_False, "(13",line);
    if (ierr == xf_NOT_FOUND)
      break;
    if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR);
    
    ierr = xf_Error(xf_GetZoneInfo(line,zone_info));
    if (ierr != xf_OK) return ierr;
    
    ierr = sscanf(zone_info,"%s %s %s %s %d\n",zone_id, hex0, hex1, hex2, &face_type);
    if (ierr != 5) return xf_Error(xf_FILE_READ_ERROR);
    //declaration
    if (xf_Hex2DecInt(zone_id) == 0) {
      nFaceTot = xf_Hex2DecInt(hex1);
      if (nFaceTot <= 0){
        xf_printf("No faces in the mesh\n");
        return xf_Error(xf_FILE_READ_ERROR);
      }
    }
    else {
      fgrp = nFaceSection;
      nFaceSection++;
      ierr = xf_Error(xf_ReAlloc((void **)&FaceGroup, nFaceSection, 
                                 sizeof(xf_FluentFaceGroup)));
      if (ierr != xf_OK) return ierr;
      
      FaceGroup[fgrp].nFace = xf_Hex2DecInt(hex1)-xf_Hex2DecInt(hex0)+1;
      FaceGroup[fgrp].First = xf_Hex2DecInt(hex0)-1;//zero base
      FaceGroup[fgrp].Last = xf_Hex2DecInt(hex1)-1;//zero base
     
      //~~~~~by Yu for recent version of fluent .msh file
      //if face_type == 0, user should specify by himself
      if(face_type == 0){
      ierr = xf_Error(xf_FluentFaceType2QBasis(2, 
                                               &FaceGroup[fgrp].QBasis));
      if (ierr != xf_OK) return ierr;
      }
      else
      {
      ierr = xf_Error(xf_FluentFaceType2QBasis(face_type, 
                                               &FaceGroup[fgrp].QBasis));
      if (ierr != xf_OK) return ierr;
      }

      FaceGroup[fgrp].ZoneId = xf_Hex2DecInt(zone_id);
      //cr
      ierr = xf_Error(xf_Alloc((void **)&FaceGroup[fgrp].cr, 
                               FaceGroup[fgrp].nFace, 
                               sizeof(int)));
      if (ierr != xf_OK) return ierr;
      //cl
      ierr = xf_Error(xf_Alloc((void **)&FaceGroup[fgrp].cl, 
                               FaceGroup[fgrp].nFace, 
                               sizeof(int)));
      if (ierr != xf_OK) return ierr;
      //nNode
      ierr = xf_Error(xf_Order2nNode(FaceGroup[fgrp].QBasis, 1, 
                                     &FaceGroup[fgrp].nNode));
      if (ierr != xf_OK) return ierr;
      //Node
      ierr = xf_Error(xf_Alloc2((void ***)&FaceGroup[fgrp].Node, 
                                FaceGroup[fgrp].nFace, 
                                FaceGroup[fgrp].nNode, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      //ZoneType
      ierr = xf_Error(xf_FluentBC2ZoneType(xf_Hex2DecInt(hex2), 
                                           &FaceGroup[fgrp].ZoneType));
      if (ierr != xf_OK) return ierr;
      if (FaceGroup[fgrp].ZoneType == xfe_FluentInterior) {
        nInterior++;
      }
      if (FaceGroup[fgrp].ZoneType == xfe_FluentPeriodic) {
        nPeriodic++;
      }
      
      FaceGroup[fgrp].Title = NULL;
      
      //finding the face section beginning
      for (i = strlen(line)-4;i < strlen(line);i++){
        if (line[i] == '(')
          break;
      }
      if (line[i] != '(')
        if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
      if (line[i] != '(')
        if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
      for(i = 0; i < FaceGroup[fgrp].nFace; i++){
        if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
        ierr = xf_Error(xf_ScanXStringAlloc(line, xf_MAXSTRLEN, &n, &info));
        if (ierr != xf_OK) return ierr;
      
        //~~~~~by Yu
        //the first entry is type index
        if(face_type == 0){
        if (n != FaceGroup[fgrp].nNode + 3) //type_index + nodes + cl + cr
          return xf_Error(xf_FILE_READ_ERROR);
        for (j = 0; j < FaceGroup[fgrp].nNode; j++){
          FaceGroup[fgrp].Node[i][j] = xf_Hex2DecInt(info[j+1])-1;//zero-based
        }
        }
        else
        {
        if (n != FaceGroup[fgrp].nNode + 2) //nodes + cl + cr
          return xf_Error(xf_FILE_READ_ERROR);
        for (j = 0; j < FaceGroup[fgrp].nNode; j++){
          FaceGroup[fgrp].Node[i][j] = xf_Hex2DecInt(info[j])-1;//zero-based
        }
        }

        FaceGroup[fgrp].cl[i] = xf_Hex2DecInt(info[n-1])-1;//zero-based
        FaceGroup[fgrp].cr[i] = xf_Hex2DecInt(info[n-2])-1;
        xf_Release2((void **)info);
      }
    }
  }
  if (nInterior == 0) {
    xf_printf("No interior faces groups was found.\n");
    return xf_Error(xf_FILE_READ_ERROR);
  }
  //Get zone names
  for (fgrp = 0; fgrp < nFaceSection; fgrp++){
    sprintf(string,"(45 (%d",FaceGroup[fgrp].ZoneId);
    ierr = xf_SetFilePositionAfterString_Serial(fid, xfe_True, string,line);
    if (ierr == xf_NOT_FOUND){ 
      xf_printf("Warning: missing header for Face Section %d\n",
                FaceGroup[fgrp].ZoneId);
      continue;
    }
    else if(ierr != xf_OK)
      xf_Error(xf_FILE_READ_ERROR);
    ierr = xf_Error(xf_GetZoneInfo(line,zone_info));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc((void **)&FaceGroup[fgrp].Title, 
                             xf_MAXSTRLEN, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    
    sscanf(zone_info, "%d %s %s", &i, bc_string, FaceGroup[fgrp].Title);
    xf_printf("Read Zone %d: %s\n",FaceGroup[fgrp].ZoneId, 
              FaceGroup[fgrp].Title);
  }
  
  //merge nodes info
  nNode2Rm = 0;
  nPeriodic2Rm = 0;
  Node2Rm = Node2RmIdx = NULL;
  //Get Periodic pairs
  if (nPeriodic > 0) {
    //Create space for periodic pairs
    ierr = xf_Error(xf_ReAlloc((void **)&NodePair, nPeriodic, 
                               sizeof(xf_FluentPeriodicNodePair)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc((void **)&FacePair, nPeriodic, 
                               sizeof(xf_FluentPeriodicFacePair)));
    if (ierr != xf_OK) return ierr;
    //scan file to find periodic pair info
    rewind(fid);
    image = shadow = -1;
    i=-1;
    while (!feof(fid)) {
      //periodic pair zone: 18
      ierr = xf_SetFilePositionAfterString_Serial(fid, xfe_False, "(18", line);
      if (ierr == xf_NOT_FOUND){
        if (i == -1){
          xf_printf("No periodic pair information was found!\n");
          return xf_Error(xf_FILE_READ_ERROR);
        }
        break;
      }
      else if (ierr != xf_OK) return xf_Error(ierr);
      
      i++;//pairs counter
      ierr = xf_Error(xf_GetZoneInfo(line, zone_info));
      if (ierr!=xf_OK) return ierr;
      
      ierr = xf_Error(xf_ScanXStringAlloc(zone_info, xf_MAXSTRLEN, &n, &info));
      if (ierr != xf_OK) return ierr;
      if (n != 4) 
        return xf_Error(xf_FILE_READ_ERROR);
      nface = xf_Hex2DecInt(info[1])-xf_Hex2DecInt(info[0])+1;
      image = xf_Hex2DecInt(info[2]);
      shadow = xf_Hex2DecInt(info[3]);
      xf_Release2((void **)info);
      //find face pairs and store nodal periodicity
      sh = im = -1;
      for (n = 0; n < nFaceSection; n++) {
        if (FaceGroup[n].ZoneId == shadow) 
          sh = n;
        if (FaceGroup[n].ZoneId == image) 
          im = n;
      }
      if (sh < 0 || im < 0)
        return xf_Error(xf_MESH_ERROR);
      
      choice = -1;
      //ask user for periodicity type
      xf_printf("Please specify the periodicity type for the zone pair: \"%s\" and \"%s\"\n",
                FaceGroup[im].Title,FaceGroup[sh].Title);
      xf_printf("0 = Translational\n1 = Rotational\n");
      while((scanf("%d",&choice)!= 1) || (choice != 0 && choice != 1)){
        if (choice == 0) {
          NodePair[i].Periodicity = xfe_PeriodicityTranslational;
        }
        else if (choice == 1){
          NodePair[i].Periodicity = xfe_PeriodicityRotational;
        }
        else {
          xf_printf("0 = Translational\n1 = Rotational\n");
          xf_printf("Please type 0 or 1\n");
          continue;
        }
      }
      
      //allocate space for node pairs
      average_distance = 0.0;
      NodePair[i].nP = 0;
      NodePair[i].nPprealloc = FaceGroup[im].nFace*FaceGroup[im].nNode;
      ierr = xf_Error(xf_Alloc((void **)&NodePair[i].im_node, 
                               NodePair[i].nPprealloc, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      ierr = xf_Error(xf_Alloc((void **)&NodePair[i].sh_node, 
                               NodePair[i].nPprealloc, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      
      FacePair[i].nP = FaceGroup[im].nFace;
      ierr = xf_Error(xf_Alloc((void **)&FacePair[i].sh_fgrp, 
                               FacePair[i].nP, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      ierr = xf_Error(xf_Alloc((void **)&FacePair[i].sh_face, 
                               FacePair[i].nP, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      ierr = xf_Error(xf_Alloc((void **)&FacePair[i].im_fgrp, 
                               FacePair[i].nP, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      ierr = xf_Error(xf_Alloc((void **)&FacePair[i].im_face, 
                               FacePair[i].nP, sizeof(int)));
      if (ierr!=xf_OK) return ierr;
      
      //find beginning of section
      for (n = strlen(line)-4;n < strlen(line);n++){
        if (line[n] == '(')
          break;
      }
      if (line[n] != '(')
        if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
      //loop through faces and pair nodes
      for (face = 0; face < nface; face++) {
        if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
        ierr = xf_Error(xf_ScanXStringAlloc(line, xf_MAXSTRLEN, &n, &info));
        if (ierr != xf_OK) return ierr;
        if (n != 2)//expecting image and shadow 
          return xf_Error(xf_FILE_READ_ERROR);
        image_face = xf_Hex2DecInt(info[0])-1-FaceGroup[im].First;
        shadow_face = xf_Hex2DecInt(info[1])-1-FaceGroup[sh].First;
        xf_Release2((void **)info);
        FacePair[i].im_fgrp[face] = im;
        FacePair[i].im_face[face] = image_face;
        FacePair[i].sh_fgrp[face] = sh;
        FacePair[i].sh_face[face] = shadow_face;
        for (j = 0; j < FaceGroup[im].nNode; j++) {
          im_node = FaceGroup[im].Node[image_face][j];
          sh_node = FaceGroup[sh].Node[shadow_face][FaceGroup[sh].nNode-j-1];
          temp = 0.0;
          for (d = 0; d < Mesh->Dim; d++) {
            temp += pow(Mesh->Coord[im_node][d]-Mesh->Coord[sh_node][d],2.0);
          }
          average_distance += sqrt(temp);
          ierr = xf_Error(xf_AddPeriodicNodePair(NodePair + i, im_node, sh_node));
          if (ierr != xf_OK) return ierr;
        }
      }
      average_distance /= (nface*FaceGroup[im].nNode);
      NodePair[i].average_distance = average_distance;
      xf_printf("average distance between periodic nodes in %s and %s: %1.15E\n",
                average_distance,FaceGroup[im].Title,FaceGroup[sh].Title);
      if (NodePair[i].average_distance < MEPS){
        NodePair[i].MergeNodes = xfe_True;
        nPeriodic2Rm++;
        //mark shadow nodes to remove
        for (j = 0; j < NodePair[i].nP; j++){
          ierr = xf_Error(xf_Add2OrderedSet(NodePair[i].sh_node[j], &nNode2Rm, &Node2Rm, 
                                            &Node2RmIdx));
          if (ierr != xf_OK) return ierr;
        }
      }
      else
        NodePair[i].MergeNodes = xfe_False;
    }//feof(fid)
  }//nPeriodic
  
  if (nPeriodic2Rm > 1)
    return xf_Error(xf_NOT_SUPPORTED);
  //~~~~~~by Yu
  //if (nPeriodic > 1)
  //  return xf_Error(xf_NOT_SUPPORTED);
  
  //allocate element groups
  ierr = xf_Error(xf_ReAlloc((void **) &Mesh->ElemGroup, Mesh->nElemGroup, 
                             sizeof(xf_ElemGroup)));
  if (ierr!=xf_OK) return ierr;
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
    //initialize element group (takes care of CutFlag and CutElemData)
    xf_InitElemGroup(Mesh->ElemGroup + egrp);
    //QOrder
    Mesh->ElemGroup[egrp].QOrder = 1;
    //QBasis
    Mesh->ElemGroup[egrp].QBasis = CellGroup[egrp].QBasis;
    //nElem
    Mesh->ElemGroup[egrp].nElem = CellGroup[egrp].Last-CellGroup[egrp].First+1;
    //Allocate nFace
    ierr = xf_Error(xf_Alloc((void **) &Mesh->ElemGroup[egrp].nFace, 
                             Mesh->ElemGroup[egrp].nElem, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    //nNode
    ierr = xf_Error(xf_Order2nNode(Mesh->ElemGroup[egrp].QBasis, 
                                   Mesh->ElemGroup[egrp].QOrder, 
                                   &Mesh->ElemGroup[egrp].nNode));
    if (ierr != xf_OK) return ierr;
    //Node
    ierr = xf_Error(xf_Alloc2((void ***) &Mesh->ElemGroup[egrp].Node, 
                              Mesh->ElemGroup[egrp].nElem, 
                              Mesh->ElemGroup[egrp].nNode, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    //nFace and initialize nodes
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      Mesh->ElemGroup[egrp].nFace[elem] = CellGroup[egrp].nFace;
      for (i = 0; i < Mesh->ElemGroup[egrp].nNode; i++)
        Mesh->ElemGroup[egrp].Node[elem][i] = -1;//this facilitates stitching later
    }
    //Face
    ierr = xf_Error(xf_VAlloc2((void ***) &Mesh->ElemGroup[egrp].Face, 
                               Mesh->ElemGroup[egrp].nElem, 
                               Mesh->ElemGroup[egrp].nFace, sizeof(xf_Face)));
    if (ierr!=xf_OK) return ierr;
  }
 
  printf("here now 1~\n");getchar();
  //Loop through face groups and stitch elements
  Mesh->nIFace = 0;
  for (fgrp = 0; fgrp < nFaceSection; fgrp++){
    switch (FaceGroup[fgrp].ZoneType) {
      case xfe_FluentInterior:
        Mesh->nIFace += FaceGroup[fgrp].nFace;
        ierr = xf_Error(xf_ReAlloc((void **) &Mesh->IFace, Mesh->nIFace, 
                                   sizeof(xf_IFace)));
        if (ierr!=xf_OK) return ierr;
        for (face = 0; face < FaceGroup[fgrp].nFace; face++){
          //find element group
          for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
            if (FaceGroup[fgrp].cl[face] >= CellGroup[egrp].First &&
                FaceGroup[fgrp].cl[face] <= CellGroup[egrp].Last) {
              egL = egrp;
            }
            if (FaceGroup[fgrp].cr[face] >= CellGroup[egrp].First &&
                FaceGroup[fgrp].cr[face] <= CellGroup[egrp].Last) {
              egR = egrp;
            }
          }
          eL = FaceGroup[fgrp].cl[face]-CellGroup[egL].First;
          eR = FaceGroup[fgrp].cr[face]-CellGroup[egR].First;
          //stitch elements
          //left cell
          ierr = xf_Error(xf_StitchElement(Mesh, egL, eL, 
                                           FaceGroup[fgrp].Node[face], 
                                           xfe_True, &fL));
          if (ierr!=xf_OK) return ierr;
          //right cell
          ierr = xf_Error(xf_StitchElement(Mesh, egR, eR, 
                                           FaceGroup[fgrp].Node[face], 
                                           xfe_False, &fR));
          if (ierr!=xf_OK) return ierr;
          n = face+Mesh->nIFace-FaceGroup[fgrp].nFace;
          Mesh->IFace[n].ElemGroupL = egL;
          Mesh->IFace[n].ElemL = eL;
          Mesh->IFace[n].FaceL = fL;
          Mesh->IFace[n].OrientL = 1;
          
          Mesh->IFace[n].ElemGroupR = egR;
          Mesh->IFace[n].ElemR = eR;
          Mesh->IFace[n].FaceR = fR;          
          Mesh->IFace[n].OrientR = -1;
          
          Mesh->IFace[n].HangNumber = 0;
          Mesh->IFace[n].CutFaceData = NULL;
        }
        break;
      case xfe_FluentOther:
        //considered boundary
        ierr = xf_Error(xf_ReAlloc((void **) &Mesh->BFaceGroup, Mesh->nBFaceGroup+1, 
                                   sizeof(xf_BFaceGroup)));
        if (ierr!=xf_OK) return ierr;
        Mesh->BFaceGroup[Mesh->nBFaceGroup].nBFace = FaceGroup[fgrp].nFace;
        ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup[Mesh->nBFaceGroup].BFace, 
                                 Mesh->BFaceGroup[Mesh->nBFaceGroup].nBFace, 
                                 sizeof(xf_BFace)));
        if (ierr!=xf_OK) return ierr;
        ierr = xf_Error(xf_AllocString(&Mesh->BFaceGroup[Mesh->nBFaceGroup].Title, 
                                       xf_MAXSTRLEN, FaceGroup[fgrp].Title));
        if (ierr != xf_OK) return ierr;
        for (face = 0; face < FaceGroup[fgrp].nFace; face++) {
          //find element group
          for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
            if (FaceGroup[fgrp].cr[face] >= CellGroup[egrp].First &&
                FaceGroup[fgrp].cr[face] <= CellGroup[egrp].Last) {
              egR = egrp;
            }
          }
          eR = FaceGroup[fgrp].cr[face]-CellGroup[egR].First;
          //right cell
          ierr = xf_Error(xf_StitchElement(Mesh, egR, eR, 
                                           FaceGroup[fgrp].Node[face], 
                                           xfe_False, &fR));
          if (ierr!=xf_OK) return ierr;
          
          Mesh->BFaceGroup[Mesh->nBFaceGroup].BFace[face].ElemGroup = egR;
          Mesh->BFaceGroup[Mesh->nBFaceGroup].BFace[face].Elem = eR;
          Mesh->BFaceGroup[Mesh->nBFaceGroup].BFace[face].Face = fR;          
          Mesh->BFaceGroup[Mesh->nBFaceGroup].BFace[face].Orient = -1;
          
          Mesh->BFaceGroup[Mesh->nBFaceGroup].BFace[face].CutFaceData = NULL;
        }
        Mesh->nBFaceGroup++;
        break;
      default:
        break;
    }
  }
  
  printf("here now 2~\n");getchar();
  /* convert periodic faces into interior faces and 
   write node periodicity to memory */
  if (nPeriodic > 0){
    Mesh->nPeriodicGroup = nPeriodic;
    ierr = xf_Error(xf_Alloc((void **) &Mesh->PeriodicGroup, 
                             Mesh->nPeriodicGroup, sizeof(xf_PeriodicGroup)));
    if (ierr!=xf_OK) return ierr;
  }

  for (i = 0; i < nPeriodic; i++){
    Mesh->nIFace += FacePair[i].nP;
    ierr = xf_Error(xf_ReAlloc((void **) &Mesh->IFace, Mesh->nIFace, 
                               sizeof(xf_IFace)));
    if (ierr!=xf_OK) return ierr;
    for (j = 0; j < FacePair[i].nP; j++) {
      image = FacePair[i].im_fgrp[j];
      image_face = FacePair[i].im_face[j];
      shadow = FacePair[i].sh_fgrp[j];
      shadow_face = FacePair[i].sh_face[j];
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
        if (FaceGroup[image].cr[image_face] >= CellGroup[egrp].First &&
            FaceGroup[image].cr[image_face] <= CellGroup[egrp].Last) {
          egR = egrp;
        }
        if (FaceGroup[shadow].cr[shadow_face] >= CellGroup[egrp].First &&
            FaceGroup[shadow].cr[shadow_face] <= CellGroup[egrp].Last) {
          egL = egrp;
        }
      }
      eL = FaceGroup[shadow].cr[shadow_face]-CellGroup[egL].First;
      eR = FaceGroup[image].cr[image_face]-CellGroup[egR].First;
      //stitch elements      
      //left cell
      ierr = xf_Error(xf_StitchElement(Mesh, egL, eL, 
                                       FaceGroup[shadow].Node[shadow_face], 
                                       xfe_False, &fL));
      if (ierr!=xf_OK) return ierr;
      //right cell
      ierr = xf_Error(xf_StitchElement(Mesh, egR, eR, 
                                       FaceGroup[image].Node[image_face], 
                                       xfe_False, &fR));
      if (ierr!=xf_OK) return ierr;
      n = j+Mesh->nIFace-FacePair[i].nP;
      Mesh->IFace[n].ElemGroupL = egL;
      Mesh->IFace[n].ElemL = eL;
      Mesh->IFace[n].FaceL = fL;
      Mesh->IFace[n].OrientL = 1;
      
      Mesh->IFace[n].ElemGroupR = egR;
      Mesh->IFace[n].ElemR = eR;
      Mesh->IFace[n].FaceR = fR;          
      Mesh->IFace[n].OrientR = -1;
      
      Mesh->IFace[n].HangNumber = 0;
      Mesh->IFace[n].CutFaceData = NULL;
    }
    //check if we are merging the nodes
    Mesh->PeriodicGroup[i].Periodicity = NodePair[i].Periodicity;
    Mesh->PeriodicGroup[i].nPeriodicNode = NodePair[i].nP;
    ierr = xf_Error(xf_Alloc((void **) &Mesh->PeriodicGroup[i].PeriodicNode, 
                             Mesh->PeriodicGroup[i].nPeriodicNode*2, 
                             sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    for (j = 0; j < NodePair[i].nP; j++){
      Mesh->PeriodicGroup[i].PeriodicNode[2*j] = NodePair[i].im_node[j];
      Mesh->PeriodicGroup[i].PeriodicNode[2*j+1] = NodePair[i].sh_node[j];
    }
  }
  
  //correct node numbers if we merged nodes
  if (nPeriodic2Rm) {
    xf_printf("Merging nodes...");fflush(stdout);
    for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
      for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++) {
        for (n = 0; n < Mesh->ElemGroup[egrp].nNode; n++) {
          /*loop through nodes and find their rank with respect to 
          the nodes will be removed*/
          ierr = xf_BinSearch(Mesh->ElemGroup[egrp].Node[elem][n], 
                              Node2Rm, 0, nNode2Rm-1, &rank);
          if (ierr == xf_NOT_FOUND){
            //do not remove it, just renumber it
            //if rank == -1 no shift is necessary
            if (rank != -1 )
              if (rank <= nNode2Rm)
                Mesh->ElemGroup[egrp].Node[elem][n]-=rank;
              else return xf_Error(xf_CODE_LOGIC_ERROR);
          }
          else if (ierr == xf_OK){
            //replace it by its image-node
            im_node = NodePair[0].im_node[Node2RmIdx[rank]];
            ierr = xf_BinSearch(im_node, Node2Rm, 0, nNode2Rm-1, &rank);
            if (ierr == xf_OK) return xf_Error(xf_CODE_LOGIC_ERROR);//it should not be found
            if (ierr != xf_NOT_FOUND) return xf_Error(ierr);//it should not be found
            if (rank <= nNode2Rm && rank != -1)
              im_node-=rank;
            Mesh->ElemGroup[egrp].Node[elem][n]=im_node;
          }
          else return xf_Error(ierr);
        }
      }
    }
    //shift and resize the coordinate array
    for (i=0; i<nNode2Rm; i++) {
      src = Node2Rm[i]-i+1;
      dest = Node2Rm[i]-i;
      if (memmove((*Mesh->Coord)+dest*Mesh->Dim,
                  (*Mesh->Coord)+src*Mesh->Dim,
                  Mesh->Dim*(Mesh->nNode-src)*sizeof(real))==NULL)
        return xf_Error(xf_MEMORY_ERROR);
    }
    Mesh->nNode-=nNode2Rm;
    ierr = xf_Error(xf_ReAlloc((void **)&Mesh->Coord, 
                                Mesh->nNode*Mesh->Dim, 
                                sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    Mesh->nPeriodicGroup = 0;
    ierr = xf_Error(xf_DestroyPeriodicGroup(Mesh->PeriodicGroup+0));
    if (ierr != xf_OK) return ierr;
    
    xf_printf("done.\n\n");fflush(stdout);
  }
  
  //clean up
  xf_Release((void *)CellGroup);
  for (i = 0; i < nFaceSection; i++) {
    xf_Release2((void **)FaceGroup[i].Node);
    xf_Release((void *)FaceGroup[i].cr);
    xf_Release((void *)FaceGroup[i].cl);
    xf_Release((void *)FaceGroup[i].Title);
  }
  xf_Release((void *)FaceGroup);
  for (i = 0; i < nPeriodic; i++) {
    xf_Release((void *)NodePair[i].im_node);
    xf_Release((void *)NodePair[i].sh_node);
    xf_Release((void *)FacePair[i].im_fgrp);
    xf_Release((void *)FacePair[i].im_face);
    xf_Release((void *)FacePair[i].sh_fgrp);
    xf_Release((void *)FacePair[i].sh_face);
  }
  xf_Release((void *)NodePair);
  xf_Release((void *)FacePair);
  xf_Release((void *)Node2Rm);
  xf_Release((void *)Node2RmIdx);
  
  return xf_OK;
}