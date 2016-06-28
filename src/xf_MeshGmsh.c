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
 FILE:  xf_MeshGmsh.c
 
 This file contains functions for working with the Gmsh format.
 
 */

#include "xf_MeshStruct.h"

static const int NTYPES = 93;
static const char GMSHVERSION[] = "2.2";

typedef struct
{
  int nNode;
  /* Number of nodes */
  
  enum xfe_BasisType QBasis;
  /* Geometrical basis */
  
  int QOrder;
  /* Geometrical basis order */
  
  enum xfe_ShapeType Shape; 
  /* Geometrical shape */
  
  enum xfe_Bool Supported;
  /* True is element-type is supported in Xflow */
  
  int *NodeOrder;
  /* Index to array is XFlow numbering, and content of list is Gmsh numbering */
}
xf_GmshElem; //Please refer to http://www.geuz.org/gmsh for information

typedef struct
{
  int dim;
  int Group;//xflow-based group numbering (-1 corresponds to unused group or element group)
  int tag;  //tag to link physical name with grouped element; extend by Yu
  char *Title;
}
xf_GmshPEntity;


/******************************************************************/
//   FUNCTION Definition: xf_CreateGmshElemInfo
static int
xf_CreateGmshElemInfo(xf_GmshElem **pelem_type)
{
  int ierr, i;
  const int ntypes = NTYPES;
  xf_GmshElem *elem_type;
  
  /* Element order - DO NOT CHANGE */
  static int elem_order1[] = {0, 1};
  static int elem_order2[] = {0, 1, 2};
  static int elem_order3[] = {0, 1, 3, 2};
  static int elem_order4[] = {0, 1, 2, 3};
  static int elem_order5[] = {0, 1, 3, 2, 4, 5, 7, 6};
  static int elem_order8[] = {0, 2, 1};
  static int elem_order9[] = {0, 3, 1, 5, 4, 2};
  static int elem_order10[] = {0, 4, 1, 7, 8, 5, 3, 6, 2};
  static int elem_order11[] = {0, 4, 1, 6, 5, 2, 7, 9, 8, 3};
  //added by Yu; static int elem_order12[] = {-1};
  static int elem_order12[] = {0, 8, 1, 9, 20, 11, 3, 13, 2, 10, 21, 12,
                               22, 26, 23, 15, 24, 14, 4, 16, 5, 17, 25,
                               18, 7, 19, 6};
  static int elem_order21[] = {0, 3, 4, 1, 8, 9, 5, 7, 6, 2};
  static int elem_order23[] = {0, 3, 4, 5, 1, 11, 12, 13, 6, 
    10, 14, 7, 9, 8, 2};
  static int elem_order26[] = {0, 2, 3, 1};
  static int elem_order27[] = {0, 2, 3, 4, 1};
  static int elem_order28[] = {0, 2, 3, 4, 5, 1};
  static int elem_order36[] = {0, 4, 5, 1, 11, 12, 13, 6, 10, 15, 14, 
    7, 3, 9, 8, 2};
  static int elem_order37[] = {0, 4, 5, 6, 1, 15, 16, 20, 17, 7,
    14, 23, 24, 21, 8, 13, 19, 22, 18, 9,
    3, 12, 11, 10, 2};
  static int elem_order93[] = {0, 8, 9, 10, 1, 11, 44, 51, 47, 17,
    12, 48, 52, 50, 18, 13, 45, 49, 46, 19,
    3, 25, 24, 23, 2,
    14, 53, 57, 54, 20, 62, 98, 106, 99, 71,
    69, 107, 118, 109, 75, 65, 101, 111, 100, 72,
    29, 80, 87, 83, 26,
    15, 60, 61, 58, 21, 66, 108, 119, 110, 78,
    70, 120, 124, 121, 79, 68, 113, 122, 112, 76,
    30, 84, 88, 86, 27,
    16, 56, 59, 55, 22, 63, 102, 114, 103, 74,
    67, 115, 123, 116, 77, 64, 105, 117, 104, 73,
    31, 81, 85, 82, 28,
    4, 32, 33, 34, 5, 35, 89, 93, 90, 38,
    36, 96, 97, 94, 39, 37, 92, 95, 91, 40,
    7, 43, 42, 41, 6};
  
  ierr = xf_Error(xf_Alloc((void **) pelem_type, ntypes + 1, sizeof(xf_GmshElem)));
  if (ierr != xf_OK) return ierr;
  elem_type = (*pelem_type);
  
  /* Default is to not support the element type */
  for (i = 1; i <= ntypes; i++) {
    elem_type[i].nNode = -1;
    elem_type[i].QBasis = xfe_BasisLast;
    elem_type[i].QOrder = -1;
    elem_type[i].Shape = xfe_ShapeLast;
    elem_type[i].Supported = xfe_False;
    elem_type[i].NodeOrder = NULL;
  }
  
  /* Fill the element information */
  elem_type[1].nNode = 2; 
  elem_type[1].QOrder = 1; 
  elem_type[1].QBasis = xfe_SegLagrange; 
  elem_type[1].Shape = xfe_Segment; 
  elem_type[1].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&elem_type[1].NodeOrder, elem_type[1].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[1].nNode; i++)
    elem_type[1].NodeOrder[i] = elem_order1[i];
  
  elem_type[2].nNode = 3;
  elem_type[2].QOrder = 1;
  elem_type[2].QBasis = xfe_TriLagrange;
  elem_type[2].Shape = xfe_Triangle; 
  elem_type[2].Supported = xfe_True; 
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[2].NodeOrder), elem_type[2].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[2].nNode; i++)
    elem_type[2].NodeOrder[i] = elem_order2[i];
  
  elem_type[3].nNode = 4;
  elem_type[3].QOrder = 1;
  elem_type[3].QBasis = xfe_QuadLagrange;
  elem_type[3].Shape = xfe_Quadrilateral;
  elem_type[3].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[3].NodeOrder), elem_type[3].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[3].nNode; i++)
    elem_type[3].NodeOrder[i] = elem_order3[i];
  
  elem_type[4].nNode = 4; 
  elem_type[4].QBasis = xfe_TetLagrange; 
  elem_type[4].Shape = xfe_Tetrahedron; 
  elem_type[4].QOrder = 1; 
  elem_type[4].Supported = xfe_True; 
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[4].NodeOrder), elem_type[4].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[4].nNode; i++)
    elem_type[4].NodeOrder[i] = elem_order4[i];
  
  elem_type[5].nNode = 8; 
  elem_type[5].QBasis = xfe_HexLagrange; 
  elem_type[5].QOrder = 1; 
  elem_type[5].Shape = xfe_Hexahedron; 
  elem_type[5].Supported = xfe_True; 
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[5].NodeOrder), elem_type[5].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[5].nNode; i++)
    elem_type[5].NodeOrder[i] = elem_order5[i];
  
  
  /* 6-node prism (not supported) */
  elem_type[6].nNode = 6;
  
  /* 5-node pyramid (not supported) */
  elem_type[7].nNode = 5;
  
  elem_type[8].nNode = 3;
  elem_type[8].QBasis = xfe_SegLagrange;
  elem_type[8].QOrder = 2;
  elem_type[8].Shape = xfe_Segment;
  elem_type[8].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[8].NodeOrder), elem_type[8].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[8].nNode; i++)
    elem_type[8].NodeOrder[i] = elem_order8[i];
  
  elem_type[9].nNode = 6;
  elem_type[9].QOrder = 2;
  elem_type[9].QBasis = xfe_TriLagrange;
  elem_type[9].Shape = xfe_Triangle;
  elem_type[9].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[9].NodeOrder), elem_type[9].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[9].nNode; i++)
    elem_type[9].NodeOrder[i] = elem_order9[i];
  
  elem_type[10].nNode = 9; 
  elem_type[10].QOrder = 2;
  elem_type[10].QBasis = xfe_QuadLagrange;
  elem_type[10].Shape = xfe_Quadrilateral;
  elem_type[10].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[10].NodeOrder), elem_type[10].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[10].nNode; i++)
    elem_type[10].NodeOrder[i] = elem_order10[i];
  
  elem_type[11].nNode = 10; 
  elem_type[11].QOrder = 2;
  elem_type[11].QBasis = xfe_TetLagrange;
  elem_type[11].Shape = xfe_Tetrahedron;
  elem_type[11].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[11].NodeOrder), elem_type[11].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[11].nNode; i++)
    elem_type[11].NodeOrder[i] = elem_order11[i];
  
  elem_type[12].nNode = 27; 
  elem_type[12].QOrder = 2;
  elem_type[12].QBasis = xfe_HexLagrange;
  elem_type[12].Shape = xfe_Hexahedron;
  elem_type[12].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[12].NodeOrder), elem_type[12].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[12].nNode; i++)
    elem_type[12].NodeOrder[i] = elem_order12[i];
  
  /* 18-node q=2 prism (not supported) */
  elem_type[13].nNode = 18;
  
  /* 14-node q=2 prism (not supported) */
  elem_type[14].nNode = 14;
  
  /* point (q=0) */
  elem_type[15].nNode = 1;
  elem_type[15].QOrder =  0;
  elem_type[15].Shape = xfe_ShapeLast;
  elem_type[15].QBasis = xfe_BasisLast;
  elem_type[15].Supported = xfe_True;
  
  elem_type[16].nNode = 8; 
  elem_type[16].QOrder = 2;
  elem_type[16].Shape = xfe_Quadrilateral;
  elem_type[16].QBasis = xfe_QuadLagrange;
  elem_type[16].Supported = xfe_False;
  
  elem_type[17].nNode = 20; 
  elem_type[17].QOrder = 2;
  elem_type[17].Shape = xfe_Hexahedron;
  elem_type[17].QBasis = xfe_HexLagrange;
  elem_type[17].Supported = xfe_False;
  
  elem_type[18].nNode = 15; 
  elem_type[18].QOrder = -1;
  elem_type[18].Shape = xfe_ShapeLast;
  elem_type[18].QBasis = xfe_BasisLast;
  elem_type[18].Supported = xfe_False;
  
  elem_type[19].nNode = 13; 
  elem_type[19].QOrder = -1;
  elem_type[19].Shape = xfe_ShapeLast;
  elem_type[19].QBasis = xfe_BasisLast;
  elem_type[19].Supported = xfe_False;
  
  elem_type[20].nNode = 9; 
  elem_type[20].QOrder = 3;
  elem_type[20].Shape = xfe_Triangle;
  elem_type[20].QBasis = xfe_TriLagrange;
  elem_type[20].Supported = xfe_False;
  
  
  elem_type[21].nNode = 10; 
  elem_type[21].QOrder = 3;
  elem_type[21].Shape = xfe_Triangle;
  elem_type[21].QBasis = xfe_TriLagrange;
  elem_type[21].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[21].NodeOrder), elem_type[21].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[21].nNode; i++)
    elem_type[21].NodeOrder[i] = elem_order21[i];
  
  
  /* 12-node q=4 incomplete triangle (not supported) */
  elem_type[22].nNode = 12; 
  elem_type[22].QOrder = 4;
  elem_type[22].QBasis = xfe_TriLagrange;
  elem_type[22].Shape = xfe_Triangle;
  elem_type[22].Supported = xfe_False;
  
  elem_type[23].nNode = 15; 
  elem_type[23].QOrder = 4;
  elem_type[23].QBasis = xfe_TriLagrange;
  elem_type[23].Shape = xfe_Triangle;
  elem_type[23].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[23].NodeOrder), elem_type[23].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[23].nNode; i++)
    elem_type[23].NodeOrder[i] = elem_order23[i];
  
  /* 15-node q=5 incomplete triangle (not supported) */
  elem_type[24].nNode = 15; 
  elem_type[24].QOrder = 5;
  elem_type[24].QBasis = xfe_TriLagrange;
  elem_type[24].Shape = xfe_Triangle;
  elem_type[24].Supported = xfe_False;
  
  elem_type[25].nNode = 21; 
  elem_type[25].QOrder = 5;
  elem_type[25].QBasis = xfe_TriLagrange;
  elem_type[25].Shape = xfe_Triangle;
  elem_type[25].Supported = xfe_True;
  /* need node list */
  
  elem_type[26].nNode = 4; 
  elem_type[26].QOrder = 3;
  elem_type[26].QBasis = xfe_SegLagrange;
  elem_type[26].Shape = xfe_Segment;
  elem_type[26].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[26].NodeOrder), elem_type[26].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[26].nNode; i++)
    elem_type[26].NodeOrder[i] = elem_order26[i];
  
  elem_type[27].nNode = 5; 
  elem_type[27].QOrder = 4;
  elem_type[27].QBasis = xfe_SegLagrange;
  elem_type[27].Shape = xfe_Segment;
  elem_type[27].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[27].NodeOrder), elem_type[27].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[27].nNode; i++)
    elem_type[27].NodeOrder[i] = elem_order27[i];
  
  elem_type[28].nNode = 6; 
  elem_type[28].QOrder = 5;
  elem_type[28].QBasis = xfe_SegLagrange;
  elem_type[28].Shape = xfe_Segment;
  elem_type[28].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[28].NodeOrder), elem_type[28].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[28].nNode; i++)
    elem_type[28].NodeOrder[i] = elem_order28[i];
  
  elem_type[29].nNode = 20; 
  elem_type[29].QOrder = 3;
  elem_type[29].QBasis = xfe_TetLagrange;
  elem_type[29].Shape = xfe_Tetrahedron;
  elem_type[29].Supported = xfe_True;
  /* need node list */
  
  elem_type[30].nNode = 35; 
  elem_type[30].QOrder = 4;
  elem_type[30].QBasis = xfe_TetLagrange;
  elem_type[30].Shape = xfe_Tetrahedron;
  elem_type[30].Supported = xfe_True;
  /* need node list */
  
  elem_type[31].nNode = 56;
  elem_type[31].QOrder = 5;
  elem_type[31].QBasis = xfe_TetLagrange;
  elem_type[31].Shape = xfe_Tetrahedron;
  elem_type[31].Supported = xfe_True;
  /* need node list */
 
  //new added May 2014
  elem_type[36].nNode = 16;
  elem_type[36].QOrder = 3;
  elem_type[36].QBasis = xfe_QuadLagrange;
  elem_type[36].Shape = xfe_Quadrilateral;
  elem_type[36].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[36].NodeOrder), elem_type[36].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[36].nNode; i++)
    elem_type[36].NodeOrder[i] = elem_order36[i];

  elem_type[37].nNode = 25;
  elem_type[37].QOrder = 4;
  elem_type[37].QBasis = xfe_QuadLagrange;
  elem_type[37].Shape = xfe_Quadrilateral;
  elem_type[37].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[37].NodeOrder), elem_type[37].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[37].nNode; i++)
    elem_type[37].NodeOrder[i] = elem_order37[i];
  
  elem_type[93].nNode = 125;
  elem_type[93].QOrder = 4;
  elem_type[93].QBasis = xfe_HexLagrange;
  elem_type[93].Shape = xfe_Hexahedron;
  elem_type[93].Supported = xfe_True;
  ierr = xf_Error(xf_Alloc((void **)&(elem_type[93].NodeOrder), elem_type[93].nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i < elem_type[93].nNode; i++)
    elem_type[93].NodeOrder[i] = elem_order93[i];
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyGmshElemInfo
static int
xf_DestroyGmshElemInfo(xf_GmshElem *elem_type)
{
  /*
   PURPOSE: Destroys gmsh elem_type lookup structure
   INPUTS: elem_type = gmsh element type structure
   OUTPUS: none
   RETURN: error code
   */
  int i, ntypes = NTYPES;
  
  /* Default is to not support the element type */
  for (i = 1; i <= ntypes; i++) 
    if (elem_type[i].NodeOrder != NULL) xf_Release( (void *) elem_type[i].NodeOrder);
  
  xf_Release( (void *) elem_type);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_NodalOrderGmsh2Gri
static int
xf_NodalOrderGmsh2Gri(xf_GmshElem elem_type, int *nvec)
{
  int i, ierr;
  int nnode, *buf, index;
  
  if (!elem_type.Supported)
    return xf_NOT_SUPPORTED;
  else if (elem_type.NodeOrder == NULL){
    xf_printf("Element type is supported, but the order does not yet exist\n");
    return xf_NOT_SUPPORTED;
  }
  
  ierr = xf_Error(xf_Order2nNode(elem_type.QBasis, elem_type.QOrder, 
                                 &nnode));
  if (ierr != xf_OK) return ierr;
 
  ierr = xf_Error(xf_Alloc((void **)&buf, nnode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  /* Put order into buffer */
  for (i = 0; i < nnode; i++) buf[i] = nvec[elem_type.NodeOrder[i]];
  
  /* Copy back to array */
  for (i = 0; i < nnode; i++) nvec[i] = buf[i];
  
  xf_Release((void *)buf);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ReadGmshFile
int
xf_ReadGmshFile( const char *InputFile, xf_Mesh *Mesh)
{
  int ierr, nNodes, ftype, nbytes, n, i, j, k, elm_number;
  int elm_type, ntags, np_entity, nfnode, *nvec, *edata;
  int dim, p_entity, nItems, fgrp, bface, egrp, elem; 
  int nnode, nmax, nface, *fvec, face, nElemTot, nLocFaceMax, nIFaceMax;
  char line[xf_MAXLINELEN], ver[xf_MAXSTRLEN], buffer[xf_MAXSTRLEN];
  enum xfe_Bool found, exists, need_interpolation;
  real **nodes, dmin[3], dmax[3], TV[3], x,y,z;
  xf_FaceHash **Node2Hash, *hinfo;
  xf_GmshElem *GmshTypesInfo;
  xf_GmshPEntity  *physical_entity;
  xf_ElemGroup *DummyEGroup;
  xf_BFaceGroup *DummyBFGroup;
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
    xf_printf("Reading .msh file.\n");
  }
  
  //Get MeshFormat
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_True, "$MeshFormat",NULL));
  if (ierr != xf_OK) xf_Error(xf_FILE_READ_ERROR);
  // version file-type nbytes
  if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  ierr = sscanf(line, "%s %d %d\n", ver, &ftype, &nbytes);
  if (ierr != 3)
    return xf_Error(xf_FILE_READ_ERROR);
  
  if (strcmp(GMSHVERSION, ver) != 0) 
    xf_printf("Warning: Gmsh version different from supported %s version!\n",GMSHVERSION);
  //search for end of section for verification purposes
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_False, "$EndMeshFormat",NULL));
  if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR); 
  
  //Nodes
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_True, "$Nodes",NULL));
  if (ierr != xf_OK) xf_Error(xf_FILE_READ_ERROR);
  if (fscanf(fid, "%d\n",&nNodes) != 1)
    return xf_Error(xf_FILE_READ_ERROR);
  //read nodes into buffer
  ierr = xf_Error(xf_Alloc2((void ***)&nodes, 3, nNodes, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  for (n = 0; n < nNodes; n++) {
    if (fscanf(fid, "%d %lf %lf %lf\n",&j, &x, &y, &z) != 4)
      return xf_Error(xf_FILE_READ_ERROR);
    
    if (j > nNodes)//file error no support for sparse nodelist
      return xf_Error(xf_FILE_READ_ERROR);
    
    nodes[0][j-1] = x;//j-1 for C numbering
    nodes[1][j-1] = y;
    nodes[2][j-1] = z;
    
    //get maximum variation in each coordinate direction
    if (n == 0){
      for (i = 0; i < 3; i++){
        dmin[i] = dmax[i] = nodes[i][j-1];
      }
    }
    else {
      for (i = 0; i < 3; i++){
        if (nodes[i][j-1] < dmin[i])
          dmin[i] = nodes[i][j-1];
        if (nodes[i][j-1] > dmax[i])
          dmax[i] = nodes[i][j-1];
        TV[i] = dmax[i]-dmin[i];
      }
    }
  }//nNodes
  
  //search for end of section for verification purposes
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_False, "$EndNodes",NULL));
  if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR); 
  
  Mesh->Dim = 3;
  for (i = 0; i < 3; i++){
    if (TV[i] <= MEPS)
      Mesh->Dim--;
  }
  Mesh->nNode = nNodes;
  //alloc node coordinates
  ierr = xf_Error(xf_Alloc2((void ***)&Mesh->Coord, Mesh->nNode, 
                            Mesh->Dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  n = 0;
  for (i = 0; i < 3; i++) {
    if (TV[i] > MEPS){
      for (j = 0; j < nNodes; j++) {
        Mesh->Coord[j][n] = nodes[i][j];
      }
      n++;
    }
  }
  xf_Release2((void **)nodes);
  
  //Now, let's get the physical entities so we can relate the corresponding elements
  //$PhysicalNames
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_True, "$PhysicalNames",NULL));
  if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR);
  //number of names
  if (fscanf(fid, "%d\n",&np_entity) != 1)
    return xf_Error(xf_FILE_READ_ERROR);
  //store physical entities for indexing purposes
  ierr = xf_Error(xf_Alloc((void **)&physical_entity, np_entity, sizeof(xf_GmshPEntity)));
  if (ierr != xf_OK) return ierr;
 
  //index bug correct on Dec 2013, by Yu
  for (i = 0; i < np_entity; i++) {
    if (fscanf(fid, "%d %d %s\n",&dim, &ntags, line) != 3)
      return xf_Error(xf_FILE_READ_ERROR);
    physical_entity[i].dim = dim;      //element dimension
    physical_entity[i].Group = -1;  
    physical_entity[i].tag = ntags;    //element "physical" name tag
    //get names (without quotes--check Gmsh manual)
    for (n = 1; n < strlen(line)-1; n++) {
      buffer[n-1] = line[n];
    }
    buffer[n-1] = '\0';
    ierr = xf_Error(xf_AllocString(&physical_entity[i].Title, 
                                   xf_MAXSTRLEN, buffer));
    if (ierr != xf_OK) return ierr;
  }
  //search for end of section for verification purposes
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_False, "$EndPhysicalNames",NULL));
  if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR);
  
  //This structure stores the element types for reference
  ierr = xf_Error(xf_CreateGmshElemInfo(&GmshTypesInfo));
  if (ierr != xf_OK) return ierr;
  
  //$Elements
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_True, "$Elements",NULL));
  if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR);
  
  //number of items(includes cells, faces and edges)
  if (fscanf(fid, "%d\n",&nItems) != 1)
    return xf_Error(xf_FILE_READ_ERROR);
  
  DummyEGroup = NULL;
  DummyBFGroup = NULL;
  Mesh->nBFaceGroup = 0;
  Mesh->nElemGroup = 0;
  for (n = 0; n < nItems; n++) {
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (sscanf(line, "%d %d %d %d",&elm_number, &elm_type, &ntags, &p_entity) != 4)
      return xf_Error(xf_FILE_READ_ERROR);

    //match "physical" name tage to find index
    for(i = 0; i < np_entity; i++)
       if(physical_entity[i].tag == p_entity)
          break;

    p_entity = i; //was a bug before; fixed by Yu Lv

    //check element dimension
    if (physical_entity[p_entity].dim == Mesh->Dim) {//"volumetric" element
      //check if Element group exists
      found = xfe_False;
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
        if (GmshTypesInfo[elm_type].QOrder == DummyEGroup[egrp].QOrder &&
            GmshTypesInfo[elm_type].QBasis == DummyEGroup[egrp].QBasis) {
          found = xfe_True;
          break;//egrp is set to the correct value automatically
         }
       }
      if (!found) {//new combination of QOrder and QBasis
        Mesh->nElemGroup++;
        ierr = xf_Error(xf_ReAlloc((void **) &DummyEGroup, Mesh->nElemGroup, 
                                   sizeof(xf_ElemGroup)));
        if (ierr != xf_OK)  return ierr;
        //initialize new element group
        xf_InitElemGroup(DummyEGroup + egrp);
        
        DummyEGroup[egrp].QOrder = GmshTypesInfo[elm_type].QOrder;
        DummyEGroup[egrp].QBasis = GmshTypesInfo[elm_type].QBasis;
        ierr = xf_Error(xf_Order2nNode(DummyEGroup[egrp].QBasis, 
                                       DummyEGroup[egrp].QOrder, 
                                       &DummyEGroup[egrp].nNode));
        if (ierr != xf_OK) return ierr;
        
        //assuming no hanging nodes/faces
        ierr = xf_Error(xf_Alloc((void **) &DummyEGroup[egrp].nFace, 
                                 1, sizeof(int)));
        if (ierr!=xf_OK) return ierr;
        // number of faces for this element group
        ierr = xf_Error(xf_Basis2nFace(DummyEGroup[egrp].QBasis, 
                                       &DummyEGroup[egrp].nFace[0]));
        if (ierr != xf_OK) return ierr;
       }
      DummyEGroup[egrp].nElem++;
    }
    else if (physical_entity[p_entity].dim < Mesh->Dim){//Boundary
      //check if BFace group exists
      found = xfe_False;
      for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++) {
        if (strcmp(DummyBFGroup[fgrp].Title, 
                   physical_entity[p_entity].Title) == 0) {
          found = xfe_True;
          break;//fgrp is set to the correct value automatically
        }
      }
      if (!found) {//new BFGroup
        Mesh->nBFaceGroup++;
        ierr = xf_Error(xf_ReAlloc((void **) &DummyBFGroup, Mesh->nBFaceGroup, 
                                   sizeof(xf_BFaceGroup)));
        if (ierr != xf_OK)  return ierr;
        //store Title
        ierr = xf_Error(xf_AllocString(&DummyBFGroup[fgrp].Title, xf_MAXSTRLEN, 
                                       physical_entity[p_entity].Title));
        if (ierr != xf_OK) return ierr;
        DummyBFGroup[fgrp].nBFace = 0;
      }
      DummyBFGroup[fgrp].nBFace++;
      physical_entity[p_entity].Group = fgrp;
    }
    else 
      return xf_Error(xf_MESH_ERROR);
  }
  
  //check with output
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++) {
    printf("fgrp: %d nBFace: %d Title: %s\n",fgrp,DummyBFGroup[fgrp].nBFace,DummyBFGroup[fgrp].Title);
  }

  //verification
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_False, "$EndElements",NULL));
  if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR);
  
  //allocate Mesh structures
  //Boundaries
  ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup, Mesh->nBFaceGroup, 
                           sizeof(xf_BFaceGroup)));
  if (ierr != xf_OK) return ierr;
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++) {
    //nBFace
    Mesh->BFaceGroup[fgrp].nBFace = DummyBFGroup[fgrp].nBFace;
    DummyBFGroup[fgrp].nBFace = 0; //use this as counter later
    //Title
    ierr = xf_Error(xf_AllocString(&Mesh->BFaceGroup[fgrp].Title, xf_MAXSTRLEN, 
                                   DummyBFGroup[fgrp].Title));
    if (ierr != xf_OK) return ierr;
    //BFace
    ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup[fgrp].BFace, 
                             Mesh->BFaceGroup[fgrp].nBFace, sizeof(xf_BFace)));
    if (ierr != xf_OK)  return ierr;
  }
  //Element groups
  nElemTot = 0;
  ierr = xf_Error(xf_ReAlloc((void **) &Mesh->ElemGroup, Mesh->nElemGroup, 
                             sizeof(xf_ElemGroup)));
  if (ierr!=xf_OK) return ierr;
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
    //initialize element group (takes care of CutFlag and CutElemData)
    xf_InitElemGroup(Mesh->ElemGroup + egrp);
    //QOrder
    Mesh->ElemGroup[egrp].QOrder = DummyEGroup[egrp].QOrder;
    //QBasis
    Mesh->ElemGroup[egrp].QBasis = DummyEGroup[egrp].QBasis;
    //nElem
    Mesh->ElemGroup[egrp].nElem = DummyEGroup[egrp].nElem;
    DummyEGroup[egrp].nElem = 0;//we will use this as counter later
    nElemTot += Mesh->ElemGroup[egrp].nElem;
    //Allocate nFace
    ierr = xf_Error(xf_Alloc((void **) &Mesh->ElemGroup[egrp].nFace, 
                             Mesh->ElemGroup[egrp].nElem, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
    //nFace
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++) 
      Mesh->ElemGroup[egrp].nFace[elem] = DummyEGroup[egrp].nFace[0];
    //Face
    ierr = xf_Error(xf_VAlloc2((void ***) &Mesh->ElemGroup[egrp].Face, 
                               Mesh->ElemGroup[egrp].nElem, 
                               Mesh->ElemGroup[egrp].nFace, sizeof(xf_Face)));
    if (ierr!=xf_OK) return ierr;
    //nNode
    Mesh->ElemGroup[egrp].nNode = DummyEGroup[egrp].nNode;
    //Node
    ierr = xf_Error(xf_Alloc2((void ***) &Mesh->ElemGroup[egrp].Node, 
                              Mesh->ElemGroup[egrp].nElem, 
                              Mesh->ElemGroup[egrp].nNode, sizeof(int)));
    if (ierr!=xf_OK) return ierr;
  }

  //IFace
  /* Allocate IFace array (more than necessary, will resize at end) */
  ierr = xf_Error(xf_GetLocFaceMax(&nLocFaceMax));
  if (ierr != xf_OK) return ierr;
  
  // sanity check
  if (nLocFaceMax > 20) 
    xf_printf("Warning, nLocFaceMax=%d; remove this warning if this is ok.\n", nLocFaceMax);
  
  nIFaceMax = nElemTot * nLocFaceMax; // not dividing by 2 to allow non-conforming case
  ierr = xf_Error(xf_Alloc( (void **) &Mesh->IFace, nIFaceMax, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;
  
  Mesh->nIFace = 0; // running total of number of interior faces
  
  /* Create a node-based hash for looking up faces given nodes */
  ierr = xf_Error(xf_CreateNodeHash(Mesh->nNode, &Node2Hash));
  if (ierr != xf_OK) return ierr;
  //Loop through $Elements again and write info into memory
  nvec = NULL;
  edata = NULL;
  fvec = NULL;
  //rewind file and set position in the elements section
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_True, "$Elements",NULL));
  if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR);
  
  //number of items(includes cells, faces and edges)
  if (fscanf(fid, "%d\n",&nItems) != 1)
    return xf_Error(xf_FILE_READ_ERROR);
 
  for (n = 0; n < nItems; n++) {
    if (fgets(line, xf_MAXLINELEN, fid) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (sscanf(line, "%d %d %d %d",&elm_number, &elm_type,&ntags,&p_entity) != 4)
      return xf_Error(xf_FILE_READ_ERROR);
  
    //match tag with index
    for(i = 0; i < np_entity; i++)
       if(physical_entity[i].tag == p_entity)
          break;

    p_entity = i;
   
    //check support
    if (!GmshTypesInfo[elm_type].Supported)
      return xf_Error(xf_NOT_SUPPORTED);
    if (ntags < 1) return xf_Error(xf_FILE_READ_ERROR);
    //get nodes
    ierr = xf_Error(xf_ReAlloc((void **) &nvec, GmshTypesInfo[elm_type].nNode, 
                               sizeof(int)));
    if (ierr != xf_OK)  return ierr;
  
    /* 3 here is for the number of reads before tags */
    ierr = xf_Error(xf_ReAlloc((void **) &edata, 
                               3+ntags+GmshTypesInfo[elm_type].nNode, 
                               sizeof(int)));
    if (ierr != xf_OK)  return ierr;
    //scan values
    ierr = xf_Error(xf_ScanInt(line, 3+ntags+GmshTypesInfo[elm_type].nNode, edata));
    if (ierr != xf_OK) return ierr;
    
    for (j = 0; j < GmshTypesInfo[elm_type].nNode; j++)
      nvec[j] = edata[3+ntags+j]-1;//converting to C numbering
   
    if (physical_entity[p_entity].Group >= 0){//Boundary
      fgrp = physical_entity[p_entity].Group;
      bface = DummyBFGroup[fgrp].nBFace;
      DummyBFGroup[fgrp].nBFace++;
      //get number of linear nodes
      ierr = xf_Error(xf_Order2nNode(GmshTypesInfo[elm_type].QBasis, 1, &nfnode));
      if (ierr != xf_OK) return ierr;
      //add linear nodes to hash list (first nfnodes are the vertices)
      ierr = xf_Error(xf_AddFaceToHash(Node2Hash, nfnode, nvec, xfe_True, fgrp, 
                                       -1, bface, &hinfo, &exists));
      if (ierr != xf_OK) return ierr;
    }
    else if (physical_entity[p_entity].Group == -1){
      /* obviously used because it has an element in it */
      found = xfe_False;
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
        if (GmshTypesInfo[elm_type].QOrder == Mesh->ElemGroup[egrp].QOrder &&
            GmshTypesInfo[elm_type].QBasis == Mesh->ElemGroup[egrp].QBasis) {
          found = xfe_True;
          break;//egrp is set to the correct value automatically
        }
      }
      if (!found)
        return xf_Error(xf_CODE_LOGIC_ERROR);
      elem = DummyEGroup[egrp].nElem;
      DummyEGroup[egrp].nElem++;
     
      // number of faces for this element group
      ierr = xf_Error(xf_Basis2nFace(Mesh->ElemGroup[egrp].QBasis, &nface));
      if (ierr != xf_OK) return ierr;
      
      ierr = xf_Error(xf_ReAlloc((void **) &fvec, nface, sizeof(int)));
      if (ierr != xf_OK)  return ierr;
      
      // number of nodes for this element group (assuming complete set)
      ierr = xf_Error(xf_Order2nNode(Mesh->ElemGroup[egrp].QBasis, 
                                     Mesh->ElemGroup[egrp].QOrder, &nnode));
      if (ierr != xf_OK) return ierr;
      if (nnode > GmshTypesInfo[elm_type].nNode)
        need_interpolation = xfe_True;
      else if (nnode == GmshTypesInfo[elm_type].nNode)
        need_interpolation = xfe_False;
      else 
        return xf_Error(xf_MESH_ERROR);
     
      //convert nodal order
      ierr = xf_Error(xf_NodalOrderGmsh2Gri(GmshTypesInfo[elm_type], nvec));
      if (ierr != xf_OK)  return ierr;
      
      //interpolate interior nodes if necessary
      if (need_interpolation) {
        //temporary
        xf_printf("Interior nodes interpolation not implemented yet.\n");
        return xf_Error(xf_NOT_SUPPORTED);
      }
      
      // Add nodes
      for (j = 0; j<nnode; j++) Mesh->ElemGroup[egrp].Node[elem][j] = nvec[j];

    }//physical_entity[p_entity].Group == -1
    else 
      return xf_Error(xf_CODE_LOGIC_ERROR);

  }//nItems

  //The following loop is separate because of disordered writing into memory
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++) {
      /* loop over the faces and check the hash table */
      for (face=0; face< Mesh->ElemGroup[egrp].nFace[elem]; face++){
        
        // local nodes on face
        ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, 
                                         Mesh->ElemGroup[egrp].QOrder, 
                                         face, &nfnode, fvec));
        if (ierr != xf_OK) return ierr;
        
        // convert to global nodes
        for (j=0; j<nfnode; j++) fvec[j] = Mesh->ElemGroup[egrp].Node[elem][fvec[j]]; 
        
        /* Attempt to add face (fvec nodes) to hash list */
        ierr = xf_Error(xf_AddFaceToHash(Node2Hash, nfnode, fvec, xfe_False, 
                                         egrp, elem, face, &hinfo, &exists));
        if (ierr != xf_OK) return ierr;
        
        if (exists){ // face already exists in hash list
          if (hinfo->nVisit != 1){
            xf_printf("Error, more than two elements share a face, or a\n");
            xf_printf("boundary face is referenced by more than one element.\n");
            return xf_Error(xf_FILE_READ_ERROR);
          }
          
          // link elem to bface or to iface
          
          if (hinfo->BFlag == xfe_True){ // boundary face
            xf_InitBFace(Mesh->BFaceGroup[hinfo->Group].BFace + hinfo->Face);
            Mesh->BFaceGroup[hinfo->Group].BFace[hinfo->Face].ElemGroup = egrp;
            Mesh->BFaceGroup[hinfo->Group].BFace[hinfo->Face].Elem      = elem;
            Mesh->BFaceGroup[hinfo->Group].BFace[hinfo->Face].Face      = face;
            Mesh->ElemGroup[egrp].Face[elem][face].Group  = hinfo->Group;
            Mesh->ElemGroup[egrp].Face[elem][face].Number = hinfo->Face;
          }
          else{ // interior face
            xf_AddInteriorFace(Mesh, hinfo->Group, hinfo->Elem, hinfo->Face, egrp, elem, face);
          }
          
          ierr = xf_Error(xf_DeleteFaceFromHash(Node2Hash, nfnode, fvec));
          if (ierr != xf_OK) return ierr;
        }//face exists
        
      } // face
    }
  }
  
  // Make sure no faces are left in the hash
  for (i=0, face=0; i<Mesh->nNode; i++){
    if (Node2Hash[i] != NULL){
      hinfo = Node2Hash[i];
      while (hinfo != NULL){
        for (k=0; k<hinfo->nfnode; k++)
          xf_printf("%d ", hinfo->svec[k]+1);
        xf_printf("\n");
        face++;
        hinfo = hinfo->Next;
      }
    }
  } // i
  if (face != 0){
    xf_printf("Mesh connectivity error: the above %d face(s) remain(s) in the hash.\n\n",
              face);
    return xf_Error(xf_FILE_READ_ERROR);
  }
  
  // Resize IFace
  if (Mesh->nIFace > nIFaceMax) return xf_Error(xf_FILE_READ_ERROR);
  ierr = xf_Error(xf_ReAlloc( (void **) &Mesh->IFace, Mesh->nIFace, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;
  
  // update face orientation info  
  ierr = xf_Error(xf_UpdateFaceOrient(Mesh));
  if (ierr != xf_OK) return ierr;
  
  
  //search for end of section for verification purposes
  ierr = xf_Error(xf_SetFilePositionAfterString_Serial(fid, xfe_False, "$EndElements",NULL));
  if (ierr != xf_OK) return xf_Error(xf_FILE_READ_ERROR);
  
  //let's clean-up this mess!
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
    ierr = xf_Error(xf_DestroyElemGroup(DummyEGroup + egrp));
    if (ierr != xf_OK) return ierr;
  }
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++) {
    DummyBFGroup[fgrp].nBFace = 0;
    DummyBFGroup[fgrp].BFace = NULL;
    ierr = xf_Error(xf_DestroyBFaceGroup(DummyBFGroup + fgrp));
    if (ierr != xf_OK) return ierr;
  }
  
  // destroy lookup structure
  ierr = xf_Error(xf_DestroyGmshElemInfo(GmshTypesInfo));
  if (ierr != xf_OK) return ierr;
  
  
  xf_Release((void *)DummyEGroup);
  xf_Release((void *)DummyBFGroup);
  for (n = 0; n < np_entity; n++) {
    xf_Release((void*)physical_entity[n].Title);
  }
  xf_Release((void *)physical_entity);
  xf_Release((void *)edata);
  xf_Release((void *)fvec);
  xf_Release((void *)nvec);
  xf_Release((void *)Node2Hash);
  fclose(fid);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_NodalOrderGri2Gmsh
static int
xf_NodalOrderGri2Gmsh(xf_GmshElem elem_type, int *nvec)
{
  // on return nvec[gmsh index] = gri index
  int i, ierr;
  
  for (i = 0; i < elem_type.nNode; i++) nvec[elem_type.NodeOrder[i]] = i;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_GmshElemType
static int
xf_GmshElemType(xf_Mesh *Mesh, int egrp, int face, int *pelm_type)
{
  int ierr, QOrder;
  int maxorder = 5;     // q =   0  1  2  3  4  5
  int igmsh_Segment[6]       = {-1, 1, 8,26,27,28};
  int igmsh_Triangle[6]      = {-1, 2, 9,21,23,25};
  int igmsh_Tetrahedron[6]   = {-1, 4,11,29,30,31};
  int igmsh_Quadrilateral[6] = {-1, 3,10,-1,37,-1};
  int igmsh_Hexahedron[6]    = {-1, 5,-1,-1,93,-1};
  enum xfe_BasisType QBasis; 
  enum xfe_ShapeType Shape;
  
  // element geometry basis and order
  QBasis = Mesh->ElemGroup[egrp].QBasis;  
  QOrder = Mesh->ElemGroup[egrp].QOrder;
  
  // pull off element shape
  ierr = xf_Error(xf_Basis2Shape(QBasis, &Shape));
  if (ierr != xf_OK) return ierr;
  
  // convert element shape to face shape if face provided
  if (face >= 0){ // requesting face elm_type
    ierr = xf_Error(xf_FaceShape(Shape, face, &Shape));
    if (ierr != xf_OK) return ierr;
  }
  
  // gmsh only supports up to a certain order
  if (QOrder > maxorder) return xf_Error(xf_NOT_SUPPORTED);
  
  // do the conversion
  switch (Shape) {
    case xfe_Segment:       (*pelm_type) = igmsh_Segment[QOrder];       break;
    case xfe_Triangle:      (*pelm_type) = igmsh_Triangle[QOrder];      break;
    case xfe_Tetrahedron:   (*pelm_type) = igmsh_Tetrahedron[QOrder];   break;
    case xfe_Quadrilateral: (*pelm_type) = igmsh_Quadrilateral[QOrder]; break;
    case xfe_Hexahedron:    (*pelm_type) = igmsh_Hexahedron[QOrder];    break;
    default: return xf_Error(xf_NOT_SUPPORTED); break;
  }
  
  // some types are not yet supported; indicated by -1 elem type
  if ((*pelm_type) < 0) return xf_Error(xf_NOT_SUPPORTED);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_WriteGmshFile
int
xf_WriteGmshFile(xf_Mesh *Mesh, char *OutputFile)
{
  int ierr, i, j, elem_tot, nElemTot, nBFaceTot, elm_type;
  int egrp, elem, node, fgrp, face, faceref, *nvec, *fvec, nfvec, count, n_tags;
  int p_entity, g_entity, nfnode;
  xf_GmshElem *GmshTypesInfo;
  FILE *fid;
  
  if (OutputFile == NULL) 
    return xf_Error(xf_INPUT_ERROR);
  
  if ((fid =fopen(OutputFile ,"w")) == NULL) 
    return xf_Error(xf_FILE_WRITE_ERROR);
  
  // create lookup table for gmsh element types
  ierr = xf_Error(xf_CreateGmshElemInfo(&GmshTypesInfo));
  if (ierr != xf_OK) return ierr;
  
  //$MeshFormat
  fprintf(fid,"$MeshFormat\n");
  fprintf(fid,"%s 0 %zu\n",GMSHVERSION,sizeof(real));
  fprintf(fid,"$EndMeshFormat\n");
  
  //$Nodes
  fprintf(fid,"$Nodes\n");
  fprintf(fid,"%d\n",Mesh->nNode);
  for (i = 0; i < Mesh->nNode; i++) {
    fprintf(fid, "%d ",i+1);
    for (j = 0; j < Mesh->Dim; j++) {
      fprintf(fid, "%.15E ",Mesh->Coord[i][j]);
    }
    //complete dimentions with zeros
    for (j = Mesh->Dim; j < 3; j++) {
      fprintf(fid, "0 ");
    }
    fprintf(fid, "\n");
  }
  fprintf(fid,"$EndNodes\n");
  
  //$Elements
  //count the total number of elements to write
  nElemTot = 0;
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++)
    nElemTot += Mesh->ElemGroup[egrp].nElem;
  nBFaceTot = 0;
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++)
    nBFaceTot += Mesh->BFaceGroup[fgrp].nBFace;
  elem_tot = nElemTot+nBFaceTot;
  fprintf(fid, "$Elements\n");
  fprintf(fid, "%d\n",elem_tot);
  count = 0;
  nvec = NULL;
  fvec = NULL;
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++) {
    // get gmsh element-type
    ierr = xf_Error(xf_GmshElemType(Mesh, egrp, -1, &elm_type));
    if (ierr != xf_OK) return ierr;
    //verification
    if (GmshTypesInfo[elm_type].nNode != Mesh->ElemGroup[egrp].nNode){
      xf_printf("Number of nodes incompatible between element types.\n");
      return xf_Error(xf_CODE_LOGIC_ERROR);
    }
    // reallocate nvec buffer
    ierr = xf_Error(xf_ReAlloc((void **) &nvec, Mesh->ElemGroup[egrp].nNode, 
                               sizeof(int)));
    if (ierr != xf_OK)  return ierr;
    
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++) {
      // element number:
      fprintf(fid, "%d ",count+1);
      n_tags = 2;//look at gmsh manual
      // by default the physical entity for the "volumetric" elements will be set to 0
      p_entity = 0; 
      g_entity = 0; //not used
      fprintf(fid, "%d %d %d %d ",elm_type, n_tags, p_entity, g_entity);
      
      // convert gri ordering to gmsh format
      ierr = xf_Error(xf_NodalOrderGri2Gmsh(GmshTypesInfo[elm_type], nvec));
      if (ierr != xf_OK)  return ierr;
      
      // write out nodes
      for (node = 0; node < Mesh->ElemGroup[egrp].nNode; node++)
        fprintf(fid,"%d ",Mesh->ElemGroup[egrp].Node[elem][nvec[node]]+1);
      
      //terminate line
      fprintf(fid, "\n");
      
      count++;
    }//elem
  }//egrp
  
  p_entity++;
  
  nfvec = -1;
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++) {
    for (face = 0; face < Mesh->BFaceGroup[fgrp].nBFace; face++) {
      
      // element information
      elem = Mesh->BFaceGroup[fgrp].BFace[face].Elem;
      egrp = Mesh->BFaceGroup[fgrp].BFace[face].ElemGroup;
      faceref = Mesh->BFaceGroup[fgrp].BFace[face].Face;
      
      // re-allocate nvec, fvec if necessary
      if (Mesh->ElemGroup[egrp].nNode > nfvec){
        nfvec = Mesh->ElemGroup[egrp].nNode;
        ierr = xf_Error(xf_ReAlloc((void **) &nvec, nfvec, sizeof(int)));
        if (ierr != xf_OK)  return ierr;
        ierr = xf_Error(xf_ReAlloc((void **) &fvec, nfvec, sizeof(int)));
        if (ierr != xf_OK)  return ierr;
      }
      
      // get nodes on face (reference coordinate)
      ierr = xf_Error(xf_NodesOnFace(Mesh->ElemGroup[egrp].QBasis, 
                                     Mesh->ElemGroup[egrp].QOrder, faceref, 
                                     &nfnode, fvec));
      if (ierr != xf_OK)  return ierr;
      
      //get gmsh element-type
      ierr = xf_Error(xf_GmshElemType(Mesh, egrp, faceref, &elm_type));
      if (ierr != xf_OK) return ierr;
      
      //reorder nodes in gmsh standard
      ierr = xf_Error(xf_NodalOrderGri2Gmsh(GmshTypesInfo[elm_type], nvec));
      if (ierr != xf_OK)  return ierr;
      
      n_tags = 2;//look at gmsh manual
      g_entity = 0; //not used
      fprintf(fid, "%d %d %d %d %d ", count+1, elm_type, n_tags, p_entity, g_entity);
      for (node = 0; node < nfnode; node++)
        fprintf(fid,"%d ",Mesh->ElemGroup[egrp].Node[elem][fvec[nvec[node]]]+1);
      //terminate line
      fprintf(fid, "\n");
      count++;
    }//face
    p_entity++;
  }
  
  //verification
  if (count != elem_tot)
    return xf_Error(xf_CODE_LOGIC_ERROR);
  fprintf(fid, "$EndElements\n");
  
  //$PhysicalNames
  fprintf(fid, "$PhysicalNames\n");
  fprintf(fid, "%d\n",p_entity);
  fprintf(fid, "%d 0 \"MeshInterior\"\n",Mesh->Dim);
  for (fgrp=0; fgrp < Mesh->nBFaceGroup; fgrp++) {
    fprintf(fid, "%d %d \"%s\"\n",Mesh->Dim-1, fgrp+1, Mesh->BFaceGroup[fgrp].Title);
  }
  fprintf(fid, "$EndPhysicalNames\n");
  
  fclose(fid);
  
  
  // destroy lookup structure
  ierr = xf_Error(xf_DestroyGmshElemInfo(GmshTypesInfo));
  if (ierr != xf_OK) return ierr;
  
  xf_Release((void *)nvec);
  xf_Release((void *)fvec);
  
  return xf_OK;
}
