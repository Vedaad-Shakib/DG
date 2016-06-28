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
  FILE:  xf_MeshBamg.c

  This file contains functions for interfacing with BAMG -- the
  Bi-dimensional Anisotropic Mesh Generator from Inria (now in
  FreeFEM++).

*/

#include <stdlib.h>


/******************************************************************/
//   FUNCTION Definition: xf_CreateBamgGeomPoints
int 
xf_CreateBamgGeomPoints(xf_BamgGeomPoints **pBamgPoints)
{
  int ierr;

  ierr = xf_Error(xf_Alloc( (void **) pBamgPoints, 1, 
			    sizeof(xf_BamgGeomPoints) ));
  if (ierr != xf_OK) return ierr;

  (*pBamgPoints)->nVert   = 0;
  (*pBamgPoints)->Vert    = NULL;
  (*pBamgPoints)->Tangent = NULL;
  (*pBamgPoints)->nEdge   = 0;
  (*pBamgPoints)->Edge    = NULL;
  (*pBamgPoints)->nCorner = 0;
  (*pBamgPoints)->Corner  = NULL;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyBamgGeomPoints
int 
xf_DestroyBamgGeomPoints(xf_BamgGeomPoints *BamgPoints)
{
  int ierr;

  xf_Release( (void *) BamgPoints->Vert);
  xf_Release( (void *) BamgPoints->Tangent);
  xf_Release( (void *) BamgPoints->Edge);
  xf_Release( (void *) BamgPoints->Corner);

  xf_Release( (void *) BamgPoints);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ReadBamgFile
int 
xf_ReadBamgFile( const char *InputFile, char **InputStrings, xf_Mesh *Mesh)
{
  int ierr, nNode, nelemtot, dim;
  int i, j, k, fi, nBFG, iBFG, nBFace, nfnode;
  int nbfacetot, *vnbface = NULL;
  int nElem, nface, nnode;
  int egrp, elem, face, iBFace, nIFaceMax;
  int iString = 0;
  int fvec[4], ivec[4];
  char *line, line0[xf_MAXLONGLINELEN];
  char title[xf_MAXSTRLEN];
  FILE *fbamg = NULL;
  enum xfe_Bool exists;
  xf_FaceHash **Node2Hash, *hinfo;
  real Volume;
  
  // Check input
  if (((InputFile == NULL) && (InputStrings == NULL)) || 
      ((InputFile != NULL) && (InputStrings != NULL))) 
    return xf_Error(xf_INPUT_ERROR);

  if (InputFile != NULL){
    // Open InputFile
    fbamg = fopen(InputFile, "r");
    if(!fbamg){
      xf_printf("File %s not found.\n\n", InputFile);
      return xf_Error(xf_FILE_READ_ERROR);
    }
    xf_printf("Reading Bamg file.\n");
  }
  else{
    // Start at the first string
    iString = 0;
  }
  
  /* Set dimension */
  Mesh->Dim = dim = 2;

  /* Vertices (nodes) */
  nNode = -1;
  do{
    /* Read in line */
    ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
    if (strncmp(line, "Vertices", 8) == 0){
      ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ScanInt(line, 1, &nNode));
      if (ierr != xf_OK) return ierr;
    }
  } while (nNode == -1);

  /* Read in nodes */
  Mesh->nNode = nNode;
    
  ierr = xf_Error(xf_Alloc2((void ***) &Mesh->Coord, nNode, dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  for(i=0; i<nNode; i++){
    ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_ScanReal(line, dim, Mesh->Coord[i]));
    if (ierr != xf_OK) return ierr;
  }

  /* Create a node-based hash for looking up faces given nodes */
  ierr = xf_Error(xf_CreateNodeHash(nNode, &Node2Hash));
  if (ierr != xf_OK) return ierr;


  /* Edges (boundary face groups) */
  ierr = xf_Error(xf_RewindFileOrStrings(fbamg, &iString));
  if (ierr != xf_OK) return ierr;
  nbfacetot = -1;
  do{
    /* Read in line */
    ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
    if (strncmp(line, "Edges", 5) == 0){
      ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ScanInt(line, 1, &nbfacetot));
      if (ierr != xf_OK) return ierr;
    }
  } while (nbfacetot == -1);

  /* First determine number of boundary face groups and bfaces in each group */
  nBFG    = 0;
  vnbface = NULL;
  for (fi=0; fi<nbfacetot; fi++){
    /* Read in line */
    ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ScanInt(line, 3, ivec));
    if (ierr != xf_OK) return ierr;
    if ((iBFG = ivec[2]) < 0) return xf_Error(xf_FILE_READ_ERROR);
    if (iBFG > nBFG){ // note, iBFG is 1-based in bamg file
      ierr = xf_Error(xf_ReAlloc((void **) &vnbface, iBFG, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      for (j=nBFG; j<iBFG; j++) vnbface[j] = 0;
      nBFG = iBFG;
    }
    vnbface[iBFG-1] += 1;
  }

  if (nBFG == 0) return xf_Error(xf_FILE_READ_ERROR);
  
  Mesh->BFaceGroup = NULL;
 

  /* Allocate nBFG bounday face groups */
  Mesh->nBFaceGroup = nBFG;
  ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup, nBFG, sizeof(xf_BFaceGroup)));
  if (ierr != xf_OK) return ierr;
  for (iBFG=0; iBFG<nBFG; iBFG++){
    Mesh->BFaceGroup[iBFG].nBFace = 0; // will be counter on second read
    sprintf(title, "BoundaryGroup%d\0", iBFG);
     /* Create name for bfgroup */
    ierr = xf_Error(xf_AllocString(&Mesh->BFaceGroup[iBFG].Title, xf_MAXSTRLEN, title));
    if (ierr != xf_OK) return ierr;
    /* Initialize .BFace */
    ierr = xf_Error(xf_Alloc((void **) &Mesh->BFaceGroup[iBFG].BFace, vnbface[iBFG], sizeof(xf_BFace)));
    if (ierr != xf_OK)  return ierr;
  }


  /* Read boundary faces, add to hash list */
  ierr = xf_Error(xf_RewindFileOrStrings(fbamg, &iString));
  if (ierr != xf_OK) return ierr;
  nbfacetot = -1;
  do{
    /* Read in line */
    ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
    if (strncmp(line, "Edges", 5) == 0){
      ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ScanInt(line, 1, &nbfacetot));
      if (ierr != xf_OK) return ierr;
    }
  } while (nbfacetot == -1);

  
  for (fi=0; fi<nbfacetot; fi++){
    /* Read in line */
    ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ScanInt(line, 3, ivec));
    if (ierr != xf_OK) return ierr;
    if ((iBFG = (ivec[2]-1)) < 0) return xf_Error(xf_FILE_READ_ERROR);
    
    nfnode = 2; // number of nodes on a face; we're dealing with triangles

    for (i=0; i<nfnode; i++){
      ivec[i]--;  // convert to 0-starting index
      if ((ivec[i] < 0) || (ivec[i] >= nNode)){
	xf_printf("Node index = %d out of range when reading boundary edges.\n", ivec[i]+1);
	xf_printf("Note, the indexing should start at 1.\n");
	return xf_Error(xf_FILE_READ_ERROR);
      }
    }

    iBFace = Mesh->BFaceGroup[iBFG].nBFace++;

    ierr = xf_Error(xf_AddFaceToHash(Node2Hash, nfnode, ivec, xfe_True, iBFG, 
				     -1, iBFace, &hinfo, &exists));
    if (ierr != xf_OK) return ierr;
      
    if (exists){ // face should not exist in hash list
      xf_printf("Error, boundary face repeated:\n");
      for (i=0; i<nfnode; i++) xf_printf("%d ", ivec[i]+1);
      xf_printf("\n (1-base numbering)\n");
      return xf_Error(xf_FILE_READ_ERROR);
    }
    
  } // fi


  // boundary faces consistency check
  for (iBFG=0; iBFG<nBFG; iBFG++)
    if (Mesh->BFaceGroup[iBFG].nBFace != vnbface[iBFG]) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  xf_Release( (void *) vnbface);

  /* Triangles (elements) */
  ierr = xf_Error(xf_RewindFileOrStrings(fbamg, &iString));
  if (ierr != xf_OK) return ierr;
  nelemtot = -1;
  do{
    /* Read in line */
    ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
    if (strncmp(line, "Triangles", 9) == 0){
      ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_ScanInt(line, 1, &nelemtot));
      if (ierr != xf_OK) return ierr;
    }
  } while (nelemtot == -1);


  // put all elements in one element group
  Mesh->nElemGroup = 1;
  ierr = xf_Error(xf_ReAlloc((void **) &Mesh->ElemGroup, Mesh->nElemGroup, 
			     sizeof(xf_ElemGroup)));
  if (ierr!=xf_OK) return ierr;

  egrp  = 0;
  nElem = nelemtot;
  nface = 3;
  nnode = 3;
  
  xf_InitElemGroup(Mesh->ElemGroup + egrp);

  Mesh->ElemGroup[egrp].QBasis = xfe_TriLagrange;
  Mesh->ElemGroup[egrp].QOrder = 1;
  Mesh->ElemGroup[egrp].nElem  = nElem;
  ierr = xf_Error(xf_Alloc((void **) &Mesh->ElemGroup[egrp].nFace, nElem, sizeof(int)));
  if (ierr!=xf_OK) return ierr;
  for (elem=0; elem<nElem; elem++) Mesh->ElemGroup[egrp].nFace[elem] = nface;
  Mesh->ElemGroup[egrp].nNode = nnode;
  ierr = xf_Error(xf_VAlloc2((void ***) &Mesh->ElemGroup[egrp].Face, nElem, 
			     Mesh->ElemGroup[egrp].nFace, sizeof(xf_Face)));
  if (ierr!=xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc2((void ***) &Mesh->ElemGroup[egrp].Node, nElem, nnode, sizeof(int)));
  if (ierr!=xf_OK) return ierr;

  
  /* Allocate IFace array (exactly the right amount) */
  nIFaceMax = (nelemtot*3-nbfacetot)/2; 
  ierr = xf_Error(xf_Alloc( (void **) &Mesh->IFace, nIFaceMax, sizeof(xf_IFace)));
  if (ierr != xf_OK) return ierr;

  Mesh->nIFace = 0; // running total of number of interior faces
  
  // Loop over elements and read them in, node by node
  for (elem=0; elem<nelemtot; elem++){
  
    // read in nodes for an element
    ierr = xf_Error(xf_LineFromFileOrStrings(fbamg, InputStrings, &iString, line0, &line));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ScanInt(line, nnode, ivec));
    if (ierr != xf_OK) return ierr;

    for (i=0; i<nnode; i++){
      ivec[i] -= 1;  // to C numbering
      if ((ivec[i] < 0) || (ivec[i] >= nNode)){
	xf_printf("Node index = %d out of range when reading elements.\n", ivec[i]+1);
	xf_printf("Note, the indexing should start at 1.\n");
	return xf_Error(xf_FILE_READ_ERROR);
      }
    }
    
    // make sure element volume (jacobian at quad points) is positive
    ierr = xf_LinearElemJacobian(xfe_TriLagrange, 1, ivec, Mesh->Coord, &Volume);
    if (ierr != xf_OK) return xf_Error(ierr);
    if ((ierr == xf_OK) && (Volume <= 0.)){
      xf_printf("Swapping nodes to fix negative volume in element = %d\n", elem);
      swap(ivec[0], ivec[1], k);
    }
    
    // Add nodes
    for (k = 0; k<nnode; k++) Mesh->ElemGroup[egrp].Node[elem][k] = ivec[k];
    
    /* loop over the faces and check the hash table */
    for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){

      // local nodes on face
      ierr = xf_Error(xf_Q1NodesOnFace(xfe_TriLagrange, 1, face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;
      
      // convert to global nodes
      for (k=0; k<nfnode; k++) fvec[k] = ivec[fvec[k]]; 
            
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
      }
      
    } // face
    
  } // elem
 
      
  // Make sure no faces are left in the hash; print out remaining faces if any
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

  // Check number of interior faces
  if (Mesh->nIFace != nIFaceMax) return xf_Error(xf_FILE_READ_ERROR);

  // update face orientation info  
  ierr = xf_Error(xf_UpdateFaceOrient(Mesh));
  if (ierr != xf_OK) return ierr;

  xf_Release((void *) Node2Hash);

  if (fbamg != NULL) fclose(fbamg);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_WriteBamgGeometry_Points
int 
xf_WriteBamgGeometry_Points( xf_BamgGeomPoints *BP, char *OutputFile )
{
  int i;
  FILE *fbg;

  // open file
  if ((fbg = fopen(OutputFile ,"w"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);

  // header info
  fprintf(fbg, "MeshVersionFormatted 0\n\n");
  fprintf(fbg, "Dimension 2\n\n");

  // write vertices
  fprintf(fbg, "\nVertices\n%d\n", BP->nVert);
  for (i=0; i<BP->nVert; i++)
    fprintf(fbg, "%.15E %.15E %d\n", BP->Vert[2*i+0], BP->Vert[2*i+1], 1);
  
  // write edges
  fprintf(fbg, "\nEdges\n%d\n", BP->nEdge);
  for (i=0; i<BP->nEdge; i++)
    fprintf(fbg, "%d %d %d\n", BP->Edge[3*i+0], BP->Edge[3*i+1], BP->Edge[3*i+2]);

  // domain interior is always on left of boundary loops
  fprintf(fbg, "\nSubDomain\n%d\n", 1);
  fprintf(fbg, "%d %d %d %d\n", 2, 1, 1, 1);

  // write corners
  if (BP->nCorner > 0){
    fprintf(fbg, "\nCorners\n%d\n", BP->nCorner);
    for (i=0; i<BP->nCorner; i++)
      fprintf(fbg, "%d\n", BP->Corner[i]);
  }
  
  // close file
  fclose(fbg);

  return xf_OK;

}


/******************************************************************/
//   FUNCTION Definition: xf_WriteBamgFile
int 
xf_WriteBamgFile( xf_Mesh *Mesh, enum xfe_Bool AnisotropicFlag,
		  real *NodeMetric, char *OutputFile )
{
  int ierr, i, k, dim, dim2;
  int nelemtot, nbfacetot, nNode, nnode;
  int ibfgrp, ibface, egrp, elem, face, node;
  int nfnode, nvec[3], fvec[3];
  int rnk, IM[3] = {0,1,3};
  int *Node, *NodeFlag = NULL;
  char MetricFile[xf_MAXSTRLEN];
  //char gfile[xf_MAXSTRLEN];
  xf_ElemGroup *EG;
  xf_BFaceGroup *BFGroup;
  FILE *fbamg, *fmetric;


  dim  = 2;
  dim2 = dim*dim;

  // first, write a geometry file
  //sprintf(gfile, "%s.geom", OutputFile);
  //return xf_Error(xf_NOT_SUPPORTED);

  // write mesh file, referencing the geometry
  if ((fbamg = fopen(OutputFile ,"w"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);

  // header info
  fprintf(fbamg, "MeshVersionFormatted 0\n\n");
  fprintf(fbamg, "Dimension 2\n\n");
  //fprintf(fbamg, "Geometry\n\"%s\"", gfile); // reference to geometry

  // create a node flag, in case we have high-order (geom) elements
  ierr = xf_Error(xf_Alloc((void **) &NodeFlag, Mesh->nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<Mesh->nNode; i++) NodeFlag[i] = -1;
  

  /** print out node coordinates (all except Q>1 nodes) **/
  nNode = 0;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0, EG=Mesh->ElemGroup+egrp; elem<EG->nElem; elem++){
      if (EG->QBasis != xfe_TriLagrange) return xf_Error(xf_NOT_SUPPORTED);
      ierr = xf_Error(xf_Q1Nodes(EG->QBasis, EG->QOrder, &nnode, nvec));
      if (ierr != xf_OK) return ierr;
      if (nnode != 3) return xf_Error(xf_CODE_LOGIC_ERROR);
      for (k=0; k<nnode; k++){
	node = EG->Node[elem][nvec[k]];
	if (NodeFlag[node] == -1){
	  NodeFlag[node] = 1;
	  nNode++;
	}
      }
    } // elem
  
  fprintf(fbamg, "\nVertices\n%d\n", nNode);
  
  nNode=0;
  for (i=0; i<Mesh->nNode; i++)
    if (NodeFlag[i] == 1){
      fprintf(fbamg, "%.15E %.15E %d\n", Mesh->Coord[i][0], Mesh->Coord[i][1], 1);
      NodeFlag[i] = nNode++;
    }
      
  
  /** print out boundary edges **/
  
  for (ibfgrp=0, nbfacetot=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++)
    nbfacetot += Mesh->BFaceGroup[ibfgrp].nBFace;

  fprintf(fbamg, "\nEdges\n%d\n", nbfacetot);

  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    BFGroup = Mesh->BFaceGroup + ibfgrp;
    for (ibface=0; ibface<BFGroup->nBFace; ibface++){
      egrp = BFGroup->BFace[ibface].ElemGroup;
      elem = BFGroup->BFace[ibface].Elem;
      face = BFGroup->BFace[ibface].Face;
      ierr = xf_Error(xf_Q1NodesOnFace(Mesh->ElemGroup[egrp].QBasis, 
				       Mesh->ElemGroup[egrp].QOrder, face, 
				       &nfnode, nvec));
      if (ierr != xf_OK) return ierr;
      Node = Mesh->ElemGroup[egrp].Node[elem];
      if (nfnode != 2) return xf_Error(xf_NOT_SUPPORTED);
      // +1 required to get from C to Fortran indexing
      fprintf(fbamg, "%d %d %d\n", NodeFlag[Node[nvec[0]]]+1, NodeFlag[Node[nvec[1]]]+1, ibfgrp+1);
    } // ibface
  } // ibfgrp

  /** print out triangles **/

  // total number of elements
  ierr = xf_Error(xf_GetnElem(Mesh, NULL, &nelemtot));
  if (ierr != xf_OK) return ierr;

  fprintf(fbamg, "\nTriangles\n%d\n", nelemtot);
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0, EG=Mesh->ElemGroup+egrp; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      ierr = xf_Error(xf_Q1Nodes(EG->QBasis, EG->QOrder, &nnode, nvec));
      if (ierr != xf_OK) return ierr;
      for (k=0; k<3; k++)
	fprintf(fbamg, "%d ", NodeFlag[EG->Node[elem][nvec[k]]]+1);
      fprintf(fbamg, "1\n");
    }

  /** SubDomainFromGeom **/

  fprintf(fbamg, "\nSubDomainFromGeom %d\n", 1);
  fprintf(fbamg, "%d %d %d %d\n", 2, 1, 1, 1);

  fprintf(fbamg, "\nEnd\n");

  fclose(fbamg);


  // Write Metric if requested
  if (NodeMetric != NULL){
    rnk = ((AnisotropicFlag) ? 3 : 1);
    
    // open metric file
    sprintf(MetricFile, "%s.metric", OutputFile);
    if ((fmetric = fopen(MetricFile ,"w"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);

    // write out header
    fprintf(fmetric, "%d %d\n", nNode, rnk);
    
    // write out metric entries
    for (i=0; i<Mesh->nNode; i++)
      if (NodeFlag[i] != -1){
	for (k=0; k<rnk; k++)
	  fprintf(fmetric, "%.15E ", NodeMetric[i*dim2+IM[k]]);
	fprintf(fmetric, "\n");
      }

    fclose(fmetric);
  }
 
  // release memory
  xf_Release( (void *) NodeFlag);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MetricElem2Node
static int 
xf_MetricElem2Node( xf_Mesh *Mesh, xf_Vector *EM, real *NM )
{
  int ierr;
  int egrp, elem, j, k, node;
  int dim, dim2;
  int *counter = NULL;
  real *E;

  dim  = Mesh->Dim;
  dim2 = dim*dim;

  // create and zero out a counter
  ierr = xf_Error(xf_Alloc( (void **) &counter, Mesh->nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<Mesh->nNode; k++) counter[k] = 0;

  // zero out NM
  for (k=0; k<Mesh->nNode*dim2; k++) NM[k] = 0.;

  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      for (j=0; j<Mesh->ElemGroup[egrp].nNode; j++){
	node = Mesh->ElemGroup[egrp].Node[elem][j];
	E = EM->GenArray[egrp].rValue[elem];
	for (k=0; k<dim2; k++) NM[node*dim2+k] += E[k];
	counter[node] += 1;
      }
  
  // set NM to average of adjacent elements (above was just sum)
  for (node=0; node<Mesh->nNode; node++){
    if (counter[node] != 0) 
      for (k=0; k<dim2; k++) NM[node*dim2+k] /= ((real) counter[node]);
  }

  // release memory
  xf_Release( (void *) counter);
  

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_RefineMeshBamg
int 
xf_RefineMeshBamg( xf_All *All, enum xfe_Bool AnisotropicFlag,
		   xf_Vector *EM, real *NMin, xf_Mesh **pOldMesh )
{
  int ierr, j, nNode;
  int dim, dim2;
  int ibfgrp;

  // hardcoded default values for calling BAMG
  char BamgIn[]  = "BamgIn.bamg";   // BAMG gets this file as input
  char BamgOut[] = "BamgOut.bamg";  // BAMG outputs this file
  int  MaxVertices = 100000;        // max # of vertices
  int  NbSmooth    = 3;             // number of smoothing iterations
  real BamgRatio   = 0.0;           // < 1.1 means no smoothing
  real MaxAniso    = 1e6;           // max level of anisotropy  

  char cmd[xf_MAXLINELEN];

  real *NM = NULL;  // node-based metric
  real h, hmin, hmax;
  real a, b, c, l1, l2;

  xf_Mesh *Mesh, *NewMesh;
  
  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  if (dim != 2) return xf_Error(xf_INPUT_ERROR);
  dim2 = dim*dim;

  // print out parameters
  xf_printf("AnisotropicFlag = %d\n",    AnisotropicFlag);
  xf_printf("BamgRatio       = %.15e\n", BamgRatio);
  xf_printf("MaxAniso        = %.15e\n", MaxAniso);
  xf_printf("MaxVertices     = %d\n",    MaxVertices);
  xf_printf("NbSmooth        = %d\n",    NbSmooth);

  // number of nodes
  nNode = Mesh->nNode;

  // set NM = node-based metric
  if (EM != NULL){
    ierr = xf_Error(xf_Alloc( (void **) &NM, nNode*dim2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    // calculate node-based metric from element-based metric
    ierr = xf_Error(xf_MetricElem2Node(Mesh, EM, NM));
    if (ierr != xf_OK) return ierr;
  }
  else NM = NMin;

  // Write out BAMG file with metric
  ierr = xf_Error(xf_WriteBamgFile(Mesh, AnisotropicFlag, NM, BamgIn));
  if (ierr != xf_OK) return ierr;

  /* determine hmin and hmax for BAMG input */
  hmin = 1e30;
  hmax = -1;
  if (AnisotropicFlag){
    for (j=0; j < nNode; j++){
      a = NM[j*dim2+0];  b = NM[j*dim2+1];  c = NM[j*dim2+3];
      if ((a > 0.0) || (c > 0.0)){
	l1 = 0.5*(a+c - sqrt(4.0*b*b + (a-c)*(a-c)) );
	l2 = 0.5*(a+c + sqrt(4.0*b*b + (a-c)*(a-c)) );
	if (l1 > 0.) hmax = max(hmax, 1.0/sqrt(l1));
	if (l2 > 0.) hmin = min(hmin, 1.0/sqrt(l2));
      }
    }
  }
  else{ // isotropic
    for (j=0; j < nNode; j++){
      h = NM[j*dim2];
      if (h > 0.){
	hmin = min(hmin, h);
	hmax = max(hmax, h);
      }
    }
  }

  /* check to see whether command processor exists */
  if (system(NULL) == 0){
    xf_printf("Error, command processor does not exist.  Cannot call BAMG.\n");
    return xf_Error(xf_SYSTEM_ERROR);
  }
  

  /* System call to BAMG */
  sprintf(cmd, "./bamg -b %s -M %s.metric -splitpbedge -anisomax %g -hmin %.12E -hmax %.12E -ratio %g -nbv %d -NbSmooth %d -noKeepBackVertices -o %s", 
	  BamgIn, BamgIn, MaxAniso, hmin, hmax, BamgRatio, MaxVertices, NbSmooth, BamgOut);
  xf_printf("Calling BAMG with the command:\n%s\n\n", cmd);
  xf_printf("----------------------------------------------------\n");
  ierr = system(cmd);
  xf_printf("----------------------------------------------------\n\n");
  xf_printf("BAMG finished.\n");
  
  if (ierr != 0){
    xf_printf("The bamg call could not be executed.  Make sure bamg is in working directory.\n");
    return xf_Error(xf_SYSTEM_ERROR);
  }


  // Create blank-slate mesh for reading BAMG output
  ierr = xf_Error(xf_CreateMesh(&NewMesh));
  if (ierr != xf_OK) return ierr;

  // Read BAMG output mesh
  ierr = xf_Error( xf_ReadBamgFile(BamgOut, NULL, NewMesh) );
  if (ierr != xf_OK) return ierr;

  // Rename boundaries in NewMesh
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    // release existing title in NewMesh
    xf_Release((void *) NewMesh->BFaceGroup[ibfgrp].Title);
    // allocate new title based on previous one
    ierr = xf_Error(xf_AllocString(&NewMesh->BFaceGroup[ibfgrp].Title, xf_MAXSTRLEN, 
				   Mesh->BFaceGroup[ibfgrp].Title));
    if (ierr != xf_OK) return ierr;
  }  

  // Destroy existing mesh if no request to keep it
  if (pOldMesh == NULL){
    ierr = xf_Error(xf_DestroyMesh(All->Mesh));
    if (ierr != xf_OK) return ierr;
  }
  else (*pOldMesh) = All->Mesh;

  // store new mesh in All
  Mesh = All->Mesh = NewMesh;


  // TODO: note HO boundaries in input and curve output boundaries accordingly

  // destroy NM if created it
  if (EM != NULL) xf_Release((void *) NM);

  return xf_OK;
}

