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
  FILE:  xf_StructuredHO.c

  This program makes a high-order mesh from a structured linear mesh

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_Arg.h"
#include "xf_MeshTools.h"
#include "xf_Math.h"

// maximum "Q" for high-order
#define MQ 4

// threshold for corner eval, cos(angle)
#define DP_THRESH 0.7

/******************************************************************/
//  FUNCTION Definition: xf_OppositeFace
static int 
xf_OppositeFace(enum xfe_BasisType Basis, int f0, int *f1)
{
/*

PURPOSE: 

  Determines face opposite f0 on elem of type Basis
  
INPUTS:

  QBasis : element basis type
  f0 : face number on element

OUTPUTS:

  f1 : face number opposite f0

RETURNS: Error Code

*/
  int FHex[6] = {5,3,4,1,2,0};

  (*f1) = -1;

  if (Basis == xfe_QuadLagrange){
    (*f1) = (f0+2)%4;
  }
  else if (Basis == xfe_HexLagrange){
    (*f1) = FHex[f0];
  }
  else return xf_Error(xf_NOT_SUPPORTED);

  return xf_OK;
}



/******************************************************************/
//  FUNCTION Definition: xf_NormalizeVector
static void
xf_NormalizeVector(int n, real *x)
{
  int i;
  real xx;

  for (i=0, xx=0.; i<n; i++) xx += x[i]*x[i];
  xx = sqrt(xx);
  for (i=0; i<n; i++) x[i] /= xx;

}


/******************************************************************/
//  FUNCTION Definition: xf_PrintVector
static void
xf_PrintVector(int n, real *x)
{
  int i;

  for (i=0; i<n; i++) xf_printf("%.3E ", x[i]);
  xf_printf("\n");

}

/******************************************************************/
//  FUNCTION Definition: xf_OrthoFaces
static int 
xf_OrthoFaces(enum xfe_BasisType Basis, int *F, int *OF)
{
/*

PURPOSE: 

  Determines two faces (opposite each other) orthogonal to the (dim-1)
  faces in F.  The faces in F must not be opposite each other.
  
INPUTS:

  QBasis : element basis type
  F : (dim-1) faces (not opposite each other)

OUTPUTS:

  OF : face pair orthogonal to F

RETURNS: Error Code

*/
  int ierr;
  int nface, dim, i;
  int nleft, f;
  int M[6];
  int FHex[6] = {5,3,4,1,2,0};

  nface  = ((dim == 2) ?     4 :        6);   // # faces per elem
  dim = ((Basis == xfe_QuadLagrange) ? 2 : 3);

  for (i=0; i<nface; i++) M[i] = 0;
  for (i=0; i<(dim-1); i++){
    ierr = xf_Error(xf_OppositeFace(Basis, F[i], &f));
    if (ierr != xf_OK) return ierr;
    M[F[i]] = M[f] = 1;
  }

  nleft=0;
  for (i=0; i<nface; i++)
    if (M[i] == 0) OF[nleft++] = i;

  if (nleft != 2) return xf_Error(xf_INPUT_ERROR);
    
  return xf_OK;
}




/******************************************************************/
//  FUNCTION Definition: xf_StartingElem
static int 
xf_StartingElem(xf_Mesh *Mesh, int *ElemFlag, int *estart, int *BF)
{
/*

PURPOSE: 

  Calculates the starting element and direction for a high q agglomeration traversal
  
INPUTS:

  Mesh : mesh structure
  ElemFlag : integer flag over elements, for working storage only

OUTPUTS:

  estart : starting element
  BF : vector of 3 indices indicating local faces to start with
       opposite BF[0] will be the fastest running index (i)
       opposite BF[1] will be the medium  running index (j)
       opposite BF[2] will be the slowest running index (k)

RETURNS: Error Code

*/
  int ierr, dim, i, j;
  int nface, face;
  int elem, ie;
  int nbfgrp, ibfgrp, ibface;
  int e0, e1, ne, egR, f1;
  int nBF, nBN, found;
  int in = 0;
  int OF[2], BN[2];
  int *elist;
  real dpmin, dp0, dp1, dpmin0, dpmin1;
  real n0[2][3], n1[2][3];
  real x0[3] = {0,0,0};
  xf_Face Face;
  xf_BFace BFace;
  xf_ElemGroup *EG;
  xf_NormalData *NData = NULL;
  
  dim    = Mesh->Dim;
  nbfgrp = Mesh->nBFaceGroup;
  EG = Mesh->ElemGroup + 0;

  nface  = ((dim == 2) ?     4 :        6);   // # faces per elem

  // set ElemFlag = # bfaces on each elem
  for (elem=0; elem<EG->nElem; elem++) ElemFlag[elem] = 0;
  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++)
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++)
      ElemFlag[Mesh->BFaceGroup[ibfgrp].BFace[ibface].Elem]++;
  
  // starting elem = one with dim adjacent boundaries
  for (elem=0, e0=-1; elem<EG->nElem; elem++)
    if (ElemFlag[elem] == dim){
      e0 = elem; 
      break;
    }

  if (e0 < 0){ // no such starting element (corner) exists

    // this is still ok; pick elements with dim-1 bfaces
    for (elem=0, ne=0; elem<EG->nElem; elem++)
      if (ElemFlag[elem] == dim-1) ne++;

    if (ne <= 0) return xf_Error(xf_INPUT_ERROR); // mesh not correct

    ierr = xf_Error(xf_Alloc( (void **) &elist, ne, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    for (elem=0, ne=0; elem<EG->nElem; elem++)
      if (ElemFlag[elem] == dim-1) elist[ne++] = elem;

    // loop through elist, find elem whose neighbors have one fewer bfaces
    found = 0;
    for (ie=0; ie<ne; ie++){
      elem = elist[ie];

      // determine (dim-1) local bfaces BF[0:(dim-2)]
      for (i=0; i<dim; i++) BF[i] = -1;
      for (face=0, nBF=0; face<nface; face++)
	if (EG->Face[elem][face].Group >= 0)
	  BF[nBF++] = face;
      if (nBF != dim-1) return xf_Error(xf_INPUT_ERROR);
      
      // determine two orthogonal faces to BF, OF[0], OF[1]
      ierr = xf_Error(xf_OrthoFaces(EG->QBasis, BF, OF));
      if (ierr != xf_OK) return ierr;

      // if ortho neighbor makes a sharp turn, done
      for (in=0; in<2; in++){
	// neighbor element, e1
	ierr = xf_Error(xf_NeighborAcrossFace(Mesh, 0, elem, OF[in], &egR, &e1, &f1));
	if (ierr != xf_OK) return ierr;

	if (e1 < 0)  return xf_Error(xf_INPUT_ERROR); 

	for (face=0, nBN=0; face<nface; face++)
	  if (EG->Face[e1][face].Group >= 0) BN[nBN++] = face;
	if (nBN == dim-2){  // neighbor has one fewer bface, done
	  found = 1;
	  break;
	}

	if (nBN != dim-1) return xf_Error(xf_MESH_ERROR);

	// determine normals
	for (i=0; i<dim-1; i++){
	  // elem
	  Face = EG->Face[elem][BF[i]];
	  BFace = Mesh->BFaceGroup[Face.Group].BFace[Face.Number];
	  ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, 1, x0, &NData, NULL));
	  if (ierr != xf_OK) return ierr;
	  xf_NormalizeVector(3, NData->n);
	  for (j=0; j<dim; j++) n0[i][j] = NData->n[j];
	  // e1
	  Face = EG->Face[e1][BN[i]];
	  BFace = Mesh->BFaceGroup[Face.Group].BFace[Face.Number];
	  ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, 1, x0, &NData, NULL));
	  if (ierr != xf_OK) return ierr;
	  xf_NormalizeVector(3, NData->n);
	  for (j=0; j<dim; j++) n1[i][j] = NData->n[j];
	} // i
	
	if (dim == 2){
	  dpmin = n0[0][0]*n1[0][0] + n0[0][1]*n1[0][1];
	}
	else{
	  for (i=0, dp0=0.; i<dim; i++) dp0 += n0[0][i]*n1[0][i];
	  for (i=0, dp1=0.; i<dim; i++) dp1 += n0[1][i]*n1[1][i];

	  dpmin0 = min(dp0,dp1);
	  for (i=0, dp0=0.; i<dim; i++) dp0 += n0[0][i]*n1[1][i];
	  for (i=0, dp1=0.; i<dim; i++) dp1 += n0[1][i]*n1[0][i];

	  dpmin1 = min(dp0,dp1);
	  dpmin = max(dpmin0, dpmin1);
	}
       
	
	// did we hit a corner?
	if (dpmin < DP_THRESH){
	  found = 1;
	  break;
	}

      } // in

      if (found) break;
    } // ie
    
    if (!found) return xf_Error(xf_INPUT_ERROR); // mesh not correct
    // note, can arbitrarily choose one of the elements in elist here
    
    // corner directions are BF[0], BF[1], OF[in]
    BF[2] = OF[in];
    e0 = elem;

    xf_Release( (void *) elist);
    
  } // if corner does not exist
  else{
    // identify dim boundaries on ecorner
    for (i=0; i<dim; i++) BF[i] = -1;
    for (face=0, nBF=0; face<nface; face++)
      if (EG->Face[e0][face].Group >= 0)
	BF[nBF++ ] = face;
    if (nBF != dim) return xf_Error(xf_INPUT_ERROR); // ecorner not a corner
  }

  (*estart) = e0;

  /* Destroy Normal Data */
  ierr = xf_Error(xf_DestroyNormalData(NData));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}



/******************************************************************/
//  FUNCTION Definition: xf_StartingFront
static int 
xf_StartingFront(xf_Mesh *Mesh, int Q, int ecorner, int *BF, int *Front)
{
/*

PURPOSE: 

  Calculates the starting Front (elems and faces) based on a corner
  elem.
  
INPUTS:

  Mesh : mesh structure
  Q : desired geom order; Front will have Q^(dim-1) elems
  ecorner : corner elem
  BF : directions to take from corner

OUTPUTS:

  Front : elem,face pairs defining the front

RETURNS: Error Code

*/
  int ierr;
  int dim, i, j, k;
  int face, nface;
  int nlayer;
  int e0, e1, egR, eR, fR;
  int f0, f1;
  int ibfgrp;
  int match;
  int iQ, iQ0, iQ1;
  int T[3];
  // corner to face maps, sorted (HS) and following required right-hand rule (HC)
  int HS[8][3] = {{0,1,2},{0,1,4},{0,2,3},{0,3,4},{1,2,5},{1,4,5},{2,3,5},{3,4,5}};
  int HC[8][3] = {{0,1,2},{0,4,1},{0,2,3},{0,3,4},{1,5,2},{4,5,1},{2,5,3},{3,5,4}};
  xf_ElemGroup *EG;

  dim = Mesh->Dim;
  EG = Mesh->ElemGroup + 0;

  nface  = ((dim == 2) ?     4 :        6);   // # faces per elem

  // sort BF so that right-hand rule positivity is maintained
  if (dim == 2){
    if ( ((BF[0]+1)%4) != BF[1]) swap(BF[0], BF[1], i);
    if ( ((BF[0]+1)%4) != BF[1]) return xf_Error(xf_INPUT_ERROR);
  }
  else{
    
    // sort BF
    for (i=0; i<3; i++)
      for (j=0; j<2; j++)
	if (BF[j] > BF[j+1]) swap(BF[j], BF[j+1], k);

    // match to known cases
    match = 0;
    for (i=0; i<8; i++){
      for (j=0, match=1; j<3; j++)
	if (BF[j] != HS[i][j]) match =0;
      if (match) break;
    }
    if (!match) return xf_Error(xf_INPUT_ERROR);
    for (j=0; j<3; j++) BF[j] = HC[i][j];
  }

  // ensure that BF[1] and BF[2] are boundary faces (BF[0] can be interior)
  for (i=0; i<dim; i++) 
    if (EG->Face[ecorner][BF[i]].Group < 0){
      for (j=0; j<dim; j++) T[j] = BF[(i+j)%dim];
      for (j=0; j<dim; j++) BF[j] = T[j];
      break;
    }  
  for (i=1; i<dim; i++)
    if (EG->Face[ecorner][BF[i]].Group < 0) return xf_Error(xf_INPUT_ERROR);

  // f0 = first local face
  f0 = BF[0];

  // traverse Q-1 times opposite f0 from ecorner -> first Q elems in Front
  Front[0] = e0 = ecorner;
  for (iQ=1; iQ<Q; iQ++){
    ierr = xf_Error(xf_OppositeFace(EG->QBasis, f0, &f1));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_NeighborAcrossFace(Mesh, 0, e0, f1, &egR, &eR, &fR));
    if (ierr != xf_OK) return ierr;
    if (eR < 0) return xf_Error(xf_MESH_ERROR);
    Front[2*iQ] = eR;
    e0 = eR;
    f0 = fR;
  } // iQ
  
  // if 3d, extend layer orthogonal to BF[1]
  if (dim == 3){
    ibfgrp = EG->Face[ecorner][BF[1]].Group; // bface[BF[1]]
    for (iQ0=0; iQ0<Q; iQ0++){
      e0 = Front[2*iQ0]; // starting elem for extension
      // f0 = local face on e0 adjacent to ibfgrp
      for (face=0, f0=-1; face<nface; face++)
	if (EG->Face[e0][face].Group == ibfgrp) f0 = face;
      if (f0 < 0) return xf_Error(xf_MESH_ERROR);
      for (iQ1=1; iQ1<Q; iQ1++){
	ierr = xf_Error(xf_OppositeFace(EG->QBasis, f0, &f1));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_NeighborAcrossFace(Mesh, 0, e0, f1, &egR, &eR, &fR));
	if (ierr != xf_OK) return ierr;
	if (eR < 0) return xf_Error(xf_MESH_ERROR);
	Front[2*(iQ0+Q*iQ1)] = eR;
	e0 = eR;
	f0 = fR;
      } // iQ1
    } // iQ0
  }

  // use bface(BF[dim-1]) to define layer direction (local faces)
  nlayer = ((dim == 2) ?     Q :      Q*Q);
  ibfgrp = EG->Face[ecorner][BF[dim-1]].Group;
  for (i=0; i<nlayer; i++){
    e0 = Front[2*i]; // ith elem in layer
    // f0 = local face on e0 adjacent to ibfgrp
    for (face=0, f0=-1; face<nface; face++)
      if (EG->Face[e0][face].Group == ibfgrp) f0 = face;
    if (f0 < 0) return xf_Error(xf_MESH_ERROR);
    // opposite f0 is the face we want in the front
    ierr = xf_Error(xf_OppositeFace(EG->QBasis, f0, Front + 2*i+1));
    if (ierr != xf_OK) return ierr;
  } // i

  return xf_OK;
}


/******************************************************************/
//  FUNCTION Definition: xf_NextLayer
static int 
xf_NextLayer(xf_Mesh *Mesh, int Q, const int *Front, int *Next)
{
/*

PURPOSE: 

  Calculates the next layer of elements and faces based on a Front
  (elems and faces).
  
INPUTS:

  Mesh : mesh structure
  Q : desired geometry order
  Front : elem,face pairs defining the front

OUTPUTS:
  Next : elem,face pairs defining next front, adjacent to the input Front

RETURNS: Error Code

*/
  int ierr, dim;
  int i, egR, f1;
  int nlayer;
  enum xfe_BasisType QBasis;

  dim = Mesh->Dim;
  QBasis = Mesh->ElemGroup[0].QBasis;

  nlayer = ((dim == 2) ?     Q :      Q*Q);   // # elems per layer

  for (i=0; i<nlayer; i++){
    ierr = xf_Error(xf_NeighborAcrossFace(Mesh, 0, Front[2*i+0], Front[2*i+1], 
					  &egR, Next+2*i, &f1));
    if (ierr != xf_OK) return ierr;
    if (egR < 0){
      xf_printf("\nHit boundary when calculating next layer (e=%d, f=%d).\n",
		Front[2*i+0], Front[2*i+1]);
      return xf_Error(xf_OUT_OF_BOUNDS);
    }
    // get opposite face to establish front
    ierr = xf_Error(xf_OppositeFace(QBasis, f1, Next+2*i+1));
    if (ierr != xf_OK) return ierr;

  }

  
  return xf_OK;
}


/******************************************************************/
//  FUNCTION Definition: xf_SubList2Node
static int 
xf_SubList2Node(xf_Mesh *Mesh, int Q, const int *SubList, int *NewNode)
{
/*

PURPOSE: 

  Calculates the q=Q nodes associated with a q=1 subelement list.
  
INPUTS:

  Mesh : mesh structure
  Q : desired geom order
  SubList : ordered list of elements (Q^dim)

OUTPUTS:  

  Node : list of q=Q nodes

RETURNS: Error Code

*/

  int dim, i, j, k, nsub;
  int inode, elem, nnode;
  int ie, je, ke, in0;
  int nadj, cnt;
  int found;
  int Qp, neQ, nnQ1;
  int *Node;
  int plist[2*(MQ+1)*(MQ+1)*(MQ+1)];
  int alist[2*(MQ+1)*(MQ+1)*(MQ+1)][8];
  int II[8] = {0,1,0,1,0,1,0,1};
  int JJ[8] = {0,0,1,1,0,0,1,1};
  int KK[8] = {0,0,0,0,1,1,1,1};
  static int alist0_Q = 0;
  static int alist0[2*(MQ+1)*(MQ+1)*(MQ+1)][9];

  dim = Mesh->Dim;
  Qp     = Q + 1;
  neQ    = ((dim == 2) ? Qp*Qp : Qp*Qp*Qp);   // # q=Q nodes per elem
  nsub   = ((dim == 2) ?   Q*Q :    Q*Q*Q);   // # subelements
  nnQ1   = ((dim == 2) ?     4 :        8);   // # q=1 nodes per elem

  // (re)build reference adjacency list if necessary
  if (alist0_Q != Q){ 
    for (elem=0; elem<nsub; elem++){
      // element indices
      ie =  elem%Q;
      je = (elem/Q)%Q;
      ke = elem/(Q*Q);
      in0 = ke*Qp*Qp + je*Qp + ie;
      for (i=0; i<nnQ1; i++){
	inode = in0 + KK[i]*Qp*Qp + JJ[i]*Qp + II[i];
	nadj = alist0[inode][0]++;
	alist0[inode][nadj+1] = elem;
      } // i
    } // elem
    alist0_Q = Q;

  }

  // build list of global nodes and adjacent elements
  for (i=0; i<neQ; i++){
    plist[2*i+0] = -1;   // glob node number
    plist[2*i+1] = 0;    // adjacent elem counter
  }
  for (i=0; i<neQ; i++)
    for (j=0; j<8; j++) 
      alist[i][j] = -1;  // adjacent element list

  nnode = 0;
  for (elem=0; elem<nsub; elem++){
    Node = Mesh->ElemGroup[0].Node[SubList[elem]];
    for (inode=0; inode<nnQ1; inode++){
      k = -1;
      for (i=0; i<nnode; i++)
	if (plist[2*i] == Node[inode]){
	  k = i;
	  break;
	}
      if (k == -1) k = nnode++;
      plist[2*k+0] = Node[inode];
      nadj = plist[2*k+1]++;
      alist[k][nadj] = elem;
    } // inode
  } // elem

  if (nnode != neQ) return xf_Error(xf_OUT_OF_BOUNDS);


  // Determine new nodes by comparing element adjacencies  
  for (inode=0; inode<nnode; inode++){
    k = -1;
    for (i=0; i<nnode; i++)
      if (plist[2*i+1] == (cnt=alist0[inode][0])){ // counters match
	for (j=0, found=1; j<cnt; j++)
	  if (alist[i][j] != alist0[inode][1+j]) found = 0;
	if (found){ // element adjacencies match
	  k = i;
	  break;
	}
      }
    if (k < 0) return xf_Error(xf_NOT_FOUND);
    NewNode[inode] = plist[2*k]; // store global node number of inode
  } // inode
  
  return xf_OK;
}


/******************************************************************/
//  FUNCTION Definition: xf_MoveEdgeDist2D
static int 
xf_MoveEdgeDist2D(int Q, int dim, real **x, real **xout)
{
/*

PURPOSE: 

  Moves 2D-face nodes (from 3D elements) to high order positions,
  using existing (Q+1)x(Q+1) node positions.  Only interior nodes are
  moved; edges are left alone.
  
INPUTS:

  Q : geometry order
  x    : x[i][d] = coordinate d of node i
         0 <= d < dim,  0 <= i < (Q+1)*(Q+1)

OUTPUTS:  

  xout : moved nodes

RETURNS: Error code

*/
  int ierr;
  int i, j, k, n, Qp;
  int ij, d;
  int P[(MQ+1)*(MQ+1)];
  real dist;
  real XD[(MQ+1)*(MQ+1)];
  real YD[(MQ+1)*(MQ+1)];
  real X[(MQ+1)*(MQ+1)];
  real xnew[(MQ+1)*(MQ+1)*3];
  xf_BasisData *PhiData = NULL;

  Qp = Q+1;

  // evaluate reference coords (XD,YD) of every node, using distances
  for (j=0; j<Q+1; j++){ // XD
    n = j*(Q+1);
    XD[n] = 0.;
    for (i=1; i<Q+1; i++){
      for (d=0, dist=0; d<dim; d++)
	dist +=(x[n+i][d]-x[n+i-1][d])*(x[n+i][d]-x[n+i-1][d]);
      XD[n+i] = XD[n+i-1] + sqrt(dist);
    }
    dist = XD[n+Q];
    for (i=0; i<Q+1; i++) XD[n+i] /= dist;
  }
/*   // TEMPORARY */
/*   xf_printf("XD = [\n"); */
/*   for (j=0; j<Q+1; j++){ */
/*     for (i=0; i<Q+1; i++) */
/*       xf_printf("%.12f ", XD[j*(Q+1)+i]); */
/*     xf_printf("\n"); */
/*   } */
/*   xf_printf("]\n"); */
      
  for (i=0; i<Q+1; i++){ // YD
    YD[i] = 0.;
    for (j=1; j<Q+1; j++){
      for (d=0, dist=0; d<dim; d++)
	dist +=(x[Qp*j+i][d]-x[Qp*(j-1)+i][d])*(x[Qp*j+i][d]-x[Qp*(j-1)+i][d]);
      YD[Qp*j+i] = YD[Qp*(j-1)+i] + sqrt(dist);
    }
    dist = YD[Qp*Q+i];
    for (j=0; j<Q+1; j++) YD[Qp*j+i] /= dist;
  }
/*   // TEMPORARY */
/*   xf_printf("YD = [\n"); */
/*   for (j=0; j<Q+1; j++){ */
/*     for (i=0; i<Q+1; i++) */
/*       xf_printf("%.12f ", YD[j*(Q+1)+i]); */
/*     xf_printf("\n"); */
/*   } */
/*   xf_printf("]\n"); */

  // lump XD and YD into one X
  for (k=0; k<Qp*Qp; k++){
    X[2*k+0] = XD[k];
    X[2*k+1] = YD[k];
  }

  // evaluate tensor-product Lagrange basis at nodes -> build matrix
  ierr = xf_Error(xf_EvalBasis(xfe_QuadLagrange, Q, xfe_True, Qp*Qp, X, 
			       xfb_Phi, &PhiData));
  if (ierr != xf_OK) return ierr;

  // Solve (Phi)*xnew = x   (dim rhs vectors)
  ierr = xf_Error(xf_ComputePLU(PhiData->Phi, Qp*Qp, P));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<Qp*Qp; k++) 
    for (d=0; d<dim; d++) xnew[dim*k+d] = x[k][d];
  ierr = xf_Error(xf_SolvePLU_Matrix(PhiData->Phi, P, Qp*Qp, dim, xnew));
  if (ierr != xf_OK) return ierr;

  // set interior nodes
  for (j=1; j<Q; j++)
    for (i=1; i<Q; i++){
      ij = j*(Q+1) + i;
      for (d=0; d<dim; d++) xout[ij][d] = xnew[dim*ij+d];
    }

  // destroy BasisData
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;


  return xf_OK;
}


/******************************************************************/
//  FUNCTION Definition: xf_MoveBFaceNodes
static int 
xf_MoveFaceNodes(xf_Mesh *Mesh, int Q,const int *NewNode, int *NodeFlag,
                 real **NewCoord)
{
/*

PURPOSE: 

  Moves face nodes in 3D, using interior point info
  
INPUTS:

  Mesh : mesh structure
  Q : geometry order
  NewNode : list of q=Q nodes
  NodeFlag : list of flags indicating which nodes have been moved

OUTPUTS:  

RETURNS: Error Code

*/
  int ierr, dim, k;
  int i, j, ij;
  int nface, iface;
  enum xfe_Bool moved;
  int FN[(MQ+1)*(MQ+1)];
  real *xf[(MQ+1)*(MQ+1)*(MQ+1)];
  real *xfnew[(MQ+1)*(MQ+1)*(MQ+1)];
  enum xfe_BasisType QBasis;

  QBasis = Mesh->ElemGroup[0].QBasis;

  if (Mesh->Dim != 3) return xf_Error(xf_INPUT_ERROR); 
  
  nface = 6;

  for (iface=0; iface<nface; iface++){
    
    ierr = xf_Error(xf_NodesOnFace(QBasis, Q, iface, &k, FN));
    if (ierr != xf_OK) return ierr;
    if (k != (Q+1)*(Q+1)) return xf_Error(xf_OUT_OF_BOUNDS);
    
    // check if moved; if so, continue
    moved = xfe_False;
    for (j=1; j<Q; j++)    
      for (i=1; i<Q; i++){
	ij = j*(Q+1) + i;
	if (NodeFlag[NewNode[FN[ij]]] != -1) moved = xfe_True;
      }
    if (moved) continue;

    /* use interior point info in face movement */

    // first, copy over node coord pointers into xf
    for (k=0; k<(Q+1)*(Q+1); k++) xf[k] = Mesh->Coord[NewNode[FN[k]]];
    for (k=0; k<(Q+1)*(Q+1); k++) xfnew[k] = NewCoord[NewNode[FN[k]]];
    // next, move nodes using interior point info (only interior nodes moved)
    ierr = xf_Error(xf_MoveEdgeDist2D(Q, 3, xf, xfnew));
    if (ierr != xf_OK) return ierr;
    // mark interior nodes as moved (should not have been moved before)
    for (j=1; j<Q; j++)    
      for (i=1; i<Q; i++){
	ij = j*(Q+1) + i;
	if (NodeFlag[NewNode[FN[ij]]] != -1) return xf_Error(xf_CODE_LOGIC_ERROR);
	NodeFlag[NewNode[FN[ij]]] = 1;
      }

  } // iface
  return xf_OK;
}

/******************************************************************/
//  FUNCTION Definition: xf_GetDistance
static real
xf_GetDistance(int dim, real *a, real *b, real *dist){
  int d;
  (*dist) = 0;
  for (d=0; d<dim; d++) (*dist) += (a[d]-b[d])*(a[d]-b[d]);
  (*dist) = sqrt((*dist));
}

/******************************************************************/
//  FUNCTION Definition: xf_ProjectToLine
static real
xf_ProjectToLine(int dim, real *a, real *b, real *p, real *q){
  // p is projected onto the line through a,b; result is point q
  // idea: q = a + s*(b-a);  (q-p) dot (b-a) = 0;
  int d;
  real num, den, s;

  for (d=0, num=0.; d<dim; d++) num += (p[d]-a[d])*(b[d]-a[d]);
  for (d=0, den=0.; d<dim; d++) den += (b[d]-a[d])*(b[d]-a[d]);

  s = num/den;

  for (d=0; d<dim; d++) q[d] = a[d] + s*(b[d]-a[d]);
}


/******************************************************************/
//  FUNCTION Definition: xf_MoveNodes
static int 
xf_MoveNodes(xf_Mesh *Mesh, int Q, const int *NewNode, int *NodeFlag,
             real **NewCoord, real SnapFactor)
{
/*

PURPOSE: 

  Moves unmoved nodes into parametric mapping position
  
INPUTS:

  Mesh : mesh structure
  Q : geometry order
  NewNode : list of q=Q nodes
  NodeFlag : list of flags indicating which nodes have been moved
  SnapFactor : if > 0, dictates tolerance below which curved edges will
               be snapped to linear

OUTPUTS:  

RETURNS: Error Code

*/
  int ierr, dim, d, dn=0;
  int i, j, k, l;
  int ii, jj, kk;
  int ij, ijk;
  int i0, Qp2, id;
  int iedge, nedge;
  int iface, nface;
  int EN[MQ+1];
  int ivec[3], jvec[3], kvec[3];
  int FN[(MQ+1)*(MQ+1)];
  enum xfe_Bool moved, CanSnap = xfe_False;
  real dist, s, svec[MQ+1];
  real *xe[MQ+1], *xf[(MQ+1)*(MQ+1)*(MQ+1)];
  real *xfnew[(MQ+1)*(MQ+1)*(MQ+1)];
  real xnew[MQ+1][3];
  real phi[MQ+1];
  real xref[3] = {0};
  real xproj[3], delta, maxdelta;
  enum xfe_BasisType QBasis;

  dim = Mesh->Dim;
  QBasis = Mesh->ElemGroup[0].QBasis;

  nedge = ((dim == 2) ? 4 : 12);
  
  /*-----------------*/
  /* Move edge nodes */
  /*-----------------*/

  // Move edge nodes
  for (iedge=0; iedge<nedge; iedge++){

    // pull off edge nodes
    if (dim == 2){
      switch (iedge){
      case 0: i0 = 0;       d =    1; dn =  Q*(Q+1); break;
      case 1: i0 = Q;       d =  Q+1; dn =       -Q; break;
      case 2: i0 = Q*(Q+2); d =   -1; dn = -Q*(Q+1); break;
      case 3: i0 = Q*(Q+1); d = -Q-1; dn =        Q; break;
      }
    }
    else{
      Qp2 = (Q+1)*(Q+1); // vertical layer diff
      id  =  Qp2*Q;      // diff from bottom to top
      switch (iedge){
	// edges on bottom
      case  0: i0 = 0;          d =    1; break;
      case  1: i0 = Q;          d =  Q+1; break;
      case  2: i0 = Q*(Q+2);    d =   -1; break;
      case  3: i0 = Q*(Q+1);    d = -Q-1; break;
	// edges on sides
      case  4: i0 = 0;          d =  Qp2; break;
      case  5: i0 = Q;          d =  Qp2; break;
      case  6: i0 = Q*(Q+2);    d =  Qp2; break;
      case  7: i0 = Q*(Q+1);    d =  Qp2; break;
	// edges on top
      case  8: i0 = 0      +id; d =    1; break;
      case  9: i0 = Q      +id; d =  Q+1; break;
      case 10: i0 = Q*(Q+2)+id; d =   -1; break;
      case 11: i0 = Q*(Q+1)+id; d = -Q-1; break;
      }
      CanSnap = xfe_False; // for now
    }
    for (k=0, i=i0; k<Q+1; k++, i+=d) EN[k] = i;

    if (NodeFlag[NewNode[EN[1]]] != -1) continue; // edge nodes already moved

    for (k=0; k<Q+1; k++) // edge point coords
      xe[k] = NewCoord[NewNode[EN[k]]];

    // check if we can snap edge to linear
    if ((dim == 2) && (SnapFactor > 0.)){
      CanSnap = xfe_True;
      for (k=1, i=i0+d; ((k<Q) && (CanSnap)); k++, i+=d){
	xf_ProjectToLine(dim, xe[0], xe[Q], xe[k], xproj);
	xf_GetDistance(dim, xe[k], xproj, &delta);
	xf_GetDistance(dim, xe[k], Mesh->Coord[NewNode[i+dn]], &maxdelta);
	if (delta > SnapFactor*maxdelta) CanSnap = xfe_False;
      }
      if (CanSnap){ // perform the snap
	for (k=1; k<Q; k++){
	  xf_ProjectToLine(dim, xe[0], xe[Q], xe[k], xproj);
	  for (i=0;i<dim;i++) xe[k][i] = xproj[i];
	}
      }
    }
    
    // ref pos vec determined by inter-node distances
    svec[0] = 0.0;
    for (k=0; k<Q; k++){
      for (d=0, dist=0.; d<dim; d++)
	dist += (xe[k+1][d]-xe[k][d])*(xe[k+1][d]-xe[k][d]);
      svec[k+1] = svec[k]+sqrt(dist);
    }

    // normalize svec
    for (k=1; k<Q; k++) svec[k] /= svec[Q];
    svec[Q] = 1.0;

    // loop over middle nodes and interpolate new positions
    for (i=1; i<Q; i++){
      
      // error if this node was already moved
      if (NodeFlag[NewNode[EN[i]]] != -1) return xf_Error(xf_CODE_LOGIC_ERROR);
      NodeFlag[NewNode[EN[i]]] = 1; // flag as moved
      
      s = ((real) i) / ((real) Q); // ref position

      // obtain 1D Lagrange basis fcns at s, based on svec nodes
      ierr = xf_Error(xf_BasisLagrange1D(s, svec, Q+1, phi, NULL, NULL));
      if (ierr != xf_OK) return ierr;

      // interpolate new node position
      for (d=0; d<dim; d++)
	for (k=0, xnew[i][d]=0.; k<Q+1; k++)
	  xnew[i][d] += phi[k]*xe[k][d];
      
    } // i

    // set new coords
    for (i=1; i<Q; i++)
      for (d=0; d<dim; d++)
	xe[i][d] = xnew[i][d];

  } // iedge


  nface = ((dim == 2) ? 1 : 6);
  
  /*-------------------------------------------------*/
  /* Move face nodes (not transfinite interpolation) */
  /*-------------------------------------------------*/

  for (iface=0; iface<nface; iface++){

    if (dim == 3){
      ierr = xf_Error(xf_NodesOnFace(QBasis, Q, iface, &k, FN));
      if (ierr != xf_OK) return ierr;
      if (k != (Q+1)*(Q+1)) return xf_Error(xf_OUT_OF_BOUNDS);

      // check if moved; if so, continue;
      moved = xfe_False;
      for (j=1; j<Q; j++)    
	for (i=1; i<Q; i++){
	  ij = j*(Q+1) + i;
	  if (NodeFlag[NewNode[FN[ij]]] != -1) moved = xfe_True;
	}
      if (moved) continue;

      // first, copy over node coord pointers into xf
      for (k=0; k<(Q+1)*(Q+1); k++) xf[k] = Mesh->Coord[NewNode[FN[k]]];
      for (k=0; k<(Q+1)*(Q+1); k++) xfnew[k] = NewCoord[NewNode[FN[k]]];
      // next, move nodes using interior point info (only interior nodes moved)
      ierr = xf_Error(xf_MoveEdgeDist2D(Q, 3, xf, xfnew)); 
      if (ierr != xf_OK) return ierr;
      // mark interior nodes as moved (should not have been moved before)
      for (j=1; j<Q; j++)    
	for (i=1; i<Q; i++){
	  ij = j*(Q+1) + i;
	  if (NodeFlag[NewNode[FN[ij]]] != -1) return xf_Error(xf_CODE_LOGIC_ERROR);
	  NodeFlag[NewNode[FN[ij]]] = 1;
	}
      continue;
    }
    else{
      for (k=0; k<(Q+1)*(Q+1); k++) FN[k] = k;
      // do not do transfinite interpolation even in 2D; instead ...
      
      // copy over node coord pointers into xf
      for (k=0; k<(Q+1)*(Q+1); k++) xf[k] = Mesh->Coord[NewNode[FN[k]]];
      for (k=0; k<(Q+1)*(Q+1); k++) xfnew[k] = NewCoord[NewNode[FN[k]]];
      // next, move nodes using interior point info (only interior nodes moved)
      ierr = xf_Error(xf_MoveEdgeDist2D(Q, 2, xf, xfnew)); 
      if (ierr != xf_OK) return ierr;
      // mark interior nodes as moved (should not have been moved before)
      for (j=1; j<Q; j++)    
	for (i=1; i<Q; i++){
	  ij = j*(Q+1) + i;
	  if (NodeFlag[NewNode[FN[ij]]] != -1) return xf_Error(xf_CODE_LOGIC_ERROR);
	  NodeFlag[NewNode[FN[ij]]] = 1;
	}
      continue;
    }

    // loop over interior points
    for (j=1; j<Q; j++)    
      for (i=1; i<Q; i++){

	// if face node was already moved, continue
	ij = j*(Q+1) + i;
	if (NodeFlag[NewNode[FN[ij]]] != -1) continue;

        // should never get here (not using trasfinite 2d interpolation any more)
        return xf_Error(xf_CODE_LOGIC_ERROR);

	// ref coords (equispaced)
	xref[0] = ((real) i) / ((real) Q);
	xref[1] = ((real) j) / ((real) Q);

	jvec[0]=0; jvec[1]=j; jvec[2]=Q;
	ivec[0]=0; ivec[1]=i; ivec[2]=Q;

	// get coord list for transfinite interpolation
	k = 0;
	for (jj=0; jj<3; jj++)
	  for (ii=0; ii<3; ii++){
	    ij = jvec[jj]*(Q+1) + ivec[ii];
	    xfnew[k++] = NewCoord[NewNode[FN[ij]]];
	  }

	// move interior point
	xf_Transfinite2D(dim, xref, xfnew);

	// flag interior node as moved
	ij = j*(Q+1) + i;
	NodeFlag[NewNode[FN[ij]]] = 1;

      } // i
  } // iface


  /*----------------------------------------------------*/
  /* (3D) Move volume nodes (transfinite interpolation) */
  /*----------------------------------------------------*/
 
  if (dim == 3){

    // loop over interior points
    for (k=1; k<Q; k++)
      for (j=1; j<Q; j++)    
	for (i=1; i<Q; i++){

	  // ref coords (equispaced)
	  xref[0] = ((real) i) / ((real) Q);
	  xref[1] = ((real) j) / ((real) Q);
	  xref[2] = ((real) k) / ((real) Q);

	  kvec[0]=0; kvec[1]=k; kvec[2]=Q;
	  jvec[0]=0; jvec[1]=j; jvec[2]=Q;
	  ivec[0]=0; ivec[1]=i; ivec[2]=Q;

	  // get coord list for transfinite interpolation
	  l = 0;
	  for (kk=0; kk<3; kk++)
	    for (jj=0; jj<3; jj++)
	      for (ii=0; ii<3; ii++){
		ijk = kvec[kk]*(Q+1)*(Q+1) + jvec[jj]*(Q+1) + ivec[ii];
		xf[l++] = NewCoord[NewNode[ijk]];
	      }
	  
	  // move interior point
	  xf_Transfinite3D(dim, xref, xf);

	  // flag interior node as moved
	  ijk = k*(Q+1)*(Q+1) + j*(Q+1) + i;
	  NodeFlag[NewNode[ijk]] = 1;
	  
	} // i
  } // if dim == 3

  return xf_OK;
}


/******************************************************************/
//  FUNCTION Definition: xf_SubList2FaceLoc
static int 
xf_SubList2FaceLoc(xf_Mesh *Mesh, int Q, const int *SubList, 
		   int faceloc, int *ElemFaceList)
{
/*

PURPOSE: 

  Calculates the subelems and faces on a local face of the HO elem
  
INPUTS:

  Mesh : mesh structure
  Q : desired geom order
  SubList : ordered list of elements (Q^dim)
  faceloc : local face of interest

OUTPUTS:  

  ElemFaceList : elem,face pairs for the faceloc of interest

RETURNS: Error Code

*/
  int ierr, i, j, dim;
  int e0, e1, f0, f1;
  int p, i0, d, d1, dd;
  int itot, jmax;
  enum xfe_BasisType QBasis;
  
  dim = Mesh->Dim;
  QBasis = Mesh->ElemGroup[0].QBasis;

  p = Q-1; // code taken from xf_NodesOnFace; we're interested in elems
  if (dim==2){
    switch (faceloc){
    case 0: i0 = 0;       d =    1; dd =  Q; break;
    case 1: i0 = p;       d =  p+1; dd = -1; break;
    case 2: i0 = p*(p+2); d =   -1; dd = -Q; break;
    case 3: i0 = p*(p+1); d = -p-1; dd =  1; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
    d1 = 0;
  }
  else{
    switch (faceloc){
    case 0: i0 = 0;             d =   p+1 ; d1 =          1 ; dd =  Q*Q; break;
    case 1: i0 = 0;             d =     1 ; d1 = (p+1)*(p+1); dd =    Q; break;
    case 2: i0 = p;             d =   p+1 ; d1 = (p+1)*(p+1); dd =   -1; break;
    case 3: i0 = (p+1)*(p+1)-1; d =    -1 ; d1 = (p+1)*(p+1); dd =   -Q; break;
    case 4: i0 = p*(p+1);       d = -(p+1); d1 = (p+1)*(p+1); dd =    1; break;
    case 5: i0 = p*(p+1)*(p+1); d =     1 ; d1 =        p+1 ; dd = -Q*Q; break;
    default: return xf_Error(xf_OUT_OF_BOUNDS); break;
    }
  }

  itot = ((dim == 2) ? Q : Q*Q); // total # elems on a face
  jmax = ((dim == 2) ? 1 :   Q); // # iters in j for face traversal

  for (j=0; j<jmax; j++)
    for (i=0; i<Q; i++){
      e0 = i0 + i*d + j*d1; // loc subelem on face
      e1 = e0 + dd;         // adjacent loc subelem one layer interior
      if ((e0<0) || (e1<0)) return xf_Error(xf_CODE_LOGIC_ERROR);

      // get actual elem numbers (global)
      e0 = SubList[e0];
      e1 = SubList[e1];

      ierr = xf_Error(xf_CommonFace(Mesh, 0, e0, 0, e1, &f0, &f1));
      if (ierr != xf_OK) return ierr;

      ierr = xf_Error(xf_OppositeFace(QBasis, f0, &f0));
      if (ierr != xf_OK) return ierr;

      itot--;
      ElemFaceList[2*itot+0] = e0;
      ElemFaceList[2*itot+1] = f0;
    } // i

  if (itot != 0) return xf_Error(xf_CODE_LOGIC_ERROR);

  return xf_OK;
}

/******************************************************************/
//  FUNCTION Definition: xf_ConvertHO
static int 
xf_ConvertHO(xf_All *All, int Q, char *outFile, 
	     enum xfe_Bool RemapFlag, enum xfe_Bool UseInt,
	     real SnapFactor)
{
/*

PURPOSE: 

  Converts All (structured) to a HO .gri.
  
INPUTS:

  All : All structure
  Q : desired geometry order
  outFile : name of output .gri file
  RemapFlag : true to use transfinite interpolation to remap interior nodes

OUTPUTS:  (None, outfile is written)

RETURNS: Error Code

*/

  int ierr, i, k, d, dim;
  int flag, egR, eR, faceR;
  int nstack, nlayer, nsub;
  int inode, nnode, nface;
  int faceloc, StackSize;
  int ibfgrp, nbfgrp, ibface;
  int elem, face;
  int e0, f0, e1, f1;
  int iElem, nElem;
  int Qp, neQ, nf1, iQ;
  int ipercent;
  int Front[2*MQ*MQ];
  int Next[2*MQ*MQ];
  int SubList[MQ*MQ*MQ];
  int fvec[4];
  int NewNode[(MQ+1)*(MQ+1)*(MQ+1)];
  int ElemFaceList[2*MQ*MQ];
  int BF[3];
  int *NodeFlag = NULL;
  int *ElemFlag = NULL;
  int *ElemList = NULL;
  int *nBFace = NULL;
  int **ElemStack = NULL;
  int **BFaceList = NULL;
  enum xfe_Bool bflag, DoAgglomerate;
  real **NewCoord = NULL;
  xf_ElemGroup *EG;
  FILE *fgri;
  xf_Mesh *Mesh;

  Mesh   = All->Mesh;
  dim    = Mesh->Dim;
  nbfgrp = Mesh->nBFaceGroup;

  Qp     = Q + 1;
  neQ    = ((dim == 2) ? Qp*Qp : Qp*Qp*Qp);   // # q=Q nodes per element
  nf1    = ((dim == 2) ?     2 :        4);   // # q=1 nodes per face
  nlayer = ((dim == 2) ?     Q :      Q*Q);   // # elems per layer
  nsub   = Q*nlayer;                          // # subelems in a HO elem
  nface  = ((dim == 2) ?     4 :        6);   // # faces per elem

  // input checks
  if ((Mesh->nElemGroup != 1) ||
      ((Mesh->ElemGroup[0].QOrder != 1) && (Mesh->ElemGroup[0].QOrder != Q)) ||
      ((Mesh->ElemGroup[0].QBasis != xfe_QuadLagrange) &&
       (Mesh->ElemGroup[0].QBasis != xfe_HexLagrange)))
    return xf_Error(xf_INPUT_ERROR);
    
  // Q must be between 2 and MQ inclusive
  if ((Q < 2) || (Q > MQ)){
    xf_printf("Q given as %d: but it must be between 2 and %d inclusive.\n", Q, MQ);
    xf_printf("Increase MQ in define if necessary\n");
    return xf_Error(xf_OUT_OF_BOUNDS);
  }

  EG = Mesh->ElemGroup + 0;

  // agglomerate only if order difference in Q
  DoAgglomerate = (Mesh->ElemGroup[0].QOrder != Q);

  if (DoAgglomerate){


    xf_printf("Converting to Q = %d\n", Q);

    if ((EG->nElem % nsub) != 0) return xf_Error(xf_INPUT_ERROR);

    for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++)
      if ((Mesh->BFaceGroup[ibfgrp].nBFace % nlayer) != 0)
	return xf_Error(xf_INPUT_ERROR);

    // sizes for new mesh
    nElem = EG->nElem/nsub;  // number of new elements


    // Allocate memory
    StackSize = EG->nElem;
    ierr = xf_Error(xf_Alloc2( (void ***) &ElemStack, StackSize, 2*nlayer, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &ElemList, nElem*neQ, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &ElemFlag, EG->nElem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &nBFace, nbfgrp, sizeof(int)));
    if (ierr != xf_OK) return ierr;
		  
    for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++)
      nBFace[ibfgrp] = Mesh->BFaceGroup[ibfgrp].nBFace/nlayer * nf1;
  
    ierr = xf_Error(xf_VAlloc2( (void ***) &BFaceList, nbfgrp, nBFace, sizeof(int)));
    if (ierr != xf_OK) return ierr;
  
    // zero out counters
    for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++) nBFace[ibfgrp] = 0;
  

    /*--------------------*/
    /* Find starting elem */
    /*--------------------*/
  
    xf_printf("Finding starting element and front.\n");

    // determine starting element and direction
    ierr = xf_Error(xf_StartingElem(Mesh, ElemFlag, &e0, BF));
    if (ierr != xf_OK) return ierr;

    // determine starting front and push onto stack
    ierr = xf_Error(xf_StartingFront(Mesh, Q, e0, BF, ElemStack[0]));
    if (ierr != xf_OK) return ierr;

    // set ElemFlag to -1
    for (elem=0; elem<EG->nElem; elem++) ElemFlag[elem] = -1;

    /*---------------------------*/
    /* Loop until stack is empty */
    /*---------------------------*/

    xf_printf("Begining conversion loop over elements.\n");

    iElem = 0; // counter for new added elements
    nstack = 1;
    ipercent = 0;
    while (iElem < nElem){

      // track progress
      if ((((real) iElem)/((real) nElem)) > ipercent*.01){
	xf_printf("  %2d%% elements converted\n", ipercent);
	ipercent += 10;
      }

      if (nstack <= 0) return xf_Error(xf_OUT_OF_BOUNDS);

      // pop element front
      nstack--;
      for (i=0; i<2*nlayer; i++) Front[i] = ElemStack[nstack][i];

      // loop over Q layers and assemble rest of subelement list
      for (iQ=0; iQ<Q-1; iQ++){
	// add front elems to subelement list
	for (i=0; i<nlayer; i++) SubList[iQ*nlayer+i] = Front[2*i];
	// calculate next layer
	ierr = xf_Error(xf_NextLayer(Mesh, Q, Front, Next));
	if (ierr != xf_OK) return ierr;
	// set Front == Next
	for (i=0; i<2*nlayer; i++) Front[i] = Next[i];
      } // iQ

      // add last front elems to subelement list
      for (i=0; i<nlayer; i++)
	SubList[(Q-1)*nlayer+i] = Front[2*i];
    

      /* check that the subelement list is not flagged; allow all
	 elements to be flagged, in which case move to next on stack
	 (this means we put the same coarse elem on the stack from
	 different sides) */
      for (i=0, k=0; i<nsub; i++) 
	if (ElemFlag[SubList[i]] >= 0) k++;
      if (k==nsub) continue;
      if (k!=   0) return xf_Error(xf_OUT_OF_BOUNDS);

      // flag elements as taken by iElem
      for (i=0; i<nsub; i++) ElemFlag[SubList[i]] = iElem;
    

      // calculate new elem node list from SubList
      ierr = xf_Error(xf_SubList2Node(Mesh, Q, SubList, NewNode));
      if (ierr != xf_OK) return ierr;


      // loop over faces of macro element, determine which are boundary
      for (faceloc=0; faceloc<nface; faceloc++){
	// calculate elem and face list for faceloc
	ierr = xf_Error(xf_SubList2FaceLoc(Mesh, Q, SubList, faceloc, ElemFaceList));
	if (ierr != xf_OK) return ierr;
	e0 = ElemFaceList[0];
	f0 = ElemFaceList[1];
      } // faceloc

      // store NewNodes in ElemList
      for (i=0; i<neQ; i++) ElemList[iElem*neQ+i] = NewNode[i];

      // loop over faces, add faces to boundary or fronts to stack
      for (faceloc=0; faceloc<nface; faceloc++){

      
	// calculate elem and face list for faceloc
	ierr = xf_Error(xf_SubList2FaceLoc(Mesh, Q, SubList, faceloc, ElemFaceList));
	if (ierr != xf_OK) return ierr;

	e0 = ElemFaceList[0];
	f0 = ElemFaceList[1];

	// if boundary, check all same and add bface (use Q1 nodes from NewNode)
	if ((ibfgrp = EG->Face[e0][f0].Group) >= 0){
	  for (i=0; i<nlayer; i++){
	    e1 = ElemFaceList[2*i+0];
	    f1 = ElemFaceList[2*i+1];
	    if (EG->Face[e1][f1].Group != ibfgrp) return xf_Error(xf_OUT_OF_BOUNDS);
	  } // i
	  ierr = xf_Error(xf_Q1NodesOnFace(EG->QBasis, Q, faceloc, &nnode, fvec));
	  if (ierr != xf_OK) return ierr;
	  if (nnode != nf1) return xf_Error(xf_OUT_OF_BOUNDS);
	  for (k=0; k<nf1; k++) 
	    BFaceList[ibfgrp][nf1*nBFace[ibfgrp]+k] = NewNode[fvec[k]];
	  nBFace[ibfgrp]++;
	}
	else{
	  // interior face
	  ierr = xf_Error(xf_NeighborAcrossFace(Mesh, 0, e0, f0, &egR, &eR, &faceR));
	  if (ierr != xf_OK) return ierr;

	  // check that all elems across faceloc are flagged the same
	  flag = ElemFlag[eR];
	  for (i=0; i<nlayer; i++){
	    e1 = ElemFaceList[2*i+0];
	    f1 = ElemFaceList[2*i+1];
	    if (EG->Face[e1][f1].Group >= 0) return xf_Error(xf_OUT_OF_BOUNDS);
	    ierr = xf_Error(xf_NeighborAcrossFace(Mesh, 0, e1, f1, &egR, &eR, &faceR));
	    if (ierr != xf_OK) return ierr;
	    if (ElemFlag[eR] != flag) return xf_Error(xf_OUT_OF_BOUNDS);
	  } // i

	  // if flagged as taken, continue
	  if (flag >= 0) continue;

	  // if not flagged, add front to stack
	  for (i=0; i<nlayer; i++){
	    ierr = xf_Error(xf_NeighborAcrossFace(Mesh, 0, ElemFaceList[2*i+0], 
						  ElemFaceList[2*i+1], &egR, &e1, &f1));
	    if (ierr != xf_OK) return ierr;

	    // get opposite face
	    ierr = xf_Error(xf_OppositeFace(EG->QBasis, f1, &f1));
	    if (ierr != xf_OK) return ierr;

	    // add elem to stack
	    ElemStack[nstack][2*i+0] = e1; 
	    ElemStack[nstack][2*i+1] = f1;
	  } // i
	  nstack++;

	  if (nstack > StackSize){
	    xf_printf("Stack overflow; increase StackSize or fix mesh.\n");
	    return xf_Error(xf_OUT_OF_BOUNDS);
	  }

	} // end else interior

      } // faceloc

      iElem++;

    } // end while istack > 0

    
    // Release intermediate memory
    xf_Release ( (void  *) ElemFlag);
    xf_Release2( (void **) ElemStack);

  }
  else{
    // no agglomeration, but still need ElemList for remap
    nElem = EG->nElem;

    // Allocate memory
    ierr = xf_Error(xf_Alloc( (void **) &ElemList, nElem*neQ, sizeof(int)));
    if (ierr != xf_OK) return ierr;

    // ElemList just stores same info as .Node structure
    for (elem=0; elem<nElem; elem++)
      for (i=0; i<neQ; i++) ElemList[elem*neQ+i] = Mesh->ElemGroup[0].Node[elem][i];
  }

  
  // if remap requested, loop over elements and move nodes
  if (RemapFlag){

    // allocate and zero out nodeflag
    ierr = xf_Error(xf_Alloc( (void **) &NodeFlag, Mesh->nNode, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    for (inode=0; inode<Mesh->nNode; inode++) NodeFlag[inode] = -1;


    // allocate and set new node coordinates
    ierr = xf_Error(xf_Alloc2( (void ***) &NewCoord, Mesh->nNode, Mesh->Dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    for (inode=0; inode<Mesh->nNode; inode++) 
      for (i=0; i<Mesh->Dim; i++) NewCoord[inode][i] = Mesh->Coord[inode][i];

    xf_printf("Remapping nodes ... ");

    // first move 3D faces using interior info, if 3D and UseInt requested
    if ((dim == 3) && (UseInt))
      for (elem=0; elem<nElem; elem++){
	ierr = xf_Error(xf_MoveFaceNodes(Mesh, Q, ElemList+neQ*elem, NodeFlag, NewCoord));
	if (ierr != xf_OK) return ierr;
        for (inode=0; inode<Mesh->nNode; inode++) 
          for (i=0; i<Mesh->Dim; i++) Mesh->Coord[inode][i] = NewCoord[inode][i];
      } // elem
    // next, move all other nodes into parametric mapping position
    for (elem=0; elem<nElem; elem++){
      ierr = xf_Error(xf_MoveNodes(Mesh, Q, ElemList+neQ*elem, NodeFlag, NewCoord, SnapFactor));
      if (ierr != xf_OK) return ierr;
      
    } // elem
    xf_printf("done.\n");

    // set new coords
    xf_Release2( (void **) Mesh->Coord);
    Mesh->Coord = NewCoord;

    xf_Release ( (void  *) NodeFlag);
  }



  /*-----------------*/
  /* Write .gri file */
  /*-----------------*/

  xf_printf("Writing .gri output file.\n");

  if (DoAgglomerate){

    // open outFile for writing
    if ((fgri =fopen(outFile ,"w"))==NULL) return xf_Error(xf_FILE_WRITE_ERROR);

    // Write header and nodes
    fprintf(fgri, "%d %d %d\n", Mesh->nNode, nElem, dim);
    for (inode=0; inode<Mesh->nNode; inode++){
      for (d=0; d<dim; d++) fprintf(fgri, "%.15E ", Mesh->Coord[inode][d]);
      fprintf(fgri, "\n");
    }

    // Write boundary faces
    fprintf(fgri, "%d\n", Mesh->nBFaceGroup);

    for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){
      fprintf(fgri, "%d %d %s\n", nBFace[ibfgrp], nf1, Mesh->BFaceGroup[ibfgrp].Title);
      for (ibface=0; ibface<nBFace[ibfgrp]; ibface++){
	for (k=0; k<nf1; k++)
	  fprintf(fgri, "%d ", BFaceList[ibfgrp][nf1*ibface + k]+1);
	fprintf(fgri, "\n");
      }
    } // ibfgrp
  

    // Write elements
    if (dim == 2)
      fprintf(fgri, "%d %d QuadLagrange\n", nElem, Q);
    else
      fprintf(fgri, "%d %d HexLagrange\n", nElem, Q);
    for (elem=0; elem<nElem; elem++){
      for (k=0; k<neQ; k++)
	fprintf(fgri, "%d ", ElemList[neQ*elem + k]+1);
      fprintf(fgri, "\n");
    } // elem

    fclose(fgri);


    xf_printf("Wrote %s\n", outFile);

    // Release memory
    xf_Release2( (void **) BFaceList);
    xf_Release ( (void  *) nBFace);

  }
  else{
    // no agglomeration took place, just write out standard .gri file
    ierr = xf_Error(xf_WriteGriFile(Mesh, outFile));
    if (ierr != xf_OK) return ierr;
  
  }

  xf_Release ( (void  *) ElemList);

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, len, Q;
  char *ArgIn[] = {"in", "NULL", "input file (with extension)",
		   "out", "NULL", "output .gri file (with extension)",
		   "Q", "2", "desired geometry order",
		   "remap", "True", "transfinite remap (use for skewed meshes)",
		   "volcheck", "False", "True to check for negative volumes",
		   "useint", "True", "True to use interior face nodes in 3d",
		   "snapfactor", "-1", ">0 to snap edges to linear (e.g. 0.1)",
		   "\0"};
  char inFile[xf_MAXSTRLEN];  
  char outFile[xf_MAXSTRLEN];
  char *inExt, *outExt;
  enum xfe_Bool RemapFlag;
  enum xfe_Bool VolCheck;
  enum xfe_Bool UseInt;
  real SnapFactor;
  xf_KeyValue KeyValue;
  xf_All *All;

  xf_printf("\n");
  xf_printf("=== xf_StructuredHO: Linear to high-order mesh conversion ===\n");
  xf_printf("\n");
    
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  xf_printf("nKey = %d\n", KeyValue.nKey);
  for (i=0; i<KeyValue.nKey; i++)
    xf_printf("%d : Key = %s, Value = %s\n", i, KeyValue.Key[i], KeyValue.Value[i]);
    
  // Get inFile
  ierr = xf_GetKeyValue(KeyValue, "in", inFile);
  if (ierr != xf_OK) return ierr;

  // Get outFile
  ierr = xf_GetKeyValue(KeyValue, "out", outFile);
  if (ierr != xf_OK) return ierr;

   // Get Q value
  ierr = xf_GetKeyValueInt(KeyValue, "Q", &Q);
  if (ierr != xf_OK) return ierr;

  // Get RemapFlag
  ierr = xf_GetKeyValueBool(KeyValue, "remap", &RemapFlag);
  if (ierr != xf_OK) return ierr;

  // Get VolCheck
  ierr = xf_GetKeyValueBool(KeyValue, "volcheck", &VolCheck);
  if (ierr != xf_OK) return ierr;

  // Get UseInt
  ierr = xf_GetKeyValueBool(KeyValue, "useint", &UseInt);
  if (ierr != xf_OK) return ierr;

  // Get SnapFactor
  ierr = xf_GetKeyValueReal(KeyValue, "snapfactor", &SnapFactor);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  
  // extensions are required
  if ((len = strlen(inFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
  inExt = inFile + len - 4; // pointer to extension

  if ((len = strlen(outFile)) < 4) return xf_Error(xf_FILE_READ_ERROR);
  outExt = outFile + len - 4; // pointer to extension
  
  // only .gri output files are supported
  if (strncmp(outExt, ".gri", 4) != 0) return xf_Error(xf_NOT_SUPPORTED);


  /*-------------*/
  /* Read inFile */
  /*-------------*/

  ierr = xf_Error(xf_ReadAllInputFile(inFile, NULL, xfe_False, &All));
  if (ierr != xf_OK) return ierr;


  /*------------------------*/
  /* Write out HO .gri file */
  /*------------------------*/
  ierr = xf_Error(xf_ConvertHO(All, Q, outFile, RemapFlag, UseInt, SnapFactor));
  if (ierr != xf_OK) return ierr;


  /*--------------------------*/
  /* Check volumes if desired */
  /*--------------------------*/

  if (VolCheck){
    xf_printf("Checking volumes.\n", Q);
    
    // do not error out on check
    ierr = xf_Error(xf_CheckVolumes(All->Mesh, NULL, xfe_False, 0, NULL, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;
  
  xf_printf("xf_StructuredHO finished.\n");

  return xf_OK;
}
