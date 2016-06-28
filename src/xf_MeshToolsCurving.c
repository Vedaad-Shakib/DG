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
 FILE:  xf_MeshToolsCurving.c
 
 This file contains functions for curving meshes to high-order
 geometry (Q > 1).

*/

#include "xf_Geom.h"
#include "xf_String.h"

/******************************************************************/
//   FUNCTION Definition: xf_SplitNewElem
static void
xf_SplitNewElem(int *elemflag, int egrp, int egrp0, int egrp1, 
		int *pegrp, int *pelem)
{
  int flag;
  if ((*pegrp) == egrp){
    flag = elemflag[(*pelem)];
    (*pegrp) = (flag < 0) ?   egrp0  : egrp1;
    (*pelem) = (flag < 0) ?  -1-flag : flag; 
  }
  
}

/******************************************************************/
//   FUNCTION Definition: xf_SplitElemGroup
static int 
xf_SplitElemGroup(xf_Mesh *Mesh, int egrp, int *elemflag, int *pegrpnew)
{
  int ierr;
  int i0, i1, inew, flag;
  int i, j, nelem, nelem0, nelem1;
  int egrp0, egrp1;
  int iiface, ibfgrp, ibface;
  xf_ElemGroup *EG, EG0, EG1, *EGnew, EGtemp;
  xf_IFace *IFace;
  xf_BFace *BFace;

  // for convenience
  EG = Mesh->ElemGroup+egrp; // pointer to element group
  nelem = EG->nElem;         // number of elements

  // count number of elements and check trivial cases
  for (i=0, nelem0=nelem1=0; i<nelem; i++)
    if (elemflag[i] < 0) nelem0++; else nelem1++;
  if (nelem0==0){
    if (pegrpnew != NULL) (*pegrpnew) = egrp;
    return xf_OK;
  }
  if (nelem1==0) return xf_Error(xf_INPUT_ERROR);

  // initialize new elem groups to orig
  EG0 = EG1 = *EG;
  
  // elem group 0
  EG0.nElem = nelem0;
  ierr = xf_Error(xf_Alloc( (void **) &EG0.nFace, nelem0, sizeof(int) ));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc2( (void ***) &EG0.Node, nelem0, EG->nNode, sizeof(int) ));
  if (ierr != xf_OK) return ierr;

  // elem group 1
  EG1.nElem = nelem1;
  ierr = xf_Error(xf_Alloc( (void **) &EG1.nFace, nelem1, sizeof(int) ));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc2( (void ***) &EG1.Node, nelem1, EG->nNode, sizeof(int) ));
  if (ierr != xf_OK) return ierr;

  // store numbering of elems in new groups in elemflag
  for (i=i0=i1=0; i<nelem; i++){
    if (elemflag[i] < 0) elemflag[i] = -(i0++)-1;
    else elemflag[i] = i1++;
  }

  // loop over orig elems: set nFace, Node
  for (i=0; i<nelem; i++){
    flag = elemflag[i];
    EGnew = (flag < 0) ?    &EG0 : &EG1;
    inew  = (flag < 0) ? -1-flag : flag;
    EGnew->nFace[inew] = EG->nFace[i];
    for (j=0; j<EG->nNode; j++) EGnew->Node[inew][j] = EG->Node[i][j];
  }

  // VAlloc Face on each allocated elem group
  ierr = xf_Error(xf_VAlloc2( (void ***) &EG0.Face, nelem0, EG0.nFace, sizeof(xf_Face)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VAlloc2( (void ***) &EG1.Face, nelem1, EG1.nFace, sizeof(xf_Face)));
  if (ierr != xf_OK) return ierr;
  
  // copy over face info
  for (i=i0=i1=0; i<nelem; i++){
    flag = elemflag[i];
    EGnew = (flag < 0) ?    &EG0 : &EG1;
    inew  = (flag < 0) ? -1-flag : flag;
    for (j=0; j<EG->nFace[i]; j++) EGnew->Face[inew][j] = EG->Face[i][j];
  }

  // reallocate Mesh->ElemGroup
  Mesh->nElemGroup++;
  ierr = xf_Error(xf_ReAlloc((void **) &Mesh->ElemGroup, Mesh->nElemGroup, 
			     sizeof(xf_ElemGroup)));
  if (ierr!=xf_OK) return ierr;

  // swap groups egrp and EG0
  swap(Mesh->ElemGroup[egrp], EG0, EGtemp);

  // store EG1 as last group in Mesh->ElemGroup
  Mesh->ElemGroup[Mesh->nElemGroup-1] = EG1;
  
  // old and new groups
  egrp0 = egrp;
  egrp1 = Mesh->nElemGroup-1;

  // set pointer to newly-created element group
  if (pegrpnew != NULL) (*pegrpnew) = egrp1;

  // point interior faces to correct groups
  for (iiface=0; iiface<Mesh->nIFace; iiface++){
    IFace = Mesh->IFace+iiface;
    xf_SplitNewElem(elemflag, egrp, egrp0, egrp1, &IFace->ElemGroupL, &IFace->ElemL);
    xf_SplitNewElem(elemflag, egrp, egrp0, egrp1, &IFace->ElemGroupR, &IFace->ElemR);
  }

  // point boundary faces to correct groups
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++)
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      BFace = Mesh->BFaceGroup[ibfgrp].BFace + ibface;
      xf_SplitNewElem(elemflag, egrp, egrp0, egrp1, &BFace->ElemGroup, &BFace->Elem);
    } // ibface

  // Destroy original group (now in EG0)
  xf_Release( (void  *) EG0.nFace);
  xf_Release2((void **) EG0.Face);
  xf_Release2((void **) EG0.Node);
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ChangeQElemGroup
static int 
xf_ChangeQElemGroup(xf_Mesh *Mesh, int egrp, int Q)
{
  int ierr, dim, i, j;
  int nelem, Q0;
  int nn, nn0;
  int nextra, ntot;
  int elem;
  int nl0, nl;
  int i0, ig;
  int nlv0[xf_MAXQ1NODE], nlv[xf_MAXQ1NODE];
  int *Node, *Node0 = NULL;
  enum xfe_Bool skip;
  real *xref = NULL, *xglob = NULL;
  xf_ElemGroup *EG = NULL;

  // for convenience
  EG = Mesh->ElemGroup+egrp; // pointer to element group
  nelem = EG->nElem;         // number of elements
  dim = Mesh->Dim;           // dim

  // store current Q
  Q0 = EG->QOrder;

  // error if Q is less than current geometry order
  if (Q  < Q0) return xf_Error(xf_INPUT_ERROR);
  if (Q == Q0) return xf_OK; // nothing to do

  // number of nodes, old and new
  ierr = xf_Error(xf_Order2nNode(EG->QBasis, Q0, &nn0));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Order2nNode(EG->QBasis, Q, &nn));
  if (ierr != xf_OK) return ierr;
  
  // determine number of new nodes, globally
  nextra = (nn-nn0)*nelem;
  if (nextra < 0) return xf_Error(xf_CODE_LOGIC_ERROR);

  // reallocate memory for global node coordinates
  ierr = xf_Error(xf_ReAllocCopy2( (void ***) &Mesh->Coord, 
				   Mesh->nNode, dim, Mesh->nNode+nextra, 
				   dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ntot = Mesh->nNode;    // original number of nodes
  Mesh->nNode += nextra; // new number of nodes
  
  // reallocate memory for Node structure
  ierr = xf_Error(xf_ReAllocCopy2((void ***) &EG->Node, nelem, 
				  nn0, nelem, nn, sizeof(int)));
  if (ierr!=xf_OK) return ierr;

  // allocate memory for ref and glob coord storage (for one elem)
  ierr = xf_Error(xf_Alloc( (void **) &xref, 2*nn*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  xglob = xref + nn*dim;

  // allocate memory for node index storage
  ierr = xf_Error(xf_Alloc( (void **) &Node0, nn0, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  
  // calculate ref coordinates for order Q geom approx
  ierr = xf_Error(xf_LagrangeNodes(EG->QBasis, Q, &i, xref, NULL));
  if (ierr != xf_OK) return ierr;
  if (i != nn) return xf_Error(xf_CODE_LOGIC_ERROR);

  // loop over elements
  for (elem=0; elem<nelem; elem++){
    
    // determine global positions of order Q geom approx nodes
    ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True, 
				    nn, xref, xglob));
    if (ierr != xf_OK) return ierr;

    // store away global indices of order Q0 (current) nodes
    for (i=0; i<nn0; i++) Node0[i] = EG->Node[elem][i];
    Node = EG->Node[elem];

    // identify local indices of Q0 (current) linear nodes
    ierr = xf_Error(xf_Q1Nodes(EG->QBasis, Q0, &nl0, nlv0));
    if (ierr != xf_OK) return ierr;
    
    // identify local indices of new linear nodes
    ierr = xf_Error(xf_Q1Nodes(EG->QBasis, Q, &nl, nlv));
    if (ierr != xf_OK) return ierr;

    // sanity check
    if (nl0 != nl) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    // set linear nodes in new Node list to same as in old list
    for (i=0; i<nl; i++) Node[nlv[i]] = Node0[nlv0[i]];
    for (i=0; i<nl; i++) Node0[nlv0[i]] = -1; // means already used

    // set all other nodes
    for (i=0, i0=0; i<nn; i++){

      // is this a Q1 node?  If so, skip.
      for (j=0, skip=xfe_False; (j<nl) && (!skip); j++) if (nlv[j] == i) skip = xfe_True;
      if (skip) continue;
      
      // identify next available global node index
      while ((i0<nn0) && (Node0[i0]<0)) i0++;
      if (i0 >= nn0) ig = ntot++;
      else ig = Node0[i0];

      // intermediate sanity check
      if (ig >= Mesh->nNode) return xf_Error(xf_CODE_LOGIC_ERROR);

      // set new node index and coordinate
      Node[i] = ig;
      for (j=0; j<dim; j++) Mesh->Coord[Node[i]][j] = xglob[i*dim+j];
      
      // increment i0 = index into available nodes from old .Node
      i0++;

    } // i

  } // elem

  // consistency check for total number of nodes
  if (ntot != Mesh->nNode) return xf_Error(xf_CODE_LOGIC_ERROR);

  // set new geom order
  EG->QOrder = Q;

  // set new nNode for elem group
  if (EG->nNode != nn0) return xf_Error(xf_CODE_LOGIC_ERROR);
  EG->nNode = nn;


  xf_Release( (void *) xref);
  xf_Release( (void *) Node0);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_EdgeNodes
static int 
xf_EdgeNodes(enum xfe_BasisType QBasis, int Q, int *pnedge, int **VEN)
{
/*

PURPOSE: 

  Determines local nodes on all edges of an element.
  
INPUTS:

  QBasis : basis describing element geometry
  Q : geometry order

OUTPUTS:  

  (*pnedge) : number of edges
  VEN: vector of high-order nodes on all edges: [nedge, Q+1]

RETURNS: Error Code

*/

  int ierr;
  int nedge, Qp2, id, iedge;
  int i0, d, i, k;
  enum xfe_ShapeType Shape;

  // determine Shape of element
  ierr = xf_Error(xf_Basis2Shape(QBasis, &Shape));
  if (ierr != xf_OK) return ierr;

  (*pnedge) = 0;

  // switch based on element shape
  switch (Shape){
  case xfe_Hexahedron:
    nedge = (*pnedge) = 12;
    Qp2 = (Q+1)*(Q+1); // vertical layer diff
    id  =  Qp2*Q;      // diff from bottom to top
    for (iedge=0; iedge<nedge; iedge++){
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
      default: return xf_Error(xf_OUT_OF_BOUNDS); break;
      }
      for (k=0, i=i0; k<Q+1; k++, i+=d) VEN[iedge][k] = i;
    } // iedge
    break;
  default: 
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemInterpolateDelta_Tri
static int 
xf_ElemInterpolateDelta_Tri(xf_Mesh *Mesh, int egrp, int elem, int *NodeFlag, 
			    real *NodeDelta)
{
  int ierr;
  int i, j, k, ip, dim, Q, gn;
  int nfnode, nf0, nf1, nQ1;
  int edge, nedge, nn, face;
  int nlv[xf_MAXQ1NODE];
  int *Node;
  int *fvec = NULL;
  enum xfe_Bool skip;
  real s, ds, tol, w, ws;
  real *xref = NULL, *X;
  real *delta = NULL;
  real phi[3], delQ1[6], del[2];
  real *dx0, *dx1;
  xf_ElemGroup *EG;

  dim   = Mesh->Dim;
  EG    = Mesh->ElemGroup+egrp;
  Q     = EG->QOrder;
  Node  = EG->Node[elem];


  /** First interpolate face interiors **/
  
  // allocate memory
  nn = (Q+1);
  ierr = xf_Error(xf_Alloc((void **) &fvec, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // loop over faces
  for (face=0; face<3; face++){
    
    // check if face interior needs to be moved based on Q1 nodes
    ierr = xf_Error(xf_Q1NodesOnFace(EG->QBasis, Q, face, &nfnode, fvec));
    if (ierr != xf_OK) return ierr;
    for (k=0, skip=xfe_True; (k<nfnode)&&(skip); k++)
      if (NodeFlag[Node[fvec[k]]] != -1) skip=xfe_False;
    if (skip) continue;

    // pull off all nodes on face
    ierr = xf_Error(xf_NodesOnFace(EG->QBasis, Q, face, &nfnode, fvec));
    if (ierr != xf_OK) return ierr;
    if (nfnode != nn) return xf_Error(xf_CODE_LOGIC_ERROR);

    // loop over interior nodes on face
    for (i=1; i<Q; i++){
      gn = Node[fvec[i]];
      if (NodeFlag[gn] != -1) continue; // already moved
      dx0 = NodeDelta + Node[fvec[0]]*dim;
      dx1 = NodeDelta + Node[fvec[Q]]*dim;
      s = ((real) i) / ((real) Q); // ref x pos on face
      for (k=0; k<dim; k++) NodeDelta[dim*gn+k] = (1.-s)*dx0[k] + s*dx1[k];
      NodeFlag[gn] = -2; // indicates moved, but not originally flagged
    } // i

  } // face


  /** Last, interpolate element interior **/

  // number of nodes in triangle
  if ((nn=EG->nNode) != (Q+1)*(Q+2)/2) return xf_Error(xf_CODE_LOGIC_ERROR);

  // copy over deltas to a local structure, delta
  ierr = xf_Error(xf_Alloc( (void **) &delta, nn*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nn; i++)
    for (k=0; k<dim; k++)
      delta[i*dim+k] = NodeDelta[Node[i]*dim+k];

  // save Q1 deltas
  ierr = xf_Error(xf_Q1Nodes(EG->QBasis, Q, &nQ1, nlv));
  if (ierr != xf_OK) return ierr;
  if (nQ1 != 3) return xf_Error(xf_INPUT_ERROR);
  for (i=0; i<3; i++)
    for (k=0; k<dim; k++)
      delQ1[i*dim+k] = delta[nlv[i]*dim+k];
  
  // allocate memory for ref coord storage
  ierr = xf_Error(xf_Alloc( (void **) &xref, dim*EG->nNode, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  
  // calculate ref coordinates for order Q geom approx
  ierr = xf_Error(xf_LagrangeNodes(EG->QBasis, Q, &i, xref, NULL));
  if (ierr != xf_OK) return ierr;
  if (i != EG->nNode) return xf_Error(xf_CODE_LOGIC_ERROR);

  // subtract Q1 deltas from all deltas
  for (i=0; i<nn; i++){
    // X = ref-space coord of node i
    X = xref + dim*i; 
    // Q1 basis functions
    phi[0] = 1-X[0]-X[1];
    phi[1] = X[0];
    phi[2] = X[1];
    // interpolate delta from Q1 nodes
    for (k=0; k<dim; k++)
      for (j=0; j<3; j++)
	delta[i*dim+k] -= phi[j]*delQ1[dim*j+k];
  } // i

  // loop over all nodes
  for (i=0; i<nn; i++){
    gn = Node[i];
    // skip if node moved
    if (NodeFlag[gn] != -1) continue;
    // X = ref-space coord of node i
    X = xref + dim*i; 
    // Q1 basis functions
    phi[0] = 1-X[0]-X[1];
    phi[1] = X[0];
    phi[2] = X[1];
    // skip if not interior node
    tol = 0.25/((real) Q); // floating point check is easy and accurate here
    for (k=0, j=0; k<nQ1; k++) j += (phi[k] < tol);
    if (j > 0) continue;
    // set node delta to interpolated Q1 delta
    for (k=0; k<dim; k++)
      for (j=0; j<3; j++)
	NodeDelta[gn*dim+k] += phi[j]*delQ1[dim*j+k];

    // add effect of 3 faces, in turn
    for (face=0; face<3; face++){

      // get indices of nodes on face
      ierr = xf_Error(xf_NodesOnFace(EG->QBasis, Q, face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;

      switch (face){
      case 0: w = X[0]+X[1]; s = X[1]/w; break;
      case 1: w = 1.0 -X[0]; s = X[1]/w; break;
      case 2: w = 1.0 -X[1]; s = X[0]/w; break;
      default: return xf_Error(xf_OUT_OF_BOUNDS); break;
      }
      
      // using s and deltas on face, set local delta, del (linear interp)
      ds = 1.0/((real) Q);
      ip = s/ds;
      ws = 1. - (s-ip*ds);
      for (k=0; k<dim; k++) 
	del[k] = ws*delta[dim*fvec[ip]+k] + (1.0-ws)*delta[dim*fvec[ip+1]+k];

      // multiply del by weight, add to NodeDelta
      for (k=0; k<dim; k++) NodeDelta[gn*dim+k] += del[k]*w;

    } // face
    
    // Indicates node has moved, but not originally flagged
    NodeFlag[gn] = -2; 
  }

  xf_Release( (void *) delta);
  xf_Release( (void *) xref);
 
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemInterpolateDelta_Quad
static int 
xf_ElemInterpolateDelta_Quad(xf_Mesh *Mesh, int egrp, int elem, int *NodeFlag, 
			     real *NodeDelta)
{
  int ierr;
  int i, j, k, l, dim, Q, gn;
  int ii, jj, ij, iijj;
  int nfnode, nf0, nf1;
  int edge, nedge, nn, face;
  int ivec[3], jvec[3];
  int *Node;
  int *fvec = NULL;
  enum xfe_Bool skip;
  real s;
  real xref[2];
  real *delta[9];
  real *dx0, *dx1;
  xf_ElemGroup *EG;

  dim   = Mesh->Dim;
  EG    = Mesh->ElemGroup+egrp;
  Q     = EG->QOrder;
  Node  = EG->Node[elem];


  /** First interpolate face interiors **/
  
  // allocate memory
  nn = (Q+1);
  ierr = xf_Error(xf_Alloc((void **) &fvec, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // loop over faces
  for (face=0; face<4; face++){
    
    // check if face interior needs to be moved based on Q1 nodes
    ierr = xf_Error(xf_Q1NodesOnFace(EG->QBasis, Q, face, &nfnode, fvec));
    if (ierr != xf_OK) return ierr;
    for (k=0, skip=xfe_True; (k<nfnode)&&(skip); k++)
      if (NodeFlag[Node[fvec[k]]] != -1) skip=xfe_False;
    if (skip) continue;

    // pull off all nodes on face
    ierr = xf_Error(xf_NodesOnFace(EG->QBasis, Q, face, &nfnode, fvec));
    if (ierr != xf_OK) return ierr;
    if (nfnode != nn) return xf_Error(xf_CODE_LOGIC_ERROR);

    // loop over interior nodes on face
    for (i=1; i<Q; i++){
      gn = Node[fvec[i]];
      if (NodeFlag[gn] != -1) continue; // already moved
      dx0 = NodeDelta + Node[fvec[0]]*dim;
      dx1 = NodeDelta + Node[fvec[Q]]*dim;
      s = ((real) i) / ((real) Q); // ref x pos on face
      for (k=0; k<dim; k++) NodeDelta[dim*gn+k] = (1.-s)*dx0[k] + s*dx1[k];
      NodeFlag[gn] = -2; // indicates moved, but not originally flagged
    } // i

  } // face

  // release memory
  xf_Release( (void *) fvec);

  /** Last, interpolate element interior **/
  for (j=1; j<Q; j++)    
    for (i=1; i<Q; i++){
	
      // ref coords (equispaced)
      xref[0] = ((real) i) / ((real) Q);
      xref[1] = ((real) j) / ((real) Q);

      // build delta for transfinite interpolation
      jvec[0]=0; jvec[1]=j; jvec[2]=Q;
      ivec[0]=0; ivec[1]=i; ivec[2]=Q;
      for (jj=0, k=0; jj<3; jj++)
	for (ii=0; ii<3; ii++){
	  iijj = jvec[jj]*(Q+1) + ivec[ii];
	  delta[k++] = NodeDelta + Node[iijj]*dim;
	}
	
      // call transfinite interpolation to obtain delta for face-interior node
      xf_Transfinite2D(dim, xref, delta);

      // flag node as moved, but not originally flagged
      ij = j*(Q+1) + i;
      NodeFlag[Node[ij]] = -2;
	
    } // i
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_ElemInterpolateDelta_Hex
static int 
xf_ElemInterpolateDelta_Hex(xf_Mesh *Mesh, int egrp, int elem, int *NodeFlag, 
			    real *NodeDelta)
{
  int ierr;
  int i, j, k, l, dim, Q, gn;
  int ii, jj, kk, ij, ijk, iijj;
  int nfnode, nf0, nf1;
  int edge, nedge, nn, face;
  int ivec[3], jvec[3], kvec[3];
  int *Node;
  int *fvec = NULL;
  int **VEN = NULL, *EN=NULL;
  enum xfe_Bool skip;
  real s;
  real xref[3];
  real *delta[27];
  real *dx0, *dx1;
  xf_ElemGroup *EG;

  dim   = Mesh->Dim;
  EG    = Mesh->ElemGroup+egrp;
  Q     = EG->QOrder;
  Node  = EG->Node[elem];

  /** First interpolate edges **/

  // number of edges (nedge) and nodes on each edge (EN[*,*])
  ierr = xf_Error(xf_Alloc2((void ***) &VEN, xf_MAXEDGES, Q+1, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_EdgeNodes(EG->QBasis, Q, &nedge, VEN));
  if (ierr != xf_OK) return ierr;

  for (edge=0; edge<nedge; edge++){
    EN = VEN[edge];
    nf0 = NodeFlag[Node[EN[0]]]; // flag on first node of edge
    nf1 = NodeFlag[Node[EN[Q]]]; // flag on last node of edge
    dx0 = NodeDelta + dim*Node[EN[0]]; // delta on first node
    dx1 = NodeDelta + dim*Node[EN[Q]]; // delta on last node
    if ((nf0>-1) && (nf1>-1)){ // both ends flagged
      for (j=1,ierr=0; j<Q; j++) ierr += (NodeFlag[Node[EN[j]]]==-1);
      if (ierr != 0) return xf_Error(xf_CODE_LOGIC_ERROR); // in-between should be flagged
    }
    else if (nf0*nf1 <= 0){ // one flagged, other not
      for (j=1; j<Q; j++){
	gn = Node[EN[j]]; // global node number
	if (NodeFlag[gn] != -1) continue; // already moved
	s = ((real) j) / ((real) Q); // dist along edge, for linear interp of delta
	for (k=0; k<dim; k++) NodeDelta[dim*gn+k] = (1.-s)*dx0[k] + s*dx1[k];
	NodeFlag[gn] = -2; // indicates moved, but not originally flagged
      }
    }
  } // edge

  xf_Release2( (void **) VEN);

  
  /** Next interpolate face interiors **/
  
  // allocate memory
  nn = (Q+1)*(Q+1);
  ierr = xf_Error(xf_Alloc((void **) &fvec, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // loop over faces
  for (face=0; face<6; face++){
    
    // check if face interior needs to be moved based on Q1 nodes
    ierr = xf_Error(xf_Q1NodesOnFace(EG->QBasis, Q, face, &nfnode, fvec));
    if (ierr != xf_OK) return ierr;
    for (k=0, skip=xfe_True; (k<nfnode)&&(skip); k++)
      if (NodeFlag[Node[fvec[k]]] != -1) skip=xfe_False;
    if (skip) continue;

    // pull off all nodes on face
    ierr = xf_Error(xf_NodesOnFace(EG->QBasis, Q, face, &nfnode, fvec));
    if (ierr != xf_OK) return ierr;
    if (nfnode != nn) return xf_Error(xf_CODE_LOGIC_ERROR);

    // loop over interior nodes on face
    for (j=1; j<Q; j++){
      for (i=1; i<Q; i++){
	ij = j*(Q+1) + i;
	if (NodeFlag[Node[fvec[ij]]] != -1) continue; // already moved
	xref[0] = ((real) i) / ((real) Q); // ref x pos on face
	xref[1] = ((real) i) / ((real) Q); // ref y pos on edge

	// build delta for transfinite interpolation
	jvec[0]=0; jvec[1]=j; jvec[2]=Q;
	ivec[0]=0; ivec[1]=i; ivec[2]=Q;
	for (jj=0, k=0; jj<3; jj++)
	  for (ii=0; ii<3; ii++){
	    iijj = jvec[jj]*(Q+1) + ivec[ii];
	    delta[k++] = NodeDelta + Node[fvec[iijj]]*dim;
	  }
	
	// call transfinite interpolation to obtain delta for face-interior node
	xf_Transfinite2D(dim, xref, delta);

	// flag node as moved, but not originally flagged
	NodeFlag[Node[fvec[ij]]] = -2;

      } // j
    } // i

  } // face

  // release memory
  xf_Release( (void *) fvec);

  /** Last, interpolate element interior **/
  for (k=1; k<Q; k++)
    for (j=1; j<Q; j++)    
      for (i=1; i<Q; i++){
	
	// ref coords (equispaced)
	xref[0] = ((real) i) / ((real) Q);
	xref[1] = ((real) j) / ((real) Q);
	xref[2] = ((real) k) / ((real) Q);

	// get coord list for transfinite interpolation	
	kvec[0]=0; kvec[1]=k; kvec[2]=Q;
	jvec[0]=0; jvec[1]=j; jvec[2]=Q;
	ivec[0]=0; ivec[1]=i; ivec[2]=Q;
	for (kk=0, l=0; kk<3; kk++)
	  for (jj=0; jj<3; jj++)
	    for (ii=0; ii<3; ii++){
	      ijk = kvec[kk]*(Q+1)*(Q+1) + jvec[jj]*(Q+1) + ivec[ii];
	      delta[l++] = NodeDelta + Node[ijk]*dim;
	    }
	
	// call transfinite interpolation to obtain delta for elem-interior node
	xf_Transfinite3D(dim, xref, delta);
	
	// flag interior node as moved, but not originally flagged
	ijk = k*(Q+1)*(Q+1) + j*(Q+1) + i;
	if (NodeFlag[Node[ijk]] != -1) return xf_Error(xf_CODE_LOGIC_ERROR);
	NodeFlag[Node[ijk]] = -2;
	
      } // i
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ElemInterpolateDelta
static int 
xf_ElemInterpolateDelta(xf_Mesh *Mesh, int egrp, int elem, int *NodeFlag, 
			real *NodeDelta)
{
  int ierr;
  enum xfe_ShapeType Shape;

  // determine Shape of element
  ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape));
  if (ierr != xf_OK) return ierr;

  // switch based on element shape
  switch (Shape){
  case xfe_Triangle:
    ierr = xf_Error(xf_ElemInterpolateDelta_Tri(Mesh, egrp, elem, NodeFlag, NodeDelta));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_Quadrilateral:
    ierr = xf_Error(xf_ElemInterpolateDelta_Quad(Mesh, egrp, elem, NodeFlag, NodeDelta));
    if (ierr != xf_OK) return ierr;
    break;
  case xfe_Hexahedron:
    ierr = xf_Error(xf_ElemInterpolateDelta_Hex(Mesh, egrp, elem, NodeFlag, NodeDelta));
    if (ierr != xf_OK) return ierr;
    break;
  default: 
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_NodesElemsOnBoundaries
static int 
xf_NodesElemsOnBoundaries(xf_Mesh *Mesh, const char *BFGTitles, 
			  int **pNodeFlag, int ***pElemFlag)
{
  int ierr, dim, i, j, nn;
  int negrp0, nelemtot;
  int egrp, egrp0, elem, face;
  int ibfgrp, ibface;
  int count, nBFG, nfnode, node;
  int nedge, edge, nf0, nf1;
  int *nElem = NULL;
  int *fvec = NULL;
  int *NodeFlag = NULL;
  int **ElemFlag = NULL;
  int **VEN = NULL, *EN;
  char **BFGTitleArray = NULL;
  enum xfe_Bool CurveThis = xfe_False;
  xf_BFace *BFace;
  xf_ElemGroup *EG;

  dim = Mesh->Dim;
  
  // number of original element groups
  negrp0 = Mesh->nElemGroup;

  // determine bface groups with which we are working
  ierr = xf_Error(xf_ScanXStringAlloc(BFGTitles, xf_MAXSTRLEN, &nBFG, &BFGTitleArray));
  if (ierr != xf_OK) return ierr;

  // create and initialize NodeFlag
  ierr = xf_Error(xf_Alloc( (void **) &NodeFlag, Mesh->nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<Mesh->nNode; i++) NodeFlag[i] = -1;

  // over-allocate fvec using max # nodes per element
  for (egrp=0, nn=-1; egrp<Mesh->nElemGroup; egrp++)
    nn = max(nn, Mesh->ElemGroup[egrp].nNode);
  ierr = xf_Error(xf_Alloc((void **) &fvec, nn, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // allocate ElemFlag, initialize to -1
  ierr = xf_Error(xf_GetnElem(Mesh, &nElem, &nelemtot));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_VAlloc2((void ***) &ElemFlag, negrp0, nElem, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nelemtot; i++) ElemFlag[0][i] = -1; // initialize to -1
  xf_Release( (void *) nElem);

  // mark nodes by looping over boundary groups and looking at adjacent elements
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){
    // are we curving this boundary?  Check input
    for (i=0, CurveThis=xfe_False; i<nBFG; i++)
      if (strcmp(Mesh->BFaceGroup[ibfgrp].Title, BFGTitleArray[i]) == 0) CurveThis = xfe_True;
    if (!CurveThis) continue;
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      BFace = Mesh->BFaceGroup[ibfgrp].BFace+ibface;
      EG = Mesh->ElemGroup + BFace->ElemGroup;
      ierr = xf_Error(xf_NodesOnFace(EG->QBasis, EG->QOrder, BFace->Face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;
      // mark nodes as on boundary using boundary group number, ibfgrp
      for (j=0; j<nfnode; j++) NodeFlag[EG->Node[BFace->Elem][fvec[j]]] = ibfgrp; 
      // mark element by setting ElemFlag
      ElemFlag[BFace->ElemGroup][BFace->Elem] = 1;
    }
  } // ibfgrp

  // release some memory
  xf_Release2( (void **) BFGTitleArray);
  xf_Release ( (void  *) fvec);

  // in 3D, there may be elements with only edges on boundaries -- find these
  if (dim == 3){
    VEN = NULL;
    for (egrp=0; egrp<negrp0; egrp++){
      EG = Mesh->ElemGroup + egrp;
      // Determine number of edges (nedge) and nodes on each edge (VEN[*,*])
      ierr = xf_Error(xf_ReAlloc2((void ***) &VEN, xf_MAXEDGES, EG->QOrder+1, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_EdgeNodes(EG->QBasis, EG->QOrder, &nedge, VEN));
      if (ierr != xf_OK) return ierr;

      for (elem=0; elem<EG->nElem; elem++){  // loop over elements
	if (ElemFlag[egrp][elem] != -1) continue; // elem already flagged
	for (edge=0; edge<nedge; edge++){ // loop over edges
	  EN = VEN[edge];
	  nf0 = NodeFlag[EG->Node[elem][EN[0]]]; // flag on first node of edge
	  nf1 = NodeFlag[EG->Node[elem][EN[EG->QOrder]]]; // flag on last node of edge
	  if ((nf0 >= 0) && (nf0==nf1)){  // flags are nonzero and agree
	    for (j=1; j<EG->QOrder; j++)
	      if (NodeFlag[EG->Node[elem][EN[j]]] == -1) 
		NodeFlag[EG->Node[elem][EN[j]]] = nf0; // mark high-order edge nodes
	    ElemFlag[egrp][elem] = 1; // mark element
	  }
	  // error if edge straddles two boundaries
	  else if ((nf0 >= 0) && (nf1 >= 0)) return xf_OUT_OF_BOUNDS;
	}
      } // elem
    } // egrp

    // release some memory
    xf_Release2( (void **) VEN);
  }

  if (pNodeFlag != NULL) (*pNodeFlag) = NodeFlag;
  else xf_Release( (void  *) NodeFlag);

  if (pElemFlag != NULL) (*pElemFlag) = ElemFlag;
  else xf_Release2((void **) ElemFlag);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CurveMeshBoundary
int 
xf_CurveMeshBoundary(xf_Mesh *Mesh, xf_Geom *Geom, 
		     const char *BFGTitles, int Q)
{
  int ierr, dim, i, j, nn;
  int negrp0, nelemtot;
  int egrp, egrp0, egrp1, elem, face;
  int ibfgrp, ibface;
  int count, nBFG, nfnode, node;
  int nedge, edge, nf0, nf1;
  int *nElem = NULL;
  int *fvec = NULL;
  int *NodeFlag = NULL;
  int **ElemFlag = NULL;
  int **VEN = NULL, *EN;
  char **BFGTitleArray = NULL;
  enum xfe_Bool split = xfe_False;
  enum xfe_Bool CurveThis = xfe_False;
  enum xfe_Bool CurvingFlag = xfe_False;
  real *NodeDelta = NULL;
  real *delta = NULL;
  real x0[3];
  xf_BFace *BFace;
  xf_ElemGroup *EG;

  dim = Mesh->Dim;
  
  // number of original element groups
  negrp0 = Mesh->nElemGroup;

  // determine bface groups with which we are working
  ierr = xf_Error(xf_ScanXStringAlloc(BFGTitles, xf_MAXSTRLEN, &nBFG, &BFGTitleArray));
  if (ierr != xf_OK) return ierr;


  /*--------------------------------*/
  /*  Identify nodes on boundaries  */
  /*--------------------------------*/

  ierr = xf_Error(xf_NodesElemsOnBoundaries(Mesh, BFGTitles, &NodeFlag, &ElemFlag));
  if (ierr != xf_OK) return ierr;


  /*-------------------------------------------*/
  /*  Separate off high Q groups if requested  */
  /*-------------------------------------------*/

  if (Q > 0){
    
    /** split off flagged elements from groups if necessary **/
    for (egrp0=0; egrp0<negrp0; egrp0++){
      // should we split this element group?
      for (elem=0, split=xfe_False, count=0; elem<Mesh->ElemGroup[egrp0].nElem; elem++)
	if (ElemFlag[egrp0][elem] > -1){
	  split = xfe_True;
	  count++;
	}
      // do not split if all elements flagged
      if (count == Mesh->ElemGroup[egrp0].nElem) split == xfe_False;
	
      if (split){
	// Split elements off to a new group
	ierr = xf_Error(xf_SplitElemGroup(Mesh, egrp0, ElemFlag[egrp0], &egrp1));
	if (ierr != xf_OK) return ierr;
	xf_printf("Splitting egrp0=%d; created egrp1=%d\n", egrp0, egrp1);
      }

      // make new group HO if Q is higher than current one
      if (Q > Mesh->ElemGroup[egrp0].QOrder){
	// increase Q of the newly-created group, egrp1
	ierr = xf_Error(xf_ChangeQElemGroup(Mesh, egrp1, Q));
	if (ierr != xf_OK) return ierr;
      }
      
    } // egrp0
    
    CurvingFlag = xfe_True;

    /*  Identify nodes on boundaries again, since element groups/nodes changed */
    xf_Release2( (void **) ElemFlag);
    xf_Release ( (void  *) NodeFlag);
    ierr = xf_Error(xf_NodesElemsOnBoundaries(Mesh, BFGTitles, &NodeFlag, &ElemFlag));
    if (ierr != xf_OK) return ierr;

  } // end if asking for high Q


  // If geometry not available, cannot do any more
  if (Geom == NULL){
    xf_Release2( (void **) ElemFlag);
    xf_Release ( (void  *) NodeFlag);
    return xf_OK;
  }
  
  /*-------------------------------------------------------------*/
  /*  Compute a delta vector for each node requiring projection  */
  /*-------------------------------------------------------------*/

  xf_printf("Projecting nodes to boundaries.\n");

  // create a delta array over all nodes
  ierr = xf_Error(xf_Alloc( (void **) &NodeDelta, Mesh->nNode*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<Mesh->nNode*dim; i++) NodeDelta[i] = 0.;

  // project boundary nodes to geometry, get deltas
  for (i=0; i<Mesh->nNode; i++){
    if ((ibfgrp=NodeFlag[i]) > -1){
      for (j=0; j<dim; j++) x0[j] = Mesh->Coord[i][j];
      ierr = xf_Error(xf_ProjectToGeom(Geom, -1, Mesh->BFaceGroup[ibfgrp].Title,
				       dim, 1, x0));
      if (ierr != xf_OK) return ierr;
      for (j=0; j<dim; j++) NodeDelta[i*dim+j] = x0[j] - Mesh->Coord[i][j];
    }
  } // i


  /*------------------------------------------*/
  /* Extend delta vector to affected elements */
  /*------------------------------------------*/

  xf_printf("Extending delta vector to element interiors.\n");

  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      if (ElemFlag[egrp][elem] > -1){ // element flagged for curving

	// Determine deltas for nodes inside (and on faces of) this element
	ierr = xf_Error(xf_ElemInterpolateDelta(Mesh, egrp, elem, NodeFlag, NodeDelta));
	if (ierr != xf_OK) return ierr;

      }
    } // elem
  } // egrp


  /*----------------------*/
  /*  Apply delta vector  */
  /*----------------------*/
  
  for (i=0; i<Mesh->nNode; i++)
    for (j=0; j<dim; j++) Mesh->Coord[i][j] += NodeDelta[i*dim+j];
  

  // Release memory
  xf_Release2( (void **) ElemFlag);
  xf_Release ( (void  *) NodeDelta);
  xf_Release ( (void  *) NodeFlag);
 

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_SnapToGeom
int 
xf_SnapToGeom(xf_All *All)
{
  int ierr, iC;
  xf_Mesh *Mesh;
  xf_Geom *Geom;
  char BFGTitles[xf_MAXSTRLEN];
  char instr[xf_MAXSTRLEN];
  char *s;
  
  // nothing to do if Mesh or Geom do not exist
  if ((Mesh=All->Mesh) == NULL) return xf_OK;
  if ((Geom=All->Geom) == NULL) return xf_OK;

  // loop over components, pull off all BFGTitles and store into one string
  sprintf(BFGTitles, "\0");
  for (iC=0; iC<Geom->nComp; iC++){
    if ((s = Geom->Comp[iC].BFGTitle) != NULL){
      sprintf(instr, "%s %s", BFGTitles, s);
      strcpy(BFGTitles, instr);
    }
  } // iC
  
  // nothing to do if no BFGs selected
  if (strlen(BFGTitles) == 0) return xf_OK;

  // curve mesh, do not change Q (pass in -1)
  xf_printf("Snapping the following boundaries: %s\n", BFGTitles);
  ierr = xf_Error(xf_CurveMeshBoundary(Mesh, Geom, BFGTitles, -1));
  if (ierr != xf_OK) return ierr;

  return xf_OK;

}




/*-------------*/
/*   LEGACY    */
/*-------------*/


/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_CurveElementInterior_TriQuad */
/* static int  */
/* xf_CurveElementInterior_TriQuad(xf_Mesh *Mesh, int egrp, int elem, int face, real *delta) */
/* { */
/*   int ierr, i, j, k, dim; */
/*   int ip, Q, nQ1, nNode, nfnode; */
/*   int *Node; */
/*   int *fvec = NULL; */
/*   int nlv[xf_MAXQ1NODE]; */
/*   enum xfe_Bool skip; */
/*   xf_ElemGroup *EG = NULL; */
/*   real s, ds, tol, w, ws; */
/*   real *xref = NULL, *X, *x; */
/*   real phi[4], del[2], deltaQ1[8]; */

/*   // for convenience */
/*   EG    = Mesh->ElemGroup+egrp; */
/*   Q     = EG->QOrder; */
/*   nNode = EG->nNode;               // number of nodes */
/*   Node  = EG->Node[elem];          // pointer to global node list */

/*   // spatial dimension should be 2 */
/*   if ((dim=Mesh->Dim) != 2) return xf_Error(xf_INPUT_ERROR); */

/*   // get indices of nodes on face */
/*   ierr = xf_Error(xf_Alloc((void **) &fvec, nNode, sizeof(int))); // over-allocate */
/*   if (ierr != xf_OK) return ierr; */
/*   ierr = xf_Error(xf_NodesOnFace(EG->QBasis, Q, face, &nfnode, fvec)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // get Q1 nodes */
/*   ierr = xf_Error(xf_Q1Nodes(EG->QBasis, Q, &nQ1, nlv)); */
/*   if (ierr != xf_OK) return ierr; */
/*   if ((nQ1 != 3) && (nQ1 != 4)) return xf_Error(xf_INPUT_ERROR); */

/*   // allocate memory for ref coord storage */
/*   ierr = xf_Error(xf_Alloc( (void **) &xref, dim*nNode, sizeof(real))); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   // calculate ref coordinates for order Q geom approx */
/*   ierr = xf_Error(xf_LagrangeNodes(EG->QBasis, Q, &i, xref, NULL)); */
/*   if (ierr != xf_OK) return ierr; */
/*   if (i != nNode) return xf_Error(xf_CODE_LOGIC_ERROR); */

/*   // subtract edge endpoint change from edge interior changes (will be superimposing effects) */
/*   for (i=1; i<(nfnode-1); i++){ */
/*     // apply perturbation now, as edge nodes will be skipped in general loop over all nodes */
/*     x = Mesh->Coord[Node[fvec[i]]]; */
/*     for (k=0; k<dim; k++) x[k] += delta[fvec[i]*dim+k]; */

/*     s = ((real) i)/ ((real) (nfnode-1.0)); */
/*     for (k=0; k<dim; k++) */
/*       delta[fvec[i]*dim+k] -= (delta[fvec[0]*dim+k]*(1.-s) + delta[fvec[nfnode-1]*dim+k]*s); */
/*   }     */

/*   // store Q1 deltas and set values to zero in the delta vector */
/*   for (k=0; k<dim; k++) */
/*     for (j=0; j<nQ1; j++){ */
/*       deltaQ1[j*dim+k] = delta[dim*nlv[j]+k]; */
/*       delta[dim*nlv[j]+k] = 0.; */
/*     } */

/*   // loop over all nodes */
/*   for (i=0; i<nNode; i++){ */
    
/*     // skip Q1 nodes and edge-interior nodes from face in question */
/*     skip = xfe_False; */
/*     for (j=0; (j<nQ1) && (!skip); j++) if (nlv[j]==i) skip = xfe_True; */
/*     for (j=1; (j<(nfnode-1)) && (!skip); j++) if (fvec[j]==i) skip = xfe_True; */
/*     if (skip) continue; */
    
/*     // X = ref-space coord of node i */
/*     X = xref + dim*i; */

/*     // Q1 basis functions */
/*     if (nQ1 == 3){ // Triangle */
/*       phi[0] = 1-X[0]-X[1]; */
/*       phi[1] = X[0]; */
/*       phi[2] = X[1]; */
/*     } */
/*     else{ // Quadrilateral */
/*       phi[0] = (1-X[0])*(1-X[1]); */
/*       phi[1] = X[0]*(1-X[1]); */
/*       phi[2] = (1-X[0])*X[1]; */
/*       phi[3] = X[0]*X[1]; */
/*     } */

/*     // interpolate delta from Q1 nodes to get local delta */
/*     for (k=0; k<dim; k++) */
/*       for (j=0, del[k]=0.; j<nQ1; j++) */
/* 	del[k] += phi[j]*deltaQ1[dim*j+k]; */

/*     // add to position of global node */
/*     x = Mesh->Coord[Node[i]]; */
/*     for (k=0; k<dim; k++) x[k] += del[k]; */

/*     // At this point, we are done if not an interior node */
/*     tol = 0.25/((real) Q); // floating point check is easy and accurate here */
/*     for (k=0, j=0; k<nQ1; k++) j += (phi[k] < tol); */
/*     if (j > 0) continue; */

/*     // on an interior node, we also need to account for effect of edge motion */

/*     if (nQ1 == 3){ // Triangle */
/*       // w = weight = one minus distance from face (X+Y, 1-X, 1-Y) */
/*       // s = position proj onto face ( Y/(X+Y), Y/(1-X), X/(1-Y),  */
/*       switch (face){ */
/*       case 0: w = X[0]+X[1]; s = X[1]/w; break; */
/*       case 1: w = 1.0 -X[0]; s = X[1]/w; break; */
/*       case 2: w = 1.0 -X[1]; s = X[0]/w; break; */
/*       default: return xf_Error(xf_OUT_OF_BOUNDS); break; */
/*       } */
      
/*       // using s and deltas on face, set local delta, del (linear interp) */
/*       ds = 1.0/((real) Q); */
/*       ip = s/ds; */
/*       ws = 1. - (s-ip*ds); */
/*       for (k=0; k<dim; k++)  */
/* 	del[k] = ws*delta[dim*fvec[ip]+k] + (1.0-ws)*delta[dim*fvec[ip+1]+k]; */
      
/*     } */
/*     else{ // Quadrilateral */
/*       switch (face){ */
/*       case 0: w = 1-X[1]; s =   X[0]; break; */
/*       case 1: w =   X[0]; s =   X[1]; break; */
/*       case 2: w =   X[1]; s = 1-X[0]; break; */
/*       case 3: w = 1-X[0]; s = 1-X[1]; break; */
/*       default: return xf_Error(xf_OUT_OF_BOUNDS); break; */
/*       } */
      
/*       // using s and deltas on face, set local delta, del */
/*       ds = 1.0/((real) Q); */
/*       ip = (s+0.1*ds)/ds; */
/*       for (k=0; k<dim; k++) del[k] = delta[dim*fvec[ip]+k];       */

/*     } */

/*     // multiply del by weight */
/*     for (k=0; k<dim; k++) del[k] *= w; */

/*     // change global coordinate using del */
/*     x = Mesh->Coord[Node[i]]; */
/*     for (k=0; k<dim; k++) x[k] += del[k]; */

/*   } // i */

/*   xf_Release( (void *) fvec); */
/*   xf_Release( (void *) xref); */

/*   return xf_OK; */
/* } */



/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_CurveElementInterior */
/* static int  */
/* xf_CurveElementInterior(xf_Mesh *Mesh, int egrp, int elem, int face,  real *delta) */
/* { */
/*   int ierr; */
/*   enum xfe_ShapeType Shape; */

/*   // determine Shape of element */
/*   ierr = xf_Error(xf_Basis2Shape(Mesh->ElemGroup[egrp].QBasis, &Shape)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // switch based on element shape */
/*   switch (Shape){ */
/*   case xfe_Triangle: */
/*   case xfe_Quadrilateral: */
/*     ierr = xf_Error(xf_CurveElementInterior_TriQuad(Mesh, egrp, elem, face, delta)); */
/*     if (ierr != xf_OK) return ierr; */
/*     break; */
/*   default:  */
/*     return xf_Error(xf_NOT_SUPPORTED); */
/*     break; */
/*   } */

/*   return xf_OK; */
/* } */


/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_CurveMeshBoundaryOld */
/* int  */
/* xf_CurveMeshBoundaryOld(xf_Mesh *Mesh, xf_Geom *Geom,  */
/* 		     const char *BFGTitles, int Q) */
/* { */
/*   int ierr, dim, i, j, nn; */
/*   int negrp0, nelemtot; */
/*   int egrp, egrp0, egrp1, elem, face; */
/*   int ibfgrp, ibface; */
/*   int count, nBFG, nfnode, node; */
/*   int *nElem = NULL; */
/*   int *fvec = NULL; */
/*   int **ElemInd = NULL; */
/*   char **BFGTitleArray = NULL; */
/*   enum xfe_Bool split = xfe_False; */
/*   enum xfe_Bool CurveThis = xfe_False; */
/*   enum xfe_Bool CurvingFlag = xfe_False; */
/*   real *NodeDelta = NULL; */
/*   real *delta = NULL; */
/*   real *xvec = NULL; */
/*   xf_BFace *BFace; */
/*   xf_ElemGroup *EG; */

/*   dim = Mesh->Dim; */
  
/*   // number of original element groups */
/*   negrp0 = Mesh->nElemGroup; */

/*   // determine bface groups with which we are working */
/*   ierr = xf_Error(xf_ScanXStringAlloc(BFGTitles, xf_MAXSTRLEN, &nBFG, &BFGTitleArray)); */
/*   if (ierr != xf_OK) return ierr; */

  
/*   /\** First, ensure that Q1 nodes are on boundary **\/ */

/*   /\** identify elements for curving, place into separate element groups **\/ */

/*   CurvingFlag = xfe_False; */

/*   if (Q > 0){ */

/*     // first create an indicator over elements */
/*     ierr = xf_Error(xf_GetnElem(Mesh, &nElem, &nelemtot)); */
/*     if (ierr != xf_OK) return ierr; */
/*     ierr = xf_Error(xf_VAlloc2((void ***) &ElemInd, negrp0, nElem, sizeof(int))); */
/*     if (ierr != xf_OK) return ierr; */
/*     for (i=0; i<nelemtot; i++) ElemInd[0][i] = -1; // initialize to -1 */
    
/*     // loop over boundary groups, fill in ElemInd */
/*     for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){ */
/*       // are we curving this boundary?  Check input */
/*       for (i=0, CurveThis=xfe_False; i<nBFG; i++) */
/* 	if (strcmp(Mesh->BFaceGroup[ibfgrp].Title, BFGTitleArray[i]) == 0) CurveThis = xfe_True; */
/*       if (!CurveThis) continue; */
/*       for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){ */
/* 	BFace = Mesh->BFaceGroup[ibfgrp].BFace+ibface; */
/* 	ElemInd[BFace->ElemGroup][BFace->Elem] = ibfgrp; */
/*       } */
/*     } */
    
/*     /\** split off flagged elements from groups if necessary **\/ */
/*     for (egrp0=0; egrp0<negrp0; egrp0++){ */
/*       // should we split this element group? */
/*       for (elem=0, split=xfe_False, count=0; elem<Mesh->ElemGroup[egrp0].nElem; elem++) */
/* 	if (ElemInd[egrp0][elem] > -1){ */
/* 	  split = xfe_True; */
/* 	  count++; */
/* 	} */
/*       // do not split if all elements flagged */
/*       if (count == Mesh->ElemGroup[egrp0].nElem) split == xfe_False; */
	
/*       if (split){ */
/* 	// Split elements off to a new group */
/* 	ierr = xf_Error(xf_SplitElemGroup(Mesh, egrp0, ElemInd[egrp0], &egrp1)); */
/* 	if (ierr != xf_OK) return ierr; */
/* 	xf_printf("Splitting egrp0=%d; created egrp1=%d\n", egrp0, egrp1); */
/*       } */

/*       // make new group HO if Q is higher than current one */
/*       if (Q > Mesh->ElemGroup[egrp0].QOrder){ */
/* 	// increase Q of the newly-created group, egrp1 */
/* 	ierr = xf_Error(xf_ChangeQElemGroup(Mesh, egrp1, Q)); */
/* 	if (ierr != xf_OK) return ierr; */
/*       } */
      
/*     } // egrp0 */
    
/*     // release memory */
/*     xf_Release( (void *) nElem); */
/*     xf_Release2( (void **) ElemInd); */
    
/*     CurvingFlag = xfe_True; */

/*   } // end if asking for high Q */

/*   // project elements adjacent to boundaries to geometry if available */
/*   if (Geom != NULL){ */

/*     // create a delta array over all nodes */
/*     ierr = xf_Error(xf_Alloc( (void **) &NodeDelta, Mesh->nNode*dim, sizeof(real))); */
/*     if (ierr != xf_OK) return ierr; */
/*     for (i=0; i<Mesh->nNode*dim; i++) NodeDelta[i] = 0.; */

/*     // get max # nodes per element */
/*     for (egrp=0, nn=-1; egrp<Mesh->nElemGroup; egrp++) */
/*       nn = max(nn, Mesh->ElemGroup[egrp].nNode); */
    
/*     // fvec, xvec, delta */
/*     ierr = xf_Error(xf_Alloc((void **) &fvec, nn, sizeof(int))); */
/*     if (ierr != xf_OK) return ierr; */
/*     ierr = xf_Error(xf_Alloc((void **) &xvec, nn*dim, sizeof(real))); */
/*     if (ierr != xf_OK) return ierr; */
/*     ierr = xf_Error(xf_Alloc((void **) &delta, nn*dim, sizeof(real))); */
/*     if (ierr != xf_OK) return ierr; */

/*     // loop over boundary groups */
/*     for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++){ */
      
/*       // are we curving this boundary?  Check input */
/*       for (i=0, CurveThis=xfe_False; i<nBFG; i++) */
/* 	if (strcmp(Mesh->BFaceGroup[ibfgrp].Title, BFGTitleArray[i]) == 0) CurveThis = xfe_True; */
/*       if (!CurveThis) continue; */

/*       xf_printf("%s ibfgrp = %d, title = %s\n",  */
/* 		(CurvingFlag) ? "Curving" : "Snapping",  */
/* 		ibfgrp, Mesh->BFaceGroup[ibfgrp].Title); */

/*       // loop over boundary faces */
/*       for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){ */

/* 	// pull off adjacent element group, number, face */
/* 	BFace = Mesh->BFaceGroup[ibfgrp].BFace+ibface; */

/* 	egrp = BFace->ElemGroup; */
/* 	elem = BFace->Elem; */
/* 	face = BFace->Face; */

/* 	EG = Mesh->ElemGroup + egrp; */

/* 	// get nodes on face */
/* 	ierr = xf_Error(xf_NodesOnFace(EG->QBasis, EG->QOrder, face, &nfnode, fvec)); */
/* 	if (ierr != xf_OK) return ierr; */

/* 	// fill in xvec */
/* 	for (i=0; i<nfnode; i++) */
/* 	  for (j=0; j<dim; j++) */
/* 	    xvec[i*dim+j] = Mesh->Coord[EG->Node[elem][fvec[i]]][j]; */

/* 	// project boundary nodes to geometry, get deltas */
/* 	ierr = xf_Error(xf_ProjectToGeom(Geom, -1, Mesh->BFaceGroup[ibfgrp].Title, */
/* 					 dim, nfnode, xvec)); */
/* 	if (ierr != xf_OK) return ierr; */

/* 	// fill in projection delta */
/* 	for (i=0; i<EG->nNode*dim; i++) delta[i] = 0.; */
/* 	for (i=0; i<nfnode; i++){ */
/* 	  for (j=0, node=fvec[i]; j<dim; j++) */
/* 	    delta[node*dim+j] = xvec[i*dim+j] - Mesh->Coord[EG->Node[elem][node]][j]; */
/* /\* 	  xf_printf("ibface=%d, egrp=%d, elem=%d, i=%d, node=%d, delta = %.5f %.5f\n", *\/ */
/* /\* 		    ibface, egrp, elem, i, node, delta[node*dim+0], delta[node*dim+1]); *\/ */
/* 	} */
	
/* 	// store deltas for Q1 nodes */
/* 	for (i=0; i<nfnode; i+=(nfnode-1)) */
/* 	  for (j=0, node=fvec[i]; j<dim; j++) */
/* 	    NodeDelta[EG->Node[elem][node]*dim+j] = delta[node*dim+j]; */

/* 	// propagate deltas from face into element interior (Q1 nodes not moved) */
/* 	ierr = xf_Error(xf_CurveElementInterior(Mesh, egrp, elem, face, delta)); */
/* 	if (ierr != xf_OK) return ierr; */

	
/*       } // ibface */
/*     } // ibfgrp */

/*     // apply NodeDelta to move Q1 nodes */
/*     for (i=0; i<Mesh->nNode; i++) */
/*       for (j=0; j<dim; j++) Mesh->Coord[i][j] += NodeDelta[i*dim+j]; */
    
/*     xf_Release( (void *) fvec); */
/*     xf_Release( (void *) xvec); */
/*     xf_Release( (void *) delta); */
/*     xf_Release( (void *) NodeDelta); */

/*   } // end if Geom != NULL */

  
/*   // release memory */
/*   xf_Release2( (void **) BFGTitleArray); */

/*   return xf_OK; */
/* } */



