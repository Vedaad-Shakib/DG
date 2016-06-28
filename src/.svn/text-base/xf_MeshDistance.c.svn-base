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
  FILE:  xf_MeshDistance.c

  Functions for calculating distance functions

*/


#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_Basis.h"
#include "xf_BasisFcn.h"
#include "xf_Memory.h"
#include "xf_Data.h"
#include "xf_Math.h"
#include "xf_Param.h"
#include "xf_MathLapack.h"
#include "xf_Quad.h"
#include "xf_MeshToolsStruct.h"
#include "xf_MeshTools.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_DataMath.h"
#include "xf_Residual.h"
#include "xf_EqnSetHook.h"


/******************************************************************/
//   FUNCTION Definition: xf_BCIsWall
static int 
xf_BCIsWall(xf_BC *BC, enum xfe_Bool *IsWall)
{
/*
PURPOSE: 

  Determines if boundary condition BC corresponds to a wall, e.g. for
  wall distance calculation purposes.

INPUTS: 
 
  BC : boundary condition in question

OUTPUTS:

  IsWall : true if BC is a wall
  
RETURNS: Error code

*/
  int ierr;

  // try calling an equation-specific function
  ierr = xf_EqnSetBCIsWall(BC, IsWall);
  // if no eqnset loaded, assume BC is a wall
  if (ierr == xf_DYNAMIC_LIBRARY_ERROR) (*IsWall) = xfe_True;
  // error out if some other error occured
  else if (ierr != xf_OK) return xf_Error(ierr);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_Dist2Segment
static int
xf_Dist2Segment( real *xseg, real *x0, real *pdist, real *xproj)
{
/* 
PURPOSE

  Computes distance to a segment in x,y space
  
INPUTS:

  xseg : coordinates of segment, unrolled
         [x0,y0, x1,y1]
  x0 : point from which we are computing the distance
  
OUTPUTS:

  (*pdist) : minimum distance to the segment
  xproj    : point on segment closest to x0 (optional)

RETURNS:  Error code

*/
  real n[2], xp[2];
  real NN;
  real d, d1, d2;

  // get line into form a*x + b*y + c = 0
  // (a,b) = normal vector components
  // c = normal dot (-x1)
  n[0] = -xseg[3] + xseg[1];
  n[1] =  xseg[2] - xseg[0];
  NN = sqrt(n[0]*n[0]+n[1]*n[1]);
  if (NN == 0.) return xf_Error(xf_INPUT_ERROR);
  n[0] /= NN; n[1] /= NN;

  // distance from line
  d = n[0]*(x0[0]-xseg[0]) + n[1]*(x0[1]-xseg[1]);
  d = fabs(d);

  // check distance to endpoint 1
  n[0] = xseg[2] - xseg[0];
  n[1] = xseg[3] - xseg[1];
  if (NN == 0.) return xf_Error(xf_INPUT_ERROR);
  NN = sqrt(n[0]*n[0]+n[1]*n[1]);
  n[0] /= NN; n[1] /= NN;
  d1 = n[0]*(x0[0]-xseg[0]) + n[1]*(x0[1]-xseg[1]);

  // check distance to endpoint 2
  n[0] = -n[0];
  n[1] = -n[1];
  d2 = n[0]*(x0[0]-xseg[2]) + n[1]*(x0[1]-xseg[3]);

  if (d1 < 0.){
    (*pdist) = xf_Distance(x0, xseg, 2);
    xp[0] = xseg[0]; xp[1] = xseg[1];
  }
  else if (d2 < 0.){
    (*pdist) = xf_Distance(x0, xseg+2, 2);
    xp[0] = xseg[2]; xp[1] = xseg[3];
  }
  else{
    (*pdist) = d;
    xp[0] = xseg[2] + n[0]*d2; xp[1] = xseg[3] + n[1]*d2; // note n is 2->1
  }

  if (xproj != NULL){
    xproj[0] = xp[0];
    xproj[1] = xp[1];
  }

  return xf_OK;

}


/******************************************************************/
//   FUNCTION Definition: xf_Dist2Segment3d
static int
xf_Dist2Segment3d( real *xseg, real *x0, real *pdist, real *xproj)
{
/* 
PURPOSE

  Computes distance to a segment in x,y,z space
  
INPUTS:

  xseg : coordinates of segment, unrolled
         [x0,y0,z0, x1,y1,z1]
  x0 : point from which we are computing the distance
  
OUTPUTS:

  (*pdist) : minimum distance to the segment
  xproj    : point on segment closest to x0 (optional)

RETURNS:  Error code

*/
  int i;
  real *x[2];
  real n[3], xp[3];
  real cp[3];
  real v0[3], v1[3];
  real NN;
  real dp0, dp1;

  // segment endpoints
  x[0] = xseg; x[1] = xseg+3;

  // vector pointing along segment, from node 0 to node 1
  for (i=0; i<3; i++) n[i] = x[1][i] - x[0][i];  
  NN = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if (NN == 0.) return xf_Error(xf_INPUT_ERROR);
  n[0] /= NN; n[1] /= NN; n[2] /= NN;
  
  // is point closest to node 0 or node 1?
  for (i=0; i<3; i++) v0[i] = x0[i] - x[0][i];
  for (i=0, dp0=0.; i<3; i++) dp0 += v0[i]*n[i];
  for (i=0; i<3; i++) v1[i] = x0[i] - x[1][i];
  for (i=0, dp1=0.; i<3; i++) dp1 -= v1[i]*n[i];
  if (dp0 <= 0.){
    (*pdist) = xf_Distance(x0, x[0], 3);
    for (i=0; i<3; i++) xp[i] = x[0][i];
  }
  else if (dp1 <= 0.){
    (*pdist) = xf_Distance(x0, x[1], 3);
    for (i=0; i<3; i++) xp[i] = x[1][i];
  }
  else{
    // point is closest to segment, have to project
    xf_CrossProduct(n, v0, cp);
    (*pdist) = sqrt(cp[0]*cp[0] + cp[1]*cp[1] +cp[2]*cp[2]);
    for (i=0; i<3; i++) xp[i] = x[0][i] + dp0*n[i];
  }
  
  if (xproj != NULL){
    xproj[0] = xp[0];
    xproj[1] = xp[1];
    xproj[2] = xp[2];
  }


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Dist2Triangle
static int
xf_Dist2Triangle( real *xtri, real *x0, real *pdist, real *xproj)
{
/* 
PURPOSE

  Computes distance to a triangle in x,y,z space
  
INPUTS:

  xtri : coordinates of triangle, unrolled
         [x0,y0,z0, x1,y1,z1, x2,y2,z2]
  x0 : point from which we are computing the distance
  
OUTPUTS:

  (*pdist) : minimum distance to the triangle
  xproj    : point on segment closest to x0 (optional)

RETURNS:  Error code

*/
  int ierr;
  int i, e;
  int i1, i2;
  real *x[3];
  real v1[3], v2[3], ve[3];
  real n[3], ne[3], xp[3];
  real xseg[6];
  real NN;
  real d, de;
  
  // node coordinates
  x[0] = xtri; x[1] = xtri+3; x[2] = xtri+6;
  
  // normal to plane of triangle
  for (i=0; i<3; i++) v1[i] = x[1][i] - x[0][i];
  for (i=0; i<3; i++) v2[i] = x[2][i] - x[0][i];
  xf_CrossProduct(v1, v2, n);
  NN = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if (NN == 0.) return xf_Error(xf_INPUT_ERROR);
  n[0] /= NN; n[1] /= NN; n[2] /= NN;
  
  // distance from plane
  for (i=0, d=0.; i<3; i++) d += n[i]*(x0[i]-x[0][i]);
  (*pdist) = fabs(d);

  // position on plane, if requested
  if (xproj != NULL) for (i=0; i<3; i++) xproj[i] = x0[i] - d*n[i];

  // check if projection is within triangle
  for (e=0; e<3; e++){ // loop over edges
    // indices of the two nodes on the edge
    i1 = (e+1)%3;
    i2 = (e+2)%3;
    // calculate ne = outward-pointing edge normal
    for (i=0; i<3; i++) ve[i] = x[i2][i] - x[i1][i];
    xf_CrossProduct(ve, n, ne);
    NN = sqrt(ne[0]*ne[0]+ne[1]*ne[1]+ne[2]*ne[2]);
    ne[0] /= NN; ne[1] /= NN; ne[2] /= NN;
    // check dot product to see on which side x0 is on
    for (i=0, de=0.; i<3; i++) de += ne[i]*(x0[i]-x[i1][i]);
    if (de >= 0.){
      // x0 is outside of triangle, take distance to edge segment
      for (i=0; i<3; i++) xseg[  i] = x[i1][i];
      for (i=0; i<3; i++) xseg[3+i] = x[i2][i];
      ierr = xf_Error(xf_Dist2Segment3d(xseg, x0, pdist, xproj));
      if (ierr != xf_OK) return ierr;
      break; // nothing else to check
    }
  }
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_Dist2Face
int
xf_Dist2Face(int dim, const real *x0in, enum xfe_BasisType FBasis,
	     int Order, int nf, real *X, real *pdist, real *xproj)
{
  int ierr;
  int i, k, dm1;
  int refine, nnode;
  int nsplit, isplit, minisplit;
  int *iv;
  int *vsplit = NULL;
  enum xfe_ShapeType FShape; 
  real mindist, dist;
  real x0[3], xp0[3];
  real xtri[9];
  real xseg[4];
  real xref[3], A[9], x01[3], x02[3], rhs[3], khat[3], iA[9], xi, eta;
  real *x = NULL, *xp = NULL;
  real *coord = NULL;
  xf_BasisData *PhiData = NULL;

  dm1 = dim-1; // face dimension

  if ((dm1 != 1) && (dm1 != 2)) return xf_Error(xf_NOT_SUPPORTED);

  if (nf <= 1) return xf_Error(xf_INPUT_ERROR);

  // initialize xproj and xp for projections
  if (xproj != NULL){
    for (i=0; i<dim; i++) xproj[i] = x0[i];
    xp = xp0;
  }
  else xp = NULL;

  // set x0 to input value (in case have to modify later)
  for (i=0; i<dim; i++) x0[i] = x0in[i];

  // initialize mindist
  mindist = 1e30;

  // get face shape
  ierr = xf_Error(xf_Basis2Shape(FBasis, &FShape));
  if (ierr != xf_OK) return ierr;


  // split up shape into segments or triangles
  refine = 2*Order+2; // refinement index
  if (Order == 1) refine = 1; // linear geometry does not need subdivision
  ierr = xf_Error(xf_GetRefineCoords(FShape, refine, &nnode, &coord,
				     &nsplit, &vsplit, NULL, NULL));
  if (ierr != xf_OK) return ierr;
  
  // evaluate basis functions at coord
  PhiData = NULL;
  ierr = xf_Error(xf_EvalBasis(FBasis, Order, xfe_True, nnode,
			       coord, xfb_Phi, &PhiData));
  if (ierr != xf_OK) return ierr;

  // allocate memory
  ierr = xf_Error(xf_ReAlloc((void **) &x, nnode*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // interpolate x(coord)
  xf_MxM_Set(PhiData->Phi, X, nnode, nf, dim,  x);

  // loop over subshapes = segments/triangles
  minisplit = -1;
  for (isplit=0; isplit<nsplit; isplit++){

    // calculate distance to subshape
    if (dm1 == 2){
      // distance to triangle
      iv = vsplit + 3*isplit;
      for (k=0; k<3; k++) 
	for (i=0; i<3; i++) xtri[k*3+i] = x[dim*iv[k]+i];
      ierr = xf_Error(xf_Dist2Triangle(xtri, x0, &dist, xp));
      if (ierr != xf_OK) return ierr;
    }
    else{
      // distance to line
      iv = vsplit + 2*isplit;
      for (k=0; k<2; k++) 
	for (i=0; i<2; i++) xseg[k*2+i] = x[dim*iv[k]+i];
      ierr = xf_Error(xf_Dist2Segment(xseg, x0, &dist, xp));
      if (ierr != xf_OK) return ierr;
    }

    // set mindist and xproj appropriately
    if (dist < mindist){
      minisplit = isplit;
      mindist = dist;
      if (xproj != NULL) for (k=0; k<dim; k++) xproj[k] = xp[k];
    }
  } // isplit
  
  // snap xp to true curved geometry
  if ((minisplit >= 0) && (xproj != NULL) && (Order > 1)){
    // determine xref based on current xproj and minisplit
    isplit = minisplit;
    if (dm1 == 2){ // dealing with a triangle
      iv = vsplit + 3*isplit;
      for (k=0; k<3; k++) for (i=0; i<3; i++) xtri[k*3+i] = x[dim*iv[k]+i];
      for (i=0; i<3; i++) x01[i] = xtri[1*3+i]-xtri[0*3+i];
      for (i=0; i<3; i++) x02[i] = xtri[2*3+i]-xtri[0*3+i];
      xf_CrossProduct(x01, x02, khat); // khat = perp to tri
      // build a 3x3 system to solve for xi, eta, delta
      // note, delta should be zero but we need a fully-determined system
      for (i=0; i<3; i++){A[3*i+0] = x01[i]; A[3*i+1] = x02[i]; A[3*i+2] = khat[i];}
      for (i=0; i<3; i++) rhs[i] = xproj[i] - xtri[0*3+i];
      ierr = xf_Error(xf_MatDetInv(A, 3, NULL, iA));
      if (ierr != xf_OK) return ierr;
      for (i=0, xi =0.; i<3; i++) xi  += iA[0*3+i]*rhs[i];
      for (i=0, eta=0.; i<3; i++) eta += iA[1*3+i]*rhs[i];
      for (i=0; i<2; i++) // convert (xi,eta) of isplit into ref coord on face
	xref[i] = coord[2*iv[0]+i]*(1.-xi-eta) + coord[2*iv[1]+i]*xi + coord[2*iv[2]+i]*eta;
    }
    else{ // dealing with a line
      iv = vsplit + 2*isplit;
      for (k=0; k<2; k++) for (i=0; i<2; i++) xseg[k*2+i] = x[dim*iv[k]+i];
      xi = xf_Distance(xseg+2*0, xproj, 2)/xf_Distance(xseg+2*0, xseg+2*1, 2);
      xref[0] = coord[iv[0]]*(1.-xi) + coord[iv[1]]*xi;
    }
    
    // evaluate basis functions of face interpolation
    ierr = xf_Error(xf_EvalBasis(FBasis, Order, xfe_True, 1,
				 xref, xfb_Phi, &PhiData));
    if (ierr != xf_OK) return ierr;

    // interpolate face coords -> xproj
    xf_MxM_Set(PhiData->Phi, X, 1, nf, dim, xproj);

  } // end snapping to true geom


  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  // Release memory
  xf_Release((void *) coord);
  xf_Release((void *) vsplit);
  xf_Release((void *) x);

  // return minimum distance
  (*pdist) = mindist;


  return xf_OK;

}



/******************************************************************/
//   FUNCTION Definition: xf_CalcDistFcnElems
static int 
xf_CalculateDistFcnElems(xf_All *All, int nface, int *IData, real *XFace, 
			 int nf0, int nneigh, xf_Vector *DistFcn)
{
/*
PURPOSE: 

  Within each non-halo element, calculates distances to the set of
  faces in IData/XFace.  The distance function on each element is then
  defined as the minimum distance to one of the faces.

  Some heuristics are employed in this calculation, and the distance
  may not be the minimum distance if some very highly nonlinear,
  high-order geometry elements are present.  Anisotropic/skewed
  elements are fine.

INPUTS: 
 
  All: All structure
  nface : number of faces to consider in distance calculation
  IData : integer list containing FBasis, Order, and nf for each face
          where nf = # Lagrange nodes to interpolate the face
  XFace : coordinates of each of the Lagrange nodes for each face
  nf0 : max number of nodes per face (for indexing of XFace)
  nneigh : max number of neighbors around one node

OUTPUTS:

  DistFcn : distance function on each element

RETURNS: Error code

*/

  int ierr, dim, i;
  int iegrp, ielem;
  int ipoint, npoint;
  int iface, nf;
  int Order;
  int ilarge, ineigh;
  int *ivec = NULL;
  enum xfe_Bool BFlag, changed;
  enum xfe_ShapeType Shape, FShape;
  enum xfe_BasisType FBasis;
  real mindist, dist, dmin;
  real *xpoint = NULL;
  real *xref  = NULL, *xglob = NULL;
  real *X = NULL;
  real *dvec = NULL;
  xf_BasisData *GPhiData = NULL; // for ref2glob transformation
  xf_Mesh *Mesh;

  Mesh  = All->Mesh;
  dim   = Mesh->Dim;

  nneigh = min(2*nneigh, nface); // overestimate of number of neighbor faces per Q1 node

  // allocate memory
  ierr = xf_Error(xf_Alloc((void **) &dvec, nneigh, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc((void **) &ivec, nneigh, sizeof(int)));
  if (ierr != xf_OK) return ierr;


  /*-------------------------*/
  /* Main loop over elements */
  /*-------------------------*/

  // loop over all element groups
  for (iegrp=0; iegrp<Mesh->nElemGroup; iegrp++){

    // reference-space coords of Lagrange nodes
    ierr = xf_Error(xf_LagrangeNodes(DistFcn->Basis[iegrp], DistFcn->Order[iegrp], &npoint, NULL, &xref));
    if (ierr != xf_OK) return ierr;

    // reallocate xglob to store global coords
    ierr = xf_Error(xf_ReAlloc((void **) &xglob, dim*npoint, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    // loop over elements
    for (ielem=0; ielem<Mesh->ElemGroup[iegrp].nElem; ielem++){
      changed = xfe_False;

      // ref2glob -> global coordinates of the points
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, iegrp, ielem, &GPhiData, xfe_True,
				      npoint, xref, xglob));
      if (ierr != xf_OK) return ierr;

      // loop over points
      for (ipoint=0; ipoint<npoint; ipoint++){

	// point coordinate
	xpoint = xglob + dim*ipoint;

	// determine closest nneigh faces based on node coordinates
	ilarge = 0;
	for (i=0; i<nneigh; i++) dvec[i] = 1e30;
	for (i=0; i<nneigh; i++) ivec[i] = -1;
      
	for (iface=0; iface<nface; iface++){
	  nf     = IData[3*iface+2];
	  X      = XFace + iface*nf0*dim;
	  dist   = xf_Distance(xpoint, X+0, dim);
	  for (i=1; i<nf; i++) dist = min(dist, xf_Distance(xpoint, X+dim*i, dim));
	  if (dist < dvec[ilarge]){
	    dvec[ilarge] = dist;
	    ivec[ilarge] = iface;
	    for (i=1, ilarge=0; i<nneigh; i++) // reset ilarge to max dist index
	      if (dvec[i] > dvec[ilarge]) ilarge = i;
	  }
	}
	
	if (ivec[ilarge] < 0) return xf_Error(xf_CODE_LOGIC_ERROR);

	// minimum distance to nodes
	dmin = dvec[0];
	for (i=1; i<nneigh; i++) dmin = min(dmin, dvec[i]);

	mindist = dmin; // start off with closest node distance
	
	// loop over closest faces to get more accurate distances
	for (ineigh=0; ineigh<nneigh; ineigh++){

	  iface = ivec[ineigh];

	  if (iface < 0) return xf_Error(xf_CODE_LOGIC_ERROR);

	  FBasis = IData[3*iface+0];
	  Order  = IData[3*iface+1];
	  nf     = IData[3*iface+2];
	  
	  X      = XFace + iface*nf0*dim;
	  
	  // calculate min distance to face, including dv
	  ierr = xf_Error(xf_Dist2Face(dim, xpoint, FBasis, Order, nf, X, &dist, NULL));
	  if (ierr != xf_OK) return ierr;

	  // set minimum distance
	  mindist = min(mindist, dist);
	  
	} // fi

	// check if mindist is smaller than current distance
	if (mindist < DistFcn->GenArray[iegrp].rValue[ielem][ipoint])
	  DistFcn->GenArray[iegrp].rValue[ielem][ipoint] = mindist;

      } // ipoint

    } // ielem
  } // iegrp


  // destroy GPhiData
  ierr = xf_Error(xf_DestroyBasisData(GPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release ( (void  *) xref);
  xf_Release ( (void  *) xglob);
  xf_Release ( (void  *) dvec);
  xf_Release ( (void  *) ivec);

 
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MaxBFacesPerNode
static int 
xf_MaxBFacesPerNode(xf_All *All, int *pnfmax)
{
  int ierr, i, j;
  int nNode;
  int ibfgrp, ibface;
  int egrp, elem, face;
  int nfnode;
  int nmax;
  int fvec[xf_MAXQ1FACENODE];
  int *Node2nFace = NULL;
  int *Node;
  xf_BFace BFace;
  xf_ElemGroup *EG;
  xf_Mesh *Mesh;

  Mesh   = All->Mesh;
  nNode  = Mesh->nNode;

  // allocate an integer counter at nodes, of adjacent faces
  ierr = xf_Error(xf_Alloc( (void **) &Node2nFace, nNode, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (i=0; i<nNode; i++) Node2nFace[i] = 0;

  // loop over boundary faces
  for (ibfgrp=0; ibfgrp<Mesh->nBFaceGroup; ibfgrp++)
    for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
      BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
      egrp = BFace.ElemGroup;
      elem = BFace.Elem;
      face = BFace.Face;
            
      // local Q1 nodes on face of elem
      EG = Mesh->ElemGroup + egrp;
      ierr = xf_Error(xf_Q1NodesOnFace(EG->QBasis, EG->QOrder, face, &nfnode, fvec));
      if (ierr != xf_OK) return ierr;
    
      Node = Mesh->ElemGroup[egrp].Node[elem];
      for (j=0; j<nfnode; j++) Node2nFace[Node[fvec[j]]] += 1;
    } // ibface

  nmax = -1;
  for (i=0; i<nNode; i++) nmax = max(nmax, Node2nFace[i]);
  
  (*pnfmax) = nmax;
  
  xf_Release ( (void  *) Node2nFace);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CalculateDistFcn
int 
xf_CalculateDistFcn(xf_All *All)
{
  int ierr, i, d, dim;
  int myRank, nProc;
  int negrp, egrp;
  int elem, face;
  int DistFcnOrder;
  int ibface, negrphalo;
  int nface, nfacetot;
  int ibfgrp, nbfgrp;
  int iface, nf, nf0;
  int nfnode;
  int nneigh;
  int nWall;
  int *OrderVec = NULL;
  int *IData = NULL, *IDataGlob = NULL;
  int *fv = NULL;
  int *nFace = NULL;
  int *Node;
  int *ibuf = NULL;

  char DistFcnWallBoundaries[xf_MAXSTRLEN];
  char **WallNames = NULL;
  enum xfe_ShapeType Shape, FShape;
  enum xfe_BasisType FBasis;
  enum xfe_BasisType *BasisVec = NULL;
  enum xfe_Bool ParallelFlag;
  enum xfe_Bool IsWall;
  enum xfe_Bool Found, done;
  enum xfe_Bool EqnSetValid;
  enum xfe_Bool WriteDistFcn;
  enum xfe_Bool OverrideFlag;

  real *XFace = NULL, *XFaceGlob = NULL;
  real *rbuf = NULL;
  real *X = NULL;
  xf_BFace   BFace;
  xf_Face   *FaceList;
  xf_Data   *D;
  xf_Vector *DistFcn;
  xf_Mesh   *Mesh;
  xf_EqnSet *EqnSet;
  xf_ElemGroup *EG;

  Mesh   = All->Mesh;
  EqnSet = All->EqnSet;

  dim    = Mesh->Dim;

  xf_printf("Calculating the distance function ... "); fflush(stdout);

  // obtain number of processors
  ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
  if (ierr != xf_OK) return ierr;

  ParallelFlag = (nProc > 1);

  // get order of the distance function
  ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "DistFcnOrder", &DistFcnOrder));
  if (ierr != xf_OK) return ierr;

  // Check if DistFcn should be made writable
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "WriteDistFcn", &WriteDistFcn));
  if (ierr != xf_OK) return ierr;

  // allocate Basis and Order vectors for vector search
  negrp = Mesh->nElemGroup;
  negrphalo = (ParallelFlag ? 2*negrp : negrp);
  ierr = xf_Error(xf_Alloc( (void **) &BasisVec, negrphalo, sizeof(enum xfe_BasisType)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &OrderVec, negrphalo, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // request Lagrange basis of order DistFcnOrder
  for (egrp=0; egrp<negrphalo; egrp++){
    ierr = xf_Error(xf_Basis2UniformLagrange(Mesh->ElemGroup[egrp].QBasis, BasisVec+egrp));
    if (ierr != xf_OK) return ierr;
    OrderVec[egrp] = DistFcnOrder;
  }

  
  // look for DistFcn vector (do not make parallel)
  ierr = xf_Error(xf_FindVector(All, "WallDistance", xfe_LinkageGlobElem, 1, 
				NULL, 0, 0, BasisVec, OrderVec, NULL, NULL, NULL, 
				xfe_SizeReal, xfe_False, xfe_True, &D, &DistFcn, &Found));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = WriteDistFcn;

  xf_Release( (void  *) BasisVec);
  xf_Release( (void  *) OrderVec);

  // consider returning immediately if DistFcn was found

  // set distances to a very large number (will take minimum)
  ierr = xf_Error(xf_SetConstVector(DistFcn, 0, 1e30));
  if (ierr != xf_OK) return ierr;

  // number of bface groups
  nbfgrp = Mesh->nBFaceGroup;


  // allocate a FaceList vector for storing faces to which we need distances
  nface=0;
  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++) nface += Mesh->BFaceGroup[ibfgrp].nBFace;
  ierr = xf_Error(xf_Alloc((void **) &FaceList, nface, sizeof(xf_Face)));
  if (ierr != xf_OK) return ierr;



  /*------------------------------*/
  /*  Build a list of wall faces  */
  /*------------------------------*/

  // get list of any wall boundaries
  ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "DistFcnWallBoundaries", DistFcnWallBoundaries));
  if (ierr != xf_OK) return ierr;

  // Check if "DistFcnWallBoundaries" is not null -- use this first
  OverrideFlag = xfe_False;
  if (xf_NotNull(DistFcnWallBoundaries)){
    ierr = xf_Error(xf_ScanXStringAlloc(DistFcnWallBoundaries, xf_MAXSTRLEN, &nWall, &WallNames));
    if (ierr != xf_OK) return ierr;
    
    if (nWall > 0) OverrideFlag = xfe_True;
  }

  nface = 0;
  // if EqnSet is not initialized, set EqnSetValid=False
  EqnSetValid = (EqnSet->BCs != NULL);
  if (EqnSetValid){
    ierr = xf_Error(xf_SortEqnSetBCs(Mesh, EqnSet->BCs+0));
    if (ierr != xf_OK) return ierr;
  }

  for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){

    if (OverrideFlag){          // first check if overriding eqnset specification
      IsWall = xfe_False;
      for (i=0; i<nWall; i++){
	if (strcmp(WallNames[i], Mesh->BFaceGroup[ibfgrp].Title) == 0){
	  IsWall = xfe_True;
	  break;
	}
      } // i
    }
    else if (EqnSetValid){  // use eqnset if valid
      ierr = xf_Error(xf_BCIsWall(EqnSet->BCs[0].BC+ibfgrp, &IsWall));
      if (ierr != xf_OK) return ierr;
    }
    else IsWall = xfe_True;

    if (IsWall){ // add to FaceList
      for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){
	FaceList[nface+ibface].Group  = ibfgrp;
	FaceList[nface+ibface].Number = ibface;
      }
      nface += Mesh->BFaceGroup[ibfgrp].nBFace;
    }
  } // ibfgrp

  xf_Release2( (void **) WallNames);

  // error if no faces flagged (i.e. no walls)
  nfacetot = nface;
  ierr = xf_Error(xf_MPI_Allreduce(&nfacetot, 1, xfe_SizeInt, xfe_MPI_SUM));
  if (ierr != xf_OK) return ierr;

  if (nfacetot <= 0){
    xf_printf("Error, no wall boundary faces identified in DistFcn calculation.\n");
    return xf_Error(xf_NOT_FOUND);
  }


  // create vector containing face information
  ierr = xf_Error(xf_Alloc( (void **) &IData, nface*3, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  // determine max number of Lagrange nodes necessary to describe the faces
  nf0 = -1;
  for (iface=0; iface<nface; iface++){
    ibfgrp = FaceList[iface].Group;
    ibface = FaceList[iface].Number;

    if (ibfgrp < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
    egrp = BFace.ElemGroup;
    elem = BFace.Elem;
    face = BFace.Face;

    EG = Mesh->ElemGroup + egrp;

    // get element shape
    ierr = xf_Error(xf_Basis2Shape(EG->QBasis, &Shape));
    if (ierr != xf_OK) return ierr;
    
    // get face shape
    ierr = xf_Error(xf_FaceShape(Shape, face, &FShape));
    if (ierr != xf_OK) return ierr;
    
    // get default Lagrange basis for face
    ierr = xf_Error(xf_Shape2UniformLagrange(FShape, &FBasis));
    if (ierr != xf_OK) return ierr;
    
    // number of Lagrange nodes on face, nf
    ierr = xf_Error(xf_Order2nNode(FBasis, EG->QOrder, &nf));
    if (ierr != xf_OK) return ierr;

    IData[3*iface + 0] = FBasis;
    IData[3*iface + 1] = EG->QOrder;
    IData[3*iface + 2] = nf;

    nf0 = max(nf0, nf);
  }

  // reduce-max number of nodes per face
  ierr = xf_Error(xf_MPI_Allreduce(&nf0, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  if (nf0 <= 0) return xf_Error(xf_CODE_LOGIC_ERROR);

  // construct vector of node coordinates
  ierr = xf_Error(xf_Alloc( (void **) &XFace, nface*nf0*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_Alloc( (void **) &fv, nf0, sizeof(int)));
  if (ierr != xf_OK) return ierr;

  for (iface=0; iface<nface; iface++){
    ibfgrp = FaceList[iface].Group;
    ibface = FaceList[iface].Number;

    if (ibfgrp < 0) return xf_Error(xf_CODE_LOGIC_ERROR);
    
    BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
    egrp = BFace.ElemGroup;
    elem = BFace.Elem;
    face = BFace.Face;

    EG = Mesh->ElemGroup + egrp;

    // take face nodes from self elem
    ierr = xf_Error(xf_NodesOnFace(EG->QBasis, EG->QOrder, face, &nfnode, fv));
    if (ierr != xf_OK) return ierr;
    if (nfnode != IData[3*iface+2]) return xf_Error(xf_CODE_LOGIC_ERROR);

    // set XFace using fv
    X = XFace + iface*nf0*dim;
    for (i=0; i<nf0; i++) X[i] = 0.; // to avoid sending uninitialized data
    Node = EG->Node[elem];
    for (i=0; i<nfnode; i++)
      for (d=0; d<dim; d++)
	X[dim*i + d] = Mesh->Coord[Node[fv[i]]][d];
  }


  // release memory
  xf_Release( (void *) fv);
  xf_Release( (void *) FaceList);

  
  if (ParallelFlag){  // parallel specific : each proc gets all the bfaces

    // allocate space for consolidated versions of IData and XFace
    ierr = xf_Error(xf_Alloc( (void **) &IDataGlob, nfacetot*3, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &XFaceGlob, nfacetot*nf0*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    // gather number of faces vector to root
    if (myRank == 0){
      ierr = xf_Error(xf_Alloc( (void **) &nFace, nProc, sizeof(int)));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_MPI_Gather((void *) &nface, (void *) nFace, 1*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    
    // send IData and XFace to root
    ibuf = (myRank == 0) ? IDataGlob : NULL;
    ierr = xf_Error(xf_MPI_Gatherv((void *) IData, nface, (void *) ibuf, nFace, 3*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    rbuf = (myRank == 0) ? XFaceGlob : NULL;
    ierr = xf_Error(xf_MPI_Gatherv((void *) XFace, nface, (void *) rbuf, nFace, nf0*dim*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;
    

    // Broadcast data to all procs
    ierr = xf_Error(xf_MPI_Bcast( (void *) IDataGlob, nfacetot*3*sizeof(int), 0));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_MPI_Bcast( (void *) XFaceGlob, nfacetot*nf0*dim*sizeof(real), 0));
    if (ierr != xf_OK) return ierr;

    // destroy previous local versions
    xf_Release( (void *) IData);
    xf_Release( (void *) XFace);
    xf_Release( (void *) nFace);

    IData = IDataGlob;
    XFace = XFaceGlob;
    
  } // end if parallel 

  nface = nfacetot;

  // maximum number of bfaces per Q1 node
  ierr = xf_Error(xf_MaxBFacesPerNode(All, &nneigh));
  if (ierr != xf_OK) return ierr;
  
  // reduce-max nneigh
  ierr = xf_Error(xf_MPI_Allreduce(&nneigh, 1, xfe_SizeInt, xfe_MPI_MAX));
  if (ierr != xf_OK) return ierr;
  if (nneigh <= 0) return xf_Error(xf_CODE_LOGIC_ERROR);

  // calculate distances for each elem to faces in IData,XFace
  ierr = xf_Error(xf_CalculateDistFcnElems(All, nface, IData, XFace, nf0, nneigh, DistFcn));
  if (ierr != xf_OK) return ierr;


  xf_printf("done.\n");

  // release memory
  xf_Release( (void *) IData);
  xf_Release( (void *) XFace);

  return xf_OK;
}








/*----------------------*/
/*   LEGACY FUNCTIONS   */
/*----------------------*/





/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_FaceElementHalo */
/* static int  */
/* xf_FaceElementHalo(xf_Mesh *Mesh, int ibfgrp, int ibface,  */
/* 		   enum xfe_Bool HaloFlag, int *egrp, int *elem, int *face) */
/* { */
/* /\* */
/* PURPOSE:  */

/*   Determines halo element adjacent to face (ibfgrp, iiface), if */
/*   HaloFlag is true.  If HaloFlag is false, returns the non-halo */
/*   element.  There must be one halo and one non-halo element adjacent */
/*   to the face.  The only exception is if a boundary face is passed in, */
/*   in which case "-1" is returned if a halo element is desired. */

/* INPUTS:  */
 
/*   Mesh : mesh structure */
/*   ibfgrp, ibface : group and number of face in question */
/*   HaloFlag : True if want halo element */

/* OUTPUTS: */

/*   egrp, elem, face : element in question */

/* RETURNS: Error code */

/* *\/ */
/*   int k, negrp; */
/*   int egrpS, elemS, faceS; */
/*   int egrpH, elemH, faceH; */
  
/*   negrp = Mesh->nElemGroup; */

/*   // pull off self (S) and halo (H) elements */
/*   xf_FaceElements(Mesh, ibfgrp, ibface, &egrpS, &elemS, &faceS,  */
/* 		  &egrpH, &elemH, &faceH); */

/*   if (ibfgrp >= 0){ // boundary face */
/*     if (HaloFlag){ */
/*       if (egrp != NULL) (*egrp) = -1; */
/*       if (elem != NULL) (*elem) = -1; */
/*       if (face != NULL) (*face) = -1; */
/*     } */
/*     else{ */
/*       if (egrp != NULL) (*egrp) = egrpS; */
/*       if (elem != NULL) (*elem) = elemS; */
/*       if (face != NULL) (*face) = faceS; */
/*     } */
/*   } */
/*   else{ // interior face */
    
/*     // Make sure one is halo, store in H */
/*     if (egrpS >= negrp){ */
/*       swap(egrpS, egrpH, k); */
/*       swap(elemS, elemH, k); */
/*       swap(faceS, faceH, k); */
/*     } */
/*     if (egrpH < negrp) return xf_Error(xf_PARALLEL_ERROR); // one should be a halo */

/*     if (HaloFlag){ */
/*       if (egrp != NULL) (*egrp) = egrpH; */
/*       if (elem != NULL) (*elem) = elemH; */
/*       if (face != NULL) (*face) = faceH; */
/*     } */
/*     else{ */
/*       if (egrp != NULL) (*egrp) = egrpS; */
/*       if (elem != NULL) (*elem) = elemS; */
/*       if (face != NULL) (*face) = faceS; */
/*     } */
    
/*   }   */

/*   return xf_OK; */
/* } */



/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_Transpose */
/* static void */
/* xf_Transpose(const real *A, int rA, int cA, real *B) */
/* { */

/* /\* */
/* PURPOSE:  */

/*   Transposes A and stores the result in B */

/* INPUTS:  */
 
/*   A : matrix to be transposed */
/*   rA : number of rows in A */
/*   cA : number of columns in A */
 
/* OUTPUTS: */
 
/*   B : cA by rA matrix = transpose(A) */
 
/* RETURNS: None */

/* *\/ */

/*   int i, j; */

/*   for (i=0; i<rA; i++) */
/*     for (j=0; j<cA; j++) */
/*       B[j*rA+i] = A[i*cA+j]; */
/* } */



/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_Dist2FaceResidual */
/* static int  */
/* xf_Dist2FaceResidual(int dim, const real *x0, int nf, real *X, real *D, */
/* 		     xf_BasisData *PhiData, real *R, real *R_xi, real *dist) */
/* { */
/* /\* */
/* PURPOSE:  */

/*   Calculates residual and linearization and distance to face */
/*   corresponding to the possibly nonlinear minimization problem */
/*   specified in xf_Dist2Face.   */

/* INPUTS:  */
 
/*   dim  : dimension of x0, X */
/*   x0   : point from which to compute distance */
/*   nf : number of Lagrange nodes interpolating face */
/*   X  : coords of face Lagrange points */
/*   D  : distance function values at Lagrange points */
/*   PhiData : basis function values, gradients, Hessians at xi */
  
/* OUTPUTS: */

/*   R : residual for distance minimization; size dim(xi) */
/*   R_xi : linearization of R w.r.t xi;     size dim(xi) x dim(xi) */
/*   (*dist)  : distance between x0 and x(xi) */
 
/* RETURNS: Error code */

/* *\/ */
/*   int i, j, k; */
/*   int fdim, fdim2; */
/*   real x[3]; */
/*   real d; */
/*   real r[3]; */
/*   real rmag; */
/*   real rmag_xi[3]; */
/*   real T[27]; */
/*   real x_xi[9]; */
/*   real d_xi[3]; */
/*   real x_xi2[27]; */
/*   real d_xi2[9]; */
/*   real t, r2; */

/*   // input checks */
/*   if (nf != PhiData->nn) return xf_Error(xf_INPUT_ERROR); */
/*   if ((fdim=PhiData->dim) > 2) return xf_Error(xf_INPUT_ERROR); */

/*   // interpolate x(xi) and d(xi) */
/*   xf_MxM_Set(PhiData->Phi, X, 1, nf, dim,  x); */
/*   xf_MxM_Set(PhiData->Phi, D, 1, nf,   1, &d); */

/*   // calculate r = x(xi) - x0, and mag of r */
/*   for (i=0; i<dim; i++) r[i] = x[i] - x0[i]; */
/*   for (i=0, rmag=0; i<dim; i++) rmag += r[i]*r[i]; */
/*   rmag = sqrt(rmag); */
  
/*   /\*  xf_printf("  x = %6E %.6E %.6E, x0 = %.6E %.6E %.6E\n", x[0], x[1], x[2], x0[0], x0[1], x0[2]); *\/ */
/*   /\*  xf_printf("  r = %6E %.6E %.6E, rmag = %.10E\n", r[0], r[1], r[2], rmag); *\/ */

/*   // distance = rmag + d */
/*   if (dist != NULL) (*dist) = rmag+d; */

/*   if (rmag == 0.0){ // indicates x0 is on face, at xi */
/*     for (i=0; i<fdim; i++) R[i] = 0.0; */
/*     for (i=0; i<fdim*fdim; i++) R_xi[i] = 0.0; */
/*     for (i=0; i<fdim; i++) R_xi[i*(fdim+1)] = 1.0; // R_xi set to I */
/*     return xf_OK; */
/*   } */
  

/*   // interpolate x_xi(xi) and d_xi(xi) */
/*   xf_MxM_Set(PhiData->GPhi, X, fdim, nf, dim, T); */
/*   xf_Transpose(T, fdim, dim, x_xi); */
/*   xf_MxM_Set(PhiData->GPhi, D, fdim, nf, 1, d_xi); */

/*   /\*  for (i=0; i<nf; i++) *\/ */
/*   /\*     for (j=0; j<fdim; j++) *\/ */
/*   /\*       xf_printf("X[%d,%d] = %.6E\n", i, j, X[i*fdim+j]); *\/ */
/*   /\*   for (j=0; j<fdim; j++) *\/ */
/*   /\*     for (i=0; i<nf; i++) *\/ */
/*   /\*       xf_printf("GPhi[%d, %d] = %.6E\n", j, i, PhiData->GPhi[j*nf+i]); *\/ */
  
/*   // construct rmag_xi{i} = 1/rmag * r{k} * r_xi{k,i} */
/*   // note, r_xi{k,i} = x_xi{k,i} */
/*   for (i=0; i<fdim; i++){ */
/*     rmag_xi[i] = 0.; */
/*     for (k=0; k<dim; k++) */
/*       rmag_xi[i] += 1/rmag * r[k]*x_xi[k*fdim+i]; */
/*   } */

/*   // construct R{i} = rmag_xi{i} + d_xi{i} */
/*   for (i=0; i<fdim; i++) R[i] = rmag_xi[i] + d_xi[i]; */
  
/*   // interpolate x_xi2(xi) and d_xi2 */
/*   xf_MxM_Set(PhiData->HPhi, X, fdim*fdim, nf, dim, T); */
/*   xf_Transpose(T, fdim*fdim, dim, x_xi2); */
/*   xf_MxM_Set(PhiData->HPhi, D, fdim*fdim, nf, 1, d_xi2); */


/*   /\* construct: */

/*      R_xi{i,j} = - 1/rmag^2 *rmag_xi{j} * r_xi{k,i} */
/*                  + 1/rmag * r_xi{k,j} * r_xi{k,i} */
/* 		 + 1/rmag * r{k} * r_xi2{k,i,j} */
/* 		 + d_xi2{i,j} */
/*      Note,  r_xi2 = x_xi2 */
/*   *\/ */
/*   r2 = rmag*rmag; */
/*   fdim2 = fdim*fdim; */
/*   for (i=0; i<fdim; i++){ */
/*     for (j=0; j<fdim; j++){ */
/*       t = 0.; */
/*       for (k=0; k<dim; k++){ */
/* 	t -= rmag_xi[j]/r2 * r[k]*x_xi[k*fdim+i]; */
/* 	t += 1/rmag * x_xi[k*fdim+i]*x_xi[k*fdim+j]; */
/* 	t +=  r[k]/rmag * x_xi2[k*fdim2 + i*fdim+j]; */
/*       } */
/*       t += d_xi2[i*fdim+j]; */
/*       R_xi[i*fdim+j] = t; */
/*     } // j */
/*   } // i */

/*   return xf_OK; */
/* } */


/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_Dist2FaceNewton */
/* static int  */
/* xf_Dist2FaceNewton(int dim, const real *x0in, enum xfe_BasisType FBasis,  */
/* 	     int Order, int nf, real *X, real *D, real *pdist) */
/* { */
/* /\* */
/* PURPOSE:  */

/*   Calculates minimum: */

/*     (distance between x0 and a Lagrange-interpolated face) */
/*   + (distance function interpolated at face) */

/*   i.e. */

/*      (*pdist) = min_(xi{i}) [(r{k} * r{k})^.5 + d] */
  
/*   where  r{k} = x{k} - x0{k} */
/*          x{k} = X{n,k}*Phi{n}(xi{i}) */
/*          d    = D{n}*Phi{n}(xi{i}) */

/*   Note, the face geometry and the current distance function are */
/*   interpolated using the same Lagrange basis and order. */

/*   Also, the Q1 nodes of the face are not considered, as they are */
/*   assumed to be taken care of outside this function. */
  

/* INPUTS:  */
 
/*   dim  : dimension of x0 */
/*   x0   : point from which to compute distance */
/*   FBasis : Lagrange basis of face */
/*   Order : Order with which face is interpolated */
/*   nf : number of Lagrange nodes interpolating face */
/*   X  : coords of face Lagrange points */
/*   D  : distance function values at Lagrange points */
  
/* OUTPUTS: */

/*   (*pdist)  : calculated minimum distance, as described above */
 
/* RETURNS: Error code */

/* *\/ */
/*   int ierr; */
/*   int i, d, dm1; */
/*   int iedge, nedge, ne; */
/*   int *ev = NULL; */
/*   enum xfe_Bool converged, inside, singular; */
/*   enum xfe_ShapeType FShape, EShape; */
/*   enum xfe_BasisType EBasis; */
/*   int isafety; */
/*   int iNewton; */
/*   int MaxNewton = 20;    // maximum number of Newton iterations */
/*   real MaxStep = 0.1;   // maximum step size in ref space */
/*   real omega;            // Newton under-relax factor */
/*   real tol = 1e-10;      // tolerance on Newton */
/*   real tolinside = 1e-5; // tolerance on inside ref elem */
/*   real mindist, dist; */
/*   real dximag, det; */
/*   real Rnorm, Rnorm0; */
/*   real xi[2], dxi[2]; */
/*   real R[2];  */
/*   real R_xi[4], iR_xi[4]; */
/*   real x0[3], xmin[3], xmax[3]; */
/*   real hmax; */
/*   real *XE = NULL, *DE = NULL; */
/*   xf_BasisData *PhiData; */

/*   dm1 = dim-1; // face dimension */

/*   if ((dm1 != 1) && (dm1 != 2)) return xf_Error(xf_NOT_SUPPORTED); */

/*   if (nf <= 1) return xf_Error(xf_INPUT_ERROR); */

/*   // set x0 to input value (in case have to modify later) */
/*   for (d=0; d<dim; d++) x0[d] = x0in[d]; */

/*   // initialize mindist */
/*   mindist = 1e30; */

/*   // get face shape */
/*   ierr = xf_Error(xf_Basis2Shape(FBasis, &FShape)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // set starting xi in (approximate) center of face */
/*   switch (FShape){ */
/*   case xfe_Segment: */
/*     xi[0] = 0.5; */
/*     break; */
/*   case xfe_Triangle: */
/*     xi[0] = 0.3; xi[1] = 0.3; */
/*     break; */
/*   case xfe_Quadrilateral: */
/*     xi[0] = 0.5; xi[1] = 0.5; */
/*     break; */
/*   default: */
/*     return xf_Error(xf_NOT_SUPPORTED); */
/*     break; */
/*   } */

/*   /\*   xf_printf("In Dist2Face.\n"); *\/ */
/*   /\*   xf_printf("x0 = %.6E %.6E\n", x0[0], x0[1]); *\/ */
/*   /\*   for (i=0; i<nf; i++) *\/ */
/*   /\*     for (d=0; d<dim; d++) *\/ */
/*   /\*       xf_printf("X[%d,%d] = %.6E\n", i, d, X[i*dim+d]); *\/ */
  
/*   /\*-------------*\/ */
/*   /\* Face Newton *\/ */
/*   /\*-------------*\/ */

/*   // allocate and initialize variables before loop */
/*   PhiData   = NULL; */
/*   converged = xfe_False; */
/*   dist      = mindist; */

/*   // Begin Newton-Raphson iteration */
/*   isafety  = 0; */
/*   singular = xfe_True; */
/*   while ((singular) && (isafety++ < 3)){ */
    
/*     singular = xfe_False; */
    
/*     for (iNewton=0; iNewton<MaxNewton; iNewton++){ */

/*       // compute basis shape, grad, hess at xi */
/*       ierr = xf_Error(xf_EvalBasis(FBasis, Order, xfe_True, 1, xi, xfb_All, &PhiData)); */
/*       if (ierr != xf_OK) return ierr; */
      
/*       // compute residual, R{i}, linearization, R{i}_xi{j}, and distance, dist */
/*       ierr = xf_Error(xf_Dist2FaceResidual(dim, x0, nf, X, D, PhiData, R, R_xi, &dist)); */
/*       if (ierr != xf_OK) return ierr; */
      
/*       /\*   if ((DEBUGFLAG) || ((iNewton > 10) && (xi[0] > 0.0) && (xi[0] < 1.0)) ){ *\/ */
/*       /\* 	xf_printf("iNewton=%d, xi = %.3E %.3E, R = %.3E %.3E, dist = %.3E\n", *\/ */
/*       /\* 		  iNewton, xi[0], xi[1], R[0], R[1], dist); *\/ */
/*       /\* 	xf_printf("  R_xi = %.3E %.3E %.3E %.3E\n", R_xi[0], R_xi[1], R_xi[2], R_xi[3]); *\/ */
/*       /\*       } *\/ */
      
/*       // convergence check */
/*       for (d=0, Rnorm=0.; d<dm1; d++) Rnorm += R[d]*R[d]; */
/*       Rnorm = sqrt(Rnorm); */
/*       if (iNewton == 0) Rnorm0 = Rnorm; */
/*       else converged = ( ((Rnorm/Rnorm0<tol) || (Rnorm<10.0*MEPS)) && (iNewton>4) ); */
/*       if (Rnorm == 0.0) converged = xfe_True; */
      
/*       /\*       xf_printf(" Rnorm0 = %.10E, Rnorm = %.15E, converged = %d\n", Rnorm0, Rnorm, converged); *\/ */
      
/*       if (converged) break; */
      
/*       if (iNewton != MaxNewton-1){ */
/* 	// calculate dxi */
/* 	if (dm1 == 1){ */
/* 	  if (R_xi[0] == 0.0){ */
/* 	    singular = xfe_True; */
/* 	    break; */
/* 	  } */
/* 	  dxi[0] = -R[0]/R_xi[0]; */
/* 	  dximag = fabs(dxi[0]); */
/* 	} */
/* 	else{ */
/* 	  ierr = xf_Error(xf_MatDetInv(R_xi, 2, &det, iR_xi)); */
/* 	  if (ierr != xf_OK) return ierr; */
/* 	  if (det == 0.0){ */
/* 	    singular = xfe_True; */
/* 	    break; */
/* 	  } */
/* 	  xf_MxV_Neg(iR_xi, R, 2, 2, dxi); */
/* 	  dximag = sqrt(dxi[0]*dxi[0] + dxi[1]*dxi[1]); */
/* 	} */
	
/* 	// increment xi, with damping if necessary */
/* 	omega = ((dximag > MaxStep) ? MaxStep/dximag : 1.0); */
/* 	if (iNewton < 3) omega *= 0.5; */
/* 	if (iNewton > 7) omega *= 0.5;//1.0/( (real) iNewton-6.0); */
/* 	for (d=0; d<dm1; d++) xi[d] += omega*dxi[d]; */
	
/*       } */
    
/*     }  // end Newton-Raphson on face */
    
/*     if (singular){ */
/*       //xf_printf("Warning, 1 singularity encountered.  Continuing.\n"); */
/*       // nudge x0in -> x0 to prevent singularity errors */
/*       for (d=0; d<dim; d++) xmin[d] = xmax[d] = X[d]; */
/*       for (i=1; i<nf; i++) */
/* 	for (d=0; d<dim; d++){ */
/* 	  xmin[d] = min(xmin[d], X[i*dim+d]); */
/* 	  xmax[d] = max(xmax[d], X[i*dim+d]); */
/* 	} */
/*       for (d=0, hmax=0.; d<dim; d++) hmax = max(xmax[d]-xmin[d], hmax); */
/*       for (d=0; d<dim; d++) x0[d] += 1.7*PI*MEPS*hmax*(d+1.0); */
/*       //xf_printf("  New x0 = %.10E %.10E\n", x0[0], x0[1]); */
/*     }  */
/*   }     */
	
/*   // warn if singularity could not be resolved */
/*   if (singular){ */
/*     xf_printf("Warning, singularity encountered in face distance calculation.  Continuing.\n"); */
/*   } */
/*   else{  // otherwise set distance */
    
/*     // if xi inside face, set mindist */
/*     ierr = xf_Error(xf_InsideShape(FShape, xi, tolinside, &inside)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     if (inside) mindist = min(mindist, dist); */
/*   } */

/*   /\*-------------*\/ */
/*   /\* Edge Newton *\/ */
/*   /\*-------------*\/ */

/*   // if dim == 3, consider edges as well, via Newton-Rapshon */
/*   if (dim == 3){ */

/*     //xf_printf("\nEdge iterations:\n"); */

/*     // get number of edges */
/*     ierr = xf_Error(xf_Shape2nFace(FShape, &nedge)); */
/*     if (ierr != xf_OK) return ierr; */

/*     // (over)allocate ev, used for edge node vector */
/*     ierr = xf_Error(xf_Alloc( (void **) &ev, nf, sizeof(int))); */
/*     if (ierr != xf_OK) return ierr; */

/*     // also (over)allocate XE, DE */
/*     ierr = xf_Error(xf_Alloc( (void **) &XE, nf*dim, sizeof(real))); */
/*     if (ierr != xf_OK) return ierr; */
/*     ierr = xf_Error(xf_Alloc( (void **) &DE, nf    , sizeof(real))); */
/*     if (ierr != xf_OK) return ierr; */

/*     // loop over edges */
/*     for (iedge=0; iedge<nedge; iedge++){ */

/*       //xf_printf("Edge %d\n", iedge); */
    
/*       // get edge shape */
/*       ierr = xf_Error(xf_FaceShape(FShape, iedge, &EShape)); */
/*       if (ierr != xf_OK) return ierr; */

/*       // get default Lagrange basis for edge */
/*       ierr = xf_Error(xf_Shape2UniformLagrange(EShape, &EBasis)); */
/*       if (ierr != xf_OK) return ierr; */
    
/*       // get nodes on edge */
/*       ierr = xf_Error(xf_NodesOnFace(FBasis, Order, iedge, &ne, ev)); */
/*       if (ierr != xf_OK) return ierr; */

/*       // pull off locations and distances at edge nodes */
/*       for (i=0; i<ne; i++){ */
/* 	for (d=0; d<dim; d++) XE[dim*i+d] = X[dim*ev[i]+d]; */
/* 	DE[i] = D[ev[i]]; */
/*       } */

/*       // set initial s (edge param) */
/*       xi[0] = 0.5; */
      
/*       // initialize variables before Newton-Raphson */
/*       converged = xfe_False; */
/*       dist      = mindist; */

/*       // Begin Newton-Raphson iteration */
/*       isafety  = 0; */
/*       singular = xfe_True; */
/*       while ((singular) && (isafety++ < 3)){ */
	
/* 	singular = xfe_False; */

/* 	for (iNewton=0; iNewton<MaxNewton; iNewton++){ */
	  
/* 	  // compute basis shape, grad, hess at xi */
/* 	  ierr = xf_Error(xf_EvalBasis(EBasis, Order, xfe_True, 1, xi, xfb_All, &PhiData)); */
/* 	  if (ierr != xf_OK) return ierr; */
	  
/* 	  // compute residual, R{i}, linearization, R{i}_xi{j}, and distance, dist */
/* 	  ierr = xf_Error(xf_Dist2FaceResidual(dim, x0, ne, XE, DE, PhiData, R, R_xi, &dist)); */
/* 	  if (ierr != xf_OK) return ierr; */
	  
/* 	  /\* 	  xf_printf("iNewton=%d, xi = %.3E, R = %.3E, dist = %.3E\n", iNewton, xi[0], R[0], dist); *\/ */
/* 	  /\* 	  xf_printf("  R_xi = %.3E\n", R_xi[0]); *\/ */
	  
/* 	  // convergence check */
/* 	  Rnorm = sqrt(R[0]*R[0]); */
/* 	  if (iNewton == 0) Rnorm0 = Rnorm; */
/* 	  else converged = ( ((Rnorm/Rnorm0<tol) || (Rnorm<10.0*MEPS)) && (iNewton>2) ); */
/* 	  if (Rnorm == 0.0) converged = xfe_True; */
	  
/* 	  /\* 	  xf_printf(" Rnorm0 = %.10E, Rnorm = %.15E, converged = %d\n", Rnorm0, Rnorm, converged); *\/ */
	  
/* 	  if (converged) break; */
	  
/* 	  if ((!converged) && (iNewton != MaxNewton-1)){ */
/* 	    // calculate dxi */
/* 	    if (R_xi[0] == 0.0) { */
/* 	      singular = xfe_True; */
/* 	      break; */
/* 	    } */
/* 	    dxi[0] = -R[0]/R_xi[0]; */
/* 	    dximag = fabs(dxi[0]); */
	    
/* 	    // increment xi, with damping if necessary */
/* 	    omega = ((dximag > MaxStep) ? MaxStep/dximag : 1.0); */
/* 	    xi[0] += omega*dxi[0];  */
/* 	  } */
/* 	}  // end Newton-Raphson */
	
/* 	if (singular){ */
/* 	  // nudge x0in -> x0 to prevent singularity errors */
/* 	  for (d=0; d<dim; d++) xmin[d] = xmax[d] = X[d]; */
/* 	  for (i=1; i<nf; i++) */
/* 	    for (d=0; d<dim; d++){ */
/* 	      xmin[d] = min(xmin[d], X[i*dim+d]); */
/* 	      xmax[d] = max(xmax[d], X[i*dim+d]); */
/* 	    } */
/* 	  for (d=0, hmax=0.; d<dim; d++) hmax = max(xmax[d]-xmin[d], hmax); */
/* 	  for (d=0; d<dim; d++) x0[d] = x0in[d] + 1.7*PI*MEPS*hmax; */
/* 	}  */
	
/*       } */

/*       // warn if singularity could not be resolved */
/*       if (singular){ */
/* 	xf_printf("Warning, singularity encountered in edge distance calculation.  Continuing.\n"); */
/*       } */
/*       else{  // otherwise set distance */
/* 	// if xi inside edge, set mindist */
/* 	ierr = xf_Error(xf_InsideShape(EShape, xi, tolinside, &inside)); */
/* 	if (ierr != xf_OK) return ierr; */
	
/* 	if (inside) mindist = min(mindist, dist); */
/*       } */
      
/*     } // iedge */

/*     xf_Release( (void *) ev); */
/*     xf_Release( (void *) XE); */
/*     xf_Release( (void *) DE); */

/*   } // end edge consideration */



/*   /\* Destroy Basis Data *\/ */
/*   ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */


/*   // return minimum distance */
/*   (*pdist) = mindist; */


/*   return xf_OK; */

/* } */


/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_Dist2Face */
/* static int */
/* xf_Dist2Face(int dim, const real *x0in, enum xfe_BasisType FBasis, */
/* 	     int Order, int nf, real *X, real *D, real *pdist) */
/* { */
/* /\* */
/* PURPOSE: */

/*   Calculates minimum: */

/*     (distance between x0 and a Lagrange-interpolated face) */
/*   + (distance function interpolated at face) */

/*   i.e. */

/*      (*pdist) = min_(xi{i}) [(r{k} * r{k})^.5 + d] */
  
/*   where  r{k} = x{k} - x0{k} */
/*          x{k} = X{n,k}*Phi{n}(xi{i}) */
/*          d    = D{n}*Phi{n}(xi{i}) */

/*   Note, the face geometry and the current distance function are */
/*   interpolated using the same Lagrange basis and order.  This function */
/*   uses subdivision of the face into subfaces, and hence is approximate */
/*   for curved faces.  It is also approximate if the interpolated */
/*   distance function varies across the face (occurs in parallel only). */
/*   On the plus side, this function is more robust than the Newton */
/*   version. */
  

/* INPUTS: */
 
/*   dim  : dimension of x0 */
/*   x0   : point from which to compute distance */
/*   FBasis : Lagrange basis of face */
/*   Order : Order with which face is interpolated */
/*   nf : number of Lagrange nodes interpolating face */
/*   X  : coords of face Lagrange points */
/*   D  : distance function values at Lagrange points */
  
/* OUTPUTS: */

/*   (*pdist)  : calculated minimum distance, as described above */
 
/* RETURNS: Error code */

/* *\/ */
/*   int ierr; */
/*   int i, k, dm1; */
/*   int refine, nnode; */
/*   int nsplit, isplit; */
/*   int *iv; */
/*   int *vsplit = NULL; */
/*   enum xfe_ShapeType FShape;  */
/*   real mindist, dist; */
/*   real x0[3]; */
/*   real xtri[9]; */
/*   real xseg[4]; */
/*   real *x = NULL; */
/*   real *d = NULL; */
/*   real *coord = NULL; */
/*   xf_BasisData *PhiData = NULL; */

/*   dm1 = dim-1; // face dimension */

/*   if ((dm1 != 1) && (dm1 != 2)) return xf_Error(xf_NOT_SUPPORTED); */

/*   if (nf <= 1) return xf_Error(xf_INPUT_ERROR); */

/*   // set x0 to input value (in case have to modify later) */
/*   for (i=0; i<dim; i++) x0[i] = x0in[i]; */

/*   // initialize mindist */
/*   mindist = 1e30; */

/*   // get face shape */
/*   ierr = xf_Error(xf_Basis2Shape(FBasis, &FShape)); */
/*   if (ierr != xf_OK) return ierr; */


/*   // split up shape into segments or triangles */
/*   refine = 2*Order+2; // refinement index */
/*   ierr = xf_Error(xf_GetRefineCoords(FShape, refine, &nnode, &coord, */
/* 				     &nsplit, &vsplit, NULL, NULL)); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   // evaluate basis functions at coord */
/*   PhiData = NULL; */
/*   ierr = xf_Error(xf_EvalBasis(FBasis, Order, xfe_True, nnode, */
/* 			       coord, xfb_Phi, &PhiData)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // allocate memory */
/*   ierr = xf_Error(xf_ReAlloc((void **) &x, nnode*dim, sizeof(real))); */
/*   if (ierr != xf_OK) return ierr; */
/*   ierr = xf_Error(xf_ReAlloc((void **) &d, nnode    , sizeof(real))); */
/*   if (ierr != xf_OK) return ierr; */

/*   // interpolate x(coord) and d(coord) */
/*   xf_MxM_Set(PhiData->Phi, X, nnode, nf, dim,  x); */
/*   xf_MxM_Set(PhiData->Phi, D, nnode, nf,   1,  d); */

/*   // loop over subshapes = segments/triangles */
/*   for (isplit=0; isplit<nsplit; isplit++){ */

/*     // calculate distance to subshape */
/*     if (dm1 == 2){ */
/*       // distance to triangle */
/*       iv = vsplit + 3*isplit; */
/*       for (k=0; k<3; k++)  */
/* 	for (i=0; i<3; i++) xtri[k*3+i] = x[dim*iv[k]+i]; */
/*       ierr = xf_Error(xf_Dist2Triangle(xtri, x0, &dist)); */
/*       if (ierr != xf_OK) return ierr; */
/*       for (k=0; k<3; k++) dist += d[iv[k]]/3.0; // only approximate */
/*     } */
/*     else{ */
/*       // distance to line */
/*       iv = vsplit + 2*isplit; */
/*       for (k=0; k<2; k++)  */
/* 	for (i=0; i<2; i++) xseg[k*2+i] = x[dim*iv[k]+i]; */
/*       ierr = xf_Error(xf_Dist2Segment(xseg, x0, &dist)); */
/*       if (ierr != xf_OK) return ierr; */
/*       for (k=0; k<2; k++) dist += d[iv[k]]/2.0; // only approximate */
/*     } */

/*     // set mindist appropriately */
/*     mindist = min(mindist, dist); */
/*   } // isplit */


/*   /\* Destroy Basis Data *\/ */
/*   ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // Release memory */
/*   xf_Release((void *) coord); */
/*   xf_Release((void *) vsplit); */
/*   xf_Release((void *) x); */
/*   xf_Release((void *) d); */

/*   // return minimum distance */
/*   (*pdist) = mindist; */


/*   return xf_OK; */

/* } */




/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_CalcDistFcnElems */
/* static int  */
/* xf_CalculateDistFcnElems(xf_All *All, int nface, xf_Face *FaceList, */
/* 			 xf_Vector *DistFcn, xf_Vector *ElemChanged) */
/* { */
/* /\* */
/* PURPOSE:  */

/*   Within each non-halo element, calculates distances to the set of */
/*   faces in FaceList.  The distance function on each element is then */
/*   defined as the minimum total distance, where the total distance */
/*   consists of the distance to a face plus the distance function at */
/*   that face.  The distance function at boundary faces is assumed to be */
/*   zero, while the distance function on interior halo-adjacent faces is */
/*   taken from the halo element.  Note, if not parallel, no interior */
/*   faces should be present in FaceList.  The distance function is */
/*   interpolated on each element by the order specified in the DistFcn */
/*   vector (Lagrange interpolation).  If the calculated minimum distance */
/*   is smaller than the one already present, ElemChanged is set to 1. */

/*   Some heuristics are employed in this calculation, and the distance */
/*   may not be the minimum distance if some very highly nonlinear, */
/*   high-order geometry elements are present.  Anisotropic/skewed */
/*   elements are fine. */

/* INPUTS:  */
 
/*   All: All structure */
/*   nface : number of faces to consider in distance calculation */
/*   FaceList : list of faces to consider */

/* OUTPUTS: */

/*   DistFcn : distance function on each element (valid on input, but updated) */
/*   ElemChanged: elements with new distance function values have this set to 1 */

/* RETURNS: Error code */

/* *\/ */

/*   int ierr, i, j, dim, d; */
/*   int ibfgrp, ibface; */
/*   int iegrp, ielem; */
/*   int egrp, elem, face; */
/*   int nfnode, nNode, nnode, inode; */
/*   int ipoint, npoint; */
/*   int pegrp, fi, nn; */
/*   int nf, pnf, pnn; */
/*   int Order; */
/*   int faceorient; */
/*   int fvec[xf_MAXQ1FACENODE]; */
/*   int *Node; */
/*   int *Node2nFace; */
/*   int *NodeList; */
/*   int *fv = NULL; */
/*   enum xfe_Bool BFlag, changed; */
/*   enum xfe_ShapeType Shape, FShape; */
/*   enum xfe_BasisType FBasis; */
/*   real mindist, dist, currdist; */
/*   real dtol = 0.99;              // tolerance for alerting about a minimum */
/*   real *NodeDist; */
/*   real *xpoint; */
/*   real *xref  = NULL, *xface = NULL, *xglob = NULL; */
/*   real *xv    = NULL, *dv    = NULL, *Xv    = NULL; */
/*   real *xelem = NULL; */
/*   xf_Face  **Node2Face; */
/*   xf_Face  Face; */
/*   xf_IFace IFace; */
/*   xf_ElemGroup *EG; */
/*   xf_BasisData *PhiData  = NULL; // for q1 node interpolation  */
/*   xf_BasisData *GPhiData = NULL; // for ref2glob transformation */
/*   xf_BasisData *QPhiData = NULL; // for interpolating geometry */
/*   xf_BasisData *DPhiData = NULL; // for interpolating DistFcn */
/*   xf_Mesh *Mesh; */

/*   Mesh  = All->Mesh; */
/*   dim   = Mesh->Dim; */
/*   nNode = Mesh->nNode; */

/*   // return immediately if nothing to do */
/*   if (nface <= 0) return xf_OK; */

/*   // allocate an integer counter at nodes, of adjacent faces */
/*   ierr = xf_Error(xf_Alloc( (void **) &Node2nFace, nNode, sizeof(int))); */
/*   if (ierr != xf_OK) return ierr; */
/*   for (i=0; i<nNode; i++) Node2nFace[i] = 0; */
  

/*   // flag all Q1 nodes involved (i.e. with an adjacent face in FaceList) */
/*   for (i=0; i<nface; i++){ */
/*     // elements adjacent to FaceList[i] */
/*     xf_FaceElements(Mesh, FaceList[i].Group, FaceList[i].Number,  */
/* 		    &egrp, &elem, &face, NULL, NULL, NULL); */
    
/*     // local Q1 nodes on face of elem */
/*     EG = Mesh->ElemGroup + egrp; */
/*     ierr = xf_Error(xf_Q1NodesOnFace(EG->QBasis, EG->QOrder, face, &nfnode, fvec)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     Node = Mesh->ElemGroup[egrp].Node[elem]; */
/*     for (j=0; j<nfnode; j++) Node2nFace[Node[fvec[j]]] += 1; */
/*   } // i */

  

/*   /\* Create Node2Face = list of faces adjacent to the involved Q1 nodes *\/ */

/*   // allocate Node2Face, variable length */
/*   ierr = xf_Error(xf_VAlloc2( (void ***) &Node2Face, nNode, Node2nFace, sizeof(xf_Face))); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   // fill in Node2Face */
/*   for (i=0; i<nNode; i++) Node2nFace[i] = 0; */
/*   for (i=0; i<nface; i++){ */
/*     // elements adjacent to FaceList[i] */
/*     xf_FaceElements(Mesh, FaceList[i].Group, FaceList[i].Number,  */
/* 		    &egrp, &elem, &face, NULL, NULL, NULL); */
    
/*     // local Q1 nodes on face of elem */
/*     EG = Mesh->ElemGroup + egrp; */
/*     ierr = xf_Error(xf_Q1NodesOnFace(EG->QBasis, EG->QOrder, face, &nfnode, fvec)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     Node = Mesh->ElemGroup[egrp].Node[elem]; */
/*     for (j=0; j<nfnode; j++){ */
/*       inode = Node[fvec[j]]; */
/*       Node2Face[inode][Node2nFace[inode]++] = FaceList[i]; */
/*     } */
/*   } // i */

  
/*   // allocate NodeList = list of Q1 nodes involved */
/*   for (i=0, nnode=0; i<nNode; i++) nnode += (Node2nFace[i] > 0); */
/*   if (nnode <= 0) return xf_Error(xf_MESH_ERROR); // faces must have adjacent Q1 nodes */
/*   ierr = xf_Error(xf_Alloc( (void **) &NodeList, nnode, sizeof(int))); */
/*   if (ierr != xf_OK) return ierr; */

/*   // fill in NodeList */
/*   for (i=0, nnode=0; i<nNode; i++)  */
/*     if (Node2nFace[i] > 0) NodeList[nnode++] = i; */



/*   /\* Create NodeDist = current distance function at involved Q1 nodes *\/ */

/*   // allocate NodeDist */
/*   ierr = xf_Error(xf_Alloc( (void **) &NodeDist, nNode, sizeof(real))); */
/*   if (ierr != xf_OK) return ierr; */
/*   for (i=0; i<nNode; i++) NodeDist[i] = 1e30; */

/*   // fill in NodeDist */
/*   pegrp   = -1; */
/*   for (i=0; i<nface; i++){ */
/*     // are we on a boundary face? */
/*     BFlag = (FaceList[i].Group >= 0); */

/*     // halo element adjacent to FaceList[i] (or self if BFlag == True) */
/*     ierr = xf_Error(xf_FaceElementHalo(Mesh, FaceList[i].Group, FaceList[i].Number,  */
/* 				       !BFlag, &egrp, &elem, &face)); */
/*     if (ierr != xf_OK) return ierr; */

/*     if (egrp < 0){ // indicates a boundary face */
/*       ierr = xf_Error(xf_FaceElementHalo(Mesh,  FaceList[i].Group, FaceList[i].Number,  */
/* 					 xfe_False, &egrp, &elem, &face)); */
/*       if (ierr != xf_OK) return ierr; */
/*     } */
      

/*     // local Q1 nodes on face of elem */
/*     EG = Mesh->ElemGroup + egrp; */
/*     ierr = xf_Error(xf_Q1NodesOnFace(EG->QBasis, EG->QOrder, face, &nfnode, fvec)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     Node = Mesh->ElemGroup[egrp].Node[elem]; */

/*     if (BFlag){ // distance is zero for boundary faces */
/*       for (j=0; j<nfnode; j++){ */
/* 	inode = Node[fvec[j]]; */
/* 	NodeDist[inode] =0.0; */
/*       } */
/*     } */
/*     else{ // interpolate distance from halo */

/*       if (egrp != pegrp){ */
/* 	// reference-space coords of all geometry nodes */
/* 	ierr = xf_Error(xf_LagrangeNodes(EG->QBasis, EG->QOrder, &nn, NULL, &xref)); */
/* 	if (ierr != xf_OK) return ierr; */
/* 	// compute DistFcn basis functions at geometry nodes */
/* 	ierr = xf_Error(xf_EvalBasis(DistFcn->Basis[egrp], DistFcn->Order[egrp],  */
/* 				     xfe_True, nn, xref, xfb_Phi, &PhiData)); */
/* 	if (ierr != xf_OK) return ierr; */
/* 	// realloc dv = distance vector */
/* 	ierr = xf_Error(xf_ReAlloc( (void **) &dv, nn, sizeof(real))); */
/* 	if (ierr != xf_OK) return ierr; */
/* 	pegrp = egrp; */
/*       } */
      
/*       if (nn != PhiData->nq) return xf_Error(xf_CODE_LOGIC_ERROR); */

/*       // interpolate DistFcn at all geometry nodes */
/*       xf_MxM_Set(PhiData->Phi, DistFcn->GenArray[egrp].rValue[elem], nn, PhiData->nn, 1, dv); */
      
/*       // set distance at nodes to minimum of current and interpolated values */
/*       for (j=0; j<nfnode; j++){ */
/* 	inode = Node[fvec[j]]; // this is a Q1 node */
/* 	NodeDist[inode] = min(NodeDist[inode], dv[fvec[j]]); */
/*       } */
/*     } */
     
/*   } // i */

/*   // destroy PhiData */
/*   ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */



/*   /\*-------------------------*\/ */
/*   /\* Main loop over elements *\/ */
/*   /\*-------------------------*\/ */

/*   pnf = -1; // previous number of face Lagrange points, for allocation */
/*   pnn = -1; // previous number of elem Lagrange points, for allocation */

/*   // loop over all element groups */
/*   for (iegrp=0; iegrp<Mesh->nElemGroup; iegrp++){ */

/*     // reference-space coords of Lagrange nodes */
/*     ierr = xf_Error(xf_LagrangeNodes(DistFcn->Basis[iegrp], DistFcn->Order[iegrp], &npoint, NULL, &xref)); */
/*     if (ierr != xf_OK) return ierr; */

/*     // reallocate xglob to store global coords */
/*     ierr = xf_Error(xf_ReAlloc((void **) &xglob, dim*npoint, sizeof(real))); */
/*     if (ierr != xf_OK) return ierr; */

/*     // loop over elements */
/*     for (ielem=0; ielem<Mesh->ElemGroup[iegrp].nElem; ielem++){ */
/*       changed = xfe_False; */

/*       // ref2glob -> global coordinates of the points */
/*       ierr = xf_Error(xf_Ref2GlobElem(Mesh, iegrp, ielem, &GPhiData, xfe_True, */
/* 				      npoint, xref, xglob)); */
/*       if (ierr != xf_OK) return ierr; */

/*       // loop over points */
/*       for (ipoint=0; ipoint<npoint; ipoint++){ */
	
/* 	xpoint = xglob + dim*ipoint; */

/* 	// calculate distance to closest Q1 node */
/* 	inode = NodeList[0]; */
/* 	mindist = NodeDist[inode] + xf_Distance(xpoint, Mesh->Coord[inode], dim); */
/* 	for (i=1; i<nnode; i++){ */
/* 	  j = NodeList[i]; */
/* 	  dist = NodeDist[j] + xf_Distance(xpoint, Mesh->Coord[j], dim); */
/* 	  if (dist < mindist){ */
/* 	    mindist = dist; */
/* 	    inode   = j; */
/* 	  } */
/* 	} */
	
/* 	/\* 	// TEMPORARY *\/ */
/* 	/\* 	if ( (fabs(xpoint[0]-5.2716449019E-02) < .0005) && (fabs(xpoint[1]-3.8951786423E-02) < 0.0005)){ *\/ */
/* 	/\* 	  xf_printf("xpoint = %.10E %.10E\n", xpoint[0], xpoint[1]); *\/ */
/* 	/\* 	  xf_printf("mindist = %.10E\n", mindist); *\/ */
/* 	/\* 	  xf_printf("inode = %d\n", inode); *\/ */
/* 	/\* 	  xf_printf("xnode = %.10E %.10E\n", Mesh->Coord[inode][0], Mesh->Coord[inode][1]); *\/ */
/* 	/\* 	  DEBUGFLAG = xfe_True; *\/ */
/* 	/\* 	} *\/ */
/* 	/\* 	else *\/ */
/* 	/\* 	  DEBUGFLAG = xfe_False; *\/ */
	
/* 	// loop over faces adjacent to inode */
/* 	for (fi=0; fi<Node2nFace[inode]; fi++){ */
/* 	  Face = Node2Face[inode][fi]; */

/* 	  BFlag = (Face.Group >= 0); // on boundary? */

/* 	  // halo element adjacent to Face (or self if on boundary) */
/* 	  ierr = xf_Error(xf_FaceElementHalo(Mesh, Face.Group, Face.Number, */
/* 					     !BFlag, &egrp, &elem, &face)); */
/* 	  if (ierr != xf_OK) return ierr; */

/* 	  EG = Mesh->ElemGroup + egrp; */

	 
/* 	  if (BFlag) */
/* 	    // set Order to geometry order */
/* 	    Order = Mesh->ElemGroup[egrp].QOrder; */
/* 	  else */
/* 	    // set Order to max of geometry order and DistFcn interpolation order */
/* 	    Order = max(Mesh->ElemGroup[egrp].QOrder, DistFcn->Order[egrp]); */
	  

/* 	  // get element shape */
/* 	  ierr = xf_Error(xf_Basis2Shape(EG->QBasis, &Shape)); */
/* 	  if (ierr != xf_OK) return ierr; */
	  
/* 	  // get face shape */
/* 	  ierr = xf_Error(xf_FaceShape(Shape, face, &FShape)); */
/* 	  if (ierr != xf_OK) return ierr; */
	  
/* 	  // get default Lagrange basis for face */
/* 	  ierr = xf_Error(xf_Shape2UniformLagrange(FShape, &FBasis)); */
/* 	  if (ierr != xf_OK) return ierr; */
	  
/* 	  // number of Lagrange nodes on face, nf */
/* 	  ierr = xf_Error(xf_Order2nNode(FBasis, Order, &nf)); */
/* 	  if (ierr != xf_OK) return ierr; */

/* 	  if (BFlag){ */

/* 	    // reallocate memory depending on Order */
/* 	    if (pnf < nf){ */
/* 	      ierr = xf_Error(xf_ReAlloc( (void **) &fv   , nf    , sizeof(int))); */
/* 	      if (ierr != xf_OK) return ierr; */

/* 	      ierr = xf_Error(xf_ReAlloc( (void **) &xv   , nf*dim, sizeof(real))); */
/* 	      if (ierr != xf_OK) return ierr; */

/* 	      ierr = xf_Error(xf_ReAlloc( (void **) &dv   , nf    , sizeof(real))); */
/* 	      if (ierr != xf_OK) return ierr; */

/* 	      pnf = nf; */
/* 	    } */

/* 	    // take face nodes from self elem */
/* 	    ierr = xf_Error(xf_NodesOnFace(EG->QBasis, EG->QOrder, face, &nfnode, fv)); */
/* 	    if (ierr != xf_OK) return ierr; */
/* 	    if (nfnode != nf) return xf_Error(xf_CODE_LOGIC_ERROR); */

/* 	    // set xv using fv */
/* 	    Node = Mesh->ElemGroup[egrp].Node[elem]; */
/* 	    for (i=0; i<nfnode; i++) */
/* 	      for (d=0; d<dim; d++) */
/* 		xv[dim*i + d] = Mesh->Coord[Node[fv[i]]][d]; */

/* 	    // dv = 0 */
/* 	    for (i=0; i<nfnode; i++) dv[i] = 0.; */

/* 	  } */
/* 	  else{ */
/* 	    // interior face */
/* 	    IFace = Mesh->IFace[Face.Number]; */

/* 	    // face orientation w.r.t halo */
/* 	    faceorient = ((egrp == IFace.ElemGroupL) ? IFace.OrientL : IFace.OrientR); */

/* 	    // reallocate memory depending on Order */
/* 	    if (pnf < nf){ */
/* 	      ierr = xf_Error(xf_ReAlloc( (void **) &dv   , nf    , sizeof(real))); */
/* 	      if (ierr != xf_OK) return ierr; */

/* 	      ierr = xf_Error(xf_ReAlloc( (void **) &xv   , nf*dim, sizeof(real))); */
/* 	      if (ierr != xf_OK) return ierr; */

/* 	      ierr = xf_Error(xf_ReAlloc( (void **) &xface, nf*dim, sizeof(real))); */
/* 	      if (ierr != xf_OK) return ierr; */

/* 	      ierr = xf_Error(xf_ReAlloc( (void **) &xelem, nf*dim, sizeof(real))); */
/* 	      if (ierr != xf_OK) return ierr; */

/* 	      pnf = nf; */
/* 	    } */

/* 	    // pull off face Lagrange points for Order, xface */
/* 	    ierr = xf_Error(xf_LagrangeNodesEqual(FShape, Order, xface)); */
/* 	    if (ierr != xf_OK) return ierr; */

/* 	    // convert xface to xelem */
/* 	    ierr = xf_Error(xf_RefFace2Interpol(Mesh, egrp, elem, face, faceorient, */
/* 						nf, xface, xelem)); */
/* 	    if (ierr != xf_OK) return ierr; */

/* 	    // Interpolate geometry at xelem, xv */
/* 	    ierr = xf_Error(xf_EvalBasis(EG->QBasis, EG->QOrder, */
/* 					 xfe_True, nf, xelem, xfb_Phi, &QPhiData)); */
/* 	    if (ierr != xf_OK) return ierr; */

/* 	    if (QPhiData->nn > pnn){ */
/* 	      ierr = xf_Error(xf_ReAlloc( (void **) &Xv, QPhiData->nn*dim, sizeof(real))); */
/* 	      if (ierr != xf_OK) return ierr; */
/* 	      pnn = QPhiData->nn; */
/* 	    } */
/* 	    Node = Mesh->ElemGroup[egrp].Node[elem]; */
/* 	    for (i=0; i<QPhiData->nn; i++) */
/* 	      for (d=0; d<dim; d++) */
/* 		Xv[dim*i + d] = Mesh->Coord[Node[i]][d]; */

/* 	    xf_MxM_Set(QPhiData->Phi, Xv, nf, QPhiData->nn, dim, xv); */

/* 	    // Interpolate DistFcn at xelem, dv */
/* 	    ierr = xf_Error(xf_EvalBasis(DistFcn->Basis[egrp], DistFcn->Order[egrp], */
/* 					 xfe_True, nf, xelem, xfb_Phi, &DPhiData)); */
/* 	    if (ierr != xf_OK) return ierr; */
/* 	    xf_MxM_Set(DPhiData->Phi, DistFcn->GenArray[egrp].rValue[elem], nf, */
/* 		       DPhiData->nn, 1, dv); */
/* 	  } */

	  
/* 	  // calculate min distance to face, including dv */
/* 	  ierr = xf_Error(xf_Dist2Face(dim, xpoint, FBasis, Order, nf, xv, dv, &dist)); */
/* 	  if (ierr != xf_OK) return ierr; */

/* 	  // set minimum distance */
/* 	  mindist = min(mindist, dist); */
	  
/* 	} // fi */
	
/* 	// check if mindist is smaller than current distance */
/* 	if (mindist < (currdist = DistFcn->GenArray[iegrp].rValue[ielem][ipoint])){ */
/* 	  changed = (mindist < dtol*currdist); // tolerance to avoid excessive iterations */
/* 	  DistFcn->GenArray[iegrp].rValue[ielem][ipoint] = mindist; */
/* 	} */

/*       } // ipoint */


/*       // set ElemChanged if one of the distances changed */
/*       ElemChanged->GenArray[iegrp].iValue[ielem][0] = changed; */
/*     } // ielem */
/*   } // iegrp */


/*   // destroy GPhiData */
/*   ierr = xf_Error(xf_DestroyBasisData(GPhiData, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */
/*   // destroy QPhiData */
/*   ierr = xf_Error(xf_DestroyBasisData(QPhiData, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */
/*   // destroy DPhiData */
/*   ierr = xf_Error(xf_DestroyBasisData(DPhiData, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */



/*   // release memory */
/*   xf_Release ( (void  *) Node2nFace); */
/*   xf_Release2( (void **) Node2Face); */
/*   xf_Release ( (void  *) NodeList); */
/*   xf_Release ( (void  *) NodeDist); */
/*   xf_Release ( (void  *) xref); */
/*   xf_Release ( (void  *) xglob); */
/*   xf_Release ( (void  *) xface); */
/*   xf_Release ( (void  *) xelem); */
/*   xf_Release ( (void  *) fv); */
/*   xf_Release ( (void  *) dv); */
/*   xf_Release ( (void  *) xv); */
/*   xf_Release ( (void  *) Xv); */

/*   return xf_OK; */
/* } */



/* /\******************************************************************\/ */
/* //   FUNCTION Definition: xf_CalculateDistFcn */
/* int  */
/* xf_CalculateDistFcn(xf_All *All) */
/* { */
/*   int ierr; */
/*   int myRank, nProc; */
/*   int negrp, egrp; */
/*   int DistFcnOrder; */
/*   int ibface, negrphalo; */
/*   int nface, nfacetot; */
/*   int ibfgrp, nbfgrp; */
/*   int isafety, nsafety; */
/*   int nIFaceRegular; */
/*   int iiface; */
/*   int egrpH, elemH, faceH; */
/*   int nWall, i; */
/*   int *OrderVec = NULL; */
/*   char DistFcnWallBoundaries[xf_MAXSTRLEN]; */
/*   char **WallNames = NULL; */
/*   enum xfe_BasisType *BasisVec = NULL; */
/*   enum xfe_Bool ParallelFlag; */
/*   enum xfe_Bool IsWall; */
/*   enum xfe_Bool Found, done; */
/*   enum xfe_Bool EqnSetValid; */
/*   enum xfe_Bool WriteDistFcn; */
/*   enum xfe_Bool OverrideFlag; */
/*   xf_Face   *FaceList; */
/*   xf_Data   *D; */
/*   xf_Vector *DistFcn; */
/*   xf_Vector *ElemChanged; */
/*   xf_Mesh   *Mesh; */
/*   xf_EqnSet *EqnSet; */

/*   Mesh   = All->Mesh; */
/*   EqnSet = All->EqnSet; */

/*   // obtain number of processors */
/*   ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc)); */
/*   if (ierr != xf_OK) return ierr; */

/*   ParallelFlag = (nProc > 1); */

/*   // get order of the distance function */
/*   ierr = xf_Error(xf_GetKeyValueInt(All->Param->KeyValue, "DistFcnOrder", &DistFcnOrder)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // Check if DistFcn should be made writable */
/*   ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "WriteDistFcn", &WriteDistFcn)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // allocate Basis and Order vectors for vector search */
/*   negrp = Mesh->nElemGroup; */
/*   ierr = xf_Error(xf_Alloc( (void **) &BasisVec, negrp, sizeof(enum xfe_BasisType))); */
/*   if (ierr != xf_OK) return ierr; */
/*   ierr = xf_Error(xf_Alloc( (void **) &OrderVec, negrp, sizeof(int))); */
/*   if (ierr != xf_OK) return ierr; */

/*   // request Lagrange basis of order DistFcnOrder */
/*   for (egrp=0; egrp<negrp; egrp++){ */
/*     ierr = xf_Error(xf_Basis2UniformLagrange(Mesh->ElemGroup[egrp].QBasis, BasisVec+egrp)); */
/*     if (ierr != xf_OK) return ierr; */
/*     OrderVec[egrp] = DistFcnOrder; */
/*   } */

  
/*   // look for DistFcn vector (make parallel so can exchange halo) */
/*   ierr = xf_Error(xf_FindVector(All, "WallDistance", xfe_LinkageGlobElem, 1,  */
/* 				NULL, 0, 0, BasisVec, OrderVec, NULL, xfe_SizeReal,  */
/* 				xfe_True, xfe_True, &D, &DistFcn, &Found)); */
/*   if (ierr != xf_OK) return ierr; */
/*   D->ReadWrite = WriteDistFcn; */

/*   xf_Release( (void  *) BasisVec); */
/*   xf_Release( (void  *) OrderVec); */

/*   // consider returning immediately if DistFcn was found */

/*   // set distances to a very large number (will take minimum) */
/*   ierr = xf_Error(xf_SetConstVector(DistFcn, 0, 1e30)); */
/*   if (ierr != xf_OK) return ierr; */
  
/*   // look for ElemChanged vector (make parallel so can exchange halo) */
/*   ierr = xf_Error(xf_FindVector(All, "ElemChanged", xfe_LinkageGlobElem, 1,  */
/* 				NULL, 0, 0, NULL, NULL, NULL, xfe_SizeInt,  */
/* 				xfe_True, xfe_False, NULL, &ElemChanged, &Found)); */
/*   if (ierr != xf_OK) return ierr; */
/*   if (Found) return xf_Error(xf_CODE_LOGIC_ERROR); // should not exist */


/*   // number of bface groups */
/*   nbfgrp = Mesh->nBFaceGroup; */

/*   // total number of groups, including halo */
/*   negrphalo = (ParallelFlag ? 2*negrp : negrp); */

/*   // allocate a FaceList vector for storing faces to which we need distances */
/*   nface=0; */
/*   for (egrp=negrp; egrp<negrphalo; egrp++) nface += Mesh->ElemGroup[egrp].nElem; */
/*   for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++) nface += Mesh->BFaceGroup[ibfgrp].nBFace; */

/*   ierr = xf_Error(xf_Alloc((void **) &FaceList, nface, sizeof(xf_Face))); */
/*   if (ierr != xf_OK) return ierr; */
  

/*   // get list of any wall boundaries */
/*   ierr = xf_Error(xf_GetKeyValue(All->Param->KeyValue, "DistFcnWallBoundaries", DistFcnWallBoundaries)); */
/*   if (ierr != xf_OK) return ierr; */

/*   // Check if "DistFcnWallBoundaries" is not null -- use this first */
/*   OverrideFlag = xfe_False; */
/*   if (xf_NotNull(DistFcnWallBoundaries)){ */
/*     ierr = xf_Error(xf_ScanXStringAlloc(DistFcnWallBoundaries, xf_MAXSTRLEN, &nWall, &WallNames)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     if (nWall > 0) OverrideFlag = xfe_True; */
/*   } */


/*   // build a list of wall faces */
/*   nface = 0; */
/*   // if EqnSet is not initialized, set EqnSetValid=False */
/*   EqnSetValid = (EqnSet->BCs != NULL); */
/*   if (EqnSetValid){ */
/*     ierr = xf_Error(xf_SortEqnSetBCs(Mesh, EqnSet->BCs+0)); */
/*     if (ierr != xf_OK) return ierr; */
/*   } */


/*   for (ibfgrp=0; ibfgrp<nbfgrp; ibfgrp++){ */

/*     if (OverrideFlag){          // first check if overriding eqnset specification */
/*       IsWall = xfe_False; */
/*       for (i=0; i<nWall; i++){ */
/* 	if (strcmp(WallNames[i], Mesh->BFaceGroup[ibfgrp].Title) == 0){ */
/* 	  IsWall = xfe_True; */
/* 	  break; */
/* 	} */
/*       } // i */
/*     } */
/*     else if (EqnSetValid){  // use eqnset if valid */
/*       ierr = xf_Error(xf_BCIsWall(EqnSet->BCs[0].BC+ibfgrp, &IsWall)); */
/*       if (ierr != xf_OK) return ierr; */
/*     } */
/*     else IsWall = xfe_True; */

/*     if (IsWall){ // add to FaceList */
/*       for (ibface=0; ibface<Mesh->BFaceGroup[ibfgrp].nBFace; ibface++){ */
/* 	FaceList[nface+ibface].Group  = ibfgrp; */
/* 	FaceList[nface+ibface].Number = ibface; */
/*       } */
/*       nface += Mesh->BFaceGroup[ibfgrp].nBFace; */
/*     } */
/*   } // ibfgrp */

/*   xf_Release2( (void **) WallNames); */

/*   // error if no faces flagged (i.e. no walls) */
/*   nfacetot = nface; */
/*   ierr = xf_Error(xf_MPI_Allreduce(&nfacetot, 1, xfe_SizeInt, xfe_MPI_SUM)); */
/*   if (ierr != xf_OK) return ierr; */

/*   if (nfacetot <= 0){ */
/*     xf_printf("Error, no wall boundary faces identified in DistFcn calculation.\n"); */
/*     return xf_Error(xf_NOT_FOUND); */
/*   } */

  

/*   // begin main distance calculation loop (intended for parallel) */
/*   done    = xfe_False; */
/*   isafety = 0; */
/*   nsafety = 2*nProc; */
/*   while ((!done) && (isafety++ < nsafety)){ */

/*     xf_printf("iteration = %d\n", isafety); */

/*     // done if no more faces to which we need distances */
/*     nfacetot = nface; */
/*     ierr = xf_Error(xf_MPI_Allreduce(&nfacetot, 1, xfe_SizeInt, xfe_MPI_SUM)); */
/*     if (ierr != xf_OK) return ierr; */

/*     if (nfacetot == 0){ */
/*       done = xfe_True; */
/*       break; */
/*     } */

/*     // recalculate distances for each elem to faces in FaceList */
/*     ierr = xf_Error(xf_CalculateDistFcnElems(All, nface, FaceList, DistFcn, ElemChanged)); */
/*     if (ierr != xf_OK) return ierr; */
    
/*     // communicate ElemChanged */
/*     ierr = xf_Error(xf_HaloExchangeVectorBegin(ElemChanged)); */
/*     if (ierr != xf_OK) return ierr; */
/*     ierr = xf_Error(xf_HaloExchangeVectorEnd(ElemChanged)); */
/*     if (ierr != xf_OK) return ierr; */

/*     // communicate DistFcn */
/*     ierr = xf_Error(xf_HaloExchangeVectorBegin(DistFcn)); */
/*     if (ierr != xf_OK) return ierr; */
/*     ierr = xf_Error(xf_HaloExchangeVectorEnd(DistFcn)); */
/*     if (ierr != xf_OK) return ierr; */


/*     // create FaceList based on halo ElemChanged */
/*     nface = 0; */
/*     nIFaceRegular = ((ParallelFlag) ? Mesh->ParallelInfo->nIFaceRegular : Mesh->nIFace); */
/*     for (iiface=nIFaceRegular; iiface<Mesh->nIFace; iiface++){ // loop over halo ifaces */

/*       // pull off halo element adjacent to iiface */
/*       ierr = xf_Error(xf_FaceElementHalo(Mesh, -1, iiface, xfe_True, &egrpH, &elemH, &faceH)); */
/*       if (ierr != xf_OK) return ierr; */

/*       // if halo distance changed, mark face */
/*       if (ElemChanged->GenArray[egrpH].iValue[elemH][0] == 1){ */
/* 	FaceList[nface  ].Group  = -1; */
/* 	FaceList[nface++].Number = iiface; */
/*       } */

/*     } // iiface */

/*     // clear ElemChanged (clears halos too) */
/*     ierr = xf_Error(xf_SetZeroVector(ElemChanged)); */
/*     if (ierr != xf_OK) return ierr; */


/*   } // end while !done */


/*   // warning if did not finish iterations */
/*   if (!done) */
/*     xf_printf("Warning, maximum number of iterations taken in DistFcn calculation.\n"); */


/*   // destroy temporary vectors */
/*   ierr = xf_Error(xf_DestroyVector(ElemChanged, xfe_True)); */
/*   if (ierr != xf_OK) return ierr; */


/*   // release memory */
/*   xf_Release( (void *) FaceList); */


/*   return xf_OK; */
/* } */




#if( UNIT_TEST==1 )
#include "xf_MeshDistance.test.in"
#endif

