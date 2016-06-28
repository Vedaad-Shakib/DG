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


#include "xf_Unit.h"
#include "xf_All.h"
#include "xf_Mesh.h"


TEST_xf_InsideShape()
{
  int i, ierr;
  enum xfe_Bool inside;
  real xSeg[4][1]  = {{0.5},{-1e-6},{1.0},{2.0}};
  int  iSeg[4]     = {1,1,1,0};
  real xQuad[4][2] = {{0.5,0.5},{-1e-6,1e-6},{1.0,0.0},{2.0,0.5}};
  int  iQuad[4]    = {1,1,1,0};
  real xTri[4][2]  = {{0.2,0.2},{-1e-6,1e-6},{1.0,0.0},{0.5,0.6}};
  int  iTri[4]     = {1,1,1,0};
  real xTet[5][3]  = {{0.2,0.2,0.2},{-1e-6,1e-6,1e-6},{1.0,0.0,0.0},{0.3,0.4,0.4},
		      {0.6,0.1,0.6}};
  int  iTet[5]     = {1,1,1,0,0};
  real xHex[4][3]  = {{0.6,0.6,0.6},{-1e-6,1e-6,1e-6},{1.0,1.0,1.0},{0.4,0.4,1.5}};
  int  iHex[4]     = {1,1,1,0};
  real tol = 1e-5;
  
  // Seg
  for (i=0; i<4; i++){
    ierr = xf_Error(xf_InsideShape(xfe_Segment, xSeg[i], tol, &inside));
    if (ierr != xf_OK) return ierr;
    xf_AssertEqual(inside, iSeg[i]);
  }
 
  // Quad
  for (i=0; i<4; i++){
    ierr = xf_Error(xf_InsideShape(xfe_Quadrilateral, xQuad[i], tol, &inside));
    if (ierr != xf_OK) return ierr;
    xf_AssertEqual(inside, iQuad[i]);
  }
  
  // Tri
  for (i=0; i<4; i++){
    ierr = xf_Error(xf_InsideShape(xfe_Triangle, xTri[i], tol, &inside));
    if (ierr != xf_OK) return ierr;
    xf_AssertEqual(inside, iTri[i]);
  }

  // Tet
  for (i=0; i<5; i++){
    ierr = xf_Error(xf_InsideShape(xfe_Tetrahedron, xTet[i], tol, &inside));
    if (ierr != xf_OK) return ierr;
    xf_AssertEqual(inside, iTet[i]);
  }
 
  // Hex
  for (i=0; i<4; i++){
    ierr = xf_Error(xf_InsideShape(xfe_Hexahedron, xHex[i], tol, &inside));
    if (ierr != xf_OK) return ierr;
    xf_AssertEqual(inside, iHex[i]);
  }

  return xf_OK;
}


TEST_xf_NodesOnFace_Tri()
{
  int ierr, nfnode, fvec[4];
  int fvec0[] = {3,6,8,9};
  int fvec1[] = {9,7,4,0};
  int fvec2[] = {0,1,2,3};

  ierr = xf_Error(xf_NodesOnFace(xfe_TriLagrange, 3, 0, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nfnode, 4);
  xf_AssertIntVectorEqual(fvec, fvec0, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_TriLagrange, 3, 1, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec1, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_TriLagrange, 3, 2, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec2, nfnode);

  return xf_OK;  
}

TEST_xf_NodesOnFace_Tet()
{
  int ierr, nfnode, fvec[10];
  int fvec0[] = {3,6,8,9,12,14,15,17,18,19};
  int fvec1[] = {0,10,16,19,4,13,18,7,15,9};
  int fvec2[] = {0,1,2,3,10,11,12,16,17,19};
  int fvec3[] = {0,4,7,9,1,5,8,2,6,3};

  ierr = xf_Error(xf_NodesOnFace(xfe_TetLagrange, 3, 0, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nfnode, 10);
  xf_AssertIntVectorEqual(fvec, fvec0, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_TetLagrange, 3, 1, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec1, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_TetLagrange, 3, 2, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec2, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_TetLagrange, 3, 3, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec3, nfnode);

  return xf_OK;  
}


TEST_xf_NodesOnFace_Hex()
{
  int ierr, nfnode, fvec[9];
  int fvec0[] = {0,3,6,1,4,7,2,5,8};
  int fvec1[] = {0,1,2,9,10,11,18,19,20};
  int fvec2[] = {2,5,8,11,14,17,20,23,26};
  int fvec3[] = {8,7,6,17,16,15,26,25,24};
  int fvec4[] = {6,3,0,15,12,9,24,21,18};  
  int fvec5[] = {18,19,20,21,22,23,24,25,26};


  ierr = xf_Error(xf_NodesOnFace(xfe_HexLagrange, 2, 0, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(nfnode, 9);
  xf_AssertIntVectorEqual(fvec, fvec0, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_HexLagrange, 2, 1, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec1, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_HexLagrange, 2, 2, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec2, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_HexLagrange, 2, 3, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec3, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_HexLagrange, 2, 4, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec4, nfnode);

  ierr = xf_Error(xf_NodesOnFace(xfe_HexLagrange, 2, 5, &nfnode, fvec));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertIntVectorEqual(fvec, fvec5, nfnode);

  return xf_OK;  
}


TEST_xf_RotateNodes()
{
  int ierr, i, nfnode, fvec[9];
  int nvec1[4] = {1,3,7,5};
  int nvecf[4], nvec2[4];
  int nvec2_true[4] = {1,5,7,3};

  ierr = xf_Error(xf_RotateNodes(xfe_Quadrilateral, 1, 4, nvec1, nvecf));
  xf_AssertEqual(ierr, xf_OK);
  ierr = xf_Error(xf_InvRotateNodes(xfe_Quadrilateral, 7, 4, nvecf, nvec2));
  xf_AssertEqual(ierr, xf_OK);

  xf_AssertIntVectorEqual(nvec2, nvec2_true, 4);

  return xf_OK;  
}



TEST_xf_BasisLagrange1D()
{
  int ierr;
  real xnode[3] = {-1, 0, 1};
  real x = 0;
  real phi[3], gphi[3], hphi[3];
  real  phi_true[3] = {0, 1, 0};
  real gphi_true[3] = {-0.5, 0, 0.5};
  real hphi_true[3] = {1, -2, 1};

  ierr = xf_Error(xf_BasisLagrange1D(x, xnode, 3, phi, gphi, hphi));
  if (ierr != xf_OK) return ierr;

  xf_AssertRealVectorWithin( phi,  phi_true, 3, UTOL2);
  xf_AssertRealVectorWithin(gphi, gphi_true, 3, UTOL2);
  xf_AssertRealVectorWithin(hphi, hphi_true, 3, UTOL2);

  return xf_OK;  
}


TEST_xf_Shape_Grad_TensorLagrange()
{
  int ierr, i[3], itype, o, d, dim, iq, nq, k, nn;
  int d0, d1, irun;
  enum xfe_BasisType Basis;
  enum xfe_BasisType bases[6] = {xfe_SegLagrange, xfe_QuadLagrange, xfe_HexLagrange, 
				 xfe_QuadLagrangeGauss, xfe_HexLagrangeGauss};
  int ntype[6] = {0,0,0,1,1};
  real *xq, *wq, *x, s=0, t, gt[3], ht[9];
  real phi[216], gphi[648], hphi[1944], F[216], xn[648];

  xq = NULL;
  wq = NULL;

  for (irun=0; irun<2; irun++){

    Basis = bases[irun];
    itype = ntype[irun];

    // dim
    ierr = xf_Error(xf_Basis2Dim(Basis, &dim));
    if (ierr != xf_OK) return ierr;
    
    for (o=0; o<=5; o++){
      if (dim == 1)
	ierr = xf_Error(xf_QuadLine(o, &nq, &xq, &wq));
      else if (dim == 2)
	ierr = xf_Error(xf_QuadQuadrilateral(o, &nq, &xq, &wq));
      else
	ierr = xf_Error(xf_QuadHexahedron(o, &nq, &xq, &wq));
      xf_AssertEqual(ierr, xf_OK);

      for (i[0]=0; i[0]<=o; i[0]++){
	for (i[1]=0; i[1]<=o*(dim>=2); i[1]++){
	  for (i[2]=0; i[2]<=o*(dim>=3); i[2]++){
	    
	    // test shape F at quad points, where F = x^i0 * y^i1 *z^i2
	    ierr = xf_Error(xf_LagrangeNodes(Basis, o, NULL, xn, NULL));
	    xf_AssertEqual(ierr, xf_OK);
	    
	    for (d=0, nn=1; d<dim; d++) nn *= (o+1);
	    for (k=0; k<nn; k++)
	      for (d=0, F[k]=1.; d<dim; d++) F[k] *= xf_PowInt(xn[dim*k+d],i[d]);
	    
	    for (iq=0; iq<nq; iq++){

	      // shape
	      x = xq+dim*iq;
	      if (itype == 0)
		ierr = xf_Error(xf_Shape_TensorLagrange(dim, o, x, phi));
	      else
		ierr = xf_Error(xf_Shape_TensorLagrangeGauss(dim, o, x, phi));
	      xf_AssertEqual(ierr, xf_OK);
	      for (k=0, t=0.; k<nn; k++) t += F[k] * phi[k];
	      for (d=0, s=1.; d<dim; d++) s *= xf_PowInt(x[d],i[d]);	      
	      xf_AssertWithin(t, s, UTOL2);

	      // grad
	      x = xq+dim*iq;
	      if (itype == 0)
		ierr = xf_Error(xf_Grad_TensorLagrange(dim, o, x, gphi, nn));
	      else
		ierr = xf_Error(xf_Grad_TensorLagrangeGauss(dim, o, x, gphi, nn));
	      xf_AssertEqual(ierr, xf_OK);
	      for (d=0; d<dim; d++) gt[d] = 0;
	      for (k=0; k<nn; k++) 
		for (d=0; d<dim; d++) gt[d] += F[k] * gphi[nn*d+k];
	      for (d=0; d<dim; d++){
		for (k=0, s=i[d]; k<dim; k++)
		  s *= xf_PowInt(x[k],i[k]-(k==d));
		xf_AssertWithin(gt[d], s, UTOL2);
	      }

	      // Hess
	      x = xq+dim*iq;
	      if (itype == 0)
		ierr = xf_Error(xf_Hess_TensorLagrange(dim, o, x, hphi, nn));
	      else
		ierr = xf_Error(xf_Hess_TensorLagrangeGauss(dim, o, x, hphi, nn));
	      xf_AssertEqual(ierr, xf_OK);
	      for (d=0; d<dim*dim; d++) ht[d] = 0;
	      for (k=0; k<nn; k++) 
		for (d0=0; d0<dim; d0++)
		  for (d1=0; d1<dim; d1++) ht[d0*dim+d1] += F[k] * hphi[nn*(d0*dim+d1)+k];
	      for (d0=0; d0<dim; d0++)
		for (d1=0; d1<dim; d1++){
		  if (d0 != d1){
		    s = i[d0]*i[d1];
		    for (k=0; k<dim; k++)
		      s *= xf_PowInt(x[k],i[k]-((k==d0) || (k==d1)) );
		  }
		  else{
		    s = i[d0]*(i[d0]-1);
		    for (k=0; k<dim; k++)
		      s *= xf_PowInt(x[k],i[k]-2*(k==d0) );
		  }
		  d = d0*dim+d1;
		  xf_AssertWithin(ht[d], s, UTOL4);
		}

	    } //iq
	  } // i[2]
	} // i[1]
      } // i[0]
    } //o
  } // irun

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}


TEST_xf_EvalPhysicalGrad_Q1TriLagrange()
{
  int ierr;
  real xq[4] = {0.3, 0.3, 0.2, 0.2};
  real gPhi_true[12] = {-1,1,0,-1,1,0, -1,0,1,-1,0,1}; 
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;
  
  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, 0);

  // compute basis functions (and grads) if quad or basis or order changed
  PhiData  = NULL;
  ierr = xf_Error(xf_EvalBasis(xfe_TriLagrange, 1, xfe_True, 2, xq, xfb_All, &PhiData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(PhiData->nn, 3);
  xf_AssertEqual(PhiData->nq, 2);
  
  /* Compute geometry Jacobian */
  JData    = NULL;
  ierr = xf_Error(xf_ElemJacobian(Mesh, 0, 0, 2, xq, xfb_detJ | xfb_iJ, xfe_True, &JData));
  xf_AssertEqual(ierr, xf_OK);

  /* Evaluate physical gradients */
  ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(PhiData->gPhi, gPhi_true, 12, UTOL1);

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy mesh */
  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, 0);

  return xf_OK;  
}


TEST_xf_EvalPhysicalGrad_Q1TriLagrange_Scaled()
{
  int ierr, k, d;
  real xq[4] = {0.3, 0.3, 0.2, 0.2};
  real gPhi_true[12] = {0,2,-2,0,2,-2, -2,2,0,-2,2,0};
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;
  
  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, 0);

  for (k=0; k<Mesh->nNode; k++)
    for (d=0; d<Mesh->Dim; d++)
      Mesh->Coord[k][d] *= 0.5;

  // compute basis functions (and grads) if quad or basis or order changed
  PhiData  = NULL;
  ierr = xf_Error(xf_EvalBasis(xfe_TriLagrange, 1, xfe_True, 2, xq, xfb_All, &PhiData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(PhiData->nn, 3);
  xf_AssertEqual(PhiData->nq, 2);
  
  /* Compute geometry Jacobian */
  JData    = NULL;
  ierr = xf_Error(xf_ElemJacobian(Mesh, 0, 1, 2, xq, xfb_detJ | xfb_iJ, xfe_True, &JData));
  xf_AssertEqual(ierr, xf_OK);

  /* Evaluate physical gradients */
  ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(PhiData->gPhi, gPhi_true, 12, UTOL1);

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy mesh */
  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, 0);

  return xf_OK;
}

TEST_xf_EvalPhysicalGrad_Q2TriLagrange()
{
  int ierr;
  real xq[4] = {0.3, 0.3, 0.2, 0.2};
  real gPhi_true[12] = {-1,1,0,-1,1,0, -1,0,1,-1,0,1}; 
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;
  
  ierr = xf_Error(xf_UnitTwoQ2TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, 0);

  // compute basis functions (and grads) if quad or basis or order changed
  PhiData  = NULL;
  ierr = xf_Error(xf_EvalBasis(xfe_TriLagrange, 1, xfe_True, 2, xq, xfb_All, &PhiData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(PhiData->nn, 3);
  xf_AssertEqual(PhiData->nq, 2);
  
  /* Compute geometry Jacobian */
  JData    = NULL;
  ierr = xf_Error(xf_ElemJacobian(Mesh, 0, 0, 2, xq, xfb_detJ | xfb_iJ, xfe_True, &JData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(JData->nq, 2);

  /* Evaluate physical gradients */
  ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertRealVectorWithin(PhiData->gPhi, gPhi_true, 12, UTOL1);

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy mesh */
  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, 0);

  return xf_OK;  
}


TEST_xf_EvalPhysicalGrad_Q2TriLagrange_Scaled()
{
  int ierr, k, d;
  real xq[4] = {0.3, 0.3, 0.2, 0.2};
  real gPhi_true[12] = {0,2,-2,0,2,-2, -2,2,0,-2,2,0};
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;
  
  ierr = xf_Error(xf_UnitTwoQ2TriangleMesh(&Mesh, xfe_True));
  xf_AssertEqual(ierr, 0);

  for (k=0; k<Mesh->nNode; k++)
    for (d=0; d<Mesh->Dim; d++)
      Mesh->Coord[k][d] *= 0.5;

  // compute basis functions (and grads) if quad or basis or order changed
  PhiData  = NULL;
  ierr = xf_Error(xf_EvalBasis(xfe_TriLagrange, 1, xfe_True, 2, xq, xfb_All, &PhiData));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(PhiData->nn, 3);
  xf_AssertEqual(PhiData->nq, 2);
  
  /* Compute geometry Jacobian */
  JData    = NULL;
  ierr = xf_Error(xf_ElemJacobian(Mesh, 0, 1, 2, xq, xfb_detJ | xfb_iJ, xfe_True, &JData));
  xf_AssertEqual(ierr, xf_OK);

  /* Evaluate physical gradients */
  ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
  xf_AssertEqual(ierr, xf_OK);
  
  xf_AssertRealVectorWithin(PhiData->gPhi, gPhi_true, 12, UTOL1);

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  xf_AssertEqual(ierr, xf_OK);

  /* Destroy mesh */
  ierr = xf_Error(xf_DestroyMesh(Mesh));
  xf_AssertEqual(ierr, 0);

  return xf_OK;
}



TEST_xf_ElemMassMatrix()
{
  int ierr;
  real *MM, *iMM, fac;
  real MM0_true[1] = {0.5};
  real MM1_true[9] = {1./12., 1./24., 1./24., 1./24., 1./12., 1./24., 1./24., 1./24., 1./12.};
  real MM2_true[6] = {1./60., 0.0, -1./360., 0.0, -1./90., -1./360.}; // just first terms
  real iMM0_true[1] = {2.0};
  real iMM1_true[9] = {18,-6,-6,-6,18,-6,-6,-6,18};
  real iMM2_true[6] = {72,-3,12,-3,12,12};

  xf_All *All;
  
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  xf_AssertEqual(ierr, 0);

  ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&All->Mesh, xfe_False));
  xf_AssertEqual(ierr, 0);

  // p = 0
  ierr = xf_Error(xf_ElemMassMatrix(All, 0, 0, xfe_TriLagrange, 0, NULL, NULL, &MM, &fac));
  xf_AssertEqual(ierr, 0);
  xf_AssertWithin(fac, 1.0, UTOL1);
  xf_AssertRealVectorWithin(MM, MM0_true, 1, UTOL1);
  
  ierr = xf_Error(xf_ElemInvMassMatrix(All, 0, 0, xfe_TriLagrange, 0, NULL, NULL, &iMM, &fac));
  xf_AssertEqual(ierr, 0);
  xf_AssertWithin(fac, 1.0, UTOL1);
  xf_AssertRealVectorWithin(iMM, iMM0_true, 1, UTOL3);


  // p = 1
  ierr = xf_Error(xf_ElemMassMatrix(All, 0, 0, xfe_TriLagrange, 1, NULL, NULL, &MM, &fac));
  xf_AssertEqual(ierr, 0);
  xf_AssertWithin(fac, 1.0, UTOL1);
  xf_AssertRealVectorWithin(MM, MM1_true, 9, UTOL1);

  ierr = xf_Error(xf_ElemInvMassMatrix(All, 0, 0, xfe_TriLagrange, 1, NULL, NULL, &iMM, &fac));
  xf_AssertEqual(ierr, 0);
  xf_AssertWithin(fac, 1.0, UTOL1);
  xf_AssertRealVectorWithin(iMM, iMM1_true, 9, UTOL3);

  // p = 2
  ierr = xf_Error(xf_ElemMassMatrix(All, 0, 0, xfe_TriLagrange, 2, NULL, NULL, &MM, &fac));
  xf_AssertEqual(ierr, 0);
  xf_AssertWithin(fac, 1.0, UTOL1);
  xf_AssertRealVectorWithin(MM, MM2_true, 6, UTOL1);

  ierr = xf_Error(xf_ElemInvMassMatrix(All, 0, 0, xfe_TriLagrange, 2, NULL, NULL, &iMM, &fac));
  xf_AssertEqual(ierr, 0);
  xf_AssertWithin(fac, 1.0, UTOL1);
  xf_AssertRealVectorWithin(iMM, iMM2_true, 6, UTOL3);

  ierr = xf_Error(xf_DestroyAll(All));
  xf_AssertEqual(ierr, 0);

  return xf_OK;  
}


TEST_xf_ComputeTransferMatrix()
{
  int ierr;
  int k, p, nn;
  real T[3600], iT[3600];
  real U0[100], U1[100], U2[100];
  
  // Lagrange to Hierarchical tri basis, various orders
  for (p=1; p<=5; p++){
    nn = (p+1)*(p+2)/2; // rank
    // forward transfer matrix
    ierr = xf_Error(xf_ComputeTransferMatrix(xfe_TriLagrange, p, xfe_TriHierarch, 
					     p, xfe_Triangle, 0,  T));
    xf_AssertEqual(ierr, xf_OK);
    // inverse transfer matrix
    ierr = xf_Error(xf_ComputeTransferMatrix(xfe_TriHierarch, p, xfe_TriLagrange, 
					     p, xfe_Triangle, 0, iT));
    xf_AssertEqual(ierr, xf_OK);
    // random discrete data
    for (k=0; k<nn; k++) U0[k] = 1.0+k/100.0-k*k/10000.0;
    // apply T
    xf_MxV_Set(T,  U0, nn, nn, U1);
    xf_MxV_Set(iT, U1, nn, nn, U2);
    // check answer (U0 should be U2)
    xf_AssertRealVectorWithin(U0, U2, nn, UTOL3);

  } // p

  return xf_OK;  
}




TEST_xf_ProjectOnElemQR()
{
  int ierr, i[3], nq, iq, k, nn, d, o;
  int nBasis, itest, iBasis, itotal;
  enum xfe_BasisType *BasisVec;
  enum xfe_BasisType BTri[2]  = {xfe_TriLagrange,  xfe_TriHierarch};
  enum xfe_BasisType BTet[2]  = {xfe_TetLagrange,  xfe_TetHierarch};
  enum xfe_BasisType BQuad[2] = {xfe_QuadLagrange, xfe_QuadLagrangeGauss};
  enum xfe_BasisType BHex[2]  = {xfe_HexLagrange,  xfe_HexLagrangeGauss};
  enum xfe_BasisType Basis;
  enum xfe_Bool QuadChanged, TensorFlag;
  int Order;
  int egrp, elem, dim, sr;
  int *P;
  real *u, *ue, *xq, *U;
  real *Q, *R;
  xf_BasisData *PhiData, *GeomPhiData;
  xf_QuadData *QuadData;
  xf_Mesh *Mesh;
  xf_All *All;
  

  for (itest=0; itest<4; itest++){

    ierr = xf_Error(xf_CreateAll(&All, xfe_False));
    xf_AssertEqual(ierr, 0);

    switch (itest){
    case 0:
      ierr = xf_Error(xf_UnitTwoQ1TriangleMesh(&All->Mesh, xfe_False));
      xf_AssertEqual(ierr, 0);
      BasisVec = BTri; nBasis = 2; TensorFlag = xfe_False;
      break;
    case 1:
      ierr = xf_Error(xf_UnitQ2TetMesh(&All->Mesh, xfe_False));
      xf_AssertEqual(ierr, 0);
      BasisVec = BTet; nBasis = 2; TensorFlag = xfe_False;
      break;
    case 2:
      ierr = xf_Error(xf_UnitQ1QuadMesh(&All->Mesh, xfe_False));
      xf_AssertEqual(ierr, 0);
      BasisVec = BQuad; nBasis = 2; TensorFlag = xfe_True;
      break;
    case 3:
      ierr = xf_Error(xf_UnitQ2HexMesh(&All->Mesh, xfe_False));
      xf_AssertEqual(ierr, 0);
      BasisVec = BHex; nBasis = 2; TensorFlag = xfe_True;
      break;
    default:
      return xf_Error(xf_NOT_SUPPORTED);
      break;
    }

    Mesh = All->Mesh;
    dim = Mesh->Dim;
    sr  = 1;
    Order = 3;

    U           = NULL;
    PhiData     = NULL;
    QuadData    = NULL;
    u           = NULL;
    ue          = NULL;
    Q           = NULL;
    R           = NULL;
    P           = NULL;
    
    // loop over the bases
    for (iBasis=0; iBasis<nBasis; iBasis++){

      Basis = BasisVec[iBasis];

      for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

	// determine nn = # unknowns for elements in this group
	ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
	xf_AssertEqual(ierr, xf_OK);

	ierr = xf_Error(xf_ReAlloc( (void **) &U, nn*sr, sizeof(real)));
	xf_AssertEqual(ierr, xf_OK);

	for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
	  // max quad points on elem
	  ierr = xf_Error(xf_QuadElemAtLeast(Mesh, egrp, elem, 2*Order, nn, &QuadData, &QuadChanged));
	  xf_AssertEqual(ierr, xf_OK);
	  nq = QuadData->nquad;
	  xq = QuadData->xquad;
      
	  // reallocate memory
	  ierr = xf_Error(xf_ReAlloc( (void **) &u, nq*sr, sizeof(real)));
	  xf_AssertEqual(ierr, xf_OK);
      
	  ierr = xf_Error(xf_ReAlloc( (void **) &ue, nq*sr, sizeof(real)));
	  xf_AssertEqual(ierr, xf_OK);


	  // compute basis functions
	  ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged,  nq, xq, xfb_Phi, &PhiData));
	  xf_AssertEqual(ierr, xf_OK);

	  // set ue to a function of the appropriate order and check its projection
	  o = Order;
	  itotal = 0;
	  for (i[0]=0; i[0]<=o; i[0]++){
	    for (i[1]=0; i[1]<=o-i[0]*(!TensorFlag); i[1]++){
	      for (i[2]=0; i[2]<=(o-(i[0]+i[1])*(!TensorFlag))*(dim==3); i[2]++){
	    
		// set ue to a monomial
		for (iq=0; iq<nq; iq++)
		  for (d=0, ue[iq]=1.; d<dim; d++) ue[iq] *= xf_PowInt(xq[dim*iq+d],i[d]);

		// Project ue to get coefficients U
		ierr = xf_Error(xf_ProjectOnElemQR(Mesh, egrp, elem, sr, QuadData, 
						   Basis, Order, (itotal==0), &Q, &R,
						   &P, ue, U));
		xf_AssertEqual(ierr, xf_OK);

		// Calculate U at xq -> u
		xf_MxM_Set(PhiData->Phi, U, nq, PhiData->nn, sr, u);
	    
		// check if got IC exactly
		for (iq=0; iq<nq; iq++)
		  xf_AssertWithin(u[iq], ue[iq], UTOL3);
	        
		itotal++;
	      } // i[2]
	    } // i[1]
	  } // i[0]
	} // elem
      } // egrp
    } // iBasis


    // Only destroy QuadData if points are generic
    ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
    xf_AssertEqual(ierr, xf_OK);

    /* Destroy Basis Data */
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    xf_AssertEqual(ierr, xf_OK);

    xf_Release((void *) U);
    xf_Release((void *) u);
    xf_Release((void *) ue);
    xf_Release((void *) Q);
    xf_Release((void *) R);
    xf_Release((void *) P);

  
    ierr = xf_Error(xf_DestroyAll(All));
    xf_AssertEqual(ierr, 0);
  } // itest

  return xf_OK;  
}

