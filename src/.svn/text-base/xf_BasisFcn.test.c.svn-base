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
#include "xf_AllStruct.h"
#include "xf_QuadRule.h"
#include "xf_Memory.h"
#include "xf_Basis.h"
#include "xf_Math.h"


TEST_xf_Shape_TriLagrange()
{
  // Tests whether U(Order) -> U(Order+1) -> U(Order) is the identity
  int ierr, o, k, nn, nn1;
  int iBasis;
  enum xfe_BasisType Basis, VBasis[2] = {xfe_TriLagrange, xfe_TriHierarch};
  real F0[36], F[36], F1[49], T[1764];

  for (iBasis=0; iBasis<2; iBasis++){
    
    Basis = VBasis[iBasis];

    for (o=0; o<=6; o++){

      // number of dofs at Order o
      nn = (o+1)*(o+2)/2;

      // number of dofs at Order o+
      nn1 = (o+2)*(o+3)/2;

      // Initialize F0 to some pseudo-random values
      for (k=0; k<nn; k++) F0[k] = (k*k + k%3)%17 - 7;

      // compute transfer matrix from o to o+1
      ierr = xf_Error(xf_ComputeTransferMatrix(Basis, o, Basis, o+1, xfe_Triangle, 0, T));
      xf_AssertEqual(ierr, xf_OK);

      // transfer F0 from o to o+1 -> F1
      xf_MxV_Set(T, F0, nn1, nn, F1);

      // compute transfer matrix from o+1 to o
      ierr = xf_Error(xf_ComputeTransferMatrix(Basis, o+1, Basis, o, xfe_Triangle, 0, T));
      xf_AssertEqual(ierr, xf_OK);
      
      // transfer F1 from o+1 to o -> F
      xf_MxV_Set(T, F1, nn, nn1, F);

      // F0 and F should be the same
      xf_AssertRealVectorWithin(F, F0, nn, 100.0*UTOL5);
      
    } // o

  } // iBasis
  
  return xf_OK;  
}



TEST_xf_Shape_TriLagrangeHierarch()
{
  int ierr, o, iq, nq, k, nn, a, b;
  int IsHierarch;
  real *xq, *wq, *x, t;
  real phi[49], F[49], F1[49], T[2401], xn[98];

  xq = NULL;
  wq = NULL;

  for (IsHierarch=0; IsHierarch<2; IsHierarch++){

    for (o=0; o<=6; o++){

      ierr = xf_Error(xf_QuadTriangle(2*o, &nq, &xq, &wq));
      xf_AssertEqual(ierr, xf_OK);
      
      if (IsHierarch){
        // compute transfer matrix from Lagrange to Hierarch
        ierr = xf_Error(xf_ComputeTransferMatrix(xfe_TriLagrange, o, xfe_TriHierarch, 
                                                 o, xfe_Triangle, 0, T));
        xf_AssertEqual(ierr, xf_OK);
      }


      for (a=0; a<=o; a++){
	b = o-a;

	// test shape F at quad points, where F = x^a * y^b
	ierr = xf_Error(xf_LagrangeNodes(xfe_TriLagrange, o, NULL, xn, NULL));
	xf_AssertEqual(ierr, xf_OK);

	// fill in Lagrange values
	nn = (o+1)*(o+2)/2;
	for (k=0; k<nn; k++) 
	  F[k] = xf_PowInt(xn[2*k+0],a)*xf_PowInt(xn[2*k+1],b);

	// map to hierarchical basis if necessary
	if (IsHierarch){
	  xf_MxV_Set(T, F, nn, nn, F1);
	  for (k=0; k<nn; k++) F[k] = F1[k];
	}

	for (iq=0; iq<nq; iq++){
	
	  if (IsHierarch){ 
	    ierr = xf_Error(xf_Shape_TriHierarch(o, xq+2*iq, phi));
	    xf_AssertEqual(ierr, xf_OK);
	  }
	  else{
	    ierr = xf_Error(xf_Shape_TriLagrange(o, xq+2*iq, phi));
	    xf_AssertEqual(ierr, xf_OK);
	  }

	  t = 0.;
	  for (k=0; k<nn; k++) t += F[k] * phi[k];
	  x = xq+2*iq;
	  xf_AssertWithin(t, xf_PowInt(x[0],a) * xf_PowInt(x[1],b), UTOL2);
	}
      } // a
    } //o

  } // IsHierarch

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}


TEST_xf_Shape_TetLagrangeHierarch()
{
  int ierr, o, iq, nq, k, nn, a, b, c;
  int IsHierarch;
  real *xq, *wq, *x, t;
  real phi[56], F[56], F1[56], T[3136], xn[168];

  xq = NULL;
  wq = NULL;

  for (IsHierarch=0; IsHierarch<2; IsHierarch++){

    for (o=0; o<=5; o++){

      ierr = xf_Error(xf_QuadTetrahedron(o, &nq, &xq, &wq));
      xf_AssertEqual(ierr, xf_OK);

      if (IsHierarch){
	// compute transfer matrix from Lagrange to Hierarch
	ierr = xf_Error(xf_ComputeTransferMatrix(xfe_TetLagrange, o, xfe_TetHierarch, 
						 o, xfe_Triangle, 0, T));
	xf_AssertEqual(ierr, xf_OK);
      }

      for (a=0; a<=o; a++){
	for (b=0; b<=(o-a); b++){
	  c = o-a-b;

	  // test shape F at quad points, where F = x^a * y^b *z^c
	  ierr = xf_Error(xf_LagrangeNodes(xfe_TetLagrange, o, NULL, xn, NULL));
	  xf_AssertEqual(ierr, xf_OK);

	  // fill in Lagrange values
	  nn = (o+1)*(o+2)*(o+3)/6;	
	  for (k=0; k<nn; k++) 
	    F[k] = xf_PowInt(xn[3*k+0],a)*xf_PowInt(xn[3*k+1],b)*xf_PowInt(xn[3*k+2],c);

	  // map to hierarchical basis if necessary
	  if (IsHierarch){
	    xf_MxV_Set(T, F, nn, nn, F1);
	    for (k=0; k<nn; k++) F[k] = F1[k];
	  }

	  for (iq=0; iq<nq; iq++){

	    if (IsHierarch){ 
	      ierr = xf_Error(xf_Shape_TetHierarch(o, xq+3*iq, phi));
	      xf_AssertEqual(ierr, xf_OK);
	    }
	    else{
	      ierr = xf_Error(xf_Shape_TetLagrange(o, xq+3*iq, phi));
	      xf_AssertEqual(ierr, xf_OK);
	    }
	    
	    t = 0.;
	    for (k=0; k<nn; k++) t += F[k] * phi[k];
	    x = xq+3*iq;
	    xf_AssertWithin(t, xf_PowInt(x[0],a) * xf_PowInt(x[1],b) * xf_PowInt(x[2],c), UTOL3);
	  } //iq

	} // b
      } //a
    } //o

  } // IsHierarch
    
  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}



TEST_xf_Grad_LineLagrange()
{
  int ierr, o, iq, nq, i, nn;
  real sum, expected, *xq, *wq;
  real gphi[6], t, F[6];

  xq = NULL;
  wq = NULL;

  for (o=1; o<=5; o++){

    ierr = xf_Error(xf_QuadLine(o, &nq, &xq, &wq));
    xf_AssertEqual(ierr, xf_OK);

    // test int (dF/dx), where F = x^o
    expected = 1.0;
    
    nn = o + 1;
    for (i=0; i<nn; i++){  // Lagrange basis
      t = ( (real) i) / ( (real) (nn-1));
      F[i] = xf_PowInt(t, o);
    } // i

    sum = 0.;
    for (iq=0; iq<nq; iq++){
      ierr = xf_Error(xf_Grad_LineLagrange(o, xq[iq], gphi));
      xf_AssertEqual(ierr, xf_OK);
      for (i=0; i<nn; i++) sum += wq[iq]*F[i] * gphi[i];
    }

    xf_AssertWithin(sum, expected, UTOL1);

  } //o

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}


TEST_xf_Grad_TriLagrangeHierarch()
{
  int ierr, o, iq, nq, k, nn, a, b;
  int IsHierarch;
  real *xq, *wq, *x, tx, ty;
  real gphi[42], F[21], F1[21], T[441], xn[42];

  xq = NULL;
  wq = NULL;

  for (IsHierarch=0; IsHierarch<2; IsHierarch++){

    for (o=0; o<=5; o++){
    
      ierr = xf_Error(xf_QuadTriangle(2*o, &nq, &xq, &wq));
      xf_AssertEqual(ierr, xf_OK);

      if (IsHierarch){
	// compute transfer matrix from Lagrange to Hierarch
	ierr = xf_Error(xf_ComputeTransferMatrix(xfe_TriLagrange, o, xfe_TriHierarch, 
						 o, xfe_Triangle, 0, T));
	xf_AssertEqual(ierr, xf_OK);
      }

      for (a=0; a<=o; a++){
	b = o-a;

	// test grad F at quad points, where F = x^a * y^b
	ierr = xf_Error(xf_LagrangeNodes(xfe_TriLagrange, o, NULL, xn, NULL));
	xf_AssertEqual(ierr, xf_OK);

	// fill in Lagrange values
	nn = (o+1)*(o+2)/2;
	for (k=0; k<nn; k++) 
	  F[k] = xf_PowInt(xn[2*k+0],a)*xf_PowInt(xn[2*k+1],b);

	// map to hierarchical basis if necessary
	if (IsHierarch){
	  xf_MxV_Set(T, F, nn, nn, F1);
	  for (k=0; k<nn; k++) F[k] = F1[k];
	}

	for (iq=0; iq<nq; iq++){

	  if (IsHierarch){ 
	    ierr = xf_Error(xf_Grad_TriHierarch(o, xq+2*iq, gphi, nn));
	    xf_AssertEqual(ierr, xf_OK);
	  }
	  else{
	    ierr = xf_Error(xf_Grad_TriLagrange(o, xq+2*iq, gphi, nn));
	    xf_AssertEqual(ierr, xf_OK);
	  }

	  tx = ty = 0.;
	  for (k=0; k<nn; k++){
	    tx += F[k] * gphi[ 0+k];
	    ty += F[k] * gphi[nn+k];
	  }
	  x = xq+2*iq;
	  xf_AssertWithin(tx, a*xf_PowInt(x[0],a-1) *   xf_PowInt(x[1],b  ), UTOL3);
	  xf_AssertWithin(ty,   xf_PowInt(x[0],a  ) * b*xf_PowInt(x[1],b-1), UTOL3);      
	}
      } // a
    } //o

  }

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}


TEST_xf_Grad_TetLagrangeHierarch()
{
  int ierr, o, iq, nq, k, nn, a, b, c;
  int IsHierarch;
  real *xq, *wq, *x, tx, ty, tz, t;
  real gphi[168], F[56], F1[56], T[3136], xn[168];

  xq = NULL;
  wq = NULL;

  for (IsHierarch=0; IsHierarch<2; IsHierarch++){

    for (o=0; o<=5; o++){

      ierr = xf_Error(xf_QuadTetrahedron(o, &nq, &xq, &wq));
      xf_AssertEqual(ierr, xf_OK);

      if (IsHierarch){
	// compute transfer matrix from Lagrange to Hierarch
	ierr = xf_Error(xf_ComputeTransferMatrix(xfe_TetLagrange, o, xfe_TetHierarch, 
						 o, xfe_Triangle, 0, T));
	xf_AssertEqual(ierr, xf_OK);
      }
    
      for (a=0; a<=o; a++){
	for (b=0; b<=(o-a); b++){
	  c = o-a-b;

	  // test grad F at quad points, where F = x^a * y^b *z^c
	  ierr = xf_Error(xf_LagrangeNodes(xfe_TetLagrange, o, NULL, xn, NULL));
	  xf_AssertEqual(ierr, xf_OK);

	  // fill in Lagrange values
	  nn = (o+1)*(o+2)*(o+3)/6;	
	  for (k=0; k<nn; k++) 
	    F[k] = xf_PowInt(xn[3*k+0],a)*xf_PowInt(xn[3*k+1],b)*xf_PowInt(xn[3*k+2],c);

	  // map to hierarchical basis if necessary
	  if (IsHierarch){
	    xf_MxV_Set(T, F, nn, nn, F1);
	    for (k=0; k<nn; k++) F[k] = F1[k];
	  }

	  for (iq=0; iq<nq; iq++){
	    if (IsHierarch){
	      ierr = xf_Error(xf_Grad_TetHierarch(o, xq+3*iq, gphi, nn));
	      xf_AssertEqual(ierr, xf_OK);
	    }
	    else{
	      ierr = xf_Error(xf_Grad_TetLagrange(o, xq+3*iq, gphi, nn));
	      xf_AssertEqual(ierr, xf_OK);
	    }

	    tx = ty = tz = 0.;
	    for (k=0; k<nn; k++){
	      tx += F[k] * gphi[   0+k];
	      ty += F[k] * gphi[  nn+k];
	      tz += F[k] * gphi[2*nn+k];
	    }
	    x = xq+3*iq;
	    xf_AssertWithin(tx, a*xf_PowInt(x[0],a-1) *   xf_PowInt(x[1],b  ) *   xf_PowInt(x[2],c  ), UTOL4);
	    xf_AssertWithin(ty,   xf_PowInt(x[0],a  ) * b*xf_PowInt(x[1],b-1) *   xf_PowInt(x[2],c  ), UTOL4); 
	    xf_AssertWithin(tz,   xf_PowInt(x[0],a  ) *   xf_PowInt(x[1],b  ) * c*xf_PowInt(x[2],c-1), UTOL4); 
	  } //iq

	} // b
      } //a
    } //o

  }
    
  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}





TEST_xf_Hess_TriLagrange()
{
  int ierr, o, iq, nq, k, nn, a, b;
  real *xq, *wq, *x, txx, txy, tyx, tyy;
  real hphi[84], F[21], xn[42];

  xq = NULL;
  wq = NULL;

  for (o=0; o<=5; o++){

    ierr = xf_Error(xf_QuadTriangle(2*o, &nq, &xq, &wq));
    xf_AssertEqual(ierr, xf_OK);

    for (a=0; a<=o; a++){
      b = o-a;

      // test hessian of F at quad points, where F = x^a * y^b
      ierr = xf_Error(xf_LagrangeNodes(xfe_TriLagrange, o, NULL, xn, NULL));
      xf_AssertEqual(ierr, xf_OK);

      nn = (o+1)*(o+2)/2;
      for (k=0; k<nn; k++) 
	F[k] = xf_PowInt(xn[2*k+0],a)*xf_PowInt(xn[2*k+1],b);

      for (iq=0; iq<nq; iq++){
	ierr = xf_Error(xf_Hess_TriLagrange(o, xq+2*iq, hphi, nn));
	xf_AssertEqual(ierr, xf_OK);
	txx = txy = tyx = tyy = 0.;
	for (k=0; k<nn; k++){
	  txx += F[k] * hphi[   0+k];
	  txy += F[k] * hphi[1*nn+k];
	  tyx += F[k] * hphi[2*nn+k];
	  tyy += F[k] * hphi[3*nn+k];
	}
	x = xq+2*iq;
	xf_AssertWithin(txx, a*(a-1)*xf_PowInt(x[0],a-2) *         xf_PowInt(x[1],b  ), UTOL3);
	xf_AssertWithin(txy,       a*xf_PowInt(x[0],a-1) *       b*xf_PowInt(x[1],b-1), UTOL3);
	xf_AssertWithin(tyx,                                                       txy, UTOL3);
	xf_AssertWithin(tyy,         xf_PowInt(x[0],a  ) * b*(b-1)*xf_PowInt(x[1],b-2), UTOL3);      
      }
    } // a
  } //o

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}



TEST_xf_Hess_TetLagrange()
{
  int ierr, o, iq, nq, k, nn, a, b, c;
  real *xq, *wq, *x;
  real txx, txy, txz, tyx, tyy, tyz, tzx, tzy, tzz;
  real hphi[504], F[56], xn[168];

  xq = NULL;
  wq = NULL;

  for (o=0; o<=5; o++){

    ierr = xf_Error(xf_QuadTetrahedron(o, &nq, &xq, &wq));
    xf_AssertEqual(ierr, xf_OK);

    for (a=0; a<=o; a++){
      for (b=0; b<=(o-a); b++){
	c = o-a-b;

	// test shape F at quad points, where F = x^a * y^b *z^c
	ierr = xf_Error(xf_LagrangeNodes(xfe_TetLagrange, o, NULL, xn, NULL));
	xf_AssertEqual(ierr, xf_OK);

	nn = (o+1)*(o+2)*(o+3)/6;	
	for (k=0; k<nn; k++) 
	  F[k] = xf_PowInt(xn[3*k+0],a)*xf_PowInt(xn[3*k+1],b)*xf_PowInt(xn[3*k+2],c);

	for (iq=0; iq<nq; iq++){
	  ierr = xf_Error(xf_Hess_TetLagrange(o, xq+3*iq, hphi, nn));
	  xf_AssertEqual(ierr, xf_OK);
	  txx = txy = txz = 0.;
	  tyx = tyy = tyz = 0.;
	  tzx = tzy = tzz = 0.;
	  for (k=0; k<nn; k++){
	    txx += F[k] * hphi[   0+k];
	    txy += F[k] * hphi[1*nn+k];
	    txz += F[k] * hphi[2*nn+k];
	    tyx += F[k] * hphi[3*nn+k];
	    tyy += F[k] * hphi[4*nn+k];
	    tyz += F[k] * hphi[5*nn+k];
	    tzx += F[k] * hphi[6*nn+k];
	    tzy += F[k] * hphi[7*nn+k];
	    tzz += F[k] * hphi[8*nn+k];
	  }
	  x = xq+3*iq;
	  xf_AssertWithin(txx, a*(a-1)*xf_PowInt(x[0],a-2) *         xf_PowInt(x[1],b  ) *         xf_PowInt(x[2],c  ), UTOL3);
	  xf_AssertWithin(txy,       a*xf_PowInt(x[0],a-1) *       b*xf_PowInt(x[1],b-1) *         xf_PowInt(x[2],c  ), UTOL3); 
	  xf_AssertWithin(txz,       a*xf_PowInt(x[0],a-1) *         xf_PowInt(x[1],b  ) *       c*xf_PowInt(x[2],c-1), UTOL3); 
	  xf_AssertWithin(tyx,                                                                                     txy, UTOL3); 
	  xf_AssertWithin(tyy,         xf_PowInt(x[0],a  ) * b*(b-1)*xf_PowInt(x[1],b-2) *         xf_PowInt(x[2],c  ), UTOL3);
	  xf_AssertWithin(tyz,         xf_PowInt(x[0],a  ) *       b*xf_PowInt(x[1],b-1) *       c*xf_PowInt(x[2],c-1), UTOL3); 
	  xf_AssertWithin(tzx,                                                                                     txz, UTOL3); 
	  xf_AssertWithin(tzy,                                                                                     tyz, UTOL3); 
	  xf_AssertWithin(tzz,         xf_PowInt(x[0],a  ) *         xf_PowInt(x[1],b  ) * c*(c-1)*xf_PowInt(x[2],c-2), UTOL3); 
	} //iq

      } // b
    } //a
  } //o

    
  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}

