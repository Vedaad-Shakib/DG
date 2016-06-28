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


/* static helper functions */
static real
mypow(real x, int n){
  int k;
  real t;
  t = 1.0;
  for (k=0; k<n; k++) t *= x;
  return t;
}


TEST_xf_QuadLine()
{
  int ierr, o, iq, nq;
  real sum, expected, *xq, *wq;

  xq = NULL;
  wq = NULL;

  for (o=0; o<40; o++){

    ierr = xf_Error(xf_QuadLine(o, &nq, &xq, &wq));
    xf_AssertEqual(ierr, xf_OK);

    // int x^n
    expected = 1.0/((real) o + 1.0);

    sum = 0.;
    for (iq=0; iq<nq; iq++)
      sum += wq[iq] * mypow(xq[iq],o);

    xf_AssertWithin(sum, expected, UTOL1);

  } //o

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}


TEST_xf_QuadTriangle()
{
  int ierr, o, iq, nq;
  real sum, expected, *xq, *wq;

  xq = NULL;
  wq = NULL;

  for (o=1; o<42; o++){

    ierr = xf_Error(xf_QuadTriangle(o, &nq, &xq, &wq));
    xf_AssertEqual(ierr, xf_OK);

    // int(int(x*y^(o-1),0,1-y), y, 0, 1)
    expected = 0.5/((real) o) - 1.0/((real) o + 1.0) + 0.5/((real) o + 2.0);

    sum = 0.;
    for (iq=0; iq<nq; iq++)
      sum += wq[iq] * xq[2*iq+0] * mypow(xq[2*iq+1],o-1);

    xf_AssertWithin(sum, expected, UTOL1); 
    
    sum = 0.;
    for (iq=0; iq<nq; iq++)
      sum += wq[iq] * xq[2*iq+1] * mypow(xq[2*iq+0],o-1);
    xf_AssertWithin(sum, expected, UTOL1);

  } //o

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}


TEST_xf_QuadTetrahedron()
{
  int ierr, o, iq, nq;
  real r, sum, expected, *xq, *wq;

  xq = NULL;
  wq = NULL;

  for (o=2; o<11; o++){

    ierr = xf_Error(xf_QuadTetrahedron(o, &nq, &xq, &wq));
    xf_AssertEqual(ierr, xf_OK);

    // int(int(int(x*y*z^(o-2),0,1-y-z), y, 0, 1-z),z,0,1)
    r = ((real) o) - 1.0;
    expected = 1./24.*(1./r-4./(r+1.)+6./(r+2.)-4./(r+3.)+1./(r+4.));

    sum = 0.;
    for (iq=0; iq<nq; iq++)
      sum += wq[iq] * xq[3*iq+0] * xq[3*iq+1] * mypow(xq[3*iq+2],o-2);

    xf_AssertWithin(sum, expected, UTOL1);

  } //o

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}

TEST_xf_QuadQuadrilateral()
{
  int ierr, o, iq, nq;
  real sum, expected, *xq, *wq;

  xq = NULL;
  wq = NULL;

  for (o=1; o<40; o++){

    ierr = xf_Error(xf_QuadQuadrilateral(o, &nq, &xq, &wq));
    xf_AssertEqual(ierr, xf_OK);

    // int(int(x^o*y^o,0,1), y, 0, 1)
    expected = (1.0/((real) o + 1.0));
    expected *= expected;

    for (iq=0, sum=0.; iq<nq; iq++)
      sum += wq[iq] * mypow(xq[2*iq+0],o) * mypow(xq[2*iq+1],o);
    xf_AssertWithin(sum, expected, UTOL1);

  } //o

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}

TEST_xf_QuadHexahedron()
{
  int ierr, o, iq, nq;
  real sum, expected, fac, *xq, *wq;

  xq = NULL;
  wq = NULL;

  for (o=1; o<40; o++){

    ierr = xf_Error(xf_QuadHexahedron(o, &nq, &xq, &wq));
    xf_AssertEqual(ierr, xf_OK);

    // int(int(int(x^o*y^o*z^o,0,1), y, 0, 1),z,0,1)
    fac = (1.0/((real) o + 1.0));
    expected = fac*fac*fac;

    for (iq=0, sum=0.; iq<nq; iq++)
      sum += wq[iq] * mypow(xq[3*iq+0],o) * mypow(xq[3*iq+1],o) * mypow(xq[3*iq+2],o);
    xf_AssertWithin(sum, expected, UTOL1);

  } //o

  xf_Release( (void *) xq);
  xf_Release( (void *) wq);
  
  return xf_OK;  
}
