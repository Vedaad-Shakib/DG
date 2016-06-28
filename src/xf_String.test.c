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

TEST_xf_ScanXInt()
{
  int ierr, n;
  char line[80] = "1 -5 0 ";
  char line2[80] = "";
  int v[3];
  int v_true[3] = {1, -5, 0};


  ierr = xf_Error(xf_ScanXInt(line, &n, v));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(n, 3);
  xf_AssertIntVectorEqual(v, v_true, 3);
  
  ierr = xf_Error(xf_ScanXInt(line2, &n, v));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(n, 0);

  return xf_OK;
}


TEST_xf_ScanXReal()
{
  int ierr, n;
  char line[80] = "0.1 1e-2 3.4 2";
  char line2[80] = "";
  char line3[80] = "0.18 0.04 18.849555921538759 0.0";
  real v[4];
  real v_true[4] = {0.1, 0.01, 3.4, 2.0};
  real v_true3[4] = {0.18, 0.04, 18.849555921538759, 0.0};


  ierr = xf_Error(xf_ScanXReal(line, &n, v));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(n, 4);
  xf_AssertRealVectorWithin(v, v_true, 4, UTOL0);
  
  ierr = xf_Error(xf_ScanXReal(line2, &n, v));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(n, 0);

  ierr = xf_Error(xf_ScanXReal(line3, &n, v));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(n, 4);
  xf_AssertRealVectorWithin(v, v_true3, 4, UTOL0);

  return xf_OK;
}


TEST_xf_ScanXStringAlloc()
{
  int ierr, n;
  char line[80] = "  this is  only a test  123 ";
  char **Strings;

  ierr = xf_Error(xf_ScanXStringAlloc(line, 20, &n, &Strings));
  xf_AssertEqual(ierr, xf_OK);
  xf_AssertEqual(n, 6);
  
  xf_AssertEqual(strncmp(Strings[0], "this", 4), 0);
  xf_AssertEqual(strncmp(Strings[1], "is",   2), 0);
  xf_AssertEqual(strncmp(Strings[2], "only", 4), 0);
  xf_AssertEqual(strncmp(Strings[3], "a",    1), 0);
  xf_AssertEqual(strncmp(Strings[4], "test", 4), 0);
  xf_AssertEqual(strncmp(Strings[5], "123",  3), 0);
  
  xf_Release2( (void **) Strings);

  return xf_OK;
}


TEST_xf_PickRealsFromLine()
{
  int ierr;
  char line[80] = "  a blah cat dog";
  char keys[80] = " dog blah  ";
  real data[4] = {1., 2., 3., 4.};
  real values0[2] = {4., 2.};
  real values[2];

  ierr = xf_Error(xf_PickRealsFromLine(line, keys, data, values));
  xf_AssertEqual(ierr, xf_OK);

  xf_AssertRealVectorWithin(values, values0, 2, UTOL0);

  return xf_OK;
}
