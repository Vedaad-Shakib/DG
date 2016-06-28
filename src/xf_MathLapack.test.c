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

TEST_xf_EigSymTriDiag()
{
  int ierr, k;
  real D[5] = {1, 2, 1, 2, 1};
  real E[4] = {2, -1, 3, 1};
  real Z[25];
  real D0[5] = { -1.913960286522578, // From Matlab
		 -0.460835322521163,
		 1.000000000000000,
		 3.460835322521164,
		 4.913960286522576};
  
  ierr = xf_EigSymTriDiag(5, D, E, Z);
  xf_AssertEqual(ierr, 0);
  
  xf_AssertRealVectorWithin(D, D0, 5, UTOL1);

  /*   for (k=0; k<25; k++){ */
  /*     if (k % 5 == 0) xf_printf("\n"); */
  /*     xf_printf("%.6E ", Z[k]); */
  /*   } */
  /*   xf_printf("\n"); */
  
  return xf_OK;
}


TEST_xf_EigSym()
{
  int ierr, k;
  real A[9] = {4,2,1, 2,4,1, 1,1,9};
  real E[3];
  real E0[3] = {2.000000000000000,
		5.438447187191166,
		9.561552812808831};
  
  ierr = xf_EigSym(3, A, E);
  xf_AssertEqual(ierr, 0);
  
  xf_AssertRealVectorWithin(E, E0, 3, UTOL1);

  /*   for (k=0; k<9; k++){ */
  /*     if (k % 3 == 0) xf_printf("\n"); */
  /*     xf_printf("%.10E ", A[k]); */
  /*   } */
  /*   xf_printf("\n"); */
  
  return xf_OK;
}


TEST_xf_CholDecomp()
{
  int ierr, k;
  real A[9] = {4,2,1, 2,4,1, 1,1,9};
  real U[9] = {2.000000000000000,1.000000000000000,0.500000000000000,
               2,   1.732050807568877,   0.288675134594813,
               1,                   1,   2.943920288775949}; 


  ierr = xf_CholDecomp(3, A, xfe_True);
  xf_AssertEqual(ierr, 0);
  
  xf_AssertRealVectorWithin(A, U, 9, UTOL1);

  /* for (k=0; k<9; k++){ */
  /*   if (k % 3 == 0) xf_printf("\n"); */
  /*   xf_printf("%.10E ", A[k]); */
  /* } */
  /* xf_printf("\n"); */
  
  return xf_OK;
}

