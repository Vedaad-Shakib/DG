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


TEST_xf_MxM()
{
  real A[6] = {1, 2, 3, 4, 5, 6};
  real B[6] = {1, 2, 1, 3, 1, 2};
  real C[4] = {2, 3, -1, -10};
  real C0[4] = {-4, -11, -16, -45};
  real C1[4] = {2, 3, -1, -10};
  real C2[4] = {6, 14, 15, 35};
  real C3[4] = {-6, -14, -15, -35};

  xf_MxM_Sub(A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C0, 4, UTOL0);

  xf_MxM_Add(A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C1, 4, UTOL0);

  xf_MxM_Set(A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C2, 4, UTOL0);

  xf_MxM_Neg(A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C3, 4, UTOL0);

  return xf_OK;  
}


TEST_xf_MTxM_Add()
{
  real A[6] = {1, 4, 2, 5, 3, 6};
  real B[6] = {1, 2, 2, 1, 3, 1};
  real C[4] = {1, 1, 0, 0};
  real C_true[4] = {15, 8, 32, 19};

  xf_MTxM_Add(A, B, 2, 3, 2, C);
  xf_AssertRealVectorWithin(C, C_true, 4, UTOL0);

  return xf_OK;  
}

TEST_xf_MTxM_Sub()
{
  real A[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  real B[12] = {2, 1, 3, -1, 2, 4, 3, 2, 1, -1, -3, -2};
  real C[6] = {1, 2, 3, 4, 5, 6};
  real C_true[6] = {-6, 6, -3, -6, 7, -6};

  xf_MTxM_Sub(A, B, 2, 4, 3, C);
  xf_AssertRealVectorWithin(C, C_true, 6, UTOL0);

  return xf_OK;  
}


TEST_xf_MTxM_Set()
{
  real A[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  real B[12] = {2, 1, 3, -1, 2, 4, 3, 2, 1, -1, -3, -2};
  real C[6] = {1, 2, 3, 4, 5, 6};
  real C_true[6] = {7, -4, 6, 10, -2, 12};

  xf_MTxM_Set(A, B, 2, 4, 3, C);
  xf_AssertRealVectorWithin(C, C_true, 6, UTOL0);

  return xf_OK;  
}

TEST_xf_MTxM_Neg()
{
  real A[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  real B[12] = {2, 1, 3, -1, 2, 4, 3, 2, 1, -1, -3, -2};
  real C[6] = {1, 2, 3, 4, 5, 6};
  real C_true[6] = {-7, 4, -6, -10, 2, -12};

  xf_MTxM_Neg(A, B, 2, 4, 3, C);
  xf_AssertRealVectorWithin(C, C_true, 6, UTOL0);

  return xf_OK;  
}

