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



TEST_xf_EigIterLanczos()
{
  int ierr, i, k;
  int nev = 5, ncv = 10, r = 16;
  int nelem[1] = {1};
  enum xfe_Bool done;
  enum xfe_EigStatusType Status;
  xf_Vector *V, *W;
  xf_VectorSet *VS, *EV;
  xf_EigSolverData *EigSolverData;
  real *rV, *rW;
  real E[11];
  real E0[5] = {3.205269272758331E+00, 3.478017834441318E+00,
		3.700434271459227E+00, 3.864944458808710E+00,
		3.965946199367802E+00}; // from Matlab
  real EV0[16] = {2.737176511998747E-01,  3.299034737043139E-01,
		  1.239048691826215E-01,  -1.805647421315387E-01,
		  -3.415340048050505E-01,  -2.310756993106725E-01,
		  6.302556436317688E-02,  3.070384761888058E-01,
		  3.070384756288038E-01,  6.302556503096360E-02,
		  -2.310756998114440E-01,  -3.415340048469053E-01,
		  -1.805647419061531E-01,  1.239048692632437E-01,
		  3.299034737712969E-01,  2.737176507833078E-01};    
/*   2.737176118030529E-01, 3.299032959518075E-01, */
/* 		  1.239048788823580E-01, -1.805645312063273E-01, */
/* 		  -3.415339241533843E-01,-2.310758234041150E-01, */
/* 		  6.302537068805480E-02, 3.070385627146517E-01, */
/* 		  3.070385623131781E-01, 6.302567508542938E-02, */
/* 		  -2.310758297424565E-01, -3.415340503347191E-01, */
/* 		  -1.805646934295058E-01, 1.239048340564260E-01, */
/* 		  3.299034791083640E-01, 2.737177350382031E-01}; */
  // EV0 is the eigenvector corresponding to eig 3.205...

  // create a vectorset, VS to be used for the Lanczos vectors
  ierr = xf_Error(xf_UnitrVectorSet(ncv+1, 1, nelem, &r, &VS));
  xf_AssertEqual(ierr, xf_OK);

  
  // create a vectorset, EV for the eigenvectors
  ierr = xf_Error(xf_UnitrVectorSet(nev, 1, nelem, &r, &EV));
  xf_AssertEqual(ierr, xf_OK);

  EigSolverData = NULL;
  done = xfe_False;
  while (!done){
    ierr = xf_Error(xf_EigIterLanczos(NULL, nev, ncv, VS, UTOL5, xfe_VerbosityLow, 
				      0, xfe_False, xfe_False, E, EV, &V, &W, 
				      &EigSolverData, &Status));
    xf_AssertEqual(ierr, xf_OK);
    done = (Status == xfe_EigConverged);
    if (!done){
      // W = A*V
      rV = V->GenArray[0].rValue[0];
      rW = W->GenArray[0].rValue[0];
      for (k=0; k<r; k++){
	rW[k] = 2.0*rV[k];
	if (k>0  ) rW[k] += rV[k-1];
	if (k<r-1) rW[k] += rV[k+1];
      }
      xf_AssertEqual((EigSolverData->iIter < 100), xfe_True);
    }
  }

  xf_AssertRealVectorWithin(E, E0, 5, UTOL5);
  rV = EV->Vector[0].GenArray[0].rValue[0];
  xf_AssertRealVectorWithin(rV, EV0, 16, UTOL5);

 /*  for (k=0; k<nev; k++){ */
/*     xf_printf("E[%d] = %.15E\n", k, E[k]); */
/*     rV = EV->Vector[k].GenArray[0].rValue[0]; */
/*     for (i=0; i<r; i++){ */
/*       xf_printf("  %.15E\n", rV[i]); */
/*     } */
/*   } */

  ierr = xf_Error(xf_DestroyVectorSet(VS));
  xf_AssertEqual(ierr, xf_OK);

  ierr = xf_Error(xf_DestroyVectorSet(EV));
  xf_AssertEqual(ierr, xf_OK);


  return xf_OK;  
}
