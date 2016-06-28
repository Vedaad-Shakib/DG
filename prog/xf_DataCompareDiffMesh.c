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
  FILE:  xf_DataCompareDiffMesh.c

  This program compares two data solutions on different meshes (same
  domain geometry)

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_MeshTools.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Math.h"
#include "xf_Quad.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Arg.h"
#include "xf_EqnSetHook.h"
#include "xf_Solver.h"
#include "xf_Residual.h"


/******************************************************************/
//   FUNCTION Definition: xf_CalculateErrorNorm
static int 
xf_CalculateErrorNorm(xf_All *All1, xf_All *All2, xf_Vector *U1, 
		      xf_Vector *U2, xf_Vector *dU1, const char *scalar)
{
  int ierr, dim, sr, k, d, iq, nq, pnq;
  int egrp1, elem1, egrp2, elem2;
  int Order1, Order2, QuadOrder;
  int nelemtot1, counter, skipped;
  int *IParam = NULL;
  real *RParam = NULL;
  enum xfe_BasisType Basis1, Basis2;
  enum xfe_Bool QuadChanged, skipelem;
  real *EU1, *EU2, *xq, *u1, *u2, *wq, *s1, *s2;
  real *xglob, xref[3];
  real e1, e2, enorm2, enormA, snorm2;
  real integral1, integral2;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData1;
  xf_BasisData *PhiData2;
  xf_BasisData *GeomPhiData;
  xf_JacobianData *JData;
  xf_ElemSearchStruct ESS;
  xf_Mesh *Mesh1;
  xf_Mesh *Mesh2;


  Mesh1 = All1->Mesh;
  Mesh2 = All2->Mesh;
  dim  = Mesh1->Dim;
  if (dim != Mesh2->Dim) return xf_Error(xf_INPUT_ERROR);

  // both vectors should have the same state rank
  sr   = U1->StateRank;
  if (sr != U2->StateRank) return xf_Error(xf_INPUT_ERROR);

  // number of elements in Mesh1
  ierr = xf_Error(xf_GetnElem(Mesh1, NULL, &nelemtot1));
  if (ierr != xf_OK) return ierr;

  // build element search structure on All2
  ierr = xf_Error(xf_BuildElemSearchStructure(All2, &ESS));
  if (ierr != xf_OK) return ierr;

  // pull off IParam and RParam
  if (xf_NotNull(scalar)){
    ierr = xf_Error(xf_RetrieveFcnParams(All1, All1->EqnSet, &IParam, &RParam, NULL, NULL));
    if (ierr != xf_OK) return ierr;
  }
  
  if (dU1 != NULL) return xf_Error(xf_NOT_SUPPORTED); // for now

  QuadData = NULL;
  PhiData1 = NULL;
  PhiData2 = NULL;
  JData    = NULL;
  u1       = NULL;
  u2       = NULL;
  wq       = NULL;
  xglob    = NULL;
  s1       = NULL;
  s2       = NULL;
  GeomPhiData = NULL;
  snorm2   = 0;
  enorm2   = 0;
  enormA   = 0;
  e1 = e2  = 0;
  integral1 = integral2 = 0.;
  pnq = -1;
  counter = 0;
  skipped = 0;
  for (egrp1=0; egrp1<Mesh1->nElemGroup; egrp1++){

    Basis1 = U1->Basis[egrp1]; Order1 = U1->Order[egrp1];

    for (elem1=0; elem1<Mesh1->ElemGroup[egrp1].nElem; elem1++){
      
      // print progress
      if (++counter%100 == 0)
	xf_printf(" %d/%d elems processed; skipped %d so far\n", 
		  counter, nelemtot1, skipped);

      /* Pull off quad points for the element; will not recalculate in generic case */
      QuadOrder = 2*Order1; // assume Order2 will be coarser or same order
      ierr = xf_Error(xf_QuadElem(Mesh1, egrp1, elem1, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;

      // compute basis functions for U1 vector
      ierr = xf_Error(xf_EvalBasis(Basis1, Order1, QuadChanged, 
				   nq, xq, xfb_Phi, &PhiData1));
      if (ierr != xf_OK) return ierr;

      // element Jacobian
      ierr = xf_Error(xf_ElemJacobian(Mesh1, egrp1, elem1, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;

      // re-allocate memory if quad points increased
      if (nq > pnq){
	pnq = nq;
	ierr = xf_Error(xf_ReAlloc( (void **)  &u1, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **)  &u2, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	if (xf_NotNull(scalar)){
	  ierr = xf_Error(xf_ReAlloc( (void **)  &s1, nq*1, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **)  &s2, nq*1, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}
      }

      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      // interpolate u1
      EU1 = U1->GenArray[egrp1].rValue[elem1];
      xf_MxM_Set(PhiData1->Phi, EU1, nq, PhiData1->nn, sr, u1);

      // obtain global coords of quad points
      ierr = xf_Error(xf_Ref2GlobElem(Mesh1, egrp1, elem1, &GeomPhiData, QuadChanged, 
				      nq, xq, xglob));
      if (ierr != xf_OK) return ierr;
      
      // calculate U2 at the global points -> u2
      skipelem = xfe_False;
      for (iq=0; iq<nq; iq++){
	// find element that contains xglob + dim*iq
	ierr = xf_FindElemUsingSearchStructure(All2, xglob+dim*iq, &ESS, &egrp2, 
					       &elem2, xref);
	if (ierr == xf_NOT_FOUND) skipelem = xfe_True;
	else if (ierr != xf_OK) return ierr;
	skipelem = (skipelem || (egrp2 < 0) || (elem2 < 0));
	if (skipelem) break;

	Basis2 = U2->Basis[egrp2]; Order2 = U2->Order[egrp2];
	EU2 = U2->GenArray[egrp2].rValue[elem2];

	// compute basis functions for U2 vector
	ierr = xf_Error(xf_EvalBasis(Basis2, Order2, xfe_True, 
				     1, xref, xfb_Phi, &PhiData2));
	if (ierr != xf_OK) return ierr;
	
	// interpolate u2
	xf_MxM_Set(PhiData2->Phi, EU2, 1, PhiData2->nn, sr, u2 + sr*iq);
      } // iq

      // should we skip this element?
      if (skipelem){
	xf_printf("Warning, Skipping egrp1=%d, elem1=%d.\n", egrp1, elem1);
	skipped++;
	continue;
      }
 
      // calculate scalar if requested
      if (xf_NotNull(scalar)){
	ierr = xf_Error(xf_EqnSetScalar(All1->EqnSet, scalar, IParam, RParam, nq, u1, 
					NULL, s1, NULL, NULL, NULL));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_EqnSetScalar(All1->EqnSet, scalar, IParam, RParam, nq, u2, 
					NULL, s2, NULL, NULL, NULL));
	if (ierr != xf_OK) return ierr;
	for (iq=0; iq<nq; iq++)
	  snorm2 += (s1[iq]-s2[iq])*(s1[iq]-s2[iq])*wq[iq];
      }

      // add to square error norm via quadrature
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  enorm2 += (u1[iq*sr+k]-u2[iq*sr+k])*(u1[iq*sr+k]-u2[iq*sr+k])*wq[iq];
   
	
      // add to integrals of data1 and data2
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  integral1 += (u1[iq*sr+k])*wq[iq];
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  integral2 += (u2[iq*sr+k])*wq[iq];
	  

      // add to square norms of data1 and data2
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  e1 += (u1[iq*sr+k])*(u1[iq*sr+k])*wq[iq];
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  e2 += (u2[iq*sr+k])*(u2[iq*sr+k])*wq[iq];

    } // elem1
  } // egrp1
  
  xf_printf("Integral of data1 = %.15E\n", integral1);
  xf_printf("Integral of data2 = %.15E\n", integral2);
  xf_printf("Norm of data1 = %.15E\n", sqrt(e1));
  xf_printf("Norm of data2 = %.15E\n", sqrt(e2));
  if (xf_NotNull(scalar)) 
    xf_printf("%s L2 error norm = %.15E\n", scalar, sqrt(snorm2));
  xf_printf("Continuous L2 error norm = %.15E\n", sqrt(enorm2));


  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

  // Destroy Basis Data
  ierr = xf_Error(xf_DestroyBasisData(PhiData1, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyBasisData(PhiData2, xfe_True));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  // destroy element search structure
  ierr = xf_Error(xf_DestroyElemSearchStructure(&ESS));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) u1);
  xf_Release( (void *) u2);
  xf_Release( (void *) wq);
  xf_Release( (void *) xglob);
  xf_Release( (void *) s1);
  xf_Release( (void *) s2);
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i;
  enum xfe_Bool TakeDiff, found;
  char *ArgIn[] = {"xfa1", "NULL", "first (fine) .xfa file",
		   "title1", "State", "title of first data",
		   "xfa2", "NULL", "second .xfa file",
		   "title2", "State", "title of second data",
		   "diff", "False", "Write out difference to diff.data (fine mesh)",
		   "scalar", "None", "Name of scalar for L2 error computation",
		   "\0"};
  char xfa1[xf_MAXSTRLEN];
  char xfa2[xf_MAXSTRLEN];
  char title1[xf_MAXSTRLEN];
  char title2[xf_MAXSTRLEN];
  char scalar[xf_MAXSTRLEN];
  xf_KeyValue KeyValue;
  xf_Vector *U1, *U2, *dU1;
  xf_Data *D;
  xf_All *All1, *All2;
  
  xf_printf("\n");
  xf_printf("=== Data Comparison (different meshes) ===\n");
  xf_printf("\n");
    
      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  // Get xfa1
  ierr = xf_GetKeyValue(KeyValue, "xfa1", xfa1);
  if (ierr != xf_OK) return ierr;

  // Get xfa2
  ierr = xf_GetKeyValue(KeyValue, "xfa2", xfa2);
  if (ierr != xf_OK) return ierr;

  // Get title1
  ierr = xf_GetKeyValue(KeyValue, "title1", title1);
  if (ierr != xf_OK) return ierr;

  // Get title2
  ierr = xf_GetKeyValue(KeyValue, "title2", title2);
  if (ierr != xf_OK) return ierr;

  /* TakeDiff? */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValue, "diff", &TakeDiff));
  if (ierr != xf_OK) return ierr;

  // Get scalar
  ierr = xf_GetKeyValue(KeyValue, "scalar", scalar);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  /* Create .xfa structures */
  ierr = xf_Error(xf_CreateAll(&All1, xfe_False));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_CreateAll(&All2, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Read .xfa files */
  xf_printf("Reading %s\n", xfa1);
  ierr = xf_Error(xf_ReadAllBinary(xfa1, All1));
  if (ierr!=xf_OK) return ierr;
  xf_printf("Reading %s\n", xfa2);
  ierr = xf_Error(xf_ReadAllBinary(xfa2, All2));
  if (ierr!=xf_OK) return ierr;


  if (xf_NotNull(scalar)){
    // load dynamic library
    ierr = xf_LoadEqnSetLibrary(All1->EqnSet->EqnSetLibrary);
    if (ierr != xf_OK) return ierr;
    // register equation set
    ierr = xf_Error(xf_EqnSetRegister(All1->EqnSet));
    if (ierr != xf_OK) return ierr;
  }


  /* Locate title1 data */
  found = xfe_False;
  D = All1->DataSet->Head;
  while (D != NULL){
    if ((strcmp(D->Title, title1) == 0)){
      found = xfe_True;
      break;
    }
    D = D->Next;
  }
  if (!found){
    xf_printf("Could not find title = %s in Dataset of %s.\n", title1, xfa1);
    return xf_Error(xf_NOT_FOUND);
  }
  if (D->Type != xfe_Vector) return xf_Error(xf_INPUT_ERROR);
  U1 = (xf_Vector *) D->Data;

  /* Locate title2 data */
  found = xfe_False;
  D = All2->DataSet->Head;
  while (D != NULL){
    if ((strcmp(D->Title, title2) == 0)){
      found = xfe_True;
      break;
    }
    D = D->Next;
  }
  if (!found){
    xf_printf("Could not find title = %s in Dataset of %s.\n", title2, xfa2);
    return xf_Error(xf_NOT_FOUND);
  }
  if (D->Type != xfe_Vector) return xf_Error(xf_INPUT_ERROR);
  U2 = (xf_Vector *) D->Data;

  /* Create a similar vector, dU1, if taking difference */
  if (TakeDiff){
    // locate temporary state vector, dU
    ierr = xf_Error(xf_FindSimilarVector(All1, U1, "dU", xfe_True, xfe_True,
					 NULL, &dU1, NULL));
    if (ierr != xf_OK) return ierr;

  }

  /* Compute (and print out) L2 error norm */
  ierr = xf_Error(xf_CalculateErrorNorm(All1, All2, U1, U2, ((TakeDiff) ? dU1 : NULL), 
					scalar));
  if (ierr != xf_OK) return ierr;


  // take a difference if desired
  if (TakeDiff){
    ierr = xf_Error(xf_DumpVectorBinary(All1->Mesh, "Vector", dU1, "diff.data"));
    if (ierr != xf_OK) return ierr;
  }
    
  /* Destroy .xfa structures */
  ierr = xf_Error(xf_DestroyAll(All1));
  if (ierr!=xf_OK) return ierr;
  ierr = xf_Error(xf_DestroyAll(All2));
  if (ierr!=xf_OK) return ierr;


  xf_printf("xf_DataCompareDiffMesh finished.\n");

  return xf_OK;
}
