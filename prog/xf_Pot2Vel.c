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
  FILE:  xf_Pot2Vel.c

  This program reads in a scalar potential field and writes out a
  velocity field obtained by taking the (minus) gradient of the
  potential.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_Basis.h"
#include "xf_MeshTools.h"
#include "xf_Math.h"
#include "xf_Param.h"
#include "xf_Data.h"
#include "xf_Arg.h"

/******************************************************************/
//   FUNCTION Definition: xf_Potential2Velocity
static int 
xf_Potential2Velocity(xf_All *All, xf_Vector *U)
{
  int ierr, dim, d, i, nn, nnU;
  int egrp, elem, negrp;
  int *Order, *OrderV;
  enum xfe_BasisType *Basis, *BasisV;
  char *VelocityName[] = {"XVelocity", "YVelocity", "ZVelocity"};
  real *EU, *EV, *xn, *gu;
  xf_BasisData *PhiData;
  xf_JacobianData *JData;
  xf_Vector *V;
  xf_Data *D;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;
  dim  = Mesh->Dim;
  negrp = Mesh->nElemGroup;
  
  if (U->StateRank != 1) return xf_Error(xf_INPUT_ERROR);
  if ((Basis = U->Basis) == NULL) return xf_Error(xf_INPUT_ERROR);
  if ((Order = U->Order) == NULL) return xf_Error(xf_INPUT_ERROR);
  
  // BasisV = basis used for velocity interpolation
  ierr = xf_Error(xf_Alloc( (void **) &BasisV, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (egrp=0; egrp<negrp; egrp++){
    ierr = xf_Error(xf_Basis2Lagrange(Basis[egrp], BasisV+egrp));
    if (ierr != xf_OK) return ierr;
  }

  // OrderV = order used for velocity interpolation
  ierr = xf_Error(xf_Alloc( (void **) &OrderV, negrp, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (egrp=0; egrp<negrp; egrp++){
    if (Order[egrp] == 0)
      xf_printf("Warning, p=0 potential will give a zero velocity.\n");
    OrderV[egrp] = max(Order[egrp]-1,0);
  }

  // Create a velocity vector
  ierr = xf_Error(xf_FindVector(All, "VelocityField", xfe_LinkageGlobElem, dim, VelocityName, 
				0, 0, Basis, OrderV, NULL, NULL, NULL, xfe_SizeReal, xfe_False, 
				xfe_True, &D, &V, NULL));
  if (ierr != xf_OK) return ierr;
  D->ReadWrite = xfe_True; // make velocity writable

  // no longer need BasisV or OrderV
  xf_Release( (void *) BasisV);
  xf_Release( (void *) OrderV);

/*   // set StateName for velocity */
/*   ierr = xf_Error(xf_Alloc2((void ***)&V->StateName, dim, xf_MAXSTRLEN, sizeof(char))); */
/*   if (ierr != xf_OK) return ierr; */
/*   for (d=0; d<dim; d++) */
/*     strcpy(V->StateName[d], VelocityName[d]); */


  // initialize variables to NULL
  xn          = NULL;
  gu          = NULL;
  PhiData     = NULL;
  JData       = NULL;
  
  // Loop over elements, take gradient of potential
  for (egrp=0; egrp<negrp; egrp++){

    // number of nodes for velocity interpolation
    ierr = xf_Error(xf_Order2nNode(V->Basis[egrp], V->Order[egrp], &nn));
    if (ierr != xf_OK) return ierr;

    // determine lagrange nodes for velocity interpolation
    ierr = xf_Error(xf_LagrangeNodes(V->Basis[egrp], V->Order[egrp], NULL, NULL, &xn));
    if (ierr != xf_OK) return ierr;

    // re-allocate gu
    ierr = xf_Error(xf_ReAlloc( (void **) &gu, nn*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      // compute basis functions for U vector
      ierr = xf_Error(xf_EvalBasis(U->Basis[egrp], U->Order[egrp], xfe_True, 
				   nn, xn, xfb_All, &PhiData));
      if (ierr != xf_OK) return ierr;

      nnU = PhiData->nn;

      /* Compute geometry Jacobian; if not constant, compute at quad
         points.  Note if jacobian is constant, only one Jacobian will
         be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nn, xn, xfb_iJ, 
				      xfe_True, &JData));
      if (ierr != xf_OK) return ierr;

      // convert reference basis grads (GPhi) to physical grads, gPhi
      ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
      if (ierr != xf_OK) return ierr;

      // interpolate gradient of U
      EU = U->GenArray[egrp].rValue[elem]; // U on elem [nnU*1]
      for (d=0; d<dim; d++)
	xf_MxM_Set(PhiData->gPhi+nnU*nn*d, EU, nn, nnU, 1, gu+nn*d);

      // copy gu to EV (minus sign included here since V = -grad(Potential)
      EV = V->GenArray[egrp].rValue[elem];
      for (i=0; i<nn; i++)
	for (d=0; d<dim; d++)
	  EV[i*dim+d] = -gu[d*nn+i];

    } // elem
  } // egrp

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  xf_Release( (void *) xn);
  xf_Release( (void *) gu);

  return xf_OK;

}


/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, len;
  char *ArgIn[] = {"in", "NULL", "input .xfa containing scalar potential",
		   "out", "NULL", "name for output .data with velocity",
		   "\0"};
  char inFile[xf_MAXSTRLEN];  
  char outFile[xf_MAXSTRLEN];
  xf_KeyValue KeyValue;
  xf_Data *D;
  xf_Vector *U;
  xf_All *All;
  
  xf_printf("\n");
  xf_printf("=== xf_Pot2Vel: Potential Field To Velocity  ===\n");
  xf_printf("\n");
    
      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);
    
  xf_printf("nKey = %d\n", KeyValue.nKey);
  for (i=0; i<KeyValue.nKey; i++)
    xf_printf("%d : Key = %s, Value = %s\n", i, KeyValue.Key[i], KeyValue.Value[i]);
    
  // Get inFile
  ierr = xf_GetKeyValue(KeyValue, "in", inFile);
  if (ierr != xf_OK) return ierr;

  // Get outFile
  ierr = xf_GetKeyValue(KeyValue, "out", outFile);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Read .xfa file*/
  ierr = xf_Error(xf_ReadAllBinary(inFile, All));
  if (ierr!=xf_OK) return ierr;

  
  // Pull off primal state
  ierr = xf_Error(xf_FindPrimalState(All->DataSet, 0, &D, NULL));
  if (ierr != xf_OK) return ierr;
  U = (xf_Vector *) D->Data;

  // Make all data in All->DataSet not writable
  D = All->DataSet->Head;
  while (D != NULL){
    D->ReadWrite = xfe_False;
    D = D->Next;
  }

  // Determine velocity (made writeable in All->DataSet)
  ierr = xf_Error(xf_Potential2Velocity(All, U));
  if (ierr!=xf_OK) return ierr;


  // Write velocity field (All->DataSet as .data)
  ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, All->DataSet, NULL, outFile));
  if (ierr!=xf_OK) return ierr;
  

  xf_printf("xf_Pot2Vel finished.\n");

  return xf_OK;
}
