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
  FILE:  xf_UnitRun.c

  This file contains functions for setting up and running cases from
  within unit tests.

*/

#include "xf_Output.h"
#include "xf_EqnSet.h"
#include "xf_EqnSetHook.h"

/******************************************************************/
//   FUNCTION Definition: xf_InitializeTestRun
static int
xf_InitializeTestRun(xf_All *All, enum xfe_BasisType Basis, int Order,
		     xf_Vector **pU)
{
/*
PURPOSE: 

  Load dynamic library, register eqnset, initialize solution

INPUTS:

  All : All structure that already exists
  Basis, Order : desired solution interpolation
                 Order < 0 means use variable order up to max |Order|
  
OUTPUTS:

  (*pU) : Primal state vector

RETURN:

  None
*/

  int ierr, i, j, k;
  xf_Vector *VOrder = NULL;
  
  // Set dimension in EqnSet
  All->EqnSet->Dim = All->Mesh->Dim;

  // Load Dynamic Library
  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;

  // Register EqnSet
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;

  // Initialize a solution vector
  ierr = xf_Error(xf_NewGREVector(All, Basis, (Order<0) ? -Order: Order, All->EqnSet->StateRank,  
				  xfe_True, pU));
  if (ierr != xf_OK) return ierr;
  (*pU)->SolverRole = xfe_SolverRolePrimalState;

  // project (*pU) if Order < 0
  if (Order < 0){
    // locate a vector for specifying desired order
    ierr = xf_Error(xf_FindVector(All, "VOrder", xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				  NULL, NULL, NULL, NULL, NULL, xfe_SizeInt, xfe_False,
				  xfe_False, NULL, &VOrder, NULL));
    if (ierr != xf_OK) return ierr;
    
    for (i=0,k=0; i<VOrder->nArray; i++)
      for (j=0; j<VOrder->GenArray[i].n; j++, k++)
	VOrder->GenArray[i].iValue[j][0] = k % (1-Order);
    
    // project (*pU)
    ierr = xf_Error(xf_ProjectVectorInPlace_VOrder(All->Mesh, All->DataSet, (*pU), NULL,
						   xfe_BasisLast, VOrder));
    if (ierr != xf_OK) return ierr;

    // Destroy VOrder
    ierr = xf_Error(xf_DestroyVector(VOrder, xfe_True));
    if (ierr != xf_OK) return ierr;
  }

  
  // store vector in a new data structure with Title="State"
  ierr = xf_Error(xf_DataSetAdd(All->DataSet, "State", xfe_Vector,
				xfe_True, (void *) (*pU), NULL));
  if (ierr != xf_OK) return ierr;
  
  /* Initialize state */
  ierr = xf_Error(xf_InitState(All, (*pU)));
  if (ierr != xf_OK) return ierr;

  // Set key values for a default run
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "nIterNonlinear", "50"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "nIterLinear", "50"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "Verbosity", "Low"));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_SetKeyValue(All->Param->KeyValue, "CFL", "1e30"));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}
