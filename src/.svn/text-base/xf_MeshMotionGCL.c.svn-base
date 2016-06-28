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
  FILE:  xf_MeshMotionGCL.c

  This file contains functions that implement the Geometric
  Conservation Law for mesh motion.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_Math.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xf_Quad.h"
#include "xf_Basis.h"
#include "xf_MeshTools.h"
#include "xf_MeshMotion.h"
#include "xf_SolverTools.h"
#include "xf_Solver.h"
#include "xf_LinearSolver.h"
#include "xf_Residual.h"


/******************************************************************/
//   FUNCTION Definition: xf_FindMeshMotionGCLVector
int 
xf_FindMeshMotionGCLVector( xf_All *All, xf_Vector **pGCL)
{
  int ierr;
  xf_Data *D = NULL;

  // find GCL vector by title
  ierr = xf_Error(xf_FindDataByTitle(All->DataSet, GCLVectorTitle, 
                                     xfe_Vector, &D));
  if (ierr != xf_OK) return ierr;
 
  // return pointer to GCL vector
  if (pGCL != NULL) (*pGCL) = (xf_Vector *) D->Data;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_InitMeshMotionGCLVector
int 
xf_InitMeshMotionGCLVector( xf_All *All, xf_Vector *U, int ind, 
                            enum xfe_Bool InitFlag, xf_Vector **pGCL)
{
  int ierr;
  enum xfe_Bool Found;
  char Title[xf_MAXSTRLEN];
  xf_Data *D = NULL;
  xf_Vector *GCL = NULL;

  // try finding it; if does not exist, create 
  if (ind < 0)
    sprintf(Title, "%s", GCLVectorTitle);
  else
    sprintf(Title, "%s_%d", GCLVectorTitle, ind);
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, 
                                U->Basis, U->Order, U->nComp, U->vOrder, NULL, 
				xfe_SizeReal, xfe_True, xfe_True, &D, &GCL, &Found));
  if (ierr != xf_OK) return ierr;
  // make writeable, because it needs to accompany the State
  D->ReadWrite = xfe_True; 
  // set solver role so that GCL data is kept during re-partitioning
  GCL->SolverRole = xfe_SolverRoleOther;

  // initialize to 1 if just created or if InitFlag is True
  if ( (!Found) || (InitFlag)){
    // initialize to 1 as if this were a Lagrange-interpolated vector
    ierr = xf_Error(xf_SetConstVector(GCL, 1, 1.0));
    if (ierr != xf_OK) return ierr;
    // convert from Lagrange basis to actual Basis    
    ierr = xf_Error(xf_ConvertVectorFromLagrange(GCL));
    if (ierr != xf_OK) return ierr;
  }
  
  // return GCL if asked for it (vector is already in All)
  if (pGCL != NULL) (*pGCL) = GCL;

  return xf_OK;
} 


/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionGCL_LinTitle
static void
xf_MeshMotionGCL_LinTitle( const char *Root, int ind, char *Title)
{
  /* Provides a consistent name for linearization of "Root" with
     respect to the GCL state */
  if (ind >= 0)
    sprintf(Title, "%s_GCL_%d", Root, ind);
  else
    sprintf(Title, "%s_GCL", Root);
}


/******************************************************************/
//   FUNCTION Definition: xf_FindMeshMotionGCLLinearization
int 
xf_FindMeshMotionGCLLinearization( xf_All *All, const char *Root, int ind,
				   xf_Vector **pRoot_GCL)
{
  int ierr;
  char Title[xf_MAXSTRLEN];
  xf_Data *D = NULL;

  // find desired vector by title
  xf_MeshMotionGCL_LinTitle(Root, ind, Title);
  ierr = xf_FindDataByTitle(All->DataSet, Title, xfe_Vector, &D);
  if (pRoot_GCL != NULL)
    (*pRoot_GCL) = (ierr == xf_OK) ? (xf_Vector *) D->Data : NULL;

  return ierr;
}



/******************************************************************/
//   FUNCTION Definition: xf_InitMeshMotionGCLLinearization
int 
xf_InitMeshMotionGCLLinearization( xf_All *All, xf_Vector *U, const char *OutputName,
                                   int ind, enum xfe_Bool InitFlag, xf_Vector **pJ_GCL)
{
  int ierr, len;
  char Title[xf_MAXSTRLEN];
  enum xfe_Bool Found;
  xf_Vector *J_GCL = NULL;

  // look for or create a linearization vector, J_GCL
  xf_MeshMotionGCL_LinTitle(OutputName, ind, Title);
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				U->Basis, U->Order, U->nComp, U->vOrder, NULL, 
				xfe_SizeReal, xfe_True, xfe_True, NULL, &J_GCL, &Found));
  if (ierr != xf_OK) return ierr;

  // J_GCL = 0 if necessary
  if ( (!Found) || (InitFlag)){
    ierr = xf_Error(xf_SetZeroVector(J_GCL));
    if (ierr != xf_OK) return ierr;
  }

  if (pJ_GCL != NULL) (*pJ_GCL) = J_GCL;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_InitMeshMotionGCLAdjoint
int 
xf_InitMeshMotionGCLAdjoint( xf_All *All, xf_Vector *Psi, int ind,
                             enum xfe_Bool InitFlag, xf_Vector **pPsiGCL)
{
  int ierr, len;
  char Title[xf_MAXSTRLEN];
  enum xfe_Bool Found;
  xf_Vector *PsiGCL = NULL;
  xf_Vector *J_GCL = NULL;

  // try finding it; if does not exist, create 
  if (ind < 0)
    sprintf(Title, "%s_%s", GCLAdjointTitle, Psi->OutputName);
  else
    sprintf(Title, "%s_%s_%d", GCLAdjointTitle, Psi->OutputName, ind);
  ierr = xf_Error(xf_FindVector(All, Title, xfe_LinkageGlobElem, 1, NULL, 0, 0, 
				Psi->Basis, Psi->Order, Psi->nComp, Psi->vOrder, NULL, 
				xfe_SizeReal, xfe_True, xfe_True, NULL, &PsiGCL, &Found));
  if (ierr != xf_OK) return ierr;
  
  // Set OutputName (if not found)
  if (!Found){
    len = strlen(Psi->OutputName)+1;
    ierr = xf_Error(xf_ReAlloc( (void **) &PsiGCL->OutputName, len, sizeof(char)));
    if (ierr != xf_OK) return ierr;
    strcpy(PsiGCL->OutputName, Psi->OutputName);
  }      
    
  // Initialize adjoint to zero if requested or if just created
  if ( (!Found) || (InitFlag)){
    ierr = xf_Error(xf_SetZeroVector(PsiGCL));
    if (ierr != xf_OK) return ierr;
  }

  // return PsiGCL if asked for it (vector is already in All)
  if (pPsiGCL != NULL) (*pPsiGCL) = PsiGCL;

  // also make sure that a linearization vector, J_GCL, exists for this output  
  ierr = xf_Error(xf_InitMeshMotionGCLLinearization(All, Psi, Psi->OutputName, 
                                                    ind, InitFlag, &J_GCL));
  if (ierr != xf_OK) return ierr;


  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_GCLFlux
static void
xf_GCLFlux( xf_MotionData *MD, int nq, int dim, real *F)
{      
  int iq, d, dim2;

  dim2 = dim*dim;

  // first set F{d,q} = vg{q,d}
  for (iq=0; iq<nq; iq++) 
    for (d=0; d<dim; d++)
      F[d*nq+iq] = MD->vg[dim*iq+d];
  
  // then multiply F{d,q} by g*Ginv
  for (iq=0; iq<nq; iq++) 
    xf_dMxM(MD->Ginv+iq*dim2, MD->g[iq], dim, 1, nq, F+iq);

}



/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionGCLResidual
int 
xf_MeshMotionGCLResidual( xf_All *All, xf_Vector *GCL, enum xfe_Bool ZeroFlag, xf_Vector *RGCL)
{
  int ierr, d;
  int egrp, elem, face;
  int iq, nq, nn, pnq, dim;
  int Order, QuadOrder;
  int Orient;
  enum xfe_Bool QuadChanged = xfe_False;
  enum xfe_BasisType Basis;
  real Time;
  real *xq = NULL;
  real *g = NULL, *wq = NULL, *xglob = NULL;
  real *xelem = NULL, *wn = NULL, *F = NULL;
  real *ER = NULL;
  xf_QuadData *QuadData = NULL;
  xf_JacobianData *JData = NULL;
  xf_BasisData *PhiData = NULL;
  xf_BasisData *GeomPhiData = NULL;
  xf_BasisTable *PhiTable = NULL;
  xf_NormalData *NData = NULL;
  xf_MotionData *MD = NULL;
  xf_Mesh *Mesh;
  xf_EqnSet *EqnSet;

  EqnSet = All->EqnSet; // so that we use the same quad rules as state integration
  Mesh = All->Mesh;
  dim = Mesh->Dim;

  if (ZeroFlag){
    // RGCL = 0
    ierr = xf_Error(xf_SetZeroVector(RGCL));
    if (ierr != xf_OK) return ierr;
  }

  // create mesh motion data
  MD = NULL;
  ierr = xf_Error(xf_CreateMotionData(All, &MD));
  if (ierr != xf_OK) return ierr;
  // specify what we want allocated in mesh motion data
  ierr = xf_Error(xf_AllocMotionData(xfb_MD_vg | xfb_MD_g | xfb_MD_Ginv, 1, dim, MD));
  if (ierr != xf_OK) return ierr;

  //determine Time
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
  if (ierr != xf_OK) return ierr;

  // --- Element-interior contribution to RGCL ---
  pnq = -1;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // on each element, do an element-interior calculation and a face calculation

      // Determine Basis and Order from the state, GCL
      Basis = GCL->Basis[egrp];
      Order = xf_InterpOrder(GCL, egrp, elem);

      // determine required interior integration order
      ierr = xf_Error(xf_GetQuadOrderElem(Mesh, EqnSet, egrp, Order, &QuadOrder));
      if (ierr != xf_OK) return ierr;
      
      /* Pull off quad points for the element; will not recalculate if
	 Basis/Order have not changed. */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;

      // compute basis functions (and grads) if quad or basis or order changed
      ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, 
				   xfb_Phi | xfb_GPhi | xfb_gPhi, &PhiData));
      if (ierr != xf_OK) return ierr;
     
      /* Compute geometry Jacobian; if not constant, compute at quad
	 points.  Note if jacobian is constant, only one Jacobian will
	 be computed/returned. */
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ | xfb_iJ, 
				      QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;
      
      // convert reference basis grads (GPhi) to physical grads, gPhi
      ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
      if (ierr != xf_OK) return ierr;
      
      nn = PhiData->nn; // number of interpolation nodes
      
      // re-allocate data if quad points increased
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &F, nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      // obtain global coords of quad points
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
				      nq, xq, xglob));
      if (ierr != xf_OK) return ierr;

      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      // obtain mesh motion transformation map (analytical)
      ierr = xf_Error(xf_MeshMotionMap( -1, -1, NULL, All->Mesh->Motion, 
					nq, dim, Time, xglob, MD));
      if (ierr != xf_OK) return ierr;
  
      // GCL residual on elem [nn*1]
      ER = RGCL->GenArray[egrp].rValue[elem]; 
  
      // calculate F at quad points, using MD
      xf_GCLFlux(MD, nq, dim, F);
	
      // multiply F by quad weights*J
      for (d=0; d<dim; d++)
	xf_ColMult(F+nq*d, wq, nq, 1, 1); // F is modified here
      
      // Add to R: ER{n} += sum_d sum_q gPhi{d,q,n}^T * F{i,q}*wq{q}
      for (d=0; d<dim; d++)
	xf_MTxM_Add(PhiData->gPhi+nn*nq*d, F+nq*d, nn, nq, 1, ER);

    } // elem
  } // egrp

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;
  QuadData = NULL;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  PhiData = NULL;

  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  GeomPhiData = NULL;

  // Release memory
  xf_Release( (void *) xglob);  xglob = NULL;
  xf_Release( (void *) F);      F     = NULL;
  xf_Release( (void *) wq);     wq    = NULL;


  // --- Face contribution to RGCL ---

  // create a basis look-up table
  ierr = xf_Error(xf_CreateBasisTable(&PhiTable));
  if (ierr != xf_OK) return ierr;
  
  pnq = -1;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      
      // Determine Basis and Order from the state, GCL
      Basis = GCL->Basis[egrp];
      Order = xf_InterpOrder(GCL, egrp, elem);

      // GCL residual on elem [nn*1]
      ER = RGCL->GenArray[egrp].rValue[elem]; 

      // loop over faces
      for (face=0; face<Mesh->ElemGroup[egrp].nFace[elem]; face++){

	// determine required integration order
	ierr = xf_Error(xf_GetQuadOrderGeneralFace(Mesh, EqnSet, egrp, Order, &QuadOrder));
	if (ierr != xf_OK) return ierr;

	/* Pull off quad points for the iface; will not recalculate if
	   Basis/Order have not changed. */
	ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face, QuadOrder, &QuadData, &QuadChanged));
	if (ierr != xf_OK) return ierr;
    
	nq = QuadData->nquad;
	xq = QuadData->xquad;
	wq = QuadData->wquad;

	// determine face orientation (should not matter; doing one elem at a time)
	Orient = 0;
    
	// compute basis functions if quad or basis or order changed
	ierr = xf_Error(xf_EvalBasisOnFaceUsingTable(Mesh, egrp, elem, face, Orient,
						     Basis, Order, QuadChanged, nq, xq, 
						     xfb_Phi, &PhiData, PhiTable, &xelem));
	if (ierr != xf_OK) return ierr;

	/* Compute normal(s) at quad points.  If face is straight, only
	   one normal will be computed/returned. */
	ierr = xf_Error(xf_ElemNormal(Mesh, egrp, elem, face, Orient, nq, xq, &NData));
	if (ierr != xf_OK) return ierr;
    
	nn = PhiData->nn;

	// re-allocate data if quad points increased
	if (nq > pnq){
	  ierr = xf_Error(xf_ReAlloc( (void **) &wn, nq*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	  ierr = xf_Error(xf_ReAlloc( (void **) &F, nq*dim, sizeof(real)));
	  if (ierr != xf_OK) return ierr;
	}

	// construct wn = weighted normals = normals multiplied by quad weights
	for (d=0; d<dim; d++)
	  for (iq=0;iq<nq; iq++) 
	    wn[iq*dim+d] = NData->n[iq*dim*(NData->nq!=1)+d]*wq[iq];

	// obtain global coords of quad points
	ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, 
					xfe_True, nq, xelem, xglob));
	if (ierr != xf_OK) return ierr;

	// obtain transformation map if doing mesh motion
	ierr = xf_Error(xf_MeshMotionMap( -1, -1, NULL, All->Mesh->Motion, 
					  nq, dim, Time, xglob, MD));
	if (ierr != xf_OK) return ierr;
	
	// calculate GCL fluxes
	xf_GCLFlux(MD, nq, dim, F);
	
	// dot F with normal (includes weights) store in first part of F
	xf_ColMult(F+0*nq, wn+0, nq, 1, dim);
	for (d=1; d<dim; d++) xf_ColMult_Add(F+d*nq, wn+d, nq, 1, dim, F+0*nq);

	// Subtract from R: ER{n} -= sum_q Phi{q,n}^T * F{0,q}
	xf_MTxM_Sub(PhiData->Phi, F, nn, nq, 1, ER);

      } // face
    } // elem
  } // egrp


  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_False));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Table */
  ierr = xf_Error(xf_DestroyBasisTable(PhiTable));
  if (ierr != xf_OK) return ierr;

  /* Destroy Geometry Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;
  
  /* Destroy Normal data */
  ierr = xf_Error(xf_DestroyNormalData(NData));
  if (ierr != xf_OK) return ierr;

  // Destroy mesh motion data
  xf_DestroyMotionData(MD);

  // Release memory
  xf_Release( (void *) xglob);
  xf_Release( (void *) xelem);
  xf_Release( (void *) F);
  xf_Release( (void *) wn);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionGCLSolveSystem
int 
xf_MeshMotionGCLSolveSystem( xf_All *All, real c, xf_Vector *S, xf_Vector *GCL)
{
  int ierr;
  
  // Compute GCL Residual vector, add to S
  ierr = xf_MeshMotionGCLResidual(All, GCL, xfe_False, S);
  if (ierr != xf_OK) return ierr;
  
  // multiply S by -1.0/c*inv(M)
  ierr = xf_Error(xf_MultInvMassMatrix(All, -1.0/c, NULL, S));
  if (ierr != xf_OK) return ierr;

  // Set GCL = S
  ierr = xf_Error(xf_SetVector(S, xfe_Set, GCL));
  if (ierr != xf_OK) return ierr;

  // communicate halo data in GCL
  ierr = xf_Error(xf_HaloExchangeVectorBegin(GCL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_HaloExchangeVectorEnd(GCL));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionMap_GCL
int 
xf_MeshMotionMap_GCL( xf_Vector *GCL, int egrp, int elem, xf_BasisData *PhiData, 
		      int npoint, int dim, real *gb, real *gbigb_X)
{
  int ierr, d;
  int i, nn, sr;
  real val;
  real *G;

  nn = PhiData->nn;   // number of basis functions
  sr = 1;             // GCL state rank

  // Need gb ("gbar")
  if (gb == NULL) return xf_Error(xf_INPUT_ERROR);

  // GCL on elem egrp, elem [nn*sr]
  G = GCL->GenArray[egrp].rValue[elem];

  /* interpolate gradient first, using gb as temp storage to account
     for non-standard storage of gbigb_X (dim is fastest running
     unrolled index) */
  if ((gbigb_X != NULL) && (PhiData->gPhi != NULL))
    for (d=0; d<dim; d++){
      xf_MxM_Set(PhiData->gPhi+nn*npoint*d, G, npoint, nn, sr, gb);
      for (i=0; i<npoint; i++) gbigb_X[dim*i + d] = gb[i];
    }

  // interpolate state at given points
  xf_MxM_Set(PhiData->Phi, G, npoint, nn, sr, gb);

  // multiply gradient by 1/gbar
  if (gbigb_X != NULL)
    for (i=0; i<npoint; i++)
      for (d=0, val=gb[i]; d<dim; d++)
        gbigb_X[dim*i+d] /= val;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DGTimeFindGCLVectors
int
xf_DGTimeFindGCLVectors(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                        xf_Vector ***pGCLj, xf_Data ***pDj)
{
  int ierr, OrderTime, j;
  char Title[xf_MAXSTRLEN];
  xf_Data *D;
 
  // Determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;
  
  // Allocate space for pointers to GCLTime vectors, which store the computed GCL at each time node
  if (pGCLj != NULL){
    ierr = xf_Error(xf_Alloc( (void **) pGCLj, OrderTime+1, sizeof(xf_Vector *) ));
    if(ierr != xf_OK) return ierr;
  }

  // Allocate space for pointers to parent data structures, if requested
  if (pDj != NULL){
    ierr = xf_Error(xf_Alloc( (void **) pDj, OrderTime+1, sizeof(xf_Data *) ));
    if (ierr != xf_OK) return ierr;
  }
  
  // Find already-computed GCLTime vectors
  for(j=0;j<=OrderTime;j++){
    sprintf(Title, "%s_%d", GCLVectorTitle, j);
    ierr = xf_Error(xf_FindDataByTitle(All->DataSet, Title, xfe_Vector, &D));
    if (ierr != xf_OK) return ierr;
    if (pGCLj != NULL) (*pGCLj)[j] = (xf_Vector *) D->Data;
    if (pDj != NULL) (*pDj)[j] = D;
  }

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DGTimeInterpolateGCL
int
xf_DGTimeInterpolateGCL(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                        real tq, real *phiout)
{
  int ierr;
  xf_Vector **GCLj;
  xf_Vector *GCL;
  
  // Find primary GCL vector stored in All->DataSet
  ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
  if (ierr != xf_OK) return ierr;
  
  // Find already-computed GCLTime vectors
  ierr = xf_Error(xf_DGTimeFindGCLVectors(All, TimeScheme, &GCLj, NULL));
  if(ierr != xf_OK) return ierr;
  
  // Interpolate GCLTime to GCL(tq) 
  ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, GCLj, 1, &tq, &GCL, phiout));
  if(ierr != xf_OK) return ierr;
    
  xf_Release( (void *) GCLj);

  return 0;
}


/******************************************************************/
//   FUNCTION Definition: xf_DGTimeFillSourceGCL
static int
xf_DGTimeFillSourceGCL(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                       xf_Vector **GCLj, xf_Vector *GCLprev, 
                       real Time, real TimeStep, xf_Vector **SGCLi)
{
 /*
 PURPOSE: 

   Compute int_{Time}^{Time+TimeStep} Phi^n(t)RGCL(x,t), where
   RGCL(x,t) is the spatial residual term in the GCL equation and
   "n" is the temporal node that Phi is associated with. E.g. for
   the first basis function on the timeslab, n=0, while for the last
   basis in a DG2 scheme, n=2.
   
   Also includes contribution of GCL contribution from previous time
   slab, GCLprev.

 INPUTS:
 
   All: All struct
   TimeScheme: Timescheme, e.g. DG1 or DG2
   GCLj : GCL vectors on current time slab
   GCLprev : GCL at end of previous time slab
   Time : time at start of slab
   TimeStep: current timestep
   
 OUTPUTS:
 
   SGCLi: Filled source at each temporal node
   
 RETURNS:
   
   Error Code
 */
    
  int ierr, OrderTime, iOrder, it, nq, n;
  char Title[xf_MAXSTRLEN];
  real phi[xf_MAXDGTIMENODE]; //stores temporal bases evaluated at a given quad point
  real tq[xf_MAXDGTIMENODE]; //stores times at quad points
  real wq[xf_MAXDGTIMENODE]; //stores temporal quad weights
  real Vprev[xf_MAXDGTIMENODE];
  xf_Vector *GCL;
  xf_Vector *RGCL; //temporary vector storing GCL residual at each temporal quad point

  
  // Get temporal order from input TimeScheme
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;
  
  // Get variables for TimeScheme
  ierr = xf_Error(xf_DGTimeSchemeVars(TimeScheme, NULL, Vprev, &nq, tq, wq, NULL));
  if (ierr != xf_OK) return ierr;


  // allocate a temporary residual
  sprintf(Title, "ResidualGCL");
  ierr = xf_Error(xf_FindSimilarVector(All, GCLj[0], "ResidualGCL", xfe_True, xfe_True, NULL, &RGCL, NULL));
  if(ierr != xf_OK) return ierr;

  // initialize SGCLi to zero
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    ierr = xf_Error(xf_SetZeroVector(SGCLi[iOrder]));
    if (ierr != xf_OK) return ierr;
  }

  // Find primary GCL vector stored in All->DataSet
  ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
  if (ierr != xf_OK) return ierr;

  //Loop over all quad points to perform integral of RGCL in time and build up source terms
  for(it=0;it<nq;it++){
  
    // Time must be modified here since xf_MeshMotionGCLResidual does not take Time as input
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time + tq[it]*TimeStep));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DGTimeInterpolateGCL(All, TimeScheme, tq[it], phi));
    if (ierr != xf_OK) return ierr;

    // Compute RGCL over whole domain at current temporal quad point
    xf_MeshMotionGCLResidual(All, GCL, xfe_True, RGCL);
    
    // Fill Source vectors associated with each temporal basis
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      // Compute Source[iOrder] = int_{Time}^{Time + TimeStep} Phi^iOrder * R(x,t) 
      ierr = xf_Error(xf_VectorMultSet(RGCL, phi[iOrder]*wq[it]*TimeStep, xfe_Add, SGCLi[iOrder]));
      if (ierr != xf_OK) return ierr; 
    }
  } // end loop over quad points
  
  // add contribution from GCLprev
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    if (Vprev[iOrder] != 0.){
      // add Vprev*M'*Uprev to source
      ierr = xf_Error(xf_AddMassMatrix(All, Vprev[iOrder], NULL, GCLprev, SGCLi[iOrder], NULL, NULL));
      if (ierr != xf_OK) return ierr;
    }
  }

  // Set Time back to beginning of time step (just in case)
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DGTimeUnsteadyResidualGCL
int
xf_DGTimeUnsteadyResidualGCL(xf_All *All, enum xfe_TimeSchemeType TimeScheme, 
                             xf_Vector **GCLj, xf_Vector *GCLprev, xf_Vector *GCLtemp,
                             real Time, real TimeStep, xf_Vector **RGCLi)
{    
  int ierr;
  int iOrder, jOrder, OrderTime;
  real Atime[xf_MAXDGTIMENODE*xf_MAXDGTIMENODE];

  // Get temporal order from input TimeScheme
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;
  
  // Get variables for TimeScheme
  ierr = xf_Error(xf_DGTimeSchemeVars(TimeScheme, Atime, NULL, NULL, NULL, NULL, NULL));
  if (ierr != xf_OK) return ierr;

  // initialize RGCLi to zero
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    ierr = xf_Error(xf_SetZeroVector(RGCLi[iOrder]));
    if (ierr != xf_OK) return ierr;
  }
  
  // Add GCL "source vector" to RGCL: i.e. residual without temporal stiffness matrix term
  ierr = xf_Error(xf_DGTimeFillSourceGCL(All, TimeScheme, GCLj, GCLprev, Time, TimeStep, RGCLi));
  if(ierr != xf_OK) return ierr;

  // Take care of temporal stiffness matrix
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    // temporal stiffness matrix
    for (jOrder = 0; jOrder <= OrderTime; jOrder++){
      ierr = xf_Error(xf_VectorMultSet(GCLj[jOrder],  Atime[iOrder*(OrderTime+1)+jOrder], 
				       ((jOrder==0) ? xfe_Set : xfe_Add), GCLtemp));
      if (ierr != xf_OK) return ierr;
    }
    // add M*GCLtemp to temporal residual RGCLi[iOrder]
    ierr = xf_Error(xf_AddMassMatrix(All, 1.0, NULL, GCLtemp, RGCLi[iOrder], NULL, NULL));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_DGTimeGCLSolve
int 
xf_DGTimeGCLSolve(xf_All *All, enum xfe_TimeSchemeType TimeScheme, real Time, 
                  real TimeStep, xf_SolverData *SolverData, xf_Vector *GCLprev,
                  xf_Vector **GCLj)
{

  int OrderTime, ierr, i, n, nq;
  char Title[xf_MAXSTRLEN];
  real iAT[xf_MAXDGTIMENODE*xf_MAXDGTIMENODE]; //inverse temporal stiffness matrix
  xf_Vector **SGCLi = NULL; //Source terms in GCL equation associated with each temporal basis
  
  // Determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;
  
  // Obtain number of temporal quad points
  ierr = xf_Error(xf_DGTimeSchemeVars(TimeScheme, NULL, NULL, &nq, NULL, NULL, iAT));
  if (ierr != xf_OK) return ierr;
  
  // Allocate pointers to Source GCL vector of vectors
  ierr = xf_Error(xf_Alloc( (void **) &SGCLi, OrderTime+1, sizeof(xf_Vector *)));
  if(ierr != xf_OK) return ierr;
    
  // Allocate SourceGCL for each temporal node and initialize both it and GCLTime to zero
  for(i=0;i<OrderTime+1;i++){
    sprintf(Title, "SourceGCL_%d", i);
    ierr = xf_Error(xf_FindSimilarVector(All, GCLj[0], Title, xfe_True, xfe_True, NULL, SGCLi+i, NULL));
    if(ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_SetZeroVector(SGCLi[i]));
    if(ierr != xf_OK) return ierr;
  }
  
  // Compute source vector associated with each temporal basis
  ierr = xf_Error(xf_DGTimeFillSourceGCL(All, TimeScheme, GCLj, GCLprev, Time, TimeStep, SGCLi));
  if(ierr != xf_OK) return ierr;
  
  // Multiply source vectors by negative inverse mass matrices (Source[n] = -M^{-1}Source[n])
  for(i=0;i<=OrderTime;i++){
    ierr = xf_Error(xf_MultInvMassMatrix(All, -1.0, NULL, SGCLi[i]));
    if (ierr != xf_OK) return ierr;
  }
    
  // Multiply modified source vectors by inverse temporal stiffness matrix entries to obtain GCLTime vectors
  for(n=0;n<=OrderTime;n++){
    for(i=0;i<=OrderTime;i++){
      ierr = xf_Error(xf_VectorMultSet(SGCLi[i], iAT[i+n*(OrderTime+1)],
                                       (i==0) ? xfe_Set : xfe_Add, GCLj[n]));
      if (ierr != xf_OK) return ierr;
    }
  }

  //Memory cleanup
  xf_Release( (void *) SGCLi);
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DGinTimeStepGCLAdjoint
int
xf_DGinTimeStepGCLAdjoint(xf_All *All, enum xfe_TimeSchemeType TimeScheme,
                          real Time, real TimeStep, int iSlab, int nSlab,
                          xf_SolverData *SolverData, xf_Vector **Uj,
                          int nPsi,  xf_Vector ***Psii, xf_Vector **PsiGCLnext, 
                          xf_Vector ***PsiGCLi)
{
/*
PURPOSE:

  Takes a step of a DG in time unsteady adjoint solver.  On input, the
  latest adjoint (left node of future adjacent time slab) is stored in
  Psinext.  On output, the entire adjoint in the current time slab is
  computed and stored in Psii.
  
INPUTS:

  All         : all structure
  TimeScheme  : what time scheme to use
  Time        : current time
  TimeStep    : delta t
  iSlab       : current time slab
  nSlab       : number of slabs
  SolverData  : solver data structure (contains CFL)
  Uj          : state vectors on current time slab
  nPsi        : number of adjoints
  Psi         : generic adjoint array storing OutputName
  Psinext     : adjoint on left node of next time slab
  
OUTPUTS: 

  Psii        : pointer to ALL the adjoints on the current time slab (0=left)
  Redo        : on error, this flag is set to True to redo the time step

RETURN: Error code

*/

  int ierr, nIter, iAdjoint, iq, nq;
  int iOrder, jOrder, OrderTime;
  enum xfe_Verbosity Verbosity;
  char *OutputName;
  char Title[xf_MAXSTRLEN];
  real tq[xf_MAXDGTIMENODE], wq[xf_MAXDGTIMENODE]; 
  real Vprev[xf_MAXDGTIMENODE];
  real iAT[xf_MAXDGTIMENODE*xf_MAXDGTIMENODE];
  real phi[xf_MAXDGTIMENODE];
  xf_Vector **Adji, *Utemp, *Adjtemp;;
  xf_Vector **SGCLj, **AdjGCLi, *AdjGCLnext, *J_GCL, *GCLtemp;
  xf_Vector *GCL;


  // determine order in time
  ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
  if (ierr != xf_OK) return ierr;

  // pull off variables for this time scheme
  ierr = xf_Error(xf_DGTimeSchemeVars(TimeScheme, NULL, Vprev, &nq, tq, wq, iAT));
  if (ierr != xf_OK) return ierr;

  // check input args
  if (nPsi <= 0) return xf_Error(xf_INPUT_ERROR);

  // determine verbosity level
  ierr = xf_Error(xf_GetKeyValueEnum(All->Param->KeyValue, "Verbosity", 
                                     xfe_VerbosityName, (int ) xfe_VerbosityLast, 
                                     (int *) &Verbosity));
  if (ierr != xf_OK) return ierr;
  
  // Find GCL vector
  ierr = xf_Error(xf_FindMeshMotionGCLVector(All, &GCL));
  if (ierr != xf_OK) return ierr;


  // Find source vector for GCL adjoint
  ierr = xf_Error(xf_Alloc( (void **) &SGCLj, OrderTime+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (iOrder = 0; iOrder <= OrderTime; iOrder++){
    sprintf(Title, "GCLAdjointSource_%d", iOrder);
    ierr = xf_Error(xf_FindSimilarVector(All, PsiGCLi[0][0], Title, xfe_True, 
                                         xfe_True, NULL, SGCLj + iOrder, NULL));
    if (ierr != xf_OK) return ierr;
  }

  // locate temporary state vector
  ierr = xf_Error(xf_FindSimilarVector(All, Uj[0], "Utemp_DG", xfe_True, xfe_True, 
				       NULL, &Utemp, NULL));
  if (ierr != xf_OK) return ierr;

  // locate temporary GCL vector
  ierr = xf_Error(xf_FindSimilarVector(All, GCL, "GCLtemp_DG", xfe_True, xfe_True, 
				       NULL, &GCLtemp, NULL));
  if (ierr != xf_OK) return ierr;

  // locate temporary non-GCL adjoint vector
  ierr = xf_Error(xf_FindSimilarVector(All, Psii[0][0], "Adjtemp_DG", xfe_True, xfe_True, 
				       NULL, &Adjtemp, NULL));
  if (ierr != xf_OK) return ierr;


  /***  Loop over adjoints  ***/
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    
    // current adjoint vectors
    Adji       = Psii[iAdjoint];
    AdjGCLi    = PsiGCLi[iAdjoint]; 
    AdjGCLnext = PsiGCLnext[iAdjoint];
    
    // associated output name (make sure not NULL)
    if ((OutputName = AdjGCLi[0]->OutputName) == NULL) 
      return xf_Error(xf_INPUT_ERROR);
    
    // no need to initialize AdjGCLi, since we'll be doing a direct solve
        
    /*** Calculate RHS (on left) for adjoint equations  ***/
      
    // initialize source to zero
    for (iOrder = 0; iOrder <= OrderTime; iOrder++){
      ierr = xf_Error(xf_SetZeroVector(SGCLj[iOrder]));
      if (ierr != xf_OK) return ierr;
    }      

    /* add contributions from stiffness matrix terms and influence of Adjnext*/
    for (jOrder = 0; jOrder <= OrderTime; jOrder++){

      /* Output linearization J_GCL should already be computed from
         regular adjoint call to IncrementUnsteadyOutputs.  */
      ierr = xf_Error(xf_FindMeshMotionGCLLinearization(All, OutputName, jOrder, &J_GCL));
      if (ierr != xf_OK) return ierr;

      // set SGCLj[jOrder] = J_GCL
      ierr = xf_Error(xf_SetVector(J_GCL, xfe_Set, SGCLj[jOrder]));
      if (ierr != xf_OK) return ierr;

      // influence of AdjGCLnext
      if ((iSlab != nSlab-1) && (Vprev[OrderTime-jOrder] != 0.)){
        // add M'*AdjGCLnext to SGCLj[jOrder]
        ierr = xf_Error(xf_AddMassMatrix(All, Vprev[OrderTime-jOrder], NULL, AdjGCLnext,
                                         SGCLj[jOrder], NULL, NULL));
        if (ierr != xf_OK) return ierr;
      }
    } // jOrder

    /* add contributions from quadrature points */
    for (iq=0; iq<nq; iq++){ // loop over the quadrature points
      
      // Set Time to that at the quadrature point, Time + tq[iq]*TimeStep
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time", Time+tq[iq]*TimeStep));
      if (ierr != xf_OK) return ierr;
	
      // Calculate U at the quad point Utemp = U(tq), and phii(tq); GCL is automatically interpolated as well
      ierr = xf_Error(xf_DGTimeInterpolateState(All,TimeScheme, Uj, tq+iq, &Utemp, phi));
      if (ierr != xf_OK) return ierr;
	
      // Calculate Adjoint at the quad point, Adjtemp = Adj(tq)
      ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Adji, 1, tq+iq, &Adjtemp, NULL));
      if (ierr != xf_OK) return ierr;
      
      // GCLtemp = 0
      ierr = xf_Error(xf_SetZeroVector(GCLtemp));
      if (ierr != xf_OK) return ierr;

      // Calculate Residual_gbar contribution to GCL adjoint source -> GCLtemp
      xf_printf("Calculating PsiTxR_G using finite differences at time quad point %d out of %d\n", iq, nq);
      xf_printf(" (this may take a while)\n");
      ierr = xf_Error(xf_MeshMotionGCL_PsiTxR_G(All, Utemp, 1, &Adjtemp, &GCLtemp));
      if (ierr != xf_OK) return ierr;
      xf_printf("Done calculating PsiTxR_G\n");
      
      // SGCLj[j] += wq(iq)*phij(tq)*TimeStep*GCLtemp
      for (jOrder = 0; jOrder <= OrderTime; jOrder++){
        ierr = xf_Error(xf_VectorMultSet(GCLtemp, wq[iq]*phi[jOrder]*TimeStep,
                                         xfe_Add, SGCLj[jOrder]));
        if (ierr != xf_OK) return ierr;
      }
    } // iq

    
    /* Adjoint solve = transpose mass matrix inversion */
    
    // Multiply source vectors by inverse mass matrices
    for (jOrder=0; jOrder<=OrderTime; jOrder++){
      ierr = xf_Error(xf_MultInvMassMatrix(All, 1.0, NULL, SGCLj[jOrder]));
      if (ierr != xf_OK) return ierr;
    }

    /* Multiply modified source vectors by neg. inverse temporal
       stiffness matrix transpose (-iAT)^T. Store result in AdjGCLi*/
    for (iOrder=0; iOrder<=OrderTime; iOrder++){
      for (jOrder=0; jOrder<=OrderTime; jOrder++){
        ierr = xf_Error(xf_VectorMultSet(SGCLj[jOrder], iAT[iOrder+jOrder*(OrderTime+1)], 
                                         ((jOrder==0) ? xfe_Neg : xfe_Sub), AdjGCLi[iOrder]));
        if (ierr != xf_OK) return ierr;
      }
    }
 
  } // iAdjoint

  // memory cleanup
  xf_Release( (void *) SGCLj);

  return xf_OK;

}





/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionGCLSolveAdjoints
int
xf_MeshMotionGCLSolveAdjoints(xf_All *All, real c, real d, int nPsi, 
			      xf_Vector **SGCL, xf_Vector **PsiGCL)
{
  int ierr, iAdjoint;
  xf_Vector *AdjGCL, *J_GCL;
  
  // check input args
  if (nPsi <= 0) return xf_Error(xf_INPUT_ERROR);

  // need c nonzero!
  if (c == 0.) return xf_Error(xf_INPUT_ERROR);
  
  // Loop over adjoints
  for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){
    
    // current adjoint vector
    AdjGCL = PsiGCL[iAdjoint]; 
    
    // need an associated output name
    if (AdjGCL->OutputName == NULL) return xf_Error(xf_INPUT_ERROR);
    
    /* Output linearization J_GCL should already be computed from
       regular adjoint call to CalculateOutput.  */
    ierr = xf_Error(xf_FindMeshMotionGCLLinearization(All, AdjGCL->OutputName, -1, &J_GCL));
    if (ierr != xf_OK) return ierr;

    // AdjGCL = d*J_GCL
    ierr = xf_Error(xf_VectorMultSet(J_GCL, d, xfe_Set, AdjGCL));
    if (ierr != xf_OK) return ierr;
    
    // add source to AdjGCL
    if (SGCL != NULL){
      ierr = xf_Error(xf_SetVector(SGCL[iAdjoint], xfe_Add, AdjGCL));
      if (ierr != xf_OK) return ierr;
    }
    
    /* Solve linear system: RGCL_GCL^T*AdjGCL + J_GCL = 0, where
       RGCL_GCL is just c*M, since the spatial Jacobian of the GCL
       equation is zero.  Note that we are also storing J_GCL in
       AdjGCL.

       Therefore, AdjGCL = (-1/c)*M^{-1} * AdjGCL
    */
    ierr = xf_Error(xf_MultInvMassMatrix(All, (-1./c), NULL, AdjGCL));
    if (ierr != xf_OK) return ierr;
    
  } // iAdjoint
  
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionGCL_PsiTxR_GOld
int
xf_MeshMotionGCL_PsiTxR_GOld(xf_All *All, xf_Vector *U, xf_Vector *GCL, 
			  int nPsi, xf_Vector **Psi, xf_Vector **SGCL)
{
  // Calculates Psi^T * (Residual_gbar) using finite difference approximations
  // This is a slow version because it does not save "static data" for residual eval
  int ierr, k, sr;
  int egrp, elem, face, nface;
  int rmax, nn, r;
  int nN, nnN;
  int egN, eN;
  int iAdjoint;
  enum xfe_Bool Found;
  real ep = 1e-5; // epsilon for finite-differencing gbar
  real *ER0, *ER;
  real *EGCLN, *EPsi, *ESGCLN;
  xf_JacobianMatrix *R_U;
  xf_SolverData *SolverData;
  real dp;
  xf_Mesh *Mesh;

  Mesh = All->Mesh;

  // locate Jacobian matrix
  ierr = xf_Error(xf_FindJacobianMatrix(All->Mesh, All->DataSet, U, NULL,
                                        xfe_False, NULL, &R_U, &Found));
  if (ierr != xf_OK) return ierr;
  if (!Found) return xf_Error(xf_NOT_FOUND);

  // state rank
  sr = R_U->StateRank;
  
  // get rmax = max(nn*sr)
  rmax=0;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++)
      rmax = max(  sr*xf_Jacobian_n(R_U,egrp,elem), rmax);

  // allocate ER and ER0 = residual vectors on a single element
  ierr = xf_Error(xf_Alloc( (void **) &ER0, 2*rmax, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ER = ER0 + rmax;

  // create/allocate SolverData
  ierr = xf_Error(xf_CreateSolverData(All, &SolverData));
  if (ierr != xf_OK) return ierr;

  // loop over element groups
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
   
    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

      nn = ((R_U->vnvec == NULL) ? R_U->nvec[egrp] : R_U->vnvec[egrp][elem]);
      r  = nn*sr;
      nface = Mesh->ElemGroup[egrp].nFace[elem];

      // calculate ER0 = spatial residual on elem
      ierr = xf_Error(xf_CalculateResidualLeanElem(All, egrp, elem, U, ER0, NULL, 
                                                   NULL, NULL, R_U, SolverData));
      if (ierr != xf_OK) return ierr;

      // loop over self and neighbors
      for (face=-1; face<nface; face++){
        if (face==-1){ // self
          egN = egrp;
          eN  = elem;
        }
        else{          // neighbor
          egN = R_U->egrpN[egrp][elem][face];
          eN  = R_U->elemN[egrp][elem][face];
        }
       
        if (egN < 0) continue; // boundary face
       
        // # of basis functions in adjacent element
        nnN = ((R_U->vnvec == NULL) ? R_U->nvec[egN] : R_U->vnvec[egN][eN]);

        // GCL on neighbor
        EGCLN = GCL->GenArray[egN].rValue[eN];
       		
        // loop over basis functions on adjacent element
        for (nN = 0; nN < nnN; nN++){
         
          // increment GCL on adjacent element
          EGCLN[nN] += ep;

          // calculate ER = spatial residual on elem
          ierr = xf_Error(xf_CalculateResidualLeanElem(All, egrp, elem, U, ER, NULL, 
                                                       NULL, NULL, R_U, SolverData));
          if (ierr != xf_OK) return ierr;
         
          // Set ER{k} = (ER{k}-ER0{k})/ep
          for (k=0; k<r; k++) ER[k] = (ER[k]-ER0[k])/ep;

          // Add ER{k} dot Psi[iAdjoint]{k} to ESGCLN{nN}
          for (iAdjoint=0; iAdjoint<nPsi; iAdjoint++){

            // Adjoint, Psi, on current element
            EPsi = Psi[iAdjoint]->GenArray[egrp].rValue[elem];

            // SGCL on neighbor element
            ESGCLN = SGCL[iAdjoint]->GenArray[egN].rValue[eN];

            // take dot product
            for (k=0, dp=0.; k<r; k++) dp += EPsi[k]*ER[k];

            // increment SGCL on neighbor
            ESGCLN[nN] += dp;

          } // iAdjoint

          // decrement GCL on adjacent element (return to orig value)
          EGCLN[nN] -= ep;
         
        } // nN

      } // face

    } // elem

  } // egrp


  // destroy SolverData
  ierr = xf_Error(xf_DestroySolverData(SolverData));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) ER0);

  return xf_OK;

}






