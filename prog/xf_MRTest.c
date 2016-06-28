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
  FILE:  xf_Test.c

  This program is used for testing of model reduction

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_Memory.h"
#include "xf_All.h"
#include "xf_Mesh.h"
#include "xf_MeshTools.h"
#include "xf_Param.h"
#include "xf_Basis.h"
#include "xf_EqnSetHook.h"
#include "xf_Data.h"
#include "xf_EqnSet.h"
#include "xf_Residual.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Quad.h"
#include "xf_Arg.h"
#include "xf_MRStruct.h"
#include "xf_MRCommon.h"


// for storing the run list
typedef struct
{
  int iterm; // residual term to set
  char resKey[xf_MAXSTRLEN];   // residual key to set
  char resValue[xf_MAXSTRLEN]; // residual value to set
}
xf_RunListType;

/******************************************************************/
//   FUNCTION Definition: xf+ComputeNonlinearError
static int 
xf_ComputeNonlinearError(xf_All *All, xf_Vector *U, xf_ReducedModel *RM,
			 xf_VectorSet *CardSet, real *peM)
{
  int ierr, M, m, sr, dim;
  int egrp, elem, nq, pnq, iq;
  int *QuadOrder, *IParam;
  enum xfe_Bool QuadChanged;
  real snorm, enorm, *xref, *RParam, u0[1];
  real *sz, *u, *s, *sM, *wq, *xq;
  real *EU, *EPsi;
  xf_IPointType *z;
  xf_Mesh *Mesh;
  xf_QuadData *QuadData;
  xf_JacobianData *JData;
  xf_BasisData *PhiData;
  xf_EqnSet *EqnSet;

  Mesh = All->Mesh;
  dim = Mesh->Dim;
  EqnSet = All->EqnSet;
  sr = EqnSet->StateRank;

  if (sr != 1) return xf_Error(xf_INPUT_ERROR);
  
  // build eqnset-desired parameter lists for passing into functions
  ierr = xf_Error(xf_RetrieveFcnParams(NULL, EqnSet, &IParam, &RParam, NULL, NULL));
  if (ierr != xf_OK) return ierr;


  snorm = 0.0; // L2 norm of nonlinear term calculated from exact U
  enorm = 0.0; // L2 norm of nonlinear term calculated from the
	       // coefficient function expansion

  // pull off QuadOrder, make sure not NULL
  if ((QuadOrder = CardSet->Vector[0].Order) == NULL) 
    return xf_Error(xf_INPUT_ERROR);

  // pull of z, M, etc.
  z = RM->z;
  M = z->M;
  
  // allocate sz
  ierr = xf_Error(xf_Alloc( (void **) &sz, M, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // calculate sz at all interpolation points z
  PhiData = NULL;
  for (m=0; m<M; m++){ // loop over the M points
    egrp = z->egrp[m];
    elem = z->elem[m];
    xref = z->xref + m*dim;

    ierr = xf_Error(xf_EvalBasis(U->Basis[egrp], U->Order[egrp], xfe_True, 
				 1, xref, xfb_Phi, &PhiData));
    if (ierr != xf_OK) return ierr;

    // interpolate u at xref
    EU = U->GenArray[egrp].rValue[elem];
    xf_MxM_Set(PhiData->Phi, EU, 1, PhiData->nn, sr, u0);

    // evaluate the nonlinear term
    ierr = xf_Error(xf_EqnSetSourceS(EqnSet, RM->ResTerm, 1, IParam, RParam, 
				     1, u0, NULL, NULL, NULL, NULL, sz + m, 
				     NULL, NULL, NULL));
    if (ierr != xf_OK) return ierr;
  } // m


  // loop over element groups
  QuadData = NULL;
  JData    = NULL;
  u        = NULL;
  s        = NULL;
  sM       = NULL;
  wq       = NULL;
  pnq = -1;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    // loop over elements
    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
  
      /* Pull off quad points for the element; will not recalculate in generic case */
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder[egrp], &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;

      // compute basis functions for U vector
      ierr = xf_Error(xf_EvalBasis(U->Basis[egrp], U->Order[egrp], QuadChanged, 
				   nq, xq, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
      
      // element Jacobian
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;


      // re=allocate space for u and s
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc((void **) &u, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc((void **) &s, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &sM, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }
      
      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      

      // interpolate u at quad points
      EU = U->GenArray[egrp].rValue[elem];
      xf_MxM_Set(PhiData->Phi, EU, nq, PhiData->nn, sr, u); 

      // calculate s(U) at quad points
      ierr = xf_Error(xf_EqnSetSourceS(EqnSet, RM->ResTerm, 1, IParam, RParam, 
				       nq, u, NULL, NULL, NULL, NULL, s, NULL, 
				       NULL, NULL));
      if (ierr != xf_OK) return ierr;

      // calculate sM at quad points using CardSet and sz
      for (iq=0; iq<nq; iq++) sM[iq] = 0.;
      for (m=0; m<M; m++){
	EPsi = CardSet->Vector[m].GenArray[egrp].rValue[elem];
	if (nq != CardSet->Vector[m].GenArray[egrp].r) return xf_Error(xf_INCOMPATIBLE);
	for (iq=0; iq<nq; iq++)
	  sM[iq] += EPsi[iq]*sz[m];
      }

      // add to snorm and enorm
      for (iq=0; iq<nq; iq++) snorm += wq[iq]*s[iq]*s[iq];
      for (iq=0; iq<nq; iq++) enorm += wq[iq]*(s[iq]-sM[iq])*(s[iq]-sM[iq]);

      pnq = nq;
    } // elem

  } // egrp


  if (snorm < MEPS) return xf_Error(xf_OUT_OF_BOUNDS);

  snorm = sqrt(snorm);
  enorm = sqrt(enorm);
  //xf_printf("  enorm = %.10E, snorm = %.10E\n", enorm, snorm);

  (*peM) = enorm/snorm;

  // Only destroy QuadData if points are generic
  ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
  if (ierr != xf_OK) return ierr;

   /* Destroy geometry Jacobian Data */
  ierr = xf_Error(xf_DestroyJacobianData(JData));
  if (ierr != xf_OK) return ierr;

  /* Destroy Basis Data */
  ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
  if (ierr != xf_OK) return ierr;

  // release memory
  xf_Release( (void *) IParam);
  xf_Release( (void *) RParam);
  xf_Release( (void *) sz);
  xf_Release( (void *) sM);
  xf_Release( (void *) u);
  xf_Release( (void *) s);
  xf_Release( (void *) wq);

  return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i, iRun, nRun;
  enum xfe_Bool BatchMode;
  char *ArgIn[] = {"xfa", "NULL", ".xfa file name to read (contains mesh)",
		   "exact", "NULL", "truth solution .data file root",
		   "basis", "NULL", "basis .data file",
		   "rom", "NULL", "Reduced Model (.rom) file name",
		   "inp", "online.inp", "online input file with run info",
		   "\0"};
  char xfaFile[xf_MAXSTRLEN];  
  char exactRoot[xf_MAXSTRLEN];  
  char exactFile[xf_MAXSTRLEN];  
  char basisFile[xf_MAXSTRLEN];
  char romFile[xf_MAXSTRLEN];
  char inpFile[xf_MAXSTRLEN];
  char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
  char line[xf_MAXLINELEN];
  FILE *fid, *finp;
  real Value, eM, eMmax;
  xf_DataSet *DataSet = NULL;
  xf_Data *D;
  xf_VectorSet *CardSet;
  xf_Vector *U;
  xf_KeyValue KeyValue;
  xf_ReducedModel *RM;
  xf_RunListType *RunList;
  xf_All *All;
  
  xf_printf("\n");
  xf_printf("=== xf_Test: MR Testing  ===\n");
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
    
  // Get xfaFile
  ierr = xf_GetKeyValue(KeyValue, "xfa", xfaFile);
  if (ierr != xf_OK) return ierr;

  // Get exactRoot
  ierr = xf_GetKeyValue(KeyValue, "exact", exactRoot);
  if (ierr != xf_OK) return ierr;

  // Get basisFile
  ierr = xf_GetKeyValue(KeyValue, "basis", basisFile);
  if (ierr != xf_OK) return ierr;

  // Get romFile
  ierr = xf_GetKeyValue(KeyValue, "rom", romFile);
  if (ierr != xf_OK) return ierr;

  // Get inpFile 
  ierr = xf_GetKeyValue(KeyValue, "inp", inpFile);
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;

  
  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  /* Read .xfa file*/
  ierr = xf_Error(xf_ReadAllBinary(xfaFile, All));
  if (ierr!=xf_OK) return ierr;

  
  // create reduced model
  ierr = xf_Error(xf_CreateReducedModel(&RM));
  if (ierr != xf_OK) return ierr;

  // read reduced model
  if ((fid = fopen(romFile, "rb")) == NULL)  return xf_Error(xf_FILE_READ_ERROR);
  ierr = xf_Error(xf_ReadReducedModelBinary(fid, RM));
  if (ierr != xf_OK) return ierr;
  if (fclose(fid)!= 0) return xf_Error(xf_FILE_READ_ERROR);


  // Dynamically load eqnset library
  ierr = xf_Error(xf_LoadEqnSetLibrary(All->EqnSet->EqnSetLibrary));
  if (ierr != xf_OK) return ierr;

  // Register the equation set
  ierr = xf_Error(xf_EqnSetRegister(All->EqnSet));
  if (ierr != xf_OK) return ierr;
 
 
  // read basis data set
  ierr = xf_Error(xf_CreateDataSet(&DataSet));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, basisFile, DataSet));
  if (ierr!=xf_OK) return ierr;


  // DataSet should contain vectorset named "SSnapSet" which contains the cardinal set
  D = DataSet->Head;
  CardSet = NULL;
  while (D != NULL){
    if ((D->Type == xfe_VectorSet) && (strcmp(D->Title, "SSnapSet") == 0)){
      xf_printf("Found Cardinal set.\n");
      CardSet = (xf_VectorSet *) D->Data;
      D->Data = NULL;
      break;
    }
    D = D->Next;
  }
  if (CardSet == NULL) return xf_Error(xf_INPUT_ERROR);

  /*------------------*/
  /* Read input file  */
  /*------------------*/

  // open input file
  if ((finp = fopen(inpFile, "r")) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  // read first line : nRun
  if (fgets(line, xf_MAXLINELEN, finp) == NULL) return xf_Error(xf_FILE_READ_ERROR);
  if (sscanf(line, "%d", &nRun) != 1) return xf_Error(xf_FILE_READ_ERROR);
  
  // allocate RunList
  ierr = xf_Error(xf_Alloc((void **) &RunList, nRun, sizeof(xf_RunListType)));
  if (ierr != xf_OK) return ierr;
  
  for (iRun=0; iRun<nRun; iRun++){
    
    /* iterm */
    if (fgets(line, xf_MAXLINELEN, finp) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (strncmp(line, "iterm", 5)!=0) return xf_Error(xf_FILE_READ_ERROR);
    ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
    if (ierr != xf_OK) return ierr;
    sscanf(value, "%d", &RunList[iRun].iterm);

    /* key */
    if (fgets(line, xf_MAXLINELEN, finp) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (strncmp(line, "key", 3)!=0) return xf_Error(xf_FILE_READ_ERROR);
    ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
    if (ierr != xf_OK) return ierr;
    strcpy(RunList[iRun].resKey, value);

    /* value */
    if (fgets(line, xf_MAXLINELEN, finp) == NULL) return xf_Error(xf_FILE_READ_ERROR);
    if (strncmp(line, "value", 5)!=0) return xf_Error(xf_FILE_READ_ERROR);
    ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
    if (ierr != xf_OK) return ierr;
    strcpy(RunList[iRun].resValue, value);
    
    if (feof(finp)) return xf_Error(xf_FILE_READ_ERROR);
  } // iRun

  fclose(finp);


  // loop over runs in online.inp, read exact data, calculate max error

  // loop over runs
  eMmax = 0.0;
  for (iRun=0; iRun<nRun; iRun++){

    // read exact solution
    sprintf(exactFile, "%s_%d.data", exactRoot, iRun);
    ierr = xf_Error(xf_DestroyDataSet(DataSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_CreateDataSet(&DataSet));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ReadDataSetBinary(NULL, NULL, exactFile, DataSet));
    if (ierr!=xf_OK) return ierr;
     
    // find state in exact solution
    ierr = xf_Error(xf_FindPrimalState(DataSet, 0, &D, NULL));
    if (ierr != xf_OK) return ierr;
    U = (xf_Vector *) D->Data;

    // calculate nonlinear error
    ierr = xf_Error(xf_ComputeNonlinearError(All, U, RM, CardSet, &eM));
    if (ierr != xf_OK) return ierr;

    xf_printf(" %d  eM = %.10E\n", iRun, eM);

    eMmax = max(eMmax, eM);

  } // iRun

  xf_printf("eMmax = %.10E\n", eMmax);

  xf_Release( (void *) RunList);

  // destroy CardSet
  ierr = xf_Error(xf_DestroyVectorSet(CardSet));
  if (ierr != xf_OK) return ierr; 

  /* Destroy Reduced Model */
  ierr = xf_Error(xf_DestroyReducedModel(RM, xfe_True));
  if (ierr != xf_OK) return ierr;

  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;

  /* Destroy dataset */
  ierr = xf_Error(xf_DestroyDataSet(DataSet));
  if (ierr != xf_OK) return ierr;


  xf_printf("xf_Test finished.\n");

  return xf_OK;
}
