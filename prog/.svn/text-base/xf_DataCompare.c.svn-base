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
  FILE:  xf_DataCompare.c

  This program compares two data solutions on the same mesh

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
#include "xf_Solver.h"
#include "xf_Quad.h"
#include "xf_QuadRule.h"



/******************************************************************/
//   FUNCTION Definition: xf_AnalyticalFcn
static int 
xf_AnalyticalFcn(real *xglob, real *F)
{
  real x, y;
  
  // various hardcoded analytical functions can go here
  x = xglob[0];
  y = xglob[1];

  // SCALAR
  //(*F) = exp(.1*sin(5.1*x-6.2*y)+0.3*cos(4.3*x+3.4*y));

  // FORCED POISEUILLE
/*   F[0] = 1.0; */
/*   F[1] = 0.05*sin(2.0*y); */
/*   F[2] = 0.0; */
/*   F[3] = 1/.4 + .5*F[1]*F[1]/F[0]; */
  
  // POISEUILLE
 /*  F[0] = 1.0; */
/*   F[1] = .5*.05*(1-y*y); */
/*   F[2] = 0.0; */
/*   F[3] = (1-.05*x)/.4 + .5*F[1]*F[1]/F[0]; */

  // MANUFACTURED
  F[0] = 1.0 + .05*sin(4*x+3*y);
  F[1] = F[0]*(.1-.1*cos(2*x+4*y));
  F[2] = F[0]*(.05+.02*cos(3*x+6*y));
  F[3] = (1.0+.04*sin(5*x-7*y))/.4 + .5*F[1]*F[1]/F[0] + .5*F[2]*F[2]/F[0];


  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CalculateErrorNorm
static int 
xf_CalculateErrorNorm(xf_All *All, xf_Vector *U1, xf_Vector *U2, 
		      enum xfe_Bool Analytical, enum xfe_Bool PrintStats,
		      real *penorm2)
{
  int ierr, dim, sr, k, d, iq, nq, pnq;
  int egrp, elem, nu, ns;
  int Order1, Order2, QuadOrder;
  enum xfe_BasisType Basis1, Basis2;
  enum xfe_Bool QuadChanged;
  real *EU1, *EU2, *xq, *u1, *u2, *wq, enorm2;
  real *uA = NULL, *xglob;
  real xcenter1[3] = {0., 0., 0.};
  real xcenter2[3] = {0., 0., 0.};
  real min1=1e30, max1=-1e30, min2=1e30, max2=-1e30;
  real xmin1[3], xmax1[3], xmin2[3], xmax2[3];
  real e1, e2, enormA;
  real integral1, integral2;
  xf_QuadData *QuadData;
  xf_BasisData *PhiData1;
  xf_BasisData *PhiData2;
  xf_BasisData *GeomPhiData;
  xf_JacobianData *JData;
  xf_Mesh *Mesh;


  Mesh = All->Mesh;
  dim  = Mesh->Dim;

  // both vectors should have the same state rank
  sr   = U1->StateRank;
  if (sr != U2->StateRank) return xf_Error(xf_INPUT_ERROR);

  ierr = xf_Error(xf_Alloc( (void **) &uA, sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  QuadData = NULL;
  PhiData1 = NULL;
  PhiData2 = NULL;
  JData    = NULL;
  u1       = NULL;
  u2       = NULL;
  wq       = NULL;
  xglob    = NULL;
  GeomPhiData = NULL;
  enorm2   = 0;
  enormA   = 0;
  e1 = e2  = 0;
  integral1 = integral2 = 0.;
  pnq = -1;
  for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

    Basis1 = U1->Basis[egrp]; Order1 = U1->Order[egrp];
    Basis2 = U2->Basis[egrp]; Order2 = U2->Order[egrp];

    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
      /* Pull off quad points for the element; will not recalculate in generic case */
      QuadOrder = Order1 + Order2;
      ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
      if (ierr != xf_OK) return ierr;

      nq = QuadData->nquad;
      xq = QuadData->xquad;

      // compute basis functions for U1 vector
      ierr = xf_Error(xf_EvalBasis(Basis1, Order1, QuadChanged, 
				   nq, xq, xfb_Phi, &PhiData1));
      if (ierr != xf_OK) return ierr;

      // compute basis functions for U2 vector
      ierr = xf_Error(xf_EvalBasis(Basis2, Order2, QuadChanged, 
				   nq, xq, xfb_Phi, &PhiData2));
      if (ierr != xf_OK) return ierr;


      // element Jacobian
      ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
      if (ierr != xf_OK) return ierr;


      // re-allocate memory if quad points increased
      if (nq > pnq){
	ierr = xf_Error(xf_ReAlloc( (void **)  &u1, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **)  &u2, nq*sr, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
	if (ierr != xf_OK) return ierr;
	ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
	if (ierr != xf_OK) return ierr;
      }

      // form detJ-multiplied quad weight vector, wq
      for (iq=0; iq<nq; iq++) 
	wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
      
      EU1 = U1->GenArray[egrp].rValue[elem];
      EU2 = U2->GenArray[egrp].rValue[elem];

      xf_MxM_Set(PhiData1->Phi, EU1, nq, PhiData1->nn, sr, u1); // interpolate u1
      xf_MxM_Set(PhiData2->Phi, EU2, nq, PhiData2->nn, sr, u2); // interpolate u2

      // add to square error norm via quadrature
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  enorm2 += (u1[iq*sr+k]-u2[iq*sr+k])*(u1[iq*sr+k]-u2[iq*sr+k])*wq[iq];
      
      // obtain global coords of quad points
      ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged, 
				      nq, xq, xglob));
      if (ierr != xf_OK) return ierr;
      
      // Analytical error if requested
      if (Analytical){
	//for (iq=0; iq<nq; iq++)
	for (iq=0; iq<nq; iq++){
	  ierr = xf_Error(xf_AnalyticalFcn(xglob+iq*dim, uA));
	  if (ierr != xf_OK) return ierr;

	  for (k=0; k<sr; k++)
	    enormA += (u1[iq*sr+k]-uA[k])*(u1[iq*sr+k]-uA[k])*wq[iq];
	  
	  //xf_printf("%.10E %.10E %.10E\n", xglob[0], xglob[1], u1[iq*sr+1] - uA[1]);
	}
	
      }

      // add to integrals of data1 and data2
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  integral1 += (u1[iq*sr+k])*wq[iq];
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  integral2 += (u2[iq*sr+k])*wq[iq];

      // add to integrals for xcenter1 and xcenter2
      for (d=0; d<dim; d++)
	for (iq=0; iq<nq; iq++)
	  for (k=0; k<sr; k++)
	    xcenter1[d] += (u1[iq*sr+k])*wq[iq]*xglob[dim*iq+d];
      for (d=0; d<dim; d++)
	for (iq=0; iq<nq; iq++)
	  for (k=0; k<sr; k++)
	    xcenter2[d] += (u2[iq*sr+k])*wq[iq]*xglob[dim*iq+d];
      
      // take min and max
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++){
	  if (u1[iq*sr+k] < min1) for(d=0;d<dim;d++) xmin1[d]=xglob[dim*iq+d];
	  if (u1[iq*sr+k] > max1) for(d=0;d<dim;d++) xmax1[d]=xglob[dim*iq+d];
	  min1 = min(u1[iq*sr+k], min1);
	  max1 = max(u1[iq*sr+k], max1);
	}
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++){
	  if (u2[iq*sr+k] < min2) for(d=0;d<dim;d++) xmin2[d]=xglob[dim*iq+d];
	  if (u2[iq*sr+k] > max2) for(d=0;d<dim;d++) xmax2[d]=xglob[dim*iq+d];
	  min2 = min(u2[iq*sr+k], min2);
	  max2 = max(u2[iq*sr+k], max2);
	}

	  

      // add to square norms of data1 and data2
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  e1 += (u1[iq*sr+k])*(u1[iq*sr+k])*wq[iq];
      for (iq=0; iq<nq; iq++)
	for (k=0; k<sr; k++)
	  e2 += (u2[iq*sr+k])*(u2[iq*sr+k])*wq[iq];

      pnq = nq;
    } // elem
  } // egrp
  
  if (PrintStats){
    xf_printf("Centroid of data1 = ");
    if (integral1 == 0.) xf_printf("UNDEFINED\n");
    else for (d=0; d<dim; d++) xf_printf(" %.8E", xcenter1[d]/integral1);
    xf_printf("\n");
    xf_printf("Centroid of data2 = ");
    if (integral2 == 0.) xf_printf("UNDEFINED\n");
    else for (d=0; d<dim; d++) xf_printf(" %.8E", xcenter2[d]/integral2);
    xf_printf("\n");
    xf_printf("Min of data1 = %.10E at ", min1);
    for (d=0; d<dim; d++) xf_printf(" %.6E", xmin1[d]);
    xf_printf("\nMax of data1 = %.10E at ", max1);
    for (d=0; d<dim; d++) xf_printf(" %.6E", xmax1[d]);
    xf_printf("\nMin of data2 = %.10E at ", min2);
    for (d=0; d<dim; d++) xf_printf(" %.6E", xmin2[d]);
    xf_printf("\nMax of data2 = %.10E at ", max2);
    for (d=0; d<dim; d++) xf_printf(" %.6E", xmax2[d]);
    xf_printf("\n");
  }
  if (penorm2 != NULL) (*penorm2) = enorm2;
  else{
    xf_printf("Integral of data1 = %.15E\n", integral1);
    xf_printf("Integral of data2 = %.15E\n", integral2);
    xf_printf("Norm of data1 = %.15E\n", sqrt(e1));
    xf_printf("Norm of data2 = %.15E\n", sqrt(e2));
    xf_printf("Continuous L2 error norm = %.15E\n", sqrt(enorm2));
    if (Analytical)
      xf_printf("Analytical L2 error norm:\n%.15E\n", sqrt(enormA));
  }
  


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


  // release memory
  xf_Release( (void *) u1);
  xf_Release( (void *) u2);
  xf_Release( (void *) wq);
  xf_Release( (void *) xglob);
  xf_Release( (void *) uA);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_UnsteadyErrorNorm
static int 
xf_UnsteadyErrorNorm(xf_All *All, xf_TimeHistData *TimeHistData,
		     const char *data1File, const char *data2File,
		     int *itoffset, enum xfe_Bool TakeDiff)
{
  int ierr;
  int nTime, iTime;
  int iOrder;
  int OrderTime, QuadOrder, iq, nq;
  int nvector1, nvector2;
  enum xfe_TimeSchemeType TimeScheme;
  char dataFile[xf_MAXSTRLEN];
  real TimeStep, Time;
  real enorm2, enorm2loc, wtemp;
  real *tq = NULL, *wq = NULL;
  xf_Vector *Ui1[xf_MAXDGTIMENODE], *Ui2[xf_MAXDGTIMENODE];
  xf_Vector *U1, *U2;
  xf_DataSet *DataSet1, *DataSet2;
  xf_Data *D;


  nTime = TimeHistData->nTime; // number of slabs
  
  enorm2 = 0.0;
  
  for (iTime=0; iTime<nTime; iTime++){

    // load data1
    ierr = xf_Error(xf_CreateDataSet(&DataSet1));
    if (ierr != xf_OK) return ierr;
    sprintf(dataFile, "%s%d.data", data1File, iTime + itoffset[0]);
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, dataFile, DataSet1));
    if (ierr!=xf_OK) return ierr;
    // use data vectors in order of data
    for (iOrder=0, D=DataSet1->Head; D!=NULL; iOrder++, D=D->Next){
      Ui1[iOrder] = (xf_Vector *) D->Data;
      xf_printf("U1: Title = %s\n", D->Title);
    }
    nvector1 = iOrder;

    // load data2
    ierr = xf_Error(xf_CreateDataSet(&DataSet2));
    if (ierr != xf_OK) return ierr;
    sprintf(dataFile, "%s%d.data", data2File, iTime + itoffset[1]);
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, dataFile, DataSet2));
    if (ierr!=xf_OK) return ierr;
    // use data vectors in order of data
    for (iOrder=0, D=DataSet2->Head; D!=NULL; iOrder++, D=D->Next){
      Ui2[iOrder] = (xf_Vector *) D->Data;
      xf_printf("U2: Title = %s\n", D->Title);
    }
    nvector2 = iOrder;

    // find temporary vectors
    ierr = xf_Error(xf_FindSimilarVector(All, Ui1[0], "Utemp1", 
					 xfe_True, xfe_True, 
					 NULL, &U1, NULL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_FindSimilarVector(All, Ui2[0], "Utemp2", 
					 xfe_True, xfe_True, 
					 NULL, &U2, NULL));
    if (ierr != xf_OK) return ierr;


    TimeScheme = TimeHistData->TimeScheme[iTime];
    TimeStep   = TimeHistData->TimeStep[iTime];
    Time       = TimeHistData->Time[iTime];

    // determine order in time
    ierr = xf_Error(xf_DGTimeScheme2Order(TimeScheme, &OrderTime));
    if (ierr != xf_OK) return ierr;

    // verify that have the correct number of vectors
    if ((nvector1 != nvector2) || (nvector1 != (OrderTime+1)))
      return xf_Error(xf_INPUT_ERROR);

    QuadOrder = 2*OrderTime+1; // hardcoded
      
    // quadrature points on time slab
    tq = NULL;
    wq = NULL;
    ierr = xf_Error(xf_QuadLine(QuadOrder, &nq, &tq, &wq));
    if (ierr != xf_OK) return ierr;

    // loop over quad points
    for (iq=0; iq<nq; iq++){
      // quadrature weight at this point
      wtemp =  wq[iq]*TimeStep;

      // Interpolate states
      ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Ui1, 1, tq+iq, &U1, NULL));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_DGTimeInterpolate(TimeScheme, Ui2, 1, tq+iq, &U2, NULL));
      if (ierr != xf_OK) return ierr;

      // calculate error
      ierr = xf_Error(xf_CalculateErrorNorm(All, U1, U2, xfe_False, xfe_False, &enorm2loc));
      if (ierr != xf_OK) return ierr;

      xf_printf("iTime = %d, iq = %d, tq = %.4f, enorm2loc = %.4e\n", 
		iTime, iq, tq[iq], enorm2loc);

      // add to running total
      enorm2 += wtemp*enorm2loc;

    } // iq
      
    xf_Release( (void *) tq);
    xf_Release( (void *) wq);

    // write out difference if requested
    if (TakeDiff){
      // Ui1 -= Ui2
      for (iOrder=0; iOrder<=OrderTime; iOrder++){
	ierr = xf_Error(xf_SetVector(Ui2[iOrder], xfe_Sub, Ui1[iOrder]));
	if (ierr != xf_OK) return ierr;
      }
      // write out Ui1
      sprintf(dataFile, "diff%d.data", iTime);
      ierr = xf_Error(xf_WriteDataSetBinary(All->Mesh, DataSet1, NULL, dataFile));
      if (ierr != xf_OK) return ierr;
    }


    // destroy DataSet1
    ierr = xf_Error(xf_DestroyDataSet(DataSet1));
    if (ierr != xf_OK) return ierr;
      
    // destroy DataSet2
    ierr = xf_Error(xf_DestroyDataSet(DataSet2));
    if (ierr != xf_OK) return ierr;
      
  } // iTime

  xf_printf("Space-time continuous L2 error norm = %.10E\n", sqrt(enorm2));

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: main
int 
main(int argc, char *argv[])
{
  int ierr, i;
  int itoffset[2];
  enum xfe_Bool Analytical, TakeDiff, PrintStats, doUnsteady;
  char *ArgIn[] = {"data1", "NULL", "first .data file",
		   "data2", "NULL", "second .data file",
		   "xfa", "NULL", ".xfa file containing the mesh",
		   "diff", "False", "Write out difference to diff.data",
		   "stats", "False", "if true, statistics of states will be printed",
		   "Analytical", "False", "is an analytical error norm desired?",
		   "TimeHist", "None", "time history file for unsteady L2 error calc",
		   "itoffset1", "0", "offset relative to ind in TimeHist",
		   "itoffset2", "0", "offset relative to ind in TimeHist",
		   "\0"};
  char data1File[xf_MAXSTRLEN];
  char data2File[xf_MAXSTRLEN];
  char xfaFile[xf_MAXSTRLEN];
  char TimeHistFile[xf_MAXSTRLEN];
  char line[xf_MAXLINELEN];
  FILE *fid;
  xf_KeyValue KeyValue;
  xf_TimeHistData *TimeHistData = NULL;
  xf_Vector *U1, *U2;
  xf_DataSet *DataSet1, *DataSet2;
  xf_Data *D1, *D2;
  xf_All *All;
  
  xf_printf("\n");
  xf_printf("=== Data Comparison ===\n");
  xf_printf("\n");
    
      
  // initialize key-value
  ierr = xf_Error(xf_InitKeyValue(&KeyValue));
  if (ierr != xf_OK) return ierr;

  // parse arguments
  ierr = xf_ParseArg(ArgIn, argc, argv, &KeyValue);
  if (ierr == xf_FORCE_QUIT) return xf_OK;
  if (ierr != xf_OK) return xf_Error(ierr);

  // Get data1
  ierr = xf_GetKeyValue(KeyValue, "data1", data1File);
  if (ierr != xf_OK) return ierr;

  // Get data2
  ierr = xf_GetKeyValue(KeyValue, "data2", data2File);
  if (ierr != xf_OK) return ierr;

  // Get xfaFile
  ierr = xf_GetKeyValue(KeyValue, "xfa", xfaFile);
  if (ierr != xf_OK) return ierr;

  /* TakeDiff? */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValue, "diff", &TakeDiff));
  if (ierr != xf_OK) return ierr;

  /* Print stats? */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValue, "stats", &PrintStats));
  if (ierr != xf_OK) return ierr;

  /* Analytical */
  ierr = xf_Error(xf_GetKeyValueBool(KeyValue, "Analytical", &Analytical));
  if (ierr != xf_OK) return ierr;

  // Get TimeHistFile
  ierr = xf_GetKeyValue(KeyValue, "TimeHist", TimeHistFile);
  if (ierr != xf_OK) return ierr;

  /* itoffset1 */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValue, "itoffset1", itoffset+0));
  if (ierr != xf_OK) return ierr;

  /* itoffset2 */
  ierr = xf_Error(xf_GetKeyValueInt(KeyValue, "itoffset2", itoffset+1));
  if (ierr != xf_OK) return ierr;

  // destroy key-value
  ierr = xf_Error(xf_DestroyKeyValue(&KeyValue));
  if (ierr!=xf_OK) return ierr;


  /* Create .xfa structure */
  ierr = xf_Error(xf_CreateAll(&All, xfe_False));
  if (ierr != xf_OK) return ierr;
  
  /* Read .xfa file */
  xf_printf("Reading .xfa file\n");
  ierr = xf_Error(xf_ReadAllBinary(xfaFile, All));
  if (ierr!=xf_OK) return ierr;

  // read input time history if specified
  doUnsteady = xfe_False;
  if (xf_NotNull(TimeHistFile)){
    doUnsteady = xfe_True;
    ierr = xf_Error(xf_ReadTimeHistData(TimeHistFile, NULL, &TimeHistData));
    if (ierr != xf_OK) return ierr;
  }


  // Treat unsteady mode separately
  if (doUnsteady){

    ierr = xf_Error(xf_UnsteadyErrorNorm(All, TimeHistData, data1File, data2File,
					 itoffset, TakeDiff));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_DestroyTimeHistData(TimeHistData));
    if (ierr != xf_OK) return ierr;

  }
  else{
    /* Read data1 */
    xf_printf("Reading first .data file\n");
    ierr = xf_Error(xf_CreateDataSet(&DataSet1));
    if (ierr != xf_OK) return ierr;
  
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, data1File, DataSet1));
    if (ierr!=xf_OK) return ierr;
  
    ierr = xf_FindPrimalState(DataSet1, 0, &D1, NULL);
    if (ierr == xf_OK){
      U1 = (xf_Vector *) D1->Data;
    }
    else if (ierr == xf_NOT_FOUND){
      U1 = (xf_Vector *) DataSet1->Head->Data;
    }
    else if (ierr != xf_OK) return xf_Error(ierr);


    /* Read data2 */
    xf_printf("Reading second .data file\n");
    ierr = xf_Error(xf_CreateDataSet(&DataSet2));
    if (ierr != xf_OK) return ierr;
  
    ierr = xf_Error(xf_ReadDataSetBinary(All->Mesh, NULL, data2File, DataSet2));
    if (ierr!=xf_OK) return ierr;
  
    ierr = xf_FindPrimalState(DataSet2, 0, &D2, NULL);
    if (ierr == xf_OK){
      U2 = (xf_Vector *) D2->Data;
    }
    else if (ierr == xf_NOT_FOUND){
      U2 = (xf_Vector *) DataSet2->Head->Data;
    }
    else if (ierr != xf_OK) return xf_Error(ierr);


    /* Compute (and print out) L2 error norm */
    ierr = xf_Error(xf_CalculateErrorNorm(All, U1, U2, Analytical, PrintStats, NULL));
    if (ierr != xf_OK) return ierr;


    // take a difference if desired
    if (TakeDiff){
      ierr = xf_Error(xf_VectorMultSet(U2, 1.0, xfe_Sub, U1));
      if (ierr != xf_OK) return ierr;
    
      ierr = xf_Error(xf_DumpVectorBinary(All->Mesh, "Vector", U1, "diff.data"));
      if (ierr != xf_OK) return ierr;
    }

    // destroy DataSet1
    ierr = xf_Error(xf_DestroyDataSet(DataSet1));
    if (ierr != xf_OK) return ierr;
    
    // destroy DataSet2
    ierr = xf_Error(xf_DestroyDataSet(DataSet2));
    if (ierr != xf_OK) return ierr;
  }
    
  /* Destroy .xfa structure */
  ierr = xf_Error(xf_DestroyAll(All));
  if (ierr!=xf_OK) return ierr;

  xf_printf("xf_DataCompare finished.\n");

  return xf_OK;
}
