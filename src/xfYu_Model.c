/* account for how to implement my model */
#include "xf.h"
#include "xf_AllStruct.h"
#include "xf_MPI.h"
#include "xf_Memory.h"
#include "xf_Basis.h"
#include "xf_All.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xfYu_Model.h"
#include "xf_Math.h"
#include "xf_Quad.h"
#include "xf_String.h"
#include "xf_Param.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "reaction_routine.h"
//#include "homo0D.h"
#include "xfYu_EntropyBounding.h"

//distributed source code
#include "./xfYu_Initialization.c"
#include "./xfYu_RiemannSolver.c"

static int init_b();

/********************************************************************************************/
//function: for gamma evaluation
int
Yu_GammaValueEvaluate(const real p, const real rho, real *Y, const int nspecies, real *GammaEval,
                      enum xfe_Bool DetailChem)
{
   //Y1: air (gamma 1.4, Cv 0.72)
   //Y2: He+0.28air (gamma 1.648, Cv 2.44)

   real Cp[9], Cv[9];

   if(DetailChem == xfe_True)
   {
   //   chemkin_backoutGamma(rho, p, Y, nspecies, GammaEval);
   }
   else
   {
      Cv[0] = 0.72;      Cv[1] = 2.44;
      Cp[0] = 1.4*Cv[0]; 
      Cp[1] = 1.648*Cv[1];

      if(Y[0]<0.0) { Y[0] = 0.0; Y[1] = 1.0; }
      if(Y[0]>1.0) { Y[0] = 1.0; Y[1] = 0.0; }

      (*GammaEval) = (Y[0]*Cp[0] + Y[1]*Cp[1])/(Y[0]*Cv[0] + Y[1]*Cv[1]);
      // (*GammaEval) = 1.4;
   }

   return xf_OK;
}
/********************************************************************************************/
//function: include user-defined vector in ouput data
int
Yu_OutputVectorUpdate(xf_All *All, Yu_Model *Model, xf_Vector *Gamma, xf_Vector *OutVec, 
                      xf_Vector *State)
{
    int ierr, dim, i;
    int iegrp, ielem;
    real rho, u[3], p, *Y, sumY;
    real *EU, *EGamma, *pOutVec;
    enum xfe_BasisType Basis;
    xf_Mesh *Mesh;
    xf_Data *D;
    xf_Vector *P0State;
     
    Mesh = All->Mesh;
    dim  = Mesh->Dim;

    //project state vector down to P0
    ierr = xf_Error(xf_FindSimilarVector(All, State, "GammaComput", xfe_False, xfe_True, NULL, &P0State, NULL));
    if(ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_SetVector(State, xfe_Set, P0State));
    if(ierr != xf_OK) return ierr;

    //find proper basis; keep identical;
    Basis = State->Basis[0];
    
    ierr = xf_Error(xf_ProjectVectorInPlace(Mesh, All->DataSet, P0State, NULL, Basis, NULL, NULL, 0, 0)); 
    if(ierr != xf_OK) return ierr;
    
    //loop over each element to compute gamma out
    for (iegrp=0; iegrp<Mesh->nElemGroup; iegrp++){
        
        for (ielem=0; ielem<Mesh->ElemGroup[iegrp].nElem; ielem++){
            
            EU      = P0State->GenArray[iegrp].rValue[ielem];  //one set of variable on each elem since P0
            EGamma  = Gamma->GenArray[iegrp].rValue[ielem];
            pOutVec = OutVec->GenArray[iegrp].rValue[ielem];

            rho = EU[0];
            if(dim == 2)
            { 
               u[0] = EU[1]/EU[0]; u[1] = EU[2]/EU[0];
               p = (EGamma[0] - 1.0)*(EU[3] - 0.5*rho*(u[0]*u[0] + u[1]*u[1])); 
            }
            else if(dim == 3)
            { 
               u[0] = EU[1]/EU[0]; u[1] = EU[2]/EU[0];  u[2] = EU[3]/EU[0]; 
               p = (EGamma[0] - 1.0)*(EU[4] - 0.5*rho*(u[0]*u[0] + u[1]*u[1] + u[2]*u[2])); 
            }
            else
            {
            return xf_NOT_SUPPORTED;
            }

            //specify algorithms for output vector
            if(pOutVec[0] < p)
               pOutVec[0] = p;
  
         }//ielem
      }//iegrp
 
    return xf_OK;
}
/********************************************************************************************/
//function: to update heat capacity ratio when needed
int
Yu_GammaVectorUpdate(xf_All *All, Yu_Model *Model, xf_Vector *Gamma, xf_Vector *State)
{
    int ierr, dim, i;
    int iegrp, ielem;
    real rho, u, v, p, *Y, sumY, GammaEval;
    real *EU, *EGamma;
    xf_Mesh *Mesh;
    xf_Data *D;
    xf_Vector *P0State;
    
    Mesh = All->Mesh;
    dim  = Mesh->Dim;
    
    //only solve for N-1 species and evaluate extra one using 1.0-SumY[i]
    ierr = xf_Error(xf_Alloc( (void **) &Y, (Model->nVars - 4) + 1, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    //project state vector down to P0
    ierr = xf_Error(xf_FindSimilarVector(All, State, "GammaComput", xfe_False, xfe_True, NULL, &P0State, NULL));
    if(ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_SetVector(State, xfe_Set, P0State));
    if(ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_ProjectVectorInPlace(Mesh, All->DataSet, P0State, NULL, xfe_QuadLagrange, NULL, NULL, 0, 0)); 
    if(ierr != xf_OK) return ierr;
    
    //loop over each element to compute gamma out
    for (iegrp=0; iegrp<Mesh->nElemGroup; iegrp++){
        
        for (ielem=0; ielem<Mesh->ElemGroup[iegrp].nElem; ielem++){
            
            EU     = P0State->GenArray[iegrp].rValue[ielem];  //one set of variable on each elem since P0
            EGamma = Gamma->GenArray[iegrp].rValue[ielem]; 
            
            rho = EU[0];
            u   = EU[1]/EU[0];
            v   = EU[2]/EU[0];
            p   = (EGamma[0] - 1.0)*(EU[3] - 0.5*(EU[1]*EU[1]+EU[2]*EU[2])/EU[0]);
            sumY= 0.0;
            for(i=0; i<(Model->nVars - 4); i++)
            {
                Y[i]  = EU[4+i]/EU[0];
                sumY += Y[i];
            }
            Y[Model->nVars-4] = 1.0 - sumY; 
            
            //put the gamma compute routine here
            Yu_GammaValueEvaluate(p, rho, Y, (Model->nVars - 4) + 1, &GammaEval, Model->DetailChem);
            //only test now
            EGamma[0] = GammaEval;
            
        } // ielem
    } // iegrp
    
    
    xf_Release( (void *) Y);
    
    return xf_OK;
}

/********************************************************************************************/
//function: for energy update according to pressure
int 
Yu_EnergyCorrection(xf_All *All, Yu_Model *Model, xf_Vector *State)
{
   int ierr, i, j, k, dim, nn, egrp, elem, sr;
   int iq, nq, pnq, pOrder, Order, QuadOrder;
   int *P;
   real *Q, *R;
   real *xq, *xglob, x0[3] = {0.}, *qU, *U, *EU;
   real GammaNew, GammaOld;
   real rho, u, v, p;
   enum xfe_Bool QuadChanged;
   enum xfe_BasisType Basis;
   xf_BasisData *PhiData, *GeomPhiData;
   xf_QuadData *QuadData;
   xf_Mesh *Mesh;
   xf_Vector *GammaOldVec, *GammaNewVec;
   xf_Data   *GammaDat;

   Mesh = All->Mesh;
   dim  = Mesh->Dim;
   sr   = Model->nVars;

   //find Gamma Vector
   ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
   if(ierr == xf_NOT_FOUND)
   {
      xf_printf("Cannot find heat capacity ratio...\n");
      return ierr;
   }
   else
      GammaNewVec = (xf_Vector *) GammaDat->Data;

   ierr = xf_Error(xf_FindSimilarVector(All, GammaNewVec, "HeatCapacityRatio_Temp", xfe_False, xfe_True, NULL, &GammaOldVec, NULL));
   if(ierr != xf_OK) return ierr;
   
   ierr = xf_Error(xf_SetVector(GammaNewVec, xfe_Set, GammaOldVec));
   if(ierr != xf_OK) return ierr;
    
   //update Gamma Vector
   ierr = xf_Error(Yu_GammaVectorUpdate(All, Model, GammaNewVec, State));
   if(ierr != xf_OK) return ierr;
    
   if(State->Linkage != xfe_LinkageGlobElem) return xf_Error(xf_NOT_SUPPORTED);

   /* Here we perform least-squares projection for energy*/
   QuadData    = NULL;
   PhiData     = NULL;
   qU          = NULL;
   xglob       = NULL;
   Q           = NULL;
   R           = NULL;
   P           = NULL;

   for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      Basis = State->Basis[egrp];

      pOrder = -1;
      pnq    = -1; // here so that we reallocate at least once for every element group

      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

         EU = State->GenArray[egrp].rValue[elem];

         GammaOld = GammaOldVec->GenArray[egrp].rValue[elem][0];
         GammaNew = GammaNewVec->GenArray[egrp].rValue[elem][0];
        
        

         if(fabs(GammaOld - GammaNew) < Model->Gammathreshold)
            continue;

         //get interpolation order
         Order = xf_InterpOrder(State, egrp, elem);

         if (Order != pOrder){
         
            pOrder = Order;

            // determine nn = # unknowns for elements in this group
            ierr = xf_Error(xf_Order2nNode(Basis, Order, &nn));
            if (ierr != xf_OK) return ierr;

            // determine quadrature order (not using eqnset)
            ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order, &QuadOrder));
            if (ierr != xf_OK) return ierr;
         }

         /* Pull off quad points for the element; make sure have enough
         quad points (more than nn). */
         ierr = xf_Error(xf_QuadElemAtLeast(Mesh, egrp, elem, QuadOrder, nn, &QuadData, &QuadChanged));
         if (ierr != xf_OK) return ierr;
         nq = QuadData->nquad;
         xq = QuadData->xquad;

         // request basis fcn for each interpolate points
         ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi, &PhiData));
         if (ierr != xf_OK) return ierr;

         // re-allocate data if quad points increased
         if (nq > pnq){
            ierr = xf_Error(xf_ReAlloc( (void **) &qU, nq*sr, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
            if (ierr != xf_OK) return ierr;
         }

         // interpolate state at interpolate points
         xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, qU);

         //specify how energy will change
          for(k=0; k<nq; k++)
          {
              U = qU + k*sr;
              
              rho = U[0];
              u   = U[1]/U[0];
              v   = U[2]/U[0];
              p   = (GammaOld - 1.0)*(U[3] - 0.5*rho*(u*u + v*v));
              
              U[3] = p/(GammaNew - 1.0) + 0.5*rho*(u*u + v*v);
          }

         // Project u to get coefficients EU
         ierr = xf_Error(xf_ProjectOnElemQR_OnlyEnergy(Mesh, egrp, elem, sr, QuadData,
                                                       Basis, Order, QuadChanged, &Q, &R, &P,
                                                       qU, State->GenArray[egrp].rValue[elem]));
         //ierr = xf_Error(xf_ProjectOnElemLeastSquares_OnlyEnergy(Mesh, egrp, elem, sr, QuadData,
         //                                                        Basis, Order, QuadChanged, &Q, &R, &P,
         //                                                        qU, State->GenArray[egrp].rValue[elem]));
         if (ierr != xf_OK) return ierr;

         pnq = nq;
      } // elem
   } //egrp

   // Only destroy QuadData if points are generic
   ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
   if (ierr != xf_OK) return ierr;

   /* Destroy Basis Data */
   ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
   if (ierr != xf_OK) return ierr;

   xf_Release((void *) qU);
   xf_Release((void *) xglob);
   xf_Release((void *) Q);
   xf_Release((void *) R);
   xf_Release((void *) P);

   return xf_OK;
}

/******************************************************************/
// functions regarding double flux method
// (gamma is piecewise constant serving as a supporting vector with P0 now)
int
Yu_GammaVectorCreate(xf_All *All, Yu_Model *Model, xf_Vector *State, const char VectorName[], const real ConstInitVal)
{
   int ierr, i, d, dim, sr;
   int myRank, nProc, nquans;
   int negrp, egrp;
   int elem, face;
   int GammaOrder, len;
   int ibface, negrphalo;
   int nface, nfacetot;
   int ibfgrp, nbfgrp;
   int iface;
   int *OrderVec = NULL;

   enum xfe_ShapeType Shape, FShape;
   enum xfe_BasisType FBasis;
   enum xfe_Bool ParallelFlag;
   enum xfe_Bool Found, InitFlag;
   enum xfe_Bool Writable;
   enum xfe_BasisType *BasisVec = NULL;

   xf_Data      *D;
   xf_Vector    *GammaTemp;
   xf_Mesh      *Mesh;
   xf_ElemGroup *EG;

   Mesh = All->Mesh;
   dim  = Mesh->Dim;
   sr   = Model->nVars;

   //working with the specified vectors
   xf_printf("Calculating %s ... ", VectorName); fflush(stdout);

   //obtain number of processors
   ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
   if (ierr != xf_OK) return ierr;

   ParallelFlag = (nProc > 1);

   // get order of gamma variable (now piecewise constant)
   GammaOrder = 0;
   // get whether gamma is writable
   Writable = xfe_True;

   // allocate Basis and Order vectors for vector search
   negrp = Mesh->nElemGroup;
   negrphalo = (ParallelFlag ? 2*negrp : negrp);
   ierr = xf_Error(xf_Alloc( (void **) &BasisVec, negrphalo, sizeof(enum xfe_BasisType)));
   if(ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc( (void **) &OrderVec, negrphalo, sizeof(int)));
   if(ierr != xf_OK) return ierr;

   // get Lagrange basis with order of GammaOrder
   for (egrp=0; egrp<negrphalo; egrp++){
      ierr = xf_Error(xf_Basis2UniformLagrange(Mesh->ElemGroup[egrp].QBasis, BasisVec+egrp));
      if (ierr != xf_OK) return ierr;
      OrderVec[egrp] = GammaOrder;
   }

   Found = xfe_False;
   // look for whether it exists already
   nquans = 1;
   len = strlen(VectorName);

   //the first rank: cellwise activiation flag;
   //the 2-sr+1 rank: the AVmodel viscosity value;
   //the last two ranks: entropy lower and upper bound;
   if(strncmp(VectorName, "AVmodel", len) == 0) //needs two quantities for AV model
      nquans = sr+5; 
   
   ierr = xf_Error(xf_FindVector(All, VectorName, xfe_LinkageGlobElem, nquans, 
                                 NULL, 0, 0, BasisVec, OrderVec, NULL, NULL, NULL,
                                 xfe_SizeReal, xfe_True, xfe_True, &D, &GammaTemp, &Found));
   
   if (ierr != xf_OK) return ierr;
   D->ReadWrite = Writable; 

   xf_Release( (void *) BasisVec);
   xf_Release( (void *) OrderVec);

   // set gamma to a constant value now
   if(Found)  //means restart from data merge
   { 
      xf_printf("Vector %s is successfully loaded from restart file~\n", VectorName);
      
   }
   else
   {
      //default zero value setting
      ierr = xf_Error(xf_SetConstVector(GammaTemp, 0, ConstInitVal));
      if(ierr != xf_OK) return ierr;
   }

      len = strlen(VectorName);
      if(Model->EntropyBdFlag && strncmp(VectorName, "MaxPressure", len) == 0)
      {
         //use "MaxPressure" as "MinEntropy" in the cell when entropy
         //bounding is activated; reloaded vector have to be re-initialized;
         //it cannot be re-used.
         if(!Found){
                xf_printf("Entropy Bounding DG is activiated; low bound is intialized!");
                ierr = xf_Error(xf_SetConstVector(GammaTemp, 0, Model->LimM));
                if(ierr != xf_OK) return ierr;
         }
                //finally set the vector to Model struct
                Model->EntropyVec = GammaTemp;
     
         //this is for initialization with flag set to True
         //InitFlag = xfe_True;
         //ierr = xf_Error(Yu_ConductEntropyBounding(All, Model, &State, InitFlag));
         //if(ierr != xf_OK) return ierr;
      }
      // if(Model->GammaVaryFlag)
      // {
      //    ierr = xf_Error(Yu_GammaVectorUpdate(All, Model, GammaTemp, State));
      //    if(ierr != xf_OK) return ierr;
      // }


      //if adaptive time stepping is applied
      len = strlen(VectorName);
      if(strncmp(VectorName, "ElemMaxCharSpeed", len) == 0)
         Model->MaxCharSpeed = GammaTemp;
   
      if(strncmp(VectorName, "ElemMinFaceLen", len) == 0)
         Model->MinFaceLen = GammaTemp;

      if(strncmp(VectorName, "AVmodel", len) == 0)
         Model->AVmodel_data = GammaTemp;

      if(strncmp(VectorName, "RiemannIndicator", len) == 0)
         Model->RiemannIndicator = GammaTemp; 

   xf_printf("done.\n");
   return xf_OK;
}

/******************************************************************/
// read model from input file
// rewritten by Yu for descent clarity; not need for ordering input
static 
enum xfe_Bool TrueOrFalse(const char *value)
{
   if(strcmp(value, "True") == 0)
      return xfe_True;
   if(strcmp(value, "False") == 0)
      return xfe_False;

   xf_printf("Error in model.yu;\nStop run and check!\n");
   getchar();
}
int
PullinModel(Yu_Model *Model)
{
   int i, j, k, tmp, ierr;
   int Num_elem, Num_spe, Num_reac;
   char key[xf_MAXSTRLEN], value[xf_MAXSTRLEN];
   char dummy[xf_MAXSTRLEN], *line;
   char outputkey[xf_MAXSTRLEN], outputvalue[xf_MAXSTRLEN], outputdummy[xf_MAXSTRLEN];
   FILE *fp, *fdata;
   
   //if model.yu file is not exist
   fp = fopen("model.yu","r");
   if(fp == NULL){
      xf_printf("model.yu file does not exist!\n");
      return xf_Error(xf_FILE_READ_ERROR);
   }
   else
      xf_printf("Reading model.yu!\n");

   //sponge layer specified later
   Model->Sponge = NULL;
   /*-------------------------------------------*/
   /* read key and handle its input individually*/
   /*-------------------------------------------*/

   do{

      //not read anything
      if (fgets(dummy, xf_MAXLINELEN, fp) == NULL)
         continue;
      line = dummy; //set pointer

      // if blank line or comment, continue to next line 
      if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;

      //if we read the key; capture the input value
      ierr = xf_Error(xf_ReadKey(line, "=", key, value, xf_MAXSTRLEN));
      if (ierr != xf_OK) return ierr;

      //handle different key in different way:
      if(strcmp(key, "ConvergStudy") == 0){
         Model->ConvergStudyFlag = TrueOrFalse(value);
         continue;
      }

      if(strcmp(key, "dim") == 0){
         Model->dim = atoi(value);
         continue;
      }

      //mesh is curved or not
      if(strcmp(key, "CurvedMesh") == 0){
         Model->TwistFlag = TrueOrFalse(value);
         continue;
      }

      if(strcmp(key, "InterpolateIC") == 0){
         Model->InterpolateInitData = atoi(value);
         continue;
      }

      if(strcmp(key, "basis") == 0){
         strcpy(Model->basis, value);
         xf_printf("Basis %s used!\n", Model->basis);
         continue;
      }

      if(strcmp(key, "timescheme") == 0){
         if(strcmp(value, "FE")==0)           Model->typeTimeScheme = FE;
         else if(strcmp(value, "RK45")==0)    Model->typeTimeScheme = RK45;
         else if(strcmp(value, "SSPRK23")==0) Model->typeTimeScheme = SSPRK23;
         else if(strcmp(value, "SSPRK34")==0) Model->typeTimeScheme = SSPRK34;
         else 
         { printf("time stepping scheme is not defined\n");   return xf_Error(xf_NOT_SUPPORTED); }
   
         printf("Time stepping scheme is %d\n", Model->typeTimeScheme);
         continue;
      }

      if(strcmp(key, "timefactor") == 0){
         Model->time_fac = atof(value);
         continue;
      }

      //polynomial order
      if(strcmp(key, "order") == 0){
         Model->order = atoi(value);
         continue;
      }

      //variable setup
      if(strcmp(key, "Vars") == 0){
         Model->nVars = atoi(value);

         Model->nameVars = (char **)malloc(sizeof(char *) * Model->nVars);
         for(i=0; i<Model->nVars; i++){
            Model->nameVars[i] = (char *)malloc(sizeof(char) * 30);
            fscanf(fp, "%s", Model->nameVars[i]);
         }
         Model->initVars = (real *)malloc(sizeof(real) * Model->nVars);
         for(i=0; i<Model->nVars; i++){
            fscanf(fp, "%lf", &Model->initVars[i]);
         }
         
         //jump one line
         fgets(dummy, xf_MAXLINELEN, fp);
         
         continue;
      }

      //solution limiting setup
      if(strcmp(key, "WhetherLimiting") == 0){
         tmp = atoi(value);
   
         if(tmp == 1){
            Model->LimiterFlag = xfe_True;
            Model->EntropyBdFlag = xfe_False;
         }
         else if(tmp == 2){
            Model->LimiterFlag = xfe_False;
            Model->EntropyBdFlag = xfe_True;
            Model->ConstEntropyBdFlag = xfe_True;
         }
         else if(tmp == 3){
            Model->LimiterFlag = xfe_False;
            Model->EntropyBdFlag = xfe_True;
            Model->ConstEntropyBdFlag = xfe_False;
         }
         else if(tmp == 0){
            Model->LimiterFlag = xfe_False;
            Model->EntropyBdFlag = xfe_False;
         }
         else{
            xf_printf("Parameter Input Error!\n");
            return xf_NOT_SUPPORTED;
         }
         continue;
      }

      //entropy bound
      if(strcmp(key, "LimiterEntropy") == 0){
         Model->LimM = atof(value);
         xf_printf("Entropy bound:%lf\n", Model->LimM);
         continue;
      }

      //whether diffusion is actived
      if(strcmp(key, "DiffFlag") == 0){
         tmp = atoi(value);
            if(tmp == 1)
            {
               Model->DiffFlag   = xfe_True;
               Model->Sutherland = xfe_False;
               Model->AVmodel    = xfe_False;
            }
            else if(tmp == 2)
            {
               Model->DiffFlag   = xfe_True;
               Model->Sutherland = xfe_True;
               Model->AVmodel    = xfe_False;
            }
            else if(tmp == 3)
            {
               Model->DiffFlag   = xfe_True;
               Model->Sutherland = xfe_False;
               Model->AVmodel    = xfe_True;
            }
            else if(tmp == 4)
            {
               Model->DiffFlag   = xfe_True;
               Model->Sutherland = xfe_True;
               Model->AVmodel    = xfe_True;
            }
            else
            {
               Model->DiffFlag   = xfe_False;
               Model->Sutherland = xfe_False;
               Model->AVmodel    = xfe_False;
            }

            //initialize required data structure for AV
            //if(Model->AVmodel){
            //   ierr = xf_Error(Yu_GammaVectorCreate(All, &Model, State, "AVmodel", 0.0));
            //   if(ierr != xf_OK) return ierr;
            //}

            continue;
      }

      //diffusion parameter
       if(strcmp(key, "Ru") == 0){
           //universial gas constant
           Model->Ru =  atof(value);
           continue;
       }
       if(strcmp(key, "Viscosity") == 0){
           Model->mu_c = atof(value);
           continue;
       }
       if(strcmp(key, "Conductivity") == 0){
           Model->kappa_c = atof(value);
           continue;
       }
       if(strcmp(key, "Diffusivity") == 0){
           Model->Diff_c = atof(value);
           continue;
       }
       if(strcmp(key, "MolecularWeight") == 0){
           Model->moleW[0] = atof(value);
           continue;
       }
       
       //read whether we use variable gamma
       if(strcmp(key, "VariableGamma") == 0){
           Model->GammaVaryFlag = TrueOrFalse(value);
           continue;
       }
       
       //read the initial value of gamma
       if(strcmp(key, "InitGammaValue") == 0){
           Model->GammaInit = atof(value);
           continue;
       }
       
       //gamma value threshold at interface
       if(strcmp(key, "Gammathreshold") == 0){
           Model->Gammathreshold = atof(value);
           continue;
       }
       
       if(strcmp(key, "IncludeReaction") == 0){
           tmp = atoi(value);

           Model->ChemSource = xfe_False;
           Model->DetailChem = xfe_False;
           Model->MMSSource  = xfe_False;

           if(tmp == 0)
           {
               //i=0 without chemistry
               Model->ChemSource = xfe_False;
               Model->DetailChem = xfe_False;
               xf_printf("Chemistry is disabled!\n");
           }
           else if(tmp == 1)
           {
               //i=1 with simplified chemistry
               Model->ChemSource = xfe_True;
               Model->DetailChem = xfe_False;
           }
           else if(tmp == 2)
           {
               //i=2 with detailed chemistry
               Model->ChemSource = xfe_True;
               Model->DetailChem = xfe_True;
           }
           else if(tmp == 3)
           {
               //special source term 
               Model->MMSSource = xfe_True;
           }
           else
               return xf_NOT_SUPPORTED;
           
           continue;
       }

       //deal with boundary conditions
       if(strcmp(key, "nBC") == 0){
           Model->nBCs = atoi(value);
           if(Model->nBCs != 0){
               Model->nameBCs = (char **) malloc(sizeof(char *) * Model->nBCs);
               Model->typeBCs = (int *)   malloc(sizeof(int)    * Model->nBCs);
               Model->paraBCs = (real **) malloc(sizeof(real *) * Model->nBCs);
           }
           else 
           {
              Model->nameBCs = NULL;
              Model->typeBCs = NULL;
              Model->paraBCs = NULL;
           }

           //store data for each boundary
           Model->BCFunSpf   = xfe_False;
           for(i=0; i<Model->nBCs; i++){
               Model->nameBCs[i] = (char *)malloc(sizeof(char) * 30);
               fscanf(fp, "%s%s%d", Model->nameBCs[i], value, &tmp);
               
               if(strcmp(value, "fullBC") == 0)              Model->typeBCs[i] = fullBC;
               else if(strcmp(value, "noneBC") == 0)         Model->typeBCs[i] = noneBC;
               else if(strcmp(value, "subinflowBC") == 0)    Model->typeBCs[i] = subinflowBC;
               else if(strcmp(value, "staginflowBC") == 0)   Model->typeBCs[i] = staginflowBC;
               else if(strcmp(value, "staticPBC") == 0)      Model->typeBCs[i] = staticPBC;
               else if(strcmp(value, "noslipwallBC") == 0)   Model->typeBCs[i] = noslipwallBC;
               else if(strcmp(value, "isothermwallBC") == 0) Model->typeBCs[i] = isothermwallBC;
               else if(strcmp(value, "slipwallBC") == 0)     Model->typeBCs[i] = slipwallBC;
               else if(strcmp(value, "farPBC") == 0)         Model->typeBCs[i] = farPBC;
               else if(strcmp(value, "exactBC") == 0)        Model->typeBCs[i] = exactBC;
               else
               { xf_printf("boundary condition is not defined\n");   return xf_Error(xf_NOT_SUPPORTED); }
           
               if(tmp > 0 && tmp != 999) {
                   //Model->BCFunSpf   = xfe_False;
                   Model->paraBCs[i] = (real *) malloc(sizeof(real) * tmp);
                   for(j=0; j<tmp; j++)
                      fscanf(fp, "%s%lf", value, &Model->paraBCs[i][j]);
                }
               else if(tmp == 999)   //user defined boundary condition in code
               {
                  Model->BCFunSpf   = xfe_True;
                  ierr = xf_Error(init_b());
                  if (ierr != xf_OK) return ierr;
                   
                  Model->paraBCs[i] = NULL;
               }
               else
                   Model->paraBCs[i] = NULL;

            }
       
          //jump one line
          fgets(dummy, xf_MAXLINELEN, fp);
          
           xf_printf("Boundary condition read!\n");
           continue;
       } 

       if(strcmp(key, "nOutput") == 0){
       
          Model->nOutput = atoi(value);
if(Model->nOutput > 0){ 
          //allocation 
          ierr = xf_Error(xf_Alloc((void **) &Model->Output, Model->nOutput, sizeof(Yu_Output)));
          if (ierr != xf_OK) return ierr;

          for(i=0; i<Model->nOutput; i++)
          {
         
             //initialize
             ierr = xf_Error(Yu_InitOutput(&Model->Output[i]));
             if (ierr != xf_OK) return ierr;

             //read necessary specification
             while(1){
             fgets(outputdummy, xf_MAXLINELEN, fp);
      
             ierr = xf_Error(xf_ReadKey(outputdummy, "=", outputkey, outputvalue, xf_MAXSTRLEN));
             if (ierr != xf_OK) return ierr;

             if(strcmp(outputvalue, "EndOutputDefine") == 0)
                break;
             else
             {
                if(strcmp(outputkey, "OutputName") == 0)
                {   sprintf(Model->Output[i].Name, "%s", outputvalue); 
                   continue;}

                if(strcmp(outputkey, "OutputType") == 0)
                {   
                   if(strcmp(outputvalue, "DomainIntegral") == 0)
                      Model->Output[i].Type = xfe_DomainIntegral;
                   else if(strcmp(outputvalue, "BoundaryIntegral") == 0)
                   {
                      Model->Output[i].Type = xfe_BoundaryIntegral;
                      Model->Output[i].UsesFlux = xfe_True;  //only support the flux associated output
                   }
                   else if(strcmp(outputvalue, "PointValue") == 0)
                      Model->Output[i].Type = xfe_PointValue;
                   else
                      return xf_NOT_SUPPORTED;
                   continue;
                }

                if(strcmp(outputkey, "nOutQuantities") == 0)
                {
                   Model->Output[i].nVars = atoi(outputvalue);
                   ierr = xf_Error(xf_Alloc((void **) &Model->Output[i].IVars, Model->Output[i].nVars, sizeof(int)));
                   if (ierr != xf_OK) return ierr;
                   //read all quantities specification
                   for(j=0; j<Model->Output[i].nVars; j++)
                   {
                      fscanf(fp, "%s", value);
                      for(k=0; k<OutputQuantity_last; k++)
                         if(strcmp(value, OutputQuanityName[k]) == 0)
                         {Model->Output[i].IVars[j] = k; break;}
                      if(k == OutputQuantity_last) 
                         return xf_NOT_SUPPORTED; //quan is not defined
                      
                   }//j
                      //jump for a line
                      fgets(outputdummy, xf_MAXLINELEN, fp);
                      continue;
                }

             
                //some parameters for boundary integral output
                if(strcmp(outputkey, "BoundaryName") == 0)
                {  sprintf(Model->Output[i].BFGTitles, "%s", outputvalue); continue;}

                if(strcmp(outputkey, "nFluxComponent") == 0)
                {  Model->Output[i].nFluxComponent = atoi(outputvalue); 
                   ierr = xf_Error(xf_Alloc((void **) &Model->Output[i].FluxComponentRanks, 
                                   Model->Output[i].nFluxComponent, sizeof(int)));
                   if (ierr != xf_OK) return ierr;
                   ierr = xf_Error(xf_Alloc((void **) &Model->Output[i].FluxComponentWeights, 
                                   Model->Output[i].nFluxComponent, sizeof(real)));
                   if (ierr != xf_OK) return ierr;

                   continue;
                }

                if(strcmp(outputkey, "nFluxRanks") == 0)
                {
                   for(j=0; j<Model->Output[i].nFluxComponent; j++)
                      sscanf(outputvalue, "%d", &Model->Output[i].FluxComponentRanks[j]);
                   printf("%d", Model->Output[0].FluxComponentRanks[0]);
                   continue;
                }

                if(strcmp(outputkey, "nFluxWeights") == 0)
                {
                   for(j=0; j<Model->Output[i].nFluxComponent; j++)
                      sscanf(outputvalue, "%lf", &Model->Output[i].FluxComponentWeights[j]);
                   printf("%lf", Model->Output[0].FluxComponentWeights[0]);
                   continue;
                }

                if(strcmp(outputkey, "RestartPointStat") == 0)
                {  Model->Output[i].RestartPointStat = TrueOrFalse(outputvalue);
                   continue;}

                if(strcmp(outputkey, "OutputSampleInterval") == 0)
                {  Model->Output[i].SampleTimeInv = atof(outputvalue);
                   continue;}

                if(strcmp(outputkey, "OutputFileWriteFreq") == 0)
                {  Model->Output[i].OutputFile_Freq_Ratio_DataFile = atoi(outputvalue);
                   continue;}

                if(strcmp(outputkey, "SequanceDump") == 0)
                {  Model->Output[i].SequanceDump = TrueOrFalse(outputvalue);
                   continue;}

                if(strcmp(outputkey, "WriteOffset") == 0)
                {  Model->Output[i].File_Write_Offset = atoi(outputvalue);
                   continue;}
             }
             }

          }//i

          //jump one line
          fgets(dummy, xf_MAXLINELEN, fp);
          
          xf_printf("Output specified!\n");
}
else
{
   Model->Output = NULL;
}
          continue;
       }
       //if we run mesh adaptation
       if(strcmp(key, "MeshAdaptation") == 0){
          Model->Stat_h_Adapt = TrueOrFalse(value);
          continue;
       }
       if(strcmp(key, "DynamicPAdaptation") == 0){
          Model->Dyn_p_Adapt = TrueOrFalse(value);
          
          if(Model->Dyn_p_Adapt)//initialize the param struct
          {
             xf_Error(xf_Alloc((void **) &Model->Dyn_p_Adapt_param, 1, sizeof(Yu_Dyn_p_Adapt_param)));
             if (ierr != xf_OK) return ierr;
             
             //read in parameters
             while(1){
                fgets(outputdummy, xf_MAXLINELEN, fp);
             
                ierr = xf_Error(xf_ReadKey(outputdummy, "=", outputkey, outputvalue, xf_MAXSTRLEN));
                if (ierr != xf_OK) return ierr;
             
                if(strcmp(outputvalue, "EndDynamicPAdaptation") == 0)
                   break;
                else
                {
                   if(strcmp(outputkey, "MaxOrder") == 0)
                      Model->Dyn_p_Adapt_param->MaxOrder = atoi(outputvalue);
                   if(strcmp(outputkey, "MinOrder") == 0)
                      Model->Dyn_p_Adapt_param->MinOrder = atoi(outputvalue);
                   if(strcmp(outputkey, "spongeOrder") == 0)
                      Model->Dyn_p_Adapt_param->spongeOrder = atoi(outputvalue);
                   if(strcmp(outputkey, "refinefrac") == 0)
                      Model->Dyn_p_Adapt_param->refinefrac = atof(outputvalue);
                   if(strcmp(outputkey, "coarsenfrac") == 0)
                      Model->Dyn_p_Adapt_param->coarsenfrac = atof(outputvalue);
                   if(strcmp(outputkey, "Yu_adapt_time_size") == 0)
                      Model->Dyn_p_Adapt_param->Yu_adapt_time_size = atof(outputvalue);
                   if(strcmp(outputkey, "xfa_file_write_inv") == 0)
                      Model->Dyn_p_Adapt_param->xfa_file_write_inv = atof(outputvalue);
                   if(strcmp(outputkey, "Yu_load_balance_time_size") == 0)
                      Model->Dyn_p_Adapt_param->Yu_load_balance_time_size = atof(outputvalue);
                }
          }
             //jump one line
             fgets(dummy, xf_MAXLINELEN, fp);
          }//
          continue;
       }


   }while(feof(fp) == 0);
 
   //consistency check for AV model
   if(Model->AVmodel && ! Model->EntropyBdFlag)
   { xf_printf("AV always goes with entropy bounding!\n"); return xf_NOT_SUPPORTED;}

   //for negative pressure checking
   //for legacy purpose; keep them there
   Model->Num_negPckpnt = 0;
   Model->Coord_negPckpnt = NULL;
   Model->Phi_negPckpnt = NULL;

   //close model.yu
   fclose(fp);

   if(Model->DetailChem == xfe_True)
   {
      //initialize CHEMKIN if necessary
      //distabled for now
      /*
      if(CHEMKIN_Init()!=0)
         printf("Error in CK initialization!~\n");
      if(Chemkin_DoubleFlux()!= 0)
         printf("Error in Introducing Thermal Model!~\n");
   
      CK_Basic_Info(&Num_elem, &Num_spe, &Num_reac);
      xf_printf("+++++++++CHEMKIN is initialized+++++++++\n");
      xf_printf("Include: %d elements; %d species; %d reactions\n", Num_elem, Num_spe, Num_reac);
      
      //do not forget
      CK_Moleweight_Species(Model->moleW);
      */
   }
    else
    {
        for(i=1; i<(Model->nVars - Model->dim - 2)+1; i++)
            Model->moleW[i] = Model->moleW[0];
    }

   
    line = NULL;
   return xf_OK;
}

/******************************************************************/
// free the allocation for data attached to Model
int
DestroyModel(Yu_Model *Model)
{
   int i, j;

   if(Model == NULL) return xf_OK;

   //for negative pressure checking
   if(Model->Coord_negPckpnt != NULL)
      xf_Release(Model->Coord_negPckpnt);
   if(Model->Phi_negPckpnt != NULL)
      xf_Release(Model->Phi_negPckpnt);

   //free variable name
   for(i=0; i<Model->nVars; i++)
      if(Model->nameVars != NULL)
      xf_Release(Model->nameVars[i]);
   if(Model->nameVars != NULL)
   xf_Release(Model->nameVars);

   //free initial state
   if(Model->initVars != NULL)
   xf_Release(Model->initVars);

   //free boundary name
   for(i=0; i<Model->nBCs; i++)
      if(Model->nameBCs[i] != NULL)
      xf_Release(Model->nameBCs[i]);
   if(Model->nameBCs != NULL)
   xf_Release(Model->nameBCs);

   //free boundary type
   if(Model->typeBCs != NULL)
      xf_Release(Model->typeBCs);

   //free boundary parameters
   for(i=0; i<Model->nBCs; i++)
   if(Model->paraBCs[i] != NULL)
      free(Model->paraBCs[i]);
   if(Model->paraBCs != NULL)
      free(Model->paraBCs);

   //free CHEMKIN allocation if necessary
   //if(Model->DetailChem == xfe_True)
   //   CHEMKIN_Free();

   //if entropy bound is used
   if(Model->EntropyBdFlag){
      DestroyEntropyBoundStruct();
   }

   //set to NULL
   Model = NULL;

   return xf_OK;
}

/******************************************************************/
//provide CFL number according to Yu's talented idea
//currently not able to handle geomtric info of each cell
static int
CFLnumber(const int Order, enum xfe_BasisType Basis, real *CFL)
{
   int ierr;
   enum xfe_ShapeType Shape;


   ierr = xf_Error(xf_Basis2Shape(Basis, &Shape));
   if(ierr != xf_OK) return ierr;

   switch(Shape){
      case xfe_Segment:
         if(Order == 0)
            (*CFL) = 1.0;
         else if(Order == 1)
            (*CFL) = 0.5;
         else if(Order == 2)
            (*CFL) = 0.167;
         else if(Order == 3)
            (*CFL) = 0.123;
         else if(Order == 4)
            (*CFL) = 0.073;
         else
            return xf_NOT_SUPPORTED;
         break;
     case xfe_Triangle:
         if(Order == 0)
            (*CFL) = 1.0;
         else if(Order == 1)
            (*CFL) = 0.135;
         else if(Order == 2)
            (*CFL) = 0.067;
         else if(Order == 3)
            (*CFL) = 0.058;
         else if(Order == 4)
            (*CFL) = 0.033;
         else
            return xf_NOT_SUPPORTED;
         break;
     case xfe_Quadrilateral:
         if(Order == 0)
            (*CFL) = 1.0;
         else if(Order == 1)
            (*CFL) = 0.25;
         else if(Order == 2)
            (*CFL) = 0.083;
         else if(Order == 3)
            (*CFL) = 0.062;
         else if(Order == 4)
            (*CFL) = 0.036;
         else
            return xf_NOT_SUPPORTED;
         break;
     case xfe_Tetrahedron:
         if(Order == 0)
            (*CFL) = 1.0;
         else if(Order == 1)
            (*CFL) = 0.066;
         else if(Order == 2)
            (*CFL) = 0.035;
         else if(Order == 3)
            (*CFL) = 0.015;
         else if(Order == 4)
            (*CFL) = 0.013;
         else
            return xf_NOT_SUPPORTED;
         break;
     case xfe_Hexahedron:
         if(Order == 0)
            (*CFL) = 1.0;
         else if(Order == 1)
            (*CFL) = 0.167;
         else if(Order == 2)
            (*CFL) = 0.056;
         else if(Order == 3)
            (*CFL) = 0.041;
         else if(Order == 4)
            (*CFL) = 0.024;
         else
            return xf_NOT_SUPPORTED;
         break;

      default:
         return xf_Error(xf_UNKNOWN_SHAPE);
         break;
   }

   return xf_OK;
}

/******************************************************************/
//estimate time step with
int
Yu_EstimateTimeStep(xf_All *All, Yu_Model *Model, xf_Vector *State, real *timestep)
{
    int iegrp, ielem, ierr, Order;
    real tmp, CFL, fac, advDt, diffDt;
    real MaxCharSpeed_elem, MinFaceLen_elem;
    xf_Vector *MaxCharSpeed, *MinFaceLen;
    enum xfe_BasisType Basis;
    xf_Mesh *Mesh;

    Mesh = All->Mesh;
    
    //catch data pointer for Model
    MaxCharSpeed = Model->MaxCharSpeed;
    MinFaceLen   = Model->MinFaceLen;
    
    //loop over each element to compute the time step
    tmp = 1.0e+16;
    fac = 1.0;
    for (iegrp=0; iegrp<Mesh->nElemGroup; iegrp++){
        
        Basis = State->Basis[iegrp];
 
        for (ielem=0; ielem<Mesh->ElemGroup[iegrp].nElem; ielem++){
        
            //need to know the basis for CFL number
            Order = xf_InterpOrder(State, iegrp, ielem); 
            //ierr = xf_Error(CFLnumber(Model->order, Basis, &CFL));
            ierr = xf_Error(CFLnumber(Order, Basis, &CFL));
            if(ierr != xf_OK) return ierr;
            
            MaxCharSpeed_elem = MaxCharSpeed->GenArray[iegrp].rValue[ielem][0];
            MinFaceLen_elem   = MinFaceLen->GenArray[iegrp].rValue[ielem][0];
           

            advDt = MinFaceLen_elem / MaxCharSpeed_elem * CFL  * 0.8;
            diffDt = 1.e+16;
            if(Model->DiffFlag){
               // we need consider the influence of viscous part
               //factor for diffusion time step criterion
               if(Model->order <=1)
                  fac = 2.8;
               else if(Model->order == 2)
                  fac = 0.86;
               else if(Model->order == 3)
                  fac = 0.4;
               else if(Model->order == 4)
                  fac = 0.24;
               else if(Model->order == 5)
                  fac = 0.16;
               else
                  return xf_NOT_SUPPORTED;

               //test with scalar adv-diff
               if(LinearCase)
                  diffDt = fac*pow(MinFaceLen_elem, 2.0)/LapVis/pow((real)2*Model->order+1,2.);
               else
               {
                  diffDt = fac*pow(MinFaceLen_elem, 2.0)/Model->mu_c/pow((real)2*Model->order+1,2.);
                  //xf_printf("The diffustion CFL is not implemented for NS type eqn!\n");
                  //return xf_NOT_SUPPORTED;
               }
           
               //if(diffDt< tmp)
               //   tmp = diffDt;
            }

            if(Model->DiffFlag && !Model->AVmodel)
            advDt = 1.0 * pow((pow(advDt, -1.) + pow(diffDt, -1.)),-1);
            else
               advDt *= 0.7;

            if(advDt < tmp)
               tmp = advDt; 

            if(tmp < 1.e-8)
                xf_printf("Warning, time step is too close to zero!\n");
        }
    }
    
    
    //return value
    (*timestep) = Model->time_fac * tmp;

    //for parallel implementation
    ierr = xf_Error(xf_MPI_Allreduce(timestep, 1, xfe_SizeReal, xfe_MPI_MIN));
    if (ierr != xf_OK) return ierr;

    //set in to structure
    Model->dt_size = (*timestep);

    //delete the info of previous step
    ierr = xf_Error(xf_SetConstVector(MaxCharSpeed, 0, 0.0));
    if(ierr != xf_OK) return ierr;

   return xf_OK;
}

/******************************************************************/
// Model based quadrature order requirement
int
ElemQuadOrder(int Order, int *QuadOrder)
{ 
   (*QuadOrder) = 2*Order + 1;

   return xf_OK;
}

/******************************************************************/
// see above
int 
FaceQuadOrder(int Order, int *QuadOrder)
{
   (*QuadOrder) = 2*Order + 1;

   return xf_OK;
}

/******************************************************************/
// Evaluate Sponge Source Term
int
SpongeSource(int nq, int sr, int dim, real *xglob, real *U, real *S, 
             real gamma, real dt, real tcurr, real SpongeMeasure, Yu_Sponge *SP)
{
   int ierr, i, j, tmp;
   real *Un, *Sn, *xn, pt[2], r;
   real Inf[10], omg;
   real rho, u, v, w, p, perbE, omega;

   //load in far field state variables
   if(dim == 2)
   {
      Inf[0] = SP->State_Inf[0];
      Inf[1] = SP->State_Inf[0] * SP->State_Inf[1];
      Inf[2] = SP->State_Inf[0] * SP->State_Inf[2];
      Inf[3] = (SP->State_Inf[3])/(gamma-1.) + 0.5 * (Inf[1]*Inf[1]+Inf[2]*Inf[2])/Inf[0];
   
   }

   if(dim == 3){
      Inf[0] = SP->State_Inf[0];
      Inf[1] = SP->State_Inf[0] * SP->State_Inf[1];
      Inf[2] = SP->State_Inf[0] * SP->State_Inf[2];
      Inf[3] = SP->State_Inf[0] * SP->State_Inf[3];
      Inf[4] = (SP->State_Inf[4])/(gamma-1.) + 0.5 * (Inf[1]*Inf[1]+Inf[2]*Inf[2]+Inf[3]*Inf[3])/Inf[0];
   }

   //damping factor
   //omg = sqrt(1.4) * 5. / 10.;
   //center[0] = 5.0; center[1] = 5.0;
   for(i=0; i<nq; i++)
   {
      Un = U + i*sr;
      Sn = S + i*sr;
      xn = xglob + i*dim;

      //r = sqrt( (xn[0]-center[0]) * (xn[0]-center[0]) + (xn[1]-center[1])*(xn[1]-center[1]) ) - 5.0;
      //r = fabs(xn[0] - center[0]) - 5.0; 
      //if(r>0. && SpongeMeasure > 0.){
      for(j=0; j<sr; j++)
         Sn[j] = 0.;
      
      if(SpongeMeasure > 0.){
     
         pt[0] = xn[0]; pt[1] = xn[1];
         if(dim == 3)
            pt[1] = sqrt(xn[1]*xn[1] + xn[2]*xn[2]);
     
    //    if(SP->WheInflowSponge)
    //    {
           //specify inflow sponge 
    //       r = fabs(pt[0]);
    //    }
    //    else
    //  if we need to specify a mult-state sponge; it has to be done in the code; here

        {
           //this is a zonal based sponge term specification
           //         |                        |
           //  1001   |            1           |   11
           //---------|------------------------|---------
           //         |                        |
           //  1000   |                        |   10
           //         |                        |
           //---------|------------------------|---------
           //  1100   |          100           |  110
           //         |                        |  
         switch((int) SpongeMeasure){
            case 1:
               //beyond frist boundary
               r = fabs(SP->line[0*3+0]*pt[0]+SP->line[0*3+1]*pt[1]+SP->line[0*3+2])
                   /sqrt(SP->line[0*3+0]*SP->line[0*3+0]+SP->line[0*3+1]*SP->line[0*3+1]); 
               break;

            case 10:
               //second
               r = fabs(SP->line[1*3+0]*pt[0]+SP->line[1*3+1]*pt[1]+SP->line[1*3+2])
                   /sqrt(SP->line[1*3+0]*SP->line[1*3+0]+SP->line[1*3+1]*SP->line[1*3+1]); 
               break;

            case 100:
               //third
               r = fabs(SP->line[2*3+0]*pt[0]+SP->line[2*3+1]*pt[1]+SP->line[2*3+2])
                   /sqrt(SP->line[2*3+0]*SP->line[2*3+0]+SP->line[2*3+1]*SP->line[2*3+1]); 
               break;

            case 1000:
               //forth
               r = fabs(SP->line[3*3+0]*pt[0]+SP->line[3*3+1]*pt[1]+SP->line[3*3+2])
                   /sqrt(SP->line[3*3+0]*SP->line[3*3+0]+SP->line[3*3+1]*SP->line[3*3+1]); 
               break;

            case 11:
               //to intersect first and second lines (point 2)
               r = sqrt((pt[0]-SP->point[1*3+0])*(pt[0]-SP->point[1*3+0]) 
                   + (pt[1]-SP->point[1*3+1])*(pt[1]-SP->point[1*3+1]));
               break;

            case 110:
               //to intersect second and third (point 3)
               r = sqrt((pt[0]-SP->point[2*3+0])*(pt[0]-SP->point[2*3+0]) 
                   + (pt[1]-SP->point[2*3+1])*(pt[1]-SP->point[2*3+1]));
               break;

            case 1100:
               //to intersect third and forth (point 4)
               r = sqrt((pt[0]-SP->point[3*3+0])*(pt[0]-SP->point[3*3+0]) 
                   + (pt[1]-SP->point[3*3+1])*(pt[1]-SP->point[3*3+1]));
               break;

            case 1001:
               //to intersect forth and first (point 1)
               r = sqrt((pt[0]-SP->point[0*3+0])*(pt[0]-SP->point[0*3+0]) 
                   + (pt[1]-SP->point[0*3+1])*(pt[1]-SP->point[0*3+1]));
               break;

         }
        }
          
        omg = 20.0 * r * r;
          
        //can be specified individually for different zones
        if(SP->WheMultiStateSponge)
        {
            tmp = (int) SpongeMeasure;
            if(tmp == 1001 || tmp == 1000 || tmp == 1100){
                  
                u = 3./2.0 + 1./2.0 * tanh(2.0 * xn[1] / 1.0);
                v = 0.0;
                rho = 1.0;
                p = 11.4285714286;
                
                omg = 2.5/10./10. * r * r;
            }
              
            if(tmp == 1 || tmp == 11){
                u = 2.; v = 0.; rho = 1.; p = 11.4285714286;
                omg = 2.5/20./20. * r * r;
            }
              
            if(tmp == 100 || tmp == 110){
                u = 1.; v = 0.; rho = 1.; p = 11.4285714286;
                omg = 2.5/20./20. * r * r;
            }
              
            if(tmp == 10){
                u = 3./2.0 + 1./2.0 * tanh(2.0 * xn[1] / 40.0);
                v = 0.; rho = 1.; p = 11.4285714286;
                omg = 2.5/800./800. * r * r;
            }
              
      
            //perturbation can be added out of sponge
            perbE = 0.75 * exp(-log(2.) * ((xn[0]-20.)*(xn[0]-20.) + xn[1]*xn[1])/0.16/0.16);
            omega = 2. * M_PI * 0.05;
            r = 0.001*sin(omega * tcurr) + 0.0005*sin(omega*tcurr + M_PI/2.);
            u += perbE * (xn[1])/0.16 * r; 
            v -= perbE * (xn[0]-20.)/0.16 * r;
            
            Inf[0] = rho; Inf[1] = rho*u; Inf[2] = rho*v; Inf[3] = p/(1.4-1.) + 0.5*rho*(u*u+v*v);
        }


        //can be hidden for some purposes
      //for(j=0; j<dim+2; j++)
      //   Sn[j] = (1. - exp(-omg*dt))*(Un[j] - Inf[j]);
  
      }

   }
   
   return xf_OK;
}
/******************************************************************/
// Evaluate Chemical Source Term
int
ChemicalSource(const int nq, const int sr, const int dim, 
               const real *U, real *xglob, real *S, const real Gamma,
               enum xfe_Bool DetailChem, const real dt_size)
{
   int ierr, i, j, k;
   int Nspe;
   real *Un, *Sn, *xn;
   real rho, u, v, w, tmp_p, p, Y[50], T, sum;
   real R, E_R, M, Q, A;


   if(DetailChem == xfe_True)
   {
      //number of species
      Nspe = sr - 4 + 1;

      for(i=0; i<nq; i++)
      {
         Un = U + i*sr;
         Sn = S + i*sr;
         
         rho = Un[0];
         u   = Un[1]/Un[0];
         v   = Un[2]/Un[0];
         p   = (Gamma - 1.0)*(Un[3] - 0.5*rho*(u*u + v*v));
      
         if(sr > 4)
         {
            sum = 0.0;
            for(k=4; k<sr; k++)
            {
              //if out of range, just ignor it
              //lower bound
              Y[k-4] = Un[k]/Un[0];
              if(Y[k-4] < MEPS)
                 Y[k-4] = MEPS;

              //higher bound
              if(Y[k-4] > 1.0)
                 Y[k-4]= 1.0;

              sum += Y[k-4];
            }
            Y[sr-4] = 1.0 - sum;
         }
         else
         {
            xf_printf("Scalar is not specified; model fails!");
            return xf_OUT_OF_BOUNDS;
         }

         //call stiff solver
         tmp_p = p;
         //one-step chemistry model
         //gamma = 1.24; Q_RT0 = 26.57; E_RT0 = 40; A = 2703.3;
         //stiffsolver(Y, Nspe, &p, rho, Gamma, dt_size);

         //figure out source term
         Sn[0]  = 0.0;
         Sn[1]  = 0.0;
         Sn[2]  = 0.0;
     
         //be careful with the sign in front
         //for(k=4; k<sr; k++)
         //   Sn[k] = 0.0;

         //reactive scalar 
         Sn[4] = - dt_size * rho * 2703.3 * (1.0-Y[0]) * exp(-40.0/ (p/rho));

         //total non-chemical energy
         Sn[3] = 26.5714 * Sn[4];
         //Sn[3]  = -(p - tmp_p)/(Gamma - 1.0);
         //for(k=4; k<sr; k++)
         //   Sn[k] = -(Y[k-4] - Un[k]);

      }
   }
   else
   {
      //Yu's one step detonation model
      //H2/O2/Ar = 2/1/7 at T = 293 and Pa = 26.7kPa
      //E_R = 1.012969108494726e+04;
      //M = 0.001 * 27.97;
      //R = 8.31446/M;
      //Q = 4.367523837442915e+03 * R;
      //A = 1.82e+08;
   
      for(i=0; i<nq; i++)
      {
         Un = U + i*sr;
         Sn = S + i*sr;
         xn = xglob + i*dim;
   
         rho = Un[0];
         u   = Un[1]/Un[0];
  /*       v   = Un[2]/Un[0];
         p   = (Gamma - 1.0)*(Un[3] - 0.5*rho*(u*u + v*v));
         T   = p/rho/R;
         if(sr > 4)
            Y[0] = Un[4]/Un[0];
         else
         {
            xf_printf("Scalar is not specified; model fails!");
            return xf_OUT_OF_BOUNDS;
         }
         
         //specify source terms
         Sn[0]  = 0.0;
         Sn[1]  = 0.0;
         Sn[2]  = 0.0;
         if(Y[0] > 1.0)
         { 
            Sn[3] = 0.0;
            Sn[4] = 0.0;
         }
         else
         {
            //Source contribute to Residuel on the left side of the eqn
            Sn[3]  = -Q*rho*A*(1.0 - Y[0])*exp(-E_R/T);
            Sn[4]  = -rho*A*(1.0 - Y[0])*exp(-E_R/T);
         }
   
         for(j=5; j<sr; j++)
            Sn[5] = 0.0;
            */
         for(j=0; j<sr; j++)
            Sn[j] = 0.;
         
         //below is for channel flow
         //channel flow source on x-axis momentum
         Sn[1] = -1.0;
         //channel flow source on energy equation
         Sn[4] = -1. * u; 

         //below is for cylinder case shedding excitation
         //if(fabs(xn[0])<=2. && fabs(xn[1])<=2.)
         //{
         //   Sn[1] = 0.05 * sin(4. * M_PI * xn[0]) * cos(2. * M_PI * xn[1]);
         //   Sn[2] = 0.05 * cos(4. * M_PI * xn[0]) * sin(2. * M_PI * xn[1]);
         //}
      }
   }

   return xf_OK;
}



/******************************************************************/
int
Yu_GhostStateConstruct(const int nq, const int sr, const int dim, const real *quL, const real *quR, 
                       real *quLghost, real *quRghost, const real GammaL, const real GammaR)
{
   int ierr, i, j, sr2;
   const real *UL, *UR;
   real *ULghost, *URghost;
   real pl, pr, rhol, rhor, ul, ur;

   sr2 = sr*sr;

   for(i=0; i<nq; i++){
   
      UL  = quL + i*sr;
      UR  = quR + i*sr;
      ULghost = quLghost + i*sr;
      URghost = quRghost + i*sr;

      for(j=0; j<sr; j++)
      {
         ULghost[j] = UL[j];
         URghost[j] = UR[j];
      }

      if(dim == 2){
      pl = (GammaL - 1.0) * (UL[3] - 0.5*(UL[1]*UL[1] + UL[2]*UL[2])/UL[0]);
      ULghost[3] = pl/(GammaR - 1.0) + 0.5*(UL[1]*UL[1] + UL[2]*UL[2])/UL[0];
      pr = (GammaR - 1.0) * (UR[3] - 0.5*(UR[1]*UR[1] + UR[2]*UR[2])/UR[0]);
      URghost[3] = pr/(GammaL - 1.0) + 0.5*(UR[1]*UR[1] + UR[2]*UR[2])/UR[0];
      }
      else 
         return xf_Error(xf_OUT_OF_BOUNDS);
   }

   return xf_OK;
}
/******************************************************************/
// Evaluate Boundary State for Slip-Wall or Symmetrical BC
static int
SlipWallBC(const real *UI, const real*n, real *UB, const int Num_Vars, const int dim)
{
       int d, k;
       double rVn;
               
       for (d=0, rVn=0.0; d<dim; d++) rVn += UI[d+1]*n[d];
                   
       for (k=0; k<Num_Vars; k++) UB[k] = UI[k];
         
       for (d=0; d<dim; d++) UB[d+1] -= n[d]*rVn;
       
       return xf_OK;
}
/******************************************************************/
// Evaluate Boundary State for No-slip wall (make sense for NS)
static int
NoSlipWallBC(const real *UI, const real *n, real *UB, const int Num_Vars, const int dim)
{
   int d, k;

   for (k=0; k<Num_Vars; k++) UB[k] = UI[k];

   for(d=0; d<dim; d++) UB[d+1] = 0.0;

   return xf_OK;
}
/******************************************************************/
// Evaluate Boundary State for No-slip wall (make sense for NS)
static int
NoSlipIsoThermalWallBC(const real *UI, const real *n, real *UB, const int Num_Vars, const int dim, 
                       const real gimoR, const real *localBCpara)
{
   int d, k;
   real wallT;

   wallT = localBCpara[0];

   for (k=0; k<Num_Vars; k++) UB[k] = UI[k];

   //set noslip condition
   for(d=0; d<dim; d++) UB[d+1] = 0.0;

   //conserve total energy
   UB[dim+1] = UI[dim+1];

   //density change according to the given wall temperature
   UB[0] = UI[dim+1] * gimoR/ wallT; 

   return xf_OK;
}
/******************************************************************/
// Evaluate Boundary State for Characteristic non-reflecting Boundary
// Input from model.yu:
//   localBCpara[i] = rInf, uInf[dim], pInf, Yinf[sr-dim-2]
static int
FarFieldPressureBC(const int dim, const int sr, const real *x, const real *UI, const real *localBCpara, 
                   const real *n, real *UB, const real gamma)
{
   int i, j, k;
   real cI, rI, uI, pI, VnI, rVI2, JpI;
   real Vninf, cinf, Jminf, VnB, cB;
   real S, VB2, v;
   //infinity far state
   real pInf, rInf, uInf[3], YInf[50];

   //pull in the parameter for BC input
   rInf = localBCpara[0];
   
   for(j=0; j<dim; j++) uInf[j] = localBCpara[1+j];
   
   //for test !!!!
   //intro a velocity profile
   //uInf[1]= 0.0;
   //uInf[0]=  0.2 + 0.1*tanh((x[1]-0.5) / 0.2);
   
   pInf = localBCpara[dim+1];
   for(j=dim+2; j<sr; j++) YInf[j-dim-2] = localBCpara[j];

   // calculate interior normal velocity, VI dot n
   rI = UI[0];
   for (i=0, VnI=0.; i<dim; i++)  VnI += UI[i+1]/rI*n[i];

   // compute interior pressure
   for (i=0, rVI2=0.0; i<dim; i++) rVI2 += UI[i+1]*UI[i+1]/rI;
   pI = (gamma-1.)*(UI[1+dim] - 0.5*rVI2);

   // computer interior speed of sound
   cI = sqrt(gamma*pI/rI);

   //computer interior Rieman invariant, JpI = VnI + 2*cI/(gamma-1.)
   JpI = VnI + 2.*cI/(gamma-1.);

   // uInf dot n
   for (i=0, Vninf=0.; i<dim; i++) Vninf += uInf[i]*n[i];

   //speed of sound for infinity state
   cinf = sqrt(gamma*pInf/rInf);

   // J-Riemann invariant for infinity state
   Jminf = Vninf - 2.*cinf/(gamma-1.);

   // normal velocity for boundary state
   VnB = 0.5*(JpI + Jminf);

   // speed of sound for boundary state
   cB = (gamma-1.)*0.25*(JpI - Jminf);

   if(VnI < 0.0) { //inflow
   
      //entropy for boundary state, taken from infinity state
      S = pInf/pow(rInf, gamma);

      //density for boundary state
      UB[0] = pow(cB*cB/(S*gamma), 1./(gamma-1.));

      //momentum for boundary state
      for (i=0, VB2=0.; i<dim; i++){
         v = uInf[i] + (VnB - Vninf)*n[i];
         UB[1+i] = UB[0]*v;
         VB2 += v*v;
      }

      //species density
      for(j=dim+2; j<sr; j++)
         //something weird with this BC
         //UB[j] = UB[0] * YInf[j];
         UB[j] = UI[j];
   }
   else { //outflow

      //entropy for boundary state, taken from interior state
      S = pI/pow(rI, gamma);

      //density for boundary state
      UB[0] = pow(cB*cB/(S*gamma), 1./(gamma-1.));

      //momentum for boundary state
      for (i=0, VB2=0.; i<dim; i++){
         v = UI[1+i]/rI + (VnB - VnI)*n[i];
         UB[1+i] = UB[0]*v;
         VB2 += v*v;
      }

      //species density
      for(j=dim+2; j<sr; j++)
         //something weird with this BC
         //UB[j] = UB[0] * UI[j]/rI;
         //something weird with this BC, too
         //UB[j] = UB[0] * YInf[j];
         UB[j] = UI[j];
   }

   //energy for boundary state
   UB[1+dim] = (1./gamma)*(1./(gamma-1.))*cB*cB*UB[0] + 0.5*UB[0]*VB2;

   //check far field BC
   //for(j=0; j<sr; j++)
   //   printf("%.12lf %.12lf\n", UI[j], UB[j]);
   //printf("pressure %.12lf %.12lf\n", pI, pInf);
   //getchar();

   return xf_OK;
}
/******************************************************************/
// Evaluate Boundary State for Subsonic Inflow using outgoing Characteristics
// given inflow velocity and density
static real A[6], St[6], phi[6], psi[6];
static real signA[6], signSt[6];
static int
init_b()
{
   int m, n;
   
   for(m=0; m<3; m++)
      for(n=0; n<2; n++)
      {
         A[m*2+n] = 0.01 + 0.01 * (real)(m*2+n);
         St[m*2+n] = 0.1 + 0.1 * (real)(m*2+n);
         signA[m*2+n] = 1.0;
         signSt[m*2+n] = 1.0;
         phi[m*2+n] = 0.0;
         psi[m*2+n] = 0.0;
      }

   return xf_OK;
}

int
Yu_UserDefinedBCparamsUpdate(const real t, const real dt_size)
{
   int myRank, ierr, chance, r, m, n;
   real dt_ref, tmp;

   ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
   if (ierr != xf_OK) return ierr;

if(myRank==0)
{
   dt_ref = 0.0085/1.613;
   tmp = dt_size/dt_ref/20.0;
   tmp = 1./tmp;
   
   chance = (int) tmp;
   srand(time(NULL));
   
   //update A
   r = rand() % chance;
   if(r==0)
      for(m=0; m<3; m++)
         for(n=0; n<2; n++)
         {
            if(A[m*2+n]<0.01 || A[m*2+n]>0.07)
               signA[m*2+n] *=-1;
            
            A[m*2+n] += signA[m*2+n]*0.0001;
         }

   //update St 
   r = rand() * chance;
   if(r==0)
      for(m=0; m<3; m++)
         for(n=0; n<2; n++)
         {
            if(St[m*2+n]<0.1 || St[m*2+n]>0.7)
               signSt[m*2+n] *=-1;
            
            St[m*2+n] += signSt[m*2+n]*0.00085;
         }

   //update phi
   r = rand() % chance;
   if(r==0)
      for(m=0; m<3; m++)
         for(n=0; n<2; n++)
            phi[m*2+n] += 0.00085;

   //update psi
   r = rand() % chance;
   if(r==0)
      for(m=0; m<3; m++)
         for(n=0; n<2; n++)
            psi[m*2+n] += 0.00085;
}
   //broadcast the data out
   ierr = xf_Error(xf_MPI_Bcast((void *) A, 6*sizeof(real), 0));
   if(ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_MPI_Bcast((void *) St, 6*sizeof(real), 0));
   if(ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_MPI_Bcast((void *) phi, 6*sizeof(real), 0));
   if(ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_MPI_Bcast((void *) psi, 6*sizeof(real), 0));
   if(ierr != xf_OK) return ierr;

   return xf_OK;
}
static int
compute_b(const real theta, const real t, real *b)
{
   int m, n;

   //assemble b
   (*b) = 12.5;
   for(m=0; m<3; m++)
      for(n=0; n<2; n++)
      {
         (*b) += A[m*2+n]*cos(St[m*2+n]*t/2.0 + phi[m*2+n])*cos(m*theta+psi[m*2+n]);
      }

   return xf_OK;
}
static int
SubsonicInflowBC(const int dim, const int sr, const real *UI, const real *uInput, 
                 const real ptInput, const real *YInput, const real *n, real *UB, const real gamma)
{
   int i, k, j;
   real dI, uI[3], pI, cI, dB, uB[3], pB, cB;
   real VnI, dV2I, VnB, dV2B, Ma, MaI, Jp, tmp;
   real a, b, c, Ht, sol1, sol2, sB;
   real TB, TI;

   //initialize UB
   for(k=0; k<sr; k++) UB[k] = UI[k];

   //calculate interior normal velocity, VI dot N
   dI = UI[0];
   for(i=0, VnI=0.0; i<dim; i++) 
      VnI += UI[1+i]/dI*n[i];

   //if (VnI > 0.)
   //   xf_printf("Warning: reversed flow on Inflow boundary.\n");

   //compute interior pressure and speed of sound
   for (i=0, dV2I=0.; i<dim; i++) dV2I += UI[1+i]*UI[1+i]/dI;
   pI = (gamma - 1.) * (UI[1+dim] - 0.5*dV2I);

   cI = sqrt(gamma*pI/dI);

   //computer interior Riemann invariant
   Jp = VnI + 2.*cI/(gamma - 1.);
   //Jp = - Jp;

   //computer interior total enthalpy
   //Ht = pI/dI*gamma/(gamma-1.) + 0.5*dV2I/dI;
   //provided
   //Ht = ptInput;

   //compute exterior normal velocity, uInput dot N and dV2B
   for(i=0, VnB=0.0, dV2B=0.0; i<dim; i++) {
      VnB += uInput[i]*n[i];
      dV2B += uInput[i]*uInput[i];
   }
   
   //a = 1. + 2./(gamma-1.);
   //b = 2.*Jp;
   //c = (gamma-1.)/2.*(Jp*Jp - 2.*Ht);
   //sol1 = -b/2./a + sqrt(b*b - 4.*a*c)/2./a;
   //sol2 = -b/2./a - sqrt(b*b - 4.*a*c)/2./a;


   //computer exterior speed of sound
   cB = (gamma - 1.)*(Jp - VnB)/2.;
   Ma = fabs(VnB) / cB;
   MaI = fabs(VnI) / cI; 
  // if(Ma > 1.0) return xf_NON_PHYSICAL;

   //!!!!~~~~~~~
   //new way by specifying the Mach number 
   // this way is numerically unstable
   //Ma = uInput[0];
   //cB = Jp / (-Ma + 2./(gamma-1.)) ;
/*  
   if(cB < 0.001)
   {
      //need a trick here
      cB = 0.001;
   }
*/
   //cB = sol1;
   //if (sol2 > cB)
   //   cB = sol2;
   //if(cB < 0.) 
   //{
   //   return xf_Error(xf_NON_PHYSICAL);
   //   getchar();
   //}

   //computer exterior Mach number
   //VnB = 2.*cB/(gamma-1.) + Jp;
   //Ma = uInput[0] / cB;
   //Ma = VnB / cB;
   //if(Ma > 0.) return xf_Error(xf_NON_PHYSICAL);

   //tmp = 1. + (gamma-1.)/2.*Ma*Ma;
   //pB = uInput[0] * pow(tmp, -gamma/(gamma-1.));
   //dB = uInput[1] * pow(tmp, -1./(gamma-1.));
   //tmp = pB/Ht*gamma/(gamma-1.);
   //dB = tmp * pow(tmp, -1./(gamma-1.));
   //printf("%lf %lf %lf %lf\n", uInput[0], uInput[1], pB, dB);
   //getchar();
   //back out exterior static pressure 
   //dB = ptInput/pow(1.+0.5*(gamma-1.)*Ma*Ma, 1./(gamma-1.));
   //dB = ptInput;
   //pB = cB*cB*dB/gamma;
   //try to set entropy
   sB = ptInput;
   //sB = 1.04406635586; 
   //sB = 1. / pow(dI, gamma - 1.);
   dB = pow((cB*cB/sB/gamma), 1./(gamma-1.));
   pB = cB*cB*dB/gamma; //ptInput;// * pow(tmp, -gamma/(gamma-1.)); 

   Ma = pow(1. + (gamma-1.)/2. * Ma * Ma, gamma/(gamma-1.));
   MaI = pow(1. + (gamma-1.)/2. * MaI * MaI, gamma/(gamma-1.));
   
   //pB = pI * MaI / Ma;
   pB = pI; 
   dB = pow(pB / sB, 1./gamma);
   //dB = 1.0;
   //pB = cB*cB*dB/gamma;

   //TI = pI / dI;
   //TB = cB * cB / gamma; 
   //prescrible pressure
   //dB = TB*dI/1.; 
   //pB = 1.0;
   //dB = pB/TB;
   //dB = pB/TB;
   //dB = dI * fabs((1. - TI) / ( TB - TI));
   //pB = cB*cB*dB/gamma; //ptInput;// * pow(tmp, -gamma/(gamma-1.)); 

   //printf("%lf %lf %lf %lf\n", TI, TB, dB, pB); getchar();
   //dB = dI * (TI - TB) / (1. - TB);
   //pB = TB / dB; 

   //printf("%lf %lf %lf %lf\n", dI, pI, dB, pB);
   //dB = 0.687044203472957;
   //dB = 1.0;
   //pB = cB*cB*dB/gamma;
 
 //uB[0] = 0.12 * cB;
   //uB[1] = 0.;
   //dV2B = uB[0] * uB[0];
   //back out exterior density
   //dB = gamma*pB/cB/cB;
 
   //Consolidate exterior state
   UB[0] = dB;
   //for(k=0; k<dim; k++)
   //{
   //   UB[1+k] = dB*uInput[k];
   //   printf("%lf\n", UB[1+k]);
   //}
   UB[1] = uInput[0] * dB;
   //UB[1] = Ma * cB * dB;
   UB[2] = uInput[1] * dB;
   //UB[2] = UI[2];
   if(dim == 3)
      UB[3] = uInput[2] * dB;

   UB[1+dim] = pB/(gamma - 1.) + 0.5*(UB[1]*UB[1] + UB[2] * UB[2])/UB[0];

   //for passive scalar
   for (i=2+dim; i<sr; i++)
      UB[i] = UB[0]*YInput[i-2-dim];

   return xf_OK;
}

/******************************************************************/
// Evaluate Boundary State for Subsonic Outflow with Fixed Static Pressure
static int
StaticPressureBC(const int dim, const int sr, const real *UI, const real p,
                 const real *n, real *UB, const real gamma)
{
   int i, k, j;
   real dI, uI[3], pI, cI, dB, uB[3], pB, cB;
   real VnI, dV2I, VnB, dV2B, Ma;

   //initialize UB
   pB = p;
   for(k=0; k<sr; k++) UB[k] = UI[k];

   //calculate interior normal velocity, VI dot N
   dI = UI[0];
   for(i=0, VnI=0.0; i<dim; i++) VnI += UI[1+i]/dI*n[i];
#ifdef test
   if (VnI < 0.0)
     xf_printf("Warning: reversed flow on Outflow boundary.\n");
#endif
   //compute interior pressure and speed of sound
   for (i=0, dV2I=0.; i<dim; i++) dV2I += UI[1+i]*UI[1+i]/dI;
   pI = (gamma - 1.) * (UI[1+dim] - 0.5*dV2I);

   cI = sqrt(gamma*pI/dI);
   
   //if Normal Mach number is supersonic, just use noneBC
   Ma = VnI/cI;
   if(Ma >= 1.0) return xf_OK;

   //set exterior density based on enstropic expansion
   dB = dI*pow(pB/pI, 1./gamma);
 
   //compute the exterior speed of sound
   cB = sqrt(gamma*pB/dB);
 
   //evaluate exterior normal velocity
   VnB = 2./(gamma-1.)*(cI - cB);
   for (i=0; i<dim; i++)
     UB[1+i] = dB*VnB*n[i] + dB*UI[1+i]/dI;

   //set exterior energy
   for (i=0, dV2B = 0.; i<dim; i++)
     dV2B += UB[1+i]*UB[1+i]/dB;  
   UB[1+dim] = pB/(gamma-1.) + 0.5*dV2B;

   //for passive scalar
   for (i=2+dim; i<sr; i++)
      UB[i] = dB/dI * UI[i];

   return xf_OK;
}
/******************************************************************/
//function: stagnition inflow boundary
static int
StagInflowBC(const int dim, const int sr, const real *UI,
                   const real Tt, const real pt, const real *angle,
                   const real *n, real *UB, const real gamma)
{
    int i, k, j;
    int ir, irV[3], irE, *PIS;
    real R, rB, cB, pB, TB, rI, pI, cI, VnI, rVI2;
    real a, b, c, disc;
    real MB1, MB2, MB;
    real Mfac;
    real NB[3], dn, Jp;
    real VB[3], VB2, fac;
    real gam, igam, gmi, igmi;
    
    //sr  = EqnSet->StateRank;
    //sr2 = sr*sr;
    //dim = EqnSet->Dim;
    //PIS = EqnSet->PosInState;
    
    // state indices
    ir  = 0; //PIS[xfe_Density];
    //xf_MomentumPIS(PIS, dim, irV);
    //irE = PIS[xfe_Energy];
    if(dim == 2){
        irV[0] = 1; irV[1] = 2; irE = 3;
    }
    if(dim == 3){
        irV[0] = 1; irV[1] = 2; irV[2] = 3; irE = 4;
    }
    
    // pull off parameters
    //R   = RParam[xfe_GasConstant];
    R = 1.0;
    gam = gamma;
    igam = 1./gam;
    gmi  = (gam-1.);
    igmi = 1./gmi;
    
    // initialize UB, and UB_UI to interior, identity, respectively
    for (k=0; k<sr; k++) UB[k] = UI[k];
    
    // calculate interior normal velocity, VI dot N (will be negative, as N points out)
    rI = UI[ir];
    for (i=0, VnI=0.; i<dim; i++) VnI += UI[irV[i]]/rI*n[i];
    
    //be quiet here;
    //if (VnI > 0.){
    //    xf_printf("Warning: reversed flow on Inflow boundary.\n");
    //}
    
    // compute interior pressure
    for (i=0, rVI2=0.; i<dim; i++) rVI2 += UI[irV[i]]*UI[irV[i]]/rI;
    pI = gmi*(UI[irE] - 0.5*rVI2);

    if (pI < 0.) return xf_Error(xf_NON_PHYSICAL);
    
    // compute interior speed of sound
    cI = sqrt(gam*pI/rI);
    
    // compute exterior direction, NB, using angles angle[0], angle[1]
    NB[0] = cos(angle[0]);
    NB[1] = sin(angle[0]);
    if (dim == 3){
        NB[0] *= cos(angle[1]);
        NB[1] *= cos(angle[1]);
        NB[2]  = sin(angle[1]);
    }
    
    for (i=0, dn=0.; i<dim; i++) dn += NB[i]*n[i];  // dn = NB dot N
    
    // compute interior Riemann invariant, Jp = VnI + 2*cI/(gam-1)
    Jp = VnI + 2.*cI*igmi;
    
    // solve for MB = exterior Mach number
    a = gam*R*Tt*dn*dn - 0.5*gmi*Jp*Jp;   // MB^2 coeff
    b = 4.*gam*R*Tt*dn*igmi;              // MB^1 coeff
    c = 4.*gam*R*Tt*igmi*igmi - Jp*Jp;    // MB^0 coeff
    disc = b*b-4.*a*c;                    // discriminant
    
    // error if no solution
    if (disc <= 0.0){
        //xf_printf("No (or degenerate) solution for MB in Inflow state calculation.\n");
        //return xf_Error(xf_NON_PHYSICAL);
    
        //propose a trick to fix this
        MB = 0.001;
    }
    else{
        // two possible solutions, but will choose smaller positive one
        MB1 = 0.5*(-b-sqrt(disc))/a;
        MB2 = 0.5*(-b+sqrt(disc))/a;
    
        if ((MB1 < 0.) && (MB2 < 0.)){
            //xf_printf("Error, negative mach number, MB, at inflow.\n");
            //return xf_Error(xf_NON_PHYSICAL);
        
            //the same trick to fix this
            MB = 0.001;
        }
        else if (MB1 < 0.){ // MB2 is non-negative
            MB = MB2; fac =  1.0;
        }
        else if (MB2 < 0.){ // MB1 is non-negative
            MB = MB1; fac = -1.0;
        }
        else{               // both non-negative
            MB  = ((MB1 < MB2) ? MB1 : MB2);
            fac = ((MB1 < MB2) ? -1. :  1.);
        }
    }
    
    // Mfac = 1 + (gamma-1)/2 * MB^2 is used below
    Mfac = 1. + .5*gmi*MB*MB;
    
    // compute exterior temperature using MB and Tt
    TB = Tt/Mfac;
    
    // compute exterior pressure using MB and pt
    pB = pt*pow(Mfac, -gam*igmi);
    
    // compute exterior density using pB, TB, R
    rB = pB/(R*TB);
    
    // compute exterior speed of sound using pB, rB
    cB = sqrt(gam*pB/rB);
    
    // compute exterior velocity using NB (direction) and cB, MB
    for (i=0; i<dim; i++) VB[i] = MB*cB*NB[i];
    for (i=0, VB2=0.; i<dim; i++) VB2 += VB[i]*VB[i];
    
    // set exterior state and derivatives
    UB[ir] = rB;
    for (i=0; i<dim; i++) UB[irV[i]] = rB*VB[i];
    UB[irE] = pB*igmi + 0.5*rB*VB2;
    
    return xf_OK;
}


/******************************************************************/
//function: specify visous flux at boundary
int
VisFluxBoundaryState(const int nq, const int sr, const int dim, const real *wn,
                     const real *xglob, const real *qgUI, real *qgUB, enum xfe_Bool *VisFluxFlag,
                     int *SetIndex, const int localBCtype, const real *localBCpara, enum xfe_Bool AVFlag)
{
   int i, j, k, ierr, iq;
   real *gUI, *gUB;

  // printf("%d\n", localBCtype);
  // getchar();

   for(iq=0; iq<nq; iq++){

      //set pointers
      gUI = qgUI + iq*sr;
      gUB = qgUB + iq*sr;
 
 if(!AVFlag)     
 {     if(localBCtype == fullBC || localBCtype == exactBC)
      {
         if(iq == 0) {
            (*VisFluxFlag) = xfe_True;
         for(k=0; k<sr; k++)
            SetIndex[k] = 0;
         }
      }
      else if(localBCtype == noneBC || localBCtype == isothermwallBC)
      {
         gUB[0] = 0.;

         if(iq == 0) {
            (*VisFluxFlag) = xfe_True;
         for(k=0; k<sr; k++)
            SetIndex[k] = 0;
         
         //previous wrong setting
         //SetIndex[0] = 1;
         }

      }
      else if(localBCtype == slipwallBC || localBCtype == noslipwallBC)
      {
         //specific boundary flux
         //only adiabatic is used
         gUB[1+dim] = 0.0;
         
         if(iq == 0) {
            (*VisFluxFlag) = xfe_True;
         for(k=0; k<sr; k++)
            SetIndex[k] = 0;
         
         SetIndex[1+dim] = 1;
         }
/*if(AVFlag) 
{
    gUB[0] = 0.0;
    for (j=0; j<dim; j++) gUB[1+j] = 0.0;
    gUB[1+dim] = 0.0;

    if(iq == 0) {
       (*VisFluxFlag) = xfe_True;
       for(k=0; k<sr; k++)
          SetIndex[k] = 0;

       SetIndex[0] = 1;
       for (j=0; j<dim; j++) SetIndex[1+j] = 1;
       SetIndex[1+dim] = 1;
    }
}*/
      }
      else if(localBCtype == staticPBC || localBCtype == farPBC || localBCtype == subinflowBC || localBCtype == staginflowBC)
      {
         //zero out momentum flux
         for (j=0; j<dim; j++) gUB[1+j] = 0.0;

         //zero out species density transport flux
         //do not know whether affect
         for (k=2+dim; k<sr; k++) gUB[k] = 0.0;

         //zero out energy flux
         gUB[1+dim] = 0.0;
         
         if(iq == 0) {
            (*VisFluxFlag) = xfe_True;
         for(k=0; k<sr; k++)
            SetIndex[k] = 0;
      
         for (j=0; j<dim; j++) SetIndex[1+j] = 1;
         //do not know whether affect
         for (k=2+dim; k<sr; k++) SetIndex[k] = 1;
         
         SetIndex[1+dim] = 1;
         }
      }
      else
          return xf_Error(xf_OUT_OF_BOUNDS);
 } 
 else
 { 
    //AVmodel is actived; not jointly use for setting physical viscous boundary
    //zero out all the viscous flux for numerical stabilty
    gUB[0] = 0.0;
    for (j=0; j<dim; j++) gUB[1+j] = 0.0;
    gUB[1+dim] = 0.0;

    if(iq == 0) {
       (*VisFluxFlag) = xfe_True;
       for(k=0; k<sr; k++)
          SetIndex[k] = 0;

       SetIndex[0] = 1;
       for (j=0; j<dim; j++) SetIndex[1+j] = 1;
       SetIndex[1+dim] = 1;
    }
 }
   }


   return xf_OK;
}

/******************************************************************/
// Entrance to Evaluate Boundary States for Various BC types
int
ConvFluxBoundaryState(Yu_Model *Model, const int nq, const real *qUI, real *qUB, 
                      const real *qn, const real *xglob, const real Time, const int BCindex,
                      const real Gamma)
{
   int iq, ierr, d, j, sr, dim;
   int localBCtype;
   const real *UI, *n, *xg;
   real *UB, tmp, b, *localBCpara, gimoR;
   real NN, N[3], uInput[3], YInput[50];
   real eps = 1.e-13, dt, theta;
   enum xfe_Bool UserSpfBCFun;

   //abstrat model data from model
   sr = Model->nVars;
   dim = Model->dim;
   UserSpfBCFun = Model->BCFunSpf;
   dt = Model->dt_size;

   //this boundary info
   localBCtype = Model->typeBCs[BCindex];
   localBCpara = Model->paraBCs[BCindex];

   for(iq=0; iq<nq; iq++){
       
       //set pointers
       UI    = qUI + iq*sr;
       UB    = qUB + iq*sr;
       xg    = xglob + iq*dim;
 
       // compute normalized normal vector
       n = qn+dim*iq;
       for (d=0, NN=0.0; d<dim; d++) NN += n[d]*n[d];
       NN = sqrt(NN);
       for (d=0; d<dim; d++) N[d] = n[d]/NN;
       
       //set state
       if(localBCtype == fullBC)
       {
          if(UserSpfBCFun)
          {
    /*         UB[0] = 1.0;
             UB[1] = 2.36643273 * tanh(xg[1]/0.01);
             UB[2] = 0.0;
             UB[3] = 2.5 + 0.5 * UB[1] * UB[1] / UB[0];
             UB[4] = 1.0;
    */
             UB[0] = 1.;
             UB[1] = 3./2.0 + 1./2.0 * tanh(2.0 * xg[1] / 1.0);
             UB[2] = 0.;
             UB[3] = 11.4285714286/(1.4 - 1.) + 0.5 * (UB[1]*UB[1] + UB[2]*UB[2]) / UB[0];

          }
          else
          for(j=0; j<sr; j++)
              UB[j] = localBCpara[j];
       }
       else if(localBCtype == slipwallBC)
          SlipWallBC(UI, N, UB, sr, dim);
       else if(localBCtype == noslipwallBC)
          NoSlipWallBC(UI, N, UB, sr, dim);
       else if (localBCtype == isothermwallBC)
       {
          gimoR = (Model->GammaInit - 1.)/(Model->Ru /Model->moleW[0]) ;
          NoSlipIsoThermalWallBC(UI, N, UB, sr, dim, gimoR, localBCpara);
       }
       else if(localBCtype == noneBC)
       {          
          for(j=0; j<sr; j++)
              UB[j] = UI[j];
      
          if(Model->InterpolateInitData== 20)
          {
             tmp = (Gamma - 1.) * (UB[3] - 0.5*(UB[2]*UB[2] + UB[1]*UB[1])/UB[0]); 
             if(n[1] > eps) //upper 
                UB[2] = 0.1 * UB[0];

             if(n[1] < eps && xg[0] <= 0.0) //lower left
                UB[2] = 0.1 * UB[0];
             if(n[1] < eps && xg[0] > 0.0) //lower right
                UB[2] = -0.6259 * UB[0];
             
             if(n[0] > eps) //right
                UB[1] = 0.1 * UB[0];
             if(n[0] < eps && xg[1] <= 0.0) // left lower
                UB[1] = 0.1 * UB[0];
             if(n[0] < eps && xg[1] > 0.0) // left upper
                UB[1] = -0.6259 * UB[0];
          
             UB[3] = tmp/(Gamma-1.) + 0.5*(UB[2]*UB[2] + UB[1]*UB[1])/UB[0];
          }
       
       }
       else if(localBCtype == subinflowBC)
       {
   
          if(UserSpfBCFun)
          {
    //for turbulent jet   
  
             //user specification
   /*          tmp = sqrt(xg[0]*xg[0] + xg[1]*xg[1] + xg[2]*xg[2]);
             
             //if the perturbation is not invoked.
             b = 12.5;
             
             //compute the perturbed parameter b 
             //the thichness of the laminar layer
             //the implementation refers to J.B.Fruend's JFM 2001

             
             theta = acos(xg[1]/tmp);
             ierr = xf_Error(compute_b(theta, Time, &b));
             if (ierr != xf_OK) return ierr;
            
             uInput[2] =  0.5*(1.0 - tanh(b*(tmp - 1./tmp))) +  0.0001;
             uInput[1] = 0.0; uInput[0] = 0.0;

             tmp = 0.88183421516; 
             for(j=0; j<sr-2-dim; j++)
                YInput[j] = 1.0;
// */
             //following is for 2d mixing layer
  /*           uInput[0] = 3./2.0 + 1./2.0 * tanh(2.0 * xg[1] / 1.0);
             uInput[1] = 0.0;
             tmp = 11.4285714286; 
*/             
/*             for(j=0; j<sr-2-dim; j++)
                YInput[j] = 1.0;

             tmp = sqrt(xg[1]*xg[1] + xg[2]*xg[2]);
             tmp = 1. - tmp;
             if(tmp > 0.3)
             uInput[0] = 0.9;
             else
             {
                b = tmp / 0.05;
                uInput[0] = 0.9 * b * (2.-2.*b*b + b*b*b);
             }
*/             uInput[0] = 0.1;
             uInput[1] = 0.0;
             uInput[2] = 0.0;
             tmp = 0.71428571428;
          }
          else
          {
             //load params from 
             for(j=0; j<dim; j++)
                uInput[j] = localBCpara[j];
             //totol pressure
             tmp = localBCpara[dim];
             //passive scalar
             for(j=0; j<sr-2-dim; j++)
                YInput[j] = localBCpara[dim+1+j];
          }

          //call BC routine
          ierr = xf_Error(SubsonicInflowBC(dim, sr, UI, uInput, tmp, YInput, N, UB, Gamma));
          if (ierr != xf_OK) return ierr;
       }
       else if(localBCtype == staginflowBC)
       {
          tmp = localBCpara[0];
          b   = localBCpara[1]; 
          for(j=0; j<dim; j++)
              uInput[j] = localBCpara[2+j];
          ierr = xf_Error(StagInflowBC(dim, sr, UI, tmp, b, uInput, N, UB, Gamma));
          if (ierr != xf_OK) return ierr;
       }
       else if(localBCtype == staticPBC)
       {
          ierr = xf_Error(StaticPressureBC(dim, sr, UI, localBCpara[0], N, UB, Gamma));
          if (ierr != xf_OK) return ierr;
       }
       else if(localBCtype == farPBC)
       {
          if(UserSpfBCFun && strcmp(Model->nameBCs[BCindex],"right") == 0)
          {
             //specified for individul case
             localBCpara[1] = 3./2.0 + 1./2.0 * tanh(2.0 * xg[1] / 40.0); 
          }
          ierr = xf_Error(FarFieldPressureBC(dim, sr, xg, UI, localBCpara, N, UB, Gamma));
          if (ierr != xf_OK) return ierr;
       }
       else if(localBCtype == exactBC)
       {

if(LinearCase)
{
   //for AdvDiffRectChannel
   for(j=0; j<sr; j++)
      UB[j] =  0.;
   //left inflow
   if(fabs(xg[0]-0.0) < eps)
   {
      UB[0] = sin(PI*xg[1]);
   }
   else if(fabs(xg[1]-1.0) < eps) //top wall
   {
      UB[0] = 0.0;
   }
   else if(fabs(xg[1]-0.0) < eps) //bottom wall
   {
      UB[0] = 0.0;
   }
   else
      return xf_Error(xf_NOT_SUPPORTED); 
}
else
{
          //Double Mach Reflection setup
          if(fabs(xg[1]-1.0) < eps)  //top BC
          {
             tmp = 1./6. + (1. + 20.*Time)/sqrt(3.);
             if(xg[0]  < tmp)
             {
                UB[0] = 8.0;
                UB[1] = 57.15768;
                UB[2] = -33.0;
                UB[3] = 563.500024;
                UB[4] = 8.0;
             }
             else
             {
                UB[0] = 1.4;
                UB[1] = 0.0;
                UB[2] = 0.0;
                UB[3] = 2.5;
                UB[4] = 1.4;
             }

          }
          else if(fabs(xg[1]-0.0) < eps) //bottom BC
          {
             if(xg[0] < 1./6.)
             {
                UB[0] = 8.0;
                UB[1] = 57.15768;
                UB[2] = -33.0;
                UB[3] = 563.500024;
                UB[4] = 8.0;
             }
             else
             {
                SlipWallBC(UI, N, UB, sr, dim);
             }
          }
          else
             return xf_Error(xf_NOT_SUPPORTED); 
       
       
}//if LinearCase
       }
       else 
          return xf_Error(xf_OUT_OF_BOUNDS);
   
   }

   return xf_OK;
}

/******************************************************************/
// Evaluate Riemann Flux at the Boundary
/*int 
ConvFluxBoundaryFace(const int nq, const int sr, const int dim, const real *qUI, const real *qUB, 
                     const real *qn, real *qF, const real Gamma, real *MaxCharSpeed)
{
   int iq, ierr;
   const real *UI, *UB, *n;
   real *F;
   real NN, N[3];

   for(iq=0; iq<nq; iq++){

      //set pointers
      UI    = qUI + iq*sr;
      UB    = qUB + iq*sr;
      n     = qn  + iq*dim;
      F     = qF  + iq*sr;

       if(dim == 2){
           ierr = xf_Error(ConvFluxFace2D(sr, UI, UB, n, F, Gamma, MaxCharSpeed));
           if (ierr != xf_OK) return ierr;
       }
       else if(dim == 3){
           ierr = xf_Error(ConvFluxFace3D(sr, UI, UB, n, F, Gamma, MaxCharSpeed));
           if (ierr != xf_OK) return ierr;
       }
       else
           return xf_Error(xf_OUT_OF_BOUNDS);
   }

   return xf_OK;
}*/

/*******************************************************************/
//data struture for data reading from pre-defined file
/*
typedef struct
{
   real *DataPointer;
   int  Var_len;
   int  Var_num;
   real Offset;
}
FileData;
static FileData Init_from_File={NULL, 0, 0, 0.0};

//read 1D data from user-defined file for afterward initialization
int
ReadDataFromUserDefineFile(Yu_Model *Model, const char *FileName)
{
   int i, j, ierr;
   int Index_i, Index_j;
   real *DataU, x, dummy;
   FILE *fp;

   fp = fopen(FileName, "r");

   //read index
   Index_i = Model->DataFileNumRow;
   Index_j = Model->DataFileNumCol;
   x = Model->DataFileLocus;

   //allocate
   ierr = xf_Error(xf_Alloc((void **) &DataU, (Index_i)*(Index_j), sizeof(real)));
   if (ierr != xf_OK) return ierr;

   //skip the first line
   char ignore[1024];
   fgets(ignore, 1023, fp);

   //load data for file
   for(i=0; i<Index_i; i++)
   {
      for(j=0; j<Index_j; j++)
     { 
         fscanf(fp, "%lf", &DataU[i*Index_j+j]);
     }
   }

   fclose(fp);

   //set up data structure
   Init_from_File.Offset  = x;
   Init_from_File.Var_len = Index_i;
   Init_from_File.Var_num = Index_j;
   Init_from_File.DataPointer = DataU;

   return xf_OK;
}
void
DeleteFileDataAllocation()
{
   xf_Release((void *) Init_from_File.DataPointer);
}
*/

/*******************************************************************/
//Make the life easier; Please specify analytic expression here
//for initialization and convergence study
static int
AnalyticalExpression(const int dim, const real *xglob, real *qu)
{
   real para, phi, alpha, gam;
   real rho, u, v, w, p;
   real r, x0, y0, uinf, vinf, x, y;

   gam = 1.4;
   para = 1./(1.4-1.0);
   phi = 1.0;
   alpha = 4.0;
   x0 = 5.;  y0 = 5.;
   x = xglob[0]; y = xglob[1];
   r = sqrt(pow(x-x0,2.)+pow(y-y0,2.));

   //in the same rank
   rho = pow((1.-alpha*alpha*(gam-1.)/16./phi/gam/PI/PI*exp(2.*phi*(1.-r*r))), para);
   u = 1.-alpha/2./PI*(y-y0)*exp(phi*(1.-r*r));
   v = 1.+alpha/2./PI*(x-x0)*exp(phi*(1.-r*r));
   if(dim == 3)
      w = 0.0;
   else
      w = 0.0;
   p = pow(rho, gam);

   //conservative variables
   qu[0] = rho;
   qu[1] = rho*u;
   qu[2] = rho*v;
   if(dim == 3)
      qu[3] = rho*w;
   qu[dim+1] = p/(gam-1.) + 0.5*rho*(u*u+v*v+w*w);
  
   return xf_OK;

}

/*******************************************************************/
// Specify initial conditions
// Usually need to modify and re-compile the code before running
int
InitializeDataInterpolation(Yu_Model *Model, const int nq, const int sr, const int dim, 
                            const real *qxglob, real *qU, const real Gamma)
{
   int ierr, iq, k, Nspe, i, j;
   real r, theta, tmp, meanW;
   real rho, p, u, v, Y, x0, y0, Velmov;
   const real *xglob;
   real *U, pref, rhoref, uref, unit;
   enum xfe_Bool found;

   Nspe   = sr - 4 + 1;
   pref   = 101325.0;
   rhoref = 1.29;
   uref   = sqrt(pref/rhoref);

   // loop over nq
   for (iq=0; iq<nq; iq++){

      U = qU + iq*sr;
      xglob = qxglob + iq*dim;


      //specify initial condition here
      //Shock bubble case (Marquina, Mulet, JCP, 2003)
      //1.22 Mach shock
/*      if(xglob[0]<0.225)
      U[0] = 1.0*rhoref;
      else 
         U[0] = 1.3764*rhoref;

     if(xglob[0]<0.225)
     U[1] = 0.0;
     else
        U[1] = U[0]*(-0.3336)*uref;

      U[2] = 0.0;

     if(xglob[0]<0.225)
        U[3] = (pref/1.4)/(Gamma - 1.0);
     else
        U[3] = (1.5698*pref/1.4)/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0]; 
 
        U[4] = U[0]*1.0;
  
     //make the cylinder
     r = sqrt(pow((xglob[0] - 0.195), 2.0) + pow((xglob[1]), 2.0));
     if(r<0.025)
     {
        U[0] = 0.1819*rhoref;
        U[1] = 0.0;
        U[2] = 0.0;
        U[3] = (pref/1.4)/(Gamma - 1.0);
        U[4] = U[0]*0.0;
     }
*/
      //for shock diffraction case with 2.4 Mach shock
      //Aug 2013, proceeding work mesh test
/*
      if(xglob[0] < 0.01){
         U[0] = 3.212;
         U[1] = 6.281461;
         U[2] = 0.0;
         U[3] = 22.52601;
         U[4] = 3.212;
      }
      else
      {
         U[0] = 1.0;
         U[1] = 0.0;
         U[2] = 0.0;
         U[3] = 2.5;
         U[4] = 0.0;
      }
*/
      //Aug 2013, proceeding work one-step chemistry
      //might need finer mesh      
/*      if(xglob[0] < 0.005){
         U[0] = 0.230891;
         U[1] = U[0] * 1.263967820547646e+03;
         U[2] = 0.0;
         U[3] = 462886.0/(Gamma - 1.) + 0.5 * U[1] * U[1] / U[0];
         U[4] = 0.0;
         U[5] = 0.0;
         U[6] = U[0];
      }
      else
      {
         U[0] = 0.125189;
         U[1] = 0.0;
         U[2] = 0.0;
         U[3] = 26000.0/(Gamma - 1.0);
         U[4] = 0.111898 * U[0];
         U[5] = 0.888102 * U[0];
         U[6] = 0.0;
      }
*/
      //Aug 2013, proceeding work detailed chemistry
/*       unit = 1.;
       if(Init_from_File.DataPointer != NULL)  //mean data is loaded
       {
           found = xfe_False;
           
           for(i=0; i<Init_from_File.Var_len - 1; i++)
               if(Init_from_File.DataPointer[i*Init_from_File.Var_num] >= (xglob[0]+Init_from_File.Offset)*unit  &&
                  Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num] <= (xglob[0]+Init_from_File.Offset)*unit)
               {
                   found = xfe_True;
                   break;
               }
           
           if(found && Model->DetailChem)
           {
               tmp = ((xglob[0]+Init_from_File.Offset)*unit - Init_from_File.DataPointer[i*Init_from_File.Var_num])
               / (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num] - Init_from_File.DataPointer[i*Init_from_File.Var_num]);
               
               rho = Init_from_File.DataPointer[i*Init_from_File.Var_num + 4]; //look at the header of input file
              
               u   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 1];
               //u   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 4];
               p   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 5];
               
               U[0] = tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 4] - rho) + rho;
               U[1] = U[0] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 1] - u) + u);
               //if(xglob[0] < (-Init_from_File.Offset+0.003) && xglob[0] > (-Init_from_File.Offset-0.001))
               //    U[1] = 2.0 * sin(3.14*xglob[1]/0.001);
               //else
               //    U[1] = 0.;
               if(xglob[0] < 0.006 && xglob[0] > 0.004)
                  U[1] += 0.1 * U[1] * sin(3.14*xglob[1] / 0.0005);
               U[2] = 0.0;
               U[3] = (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 5] - p) + p)/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0];
               
               meanW = 0.;
               for(j=0; j<Nspe; j++)
               {
                   Y      = Init_from_File.DataPointer[i*Init_from_File.Var_num + 6 + j];
                   meanW += Model->moleW[j] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 6 + j] - Y) + Y);
               }
               
               for(j=0; j<Nspe-1; j++)
               {
                   Y      = Init_from_File.DataPointer[i*Init_from_File.Var_num + 6 + j];
                   U[4+j] = Model->moleW[j] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 6 + j] - Y) + Y);
                   U[4+j] = U[0] * U[4+j] / meanW;
               }
           }
           else if(found && !Model->DetailChem)
           {
               tmp = ((xglob[0]+Init_from_File.Offset)*unit - Init_from_File.DataPointer[i*Init_from_File.Var_num])
               / (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num] - Init_from_File.DataPointer[i*Init_from_File.Var_num]);
               rho = Init_from_File.DataPointer[i*Init_from_File.Var_num + 3]; //look at the header of input file
               u   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 1];
               p   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 4];
               Y   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 5];
               
               U[0] = tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 3] - rho) + rho;
               U[1] = U[0] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 1] - u) + u);
               
               if(xglob[0] < 0.006 && xglob[0] > 0.004)
               U[1] += 0.1 * U[1] * sin(3.14*xglob[1] / 0.0005);
               
               U[2] = 0.0;
               U[3] = (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 4] - p) + p)/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0];
               U[4] = U[0] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 5] - Y) + Y);
               
           }
           else
               return xf_Error(xf_NOT_SUPPORTED);
       }
*/
       //Interpolate IC with analytical expression
      //ierr = xf_Error(AnalyticalExpression(dim, xglob, U));
      //if(ierr != xf_OK) return ierr;
      //U[0] = 1.4; U[1] = 4.2; U[2] = 0.0; U[3]= 8.8; U[4] = 1.4;
      //if(xglob[0] < 0.6 && xglob[0]> 0.59 && xglob[1] < 0.2)
      //{ U[1] = 0.0; U[3] = 2.5;}

      //if initialization is conducted using initial from file
      if(Model->InterpolateInitData == 2)
         Yu_FlowFieldInit_from_datafile(Model->init_datafile, sr, Model->dim, xglob, U);
      else 
      //Initialization for Double Mach reflection
         Yu_FlowFieldInit(Model->InterpolateInitData, sr, Model->dim, xglob, U);

      //order of species: H        H2       O        OH
      //H2O      O2       HO2      H2O2     N2     
      //for shock flame interaction
      //if(xglob[0] > 0.0077)
/*      if(xglob[0] > 0.007)
      {
         U[0] = 0.849057028866;
         U[1] = 0.0;
         U[2] = 0.0;
         U[3] = 253312.5;
         U[4] = 0.0;  U[5] = 0.028301887 * U[0]; U[6] = 0.0; U[7] = 0.0;
         U[8] = 0.0;  U[9]= 0.2264151 * U[0];    U[10] = 0.0; U[11] = 0.0;
      }
      else
      {
         U[0] = 3.1724782019;
         U[1] = 838.46 * U[0];
         U[2] = 0.0;
         U[3] = 9.1635250106e+5/0.4 + 0.5 * U[1] * U[1] / U[0];
         U[4] = 0.0;  U[5] = 0.028301887 * U[0]; U[6] = 0.0; U[7] = 0.0;
         U[8] = 0.0;  U[9]= 0.2264151 * U[0];    U[10] = 0.0; U[11] = 0.0;
      }

      //CEA output: T 2379.34 K; rho 1.2431e-1; mass fraction: 
      //H 7.3268e-5; H2 1.2584e-3; H2O 2.4020e-1; O 3.4621e-4;
      //OH 5.2329e-3; O2 6.1279e-3
     r = sqrt(pow((xglob[0] - 0.0111), 2.0) + pow((xglob[1]), 2.0));
     if(r<0.003)
     {
        U[0] = 101325.0/(8.31446/24.271e-3)/2380.0;
        U[1] = 0.0;
        U[2] = 0.0;
        U[3] = 101325.0/0.4;
        U[4] = U[0]*7.3268e-5;  U[5] = U[0]*1.2584e-3; U[6] = U[0]*3.4621e-4; U[7] = U[0]*5.2329e-3;
        U[8] = U[0]*2.4020e-1;  U[9]=  U[0]*6.1279e-3; U[10] = U[0]*1.4955e-6; U[11] = U[0] * 1.6851e-7;
     }
 */     
      //
       
/*     
   //simple test (normal shock)
   U[0] = 1.0;

//   r = xglob[0]+xglob[1];
   if(xglob[0]+0.25*xglob[1]<0.1)
//   if(xglob[0]<0.1)
      U[1] = 40.0;
   else
      U[1] = 0.0;

  //    U[1] = 20.0 - 20.0*tanh((r-0.2)/0.005);

   U[2] = 0.0;

   if(xglob[0]+0.25*xglob[1]<0.1)
//   if(xglob[0]<0.1)
      U[3] = 1000.0;
   else
      U[3] = 200.0;

//   U[3] = 600.0 - 400.0*tanh((r-0.2)/0.005);
*/
/*
      U[0] = 10.0*exp(-(pow((xglob[0]-0.005)/0.001, 2.0) + pow((xglob[1]-0.0017)/0.001, 2.0))) + 1.0;
      U[1] = 10.0;
      U[2] = 10.0;
      U[3] = 1000.0;
      U[4] = 1.0;
*/
      //specify for detonation with simple chemistry
      //Velmov = -400.0;
/*      
      if(xglob[0]<0.003)
         U[0] = 0.84;// + 0.2*sin(xglob[1]/0.003*6*3.14);
      else
         U[0] = 0.493;

      //U[1] = (-2845.0 + Velmov ) * U[0];// + 100.0 * sin(xglob[1]/0.003*6*3.14))*U[0];
      
         U[1] = 0.0;
      //add small disturbance
      if(xglob[0]>0.0004 && xglob[0]<0.0007)
         U[2] = 1.0*sin(xglob[1]/0.003*6*3.14);
      else 
         U[2] = 0.0;

      if(xglob[0]<0.0005)
         U[3] = 34.0e+5/(Gamma - 1.0);
      else
         U[3] = 1.0e+5/(Gamma - 1.0);

      U[3] = U[3] + 0.5*(U[1]*U[1] + U[2]*U[2])/U[0];

      if(xglob[0]<0.0005)
         U[4] = U[0];
      else
         U[4] = 0.0;
*/
      //test the far field pressure boundary condition
/*      U[0] = 1.0;

      U[1] = 0.0;

      if(xglob[0]>0.7 && xglob[0]<1.2)
         U[2] = 1.0*sin(xglob[1]/0.2*6*3.14);
      else 
         U[2] = 0.0;

      U[3] = 101325.0/(0.4) + 0.5*(U[1]*U[1]+U[2]*U[2])/U[0];

      U[4] = U[0];
*/
     /* sb_rho  += U[0];
      sb_rhou += U[1];
      sb_rhov += U[2];
      sb_rhoE += U[3];
      sb_air  += U[4];
      *///test with chemkin on March 5
/*         U[0] = 1.22;
         U[1] = 0.0;
         U[2] = 0.0;

        if(xglob[0]<0.003)
           U[3] = 3.0e+5 / (Gamma - 1.0);
        else
            U[3] = 1.0e+5/(Gamma - 1.0);

        if(xglob[0]<0.003)
           U[4] = 0.0;
        else
           U[4] = 1.22 * 0.033;

         U[5] = 1.22 * 0.21;

         U[6] = 0.0;
*/
      //initialization using user-defined data file
      //for perturbed premixed flame H2/O2 at phi = 0.4 with 70% Ar dillution
      /*
      unit = 100.0;
      if(Init_from_File.DataPointer != NULL)  //mean data is loaded
      {
         found = xfe_False;

         for(i=0; i<Init_from_File.Var_len - 1; i++)
            if(Init_from_File.DataPointer[i*Init_from_File.Var_num] <= (xglob[0]+Init_from_File.Offset)*unit  &&
               Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num] >= (xglob[0]+Init_from_File.Offset)*unit)
            {
               found = xfe_True;
               break;
            }

         if(found)
         {
            tmp = ((xglob[0]+Init_from_File.Offset)*unit - Init_from_File.DataPointer[i*Init_from_File.Var_num])
                / (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num] - Init_from_File.DataPointer[i*Init_from_File.Var_num]);

            rho = Init_from_File.DataPointer[i*Init_from_File.Var_num + 2]; //look at the header of input file
            u   = 0.0;
            //u   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 4];
            p   = 101325.0;
            // p   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 3];


            U[0] = tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 2] - rho) + rho;
            U[0] *= 1000.0;
            //U[1] = U[0] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 2] - u  ) + u);
            if(xglob[0] < (-Init_from_File.Offset+0.003) && xglob[0] > (-Init_from_File.Offset-0.001))
            U[1] = 2.0 * sin(3.14*xglob[1]/0.001);
            else
               U[1] = 0.;
            U[2] = 0.0;
            U[3] = p/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0];;
            //U[3] = (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 3] - p) + p)/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0];

            meanW = 0.;
            for(j=0; j<Nspe; j++)
            {
               Y      = Init_from_File.DataPointer[i*Init_from_File.Var_num + 5 + j];
               meanW += Model->moleW[j] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 5 + j] - Y) + Y);
            }           

            for(j=0; j<Nspe-1; j++)
            {
               Y      = Init_from_File.DataPointer[i*Init_from_File.Var_num + 5 + j];
               U[4+j] = Model->moleW[j] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 5 + j] - Y) + Y);
               U[4+j] = U[0] * U[4+j] / meanW;
            }
         }
         else
            return xf_Error(xf_NOT_SUPPORTED);
      }
      */
      //test 1D premixed flame of H2/O2/Ar with different phi and dillution
/*      unit = 1.;
      if(Init_from_File.DataPointer != NULL)  //mean data is loaded
      {
         found = xfe_False;

         for(i=0; i<Init_from_File.Var_len - 1; i++)
            if(Init_from_File.DataPointer[i*Init_from_File.Var_num] <= (xglob[0]+Init_from_File.Offset)*unit  &&
               Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num] >= (xglob[0]+Init_from_File.Offset)*unit)
            {
               found = xfe_True;
               break;
            }

         if(found)
         {
            tmp = ((xglob[0]+Init_from_File.Offset)*unit - Init_from_File.DataPointer[i*Init_from_File.Var_num])
                / (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num] - Init_from_File.DataPointer[i*Init_from_File.Var_num]);

            rho = Init_from_File.DataPointer[i*Init_from_File.Var_num + 4]; //look at the header of input file
            u   = 0.0;
            //u   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 4];
            p   = 101325.0;
            // p   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 3];


            U[0] = tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 4] - rho) + rho;
            U[1] = 0.;
            U[2] = 0.0;
            U[3] = p/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0];;
            //U[3] = (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 3] - p) + p)/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0];

            meanW = 0.;
            for(j=0; j<Nspe; j++)
            {
               Y      = Init_from_File.DataPointer[i*Init_from_File.Var_num + 5 + j];
               meanW += Model->moleW[j] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 5 + j] - Y) + Y);
            }           

            for(j=0; j<Nspe-1; j++)
            {
               Y      = Init_from_File.DataPointer[i*Init_from_File.Var_num + 5 + j];
               U[4+j] = Model->moleW[j] * (tmp * (Init_from_File.DataPointer[(i+1)*Init_from_File.Var_num + 5 + j] - Y) + Y);
               U[4+j] = U[0] * U[4+j] / meanW;
            }
         }
         else
            return xf_Error(xf_NOT_SUPPORTED);
      }
*/
            /*
         else
         {
            rho = Init_from_File.DataPointer[i*Init_from_File.Var_num + 1];
            u   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 2];
            p   = Init_from_File.DataPointer[i*Init_from_File.Var_num + 3];

            U[0] = rho;
            U[1] = rho*u;
            U[2] = 0.0;
            U[3] = p/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0];

            for(j=0; j<Nspe-1; j++)
            {
               Y      = Init_from_File.DataPointer[i*Init_from_File.Var_num + 7 + j];
               U[4+j] = Y * U[0];
            }
         }
         */
         //add perturbation according to Houim
         //fresh unburnt gas exsiting at a box right behind the shock front
         //P = 47kPa; T = 2084K; rho = 1.0076 * rho_unburnt
/*         if(xglob[0]>=(0.016 - 0.003) && xglob[0] <= 0.016 &&
            xglob[1]>=(0.005 - 0.0015) && xglob[1] <= (0.005 + 0.0015))
         {
            rho = 1.0076 * Init_from_File.DataPointer[(Init_from_File.Var_len - 1)*Init_from_File.Var_num + 1];
            u   = Init_from_File.DataPointer[(Init_from_File.Var_len - 1)*Init_from_File.Var_num + 2];
            p   = 47000.0;

            U[0] = rho;
            U[1] = rho*u;
            U[2] = 0.0;
            U[3] = p/(Gamma - 1.0) + 0.5*U[1]*U[1]/U[0];

            for(j=0; j<Nspe-1; j++)
            {
               Y      = Init_from_File.DataPointer[(Init_from_File.Var_len - 1)*Init_from_File.Var_num + 7 + j];
               U[4+j] = Y * U[0];
            }
         }

      }
*/
   }

   return xf_OK;
}

//*******************************************************
//Error evaluation: L1 or L2 norm with analytical expression
//for the whole computational domain
//subject to change
int
EvaluateError(xf_All *All, xf_Vector *U, enum xfe_Bool PrintFlag)
{
    int ierr, dim, sr, i, k, d, iq, nq, pnq;
    int egrp, elem, nu, ns;
    int Order, QuadOrder, myRank, nProc;
    enum xfe_BasisType Basis;
    enum xfe_Bool QuadChanged;
    real *EU, *xq, *u, *wq, enorm2;
    real *uA = NULL, *xglob;
    real xcenter[3] = {0., 0., 0.};
    real min1=1e30, max1=-1e30, min2=1e30, max2=-1e30;
    real xmin1[3], xmax1[3], xmin2[3], xmax2[3];
    real e, enormA[5], integral, tmp1, tmp2;
    FILE *fout;
    xf_QuadData *QuadData;
    xf_BasisData *PhiData;
    xf_BasisData *GeomPhiData;
    xf_JacobianData *JData;
    xf_Mesh *Mesh;
    
    Mesh = All->Mesh;
    dim  = Mesh->Dim;
    
    // both vectors should have the same state rank
    sr   = U->StateRank;
    
    ierr = xf_Error(xf_Alloc( (void **) &uA, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    QuadData = NULL;
    PhiData = NULL;
    JData    = NULL;
    u       = NULL;
    wq       = NULL;
    xglob    = NULL;
    GeomPhiData = NULL;
    enorm2   = 0;
    for(i=0; i<5; i++)
       enormA[i] = 0;
    e       = 0;
    integral = 0.;
    pnq = -1;
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
        
        Basis = U->Basis[egrp]; 
        //Order = U->Order[egrp];
        
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
            /* Pull off quad points for the element; will not recalculate in generic case */
            //QuadOrder = Order;

            //modified to handle variable order elements
            Order = xf_InterpOrder(U, egrp, elem);
            QuadOrder = Order;

            ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
            if (ierr != xf_OK) return ierr;
            
            nq = QuadData->nquad;
            xq = QuadData->xquad;
            
            // compute basis functions for U vector
            ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged,
                                         nq, xq, xfb_Phi, &PhiData));
            if (ierr != xf_OK) return ierr;
            
            // element Jacobian
            ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ, QuadChanged, &JData));
            if (ierr != xf_OK) return ierr;
            
            // re-allocate memory if quad points increased
            if (nq > pnq){
                ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
                if (ierr != xf_OK) return ierr;
                ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
                if (ierr != xf_OK) return ierr;
                ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
                if (ierr != xf_OK) return ierr;
            }
            
            // form detJ-multiplied quad weight vector, wq
            // embed the element info (element size, etc.)
            for (iq=0; iq<nq; iq++)
                wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
            
            EU = U->GenArray[egrp].rValue[elem];
            
            xf_MxM_Set(PhiData->Phi, EU, nq, PhiData->nn, sr, u); // interpolate u
            
            // obtain global coords of quad points
            ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged,
                                            nq, xq, xglob));
            if (ierr != xf_OK) return ierr;
            
            // Analytical error if requested
            //if (Analytical)
            //{
            //for different variables
            for (iq=0; iq<nq; iq++){
                    
               //!!!!error evaluation by refering to analytical expression
               //for (k=0; k<dim+2; k++){
                    
                    //ierr = xf_Error(AnalyticalExpression(dim, xglob+iq*dim, uA));
                    ierr = xf_Error(Yu_FlowFieldInit(12, sr, dim, xglob+iq*dim, uA));
                    if (ierr != xf_OK) return ierr;
                    
                    //for (k=0; k<sr; k++)
                    k=0;
                    enormA[k] += (u[iq*sr+k]-uA[k])*(u[iq*sr+k]-uA[k])*wq[iq];
                    
                //}
            
               //!!!! error evaluation by entropy generation
               /*
               tmp1 = u[iq*sr+0]; tmp2 = 0.0;
               for (k=0; k<dim; k++)
                  tmp2 += pow(u[iq*sr+k+1],2.0);
               tmp2 = (1.4 - 1.0) * (u[iq*sr+dim+1] - 0.5*tmp2/tmp1);
               tmp1 = tmp2 / pow(tmp1, 1.4) / 0.62433941058;
               enormA[0] += (tmp1 - 1.0)*(tmp1 - 1.0) * wq[iq];
               */
            }
            //}
            
            // add to integrals of data1 and data2
            //for (iq=0; iq<nq; iq++)
            //    for (k=0; k<sr; k++)
            //        integral1 += (u1[iq*sr+k])*wq[iq];
            
            // add to integrals for xcenter1 and xcenter2
            //for (d=0; d<dim; d++)
            //    for (iq=0; iq<nq; iq++)
            //        for (k=0; k<sr; k++)
            //            xcenter1[d] += (u1[iq*sr+k])*wq[iq]*xglob[dim*iq+d];
            
            // take min and max
            //for (iq=0; iq<nq; iq++)
            //    for (k=0; k<sr; k++){
            //        if (u1[iq*sr+k] < min1) for(d=0;d<dim;d++) xmin1[d]=xglob[dim*iq+d];
            //        if (u1[iq*sr+k] > max1) for(d=0;d<dim;d++) xmax1[d]=xglob[dim*iq+d];
            //        min1 = min(u1[iq*sr+k], min1);
            //       max1 = max(u1[iq*sr+k], max1);
            //    }
            
            // add to square norms of data1 and data2
            //for (iq=0; iq<nq; iq++)
            //    for (k=0; k<sr; k++)
            //        e1 += (u1[iq*sr+k])*(u1[iq*sr+k])*wq[iq];
            
            pnq = nq;
        } // elem
    } // egrp
  
    //Determine myRank
    ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
    if (ierr != xf_OK) return ierr;

    //because of parallel
    ierr = xf_Error(xf_MPI_Allreduce(&enormA[0], 1, xfe_SizeReal, xfe_MPI_SUM));
    //ierr = xf_Error(xf_MPI_Reduce(&enormA[0], &tmp1,  1, xfe_SizeReal, xfe_MPI_SUM, 0));
    if (ierr != xf_OK) return ierr;

    if(PrintFlag && myRank == 0)
    {
       //user specification
        xf_printf("L2 error norm = %.15E\n", sqrt(enormA[0]));
        fout = fopen("Err","w");
        for(i=0; i<dim+2; i++)
           fprintf(fout, "%.15E\n", sqrt(enormA[i]));
        fclose(fout);

    }
    
    /* Destroy geometry Jacobian Data */
    ierr = xf_Error(xf_DestroyJacobianData(JData));
    if (ierr != xf_OK) return ierr;
    
    // Only destroy QuadData if points are generic
    ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
    if (ierr != xf_OK) return ierr;
    
    // Destroy Basis Data
    ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
    if (ierr != xf_OK) return ierr;
    
    // release memory
    xf_Release( (void *) u);
    xf_Release( (void *) wq);
    xf_Release( (void *) xglob);
    xf_Release( (void *) uA);
    
   return xf_OK;
}
