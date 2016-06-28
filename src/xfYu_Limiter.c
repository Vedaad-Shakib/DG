/*
 * Limiter implementation by Yu
 * The detail could be found in Cockburn and Shu's earlier work
 * Date: Oct, 2012
 * Email: lvyu@umich.edu
 */
#include "stdlib.h"
#include "stdio.h"
#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_MPI.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Memory.h"
#include "xf_Math.h"
#include "xf_DataMath.h"
#include "xf_Quad.h"
#include "xf_Basis.h"
#include "xf_QuadRule.h"
#include "xf_MeshTools.h"
#include "xf_All.h"
#include "xf_Solver.h"
#include "xfYu_Model.h"
#include "xf_MathLapack.h"

#define thredhold 1.0e-11

static int
FillinCharTransMatrix(real *iR_x, real *R_x, real *iR_y, real *R_y, Yu_Model *Model, const real gamma, const real *U);
static real
minmod(real a, real b, real c);
static void
LeftMulti(real *R, real *U, int sr);
static int
TransposeMatrix(real *M, int sr);
static int
Evaluate_Derivative_Coeff(real *c, const int egrp, const int elem, xf_Vector *U, const int order, const int sr);
static int
Evaluate_Derivative_Coeff_VersionII(real *c, real *U, const int order, const int sr);
static int
BackOut_BasisFcn_Coeff(const real *c, const int egrp, const int elem, xf_Vector *U, const int order, const int sr, const int rank);
static int
Convert_Derivative_to_PhysicsSpace(real *moment, const real *Trans_sign, const int order, const int sr);
static int
Convert_Derivative_to_ReferenceSpace(real *moment, const real *Trans_sign, const int order, const int sr, const int rank);
static int
Finish_Limiting(const int *Limited, const int sr);
static int
NegativeDensityChecking(Yu_Model *Model, const real gamma, const int egrp, const int elem, xf_Vector *U,
                        const real avgrho, const real *avgrhoY, real *facrho, real *facrhoY);
static int
NegativePressureChecking(Yu_Model *Model, const real gamma, const int egrp, const int elem, xf_Vector *U, 
                         const real avgrho, const real avgp, real *facp);
static int 
WhetherInTrouble(int TroubleMarker[], const int sr);
static int 
PutToTroubleStore(real *Poly, real *TroubleStore, const int order, const int sr, const int k);


//struct to rearrange the element sequance for parallel speed-up
typedef struct{
    int NumElem;   //number of element in block of the current cpu
    int NumBrk;    //number of element without neighbor in halo
    //re-arrange the element: interior element (without neighbor in halo) put first (=0)
    //elements with neighbor in halo put last (=-1)
    int *RearrangedEgrg;
    int *RearrangedElem;
    
    //allocation for smoothness indicator
    real *SmoothIndtor_wq, *SmoothIndtor_gu, *SmoothIndtor_gu2;
    
    //allocation for postivity check
    real *PostCheckState;
}
RearrangedElemList;
//global case
static RearrangedElemList Relist;

int
xf_FullinLimiterStruct(xf_All *All, Yu_Model *Model, Yu_Limiter ***pLimiter)
{

   int ierr, i, j, k, d;
   int face, sr, nn, nq, nq_ref;
   int egrp, elem, dim, localBCflag, BCtype;
   int grpIndx, elmIndx, bgnIndx, endIndx;
   int ibfgrp, iiface, ibface;
   real eff=1.0e-6;
   real *xq, nvec[3], xq_ref[8], *xglob, a[4], b[4];
   int QdLin_nq, QdSrf_nq;
   real *QdLin_xq=NULL, *QdLin_wq=NULL;
   real *QdSrf_xq=NULL, *QdSrf_wq=NULL;
   enum xfe_Bool QuadChanged, IamL, flag;
   xf_JacobianData *JData;
   xf_BasisData *GeomPhiData;
   xf_IFace IFace;
   xf_BFace BFace;
   xf_QuadData *QuadData;
   xf_NormalData *NData;
   xf_Mesh *Mesh;
   Yu_Limiter **Limiter;

   Mesh = All->Mesh;
   sr   = Model->nVars;
   dim  = All->Mesh->Dim;
    nn   = (Model->order+1)*(Model->order+1);

   //check dim first before trying to do anything regarding limiter
   if(dim != 2)
   {
      printf("3D limiter is not provided!"); 
      return xf_Error(xf_NOT_SUPPORTED);
   }

   if(Model->order > 2)
   {
      printf("Higher order limiter is not supported!");
      return xf_Error(xf_NOT_SUPPORTED);
   }

   if(strcmp(Model->basis,"QuadLagrange") !=0)
   {
      printf("Limiter is developed only for Quads mesh!");
      return xf_Error(xf_NOT_SUPPORTED);
   }

   QuadData = NULL;
   NData    = NULL;
   GeomPhiData = NULL;
   xglob       = NULL;
   JData       = NULL;

   //remember here only for linear quads_limiting
   nq_ref = 4;
   xq_ref[0] = 0.0; xq_ref[1] = 0.0;
   xq_ref[2] = 1.0; xq_ref[3] = 0.0;
   xq_ref[4] = 1.0; xq_ref[5] = 1.0;
   xq_ref[6] = 0.0; xq_ref[7] = 1.0;
   
   //obtain Gauss points for negative pressure checking
   ierr = xf_Error(xf_QuadLine(2*Model->order + 1, &QdLin_nq, &QdLin_xq, &QdLin_wq));
   if(ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_QuadQuadrilateral(2*Model->order + 1, &QdSrf_nq, &QdSrf_xq, &QdSrf_wq));
   if(ierr != xf_OK) return ierr;

   Model->Num_negPckpnt = QdLin_nq*4 + QdSrf_nq;

   ierr = xf_Error(xf_Alloc( (void **) &Model->Coord_negPckpnt, Model->Num_negPckpnt * 2, sizeof(real)));
   if(ierr != xf_OK) return ierr;
    
    for(k=0; k<2*QdSrf_nq; k++)
        Model->Coord_negPckpnt[k] = QdSrf_xq[k]; 
    for(k=0; k<QdLin_nq; k++){
        //edge y=0
        Model->Coord_negPckpnt[2*k+2*QdSrf_nq]              = QdLin_xq[k];
        Model->Coord_negPckpnt[2*k+1+2*QdSrf_nq]            = 0.0;
        //edge y=1
        Model->Coord_negPckpnt[2*k+2*QdSrf_nq+2*QdLin_nq]   = QdLin_xq[k];
        Model->Coord_negPckpnt[2*k+1+2*QdSrf_nq+2*QdLin_nq] = 1.0;
        //edge x=0
        Model->Coord_negPckpnt[2*k+2*QdSrf_nq+4*QdLin_nq]   = 0.0;
        Model->Coord_negPckpnt[2*k+1+2*QdSrf_nq+4*QdLin_nq] = QdLin_xq[k];
        //edge x=1
        Model->Coord_negPckpnt[2*k+2*QdSrf_nq+6*QdLin_nq]   = 1.0;
        Model->Coord_negPckpnt[2*k+1+2*QdSrf_nq+6*QdLin_nq] = QdLin_xq[k];
    }
    
    ierr = xf_Error(xf_EvalBasis(xfe_QuadLagrange, Model->order, QuadChanged, Model->Num_negPckpnt, 
                                 Model->Coord_negPckpnt, xfb_Phi, &Model->Phi_negPckpnt));
    if (ierr != xf_OK) return ierr;
/*    
   for(k=0; k<QdLin_nq; k++){
       //edge y=0
       Model->Coord_negPckpnt[2*k]              = QdLin_xq[k];
       Model->Coord_negPckpnt[2*k+1]            = 0.0;
       //edge y=1
       Model->Coord_negPckpnt[2*k+2*QdLin_nq]   = QdLin_xq[k];
       Model->Coord_negPckpnt[2*k+1+2*QdLin_nq] = 1.0;
       //edge x=0
       Model->Coord_negPckpnt[2*k+4*QdLin_nq]   = 0.0;
       Model->Coord_negPckpnt[2*k+1+4*QdLin_nq] = QdLin_xq[k];
       //edge x=1
       Model->Coord_negPckpnt[2*k+6*QdLin_nq]   = 1.0;
       Model->Coord_negPckpnt[2*k+1+6*QdLin_nq] = QdLin_xq[k];
   }

   

   for(k=0; k<12; k++)
      printf("%lf %lf\n", Model->Coord_negPckpnt[2*k], Model->Coord_negPckpnt[2*k+1]);
   return xf_Error(xf_NOT_SUPPORTED);

   xf_BasisData *PhiData = NULL;
   ierr = xf_Error(xf_EvalBasis(xfe_QuadLagrange, 2, QuadChanged, nq_ref, xq_ref,
                                xfb_Phi, &PhiData));
   if (ierr != xf_OK) return ierr;
   for(k=0; k<36; k++)
      printf("%lf ", PhiData->Phi[k]);
   return xf_Error(xf_NOT_SUPPORTED);
*/
   //allocate for pointer to Limiter
   if (( (*pLimiter) = (Yu_Limiter **)malloc(Mesh->nElemGroup * sizeof(Yu_Limiter *))) == NULL)
   return xf_Error(xf_MEMORY_ERROR);
   Limiter = (*pLimiter);

   //count how many elements are on this processor
    Relist.NumElem = 0;
    
   //loop over element groups
   for (egrp=0; egrp<Mesh->nElemGroup; egrp++){

      //assign pointer to Limiter
      if ((Limiter[egrp] = (Yu_Limiter *)malloc(Mesh->ElemGroup[egrp].nElem * sizeof(Yu_Limiter))) == NULL)
      return xf_Error(xf_MEMORY_ERROR);

      // loop over elements
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){

          //count element number
          Relist.NumElem++;
          
         //deal with the orientation stuff
         ierr = xf_Error(xf_Alloc( (void **) &xglob, 4*dim, sizeof(real)));
         if (ierr != xf_OK) return ierr;

         ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged,
                                         nq_ref, xq_ref, xglob));
         if (ierr != xf_OK) return ierr;

         ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq_ref, xq_ref, xfb_detJ | xfb_iJ,
                                         QuadChanged, &JData));
         if (ierr != xf_OK) return ierr;

         //use iJ to determine the direction of reference coordinates w.r.t the physical space
         if(JData[0].iJ[0]>eff && fabs(JData[0].iJ[1])<eff && fabs(JData[0].iJ[2]) < eff && JData[0].iJ[3] > eff)
         {
            //origin on node 1
            Limiter[egrp][elem].deltax = 1.0 / JData[0].iJ[0];
            Limiter[egrp][elem].deltay = 1.0 / JData[0].iJ[3];
            a[0] = 1.0; a[1] = 0.0; a[2] = 0.0; a[3] = 1.0;
            b[0] = 1.0; b[1] = 0.0; b[2] = 0.0; b[3] = 1.0;
            //Limiter[egrp][elem].signTranMatr[4] = {1.0, 0.0, 0.0, 1.0};
         }
         else if(fabs(JData[0].iJ[0])<eff && JData[0].iJ[1]>eff && JData[0].iJ[2] < -eff && fabs(JData[0].iJ[3]) < eff)
         {
            //origin on node 2
            Limiter[egrp][elem].deltax = 1.0 / fabs(JData[0].iJ[2]);
            Limiter[egrp][elem].deltay = 1.0 / JData[0].iJ[1];
            a[0] = 0.0; a[1] = 1.0; a[2] = -1.0; a[3] = 0.0;
            b[0] = 0.0; b[1] = -1.0; b[2] = 1.0; b[3] = 0.0;
         }
         else if(JData[0].iJ[0]<-eff && fabs(JData[0].iJ[1]) <eff && fabs(JData[0].iJ[2]) < eff && JData[0].iJ[3] < -eff)
         {
            //origin on node 3
            Limiter[egrp][elem].deltax = 1.0 / fabs(JData[0].iJ[0]);
            Limiter[egrp][elem].deltay = 1.0 / fabs(JData[0].iJ[3]);
            a[0] = -1.0; a[1] = 0.0; a[2] = 0.0; a[3] = -1.0;
            b[0] = -1.0; b[1] = 0.0; b[2] = 0.0; b[3] = -1.0;
         }
         else if(fabs(JData[0].iJ[0]) < eff && JData[0].iJ[1] <-eff && JData[0].iJ[2] > eff && fabs(JData[0].iJ[3]) < eff)
         {
            //origin on node 4
            Limiter[egrp][elem].deltax = 1.0 / JData[0].iJ[2];
            Limiter[egrp][elem].deltay = 1.0 / fabs(JData[0].iJ[1]);
            a[0] = 0.0; a[1] = -1.0; a[2] = 1.0; a[3] = 0.0;
            b[0] = 0.0; b[1] = 1.0; b[2] = -1.0; b[3] = 0.0;
         }
         else 
            return xf_Error(xf_OUT_OF_BOUNDS);

         //set orientation transformation matrix
         for(k=0; k<4; k++){
             Limiter[egrp][elem].signTranMatr[k] = a[k];
             Limiter[egrp][elem].signBackMatr[k] = b[k];
         }

         //number of face
         Limiter[egrp][elem].nface = Mesh->ElemGroup[egrp].nFace[elem];   

         for(face=0; face<Limiter[egrp][elem].nface; face++) {
         
            /** Interior Face **/
            if ((ibfgrp = Mesh->ElemGroup[egrp].Face[elem][face].Group) == xf_INTERIORFACE){
            
               iiface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
               IFace  = Mesh->IFace[iiface];

               //find the normal of the elem on site
               ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face,
                                           1, &QuadData, &QuadChanged));
               if (ierr != xf_OK) return ierr;
               nq = QuadData->nquad;
               xq = QuadData->xquad;

               ierr = xf_Error(xf_IFaceNormal(Mesh, IFace, nq, xq, &NData));
               if (ierr != xf_OK) return ierr;

               for(d=0; d<dim; d++)
                  nvec[d] = NData->n[d];

               //check which side elem is on
               ierr = xf_Error(xf_IsElemOnLeft(IFace, egrp, elem, &IamL));
               if(!IamL) 
               {
                  //reverse the normal direction
                  for(d=0; d<dim; d++)
                     nvec[d] *= (-1.0);
               }

               //use normals to find the location w.r.t. the elem on site
               grpIndx = ((IamL) ? IFace.ElemGroupR : IFace.ElemGroupL);
               elmIndx = ((IamL) ? IFace.ElemR      : IFace.ElemL);
               BCtype  = 0;

               if(fabs(nvec[0])<eff&&nvec[1]>eff)        //up
               {   
                   Limiter[egrp][elem].adjEgrp[2] = grpIndx;  
                   Limiter[egrp][elem].adjEelm[2] = elmIndx;  
                   Limiter[egrp][elem].adjEtyp[2] = BCtype;
               }
               else if(fabs(nvec[0])<eff&&nvec[1]<-eff)   //down
               {   
                   Limiter[egrp][elem].adjEgrp[3] = grpIndx;
                   Limiter[egrp][elem].adjEelm[3] = elmIndx; 
                   Limiter[egrp][elem].adjEtyp[3] = BCtype;
               }
               else if(fabs(nvec[1])<eff&&nvec[0]>eff)   //right
               {    
                   Limiter[egrp][elem].adjEgrp[1] = grpIndx; 
                   Limiter[egrp][elem].adjEelm[1] = elmIndx; 
                   Limiter[egrp][elem].adjEtyp[1] = BCtype; 
               }
               else if(fabs(nvec[1])<eff&&nvec[0]<-eff)   //left
               {   
                   Limiter[egrp][elem].adjEgrp[0] = grpIndx;
                   Limiter[egrp][elem].adjEelm[0] = elmIndx;  
                   Limiter[egrp][elem].adjEtyp[0] = BCtype; 
               }
               else
                  return xf_Error(xf_NOT_SUPPORTED);
             
               //free allocation
               ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
               if (ierr != xf_OK) return ierr;
               ierr = xf_Error(xf_DestroyNormalData(NData));
               if (ierr != xf_OK) return ierr;

               QuadData = NULL;
               NData    = NULL;
            }
            else if(ibfgrp >= 0)     /** Boundary Face **/
            {
               
               ibface = Mesh->ElemGroup[egrp].Face[elem][face].Number;
               BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];

               // skip zero-measure boundary groups
               if (strncmp(Mesh->BFaceGroup[ibfgrp].Title, "ZeroMeasure", 11) != 0){

               //find the normal first
               ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, face, 
                                           1, &QuadData, &QuadChanged));
               if (ierr != xf_OK) return ierr;

               nq = QuadData->nquad;
               xq = QuadData->xquad;

               ierr = xf_Error(xf_BFaceNormal(Mesh, BFace, nq, xq, &NData, NULL));
               if (ierr != xf_OK) return ierr;

               for(d=0; d<dim; d++)
                  nvec[d] = NData->n[d];

               //go to find the boundary condition type
               for(i=0; i<Model->nBCs; i++)
                  if(strcmp(Mesh->BFaceGroup[ibfgrp].Title, Model->nameBCs[i]) == 0)
                  { localBCflag = i; break; }
               if(i==Model->nBCs) {
                  printf("Boundary name does not match model defining.\n");
                  xf_Error(xf_BOUNDARY_CONDITION_ERROR);
               }

               if(Model->typeBCs[localBCflag] == fullBC)        BCtype = 2;      //supersonic inflow
               else if(Model->typeBCs[localBCflag] == noneBC)   BCtype = 3;      //supersonic outflow
               else if(Model->typeBCs[localBCflag] == staticPBC)   BCtype = 4;   //subsonic outflow
               else if(Model->typeBCs[localBCflag] == slipwallBC || Model->typeBCs[localBCflag] == noslipwallBC) BCtype = 1;    //slip wall
               else   return xf_Error(xf_NOT_SUPPORTED);

               //find the location w.r.t the elem on site
               grpIndx = -1;
               elmIndx = -1;
             
               if(fabs(nvec[0])<eff&&nvec[1]>eff)        //up
               {   
                   Limiter[egrp][elem].adjEgrp[2] = grpIndx;  
                   Limiter[egrp][elem].adjEelm[2] = elmIndx;  
                   Limiter[egrp][elem].adjEtyp[2] = BCtype;
               }
               else if(fabs(nvec[0])<eff&&nvec[1]<-eff)   //down
               {   
                   Limiter[egrp][elem].adjEgrp[3] = grpIndx;
                   Limiter[egrp][elem].adjEelm[3] = elmIndx; 
                   Limiter[egrp][elem].adjEtyp[3] = BCtype;
               }
               else if(fabs(nvec[1])<eff&&nvec[0]>eff)   //right
               {   
                   Limiter[egrp][elem].adjEgrp[1] = grpIndx; 
                   Limiter[egrp][elem].adjEelm[1] = elmIndx; 
                   Limiter[egrp][elem].adjEtyp[1] = BCtype; 
               }
               else if(fabs(nvec[1])<eff&&nvec[0]<-eff)   //left
               {   
                   Limiter[egrp][elem].adjEgrp[0] = grpIndx;
                   Limiter[egrp][elem].adjEelm[0] = elmIndx;  
                   Limiter[egrp][elem].adjEtyp[0] = BCtype; 
               }
               else
                  return xf_Error(xf_NOT_SUPPORTED);
               
               //free allocation
               ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
               if (ierr != xf_OK) return ierr;
               ierr = xf_Error(xf_DestroyNormalData(NData));
               if (ierr != xf_OK) return ierr;

               QuadData = NULL;
               NData    = NULL;
            
               }
            }
            // do nothing for NULL faces
         }
            
            // free orientation stuff
            ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_DestroyJacobianData(JData));
            if (ierr != xf_OK) return ierr;

            GeomPhiData = NULL;
            JData       = NULL;
            xf_Release( (void *) xglob);
      }
   }
    
    //allocate memory for relisting
    ierr = xf_Error(xf_Alloc((void **)&Relist.RearrangedEgrg, Relist.NumElem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **)&Relist.RearrangedElem, Relist.NumElem, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    bgnIndx = 0;
    endIndx = Relist.NumElem - 1;
    
    //allocate memory for smoothness indicator
    //preparation for smooth indicator calculation
    QuadChanged = xfe_True;
    
    // get required order for quadrature
    int QuadOrder;
    ierr = xf_Error(ElemQuadOrder(Model->order, &QuadOrder));
    if (ierr != xf_OK) return ierr;
    
    // get quadrature data
    ierr = xf_Error(xf_QuadElem(Mesh, 0, 0, QuadOrder, &QuadData, &QuadChanged));
    if (ierr != xf_OK) return ierr;
    
    nq = QuadData->nquad;
    
    //allocate memory for postivity check
    ierr = xf_Error(xf_Alloc( (void **) &Relist.PostCheckState, Model->Num_negPckpnt*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc( (void **) &Relist.SmoothIndtor_wq, nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &Relist.SmoothIndtor_gu, nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &Relist.SmoothIndtor_gu2, nq, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    //allocate memory here
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
        
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
            
            //fill in the moments of self element
            ierr = xf_Error(xf_Alloc((void **)&Limiter[egrp][elem].moment, nn*sr, sizeof(real)));
            if (ierr != xf_OK) return ierr;
        
            //fill the rearranged element list
            flag = xfe_False;
            for(face=0; face<Limiter[egrp][elem].nface; face++) {
                
                //judge whether this element has neighbor in halogroup
                if(Limiter[egrp][elem].adjEgrp[face] >= Mesh->nElemGroup)
                    flag = xfe_True;
            }
            
            if(flag)
            {
                Relist.RearrangedEgrg[endIndx] = egrp;
                Relist.RearrangedElem[endIndx] = elem;
                endIndx--;
            }
            else
            {
                Relist.RearrangedEgrg[bgnIndx] = egrp;
                Relist.RearrangedElem[bgnIndx] = elem;
                bgnIndx++;
            }
            
        } //elem
    }//egrp
    
    //check consistency
    if( (bgnIndx - endIndx) != 1)
    {
        xf_printf("Error in re-arrange element for parallel limiting");
        return xf_OUT_OF_BOUNDS;
    }
    else
        Relist.NumBrk = bgnIndx;

   return xf_OK;
}

int
DestroyLimiterStruct(xf_All *All, Yu_Limiter **Limiter)
{
    int ierr, i, d, egrp, elem;
    xf_Mesh *Mesh;
    
    Mesh = All->Mesh;
    
    //delete stuff for re-arrangement the element (useful for parallel)
    if(Relist.RearrangedEgrg != NULL)
        xf_Release(Relist.RearrangedEgrg);
    if(Relist.RearrangedElem != NULL)
        xf_Release(Relist.RearrangedElem);
    
    // destroy memory for smooth indicator
    xf_Release(Relist.SmoothIndtor_wq);  Relist.SmoothIndtor_wq = NULL;
    xf_Release(Relist.SmoothIndtor_gu);  Relist.SmoothIndtor_gu = NULL;
    xf_Release(Relist.SmoothIndtor_gu2); Relist.SmoothIndtor_gu2 = NULL;
    
    // destory memory for postitivity check
    xf_Release(Relist.PostCheckState);  Relist.PostCheckState = NULL;
    
    //free the allocation during limiting process
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
        
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
            
            if(Limiter[egrp][elem].moment != NULL)
                xf_Release(Limiter[egrp][elem].moment);
            
            Limiter[egrp][elem].moment = NULL;
            
        } //elem
    } //egrp
    
    //loop over element groups
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
        
        free(Limiter[egrp]);
        Limiter[egrp] = NULL;
    }
    
    free(Limiter);
    Limiter = NULL;
        
    return xf_OK;
}
static int
LocalValueFetch(real *Uval, const int sr, const int n, xf_BasisData *PhiData, const int egrp, const int elem, xf_Vector *U)
{
    int ierr;
    real *EU;
    EU = U->GenArray[egrp].rValue[elem];
    
    xf_MxM_Set(PhiData->Phi, EU, n, PhiData->nn, sr, Uval);
    
    return xf_OK;
}
static int
SmoothIndicator(real *Beta, real *TroubleStore, xf_JacobianData *JData, xf_BasisData *QuadPhiData, xf_QuadData *QuadData,
                const int nq, const int nn, const int dim_rank, const real h)
{
    int ierr, i, dim, d, iq, dim2_rank;
   
    if(dim_rank == 0)
       dim2_rank = 0;
    else if(dim_rank == 1)
       dim2_rank = 3;
    else
       xf_NOT_SUPPORTED;

    //allocate memory
    //ierr = xf_Error(xf_Alloc( (void **) &Relist.SmoothIndtor_wq, nq, sizeof(real)));
    //if (ierr != xf_OK) return ierr;
    //ierr = xf_Error(xf_Alloc( (void **) &Relist.SmoothIndtor_gu, nq, sizeof(real)));
    //if (ierr != xf_OK) return ierr;
    //ierr = xf_Error(xf_Alloc( (void **) &Relist.SmoothIndtor_gu2, nq, sizeof(real)));
    //if (ierr != xf_OK) return ierr;
       
    //get gradient
    xf_MxM_Set(QuadPhiData->GPhi+nn*nq*dim_rank,  TroubleStore, nq, nn, 1, Relist.SmoothIndtor_gu);
    xf_MxM_Set(QuadPhiData->HPhi+nn*nq*dim2_rank, TroubleStore, nq, nn, 1, Relist.SmoothIndtor_gu2);

    for (iq=0; iq<nq; iq++)
        //wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)];
        Relist.SmoothIndtor_wq[iq] = QuadData->wquad[iq];
   
    (*Beta) = 0.0;
    for (iq=0; iq<nq; iq++)
        //(*Beta) += wq[iq] * (1.0/h) * (pow(gu[iq], 2.0) + pow(gu2[iq], 2.0));
        (*Beta) += Relist.SmoothIndtor_wq[iq] * (pow(Relist.SmoothIndtor_gu[iq], 2.0) + pow(Relist.SmoothIndtor_gu2[iq], 2.0));

    // destroy memory
    //xf_Release(wq);  wq = NULL;
    //xf_Release(gu);  gu = NULL;
    //xf_Release(gu2); gu2 = NULL;

    return xf_OK;
}

int
Yu_ConductLimiting(xf_All *All, Yu_Model *Model, Yu_Limiter **Limiter, xf_Vector **pU)
{
    
    int ierr, i, j, k, d, nn, nq, list;
    int egrp, elem, sr, sr2;
    int *TroubleMarker;
    real *EU, tmp1, tmp2, scaling, *gammatmp, *Xpoly, *Ypoly;
    real BetaL, Beta, BetaR, WL, W, WR, sum, *NeightElem[4];
    real avgrho, avgu, avgv, avgp, avgrhoY[50], facrho, facp, facrhoY[50];
    enum xfe_Bool MeshIsParallel;
    enum xfe_BasisType Basis;
    real *EigVL_x, *EigVR_x, *EigVL_y, *EigVR_y;
    real *dul, *dur, *dut, *dub, *uval, *tmpp, *uvalP2;
    real *TroubleStore, *xq, psd, psd_rho, psd_u, psd_v, psd_p;
    xf_Mesh *Mesh;
    xf_Data *D;
    xf_Vector *P0U, *U, *GammaVec;

    //added Jan for Smooth Indicator Evaluation
    enum xfe_Bool QuadChanged, flag;
    int QuadOrder;
    xf_QuadData *QuadData=NULL;
    xf_JacobianData *JData=NULL;
    xf_BasisData *QuadPhiData=NULL;
    
    Mesh = All->Mesh;
    sr   = Model->nVars;
    sr2  = sr*sr;
    nn   = (Model->order+1)*(Model->order+1);
    ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &D);
    if(ierr == xf_NOT_FOUND) {xf_printf("Cannot find heat capacity ratio...\n"); return ierr;}
    else
        GammaVec = (xf_Vector *) D->Data;
    // get required order for quadrature
    ierr = xf_Error(ElemQuadOrder(Model->order, &QuadOrder));
    if (ierr != xf_OK) return ierr;
    
    //Mar 5; energy projection is implemented
    int *P, Enr_nq;
    real p;
    real *Q, *R, *Enr_qU, *Enr_qU2, *Enr_U, *Enr_xq;
    xf_BasisData *Enr_PhiData;
    xf_QuadData *Enr_QuadData;
    P       = NULL;
    Q       = NULL;
    R       = NULL;
    Enr_PhiData   = NULL;
    Enr_QuadData  = NULL;
    /* Pull off quad points for the element; make sure have enough
     * quad points (more than nn). */
    QuadChanged = xfe_True;
    
    ierr = xf_Error(xf_QuadElemAtLeast(Mesh, 0, 0, QuadOrder, nn, &Enr_QuadData, &QuadChanged));
    if (ierr != xf_OK) return ierr;
    Enr_nq = Enr_QuadData->nquad;
    Enr_xq = Enr_QuadData->xquad;

    ierr = xf_Error(xf_EvalBasis(xfe_QuadLagrange, Model->order, QuadChanged, Enr_nq, Enr_xq, xfb_Phi, &Enr_PhiData));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(xf_Alloc( (void **) &Enr_qU, Enr_nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_Alloc( (void **) &Enr_qU2, Enr_nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    /************do parallel communication right now************/
    U = (*pU);
    MeshIsParallel = (Mesh->ParallelInfo != NULL);
    if(MeshIsParallel){
        ierr = xf_Error(xf_HaloExchangeVectorBegin(U));
        if (ierr != xf_OK) return ierr;
        ierr = xf_Error(xf_HaloExchangeVectorBegin(GammaVec));
        if (ierr != xf_OK) return ierr;
    }
    
    /***********************************************************/
    //preparation for smooth indicator calculation
    QuadChanged = xfe_True;
    
    // get quadrature data
    ierr = xf_Error(xf_QuadElem(Mesh, 0, 0, QuadOrder, &QuadData, &QuadChanged));
    if (ierr != xf_OK) return ierr;
    
    nq = QuadData->nquad;
    xq = QuadData->xquad;
    
    // compute basis functions (and grads) if quad or basis or order changed
    ierr = xf_Error(xf_EvalBasis(xfe_QuadLagrange, Model->order, QuadChanged, nq, xq,
                                 xfb_Phi | xfb_GPhi | xfb_gPhi | xfb_HPhi, &QuadPhiData));
    if (ierr != xf_OK) return ierr;
    
    // compute geometry Jacobians
    ierr = xf_Error(xf_ElemJacobian(Mesh, 0, 0, nq, xq, xfb_detJ | xfb_iJ,
                                    QuadChanged, &JData));
    if (ierr != xf_OK) return ierr;

    // convert reference basis grads (GPhi) to physical grads, gPhi
    //ierr = xf_Error(xf_EvalPhysicalGrad(QuadPhiData, JData));
    //if (ierr != xf_OK) return ierr;

    /***********************************************************/
    
    
    
    // allocate for eigenvectors
    ierr = xf_Error(xf_Alloc((void **) &EigVL_x, sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &EigVR_x, sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &EigVL_y, sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &EigVR_y, sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &Xpoly, nn*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &Ypoly, nn*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    
    // allocate for temporal momery
    ierr = xf_Error(xf_Alloc((void **) &dul, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &dur, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &dut, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &dub, sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &uval,  2*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &uvalP2,  2*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &TroubleMarker,  sr, sizeof(int)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &TroubleStore, nn, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    
    
    for(i=0; i<9; i++)
        if(strcmp(Model->basis, xfe_BasisName[i])==0)
        {  Basis = i; break;}
    
    if(i == 9) xf_Error(xf_OUT_OF_BOUNDS);
     
    //fill in the data in Limiter Struct for preparation
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
        
        for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
            
            //convert to full-order basis on physical space (note:support P1 and P2 now)
            Evaluate_Derivative_Coeff(Limiter[egrp][elem].moment, egrp, elem, U, Model->order, sr);
            
        }//elem
    }//egrp
  
    //return xf_NOT_SUPPORTED;
    //only allocate once is enough.
    for(i=0; i<4; i++)
    {
       ierr = xf_Error(xf_Alloc((void **)&NeightElem[i], nn*sr, sizeof(real)));
       if (ierr != xf_OK) return ierr;
    }

    //loop again to fill in the data of neighbors
    //for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
    
    //    for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
    for (list=0; list<Relist.NumElem; list++)
    {
        //fill in corresponding group and element index
        egrp = Relist.RearrangedEgrg[list];
        elem = Relist.RearrangedElem[list];
        
        if(list == (Relist.NumBrk-1))
        {
            //check whether finish the parallel communication
            if(MeshIsParallel){
                ierr = xf_Error(xf_HaloExchangeVectorEnd(U));
                if (ierr != xf_OK) return ierr;
                ierr = xf_Error(xf_HaloExchangeVectorEnd(GammaVec));
                if (ierr != xf_OK) return ierr;
            }
        }
        
            for (i=0; i<Limiter[egrp][elem].nface; i++)
            {
                if(Limiter[egrp][elem].adjEgrp[i] == -1 && Limiter[egrp][elem].adjEelm[i] == -1)
                {
                    //boundary ghost element
                    for(k=0; k<nn*sr; k++)
                        NeightElem[i][k] = Limiter[egrp][elem].moment[k];
                    
                    if(Limiter[egrp][elem].adjEtyp[i] == 1)   //wallBC
                    {
                        //left & right orientation
                        if(i == 0 || i == 1)  
                           NeightElem[i][1] *= (-1.0); 
                        //top & bottom orientation
                        else if(i == 2 || i == 3) 
                           NeightElem[i][2] *= (-1.0); 
                        else
                            return xf_Error(xf_OUT_OF_BOUNDS);
                        
                    }
                    else if(Limiter[egrp][elem].adjEtyp[i] == 2) //inletBC
                    {
                        //for(k=0; k<sr; k++)
                        //    NeightElem[i][k] = 2.0 * Model->limiterVars[k] - Limiter[egrp][elem].moment[k];
                        //set higher order moments (p>0) to zeros
                        for(k=sr; k<nn*sr; k++)
                            NeightElem[i][k] = 0.0;
                    }
                    else  //outletBC or others
                    {
                        //Neumann condition is imposed
                        //set higher order moments (p>0) to zeros
                        for(k=sr; k<nn*sr; k++) {
                            NeightElem[i][k] = 0.0;
                        }
                    }
                }
                else
                {
                    //handle variable thermal property
                    tmp1 = GammaVec->GenArray[egrp].rValue[elem][0];
                    tmp2 = GammaVec->GenArray[Limiter[egrp][elem].adjEgrp[i]].rValue[Limiter[egrp][elem].adjEelm[i]][0];
                  
                    if(Model->GammaVaryFlag && fabs(tmp1 - tmp2) > Model->Gammathreshold)
                    {
                       // find interpolate state for least square projection
                       xf_MxM_Set(Enr_PhiData->Phi, U->GenArray[Limiter[egrp][elem].adjEgrp[i]].rValue[Limiter[egrp][elem].adjEelm[i]], Enr_nq, nn, sr, Enr_qU);

                       //specify how energy will change
                       for(k=0; k<Enr_nq; k++)
                       {
                          Enr_U = Enr_qU + k*sr;
                          p   = (tmp2 - 1.0)*(Enr_U[3] - 0.5*(Enr_U[1]*Enr_U[1] + Enr_U[2]*Enr_U[2])/Enr_U[0]);
                          Enr_U[3] = p/(tmp1 - 1.0) + 0.5*(Enr_U[1]*Enr_U[1] + Enr_U[2]*Enr_U[2])/Enr_U[0];
                    
                       }

                       EU = U->GenArray[Limiter[egrp][elem].adjEgrp[i]].rValue[Limiter[egrp][elem].adjEelm[i]];
                       for(k=0; k<nn*sr; k++)
                          Enr_qU2[k] = EU[k];

                       //ierr = xf_Error(xf_ProjectOnElemLeastSquares_OnlyEnergy(Mesh, 0, 0, sr, Enr_QuadData, xfe_QuadLagrange, Model->order, QuadChanged,
                       //                                                        &Q, &R, &P, Enr_qU, Enr_qU2));
                       ierr = xf_Error(xf_ProjectOnElemQR_OnlyEnergy(Mesh, 0, 0, sr, Enr_QuadData, xfe_QuadLagrange, Model->order, QuadChanged,
                                                                               &Q, &R, &P, Enr_qU, Enr_qU2));
                       if (ierr != xf_OK) return ierr;

                       Evaluate_Derivative_Coeff_VersionII(NeightElem[i], Enr_qU2, Model->order, sr);

                    }
                    else
                    {
                       //interior element
                       //not in halo group
                       if(Limiter[egrp][elem].adjEgrp[i] < Mesh->nElemGroup)
                       {
                          for(k=0; k<nn*sr; k++)
                              NeightElem[i][k] = Limiter[Limiter[egrp][elem].adjEgrp[i]][Limiter[egrp][elem].adjEelm[i]].moment[k];
                       }
                       else                       
                       //in halo group (indeed parallel)
                       {
                           Evaluate_Derivative_Coeff(NeightElem[i], Limiter[egrp][elem].adjEgrp[i], Limiter[egrp][elem].adjEelm[i], U, Model->order, sr);
                       }
                    }

               /*     //handle variable thermal property
                    if(Model->GammaVaryFlag && fabs(tmp1 - tmp2) > Model->Gammathreshold)
                    {
                       psd_rho = NeightElem[i][0];
                       psd_u   = NeightElem[i][1]/NeightElem[i][0];
                       psd_v   = NeightElem[i][2]/NeightElem[i][0];
                       
                       //reconstruction a kind of psuedo state for energy at neighbor
                       for(k=0; k<nn; k++)
                       {
                          psd   = - 0.5*(psd_u*psd_u + psd_v*psd_v)*NeightElem[i][k*sr] + psd_u*NeightElem[i][k*sr+1]
                                  + psd_v*NeightElem[i][k*sr+2];
                          psd_p = (tmp2 - 1.0) * (NeightElem[i][k*sr+3] - psd);
                          
                          NeightElem[i][k*sr+3] = psd_p/(tmp1 - 1.0) + psd;
                       }
                    }
               */     
                }
                
            }//nfac

            //whether check negative value
            flag = xfe_False;
            //diagnolize Flux Jacobian
            gammatmp = GammaVec->GenArray[egrp].rValue[elem];
            ierr = xf_Error(FillinCharTransMatrix(EigVL_x, EigVR_x, EigVL_y, EigVR_y, Model, gammatmp[0], Limiter[egrp][elem].moment));
            if (ierr != xf_OK) return ierr;
            
            //fill Xpoly and Ypoly before go limiting
            for(k=0; k<nn*sr; k++)
            {
                Xpoly[k] = Limiter[egrp][elem].moment[k];
                Ypoly[k] = Limiter[egrp][elem].moment[k];
            }
            
            //Prepare for Trouble detection
            for(k=0; k<sr; k++){
                dul[k] = Limiter[egrp][elem].moment[k] - NeightElem[0][k];
                dur[k] = NeightElem[1][k] - Limiter[egrp][elem].moment[k];
                dut[k] = NeightElem[2][k] - Limiter[egrp][elem].moment[k];
                dub[k] = Limiter[egrp][elem].moment[k] - NeightElem[3][k];
                
                //left-right
                uval[k] = Limiter[egrp][elem].moment[k+sr];
                //top-bottonm
                uval[sr+k] = Limiter[egrp][elem].moment[k+2*sr];
                
                if(Model->order == 2)
                {
                    uvalP2[k]    = 1.5*Limiter[egrp][elem].moment[k+4*sr];
                    uvalP2[sr+k] = 1.5*Limiter[egrp][elem].moment[k+5*sr];
                }
            }
           
            tmpp = uval;
            LeftMulti(EigVL_x, tmpp,  sr);
            LeftMulti(EigVL_x, dul,   sr);
            LeftMulti(EigVL_x, dur,   sr);
            
            
            tmpp = uval + sr;
            LeftMulti(EigVL_y, tmpp,  sr);
            LeftMulti(EigVL_y, dut,   sr);
            LeftMulti(EigVL_y, dub,   sr);
        
            if(Model->order == 2)
            {
                tmpp = uvalP2;
                LeftMulti(EigVL_x, tmpp,  sr);
                tmpp = uvalP2 + sr;
                LeftMulti(EigVL_y, tmpp,  sr);
            }
            
            //do x-direction
            //evalute at point (0, 0.5) and (1, 0.5)
            //Trouble Detection
            for(k=0; k<sr; k++)
            {
                TroubleMarker[k] = 1;
                if(fabs(uval[k]) > 0.5*fabs(dul[k] - dur[k]))
                   tmp1 = minmod(uval[k], dul[k], dur[k]);
                else
                    tmp1 = uval[k];
                
                if(fabs(tmp1 - uval[k])<thredhold)
                    TroubleMarker[k] = 0;

                    //additional check for higher-order
                    if(Model->order == 2 && 
                       fabs(uvalP2[k]) > 2.0*tmp1) 
                    flag = xfe_True;

                //additional constrain for positivity check
                tmp1 = minmod(uval[k], dul[k], dur[k]);
                if(fabs(tmp1 - uval[k])>thredhold)
                    flag = xfe_True;

            }
            
            
            if(WhetherInTrouble(TroubleMarker,sr) == 1)
            {
                //will check negative value
                flag = xfe_True;
                
                //transform self, left and right neighbors to Characteristic space
                for(j=0; j<nn; j++)
                {
                    tmpp = Xpoly;
                    tmpp = tmpp + j*sr;
                    LeftMulti(EigVL_x, tmpp, sr);
                    
                    tmpp = NeightElem[0];
                    tmpp = tmpp + j*sr;
                    LeftMulti(EigVL_x, tmpp, sr);
                    
                    tmpp = NeightElem[1];
                    tmpp = tmpp + j*sr;
                    LeftMulti(EigVL_x, tmpp, sr);
                }
                
                //adjust neighbors to hold the mean of self
                for(k=0; k<sr; k++)
                {
                    NeightElem[0][k] = Xpoly[k];
                    NeightElem[1][k] = Xpoly[k];
                }
                
                //computation the smoothness indicator
                for(k=0; k<sr; k++)
                if(TroubleMarker[k]==1)
                {
                    PutToTroubleStore(NeightElem[0], TroubleStore, Model->order, sr, k);
                    SmoothIndicator(&BetaL, TroubleStore, JData, QuadPhiData, QuadData,
                                    nq, nn, 0, Limiter[egrp][elem].deltax);
                    //BetaL = BetaL/Limiter[egrp][elem].deltay;
                    
                    PutToTroubleStore(NeightElem[1], TroubleStore, Model->order, sr, k);
                    SmoothIndicator(&BetaR, TroubleStore, JData, QuadPhiData, QuadData,
                                    nq, nn, 0, Limiter[egrp][elem].deltax);
                    //BetaR = BetaR/Limiter[egrp][elem].deltay;
                    
                    PutToTroubleStore(Xpoly, TroubleStore, Model->order, sr, k);
                    SmoothIndicator(&Beta, TroubleStore, JData, QuadPhiData, QuadData,
                                    nq, nn, 0, Limiter[egrp][elem].deltax);
                    //Beta = Beta/Limiter[egrp][elem].deltay;

                    WL = 0.001/pow(1.0e-6 + BetaL, 2.0);
                    WR = 0.001/pow(1.0e-6 + BetaR, 2.0);
                    W  = 0.998/pow(1.0e-6 + Beta , 2.0);
                    sum = WL + W + WR;
                    
                    WL = WL/sum; W = W/sum; WR = WR/sum;

                    //reconstruct a smoother polynomial
                    for(j=0; j<nn; j++)
                        Xpoly[j*sr+k] =
                                        WL * NeightElem[0][j*sr+k]
                                      + W  * Xpoly[j*sr+k]
                                      + WR * NeightElem[1][j*sr+k];
                    
                }
                
                for(j=0; j<nn; j++)
                {
                    tmpp = Xpoly;
                    tmpp = tmpp + j*sr;
                    LeftMulti(EigVR_x, tmpp, sr);
                }
                
            }//x-direction finished
            
            //do y-direction
            //evaluate at point (0.5, 0) and (0.5, 1)
            //Trouble Detection
            for(k=0; k<sr; k++)
            {
                TroubleMarker[k] = 1;
                if(fabs(uval[sr+k]) > 0.5*fabs(dub[k] - dut[k]))
                tmp1 = minmod(uval[sr+k], dut[k], dub[k]);
                else
                    tmp1 = uval[sr+k];
                
                if(fabs(tmp1 - uval[sr+k])<thredhold)
                   TroubleMarker[k] = 0;
                
                   //additional check for higher-order
                   if(Model->order == 2 
                      && fabs(uvalP2[sr+k]) > 2.0*tmp1)
                    flag = xfe_True;
           
                //additional constrain for positivity check
                tmp1 = minmod(uval[sr+k], dut[k], dub[k]);
                if(fabs(tmp1 - uval[sr+k])>thredhold)
                    flag = xfe_True;
                  
            }

            if(WhetherInTrouble(TroubleMarker,sr) == 1)
            {
                //will check negative value
                flag = xfe_True;
                
                //transform self, top and bottom neighbors to Characteristic space
                for(j=0; j<nn; j++)
                {
                    tmpp = Ypoly;
                    tmpp = tmpp + j*sr;
                    LeftMulti(EigVL_y, tmpp, sr);
                    
                    tmpp = NeightElem[2];
                    tmpp = tmpp + j*sr;
                    LeftMulti(EigVL_y, tmpp, sr);
                    
                    tmpp = NeightElem[3];
                    tmpp = tmpp + j*sr;
                    LeftMulti(EigVL_y, tmpp, sr);
                }
                
                //adjust neighbors to hold the mean of self
                for(k=0; k<sr; k++)
                {
                    NeightElem[2][k] = Ypoly[k];
                    NeightElem[3][k] = Ypoly[k];
                }
                
                //computation the smoothness indicator
                for(k=0; k<sr; k++)
                if(TroubleMarker[k]==1)
                {
                    PutToTroubleStore(NeightElem[2], TroubleStore, Model->order, sr, k);
                    SmoothIndicator(&BetaL, TroubleStore, JData, QuadPhiData, QuadData,
                                    nq, nn, 1, Limiter[egrp][elem].deltay);
                    //BetaL = BetaL/Limiter[egrp][elem].deltax;
                    
                    PutToTroubleStore(NeightElem[3], TroubleStore, Model->order, sr, k);
                    SmoothIndicator(&BetaR, TroubleStore, JData, QuadPhiData, QuadData,
                                    nq, nn, 1, Limiter[egrp][elem].deltay);
                    //BetaR = BetaR/Limiter[egrp][elem].deltax;
                    
                    PutToTroubleStore(Ypoly, TroubleStore, Model->order, sr, k);
                    SmoothIndicator(&Beta, TroubleStore, JData, QuadPhiData, QuadData,
                                    nq, nn, 1, Limiter[egrp][elem].deltay);
                    //Beta = Beta/Limiter[egrp][elem].deltax;
                        
                    WL = 0.001/pow(1.0e-6 + BetaL, 2.0);
                    WR = 0.001/pow(1.0e-6 + BetaR, 2.0);
                    W  = 0.998/pow(1.0e-6 + Beta , 2.0);
                    sum = WL + W + WR;
                        
                    WL = WL/sum; W = W/sum; WR = WR/sum;
                        
                    //reconstruct a smoother polynomial
                    for(j=0; j<nn; j++)
                        Ypoly[j*sr+k] = WL * NeightElem[2][j*sr+k]
                                      + W  * Ypoly[j*sr+k]
                                      + WR * NeightElem[3][j*sr+k];
                
                }
                
                for(j=0; j<nn; j++)
                {
                    tmpp = Ypoly;
                    tmpp = tmpp + j*sr;
                    LeftMulti(EigVR_y, tmpp, sr);
                }
            }//y-direction finished
    
            //for time saving, no need to modified the original polynomial
            //we could directly go to next, if limiting is not activated.
            if(!flag)
                continue;
            
            //finally, combine X & Y direction to form limited polynomial
            for(k=0; k<nn*sr; k++)
                Limiter[egrp][elem].moment[k] = 0.5*(Xpoly[k] + Ypoly[k]);
           
            for(k=0; k<sr; k++)
                BackOut_BasisFcn_Coeff(Limiter[egrp][elem].moment, egrp, elem, U, Model->order, sr, k);
            
            //very important step for troubled element
            //positivity preserving limiter
            if(flag)
            {
                //compute mean values of density and pressure
                avgrho = Limiter[egrp][elem].moment[0];
                avgu   = Limiter[egrp][elem].moment[1]/Limiter[egrp][elem].moment[0];
                avgv   = Limiter[egrp][elem].moment[2]/Limiter[egrp][elem].moment[0];
                avgp   = (gammatmp[0] - 1.0) * (Limiter[egrp][elem].moment[3] - 0.5*avgrho*(pow(avgu, 2.0) + pow(avgv, 2.0)));
                //if(avgp < Model->PressureThredhold)
                if(avgp < 0.9*Model->LimM * pow(avgrho, gammatmp[0]))
                { 
                   printf("Solution is not TVD; suggest reduce time step! %d\n", elem); 
                   printf("avgp %lf; avgrho %lf\n", avgp, avgrho);
                   return xf_OUT_OF_BOUNDS; 
                }

                for(j=0; j<(sr-4); j++)
                    avgrhoY[j] = Limiter[egrp][elem].moment[4+j];
                
                //positivity guarantee for density
                NegativeDensityChecking(Model, gammatmp[0], egrp, elem, U, avgrho, avgrhoY, &facrho, facrhoY);

                for(k=1; k<nn; k++)
                {
                    Limiter[egrp][elem].moment[sr*k] = Limiter[egrp][elem].moment[sr*k] * facrho;
                    for(j=0; j<(sr-4); j++)
                        Limiter[egrp][elem].moment[sr*k+4+j] = Limiter[egrp][elem].moment[sr*k+4+j] * facrhoY[j];
                }
                
                for(k=0; k<sr; k++)
                    BackOut_BasisFcn_Coeff(Limiter[egrp][elem].moment, egrp, elem, U, Model->order, sr, k);
                
                //positivity guarantee for pressure
                NegativePressureChecking(Model, gammatmp[0], egrp, elem, U, avgrho, avgp, &facp);
                for(k=sr; k<nn*sr; k++)
                    Limiter[egrp][elem].moment[k] = Limiter[egrp][elem].moment[k] * facp;
                for(k=0; k<sr; k++)
                    BackOut_BasisFcn_Coeff(Limiter[egrp][elem].moment, egrp, elem, U, Model->order, sr, k);
               
            }
            
            
            
    }//list         
 //       } //elem
 //   }//egrp

            //release all the smooth indicator stuff
            ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
            if (ierr != xf_OK) return ierr;
            
            ierr = xf_Error(xf_DestroyBasisData(QuadPhiData, xfe_False));
            if (ierr != xf_OK) return ierr;
            
            ierr = xf_Error(xf_DestroyJacobianData(JData));
            if (ierr != xf_OK) return ierr;

    //release temporary memory
    xf_Release(EigVL_x); EigVL_x = NULL;
    xf_Release(EigVR_x); EigVR_x = NULL;
    xf_Release(EigVL_y); EigVL_y = NULL;
    xf_Release(EigVR_y); EigVR_y = NULL;
    xf_Release(Xpoly); Xpoly = NULL;
    xf_Release(Ypoly); Ypoly = NULL;
    xf_Release(dul); dul = NULL;
    xf_Release(dur); dur = NULL;
    xf_Release(dut); dut = NULL;
    xf_Release(dub); dub = NULL;
    xf_Release(uval); uval = NULL;
    xf_Release(uvalP2); uvalP2 = NULL;
    xf_Release(TroubleMarker); TroubleMarker = NULL;
    xf_Release(TroubleStore); TroubleStore = NULL;

    //Mar 5; Free allocation for Energy projection
    xf_Release((void *) Q);
    xf_Release((void *) R);
    xf_Release((void *) P);
    xf_Release((void *) Enr_qU);
    xf_Release((void *) Enr_qU2);

    ierr = xf_Error(xf_DestroyGenericQuadData(Enr_QuadData));
    if (ierr != xf_OK) return ierr;

    ierr = xf_Error(xf_DestroyBasisData(Enr_PhiData, xfe_True));
    if (ierr != xf_OK) return ierr;


    for(i=0; i<4; i++)
    {
       if(NeightElem[i] != NULL)
          xf_Release(NeightElem[i]);
       NeightElem[i] = NULL;
    }

   return xf_OK;
}

real minmod(real a, real b, real c)
{
   real sgn, m;
           
   if(a * b <= 0.0 || b * c <= 0.0)
      return 0.0;
               
    sgn = (a > 0.0) ? 1.0 : -1.0;
      a = fabs(a);
      b = fabs(b);
      c = fabs(c);
      m = (a < b) ? a : b;
      m = (c < m) ? c : m;
     
      return sgn * m;
}

//left multiply Matrix R(sr*sr) on U(sr); note R is column priority 
static void
LeftMulti(real *R, real *U, int sr)
{ 
   int ierr, i, j;
   real *tmpU;

   ierr = xf_Error(xf_Alloc((void **) &tmpU, sr, sizeof(real)));

   for(i=0; i<sr; i++)
      tmpU[i] = U[i];

   for(i=0; i<sr; i++)
   {
      U[i] = 0.0;
      for(j=0; j<sr; j++)
         U[i] += R[j + i*sr] * tmpU[j];
   }

   xf_Release(tmpU);
}

static int
TransposeMatrix(real *M, int sr)
{
   int ierr, i, j;
   real *tmpU;

   ierr = xf_Error(xf_Alloc((void **) &tmpU, sr*sr, sizeof(real)));
   if (ierr != xf_OK) return ierr;

   for(i=0; i<sr*sr; i++)
      tmpU[i] = M[i];
   
   for(i=0; i<sr; i++)
      for(j=0; j<sr; j++)
         M[i*sr + j] = tmpU[i + sr*j];
   
   xf_Release(tmpU); tmpU = NULL;
   
   return xf_OK;
}

int
DiagonalizeBasedonIDOLoveCFD(const real *n, real *Ln, real *Rn, real gamma, const real *U, const int sr)
{
    int i, j, nscalars;
    real l[3], qn, ql, q2, c, H;
    real rho, u, v, p, E;

    //only support 2d computation
    nscalars = sr - 4;
    
    rho = U[0];
    u   = U[1]/U[0];
    v   = U[2]/U[0];
    E   = U[3]/U[0];
    p   = (gamma - 1.0)*(U[3] - 0.5*rho*(u*u + v*v));
    
    if(rho <= 0 || p<=0 )
    {
       printf("%lf %lf\n", rho, p);
       return xf_Error(xf_NON_PHYSICAL);
    }
    
    l[0] = -n[1]; l[1] = n[0];
    qn = u*n[0] + v*n[1];
    ql = u*l[0] + v*l[1];
    q2  = u*u + v*v;
    c = sqrt(gamma * p/rho);
    H = E + p/rho;

    //first row 
    Ln[0] = 0.5*((gamma - 1.0)*q2/2.0/c/c + qn /c);
    Ln[1] = -0.5*((gamma - 1.0)*u/c/c + n[0]/c);
    Ln[2] = -0.5*((gamma - 1.0)*v/c/c + n[1]/c);
    Ln[3] = (gamma - 1.0)/2.0/c/c;
    for(i=0; i<nscalars; i++)
       Ln[4+i] = 0.0;

    //second row
    Ln[sr] = 1.0 - (gamma - 1.0)*q2/2.0/c/c;
    Ln[sr+1] = (gamma - 1.0)*u/c/c;
    Ln[sr+2] = (gamma - 1.0)*v/c/c;
    Ln[sr+3] = -(gamma - 1.0)/c/c;
    for(i=0; i<nscalars; i++)
       Ln[sr+4+i] = 0.0;

    //third row
    Ln[2*sr] = 0.5*((gamma - 1.0)*q2/2.0/c/c - qn/c);
    Ln[2*sr+1] = -0.5*((gamma - 1.0)*u/c/c - n[0]/c);
    Ln[2*sr+2] = -0.5*((gamma - 1.0)*v/c/c - n[1]/c);
    Ln[2*sr+3] = (gamma - 1.0)/2.0/c/c;
    for(i=0; i<nscalars; i++)
       Ln[2*sr+4+i] = 0.0;
    
    //forth row 
    Ln[3*sr] = - ql;
    Ln[3*sr+1] = l[0];
    Ln[3*sr+2] = l[1];
    Ln[3*sr+3] = 0.0;
    for(i=0; i<nscalars; i++)
       Ln[3*sr+4+i] = 0.0;
   
    //fifth to end rows
    for(i=0; i<nscalars; i++)
    {
       Ln[(4+i)*sr]   = -U[4+i]/U[0]/U[0];
       Ln[(4+i)*sr+1] = 0.0;
       Ln[(4+i)*sr+2] = 0.0;
       Ln[(4+i)*sr+3] = 0.0;
       for(j=0; j<nscalars; j++)
       if(i == j)
          Ln[(4+i)*sr+4+j] = 1.0/U[0];
       else
          Ln[(4+i)*sr+4+j] = 0.0;
    }

    //first row
    Rn[0] = 1.0; Rn[1] = 1.0; Rn[2] = 1.0; Rn[3] = 0.0;
    for(i=0; i<nscalars; i++)
       Rn[4+i] = 0.0;
    
    //second row
    Rn[sr] = u - c*n[0]; Rn[sr+1] = u; Rn[sr+2] = u + c*n[0]; Rn[sr+3] = l[0];
    for(i=0; i<nscalars; i++)
       Rn[sr+4+i] = 0.0;
    
    //third row
    Rn[2*sr] = v - c*n[1]; Rn[2*sr+1] = v; Rn[2*sr+2] = v + c*n[1]; Rn[2*sr+3] = l[1];
    for(i=0; i<nscalars; i++)
       Rn[2*sr+4+i] = 0.0;
    
    //forth row
    Rn[3*sr] = H - qn*c; Rn[3*sr+1] = q2/2.0; Rn[3*sr+2] = H + qn*c; Rn[3*sr+3] = ql;
    for(i=0; i<nscalars; i++)
       Rn[3*sr+4+i] = 0.0;
  
    //fifth to end rows
    for(i=0; i<nscalars; i++)
    {
       Rn[(4+i)*sr]   = U[4+i]/U[0];
       Rn[(4+i)*sr+1] = U[4+i]/U[0];
       Rn[(4+i)*sr+2] = U[4+i]/U[0];
       Rn[(4+i)*sr+3] = 0.0;
       for(j=0; j<nscalars; j++)
       if(i == j)
          Rn[(4+i)*sr+4+j] = U[0];
       else
          Rn[(4+i)*sr+4+j] = 0.0;
    }

    return xf_OK;
}

int
FillinCharTransMatrix(real *iR_x, real *R_x, real *iR_y, real *R_y, Yu_Model *Model, const real gamma, const real *U)
{
   int ierr, i, j, k, sr, sr2;
   real u, v, p, E, n[3];
   real *Fx_U, *Fy_U;

   sr = Model->nVars;
   sr2 = sr*sr;

    //diagonalize Jacobi analytically
    //x-direction
    n[0] = 1.0; n[1] = 0.0;
    DiagonalizeBasedonIDOLoveCFD(n, iR_x, R_x, gamma, U, sr);
    
    //y-direction
    n[0] = 0.0; n[1] = 1.0;
    DiagonalizeBasedonIDOLoveCFD(n, iR_y, R_y, gamma, U, sr);

   return xf_OK;
}
static int
Evaluate_Derivative_Coeff(real *c, const int egrp, const int elem, xf_Vector *U, const int order, const int sr)
{
    int k, j;
    real *EU, tmp[9];
    real tranMP1[16] = { 0.25,  0.25,  0.25,  0.25,
                        -0.25,  0.25, -0.25,  0.25,
                        -0.25, -0.25,  0.25,  0.25,
                         0.25, -0.25, -0.25,  0.25 };
    real tranMP2[81] = {1.0/36.0,  1.0/9.0,  1.0/36.0,  1.0/9.0,  4.0/9.0,  1.0/9.0, 1.0/36.0,  1.0/9.0, 1.0/36.0,
                        -1.0/12.0,  0.0,     1.0/12.0,  -1.0/3.0,  0.0,      1.0/3.0, -1.0/12.0,  0.0,    1.0/12.0,
                        1.0/18.0,  -1.0/9.0, 1.0/18.0,  2.0/9.0,  -4.0/9.0, 2.0/9.0, 1.0/18.0, -1.0/9.0, 1.0/18.0,
                        -1.0/12.0, -1.0/3.0, -1.0/12.0, 0.0,      0.0,         0.0,  1.0/12.0,  1.0/3.0, 1.0/12.0,
                        1.0/4.0,      0.0,   -1.0/4.0,  0.0,      0.0,         0.0,  -1.0/4.0,  0.0,     1.0/4.0,
                        -1.0/6.0,   1.0/3.0, -1.0/6.0,  0.0,      0.0,         0.0,  1.0/6.0,  -1.0/3.0, 1.0/6.0,
                        1.0/18.0,   2.0/9.0, 1.0/18.0,  -1.0/9.0, -4.0/9.0, -1.0/9.0, 1.0/18.0, 2.0/9.0, 1.0/18.0,
                        -1.0/6.0,     0.0,   1.0/6.0,   1.0/3.0,  0.0,      -1.0/3.0, -1.0/6.0, 0.0,     1.0/6.0,
                        1.0/9.0,   -2.0/9.0, 1.0/9.0,   -2.0/9.0, 4.0/9.0,  -2.0/9.0, 1.0/9.0, -2.0/9.0, 1.0/9.0};
    real LetinOrder[81] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]

    if(order == 1)
    for(k=0; k<sr; k++) {
    
       for(j=0; j<4; j++)
           tmp[j] = EU[j*sr+k];
            
       LeftMulti(tranMP1, tmp, 4);
            
       for(j=0; j<4; j++)
           c[j*sr+k] = tmp[j];
    }
    else if(order == 2)
    for(k=0; k<sr; k++) {
       
       for(j=0; j<9; j++)
           tmp[j] = EU[j*sr+k];

       //note: the order is not in order "00 10 20 01 11 21 02 12 22"        
       LeftMulti(tranMP2, tmp, 9);

       //let's put everything in order
       LeftMulti(LetinOrder, tmp, 9);

       for(j=0; j<9; j++)
           c[j*sr+k] = tmp[j];
    }
    else
       return xf_Error(xf_NOT_SUPPORTED);
    
    return xf_OK;
}
static int
Evaluate_Derivative_Coeff_VersionII(real *c, real *EU, const int order, const int sr)
{
    int k, j;
    real tmp[9];
    real tranMP1[16] = { 0.25,  0.25,  0.25,  0.25,
        -0.25,  0.25, -0.25,  0.25,
        -0.25, -0.25,  0.25,  0.25,
        0.25, -0.25, -0.25,  0.25 };
    real tranMP2[81] = {1.0/36.0,  1.0/9.0,  1.0/36.0,  1.0/9.0,  4.0/9.0,  1.0/9.0, 1.0/36.0,  1.0/9.0, 1.0/36.0,
        -1.0/12.0,  0.0,     1.0/12.0,  -1.0/3.0,  0.0,      1.0/3.0, -1.0/12.0,  0.0,    1.0/12.0,
        1.0/18.0,  -1.0/9.0, 1.0/18.0,  2.0/9.0,  -4.0/9.0, 2.0/9.0, 1.0/18.0, -1.0/9.0, 1.0/18.0,
        -1.0/12.0, -1.0/3.0, -1.0/12.0, 0.0,      0.0,         0.0,  1.0/12.0,  1.0/3.0, 1.0/12.0,
        1.0/4.0,      0.0,   -1.0/4.0,  0.0,      0.0,         0.0,  -1.0/4.0,  0.0,     1.0/4.0,
        -1.0/6.0,   1.0/3.0, -1.0/6.0,  0.0,      0.0,         0.0,  1.0/6.0,  -1.0/3.0, 1.0/6.0,
        1.0/18.0,   2.0/9.0, 1.0/18.0,  -1.0/9.0, -4.0/9.0, -1.0/9.0, 1.0/18.0, 2.0/9.0, 1.0/18.0,
        -1.0/6.0,     0.0,   1.0/6.0,   1.0/3.0,  0.0,      -1.0/3.0, -1.0/6.0, 0.0,     1.0/6.0,
        1.0/9.0,   -2.0/9.0, 1.0/9.0,   -2.0/9.0, 4.0/9.0,  -2.0/9.0, 1.0/9.0, -2.0/9.0, 1.0/9.0};
    real LetinOrder[81] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    
    if(order == 1)
        for(k=0; k<sr; k++) {
            
            for(j=0; j<4; j++)
                tmp[j] = EU[j*sr+k];
            
            LeftMulti(tranMP1, tmp, 4);
            
            for(j=0; j<4; j++)
                c[j*sr+k] = tmp[j];
        }
    else if(order == 2)
        for(k=0; k<sr; k++) {
            
            for(j=0; j<9; j++)
                tmp[j] = EU[j*sr+k];
            
            //note: the order is not in order "00 10 20 01 11 21 02 12 22"        
            LeftMulti(tranMP2, tmp, 9);
            
            //let's put everything in order
            LeftMulti(LetinOrder, tmp, 9);
            
            for(j=0; j<9; j++)
                c[j*sr+k] = tmp[j];
        }
    else
        return xf_Error(xf_NOT_SUPPORTED);
    
    return xf_OK;
}

static int PutToTroubleStore(real *Poly, real *TroubleStore, const int order, const int sr, const int rank)
{
    int i, j, k;
    real tmp[9];
    real tranMP1[16]={ 1.0, -1.0, -1.0,  1.0,
                       1.0,  1.0, -1.0, -1.0,
                       1.0, -1.0,  1.0, -1.0,
                       1.0,  1.0,  1.0,  1.0 };
    real tranMP2[81]={ 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0,
                       1.0,  0.0, -0.5, -1.0,  0.0,  0.5,  1.0,  0.0, -0.5,
                       1.0,  1.0,  1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,
                       1.0, -1.0,  1.0,  0.0,  0.0,  0.0, -0.5,  0.5, -0.5,
                       1.0,  0.0, -0.5,  0.0,  0.0,  0.0, -0.5,  0.0,  0.25,
                       1.0,  1.0,  1.0,  0.0,  0.0,  0.0, -0.5, -0.5, -0.5,
                       1.0, -1.0,  1.0,  1.0, -1.0,  1.0,  1.0, -1.0,  1.0,
                       1.0,  0.0, -0.5,  1.0,  0.0, -0.5,  1.0,  0.0, -0.5,
                       1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 };
    real LetinOrder[81] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
    
    i = rank;
    {
        if(order == 2){
            for(j=0; j<9; j++)
                tmp[j] = Poly[j*sr+i];
            
            //put back into original order
            LeftMulti(LetinOrder, tmp, 9);
            LeftMulti(tranMP2, tmp, 9);
            
            for(j=0; j<9; j++)
                TroubleStore[j] = tmp[j];
        }
        else if(order == 1){
            for(j=0; j<4; j++)
                tmp[j] = Poly[j*sr+i];
            
            LeftMulti(tranMP1, tmp, 4);
            
            for(j=0; j<4; j++)
                TroubleStore[j] = tmp[j];
            
        }
        else
            return xf_Error(xf_NOT_SUPPORTED);
    }
    
    return xf_OK;
}
static int
BackOut_BasisFcn_Coeff(const real *c, const int egrp, const int elem, xf_Vector *U, const int order, const int sr, const int rank)
{
    int i, j, k;
    real tmp[9];
    real tranMP1[16]={ 1.0, -1.0, -1.0,  1.0,
                       1.0,  1.0, -1.0, -1.0,
                       1.0, -1.0,  1.0, -1.0,
                       1.0,  1.0,  1.0,  1.0 };

    real tranMP2[81]={ 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0,
                       1.0,  0.0, -0.5, -1.0,  0.0,  0.5,  1.0,  0.0, -0.5,
                       1.0,  1.0,  1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,
                       1.0, -1.0,  1.0,  0.0,  0.0,  0.0, -0.5,  0.5, -0.5,
                       1.0,  0.0, -0.5,  0.0,  0.0,  0.0, -0.5,  0.0,  0.25,
                       1.0,  1.0,  1.0,  0.0,  0.0,  0.0, -0.5, -0.5, -0.5,
                       1.0, -1.0,  1.0,  1.0, -1.0,  1.0,  1.0, -1.0,  1.0,
                       1.0,  0.0, -0.5,  1.0,  0.0, -0.5,  1.0,  0.0, -0.5,
                       1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 };
    real LetinOrder[81] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
                           

    real *EU;
    EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
    
    //for(i=0; i<sr; i++) 
    i = rank;
    {
        if(order == 2){
            for(j=0; j<9; j++)
                tmp[j] = c[j*sr+i];
       
            //put back into original order
            LeftMulti(LetinOrder, tmp, 9);
            LeftMulti(tranMP2, tmp, 9);
        
            for(j=0; j<9; j++)
                EU[j*sr+i] = tmp[j];
        }
        else if(order == 1){
            for(j=0; j<4; j++)
                tmp[j] = c[j*sr+i];
            
            LeftMulti(tranMP1, tmp, 4);
            
            for(j=0; j<4; j++)
                EU[j*sr+i] = tmp[j];
           
        }
        else
            return xf_Error(xf_NOT_SUPPORTED);
    }
    
    return xf_OK;
}
static int
Finish_Limiting(const int *Limited, const int sr)
{
   int i, flag;
   
   flag = 1;
   for(i=0; i<sr; i++)
   if(Limited[i] != 1)
   { 
      flag = 0;
      break;
   }

   if(flag ==0)
      return 0;
   else
      return 1;
} 
static int
WhetherInTrouble(int TroubleMarker[], const int sr)
{ 
    int j;
    int flag;

    flag = 0;
    for(j=0; j<sr; j++)
    if(TroubleMarker[j] == 1)
    { 
        flag = 1;
        break;
    }
        
    if(flag == 1)
        return 1;
    else
        return 0;
}
static int
findminimum(const real *val, const int N, real *min)
{
    int i;
    
    (*min) = val[0];
    for(i=1; i<N; i++)
       if(val[i]<(*min))
       (*min) = val[i];
    
    return xf_OK;
}
static int
NegativeDensityChecking(Yu_Model *Model, const real gamma, const int egrp, const int elem, xf_Vector *U,
                        const real avgrho, const real *avgrhoY, real *facrho, real *facrhoY)
{
    int ierr, k, spe, flag, sr;
    real rho, u, v, p, E, min_rho, rhoY[50], min_rhoY[50];
    real *EU, *State;
    EU = U->GenArray[egrp].rValue[elem];
    sr = Model->nVars;
    
    //ierr = xf_Error(xf_Alloc( (void **) &allState, Model->Num_negPckpnt*sr, sizeof(real)));
    //if (ierr != xf_OK) return ierr;
          
    xf_MxM_Set(Model->Phi_negPckpnt->Phi, EU, Model->Num_negPckpnt, Model->Phi_negPckpnt->nn, sr, Relist.PostCheckState);
    
    min_rho = 1.0e+12;
    for(spe=0; spe<(sr-4); spe++)
    min_rhoY[spe] = 1.0e+12;

    for(k=0; k<Model->Num_negPckpnt; k++)
    {
        State = Relist.PostCheckState + k*sr;
        rho = State[0];
        if(rho<min_rho)
            min_rho = rho;       

        for(spe=0; spe<(sr-4); spe++)
        {
           rhoY[spe] = State[4+spe];
           if(rhoY[spe] < min_rhoY[spe])
              min_rhoY[spe] = rhoY[spe];
        }
    }
   
    (*facrho) = fabs((avgrho - 1.0e-13)/(avgrho - min_rho));
    if((*facrho)>1.0 || min_rho > 0.0)
        (*facrho) = 1.0;

    for(spe=0; spe<(sr-4); spe++)
    {
       facrhoY[spe] = fabs((avgrhoY[spe] - 1.0e-13)/(avgrhoY[spe] - min_rhoY[spe]));
       if(facrhoY[spe]>1.0 || min_rhoY[spe]>0.0)
           facrhoY[spe] = 1.0;
    }
    
    //xf_Release(allState);
    
    return xf_OK;
}
static int
NegativePressureChecking(Yu_Model *Model, const real gamma, const int egrp, const int elem, xf_Vector *U, const real avgrho, 
                         const real avgp, real *facp)
{
   int ierr, k, sr, fac;
   real rho, u, v, p, E;
   real *EU, *State;
   EU = U->GenArray[egrp].rValue[elem];
   sr = Model->nVars;

   //ierr = xf_Error(xf_Alloc( (void **) &allState, Model->Num_negPckpnt*sr, sizeof(real)));
   //if (ierr != xf_OK) return ierr;

   xf_MxM_Set(Model->Phi_negPckpnt->Phi, EU, Model->Num_negPckpnt, Model->Phi_negPckpnt->nn, sr, Relist.PostCheckState);

    (*facp) = 1.0;
   for(k=0; k<Model->Num_negPckpnt; k++)
   {
      State = Relist.PostCheckState + k*sr;
      rho = State[0];
      u   = State[1]/State[0];
      v   = State[2]/State[0];
      p   = (gamma - 1.0)*(State[3] - 0.5*rho*(u*u + v*v));

      //if(p<Model->PressureThredhold)
      //replace by Yu's entropy bounding idea intead of positivity preserving
      //implemented with little relaxation to avoid too much numerical dissipation
      if(p < 0.9*(Model->LimM * pow(rho, gamma))) 
      {
          //this seems a wrong expression; changed Apr 2013
          //fac = avgp/(avgp - p + Model->PressureThredhold);
          //fac = (avgp - Model->PressureThredhold)/(avgp - p);
          //changed again with entropy bound
          fac = (avgp - Model->LimM * pow(avgrho, gamma)) / (avgp - Model->LimM * pow(avgrho, gamma) + Model->LimM * pow(rho, gamma) - p);

          if(fac<(*facp))
              (*facp) = fac;
      }
   }

   //xf_Release(allState);

   return xf_OK;
}
static int
Convert_Derivative_to_PhysicsSpace(real *moment, const real *Trans_sign, const int order, const int sr)
{
    int i, j, k;
    real tmp[9];
    real dchdx = Trans_sign[0];
    real detdx = Trans_sign[2];
    real dchdy = Trans_sign[1];
    real detdy = Trans_sign[3];
    
    if(order == 1)
        for(k=0; k<sr; k++)
        {
            //partial_x
            tmp[1] = moment[k+sr]*dchdx + moment[k+2*sr]*detdx;
            //partial_y
            tmp[2] = moment[k+sr]*dchdy + moment[k+2*sr]*detdy;
            //partial_x & partial_y
            tmp[3] = moment[k+3*sr]*(dchdy*detdx + detdy*detdx);
            
            //restore
            for(j=1; j<=3; j++)
                moment[k+j*sr] = tmp[j];
        }
    else if(order == 2)
        for(k=0; k<sr; k++)
        {
            //partial_x
            tmp[1] = moment[k+sr]*dchdx + moment[k+2*sr]*detdx;
            //partial_y
            tmp[2] = moment[k+sr]*dchdy + moment[k+2*sr]*detdy;
            //partial_x & partial_y
            tmp[3] = moment[k+3*sr]*(dchdy*detdx + detdy*detdx);
            //partial_x & partial_x
            tmp[4] = moment[k+4*sr]*dchdx*dchdx + moment[k+5*sr]*detdx*detdx;
            //partial_y & partial_y
            tmp[5] = moment[k+4*sr]*dchdy*dchdy + moment[k+5*sr]*detdy*detdy;
            //partial_x & partial_x & partial_y
            tmp[6] = moment[k+6*sr]*dchdx*dchdx*detdy + moment[k+7*sr]*dchdy*detdx*detdx;
            //partial_x & partial_y & partial_y
            tmp[7] = moment[k+6*sr]*dchdy*dchdy*detdx + moment[k+7*sr]*dchdx*detdy*detdy;
            //partial_x & partial_x & partial_y & partial_y
            tmp[8] = moment[k+8*sr]*(dchdx*dchdx*detdy*detdy + dchdy*dchdy*detdx*detdx);
            
            //restore
            for(j=1; j<=8; j++)
                moment[k+j*sr] = tmp[j];
        }
    else 
        return xf_Error(xf_NOT_SUPPORTED);
    
    return xf_OK;
}
static int
Convert_Derivative_to_ReferenceSpace(real *moment, const real *Trans_sign, const int order, const int sr, const int rank)
{
    int i, j, k;
    real tmp[9];
    real dxdch = Trans_sign[0];
    real dxdet = Trans_sign[2];
    real dydch = Trans_sign[1];
    real dydet = Trans_sign[3];
    
    if(order == 1)
       // for(k=0; k<sr; k++)    
    {  
       k = rank;

            //partial_x
            tmp[1] = moment[k+sr]*dxdch + moment[k+2*sr]*dydch;
            //partial_y
            tmp[2] = moment[k+sr]*dxdet + moment[k+2*sr]*dydet;
            //partial_x & partial_y
            tmp[3] = moment[k+3*sr]*(dxdch*dydet + dxdet*dydch);
            
            //restore
            for(j=1; j<=3; j++)
                moment[k+j*sr] = tmp[j];
    }
    else if(order == 2)
        //for(k=0; k<sr; k++)
    {
       k = rank;

            //partial_ch
            tmp[1] = moment[k+sr]*dxdch + moment[k+2*sr]*dydch;
            //partial_et
            tmp[2] = moment[k+sr]*dxdet + moment[k+2*sr]*dydet;
            //partial_ch & partial_et
            tmp[3] = moment[k+3*sr]*(dxdch*dydet + dxdet*dydch);
            //partial_ch & partial_ch
            tmp[4] = moment[k+4*sr]*dxdch*dxdch + moment[k+5*sr]*dydch*dydch;
            //partial_et & partial_et
            tmp[5] = moment[k+4*sr]*dxdet*dxdet + moment[k+5*sr]*dydet*dydet;
            //partial_ch & partial_ch & partial_et
            tmp[6] = moment[k+6*sr]*dxdch*dxdch*dydet + moment[k+7*sr]*dydch*dydch*dxdet;
            //partial_et & partial_et & partial_ch
            tmp[7] = moment[k+6*sr]*dxdet*dxdet*dydch + moment[k+7*sr]*dydet*dydet*dxdch;
            //partial_ch & partial_ch & partial_et & partial_et
            tmp[8] = moment[k+8*sr]*(dxdch*dxdch*dydet*dydet + dxdet*dxdet*dydch*dydch);
        
            //restore
            for(j=1; j<=8; j++)
                moment[k+j*sr] = tmp[j];
    }
    else 
        return xf_Error(xf_NOT_SUPPORTED);
    
    return xf_OK;
}
