//------------------legacy code with novel implementation---------------------//
//this .c file includes the relevant implementation for Yu's entropy bounding
//idea for shock capturing based on DG. The background knowledge will be published.
//Author: Yu Lv
//Email: ylv@stanford.edu
//Note that the idea is expected to apply on 2D and 3D
//Current phase, it is not expected to be jointly used with double flux
//only used with the following basis function:
/*
 * Triangle: TriHierarch
 * Tetrahedron: TetHierarch
 * Quadrilateral: QuadLegendre
 * Hexahedron: HexLegendre
 */
#include "xf.h"
#include "xf_AllStruct.h"
#include "xf_MPI.h"
#include "xf_Memory.h"
#include "xf_Basis.h"
#include "xf_All.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xfYu_Model.h"
#include "xfYu_ModelStruct.h"
#include "xf_Math.h"
#include "xf_Quad.h"
#include "xf_MeshTools.h"
#include "xf_Adapt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//data structure for re-arrange element list for parallel
typedef struct{
   int NumElem;   //number of element in block of the current cpu
   int NumBrk;    //number of element without neighbor in halo
   int *RearrangedEgrg;
   int *RearrangedElem;
}
RearrangedElemList;
static RearrangedElemList Relist;


//thredhold for density positivity
#define eps 1.0e-12
static int
Refresh(xf_QuadData **QuadDataElem, xf_QuadData **QuadDataFace, xf_JacobianData **JData,
        xf_BasisData **PhiData, real **Quadx, real **QuadStat, real **QuadGrad)
{
   //default setting: maximum number of face = 6
   
   int ierr, i;

    if((*QuadDataElem) != NULL){
      ierr = xf_Error(xf_DestroyGenericQuadData(*QuadDataElem));
      if (ierr != xf_OK) return ierr;
   }

   if((*QuadDataFace) != NULL){
      ierr = xf_Error(xf_DestroyGenericQuadData(*QuadDataFace));
      if (ierr != xf_OK) return ierr;
    }

   if((*JData) != NULL){
      ierr = xf_Error(xf_DestroyJacobianData(*JData));
      if (ierr != xf_OK) return ierr;
   }

   if ((*PhiData) != NULL){
      ierr = xf_Error(xf_DestroyBasisData(*PhiData, xfe_True));
      if (ierr != xf_OK) return ierr;
   }

   xf_Release( (void *) *Quadx);
   xf_Release( (void *) *QuadStat);

   *QuadDataElem = NULL;
   *QuadDataFace = NULL;
   *JData        = NULL;
   *PhiData      = NULL;
   *Quadx        = NULL;
   *QuadStat     = NULL;

   if(QuadGrad != NULL)
   {
      xf_Release( (void *) *QuadGrad);
      *QuadGrad = NULL;
   }

   return xf_OK;
}

/*************************************************************************************/
//evalute the minimum length scale 
int
Yu_MinFaceLengthScale(xf_All *All, Yu_Model *Model, xf_Vector *U)
{
   int ierr, i, j, k, iq;
   int elem, egrp, sr, Order, pOrder, QuadOrder;
   int nq, nn, dim, nface, endIndx, bgnIndx;
   int Neigh_egrp, Neigh_elem, Neigh_face;
   real *EU, *MinFaceLen;
   real *xq, Velem, FaceJ, MaxFS;
   enum xfe_ShapeType Shape;
   enum xfe_Bool QuadChanged, flag;
   enum xfe_BasisType Basis;
   xf_QuadData *QuadDataElem, *QuadDataFace;
   xf_JacobianData *JData;
   xf_NormalData *NData[6];      //only consider six-face element at maximum
   xf_Mesh *Mesh;

   //////////////////////
   Mesh = All->Mesh;
   sr   = Model->nVars;
   dim  = Mesh->Dim;

   //////////////////////
   QuadDataFace = NULL;
   QuadDataElem = NULL;
   JData        = NULL;
   for(i=0; i<6; i++)
      NData[i] = NULL;

   if(Model->EntropyBdFlag)
      Relist.NumElem = 0;
   
   //loopping over all elements
   for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
 
         EU = U->GenArray[egrp].rValue[elem]; //state vector on element
         Order = xf_InterpOrder(U, egrp, elem);
         
         //quadrature points inside element
         ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order, &QuadOrder));
         if (ierr != xf_OK) return ierr;
     
         ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadDataElem, &QuadChanged));
         if (ierr != xf_OK) return ierr;
         
         //element volume Jacobian information
         QuadChanged = xfe_False;
         ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, QuadDataElem->nquad,
                                         QuadDataElem->xquad, xfb_detJ, QuadChanged, &JData));
         if (ierr != xf_OK) return ierr;
        
         //compute the element volume
         Velem = 0.0;
         for(k=0; k<QuadDataElem->nquad; k++)
            Velem += QuadDataElem->wquad[k]*JData->detJ[k*(JData->nq!=1)];
     
         //quadrature points on element faces (the same on one element)
         Shape = QuadDataElem->Shape;
         ierr = xf_Error(xf_Shape2nFace(Shape, &nface));
         if (ierr != xf_OK) return ierr;

         ierr = xf_Error(xf_GetQuadOrderGeneralFace(Mesh, NULL, egrp, Order, &QuadOrder));
         if (ierr != xf_OK) return ierr;

         //find the maximum length scale on element face
         MaxFS = 0.0;
         for (i=0; i<nface; i++){
             
            //access the quadrature rule for each face
            ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, i, QuadOrder, &QuadDataFace,
                                        &QuadChanged));
            if (ierr != xf_OK) return ierr;
             
            nq = QuadDataFace->nquad;
            xq = QuadDataFace->xquad;
             
            ierr = xf_Error(xf_ElemNormal(Mesh, egrp, elem, i, 0, nq, xq, &NData[i]));
            if (ierr != xf_OK) return ierr;
          
            FaceJ = 0.0;
            for(iq=0; iq<NData[i]->nq; iq++)
             {
               for(k=0; k<NData[i]->dim; k++)
                  FaceJ += pow(NData[i]->n[iq*(NData[i]->dim)*(NData[i]->nq!=1)+k], 2.0);
               FaceJ = sqrt(FaceJ);

               if(FaceJ < 1.0e-10) 
                  xf_printf("Warning, the face length is close to ZERO!\n");

               if(MaxFS < FaceJ)
                  MaxFS = FaceJ;
            }

         }

         //make record for minimum length scale for time step evaluation
         MinFaceLen = Model->MinFaceLen->GenArray[egrp].rValue[elem];
         MinFaceLen[0] = Velem / MaxFS;

         //renew all the allocation
         if(QuadDataElem != NULL){
            ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataElem));
            if (ierr != xf_OK) return ierr;
         }

         if(QuadDataFace != NULL){
            ierr = xf_Error(xf_DestroyGenericQuadData(QuadDataFace));
            if (ierr != xf_OK) return ierr;
         }
   
         if(JData != NULL){
            ierr = xf_Error(xf_DestroyJacobianData(JData));
            if (ierr != xf_OK) return ierr;
         }

         for(i=0; i<6; i++)
             if(NData[i] != NULL){
                ierr = xf_Error(xf_DestroyNormalData(NData[i]));
                if (ierr != xf_OK) return ierr;
                NData[i] = NULL;
             }

         QuadDataFace = NULL;
         QuadDataElem = NULL;
         JData        = NULL;
         
         if(Model->EntropyBdFlag)
            Relist.NumElem ++;
      }//elem
   }//egrp
   
   //re-shuffle the element list
   if(Model->EntropyBdFlag){
      Relist.RearrangedEgrg = NULL;
      Relist.RearrangedElem = NULL;
      
      ierr = xf_Error(xf_Alloc((void **)&Relist.RearrangedEgrg, Relist.NumElem, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_Alloc((void **)&Relist.RearrangedElem, Relist.NumElem, sizeof(int)));
      if (ierr != xf_OK) return ierr;
      bgnIndx = 0;  endIndx = Relist.NumElem - 1;
      
      //loop again in order to fill Relist
      for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
         for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
            
            flag = xfe_False;
            for(i = 0; i < Mesh->ElemGroup[egrp].nFace[elem]; i++)
            {
               ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, i, &Neigh_egrp, &Neigh_elem, &Neigh_face));
               if(ierr != xf_OK) return ierr;
               
               //if face is no boundary; never mind
               if(Neigh_egrp==-1 || Neigh_elem==-1 ||  Neigh_face==-1)
                  continue;
               
               //adjacent to halo interface
               if(Neigh_egrp >= Mesh->nElemGroup)
               {
                  //printf("%d %d --- %d %d\n", egrp, elem, Neigh_egrp, Neigh_elem);
                  //xf_printf("...%d %d --- %d %d\n", egrp, elem, Neigh_egrp, Neigh_elem);
                  flag = xfe_True;break;}
            }//nface
            
            if(flag){
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
            
         }//egrp & elem
     
      //check consistency
      if( (bgnIndx - endIndx) != 1) return xf_OUT_OF_BOUNDS;
      else Relist.NumBrk = bgnIndx;
   }

   return xf_OK;
}
int DestroyEntropyBoundStruct()
{
   if(Relist.RearrangedEgrg != NULL)
      xf_Release(Relist.RearrangedEgrg);
   if(Relist.RearrangedElem != NULL)
      xf_Release(Relist.RearrangedElem);
   
   return xf_OK;
}

/*************************************************************************************/
static const real s0 = 30.0; // accounts for s entropy offset (to ensure positivity)
//function: compute entropy value given state U
static void evalEntropy(real *U, real gamma, int dim, real *entropy)
{
   real d, p;

   d = U[0];
   if(dim == 2)
      p = (gamma-1.) * (U[3] - 0.5 * (U[1]*U[1] + U[2]*U[2]) / U[0]);
   
   if(dim == 3)
      p = (gamma-1.) * (U[4] - 0.5 * (U[1]*U[1] + U[2]*U[2] + U[3]*U[3]) / U[0]);

   (*entropy) = d * (log(p) - gamma * log(d) + s0);

}




/*************************************************************************************/
//function: compute divergence of entropy flux given state U and gradient gU
//note: \p (u rho s) / \p (x) + \p (v rho s) / \p (y) + \p (w rho s) / \p (z)
/*
\p (s) = \p (log（p) - gamma*log(rho))
1/p \p p = (gamma-1.0)/p \p (rhoE - 0.5* (rho u)^2 / rho - 0.5* (rho v)^2 / rho)
*/
static void linearizeEntropy(real *qU, real *qgU, xf_QuadData *QuadDataElem, xf_JacobianData *JData, 
                             real gamma, int dim, const int sr, real *entropy, real *diverg,
                             real *pmin, real *pmax, Yu_Model *Model, real *DenGrad, real *VorxMag)
{
   int i, j, k, nq, indx1, indx2, indx3;
   real coeff[10], *U, gU[50*3];
   real tmp = 0.0, rho, p, u, v, w;
   real u2, c_log_p, wq;
   real du_dx, dv_dx, dw_dx;
   real du_dy, dv_dy, dw_dy;
   real du_dz, dv_dz, dw_dz;
   real tau12, tau23, tau31, nablau, tau_nablau;

   //setting 
   nq = QuadDataElem->nquad; 

   *entropy = 0.0;
   *diverg = 0.0;
   if(pmin!=NULL){
   *pmin = 1.e+16;
   }
   if(pmax!=NULL){
   *pmax = -1.e+16;
   }

   if(DenGrad != NULL){
      *DenGrad = 0.;
   }

   if(VorxMag != NULL){
      *VorxMag = 0.;
   }
   for(i=0; i<nq; i++)
   {
        U = qU + i*sr;
        for(j=0; j<dim; j++)
           for(k=0; k<sr; k++)
              gU[j*sr+k] = qgU[j*nq*sr + i*sr + k];

        //account for quadrature weight
        wq = QuadDataElem->wquad[i]*JData->detJ[i*(JData->nq!=1)];

   rho = U[0];
      
   u = U[1]/U[0]; v = U[2]/U[0];
      
   if(dim == 3)
      w = U[3]/U[0];
   else
      w = 0.0;

   u2 = u*u + v*v + w*w;
   p = (gamma-1.) * (U[dim+1] - 0.5 * rho * u2);

   if(pmin!=NULL && p<(*pmin)) *pmin = p;
   if(pmax!=NULL && p>(*pmax)) *pmax = p;
   
   if(dim == 2){
      // x-direction
      c_log_p = (gamma-1.) * rho * u / p;
      coeff[0] = c_log_p * u2 /2. - gamma * u;
      coeff[1] = (log(p) - gamma*log(rho) + s0 - c_log_p * u);
      coeff[2] = - c_log_p * v;
      coeff[3] = c_log_p;
      for(j=0; j<4; j++)
         *diverg += coeff[j] * gU[0*sr + j] * wq;
      // y-direction
      c_log_p = (gamma-1.) * rho * v / p;
      coeff[0] = c_log_p * u2 /2. - gamma * v;
      coeff[1] = - c_log_p * u;
      coeff[2] = (log(p) - gamma*log(rho) + s0 - c_log_p * v);
      coeff[3] = c_log_p;

      for(j=0; j<4; j++)
         *diverg += coeff[j] * gU[1*sr + j] * wq;

      if(Model->DiffFlag && !Model->AVmodel)
      {
         indx1 = 0*sr + 0; indx2 = 1*sr + 0;
         du_dx = (gU[0*sr + 1] - u * gU[indx1]) / rho;
         du_dy = (gU[1*sr + 1] - u * gU[indx2]) / rho;
         dv_dx = (gU[0*sr + 2] - v * gU[indx1]) / rho;
         dv_dy = (gU[1*sr + 2] - v * gU[indx2]) / rho;
         tau12 = (du_dy + dv_dx);
         nablau = du_dx+dv_dy;
         tau_nablau = tau12*tau12 + du_dx*du_dx + dv_dy*dv_dy - 2./3.*nablau*nablau;
         tau_nablau *= (Model->mu_c);
         *diverg -= (gamma-1.) * rho / p * tau_nablau * wq;
      }

      if(DenGrad != NULL)
         *DenGrad += wq * sqrt(gU[0*sr+0]*gU[0*sr+0] + gU[1*sr+0]*gU[1*sr+0]);
      if(VorxMag != NULL)
      {
         indx1 = 0*sr + 0; indx2 = 1*sr + 0;
         dv_dx = (gU[0*sr + 2] - v * gU[indx1]) / rho;
         du_dy = (gU[1*sr + 1] - u * gU[indx2]) / rho;
         *VorxMag += wq * fabs(dv_dx - du_dy);
      }
   }
   if(dim == 3){
      // x-direction
      c_log_p = (gamma-1.) * rho * u / p;
      coeff[0] = c_log_p * u2 /2. - gamma * u;
      coeff[1] = (log(p) - gamma*log(rho) + s0 - c_log_p * u);
      coeff[2] = - c_log_p * v;
      coeff[3] = - c_log_p * w;
      coeff[4] = c_log_p;
      for(j=0; j<5; j++)
         *diverg += coeff[j] * gU[0*sr + j] * wq;
      // y-direction
      c_log_p = (gamma-1.) * rho * v / p;
      coeff[0] = c_log_p * u2 /2. - gamma * v;
      coeff[1] = -c_log_p * u;
      coeff[2] = (log(p) - gamma*log(rho) + s0 - c_log_p * v);
      coeff[3] = -c_log_p * w;
      coeff[4] = c_log_p; 
      for(j=0; j<5; j++)
         *diverg += coeff[j] * gU[1*sr + j] * wq;
      // z-direction
      c_log_p = (gamma-1.) * rho * w / p;
      coeff[0] = c_log_p * u2 /2. - gamma * w;
      coeff[1] = c_log_p * u;
      coeff[2] = c_log_p * v;
      coeff[3] = (log(p) - gamma*log(rho) + s0 - c_log_p * w);
      coeff[4] = c_log_p;
      for(j=0; j<5; j++)
         *diverg += coeff[j] * gU[2*sr + j] * wq;

      if(Model->DiffFlag && !Model->AVmodel)
      {
         indx1=0*sr + 0; indx2=1*sr + 0; indx3=2*sr + 0;
         du_dx = (gU[0*sr + 1] - u * gU[indx1]) / rho;
         du_dy = (gU[1*sr + 1] - u * gU[indx2]) / rho;
         du_dz = (gU[2*sr + 1] - u * gU[indx3]) / rho;
         dv_dx = (gU[0*sr + 2] - v * gU[indx1]) / rho;
         dv_dy = (gU[1*sr + 2] - v * gU[indx2]) / rho;
         dv_dz = (gU[2*sr + 2] - v * gU[indx3]) / rho;
         dw_dx = (gU[0*sr + 3] - w * gU[indx1]) / rho;
         dw_dy = (gU[1*sr + 3] - w * gU[indx2]) / rho;
         dw_dz = (gU[2*sr + 3] - w * gU[indx3]) / rho;

         tau12 = du_dy+dv_dx; tau23 = dv_dz+dw_dy; tau31=dw_dx+du_dz;
         nablau = du_dx + dv_dy + dw_dz;
         tau_nablau = tau12*tau12 + tau23*tau23 + tau31*tau31 + du_dx*du_dx
                    + dv_dy*dv_dy + dw_dz*dw_dz - 2./3.*nablau* nablau;

         tau_nablau *= (Model->mu_c);
         *diverg -= (gamma-1.) * rho / p * tau_nablau * wq;

         //for 3D channel flow
         //if(Model->MMSSource)
         //   *diverg += (gamma-1.) * rho * u / p * wq;
      }
      
      if(DenGrad != NULL)
         *DenGrad += wq * sqrt(gU[0*sr+0]*gU[0*sr+0] + gU[1*sr+0]*gU[1*sr+0] + gU[2*sr+0]*gU[2*sr+0]);
      
      if(VorxMag != NULL)
      {
         indx1 = 0*sr + 0; indx2 = 1*sr + 0; indx3 = 2*sr + 0;
         dv_dx = (gU[0*sr + 2] - v * gU[indx1]) / rho;
         du_dy = (gU[1*sr + 1] - u * gU[indx2]) / rho;
         
         dw_dy = (gU[1*sr + 3] - w * gU[indx2]) / rho;
         dv_dz = (gU[2*sr + 2] - v * gU[indx3]) / rho;
         
         du_dz = (gU[2*sr + 1] - u * gU[indx3]) / rho;
         dw_dx = (gU[0*sr + 3] - w * gU[indx1]) / rho;
         tau12 = (dv_dx - du_dy) * ( dv_dx - du_dy);
         tau23 = (dw_dy - dv_dz) * ( dw_dy - dv_dz);
         tau31 = (du_dz - dw_dx) * ( du_dz - dw_dx);
         *VorxMag += wq * sqrt(tau12 + tau23 + tau31);
      }
   
   }

   *entropy += rho * (log(p) - gamma*log(rho) + s0) * wq;
   
   //weighted by quadrature weights
   //(*entropy) *= wq;
   //(*diverg)  *= wq; 
   }//i
}

/*************************************************************************************/
// subroutine: evaluate element-wise artificial viscosity
static int 
//Yu_EvaluateArtViscs(const int nq, const int nn, enum xfe_BasisType Basis, const real *U,
//                    const real *Phidata, const real *wq, xf_JacobianData *JData,
//                    const real l, Yu_Model *Model, real *AV)
Yu_EvaluateArtViscs(enum xfe_Bool InitFlag, const real l, const real wmax, Yu_Model *Model, real *AV, real *Decay_eff)
{
   int i, j, k, ierr, sr, Order, dim;
   real int_num, int_denom, so, se, eo, kappa;
   real alp, c1, c2, c3, c4;
   //real *Utrunc, *Utrunc_quad, *U_quad;

   sr = Model->nVars;
   Order = Model->order;
   dim = Model->dim;

   //clear up stored stuff
   for(i=1; i<sr+1; i++)
      AV[i] = 0.0;

//   ierr = xf_Error(xf_Alloc((void **) &Utrunc, nn*sr, sizeof(real)));
//   if (ierr != xf_OK) return ierr;
//   ierr = xf_Error(xf_Alloc((void **) &Utrunc_quad, nq*sr, sizeof(real)));
//   if (ierr != xf_OK) return ierr;
//   ierr = xf_Error(xf_Alloc((void **) &U_quad, nq*sr, sizeof(real)));
//   if (ierr != xf_OK) return ierr;

 //  switch(Basis){
 //     case xfe_TriHierarch: 
 //        k = (Order*Order+Order)/2;
 //        break;
 //     default:
 //        k = nn - 1;
 //        break;
 //  }

   if (Order == 0)
   {   c2 = 0.0; c4 = 1.0; alp = 1.0;}
   else if (Order == 1)
   {alp = 0.5;
    c1 = 36.0; c2 = 6.0; c4=20.5;}
   else if (Order == 2)
   {alp = 0.167; alp *= 0.8;
    c1 = 147.3; c2 = 11.8; c4= 74.0;}
   else if (Order == 3)
   {alp = 0.123; alp *= 0.8;
    c1 = 427.0;  c2 = 19.1; c4 = 173.0;}
   else if (Order == 4)
   {alp = 0.073; alp *= 0.8;
    c1 = 1000.14; c2 = 27.8; c4 = 362.3;}
   else
      return xf_NOT_SUPPORTED;
   
   //c3 = alp * c2 * 0.5/ c4 / (real) dim;
   c3 = alp * c2 * 0.5/ c4;
   *Decay_eff =  c3 * c4/alp;  
   //*Decay_eff = 20.0 * c3 * c4;  

   //truncated polynomial 
   //for (i=0; i<nn*sr; i++)
   //   Utrunc[i] = 0.2;

   //for (i=k; i<nn; i++)
   //   for (j=0; j<sr; j++)
   //      Utrunc[i*sr+j] = U[i*sr+j];

   //compute inner-product using quadrature
   //xf_MxM_Set(Phidata, U     , nq, nn, sr, U_quad);
   //xf_MxM_Set(Phidata, Utrunc, nq, nn, sr, Utrunc_quad);
 
   for (j=0; j<sr; j++) {
   //   int_num = 0.0;
   //   int_denom = 0.0;

   //   for(i=0; i<nq; i++){
   //      k = i*sr + j;
   //      int_num  += Utrunc_quad[k]*Utrunc_quad[k]*wq[i]*JData->detJ[i*(JData->nq!=1)];
   //      int_denom+= U_quad[k]*U_quad[k]*wq[i]*JData->detJ[i*(JData->nq!=1)];
   //   }//i

   //   if ( fabs(int_denom) < 1e-10 )
   //      se = -1000.0;
   //   else
   //      se = log10(int_num/int_denom);

      if (Model->dt_size < 1.0e-10)
         eo = 0.0;
      else
        // eo = c3 * l * l / Model->dt_size *4.0;
         //eo = c3 * l * wmax /alp;
         eo = c3 * l /alp;

    //  so    = log10(1.0/ (Order*Order*Order*Order));
    //  kappa = 10.0*fabs(so);

    //  if ( se < so-kappa )
    //     *(AV+1+j) = 0.0;
    //  else if ( se <= so+kappa )
    //     *(AV+1+j) = eo*0.5*(1.0 + sin(M_PI*(se-so)*0.5/kappa));
    //  else
    //     *(AV+1+j) = eo;

   if(!InitFlag)
      *(AV+1+j) = eo;
   else
      *(AV+1+j) = 0.0;
   }//j


  // xf_Release(Utrunc);
  // xf_Release(Utrunc_quad);
  // xf_Release(U_quad);

   return xf_OK;
}

typedef struct
{
   xf_QuadData     *QuadDataElem;
   xf_QuadData     *QuadDataFace;
   xf_BasisData    *PhiData;
   xf_JacobianData *JData;
   xf_BasisData    *QuadPhiData;
   real            *Quadx;
   real            *QuadStat;
   real            *QuadGrad;
   int             pnq, pnq_entropy;
}
Yu_StaticDataElem;

static int
Yu_CreateStaticDataElem(Yu_StaticDataElem **pSD)
{
   int ierr;
   Yu_StaticDataElem *SD;
   
   ierr = xf_Error(xf_Alloc( (void **) pSD, 1, sizeof(Yu_StaticDataElem)));
   if (ierr != xf_OK) return ierr;
   
   SD = (*pSD);
   SD->QuadDataElem = NULL;
   SD->QuadDataFace = NULL;
   SD->JData        = NULL;
   SD->PhiData      = NULL;
   SD->QuadPhiData  = NULL;
   SD->Quadx        = NULL;
   SD->QuadStat     = NULL;
   SD->QuadGrad     = NULL;
   SD->pnq          = -1;
   SD->pnq_entropy  = -1;
   return xf_OK;
}

static int
Yu_DestroyStaticDataElem(Yu_StaticDataElem *SD)
{
   int ierr, i;
   
   if (SD == NULL) return xf_OK;
   
   //if((*QuadDataElem) != NULL){
      ierr = xf_Error(xf_DestroyGenericQuadData(SD->QuadDataElem));
      if (ierr != xf_OK) return ierr;
   //}
   
   //if((*QuadDataFace) != NULL){
      ierr = xf_Error(xf_DestroyGenericQuadData(SD->QuadDataFace));
      if (ierr != xf_OK) return ierr;
   //}
   
   //if((*JData) != NULL){
      ierr = xf_Error(xf_DestroyJacobianData(SD->JData));
      if (ierr != xf_OK) return ierr;
   //}
   
   //if ((*PhiData) != NULL){
      ierr = xf_Error(xf_DestroyBasisData(SD->PhiData, xfe_True));
      if (ierr != xf_OK) return ierr;
   //}
   
   //if (QuadPhiData != NULL){
      ierr = xf_Error(xf_DestroyBasisData(SD->QuadPhiData, xfe_True));
      if (ierr != xf_OK) return ierr;
   //   QuadPhiData = NULL;
   //}
   
   xf_Release( (void *) SD->Quadx);
   xf_Release( (void *) SD->QuadStat);
   
   SD->QuadDataElem = NULL;
   SD->QuadDataFace = NULL;
   SD->JData        = NULL;
   SD->PhiData      = NULL;
   SD->QuadPhiData  = NULL;
   SD->Quadx        = NULL;
   SD->QuadStat     = NULL;
   
   if(SD->QuadGrad != NULL)
   {
      xf_Release( (void *) SD->QuadGrad);
   
      SD->QuadGrad = NULL;
   }
   
   
   xf_Release( (void *) SD);
   
   return xf_OK;
}
/*************************************************************************************/
//the side product of this subroutine is a robust estimate of dt
int
Yu_ConductEntropyBounding(xf_All *All, Yu_Model *Model, xf_Vector **pU, enum xfe_Bool* CtrlFlags)
{
   int ierr, i, j, k, iq, ginx, list, pnq;
   int elem, egrp, sr, Order, pOrder, QuadOrder;
   int MaxOrder, pNumQuadPtn, PtnInFace[40], QuadOrderFace[40]; 
   int nface, nface2, NumQuadPtn, nq, nn, dim, *Node;
   int Neigh_egrp, Neigh_elem, Neigh_face, TotelElem;
   real EntrpBound, *QuadStat, *Quadx, *EU, *xq, *MinEntropyVec;
   real avg_d, fac_d, avg_p, pmin, pmax, fac_p, avg_u, Utmp[5], denom;
   real avg_Y, fac_Y, eps_Y=1.e-4; //if we consider one-step reaction
   real Gam, p, d, s, dtmp, *EntropyInCell, *EntropyBound;
   enum xfe_ShapeType Shape;
   enum xfe_Bool QuadChanged, LocalMinFlag, MeshIsParallel;
   enum xfe_BasisType Basis;
   xf_QuadData *QuadDataElem, *QuadDataFace;
   xf_BasisData *PhiData, *QuadPhiData;
   xf_JacobianData *JData;
   xf_Mesh *Mesh;
   xf_Vector *U, *AdaptIndicator;
   Yu_StaticDataElem *StaticData = NULL;
   ///////for AV model
   xf_Vector *MinFaceLen, *AVmodel_data;
   real *MinFaceLen_elem, *AVmodel_data_elem, *AVmodel_data_neigh;
   //enum xfe_Bool AVmodel_Flag=xfe_False;
   real AV_s, AV_fd, vol, s_ref, thres, pcut, Decay_eff, wmax;
   real *QuadGrad, fromneigh[20];
   int myRank;

   //control flag for AV and p-adaptation
   real DenGrad, VorxMag;
   enum xfe_Bool InitFlag, RecdFlag, EvalFlag;

   //adaption control flag
   enum xfe_Bool p_Adapt_active;

   //xf_Data *limiting_data; 
   //xf_Vector *limiting; 


      
   //ierr = xf_FindDataByTitle(All->DataSet, "Limit_Factor", xfe_Vector, &limiting_data);
   //if(ierr == xf_NOT_FOUND)
   //{                
   //   xf_printf("Cannot find limiting indicator...\n");
   //   return ierr;
   //}    
   //else     
   //   limiting = (xf_Vector *) limiting_data->Data;

   InitFlag = CtrlFlags[0];
   RecdFlag = CtrlFlags[1];
   EvalFlag = CtrlFlags[2];

   p_Adapt_active = Model->Dyn_p_Adapt && (InitFlag || RecdFlag || EvalFlag);

   ierr = xf_Error(xf_MPI_GetRank(&myRank, NULL));
   if(ierr != xf_OK) return ierr;

   //////////////////////
   Mesh = All->Mesh;
   sr   = Model->nVars;
   dim  = Mesh->Dim;
   Gam  = Model->GammaInit;
   U    = *pU;
   EntrpBound = Model->LimM;
   MeshIsParallel = (Mesh->ParallelInfo != NULL);
   
   if(Model->AVmodel){
      MinFaceLen = Model->MinFaceLen;
      AVmodel_data = Model->AVmodel_data;
    
      //clear up the previous data
      //ierr = xf_Error(xf_SetZeroVector(AVmodel_data));
      //if(ierr != xf_OK) return ierr;
   }
   
   //////////////////////
   Quadx = NULL;
   QuadStat = NULL;
   EntropyInCell = NULL;
   EntropyBound  = NULL;

   /////halo exchange to ensure element
   /////entropy change of the neighbors; only use if entropy bound is not constant 
   if(!Model->ConstEntropyBdFlag && MeshIsParallel){
      ierr = xf_Error(xf_HaloExchangeVectorBegin(Model->EntropyVec));
      if(ierr != xf_OK) return ierr;
   }

   //find error indicator for adaptation
   if(Model->Stat_h_Adapt){
      ierr = xf_Error(Yu_CreateOrFindAdaptIndicator(All, xfe_False, 1, &AdaptIndicator));
      if (ierr != xf_OK) return ierr;
   }

   if(p_Adapt_active){
      ierr = xf_Error(Yu_CreateOrFindAdaptIndicator(All, xfe_False, NUM_VAR_DYM_P, &AdaptIndicator));
      if (ierr != xf_OK) return ierr;
   }
   

    //allocation for minimum entropy
    TotelElem = 0;
    for (egrp=0; egrp<Mesh->nElemGroup; egrp++)
        TotelElem += Mesh->ElemGroup[egrp].nElem;
    ierr = xf_Error(xf_Alloc((void **) &EntropyInCell, TotelElem, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc((void **) &EntropyBound, TotelElem, sizeof(real)));
    if (ierr != xf_OK) return ierr;


   //clear the AVmodel flag in order to accommandate the adjacent effect
   //for(list = 0; list < Relist.NumElem; list++){
      
      //entries
      //egrp = Relist.RearrangedEgrg[list];
      //elem = Relist.RearrangedElem[list];
      
            
      //if(Model->AVmodel){
      //   AVmodel_data_elem = AVmodel_data->GenArray[egrp].rValue[elem];
      //   AVmodel_data_elem[0] = 0.0;
      //}

   //}//list

   //create tempory memory allocation
   Yu_CreateStaticDataElem(&StaticData);
   ierr = xf_Error(xf_Alloc((void **) &Quadx, 600*dim, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc((void **) &QuadStat,600*sr, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   
   
   //first, loop over all quadrature points to find beta
   pOrder = -1;
   pNumQuadPtn = 0;
   ginx = 0;
   //for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
   //   for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
   for(list = 0; list < Relist.NumElem; list++){{
      
      //entries
      egrp = Relist.RearrangedEgrg[list];
      elem = Relist.RearrangedElem[list];
      
      //check if halo exchange is finished;
      if(!Model->ConstEntropyBdFlag && list >= (Relist.NumBrk-1) && MeshIsParallel){
         ierr = xf_Error(xf_HaloExchangeVectorEnd(Model->EntropyVec));
         if (ierr != xf_OK) return ierr;
      }
         //find proper basis type
         Basis = U->Basis[0];

         EU = U->GenArray[egrp].rValue[elem]; //state vector on element
         Order = xf_InterpOrder(U, egrp, elem);
         //for first-order FVM, the entropy principle is automatically satisfied
         if(Order == 0 && !p_Adapt_active) continue;

         //logic is not general; need to be revised????
  //   if(Order != pOrder )
        {
/*           ////refresh memory if not NULL
           if(Model->AVmodel || p_Adapt_active){
              ierr = xf_Error(Refresh(&QuadDataElem, &QuadDataFace, &JData, &PhiData, &Quadx, &QuadStat, &QuadGrad));
              if (ierr != xf_OK) return ierr;
           }
           else{
              ierr = xf_Error(Refresh(&QuadDataElem, &QuadDataFace, &JData, &PhiData, &Quadx, &QuadStat, NULL));
              if (ierr != xf_OK) return ierr;
           }
   //      ierr = xf_Error(Refresh(&QuadDataElem, &QuadDataFace, &JData, &PhiData, &Quadx, &QuadStat));
   //      if (ierr != xf_OK) return ierr;
        
         if (QuadPhiData != NULL){
            ierr = xf_Error(xf_DestroyBasisData(QuadPhiData, xfe_True));
            if (ierr != xf_OK) return ierr;
            QuadPhiData = NULL;
         }
      */      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //following way is used to avoid re-allocation
            QuadDataElem = StaticData->QuadDataElem;
            QuadDataFace = StaticData->QuadDataFace;
            JData        = StaticData->JData;
            PhiData      = StaticData->PhiData;
            QuadPhiData  = StaticData->QuadPhiData;
           // Quadx        = StaticData->Quadx;
           // QuadStat     = StaticData->QuadStat;
            QuadGrad     = StaticData->QuadGrad;
            pnq          = StaticData->pnq;
      

         NumQuadPtn = 0;
         //quadrature points inside element
         ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, Order, &QuadOrder));
         if (ierr != xf_OK) return ierr;
     
         QuadChanged = xfe_True;
         ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadDataElem, &QuadChanged));
         if (ierr != xf_OK) return ierr;

         NumQuadPtn += QuadDataElem->nquad;

         //quadrature points on element faces (the same on one element)
         Shape = QuadDataElem->Shape;
         ierr = xf_Error(xf_Shape2nFace(Shape, &nface));
         if (ierr != xf_OK) return ierr;

        
         nface2 = Mesh->ElemGroup[egrp].nFace[elem];
         
         //looping over all faces to find proper quadrature rule
  /*
         for(i=0; i<nface2; i++){
             //first, what is the neighbor
             ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, i, &Neigh_egrp, &Neigh_elem, &Neigh_face));
             if(ierr != xf_OK) return ierr;
             
             if(Neigh_egrp==-1 || Neigh_elem==-1 ||  Neigh_face==-1)
             {Neigh_egrp = egrp; Neigh_elem = elem;}

             //second, dispatch proper quadrature rule
             MaxOrder = xf_InterpOrder(U, Neigh_egrp, Neigh_elem);
             if(Order > MaxOrder)
                 MaxOrder = Order;
             
             ierr = xf_Error(xf_GetQuadOrderAcrossFace(Mesh, NULL, egrp, Neigh_egrp, MaxOrder, &QuadOrder));
             if (ierr != xf_OK) return ierr;
            
             QuadChanged = xfe_True;
             ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, i, QuadOrder, &QuadDataFace,
                                         &QuadChanged));
             if (ierr != xf_OK) return ierr;
             
             nq = QuadDataFace->nquad;
             NumQuadPtn += nq;
             PtnInFace[i] = nq;
             QuadOrderFace[i] = QuadOrder;
         }
      
         ////memory allocation
        // if(NumQuadPtn > pNumQuadPtn)
        if(NumQuadPtn > pnq){
             //ierr = xf_Error(xf_ReAlloc((void **) &Quadx, NumQuadPtn*dim, sizeof(real)));
             //if (ierr != xf_OK) return ierr;
             //ierr = xf_Error(xf_ReAlloc((void **) &QuadStat, NumQuadPtn*sr, sizeof(real)));
             //if (ierr != xf_OK) return ierr;
             ierr = xf_Error(xf_ReAlloc((void **) &Quadx, NumQuadPtn*dim, sizeof(real)));
             if (ierr != xf_OK) return ierr;
             ierr = xf_Error(xf_ReAlloc((void **) &QuadStat, NumQuadPtn*sr, sizeof(real)));
             if (ierr != xf_OK) return ierr;
        }
  */
      
           if(Model->AVmodel || p_Adapt_active)
           {
              if(QuadDataElem->nquad > StaticData->pnq_entropy)
              //request gradiant information
              ierr = xf_Error(xf_ReAlloc((void **) &QuadGrad, sr*dim*QuadDataElem->nquad, sizeof(real)));
              if (ierr != xf_OK) return ierr;
            }
   
         //coordinates of all quadrature points inside element
         xq = QuadDataElem->xquad;
         for(k=0; k<dim*QuadDataElem->nquad; k++){
            Quadx[k] = xq[k];
         }

         j = QuadDataElem->nquad;
          for (i=0; i<nface2; i++){
             //first, what is the neighbor
             ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, i, &Neigh_egrp, &Neigh_elem, &Neigh_face));
             if(ierr != xf_OK) return ierr;
             
             if(Neigh_egrp==-1 || Neigh_elem==-1 ||  Neigh_face==-1)
             {Neigh_egrp = egrp; Neigh_elem = elem;}
             
             //second, dispatch proper quadrature rule
             MaxOrder = xf_InterpOrder(U, Neigh_egrp, Neigh_elem);
             if(Order > MaxOrder)
                MaxOrder = Order;
             
             ierr = xf_Error(xf_GetQuadOrderAcrossFace(Mesh, NULL, egrp, Neigh_egrp, MaxOrder, &QuadOrder));
             if (ierr != xf_OK) return ierr;
             
             QuadChanged = xfe_True;
            ierr = xf_Error(xf_QuadFace(Mesh, egrp, elem, i, QuadOrder, &QuadDataFace,
                                         &QuadChanged));
             if (ierr != xf_OK) return ierr;

             nq = QuadDataFace->nquad;
             xq = QuadDataFace->xquad;
             
             
             NumQuadPtn += nq;
             
            //fetch basis fcn for quadrature points on faces
            ////face orientation should not matter (set to zero)
            ierr = xf_Error(xf_RefFace2Interpol(Mesh, egrp, elem, i, 0, nq, xq,
                                               Quadx+dim*j));
            if (ierr != xf_OK) return ierr;
             
             j += nq;
          }//
          //check if momery overflow
           if(NumQuadPtn > 600) return xf_OUT_OF_BOUNDS;
           
         //state vector of all quadrature points
           QuadChanged = xfe_True;
         ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, NumQuadPtn, Quadx, 
                                      xfb_Phi, &PhiData));
         if (ierr != xf_OK) return ierr;

         nn = PhiData->nn;
         
         /////////////////////////////////////
         pOrder = Order;
         pNumQuadPtn = NumQuadPtn;

      }//if(pOrder != Order)

         //!!!Jacobian informaton has to be updated for different elements
         //QuadChanged = xfe_False;
         if(Model->AVmodel || p_Adapt_active){
            QuadChanged = xfe_True;
           ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, QuadDataElem->nquad,
                                           QuadDataElem->xquad, xfb_detJ | xfb_iJ, QuadChanged, &JData));
           if (ierr != xf_OK) return ierr;
         }
         else{
           ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, QuadDataElem->nquad,
                                           QuadDataElem->xquad, xfb_detJ, QuadChanged, &JData));
           if (ierr != xf_OK) return ierr;
         }

         //take state values at all quadature points
         xf_MxM_Set(PhiData->Phi, EU, NumQuadPtn, nn, sr, QuadStat);

         if(Model->AVmodel){
         
            AVmodel_data_elem = AVmodel_data->GenArray[egrp].rValue[elem];
            if(InitFlag)
            {
               AVmodel_data_elem[0] = 0.0; 
               AVmodel_data_elem[sr+1] = 0.0; 
               AVmodel_data_elem[sr+2] = 0.0;
            }
         }
         
         //Apply Yu's entropy bounding idea
         //For initialization purpose, initial entropy has to be provided.
         if(InitFlag)
         {
            //initialization
            EntrpBound = Model->LimM;

            //if dynamic_p_adaptation is invoked 
            if(Model->Dyn_p_Adapt)
            {
               AdaptIndicator->GenArray[egrp].rValue[elem][0] = 0.0; //simultaneous entropy residual
               AdaptIndicator->GenArray[egrp].rValue[elem][1] = 0.0; //stored entropy value
               AdaptIndicator->GenArray[egrp].rValue[elem][2] = 0.0; //stored entropy flux term (currently only consider convection)
               AdaptIndicator->GenArray[egrp].rValue[elem][3] = 1.0e+16; //stored entropy flux term (currently only consider convection)
            }
         }
         else
              EntrpBound = Model->EntropyVec->GenArray[egrp].rValue[elem][0];

         ////1. compute mean density and pressure with nonlinear element
         denom = 0.;
         for (j=0; j<dim+2; j++)
            Utmp[j] = 0.;
         //~~~~one-step reaction
         if(Model->ChemSource) 
            avg_Y = 0;
         for (i=0; i<QuadDataElem->nquad; i++){
            dtmp = QuadDataElem->wquad[i]*JData->detJ[i*(JData->nq!=1)];
            denom += dtmp;
            for (j=0; j<dim+2; j++)
               Utmp[j] += dtmp * QuadStat[i*sr+j]; 
            if(Model->ChemSource) 
               avg_Y += dtmp * QuadStat[i*sr+dim+2];
         }//i
         vol = denom;
     
            
            //averaged state variable 
            for (j=0; j<dim+2; j++)
               Utmp[j] = Utmp[j]/denom;
          
            //if do averaging in this way, it violates the conservation law
            //for (j=0; j<dim+2; j++)
            //   Utmp[j] = EU[j];

            //~~~~one-step reaction
            //averaged rhoY
            if(Model->ChemSource) 
            {
               avg_Y = avg_Y/denom;
               //positivity check
               //we should reduce time-step; or just skip
               /*
               if(avg_Y<-eps_Y) 
               {
                 // printf("sr %d %lf %lf \n",sr, avg_Y, Utmp[0]);
                  for (i=0; i<QuadDataElem->nquad; i++)
                     printf("%.12lf\n", QuadStat[i*sr+dim+2]);
                  getchar();
                  return xf_NON_PHYSICAL;
               }*/
            }

            //averaged density
            avg_d = Utmp[0];
            dtmp = Utmp[1]*Utmp[1]+Utmp[2]*Utmp[2];
            if(dim == 3)
               dtmp += Utmp[3]*Utmp[3];
            //averaged pressure
            avg_p = (Gam - 1.)*(Utmp[dim+1] - 0.5*dtmp/avg_d);
 
            avg_u = sqrt(dtmp/avg_d/avg_d);

            //indicator for use of Riemann solvers
            if(avg_p< 1.2 * (Model->LimM) * pow(avg_d, Gam))
               Model->RiemannIndicator->GenArray[egrp].rValue[elem][0] = 1.;
            else
               Model->RiemannIndicator->GenArray[egrp].rValue[elem][0] = 0.;


            //////logistic check
            if(avg_d<0. || avg_p<EntrpBound * pow(avg_d, Gam)-1.0e-6){
               Node = Mesh->ElemGroup[egrp].Node[elem];
               printf("\nerror element location %lf, %lf\n", Mesh->Coord[Node[0]][0],
                         Mesh->Coord[Node[0]][1]);
               printf("elem %d, density=%.12lf, pressure=%.12lf,Minentropy=%.12lf\n", elem, avg_d, avg_p, EntrpBound);
               printf("mean density or pressure goes to non-positive values!\n");
               
               //play a trick before quit; do not give up
               dtmp = EntrpBound;
               for(i=0; i<nface; i++)
               {
                  ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, i, &Neigh_egrp, &Neigh_elem, &Neigh_face));
                  if(ierr != xf_OK) return ierr;
                  
                  //if face is no boundary; never mind
                  if(Neigh_egrp==-1 || Neigh_elem==-1 ||  Neigh_face==-1)
                     continue;
                  
                  MinEntropyVec = Model->EntropyVec->GenArray[Neigh_egrp].rValue[Neigh_elem];
                  //printf("neight entropy bound %lf\n", MinEntropyVec[0]);
                  if(MinEntropyVec[0] < dtmp)
                     dtmp = MinEntropyVec[0];
                  
               }

               if(dtmp >= EntrpBound || dtmp > avg_p / pow(avg_d, Gam))
                 if(dtmp < Model->LimM)
                  return xf_NON_PHYSICAL;
                 else
                  dtmp = Model->LimM; 
               else
                  EntrpBound = avg_p / pow(avg_d, Gam) - 1.0e-6;
            }

            
            /////////////////////////////
            //if AV model is activiated//
            /////////////////////////////
            if(Model->AVmodel){
               //state vector of interior quadrature points
               //note: the number of points and point coordinates are consistent to xf_ElemJacobian
               ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, QuadDataElem->nquad, 
                                            QuadDataElem->xquad, xfb_Phi | xfb_GPhi | xfb_gPhi, &QuadPhiData));
               if (ierr != xf_OK) return ierr;
      
               // convert reference basis grads (GPhi) to physical grads, gPhi
               ierr = xf_Error(xf_EvalPhysicalGrad(QuadPhiData, JData));
               if (ierr != xf_OK) return ierr;

               //MinFaceLen_elem   = MinFaceLen->GenArray[egrp].rValue[elem];
               AVmodel_data_elem = AVmodel_data->GenArray[egrp].rValue[elem];

               if(avg_p < 1.e-8) avg_p = 1.e-8; 
               if(avg_d < 1.e-8) avg_d = 1.e-8;
               wmax = sqrt(Gam*avg_p / avg_d) + avg_u;
               //set the amount of viscosity; if initialization, set to zero
               //ierr = xf_Error(Yu_EvaluateArtViscs(QuadDataElem->nquad, nn, Basis, EU, QuadPhiData->Phi,
               //                QuadDataElem->wquad, JData, MinFaceLen_elem[0], Model, AVmodel_data_elem));
               //if (ierr != xf_OK) return ierr;
               if(EvalFlag){
               //ierr = xf_Error(Yu_EvaluateArtViscs(InitFlag, MinFaceLen_elem[0], wmax, Model, AVmodel_data_elem, &Decay_eff));
               ierr = xf_Error(Yu_EvaluateArtViscs(InitFlag, pow(vol*2.0, 1./(real)dim), wmax, Model, AVmodel_data_elem, &Decay_eff));
               if (ierr != xf_OK) return ierr;
               }

               //for DT
               if(!InitFlag){  //avoid invoking for initialization
               dtmp = avg_p;
               //pcut = avg_p;
               for(i=0; i<nface; i++)
               {
                  ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, i, &Neigh_egrp, &Neigh_elem, &Neigh_face));
                  if(ierr != xf_OK) return ierr;
                  
                  //if face is no boundary; never mind
                  if(Neigh_egrp==-1 || Neigh_elem==-1 ||  Neigh_face==-1)
                     continue;
                  
                  AVmodel_data_neigh = Model->AVmodel_data->GenArray[Neigh_egrp].rValue[Neigh_elem];
                  //printf("neight entropy bound %lf\n", MinEntropyVec[0]);
                  if(AVmodel_data_neigh[sr+4] > dtmp)
                     dtmp = AVmodel_data_neigh[sr+4];
                  //if(fabs(AVmodel_data_neigh[sr+4] - avg_p) > fabs(pcut - avg_p))
                  //   pcut = AVmodel_data_neigh[sr+4];
               }
               }

               //tolerate threshold
               //AVmodel_data_elem[sr+4] = avg_p;
               AV_s = (1.05*1.05 - 1.) * 2. * Gam /(Gam+1.) * avg_p + avg_p;
               AV_fd = 2./(Gam-1.) * (AV_s/avg_p - 1.) / (AV_s/avg_p + (Gam+1.)/(Gam-1.));
               thres = AV_fd * (1.05*sqrt(Gam*avg_p / avg_d) + avg_u);

               wmax = 1.1*sqrt(Gam*avg_p / avg_d) + avg_u;
               //threshold with neighbor info
               if(!InitFlag){
//               pcut = dtmp;
               AV_s = 2./(Gam-1.) * (dtmp/avg_p - 1.) / (dtmp/avg_p + (Gam+1.)/(Gam-1.));
               AV_fd = sqrt(Gam * avg_p / avg_d) * sqrt(1.+ (Gam+1.)/2./Gam  * (dtmp - avg_p)/avg_p);
               if(AV_fd + avg_u > wmax)
                  wmax = AV_fd + avg_u;
               
               AV_fd = AV_s * (AV_fd + avg_u);
               if(AV_fd > thres) 
                  thres = AV_fd;
               }

               //!!!!it is risky to conduct shock detection here since negative density and pressure 
               //have not been ruled out.
            }

         ////2. ensure positivity of density
         fac_d = 1.;
         fac_Y = 1.;
         for (i=0; i<NumQuadPtn; i++){
            dtmp = fabs((avg_d - eps)/(avg_d - QuadStat[i*sr+0]));
               if(dtmp < fac_d)
                  fac_d = dtmp;
            
            //~~~~one-step reaction
            if(Model->ChemSource && avg_Y > eps_Y){
               dtmp = fabs((avg_Y - eps_Y)/(avg_Y - QuadStat[i*sr+dim+2]));
               if(dtmp < fac_Y)
                  fac_Y = dtmp;
            }
            else
               fac_Y = 0.; //if -eps <= avg_Y <= eps; set to uniform profile
         }//i


         //set detector to zero 
         if(Model->AVmodel && EvalFlag)
            AVmodel_data_elem[0] = 0.0;

         if(fac_d < 1.) {
            //////Legendre basis used; modify density state
            //!!!!conservation has to be preserved
            //rank of density variable is at first
            switch(Basis){
             case xfe_TriHierarch:
             case xfe_TetHierarch:
                for(i=0; i<dim+1; i++)
                   EU[i*sr+0] = (1.-fac_d)*avg_d + fac_d*EU[i*sr+0];
                for(i=dim+1; i<nn; i++)
                   EU[i*sr+0] *= fac_d;
                break;

             case xfe_QuadLegendre:
             case xfe_HexLegendre:
                EU[0*sr+0] = (1.-fac_d)*avg_d + fac_d*EU[0*sr+0];
                for (i=1; i<nn; i++)
                   EU[i*sr+0] *= fac_d;
                break;

             default:
               return xf_Error(xf_UNKNOWN_BASIS);
               break;
            }

            if(Model->Stat_h_Adapt)
               AdaptIndicator->GenArray[egrp].rValue[elem][0] = 1.0;

            //xf_printf("Bounding takes place!\n");
            xf_MxM_Set(PhiData->Phi, EU, NumQuadPtn, nn, sr, QuadStat);
         
           if(Model->AVmodel && EvalFlag)
              AVmodel_data_elem[0] = 1.0;
         }

         //~~~~one-step reaction
         if(Model->ChemSource && fac_Y < 1.){
            switch(Basis){
               case xfe_TriHierarch:
               case xfe_TetHierarch:
                   for(i=0; i<dim+1; i++)
                      EU[i*sr+dim+2] = (1.-fac_Y)*avg_Y + fac_Y*EU[i*sr+dim+2];
                   for(i=dim+1; i<nn; i++)
                      EU[i*sr+dim+2] *= fac_Y;
               break;

               case xfe_QuadLegendre:
               case xfe_HexLegendre:
                   EU[0*sr+dim+2] = (1.-fac_Y)*avg_Y + fac_Y*EU[0*sr+dim+2];
                   for (i=1; i<nn; i++)
                      EU[i*sr+dim+2] *= fac_Y;
               break;
            }
         }

         ////3. ensure entropy minimum principle
         fac_p = 1.;
         denom = avg_p - EntrpBound * pow(avg_d, Gam);
         if(fabs(denom) > 1.0e-8)  //relaxation
            for (i=0; i<NumQuadPtn; i++){
               dtmp = QuadStat[i*sr+1]*QuadStat[i*sr+1] + QuadStat[i*sr+2]*QuadStat[i*sr+2];
               if(dim == 3)
                  dtmp += QuadStat[i*sr+3]*QuadStat[i*sr+3];
               p = (Gam - 1.)*(QuadStat[i*sr+dim+1] - 0.5*dtmp/QuadStat[i*sr+0]);

               if(p < EntrpBound*pow(QuadStat[i*sr+0], Gam) - eps){
                  dtmp = denom/(denom + EntrpBound*pow(QuadStat[i*sr+0], Gam)-p);
                  if(dtmp < fac_p && dtmp > 0.)
                     fac_p = dtmp;
               }//if
         
            }//i

         //debug check
         //limiting->GenArray[egrp].rValue[elem][0] = 0.;
         if(fac_p < 1.)
         {   
            //xf_printf("elem %d; Pressure %lf; EntropyMin %lf\n", elem, fac_p, EntrpBound);
            Model->Num_negPckpnt=1;

            if(fac_p < 0.0)
               xf_printf("Warning, something wrong with pressure scaling!\n");
            //////modifying the state vector
            //!!!!conservation has to be preserved
            switch(Basis) {
               case xfe_QuadLegendre:
               case xfe_HexLegendre:
                  for(j=0; j<sr; j++)
                     EU[0*sr+j] = (1.-fac_p)*Utmp[j] + fac_p*EU[0*sr+j];
                  for(i=1; i<nn; i++)
                     for(j=0; j<sr; j++)
                        EU[i*sr+j] *= fac_p;
                  break;

               case xfe_TriHierarch:
               case xfe_TetHierarch:
                  for(i=0; i<dim+1; i++)
                     for(j=0; j<sr; j++)
                        EU[i*sr+j] = (1.-fac_p)*Utmp[j] + fac_p*EU[i*sr+j];
                  for(i=dim+1; i<nn; i++)
                     for(j=0; j<sr; j++)
                        EU[i*sr+j] *= fac_p;
                  break;
       
               default:
                  return xf_Error(xf_UNKNOWN_BASIS);
                  break;
            }
        
            if(Model->Stat_h_Adapt)
               AdaptIndicator->GenArray[egrp].rValue[elem][0] = 1.0;
            
           if(Model->AVmodel && EvalFlag)
              AVmodel_data_elem[0] = 1.0;
            //xf_printf("Bounding takes place!\n");
            //not necessary
            //xf_MxM_Set(PhiData->Phi, EU, NumQuadPtn, nn, sr, QuadStat);
         
            //limiting->GenArray[egrp].rValue[elem][0] = fac_p;
         }         

         //////////////////////////////////////////////////////
             //implementation for p-adaptation
             if(p_Adapt_active)
             {
                //retrieve the state and grad info
                QuadChanged = xfe_True;
                ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, QuadDataElem->nquad, 
                                             QuadDataElem->xquad, xfb_Phi | xfb_GPhi | xfb_gPhi, &QuadPhiData));
                if (ierr != xf_OK) return ierr;
      
                // convert reference basis grads (GPhi) to physical grads, gPhi
                ierr = xf_Error(xf_EvalPhysicalGrad(QuadPhiData, JData));
                if (ierr != xf_OK) return ierr;

                xf_MxM_Set(PhiData->Phi, EU, NumQuadPtn, nn, sr, QuadStat);
                for(k=0; k<dim; k++)
                   xf_MxM_Set(QuadPhiData->gPhi+nn*(QuadDataElem->nquad)*k, EU, QuadDataElem->nquad, 
                              nn, sr, QuadGrad+(QuadDataElem->nquad)*sr*k);

                //compute element-wise entropy value and flux-divergence
                //linearizeEntropy(QuadStat, QuadGrad, QuadDataElem, JData, Gam, dim, sr, &AV_s, &AV_fd, NULL, NULL, Model, &DenGrad, &VorxMag);
                linearizeEntropy(QuadStat, QuadGrad, QuadDataElem, JData, Gam, dim, sr, &AV_s, &AV_fd, NULL, NULL, Model, NULL, NULL);

                //form the entropy residual
                AV_s  /= vol; 
                AV_fd /= vol;
                //DenGrad /= vol;
                //VorxMag /= vol;

                //implementation for entropy residual
                /*if(RecdFlag && (!InitFlag))
                {
                   AdaptIndicator->GenArray[egrp].rValue[elem][1] = AV_s;
                   AdaptIndicator->GenArray[egrp].rValue[elem][2] = AV_fd;
                }

                if(EvalFlag && (!InitFlag))
                {
                   AV_s = AV_s - AdaptIndicator->GenArray[egrp].rValue[elem][1];
                   AV_fd = AV_fd + AdaptIndicator->GenArray[egrp].rValue[elem][2];

                   AdaptIndicator->GenArray[egrp].rValue[elem][0] += fabs(AV_s + Model->dt_size * 0.5 * AV_fd);
                }*/
                
                //implementation for entropy
 /*
                if(RecdFlag && (!InitFlag))
                {
                   AdaptIndicator->GenArray[egrp].rValue[elem][1] = AV_s;
                   AdaptIndicator->GenArray[egrp].rValue[elem][2] = AV_fd;
                }

                if(EvalFlag && (!InitFlag))
                {
                   AV_s = AV_s - AdaptIndicator->GenArray[egrp].rValue[elem][1];
                   AV_fd = AV_fd + AdaptIndicator->GenArray[egrp].rValue[elem][2];

                   AdaptIndicator->GenArray[egrp].rValue[elem][0] += fabs(AV_s); // + Model->dt_size * 0.5 * AV_fd);
                }
  //              */
                //implementation for density gradient
                //if(EvalFlag && (!InitFlag))
                //AdaptIndicator->GenArray[egrp].rValue[elem][0] += DenGrad;
                //implementation for vortex magnitude
                //if(EvalFlag && (!InitFlag))
                //AdaptIndicator->GenArray[egrp].rValue[elem][0] += VorxMag;

             }
             
             
             //////////////////////////////////////////////////////
             //new detection mechanism: based on the magnitude of entropy residual
             if(Model->AVmodel && (RecdFlag || EvalFlag))
             {
                xf_MxM_Set(PhiData->Phi, EU, NumQuadPtn, nn, sr, QuadStat);
                //take gradient information 
                for(k=0; k<dim; k++)
                   xf_MxM_Set(QuadPhiData->gPhi+nn*(QuadDataElem->nquad)*k, EU, QuadDataElem->nquad, 
                              nn, sr, QuadGrad+(QuadDataElem->nquad)*sr*k);

                //compute \p (rho s) / \p t + \nabla (u rho s)
                linearizeEntropy(QuadStat, QuadGrad, QuadDataElem, JData, Gam, dim, sr, &AV_s, &AV_fd, &pmin, &pmax, Model, NULL, NULL);
                
                //form the entropy residual
                AV_s  /= vol;
                AV_fd /= vol;

                if(EvalFlag){
                //s_ref = fabs(AVmodel_data_elem[sr+1]);
                s_ref = 0.5 * (s_ref + fabs(AV_s));
                //if(s_ref < fabs(AV_s)) s_ref = fabs(AV_s);

                dtmp = AV_s - AVmodel_data_elem[sr+1] + Model->dt_size * 0.5*(AV_fd + AVmodel_data_elem[sr+2]);
             
                //!!!!!! a length factor is missing!!!!!!
                //use static thresholding
                //thres = s_ref*(Model->dt_size) * 0.01; 
                
                //use DT
                thres = s_ref*(Model->dt_size) * thres * 5.0 ; /// pow(2.*vol, 1./(real)dim)/10.0 ;
                }

                //store data for next time step
                if(RecdFlag){
                AVmodel_data_elem[sr+1] = AV_s;
                AVmodel_data_elem[sr+2] = AV_fd;
                AVmodel_data_elem[sr+4] = avg_p;
                }
                if(EvalFlag){
                AVmodel_data_elem[sr+3] = fabs(dtmp)/thres; //save the present residual 

                //introduce log-scaling argument
                s_ref = 1.0;
                if(fabs(dtmp) > thres && fabs(pmax - pmin) > 0.02*avg_p)
                {
                   AVmodel_data_elem[0] = 1.0;
             //      if(log(fabs(dtmp)/thres) < Decay_eff)
             //         s_ref = log(fabs(dtmp)/thres)/Decay_eff;
                }

                //finialize viscosity; whether activiated and scaling
                for(k=1; k<sr+1; k++)
                   AVmodel_data_elem[k] *= s_ref*AVmodel_data_elem[0] * wmax; 
                   
                //linear argument; performance not good
                //   if(fabs(dtmp) < 1000.0 * (Model->dt_size))
                //     AVmodel_data_elem[k] *= fabs(dtmp) / (1000.0 * (Model->dt_size));


                }
             }
     
         //////////////////////////////////////////
         //Very important part
         //compute the Min Entropy in every element
         //adaptive Entropy calculation
         if(InitFlag != xfe_True)
         {
             
            //first, esitmate the in-cell minimum entropy (s)
             //data, avg_d in-cell entropy min; avg_p in-cell entropy max
             //dtmp disctance from min to its closest neighbor.
             xf_MxM_Set(PhiData->Phi, EU, NumQuadPtn, nn, sr, QuadStat);
             avg_d = 1.0e+30;
             avg_p =-1.0e+30;
             j = 0;
             for(i=0; i<NumQuadPtn; i++){
                 d = QuadStat[i*sr+0];
                 dtmp = QuadStat[i*sr+1]*QuadStat[i*sr+1] + QuadStat[i*sr+2]*QuadStat[i*sr+2];
                 if(dim == 3)
                     dtmp += QuadStat[i*sr+3]*QuadStat[i*sr+3];
                 p = (Gam - 1.)*(QuadStat[i*sr+dim+1] - 0.5*dtmp/QuadStat[i*sr+0]);
                 s = p/pow(d, Gam);
                 if(s<avg_d){
                     j = i;
                     avg_d = s;
                 }
                 
                 if(s>avg_p)
                     avg_p = s;

             }
                 
             dtmp = 1.0e+30;
             for(i=0; i<NumQuadPtn; i++){
                 if(j == i)
                     continue;
                 
                 denom = 0.0;
                 for(k=0; k<dim; k++){
                     denom += pow(Quadx[j*dim + k]  - Quadx[i*dim + k], 2.0);
                 }
                 
                 if(dtmp > sqrt(denom))
                     dtmp =  sqrt(denom);
             }
             //first order approximation
             avg_d = avg_d - (avg_p - avg_d)/(real) (Order +1 ); 
          
            //second, look at the neighbors for constranting (avg_p)
             avg_p = EntrpBound; 
             if(!Model->ConstEntropyBdFlag)
             for(i=0; i<nface; i++)
             {
                 ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, i, &Neigh_egrp, &Neigh_elem, &Neigh_face));
                 if(ierr != xf_OK) return ierr;
                
                 //if face is no boundary; never mind
                 if(Neigh_egrp==-1 || Neigh_elem==-1 ||  Neigh_face==-1)
                    continue;
               
                 MinEntropyVec = Model->EntropyVec->GenArray[Neigh_egrp].rValue[Neigh_elem];
                 if(MinEntropyVec[0] < avg_p)
                     avg_p = MinEntropyVec[0];
                  
             }

             //third, save data for further assignment
             EntropyInCell[ginx] = avg_d;
             EntropyBound[ginx]  = avg_p;
         }
          else //for initialization purpose
          {
              if(EntropyBound[ginx] < EntrpBound)
                  EntropyBound[ginx] = EntrpBound;
          }
          
          ginx++;

      
      StaticData->QuadDataElem = QuadDataElem;
      StaticData->QuadDataFace = QuadDataFace;
      StaticData->JData = JData;
      StaticData->PhiData = PhiData;
      StaticData->QuadPhiData = QuadPhiData;
     // StaticData->Quadx = Quadx;
     // StaticData->QuadStat = QuadStat;
      StaticData->QuadGrad = QuadGrad;
      StaticData->pnq = pnq;
      StaticData->pnq_entropy = QuadDataElem->nquad;
      }//elem
   }//egrp
 
   //trouble proprogation procedure to neighbors
   //ierr = xf_Error(xf_MPI_Allreduce(AVmodel_data->GenArray[egrp].rValue[elem], 1, xfe_SizeReal, xfe_MPI_SUM));
   //if (ierr != xf_OK) return ierr;

    //final bounded entropy assignment
   ginx=0;
   if(!InitFlag) 
   // for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
   //      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
   
//      if(MeshIsParallel && Model->AVmodel && EvalFlag){
//         ierr = xf_Error(xf_HaloExchangeVectorBegin(Model->AVmodel_data));
//         if(ierr != xf_OK) return ierr;
//      }

      for(list = 0; list < Relist.NumElem; list++){{
         
         //entries
         egrp = Relist.RearrangedEgrg[list];
         elem = Relist.RearrangedElem[list];

               
//         if(list >= (Relist.NumBrk-1) && MeshIsParallel && Model->AVmodel && EvalFlag){
//            ierr = xf_Error(xf_HaloExchangeVectorEnd(Model->AVmodel_data));
//            if (ierr != xf_OK) return ierr;
//         }
           
            MinEntropyVec= Model->EntropyVec->GenArray[egrp].rValue[elem];
            
            if(EntropyInCell[ginx] < EntropyBound[ginx])
                MinEntropyVec[0] = EntropyBound[ginx];
            else
                MinEntropyVec[0] = EntropyInCell[ginx];
 /*         
            if(Model->AVmodel && EvalFlag){
               AVmodel_data_elem = Model->AVmodel_data->GenArray[egrp].rValue[elem];
             for(i=0; i<nface; i++)
             {
                 ierr = xf_Error(xf_NeighborAcrossFace(Mesh, egrp, elem, i, &Neigh_egrp, &Neigh_elem, &Neigh_face));
                 if(ierr != xf_OK) return ierr;
                
                 //if face is no boundary; never mind
                 if(Neigh_egrp==-1 || Neigh_elem==-1 ||  Neigh_face==-1)
                 {
                      continue;
                 }
               
                 AVmodel_data_neigh = Model->AVmodel_data->GenArray[Neigh_egrp].rValue[Neigh_elem];
                 
                  if(AVmodel_data_neigh[0]>1.e-3 && AVmodel_data_elem[0]<1.e-6){ 
                     for(j=1; j<sr+1; j++)
                        AVmodel_data_elem[j] = AVmodel_data_neigh[j];
                     break;
                  }
             }
            }
*/
     //      xf_printf("bound %lf, In-cell %lf", EntropyBound[ginx], EntropyInCell[ginx]);
     //      if(!InitFlag)
     //      getchar();
     //
            //relaxtion
            //Way Number One
            //if(MinEntropyVec[0] > Model->LimM && MinEntropyVec[0] < 1.01 * Model->LimM)
            //MinEntropyVec[0] -= 1.0e-2;
            //if(MinEntropyVec[0] < 1.01 * Model->LimM)
            //   MinEntropyVec[0] = Model->LimM;

            //a further relaxation to avoid machine error
            MinEntropyVec[0] *= 0.99;
            if(MinEntropyVec[0] < Model->LimM)
               MinEntropyVec[0] = Model->LimM;
            
            MinEntropyVec[0] = Model->LimM;

            ginx++;
        }//elem
     }//egrp
    
  
   //let's get ride of this mess
/*   if(Model->AVmodel || p_Adapt_active){
      ierr = xf_Error(Refresh(&QuadDataElem, &QuadDataFace, &JData, &PhiData, &Quadx, &QuadStat, &QuadGrad));
      if (ierr != xf_OK) return ierr;
   }
   else{
      ierr = xf_Error(Refresh(&QuadDataElem, &QuadDataFace, &JData, &PhiData, &Quadx, &QuadStat, NULL));
      if (ierr != xf_OK) return ierr;
   }
         
   if (QuadPhiData != NULL){
      ierr = xf_Error(xf_DestroyBasisData(QuadPhiData, xfe_True));
      if (ierr != xf_OK) return ierr;
      QuadPhiData = NULL;
   }
*/

   Yu_DestroyStaticDataElem(StaticData);
   xf_Release( (void *) Quadx);
   xf_Release( (void *) QuadStat);
   xf_Release( (void *) EntropyInCell);
   xf_Release( (void *) EntropyBound);
    
    EntropyInCell= NULL;
    EntropyBound = NULL;

   return xf_OK;
}
