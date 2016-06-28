/*
  FILE:  Yu_DiffusDiscret.c

  This file contains the residual-calculation functions specific to
  diffusion terms.

*/

#include "xf_AllStruct.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_SolverStruct.h"
#include "xf_Data.h"
#include "xf_Basis.h"
#include "xf_Param.h"
#include "xf_Memory.h"
#include "xf_Math.h"
#include "xf_MeshTools.h"
#include "xf_EqnSetHook.h"
#include "xf_ResidualStab.h"
#include "xf_MeshMotionStruct.h"
#include "xf_LinQueue.h"
#include "xfYu_Model.h"
#include "xfYu_ModelStruct.h"
//#include "reaction_routine.h"

//cannot handle variables more than 50
#define MAXSR 50

/******************************************************************/
//   FUNCTION Definition: xf_JumpPenaltyEta
static int
xf_JumpPenaltyEta(//enum xfe_DiffDiscType DiffDisc,
                  enum xfe_BasisType Basis, int Order, int nFace,
                  real *eta)
{
    /* Coefficient for stabilization */
    
    int ierr;
    enum xfe_ShapeType Shape;
    int MaxOrderTriIP = 5;
    int MaxOrderTetIP = 5;
    int MaxOrderQuadIP = 8;
    int MaxOrderHexIP = 8;
    // real etaTriIP[6] = {1.0, 4.0, 11.0, 14.0, 17.0, 20.0};
    real etaTriIP[6] = {1.0, 4.0, 11.0, 20.0, 30.0, 40.0};
    //real etaQuadIP[9] = {1.0, 1.5, 3.0, 4.0, 8.0, 10.0, 20.0, 24.0, 50.0};
    //real etaQuadIP[9] = {1.0, 4.0, 8.0, 16.0, 20.0, 30.0, 35.0, 45.0, 50.0};
    real etaQuadIP[9] = {1.0, 4.0, 12.0, 12.0, 20.0, 30.0, 35.0, 45.0, 50.0};
    real etaTetIP[6] = {1.0, 5.0, 20.0, 30.0, 40.0, 50.0};
    //real etaHexIP[9] = {1.0, 2.5, 5.0, 8.0, 10.0, 14.0, 20.0, 24.0, 50.0};
    
    real etaHexIP[9] = {1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0};
    
    ierr = xf_Error(xf_Basis2Shape(Basis, &Shape));
    if (ierr != xf_OK) return ierr;
    
    //switch (DiffDisc){
    //    case xfe_DiffDiscIP:
            switch (Shape){
                case xfe_Triangle:
                    (*eta) = etaTriIP[min(Order,MaxOrderTriIP)]; break;
                case xfe_Quadrilateral:
                    (*eta) = etaQuadIP[min(Order,MaxOrderQuadIP)]; break;
                case xfe_Tetrahedron:
                    (*eta) = etaTetIP[min(Order,MaxOrderTetIP)]; break;
                case xfe_Hexahedron:
                    (*eta) = etaHexIP[min(Order,MaxOrderHexIP)]; break;
                default: 
                    return xf_Error(xf_NOT_SUPPORTED); break;
            }
   //         break;
   //     default:
   //         return xf_Error(xf_NOT_SUPPORTED);
   //         break;
   // }
    
    (*eta) *= nFace;
    
    return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: LYDG_CreateDiffJumpData
int
LYDG_CreateDiffJumpData(xf_DiffJumpData **pDData)
{
  int ierr;
  xf_DiffJumpData *DData;
  
  ierr = xf_Error(xf_Alloc( (void **) pDData, sizeof(xf_DiffJumpData), 1));
  if (ierr != xf_OK) return ierr;

  DData = (*pDData);

  DData->Need_Grad = xfe_True;
  DData->ConstAL   = xfe_False;
  DData->ConstAR   = xfe_False;
  DData->alpha     = 0.0;
  DData->duL_uL = DData->duL_uR = 0.0;
  DData->duR_uL = DData->duR_uR = 0.0;

  DData->dunL     = NULL; DData->dunR     = NULL;
  DData->ALguL    = NULL; DData->ARguR    = NULL;
  DData->ALdunL   = NULL; DData->ARdunR   = NULL;
  DData->AL       = NULL; DData->AR       = NULL;
  DData->A_uLguL  = NULL; DData->A_uRguR  = NULL;
  DData->A_uLdunL = NULL; DData->A_uRdunR = NULL;
  DData->N        = NULL; 
  DData->Qn       = NULL; 
  DData->Qn_dp    = NULL; 
  DData->Qn_u     = NULL;
  DData->Qn_gu    = NULL;
  DData->AwStabL  = NULL; DData->AwStabR  = NULL;
  DData->T        = NULL;
  DData->Aw       = NULL;
  DData->DnL      = NULL; DData->DnR      = NULL;
  DData->SLL      = NULL; DData->SLR      = NULL;
  DData->SRL      = NULL; DData->SRR      = NULL;
  DData->PL       = NULL; DData->PR       = NULL;

  DData->RL_Fholder = NULL;
  DData->RR_Fholder = NULL;

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: LYDG_ReAllocDiffJumpData
int
LYDG_ReAllocDiffJumpData(int sr, int nq, int dim, int nn, 
		       enum xfe_Bool Need_Grad, xf_DiffJumpData *DData)
{
  int ierr;
  int sr2, Tsize;

  sr2 = sr*sr;
  DData->Need_Grad = Need_Grad;

  ierr = xf_Error(xf_ReAlloc( (void **) &DData->dunL, dim*nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->dunR, dim*nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ALguL, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ARguR, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ALdunL, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ARdunR, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn, sr*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_dp, sr*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->N , nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->DnL, 2*dim*nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  DData->DnR = DData->DnL + dim*nn*sr;

  // BR2 needs these even without gradient
  Tsize = nn*nq;
  Tsize = max(Tsize, nn*sr);
  Tsize = max(Tsize, nn*nn);
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->T, Tsize, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->SLL, 4*dim*nn*nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  DData->SLR = DData->SLL +   dim*nn*nn;
  DData->SRL = DData->SLL + 2*dim*nn*nn;
  DData->SRR = DData->SLL + 3*dim*nn*nn;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->PL, 2*dim*nn*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  DData->PR = DData->PL + dim*nn*nq;

  //added July 2015
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->RL_Fholder, nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->RR_Fholder, nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  if (Need_Grad){
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AL, sr2*nq*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AR, sr2*nq*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uLguL , sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uRguR , sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uLdunL, sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uRdunR, sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_u, nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_gu, dim*nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AwStabL, dim*nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AwStabR, dim*nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Aw, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: LYDG_DestroyDiffJumpData
void
LYDG_DestroyDiffJumpData(xf_DiffJumpData *DData)
{
  xf_Release( (void *) DData->dunL);     xf_Release( (void *) DData->dunR);
  xf_Release( (void *) DData->AL);       xf_Release( (void *) DData->AR);
  xf_Release( (void *) DData->ALguL);    xf_Release( (void *) DData->ARguR);
  xf_Release( (void *) DData->ALdunL);   xf_Release( (void *) DData->ARdunR);
  xf_Release( (void *) DData->A_uLguL);  xf_Release( (void *) DData->A_uRguR);
  xf_Release( (void *) DData->A_uLdunL); xf_Release( (void *) DData->A_uRdunR);
  xf_Release( (void *) DData->N);
  xf_Release( (void *) DData->Qn);
  xf_Release( (void *) DData->Qn_dp);
  xf_Release( (void *) DData->Qn_u);
  xf_Release( (void *) DData->Qn_gu);
  xf_Release( (void *) DData->AwStabL);  xf_Release( (void *) DData->AwStabR);
  xf_Release( (void *) DData->T);
  xf_Release( (void *) DData->Aw);
  xf_Release( (void *) DData->DnL);
  xf_Release( (void *) DData->SLL);
  xf_Release( (void *) DData->PL);
  xf_Release( (void *) DData->RL_Fholder);
  xf_Release( (void *) DData->RR_Fholder);

  xf_Release( (void *) DData);
}
/******************************************************************/
//Calculate transport property according to "Model" specification
static int 
Yu_TransportProperty(Yu_Model *Model, const real p, const real rho, const real meanW,
                     const real *X, const real gamma, real *T_out, real *mu,
                     real *kappa, real *Diff)
{
   int i, sr, dim, Nspe;
   real T, R;

   sr = Model->nVars;
   dim = Model->dim;
   Nspe = sr - 2 - dim;
   
   if(Model->DetailChem)
   {
      //CK_Massfrac_to_Molefrac(Y, X);
      //CK_Moleweight_Mixture(Y, &meanW);
      //gas constant of the mixture
      R = Model->Ru / meanW;
      T = p/rho/R;
      
      //prvide transport property(vicosity, diffusion, conductivity)
      //CK_Viscous_Mixture(T, X, mu);
      //CK_Conduct_Mixture(T, X, kappa);
      //CK_Diffus_Species_Mixture(p, T, X, Diff);
      
   }
   else
   {
      R = Model->Ru / meanW;  //this is specific universal gas const
      
      if( Model->Sutherland == xfe_False)
      {
         //use default value for Model definition
         (*mu)    = Model->mu_c;
         (*kappa) = Model->kappa_c;
         for(i=0; i<=Nspe; i++)
            Diff[i] = Model->Diff_c;
         
         T = p/rho/R;

      }
      else
      {
         
         //back to dimensional form
         T = p / rho /R;
         (*mu) = (Model->mu_c) * pow((1.4-1.0) * T, 0.666);
         (*kappa) = (Model->kappa_c) * pow((1.4-1.0) * T, 0.666);
         //(*mu) = (Model->mu_c) * pow(T/214.2857144, 0.666);
         //(*kappa) = (Model->kappa_c) * pow(T/214.2857144, 0.666);
         for(i=0; i<=Nspe; i++)
            Diff[i] = Model->Diff_c;
       /*  T = (p*101325.0)/(rho*1.29)/R;
         
         //if Sutherland's Law is used for air only
         (*mu)  = 1.7894e-5 * pow(T/288.15, 1.5) * ((288.15 + 110.0)/(T + 110.0));
         (*mu)  = (*mu) / 1.29 / 1.0e-5 / 280.261506;
         (*kappa) = (*mu) * (gamma * R / (gamma - 1.0)) / 0.75;
         for(i=0; i<=Nspe; i++)
            //    Diff[i] = mu / rho;  //basically assume Lewis number equal to one
            Diff[i] = 0.0;
         
         //back to non-dimensional form
         
         T = T / 280.261506/ 280.261506;
         */
      }
   }
   
   (*T_out) = T;

   
   return xf_OK;
}

static int
MaxViscosity(Yu_Model *Model, const int nq, real *uL, real *uR, real *max_vis)
{
    int i, j, k, dim, sr, sr2;
    int Nspe, ierr;
    real rho, rhou, rhov, rhow, rhoE, rhoY[50], Y[50];
    real u, v, w, p, tmp, meanW, X[50];
    real T, mu, gamma, kappa, Diff[50];
    real *U_i;
    
    dim = Model->dim;
    sr  = Model->nVars;
    gamma = Model->GammaInit; 
     
    sr2  = sr*sr;
    Nspe = sr - 2 - dim;

    (*max_vis) = 0.;
    for(j=0; j<2; j++)
    for(k=0; k<nq; k++)
    {
        if(j==0)
            U_i = uL + k*sr;
        else
            U_i = uR + k*sr;
        //pull in variables
        //may have problem; double check
        rho  = U_i[0];
        u    = U_i[1]/rho;
        v    = U_i[2]/rho;
        if(dim == 3)
        w    = U_i[3]/rho;
        else
           w = 0.;
        
        rhoE = U_i[dim+1];
        for(i=2+dim; i<sr; i++)
            rhoY[i-2-dim] = U_i[i];
        
        p = (gamma-1.)*(rhoE - 0.5*rho*(u*u + v*v + w*w));
        tmp   = 0.;
        meanW = 0.;
        for(i=0; i<Nspe; i++)
        {
            Y[i]   = rhoY[i]/rho;
            tmp   += Y[i];
            meanW += Y[i]/(Model->moleW[i]);
        }
        Y[Nspe] = 1.0 - tmp;
        meanW  += Y[Nspe]/(Model->moleW[Nspe]);
        meanW   = 1./meanW;
        
        //mole fraction
        for(i=0; i<=Nspe; i++)
            X[i] = meanW*Y[i]/(Model->moleW[i]);
        
        //pull in thermal transport parameters; modulized
        ierr = xf_Error(Yu_TransportProperty(Model, p, rho, meanW, X, gamma,
                                             &T, &mu, &kappa, Diff));
    
    
        if(mu > *max_vis) *max_vis = mu;
    }

    return xf_OK;
}
/******************************************************************/
//Diffusion discretization two dimensional version
static int
LYDG_PhysDiffModel2d(const int sr, const int dim, const int off,
                   const int off2, const real *U_i, const real *gU_i, 
                   const real *W_i, real *Aw_i, real *A_i, const real gamma, 
                   Yu_Model *Model)
{
    int ierr, sr2, i, j, k, Nspe;
    real mu, kappa, moleWn, moleWi; 
    real rho, rhou, rhov, rhoY[MAXSR], rhoE, u, v, p;
    real Y[MAXSR], X[MAXSR], Diff[MAXSR], meanW, R;
    real T, T_U[MAXSR+50], tmp;
    real A0[MAXSR*MAXSR];
    real *A00, *A01, *A10, *A11;
    real V2, c1, c2;
    const real *W0, *W1;
    
    sr2  = sr*sr;
    Nspe = sr - 2 - dim;
    
    //pull in variables
    //may have problem; double check
    rho  = U_i[0];
    rhou = U_i[1];
    rhov = U_i[2];
    rhoE = U_i[3];
    for(i=2+dim; i<sr; i++)
        rhoY[i-2-dim] = U_i[i];
    
    //back out primitive variables
    u = rhou/rho;
    v = rhov/rho;
    p = (gamma-1.)*(rhoE - 0.5*rho*(u*u + v*v));
    tmp   = 0.;
    meanW = 0.;
    for(i=0; i<Nspe; i++)
    {
        Y[i]   = rhoY[i]/rho;
        tmp   += Y[i];
        meanW += Y[i]/(Model->moleW[i]);
    }
    Y[Nspe] = 1.0 - tmp;
    meanW  += Y[Nspe]/(Model->moleW[Nspe]);
    meanW   = 1./meanW;
   
    //mole fraction
    for(i=0; i<=Nspe; i++)
       X[i] = meanW*Y[i]/(Model->moleW[i]);

    //pull in thermal transport parameters; modulized
    ierr = xf_Error(Yu_TransportProperty(Model, p, rho, meanW, X, gamma,
                                         &T, &mu, &kappa, Diff));
    if (ierr != xf_OK) return ierr;
   
    //set pointer from A
    if (A_i != NULL){
        A00 = A_i       ;
        A01 = A_i+  off2;
        A10 = A_i+2*off2;
        A11 = A_i+3*off2;
    }
    else
    {
        A00 = A01 = A10 = A11 = A0;
    }
    
    if(W_i != NULL){
        W0 = W_i;
        W1 = W_i + off;
    }
    
    //provide linearization for T
    R = (Model->Ru)/meanW;
    V2 = u*u + v*v;
    moleWn =  Model->moleW[Nspe];
    //rho
    T_U[0] = 0.5*(gamma-1.)/rho/R*V2 - T/rho/R*Model->Ru/moleWn;
    //rhou
    T_U[1] = (gamma-1.)/rho/R*(-u);
    //rhov
    T_U[2] = (gamma-1.)/rho/R*(-v);
    //rhoE
    T_U[3] = (gamma-1.)/rho/R;
    //rhoY
    for(i=0; i<Nspe; i++) {
        moleWi = Model->moleW[i];
        T_U[4+i] = -T/rho/R*Model->Ru*(1./moleWi - 1./moleWn);
    }
   
    //constants for use below
    c1 = 4./3.;
    c2 = 2./3.;
    
    /* A00 */
    for(k=0; k<sr2; k++) A0[k] = 0.;
    A0[sr*1 + 0]   =    -c1*mu/rho*u;
    A0[sr*1 + 1]   =     c1*mu/rho;
    A0[sr*2 + 0]   =    -mu/rho*v;
    A0[sr*2 + 2]   =     mu/rho;
    A0[sr*3 + 0]   =    -mu/rho*(c1*u*u + v*v) + kappa*T_U[0];
    A0[sr*3 + 1]   =     c1*mu/rho*u + kappa*T_U[1];
    A0[sr*3 + 2]   =     mu/rho*v + kappa*T_U[2];
    A0[sr*3 + 3]   =     kappa*T_U[3];
    //for conduction in energy
    for(i=4; i<sr; i++)
        A0[sr*3 + i] = kappa*T_U[i];
    //use the simplest model for species diffusion
    for(i=4; i<sr; i++) {
        A0[sr*i + 0]  =  -Diff[i-4]*Y[i-4];
        A0[sr*i + i]  = Diff[i-4];
    }
    if(A_i != NULL)
        xf_V_Add(A0, sr2, xfe_Add, A00);
    if(Aw_i != NULL)
        xf_MxV(A0, W0, sr, sr, xfe_Add, Aw_i+0*off);
    
    /* A01 */
    for(k=0; k<sr2; k++) A0[k] = 0.;
    A0[sr*1 + 0]  =    c2*mu/rho*v;
    A0[sr*1 + 2]  =   -c2*mu/rho;
    A0[sr*2 + 0]  =   -mu/rho*u;
    A0[sr*2 + 1]  =    mu/rho;
    A0[sr*3 + 0]  =   (c2-1.)*mu/rho*u*v;
    A0[sr*3 + 1]  =    mu/rho*v;
    A0[sr*3 + 2]  =   -c2*mu/rho*u;
    //no cross-differentiation term for conduction and diffusion
    if(A_i != NULL)
        xf_V_Add(A0, sr2, xfe_Add, A01);
    if(Aw_i != NULL)
        xf_MxV_Add(A0, W1, sr, sr, Aw_i+0*off);
    
    /* A01 */
    for(k=0; k<sr2; k++) A0[k] = 0.;
    A0[sr*1 + 0]   =  -mu/rho*v;
    A0[sr*1 + 2]   =   mu/rho;
    A0[sr*2 + 0]   =   c2*mu/rho*u;
    A0[sr*2 + 1]   =  -c2*mu/rho;
    A0[sr*3 + 0]   =  (c2-1.)*mu/rho*u*v;
    A0[sr*3 + 1]   =  -c2*mu/rho*v;
    A0[sr*3 + 2]   =   mu/rho*u;
    //no cross-differentiation term for conduction and diffusion
    if(A_i != NULL)
        xf_V_Add(A0, sr2, xfe_Add, A10);
    if(Aw_i != NULL)
        xf_MxV(A0, W0, sr, sr, xfe_Add, Aw_i+1*off);
    
    /* A11 */
    for(k=0; k<sr2; k++) A0[k] = 0.;
    A0[sr*1 + 0]  =  -mu/rho*u;
    A0[sr*1 + 1]  =   mu/rho;
    A0[sr*2 + 0]  =  -c1*mu/rho*v;
    A0[sr*2 + 2]  =   c1*mu/rho;
    A0[sr*3 + 0]  = -mu/rho*(u*u + c1*v*v) + kappa*T_U[0];
    A0[sr*3 + 1]  =  mu/rho*u + kappa*T_U[1];
    A0[sr*3 + 2]  =  c1*mu/rho*v + kappa*T_U[2];
    A0[sr*3 + 3]  =  kappa*T_U[3];
    //for conduction in energy
    for(i=4; i<sr; i++)
        A0[sr*3 + i] = kappa*T_U[i];
    //use the simple model for species diffusion
    for(i=4; i<sr; i++) {
        A0[sr*i + 0]  =  -Diff[i-4]*Y[i-4];
        A0[sr*i + i]  =   Diff[i-4];
    }
    if(A_i != NULL)
        xf_V_Add(A0, sr2, xfe_Add, A11);
    if(Aw_i != NULL)
        xf_MxV_Add(A0, W1, sr, sr, Aw_i+1*off);
    
    return xf_OK;
}

/******************************************************************/
//Diffusion discretization two dimensional version
static int
LYDG_PhysDiffModel3d(const int sr, const int dim, const int off,
                     const int off2, const real *U_i, const real *gU_i,
                     const real *W_i, real *Aw_i, real *A_i, const real gamma,
                     Yu_Model *Model)
{
   int ierr, sr2, i, j, k, Nspe, speoff;
   real mu, kappa, nu, moleWn, moleWi;
   real rho, rhou, rhov, rhow, rhoY[50], rhoE, u, v, w, p;
   real Y[50], X[50], Diff[50], meanW, R;
   real T, T_U[MAXSR+50], tmp;
   real A0[MAXSR*MAXSR];
   real *A00, *A01, *A02;
   real *A10, *A11, *A12;
   real *A20, *A21, *A22;
   real V2, c1, c2, eps;
   const real *W0, *W1, *W2;
   
   sr2  = sr*sr;
   Nspe = sr - 2 - dim;
   eps  = 1.0-14;  //avoid zero denominator
   
   //pull in variables
   //may have problem; double check
   rho  = U_i[0];
   rhou = U_i[1];
   rhov = U_i[2];
   rhow = U_i[3];
   rhoE = U_i[4];
   for(i=2+dim; i<sr; i++)
      rhoY[i-2-dim] = U_i[i];
   
   //back out primitive variables
   u = rhou/rho;
   v = rhov/rho;
   w = rhow/rho;
   p = (gamma-1.)*(rhoE - 0.5*rho*(u*u + v*v + w*w));
   tmp   = 0.;
   meanW = 0.;
   for(i=0; i<Nspe; i++)
   {
      Y[i]   = rhoY[i]/rho;
      tmp   += Y[i];
      meanW += Y[i]/(Model->moleW[i]);
   }
   Y[Nspe] = 1.0 - tmp;
   meanW  += Y[Nspe]/(Model->moleW[Nspe]);
   meanW   = 1./meanW;
   
   //mole fraction
   for(i=0; i<=Nspe; i++)
      X[i] = meanW*Y[i]/(Model->moleW[i]);
   
   //pull in thermal transport parameters
   ierr = xf_Error(Yu_TransportProperty(Model, p, rho, meanW, X, gamma,
                                        &T, &mu, &kappa, Diff));
   if (ierr != xf_OK) return ierr;

   //set pointer from A
   if (A_i != NULL){
      A00 = A_i       ;  A01 = A_i+  off2; A02 = A_i+2*off2;
      A10 = A_i+3*off2;  A11 = A_i+4*off2; A12 = A_i+5*off2;
      A20 = A_i+6*off2;  A21 = A_i+7*off2; A22 = A_i+8*off2;
   }
   else{
      A00 = A01 = A02 = A10 = A11 = A12 = A20 = A21 = A22 = A0;
   }
   
   if (W_i != NULL){
      W0 = W_i;
      W1 = W_i +  off;
      W2 = W_i + 2*off;
   }
   
   //provide linearization for T
   R = (Model->Ru)/meanW;
   V2 = u*u + v*v + w*w;
   moleWn =  Model->moleW[Nspe];
   //rho
   T_U[0] = 0.5*(gamma-1.)/rho/R*V2 - T/rho/R*Model->Ru/moleWn;
   //rhou
   T_U[1] = (gamma-1.)/rho/R*(-u);
   //rhov
   T_U[2] = (gamma-1.)/rho/R*(-v);
   //rhow
   T_U[3] = (gamma-1.)/rho/R*(-w);
   //rhoE
   T_U[4] = (gamma-1.)/rho/R;
   //rhoY
   for(i=0; i<Nspe; i++) {
      moleWi = Model->moleW[i];
      T_U[2+dim+i] = -T/rho/R*Model->Ru*(1./moleWi - 1./moleWn);
   }
   
   //constants for use below
   c1 = 4./3.;
   c2 = 2./3.;
   speoff = 2+dim;
   nu = mu/rho;
   
   /* A00 */
   for(k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*1 + 0]   =    -c1*nu*u;
   A0[sr*1 + 1]   =     c1*nu;
   A0[sr*2 + 0]   =    -nu*v;
   A0[sr*2 + 2]   =     nu;
   A0[sr*3 + 0]   =    -nu*w;
   A0[sr*3 + 3]   =     nu;
   A0[sr*4 + 0]   =    -nu*(c1*u*u + v*v + w*w) + kappa*T_U[0];
   A0[sr*4 + 1]   =     c1*nu*u + kappa*T_U[1];
   A0[sr*4 + 2]   =     nu*v + kappa*T_U[2];
   A0[sr*4 + 3]   =     nu*w + kappa*T_U[3];
   A0[sr*4 + 4]   =     kappa*T_U[4];
   //for conduction in energy
   for(i=speoff; i<sr; i++)
      A0[sr*4 + i] = kappa*T_U[i];
   //use the simplest model for species diffusion
   for(i=speoff; i<sr; i++) {
      A0[sr*i + 0]  =  -Diff[i-speoff]*Y[i-speoff];
      A0[sr*i + i]  = Diff[i-speoff];
   }
   if(A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A00);
   if(Aw_i != NULL)
      xf_MxV(A0, W0, sr, sr, xfe_Add, Aw_i+0*off);
   
   /* A01 */
   for(k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*1 + 0]  =    c2*nu*v;
   A0[sr*1 + 2]  =   -c2*nu;
   A0[sr*2 + 0]  =   -nu*u;
   A0[sr*2 + 1]  =    nu;
   A0[sr*4 + 0]  =   (c2-1.)*nu*u*v;
   A0[sr*4 + 1]  =    nu*v;
   A0[sr*4 + 2]  =   -c2*nu*u;
   //no cross-differentiation term for conduction and diffusion
   if(A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A01);
   if(Aw_i != NULL)
      xf_MxV_Add(A0, W1, sr, sr, Aw_i+0*off);
   
   /* A02 */
   for (k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*1 + 0] =  c2*nu*w;
   A0[sr*1 + 3] = -c2*nu;
   A0[sr*3 + 0] = -nu*u;
   A0[sr*3 + 1] =  nu;
   A0[sr*4 + 0] = (c2-1.)*nu*u*w;
   A0[sr*4 + 1] =  nu*w;
   A0[sr*4 + 3] = -c2*nu*u;
   if (A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A02);
   if (Aw_i != NULL)
      xf_MxV_Add(A0, W2, sr, sr, Aw_i+0*off);
   
   /* A10*/
   for(k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*1 + 0]   =  -nu*v;
   A0[sr*1 + 2]   =   nu;
   A0[sr*2 + 0]   =   c2*nu*u;
   A0[sr*2 + 1]   =  -c2*nu;
   A0[sr*4 + 0]   =  (c2-1.)*nu*u*v;
   A0[sr*4 + 1]   =  -c2*nu*v;
   A0[sr*4 + 2]   =   nu*u;
   //no cross-differentiation term for conduction and diffusion
   if(A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A10);
   if(Aw_i != NULL)
      xf_MxV(A0, W0, sr, sr, xfe_Add, Aw_i+1*off);
   
   /* A11 */
   for(k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*1 + 0]  =  -nu*u;
   A0[sr*1 + 1]  =   nu;
   A0[sr*2 + 0]  =  -c1*nu*v;
   A0[sr*2 + 2]  =   c1*nu;
   A0[sr*3 + 0]  =  -nu*w;
   A0[sr*3 + 3]  =   nu;
   A0[sr*4 + 0]  = -nu*(u*u + c1*v*v + w*w) + kappa*T_U[0];
   A0[sr*4 + 1]  =  nu*u + kappa*T_U[1];
   A0[sr*4 + 2]  =  c1*nu*v + kappa*T_U[2];
   A0[sr*4 + 3]  =  nu*w + kappa*T_U[3];
   A0[sr*4 + 4]  =  kappa*T_U[4];
   //for conduction in energy
   for(i=speoff; i<sr; i++)
      A0[sr*4 + i] = kappa*T_U[i];
   //use the simple model for species diffusion
   for(i=speoff; i<sr; i++) {
      A0[sr*i + 0]  =  -Diff[i-speoff]*Y[i-speoff];
      A0[sr*i + i]  =   Diff[i-speoff];
   }
   if(A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A11);
   if(Aw_i != NULL)
      xf_MxV_Add(A0, W1, sr, sr, Aw_i+1*off);
   
   /* A12 */
   for (k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*2 + 0] =  c2*nu*w;
   A0[sr*2 + 3] = -c2*nu;
   A0[sr*3 + 0] = -nu*v;
   A0[sr*3 + 2] =  nu;
   A0[sr*4 + 0] = (c2-1.)*nu*v*w;
   A0[sr*4 + 2] =  nu*w;
   A0[sr*4 + 3] = -c2*nu*v;
   if (A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A12);
   if (Aw_i != NULL)
      xf_MxV_Add(A0, W2, sr, sr, Aw_i+1*off);
   
   /* A20 */
   for (k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*1 + 0] = -nu*w;
   A0[sr*1 + 3] =  nu;
   A0[sr*3 + 0] =  c2*nu*u;
   A0[sr*3 + 1] = -c2*nu;
   A0[sr*4 + 0] = (c2-1.)*nu*u*w;
   A0[sr*4 + 1] = -c2*nu*w;
   A0[sr*4 + 3] =  nu*u;
   if (A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A20);
   if (Aw_i != NULL)
      xf_MxV(A0, W0, sr, sr, xfe_Add, Aw_i+2*off);
   
   /* A21 */
   for (k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*2 + 0] = -nu*w;
   A0[sr*2 + 3] =  nu;
   A0[sr*3 + 0] =  c2*nu*v;
   A0[sr*3 + 2] = -c2*nu;
   A0[sr*4 + 0] = (c2-1.)*nu*v*w;
   A0[sr*4 + 2] = -c2*nu*w;
   A0[sr*4 + 3] =  nu*v;
   if (A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A21);
   if (Aw_i != NULL)
      xf_MxV_Add(A0, W1, sr, sr, Aw_i+2*off);
   
   /* A22 */
   for (k=0; k<sr2; k++) A0[k] = 0.;
   A0[sr*1 + 0] = -nu*u;
   A0[sr*1 + 1] =  nu;
   A0[sr*2 + 0] = -nu*v;
   A0[sr*2 + 2] =  nu;
   A0[sr*3 + 0] = -c1*nu*w;
   A0[sr*3 + 3] =  c1*nu;
   A0[sr*4 + 0] = -nu*(u*u + v*v + c1*w*w) + kappa*T_U[0];
   A0[sr*4 + 1] =  nu*u + kappa*T_U[1];
   A0[sr*4 + 2] =  nu*v + kappa*T_U[2];
   A0[sr*4 + 3] =  c1*nu*w + kappa*T_U[3];
   A0[sr*4 + 4] =  kappa*T_U[4];
   //for conduction in energy
   for(i=speoff; i<sr; i++)
      A0[sr*4 + i] = kappa*T_U[i];
   //use the simple model for species diffusion
   for(i=speoff; i<sr; i++) {
      A0[sr*i + 0]  =  -Diff[i-speoff]*Y[i-speoff];
      A0[sr*i + i]  =   Diff[i-speoff];
   }
   if (A_i != NULL)
      xf_V_Add(A0, sr2, xfe_Add, A22);
   if (Aw_i != NULL)
      xf_MxV_Add(A0, W2, sr, sr, Aw_i+2*off);
   
   return xf_OK;
}

static int
LYDG_LaplaceDiffModel(const int sr, const int dim, const int off, const int off2, 
                      const real *U_i, const real *gU_i,
                      const real *W_i, real *Aw_i, real *A_i, const real gamma,
                      Yu_Model *Model, const real *AV)
{
   int sr2, i, j, k, l;
   real u, v, w, E, rij, r, p, c;
   real R, Pr, gmi, V2, Y[50];
   real A0[MAXSR*MAXSR];
   enum xfe_AddType AddFlag, AddFlag2, CurAddFlag;
   const real *Wj;

   sr2 = sr*sr;
/*
   r = U_i[0];
   u = U_i[1]/U_i[0];
   v = U_i[2]/U_i[0];
   w = (dim == 2) ? 0.0 : U_i[dim]/U_i[0];
   E = U_i[dim+1]/U_i[0];
   for(i=dim+2; i<sr; i++)
      Y[i-dim-2] = U_i[i]/U_i[0];
  
   V2 = u*u + v*v + w*w;
   gmi = gamma - 1.0;
*/
   AddFlag = xfe_Add;
   AddFlag2 = xf_GetAddFlag2(AddFlag);

   //Aij
   for (i=0; i<dim; i++){
      for (j=0; j<dim; j++){
         
         //Evaluate ResMetric for AV 
         //rij = ((ResMetric == NULL) ? 1.0 : ResMetric[i*dim+j]);
         //here is for linear case test
         //now incorporate the AV model
         if(i == j){ 
         //rij = AV[0]; //LapVis;
         rij = 1.0;
         }
         else
            rij = 0.0;

         //zero out A0
         for (k=0; k<sr2; k++) A0[k] = 0.;

         for (k=0; k<sr; k++) A0[sr*k + k] += rij*AV[1+k]; // Laplace portion
         
            CurAddFlag = (j==0) ? AddFlag : AddFlag2;

         if(A_i  != NULL)
            xf_V_Add(A0, sr2, AddFlag, A_i + off2*(i*dim+j));
         if(Aw_i != NULL)
            xf_MxV(A0, W_i+j*off, sr, sr, CurAddFlag, Aw_i+i*off);

      }//i
   }//j

   return xf_OK;


}
/******************************************************************/
//   FUNCTION Definition: LYDG_InteriorViscTerm
int
LYDG_InteriorViscTerm(Yu_Model *Model, int nDiff, int nq, real *u,
                      real *gu, real *w, real *Aw, real *A, const real gamma, 
                      const real *AV)
{
  int ierr, i, j, d, dim, sr, sr2, iDiff, k, iq;
  int off, off2;
  const real *U_i, *W_i;
  real *gU_i, gU0[3*MAXSR];
  real *Aw_i, *A_i, mu;
  
  dim = Model->dim;
  sr  = Model->nVars;
  sr2 = sr*sr;

  off = sr *nq;
  off2= sr2*nq;

  if (nDiff <= 0) return xf_Error(xf_CODE_LOGIC_ERROR);
  
  // Initialize Aw, A, etc. to zero
  if (Aw     != NULL) for (k=0; k<dim*nq*sr     ; k++) Aw[k]     = 0.0;
  if (A      != NULL) for (k=0; k<dim*dim*nq*sr2; k++) A[k]      = 0.0;

  //modify here to add more specification 
//  ierr = xf_Error(xf_EqnSetDiffA(EqnSet, ResTerm+iDiff, IParam, RParam,
//				 ResMetric, nq, u, gu, w, Aw, A, A_uw, 
//				 &ConstATemp, CF));
//  if (ierr != xf_OK) return ierr;

  //loop over each quadrature point
  for (i=0; i<nq; i++){
     U_i = u + i*sr;

     //pull off state gradient if not NULL
     if(gu == NULL) gU_i = NULL;
     else {
        gU_i = gU0;
        for(j=0; j<dim; j++)
           for(k=0; k<sr; k++)
              gU_i[j*sr+k] = gu[j*nq*sr + i*sr + k];
     }

     W_i  = ((w   == NULL) ? NULL : w     + i*sr );
     Aw_i = ((Aw  == NULL) ? NULL : Aw    + i*sr );
     A_i  = ((A   == NULL) ? NULL : A     + i*sr2);

if(LinearCase)
{
   //make sense to use Laplacian diffusion model 
   ierr = xf_Error(LYDG_LaplaceDiffModel(sr, dim, off, off2, U_i, gU_i, W_i, Aw_i, A_i, gamma, Model, AV));
   if(ierr != xf_OK) return ierr;
}
else
{
     if(dim == 2)
     {
        if(Model->AVmodel)
        {
           ierr = xf_Error(LYDG_LaplaceDiffModel(sr, dim, off, off2, U_i, gU_i, W_i, Aw_i, A_i, gamma, Model, AV));
           if(ierr != xf_OK) return ierr;
        }
        else
        {
           ierr = xf_Error(LYDG_PhysDiffModel2d(sr, dim, off, off2, U_i, gU_i, W_i, Aw_i, A_i, gamma, Model));
           if(ierr != xf_OK) return ierr;
        }
     }
     else if(dim == 3)
     {
        ierr = xf_Error(LYDG_PhysDiffModel3d(sr, dim, off, off2, U_i, gU_i, W_i, Aw_i,
                                             A_i, gamma, Model));
        if(ierr != xf_OK) return ierr;
     }
     else
        return xf_NOT_SUPPORTED;
}
   }//i


   return xf_OK;
}



/******************************************************************/
/* we provide the linearization for diffusion term in Navier-Stokes
 * Q{k, i} = [0;
 *            tau{i,0};
 *            tau{i,1};
 *            tau{i,2};
 *            u{j}*tau{i,j} + kappa*g{i}T
 *            ]
 where, g{i} = gradient operator in i'th spatial direction
 tau{i,j} = mu*(g{i}u{j} + g{j}u{i} - 2/3*g{k}u{k} delta {i,j})
 delta{i,j} = 1 if i==j, 0 otherwise
 u{i} = velocity component in i'th spatial direction
 mu = visocity, which depends on the entire state vector

 Writing Q{k,i} = A{k,i,l,j}*g{j}U{l}, where U{l} is the conservative state means

 A{k,0,l,0}*g{0}U{l} = 
 [ 0;
   4/3*mu*g{0}u{0};
   mu*g{0}u{1};
   mu*g{0}u{2};
   mu*(4/3*g{0}u{0}*u{0} + g{0}u{1}*u{1} + g{0}u{2}*u{2}) + kappa*g{0}T
 ]

We can convert the gradients above to gradients of the conservative state variable
using:

g{i}u{j} = -(1/r)*u{j}*g{i}r + (1/r)*g{i}ru{j}
where,
r = density

Also, we write g{i}T = T_U{l}*g{i}U{l}
where,
T = (gam-1)/R * (E - .5*u{k}*u{k})
so that
T_U{l} = (gam-1)/(R*r) * [-E + u{k}*u{k}, -u{0}, -u{1}, -u{2}, 1]

Thus, we have
A{k,0,l,0} = 
[
      0         |     0       |    0     |    0      |    0     |
-----------------------------------------------------------------
-4/3*mu/r*u{0}  | 4/3*mu/r    |    0     |    0      |    0     |
-----------------------------------------------------------------
-mu/r*u{1}      |     0       |  mu/r    |    0      |    0     |
-----------------------------------------------------------------
-mu/r*u{2}      |     0       |    0     |  mu/r     |    0     |
-----------------------------------------------------------------
-4/3*mu/r*u{0}^2|4/3*mu/r*u{0}|mu/r*u{1}+|mu/r*u{2}  |kap*T_U{4}|
 -mu/r*u{1}^2   |+kap*T_U{1}  |kap*T_U{2}|+kap*T_U{3}|
 -mu/r*u{2}^2   |
 +kap*T_U{0}    |

Similarly,
A{k,0,l,1}*g{1}u{l} = 
[ 0;
  -2/3*mu*g{1}u{1};
  mu*g{1}u{0};
  0;
  mu*(-2/3*g{1}u{1}*u{0} + g{1}u{0}*u{1})
]
so that A{k,0,l,1} = 
{
       0        |      0      |     0      |      0     |     0
----------------|-------------|------------|------------|-----------
  2/3*mu/r*u{1} |      0      |  -2/3*mu/r |      0     |     0
----------------|-------------|------------|------------|-----------
    -mu/r*u{0}  |   mu/r      |     0      |      0     |     0
 ---------------|-------------|------------|------------|-----------
       0        |      0      |     0      |      0     |     0
 ---------------|-------------|------------|------------|-----------
2/3*mu/r*u{1}*  |  mu/r*u{1}  |-2/3*mu/r   |      0     |     0
u{0}-mu/r*u{0}*u{1}|          |*u{0}
]

A{k,0,l,2}*g{2}U{l} = 
[ 0;
  -2/3*mu*g{2}u{2};
  0;
  mu*g{2}u{0};
  mu*(-2/3*g{2}u{2}*u{0} + g{2}u{0}*u{2})
]
so that A{k,0,l,2} = 
[
       0        |      0      |     0      |      0     |     0
----------------|-------------|------------|------------|-----------
  2/3*mu/r*u{2} |      0      |  -2/3*mu/r |      0     |     0
----------------|-------------|------------|------------|-----------
       0        |      0      |     0      |      0     |     0
 ---------------|-------------|------------|------------|-----------
    -mu*u{0}    |    mu/r     |     0      |      0     |     0
 ---------------|-------------|------------|------------|-----------
2/3*mu/r*u{2}*  |  mu/r*u{2}  |     0      |-2/3*mu/r   |     0
u{0}-mu/r*u{0}*u{2}|                       |*u{0}
}

A{k,1,l,0}*g{0}U{l} = 
[  0;
  mu*g{0}u{1};
  -2/3*mu*g{0}u{0};
   0;
  mu*(g{0}u{1}*u{0} - 2/3*g{0}u{0}*u{1})
]
so that A{k,1,l,0} = 
[
       0        |      0      |     0      |      0     |     0
----------------|-------------|------------|------------|-----------
   -mu/r*u{1}   |      0      |   mu/r     |      0     |     0
----------------|-------------|------------|------------|-----------
  2/3*mu/r*u{0} | -2/3*mu/r   |     0      |      0     |     0
 ---------------|-------------|------------|------------|-----------
       0        |      0      |     0      |      0     |     0
 ---------------|-------------|------------|------------|-----------
-mu/r*u{1}*u{0} |-2/3*mu/r*u{1}|mu/r*u{0}  |      0     |     0
+2/3*mu/r*u{0}*u{1}|                     
]

A{k,1,l,1}*g{1}U{l} = 
[  0;
  mu*g{1}u{0};
  4/3*mu*g{1}u{1};
  mu*g{1}u{2};
  mu*(g{1}u{0}*u{0} + 4/3*g{1}u{1}*u{1} + g{1}u{2}*u{2}) + kappa*g{1}T
]
so that A{k,1,l,1} = 
[
        0       |      0      |     0     |      0     |     0
----------------|-------------|-----------|------------|-----------
   -mu/r*u{0}   |    mu/r     |     0     |      0     |     0
----------------|-------------|-----------|------------|-----------
 -4/3*mu/r*u{1} |      0      |  4/3*mu/r |      0     |     0
----------------|-------------|-----------|------------|-----------
   -mu/r*u{2}   |      0      |     0     |   mu/r     |     0
----------------|-------------|-----------|------------|-----------
 -mu/r*u{0}^2   | mu/r*u{0}   |4/3*mu/r*u{1}|mu/r*u{2} | kap*T_U{4}
-4/3*mu/r*u{1}^2| +kap*T_U{1} |+kap*T_U{2}| +kap*T_U{3}|
  -mu/r*u{2}^2  |
   +kap*T_U{0}  |
]

A{k,1,l,2}*g{2}U{l} = 
[ 0;
  0;
  -2/3*mu*g{2}u{2};
 mu*g{2}u{1};
 mu*(-2/3*g{2}u{2}*u{1} + g{2}u{1}*u{2})
]
so that A{k,1,l,2} = 
[
      0        |      0      |     0      |      0     |     0
---------------|-------------|------------|------------|-----------
      0        |      0      |     0      |      0     |     0
---------------|-------------|------------|------------|-----------
 2/3*mu/r*u{2} |      0      |     0      |   -2/3*mu/r|     0
---------------|-------------|------------|------------|-----------
   -mu/r*u{1}  |      0      |   mu/r     |      0     |     0
---------------|-------------|------------|------------|-----------
2/3*mu/r*u{2}*u{1}|   0      | mu/r*u{2}  |-2/3*mu/r*u{1}|     0
  -mu/r*u{1}*u{2}| 
]

A{k,2,l,0}*g{0}U{l} = 
[ 0;
 mu*g{0}u{2};
  0;
 -2/3*mu*g{0}u{0};
 mu*(g{0}u{2}*u{0} - 2/3*g{0}u{0}*u{2})
]

so that A{k,2,l,0} = 
[
      0        |      0      |     0      |      0     |     0
---------------|-------------|------------|------------|-----------
  -mu/r*u{2}   |      0      |     0      |    mu/r    |     0
---------------|-------------|------------|------------|-----------
      0        |      0      |     0      |      0     |     0
---------------|-------------|------------|------------|-----------
 2/3*mu/r*u{0} |  -2/3*mu/r  |     0      |      0     |     0
---------------|-------------|------------|------------|-----------
-mu/r*u{2}*u{0}| -2/3*mu/r*u{2}|     0      | mu/r*u{0}|     0
+2/3*mu/r*u{0}*u{2}|            
]

A{k,2,l,1}*g{1}U{l} = 
[ 0;
  0;
  mu*g{1}u{2};
  -2/3*mu*g{1}u{1};
  mu*(g{1}u{2}*u{1} - 2/3*g{1}u{1}*u{0})
]

so that A{k,2,l,1} = 
[
      0        |      0      |     0      |      0     |     0
---------------|-------------|------------|------------|-----------
      0        |      0      |     0      |      0     |     0
---------------|-------------|------------|------------|-----------
    -mu/r*u{2} |      0      |     0      |   mu/r     |     0
---------------|-------------|------------|------------|-----------
 2/3*mu/r*u{1} |      0      |  -2/3*mu/r |      0     |     0
---------------|-------------|------------|------------|-----------
  -mu*u{2}*u{1}|      0      |-2/3*mu/r*u{2}|mu/r*u{1} |     0
+2/3*mu/r*u{1}*u{2}|
]

A{k,2,l,2}*g{2}U{l} = 
[ 0;
  mu*g{2}u{0};
  mu*g{2}u{1};
  4/3*mu*g{2}u{1};
  mu*(g{2}u{0}*u{0} + g{2}u{1}*u{1} + 4/3*g{2}u{2}*u{2}) + kap*g{1}T
]

so that A{k,2,l,2} = 
[
      0       |      0      |     0     |      0     |     0
--------------|-------------|-----------|------------|-----------
  -mu/r*u{0}  |   mu/r      |     0     |      0     |     0
--------------|-------------|-----------|------------|-----------
  -mu/r*u{1}  |      0      |  mu/r     |      0     |     0
--------------|-------------|-----------|------------|-----------
-4/3*mu/r*u{2}|      0      |     0     | 4/3*mu/r   |     0
--------------|-------------|-----------|------------|-----------
-mu/r*u{0}^2  | mu/r*u{0}   | mu/r*u{1} |4/3*mu/r*u{2}| kap*T_U{4}
-mu/r*u{1}^2  | +kap*T_U{1} |+kap*T_U{2}| +kap*T_U{3}|
-4/3*mu/r*u{2}^2  |
  +kap*T_U{0} |
]

*/

/******************************************************************/
//   FUNCTION Definition: LYDG_FaceViscTerm
//   BR2 DG duffusion discretization is employed
//  July 2015; incorporate Yu's modification for reducing numerical dissipation
   
int
LYDG_FaceViscTerm(xf_All *All, Yu_Model *Model, const int iiface, const xf_BasisData *PhiDataL,
                const xf_BasisData *PhiDataR, const xf_BasisData *ResPhiDataL,
                const xf_BasisData *ResPhiDataR, real ElemVolL, real ElemVolR,
                enum xfe_BasisType BasisL, enum xfe_BasisType BasisR, int OrderL,
                int OrderR, int nq, const real *wn, real *uL, real *uR, real *guL,
                real *guR, real *RL, real *RR, const real gammaL, const real gammaR,
                xf_DiffJumpData *DData)
{
  int ierr, i, j, k, iq, dim, ii, sr, sr2;
  int egrpL, egrpR, elemL, elemR, nFace, Order;
  int nL, nR, n, n2;
  int ResOrderL, ResOrderR, RnL, RnR;
  real *uLghost, *uRghost, max_vis;
  real *T, nval, eta, FaceArea, ihL, ihR;
  real *iMML, *iMMR, facL, facR;
  real *wL = NULL, *wR = NULL;
  real *AVelemL, *AVelemR, AVMid[50];
  real *RL_Fholder, *RR_Fholder;
  enum xfe_BasisType Basis;
  xf_Mesh *Mesh;
  xf_IFace IFace;
  
  Mesh   = All->Mesh;
  IFace  = Mesh->IFace[iiface];

  dim  = Model->dim;
  sr   = Model->nVars;
  sr2  = sr*sr;

  nL = PhiDataL->nn;
  nR = PhiDataR->nn;

  // residual orders (possibly different from state orders)
  // residual order should be identical to state order
  ResOrderL = ResPhiDataL->Order;
  ResOrderR = ResPhiDataR->Order;
  RnL = ResPhiDataL->nn;
  RnR = ResPhiDataR->nn;
  
  //element group info
  egrpL = IFace.ElemGroupL; egrpR = IFace.ElemGroupR;
  elemL = IFace.ElemL;      elemR = IFace.ElemR;

  // pull off temporary matrices
  T = DData->T;

  // wL is the vector that AL multiplies
  wL = DData->dunL;  // using this space as temporary storage
  for (k=0;k<dim*nq*sr;k++) wL[k] = guL[k]; // set wL = gradient

  if(Model->AVmodel)
  {
     AVelemL = Model->AVmodel_data->GenArray[egrpL].rValue[elemL];
     AVelemR = Model->AVmodel_data->GenArray[egrpR].rValue[elemR];
     for(k=0; k<sr+4; k++)
        AVMid[k] = 0.5 *(AVelemL[k] + AVelemR[k]);
  }
  else
  {
     AVelemL = NULL;
     AVelemR = NULL;
  }

  ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uL, guL, wL, DData->ALguL, DData->AL, gammaL, AVelemL));
  if (ierr != xf_OK) return ierr;
  //ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uL, guL, wL, DData->ALguL, DData->AL, gammaL, AVMid));
  //if (ierr != xf_OK) return ierr;

  // wR is the vector that AR multiplies
  wR = DData->dunR;
  for (k=0;k<dim*nq*sr;k++) wR[k] = guR[k]; // set wR = gradient

  ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uR, guR, wR, DData->ARguR, DData->AR, gammaR, AVelemR));
  if (ierr != xf_OK) return ierr;
  //ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uR, guR, wR, DData->ARguR, DData->AR, gammaR, AVMid));
  //if (ierr != xf_OK) return ierr;

  //!!no need to change for ghost state model
  // flux term first: Qn{q,k} = 0.5*(ALguL{i,q,k}+ARguR{i,q,k})*wn{q,i}
  xf_ColcMult_Set(DData->ALguL+0*nq*sr, wn+0, nq, sr, dim, 0.5, DData->Qn);
  xf_ColcMult_Add(DData->ARguR+0*nq*sr, wn+0, nq, sr, dim, 0.5, DData->Qn);
  for(i=1; i<dim; i++){
    xf_ColcMult_Add(DData->ALguL+i*nq*sr, wn+i, nq, sr, dim, 0.5, DData->Qn);
    xf_ColcMult_Add(DData->ARguR+i*nq*sr, wn+i, nq, sr, dim, 0.5, DData->Qn);
  }

  /* Derivatives of flux term */

  /*
   * Qn_guL{j,q,k;a} = 0.5*AL{i,j,q,k,a}*wn{i,q}
   * Qn_uL{q,k;a}    = 0.5*A_uLguL{i,q,k,a}*wn{i,q}
   */

  // Need u-flux: uhat = 0.5*(uL+uR) + alpha*(uL-uR)
  //ierr = xf_Error(xf_DiffFluxUJump(DiffDisc, &DData->alpha));
  //if (ierr != xf_OK) return ierr;
  DData->alpha = 0.0;

  /* Need jumps, so calculate:
     dunL = (uL - uhat)*nL = (0.5-alpha) * (uL - uR) * nL 
     dunR = (uR - uhat)*nR = (0.5+alpha) * (uR - uL) * nR = (0.5+alpha) * (uL - uR) * nL
     Note: quad weights are included

     Derivatives defined according to:
     dunL_uL = duL_uL*nL
     dunL_uR = duL_uR*nL
     dunR_uL = duR_uL*nL (yes, this is nL)
     dunR_uR = duR_uR*nL (same here)
   */

  //construct Puesdo-state for opposite side
  if(Model->GammaVaryFlag){
  //--------------------------------------------------------------------------------------------//
  ierr = xf_Error(xf_Alloc( (void **) &uLghost, nq*sr, sizeof(real)));
  if(ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &uRghost, nq*sr, sizeof(real)));
  if(ierr != xf_OK) return ierr;
  ierr = xf_Error(Yu_GhostStateConstruct(nq, sr, dim, uL, uR, uLghost, uRghost, gammaL, gammaR));
  if(ierr != xf_OK) return ierr;
  //--------------------------------------------------------------------------------------------//
  }

  DData->duL_uL =  (0.5-DData->alpha); DData->duL_uR = -(0.5-DData->alpha);
  DData->duR_uR = -(0.5+DData->alpha); DData->duR_uL =  (0.5+DData->alpha);
  for(i=0; i<dim; i++)
     for(iq=0; iq<nq; iq++)
        for (k=0, nval=wn[dim*iq+i]; k<sr; k++){
           
           if(Model->GammaVaryFlag){
              DData->dunL[(nq*i+iq)*sr+k] = DData->duL_uL*(uL[iq*sr+k]-uRghost[iq*sr+k])*nval;
              DData->dunR[(nq*i+iq)*sr+k] = DData->duR_uL*(uLghost[iq*sr+k]-uR[iq*sr+k])*nval;
           }
           else
           {
              DData->dunL[(nq*i+iq)*sr+k] = DData->duL_uL*(uL[iq*sr+k]-uR[iq*sr+k])*nval;
              DData->dunR[(nq*i+iq)*sr+k] = DData->duR_uL*(uL[iq*sr+k]-uR[iq*sr+k])*nval;
           }
        }

  // calculate ALdu [and A_uLdu]
  //ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uL, guL, DData->dunL, DData->ALdunL, NULL, gammaL, AVelemL));
  //if(ierr != xf_OK) return ierr;
  ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uL, guL, DData->dunL, DData->ALdunL, NULL, gammaL, AVMid));
  if(ierr != xf_OK) return ierr;

  //ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uR, guR, DData->dunR, DData->ARdunR, NULL, gammaR, AVelemR));
  //if(ierr != xf_OK) return ierr;
  ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uR, guR, DData->dunR, DData->ARdunR, NULL, gammaR, AVMid));
  if(ierr != xf_OK) return ierr;

  //what we do here is to compuate the maximum viscosity
  ierr = xf_Error(MaxViscosity(Model, nq, uL, uR, &max_vis));
  if(ierr != xf_OK) return ierr;
 
  // Calculate DData->N = normalized wn; and FaceArea
  for (iq=0, FaceArea=0.; iq<nq; iq++){
     for (i=0, nval=0.; i<dim; i++) nval += wn[dim*iq+i]*wn[dim*iq+i];
     FaceArea += (nval = sqrt(nval));
     for (i=0; i<dim; i++) DData->N[dim*iq+i] = wn[dim*iq+i]/nval;
  }

if(xfe_True)
{
    /*------------------------*/
    /* Viscous Discretization */
    /* Interior Penalty       */
    /*------------------------*/
    
    // need eta for IP stabilization
    Basis = min(Mesh->ElemGroup[egrpL].QBasis, Mesh->ElemGroup[egrpR].QBasis);
    Order = max(ResOrderL, ResOrderR); // accounts for possible residual p-dependence
    nFace = max(Mesh->ElemGroup[egrpL].nFace[elemL], Mesh->ElemGroup[egrpR].nFace[elemR]);
    //ierr = xf_Error(xf_JumpPenaltyEta(DiffDisc, Basis, Order, nFace, &eta));
    ierr = xf_Error(xf_JumpPenaltyEta(Basis, Order, nFace, &eta));
    if (ierr != xf_OK) return ierr;
 
    // ih = average (1/h) normal to face
    ihL = FaceArea/ElemVolL;
    ihR = FaceArea/ElemVolR;
    
    //  Qn{q,k} (see above)
/*    for (i=0; i<dim; i++){
        xf_ColcMult_Add(DData->ALdunL+sr*nq*i, DData->N+i, nq, sr, dim,
                        -eta*0.5*ihL, DData->Qn);
        xf_ColcMult_Add(DData->ARdunR+sr*nq*i, DData->N+i, nq, sr, dim,
                        -eta*0.5*ihR, DData->Qn);
    }
*/
    //compare Qn with inviscid dissipation coefficient
    ierr = xf_Error(DissipationCorrect(nq, sr, dim, uL, uR, wn, gammaR, DData->Qn, DData->Qn_dp, max_vis*(real)eta*max(ihL, ihR)));
    if (ierr != xf_OK) return ierr;
    
    /* Add Q term to R:
     RL{n,k} -= PhiL{n,q}*Qn{q,k},  sum over q
     RR{n,k} += PhiR{n,q}*Qn{q,k},  sum over q (+ because nR = -nL)
     */
    if (RL != NULL)
        xf_MTxM_Sub(PhiDataL->Phi, DData->Qn, nL, nq, sr, RL);
    if (RR != NULL)
        xf_MTxM_Add(PhiDataR->Phi, DData->Qn, nR, nq, sr, RR);
}
    
    //temporarily clip this for turbulence test
if(xfe_False)
    {
  /*------------------------*/
  /*Viscous Discretization  */
  /* BR2 formulation        */
  /*------------------------*/

  /* BR2 = Second form of Bassi & Rebay discretization
   *
   * Stabilization delta terms:
   *
   * Qn{q,k}   -= eta*dn{q,k} 
   * dn{q,k}    = 0.5*(dnL{q,k} + dnR{q,k})
   * dnL{q,k}   = ResPhiL{q,n}*DnL{i,n,k}*wn{i,q}
   * dnR{q,k}   = ResPhiR{q,n}*DnR{i,n,k}*wn{i,q}
   * DnL{i,n,k} = ResiML{n,m} * 0.5*ResPhiL{g,m}*ALdunL{i,g,k}
   * DnR{i,n,k} = ResiMR{n,m} * 0.5*ResPhiR{g,m}*ARdunR{i,g,k}
   */

  /* First add existing Qn to RL and RR, since we will be dealing directly with RL and RR.
   * Note that derivatives like Qn_uL, Qn_guL, etc. have already been added above */

  //July 2015 change 
  //if (RL != NULL)
  //   xf_MTxM_Sub(PhiDataL->Phi, DData->Qn, nL, nq, sr, RL);
  //if (RR != NULL)
  //   xf_MTxM_Add(PhiDataR->Phi, DData->Qn, nR, nq, sr, RR);
  RL_Fholder = DData->RL_Fholder;
  RR_Fholder = DData->RR_Fholder;
  for(i=0; i<nL*sr; i++) RL_Fholder[i] = 0.;
  xf_MTxM_Sub(PhiDataL->Phi, DData->Qn, nL, nq, sr, RL_Fholder);
  for(i=0; i<nR*sr; i++) RR_Fholder[i] = 0.;
  xf_MTxM_Add(PhiDataR->Phi, DData->Qn, nR, nq, sr, RR_Fholder);

  /* Now can deal with delta terms. First pull off inverse mass matrices,
   * iML, iMR. These are at the residual order. */
  ierr = xf_Error(xf_ElemInvMassMatrix(All, egrpL, elemL, BasisL,
                                       ResOrderL, NULL, NULL, &iMML, &facL));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ElemInvMassMatrix(All, egrpR, elemR, BasisR,
                                       ResOrderR, NULL, NULL, &iMMR, &facR));
  if (ierr != xf_OK) return ierr;

  /* Calculate delta coefficients DnL, DnR (note use of facL and facR for mass matrix scaling) */
  for (i=0; i<dim; i++){
     xf_MTxM_Set(ResPhiDataL->Phi, DData->ALdunL+i*nq*sr, RnL, nq, sr, T);
     xf_MxM_Set(iMML, T, RnL, RnL, sr, DData->DnL+i*RnL*sr);
  }
  for (k=0; k<dim*RnL*sr; k++) DData->DnL[k] *= 0.5*facL;
  for (i=0; i<dim; i++){
     xf_MTxM_Set(ResPhiDataR->Phi, DData->ARdunR+i*nq*sr, RnR, nq, sr, T);
     xf_MxM_Set(iMMR, T, RnR, RnR, sr, DData->DnR+i*RnR*sr);
  }
  for (k=0; k<dim*RnR*sr; k++) DData->DnR[k] *= 0.5*facR;

  //Determine eta
  nFace = max(Mesh->ElemGroup[egrpL].nFace[elemL], Mesh->ElemGroup[egrpR].nFace[elemR]);
  eta = ((real) nFace);
  // ihL and ihR used to set coupling
  ihL = FaceArea/ElemVolL;
  ihR = FaceArea/ElemVolR;

  /* Calculate
       SLL{i,n,m} = 0.5*eta*PhiL{n,q}*ResPhiL{q,m}*wn{i,q}
       SLR{i,n,m} = 0.5*eta*PhiL{n,q}*ResPhiR{q,m}*wn{i,q}
       SRL{i,n,m} = 0.5*eta*PhiR{n,q}*ResPhiL{q,m}*wn{i,q}
       SRR{i,n,m} = 0.5*eta*PhiR{n,q}*ResPhiR{q,m}*wn{i,q}
   */

  n2 = max(nL, nR)*max(RnL, RnR);  //for ease of indexing into S**
  for (i=0; i<dim; i++){
     xf_ColcMult_Set(PhiDataL->Phi, wn+i, nq,  nL, dim, 0.5*eta, T); // row
     xf_MTxM_Set(T, ResPhiDataL->Phi, nL, nq, RnL, DData->SLL+i*n2); // col

     xf_ColcMult_Set(PhiDataL->Phi, wn+i, nq,  nL, dim, 0.5*eta, T); // row
     xf_MTxM_Set(T, ResPhiDataR->Phi, nL, nq, RnR, DData->SLR+i*n2); // col

     xf_ColcMult_Set(PhiDataR->Phi, wn+i, nq,  nR, dim, 0.5*eta, T); // row
     xf_MTxM_Set(T, ResPhiDataL->Phi, nR, nq, RnL, DData->SRL+i*n2); // col

     xf_ColcMult_Set(PhiDataR->Phi, wn+i, nq,  nR, dim, 0.5*eta, T); // row
     xf_MTxM_Set(T, ResPhiDataR->Phi, nR, nq, RnR, DData->SRR+i*n2); // col
  }

  /* Add to RL and RR
   * RL{n,k} += SLL{i,n,m}*DnL{i,m,k} + SLR{i,n,m}*DnR{i,m,k}
   * RR{n,k} -= SRL{i,n,m}*DnL{i,m,k} + SRR{i,n,m}*DnR{i,m,k}
   */

 // if(RL != NULL){
 //    for(i=0; i<dim; i++){
 //       xf_MxM_Add(DData->SLL+i*n2, DData->DnL+i*RnL*sr, nL, RnL, sr, RL);
 //       xf_MxM_Add(DData->SLR+i*n2, DData->DnR+i*RnR*sr, nL, RnR, sr, RL);
 //    }
 // }

 // if(RR != NULL){
 //    for (i=0; i<dim; i++){
 //       xf_MxM_Sub(DData->SRL+i*n2, DData->DnL+i*RnL*sr, nR, RnL, sr, RR);
 //       xf_MxM_Sub(DData->SRR+i*n2, DData->DnR+i*RnR*sr, nR, RnR, sr, RR);
 //    }
 // }
  for(i=0; i<dim; i++){
     xf_MxM_Add(DData->SLL+i*n2, DData->DnL+i*RnL*sr, nL, RnL, sr, RL_Fholder);
     xf_MxM_Add(DData->SLR+i*n2, DData->DnR+i*RnR*sr, nL, RnR, sr, RL_Fholder);
     xf_MxM_Sub(DData->SRL+i*n2, DData->DnL+i*RnL*sr, nR, RnL, sr, RR_Fholder);
     xf_MxM_Sub(DData->SRR+i*n2, DData->DnR+i*RnR*sr, nR, RnR, sr, RR_Fholder);
  }

  //one-time addition to element residual
  if (RL != NULL) {
     xf_V_Add(RL_Fholder, nL*sr, xfe_Add, RL);
  }

  if (RR != NULL) {
     xf_V_Add(RR_Fholder, nR*sr, xfe_Add, RR);
  }
  //check flux conservation
}
  if(Model->GammaVaryFlag){
     xf_Release( (void *) uLghost);
     xf_Release( (void *) uRghost);
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: LYDG_CreateDiffBCData
int
LYDG_CreateDiffBCData(xf_DiffBCData **pDData)
{
  int ierr;
  xf_DiffBCData *DData;
  
  ierr = xf_Error(xf_Alloc( (void **) pDData, sizeof(xf_DiffBCData), 1));
  if (ierr != xf_OK) return ierr;

  DData = (*pDData);

  DData->Need_Grad = xfe_True;
  DData->ConstAB   = xfe_False;
  DData->alpha     = 0.0;

  DData->dun       = NULL; 
  DData->du_uI     = NULL;
  DData->du_gb     = NULL;
  DData->ABgu      = NULL; 
  DData->ABdun     = NULL; 
  DData->AB        = NULL; 
  DData->A_uBgu    = NULL; 
  DData->A_uBdun   = NULL; 
  DData->AB_uIdun  = NULL;
  DData->N         = NULL; 
  DData->Qn        = NULL; 
  DData->Qn_uI     = NULL; 
  DData->Qn_guI    = NULL; 
  DData->Qn_gb     = NULL; 
  DData->Qn_ggb    = NULL;
  DData->AwStab    = NULL; 
  DData->guB       = NULL; 
  DData->guB_guI   = NULL; 
  DData->Aw        = NULL; 
  DData->T         = NULL; 
  DData->Dn        = NULL; 
  DData->S         = NULL; 
  DData->P         = NULL; 
  DData->Vw        = NULL; 

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: LYDG_ReAllocDiffBCData
int
LYDG_ReAllocDiffBCData(int sr, int nq, int dim, int nn, 
		     enum xfe_Bool Need_Grad, xf_DiffBCData *DData)
{
  int ierr;
  int sr2, Tsize;

  sr2 = sr*sr;
  DData->Need_Grad = Need_Grad;

  ierr = xf_Error(xf_ReAlloc( (void **) &DData->dun, dim*nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ABgu, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->ABdun, sr*nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn, sr*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->guB, dim*nq*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->N , nq*dim, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Dn, dim*nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // BR2 needs these even without R_U request
  Tsize = dim*nq*sr2;
  Tsize = max(Tsize, nn*nq);
  Tsize = max(Tsize, nn*sr2);
  Tsize = max(Tsize, nn*nn);
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->T, Tsize, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->S, dim*nn*nn, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->P, dim*nn*nq, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // for output calculation
  ierr = xf_Error(xf_ReAlloc( (void **) &DData->Vw, dim*nn*sr, sizeof(real)));
  if (ierr != xf_OK) return ierr;

	
  if (Need_Grad){ 
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->du_uI, nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->du_gb, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AB, sr2*nq*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uBgu , sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->A_uBdun, sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AB_uIdun, sr2*nq*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_uI, nq*nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_guI, dim*nq*sr2, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_gb, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Qn_ggb, dim*nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->AwStab, dim*nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->guB_guI, sr2*nq*dim*dim, sizeof(real)));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_ReAlloc( (void **) &DData->Aw, nq*sr, sizeof(real)));
    if (ierr != xf_OK) return ierr;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyDiffBCData
void
LYDG_DestroyDiffBCData(xf_DiffBCData *DData){

  xf_Release( (void *) DData->dun);
  xf_Release( (void *) DData->du_uI);     
  xf_Release( (void *) DData->du_gb);      
  xf_Release( (void *) DData->ABgu);   
  xf_Release( (void *) DData->ABdun);    
  xf_Release( (void *) DData->AB);  
  xf_Release( (void *) DData->A_uBgu); 
  xf_Release( (void *) DData->A_uBdun);
  xf_Release( (void *) DData->AB_uIdun);
  xf_Release( (void *) DData->N);
  xf_Release( (void *) DData->Qn);
  xf_Release( (void *) DData->Qn_uI);   
  xf_Release( (void *) DData->Qn_guI);
  xf_Release( (void *) DData->Qn_gb);   
  xf_Release( (void *) DData->Qn_ggb);
  xf_Release( (void *) DData->AwStab);  
  xf_Release( (void *) DData->guB);
  xf_Release( (void *) DData->guB_guI);
  xf_Release( (void *) DData->Aw);
  xf_Release( (void *) DData->T);
  xf_Release( (void *) DData->Dn);
  xf_Release( (void *) DData->S);
  xf_Release( (void *) DData->P);
  xf_Release( (void *) DData->Vw);

  xf_Release( (void *) DData);
}

/******************************************************************/
//   FUNCTION Definition:  LYDG_SetVisFluxBC
static int
LYDG_SetVisFluxBC(int *SetIndex, int nq, int sr, int dim, const real *guB,
                  real *Qn, real *Aw, enum xfe_Bool ZeroFlag)
{
   //account for BC-imposed values of viscous flux
   //Zero flag means we just want to zero out components that are set up
   int k, iq, i, l, ii, sr2;

   sr2 = sr*sr;

   for (k=0; k<sr; k++)
      if (SetIndex[k]){
         for (iq=0;iq<nq; iq++){
           if (Qn != NULL){
              if (ZeroFlag) Qn[iq*sr+k] = 0.0;
              else Qn[iq*sr+k] = guB[iq*sr+k];
           }
           if(Aw != NULL)
              Aw[iq*sr+k] = 0.0;
         }
      }

   return xf_OK;
}
/******************************************************************/
//   FUNCTION Definition:  LYDG_BoundaryViscTerm
//   remark: current no need for special treat on variable gamma
int
LYDG_BoundaryViscTerm(xf_All *All, Yu_Model *Model, int ibfgrp, int ibface, 
                      const xf_BasisData *PhiData, const xf_BasisData *ResPhiData, 
                      real ElemVol, enum xfe_BasisType Basis, int Order, int nq, 
                      const real *wn, const real *xglob, real *uI, real *uB, 
                      real *guI, real *ER, xf_DiffBCData *DData, const real gamma,
                      xf_OutputEvalData *OutputEval)
{ 
  int ierr, i, j, k, l, dim, ii, jj, sr, sr2;
  int egrp, elem, iq, nn, n, n2;
  int ResOrder, Rnn, *SetIndex, localBCflag;
  enum xfe_Bool VisFluxFlag, Identity_gu_guI = xfe_False;
  real *T, *gu, *gu_guI, *iMM, *w=NULL;
  real eta, fac, nval, FaceArea, ih, val, t; 
  xf_BFace BFace;
  xf_Mesh *Mesh;
  //for output post-processing
  ////////////////////////////
  real xfac;
  int *FluxMoments = NULL;
  real *FluxWeights = NULL, *EV_U = NULL, *EV_G = NULL, *Value = NULL;
  real *AVelem;
  ////////////////////////////
   
  if (OutputEval != NULL){
     Value       = OutputEval->Value;
     EV_U        = OutputEval->EV_U;
     EV_G        = OutputEval->EV_G;
     FluxWeights = OutputEval->FluxWeights;
     FluxMoments = OutputEval->FluxMoments;
  }
   
  Mesh = All->Mesh;
  dim  = Model->dim;
  sr   = Model->nVars;
  sr2  = sr*sr;

  //element group info
  BFace = Mesh->BFaceGroup[ibfgrp].BFace[ibface];
  egrp  = BFace.ElemGroup;
  elem  = BFace.Elem;

  //pull off temporary matrices
  T = DData->T;
  nn = PhiData->nn;

  //residual order (possibly different from state order)
  ResOrder = ResPhiData->Order;
  Rnn      = ResPhiData->nn;

  //check out BC type before setting viscous flux
  for(i=0; i<Model->nBCs; i++)
     if(strcmp(Mesh->BFaceGroup[ibfgrp].Title, Model->nameBCs[i]) == 0)
     {  localBCflag = i; break;  }

  if(i==Model->nBCs) {
     printf("Boundary name does not match model defining.\n");
     xf_Error(xf_BOUNDARY_CONDITION_ERROR);
  }

  //allocate whether set BC visous flux
  ierr = xf_Error(xf_Alloc((void **) &SetIndex, sr, sizeof(int)));
  if (ierr != xf_OK) return ierr;
  for (k=0; k<sr; k++) SetIndex[k] = 0;

  // obtain viscous BC flux or grad u, using eqnset
  // TODO: Normal passed in should be transformed, before and after, when Motion is on
  //       (OK for now because currently no visc flux BCs explicitly use the normal)
  ierr = xf_Error(VisFluxBoundaryState(nq, sr, dim, wn, xglob, guI, DData->guB, &VisFluxFlag, 
                                       SetIndex, Model->typeBCs[localBCflag], Model->paraBCs[localBCflag], Model->AVmodel));
  if (ierr != xf_OK) return ierr;

  //VisFluxFlag == True means DData->guB is directly set on viscous flux, but gu is not
  //identical to guB; it might only be set on certain variable component;
  //VisFluxFlag == False means gu is directly set using DData->guB;

  if (VisFluxFlag){
    gu = guI;
    Identity_gu_guI = xfe_True;
   }
  else{
    Identity_gu_guI = xfe_False;
    gu = DData->guB;
    gu_guI = DData->guB_guI;
  }

  // w is the vector that AB multiplies
  w = DData->dun;  // using this space as temporary storage
  for (k=0;k<dim*nq*sr;k++) w[k] = gu[k]; // set w = gradient
  
  // Calculate AB*gu, A_uB*gu ; uB is now physical
  //ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm, nDiff, IParam, RParam, ResMetric, 
  //				  StabVisc, nq, uB, guI, w, xfe_False, MD, DData->ABgu, 
  //				  DData->AB, DData->A_uBgu, DData->AwStab, 
  //				  &DData->ConstAB, CF));
  if(Model->AVmodel)
  {   
     AVelem = Model->AVmodel_data->GenArray[egrp].rValue[elem];
  }
  else
     AVelem = NULL;

  ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uB, guI, w, DData->ABgu, 
                  DData->AB, gamma, AVelem));
  if (ierr != xf_OK) return ierr;

  // flux term first: Qn{q,k} = ABgu{i,q,k}*wn{q,i}
  xf_ColMult_Set(DData->ABgu+0*nq*sr, wn+0, nq, sr, dim, DData->Qn);
  for (i=1; i<dim; i++)
    xf_ColMult_Add(DData->ABgu+i*nq*sr, wn+i, nq, sr, dim, DData->Qn);
  
  // need bc flux in u: uhat = uB + alpha*(uI-uB);
  //ierr = xf_Error(xf_DiffFluxUBC(DiffDisc, &DData->alpha));
  //if (ierr != xf_OK) return ierr;
  DData->alpha = 0.0;

  // need jump in u: dun = (uI - uhat)*wn = (1-alpha)*(uI-uB)*wn
  fac = 1.-DData->alpha;
  for (i=0; i<dim; i++)
    for (iq=0;iq<nq; iq++)
      for (k=0, nval=wn[dim*iq+i]; k<sr; k++)
	DData->dun[(nq*i+iq)*sr+k] = fac*(uI[iq*sr+k]-uB[iq*sr+k])*nval;
  
  // Calculate AB, AB*dun, A_uB*dun ; uB is now physical
  //ierr = xf_Error(xf_ComputeDiffA(EqnSet, ResTerm, nDiff, IParam, RParam, ResMetric,
  //				  StabVisc, nq, uB, guI, DData->dun, xfe_False, MD, 
  //				  DData->ABdun, NULL, DData->A_uBdun, DData->AwStab, 
  //				  &DData->ConstAB, NULL));
  ierr = xf_Error(LYDG_InteriorViscTerm(Model, 1, nq, uB, guI, DData->dun, DData->ABdun, 
                  NULL, gamma, AVelem));
  if (ierr != xf_OK) return ierr;
 
  // Calculate DData->N = normalized wn; and FaceArea
  for (iq=0, FaceArea=0.; iq<nq; iq++){
    for (i=0, nval=0.; i<dim; i++) nval += wn[dim*iq+i]*wn[dim*iq+i];
    FaceArea += (nval = sqrt(nval));
    for (i=0; i<dim; i++) DData->N[dim*iq+i] = wn[dim*iq+i]/nval;
  }

  /********************/
  /*   BR2 Diffusion  */
  /********************/
  
    /* BR2 = Second form of Bassi & Rebay discretization

       Stabilization delta terms:

       Qn{q,k}   -= eta*dn{q,k}
       dn{q,k}    = ResPhi{q,n}*Dn{i,n,k}*wn{i,q}
       Dn{i,n,k}  = ResiM{n,m} * ResPhi{q,m}*ABdun{i,q,k}
       Dn_UI{i,n,k;o,a} = ResiM{n,m} * ResPhi{q,m}*(AB_uIdun{i,q,k;a} + AB{i,j,q,k,l}*du_uI{q,l;a}*wn{q,j})*Phi{q,o}
                                                    \_____________ ABdun_uI{i,q,k;a} ____________________/
       Since,
       ER{n,k} += Phi{n,q}*eta*ResPhi{q,m}*Dn{i,m,k}*wn{i,q} = S{i,n,m}*Dn{i,m,k}
       then
       ER_UI{n,k;o,a} += S{i,n,m}*Dn_UI{i,m,k;o,a}

       N{i,q} = normalized wn{i,q}
       Note, in using the mass matrix, iMM must be multiplied by fac

       Note, ResPhi and ResiM refer to quantities computed at a
       possibly different order than the state approximation.  This is
       useful for error estimation when the order of the residual
       needs to remain constant.

     */

    /* First add existing Qn to ER, since we will be dealing
       directly with ER.  Note that derivatives like Qn_uI
       Qn_guI, etc. have already been added above. */

    // Store flux components set by BC in Qn
    if (VisFluxFlag){
      ierr = xf_Error(LYDG_SetVisFluxBC(SetIndex, nq, sr, dim, DData->guB, DData->Qn, NULL,
					xfe_False));
      if (ierr != xf_OK) return ierr;
    }

    /* Add Q term to R:
       ER{n,k} -= Phi{n,q}*Qn{q,k},  sum over q
    */
    if (ER != NULL)
      xf_MTxM_Sub(PhiData->Phi, DData->Qn, nn, nq, sr, ER);

    //deal with output for logging or post-processing
    /////////////////////////////////
    if (Value != NULL){
       for (iq=0; iq<nq; iq++)
          for (k=0; k<sr; k++){
             xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
             (*Value) -= FluxWeights[k]*DData->Qn[iq*sr+k]*xfac;
          }
    }
   
    /* Now can deal with delta terms.  First pull off inverse mass
       matrix, iMM. This is at the residual order. */
    ierr = xf_Error(xf_ElemInvMassMatrix(All, egrp, elem, Basis, ResOrder, 
                                         NULL, NULL, &iMM, &fac));
    if (ierr != xf_OK) return ierr;
    
    /* Calculate delta coefficients Dn (note use of fac for mass matrix scaling) */
    for (i=0; i<dim; i++){
      xf_MTxM_Set(ResPhiData->Phi, DData->ABdun+i*nq*sr, Rnn, nq, sr, T);
      xf_MxM_Set(iMM, T, Rnn, Rnn, sr, DData->Dn+i*Rnn*sr);
    }
    for (k=0; k<dim*Rnn*sr; k++) DData->Dn[k] *= fac;

    // delta is zero for flux components prescribed by BC -> zero out appropriate cols of Dn
    if (VisFluxFlag){
      ierr = xf_Error(LYDG_SetVisFluxBC(SetIndex, dim*Rnn, sr, dim, NULL, DData->Dn, NULL, 
                                          xfe_True));
      if (ierr != xf_OK) return ierr;
    }
    
    // Determine eta
    eta = 0.5*((real) Mesh->ElemGroup[egrp].nFace[elem]);
    // ih used to set coupling
    ih = FaceArea/ElemVol;

    /* Calculate:
         S{i,n,m} = eta*Phi{n,q}*ResPhi{q,m}*wn{i,q}
    */

    n2 = nn*Rnn;
    for (i=0; i<dim; i++){
      xf_ColcMult_Set(PhiData->Phi, wn+i, nq,  nn, dim, eta, T);   // row
      xf_MTxM_Set(T, ResPhiData->Phi, nn, nq, Rnn, DData->S+i*n2); // col
    }    

    /* Add to Qn (so that the returned Qn is actually a valid flux):
         Qn{q,k}   -= eta*dn{q,k}
         dn{q,k}    = ResPhi{q,n}*Dn{i,n,k}*wn{i,q}
    */
    for (i=0; i<dim; i++){
      // T{q,k}  = ResPhi{q,n}*Dn{i,n,k}
      xf_MxM_Set(ResPhiData->Phi, DData->Dn+i*Rnn*sr, nq, Rnn, sr, T);
      // Qn{q,k} += -etaT{q,k}*wn{i,q}
      xf_ColcMult_Add(T, wn+i, nq, sr, dim, -eta, DData->Qn);
    } 

    /* Add to ER:
         ER{n,k} += S{i,n,m}*Dn{i,m,k}
    */
    if (ER != NULL) 
      for (i=0; i<dim; i++)
      	xf_MxM_Add(DData->S+i*n2, DData->Dn+i*Rnn*sr, nn, Rnn, sr, ER);

    //deal with output for logging or post-processing
    /////////////////////////////////
    if (Value != NULL){
       for (i=0; i<dim; i++){
          for (iq=0; iq<nq; iq++)
             for (k=0; k<sr; k++){
                xfac = ((FluxMoments[k] < 0) ? 1.0 : xglob[iq*dim+FluxMoments[k]]);
                T[iq*sr + k] = wn[iq*dim+i]*xfac*eta*FluxWeights[k];
             }
          xf_MTxM_Set(PhiData->Phi, T, Rnn, nq, sr, DData->Vw+i*Rnn*sr);
       }
       
       xf_DotProduct(DData->Vw, DData->Dn, dim*Rnn*sr, &val);
       (*Value) += val;
    }
      
        
    xf_Release((void*) SetIndex);

  return xf_OK;
}
