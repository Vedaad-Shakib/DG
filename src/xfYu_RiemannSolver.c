//include all the Riemann solver in this routine

#define MAXSR 50


/******************************************************************/
// Euler Flux Evaluate for 2D setup
static int
ConvFluxInterior2D(const int sr, const real *U, real *F, real *G, const real Gamma)
{
   int i, j;
   real rho, u, v, p;
   
   rho = U[0];
   u   = U[1] / U[0];
   v   = U[2] / U[0];
   p   = (U[3] - 0.5 * (U[1]*U[1] + U[2]*U[2])/U[0]) * (Gamma - 1.0);
   
   if(LinearCase)
   {
      for(i=0; i<sr; i++)
      { F[i] = 0.; G[i] = 0.;}
      
      //need define advection speed; default set to unit
      F[0] = AdvSpd*rho; 
      
   }
   else{
      
      if(rho <0)
      {
         printf("Density goes to negative!\n");
         return xf_Error(xf_NON_PHYSICAL);
      }
      
      if(p <0)
      {
         printf("Pressure goes to negative!\n");
         return xf_Error(xf_NON_PHYSICAL);
      }
      
      F[0] = rho*u;
      F[1] = rho*u*u + p;
      F[2] = rho*u*v;
      F[3] = u*(U[3] + p);
      
      G[0] = rho*v;
      G[1] = rho*v*u;
      G[2] = rho*v*v + p;
      G[3] = v*(U[3] + p);
      
      for(i=4; i<sr; i++)
      {
         F[i] = U[i]*u;
         G[i] = U[i]*v;
      }
   }
   return xf_OK;
}

/******************************************************************/
// Euler Flux Evaluate for 3D setup
static int
ConvFluxInterior3D(const int sr, const real *U, real *F, real *G, real *H, const real Gamma)
{
   int i, j;
   real rho, u, v, w, p;
   
   rho = U[0];
   u   = U[1] / U[0];
   v   = U[2] / U[0];
   w   = U[3] / U[0];
   p   = (U[4] - 0.5 * (U[1]*U[1] + U[2]*U[2] + U[3]*U[3])/U[0]) * (Gamma - 1.0);
   
   if(LinearCase)
   {
      for(i=0; i<sr; i++)
      { F[i] = 0.; G[i] = 0.; H[i] = 0.;}
      
      //need define advection speed; default set to unit
      F[0] = AdvSpd*rho; 
      
   }
   else{
      
      F[0] = rho*u;
      F[1] = rho*u*u + p;
      F[2] = rho*u*v;
      F[3] = rho*u*w;
      F[4] = u*(U[4] + p);
      
      G[0] = rho*v;
      G[1] = rho*v*u;
      G[2] = rho*v*v + p;
      G[3] = rho*v*w;
      G[4] = v*(U[4] + p);
      
      H[0] = rho*w;
      H[1] = rho*w*u;
      H[2] = rho*w*v;
      H[3] = rho*w*w + p;
      H[4] = w*(U[4] + p);
      
      for(i=5; i<sr; i++)
      {
         F[i] = U[i]*u;
         G[i] = U[i]*v;
         H[i] = U[i]*w;
      }
   }
   
   return xf_OK;
}

/******************************************************************/
// Entrance for Flux Evaluation
int
ConvFluxInterior(const int nq, const int sr, const int dim,
                 const real *U, real *F, const real Gamma)
{
   int ierr, i, off, off2, sr2;
   sr2 = sr * sr;
   off = sr * nq;
   off2 = sr2 * nq;
   
   if (dim == 2) {
      for (i=0; i<nq; i++){
         ierr = xf_Error(ConvFluxInterior2D(sr, U+i*sr, F+i*sr, F+off+i*sr, Gamma));
         if (ierr != xf_OK) return ierr;
      }
   }
   else if(dim == 3){
      for (i=0; i<nq; i++){
         ierr = xf_Error(ConvFluxInterior3D(sr, U+i*sr, F+i*sr, F+off+i*sr, F+2*off+i*sr, Gamma));
         if (ierr != xf_OK) return ierr;
      }
   }
   else
      return xf_Error(xf_OUT_OF_BOUNDS);
   
   
   return xf_OK;
}


/******************************************************************/
// Riemann Flux for 2D setup
static int
ConvFluxFace2D(const int sr, const real *UL, const real *UR, const real *n,
               real *F, const real Gamma, real *MaxCharSpeed)
{
   int i, j, k;
   real rhol, ul, vl, pl, al, unl, Fl[4];
   real rhor, ur, vr, pr, ar, unr, Fr[4];
   real cl, sssul, ssstl, c1l, c2l;
   real cr, sssur, ssstr, c1r, c2r;
   real rhom, am, ppvrs, ps, pspl, pspr, ql, qr, ssl, ssr, sss;
   real NN, NN1, n1[3];
   real fac, Machl, Machr;
   
   NN = sqrt(n[0]*n[0] + n[1]*n[1]);
   
   NN1 = 1/NN;
   
   for(i=0; i<2; i++)
      n1[i] = n[i]/NN;

if(LinearCase)
{
   for(i=0; i<sr; i++)
      F[i] = 0.;
   
   unl = AdvSpd*n1[0] + 0.*n1[1];
   unr = AdvSpd*n1[0] + 0.*n1[1];

   Fl[0] = UL[0] * unl;
   Fr[0] = UR[0] * unr;

   F[0] = 0.5*NN*((Fl[0] + Fr[0]) - AdvSpd*(UR[0] - UL[0]));
   *MaxCharSpeed = AdvSpd;
}
else
{
   //Left state
   rhol = UL[0];
   ul   = UL[1] / UL[0];
   vl   = UL[2] / UL[0];
   pl   = (UL[3] - 0.5*(UL[1]*UL[1] + UL[2]*UL[2])/UL[0])*(Gamma - 1.0);
   al   = sqrt(Gamma*pl/rhol);
   unl  = (ul*n1[0] + vl*n1[1]);
   Machl = fabs(unl) / al;
   
   Fl[0]   = rhol*unl;
   Fl[1]   = n1[0]*pl + rhol*ul*unl;
   Fl[2]   = n1[1]*pl + rhol*vl*unl;
   Fl[3]   = (UL[3] + pl)*unl;
   
   //Right state
   rhor = UR[0];
   ur   = UR[1] / UR[0];
   vr   = UR[2] / UR[0];
   pr   = (UR[3] - 0.5*(UR[1]*UR[1] + UR[2]*UR[2])/UR[0])*(Gamma - 1.0);
   ar   = sqrt(Gamma*pr/rhor);
   unr  = (ur*n1[0] + vr*n1[1]);
   Machr = fabs(unr) / ar;
   
   //Take maximum characteristic speed
   //(*MaxCharSpeed) = fabs(unl) + al;
   //if((*MaxCharSpeed) < (fabs(unr) + ar))
   //   (*MaxCharSpeed) = (fabs(unr) + ar);
   (*MaxCharSpeed) = sqrt(ul*ul + vl*vl) + al;
   if((*MaxCharSpeed) < (sqrt(ur*ur + vr*vr) + ar))
      (*MaxCharSpeed) = (sqrt(ur*ur + vr*vr) + ar);
   
   if(rhol <0 || rhor <0 || pl <0 || pr <0)
   {
      if(rhol <0 || rhor <0)
         printf("Density goes negative!");
      
      if(pl <0 || pr <0)
         printf("Pressure goes negative!");
      
      return xf_Error(xf_NON_PHYSICAL);
   }
   
   Fr[0]   = rhor*unr;
   Fr[1]   = n1[0]*pr + rhor*ur*unr;
   Fr[2]   = n1[1]*pr + rhor*vr*unr;
   Fr[3]   = (UR[3] + pr)*unr;
   
   rhom = 0.5*(rhol + rhor);
   am   = 0.5*(al + ar);
   
   // step1: estimate pressure at star region
   ppvrs = 0.5*((pr+pl) - (unr-unl)*rhom*am);
   
   if(ppvrs > 0.0)
      ps = ppvrs;
   else
      ps = 0.0;
   
   pspl = ps/pl;
   pspr = ps/pr;
   
   // step2: left and right wave speed
   ql = 1.0;
   if(pspl > 1.0)
      ql = sqrt(1.0+(Gamma+1.0)/(2.0*Gamma)*(pspl-1.0));
   ssl = unl - al*ql;
   
   qr = 1.0;
   if(pspr > 1.0)
      qr = sqrt(1.0+(Gamma+1.0)/(2.0*Gamma)*(pspr-1.0));
   ssr = unr + ar*qr;
   
   //step3: shear wave speed
   sss = 0.5*(unl + unr) + 0.5*(pl-pr)/(rhom*am);
   
   //step4: compute flux
   if(ssl >= 0.0){
      F[0] = NN*Fl[0];
      F[1] = NN*Fl[1];
      F[2] = NN*Fl[2];
      F[3] = NN*Fl[3];
      for(i=4; i<sr; i++)
         F[i] = NN*(UL[i]*unl);
   }
   else if(ssr<=0.0){
      F[0] = NN*Fr[0];
      F[1] = NN*Fr[1];
      F[2] = NN*Fr[2];
      F[3] = NN*Fr[3];
      for(i=4; i<sr; i++)
         F[i] = NN*(UR[i]*unr);
   }
   else if((ssl<=0.0)&&(sss>=0.0)){
      cl = (ssl - unl)/(ssl - sss);
      sssul = sss - unl;
      ssstl = sss + pl/rhol/(ssl - unl);
      c1l = rhol*cl*sssul;
      c2l = rhol*cl*sssul*ssstl;
      
      F[0] = NN*(Fl[0] + ssl*(rhol*(cl - 1.0)));
      F[1] = NN*(Fl[1] + ssl*(rhol*ul*(cl - 1.0) + c1l*n1[0]));
      F[2] = NN*(Fl[2] + ssl*(rhol*vl*(cl - 1.0) + c1l*n1[1]));
      F[3] = NN*(Fl[3] + ssl*(UL[3]*(cl - 1.0) + c2l));
      for(i=4; i<sr; i++)
         F[i] = NN*(UL[i]*unl + ssl*UL[i]*(cl - 1.0));
      
      //try to blend the central flux to reduce dissipation at Low Mach regime
      /*     if(Machl > Machr)
       fac = Machl;
       else
       fac = Machr;
       
       if(fac < 1.)
       {
       //hard-coded using AUSM idea (cut-off Mach number)
       if(fac < 0.2) fac = 0.2;
       fac = fac * (2. - fac) * 2.;
       
       for(i=0; i<4; i++)
       F[i] = (1.- fac) * 0.5 * NN * (Fl[i] + Fr[i]) + fac * F[i];
       
       for(i=4; i<sr; i++)
       F[i] = (1.- fac) * 0.5 * NN * (UL[i]*unl + UR[i]*unr) + fac * F[i];
       }
       
       */
   }
   else if((sss<=0.0)&&(ssr>=0.0)){
      cr = (ssr - unr)/(ssr - sss);
      sssur = sss - unr;
      ssstr = sss + pr/rhor/(ssr - unr);
      c1r = rhor*cr*sssur;
      c2r = rhor*cr*sssur*ssstr;
      
      F[0] = NN*(Fr[0] + ssr*(rhor*(cr - 1.0)));
      F[1] = NN*(Fr[1] + ssr*(rhor*ur*(cr - 1.0) + c1r*n1[0]));
      F[2] = NN*(Fr[2] + ssr*(rhor*vr*(cr - 1.0) + c1r*n1[1]));
      F[3] = NN*(Fr[3] + ssr*(UR[3]*(cr - 1.0) + c2r));
      for(i=4; i<sr; i++)
         F[i] = NN*(UR[i]*unr + ssr*UR[i]*(cr - 1.0));
      
      //try to blend the central flux to reduce dissipation at Low Mach regime
      //currently not used for numerical test
      /*
       if(Machl > Machr)
       fac = Machl;
       else
       fac = Machr;
       
       if(fac < 1.)
       {
       //hard-coded using AUSM idea (cut-off Mach number)
       if(fac < 0.2) fac = 0.2;
       fac = fac * (2. - fac) * 2.;
       for(i=0; i<4; i++)
       F[i] = (1.- fac) * 0.5 * NN * (Fl[i] + Fr[i]) + fac * F[i];
       
       for(i=4; i<sr; i++)
       F[i] = (1.- fac) * 0.5 * NN * (UL[i]*unl + UR[i]*unr) + fac * F[i];
       }
       */
   }
   else
   {
      printf("%lf %lf %lf %lf\n", rhor, rhol, pr, pl);
      return xf_Error(xf_OUT_OF_BOUNDS);
   }
}   
   //special case for noslipwall
   /*if(fabs(ssr + ssl) < 1.e-15)
    {
    F[0] = 0.0; F[1] = 0.0; F[2] = 0.0; F[3] = 0.0;
    for(i=4; i<sr; i++)
    F[i] = 0.0;
    }*/
   
   return xf_OK;
}
/******************************************************************/
// Riemann Flux for 3D setup
static int
ConvFluxFace3D(const int sr, const real *UL, const real *UR, const real *n,
               real *F, const real Gamma, real *MaxCharSpeed)
{
   
   int i, j, k;
   real rhol, ul, vl, wl, pl, al, unl, Fl[5];
   real rhor, ur, vr, wr, pr, ar, unr, Fr[5];
   real cl, sssul, ssstl, c1l, c2l;
   real cr, sssur, ssstr, c1r, c2r;
   real rhom, am, ppvrs, ps, pspl, pspr, ql, qr, ssl, ssr, sss;
   real NN, NN1, n1[3];
   
   NN = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   NN1 = 1/NN;
   
   for(i=0; i<3; i++)
      n1[i] = n[i]/NN;
   
if(LinearCase)
{
   for(i=0; i<sr; i++)
      F[i] = 0.;
   
   unl = AdvSpd*n1[0] + 0.*n1[1] + 0.*n1[2];
   unr = AdvSpd*n1[0] + 0.*n1[1] + 0.*n1[2];

   Fl[0] = UL[0] * unl;
   Fr[0] = UR[0] * unr;

   F[0] = 0.5*NN*((Fl[0] + Fr[0]) - AdvSpd*(UR[0] - UL[0]));
   *MaxCharSpeed = AdvSpd;
}
else
{
   //Left state
   rhol = UL[0];
   ul   = UL[1] / UL[0];
   vl   = UL[2] / UL[0];
   wl   = UL[3] / UL[0];
   pl   = (UL[4] - 0.5*(UL[1]*UL[1] + UL[2]*UL[2] + UL[3]*UL[3])/UL[0])*(Gamma - 1.0);
   al   = sqrt(Gamma*pl/rhol);
   unl  = (ul*n1[0] + vl*n1[1] + wl*n1[2]);
   
   Fl[0]   = rhol*unl;
   Fl[1]   = n1[0]*pl + rhol*ul*unl;
   Fl[2]   = n1[1]*pl + rhol*vl*unl;
   Fl[3]   = n1[2]*pl + rhol*wl*unl;
   Fl[4]   = (UL[4] + pl)*unl;
   
   //Right state
   rhor = UR[0];
   ur   = UR[1] / UR[0];
   vr   = UR[2] / UR[0];
   wr   = UR[3] / UR[0];
   pr   = (UR[4] - 0.5*(UR[1]*UR[1] + UR[2]*UR[2] + UR[3]*UR[3])/UR[0])*(Gamma - 1.0);
   ar   = sqrt(Gamma*pr/rhor);
   unr  = (ur*n1[0] + vr*n1[1] + wr*n1[2]);
   
   Fr[0]   = rhor*unr;
   Fr[1]   = n1[0]*pr + rhor*ur*unr;
   Fr[2]   = n1[1]*pr + rhor*vr*unr;
   Fr[3]   = n1[2]*pr + rhor*wr*unr;
   Fr[4]   = (UR[4] + pr)*unr;
   
   //(*MaxCharSpeed) = fabs(unl) + al;
   //if((*MaxCharSpeed) < (fabs(unr) + ar))
   //   (*MaxCharSpeed) = fabs(unr) + ar;
   
   (*MaxCharSpeed) = sqrt(ul*ul + vl*vl + wl*wl) + al;
   am = sqrt(ur*ur + vr*vr + wr*wr) + ar;
   if((*MaxCharSpeed) < am) 
      (*MaxCharSpeed) = am; 
   
   //physical sense check
   if(rhol <0 || rhor <0 || pl <0 || pr <0)
   {
      if(rhol <0 || rhor <0)
         printf("Density goes negative!");
      
      if(pl <0 || pr <0)
         printf("Pressure goes negative!");
      
      return xf_Error(xf_NON_PHYSICAL);
   }
   
   rhom = 0.5*(rhol + rhor);
   am   = 0.5*(al + ar);
   
   // step1: estimate pressure at star region
   ppvrs = 0.5*((pr+pl) - (unr-unl)*rhom*am);
   
   if(ppvrs > 0.0)
      ps = ppvrs;
   else
      ps = 0.0;
   
   pspl = ps/pl;
   pspr = ps/pr;
   
   // step2: left and right wave speed
   ql = 1.0;
   if(pspl > 1.0)
      ql = sqrt(1.0+(Gamma+1.0)/(2.0*Gamma)*(pspl-1.0));
   ssl = unl - al*ql;
   
   qr = 1.0;
   if(pspr > 1.0)
      qr = sqrt(1.0+(Gamma+1.0)/(2.0*Gamma)*(pspr-1.0));
   ssr = unr + ar*qr;
   
   //step3: shear wave speed
   sss = 0.5*(unl + unr) + 0.5*(pl-pr)/(rhom*am);
   
   //step4: compute flux
   if(ssl >= 0.0){
      F[0] = NN*Fl[0];
      F[1] = NN*Fl[1];
      F[2] = NN*Fl[2];
      F[3] = NN*Fl[3];
      F[4] = NN*Fl[4];
      for(i=5; i<sr; i++)
         F[i] = NN*(UL[i]*unl);
   }
   else if(ssr<=0.0){
      F[0] = NN*Fr[0];
      F[1] = NN*Fr[1];
      F[2] = NN*Fr[2];
      F[3] = NN*Fr[3];
      F[4] = NN*Fr[4];
      for(i=5; i<sr; i++)
         F[i] = NN*(UR[i]*unr);
   }
   else if((ssl<=0.0)&&(sss>=0.0)){
      cl = (ssl - unl)/(ssl - sss);
      sssul = sss - unl;
      ssstl = sss + pl/rhol/(ssl - unl);
      c1l = rhol*cl*sssul;
      c2l = rhol*cl*sssul*ssstl;
      
      F[0] = NN*(Fl[0] + ssl*(rhol*(cl - 1.0)));
      F[1] = NN*(Fl[1] + ssl*(rhol*ul*(cl - 1.0) + c1l*n1[0]));
      F[2] = NN*(Fl[2] + ssl*(rhol*vl*(cl - 1.0) + c1l*n1[1]));
      F[3] = NN*(Fl[3] + ssl*(rhol*wl*(cl - 1.0) + c1l*n1[2]));
      F[4] = NN*(Fl[4] + ssl*(UL[4]*(cl - 1.0) + c2l));
      for(i=5; i<sr; i++)
         F[i] = NN*(UL[i]*unl + ssl*UL[i]*(cl - 1.0));
   }
   else if((sss<=0.0)&&(ssr>=0.0)){
      cr = (ssr - unr)/(ssr - sss);
      sssur = sss - unr;
      ssstr = sss + pr/rhor/(ssr - unr);
      c1r = rhor*cr*sssur;
      c2r = rhor*cr*sssur*ssstr;
      
      F[0] = NN*(Fr[0] + ssr*(rhor*(cr - 1.0)));
      F[1] = NN*(Fr[1] + ssr*(rhor*ur*(cr - 1.0) + c1r*n1[0]));
      F[2] = NN*(Fr[2] + ssr*(rhor*vr*(cr - 1.0) + c1r*n1[1]));
      F[3] = NN*(Fr[3] + ssr*(rhor*wr*(cr - 1.0) + c1r*n1[2]));
      F[4] = NN*(Fr[4] + ssr*(UR[4]*(cr - 1.0) + c2r));
      for(i=5; i<sr; i++)
         F[i] = NN*(UR[i]*unr + ssr*UR[i]*(cr - 1.0));
   }
   else
   {
      printf("%lf %lf %lf %lf\n", rhor, rhol, pr, pl);
      return xf_Error(xf_OUT_OF_BOUNDS);
   }
}

   return xf_OK;
}

/******************************************************************/
//Roe Riemann flux for 2D setup
static int
ConvFluxFace2D_Roe(const int sr, const real *UL, const real *UR,
                   const real *n, const real *vg, const real *RParam, real *F,
                   const real Gamma, real *MaxCharSpeed)
{
   /*
    PURPOSE:
    
    Evaluates 2d Roe-averaged approximate Riemann flux.
    
    INPUTS:
    
    PIS: Position In State vector;: maps enumerated type to state index
    UL: left state vector
    UR: left state vector
    n: normal vector (not necessarily of unit magnitude)
    vg : interface velocity
    RParam: function real parameters
    
    OUTPUTS:
    
    F : flux dotted with n
    
    RETURN:
    
    Error Code
    */
   real rhol1, rhor1;
   real ul, vl;
   real ur, vr;
   real u2l;
   real u2r;
   real unl;
   real unr;
   real pl;
   real hl;
   real pr;
   real hr;
   real FL[MAXSR];
   real FR[MAXSR];
   
   real di;
   real d1;
   real ui;
   real vi;
   real hi;
   real af;
   real ucp;
   real c2;
   real ci;
   
   real du[MAXSR], ci1, c21, el[3], ep;
   
   real lam[3];
   
   real eps;
   real eps1;
   
   real s1;
   real s2;
   real l3;
   real G1;
   real G2;
   real C1;
   real C2;
   
   real al, ar;
   real Mar, Mal, Ma;
   enum xfe_Bool LM_Fix_flag;
   //!!!!!!!!!!
   LM_Fix_flag = xfe_True;
   
   real gmi, gmi1;
   real NN, NN1, n1[3];
   
   int i, ir, iru, irv, irE;
   int iP, nPass, iPass[MAXSR];
   
   real kl, kr, drk, ki, vgn;
   
   //variable index
   //ir  = PIS[xfe_Density];
   //iru = PIS[xfe_XMomentum];
   //irv = PIS[xfe_YMomentum];
   //irE = PIS[xfe_Energy];
   ir = 0; iru = 1; irv = 2; irE = 3;
   
   gmi = Gamma - 1.0;
   gmi1 = 1/gmi;
   
   NN = sqrt(n[0]*n[0] + n[1]*n[1]);
   NN1 = 1/NN;
   
   for (i=0; i<2; i++)
      n1[i] = n[i]*NN1;
   
   if((UL[ir] <= 0)||(UR[ir]<0))
      return xf_Error(xf_NON_PHYSICAL);
   
   // interface normal velocity
   vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1];
   
   // Left State
   rhol1 = 1/UL[ir];
   
   ul    = UL[iru]*rhol1;
   vl    = UL[irv]*rhol1;
   u2l   = (ul*ul + vl*vl)*UL[ir];
   unl   = (ul*n1[0] + vl*n1[1]);
   pl    = (UL[irE] - 0.5*u2l  )*gmi;
   hl    = (UL[irE] + pl  )*rhol1;
   al    = sqrt(pl*Gamma*rhol1);
   Mal   =  (fabs(ul)+fabs(ur)) / al;
   
   FL[ir]     = UL[ir]*unl;
   FL[iru]    = n1[0]*pl   + UL[iru]*unl;
   FL[irv]    = n1[1]*pl   + UL[irv]*unl;
   FL[irE]    = (pl   + UL[irE])*unl;
   
   // Right State
   rhor1 = 1/UR[ir];
   
   ur    = UR[iru]*rhor1;
   vr    = UR[irv]*rhor1;
   u2r   = (ur*ur + vr*vr)*UR[ir];
   unr   = (ur*n1[0] + vr*n1[1]);
   pr    = (UR[irE] - 0.5*u2r  )*gmi;
   hr    = (UR[irE] + pr  )*rhor1;
   ar    = sqrt(pr*Gamma*rhor1);
   Mar   = (fabs(ur)+fabs(ul)) / ar;
   
   FR[ir ]    = UR[ir]*unr;
   FR[iru]    = n1[0]*pr   + UR[iru]*unr;
   FR[irv]    = n1[1]*pr   + UR[irv]*unr;
   FR[irE]    = (pr   + UR[irE])*unr;
   
   //estimate maximum wave speed
   (*MaxCharSpeed) = sqrt(ul*ul + vl*vl) + al;
   if((*MaxCharSpeed) < (sqrt(ur*ur + vr*vr) + ar))
      (*MaxCharSpeed) = (sqrt(ur*ur + vr*vr) + ar);
   
   // Average state
   di     = sqrt(UR[ir]*rhol1);
   d1     = 1.0/(1.0+di);
   
   ui     = (di*ur + ul)*d1;
   vi     = (di*vr + vl)*d1;
   hi     = (di*hr+hl)*d1;
   
   af     = 0.5*(ui*ui   +vi*vi );
   ucp    = ui*n1[0] + vi*n1[1];
   c2     = gmi*(hi   -af   );
   
   if(c2 <= 0)
      return xf_Error(xf_NON_PHYSICAL);
   
   ci    = sqrt(c2);
   ci1   = 1/ci;
   Ma = (fabs(ui) + fabs(vi)) / ci;
   
   // du = UR-UL
   du[ir ] = UR[ir ] - UL[ir ];
   du[iru] = UR[iru] - UL[iru];
   du[irv] = UR[irv] - UL[irv];
   du[irE] = UR[irE] - UL[irE];
   
   // eigenvalues
   lam[0] = ucp-vgn+ ci;
   lam[1] = ucp-vgn- ci;
   lam[2] = ucp-vgn;
   
   // Entropy fix
   ep     = 1e-2;
   eps    = ep*ci;
   
   for (i=0; i<3; i++)
      if ((lam[i]<eps)&&(lam[i]>-eps)){
         eps1 = 1/eps;
         
         lam[i] = 0.5*(eps+lam[i]*lam[i]*eps1);
      }
   
   // define el = sign(lam[i])
   for (i=0; i<3; i++)
      if (lam[i]<0)
         el[i] = -1;
      else
         el[i] =  1;
   
   // average and half-difference of 1st and 2nd eigs
   s1    = 0.5*(el[0]*lam   [0]+el[1]*lam   [1]);
   s2    = 0.5*(el[0]*lam   [0]-el[1]*lam   [1]);
   
   // third eigenvalue, absolute value
   l3    = el[2]*lam[2];
   
   // left eigenvector product generators (see Theory guide)
   G1    = gmi*(af*du[ir] - ui*du[iru] - vi*du[irv] + du[irE]);
   G2    = -ucp*du[ir]+du[iru]*n1[0]+du[irv]*n1[1];
   
   if(LM_Fix_flag){
   if(Ma>=1.0 || Mal>=1.0 || Mar>=1.0)
      Ma = 1.0;
   G2 *= Ma;
   }
   
   // required functions of G1 and G2 (again, see Theory guide)
   C1    = G1*(s1-l3)*ci1*ci1 + G2*s2*ci1;
   C2    = G1*s2*ci1          + G2*(s1-l3);
  
   // flux assembly
   F[ir ]    = NN*(0.5*(FL[ir ]+FR[ir ])-0.5*(l3*du[ir ] + C1   ));
   F[iru]    = NN*(0.5*(FL[iru]+FR[iru])-0.5*(l3*du[iru] + C1*ui + C2*n1[0]));
   F[irv]    = NN*(0.5*(FL[irv]+FR[irv])-0.5*(l3*du[irv] + C1*vi + C2*n1[1]));
   F[irE]    = NN*(0.5*(FL[irE]+FR[irE])-0.5*(l3*du[irE] + C1*hi + C2*ucp  ));
   
   
   /*** Take care of passive scalars ***/
   //xf_PassiveEulerScalars(PIS, &nPass, iPass);
   
   for (iP=4; iP<sr; iP++){
      
      // state rank index of passive scalar
      //iP = iPass[i];
      
      // left and right scalars
      kl  = UL[iP]*rhol1;
      kr  = UR[iP]*rhor1;
      drk = UR[iP] - UL[iP];
      
      // Roe-averaged scalar
      ki     = (di*kr+kl)*d1;
      
      // flux assembly
      F[iP] = NN*( 0.5*(UL[iP]*unl+UR[iP]*unr)-0.5*(l3*drk + C1*ki));
      
   } // i
   
   // account for mesh motion
   //if (vg != NULL) xf_AddMotionToFlux0(PIS, NN*0.5*vgn, UL, UR, F);
   
   return xf_OK;
   
}

/******************************************************************/
// Roe Riemann flux for 3d setup
static int
ConvFluxFace3D_Roe(const int sr, const real *UL, const real *UR,
                  const real *n, const real *vg, const real *RParam, real *F, 
                  const real Gamma, real *MaxCharSpeed)
{
   /*
    PURPOSE:
    
    Evaluates 3d Roe-averaged approximate Riemann flux.
    
    INPUTS:
    
    PIS: Position In State vector;: maps enumerated type to state index
    UL: left state vector
    UR: left state vector
    n: normal vector (not necessarily of unit magnitude)
    vg : interface velocity
    RParam: function real parameters
    
    OUTPUTS:
    
    F : flux dotted with n
    
    RETURN:
    
    Error Code
    */
   
   real rhol1, rhor1;
   real ul, vl, wl;
   real ur, vr, wr;
   real u2l;
   real u2r;
   real unl;
   real unr;
   real pl;
   real hl;
   real pr;
   real hr;
   real FL[MAXSR];
   real FR[MAXSR];
   real al, ar;

   real di;
   real d1;
   real ui;
   real vi;
   real wi;
   real hi;
   real af;
   real ucp;
   real c2;
   real ci;
   
   real du[MAXSR], ci1, c21, el[3], ep;
   
   real lam[3];
   
   real eps;
   real eps1;
   
   real s1;
   real s2;
   real l3;
   real G1;
   real G2;
   real ct1;
   real ct2;
   real C1;
   real C2;
   
   real Mar, Mal, Ma;
   enum xfe_Bool LM_Fix_flag;
   //!!!!!!!!!!
   LM_Fix_flag = xfe_True;
   
   real gmi, gmi1;
   real NN, NN1, n1[3];
   
   int i, ir, iru, irv, irw, irE;
   int iP, nPass, iPass[MAXSR];
   
   real kl, kr, drk, ki, vgn;
   
   //variable index
   //ir  = PIS[xfe_Density];
   //iru = PIS[xfe_XMomentum];
   //irv = PIS[xfe_YMomentum];
   //irw = PIS[xfe_ZMomentum];
   //irE = PIS[xfe_Energy];
   ir = 0; iru = 1; irv = 2; irw = 3; irE = 4;
   
   
   gmi = Gamma - 1.0;
   gmi1 = 1/gmi;
   
   NN = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   NN1 = 1/NN;
   
   for (i=0; i<3; i++) n1[i] = n[i]*NN1;
   
   if((UL[ir] <= 0)||(UR[ir]<0)) return xf_Error(xf_NON_PHYSICAL);
   
   // interface normal velocity
   //vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1]+vg[2]*n1[2];
   vgn = 0.0;

   // Left State
   rhol1 = 1/UL[ir];
   
   ul    = UL[iru]*rhol1;
   vl    = UL[irv]*rhol1;
   wl    = UL[irw]*rhol1;
   u2l   = (ul*ul + vl*vl + wl*wl)*UL[ir];
   unl   = (ul  *n1[0] + vl  *n1[1] + wl  *n1[2]);
   pl    = (UL[irE] - 0.5*u2l  )*gmi;
   hl    = (UL[irE] + pl  )*rhol1;
   al    = sqrt(Gamma*pl/UL[ir]);
   Mal   =  (fabs(ul)+fabs(vl)+fabs(wl)) / al;

   FL[ir]     = UL[ir]*unl;
   FL[iru]    = n1[0]*pl   + UL[iru]*unl;
   FL[irv]    = n1[1]*pl   + UL[irv]*unl;
   FL[irw]    = n1[2]*pl   + UL[irw]*unl;
   FL[irE]    = (pl   + UL[irE])*unl;
   
   // Right State
   rhor1 = 1/UR[ir];
   
   ur    = UR[iru]*rhor1;
   vr    = UR[irv]*rhor1;
   wr    = UR[irw]*rhor1;
   u2r   = (ur*ur + vr*vr + wr*wr)*UR[ir];
   unr   = (ur  *n1[0] + vr  *n1[1] + wr  *n1[2]);
   pr    = (UR[irE] - 0.5*u2r  )*gmi;
   hr    = (UR[irE] + pr  )*rhor1;
   ar    = sqrt(Gamma*pr/UR[ir]);
   Mar   = (fabs(ur)+fabs(vr)+fabs(wr)) / ar;
   
   FR[ir ]    = UR[ir]*unr;
   FR[iru]    = n1[0]*pr   + UR[iru]*unr;
   FR[irv]    = n1[1]*pr   + UR[irv]*unr;
   FR[irw]    = n1[2]*pr   + UR[irw]*unr;
   FR[irE]    = (pr   + UR[irE])*unr;
   
   (*MaxCharSpeed) = sqrt(ul*ul + vl*vl + wl*wl) + al;
   di = sqrt(ur*ur + vr*vr + wr*wr) + ar;
   if((*MaxCharSpeed) < di) 
      (*MaxCharSpeed) = di; 
   
   // Average state
   di     = sqrt(UR[ir]*rhol1);
   d1     = 1.0/(1.0+di);
   
   ui     = (di   *ur+   ul  )*d1;
   vi     = (di   *vr+   vl  )*d1;
   wi     = (di   *wr+   wl  )*d1;
   hi     = (di*hr+hl)*d1;
   
   af     = 0.5*(ui*ui   +vi*vi   +wi*wi  );
   ucp    = ui   *n1[0]+vi   *n1[1]+wi   *n1[2];
   c2     = gmi*(hi   -af   );
   
   if(c2 <= 0) return xf_Error(xf_NON_PHYSICAL);
   
   ci    = sqrt(c2);
   ci1   = 1/ci;
   
   // du = UR-UL
   du[ir ] = UR[ir ] - UL[ir ];
   du[iru] = UR[iru] - UL[iru];
   du[irv] = UR[irv] - UL[irv];
   du[irw] = UR[irw] - UL[irw];
   du[irE] = UR[irE] - UL[irE];
   
   // eigenvalues
   lam[0] = ucp-vgn+ ci;
   lam[1] = ucp-vgn- ci;
   lam[2] = ucp-vgn;
   
   // Entropy fix
   ep     = 1e-2;
   eps    = ep*ci;
   
   for (i=0; i<3; i++)
      if ((lam[i]<eps)&&(lam[i]>-eps)){
         eps1 = 1/eps;
         
         lam[i] = 0.5*(eps+lam[i]*lam[i]*eps1);
      }
   
   for (i=0; i<3; i++)
      if (lam[i]<0)
         el[i] = -1;
      else
         el[i] =  1;
   
   // The following parameters are described in the theory guide
   s1    = 0.5*(el[0]*lam   [0]+el[1]*lam   [1]);
   s2    = 0.5*(el[0]*lam   [0]-el[1]*lam   [1]);
   
   // third eigenvalue, absolute value
   l3     = el[2]*lam    [2];
   
   G1    = gmi*(af*du[ir] - ui*du[iru] - vi*du[irv] - wi*du[irw]+du[irE]);
   G2    = -ucp*du[ir]+du[iru]*n1[0]+du[irv]*n1[1]+du[irw]*n1[2];
   
   if(LM_Fix_flag){
      Ma = max(Mal, Mar);
      if(Ma>=1.0)
         Ma = 1.0;
      G2 *= Ma;
   }
   
   c21   = ci1*ci1;
   ct1   = (s2   *G1            *ci1);
   ct2   = (s2   *G2            *ci1);
   C1    = ct2   +(s1 - l3   )*G1*c21;
   C2    = ct1   +(s1 - l3   )*G2;
  
   //l3 = 0.; C1 = 0.; C2 = 0.;
   F[ir ]    = NN*(0.5*(FL[ir ]+FR[ir ])-0.5*(l3*du[ir ]+C1   ));
   F[iru]    = NN*(0.5*(FL[iru]+FR[iru])-0.5*(l3*du[iru]+C1   *ui+C2   *n1[0]));
   F[irv]    = NN*(0.5*(FL[irv]+FR[irv])-0.5*(l3*du[irv]+C1   *vi+C2   *n1[1]));
   F[irw]    = NN*(0.5*(FL[irw]+FR[irw])-0.5*(l3*du[irw]+C1   *wi+C2   *n1[2]));
   F[irE]    = NN*(0.5*(FL[irE]+FR[irE])-0.5*(l3*du[irE]+C1   *hi+C2   *ucp  ));
   
   
   /*** Take care of passive scalars ***/
   //xf_PassiveEulerScalars(PIS, &nPass, iPass);
   
   for (iP=5; iP<sr; iP++){
      
      // state rank index of passive scalar
      //iP = iPass[i];
      
      // left and right scalars
      kl  = UL[iP]*rhol1;
      kr  = UR[iP]*rhor1;
      drk = UR[iP] - UL[iP];
      
      // Roe-averaged scalar
      ki     = (di*kr+kl)*d1;
      
      // flux assembly
      F[iP] = NN*( 0.5*(UL[iP]*unl+UR[iP]*unr)-0.5*(l3*drk + C1*ki));
      
   } // i
   
   // account for mesh motion
   //if (vg != NULL) xf_AddMotionToFlux0(PIS, NN*0.5*vgn, UL, UR, F);
   
   return xf_OK;
   
}

/******************************************************************/
// Central flux for 2d setup; wrong for gas-dynamics but might work
// well for turbulence
static int
ConvFluxFace2D_Central(const int sr, const real *UL, const real *UR,
                       const real *n, real *F,
                       const real Gamma, real *MaxCharSpeed)
{
    /*
     PURPOSE:
     
     Evaluates 2d Rusanov-averaged approximate Riemann flux, and derivatives.
     
     INPUTS:
     
     PIS: Position In State vector;: maps enumerated type to state index
     sr : state rank
     UL: left state vector
     UR: left state vector
     n: normal vector (not necessarily of unit magnitude)
     vg : interface velocity
     RParam: function real parameters
     
     OUTPUTS:
     
     F : flux dotted with n
     F_UL, F_UR : flux derivatives w.r.t UL and UR
     
     RETURN:
     
     Error Code
     */
    
    real rhol1;
    real ul, vl;
    real rhor1;
    real ur, vr;
    real u2l, u2r, unl, unr;
    real pl, hl, pr, hr;
    real FL[MAXSR];
    real FR[MAXSR];
    
    real di, d1, ui, vi, hi, af;
    real ucp, c2, ci1, c, ci;
    real du[MAXSR];
    
    real lam[3], wl, wr;
    
    real gmi, gmi1;
    real NN, NN1, n1[2];
    
    int ic;
    int i, ierr;
    int ind, ir, iru, irv, irE;
    int j, iP, nPass, iPass[MAXSR];
    
    real drk, vgn;
    
    //variable index
    //ir  = PIS[xfe_Density];
    //iru = PIS[xfe_XMomentum];
    //irv = PIS[xfe_YMomentum];
    //irE = PIS[xfe_Energy];
    ir = 0; iru = 1; irv = 2; irE = 3;
    
    gmi = Gamma - 1.0;
    gmi1 = 1/gmi;
    
    NN = sqrt(n[0]*n[0] + n[1]*n[1]);
    NN1 = 1/NN;
    
    for (i=0; i<2; i++) n1[i] = n[i]*NN1;
    
    if((UL[ir] <= 0)||(UR[ir]<0)) return xf_Error(xf_NON_PHYSICAL);
    
    // interface normal velocity
    //vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1];
    vgn = 0.0;
    
    // Left State
    rhol1 = 1/UL[ir];
    ul    = UL[iru]*rhol1;
    vl    = UL[irv]*rhol1;
    u2l   = (ul*ul + vl*vl)*UL[ir];
    unl   = (ul   *n1[0] + vl   *n1[1]);
    pl    = (UL[irE] - 0.5*u2l   )*gmi;
    wl    = sqrt(ul*ul + vl*vl) + sqrt(pl*Gamma*rhol1);
    
    if (pl <= 0.0) return xf_Error(xf_NON_PHYSICAL);
    
    hl    = (UL[irE] + pl  )*rhol1;
    FL   [ir]      = UL[ir]*unl;
    FL   [iru]      = n1[0]*pl    + UL[iru]*unl;
    FL   [irv]      = n1[1]*pl    + UL[irv]*unl;
    FL   [irE]      = (pl   + UL[irE])*unl;
    
    // Right State
    rhor1 = 1/UR[ir];
    ur     = UR[iru]*rhor1;
    vr     = UR[irv]*rhor1;
    u2r    = (ur*ur + vr*vr)*UR[ir];
    unr    = (ur   *n1[0] + vr   *n1[1]);
    pr     = (UR[irE] - 0.5*u2r   )*gmi;
    wr     = sqrt(ur*ur + vr*vr) + sqrt(pr*Gamma*rhor1);
    
    if (pr <= 0.0){
        xf_printf("UR = %.10E %.10E %.10E %.10E\n", UR[ir], UR[iru], UR[irv], UR[irE]);
        return xf_Error(xf_NON_PHYSICAL);
    }
    
    hr     = (UR[irE] + pr   )*rhor1;
    FR   [ir]      = UR[ir]*unr;
    FR   [iru]      = n1[0]*pr    + UR[iru]*unr;
    FR   [irv]      = n1[1]*pr    + UR[irv]*unr;
    FR   [irE]      = (pr   + UR[irE])*unr;
    
    F   [ir       ] = NN*0.5*(FL[ir]+FR[ir ]);
    F   [iru      ] = NN*0.5*(FL[iru]+FR[iru]);
    F   [irv      ] = NN*0.5*(FL[irv]+FR[irv]);
    F   [irE      ] = NN*0.5*(FL[irE]+FR[irE]);
    
    /*** Take care of passive scalars ***/
    //xf_PassiveEulerScalars(PIS, &nPass, iPass);
    
    for (iP=4; iP<sr; iP++){
        
        // flux assembly for passive scalar
        F[iP] = NN*0.5*(UL[iP]*unl+UR[iP]*unr);
        
    } // i
    
    // account for mesh motion
    //if (vg != NULL) xf_AddMotionToFlux(PIS, NN*0.5*vgn, sr, UL, UR, F, F_UL, F_UR);
    
    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition:  xf_ConvFluxRusanov3d
static int
ConvFluxFace3D_Central(const int sr, const real *UL, const real *UR,
                       const real *n, real *F,
                       const real Gamma, real *MaxCharSpeed)
{
    /*
     PURPOSE:
     
     Evaluates 3d Rusanov-averaged approximate Riemann flux, and derivatives.
     
     INPUTS:
     
     PIS: Position In State vector;: maps enumerated type to state index
     sr : state rank
     UL: left state vector
     UR: left state vector
     n: normal vector (not necessarily of unit magnitude)
     vg : interface velocity
     RParam: function real parameters
     
     OUTPUTS:
     
     F : flux dotted with n
     F_UL, F_UR : flux derivatives w.r.t UL and UR
     
     RETURN:
     
     Error Code
     */
    
    
    real rhol1;
    real ul, vl, wl;
    real rhor1;
    real ur, vr, wr;
    real u2l, u2r, unl, unr;
    real pl, hl, pr, hr;
    real FL[MAXSR];
    real FR[MAXSR];
    real al, ar;
    
    real di, d1, ui, vi, wi, hi;
    real af, ucp;
    real c2, ci1, ci;
    
    real du[MAXSR];
    
    real c;
    
    real lam[3];
    
    real gmi, gmi1;
    real NN, NN1, n1[3];
    
    int i, ierr;
    int ic;
    int ind, ir, iru, irv, irw, irE;
    int j, iP, nPass, iPass[MAXSR];
    
    real drk, vgn;
    
    //variable index
    //ir  = PIS[xfe_Density];
    //iru = PIS[xfe_XMomentum];
    //irv = PIS[xfe_YMomentum];
    //irw = PIS[xfe_ZMomentum];
    //irE = PIS[xfe_Energy];
    ir = 0; iru = 1; irv = 2; irw = 3; irE = 4;
    
    gmi = Gamma - 1.0;
    gmi1 = 1/gmi;
    
    NN = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    NN1 = 1/NN;
    
    for (i=0; i<3; i++) n1[i] = n[i]*NN1;
    
    if((UL[ir] <= 0)||(UR[ir]<0)) return xf_Error(xf_NON_PHYSICAL);
    
    // interface normal velocity
    //vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1]+vg[2]*n1[2];
    
    // Left State
    rhol1 = 1/UL[ir];
    ul    = UL[iru]*rhol1;
    vl    = UL[irv]*rhol1;
    wl    = UL[irw]*rhol1;
    u2l   = (ul*ul + vl*vl + wl*wl)*UL[ir];
    
    unl   = (ul  *n1[0] + vl  *n1[1] + wl  *n1[2]);
    pl    = (UL[irE] - 0.5*u2l  )*gmi;
    hl    = (UL[irE] + pl  )*rhol1;
    al    = sqrt(Gamma * pl * rhol1);
    FL   [ir]    = UL[ir]*unl;
    FL   [iru]    = n1[0]*pl   + UL[iru]*unl;
    FL   [irv]    = n1[1]*pl   + UL[irv]*unl;
    FL   [irw]    = n1[2]*pl   + UL[irw]*unl;
    FL   [irE]    = (pl   + UL[irE])*unl;
    
    
    // Right State
    rhor1 = 1/UR[ir];
    ur    = UR[iru]*rhor1;
    vr    = UR[irv]*rhor1;
    wr    = UR[irw]*rhor1;
    u2r   = (ur*ur + vr*vr + wr*wr)*UR[ir];
    
    unr   = (ur  *n1[0] + vr  *n1[1] + wr  *n1[2]);
    pr    = (UR[irE] - 0.5*u2r  )*gmi;
    hr    = (UR[irE] + pr  )*rhor1;
    ar    = sqrt(Gamma * pr * rhor1);
    FR   [ir]    = UR[ir]*unr;
    FR   [iru]    = n1[0]*pr   + UR[iru]*unr;
    FR   [irv]    = n1[1]*pr   + UR[irv]*unr;
    FR   [irw]    = n1[2]*pr   + UR[irw]*unr;
    FR   [irE]    = (pr   + UR[irE])*unr;
    
    (*MaxCharSpeed) = sqrt(ul*ul + vl*vl + wl*wl) + al;
    if((*MaxCharSpeed) < (sqrt(ur*ur + vr*vr + wr*wr) + ar))
        (*MaxCharSpeed) = (sqrt(ur*ur + vr*vr + wr*wr) + ar);
    
    F   [ir       ] = NN*0.5*(FL[ir]+FR[ir]);
    F   [iru      ] = NN*0.5*(FL[iru]+FR[iru]);
    F   [irv      ] = NN*0.5*(FL[irv]+FR[irv]);
    F   [irw      ] = NN*0.5*(FL[irw]+FR[irw]);
    F   [irE      ] = NN*0.5*(FL[irE]+FR[irE]);
    
    
    /*** Take care of passive scalars ***/
    //xf_PassiveEulerScalars(PIS, &nPass, iPass);
    
    for (iP=5; iP<sr; iP++){
        
        // flux assembly for passive scalar
        F[iP] = NN*0.5*(UL[iP]*unl+UR[iP]*unr);
        
    } // i
    
    // account for mesh motion
    //if (vg != NULL) xf_AddMotionToFlux(PIS, NN*0.5*vgn, sr, UL, UR, F, F_UL, F_UR);
    
    return xf_OK;
    
}

/******************************************************************/
// Rusanov/local Lax-Friedrichs Riemann flux for 2d setup
static int
ConvFluxFace2D_Rusanov(const int sr, const real *UL, const real *UR,
                     const real *n, real *F, 
                     const real Gamma, real *MaxCharSpeed)
{
   /*
    PURPOSE:
    
    Evaluates 2d Rusanov-averaged approximate Riemann flux, and derivatives.
    
    INPUTS:
    
    PIS: Position In State vector;: maps enumerated type to state index
    sr : state rank
    UL: left state vector
    UR: left state vector
    n: normal vector (not necessarily of unit magnitude)
    vg : interface velocity
    RParam: function real parameters
    
    OUTPUTS:
    
    F : flux dotted with n
    F_UL, F_UR : flux derivatives w.r.t UL and UR
    
    RETURN:
    
    Error Code
    */
   
   real rhol1;
   real ul, vl;
   real rhor1;
   real ur, vr;
   real u2l, u2r, unl, unr;
   real pl, hl, pr, hr;
   real FL[MAXSR];
   real FR[MAXSR];
   
   real di, d1, ui, vi, hi, af;
   real ucp, c2, ci1, c, ci;
   real du[MAXSR];
   
   real lam[3], wl, wr;
   
   real gmi, gmi1;
   real NN, NN1, n1[2];
   
   int ic;
   int i, ierr;
   int ind, ir, iru, irv, irE;
   int j, iP, nPass, iPass[MAXSR];
   
   real drk, vgn;
   
   //variable index
   //ir  = PIS[xfe_Density];
   //iru = PIS[xfe_XMomentum];
   //irv = PIS[xfe_YMomentum];
   //irE = PIS[xfe_Energy];
   ir = 0; iru = 1; irv = 2; irE = 3;
   
   gmi = Gamma - 1.0;
   gmi1 = 1/gmi;
   
   NN = sqrt(n[0]*n[0] + n[1]*n[1]);
   NN1 = 1/NN;
   
   for (i=0; i<2; i++) n1[i] = n[i]*NN1;
   
   if((UL[ir] <= 0)||(UR[ir]<0)) return xf_Error(xf_NON_PHYSICAL);
   
   // interface normal velocity
   //vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1];
   vgn = 0.0;

   // Left State
   rhol1 = 1/UL[ir];
   ul    = UL[iru]*rhol1;
   vl    = UL[irv]*rhol1;
   u2l   = (ul*ul + vl*vl)*UL[ir];
   unl   = (ul   *n1[0] + vl   *n1[1]);
   pl    = (UL[irE] - 0.5*u2l   )*gmi;
   wl    = sqrt(ul*ul + vl*vl) + sqrt(pl*Gamma*rhol1);

   if (pl <= 0.0) return xf_Error(xf_NON_PHYSICAL);
   
   hl    = (UL[irE] + pl  )*rhol1;
   FL   [ir]      = UL[ir]*unl;
   FL   [iru]      = n1[0]*pl    + UL[iru]*unl;
   FL   [irv]      = n1[1]*pl    + UL[irv]*unl;
   FL   [irE]      = (pl   + UL[irE])*unl;
   
   // Right State
   rhor1 = 1/UR[ir];   
   ur     = UR[iru]*rhor1;
   vr     = UR[irv]*rhor1;
   u2r    = (ur*ur + vr*vr)*UR[ir];
   unr    = (ur   *n1[0] + vr   *n1[1]);
   pr     = (UR[irE] - 0.5*u2r   )*gmi;
   wr     = sqrt(ur*ur + vr*vr) + sqrt(pr*Gamma*rhor1);
   
   if (pr <= 0.0){
      xf_printf("UR = %.10E %.10E %.10E %.10E\n", UR[ir], UR[iru], UR[irv], UR[irE]);
      return xf_Error(xf_NON_PHYSICAL);
   }
   
   hr     = (UR[irE] + pr   )*rhor1;
   FR   [ir]      = UR[ir]*unr;
   FR   [iru]      = n1[0]*pr    + UR[iru]*unr;
   FR   [irv]      = n1[1]*pr    + UR[irv]*unr;
   FR   [irE]      = (pr   + UR[irE])*unr;
   
   
   // Average state
   di     = sqrt(UR[ir]*rhol1);
   d1     = 1.0/(1.0+di);
   ui     = (di   *ur+   ul  )*d1;
   vi     = (di   *vr+   vl  )*d1;
   hi      = (di*hr+hl)*d1;
   af      = 0.5*(ui*ui   +vi*vi);
   ucp     = ui    *n1[0]+vi    *n1[1];
   c2      = gmi*(hi   -af   );

   
   if(c2 <= 0)
      return xf_Error(xf_NON_PHYSICAL);
   
   ci     = sqrt(c2);
   ci1    = 1/ci;

   
   // du = UR-UL
   du[ir ] = UR[ir ] - UL[ir ];
   du[iru] = UR[iru] - UL[iru];
   du[irv] = UR[irv] - UL[irv];
   du[irE] = UR[irE] - UL[irE];
   
   // eigenvalues
   //lam   [0]  = ucp-vgn+ ci;
   //lam   [1]  = ucp-vgn- ci;
   //lam   [2]  = ucp-vgn;

   
   //ic = (lam[0] > lam[1]) ?
   //((lam[0] > lam[2]) ? 0 : 2) :
   //((lam[1] > lam[2]) ? 1 : 2);
   
   //c     = lam[ic];
   c  = wl;
   if(c < wr) c = wr;

   //maximum wave speed for time step estimate
   *MaxCharSpeed = c;

   F   [ir       ] = NN*(0.5*(FL[ir]+FR[ir ])-0.5*c*du[ir]);
   F   [iru      ] = NN*(0.5*(FL[iru]+FR[iru])-0.5*c*du[iru]);
   F   [irv      ] = NN*(0.5*(FL[irv]+FR[irv])-0.5*c*du[irv]);
   F   [irE      ] = NN*(0.5*(FL[irE]+FR[irE])-0.5*c*du[irE]);

   /*** Take care of passive scalars ***/
   //xf_PassiveEulerScalars(PIS, &nPass, iPass);
   
   for (iP=4; iP<sr; iP++){
      
      // state rank index of passive scalar
      //iP = iPass[i];
      
      drk = UR[iP] - UL[iP];
     
      // flux assembly
      F[iP] = NN*( 0.5*(UL[iP]*unl+UR[iP]*unr)-0.5*c*drk);
     
   } // i
   
   // account for mesh motion
   //if (vg != NULL) xf_AddMotionToFlux(PIS, NN*0.5*vgn, sr, UL, UR, F, F_UL, F_UR);
   
   return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition:  xf_ConvFluxRusanov3d
static int
ConvFluxFace3D_Rusanov(const int sr, const real *UL, const real *UR,
                     const real *n, real *F,
                     const real Gamma, real *MaxCharSpeed)
{
   /*
    PURPOSE:
    
    Evaluates 3d Rusanov-averaged approximate Riemann flux, and derivatives.
    
    INPUTS:
    
    PIS: Position In State vector;: maps enumerated type to state index
    sr : state rank
    UL: left state vector
    UR: left state vector
    n: normal vector (not necessarily of unit magnitude)
    vg : interface velocity
    RParam: function real parameters
    
    OUTPUTS:
    
    F : flux dotted with n
    F_UL, F_UR : flux derivatives w.r.t UL and UR
    
    RETURN:
    
    Error Code
    */
   
   
   real rhol1;
   real ul, vl, wl;
   real rhor1;
   real ur, vr, wr;
   real u2l, u2r, unl, unr;
   real pl, hl, pr, hr;
   real FL[MAXSR];
   real FR[MAXSR];
   real al, ar;

   real di, d1, ui, vi, wi, hi;
   real af, ucp;
   real c2, ci1, ci;
   
   real du[MAXSR];
   
   real c;
   
   real lam[3];
   
   real gmi, gmi1;
   real NN, NN1, n1[3];
   
   int i, ierr;
   int ic;
   int ind, ir, iru, irv, irw, irE;
   int j, iP, nPass, iPass[MAXSR];
   
   real drk, vgn;
   
   //variable index
   //ir  = PIS[xfe_Density];
   //iru = PIS[xfe_XMomentum];
   //irv = PIS[xfe_YMomentum];
   //irw = PIS[xfe_ZMomentum];
   //irE = PIS[xfe_Energy];
   ir = 0; iru = 1; irv = 2; irw = 3; irE = 4;
   
   gmi = Gamma - 1.0;
   gmi1 = 1/gmi;
   
   NN = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   NN1 = 1/NN;
   
   for (i=0; i<3; i++) n1[i] = n[i]*NN1;
   
   if((UL[ir] <= 0)||(UR[ir]<0)) return xf_Error(xf_NON_PHYSICAL);
   
   // interface normal velocity
   //vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1]+vg[2]*n1[2];
   
   // Left State
   rhol1 = 1/UL[ir];
   ul    = UL[iru]*rhol1;
   vl    = UL[irv]*rhol1;
   wl    = UL[irw]*rhol1;
   u2l   = (ul*ul + vl*vl + wl*wl)*UL[ir];
   
   unl   = (ul  *n1[0] + vl  *n1[1] + wl  *n1[2]);
   pl    = (UL[irE] - 0.5*u2l  )*gmi;
   hl    = (UL[irE] + pl  )*rhol1;
   al    = sqrt(Gamma * pl * rhol1);
   FL   [ir]    = UL[ir]*unl;
   FL   [iru]    = n1[0]*pl   + UL[iru]*unl;
   FL   [irv]    = n1[1]*pl   + UL[irv]*unl;
   FL   [irw]    = n1[2]*pl   + UL[irw]*unl;
   FL   [irE]    = (pl   + UL[irE])*unl;

   
   // Right State
   rhor1 = 1/UR[ir];
   ur    = UR[iru]*rhor1;
   vr    = UR[irv]*rhor1;
   wr    = UR[irw]*rhor1;
   u2r   = (ur*ur + vr*vr + wr*wr)*UR[ir];

   unr   = (ur  *n1[0] + vr  *n1[1] + wr  *n1[2]);
   pr    = (UR[irE] - 0.5*u2r  )*gmi;
   hr    = (UR[irE] + pr  )*rhor1;
   ar    = sqrt(Gamma * pr * rhor1);
   FR   [ir]    = UR[ir]*unr;
   FR   [iru]    = n1[0]*pr   + UR[iru]*unr;
   FR   [irv]    = n1[1]*pr   + UR[irv]*unr;
   FR   [irw]    = n1[2]*pr   + UR[irw]*unr;
   FR   [irE]    = (pr   + UR[irE])*unr;
     
   (*MaxCharSpeed) = sqrt(ul*ul + vl*vl + wl*wl) + al;
   if((*MaxCharSpeed) < (sqrt(ur*ur + vr*vr + wr*wr) + ar)) 
      (*MaxCharSpeed) = (sqrt(ur*ur + vr*vr + wr*wr) + ar); 
 /*  
   // Average state
   di     = sqrt(UR[ir]*rhol1);
   d1     = 1.0/(1.0+di);
   ui     = (di   *ur+   ul  )*d1;
   vi     = (di   *vr+   vl  )*d1;
   wi     = (di   *wr+   wl  )*d1;
   hi     = (di*hr+hl)*d1;
   
   af     = 0.5*(ui*ui   +vi*vi   +wi*wi  );
   ucp    = ui   *n1[0]+vi   *n1[1]+wi   *n1[2];
   c2     = gmi*(hi   -af   );
   
   if(c2 <= 0)
      return xf_Error(xf_NON_PHYSICAL);
   
   ci    = sqrt(c2);
   ci1   = 1/ci;
*/   
   // du = UR-UL
   du[ir ] = UR[ir ] - UL[ir ];
   du[iru] = UR[iru] - UL[iru];
   du[irv] = UR[irv] - UL[irv];
   du[irw] = UR[irw] - UL[irw];
   du[irE] = UR[irE] - UL[irE];
/*   
   // eigenvalues
   lam   [0] = ucp-vgn+ ci;
   lam   [1] = ucp-vgn- ci;
   lam   [2] = ucp-vgn ;
   
   ic = (lam[0] > lam[1]) ?
   ((lam[0] > lam[2]) ? 0 : 2) :
   ((lam[1] > lam[2]) ? 1 : 2);
   
   c     = lam[ic];
 */  
   c     = (*MaxCharSpeed);
   
   F   [ir       ] = NN*(0.5*(FL[ir]+FR[ir])-0.5*c*du[ir]);
   F   [iru      ] = NN*(0.5*(FL[iru]+FR[iru])-0.5*c*du[iru]);
   F   [irv      ] = NN*(0.5*(FL[irv]+FR[irv])-0.5*c*du[irv]);
   F   [irw      ] = NN*(0.5*(FL[irw]+FR[irw])-0.5*c*du[irw]);
   F   [irE      ] = NN*(0.5*(FL[irE]+FR[irE])-0.5*c*du[irE]);

   
   /*** Take care of passive scalars ***/
   //xf_PassiveEulerScalars(PIS, &nPass, iPass);
   
   for (iP=5; iP<sr; iP++){
      
      // state rank index of passive scalar
      //iP = iPass[i];
      
      drk = UR[iP] - UL[iP];
      
      // flux assembly
      F[iP] = NN*( 0.5*(UL[iP]*unl+UR[iP]*unr)-0.5*c*drk);
      
   } // i
   
   // account for mesh motion
   //if (vg != NULL) xf_AddMotionToFlux(PIS, NN*0.5*vgn, sr, UL, UR, F, F_UL, F_UR);
   
   return xf_OK;
   
}

/******************************************************************/
// Entrance for Riemann Flux Evaluation
int
ConvFluxInteriorFace(const int nq, const int sr, const int dim, const real *qUL, const real *qUR,
                     const real *qn, real *qF, const real Gamma, real *MaxCharSpeed)
{
   int ierr, i, sr2;
   const real *UL, *UR, *n;
   real *F;
   
   sr2 = sr * sr;
   
   for (i=0; i<nq; i++){
      
      UL   = qUL+i*sr;
      UR   = qUR+i*sr;
      n    = qn +i*dim;
      F    = qF +i*sr;
      
      
      if(dim == 2){
         //ierr = xf_Error(ConvFluxFace2D(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
         //ierr = xf_Error(ConvFluxFace2D_Rusanov(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
         ierr = xf_Error(ConvFluxFace2D_Roe(sr, UL, UR, n, NULL, NULL, F, Gamma, MaxCharSpeed));
         if (ierr != xf_OK) return ierr;
      }
      else if(dim == 3){
         //ierr = xf_Error(ConvFluxFace3D(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
         //if (ierr != xf_OK) return ierr;
      //   ierr = xf_Error(ConvFluxFace3D_Rusanov(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
      //   if (ierr != xf_OK) return ierr;
      //    ierr = xf_Error(ConvFluxFace3D_Roe(sr, UL, UR, n, NULL, NULL, F, Gamma, MaxCharSpeed));
      //    if (ierr != xf_OK) return ierr;
         ierr = xf_Error(ConvFluxFace3D_Central(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
         if (ierr != xf_OK) return ierr;
      }
      else
         return xf_Error(xf_OUT_OF_BOUNDS);
      
   }
   
   return xf_OK;
}

/******************************************************************/
// Entrance for Riemann Flux Evaluation
int
ConvFluxBoundaryFace(const int nq, const int sr, const int dim, const real *qUL, const real *qUR,
                     const real *qn, real *qF, const real Gamma, real *MaxCharSpeed)
{
    int ierr, i, sr2;
    const real *UL, *UR, *n;
    real *F;
    
    sr2 = sr * sr;
    
    for (i=0; i<nq; i++){
        
        UL   = qUL+i*sr;
        UR   = qUR+i*sr;
        n    = qn +i*dim;
        F    = qF +i*sr;
        
        
        if(dim == 2){
            //ierr = xf_Error(ConvFluxFace2D(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            ierr = xf_Error(ConvFluxFace2D_Rusanov(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            //ierr = xf_Error(ConvFluxFace2D_Roe(sr, UL, UR, n, NULL, NULL, F, Gamma, MaxCharSpeed));
            if (ierr != xf_OK) return ierr;
        }
        else if(dim == 3){
            //ierr = xf_Error(ConvFluxFace3D(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            //if (ierr != xf_OK) return ierr;
            //   ierr = xf_Error(ConvFluxFace3D_Rusanov(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            //   if (ierr != xf_OK) return ierr;
            ierr = xf_Error(ConvFluxFace3D_Roe(sr, UL, UR, n, NULL, NULL, F, Gamma, MaxCharSpeed));
            if (ierr != xf_OK) return ierr;
            //ierr = xf_Error(ConvFluxFace3D_Central(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            //if (ierr != xf_OK) return ierr;
        }
        else
            return xf_Error(xf_OUT_OF_BOUNDS);
        
    }
    
    return xf_OK;
}

/******************************************************************/
// Roe Riemann flux for 3d setup
static int
ConvFluxFace3D_Dissip_Roe(const int sr, const real *UL, const real *UR,
                          const real *n, const real Gamma, real *Qn, 
                          const real eta)
{
    /*
     PURPOSE:
     
     Evaluates 3d Roe-averaged approximate Riemann flux.
     
     INPUTS:
     
     PIS: Position In State vector;: maps enumerated type to state index
     UL: left state vector
     UR: left state vector
     n: normal vector (not necessarily of unit magnitude)
     vg : interface velocity
     RParam: function real parameters
     
     OUTPUTS:
     
     F : flux dotted with n
     
     RETURN:
     
     Error Code
     */
    
    real rhol1, rhor1;
    real ul, vl, wl;
    real ur, vr, wr;
    real u2l;
    real u2r;
    real unl;
    real unr;
    real pl;
    real hl;
    real pr;
    real hr;
//    real FL[MAXSR];
//    real FR[MAXSR];
    real al, ar;
    
    real di;
    real d1;
    real ui;
    real vi;
    real wi;
    real hi;
    real af;
    real ucp;
    real c2;
    real ci;
    
    real du[MAXSR], ci1, c21, el[3], ep;
    
    real lam[3];
    
    real eps;
    real eps1;
    
    real s1;
    real s2;
    real l3;
    real G1;
    real G2;
    real ct1;
    real ct2;
    real C1;
    real C2;
    
    real Mar, Mal, Ma;
    enum xfe_Bool LM_Fix_flag;
    //!!!!!!!!!!
    LM_Fix_flag = xfe_True;
    
    real gmi, gmi1;
    real NN, NN1, n1[3];
    
    int i, ir, iru, irv, irw, irE;
    int iP, nPass, iPass[MAXSR];
    
    real kl, kr, drk, ki, vgn;
    
    //variable index
    //ir  = PIS[xfe_Density];
    //iru = PIS[xfe_XMomentum];
    //irv = PIS[xfe_YMomentum];
    //irw = PIS[xfe_ZMomentum];
    //irE = PIS[xfe_Energy];
    ir = 0; iru = 1; irv = 2; irw = 3; irE = 4;
    
    
    gmi = Gamma - 1.0;
    gmi1 = 1/gmi;
    
    NN = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    NN1 = 1/NN;
    
    for (i=0; i<3; i++) n1[i] = n[i]*NN1;
    
    if((UL[ir] <= 0)||(UR[ir]<0)) return xf_Error(xf_NON_PHYSICAL);
    if((UL[irE] <= 0)||(UR[irE]<0)) return xf_Error(xf_NON_PHYSICAL);
    
    // interface normal velocity
    //vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1]+vg[2]*n1[2];
    vgn = 0.0;
    
    // Left State
    rhol1 = 1/UL[ir];
    
    ul    = UL[iru]*rhol1;
    vl    = UL[irv]*rhol1;
    wl    = UL[irw]*rhol1;
    u2l   = (ul*ul + vl*vl + wl*wl)*UL[ir];
    unl   = (ul  *n1[0] + vl  *n1[1] + wl  *n1[2]);
    pl    = (UL[irE] - 0.5*u2l  )*gmi;
    hl    = (UL[irE] + pl  )*rhol1;
    al    = sqrt(Gamma*pl/UL[ir]);
    Mal   = sqrt(ul*ul + vl*vl + wl*wl)/ al;
    
    // Right State
    rhor1 = 1/UR[ir];
    
    ur    = UR[iru]*rhor1;
    vr    = UR[irv]*rhor1;
    wr    = UR[irw]*rhor1;
    u2r   = (ur*ur + vr*vr + wr*wr)*UR[ir];
    unr   = (ur  *n1[0] + vr  *n1[1] + wr  *n1[2]);
    pr    = (UR[irE] - 0.5*u2r  )*gmi;
    hr    = (UR[irE] + pr  )*rhor1;
    ar    = sqrt(Gamma*pr/UR[ir]);
    Mar   = sqrt(ur*ur + vr*vr + wr*wr) / ar;
    
    // du = UR-UL
    du[ir ] = UR[ir ] - UL[ir ];
    du[iru] = UR[iru] - UL[iru];
    du[irv] = UR[irv] - UL[irv];
    du[irw] = UR[irw] - UL[irw];
    du[irE] = UR[irE] - UL[irE];
    
    // Average state
    di     = sqrt(UR[ir]*rhol1);
    d1     = 1.0/(1.0+di);
    
    ui     = (di   *ur+   ul  )*d1;
    vi     = (di   *vr+   vl  )*d1;
    wi     = (di   *wr+   wl  )*d1;
    hi     = (di*hr+hl)*d1;
    
    af     = 0.5*(ui*ui   +vi*vi   +wi*wi  );
    ucp    = ui   *n1[0]+vi   *n1[1]+wi   *n1[2];
    c2     = gmi*(hi   -af   );
    
    
    
    //if(c2 <= 0)
    if(xfe_True)
    {   //never give up return xf_Error(xf_NON_PHYSICAL);
        //this situation means Roe average does not exist
        vgn = max(Mal * al + al, Mar * ar + ar);
        
        //Opt 1: simple trick
        vgn = max(0.5*vgn, eta);      
        
        Qn   [ir       ] = NN*(vgn*du[ir]);
        Qn   [iru      ] = NN*(vgn*du[iru]);
        Qn   [irv      ] = NN*(vgn*du[irv]);
        Qn   [irw      ] = NN*(vgn*du[irw]);
        Qn   [irE      ] = NN*(vgn*du[irE]);
        
        
        /*** Take care of passive scalars ***/
        //xf_PassiveEulerScalars(PIS, &nPass, iPass);
        
        for (iP=5; iP<sr; iP++){
            
            // state rank index of passive scalar
            //iP = iPass[i];
            
            drk = UR[iP] - UL[iP];
            
            // flux assembly
            Qn[iP] = NN*(vgn*drk);
            
        } // i
        
    }
    else
    {
    ci    = sqrt(c2);
    ci1   = 1/ci;
    
    // eigenvalues
    lam[0] = ucp + ci;
    lam[1] = ucp - ci;
    lam[2] = ucp;
    
    // Entropy fix
 /*   ep     = 1e-2;
    eps    = ep*ci;
    
    for (i=0; i<3; i++)
        if ((lam[i]<eps)&&(lam[i]>-eps)){
            eps1 = 1/eps;
            
            lam[i] = 0.5*(eps+lam[i]*lam[i]*eps1);
        }
 */   
    for (i=0; i<3; i++)
        if (lam[i]<0)
            el[i] = -1;
        else
            el[i] =  1;
    
    // The following parameters are described in the theory guide
    s1    = 0.5*(el[0]*lam   [0]+el[1]*lam   [1]);
    s2    = 0.5*(el[0]*lam   [0]-el[1]*lam   [1]);
    
    // third eigenvalue, absolute value
    l3     = el[2]*lam    [2];
    
    G1    = gmi*(af*du[ir] - ui*du[iru] - vi*du[irv] - wi*du[irw]+du[irE]);
    G2    = -ucp*du[ir]+du[iru]*n1[0]+du[irv]*n1[1]+du[irw]*n1[2];
    
    if(LM_Fix_flag){
        Ma = max(Mal, Mar);
        if(Ma>=1.0)
            Ma = 1.0;
        G2 *= Ma;
    }
    
    c21   = ci1*ci1;
    ct1   = (s2   *G1            *ci1);
    ct2   = (s2   *G2            *ci1);
    C1    = ct2   +(s1 - l3   )*G1*c21;
    C2    = ct1   +(s1 - l3   )*G2;
    
/*    Qn[ir]   = max( Qn[ir],  NN*0.5*(l3*du[ir ]+C1) );
    Qn[iru]  = max( Qn[iru], NN*0.5*(l3*du[iru]+C1   *ui+C2   *n1[0]));
    Qn[irv]  = max( Qn[irv], NN*0.5*(l3*du[irv]+C1   *vi+C2   *n1[1]));
    Qn[irw]  = max( Qn[irw], NN*0.5*(l3*du[irw]+C1   *wi+C2   *n1[2]));
    Qn[irE]  = max( Qn[irE], NN*0.5*(l3*du[irE]+C1   *hi+C2   *ucp  ));
*/    
    Qn[ir]   = ( NN*0.5*(l3*du[ir ]+C1) );
    Qn[iru]  = ( NN*0.5*(l3*du[iru]+C1   *ui+C2   *n1[0]));
    Qn[irv]  = ( NN*0.5*(l3*du[irv]+C1   *vi+C2   *n1[1]));
    Qn[irw]  = ( NN*0.5*(l3*du[irw]+C1   *wi+C2   *n1[2]));
    Qn[irE]  = ( NN*0.5*(l3*du[irE]+C1   *hi+C2   *ucp  ));
    
    
    /*** Take care of passive scalars ***/
    //xf_PassiveEulerScalars(PIS, &nPass, iPass);
    
    for (iP=5; iP<sr; iP++){
        
        // state rank index of passive scalar
        //iP = iPass[i];
        
        // left and right scalars
        kl  = UL[iP]*rhol1;
        kr  = UR[iP]*rhor1;
        drk = UR[iP] - UL[iP];
        
        // Roe-averaged scalar
        ki     = (di*kr+kl)*d1;
        
        // flux assembly
        //Qn[iP] = max(Qn[iP], NN*0.5*(l3*drk + C1*ki));
        Qn[iP] = ( NN*0.5*(l3*drk + C1*ki));
        
    } // i
    
    }
    
    return xf_OK;
    
}

/******************************************************************/
// Roe Riemann flux for 2d setup
static int
ConvFluxFace2D_Adaptive_Roe(const int sr, const real *UL, const real *UR,
                          const real *n, const real Gamma, real *Qn,
                          const real eta)
{
    /*
     PURPOSE:
     
     Evaluates 2d Roe-averaged approximate Riemann flux.
     
     INPUTS:
     
     PIS: Position In State vector;: maps enumerated type to state index
     UL: left state vector
     UR: left state vector
     n: normal vector (not necessarily of unit magnitude)
     vg : interface velocity
     RParam: function real parameters
     
     OUTPUTS:
     
     F : flux dotted with n
     
     RETURN:
     
     Error Code
     */
    
    real rhol1, rhor1;
    real ul, vl, wl;
    real ur, vr, wr;
    real u2l;
    real u2r;
    real unl;
    real unr;
    real pl;
    real hl;
    real pr;
    real hr;
    //    real FL[MAXSR];
    //    real FR[MAXSR];
    real al, ar;
    
    real di;
    real d1;
    real ui;
    real vi;
    real wi;
    real hi;
    real af;
    real ucp;
    real c2;
    real ci;
    
    real du[MAXSR], ci1, c21, el[3], ep;
    
    real lam[3];
    
    real eps;
    real eps1;
    real tmp, max_char;
    
    real s1;
    real s2;
    real l3;
    real G1;
    real G2;
    real ct1;
    real ct2;
    real C1;
    real C2;
    
    real Mar, Mal, Ma;
    real alpha[MAXSR];
    enum xfe_Bool LM_Fix_flag;
    //!!!!!!!!!!
    LM_Fix_flag = xfe_True;
    
    real gmi, gmi1;
    real NN, NN1, n1[3];
    
    int i, ir, iru, irv, irw, irE;
    int iP, nPass, iPass[MAXSR];
    
    real kl, kr, drk, ki, vgn;
    real enpl, enpr;
    
    //variable index
    //ir  = PIS[xfe_Density];
    //iru = PIS[xfe_XMomentum];
    //irv = PIS[xfe_YMomentum];
    //irE = PIS[xfe_Energy];
    ir = 0; iru = 1; irv = 2; irE = 3;
    
    
    gmi = Gamma - 1.0;
    gmi1 = 1/gmi;
    
    NN = sqrt(n[0]*n[0] + n[1]*n[1]);
    NN1 = 1/NN;
    
    for (i=0; i<2; i++) n1[i] = n[i]*NN1;
    
    if((UL[ir] <= 0)||(UR[ir]<0)) return xf_Error(xf_NON_PHYSICAL);
    if((UL[irE] <= 0)||(UR[irE]<0)) return xf_Error(xf_NON_PHYSICAL);
    
    // interface normal velocity
    //vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1]+vg[2]*n1[2];
    vgn = 0.0;
    
    // Left State
    rhol1 = 1/UL[ir];
    
    ul    = UL[iru]*rhol1;
    vl    = UL[irv]*rhol1;
    u2l   = (ul*ul + vl*vl)*UL[ir];
    unl   = (ul  *n1[0] + vl  *n1[1]);
    pl    = (UL[irE] - 0.5*u2l  )*gmi;
    hl    = (UL[irE] + pl  )*rhol1;
    al    = sqrt(Gamma*pl/UL[ir]);
    Mal   = sqrt(ul*ul + vl*vl)/ al;
    enpl  = pl/pow(UL[ir], Gamma);
    
    // Right State
    rhor1 = 1/UR[ir];
    
    ur    = UR[iru]*rhor1;
    vr    = UR[irv]*rhor1;
    u2r   = (ur*ur + vr*vr)*UR[ir];
    unr   = (ur  *n1[0] + vr  *n1[1]);
    pr    = (UR[irE] - 0.5*u2r  )*gmi;
    hr    = (UR[irE] + pr  )*rhor1;
    ar    = sqrt(Gamma*pr/UR[ir]);
    Mar   = sqrt(ur*ur + vr*vr) / ar;
    enpr  = pr/pow(UR[ir], Gamma);
    
    // du = UR-UL
    du[ir ] = UR[ir ] - UL[ir ];
    du[iru] = UR[iru] - UL[iru];
    du[irv] = UR[irv] - UL[irv];
    du[irE] = UR[irE] - UL[irE];
    
    // Average state
    di     = sqrt(UR[ir]*rhol1);
    d1     = 1.0/(1.0+di);
    
    ui     = (di   *ur+   ul  )*d1;
    vi     = (di   *vr+   vl  )*d1;
    hi     = (di*hr+hl)*d1;
    
    af     = 0.5*(ui*ui   +vi*vi  );
    ucp    = ui   *n1[0]+vi   *n1[1];
    c2     = gmi*(hi   -af   );
    
    
    //no matter how to adapt the Riemann;
    //the stabilization for diffusion flux term is respected
    vgn = eta;
    
    Qn   [ir       ] += NN*(vgn*du[ir]);
    Qn   [iru      ] += NN*(vgn*du[iru]);
    Qn   [irv      ] += NN*(vgn*du[irv]);
    Qn   [irE      ] += NN*(vgn*du[irE]);
    for (iP=4; iP<sr; iP++){
        
        // state rank index of passive scalar
        //iP = iPass[i];
        
        return xf_NOT_SUPPORTED;
        drk = UR[iP] - UL[iP];
        
        // flux assembly
        Qn[iP] += NN*(vgn*drk);
        
    } // i
    
  /*
    ci    = sqrt(c2);
    ci1   = 1/ci;
    gmic2 = (Gamma-1.)*0.5*ci1*ci1;
    ucg   = ucp * ci * gmi1;
    
    G1    = gmi*(af*du[ir] - ui*du[iru] - vi*du[irv] - wi*du[irw]+du[irE]);
    G2    = -ucp*du[ir]+du[iru]*n1[0]+du[irv]*n1[1]+du[irw]*n1[2];
    
    //trick starts: step 1-project Qn onto char-space
    Qn_u =
    alpha[0] = gmic2 * ((af-ucg) * Qn[0] + ci*gmi1 * );
    */
    
    //if(c2 <= 0)
    if(xfe_False)
    {   //never give up return xf_Error(xf_NON_PHYSICAL);
        //this situation means Roe average does not exist
        vgn = 0.5*max(Mal * al + al, Mar * ar + ar);
 /*       if(enpr < 0.5 || enpl < 0.5)
        {
            Qn   [ir       ] += NN*(vgn*du[ir]);
            Qn   [iru      ] += NN*(vgn*du[iru]);
            Qn   [irv      ] += NN*(vgn*du[irv]);
            Qn   [irw      ] += NN*(vgn*du[irw]);
            Qn   [irE      ] += NN*(vgn*du[irE]);
        }
  */      
        
        
        
        //Opt 1: simple trick
        eps = 1.e-13;
        
        //test sign consistency
        max_char = 0.;
        for(iP=0; iP<sr; iP++)
        {
            if(fabs(Qn[iP]) < eps)
            {
                alpha[iP] = 0.;
                continue;
            }
            
            tmp = Qn[iP]/( (real)sign(du[iP]) * (fabs(du[iP]) + eps));
            
            if(tmp< 0.){
             //sorry, this idea does not work; just add dissiptation as usual
               vgn = 0.5*max(Mal * al + al, Mar * ar + ar);
                
                Qn   [ir       ] += NN*(vgn*du[ir]);
                Qn   [iru      ] += NN*(vgn*du[iru]);
                Qn   [irv      ] += NN*(vgn*du[irv]);
                Qn   [irE      ] += NN*(vgn*du[irE]);
                break;
            }
            else
            {
                alpha[iP] = tmp;
                if(tmp > max_char)
                    max_char = tmp;
            }
        }
        
        if(iP == sr){
            //this case means all signs are consistent
            vgn = 0.5*max(Mal * al + al, Mar * ar + ar);
            for(i=0; i<sr; i++)
                Qn[i] = NN * vgn *  du[i];
        }


  /*
        if(iP == sr){
            //this case means all signs are consistent
            vgn = 0.5*max(Mal * al + al, Mar * ar + ar);
            if(max_char < vgn){
                for(i=0; i<sr; i++)
                    Qn[i] = NN*(vgn * du[i]);
            
                xf_printf("I am good\n");
            }
            else
            {
                if(max_char > 2. *vgn)
                {
                    //give up
                    for(i=0; i<sr; i++)
                        Qn[i] += NN*(vgn * du[i]);
                }
                else
                {
                    //not give up
                    for(i=0; i<sr; i++)
                        Qn[i] = NN*(max_char * du[i]);
                }
            }
        }
    */
    
        
    }
    
    if(xfe_True)
    {
        ci    = sqrt(c2);
        ci1   = 1/ci;
        
        // eigenvalues
        //lam[0] = ucp + ci;
        //lam[1] = ucp - ci;

        if(LM_Fix_flag){
            Ma = max(Mal, Mar);
            if(Ma>=0.8)
                Ma = 1.0;

            //cut-off
            if(Ma<0.01)
                Ma = 0.01;
        }

        //Ma = 1.;
        lam[0] = ucp + ci * Ma;
        lam[1] = ucp - ci * Ma;
        lam[2] = ucp;
        
        // Entropy fix
         ep     = 1e-2;
         eps    = ep*ci;
         
         for (i=0; i<3; i++)
         if ((lam[i]<eps)&&(lam[i]>-eps)){
         eps1 = 1/eps;
         
         lam[i] = 0.5*(eps+lam[i]*lam[i]*eps1);
         }
         
        for (i=0; i<3; i++)
            if (lam[i]<0)
                el[i] = -1;
            else
                el[i] =  1;
        
        // The following parameters are described in the theory guide
        s1    = 0.5*(el[0]*lam   [0]+el[1]*lam   [1]);
        s2    = 0.5*(el[0]*lam   [0]-el[1]*lam   [1]);
        
        // third eigenvalue, absolute value
        l3     = el[2]*lam    [2];
        
        G1    = gmi*(af*du[ir] - ui*du[iru] - vi*du[irv] +du[irE]);
        G2    = -ucp*du[ir]+du[iru]*n1[0]+du[irv]*n1[1];
        
    /*    if(LM_Fix_flag){
            Ma = max(Mal, Mar);
            if(Ma>=1.0)
                Ma = 1.0;
            G2 *= Ma;
        }
    */    
        c21   = ci1*ci1;
        ct1   = (s2   *G1            *ci1);
        ct2   = (s2   *G2            *ci1);
        C1    = ct2   +(s1 - l3   )*G1*c21;
        C2    = ct1   +(s1 - l3   )*G2;
        
        /*    Qn[ir]   = max( Qn[ir],  NN*0.5*(l3*du[ir ]+C1) );
         Qn[iru]  = max( Qn[iru], NN*0.5*(l3*du[iru]+C1   *ui+C2   *n1[0]));
         Qn[irv]  = max( Qn[irv], NN*0.5*(l3*du[irv]+C1   *vi+C2   *n1[1]));
         Qn[irw]  = max( Qn[irw], NN*0.5*(l3*du[irw]+C1   *wi+C2   *n1[2]));
         Qn[irE]  = max( Qn[irE], NN*0.5*(l3*du[irE]+C1   *hi+C2   *ucp  ));
         */
        Qn[ir]   += ( NN*0.5*(l3*du[ir ]+C1) );
        Qn[iru]  += ( NN*0.5*(l3*du[iru]+C1   *ui+C2   *n1[0]));
        Qn[irv]  += ( NN*0.5*(l3*du[irv]+C1   *vi+C2   *n1[1]));
        Qn[irE]  += ( NN*0.5*(l3*du[irE]+C1   *hi+C2   *ucp  ));
        
        
        /*** Take care of passive scalars ***/
        //xf_PassiveEulerScalars(PIS, &nPass, iPass);
        
        for (iP=4; iP<sr; iP++){
            
            // state rank index of passive scalar
            //iP = iPass[i];
            
            // left and right scalars
            kl  = UL[iP]*rhol1;
            kr  = UR[iP]*rhor1;
            drk = UR[iP] - UL[iP];
            
            // Roe-averaged scalar
            ki     = (di*kr+kl)*d1;
            
            // flux assembly
            //Qn[iP] = max(Qn[iP], NN*0.5*(l3*drk + C1*ki));
            Qn[iP] += ( NN*0.5*(l3*drk + C1*ki));
            
        } // i
        
    }
    
    return xf_OK;
    
}

/******************************************************************/
// Roe Riemann flux for 3d setup
static int
ConvFluxFace3D_Adaptive_Roe(const int sr, const real *UL, const real *UR,
                          const real *n, const real Gamma, real *Qn,
                          const real eta)
{
    /*
     PURPOSE:
     
     Evaluates 3d Roe-averaged approximate Riemann flux.
     
     INPUTS:
     
     PIS: Position In State vector;: maps enumerated type to state index
     UL: left state vector
     UR: left state vector
     n: normal vector (not necessarily of unit magnitude)
     vg : interface velocity
     RParam: function real parameters
     
     OUTPUTS:
     
     F : flux dotted with n
     
     RETURN:
     
     Error Code
     */
    
    real rhol1, rhor1;
    real ul, vl, wl;
    real ur, vr, wr;
    real u2l;
    real u2r;
    real unl;
    real unr;
    real pl;
    real hl;
    real pr;
    real hr;
    //    real FL[MAXSR];
    //    real FR[MAXSR];
    real al, ar;
    
    real di;
    real d1;
    real ui;
    real vi;
    real wi;
    real hi;
    real af;
    real ucp;
    real c2;
    real ci;
    
    real du[MAXSR], ci1, c21, el[3], ep;
    
    real lam[3];
    
    real eps;
    real eps1;
    real tmp, max_char;
    
    real s1;
    real s2;
    real l3;
    real G1;
    real G2;
    real ct1;
    real ct2;
    real C1;
    real C2;
    
    real Mar, Mal, Ma;
    real alpha[MAXSR];
    enum xfe_Bool LM_Fix_flag;
    //!!!!!!!!!!
    LM_Fix_flag = xfe_True;
    
    real gmi, gmi1;
    real NN, NN1, n1[3];
    
    int i, ir, iru, irv, irw, irE;
    int iP, nPass, iPass[MAXSR];
    
    real kl, kr, drk, ki, vgn;
    real enpl, enpr;
    
    //variable index
    //ir  = PIS[xfe_Density];
    //iru = PIS[xfe_XMomentum];
    //irv = PIS[xfe_YMomentum];
    //irw = PIS[xfe_ZMomentum];
    //irE = PIS[xfe_Energy];
    ir = 0; iru = 1; irv = 2; irw = 3; irE = 4;
    
    
    gmi = Gamma - 1.0;
    gmi1 = 1/gmi;
    
    NN = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    NN1 = 1/NN;
    
    for (i=0; i<3; i++) n1[i] = n[i]*NN1;
    
    if((UL[ir] <= 0)||(UR[ir]<0)) return xf_Error(xf_NON_PHYSICAL);
    if((UL[irE] <= 0)||(UR[irE]<0)) return xf_Error(xf_NON_PHYSICAL);
    
    // interface normal velocity
    //vgn = (vg == NULL) ? 0. : vg[0]*n1[0]+vg[1]*n1[1]+vg[2]*n1[2];
    vgn = 0.0;
    
    // Left State
    rhol1 = 1/UL[ir];
    
    ul    = UL[iru]*rhol1;
    vl    = UL[irv]*rhol1;
    wl    = UL[irw]*rhol1;
    u2l   = (ul*ul + vl*vl + wl*wl)*UL[ir];
    unl   = (ul  *n1[0] + vl  *n1[1] + wl  *n1[2]);
    pl    = (UL[irE] - 0.5*u2l  )*gmi;
    hl    = (UL[irE] + pl  )*rhol1;
    al    = sqrt(Gamma*pl/UL[ir]);
    Mal   = sqrt(ul*ul + vl*vl + wl*wl)/ al;
    enpl  = pl/pow(UL[ir], Gamma);
    
    // Right State
    rhor1 = 1/UR[ir];
    
    ur    = UR[iru]*rhor1;
    vr    = UR[irv]*rhor1;
    wr    = UR[irw]*rhor1;
    u2r   = (ur*ur + vr*vr + wr*wr)*UR[ir];
    unr   = (ur  *n1[0] + vr  *n1[1] + wr  *n1[2]);
    pr    = (UR[irE] - 0.5*u2r  )*gmi;
    hr    = (UR[irE] + pr  )*rhor1;
    ar    = sqrt(Gamma*pr/UR[ir]);
    Mar   = sqrt(ur*ur + vr*vr + wr*wr) / ar;
    enpr  = pr/pow(UR[ir], Gamma);
    
    // du = UR-UL
    du[ir ] = UR[ir ] - UL[ir ];
    du[iru] = UR[iru] - UL[iru];
    du[irv] = UR[irv] - UL[irv];
    du[irw] = UR[irw] - UL[irw];
    du[irE] = UR[irE] - UL[irE];
    
    // Average state
    di     = sqrt(UR[ir]*rhol1);
    d1     = 1.0/(1.0+di);
    
    ui     = (di   *ur+   ul  )*d1;
    vi     = (di   *vr+   vl  )*d1;
    wi     = (di   *wr+   wl  )*d1;
    hi     = (di*hr+hl)*d1;
    
    af     = 0.5*(ui*ui   +vi*vi   +wi*wi  );
    ucp    = ui   *n1[0]+vi   *n1[1]+wi   *n1[2];
    c2     = gmi*(hi   -af   );
    
    
    //no matter how to adapt the Riemann;
    //the stabilization for diffusion flux term is respected
    vgn = eta;
    
    Qn   [ir       ] += NN*(vgn*du[ir]);
    Qn   [iru      ] += NN*(vgn*du[iru]);
    Qn   [irv      ] += NN*(vgn*du[irv]);
    Qn   [irw      ] += NN*(vgn*du[irw]);
    Qn   [irE      ] += NN*(vgn*du[irE]);
    for (iP=5; iP<sr; iP++){
        
        // state rank index of passive scalar
        //iP = iPass[i];
        
        return xf_NOT_SUPPORTED;
        drk = UR[iP] - UL[iP];
        
        // flux assembly
        Qn[iP] += NN*(vgn*drk);
        
    } // i
    
  /*
    ci    = sqrt(c2);
    ci1   = 1/ci;
    gmic2 = (Gamma-1.)*0.5*ci1*ci1;
    ucg   = ucp * ci * gmi1;
    
    G1    = gmi*(af*du[ir] - ui*du[iru] - vi*du[irv] - wi*du[irw]+du[irE]);
    G2    = -ucp*du[ir]+du[iru]*n1[0]+du[irv]*n1[1]+du[irw]*n1[2];
    
    //trick starts: step 1-project Qn onto char-space
    Qn_u =
    alpha[0] = gmic2 * ((af-ucg) * Qn[0] + ci*gmi1 * );
    */
    
    //if(c2 <= 0)
    if(xfe_False)
    {   //never give up return xf_Error(xf_NON_PHYSICAL);
        //this situation means Roe average does not exist
        vgn = 0.5*max(Mal * al + al, Mar * ar + ar);
 /*       if(enpr < 0.5 || enpl < 0.5)
        {
            Qn   [ir       ] += NN*(vgn*du[ir]);
            Qn   [iru      ] += NN*(vgn*du[iru]);
            Qn   [irv      ] += NN*(vgn*du[irv]);
            Qn   [irw      ] += NN*(vgn*du[irw]);
            Qn   [irE      ] += NN*(vgn*du[irE]);
        }
  */      
        
        
        
        //Opt 1: simple trick
        eps = 1.e-13;
        
        //test sign consistency
        max_char = 0.;
        for(iP=0; iP<sr; iP++)
        {
            if(fabs(Qn[iP]) < eps)
            {
                alpha[iP] = 0.;
                continue;
            }
            
            tmp = Qn[iP]/( (real)sign(du[iP]) * (fabs(du[iP]) + eps));
            
            if(tmp< 0.){
             //sorry, this idea does not work; just add dissiptation as usual
               vgn = 0.5*max(Mal * al + al, Mar * ar + ar);
                
                Qn   [ir       ] += NN*(vgn*du[ir]);
                Qn   [iru      ] += NN*(vgn*du[iru]);
                Qn   [irv      ] += NN*(vgn*du[irv]);
                Qn   [irw      ] += NN*(vgn*du[irw]);
                Qn   [irE      ] += NN*(vgn*du[irE]);
                break;
            }
            else
            {
                alpha[iP] = tmp;
                if(tmp > max_char)
                    max_char = tmp;
            }
        }
        
        if(iP == sr){
            //this case means all signs are consistent
            vgn = 0.5*max(Mal * al + al, Mar * ar + ar);
            for(i=0; i<sr; i++)
                Qn[i] = NN * vgn *  du[i];
        }


  /*
        if(iP == sr){
            //this case means all signs are consistent
            vgn = 0.5*max(Mal * al + al, Mar * ar + ar);
            if(max_char < vgn){
                for(i=0; i<sr; i++)
                    Qn[i] = NN*(vgn * du[i]);
            
                xf_printf("I am good\n");
            }
            else
            {
                if(max_char > 2. *vgn)
                {
                    //give up
                    for(i=0; i<sr; i++)
                        Qn[i] += NN*(vgn * du[i]);
                }
                else
                {
                    //not give up
                    for(i=0; i<sr; i++)
                        Qn[i] = NN*(max_char * du[i]);
                }
            }
        }
    */
    
        
    }
    
    if(xfe_True)
    {
        ci    = sqrt(c2);
        ci1   = 1/ci;
        
        // eigenvalues
        //lam[0] = ucp + ci;
        //lam[1] = ucp - ci;

        if(LM_Fix_flag){
            Ma = max(Mal, Mar);
            if(Ma>=0.8)
                Ma = 1.0;

            //cut-off
            if(Ma<0.01)
                Ma = 0.01;
        }

        //Ma = 1.;
        lam[0] = ucp + ci * Ma;
        lam[1] = ucp - ci * Ma;
        lam[2] = ucp;
        
        // Entropy fix
         ep     = 1e-2;
         eps    = ep*ci;
         
         for (i=0; i<3; i++)
         if ((lam[i]<eps)&&(lam[i]>-eps)){
         eps1 = 1/eps;
         
         lam[i] = 0.5*(eps+lam[i]*lam[i]*eps1);
         }
         
        for (i=0; i<3; i++)
            if (lam[i]<0)
                el[i] = -1;
            else
                el[i] =  1;
        
        // The following parameters are described in the theory guide
        s1    = 0.5*(el[0]*lam   [0]+el[1]*lam   [1]);
        s2    = 0.5*(el[0]*lam   [0]-el[1]*lam   [1]);
        
        // third eigenvalue, absolute value
        l3     = el[2]*lam    [2];
        
        G1    = gmi*(af*du[ir] - ui*du[iru] - vi*du[irv] - wi*du[irw]+du[irE]);
        G2    = -ucp*du[ir]+du[iru]*n1[0]+du[irv]*n1[1]+du[irw]*n1[2];
        
    /*    if(LM_Fix_flag){
            Ma = max(Mal, Mar);
            if(Ma>=1.0)
                Ma = 1.0;
            G2 *= Ma;
        }
    */    
        c21   = ci1*ci1;
        ct1   = (s2   *G1            *ci1);
        ct2   = (s2   *G2            *ci1);
        C1    = ct2   +(s1 - l3   )*G1*c21;
        C2    = ct1   +(s1 - l3   )*G2;
        
        /*    Qn[ir]   = max( Qn[ir],  NN*0.5*(l3*du[ir ]+C1) );
         Qn[iru]  = max( Qn[iru], NN*0.5*(l3*du[iru]+C1   *ui+C2   *n1[0]));
         Qn[irv]  = max( Qn[irv], NN*0.5*(l3*du[irv]+C1   *vi+C2   *n1[1]));
         Qn[irw]  = max( Qn[irw], NN*0.5*(l3*du[irw]+C1   *wi+C2   *n1[2]));
         Qn[irE]  = max( Qn[irE], NN*0.5*(l3*du[irE]+C1   *hi+C2   *ucp  ));
         */
        Qn[ir]   += ( NN*0.5*(l3*du[ir ]+C1) );
        Qn[iru]  += ( NN*0.5*(l3*du[iru]+C1   *ui+C2   *n1[0]));
        Qn[irv]  += ( NN*0.5*(l3*du[irv]+C1   *vi+C2   *n1[1]));
        Qn[irw]  += ( NN*0.5*(l3*du[irw]+C1   *wi+C2   *n1[2]));
        Qn[irE]  += ( NN*0.5*(l3*du[irE]+C1   *hi+C2   *ucp  ));
        
        
        /*** Take care of passive scalars ***/
        //xf_PassiveEulerScalars(PIS, &nPass, iPass);
        
        for (iP=5; iP<sr; iP++){
            
            // state rank index of passive scalar
            //iP = iPass[i];
            
            // left and right scalars
            kl  = UL[iP]*rhol1;
            kr  = UR[iP]*rhor1;
            drk = UR[iP] - UL[iP];
            
            // Roe-averaged scalar
            ki     = (di*kr+kl)*d1;
            
            // flux assembly
            //Qn[iP] = max(Qn[iP], NN*0.5*(l3*drk + C1*ki));
            Qn[iP] += ( NN*0.5*(l3*drk + C1*ki));
            
        } // i
        
    }
    
    return xf_OK;
    
}

/******************************************************************/
// Entrance for Dissipation Correct combined use with Interior Panelty method
// note: modified the diffusion flux with dissipation in Rieman flux
int
DissipationCorrect(const int nq, const int sr, const int dim, const real *qUL, const real *qUR,
                   const real *qn, const real Gamma, real *qQn, real *qFv, const real eta) 
{
    int ierr, i, j, sn;
    const real *UL, *UR, *n;
    real *Qn, *Fv;
    real Qn_sum[MAXSR], Fv_sum[MAXSR];

    for(j=0; j<sr; j++){
       Qn_sum[j] = 0.;
       Fv_sum[j] = 0.;
    }

    for (i=0; i<nq; i++){
        
        UL   = qUL + i*sr;
        UR   = qUR + i*sr;
        n    = qn  + i*dim;
        Qn   = qQn + i*sr;
        Fv   = qFv + i*sr; 
        
        if(dim == 2){
            //ierr = xf_Error(ConvFluxFace2D(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            //ierr = xf_Error(ConvFluxFace2D_Rusanov(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            //ierr = xf_Error(ConvFluxFace2D_Roe(sr, UL, UR, n, NULL, NULL, F, Gamma, MaxCharSpeed));
            //if (ierr != xf_OK) return ierr;

            ierr = xf_Error(ConvFluxFace2D_Adaptive_Roe(sr, UL, UR, n, Gamma, Qn, eta));
            if (ierr != xf_OK) return ierr;
        }
        else if(dim == 3){
            //ierr = xf_Error(ConvFluxFace3D(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            //if (ierr != xf_OK) return ierr;
            //   ierr = xf_Error(ConvFluxFace3D_Rusanov(sr, UL, UR, n, F, Gamma, MaxCharSpeed));
            //   if (ierr != xf_OK) return ierr;
            //ierr = xf_Error(ConvFluxFace3D_Dissip_Roe(sr, UL, UR, n, Gamma, Fv, eta));
            //if (ierr != xf_OK) return ierr;
            
            ierr = xf_Error(ConvFluxFace3D_Adaptive_Roe(sr, UL, UR, n, Gamma, Qn, eta));
            if (ierr != xf_OK) return ierr;
        }
        else
            return xf_Error(xf_OUT_OF_BOUNDS);
      
        for(j=0; j<sr; j++){
           Qn_sum[j] += Qn[j];
           Fv_sum[j] += Fv[j];
        }
    }
    
    //dissiption subtraction based on face-integration
    for(j=0; j<sr; j++)
          for(i=0; i<nq; i++)
          {
             //opt 2: a tricky trick   
             //*(qQn + i*sr + j) += *(qFv + i*sr + j);
                

          }
 //      if(fabs(Qn_sum[j]) < fabs(Fv_sum[j]) || sign(Qn_sum[j]) != sign(Fv_sum[j]))
       {
 /*         sn = sign(Fv_sum[j]);

          if(sn == sign(Qn_sum[j]))
          {
             if(fabs(Qn_sum[j]) < fabs(Fv_sum[j]))
             {
                for(i=0; i<nq; i++)
                   *(qQn + i*sr + j) += *(qFv + i*sr + j);
             }
          }
          else
          {
             for(i=0; i<nq; i++)
                *(qQn + i*sr + j) += *(qFv + i*sr + j);
          }
         
          for(i=0; i<nq; i++)
             if(j == 0)
                *(qQn + i*sr + j) += *(qFv + i*sr + j);
             else{
             sn = sign(*(qFv + i*sr + j)); 
             if(sn == sign(*(qQn + i*sr + j)))
             *(qQn + i*sr + j) = max(fabs(*(qQn + i*sr + j)),  fabs(*(qFv + i*sr + j))) * (real)sn;
             else
                *(qQn + i*sr + j) += *(qFv + i*sr + j);
             }
 */      }
       
    return xf_OK;
}
