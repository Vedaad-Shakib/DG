//routine used for initialization


/******************************************************************/
//when initial condition is obtained from a ASCII file, use the routine
//to read data and execute initialization
//Note: the subroutine is designed to handle ordered data; right now plot3D 
//data format is implemented; the data is loaded into a pointer of Model;
//the data is then released after initialization is finished.
int
init_data_read_from_file(Yu_datafile **pdataholder)
{
   int ierr, i, j, k, var, sr;
   char key[xf_MAXSTRLEN], *value;
   char dummy[xf_MAXSTRLEN], *line;
   FILE *fid;
   Yu_datafile *dataholder;

   ierr = xf_Error(xf_Alloc( (void **) pdataholder, 1, sizeof(Yu_datafile)));
   if (ierr != xf_OK) return ierr;
   dataholder = (*pdataholder);

   //user specify the data dimension
   sr = 8;  //number of vars in file
   dataholder->rank   = sr;

   //for HIT box
 //  /*
   dataholder->dim[0] = 129;
   dataholder->dim[1] = 129;
   dataholder->dim[2] = 129;
   for(i=0; i<3; i++)
   {
      dataholder->xyz_min_max[i][0] = -M_PI;
      dataholder->xyz_min_max[i][1] = M_PI;
   }
 //  */
   //for channel flow
  /*
   dataholder->dim[0] = 48;
   dataholder->dim[1] = 48;
   dataholder->dim[2] = 48;
      
   dataholder->xyz_min_max[0][0] = 0.;
   dataholder->xyz_min_max[0][1] = 2.*M_PI;
   dataholder->xyz_min_max[1][0] = -1.;
   dataholder->xyz_min_max[1][1] = 1.;
   dataholder->xyz_min_max[2][0] = 0.;
   dataholder->xyz_min_max[2][1] = M_PI;
   */

   //turbulence box with Peter Ma 
/*   sr = 6;  //initial data only include
   dataholder->rank = sr;
   dataholder->dim[0] = 193;
   dataholder->dim[1] = 193;
   dataholder->dim[2] = 193;
 */  for(i=0; i<3; i++)
   {
      dataholder->xyz_min_max[i][0] = -M_PI;
      dataholder->xyz_min_max[i][1] = M_PI;
   }


   dataholder->data = NULL;

   ierr = xf_Error(xf_Alloc((void **) &(dataholder->data), (dataholder->dim[0])*
            (dataholder->dim[1])*(dataholder->dim[2])*sr, sizeof(real)));
   if(ierr != xf_OK) return ierr;
      
   //check default input data file name
   fid = fopen("lydg_init.dat", "r");  

   if(fid == NULL)
      return xf_INPUT_ERROR; 

   //seem ready to start loading
   xf_printf("start loading data for initialization...\n");
   do{
      if (fgets(dummy, xf_MAXLINELEN, fid) == NULL)
         continue;
      line = dummy; //set pointer

         //if (xf_TrimAndCheckBlank(&line, xf_MAXLINELEN)) continue;
      
      value = strstr(line, "DT"); 
      if(value != NULL)
         break;
         
   }while(1);
 
      //ierr = xf_Error(xf_Alloc((void **) &(dataholder->data), (dataholder->dim[0])*
      //                         (dataholder->dim[1])*(dataholder->dim[2])*sr, sizeof(real)));
      //if(ierr != xf_OK) return ierr;

      for(k=0;k<dataholder->dim[2];k++)
         for(j=0;j<dataholder->dim[1];j++)
            for(i=0;i<dataholder->dim[0];i++)
            {
               for(var=0; var<sr; var++)
               fscanf(fid, "%lf", &(dataholder->data[k*(dataholder->dim[1])*(dataholder->dim[0])*sr +
                      j*(dataholder->dim[0])*sr + i*sr + var]));
            }
    
   fclose(fid);

      return xf_OK;
} 

/******************************************************************/
//after data is loaded into a struct in Model, it is now assigned to spatial
//location with interpolation
//might require modification by users
static int
Yu_FlowFieldInit_from_datafile(Yu_datafile *holder, const int sr, const int dim, 
                               const real *xglob, real *FV)
{
   int ierr, Iz[2], Iy[2], Ix[2];
   int i, j, k, tmp1, tmp2, tmp3;
   real z0, y0, x0, zSpan, ySpan, xSpan, coor[3];
   real start, end, interp_once[4][10], interp_twice[2][10];  

   if(holder==NULL && dim != 3)
      return xf_CODE_LOGIC_ERROR;

   coor[0] = xglob[0];
   coor[1] = xglob[1];
   coor[2] = xglob[2];
   //check the range 
   if(xglob[0]<holder->xyz_min_max[0][0] || xglob[0]>holder->xyz_min_max[0][1])
      return xf_OUT_OF_BOUNDS;
   if(coor[1]<holder->xyz_min_max[1][0] || coor[1]>holder->xyz_min_max[1][1])
      return xf_OUT_OF_BOUNDS;
   if(xglob[2]<holder->xyz_min_max[2][0] || xglob[2]>holder->xyz_min_max[2][1])
      return xf_OUT_OF_BOUNDS;

   //shifting for z-axis
   //if(coor[2] > holder->xyz_min_max[2][1]) coor[2] -= holder->xyz_min_max[2][1];
   //shifting for x-axis
   //if(coor[0] > holder->xyz_min_max[0][1]) coor[0] -= holder->xyz_min_max[0][1];

   z0 = holder->xyz_min_max[2][0];
   zSpan = (holder->xyz_min_max[2][1] - holder->xyz_min_max[2][0]) / (real) (holder->dim[2] - 1) ;
   y0 = holder->xyz_min_max[1][0];
   ySpan = (holder->xyz_min_max[1][1] - holder->xyz_min_max[1][0]) / (real) (holder->dim[1] - 1) ;
   x0 = holder->xyz_min_max[0][0];
   xSpan = (holder->xyz_min_max[0][1] - holder->xyz_min_max[0][0]) / (real) (holder->dim[0] - 1) ;

   //find upper and lower bounds for z-index
   for(k=0; k<holder->dim[2]; k++)
      if ( z0+(real)k * zSpan <= coor[2] && z0+(real)(k+1) * zSpan >= coor[2])
      { Iz[0]=k; Iz[1]=k+1; break;}
   for(j=0; j<holder->dim[1]; j++)
      if ( y0+(real)j * ySpan <= coor[1] && y0+(real)(j+1) * ySpan >= coor[1])
      { Iy[0]=j; Iy[1]=j+1; break;}
   for(i=0; i<holder->dim[0]; i++)
      if ( x0+(real)i * xSpan <= coor[0] && x0+(real)(i+1) * xSpan >= coor[0])
      { Ix[0]=i; Ix[1]=i+1; break;}

   //interpolate on x-dimension
   tmp1 = (holder->dim[1])*(holder->dim[0])*(holder->rank);
   tmp2 = (holder->dim[0])*(holder->rank);
   tmp3 = holder->rank;
  
   for(k=0; k<2; k++)
   {
      for(j=0; j<2; j++)
      {
         for(i=0; i<tmp3; i++)
         {
            start = holder->data[Iz[k]*tmp1+Iy[j]*tmp2+Ix[0]*tmp3+i];
            end   = holder->data[Iz[k]*tmp1+Iy[j]*tmp2+Ix[1]*tmp3+i];
               
            interp_once[k*2+j][i] = start - (start-end)/xSpan * (coor[0] - (x0+(real)Ix[0]*xSpan));
         }
      }
         
      for(i=0; i<tmp3; i++)
      {
         start = interp_once[k*2+0][i];  end = interp_once[k*2+1][i];
         interp_twice[k][i] = start - (start-end)/ySpan * (coor[1] - (y0+(real)Iy[0]*ySpan)); 
      }
   }

   for(i=3; i<tmp3; i++)
   {
      start = interp_twice[0][i]; end = interp_twice[1][i];

      FV[i-3] = start - (start-end)/zSpan * (coor[2] - (z0+(real)Iz[0]*zSpan));
   
      //with Peter Ma; only interpolate velocity
      //FV[i-3+1] = start - (start-end)/zSpan * (coor[2] - (z0+(real)Iz[0]*zSpan));
   }
   
   //correct internal energy 
   //FV[4] = 5785.7143/0.4 + 0.5*(FV[1]*FV[1] + FV[2]*FV[2] + FV[3]*FV[3])/FV[0];
   //FV[0] = 1.;
   //FV[4] = 214.2857144/0.4 + 0.5*(FV[1]*FV[1] + FV[2]*FV[2] + FV[3]*FV[3])/FV[0];
   
   if(sr == 6)
   //for the unused scalar
   FV[5] = 0.0;

   return xf_OK;
}

/******************************************************************/
// initialization of flow fluid for numerical test
static int
Yu_FlowFieldInit(const int FcnOpt, const int sr, const int dim, const real *xglob, real *FV)
 {
   int ierr, k, i, j;
   real x, y, z, r;
   real rho, u, v, w, p, Y;
   real Rconst = 8.314;
   real Gamma = 1.4;
   real V0, T0, L;
   real phi, alpha, x0, y0, Ms, rho2, p2, tau; 
   real mu0 = 0.; //Model->mu_c;

   //which case we are testing
   switch(FcnOpt) {
      case Yu_AdvDiffRectChannel:
         for(i=0; i<sr; i++)
            FV[i] = 0.;

         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1];

         L = 3; 
         //r1
         x0 = AdvSpd/2./LapVis + sqrt(pow(AdvSpd,2.)/4./pow(LapVis,2.) + pow(PI,2.));
         //r2
         y0 = AdvSpd/2./LapVis - sqrt(pow(AdvSpd,2.)/4./pow(LapVis,2.) + pow(PI,2.));

         FV[0] = sin(PI*y) * (y0*exp(x0*x+y0*L) - x0*exp(x0*L+y0*x))/(y0*exp(y0*L) - x0*exp(x0*L));
      break;

      case Yu_PoiseuilleFlow:

         return xf_Error(xf_INPUT_ERROR);
         //require two parameters
         V0 = 1.;
         L  = 1.;

         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1];

         for (k=0; k<sr; k++) FV[k] = 0;
         FV[0] = 1.0;
         FV[1] = -0.5/mu0*V0*(1.0 - y*y);
         p = 1.0 + V0*x;
         FV[3] = p/(Gamma-1.0) + 0.5*FV[1]*FV[1]/FV[0];

      break;

      case Yu_TaylorGreenVortex:
         //condition taken from high-order workshop  
         //from Carton et al. Int. J. Numer. Meth. Fluids
         //L = 1.; V0 = 1.; mu0 = 1./1600.;
         //rho0 = 1.;

         //required parameters
         V0 = 1.; 
         L  = 1.;

         if(dim != 3) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1]; z = xglob[2];

         p = 1.0 + 0.014 * V0 * V0 / 16.0 * 
            (cos(2.0*x/L)+cos(2.0*y/L)) * (cos(2.0*z/L)+2.0);
         FV[0] = 0.014;
         FV[1] =  FV[0] * V0 * sin(x/L)*cos(y/L)*cos(z/L);
         FV[2] = -FV[0] * V0 * cos(x/L)*sin(y/L)*cos(z/L);
         FV[3] = 0.0;
         FV[4] = p/(Gamma-1.0) + 0.5*(FV[1]*FV[1] + FV[2]*FV[2])/FV[0];

         FV[5] = 0.0;
         break;

      case Yu_DoubleMachReflection:
         for(i=0; i<sr; i++)
            FV[i] = 0.;

         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1];
         
         if (xglob[0] < 1./6. + xglob[1]/sqrt(3.0))
         {
            //post-shock state
            FV[0] = 8.0; FV[1] = 57.15768; FV[2] = -33.0;
            FV[3] = 563.500024; FV[4] = 8.0;
         }
         else
         {
            //pre-shock state
            FV[0] = 1.4; FV[1] = 0.0; FV[2] = 0.0;
            FV[3] = 2.5; FV[4] = 1.4;
         }
      break;

      case Yu_CompressVortexAdvect:
         
         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1];

         phi = 1.0; alpha = 4.0;
         x0 = 5.0; y0 = 5.0;
         r = sqrt(pow(x-x0, 2.)+pow(y-y0, 2.));

         rho = pow((1.-alpha*alpha*(Gamma-1.)/16./phi/Gamma/PI/PI*exp(2.*phi*(1.-r*r))), 1./(Gamma-1.));
         u   = 1.0-alpha/2./PI*(y-y0)*exp(phi*(1.-r*r));
         v   = 1.0+alpha/2./PI*(x-x0)*exp(phi*(1.-r*r));
         p   = pow(rho, Gamma);

         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v;
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         if(sr > dim + 2)
         FV[4] = 0.0;
      break;
         
      case Yu_acoustic_pulse:
         
         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1];
         
         phi = 1.0; alpha = 4.0;
         x0 = 5.0; y0 = 5.0;
         r = sqrt(pow(x-x0, 2.)+pow(y-y0, 2.));
         
         rho = 1.0 + 0.5 * 2. * exp(-r*r) ; //pow((1.-alpha*alpha*(Gamma-1.)/16./phi/Gamma/PI/PI*exp(2.*phi*(1.-r*r))), 1./(Gamma-1.));
         u   = 0.0; //-alpha/2./PI*(y-y0)*exp(phi*(1.-r*r));
         //v   = 0.5+alpha/2./PI*(x-x0)*exp(phi*(1.-r*r));
         v   = 0.0; //alpha/2./PI*(x-x0)*exp(phi*(1.-r*r));
         p   = 1./Gamma + 0.5 * 2. * exp(-r*r) ; //pow(rho, Gamma);
         
         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v;
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         break;

      case Yu_ShockVortexInteraction:
         
         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1];

         Ms = 1.1;
         rho2 = (Gamma+1) * Ms*Ms / (2.+(Gamma - 1.) * Ms * Ms);
         p2 = 1. + 2.*Gamma*(Ms*Ms - 1.) / (Gamma + 1.);

         phi = 0.204; alpha = 0.3; 
         x0 = 0.25; y0 = 0.5;
         r = sqrt(pow(x-x0,2.)+pow(y-y0,2.));

         tau = r/0.05;
         rho = pow((1.-alpha*alpha*(Gamma-1.)/4./phi/Gamma*exp(2.*phi*(1.-tau*tau))), 1./(Gamma-1.));
         u   = sqrt(Gamma)*Ms + alpha*tau*(y-y0)/r*exp(phi*(1.-tau*tau));
         v   = -alpha*tau*(x-x0)/r*exp(phi*(1.-tau*tau));
         p   = pow(rho, Gamma);

         if(x>0.5)
         {
            FV[0] = rho2; Ms = (1.+0.5*(Gamma-1.)*Ms*Ms) / (Gamma*Ms*Ms - 0.5*(Gamma - 1.)); Ms = sqrt(Ms);
            FV[1] = rho2*Ms*sqrt(1.4 *p2/rho2); FV[2] = 0.0;
            FV[3] = p2/(Gamma-1.) + 0.5 * FV[1] * FV[1] / FV[0];
         }
         else
         {
            FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v;
            FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
          }
         
         FV[4] = 0.0;
         break;

      case Yu_ForwardFacingStep:
         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1];

         FV[0] = 1.4; FV[1] = 4.2; FV[2] = 0.0; FV[3]= 8.8; FV[4] = 1.4;
      //   if(x < 0.6 && x > 0.59 && y < 0.2)
      //   {
      //      FV[1] = 0.0; FV[3] = 2.5;
      //   }

      break;

      case Yu_Kelvin_Helmholtz_Instability:
         
         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0];  y = xglob[1];

         if(fabs(y)>=0.25)
         { rho = 1.0; u = 0.5; v = 0.0; p = 2.5; }

         if(fabs(y)<0.25)
         { rho = 2.0; u = -0.5; v = 0.0; p = 2.5; }

         //add velocity perturbation
         u += 0.01 * sin(2.0*M_PI*x);
         v += 0.01 * sin(2.0*M_PI*x);

         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v; 
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         FV[4] = FV[0];

      break;

      case Yu_acoustic_ML:
         if(dim!=2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0]; y = xglob[1];

         u = 3./2.0 + 1./2.0 * tanh(2.0 * y / 1.0);
         v = 0.0;
         rho = 1.0;
         p = 11.4285714286;
         
         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v; 
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         FV[4] = FV[0];

      break;

      case Yu_2DRiemannProblem2:

         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0];  y = xglob[1];

         if(y>0.0 && x<=0.0)
         { rho = 0.5065; u = 0.8939; v = 0.0; p = 0.35;}
         if(y>0.0 && x>0.0)
         { rho = 1.1; u = 0.0; v = 0.0; p = 1.1;}
         if(y<=0.0 && x<=0.0)
         { rho = 1.1; u = 0.8939; v = 0.8939; p = 1.1;}
         if(y<=0.0 && x>0.0)
         { rho = 0.5065; u = 0.0; v = 0.8939; p = 0.35;}
         
         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v; 
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         FV[4] = FV[0];
      break;

      case Yu_2DRiemannProblem3:

         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0];  y = xglob[1];

         if(y>0.0 && x<=0.0)
         { rho = 2.0; u = 0.75; v = 0.5; p = 1.0;}
         if(y>0.0 && x>0.0)
         { rho = 1.0; u = 0.75; v = -0.5; p = 1.0;}
         if(y<=0.0 && x<=0.0)
         { rho = 1.0; u = -0.75; v = 0.5; p = 1.0;}
         if(y<=0.0 && x>0.0)
         { rho = 3.0; u = -0.75; v = -0.5; p = 1.0;}
         
         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v; 
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         FV[4] = FV[0];
      break;

      case Yu_2DRiemannProblem4:

         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0];  y = xglob[1];

         if(y>0.0 && x<=0.0)
         { rho = 1.0; u = -0.6259; v = 0.1; p = 1.0;}
         if(y>0.0 && x>0.0)
         { rho = 0.5197; u = 0.1; v = 0.1; p = 0.4;}
         if(y<=0.0 && x<=0.0)
         { rho = 0.8; u = 0.1; v = 0.1; p = 1.0;}
         if(y<=0.0 && x>0.0)
         { rho = 1.0; u = 0.1; v = -0.6259; p = 1.0;}
         
         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v; 
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         FV[4] = FV[0];
      break;
      
      case Yu_2DRiemannProblem5:

         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0];  y = xglob[1];

         if(y>0.0 && x<=0.0)
         { rho = 1.0; u = 0.7276; v = 0.0; p = 1.0;}
         if(y>0.0 && x>0.0)
         { rho = 0.5313; u = 0.0; v = 0.0; p = 0.4;}
         if(y<=0.0 && x<=0.0)
         { rho = 0.8; u = 0.0; v = 0.0; p = 1.0;}
         if(y<=0.0 && x>0.0)
         { rho = 1.0; u = 0.0; v = 0.7276; p = 1.0;}

         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v; 
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         FV[4] = FV[0];
      break;

      case Yu_2DSodTube:
         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0];  y = xglob[1];

         if(x < 0.5)
         { rho = 1.0; u = 0.0; v = 0.0; p = 1.0;}
         else
         { rho = 0.125; u = 0.0; v= 0.0; p = 0.1;}

         FV[0] = rho; FV[1] = rho*u; FV[2] = rho*v; 
         FV[3] = p/(Gamma-1.) + 0.5*rho*(u*u + v*v);
         FV[4] = FV[0];
      break;

      case Yu_Density_Wave_ConstVel:
        
         x = xglob[0];

         rho = 2.0 + sin(2.0*M_PI*x);
         u = 1.0;
       
         FV[0] = rho; FV[1] = rho*u; FV[2] = 0.0; 
         if(dim == 3) FV[dim] = 0.0;
         FV[dim+1] = 2.0/(Gamma-1.0) + 0.5*rho*u*u;

      break;
      
      case Yu_plane_acoustic_test:
        
         x = xglob[0];

         rho = 1.0;
         p = 1./Gamma * (1.0 + 0.000001*sin(2.0*M_PI*x * 2.0));
         u = 0.0;
       
         FV[0] = rho; FV[1] = rho*u; FV[2] = 0.0; 
         if(dim == 3) FV[dim] = 0.0;
         FV[dim+1] = p/(Gamma-1.0) + 0.5*rho*u*u;

      break;

      case Yu_quasi1D_Detonation:
         if(dim != 2) return xf_Error(xf_INPUT_ERROR);
         x = xglob[0];
         
         if(x> 2.0) 
         {
            rho = 1.0; p = 1.0; u = 0.0; Y = 0.0;
         }
         else if(x<=2.0 && x>0.5)
         {
            rho = (7. - 1.754929) /1.5 * (x - 0.5) + 1.754929;
            p = (28.0 - 14.33557) /1.5 * (x - 0.5) + 14.33557;
            u = sqrt(1.24) * 5.0 * (1. - 1./rho);
            Y = -1./1.5 * (x - 0.5) + 1.0;
         }
         else  //left inflow; also equilibirum state
         {
            p = 14.33557; rho = 1.754929;  u = sqrt(1.24) * 5.0 * (1. - 1./rho);
            Y = 1.0;
         }
         
         FV[0] = rho; 
         FV[1] = rho * u; 
         FV[2] = 0.0;
         FV[3] = p/(1.24 - 1.0)  + 0.5*rho*u*u;
         FV[4] = rho * Y;

         //change to initialize using a stationary shock
/*         Ms = 5.0;
         Gamma = 1.24;
         if(x>10.0){
         rho2 = (Gamma+1) * Ms*Ms / (2.+(Gamma - 1.) * Ms * Ms);
         p2 = 1. + 2.*Gamma*(Ms*Ms - 1.) / (Gamma + 1.);
         Ms = (1.+0.5*(Gamma-1.)*Ms*Ms) / (Gamma*Ms*Ms - 0.5*(Gamma - 1.)); Ms = sqrt(Ms);
         FV[0] = rho2;
         FV[1] = rho2*Ms*sqrt(Gamma *p2/rho2);
         FV[2] = 0.0;
         FV[3] = p2/(Gamma-1.) + 0.5 * FV[1] * FV[1] / FV[0];
         }
         else{
         rho = 1.0; p = 1.0; u = sqrt(Gamma) * Ms;
         FV[0] = rho; 
         FV[1] = rho*u; 
         FV[2] = 0.0;
         FV[3] = p/(Gamma - 1.0)  + 0.5*rho*u*u;
         }

         FV[4] = 0.0;
*/
        break;

      case Yu_channel_flow:
        x = xglob[0]; y = xglob[1]; z = xglob[2];
 
        tau = 0.03;
        rho = 1.0;
        //for Re_tau 950
        //u = 21.0 * (1. - y*y) + tau * sin(x)*cos(y*M_PI)*cos(z * 2.);
        
        //for Re_tau 180
        u = 16.5 * (1. - y*y) + tau * sin(x)*cos(y*M_PI)*cos(z * 2.);
        v = - tau * cos(x)*sin(y*M_PI)*cos(z * 2.);
        w = 0.;
        //u = 21.0 * (1. - y*y/M_PI/M_PI)  +  ((real) rand())/((real) RAND_MAX);
        //v = ((real) rand())/((real) RAND_MAX);
        //w = ((real) rand())/((real) RAND_MAX);

        //for Re_tau 950; Ma 0.1
        //p = 31500.0;
        //
        //for Re_tau 180; Ma 0.2
        p = 4861.60714286;

        FV[0] = rho; 
        FV[1] = rho * u;  
        FV[2] = rho * v;
        FV[3] = rho * w;
        FV[4] = p/(1.4 - 1.0)  + 0.5*rho*(u*u + v*v + w*w);

        break;

      case Yu_BCstudy_acoustic_pulse:
         x = xglob[0]; y = xglob[1];

         phi = 0.1; 
         alpha = 0.1;
         x0 = 0.1; y0 = 0.0;
         r = sqrt(pow(x-x0, 2.)+pow(y-y0, 2.));
         
         r /= phi;

         rho = 1.0 + alpha * exp(- r * r);
         u = 0.4;
         v = 0.0;
         p = 1./Gamma + alpha * exp(- r * r);

        FV[0] = rho; 
        FV[1] = rho * u;  
        FV[2] = rho * v;
        FV[3] = p/(1.4 - 1.0)  + 0.5*rho*(u*u + v*v);

        break;

        case Yu_BCstudy_vortex_transport:
           x = xglob[0]; y = xglob[1];
           
           phi = 0.1;
           alpha = 2.0;
           x0 = 0.1; y0 = 0.0;
           r = sqrt(pow(x-x0, 2.)+pow(y-y0, 2.));
           
           r /= phi;
           
           rho = pow((1.-alpha*alpha*(Gamma-1.)/16./phi/Gamma/PI/PI*exp(2.*phi*(1.-r*r))), 1./(Gamma-1.));
           u   = 1.0-alpha/2./PI*(y-y0)*exp(phi*(1.-r*r));
           v   = 0.0+alpha/2./PI*(x-x0)*exp(phi*(1.-r*r));
           p   = pow(rho, Gamma);
           
           //rho = 1.0 + alpha * exp(- r * r);
           //u = 0.4;
           //v = 0.0;
           //p = 1./Gamma + alpha * exp(- r * r);
           
           FV[0] = rho;
           FV[1] = rho * u;
           FV[2] = rho * v;
           FV[3] = p/(1.4 - 1.0)  + 0.5*rho*(u*u + v*v);
           
        break;
  
      case Yu_debug:
        rho = 1.0;
        x = xglob[0];
        y = xglob[1]; 
        z = xglob[2];

        u = x;
        v = y; 
        w = z;

           
        FV[0] = 1.;
        FV[1] = rho * u;
        FV[2] = rho * v;
        FV[3] = rho * w;
        FV[4] = 1.4/(1.4-1.) + 0.5 * rho * (u*u + v*v + w*w); 
       
        break;
           
      default:
        return xf_Error(xf_NOT_SUPPORTED);
        break;
   }

   return xf_OK;
}
