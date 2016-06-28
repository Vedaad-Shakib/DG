//The quantities we would need for post-processing
#include "xf.h"
#include "xf_AllStruct.h"
#include "xf_MPI.h"
#include "xf_Memory.h"
#include "xf_Basis.h"
#include "xf_All.h"
#include "xf_Data.h"
#include "xf_DataMath.h"
#include "xfYu_Statistics.h"
#include "xf_Math.h"
#include "xf_MeshTools.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/******************************************************************/
// write string into tecplot binary format   
static int
write_string_to_tecplot(FILE *fout, const char *fname, int len) 
{
      int i, tmp; 
      int *StrIntData;

      for(i=0; i<len; i++){
         tmp = (int)fname[i];
         fwrite(&tmp, sizeof(int), 1, fout);
      }

      //the last digit should be 0
      tmp = 0; 
      fwrite(&tmp, sizeof(int), 1, fout);
              
      return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_HydroVorticityMag
static int
Compute_Vorticity_Mag(const int sr, const int dim, const real *U,
                     const real *gU, real *pS)
{
   int ierr, i, j, k;
   int ij, ji, ir, irV[3];
   real gu[9]; // velocity gradient
   real V[3];
   real S, r, t;
   
   
   // state indices
   ir  = 0;
   for (i=0; i<dim; i++) irV[i] = 1 + i;
   
   // state components
   r  = U[ir ];  // density
   for (i=0; i<dim; i++) V[i] = U[irV[i]]/r; // velocity
   
   // gu = velocity gradient
   for (i=0; i<dim; i++)
      for (j=0; j<dim; j++)
         gu[i*dim+j] = gU[i*sr+irV[j]]/r - V[j]*gU[i*sr+ir]/r;
   
   // S = magnitude of vorticity
   S = 0.;
   for (i=0; i<dim; i++)
      for (j=i+1; j<dim; j++){
         ij = i*dim+j;
         ji = j*dim+i;
         t = gu[ij] - gu[ji];
         S += t*t;
      }// j
   
   if (pS != NULL) (*pS) = sqrt(S);
   
   return xf_OK;
   
}

/******************************************************************/
//   FUNCTION Definition: xf_HydroVorticityMag
static int
Compute_Divergence_Mag(const int sr, const int dim, const real *U,
                       const real *gU, real *pS)
{
    int ierr, i, j, k;
    int ij, ji, ir, irV[3];
    real gu[9]; // velocity gradient
    real V[3];
    real S, r, t;
    
    
    // state indices
    ir  = 0;
    for (i=0; i<dim; i++) irV[i] = 1 + i;
    
    // state components
    r  = U[ir ];  // density
    for (i=0; i<dim; i++) V[i] = U[irV[i]]/r; // velocity
    
    // gu = velocity gradient
    for (i=0; i<dim; i++)
        for (j=0; j<dim; j++)
            gu[i*dim+j] = gU[i*sr+irV[j]]/r - V[j]*gU[i*sr+ir]/r;
    
    // S = magnitude of vorticity
    S = 0.;
    for (i=0; i<dim; i++)
        for (j=0; j<dim; j++)
        {     
            if(i==j)
                S += gu[i*dim+j];
        }

    //if (pS != NULL) (*pS) = sqrt(S);
    if (pS != NULL) (*pS) = S;
   
    return xf_OK;
    
}

/******************************************************************/
int 
Yu_InitOutput(Yu_Output *Output)
{
   Output->egrp = 0;
   Output->elem = 0;
   Output->Value = 0.0;
   Output->Type = 0;
   Output->nVars = 1; //at least set one variable as default
   Output->UsesFlux = xfe_False;
   Output->FluxComponentNames = NULL;
   Output->FluxComponentRanks = NULL;
   Output->FluxComponentWeights = NULL;
   Output->FluxComponentMoments = NULL;

   //imp: pointwise output specification
   Output->nPoints = 0;
   Output->pretime = 0.0;
   Output->accumtime = 0.0;
   Output->RestartPointStat = xfe_False;
   Output->egrp = NULL; //global index; might not need
   Output->elem = NULL; //global index; might not need
   Output->xref = NULL;
   Output->egrpLocal = NULL;
   Output->elemLocal = NULL;
   Output->PointData = NULL;
   Output->SampleTimeInv = 0.0;
   Output->OutputFile_Freq_Ratio_DataFile = 1;
   Output->File_Write_Offset = 0;
   Output->SequanceDump = xfe_False;

   return xf_OK;
}

/******************************************************************/
//function: specify the location of sampling pointwise output
//locations should be specified in analytical way but can be multi-D
static int SpecifyPointLocation(Yu_Output * Output, enum xfe_Bool WhtChkDim, int *Dim)
{
   int i, j, k, dim, totsize, ierr;
   int xdim, ydim, zdim;  //for Cartesian system
   real xorg, yorg, zorg;
   real xlen, ylen, zlen, rad, theta, incre; 
   real *pdata;
   int rdim, adim;  //for cylindrical system

   //logics check 
   if(Output->Type != xfe_PointValue)
      return xf_CODE_LOGIC_ERROR;

   dim = 3; // for generality
   //output name as an identifier
   if(strcmp(Output->Name, "StatPlane") == 0)
   {
      xdim = 121; ydim = 1; zdim = 73;
     
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = ydim; Dim[2] = zdim;
         return xf_OK;
      }

      xlen = 12.0; ylen = 0.0; zlen = 36.0;
      xorg = -6.0; yorg = 0.0; zorg = 0.0;
      Output->nPoints = xdim * ydim * zdim;
   }

   if(strcmp(Output->Name, "PlaneAcoustic") == 0)
   {
      xdim = 36; rdim = 1; zdim = 361;
         
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = rdim; Dim[2] = zdim;
         return xf_OK;
      }

      zlen = 36.0; rad = 2.0; incre = 4.0 / zlen;
      Output->nPoints = xdim * rdim * zdim;
   }

   if(strcmp(Output->Name, "Padapt_Test") == 0)
   {
      xdim = 1; ydim = 1; zdim = 1;
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = ydim; Dim[2] = zdim;
         return xf_OK;
      }
      
      Output->nPoints = 1;
   }

   if(strcmp(Output->Name, "acoustic_ML_Stat") == 0)
   {
      xdim = 301; ydim = 21; zdim = 1;
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = ydim; Dim[2] = zdim;
         return xf_OK;
      }
      
      Output->nPoints = xdim * ydim * zdim;
   }
   
   if(strcmp(Output->Name, "acoustic_ML_pProbe") == 0)
   {
      xdim = 8; ydim = 1; zdim = 1;
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = ydim; Dim[2] = zdim;
         return xf_OK;
      }
      
      Output->nPoints = xdim * ydim * zdim;
   }

   //a special output for data pointwise data in regular domain
   if(strcmp(Output->Name, "data_dump") == 0)
   {
     //for initialization of channel flow
     // xdim = 48; ydim = 48; zdim = 48;
      xdim = 24; ydim = 64; zdim = 24;
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = ydim; Dim[2] = zdim;
         return xf_OK;
      }
      
      Output->nPoints = xdim * ydim * zdim;
   }
   
   //a special output for data pointwise data in regular domain
   if(strcmp(Output->Name, "box_HIT_dump") == 0)
   {
     //for initialization of channel flow
      xdim = 128; ydim = 128; zdim = 128;
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = ydim; Dim[2] = zdim;
         return xf_OK;
      }
      
      Output->nPoints = xdim * ydim * zdim;
   }

   if(strcmp(Output->Name, "acoustic_circle") == 0)
   {
      xdim = 300; rdim = 1; zdim = 1;
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = rdim; Dim[2] = zdim;
         return xf_OK;
      }

      Output->nPoints = xdim * rdim * zdim;
   }

   //specify cylindrical region for data dump
   if(strcmp(Output->Name, "cylind_data_dump") == 0)
   {
      xdim = 300; rdim = 300; zdim = 1;
      if(WhtChkDim)
      {
         Dim[0] = xdim; Dim[1] = rdim; Dim[2] = zdim;
         return xf_OK;
      }
     
      incre = (32. - 0.5)/(real) (xdim - 1.);
      Output->nPoints = xdim * rdim * zdim;
   }

   //allocate memory for this output
   ierr = xf_Error(xf_Alloc((void **) &Output->egrpLocal, Output->nPoints, sizeof(int)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc((void **) &Output->elemLocal, Output->nPoints, sizeof(int)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc((void **) &Output->xref, Output->nPoints*dim, sizeof(real)));
   if (ierr != xf_OK) return ierr;

   //total size of the data
   totsize = Output->nPoints*(dim + Output->nVars);
   ierr = xf_Error(xf_Alloc((void **) &Output->PointData, totsize, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   
   for(i=0; i<totsize; i++)
      Output->PointData[i] = 0.0;

   if(strcmp(Output->Name, "StatPlane") == 0)
   {
      //set coordinates for each point
      for(i=0; i<xdim; i++)
         for(j=0; j<ydim; j++)
            for(k=0; k<zdim; k++)
            {
               pdata = Output->PointData + (i*ydim*zdim + j*zdim + k)*(dim + Output->nVars);
               if(xdim == 1)
                  pdata[0] = xorg + (real)i * xlen/(real)xdim;
               else
                  pdata[0] = xorg + (real)i * xlen/(real)(xdim-1);
               if(ydim == 1)
                  pdata[1] = yorg + (real)j * ylen/(real)ydim;
               else
                  pdata[1] = yorg + (real)j * ylen/(real)(ydim-1);
               if(zdim == 1)
                  pdata[2] = zorg + (real)k * zlen/(real)zdim;
               else
                  pdata[2] = zorg + (real)k * zlen/(real)(zdim-1);
            
            }
   }
  
   if(strcmp(Output->Name, "PlaneAcoustic") == 0)
   {
      //set coordinates for each point
      for(i=0; i<xdim; i++)  //azithmul
         for(j=0; j<rdim; j++)
            for(k=0; k<zdim; k++)
            {
               pdata = Output->PointData + (i*rdim*zdim + j*zdim + k)*(dim + Output->nVars);
               theta = (real) i * 2.*M_PI/(real)xdim;
               pdata[2] = 0.0 + (real)k * zlen/(real)(zdim-1);
              
               pdata[0] = cos(theta) * (rad + incre*pdata[2]); 
               pdata[1] = sin(theta) * (rad + incre*pdata[2]);
            }
   }

      if(strcmp(Output->Name, "acoustic_circle") == 0)
      {
         for(i=0; i<xdim; i++)  //azithmul
            for(j=0; j<rdim; j++) 
               for(k=0; k<zdim; k++) 
               {
                  pdata = Output->PointData + (i*rdim*zdim + j*zdim + k)*(dim + Output->nVars);
                  theta = (real) i * 2.*M_PI/(real)xdim;
                  rad = 12.9;
                  pdata[0] = 1.86 + cos(theta) * rad; 
                  pdata[1] = sin(theta) * rad;
                  pdata[2] = 0.;
               }
      }

      if(strcmp(Output->Name, "cylind_data_dump") == 0)
      {
                     
         //set coordinates for each point
         for(i=0; i<xdim; i++)  //azithmul
            for(j=0; j<rdim; j++) 
               for(k=0; k<zdim; k++) 
               {    
                  pdata = Output->PointData + (i*rdim*zdim + j*zdim + k)*(dim + Output->nVars);
                  theta = (real) i * 2.*M_PI/(real)xdim;
                  pdata[2] = 0.0 + (real)k * zlen/(real)(zdim-1);
                  pdata[0] = cos(theta) * (rad + incre*pdata[2]); 
                  pdata[1] = sin(theta) * (rad + incre*pdata[2]);
               }    
      }


   if(strcmp(Output->Name, "Padapt_Test") == 0)
   {
      // on 2d domain [0, 10] X [0, 10]
      pdata = Output->PointData;
      pdata[0] = 5.0; pdata[1] = 5.0; pdata[2] = 0.;
   }

   if(strcmp(Output->Name, "acoustic_ML_Stat") == 0)
   {
      //set coordinates for each point
      for(i=0; i<xdim; i++)
         for(j=0; j<ydim; j++)
            for(k=0; k<zdim; k++)
            {
               pdata = Output->PointData + (i*ydim*zdim + j*zdim + k)*(dim + Output->nVars);
           
               pdata[0] = 20. + (real)i * 300.0 / (real)(xdim-1);
               pdata[1] = -10.0 + (real)j * 20.0 / (real)(ydim-1);
               pdata[2] = 0.; 
            }
   }
  
   if(strcmp(Output->Name, "data_dump") == 0)
   {
      //set coordinates for each point
      for(i=0; i<xdim; i++)
         for(j=0; j<ydim; j++)
            for(k=0; k<zdim; k++)
            {
               pdata = Output->PointData + (i*ydim*zdim + j*zdim + k)*(dim + Output->nVars);

               //for initialization of channel flow
               pdata[0] = 2. * M_PI / (real)xdim +  (real) i * 2. * M_PI / (real)xdim;
               pdata[1] = - cos(M_PI * (real) j / (real)(ydim - 1)); 
               pdata[2] = M_PI / (real)zdim + (real) k * M_PI / (real)zdim;

               //pdata[0] =  4. * M_PI/(real)xdim + (real) i * 4. * M_PI / (real)xdim;
               //pdata[1] = - cos(M_PI * (real) j / (real)(ydim - 1));
               //pdata[2] =  2. * M_PI/(real)zdim  + (real) k * 2. * M_PI/(real)zdim;
            }
   }
   
   if(strcmp(Output->Name, "box_HIT_dump") == 0)
   {
      //set coordinates for each point
      for(i=0; i<xdim; i++)
         for(j=0; j<ydim; j++)
            for(k=0; k<zdim; k++)
            {
               pdata = Output->PointData + (i*ydim*zdim + j*zdim + k)*(dim + Output->nVars);
               
               pdata[0] = -M_PI +  ((real)i + 0.5) * 2. * M_PI / (real)xdim;
               pdata[1] = -M_PI +  ((real)j + 0.5) * 2. * M_PI / (real)ydim; 
               pdata[2] = -M_PI +  ((real)k + 0.5) * 2. * M_PI / (real)zdim;
            }
   }


   if(strcmp(Output->Name, "acoustic_ML_pProbe") == 0)
   {
      for(i=0; i<xdim; i++)
      {
         pdata = Output->PointData + i*(dim + Output->nVars);
         pdata[0] = 60.0 + (real)i * 40.0;
         pdata[1] = 0.0;  pdata[2] = 0.0;
      }
   }
   return xf_OK;
}

/******************************************************************/
//call to initialize pointwise output statistics
// should called with Model structure initializaiton
int
Yu_OutputPointValueInit(xf_All *All, Yu_Output *Output)
{
   //specify the set of point for output statistics
   int ierr, i, j, k, l, dim;
   real xpoint[3], chk[3], *gridx;
   char buf[xf_MAXSTRLEN];
   FILE *fout;
   struct stat st = {0};

   //build EES struct before pointwise searching 
   xf_ElemSearchStruct *ESS = NULL;

   ierr = xf_Error(xf_Alloc((void **) &ESS, 1, sizeof(xf_ElemSearchStruct)));
   if (ierr != xf_OK) return ierr;

   ierr = xf_Error(xf_BuildElemSearchStructure(All, ESS));
   if (ierr != xf_OK) return ierr;

   //get point location info
   ierr = xf_Error(SpecifyPointLocation(Output, xfe_False, NULL));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc((void **) &gridx, 3*(Output->nPoints), sizeof(real)));
   if (ierr != xf_OK) return ierr;

   //find local coordiate of each point in DG cell
   dim = 3;
   for(i=0; i<Output->nPoints; i++)
      for(j=0; j<dim; j++)
         gridx[i*dim+j] = *(Output->PointData + i*(dim + Output->nVars) + j);
      
   ierr = xf_Error(xf_FindElemUsingSearchStructure(All, Output->nPoints, gridx, ESS, 
                   Output->egrpLocal, Output->elemLocal, Output->xref));
   if (ierr != xf_OK) return ierr;

 /*  
   for(i=0; i<Output->nPoints; i++)
   {
      xpoint[0] = *(Output->PointData + i*(dim + Output->nVars) + 0); 
      xpoint[1] = *(Output->PointData + i*(dim + Output->nVars) + 1); 
      xpoint[2] = *(Output->PointData + i*(dim + Output->nVars) + 2); 
         
      //searching point globally and store on local output struct
      //if the point is not in current block; elemLocal=egrpLocal=-1
      ierr = xf_Error(xf_FindElemUsingSearchStructure(All, xpoint, ESS, 
                      Output->egrpLocal+i, Output->elemLocal+i, Output->xref+i*dim));
      if (ierr != xf_OK) return ierr;
   }
*/
   //load data from restart file if specified
   //if just dump sequance data; no need for restart
   if(Output->RestartPointStat && !Output->SequanceDump)
   {
      //check if the target folder exists
      sprintf(buf, "%s_Statistic", Output->Name);

      if (stat(buf, &st) == -1) {
         xf_printf("the restart folder does not exist!\n");
         return xf_CODE_LOGIC_ERROR;
      }

      //logic check: (1) number of points (2) number of vars (3) point coordinates
      sprintf(buf, "%s_Statistic/statistics_%s.bac", Output->Name, Output->Name);
      fout = fopen(buf,"r");
      
      fread(&k, sizeof(int), 1, fout); if(k != Output->nPoints) return xf_FILE_READ_ERROR;
      fread(&k, sizeof(int), 1, fout); if(k != Output->nVars)   return xf_FILE_READ_ERROR;
      
      for(i=0; i<Output->nVars; i++)
      {
         fread(&k, sizeof(int), 1, fout);
         if(k != Output->IVars[i])
            return xf_FILE_READ_ERROR;
      }
      
      fread(&Output->accumtime, sizeof(double), 1, fout); 
      xf_printf("statistic sampling time: %lf\n", Output->accumtime);

      for(i=0; i<Output->nPoints; i++)
      {
         xpoint[0] = *(Output->PointData + i*(dim + Output->nVars) + 0); 
         xpoint[1] = *(Output->PointData + i*(dim + Output->nVars) + 1); 
         xpoint[2] = *(Output->PointData + i*(dim + Output->nVars) + 2); 
      
         fread(chk, sizeof(double), 3, fout); 

         if(xpoint[0]!=chk[0] || xpoint[1]!=chk[1] || xpoint[2]!=chk[2])
            return xf_FILE_READ_ERROR;
      }
      
      //everything goes; let's load the data from restart file
      for(i=0; i<Output->nPoints; i++)
         fread(Output->PointData + i*(dim + Output->nVars) + dim, sizeof(double), Output->nVars, fout);

      fclose(fout);

      xf_printf("Data sample continued ...\n");
   }
            
   //destroy EES struct
   ierr = xf_Error(xf_DestroyElemSearchStructure(ESS));
   if (ierr != xf_OK) return ierr;

   xf_Release((void *) gridx);

   return xf_OK;
}

/******************************************************************/
//pointwise calculation for output
//dimension of s = QutputQuantityActive[0];
static int
Yu_OutputScalarRaw(Yu_Model *Model, Yu_Output *Output, const int nq, const real *qU,
                const real *qgU, real *s)
{
   int i, j, k, sr, dim, ierr;
   real rho, u, v, w, p, vort;
   real *U, gU[50*3], gamma;
   
   sr = Model->nVars;
   dim = Model->dim;
   
   if(Model->GammaVaryFlag) return xf_NOT_SUPPORTED;
   gamma = Model->GammaInit;
   
   for (i=0; i<nq; i++)
   {
      U = qU + sr*i;
      if(qgU != NULL)
      for(j=0; j<dim; j++)
         for(k=0; k<sr; k++)
            gU[j*sr+k] = qgU[j*nq*sr + i*sr + k];
      
      rho = U[0];
      u = U[1]/rho; v = U[2]/rho;
      if(dim == 3)
         w = U[3]/rho;
      else
         w = 0.0;
      p = (gamma-1.)*(U[dim+1] - 0.5*rho*(u*u + v*v + w*w));
     
      //for testing
      //s[i] = 0.5*U[2]*U[2]/U[0];
      for(j=0; j<Output->nVars; j++)
         //switch(OutputQuantityActive[j]){
         switch(Output->IVars[j]){
            case rho_bar:
               s[i*Output->nVars + j] = rho;
               break;
           
            case rhou:
               s[i*Output->nVars + j] = U[1];
               break;

            case rhov:
               s[i*Output->nVars + j] = U[2];
               break;
            
            case rhow:
               s[i*Output->nVars + j] = U[3];
               break;

            case rhoE:
               s[i*Output->nVars + j] = U[dim+1];
               break;

            case rho_rms:
               s[i*Output->nVars + j] = (rho-1.)*(rho-1.);
               break;
            
            case u_bar:
               s[i*Output->nVars + j] = u;
               break;
               
            case u_rms:
               s[i*Output->nVars + j] = u*u;
               break;
              
            case v_bar:
               s[i*Output->nVars + j] = v;
               break;
               
            case v_rms:
               s[i*Output->nVars + j] = v*v;
               break;
            
            case w_bar:
               s[i*Output->nVars + j] = w;
               break;
               
            case w_rms:
               s[i*Output->nVars + j] = w*w;
               break;
          
            case RS_uv:
               s[i*Output->nVars + j] = u*v;
               break;
            
            case RS_uw:
               s[i*Output->nVars + j] = u*w;
               break;
            
            case RS_vw:
               s[i*Output->nVars + j] = v*w;
               break;

            case p_mean:
               s[i*Output->nVars + j] = p;
               break;

            case kinenergy:
               s[i*Output->nVars + j] = 0.5*rho*(u*u+v*v+w*w);
               break;

            case enstrophy:
              // ierr = xf_Error(Compute_Vorticity_Mag(sr, dim, U, gU, &vort)); 
              // if (ierr != xf_OK) return ierr;
              // s[i*Output->nVars + j] = 0.5 * pow(vort, 2.0); //def: Carton, et al Int. J Numer. Meth. Fluids.
              ierr = xf_Error(Compute_Divergence_Mag(sr, dim, U, gU, &vort));
              if (ierr != xf_OK) return ierr;
              //s[i*Output->nVars + j] = vort*vort;
              s[i*Output->nVars + j] = fabs(vort);
              break;

            default:
               return xf_Error(xf_NOT_SUPPORTED);
               break;
         }
      
   }
   
   return xf_OK;
}

/******************************************************************/
static int
Yu_userdefine_OutputScalarRaw(Yu_Model *Model, Yu_Output *Output, const int nq, const real *qU,
                   const real *qgU, real *s)
{
   int i, j, k, sr, dim, ierr;
   real rho, u, v, w, p, vort;
   real *U, gU[50*3], gamma;
   
   sr = Model->nVars;
   dim = Model->dim;
   
   if(Model->GammaVaryFlag) return xf_NOT_SUPPORTED;
   gamma = Model->GammaInit;
   
   for (i=0; i<nq; i++)
   {
      U = qU + sr*i;
      if(qgU != NULL)
         for(j=0; j<dim; j++)
            for(k=0; k<sr; k++)
               gU[j*sr+k] = qgU[j*nq*sr + i*sr + k];
      
      rho = U[0];
      u = U[1]/rho; v = U[2]/rho;
      if(dim == 3)
         w = U[3]/rho;
      else
         w = 0.0;
      p = (gamma-1.)*(U[dim+1] - 0.5*rho*(u*u + v*v + w*w));
      
      //for testing
      //s[i] = 0.5*U[2]*U[2]/U[0];
      for(j=0; j<Output->nVars; j++)
         //switch(OutputQuantityActive[j]){
         switch(Output->IVars[j]){
            case rho_bar:
               s[i*Output->nVars + j] = rho;
               break;
               
            case rho_rms:
               s[i*Output->nVars + j] = (rho-1.) * (rho-1.);
               break;
            
            case rhou:
               s[i*Output->nVars + j] = U[1];
               break;

            case rhov:
               s[i*Output->nVars + j] = U[2];
               break;
            
            case rhow:
               s[i*Output->nVars + j] = U[3];
               break;

            case rhoE:
               s[i*Output->nVars + j] = U[dim+1];
               break;

               
            case u_bar:
               s[i*Output->nVars + j] = u;
               break;
               
            case u_rms:
               s[i*Output->nVars + j] = u*u;
               break;
               
            case v_bar:
               s[i*Output->nVars + j] = v;
               break;
               
            case v_rms:
               s[i*Output->nVars + j] = v*v;
               break;
               
            case w_bar:
               s[i*Output->nVars + j] = w;
               break;
               
            case w_rms:
               s[i*Output->nVars + j] = w*w;
               break;
               
            case RS_uv:
               s[i*Output->nVars + j] = u*v;
               break;
               
            case RS_uw:
               s[i*Output->nVars + j] = u*w;
               break;
               
            case RS_vw:
               s[i*Output->nVars + j] = v*w;
               break;
               
            case p_mean:
               s[i*Output->nVars + j] = p;
               break;
               
            case kinenergy:
               s[i*Output->nVars + j] = 0.5*rho*(u*u+v*v+w*w);
               break;
               
            case enstrophy:
              // ierr = xf_Error(Compute_Vorticity_Mag(sr, dim, U, gU, &vort));
              // if (ierr != xf_OK) return ierr;
              // s[i*Output->nVars + j] = 0.5 * pow(vort, 2.0); //def: Carton, et al Int. J Numer. Meth. Fluids.
               ierr = xf_Error(Compute_Divergence_Mag(sr, dim, U, gU, &vort));
               if (ierr != xf_OK) return ierr;
               //s[i*Output->nVars + j] = vort * vort; 
               s[i*Output->nVars + j] = fabs(vort);
               break;
               
            default:
               return xf_Error(xf_NOT_SUPPORTED);
               break;
         }
      
   }
   
   return xf_OK;
}
/******************************************************************/
//raw data for statistics should be processed before dump for visualization
static int
Yu_OutputScalarMature(char *Name, const int nVars, int *IVars, const real *data, const real time, real *out)
{
   int i, j, k; 
   real tmp;

   if(strcmp(Name, "StatPlane") == 0 || strcmp(Name, "data_dump") == 0)
   {
      for(j=0; j<nVars; j++)
         switch(IVars[j]){
            case rho_bar:
            case u_bar:
            case v_bar:
            case w_bar:
            case p_mean:
               out[j] = data[j]/time;
               break;

            case rho_rms:
            case u_rms:
            case v_rms:
            case w_rms:
            case rhou:
            case rhov:
            case rhoE:
               out[j] = data[j]/time  - (data[j-1]/time)*(data[j-1]/time);
               //out[j] = data[j]/time;
               break;

            case RS_uv:   //need to make sure u_bar and v_bar are sampled along
               tmp = 1.0;
               for(k=0, i=0; k<nVars; k++)
                  if(IVars[k]==u_bar || IVars[k]==v_bar)
                  {  tmp *= data[k]/time; i++; }
               if(i != 2) return xf_CODE_LOGIC_ERROR;
               out[j] = data[j]/time - tmp;
              
               //temporary change for channel flow
               //out[j] = data[j]/time;
               break;

            case RS_uw:   //need to make sure u_bar and v_bar are sampled along
               tmp = 1.0;
               for(k=0, i=0; k<nVars; k++)
                  if(IVars[k]==u_bar || IVars[k]==w_bar)
                  {  tmp *= data[k]/time; i++; }
               if(i != 2) return xf_CODE_LOGIC_ERROR;
               out[j] = data[j]/time - tmp;
               break;

            case RS_vw:   //need to make sure u_bar and v_bar are sampled along
               tmp = 1.0;
               for(k=0, i=0; k<nVars; k++)
                  if(IVars[k]==v_bar || IVars[k]==w_bar)
                  {  tmp *= data[k]/time; i++; }
               if(i != 2) return xf_CODE_LOGIC_ERROR;
               out[j] = data[j]/time - tmp;
               break;
            
            default:
               return xf_Error(xf_NOT_SUPPORTED);
               break;
      }
   }

   if(strcmp(Name, "Padapt_Test") == 0)
   {
      for(j=0; j<nVars; j++)
         switch(IVars[j]){
            case rho_bar:
            case u_bar:
            case v_bar:
            case w_bar:
            case p_mean:
               out[j] = data[j]/time;
               break;
            
            default:
               return xf_Error(xf_NOT_SUPPORTED);
               break;
         }
   }
   
   if(strcmp(Name, "acoustic_ML_Stat") == 0)
   {
      for(j=0; j<nVars; j++)
         switch(IVars[j]){
            case rho_bar:
            case u_bar:
            case v_bar:
               out[j] = data[j]/time;
               break;
               
            case rho_rms:
            case u_rms:
            case v_rms:
               out[j] = data[j]/time  - (data[j-1]/time)*(data[j-1]/time);
               break;
               
            case RS_uv:   //need to make sure u_bar and v_bar are sampled along
               tmp = 1.0;
               for(k=0, i=0; k<nVars; k++)
                  if(IVars[k]==u_bar || IVars[k]==v_bar)
                  {  tmp *= data[k]/time; i++; }
               if(i != 2) return xf_CODE_LOGIC_ERROR;
               out[j] = data[j]/time - tmp;
               break;
            
            default:
               return xf_Error(xf_NOT_SUPPORTED);
               break;
         }
   }
               
   return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_OutputPointValue
//half borrow from XFlow; half extension done by Yu
//for shock simulation, several probes can be placed using this routine.
static int
xf_OutputPointValue( xf_All *All, Yu_Output *Output, const xf_Vector *U,
                    real weight, real *Value, xf_Vector *Value_U,
                    enum xfe_AddType AddFlag)
{
   int ierr, k, i, j, n, sr, egrp, elem;
   int Order, nn, dim, myRank, nProc;
   int *IParam;
   real Gamma, dt;
   enum xfe_BasisType Basis;
   enum xfe_Bool MotionOn = xfe_False;
   real *RParam, *xq, *wq, *u, *s, *s_u;
   real *xref, *EU, *EV_U, *data;
   real Time, xglob[3], dp;
   xf_BasisData *PhiData;
   xf_Vector *J_GCL, *GammaVec;
   xf_Data *GammaDat;
   //xf_EqnSet *EqnSet;
   Yu_Model *Model;
   xf_Mesh *Mesh;
   xf_MotionData *MD = NULL;
   
   if (Output->UsesFlux) return xf_Error(xf_INPUT_ERROR);
   
   //ierr = xf_FindDataByTitle(All->DataSet, "HeatCapacityRatio", xfe_Vector, &GammaDat);
   //if(ierr == xf_NOT_FOUND)
   //{
   //   xf_printf("Cannot find heat capacity ratio...\n");
   //   return ierr;
   //}
   //else
   //   GammaVec = (xf_Vector *) GammaDat->Data;
   
   //EqnSet = All->EqnSet;
   Model  = All->Model;
   //sr     = EqnSet->StateRank;
   sr     = Model->nVars;
   Mesh   = All->Mesh;
   dim    = Mesh->Dim;
   
   //if (Value_U != NULL){
      // zero out Value_U if requesting a set or neg
   //   if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
   //      ierr = xf_Error(xf_SetZeroVector(Value_U));
   //      if (ierr != xf_OK) return ierr;
   //   }
   //}
   
   // build eqnset-desired parameter lists for passing into functions
   //ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
   //if (ierr != xf_OK) return ierr;
   
   // determine if we need mesh motion
   //MotionOn = ((Mesh->Motion != NULL) && (Mesh->Motion->Active));
   //if (MotionOn){
   //   ierr = xf_Error(xf_CreateMotionData(All, &MD));
   //   if (ierr != xf_OK) return ierr;
   //}
   
   // locate J_GCL (zero out if AddFlag is set/neg)
   //if (Value_U != NULL){ // only if also calculating Value_U
   //   ierr = xf_Error(xf_FindOutputGCLLinearization(All, Output->Name, AddFlag, &J_GCL));
   //   if (ierr != xf_OK) return ierr;
   //}
   //else J_GCL = NULL;
   
   //will be updated after output is processed
   //Output->pretime = Time;
   
   // initialize vars to NULL
   PhiData = NULL;
   u       = NULL;
   s       = NULL;
   s_u     = NULL;
   
   /* Calculate element info.  May require a global-to-local map, but
    this is done just once. */ 
   if (Output->elemLocal == NULL){

      ierr = xf_Error(Yu_OutputPointValueInit(All, Output));
      if (ierr != xf_OK) return ierr;

      xf_printf("Output point location info is found!\n");
      // ierr = xf_Error(xf_Alloc( (void **) &Output->egrpLocal, Output->nPoints, sizeof(int)));
      //if (ierr != xf_OK) return ierr;
      //ierr = xf_Error(xf_Alloc( (void **) &Output->elemLocal, Output->nPoints, sizeof(int)));
      //if (ierr != xf_OK) return ierr;
      
      //for(i=0; i<Output->nPoints; i++){
      //   ierr = xf_Error(xf_GetLocalElem(All->Mesh, Output->egrp[i], Output->elem[i],
      //                                   Output->egrpLocal + i, Output->elemLocal + i));
      //   if (ierr != xf_OK) return ierr;
      //}
   }
   
   // determine Time & Time step
   ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
   if (ierr != xf_OK) return ierr;
   
   dt = Time - Output->pretime;
   Output->accumtime += dt;
   xf_printf("Output:%s pointwise data sample at %lf\n", Output->Name, Time);
   
   // allocate data
   ierr = xf_Error(xf_ReAlloc( (void **)  &u, sr, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_ReAlloc( (void **)  &s, Output->nVars, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   if (Value_U != NULL){
      ierr = xf_Error(xf_ReAlloc( (void **)  &s_u, sr, sizeof(real)));
      if (ierr != xf_OK) return ierr;
   }
  
   // element information
   for(i=0; i<Output->nPoints; i++){
   
      egrp = Output->egrpLocal[i];
      elem = Output->elemLocal[i];
      //imp: always use 3D coordinates
      xref = Output->xref + 3*i;
  
      if (egrp < 0){ // this is the case in parallel if loc is not on proc
         if (Value != NULL) (*Value) = 0.0;
     
         //jump over the coordinates 
         for(j=3; j<Output->nVars+3; j++)
            Output->PointData[i*(3+Output->nVars) + j] = 0.0; 
       
         //directly jump to the next sample point
         continue;
      }
   
   // need global coordinate of point if motion is on
   //if (MotionOn){
   //   ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, NULL, xfe_True, 1, xref, xglob));
   //   if (ierr != xf_OK) return ierr;
   //}
  
      // basis and order of state vector on egrp
      Basis = U->Basis[egrp];
      Order = xf_InterpOrder(U, egrp, elem);
   
   // compute basis functions at xref
      ierr = xf_Error(xf_EvalBasis(Basis, Order, xfe_True, 1, xref, xfb_Phi, &PhiData));
      if (ierr != xf_OK) return ierr;
   
      nn = PhiData->nn; // number of interpolation nodes
   
   // obtain transformation map if doing mesh motion
   //if (MotionOn){
   //   ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, Mesh->Motion,
   //                                    1, dim, Time, xglob, MD));
   //   if (ierr != xf_OK) return ierr;
   //}
   
      EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
   //Gamma = GammaVec->GenArray[egrp].rValue[elem][0];
   
   // interpolate state at quad points
      xf_MxM_Set(PhiData->Phi, EU, 1, nn, sr, u);
   
   // transform state to physical
   //if (MotionOn) xf_ModMotionPreEqnCall(1, dim, sr, MD, u, NULL);
   
   // call eqnset specific function for scalar
      ierr = xf_Error(Yu_OutputScalarRaw(Model, Output, 1, u, NULL, s));
      if (ierr != xf_OK) return ierr;
      
      //ierr = xf_Error(xf_EqnSetScalar(EqnSet, Output->ScalarName, IParam,
      //                             RParam, 1, u, NULL, s, s_u, NULL, NULL, Gamma));
      //if (ierr != xf_OK) return ierr;
   
      //update the output data
      if(!Output->SequanceDump)
         for(j=0; j<Output->nVars; j++)
            Output->PointData[i*(3+Output->nVars) + 3 + j] += dt*s[j];
      else
         for(j=0; j<Output->nVars; j++)
            Output->PointData[i*(3+Output->nVars) + 3 + j] = s[j];
   // Divide output linearization by gbar to make it wrt reference state
   //if ((MotionOn) && (Value_U != NULL) && (s_u != NULL))
   //   xf_ColDiv(s_u, MD->gb, 1, sr, 1);
   
      if (Value != NULL){
      //if (Output->DomainNorm == xfe_DomainNormL2){
      //   (*Value) = s[0]*s[0]*weight; //square of scalar
      //}
      //else if (Output->DomainNorm == xfe_DomainNormNone)
         ////////replaced/////////
         //(*Value) = s[0]*weight;
         (*Value) = 0.0;
         
      // set Value :: include weight here
      }
      
   }//i
   // Value_U{n,k} @= s_u{1,k} * Phi{1,n}
   //if (Value_U != NULL){
   //   EV_U = Value_U->GenArray[egrp].rValue[elem]; // Value_U on elem [nn*sr]
   //   for (k=0; k<sr; k++) s_u[k] *= weight; // :: include weight here
   //   xf_MTxM(PhiData->Phi, s_u, nn, 1, sr, AddFlag, EV_U);
   //   if (J_GCL != NULL){
         // J_GCL{n} = J_ubar{k} dot (-u{k}) * Phi{1,n}
         // note, u{k} is the physical state
   //      xf_DotProduct(s_u, u, sr, &dp);
   //      xf_cV_Add(PhiData->Phi, -dp, nn, AddFlag, J_GCL->GenArray[egrp].rValue[elem]);
   //   }
   //}
   
   /* Destroy Basis Data */
   ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
   if (ierr != xf_OK) return ierr;
   
   /* Destroy mesh motion data */
   xf_DestroyMotionData(MD);
  
   //xf_Release( (void *) IParam);
   //xf_Release( (void *) RParam);
   xf_Release( (void *) u);
   xf_Release( (void *) s);
   xf_Release( (void *) s_u);
   
   return xf_OK;
}




/******************************************************************/
//   FUNCTION Definition: xf_OutputDomainIntegral
static int
xf_OutputDomainIntegral( xf_All *All, Yu_Output *Output, const xf_Vector *U,
                        real weight, real *Value, xf_Vector *Value_U,
                        enum xfe_AddType AddFlag)
{
   int ierr, i, k, sr, sr2, iq, nq, pnq, dim;
   int Order, QuadOrder, QuadOrder0, IntOrder, nn;
   int egrp, elem, nelemtot, myRank, nProc;
   int *IParam;
   char buf[xf_MAXSTRLEN];
   enum xfe_BasisType Basis;
   enum xfe_Bool found, QuadChanged, Need_WD, DoOffset = xfe_False;
   enum xfe_Bool MotionOn = xfe_False;
   real *RParam, *xq, *wq, *u, *gu, *s, *s_u;
   real *EU, *EV_U, *EV_G, LocValue, Time;
   real *xglob, *f, *wd, wdf, offset;
   real DiscreteL2Sum, vol, Volume, ElemLocValue, dp;
   xf_QuadData *QuadData;
   xf_BasisData *PhiData, *GeomPhiData, *WDPhiData;
   xf_Vector *Vtemp, *WD;
   xf_Vector *J_GCL = NULL;
   xf_JacobianData *JData;
   xf_Mesh *Mesh;
   //xf_EqnSet *EqnSet;
   Yu_Model *Model;
   xf_MotionData *MD = NULL;
   FILE *fout;
   
   Mesh = All->Mesh;
   dim  = Mesh->Dim;
   
   //EqnSet = All->EqnSet;
   Model  = All->Model;
   //sr     = EqnSet->StateRank;
   sr     = Model->nVars;
   sr2    = sr*sr;
   
   // build eqnset-desired parameter lists for passing into functions
   //ierr = xf_Error(xf_RetrieveFcnParams(All, EqnSet, &IParam, &RParam, NULL, NULL));
   //if (ierr != xf_OK) return ierr;
  
   //reset
   if(Value != NULL)
      for(i=0; i<Output->nVars; i++)
         Value[i] = 0.0;

   // set LocValue to 0
   LocValue = 0.0;
   
   // used for discrete L2 error norms
   DiscreteL2Sum = 0;
   
   // volume for PerVol L2 norm
   Volume = 0.;
   
   // need a temporary vector if doing an L2 norm
   //if ( ((Output->DomainNorm == xfe_DomainNormL2) ||
   //      (Output->DomainNorm == xfe_DomainNormL2PerVol))
   //    && (Value_U != NULL)){
   //   ierr = xf_Error(xf_FindSimilarVector(All, Value_U, "Vtemp", xfe_False,
   //                                        xfe_True, NULL, &Vtemp, NULL));
   //   if (ierr != xf_OK) return ierr;
   //   ierr = xf_Error(xf_SetZeroVector(Vtemp));
   //   if (ierr != xf_OK) return ierr;
   //}
   
   // linearization should not be requested on an error norm output ...
   // not impossible to implement, but cannot think of when would need it
   //if ((Output->DomainNorm == xfe_DomainNormL2Error) && (Value_U != NULL))
   //   return xf_Error(xf_NOT_SUPPORTED);
   
   // linearization not implemented for a discrete L2 norm
   //if ((Output->DomainNorm == xfe_DomainNormL2Discrete) && (Value_U != NULL))
   //   return xf_Error(xf_NOT_SUPPORTED);
   
   // check if need a wall distance or offset
   Need_WD  = xfe_False;
   DoOffset = xfe_False;
   //if ((Output->Function != NULL) && (strncmp(Output->Function, "WallDistance", 12) == 0)){
      /*
       The wall distance is used to introduce a multiplicative factor
       into an interior domain integral.  This factor is exp(-wdf*wd^2),
       where wdf is the wall distance decay factor, prescribed in
       Output->Data.
       */
   //   if ((Output->Data != NULL) && (sscanf(Output->Data, "%lf", &wdf) == 1)){
   //      ierr = xf_Error(xf_FindSupportedVector(All, "WallDistance", &WD));
   //      if (ierr != xf_OK) return ierr;
   //      Need_WD = xfe_True;
   //   }
   //   else
   //      xf_printf("Warning: wall distance requested but decay factor not specified in Data.\n");
   //}
   //else{
      // a value in Output->Data indicates a desired offset in scalar value
   //   if ((Output->Data != NULL) && (sscanf(Output->Data, "%lf", &offset) == 1))
   //      DoOffset = xfe_True;
   //}
   
   // determine if we need mesh motion
   //MotionOn = ((Mesh->Motion != NULL) && (Mesh->Motion->Active));
   //if (MotionOn){
   //   ierr = xf_Error(xf_CreateMotionData(All, &MD));
   //   if (ierr != xf_OK) return ierr;
   //}
   
   // locate J_GCL (zero out if AddFlag is set/neg)
   //if (Value_U != NULL){ // only if also calculating Value_U
   //   ierr = xf_Error(xf_FindOutputGCLLinearization(All, Output->Name, AddFlag, &J_GCL));
   //   if (ierr != xf_OK) return ierr;
   //}
   //else J_GCL = NULL;
   
   // determine Time
   ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
   if (ierr != xf_OK) return ierr;
   
   // initialize vars to NULL
   QuadData    = NULL;
   PhiData     = NULL;
   JData       = NULL;
   u           = NULL;
   gu          = NULL;
   wq          = NULL;
   s           = NULL;
   s_u         = NULL;
   f           = NULL;
   xglob       = NULL;
   GeomPhiData = NULL;
   WDPhiData   = NULL; // for wall distance basis
   wd          = NULL; // for wall distance values
   
   pnq         = -1;  // previous number of quad points
   
   // loop over element groups
   for (egrp=0; egrp<Mesh->nElemGroup; egrp++){
      // Determine Basis and Order from the state, U
      Basis = U->Basis[egrp];
      
      ElemLocValue = LocValue;
      
      // loop over elements
      for (elem=0; elem<Mesh->ElemGroup[egrp].nElem; elem++){
         
         // get interpolation order
         Order = xf_InterpOrder(U, egrp, elem);
         
         // determine required integration order using eqn-set specific rule
         IntOrder = Order;
         ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, IntOrder, &QuadOrder));
         if (ierr != xf_OK) return ierr;
         
         /* need higher order integration if doing L2, since integrand is
          squared but do not account for eqn-set here; could lose
          accuracy if calculating the L2 norm of a highly-nonlinear
          function of the state.  If we did account for eqn-set, we
          would run into max quad order problems in some cases. */
        //might be important later; notes by Yu Lv
        // if ((Output->DomainNorm == xfe_DomainNormL2) ||
        //     (Output->DomainNorm == xfe_DomainNormL2PerVol) ||
        //     (Output->DomainNorm == xfe_DomainNormL2Error)){
        //    ierr = xf_Error(xf_GetQuadOrderElem(Mesh, NULL, egrp, 2*IntOrder, &QuadOrder0));
        //    if (ierr != xf_OK) return ierr;
        //    QuadOrder = max(QuadOrder0, QuadOrder);
        // }
         
         /* Pull off quad points for the element; will not recalculate if
          Basis/Order have not changed. */
         ierr = xf_Error(xf_QuadElem(Mesh, egrp, elem, QuadOrder, &QuadData, &QuadChanged));
         if (ierr != xf_OK) return ierr;
         
         nq = QuadData->nquad;
         xq = QuadData->xquad;
         
         // compute basis functions (and grads) if quad or basis or order changed
         ierr = xf_Error(xf_EvalBasis(Basis, Order, QuadChanged, nq, xq, xfb_Phi | xfb_GPhi | xfb_gPhi, &PhiData));
         if (ierr != xf_OK) return ierr;
         
         /* Compute geometry Jacobian; if not constant, compute at quad
          points.  Note if jacobian is constant, only one Jacobian will
          be computed/returned. */
         ierr = xf_Error(xf_ElemJacobian(Mesh, egrp, elem, nq, xq, xfb_detJ | xfb_iJ, QuadChanged, &JData));
         if (ierr != xf_OK) return ierr;
         
         nn = PhiData->nn; // number of interpolation nodes
         
         // re-allocate data if quad points increased
         if (nq > pnq){
            ierr = xf_Error(xf_ReAlloc( (void **)  &u, nq*sr, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_ReAlloc( (void **)  &gu, dim*nq*sr, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_ReAlloc( (void **) &wq, nq, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            //multiple output quantities
            ierr = xf_Error(xf_ReAlloc( (void **)  &s, nq*Output->nVars, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            ierr = xf_Error(xf_ReAlloc( (void **) &xglob, nq*dim, sizeof(real)));
            if (ierr != xf_OK) return ierr;
            if (Value_U != NULL){
               ierr = xf_Error(xf_ReAlloc( (void **)  &s_u, nq*sr, sizeof(real)));
               if (ierr != xf_OK) return ierr;
            }
            //if (Output->DomainNorm == xfe_DomainNormL2Error){
            //   ierr = xf_Error(xf_ReAlloc( (void **)  &f, nq*sr, sizeof(real)));
            //   if (ierr != xf_OK) return ierr;
            //}
         }
         
         // obtain global coords of quad points
         ierr = xf_Error(xf_Ref2GlobElem(Mesh, egrp, elem, &GeomPhiData, QuadChanged,
                                         nq, xq, xglob));
         if (ierr != xf_OK) return ierr;
         
         // obtain transformation map if doing mesh motion
         //if (MotionOn){
         //   ierr = xf_Error(xf_MeshMotionMap( egrp, elem, PhiData, Mesh->Motion,
         //                                    nq, dim, Time, xglob, MD));
         //   if (ierr != xf_OK) return ierr;
         //}
         
         EU = U->GenArray[egrp].rValue[elem]; // U on elem [nn*sr]
         
         // interpolate state at quad points
         xf_MxM_Set(PhiData->Phi, EU, nq, nn, sr, u);
         
         // transform state to physical
         //if (MotionOn) xf_ModMotionPreEqnCall(nq, dim, sr, MD, u, NULL);
         
         /* Compute geometry Jacobian */
         //ierr = xf_Error(xf_ElemJacobian(All->Mesh, egrp, elem, nq, xq,
         //                                xfb_detJ | xfb_iJ | xfb_J, QuadChanged, &JData));
         //if (ierr != xf_OK) return ierr;
         
         /* convert reference basis grads (GPhi) to physical grads, gPhi */
         ierr = xf_Error(xf_EvalPhysicalGrad(PhiData, JData));
         if (ierr != xf_OK) return ierr;
         
         // compute state gradient
         for (i=0; i<dim; i++)
            xf_MxM_Set(PhiData->gPhi+nn*nq*i, EU, nq, nn, sr, gu + nq*sr*i);
         
         // transform state gradient to physical
         //if (MotionOn) xf_ModMotionPhysGrad(nq, dim, sr, MD, u, gu);
         
         // form detJ-multiplied quad weight vector, wq, multiplied by weight
         for (iq=0; iq<nq; iq++)
            wq[iq] = QuadData->wquad[iq]*JData->detJ[iq*(JData->nq!=1)]*weight;
         
         // account for mesh motion in quad weights (elements are deformed in physical space)
         //if (MotionOn) xf_ColMult(wq, MD->g, nq, 1, 1);
         
         // interpolate wall distance if necessary and adjust quad weights
         //if (Need_WD){
         //   if (nq > pnq){
         //      ierr = xf_Error(xf_ReAlloc( (void **) &wd, nq*WD->StateRank, sizeof(real)));
         //      if (ierr != xf_OK) return ierr;
         //   }
         //   ierr = xf_Error(xf_EvalBasis(WD->Basis[egrp], WD->Order[egrp], QuadChanged,
         //                                nq, xq, xfb_Phi, &WDPhiData));
         //   if (ierr != xf_OK) return ierr;
         //   xf_MxM_Set(WDPhiData->Phi, WD->GenArray[egrp].rValue[elem], nq,
         //              WDPhiData->nn, WD->StateRank, wd);
            
            // quad weights are adjusted here to incorporate the wall distance
         //   for (iq=0; iq<nq; iq++) wq[iq] *= exp(-wdf*wd[iq]*wd[iq]);
         //}
         
         
         /** Distinguish between state error and scalar norms **/
         
         //if (Output->DomainNorm == xfe_DomainNormL2Error){
            
            // call eqnset specific function for state values
         //   ierr = xf_Error(xf_EqnSetFcnState(EqnSet, Output->Function, Output->Data, IParam,
         //                                     RParam, nq, (MotionOn) ? MD->x : xglob, &Time, f));
         //   if (ierr != xf_OK) return ierr;
            
            // set s = square error at each quad point
         //   for (iq=0; iq<nq; iq++){
         //      for (k=0, s[iq]=0.; k<sr; k++){
         //         i = iq*sr+k;
         //         s[iq] += (u[i]-f[i])*(u[i]-f[i]);
         //      }
         //   }
            
            // linearization is possible, but currently see no use for it
         //   if (Value_U != NULL) return xf_Error(xf_NOT_SUPPORTED);
            
         //}
         //else
         {
            
            // call eqnset specific function for scalar
            
            //default routine for volume intergral quantities; 
            //ierr = xf_Error(Yu_OutputScalarRaw(Model, Output, nq, u, gu, s));
            //if (ierr != xf_OK) return ierr;

            //user defined but index of quantities are borrow from default setting
            ierr = xf_Error(Yu_userdefine_OutputScalarRaw(Model, Output, nq, u, gu, s));
            if (ierr != xf_OK) return ierr;

            
            // Divide output linearization by gbar to make it wrt reference state
            //if ((MotionOn) && (Value_U != NULL) && (s_u != NULL))
            //   xf_ColDiv(s_u, MD->gb, nq, sr, 1);
            
            // subtract an offset if specified
            //if (DoOffset) for (iq=0; iq<nq; iq++) s[iq] -= offset;
            
            // subtract linearization from adjoint residual if not null
            //if (Value_U != NULL){
            //   if ((Output->DomainNorm == xfe_DomainNormL2) || (Output->DomainNorm == xfe_DomainNormL2PerVol))
            //      EV_U = Vtemp->GenArray[egrp].rValue[elem];
            //   else
            //      EV_U = Value_U->GenArray[egrp].rValue[elem]; // Value_U on elem [nn*sr]
               
               // multiply s_u by wq
            //   xf_ColMult(s_u, wq, nq, sr, 1);
               
               /*
                Need to set, add-to, subtract-from, etc. Value_U:
                
                No norm: Value   = int(s   dx)
                Value_U = int(s_U dx)
                
                L2 Norm: Value   = sqrt(int(s^2 dx))
                Value_U = int(s*s_U dx) * 1/Value
                
                On each elem, s_U{n,k} = s_u{q,k} * Phi{q,n}
                */
               
             //  if ((Output->DomainNorm == xfe_DomainNormL2) || (Output->DomainNorm == xfe_DomainNormL2PerVol))
             //     xf_ColMult(s_u, s, nq, sr, 1); // multiply s_u by s
             //  else if (Output->DomainNorm != xfe_DomainNormNone)
             //     return xf_Error(xf_NOT_SUPPORTED);
               
               // Value_U{n,k} @= int(s_u{q,k} * Phi{q,n})
             //  xf_MTxM(PhiData->Phi, s_u, nn, nq, sr, AddFlag, EV_U);
               
               // Value_G{n} @= int(s_u{q,k} * (-u{q,k}) * Phi{q,n}) -- u is physical here
             //  if (J_GCL != NULL){
             //     EV_G = J_GCL->GenArray[egrp].rValue[elem]; // J_GCL on elem [nn*sr]
             //     for (iq=0; iq<nq; iq++){
             //        xf_DotProduct(s_u+iq*sr, u+iq*sr, sr, &dp);
             //        s_u[iq] = -dp; // s_u is overwritten here
             //     }
             //     xf_MTxM(PhiData->Phi, s_u, nn, nq, 1, AddFlag, EV_G);
             //  }
               
               
            //}
            
            // modify scalar depending on norm
            //if (Output->DomainNorm == xfe_DomainNormL1)
            //   for (iq=0; iq<nq; iq++) s[iq] = fabs(s[iq]);
            //else if ((Output->DomainNorm == xfe_DomainNormL2) || (Output->DomainNorm == xfe_DomainNormL2PerVol))
            //   for (iq=0; iq<nq; iq++) s[iq] *= s[iq];
            // else do nothing (L2Discrete norm just uses s[iq])
            
         }
         
         // sum scalar*wq over quad points, add to Value
         //for (iq=0; iq<nq; iq++) LocValue += wq[iq]*s[iq];
         for(iq=0; iq<nq; iq++)
            for(k=0; k<Output->nVars; k++)
               Value[k] += wq[iq] * s[iq*Output->nVars+k];

         // sum quad weights to get volume
         for (iq=0, vol=0; iq<nq; iq++) vol += wq[iq];
         Volume += vol; // increment running total
         
         // discrete L2 norm requires element averages
         //if (Output->DomainNorm == xfe_DomainNormL2Discrete){
         //   ElemLocValue = LocValue - ElemLocValue;
         //   ElemLocValue /= vol;
         //   DiscreteL2Sum += ElemLocValue*ElemLocValue; // sum of square averages
         //}
         
         pnq = nq;
      } // elem
      
   } // egrp
   
   //if (Output->DomainNorm == xfe_DomainNormL2Discrete){
      // special sum/reduce for discrete L2 norm
   //   ierr = xf_Error(xf_MPI_Allreduce(&DiscreteL2Sum, 1, xfe_SizeReal, xfe_MPI_SUM));
   //   if (ierr != xf_OK) return ierr;
      
   //   nelemtot = Mesh->ElemGroup[egrp].nElem;
   //   ierr = xf_Error(xf_MPI_Allreduce(&nelemtot, 1, xfe_SizeInt, xfe_MPI_SUM));
   //   if (ierr != xf_OK) return ierr;
      
   //   LocValue = DiscreteL2Sum / ( (real) nelemtot);
   //}
   //else
   {
      // sum reduce LocValue
      //ierr = xf_Error(xf_MPI_Allreduce(&LocValue, 1, xfe_SizeReal, xfe_MPI_SUM));
      ierr = xf_Error(xf_MPI_Allreduce(Value, Output->nVars, xfe_SizeReal, xfe_MPI_SUM));
      if (ierr != xf_OK) return ierr;
         xf_printf("%lf\n", Value[3]);
   }
   
   // take sqrt of Value if L2 norm
   //if ((Output->DomainNorm == xfe_DomainNormL2) ||
   //    (Output->DomainNorm == xfe_DomainNormL2Error) ||
   //    (Output->DomainNorm == xfe_DomainNormL2PerVol) ||
   //    (Output->DomainNorm == xfe_DomainNormL2Discrete)){
      
   //   if (Volume == 0.) return xf_Error(xf_OUT_OF_BOUNDS);
   //   if (Output->DomainNorm == xfe_DomainNormL2PerVol) weight *= 1./Volume; // per volume comes in here
      
   //   LocValue = sqrt(fabs(LocValue*weight));
   //   if (weight < 0) LocValue = -LocValue;
      
      // need to multiply linearization by 1/Value (see above description)
   //   if ((weight != 0.) && (Value_U != NULL)){
   //      ierr = xf_Error(xf_VectorMultSet(Vtemp, fabs(weight/LocValue), AddFlag, Value_U));
   //      if (ierr != xf_OK) return ierr;
   //   }
   //}
   
   //if (Value != NULL) (*Value) = LocValue;
  
   ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
   if (ierr != xf_OK) return ierr;

   //dump the output to the specific file
   if(myRank == 0){
      sprintf(buf, "DomainIntg_%s.dat", Output->Name);
      fout = fopen(buf, "a");  //for appending

      //sampling time
      fprintf(fout, "%.10f ", Time);

      for(i=0; i<Output->nVars; i++)
         fprintf(fout, "%.10lf  ", Value[i]);
      fprintf(fout, "\n");
      fclose(fout);
   }
   
   // Only destroy QuadData if points are generic
   ierr = xf_Error(xf_DestroyGenericQuadData(QuadData));
   if (ierr != xf_OK) return ierr;
   
   /* Destroy Basis Data */
   ierr = xf_Error(xf_DestroyBasisData(PhiData, xfe_True));
   if (ierr != xf_OK) return ierr;
   
   /* Destroy Geometry Basis Data */
   ierr = xf_Error(xf_DestroyBasisData(GeomPhiData, xfe_True));
   if (ierr != xf_OK) return ierr;
   
   /* Destroy Wall Distance Basis Data */
   ierr = xf_Error(xf_DestroyBasisData(WDPhiData, xfe_True));
   if (ierr != xf_OK) return ierr;
   
   /* Destroy geometry Jacobian Data */
   ierr = xf_Error(xf_DestroyJacobianData(JData));
   if (ierr != xf_OK) return ierr;
   
   /* Destroy mesh motion data */
   xf_DestroyMotionData(MD);
   
   xf_Release( (void  *) IParam);
   xf_Release( (void  *) RParam);
   xf_Release( (void *) u);
   xf_Release( (void *) gu);
   xf_Release( (void *) wq);
   xf_Release( (void *) s);
   xf_Release( (void *) s_u);
   xf_Release( (void *) f);
   xf_Release( (void *) xglob);
   xf_Release( (void *) wd);
   
   return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_OutputBoundaryIntegral
static int
xf_OutputBoundaryIntegral( xf_All *All, Yu_Output *Output, const xf_Vector *U,
                          real weight, real *Value, xf_Vector *Value_U,
                          enum xfe_AddType AddFlag)
{
   int ierr, k, i, sr, nset, nBFG, *BFGs;
   int myRank, nProc;
   enum xfe_Bool NegativeFlag = xfe_False;
   char DumpFile[xf_MAXSTRLEN];
   int  *FluxMoments;
   real *FluxWeights, Time;
   xf_Vector *J_GCL = NULL;
   //xf_EqnSet *EqnSet;
   Yu_Model *Model;
   FILE *fout;
   xf_OutputEvalData OutputEval;
   
   //EqnSet = All->EqnSet;
   Model = All->Model;
   //sr = EqnSet->StateRank;
   sr = Model->nVars;
   
   // If not using a flux, can use a vector dot n, or a scalar
   if (!Output->UsesFlux){
      //return xf_Error(xf_OutputBoundaryIntegral_VectorScalar(All, Output, Output->VectorName != NULL,
      //                                                       U, weight, Value, Value_U, AddFlag));
      xf_printf("Direct integral on vector or scalar is not supported!\n");
      return xf_NOT_SUPPORTED;
   }
   
   if (Output->nFluxComponent <= 0) return xf_Error(xf_INPUT_ERROR);
   if (Output->FluxComponentWeights == NULL) return xf_Error(xf_INPUT_ERROR);
   
   // allocate a flux weight vector
   ierr = xf_Error(xf_Alloc((void **) &FluxWeights, sr, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc((void **) &FluxMoments, sr, sizeof(int)));
   if (ierr != xf_OK) return ierr;
   
   for (k=0; k<sr; k++) FluxWeights[k] =  0.;
   for (k=0; k<sr; k++) FluxMoments[k] = -1;
   
   // set weights according to request; include output weight here
   nset = 0;
   for (i=0; i<Output->nFluxComponent; i++){
      for (k=0; k<sr; k++)
        // if (strcmp(Output->FluxComponentNames[i], EqnSet->StateName[k]) == 0){
         if(k == Output->FluxComponentRanks[i]){
            FluxWeights[k] = Output->FluxComponentWeights[i]*weight;
            if (Output->FluxComponentMoments != NULL)
               FluxMoments[k] = Output->FluxComponentMoments[i];
            nset++;
         }
   }
   if (nset != Output->nFluxComponent) return xf_Error(xf_OUT_OF_BOUNDS);
   
   // determine which boundary face groups to integrate over
   //nBFG = Output->nBFG;
   nBFG = 1; //default: the integral is only conducted on one specific boundary
   ierr = xf_Error(xf_Alloc((void **) &BFGs, nBFG, sizeof(int)));
   if (ierr != xf_OK) return ierr;
   
   nset = 0;
   for (i=0; i<nBFG; i++){
      for (k=0; k<All->Mesh->nBFaceGroup; k++)
         if (strcmp(Output->BFGTitles, All->Mesh->BFaceGroup[k].Title) == 0){
            BFGs[i] = k;
            nset++;
         }
   }
   if ((nset != nBFG) || (nBFG > All->Mesh->nBFaceGroup))
      return xf_Error(xf_OUT_OF_BOUNDS);
   
   /*if (Value_U != NULL){
      // zero out Value_U if requesting a set or neg
      if ((AddFlag == xfe_Set) || (AddFlag == xfe_Neg)){
         ierr = xf_Error(xf_SetZeroVector(Value_U));
         if (ierr != xf_OK) return ierr;
      }
      // we will use a function that adds to Value_U; so set weights if want neg
      if ((AddFlag == xfe_Sub) || (AddFlag == xfe_Neg)){
         NegativeFlag = xfe_True;
         for (k=0; k<sr; k++) FluxWeights[k] *= -1.;
      }
   }*/
   
   // locate J_GCL (zero out if AddFlag is set/neg)
   //if (Value_U != NULL){ // only if also calculating Value_U
   //   ierr = xf_Error(xf_FindOutputGCLLinearization(All, Output->Name, AddFlag, &J_GCL));
   //   if (ierr != xf_OK) return ierr;
   //}
   //else
   J_GCL = NULL;
   
   // roll output evaluation inputs into a structure
   OutputEval.Value       = Value;
   OutputEval.Value_U     = Value_U;
   OutputEval.Value_G     = J_GCL;
   OutputEval.FluxWeights = FluxWeights;
   OutputEval.FluxMoments = FluxMoments;
   OutputEval.fidDump     = NULL;
   // DumpFile
   /*if (Output->DumpFile != NULL){
      // check if parallel
      ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
      if (ierr != xf_OK) return ierr;
      
      // use proc number to write to different files in parallel
      strcpy(DumpFile, Output->DumpFile);
      if (nProc > 1) sprintf(DumpFile, "%s.%d\0", Output->DumpFile, myRank);
      
      OutputEval.fidDump = fopen(DumpFile, "w");
      if (OutputEval.fidDump == NULL)
         xf_printf("Warning, could not open the file %s for output.\n", DumpFile);
   }*/
   
   // call boundary residual routine to evaluate the output
   ierr = xf_Error(xf_CalculateResidualBFaces(All, U, NULL, NULL, &OutputEval, nBFG, BFGs, NULL));
   if (ierr != xf_OK) return ierr;

   // close DumpFile
   //if (Output->DumpFile != NULL){
   //   fclose(OutputEval.fidDump);
   //}
  

   // determine Time
   ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "Time", &Time));
   if (ierr != xf_OK) return ierr;

   //need to consider contributions from different processors
   ierr = xf_Error(xf_MPI_Allreduce(&OutputEval.Value[0], 1, xfe_SizeReal, xfe_MPI_SUM));
   if (ierr != xf_OK) return ierr;
   

   ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
   if (ierr != xf_OK) return ierr;
    
   if(myRank == 0)
   {
      sprintf(DumpFile, "BoundaryIntg_%s.dat", Output->Name);
      fout = fopen(DumpFile, "a");
      //only one output value for each boundary integral output
      fprintf(fout, "%.10lf %.10lf\n", Time, OutputEval.Value[0]);

      fclose(fout);
   }

   // correct Value back to positive if was using negative weights
   // do not care for now
   //if ((NegativeFlag) && (Value != NULL))
   //   (*Value) *= -1.;
   
   // include weight in output
   // do not care for now
   //if (Value != NULL) (*Value) *= weight;
   
   xf_Release( (void *) FluxWeights);
   xf_Release( (void *) FluxMoments);
   xf_Release( (void *) BFGs);
   
   return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: Yu_StatisticsOutput
//   borrowed from xflow
int
Yu_StatisticsOutput(xf_All *All, const char *Name, const xf_Vector *U,
                   real *Value, xf_Vector *Value_U, enum xfe_AddType AddFlag1)
{
   int ierr, nOutput, iOutput;
   Yu_Model *Model;
   Yu_Output *Output0, *Output;
   enum xfe_AddType AddFlag, AddFlag2;
   enum xfe_Bool found;
   real weight, LocValue=0.;
   
   //find it from Yu_Model
   Model = All->Model;
   //ierr = xf_Error(xf_FindOutput(All->EqnSet, Name, &Output0));
   //if (ierr != xf_OK) return ierr;
   //find Output from Model
   for(iOutput=0, found=xfe_False; iOutput<Model->nOutput; iOutput++)
   {
      Output = Model->Output + iOutput;
      if(strcmp(Name, Output->Name) == 0)
      {
         found = xfe_True;
         break;
      }
   } //iOutput
   if (!found) return xf_NOT_FOUND;

   //nOutput = ((Output0->Type == xfe_SumOutput) ? Output0->nSumOutput : 1);
   nOutput = 1;

   // zero out Value
   //if (Value != NULL) (*Value) = 0.0;
   
   //temporary removed; Yu; June 2015
   if(Value == NULL) return xf_NULL_OR_STALE_POINTER;

   //output has no relation; this "weight-sum" functionality is deactiviated
   // AddFlag2 is used for additional operations: always either add or sub
   //AddFlag2 = xf_GetAddFlag2(AddFlag1);
   
   for (iOutput=0; iOutput<nOutput; iOutput++){
      
      // determine whether we're setting or adding to Value_U
      //AddFlag = ((iOutput == 0) ? AddFlag1 : AddFlag2);
      
      // pull off desired output (only different if summing outputs)
      //if (Output0->Type == xfe_SumOutput){
      //   ierr = xf_Error(xf_FindOutput(All->EqnSet, Output0->SumOutputNames[iOutput],
      //                                 &Output));
      //   if (ierr != xf_OK) return ierr;
      //   weight = Output0->SumOutputWeights[iOutput];
      //}
      //else{
      //   Output = Output0;
      //   weight = 1.0;
      //}
      weight = 1.0;
      
      switch (Output->Type){
         case xfe_DomainIntegral:
            //return xf_Error(xf_NOT_SUPPORTED);
            ierr = xf_Error(xf_OutputDomainIntegral(All, Output, U, weight, Value,
                                                    Value_U, AddFlag));
            if (ierr != xf_OK) return ierr;
            break;
            
         case xfe_BoundaryIntegral:
            ierr = xf_Error(xf_OutputBoundaryIntegral(All, Output, U, weight, &LocValue,
                                                      Value_U, AddFlag));
            if (ierr != xf_OK) return ierr;
            
            //ierr = xf_Error(xf_MPI_Allreduce(&LocValue, 1, xfe_SizeReal, xfe_MPI_SUM));
            //if (ierr != xf_OK) return ierr;
            break;
            
         case xfe_CutPlaneIntegral:
         case xfe_LineIntegral:
            //ierr = xf_Error(xf_OutputLineCutPlaneIntegral(All, Output, U, weight, &LocValue,
            //                                              Value_U, AddFlag));
            //if (ierr != xf_OK) return ierr;
            
            //ierr = xf_Error(xf_MPI_Allreduce(&LocValue, 1, xfe_SizeReal, xfe_MPI_SUM));
            //if (ierr != xf_OK) return ierr;
            return xf_Error(xf_NOT_SUPPORTED);
            break;
            
         case xfe_PointValue:
            //currently not supported, since not useful
            //return xf_Error(xf_NOT_SUPPORTED);
            ierr = xf_Error(xf_OutputPointValue(All, Output, U, weight, &LocValue,
                                                Value_U, AddFlag));
            if (ierr != xf_OK) return ierr;
            
            //ierr = xf_Error(xf_MPI_Allreduce(&LocValue, 1, xfe_SizeReal, xfe_MPI_SUM));
            //if (ierr != xf_OK) return ierr;
            break;
            
         default:
            return xf_Error(xf_NOT_SUPPORTED);
            break;
      }
      
      
      // add to requested value (weight was included in individual functions)
      //if (Value != NULL) (*Value) += LocValue;
      
   } // iOutput
   
   return xf_OK;
}

/******************************************************************/
//compute pointwise derived quantities for visualizaiton purpose
//
int
Yu_VisualizationOutput(const int sr, const int dim, const int nq,
                       const real *qU, const real *qgU, const real *qout)
{
   int ierr, i, j, k;
   real tmp, vort, fac, gamma, p;
   real *U, gU[50*3], *out;
   real pu_px, pv_py, pw_pz;

   //default value; can vary for different cases
   gamma = 1.40;

   for(i=0; i<nq; i++)
   {
      out = qout + i*Yu_ScalarLast;
      U   = qU + i*sr;
      for(j=0; j<dim; j++)
         for(k=0; k<sr; k++)
            gU[j*sr+k] = qgU[j*nq*sr + i*sr + k];
     
      pu_px = (gU[0*sr+1]-U[1]/U[0]*gU[0*sr+0])/U[0];
      pv_py = (gU[1*sr+2]-U[2]/U[0]*gU[1*sr+0])/U[0];
      if(dim==3)
      pw_pz = (gU[2*sr+3]-U[3]/U[0]*gU[2*sr+0])/U[0];
      else
         pw_pz=0.;
      for(j=0, tmp=0; j<dim; j++)
         tmp += pow(U[j+1]/U[0], 2.0);
      
      p = (gamma-1.0)*(U[dim+1] - 0.5*U[0]*tmp);
      
      ierr = xf_Error(Compute_Vorticity_Mag(sr, dim, U, gU, &vort));
      if(ierr != xf_OK) return ierr;
      
      //pressure
      out[0] = p;
      //vorticity magnitude
      out[1] = vort;
      //out[1] = U[1]/U[0];
      //kinetic energy
      out[2] = 0.5 * U[0] * tmp;
      //enstropy
      out[3] = 0.5 * U[0] * pow(vort, 2.0);
      //dilation; not supported right now
      out[4] = pu_px + pv_py + pw_pz;
      //Mach number
      out[5] = sqrt(tmp)/sqrt(gamma * p/ U[0]);
      //Entropy
      out[6] = p / pow(U[0], gamma);
   }
   


   return xf_OK;
}

//*****************************************************************//
//dump restart datafile and tecplot file in a overwritting manner
int
Yu_OutputPointValueDump(xf_All *All, Yu_Output *Output)
{
   int ierr, i, j, k, l, dim;
   int myRank, nProc, Dim[3];
   real *min, *max, *data, *mature, *xholder;
   real x, y, z;
   float tmp;
   int *Noverlap;
   FILE *ftec, *fout;
   char buf[xf_MAXSTRLEN];
   struct stat st = {0};

   //logistic check
   if(Output->nPoints==0 || Output->PointData==NULL)
   {
      xf_printf("Data requesting dump has not be initiated!\n");
      xf_printf("Dump interval is smaller than sample interval\n");
      return xf_CODE_LOGIC_ERROR;
   }

   //MPI reduce for data gathering 
   ierr = xf_Error(xf_MPI_GetRank(&myRank, &nProc));
   if (ierr != xf_OK) return ierr;
   
   //default setting 
   dim = 3;
   SpecifyPointLocation(Output, xfe_True, Dim);

   //remember the original point coordinates
   ierr = xf_Error(xf_Alloc((void **) &xholder, Output->nPoints*dim, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc((void **) &Noverlap, Output->nPoints, sizeof(int)));
   if (ierr != xf_OK) return ierr;

   for(i=0; i<Output->nPoints; i++)
   {
      xholder[i*dim+0] = Output->PointData[i*(Output->nVars+dim)+0];
      xholder[i*dim+1] = Output->PointData[i*(Output->nVars+dim)+1];
      xholder[i*dim+2] = Output->PointData[i*(Output->nVars+dim)+2];

      if(Output->egrpLocal[i]<0 && Output->elemLocal[i]<0)
         Noverlap[i] = 0;
      else
         Noverlap[i] = 1;
   }

   //reduce sum for Noverlap and data
   ierr = xf_Error(xf_MPI_Allreduce(Noverlap, Output->nPoints, xfe_SizeInt, xfe_MPI_SUM));
   if (ierr != xf_OK) return ierr;
   
  //this one is problemic since the the chunck of data is very LARGE!!
   for(i=0; i<Output->nVars+dim; i++){
      ierr = xf_Error(xf_MPI_Allreduce((Output->PointData + i*(Output->nPoints)), Output->nPoints,
                                        xfe_SizeReal, xfe_MPI_SUM));
      if (ierr != xf_OK) return ierr;
   }

   for(i=0; i<Output->nPoints; i++)
   {
      if(Noverlap[i]<=0) return xf_NOT_SUPPORTED;
      for(j=dim; j<Output->nVars+dim; j++)
         Output->PointData[i*(Output->nVars+dim) + j] /= (real) Noverlap[i];
   
      //put the coordinates back in place
      Output->PointData[i*(Output->nVars+dim)+0] = xholder[i*dim+0]; 
      Output->PointData[i*(Output->nVars+dim)+1] = xholder[i*dim+1]; 
      Output->PointData[i*(Output->nVars+dim)+2] = xholder[i*dim+2]; 
   }

   xf_Release((void *) xholder);  
   xf_Release((void *) Noverlap);
  
   //let the first thread do it; others just return
   if(myRank != 0)
      return xf_OK;

   //check if the target folder exist
   sprintf(buf, "./%s_Statistic", Output->Name);
   if (stat(buf, &st) == -1) {
      //if not exsit, creat one
      mkdir(buf, 0700);
   }

   k = Output->nVars + dim;
   ierr = xf_Error(xf_Alloc( (void **) &min, k, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc( (void **) &max, k, sizeof(real)));
   if (ierr != xf_OK) return ierr;
   ierr = xf_Error(xf_Alloc( (void **) &mature, k, sizeof(real)));
   if (ierr != xf_OK) return ierr;

   //dump binary tec360 file
   if(!Output->SequanceDump)
      sprintf(buf, "./%s_Statistic/statistics_tec_%s.dat", Output->Name, Output->Name);
   else
      sprintf(buf, "./%s_Statistic/statistics_tec_%s_U%d.dat", Output->Name, Output->Name, Output->File_Write_Offset);
   ftec = fopen(buf,"w");
   //write header 
   sprintf(buf, "#!TDV111"); fwrite(buf, sizeof(char), 8, ftec);
   k=1; fwrite(&k, sizeof(int), 1, ftec);
   k=0; fwrite(&k, sizeof(int), 1, ftec);
   sprintf(buf, "TITLE");
   write_string_to_tecplot(ftec, buf, strlen(buf));
   //note the dimension is alway 3D
   k=Output->nVars+3; fwrite(&k, sizeof(int), 1, ftec);

   //variable names
     sprintf(buf,"X"); write_string_to_tecplot(ftec, buf, strlen(buf));
     sprintf(buf,"Y"); write_string_to_tecplot(ftec, buf, strlen(buf));
     sprintf(buf,"Z"); write_string_to_tecplot(ftec, buf, strlen(buf));
     for(i=0; i<Output->nVars; i++){
     k=Output->IVars[i];
     sprintf(buf,OutputQuanityName[k]);write_string_to_tecplot(ftec, buf, strlen(buf));
     }
     tmp = 299.0;fwrite(&tmp, sizeof(float), 1, ftec);
     sprintf(buf, "ZONE_%d", 1); write_string_to_tecplot(ftec, buf, strlen(buf));
     k=-1; fwrite(&k, sizeof(int), 1, ftec);
     k=-2; fwrite(&k, sizeof(int), 1, ftec);
     fwrite(&Output->pretime, sizeof(double), 1, ftec);
     k=-1; fwrite(&k, sizeof(int), 1, ftec);
     k=0; fwrite(&k, sizeof(int), 1, ftec);
     k=0; fwrite(&k, sizeof(int), 1, ftec);
     k=0; fwrite(&k, sizeof(int), 1, ftec);
     k=0; fwrite(&k, sizeof(int), 1, ftec);
     //misc user-defined face; be careful here
     k=0; fwrite(&k, sizeof(int), 1, ftec);
     k=Dim[0]; fwrite(&k, sizeof(int), 1, ftec);
     k=Dim[1]; fwrite(&k, sizeof(int), 1, ftec);
     k=Dim[2]; fwrite(&k, sizeof(int), 1, ftec);
     k=0; fwrite(&k, sizeof(int), 1, ftec);

     //header data seperation
     tmp = 357.0;fwrite(&tmp, sizeof(float), 1, ftec);
     tmp = 299.0;fwrite(&tmp, sizeof(float), 1, ftec);
     //variable precision
     for(i=0; i<Output->nVars+3; i++)
     {k=2; fwrite(&k, sizeof(int), 1, ftec);}
     k=0; fwrite(&k, sizeof(int), 1, ftec);
     k=0; fwrite(&k, sizeof(int), 1, ftec);
     k=-1; fwrite(&k, sizeof(int), 1, ftec);

     //find min and max
     for(l=0; l<Output->nVars+dim; l++)
     {  min[l] = 1.0e+16; max[l] = -1.0e+16;}
     
     //note that the coordinate ordering has been done
     for(i=0; i<Output->nPoints; i++)
     {
        data = Output->PointData + i * (Output->nVars+dim);
        mature[0] = data[0]; mature[1] = data[1]; mature[2] = data[2]; 

        //otherwise, for sequance dump;no need further treatment
        if(!Output->SequanceDump)  
        Yu_OutputScalarMature(Output->Name, Output->nVars, Output->IVars, data+dim, Output->accumtime, mature+dim);
        else
           for(l=dim; l<Output->nVars+dim; l++)
              mature[l] = data[l];

        for(l=0; l<Output->nVars+dim; l++)
        {
            if(mature[l]<min[l]) min[l] = mature[l];
            if(mature[l]>max[l]) max[l] = mature[l];
       
        }
     }
     
     for(l=0; l<Output->nVars+dim; l++)
     {  fwrite(&min[l], sizeof(double), 1, ftec);
        fwrite(&max[l], sizeof(double), 1, ftec);
     }

     for(l=0; l<Output->nVars+dim; l++)
     {
        for(k=0; k<Dim[2]; k++)
           for(j=0; j<Dim[1]; j++)
              for(i=0; i<Dim[0]; i++)
              { 
                 data = Output->PointData + (Dim[1]*Dim[2]*i + Dim[2]*j + k)*(Output->nVars + dim);

                 if(l < dim)  //coordinates
                 {
                    fwrite(&data[l], sizeof(double), 1, ftec);
                    continue;
                 }
                 else
                 {
                    if(!Output->SequanceDump)
                    {
                       Yu_OutputScalarMature(Output->Name, Output->nVars, Output->IVars, data+dim, Output->accumtime, mature+dim);
                       fwrite(&mature[l], sizeof(double), 1, ftec);
                    }
                    else
                       fwrite(&data[l], sizeof(double), 1, ftec);
                 }
              }
     }

     fclose(ftec);

     //save data for restart; save to the corresponding folder
     sprintf(buf, "./%s_Statistic/statistics_%s.bac", Output->Name, Output->Name);
     fout = fopen(buf,"w");

     fwrite(&Output->nPoints, sizeof(int), 1, fout);
     fwrite(&Output->nVars, sizeof(int), 1, fout);
      
     for(i=0; i<Output->nVars; i++)
        fwrite(&Output->IVars[i], sizeof(int), 1, fout);
     
     fwrite(&Output->accumtime, sizeof(double), 1, fout);
     
     //point coordinates
     for(i=0; i<Output->nPoints; i++)
        fwrite(Output->PointData+i*(Output->nVars+dim), sizeof(double), 3, fout);
     
     //point wise data
     for(i=0; i<Output->nPoints; i++)
        fwrite(Output->PointData+i*(Output->nVars+dim)+dim, sizeof(double), Output->nVars, fout);

     fclose(fout);

     xf_Release((void *) min);
     xf_Release((void *) max);
     xf_Release((void *) mature);

   return xf_OK;
}
