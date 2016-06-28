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
  FILE:  xf_SolverIRK.c

  This file contains functions for implicit Runge-Kutta solvers.

*/


#include "xf_MeshMotionGCL.h"

/******************************************************************/
//   FUNCTION Definition: xf_TakeImplicitIRKStep
static int
xf_TakeImplicitIRKStep(xf_All *All, const char *SavePrefix,
		       enum xfe_Bool RestartFlag, real Time, real TimeStep,
		       xf_SolverData *SolverData, xf_Vector **Ui, 
		       enum xfe_TimeSchemeType TimeScheme)
{
  /*
   PURPOSE:
   
     Takes step of an implicit Runge-Kutta method.
   
   INPUTS:
   
     All          : All data structure
     SavePrefix   : prefix for saving files
     RestartFlag  : (not used)
     Time         : physical time at the beginning of the time step
     SolverData   : structure for calling solver routines
     Ui           : Ui[0] = state at beginning of time step
                    Ui[1] = temporary storage state (must exist)
     TimeScheme   : enumerated type indicating which scheme is used
   
   OUTPUTS:

     Ui[0]        : state at end of time step
   
   RETURN: Error code
   
   */


  int ierr,i,j,ncol, iStart = 1;
  int nStage; // Number of stages
  enum xfe_Bool ReuseJacobian;
  enum xfe_Bool UseGCL;
  real CFLStart, c;
  xf_Vector *S;  
  real *a, *dt;
  char title[xf_MAXSTRLEN];
  char PreHeader[xf_MAXSTRLEN];
  xf_Vector **F;  // temporary vectors
  xf_Vector *GCL = NULL, *GCL0 = NULL, **FGCL = NULL, *SGCL = NULL; // GCL vectors
  

  // Declare appropriate coefficients for DIRK3
  real aDIRK3[9] = {0.435866521508459, 0.0, 0.0,
		0.2820667393, 0.435866521508459, 0.0,
		1.2084966492, -0.6443631707, 0.435866521508549};
  
  real dtDIRK3[3]   = {0.435866521508459, 0.7179332608, 1.0};
  
  
  // Declare appropriate coefficients for DIRK4
  real aDIRK4[25] = {1.0/4.0, 0.0, 0.0, 0.0, 0.0,
		 1.0/2.0, 1.0/4.0, 0.0, 0.0, 0.0,
		 17.0/50.0, -1.0/25.0, 1.0/4.0, 0.0, 0.0,
		 371.0/1360.0, -137.0/2720.0, 15.0/544.0, 1.0/4.0, 0.0,
		 25.0/24.0, -49.0/48.0, 125.0/16.0, -85.0/12.0, 1.0/4.0};
  
  real dtDIRK4[5]   = {1.0/4.0, 3.0/4.0, 11.0/20.0, 1.0/2.0, 1.0};
  
  
  // Declare appropriate coefficients for ESDIRK3
  real aESDIRK3[16] = {0.0, 0.0, 0.0, 0.0,
		 1767732205903.0/4055673282236.0, 1767732205903.0/4055673282236.0, 0.0, 0.0,
		 2746238789719.0/10658868560708.0, -640167445237.0/6845629431997.0, 1767732205903.0/4055673282236.0, 0.0,
		 1471266399579.0/7840856788654.0, -4482444167858.0/7529755066697.0, 11266239266428.0/11593286722821.0, 1767732205903.0/4055673282236.0};
  
  real dtESDIRK3[4]    = {0.0, 1767732205903.0/2027836641118.0, 3.0/5.0, 1.0};
  
  
  //Declare appropriate coefficients for ESDIRK4
  real aESDIRK4[36] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.25, 0.25, 0.0, 0.0, 0.0, 0.0,
		 8611.0/62500.0, -1743.0/31250.0, 0.25, 0.0, 0.0, 0.0,
		 5012029.0/34652500.0, -654441.0/2922500.0, 174375.0/388108.0, 0.25, 0.0, 0.0,
		 15267082809.0/155376265600.0, -71443401.0/120774400.0, 730878875.0/902184768.0,2285395.0/8070912.0, 0.25, 0,
		 82889.0/524892.0, 0.0, 15625.0/83664.0, 69875.0/102672.0, -2260.0/8211.0, 0.25};

  real dtESDIRK4[6] = {0.0, 0.5, 83.0/250.0, 31.0/50.0, 17.0/20.0, 1.0};

  
  // Declare appropriate coefficients for ESDIRK5
  real aESDIRK5[64] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 41.0/200.0, 41.0/200.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 41.0/400.0, -567603406766.0/11931857280679.0, 41.0/200.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 683785636431.0/9252920307686.0, 0.0, -110385047103.0/1367015193373.0, 41.0/200.0, 0.0, 0.0, 0.0, 0.0,
		 3016520224154.0/10081342136671.0, 0.0, 30586259806659.0/12414158314087.0, -22760509404356.0/11113319521817.0, 41.0/200.0, 0.0, 0.0, 0.0,
		 218866479029.0/1489978393911.0, 0.0, 638256894668.0/5436446318841.0, -1179710474555.0/5321154724896.0, -60928119172.0/8023461067671.0, 41.0/200.0, 0.0, 0.0,
		 1020004230633.0/5715676835656.0, 0.0, 25762820946817.0/25263940353407.0, -2161375909145.0/9755907335909.0, -211217309593.0/5846859502534.0, -4269925059573.0/7827059040749.0, 41.0/200.0, 0.0,
		 -872700587467.0/9133579230613.0, 0.0, 0.0, 22348218063261.0/9555858737531.0, -1143369518992.0/8141816002931.0, -39379526789629.0/19018526304540.0, 32727382324388.0/42900044865799.0, 41.0/200.0};
		   
  real dtESDIRK5[8] = {0.0, 41.0/100.0, 2935347310677.0/11292855782101.0, 1426016391358.0/7196633302097.0, 92.0/100.0, 24.0/100.0, 3.0/5.0, 1.0};
  
  // Determine the number of stages and set pointers
  switch (TimeScheme){
    case xfe_TimeSchemeDIRK3:
      a = aDIRK3;   dt = dtDIRK3;   ncol = 3; iStart = 0; nStage = 2; break;
    case xfe_TimeSchemeDIRK4:
      a = aDIRK4;   dt = dtDIRK4;   ncol = 5; iStart = 0; nStage = 4; break;
    case xfe_TimeSchemeESDIRK3:
      a = aESDIRK3; dt = dtESDIRK3; ncol = 4; iStart = 1; nStage = 3; break;
    case xfe_TimeSchemeESDIRK4:
      a = aESDIRK4; dt = dtESDIRK4; ncol = 6; iStart = 1; nStage = 5; break;
    case xfe_TimeSchemeESDIRK5:
      a = aESDIRK5; dt = dtESDIRK5; ncol = 8; iStart = 1; nStage = 7; break;
    default:
      return xf_Error(xf_NOT_SUPPORTED); 
      break;
  }
  
  // Initial CFL for artificial time stepping
  ierr = xf_Error(xf_GetKeyValueReal(All->Param->KeyValue, "CFL", &CFLStart));
  if (ierr != xf_OK) return ierr;
  
  // Do NOT reuse R_U for now
  ReuseJacobian = xfe_False;

  // determine if using a Geometric Conservation Law
  ierr = xf_Error(xf_GetKeyValueBool(All->Param->KeyValue, "UseGCL", &UseGCL));
  if (ierr != xf_OK) return ierr;
  if (UseGCL){
    ierr = xf_Error(xf_InitMeshMotionGCLVector(All, Ui[0], -1, xfe_False, &GCL));
    if (ierr != xf_OK) return ierr;
  }
  
  // locate array of required update vectors, F
  ierr = xf_Error(xf_Alloc( (void **) &F, nStage+1, sizeof(xf_Vector *)));
  if (ierr != xf_OK) return ierr;
  for (i = 0; i <= nStage; i++){
    sprintf(title, "F%d_IRK", i);
    ierr = xf_Error(xf_FindSimilarVector(All, Ui[0], title, xfe_True, xfe_True,
					 NULL, F + i, NULL));
    if (ierr != xf_OK) return ierr;
  }
  S = F[nStage]; // will serve as source vector
  
  // locate temporary vectors for GCL
  if (UseGCL){
    sprintf(title, "ImplicitTimeStepGCL0_%d", i);
    ierr = xf_Error(xf_FindSimilarVector(All, GCL, title, xfe_True, xfe_True, NULL, &GCL0, NULL));
    if (ierr != xf_OK) return ierr;
    ierr = xf_Error(xf_Alloc( (void **) &FGCL, nStage+1, sizeof(xf_Vector *)));
    if (ierr != xf_OK) return ierr;
    for (i=0; i<=nStage; i++){
      sprintf(title, "ImplicitTimeStepFGCL_%d", i);
      ierr = xf_Error(xf_FindSimilarVector(All, GCL, title, xfe_True, xfe_True, NULL, FGCL+i, NULL));
      if (ierr != xf_OK) return ierr;
    }
    SGCL = FGCL[nStage];
  }


  // Set Ui[1] = Ui[0] (the current state)    
  ierr = xf_Error(xf_SetVector(Ui[0], xfe_Set, Ui[1]));
  if (ierr != xf_OK) return ierr;

  // same for GCL
  if (UseGCL){
    // set GCL0 = GCL
    ierr = xf_Error(xf_SetVector(GCL, xfe_Set, GCL0));
    if (ierr != xf_OK) return ierr;
  }

  // begin loop over stages
  for (i=iStart; i<=nStage; i++){
    
    // calculate residual at previous state
    if (i != 0){
      // Set time to Time+dt[i-1]*TimeStep
      ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time",Time + dt[i-1]*TimeStep));
      if (ierr != xf_OK) return ierr;
      
      if (UseGCL){ // GCL comes first
	ierr = xf_MeshMotionGCLResidual(All, GCL, xfe_True, FGCL[i-1]);
	if (ierr != xf_OK) return ierr;
      }
      ierr = xf_CalculateResidual(All, Ui[0], F[i-1], NULL, NULL);
      if (ierr != xf_OK) return ierr;
    } 

    // Construct S = -1*u^{n}
    ierr = xf_Error(xf_VectorMultSet(Ui[1], -1.0, xfe_Set, S));
    if (ierr != xf_OK) return ierr;
    // Set S = M/dt * S
    ierr = xf_Error(xf_MultMassMatrix(All, 1.0/TimeStep, S));
    if (ierr != xf_OK) return ierr;

    // same with GCL
    if (UseGCL){
      ierr = xf_Error(xf_VectorMultSet(GCL0, -1.0, xfe_Set, SGCL));
      if (ierr != xf_OK) return ierr;
      ierr = xf_Error(xf_MultMassMatrix(All, 1.0/TimeStep, SGCL));
      if (ierr != xf_OK) return ierr;
    }
    
    for (j=0; j<i; j++){
      // Add previously-computed residuals to S: S = -u^{n}/dt + a[i][j]*F[j]
      ierr = xf_Error(xf_VectorMultSet(F[j],  a[i*ncol+j], xfe_Add, S));
      if (ierr != xf_OK) return ierr;
      if (UseGCL){
	ierr = xf_Error(xf_VectorMultSet(FGCL[j],  a[i*ncol+j], xfe_Add, SGCL));
	if (ierr != xf_OK) return ierr;
      }
    }

    // Set S = 1/a[i][i] * S
    ierr = xf_Error(xf_VectorMult(S, 1.0/a[i*ncol+i]));
    if (ierr != xf_OK) return ierr;
    if (UseGCL){
      ierr = xf_Error(xf_VectorMult(SGCL, 1.0/a[i*ncol+i]));
      if (ierr != xf_OK) return ierr;
    }    

    c = 1.0/(a[i*ncol+i]*TimeStep); // coefficient on M*u^{i}

    // Set initial CFL for pseudo-time stepping
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart));
    if (ierr != xf_OK) return ierr;

    // Write out header for this time step
    sprintf(PreHeader, "%s %.15E  %s stage = %d", "%% Time = ", Time+dt[i]*TimeStep, 
	    xfe_TimeSchemeName[TimeScheme], i);
    ierr = xf_Error(xf_WriteLogHeader(All, PreHeader));
    if (ierr != xf_OK) return ierr;


    // Set time to Time+dt[i]*TimeStep
    ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "Time",Time+dt[i]*TimeStep));
    if (ierr != xf_OK) return ierr;

    // Call non-linear solver
    if (UseGCL){
      // GCL first; this is a trivial linear "solve" with the mass matrix
      ierr = xf_Error(xf_MeshMotionGCLSolveSystem(All, c, SGCL, GCL));
      if (ierr != xf_OK) return ierr;
    }
    ierr = xf_Error(xf_SolveNonlinearSystem(All, c, ReuseJacobian, S, Ui+0));
    if (ierr != xf_OK){
      xf_printf("Error: nonlinear solve failed at Time = %.6E in IRK step\n", Time);
      return ierr;
    }

  } // i over stages

  // Set initial CFL or pseudo-time stepping
  ierr = xf_Error(xf_SetKeyValueReal(All->Param->KeyValue, "CFL", CFLStart));
  if (ierr != xf_OK) return ierr;

  xf_Release((void *) F);
  if (UseGCL) xf_Release( (void *) FGCL);
  
  return xf_OK;

}

