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
  FILE:  xf_MeshMotionAnalytical.c

  This file contains functions for analytical mesh motion.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_String.h"
#include "xf_Param.h"
#include "xf_Math.h"
#include "xf_MeshMotion.h"

// input/output structures
typedef struct
{
  int dim;
  real Time;
  real *X;
}
xf_AnaInData;

typedef struct
{
  real *x;
  real *H;
  real *H_X;
  real *x_t;
}
xf_AnaOutData;


/******************************************************************/
//   FUNCTION Definition: xf_InitAnaMotions
static int 
xf_InitAnaMotions( xf_AnaMotions *AnaMotions){

  int ierr;

  AnaMotions->nTerm        = 0;
  AnaMotions->AnaMotion    = NULL;
  AnaMotions->BlendType    = xfe_AnaBlend_None;
  AnaMotions->BlendRParam  = NULL;
  AnaMotions->BlendIParam  = NULL;

  /* Initialize key-value structure */
  ierr = xf_Error(xf_InitKeyValue(&AnaMotions->BlendKeyValue));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_CreateAnaMotions
int 
xf_CreateAnaMotions( xf_AnaMotions **pAnaMotions){

  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pAnaMotions, 1, sizeof(xf_AnaMotions)));
  if (ierr != xf_OK) return ierr;

  ierr = xf_Error(xf_InitAnaMotions((*pAnaMotions)));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_CreateAnaMotionsSet
int 
xf_CreateAnaMotionsSet( xf_AnaMotionsSet **pAnaMotionsSet){

  int ierr;

  ierr = xf_Error(xf_Alloc((void **) pAnaMotionsSet, 1, sizeof(xf_AnaMotionsSet)));
  if (ierr != xf_OK) return ierr;
  
  (*pAnaMotionsSet)->nMotion       = 0;
  (*pAnaMotionsSet)->AnaMotions    = NULL;
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AllocAnaMotions
int 
xf_AllocAnaMotions( xf_AnaMotions *AnaMotions, int nTerm)
{
  int ierr, i;
  xf_AnaMotion *pAnaMotion;

  ierr = xf_Error(xf_Alloc((void **) &AnaMotions->AnaMotion, nTerm, sizeof(xf_AnaMotion)));
  if (ierr != xf_OK) return ierr;

  AnaMotions->nTerm = nTerm;
  
  for (i=0; i<nTerm; i++){
    pAnaMotion               = AnaMotions->AnaMotion+i;
    pAnaMotion->Name         = NULL;
    pAnaMotion->MotionType   = xfe_AnaMotionLast;
    pAnaMotion->MotionRParam = NULL;
    pAnaMotion->MotionIParam = NULL;

    /* Initialize key-value structure */
    ierr = xf_Error(xf_InitKeyValue(&pAnaMotion->MotionKeyValue));
    if (ierr != xf_OK) return ierr;
  } // i
  
  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_AllocAnaMotionsSet
int 
xf_AllocAnaMotionsSet( xf_AnaMotionsSet *AnaMotionsSet, int nMotion)
{
  int ierr, i;
  xf_AnaMotion *pAnaMotionSet;

  ierr = xf_Error(xf_Alloc((void **) &AnaMotionsSet->AnaMotions, nMotion, sizeof(xf_AnaMotions)));
  if (ierr != xf_OK) return ierr;

  AnaMotionsSet->nMotion = nMotion;
  
  for (i=0; i<nMotion; i++){
    ierr = xf_Error(xf_InitAnaMotions(AnaMotionsSet->AnaMotions+i));
    if (ierr != xf_OK) return ierr;
  } // i
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyAnaMotion
static int 
xf_DestroyAnaMotion( xf_AnaMotion *AnaMotion){

  int ierr;

  xf_Release( (void *) AnaMotion->Name);
  xf_Release( (void *) AnaMotion->MotionRParam);
  xf_Release( (void *) AnaMotion->MotionIParam);

  /* Destroy key-value structure */
  ierr = xf_Error(xf_DestroyKeyValue(&AnaMotion->MotionKeyValue));
  if (ierr != xf_OK) return ierr;
  
  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_DestroyAnaMotions
static int 
xf_DestroyAnaMotions( xf_AnaMotions *AnaMotions, enum xfe_Bool DestroySelf){

  int ierr, i;

  if (AnaMotions == NULL) return xf_OK;
  
  for (i=0; i<AnaMotions->nTerm; i++){
    ierr = xf_Error(xf_DestroyAnaMotion(AnaMotions->AnaMotion+i));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *) AnaMotions->AnaMotion);
  
  xf_Release( (void *) AnaMotions->BlendRParam);
  xf_Release( (void *) AnaMotions->BlendIParam);

  /* Destroy key-value structure */
  ierr = xf_Error(xf_DestroyKeyValue(&AnaMotions->BlendKeyValue));
  if (ierr != xf_OK) return ierr;  

  // destroy self only if requested
  if (DestroySelf) xf_Release((void *) AnaMotions);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_DestroyAnaMotionsSet
int 
xf_DestroyAnaMotionsSet( xf_AnaMotionsSet *AnaMotionsSet){

  int ierr, i;

  if (AnaMotionsSet == NULL) return xf_OK;
  
  for (i=0; i<AnaMotionsSet->nMotion; i++){
    ierr = xf_Error(xf_DestroyAnaMotions(AnaMotionsSet->AnaMotions+i, xfe_False));
    if (ierr != xf_OK) return ierr;
  }
  xf_Release((void *) AnaMotionsSet->AnaMotions);
  
  xf_Release((void *) AnaMotionsSet);

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ApplyAnaMotion_Plunge
static int
xf_ApplyAnaMotion_Plunge(xf_AnaMotion *AnaMotion, xf_AnaInData *In, 
			 xf_AnaOutData *Out)
{
  int ierr, d, dim, nset;
  real Time, w, phi, *A, *offset;
  enum {XAmplitude, YAmplitude, ZAmplitude,  
	XOffset, YOffset, ZOffset,
	Frequency, Phase, RPLast};
  char *RPName[RPLast] = {"XAmplitude", "YAmplitude", "ZAmplitude", 
			  "XOffset", "YOffset", "ZOffset", 
			  "Frequency", "Phase"};
  real RPDef[RPLast] = {0., 0., 0., 0., 0., 0., 0., 0.};
  real *RP;

  // fill in if do not have parameters in shortcut storage
  if (AnaMotion->MotionRParam == NULL){
    xf_AllocFillRParam(AnaMotion->MotionKeyValue, RPLast, RPName, RPDef, 
		       &AnaMotion->MotionRParam, &nset);
    if (nset != (d=AnaMotion->MotionKeyValue.nKey)){
      xf_printf("Error: %d key-value(s) not understood in Plunge motion.\n", d-nset);
      return xf_Error(xf_INPUT_ERROR);
    }
  }
  RP = AnaMotion->MotionRParam;

  dim  = In->dim;
  Time = In->Time;
  w      = RP[Frequency];
  phi    = RP[Phase];
  A      = RP+XAmplitude;
  offset = RP+XOffset;

  // new position, x
  if (Out->x != NULL)
    for (d=0; d<dim; d++)
      Out->x[d] = In->X[d] + offset[d] + A[d]*sin(w*Time + phi);

  // H = x_X = Identity
  if (Out->H != NULL){
    for (d=0; d<dim*dim; d++) Out->H[d] = 0.;
    for (d=0; d<dim*dim; d+=(dim+1)) Out->H[d] = 1.;
   }
  
  // H_X = x_XX = 0
  if (Out->H_X != NULL)
    for (d=0; d<dim*dim*dim; d++) Out->H_X[d] = 0.;

  // x_t
  if (Out->x_t != NULL)
    for (d=0; d<dim; d++)
      Out->x_t[d] = A[d]*w*cos(w*Time + phi);

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyAnaMotion_Pitch
static int
xf_ApplyAnaMotion_Pitch(xf_AnaMotion *AnaMotion, xf_AnaInData *In, xf_AnaOutData *Out)
{
  int ierr, d, d0, d1, d2, dim, nset, divcount;
  real Time, w, phi, *A, *Xc, Tt;
  real theta, ct, st, dp0, dp1, wo;
  real theta_t, wn, a, b, f;
  real f_wn, a_wn, b_wn, transtol;
  real Ar, Av[3];
  real X[3], x[3];
  static real w1 =0.0, A1=0.0;

  enum {XOrigin, YOrigin, ZOrigin, 
	XAmplitude, YAmplitude, ZAmplitude,
	Frequency, Phase, TransientTime, RPLast};
  char *RPName[RPLast] = {"XOrigin", "YOrigin", "ZOrigin", 
			  "XAmplitude", "YAmplitude", "ZAmplitude", 
			  "Frequency", "Phase", "TransientTime"};
  real RPDef[RPLast] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  real *RP;

  // fill in if do not have parameters in shortcut storage
  if (AnaMotion->MotionRParam == NULL){
    xf_AllocFillRParam(AnaMotion->MotionKeyValue, RPLast, RPName, RPDef, 
		       &AnaMotion->MotionRParam, &nset);
    if (nset != (d=AnaMotion->MotionKeyValue.nKey)){
      xf_printf("Error: %d key-value(s) not understood in Pitch motion.\n", d-nset);
      return xf_Error(xf_INPUT_ERROR);
    }
  }
  RP = AnaMotion->MotionRParam;

  dim    = In->dim;
  Time   = In->Time;
  Xc     = RP + 0;
  A      = RP+XAmplitude;
  w      = RP[Frequency];
  phi    = RP[Phase];
  Tt     = RP[TransientTime];

  for (d=0,Ar=0.; d<3; d++) Ar+=A[d]*A[d]; Ar = sqrt(Ar); // amplitude of rot

  if (Ar != 0.){
    for (d=0; d<3; d++) Av[d] = A[d]/Ar; // direction of rot
  }
  else{
    for (d=0; d<3; d++) Av[d] = 0.0; //to avoid dividing by 0 if Ar=0
  }
 
  // determine blending amplitude and frequency 

  if ((Tt > 0.0) && (w1 == 0.0)){
    transtol = 1e-5; //tolerance for newton solve
    divcount = 0; 
    wo = 2*w; //initial guess for newton solve
    wn = wo; 
    f = 1.0; //initial value to get loop running

    while (fabs(f) > transtol){
      a = Ar*sin(w*Tt + phi)/(sin(wn*Tt-PI/2)+1);
      b = cos(wn*Tt-PI/2);
      a_wn = -Ar*sin(w*Tt + phi)*Tt*cos(wn*Tt-PI/2)/((sin(wn*Tt-PI/2)+1)*(sin(wn*Tt-PI/2)+1));
      b_wn = -Tt*sin(wn*Tt-PI/2);
      f = a*b*wn - Ar*w*cos(w*Tt + phi); //system to find roots of
      f_wn = a_wn*b*wn + a*b + a*b_wn*wn;
      wn = wn - f/f_wn; //update root estimate

      if ((f <= transtol) && (fabs(wn) > 2*PI/Tt)){
        divcount+= 1;
        f = 1.0;
        wn = pow((2.0/3.0),divcount)*wo; //wrong root; change initial guess and try again
      }
    }
    
    w1 = wn; //set frequency of blending sinusoid
    A1 = Ar*sin(w*Tt + phi)/(sin(w1*Tt-PI/2)+1); // set amplitude of blending sinusoid
  }

  if ((Tt > 0.) && (Time <= Tt)){
    theta = A1*(sin(w1*Time-PI/2) + 1);
    theta_t = A1*w1*cos(w1*Time-PI/2);
  }
  else{
    theta = Ar*sin(w*Time + phi);
    theta_t = Ar*w*cos(w*Time + phi); // desired sinusoidal motion
  }
    
  ct = cos(theta);
  st = sin(theta);

  // set 3D X (will use a generic formula)
  for (d=0; d<3; d++) X[d] = ((d<dim) ? In->X[d] : 0.);

  // new position, x
  if (Out->x != NULL){
    dp0 = Av[0]*Xc[0] + Av[1]*Xc[1] + Av[2]*Xc[2];
    dp1 = Av[0]*X[0] + Av[1]*X[1] + Av[2]*X[2];
    for (d=0; d<3; d++){
      d1 = (d+1)%3;
      d2 = (d+2)%3;
      x[d] = ((Xc[d]*(1. - Av[d]*Av[d]) 
	       - Av[d]*(dp0-Xc[d]*Av[d] - dp1))*(1.-ct)
	      +X[d]*ct + (Av[d1]*(X[d2]-Xc[d2])-Av[d2]*(X[d1]-Xc[d1]))*st); 
    } // d
    for (d=0; d<dim; d++)
      Out->x[d] = x[d];
  }

  // H = x_X = first derivative
  if (Out->H != NULL){

    if (dim == 3){
      Out->H[0] = Av[0]*Av[0]*(1.-ct) + ct;
      Out->H[1] = Av[0]*Av[1]*(1.-ct) - Av[2]*st;
      Out->H[2] = Av[0]*Av[2]*(1.-ct) + Av[1]*st;
      Out->H[3] = Av[1]*Av[0]*(1.-ct) + Av[2]*st;
      Out->H[4] = Av[1]*Av[1]*(1.-ct) + ct;
      Out->H[5] = Av[1]*Av[2]*(1.-ct) - Av[0]*st;
      Out->H[6] = Av[2]*Av[0]*(1.-ct) - Av[1]*st;
      Out->H[7] = Av[2]*Av[1]*(1.-ct) + Av[0]*st;
      Out->H[8] = Av[2]*Av[2]*(1.-ct) + ct;
    }
    else if (dim == 2){
      Out->H[0] = Av[0]*Av[0]*(1.-ct) + ct;
      Out->H[1] = Av[0]*Av[1]*(1.-ct) - Av[2]*st;
      Out->H[2] = Av[1]*Av[0]*(1.-ct) + Av[2]*st;
      Out->H[3] = Av[1]*Av[1]*(1.-ct) + ct;
      Out->H[4] = 0.;
      Out->H[5] = 0.;
      Out->H[6] = 0.;
      Out->H[7] = 0.;
      Out->H[8] = 0.;
    }
    
  }
  
  // H_X = x_XX = 0
  if (Out->H_X != NULL)
    for (d=0; d<27; d++) Out->H_X[d] = 0.; //second derivative is zero 

  // x_t
  if (Out->x_t != NULL){
    Out->x_t[0] = ((-st + Av[0]*Av[0]*st)*(X[0]-Xc[0]) + (Av[0]*Av[1]*st - Av[2]*ct)*(X[1]-Xc[1]) + (Av[0]*Av[2]*st + Av[1]*ct)*(X[2]-Xc[2]))*theta_t;
    Out->x_t[1] = ((Av[1]*Av[0]*st + Av[2]*ct)*(X[0]-Xc[0]) + (-st + Av[1]*Av[1]*st)*(X[1]-Xc[1]) + (Av[1]*Av[2]*st - Av[0]*ct)*(X[2]-Xc[2]))*theta_t;
    Out->x_t[2] = ((Av[2]*Av[0]*st - Av[1]*ct)*(X[0]-Xc[0]) + (Av[2]*Av[1]*st + Av[0]*ct)*(X[1]-Xc[1]) + (-st + Av[2]*Av[2]*st)*(X[2]-Xc[2]))*theta_t;
  }

  return xf_OK;

}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyAnaMotion_PerssonFS
static int
xf_ApplyAnaMotion_PerssonFS(xf_AnaMotion *AnaMotion, xf_AnaInData *In, 
			 xf_AnaOutData *Out)
{
  int ierr, d, dim, nset;
  real Time, w, phi, *A, *offset, X[3], to;
  enum {To, RPLast};
  char *RPName[RPLast] = {"To"};
  real RPDef[RPLast] = {1.0};
  real sinx0, sinx1, sint2, sint4;
  real cosx0, cosx1, cost2, cost4;
  real *RP;

  // fill in if do not have parameters in shortcut storage
  if (AnaMotion->MotionRParam == NULL){
    xf_AllocFillRParam(AnaMotion->MotionKeyValue, RPLast, RPName, RPDef, 
		       &AnaMotion->MotionRParam, &nset);
    if (nset != (d=AnaMotion->MotionKeyValue.nKey)){
      xf_printf("Error: %d key-value(s) not understood in PerssonFS motion.\n", d-nset);
      return xf_Error(xf_INPUT_ERROR);
    }
  }
  RP = AnaMotion->MotionRParam;

  dim  = In->dim;
  Time = In->Time;
  to = RP[To];

  // set X 
  for (d=0; d<2; d++) X[d] = ((d<dim) ? In->X[d] : 0.);

  // convenient pre-computations
  sinx0 = sin(2.0*PI*X[0]/20);  cosx0 = cos(2.0*PI*X[0]/20); 
  sinx1 = sin(PI*X[1]/7.5);	cosx1 = cos(PI*X[1]/7.5);    
  sint2 = sin(2*PI*Time/to);	cost2 = cos(2*PI*Time/to);   
  sint4 = sin(4*PI*Time/to);	cost4 = cos(4*PI*Time/to);

  // new position, x
  if (Out->x != NULL){
    Out->x[0] = X[0] + 2.0*sinx0*sinx1*sint2;
    Out->x[1] = X[1] + 1.5*sinx0*sinx1*sint4;
  }

  // H = x_X 
  if (Out->H != NULL){
    Out->H[0] = 1.0 + sinx1*sint2*(4*PI/20)*cosx0;
    Out->H[1] = 2.0*sinx0*sint2*cosx1*(PI/7.5);
    Out->H[2] = sinx1*sint4*(1.5*2*PI/20)*cosx0;
    Out->H[3] = 1.0 + 1.5*sinx0*sint4*(PI/7.5)*cosx1;
   }
  
  // H_X = x_XX
  if (Out->H_X != NULL){

    Out->H_X[0] = -1.0*sinx1*sint2*(4*PI/20)*(2*PI/20)*sinx0;
    Out->H_X[1] = sint2*(4*PI/20)*cosx0*(PI/7.5)*cosx1;
    Out->H_X[2] = sinx1*sint4*(1.5*2*PI/20)*(2*PI/20)*sinx0*(-1.0);
    Out->H_X[3] = (PI/7.5)*cosx1*sint4*(1.5*2*PI/20)*cosx0;
    Out->H_X[4] = sint2*(4*PI/20)*cosx0*(PI/7.5)*cosx1;
    Out->H_X[5] = (2*PI/7.5)*sinx0*sint2*sinx1*(-PI/7.5);
    Out->H_X[6] = (PI/7.5)*cosx1*sint4*(1.5*2*PI/20)*cosx0;
    Out->H_X[7] = 1.5*sinx0*sint4*(PI/7.5)*(PI/7.5)*sinx1*(-1.0);

  }

  // x_t
  if (Out->x_t != NULL){

    Out->x_t[0] = 2.0*sinx0*sinx1*(2*PI/to)*cost2;
    Out->x_t[1] = 1.5*sinx0*sinx1*(4*PI/to)*cost4;
  }

  return xf_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_ApplyAnaMotion
static int
xf_ApplyAnaMotion(xf_AnaMotion *AnaMotion, xf_AnaInData *In, xf_AnaOutData *Out){
  int ierr;

  switch (AnaMotion->MotionType){
  case xfe_AnaMotion_Plunge:
    ierr = xf_Error(xf_ApplyAnaMotion_Plunge(AnaMotion, In, Out));
    break;
  case xfe_AnaMotion_Pitch:
    ierr = xf_Error(xf_ApplyAnaMotion_Pitch(AnaMotion, In, Out));
    break;
  case xfe_AnaMotion_PerssonFS:
    ierr = xf_Error(xf_ApplyAnaMotion_PerssonFS(AnaMotion, In, Out));
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_BlendMotion
static int
xf_BlendMotion(xf_AnaMotions *Motions, int dim, const real *X, real *pb, 
	       real *b_X, real *b_XX, enum xfe_Bool *OutsideBlend)
{
  int ierr, i, j, k, nsetR, nsetI;
  enum xfe_AnaBlendType BlendType;
  xf_KeyValue BlendKeyValue;
  enum {XCenter, YCenter, ZCenter, Radius, Dist, RPLast};
  enum {ExtrudeInX, ExtrudeInY, ExtrudeInZ, IPLast};
  char *RPName[RPLast] = {"XCenter", "YCenter", "ZCenter", "Radius", "Dist"}; //Move ExtrudeMotion to IParam vector eventually, since it's an integer
  real RPDef[RPLast] = {0., 0., 0., 1., 10.0};
  real *RP, **pRP, *Xc;
  char *IPName[IPLast] = {"ExtrudeMotionInX", "ExtrudeMotionInY", "ExtrudeMotionInZ"};
  int IPDef[IPLast] = {0, 0, 0}; //default values; no extrusion
  int *IP, **pIP;
  int ExtrudeFlag[3];
  real Rc, D, d, d2, z, z2, z3, z4, r, tmp;
  real r_z, r_zz;
  real d_X[3], d_XX[9];

  BlendType     = Motions->BlendType;
  BlendKeyValue = Motions->BlendKeyValue;
  pRP           = &Motions->BlendRParam;
  pIP           = &Motions->BlendIParam;

  // initialize to zero
  (*pb) = 0.;
  if (b_X  != NULL) for (i=0; i<dim; i++) b_X[i] = 0.;
  if (b_XX != NULL) for (i=0; i<dim*dim; i++) b_XX[i] = 0.;

  // initialize OutsideBlend: assume affected by blending
  if (OutsideBlend != NULL) (*OutsideBlend) = xfe_False; 


  if (BlendType == xfe_AnaBlend_None) return xf_OK; // nothing else to do

  if (((*pRP) == NULL)||((*pIP) == NULL)){//fill both real and int params; "or" statement since we need both nsetR and nsetI for check
    xf_AllocFillRParam(BlendKeyValue, RPLast, RPName, RPDef, pRP, &nsetR);
    xf_AllocFillIParam(BlendKeyValue, IPLast, IPName, IPDef, pIP, &nsetI);     
     
    if((nsetR + nsetI) != (k=BlendKeyValue.nKey)){ //if total number of keys actually set not equal to number read from file 
      xf_printf("Error: %d key-value(s) not understood in Blend motion. \n", k-(nsetR+nsetI)); 
      return xf_Error(xf_INPUT_ERROR); 
    } 

  }

  RP = (*pRP);
  Xc = RP + XCenter;   // center of blending
  Rc = RP[Radius];     // radius of rigid body motion
  D  = RP[Dist];       // blending takes place over this distance

  IP = (*pIP);
  ExtrudeFlag[0] = IP[ExtrudeInX]; //1 if extruding mesh motion in x direction, 0 otherwise
  ExtrudeFlag[1] = IP[ExtrudeInY];
  ExtrudeFlag[2] = IP[ExtrudeInZ]; 
  
  // compute d(X) and derivatives
  //for (i=0, d2=0.; i<dim-ExtrudeInZFlag; i++) d2 += (X[i]-Xc[i])*(X[i]-Xc[i]);
 
  for(i=0, d2=0.; i<dim; i++)
    d2 += ((ExtrudeFlag[i]) ? 0. : (X[i]-Xc[i])*(X[i]-Xc[i]) );

  d = sqrt(d2);
  if (d >= Rc){
    if (d == 0.) return xf_Error(xf_SINGULAR);
    if (b_X != NULL) for (i=0; i<dim; i++) d_X[i] = (X[i]-Xc[i])/d;
    if (b_XX != NULL) 
      for (i=0; i<dim; i++){
	for (j=0, tmp=-(X[i]-Xc[i])/(d2*d); j<dim; j++) 
	  d_XX[i*dim+j] = tmp*(X[j]-Xc[j]);
	d_XX[i*dim+i] += 1./d;
      }
  }
  d = d - Rc; // actual definition, above was for derivative convenience

  // compute r, r_z, r_zz, z = d/D;
  z = d/D; z2 = z*z; z3 = z2*z; z4=z3*z;
  switch (BlendType){
  case xfe_AnaBlend_Cubic:
    r = 3.*z2 - 2.*z3;
    if (b_X != NULL) r_z = 6.*z - 6.*z2;
    if (b_XX != NULL) r_zz = 6. - 12.*z;
    break;
  case xfe_AnaBlend_Quintic:
    r = 10.*z3 - 15.*z2*z2 + 6.*z3*z2;
    if (b_X != NULL) r_z = 30.*z2 - 60.*z3 + 30.*z2*z2;
    if (b_XX != NULL) r_zz = 60.*z - 180.*z2 + 120.*z3;
    break;
  case xfe_AnaBlend_Septic:
    r = -20.*z4*z3 + 70.*z3*z3 - 84*z3*z2 + 35*z2*z2;
    if (b_X != NULL) r_z = -140.*z3*z3 + 420.*z3*z2 - 420.*z4 + 140.*z3;
    if (b_XX != NULL) r_zz = -840.*z3*z2 + 2100.*z4 - 1680.*z3 + 420.*z2;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }
  
  // set b, b_X, b_XX
  if (d < 0.)
    (*pb) = 0.;
  else if (d>D){
    (*pb) = 1.;
    // outside blending region in this case
    if (OutsideBlend != NULL) (*OutsideBlend) = xfe_True; 
  }
  else{
    (*pb) = r;
    if (b_X  != NULL){
      for (i=0; i<dim; i++){ // skip certain dimensions if extruding
	b_X[i] = ( (ExtrudeFlag[i]) ? 0. : r_z/D*d_X[i]);
      }
    }
    if (b_XX != NULL) 
      for (i=0; i<dim; i++){
        if(ExtrudeFlag[i])
          continue;
        else{
	  for (j=0; j<dim; j++){
	    b_XX[i*dim+j] = ( (ExtrudeFlag[j]) ? 0. : r_zz/(D*D)*d_X[j]*d_X[i] + r_z/D*d_XX[i*dim+j]);
          } // j
	} // else     
      } // i
  }

  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_dimMxM
static void
xf_dimMxM(const real *A, const real *B, int d, real *C)
{
  int i, j, k, id;
  // C = A*B, all matrices dxd squares, stored unrolled in row-major
  for (i=0; i<d; i++)
    for (j=0, id=i*d; j<d; j++)
      for (k=0, C[id+j]=0.; k<d; k++)
	C[id+j] += A[id+k]*B[k*d+j];
}



/******************************************************************/
//   FUNCTION Definition: xf_MapDet
static int
xf_MapDet(const real *G, const real *G_X, int dim, real *pg, real *gig_X)
{
  // computes g = det(G) and gig_X (vector) based on G_X (tensor)
  int ierr;
  int k, kd;
  real g;
  
  if (dim == 1)
    g = G[0];
  else if (dim == 2)
    g = G[0]*G[3]-G[1]*G[2];
  else if (dim == 3)
    g =  G[0]*G[4]*G[8] + G[1]*G[5]*G[6] + G[2]*G[3]*G[7]
      -  G[0]*G[7]*G[5] - G[3]*G[1]*G[8] - G[6]*G[4]*G[2];
  else return xf_Error(xf_NOT_SUPPORTED);

  if (gig_X != NULL){
    if (G_X == NULL) return xf_Error(xf_INPUT_ERROR);
    for (k=0; k<dim; k++) gig_X[k] = 0.;
    if (dim == 1)
      gig_X[0] = G_X[0];
    else if (dim == 2)
      for (k=0; k<2; k++){
	kd = k*4;
	gig_X[k] += G_X[0+kd]*G[3] + G[0]*G_X[3+kd] - G_X[1+kd]*G[2] - G[1]*G_X[2+kd];
      }
    else{
      for (k=0; k<3; k++){
	kd = k*9;	
	gig_X[k] += G_X[0+kd]*G[4]*G[8] + G[0]*G_X[4+kd]*G[8] + G[0]*G[4]*G_X[8+kd]
	  + G_X[1+kd]*G[5]*G[6] + G[1]*G_X[5+kd]*G[6] + G[1]*G[5]*G_X[6+kd]
	  + G_X[2+kd]*G[3]*G[7] + G[2]*G_X[3+kd]*G[7] + G[2]*G[3]*G_X[7+kd]
	  - G_X[0+kd]*G[7]*G[5] - G[0]*G_X[7+kd]*G[5] - G[0]*G[7]*G_X[5+kd]
	  - G_X[3+kd]*G[1]*G[8] - G[3]*G_X[1+kd]*G[8] - G[3]*G[1]*G_X[8+kd]
	  - G_X[6+kd]*G[4]*G[2] - G[6]*G_X[4+kd]*G[2] - G[6]*G[4]*G_X[2+kd];
      }
    }
    for (k=0; k<dim; k++) gig_X[k] /= g;
  }
      
  // general code, not as intuitive
/*   int  I2[2][2] = {{0,3},{1,2}}; */
/*   real s2[2]    = {1., -1.}; */
/*   int  I3[6][3] = {{0,4,8}, {1,5,6}, {2,3,7}, {0,7,5}, {3,1,8}, {6,4,2}}; */
/*   real s3[6]    = {1., 1., 1., -1. -1. -1.}; */
/*   // determinant, g */
/*   for (i=0, g=0.; i<n; i++){ */
/*     for (j=0, t=s[i].; j<dim; j++) */
/*       t *= G[I[i*dim+j]]; */
/*     g += t; */
/*   } */
  

/*   // linearization, g_X */
/*   if (g_X != NULL){ */
/*     if (G_X == NULL) return xf_Error(xf_INPUT_ERROR);; */
/*     for (k=0; k<dim; k++) */
/*       for (i=0, g_X[k]=0.; i<n; i++) */
/* 	for (j=0, t=s[i].; j<dim; j++) */
/* 	  for (l=0; l<dim; l++){ */
/* 	    if (j==l) t *= G_X[I[i*dim+j]+l*dim*dim]; */
/* 	    else t *= G[I[i*dim+j]]; */
/* 	g += t; */
/*       } */
/*   } */

  if (pg != NULL) (*pg) = g;

  return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_MeshMotionMap_Analytical
int 
xf_MeshMotionMap_Analytical( xf_AnaMotionsSet *AnaMotionsSet, int npoint, int dim, 
			     real Time, const real *qX, xf_MotionData *MData)
{
  int ierr;
  int k, n, i, ipoint;
  int ki, kj, kk, kij;
  int dim2, dim3;
  int nMotion, iMotion, nTerm;
  enum xfe_Bool NeedG, NeedG_X;
  enum xfe_Bool OutsideBlend;
  real xprev[3], Gprev[9], Gprev_X[27];
  real G0[9], G0_X[27];
  real H[9], H_xprev[27];
  real x_t[3], G2[9], Ttemp[27], vtemp[3];
  real *x, *G, *G_X, *gig_X, *vg;
  real g, b, b_X[3], b_XX[9];
  real txp, tX;
  const real *X;
  xf_AnaMotion *AnaMotion = NULL;
  xf_AnaMotions *AnaMotions = NULL;
  xf_AnaInData InData;
  xf_AnaOutData OutData;

  // error if no motions
  if ((nMotion = AnaMotionsSet->nMotion) <= 0) return xf_Error(xf_INPUT_ERROR);

  // convenient variables
  dim2 = dim*dim;
  dim3 = dim*dim2;

  // what will we need in terms of derivatives?
  NeedG_X = (MData->gbigb_X != NULL);
  NeedG   = ((MData->G != NULL) || NeedG_X);

  // set quantities to zero/one/identity
  if (MData->vg    != NULL) for (k=0; k<npoint*dim;     k++) MData->vg[k]    = 0.;
  if (MData->G     != NULL) 
    for (ipoint=0; ipoint<npoint; ipoint++){ // all G's set to identity
      for (k=0; k<dim2; k++) MData->G[ipoint*dim2+k]        = 0.;
      for (k=0; k<dim2; k+=(dim+1)) MData->G[ipoint*dim2+k] = 1.;
    }
  if (MData->g       != NULL) for (k=0; k<npoint;         k++) MData->g[k]       = 1.;
  if (MData->gbigb_X != NULL) for (k=0; k<npoint*dim;     k++) MData->gbigb_X[k] = 0.;

  // initialize MData->x to X
  if (MData->x != NULL) for (k=0; k<npoint*dim; k++) MData->x[k] = qX[k];
  else return xf_Error(xf_INPUT_ERROR);


  // loop over points
  for (ipoint=0; ipoint<npoint; ipoint++){
    
    // Set pointers to x, G, g, gig_X, vg
    x     = MData->x + ipoint*dim;
    G     = ((NeedG) ? G0 : NULL);
    G_X   = ((NeedG_X) ? G0_X : NULL);
    gig_X = ((MData->gbigb_X != NULL) ? MData->gbigb_X + ipoint*dim : NULL);
    vg    = ((MData->vg != NULL) ? MData->vg + ipoint*dim : NULL);
    
    // initialize G to identity, G_X to 0, g to 1, gig_X to 0, vg to 0
    if (NeedG){
      for (k=0; k<dim2; k++) G[k] = 0.;
      for (k=0; k<dim2; k+=(dim+1)) G[k] = 1.;
    }
    if (NeedG_X) for (k=0; k<dim3; k++) G_X[k] = 0.;
    if (vg  != NULL) for (k=0; k<dim; k++) vg[k] = 0.;
    if (gig_X != NULL) for (k=0; k<dim; k++) gig_X[k] = 0.;
    

    /* Begin loop over motions.  The motions are assumed independent in
       that the blending regions of the motions are disjoint.  So we
       just perform a simple superposition of the mappings. */
    
    for (iMotion=0; iMotion<nMotion; iMotion++){

      AnaMotions = AnaMotionsSet->AnaMotions + iMotion;

      // error if no terms
      if ((nTerm = AnaMotions->nTerm) <= 0) return xf_Error(xf_INPUT_ERROR);
      
      /* before we do anything, check if outside of blending region --
	 means this mesh motion is not active.  Perform this check
	 while computing blending function and derivatives. */
      ierr = xf_Error(xf_BlendMotion(AnaMotions, dim, qX+ipoint*dim,
				     &b, b_X, (NeedG_X) ? b_XX : NULL, &OutsideBlend));
      if (ierr != xf_OK) return ierr;
      if (OutsideBlend) continue; // nothing to do, outside blending

      /*
        Recursion formulas for combinations of motions:
      
	x{i} = x{i}(x{i-1}( ... (x{0}))), where x{0} = X = original ref coord
	
	G{i} = x{i}_X = H{i}*G{i-1}, where H{i} = x{i}_x{i-1}

	G{i}_X = H{i}_X*G{i-1} + H{i}*G{i-1}_X
               = H{i}_x{i-1}*G{i-1}*G{i-1} + H{i}*G{i-1}_X

        Also, all x{i} are potentially functions of t
	vg{i} = x{i}_x{i-1}*vg{i-1} + x{i}_t
      */


      // loop over analytical motions
      for (i=1; i<=nTerm; i++){
      
	AnaMotion = AnaMotions->AnaMotion + i-1;

	// At this point have x_{i-1}, H{i-1}, G{i-1}, G{i-1}_X, vg{i-1}

	// set xprev, Gprev, Gprev_X
	for (k=0; k<dim; k++) xprev[k] = x[k];
	if (NeedG  ) for (k=0; k<dim2; k++) Gprev[k] = G[k];
	if (NeedG_X) for (k=0; k<dim3; k++) Gprev_X[k] = G_X[k];

	// Call motion routine with x_{i-1} to obtain x{i}, H{i}, h{i}, h{i}_x{i-1}, x{i}_t
	InData.dim  = dim;
	InData.Time = Time;
	InData.X    = xprev;
	OutData.x   = x;
	OutData.H   = ((NeedG  ) ? H   : NULL);
	OutData.H_X = ((NeedG_X) ? H_xprev : NULL);
	OutData.x_t = ((vg != NULL) ? x_t : NULL);
	ierr = xf_Error(xf_ApplyAnaMotion(AnaMotion, &InData, &OutData));
	if (ierr != xf_OK) return ierr; 
      
	// G{i} = H{i}*Gprev (matrix multiplication)
	if (NeedG) xf_dimMxM(H, Gprev, dim, G);

	// G{i}_X = H{i}_X*G{i-1} + H{i}*G{i-1}_X
	//        = H{i}_x{i-1}*G{i-1}*G{i-1} + H{i}*G{i-1}_X
	if (NeedG_X){
	
	  // G{i}_X = H{i}*G{i-1}_X
	  for (k=0; k<dim; k++) xf_dimMxM(H, Gprev_X+k*dim2, dim, G_X+k*dim2);
	  // G2 = G{i-1}*G{i-1}
	  xf_dimMxM(Gprev, Gprev, dim, G2);
	  // Ttemp = H{i}_x{i-1}*G2
	  for (k=0; k<dim; k++) xf_dimMxM(H_xprev+k*dim2, G2, dim, Ttemp+k*dim2);
	  // G_X += Ttemp
	  for (k=0; k<dim3; k++) G_X[k] += Ttemp[k];
	}
      
	// vg{i} = H{i}*vg{i-1} + x{i}_t
	if (vg != NULL){
	  xf_MxV(H, vg, dim, dim, xfe_Set, vtemp);
	  for (k=0; k<dim; k++) vg[k] = vtemp[k] + x_t[k];
	}
    
      } // iTerm
  
      // set xprev, Gprev, Gprev_X
      for (k=0; k<dim; k++) xprev[k] = x[k];
      if (NeedG  ) for (k=0; k<dim2; k++) Gprev[k] = G[k];
      if (NeedG_X) for (k=0; k<dim3; k++) Gprev_X[k] = G_X[k];
      X = qX + dim*ipoint;

      // x = b*X + (1-b)*x
      for (k=0; k<dim; k++) x[k] = b*X[k] + (1.-b)*xprev[k];

      // G = x_X = X*b_X + b*I - x*b_X + (1-b)*x_X
      /* 	for (ki=0,k=0; ki<dim; ki++){ */
      /* 	  for (kj=0,txp=xprev[ki],tX=X[ki]; kj<dim; kj++,k++) */
      /* 	    G[k] = Gprev[k]*(1.-b) - txp*b_X[kj] + tX*b_X[kj]; */
      /* 	  G[ki*(dim+1)] += b; */
      /* 	} */
      if (NeedG){
	if (dim == 2){ // unrolled for speed
	  G[0] = Gprev[0]*(1.-b) + b - xprev[0]*b_X[0] + X[0]*b_X[0];
	  G[1] = Gprev[1]*(1.-b)     - xprev[0]*b_X[1] + X[0]*b_X[1];
	  G[2] = Gprev[2]*(1.-b)     - xprev[1]*b_X[0] + X[1]*b_X[0];
	  G[3] = Gprev[3]*(1.-b) + b - xprev[1]*b_X[1] + X[1]*b_X[1];
	}
	else{  // unrolled for speed
	  G[0] = Gprev[0]*(1.-b) + b - xprev[0]*b_X[0] + X[0]*b_X[0];
	  G[1] = Gprev[1]*(1.-b)     - xprev[0]*b_X[1] + X[0]*b_X[1];
	  G[2] = Gprev[2]*(1.-b)     - xprev[0]*b_X[2] + X[0]*b_X[2];
	  G[3] = Gprev[3]*(1.-b)     - xprev[1]*b_X[0] + X[1]*b_X[0];
	  G[4] = Gprev[4]*(1.-b) + b - xprev[1]*b_X[1] + X[1]*b_X[1];
	  G[5] = Gprev[5]*(1.-b)     - xprev[1]*b_X[2] + X[1]*b_X[2];
	  G[6] = Gprev[6]*(1.-b)     - xprev[2]*b_X[0] + X[2]*b_X[0];
	  G[7] = Gprev[7]*(1.-b)     - xprev[2]*b_X[1] + X[2]*b_X[1];
	  G[8] = Gprev[8]*(1.-b) + b - xprev[2]*b_X[2] + X[2]*b_X[2];
	}
      }


      // G_X = x_XX' = X_X'*b_X + X*b_XX' + b_X'*I
      //               - x_X'*b_X - x*b_XX'
      //               - b_X'*x_X + (1-b)*x_XX'
      if (NeedG_X)
	for (ki=0; ki<dim; ki++){
	  for (kj=0; kj<dim; kj++){
	    for (kk=0; kk<dim; kk++){
	      kij = ki*dim+kj;
	      k = kij + kk*dim2;
	      G_X[k] = (1.-b)*Gprev_X[k];
	      G_X[k] -= Gprev[kij]*b_X[kk];
	      G_X[k] -= xprev[ki]*b_XX[kj*dim+kk];
	      G_X[k] -= Gprev[ki*dim+kk]*b_X[kj];
	      if (ki == kj) G_X[k] += b_X[kk];
	      G_X[k] += X[ki] * b_XX[kj*dim+kk];
	      if (ki == kk) G_X[k] += b_X[kj];
	    }
	  }
	}

      // adjust velocity, vg *= (1-b)
      if (vg != NULL) for (k=0; k<dim; k++) vg[k] *= (1.-b);

    
      // calculate determinant and derivatives
      if (NeedG){
	ierr = xf_Error(xf_MapDet(G, G_X, dim, &g, gig_X));
	if (ierr != xf_OK) return ierr;
      }
    
      // set G, g in output
      if (MData->G != NULL) for (k=0; k<dim2; k++) MData->G[ipoint*dim2+k] = G0[k];
      if (MData->g != NULL) MData->g[ipoint] = g;

    } // iMotion

  } // ipoint

  return xf_OK;
}


#if( UNIT_TEST==1 )
#include "xf_MeshMotionAnalytical.test.in"
#endif




