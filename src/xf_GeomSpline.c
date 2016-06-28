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
  FILE:  xf_GeomSpline.c

  This file contains functions for working with spline geometry
  components.

*/


#include "xf_AllStruct.h"
#include "xf_Memory.h"
#include "xf_MPI.h"
#include "xf_IO.h"
#include "xf_String.h"
#include "xf_Mesh.h"
#include "xf_Param.h"
#include "xf_Math.h"
#include "xf_QuadRule.h"


/* Certain functions or portions of functions in this file were
   translated by K. Fidkowski from Fortran code written by Mark Drela:
   spline.f, (C) 2000 Mark Drela. */

/******************************************************************/
//  FUNCTION Definition: seval
static real 
seval(real SS, real *X, real *XS, real *S, int N)
{
/*

PURPOSE:  Evaluates X(SS).  XS must have been already calculated.

INPUTS
  
  SS: s-value at which to evaluate
  X, XS, S: spline info
  N: number of splines

OUTPUTS: 

RETURNS:

  Evaluated spline

*/
  int ilow, i, imid;
  real t, cx1, cx2, ds, val;
  
  ilow = 0; // 1
  i = N-1; // N

  while ((i-ilow) > 1){

    imid = (i+ilow)/2; // truncates decimal
    if (SS < S[imid])
      i = imid;
    else
      ilow = imid;
  }

  ds = S[i] - S[i-1];
  t = (SS - S[i-1]) / ds;
  cx1 = ds*XS[i-1] - X[i] + X[i-1];
  cx2 = ds*XS[i]   - X[i] + X[i-1];

  val = t*X[i] + (1.0-t)*X[i-1] + (t-t*t)*((1.0-t)*cx1 - t*cx2);

  return val;
} // end seval


/******************************************************************/
//  FUNCTION Definition: teval
static real 
teval(real SS, real *X, real *XS, real *S, int N)
{
/*

PURPOSE:  Evaluates spline tangent, dX/dS

INPUTS:

  SS: s-value at which to evaluate
  X, XS, S: spline info
  N: number of splines

OUTPUTS: 

RETURNS:

  Evaluated spline tangent

*/

  int ilow, i, imid;
  real t, cx1, cx2, ds, val;
  
  ilow = 0; // 1
  i = N-1; // N

  while ((i-ilow) > 1){

    imid = (i+ilow)/2; // truncates decimal
    if (SS < S[imid])
      i = imid;
    else
      ilow = imid;
  }

  ds = S[i] - S[i-1];
  t = (SS - S[i-1]) / ds;
  cx1 = ds*XS[i-1] - X[i] + X[i-1];
  cx2 = ds*XS[i]   - X[i] + X[i-1];

  val = ( X[i]-X[i-1] + (1-4.0*t+3.0*t*t)*cx1 - (2.0*t-3.0*t*t)*cx2 )/ds;

  return val;
} // end teval

/******************************************************************/
//  FUNCTION Definition: xf_SplineEval
static void
xf_SplineEval(real SS, xf_GeomCompSpline *GCS, real *x)
{

  x[0] = seval(SS, GCS->X, GCS->XS, GCS->S, GCS->N);
  x[1] = seval(SS, GCS->Y, GCS->YS, 
	       GCS->S, GCS->N);
}

/******************************************************************/
//  FUNCTION Definition: xf_SplineTanEval
static void
xf_SplineTanEval(real SS, xf_GeomCompSpline *GCS, real *t)
{

  t[0] = teval(SS, GCS->X, GCS->XS, 
	       GCS->S, GCS->N);
  t[1] = teval(SS, GCS->Y, GCS->YS, 
	       GCS->S, GCS->N);
}

/******************************************************************/
//  FUNCTION Definition: trisol
static int 
trisol(real *A, real *B, real *C, real *D, int kk)
{
/*-----------------------------------------  */
/*    Solves KK long, tri-diagonal system  | */
/*                                         | */
/*             A C          D              | */
/*            B A C        D               | */
/*               B A .      .              | */
/*                 . . C    .              | */
/*                   B A    D              | */
/*                                         | */
/*     The righthand side D is replaced by | */
/*     the solution.  A, C are destroyed.  | */
/*-----------------------------------------  */

  int k, km;

  for (k=1; k<kk; k++){
    km = k-1;
    C[km] = C[km] / A[km];
    D[km] = D[km] / A[km];
    A[k] = A[k] - B[k]*C[km];
    D[k] = D[k] - B[k]*D[km];
  } // k

  D[kk-1] = D[kk-1]/A[kk-1];

  for (k=(kk-2); k>=0; k--){
    D[k] = D[k] - C[k]*D[k+1];
  }
  
  return(xf_OK);
} // end trisol


/******************************************************************/
//  FUNCTION Definition: splind
static int 
splind(real *X, real *XS, real *S, int N, real XS1, real XS2)
{
/* C------------------------------------------------------- */
/* C     Calculates spline coefficients for X(S).          | */
/* C     Specified 1st derivative and/or usual zero 2nd    | */
/* C     derivative end conditions are used.               | */
/* C     To evaluate the spline at some value of S,        | */
/* C     use SEVAL and/or DEVAL.                           | */
/* C                                                       | */
/* C     S        independent variable array (input)       | */
/* C     X        dependent variable array   (input)       | */
/* C     XS       dX/dS array                (calculated)  | */
/* C     N        number of points           (input)       | */
/* C     XS1,XS2  endpoint derivatives       (input)       | */
/* C              If = 999.0, then usual zero second       | */
/* C              derivative end condition(s) are used     | */
/* C              If = -999.0, then zero third             | */
/* C              derivative end condition(s) are used     | */
/* C                                                       | */
/* C------------------------------------------------------- */

  int ierr, i;
  real dsm, dsp;
  real *A, *B, *C;


  /* allocate A, B, C */
  ierr = xf_Error(xf_Alloc((void **) &A, 3*N, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  B = A+N; C = B+N;

  for (i=1; i<N-1; i++){
    dsm = S[i] - S[i-1];
    dsp = S[i+1] - S[i];
    B[i] = dsp;
    A[i] = 2.0*(dsm+dsp);
    C[i] = dsm;
    XS[i] = 3.0*((X[i+1]-X[i])*dsm/dsp + (X[i]-X[i-1])*dsp/dsm);
  }

  if (XS1 == 999.0){
    //----- set zero second derivative end condition
    A[0] = 2.0;
    C[0] = 1.0;
    XS[0] = 3.0*(X[1]-X[0]) / (S[1]-S[0]);
  }
  else if (XS1 == -999.0){
    //----- set zero third derivative end condition
    A[0] = 1.0;
    C[0] = 1.0;
    XS[0] = 2.0*(X[1]-X[0]) / (S[1]-S[0]);
  }
  else{
    //----- set specified first derivative end condition
    A[0] = 1.0;
    C[0] = 0.0;
    XS[0] = XS1;
  }

  if (XS2 == 999.0){
    B[N-1] = 1.0;
    A[N-1] = 2.0;
    XS[N-1] = 3.0*(X[N-1]-X[N-2]) / (S[N-1]-S[N-2]);
  }
  else if (XS2 == -999.0){
    B[N-1] = 1.0;
    A[N-1] = 1.0;
    XS[N-1] = 2.0*(X[N-1]-X[N-2]) / (S[N-1]-S[N-2]);
  }
  else{
    A[N-1] = 1.0;
    B[N-1] = 0.0;
    XS[N-1] = XS2;
  }
  
  if ( (N==2) && (XS1 == -999.0) && (XS2 == -999.0)){
    B[N-1] = 1.0;
    A[N-1] = 2.0;
    XS[N-1] = 3.0*(X[N-1]-X[N-2]) / (S[N-1]-S[N-2]);
  }

  //---- solve for derivative array XS
  ierr = xf_Error(trisol(A,B,C,XS,N));
  if (ierr != xf_OK) return ierr;
 		     
  xf_Release( (void *) A);
  
  return xf_OK;
} // end splind


/******************************************************************/
//  FUNCTION Definition: scalc
static int 
scalc(real *X, real *Y, real *S, int N)
{
/* C---------------------------------------- */
/* C     Calculates the arc length array S  | */
/* C     for a 2-D array of points (X,Y).   | */
/* C---------------------------------------- */

  int i;

  S[0] = 0.0;
  for (i=1; i<N; i++){
    S[i] = S[i-1] + sqrt( (X[i]-X[i-1])*(X[i]-X[i-1]) + (Y[i]-Y[i-1])*(Y[i]-Y[i-1]) );
  }

  return xf_OK;
} // end scalc


/******************************************************************/
//  FUNCTION Definition: segspl
static int 
segspl(real *X, real *XS, real *S, int N)
{
/* C----------------------------------------------- */
/* C     Splines X(S) array just like SPLINE,      | */
/* C     but allows derivative discontinuities     | */
/* C     at segment joints.  Segment joints are    | */
/* C     defined by identical successive S values. | */
/* C----------------------------------------------- */

  int ierr, iseg, iseg0, nseg;
  real t;
  real tol = 1e-13;

  if (fabs(S[0]-S[1]) < tol){
    printf("SEGSPL:  First input point duplicated\n");
    return xf_Error(xf_INPUT_ERROR);
  }

  if (fabs(S[N-1]-S[N-2]) < tol){
    printf("SEGSPL:  Last  input point duplicated\n");
    return xf_Error(xf_INPUT_ERROR);
  }


  iseg0 = 0;
  t = -999.0;
  for (iseg = 1; iseg<N-2; iseg++){
    if (fabs(S[iseg]-S[iseg+1]) < tol){
      nseg = iseg - iseg0 + 1;
      ierr = xf_Error( splind( X+iseg0, XS+iseg0, S+iseg0, nseg, t, t) );
      if (ierr != xf_OK) return ierr;
      iseg0 = iseg+1;
    }
  } // iseg;

  nseg = N - iseg0;
  ierr = xf_Error( splind( X+iseg0, XS+iseg0, S+iseg0, nseg, t, t) );
  if (ierr != xf_OK) return ierr;

  return xf_OK;
} // segspl
  

/******************************************************************/
//  FUNCTION Definition: splnxy
static int 
splnxy(real *X, real *XS, real *Y, real *YS, real *S, int N)
{
/* C-----------------------------------------  */
/* C     Splines 2-D shape X(S), Y(S), along | */
/* C     with true arc length parameter S.   | */
/* C-----------------------------------------  */

  int ierr, npass, ipass, i, k, nq;
  real ds, dx, dy, cx1, cx2, cy1, cy2;
  real t, sint;
  real *xq, *wq, s, xs, ys, *Snew;
  real serr;
 
  npass = 10;
  
  // allocate Snew
  ierr = xf_Error(xf_Alloc( (void **) &Snew, N, sizeof(real)));
  if (ierr != xf_OK) return ierr;
    
  // first estimate of arc length parameter
  ierr = xf_Error(scalc(X,Y,S,N));
  if (ierr != xf_OK) return ierr;

  // spline X(S) and Y(S)
  ierr = xf_Error(segspl(X, XS, S, N));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(segspl(Y, YS, S, N));
  if (ierr != xf_OK) return ierr;

  
  /* Get quadrature rule (order 10 seems high enough) */
  xq = (real *) NULL;
  wq = (real *) NULL;
  ierr = xf_Error( xf_QuadLine(10, &nq, &xq, &wq) );
  if (ierr != xf_OK) return ierr;
  

  //---- re-integrate true arc length
  ipass = 0;
  while (ipass < npass){

    serr = 0.0;

    Snew[0] = S[0];

    for (i=1; i<N; i++){
      dx = X[i] - X[i-1];
      dy = Y[i] - Y[i-1];
      ds = S[i] - S[i-1];

      cx1 = ds*XS[i-1] - dx;
      cx2 = ds*XS[i  ] - dx;
      cy1 = ds*YS[i-1] - dy;
      cy2 = ds*YS[i  ] - dy;

      /* integrate using quadrature */
      sint = 0.0;
      
      for (k=0; k<nq; k++){
	s = S[i-1] + (1.0+xq[k])*0.5*(S[i]-S[i-1]);
	xs = teval(s, X, XS, S, N);
	ys = teval(s, Y, YS, S, N);
	sint += 0.5*wq[k]*sqrt(xs*xs + ys*ys)*ds;
      } // k

      if (fabs(sint-ds) > fabs(serr)) serr = sint - ds;
      
      Snew[i] = Snew[i-1] + sint;
    } // i

    serr = serr / (Snew[N-1] - Snew[0]);

    /* set S = Snew */
    for (i=0; i<N; i++) S[i] = Snew[i];

    // re-spline X(S) and Y(S)
    ierr = xf_Error(segspl(X, XS, S, N));
    if (ierr != xf_OK) return ierr;
    
    ierr = xf_Error(segspl(Y, YS, S, N));
    if (ierr != xf_OK) return ierr;

    ipass++;
  } // ipass
  
  xf_Release( (void *) Snew);
  xf_Release( (void *) xq  );
  xf_Release( (void *) wq  );

  return xf_OK;

} // end splnxy



/******************************************************************/
//  FUNCTION Definition: xf_SplinePoints
static int 
xf_SplinePoints(xf_GeomCompSpline *GCS)
{
  int ierr, N;

  // only need to spline if S is NULL (else a spline already exists)
  if (GCS->S != NULL) return xf_OK;

  // number of points
  N = GCS->N;

  // allocate memory for S, XS, and YS
  ierr = xf_Error(xf_Alloc( (void **) &GCS->S , N, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &GCS->XS, N, sizeof(real)));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_Alloc( (void **) &GCS->YS, N, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // spline X and Y simultaneously
  ierr = xf_Error(splnxy(GCS->X, GCS->XS, GCS->Y, GCS->YS, GCS->S, GCS->N));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
} 


/******************************************************************/
//   FUNCTION Definition: xf_SegmentProject
static int 
xf_SegmentProject(real *p0, real *p1, real *q, 
		 real *pc, real *dist, real eps)
{

  real dp0, dp1, d0, d1, a2, s, pc2[2];

  // if (p0==p1) this is not a segment; return distance to p0
  if ((fabs(p0[0]-p1[0]) < eps) && (fabs(p0[1]-p1[1]) < eps)){
    d0 = sqrt( (q[0]-p0[0])*(q[0]-p0[0]) + (q[1]-p0[1])*(q[1]-p0[1]) );
    (*dist) = d0;
    if (pc != NULL){
      pc[0] = p0[0];
      pc[1] = p0[1];
    }
    return xf_OK;
  }

  dp0 = (q[0]-p0[0])*(p1[0]-p0[0]) + (q[1]-p0[1])*(p1[1]-p0[1]);
  dp1 = (q[0]-p1[0])*(p0[0]-p1[0]) + (q[1]-p1[1])*(p0[1]-p1[1]);

  d0 = sqrt( (q[0]-p0[0])*(q[0]-p0[0]) + (q[1]-p0[1])*(q[1]-p0[1]) );

  if (dp0*dp1 < 0){ // projection not on segment, take min dist to nodes
    d1 = sqrt( (q[0]-p1[0])*(q[0]-p1[0]) + (q[1]-p1[1])*(q[1]-p1[1]) );
    if (d0 < d1){
      (*dist) = d0;
      if (pc != NULL){
	pc[0] = p0[0];
	pc[1] = p0[1];
      }
    }
    else{
      (*dist) = d1;
      if (pc != NULL){
	pc[0] = p1[0];
	pc[1] = p1[1];
      }
    }
  }
  else{ // projection onto segment
    a2 = (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) ;
    s = dp0/a2;

    pc2[0] = p0[0]*(1-s) + p1[0]*s;
    pc2[1] = p0[1]*(1-s) + p1[1]*s;

    if (pc != NULL){
      pc[0] = pc2[0]; pc[1] = pc2[1];
    }
    (*dist) = sqrt((pc2[0]-q[0])*(pc2[0]-q[0]) + (pc2[1]-q[1])*(pc2[1]-q[1]));
  }
  
  return xf_OK;
}


/******************************************************************/
//  FUNCTION Definition: xf_ProjectToSpline
static int 
xf_ProjectToSpline(xf_GeomCompSpline *GCS, real *xn, real *sn)
{
/*
  PURPOSE: 

    Projects point xn to cubic spline in GCS

  INPUTS:

    GCS : spline geometry component
    xn : point to projecdt
  
  OUTPUTS: 

    xn : modified (shifted) point
    sn : s-value on spline (optional)

  RETURNS:  Error code

*/

  int ierr, i, i_min, N;
  int k, nb, flag0[2], flag1[2];
  int knot_flag, knmax, kn, knmin, i0, i1;
  enum xfe_Bool flag;
  real *X, *Y, *S, *XS, *YS;
  real p0[2], p1[2], p[2];
  real t[2], v[2], mag, dot;
  real kdist[3], kp[3][2], ks[3];
  real eps, dist, dist_min;
  real d0, d1, d01, s, s0, s1, ss;

  eps = 1e-12;

  // calculate spline if have not yet done so
  if (GCS->S == NULL){
    ierr = xf_Error(xf_SplinePoints(GCS));
    if (ierr != xf_OK) return ierr;
  }
  
  N  = GCS->N;
  X  = GCS->X;
  Y  = GCS->Y;
  S  = GCS->S;
  XS = GCS->XS;
  YS = GCS->YS;

  
  flag = xfe_True;
  i_min = -1;
  dist_min = -1.0;
  /* Loop over spline knots.  Determine knot segment on spline that is
     closest to xn.  This is done by distance checks to the linear
     knot segments.  If the minimal distance is to a knot rather than
     a segment, knot_flag is set to true. */
  for (i=0; i<(N-1); i++){
    p0[0] = X[i]  ;   p0[1] = Y[i];
    p1[0] = X[i+1];   p1[1] = Y[i+1];

    /* if p0 == p1 (i.e. corner), continue */
    d01 =  sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) );
    if (d01 < eps)
      continue;

    /* check to see if xn is on one of knots -- if so, return */
    d0  = sqrt( (xn[0]-p0[0])*(xn[0]-p0[0]) + (xn[1]-p0[1])*(xn[1]-p0[1]) );
    d1  = sqrt( (xn[0]-p1[0])*(xn[0]-p1[0]) + (xn[1]-p1[1])*(xn[1]-p1[1]) );
    
    if (d0 < eps){
      if (sn != NULL) *sn = S[i];
      return xf_OK;
    }

    if (d1 < eps){
      if (sn != NULL) *sn = S[i+1];
      return xf_OK;
    }

    ierr = xf_Error(xf_SegmentProject(p0, p1, xn, NULL, &dist, eps));
    if (ierr != xf_OK)
      return ierr;
    
    if ((flag) || (dist < dist_min)){
      dist_min = dist;
      i_min = i;
      flag = xfe_False;
    }
  } // i

  knmax = 3;
  dist_min = 1e30;
  knmin = -1;
  for (kn=0; kn<knmax; kn++){

    i0 = i_min-1 + kn;
    i1 = i0 + 1;

    if ((i0 < 0) || (i1 > (N-1))){
      kdist[kn] = 1e30;
      continue;
    }

    s0 = S[i0];
    s1 = S[i1];
    

    p0[0] = X[i0];   p0[1] = Y[i0];
    p1[0] = X[i1];   p1[1] = Y[i1];

    if (fabs(s0-s1) < MEPS){
      d0 = sqrt( ( xn[0]-p0[0])*( xn[0]-p0[0]) + ( xn[1]-p0[1])*( xn[1]-p0[1]) );
      kdist[kn] = d0;
      continue;
    }


    nb = 25; // number of bisection iterations
    flag0[kn] = flag1[kn] = xfe_False;

    for (k=0; k<nb; k++){  // determine closest s via smart bisection
      s = 0.5*(s0+s1);

      xf_SplineEval(s, GCS, p);

      // get unit tangent vector at s
      xf_SplineTanEval(s, GCS, t);
      mag = sqrt(t[0]*t[0] + t[1]*t[1]);
      t[0] = t[0]/mag; t[1] = t[1]/mag;
     

      // get unit vector pointing from p to xn
      v[0] = xn[0]-p[0]; v[1] = xn[1]-p[1];
      mag = sqrt(v[0]*v[0] + v[1]*v[1]);
      
      if (mag < eps){ // xn is on spline and we hit it
	if (sn != NULL) *sn = s;
	return xf_OK;
      }
     
      v[0] = v[0]/mag; v[1] = v[1]/mag;

      dot = t[0]*v[0] + t[1]*v[1];
      
      if (dot < 0){
	s1 = s;
	p1[0] = p[0]; p1[1] = p[1];
	flag0[kn] = xfe_True;
      }
      else if (dot > 0){
	s0 = s;
	p0[0] = p[0]; p0[1] = p[1];
	flag1[kn] = xfe_True;
      }
      else{ // dot == 0, meaning p is projection, exit
	if (sn != NULL) *sn = s;
	xn[0] = p[0]; xn[1] = p[1];
	return xf_OK;
      }
      
    } // end for k

    ierr = xf_Error(xf_SegmentProject(p0, p1, xn, p, &dist, eps));
    if (ierr != xf_OK)
      return ierr;
    
    d01 = sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) );
    d0  = sqrt( ( p[0]-p0[0])*( p[0]-p0[0]) + ( p[1]-p0[1])*( p[1]-p0[1]) );
    ss = 0;
    if (d0 < eps)
      ss = 0;
    else
      ss = d0/d01;
    s = s0 + ss*(s1-s0);
    
    xf_SplineEval(s, GCS, p);

    kdist[kn] = dist;
    kp[kn][0] = p[0];
    kp[kn][1] = p[1];
    ks[kn] = s;

    if (dist < dist_min){
      dist_min = dist;
      knmin = kn;
    }
    
  } // kn
  
  if (!flag1[knmin]){
    xn[0] = X[i_min-1+knmin];
    xn[1] = Y[i_min-1+knmin];
    if (sn != NULL) *sn = S[i_min-1+knmin];
    return xf_OK;
  }
  if (!flag0[knmin]){
    xn[0] = X[i_min+knmin];
    xn[1] = Y[i_min+knmin];
    if (sn != NULL) *sn = S[i_min+knmin];
    return xf_OK;
  }
  xn[0] = kp[knmin][0];
  xn[1] = kp[knmin][1];
  if (sn != NULL) *sn = ks[knmin];
  
    return xf_OK;
}


/******************************************************************/
//   FUNCTION Definition: xf_ProjectToGeomComp_Spline
int 
xf_ProjectToGeomComp_Spline( xf_GeomCompSpline *GCS, int dim, real *x)
{

  if (dim != 2) return xf_Error(xf_INPUT_ERROR);

  if (GCS->Order == 3)
    return xf_ProjectToSpline(GCS, x, NULL);
  else
    return xf_Error(xf_NOT_SUPPORTED);


  return xf_OK;
}



/******************************************************************/
//   FUNCTION Definition: xf_PointsOnGeom_Spline
int 
xf_PointsOnGeom_Spline( xf_GeomCompSpline *GCS, int dim, 
			int np, real dl, enum xfe_GeomSpacingType Spacing, 
			int *pnout, real **px)
{
  int ierr;
  int i, nout;

/*   // if points not yet splined, do so now */
/*   if (GCS->S == NULL){ */
/*     ierr = xf_Error(xf_SplinePoints(GCS)); */
/*     if (ierr != xf_OK) return ierr; */
/*   } */

/*   // how many points do we need? */
/*   if (np > 0) nout = np; */
/*   else nout = S[GCS->N-1]/dl; */
/*   (*pnout) = nout; */

  // for now, just return input spline points, ignoring spacing requests
  (*pnout) = nout = GCS->N-1;  

  // allocate memory for points
  ierr = xf_Error(xf_Alloc( (void **) px, 2*nout, sizeof(real)));
  if (ierr != xf_OK) return ierr;

  // fill in points (assume closed spline)
  for (i=0; i<2*nout; i++) (*px)[i] = GCS->X[i];

  return xf_OK;
}

