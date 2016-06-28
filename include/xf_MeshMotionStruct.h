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

#ifndef _xf_MeshMotionStruct_h
#define _xf_MeshMotionStruct_h 1

/*
  FILE:  xf_MeshMotionStruct.h

  This file contains the xflow mesh motion structures

*/

#include "xf.h"
#include "xf_DataStruct.h"
#include "xf_ParamStruct.h"


//unified name for GCL vector
#define GCLVectorTitle "GCLVectorMain"
#define GCLAdjointTitle "GCLAdjointMain"

/*------------------------*/
/* ANALYTICAL MESH MOTION */
/*------------------------*/

/* Analytical motion types */
enum xfe_AnaMotionType{
  xfe_AnaMotion_Plunge,    /* Plunge  */
  xfe_AnaMotion_Pitch,     /* Pitch   */
  xfe_AnaMotion_PerssonFS, /* Persson's FStream Sinusoidal Mapping */
  xfe_AnaMotionLast
};
/* corresponding names */
static char *xfe_AnaMotionName[xfe_AnaMotionLast] = {
  "Plunge",
  "Pitch",
  "PerssonFS"
};


/* Analytical blending functions */
enum xfe_AnaBlendType{
  xfe_AnaBlend_None,        /* No blending */
  xfe_AnaBlend_Cubic,       /* Cubic blending */
  xfe_AnaBlend_Quintic,     /* Quintic blending */
  xfe_AnaBlend_Septic,      /* Septic (7th order) blending */
  xfe_AnaBlendLast
};
/* corresponding names */
static char *xfe_AnaBlendName[xfe_AnaBlendLast] = {
  "None",
  "Cubic",
  "Quintic",
  "Septic"
};

/* Structure for one expression in analytical motion; no blending here */
typedef struct
{
  char *Name;
  /* An arbitrary (descriptive) name given to this motion */

  enum xfe_AnaMotionType MotionType;
  /* Type of analytical motion */

  xf_KeyValue MotionKeyValue;
  /* Keys and Values understood by the motion function */

  real *MotionRParam;
  int  *MotionIParam;
  /* real/integer key values for quick access; not read or written */

}
xf_AnaMotion;


/* Analytical motions structure; blending here */
typedef struct
{
  int nTerm;
  /* Multiple analytical expressions can be superimposed */

  xf_AnaMotion *AnaMotion;
  /* Vector of analytical expressions */

  enum xfe_AnaBlendType BlendType;
  /* Blending function */
  
  xf_KeyValue BlendKeyValue;
  /* Keys and Values understood by the blending function */
  
  real *BlendRParam;
  int  *BlendIParam;
  /* real/integer key values for quick access; not read or written */

}
xf_AnaMotions;

/* Structure for set of analytical motions: can have more than one if
   have objects moving independently, each one requiring own blending
   function. */
typedef struct
{
  int nMotion;
  /* Multiple analytical expressions can be superimposed */

  xf_AnaMotions *AnaMotions;
  /* Vector of analytical motions, each one with own blending */

}
xf_AnaMotionsSet;


/*-------------------------------*/
/* GENERAL MeshMotion STRUCTURES */
/*-------------------------------*/


/* General (high-level) motion types.  For example, analytical versus
   spring analogy. */
enum xfe_MotionType{
  xfe_Motion_Analytical,    /* Prescribed analytical function */
  xfe_MotionLast,
};
/* corresponding names */
static char *xfe_MotionName[xfe_MotionLast] = {
  "Analytical"
};


/* MeshMotion structure definition */
typedef struct
{   
  enum xfe_MotionType Type; 
  /* Type of mesh motion (e.g. analytical) */

  void *Data;
  /* Data associated with this mesh motion */  

  enum xfe_Bool Active;
  /* Indicates whether the motion is turned on */
}
xf_MeshMotion;

/* Mesh Motion map data structure */
typedef struct
{

  int npoint;  /* number of points */
  
  int dim;     /* spatial dimension */
  
  real *x;     /* mapped coordinates in time-varying domain [npoint*dim] */
  
  real *vg;    /* velocity of mapping = dx/dt [npoint*dim] */

  real *G;     /* Jacobian of mapping = x_X [npoint*dim*dim] */
  
  real *g;     /* determinant of G [npoint] */

  real *gb;    /* gbar = interpolated GCL vector; points to g if GCL is off */
  
  real *gbigb_X; /* gb^{-1} * interpolated derivative of gb w.r.t. X [npoint*dim];
		    when GCL is off, g is used in place of gb */

  real *Ginv;  /* inverse of Jacobian [npoint*dim*dim] */

  xf_Vector *GCLVector;  /* approximated gbar vector if using GCL */
  
}
xf_MotionData;


// bit masks for specifying what we want allocated in mesh motion data
#define xfb_MD_x        1
#define xfb_MD_vg       2
#define xfb_MD_G        4
#define xfb_MD_g        8
#define xfb_MD_gb       16
#define xfb_MD_gbigb_X  32
#define xfb_MD_Ginv     64


#endif // end ifndef _xf_MeshMotionStruct_h
