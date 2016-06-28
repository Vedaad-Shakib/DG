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


#include "xf_Unit.h"
#include "xf_All.h"

/* TODO: add a utility function for DefaultAnaMotionsSet, with
   possibly independent mesh motions */


// utility function for creating a motion
static int
xf_DefaultAnaMotions(xf_AnaMotions **pMotions, enum xfe_AnaMotionType MotionType,
		     enum xfe_AnaBlendType BlendType)
{
  // TODO: support test of multiple motions
  int ierr;
  char *KVM_Plunge[] = {"XAmplitude", "1.0",
			"YAmplitude", "0.5",
			"ZAmplitude", "0.3",
			"Frequency", "7.0",
			"Phase", "0.1",
			"\0"};

  char *KVM_Pitch[] = {"XOrigin", "10.0", 
                       "YOrigin", "5.2",
                       "ZOrigin", "2.3",
                       "XAmplitude", "1.0",
                       "YAmplitude", "0.25",
                       "ZAmplitude", "0.45",
                       "Frequency", "1.5",
                       "Phase", "0.4",
                       "TransientTime", "0.0",
                       "\0"};
  char *KVM_PerssonFS[] = {"To", "1.0", "\0"};
  char *KVB[] = {"XCenter", "-1.7",
		"YCenter", "1.4",
		"ZCenter", "0.1",
		"Radius", "0.3",
		"Dist", "7.0",
		"\0"};
  char **KVM = NULL;
  xf_AnaMotions *Motions;

  // create/fill AnaMotions structure
  ierr = xf_Error(xf_CreateAnaMotions(pMotions));
  xf_AssertEqual(ierr, xf_OK);
  Motions = (*pMotions);
  ierr = xf_Error(xf_AllocAnaMotions(Motions, 1));
  xf_AssertEqual(ierr, xf_OK);
  Motions->AnaMotion[0].MotionType = MotionType;
  Motions->BlendType = BlendType;
  
  // choose motion type
  switch(MotionType){
  case xfe_AnaMotion_Plunge:
    KVM = KVM_Plunge;
    break;
  case xfe_AnaMotion_Pitch:
    KVM = KVM_Pitch;
    break;
  case xfe_AnaMotion_PerssonFS:
    KVM = KVM_PerssonFS;
    break;
  default:
    return xf_Error(xf_NOT_SUPPORTED);
    break;
  }

  // choose blend type (just one for now, already set)
 
  // set key values
  ierr = xf_Error(xf_AddKeyValueList(&Motions->AnaMotion[0].MotionKeyValue, 
				     KVM, xfe_False, xfe_False));
  if (ierr != xf_OK) return ierr;
  ierr = xf_Error(xf_AddKeyValueList(&Motions->BlendKeyValue, 
				     KVB, xfe_False, xfe_False));
  if (ierr != xf_OK) return ierr;

  return xf_OK;
}


// utility function for testing analytical mesh motion that has already been defined
static int
xf_TestApplyAnaMotion(xf_AnaMotion *AnaMotion, xf_AnaInData InData)
{
  int ierr, i, j, dim;
  real x0[3], H0[9], H0_X[27], x0_t[3];
  real x[3], H[9], H_X[27], x_t[3];
  xf_AnaOutData OutData;
  real fd, an;
  real eps = 1e-5, tol = 1e-8;
  
  dim = InData.dim;

  // set OutData variables
  OutData.x = x0; OutData.H = H0; OutData.H_X = H0_X; OutData.x_t = x0_t;

  // call mesh motion routine
  ierr = xf_Error(xf_ApplyAnaMotion(AnaMotion, &InData, &OutData));
  if (ierr != xf_OK) return ierr; 

  // ping X derivatives
  for (j=0; j<dim; j++){

    // perturb X
    InData.X[j] += eps;

    // call motion routine
    OutData.x = x; OutData.H = H; OutData.H_X = H_X; OutData.x_t = x_t;
    ierr = xf_Error(xf_ApplyAnaMotion(AnaMotion, &InData, &OutData));
    if (ierr != xf_OK) return ierr; 

    // check derivatives against finite difference
    // H
    for (i=0; i<dim; i++){
      fd = (x[i] - x0[i])/eps;
      an = 0.5*(H0[i*dim+j] + H[i*dim+j]);
      if (fabs(fd-an) > tol){
	xf_printf("H ping fail[i,j=%d,%d]: fd = %.15E, an = %.15E, tol = %.5E\n", 
		  i,j, fd, an, tol);
	return xf_Error(xf_PING_FAILED);
      } 
    } // j
    // H_X
    for (i=0; i<dim*dim; i++){
      fd = (H[i] - H0[i])/eps;
      an = 0.5*(H0_X[i+j*dim*dim] + H_X[i+j*dim*dim]);
      if (fabs(fd-an) > tol){
	xf_printf("H_X ping fail[i,j=%d,%d]: fd = %.15E, an = %.15E, tol = %.5E\n", 
		  i,j, fd, an, tol);
	return xf_Error(xf_PING_FAILED);
      } 
    } // i

    // unperturb X
    InData.X[j] -= eps;

  } // j

  // ping Time derivative
  InData.Time += eps;
  OutData.x = x; OutData.H = H; OutData.H_X = H_X; OutData.x_t = x_t;
  ierr = xf_Error(xf_ApplyAnaMotion(AnaMotion, &InData, &OutData));
  if (ierr != xf_OK) return ierr; 
  for (j=0; j<dim; j++){
    fd = (x[j] - x0[j])/eps;
    an = 0.5*(x0_t[j] + x_t[j]);
    if (fabs(fd-an) > tol){
      xf_printf("time ping fail[j=%d]: fd = %.15E, an = %.15E, tol = %.5E\n", 
		  j, fd, an, tol);
      return xf_Error(xf_PING_FAILED);
    } 
  } // j

  return xf_OK;
}


TEST_xf_ApplyAnaMotion_Plunge()
{
  int ierr;
  xf_AnaMotions *Motions;
  xf_AnaInData InData;
  real X[3] = {0.1, 0.3, 0.25}, x[3];


  // create Motions structure
  ierr = xf_Error(xf_DefaultAnaMotions(&Motions, xfe_AnaMotion_Plunge,
				       xfe_AnaBlend_None));
  xf_AssertEqual(ierr, xf_OK);
  
  // call motion testing routine
  InData.dim  = 3;
  InData.Time = 1.2;
  InData.X    = X;
  ierr = xf_Error(xf_TestApplyAnaMotion(Motions->AnaMotion+0, InData));
  xf_AssertEqual(ierr, xf_OK);
  
  // destroy AnaMotions
  ierr = xf_Error(xf_DestroyAnaMotions(Motions, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;

}

TEST_xf_ApplyAnaMotion_Pitch()
{
  int ierr;
  xf_AnaMotions *Motions;
  xf_AnaInData InData;
  real X[3] = {0.1, 0.3, 0.25}, x[3];


  // create Motions structure
  ierr = xf_Error(xf_DefaultAnaMotions(&Motions, xfe_AnaMotion_Pitch,
				       xfe_AnaBlend_None));
  xf_AssertEqual(ierr, xf_OK);
  
  // call motion testing routine
  InData.dim  = 3;
  InData.Time = 1.2;
  InData.X    = X;
  ierr = xf_Error(xf_TestApplyAnaMotion(Motions->AnaMotion+0, InData));
  xf_AssertEqual(ierr, xf_OK);
  
  // destroy AnaMotions
  ierr = xf_Error(xf_DestroyAnaMotions(Motions, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;

}

TEST_xf_ApplyAnaMotion_PerssonFS()
{
  int ierr;
  xf_AnaMotions *Motions;
  xf_AnaInData InData;
  real X[3] = {0.1, 0.3, 0.0}, x[3];


  // create Motions structure
  ierr = xf_Error(xf_DefaultAnaMotions(&Motions, xfe_AnaMotion_PerssonFS,
				       xfe_AnaBlend_None));
  xf_AssertEqual(ierr, xf_OK);
  
  // call motion testing routine
  InData.dim  = 2;
  InData.Time = .4;
  InData.X    = X;
  ierr = xf_Error(xf_TestApplyAnaMotion(Motions->AnaMotion+0, InData));
  xf_AssertEqual(ierr, xf_OK);
  
  // destroy AnaMotions
  ierr = xf_Error(xf_DestroyAnaMotions(Motions, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
  
  return xf_OK;

}

TEST_xf_BlendMotion()
{
  int ierr, i, j, dim;
  xf_AnaMotions *Motions;
  enum xfe_AnaBlendType BlendType;
  enum xfe_AnaMotionType MType;
  real X[3] = {0.1, 0.3, 0.25}, x[3];
  real b, b_X[3], b_XX[9];
  real b0, b0_X[3], b0_XX[9];
  real fd, an, eps = 1e-5, tol = 1e-8;

  for (MType=0; MType<xfe_AnaMotionLast; MType++){

  // create Motions structure
  ierr = xf_Error(xf_DefaultAnaMotions(&Motions, MType,
				       xfe_AnaBlend_Cubic));
  xf_AssertEqual(ierr, xf_OK);
  
  for (BlendType=0; BlendType<xfe_AnaBlendLast; BlendType++){
    
    Motions->BlendType = BlendType;
    dim = 3;

    // call blend routine
    ierr = xf_Error(xf_BlendMotion(Motions, dim, X, &b0, b0_X, b0_XX, NULL));
    xf_AssertEqual(ierr, xf_OK);

    // ping X derivatives
    for (j=0; j<dim; j++){

      // perturb X
      X[j] += eps;

      // call blend routine
      ierr = xf_Error(xf_BlendMotion(Motions, dim, X, &b, b_X, b_XX, NULL));
      xf_AssertEqual(ierr, xf_OK);

      // b_X
      fd = (b - b0)/eps;
      an = 0.5*(b0_X[j] + b_X[j]);
      if (fabs(fd-an) > tol){
	xf_printf("b_X ping fail[j=%d]: fd = %.15E, an = %.15E, tol = %.5E\n", 
		  j, fd, an, tol);
	return xf_Error(xf_PING_FAILED);
      } // j
      // b_XX
      for (i=0; i<dim; i++){
	fd = (b_X[i] - b0_X[i])/eps;
	an = 0.5*(b0_XX[i*dim+j] + b_XX[i*dim+j]);
	if (fabs(fd-an) > tol){
	  xf_printf("b_XX ping fail[i,j=%d,%d]: fd = %.15E, an = %.15E, tol = %.5E\n", 
		    i,j, fd, an, tol);
	  return xf_Error(xf_PING_FAILED);
	} 
      } // i

      // unperturb X
      X[j] -= eps;

    } // j
  }
  
  // destroy AnaMotions
  ierr = xf_Error(xf_DestroyAnaMotions(Motions, xfe_True));
  xf_AssertEqual(ierr, xf_OK);
} //end loop over motion type
  
  return xf_OK;

}

TEST_xf_dimMxM()
{
  real A[9] = {1,2,5,4,3,6,7,9,8};
  real B[9] = {3,2,1,-1,-2,-3,-2,2,2};
  real C[9];
  real C0[9] = {-9,8,5,-3,14,7,-4,12,-4};

  xf_dimMxM(A,B,3,C);
  xf_AssertRealVectorWithin(C, C0, 9, UTOL0);

  return xf_OK;
}


TEST_xf_MapDet()
{
  int ierr, i, d, dim;
  real G2[4] = {1,5,4,7};
  real G2_X[8] = {1,2,3,4, 8,7,6,5};
  real G3[9] = {1,5,4, 7,-2,3, 8,9,6};
  real G3_X[27] = {1,3,2,4,5,7,6,8,9,
		   6,8,7,3,4,2,1,9,5,
		   7,6,4,5,3,2,8,1,9};
  real *G, *G_X;
  real g, gbigb_X[3];
  real g0, gbigb0_X[3];
  real fd, an, eps=1e-5, tol=1e-8;

  for (dim=2; dim<=3; dim++){

    if (dim==2){
      G = G2; G_X = G2_X;
    }
    else{
      G = G3; G_X = G3_X;
    }

    // ping gbigb_X
    ierr = xf_Error(xf_MapDet(G,G_X,dim,&g0,gbigb0_X));
    xf_AssertEqual(ierr, xf_OK);
    
    for (i=0; i<dim;i++){
      // use G_X to perturb G
      for (d=0; d<dim*dim; d++) G[d] += G_X[i*dim*dim+d]*eps;
      // recalculate determinant and gbigb_X
      ierr = xf_Error(xf_MapDet(G,G_X,dim,&g,gbigb_X));
      xf_AssertEqual(ierr, xf_OK);
      // finite difference check
      fd = (log(g)-log(g0))/eps;
      an = 0.5*(gbigb0_X[i] + gbigb_X[i]);
      if (fabs(fd-an) > tol){
	xf_printf("ping fail[dim=%d, i=%d]: fd = %.15E, an = %.15E, tol = %.5E\n", 
		  dim, i, fd, an, tol);
	return xf_Error(xf_PING_FAILED);
      } 
      // unperturb G
      for (d=0; d<dim*dim; d++) G[d] -= G_X[i*dim*dim+d]*eps;

    } // i
  } // dim

  return xf_OK;
}


TEST_xf_TestMeshMotionMap_Analytical()
{
  int ierr, i, j, dim, npoint;
  enum xfe_AnaMotionType MType;
  real X[3] = {0.4, 0.7, -0.65};
  real x0[3], G0[9], vg0[3], gbigb0_X[3];
  real  x[3],  G[9],  vg[3],  gbigb_X[3];
  real g, g0, Time;
  xf_MotionData MData;
  xf_AnaMotions *Motions = NULL;
  xf_AnaMotionsSet AnaMotionsSet;
  real fd, an;
  real eps = 1e-5, tol = 1e-8;
  
  for (MType=0; MType<xfe_AnaMotionLast; MType++){

  // create Motions structure
  ierr = xf_Error(xf_DefaultAnaMotions(&Motions, MType,
				       xfe_AnaBlend_Cubic));
  xf_AssertEqual(ierr, xf_OK);

  npoint = 1;
  dim = 3;
  Time = 1.1;
  MData.npoint = npoint;
  if (MType == xfe_AnaMotion_PerssonFS) dim = 2; //Persson mapping only meant for 2D
  MData.dim    = dim;

  // fill a motion set super-structure
  AnaMotionsSet.nMotion = 1;
  AnaMotionsSet.AnaMotions = Motions;

  // set MData variables
  MData.x = x0; MData.vg = vg0; MData.G = G0; MData.g = &g0; MData.gbigb_X = gbigb0_X;

  // call mesh motion routine
  ierr = xf_Error(xf_MeshMotionMap_Analytical(&AnaMotionsSet, npoint, dim, Time,
					      X, &MData));
  if (ierr != xf_OK) return ierr; 

  // ping X derivatives
  for (j=0; j<dim; j++){

    // perturb X
    X[j] += eps;

    // call motion routine
    MData.x = x; MData.vg = vg; MData.G = G; MData.g = &g; MData.gbigb_X = gbigb_X;
    ierr = xf_Error(xf_MeshMotionMap_Analytical(&AnaMotionsSet, npoint, dim, Time,
						X, &MData));
    if (ierr != xf_OK) return ierr;

    // check derivatives against finite difference
    // G
    for (i=0; i<dim; i++){
      fd = (x[i] - x0[i])/eps;
      an = 0.5*(G0[i*dim+j] + G[i*dim+j]);
      if (fabs(fd-an) > tol){
	xf_printf("G ping fail[i,j=%d,%d]: fd = %.15E, an = %.15E, tol = %.5E\n",
		  i,j, fd, an, tol);
	return xf_Error(xf_PING_FAILED);
      }
    } // j
    // gbigb_X
    fd = (log(g) - log(g0))/eps;
    an = 0.5*(gbigb0_X[j] + gbigb_X[j]);
    if (fabs(fd-an) > tol){
      xf_printf("gbigb_X ping fail[j=%d]: fd = %.15E, an = %.15E, tol = %.5E\n",
		j, fd, an, tol);
      return xf_Error(xf_PING_FAILED);
    }

    // unperturb X
    X[j] -= eps;

  } // j

  // ping Time derivative
  Time += eps;
  MData.x = x; MData.vg = vg; MData.G = G; MData.g = &g; MData.gbigb_X = gbigb_X;
  ierr = xf_Error(xf_MeshMotionMap_Analytical(&AnaMotionsSet, npoint, dim, Time,
					      X, &MData));
  if (ierr != xf_OK) return ierr; 

  for (j=0; j<dim; j++){
    fd = (x[j] - x0[j])/eps;
    an = 0.5*(vg0[j] + vg[j]);
    if (fabs(fd-an) > tol){
      xf_printf("time ping fail[j=%d]: fd = %.15E, an = %.15E, tol = %.5E\n", 
		  j, fd, an, tol);
      return xf_Error(xf_PING_FAILED);
    } 
    
    // unperturb Time
    Time -= eps;
  } // j

  // destroy AnaMotions
  ierr = xf_Error(xf_DestroyAnaMotions(Motions, xfe_True));
  xf_AssertEqual(ierr, xf_OK);

} //end MotionType loop

  return xf_OK;
}
