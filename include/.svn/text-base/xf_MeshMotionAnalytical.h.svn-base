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

#ifndef _xf_MeshMotionAnalytical_h
#define _xf_MeshMotionAnalytical_h 1

/*
  FILE:  xf_MeshMotionAnalytical.h

  This file contains the headers for functions dealing with analytical
  mesh motions.

*/

/******************************************************************/
//   FUNCTION Prototype: xf_CreateAnaMotionsSet
extern int 
xf_CreateAnaMotionsSet( xf_AnaMotionsSet **pAnaMotionsSet);
/*
PURPOSE:

  Creates analytical mesh motion structure set

INPUTS:

  pAnaMotions : pointer to analytical mesh motion structure set

OUTPUTS: 
 
  (*pAnaMotionsSet) : allocated/initialized analytical mesh motion 
                      structure set

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_AllocAnaMotions
extern int 
xf_AllocAnaMotions( xf_AnaMotions *AnaMotions, int nTerm);
/*
PURPOSE:

  Allocates analytical mesh motions structure with nTerm terms

INPUTS:

  AnaMotions : analytical mesh motions structure
  nTerm      : number of analytical terms to allocate

OUTPUTS: 

  AnaMotions : modified 

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_AllocAnaMotionsSet
extern int 
xf_AllocAnaMotionsSet( xf_AnaMotionsSet *AnaMotionsSet, int nMotion);
/*
PURPOSE:

  Allocates set of nMotion analytical mesh motion structures

INPUTS:

  AnaMotionsSet : set of analytical mesh motion structures
  nTerm         : number of mesh motions to allocate

OUTPUTS: 

  AnaMotions : modified 

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_DestroyAnaMotionsSet
extern int 
xf_DestroyAnaMotionsSet( xf_AnaMotionsSet *AnaMotionsSet);
/*
PURPOSE:

  Destroys set of analytical mesh motion structures

INPUTS:

  AnaMotionsSet : set of analytical mesh motion structures

OUTPUTS: 

  None, structure is destroyed

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_MeshMotionMap_Analytical
extern int 
xf_MeshMotionMap_Analytical( xf_AnaMotionsSet *AnaMotionsSet, int npoint, int dim, 
			     real Time, const real *X,  xf_MotionData *MData);
/*
PURPOSE:

  Analytic-motion version of xf_MeshMotionMap.  Calculates mapping
  between reference (X) and moving (x) coordinates. 

INPUTS:

  AnaMotions : Analytical motions set structure that defines the motion
  npoint     : number of points to map; each variable is stored 
               unrolled sequentially
  dim        : spatial dimension
  Time       : physical Time
  X          : reference coordinates, size [dim] for each point

OUTPUTS: 
 
  MData : pointer to structure containing motion data

RETURN:

  Error Code
*/



#endif // end ifndef _xf_MeshMotionAnalytical_h
