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

#ifndef _xf_MeshMotion_h
#define _xf_MeshMotion_h 1

/*
  FILE:  xf_MeshMotion.h

  This file contains the headers for functions dealing with the MeshMotion structure.

*/


/******************************************************************/
//   FUNCTION Prototype: xf_CreateMeshMotion
extern int 
xf_CreateMeshMotion( xf_MeshMotion **pMotion);
/*
PURPOSE:

  Creates a MeshMotion structure and all of its children.  Memory is
  allocated and initial (zero) values are set.

INPUTS:

  pMotion : address of pointer to MeshMotion structure

OUTPUTS: 

  None: structure is allocated

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyMeshMotion
extern int 
xf_DestroyMeshMotion( xf_MeshMotion *Motion);
/*
PURPOSE:

  Destroys a MeshMotion structure (including self).

INPUTS:

  Motion : pointer to MeshMotion structure

OUTPUTS: 

  None: structure is allocated

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyMeshMotion
extern int 
xf_DestroyMeshMotion( xf_MeshMotion *Motion);
/*
PURPOSE:

  Destroys a MeshMotion structure (including self).

INPUTS:

  Motion : pointer to MeshMotion structure

OUTPUTS: 

  None: structure is allocated

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_InitMotionData
extern void
xf_InitMotionData( xf_MotionData *MData);
/*
PURPOSE:

  Initializes motion data structure entries to 0 or NULL.

INPUTS:

  MData: pointer to structure containing motion data

OUTPUTS:

  None, entries of MData get initialized

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_CreateMotionData
extern int 
xf_CreateMotionData( xf_All *All, xf_MotionData **pMData);
/*
PURPOSE:

  Creates and initializes motion data structure.

INPUTS:

  pMData: address of pointer to structure containing motion data

OUTPUTS:

  Created, initialized (*pMData)

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyMotionData
extern void 
xf_DestroyMotionData( xf_MotionData *MData);
/*
PURPOSE:

  Destroys motion data structure, including self.

INPUTS:

  MData: pointer to structure containing motion data

OUTPUTS:

  None, MData is destroyed

RETURN: None

*/ 

/******************************************************************/
//   FUNCTION Prototype: xf_AllocMotionData
extern int 
xf_AllocMotionData( unsigned int AllocFlag, int npoint, int dim, xf_MotionData *MData);
/*
PURPOSE:

  Allocates MData according to AllocFlag.  See the bitwise allocation
  options defined in xf_MeshMotionStruct.h.

INPUTS:

  AllocFlag  : bits of this number indicate which vectors to allocate
  npoint     : number of points for which we need to allocate MData
  dim        : spatial dimension

OUTPUTS: 
 
  MData : pointer to structure containing allocated motion data

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_MeshMotionMap
extern int 
xf_MeshMotionMap( int egrp, int elem, xf_BasisData *PhiData,
		  xf_MeshMotion *Motion, int npoint, int dim, 
		  real Time, const real *X,  xf_MotionData *MData);
/*
PURPOSE:

  Calculates mapping between reference (X) and moving (x) coordinates,
  for multiple points with given coordinates and times.

INPUTS:

  egrp, elem : element on which map is being computed (only if using GCL)
  PhiData    : basis functions/gradients at the points (only if using GCL)
  Motions    : MeshMotion structure that defines the motion
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




#endif // end ifndef _xf_MeshMotion_h
