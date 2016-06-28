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

#ifndef _xf_Adapt_h
#define _xf_Adapt_h 1

/*
  FILE:  xf_Adapt.h

  This file contains the headers for the top-level adaptation functions

*/

#include "xf_AdaptStruct.h"


/******************************************************************/
//   FUNCTION Prototype: xf_GetNumRefOpt
extern int 
xf_GetNumRefOpt(enum xfe_ShapeType Shape, int *nref);

/*
PURPOSE:

  Returns number of refinement options for a given element Shape.
  This applies to hanging node refinement.

INPUTS:

  Shape : shape to consider for refinement

OUTPUTS: 

  (*nref) : number of refinement options for this Shape

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_AdaptAll
extern int 
xf_AdaptAll(xf_All *All, int iAdapt, enum xfe_Bool *pDoneAdapt,
            xf_Vector *Indicator);
/*
PURPOSE:

  Top-level adaptation function.  Determines what kind of error
  estimation and adaptation is requested and performs it.  Should be
  called with primal and adjoint (if necessary) solutions already
  present.

INPUTS:

  All : All structure
  iAdapt : adaptation iteration (for logging purposes)
  Indicator: externally provided indicator. If NULL, this 
             function will compute one.

OUTPUTS: 

  (*pDoneAdapt) : True if adaptation has terminated (e.g. tolerance met)

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_AdaptTimeHistData
extern int 
xf_AdaptTimeHistData(xf_TimeHistData *TimeHistData, int *RefIndicatorTime,
		     xf_TimeHistData **OldTimeHistData);
/*
PURPOSE:

  Adapts Time history data structure using RefIndicatorTime.  Assumes
  a finite element in time discretization.  Coarsening and refinement
  are supported.

INPUTS:

  TimeHistData : structure to adapt
  RefIndicatorTime : temporal refinement/coarsening indicator
                      0 = nothing requested
                      1 = refinement requested
                     -1 = coarsening requested

OUTPUTS: 

  TimeHistData : structures are reallocated
  OldTimeHistData : bare-bones original time history data, storing
                    time steps and slab sizes (optional output)
		   
RETURN: Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_SplitTimeHistData
extern int 
xf_SplitTimeHistData(xf_TimeHistData *TimeHistData, int *RefIndicatorTime, int iSplit);
/*
PURPOSE:

  Splits time slab iSplit or does a uniform time slab refinement if
  iSplit == -1.  Splits potentially many time slabs if
  RefIndicatorTime is provided.

INPUTS:

  TimeHistData : structure to adapt
  RefIndicatorTime : temporal refinement/coarsening indicator
                      0 = nothing requested
                      1 = refinement requested 
                     (optional -- iSplit used if NULL)
  iSplit : this is the index of the one time slab that will be
           refined.  -1 indicates uniform refinement.

OUTPUTS: 

  TimeHistData : structures are reallocated
		   
RETURN: Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_AdaptAllUnsteady
extern int 
xf_AdaptAllUnsteady(xf_All *All, int iAdapt, const char *SavePrefixNext, 
		    xf_TimeHistData *TimeHistData, enum xfe_Bool *pDoneAdapt);
/*
PURPOSE:

  Top-level unsteady adaptation function. Adapts in space (modifies
  mesh) and in time (modifies TimeHistData).

INPUTS:

  All : All structure
  iAdapt : adaptation iteration (for logging purposes)
  TimeHistData : current time history
  (*pDoneAdapt) : must be set on input; if True, this function will
                  write out a couple files (end of adapt) and exit.

OUTPUTS: 

  TimeHistData : new time history, with prescribed time steps
  (*pDoneAdapt) : True if adaptation has terminated (e.g. tolerance met)

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Definition: xf_AdaptRobust
extern int
xf_AdaptRobust(xf_All *All, xf_Vector *U, int *piAdaptRobust, 
               enum xfe_Bool *pAdapted, xf_SolverData *SolverData);

/* 
 PURPOSE:
 
 Performs mesh adaptation in order to improve convergence.
 
 INPUTS:
 
 All : All structure
 U : Primal state vector
 piAdaptRobust: pointer to index of the adaptation cycle
 pAdapted : if True, adaptation was performed.
 SolverData: structure storing the solver parameters
 
 OUTPUTS: 
 
 None. All->Mesh gets modified and data is transferred
 
 RETURN:
 
 Error Code 
 */

//Yu's add-on
extern int 
Yu_CreateOrFindAdaptIndicator(xf_All *All, enum xfe_Bool WithRefOpt, 
                              int numrank, xf_Vector **pAdaptIndicator);

extern int
Yu_MeshAdaptation(xf_All *All, Yu_Model *Model);

extern int
Yu_CreateRefIndicator(xf_All *All, xf_Vector *AdaptIndicator,
                      xf_Vector **pRefIndicator, int *nelemref);

#endif // end ifndef _xf_Adapt_h
