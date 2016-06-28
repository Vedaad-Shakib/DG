/*------------------------------------------------------------------*/
/* XFLOW: A discontinuous Galerkin finite element software library. */
/*                                                                  */
/*                    Copyright  2007-2008                          */
/*           Krzysztof J. Fidkowski, kfid@alum.mit.edu              */
/*                                                                  */
/*                    Copyright  2008-2011                          */
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

#ifndef _xfYu_Solver_h
#define _xfYu_Solver_h 1

/*
  FILE:  xf_Solver.h

  This file contains the headers for the top-level solver functions

*/

#include "xf_SolverStruct.h"
#include "xf_SolverTools.h"
#include "xf_AdaptStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_ApplyTimeScheme
extern int
xfYu_ApplyTimeScheme(xf_All *All, Yu_Model *Model, Yu_Limiter ** Limiter, 
                    const char *SavePrefix, enum xfe_Bool RestartFlag,
                    xf_Vector *U0, xf_TimeHistData *TimeHistData);
/*
PURPOSE:

  Runs the unsteady time-stepping solver.  Pulls info from All.

  Each implicit time scheme is written as:
   
  M/dt * (c0*u^{n+1} + c1*u^n + c2*u^{n-1} + ... ) + R(u^{n+1}) = 0

  The coefficients are particular to the unsteady scheme.


INPUTS:

  All: All structure
  SavePrefix : prefix for writing unsteady write-interval files.
  RestartFlag: True if the run is restarted, in which case the appropriate 
               time-index states should be present in All->DataSet.  For non-
	       restarted runs, additional vectors are created and initialized
	       to the time-index=0 state.
  U0: initial/time-index=0 state

  TimeHistData : structure used to accumulate history of the time
                 values run and any outputs at all times.  Optional:
                 can pass in as NULL.  If not NULL, the structure must
                 have been created before the call, and the desired
                 number of Outputs must be stored in nOutput, with
                 corresponding names in OutputNames.  The Time vector
                 and the OutputValues vector are reallocated in this
                 function.

OUTPUTS:

  U0: final time state data is stored here
  TimeHistData : If given as not NULL, modified according to the above
                 description.

RETURN:

  Error code
*/







#endif // end ifndef _xf_Solver_h
