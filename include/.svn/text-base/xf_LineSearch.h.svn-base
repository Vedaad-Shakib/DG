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

#ifndef _xf_LineSearch_h
#define _xf_LineSearch_h 1

/*
 FILE:  xf_LineSearch.h
 
 This file contains the headers for the functions related 
 to the line-search
 
 */

#define xf_mu1 1e-4
#define xf_mu2 9e-1

/******************************************************************/
//   FUNCTION Definition: xf_LineSearch
extern int 
xf_LineSearch(xf_All *All, xf_Vector *U, xf_Vector *P, 
              xf_SolverData *SolverData);

/*
 PURPOSE:
 
 Calculates how much of P should me added to U
 
 INPUTS:
 
 All: All structure
 U: primal state 
 P: search direction
 SolverData: structure that stores the step length
             (SolverData->UpdateFrac)
 
 OUTPUTS: 
 
 None: SolverData get modified
 
 RETURN:
 
 Error Code
 */


#endif // end ifndef _xf_LineSearch_h