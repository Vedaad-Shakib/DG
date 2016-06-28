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

#ifndef _xf_ErrEst_h
#define _xf_ErrEst_h 1

/*
  FILE:  xf_ErrEst.h

  This file contains the headers for error estimation functions

*/

//#include "xf_AdaptStruct.h"


/******************************************************************/
//   FUNCTION Prototype: xf_PreFineSpaceSolve
extern int 
xf_PreFineSpaceSolve(xf_All *All, xf_KeyValue *pKeyValueOrig);
/*

PURPOSE: 

  Prepares parameters for a fine space solve.  Specifically, the
  values from all key-value pairs whose keys begin with "FineSpace_"
  overwrite the values of the original key-values.
  e.g. "FineSpace_nIterNonlinear" overwrites "nIterNonlinear".
  Original keys and values are stored in KeyValueOrig, which is
  initialized in this function.
  
INPUTS:

  All : All file

OUTPUTS: 

  (*pKeyValueOrig) : Original key value pairs (before they were overwritten)


RETURNS: Error Code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_PostFineSpaceSolve
extern int 
xf_PostFineSpaceSolve(xf_All *All, xf_KeyValue *pKeyValueOrig);
/*

PURPOSE: 

  This function should be called after a fine space solve.  It reverts
  all key values from KeyValueOrig to All.
  
INPUTS:

  All : All file
  KeyValueOrig : Original key value pairs (before they were overwritten)
                 Note, KeyValueOrig is destroyed in this function

OUTPUTS: None, All has its key-value pairs reverted

RETURNS: Error Code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ErrEstOutput
extern int 
xf_ErrEstOutput(xf_All *All, char *OutputName, char *VariableSet, 
		enum xfe_Bool StoreFineSpace, enum xfe_Bool ReuseFineSpace,
		xf_Vector *ElemIndicator, enum xfe_Bool ElemIndSign, 
                real *pErrOutput);
/*
PURPOSE:

  Calculates an output-based error estimate, and stores it in
  (*pOutputError).  If requested, an elemental error indicator is also
  computed and stored in ErrIndicator.  Error estimates under various
  refinement options of each element can also be returned, if
  WithRefOpt is true.

INPUTS:

  All : All structure
  OutputName : name of output for which to calculate error estimate
  VariableSet : name of entropy variable set for which to calculate error estimate
                (only one of OutputName, VariableSet can be non-null in
		 the string sense, but both should be passed in)
  StoreFineSpace : if True, fine-space primal and dual solutions will be
                   stored in All->DataSet (assuming ReuseFineSpace is False)
  ReuseFineSpace : if True, fine-space primal and dual solutions will not be
                   created; rather they will be searched for in All->DataSet.
		   Note: to delete the fine-space solutions, call with
		   ReuseFineSpace=True, but StoreFineSpace=False.
  ElemIndSign: if True, ElemIndicator is a signed value.

OUTPUTS: 

  ElemIndicator : Calculated error estimate.  This vector must exist and
                  its size must be consistent with the WithRefOpt and
		  IsotropicOnly flags.

  (*pOutputError) : unadulturated error estimate (optional)

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Definition: xf_SeparateErrEst

extern int 
xf_SeparateErrEst(xf_All *All, const char *CombinedOutputName, 
                  xf_Vector *Eh, xf_Vector *UH_h);
/*
 PURPOSE:
 
 Calculates an error estimate for each individual output bundled 
 in the CombinedOutputName through a linearization. 
 
 INPUTS:
 
 All : All structure
 CombinedOutputName : string for locating the correct output
 Eh : Solution error estimate mapped into a finer space
 UH_h : Current solution mapped into the same finer space as Eh
 
 OUTPUTS: 
 
 ErrEst for each output in CombinedOutput is calculated and 
 SumOutputWeights are calculated.
 
 RETURN:
 
 Error Code
 
 */

/******************************************************************/
//   FUNCTION Definition: xf_ErrEstSolution

extern int 
xf_ErrEstSolution(xf_All *All, xf_Vector *U, const char *CombinedOutputName, 
                  enum xfe_Bool StoreFineSpace, enum xfe_Bool ReuseFineSpace);
/*
 PURPOSE:
 
 Calculates a solution error estimate vector, and stores it in
 Eh. For framework reference see R. Hartmann SIAM 2008, 
 Vol. 31, No. 1
 
 INPUTS:
 
 All : All structure
 U: Primal state
 CombineOutputName: string with combined output name for SeparateErrEst
 StoreFineSpace : if True, fine-space state approximate solution will be
 stored in All->DataSet (assuming ReuseFineSpace is False)
 ReuseFineSpace : if True, fine-space solution will not be
    created; rather they will be searched for in All->DataSet.
 
 OUTPUTS: 
 
 None : Eh = Uh - UH->h is computed and sent to SeparateErrEst
 
 RETURN:
 
 Error Code
 
 */


#endif // end ifndef _xf_ErrEst_h
