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

#ifndef _xf_SolverTools_h
#define _xf_SolverTools_h 1

#include"xfYu_Model.h"

/*
  FILE:  xf_SolverTools.h

  This file contains the headers for the top-level solver functions

*/


/******************************************************************/
// counterpart of the fcn: xf_InitAlterState
extern int 
xfYu_InitAlterState(xf_All *All, Yu_Model *Model, char **FcnName, 
                    real *FcnParam, xf_Vector *U, const int sr); 
/******************************************************************/
//   FUNCTION Prototype:  xf_CheckUserHalt
extern enum xfe_Bool
xf_CheckUserHalt(char *fname);
/*
PURPOSE:

  Returns True if file named fname (or "STOP" if NULL is passed in)
  exists with the first line consisting of the word "STOP".

INPUTS:

  fname : file name to check for first line containing "STOP".  If
          this is passed in as NULL, the filename "STOP" is used.

OUTPUTS: None

RETURN: True if STOP file exists.

*/

/******************************************************************/
//   FUNCTION Prototype: xf_MultMassMatrix
extern int
xf_MultMassMatrix(xf_All *All, real c, xf_Vector *U);
/*
PURPOSE:
   
  Multiplies U by the Mass matrix, M:

                U = c*M*U
 
  where c is a real constant.  The mass matrix is found from/created
  in All.

INPUTS:

  All: All structure
  c : constant in front of M*U product
  U : state vector

OUTPUTS:

  U : modified state vector

RETURN:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_MultInvMassMatrix
extern int
xf_MultInvMassMatrix(xf_All *All, real c0, xf_Vector *dt, xf_Vector *U);
/*
PURPOSE:
   
  Multiplies U by the inverse Mass matrix, M^{-1}:

                U = c*M^{-1}*U
 
  where c is a real constant = c0 or c0*dt(elem) if dt is not NULL.
  The inverse mass matrix is found from/created in All.

INPUTS:

  All: All structure
  c0 : constant in front of M^{-1}*U product
  dt : vector of local time steps (optional)
  U : state vector

OUTPUTS:

  U : modified state vector

RETURN:

  Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_MultRootMassMatrix
extern int
xf_MultRootMassMatrix(xf_All *All, real c, enum xfe_Bool InverseFlag,
		      xf_VectorSet *VS);
/*
PURPOSE:
   
  Multiplies all vectors in VS by the root of the (inverse) Mass matrix:
  
         U = c*M^{.5}*U  or U = c*M^{-.5}*U

  where c is a real constant.  The (inverse) mass matrix is found
  from/created in All.

INPUTS:

  All: All structure
  c : constant in front of product
  InverseFlag : if true, M^{-1} is used instead of M
  VS : vector set; all vectors in set are affected

OUTPUTS:

  VS: modified vector set

RETURN:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_MultStiffMatrix
extern int
xf_MultStiffMatrix(xf_All *All, real c, xf_Vector *U);
/*
PURPOSE:
   
  Multiplies U by the Stiffness matrix, K:

                U = c*K*U    K_ij = int(gPhi{i}*gphi{j})
 
  where c is a real constant.  The stiffness matrix is constructed on
  the fly -- hence this function is not optimized.

INPUTS:

  All: All structure
  c : constant in front of K*U product
  U : state vector

OUTPUTS:

  U : modified state vector

RETURN:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_AddMassMatrix
extern int
xf_AddMassMatrix(xf_All *All, real c, xf_Vector *dt, xf_Vector *U, 
		 xf_Vector *R, xf_JacobianMatrix *R_U, xf_SolverData *SolverData);
/*
PURPOSE:
   
  Adds Mass matrix to residual and Jacobian according to:

    R += c*M*U, R_U += (c+1/dt)*M

  where c is a real constant.  The mass matrix is found from/created
  in All.  R, R_U, dt can be NULL;

INPUTS:

  All: All structure
  c : constant in front of M
  dt : time step vector (1 time step for each element), optional
  U : state vector

OUTPUTS:

  R : modified residual vector, optional
  R_U : modified Jacobian matrix, optional

RETURN:

  Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_InitAlterState
extern int 
xf_InitAlterState(xf_All *All, const char *FcnName, 
		  const char *FcnParam, xf_Vector *U);
/*
PURPOSE:

  General function for initializing or altering a state vector.  This
  is the core function for which xf_InitState is a wrapper.

INPUTS:

  All : All structure
  FcnName : Name of alteration function (or NULL)
  FcnParam : string of parameters for alteration function (or NULL)
  U  : state; input only for alteration

OUTPUTS: 

  U : modified or initialized state

RETURN: Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_InitState
extern int 
xf_InitState(xf_All *All, xf_Vector *U);

/*
PURPOSE:

  Initializes a state vector U according to the initial conditions
  specified in the EqnSet structure.  Performs least-squares
  projection if a function of space is specified for the initial
  condition.

INPUTS:

  All : All structure containing the EqnSet with ICs

OUTPUTS: 

  U : vector of data that gets initialized.  Must be pre-allocated.

RETURN:

  Error Code (e.g. if the initial condition type or data are not understood)
*/


/******************************************************************/
//   FUNCTION Prototype: xf_AlterState
extern int 
xf_AlterState(xf_All *All, const char *FcnName, const char *FcnParam,
	      xf_Vector *U);

/*
PURPOSE:

  Similar to InitState, but the state is altered instead of initialized.  

INPUTS:

  All : All structure containing the EqnSet with ICs
  FcnName : name of altering function
  FcnParam: parameters for function, rolled into a string
  U  : initial state for function

OUTPUTS: 

  U : vector of data that gets altered.  Must be pre-allocated.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ChangeVariableSet
extern int 
xf_ChangeVariableSet(xf_All *All, const char *VariableSet, xf_Vector *U);
/*
PURPOSE:

  Changes entries in state vector U from one variable set to another.
  Uses QR projection to determine the coefficients of projection
  functions of the same Basis,Order as currently in U.  Calls
  equation-specific function at each quadrature point.

INPUTS:

  All : All structure
  VariableSet : desired variable set (understood by the equation set)
  U : conservative state variable vector
  
OUTPUTS: 

  U : vector of transformed variables in variableSet. Overwrites conservative
      variables in U

RETURN:

  Error Code 
*/

/******************************************************************/
//   FUNCTION Prototype: xf_ScaleState
extern int 
xf_ScaleState(xf_All *All, xf_ICs *ICsOrig, xf_Vector *U);
/*
PURPOSE:

  Scales state vector U in some EqnSet specific way based on the
  original initial condition, ICsOrig, and the new IC in All->EqnSet.
  Usefule for restarting a calculation during parameter sequencing
  (e.g. ramping up the Mach number).

INPUTS:

  All : All structure
  ICsOrig : original IC that gave rise to the state U
  U : conservative state variable vector
  
OUTPUTS: 

  U : scaled state vector based on All->EqnSet->ICs and ICsOrig.

RETURN:

  Error Code 
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReconstructVector
extern int 
xf_ReconstructVector(xf_All *All, int OrderIncrement, xf_Vector *U,
		     enum xfe_Bool AllocFlag, xf_Vector **pV);
/*
PURPOSE:

  Spatially reconstructs a vector to higher order using nearest
  neighbor information.  The vector must be spatially interpolated.
  The reconstruction is performed via least-squares, with quadrature
  integration over the element volumes.  The reconstructed vector is
  stored in (*pV).

INPUTS:

  All            : All structure
  OrderIncrement : Increment of current order during reconstruction
  U              : vector to reconstruct
  AllocFlag      : if True, space for (*pV) will be allocated in the call
  
OUTPUTS: 

  (*pV) : reconstructed vector

RETURN:

  Error Code 
*/


/******************************************************************/
//   FUNCTION Prototype: xf_VectorTextOut
extern int
xf_VectorTextOut(const char *fname, enum xfe_Bool Project2Lagrange, xf_Vector *U);
/*
PURPOSE:

  Writes out a vector to a text file.

INPUTS:

  fname : file name
  Project2Lagrange: if True, vector will be projected to Lagrange before writing
  U : state vector

OUTPUTS: 

  None, files are written

RETURN: Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_DumpSystemSparse
extern int
xf_DumpSystemSparse(xf_All *All, xf_Vector *U);
/*
PURPOSE:

  Serial only.  Calculates Jacobian matrix and Mass matrix, and writes
  to A.txt and M.txt respectively, using sparse storage format (row,
  col, value).

INPUTS:

  All : all structure
  U : vector about which to calculate Jacobian R_U

OUTPUTS: 

  None, files are written

RETURN: Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_GetTimeSchemeInfo
extern int
xf_GetTimeSchemeInfo(enum xfe_TimeSchemeType TimeScheme, 
		     enum xfe_Bool *Explicit, enum xfe_Bool *MultiStep,
		     enum xfe_Bool *MultiStage);
/*
PURPOSE:

  Returns information about the requested TimeScheme

INPUTS:

  TimeScheme : the time scheme in question

OUTPUTS: 

  (*Explicit)   : True if TimeScheme is explicit   (optional)
  (*MultiStep)  : True if TimeScheme is MultiStep  (optional) 
  (*MultiStage) : True if TimeScheme is MultiStage (optional) 

RETURN: Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_CreateTimeHistData
extern int 
xf_CreateTimeHistData( xf_TimeHistData **pTimeHistData);
/*
PURPOSE:

  Creates (allocates) a Time history data structure.

INPUTS:

  pTimeHistData : pointer to structure to allocate

OUTPUTS: 

  (*pTimeHistData) : allocated structure

RETURN: Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_DestroyTimeHistData
extern int 
xf_DestroyTimeHistData( xf_TimeHistData *TimeHistData);
/*
PURPOSE:

  Destroys a Time history data structure, including all sub-arrays.

INPUTS:

  TimeHistData : pointer to structure to destroy

OUTPUTS: 

  None; TimeHistData is destroyed

RETURN: Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ReadTimeHistData
extern int 
xf_ReadTimeHistData( const char *fname, const char *LogOutput,
		     xf_TimeHistData **pTimeHistData);
/*

PURPOSE: 

  Reads an ASCII TimeHistData file into a structure.
  
INPUTS:

  fname : name of TimeHistData file to read
  LogOutput : single string containing space-separated names of
              outputs that will be logged.

OUTPUTS: 

  (*pTimeHistData) : allocated structure

RETURNS: Error Code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_CreateUniformTimeHistData
extern int 
xf_CreateUniformTimeHistData( xf_All *All, const char *LogOutput,
			      xf_TimeHistData **pTimeHistData);
/*
PURPOSE:

  Creates (allocates) a Time history data structure and fills it with
  uniform time step information.

INPUTS:

  All : all structure containing time scheme and other parameters
  LogOutput : single string containing space-separated names of
              outputs that will be logged.

OUTPUTS: 

  (*pTimeHistData) : allocated structure

RETURN: Error code

*/





#endif // end ifndef _xf_SolverTools_h
