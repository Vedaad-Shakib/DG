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

#ifndef _xf_Output_h
#define _xf_Output_h 1

#include "xf_OutputStruct.h"

/*
  FILE:  xf_Output.h

  This file contains the headers for functions in xf_Output.c

*/

/******************************************************************/
//   FUNCTION Prototype: xf_DestroyCutPlaneIntersect
extern int
xf_DestroyCutPlaneIntersect(xf_CutPlaneIntersect *CutPlaneIntersect);
/*
PURPOSE:

  Destroys cut-plane intersection data

INPUTS:

  CutPlaneIntersect : data to be destroyed

OUTPUTS: 

  None

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_OutputExactErrorNorm
extern int
xf_OutputExactErrorNorm(xf_All *All, xf_Vector *U, xf_Vector *Ue, real *Value);
/*
PURPOSE:

  Calculates an L2 error norm between U and Ue on the domain in All.

INPUTS:

  All : All structure
  Name : output name
  U  : state vector
  Ue : exact state vector

OUTPUTS: 

  Value : computed norm

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_FindOutput
extern int 
xf_FindOutput( const xf_EqnSet *EqnSet, const char *Name, xf_Output **pOutput);
/*
PURPOSE:

  Finds an output entitled Name in the EqnSet output structure.

INPUTS:

  EqnSet : equation-set structure which to search
  Name : desired name of output

OUTPUTS: 

  (*pOutput) : pointer to output, if found.

RETURN:

  Error code
  xf_NOT_FOUND if output is not found
*/


/******************************************************************/
//   FUNCTION Prototype: xf_IsOutputUnsteady
extern int 
xf_IsOutputUnsteady( const xf_EqnSet *EqnSet, const char *Name, 
		     enum xfe_Bool *IsUnsteady);
/*
PURPOSE:

  Determines whether Output entitled Name is an unsteady output.

INPUTS:

  EqnSet : equation-set structure which to search
  Name : desired name of output

OUTPUTS: 

  (*IsUnsteady) : True if output is unsteady

RETURN:

  Error code
  xf_NOT_FOUND if output is not found
*/


/******************************************************************/
//   FUNCTION Prototype: xf_ReadStoredOutput
extern int 
xf_ReadStoredOutput( const xf_EqnSet *EqnSet, const char *Name, 
		     real *Value, real *ErrEst);
/*
PURPOSE:

  Reads values stored in Output->Value and Output->ErrEst (both
  optional) for an output stored in EqnSet and entitled Name.

INPUTS:

  EqnSet : equation-set structure which to search
  Name : desired name of output

OUTPUTS: 

  (*Value)  : value of output (optional, can be NULL)
  (*ErrEst) : error estimate of output (optional, can be NULL)

RETURN:

  Error code
  xf_NOT_FOUND if output is not found
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CalculateOutput
extern int 
xf_CalculateOutput( xf_All *All, const char *Name, const xf_Vector *U, 
		    real *Value, xf_Vector *Value_U, enum xfe_AddType AddFlag);
/*
PURPOSE:

  Calculates output that matches Name.  If not NULL, Value_U is
  set/incremented with the linearization of the output w.r.t the
  state, U.  AddFlag controls the setting vs. incrementing.

INPUTS:

  All : All structure
  Name : output name
  U  : state vector
  AddFlag : Controls whether Value_U is set or incremented

OUTPUTS: 

  Value : computed output
  Value_U : linearization of output w.r.t state

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_IncrementUnsteadyOutputs
extern int 
xf_IncrementUnsteadyOutputs( xf_All *All, const char *Name, 
			     enum xfe_TimeSchemeType TimeScheme, int nU,
			     xf_Vector **Ui, xf_Vector *Utemp, real Time, 
			     real TimeStep, enum xfe_Bool AtStart, 
			     enum xfe_Bool AtEnd, real *Value, xf_Vector **Value_U );
/*
PURPOSE:

  Increments the value of (or just evaluates in some cases) an
  unsteady output.  The output in question could be an integral over
  time, and this function is only called for one time node or slab at
  a time.  For this reason, the output is incremented on each call.

  The output is specified in Name, but if this input is NULL, all
  unsteady outputs in the key-value LogOutput will be considered.

INPUTS:

  All       : All structure
  Name      : output name; if NULL all outputs in LogOutput are used
  TimeScheme: time scheme used by calling function
  nU        : number of state vectors passed in (e.g. DG1 has 2)
  Ui        : vector of state vectors
  Utemp     : a temporary state vector (storage for interpolation)
  Time      : Start time of time step/slab.
  TimeStep  : size of time step
  AtStart   : True if at start of time integration (output will be initialized)
  AtEnd     : True if at end of time integration (can apply Final TimeNorm)

OUTPUTS: 

  Value     : computed output -- if this is NULL, the Output->Value will
              not be changed
  Value_U   : linearization of output w.r.t state vectors (nU of these)
              value is set (not incremented)
  
  Also, the Value field of each Output is incremented.  Note, initalization
  happens if AtStart=True.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_PingOutput
extern int
xf_PingOutput(xf_All *All, const char *Name, xf_Vector *U);
/*
PURPOSE:

  Pings output Name to determine if linearization is correct

INPUTS:

  All : All structure
  Name : output name
  U  : state vector

OUTPUTS: 

  None : output is pinged (info is written to stdout)

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_SetOutputDependentBCs
extern int
xf_SetOutputDependentBCs(xf_All *All, xf_Vector *U, xf_KeyValue *pOutputsForBCs, 
                         enum xfe_Bool PostSolve, enum xfe_Bool *Converged, 
                         enum xfe_Bool *Trim);
/*
PURPOSE:

  Sets boundary conditions that depend on outputs.  The outputs are
  taken from (*pOutputsForBCs) if this key-value structure has some
  entries in it.  If (*pOutputsForBCs) is an empty key-value
  structure, it is filled in using U to calculate the outputs, prior
  to calling the equation-set specific function that sets the BCs.

INPUTS:

  All : All structure
  U  : state vector
  (*pOutputsForBCs) : on input, may have outputs that should be used to set the BCs

OUTPUTS: 

  (*pOutputsForBCs) : if given as empty key-value, gets filled-in with outputs
                      calculated using the state U
  (*Converged) : True or False if BC is of trim type, always True otherwise.
  (*Trim) : True if BC is of trim type, False otherwise.

RETURN:

  Error Code
*/


/******************************************************************/
//   FUNCTION Prototype: xf_CompOutputSensitivity
extern int
xf_CompOutputSensitivities(xf_All *All, xf_Vector *U, xf_Output *Output);
/*
 PURPOSE:
 
 Computes the sensitivities associated to Output. The pointers to the 
 sensitivity structures are stored in the xf_Output structure.
 
 INPUTS:
 
 All : All structure
 U  : state vector
 Output : Output for which the sensitivities will be computed
 
 OUTPUTS: 
 
 None: Output get modified.
 
 RETURN:
 
 Error Code
 */


/******************************************************************/
//   FUNCTION Prototype: xf_ErrEstBC
extern int
xf_ErrEstBC(xf_All *All, xf_Vector *U, int nPsi, xf_Vector **Psi,
            xf_KeyValue *pOutputsForBCs);
/*
PURPOSE:

  Prints out output error estimates for all outputs for which the
  adjoints (Psi) are given.  The error in question is the error due to
  output-based BCs not being set exactly.  The function sets these BCs
  using the latest outputs from the state U, computes the resulting
  residual, and weights the residual using the adjoint in Psi.  On
  exit, the BCs are reset using given original outputs in
  (*pOutputsForBCs) ... subsequently the outputs in pOutputsForBCs are
  modified to be the latest ones based on the state.

INPUTS:

  All : All structure
  U  : state vector
  nPsi : number of adjoint vectors
  Psi : adjoint vectors
  (*pOutputsForBCs) : on input, contains old outputs -- BCs will be reset using these

OUTPUTS: 

  (*pOutputsForBCs) : on exit, contains outputs calculated using the state U

RETURN:

  Error Code
*/

extern int 
xf_GetLocalElem(xf_Mesh *Mesh, int egrpglob, int elemglob, 
                      int *egrploc, int *elemloc);

#endif // end ifndef _xf_Output_h
