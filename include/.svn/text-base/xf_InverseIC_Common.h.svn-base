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

#ifndef _xf_InverseIC_Common_h
#define _xf_InverseIC_Common_h 1

/*
  FILE:  xf_InverseIC_Common.h

  This file contains function prototypes for functions in xf_InverseIC_Common.c

*/



/******************************************************************/
//   FUNCTION Prototype: xf_ReadICParamFile
extern int
xf_ReadICParamFile(const char *ICParamFile, int *pnParam, 
		     char ***pParamName, real **pParamRange, real **pParamValue);
/*
PURPOSE:

  Reads a parameter file ICParamFile and fills in the required
  vectors (see below).  Format of the param file should be:
 
  ParamName  Min  Max  Start

INPUTS:

  ICParamFile : name of file to read
  
OUTPUTS: 
  
  (*pnParam) : number of parameters
  (*pParamName) : names of parameters
  (*pParamRange) : min and max ranges for parameters:
                   (*pParamRange)[2*iParam+0] = minimum
                   (*pParamRange)[2*iParam+1] = maximum
  (*pParamValue) : starting guess values for parameters

RETURN: Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ReadICParamFileEnsemble
extern int
xf_ReadICParamFileEnsemble(const char *ICParamFileEnsemble, int *pnParam,
        char ***pParamName, real **pParamRange, real **pParamWindowRange);
/*
PURPOSE:

  Reads a parameter file ICParamFileEnsemble and fills in the required
  vectors (see below).  Format of the param file should be:
 
  ParamName  Min  Max InitMin InitMax

INPUTS:

  ICParamFileEnsemble : name of file to read
  
OUTPUTS: 
  
  (*pnParam) : number of parameters
  (*pParamName) : names of parameters
  (*pParamRange) : min and max ranges for parameters:
                   (*pParamRange)[2*iParam+0] = minimum
                   (*pParamRange)[2*iParam+1] = maximum
  (*pParamWindowRange) : min and max intiial ranges for parameters:
                   (*pParamWindowRange)[2*iParam+0] = minimum
                   (*pParamWindowRange)[2*iParam+1] = maximum


RETURN: Error code

*/

/******************************************************************/
//   FUNCTION Prototype: xf_ICParam2State
extern int
xf_ICParam2State(xf_All *All, int nParam, char **ParamName, real *ParamValue,
		 real *epsParam, xf_Vector *U, xf_Vector *U_mu);
/*
PURPOSE:

  Converts a set of initial condition parameters to the actual initial
  condition vector.  The gradient with respect to the parameters is
  also computed if requested.
 
INPUTS:

  All : all file
  nParam : number of parameters to set
  ParamName : names of parameters
  ParamValue : values of parameters
  epsParam : finite difference epsilons for U_mu calculation
  
OUTPUTS: 

  U : initial condition vector
  U_mu : gradient of initial vector w.r.t parameters (optional)

RETURN: Error code

*/


/******************************************************************/
//   FUNCTION Prototype : xf_ReadICOutputFile
extern int
xf_ReadICOutputFile(const char *OutputFile, int *pnOut, char ***pOutName,
		    int **pOutTimeIndex, real **pOutValue, real **pOutStdErr);
/*
PURPOSE:

  Reads a sensor output file, with format
 
  OutputName TimeIndex Value StdError

INPUTS:

  OutputFile : name of file to read
  
OUTPUTS: 
  
  (*pnOut) : number of outputs
  (*pOutName) : names of outputs
  (*pOutTimeIndex) : time indices indicating when outputs were taken
                     note, time step is defined in .job file
  (*pOutValue) : output readings
  (*pOutStdErr) : output measurement standard errors
  

RETURN: Error code

*/



/******************************************************************/
//   FUNCTION Prototype: xf_LoadSensitivityVectors
extern int
xf_LoadSensitivityVectors(xf_All *All, int nOut, char **OutName, 
			  int *OutTimeIndex, xf_Vector ***pJ_U0);
/*
PURPOSE:

  Check for or create output sensitivity vectors
  Also, Load these vectors into memory           

INPUTS:

  All : All structure
  nOut : number of outputs
  OutName : names of outputs
  OutTimeIndex : time indices of outputs
  
OUTPUTS: 
  
  J_U0 : sensitivity vectors

RETURN: Error code

*/


/******************************************************************/
//   FUNCTION Prototype: xf_CalculateLinearOutputs
extern int
xf_CalculateLinearOutputs(xf_All *All, int nParam, char **ParamName, real *ParamValue,
			  int nOut, xf_Vector **J_U0, real *epsParam, real *J, 
			  real *J_mu, xf_Vector **pU0);
/*
PURPOSE:

  Calculates outputs J for given initial condition parameters using
  adjoint-based sensitivity vectors J_U0.  Assumes linear U0 to J
  relationship:

  J(i) = J_U0(i)^T * U0

INPUTS:

  All : All structure
  nParam : number of parameters to set
  nOut : number of outputs
  J_U0 : nOut sensitivity vectors
  pU0  : pointer to a state vector pointer (NULL on first call)
         used internally in calculation; pass in on all future calls
  
OUTPUTS: 

  J : nOut outputs  
  (*pU0) : modified state vector 


RETURN: Error code

*/




#endif // end ifndef _xf_InverseIC_Common_h

