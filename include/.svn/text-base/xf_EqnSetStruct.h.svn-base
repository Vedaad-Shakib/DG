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

#ifndef _xf_EqnSetStruct_h
#define _xf_EqnSetStruct_h 1

/*
  FILE:  xf_EqnSetStruct.h

  This file contains the xflow equation-set structures

*/

#include "xf.h"
#include "xf_OutputStruct.h"


/* Residual term types */
enum xfe_ResTermType{
  xfe_ResTermConv,    /* Convection term  */
  xfe_ResTermDiff,    /* Diffusion term   */
  xfe_ResTermSource,  /* Source term      */
  xfe_ResTermUnknown, /* Unknown term     */
  xfe_ResTermLast
};

/* corresponding names */
static char *xfe_ResTermName[xfe_ResTermLast] = {
  "Convection",
  "Diffusion",
  "Source",
  "Unknown"
};




/* One Residual Term structure.*/
typedef struct
{
  enum xfe_ResTermType Type; 
  /* Type of residual term (conv, diff, source); the discretization is
     performed outside the equation set functions according to this
     type. */

  xf_KeyValue KeyValue;
  /* Keys and Values understood by the equation set functions for
     determining the details of this residual term (e.g. type of A
     matrix for a diffusion term, or type of flux function for a
     convection term) */

  enum xfe_Bool Active;
  /* If False, this term is not used during residual calculation.
     This flag is set to True on term creation.  Not read or
     written. */
}
xf_ResTerm;


/* Set of residual terms */
typedef struct
{
  int nResTerm;
  /* Number of residual terms */

  xf_ResTerm *ResTerm;
  /* Residual term structures */
}
xf_ResTerms;


typedef struct
{
  
  char *Type;
  /* String signifying type of initial condition; understood by
     EqnSet.  Read and written. */

  char *Function;
  /* String identifying any function associated with this IC.  NULL if
     no function.  Read and written. */

  char *Header;
  /* Header string for the IC Data.  Optional.  Read and written. */

  char *Data;
  /* Data string for this IC.  Optional.  Read and written. */

  char *AlterFunction;
  /* String identifying any alteration function associated with this
     IC.  NULL if no function.  Read and written. */

  char *AlterData;
  /* Data string for AlterFunction.  Optional.  Read and written. */

  enum xfe_Bool PriorSteadySolve;
  /* For unsteady cases, if true, the solver will be started with a
     steady-state initial condition. */

}
xf_IC;

/* Set of initial conditions */
typedef struct
{
  int nIC;
  /* Number of initial conditions.  Usually 1. */

  xf_IC *IC;
  /* Pointer to initial condition structure(s). */
}
xf_ICs;


typedef struct
{
  
  char *BFGTitle;
  /* Title of boundary face group to which this boundary condition
     applies. Read and written. */

  char *Type;
  /* String identifying the type of BC.  Understood by the equation
     set and hence equation-set specific.  Read and written */

  char *Function;
  /* String identifying any function associated with this BC.  NULL if
     no function.  Read and written. */

  char *Header;
  /* Header string for the BC Data.  Optional.  Read and written. */

  char *Data;
  /* Data string for this BC.  Optional.  Read and written. */

  char *OutputLinkage;
  /* Name of outputs on which this BC depends.  Read and written. */

  int nBCParam;
  /* Number of BC params.  Determined by equation set from parsing the
     Data list.  Not read or written. */

  real *BCParam;
  /* Vector of length nBCParam containing the BC params.  Set by
     equation set.  Not read or written. */

}
xf_BC;


/* Set of boundary conditions */
typedef struct
{
  int nBC;
  /* Number of boundary conditions.  Has to equal the number of
     boundary face groups. Read and written. */

  xf_BC *BC;
  /* Boundary condition structures.  Read and written. */
}
xf_BCs;


/* Set of outputs */
typedef struct
{
  int nOutput;
  /* Number of outputs. Read and written */

  xf_Output *Output;
  /* Output structures.  Read and written. */
}
xf_Outputs;



/* EqnSet structure definition */
typedef struct
{ 
  
  char *EqnSetLibrary;
  /* Name of the equation set library.  Read and written. */

  int Dim;
  /* Spatial dimension for this equation set, if applicable. Read and
     written. */

  xf_KeyValue KeyValue;
  /* Key-Value variables.  All equation-set-specific parameters
     are stored here; they are read and written. */

  xf_ResTerms *ResTerms;
  /* Spatial residual terms.  Read and written. */

  xf_ICs *ICs;
  /* Initial conditions.  Read and written. */

  xf_BCs *BCs;
  /* Boundary conditions.  Read and written. */

  xf_Outputs *Outputs;
  /* Outputs.  Read and written. */

  int StateRank;
  /* Rank of the state vector.  Read and written. */

  char **StateName;
  /* StateRank strings containing the state names.  Read and written. */

  int nPosInState;
  /* number of possible equation-set states (understood by the
     equation-set library).  Not read or written. */

  int *PosInState;
  /* An index used by the equation set to convert between an
     enumerated equation-set specific type and the position in the
     state vector.  For example, an equation-set function that wants
     to access xfe_Density would pull off the value corresponding
     to PosInState[xfe_Density].  This value is not read or written,
     as it can be constructed from StateName */


  int nIParam;
  /* Number of integer parameters passed to equation-set specific
     functions.  Not read or written; set during equation-set
     registration.  */

  char **IParamKey;
  /* Names of integer parameter keys.  Not read or written. */

  int nRParam;
  /* Number of real parameters passed to equation-set specific
     functions.  Not read or written. */

  char **RParamKey;
  /* Names of real parameter keys.  Not read or written.  */

  int nAuxU;
  /* Number of auxiliary vectors (interpolated) passed to equation-set
     specific functions.  Not read or written. */

  char **AuxUNames;
  /* Names of auxiliary vectors.  Not read or written.  */
  
  int *PosInAuxU;
  /* An index used by the equation set to convert between an
     enumerated equation-set specific type and the position in the
     auxiliary vector list.  For example, if an eqnset function wants
     to access xfe_AuxVectorVelocity, it would refer to the
     PosInAux[xfe_AuxVectorVelocity] component of the AuxU vector
     array.   Not read or written. */
  
}
xf_EqnSet;


#endif // end ifndef _xf_EqnSetStruct_h
