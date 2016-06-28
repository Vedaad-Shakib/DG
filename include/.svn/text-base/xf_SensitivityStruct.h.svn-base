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

#ifndef _xf_SensitivityStruct_h
#define _xf_SensitivityStruct_h 1

/*
 FILE:  xf_SensitivityStruct.h
 
 This file contains the xflow structures for storing sensitivities
 
 */

/* Sensitivity types */
enum xfe_SensitivityType {
  xfe_SensitivityEqnset,  //sensitivity w.r.t. eqnset parameters
  xfe_SensitivityMesh,    //sensitivity w.r.t. mesh parameters
  xfe_SensitivitySolver,  //sensitivity w.r.t. solver parameters
  xfe_SensitivityNone,
  xfe_SensitivityLast
};

/* corresponding names */
static char *xfe_SensitivityName[xfe_SensitivityLast] = {
  "Eqnset",
  "Mesh",
  "Solver",
  "None"
};

/* Sensitivity structure */
typedef struct
{
  enum xfe_SensitivityType Type; //type of sensitivity (defined above)
  char *ParamName; //Type of parameter (angle of attack, gas constant, node position, etc.)
  int nVar; //rank of sensitivity vector, e.g., 1 for angle-of-attack, nNode for node positions..
  real *value; //value[i=0->nVar-1] stores derivative w.r.t. variable "i"
  xf_KeyValue KeyValue;  //this stores additional parameters for computing sentitivities (e.g. plane normal for AOA)
}
xf_Sensitivity;




#endif // end ifndef _xf_SensitivityStruct_h
