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

#ifndef _xf_AllStruct_h
#define _xf_AllStruct_h 1

/*
  FILE:  xf_AllStruct.h

  This file contains the top-level xflow data structure, which
  references the mesh, data, parameters, and equation set structures.

*/


#include "xf.h"
#include "xf_MeshStruct.h"
#include "xf_GeomStruct.h"
#include "xf_DataStruct.h"
#include "xf_ParamStruct.h"
#include "xf_EqnSetStruct.h"
#include "xfYu_ModelStruct.h"

typedef struct
{

  xf_Mesh *Mesh;
  /* Mesh structure and connectivity */

  xf_Geom *Geom;
  /* Geometry information */

  xf_DataSet *DataSet;
  /* Grid and solution data set.  Temporary data can be stored here. */

  xf_Param *Param;
  /* Case and solver parameters */

  xf_EqnSet *EqnSet;
  /* Equation-set specific information */

  Yu_Model *Model;
  /* Basically replace eqnset with Yu's specification*/

} xf_All; 


#endif // end ifndef _xf_AllStruct_h
