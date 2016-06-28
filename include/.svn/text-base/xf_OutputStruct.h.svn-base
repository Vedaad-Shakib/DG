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

#ifndef _xf_OutputStruct_h
#define _xf_OutputStruct_h 1

/*
  FILE:  xf_OutputStruct.h

  This file contains the xflow Output-specific structures

*/

#include "xf_IO.h"
#include "xf_SolverStruct.h"
#include "xf_SensitivityStruct.h"

/* Output type */
enum xfe_OutputType {
  xfe_DomainIntegral, 
  xfe_BoundaryIntegral, 
  xfe_LineIntegral,
  xfe_CutPlaneIntegral,
  xfe_PointValue,
  xfe_SumOutput,
  xfe_OutputLast
};

/* corresponding names */
static char *xfe_OutputName[xfe_OutputLast] = {
  "DomainIntegral", 
  "BoundaryIntegral",
  "LineIntegral",
  "CutPlaneIntegral",
  "PointValue",
  "SumOutput"
};

/* Domain Integral norm type */
enum xfe_DomainNormType {
  xfe_DomainNormNone,       // No norm, just a straight-up integral
  xfe_DomainNormL1,         // L1 norm: absolute value inside integral
  xfe_DomainNormL2,         // L2 norm: ^2 inside integral, ^.5 outside
  xfe_DomainNormL2Error,    // L2 Error norm: needs an associated function
  xfe_DomainNormL2Discrete, // Discrete L2 norm of element averages
  xfe_DomainNormL2PerVol,   // L2 norm per unit volume
  xfe_DomainNormLast
};

/* corresponding names */
static char *xfe_DomainNormName[xfe_DomainNormLast] = {
  "None", 
  "L1",
  "L2",
  "L2Error",
  "L2Discrete",
  "L2PerVol",
};


/* Time norm type (for unsteady runs)*/
enum xfe_TimeNormType {
  xfe_TimeNormNone,           // no norm; output calculated at every time
  xfe_TimeNormFinal,          // output measured at final time
  xfe_TimeNormPoint,          // output measured at specified time
  xfe_TimeNormIntegral,       // time integral of output
  xfe_TimeNormSquareIntegral, // time integral of square of output
  xfe_TimeNormLast
};

/* corresponding names */
static char *xfe_TimeNormName[xfe_TimeNormLast] = {
  "None",
  "Final", 
  "Point",
  "Integral",
  "SquareIntegral"
};


/* CutPlaneIntersect data structure: stores result of intersecting a
   cut plane with mesh elements.  Also used for storing intersections
   of elements with lines. */
typedef struct
{
  int  nelem;   // number of elements intersected
  int  *egrp;   // group numbers of intersected elements [nelem]
  int  *elem;   // element numbers of intersected elements [nelem]
  int  *nquad;  // number of quadrature points per element [nelem]
  real **xquad; // ref-space coords of quad points [nelem x (nquad[elem]*dim)]
  real **wquad; // weights of quad points [nelem x nquad[elem]]
}
xf_CutPlaneIntersect;


/* Output structure */
typedef struct
{
  char *Name;
  /* Unique identifying name of output. Read and written. */

  enum xfe_OutputType Type;
  /* DomainIntegral, BoundaryIntegral, LineIntegral, etc. Read and
     written. */

  enum xfe_DomainNormType DomainNorm;
  /* Relevant only for DomainIntegral output.  Norm used for the
     output. */
 
  enum xfe_TimeNormType TimeNorm;
  /* Relevant only for unsteady outputs.  Describes how the output is
     measured in the time domain (point value, integral, etc.). */

  enum xfe_Bool UsesFlux;
  /* Relevant only for BoundaryIntegral output.  If true, output is
     based on boundary flux. Read and written. */

  char *ScalarName;
  /* Relevant if UsesFlux == False.  Stores the name of the scalar to
     be computed for the output.  Understood by the eqnset.  Read and
     written. */

  char *VectorName;
  /* Relevant if UsesFlux == False, and only valid for
     BoundaryIntegral outputs.  Stores the name of the vector to be
     dotted with the boundary normal and integrated to computed the
     output.  Understood by the eqnset.  Read and written. */

  char *Function;
  /* String identifying any function associated with this Output.
     NULL if no function.  Read and written. */

  char *Data;
  /* Data string, for example to be used with Function.  Read and
     written. */

  int nFluxComponent;
  /* Relevant if UsesFlux == True.  Number of flux components that
     will be combined to produce the single value output. */

  char **FluxComponentNames;
  /* Relevant if UsesFlux == True.  Array list of eqnset-understood
     state vector component names.  These flux components will be
     combined to produce the single value output. Read and written. */

  real *FluxComponentWeights;
  /* Relevant if UsesFlux == True. Vector (size nFluxComponent) of
     weights used in linearly-combining the flux components to produce
     the output. Read and written. */

  int *FluxComponentMoments;
  /* Relevant if UsesFlux == True. Vector (size nFluxComponent) of
     integers between -1 and 2 that signify which position component
     should be used to multiply the weighted flux component terms.  -1
     means no weight by position, while 0,1,2 correspond to x,y,z,
     respectively. Read and written. */

  real LineCoord[6];
  /* Relevant only for LineIntegral output: x0,x1,y0,y1,[z0,z1]. */

  real CutPlane[4];
  /* Relevant only for CutPlaneIntegral: coefficients (a,b,c,d) in
     plane equation: ax + by + cz + d = 0 */

  int nBFG;
  /* For BoundaryIntegral output, number of boundary face groups
     (BFGs) associated with this output.  Read and written. */

  char **BFGTitles;
  /* For BoundaryIntegral output, array list of strings specifying the
     titles of the boundary face groups associated with the output.
     Read and written. */

  char *DumpFile;
  /* Point distribution data from BoundaryIntegral outputs can be
     dumped to a file.  DumpFile (if not NULL) stores the name of the
     file to which the integrand point values are written.  Currently
     only supported for UsesFlux=True (which includes pressure forces,
     for example).  Read and written */

  /*----------------*/
  /* For PointValue */
  /*----------------*/

  int egrp, elem;
  /* Element information for PointValue output. */
  
  real xref[3];
  /* Element reference coordinate for PointValue output. */

  /*---------------*/
  /* For SumOutput */
  /*---------------*/

  int nSumOutput;
  /* Number of outputs summed for a combination "SumOutput". */

  char **SumOutputNames;
  /* Names of outputs included in the sum [nSumOutput].*/
  
  real *SumOutputWeights;
  /* The output sum can be weighted; these are the weights
     [nSumOutput]. */
  
  real *SumOutputErrTols;
  /* Stores the error tolerances for each output when adapting the 
   mesh based on error estimates */

  /*----------------------*/
  /* For Unsteady Outputs */
  /*----------------------*/

  real StartTime;
  /* Time at which the unsteady output measurement begins. For point
     value outputs, this is the time at which the single measurement
     is taken. */

  real EndTime;
  /* Time at which the unsteady output measurement ends.  Disregarded
     for point value outputs. */
  
  /*------------------------*/
  /* For output sensitivity */
  /*------------------------*/
  
  int nSensitivity;
  
  xf_Sensitivity *Sensitivity;

  /*---------------------*/
  /* Internal structures */
  /*---------------------*/
  
  real Value;
  /* Stores latest-calculated value for this output. Not read or
     written. */

  real ErrEst;
  /* Stores latest-calculated error estimate for this output. Not read
     or written. */

  xf_CutPlaneIntersect *CutPlaneIntersect;
  /* Stores intersection information between cut plane and mesh
     elements; used for performing integrals on the cut plane.  Not
     read or written. */

  int *elemLocal;
  /* Processor-local element information for PointValue output.
     Size-2 vector corresponding to (egrp,elem).  Not read or
     written. */
  
}
xf_Output;


/* Output evaluation structure used to pass information related to
 output calculation. */
typedef struct
{
  real      *Value;       // pointer to value to calculate
  xf_Vector *Value_U;     // pointer to linearization vector of value
  xf_Vector *Value_G;     // pointer to GCL linearization vector of value
  real      *EV_U;        // pointer to linearization on one element
  real      *EV_G;        // pointer to GCL linearization on one element
  real      *FluxWeights; // weights in summation of output values
  int       *FluxMoments; // moment indicators for output summation
  FILE      *fidDump;     // if not NULL, file to which point data should be dumped
}
xf_OutputEvalData;



#endif // end ifndef _xf_OutputStruct_h
