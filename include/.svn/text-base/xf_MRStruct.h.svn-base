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

#ifndef _xf_MRStruct_h
#define _xf_MRStruct_h 1

/*
  FILE:  xf_MRStruct.h

  This file contains structures used in Model Reduction

*/


/* Interpolation Point structure*/
typedef struct
{
  int M; // number of points
  int Dim; // dimension
  int *proc; // processor numbers of points
  int *egrp; // element groups of points
  int *elem; // element numbers of points
  int *node; // Quad node numbers of points
  real *xref; // ref coords in each elem [M*Dim]
}
xf_IPointType;


/* Reduced model structure */
typedef struct
{
  int N;  // # of basis functions
  real *A; // reduced matrix for linear terms [N by N]
  real *L; // constant term corresponding to linear reduction [N]

  int nOutput; // # of outputs
  real *F;  // matrix for constructing reduced output [nOutput x N]
  real *F0; // constant vector for output calc  [nOutput x 1]

  int nNonLinear; // # of nonlinear terms
  int *M; // # of coeff expansion functions for each nonlinear term
  real **E; // E matrices for each nonlinear term [N by M]
  real **D; // D matrices for each nonlinear term [M by N]
  xf_IPointType *z; // interpolation point sets for each nonlinear term
  xf_ResTerm *ResTerm; // Residual term for each nonlinear term

  xf_EqnSet *EqnSet; // equation set structure

  // Potentially useful for Newton restart; not required for online
  int nSnap; // # of snapshots (not read or written to binary)
  real *UN; // matrix of snapshots dotted with basis vectors [nSnap by N]
            // not read or written to binary
}
xf_ReducedModel;

 

#endif // end ifndef _xf_MRStruct_h
