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

#ifndef _xf_BasisFcn_h
#define _xf_BasisFcn_h 1

/*
  FILE:  xf_BasisFcn.h

  This file contains the headers for functions that return basis functions

*/

/******************************************************************/
//   FUNCTION Prototype: xf_Shape_TriLagrange
extern int 
xf_Shape_TriLagrange(int p, const real *xy, real *phi);
/*
PURPOSE:

  Returns TriLagrange basis of order p at point xy.

INPUTS:

  p : desired order of basis
  xy: point at which to evaluate

OUTPUTS: 

  phi : basis functions [nn(p)]

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Shape_TriHierarch
extern int 
xf_Shape_TriHierarch(int p, const real *xy, real *phi);
/*
PURPOSE:

  Returns TriHierarch basis of order p at point xy.

INPUTS:

  p : desired order of basis
  xy: point at which to evaluate

OUTPUTS: 

  phi : basis functions [nn(p)]

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Shape_TetLagrange
extern int 
xf_Shape_TetLagrange(int p, const real *xyz, real *phi);
/*
PURPOSE:

  Returns TetLagrange basis of order p at point xyz.

INPUTS:

  p : desired order of basis
  xyz: point at which to evaluate

OUTPUTS: 

  phi : basis functions [nn(p)]

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Shape_TetHierarch
extern int 
xf_Shape_TetHierarch(int p, const real *xyz, real *phi);
/*
PURPOSE:

  Returns TetHierarch basis of order p at point xyz.

INPUTS:

  p : desired order of basis
  xyz: point at which to evaluate

OUTPUTS: 

  phi : basis functions [nn(p)]

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_Grad_LineLagrange
extern int 
xf_Grad_LineLagrange(int p, const real x, real *gphi);
/*
PURPOSE:

  Returns LineLagrange basis gradients of order p at point x.

INPUTS:

  p : desired order of basis
  x: point at which to evaluate

OUTPUTS: 

  gphi : basis functions [nn(p)*1]. 

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Grad_TriLagrange
extern int 
xf_Grad_TriLagrange(int p, const real *xy, real *gphi, int n);
/*
PURPOSE:

  Returns TriLagrange basis gradients of order p at point xy.

INPUTS:

  p : desired order of basis
  xy: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim].  The dim individual gradients
         are offset by n in memory

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Grad_TriHierarch
extern int 
xf_Grad_TriHierarch(int p, const real *xy, real *gphi, int n);
/*
PURPOSE:

  Returns TriHierarch basis gradients of order p at point xy.

INPUTS:

  p : desired order of basis
  xy: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim].  The dim individual gradients
         are offset by n in memory

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Grad_TetLagrange
extern int 
xf_Grad_TetLagrange(int p, const real *xyz, real *gphi, int n);
/*
PURPOSE:

  Returns TetLagrange basis gradients of order p at point xyz.

INPUTS:

  p : desired order of basis
  xyz: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim].  The dim individual gradients
         are offset by n in memory

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Grad_TetHierarch
extern int 
xf_Grad_TetHierarch(int p, const real *xy, real *gphi, int n);
/*
PURPOSE:

  Returns TetHierarch basis gradients of order p at point xyz.

INPUTS:

  p : desired order of basis
  xyz: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim].  The dim individual gradients
         are offset by n in memory

RETURN:

  Error Code
*/



/******************************************************************/
//   FUNCTION Prototype: xf_Hess_TriLagrange
extern int 
xf_Hess_TriLagrange(int p, const real *xy, real *gphi, int n);
/*
PURPOSE:

  Returns TriLagrange basis Hessians of order p at point xy.

INPUTS:

  p : desired order of basis
  xy: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim*dim].  The dim*dim individual
         Hessians are offset by n in memory

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Hess_TriHierarch
extern int 
xf_Hess_TriHierarch(int p, const real *xy, real *gphi, int n);
/*
PURPOSE:

  Returns TriHierarch basis Hessians of order p at point xy.

INPUTS:

  p : desired order of basis
  xy: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim*dim].  The dim*dim individual
         Hessians are offset by n in memory

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Hess_TetLagrange
extern int 
xf_Hess_TetLagrange(int p, const real *xyz, real *gphi, int n);
/*
PURPOSE:

  Returns TetLagrange basis Hessians of order p at point xyz.

INPUTS:

  p : desired order of basis
  xyz: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim*dim].  The dim*dim individual
         Hessians are offset by n in memory

RETURN:

  Error Code
*/

/******************************************************************/
//   FUNCTION Prototype: xf_Hess_TetHierarch
extern int 
xf_Hess_TetHierarch(int p, const real *xy, real *gphi, int n);
/*
PURPOSE:

  Returns TetHierarch basis Hessians of order p at point xyz.

INPUTS:

  p : desired order of basis
  xyz: point at which to evaluate
  n : linear storage offset between _X, _Y, etc.

OUTPUTS: 

  gphi : basis functions [nn(p)*dim*dim].  The dim*dim individual
         Hessians are offset by n in memory

RETURN:

  Error Code
*/




#endif // end ifndef _xf_BasisFcn_h
