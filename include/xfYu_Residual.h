 #ifndef _xfYu_Residual_h
 #define _xfYu_Residual_h 1

#include "xf_SolverStruct.h"

/******************************************************************/
//   FUNCTION Prototype: xf_CalculateResidual
extern int
xfYu_CalculateResidual(xf_All *All, Yu_Model *Model, xf_Vector *U, xf_Vector *R,
                     xf_JacobianMatrix *R_U, xf_SolverData *SolverData);

#endif
