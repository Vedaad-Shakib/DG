
#ifndef _xfYu_AdaptSolver_h
#define _xfYu_AdaptSolver_h 1

#include "xf_SolverStruct.h"
#include "xf_SolverTools.h"
#include "xf_AdaptStruct.h"

extern int 
xfYu_ApplyUnsteadyAdapt(xf_All *All, Yu_Model *Model, Yu_Limiter ** Limiter, 
                        const char *SavePrefix, enum xfe_Bool RestartFlag,
                        xf_Vector *U0, xf_TimeHistData *TimeHistData);

//int
//SpongeSource(int nq, int sr, int dim, real *xglob, real *U, real *S, real dt);



#endif

