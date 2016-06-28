
#ifndef _xfYu_Limiter_h
#define _xfYu_Limiter_h 1

int
xf_FullinLimiterStruct(xf_All *All, Yu_Model *Model, Yu_Limiter ***pLimiter);

int
DestroyLimiterStruct(xf_All *All, Yu_Limiter **Limiter);

int
Yu_ConductLimiting(xf_All *All, Yu_Model *Model, Yu_Limiter **Limiter, xf_Vector **pU);

#endif
