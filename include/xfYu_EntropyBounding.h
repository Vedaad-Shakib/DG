/* header for xfYu_EntropyBounding.c */


#ifndef _xfYu_EntropyBounding_h
#define _xfYu_EntropyBounding_h 1

int 
Yu_ConductEntropyBounding(xf_All *All, Yu_Model *Model, xf_Vector **pU,
                          enum xfe_Bool *CtrlSeq);
int 
DestroyEntropyBoundStruct();

int
Yu_MinFaceLengthScale(xf_All *All, Yu_Model *Model, xf_Vector **pU);

#endif
