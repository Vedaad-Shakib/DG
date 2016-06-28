//my script for post-processing data
//
//
//

#ifndef _xfYu_Statistics_h
#define _xfYu_Statistics_h 1

enum Yu_ScalarType{
     Yu_Pressure,
     Yu_VorticityMag,
     Yu_KineticEnergy,
     Yu_Enstrophy,
     Yu_Dilatation,
     Yu_Mach,
     Yu_Entropy,
     Yu_ScalarLast
};
static char *Yu_ScalarName[Yu_ScalarLast] = {
     "Pressure",
     "VorticityMag",
     "KineticEnergy",
     "Enstrophy",
     "Dilatation",
     "Mach",
     "Entropy"
};


//member functions
int 
Yu_InitOutput(Yu_Output *Output);

int
Yu_VisualizationOutput(const int sr, const int dim, const int nq,
                       const real *qU, const real *qgU, const real *qout);

int
Yu_StatisticsOutput(xf_All *All, const char *Name, const xf_Vector *U,
                    real *Value, xf_Vector *Value_U, enum xfe_AddType AddFlag1);

int
Yu_OutputPointValueDump(xf_All *All, Yu_Output *Output);
#endif

