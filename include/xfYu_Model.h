/*header for xfYu_Model.h */

#ifndef _xfYu_Model_h
#define _xfYu_Model_h 1

#include "xf.h"
#include "xf_MeshStruct.h"
#include "xf_GeomStruct.h"
#include "xf_DataStruct.h"
#include "xf_ParamStruct.h"
#include "xf_EqnSetStruct.h"
#include "xfYu_Statistics.h"
#include "xf_AllStruct.h"
#include "xfYu_ModelStruct.h"

/*
enum BoundaryConditionType{
     fullBC,
     noneBC,
     subinflowBC,
     staticPBC,
     noslipwallBC,
     slipwallBC,
     farPBC,     //characterisic based non-reflecting boundary;or pressure at far field(fluent)
     exactBC     //useful for some test
};

enum TimeIntegrationScheme{
     FE,
     RK45,
     SSPRK23,
     SSPRK34
};

//specification for some test cases
#define LinearCase xfe_False
#define Yu_AdvDiffAnnulus 101
#define Yu_AdvDiffRectChannel 102
#define AdvSpd 1.0
#define LapVis 0.001
#define Yu_TaylorGreenVortex 10
#define Yu_DoubleMachReflection 11
#define Yu_CompressVortexAdvect 12
#define Yu_PoiseuilleFlow 13

typedef struct
{
   int  dim;

   //whether normal quads or rectangle
   enum xfe_Bool TwistFlag;

   //basis type
   char basis[50];

   //time stepping scheme
   int typeTimeScheme;

   //basis order
   int  order;

   //number of variables
   int  nVars;

   //name of variables
   char **nameVars;

   //initial state of variables
   real *initVars;

   //number of boundary conditions
   int  nBCs;

   //name of boundary conditions
   char **nameBCs;

   //type of boundary conditions
   int  *typeBCs;

   //parameters of boundary condition
   real **paraBCs;

   //whether specify inflow boundary function
   enum xfe_Bool BCFunSpf;   //user specification on code

   //whether use limiter
   enum xfe_Bool LimiterFlag;

   //whether enable entropy bounding
   enum xfe_Bool EntropyBdFlag;

   //Vector point for MinEntropy
   xf_Vector *EntropyVec;

   //Vector related to flexible time stepping
   xf_Vector *MaxCharSpeed;
   xf_Vector *MinFaceLen;

   //Initial Data need to be interpolated
   int InterpolateInitData;

   //M value for TVB limiter
   real LimM;

   //number of negative pressure checking points
   int Num_negPckpnt;
   
   //coordinate of negative pressure checking points
   real *Coord_negPckpnt;

   //basis fcn coeff of negative pressure checking points
   xf_BasisData *Phi_negPckpnt;
   
   //whether let gamma (heat capacity ratio) varying
   int GammaVaryFlag;

   //initial value of gamma
   //for unvariable gamma: init_value = always value
   real GammaInit;

   //gamma difference threshold for face flux evaluation twice
   real Gammathreshold;

   //whether deal with reaction
   enum xfe_Bool ChemSource;

   //whether use detailed chemistry
   enum xfe_Bool DetailChem;

   //time step for current time; useful with detailed chemistry
   real dt_size;

   //user-defined file for data loading and initialization
   char DataFileName[64];

   //user-defined file: number of rows
   int DataFileNumRow;

   //user-defined file: number of cols
   int DataFileNumCol;

   //user-defined file: locus of front (flame or shock)
   real DataFileLocus;

   //provide default constant value for transport
   //if detailed chemistry is provoked, we use chemkin for computation
   //flag to diffusion activation
   enum xfe_Bool DiffFlag;

   //whether use Sutherland's law
   enum xfe_Bool Sutherland;

   //universal gas constant
   real Ru;

   //viscosity mu, unit: kg/m*s
   real mu_c;

   //heat conductivity, unit: W/m*K
   real kappa_c;

   //diffusivity kept the same for different species, unit: m^2/s
   real Diff_c;

   //molecular weight for simple chemistry model
   real moleW[50];

   //flag to activate convergence study
   enum xfe_Bool ConvergStudyFlag;

   //flag for mesh adaptation
   enum xfe_Bool MeshAdapt;

   //output for post-processing or logging
   int nOutput;

   Yu_Output *Output;
} 
Yu_Model;

typedef struct{
    
   //number of faces
   int nface;

   //elem group index of element:left-right-up-down
   int adjEgrp[4];

   //elem index in corresponding group of element
   int adjEelm[4];
   
   //type of element: 0/interior; 1/wallBC; 2/inletBC; 3/outletBC
   int adjEtyp[4];

   //data for the variable states in neighborhoods
   real *adjVars[4];

   //data for the variable states in self
   real *selfVars;

   //data for variable moment in order spatial-dimen * nVars 
   //(c_00, c_10, c_01, c11 for P1 || c_ii for P2)
   real *moment;

   //grid size x direction
   real deltax;
   
   //grid size y direction
   real deltay;

   //the tranformation matrix for first moments
   real signTranMatr[4];

   //the tranformation matrix back to first moments
   real signBackMatr[4];
}
Yu_Limiter;
*/

int PullinModel(Yu_Model *Model);

int DestroyModel(Yu_Model *Model);

int
ElemQuadOrder(int Order, int *QuadOrder);

int
FaceQuadOrder(int Order, int *QuadOrder);

int 
ConvFluxInterior(const int nq, const int sr, const int dim,
                 const real *U, real *F, const real Gamma);
int
ChemicalSource(const int nq, const int sr, const int dim, 
               const real *U, real *xglob, real *S, const real Gamma,
               enum xfe_Bool DetailChem, const real dt);
int
SpongeSource(int nq, int sr, int dim, real *xglob, real *U, real *S, 
             real gamma, real dt, real tcurr, real SpongeMeasure, Yu_Sponge *SP);

int
ConvFluxInterior_WENOlimiter(const int nq, const int sr, const int dim,
                 const real *Uxlim, const real *Uylim, real *F, const real Gamma);
int
ConvFluxInteriorFace(const int nq, const int sr, const int dim, 
                     const real *qUL, const real *qUR, const real *qn, 
                     real *qF, const real Gamma, real *MaxCharSpeed);
int
ConvFluxBoundaryFace(const int nq, const int sr, const int dim, 
                     const real *qUL, const real *qUR, const real *qn, 
                     real *qF, const real Gamma, real *MaxCharSpeed);

int
ConvFluxBoundaryState(const int nq, Yu_Model *Model,  const real *qUI, real *qUB,
                      const real *qn, const real *xglob, const real Time, const int localBCflag, 
                      const real Gamma);

int
ConvFluxBoundaryFace(const int nq, const int sr, const int dim, const real *qUI, const real *qUB,
                  const real *qn, real *qF, const real Gamma, real *MaxCharSpeed);

int
InitializeDataInterpolation(Yu_Model *Model, const int nq, const int sr, const int dim, 
                            const real *qxglob, real *qU, const real Gamma);

int
Yu_GammaVectorCreate(xf_All *All, Yu_Model *Model, xf_Vector *State, const char Name[], const real InitVal);

int
Yu_GammaVectorUpdate(xf_All *All, Yu_Model *Model, xf_Vector *Gamma, xf_Vector *State);

int
Yu_GhostStateConstruct(const int nq, const int sr, const int dim, const real *quL, const real *quR,
                       real *quLghost, real *quRghost, const real GammaL, const real GammaR);

int
Yu_EnergyCorrection(xf_All *All, Yu_Model *Model, xf_Vector *S);

int
ReadDataFromUserDefineFile(Yu_Model *Model, const char *FileName);

void 
DeleteFileDataAllocation();

int
VisFluxBoundaryState(const int nq, const int sr, const int dim, const real *wn,
                     const real *xglob, const real *qgUI, real *qgUB, enum xfe_Bool *VisFluxFlag,
                     int *SetIndex, const int localBCtype, const real *localBCpara, enum xfe_Bool AVFlag);

int
EvaluateError(xf_All *All, xf_Vector *U, enum xfe_Bool PrintFlag);

int 
Yu_UserDefinedBCparamsUpdate(const real t, const real dt_size);

int
init_data_read_from_file(Yu_datafile **dataholder);

int
DissipationCorrect(const int nq, const int sr, const int dim, const real *qUL, const real *qUR,
                   const real *qn, const real Gamma, real *qQn, real *qQn_dp, 
                   const real eta);
#endif
