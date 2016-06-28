/*header for xfYu_Model.h */

#ifndef _xfYu_Model_h
#define _xfYu_Model_h 1

#include "xf.h"
#include "xf_MeshStruct.h"
#include "xf_GeomStruct.h"
#include "xf_DataStruct.h"
#include "xf_ParamStruct.h"
#include "xf_EqnSetStruct.h"

//for dynamic-p adaptation
#define NUM_VAR_DYM_P 4

enum BoundaryConditionType{
     fullBC,
     noneBC,
     subinflowBC,
     staginflowBC,
     staticPBC,
     noslipwallBC,
     slipwallBC,
     isothermwallBC, //noslip wall with fixed temperature   
     farPBC,         //characterisic based non-reflecting boundary;or pressure at far field(fluent)
     exactBC         //useful for some test
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
#define Yu_ForwardFacingStep 14
#define Yu_Kelvin_Helmholtz_Instability 15
#define Yu_Density_Wave_ConstVel 16
#define Yu_ShockVortexInteraction 17
#define Yu_2DRiemannProblem2 18
#define Yu_2DRiemannProblem3 19
#define Yu_2DRiemannProblem4 20
#define Yu_2DRiemannProblem5 21
#define Yu_2DSodTube 22
#define Yu_quasi1D_Detonation 23
#define Yu_acoustic_ML 24
#define Yu_acoustic_pulse 25
#define Yu_plane_acoustic_test 26
#define Yu_channel_flow 27
#define Yu_BCstudy_acoustic_pulse 28
#define Yu_BCstudy_vortex_transport 29
#define Yu_debug 30
//quantities for output statistics
enum OutputQuantity{
   rho_bar,
   rho_rms,
   u_bar,
   u_rms,
   v_bar,
   v_rms,
   w_bar,
   w_rms,
   p_mean,
   RS_uv,
   RS_uw,
   RS_vw,
   T_mean,
   T_rms,
   kinenergy,
   enstrophy,
   rhou,
   rhov,
   rhow,
   rhoE,
   OutputQuantity_last  
};

static char *OutputQuanityName[OutputQuantity_last] = {
   "rho_bar",
   "rho_rms",
   "u_bar",
   "u_rms",
   "v_bar",
   "v_rms",
   "w_bar",
   "w_rms",
   "p_mean",
   "RS_uv",
   "RS_uw",
   "RS_vw",
   "T_mean",
   "T_rms",
   "kinenergy",
   "enstrophy",
   "rhou",
   "rhov",
   "rhow",
   "rhoE",
};

//integral output for logging and post-processing
//no need to be sophisticated as Kfid used
//from now, only consider one output for each type
//type: domain-integral 
//      point-wise
//      boundary-integral
typedef struct
{
   char Name[30];
   /*identity of an output instance*/
   
   enum xfe_OutputType Type;
   /*see xf_OutputStruct.h for detail*/

   int nVars;     //number of different output quantities
   int *IVars;    //specified indice of variable in my code

   /*for point-type of output*/
   int nPoints;   //number of points for pointwise output
   enum xfe_Bool RestartPointStat; //read restart file for previous collection
   real pretime;  //the previous time level for statistics purpose
   real accumtime;//time accumulated at restart file
   int *egrp;     //an array for all output points; might not needed?
   int *elem;     //an array for all output points; might not needed?
   real *xref;    //coordinates array for all output points
   int *egrpLocal;//an array for all output points
   int *elemLocal;//an array for all output pionts
   real *PointData;    //holding all the data 

   real Value;
   /*store the latest output value*/
   
   //int nBFG;
   char BFGTitles[30];
   enum xfe_Bool UsesFlux;//current only support the case when this is True
   int nFluxComponent;
   char **FluxComponentNames;//can be replaced by the following entry
   int *FluxComponentRanks;
   real *FluxComponentWeights;
   int *FluxComponentMoments;//this is not useful
   /*for boundary integral type of output; see xf_OutputStruct.h for specification details*/

   real SampleTimeInv;    //time interval for data sample 
   int OutputFile_Freq_Ratio_DataFile;
   enum xfe_Bool SequanceDump;
   int File_Write_Offset;

   //3D self-similiar averaging
   int threeDSSMod;     //'1' rotating; '2' translating
   real threeDSSParam[3];    //if '1' param: org (y, z) and x-axis; if '2' param: x_str x_end x-axis
}
Yu_Output;

//class for initialization with input data file
typedef struct{
   char Name[60];

   real *data;

   int rank;
   int dim[3];             //number of data along x y z
   real xyz_min_max[3][2]; 
}
Yu_datafile; 

//dynamic p adaptation parameter specification 
typedef struct
{
   int MaxOrder;
   int MinOrder;
   int spongeOrder;
   real refinefrac;
   real coarsenfrac;
   real Yu_adapt_time_size;   //time slot size 
   real xfa_file_write_inv;   //time interval for .xfa file dumping
   real Yu_load_balance_time_size; //time slot size for checking and conducting load balance
   
   //performance record for load balancing
   real best_time;
   real best_imp;   //load imbalance metric
   int adapt_index; 

}
Yu_Dyn_p_Adapt_param;

//sponge parameters
typedef struct
{
   real State_Inf[5];
   real line[4*3];
   real point[4*3];
   enum xfe_Bool WheMultiStateSponge;  
}
Yu_Sponge;

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

   //whether enable entropy bounding
   enum xfe_Bool ConstEntropyBdFlag;
   
   //Vector point for MinEntropy
   xf_Vector *EntropyVec;

   //Vector related to flexible time stepping
   xf_Vector *MaxCharSpeed;
   xf_Vector *MinFaceLen;

   //indicator for Riemann solver
   xf_Vector *RiemannIndicator; 

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

   //whether use some special source for MMS or channel flow
   enum xfe_Bool MMSSource;

   //time step for current time; useful with detailed chemistry
   real dt_size;

   //user-defined file for data loading and initialization
   Yu_datafile *init_datafile;

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
   real moleW[30];

   //flag to activate convergence study
   enum xfe_Bool ConvergStudyFlag;

   //data for artificial viscosity
   enum xfe_Bool AVmodel;
   xf_Vector *AVmodel_data;

   //output for post-processing or logging
   int nOutput;

   Yu_Output *Output;

   //parameter for adaptation
   enum xfe_Bool Stat_h_Adapt;               //h-refinement
   enum xfe_Bool Dyn_p_Adapt;                //p-adaptation
   Yu_Dyn_p_Adapt_param *Dyn_p_Adapt_param;  //parameters wrapped in struct

   xf_Vector *Sponge;       //point to sponge specification
   Yu_Sponge *SpongeParam;  //parameters related to sponge setup

   real cputime_indicator;  //estimate the cpu time for each process

   real time_fac;  // a time factor to scale time-step-size

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
  
#endif
