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

#ifndef _xf_ParamDefault_h
#define _xf_ParamDefault_h 1

/*
  FILE:  xf_ParamDefault.h

  This file contains the default parameter keys and values

*/

/* Key definitions for key/value parameters */

static char *xf_DefaultParamList[] =
{

  /*==================*/
  /* Case Information */
  /*==================*/

  "Restart", "False",
  /* True to restart from an existing solution (.xfa file).  If False,
     any existing solution in the InputFile is overwritten. */

  "InputFile", "None",
  /* Name of the input file (e.g. .xfa, .gri, .msh). */

  "SavePrefix", "None",
  /* Root name of all files saved to disk (suffixes will be added automatically). */

  "Default", "False",
  /* If True, parameters in the input file will be ignored, and those
     parameters not specified in the job file will take on default
     values. */

  "EqnSetFile", "None",
  /* Name of the equation set file to use, including suffix. */

  "InterpBasis", "TriLagrange",
  /* Basis used for interpolating the solution.  See list of
     definitions in xf_BasisStruct.h */

  "InterpOrder", "0",
  /* Spatial order for interpolating the solution. Specify as
     "Variable" in order to use existing order information for a
     restart run. */

  "SequenceOrderAdd", "0",
  /* For order sequencing in steady runs.  The final solution order
     will be InterpOrder (or current order) + SequenceOrderAdd.  The
     solution will start at InterpOrder and work up to the final
     order. */

  "iRestart", "0",
  /* The restart index of the current case.  Incremented automatically
     upon each run. */

  "DataMerge", "None",
  /* Space-separated list of .data files (including suffix) which will
     be read and merged with the input file. */

  "ScaleStateUsingIC", "False",
  /* If True, and if Restart is True, the existing state vector will be
     rescaled in an eqnset-specific way before the calculation begins.
     This rescaling is based on the new versus original initial
     conditions.  Useful for parameter sequencing.  */



  /*============================*/
  /* General Solver Information */
  /*============================*/

  "NonlinearSolver", "Newton",
  /* Requested nonlinear solver (see xf_SolverStruct.h). */  

  "LinearSolver", "GMRES",
  /* Requested linear solver (see xf_SolverStruct.h). */
  
  "Preconditioner", "BlockJacobi",
  /* Preconditioner for the linear solver (see xf_SolverStruct.h). */

  "iIterNonlinear", "0",
  /* Current nonlinear iteration number, for restart purposes. */

  "nIterNonlinear", "0",
  /* Maximum number of nonlinear iterations requested. */
  
  "nIterLinear", "0",
  /* Maximum number of linear iterations requested.*/

  "ReusedU", "False",
  /* If True, dU is used as initial guess in the linear solver. */
  
  "nSmoothLinear", "2",
  /* Number of linear smoothing steps to take during each application
     of the preconditioner in GMRES. */
 
  "nIterGMRESOuter", "15",
  /* Number of GMRES outer restart iterations. */

  "nIterGMRESInner", "50",
  /* Number of GMRES inner iterations (= number of storage vectors
     required). */

  "CheckHaltLinear", "True",
  /* If True, checks for existence of a user-placed STOP file (a file
     named STOP, with STOP written on the first line) during linear
     iterations. */

  "ResidualTolerance", "1e-12",
  /* Execution stops when the L1 residual norm drops below this value
     (absolute).*/

  "MinUpdateFraction", "1e-2",
  /* Minimum update fraction that is still considered acceptable;
     below this value, the iteration will be re-done with a lower
     CFL.*/

  "LocalStateUpdate", "False",
  /* Allow solution to proceed even if MinUpdateFraction is violated
     on some elements; on those troublesome elements, take local small
     updates. */

  "MaxQuadFlag", "False",
  /* Set to True to use highest quad rule everywhere. Not currently
     supported. */

  "ResidualCriterionFlag", "True",
  /* If True, number of linear iterations will depend on a relative
     residual criterion. */
  
  "MinLinResDecreaseFactor", "1e-3",
  /* Minimum amount by which to decrease the linear residual if using
     ResidualCriterionFlag = True. */

  "SortLines", "True",
  /* If True, lines will be sorted using residual, for LineGS preconditioning. */

  "PingResidual", "False",
  /* If True, residual will be pinged.  This is a lengthy testing
     procedure in which the residual derivatives are checked against
     finite differences.  */

  "PingOutput", "None",
  /* If not None/Null, the specified output will be pinged -- that is,
     the derivatives will be checked via finite differences. */

  "LineViz", "False",
  /* If True, a LineID vector will be created for visualizing lines
     used in the line-based preconditioner. */

  "ProcViz", "False",
  /* If True, a ProcID vector will be created for visualizing
     processor partitions. */

  "WriteResidual", "False",
  /* If True, the residual vector will be made writable, e.g. for
     plotting at the end of the simulation. */

  "WriteUpdate", "False",
  /* If True, the update vector will be made writable. */

  "WriteLog", "True",
  /* If True, <SavePrefix>.log file is written. */

  "DumpSystem", "False",
  /* If true, system (mass and Jacobian matrices) will be written to
     text files (M.txt and A.txt).  Only supported in serial.  Should
     only use for small-ish problems. */

  "BLASLibName", "None",
  /* Name/path to BLAS library.  No longer supported. */

  "LinearFlag", "False",
  /* If true, system is assumed to be linear and the Jacobian is not
     recomputed during the "nonlinear" solve, assuming it exists --
     i.e. the Jacobian will be computed on first call, but will not be
     recomputed on subsequent calls (e.g. during time-stepping). */

  "MeshIsLinear", "False",
  /* Irrelevant for triangles and tets. Quads and hexes are generally
     considered nonlinear even for Q=1, because the Jacobian inside is
     not constant.  This flag overwrites this assumption, so a user
     can run on Cartesian meshes without storing the mass matrix. */

  "InterpolateIC", "False",
  /* If true, functional initial conditions (ICs) are interpolated via
     pointwise sampling instead of performing a QR least-squares
     minimization.  Note, for this interpolation, a Lagrange basis
     must be used.  On the other hand, QR least squares can take any
     basis functions. */

  "SteadyWriteInterval", "-1",
  /* For steady runs, how often (e.g. Newton iterations) the state
     vector should be written to disk.  If negative, state is not
     written during the course of the steady solve. */

  /*==================================*/
  /* CFL and Time-Stepping Parameters */
  /*==================================*/

  "CFL", "1.0",
  /* Courant-Friedrichs-Lewy number: used for improving robustness via
     artificial time stepping during Newton non-linear solve. */

  "CFLMax", "1e30",
  /* Maximum CFL: during the course of a regular pseudo-time
     continuation nonlinear solve, CFL is increased to take larger and
     larger artificial time steps.  This is a bound on how high CFL is
     taken. */

  "CFLMin", "1e-3",
  /* Minimum CFL: if a nonlinear solve is not going well (small update
     fractions, for example), the CFL is decreased according to a
     procedure outlined in the solver code.  If the CFL drops below
     CFLMin, the solver gives up and returns a non-convergence
     error. */
  
  "CFLEvolution", "Exp",
  /* Type of CFL evolution strategy to be used, e.g., SER-A, SER-B,
   Exponential, etc. These correspond to enumerator types in
   xf_SolverStruct.h. */

  "CFLDecreaseFactor", "10",
  /* Factor by which to decrease the CFL when problems are
     encountered during a nonlinear solve. */

  "CFLIncreaseFactor", "2",
  /* Factor by which to increase the CFL when the nonlinear solve is
     going well. */

  "LocalTimeStepping", "True",
  /* If True, the artificial time step is different for each element
     and is set according to the CFL combined with the element size
     (hydraulic diameter). If False, a uniform time step is used for
     all elements. */

  "UseConstDt", "True",
  /* For unsteady runs, UseConstDt=True means that a constant time
     step based on the start/end time and number of time steps is
     used.  If UseConstDt is False (and if the solver supports this
     capability), the time step is set using the CFL and could be
     different at each unsteady iteration.  */

  "TimeScheme", "Steady",
  /* For unsteady runs, set this to an unsteady scheme (e.g. BDF, DG1,
     DIRK3, etc.). See xf_SolverStruct.h */

  "Time", "0.0",
  /* For unsteady runs, the current time in the simulation.  When set
     in the .job file, this indicates the start time. */

  "EndTime", "1.0",
  /* For unsteady runs, the end time of the simulation. */

  "nTimeStep", "10",
  /* For unsteady runs, number of time steps in the simulation. */

  "UnsteadyWriteInterval", "-1",
  /* For unsteady runs, how often, in terms of unsteady iterations,
     the state vector should be written to disk. Less than zero means
     never write. */

  "WriteOffset", "0",
  /*the output data index often needs an offset for restarting
   */
  
  "WriteAugFac", "1",
  /*specify a larger spacing for data saving, useful when there requests
   *for output data saving simultaneously
   */

  "TimeStepDecreaseFactor", "2.0",
  /* Factor by which the TimeStep is decreased when problems are
     encountered during a nonlinear unsteady solve. */

  "InputTimeHist", "None",
  /* Time history file (text) to read, if user wants to prescribe a
     custom time history.  The alternative default is equally-spaced
     time steps. */


  /*===================*/
  /* Output Parameters */
  /*===================*/

  "LogOutput", "None",
  /* List of outputs which should be logged, listed sequentially by
     output name (space as separator).  The outputs need to be defined
     in the .eqn file. */
  
  "Verbosity", "Medium",
  /* Set to Low for little/no output to screen.  Medium/High for more
     ouput. */



  /*=================*/
  /* Adjoint Solver  */
  /*=================*/

  "AdjointOutputs", "None",
  /* Set to output name (e.g. Drag) to enable adjoint calculation for
     that output.  Multiple outputs can be listed in the same string,
     separated by spaces. */
   
  "nIterAdjoint", "200",
  /* Number of iterations for the (linear) adjoint solve. */

  "AdjointResidualTolerance", "1e-12",
  /* Residual tolerance for each linear adjoint solve. */
  
  "WriteAdjoint", "True",
  /* If True, adjoint(s) will be written with data. */



  /*=======*/
  /* GMRES */
  /*=======*/

  "GMRES_MaxInner", "40",
  /* Maximum number of inner iterations for GMRES (i.e. iterations until a restart). */

  "GMRES_MaxOuter", "10",
  /* Maximum number of outer iterations for GMRES.  */



  /*=======*/
  /*  ILU  */
  /*=======*/

  "ILUOrdering", "MDF",
  /* Ordering algorithm for Incomplete LU factorization linear solver */

  /*===================*/
  /* Coarse Correction */
  /*===================*/
  
  "CoarseCorrectionFlag", "False",
  /* If True, a coarse-grid correction will be applied as a
     preconditioner for the linear solve.  Smoothing will also be
     performed according to the Preconditioner. */

  "CoarseCorrectionOrder", "0",
  /* Interpolation order for the coarse correction. */

  "CoarseCorrectionIter", "10",
  /* Number of coarse-grid iterations. */

  "CoarseCorrectionSmooth", "2",
  /* Number of post-correction smoothing iterations. */

  "CoarseCorrectionRelax", "0.66666",
  /* Under-relaxation factor for smoothing. */

  "CoarseCorrectionIsAdditive", "False",
  /* If True, additive Schwarz is used; else it is multiplicative. */


  /*===================*/
  /* PMG = P Multigrid */
  /*===================*/
  
  "PMG_Cycle", "VCycle",
  /* Multigrid cycle; see xf_SolverStruct.h for definitions. */

  "PMG_nLevel", "2",
  /* Number of p multigrid levels, if not specifying the levels
     directly. <=0 indicates that PMG_CoarseOrders is specified. The
     coarse orders will be p, p-1, p-2, etc. */

  "PMG_CoarseOrders", "None",
  /* List of coarse orders to use for pMG, starting with one below the
     finest, down to coarsest.  None indicates that pMG_nLevel should
     be used instead. */

  "PMG_nPreSmooth", "1",
  /* Number of pre-smoothing iterations. */

  "PMG_nPostSmooth", "1",
  /* Number of post-smoothing iterations. */

  "PMG_nCoarseIter", "10",
  /* Number of smoothing iterations on the coarsest level. */

  "PMG_CoarseLinearSolver", "GMRES",
  /* Type of linear solver to use on the coarsest MG level. */

  "PMG_CoarsePreconditioner", "LineJacobi",
  /* Type of preconditioner on the coarsest level. */

  "IterativeUnderRelax", "1.0",
  /* Under-relaxation factor for iterative stepping.  Usually should
     set to < 1 for p-multigrid runs. */


  /*============*/
  /* DG in Time */
  /*============*/

  "DGTimeUnderRelax", "1.0",
  /* Under-relaxation factor for iterative DG in time. */

  /*==============================================*/
  /* Line Search and Residual-Optimization Solver */
  /*==============================================*/
  
  "PenalizeResidual", "False",
  /* If True, the residual vector is taxed by (1+P(U)) */
  
  "StateUpdateMethod", "MaxPrimitiveChange",
  /* Type of solution update. Options: MaxPrimitiveChange, LineSearch */

  "UpdateFracReduction", "0.5",
  /* Used to decrease update fraction in line search. */
  
  "LineSearchGreedy", "False",
  /* If True, the Line-search uses a greedy algorithm. */

  /*======================*/
  /* Domain Decomposition */
  /*======================*/

  "DD_DoDomainBlockJacobi", "False",
  /* If True, parallel message exchange will not happen during the
     linear steps of the preconditioner application. */
  
  "DD_LoadBalance", "True",
  /* If True, the cpu load is rebalanced using elemental p-order 
   and connectivity of the line-based preconditioner */
  
  "DD_DeleteNonEssential", "True",
  /* If True, non-essential data in All->DataSet will be deleted 
   before comunication. This reduces useless communication.*/
  
  "DD_CPUWorkReport", "False",
  /* If True, each processor will write a report representative of 
   the amount of cpu-work in their partition.*/
  
  "DD_RepartEveryNdt", "-1",
  /* If negative or zero, the domain will not be repartitioned 
   during the calculation. Else, the domain will be repartitioned 
   every N time-steps. */
  
  /*=========================*/
  /* Discretization-Specific */
  /*=========================*/

  "DiffusionDiscretization", "BR2",
  /* Type of diffusion discretization to use.  See SolverStruct.h. */

  "ViscStabFactor", "2.0",
  /* Factor multiplying stabilization term, e.g. eta in BR2.  In BR2,
     the code sets eta to the number of faces times this factor. */

  /*==================*/
  /* Error Estimation */
  /*==================*/

  "ErrEstOrderIncrement", "1",
  /* Order increment added to current order to get fine space for
     error estimation.  */

  "ErrEstUseReconstruct", "False",
  /* If True, the fine space solution is obtained by reconstructing
     the low-order solution, not by solving or iterating on the fine
     space. */

  "FineSpace_nIterNonlinear", "5",
  /* Number of iterations of the nonlinear solver on the fine space. */
  
  "FineSpace_LinearSolver", "Iterative",
  /* Type of linear solver to use for the fine space. */

  "FineSpace_nIterLinear", "1",
  /* Number of iterations of the linear solver on the fine space.
     Note, if using a lean iterative nonlinear solver on the fine
     space, set this value to 1, as there is no speed advantage of
     larger numbers, and with 1 linear iteration, the update will be
     checked on every iteration. */

  "FineSpace_Preconditioner", "BlockJacobiLean",
  /* Type of preconditioner to use for the fine space. */

  "FineSpace_nIterAdjoint", "5",
  /* Number of adjoint iterations on the fine space.  Note, the
     adjoint solver is the same as the FineSpace_LinearSolver, as is
     the preconditioner. */

  "FineSpace_CFL", "5",
  /* Starting CFL number for solving on the fine space. */

  "FineSpace_AdjointResidualTolerance", "1e-12",
  /* Residual tolerance for each fine-space adjoint solve. */
  
  "ErrEstUsePrimal", "True",
  /* If true, primal residual form of error estimate will be used. */

  "ErrEstUseDual", "False",
  /* If true, dual residual form of error estimate will be used. Note,
     if both Primal and Dual are true, a 50/50 contribution split will
     be used. */

  "VolumeSpecificResidual", "False",
  /* If True, the residual adaptive indicator gets divided 
   by the volume/area of the element.  */

  "ErrEstBoundingBox", "None",
  /* xmin xmax [ymin ymax [zmin zmax]] coordinates of a bounding box
     used to restrict error estimation/adaptation.  Contributions to
     the error and to the error indicator are not computed from
     elements outside this bounding box if one is specified. */

  "ErrEstResBreakdown", "None",
  /* Type of residual to use for error-estimation breakdown
     (contribution due to just one residual term).  Possibilities are
     None, Convection, Diffusion, Source. */

  "ErrEstHRef", "False",
  /* If True, mesh will be uniformly h-refined prior to calculation of
     the error estimate.  The adjoint/state will not be iterated on
     the new mesh, so the $Psi^T$ * R error estimate should not change.
     However, the conservative sum of indicators will likely increase
     due to additional subelement contributions (uncovering more of
     the residual) */

  "ErrEstOutputs", "None",
  /* Set to output name (e.g. Drag) to enable adjoint-weighted error 
   estimation only, without adapting, for that output.  Multiple outputs 
   can be listed in the same string, separated by spaces. 
   Note: if output-based adaptation is ON, then error estimation for 
   outputs listed in AdaptOutputs will be performed.*/

  /*===========================*/
  /* Unsteady Error Estimation */
  /*===========================*/

  "UErrEstOn", "False",
  /* If True, error estimation is performed during the unsteady
     adjoint solve.  The outputs considered are all the ones for which
     an adjoint is computed. */

  "UErrEstnSubSlab", "1",
  /* Number of time sub-slabs when h-refining the slabs for a finer
     temporal space.  1 means no slab refinement. */

  "UErrEstOrderIncrement", "0",
  /* Order added to current temporal order to get fine space for error
     estimation. */

  "UErrEstIterative", "False",
  /* If True, the FineSpace_* variables are used in a solve that
     iteratively corrects a coarse-space unsteady adjoint. */

  "UErrEstIterativeGS", "False",
  /* Applies to iterative adjoint solve during error estimation.  If
     True, the adjoint solve proceeds in a Gauss-Seidel fashion
     backwards in time (more accurate, not computationally more
     expensive, but less robust because errors due to inexact solve
     could grow quickly).  If False, the adjoint iterative solve is
     Jacobi in time (less accurate, more robust). */

  "UErrEstUseReconstruct", "False",
  /* If True, the fine space temporal adjoint solution is obtained by
     reconstructing the low-order temporal solution, not by solving or
     iterating on the fine temporal space. */

  "UErrEstIndicatorWriteInterval", "-1", 
  /* Write interval (in terms of time slabs) for spatially-localized
     unsteady error indicator. */

  "UErrEstAnisoMeasure", "Proj",
  /* Space-time anisotropy measure for unsteady adaptation.  The
     default is based on projection of adjoint between space-time
     elements.  See xf_SolverStruct.h for others. */

  "UErrEstConservative", "False",
  /* If True, a conservative form of unsteady error estimation is used
     -- basically consisting of more absolute values. */


  /*============*/
  /* Adaptation */
  /*============*/

  "AdaptIter", "0",
  /* Number of adaptation iterations to run.  Greater than zero
     indicates adaptation is requested. */

  "AdaptOn", "Output",
  /* What drives the adaptation?  Default is output-based (Output).
     Alternatives include residual-based (Residual), penalty-based
     (Penalty), based on a scalar quantity (Scalar), uniform
     refinement (Uniform), based on interpolation error (Interpol). */

  "AdaptMechanics", "HangNode",
  /* Mechanics for adaptation.  Default is hanging-node (HangNode).
     The alternative is order refinement (OrderRef). */

  "AdaptFixedGrowth", "False",
  /* If True, a fixed growth dof strategy is used.  */

  "AdaptFixedFraction", "0.3",
  /* Fraction of elements adapted.  e.g. used for hanging-node
     adaptation. */

  "AdaptFixedGrowthFactor", "2.0",
  /* Desired growth factor of number of d.o.f. for each 
   adaptation iteration.  */

  "AdaptCoarsenFraction", "0.0",
  /* Fraction of d.o.f. coarsened in each adaptation iteration.  When
     used with the growth factor, the desired AdaptFixedGrowthFactor
     is preserved.  For example, AdaptFixedGrowthFactor = 1.0 together
     with AdaptCoarsenFraction > 0 would yield coarsening and
     refinement (d.o.f. redistribution) but no net d.o.f. growth.  */

  "AdaptIsotropic", "True",
  /* If True, only isotropic refinement will be considered.  For
     hanging node, this means only one uniform refinement option per
     element will be used. */
  
  "AdaptCostMetric", "DeltaDOF",
  /* Type of cost metric to be used in optimization-based 
   refinement direction decision. */
  
  "AdaptIncludeP", "False",
  /* If True, p-order is considered as a refinement option to 
   be compared with h-refinement options. */
  
  "AdaptPmax", "9",
  /* p-refinement is allowed up to p = AdaptPmax. */
  
  "AdaptWriteVOrder", "True",
  /* If True, a integer vector data will be written out. */

  "AdaptOutput", "None",
  /* For output-based adaptation: on which output do we want to
     adapt? */

  "AdaptScalar", "None",
  /* For scalar-based adaptation: on which scalar do we want to
     adapt? */

  "AdaptVariableSet", "None",
  /* Output-based adaptation can also proceed using a different
     variable set (e.g. entropy veriables) as adjoints.  This
     parameter defines the type of variable set (understood by the
     equation-set) for the adjoint variable definition.  Should
     only specify one of AdaptOutput or AdaptVariableSet. */

  "AdaptTolerance", "0.0",
  /* If applicable, adaptation stops when error measure drops below
     this tolerance. */
  
  "VariableWeights", "False",
  /* If True, the weights for each output considered for Error-based
     adaptation is calculated based on how well the error tolerance is
     met. */

  "AdaptSmoothRef", "False",
  /* If True, mesh refinement "smoothing" is carried out on the cells
     flagged for refinement -- i.e. those with a high error
     indicator. */

  "AdaptOrderMin", "0",
  /* Minimum order for p-adaptation */
  
  "AdaptOrderMax", "20",
  /* Maximum order for p-adaptation */

  /*=======================*/
  /* Robustness Adaptation */
  /*=======================*/
  
  "AdaptRobust", "False",
  /* If True, the solver is going to adapt the mesh based on residual
     distribution as an attempt to achieve convergence. */
  
  "AdaptRobustWriteInterm", "False",
  /* If True, intermediate meshes and states will be 
     written out. */
  
  "AdaptRobustIndicator", "Penalty",
  /* Indicator used for robustness adaptation.  Possibilities are same
     as the AdaptOn flag in regular Adaptation.  */
  
  "AdaptRobustIndBdown", "False",
  /* If True, the robustness adaptive indicator will be separated 
   into contribution of each equation.  */

  "AdaptRobustMaxIter", "2",
  /* Maximum number of adaptation iterations if AdaptRobustness is
     True. */

  "iAdaptRobust", "0",
  /* Current iteration number for robustness adaptation. */

  "AdaptRobustFixedFraction", "0.10",
  /* Fraction of elements refined when AdaptRobust is True. */

  "AdaptRobustMaxCFLAmpFactor", "1e1",
  /* Amplification factor for the maximum CFL achieved in the run. */

  "nCFLReducedMax", "2",
  /* Maximum number of CFL reductions. */


  /*=====================*/
  /* Unsteady Adaptation */
  /*=====================*/

  "DynamicSpatialRef", "False",
  /* If true, spatial refinement is applied dynamically during
     unsteady adaptive simulations. */

  "DynamicSpatialRefOrderInc", "0",
  /* Spatial order increment for dynamic p-refinement. This gets added
     to the orders in the VOrder files. */

  "VOrderFile", "None",
  /* File used to specify variable orders for the state vector, upon
     initializing the state. */

  "TimeHistorySmoothFactor", "0.0",
  /* Factor used to smooth adapted time histories.  0 means no
     smoothing, 1 means full Laplace smoothing (stencil of 0.25, 0.5,
     0.25).*/

  
  /*===============*/
  /* Stabilization */
  /*===============*/

  "StabSwitch", "Square",
  /* Type of stabilization switch function.  Used to convert a
     non-dimensional jump or regularity estimate into an indicator
     function.  Possible choices are: Const, Linear, Square,
     Log (given in xf_SolverStruct.h). */

  "StabSwitchValue", ".05",
  /* When StabSwitch == Const, this is the constant value used for the
     indicator. */

  "StabSwitchFactor", "1.0",
  /* The stabilization switch is multiplied by this factor. Useful if
     want to quickly try less/more stabilization. */

  "IsotropicHMetric", "False",
  /* If true, an isotropic element size (H) metric is constructed.
     Otherwise, an anisotropic metric is used. */

  "IsotropicHMetricIsHD", "True",
  /* If true, the hydraulic diameter (HD), proportional to
     Volume/SurfArea, is used for the element length scale in the
     isotropic element size metric.  Only relevant if IsotropicHMetric
     is also true.  The alternative is that the longest stretching
     length is used as the element diameter. */



  /*============*/
  /* Mesh Tools */
  /*============*/

  "DistFcnOrder", "1",
  /* Order for distance function interpolation. */

  "WriteDistFcn", "False",
  /* If True, distance function will be made writable if it is computed. */

  "DistFcnWallBoundaries", "NULL",
  /* Space-separated list of boundary titles to treat as walls.
     Overrides any equation-set designated walls.  Should not have to
     use this except in special cases of manufactured solutions where
     no true walls exist. */


  /*=============*/
  /* Mesh Motion */
  /*=============*/

  "MeshMotionFile", "None",
  /* Name of mesh motion file to use, if running with unsteady mesh
     motion.  This file describes grid deformation in time. */

  "MeshMotionActive", "True",
  /* Only relevant if MeshMotionFile is specified (or if there is a
     motion structure in the All->Mesh structure): this flag indicates
     whether the specified motion is active. */

  "UseGCL", "False",
  /* If True, a Geometric Conservation Law will be used. */

  /*===========*/
  /* Geometry  */
  /*===========*/

  "GeometryFile", "None",
  /* Name of geometry file to use.  This file describes geometry
     components that can be used to snap/curve meshes to the
     geometry. */


  /*=============*/
  /* Snapshots   */
  /*=============*/

  "nEqnSet", "1",
  /* Number of equation sets to run.  > 1 for taking snapshots.  The
     equation set file for each snapshot must take the form
     <eqnsetroot>_#.eqn, where # is the snapshot number (0
     .. nEqnSet-1).  <eqnsetroot> is the EqnSetFile parameter. */

  "SnapshotRestart", "False",
  /* True to restart from a previous solution if taking snapshots. */

  
  /*-----------------------------------------------------*/
  "\0"  /* DO NOT DELETE: null terminating string */
};



#endif // end ifndef _xf_ParamDefault_h
