# Eqnset hook source
set( XF_HOOK_SRC ${XFDIR}/xf_EqnSetHook.c )
#set(ODESOLVER_LIB "-L${TOPDIR}/chemkin/lib -lchemkin -lsundials_cvode -lsundials_nvecserial")
#set(FORTRAN_SUPPORT "-L/usr/local/gfortran/lib -lgfortran -lquadmath")

# Model-reduction source
set(MR_SRC
xf_MRCommon.c
)

# Initial-condition inversion source
set(IC_SRC
xf_InverseIC_Common.c
)

## Macro to build tools
macro(serial_prog PROGNAME SRC)
  set(LIBS ${ARGN})
  add_executable(${PROGNAME} ${SRC} ${XF_HOOK_SRC})
#  target_link_libraries(${PROGNAME} ${DYLINK_LIBRARY} xfSerial DTLCHEMISTRY ${LIBS}
  target_link_libraries(${PROGNAME} ${DYLINK_LIBRARY} xfSerial ${LIBS}
     ${ODESOLVER_LIB}
     ${FORTRAN_SUPPORT}
     ${MATH_LIBRARY})
  install(TARGETS ${PROGNAME} DESTINATION ${BINDIR})
endmacro(serial_prog)
macro(parallel_prog PROGNAME SRC)
  set(LIBS ${ARGN})
  add_executable(${PROGNAME} ${SRC} ${XF_HOOK_SRC})
#  target_link_libraries(${PROGNAME} ${DYLINK_LIBRARY} xfParallel DTLCHEMISTRY ${PARALLEL_LIBRARY} ${LIBS}
  target_link_libraries(${PROGNAME} ${DYLINK_LIBRARY} xfParallel ${PARALLEL_LIBRARY} ${LIBS}
     ${ODESOLVER_LIB}
     ${FORTRAN_SUPPORT}
     ${MATH_LIBRARY})
  install(TARGETS ${PROGNAME} DESTINATION ${BINDIR})
endmacro(parallel_prog)

# main executables
serial_prog  (lydg "xf_XFlow.c")
if (PARALLEL_BUILD)
  parallel_prog(plydg "xf_XFlow.c")
endif (PARALLEL_BUILD)

# get info from xfa file
serial_prog(lydg_Info "xf_Info.c")

#serial_prog(lydg_TecOut "xf_TecOut.c")
serial_prog(lydg_TecOut "xf_Tecout.c")

# post-processing of outputs
serial_prog(lydg_Post "xf_Post.c;xf_Arg.c")
if ( PARALLEL_BUILD )
  parallel_prog(plydg_Post "xf_Post.c;xf_Arg.c;${MR_SRC}")
endif ( PARALLEL_BUILD )

# mesh statistics utility (xf_MeshStats)
serial_prog(lydg_MeshStats "xf_MeshStats.c;xf_Arg.c")

if ( PARALLEL_BUILD )
  parallel_prog(plydg_MeshStats "xf_MeshStats.c;xf_Arg.c")
endif ( PARALLEL_BUILD )


# mesh conversion utility (xf_Convert)
serial_prog(lydg_Convert "xf_Convert.c;xf_Arg.c")

# geometry conversion utility (xf_GeomConvert)
serial_prog(lydg_GeomConvert "xf_GeomConvert.c;xf_Arg.c")

# for converting a scalar potential to a velocity field (xf_Pot2Vel)
serial_prog(lydg_Pot2Vel "xf_Pot2Vel.c;xf_Arg.c")

# openGL plotter
if (GRAPHICS_BUILD)
  serial_prog(lydg_Plot "xf_Plot.c;xf_Arg.c" "${GRAPHICS_LIBRARY}")
endif (GRAPHICS_BUILD)

# openGL plotter
if (GRAPHICS_BUILD)
   serial_prog(lydg_FigureOut "xf_FigureOut.c;xf_Arg.c" "${GRAPHICS_LIBRARY}")
endif (GRAPHICS_BUILD)

# data comparison program (xf_DataCompare)
serial_prog(lydg_DataCompare "xf_DataCompare.c;xf_Arg.c")

# data peeling program (xf_DataPeel)
serial_prog(lydg_DataPeel "xf_DataPeel.c;xf_Arg.c")

# mesh refinement (xf_Refine)
serial_prog(lydg_Refine "xf_Refine.c;xf_Arg.c")

# high order mesh generation from structured mesh (xf_StructuredHO)
serial_prog(lydg_StructuredHO "xf_StructuredHO.c;xf_Arg.c")

# interactive all/mesh editing (xf_Edit)
serial_prog(lydg_Edit "xf_Edit.c;xf_Arg.c")

# data operation program (xf_DataOper)
serial_prog(lydg_DataOper "xf_DataOper.c;xf_Arg.c")

if ( PARALLEL_BUILD )
  parallel_prog(plydg_DataOper "xf_DataOper.c;xf_Arg.c;${MR_SRC}")
endif ( PARALLEL_BUILD )

# state alteration (xf_AlterState)
serial_prog(lydg_AlterState "xf_AlterState.c;xf_Arg.c")

# Initial condition sensitivity analysis
serial_prog(lydg_ICSensitivity "xf_ICSensitivity.c;xf_Arg.c")

# data to text conversion (xf_Data2Text)
serial_prog(lydg_Data2Text "xf_Data2Text.c;xf_Arg.c")

# operator writing
serial_prog(lydg_WriteOperators "xf_WriteOperators.c;xf_Arg.c")

# parameter continuation for steady-state solution
serial_prog(lydg_SteadyContinue "xf_SteadyContinue.c;xf_Arg.c")
if ( PARALLEL_BUILD )
  parallel_prog(plydg_SteadyContinue "xf_SteadyContinue.c;xf_Arg.c;${MR_SRC}")
endif ( PARALLEL_BUILD )

#--------------------------#
# Model Reduction programs #
#--------------------------#

# Offline (xf_MROffline)
serial_prog(lydg_MROffline "xf_MROffline;xf_Arg.c;${MR_SRC}")

if ( PARALLEL_BUILD )
  parallel_prog(plydg_MROffline "xf_MROffline.c;xf_Arg.c;${MR_SRC}")
endif ( PARALLEL_BUILD )

# Online (xf_MROnline)
serial_prog(lydg_MROnline "xf_MROnline.c;xf_Arg.c;${MR_SRC}")

# Solution reconstruction utility (xf_MRReconstruct)
serial_prog(lydg_MRReconstruct "xf_MRReconstruct.c;xf_Arg.c;${MR_SRC}")

# Tester
serial_prog(lydg_MRTest "xf_MRTest.c;xf_Arg.c;${MR_SRC}")

# Hessian-based Initial Condition calculation
if ( LAPACK_FOUND )
  serial_prog(lydg_MRHessianIC "xf_MRHessianIC.c;xf_Arg.c;${MR_SRC}")
endif ( LAPACK_FOUND )
if ( LAPACK_FOUND AND PARALLEL_BUILD )
  parallel_prog(plydg_MRHessianIC "xf_MRHessianIC.c;xf_Arg.c;${MR_SRC}")
endif ( LAPACK_FOUND AND PARALLEL_BUILD )

#-----------------------------#
# Initial Condition Inversion #
#-----------------------------#


# Inverse Initial Condition calculation
serial_prog(lydg_InverseIC "xf_InverseIC.c;xf_Arg.c;${IC_SRC}")

if ( PARALLEL_BUILD )
  parallel_prog(plydg_InverseIC "xf_InverseIC.c;xf_Arg.c;${IC_SRC}")
endif ( PARALLEL_BUILD )

# MCMC  Inverse Initial Condition calculation
serial_prog(lydg_InverseIC_MCMC "xf_InverseIC_MCMC.c;xf_Arg.c;${IC_SRC}")

if ( PARALLEL_BUILD )
  parallel_prog(plydg_InverseIC_MCMC "xf_InverseIC_MCMC.c;xf_Arg.c;${IC_SRC}")
endif ( PARALLEL_BUILD )


# Ensemble MCMC  Inverse Initial Condition calculation
serial_prog(lydg_InverseIC_MCMCEnsemble "xf_InverseIC_MCMCEnsemble.c;xf_Arg.c;${IC_SRC}")

if ( PARALLEL_BUILD )
  parallel_prog(plydg_InverseIC_MCMCEnsemble "xf_InverseIC_MCMCEnsemble.c;xf_Arg.c;${IC_SRC}")
endif ( PARALLEL_BUILD )

# Eigenvalue solver
if ( LAPACK_FOUND )
  serial_prog(lydg_EigSolve "xf_EigSolve.c;xf_Arg.c")
endif ( LAPACK_FOUND )
if ( LAPACK_FOUND AND PARALLEL_BUILD )
  parallel_prog(plydg_EigSolve "xf_EigSolve.c;xf_Arg.c")
endif ( LAPACK_FOUND AND PARALLEL_BUILD )
