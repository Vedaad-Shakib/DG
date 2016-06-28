################################################
### XFlow CMake Build Interface
################################################

### Set up the environment
cmake_minimum_required(VERSION 2.8)
project(XFlow C)
enable_language(CXX OPTIONAL)
enable_language(Fortran OPTIONAL)

## Set the module directory
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

include(XFlowDirs)
include(XFlowDeps)
include(CompilerSettings)

# Build type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Build types are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Include the header files
include_directories(${INCDIR})
include_directories(${TOPDIR}/chemkin/ODEsolver_forXFLOW)

## Where to put executables and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# RPATH link flags, this puts the library path directly into the binary, so
# LD_LIBRARY_PATH does not need to be set
# Note: If you want to use another version of the library, use the RUNPATH
# environment variable to override RPATH
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


#---------------#
# Base source 
#---------------#
# Include the header files for Xcode and other development tools
file(GLOB XF_INC ${XFDIR}/*.h)
file(GLOB CHEM_HEADER ${XFDIR}/chemkin/ODEsolver_forXFLOW/*.h)
set( XF_BASE_SRC
${XFDIR}/xf_All.c
${XFDIR}/xfYu_Model.c
${XFDIR}/xf_Mesh.c
${XFDIR}/xf_Heap.c
${XFDIR}/xf_Geom.c
${XFDIR}/xf_GeomIO.c
${XFDIR}/xf_GeomAnalytical.c
${XFDIR}/xf_GeomSpline.c
${XFDIR}/xf_GeomPanel.c
${XFDIR}/xf_Data.c
${XFDIR}/xf_Param.c
${XFDIR}/xf_LinQueue.c
${XFDIR}/xf_Memory.c
${XFDIR}/xf_String.c
${XFDIR}/xf_Math.c
${XFDIR}/xf_DataMath.c
${XFDIR}/xf_Error.c
#${XFDIR}/xf_Basis.c
#${XFDIR}/xf_BasisFcn.c
${XFDIR}/xfYu_Basis.c
${XFDIR}/xfYu_BasisFcn.c
${XFDIR}/xf_MeshTools.c
${XFDIR}/xf_MeshDistance.c
${XFDIR}/xf_MeshMotionGCL.c
${XFDIR}/xf_MeshMotion.c
${XFDIR}/xf_MeshMotionIO.c
${XFDIR}/xf_MeshMotionAnalytical.c
${XFDIR}/xf_EqnSet.c
${XFDIR}/xf_IO.c
${XFDIR}/xf_Dynamic.c
${XFDIR}/xf_Log.c
${XFDIR}/xf_Solver.c
${XFDIR}/xf_SolverTools.c
${XFDIR}/xf_Line.c
${XFDIR}/xf_LinearSolver.c
${XFDIR}/xf_EigSolver.c
${XFDIR}/xf_Residual.c
${XFDIR}/xf_ResidualDiff.c
${XFDIR}/xf_ResidualStab.c
${XFDIR}/xf_Quad.c
${XFDIR}/xf_QuadRule.c
${XFDIR}/xf_Output.c
${XFDIR}/xf_Adapt.c
${XFDIR}/xf_AdaptHang.c
${XFDIR}/xf_ErrEst.c
${XFDIR}/xf_LeanSolver.c
${XFDIR}/xf_AllPull.c
${XFDIR}/xf_Penalty.c
#${XFDIR}/xf_DoubleFlux.c
${XFDIR}/xfYu_Solver.c
${XFDIR}/xfYu_AdaptSolver.c
#${XFDIR}/xfYu_Residual.c
${XFDIR}/Yu_DiffusDiscret.c
${XFDIR}/xfYu_Limiter.c
${XFDIR}/xfYu_EntropyBounding.c
${XFDIR}/xfYu_Statistics.c
)

set(XF_HOOK_SRC
${XFDIR}/xf_EqnSetHook.c
)

##just added on Feb 20 for detailed chemistry
set(XF_CHEMKIN_SRC
#${TOPDIR}/chemkin/ODEsolver_forXFLOW/reaction_routine.c
#${TOPDIR}/chemkin/ODEsolver_forXFLOW/homo0D.c
)

## Add math to base source
if ( BLAS_FOUND )
  set(XF_BASE_SRC ${XFDIR}/xf_MathBlas.c ${XF_BASE_SRC})
else ( BLAS_FOUND )
  set(XF_BASE_SRC ${XFDIR}/xf_MathNoBlas.c ${XF_BASE_SRC})
endif ( BLAS_FOUND )

if ( LAPACK_FOUND )
  set(XF_BASE_SRC ${XFDIR}/xf_MathLapack.c ${XF_BASE_SRC})
else ( LAPACK_FOUND )
  set(XF_BASE ${XFDIR}/xf_MathNoLapack.c ${XF_BASE_SRC})
endif ( LAPACK_FOUND )

## Serial and parallel source addition
set( XF_SERIAL_SRC ${XFDIR}/xf_Serial.c )
set( XF_PARALLEL_SRC
${XFDIR}/xf_MPI.c
${XFDIR}/xf_Partition.c
)

#---------------#
# XFlow library
#---------------#
###################################################
#add_library(DTLCHEMISTRY
#${XF_CHEMKIN_SRC}
#${CHEM_HEADER}
#)
#install(TARGETS DTLCHEMISTRY DESTINATION ${LIBDIR})
###################################################   
add_library(xfSerial ${XF_SERIAL_SRC} ${XF_BASE_SRC} ${XF_INC} ${CHEM_HEADER})
install(TARGETS xfSerial DESTINATION ${LIBDIR})
if ( PARALLEL_BUILD )
   add_library(xfParallel ${XF_PARALLEL_SRC} ${XF_BASE_SRC} ${XF_INC} ${CHEM_HEADER})  
    install(TARGETS xfParallel DESTINATION ${LIBDIR})
endif ( PARALLEL_BUILD )

#---------------#
# Unit testing
#---------------#
#enable_testing()
#include (UnitTest)
# Create 'make check'
#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
#add_custom_target(check_HydroCommon COMMAND ${CMAKE_CTEST_COMMAND} -R HydroCommon)
#add_custom_target(check_CompressibleNS COMMAND ${CMAKE_CTEST_COMMAND} -R CompressibleNS)
#add_custom_target(check_RadHydro COMMAND ${CMAKE_CTEST_COMMAND} -R RadHydro)
#add_custom_target(check_ShallowWater COMMAND ${CMAKE_CTEST_COMMAND} -R ShallowWater)

# For each unit test a .test.in file should exist
#set(UTESTS
#xf_Solver
#xf_Basis
#xf_MeshTools
#xf_AdaptHang
#xf_Math
#xf_MeshDistance
#xf_MeshMotionAnalytical
#xf_Residual
#xf_ResidualStab
#xf_QuadRule
#xf_DataMath
#xf_String
#xf_BasisFcn
#xf_LinearSolver
#xf_Line
#xf_ErrEst
#xf_SolverTools
#xf_Output
#)

#if ( BLAS_FOUND )
#  set(UTESTS ${UTESTS} xf_MathBlas)
#else ( BLAS_FOUND )
#  set(UTESTS ${UTESTS} xf_MathNoBlas)
#endif ( BLAS_FOUND )

#if ( LAPACK_FOUND )
#  set(UTESTS ${UTESTS} xf_MathLapack xf_EigSolver)
#endif ( LAPACK_FOUND )

#set(SRC_UTEST_ADD
#${XF_HOOK_SRC}
#)
#set(LIB_UTEST_ADD
#)

# Core unit tests
#foreach(TEST IN LISTS UTESTS)
#  unit_test(${XFDIR} ${TEST})
#  add_dependencies(check ${TEST}.utest)
#endforeach(TEST IN LISTS UTESTS)

# Dynamic libraries
#add_subdirectory(${DYNDIR})
# Programs
add_subdirectory(${PROGDIR})
# Regression tests
#add_subdirectory(${TESTDIR})

# Print configuration status
message (STATUS "======================================================================")
message (STATUS " XFlow suite configuration settings")
message (STATUS "======================================================================")
message (STATUS "")
message (STATUS " Compiler: ${CMAKE_C_COMPILER}")
message (STATUS " CFLAGS: ${CMAKE_C_FLAGS}")
message (STATUS "")
if (PARALLEL_BUILD)
    message (STATUS " Parallel build enabled:")
    message (STATUS "  * MPI Library: ${MPI_C_LIBRARIES}")
    message (STATUS "  * Parmetis Library: ${PARMETIS_LIBRARIES}")
else (PARALLEL_BUILD)
    message (STATUS " Parallel build disabled")
endif (PARALLEL_BUILD)
message (STATUS "")
if (LAPACK_FOUND)
    message (STATUS " LAPACK enabled:")
    message (STATUS "  * Lapack Library: ${LAPACK_LIBRARIES}")
else (LAPACK_FOUND)
    message (STATUS " LAPACK disabled")
endif (LAPACK_FOUND)
message (STATUS "")
if (BLAS_FOUND)
    message (STATUS " BLAS enabled:")
    message (STATUS "  * BLAS Library: ${BLAS_LIBRARIES}")
elseif (BLAS_FOUND)
    message (STATUS " BLAS disabled")
endif (BLAS_FOUND)
message (STATUS "")
message (STATUS "======================================================================")
