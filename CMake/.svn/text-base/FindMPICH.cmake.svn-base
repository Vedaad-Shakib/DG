set(MPI_INCLUDE_DIR
    /usr/include/mpi
    /usr/include/mpich2
    $ENV{MPI_HOME}
    $ENV{MPI_HOME}/include
)

set(MPI_LINK_DIR
    $ENV{MPI_HOME}
    $ENV{MPI_HOME}/lib
)

FIND_PATH(MPICH_INC mpi.h
  HINTS ${MPI_INC_PATH}
)

FIND_LIBRARY(MPICH_LIB
  NAMES libmpich.so libmpich.so.1 libmpich.so.2 libmpich.a
  HINTS ${MPI_LIB_PATH}
)

FIND_LIBRARY(OPA_LIB
  NAMES libopa.so libopa.so.1 libopa.la libopa.a
  HINTS ${MPI_LIB_PATH}
)

FIND_LIBRARY(MPL_LIB
  NAMES libmpl.so libmpl.so.1 libmpl.a
  HINTS ${MPI_LIB_PATH}
)

if (MPICH_LIB AND MPL_LIB)
  set(MPI_C_LIBRARIES ${MPICH_LIB} ${MPL_LIB})
  set(MPI_C_INCLUDE_PATH ${MPICH_INC})
  set(MPI_C_FOUND TRUE)
endif (MPICH_LIB AND MPL_LIB)

