#
# Find the PARMETIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for 
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# PARMETIS_INCLUDE_DIR - where to find autopack.h
# PARMETIS_LIBRARIES   - List of fully qualified libraries to link against.
# PARMETIS_FOUND       - Do not attempt to use if false or undefined.

FIND_PATH(PARMETIS_INCLUDE_PATH parmetis.h
HINTS ${PARMETIS_INCLUDE_DIR} $ENV{PARMETIS_INC}
)

FIND_LIBRARY(PARMETIS_LIBRARY parmetis
HINTS ${PARMETIS_LINK_DIR} $ENV{PARMETIS_LINK}
)

FIND_LIBRARY(METIS_LIBRARY metis
HINTS ${PARMETIS_LINK_DIR} $ENV{PARMETIS_LINK} $ENV{METIS_LINK}
)

IF(PARMETIS_INCLUDE_PATH AND PARMETIS_LIBRARY AND METIS_LIBRARY)
    SET( PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
    SET( PARMETIS_FOUND TRUE )
ENDIF(PARMETIS_INCLUDE_PATH AND PARMETIS_LIBRARY AND METIS_LIBRARY)

