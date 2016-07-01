# Apple OS X computers require this link flag
if ( APPLE )
  set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
  set(CMAKE_SHARED_LINKER_FLAGS "-undefined dynamic_lookup")
endif ( APPLE )

# GCC Flags
if ( CMAKE_C_COMPILER MATCHES "gcc" AND NOT CMAKE_C_COMPILER MATCHES "pgcc")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wimplicit-int -Werror-implicit-function-declaration -fno-stack-protector -mmacosx-version-min=10.6")
endif ( CMAKE_C_COMPILER MATCHES "gcc" AND NOT CMAKE_C_COMPILER MATCHES "pgcc")

# PGCC Settings
if ( CMAKE_C_COMPILER MATCHES "pgcc" )
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Minline -Mnovect")
endif ( CMAKE_C_COMPILER MATCHES "pgcc" )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -pg")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")




