#-----------------------#
# XFlow unit test macro
#-----------------------#

macro (unit_test SRCDIR TEST)

# Set the unit test directory for objects
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/UnitTests)
set(UTESTDIR ${CMAKE_BINARY_DIR}/UnitTests)

# Copy the file over to the object directory
configure_file( ${SRCDIR}/${TEST}.test.c ${UTESTDIR}/${TEST}.test.c COPYONLY)
execute_process( COMMAND ${TESTDIR}/xf_test2utest.sh ${UTESTDIR}/${TEST})
# Add the CFLAG for testing
configure_file( ${SRCDIR}/${TEST}.c ${UTESTDIR}/${TEST}.c COPYONLY)
set_source_files_properties(${UTESTDIR}/${TEST}.c PROPERTIES COMPILE_FLAGS "-DUNIT_TEST=1")

# xf_UnitRun.c is in ${XFDIR}, so we need to add this directory
include_directories(${XFDIR})

# Build the unit test executable and link
add_executable(${TEST}.utest EXCLUDE_FROM_ALL ${UTESTDIR}/${TEST}.c ${XFDIR}/xf_Unit.c ${SRC_UTEST_ADD})
target_link_libraries(${TEST}.utest ${DYLINK_LIBRARY} xfSerial ${MATH_LIBRARY} ${LIB_UTEST_ADD})

# Add a CTest for the unit test
add_test(NAME ${TEST} COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TEST}.utest)

endmacro (unit_test)
