# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /share/apps/cmake/2.8.12.2/bin/cmake

# The command to remove a file.
RM = /share/apps/cmake/2.8.12.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /share/apps/cmake/2.8.12.2/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vshakib/LYDG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vshakib/LYDG/build

# Include any dependencies generated for this target.
include prog/CMakeFiles/plydg.dir/depend.make

# Include the progress variables for this target.
include prog/CMakeFiles/plydg.dir/progress.make

# Include the compile flags for this target's objects.
include prog/CMakeFiles/plydg.dir/flags.make

prog/CMakeFiles/plydg.dir/xf_XFlow.c.o: prog/CMakeFiles/plydg.dir/flags.make
prog/CMakeFiles/plydg.dir/xf_XFlow.c.o: ../prog/xf_XFlow.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/plydg.dir/xf_XFlow.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/plydg.dir/xf_XFlow.c.o   -c /home/vshakib/LYDG/prog/xf_XFlow.c

prog/CMakeFiles/plydg.dir/xf_XFlow.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/plydg.dir/xf_XFlow.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/prog/xf_XFlow.c > CMakeFiles/plydg.dir/xf_XFlow.c.i

prog/CMakeFiles/plydg.dir/xf_XFlow.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/plydg.dir/xf_XFlow.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/prog/xf_XFlow.c -o CMakeFiles/plydg.dir/xf_XFlow.c.s

prog/CMakeFiles/plydg.dir/xf_XFlow.c.o.requires:
.PHONY : prog/CMakeFiles/plydg.dir/xf_XFlow.c.o.requires

prog/CMakeFiles/plydg.dir/xf_XFlow.c.o.provides: prog/CMakeFiles/plydg.dir/xf_XFlow.c.o.requires
	$(MAKE) -f prog/CMakeFiles/plydg.dir/build.make prog/CMakeFiles/plydg.dir/xf_XFlow.c.o.provides.build
.PHONY : prog/CMakeFiles/plydg.dir/xf_XFlow.c.o.provides

prog/CMakeFiles/plydg.dir/xf_XFlow.c.o.provides.build: prog/CMakeFiles/plydg.dir/xf_XFlow.c.o

prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o: prog/CMakeFiles/plydg.dir/flags.make
prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o: ../src/xf_EqnSetHook.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o   -c /home/vshakib/LYDG/src/xf_EqnSetHook.c

prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/src/xf_EqnSetHook.c > CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.i

prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/src/xf_EqnSetHook.c -o CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.s

prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o.requires:
.PHONY : prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o.requires

prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o.provides: prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o.requires
	$(MAKE) -f prog/CMakeFiles/plydg.dir/build.make prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o.provides.build
.PHONY : prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o.provides

prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o.provides.build: prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o

# Object files for target plydg
plydg_OBJECTS = \
"CMakeFiles/plydg.dir/xf_XFlow.c.o" \
"CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o"

# External object files for target plydg
plydg_EXTERNAL_OBJECTS =

bin/plydg: prog/CMakeFiles/plydg.dir/xf_XFlow.c.o
bin/plydg: prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o
bin/plydg: prog/CMakeFiles/plydg.dir/build.make
bin/plydg: lib/libxfParallel.a
bin/plydg: /share/apps/mvapich2/2.1rc1-intel-15/lib/libmpi.so
bin/plydg: /share/apps/lib/parmetis/4.0.3/mvapich2-2.1rc1-intel-15/lib/libparmetis.a
bin/plydg: /share/apps/lib/parmetis/4.0.3/mvapich2-2.1rc1-intel-15/lib/libmetis.a
bin/plydg: /usr/lib64/libblas.so
bin/plydg: /usr/lib64/liblapack.so
bin/plydg: /usr/lib64/libblas.so
bin/plydg: /usr/lib64/liblapack.so
bin/plydg: prog/CMakeFiles/plydg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable ../bin/plydg"
	cd /home/vshakib/LYDG/build/prog && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/plydg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
prog/CMakeFiles/plydg.dir/build: bin/plydg
.PHONY : prog/CMakeFiles/plydg.dir/build

prog/CMakeFiles/plydg.dir/requires: prog/CMakeFiles/plydg.dir/xf_XFlow.c.o.requires
prog/CMakeFiles/plydg.dir/requires: prog/CMakeFiles/plydg.dir/__/src/xf_EqnSetHook.c.o.requires
.PHONY : prog/CMakeFiles/plydg.dir/requires

prog/CMakeFiles/plydg.dir/clean:
	cd /home/vshakib/LYDG/build/prog && $(CMAKE_COMMAND) -P CMakeFiles/plydg.dir/cmake_clean.cmake
.PHONY : prog/CMakeFiles/plydg.dir/clean

prog/CMakeFiles/plydg.dir/depend:
	cd /home/vshakib/LYDG/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vshakib/LYDG /home/vshakib/LYDG/prog /home/vshakib/LYDG/build /home/vshakib/LYDG/build/prog /home/vshakib/LYDG/build/prog/CMakeFiles/plydg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : prog/CMakeFiles/plydg.dir/depend
