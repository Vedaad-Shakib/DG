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
include prog/CMakeFiles/plydg_Post.dir/depend.make

# Include the progress variables for this target.
include prog/CMakeFiles/plydg_Post.dir/progress.make

# Include the compile flags for this target's objects.
include prog/CMakeFiles/plydg_Post.dir/flags.make

prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o: prog/CMakeFiles/plydg_Post.dir/flags.make
prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o: ../prog/xf_Post.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/plydg_Post.dir/xf_Post.c.o   -c /home/vshakib/LYDG/prog/xf_Post.c

prog/CMakeFiles/plydg_Post.dir/xf_Post.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/plydg_Post.dir/xf_Post.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/prog/xf_Post.c > CMakeFiles/plydg_Post.dir/xf_Post.c.i

prog/CMakeFiles/plydg_Post.dir/xf_Post.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/plydg_Post.dir/xf_Post.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/prog/xf_Post.c -o CMakeFiles/plydg_Post.dir/xf_Post.c.s

prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o.requires:
.PHONY : prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o.requires

prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o.provides: prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o.requires
	$(MAKE) -f prog/CMakeFiles/plydg_Post.dir/build.make prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o.provides.build
.PHONY : prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o.provides

prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o.provides.build: prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o

prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o: prog/CMakeFiles/plydg_Post.dir/flags.make
prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o: ../prog/xf_Arg.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/plydg_Post.dir/xf_Arg.c.o   -c /home/vshakib/LYDG/prog/xf_Arg.c

prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/plydg_Post.dir/xf_Arg.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/prog/xf_Arg.c > CMakeFiles/plydg_Post.dir/xf_Arg.c.i

prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/plydg_Post.dir/xf_Arg.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/prog/xf_Arg.c -o CMakeFiles/plydg_Post.dir/xf_Arg.c.s

prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o.requires:
.PHONY : prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o.requires

prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o.provides: prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o.requires
	$(MAKE) -f prog/CMakeFiles/plydg_Post.dir/build.make prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o.provides.build
.PHONY : prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o.provides

prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o.provides.build: prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o

prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o: prog/CMakeFiles/plydg_Post.dir/flags.make
prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o: ../prog/xf_MRCommon.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o   -c /home/vshakib/LYDG/prog/xf_MRCommon.c

prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/plydg_Post.dir/xf_MRCommon.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/prog/xf_MRCommon.c > CMakeFiles/plydg_Post.dir/xf_MRCommon.c.i

prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/plydg_Post.dir/xf_MRCommon.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/prog/xf_MRCommon.c -o CMakeFiles/plydg_Post.dir/xf_MRCommon.c.s

prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o.requires:
.PHONY : prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o.requires

prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o.provides: prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o.requires
	$(MAKE) -f prog/CMakeFiles/plydg_Post.dir/build.make prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o.provides.build
.PHONY : prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o.provides

prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o.provides.build: prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o

prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o: prog/CMakeFiles/plydg_Post.dir/flags.make
prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o: ../src/xf_EqnSetHook.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o   -c /home/vshakib/LYDG/src/xf_EqnSetHook.c

prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/src/xf_EqnSetHook.c > CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.i

prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/src/xf_EqnSetHook.c -o CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.s

prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o.requires:
.PHONY : prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o.requires

prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o.provides: prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o.requires
	$(MAKE) -f prog/CMakeFiles/plydg_Post.dir/build.make prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o.provides.build
.PHONY : prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o.provides

prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o.provides.build: prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o

# Object files for target plydg_Post
plydg_Post_OBJECTS = \
"CMakeFiles/plydg_Post.dir/xf_Post.c.o" \
"CMakeFiles/plydg_Post.dir/xf_Arg.c.o" \
"CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o" \
"CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o"

# External object files for target plydg_Post
plydg_Post_EXTERNAL_OBJECTS =

bin/plydg_Post: prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o
bin/plydg_Post: prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o
bin/plydg_Post: prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o
bin/plydg_Post: prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o
bin/plydg_Post: prog/CMakeFiles/plydg_Post.dir/build.make
bin/plydg_Post: lib/libxfParallel.a
bin/plydg_Post: /share/apps/mvapich2/2.1rc1-intel-15/lib/libmpi.so
bin/plydg_Post: /share/apps/lib/parmetis/4.0.3/mvapich2-2.1rc1-intel-15/lib/libparmetis.a
bin/plydg_Post: /share/apps/lib/parmetis/4.0.3/mvapich2-2.1rc1-intel-15/lib/libmetis.a
bin/plydg_Post: /usr/lib64/libblas.so
bin/plydg_Post: /usr/lib64/liblapack.so
bin/plydg_Post: /usr/lib64/libblas.so
bin/plydg_Post: /usr/lib64/liblapack.so
bin/plydg_Post: prog/CMakeFiles/plydg_Post.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable ../bin/plydg_Post"
	cd /home/vshakib/LYDG/build/prog && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/plydg_Post.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
prog/CMakeFiles/plydg_Post.dir/build: bin/plydg_Post
.PHONY : prog/CMakeFiles/plydg_Post.dir/build

prog/CMakeFiles/plydg_Post.dir/requires: prog/CMakeFiles/plydg_Post.dir/xf_Post.c.o.requires
prog/CMakeFiles/plydg_Post.dir/requires: prog/CMakeFiles/plydg_Post.dir/xf_Arg.c.o.requires
prog/CMakeFiles/plydg_Post.dir/requires: prog/CMakeFiles/plydg_Post.dir/xf_MRCommon.c.o.requires
prog/CMakeFiles/plydg_Post.dir/requires: prog/CMakeFiles/plydg_Post.dir/__/src/xf_EqnSetHook.c.o.requires
.PHONY : prog/CMakeFiles/plydg_Post.dir/requires

prog/CMakeFiles/plydg_Post.dir/clean:
	cd /home/vshakib/LYDG/build/prog && $(CMAKE_COMMAND) -P CMakeFiles/plydg_Post.dir/cmake_clean.cmake
.PHONY : prog/CMakeFiles/plydg_Post.dir/clean

prog/CMakeFiles/plydg_Post.dir/depend:
	cd /home/vshakib/LYDG/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vshakib/LYDG /home/vshakib/LYDG/prog /home/vshakib/LYDG/build /home/vshakib/LYDG/build/prog /home/vshakib/LYDG/build/prog/CMakeFiles/plydg_Post.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : prog/CMakeFiles/plydg_Post.dir/depend

