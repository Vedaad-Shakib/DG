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
include prog/CMakeFiles/lydg_GeomConvert.dir/depend.make

# Include the progress variables for this target.
include prog/CMakeFiles/lydg_GeomConvert.dir/progress.make

# Include the compile flags for this target's objects.
include prog/CMakeFiles/lydg_GeomConvert.dir/flags.make

prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o: prog/CMakeFiles/lydg_GeomConvert.dir/flags.make
prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o: ../prog/xf_GeomConvert.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o   -c /home/vshakib/LYDG/prog/xf_GeomConvert.c

prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/prog/xf_GeomConvert.c > CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.i

prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/prog/xf_GeomConvert.c -o CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.s

prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o.requires:
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o.requires

prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o.provides: prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o.requires
	$(MAKE) -f prog/CMakeFiles/lydg_GeomConvert.dir/build.make prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o.provides.build
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o.provides

prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o.provides.build: prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o

prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o: prog/CMakeFiles/lydg_GeomConvert.dir/flags.make
prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o: ../prog/xf_Arg.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o   -c /home/vshakib/LYDG/prog/xf_Arg.c

prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/prog/xf_Arg.c > CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.i

prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/prog/xf_Arg.c -o CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.s

prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o.requires:
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o.requires

prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o.provides: prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o.requires
	$(MAKE) -f prog/CMakeFiles/lydg_GeomConvert.dir/build.make prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o.provides.build
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o.provides

prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o.provides.build: prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o

prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o: prog/CMakeFiles/lydg_GeomConvert.dir/flags.make
prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o: ../src/xf_EqnSetHook.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o   -c /home/vshakib/LYDG/src/xf_EqnSetHook.c

prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/src/xf_EqnSetHook.c > CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.i

prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/src/xf_EqnSetHook.c -o CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.s

prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o.requires:
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o.requires

prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o.provides: prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o.requires
	$(MAKE) -f prog/CMakeFiles/lydg_GeomConvert.dir/build.make prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o.provides.build
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o.provides

prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o.provides.build: prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o

# Object files for target lydg_GeomConvert
lydg_GeomConvert_OBJECTS = \
"CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o" \
"CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o" \
"CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o"

# External object files for target lydg_GeomConvert
lydg_GeomConvert_EXTERNAL_OBJECTS =

bin/lydg_GeomConvert: prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o
bin/lydg_GeomConvert: prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o
bin/lydg_GeomConvert: prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o
bin/lydg_GeomConvert: prog/CMakeFiles/lydg_GeomConvert.dir/build.make
bin/lydg_GeomConvert: lib/libxfSerial.a
bin/lydg_GeomConvert: /usr/lib64/libblas.so
bin/lydg_GeomConvert: /usr/lib64/liblapack.so
bin/lydg_GeomConvert: /usr/lib64/libblas.so
bin/lydg_GeomConvert: /usr/lib64/liblapack.so
bin/lydg_GeomConvert: prog/CMakeFiles/lydg_GeomConvert.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable ../bin/lydg_GeomConvert"
	cd /home/vshakib/LYDG/build/prog && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lydg_GeomConvert.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
prog/CMakeFiles/lydg_GeomConvert.dir/build: bin/lydg_GeomConvert
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/build

prog/CMakeFiles/lydg_GeomConvert.dir/requires: prog/CMakeFiles/lydg_GeomConvert.dir/xf_GeomConvert.c.o.requires
prog/CMakeFiles/lydg_GeomConvert.dir/requires: prog/CMakeFiles/lydg_GeomConvert.dir/xf_Arg.c.o.requires
prog/CMakeFiles/lydg_GeomConvert.dir/requires: prog/CMakeFiles/lydg_GeomConvert.dir/__/src/xf_EqnSetHook.c.o.requires
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/requires

prog/CMakeFiles/lydg_GeomConvert.dir/clean:
	cd /home/vshakib/LYDG/build/prog && $(CMAKE_COMMAND) -P CMakeFiles/lydg_GeomConvert.dir/cmake_clean.cmake
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/clean

prog/CMakeFiles/lydg_GeomConvert.dir/depend:
	cd /home/vshakib/LYDG/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vshakib/LYDG /home/vshakib/LYDG/prog /home/vshakib/LYDG/build /home/vshakib/LYDG/build/prog /home/vshakib/LYDG/build/prog/CMakeFiles/lydg_GeomConvert.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : prog/CMakeFiles/lydg_GeomConvert.dir/depend

