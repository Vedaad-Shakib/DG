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
include prog/CMakeFiles/lydg_ICSensitivity.dir/depend.make

# Include the progress variables for this target.
include prog/CMakeFiles/lydg_ICSensitivity.dir/progress.make

# Include the compile flags for this target's objects.
include prog/CMakeFiles/lydg_ICSensitivity.dir/flags.make

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o: prog/CMakeFiles/lydg_ICSensitivity.dir/flags.make
prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o: ../prog/xf_ICSensitivity.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o   -c /home/vshakib/LYDG/prog/xf_ICSensitivity.c

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/prog/xf_ICSensitivity.c > CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.i

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/prog/xf_ICSensitivity.c -o CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.s

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o.requires:
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o.requires

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o.provides: prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o.requires
	$(MAKE) -f prog/CMakeFiles/lydg_ICSensitivity.dir/build.make prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o.provides.build
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o.provides

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o.provides.build: prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o: prog/CMakeFiles/lydg_ICSensitivity.dir/flags.make
prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o: ../prog/xf_Arg.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o   -c /home/vshakib/LYDG/prog/xf_Arg.c

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/prog/xf_Arg.c > CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.i

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/prog/xf_Arg.c -o CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.s

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o.requires:
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o.requires

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o.provides: prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o.requires
	$(MAKE) -f prog/CMakeFiles/lydg_ICSensitivity.dir/build.make prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o.provides.build
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o.provides

prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o.provides.build: prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o

prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o: prog/CMakeFiles/lydg_ICSensitivity.dir/flags.make
prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o: ../src/xf_EqnSetHook.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vshakib/LYDG/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o   -c /home/vshakib/LYDG/src/xf_EqnSetHook.c

prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.i"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/vshakib/LYDG/src/xf_EqnSetHook.c > CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.i

prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.s"
	cd /home/vshakib/LYDG/build/prog && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/vshakib/LYDG/src/xf_EqnSetHook.c -o CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.s

prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o.requires:
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o.requires

prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o.provides: prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o.requires
	$(MAKE) -f prog/CMakeFiles/lydg_ICSensitivity.dir/build.make prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o.provides.build
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o.provides

prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o.provides.build: prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o

# Object files for target lydg_ICSensitivity
lydg_ICSensitivity_OBJECTS = \
"CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o" \
"CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o" \
"CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o"

# External object files for target lydg_ICSensitivity
lydg_ICSensitivity_EXTERNAL_OBJECTS =

bin/lydg_ICSensitivity: prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o
bin/lydg_ICSensitivity: prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o
bin/lydg_ICSensitivity: prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o
bin/lydg_ICSensitivity: prog/CMakeFiles/lydg_ICSensitivity.dir/build.make
bin/lydg_ICSensitivity: lib/libxfSerial.a
bin/lydg_ICSensitivity: /usr/lib64/libblas.so
bin/lydg_ICSensitivity: /usr/lib64/liblapack.so
bin/lydg_ICSensitivity: /usr/lib64/libblas.so
bin/lydg_ICSensitivity: /usr/lib64/liblapack.so
bin/lydg_ICSensitivity: prog/CMakeFiles/lydg_ICSensitivity.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable ../bin/lydg_ICSensitivity"
	cd /home/vshakib/LYDG/build/prog && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lydg_ICSensitivity.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
prog/CMakeFiles/lydg_ICSensitivity.dir/build: bin/lydg_ICSensitivity
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/build

prog/CMakeFiles/lydg_ICSensitivity.dir/requires: prog/CMakeFiles/lydg_ICSensitivity.dir/xf_ICSensitivity.c.o.requires
prog/CMakeFiles/lydg_ICSensitivity.dir/requires: prog/CMakeFiles/lydg_ICSensitivity.dir/xf_Arg.c.o.requires
prog/CMakeFiles/lydg_ICSensitivity.dir/requires: prog/CMakeFiles/lydg_ICSensitivity.dir/__/src/xf_EqnSetHook.c.o.requires
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/requires

prog/CMakeFiles/lydg_ICSensitivity.dir/clean:
	cd /home/vshakib/LYDG/build/prog && $(CMAKE_COMMAND) -P CMakeFiles/lydg_ICSensitivity.dir/cmake_clean.cmake
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/clean

prog/CMakeFiles/lydg_ICSensitivity.dir/depend:
	cd /home/vshakib/LYDG/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vshakib/LYDG /home/vshakib/LYDG/prog /home/vshakib/LYDG/build /home/vshakib/LYDG/build/prog /home/vshakib/LYDG/build/prog/CMakeFiles/lydg_ICSensitivity.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : prog/CMakeFiles/lydg_ICSensitivity.dir/depend

