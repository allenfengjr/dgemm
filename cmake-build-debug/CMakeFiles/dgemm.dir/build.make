# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/haofeng-admin/dgemm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/haofeng-admin/dgemm/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/dgemm.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/dgemm.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/dgemm.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dgemm.dir/flags.make

CMakeFiles/dgemm.dir/dgemm.cpp.o: CMakeFiles/dgemm.dir/flags.make
CMakeFiles/dgemm.dir/dgemm.cpp.o: /Users/haofeng-admin/dgemm/dgemm.cpp
CMakeFiles/dgemm.dir/dgemm.cpp.o: CMakeFiles/dgemm.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/haofeng-admin/dgemm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dgemm.dir/dgemm.cpp.o"
	/opt/homebrew/bin/mpicc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dgemm.dir/dgemm.cpp.o -MF CMakeFiles/dgemm.dir/dgemm.cpp.o.d -o CMakeFiles/dgemm.dir/dgemm.cpp.o -c /Users/haofeng-admin/dgemm/dgemm.cpp

CMakeFiles/dgemm.dir/dgemm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dgemm.dir/dgemm.cpp.i"
	/opt/homebrew/bin/mpicc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/haofeng-admin/dgemm/dgemm.cpp > CMakeFiles/dgemm.dir/dgemm.cpp.i

CMakeFiles/dgemm.dir/dgemm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dgemm.dir/dgemm.cpp.s"
	/opt/homebrew/bin/mpicc $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/haofeng-admin/dgemm/dgemm.cpp -o CMakeFiles/dgemm.dir/dgemm.cpp.s

# Object files for target dgemm
dgemm_OBJECTS = \
"CMakeFiles/dgemm.dir/dgemm.cpp.o"

# External object files for target dgemm
dgemm_EXTERNAL_OBJECTS =

dgemm: CMakeFiles/dgemm.dir/dgemm.cpp.o
dgemm: CMakeFiles/dgemm.dir/build.make
dgemm: CMakeFiles/dgemm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/haofeng-admin/dgemm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dgemm"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dgemm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dgemm.dir/build: dgemm
.PHONY : CMakeFiles/dgemm.dir/build

CMakeFiles/dgemm.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dgemm.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dgemm.dir/clean

CMakeFiles/dgemm.dir/depend:
	cd /Users/haofeng-admin/dgemm/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/haofeng-admin/dgemm /Users/haofeng-admin/dgemm /Users/haofeng-admin/dgemm/cmake-build-debug /Users/haofeng-admin/dgemm/cmake-build-debug /Users/haofeng-admin/dgemm/cmake-build-debug/CMakeFiles/dgemm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dgemm.dir/depend

