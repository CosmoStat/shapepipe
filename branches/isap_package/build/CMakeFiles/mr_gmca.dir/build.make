# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/opt/core-4.0-amd64/gnu/cmake/3.5.2-gcc492/bin/cmake

# The command to remove a file.
RM = /usr/local/opt/core-4.0-amd64/gnu/cmake/3.5.2-gcc492/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build

# Include any dependencies generated for this target.
include CMakeFiles/mr_gmca.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mr_gmca.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mr_gmca.dir/flags.make

CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o: CMakeFiles/mr_gmca.dir/flags.make
CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o: ../src/mr_gmca.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o -c /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/mr_gmca.cc

CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/mr_gmca.cc > CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.i

CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/mr_gmca.cc -o CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.s

CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o.requires:

.PHONY : CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o.requires

CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o.provides: CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o.requires
	$(MAKE) -f CMakeFiles/mr_gmca.dir/build.make CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o.provides.build
.PHONY : CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o.provides

CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o.provides.build: CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o


# Object files for target mr_gmca
mr_gmca_OBJECTS = \
"CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o"

# External object files for target mr_gmca
mr_gmca_EXTERNAL_OBJECTS =

mr_gmca: CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o
mr_gmca: CMakeFiles/mr_gmca.dir/build.make
mr_gmca: libsparse2d.a
mr_gmca: libsparse1d.a
mr_gmca: libtools.a
mr_gmca: CMakeFiles/mr_gmca.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mr_gmca"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mr_gmca.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mr_gmca.dir/build: mr_gmca

.PHONY : CMakeFiles/mr_gmca.dir/build

CMakeFiles/mr_gmca.dir/requires: CMakeFiles/mr_gmca.dir/src/mr_gmca.cc.o.requires

.PHONY : CMakeFiles/mr_gmca.dir/requires

CMakeFiles/mr_gmca.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mr_gmca.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mr_gmca.dir/clean

CMakeFiles/mr_gmca.dir/depend:
	cd /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles/mr_gmca.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mr_gmca.dir/depend

