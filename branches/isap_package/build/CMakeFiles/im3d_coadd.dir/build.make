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
include CMakeFiles/im3d_coadd.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/im3d_coadd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/im3d_coadd.dir/flags.make

CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o: CMakeFiles/im3d_coadd.dir/flags.make
CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o: ../src/im3d_coadd.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o -c /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/im3d_coadd.cc

CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/im3d_coadd.cc > CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.i

CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/im3d_coadd.cc -o CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.s

CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o.requires:

.PHONY : CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o.requires

CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o.provides: CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o.requires
	$(MAKE) -f CMakeFiles/im3d_coadd.dir/build.make CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o.provides.build
.PHONY : CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o.provides

CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o.provides.build: CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o


# Object files for target im3d_coadd
im3d_coadd_OBJECTS = \
"CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o"

# External object files for target im3d_coadd
im3d_coadd_EXTERNAL_OBJECTS =

im3d_coadd: CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o
im3d_coadd: CMakeFiles/im3d_coadd.dir/build.make
im3d_coadd: libsparse2d.a
im3d_coadd: libsparse1d.a
im3d_coadd: libtools.a
im3d_coadd: CMakeFiles/im3d_coadd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable im3d_coadd"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/im3d_coadd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/im3d_coadd.dir/build: im3d_coadd

.PHONY : CMakeFiles/im3d_coadd.dir/build

CMakeFiles/im3d_coadd.dir/requires: CMakeFiles/im3d_coadd.dir/src/im3d_coadd.cc.o.requires

.PHONY : CMakeFiles/im3d_coadd.dir/requires

CMakeFiles/im3d_coadd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/im3d_coadd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/im3d_coadd.dir/clean

CMakeFiles/im3d_coadd.dir/depend:
	cd /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles/im3d_coadd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/im3d_coadd.dir/depend

