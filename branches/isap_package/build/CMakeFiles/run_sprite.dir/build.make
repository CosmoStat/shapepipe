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
include CMakeFiles/run_sprite.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/run_sprite.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/run_sprite.dir/flags.make

CMakeFiles/run_sprite.dir/src/run_sprite.cc.o: CMakeFiles/run_sprite.dir/flags.make
CMakeFiles/run_sprite.dir/src/run_sprite.cc.o: ../src/run_sprite.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/run_sprite.dir/src/run_sprite.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run_sprite.dir/src/run_sprite.cc.o -c /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/run_sprite.cc

CMakeFiles/run_sprite.dir/src/run_sprite.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run_sprite.dir/src/run_sprite.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/run_sprite.cc > CMakeFiles/run_sprite.dir/src/run_sprite.cc.i

CMakeFiles/run_sprite.dir/src/run_sprite.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run_sprite.dir/src/run_sprite.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/run_sprite.cc -o CMakeFiles/run_sprite.dir/src/run_sprite.cc.s

CMakeFiles/run_sprite.dir/src/run_sprite.cc.o.requires:

.PHONY : CMakeFiles/run_sprite.dir/src/run_sprite.cc.o.requires

CMakeFiles/run_sprite.dir/src/run_sprite.cc.o.provides: CMakeFiles/run_sprite.dir/src/run_sprite.cc.o.requires
	$(MAKE) -f CMakeFiles/run_sprite.dir/build.make CMakeFiles/run_sprite.dir/src/run_sprite.cc.o.provides.build
.PHONY : CMakeFiles/run_sprite.dir/src/run_sprite.cc.o.provides

CMakeFiles/run_sprite.dir/src/run_sprite.cc.o.provides.build: CMakeFiles/run_sprite.dir/src/run_sprite.cc.o


# Object files for target run_sprite
run_sprite_OBJECTS = \
"CMakeFiles/run_sprite.dir/src/run_sprite.cc.o"

# External object files for target run_sprite
run_sprite_EXTERNAL_OBJECTS =

run_sprite: CMakeFiles/run_sprite.dir/src/run_sprite.cc.o
run_sprite: CMakeFiles/run_sprite.dir/build.make
run_sprite: libsparse2d.a
run_sprite: libsparse1d.a
run_sprite: libtools.a
run_sprite: CMakeFiles/run_sprite.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable run_sprite"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run_sprite.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/run_sprite.dir/build: run_sprite

.PHONY : CMakeFiles/run_sprite.dir/build

CMakeFiles/run_sprite.dir/requires: CMakeFiles/run_sprite.dir/src/run_sprite.cc.o.requires

.PHONY : CMakeFiles/run_sprite.dir/requires

CMakeFiles/run_sprite.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run_sprite.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run_sprite.dir/clean

CMakeFiles/run_sprite.dir/depend:
	cd /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/CMakeFiles/run_sprite.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/run_sprite.dir/depend

