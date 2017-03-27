cmake_minimum_required(VERSION 2.6)
include(FindPkgConfig)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
project(great3)

#
# Find libcfitsio
#
# Locate cfitsio using pkg-config

find_package(Boost 1.55.0 COMPONENTS python)
find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
message("PYTHON: ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${BOOST_PYTHON_INCLUDE_DIR} ${NUMPY_INCLUDE_DIRS})

pkg_check_modules(CFITSIO REQUIRED cfitsio)
include_directories(${CFITSIO_INCLUDE_DIRS})
link_directories(${CFITSIO_LIBRARY_DIRS})

pkg_check_modules(GSL gsl)
if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    link_directories(${GSL_LIBRARY_DIRS})
endif(GSL_FOUND)
#
# Compilation flags
#
if(APPLE)
  set(CMAKE_CXX_FLAGS "-DMACOS -DNO_DISP_IO -g0 -O2 -fPIC -ffast-math -fomit-frame-pointer -fopenmp -Wno-write-strings -DNDEBUG")
#  set(CMAKE_CXX_FLAGS "-DMACOS -g -Wno-write-strings -DNDEBUG")
else(APPLE)
  set(CMAKE_CXX_FLAGS "-DNO_DISP_IO  -ggdb3 -fPIC -O2 -ffast-math -fomit-frame-pointer -fopenmp -Wno-write-strings -DNDEBUG")
endif(APPLE)

#
# Build tools library
#
FILE(GLOB src_lib1  "${PROJECT_SOURCE_DIR}/src/libtools/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/libtools")
add_library(tools STATIC ${src_lib1})
target_link_libraries(tools ${CFITSIO_LIBRARIES})

#
# Build sparse1d library
#
FILE(GLOB src_lib2  "${PROJECT_SOURCE_DIR}/src/libsparse1d/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/libsparse1d")
add_library(sparse1d STATIC ${src_lib2})
target_link_libraries(sparse1d ${CFITSIO_LIBRARIES})

#
# Build sparse2d library
#
FILE(GLOB src_lib2  "${PROJECT_SOURCE_DIR}/src/libsparse2d/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/libsparse2d")
add_library(sparse2d STATIC ${src_lib2})
target_link_libraries(sparse2d ${CFITSIO_LIBRARIES} ${GSL_LIBRARIES})

#
# Compile and link executables
#
# set(BINMR2D     mr1d_filter mr1d_trans mr1d_recons mw1d_filter  mr1d_ifilter mr1d_detect mr1d_dct mr1d_deconv run_sprite)
set(BINMR2D    mr_gmca  im3d_deconv mr_transform mr_recons mr_filter mr_deconv mw_deconv im3d_coadd run_sprite)


foreach(program ${BINMR2D})
 add_executable(${program} ${PROJECT_SOURCE_DIR}/src/${program}.cc)
 target_link_libraries(${program} sparse2d  sparse1d tools)
endforeach(program)

#
# Install (by default in the project directory)
#

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

#install sparse1d headers
FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libsparse1d/*.h")
INSTALL(FILES ${inc_lib} DESTINATION include)
FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libsparse2d/*.h")
INSTALL(FILES ${inc_lib} DESTINATION include)

# Install sparse1d library
INSTALL(TARGETS sparse1d sparse2d tools DESTINATION lib)

# install sparse1d executables
INSTALL(TARGETS ${BINMR2D} DESTINATION bin)

#
# Build python wrappers
#
#Try to generate the boost python code automatically using pyste - else use the cpp provided
find_program(PYSTE_TOOL NAMES pyste PATHS ENV PATH)
if(EXISTS ${PYSTE_TOOL} AND DEFINED PYSTE_PROC)
   message("GENERATE BOOSTPYTHON CODE USING BOOST_PYTHON")
   if(APPLE)
      add_custom_command(OUTPUT isap_great3_tools.cpp COMMAND pyste -I ${PROJECT_SOURCE_DIR}/src/libtools  -DMACOS
					--module=isap_great3_tools ${PROJECT_SOURCE_DIR}/python-wrappers/src/isap_great3_tools.pyste
                                    DEPENDS ${PROJECT_SOURCE_DIR}/python-wrappers/src/isap_great3_tools.pyste)
      add_custom_command(OUTPUT isap_great3_sparse2d.cpp COMMAND pyste -I ${PROJECT_SOURCE_DIR}/src/libsparse1d
                                                    -I ${PROJECT_SOURCE_DIR}/src/libsparse2d
                                                    -I ${PROJECT_SOURCE_DIR}/src/libtools -DMACOS
                                                  --module=isap_great3_sparse2d ${PROJECT_SOURCE_DIR}/python-wrappers/src/isap_great3_sparse2d.pyste
                                    DEPENDS ${PROJECT_SOURCE_DIR}/python-wrappers/src/isap_great3_sparse2d.pyste)
   else(APPLE)
      add_custom_command(OUTPUT isap_great3_tools.cpp COMMAND pyste -I "${PROJECT_SOURCE_DIR}/src/libtools"
					--module=isap_great3_tools ${PROJECT_SOURCE_DIR}/python-wrappers/src/isap_great3_tools.pyste
                                    DEPENDS ${PROJECT_SOURCE_DIR}/python-wrappers/src/isap_great3_tools.pyste)
      add_custom_command(OUTPUT isap_great3_sparse2d.cpp COMMAND pyste -I "${PROJECT_SOURCE_DIR}/src/libsparse1d"
                                                    -I "${PROJECT_SOURCE_DIR}/src/libsparse2d"
                                                    -I "${PROJECT_SOURCE_DIR}/src/libtools"
                                                  --module=isap_great3_sparse2d ${PROJECT_SOURCE_DIR}/python-wrappers/src/isap_great3_sparse2d.pyste
                                    DEPENDS ${PROJECT_SOURCE_DIR}/python-wrappers/src/isap_great3_sparse2d.pyste)
   endif(APPLE)
   add_library(pytools_isap_great3 SHARED isap_great3_tools.cpp)
   add_library(pysparse2d_isap_great3 SHARED isap_great3_sparse2d.cpp)
   INSTALL(FILES build/isap_great3_tools.cpp DESTINATION  python-wrappers/src)
   INSTALL(FILES build/isap_great3_sparse2d.cpp DESTINATION  python-wrappers/src)
else ()
   add_library(pytools_isap_great3 SHARED python-wrappers/src/isap_great3_tools.cpp)
   add_library(pysparse2d_isap_great3 SHARED python-wrappers/src/isap_great3_sparse2d.cpp)
endif()

#Create libraries
target_link_libraries(pytools_isap_great3 tools ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(pytools_isap_great3 PROPERTIES PREFIX "")
set_target_properties(pytools_isap_great3 PROPERTIES OUTPUT_NAME "isap_great3_tools")

target_link_libraries(pysparse2d_isap_great3 sparse2d sparse1d tools ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(pysparse2d_isap_great3 PROPERTIES PREFIX "")
set_target_properties(pysparse2d_isap_great3 PROPERTIES OUTPUT_NAME "isap_great3_sparse2d")

#
# Install python libraries
#
#string(FIND ${PYTHONLIBS_VERSION_STRING} "." minor_version REVERSE)
#string(SUBSTRING ${PYTHONLIBS_VERSION_STRING} 0 ${minor_version} major_version)
INSTALL(TARGETS pytools_isap_great3 DESTINATION lib/python2.7/site-packages)
INSTALL(TARGETS pysparse2d_isap_great3 DESTINATION lib/python2.7/site-packages)

