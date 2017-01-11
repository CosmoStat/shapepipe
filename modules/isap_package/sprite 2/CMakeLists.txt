cmake_minimum_required(VERSION 2.6)
include(FindPkgConfig)


project(test_sr)
#find_package(PythonLibs REQUIRED)
#find_package(Boost 1.45.0 COMPONENTS python)
#find_package(NumPy REQUIRED)
#SET(Boost_USE_STATIC_LIBS OFF)
#SET(Boost_USE_MULTITHREADED ON)
#SET(Boost_USE_STATIC_RUNTIME OFF)
#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS})

#
# Find libcfitsio
#
# Locate cfitsio using pkg-config
pkg_check_modules(CFITSIO REQUIRED cfitsio)
include_directories(${CFITSIO_INCLUDE_DIRS})
link_directories(${CFITSIO_LIBRARY_DIRS})

pkg_check_modules(GSL REQUIRED gsl)
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})


#set(BIN_DIR  bin)
#if (DEFINED ENV{ISAP})  
#set(BIN_DIR  $ENV{ISAP}/bin)
#endif  (DEFINED ENV{ISAP}) 

#
# Compilation flags
#
if(APPLE)
  set(CMAKE_CXX_FLAGS "-DMACOS -DNO_DISP_IO -g0 -O2 -ffast-math -fomit-frame-pointer -fopenmp -Wno-write-strings -DNDEBUG")
#  set(CMAKE_CXX_FLAGS "-DMACOS -g -Wno-write-strings -DNDEBUG")
else(APPLE)
  set(CMAKE_CXX_FLAGS "-DNO_DISP_IO  -g0 -O2 -ffast-math -fomit-frame-pointer -fopenmp -Wno-write-strings -DNDEBUG")
endif(APPLE)

#
# Build tools library
#
FILE(GLOB src_lib1  "${PROJECT_SOURCE_DIR}/src/libtools/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/libtools")
add_library(tools STATIC ${src_lib1})
target_link_libraries(tools ${CFITSIO_LIBRARIES} ${GSL_LIBRARIES})

#
# Build sparse1d library
#
FILE(GLOB src_lib2  "${PROJECT_SOURCE_DIR}/src/libsparse1d/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/libsparse1d")
add_library(sparse1d STATIC ${src_lib2})
target_link_libraries(sparse1d ${CFITSIO_LIBRARIES})

#
# Build sparse1d library
#
FILE(GLOB src_lib2  "${PROJECT_SOURCE_DIR}/src/libsparse2d/*.cc")
include_directories("${PROJECT_SOURCE_DIR}/src/libsparse2d")
add_library(sparse2d STATIC ${src_lib2})
target_link_libraries(sparse2d ${CFITSIO_LIBRARIES})

#
# Compile and link executables
#
# set(BINMR2D     mr1d_filter mr1d_trans mr1d_recons mw1d_filter  mr1d_ifilter mr1d_detect mr1d_dct mr1d_deconv)
# set(BINMR2D    test_sr mr_transform im3d_coadd run_sprite op_mat_test)
set(BINMR2D    run_sprite  )

# set(BINMR2D im_deconv  mr_abaque  mr_at_edge  mr_background    mr_comp  mr_compare   mr_contrast   mr_decomp  mr_deconv  mr_detect  mr_edge  mr_extract  mr_filter mr_fusion mr_getresol mr_hfilter mr_info mr_insert mr_lcomp mr_mrc mr_pfilter mr_psupport mr_rec_edge mr_recons mr_rfilter mr_segment mr_support mr_transform mr_upresol mr_visu mw_comb_filter mw_deconv mw_entrop mw_filter mw_proba mw_stat  wcomp wdecomp mr_linecol mr_ulc_filter)


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
