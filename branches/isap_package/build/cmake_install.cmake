# Install script for directory: /dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/Mr_FE.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/WT1D_FFT.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR3D_Obj.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/SB_Filter1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR2D1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR_Threshold.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR1D_NoiseModel.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/Atrou1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR1D_Obj.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/Atrou3D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR1D_Filter.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/Fista.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/GMCA.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/WT2D_CF.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/fractal.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR1D_Segment.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR1D_Sigma.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MeyerWT1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/MR1D_Regul.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/Mr_FewEvent.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/Pyr1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/Mr_FewEvent1d.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse1d/Filter.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_CorrNoise.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/Mr_FewEvent2d.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Rayleigh.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Filter.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/sprite.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Psupport.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/IM_Deconv.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/FCur.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Abaque.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MW_Deconv.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MeyerWT.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/SB_Filter.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Obj.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/sr_util.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Noise.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Contrast.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_SoftRegul.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/CMem.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Sigma.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/LineCol.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_NoiseModel.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MW_Filter.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Deconv.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR1D1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/src/libsparse2d/MR_Support.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/libsparse1d.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/libsparse2d.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/libtools.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/mr_gmca")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/im3d_deconv")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/mr_transform")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/mr_recons")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/mr_filter")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/mr_deconv")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/mw_deconv")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/im3d_coadd")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/run_sprite")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages" TYPE SHARED_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/isap_great3_tools.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages" TYPE SHARED_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/isap_great3_sparse2d.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/dsm/cosmo02/sparseastro/Euclid/code/cea-epfl/branches/isap_package/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
