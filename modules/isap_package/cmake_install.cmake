# Install script for directory: /dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/Mr_FE.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/WT1D_FFT.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR3D_Obj.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/SB_Filter1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR2D1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR_Threshold.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR1D_NoiseModel.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/Atrou1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR1D_Obj.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/Atrou3D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR1D_Filter.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/Fista.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/GMCA.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/WT2D_CF.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/fractal.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR1D_Segment.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR1D_Sigma.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MeyerWT1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/MR1D_Regul.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/Mr_FewEvent.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/Pyr1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/Mr_FewEvent1d.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse1d/Filter.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_CorrNoise.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/Mr_FewEvent2d.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Rayleigh.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Filter.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/sprite.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Psupport.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/IM_Deconv.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/FCur.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Abaque.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MW_Deconv.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MeyerWT.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/SB_Filter.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Obj.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/sr_util.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Noise.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Contrast.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_SoftRegul.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/CMem.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Sigma.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/LineCol.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_NoiseModel.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MW_Filter.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Deconv.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR1D1D.h"
    "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/src/libsparse2d/MR_Support.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/libsparse1d.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/libsparse2d.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/libtools.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/mr_gmca")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_gmca")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/im3d_deconv")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_deconv")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/mr_transform")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_transform")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/mr_recons")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_recons")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/mr_filter")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_filter")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/mr_deconv")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mr_deconv")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/mw_deconv")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/mw_deconv")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/im3d_coadd")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/im3d_coadd")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/run_sprite")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/run_sprite")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages" TYPE SHARED_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/isap_great3_tools.so")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_tools.so")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages" TYPE SHARED_LIBRARY FILES "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/isap_great3_sparse2d.so")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages/isap_great3_sparse2d.so")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/dsm/cosmo02/sparseastro/Euclid/code/EPFL_SAC_great3/trunk/isap-great3-v1.0/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
