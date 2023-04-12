#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/c4/bin/CMakeLists.txt
# brief  Instructions for building xthi Makefiles
# note   Copyright (C) 2017-2023 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(c4_bin CXX)

# -------------------------------------------------------------------------------------------------#
# Build package library
# -------------------------------------------------------------------------------------------------#

# xthi depends on OpenMP, don't attempt to build it if we don't have it
if(OpenMP_FOUND)
  add_component_executable(
    TARGET Exe_xthi
    TARGET_DEPS Lib_c4;OpenMP::OpenMP_CXX
    SOURCES ${PROJECT_SOURCE_DIR}/xthi.cc)
  install(
    TARGETS Exe_xthi
    EXPORT draco-targets
    DESTINATION ${DBSCFGDIR}bin)
  # rzvernal with cce@15 needs -fopenmp on the link line and cmake doesn't provide this
  # automatically.  Probably a bug, but let's just provide the work-around for now.
  if(DEFINED ENV{CRAY_ACCEL_TARGET}
     AND "$ENV{CRAY_ACCEL_TARGET}" MATCHES "amd"
     AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set_target_properties(
      Exe_xthi
      PROPERTIES
      LINK_FLAGS "-fopenmp")
  endif()
endif()

add_component_executable(
  TARGET Exe_ythi
  TARGET_DEPS Lib_c4;Threads::Threads
  SOURCES ${PROJECT_SOURCE_DIR}/ythi.cc)

# -------------------------------------------------------------------------------------------------#
# Installation instructions
# -------------------------------------------------------------------------------------------------#
install(
  TARGETS Exe_ythi
  EXPORT draco-targets
  DESTINATION ${DBSCFGDIR}bin)

# -------------------------------------------------------------------------------------------------#
# End cdi_eospac/CMakeLists.txt
# -------------------------------------------------------------------------------------------------#
