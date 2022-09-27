#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/c4/bin/CMakeLists.txt
# brief  Instructions for building xthi Makefiles
# note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved.
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
    TARGET_DEPS Lib_c4
    SOURCES ${PROJECT_SOURCE_DIR}/xthi.cc)
  install(
    TARGETS Exe_xthi
    EXPORT draco-targets
    DESTINATION ${DBSCFGDIR}bin)
endif()

add_component_executable(
  TARGET Exe_ythi
  TARGET_DEPS Lib_c4
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
