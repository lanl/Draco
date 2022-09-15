#--------------------------------------------*-cmake-*---------------------------------------------#
# file   ds++/bin/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 Nov 6
# brief  Provide a sample crash handler for Windows.
# note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(dsxx_bin CXX)

# -------------------------------------------------------------------------------------------------#
# Build and install programs
# -------------------------------------------------------------------------------------------------#
if(WIN32)
  add_component_executable(
    TARGET Exe_testCrashHandler
    TARGET_DEPS Lib_dsxx
    SOURCES ${PROJECT_SOURCE_DIR}/testCrashHandler.cc)
  install(
    TARGETS Exe_testCrashHandler
    EXPORT draco-targets
    DESTINATION ${DBSCFGDIR}bin)
endif()

add_component_executable(
  TARGET Exe_query_cpu
  TARGET_DEPS Lib_dsxx
  SOURCES ${Draco_SOURCE_DIR}/config/query_fma.cc)
install(
  TARGETS Exe_query_cpu
  EXPORT draco-targets
  DESTINATION ${DBSCFGDIR}bin)

# -------------------------------------------------------------------------------------------------#
# end ds++/bin/CMakeLists.txt
# -------------------------------------------------------------------------------------------------#
