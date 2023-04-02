# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   parser/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for parser/test.
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
project(parser_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB test_sources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

add_component_executable(
  TARGET Exe_tstConsole_Token_Stream
  SOURCES ${PROJECT_SOURCE_DIR}/tstConsole_Token_Stream.cc
  TARGET_DEPS Lib_parser NOEXPORT)
list(REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/tstConsole_Token_Stream.cc)

# These tests are scalar only
set(scalar_test_sources
    ${PROJECT_SOURCE_DIR}/tstToken.cc
    ${PROJECT_SOURCE_DIR}/tstToken_Equivalence.cc
    ${PROJECT_SOURCE_DIR}/tstutilities.cc
    ${PROJECT_SOURCE_DIR}/tstUnit.cc
    ${PROJECT_SOURCE_DIR}/tstString_Token_Stream.cc
    ${PROJECT_SOURCE_DIR}/tstParse_Table.cc
    ${PROJECT_SOURCE_DIR}/tstClass_Parser.cc)
list(REMOVE_ITEM test_sources ${scalar_test_sources})
add_scalar_tests(SOURCES "${scalar_test_sources}" DEPS "Lib_parser")

# Special case for driver4tstConsole_Token_Stream 1. This test doesn't work under Xcode 2. This test
# must have unique access to tstConsole_Token_Stream and we will use a RESOURCE_LOCK to prevent
# multiple simultaneous instances during parallel ctest.
if(NOT ${CMAKE_GENERATOR} MATCHES "Xcode")
  include(ApplicationUnitTest)
  # cmake-lint: disable=C0301
  add_app_unit_test(
    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/driver4tstConsole_Token_Stream.py
    APP $<TARGET_FILE_DIR:Exe_tstConsole_Token_Stream>/$<TARGET_FILE_NAME:Exe_tstConsole_Token_Stream>
    STDINFILE "${CMAKE_CURRENT_SOURCE_DIR}/console_test.inp"
    LABELS "nomemcheck")
endif()

# ------------------------------------------------------------------------------------------------ #
# The remaining tests are Parallel.
# ------------------------------------------------------------------------------------------------ #

# These 2 tests require extra environment variables to be set.
list(REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/tstFile_Token_Stream.cc
     ${PROJECT_SOURCE_DIR}/tstParallel_File_Token_Stream.cc)
add_parallel_tests(
  SOURCES "${PROJECT_SOURCE_DIR}/tstFile_Token_Stream.cc"
  PE_LIST "1"
  DEPS "Lib_parser" ENV "DRACO_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}")
add_parallel_tests(
  SOURCES "${PROJECT_SOURCE_DIR}/tstParallel_File_Token_Stream.cc"
  PE_LIST "1;2"
  DEPS "Lib_parser" ENV "DRACO_INCLUDE_PATH=${CMAKE_CURRENT_SOURCE_DIR}")

# Run remaining tests with 1 and 2 PEs.
add_parallel_tests(
  SOURCES "${test_sources}"
  PE_LIST "1;2"
  DEPS "Lib_parser")

# ------------------------------------------------------------------------------------------------ #
# End parser/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
