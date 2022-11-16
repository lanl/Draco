#--------------------------------------------*-cmake-*---------------------------------------------#
# file   compton_tools/test/CMakeLists.txt
# author Andrew Till
# date   11 May 2020
# brief  Generate build project files for compton_tools/test.
# note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(compton_tools_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_sources ${PROJECT_SOURCE_DIR}/tCompton_Native.cc ${PROJECT_SOURCE_DIR}/tCompton_Edep.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_parallel_tests(
  SOURCES "${test_sources}"
  PE_LIST "1"
  DEPS "Lib_compton_tools"
  MPI_PLUS_OMP)

# ------------------------------------------------------------------------------------------------ #
# Python based test runners
# ------------------------------------------------------------------------------------------------ #
include(ApplicationUnitTest)
if(TARGET Exe_CSKrw)
  add_app_unit_test(
    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tstCSKrw.py
    APP $<TARGET_FILE_DIR:Exe_CSKrw>/$<TARGET_FILE_NAME:Exe_CSKrw>
    TEST_ARGS "${CMAKE_CURRENT_SOURCE_DIR}/dummy_data"
    LABELS "nomemcheck")
endif()

# ------------------------------------------------------------------------------------------------ #
# end of compton_tools/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
