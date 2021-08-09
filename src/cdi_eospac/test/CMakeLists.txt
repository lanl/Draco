#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_eospac/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_eospac/test.
# note   Copyright (C) 2016-2021, Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(cdi_eospac_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_sources ${PROJECT_SOURCE_DIR}/tEospac.cc ${PROJECT_SOURCE_DIR}/tEospacWithCDI.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
if((DEFINED ENV{SESAMEPATH} AND EXISTS "$ENV{SESAMEPATH}")
   OR (DEFINED ENV{SESPATHU} AND EXISTS "$ENV{SESPATHU}")
   OR IS_DIRECTORY "/usr/gapps/lanl-data/eos"
   OR IS_DIRECTORY "/usr/projects/data/eos")

  # These tests require a datafile to be defined in the environment.

  add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_cdi_eospac")

  include(ApplicationUnitTest)

  # Run the QueryEospac executable reading interactive commands directly from the file
  # QueryEospacInput.dat.  Save the output to QueryEospacInput.out and compare this to the gold
  # file.
  add_app_unit_test(
    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tQueryEospac.cmake
    APP $<TARGET_FILE_DIR:Exe_QueryEospac>/$<TARGET_FILE_NAME:Exe_QueryEospac>
    STDINFILE ${PROJECT_SOURCE_DIR}/QueryEospacInput.dat
    GOLDFILE ${PROJECT_SOURCE_DIR}/QueryEospac.gold)
endif()

# Run the QueryEospac executable with options '--version' and '--help' to ensure these options print
# the expected output.
add_app_unit_test(
  DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tQueryEospac.cmake
  APP $<TARGET_FILE_DIR:Exe_QueryEospac>/$<TARGET_FILE_NAME:Exe_QueryEospac>
  TEST_ARGS "--version;--help")

# ------------------------------------------------------------------------------------------------ #
# End cdi_eospac/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
