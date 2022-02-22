#--------------------------------------------*-cmake-*---------------------------------------------#
# file   quadrature/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for quadrature/test.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(quadrature_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_lib_sources ${PROJECT_SOURCE_DIR}/quadrature_test.cc
                     ${PROJECT_SOURCE_DIR}/quadrature_test.hh)
file(GLOB test_sources *.cc)
list(REMOVE_ITEM test_sources ${test_lib_sources})
# Remove some tests from the generic list. These will be process individually with custom build
# options
list(REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/ftstquadrature_interfaces.cc)
# file( GLOB input_files *.inp )

# ------------------------------------------------------------------------------------------------ #
# Build library for test directory
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_quadrature_test
  TARGET_DEPS Lib_quadrature
  LIBRARY_NAME quadrature_test
  SOURCES "${test_lib_sources}"
  NOEXPORT)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_scalar_tests(SOURCES "${test_sources}" DEPS Lib_quadrature_test)

# ------------------------------------------------------------------------------------------------ #
# Wrap f90 tests
# ------------------------------------------------------------------------------------------------ #
if(HAVE_Fortran)
  add_scalar_tests(SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ftstquadrature_interfaces.cc"
                   DEPS "Lib_quadrature_ftest;Lib_quadrature_test;Lib_c4;Lib_dsxx")
endif()

# ------------------------------------------------------------------------------------------------ #
# End quadrature/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
