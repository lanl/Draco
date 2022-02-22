#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi/test
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(cdi_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_lib_sources
    ${PROJECT_SOURCE_DIR}/DummyEoS.cc ${PROJECT_SOURCE_DIR}/DummyEICoupling.cc
    ${PROJECT_SOURCE_DIR}/DummyGrayOpacity.cc ${PROJECT_SOURCE_DIR}/DummyMultigroupOpacity.cc)
file(GLOB test_sources *.cc)
file(GLOB test_headers *.hh)
list(REMOVE_ITEM test_sources ${test_lib_sources})

# ------------------------------------------------------------------------------------------------ #
# Build library for test directory
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_cdi_test
  TARGET_DEPS Lib_cdi
  LIBRARY_NAME cdi_test
  SOURCES "${test_lib_sources}"
  NOEXPORT)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_cdi_test")

# ------------------------------------------------------------------------------------------------ #
# End cdi/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
