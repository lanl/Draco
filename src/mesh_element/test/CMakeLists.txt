#--------------------------------------------*-cmake-*---------------------------------------------#
# file   mesh_element/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for mesh_element/test.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(mesh_element_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

set(test_lib_sources ${PROJECT_SOURCE_DIR}/mesh_element_test.cc
                     ${PROJECT_SOURCE_DIR}/mesh_element_test.hh)
file(GLOB test_sources *.cc)
list(REMOVE_ITEM test_sources ${test_lib_sources})

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_mesh_element;Lib_dsxx")

# ------------------------------------------------------------------------------------------------ #
# End mesh_element/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
