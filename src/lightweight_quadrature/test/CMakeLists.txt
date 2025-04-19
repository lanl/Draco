#--------------------------------------------*-cmake-*---------------------------------------------#
# file   lightweight_quadrature/test/CMakeLists.txt
# author Jim Warsa <warsa@lanl.gov>
# date   2025 Apr 1
# brief  Generate build project files for lightweight_quadrature/test.
# note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(lightweight_quadrature_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB test_sources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_lightweight_quadrature;Lib_diagnostics")

# ------------------------------------------------------------------------------------------------ #
# Copy vendor dll files to test directory (binary_dir)
# ------------------------------------------------------------------------------------------------ #
if(WIN32)
  provide_aux_files(FILES "${GSL_DLL_LIBRARIES}")
endif()

# ------------------------------------------------------------------------------------------------ #
# End lighweight_quadrature/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
