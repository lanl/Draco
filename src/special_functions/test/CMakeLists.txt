#--------------------------------------------*-cmake-*---------------------------------------------#
# file   special_functions/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Aug 1
# brief  Generate build project files for special_functions/test
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(special_functions_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB test_sources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_special_functions;Lib_diagnostics")

# ------------------------------------------------------------------------------------------------ #
# Copy vendor dll files to test directory (binary_dir)
# ------------------------------------------------------------------------------------------------ #
if(WIN32)
  provide_aux_files(FILES "${GSL_DLL_LIBRARIES}")
endif()

# ------------------------------------------------------------------------------------------------ #
# End special_functions/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
