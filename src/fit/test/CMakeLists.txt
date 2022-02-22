#--------------------------------------------*-cmake-*---------------------------------------------#
# file   fit/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for fit/test.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(fit_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_lib_sources ${PROJECT_SOURCE_DIR}/fit_test.cc ${PROJECT_SOURCE_DIR}/fit_test.hh)
file(GLOB test_sources *.cc)
list(REMOVE_ITEM test_sources ${test_lib_sources})
file(GLOB mesh_input_files *.in)
file(GLOB python_scripts python/*.py)
file(GLOB bench_files bench/*.0001)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_linear")

# ------------------------------------------------------------------------------------------------ #
# End fit/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
