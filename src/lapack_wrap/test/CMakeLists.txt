#--------------------------------------------------------------------------------------------------#
# file   lapack_wrap/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for lapack_wrap.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(lapack_wrap_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB test_sources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_dsxx;LAPACK::LAPACK")
target_include_directories(Ut_lapack_wrap_tstBlas_Level_1_exe
                           PRIVATE $<BUILD_INTERFACE:${lapack_wrap_BINARY_DIR}>)

# ------------------------------------------------------------------------------------------------ #
# End lapack_wrap/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
