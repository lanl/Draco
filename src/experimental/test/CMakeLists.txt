#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/experimental/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   Saturday, Dec 28, 2019, 15:50 pm
# brief  Instructions for building experimental/test level Makefile.
# note   Copyright (C) 2019-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(experimental_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

file(GLOB test_sources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build and register unit tests
# ------------------------------------------------------------------------------------------------ #

add_scalar_tests(SOURCES "${test_sources}" DEPS Lib_dsxx)
# target_include_directories( Ut_experimental_tst_mdspan_exe PUBLIC
# $<BUILD_INTERFACE:${draco_src_dir_SOURCE_DIR}> )

# ------------------------------------------------------------------------------------------------ #
# end experimental/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
