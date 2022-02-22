#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_analytic/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_analytic/test.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(cdi_analytic_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB test_sources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
set(parallel_test ${PROJECT_SOURCE_DIR}/tstPseudo_Line_Analytic_MultigroupOpacity.cc)
add_parallel_tests(
  SOURCES "${parallel_test}"
  DEPS "Lib_cdi_analytic"
  PE_LIST "1")

list(REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/tstPseudo_Line_Analytic_MultigroupOpacity.cc)

add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_cdi_analytic")

# ------------------------------------------------------------------------------------------------ #
# End cdi_analytic/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
