#--------------------------------------------*-cmake-*---------------------------------------------#
# file   diagnostics/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for diagnostics/test.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(diagnostics_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(scalar_test_sources ${PROJECT_SOURCE_DIR}/tstDiagnostics.cc ${PROJECT_SOURCE_DIR}/tstProcmon.cc
                        ${PROJECT_SOURCE_DIR}/tstTiming.cc)
set(parallel_test_sources ${PROJECT_SOURCE_DIR}/tstruntime_error.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_scalar_tests(SOURCES "${scalar_test_sources}" DEPS Lib_diagnostics)
add_parallel_tests(
  SOURCES "${parallel_test_sources}"
  DEPS Lib_diagnostics
  PE_LIST "1")

include(ApplicationUnitTest)

add_app_unit_test(
  DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tDracoInfo.py
  APP $<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
  TEST_ARGS "none;--author;--brief;--use_doxygen_formatting;--help;--version;--foo"
  LABELS nomemcheck)

# ------------------------------------------------------------------------------------------------ #
# End diagnostics/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
