#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/ds++/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building ds++/test level Makefile.
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(dsxx_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

file(GLOB test_sources *.cc)
list(REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/do_exception.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

# Run ScalarUnitTest in 2 modes: (a) without additional command line arguments and (b) a special
# test to demonstrate the optional '--version' command line argument
add_scalar_tests(
  SOURCES ${PROJECT_SOURCE_DIR}/tstScalarUnitTest.cc
  DEPS Lib_dsxx
  # (a) Run first case without arguments; (b) Run second case with the '--version' argument.
  TEST_ARGS " ;--version"
  # ~~~
  #         ^^^ this space before the semicolon is important!
  # The 2nd case does not print "PASS" so we register the following regex as PASS message.
  # ~~~
  PASS_REGEX "tstScalarUnitTest[.ex]*: version")
list(REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/tstScalarUnitTest.cc)

# The remaining tests
add_scalar_tests(SOURCES "${test_sources}" DEPS Lib_dsxx)

# (2013-11-13 KT):: tstdbc intentionally divides by zero (to test the isinfinity function). PGI
# issues a warning for this case that we don't care about.  Since PGI does not appear to use pragmas
# to suppress warnings, I am turning off warnings for Ut_dsxx_tstdbc_exe when the compiler is PGI.
#
# ~~~
# if( draco_isPGI )
#   set_source_files_properties(${PROJECT_SOURCE_DIR}/tstdbc.cc PROPERTIES
#     COMPILE_FLAGS --no_warnings)
# endif()
# ~~~

# tstData_Table.cc will fail for RelWithDebugInfo. Set extra CPP macros as a work-around:
set_target_properties(
  Ut_dsxx_tstData_Table_exe
  PROPERTIES
  COMPILE_DEFINITIONS $<UPPER_CASE:$<CONFIG>>)

# create the do_exception executable.
add_executable(Exe_do_exception ${PROJECT_SOURCE_DIR}/do_exception.cc)
target_link_libraries(Exe_do_exception Lib_dsxx)
set_target_properties(
  Exe_do_exception PROPERTIES
  FOLDER dsxx_test
  OUTPUT_NAME do_exception)
dbs_std_tgt_props(Exe_do_exception)

# This is currently broken for clang-cl under MSVC because we are picking up <float.h> from LLVM
# instead of from MSVC toolkits.
if(NOT (MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
  include(ApplicationUnitTest)
  add_app_unit_test(
    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/test_fpe_trap.py
    APP $<TARGET_FILE_DIR:Exe_do_exception>/$<TARGET_FILE_NAME:Exe_do_exception>
    TEST_ARGS "0;1;2;3" #< These 4 are: (0) no exceptions, (1) division by zero, (2) sqrt(-1.0) and
                        #< (3) overflow
    LABELS "nomemcheck")
endif()

# ------------------------------------------------------------------------------------------------ #
# end ds++/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
