#--------------------------------------------*-cmake-*---------------------------------------------#
# file   c4/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for c4/test.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(c4_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB test_sources *.cc)

# Remove some tests from the generic list. These will be process individually with custom build
# options
list(REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/tstfc4_hw.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

# tstAbort doesn't really work as a unit test because it always calls MPI_abort.  This causes an
# uncatchable error condition.  I suppose we could wrap the test in a python or cmake script - maybe
# we can do this sometime in the future.
list(REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/tstAbort.cc)

# This test has special arguments and special pass/fail criteria.
set(special_test ${PROJECT_SOURCE_DIR}/phw.cc)
add_parallel_tests(
  SOURCES "${special_test}"
  DEPS Lib_c4
  PE_LIST "1;2;4"
  TEST_ARGS "hello"
  PASS_REGEX "Hello"
  FAIL_REGEX "cruel")
list(REMOVE_ITEM test_sources ${special_test})

# This tests runs with extra arguments
set(special_test ${PROJECT_SOURCE_DIR}/tstParallelUnitTest.cc)
add_parallel_tests(
  SOURCES "${special_test}"
  DEPS Lib_c4
  PE_LIST "1;2;4"
  TEST_ARGS "--version;--dummy; "
  PASS_REGEX "tstParallelUnitTest[.ex]*: version")
list(REMOVE_ITEM test_sources ${special_test})

# This test uses MPI+OMP:
set(special_tests ${PROJECT_SOURCE_DIR}/tstOMP.cc ${PROJECT_SOURCE_DIR}/tstOMP_API_off.cc
                  ${PROJECT_SOURCE_DIR}/tstOMP_API_on.cc)
if(DEFINED ENV{TRAVIS})
  set(omp_test_pe_list "1") #< TRAVIS has severe limits on number of ranks/threads
elseif(CODE_COVERAGE)
  set(omp_test_pe_list "2") #< The 1 MPI rank version gets confused when running under code-coverage
                            #< configurations.
else()
  set(omp_test_pe_list "1;2")
endif()

foreach(omptest ${special_tests})
  list(REMOVE_ITEM test_sources ${omptest})
  if(OpenMP_FOUND)
    add_parallel_tests(
      MPI_PLUS_OMP
      SOURCES "${omptest}"
      DEPS Lib_c4
      PE_LIST "${omp_test_pe_list}")
  endif()
endforeach()
unset(omp_test_pe_list)

# Tests that expect exactly 2 PE (omit for DRACO_C4=SCALAR builds):
set(two_pe_tests ${PROJECT_SOURCE_DIR}/tstPingPong.cc)
list(APPEND two_pe_tests ${PROJECT_SOURCE_DIR}/tstC4_Status.cc)
list(REMOVE_ITEM test_sources ${two_pe_tests})
if(NOT ${DRACO_C4} STREQUAL "SCALAR")
  add_parallel_tests(
    SOURCES "${two_pe_tests}"
    DEPS Lib_c4
    PE_LIST "2")
endif()

# Tests that expect exactly 4 PE (omit for DRACO_C4=SCALAR builds):
set(four_pe_tests ${PROJECT_SOURCE_DIR}/tstQuoWrapper.cc)
list(REMOVE_ITEM test_sources ${four_pe_tests})

if(NOT ${DRACO_C4} STREQUAL "SCALAR")
  add_parallel_tests(
    SOURCES "${four_pe_tests}"
    DEPS Lib_c4
    PE_LIST "4"
    LABEL "nomemcheck")
  if(NOT DEFINED MPI_MAX_NUMPROCS_PHYSICAL)
    message(FATAL_ERROR "need max procs in c4/tests/CMakeLists.txt")
  endif()
  # Override PROCESSORS property.  Use 'mpirun -np 4 ...' but reserve half of the cores. The test
  # will actually use all of the cores, but only very briefly.
  math(EXPR half_of_all_cores "${MPI_MAX_NUMPROCS_PHYSICAL} / 2")
  set_tests_properties(
    c4_tstQuoWrapper_4
    PROPERTIES
    PROCESSORS "${half_of_all_cores}")
endif()

# Exclude tstTime from valgrind checks.
set(nomemcheck_tests ${PROJECT_SOURCE_DIR}/tstTime.cc ${PROJECT_SOURCE_DIR}/tstGlobalTimer.cc)
list(REMOVE_ITEM test_sources ${nomemcheck_tests})
add_parallel_tests(
  SOURCES "${nomemcheck_tests}"
  DEPS Lib_c4
  PE_LIST "1;2;4"
  LABEL "nomemcheck")

# The remainder of the tests
add_parallel_tests(
  SOURCES "${test_sources}"
  DEPS Lib_c4
  PE_LIST "1;2;4")

# ------------------------------------------------------------------------------------------------ #
# Python based test runners
# ------------------------------------------------------------------------------------------------ #
include(ApplicationUnitTest)
if(OpenMP_FOUND AND TARGET Exe_xthi)
  # xthi depends on OpenMP.  If OpenMP is not found, xthi will not be compiled.
  add_app_unit_test(
    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tstXthi.py
    APP $<TARGET_FILE_DIR:Exe_xthi>/$<TARGET_FILE_NAME:Exe_xthi>
    PE_LIST "1;2"
    LABELS "nomemcheck")
endif()
if(TARGET Exe_ythi)
  add_app_unit_test(
    DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tstYthi.py
    APP $<TARGET_FILE_DIR:Exe_ythi>/$<TARGET_FILE_NAME:Exe_ythi>
    TEST_ARGS 2
    PE_LIST "1;2"
    LABELS "nomemcheck")
endif()

# ------------------------------------------------------------------------------------------------ #
# Wrap f90 tests
# ------------------------------------------------------------------------------------------------ #
if(HAVE_Fortran AND DRACO_C4 STREQUAL MPI)
  add_parallel_tests(
    SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tstfc4_hw.cc"
    DEPS "Lib_c4_ftest;Lib_c4_fc4;Lib_c4;Lib_dsxx"
    PE_LIST "1;2")
endif()

# ------------------------------------------------------------------------------------------------ #
# End src/c4/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
