#--------------------------------------------*-cmake-*---------------------------------------------#
# file   rng/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng/test.
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(rng_test C CXX)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.0)
  string(APPEND CMAKE_C_FLAGS " -Wno-expansion-to-defined")
  string(APPEND CMAKE_CXX_FLAGS " -Wno-expansion-to-defined")
endif()

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_sources
    ${PROJECT_SOURCE_DIR}/tstRnd_Control_Inline.cc ${PROJECT_SOURCE_DIR}/tstSubrandom_Sequence.cc
    ${PROJECT_SOURCE_DIR}/tstCounter_RNG.cc)

# Random123 unit tests (these tests have special PASS/FAIL REGEX conditions)
set(random123_unit_tests
    ${PROJECT_SOURCE_DIR}/ut_M128.cpp ${PROJECT_SOURCE_DIR}/ut_aes.cpp
    ${PROJECT_SOURCE_DIR}/ut_ars.c ${PROJECT_SOURCE_DIR}/ut_carray.cpp
    ${PROJECT_SOURCE_DIR}/ut_gsl.c ${PROJECT_SOURCE_DIR}/ut_uniform.cpp)

# [2020-08-21 KT] The code in ut_gsl.c is not valid C++.  Because this test is a mix of C and C++,
# it causes issues. Most compilers are happy to compile this code with C. However, LLVM clang and
# Visual Studio cl.exe report fatal errors.  When marked as C++, this test fails with:
#
# "C:\projects\draco\src\rng\test\ut_gsl.c(62,1): error C2086: 'const gsl_rng_type \
# *gsl_rng_cbrng': redefinition"
#
# Indeed, reviewing the Random123 file gsl_cbrng.h reveals that that gsl_rng_type is declared twice,
# first on line 68 and then again on line 125 (with an assignment this time). We can continue
# building and running this tests on platforms that can processs it (all but MSVC right now).
if(CMAKE_CXX_COMPILER_ID MATCHES Clang OR CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  list(REMOVE_ITEM random123_unit_tests ${PROJECT_SOURCE_DIR}/ut_gsl.c)
else()
  if(R123_VERSION VERSION_GREATER 1.09)
    list(APPEND random123_unit_tests ${PROJECT_SOURCE_DIR}/ut_Engine.cpp)
  endif()
endif()

# Random123 known-answer tests
set(random123_known_answer_tests ${PROJECT_SOURCE_DIR}/kat_c.cc ${PROJECT_SOURCE_DIR}/kat_cpp.cpp)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #
add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_rng")

add_scalar_tests(
  SOURCES "${random123_unit_tests}"
  DEPS "Lib_rng"
  PASS_REGEX "OK|SUCCESS" FAIL_REGEX "failed;error;Failure")

target_include_directories(Ut_rng_ut_uniform_exe PRIVATE $<BUILD_INTERFACE:${rng_SOURCE_DIR}>)

add_scalar_tests(
  SOURCES "${random123_known_answer_tests}"
  DEPS "Lib_rng"
  PASS_REGEX "PASSED" FAIL_REGEX "FAILED")

if(USE_CUDA AND NOT CMAKE_CXX_COMPILER_ID MATCHES XL)
  # This test is proved by the Random123 project and is requires to nvcc (cuda) and gcc
  set(test_sources_cuda ${PROJECT_SOURCE_DIR}/kat_cuda.cu)
  set(test_deps_cuda "Lib_rng")
  add_scalar_tests(
    SOURCES "${test_sources_cuda}"
    DEPS "${test_deps_cuda}"
    PASS_REGEX "PASSED" FAIL_REGEX "FAILED")
endif()

add_component_executable(TARGET Exe_time_serial SOURCES ${PROJECT_SOURCE_DIR}/time_serial.cc
                                                        NOEXPORT)
target_include_directories(
  Exe_time_serial PRIVATE $<BUILD_INTERFACE:${rng_BINARY_DIR}>
                          $<BUILD_INTERFACE:${dsxx_BINARY_DIR}> ${RANDOM123_INCLUDE_DIR})

configure_file(${PROJECT_SOURCE_DIR}/kat_vectors ${PROJECT_BINARY_DIR}/kat_vectors COPYONLY)

# ------------------------------------------------------------------------------------------------ #
# end rng/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
