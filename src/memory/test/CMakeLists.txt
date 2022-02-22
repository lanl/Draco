#--------------------------------------------*-cmake-*---------------------------------------------#
# file   memory/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for memory/test.
# note   Copyright (C) 2013-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(memory_test CXX)

if(NOT "${DRACO_DIAGNOSTICS}" STREQUAL "7")
  # This test fails when fpe_trap is enabled.
  #
  add_scalar_tests(SOURCES "tstmemory.cc" DEPS Lib_memory LABEL "nomemcheck")
endif()

# ------------------------------------------------------------------------------------------------ #
# End memory/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
