#--------------------------------------------*-cmake-*---------------------------------------------#
# file   compton_interface/test/CMakeLists.txt
# author Kendra Keady
# date   2017 February 28
# brief  Generate build project files for compton/test.
# note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(compton_interface_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(interfaces_test_sources ${PROJECT_SOURCE_DIR}/tCSK_Interface_interfaces.cc)
set(dummy_data_test_sources ${PROJECT_SOURCE_DIR}/tCSK_Interface_dummy_data.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

# If the compton library was not found, then there is no Lib_compton_interface library to link
# against (and thus compton_interface/config.h isn't found). However, we can still build the test
# against the stub headers...
if(NOT TARGET Lib_compton_interface)

  add_parallel_tests(
    SOURCES "${interfaces_test_sources}"
    PE_LIST "1"
    DEPS "Lib_c4"
    LABEL "nomemcheck"
    MPI_PLUS_OMP)

  # generated include directive files (config.h)
  target_include_directories(Ut_compton_interface_tCSK_Interface_interfaces_exe
                             PUBLIC $<BUILD_INTERFACE:${compton_interface_BINARY_DIR}>)

else()

  add_parallel_tests(
    MPI_PLUS_OMP
    SOURCES "${interfaces_test_sources}"
    PE_LIST "1"
    DEPS "Lib_compton_interface")

  add_parallel_tests(
    MPI_PLUS_OMP
    SOURCES "${dummy_data_test_sources}"
    PE_LIST "1"
    DEPS "Lib_compton_interface")

endif()

# ------------------------------------------------------------------------------------------------ #
# end of compton_interface/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
