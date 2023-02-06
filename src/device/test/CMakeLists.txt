#--------------------------------------------*-cmake-*---------------------------------------------#
# file   device/test/CMakeLists.txt
# brief  Instructions for building device/test Makefile.
# note   Copyright (C) 2011-2023 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(device_test CXX ${GPU_DBS_STRING})

# Only build GPU tests if GPU is available.
if(HAVE_GPU AND USE_GPU)

  # ---------------------------------------------------------------------------------------------- #
  # Source files
  # ---------------------------------------------------------------------------------------------- #
  set(test_sources gpu_hello_rt_api.cu gpu_device_info.cu gpu_dual_call_test.cu)
  set(cuda_headers basic_kernels.hh Dual_Call.hh)
  set(cuda_sources basic_kernels.cu Dual_Call.cu)

  # ---------------------------------------------------------------------------------------------- #
  # Build Unit tests
  # ---------------------------------------------------------------------------------------------- #
  add_component_library(
    TARGET Lib_device_test
    TARGET_DEPS Lib_dsxx Lib_device
    LIBRARY_NAME device_test
    LIBRARY_TYPE STATIC
    LANGUAGE "${GPU_DBS_STRING}"
    SOURCES "${cuda_sources}"
    HEADERS "${cuda_headers}")

  # ---------------------------------------------------------------------------------------------- #
  # Register Unit tests
  # ---------------------------------------------------------------------------------------------- #
  set(test_deps Lib_device Lib_device_test)

  add_scalar_tests(SOURCES "${test_sources}" DEPS "${test_deps}" LANGUAGE ${GPU_DBS_STRING})

  set_tests_properties(
    device_gpu_hello_rt_api
    PROPERTIES
    PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
    FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED")

endif(HAVE_GPU AND USE_GPU)

# ------------------------------------------------------------------------------------------------ #
# End device/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
