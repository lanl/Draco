# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/dracoTesting.cmake
# brief  Enable ctest and point to the resource specification file.
# note   Copyright (C) 2011-2024 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

include(FeatureSummary)
option(BUILD_TESTING "Should we compile the tests?" ON)
add_feature_info(BUILD_TESTING BUILD_TESTING
                 "Turn off to prevent the compilation of unit tests (ctest).")

# enable ctest funcitons and run ctest in parallel if we have multiple cores.
if(BUILD_TESTING)
  include(CTest)
  enable_testing()
  set(CTEST_RESOURCE_SPEC_FILE
      "${PROJECT_BINARY_DIR}/ctest_resource_spec.json"
      CACHE FILEPATH "json file that details available resources.")
endif()

# ------------------------------------------------------------------------------------------------ #
# Copy the appropriate json file to the binary. This file provides details about platform resources
# that can be used when running ctest (eg. number of gpus).
#
# * Resources can be oversubscribed by declaring more available resources than actual.  Eg. You can
#   declare that there are 8 GPUs, when there is only 1. This would allow 8 tests to share the GPU
#   when running ctest.
function(provide_ctest_resources_spec)

  # These files will be at $draco/config or at $draco_install/cmake.
  if("${PROJECT_NAME}" STREQUAL "Draco")
    set(json_dir "${PROJECT_SOURCE_DIR}/config")
  elseif(EXISTS "${DRACO_DIR}")
    set(json_dir "${DRACO_DIR}/cmake") # we are building a client that uses the draco build sys.
  elseif(EXISTS "${draco_DIR}")
    set(json_dir "${draco_DIR}")
  endif()

  # Look for a platform specific resource file.
  if(DEFINED SITENAME_FAMILY)
    set(json_file "${json_dir}/ctest_resource_spec-${SITENAME_FAMILY}.json")
  endif()

  # If the platform file is not found, revert to the generic one.
  if(NOT EXISTS "${json_file}" OR (NOT DEFINED SITENAME_FAMILY))
    # * Consider creating this file dynamically during configure.  Can we autodetect number of GPUs
    #   or amount of memory?
    set(json_file "${json_dir}/ctest_resource_spec.json")
  endif()

  # Sanity check
  if(NOT EXISTS "${json_file}")
    message(FATAL_ERROR "Unable to locate required file ${json_file}.\n"
                        "You can disable testing or provide this file manually")
  endif()

  # Copy the specification file to the build directory so it can be used for scheduling test
  # execution via ctest.
  configure_file("${json_file}" "${PROJECT_BINARY_DIR}/ctest_resource_spec.json" COPYONLY)

endfunction(provide_ctest_resources_spec)

# ------------------------------------------------------------------------------------------------ #
# End dracoTesting.cmake
# ------------------------------------------------------------------------------------------------ #
