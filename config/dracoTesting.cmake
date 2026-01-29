# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/dracoTesting.cmake
# brief  Enable ctest and point to the resource specification file.
# note   Copyright (C) 2011-2026 Triad National Security, LLC., All rights reserved.
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
# Convert <value><unit> to GB (integer, rounding down) -
#
# * out_var: variable name to receive GB
# * value: integer value
# * unit: one of K/M/G/T (case-insensitive)
function(mem_to_gb out_var value unit)
  if(unit MATCHES "^[Gg]$")
    set(_gb "${value}")
  elseif(unit MATCHES "^[Mm]$")
    math(EXPR _gb "${value} / 1024")
  elseif(unit MATCHES "^[Kk]$")
    math(EXPR _gb "${value} / (1024*1024)")
  elseif(unit MATCHES "^[Tt]$")
    math(EXPR _gb "${value} * 1024")
  else()
    set(_gb "")
  endif()

  set(${out_var}
      "${_gb}"
      PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------------------------ #
# On systems that run SLURM, this function will query the system to determine how many GPUs are
# available and how much memory is available on each GPU.  This information is required to correctly
# construct ctest_resource_spec.json which is read by ctest to specify proper resource management
# when running tests in parallel.

# cmake-lint: disable=R0915
function(slurm_extract_gpu_info out_gpu_count out_gpu_mem_gb)
  # Default values
  set(_gpu_count "0")
  set(_gpu_mem_gb "0")
  set(${out_gpu_count}
      "${_gpu_count}"
      PARENT_SCOPE)
  set(${out_gpu_mem_gb}
      "${_gpu_mem_gb}"
      PARENT_SCOPE)

  # Ensure SLURMD_NODENAME exists
  set(_nodename "$ENV{SLURMD_NODENAME}")
  if(_nodename STREQUAL "")
    return()
  endif()

  # Run scontrol
  execute_process(
    COMMAND scontrol show node "${_nodename}"
    RESULT_VARIABLE _rc
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
  if(NOT _rc EQUAL 0)
    return()
  endif()

  # scontrol often prints key=val pairs across multiple lines; normalize to a single line
  string(REPLACE "\n" " " _out_one_line "${_out}")

  # Extract the AvailableFeatures=... value (up to next whitespace) Example:
  # AvailableFeatures=...,gpu1_model:Tesla_V100S-PCIE-32GB,gpu1_memory:32GB,gpu_count:1... Note: If
  # SLURM Gres is configured properly, that output can also be used to get GPU information (See
  # original logic in Draco PR !620). However, Darwin Gres config is broken as of 1/2026 so we have
  # removed the Gres-specific logic.
  string(REGEX MATCH "([ \t]|^)AvailableFeatures=([^ \t]+)" _m "${_out_one_line}")
  if(_m STREQUAL "")
    return()
  endif()
  set(_availablefeatures "${CMAKE_MATCH_2}")

  # ---- GPU COUNT ----
  # Match: - gpu_count:4
  string(REGEX MATCH "(^|,)(gpu_count(:[^,:]+)?:)([0-9]+)" _mgpu "${_availablefeatures}")
  if(_mgpu STREQUAL "")
    message(STATUS "No GPU info found in AvailableFeatures='${_availablefeatures}'."
                   "Expected something like gpu_count:4.")
    return()
  endif()
  set(_gpu_count "${CMAKE_MATCH_4}")

  # ---- GPU MEM PER GPU (GB) ----
  set(_gpu_mem_gb "")

  # Get memory from (_availablefeatures): gpu1_memory::<num><unit>
  if(_gpu_mem_gb STREQUAL "")
    # unit is one of K/M/G/T (case-insensitive)
    string(REGEX MATCH "(^|,)(gpu1_memory(:[^,:]+)?:)([0-9]+)([KkMmGgTt])" _mmem
                 "${_availablefeatures}")
    if(NOT _mmem STREQUAL "")
      set(_mem_val "${CMAKE_MATCH_4}")
      set(_mem_unit "${CMAKE_MATCH_5}")
      mem_to_gb(_gpu_mem_gb "${_mem_val}" "${_mem_unit}")
    endif()
  endif()

  if(_gpu_mem_gb STREQUAL "")
    set(_gpu_mem_gb "16")
  endif()

  # Return results
  set(${out_gpu_count}
      "${_gpu_count}"
      PARENT_SCOPE)
  set(${out_gpu_mem_gb}
      "${_gpu_mem_gb}"
      PARENT_SCOPE)

  message(STATUS "Parsed SLURM node info: GPUs=${_gpu_count}, GPU_MEM_GB=${_gpu_mem_gb}")
endfunction()

# ------------------------------------------------------------------------------------------------ #
# Copy the appropriate json file to the binary. This file provides details about platform resources
# that can be used when running ctest (eg. number of gpus).
#
# * Resources can be oversubscribed by declaring more available resources than actual.  Eg. You can
#   declare that there are 8 GPUs, when there is only 1. This would allow 8 tests to share the GPU
#   when running ctest.
# * With minor effort, we could also pull actual node memory.
# * https://cmake.org/cmake/help/latest/manual/ctest.1.html#ctest-resource-specification-file
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
  set(json_file "${json_dir}/ctest_resource_spec-${SITENAME_FAMILY}.json")

  # If the platform file is not found, create a generic one from scratch.
  if(NOT EXISTS "${json_file}" OR (NOT DEFINED SITENAME_FAMILY))
    set(json_file "${PROJECT_BINARY_DIR}/ctest_resource_spec.json")
    # fetch number of gpus and gpu memory by querying slurm.
    slurm_extract_gpu_info(SLURM_GPU_COUNT SLURM_GPU_MEM_GB)

    set(N "${SLURM_GPU_COUNT}")

    # Build JSON fragments for the GPU arrays
    set(_gpus_json "")
    set(_gpu_mem_json "")

    if(N GREATER 0)
      math(EXPR _last "${N} - 1")
      foreach(i RANGE 0 ${_last} 1)
        # Commas between elements (but not after last)
        set(_comma "")
        if(i LESS _last)
          set(_comma ",")
        endif()

        string(APPEND _gpus_json "      {\n" "        \"id\": \"${i}\",\n" "        \"slots\": 4\n"
               "      }${_comma}\n")

        string(APPEND _gpu_mem_json "      {\n" "        \"id\": \"${i}\",\n"
               "        \"slots\": ${SLURM_GPU_MEM_GB}\n" "      }${_comma}\n")
      endforeach()
    endif()

    # Write the final JSON file
    file(
      WRITE "${json_file}"
      "{\n"
      "  \"version\": {\n"
      "    \"major\": 1,\n"
      "    \"minor\": 0\n"
      "  },\n"
      "  \"local\": [\n"
      "    {\n"
      "      \"gpus\": [\n"
      "${_gpus_json}"
      "      ],\n"
      "      \"gpu_mem_gb\": [\n"
      "${_gpu_mem_json}"
      "      ],\n"
      "      \"mem_gb\": [\n"
      "        {\n"
      "          \"id\": \"0\",\n"
      "          \"slots\": 256\n"
      "        }\n"
      "      ]\n"
      "    }\n"
      "  ]\n"
      "}\n")
  endif()

  # Sanity check
  if(NOT EXISTS "${json_file}")
    message(FATAL_ERROR "Unable to locate required file ${json_file}.\n"
                        "You can disable testing or provide this file manually")
  endif()

  # Copy the specification file to the build directory so it can be used for scheduling test
  # execution via ctest.
  configure_file("${json_file}" "${PROJECT_BINARY_DIR}/ctest_resource_spec.json" COPYONLY)
  message(STATUS "CTest resource specifications written to "
                 "${PROJECT_BINARY_DIR}/ctest_resource_spec.json")

endfunction(provide_ctest_resources_spec)

# ------------------------------------------------------------------------------------------------ #
# End dracoTesting.cmake
# ------------------------------------------------------------------------------------------------ #
