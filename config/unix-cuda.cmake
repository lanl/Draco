# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-cuda.cmake
# brief  Establish flags for Unix/Linux - Cuda
# note   Copyright (C) 2020-2026 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

# cmake-lint: disable=C0301,W0106

if(NOT CMAKE_CUDA_COMPILER_ID STREQUAL NVIDIA)
  message(FATAL_ERROR "Draco only supports the Nvidia cuda compiler.")
endif()

# Notes:
#
# Useful options that could be added to aid debugging
#
# * https://devblogs.nvidia.com/tag/cuda/
# * https://devblogs.nvidia.com/building-cuda-applications-cmake/

# Discover CUDA compute capabilities.
if(NOT DEFINED CUDA_ARCHITECTURES)
  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/config)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/config)
  endif()
  set(OUTPUTFILE ${CMAKE_CURRENT_BINARY_DIR}/config/cuda_script)
  set(CUDAFILE ${CMAKE_CURRENT_SOURCE_DIR}/config/query_gpu.cu)
  execute_process(COMMAND ${CMAKE_CUDA_COMPILER} -lcuda ${CUDAFILE} -o ${OUTPUTFILE})
  execute_process(
    COMMAND ${OUTPUTFILE}
    RESULT_VARIABLE CUDA_RETURN_CODE
    OUTPUT_VARIABLE CUDA_ARCHITECTURES)
  if(NOT ${CUDA_RETURN_CODE} EQUAL 0)
    message(FATAL_ERROR "Unable to determine target Cuda arch.")
  endif()
  unset(OUTPUTFILE)
  unset(CUDAFILE)
  unset(CUDA_RETURN_CODE)
  # This value is automatically added to all CUDA targets. See cmake/component_macros.cmake
  set(CUDA_ARCHITECTURES
      ${CUDA_ARCHITECTURES}
      CACHE STRING "target architecture for gpu code.")
endif()

#
# Compiler Flags
#

if(NOT CUDA_FLAGS_INITIALIZED)

  set(CUDA_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_CUDA_FLAGS " -g --expt-relaxed-constexpr --expt-extended-lambda -Xcompiler"
         " -Wno-psabi -Wno-deprecated-gpu-targets")
  string(CONCAT CMAKE_CUDA_FLAGS_DEBUG " -O0 -Xcudafe --display_error_number"
                " -Xcudafe --diag_suppress=1427")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O2") # -dipo
  set(CMAKE_CUDA_FLAGS_MINSIZEREL "-O2")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O2 --generate-line-info")

  # Disable FMA at the compile level if desired
  if(DEFINED FMA_NEVER_HARDWARE)
    string(APPEND CMAKE_CUDA_FLAGS " --fmad=false")
  endif()
endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_CUDA_FLAGS)
force_compiler_flags_to_cache("CUDA")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-cuda.cmake
# ------------------------------------------------------------------------------------------------ #
