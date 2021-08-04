# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-cuda.cmake
# brief  Establish flags for Unix/Linux - Cuda
# note   Copyright (C) 2020-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

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
# -Xcudafe --diag_suppress=1427: suppress "warning: offsetof applied to non-POD"

if(NOT CUDA_FLAGS_INITIALIZED)

  set(CUDA_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_CUDA_FLAGS " -G --expt-relaxed-constexpr")
  string(APPEND CMAKE_CUDA_FLAGS " --expt-extended-lambda")
  if(CMAKE_CXX_COMPILER_ID MATCHES "XL")
    string(APPEND CMAKE_CUDA_FLAGS " -DCUB_IGNORE_DEPRECATED_CPP_DIALECT"
           " -DTHRUST_IGNORE_DEPRECATED_CPP_DIALECT")
    string(APPEND CMAKE_CUDA_FLAGS " -ccbin ${CMAKE_CXX_COMPILER} -Xcompiler -std=c++14")
    if(EXISTS /usr/gapps)
      # ATS-2
      string(APPEND CMAKE_CUDA_FLAGS " -Xcompiler --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1"
             " -Xcompiler -qxflag=disable__cplusplusOverride")
    elseif(EXISTS ${CMAKE_CXX_COMPILER_CONFIG_FILE})
      # Darwin
      string(APPEND CMAKE_CUDA_FLAGS " -Xcompiler -F${CMAKE_CXX_COMPILER_CONFIG_FILE}")
    endif()
  endif()
  string(CONCAT CMAKE_CUDA_FLAGS_DEBUG " -O0 -Xcudafe --display_error_number"
                " -Xcudafe --diag_suppress=1427")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O2") # -dipo
  set(CMAKE_CUDA_FLAGS_MINSIZEREL "-O2")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O2 --generate-line-info")

endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_CUDA_FLAGS)
force_compiler_flags_to_cache("CUDA")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-cuda.cmake
# ------------------------------------------------------------------------------------------------ #
