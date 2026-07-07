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

# Set CUDA compute capabilities to use host architecture (read: must configure on GPU node!)
if(NOT DEFINED CUDA_ARCHITECTURES)
  # This value is automatically added to all CUDA targets. See cmake/component_macros.cmake
  set(CUDA_ARCHITECTURES
      "native"
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
         " -Wno-psabi -Wno-deprecated-gpu-targets -Wdeprecated-declarations -fpermissive")
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
