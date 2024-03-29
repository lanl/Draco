# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-hip-clang.cmake
# brief  Establish flags for Unix/Linux - HIP
# note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

if(NOT CMAKE_HIP_COMPILER_ID STREQUAL Clang)
  message(FATAL_ERROR "Draco only supports the Clang hip compiler.")
endif()

# Notes:
#
# * Options that control floating point behavior are listed at
#   https://clang.llvm.org/docs/UsersManual.html#controlling-floating-point-behavior

# Discover CUDA compute capabilities.
# ~~~
# if(NOT DEFINED CUDA_ARCHITECTURES)
#   if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/config)
#     file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/config)
#   endif()
#   set(OUTPUTFILE ${CMAKE_CURRENT_BINARY_DIR}/config/cuda_script)
#   set(CUDAFILE ${CMAKE_CURRENT_SOURCE_DIR}/config/query_gpu.cu)
#   execute_process(COMMAND ${CMAKE_CUDA_COMPILER} -lcuda ${CUDAFILE} -o ${OUTPUTFILE})
#   execute_process(COMMAND ${OUTPUTFILE}
#     RESULT_VARIABLE CUDA_RETURN_CODE
#     OUTPUT_VARIABLE CUDA_ARCHITECTURES)
#   if(NOT ${CUDA_RETURN_CODE} EQUAL 0)
#     message(FATAL_ERROR "Unable to determine target Cuda arch.")
#   endif()
#   unset(OUTPUTFILE)
#   unset(CUDAFILE)
#   unset(CUDA_RETURN_CODE) # This value is automatically added to all CUDA targets.
#   # See cmake/component_macros.cmake
#   set(CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES} CACHE STRING "target architecture for gpu code.")
# endif()
# ~~~

#
# Compiler Flags
#

if(NOT HIP_FLAGS_INITIALIZED)

  set(HIP_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  # Discover HIP compute capabilities.
  if(NOT DEFINED HIP_ARCHITECTURES)
    # see if we are on a 250 or a 300 node
    if(SITENAME MATCHES "RZVernal")
      set(HIP_ARCHITECTURES "gfx90a")
    elseif(SITENAME MATCHES "RZAdams" OR SITENAME MATCHES "ElCapitan")
      set(HIP_ARCHITECTURES "gfx942")
    else()
      message(FATAL_ERROR "unix-hip-clang.cmake: Logic error in setting HIP_ARCHITECTURES"
                          "Unknown platform. Variable not set.")
    endif()
    set(HIP_ARCHITECTURES
        ${HIP_ARCHITECTURES}
        CACHE STRING "target architecture for gpu code.")
  endif()

  # -O0 worked with cce-15/rocm-5.4.3, but is broken with cce-16/rocm-5.7.1. A fix is coming in
  # rocm-6.1. For now, build Debug HIP code at -O1.
  string(CONCAT CMAKE_HIP_FLAGS_DEBUG " -O1")

  # As of cce-16/rocm-5.7.1, everything is passing at -O3.
  set(CMAKE_HIP_FLAGS_RELEASE "-O3")
  set(CMAKE_HIP_FLAGS_MINSIZEREL "-O3")
  set(CMAKE_HIP_FLAGS_RELWITHDEBINFO "-O3")

  # Disable FMA at the compile level if desired
  if(DEFINED FMA_NEVER_HARDWARE)
    string(APPEND CMAKE_HIP_FLAGS " -ffp-contract=off")
  endif()

endif()

string(APPEND CMAKE_CXX_FLAGS " -D__HIP_PLATFORM_HCC__= -D__HIP_PLATFORM_AMD__=")

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_HIP_FLAGS)
force_compiler_flags_to_cache("HIP")
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("CXX")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-cuda.cmake
# ------------------------------------------------------------------------------------------------ #
