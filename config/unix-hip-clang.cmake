# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-hip-clang.cmake
# brief  Establish flags for Unix/Linux - HIP
# note   Copyright (C) 2023-2024 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

if(NOT CMAKE_HIP_COMPILER_ID STREQUAL Clang)
  message(FATAL_ERROR "Draco only supports the Clang hip compiler.")
endif()

#
# Compiler Flags
#

if(NOT HIP_FLAGS_INITIALIZED)

  set(HIP_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  # Discover HIP compute capabilities.
  if(NOT DEFINED HIP_ARCHITECTURES)
    if(DEFINED CMAKE_HIP_ARCHITECTURES)
      list(GET CMAKE_HIP_ARCHITECTURES 1 GFXMODEL)
    else()
      # This logic is from Bob Robey, see
      # https://github.com/amd/HPCTrainingExamples/blob/main/tests/CMakeLists.txt
      execute_process(
        COMMAND rocminfo
        COMMAND grep -m 1 gfx
        COMMAND sed "s/Name://"
        COMMAND head -1
        COMMAND sed "s/ //g"
        OUTPUT_VARIABLE GFXMODEL
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
    set(HIP_ARCHITECTURES "${GFXMODEL}")
    set(HIP_ARCHITECTURES
        ${HIP_ARCHITECTURES}
        CACHE STRING "target architecture for gpu code.")
  endif()

  # -O0 worked with cce-15/rocm-5.4.3, but is broken with cce-16/rocm-5.7.1. A fix is coming in
  # rocm-6.1. For now, build Debug HIP code at -O1.
  string(CONCAT CMAKE_HIP_FLAGS_DEBUG " -O1")
  # -fno-cray -mllvm -amdgpu-early-inline-all=false -mllvm -amdgpu-function-calls=true")

  # As of cce-16/rocm-5.7.1, everything is passing at -O3.
  set(CMAKE_HIP_FLAGS_RELEASE "-O3")
  set(CMAKE_HIP_FLAGS_MINSIZEREL "-O3")
  set(CMAKE_HIP_FLAGS_RELWITHDEBINFO "-O3")

  # Disable FMA at the compile level if desired
  if(DEFINED FMA_NEVER_HARDWARE)
    string(APPEND CMAKE_HIP_FLAGS " -ffp-contract=off")
  endif()

  # We are seeing very slow HIP compile times when using the defaults provided by ATS-4 machines.
  # During HIP compiler detection (/opt/rocm-5.7.1/lib/cmake/hip-lang/hip-lang-config.cmake), the
  # following property is set by cmake:
  #
  # ~~~
  # set_property(TARGET hip-lang::device APPEND PROPERTY INTERFACE_COMPILE_OPTIONS
  #     $<$<COMPILE_LANGUAGE:HIP>:SHELL:-mllvm;-amdgpu-early-inline-all=true;
  #     -mllvm;-amdgpu-function-calls=false> )
  # ~~~
  #
  # For better compiler performance we need to reverse the two boolean values (as per discussion
  # with lc-hotline, email to KT, 6/25/2024) MAC 8/5/2024 - Removed because switching this flag
  # causes unexplained segfault errors in CCE/18.0

  get_target_property(current_value hip-lang::device INTERFACE_COMPILE_OPTIONS)
  # if(current_value MATCHES "amdgpu-early-inline-all=true") string(REPLACE
  # "amdgpu-early-inline-all=true" "amdgpu-early-inline-all=false" current_value "${current_value}")
  # string(REPLACE "amdgpu-function-calls=false" "amdgpu-function-calls=true" current_value
  # "${current_value}") set_target_properties(hip-lang::device PROPERTIES INTERFACE_COMPILE_OPTIONS
  # "${current_value}") endif()

endif()

string(APPEND CMAKE_CXX_FLAGS " -D__HIP_PLATFORM_HCC__= -D__HIP_PLATFORM_AMD__=")
string(APPEND CMAKE_CXX_FLAGS " -Wno-switch-default")

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_HIP_FLAGS)
force_compiler_flags_to_cache("HIP")
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("CXX")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-cuda.cmake
# ------------------------------------------------------------------------------------------------ #
