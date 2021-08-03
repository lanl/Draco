#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/unix-g++.cmake
# brief  Establish flags for Unix/Linux - Gnu C++
# note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

include_guard(GLOBAL)

# Note: In config/compilerEnv.cmake, the build system sets flags for
#
# * the language standard (C++14, C99, etc)
# * interprocedural optimization.

# Notes:
# ----------------------------------------
# Useful options that could be added to aid debugging
#
# * http://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
# * https://github.com/lefticus/cppbestpractices/blob/master/02-Use_the_Tools_Available.md#compilers

#
# Compiler flag checks
#
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
# arm --> -mcpu=thunderx2t99 ??? Darwin volta nodes?
check_c_compiler_flag("-march=native" HAS_MARCH_NATIVE)

#
# Compiler Flags

# * Consider using these optimization flags: -ffast-math -ftree-vectorize -fno-finite-math-only
#   -fno-associative-math -fsignaling-nans
# * Control FMA -ffp-contract=off

if(NOT CXX_FLAGS_INITIALIZED)
  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")
  string(APPEND CMAKE_C_FLAGS
         " -g -Wall -Wextra -pedantic -Wcast-align -Wpointer-arith -Wfloat-equal -Wunused-macros"
         " -Wshadow -Wformat=2")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7.0)
    string(APPEND CMAKE_C_FLAGS " -Wnull-dereference")
  endif()
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    # Pragmas in rng and Random123 don't seem to disable this check.
    string(APPEND CMAKE_C_FLAGS " -Wno-expansion-to-defined")
  endif()
  string(
    CONCAT CMAKE_C_FLAGS_DEBUG "-fno-inline -fno-eliminate-unused-debug-types -O0"
           " -Wundef -Wunreachable-code -fsanitize=bounds-strict -fdiagnostics-color=auto -DDEBUG")
  # GCC_COLORS="error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01"
  string(CONCAT CMAKE_C_FLAGS_RELEASE "-O3 -funroll-loops -D_FORTIFY_SOURCE=2"
                " -fno-var-tracking-assignments -DNDEBUG")
  set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")
  string(CONCAT CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -fno-eliminate-unused-debug-types"
                " -funroll-loops")

  if(NOT DEFINED ENV{TRAVIS})
    # Some options (including these) seem to confuse Travis. See
    # https://stackoverflow.com/questions/50024731/ld-unrecognized-option-push-state-no-as-needed
    string(APPEND CMAKE_C_FLAGS_DEBUG " -fsanitize=float-divide-by-zero")
    string(APPEND CMAKE_C_FLAGS_DEBUG " -fsanitize=float-cast-overflow")
  endif()

  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7.0)
    # See https://gcc.gnu.org/gcc-7/changes.html
    if(NOT DEFINED ENV{TRAVIS})
      string(APPEND CMAKE_C_FLAGS_DEBUG " -fsanitize=signed-integer-overflow")
      string(APPEND CMAKE_C_FLAGS_DEBUG " -Wduplicated-branches")
    endif()
  endif()

  # Systems running CrayPE use compile wrappers to specify this option.
  site_name(sitename)
  string(REGEX REPLACE "([A-z0-9]+).*" "\\1" sitename ${sitename})
  # Only add '-march=native' if -march and -mtune are not already requested. Spack will add -march=
  # options if spacific target is requested (ie. target Haswell CPU when building on Skylake).
  if(NOT
     ("${CC}" MATCHES "-mtune"
      OR "${CC}" MATCHES "-march"
      OR "${CMAKE_C_FLAGS}" MATCHES "-mtune"
      OR "${CMAKE_C_FLAGS}" MATCHES "-march"))
    if(HAS_MARCH_NATIVE
       AND NOT APPLE
       AND NOT CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv)
      string(APPEND CMAKE_C_FLAGS " -march=native")
    elseif(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
      string(APPEND CMAKE_C_FLAGS " -mcpu=powerpc64le -mtune=powerpc64le")
    endif()
  endif()

  string(APPEND CMAKE_CXX_FLAGS " ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Woverloaded-virtual")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")

  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0)
    # See https://gcc.gnu.org/gcc-8/changes.html
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Wold-style-cast -Wnon-virtual-dtor -Wuseless-cast")
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fdiagnostics-show-template-tree -Wnoexcept")
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Wsuggest-attribute=const")
  endif()

endif()

# --------------------------------------------------------------------------------------------------#
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("C;CXX")

# Toggle compiler flags for optional features
toggle_compiler_flag(OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX" "")

# --------------------------------------------------------------------------------------------------#
# End config/unix-g++.cmake
# --------------------------------------------------------------------------------------------------#
