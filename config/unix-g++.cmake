# --------------------------------------------*-cmake-*------------------------------------------- #
# file   config/unix-g++.cmake
# brief  Establish flags for Unix/Linux - Gnu C++
# note   Copyright (C) 2010-2026 Triad National Security, LLC., All rights reserved.
#
# Best Practices for C++:
#
# 1. Turn on modern C++ and strict diagnostics:
#    * -std=c++20 -Wall -Wextra -Wpedantic -Werror
# 2. Add high-value C++ warnings (low noise, good catches)
#    * Using:
#      * -Wconversion -Wshadow -Wold-style-cast -Wuseless-cast -Woverloaded-virtual
#        -Wsuggest-override -Wformat=2 -Wcast-align=strict -Wduplicated-branches -Wnon-virtual-dtor
#        -Wduplicated-cond
#    * Consider adding (#1508):
#      * -Wdelete-non-virtual-dtor -Wzero-as-null-pointer-constant -Wcast-qual -Wnull-dereference
# 3. Enable GCC’s static analyzer
#    * -fanalyzer -Wanalyzer-too-complex
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

if(DEFINED CMAKE_CXX_COMPILER_VERSION
   AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0
   AND NOT (MSVC OR MINGW))
  message(FATAL_ERROR "Draco requires GNU gcc/g++ version >= 9.0.")
endif()

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
check_c_compiler_flag("-march=native" HAS_MARCH_NATIVE)
# Figure out what gcc does when give '-march=native': 'gcc -### -march=native -x c -'

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
         " -g -Wall -Wextra -Wpedantic -Wcast-align -Wconversion -Wpointer-arith -Wfloat-equal"
         " -Wunused-macros -Wshadow -Wformat=2")
  # Pragmas in rng and Random123 don't seem to disable this check.
  string(APPEND CMAKE_C_FLAGS " -Wno-expansion-to-defined")
  string(CONCAT CMAKE_C_FLAGS_DEBUG "-fno-inline -fno-eliminate-unused-debug-types -O0"
                " -Wundef -Wunreachable-code -fdiagnostics-color=auto"
                " -Wnull-dereference -DDEBUG")
  # GCC_COLORS="error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01"
  string(CONCAT CMAKE_C_FLAGS_RELEASE "-O3 -funroll-loops -D_FORTIFY_SOURCE=2"
                " -fno-var-tracking-assignments -Wno-stringop-overflow -DNDEBUG")
  set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")
  string(CONCAT CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -fno-eliminate-unused-debug-types"
                " -funroll-loops -Wno-stringop-overflow")

  # Enable additional sanitizer flags for non-Cray (i.e. ccs-net) gcc builds. These can cause
  # link-time issues for Cray mixed builds (i.e. El Cap).
  if(NOT DEFINED ENV{TRAVIS} AND NOT CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv)
    # Some options (including these) seem to confuse Travis. See
    # https://stackoverflow.com/questions/50024731/ld-unrecognized-option-push-state-no-as-needed
    string(APPEND CMAKE_C_FLAGS_DEBUG " -fsanitize=bounds-strict")
    string(APPEND CMAKE_C_FLAGS_DEBUG " -fsanitize=float-divide-by-zero")
    string(APPEND CMAKE_C_FLAGS_DEBUG " -fsanitize=float-cast-overflow")
    string(APPEND CMAKE_C_FLAGS_DEBUG " -fsanitize=signed-integer-overflow")
    string(APPEND CMAKE_C_FLAGS_DEBUG " -Wduplicated-branches")
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

  # Disable FMA at the compile level if desired
  if(DEFINED FMA_NEVER_HARDWARE)
    string(APPEND CMAKE_C_FLAGS " -ffp-contract=off")
  endif()

  string(APPEND CMAKE_CXX_FLAGS " ${CMAKE_C_FLAGS} -Wno-psabi")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Woverloaded-virtual -Wsuggest-override")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")

  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Wold-style-cast -Wnon-virtual-dtor -Wuseless-cast")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fdiagnostics-show-template-tree -Wnoexcept")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Wsuggest-attribute=const -Wduplicated-cond -Wlogical-op")

endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("C;CXX")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-g++.cmake
# ------------------------------------------------------------------------------------------------ #
