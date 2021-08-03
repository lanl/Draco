# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-crayCC.cmake
# author Kelly Thompson
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Cray C/C++
# note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

#
# Compiler flag checks
#
if(CMAKE_CXX_COMPILER_VERSION LESS "8.4")
  message(FATAL_ERROR "Cray C++ prior to 8.4 does not support C++11.")
endif()

#
# Compiler Flags
#

# As of cmake-3.12.1, the cray flags for C++14 were not available in cmake.
if(NOT DEFINED CMAKE_CXX14_STANDARD_COMPILE_OPTION)
  set(CMAKE_CXX14_STANDARD_COMPILE_OPTION
      "-hstd=c++14"
      CACHE STRING "internal" FORCE)
endif()

if(NOT CXX_FLAGS_INITIALIZED)
  set(CXX_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_C_FLAGS " -g -DR123_USE_GNU_UINT128=0")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -DDEBUG")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -DNDEBUG")

  string(APPEND CMAKE_CXX_FLAGS " ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")

endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_C_FLAGS)
deduplicate_flags(CMAKE_CXX_FLAGS)
force_compiler_flags_to_cache("C;CXX")

toggle_compiler_flag(OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX" "")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-crayCC.cmake
# ------------------------------------------------------------------------------------------------ #
