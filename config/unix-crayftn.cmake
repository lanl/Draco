# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-crayftn.cmake
# author Kelly Thompson
# date   2008 May 30
# brief  Establish flags for Unix - Cray Fortran
# note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

#
# Compiler flags:
#
if(NOT Fortran_FLAGS_INITIALIZED)
  set(Fortran_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  # Find and cache the compiler version (2015-09-28 CMake-3.3.1 misses this).
  execute_process(
    COMMAND ${CMAKE_Fortran_COMPILER} -V
    ERROR_VARIABLE ftn_version_output
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE ".*Version ([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\1.\\2"
                       CMAKE_Fortran_COMPILER_VERSION "${ftn_version_output}")
  set(CMAKE_Fortran_COMPILER_VERSION
      ${CMAKE_Fortran_COMPILER_VERSION}
      CACHE STRING "Fortran compiler version string" FORCE)
  mark_as_advanced(CMAKE_Fortran_COMPILER_VERSION)

  string(APPEND CMAKE_Fortran_FLAGS " -g")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -DDEBUG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -DDEBUG")

endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_Fortran_FLAGS)
force_compiler_flags_to_cache("Fortran")

toggle_compiler_flag(OPENMP_FOUND ${OpenMP_Fortran_FLAGS} "Fortran" "")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-crayftn.cmake
# ------------------------------------------------------------------------------------------------ #
