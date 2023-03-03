# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/Linux-NVHPC-Fortran.cmake
# brief  Establish Fortran flags for Linux - NVHPC
# note   Copyright (C) 2022-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

#
# Compiler Flags
#
if(NOT Fortran_FLAGS_INITIALIZED)
  set(CMAKE_Fortran_COMPILER_VERSION
      ${CMAKE_Fortran_COMPILER_VERSION}
      CACHE STRING "Fortran compiler version string" FORCE)
  mark_as_advanced(CMAKE_Fortran_COMPILER_VERSION)
  set(Fortran_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")
  string(APPEND CMAKE_Fortran_FLAGS " -Mpreprocess -Mint128")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Mbounds -Mchkptr")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_DEBUG} -O3")
endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_Fortran_FLAGS)

# Toggle compiler flags for optional features
force_compiler_flags_to_cache("Fortran")

# ------------------------------------------------------------------------------------------------ #
# End config/Linux-NVHPC-Fortran.cmake
# ------------------------------------------------------------------------------------------------ #
