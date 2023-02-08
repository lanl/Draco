# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-pgf90.cmake
# author Kelly Thompson
# date   2011 June 7
# brief  Establish flags for Unix - PGI Fortran
# note   Copyright (C) 2011-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

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
  string(APPEND CMAKE_Fortran_FLAGS " -Mpreprocess")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Mbounds -Mchkptr")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_DEBUG} -O3")
endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_Fortran_FLAGS)

# Toggle compiler flags for optional features
toggle_compiler_flag(OpenMP_FOUND "${OpenMP_Fortran_FLAGS}" "Fortran" "")
force_compiler_flags_to_cache("Fortran")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-pgf90.cmake
# ------------------------------------------------------------------------------------------------ #
