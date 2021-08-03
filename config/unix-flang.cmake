# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-flang.cmake
# author Kelly Thompson
# date   Sunday, Apr 29, 2018, 19:56 pm
# brief  Establish flags for Unix/Linux - Gnu Fortran
# note   Copyright (C) 2018-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

if(NOT Fortran_FLAGS_INITIALIZED)
  mark_as_advanced(CMAKE_Fortran_COMPILER_VERSION)
endif()

# Compiler Flags

if(NOT Fortran_FLAGS_INITIALIZED)
  set(Fortran_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_Fortran_FLAGS " -g -cpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-gdwarf-3 -DDEBUG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-gdwarf-3 ${CMAKE_Fortran_FLAGS_RELEASE}")

endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_Fortran_FLAGS)
force_compiler_flags_to_cache("Fortran")

# Toggle compiler flags for optional features
if(OpenMP_Fortran_FLAGS)
  toggle_compiler_flag(OPENMP_FOUND ${OpenMP_Fortran_FLAGS} "Fortran" "")
endif()

# ------------------------------------------------------------------------------------------------ #
# End config/unix-flang.cmake
# ------------------------------------------------------------------------------------------------ #
