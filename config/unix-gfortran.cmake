# -----------------------------*-cmake-*---------------------------------------------------------- #
# file   config/unix-gfortran.cmake
# author Kelly Thompson
# date   2010 Sep 27
# brief  Establish flags for Unix/Linux - Gnu Fortran
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

#
# Compiler Flags
#

if(NOT Fortran_FLAGS_INITIALIZED)
  set(Fortran_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")

  string(APPEND CMAKE_Fortran_FLAGS " -g -ffree-line-length-none -cpp")
  string(CONCAT CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -frange-check"
                " -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-character=127 -DDEBUG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mtune=native -ftree-vectorize -funroll-loops -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE}")

  # Only add '-march=native' if -march and -mtune are not already requested. Spack will add -march=
  # options if spacific target is requested (ie. target Haswell CPU when building on Skylake).
  if(NOT
     ("${FC}" MATCHES "-mtune"
      OR "${FC}" MATCHES "-march"
      OR "${CMAKE_Fortran_FLAGS}" MATCHES "-mtune"
      OR "${CMAKE_Fortran_FLAGS}" MATCHES "-march"))
    if(NOT APPLE AND HAS_MARCH_NATIVE)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
    endif()
  endif()
endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_Fortran_FLAGS)

# Toggle compiler flags for optional features
if(OpenMP_Fortran_FLAGS)
  toggle_compiler_flag(OpenMP_FOUND "${OpenMP_Fortran_FLAGS}" "Fortran" "")
endif()
force_compiler_flags_to_cache("Fortran")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-gfortran.cmake
# ------------------------------------------------------------------------------------------------ #
