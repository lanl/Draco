# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/unix-xlf.cmake
# author Gabriel Rockefeller, Kelly Thompson <kgt@lanl.gov>
# date   2012 Nov 1
# brief  Establish flags for Unix - IBM XL Fortran
# note   Copyright (C) 2016-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)

# Let anyone who is interested in which FORTRAN compiler we're using switch on this macro.
set(CMAKE_Fortran_COMPILER_FLAVOR "XL")

# -qsuppress=
#
# * 1501-210 (W) command option t contains an incorrect subargument This is related to '-pthread'
#   showing up in the flags (a cmake bug)

if(NOT Fortran_FLAGS_INITIALIZED)
  set(Fortran_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")
  string(APPEND CMAKE_Fortran_FLAGS " -g -qlanglvl=2008std -qflag=i:w -qsuppress=1501-210")
  if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
    string(APPEND CMAKE_Fortran_FLAGS " -qarch=pwr9 -qtune=pwr9")
  endif()

  set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -qstrict=nans:operationprecision")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

# ------------------------------------------------------------------------------------------------ #
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_Fortran_FLAGS)
force_compiler_flags_to_cache("Fortran")

# ------------------------------------------------------------------------------------------------ #
# End config/unix-xlf.cmake
# ------------------------------------------------------------------------------------------------ #
