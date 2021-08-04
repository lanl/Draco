#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/unix-ifx.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Wednesday, May 26, 2021, 14:12 pm
# brief  Establish flags for Unix - Intel OneAPI Fortran (ifx)
# note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

#
# Compiler flags:
#
if(NOT Fortran_FLAGS_INITIALIZED)
  set(Fortran_FLAGS_INITIALIZED
      "yes"
      CACHE INTERNAL "using draco settings.")
  set(CMAKE_Fortran_COMPILER_VERSION
      ${CMAKE_Fortran_COMPILER_VERSION}
      CACHE STRING "Fortran compiler version string" FORCE)
  mark_as_advanced(CMAKE_Fortran_COMPILER_VERSION)

  string(APPEND CMAKE_Fortran_FLAGS " -g")
  string(CONCAT CMAKE_Fortran_FLAGS_DEBUG "-O0 -DDEBUG")
  if(NOT DEFINED CMAKE_CXX_COMPILER_WRAPPER AND NOT "${CMAKE_CXX_COMPILER_WRAPPER}" STREQUAL
                                                "CrayPrgEnv")
    string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -traceback")
  endif()
  string(CONCAT CMAKE_Fortran_FLAGS_RELEASE "-O2 -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -DDEBUG")

endif()

# remove -Werror if it exists; maybe use -diag-enable=error (syntax?)
if("${CMAKE_Fortran_FLAGS}" MATCHES "Werror")
  message(STATUS "Fortran option -Werror not supported by ifx. Removing it.\n"
                 "CMAKE_Fortran_FLAGS  = ${CMAKE_Fortran_FLAGS}")
  string(REPLACE "-Werror" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  message(STATUS "CMAKE_Fortran_FLAGS  = ${CMAKE_Fortran_FLAGS}")
endif()

# --------------------------------------------------------------------------------------------------#
# Ensure cache values always match current selection
deduplicate_flags(CMAKE_Fortran_FLAGS)
force_compiler_flags_to_cache("Fortran")

# Optional compiler flags
if(NOT ${SITENAME} STREQUAL "Trinitite")
  toggle_compiler_flag(ENABLE_SSE "-mia32 -axSSSE3" "Fortran" "") # sse3, ssse3
endif()
toggle_compiler_flag(OPENMP_FOUND ${OpenMP_Fortran_FLAGS} "Fortran" "")

# --------------------------------------------------------------------------------------------------#
# End config/unix-ifx.cmake
# --------------------------------------------------------------------------------------------------#
