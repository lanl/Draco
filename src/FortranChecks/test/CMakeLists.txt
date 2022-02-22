#--------------------------------------------------------------------------------------------------#
# file   FortranChecks/test/CMakeLists.txt
# author Kelly Thompson
# brief  Demo of building a C++ main that links against an F90 library and uses ISO_C_BINDING.
# note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project(FortranChecks_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(test_sources ${PROJECT_SOURCE_DIR}/cppmain.cc ${PROJECT_SOURCE_DIR}/tstDerivedTypes.cc
                 ${PROJECT_SOURCE_DIR}/tstDrel.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #

# When building with MVSE and linking Fortran libraries generated with MinGW's gfortran, the
# libquadmath-0.dll library must be found in the path.
if(MSVC)
  get_filename_component(CAFS_Fortran_DIR "${CAFS_Fortran_COMPILER}" DIRECTORY)
  find_file(mingw_libquadmath libquadmath-0.dll HINTS ${CAFS_Fortran_DIR})
  if(NOT mingw_libquadmath)
    message(FATAL_ERROR "ERROR in FortranChecks\test:
      The library libquadmath-0.dll must be found in your PATH to allow these tests to run!")
  endif()
endif()

add_scalar_tests(SOURCES "${test_sources}" DEPS "Lib_FC_f90sub")

# ------------------------------------------------------------------------------------------------ #
# end of FortranChecks/test/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
