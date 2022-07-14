# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   lapack_wrap/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for lapack_wrap.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
cmake_minimum_required(VERSION 3.9.0)
project(lapack_wrap CXX)

# ------------------------------------------------------------------------------------------------ #
# Skip this package if
#
# 1. LAPACK isn't found.
# 2. We are targeting an Xcode solution (necessary cmake logic is missing).
#
# ------------------------------------------------------------------------------------------------ #
if(HAVE_Fortran)

  if(NOT TARGET LAPACK::LAPACK OR ${CMAKE_GENERATOR} MATCHES Xcode)
    message(STATUS "NOTE: Skipping lapack_wrap component because LAPACK is not found (or maybe "
                   "you are targeting an Xcode solution).")
    return()
  endif()

  # ---------------------------------------------------------------------------------------------- #
  # Generate config.h (only occurs when cmake is run)
  # ---------------------------------------------------------------------------------------------- #

  # ---------------------------------------------------------------------------------------------- #
  # MSVC + MinGW's gfortran
  #
  # * FortranCInterface fails for this build mode because we don't have Fortran and C compilers
  #   defined at the same time.  Create the macros manually for this special case.
  if(MSVC)
    if(NOT EXISTS ${PROJECT_BINARY_DIR}/lapack_wrap/config.h)
      set(file_content
          "
#ifndef FC_HEADER_INCLUDED
#define FC_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define FC_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define FC_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define FC_MODULE(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/* Mangling for Fortran module symbols with underscores. */
#define FC_MODULE_(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

#endif
")
      file(WRITE ${PROJECT_BINARY_DIR}/lapack_wrap/config.h ${file_content})
    endif()

  else()

    # -------------------------------------------------------------------------------------------- #
    # Linux solution (main path)

    # This module automatically detects the API by which C and Fortran languages interact.
    #
    # * This fails on rznevada when mixing LLVM C/C++ with Cray Fortran because a subroutine in a
    #   Fortran module is represented as `mysubroutine_$mymodule_` but C doesn't allow a dollar
    #   symbol in the symbol name.  Ugh.  Let if fail for now.
    include(FortranCInterface)
    if(FortranCInterface_GLOBAL_FOUND AND FortranCInterface_MODULE_FOUND)
      fortrancinterface_header(
        lapack_wrap/config.h
        MACRO_NAMESPACE "FC_"
        SYMBOLS scopy
                dcopy
                saxpy
                daxpy
                sscal
                dscal
                sdot
                ddot
                snrm2
                dnrm2
        SYMBOL_NAMESPACE "FC_")
    endif()
  endif()

  # ---------------------------------------------------------------------------------------------- #
  # Omit this package if the mangling was not found by FortranCInterface.cmake
  if(EXISTS ${PROJECT_BINARY_DIR}/lapack_wrap/config.h)

    # -------------------------------------------------------------------------------------------- #
    # Source files
    # -------------------------------------------------------------------------------------------- #

    file(GLOB sources *.cc)
    file(GLOB headers *.hh)
    list(APPEND headers ${PROJECT_BINARY_DIR}/lapack_wrap/config.h)

    # -------------------------------------------------------------------------------------------- #
    # Installation instructions
    # -------------------------------------------------------------------------------------------- #

    install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/lapack_wrap)

    # -------------------------------------------------------------------------------------------- #
    # Unit tests
    # -------------------------------------------------------------------------------------------- #

    if(BUILD_TESTING)
      add_subdirectory(test)
    endif()

    # -------------------------------------------------------------------------------------------- #
    # Autodoc
    # -------------------------------------------------------------------------------------------- #

    process_autodoc_pages()

  else()
    message("   lapack_wrap - skipped because FortranCInterface failed.")
  endif()

endif(HAVE_Fortran)

# ------------------------------------------------------------------------------------------------ #
# End lapack_wrap/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
