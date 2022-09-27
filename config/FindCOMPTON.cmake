# --------------------------------------------*-cmake-*------------------------------------------- #
# file   config/FindCompton.cmake
# author Kendra Keady <keadyk@lanl.gov>
# date   2017 February 28
# brief  Instructions for discovering the Compton vendor libraries.
# note   Copyright (C) 2017-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# .rst: FindCOMPTON
# ---------
#
# Find the CSK COMPTON includes and libraries.
#
# CSK COMPTON is a project dedicated to the calculation and interpolation of Compton Scattering
# Kernel (CSK, hence the project name) values for use in radiative transfer simulations.
# https://gitlab.lanl.gov/keadyk/CSK_generator
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If COMPTON is found, this module defines the following :prop_tgt:`IMPORTED` targets::
#
# COMPTON::compton         - The COMPTON library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
# * COMPTON_FOUND        - True if COMPTON found on the local system
# * COMPTON_INCLUDE_DIRS - Location of COMPTON header files.
# * COMPTON_LIBRARIES    - The COMPTON libraries.
# * COMPTON_VERSION      - The version of the discovered COMPTON install.
#
# Hints ^^^^^
#
# Set ``COMPTON_ROOT_DIR`` to a directory that contains a COMPTON installation.
#
# This script expects to find libraries at ``$COMPTON_ROOT_DIR/lib`` and the COMPTON headers at
# ``$COMPTON_ROOT_DIR/include``.  The library directory may optionally provide Release and Debug
# folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of COMPTON installation
# discovered.  These variables may optionally be set to help this module find the correct files::
#
# * COMPTON_LIBRARY        - Location of the COMPTON library.
# * COMPTON_LIBRARY_DEBUG  - Location of the debug COMPTON library (if any).
#
# ------------------------------------------------------------------------------------------------ #

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

# Use OpenMP version if OpenMP is available.
if(NOT OpenMP_FOUND)
  find_package(OpenMP QUIET)
endif()

# ================================================================================================ #
# If the user has provided ``COMPTON_ROOT_DIR``, use it!  Choose items found at this location over
# system locations.
if(EXISTS "$ENV{COMPTON_ROOT_DIR}")
  file(TO_CMAKE_PATH "$ENV{COMPTON_ROOT_DIR}" COMPTON_ROOT_DIR)
  set(COMPTON_ROOT_DIR
      "${COMPTON_ROOT_DIR}"
      CACHE PATH "Prefix for COMPTON installation.")
endif()

# ================================================================================================ #
# Set COMPTON_INCLUDE_DIRS and COMPTON_LIBRARIES. Try to find the libraries at $COMPTON_ROOT_DIR (if
# provided) or in standard system locations.  These find_library and find_path calls will prefer
# custom locations over standard locations (HINTS).  If the requested file is not found at the HINTS
# location, standard system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path(
  COMPTON_INCLUDE_DIR
  NAMES multigroup_lib_builder.hh
  HINTS ${COMPTON_ROOT_DIR}/include ${COMPTON_INCLUDEDIR}
  PATH_SUFFIXES Release Debug)

set(COMPTON_LIBRARY_NAME Lib_compton_omp;Lib_compton;compton_omp;compton)
find_library(
  COMPTON_LIBRARY
  NAMES ${COMPTON_LIBRARY_NAME}
  HINTS ${COMPTON_ROOT_DIR}/lib ${COMPTON_LIBDIR}
  PATH_SUFFIXES Release Debug)
# Do we also have debug versions?
find_library(
  COMPTON_LIBRARY_DEBUG
  NAMES ${COMPTON_LIBRARY_NAME}
  HINTS ${COMPTON_ROOT_DIR}/lib ${COMPTON_LIBDIR}
  PATH_SUFFIXES Debug)
set(COMPTON_INCLUDE_DIRS ${COMPTON_INCLUDE_DIR})
set(COMPTON_LIBRARIES ${COMPTON_LIBRARY})

# Try to find the version.
if(NOT COMPTON_VERSION)
  if(EXISTS "${COMPTON_INCLUDE_DIRS}/compton.h")
    file(STRINGS "${COMPTON_INCLUDE_DIRS}/compton.h" compton_h_major
         REGEX "define COMPTON_VER_MAJOR")
    file(STRINGS "${COMPTON_INCLUDE_DIRS}/compton.h" compton_h_minor
         REGEX "define COMPTON_VER_MINOR")
    file(STRINGS "${COMPTON_INCLUDE_DIRS}/compton.h" compton_h_subminor
         REGEX "define COMPTON_VER_SUBMINOR")
    string(REGEX REPLACE ".*([0-9]+)" "\\1" COMPTON_MAJOR ${compton_h_major})
    string(REGEX REPLACE ".*([0-9]+)" "\\1" COMPTON_MINOR ${compton_h_minor})
    string(REGEX REPLACE ".*([0-9]+)" "\\1" COMPTON_SUBMINOR ${compton_h_subminor})
    set(COMPTON_VERSION
        "${COMPTON_MAJOR}.${COMPTON_MINOR}.${COMPTON_SUBMINOR}"
        CACHE STRING "CSK version" FORCE)
  endif()
  # We might also try scraping the directory name for a regex match "csk-X.X.X"
  if(NOT COMPTON_VERSION)
    string(REGEX REPLACE ".*csk-([0-9]+).([0-9]+).([0-9]+).*" "\\1" COMPTON_MAJOR
                         ${COMPTON_INCLUDE_DIR})
    string(REGEX REPLACE ".*csk-([0-9]+).([0-9]+).([0-9]+).*" "\\2" COMPTON_MINOR
                         ${COMPTON_INCLUDE_DIR})
    string(REGEX REPLACE ".*csk-([0-9]+).([0-9]+).([0-9]+).*" "\\3" COMPTON_SUBMINOR
                         ${COMPTON_INCLUDE_DIR})
    set(COMPTON_VERSION
        "${COMPTON_MAJOR}.${COMPTON_MINOR}.${COMPTON_SUBMINOR}"
        CACHE STRING "CSK version" FORCE)
  endif()
endif()

# ================================================================================================ #
# handle the QUIETLY and REQUIRED arguments and set COMPTON_FOUND to TRUE if all listed variables
# are TRUE.
find_package_handle_standard_args(
  COMPTON
  FOUND_VAR COMPTON_FOUND
  REQUIRED_VARS COMPTON_INCLUDE_DIR COMPTON_LIBRARY
  VERSION_VAR COMPTON_VERSION)

mark_as_advanced(
  COMPTON_ROOT_DIR
  COMPTON_VERSION
  COMPTON_LIBRARY
  COMPTON_INCLUDE_DIR
  COMPTON_LIBRARY_DEBUG
  COMPTON_USE_PKGCONFIG
  COMPTON_CONFIG)

# ================================================================================================ #
# Register imported libraries:
#
# ~~~
# 1. If we can find a Windows .dll file (or if we can find both Debug and Release libraries), we
#    will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and include directory.
# ~~~

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string(REPLACE ".lib" ".dll" COMPTON_LIBRARY_DLL "${COMPTON_LIBRARY}")
  string(REPLACE ".lib" ".dll" COMPTON_LIBRARY_DEBUG_DLL "${COMPTON_LIBRARY_DEBUG}")
endif()

if(COMPTON_FOUND AND NOT TARGET COMPTON::compton)
  if(EXISTS "${COMPTON_LIBRARY_DLL}")

    # Windows systems with dll libraries.
    add_library(COMPTON::compton SHARED IMPORTED)

    # Windows with dlls, but only Release libraries.
    set_target_properties(
      COMPTON::compton
      PROPERTIES
      IMPORTED_LOCATION_RELEASE "${COMPTON_LIBRARY_DLL}"
      IMPORTED_IMPLIB "${COMPTON_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${COMPTON_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")

    # If we have both Debug and Release libraries
    if(EXISTS "${COMPTON_LIBRARY_DEBUG_DLL}")
      set_property(
        TARGET COMPTON::compton
        APPEND
        PROPERTY IMPORTED_CONFIGURATIONS Debug)
      set_target_properties(
        COMPTON::compton PROPERTIES
        IMPORTED_LOCATION_DEBUG "${COMPTON_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG "${COMPTON_LIBRARY_DEBUG}")
    endif()

  else()

    # For all other environments (ones without dll libraries), create the imported library targets.
    add_library(COMPTON::compton UNKNOWN IMPORTED)
    set_target_properties(
      COMPTON::compton PROPERTIES
      IMPORTED_LOCATION "${COMPTON_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${COMPTON_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")
  endif()
endif()

# ------------------------------------------------------------------------------------------------ #
# End FindCOMPTON.cmake
# ------------------------------------------------------------------------------------------------ #
