#--------------------------------------------*-cmake-*---------------------------------------------#
# file   draco/config/FindMETIS.cmake
# author Kelly Thompson <kgt@lanl.gov>
# brief  Helper for finding metis.
# note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# ~~~
# .rst:
# FindMETIS
# ---------
#
# Find the METIS includes and libraries.
#
# METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes,
# and producing fill reducing orderings for sparse matrices. The algorithms implemented in METIS
# are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint
# partitioning schemes developed in our lab.
# http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If METIS is found, this module defines the following :prop_tgt:`IMPORTED` targets::
#
#  METIS::metis         - The METIS library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  METIS_FOUND          - True if METIS found on the local system
#  METIS_INCLUDE_DIRS   - Location of METIS header files.
#  METIS_LIBRARIES      - The METIS libraries.
#  METIS_VERSION        - The version of the discovered METIS install.
#
# Hints
# ^^^^^
#
# Set ``METIS_ROOT_DIR`` to a directory that contains a METIS installation.
#
# This script expects to find libraries at ``$METIS_ROOT_DIR/lib`` and the METIS headers at
# ``$METIS_ROOT_DIR/include``.  The library directory may optionally provide Release and Debug
# folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of METIS installation
# discovered.  These variables may optionally be set to help this module find the correct files::
#
#  METIS_LIBRARY        - Location of the METIS library.
#  METIS_LIBRARY_DEBUG  - Location of the debug METIS library (if any).
#
# ~~~

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

# ================================================================================================ #
# If the user has provided ``METIS_ROOT_DIR``, use it!  Choose items found at this location over
# system locations.
if(EXISTS "$ENV{METIS_ROOT_DIR}")
  file(TO_CMAKE_PATH "$ENV{METIS_ROOT_DIR}" METIS_ROOT_DIR)
  set(METIS_ROOT_DIR
      "${METIS_ROOT_DIR}"
      CACHE PATH "Prefix for METIS installation.")
endif()

# ================================================================================================ #
# Set METIS_INCLUDE_DIRS and METIS_LIBRARIES. Try to find the libraries at $METIS_ROOT_DIR (if
# provided) or in standard system locations.  These find_library and find_path calls will prefer
# custom locations over standard locations (HINTS).  If the requested file is not found at the HINTS
# location, standard system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path(
  METIS_INCLUDE_DIR
  NAMES metis.h
  HINTS ${METIS_ROOT_DIR}/include ${METIS_INCLUDEDIR}
  PATH_SUFFIXES Release Debug)
find_library(
  METIS_LIBRARY
  NAMES metis
  HINTS ${METIS_ROOT_DIR}/lib ${METIS_LIBDIR}
  # PATH_SUFFIXES Release Debug
)
# Do we also have debug versions?
find_library(
  METIS_LIBRARY_DEBUG
  NAMES metis
  HINTS ${METIS_ROOT_DIR}/lib ${METIS_LIBDIR}
  PATH_SUFFIXES Debug)
set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
set(METIS_LIBRARIES ${METIS_LIBRARY})

# Try to find the version.
if(NOT METIS_VERSION)
  if(EXISTS "${METIS_INCLUDE_DIRS}/metis.h")
    file(STRINGS "${METIS_INCLUDE_DIRS}/metis.h" metis_h_major REGEX "define METIS_VER_MAJOR")
    file(STRINGS "${METIS_INCLUDE_DIRS}/metis.h" metis_h_minor REGEX "define METIS_VER_MINOR")
    file(STRINGS "${METIS_INCLUDE_DIRS}/metis.h" metis_h_subminor REGEX "define METIS_VER_SUBMINOR")
    string(REGEX REPLACE ".*([0-9]+)" "\\1" METIS_MAJOR ${metis_h_major})
    string(REGEX REPLACE ".*([0-9]+)" "\\1" METIS_MINOR ${metis_h_minor})
    string(REGEX REPLACE ".*([0-9]+)" "\\1" METIS_SUBMINOR ${metis_h_subminor})
  endif()
  # We might also try scraping the directory name for a regex match "metis-X.X.X"
endif()

# ================================================================================================ #
# handle the QUIETLY and REQUIRED arguments and set METIS_FOUND to TRUE if all listed variables are
# TRUE.
find_package_handle_standard_args(
  METIS
  FOUND_VAR METIS_FOUND
  REQUIRED_VARS METIS_INCLUDE_DIR METIS_LIBRARY
  VERSION_VAR METIS_VERSION)

mark_as_advanced(
  METIS_ROOT_DIR
  METIS_VERSION
  METIS_LIBRARY
  METIS_INCLUDE_DIR
  METIS_LIBRARY_DEBUG
  METIS_USE_PKGCONFIG
  METIS_CONFIG)

# ================================================================================================ #
# Register imported libraries:
#
# 1 If we can find a Windows .dll file (or if we can find both Debug and Release libraries), we will
# set appropriate target properties for these. 2 However, for most systems, we will only register
# the import location and include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string(REPLACE ".lib" ".dll" METIS_LIBRARY_DLL "${METIS_LIBRARY}")
  string(REPLACE ".lib" ".dll" METIS_LIBRARY_DEBUG_DLL "${METIS_LIBRARY_DEBUG}")
endif()

if(METIS_FOUND AND NOT TARGET METIS::metis)
  if(EXISTS "${METIS_LIBRARY_DLL}")

    # Windows systems with dll libraries.
    add_library(METIS::metis SHARED IMPORTED)

    # Windows with dlls, but only Release libraries.
    set_target_properties(
      METIS::metis
      PROPERTIES IMPORTED_LOCATION_RELEASE "${METIS_LIBRARY_DLL}"
                 IMPORTED_IMPLIB "${METIS_LIBRARY}"
                 INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIRS}"
                 IMPORTED_CONFIGURATIONS Release
                 IMPORTED_LINK_INTERFACE_LANGUAGES "C")

    # If we have both Debug and Release libraries
    if(EXISTS "${METIS_LIBRARY_DEBUG_DLL}")
      set_property(
        TARGET METIS::metis
        APPEND
        PROPERTY IMPORTED_CONFIGURATIONS Debug)
      set_target_properties(
        METIS::metis PROPERTIES IMPORTED_LOCATION_DEBUG "${METIS_LIBRARY_DEBUG_DLL}"
                                IMPORTED_IMPLIB_DEBUG "${METIS_LIBRARY_DEBUG}")
    endif()

  else()

    # For all other environments (ones without dll libraries), create the imported library targets.
    add_library(METIS::metis UNKNOWN IMPORTED)
    set_target_properties(
      METIS::metis
      PROPERTIES IMPORTED_LOCATION "${METIS_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES
                                                      "${METIS_INCLUDE_DIRS}"
                 IMPORTED_LINK_INTERFACE_LANGUAGES "C")
  endif()
endif()
