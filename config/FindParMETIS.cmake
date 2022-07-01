#--------------------------------------------*-cmake-*---------------------------------------------#
# file   draco/config/FindParMETIS.cmake
# author Kelly Thompson <kgt@lanl.gov>
# brief  Helper for finding parmetis.
# note   Copyright (C) 2016-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# ~~~
# .rst:
# FindParMetis
# --------
#
# Find the ParMETIS and METIS includes and libraries.
#
# ParMETIS is an MPI-based parallel library that implements a variety of
# algorithms for partitioning unstructured graphs, meshes, and for computing
# fill-reducing orderings of sparse matrices. ParMETIS extends the functionality
# provided by METIS and includes routines that are especially suited for
# parallel AMR computations and large scale numerical simulations.  The
# algorithms implemented in ParMETIS are based on the parallel multilevel k-way
# graph-partitioning, adaptive repartitioning, and parallel multi- constrained
# partitioning schemes developed in our lab.
# http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If ParMETIS is found, this module defines the following :prop_tgt:`IMPORTED`
# targets::
#
#  ParMETIS::parmetis      - The main ParMETIS library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  ParMETIS_FOUND          - True if ParMETIS found on the local system
#  ParMETIS_INCLUDE_DIRS   - Location of ParMETIS header files.
#  ParMETIS_LIBRARIES      - The ParMETIS libraries.
#  ParMETIS_VERSION        - The version of the discovered ParMETIS install.
#
# Hints
# ^^^^^
#
# Set ``ParMETIS_ROOT_DIR`` to a directory that contains a ParMETIS
# installation.
#
# This script expects to find libraries at ``$ParMETIS_ROOT_DIR/lib`` and the
# ParMETIS headers at ``$ParMETIS_ROOT_DIR/include``.  The library directory may
# optionally provide Release and Debug folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of
# ParMETIS installation discovered.  These variables may optionally be set to
# help this module find the correct files::
#
#  ParMETIS_LIBRARY             - Location of the ParMETIS library.
#  ParMETIS_LIBRARY_DEBUG       - Location of the debug ParMETIS library
#                                 (if any).
#
# ~~~

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

# ================================================================================================ #
# Prerequisite: ParMETIS requires that METIS be found.
find_package(METIS QUIET)

# ================================================================================================ #
# If the user has provided ``ParMETIS_ROOT_DIR``, use it!  Choose items found at this location over
# system locations.
if(EXISTS "$ENV{ParMETIS_ROOT_DIR}")
  file(TO_CMAKE_PATH "$ENV{ParMETIS_ROOT_DIR}" ParMETIS_ROOT_DIR)
  set(ParMETIS_ROOT_DIR
      "${ParMETIS_ROOT_DIR}"
      CACHE PATH "Prefix for ParMETIS installation.")
endif()

# ================================================================================================ #
# Set ParMETIS_INCLUDE_DIRS and ParMETIS_LIBRARIES. Try to find the libraries at $ParMETIS_ROOT_DIR
# (if provided) or in standard system locations.  These find_library and find_path calls will prefer
# custom locations over standard locations (HINTS).  If the requested file is not found at the HINTS
# location, standard system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path(
  ParMETIS_INCLUDE_DIR
  NAMES parmetis.h
  HINTS ${ParMETIS_ROOT_DIR}/include ${ParMETIS_INCLUDEDIR}
  PATH_SUFFIXES Release Debug)

find_library(
  ParMETIS_LIBRARY
  NAMES parmetis
  HINTS ${ParMETIS_ROOT_DIR}/lib ${ParMETIS_LIBDIR}
  PATH_SUFFIXES Release Debug)
# Do we also have debug versions?
find_library(
  ParMETIS_LIBRARY_DEBUG
  NAMES parmetis
  HINTS ${ParMETIS_ROOT_DIR}/Debug/lib ${ParMETIS_LIBDIR})
set(ParMETIS_INCLUDE_DIRS ${ParMETIS_INCLUDE_DIR})
set(ParMETIS_LIBRARIES ${ParMETIS_LIBRARY})

# Try to find the version.
if(NOT ParMETIS_VERSION)
  if(EXISTS "${ParMETIS_INCLUDE_DIRS}/parmetis.h")
    file(STRINGS "${ParMETIS_INCLUDE_DIRS}/parmetis.h" parmetis_h_major
         REGEX "define PARMETIS_MAJOR")
    file(STRINGS "${ParMETIS_INCLUDE_DIRS}/parmetis.h" parmetis_h_minor
         REGEX "define PARMETIS_MINOR")
    file(STRINGS "${ParMETIS_INCLUDE_DIRS}/parmetis.h" parmetis_h_subminor
         REGEX "define PARMETIS_SUBMINOR")
    string(REGEX REPLACE ".*([0-9]+)" "\\1" ParMETIS_MAJOR ${parmetis_h_major})
    string(REGEX REPLACE ".*([0-9]+)" "\\1" ParMETIS_MINOR ${parmetis_h_minor})
    string(REGEX REPLACE ".*([0-9]+)" "\\1" ParMETIS_SUBMINOR ${parmetis_h_subminor})
  endif()
  # We might also try scraping the directory name for a regex match "parmetis-X.X.X"
endif()

# ================================================================================================ #
# handle the QUIETLY and REQUIRED arguments and set ParMETIS_FOUND to TRUE if all listed variables
# are TRUE.
find_package_handle_standard_args(
  ParMETIS
  FOUND_VAR ParMETIS_FOUND
  REQUIRED_VARS ParMETIS_INCLUDE_DIR ParMETIS_LIBRARY METIS_FOUND
  VERSION_VAR ParMETIS_VERSION)

mark_as_advanced(
  ParMETIS_ROOT_DIR
  ParMETIS_VERSION
  ParMETIS_LIBRARY
  ParMETIS_INCLUDE_DIR
  ParMETIS_LIBRARY_DEBUG
  ParMETIS_USE_PKGCONFIG
  ParMETIS_CONFIG)

# ================================================================================================ #
# Register imported libraries: 1. If we can find a Windows .dll file (or if we can find both Debug
# and Release libraries), we will set appropriate target properties for these. 2. However, for most
# systems, we will only register the import location and include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string(REPLACE ".lib" ".dll" ParMETIS_LIBRARY_DLL "${ParMETIS_LIBRARY}")
  string(REPLACE ".lib" ".dll" ParMETIS_LIBRARY_DEBUG_DLL "${ParMETIS_LIBRARY_DEBUG}")
endif()

if(ParMETIS_FOUND AND NOT TARGET ParMETIS::parmetis)
  if(WIN32)
    if(EXISTS "${ParMETIS_LIBRARY_DLL}")
      set(libtype SHARED)
      set(ilr "${ParMETIS_LIBRARY_DLL}")
      set(ild "${ParMETIS_LIBRARY_DEBUG_DLL}")
    else()
      set(libtype STATIC)
      set(ilr "${ParMETIS_LIBRARY}")
      set(ild "${ParMETIS_LIBRARY_DEBUG}")
    endif()

    add_library(ParMETIS::parmetis ${libtype} IMPORTED)

    set_target_properties(
      ParMETIS::parmetis
      PROPERTIES IMPORTED_LOCATION_RELEASE "${ilr}"
                 IMPORTED_IMPLIB "${ParMETIS_LIBRARY}"
                 INTERFACE_INCLUDE_DIRECTORIES "${ParMETIS_INCLUDE_DIRS}"
                 IMPORTED_CONFIGURATIONS Release
                 IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                 INTERFACE_LINK_LIBRARIES METIS::metis)

    # If we have both Debug and Release libraries
    if(EXISTS "${ild}")
      set_property(
        TARGET ParMETIS::parmetis
        APPEND
        PROPERTY IMPORTED_CONFIGURATIONS Debug)
      set_target_properties(
        ParMETIS::parmetis PROPERTIES IMPORTED_LOCATION_DEBUG "${ild}" IMPORTED_IMPLIB_DEBUG
                                                                       "${ParMETIS_LIBRARY_DEBUG}")
    endif()
    unset(libytpe)
    unset(ilr)
    unset(ild)

  else() # Linux

    # For all other environments (ones without dll libraries), create the imported library targets.
    add_library(ParMETIS::parmetis UNKNOWN IMPORTED)
    set_target_properties(
      ParMETIS::parmetis
      PROPERTIES IMPORTED_LOCATION "${ParMETIS_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES
                                                         "${ParMETIS_INCLUDE_DIRS}"
                 IMPORTED_LINK_INTERFACE_LANGUAGES "C" INTERFACE_LINK_LIBRARIES METIS::metis)
  endif()
endif()
