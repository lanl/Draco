#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/FindNDI.cmake
# date   2017 February 28
# brief  Instructions for discovering the NDI vendor libraries.
# note   Copyright (C) 2020-2025 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# .rst: FindNDI
# ---------
#
# Find the LANL NDI library and header files.
#
# A collection of C routines that can be used to access NDI data files.
# https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If NDI is found, this module defines the following :prop_tgt:`IMPORTED` targets::
#
# NDI::ndi        - The NDI library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
# * NDI_FOUND          - True if NDI found on the local system
# * NDI_INCLUDE_DIRS   - Location of NDI header files.
# * NDI_LIBRARIES      - The NDI libraries.
# * NDI_VERSION        - The version of the discovered NDI install.
#
# Hints ^^^^^
#
# Set ``NDI_ROOT_DIR`` to a directory that contains a NDI installation.
#
# This script expects to find libraries at ``$NDI_ROOT_DIR/lib`` and the NDI headers at
# ``$NDI_ROOT_DIR/include``.  The library directory may optionally provide Release and Debug
# folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of NDI installation
# discovered.  These variables may optionally be set to help this module find the correct files::
#
# * NDI_LIBRARY        - Location of the NDI library.
# * NDI_LIBRARY_DEBUG  - Location of the debug NDI library (if any).
#
# ------------------------------------------------------------------------------------------------ #

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

# =============================================================================
# If the user has provided ``NDI_ROOT_DIR``, use it!  Choose items found at this location over
# system locations.
if(EXISTS "$ENV{NDI_ROOT_DIR}")
  file(TO_CMAKE_PATH "$ENV{NDI_ROOT_DIR}" NDI_ROOT_DIR)
  set(NDI_ROOT_DIR
      "${NDI_ROOT_DIR}"
      CACHE PATH "Prefix for NDI installation.")
endif()

# =============================================================================
# Set NDI_INCLUDE_DIRS and NDI_LIBRARIES. Try to find the libraries at $NDI_ROOT_DIR (if provided)
# or in standard system locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found at the HINTS
# location, standard system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path(
  NDI_INCLUDE_DIR
  NAMES ndi.h
  HINTS ${NDI_ROOT_DIR}/include/ndi ${NDI_ROOT_DIR}/include
  PATH_SUFFIXES ndi)
set(NDI_LIBRARY_NAME ndipic ndi)

find_library(
  NDI_LIBRARY
  NAMES ${NDI_LIBRARY_NAME}
  PATHS ${NDI_ROOT_DIR}/lib)
# Do we also have debug versions?
find_library(
  NDI_LIBRARY_DEBUG
  NAMES ${NDI_LIBRARY_NAME}
  HINTS ${NDI_ROOT_DIR}/lib
  PATH_SUFFIXES Debug)
set(NDI_INCLUDE_DIRS ${NDI_INCLUDE_DIR})
set(NDI_LIBRARIES ${NDI_LIBRARY})

# Try to find the version.
if(NOT NDI_VERSION)
  if(EXISTS "${NDI_INCLUDE_DIR}/ndi_version.h")
    # pull data directly from the include file.
    file(STRINGS "${NDI_INCLUDE_DIR}/ndi_version.h" include_version_h REGEX "NDI_VERSION_DATA")
    string(REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([A-Za-z0-9]+).*" "\\1" NDI_MAJOR
                         ${include_version_h})
    string(REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([A-Za-z0-9]+).*" "\\2" NDI_MINOR
                         ${include_version_h})
    string(REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([A-Za-z0-9]+).*" "\\3" NDI_SUBMINOR
                         ${include_version_h})
  else()
    # We might also try scraping the directory name for a regex match "ndi-X.X.X"
    string(REGEX REPLACE ".*ndi-([0-9]+).([0-9]+).([A-Za-z0-9]+).*" "\\1" NDI_MAJOR
                         ${NDI_INCLUDE_DIR})
    string(REGEX REPLACE ".*ndi-([0-9]+).([0-9]+).([A-Za-z0-9]+).*" "\\2" NDI_MINOR
                         ${NDI_INCLUDE_DIR})
    string(REGEX REPLACE ".*ndi-([0-9]+).([0-9]+).([A-Za-z0-9]+).*" "\\3" NDI_SUBMINOR
                         ${NDI_INCLUDE_DIR})
  endif()
  if(NDI_MAJOR)
    set(NDI_VERSION "${NDI_MAJOR}.${NDI_MINOR}.${NDI_SUBMINOR}")
    set(NDI_VERSION_STRING "${NDI_VERSION}")
  endif()
endif()

# =============================================================================
# handle the QUIETLY and REQUIRED arguments and set NDI_FOUND to TRUE if all listed variables are
# TRUE.
find_package_handle_standard_args(
  NDI
  FOUND_VAR NDI_FOUND
  REQUIRED_VARS NDI_INCLUDE_DIR NDI_LIBRARY NDI_VERSION
  VERSION_VAR NDI_VERSION)

mark_as_advanced(
  NDI_ROOT_DIR
  NDI_VERSION
  NDI_LIBRARY
  NDI_INCLUDE_DIR
  NDI_LIBRARY_DEBUG
  NDI_USE_PKGCONFIG
  NDI_CONFIG)

# =============================================================================
# ~~~
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and Release libraries),
#    we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and include directory.
# ~~~

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string(REPLACE ".lib" ".dll" NDI_LIBRARY_DLL "${NDI_LIBRARY}")
  string(REPLACE ".lib" ".dll" NDI_LIBRARY_DEBUG_DLL "${NDI_LIBRARY_DEBUG}")
endif()

if(NDI_FOUND AND NOT TARGET NDI::ndi)
  if(WIN32)
    if(EXISTS "${NDI_LIBRARY_DLL}")

      # Windows systems with dll libraries.
      add_library(NDI::ndi SHARED IMPORTED)

      # Windows with dlls, but only Release libraries.
      set_target_properties(
        NDI::ndi
        PROPERTIES IMPORTED_LOCATION_RELEASE "${NDI_LIBRARY_DLL}"
                   IMPORTED_IMPLIB "${NDI_LIBRARY}"
                   INTERFACE_INCLUDE_DIRECTORIES "${NDI_INCLUDE_DIRS}"
                   IMPORTED_CONFIGURATIONS Release
                   IMPORTED_LINK_INTERFACE_LANGUAGES "C")

      # If we have both Debug and Release libraries
      if(EXISTS "${NDI_LIBRARY_DEBUG_DLL}")
        set_property(
          TARGET NDI::ndi
          APPEND
          PROPERTY IMPORTED_CONFIGURATIONS Debug)
        set_target_properties(NDI::ndi PROPERTIES IMPORTED_LOCATION_DEBUG "${NDI_LIBRARY_DEBUG_DLL}"
                                                  IMPORTED_IMPLIB_DEBUG "${NDI_LIBRARY_DEBUG}")
      endif()

    else()
      # Windows systems with static lib libraries.
      add_library(NDI::ndi STATIC IMPORTED)

      # Windows with dlls, but only Release libraries.
      set_target_properties(
        NDI::ndi
        PROPERTIES IMPORTED_LOCATION_RELEASE "${NDI_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES
                                                              "${NDI_INCLUDE_DIRS}"
                   IMPORTED_CONFIGURATIONS Release IMPORTED_LINK_INTERFACE_LANGUAGES "C")

      # If we have both Debug and Release libraries
      if(EXISTS "${NDI_LIBRARY_DEBUG}")
        set_property(
          TARGET NDI::ndi
          APPEND
          PROPERTY IMPORTED_CONFIGURATIONS Debug)
        set_target_properties(NDI::ndi PROPERTIES IMPORTED_LOCATION_DEBUG "${NDI_LIBRARY_DEBUG}")
      endif()

    endif()

  else()

    # For all other environments (ones without dll libraries), create the imported library targets.
    add_library(NDI::ndi UNKNOWN IMPORTED)
    set_target_properties(
      NDI::ndi
      PROPERTIES IMPORTED_LOCATION "${NDI_LIBRARY}" INTERFACE_INCLUDE_DIRECTORIES
                                                    "${NDI_INCLUDE_DIRS}"
                 IMPORTED_LINK_INTERFACE_LANGUAGES "C")
  endif()

  # This is a workaround until we find a better way to handle CP2011 data, which is not included in
  # the standard gendir.all file. It also appears to be patchily installed on our systems.
  if(NDI_FOUND AND "${NDI_CP2011_FOUND}notset" STREQUAL "notset")
    set(NDI_CP2011_FOUND "FALSE")
    set(NDI_CP2011_GENDIR "NOT-FOUND")
    if(${SITENAME_FAMILY} MATCHES "ATS-3") # roci, tycho, xroads
      set(NDI_CP2011_FOUND "TRUE")
      set(NDI_CP2011_GENDIR "/usr/projects/data/nuclear/ndi/special/cp2011/gendir")
    elseif(${SITENAME_FAMILY} MATCHES "CCS-NET")
      set(NDI_CP2011_FOUND "TRUE")
      set(NDI_CP2011_GENDIR "/ccs/opt/vendors-ec/ndi/special/cp2011/gendir")
    endif()
  endif()
  # Set and cache CMake variables:
  set(NDI_CP2011_FOUND
      "${NDI_CP2011_FOUND}"
      CACHE BOOL "NDI CP2011 data found?" FORCE)
  set(NDI_CP2011_GENDIR
      "${NDI_CP2011_GENDIR}"
      CACHE PATH "Path for NDI CP2011 data." FORCE)
endif()

# ------------------------------------------------------------------------------------------------ #
# End FindNDI.cmake
# ------------------------------------------------------------------------------------------------ #
