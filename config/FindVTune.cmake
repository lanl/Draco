#--------------------------------------------*-cmake-*---------------------------------------------#
# file   draco/config/FindVTune.cmake
# author Alex Long <along@lanl.gov>
# brief  Helper for finding VTune profiling api.
# note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# ~~~
# .rst:
# FindVTune
# --------
#
# Find the VTune ittnotify library, if a vtune module is loaded (profiling only)
#
# VTune is a profiler by Intel that can measure various aspects of performance, including hotpots
# vectorization and use of the memory subsystem. It includes an API that allows the user to control
# when VTune collects information. Right now the version of VTune is hardcoded in this file as the
# VTune modules don't set a uniform "VTUNE_ROOT_DIR" type variable to look for. The profiling can
# be controlled with the API calls to ittnotify like this:
#
#  __itt_resume();
#   run_important_code();
#  __itt_pause();
#
# Where these calls are in the the "ittnotify.h" file. Add the --start-paused" on the VTune command
# line for the simplest approach (alternatively you could immediately pause in your program and then
# resume around important segments).

# * Find VTune ittnotify. Defines: VTune_FOUND VTune_INCLUDE_DIRS VTune_LIBRARIES
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If VTune is found, this module defines the following :prop_tgt:`IMPORTED` targets::
#
# VTune::vtune      - The VTune ittnotify ibrary.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
# VTune_FOUND          - True if VTune found on the local system VTune_INCLUDE_DIRS   - Location of
# VTune ittnotify header files. VTune_LIBRARIES      - The location of the ittnotify VTune library.
#

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

# only try to find VTune if the module is loaded, try to avoid linking unnecessary libraries when
# not profiling an application
set(VTUNE_MODULE_LOADED "FALSE")
if(DEFINED ENV{VTUNE_PROFILER_2022_DIR} AND EXISTS $ENV{VTUNE_PROFILER_2022_DIR})
  set(VTUNE_MODULE_LOADED "TRUE")
  file(TO_CMAKE_PATH "$ENV{VTUNE_PROFILER_2022_DIR}" VTUNE_PROFILER_2022_DIR)
  set(VTUNE_ROOT_DIR
      "${VTUNE_PROFILER_2022_DIR}"
      CACHE PATH "Prefix for VTune installation.")
elseif(DEFINED ENV{VTUNE_AMPLIFIER_XE_2020_DIR} AND EXISTS $ENV{VTUNE_AMPLIFIER_XE_2020_DIR})
  set(VTUNE_MODULE_LOADED "TRUE")
  file(TO_CMAKE_PATH "$ENV{VTUNE_AMPLIFIER_XE_2020_DIR}" VTUNE_AMPLIFIER_XE_2020_DIR)
  set(VTUNE_ROOT_DIR
      "${VTUNE_AMPLIFIER_XE_2020_DIR}"
      CACHE PATH "Prefix for VTune installation.")

endif()

if(VTUNE_MODULE_LOADED STREQUAL "TRUE")

  find_path(
    VTune_INCLUDE_DIRS
    HINTS ${VTUNE_ROOT_DIR}/include ${VTUNE_ROOT_DIR}/sdk/src/ittnotify/
    NAMES ittnotify.h)

  find_library(VTune_LIBRARIES libittnotify.a HINTS ${VTUNE_ROOT_DIR}/lib64/)

  # ============================================================================================== #
  # handle the QUIETLY and REQUIRED arguments and set VTune_FOUND to TRUE if all listed variables
  # are TRUE.
  find_package_handle_standard_args(
    VTune
    FOUND_VAR VTune_FOUND
    REQUIRED_VARS VTune_INCLUDE_DIRS VTune_LIBRARIES)

  # For all other environments without dll libraries, create the imported library targets.
  if(VTune_FOUND AND NOT TARGET VTune::vtune)
    add_library(VTune::vtune UNKNOWN IMPORTED)
    set_target_properties(
      VTune::vtune
      PROPERTIES IMPORTED_LOCATION "${VTune_LIBRARIES}" INTERFACE_INCLUDE_DIRECTORIES
                                                        "${VTune_INCLUDE_DIRS}"
                 IMPORTED_LINK_INTERFACE_LANGUAGES "C")

  endif()
endif()
