# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/vendor_libraries.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 6
# brief  Look for any libraries which are required at the top level.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

include_guard(GLOBAL)
include(FeatureSummary)
include(setupMPI) # defines the macros setupMPILibrariesUnix|Windows

# ------------------------------------------------------------------------------------------------ #
# Helper macros for Python
# ------------------------------------------------------------------------------------------------ #
macro(setupPython)

  message(STATUS "Looking for Python....")
  # This module looks preferably for version 3 of Python. If not found, version 2 is searched.
  find_package(Python QUIET REQUIRED COMPONENTS Interpreter)
  # * Python_Interpreter_FOUND - Was the Python executable found
  # * Python_EXECUTABLE  - path to the Python interpreter
  if(Python_Interpreter_FOUND)
    message(STATUS "Looking for Python....found ${Python_EXECUTABLE}")
  else()
    message(STATUS "Looking for Python....not found")
  endif()
  # As of 2020-11-10, we require 'python@3.6:' for correct dictionary sorting features.  We can
  # build the code with 'python@2:' but some tests may fail.
  if(Python_VERSION_MAJOR STREQUAL "3" AND Python_VERSION_MINOR VERSION_LESS "6")
    message(
      FATAL_ERROR
        "When using python3, we require version 3.6+.  Python version "
        "${Python_VERSION} was discovered, which doesn't satisfy the compatibility requirement.")
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Helper functions for Random123
# ------------------------------------------------------------------------------------------------ #
function(setupRandom123)

  message(STATUS "Looking for Random123...")
  find_package(Random123 REQUIRED QUIET)
  mark_as_advanced(RANDOM123_FOUND)
  if(RANDOM123_FOUND)
    message(STATUS "Looking for Random123.found ${RANDOM123_INCLUDE_DIR}")
  else()
    message(STATUS "Looking for Random123.not found")
  endif()
endfunction()

# --------------------------------------------------------------------------------------------------
# Setup GSL (any)
# --------------------------------------------------------------------------------------------------
function(setupGSL)

  if(NOT TARGET GSL::gsl)

    message(STATUS "Looking for GSL...")

    # There are 3 ways to find gsl:

    # * Config mode.
    #
    # If CMAKE_PREFIX_PATH contains a GSL install prefix directory and the file gsl-config.cmake is
    # found somewhere in this installation tree, then the targets defined by gsl-config.cmake will
    # be used.
    find_package(GSL CONFIG QUIET)

  endif()

  if(NOT TARGET GSL::gsl) # if option #1 was successful, skip this.

    # * pkg-config mode (Linux)
    #
    # If GSL_ROOT_DIR isn't set, look for the binary 'gsl-config' in $PATH. If found, run it to
    # discover and set GSL_ROOT_DIR that will be used in method #3.

    if("$ENV{GSL_ROOT_DIR}x" STREQUAL "x" AND "${GSL_ROOT_DIR}x" STREQUAL "x")
      find_program(GSL_CONFIG gsl-config)
      if(EXISTS "${GSL_CONFIG}")
        exec_program(
          "${GSL_CONFIG}" "."
          ARGS --prefix
          OUTPUT_VARIABLE GSL_ROOT_DIR)
      endif()
    endif()

    # * Module mode.
    #
    # Locate GSL by using the value of GSL_ROOT_DIR or by looking in standard locations. We add
    # 'REQUIRED' here because if this fails, then we abort the built.
    find_package(GSL REQUIRED QUIET)

  endif()

  # Print a report
  if(TARGET GSL::gsl)
    if(TARGET GSL::gsl AND NOT GSL_LIBRARY)
      foreach(config NOCONFIG DEBUG RELEASE RELWITHDEBINFO)
        get_target_property(tmp GSL::gsl IMPORTED_LOCATION_${config})
        if(EXISTS ${tmp} AND NOT GSL_LIBRARY)
          set(GSL_LIBRARY ${tmp})
        endif()
      endforeach()
    endif()
    message(STATUS "Looking for GSL.......found ${GSL_LIBRARY}")
    mark_as_advanced(GSL_CONFIG_EXECUTABLE)
  else()
    message(STATUS "Looking for GSL.......not found")
  endif()

endfunction()

# ------------------------------------------------------------------------------------------------ #
# Setup ParMETIS (any)
# ------------------------------------------------------------------------------------------------ #
function(setupParMETIS)

  if(VCPKG_INSTALLED_DIR)
    set(thispkg "metis")
    set(thispkgmode CONFIG)
  else()
    set(thispkg "METIS")
    set(thispkgmode MODULE)
  endif()

  set(QUIET "QUIET")
  if(NOT TARGET METIS::metis AND NOT TARGET metis)
    message(STATUS "Looking for METIS...")
    find_package(${thispkg} ${thispkgmode} ${QUIET})
  endif()

  if(VCPKG_INSTALLED_DIR AND TARGET metis)
    add_library(METIS::metis ALIAS metis)
  endif()

  set(METIS_FOUND
      "${METIS_FOUND}"
      CACHE STRING "Was the metis library found?")

  # Set METIS_LIBRARY variable
  if(TARGET METIS::metis AND NOT METIS_LIBRARY)
    foreach(config RELEASE NOCONFIG DEBUG RELWITHDEBINFO)
      get_target_property(tmp ${thispkg} IMPORTED_LOCATION_${config})
      if(EXISTS ${tmp} AND NOT METIS_LIBRARY)
        set(METIS_LIBRARY ${tmp})
      endif()
    endforeach()
  endif()

  if(EXISTS "${METIS_LIBRARY}")
    message(STATUS "Looking for METIS.....found ${METIS_LIBRARY}")
  else()
    message(STATUS "Looking for METIS.....not found")
  endif()

  if(NOT TARGET ParMETIS::parmetis)
    message(STATUS "Looking for ParMETIS...")
    find_package(ParMETIS QUIET)
    if(ParMETIS_FOUND)
      message(STATUS "Looking for ParMETIS..found ${ParMETIS_LIBRARY}")
    else()
      message(STATUS "Looking for ParMETIS..not found")
    endif()
  endif()

  set(ParMETIS_FOUND
      "${ParMETIS_FOUND}"
      CACHE STRING "Was the parmetis library found?")

  unset(QUIET)
  unset(thispkg)
  unset(thispkgmode)

endfunction()

# ------------------------------------------------------------------------------
# Setup VTune
# ------------------------------------------------------------------------------
function(setupVTune)
  if(NOT TARGET VTune::vtune)
    message(STATUS "Looking for VTune...")
    find_package(VTune QUIET)
    if(VTune_FOUND)
      message(STATUS "Looking for VTune...found ${VTune_LIBRARIES}")
    else()
      message(STATUS "Looking for VTune...not found")
    endif()
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Setup Libquo (https://github.com/lanl/libquo
# ------------------------------------------------------------------------------
function(setupLIBQUO)
  if(NOT TARGET LIBQUO::libquo AND TARGET MPI::MPI_C)
    message(STATUS "Looking for LIBQUO...")
    find_package(Libquo QUIET)
    if(LIBQUO_FOUND)
      message(STATUS "Looking for LIBQUO....found ${LIBQUO_LIBRARY}")
    else()
      message(STATUS "Looking for LIBQUO....not found")
    endif()
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Setup Caliper (https://github.com/LLNL/Caliper)
# ------------------------------------------------------------------------------
function(setupCaliper)

  if(NOT TARGET CALIPER::caliper)
    message(STATUS "Looking for Caliper...")
    find_package(Caliper QUIET)
    if(CALIPER_FOUND)
      message(STATUS "Looking for Caliper...${CALIPER_LIBRARY}")
    else()
      message(STATUS "Looking for Caliper...not found")
    endif()
  endif()

endfunction()

# ------------------------------------------------------------------------------
# Setup Eospac (https://laws.lanl.gov/projects/data/eos.html)
# ------------------------------------------------------------------------------
function(setupEOSPAC)

  if(NOT TARGET EOSPAC::eospac)
    message(STATUS "Looking for EOSPAC...")
    find_package(EOSPAC QUIET)
    if(EOSPAC_FOUND)
      message(STATUS "Looking for EOSPAC....found ${EOSPAC_LIBRARY}")
    else()
      message(STATUS "Looking for EOSPAC....not found")
    endif()
    set(EOSPAC_FOUND
        "${EOSPAC_FOUND}"
        CACHE STRING "Did we find libeospac?")
  endif()

endfunction()

# ------------------------------------------------------------------------------
# Setup NDI (https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html)
# ------------------------------------------------------------------------------
function(setupNDI)

  if(NOT TARGET NDI::ndi)
    message(STATUS "Looking for NDI...")
    find_package(NDI QUIET)
    if(NDI_FOUND)
      message(STATUS "Looking for NDI....found ${NDI_LIBRARY}")
    else()
      message(STATUS "Looking for NDI....not found")
    endif()
  endif()
  set(NDI_FOUND
      "${NDI_FOUND}"
      CACHE STRING "Did we find the NDI library?")
  set(NDI_VERSION_STRING
      "${NDI_VERSION_STRING}"
      CACHE STRING ".")
  set(NDI_MAJOR
      "${NDI_MAJOR}"
      CACHE STRING ".")
  set(NDI_MINOR
      "${NDI_MINOR}"
      CACHE STRING ".")

endfunction()

# ------------------------------------------------------------------------------------------------ #
# This function should contain all the system libraries which are required to link the main objects.
# ------------------------------------------------------------------------------------------------ #
macro(setupVendorLibraries)

  message("\nVendor Setup:\n")

  setupmpilibraries()
  setupgsl()
  setupparmetis()
  setuprandom123()
  setupeospac()
  setupndi()
  setuppython()
  setupvtune()
  if(UNIX)
    # These are not supported on Windows yet.
    setuplibquo()
    setupcaliper()
    setuptorch()
  endif()

  # Doxygen ------------------------------------------------------------------
  message(STATUS "Looking for Doxygen...")
  find_package(Doxygen QUIET OPTIONAL_COMPONENTS dot mscgen)
  if(DOXYGEN_FOUND)
    message(STATUS "Looking for Doxygen...found version ${DOXYGEN_VERSION}")
  else()
    message(STATUS "Looking for Doxygen...not found")
  endif()
  set(DOXYGEN_FOUND
      "${DOXYGEN_FOUND}"
      CACHE STRING "Is doxygen available?")

  # Add commands to draco-config.cmake (which is installed for use by other projects), to setup
  # Draco's vendors
  string(
    APPEND
    Draco_EXPORT_TARGET_PROPERTIES
    "
macro( dbs_basic_setup )

  message(\"
Looking for Draco...\")
  message(\"Looking for Draco...\${draco_DIR}
  \")

  # Provide helper functions used by component CMakeLists.txt files
  # This block of code generated by draco/config/vendor_libraries.cmake.

  # Setup defaults, value checks, etc.
  include(buildEnv)
  dbsSetDefaults()

  # CMake macros that check the system for features like 'gethostname', etc.
  include( platform_checks )

  # CMake macros to query the availability of TPLs.
  # - Must do this before adding compiler flags. Otherwise find_package(MPI) will fail for LLVM
  #   when using flags '-Werror -Weverything'.
  # include(vendor_libraries)

  # Provide targets for MPI, Metis, etc.
  setupVendorLibraries()

  # Set compiler options
  include(compilerEnv)
  dbsSetupCxx()
  dbsSetupFortran()
  dbsSetupCuda()
  dbsSetupProfilerTools()

  # CMake macros like 'add_component_library' and 'add_component_executable'
  include(component_macros)

endmacro()
")

  message(" ")

endmacro()

# ------------------------------------------------------------------------------
# Setup LibTorch (https://github.com/pytorch/pytorch)
# ------------------------------------------------------------------------------
function(setuptorch)

  if(NOT TARGET TORCH::torch)
    message(STATUS "Looking for Torch...")
    find_package(Torch QUIET)
    if(TORCH_FOUND)
      add_definitions(-DLIBTORCH)
      message(STATUS "Looking for Torch...${TORCH_LIBRARY}")
    else()
      message(STATUS "Looking for Torch...not found")
    endif()
  endif()

endfunction()

# ------------------------------------------------------------------------------------------------ #
# End vendor_libraries.cmake
# ------------------------------------------------------------------------------------------------ #
