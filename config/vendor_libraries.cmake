# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   config/vendor_libraries.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 6
# brief  Look for any libraries which are required at the top level.
# note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved.
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
  set_package_properties(
    PythonInterp PROPERTIES
    URL "https://www.python.org"
    DESCRIPTION "Python interpreter"
    TYPE REQUIRED
    PURPOSE "Required for running tests and accessing features that rely on matplotlib.")
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
# Helper macros for Random123
# ------------------------------------------------------------------------------------------------ #
macro(setupRandom123)

  message(STATUS "Looking for Random123...")
  find_package(Random123 REQUIRED QUIET)
  mark_as_advanced(RANDOM123_FOUND)
  if(RANDOM123_FOUND)
    message(STATUS "Looking for Random123.found ${RANDOM123_INCLUDE_DIR}")
  else()
    message(STATUS "Looking for Random123.not found")
  endif()
  set_package_properties(
    Random123 PROPERTIES
    URL "http://www.deshawresearch.com/resources_random123.html"
    DESCRIPTION "a library of counter-based random number generators"
    TYPE REQUIRED
    PURPOSE "Required for building the rng component.")
endmacro()

# --------------------------------------------------------------------------------------------------
# Helper macros for LAPACK/Unix
#
# This module sets the following variables:
#
# * BLAS_FOUND   - set to true if a library implementing the BLAS interface is found
# * LAPACK_FOUND - set to true if a library implementing the LAPACK interface is found
#
# Provides import targets:
#
# * LAPACK::LAPACK
# * BLAS::BLAS
#
# Providers:
#
# * Linux - use spack to install netlib-lapack, openblas, or mkl, https://github.com/spack/spack
# * Windows - Use vcpkg, https://github.com/microsoft/vcpkg, or clone and build from sources,
#   https://github.com/KineticTheory/lapack-visualstudio-mingw-gfortran
#
# Special instructions:
#
# * cmake-3.21+ will find Fujitsu_SSL2 when using xlc++ on Power9.  As of Sep 2021, this provider
#   seems broken on darwin and rzansel. To avoid the bad behavior, we set BLA_VENDOR=Generic to help
#   cmake find netlib-lapack.
# --------------------------------------------------------------------------------------------------
function(setupLAPACKLibraries)
  message(STATUS "Looking for LAPACK {netlib, mkl, openblas}...")
  if(CMAKE_CXX_COMPILER_ID STREQUAL XLClang)
    # If xlc++, we need to use netlib-lapack to avoid wierd issues,
    # https://re-git.lanl.gov/draco/draco/-/issues/1361
    set(BLA_VENDOR "Generic")
  endif()
  if(NOT TARGET BLAS::BLAS)
    find_package(BLAS QUIET)
  endif()
  if(NOT TARGET LAPACK::LAPACK)
    find_package(LAPACK QUIET)
  endif()
  if(TARGET LAPACK::LAPACK)
    target_link_libraries(LAPACK::LAPACK INTERFACE BLAS::BLAS)
  endif()
  set(lapack_url "http://www.netlib.org/lapack")
  if(LAPACK_LIBRARIES MATCHES openblas)
    set(lapack_url "http://openblas.net")
  elseif(LAPACK_LIBRARIES MATCHES mkl)
    set(lapack_url "https://software.intel.com/en-us/intel-mkl")
  endif()
  set_package_properties(
    BLAS PROPERTIES
    URL "${lapack_url}"
    DESCRIPTION "Basic Linear Algebra Subprograms"
    TYPE OPTIONAL
    PURPOSE "Required for building the lapack_wrap component.")
  set_package_properties(
    lapack PROPERTIES
    URL "${lapack_url}"
    DESCRIPTION "Linear Algebra PACKage"
    TYPE OPTIONAL
    PURPOSE "Required for building the lapack_wrap component.")
  message(STATUS "Looking for LAPACK {netlib, mkl, openblas}...found ${LAPACK_LIBRARIES} "
                 "${BLAS_LIBRARIES}")
endfunction()

# --------------------------------------------------------------------------------------------------
# Setup QT (any)
# --------------------------------------------------------------------------------------------------
macro(setupQt)
  message(STATUS "Looking for Qt SDK....")

  if(CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv)
    option(USE_QT "Build QT support for Draco" OFF)
  else()
    option(USE_QT "Build QT support for Draco" ON)
  endif()

  if(USE_QT)
    # Find the QtWidgets library
    find_package(Qt5 COMPONENTS Widgets QUIET)

    if(Qt5Core_DIR)
      mark_as_advanced(Qt5Core_DIR Qt5Gui_DIR Qt5Gui_EGL_LIBRARY Qt5Widgets_DIR QTDIR)
      message(STATUS "Looking for Qt SDK....found ${Qt5Core_DIR}")
    else()
      message(STATUS "Looking for Qt SDK....not found.")
    endif()

    set_package_properties(
      Qt PROPERTIES
      URL "http://qt.io"
      DESCRIPTION "Qt is a comprehensive cross-platform C++ application framework."
      TYPE OPTIONAL
      PURPOSE "Only needed to demo qt version of draco_diagnostics.")

  endif()

endmacro()

# --------------------------------------------------------------------------------------------------
# Setup GSL (any)
# --------------------------------------------------------------------------------------------------
macro(setupGSL)

  if(NOT TARGET GSL::gsl)

    message(STATUS "Looking for GSL...")
    set(QUIET "QUIET")

    # There are 3 ways to find gsl:

    # * Config mode.
    #
    # If CMAKE_PREFIX_PATH contains a GSL install prefix directory and the file gsl-config.cmake is
    # found somewhere in this installation tree, then the targets defined by gsl-config.cmake will
    # be used.
    find_package(GSL CONFIG ${QUIET})

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
          "${GSL_CONFIG}" ARGS
          --prefix OUTPUT_VARIABLE
          GSL_ROOT_DIR)
      endif()
    endif()

    # * Module mode.
    #
    # Locate GSL by using the value of GSL_ROOT_DIR or by looking in standard locations. We add
    # 'REQUIRED' here because if this fails, then we abort the built.
    find_package(GSL REQUIRED ${QUIET})

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

  # If successful in finding GSL, provide some information for the vendor summary reported by
  # src/CMakeLists.txt.
  if(TARGET GSL::gsl)
    # Include some information that can be printed by the build system.
    set_package_properties(
      GSL PROPERTIES
      URL "https://www.gnu.org/software/gsl"
      DESCRIPTION
        "The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers."
      TYPE REQUIRED
      PURPOSE "Required for rng and quadrature components.")
  endif()
  unset(QUIET)

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Setup ParMETIS (any)
# ------------------------------------------------------------------------------------------------ #
macro(setupParMETIS)

  set(QUIET "QUIET")
  if(NOT TARGET METIS::metis)
    message(STATUS "Looking for METIS...")

    find_package(METIS CONFIG ${QUIET})
    if(NOT TARGET METIS::metis)
      find_package(METIS ${QUIET})
    endif()
    if(TARGET METIS::metis)
      if(TARGET METIS::metis AND NOT METIS_LIBRARY)
        foreach(config NOCONFIG DEBUG RELEASE RELWITHDEBINFO)
          get_target_property(tmp METIS::metis IMPORTED_LOCATION_${config})
          if(EXISTS ${tmp} AND NOT METIS_LIBRARY)
            set(METIS_LIBRARY ${tmp})
          endif()
        endforeach()
      endif()
      message(STATUS "Looking for METIS.....found ${METIS_LIBRARY}")
    else()
      message(STATUS "Looking for METIS.....not found")
    endif()

    # Include some information that can be printed by the build system.
    set_package_properties(
      METIS PROPERTIES
      DESCRIPTION "METIS"
      TYPE RECOMMENDED
      URL "http://glaros.dtc.umn.edu/gkhome/metis/metis/overview"
      PURPOSE
        "METIS is a set of serial programs for partitioning graphs, partitioning finite element
meshes, and producing fill reducing orderings for sparse matrices.")

  endif()

  if(NOT TARGET ParMETIS::parmetis)

    message(STATUS "Looking for ParMETIS...")

    find_package(ParMETIS QUIET)
    if(ParMETIS_FOUND)
      message(STATUS "Looking for ParMETIS..found ${ParMETIS_LIBRARY}")
    else()
      message(STATUS "Looking for ParMETIS..not found")
    endif()

    # Include some information that can be printed by the build system.
    set_package_properties(
      ParMETIS PROPERTIES
      DESCRIPTION "MPI Parallel METIS"
      TYPE OPTIONAL
      URL "http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview"
      PURPOSE "ParMETIS is an MPI-based parallel library that implements a variety of algorithms for
partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of sparse
matrices.")

  endif()
  unset(QUIET)
endmacro()

# ------------------------------------------------------------------------------
# Setup Libquo (https://github.com/lanl/libquo
# ------------------------------------------------------------------------------
macro(setupLIBQUO)

  if(NOT TARGET LIBQUO::libquo AND MPI_C_FOUND)
    message(STATUS "Looking for LIBQUO...")

    find_package(Libquo QUIET)

    if(LIBQUO_FOUND)
      message(STATUS "Looking for LIBQUO....found ${LIBQUO_LIBRARY}")
    else()
      message(STATUS "Looking for LIBQUO....not found")
    endif()

    # Include some information that can be printed by the build system.
    set_package_properties(
      Libquo PROPERTIES
      URL "https://github.com/lanl/libquo"
      DESCRIPTION
        "A runtime library that aids in accommodating thread-level heterogeneity in dynamic,
phased MPI+X appliations comprising single- and multi-threaded libraries."
      TYPE RECOMMENDED
      PURPOSE
        "Required for allowing draco-clients to switch MPI+X bindings and thread affinities when a
library is called instead of at program ivokation.")
  endif()

endmacro()

# ------------------------------------------------------------------------------
# Setup Caliper (https://github.com/LLNL/Caliper)
# ------------------------------------------------------------------------------
macro(setupCaliper)

  if(NOT TARGET CALIPER::caliper)
    message(STATUS "Looking for Caliper...")
    find_package(Caliper QUIET)
    if(CALIPER_FOUND)
      message(STATUS "Looking for Caliper...${CALIPER_LIBRARY}")
    else()
      message(STATUS "Looking for Caliper...not found")
    endif()
  endif()

endmacro()

# ------------------------------------------------------------------------------
# Setup Eospac (https://laws.lanl.gov/projects/data/eos.html)
# ------------------------------------------------------------------------------
macro(setupEOSPAC)

  if(NOT TARGET EOSPAC::eospac)
    message(STATUS "Looking for EOSPAC...")

    find_package(EOSPAC QUIET)

    if(EOSPAC_FOUND)
      message(STATUS "Looking for EOSPAC....found ${EOSPAC_LIBRARY}")
    else()
      message(STATUS "Looking for EOSPAC....not found")
    endif()

    # Include some information that can be printed by the build system.
    set_package_properties(
      EOSPAC PROPERTIES
      URL "https://laws.lanl.gov/projects/data/eos.html"
      DESCRIPTION "Access SESAME thermodynamic and transport data."
      TYPE OPTIONAL
      PURPOSE "Required for bulding the cdi_eospac component.")
  endif()

endmacro()

# ------------------------------------------------------------------------------
# Setup NDI (https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html)
# ------------------------------------------------------------------------------
macro(setupNDI)

  if(NOT TARGET NDI::ndi)
    message(STATUS "Looking for NDI...")

    find_package(NDI QUIET)

    if(NDI_FOUND)
      message(STATUS "Looking for NDI....found ${NDI_LIBRARY}")
    else()
      message(STATUS "Looking for NDI....not found")
    endif()

    # Include some information that can be printed by the build system.
    set_package_properties(
      NDI PROPERTIES
      URL "https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html"
      DESCRIPTION "Access nuclear data."
      TYPE OPTIONAL
      PURPOSE "Required for building the cdi_ndi component.")
  endif()

endmacro()

# ------------------------------------------------------------------------------
# Setup COMPTON (https://gitlab.lanl.gov/keadyk/CSK_generator)
# ------------------------------------------------------------------------------
macro(setupCOMPTON)

  if(NOT TARGET compton::compton)
    message(STATUS "Looking for COMPTON...")

    find_package(COMPTON QUIET)

    if(COMPTON_FOUND)
      message(STATUS "Looking for COMPTON...found ${COMPTON_LIBRARY}")
    else()
      message(STATUS "Looking for COMPTON...not found")
    endif()

    # Include some information that can be printed by the build system.
    set_package_properties(
      COMPTON PROPERTIES
      URL "https://gitlab.lanl.gov/CSK/CSK"
      DESCRIPTION "Access multigroup Compton scattering data."
      TYPE OPTIONAL
      PURPOSE "Required for building the Compton component.")
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Helper macros for setup_global_libraries()
# ------------------------------------------------------------------------------------------------ #
macro(SetupVendorLibrariesUnix)

  setupgsl()
  setupparmetis()
  setupcompton()
  setupeospac()
  setupndi()
  setuprandom123()
  setuppython()
  setupqt()
  setuplibquo()
  setupcaliper()

  # Doxygen ------------------------------------------------------------------
  message(STATUS "Looking for Doxygen...")
  find_package(Doxygen QUIET OPTIONAL_COMPONENTS dot mscgen)
  set_package_properties(
    Doxygen PROPERTIES
    URL "http://www.stack.nl/~dimitri/doxygen"
    DESCRIPTION "Doxygen autodoc generator"
    TYPE OPTIONAL
    PURPOSE "Required for building develop HTML documentation.")
  if(DOXYGEN_FOUND)
    message(STATUS "Looking for Doxygen...found version ${DOXYGEN_VERSION}")
  else()
    message(STATUS "Looking for Doxygen...not found")
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Vendors for building on Windows-based platforms.
# ------------------------------------------------------------------------------------------------ #
macro(SetupVendorLibrariesWindows)

  setupgsl()
  setupparmetis()
  setuprandom123()
  setupcompton()
  setupeospac()
  setupndi()
  setuppython()
  setupqt()

  # Doxygen ------------------------------------------------------------------
  message(STATUS "Looking for Doxygen...")
  find_package(Doxygen QUIET OPTIONAL_COMPONENTS dot mscgen)
  set_package_properties(
    Doxygen PROPERTIES
    URL "http://www.stack.nl/~dimitri/doxygen"
    DESCRIPTION "Doxygen autodoc generator"
    TYPE OPTIONAL
    PURPOSE "Required for building develop HTML documentation.")
  if(DOXYGEN_FOUND)
    message(STATUS "Looking for Doxygen...found version ${DOXYGEN_VERSION}")
  else()
    message(STATUS "Looking for Doxygen...not found")
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# Helper macros for setup_global_libraries() Assign here the library version to be used.
# ------------------------------------------------------------------------------------------------ #
macro(setVendorVersionDefaults)
  # Set the preferred search directories(ROOT)

  # Check that VENDOR_DIR is defined as a cache variable or as an environment variable. If defined
  # as both then take the environment variable.

  # See if VENDOR_DIR is set.  Try some defaults if it is not set.
  if(NOT DEFINED VENDOR_DIR AND IS_DIRECTORY "$ENV{VENDOR_DIR}")
    set(VENDOR_DIR $ENV{VENDOR_DIR})
  endif()
  # If needed, try some obvious places.
  if(NOT EXISTS "${VENDOR_DIR}")
    if(IS_DIRECTORY "/ccs/codes/radtran/vendors/Linux64")
      set(VENDOR_DIR "/ccs/codes/radtran/vendors/Linux64")
    endif()
    if(IS_DIRECTORY /usr/projects/draco/vendors)
      set(VENDOR_DIR /usr/projects/draco/vendors)
    endif()
    if(IS_DIRECTORY c:/vendors)
      set(VENDOR_DIR c:/vendors)
    endif()
  endif()
  # Cache the result
  if(IS_DIRECTORY "${VENDOR_DIR}")
    set(VENDOR_DIR
        ${VENDOR_DIR}
        CACHE PATH "Root directory where CCS-2 3rd party libraries are located." FORCE)
  endif()

  # Import environment variables related to vendors
  #
  # 1. Use command line variables (-DLAPACK_LIB_DIR=<path>
  # 2. Use environment variables ($ENV{LAPACK_LIB_DIR}=<path>)
  # 3. Try to find vendor in $VENDOR_DIR
  # 4. Don't set anything and let the user set a value in the#cache after failed 1st configure
  #    attempt.
  if(NOT DEFINED LAPACK_LIB_DIR AND IS_DIRECTORY $ENV{LAPACK_LIB_DIR})
    set(LAPACK_LIB_DIR $ENV{LAPACK_LIB_DIR})
    set(LAPACK_INC_DIR $ENV{LAPACK_INC_DIR})
  endif()
  if(NOT LAPACK_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/lapack-3.4.2/lib)
    set(LAPACK_LIB_DIR "${VENDOR_DIR}/lapack-3.4.2/lib")
    set(LAPACK_INC_DIR "${VENDOR_DIR}/lapack-3.4.2/include")
  endif()

  if(NOT DEFINED GSL_LIB_DIR)
    if(IS_DIRECTORY $ENV{GSL_LIB_DIR})
      set(GSL_LIB_DIR $ENV{GSL_LIB_DIR})
      set(GSL_INC_DIR $ENV{GSL_INC_DIR})
    elseif(IS_DIRECTORY ${VENDOR_DIR}/gsl/lib)
      set(GSL_LIB_DIR "${VENDOR_DIR}/gsl/lib")
      set(GSL_INC_DIR "${VENDOR_DIR}/gsl/include")
    endif()
  endif()

  if(NOT DEFINED ParMETIS_ROOT_DIR)
    if(IS_DIRECTORY $ENV{ParMETIS_ROOT_DIR})
      set(ParMETIS_ROOT_DIR $ENV{ParMETIS_ROOT_DIR})
    endif()
  endif()

  if(NOT DEFINED RANDOM123_INC_DIR AND IS_DIRECTORY $ENV{RANDOM123_INC_DIR})
    set(RANDOM123_INC_DIR $ENV{RANDOM123_INC_DIR})
  endif()
  if(NOT DEFINED RANDOM123_INC_DIR AND IS_DIRECTORY ${VENDOR_DIR}/Random123-1.08/include)
    set(RANDOM123_INC_DIR "${VENDOR_DIR}/Random123-1.08/include")
  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# This macro should contain all the system libraries which are required to link the main objects.
# ------------------------------------------------------------------------------------------------ #
macro(setupVendorLibraries)

  message("\nVendor Setup:\n")

  #
  # General settings
  #
  setvendorversiondefaults()
  if(NOT TARGET lapack)
    setuplapacklibraries()
  endif()

  setupmpilibraries()
  # System specific settings
  if(UNIX)
    setupvendorlibrariesunix()
  elseif(WIN32)
    setupvendorlibrarieswindows()
  else()
    message(
      FATAL_ERROR
        "\nI don't know how to setup global (vendor) libraries for this platform.  WIN32=0; "
        "UNIX=0;CMAKE_SYSTEM=${CMAKE_SYSTEM}; CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}")
  endif()

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
  include( vendor_libraries )

  # Provide targets for MPI, Metis, etc.
  setupVendorLibraries()

  # Set compiler options
  include( compilerEnv )
  dbsSetupCxx()
  dbsSetupFortran()
  dbsSetupCuda()
  dbsSetupProfilerTools()

  # CMake macros like 'add_component_library' and 'add_component_executable'
  include( component_macros )

endmacro()

")

  message(" ")

endmacro()

# ------------------------------------------------------------------------------------------------ #
# End vendor_libraries.cmake
# ------------------------------------------------------------------------------------------------ #
