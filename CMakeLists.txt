#--------------------------------------------*-cmake-*---------------------------------------------#
# file   draco/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building root level Makefile.
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
#
# Build notes:
#
# * https://re-git.lanl.gov/draco/draco/-/wikis/home
# * https://re-git.lanl.gov/draco/draco/-/wiki/Common_Configure_Options
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0 FATAL_ERROR)
string(CONCAT ddesc "An object-oriented component library supporting radiation transport "
              "applications.")
project(
  Draco
  DESCRIPTION ${ddesc}
  VERSION 7.15.0
  LANGUAGES CXX C)

# Export compile commands for compiler tools
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# Do not look for Fortran/CUDA/HIP for
#
# 1. XCode based Generators, or
# 2. Visual Studio IDE or NMake Generators (MSYS or CYGWIN environments will look for Fortran).
# 3. Ninja, Codeblocks or Eclipse CDT4 generators.
# 4. Unix Makefile generator when:
#
#   * Fortran only: $ENV{FC} is not set (e.g. clang on Linux)
#   * Cuda only   : $ENV{CUDADIR} is not set (e.g. cuda module not loaded)
if((${CMAKE_GENERATOR} MATCHES "Unix Makefiles" OR ${CMAKE_GENERATOR} MATCHES "Ninja") AND NOT MSVC)
  # For non-MSVC builds, always look for Fortran.  If $ENV{FC} or ${CMAKE_Fortran_COMPILER} are
  # defined, they will be used.  Otherwise cmake will search for a Fortran compiler based on its
  # internal rules.  Normally, it will look for gfortran first.
  enable_language(Fortran OPTIONAL)
  if(DEFINED ENV{CUDADIR}
     OR DEFINED ENV{CUDACXX}
     OR DEFINED ENV{CUDA_HOME})
    if((DEFINED USE_GPU AND USE_GPU) OR NOT DEFINED USE_GPU)
      enable_language(CUDA OPTIONAL)
    endif()
  elseif(
    "rocmcc" MATCHES $ENV{LOADEDMODULES}
    OR DEFINED ENV{HIPCXX}
    OR DEFINED CMAKE_HIP_COMPILER)
    if((DEFINED USE_GPU AND USE_GPU) OR NOT DEFINED USE_GPU)
      enable_language(HIP OPTIONAL)
      if(CMAKE_HIP_COMPILER MATCHES "CC")
        set(CMAKE_HIP_FLAGS "-x hip")
      endif()
    endif()
  endif()
endif()

# Build system configuration files are located here.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# Debug cmake scripts:
#
# 1. Review macros defined at config/debug_macros.cmake and config/print_target_properties.cmake.
# 2. Uncomment this feature to tell CMake to print include paths during target registration
#    set(CMAKE_DEBUG_TARGET_PROPERTIES INCLUDE_DIRECTORIES)
# 3. Optionally use --trace

#
# The Draco version number.
#
include(dracoVersion)
set_ccs2_software_version(Draco)

#
# Unit Test Setup
#
include(dracoTesting)

# Set default compile environment:

# Setup defaults, value checks, etc.
include(buildEnv)
# set defaults for BUILD_TYPE and INSTALL_PREFIX
dbssetdefaults()
# Initialize fields that define the exported target files (draco-config.cmake)
dbsinitexporttargets("Draco")
# Save config info
dbsconfiginfo()

# Platform Checks: Is HOST_NAME_MAX defined?  Is WinSock2.h available?  Is gethostname() available?
include(platform_checks)
set_draco_uname()
query_have_gethostname()
query_have_maxpathlen()
query_have_sys_headers() #< sets HAVE_UNISTD_H, etc.
query_have_restrict_keyword()
query_fma_on_hardware()

# Find any globally required libraries
include(FeatureSummary)
include(vendor_libraries)
setupvendorlibraries()

# Set compiler options
include(compilerEnv)
dbssetupcompilers()
dbssetupcxx()
dbssetupfortran()
dbssetupcuda()
dbssetupprofilertools()
dbssetupstaticanalyzers()

#
# Build Draco components:
#
add_subdirectory(src)
if(TARGET Exe_draco_info)
  add_subdirectory(autodoc) #< This must be processed after 'src'
endif()
add_subdirectory(config) #< This must be processed after 'src'

# install top level documents
install(FILES ChangeLog LICENSE.md README.md DESTINATION ${CMAKE_INSTALL_PREFIX})

# Export targets
install(
  EXPORT draco-targets
  DESTINATION cmake
  EXPORT_LINK_INTERFACE_LIBRARIES)

# ------------------------------------------------------------------------------------------------ #
# End of CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
