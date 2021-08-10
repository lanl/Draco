#--------------------------------------------*-cmake-*---------------------------------------------#
# file   draco/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building root level Makefile.
# note   Copyright (C) 2019-2021 Triad National Security, LLC., All rights reserved.
#
# Build notes:
#
# * https://rtt.lanl.gov/redmine/projects/draco/wiki/
# * https://rtt.lanl.gov/redmine/projects/draco/wiki/Common_Configure_Options
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.18.0 FATAL_ERROR)
string(CONCAT ddesc "An object-oriented component library supporting radiation transport "
              "applications.")
project(
  Draco
  DESCRIPTION ${ddesc}
  VERSION 7.11
  LANGUAGES CXX C)

# Export compile commands for compiler tools
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# Do not look for Fortran/CUDA for
#
# 1. XCode based Generators, or
# 2. Visual Studio IDE or NMake Generators (MSYS or CYGWIN environments will look for Fortran).
# 3. Ninja, Codeblocks or Eclipse CDT4 generators.
# 4. Unix Makefile generator when:
#
#   * Fortran only: $ENV{FC} is not set (e.g. clang on Linux)
#   * Cuda only   : $ENV{CUDADIR} is not set (e.g. cuda module not loaded)
if(${CMAKE_GENERATOR} MATCHES "Unix Makefiles")
  if(DEFINED ENV{FC} OR DEFINED CMAKE_Fortran_COMPILER)
    enable_language(Fortran OPTIONAL)
  endif()
  if(DEFINED ENV{CUDADIR}
     OR DEFINED ENV{CUDACXX}
     OR DEFINED ENV{CUDA_HOME})
    if((DEFINED USE_CUDA AND USE_CUDA) OR NOT DEFINED USE_CUDA)
      enable_language(CUDA OPTIONAL)
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
