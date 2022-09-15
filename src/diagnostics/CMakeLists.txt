#--------------------------------------------*-cmake-*---------------------------------------------#
# file   diagnostics/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for diagnostics.
# note   Copyright (C) 2010-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(diagnostics CXX)

# ------------------------------------------------------------------------------------------------ #
# Is Caliper available? If so, let's use it
# ------------------------------------------------------------------------------------------------ #
if(TARGET CALIPER::caliper)
  set(DRACO_CALIPER "1")
else()
  set(DRACO_CALIPER "0")
endif()
set(DRACO_CALIPER
    ${DRACO_CALIPER}
    CACHE BOOL "Does Draco require the Caliper library?" FORCE)

# ------------------------------------------------------------------------------------------------ #
# Special options for this component
# ------------------------------------------------------------------------------------------------ #

# DRACO_DIAGNOSICS is set in ds++ (ds++/config.h).

# From autoconf build system: "yes" --> 2, "no" --> 0
set(DRACO_TIMING
    "0"
    CACHE STRING "Timing detail level (0-2). Default is 0.")
# provide a constrained drop down menu for cmake-gui
set_property(CACHE DRACO_TIMING PROPERTY STRINGS 0 1 2)

# Are the memory features of getrusage to be activated?
option(USE_PROCMON "Use PROCMON process (cpu/memory) monitoring." OFF)
# If USE_PROCMON is requested but /proc/meminfo is unavailable, abort
if(USE_PROCMON)
  if(NOT EXISTS "/proc/meminfo")
    message(FATAL_ERROR "USE_PROCMON=ON, but /proc/meminfo unavailable.")
  endif()
endif()
if(NOT ${HAVE_SYS_TYPES})
  message(FATAL_ERROR "diagnostics/Procmon.hh requires sys/types.h.")
endif()

# ------------------------------------------------------------------------------------------------ #
# Generate config.h (only occurs when cmake is run)
# ------------------------------------------------------------------------------------------------ #
configure_file(config.h.in ${PROJECT_BINARY_DIR}/diagnostics/config.h)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB sources *.cc)
file(GLOB headers *.hh)
list(APPEND headers ${PROJECT_BINARY_DIR}/diagnostics/config.h)

if(CMAKE_CONFIGURATION_TYPES)
  set_source_files_properties(${PROJECT_SOURCE_DIR}/draco_info.cc
                              PROPERTIES COMPILE_DEFINITIONS CBT="${CMAKE_CFG_INTDIR}")
else()
  set_source_files_properties(${PROJECT_SOURCE_DIR}/draco_info.cc
                              PROPERTIES COMPILE_DEFINITIONS CBT="${CMAKE_BUILD_TYPE}")
endif()

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
set(deps Lib_c4)
if(TARGET CALIPER::caliper)
  list(APPEND deps CALIPER::caliper)
endif()

add_component_library(
  TARGET Lib_diagnostics
  TARGET_DEPS "${deps}"
  INCLUDE_DIRS "PUBLIC;$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>"
  LIBRARY_NAME ${PROJECT_NAME}
  SOURCES "${sources}"
  HEADERS "${headers}")

add_component_executable(
  TARGET Exe_draco_info
  TARGET_DEPS Lib_diagnostics
  EXE_NAME draco_info
  SOURCES ${PROJECT_SOURCE_DIR}/bin/draco_info_main.cc
  HEADERS ${PROJECT_SOURCE_DIR}/draco_info.hh)

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(
  TARGETS Exe_draco_info
  EXPORT draco-targets
  DESTINATION ${DBSCFGDIR}bin)
install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/diagnostics)

# ------------------------------------------------------------------------------------------------ #
# Unit tests
# ------------------------------------------------------------------------------------------------ #
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# ------------------------------------------------------------------------------------------------ #
# Autodoc
# ------------------------------------------------------------------------------------------------ #
process_autodoc_pages()

# ------------------------------------------------------------------------------------------------ #
# End diagnostics/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
