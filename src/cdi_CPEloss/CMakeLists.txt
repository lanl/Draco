#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_CPEloss/CMakeLists.txt
# author Ben. R. Ryan <brryan@lanl.gov>
# date   2019 Nov 4
# brief  Instructions for building cdi_CPEloss Makefiles.
# note   Copyright (C) 2020-2025 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(cdi_CPEloss CXX)

# Create cdi_CPEloss/config.h
configure_file(config.h.in ${PROJECT_BINARY_DIR}/cdi_CPEloss/config.h)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB sources *.cc)
file(GLOB headers *.hh)

# add config file to list of headers
list(APPEND headers_public "${PROJECT_BINARY_DIR}/cdi_CPEloss/config.h")

set(include_dirs PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>)
# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
set(target_deps Lib_cdi;Lib_units)
if(TARGET GSL::gsl)
  list(APPEND target_deps GSL::gsl)
endif()

add_component_library(
  TARGET Lib_cdi_CPEloss
  TARGET_DEPS "${target_deps}"
  INCLUDE_DIRS ${include_dirs}
  LIBRARY_NAME ${PROJECT_NAME}
  SOURCES "${sources}"
  HEADERS "${headers}")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_CPEloss)

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
# End cdi_CPEloss/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
