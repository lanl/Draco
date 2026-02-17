#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cpt_data/CMakeLists.txt
# author Kyle S. Beling
# date   2025 Nov 18
# brief  Instructions for building cpt_data Makefiles.
# note   Copyright (C) 2026 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(cpt_data CXX)

# Create cpt_data/config.h
configure_file(config.h.in ${PROJECT_BINARY_DIR}/cpt_data/config.h)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB sources *.cc)
file(GLOB headers *.hh)

# add config file to list of headers
list(APPEND headers_public "${PROJECT_BINARY_DIR}/cpt_data/config.h")

set(include_dirs PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>)
# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
set(target_deps Lib_units)

if(JSON_FOUND)
  list(APPEND include_dirs ${JSON_INCLUDE_DIR})
endif()

add_component_library(
  TARGET Lib_cpt_data
  TARGET_DEPS "${target_deps}"
  INCLUDE_DIRS ${include_dirs}
  LIBRARY_NAME ${PROJECT_NAME}
  SOURCES "${sources}"
  HEADERS "${headers}")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/cpt_data)

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
# End cpt_data/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
