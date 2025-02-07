#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_cpt/CMakeLists.txt
# author Ben. R. Ryan <brryan@lanl.gov>
# date   2024 Oct 9
# brief  Instructions for building cdi_cpt Makefiles.
# note   Copyright (C) 2025 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(cdi_cpt CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB sources *.cc)
file(GLOB headers *.hh)
list(APPEND headers ${PROJECT_BINARY_DIR}/cdi_cpt/config.h)

# ------------------------------------------------------------------------------------------------ #
# Generate config.h (only occurs when cmake is run)
# ------------------------------------------------------------------------------------------------ #
configure_file(config.h.in ${PROJECT_BINARY_DIR}/cdi_cpt/config.h)
set(include_dirs PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>)

# ------------------------------------------------------------------------------------------------ #
# Set library dependencies
# ------------------------------------------------------------------------------------------------ #
set(deps Lib_cdi Lib_units Lib_c4)

# If SHACL_PHYSICS is found, set the associated variables. Otherwise, the cdi_cpt lib will be built
# with stubbed-out routines.
if(TARGET shacl::physics)
  list(APPEND deps shacl::physics)
endif()

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_cdi_cpt
  TARGET_DEPS "${deps}"
  INCLUDE_DIRS "${include_dirs}"
  LIBRARY_NAME ${PROJECT_NAME}
  SOURCES "${sources}"
  HEADERS "${headers}")

target_compile_definitions(Lib_cdi_cpt PUBLIC AVALANCHE_DIR="${avalanche_DIR}")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_cpt)

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
# End cdi_cpt/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
