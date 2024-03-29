#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(cdi CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB sources *.cc)
file(GLOB explicit_instantiations *_pt.cc)
file(GLOB headers *.hh)

set(include_dirs PUBLIC $<BUILD_INTERFACE:${draco_src_dir_BINARY_DIR}/device>)

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
if(HAVE_GPU AND USE_GPU)
  if(EXISTS "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    string(APPEND include_dirs ";$<BUILD_INTERFACE:${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}>")
    string(APPEND include_dirs ";$<INSTALL_INTERFACE:${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}>")
  elseif(EXISTS "$ENV{ROCM_PATH}")
    string(APPEND include_dirs ";$<BUILD_INTERFACE:$ENV{ROCM_PATH}/include>")
    string(APPEND include_dirs ";$<INSTALL_INTERFACE:$ENV{ROCM_PATH}/include>")
  endif()
endif()

add_component_library(
  TARGET Lib_cdi
  TARGET_DEPS Lib_dsxx
  INCLUDE_DIRS "${include_dirs}"
  LIBRARY_NAME cdi
  HEADERS "${headers}"
  SOURCES "${sources}")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi)

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
# End cdi/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
