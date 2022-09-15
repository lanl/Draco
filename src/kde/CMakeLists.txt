#--------------------------------------------*-cmake-*---------------------------------------------#
#
# file   kde/CMakeLists.txt
# author Mathew Cleveland
# date   Nov. 10th 2020
# brief  CMake build environment for KDE
# note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(kde CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

file(GLOB sources *.cc)
file(GLOB explicit_instantiations *_pt.cc)
file(GLOB headers *.hh)

# Make the header files available in the IDE.
if(MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode)
  list(APPEND sources ${headers})
endif()

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #

add_component_library(
  TARGET Lib_kde
  TARGET_DEPS Lib_c4
  INCLUDE_DIRS ${include_dirs}
  LIBRARY_NAME ${PROJECT_NAME}
  HEADERS "${headers}"
  SOURCES "${sources}")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #

install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/kde)

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
# End kde/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
