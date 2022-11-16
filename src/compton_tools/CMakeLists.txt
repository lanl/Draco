#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/compton_tools/CMakeLists.txt
# author Andrew Till
# date   2020 May 11
# brief  Generate build project files for compton_tools
# note   Copyright (C) 2020-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(compton_tools CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
set(sources ${PROJECT_SOURCE_DIR}/Compton_Native.cc ${PROJECT_SOURCE_DIR}/Compton_Edep.cc)
set(headers ${PROJECT_SOURCE_DIR}/Compton_Native.hh ${PROJECT_SOURCE_DIR}/Compton_Edep.hh)

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_compton_tools
  TARGET_DEPS "Lib_c4;Lib_units;Lib_cdi_analytic"
  INCLUDE_DIRS "PUBLIC;$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>"
  LIBRARY_NAME compton_tools
  SOURCES "${sources}"
  HEADERS "${headers}")

add_component_executable(
  TARGET Exe_CSKrw
  TARGET_DEPS "Lib_units;Lib_compton_tools"
  SOURCES ${PROJECT_SOURCE_DIR}/cskrw.cc)

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
set(CMAKE_INSTALL_DO_STRIP NO)
install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/compton_tools)
install(
  TARGETS Exe_CSKrw
  EXPORT draco-targets
  DESTINATION ${DBSCFGDIR}bin)
install(PROGRAMS python/common_compton.py python/csk_reader.py python/merge_csk.py
        DESTINATION ${DBSCFGDIR}bin)

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
# End compton_tools/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
