#--------------------------------------------*-cmake-*---------------------------------------------#
# file   cdi_eospac/CMakeLists.txt
# brief  Instructions for building cdi_eospac Makefiles.
# note   Copyright (C) 2012-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project(cdi_eospac CXX)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "11.2.0")
  # Suppress build warning: "lto1: error: potential null pointer dereference"
  string(APPEND CMAKE_C_FLAGS " -Wno-null-dereference")
  string(APPEND CMAKE_CXX_FLAGS " -Wno-null-dereference")
endif()

if(EOSPAC_FOUND AND NOT DEFINED ENV{TRAVIS})

  # ---------------------------------------------------------------------------------------------- #
  # Source files
  # ---------------------------------------------------------------------------------------------- #
  file(GLOB sources *.cc)
  file(GLOB headers *.hh)
  list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/QueryEospac.cc)

  # ---------------------------------------------------------------------------------------------- #
  # Build package library
  # ---------------------------------------------------------------------------------------------- #
  add_component_library(
    TARGET Lib_cdi_eospac
    TARGET_DEPS "Lib_cdi;EOSPAC::eospac"
    INCLUDE_DIRS "PUBLIC;${EOSPAC_INCLUDE_DIR}"
    LIBRARY_NAME cdi_eospac
    SOURCES "${sources}"
    HEADERS "${headers}")

  add_component_executable(
    TARGET Exe_QueryEospac
    TARGET_DEPS Lib_cdi_eospac
    SOURCES ${PROJECT_SOURCE_DIR}/QueryEospac.cc)

  # Copy necessary dll files to the build directory
  if(MSVC AND DRACO_LIBRARY_TYPE STREQUAL SHARED)
    add_custom_command(
      TARGET Lib_cdi_eospac
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:EOSPAC::eospac>
              $<TARGET_FILE_DIR:Lib_cdi_eospac>
      COMMENT "")
  endif()

  # ---------------------------------------------------------------------------------------------- #
  # Installation instructions
  # ---------------------------------------------------------------------------------------------- #
  install(
    TARGETS Exe_QueryEospac
    EXPORT draco-targets
    DESTINATION ${DBSCFGDIR}bin)
  install(FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_eospac)

  # ---------------------------------------------------------------------------------------------- #
  # Unit tests
  # ---------------------------------------------------------------------------------------------- #
  if(BUILD_TESTING AND NOT DEFINED ENV{APPVEYOR})
    add_subdirectory(test)
  endif()

  # ---------------------------------------------------------------------------------------------- #
  # Autodoc
  # ---------------------------------------------------------------------------------------------- #
  process_autodoc_pages()

endif() # EOSPAC_FOUND

# ------------------------------------------------------------------------------------------------ #
# End cdi_eospac/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
