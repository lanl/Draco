#--------------------------------------------*-cmake-*---------------------------------------------#
# file   rng/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng.
# note   Copyright (C) 2010-2023 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.21.0)
project(rng CXX)

# Local variables ---------------------------------------------------------------------------------#
string(CONCAT header_prefix_source "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
              "$<INSTALL_INTERFACE:${DBSCFG_IMPORT_PREFIX}/rng>")
string(CONCAT header_prefix_binary "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>"
              "$<INSTALL_INTERFACE:${DBSCFG_IMPORT_PREFIX}>")

# Determine if Random123 version > 1.09.
if(EXISTS ${RANDOM123_INCLUDE_DIR}/Random123/u01fixedpt.h)
  set(R123_VERSION 1.13.2)
else()
  set(R123_VERSION 1.09)
endif()

# ------------------------------------------------------------------------------------------------ #
# Generate config.h (only occurs when cmake is run)
# ------------------------------------------------------------------------------------------------ #
set(R123_USE_CXX11 1)
set(R123_USE_STD_RANDOM 1)
set(R123_USE_64BIT 1)
configure_file(config.h.in ${PROJECT_BINARY_DIR}/rng/config.h)
set(pkg_config_h "${header_prefix_binary}/rng/config.h")

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(
  GLOB sources
  RELATIVE "${PROJECT_SOURCE_DIR}"
  *.cc *.c)
file(
  GLOB headers
  RELATIVE "${PROJECT_SOURCE_DIR}"
  *.hh *.h *.hpp)

# prefix header files with appropriate paths (build vs install).
foreach(file ${headers})
  list(APPEND headers_public ${header_prefix_source}/${file})
endforeach()
list(APPEND headers_public "${pkg_config_h}")

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
if(HAVE_GPU AND USE_GPU)
  set(deps Lib_dsxx Lib_device GSL::gsl)
else()
  set(deps Lib_dsxx GSL::gsl)
endif()

# remove 'src/device' when 'device/CMakeFiles.txt' is updated?
set(include_dirs PUBLIC $<BUILD_INTERFACE:${draco_src_dir_BINARY_DIR}/device>
                 $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> ${RANDOM123_INCLUDE_DIR})
add_component_library(
  TARGET Lib_rng
  TARGET_DEPS "${deps}"
  INCLUDE_DIRS "${include_dirs}"
  SOURCES "${sources}"
  HEADERS "${headers_public}"
  LIBRARY_NAME ${PROJECT_NAME})
# NOEXPORT ) # explict install instructions are below.

# Target sources should list headers as `PUBLIC` but this is broken as of cmake-3.17.0 (might
# actually be a Visual Studio issue?). Because this is broken we mark headers as `PRIVATE`.
# target_sources(Lib_rng PRIVATE "${headers_public}" PRIVATE "${sources}" )
# set_target_properties(Lib_rng PROPERTIES PUBLIC_HEADER "${headers_public}") if(TARGET Objlib_rng)
# target_sources(Objlib_rng PRIVATE "${headers_public}" PRIVATE "${sources}" )
# set_target_properties(Objlib_rng PROPERTIES PUBLIC_HEADER "${headers_public}") endif()

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
# install( TARGETS Lib_rng EXPORT draco-targets LIBRARY DESTINATION "${DBSCFGDIR_LIBRARY}" ARCHIVE
# DESTINATION "${DBSCFGDIR}lib" RUNTIME DESTINATION "${DBSCFGDIR}bin" PUBLIC_HEADER DESTINATION
# "${DBSCFGDIR}include/rng" )
install(FILES ${headers_public} DESTINATION ${DBSCFGDIR}include/rng)

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
# End rng/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
